#include "ahf_ilucg_gpu.h"
#include <cusolverDn.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
std::size_t product(std::size_t a, std::size_t b) {
    if (b && a > std::numeric_limits<std::size_t>::max() / b)
        throw std::overflow_error("AHF LU: allocation size overflow");
    return a * b;
}
void cuda_check(cudaError_t error, const char* operation) {
    if (error != cudaSuccess)
        throw std::runtime_error(std::string("AHF LU: ") + operation + ": " +
                                 cudaGetErrorString(error));
}
void solver_check(cusolverStatus_t error, const char* operation) {
    if (error != CUSOLVER_STATUS_SUCCESS)
        throw std::runtime_error(std::string("AHF LU: ") + operation +
                                 ": status=" + std::to_string(int(error)));
}

bool overlaps(const void* a, std::size_t na, const void* b, std::size_t nb) {
    const auto x = reinterpret_cast<std::uintptr_t>(a);
    const auto y = reinterpret_cast<std::uintptr_t>(b);
    return x <= y ? y - x < na : x - y < nb;
}

__global__ void csr_to_dense(int n, const int* rows, const int* columns,
                             const double* values, double* dense) {
    const std::size_t row = (std::size_t(blockIdx.x) * blockDim.x + threadIdx.x) / 32;
    const unsigned lane = threadIdx.x % 32;
    if (row >= std::size_t(n)) return;
    for (std::size_t k = std::size_t(rows[row]) + lane;
         k < std::size_t(rows[row + 1]); k += 32)
        dense[row + std::size_t(columns[k]) * std::size_t(n)] = values[k];
}

__global__ void check_finite(std::size_t na, const double* a,
                              std::size_t nb, const double* b, int* bad) {
    const std::size_t count = na > nb ? na : nb;
    int invalid = 0;
    for (std::size_t k = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
         k < count; k += std::size_t(blockDim.x) * gridDim.x) {
        if ((k < na && !isfinite(a[k])) || (k < nb && !isfinite(b[k])))
            invalid = 1;
    }
    if (invalid) atomicExch(bad, 1);
}

int finite_blocks(int count) {
    const int blocks = (count - 1) / 256 + 1;
    return blocks < 4096 ? blocks : 4096;
}
cudaError_t read_scalar(cudaStream_t stream, const int* device, int& host) {
    cudaError_t error = cudaMemcpyAsync(&host, device, sizeof(host),
                                       cudaMemcpyDeviceToHost, stream);
    if (error != cudaSuccess) return error;
    return cudaStreamSynchronize(stream);
}
} // namespace

struct GpuAhfLinearSolver {
    cudaStream_t stream;
    int n, nnz;
    const int *rows, *columns;
    std::size_t dense_bytes, vector_bytes;
    cusolverDnHandle_t handle = nullptr;
    double *dense = nullptr, *workspace = nullptr;
    int *pivots = nullptr, *info = nullptr;
    // 主机侧稀疏迭代解算；失败时回退稠密 LU（AMSS_AHF_DENSE_LU=1 可强制 LU）
    std::vector<int> h_rows, h_columns;
    std::vector<double> h_values, h_rhs, h_sol, h_diag;
    std::vector<double> h_p, h_v, h_s, h_t, h_r, h_rhat;

    GpuAhfLinearSolver(cudaStream_t st, int size, int nonzeros,
                       const int* row_offsets, const int* column_indices,
                       std::size_t matrix_bytes, std::size_t rhs_bytes)
        : stream(st), n(size), nnz(nonzeros), rows(row_offsets), columns(column_indices),
          dense_bytes(matrix_bytes), vector_bytes(rhs_bytes) {}
    ~GpuAhfLinearSolver() noexcept {
        cudaStreamSynchronize(stream);
        if (handle) cusolverDnDestroy(handle);
        if (info) cudaFree(info);
        if (pivots) cudaFree(pivots);
        if (workspace) cudaFree(workspace);
        if (dense) cudaFree(dense);
    }
    GpuAhfLinearSolver(const GpuAhfLinearSolver&) = delete;
    GpuAhfLinearSolver& operator=(const GpuAhfLinearSolver&) = delete;
};

GpuAhfLinearSolver* gpu_ahf_linear_create(cudaStream_t stream, int n, int nnz,
                                        const int* d_rows, const int* d_columns) {
    if (n <= 0 || nnz <= 0 || !d_rows || !d_columns ||
        n > std::numeric_limits<int>::max() + GPU_AHF_LINEAR_SINGULAR_BASE)
        throw std::invalid_argument("AHF LU: invalid dimensions or CSR pointers");
    const std::size_t count = product(std::size_t(n), std::size_t(n));
    const std::size_t dense_bytes = product(count, sizeof(double));
    const std::size_t vector_bytes = product(std::size_t(n), sizeof(double));
    product(std::size_t(nnz), sizeof(double));
    product(std::size_t(n) + 1, sizeof(int));
    if (std::size_t(nnz) > count)
        throw std::invalid_argument("AHF LU: nnz exceeds matrix size");
    std::unique_ptr<GpuAhfLinearSolver> s(new GpuAhfLinearSolver(
        stream, n, nnz, d_rows, d_columns, dense_bytes, vector_bytes));
    solver_check(cusolverDnCreate(&s->handle), "cusolverDnCreate");
    solver_check(cusolverDnSetStream(s->handle, stream), "cusolverDnSetStream");
    cuda_check(cudaMalloc(reinterpret_cast<void**>(&s->dense), dense_bytes), "dense allocation");
    int workspace_count = 0;
    solver_check(cusolverDnDgetrf_bufferSize(s->handle, n, n, s->dense, n,
                                            &workspace_count), "Dgetrf workspace query");
    if (workspace_count < 0)
        throw std::runtime_error("AHF LU: negative Dgetrf workspace size");
    const std::size_t workspace_bytes = product(
        std::size_t(workspace_count > 0 ? workspace_count : 1), sizeof(double));
    cuda_check(cudaMalloc(reinterpret_cast<void**>(&s->workspace), workspace_bytes), "LU workspace allocation");
    cuda_check(cudaMalloc(reinterpret_cast<void**>(&s->pivots), product(std::size_t(n), sizeof(int))), "pivot allocation");
    cuda_check(cudaMalloc(reinterpret_cast<void**>(&s->info), sizeof(int)), "info allocation");
    s->h_rows.assign(std::size_t(n) + 1, 0);
    s->h_columns.assign(std::size_t(nnz), 0);
    cuda_check(cudaMemcpy(s->h_rows.data(), d_rows, (std::size_t(n) + 1) * sizeof(int),
                          cudaMemcpyDeviceToHost), "CSR rows download");
    cuda_check(cudaMemcpy(s->h_columns.data(), d_columns, std::size_t(nnz) * sizeof(int),
                          cudaMemcpyDeviceToHost), "CSR columns download");
    s->h_values.assign(std::size_t(nnz), 0.0);
    s->h_rhs.assign(std::size_t(n), 0.0);
    s->h_sol.assign(std::size_t(n), 0.0);
    s->h_diag.assign(std::size_t(n), 1.0);
    s->h_p.assign(std::size_t(n), 0.0);
    s->h_v.assign(std::size_t(n), 0.0);
    s->h_s.assign(std::size_t(n), 0.0);
    s->h_t.assign(std::size_t(n), 0.0);
    s->h_r.assign(std::size_t(n), 0.0);
    s->h_rhat.assign(std::size_t(n), 0.0);
    return s.release();
}

void gpu_ahf_linear_destroy(GpuAhfLinearSolver* solver) noexcept {
    delete solver;
}

namespace {
void csr_ax(const int *rows, const int *columns, const double *a, const double *x,
            double *y, int n) {
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int k = rows[i]; k < rows[i + 1]; ++k) sum += a[k] * x[columns[k]];
        y[i] = sum;
    }
}
double vdot(const double *a, const double *b, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; ++i) sum += a[i] * b[i];
    return sum;
}
// Jacobi 预条件 BiCGStab；true=已解出并写回，false=回退稠密 LU。
bool host_sparse_solve(GpuAhfLinearSolver *s, const double *d_values,
                       const double *d_rhs, double *d_solution) noexcept {
    static const bool force_lu = std::getenv("AMSS_AHF_DENSE_LU") != nullptr;
    if (force_lu) return false;
    const int n = s->n;
    if (cudaMemcpyAsync(s->h_values.data(), d_values, std::size_t(s->nnz) * sizeof(double),
                        cudaMemcpyDeviceToHost, s->stream) != cudaSuccess) return false;
    if (cudaMemcpyAsync(s->h_rhs.data(), d_rhs, std::size_t(n) * sizeof(double),
                        cudaMemcpyDeviceToHost, s->stream) != cudaSuccess) return false;
    if (cudaStreamSynchronize(s->stream) != cudaSuccess) return false;
    const int *rows = s->h_rows.data(), *columns = s->h_columns.data();
    double *a = s->h_values.data(), *b = s->h_rhs.data();
    double *x = s->h_sol.data(), *diag = s->h_diag.data();
    double *p = s->h_p.data(), *v = s->h_v.data(), *sv = s->h_s.data(), *t = s->h_t.data();
    double *r = s->h_r.data(), *rhat = s->h_rhat.data();
    for (int i = 0; i < n; ++i) diag[i] = 1.0;
    for (int i = 0; i < n; ++i)
        for (int k = rows[i]; k < rows[i + 1]; ++k) {
            if (!std::isfinite(a[k])) return false;
            if (columns[k] == i) diag[i] = a[k];
        }
    double bnorm2 = 0.0;
    for (int i = 0; i < n; ++i) {
        if (!std::isfinite(b[i])) return false;
        bnorm2 += b[i] * b[i];
        const double d = diag[i];
        if (!std::isfinite(d) || d == 0.0) return false;
        x[i] = b[i] / d;
    }
    const double bnorm = std::sqrt(bnorm2) > 0.0 ? std::sqrt(bnorm2) : 1.0;
    const double tol2 = 1e-24 * bnorm * bnorm;   // 相对残差 1e-12
    csr_ax(rows, columns, a, x, t, n);
    double rho_old = 0.0, alpha = 0.0, omega = 1.0;
    for (int i = 0; i < n; ++i) { r[i] = b[i] - t[i]; rhat[i] = r[i]; p[i] = r[i]; }
    bool converged = vdot(r, r, n) <= tol2;
    const int max_iterations = 300;
    for (int iteration = 0; !converged && iteration < max_iterations; ++iteration) {
        const double rho = vdot(rhat, r, n);
        if (!std::isfinite(rho) || std::fabs(rho) < 1e-300) return false;
        if (iteration > 0) {
            const double beta = (rho / rho_old) * (alpha / omega);
            if (!std::isfinite(beta)) return false;
            for (int i = 0; i < n; ++i) p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }
        csr_ax(rows, columns, a, p, v, n);
        const double denominator = vdot(rhat, v, n);
        if (!std::isfinite(denominator) || std::fabs(denominator) < 1e-300) return false;
        alpha = rho / denominator;
        if (!std::isfinite(alpha)) return false;
        for (int i = 0; i < n; ++i) sv[i] = r[i] - alpha * v[i];
        csr_ax(rows, columns, a, sv, t, n);
        const double tt = vdot(t, t, n);
        if (!std::isfinite(tt) || tt < 1e-300) return false;
        omega = vdot(t, sv, n) / tt;
        if (!std::isfinite(omega) || std::fabs(omega) < 1e-300) {
            for (int i = 0; i < n; ++i) x[i] += alpha * p[i];
            return false;   // 停滞：部分修正后回退稠密 LU
        }
        for (int i = 0; i < n; ++i) {
            x[i] += alpha * p[i] + omega * sv[i];
            r[i] = sv[i] - omega * t[i];
        }
        rho_old = rho;
        converged = vdot(r, r, n) <= tol2;
    }
    if (!converged) return false;
    for (int i = 0; i < n; ++i) if (!std::isfinite(x[i])) return false;
    if (cudaMemcpyAsync(d_solution, x, std::size_t(n) * sizeof(double),
                        cudaMemcpyHostToDevice, s->stream) != cudaSuccess) return false;
    return cudaStreamSynchronize(s->stream) == cudaSuccess;
}
} // namespace

int gpu_ahf_linear_solve(GpuAhfLinearSolver* s, const double* d_values,
                         const double* d_rhs, double* d_solution) noexcept {
    if (!s || !d_values || !d_rhs || !d_solution)
        return GPU_AHF_LINEAR_INVALID_ARGUMENT;
    if (overlaps(d_solution, s->vector_bytes, d_rhs, s->vector_bytes) ||
        overlaps(d_solution, s->vector_bytes, d_values, std::size_t(s->nnz) * sizeof(double)) ||
        overlaps(d_solution, s->vector_bytes, s->rows, (std::size_t(s->n) + 1) * sizeof(int)) ||
        overlaps(d_solution, s->vector_bytes, s->columns, std::size_t(s->nnz) * sizeof(int)))
        return GPU_AHF_LINEAR_INVALID_ARGUMENT;
    if (host_sparse_solve(s, d_values, d_rhs, d_solution))
        return GPU_AHF_LINEAR_SUCCESS;

#define AHF_CUDA_RETURN(call) do { if ((call) != cudaSuccess) return GPU_AHF_LINEAR_CUDA_FAILURE; } while (0)
    AHF_CUDA_RETURN(cudaMemsetAsync(s->info, 0, sizeof(int), s->stream));
    check_finite<<<finite_blocks(s->nnz > s->n ? s->nnz : s->n), 256, 0, s->stream>>>(
        s->nnz, d_values, s->n, d_rhs, s->info);
    AHF_CUDA_RETURN(cudaGetLastError());
    int info = 0;
    AHF_CUDA_RETURN(read_scalar(s->stream, s->info, info));
    if (info) return GPU_AHF_LINEAR_NONFINITE_INPUT;

    AHF_CUDA_RETURN(cudaMemsetAsync(s->dense, 0, s->dense_bytes, s->stream));
    csr_to_dense<<<(s->n - 1) / 8 + 1, 256, 0, s->stream>>>(
        s->n, s->rows, s->columns, d_values, s->dense);
    AHF_CUDA_RETURN(cudaGetLastError());
    AHF_CUDA_RETURN(cudaMemcpyAsync(d_solution, d_rhs, s->vector_bytes,
                                    cudaMemcpyDeviceToDevice, s->stream));
    if (cusolverDnDgetrf(s->handle, s->n, s->n, s->dense, s->n,
                         s->workspace, s->pivots, s->info) != CUSOLVER_STATUS_SUCCESS)
        return GPU_AHF_LINEAR_GETRF_FAILURE;
    AHF_CUDA_RETURN(read_scalar(s->stream, s->info, info));
    if (info < 0 || info > s->n) return GPU_AHF_LINEAR_GETRF_ARGUMENT;
    if (info > 0) return GPU_AHF_LINEAR_SINGULAR_BASE - info;

    if (cusolverDnDgetrs(s->handle, CUBLAS_OP_N, s->n, 1, s->dense, s->n,
                         s->pivots, d_solution, s->n, s->info) != CUSOLVER_STATUS_SUCCESS)
        return GPU_AHF_LINEAR_GETRS_FAILURE;
    AHF_CUDA_RETURN(read_scalar(s->stream, s->info, info));
    if (info != 0) return GPU_AHF_LINEAR_GETRS_INFO;

    AHF_CUDA_RETURN(cudaMemsetAsync(s->info, 0, sizeof(int), s->stream));
    check_finite<<<finite_blocks(s->n), 256, 0, s->stream>>>(s->n, d_solution, 0, nullptr, s->info);
    AHF_CUDA_RETURN(cudaGetLastError());
    AHF_CUDA_RETURN(read_scalar(s->stream, s->info, info));
    return info ? GPU_AHF_LINEAR_NONFINITE_OUTPUT : GPU_AHF_LINEAR_SUCCESS;
#undef AHF_CUDA_RETURN
}
