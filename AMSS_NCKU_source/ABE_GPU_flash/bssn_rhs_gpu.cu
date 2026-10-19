#include "bssn_rhs.h"
#include "bssn_rhs_gpu.cuh"

#include "fmisc_gpu.cuh"
#include "diff_new_gpu.cuh"
#include "kodiss_gpu.cuh"
#include "lopsidediff_gpu.cuh"
#include "gpu_runtime_coordinator.h"

#include <cuda_runtime.h>
#include <math.h>
#include <device_launch_parameters.h>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <vector>

constexpr double SYM = 1.0;
constexpr double ANTI = -1.0;

constexpr double ONE = 1.0;
constexpr double TWO = 2.0;
constexpr double FOUR_PI  = 4.0 * M_PI;
constexpr double EIGHT_PI = 8.0 * M_PI;
constexpr double F16_PI   = 16.0 * M_PI;

constexpr double F1o3 = 1.0 / 3.0;
constexpr double F2o3 = 2.0 / 3.0;
constexpr double F3o2 = 1.5;
constexpr double HALF = 0.5;
constexpr double FF = 0.75;
constexpr double eta = 2.0;
constexpr double Z4C_KAPPA1 = 2.0e-2;
constexpr double Z4C_KAPPA2 = 0.0;

#if (GAUGE == 6 || GAUGE == 7)
extern int BHN;
extern double Mass[3];
extern double PBH[9];

__constant__ int d_bssn_gauge_bh_count;
__constant__ double d_bssn_gauge_mass[3];
__constant__ double d_bssn_gauge_pbh[9];
#endif

namespace {

constexpr int kBssnBlockX = 8;
constexpr int kBssnBlockY = 8;
constexpr int kBssnBlockZ = 4;

#if ghost_width == 5
#define BSSN_CENTERED_SMEM_TYPED(type, name) __shared__ type name[12][16][17]
#define BSSN_TRANSPORT_SMEM_TYPED(type, name) __shared__ type name[14][18][19]
#elif ghost_width == 4
#define BSSN_CENTERED_SMEM_TYPED(type, name) __shared__ type name[10][14][15]
#define BSSN_TRANSPORT_SMEM_TYPED(type, name) __shared__ type name[14][18][19]
#else
#define BSSN_CENTERED_SMEM_TYPED(type, name) __shared__ type name[8][12][13]
#define BSSN_TRANSPORT_SMEM_TYPED(type, name) __shared__ type name[10][14][15]
#endif
#define BSSN_CENTERED_SMEM(name) BSSN_CENTERED_SMEM_TYPED(double, name)
#define BSSN_TRANSPORT_SMEM(name) BSSN_TRANSPORT_SMEM_TYPED(double, name)

struct StencilBounds {
    int imin;
    int jmin;
    int kmin;
    int imax;
    int jmax;
    int kmax;
};

__device__ __forceinline__ int flatten_bssn_index(int i, int j, int k, int ex0, int ex1) {
    return i + ex0 * (j + ex1 * k);
}

__device__ __forceinline__ StencilBounds make_stencil_bounds(
    int ex0,
    int ex1,
    int ex2,
    int symmetry,
    const double* X,
    const double* Y,
    const double* Z,
    double dX,
    double dY,
    double dZ,
    int boundary_radius
) {
    StencilBounds bounds{0, 0, 0, ex0 - 1, ex1 - 1, ex2 - 1};

    if (symmetry > 0 && fabs(Z[0]) < dZ) bounds.kmin = -boundary_radius;
    if (symmetry > 1 && fabs(X[0]) < dX) bounds.imin = -boundary_radius;
    if (symmetry > 1 && fabs(Y[0]) < dY) bounds.jmin = -boundary_radius;

    return bounds;
}

#if ghost_width == 4 || ghost_width == 5
template <typename Real>
__device__ __forceinline__ void accumulate_fd6_transport_axis(
    Real velocity,
    Real d60,
    Real d12,
    Real d2,
    int center,
    int lower_bound,
    int upper_bound,
    Real* W
) {
    const Real zero = Real(0.0);
    const Real v60 = velocity * d60;

    if (velocity > zero) {
        if (center + 4 <= upper_bound && center - 2 >= lower_bound) {
            W[3] += Real(2.0) * v60;   W[4] -= Real(24.0) * v60; W[5] -= Real(35.0) * v60;
            W[6] += Real(80.0) * v60;  W[7] -= Real(30.0) * v60; W[8] += Real(8.0) * v60;
            W[9] -= v60;
            return;
        }
        if (center + 5 <= upper_bound && center - 1 >= lower_bound) {
            W[4] -= Real(10.0) * v60;  W[5] -= Real(77.0) * v60; W[6] += Real(150.0) * v60;
            W[7] -= Real(100.0) * v60; W[8] += Real(50.0) * v60; W[9] -= Real(15.0) * v60;
            W[10] += Real(2.0) * v60;
            return;
        }
    } else if (velocity < zero) {
        const Real n60 = -v60;
        if (center - 4 >= lower_bound && center + 2 <= upper_bound) {
            W[7] += Real(2.0) * n60;   W[6] -= Real(24.0) * n60; W[5] -= Real(35.0) * n60;
            W[4] += Real(80.0) * n60;  W[3] -= Real(30.0) * n60; W[2] += Real(8.0) * n60;
            W[1] -= n60;
            return;
        }
        if (center - 5 >= lower_bound && center + 1 <= upper_bound) {
            W[6] -= Real(10.0) * n60;  W[5] -= Real(77.0) * n60; W[4] += Real(150.0) * n60;
            W[3] -= Real(100.0) * n60; W[2] += Real(50.0) * n60; W[1] -= Real(15.0) * n60;
            W[0] += Real(2.0) * n60;
            return;
        }
    }

    if (center + 3 <= upper_bound && center - 3 >= lower_bound) {
        W[2] -= v60;           W[3] += Real(9.0) * v60;  W[4] -= Real(45.0) * v60;
        W[6] += Real(45.0) * v60;    W[7] -= Real(9.0) * v60;  W[8] += v60;
    } else if (center + 2 <= upper_bound && center - 2 >= lower_bound) {
        const Real v12 = velocity * d12;
        W[3] += v12;           W[4] -= Real(8.0) * v12;  W[6] += Real(8.0) * v12;
        W[7] -= v12;
    } else if (center + 1 <= upper_bound && center - 1 >= lower_bound) {
        const Real v2 = velocity * d2;
        W[4] -= v2;            W[6] += v2;
    }
}
#endif

#if ghost_width == 5
template <typename Real>
__device__ __forceinline__ void accumulate_fd8_transport_axis(
    Real velocity,
    Real d840,
    Real d60,
    Real d12,
    Real d2,
    int center,
    int lower_bound,
    int upper_bound,
    Real* W
) {
    const Real zero = Real(0.0);
    const Real v840 = velocity * d840;

    if (velocity > zero) {
        if (center + 5 <= upper_bound && center - 3 >= lower_bound) {
            W[2] -= Real(5.0) * v840;    W[3] += Real(60.0) * v840;   W[4] -= Real(420.0) * v840;
            W[5] -= Real(378.0) * v840;  W[6] += Real(1050.0) * v840; W[7] -= Real(420.0) * v840;
            W[8] += Real(140.0) * v840;  W[9] -= Real(30.0) * v840;   W[10] += Real(3.0) * v840;
            return;
        }
    } else if (velocity < zero) {
        const Real n840 = -v840;
        if (center - 5 >= lower_bound && center + 3 <= upper_bound) {
            W[0] += Real(3.0) * n840;    W[1] -= Real(30.0) * n840;   W[2] += Real(140.0) * n840;
            W[3] -= Real(420.0) * n840;  W[4] += Real(1050.0) * n840; W[5] -= Real(378.0) * n840;
            W[6] -= Real(420.0) * n840;  W[7] += Real(60.0) * n840;   W[8] -= Real(5.0) * n840;
            return;
        }
    }

    if (center + 4 <= upper_bound && center - 4 >= lower_bound) {
        W[1] += Real(3.0) * v840;      W[2] -= Real(32.0) * v840;   W[3] += Real(168.0) * v840;
        W[4] -= Real(672.0) * v840;    W[6] += Real(672.0) * v840;  W[7] -= Real(168.0) * v840;
        W[8] += Real(32.0) * v840;     W[9] -= Real(3.0) * v840;
    } else if (center + 3 <= upper_bound && center - 3 >= lower_bound) {
        const Real v60 = velocity * d60;
        W[2] -= v60;             W[3] += Real(9.0) * v60;     W[4] -= Real(45.0) * v60;
        W[6] += Real(45.0) * v60;      W[7] -= Real(9.0) * v60;     W[8] += v60;
    } else if (center + 2 <= upper_bound && center - 2 >= lower_bound) {
        const Real v12 = velocity * d12;
        W[3] += v12;             W[4] -= Real(8.0) * v12;     W[6] += Real(8.0) * v12;
        W[7] -= v12;
    } else if (center + 1 <= upper_bound && center - 1 >= lower_bound) {
        const Real v2 = velocity * d2;
        W[4] -= v2;              W[6] += v2;
    }
}
#endif

inline dim3 bssn_block_shape() {
    return dim3(kBssnBlockX, kBssnBlockY, kBssnBlockZ);
}

inline dim3 bssn_grid_shape(const int* ex, const dim3& block) {
    return dim3(
        (ex[0] + block.x - 1) / block.x,
        (ex[1] + block.y - 1) / block.y,
        (ex[2] + block.z - 1) / block.z
    );
}

struct RhsTraceField {
    const char* name;
    const double* ptr;
};

int rhs_trace_env_int(const char* name, int fallback) {
    const char* value = getenv(name);
    return value ? atoi(value) : fallback;
}

bool bssn_env_value_enabled(const char* name) {
    const char* value = getenv(name);
    if (!value) return false;

    return std::strcmp(value, "0") != 0
        && std::strcmp(value, "false") != 0
        && std::strcmp(value, "FALSE") != 0
        && std::strcmp(value, "off") != 0
        && std::strcmp(value, "OFF") != 0
        && std::strcmp(value, "no") != 0
        && std::strcmp(value, "NO") != 0;
}

bool bssn_fp32_rhs_enabled() {
#if ABEtype == 0
    return bssn_env_value_enabled("AMSS_NCKU_ALLOW_EXPERIMENTAL_FP32_RHS")
        && bssn_env_value_enabled("AMSS_NCKU_GPU_FP32_RHS");
#else
    return false;
#endif
}

void report_bssn_fp32_rhs_once() {
    static bool reported = false;
    if (reported) return;
    reported = true;

    std::cout << " BSSN GPU RHS: using FP32 float-shadow RHS kernels"
              << " (AMSS_NCKU_ALLOW_EXPERIMENTAL_FP32_RHS=1,"
              << " AMSS_NCKU_GPU_FP32_RHS=1)" << std::endl;
}

__global__ void bssn_double_to_float_fields_kernel(
    const double* const* __restrict__ sources,
    float* const* __restrict__ targets,
    int field_count,
    int point_count
) {
    int point = blockIdx.x * blockDim.x + threadIdx.x;
    int field = blockIdx.y;
    if (field >= field_count || point >= point_count) return;

    targets[field][point] = static_cast<float>(sources[field][point]);
}

__global__ void bssn_float_to_double_fields_kernel(
    const float* const* __restrict__ sources,
    double* const* __restrict__ targets,
    int field_count,
    int point_count
) {
    int point = blockIdx.x * blockDim.x + threadIdx.x;
    int field = blockIdx.y;
    if (field >= field_count || point >= point_count) return;

    targets[field][point] = static_cast<double>(sources[field][point]);
}

__global__ void bssn_zero_float_fields_kernel(
    float* const* __restrict__ targets,
    int field_count,
    int point_count
) {
    int point = blockIdx.x * blockDim.x + threadIdx.x;
    int field = blockIdx.y;
    if (field >= field_count || point >= point_count) return;

    targets[field][point] = 0.0f;
}

struct BssnRhsFp32Scratch {
    size_t element_capacity = 0;
    int pointer_capacity = 0;
    std::vector<float*> buffers;
    const double** d_double_sources = nullptr;
    double** d_double_targets = nullptr;
    float** d_float_sources = nullptr;
    float** d_float_targets = nullptr;

    template <typename Pointer>
    static void allocate_pointer_array(Pointer** target, int count) {
        void* raw = nullptr;
        CUDA_CHECK(cudaMalloc(&raw, count * sizeof(Pointer)));
        *target = static_cast<Pointer*>(raw);
    }

    void reserve(size_t element_count, int buffer_count, int pointer_count) {
        if (static_cast<int>(buffers.size()) < buffer_count) {
            const int old_count = static_cast<int>(buffers.size());
            const size_t allocation_count = element_capacity > element_count ? element_capacity : element_count;
            buffers.resize(buffer_count, nullptr);
            for (int field = old_count; field < buffer_count; ++field) {
                void* raw = nullptr;
                CUDA_CHECK(cudaMalloc(&raw, allocation_count * sizeof(float)));
                buffers[field] = static_cast<float*>(raw);
            }
        }

        if (element_capacity < element_count) {
            for (float*& buffer : buffers) {
                if (buffer) CUDA_CHECK(cudaFree(buffer));
                void* raw = nullptr;
                CUDA_CHECK(cudaMalloc(&raw, element_count * sizeof(float)));
                buffer = static_cast<float*>(raw);
            }
            element_capacity = element_count;
        } else if (element_capacity == 0) {
            element_capacity = element_count;
        }

        if (pointer_capacity < pointer_count) {
            if (d_double_sources) CUDA_CHECK(cudaFree(d_double_sources));
            if (d_double_targets) CUDA_CHECK(cudaFree(d_double_targets));
            if (d_float_sources) CUDA_CHECK(cudaFree(d_float_sources));
            if (d_float_targets) CUDA_CHECK(cudaFree(d_float_targets));

            allocate_pointer_array(&d_double_sources, pointer_count);
            allocate_pointer_array(&d_double_targets, pointer_count);
            allocate_pointer_array(&d_float_sources, pointer_count);
            allocate_pointer_array(&d_float_targets, pointer_count);
            pointer_capacity = pointer_count;
        }
    }
};

BssnRhsFp32Scratch& bssn_rhs_fp32_scratch() {
    static thread_local BssnRhsFp32Scratch scratch;
    return scratch;
}

void bssn_convert_double_to_float_fields(
    cudaStream_t stream,
    const std::vector<double*>& sources,
    const std::vector<float*>& targets,
    int point_count,
    BssnRhsFp32Scratch& scratch
) {
    if (sources.empty()) return;
    CUDA_CHECK(cudaMemcpyAsync(
        scratch.d_double_sources,
        sources.data(),
        sources.size() * sizeof(double*),
        cudaMemcpyHostToDevice,
        stream
    ));
    CUDA_CHECK(cudaMemcpyAsync(
        scratch.d_float_targets,
        targets.data(),
        targets.size() * sizeof(float*),
        cudaMemcpyHostToDevice,
        stream
    ));

    dim3 block(256);
    dim3 grid((point_count + block.x - 1) / block.x, static_cast<unsigned int>(sources.size()));
    bssn_double_to_float_fields_kernel<<<grid, block, 0, stream>>>(
        scratch.d_double_sources,
        scratch.d_float_targets,
        static_cast<int>(sources.size()),
        point_count
    );
}

void bssn_convert_float_to_double_fields(
    cudaStream_t stream,
    const std::vector<float*>& sources,
    const std::vector<double*>& targets,
    int point_count,
    BssnRhsFp32Scratch& scratch
) {
    if (sources.empty()) return;
    CUDA_CHECK(cudaMemcpyAsync(
        scratch.d_float_sources,
        sources.data(),
        sources.size() * sizeof(float*),
        cudaMemcpyHostToDevice,
        stream
    ));
    CUDA_CHECK(cudaMemcpyAsync(
        scratch.d_double_targets,
        targets.data(),
        targets.size() * sizeof(double*),
        cudaMemcpyHostToDevice,
        stream
    ));

    dim3 block(256);
    dim3 grid((point_count + block.x - 1) / block.x, static_cast<unsigned int>(sources.size()));
    bssn_float_to_double_fields_kernel<<<grid, block, 0, stream>>>(
        scratch.d_float_sources,
        scratch.d_double_targets,
        static_cast<int>(sources.size()),
        point_count
    );
}

void bssn_zero_float_fields(
    cudaStream_t stream,
    const std::vector<float*>& targets,
    int point_count,
    BssnRhsFp32Scratch& scratch
) {
    if (targets.empty()) return;
    CUDA_CHECK(cudaMemcpyAsync(
        scratch.d_float_targets,
        targets.data(),
        targets.size() * sizeof(float*),
        cudaMemcpyHostToDevice,
        stream
    ));

    dim3 block(256);
    dim3 grid((point_count + block.x - 1) / block.x, static_cast<unsigned int>(targets.size()));
    bssn_zero_float_fields_kernel<<<grid, block, 0, stream>>>(
        scratch.d_float_targets,
        static_cast<int>(targets.size()),
        point_count
    );
}

template<typename Real, int SMEM_Z, int SMEM_Y, int SMEM_X>
__device__ __forceinline__ void load_field_to_smem_rad3_mixed(
    Real (&smem)[SMEM_Z][SMEM_Y][SMEM_X], const Real* f,
    const int ex[3], int block_i, int block_j, int block_k,
    double SYM1, double SYM2, double SYM3
) {
    int tid = threadIdx.z * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x;
    double SoA[3] = {SYM1, SYM2, SYM3};
    const int tile_x = blockDim.x + 2 * kBssnTransportLoadRadius;
    const int tile_y = blockDim.y + 2 * kBssnTransportLoadRadius;
    const int tile_z = blockDim.z + 2 * kBssnTransportLoadRadius;
    const int tile_xy = tile_x * tile_y;
    const int tile_size = tile_xy * tile_z;

    for (int idx = tid; idx < tile_size; idx += 256) {
        int loc_k = idx / tile_xy;
        int rem = idx % tile_xy;
        int loc_j = rem / tile_x;
        int loc_i = rem % tile_x;

        int glob_i = block_i + loc_i - kBssnTransportLoadRadius;
        int glob_j = block_j + loc_j - kBssnTransportLoadRadius;
        int glob_k = block_k + loc_k - kBssnTransportLoadRadius;

        smem[loc_k][loc_j][loc_i] = static_cast<Real>(d_symmetry_bd_0b_t(
            kBssnTransportBoundaryRadius,
            ex[0], ex[1], ex[2], f, glob_i, glob_j, glob_k, SoA[0], SoA[1], SoA[2]
        ));
    }
}

bool rhs_nonfinite_trace_in_scope(int lev, int co) {
    if (!getenv("AMSS_TRAJECTORY_PREFIX_TRACE_RHS_NONFINITE"))
        return false;
    const int focus_lev = rhs_trace_env_int("AMSS_TRAJECTORY_PREFIX_TRACE_RHS_NONFINITE_LEVEL", -1);
    if (focus_lev >= 0 && focus_lev != lev)
        return false;
    const int focus_co = rhs_trace_env_int("AMSS_TRAJECTORY_PREFIX_TRACE_RHS_NONFINITE_CO", -1);
    if (focus_co >= 0 && focus_co != co)
        return false;
    return true;
}

void trace_rhs_nonfinite_summary(
    cudaStream_t stream,
    const char* phase,
    const int* ex,
    const double* X,
    const double* Y,
    const double* Z,
    int lev,
    int co,
    const RhsTraceField* fields,
    int field_count
) {
    if (!rhs_nonfinite_trace_in_scope(lev, co))
        return;

    const int nx = ex[0];
    const int ny = ex[1];
    const int nz = ex[2];
    const int count = nx * ny * nz;
    std::vector<double> host(count);

    int total_nonfinite = 0;
    int variables_with_nonfinite = 0;
    const char* first_var = "none";
    int first_field_index = -1;
    int first_i = -1;
    int first_j = -1;
    int first_k = -1;
    double first_value = 0.0;

    for (int field_index = 0; field_index < field_count; ++field_index) {
        if (!fields[field_index].ptr)
            continue;
        CUDA_CHECK(cudaMemcpyAsync(
            host.data(),
            fields[field_index].ptr,
            static_cast<size_t>(count) * sizeof(double),
            cudaMemcpyDeviceToHost,
            stream));
        CUDA_CHECK(cudaStreamSynchronize(stream));

        int field_nonfinite = 0;
        for (int idx = 0; idx < count; ++idx) {
            const double value = host[idx];
            if (std::isfinite(value))
                continue;
            ++field_nonfinite;
            if (first_field_index < 0) {
                first_field_index = field_index;
                first_var = fields[field_index].name;
                first_i = idx % nx;
                first_j = (idx / nx) % ny;
                first_k = idx / (nx * ny);
                first_value = value;
            }
        }
        if (field_nonfinite > 0) {
            ++variables_with_nonfinite;
            total_nonfinite += field_nonfinite;
        }
    }

    double first_x = 0.0;
    double first_y = 0.0;
    double first_z = 0.0;
    if (first_field_index >= 0) {
        CUDA_CHECK(cudaMemcpyAsync(&first_x, X + first_i, sizeof(double), cudaMemcpyDeviceToHost, stream));
        CUDA_CHECK(cudaMemcpyAsync(&first_y, Y + first_j, sizeof(double), cudaMemcpyDeviceToHost, stream));
        CUDA_CHECK(cudaMemcpyAsync(&first_z, Z + first_k, sizeof(double), cudaMemcpyDeviceToHost, stream));
        CUDA_CHECK(cudaStreamSynchronize(stream));
    }

    std::cout << " rhs_nonfinite"
              << " phase=" << phase
              << " lev=" << lev
              << " co=" << co
              << " shape=(" << nx << "," << ny << "," << nz << ")"
              << " total_nonfinite=" << total_nonfinite
              << " variables_with_nonfinite=" << variables_with_nonfinite
              << " first_var=" << first_var
              << " first_ijk=(" << first_i << "," << first_j << "," << first_k << ")"
              << " first_xyz=(" << first_x << "," << first_y << "," << first_z << ")"
              << " first_value=" << first_value
              << std::endl;
}

void trace_metric_det_summary(
    cudaStream_t stream,
    const char* phase,
    const int* ex,
    const double* X,
    const double* Y,
    const double* Z,
    int lev,
    int co,
    const double* dxx,
    const double* gxy,
    const double* gxz,
    const double* dyy,
    const double* gyz,
    const double* dzz
) {
    if (!rhs_nonfinite_trace_in_scope(lev, co))
        return;

    const int nx = ex[0];
    const int ny = ex[1];
    const int nz = ex[2];
    const int count = nx * ny * nz;
    const double det_eps = getenv("AMSS_TRAJECTORY_PREFIX_TRACE_RHS_DET_EPS")
        ? atof(getenv("AMSS_TRAJECTORY_PREFIX_TRACE_RHS_DET_EPS"))
        : 1.0e-300;

    std::vector<double> h_dxx(count), h_gxy(count), h_gxz(count);
    std::vector<double> h_dyy(count), h_gyz(count), h_dzz(count);
    CUDA_CHECK(cudaMemcpyAsync(h_dxx.data(), dxx, static_cast<size_t>(count) * sizeof(double), cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaMemcpyAsync(h_gxy.data(), gxy, static_cast<size_t>(count) * sizeof(double), cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaMemcpyAsync(h_gxz.data(), gxz, static_cast<size_t>(count) * sizeof(double), cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaMemcpyAsync(h_dyy.data(), dyy, static_cast<size_t>(count) * sizeof(double), cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaMemcpyAsync(h_gyz.data(), gyz, static_cast<size_t>(count) * sizeof(double), cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaMemcpyAsync(h_dzz.data(), dzz, static_cast<size_t>(count) * sizeof(double), cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaStreamSynchronize(stream));

    int bad_det_count = 0;
    int nonfinite_input_count = 0;
    int first_i = -1;
    int first_j = -1;
    int first_k = -1;
    double first_det = 0.0;
    double first_gxx = 0.0;
    double first_gxy = 0.0;
    double first_gxz = 0.0;
    double first_gyy = 0.0;
    double first_gyz = 0.0;
    double first_gzz = 0.0;
    double min_abs_det = 0.0;
    int min_i = -1;
    int min_j = -1;
    int min_k = -1;

    for (int idx = 0; idx < count; ++idx) {
        const double l_gxx = h_dxx[idx] + ONE;
        const double l_gxy = h_gxy[idx];
        const double l_gxz = h_gxz[idx];
        const double l_gyy = h_dyy[idx] + ONE;
        const double l_gyz = h_gyz[idx];
        const double l_gzz = h_dzz[idx] + ONE;
        const bool inputs_finite = std::isfinite(l_gxx) && std::isfinite(l_gxy) &&
            std::isfinite(l_gxz) && std::isfinite(l_gyy) && std::isfinite(l_gyz) &&
            std::isfinite(l_gzz);
        if (!inputs_finite)
            ++nonfinite_input_count;
        const double det = l_gxx * l_gyy * l_gzz + l_gxy * l_gyz * l_gxz
                         + l_gxz * l_gxy * l_gyz
                         - l_gxz * l_gyy * l_gxz
                         - l_gxy * l_gxy * l_gzz
                         - l_gxx * l_gyz * l_gyz;
        const double abs_det = fabs(det);
        if (min_i < 0 || abs_det < min_abs_det) {
            min_abs_det = abs_det;
            min_i = idx % nx;
            min_j = (idx / nx) % ny;
            min_k = idx / (nx * ny);
        }
        if (!std::isfinite(det) || abs_det <= det_eps) {
            ++bad_det_count;
            if (first_i < 0) {
                first_i = idx % nx;
                first_j = (idx / nx) % ny;
                first_k = idx / (nx * ny);
                first_det = det;
                first_gxx = l_gxx;
                first_gxy = l_gxy;
                first_gxz = l_gxz;
                first_gyy = l_gyy;
                first_gyz = l_gyz;
                first_gzz = l_gzz;
            }
        }
    }

    double first_x = 0.0;
    double first_y = 0.0;
    double first_z = 0.0;
    if (first_i >= 0) {
        CUDA_CHECK(cudaMemcpyAsync(&first_x, X + first_i, sizeof(double), cudaMemcpyDeviceToHost, stream));
        CUDA_CHECK(cudaMemcpyAsync(&first_y, Y + first_j, sizeof(double), cudaMemcpyDeviceToHost, stream));
        CUDA_CHECK(cudaMemcpyAsync(&first_z, Z + first_k, sizeof(double), cudaMemcpyDeviceToHost, stream));
        CUDA_CHECK(cudaStreamSynchronize(stream));
    }

    std::cout << " rhs_metric_det"
              << " phase=" << phase
              << " lev=" << lev
              << " co=" << co
              << " shape=(" << nx << "," << ny << "," << nz << ")"
              << " det_eps=" << det_eps
              << " bad_det_count=" << bad_det_count
              << " nonfinite_input_count=" << nonfinite_input_count
              << " min_abs_det=" << min_abs_det
              << " min_ijk=(" << min_i << "," << min_j << "," << min_k << ")"
              << " first_ijk=(" << first_i << "," << first_j << "," << first_k << ")"
              << " first_xyz=(" << first_x << "," << first_y << "," << first_z << ")"
              << " first_det=" << first_det
              << " first_metric=(" << first_gxx << "," << first_gxy << "," << first_gxz
              << "," << first_gyy << "," << first_gyz << "," << first_gzz << ")"
              << std::endl;
}

template <typename Real>
__device__ __forceinline__ void compute_inverse_metric_t(
    Real gxx,
    Real gxy,
    Real gxz,
    Real gyy,
    Real gyz,
    Real gzz,
    Real* gupxx,
    Real* gupxy,
    Real* gupxz,
    Real* gupyy,
    Real* gupyz,
    Real* gupzz
) {
    Real det = gxx * gyy * gzz + gxy * gyz * gxz
               + gxz * gxy * gyz
               - gxz * gyy * gxz
               - gxy * gxy * gzz
               - gxx * gyz * gyz;
    Real inv_det = Real(1.0) / det;

    *gupxx = (gyy * gzz - gyz * gyz) * inv_det;
    *gupxy = (gyz * gxz - gxy * gzz) * inv_det;
    *gupxz = (gxy * gyz - gyy * gxz) * inv_det;
    *gupyy = (gxx * gzz - gxz * gxz) * inv_det;
    *gupyz = (gxy * gxz - gxx * gyz) * inv_det;
    *gupzz = (gxx * gyy - gxy * gxy) * inv_det;
}

__device__ __forceinline__ void compute_inverse_metric(
    double gxx,
    double gxy,
    double gxz,
    double gyy,
    double gyz,
    double gzz,
    double* gupxx,
    double* gupxy,
    double* gupxz,
    double* gupyy,
    double* gupyz,
    double* gupzz
) {
    compute_inverse_metric_t<double>(
        gxx, gxy, gxz, gyy, gyz, gzz,
        gupxx, gupxy, gupxz, gupyy, gupyz, gupzz
    );
}

template <typename Real>
__device__ __forceinline__ void compute_raised_extrinsic_curvature_t(
    Real Axx,
    Real Axy,
    Real Axz,
    Real Ayy,
    Real Ayz,
    Real Azz,
    Real gupxx,
    Real gupxy,
    Real gupxz,
    Real gupyy,
    Real gupyz,
    Real gupzz,
    Real* Au_d_xx,
    Real* Au_d_xy,
    Real* Au_d_xz,
    Real* Au_d_yx,
    Real* Au_d_yy,
    Real* Au_d_yz,
    Real* Au_d_zx,
    Real* Au_d_zy,
    Real* Au_d_zz,
    Real* Aup_xx,
    Real* Aup_xy,
    Real* Aup_xz,
    Real* Aup_yy,
    Real* Aup_yz,
    Real* Aup_zz
) {
    *Au_d_xx = gupxx * Axx + gupxy * Axy + gupxz * Axz;
    *Au_d_xy = gupxx * Axy + gupxy * Ayy + gupxz * Ayz;
    *Au_d_xz = gupxx * Axz + gupxy * Ayz + gupxz * Azz;
    *Au_d_yx = gupxy * Axx + gupyy * Axy + gupyz * Axz;
    *Au_d_yy = gupxy * Axy + gupyy * Ayy + gupyz * Ayz;
    *Au_d_yz = gupxy * Axz + gupyy * Ayz + gupyz * Azz;
    *Au_d_zx = gupxz * Axx + gupyz * Axy + gupzz * Axz;
    *Au_d_zy = gupxz * Axy + gupyz * Ayy + gupzz * Ayz;
    *Au_d_zz = gupxz * Axz + gupyz * Ayz + gupzz * Azz;

    *Aup_xx = *Au_d_xx * gupxx + *Au_d_xy * gupxy + *Au_d_xz * gupxz;
    *Aup_xy = *Au_d_xx * gupxy + *Au_d_xy * gupyy + *Au_d_xz * gupyz;
    *Aup_xz = *Au_d_xx * gupxz + *Au_d_xy * gupyz + *Au_d_xz * gupzz;
    *Aup_yy = *Au_d_yx * gupxy + *Au_d_yy * gupyy + *Au_d_yz * gupyz;
    *Aup_yz = *Au_d_yx * gupxz + *Au_d_yy * gupyz + *Au_d_yz * gupzz;
    *Aup_zz = *Au_d_zx * gupxz + *Au_d_zy * gupyz + *Au_d_zz * gupzz;
}

__device__ __forceinline__ void compute_raised_extrinsic_curvature(
    double Axx,
    double Axy,
    double Axz,
    double Ayy,
    double Ayz,
    double Azz,
    double gupxx,
    double gupxy,
    double gupxz,
    double gupyy,
    double gupyz,
    double gupzz,
    double* Au_d_xx,
    double* Au_d_xy,
    double* Au_d_xz,
    double* Au_d_yx,
    double* Au_d_yy,
    double* Au_d_yz,
    double* Au_d_zx,
    double* Au_d_zy,
    double* Au_d_zz,
    double* Aup_xx,
    double* Aup_xy,
    double* Aup_xz,
    double* Aup_yy,
    double* Aup_yz,
    double* Aup_zz
) {
    compute_raised_extrinsic_curvature_t<double>(
        Axx, Axy, Axz, Ayy, Ayz, Azz,
        gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
        Au_d_xx, Au_d_xy, Au_d_xz, Au_d_yx, Au_d_yy, Au_d_yz, Au_d_zx, Au_d_zy, Au_d_zz,
        Aup_xx, Aup_xy, Aup_xz, Aup_yy, Aup_yz, Aup_zz
    );
}

__device__ __forceinline__ void compute_escalar_frpotential(
    double Sphi,
    double a2,
    double* V,
    double* dVdSphi
) {
#if (EScalar_CC == 1 || EScalar_CC == 4 || EScalar_CC == 5)
    *V = 0.0;
    *dVdSphi = 0.0;
#elif (EScalar_CC == 2 || EScalar_CC == 3)
    const double fac = sqrt(M_PI / 3.0);
    const double e1 = exp(-8.0 * fac * Sphi);
    const double e2 = exp(4.0 * fac * Sphi);
    *V = e1 * (1.0 - e2) * (1.0 - e2) / (32.0 * M_PI * a2);
    *dVdSphi = (1.0 / a2 / 12.0) * sqrt(3.0 / M_PI) * e1 * (-1.0 + e2);
#else
#error Unsupported EScalar_CC in compute_escalar_frpotential.
#endif
}

}  // namespace

__launch_bounds__(256, 2)
__global__ void evolve_bssn_escalar_matter_kernel(
    const int ex0, const int ex1, const int ex2,
    const int symmetry, const int lev,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    const double* __restrict__ chi, const double* __restrict__ trK,
    const double* __restrict__ dxx, const double* __restrict__ gxy, const double* __restrict__ gxz,
    const double* __restrict__ dyy, const double* __restrict__ gyz, const double* __restrict__ dzz,
    const double* __restrict__ Gamx, const double* __restrict__ Gamy, const double* __restrict__ Gamz,
    const double* __restrict__ Lap,
    const double* __restrict__ betax, const double* __restrict__ betay, const double* __restrict__ betaz,
    const double* __restrict__ Sphi, const double* __restrict__ Spi,
    double* __restrict__ Sphi_rhs, double* __restrict__ Spi_rhs,
    double* __restrict__ rho, double* __restrict__ Sx, double* __restrict__ Sy, double* __restrict__ Sz,
    double* __restrict__ Sxx, double* __restrict__ Sxy, double* __restrict__ Sxz,
    double* __restrict__ Syy, double* __restrict__ Syz, double* __restrict__ Szz,
    double a2,
    double eps
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    int dims[3] = {ex0, ex1, ex2};
    double dX = X[1] - X[0];
    double dY = Y[1] - Y[0];
    double dZ = Z[1] - Z[0];
    StencilBounds bounds = make_stencil_bounds(ex0, ex1, ex2, symmetry, X, Y, Z, dX, dY, dZ, kBssnCenteredBoundaryRadius);

    int block_i = blockIdx.x * blockDim.x;
    int block_j = blockIdx.y * blockDim.y;
    int block_k = blockIdx.z * blockDim.z;

    BSSN_CENTERED_SMEM(smem);

    double chix, chiy, chiz;
    load_field_to_smem(smem, chi, dims, block_i, block_j, block_k, SYM, SYM, SYM);
    __syncthreads();
    compute_first_derivs_smem(smem, &chix, &chiy, &chiz, dX, dY, dZ, i, j, k,
                              bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax);
    __syncthreads();

    double Lapx, Lapy, Lapz;
    load_field_to_smem(smem, Lap, dims, block_i, block_j, block_k, SYM, SYM, SYM);
    __syncthreads();
    compute_first_derivs_smem(smem, &Lapx, &Lapy, &Lapz, dX, dY, dZ, i, j, k,
                              bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax);
    __syncthreads();

    double Kx, Ky, Kz;
    double fxx, fxy, fxz, fyy, fyz, fzz;
    load_field_to_smem(smem, Sphi, dims, block_i, block_j, block_k, SYM, SYM, SYM);
    __syncthreads();
    compute_all_derivs_smem(smem, &Kx, &Ky, &Kz, &fxx, &fxy, &fxz, &fyy, &fyz, &fzz,
                            dX, dY, dZ, i, j, k,
                            bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax);

    if (i >= ex0 || j >= ex1 || k >= ex2) return;
    int idx = flatten_bssn_index(i, j, k, ex0, ex1);

    double alpn1 = Lap[idx] + ONE;
    double chin1 = chi[idx] + ONE;
    double inv_chin1 = ONE / chin1;

    double pgxx = dxx[idx] + ONE;
    double pgxy = gxy[idx];
    double pgxz = gxz[idx];
    double pgyy = dyy[idx] + ONE;
    double pgyz = gyz[idx];
    double pgzz = dzz[idx] + ONE;

    double gupxx, gupxy, gupxz, gupyy, gupyz, gupzz;
    compute_inverse_metric(pgxx, pgxy, pgxz, pgyy, pgyz, pgzz, &gupxx, &gupxy, &gupxz, &gupyy, &gupyz, &gupzz);

    double chix_conn = (gupxx * chix + gupxy * chiy + gupxz * chiz) * HALF * inv_chin1;
    double chiy_conn = (gupxy * chix + gupyy * chiy + gupyz * chiz) * HALF * inv_chin1;
    double chiz_conn = (gupxz * chix + gupyz * chiy + gupzz * chiz) * HALF * inv_chin1;

    double lap_sphi = gupxx * fxx + gupyy * fyy + gupzz * fzz
        + TWO * (gupxy * fxy + gupxz * fxz + gupyz * fyz);
    double conn_term = (Gamx[idx] + chix_conn) * Kx
        + (Gamy[idx] + chiy_conn) * Ky
        + (Gamz[idx] + chiz_conn) * Kz;
    double rhs_spi = lap_sphi - conn_term;
    rhs_spi = rhs_spi * alpn1
        + gupxx * Lapx * Kx + gupxy * Lapx * Ky + gupxz * Lapx * Kz
        + gupxy * Lapy * Kx + gupyy * Lapy * Ky + gupyz * Lapy * Kz
        + gupxz * Lapz * Kx + gupyz * Lapz * Ky + gupzz * Lapz * Kz;

    double potential = 0.0;
    double dVdSphi = 0.0;
    compute_escalar_frpotential(Sphi[idx], a2, &potential, &dVdSphi);
    rhs_spi = rhs_spi * chin1 + alpn1 * (trK[idx] * Spi[idx] - dVdSphi);

    double rhs_sphi = alpn1 * Spi[idx];
    const int ext[3] = {ex0, ex1, ex2};
    rhs_sphi += d_lopsided_point(ext, Sphi, Sphi_rhs, betax, betay, betaz, X, Y, Z, symmetry, SYM, SYM, SYM, i, j, k);
    rhs_spi += d_lopsided_point(ext, Spi, Spi_rhs, betax, betay, betaz, X, Y, Z, symmetry, SYM, SYM, SYM, i, j, k);

    if (eps > 0.0) {
        rhs_sphi += d_kodis_point(ext, Sphi, X, Y, Z, SYM, SYM, SYM, symmetry, eps, i, j, k);
        rhs_spi += d_kodis_point(ext, Spi, X, Y, Z, SYM, SYM, SYM, symmetry, eps, i, j, k);
    }

    Sphi_rhs[idx] = rhs_sphi;
    Spi_rhs[idx] = rhs_spi;

    double kinetic = HALF * (gupxx * Kx * Kx + gupyy * Ky * Ky + gupzz * Kz * Kz)
        + gupxy * Kx * Ky + gupxz * Kx * Kz + gupyz * Ky * Kz;
    double lrho = chin1 * kinetic + HALF * Spi[idx] * Spi[idx] + potential;
    rho[idx] = lrho;
    Sx[idx] = -Spi[idx] * Kx;
    Sy[idx] = -Spi[idx] * Ky;
    Sz[idx] = -Spi[idx] * Kz;

    double stress_trace = (lrho - Spi[idx] * Spi[idx]) * inv_chin1;
    Sxx[idx] = Kx * Kx - stress_trace * pgxx;
    Sxy[idx] = Kx * Ky - stress_trace * pgxy;
    Sxz[idx] = Kx * Kz - stress_trace * pgxz;
    Syy[idx] = Ky * Ky - stress_trace * pgyy;
    Syz[idx] = Ky * Kz - stress_trace * pgyz;
    Szz[idx] = Kz * Kz - stress_trace * pgzz;

    (void)lev;
}

__launch_bounds__(256, 2)
__global__ void evolve_bssn_em_matter_kernel(
    const int ex0, const int ex1, const int ex2,
    const int symmetry, const int lev,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    const double* __restrict__ chi, const double* __restrict__ trK,
    const double* __restrict__ dxx, const double* __restrict__ gxy, const double* __restrict__ gxz,
    const double* __restrict__ dyy, const double* __restrict__ gyz, const double* __restrict__ dzz,
    const double* __restrict__ Lap,
    const double* __restrict__ betax, const double* __restrict__ betay, const double* __restrict__ betaz,
    const double* __restrict__ Ex, const double* __restrict__ Ey, const double* __restrict__ Ez,
    const double* __restrict__ Bx, const double* __restrict__ By, const double* __restrict__ Bz,
    const double* __restrict__ Kpsi, const double* __restrict__ Kphi,
    const double* __restrict__ Jx, const double* __restrict__ Jy, const double* __restrict__ Jz,
    const double* __restrict__ qchar,
    double* __restrict__ Ex_rhs, double* __restrict__ Ey_rhs, double* __restrict__ Ez_rhs,
    double* __restrict__ Bx_rhs, double* __restrict__ By_rhs, double* __restrict__ Bz_rhs,
    double* __restrict__ Kpsi_rhs, double* __restrict__ Kphi_rhs,
    double* __restrict__ rho, double* __restrict__ Sx, double* __restrict__ Sy, double* __restrict__ Sz,
    double* __restrict__ Sxx, double* __restrict__ Sxy, double* __restrict__ Sxz,
    double* __restrict__ Syy, double* __restrict__ Syz, double* __restrict__ Szz,
    double eps
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    int dims[3] = {ex0, ex1, ex2};
    double dX = X[1] - X[0];
    double dY = Y[1] - Y[0];
    double dZ = Z[1] - Z[0];
    StencilBounds bounds = make_stencil_bounds(ex0, ex1, ex2, symmetry, X, Y, Z, dX, dY, dZ, kBssnCenteredBoundaryRadius);

    int block_i = blockIdx.x * blockDim.x;
    int block_j = blockIdx.y * blockDim.y;
    int block_k = blockIdx.z * blockDim.z;

    BSSN_CENTERED_SMEM(smem);

#define LOAD_FIRST_DERIVS(field, fx, fy, fz, soa1, soa2, soa3) \
    load_field_to_smem(smem, field, dims, block_i, block_j, block_k, soa1, soa2, soa3); \
    __syncthreads(); \
    compute_first_derivs_smem(smem, &fx, &fy, &fz, dX, dY, dZ, i, j, k, bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax); \
    __syncthreads()

    double Lapx, Lapy, Lapz;
    double betaxx, betaxy, betaxz, betayx, betayy, betayz, betazx, betazy, betazz;
    double chix, chiy, chiz;
    double gxxx, gxxy, gxxz, gxyx, gxyy, gxyz, gxzx, gxzy, gxzz;
    double gyyx, gyyy, gyyz, gyzx, gyzy, gyzz, gzzx, gzzy, gzzz;
    double Kpsix, Kpsiy, Kpsiz, Kphix, Kphiy, Kphiz;
    double Exx, Exy, Exz, Eyx, Eyy, Eyz, Ezx, Ezy, Ezz;
    double Bxx, Bxy, Bxz, Byx, Byy, Byz, Bzx, Bzy, Bzz;

    LOAD_FIRST_DERIVS(Lap, Lapx, Lapy, Lapz, SYM, SYM, SYM);
    LOAD_FIRST_DERIVS(betax, betaxx, betaxy, betaxz, ANTI, SYM, SYM);
    LOAD_FIRST_DERIVS(betay, betayx, betayy, betayz, SYM, ANTI, SYM);
    LOAD_FIRST_DERIVS(betaz, betazx, betazy, betazz, SYM, SYM, ANTI);
    LOAD_FIRST_DERIVS(chi, chix, chiy, chiz, SYM, SYM, SYM);
    LOAD_FIRST_DERIVS(dxx, gxxx, gxxy, gxxz, SYM, SYM, SYM);
    LOAD_FIRST_DERIVS(gxy, gxyx, gxyy, gxyz, ANTI, ANTI, SYM);
    LOAD_FIRST_DERIVS(gxz, gxzx, gxzy, gxzz, ANTI, SYM, ANTI);
    LOAD_FIRST_DERIVS(dyy, gyyx, gyyy, gyyz, SYM, SYM, SYM);
    LOAD_FIRST_DERIVS(gyz, gyzx, gyzy, gyzz, SYM, ANTI, ANTI);
    LOAD_FIRST_DERIVS(dzz, gzzx, gzzy, gzzz, SYM, SYM, SYM);
    LOAD_FIRST_DERIVS(Kpsi, Kpsix, Kpsiy, Kpsiz, SYM, SYM, SYM);
    LOAD_FIRST_DERIVS(Kphi, Kphix, Kphiy, Kphiz, SYM, SYM, SYM);
    LOAD_FIRST_DERIVS(Ex, Exx, Exy, Exz, ANTI, SYM, SYM);
    LOAD_FIRST_DERIVS(Ey, Eyx, Eyy, Eyz, SYM, ANTI, SYM);
    LOAD_FIRST_DERIVS(Ez, Ezx, Ezy, Ezz, SYM, SYM, ANTI);
    LOAD_FIRST_DERIVS(Bx, Bxx, Bxy, Bxz, SYM, ANTI, ANTI);
    LOAD_FIRST_DERIVS(By, Byx, Byy, Byz, ANTI, SYM, ANTI);
    LOAD_FIRST_DERIVS(Bz, Bzx, Bzy, Bzz, ANTI, ANTI, SYM);

#undef LOAD_FIRST_DERIVS

    if (i >= ex0 || j >= ex1 || k >= ex2) return;
    int idx = flatten_bssn_index(i, j, k, ex0, ex1);

    double alpn1 = Lap[idx] + ONE;
    double chin1 = chi[idx] + ONE;
    double chi3o2 = sqrt(chin1) * sqrt(chin1) * sqrt(chin1);

    double pgxx = dxx[idx] + ONE;
    double pgxy = gxy[idx];
    double pgxz = gxz[idx];
    double pgyy = dyy[idx] + ONE;
    double pgyz = gyz[idx];
    double pgzz = dzz[idx] + ONE;

    pgxx /= chin1; pgxy /= chin1; pgxz /= chin1;
    pgyy /= chin1; pgyz /= chin1; pgzz /= chin1;

    gxxx = (gxxx - pgxx * chix) / chin1;
    gxxy = (gxxy - pgxx * chiy) / chin1;
    gxxz = (gxxz - pgxx * chiz) / chin1;
    gxyx = (gxyx - pgxy * chix) / chin1;
    gxyy = (gxyy - pgxy * chiy) / chin1;
    gxyz = (gxyz - pgxy * chiz) / chin1;
    gxzx = (gxzx - pgxz * chix) / chin1;
    gxzy = (gxzy - pgxz * chiy) / chin1;
    gxzz = (gxzz - pgxz * chiz) / chin1;
    gyyx = (gyyx - pgyy * chix) / chin1;
    gyyy = (gyyy - pgyy * chiy) / chin1;
    gyyz = (gyyz - pgyy * chiz) / chin1;
    gyzx = (gyzx - pgyz * chix) / chin1;
    gyzy = (gyzy - pgyz * chiy) / chin1;
    gyzz = (gyzz - pgyz * chiz) / chin1;
    gzzx = (gzzx - pgzz * chix) / chin1;
    gzzy = (gzzy - pgzz * chiy) / chin1;
    gzzz = (gzzz - pgzz * chiz) / chin1;

    double gupxx, gupxy, gupxz, gupyy, gupyz, gupzz;
    compute_inverse_metric(pgxx, pgxy, pgxz, pgyy, pgyz, pgzz, &gupxx, &gupxy, &gupxz, &gupyy, &gupyz, &gupzz);

    double lEx = Ex[idx], lEy = Ey[idx], lEz = Ez[idx];
    double lBx = Bx[idx], lBy = By[idx], lBz = Bz[idx];
    double ltrK = trK[idx];

    double rhs_Ex = alpn1 * ltrK * lEx - (lEx * betaxx + lEy * betaxy + lEz * betaxz)
        - FOUR_PI * alpn1 * Jx[idx] - alpn1 * (gupxx * Kpsix + gupxy * Kpsiy + gupxz * Kpsiz)
        + chi3o2 * (((pgxz * lBx + pgyz * lBy + pgzz * lBz) * Lapy + alpn1 * (pgxz * Bxy + pgyz * Byy + pgzz * Bzy) + alpn1 * (lBx * gxzy + lBy * gyzy + lBz * gzzy))
        - ((pgxy * lBx + pgyy * lBy + pgyz * lBz) * Lapz + alpn1 * (pgxy * Bxz + pgyy * Byz + pgyz * Bzz) + alpn1 * (lBx * gxyz + lBy * gyyz + lBz * gyzz)));
    double rhs_Ey = alpn1 * ltrK * lEy - (lEx * betayx + lEy * betayy + lEz * betayz)
        - FOUR_PI * alpn1 * Jy[idx] - alpn1 * (gupxy * Kpsix + gupyy * Kpsiy + gupyz * Kpsiz)
        + chi3o2 * (((pgxx * lBx + pgxy * lBy + pgxz * lBz) * Lapz + alpn1 * (pgxx * Bxz + pgxy * Byz + pgxz * Bzz) + alpn1 * (lBx * gxxz + lBy * gxyz + lBz * gxzz))
        - ((pgxz * lBx + pgyz * lBy + pgzz * lBz) * Lapx + alpn1 * (pgxz * Bxx + pgyz * Byx + pgzz * Bzx) + alpn1 * (lBx * gxzx + lBy * gyzx + lBz * gzzx)));
    double rhs_Ez = alpn1 * ltrK * lEz - (lEx * betazx + lEy * betazy + lEz * betazz)
        - FOUR_PI * alpn1 * Jz[idx] - alpn1 * (gupxz * Kpsix + gupyz * Kpsiy + gupzz * Kpsiz)
        + chi3o2 * (((pgxy * lBx + pgyy * lBy + pgyz * lBz) * Lapx + alpn1 * (pgxy * Bxx + pgyy * Byx + pgyz * Bzx) + alpn1 * (lBx * gxyx + lBy * gyyx + lBz * gyzx))
        - ((pgxx * lBx + pgxy * lBy + pgxz * lBz) * Lapy + alpn1 * (pgxx * Bxy + pgxy * Byy + pgxz * Bzy) + alpn1 * (lBx * gxxy + lBy * gxyy + lBz * gxzy)));

    double rhs_Bx = alpn1 * ltrK * lBx - (lBx * betaxx + lBy * betaxy + lBz * betaxz)
        - alpn1 * (gupxx * Kphix + gupxy * Kphiy + gupxz * Kphiz)
        - chi3o2 * (((pgxz * lEx + pgyz * lEy + pgzz * lEz) * Lapy + alpn1 * (pgxz * Exy + pgyz * Eyy + pgzz * Ezy) + alpn1 * (lEx * gxzy + lEy * gyzy + lEz * gzzy))
        - ((pgxy * lEx + pgyy * lEy + pgyz * lEz) * Lapz + alpn1 * (pgxy * Exz + pgyy * Eyz + pgyz * Ezz) + alpn1 * (lEx * gxyz + lEy * gyyz + lEz * gyzz)));
    double rhs_By = alpn1 * ltrK * lBy - (lBx * betayx + lBy * betayy + lBz * betayz)
        - alpn1 * (gupxy * Kphix + gupyy * Kphiy + gupyz * Kphiz)
        - chi3o2 * (((pgxx * lEx + pgxy * lEy + pgxz * lEz) * Lapz + alpn1 * (pgxx * Exz + pgxy * Eyz + pgxz * Ezz) + alpn1 * (lEx * gxxz + lEy * gxyz + lEz * gxzz))
        - ((pgxz * lEx + pgyz * lEy + pgzz * lEz) * Lapx + alpn1 * (pgxz * Exx + pgyz * Eyx + pgzz * Ezx) + alpn1 * (lEx * gxzx + lEy * gyzx + lEz * gzzx)));
    double rhs_Bz = alpn1 * ltrK * lBz - (lBx * betazx + lBy * betazy + lBz * betazz)
        - alpn1 * (gupxz * Kphix + gupyz * Kphiy + gupzz * Kphiz)
        - chi3o2 * (((pgxy * lEx + pgyy * lEy + pgyz * lEz) * Lapx + alpn1 * (pgxy * Exx + pgyy * Eyx + pgyz * Ezx) + alpn1 * (lEx * gxyx + lEy * gyyx + lEz * gyzx))
        - ((pgxx * lEx + pgxy * lEy + pgxz * lEz) * Lapy + alpn1 * (pgxx * Exy + pgxy * Eyy + pgxz * Ezy) + alpn1 * (lEx * gxxy + lEy * gxyy + lEz * gxzy)));

    double rhs_Kpsi = FOUR_PI * alpn1 * qchar[idx] - alpn1 * Kpsi[idx]
        - alpn1 * (Exx + Eyy + Ezz - F3o2 / chin1 * (chix * lEx + chiy * lEy + chiz * lEz));
    double rhs_Kphi = -alpn1 * Kphi[idx]
        - alpn1 * (Bxx + Byy + Bzz - F3o2 / chin1 * (chix * lBx + chiy * lBy + chiz * lBz));

    const int ext[3] = {ex0, ex1, ex2};
    rhs_Kpsi += d_lopsided_point(ext, Kpsi, Kpsi_rhs, betax, betay, betaz, X, Y, Z, symmetry, SYM, SYM, SYM, i, j, k);
    rhs_Kphi += d_lopsided_point(ext, Kphi, Kphi_rhs, betax, betay, betaz, X, Y, Z, symmetry, SYM, SYM, SYM, i, j, k);
    rhs_Ex += d_lopsided_point(ext, Ex, Ex_rhs, betax, betay, betaz, X, Y, Z, symmetry, ANTI, SYM, SYM, i, j, k);
    rhs_Ey += d_lopsided_point(ext, Ey, Ey_rhs, betax, betay, betaz, X, Y, Z, symmetry, SYM, ANTI, SYM, i, j, k);
    rhs_Ez += d_lopsided_point(ext, Ez, Ez_rhs, betax, betay, betaz, X, Y, Z, symmetry, SYM, SYM, ANTI, i, j, k);
    rhs_Bx += d_lopsided_point(ext, Bx, Bx_rhs, betax, betay, betaz, X, Y, Z, symmetry, SYM, ANTI, ANTI, i, j, k);
    rhs_By += d_lopsided_point(ext, By, By_rhs, betax, betay, betaz, X, Y, Z, symmetry, ANTI, SYM, ANTI, i, j, k);
    rhs_Bz += d_lopsided_point(ext, Bz, Bz_rhs, betax, betay, betaz, X, Y, Z, symmetry, ANTI, ANTI, SYM, i, j, k);

    if (eps > 0.0) {
        rhs_Kpsi += d_kodis_point(ext, Kpsi, X, Y, Z, SYM, SYM, SYM, symmetry, eps, i, j, k);
        rhs_Kphi += d_kodis_point(ext, Kphi, X, Y, Z, SYM, SYM, SYM, symmetry, eps, i, j, k);
        rhs_Ex += d_kodis_point(ext, Ex, X, Y, Z, ANTI, SYM, SYM, symmetry, eps, i, j, k);
        rhs_Ey += d_kodis_point(ext, Ey, X, Y, Z, SYM, ANTI, SYM, symmetry, eps, i, j, k);
        rhs_Ez += d_kodis_point(ext, Ez, X, Y, Z, SYM, SYM, ANTI, symmetry, eps, i, j, k);
        rhs_Bx += d_kodis_point(ext, Bx, X, Y, Z, SYM, ANTI, ANTI, symmetry, eps, i, j, k);
        rhs_By += d_kodis_point(ext, By, X, Y, Z, ANTI, SYM, ANTI, symmetry, eps, i, j, k);
        rhs_Bz += d_kodis_point(ext, Bz, X, Y, Z, ANTI, ANTI, SYM, symmetry, eps, i, j, k);
    }

    Ex_rhs[idx] = rhs_Ex; Ey_rhs[idx] = rhs_Ey; Ez_rhs[idx] = rhs_Ez;
    Bx_rhs[idx] = rhs_Bx; By_rhs[idx] = rhs_By; Bz_rhs[idx] = rhs_Bz;
    Kpsi_rhs[idx] = rhs_Kpsi; Kphi_rhs[idx] = rhs_Kphi;

    double lrho = (pgxx * (lEx * lEx + lBx * lBx) + pgyy * (lEy * lEy + lBy * lBy) + pgzz * (lEz * lEz + lBz * lBz)
        + TWO * (pgxy * (lEx * lEy + lBx * lBy) + pgxz * (lEx * lEz + lBx * lBz) + pgyz * (lEy * lEz + lBy * lBz))) / EIGHT_PI;
    rho[idx] = lrho;
    Sx[idx] = (lEy * lBz - lEz * lBy) / FOUR_PI / chi3o2;
    Sy[idx] = (lEz * lBx - lEx * lBz) / FOUR_PI / chi3o2;
    Sz[idx] = (lEx * lBy - lEy * lBx) / FOUR_PI / chi3o2;

    double llEx = pgxx * lEx + pgxy * lEy + pgxz * lEz;
    double llEy = pgxy * lEx + pgyy * lEy + pgyz * lEz;
    double llEz = pgxz * lEx + pgyz * lEy + pgzz * lEz;
    double llBx = pgxx * lBx + pgxy * lBy + pgxz * lBz;
    double llBy = pgxy * lBx + pgyy * lBy + pgyz * lBz;
    double llBz = pgxz * lBx + pgyz * lBy + pgzz * lBz;
    Sxx[idx] = lrho * pgxx - (llEx * llEx + llBx * llBx) / FOUR_PI;
    Sxy[idx] = lrho * pgxy - (llEx * llEy + llBx * llBy) / FOUR_PI;
    Sxz[idx] = lrho * pgxz - (llEx * llEz + llBx * llBz) / FOUR_PI;
    Syy[idx] = lrho * pgyy - (llEy * llEy + llBy * llBy) / FOUR_PI;
    Syz[idx] = lrho * pgyz - (llEy * llEz + llBy * llBz) / FOUR_PI;
    Szz[idx] = lrho * pgzz - (llEz * llEz + llBz * llBz) / FOUR_PI;
}

template <typename FieldReal>
__launch_bounds__(256, 2)
__global__ void assemble_conformal_geometry_kernel(
    const int ex0, const int ex1, const int ex2, 
    const int symmetry, const int lev,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    const FieldReal* __restrict__ dxx, const FieldReal* __restrict__ gxy, const FieldReal* __restrict__ gxz,
    const FieldReal* __restrict__ dyy, const FieldReal* __restrict__ gyz, const FieldReal* __restrict__ dzz,
    const FieldReal* __restrict__ Gamx, const FieldReal* __restrict__ Gamy, const FieldReal* __restrict__ Gamz,
    FieldReal* __restrict__ Gamxxx_out, FieldReal* __restrict__ Gamxxy_out, FieldReal* __restrict__ Gamxxz_out,
    FieldReal* __restrict__ Gamxyy_out, FieldReal* __restrict__ Gamxyz_out, FieldReal* __restrict__ Gamxzz_out,
    FieldReal* __restrict__ Gamyxx_out, FieldReal* __restrict__ Gamyxy_out, FieldReal* __restrict__ Gamyxz_out,
    FieldReal* __restrict__ Gamyyy_out, FieldReal* __restrict__ Gamyyz_out, FieldReal* __restrict__ Gamyzz_out,
    FieldReal* __restrict__ Gamzxx_out, FieldReal* __restrict__ Gamzxy_out, FieldReal* __restrict__ Gamzxz_out,
    FieldReal* __restrict__ Gamzyy_out, FieldReal* __restrict__ Gamzyz_out, FieldReal* __restrict__ Gamzzz_out,
    FieldReal* __restrict__ Rxx_out, FieldReal* __restrict__ Rxy_out, FieldReal* __restrict__ Rxz_out,
    FieldReal* __restrict__ Ryy_out, FieldReal* __restrict__ Ryz_out, FieldReal* __restrict__ Rzz_out
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    int dims[3] = {ex0, ex1, ex2};
    FieldReal dX = X[1] - X[0]; FieldReal dY = Y[1] - Y[0]; FieldReal dZ = Z[1] - Z[0];
    StencilBounds bounds = make_stencil_bounds(ex0, ex1, ex2, symmetry, X, Y, Z, dX, dY, dZ, kBssnCenteredBoundaryRadius);

    int block_i = blockIdx.x * blockDim.x;
    int block_j = blockIdx.y * blockDim.y;
    int block_k = blockIdx.z * blockDim.z;

    BSSN_CENTERED_SMEM_TYPED(FieldReal, smem);

    FieldReal dgxx_x, dgxx_y, dgxx_z, dxx_xx, dxx_xy, dxx_xz, dxx_yy, dxx_yz, dxx_zz;
    load_field_to_smem(smem, dxx, dims, block_i, block_j, block_k, SYM, SYM, SYM);
    __syncthreads(); compute_all_derivs_smem(smem, &dgxx_x, &dgxx_y, &dgxx_z, &dxx_xx, &dxx_xy, &dxx_xz, &dxx_yy, &dxx_yz, &dxx_zz, dX, dY, dZ, i, j, k, bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax); __syncthreads();

    FieldReal dgxy_x, dgxy_y, dgxy_z, gxy_xx, gxy_xy, gxy_xz, gxy_yy, gxy_yz, gxy_zz;
    load_field_to_smem(smem, gxy, dims, block_i, block_j, block_k, ANTI, ANTI, SYM);
    __syncthreads(); compute_all_derivs_smem(smem, &dgxy_x, &dgxy_y, &dgxy_z, &gxy_xx, &gxy_xy, &gxy_xz, &gxy_yy, &gxy_yz, &gxy_zz, dX, dY, dZ, i, j, k, bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax); __syncthreads();

    FieldReal dgxz_x, dgxz_y, dgxz_z, gxz_xx, gxz_xy, gxz_xz, gxz_yy, gxz_yz, gxz_zz;
    load_field_to_smem(smem, gxz, dims, block_i, block_j, block_k, ANTI, SYM, ANTI);
    __syncthreads(); compute_all_derivs_smem(smem, &dgxz_x, &dgxz_y, &dgxz_z, &gxz_xx, &gxz_xy, &gxz_xz, &gxz_yy, &gxz_yz, &gxz_zz, dX, dY, dZ, i, j, k, bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax); __syncthreads();

    FieldReal dgyy_x, dgyy_y, dgyy_z, dyy_xx, dyy_xy, dyy_xz, dyy_yy, dyy_yz, dyy_zz;
    load_field_to_smem(smem, dyy, dims, block_i, block_j, block_k, SYM, SYM, SYM);
    __syncthreads(); compute_all_derivs_smem(smem, &dgyy_x, &dgyy_y, &dgyy_z, &dyy_xx, &dyy_xy, &dyy_xz, &dyy_yy, &dyy_yz, &dyy_zz, dX, dY, dZ, i, j, k, bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax); __syncthreads();

    FieldReal dgyz_x, dgyz_y, dgyz_z, gyz_xx, gyz_xy, gyz_xz, gyz_yy, gyz_yz, gyz_zz;
    load_field_to_smem(smem, gyz, dims, block_i, block_j, block_k, SYM, ANTI, ANTI);
    __syncthreads(); compute_all_derivs_smem(smem, &dgyz_x, &dgyz_y, &dgyz_z, &gyz_xx, &gyz_xy, &gyz_xz, &gyz_yy, &gyz_yz, &gyz_zz, dX, dY, dZ, i, j, k, bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax); __syncthreads();

    FieldReal dgzz_x, dgzz_y, dgzz_z, dzz_xx, dzz_xy, dzz_xz, dzz_yy, dzz_yz, dzz_zz;
    load_field_to_smem(smem, dzz, dims, block_i, block_j, block_k, SYM, SYM, SYM);
    __syncthreads(); compute_all_derivs_smem(smem, &dgzz_x, &dgzz_y, &dgzz_z, &dzz_xx, &dzz_xy, &dzz_xz, &dzz_yy, &dzz_yz, &dzz_zz, dX, dY, dZ, i, j, k, bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax); __syncthreads();

    FieldReal dGamxx, dGamxy, dGamxz;
    load_field_to_smem(smem, Gamx, dims, block_i, block_j, block_k, ANTI, SYM, SYM);
    __syncthreads(); compute_first_derivs_smem(smem, &dGamxx, &dGamxy, &dGamxz, dX, dY, dZ, i, j, k, bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax); __syncthreads();

    FieldReal dGamyx, dGamyy, dGamyz;
    load_field_to_smem(smem, Gamy, dims, block_i, block_j, block_k, SYM, ANTI, SYM);
    __syncthreads(); compute_first_derivs_smem(smem, &dGamyx, &dGamyy, &dGamyz, dX, dY, dZ, i, j, k, bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax); __syncthreads();

    FieldReal dGamzx, dGamzy, dGamzz;
    load_field_to_smem(smem, Gamz, dims, block_i, block_j, block_k, SYM, SYM, ANTI);
    __syncthreads(); compute_first_derivs_smem(smem, &dGamzx, &dGamzy, &dGamzz, dX, dY, dZ, i, j, k, bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax);

    if (i >= ex0 || j >= ex1 || k >= ex2) return;
    int idx = flatten_bssn_index(i, j, k, ex0, ex1);

    const FieldReal HALF = FieldReal(0.5);
    const FieldReal TWO = FieldReal(2.0);
    const FieldReal ONE = FieldReal(1.0);

    FieldReal l_gxx = dxx[idx] + ONE; FieldReal l_gxy = gxy[idx]; FieldReal l_gxz = gxz[idx];
    FieldReal l_gyy = dyy[idx] + ONE; FieldReal l_gyz = gyz[idx]; FieldReal l_gzz = dzz[idx] + ONE;

    FieldReal det = l_gxx * l_gyy * l_gzz + l_gxy * l_gyz * l_gxz
               + l_gxz * l_gxy * l_gyz
               - l_gxz * l_gyy * l_gxz
               - l_gxy * l_gxy * l_gzz
               - l_gxx * l_gyz * l_gyz;
    FieldReal inv_det = ONE / det;

    FieldReal gupxx = (l_gyy * l_gzz - l_gyz * l_gyz) * inv_det;
    FieldReal gupxy = (l_gyz * l_gxz - l_gxy * l_gzz) * inv_det;
    FieldReal gupxz = (l_gxy * l_gyz - l_gyy * l_gxz) * inv_det;
    FieldReal gupyy = (l_gxx * l_gzz - l_gxz * l_gxz) * inv_det;
    FieldReal gupyz = (l_gxy * l_gxz - l_gxx * l_gyz) * inv_det;
    FieldReal gupzz = (l_gxx * l_gyy - l_gxy * l_gxy) * inv_det;

    FieldReal G1_xxx = HALF * dgxx_x;
    FieldReal G1_xxy = HALF * dgxx_y;
    FieldReal G1_xxz = HALF * dgxx_z;
    FieldReal G1_xyy = HALF * (TWO * dgxy_y - dgyy_x);
    FieldReal G1_xyz = HALF * (dgxy_z + dgxz_y - dgyz_x);
    FieldReal G1_xzz = HALF * (TWO * dgxz_z - dgzz_x);

    FieldReal G1_yxx = HALF * (TWO * dgxy_x - dgxx_y);
    FieldReal G1_yxy = HALF * dgyy_x;
    FieldReal G1_yxz = HALF * (dgxy_z + dgyz_x - dgxz_y);
    FieldReal G1_yyy = HALF * dgyy_y;
    FieldReal G1_yyz = HALF * dgyy_z;
    FieldReal G1_yzz = HALF * (TWO * dgyz_z - dgzz_y);

    FieldReal G1_zxx = HALF * (TWO * dgxz_x - dgxx_z);
    FieldReal G1_zxy = HALF * (dgxz_y + dgyz_x - dgxy_z);
    FieldReal G1_zxz = HALF * dgzz_x;
    FieldReal G1_zyy = HALF * (TWO * dgyz_y - dgyy_z);
    FieldReal G1_zyz = HALF * dgzz_y;
    FieldReal G1_zzz = HALF * dgzz_z;

    FieldReal G2_xxx = gupxx * G1_xxx + gupxy * G1_yxx + gupxz * G1_zxx;
    FieldReal G2_xxy = gupxx * G1_xxy + gupxy * G1_yxy + gupxz * G1_zxy;
    FieldReal G2_xxz = gupxx * G1_xxz + gupxy * G1_yxz + gupxz * G1_zxz;
    FieldReal G2_xyy = gupxx * G1_xyy + gupxy * G1_yyy + gupxz * G1_zyy;
    FieldReal G2_xyz = gupxx * G1_xyz + gupxy * G1_yyz + gupxz * G1_zyz;
    FieldReal G2_xzz = gupxx * G1_xzz + gupxy * G1_yzz + gupxz * G1_zzz;

    FieldReal G2_yxx = gupxy * G1_xxx + gupyy * G1_yxx + gupyz * G1_zxx;
    FieldReal G2_yxy = gupxy * G1_xxy + gupyy * G1_yxy + gupyz * G1_zxy;
    FieldReal G2_yxz = gupxy * G1_xxz + gupyy * G1_yxz + gupyz * G1_zxz;
    FieldReal G2_yyy = gupxy * G1_xyy + gupyy * G1_yyy + gupyz * G1_zyy;
    FieldReal G2_yyz = gupxy * G1_xyz + gupyy * G1_yyz + gupyz * G1_zyz;
    FieldReal G2_yzz = gupxy * G1_xzz + gupyy * G1_yzz + gupyz * G1_zzz;

    FieldReal G2_zxx = gupxz * G1_xxx + gupyz * G1_yxx + gupzz * G1_zxx;
    FieldReal G2_zxy = gupxz * G1_xxy + gupyz * G1_yxy + gupzz * G1_zxy;
    FieldReal G2_zxz = gupxz * G1_xxz + gupyz * G1_yxz + gupzz * G1_zxz;
    FieldReal G2_zyy = gupxz * G1_xyy + gupyz * G1_yyy + gupzz * G1_zyy;
    FieldReal G2_zyz = gupxz * G1_xyz + gupyz * G1_yyz + gupzz * G1_zyz;
    FieldReal G2_zzz = gupxz * G1_xzz + gupyz * G1_yzz + gupzz * G1_zzz;

    FieldReal Gamxa = gupxx * G2_xxx + gupyy * G2_xyy + gupzz * G2_xzz + TWO * (gupxy * G2_xxy + gupxz * G2_xxz + gupyz * G2_xyz);
    FieldReal Gamya = gupxx * G2_yxx + gupyy * G2_yyy + gupzz * G2_yzz + TWO * (gupxy * G2_yxy + gupxz * G2_yxz + gupyz * G2_yyz);
    FieldReal Gamza = gupxx * G2_zxx + gupyy * G2_zyy + gupzz * G2_zzz + TWO * (gupxy * G2_zxy + gupxz * G2_zxz + gupyz * G2_zyz);

    FieldReal l_Rxx = gupxx * dxx_xx + gupyy * dxx_yy + gupzz * dxx_zz + (gupxy * dxx_xy + gupxz * dxx_xz + gupyz * dxx_yz) * TWO;
    FieldReal l_Ryy = gupxx * dyy_xx + gupyy * dyy_yy + gupzz * dyy_zz + (gupxy * dyy_xy + gupxz * dyy_xz + gupyz * dyy_yz) * TWO;
    FieldReal l_Rzz = gupxx * dzz_xx + gupyy * dzz_yy + gupzz * dzz_zz + (gupxy * dzz_xy + gupxz * dzz_xz + gupyz * dzz_yz) * TWO;
    FieldReal l_Rxy = gupxx * gxy_xx + gupyy * gxy_yy + gupzz * gxy_zz + (gupxy * gxy_xy + gupxz * gxy_xz + gupyz * gxy_yz) * TWO;
    FieldReal l_Rxz = gupxx * gxz_xx + gupyy * gxz_yy + gupzz * gxz_zz + (gupxy * gxz_xy + gupxz * gxz_xz + gupyz * gxz_yz) * TWO;
    FieldReal l_Ryz = gupxx * gyz_xx + gupyy * gyz_yy + gupzz * gyz_zz + (gupxy * gyz_xy + gupxz * gyz_xz + gupyz * gyz_yz) * TWO;

    l_Rxx = -HALF * l_Rxx + l_gxx * dGamxx + l_gxy * dGamyx + l_gxz * dGamzx + Gamxa * G1_xxx + Gamya * G1_xxy + Gamza * G1_xxz + 
          gupxx * (TWO*(G2_xxx*G1_xxx + G2_yxx*G1_xxy + G2_zxx*G1_xxz) + G2_xxx*G1_xxx + G2_yxx*G1_yxx + G2_zxx*G1_zxx) +
          gupxy * (TWO*(G2_xxx*G1_xxy + G2_yxx*G1_xyy + G2_zxx*G1_xyz + G2_xxy*G1_xxx + G2_yxy*G1_xxy + G2_zxy*G1_xxz) + G2_xxy*G1_xxx + G2_yxy*G1_yxx + G2_zxy*G1_zxx + G2_xxx*G1_xxy + G2_yxx*G1_yxy + G2_zxx*G1_zxy) + 
          gupxz * (TWO*(G2_xxx*G1_xxz + G2_yxx*G1_xyz + G2_zxx*G1_xzz + G2_xxz*G1_xxx + G2_yxz*G1_xxy + G2_zxz*G1_xxz) + G2_xxz*G1_xxx + G2_yxz*G1_yxx + G2_zxz*G1_zxx + G2_xxx*G1_xxz + G2_yxx*G1_yxz + G2_zxx*G1_zxz) + 
          gupyy * (TWO*(G2_xxy*G1_xxy + G2_yxy*G1_xyy + G2_zxy*G1_xyz) + G2_xxy*G1_xxy + G2_yxy*G1_yxy + G2_zxy*G1_zxy) + 
          gupyz * (TWO*(G2_xxy*G1_xxz + G2_yxy*G1_xyz + G2_zxy*G1_xzz + G2_xxz*G1_xxy + G2_yxz*G1_xyy + G2_zxz*G1_xyz) + G2_xxz*G1_xxy + G2_yxz*G1_yxy + G2_zxz*G1_zxy + G2_xxy*G1_xxz + G2_yxy*G1_yxz + G2_zxy*G1_zxz) + 
          gupzz * (TWO*(G2_xxz*G1_xxz + G2_yxz*G1_xyz + G2_zxz*G1_xzz) + G2_xxz*G1_xxz + G2_yxz*G1_yxz + G2_zxz*G1_zxz);

    l_Ryy = -HALF * l_Ryy + l_gxy * dGamxy + l_gyy * dGamyy + l_gyz * dGamzy + Gamxa * G1_yxy + Gamya * G1_yyy + Gamza * G1_yyz + 
          gupxx * (TWO*(G2_xxy*G1_yxx + G2_yxy*G1_yxy + G2_zxy*G1_yxz) + G2_xxy*G1_xxy + G2_yxy*G1_yxy + G2_zxy*G1_zxy) + 
          gupxy * (TWO*(G2_xxy*G1_yxy + G2_yxy*G1_yyy + G2_zxy*G1_yyz + G2_xyy*G1_yxx + G2_yyy*G1_yxy + G2_zyy*G1_yxz) + G2_xyy*G1_xxy + G2_yyy*G1_yxy + G2_zyy*G1_zxy + G2_xxy*G1_xyy + G2_yxy*G1_yyy + G2_zxy*G1_zyy) + 
          gupxz * (TWO*(G2_xxy*G1_yxz + G2_yxy*G1_yyz + G2_zxy*G1_yzz + G2_xyz*G1_yxx + G2_yyz*G1_yxy + G2_zyz*G1_yxz) + G2_xyz*G1_xxy + G2_yyz*G1_yxy + G2_zyz*G1_zxy + G2_xxy*G1_xyz + G2_yxy*G1_yyz + G2_zxy*G1_zyz) + 
          gupyy * (TWO*(G2_xyy*G1_yxy + G2_yyy*G1_yyy + G2_zyy*G1_yyz) + G2_xyy*G1_xyy + G2_yyy*G1_yyy + G2_zyy*G1_zyy) + 
          gupyz * (TWO*(G2_xyy*G1_yxz + G2_yyy*G1_yyz + G2_zyy*G1_yzz + G2_xyz*G1_yxy + G2_yyz*G1_yyy + G2_zyz*G1_yyz) + G2_xyz*G1_xyy + G2_yyz*G1_yyy + G2_zyz*G1_zyy + G2_xyy*G1_xyz + G2_yyy*G1_yyz + G2_zyy*G1_zyz) + 
          gupzz * (TWO*(G2_xyz*G1_yxz + G2_yyz*G1_yyz + G2_zyz*G1_yzz) + G2_xyz*G1_xyz + G2_yyz*G1_yyz + G2_zyz*G1_zyz);

    l_Rzz = -HALF * l_Rzz + l_gxz * dGamxz + l_gyz * dGamyz + l_gzz * dGamzz + Gamxa * G1_zxz + Gamya * G1_zyz + Gamza * G1_zzz + 
          gupxx * (TWO*(G2_xxz*G1_zxx + G2_yxz*G1_zxy + G2_zxz*G1_zxz) + G2_xxz*G1_xxz + G2_yxz*G1_yxz + G2_zxz*G1_zxz) + 
          gupxy * (TWO*(G2_xxz*G1_zxy + G2_yxz*G1_zyy + G2_zxz*G1_zyz + G2_xyz*G1_zxx + G2_yyz*G1_zxy + G2_zyz*G1_zxz) + G2_xyz*G1_xxz + G2_yyz*G1_yxz + G2_zyz*G1_zxz + G2_xxz*G1_xyz + G2_yxz*G1_yyz + G2_zxz*G1_zyz) + 
          gupxz * (TWO*(G2_xxz*G1_zxz + G2_yxz*G1_zyz + G2_zxz*G1_zzz + G2_xzz*G1_zxx + G2_yzz*G1_zxy + G2_zzz*G1_zxz) + G2_xzz*G1_xxz + G2_yzz*G1_yxz + G2_zzz*G1_zxz + G2_xxz*G1_xzz + G2_yxz*G1_yzz + G2_zxz*G1_zzz) + 
          gupyy * (TWO*(G2_xyz*G1_zxy + G2_yyz*G1_zyy + G2_zyz*G1_zyz) + G2_xyz*G1_xyz + G2_yyz*G1_yyz + G2_zyz*G1_zyz) + 
          gupyz * (TWO*(G2_xyz*G1_zxz + G2_yyz*G1_zyz + G2_zyz*G1_zzz + G2_xzz*G1_zxy + G2_yzz*G1_zyy + G2_zzz*G1_zyz) + G2_xzz*G1_xyz + G2_yzz*G1_yyz + G2_zzz*G1_zyz + G2_xyz*G1_xzz + G2_yyz*G1_yzz + G2_zyz*G1_zzz) + 
          gupzz * (TWO*(G2_xzz*G1_zxz + G2_yzz*G1_zyz + G2_zzz*G1_zzz) + G2_xzz*G1_xzz + G2_yzz*G1_yzz + G2_zzz*G1_zzz);

    l_Rxy = HALF * ( - l_Rxy + l_gxx * dGamxy + l_gxy * dGamyy + l_gxz * dGamzy + l_gxy * dGamxx + l_gyy * dGamyx + l_gyz * dGamzx + Gamxa * G1_xxy + Gamya * G1_xyy + Gamza * G1_xyz + Gamxa * G1_yxx + Gamya * G1_yxy + Gamza * G1_yxz) + 
          gupxx * (G2_xxx*G1_yxx + G2_yxx*G1_yxy + G2_zxx*G1_yxz + G2_xxy*G1_xxx + G2_yxy*G1_xxy + G2_zxy*G1_xxz + G2_xxx*G1_xxy + G2_yxx*G1_yxy + G2_zxx*G1_zxy) + 
          gupxy * (G2_xxx*G1_yxy + G2_yxx*G1_yyy + G2_zxx*G1_yyz + G2_xxy*G1_xxy + G2_yxy*G1_xyy + G2_zxy*G1_xyz + G2_xxy*G1_xxy + G2_yxy*G1_yxy + G2_zxy*G1_zxy + G2_xxy*G1_yxx + G2_yxy*G1_yxy + G2_zxy*G1_yxz + G2_xyy*G1_xxx + G2_yyy*G1_xxy + G2_zyy*G1_xxz + G2_xxx*G1_xyy + G2_yxx*G1_yyy + G2_zxx*G1_zyy) + 
          gupxz * (G2_xxx*G1_yxz + G2_yxx*G1_yyz + G2_zxx*G1_yzz + G2_xxy*G1_xxz + G2_yxy*G1_xyz + G2_zxy*G1_xzz + G2_xxz*G1_xxy + G2_yxz*G1_yxy + G2_zxz*G1_zxy + G2_xxz*G1_yxx + G2_yxz*G1_yxy + G2_zxz*G1_yxz + G2_xyz*G1_xxx + G2_yyz*G1_xxy + G2_zyz*G1_xxz + G2_xxx*G1_xyz + G2_yxx*G1_yyz + G2_zxx*G1_zyz) + 
          gupyy * (G2_xxy*G1_yxy + G2_yxy*G1_yyy + G2_zxy*G1_yyz + G2_xyy*G1_xxy + G2_yyy*G1_xyy + G2_zyy*G1_xyz + G2_xxy*G1_xyy + G2_yxy*G1_yyy + G2_zxy*G1_zyy) + 
          gupyz * (G2_xxy*G1_yxz + G2_yxy*G1_yyz + G2_zxy*G1_yzz + G2_xyy*G1_xxz + G2_yyy*G1_xyz + G2_zyy*G1_xzz + G2_xxz*G1_xyy + G2_yxz*G1_yyy + G2_zxz*G1_zyy + G2_xxz*G1_yxy + G2_yxz*G1_yyy + G2_zxz*G1_yyz + G2_xyz*G1_xxy + G2_yyz*G1_xyy + G2_zyz*G1_xyz + G2_xxy*G1_xyz + G2_yxy*G1_yyz + G2_zxy*G1_zyz) + 
          gupzz * (G2_xxz*G1_yxz + G2_yxz*G1_yyz + G2_zxz*G1_yzz + G2_xyz*G1_xxz + G2_yyz*G1_xyz + G2_zyz*G1_xzz + G2_xxz*G1_xyz + G2_yxz*G1_yyz + G2_zxz*G1_zyz);

    l_Rxz = HALF * ( - l_Rxz + l_gxx * dGamxz + l_gxy * dGamyz + l_gxz * dGamzz + l_gxz * dGamxx + l_gyz * dGamyx + l_gzz * dGamzx + Gamxa * G1_xxz + Gamya * G1_xyz + Gamza * G1_xzz + Gamxa * G1_zxx + Gamya * G1_zxy + Gamza * G1_zxz) + 
          gupxx * (G2_xxx*G1_zxx + G2_yxx*G1_zxy + G2_zxx*G1_zxz + G2_xxz*G1_xxx + G2_yxz*G1_xxy + G2_zxz*G1_xxz + G2_xxx*G1_xxz + G2_yxx*G1_yxz + G2_zxx*G1_zxz) + 
          gupxy * (G2_xxx*G1_zxy + G2_yxx*G1_zyy + G2_zxx*G1_zyz + G2_xxz*G1_xxy + G2_yxz*G1_xyy + G2_zxz*G1_xyz + G2_xxy*G1_xxz + G2_yxy*G1_yxz + G2_zxy*G1_zxz + G2_xxy*G1_zxx + G2_yxy*G1_zxy + G2_zxy*G1_zxz + G2_xyz*G1_xxx + G2_yyz*G1_xxy + G2_zyz*G1_xxz + G2_xxx*G1_xyz + G2_yxx*G1_yyz + G2_zxx*G1_zyz) + 
          gupxz * (G2_xxx*G1_zxz + G2_yxx*G1_zyz + G2_zxx*G1_zzz + G2_xxz*G1_xxz + G2_yxz*G1_xyz + G2_zxz*G1_xzz + G2_xxz*G1_xxz + G2_yxz*G1_yxz + G2_zxz*G1_zxz + G2_xxz*G1_zxx + G2_yxz*G1_zxy + G2_zxz*G1_zxz + G2_xzz*G1_xxx + G2_yzz*G1_xxy + G2_zzz*G1_xxz + G2_xxx*G1_xzz + G2_yxx*G1_yzz + G2_zxx*G1_zzz) + 
          gupyy * (G2_xxy*G1_zxy + G2_yxy*G1_zyy + G2_zxy*G1_zyz + G2_xyz*G1_xxy + G2_yyz*G1_xyy + G2_zyz*G1_xyz + G2_xxy*G1_xyz + G2_yxy*G1_yyz + G2_zxy*G1_zyz) + 
          gupyz * (G2_xxy*G1_zxz + G2_yxy*G1_zyz + G2_zxy*G1_zzz + G2_xyz*G1_xxz + G2_yyz*G1_xyz + G2_zyz*G1_xzz + G2_xxz*G1_xyz + G2_yxz*G1_yyz + G2_zxz*G1_zyz + G2_xxz*G1_zxy + G2_yxz*G1_zyy + G2_zxz*G1_zyz + G2_xzz*G1_xxy + G2_yzz*G1_xyy + G2_zzz*G1_xyz + G2_xxy*G1_xzz + G2_yxy*G1_yzz + G2_zxy*G1_zzz) + 
          gupzz * (G2_xxz*G1_zxz + G2_yxz*G1_zyz + G2_zxz*G1_zzz + G2_xzz*G1_xxz + G2_yzz*G1_xyz + G2_zzz*G1_xzz + G2_xxz*G1_xzz + G2_yxz*G1_yzz + G2_zxz*G1_zzz);

    l_Ryz = HALF * ( - l_Ryz + l_gxy * dGamxz + l_gyy * dGamyz + l_gyz * dGamzz + l_gxz * dGamxy + l_gyz * dGamyy + l_gzz * dGamzy + Gamxa * G1_yxz + Gamya * G1_yyz + Gamza * G1_yzz + Gamxa * G1_zxy + Gamya * G1_zyy + Gamza * G1_zyz) + 
          gupxx * (G2_xxy*G1_zxx + G2_yxy*G1_zxy + G2_zxy*G1_zxz + G2_xxz*G1_yxx + G2_yxz*G1_yxy + G2_zxz*G1_yxz + G2_xxy*G1_xxz + G2_yxy*G1_yxz + G2_zxy*G1_zxz) + 
          gupxy * (G2_xxy*G1_zxy + G2_yxy*G1_zyy + G2_zxy*G1_zyz + G2_xxz*G1_yxy + G2_yxz*G1_yyy + G2_zxz*G1_yyz + G2_xyy*G1_xxz + G2_yyy*G1_yxz + G2_zyy*G1_zxz + G2_xyy*G1_zxx + G2_yyy*G1_zxy + G2_zyy*G1_zxz + G2_xyz*G1_yxx + G2_yyz*G1_yxy + G2_zyz*G1_yxz + G2_xxy*G1_xyz + G2_yxy*G1_yyz + G2_zxy*G1_zyz) + 
          gupxz * (G2_xxy*G1_zxz + G2_yxy*G1_zyz + G2_zxy*G1_zzz + G2_xxz*G1_yxz + G2_yxz*G1_yyz + G2_zxz*G1_yzz + G2_xyz*G1_xxz + G2_yyz*G1_yxz + G2_zyz*G1_zxz + G2_xyz*G1_zxx + G2_yyz*G1_zxy + G2_zyz*G1_zxz + G2_xzz*G1_yxx + G2_yzz*G1_yxy + G2_zzz*G1_yxz + G2_xxy*G1_xzz + G2_yxy*G1_yzz + G2_zxy*G1_zzz) + 
          gupyy * (G2_xyy*G1_zxy + G2_yyy*G1_zyy + G2_zyy*G1_zyz + G2_xyz*G1_yxy + G2_yyz*G1_yyy + G2_zyz*G1_yyz + G2_xyy*G1_xyz + G2_yyy*G1_yyz + G2_zyy*G1_zyz) + 
          gupyz * (G2_xyy*G1_zxz + G2_yyy*G1_zyz + G2_zyy*G1_zzz + G2_xyz*G1_yxz + G2_yyz*G1_yyz + G2_zyz*G1_yzz + G2_xyz*G1_xyz + G2_yyz*G1_yyz + G2_zyz*G1_zyz + G2_xyz*G1_zxy + G2_yyz*G1_zyy + G2_zyz*G1_zyz + G2_xzz*G1_yxy + G2_yzz*G1_yyy + G2_zzz*G1_yyz + G2_xyy*G1_xzz + G2_yyy*G1_yzz + G2_zyy*G1_zzz) + 
          gupzz * (G2_xyz*G1_zxz + G2_yyz*G1_zyz + G2_zyz*G1_zzz + G2_xzz*G1_yxz + G2_yzz*G1_yyz + G2_zzz*G1_yzz + G2_xyz*G1_xzz + G2_yyz*G1_yzz + G2_zyz*G1_zzz);

    Gamxxx_out[idx] = G2_xxx; Gamxxy_out[idx] = G2_xxy; Gamxxz_out[idx] = G2_xxz;
    Gamxyy_out[idx] = G2_xyy; Gamxyz_out[idx] = G2_xyz; Gamxzz_out[idx] = G2_xzz;
    Gamyxx_out[idx] = G2_yxx; Gamyxy_out[idx] = G2_yxy; Gamyxz_out[idx] = G2_yxz;
    Gamyyy_out[idx] = G2_yyy; Gamyyz_out[idx] = G2_yyz; Gamyzz_out[idx] = G2_yzz;
    Gamzxx_out[idx] = G2_zxx; Gamzxy_out[idx] = G2_zxy; Gamzxz_out[idx] = G2_zxz;
    Gamzyy_out[idx] = G2_zyy; Gamzyz_out[idx] = G2_zyz; Gamzzz_out[idx] = G2_zzz;

    Rxx_out[idx] = l_Rxx; Ryy_out[idx] = l_Ryy; Rzz_out[idx] = l_Rzz;
    Rxy_out[idx] = l_Rxy; Rxz_out[idx] = l_Rxz; Ryz_out[idx] = l_Ryz;
}

template <typename FieldReal>
__launch_bounds__(256, 2)
__global__ void evolve_bssn_metric_fields_kernel(
    const int ex0, const int ex1, const int ex2,
    const int symmetry,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    const FieldReal* __restrict__ chi, const FieldReal* __restrict__ trK,
    const FieldReal* __restrict__ dxx, const FieldReal* __restrict__ gxy, const FieldReal* __restrict__ gxz,
    const FieldReal* __restrict__ dyy, const FieldReal* __restrict__ gyz, const FieldReal* __restrict__ dzz,
    const FieldReal* __restrict__ Axx, const FieldReal* __restrict__ Axy, const FieldReal* __restrict__ Axz,
    const FieldReal* __restrict__ Ayy, const FieldReal* __restrict__ Ayz, const FieldReal* __restrict__ Azz,
    const FieldReal* __restrict__ Lap,
    const FieldReal* __restrict__ betax, const FieldReal* __restrict__ betay, const FieldReal* __restrict__ betaz,
    FieldReal* __restrict__ chi_rhs,
    FieldReal* __restrict__ gxx_rhs, FieldReal* __restrict__ gxy_rhs, FieldReal* __restrict__ gxz_rhs,
    FieldReal* __restrict__ gyy_rhs, FieldReal* __restrict__ gyz_rhs, FieldReal* __restrict__ gzz_rhs
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    int dims[3] = {ex0, ex1, ex2};
    FieldReal dX = X[1] - X[0];
    FieldReal dY = Y[1] - Y[0];
    FieldReal dZ = Z[1] - Z[0];
    StencilBounds bounds = make_stencil_bounds(ex0, ex1, ex2, symmetry, X, Y, Z, dX, dY, dZ, kBssnCenteredBoundaryRadius);

    int block_i = blockIdx.x * blockDim.x;
    int block_j = blockIdx.y * blockDim.y;
    int block_k = blockIdx.z * blockDim.z;

    BSSN_CENTERED_SMEM_TYPED(FieldReal, smem);

    FieldReal betaxx, betaxy, betaxz;
    load_field_to_smem(smem, betax, dims, block_i, block_j, block_k, ANTI, SYM, SYM);
    __syncthreads(); compute_first_derivs_smem(smem, &betaxx, &betaxy, &betaxz, dX, dY, dZ, i, j, k, bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax); __syncthreads();

    FieldReal betayx, betayy, betayz;
    load_field_to_smem(smem, betay, dims, block_i, block_j, block_k, SYM, ANTI, SYM);
    __syncthreads(); compute_first_derivs_smem(smem, &betayx, &betayy, &betayz, dX, dY, dZ, i, j, k, bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax); __syncthreads();

    FieldReal betazx, betazy, betazz;
    load_field_to_smem(smem, betaz, dims, block_i, block_j, block_k, SYM, SYM, ANTI);
    __syncthreads(); compute_first_derivs_smem(smem, &betazx, &betazy, &betazz, dX, dY, dZ, i, j, k, bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax); __syncthreads();

    if (i >= ex0 || j >= ex1 || k >= ex2) return;
    int idx = flatten_bssn_index(i, j, k, ex0, ex1);

    const FieldReal ONE = FieldReal(1.0);
    const FieldReal TWO = FieldReal(2.0);
    const FieldReal F1o3 = FieldReal(1.0) / FieldReal(3.0);
    const FieldReal F2o3 = FieldReal(2.0) / FieldReal(3.0);

    FieldReal val_Lap = Lap[idx];
    FieldReal val_chi = chi[idx];
    FieldReal val_trK = trK[idx];
    FieldReal alpn1 = val_Lap + ONE;
    FieldReal chin1 = val_chi + ONE;

    FieldReal l_gxx = dxx[idx] + ONE;
    FieldReal l_gxy = gxy[idx];
    FieldReal l_gxz = gxz[idx];
    FieldReal l_gyy = dyy[idx] + ONE;
    FieldReal l_gyz = gyz[idx];
    FieldReal l_gzz = dzz[idx] + ONE;
    FieldReal l_Axx = Axx[idx];
    FieldReal l_Axy = Axy[idx];
    FieldReal l_Axz = Axz[idx];
    FieldReal l_Ayy = Ayy[idx];
    FieldReal l_Ayz = Ayz[idx];
    FieldReal l_Azz = Azz[idx];

    FieldReal div_beta = betaxx + betayy + betazz;

    chi_rhs[idx] = F2o3 * chin1 * (alpn1 * val_trK - div_beta);

    gxx_rhs[idx] = -TWO * alpn1 * l_Axx - F2o3 * l_gxx * div_beta + TWO * (l_gxx * betaxx + l_gxy * betayx + l_gxz * betazx);
    gyy_rhs[idx] = -TWO * alpn1 * l_Ayy - F2o3 * l_gyy * div_beta + TWO * (l_gxy * betaxy + l_gyy * betayy + l_gyz * betazy);
    gzz_rhs[idx] = -TWO * alpn1 * l_Azz - F2o3 * l_gzz * div_beta + TWO * (l_gxz * betaxz + l_gyz * betayz + l_gzz * betazz);
    gxy_rhs[idx] = -TWO * alpn1 * l_Axy + F1o3 * l_gxy * div_beta + l_gxx * betaxy + l_gxz * betazy + l_gyy * betayx + l_gyz * betazx - l_gxy * betazz;
    gyz_rhs[idx] = -TWO * alpn1 * l_Ayz + F1o3 * l_gyz * div_beta + l_gxy * betaxz + l_gyy * betayz + l_gxz * betaxy + l_gzz * betazy - l_gyz * betaxx;
    gxz_rhs[idx] = -TWO * alpn1 * l_Axz + F1o3 * l_gxz * div_beta + l_gxx * betaxz + l_gxy * betayz + l_gyz * betayx + l_gzz * betazx - l_gxz * betayy;
}

template <typename FieldReal>
__launch_bounds__(256, 2)
__global__ void evolve_bssn_curvature_fields_kernel(
    const int ex0, const int ex1, const int ex2,
    const int symmetry,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    const FieldReal* __restrict__ chi, const FieldReal* __restrict__ trK,
    const FieldReal* __restrict__ dxx, const FieldReal* __restrict__ gxy, const FieldReal* __restrict__ gxz,
    const FieldReal* __restrict__ dyy, const FieldReal* __restrict__ gyz, const FieldReal* __restrict__ dzz,
    const FieldReal* __restrict__ Axx, const FieldReal* __restrict__ Axy, const FieldReal* __restrict__ Axz,
    const FieldReal* __restrict__ Ayy, const FieldReal* __restrict__ Ayz, const FieldReal* __restrict__ Azz,
    const FieldReal* __restrict__ Lap,
    const FieldReal* __restrict__ betax, const FieldReal* __restrict__ betay, const FieldReal* __restrict__ betaz,
    const FieldReal* __restrict__ rho,
    const FieldReal* __restrict__ Sxx, const FieldReal* __restrict__ Sxy, const FieldReal* __restrict__ Sxz,
    const FieldReal* __restrict__ Syy, const FieldReal* __restrict__ Syz, const FieldReal* __restrict__ Szz,
    const FieldReal* __restrict__ Gamxxx_in, const FieldReal* __restrict__ Gamxxy_in, const FieldReal* __restrict__ Gamxxz_in,
    const FieldReal* __restrict__ Gamxyy_in, const FieldReal* __restrict__ Gamxyz_in, const FieldReal* __restrict__ Gamxzz_in,
    const FieldReal* __restrict__ Gamyxx_in, const FieldReal* __restrict__ Gamyxy_in, const FieldReal* __restrict__ Gamyxz_in,
    const FieldReal* __restrict__ Gamyyy_in, const FieldReal* __restrict__ Gamyyz_in, const FieldReal* __restrict__ Gamyzz_in,
    const FieldReal* __restrict__ Gamzxx_in, const FieldReal* __restrict__ Gamzxy_in, const FieldReal* __restrict__ Gamzxz_in,
    const FieldReal* __restrict__ Gamzyy_in, const FieldReal* __restrict__ Gamzyz_in, const FieldReal* __restrict__ Gamzzz_in,
    const FieldReal* __restrict__ Rxx_in, const FieldReal* __restrict__ Rxy_in, const FieldReal* __restrict__ Rxz_in,
    const FieldReal* __restrict__ Ryy_in, const FieldReal* __restrict__ Ryz_in, const FieldReal* __restrict__ Rzz_in,
    FieldReal* __restrict__ trK_rhs,
    FieldReal* __restrict__ Axx_rhs, FieldReal* __restrict__ Axy_rhs, FieldReal* __restrict__ Axz_rhs,
    FieldReal* __restrict__ Ayy_rhs, FieldReal* __restrict__ Ayz_rhs, FieldReal* __restrict__ Azz_rhs,
    bool use_escalar_trace_free
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    int dims[3] = {ex0, ex1, ex2};
    FieldReal dX = X[1] - X[0];
    FieldReal dY = Y[1] - Y[0];
    FieldReal dZ = Z[1] - Z[0];
    StencilBounds bounds = make_stencil_bounds(ex0, ex1, ex2, symmetry, X, Y, Z, dX, dY, dZ, kBssnCenteredBoundaryRadius);

    int block_i = blockIdx.x * blockDim.x;
    int block_j = blockIdx.y * blockDim.y;
    int block_k = blockIdx.z * blockDim.z;

    BSSN_CENTERED_SMEM_TYPED(FieldReal, smem);

    FieldReal betaxx, betaxy, betaxz;
    load_field_to_smem(smem, betax, dims, block_i, block_j, block_k, ANTI, SYM, SYM);
    __syncthreads(); compute_first_derivs_smem(smem, &betaxx, &betaxy, &betaxz, dX, dY, dZ, i, j, k, bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax); __syncthreads();

    FieldReal betayx, betayy, betayz;
    load_field_to_smem(smem, betay, dims, block_i, block_j, block_k, SYM, ANTI, SYM);
    __syncthreads(); compute_first_derivs_smem(smem, &betayx, &betayy, &betayz, dX, dY, dZ, i, j, k, bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax); __syncthreads();

    FieldReal betazx, betazy, betazz;
    load_field_to_smem(smem, betaz, dims, block_i, block_j, block_k, SYM, SYM, ANTI);
    __syncthreads(); compute_first_derivs_smem(smem, &betazx, &betazy, &betazz, dX, dY, dZ, i, j, k, bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax); __syncthreads();

    FieldReal Lapx, Lapy, Lapz, fxx_Lap, fxy_Lap, fxz_Lap, fyy_Lap, fyz_Lap, fzz_Lap;
    load_field_to_smem(smem, Lap, dims, block_i, block_j, block_k, SYM, SYM, SYM);
    __syncthreads(); compute_all_derivs_smem(smem, &Lapx, &Lapy, &Lapz, &fxx_Lap, &fxy_Lap, &fxz_Lap, &fyy_Lap, &fyz_Lap, &fzz_Lap, dX, dY, dZ, i, j, k, bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax); __syncthreads();

    FieldReal chix, chiy, chiz, chixx, chixy, chixz, chiyy, chiyz, chizz;
    load_field_to_smem(smem, chi, dims, block_i, block_j, block_k, SYM, SYM, SYM);
    __syncthreads(); compute_all_derivs_smem(smem, &chix, &chiy, &chiz, &chixx, &chixy, &chixz, &chiyy, &chiyz, &chizz, dX, dY, dZ, i, j, k, bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax);

    if (i >= ex0 || j >= ex1 || k >= ex2) return;
    int idx = flatten_bssn_index(i, j, k, ex0, ex1);

    const FieldReal ONE = FieldReal(1.0);
    const FieldReal TWO = FieldReal(2.0);
    const FieldReal HALF = FieldReal(0.5);
    const FieldReal F1o3 = FieldReal(1.0) / FieldReal(3.0);
    const FieldReal F2o3 = FieldReal(2.0) / FieldReal(3.0);
    const FieldReal F3o2 = FieldReal(1.5);
    const FieldReal FOUR_PI = FieldReal(4.0) * FieldReal(M_PI);
    const FieldReal EIGHT_PI = FieldReal(8.0) * FieldReal(M_PI);
    const FieldReal F16_PI = FieldReal(16.0) * FieldReal(M_PI);

    FieldReal val_Lap = Lap[idx];
    FieldReal val_chi = chi[idx];
    FieldReal val_trK = trK[idx];
    FieldReal alpn1 = val_Lap + ONE;
    FieldReal chin1 = val_chi + ONE;
    FieldReal inv_chin1 = ONE / chin1;

    FieldReal l_gxx = dxx[idx] + ONE;
    FieldReal l_gxy = gxy[idx];
    FieldReal l_gxz = gxz[idx];
    FieldReal l_gyy = dyy[idx] + ONE;
    FieldReal l_gyz = gyz[idx];
    FieldReal l_gzz = dzz[idx] + ONE;
    FieldReal l_Axx = Axx[idx];
    FieldReal l_Axy = Axy[idx];
    FieldReal l_Axz = Axz[idx];
    FieldReal l_Ayy = Ayy[idx];
    FieldReal l_Ayz = Ayz[idx];
    FieldReal l_Azz = Azz[idx];

    FieldReal div_beta = betaxx + betayy + betazz;

    FieldReal gupxx, gupxy, gupxz, gupyy, gupyz, gupzz;
    compute_inverse_metric_t<FieldReal>(l_gxx, l_gxy, l_gxz, l_gyy, l_gyz, l_gzz, &gupxx, &gupxy, &gupxz, &gupyy, &gupyz, &gupzz);

    FieldReal Au_d_xx, Au_d_xy, Au_d_xz, Au_d_yx, Au_d_yy, Au_d_yz, Au_d_zx, Au_d_zy, Au_d_zz;
    FieldReal Aup_xx, Aup_xy, Aup_xz, Aup_yy, Aup_yz, Aup_zz;
    compute_raised_extrinsic_curvature_t<FieldReal>(
        l_Axx, l_Axy, l_Axz, l_Ayy, l_Ayz, l_Azz,
        gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
        &Au_d_xx, &Au_d_xy, &Au_d_xz, &Au_d_yx, &Au_d_yy, &Au_d_yz, &Au_d_zx, &Au_d_zy, &Au_d_zz,
        &Aup_xx, &Aup_xy, &Aup_xz, &Aup_yy, &Aup_yz, &Aup_zz
    );

    FieldReal l_Gamxxx = Gamxxx_in[idx];
    FieldReal l_Gamxxy = Gamxxy_in[idx];
    FieldReal l_Gamxxz = Gamxxz_in[idx];
    FieldReal l_Gamxyy = Gamxyy_in[idx];
    FieldReal l_Gamxyz = Gamxyz_in[idx];
    FieldReal l_Gamxzz = Gamxzz_in[idx];
    FieldReal l_Gamyxx = Gamyxx_in[idx];
    FieldReal l_Gamyxy = Gamyxy_in[idx];
    FieldReal l_Gamyxz = Gamyxz_in[idx];
    FieldReal l_Gamyyy = Gamyyy_in[idx];
    FieldReal l_Gamyyz = Gamyyz_in[idx];
    FieldReal l_Gamyzz = Gamyzz_in[idx];
    FieldReal l_Gamzxx = Gamzxx_in[idx];
    FieldReal l_Gamzxy = Gamzxy_in[idx];
    FieldReal l_Gamzxz = Gamzxz_in[idx];
    FieldReal l_Gamzyy = Gamzyy_in[idx];
    FieldReal l_Gamzyz = Gamzyz_in[idx];
    FieldReal l_Gamzzz = Gamzzz_in[idx];

    FieldReal l_Rxx = Rxx_in[idx];
    FieldReal l_Rxy = Rxy_in[idx];
    FieldReal l_Rxz = Rxz_in[idx];
    FieldReal l_Ryy = Ryy_in[idx];
    FieldReal l_Ryz = Ryz_in[idx];
    FieldReal l_Rzz = Rzz_in[idx];

    FieldReal fxx = chixx - l_Gamxxx * chix - l_Gamyxx * chiy - l_Gamzxx * chiz;
    FieldReal fxy = chixy - l_Gamxxy * chix - l_Gamyxy * chiy - l_Gamzxy * chiz;
    FieldReal fxz = chixz - l_Gamxxz * chix - l_Gamyxz * chiy - l_Gamzxz * chiz;
    FieldReal fyy = chiyy - l_Gamxyy * chix - l_Gamyyy * chiy - l_Gamzyy * chiz;
    FieldReal fyz = chiyz - l_Gamxyz * chix - l_Gamyyz * chiy - l_Gamzyz * chiz;
    FieldReal fzz = chizz - l_Gamxzz * chix - l_Gamyzz * chiy - l_Gamzzz * chiz;

    FieldReal f_scalar = gupxx * (fxx - F3o2 * inv_chin1 * chix * chix) + gupyy * (fyy - F3o2 * inv_chin1 * chiy * chiy) + gupzz * (fzz - F3o2 * inv_chin1 * chiz * chiz) +
        TWO * (gupxy * (fxy - F3o2 * inv_chin1 * chix * chiy) + gupxz * (fxz - F3o2 * inv_chin1 * chix * chiz) + gupyz * (fyz - F3o2 * inv_chin1 * chiy * chiz));

    l_Rxx += (fxx - chix * chix * inv_chin1 * HALF + l_gxx * f_scalar) * inv_chin1 * HALF;
    l_Ryy += (fyy - chiy * chiy * inv_chin1 * HALF + l_gyy * f_scalar) * inv_chin1 * HALF;
    l_Rzz += (fzz - chiz * chiz * inv_chin1 * HALF + l_gzz * f_scalar) * inv_chin1 * HALF;
    l_Rxy += (fxy - chix * chiy * inv_chin1 * HALF + l_gxy * f_scalar) * inv_chin1 * HALF;
    l_Rxz += (fxz - chix * chiz * inv_chin1 * HALF + l_gxz * f_scalar) * inv_chin1 * HALF;
    l_Ryz += (fyz - chiy * chiz * inv_chin1 * HALF + l_gyz * f_scalar) * inv_chin1 * HALF;

    FieldReal gx_phy = (gupxx * chix + gupxy * chiy + gupxz * chiz) * inv_chin1;
    FieldReal gy_phy = (gupxy * chix + gupyy * chiy + gupyz * chiz) * inv_chin1;
    FieldReal gz_phy = (gupxz * chix + gupyz * chiy + gupzz * chiz) * inv_chin1;
    FieldReal chi_dot_Lap = gx_phy * Lapx + gy_phy * Lapy + gz_phy * Lapz;

    FieldReal D2_Lap_xx = fxx_Lap - l_Gamxxx * Lapx - l_Gamyxx * Lapy - l_Gamzxx * Lapz;
    FieldReal D2_Lap_yy = fyy_Lap - l_Gamxyy * Lapx - l_Gamyyy * Lapy - l_Gamzyy * Lapz;
    FieldReal D2_Lap_zz = fzz_Lap - l_Gamxzz * Lapx - l_Gamyzz * Lapy - l_Gamzzz * Lapz;
    FieldReal D2_Lap_xy = fxy_Lap - l_Gamxxy * Lapx - l_Gamyxy * Lapy - l_Gamzxy * Lapz;
    FieldReal D2_Lap_xz = fxz_Lap - l_Gamxxz * Lapx - l_Gamyxz * Lapy - l_Gamzxz * Lapz;
    FieldReal D2_Lap_yz = fyz_Lap - l_Gamxyz * Lapx - l_Gamyyz * Lapy - l_Gamzyz * Lapz;

    FieldReal fxx_Lap_phy = D2_Lap_xx + chix * Lapx * inv_chin1 - HALF * l_gxx * chi_dot_Lap;
    FieldReal fyy_Lap_phy = D2_Lap_yy + chiy * Lapy * inv_chin1 - HALF * l_gyy * chi_dot_Lap;
    FieldReal fzz_Lap_phy = D2_Lap_zz + chiz * Lapz * inv_chin1 - HALF * l_gzz * chi_dot_Lap;
    FieldReal fxy_Lap_phy = D2_Lap_xy + HALF * (chix * Lapy + chiy * Lapx) * inv_chin1 - HALF * l_gxy * chi_dot_Lap;
    FieldReal fxz_Lap_phy = D2_Lap_xz + HALF * (chix * Lapz + chiz * Lapx) * inv_chin1 - HALF * l_gxz * chi_dot_Lap;
    FieldReal fyz_Lap_phy = D2_Lap_yz + HALF * (chiy * Lapz + chiz * Lapy) * inv_chin1 - HALF * l_gyz * chi_dot_Lap;

    FieldReal conf_trK_rhs = gupxx * D2_Lap_xx + gupyy * D2_Lap_yy + gupzz * D2_Lap_zz + TWO * (gupxy * D2_Lap_xy + gupxz * D2_Lap_xz + gupyz * D2_Lap_yz);
    FieldReal trK_rhs_val = conf_trK_rhs - HALF * chi_dot_Lap;

    FieldReal l_Sxx = Sxx[idx];
    FieldReal l_Sxy = Sxy[idx];
    FieldReal l_Sxz = Sxz[idx];
    FieldReal l_Syy = Syy[idx];
    FieldReal l_Syz = Syz[idx];
    FieldReal l_Szz = Szz[idx];
    FieldReal S = chin1 * (gupxx * l_Sxx + gupyy * l_Syy + gupzz * l_Szz + TWO * (gupxy * l_Sxy + gupxz * l_Sxz + gupyz * l_Syz));

    FieldReal term_xx = l_Axx * Au_d_xx + l_Axy * Au_d_yx + l_Axz * Au_d_zx;
    FieldReal term_xy = l_Axx * Au_d_xy + l_Axy * Au_d_yy + l_Axz * Au_d_zy;
    FieldReal term_xz = l_Axx * Au_d_xz + l_Axy * Au_d_yz + l_Axz * Au_d_zz;
    FieldReal term_yy = l_Axy * Au_d_xy + l_Ayy * Au_d_yy + l_Ayz * Au_d_zy;
    FieldReal term_yz = l_Axy * Au_d_xz + l_Ayy * Au_d_yz + l_Ayz * Au_d_zz;
    FieldReal term_zz = l_Axz * Au_d_xz + l_Ayz * Au_d_yz + l_Azz * Au_d_zz;

    FieldReal src_xx = alpn1 * (l_Rxx - EIGHT_PI * l_Sxx) - fxx_Lap_phy;
    FieldReal src_yy = alpn1 * (l_Ryy - EIGHT_PI * l_Syy) - fyy_Lap_phy;
    FieldReal src_zz = alpn1 * (l_Rzz - EIGHT_PI * l_Szz) - fzz_Lap_phy;
    FieldReal src_xy = alpn1 * (l_Rxy - EIGHT_PI * l_Sxy) - fxy_Lap_phy;
    FieldReal src_xz = alpn1 * (l_Rxz - EIGHT_PI * l_Sxz) - fxz_Lap_phy;
    FieldReal src_yz = alpn1 * (l_Ryz - EIGHT_PI * l_Syz) - fyz_Lap_phy;

    FieldReal trA2 = l_Axx * Aup_xx + l_Ayy * Aup_yy + l_Azz * Aup_zz + TWO * (l_Axy * Aup_xy + l_Axz * Aup_xz + l_Ayz * Aup_yz);
    FieldReal bssn_trace_aux = F2o3 * val_trK * val_trK - trA2 - F16_PI * rho[idx] + EIGHT_PI * S;
    FieldReal bssn_f_trace = -F1o3 * (trK_rhs_val + alpn1 * inv_chin1 * bssn_trace_aux);
    FieldReal escalar_f_trace = F1o3 * (gupxx * src_xx + gupyy * src_yy + gupzz * src_zz
        + TWO * (gupxy * src_xy + gupxz * src_xz + gupyz * src_yz));
    FieldReal f_trace = use_escalar_trace_free ? escalar_f_trace : bssn_f_trace;

    Axx_rhs[idx] = chin1 * (src_xx - l_gxx * f_trace) + alpn1 * (val_trK * l_Axx - TWO * term_xx) + TWO * (l_Axx * betaxx + l_Axy * betayx + l_Axz * betazx) - F2o3 * l_Axx * div_beta;
    Ayy_rhs[idx] = chin1 * (src_yy - l_gyy * f_trace) + alpn1 * (val_trK * l_Ayy - TWO * term_yy) + TWO * (l_Axy * betaxy + l_Ayy * betayy + l_Ayz * betazy) - F2o3 * l_Ayy * div_beta;
    Azz_rhs[idx] = chin1 * (src_zz - l_gzz * f_trace) + alpn1 * (val_trK * l_Azz - TWO * term_zz) + TWO * (l_Axz * betaxz + l_Ayz * betayz + l_Azz * betazz) - F2o3 * l_Azz * div_beta;
    Axy_rhs[idx] = chin1 * (src_xy - l_gxy * f_trace) + alpn1 * (val_trK * l_Axy - TWO * term_xy) + l_Axx * betaxy + l_Axz * betazy + l_Ayy * betayx + l_Ayz * betazx - l_Axy * betazz + F1o3 * l_Axy * div_beta;
    Ayz_rhs[idx] = chin1 * (src_yz - l_gyz * f_trace) + alpn1 * (val_trK * l_Ayz - TWO * term_yz) + l_Axy * betaxz + l_Ayy * betayz + l_Axz * betaxy + l_Azz * betazy - l_Ayz * betaxx + F1o3 * l_Ayz * div_beta;
    Axz_rhs[idx] = chin1 * (src_xz - l_gxz * f_trace) + alpn1 * (val_trK * l_Axz - TWO * term_xz) + l_Axx * betaxz + l_Axy * betayz + l_Ayz * betayx + l_Azz * betazx - l_Axz * betayy + F1o3 * l_Axz * div_beta;

    trK_rhs[idx] = -chin1 * trK_rhs_val + alpn1 * (F1o3 * val_trK * val_trK + trA2 + FOUR_PI * (rho[idx] + S));
}

template <typename FieldReal>
__launch_bounds__(256, 2)
__global__ void evolve_bssn_gauge_fields_kernel(
    const int ex0, const int ex1, const int ex2,
    const int symmetry,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    const FieldReal* __restrict__ chi, const FieldReal* __restrict__ trK,
    const FieldReal* __restrict__ dxx, const FieldReal* __restrict__ gxy, const FieldReal* __restrict__ gxz,
    const FieldReal* __restrict__ dyy, const FieldReal* __restrict__ gyz, const FieldReal* __restrict__ dzz,
    const FieldReal* __restrict__ Axx, const FieldReal* __restrict__ Axy, const FieldReal* __restrict__ Axz,
    const FieldReal* __restrict__ Ayy, const FieldReal* __restrict__ Ayz, const FieldReal* __restrict__ Azz,
    const FieldReal* __restrict__ Gamx, const FieldReal* __restrict__ Gamy, const FieldReal* __restrict__ Gamz,
    const FieldReal* __restrict__ Lap,
    const FieldReal* __restrict__ betax, const FieldReal* __restrict__ betay, const FieldReal* __restrict__ betaz,
    const FieldReal* __restrict__ dtSfx, const FieldReal* __restrict__ dtSfy, const FieldReal* __restrict__ dtSfz,
    const FieldReal* __restrict__ Sx, const FieldReal* __restrict__ Sy, const FieldReal* __restrict__ Sz,
    const FieldReal* __restrict__ Gamxxx_in, const FieldReal* __restrict__ Gamxxy_in, const FieldReal* __restrict__ Gamxxz_in,
    const FieldReal* __restrict__ Gamxyy_in, const FieldReal* __restrict__ Gamxyz_in, const FieldReal* __restrict__ Gamxzz_in,
    const FieldReal* __restrict__ Gamyxx_in, const FieldReal* __restrict__ Gamyxy_in, const FieldReal* __restrict__ Gamyxz_in,
    const FieldReal* __restrict__ Gamyyy_in, const FieldReal* __restrict__ Gamyyz_in, const FieldReal* __restrict__ Gamyzz_in,
    const FieldReal* __restrict__ Gamzxx_in, const FieldReal* __restrict__ Gamzxy_in, const FieldReal* __restrict__ Gamzxz_in,
    const FieldReal* __restrict__ Gamzyy_in, const FieldReal* __restrict__ Gamzyz_in, const FieldReal* __restrict__ Gamzzz_in,
    FieldReal* __restrict__ Gamx_rhs, FieldReal* __restrict__ Gamy_rhs, FieldReal* __restrict__ Gamz_rhs,
    FieldReal* __restrict__ Lap_rhs,
    FieldReal* __restrict__ betax_rhs, FieldReal* __restrict__ betay_rhs, FieldReal* __restrict__ betaz_rhs,
    FieldReal* __restrict__ dtSfx_rhs, FieldReal* __restrict__ dtSfy_rhs, FieldReal* __restrict__ dtSfz_rhs
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    int dims[3] = {ex0, ex1, ex2};
    FieldReal dX = X[1] - X[0];
    FieldReal dY = Y[1] - Y[0];
    FieldReal dZ = Z[1] - Z[0];
    StencilBounds bounds = make_stencil_bounds(ex0, ex1, ex2, symmetry, X, Y, Z, dX, dY, dZ, kBssnCenteredBoundaryRadius);

    int block_i = blockIdx.x * blockDim.x;
    int block_j = blockIdx.y * blockDim.y;
    int block_k = blockIdx.z * blockDim.z;

    BSSN_CENTERED_SMEM_TYPED(FieldReal, smem);

    FieldReal betaxx, betaxy, betaxz, bx_gxxx, bx_gxyx, bx_gxzx, bx_gyyx, bx_gyzx, bx_gzzx;
    load_field_to_smem(smem, betax, dims, block_i, block_j, block_k, ANTI, SYM, SYM);
    __syncthreads(); compute_all_derivs_smem(smem, &betaxx, &betaxy, &betaxz, &bx_gxxx, &bx_gxyx, &bx_gxzx, &bx_gyyx, &bx_gyzx, &bx_gzzx, dX, dY, dZ, i, j, k, bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax); __syncthreads();

    FieldReal betayx, betayy, betayz, by_gxxy, by_gxyy, by_gxzy, by_gyyy, by_gyzy, by_gzzy;
    load_field_to_smem(smem, betay, dims, block_i, block_j, block_k, SYM, ANTI, SYM);
    __syncthreads(); compute_all_derivs_smem(smem, &betayx, &betayy, &betayz, &by_gxxy, &by_gxyy, &by_gxzy, &by_gyyy, &by_gyzy, &by_gzzy, dX, dY, dZ, i, j, k, bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax); __syncthreads();

    FieldReal betazx, betazy, betazz, bz_gxxz, bz_gxyz, bz_gxzz, bz_gyyz, bz_gyzz, bz_gzzz;
    load_field_to_smem(smem, betaz, dims, block_i, block_j, block_k, SYM, SYM, ANTI);
    __syncthreads(); compute_all_derivs_smem(smem, &betazx, &betazy, &betazz, &bz_gxxz, &bz_gxyz, &bz_gxzz, &bz_gyyz, &bz_gyzz, &bz_gzzz, dX, dY, dZ, i, j, k, bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax); __syncthreads();

    FieldReal Lapx, Lapy, Lapz;
    load_field_to_smem(smem, Lap, dims, block_i, block_j, block_k, SYM, SYM, SYM);
    __syncthreads(); compute_first_derivs_smem(smem, &Lapx, &Lapy, &Lapz, dX, dY, dZ, i, j, k, bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax); __syncthreads();

    FieldReal chix, chiy, chiz;
    load_field_to_smem(smem, chi, dims, block_i, block_j, block_k, SYM, SYM, SYM);
    __syncthreads(); compute_first_derivs_smem(smem, &chix, &chiy, &chiz, dX, dY, dZ, i, j, k, bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax); __syncthreads();

    FieldReal Kx, Ky, Kz;
    load_field_to_smem(smem, trK, dims, block_i, block_j, block_k, SYM, SYM, SYM);
    __syncthreads(); compute_first_derivs_smem(smem, &Kx, &Ky, &Kz, dX, dY, dZ, i, j, k, bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax);

    if (i >= ex0 || j >= ex1 || k >= ex2) return;
    int idx = flatten_bssn_index(i, j, k, ex0, ex1);

    const FieldReal ONE = FieldReal(1.0);
    const FieldReal TWO = FieldReal(2.0);
    const FieldReal F1o3 = FieldReal(1.0) / FieldReal(3.0);
    const FieldReal F2o3 = FieldReal(2.0) / FieldReal(3.0);
    const FieldReal F3o2 = FieldReal(1.5);
    const FieldReal EIGHT_PI = FieldReal(8.0) * FieldReal(M_PI);
    const FieldReal FF = FieldReal(0.75);
    const FieldReal eta = FieldReal(2.0);

    FieldReal val_Lap = Lap[idx];
    FieldReal val_trK = trK[idx];
    FieldReal alpn1 = val_Lap + ONE;
    FieldReal chin1 = chi[idx] + ONE;
    FieldReal inv_chin1 = ONE / chin1;

    FieldReal l_gxx = dxx[idx] + ONE;
    FieldReal l_gxy = gxy[idx];
    FieldReal l_gxz = gxz[idx];
    FieldReal l_gyy = dyy[idx] + ONE;
    FieldReal l_gyz = gyz[idx];
    FieldReal l_gzz = dzz[idx] + ONE;
    FieldReal l_Axx = Axx[idx];
    FieldReal l_Axy = Axy[idx];
    FieldReal l_Axz = Axz[idx];
    FieldReal l_Ayy = Ayy[idx];
    FieldReal l_Ayz = Ayz[idx];
    FieldReal l_Azz = Azz[idx];

    FieldReal gupxx, gupxy, gupxz, gupyy, gupyz, gupzz;
    compute_inverse_metric_t<FieldReal>(l_gxx, l_gxy, l_gxz, l_gyy, l_gyz, l_gzz, &gupxx, &gupxy, &gupxz, &gupyy, &gupyz, &gupzz);

    FieldReal Au_d_xx, Au_d_xy, Au_d_xz, Au_d_yx, Au_d_yy, Au_d_yz, Au_d_zx, Au_d_zy, Au_d_zz;
    FieldReal Aup_xx, Aup_xy, Aup_xz, Aup_yy, Aup_yz, Aup_zz;
    compute_raised_extrinsic_curvature_t<FieldReal>(
        l_Axx, l_Axy, l_Axz, l_Ayy, l_Ayz, l_Azz,
        gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
        &Au_d_xx, &Au_d_xy, &Au_d_xz, &Au_d_yx, &Au_d_yy, &Au_d_yz, &Au_d_zx, &Au_d_zy, &Au_d_zz,
        &Aup_xx, &Aup_xy, &Aup_xz, &Aup_yy, &Aup_yz, &Aup_zz
    );

    FieldReal l_Gamxxx = Gamxxx_in[idx];
    FieldReal l_Gamxxy = Gamxxy_in[idx];
    FieldReal l_Gamxxz = Gamxxz_in[idx];
    FieldReal l_Gamxyy = Gamxyy_in[idx];
    FieldReal l_Gamxyz = Gamxyz_in[idx];
    FieldReal l_Gamxzz = Gamxzz_in[idx];
    FieldReal l_Gamyxx = Gamyxx_in[idx];
    FieldReal l_Gamyxy = Gamyxy_in[idx];
    FieldReal l_Gamyxz = Gamyxz_in[idx];
    FieldReal l_Gamyyy = Gamyyy_in[idx];
    FieldReal l_Gamyyz = Gamyyz_in[idx];
    FieldReal l_Gamyzz = Gamyzz_in[idx];
    FieldReal l_Gamzxx = Gamzxx_in[idx];
    FieldReal l_Gamzxy = Gamzxy_in[idx];
    FieldReal l_Gamzxz = Gamzxz_in[idx];
    FieldReal l_Gamzyy = Gamzyy_in[idx];
    FieldReal l_Gamzyz = Gamzyz_in[idx];
    FieldReal l_Gamzzz = Gamzzz_in[idx];

    FieldReal val_Gamx_rhs = -TWO * (Lapx * Aup_xx + Lapy * Aup_xy + Lapz * Aup_xz) +
        TWO * alpn1 * (-F3o2 * inv_chin1 * (chix * Aup_xx + chiy * Aup_xy + chiz * Aup_xz) - gupxx * (F2o3 * Kx + EIGHT_PI * Sx[idx]) - gupxy * (F2o3 * Ky + EIGHT_PI * Sy[idx]) - gupxz * (F2o3 * Kz + EIGHT_PI * Sz[idx]) +
        l_Gamxxx * Aup_xx + l_Gamxyy * Aup_yy + l_Gamxzz * Aup_zz + TWO * (l_Gamxxy * Aup_xy + l_Gamxxz * Aup_xz + l_Gamxyz * Aup_yz));

    FieldReal val_Gamy_rhs = -TWO * (Lapx * Aup_xy + Lapy * Aup_yy + Lapz * Aup_yz) +
        TWO * alpn1 * (-F3o2 * inv_chin1 * (chix * Aup_xy + chiy * Aup_yy + chiz * Aup_yz) - gupxy * (F2o3 * Kx + EIGHT_PI * Sx[idx]) - gupyy * (F2o3 * Ky + EIGHT_PI * Sy[idx]) - gupyz * (F2o3 * Kz + EIGHT_PI * Sz[idx]) +
        l_Gamyxx * Aup_xx + l_Gamyyy * Aup_yy + l_Gamyzz * Aup_zz + TWO * (l_Gamyxy * Aup_xy + l_Gamyxz * Aup_xz + l_Gamyyz * Aup_yz));

    FieldReal val_Gamz_rhs = -TWO * (Lapx * Aup_xz + Lapy * Aup_yz + Lapz * Aup_zz) +
        TWO * alpn1 * (-F3o2 * inv_chin1 * (chix * Aup_xz + chiy * Aup_yz + chiz * Aup_zz) - gupxz * (F2o3 * Kx + EIGHT_PI * Sx[idx]) - gupyz * (F2o3 * Ky + EIGHT_PI * Sy[idx]) - gupzz * (F2o3 * Kz + EIGHT_PI * Sz[idx]) +
        l_Gamzxx * Aup_xx + l_Gamzyy * Aup_yy + l_Gamzzz * Aup_zz + TWO * (l_Gamzxy * Aup_xy + l_Gamzxz * Aup_xz + l_Gamzyz * Aup_yz));

    FieldReal fxx_beta = bx_gxxx + by_gxyy + bz_gxzz;
    FieldReal fxy_beta = bx_gxyx + by_gyyy + bz_gyzz;
    FieldReal fxz_beta = bx_gxzx + by_gyzy + bz_gzzz;

    FieldReal div_beta = betaxx + betayy + betazz;
    FieldReal Gamxa = gupxx * l_Gamxxx + gupyy * l_Gamxyy + gupzz * l_Gamxzz + TWO * (gupxy * l_Gamxxy + gupxz * l_Gamxxz + gupyz * l_Gamxyz);
    FieldReal Gamya = gupxx * l_Gamyxx + gupyy * l_Gamyyy + gupzz * l_Gamyzz + TWO * (gupxy * l_Gamyxy + gupxz * l_Gamyxz + gupyz * l_Gamyyz);
    FieldReal Gamza = gupxx * l_Gamzxx + gupyy * l_Gamzyy + gupzz * l_Gamzzz + TWO * (gupxy * l_Gamzxy + gupxz * l_Gamzxz + gupyz * l_Gamzyz);

    val_Gamx_rhs += F2o3 * Gamxa * div_beta - (Gamxa * betaxx + Gamya * betaxy + Gamza * betaxz) + F1o3 * (gupxx * fxx_beta + gupxy * fxy_beta + gupxz * fxz_beta) + gupxx * bx_gxxx + gupyy * bx_gyyx + gupzz * bx_gzzx + TWO * (gupxy * bx_gxyx + gupxz * bx_gxzx + gupyz * bx_gyzx);
    val_Gamy_rhs += F2o3 * Gamya * div_beta - (Gamxa * betayx + Gamya * betayy + Gamza * betayz) + F1o3 * (gupxy * fxx_beta + gupyy * fxy_beta + gupyz * fxz_beta) + gupxx * by_gxxy + gupyy * by_gyyy + gupzz * by_gzzy + TWO * (gupxy * by_gxyy + gupxz * by_gxzy + gupyz * by_gyzy);
    val_Gamz_rhs += F2o3 * Gamza * div_beta - (Gamxa * betazx + Gamya * betazy + Gamza * betazz) + F1o3 * (gupxz * fxx_beta + gupyz * fxy_beta + gupzz * fxz_beta) + gupxx * bz_gxxz + gupyy * bz_gyyz + gupzz * bz_gzzz + TWO * (gupxy * bz_gxyz + gupxz * bz_gxzz + gupyz * bz_gyzz);

    FieldReal l_dtSfx = dtSfx[idx];
    FieldReal l_dtSfy = dtSfy[idx];
    FieldReal l_dtSfz = dtSfz[idx];
    Lap_rhs[idx] = -TWO * alpn1 * val_trK;

#if (GAUGE == 0)
    betax_rhs[idx] = FF * l_dtSfx;
    betay_rhs[idx] = FF * l_dtSfy;
    betaz_rhs[idx] = FF * l_dtSfz;
    dtSfx_rhs[idx] = val_Gamx_rhs - eta * l_dtSfx;
    dtSfy_rhs[idx] = val_Gamy_rhs - eta * l_dtSfy;
    dtSfz_rhs[idx] = val_Gamz_rhs - eta * l_dtSfz;
#elif (GAUGE == 1)
    betax_rhs[idx] = Gamx[idx] - eta * betax[idx];
    betay_rhs[idx] = Gamy[idx] - eta * betay[idx];
    betaz_rhs[idx] = Gamz[idx] - eta * betaz[idx];
    dtSfx_rhs[idx] = FieldReal(0.0);
    dtSfy_rhs[idx] = FieldReal(0.0);
    dtSfz_rhs[idx] = FieldReal(0.0);
#elif (GAUGE == 2 || GAUGE == 3 || GAUGE == 4 || GAUGE == 5)
    FieldReal reta = gupxx * chix * chix + gupyy * chiy * chiy + gupzz * chiz * chiz
                + TWO * (gupxy * chix * chiy + gupxz * chix * chiz + gupyz * chiy * chiz);
#if (GAUGE == 2 || GAUGE == 4)
    FieldReal reta_denom = ONE - sqrt(chin1);
#else
    FieldReal reta_denom = ONE - chin1;
#endif
    reta = FieldReal(1.31) / FieldReal(2.0) * sqrt(reta / chin1) / (reta_denom * reta_denom);
#if (GAUGE == 2 || GAUGE == 3)
    betax_rhs[idx] = FF * l_dtSfx;
    betay_rhs[idx] = FF * l_dtSfy;
    betaz_rhs[idx] = FF * l_dtSfz;
    dtSfx_rhs[idx] = val_Gamx_rhs - reta * l_dtSfx;
    dtSfy_rhs[idx] = val_Gamy_rhs - reta * l_dtSfy;
    dtSfz_rhs[idx] = val_Gamz_rhs - reta * l_dtSfz;
#else
    betax_rhs[idx] = FF * Gamx[idx] - reta * betax[idx];
    betay_rhs[idx] = FF * Gamy[idx] - reta * betay[idx];
    betaz_rhs[idx] = FF * Gamz[idx] - reta * betaz[idx];
    dtSfx_rhs[idx] = FieldReal(0.0);
    dtSfy_rhs[idx] = FieldReal(0.0);
    dtSfz_rhs[idx] = FieldReal(0.0);
#endif
#elif (GAUGE == 6 || GAUGE == 7)
    FieldReal reta = eta;
    if (d_bssn_gauge_bh_count == 2) {
        const FieldReal total_mass = d_bssn_gauge_mass[0] + d_bssn_gauge_mass[1];
        const FieldReal A = TWO / total_mass;
        const FieldReal w1 = 12.0;
        const FieldReal w2 = w1;
        const FieldReal C1 = ONE / d_bssn_gauge_mass[0] - A;
        const FieldReal C2 = ONE / d_bssn_gauge_mass[1] - A;
        const FieldReal x = X[i];
        const FieldReal y = Y[j];
        const FieldReal z = Z[k];
        const FieldReal dx12 = d_bssn_gauge_pbh[0] - d_bssn_gauge_pbh[3];
        const FieldReal dy12 = d_bssn_gauge_pbh[1] - d_bssn_gauge_pbh[4];
        const FieldReal dz12 = d_bssn_gauge_pbh[2] - d_bssn_gauge_pbh[5];
        const FieldReal sep2 = dx12 * dx12 + dy12 * dy12 + dz12 * dz12;
        const FieldReal dx1 = d_bssn_gauge_pbh[0] - x;
        const FieldReal dy1 = d_bssn_gauge_pbh[1] - y;
        const FieldReal dz1 = d_bssn_gauge_pbh[2] - z;
        const FieldReal dx2 = d_bssn_gauge_pbh[3] - x;
        const FieldReal dy2 = d_bssn_gauge_pbh[4] - y;
        const FieldReal dz2 = d_bssn_gauge_pbh[5] - z;
        const FieldReal r1 = (dx1 * dx1 + dy1 * dy1 + dz1 * dz1) / sep2;
        const FieldReal r2 = (dx2 * dx2 + dy2 * dy2 + dz2 * dz2) / sep2;
#if (GAUGE == 6)
        reta = A + C1 / (ONE + w1 * r1) + C2 / (ONE + w2 * r2);
#else
        reta = A + C1 * exp(-w1 * r1) + C2 * exp(-w2 * r2);
#endif
    }
    betax_rhs[idx] = FF * l_dtSfx;
    betay_rhs[idx] = FF * l_dtSfy;
    betaz_rhs[idx] = FF * l_dtSfz;
    dtSfx_rhs[idx] = val_Gamx_rhs - reta * l_dtSfx;
    dtSfy_rhs[idx] = val_Gamy_rhs - reta * l_dtSfy;
    dtSfz_rhs[idx] = val_Gamz_rhs - reta * l_dtSfz;
#else
#error Unsupported BSSN GPU gauge; expected GAUGE == 0..7.
#endif

    Gamx_rhs[idx] = val_Gamx_rhs;
    Gamy_rhs[idx] = val_Gamy_rhs;
    Gamz_rhs[idx] = val_Gamz_rhs;
}

template <typename Real>
__launch_bounds__(256, 2)
__global__ void apply_transport_and_dissipation_kernel(
    const int ex0, const int ex1, const int ex2,
    const int symmetry, const int lev, const double eps,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    const Real* __restrict__ betax, const Real* __restrict__ betay, const Real* __restrict__ betaz,
    const Real* __restrict__ dxx, const Real* __restrict__ gxy, const Real* __restrict__ gxz,
    const Real* __restrict__ dyy, const Real* __restrict__ gyz, const Real* __restrict__ dzz,
    const Real* __restrict__ Axx, const Real* __restrict__ Axy, const Real* __restrict__ Axz,
    const Real* __restrict__ Ayy, const Real* __restrict__ Ayz, const Real* __restrict__ Azz,
    const Real* __restrict__ chi, const Real* __restrict__ trK,
    const Real* __restrict__ Gamx, const Real* __restrict__ Gamy, const Real* __restrict__ Gamz,
    const Real* __restrict__ Lap,
    const Real* __restrict__ dtSfx, const Real* __restrict__ dtSfy, const Real* __restrict__ dtSfz,
    Real* __restrict__ gxx_rhs, Real* __restrict__ gxy_rhs, Real* __restrict__ gxz_rhs,
    Real* __restrict__ gyy_rhs, Real* __restrict__ gyz_rhs, Real* __restrict__ gzz_rhs,
    Real* __restrict__ Axx_rhs, Real* __restrict__ Axy_rhs, Real* __restrict__ Axz_rhs,
    Real* __restrict__ Ayy_rhs, Real* __restrict__ Ayz_rhs, Real* __restrict__ Azz_rhs,
    Real* __restrict__ chi_rhs, Real* __restrict__ trK_rhs,
    Real* __restrict__ Gamx_rhs, Real* __restrict__ Gamy_rhs, Real* __restrict__ Gamz_rhs,
    Real* __restrict__ Lap_rhs,
    Real* __restrict__ betax_rhs, Real* __restrict__ betay_rhs, Real* __restrict__ betaz_rhs,
    Real* __restrict__ dtSfx_rhs, Real* __restrict__ dtSfy_rhs, Real* __restrict__ dtSfz_rhs
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    int dims[3] = {ex0, ex1, ex2};
    int block_i = blockIdx.x * blockDim.x;
    int block_j = blockIdx.y * blockDim.y;
    int block_k = blockIdx.z * blockDim.z;

    const double dX_native = X[1] - X[0];
    const double dY_native = Y[1] - Y[0];
    const double dZ_native = Z[1] - Z[0];
    Real dX = static_cast<Real>(dX_native);
    Real dY = static_cast<Real>(dY_native);
    Real dZ = static_cast<Real>(dZ_native);
#if ghost_width == 5
    Real d840dx = Real(1.0) / Real(840.0) / dX;
    Real d840dy = Real(1.0) / Real(840.0) / dY;
    Real d840dz = Real(1.0) / Real(840.0) / dZ;
    Real d60dx = Real(1.0) / Real(60.0) / dX;
    Real d60dy = Real(1.0) / Real(60.0) / dY;
    Real d60dz = Real(1.0) / Real(60.0) / dZ;
    Real d12dx = Real(1.0) / Real(12.0) / dX;
    Real d12dy = Real(1.0) / Real(12.0) / dY;
    Real d12dz = Real(1.0) / Real(12.0) / dZ;
    Real d2dx = Real(1.0) / Real(2.0) / dX;
    Real d2dy = Real(1.0) / Real(2.0) / dY;
    Real d2dz = Real(1.0) / Real(2.0) / dZ;
#elif ghost_width == 4
    Real d60dx = Real(1.0) / Real(60.0) / dX;
    Real d60dy = Real(1.0) / Real(60.0) / dY;
    Real d60dz = Real(1.0) / Real(60.0) / dZ;
    Real d12dx = Real(1.0) / Real(12.0) / dX;
    Real d12dy = Real(1.0) / Real(12.0) / dY;
    Real d12dz = Real(1.0) / Real(12.0) / dZ;
    Real d2dx = Real(1.0) / Real(2.0) / dX;
    Real d2dy = Real(1.0) / Real(2.0) / dY;
    Real d2dz = Real(1.0) / Real(2.0) / dZ;
#elif ghost_width == 3
    Real d12dx = Real(1.0) / Real(12.0) / dX;
    Real d12dy = Real(1.0) / Real(12.0) / dY;
    Real d12dz = Real(1.0) / Real(12.0) / dZ;
#elif ghost_width == 2
    Real d2dx = Real(1.0) / Real(2.0) / dX;
    Real d2dy = Real(1.0) / Real(2.0) / dY;
    Real d2dz = Real(1.0) / Real(2.0) / dZ;
#else
#error Unsupported ghost_width in apply_transport_and_dissipation_kernel; expected ghost_width == 2, 3, 4, or 5.
#endif
    StencilBounds bounds = make_stencil_bounds(ex0, ex1, ex2, symmetry, X, Y, Z, dX_native, dY_native, dZ_native, kBssnTransportBoundaryRadius);
    const int imin = bounds.imin;
    const int jmin = bounds.jmin;
    const int kmin = bounds.kmin;
    const int imax = bounds.imax;
    const int jmax = bounds.jmax;
    const int kmax = bounds.kmax;

    Real vx = Real(0.0), vy = Real(0.0), vz = Real(0.0);
    int idx = 0;
    bool is_valid_point = (i < imax && j < jmax && k < kmax);

    if (is_valid_point) {
        idx = flatten_bssn_index(i, j, k, ex0, ex1);
        vx = static_cast<Real>(__ldg(&betax[idx]));
        vy = static_cast<Real>(__ldg(&betay[idx]));
        vz = static_cast<Real>(__ldg(&betaz[idx]));
    }

#if ghost_width == 4 || ghost_width == 5
    Real Wx[11] = {0}, Wy[11] = {0}, Wz[11] = {0};
#else
    Real Wx[7] = {0}, Wy[7] = {0}, Wz[7] = {0};
#endif

    if (is_valid_point) {
#if ghost_width == 5
        accumulate_fd8_transport_axis(vx, d840dx, d60dx, d12dx, d2dx, i, imin, imax, Wx);
        accumulate_fd8_transport_axis(vy, d840dy, d60dy, d12dy, d2dy, j, jmin, jmax, Wy);
        accumulate_fd8_transport_axis(vz, d840dz, d60dz, d12dz, d2dz, k, kmin, kmax, Wz);

        if (eps > 0.0 && i - 5 >= imin && i + 5 <= imax && j - 5 >= jmin && j + 5 <= jmax && k - 5 >= kmin && k + 5 <= kmax) {
            Real edx = static_cast<Real>(eps) / Real(1024.0) / dX;
            Real edy = static_cast<Real>(eps) / Real(1024.0) / dY;
            Real edz = static_cast<Real>(eps) / Real(1024.0) / dZ;
            Wx[0] += edx;             Wy[0] += edy;             Wz[0] += edz;
            Wx[1] -= Real(10.0) * edx;  Wy[1] -= Real(10.0) * edy;  Wz[1] -= Real(10.0) * edz;
            Wx[2] += Real(45.0) * edx;  Wy[2] += Real(45.0) * edy;  Wz[2] += Real(45.0) * edz;
            Wx[3] -= Real(120.0) * edx; Wy[3] -= Real(120.0) * edy; Wz[3] -= Real(120.0) * edz;
            Wx[4] += Real(210.0) * edx; Wy[4] += Real(210.0) * edy; Wz[4] += Real(210.0) * edz;
            Wx[5] -= Real(252.0) * edx; Wy[5] -= Real(252.0) * edy; Wz[5] -= Real(252.0) * edz;
            Wx[6] += Real(210.0) * edx; Wy[6] += Real(210.0) * edy; Wz[6] += Real(210.0) * edz;
            Wx[7] -= Real(120.0) * edx; Wy[7] -= Real(120.0) * edy; Wz[7] -= Real(120.0) * edz;
            Wx[8] += Real(45.0) * edx;  Wy[8] += Real(45.0) * edy;  Wz[8] += Real(45.0) * edz;
            Wx[9] -= Real(10.0) * edx;  Wy[9] -= Real(10.0) * edy;  Wz[9] -= Real(10.0) * edz;
            Wx[10] += edx;        Wy[10] += edy;        Wz[10] += edz;
        }
#elif ghost_width == 4
        accumulate_fd6_transport_axis(vx, d60dx, d12dx, d2dx, i, imin, imax, Wx);
        accumulate_fd6_transport_axis(vy, d60dy, d12dy, d2dy, j, jmin, jmax, Wy);
        accumulate_fd6_transport_axis(vz, d60dz, d12dz, d2dz, k, kmin, kmax, Wz);

        if (eps > 0.0 && i - 4 >= imin && i + 4 <= imax && j - 4 >= jmin && j + 4 <= jmax && k - 4 >= kmin && k + 4 <= kmax) {
            Real edx = static_cast<Real>(eps) / Real(256.0) / dX;
            Real edy = static_cast<Real>(eps) / Real(256.0) / dY;
            Real edz = static_cast<Real>(eps) / Real(256.0) / dZ;
            Wx[1] -= edx;        Wy[1] -= edy;        Wz[1] -= edz;
            Wx[2] += Real(8.0) * edx;  Wy[2] += Real(8.0) * edy;  Wz[2] += Real(8.0) * edz;
            Wx[3] -= Real(28.0) * edx; Wy[3] -= Real(28.0) * edy; Wz[3] -= Real(28.0) * edz;
            Wx[4] += Real(56.0) * edx; Wy[4] += Real(56.0) * edy; Wz[4] += Real(56.0) * edz;
            Wx[5] -= Real(70.0) * edx; Wy[5] -= Real(70.0) * edy; Wz[5] -= Real(70.0) * edz;
            Wx[6] += Real(56.0) * edx; Wy[6] += Real(56.0) * edy; Wz[6] += Real(56.0) * edz;
            Wx[7] -= Real(28.0) * edx; Wy[7] -= Real(28.0) * edy; Wz[7] -= Real(28.0) * edz;
            Wx[8] += Real(8.0) * edx;  Wy[8] += Real(8.0) * edy;  Wz[8] += Real(8.0) * edz;
            Wx[9] -= edx;        Wy[9] -= edy;        Wz[9] -= edz;
        }
#elif ghost_width == 3
        Real v_dx = vx * d12dx;
        if (vx > Real(0.0)) {
            if (i + 3 <= imax) { Wx[2]-=3*v_dx; Wx[3]-=10*v_dx; Wx[4]+=18*v_dx; Wx[5]-=6*v_dx; Wx[6]+=v_dx; }
            else if (i + 2 <= imax) { Wx[1]+=v_dx; Wx[2]-=8*v_dx; Wx[4]+=8*v_dx; Wx[5]-=v_dx; }
            else if (i + 1 <= imax) { Wx[4]+=3*v_dx; Wx[3]+=10*v_dx; Wx[2]-=18*v_dx; Wx[1]+=6*v_dx; Wx[0]-=v_dx; }
        } else if (vx < Real(0.0)) {
            if (i - 3 >= imin) { Wx[4]+=3*v_dx; Wx[3]+=10*v_dx; Wx[2]-=18*v_dx; Wx[1]+=6*v_dx; Wx[0]-=v_dx; }
            else if (i - 2 >= imin) { Wx[1]+=v_dx; Wx[2]-=8*v_dx; Wx[4]+=8*v_dx; Wx[5]-=v_dx; }
            else if (i - 1 >= imin) { Wx[2]-=3*v_dx; Wx[3]-=10*v_dx; Wx[4]+=18*v_dx; Wx[5]-=6*v_dx; Wx[6]+=v_dx; }
        }

        Real v_dy = vy * d12dy;
        if (vy > Real(0.0)) {
            if (j + 3 <= jmax) { Wy[2]-=3*v_dy; Wy[3]-=10*v_dy; Wy[4]+=18*v_dy; Wy[5]-=6*v_dy; Wy[6]+=v_dy; }
            else if (j + 2 <= jmax) { Wy[1]+=v_dy; Wy[2]-=8*v_dy; Wy[4]+=8*v_dy; Wy[5]-=v_dy; }
            else if (j + 1 <= jmax) { Wy[4]+=3*v_dy; Wy[3]+=10*v_dy; Wy[2]-=18*v_dy; Wy[1]+=6*v_dy; Wy[0]-=v_dy; }
        } else if (vy < Real(0.0)) {
            if (j - 3 >= jmin) { Wy[4]+=3*v_dy; Wy[3]+=10*v_dy; Wy[2]-=18*v_dy; Wy[1]+=6*v_dy; Wy[0]-=v_dy; }
            else if (j - 2 >= jmin) { Wy[1]+=v_dy; Wy[2]-=8*v_dy; Wy[4]+=8*v_dy; Wy[5]-=v_dy; }
            else if (j - 1 >= jmin) { Wy[2]-=3*v_dy; Wy[3]-=10*v_dy; Wy[4]+=18*v_dy; Wy[5]-=6*v_dy; Wy[6]+=v_dy; }
        }

        Real v_dz = vz * d12dz;
        if (vz > Real(0.0)) {
            if (k + 3 <= kmax) { Wz[2]-=3*v_dz; Wz[3]-=10*v_dz; Wz[4]+=18*v_dz; Wz[5]-=6*v_dz; Wz[6]+=v_dz; }
            else if (k + 2 <= kmax) { Wz[1]+=v_dz; Wz[2]-=8*v_dz; Wz[4]+=8*v_dz; Wz[5]-=v_dz; }
            else if (k + 1 <= kmax) { Wz[4]+=3*v_dz; Wz[3]+=10*v_dz; Wz[2]-=18*v_dz; Wz[1]+=6*v_dz; Wz[0]-=v_dz; }
        } else if (vz < Real(0.0)) {
            if (k - 3 >= kmin) { Wz[4]+=3*v_dz; Wz[3]+=10*v_dz; Wz[2]-=18*v_dz; Wz[1]+=6*v_dz; Wz[0]-=v_dz; }
            else if (k - 2 >= kmin) { Wz[1]+=v_dz; Wz[2]-=8*v_dz; Wz[4]+=8*v_dz; Wz[5]-=v_dz; }
            else if (k - 1 >= kmin) { Wz[2]-=3*v_dz; Wz[3]-=10*v_dz; Wz[4]+=18*v_dz; Wz[5]-=6*v_dz; Wz[6]+=v_dz; }
        }

        if (eps > 0.0 && i - 3 >= imin && i + 3 <= imax && j - 3 >= jmin && j + 3 <= jmax && k - 3 >= kmin && k + 3 <= kmax) {
            Real edx = static_cast<Real>(eps) / Real(64.0) / dX; 
            Real edy = static_cast<Real>(eps) / Real(64.0) / dY; 
            Real edz = static_cast<Real>(eps) / Real(64.0) / dZ;
            Wx[0]+=edx; Wx[1]-=6*edx; Wx[2]+=15*edx; Wx[3]-=20*edx; Wx[4]+=15*edx; Wx[5]-=6*edx; Wx[6]+=edx;
            Wy[0]+=edy; Wy[1]-=6*edy; Wy[2]+=15*edy; Wy[3]-=20*edy; Wy[4]+=15*edy; Wy[5]-=6*edy; Wy[6]+=edy;
            Wz[0]+=edz; Wz[1]-=6*edz; Wz[2]+=15*edz; Wz[3]-=20*edz; Wz[4]+=15*edz; Wz[5]-=6*edz; Wz[6]+=edz;
        }
#elif ghost_width == 2
        Real v_dx = vx * d2dx;
        if (vx > Real(0.0)) {
            if (i + 2 <= imax && i >= imin) { Wx[3] -= 3.0 * v_dx; Wx[4] += 4.0 * v_dx; Wx[5] -= v_dx; }
            else if (i + 1 <= imax && i >= imin) { Wx[3] -= v_dx; Wx[4] += v_dx; }
        } else if (vx < Real(0.0)) {
            if (i - 2 >= imin && i <= imax) { Wx[1] += v_dx; Wx[2] -= 4.0 * v_dx; Wx[3] += 3.0 * v_dx; }
            else if (i - 1 >= imin && i <= imax) { Wx[2] -= v_dx; Wx[3] += v_dx; }
        }

        Real v_dy = vy * d2dy;
        if (vy > Real(0.0)) {
            if (j + 2 <= jmax && j >= jmin) { Wy[3] -= 3.0 * v_dy; Wy[4] += 4.0 * v_dy; Wy[5] -= v_dy; }
            else if (j + 1 <= jmax && j >= jmin) { Wy[3] -= v_dy; Wy[4] += v_dy; }
        } else if (vy < Real(0.0)) {
            if (j - 2 >= jmin && j <= jmax) { Wy[1] += v_dy; Wy[2] -= 4.0 * v_dy; Wy[3] += 3.0 * v_dy; }
            else if (j - 1 >= jmin && j <= jmax) { Wy[2] -= v_dy; Wy[3] += v_dy; }
        }

        Real v_dz = vz * d2dz;
        if (vz > Real(0.0)) {
            if (k + 2 <= kmax && k >= kmin) { Wz[3] -= 3.0 * v_dz; Wz[4] += 4.0 * v_dz; Wz[5] -= v_dz; }
            else if (k + 1 <= kmax && k >= kmin) { Wz[3] -= v_dz; Wz[4] += v_dz; }
        } else if (vz < Real(0.0)) {
            if (k - 2 >= kmin && k <= kmax) { Wz[1] += v_dz; Wz[2] -= 4.0 * v_dz; Wz[3] += 3.0 * v_dz; }
            else if (k - 1 >= kmin && k <= kmax) { Wz[2] -= v_dz; Wz[3] += v_dz; }
        }

        if (eps > 0.0 && i - 2 >= imin && i + 2 <= imax && j - 2 >= jmin && j + 2 <= jmax && k - 2 >= kmin && k + 2 <= kmax) {
            Real edx = static_cast<Real>(eps) / Real(16.0) / dX;
            Real edy = static_cast<Real>(eps) / Real(16.0) / dY;
            Real edz = static_cast<Real>(eps) / Real(16.0) / dZ;
            Wx[1] -= edx; Wy[1] -= edy; Wz[1] -= edz;
            Wx[2] += 4.0 * edx; Wy[2] += 4.0 * edy; Wz[2] += 4.0 * edz;
            Wx[3] -= 6.0 * edx; Wy[3] -= 6.0 * edy; Wz[3] -= 6.0 * edz;
            Wx[4] += 4.0 * edx; Wy[4] += 4.0 * edy; Wz[4] += 4.0 * edz;
            Wx[5] -= edx; Wy[5] -= edy; Wz[5] -= edz;
        }
#else
#error Unsupported ghost_width in apply_transport_and_dissipation_kernel; expected ghost_width == 2, 3, 4, or 5.
#endif
    }

#if ghost_width == 4 || ghost_width == 5
    Real W0_center = Wx[5] + Wy[5] + Wz[5];
#else
    Real W0_center = Wx[3] + Wy[3] + Wz[3];
#endif

    BSSN_TRANSPORT_SMEM_TYPED(Real, smem);

    auto process_field = [&](const Real* f, Real* f_rhs, double SYM1, double SYM2, double SYM3) {
        load_field_to_smem_rad3_mixed(smem, f, dims, block_i, block_j, block_k, SYM1, SYM2, SYM3);
        __syncthreads(); 

        if (is_valid_point) {
            int tx = threadIdx.x + kBssnTransportLoadRadius;
            int ty = threadIdx.y + kBssnTransportLoadRadius;
            int tz = threadIdx.z + kBssnTransportLoadRadius;

#if ghost_width == 4 || ghost_width == 5
            Real rhs_add = W0_center * smem[tz][ty][tx] +
                Wx[0]*smem[tz][ty][tx-5] + Wx[1]*smem[tz][ty][tx-4] + Wx[2]*smem[tz][ty][tx-3] + Wx[3]*smem[tz][ty][tx-2] + Wx[4]*smem[tz][ty][tx-1] +
                Wx[6]*smem[tz][ty][tx+1] + Wx[7]*smem[tz][ty][tx+2] + Wx[8]*smem[tz][ty][tx+3] + Wx[9]*smem[tz][ty][tx+4] + Wx[10]*smem[tz][ty][tx+5] +
                Wy[0]*smem[tz][ty-5][tx] + Wy[1]*smem[tz][ty-4][tx] + Wy[2]*smem[tz][ty-3][tx] + Wy[3]*smem[tz][ty-2][tx] + Wy[4]*smem[tz][ty-1][tx] +
                Wy[6]*smem[tz][ty+1][tx] + Wy[7]*smem[tz][ty+2][tx] + Wy[8]*smem[tz][ty+3][tx] + Wy[9]*smem[tz][ty+4][tx] + Wy[10]*smem[tz][ty+5][tx] +
                Wz[0]*smem[tz-5][ty][tx] + Wz[1]*smem[tz-4][ty][tx] + Wz[2]*smem[tz-3][ty][tx] + Wz[3]*smem[tz-2][ty][tx] + Wz[4]*smem[tz-1][ty][tx] +
                Wz[6]*smem[tz+1][ty][tx] + Wz[7]*smem[tz+2][ty][tx] + Wz[8]*smem[tz+3][ty][tx] + Wz[9]*smem[tz+4][ty][tx] + Wz[10]*smem[tz+5][ty][tx];
#elif ghost_width == 3
            Real rhs_add = W0_center * smem[tz][ty][tx] +
                Wx[0]*smem[tz][ty][tx-3] + Wx[1]*smem[tz][ty][tx-2] + Wx[2]*smem[tz][ty][tx-1] + Wx[4]*smem[tz][ty][tx+1] + Wx[5]*smem[tz][ty][tx+2] + Wx[6]*smem[tz][ty][tx+3] +
                Wy[0]*smem[tz][ty-3][tx] + Wy[1]*smem[tz][ty-2][tx] + Wy[2]*smem[tz][ty-1][tx] + Wy[4]*smem[tz][ty+1][tx] + Wy[5]*smem[tz][ty+2][tx] + Wy[6]*smem[tz][ty+3][tx] +
                Wz[0]*smem[tz-3][ty][tx] + Wz[1]*smem[tz-2][ty][tx] + Wz[2]*smem[tz-1][ty][tx] + Wz[4]*smem[tz+1][ty][tx] + Wz[5]*smem[tz+2][ty][tx] + Wz[6]*smem[tz+3][ty][tx];
#elif ghost_width == 2
            Real rhs_add = W0_center * smem[tz][ty][tx] +
                Wx[1]*smem[tz][ty][tx-2] + Wx[2]*smem[tz][ty][tx-1] + Wx[4]*smem[tz][ty][tx+1] + Wx[5]*smem[tz][ty][tx+2] +
                Wy[1]*smem[tz][ty-2][tx] + Wy[2]*smem[tz][ty-1][tx] + Wy[4]*smem[tz][ty+1][tx] + Wy[5]*smem[tz][ty+2][tx] +
                Wz[1]*smem[tz-2][ty][tx] + Wz[2]*smem[tz-1][ty][tx] + Wz[4]*smem[tz+1][ty][tx] + Wz[5]*smem[tz+2][ty][tx];
#else
#error Unsupported ghost_width in transport shared-memory accumulation; expected ghost_width == 2, 3, 4, or 5.
#endif

            f_rhs[idx] += rhs_add;
        }
        __syncthreads(); 
    };

    process_field(dxx, gxx_rhs, SYM, SYM, SYM);
    process_field(gxy, gxy_rhs, ANTI, ANTI, SYM);
    process_field(gxz, gxz_rhs, ANTI, SYM, ANTI);
    process_field(dyy, gyy_rhs, SYM, SYM, SYM);
    process_field(gyz, gyz_rhs, SYM, ANTI, ANTI);
    process_field(dzz, gzz_rhs, SYM, SYM, SYM);
    
    process_field(Axx, Axx_rhs, SYM, SYM, SYM);
    process_field(Axy, Axy_rhs, ANTI, ANTI, SYM);
    process_field(Axz, Axz_rhs, ANTI, SYM, ANTI);
    process_field(Ayy, Ayy_rhs, SYM, SYM, SYM);
    process_field(Ayz, Ayz_rhs, SYM, ANTI, ANTI);
    process_field(Azz, Azz_rhs, SYM, SYM, SYM);

    process_field(chi, chi_rhs, SYM, SYM, SYM);
    process_field(trK, trK_rhs, SYM, SYM, SYM);
    
    process_field(Gamx, Gamx_rhs, ANTI, SYM, SYM);
    process_field(Gamy, Gamy_rhs, SYM, ANTI, SYM);
    process_field(Gamz, Gamz_rhs, SYM, SYM, ANTI);
    
    process_field(Lap, Lap_rhs, SYM, SYM, SYM);
    
    process_field(betax, betax_rhs, ANTI, SYM, SYM);
    process_field(betay, betay_rhs, SYM, ANTI, SYM);
    process_field(betaz, betaz_rhs, SYM, SYM, ANTI);

#if (GAUGE == 0 || GAUGE == 2 || GAUGE == 3 || GAUGE == 6 || GAUGE == 7)
    process_field(dtSfx, dtSfx_rhs, ANTI, SYM, SYM);
    process_field(dtSfy, dtSfy_rhs, SYM, ANTI, SYM);
    process_field(dtSfz, dtSfz_rhs, SYM, SYM, ANTI);
#elif (ABEtype == 2)
    process_field(dtSfx, dtSfx_rhs, ANTI, SYM, SYM);
    process_field(dtSfy, dtSfy_rhs, SYM, ANTI, SYM);
    process_field(dtSfz, dtSfz_rhs, SYM, SYM, ANTI);
#endif
}

template <typename FieldReal>
__launch_bounds__(256, 2)
__global__ void evaluate_bssn_constraints_kernel(
    const int ex0, const int ex1, const int ex2,
    const int symmetry, const int lev,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    const FieldReal* __restrict__ chi, const FieldReal* __restrict__ trK,
    const FieldReal* __restrict__ dxx, const FieldReal* __restrict__ gxy, const FieldReal* __restrict__ gxz,
    const FieldReal* __restrict__ dyy, const FieldReal* __restrict__ gyz, const FieldReal* __restrict__ dzz,
    const FieldReal* __restrict__ Axx, const FieldReal* __restrict__ Axy, const FieldReal* __restrict__ Axz,
    const FieldReal* __restrict__ Ayy, const FieldReal* __restrict__ Ayz, const FieldReal* __restrict__ Azz,
    const FieldReal* __restrict__ Gamx, const FieldReal* __restrict__ Gamy, const FieldReal* __restrict__ Gamz,
    const FieldReal* __restrict__ rho,
    const FieldReal* __restrict__ Sx, const FieldReal* __restrict__ Sy, const FieldReal* __restrict__ Sz,
    const FieldReal* __restrict__ Gamxxx_in, const FieldReal* __restrict__ Gamxxy_in, const FieldReal* __restrict__ Gamxxz_in,
    const FieldReal* __restrict__ Gamxyy_in, const FieldReal* __restrict__ Gamxyz_in, const FieldReal* __restrict__ Gamxzz_in,
    const FieldReal* __restrict__ Gamyxx_in, const FieldReal* __restrict__ Gamyxy_in, const FieldReal* __restrict__ Gamyxz_in,
    const FieldReal* __restrict__ Gamyyy_in, const FieldReal* __restrict__ Gamyyz_in, const FieldReal* __restrict__ Gamyzz_in,
    const FieldReal* __restrict__ Gamzxx_in, const FieldReal* __restrict__ Gamzxy_in, const FieldReal* __restrict__ Gamzxz_in,
    const FieldReal* __restrict__ Gamzyy_in, const FieldReal* __restrict__ Gamzyz_in, const FieldReal* __restrict__ Gamzzz_in,
    const FieldReal* __restrict__ Rxx_in, const FieldReal* __restrict__ Rxy_in, const FieldReal* __restrict__ Rxz_in,
    const FieldReal* __restrict__ Ryy_in, const FieldReal* __restrict__ Ryz_in, const FieldReal* __restrict__ Rzz_in,
    FieldReal* __restrict__ ham_Res,
    FieldReal* __restrict__ movx_Res, FieldReal* __restrict__ movy_Res, FieldReal* __restrict__ movz_Res,
    FieldReal* __restrict__ Gmx_Res, FieldReal* __restrict__ Gmy_Res, FieldReal* __restrict__ Gmz_Res
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    int dims[3] = {ex0, ex1, ex2};
    FieldReal dX = X[1] - X[0]; FieldReal dY = Y[1] - Y[0]; FieldReal dZ = Z[1] - Z[0];
    StencilBounds bounds = make_stencil_bounds(ex0, ex1, ex2, symmetry, X, Y, Z, dX, dY, dZ, kBssnCenteredBoundaryRadius);

    int block_i = blockIdx.x * blockDim.x;
    int block_j = blockIdx.y * blockDim.y;
    int block_k = blockIdx.z * blockDim.z;

    BSSN_CENTERED_SMEM_TYPED(FieldReal, smem);

    FieldReal d_Axx_x, d_Axx_y, d_Axx_z;
    load_field_to_smem(smem, Axx, dims, block_i, block_j, block_k, SYM, SYM, SYM);
    __syncthreads(); compute_first_derivs_smem(smem, &d_Axx_x, &d_Axx_y, &d_Axx_z, dX, dY, dZ, i, j, k, bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax); __syncthreads();

    FieldReal d_Axy_x, d_Axy_y, d_Axy_z;
    load_field_to_smem(smem, Axy, dims, block_i, block_j, block_k, ANTI, ANTI, SYM);
    __syncthreads(); compute_first_derivs_smem(smem, &d_Axy_x, &d_Axy_y, &d_Axy_z, dX, dY, dZ, i, j, k, bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax); __syncthreads();

    FieldReal d_Axz_x, d_Axz_y, d_Axz_z;
    load_field_to_smem(smem, Axz, dims, block_i, block_j, block_k, ANTI, SYM, ANTI);
    __syncthreads(); compute_first_derivs_smem(smem, &d_Axz_x, &d_Axz_y, &d_Axz_z, dX, dY, dZ, i, j, k, bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax); __syncthreads();

    FieldReal d_Ayy_x, d_Ayy_y, d_Ayy_z;
    load_field_to_smem(smem, Ayy, dims, block_i, block_j, block_k, SYM, SYM, SYM);
    __syncthreads(); compute_first_derivs_smem(smem, &d_Ayy_x, &d_Ayy_y, &d_Ayy_z, dX, dY, dZ, i, j, k, bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax); __syncthreads();

    FieldReal d_Ayz_x, d_Ayz_y, d_Ayz_z;
    load_field_to_smem(smem, Ayz, dims, block_i, block_j, block_k, SYM, ANTI, ANTI);
    __syncthreads(); compute_first_derivs_smem(smem, &d_Ayz_x, &d_Ayz_y, &d_Ayz_z, dX, dY, dZ, i, j, k, bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax); __syncthreads();

    FieldReal d_Azz_x, d_Azz_y, d_Azz_z;
    load_field_to_smem(smem, Azz, dims, block_i, block_j, block_k, SYM, SYM, SYM);
    __syncthreads(); compute_first_derivs_smem(smem, &d_Azz_x, &d_Azz_y, &d_Azz_z, dX, dY, dZ, i, j, k, bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax); __syncthreads();

    FieldReal chix, chiy, chiz, chixx, chixy, chixz, chiyy, chiyz, chizz;
    load_field_to_smem(smem, chi, dims, block_i, block_j, block_k, SYM, SYM, SYM);
    __syncthreads(); compute_all_derivs_smem(smem, &chix, &chiy, &chiz, &chixx, &chixy, &chixz, &chiyy, &chiyz, &chizz, dX, dY, dZ, i, j, k, bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax); __syncthreads();

    FieldReal Kx, Ky, Kz;
    load_field_to_smem(smem, trK, dims, block_i, block_j, block_k, SYM, SYM, SYM);
    __syncthreads(); compute_first_derivs_smem(smem, &Kx, &Ky, &Kz, dX, dY, dZ, i, j, k, bounds.imin, bounds.jmin, bounds.kmin, bounds.imax, bounds.jmax, bounds.kmax); 

    if (i >= ex0 || j >= ex1 || k >= ex2) return;
    int idx = flatten_bssn_index(i, j, k, ex0, ex1);

    const FieldReal ONE = FieldReal(1.0);
    const FieldReal TWO = FieldReal(2.0);
    const FieldReal F2o3 = FieldReal(2.0) / FieldReal(3.0);
    const FieldReal EIGHT = FieldReal(8.0);
    const FieldReal F16 = FieldReal(16.0);
    const FieldReal PI = FieldReal(M_PI);

    FieldReal chin1 = chi[idx] + ONE; FieldReal val_trK = trK[idx];
    FieldReal inv_chin1 = ONE / chin1;

    FieldReal l_gxx = dxx[idx] + ONE; FieldReal l_gxy = gxy[idx]; FieldReal l_gxz = gxz[idx];
    FieldReal l_gyy = dyy[idx] + ONE; FieldReal l_gyz = gyz[idx]; FieldReal l_gzz = dzz[idx] + ONE;
    FieldReal l_Axx = Axx[idx]; FieldReal l_Axy = Axy[idx]; FieldReal l_Axz = Axz[idx];
    FieldReal l_Ayy = Ayy[idx]; FieldReal l_Ayz = Ayz[idx]; FieldReal l_Azz = Azz[idx];

    FieldReal det = l_gxx * l_gyy * l_gzz + l_gxy * l_gyz * l_gxz
               + l_gxz * l_gxy * l_gyz
               - l_gxz * l_gyy * l_gxz
               - l_gxy * l_gxy * l_gzz
               - l_gxx * l_gyz * l_gyz;
    FieldReal inv_det = ONE / det;

    FieldReal gupxx = (l_gyy * l_gzz - l_gyz * l_gyz) * inv_det;
    FieldReal gupxy = (l_gyz * l_gxz - l_gxy * l_gzz) * inv_det;
    FieldReal gupxz = (l_gxy * l_gyz - l_gyy * l_gxz) * inv_det;
    FieldReal gupyy = (l_gxx * l_gzz - l_gxz * l_gxz) * inv_det;
    FieldReal gupyz = (l_gxy * l_gxz - l_gxx * l_gyz) * inv_det;
    FieldReal gupzz = (l_gxx * l_gyy - l_gxy * l_gxy) * inv_det;

    const FieldReal l_Gamxxx = Gamxxx_in[idx]; const FieldReal l_Gamxxy = Gamxxy_in[idx]; const FieldReal l_Gamxxz = Gamxxz_in[idx];
    const FieldReal l_Gamxyy = Gamxyy_in[idx]; const FieldReal l_Gamxyz = Gamxyz_in[idx]; const FieldReal l_Gamxzz = Gamxzz_in[idx];
    const FieldReal l_Gamyxx = Gamyxx_in[idx]; const FieldReal l_Gamyxy = Gamyxy_in[idx]; const FieldReal l_Gamyxz = Gamyxz_in[idx];
    const FieldReal l_Gamyyy = Gamyyy_in[idx]; const FieldReal l_Gamyyz = Gamyyz_in[idx]; const FieldReal l_Gamyzz = Gamyzz_in[idx];
    const FieldReal l_Gamzxx = Gamzxx_in[idx]; const FieldReal l_Gamzxy = Gamzxy_in[idx]; const FieldReal l_Gamzxz = Gamzxz_in[idx];
    const FieldReal l_Gamzyy = Gamzyy_in[idx]; const FieldReal l_Gamzyz = Gamzyz_in[idx]; const FieldReal l_Gamzzz = Gamzzz_in[idx];

    FieldReal Gamxa = gupxx * l_Gamxxx + gupyy * l_Gamxyy + gupzz * l_Gamxzz + TWO * (gupxy * l_Gamxxy + gupxz * l_Gamxxz + gupyz * l_Gamxyz);
    FieldReal Gamya = gupxx * l_Gamyxx + gupyy * l_Gamyyy + gupzz * l_Gamyzz + TWO * (gupxy * l_Gamyxy + gupxz * l_Gamyxz + gupyz * l_Gamyyz);
    FieldReal Gamza = gupxx * l_Gamzxx + gupyy * l_Gamzyy + gupzz * l_Gamzzz + TWO * (gupxy * l_Gamzxy + gupxz * l_Gamzxz + gupyz * l_Gamzyz);

    Gmx_Res[idx] = Gamx[idx] - Gamxa; Gmy_Res[idx] = Gamy[idx] - Gamya; Gmz_Res[idx] = Gamz[idx] - Gamza;

    FieldReal Au_d_xx = gupxx * l_Axx + gupxy * l_Axy + gupxz * l_Axz;
    FieldReal Au_d_xy = gupxx * l_Axy + gupxy * l_Ayy + gupxz * l_Ayz;
    FieldReal Au_d_xz = gupxx * l_Axz + gupxy * l_Ayz + gupxz * l_Azz;
    FieldReal Au_d_yx = gupxy * l_Axx + gupyy * l_Axy + gupyz * l_Axz;
    FieldReal Au_d_yy = gupxy * l_Axy + gupyy * l_Ayy + gupyz * l_Ayz;
    FieldReal Au_d_yz = gupxy * l_Axz + gupyy * l_Ayz + gupyz * l_Azz;
    FieldReal Au_d_zx = gupxz * l_Axx + gupyz * l_Axy + gupzz * l_Axz;
    FieldReal Au_d_zy = gupxz * l_Axy + gupyz * l_Ayy + gupzz * l_Ayz;
    FieldReal Au_d_zz = gupxz * l_Axz + gupyz * l_Ayz + gupzz * l_Azz;

    FieldReal trA2 = Au_d_xx*Au_d_xx + Au_d_yy*Au_d_yy + Au_d_zz*Au_d_zz + TWO*(Au_d_xy*Au_d_yx + Au_d_xz*Au_d_zx + Au_d_yz*Au_d_zy);

    FieldReal conf_R = gupxx * Rxx_in[idx] + gupyy * Ryy_in[idx] + gupzz * Rzz_in[idx] + TWO * (gupxy * Rxy_in[idx] + gupxz * Rxz_in[idx] + gupyz * Ryz_in[idx]);
    
    FieldReal fxx = chixx - l_Gamxxx * chix - l_Gamyxx * chiy - l_Gamzxx * chiz;
    FieldReal fyy = chiyy - l_Gamxyy * chix - l_Gamyyy * chiy - l_Gamzyy * chiz;
    FieldReal fzz = chizz - l_Gamxzz * chix - l_Gamyzz * chiy - l_Gamzzz * chiz;
    FieldReal fxy = chixy - l_Gamxxy * chix - l_Gamyxy * chiy - l_Gamzxy * chiz;
    FieldReal fxz = chixz - l_Gamxxz * chix - l_Gamyxz * chiy - l_Gamzxz * chiz;
    FieldReal fyz = chiyz - l_Gamxyz * chix - l_Gamyyz * chiy - l_Gamzyz * chiz;

    FieldReal D2_chi = gupxx * fxx + gupyy * fyy + gupzz * fzz + TWO * (gupxy * fxy + gupxz * fxz + gupyz * fyz);
    FieldReal grad_chi_sq = gupxx * chix*chix + gupyy * chiy*chiy + gupzz * chiz*chiz + TWO * (gupxy * chix*chiy + gupxz * chix*chiz + gupyz * chiy*chiz);
    
    FieldReal phys_R = chin1 * conf_R + TWO * D2_chi - FieldReal(2.5) * inv_chin1 * grad_chi_sq;
    ham_Res[idx] = phys_R + F2o3 * val_trK * val_trK - trA2 - F16 * PI * rho[idx];

    FieldReal div_A_x = gupxx * d_Axx_x + gupyy * d_Axy_y + gupzz * d_Axz_z + gupxy * (d_Axx_y + d_Axy_x) + gupxz * (d_Axx_z + d_Axz_x) + gupyz * (d_Axy_z + d_Axz_y);
    FieldReal div_A_y = gupxx * d_Axy_x + gupyy * d_Ayy_y + gupzz * d_Ayz_z + gupxy * (d_Axy_y + d_Ayy_x) + gupxz * (d_Axy_z + d_Ayz_x) + gupyz * (d_Ayy_z + d_Ayz_y);
    FieldReal div_A_z = gupxx * d_Axz_x + gupyy * d_Ayz_y + gupzz * d_Azz_z + gupxy * (d_Axz_y + d_Ayz_x) + gupxz * (d_Axz_z + d_Azz_x) + gupyz * (d_Ayz_z + d_Azz_y);

    FieldReal Gam_A_x = Gamxa * l_Axx + Gamya * l_Axy + Gamza * l_Axz;
    FieldReal Gam_A_y = Gamxa * l_Axy + Gamya * l_Ayy + Gamza * l_Ayz;
    FieldReal Gam_A_z = Gamxa * l_Axz + Gamya * l_Ayz + Gamza * l_Azz;

    FieldReal A_Gam_x = Au_d_xx * l_Gamxxx + Au_d_yx * l_Gamxxy + Au_d_zx * l_Gamxxz + Au_d_xy * l_Gamyxx + Au_d_yy * l_Gamyxy + Au_d_zy * l_Gamyxz + Au_d_xz * l_Gamzxx + Au_d_yz * l_Gamzxy + Au_d_zz * l_Gamzxz;
    FieldReal A_Gam_y = Au_d_xx * l_Gamxxy + Au_d_yx * l_Gamxyy + Au_d_zx * l_Gamxyz + Au_d_xy * l_Gamyxy + Au_d_yy * l_Gamyyy + Au_d_zy * l_Gamyyz + Au_d_xz * l_Gamzxy + Au_d_yz * l_Gamzyy + Au_d_zz * l_Gamzyz;
    FieldReal A_Gam_z = Au_d_xx * l_Gamxxz + Au_d_yx * l_Gamxyz + Au_d_zx * l_Gamxzz + Au_d_xy * l_Gamyxz + Au_d_yy * l_Gamyyz + Au_d_zy * l_Gamyzz + Au_d_xz * l_Gamzxz + Au_d_yz * l_Gamzyz + Au_d_zz * l_Gamzzz;

    FieldReal chi_A_x = FieldReal(1.5) * inv_chin1 * (Au_d_xx * chix + Au_d_yx * chiy + Au_d_zx * chiz);
    FieldReal chi_A_y = FieldReal(1.5) * inv_chin1 * (Au_d_xy * chix + Au_d_yy * chiy + Au_d_zy * chiz);
    FieldReal chi_A_z = FieldReal(1.5) * inv_chin1 * (Au_d_xz * chix + Au_d_yz * chiy + Au_d_zz * chiz);

    movx_Res[idx] = div_A_x - Gam_A_x - A_Gam_x - chi_A_x - F2o3 * Kx - EIGHT * PI * Sx[idx];
    movy_Res[idx] = div_A_y - Gam_A_y - A_Gam_y - chi_A_y - F2o3 * Ky - EIGHT * PI * Sy[idx];
    movz_Res[idx] = div_A_z - Gam_A_z - A_Gam_z - chi_A_z - F2o3 * Kz - EIGHT * PI * Sz[idx];
}

__launch_bounds__(256, 2)
__global__ void apply_z4c_algebraic_corrections_kernel(
    const int ex0, const int ex1, const int ex2,
    const int symmetry,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    const double* __restrict__ chi, const double* __restrict__ trK, const double* __restrict__ TZ,
    const double* __restrict__ dxx, const double* __restrict__ gxy, const double* __restrict__ gxz,
    const double* __restrict__ dyy, const double* __restrict__ gyz, const double* __restrict__ dzz,
    const double* __restrict__ Axx, const double* __restrict__ Axy, const double* __restrict__ Axz,
    const double* __restrict__ Ayy, const double* __restrict__ Ayz, const double* __restrict__ Azz,
    const double* __restrict__ Lap,
    const double* __restrict__ rho,
    const double* __restrict__ Sxx, const double* __restrict__ Sxy, const double* __restrict__ Sxz,
    const double* __restrict__ Syy, const double* __restrict__ Syz, const double* __restrict__ Szz,
    double* __restrict__ chi_rhs, double* __restrict__ trK_rhs,
    double* __restrict__ Axx_rhs, double* __restrict__ Axy_rhs, double* __restrict__ Axz_rhs,
    double* __restrict__ Ayy_rhs, double* __restrict__ Ayz_rhs, double* __restrict__ Azz_rhs,
    double* __restrict__ TZ_rhs,
    double* __restrict__ ham_Res,
    double* __restrict__ movx_Res, double* __restrict__ movy_Res, double* __restrict__ movz_Res
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    if (i >= ex0 || j >= ex1 || k >= ex2) return;

    int idx = flatten_bssn_index(i, j, k, ex0, ex1);
    int dims[3] = {ex0, ex1, ex2};

    double val_chi = chi[idx];
    double chin1 = val_chi + ONE;
    double inv_chin1 = ONE / chin1;
    double val_trK = trK[idx];
    double val_TZ = TZ[idx];
    double trKd = val_trK + TWO * val_TZ;
    double alpn1 = Lap[idx] + ONE;

    double l_gxx = dxx[idx] + ONE;
    double l_gxy = gxy[idx];
    double l_gxz = gxz[idx];
    double l_gyy = dyy[idx] + ONE;
    double l_gyz = gyz[idx];
    double l_gzz = dzz[idx] + ONE;

    double l_Axx = Axx[idx];
    double l_Axy = Axy[idx];
    double l_Axz = Axz[idx];
    double l_Ayy = Ayy[idx];
    double l_Ayz = Ayz[idx];
    double l_Azz = Azz[idx];

    double gupxx, gupxy, gupxz, gupyy, gupyz, gupzz;
    compute_inverse_metric(l_gxx, l_gxy, l_gxz, l_gyy, l_gyz, l_gzz,
                           &gupxx, &gupxy, &gupxz, &gupyy, &gupyz, &gupzz);

    double Au_d_xx, Au_d_xy, Au_d_xz, Au_d_yx, Au_d_yy, Au_d_yz, Au_d_zx, Au_d_zy, Au_d_zz;
    double Aup_xx, Aup_xy, Aup_xz, Aup_yy, Aup_yz, Aup_zz;
    compute_raised_extrinsic_curvature(
        l_Axx, l_Axy, l_Axz, l_Ayy, l_Ayz, l_Azz,
        gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
        &Au_d_xx, &Au_d_xy, &Au_d_xz, &Au_d_yx, &Au_d_yy, &Au_d_yz, &Au_d_zx, &Au_d_zy, &Au_d_zz,
        &Aup_xx, &Aup_xy, &Aup_xz, &Aup_yy, &Aup_yz, &Aup_zz
    );

    double trA2 = l_Axx * Aup_xx + l_Ayy * Aup_yy + l_Azz * Aup_zz
                + TWO * (l_Axy * Aup_xy + l_Axz * Aup_xz + l_Ayz * Aup_yz);
    double trace_S = gupxx * Sxx[idx] + gupyy * Syy[idx] + gupzz * Szz[idx]
                   + TWO * (gupxy * Sxy[idx] + gupxz * Sxz[idx] + gupyz * Syz[idx]);
    double S = chin1 * trace_S;

    double d2_lap = (-trK_rhs[idx]
        + alpn1 * (F1o3 * val_trK * val_trK + trA2 + FOUR_PI * (rho[idx] + S))) * inv_chin1;
    double bssn_trace_aux = F2o3 * val_trK * val_trK - trA2 - F16_PI * rho[idx] + EIGHT_PI * S;
    double f_trace_bssn = -F1o3 * (d2_lap + alpn1 * inv_chin1 * bssn_trace_aux);

    double ham_bssn = ham_Res[idx];
    double ricci_trace = (ham_bssn - F2o3 * val_trK * val_trK + trA2 + F16_PI * rho[idx]) * inv_chin1;
    double f_trace_z4c = F1o3 * (alpn1 * ricci_trace - (EIGHT_PI * alpn1 * trace_S + d2_lap));
    double trace_delta = chin1 * (f_trace_bssn - f_trace_z4c);
    double trKd_delta = alpn1 * (trKd - val_trK);

    chi_rhs[idx] += F2o3 * chin1 * alpn1 * (trKd - val_trK);

    Axx_rhs[idx] += trace_delta * l_gxx + trKd_delta * l_Axx;
    Ayy_rhs[idx] += trace_delta * l_gyy + trKd_delta * l_Ayy;
    Azz_rhs[idx] += trace_delta * l_gzz + trKd_delta * l_Azz;
    Axy_rhs[idx] += trace_delta * l_gxy + trKd_delta * l_Axy;
    Ayz_rhs[idx] += trace_delta * l_gyz + trKd_delta * l_Ayz;
    Axz_rhs[idx] += trace_delta * l_gxz + trKd_delta * l_Axz;

    trK_rhs[idx] += alpn1 * F1o3 * (trKd * trKd - val_trK * val_trK);

    double ham_z4c = ham_bssn + F2o3 * (trKd * trKd - val_trK * val_trK);
    TZ_rhs[idx] = HALF * alpn1 * ham_z4c;

    (void)dims;
    (void)movx_Res;
    (void)movy_Res;
    (void)movz_Res;
}

__launch_bounds__(256, 2)
__global__ void apply_z4c_gauge_corrections_kernel(
    const int ex0, const int ex1, const int ex2,
    const int symmetry,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    const double* __restrict__ TZ,
    const double* __restrict__ dxx, const double* __restrict__ gxy, const double* __restrict__ gxz,
    const double* __restrict__ dyy, const double* __restrict__ gyz, const double* __restrict__ dzz,
    const double* __restrict__ Gamx, const double* __restrict__ Gamy, const double* __restrict__ Gamz,
    const double* __restrict__ Lap,
    const double* __restrict__ betax, const double* __restrict__ betay, const double* __restrict__ betaz,
    double* __restrict__ Gamx_rhs, double* __restrict__ Gamy_rhs, double* __restrict__ Gamz_rhs,
    double* __restrict__ betax_rhs, double* __restrict__ betay_rhs, double* __restrict__ betaz_rhs,
    double* __restrict__ dtSfx_rhs, double* __restrict__ dtSfy_rhs, double* __restrict__ dtSfz_rhs
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    if (i >= ex0 || j >= ex1 || k >= ex2) return;

    int idx = flatten_bssn_index(i, j, k, ex0, ex1);
    int dims[3] = {ex0, ex1, ex2};

    double l_gxx = dxx[idx] + ONE;
    double l_gxy = gxy[idx];
    double l_gxz = gxz[idx];
    double l_gyy = dyy[idx] + ONE;
    double l_gyz = gyz[idx];
    double l_gzz = dzz[idx] + ONE;

    double gupxx, gupxy, gupxz, gupyy, gupyz, gupzz;
    compute_inverse_metric(l_gxx, l_gxy, l_gxz, l_gyy, l_gyz, l_gzz,
                           &gupxx, &gupxy, &gupxz, &gupyy, &gupyz, &gupzz);

    double TZx, TZy, TZz;
    d_fderivs_point(dims, TZ, &TZx, &TZy, &TZz, X, Y, Z, SYM, SYM, SYM, symmetry, 0, i, j, k);

    double alpn1 = Lap[idx] + ONE;
    double delta_x = -F2o3 * alpn1 * (gupxx * TZx + gupxy * TZy + gupxz * TZz);
    double delta_y = -F2o3 * alpn1 * (gupxy * TZx + gupyy * TZy + gupyz * TZz);
    double delta_z = -F2o3 * alpn1 * (gupxz * TZx + gupyz * TZy + gupzz * TZz);

    Gamx_rhs[idx] += delta_x;
    Gamy_rhs[idx] += delta_y;
    Gamz_rhs[idx] += delta_z;

#if (GAUGE == 0)
    dtSfx_rhs[idx] += delta_x;
    dtSfy_rhs[idx] += delta_y;
    dtSfz_rhs[idx] += delta_z;
#elif (GAUGE == 1)
    betax_rhs[idx] = Gamx[idx] - eta * betax[idx];
    betay_rhs[idx] = Gamy[idx] - eta * betay[idx];
    betaz_rhs[idx] = Gamz[idx] - eta * betaz[idx];
    dtSfx_rhs[idx] = 0.0;
    dtSfy_rhs[idx] = 0.0;
    dtSfz_rhs[idx] = 0.0;
#else
#if (ABEtype == 2)
#error Unsupported Z4C GPU gauge; expected GAUGE == 0 or GAUGE == 1.
#endif
#endif
}

__launch_bounds__(256, 2)
__global__ void apply_z4c_transport_damping_kernel(
    const int ex0, const int ex1, const int ex2,
    const int symmetry, const double eps,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    const double* __restrict__ TZ, const double* __restrict__ Lap,
    const double* __restrict__ betax, const double* __restrict__ betay, const double* __restrict__ betaz,
    const double* __restrict__ Gmx_Res, const double* __restrict__ Gmy_Res, const double* __restrict__ Gmz_Res,
    double* __restrict__ trK_rhs,
    double* __restrict__ Gamx_rhs, double* __restrict__ Gamy_rhs, double* __restrict__ Gamz_rhs,
    double* __restrict__ TZ_rhs
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    if (i >= ex0 || j >= ex1 || k >= ex2) return;

    int idx = flatten_bssn_index(i, j, k, ex0, ex1);
    int dims[3] = {ex0, ex1, ex2};

    TZ_rhs[idx] += d_lopsided_point(
        dims, TZ, TZ_rhs, betax, betay, betaz, X, Y, Z,
        symmetry, SYM, SYM, SYM, i, j, k
    );
    if (eps > 0.0) {
        TZ_rhs[idx] += d_kodis_point(
            dims, TZ, X, Y, Z, SYM, SYM, SYM, symmetry, eps, i, j, k
        );
    }

    double alpn1 = Lap[idx] + ONE;
    TZ_rhs[idx] -= alpn1 * (TWO + Z4C_KAPPA2) * Z4C_KAPPA1 * TZ[idx];
    trK_rhs[idx] += alpn1 * Z4C_KAPPA1 * (ONE - Z4C_KAPPA2) * TZ[idx];
    Gamx_rhs[idx] -= TWO * alpn1 * Z4C_KAPPA1 * Gmx_Res[idx];
    Gamy_rhs[idx] -= TWO * alpn1 * Z4C_KAPPA1 * Gmy_Res[idx];
    Gamz_rhs[idx] -= TWO * alpn1 * Z4C_KAPPA1 * Gmz_Res[idx];
}

__launch_bounds__(256, 2)
__global__ void finalize_z4c_constraint_outputs_kernel(
    const int ex0, const int ex1, const int ex2,
    const int symmetry,
    const double* __restrict__ X, const double* __restrict__ Y, const double* __restrict__ Z,
    const double* __restrict__ chi, const double* __restrict__ trK,
    const double* __restrict__ dxx, const double* __restrict__ gxy, const double* __restrict__ gxz,
    const double* __restrict__ dyy, const double* __restrict__ gyz, const double* __restrict__ dzz,
    const double* __restrict__ Axx, const double* __restrict__ Axy, const double* __restrict__ Axz,
    const double* __restrict__ Ayy, const double* __restrict__ Ayz, const double* __restrict__ Azz,
    const double* __restrict__ Gamx, const double* __restrict__ Gamy, const double* __restrict__ Gamz,
    const double* __restrict__ Sx, const double* __restrict__ Sy, const double* __restrict__ Sz,
    const double* __restrict__ Gamxxx_in, const double* __restrict__ Gamxxy_in, const double* __restrict__ Gamxxz_in,
    const double* __restrict__ Gamxyy_in, const double* __restrict__ Gamxyz_in, const double* __restrict__ Gamxzz_in,
    const double* __restrict__ Gamyxx_in, const double* __restrict__ Gamyxy_in, const double* __restrict__ Gamyxz_in,
    const double* __restrict__ Gamyyy_in, const double* __restrict__ Gamyyz_in, const double* __restrict__ Gamyzz_in,
    const double* __restrict__ Gamzxx_in, const double* __restrict__ Gamzxy_in, const double* __restrict__ Gamzxz_in,
    const double* __restrict__ Gamzyy_in, const double* __restrict__ Gamzyz_in, const double* __restrict__ Gamzzz_in,
    double* __restrict__ movx_Res, double* __restrict__ movy_Res, double* __restrict__ movz_Res,
    double* __restrict__ Gmx_Res, double* __restrict__ Gmy_Res, double* __restrict__ Gmz_Res
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    if (i >= ex0 || j >= ex1 || k >= ex2) return;

    int idx = flatten_bssn_index(i, j, k, ex0, ex1);
    int dims[3] = {ex0, ex1, ex2};

    double l_gxx = dxx[idx] + ONE;
    double l_gxy = gxy[idx];
    double l_gxz = gxz[idx];
    double l_gyy = dyy[idx] + ONE;
    double l_gyz = gyz[idx];
    double l_gzz = dzz[idx] + ONE;

    double gupxx, gupxy, gupxz, gupyy, gupyz, gupzz;
    compute_inverse_metric(l_gxx, l_gxy, l_gxz, l_gyy, l_gyz, l_gzz,
                           &gupxx, &gupxy, &gupxz, &gupyy, &gupyz, &gupzz);

    double gxxx, gxxy, gxxz;
    double gxyx, gxyy, gxyz;
    double gxzx, gxzy, gxzz;
    double gyyx, gyyy, gyyz;
    double gyzx, gyzy, gyzz;
    double gzzx, gzzy, gzzz;
    d_fderivs_point(dims, dxx, &gxxx, &gxxy, &gxxz, X, Y, Z, SYM, SYM, SYM, symmetry, 0, i, j, k);
    d_fderivs_point(dims, gxy, &gxyx, &gxyy, &gxyz, X, Y, Z, ANTI, ANTI, SYM, symmetry, 0, i, j, k);
    d_fderivs_point(dims, gxz, &gxzx, &gxzy, &gxzz, X, Y, Z, ANTI, SYM, ANTI, symmetry, 0, i, j, k);
    d_fderivs_point(dims, dyy, &gyyx, &gyyy, &gyyz, X, Y, Z, SYM, SYM, SYM, symmetry, 0, i, j, k);
    d_fderivs_point(dims, gyz, &gyzx, &gyzy, &gyzz, X, Y, Z, SYM, ANTI, ANTI, symmetry, 0, i, j, k);
    d_fderivs_point(dims, dzz, &gzzx, &gzzy, &gzzz, X, Y, Z, SYM, SYM, SYM, symmetry, 0, i, j, k);

    Gmx_Res[idx] = Gamx[idx] - (
        gupxx*(gupxx*gxxx + gupxy*gxyx + gupxz*gxzx) +
        gupxy*(gupxx*gxyx + gupxy*gyyx + gupxz*gyzx) +
        gupxz*(gupxx*gxzx + gupxy*gyzx + gupxz*gzzx) +
        gupxx*(gupxy*gxxy + gupyy*gxyy + gupyz*gxzy) +
        gupxy*(gupxy*gxyy + gupyy*gyyy + gupyz*gyzy) +
        gupxz*(gupxy*gxzy + gupyy*gyzy + gupyz*gzzy) +
        gupxx*(gupxz*gxxz + gupyz*gxyz + gupzz*gxzz) +
        gupxy*(gupxz*gxyz + gupyz*gyyz + gupzz*gyzz) +
        gupxz*(gupxz*gxzz + gupyz*gyzz + gupzz*gzzz));
    Gmy_Res[idx] = Gamy[idx] - (
        gupxx*(gupxy*gxxx + gupyy*gxyx + gupyz*gxzx) +
        gupxy*(gupxy*gxyx + gupyy*gyyx + gupyz*gyzx) +
        gupxz*(gupxy*gxzx + gupyy*gyzx + gupyz*gzzx) +
        gupxy*(gupxy*gxxy + gupyy*gxyy + gupyz*gxzy) +
        gupyy*(gupxy*gxyy + gupyy*gyyy + gupyz*gyzy) +
        gupyz*(gupxy*gxzy + gupyy*gyzy + gupyz*gzzy) +
        gupxy*(gupxz*gxxz + gupyz*gxyz + gupzz*gxzz) +
        gupyy*(gupxz*gxyz + gupyz*gyyz + gupzz*gyzz) +
        gupyz*(gupxz*gxzz + gupyz*gyzz + gupzz*gzzz));
    Gmz_Res[idx] = Gamz[idx] - (
        gupxx*(gupxz*gxxx + gupyz*gxyx + gupzz*gxzx) +
        gupxy*(gupxz*gxyx + gupyz*gyyx + gupzz*gyzx) +
        gupxz*(gupxz*gxzx + gupyz*gyzx + gupzz*gzzx) +
        gupxy*(gupxz*gxxy + gupyz*gxyy + gupzz*gxzy) +
        gupyy*(gupxz*gxyy + gupyz*gyyy + gupzz*gyzy) +
        gupyz*(gupxz*gxzy + gupyz*gyzy + gupzz*gzzy) +
        gupxz*(gupxz*gxxz + gupyz*gxyz + gupzz*gxzz) +
        gupyz*(gupxz*gxyz + gupyz*gyyz + gupzz*gyzz) +
        gupzz*(gupxz*gxzz + gupyz*gyzz + gupzz*gzzz));

    double chix, chiy, chiz;
    d_fderivs_point(dims, chi, &chix, &chiy, &chiz, X, Y, Z, SYM, SYM, SYM, symmetry, 0, i, j, k);
    double chin1 = chi[idx] + ONE;
    double inv_chin1 = ONE / chin1;
    double gx_phy = (gupxx * chix + gupxy * chiy + gupxz * chiz) * inv_chin1;
    double gy_phy = (gupxy * chix + gupyy * chiy + gupyz * chiz) * inv_chin1;
    double gz_phy = (gupxz * chix + gupyz * chiy + gupzz * chiz) * inv_chin1;

    double Gamxxx = Gamxxx_in[idx] - ((chix + chix) * inv_chin1 - l_gxx * gx_phy) * HALF;
    double Gamyxx = Gamyxx_in[idx] - (                         - l_gxx * gy_phy) * HALF;
    double Gamzxx = Gamzxx_in[idx] - (                         - l_gxx * gz_phy) * HALF;
    double Gamxyy = Gamxyy_in[idx] - (                         - l_gyy * gx_phy) * HALF;
    double Gamyyy = Gamyyy_in[idx] - ((chiy + chiy) * inv_chin1 - l_gyy * gy_phy) * HALF;
    double Gamzyy = Gamzyy_in[idx] - (                         - l_gyy * gz_phy) * HALF;
    double Gamxzz = Gamxzz_in[idx] - (                         - l_gzz * gx_phy) * HALF;
    double Gamyzz = Gamyzz_in[idx] - (                         - l_gzz * gy_phy) * HALF;
    double Gamzzz = Gamzzz_in[idx] - ((chiz + chiz) * inv_chin1 - l_gzz * gz_phy) * HALF;
    double Gamxxy = Gamxxy_in[idx] - ( chiy * inv_chin1 - l_gxy * gx_phy) * HALF;
    double Gamyxy = Gamyxy_in[idx] - ( chix * inv_chin1 - l_gxy * gy_phy) * HALF;
    double Gamzxy = Gamzxy_in[idx] - (                  - l_gxy * gz_phy) * HALF;
    double Gamxxz = Gamxxz_in[idx] - ( chiz * inv_chin1 - l_gxz * gx_phy) * HALF;
    double Gamyxz = Gamyxz_in[idx] - (                  - l_gxz * gy_phy) * HALF;
    double Gamzxz = Gamzxz_in[idx] - ( chix * inv_chin1 - l_gxz * gz_phy) * HALF;
    double Gamxyz = Gamxyz_in[idx] - (                  - l_gyz * gx_phy) * HALF;
    double Gamyyz = Gamyyz_in[idx] - ( chiz * inv_chin1 - l_gyz * gy_phy) * HALF;
    double Gamzyz = Gamzyz_in[idx] - ( chiy * inv_chin1 - l_gyz * gz_phy) * HALF;

    double Axxx, Axxy, Axxz;
    double Axyx, Axyy, Axyz;
    double Axzx, Axzy, Axzz;
    double Ayyx, Ayyy, Ayyz;
    double Ayzx, Ayzy, Ayzz;
    double Azzx, Azzy, Azzz;
    d_fderivs_point(dims, Axx, &Axxx, &Axxy, &Axxz, X, Y, Z, SYM, SYM, SYM, symmetry, 0, i, j, k);
    d_fderivs_point(dims, Axy, &Axyx, &Axyy, &Axyz, X, Y, Z, ANTI, ANTI, SYM, symmetry, 0, i, j, k);
    d_fderivs_point(dims, Axz, &Axzx, &Axzy, &Axzz, X, Y, Z, ANTI, SYM, ANTI, symmetry, 0, i, j, k);
    d_fderivs_point(dims, Ayy, &Ayyx, &Ayyy, &Ayyz, X, Y, Z, SYM, SYM, SYM, symmetry, 0, i, j, k);
    d_fderivs_point(dims, Ayz, &Ayzx, &Ayzy, &Ayzz, X, Y, Z, SYM, ANTI, ANTI, symmetry, 0, i, j, k);
    d_fderivs_point(dims, Azz, &Azzx, &Azzy, &Azzz, X, Y, Z, SYM, SYM, SYM, symmetry, 0, i, j, k);

    double l_Axx = Axx[idx], l_Axy = Axy[idx], l_Axz = Axz[idx];
    double l_Ayy = Ayy[idx], l_Ayz = Ayz[idx], l_Azz = Azz[idx];

    double Dxxx = Axxx - (Gamxxx*l_Axx + Gamyxx*l_Axy + Gamzxx*l_Axz
                        + Gamxxx*l_Axx + Gamyxx*l_Axy + Gamzxx*l_Axz) - chix*l_Axx*inv_chin1;
    double Dxyx = Axyx - (Gamxxy*l_Axx + Gamyxy*l_Axy + Gamzxy*l_Axz
                        + Gamxxx*l_Axy + Gamyxx*l_Ayy + Gamzxx*l_Ayz) - chix*l_Axy*inv_chin1;
    double Dzxx = Axzx - (Gamxxz*l_Axx + Gamyxz*l_Axy + Gamzxz*l_Axz
                        + Gamxxx*l_Axz + Gamyxx*l_Ayz + Gamzxx*l_Azz) - chix*l_Axz*inv_chin1;
    double Dyyx = Ayyx - (Gamxxy*l_Axy + Gamyxy*l_Ayy + Gamzxy*l_Ayz
                        + Gamxxy*l_Axy + Gamyxy*l_Ayy + Gamzxy*l_Ayz) - chix*l_Ayy*inv_chin1;
    double Dyzx = Ayzx - (Gamxxz*l_Axy + Gamyxz*l_Ayy + Gamzxz*l_Ayz
                        + Gamxxy*l_Axz + Gamyxy*l_Ayz + Gamzxy*l_Azz) - chix*l_Ayz*inv_chin1;
    double Dzzx = Azzx - (Gamxxz*l_Axz + Gamyxz*l_Ayz + Gamzxz*l_Azz
                        + Gamxxz*l_Axz + Gamyxz*l_Ayz + Gamzxz*l_Azz) - chix*l_Azz*inv_chin1;

    double Dxxy = Axxy - (Gamxxy*l_Axx + Gamyxy*l_Axy + Gamzxy*l_Axz
                        + Gamxxy*l_Axx + Gamyxy*l_Axy + Gamzxy*l_Axz) - chiy*l_Axx*inv_chin1;
    double Dxyy = Axyy - (Gamxyy*l_Axx + Gamyyy*l_Axy + Gamzyy*l_Axz
                        + Gamxxy*l_Axy + Gamyxy*l_Ayy + Gamzxy*l_Ayz) - chiy*l_Axy*inv_chin1;
    double Dzxy = Axzy - (Gamxyz*l_Axx + Gamyyz*l_Axy + Gamzyz*l_Axz
                        + Gamxxy*l_Axz + Gamyxy*l_Ayz + Gamzxy*l_Azz) - chiy*l_Axz*inv_chin1;
    double Dyyy = Ayyy - (Gamxyy*l_Axy + Gamyyy*l_Ayy + Gamzyy*l_Ayz
                        + Gamxyy*l_Axy + Gamyyy*l_Ayy + Gamzyy*l_Ayz) - chiy*l_Ayy*inv_chin1;
    double Dyzy = Ayzy - (Gamxyz*l_Axy + Gamyyz*l_Ayy + Gamzyz*l_Ayz
                        + Gamxyy*l_Axz + Gamyyy*l_Ayz + Gamzyy*l_Azz) - chiy*l_Ayz*inv_chin1;
    double Dzzy = Azzy - (Gamxyz*l_Axz + Gamyyz*l_Ayz + Gamzyz*l_Azz
                        + Gamxyz*l_Axz + Gamyyz*l_Ayz + Gamzyz*l_Azz) - chiy*l_Azz*inv_chin1;

    double Dxxz = Axxz - (Gamxxz*l_Axx + Gamyxz*l_Axy + Gamzxz*l_Axz
                        + Gamxxz*l_Axx + Gamyxz*l_Axy + Gamzxz*l_Axz) - chiz*l_Axx*inv_chin1;
    double Dxyz = Axyz - (Gamxyz*l_Axx + Gamyyz*l_Axy + Gamzyz*l_Axz
                        + Gamxxz*l_Axy + Gamyxz*l_Ayy + Gamzxz*l_Ayz) - chiz*l_Axy*inv_chin1;
    double Dzxz = Axzz - (Gamxzz*l_Axx + Gamyzz*l_Axy + Gamzzz*l_Axz
                        + Gamxxz*l_Axz + Gamyxz*l_Ayz + Gamzxz*l_Azz) - chiz*l_Axz*inv_chin1;
    double Dyyz = Ayyz - (Gamxyz*l_Axy + Gamyyz*l_Ayy + Gamzyz*l_Ayz
                        + Gamxyz*l_Axy + Gamyyz*l_Ayy + Gamzyz*l_Ayz) - chiz*l_Ayy*inv_chin1;
    double Dyzz = Ayzz - (Gamxzz*l_Axy + Gamyzz*l_Ayy + Gamzzz*l_Ayz
                        + Gamxyz*l_Axz + Gamyyz*l_Ayz + Gamzyz*l_Azz) - chiz*l_Ayz*inv_chin1;
    double Dzzz = Azzz - (Gamxzz*l_Axz + Gamyzz*l_Ayz + Gamzzz*l_Azz
                        + Gamxzz*l_Axz + Gamyzz*l_Ayz + Gamzzz*l_Azz) - chiz*l_Azz*inv_chin1;

    double Kx, Ky, Kz;
    d_fderivs_point(dims, trK, &Kx, &Ky, &Kz, X, Y, Z, SYM, SYM, SYM, symmetry, 0, i, j, k);

    movx_Res[idx] = gupxx*Dxxx + gupyy*Dxyy + gupzz*Dzxz
                  + gupxy*Dxyx + gupxz*Dzxx + gupyz*Dzxy
                  + gupxy*Dxxy + gupxz*Dxxz + gupyz*Dxyz
                  - F2o3*Kx - EIGHT_PI*Sx[idx];
    movy_Res[idx] = gupxx*Dxyx + gupyy*Dyyy + gupzz*Dyzz
                  + gupxy*Dyyx + gupxz*Dyzx + gupyz*Dyzy
                  + gupxy*Dxyy + gupxz*Dxyz + gupyz*Dyyz
                  - F2o3*Ky - EIGHT_PI*Sy[idx];
    movz_Res[idx] = gupxx*Dzxx + gupyy*Dyzy + gupzz*Dzzz
                  + gupxy*Dyzx + gupxz*Dzzx + gupyz*Dzzy
                  + gupxy*Dzxy + gupxz*Dzxz + gupyz*Dyzz
                  - F2o3*Kz - EIGHT_PI*Sz[idx];
}

__launch_bounds__(256, 2)
__global__ void compute_constraint_fr_kernel(
    const int ex0, const int ex1, const int ex2,
    const double* __restrict__ chi,
    const double* __restrict__ trK,
    const double* __restrict__ rho,
    const double* __restrict__ Sphi,
    const double* __restrict__ dxx,
    const double* __restrict__ gxy,
    const double* __restrict__ gxz,
    const double* __restrict__ dyy,
    const double* __restrict__ gyz,
    const double* __restrict__ dzz,
    const double* __restrict__ Axx,
    const double* __restrict__ Axy,
    const double* __restrict__ Axz,
    const double* __restrict__ Ayy,
    const double* __restrict__ Ayz,
    const double* __restrict__ Azz,
    const double* __restrict__ Rxx,
    const double* __restrict__ Rxy,
    const double* __restrict__ Rxz,
    const double* __restrict__ Ryy,
    const double* __restrict__ Ryz,
    const double* __restrict__ Rzz,
    const double* __restrict__ Sxx,
    const double* __restrict__ Sxy,
    const double* __restrict__ Sxz,
    const double* __restrict__ Syy,
    const double* __restrict__ Syz,
    const double* __restrict__ Szz,
    double* __restrict__ Cons_fR,
    const double a2
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    if (i >= ex0 || j >= ex1 || k >= ex2) return;

    const int idx = flatten_bssn_index(i, j, k, ex0, ex1);

    const double l_gxx = dxx[idx] + 1.0;
    const double l_gxy = gxy[idx];
    const double l_gxz = gxz[idx];
    const double l_gyy = dyy[idx] + 1.0;
    const double l_gyz = gyz[idx];
    const double l_gzz = dzz[idx] + 1.0;
    const double chipn1 = chi[idx] + 1.0;

    double gupxx, gupxy, gupxz, gupyy, gupyz, gupzz;
    compute_inverse_metric(l_gxx, l_gxy, l_gxz, l_gyy, l_gyz, l_gzz,
                           &gupxx, &gupxy, &gupxz, &gupyy, &gupyz, &gupzz);

    const double l_Axx = Axx[idx];
    const double l_Axy = Axy[idx];
    const double l_Axz = Axz[idx];
    const double l_Ayy = Ayy[idx];
    const double l_Ayz = Ayz[idx];
    const double l_Azz = Azz[idx];

    const double trA2 =
        (gupxx * (
             gupxx * l_Axx * l_Axx + gupyy * l_Axy * l_Axy + gupzz * l_Axz * l_Axz +
             2.0 * (gupxy * l_Axx * l_Axy + gupxz * l_Axx * l_Axz + gupyz * l_Axy * l_Axz)) +
         gupyy * (
             gupxx * l_Axy * l_Axy + gupyy * l_Ayy * l_Ayy + gupzz * l_Ayz * l_Ayz +
             2.0 * (gupxy * l_Axy * l_Ayy + gupxz * l_Axy * l_Ayz + gupyz * l_Ayy * l_Ayz)) +
         gupzz * (
             gupxx * l_Axz * l_Axz + gupyy * l_Ayz * l_Ayz + gupzz * l_Azz * l_Azz +
             2.0 * (gupxy * l_Axz * l_Ayz + gupxz * l_Axz * l_Azz + gupyz * l_Ayz * l_Azz)) +
         2.0 * (
             gupxy * (
                 gupxx * l_Axx * l_Axy + gupyy * l_Axy * l_Ayy + gupzz * l_Axz * l_Ayz +
                 gupxy * (l_Axx * l_Ayy + l_Axy * l_Axy) +
                 gupxz * (l_Axx * l_Ayz + l_Axz * l_Axy) +
                 gupyz * (l_Axy * l_Ayz + l_Axz * l_Ayy)) +
             gupxz * (
                 gupxx * l_Axx * l_Axz + gupyy * l_Axy * l_Ayz + gupzz * l_Axz * l_Azz +
                 gupxy * (l_Axx * l_Ayz + l_Axy * l_Axz) +
                 gupxz * (l_Axx * l_Azz + l_Axz * l_Axz) +
                 gupyz * (l_Axy * l_Azz + l_Axz * l_Ayz)) +
             gupyz * (
                 gupxx * l_Axy * l_Axz + gupyy * l_Ayy * l_Ayz + gupzz * l_Ayz * l_Azz +
                 gupxy * (l_Axy * l_Ayz + l_Ayy * l_Axz) +
                 gupxz * (l_Axy * l_Azz + l_Ayz * l_Axz) +
                 gupyz * (l_Ayy * l_Azz + l_Ayz * l_Ayz))));

    const double ricci_trace =
        gupxx * Rxx[idx] + gupyy * Ryy[idx] + gupzz * Rzz[idx] +
        2.0 * (gupxy * Rxy[idx] + gupxz * Rxz[idx] + gupyz * Ryz[idx]);
    const double stress_trace =
        gupxx * Sxx[idx] + gupyy * Syy[idx] + gupzz * Szz[idx] +
        2.0 * (gupxy * Sxy[idx] + gupxz * Sxz[idx] + gupyz * Syz[idx]);

    const double RR =
        trA2 - (2.0 / 3.0) * trK[idx] * trK[idx]
        - ricci_trace * chipn1
        - 8.0 * M_PI * (3.0 * rho[idx] - stress_trace * chipn1);

#if (EScalar_CC == 1 || EScalar_CC == 4 || EScalar_CC == 5)
    const double fprim = 1.0;
#elif (EScalar_CC == 2 || EScalar_CC == 3)
    const double fprim = 1.0 + 2.0 * a2 * RR;
#else
#error Unsupported EScalar_CC in compute_constraint_fr_kernel.
#endif
    Cons_fR[idx] = Sphi[idx] - sqrt(3.0 / M_PI) / 4.0 * log(fprim);
}

void gpu_compute_rhs_bssn_launch(
    cudaStream_t &stream, int* ex, double T, double* X, double* Y, double* Z,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Gamx, double* Gamy, double* Gamz, double* Lap,
    double* betax, double* betay, double* betaz,
    double* dtSfx, double* dtSfy, double* dtSfz,
    double* chi_rhs, double* trK_rhs,
    double* gxx_rhs, double* gxy_rhs, double* gxz_rhs, double* gyy_rhs, double* gyz_rhs, double* gzz_rhs,
    double* Axx_rhs, double* Axy_rhs, double* Axz_rhs, double* Ayy_rhs, double* Ayz_rhs, double* Azz_rhs,
    double* Gamx_rhs, double* Gamy_rhs, double* Gamz_rhs, double* Lap_rhs,
    double* betax_rhs, double* betay_rhs, double* betaz_rhs,
    double* dtSfx_rhs, double* dtSfy_rhs, double* dtSfz_rhs,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Sxx, double* Sxy, double* Sxz, double* Syy, double* Syz, double* Szz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz, double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz, double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz, double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* ham_Res, double* movx_Res, double* movy_Res, double* movz_Res,
    double* Gmx_Res, double* Gmy_Res, double* Gmz_Res,
    int symmetry, int lev, double eps, int co, bool use_escalar_trace_free
) {
    dim3 block = bssn_block_shape();
    dim3 grid = bssn_grid_shape(ex, block);

#if (GAUGE == 6 || GAUGE == 7)
    cudaMemcpyToSymbolAsync(d_bssn_gauge_bh_count, &BHN, sizeof(int), 0, cudaMemcpyHostToDevice, stream);
    cudaMemcpyToSymbolAsync(d_bssn_gauge_mass, Mass, 3 * sizeof(double), 0, cudaMemcpyHostToDevice, stream);
    cudaMemcpyToSymbolAsync(d_bssn_gauge_pbh, PBH, 9 * sizeof(double), 0, cudaMemcpyHostToDevice, stream);
#endif

    trace_metric_det_summary(stream, "entry_metric", ex, X, Y, Z, lev, co,
                             dxx, gxy, gxz, dyy, gyz, dzz);

    assemble_conformal_geometry_kernel<double><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], symmetry, lev, X, Y, Z,
        dxx, gxy, gxz, dyy, gyz, dzz, Gamx, Gamy, Gamz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz
    );
    const RhsTraceField geometry_fields[] = {
        {"Gamxxx", Gamxxx}, {"Gamxxy", Gamxxy}, {"Gamxxz", Gamxxz},
        {"Gamxyy", Gamxyy}, {"Gamxyz", Gamxyz}, {"Gamxzz", Gamxzz},
        {"Gamyxx", Gamyxx}, {"Gamyxy", Gamyxy}, {"Gamyxz", Gamyxz},
        {"Gamyyy", Gamyyy}, {"Gamyyz", Gamyyz}, {"Gamyzz", Gamyzz},
        {"Gamzxx", Gamzxx}, {"Gamzxy", Gamzxy}, {"Gamzxz", Gamzxz},
        {"Gamzyy", Gamzyy}, {"Gamzyz", Gamzyz}, {"Gamzzz", Gamzzz},
        {"Rxx", Rxx}, {"Rxy", Rxy}, {"Rxz", Rxz},
        {"Ryy", Ryy}, {"Ryz", Ryz}, {"Rzz", Rzz}
    };
    trace_rhs_nonfinite_summary(stream, "after_assemble", ex, X, Y, Z, lev, co,
                                geometry_fields, sizeof(geometry_fields) / sizeof(geometry_fields[0]));

    evolve_bssn_metric_fields_kernel<double><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], symmetry, X, Y, Z,
        chi, trK, dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz,
        Lap, betax, betay, betaz,
        chi_rhs, gxx_rhs, gxy_rhs, gxz_rhs, gyy_rhs, gyz_rhs, gzz_rhs
    );
    const RhsTraceField metric_rhs_fields[] = {
        {"phi_rhs", chi_rhs}, {"gxx_rhs", gxx_rhs}, {"gxy_rhs", gxy_rhs},
        {"gxz_rhs", gxz_rhs}, {"gyy_rhs", gyy_rhs}, {"gyz_rhs", gyz_rhs},
        {"gzz_rhs", gzz_rhs}
    };
    trace_rhs_nonfinite_summary(stream, "after_metric", ex, X, Y, Z, lev, co,
                                metric_rhs_fields, sizeof(metric_rhs_fields) / sizeof(metric_rhs_fields[0]));

    evolve_bssn_curvature_fields_kernel<double><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], symmetry, X, Y, Z,
        chi, trK, dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz,
        Lap, betax, betay, betaz,
        rho, Sxx, Sxy, Sxz, Syy, Syz, Szz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
        trK_rhs,
        Axx_rhs, Axy_rhs, Axz_rhs, Ayy_rhs, Ayz_rhs, Azz_rhs,
        use_escalar_trace_free
    );
    const RhsTraceField curvature_rhs_fields[] = {
        {"trK_rhs", trK_rhs}, {"Axx_rhs", Axx_rhs}, {"Axy_rhs", Axy_rhs},
        {"Axz_rhs", Axz_rhs}, {"Ayy_rhs", Ayy_rhs}, {"Ayz_rhs", Ayz_rhs},
        {"Azz_rhs", Azz_rhs}
    };
    trace_rhs_nonfinite_summary(stream, "after_curvature", ex, X, Y, Z, lev, co,
                                curvature_rhs_fields, sizeof(curvature_rhs_fields) / sizeof(curvature_rhs_fields[0]));

    evolve_bssn_gauge_fields_kernel<double><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], symmetry, X, Y, Z,
        chi, trK, dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz,
        Gamx, Gamy, Gamz, Lap, betax, betay, betaz, dtSfx, dtSfy, dtSfz,
        Sx, Sy, Sz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Gamx_rhs, Gamy_rhs, Gamz_rhs, Lap_rhs,
        betax_rhs, betay_rhs, betaz_rhs,
        dtSfx_rhs, dtSfy_rhs, dtSfz_rhs
    );
    const RhsTraceField gauge_rhs_fields[] = {
        {"Gamx_rhs", Gamx_rhs}, {"Gamy_rhs", Gamy_rhs}, {"Gamz_rhs", Gamz_rhs},
        {"Lap_rhs", Lap_rhs}, {"Sfx_rhs", betax_rhs}, {"Sfy_rhs", betay_rhs},
        {"Sfz_rhs", betaz_rhs}, {"dtSfx_rhs", dtSfx_rhs},
        {"dtSfy_rhs", dtSfy_rhs}, {"dtSfz_rhs", dtSfz_rhs}
    };
    trace_rhs_nonfinite_summary(stream, "after_gauge", ex, X, Y, Z, lev, co,
                                gauge_rhs_fields, sizeof(gauge_rhs_fields) / sizeof(gauge_rhs_fields[0]));

    apply_transport_and_dissipation_kernel<double><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], symmetry, lev, eps, X, Y, Z,
        betax, betay, betaz, dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz, chi, trK, Gamx, Gamy, Gamz, Lap, dtSfx, dtSfy, dtSfz,
        gxx_rhs, gxy_rhs, gxz_rhs, gyy_rhs, gyz_rhs, gzz_rhs,
        Axx_rhs, Axy_rhs, Axz_rhs, Ayy_rhs, Ayz_rhs, Azz_rhs,
        chi_rhs, trK_rhs, Gamx_rhs, Gamy_rhs, Gamz_rhs, Lap_rhs,
        betax_rhs, betay_rhs, betaz_rhs, dtSfx_rhs, dtSfy_rhs, dtSfz_rhs
    );
    const RhsTraceField all_rhs_fields[] = {
        {"gxx_rhs", gxx_rhs}, {"gxy_rhs", gxy_rhs}, {"gxz_rhs", gxz_rhs},
        {"gyy_rhs", gyy_rhs}, {"gyz_rhs", gyz_rhs}, {"gzz_rhs", gzz_rhs},
        {"Axx_rhs", Axx_rhs}, {"Axy_rhs", Axy_rhs}, {"Axz_rhs", Axz_rhs},
        {"Ayy_rhs", Ayy_rhs}, {"Ayz_rhs", Ayz_rhs}, {"Azz_rhs", Azz_rhs},
        {"phi_rhs", chi_rhs}, {"trK_rhs", trK_rhs},
        {"Gamx_rhs", Gamx_rhs}, {"Gamy_rhs", Gamy_rhs}, {"Gamz_rhs", Gamz_rhs},
        {"Lap_rhs", Lap_rhs}, {"Sfx_rhs", betax_rhs}, {"Sfy_rhs", betay_rhs},
        {"Sfz_rhs", betaz_rhs}, {"dtSfx_rhs", dtSfx_rhs},
        {"dtSfy_rhs", dtSfy_rhs}, {"dtSfz_rhs", dtSfz_rhs}
    };
    trace_rhs_nonfinite_summary(stream, "after_transport", ex, X, Y, Z, lev, co,
                                all_rhs_fields, sizeof(all_rhs_fields) / sizeof(all_rhs_fields[0]));

    if (co == 0) {
        evaluate_bssn_constraints_kernel<double><<<grid, block, 0, stream>>>(
            ex[0], ex[1], ex[2], symmetry, lev, X, Y, Z,
            chi, trK, dxx, gxy, gxz, dyy, gyz, dzz,
            Axx, Axy, Axz, Ayy, Ayz, Azz, Gamx, Gamy, Gamz, rho, Sx, Sy, Sz,
            Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
            Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
            Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
            Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
            ham_Res, movx_Res, movy_Res, movz_Res,
            Gmx_Res, Gmy_Res, Gmz_Res
        );
    }
}

void gpu_compute_rhs_bssn_fp32_launch(
    cudaStream_t &stream, int* ex, double T, double* X, double* Y, double* Z,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Gamx, double* Gamy, double* Gamz, double* Lap,
    double* betax, double* betay, double* betaz,
    double* dtSfx, double* dtSfy, double* dtSfz,
    double* chi_rhs, double* trK_rhs,
    double* gxx_rhs, double* gxy_rhs, double* gxz_rhs, double* gyy_rhs, double* gyz_rhs, double* gzz_rhs,
    double* Axx_rhs, double* Axy_rhs, double* Axz_rhs, double* Ayy_rhs, double* Ayz_rhs, double* Azz_rhs,
    double* Gamx_rhs, double* Gamy_rhs, double* Gamz_rhs, double* Lap_rhs,
    double* betax_rhs, double* betay_rhs, double* betaz_rhs,
    double* dtSfx_rhs, double* dtSfy_rhs, double* dtSfz_rhs,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Sxx, double* Sxy, double* Sxz, double* Syy, double* Syz, double* Szz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz, double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz, double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz, double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* ham_Res, double* movx_Res, double* movy_Res, double* movz_Res,
    double* Gmx_Res, double* Gmy_Res, double* Gmz_Res,
    int symmetry, int lev, double eps, int co
) {
    (void)T;
    (void)rho; (void)Sx; (void)Sy; (void)Sz;
    (void)Sxx; (void)Sxy; (void)Sxz; (void)Syy; (void)Syz; (void)Szz;

    report_bssn_fp32_rhs_once();

    dim3 block = bssn_block_shape();
    dim3 grid = bssn_grid_shape(ex, block);
    const int point_count = ex[0] * ex[1] * ex[2];

    constexpr int kStateFieldCount = 24;
    constexpr int kAuxFieldCount = 24;
    constexpr int kRhsFieldCount = 24;
    constexpr int kConstraintFieldCount = 7;
    constexpr int kZeroFieldCount = 1;
    constexpr int kTotalFloatBuffers =
        kStateFieldCount + kAuxFieldCount + kRhsFieldCount + kConstraintFieldCount + kZeroFieldCount;

    BssnRhsFp32Scratch& scratch = bssn_rhs_fp32_scratch();
    scratch.reserve(point_count, kTotalFloatBuffers, kStateFieldCount);

    int slot = 0;
    auto take_buffer = [&]() -> float* { return scratch.buffers[slot++]; };

    std::vector<float*> state_fields;
    state_fields.reserve(kStateFieldCount);
    for (int field = 0; field < kStateFieldCount; ++field) state_fields.push_back(take_buffer());

    std::vector<float*> aux_fields;
    aux_fields.reserve(kAuxFieldCount);
    for (int field = 0; field < kAuxFieldCount; ++field) aux_fields.push_back(take_buffer());

    std::vector<float*> rhs_fields;
    rhs_fields.reserve(kRhsFieldCount);
    for (int field = 0; field < kRhsFieldCount; ++field) rhs_fields.push_back(take_buffer());

    std::vector<float*> constraint_fields;
    constraint_fields.reserve(kConstraintFieldCount);
    for (int field = 0; field < kConstraintFieldCount; ++field) constraint_fields.push_back(take_buffer());

    float* zero_matter = take_buffer();

    std::vector<double*> state_sources = {
        chi, trK,
        dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz,
        Gamx, Gamy, Gamz, Lap,
        betax, betay, betaz,
        dtSfx, dtSfy, dtSfz
    };
    bssn_convert_double_to_float_fields(stream, state_sources, state_fields, point_count, scratch);
    std::vector<float*> zero_fields = {zero_matter};
    bssn_zero_float_fields(stream, zero_fields, point_count, scratch);

    float* f_chi = state_fields[0];
    float* f_trK = state_fields[1];
    float* f_dxx = state_fields[2];
    float* f_gxy = state_fields[3];
    float* f_gxz = state_fields[4];
    float* f_dyy = state_fields[5];
    float* f_gyz = state_fields[6];
    float* f_dzz = state_fields[7];
    float* f_Axx = state_fields[8];
    float* f_Axy = state_fields[9];
    float* f_Axz = state_fields[10];
    float* f_Ayy = state_fields[11];
    float* f_Ayz = state_fields[12];
    float* f_Azz = state_fields[13];
    float* f_Gamx = state_fields[14];
    float* f_Gamy = state_fields[15];
    float* f_Gamz = state_fields[16];
    float* f_Lap = state_fields[17];
    float* f_betax = state_fields[18];
    float* f_betay = state_fields[19];
    float* f_betaz = state_fields[20];
    float* f_dtSfx = state_fields[21];
    float* f_dtSfy = state_fields[22];
    float* f_dtSfz = state_fields[23];

    float* f_Gamxxx = aux_fields[0];
    float* f_Gamxxy = aux_fields[1];
    float* f_Gamxxz = aux_fields[2];
    float* f_Gamxyy = aux_fields[3];
    float* f_Gamxyz = aux_fields[4];
    float* f_Gamxzz = aux_fields[5];
    float* f_Gamyxx = aux_fields[6];
    float* f_Gamyxy = aux_fields[7];
    float* f_Gamyxz = aux_fields[8];
    float* f_Gamyyy = aux_fields[9];
    float* f_Gamyyz = aux_fields[10];
    float* f_Gamyzz = aux_fields[11];
    float* f_Gamzxx = aux_fields[12];
    float* f_Gamzxy = aux_fields[13];
    float* f_Gamzxz = aux_fields[14];
    float* f_Gamzyy = aux_fields[15];
    float* f_Gamzyz = aux_fields[16];
    float* f_Gamzzz = aux_fields[17];
    float* f_Rxx = aux_fields[18];
    float* f_Rxy = aux_fields[19];
    float* f_Rxz = aux_fields[20];
    float* f_Ryy = aux_fields[21];
    float* f_Ryz = aux_fields[22];
    float* f_Rzz = aux_fields[23];

    float* f_chi_rhs = rhs_fields[0];
    float* f_trK_rhs = rhs_fields[1];
    float* f_gxx_rhs = rhs_fields[2];
    float* f_gxy_rhs = rhs_fields[3];
    float* f_gxz_rhs = rhs_fields[4];
    float* f_gyy_rhs = rhs_fields[5];
    float* f_gyz_rhs = rhs_fields[6];
    float* f_gzz_rhs = rhs_fields[7];
    float* f_Axx_rhs = rhs_fields[8];
    float* f_Axy_rhs = rhs_fields[9];
    float* f_Axz_rhs = rhs_fields[10];
    float* f_Ayy_rhs = rhs_fields[11];
    float* f_Ayz_rhs = rhs_fields[12];
    float* f_Azz_rhs = rhs_fields[13];
    float* f_Gamx_rhs = rhs_fields[14];
    float* f_Gamy_rhs = rhs_fields[15];
    float* f_Gamz_rhs = rhs_fields[16];
    float* f_Lap_rhs = rhs_fields[17];
    float* f_betax_rhs = rhs_fields[18];
    float* f_betay_rhs = rhs_fields[19];
    float* f_betaz_rhs = rhs_fields[20];
    float* f_dtSfx_rhs = rhs_fields[21];
    float* f_dtSfy_rhs = rhs_fields[22];
    float* f_dtSfz_rhs = rhs_fields[23];

#if (GAUGE == 6 || GAUGE == 7)
    cudaMemcpyToSymbolAsync(d_bssn_gauge_bh_count, &BHN, sizeof(int), 0, cudaMemcpyHostToDevice, stream);
    cudaMemcpyToSymbolAsync(d_bssn_gauge_mass, Mass, 3 * sizeof(double), 0, cudaMemcpyHostToDevice, stream);
    cudaMemcpyToSymbolAsync(d_bssn_gauge_pbh, PBH, 9 * sizeof(double), 0, cudaMemcpyHostToDevice, stream);
#endif

    assemble_conformal_geometry_kernel<float><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], symmetry, lev, X, Y, Z,
        f_dxx, f_gxy, f_gxz, f_dyy, f_gyz, f_dzz, f_Gamx, f_Gamy, f_Gamz,
        f_Gamxxx, f_Gamxxy, f_Gamxxz, f_Gamxyy, f_Gamxyz, f_Gamxzz,
        f_Gamyxx, f_Gamyxy, f_Gamyxz, f_Gamyyy, f_Gamyyz, f_Gamyzz,
        f_Gamzxx, f_Gamzxy, f_Gamzxz, f_Gamzyy, f_Gamzyz, f_Gamzzz,
        f_Rxx, f_Rxy, f_Rxz, f_Ryy, f_Ryz, f_Rzz
    );

    evolve_bssn_metric_fields_kernel<float><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], symmetry, X, Y, Z,
        f_chi, f_trK, f_dxx, f_gxy, f_gxz, f_dyy, f_gyz, f_dzz,
        f_Axx, f_Axy, f_Axz, f_Ayy, f_Ayz, f_Azz,
        f_Lap, f_betax, f_betay, f_betaz,
        f_chi_rhs, f_gxx_rhs, f_gxy_rhs, f_gxz_rhs, f_gyy_rhs, f_gyz_rhs, f_gzz_rhs
    );

    evolve_bssn_curvature_fields_kernel<float><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], symmetry, X, Y, Z,
        f_chi, f_trK, f_dxx, f_gxy, f_gxz, f_dyy, f_gyz, f_dzz,
        f_Axx, f_Axy, f_Axz, f_Ayy, f_Ayz, f_Azz,
        f_Lap, f_betax, f_betay, f_betaz,
        zero_matter, zero_matter, zero_matter, zero_matter, zero_matter, zero_matter, zero_matter,
        f_Gamxxx, f_Gamxxy, f_Gamxxz, f_Gamxyy, f_Gamxyz, f_Gamxzz,
        f_Gamyxx, f_Gamyxy, f_Gamyxz, f_Gamyyy, f_Gamyyz, f_Gamyzz,
        f_Gamzxx, f_Gamzxy, f_Gamzxz, f_Gamzyy, f_Gamzyz, f_Gamzzz,
        f_Rxx, f_Rxy, f_Rxz, f_Ryy, f_Ryz, f_Rzz,
        f_trK_rhs,
        f_Axx_rhs, f_Axy_rhs, f_Axz_rhs, f_Ayy_rhs, f_Ayz_rhs, f_Azz_rhs,
        false
    );

    evolve_bssn_gauge_fields_kernel<float><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], symmetry, X, Y, Z,
        f_chi, f_trK, f_dxx, f_gxy, f_gxz, f_dyy, f_gyz, f_dzz,
        f_Axx, f_Axy, f_Axz, f_Ayy, f_Ayz, f_Azz,
        f_Gamx, f_Gamy, f_Gamz, f_Lap, f_betax, f_betay, f_betaz, f_dtSfx, f_dtSfy, f_dtSfz,
        zero_matter, zero_matter, zero_matter,
        f_Gamxxx, f_Gamxxy, f_Gamxxz, f_Gamxyy, f_Gamxyz, f_Gamxzz,
        f_Gamyxx, f_Gamyxy, f_Gamyxz, f_Gamyyy, f_Gamyyz, f_Gamyzz,
        f_Gamzxx, f_Gamzxy, f_Gamzxz, f_Gamzyy, f_Gamzyz, f_Gamzzz,
        f_Gamx_rhs, f_Gamy_rhs, f_Gamz_rhs, f_Lap_rhs,
        f_betax_rhs, f_betay_rhs, f_betaz_rhs,
        f_dtSfx_rhs, f_dtSfy_rhs, f_dtSfz_rhs
    );

    apply_transport_and_dissipation_kernel<float><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], symmetry, lev, eps, X, Y, Z,
        f_betax, f_betay, f_betaz, f_dxx, f_gxy, f_gxz, f_dyy, f_gyz, f_dzz,
        f_Axx, f_Axy, f_Axz, f_Ayy, f_Ayz, f_Azz, f_chi, f_trK, f_Gamx, f_Gamy, f_Gamz, f_Lap, f_dtSfx, f_dtSfy, f_dtSfz,
        f_gxx_rhs, f_gxy_rhs, f_gxz_rhs, f_gyy_rhs, f_gyz_rhs, f_gzz_rhs,
        f_Axx_rhs, f_Axy_rhs, f_Axz_rhs, f_Ayy_rhs, f_Ayz_rhs, f_Azz_rhs,
        f_chi_rhs, f_trK_rhs, f_Gamx_rhs, f_Gamy_rhs, f_Gamz_rhs, f_Lap_rhs,
        f_betax_rhs, f_betay_rhs, f_betaz_rhs, f_dtSfx_rhs, f_dtSfy_rhs, f_dtSfz_rhs
    );

    if (co == 0) {
        evaluate_bssn_constraints_kernel<float><<<grid, block, 0, stream>>>(
            ex[0], ex[1], ex[2], symmetry, lev, X, Y, Z,
            f_chi, f_trK, f_dxx, f_gxy, f_gxz, f_dyy, f_gyz, f_dzz,
            f_Axx, f_Axy, f_Axz, f_Ayy, f_Ayz, f_Azz, f_Gamx, f_Gamy, f_Gamz,
            zero_matter, zero_matter, zero_matter, zero_matter,
            f_Gamxxx, f_Gamxxy, f_Gamxxz, f_Gamxyy, f_Gamxyz, f_Gamxzz,
            f_Gamyxx, f_Gamyxy, f_Gamyxz, f_Gamyyy, f_Gamyyz, f_Gamyzz,
            f_Gamzxx, f_Gamzxy, f_Gamzxz, f_Gamzyy, f_Gamzyz, f_Gamzzz,
            f_Rxx, f_Rxy, f_Rxz, f_Ryy, f_Ryz, f_Rzz,
            constraint_fields[0], constraint_fields[1], constraint_fields[2], constraint_fields[3],
            constraint_fields[4], constraint_fields[5], constraint_fields[6]
        );
    }

    std::vector<double*> rhs_targets = {
        chi_rhs, trK_rhs,
        gxx_rhs, gxy_rhs, gxz_rhs, gyy_rhs, gyz_rhs, gzz_rhs,
        Axx_rhs, Axy_rhs, Axz_rhs, Ayy_rhs, Ayz_rhs, Azz_rhs,
        Gamx_rhs, Gamy_rhs, Gamz_rhs, Lap_rhs,
        betax_rhs, betay_rhs, betaz_rhs,
        dtSfx_rhs, dtSfy_rhs, dtSfz_rhs
    };
    bssn_convert_float_to_double_fields(stream, rhs_fields, rhs_targets, point_count, scratch);

    std::vector<double*> aux_targets = {
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz
    };
    bssn_convert_float_to_double_fields(stream, aux_fields, aux_targets, point_count, scratch);

    if (co == 0) {
        std::vector<double*> constraint_targets = {
            ham_Res, movx_Res, movy_Res, movz_Res,
            Gmx_Res, Gmy_Res, Gmz_Res
        };
        bssn_convert_float_to_double_fields(stream, constraint_fields, constraint_targets, point_count, scratch);
    }
}

void gpu_compute_rhs_bssn_launch(
    cudaStream_t &stream, int* ex, double T, double* X, double* Y, double* Z,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Gamx, double* Gamy, double* Gamz, double* Lap,
    double* betax, double* betay, double* betaz,
    double* dtSfx, double* dtSfy, double* dtSfz,
    double* chi_rhs, double* trK_rhs,
    double* gxx_rhs, double* gxy_rhs, double* gxz_rhs, double* gyy_rhs, double* gyz_rhs, double* gzz_rhs,
    double* Axx_rhs, double* Axy_rhs, double* Axz_rhs, double* Ayy_rhs, double* Ayz_rhs, double* Azz_rhs,
    double* Gamx_rhs, double* Gamy_rhs, double* Gamz_rhs, double* Lap_rhs,
    double* betax_rhs, double* betay_rhs, double* betaz_rhs,
    double* dtSfx_rhs, double* dtSfy_rhs, double* dtSfz_rhs,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Sxx, double* Sxy, double* Sxz, double* Syy, double* Syz, double* Szz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz, double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz, double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz, double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* ham_Res, double* movx_Res, double* movy_Res, double* movz_Res,
    double* Gmx_Res, double* Gmy_Res, double* Gmz_Res,
    int symmetry, int lev, double eps, int co
) {
#if ABEtype == 0
    if (bssn_fp32_rhs_enabled()) {
        gpu_compute_rhs_bssn_fp32_launch(
            stream,
            ex, T, X, Y, Z,
            chi, trK,
            dxx, gxy, gxz, dyy, gyz, dzz,
            Axx, Axy, Axz, Ayy, Ayz, Azz,
            Gamx, Gamy, Gamz, Lap,
            betax, betay, betaz,
            dtSfx, dtSfy, dtSfz,
            chi_rhs, trK_rhs,
            gxx_rhs, gxy_rhs, gxz_rhs, gyy_rhs, gyz_rhs, gzz_rhs,
            Axx_rhs, Axy_rhs, Axz_rhs, Ayy_rhs, Ayz_rhs, Azz_rhs,
            Gamx_rhs, Gamy_rhs, Gamz_rhs, Lap_rhs,
            betax_rhs, betay_rhs, betaz_rhs,
            dtSfx_rhs, dtSfy_rhs, dtSfz_rhs,
            rho, Sx, Sy, Sz,
            Sxx, Sxy, Sxz, Syy, Syz, Szz,
            Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
            Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
            Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
            Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
            ham_Res, movx_Res, movy_Res, movz_Res,
            Gmx_Res, Gmy_Res, Gmz_Res,
            symmetry, lev, eps, co
        );
        return;
    }
#endif

    gpu_compute_rhs_bssn_launch(
        stream,
        ex, T, X, Y, Z,
        chi, trK,
        dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz,
        Gamx, Gamy, Gamz, Lap,
        betax, betay, betaz,
        dtSfx, dtSfy, dtSfz,
        chi_rhs, trK_rhs,
        gxx_rhs, gxy_rhs, gxz_rhs, gyy_rhs, gyz_rhs, gzz_rhs,
        Axx_rhs, Axy_rhs, Axz_rhs, Ayy_rhs, Ayz_rhs, Azz_rhs,
        Gamx_rhs, Gamy_rhs, Gamz_rhs, Lap_rhs,
        betax_rhs, betay_rhs, betaz_rhs,
        dtSfx_rhs, dtSfy_rhs, dtSfz_rhs,
        rho, Sx, Sy, Sz,
        Sxx, Sxy, Sxz, Syy, Syz, Szz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
        ham_Res, movx_Res, movy_Res, movz_Res,
        Gmx_Res, Gmy_Res, Gmz_Res,
        symmetry, lev, eps, co, false
    );
}

void gpu_compute_rhs_bssn_escalar_launch(
    cudaStream_t &stream, int* ex, double T, double* X, double* Y, double* Z,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Gamx, double* Gamy, double* Gamz, double* Lap,
    double* betax, double* betay, double* betaz,
    double* dtSfx, double* dtSfy, double* dtSfz,
    double* Sphi, double* Spi,
    double* chi_rhs, double* trK_rhs,
    double* gxx_rhs, double* gxy_rhs, double* gxz_rhs, double* gyy_rhs, double* gyz_rhs, double* gzz_rhs,
    double* Axx_rhs, double* Axy_rhs, double* Axz_rhs, double* Ayy_rhs, double* Ayz_rhs, double* Azz_rhs,
    double* Gamx_rhs, double* Gamy_rhs, double* Gamz_rhs, double* Lap_rhs,
    double* betax_rhs, double* betay_rhs, double* betaz_rhs,
    double* dtSfx_rhs, double* dtSfy_rhs, double* dtSfz_rhs,
    double* Sphi_rhs, double* Spi_rhs,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Sxx, double* Sxy, double* Sxz, double* Syy, double* Syz, double* Szz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz, double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz, double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz, double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* ham_Res, double* movx_Res, double* movy_Res, double* movz_Res,
    double* Gmx_Res, double* Gmy_Res, double* Gmz_Res,
    int symmetry, int lev, double eps, int co, double a2
) {
    dim3 block = bssn_block_shape();
    dim3 grid = bssn_grid_shape(ex, block);

    evolve_bssn_escalar_matter_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], symmetry, lev,
        X, Y, Z,
        chi, trK,
        dxx, gxy, gxz, dyy, gyz, dzz,
        Gamx, Gamy, Gamz, Lap,
        betax, betay, betaz,
        Sphi, Spi,
        Sphi_rhs, Spi_rhs,
        rho, Sx, Sy, Sz,
        Sxx, Sxy, Sxz, Syy, Syz, Szz,
        a2,
        eps
    );

    gpu_compute_rhs_bssn_launch(
        stream,
        ex, T, X, Y, Z,
        chi, trK,
        dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz,
        Gamx, Gamy, Gamz, Lap,
        betax, betay, betaz,
        dtSfx, dtSfy, dtSfz,
        chi_rhs, trK_rhs,
        gxx_rhs, gxy_rhs, gxz_rhs, gyy_rhs, gyz_rhs, gzz_rhs,
        Axx_rhs, Axy_rhs, Axz_rhs, Ayy_rhs, Ayz_rhs, Azz_rhs,
        Gamx_rhs, Gamy_rhs, Gamz_rhs, Lap_rhs,
        betax_rhs, betay_rhs, betaz_rhs,
        dtSfx_rhs, dtSfy_rhs, dtSfz_rhs,
        rho, Sx, Sy, Sz,
        Sxx, Sxy, Sxz, Syy, Syz, Szz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
        ham_Res, movx_Res, movy_Res, movz_Res,
        Gmx_Res, Gmy_Res, Gmz_Res,
        symmetry, lev, eps, co, true
    );
}

void gpu_compute_constraint_fr_launch(
    cudaStream_t &stream,
    int* ex,
    double* chi,
    double* trK,
    double* rho,
    double* Sphi,
    double* dxx,
    double* gxy,
    double* gxz,
    double* dyy,
    double* gyz,
    double* dzz,
    double* Axx,
    double* Axy,
    double* Axz,
    double* Ayy,
    double* Ayz,
    double* Azz,
    double* Rxx,
    double* Rxy,
    double* Rxz,
    double* Ryy,
    double* Ryz,
    double* Rzz,
    double* Sxx,
    double* Sxy,
    double* Sxz,
    double* Syy,
    double* Syz,
    double* Szz,
    double* Cons_fR,
    double a2
) {
    dim3 block = bssn_block_shape();
    dim3 grid = bssn_grid_shape(ex, block);

    compute_constraint_fr_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2],
        chi, trK, rho, Sphi,
        dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
        Sxx, Sxy, Sxz, Syy, Syz, Szz,
        Cons_fR, a2
    );
}

void gpu_compute_rhs_bssn_em_launch(
    cudaStream_t &stream, int* ex, double T, double* X, double* Y, double* Z,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Gamx, double* Gamy, double* Gamz, double* Lap,
    double* betax, double* betay, double* betaz,
    double* dtSfx, double* dtSfy, double* dtSfz,
    double* Ex, double* Ey, double* Ez,
    double* Bx, double* By, double* Bz,
    double* Kpsi, double* Kphi,
    double* Jx, double* Jy, double* Jz, double* qchar,
    double* chi_rhs, double* trK_rhs,
    double* gxx_rhs, double* gxy_rhs, double* gxz_rhs, double* gyy_rhs, double* gyz_rhs, double* gzz_rhs,
    double* Axx_rhs, double* Axy_rhs, double* Axz_rhs, double* Ayy_rhs, double* Ayz_rhs, double* Azz_rhs,
    double* Gamx_rhs, double* Gamy_rhs, double* Gamz_rhs, double* Lap_rhs,
    double* betax_rhs, double* betay_rhs, double* betaz_rhs,
    double* dtSfx_rhs, double* dtSfy_rhs, double* dtSfz_rhs,
    double* Ex_rhs, double* Ey_rhs, double* Ez_rhs,
    double* Bx_rhs, double* By_rhs, double* Bz_rhs,
    double* Kpsi_rhs, double* Kphi_rhs,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Sxx, double* Sxy, double* Sxz, double* Syy, double* Syz, double* Szz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz, double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz, double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz, double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* ham_Res, double* movx_Res, double* movy_Res, double* movz_Res,
    double* Gmx_Res, double* Gmy_Res, double* Gmz_Res,
    int symmetry, int lev, double eps, int co
) {
    dim3 block = bssn_block_shape();
    dim3 grid = bssn_grid_shape(ex, block);

    evolve_bssn_em_matter_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], symmetry, lev,
        X, Y, Z, chi, trK,
        dxx, gxy, gxz, dyy, gyz, dzz,
        Lap, betax, betay, betaz,
        Ex, Ey, Ez, Bx, By, Bz, Kpsi, Kphi,
        Jx, Jy, Jz, qchar,
        Ex_rhs, Ey_rhs, Ez_rhs, Bx_rhs, By_rhs, Bz_rhs, Kpsi_rhs, Kphi_rhs,
        rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz,
        eps
    );

    gpu_compute_rhs_bssn_launch(
        stream,
        ex, T, X, Y, Z,
        chi, trK,
        dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz,
        Gamx, Gamy, Gamz, Lap,
        betax, betay, betaz,
        dtSfx, dtSfy, dtSfz,
        chi_rhs, trK_rhs,
        gxx_rhs, gxy_rhs, gxz_rhs, gyy_rhs, gyz_rhs, gzz_rhs,
        Axx_rhs, Axy_rhs, Axz_rhs, Ayy_rhs, Ayz_rhs, Azz_rhs,
        Gamx_rhs, Gamy_rhs, Gamz_rhs, Lap_rhs,
        betax_rhs, betay_rhs, betaz_rhs,
        dtSfx_rhs, dtSfy_rhs, dtSfz_rhs,
        rho, Sx, Sy, Sz,
        Sxx, Sxy, Sxz, Syy, Syz, Szz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
        ham_Res, movx_Res, movy_Res, movz_Res,
        Gmx_Res, Gmy_Res, Gmz_Res,
        symmetry, lev, eps, co, false
    );
}

void gpu_compute_rhs_z4c_launch(
    cudaStream_t &stream, int* ex, double T, double* X, double* Y, double* Z,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Gamx, double* Gamy, double* Gamz, double* Lap,
    double* betax, double* betay, double* betaz,
    double* dtSfx, double* dtSfy, double* dtSfz,
    double* TZ,
    double* chi_rhs, double* trK_rhs,
    double* gxx_rhs, double* gxy_rhs, double* gxz_rhs, double* gyy_rhs, double* gyz_rhs, double* gzz_rhs,
    double* Axx_rhs, double* Axy_rhs, double* Axz_rhs, double* Ayy_rhs, double* Ayz_rhs, double* Azz_rhs,
    double* Gamx_rhs, double* Gamy_rhs, double* Gamz_rhs, double* Lap_rhs,
    double* betax_rhs, double* betay_rhs, double* betaz_rhs,
    double* dtSfx_rhs, double* dtSfy_rhs, double* dtSfz_rhs,
    double* TZ_rhs,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Sxx, double* Sxy, double* Sxz, double* Syy, double* Syz, double* Szz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz, double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz, double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz, double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* ham_Res, double* movx_Res, double* movy_Res, double* movz_Res,
    double* Gmx_Res, double* Gmy_Res, double* Gmz_Res,
    int symmetry, int lev, double eps, int co
) {
    dim3 block = bssn_block_shape();
    dim3 grid = bssn_grid_shape(ex, block);

    trace_metric_det_summary(stream, "z4c_entry_metric", ex, X, Y, Z, lev, co,
                             dxx, gxy, gxz, dyy, gyz, dzz);

    assemble_conformal_geometry_kernel<double><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], symmetry, lev, X, Y, Z,
        dxx, gxy, gxz, dyy, gyz, dzz, Gamx, Gamy, Gamz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz
    );

    evolve_bssn_metric_fields_kernel<double><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], symmetry, X, Y, Z,
        chi, trK, dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz,
        Lap, betax, betay, betaz,
        chi_rhs, gxx_rhs, gxy_rhs, gxz_rhs, gyy_rhs, gyz_rhs, gzz_rhs
    );

    evolve_bssn_curvature_fields_kernel<double><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], symmetry, X, Y, Z,
        chi, trK, dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz,
        Lap, betax, betay, betaz,
        rho, Sxx, Sxy, Sxz, Syy, Syz, Szz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
        trK_rhs,
        Axx_rhs, Axy_rhs, Axz_rhs, Ayy_rhs, Ayz_rhs, Azz_rhs,
        false
    );

    evolve_bssn_gauge_fields_kernel<double><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], symmetry, X, Y, Z,
        chi, trK, dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz,
        Gamx, Gamy, Gamz, Lap, betax, betay, betaz, dtSfx, dtSfy, dtSfz,
        Sx, Sy, Sz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Gamx_rhs, Gamy_rhs, Gamz_rhs, Lap_rhs,
        betax_rhs, betay_rhs, betaz_rhs,
        dtSfx_rhs, dtSfy_rhs, dtSfz_rhs
    );

    evaluate_bssn_constraints_kernel<double><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], symmetry, lev, X, Y, Z,
        chi, trK, dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz, Gamx, Gamy, Gamz, rho, Sx, Sy, Sz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
        ham_Res, movx_Res, movy_Res, movz_Res,
        Gmx_Res, Gmy_Res, Gmz_Res
    );

    apply_z4c_algebraic_corrections_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], symmetry, X, Y, Z,
        chi, trK, TZ, dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz, Lap,
        rho, Sxx, Sxy, Sxz, Syy, Syz, Szz,
        chi_rhs, trK_rhs,
        Axx_rhs, Axy_rhs, Axz_rhs, Ayy_rhs, Ayz_rhs, Azz_rhs,
        TZ_rhs, ham_Res, movx_Res, movy_Res, movz_Res
    );

    apply_z4c_gauge_corrections_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], symmetry, X, Y, Z, TZ,
        dxx, gxy, gxz, dyy, gyz, dzz,
        Gamx, Gamy, Gamz, Lap, betax, betay, betaz,
        Gamx_rhs, Gamy_rhs, Gamz_rhs,
        betax_rhs, betay_rhs, betaz_rhs,
        dtSfx_rhs, dtSfy_rhs, dtSfz_rhs
    );

    apply_transport_and_dissipation_kernel<double><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], symmetry, lev, eps, X, Y, Z,
        betax, betay, betaz, dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz, chi, trK, Gamx, Gamy, Gamz, Lap, dtSfx, dtSfy, dtSfz,
        gxx_rhs, gxy_rhs, gxz_rhs, gyy_rhs, gyz_rhs, gzz_rhs,
        Axx_rhs, Axy_rhs, Axz_rhs, Ayy_rhs, Ayz_rhs, Azz_rhs,
        chi_rhs, trK_rhs, Gamx_rhs, Gamy_rhs, Gamz_rhs, Lap_rhs,
        betax_rhs, betay_rhs, betaz_rhs, dtSfx_rhs, dtSfy_rhs, dtSfz_rhs
    );

    apply_z4c_transport_damping_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], symmetry, eps, X, Y, Z,
        TZ, Lap, betax, betay, betaz, Gmx_Res, Gmy_Res, Gmz_Res,
        trK_rhs, Gamx_rhs, Gamy_rhs, Gamz_rhs, TZ_rhs
    );

    finalize_z4c_constraint_outputs_kernel<<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2], symmetry, X, Y, Z,
        chi, trK, dxx, gxy, gxz, dyy, gyz, dzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz,
        Gamx, Gamy, Gamz, Sx, Sy, Sz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        movx_Res, movy_Res, movz_Res,
        Gmx_Res, Gmy_Res, Gmz_Res
    );

    const RhsTraceField z4c_rhs_fields[] = {
        {"gxx_rhs", gxx_rhs}, {"gxy_rhs", gxy_rhs}, {"gxz_rhs", gxz_rhs},
        {"gyy_rhs", gyy_rhs}, {"gyz_rhs", gyz_rhs}, {"gzz_rhs", gzz_rhs},
        {"Axx_rhs", Axx_rhs}, {"Axy_rhs", Axy_rhs}, {"Axz_rhs", Axz_rhs},
        {"Ayy_rhs", Ayy_rhs}, {"Ayz_rhs", Ayz_rhs}, {"Azz_rhs", Azz_rhs},
        {"phi_rhs", chi_rhs}, {"trK_rhs", trK_rhs}, {"TZ_rhs", TZ_rhs},
        {"Gamx_rhs", Gamx_rhs}, {"Gamy_rhs", Gamy_rhs}, {"Gamz_rhs", Gamz_rhs}
    };
    trace_rhs_nonfinite_summary(stream, "z4c_after_transport", ex, X, Y, Z, lev, co,
                                z4c_rhs_fields, sizeof(z4c_rhs_fields) / sizeof(z4c_rhs_fields[0]));
}
