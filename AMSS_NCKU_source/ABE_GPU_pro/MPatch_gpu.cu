#include "MPatch.h"

#include <cuda_runtime.h>
#include "fmisc.h"
#include "macrodef.h"
#include "gpu_info.h"

using namespace std;

// ============================================================================
// normalize_shellf_kernel + gpu_normalize_shellf_launch
// ============================================================================

__global__ void normalize_shellf_kernel(
    int NN, int num_var,
    double* d_Shellf,
    const int* d_Weight,
    int* d_err_flag
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= NN) return;

    int w = d_Weight[i];
    if (w > 1) {
        for (int j = 0; j < num_var; j++) {
            d_Shellf[j + i * num_var] /= (double)w;
        }
        // printf("\33[1;33mWARNING: normalize_shellf_kernel meets multiple weight at point %d\33[0m\n", i);
    } else if (w == 0) {
        *d_err_flag = 1;
    }
}

void gpu_normalize_shellf_launch(
    cudaStream_t stream,
    int NN, int num_var,
    double* d_Shellf,
    const int* d_Weight,
    int* d_err_flag
) {
    int blockSize = 256;
    int gridSize = (NN + blockSize - 1) / blockSize;
    normalize_shellf_kernel<<<gridSize, blockSize, 0, stream>>>(
        NN, num_var, d_Shellf, d_Weight, d_err_flag
    );
}

// ============================================================================
// find_maximum_kernel + gpu_find_maximum_launch
// ============================================================================

__global__ void find_maximum_kernel(
    int shape_0, int shape_1, int shape_2,
    int pshape_0, int pshape_1, int pshape_2,
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_fun,
    int llb_0, int llb_1, int llb_2,
    int uub_0, int uub_1, int uub_2,
    double* d_max_val, double* d_max_pos_x, double* d_max_pos_y, double* d_max_pos_z
) {
    // Grid-stride loop with single block of 256 threads
    double thread_max = 0.0;
    int thread_max_i = llb_0;
    int thread_max_j = llb_1;
    int thread_max_k = llb_2;

    int total = (shape_0 - uub_0 - llb_0) *
                (shape_1 - uub_1 - llb_1) *
                (shape_2 - uub_2 - llb_2);

    for (int idx = threadIdx.x; idx < total; idx += blockDim.x) {
        int rem = idx;
        int di = rem % (shape_0 - uub_0 - llb_0);
        rem /= (shape_0 - uub_0 - llb_0);
        int dj = rem % (shape_1 - uub_1 - llb_1);
        rem /= (shape_1 - uub_1 - llb_1);
        int dk = rem;

        int i = di + llb_0;
        int j = dj + llb_1;
        int k = dk + llb_2;

        double val = fabs(d_fun[i + j * pshape_0 + k * pshape_0 * pshape_1]);
        if (val > thread_max) {
            thread_max = val;
            thread_max_i = i;
            thread_max_j = j;
            thread_max_k = k;
        }
    }

    // Block-level reduction using shared memory
    __shared__ double smem_val[256];
    __shared__ int smem_i[256];
    __shared__ int smem_j[256];
    __shared__ int smem_k[256];

    smem_val[threadIdx.x] = thread_max;
    smem_i[threadIdx.x] = thread_max_i;
    smem_j[threadIdx.x] = thread_max_j;
    smem_k[threadIdx.x] = thread_max_k;
    __syncthreads();

    for (int offset = blockDim.x / 2; offset > 0; offset >>= 1) {
        if (threadIdx.x < offset) {
            if (smem_val[threadIdx.x + offset] > smem_val[threadIdx.x]) {
                smem_val[threadIdx.x] = smem_val[threadIdx.x + offset];
                smem_i[threadIdx.x] = smem_i[threadIdx.x + offset];
                smem_j[threadIdx.x] = smem_j[threadIdx.x + offset];
                smem_k[threadIdx.x] = smem_k[threadIdx.x + offset];
            }
        }
        __syncthreads();
    }

    // Thread 0 writes the result
    if (threadIdx.x == 0) {
        *d_max_val = smem_val[0];
        *d_max_pos_x = d_X[smem_i[0]];
        *d_max_pos_y = d_Y[smem_j[0]];
        *d_max_pos_z = d_Z[smem_k[0]];
    }
}

void gpu_find_maximum_launch(
    cudaStream_t stream,
    int shape_0, int shape_1, int shape_2,
    int pshape_0, int pshape_1, int pshape_2,
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_fun,
    int llb_0, int llb_1, int llb_2,
    int uub_0, int uub_1, int uub_2,
    double* d_max_val, double* d_max_pos_x, double* d_max_pos_y, double* d_max_pos_z
) {
    find_maximum_kernel<<<1, 256, 0, stream>>>(
        shape_0, shape_1, shape_2,
        pshape_0, pshape_1, pshape_2,
        d_X, d_Y, d_Z, d_fun,
        llb_0, llb_1, llb_2,
        uub_0, uub_1, uub_2,
        d_max_val, d_max_pos_x, d_max_pos_y, d_max_pos_z
    );
}

// ============================================================================
// Interp_Points_GPU (MPI_COMM_WORLD) -- delegates to MPI_Comm version
// ============================================================================

void Patch::Interp_Points_GPU(
    cudaStream_t stream,
    MyList<var> *VarList,
    int NN, double **d_XX,
    double *d_Shellf, int Symmetry)
{
    Interp_Points_GPU(stream, VarList, NN, d_XX, d_Shellf, Symmetry, MPI_COMM_WORLD);
}

// ============================================================================
// Interp_Points_GPU (MPI_Comm version)
// ============================================================================

void Patch::Interp_Points_GPU(
    cudaStream_t stream,
    MyList<var> *VarList,
    int NN, double **d_XX,
    double *d_Shellf, int Symmetry, MPI_Comm Comm_here)
{
    int myrank, lmyrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_rank(Comm_here, &lmyrank);

    int ordn = 2 * ghost_width;
    MyList<var> *varl = VarList;
    int num_var = 0;
    while (varl) {
        num_var++;
        varl = varl->next;
    }

    double *DH = new double[dim];
    for (int i = 0; i < dim; i++) DH[i] = getdX(i);
    double *llb = new double[dim];
    double *uub = new double[dim];

    double *d_local_shellf = GPUManager::getInstance().allocate_comm_buffer<double>(NN * num_var);
    int *d_local_weight = GPUManager::getInstance().allocate_comm_buffer<int>(NN);
    GPUManager::getInstance().memset_device_memory(d_local_shellf, NN * num_var);
    GPUManager::getInstance().memset_device_memory(d_local_weight, NN);

    int *d_global_weight = GPUManager::getInstance().allocate_comm_buffer<int>(NN, stream);
    GPUManager::getInstance().memset_device_memory(d_global_weight, NN, stream);
    GPUManager::getInstance().synchronize_memory();

    MyList<Block> *Bp = blb;
    while (Bp) {
        Block *BP = Bp->data;
        if (myrank == BP->rank) {
            for (int i = 0; i < dim; i++) {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
                llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + (ghost_width - 0.5) * DH[i];
                uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - (ghost_width - 0.5) * DH[i];
#else
#ifdef Cell
                llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + ghost_width * DH[i];
                uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - ghost_width * DH[i];
#else
#error Neither Vertex nor Cell is defined
#endif
#endif
            }

            int shape_0 = BP->shape[0];
            int shape_1 = (dim > 1) ? BP->shape[1] : 1;
            int shape_2 = (dim > 2) ? BP->shape[2] : 1;

            int shape_pad_0 = BP->shape_pad[0];
            int shape_pad_1 = (dim > 1) ? BP->shape_pad[1] : 1;
            int shape_pad_2 = (dim > 2) ? BP->shape_pad[2] : 1;

            double llb_0 = llb[0], llb_1 = (dim > 1) ? llb[1] : 0.0, llb_2 = (dim > 2) ? llb[2] : 0.0;
            double uub_0 = uub[0], uub_1 = (dim > 1) ? uub[1] : 0.0, uub_2 = (dim > 2) ? uub[2] : 0.0;

            double DH_0 = DH[0];
            double DH_1 = (dim > 1) ? DH[1] : 0.0;
            double DH_2 = (dim > 2) ? DH[2] : 0.0;
            varl = VarList;
            int k = 0;
            while (varl) {
                gpu_global_interp_launch(
                    BP->stream,
                    NN, dim,
                    d_XX[0], d_XX[1], d_XX[2],
                    shape_0, shape_1, shape_2,
                    shape_pad_0, shape_pad_1, shape_pad_2,
                    BP->d_X[0], BP->d_X[1], BP->d_X[2],
                    BP->d_fgfs[varl->data->sgfn],
                    llb_0, llb_1, llb_2,
                    uub_0, uub_1, uub_2,
                    DH_0, DH_1, DH_2,
                    ordn, varl->data->SoA[0], varl->data->SoA[1], varl->data->SoA[2],
                    Symmetry, k, num_var, d_local_shellf, d_local_weight
                );
                varl = varl->next;
                k++;
            }
        }
        if (Bp == ble) break;
        Bp = Bp->next;
    }

    GPUManager::getInstance().synchronize_all();

#if MPI_CUDA_AWARE
    MPI_Allreduce(d_local_shellf, d_Shellf, NN * num_var, MPI_DOUBLE, MPI_SUM, Comm_here);
    MPI_Allreduce(d_local_weight, d_global_weight, NN, MPI_INT, MPI_SUM, Comm_here);
#else
    double *h_local_shellf = new double[NN * num_var];
    int    *h_local_weight = new int[NN];
    GPUManager::getInstance().sync_to_cpu(h_local_shellf, d_local_shellf, NN * num_var);
    GPUManager::getInstance().sync_to_cpu(h_local_weight, d_local_weight, NN);

    double *h_global_shellf = new double[NN * num_var];
    int    *h_global_weight = new int[NN];

    GPUManager::getInstance().synchronize_memory();

    MPI_Allreduce(h_local_shellf, h_global_shellf, NN * num_var, MPI_DOUBLE, MPI_SUM, Comm_here);
    MPI_Allreduce(h_local_weight, h_global_weight, NN, MPI_INT, MPI_SUM, Comm_here);

    GPUManager::getInstance().sync_to_gpu(h_global_shellf, d_Shellf, NN * num_var, stream);
    GPUManager::getInstance().sync_to_gpu(h_global_weight, d_global_weight, NN, stream);

    delete[] h_local_shellf; delete[] h_local_weight;
    delete[] h_global_shellf; delete[] h_global_weight;
#endif

    int *d_err_flag = GPUManager::getInstance().allocate_device_memory<int>(1, stream);
    GPUManager::getInstance().memset_device_memory(d_err_flag, 1, stream);
    gpu_normalize_shellf_launch(stream, NN, num_var, d_Shellf, d_global_weight, d_err_flag);

    int h_err_flag = 0;
    GPUManager::getInstance().sync_to_cpu(&h_err_flag, d_err_flag, 1, stream);
    GPUManager::synchronize_stream(stream);

    if (h_err_flag > 0) {
        if (lmyrank == 0) {
            int *h_temp_weight = new int[NN];
            GPUManager::getInstance().sync_to_cpu(h_temp_weight, d_global_weight, NN);

            double *h_XX[3];
            for (int i = 0; i < dim; i++) {
                h_XX[i] = new double[NN];
                GPUManager::getInstance().sync_to_cpu(h_XX[i], d_XX[i], NN);
            }
            GPUManager::synchronize_stream(stream);

            for (int i = 0; i < NN; i++) {
                if (h_temp_weight[i] > 1) {
                    cout << "WARNING: Patch::Interp_Points meets multiple weight at point " << i << endl;
                }
                else if (h_temp_weight[i] == 0) {
                    cout << "ERROR: Patch::Interp_Points fails to find point (";
                    for (int j = 0; j < dim; j++) {
                        cout << h_XX[j][i];
                        if (j < dim - 1) cout << ",";
                        else cout << ") ! Out of domain bounds." << endl;
                    }
                }
            }
            delete[] h_temp_weight;
            for (int i = 0; i < dim; i++) delete[] h_XX[i];
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    delete[] DH; delete[] llb; delete[] uub;
    GPUManager::getInstance().free_comm_buffer(d_local_shellf, NN * num_var);
    GPUManager::getInstance().free_comm_buffer(d_local_weight, NN);
    GPUManager::getInstance().free_comm_buffer(d_global_weight, NN, stream);
    GPUManager::getInstance().free_device_memory(d_err_flag, 1, stream);
}

// ============================================================================
// Interp_N_Points_GPU (existing, unchanged)
// ============================================================================

bool Patch::Interp_N_Points_GPU(
    MyList<var> *VarList,
    int NN, double *d_XX_0, double *d_XX_1, double *d_XX_2,
    double *d_shellf, int *d_weight, int Symmetry
) {
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int ordn = 2 * ghost_width;
    MyList<var> *varl = VarList;
    int num_var = 0;
    while (varl) {
        num_var++;
        varl = varl->next;
    }

    double DH[3];
    for (int i = 0; i < dim; i++) {
        DH[i] = getdX(i);
    }

    MyList<Block> *Bp = blb;
    while (Bp) {
        Block *BP = Bp->data;
        if (myrank == BP->rank) {
            double llb[3], uub[3];
            for (int i = 0; i < dim; i++) {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
                llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + (ghost_width - 0.5) * DH[i];
                uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - (ghost_width - 0.5) * DH[i];
#else
#ifdef Cell
                llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + ghost_width * DH[i];
                uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - ghost_width * DH[i];
#else
#error Neither Vertex nor Cell is defined
#endif
#endif
            }

            double DH_0 = DH[0];
            double DH_1 = (dim > 1) ? DH[1] : 0.0;
            double DH_2 = (dim > 2) ? DH[2] : 0.0;
            varl = VarList;
            int k = 0;
            while (varl) {
                gpu_global_interp_launch(
                    BP->stream, NN, dim,
                    d_XX_0, d_XX_1, d_XX_2,
                    BP->shape[0], BP->shape[1], BP->shape[2],
                    BP->shape_pad[0], BP->shape_pad[1], BP->shape_pad[2],
                    BP->d_X[0], BP->d_X[1], BP->d_X[2],
                    BP->d_fgfs[varl->data->sgfn],
                    llb[0], llb[1], llb[2], uub[0], uub[1], uub[2],
                    DH_0, DH_1, DH_2,
                    ordn, varl->data->SoA[0], varl->data->SoA[1], varl->data->SoA[2],
                    Symmetry, k, num_var, d_shellf, d_weight
                );
                varl = varl->next;
                k++;
            }
        }
        Bp = Bp->next;
    }

    GPUManager::getInstance().synchronize_all();

#if MPI_CUDA_AWARE
    MPI_Allreduce(MPI_IN_PLACE, d_shellf, NN * num_var, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, d_weight, NN, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
    double *h_shellf_local = new double[NN * num_var];
    double *h_shellf_global = new double[NN * num_var];
    int *h_weight_local = new int[NN];
    int *h_weight_global = new int[NN];

    GPUManager::getInstance().sync_to_cpu(h_shellf_local, d_shellf, NN * num_var);
    GPUManager::getInstance().sync_to_cpu(h_weight_local, d_weight, NN);
    GPUManager::getInstance().synchronize_memory();

    MPI_Allreduce(h_shellf_local, h_shellf_global, NN * num_var, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(h_weight_local, h_weight_global, NN, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    GPUManager::getInstance().sync_to_gpu(h_shellf_global, d_shellf, NN * num_var);
    GPUManager::getInstance().sync_to_gpu(h_weight_global, d_weight, NN);
    GPUManager::getInstance().synchronize_memory();

    delete[] h_shellf_local;
    delete[] h_shellf_global;
    delete[] h_weight_local;
    delete[] h_weight_global;
#endif

    return true;
}

// ============================================================================
// Interp_ONE_Point_GPU (MPI_COMM_WORLD) -- delegates to MPI_Comm version
// ============================================================================

bool Patch::Interp_ONE_Point_GPU(
    cudaStream_t stream,
    MyList<var> *VarList,
    double *d_XX_0, double *d_XX_1, double *d_XX_2,
    double *d_Shellf, int Symmetry)
{
    return Interp_ONE_Point_GPU(stream, VarList, d_XX_0, d_XX_1, d_XX_2, d_Shellf, Symmetry, MPI_COMM_WORLD);
}

// ============================================================================
// Interp_ONE_Point_GPU (MPI_Comm version)
// ============================================================================

bool Patch::Interp_ONE_Point_GPU(
    cudaStream_t stream,
    MyList<var> *VarList,
    double *d_XX_0, double *d_XX_1, double *d_XX_2,
    double *d_Shellf, int Symmetry, MPI_Comm Comm_here)
{
    int myrank, lmyrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_rank(Comm_here, &lmyrank);

    int ordn = 2 * ghost_width;
    MyList<var> *varl = VarList;
    int num_var = 0;
    while (varl) {
        num_var++;
        varl = varl->next;
    }

    double *DH = new double[dim];
    for (int i = 0; i < dim; i++) DH[i] = getdX(i);
    double *llb = new double[dim];
    double *uub = new double[dim];

    int NN = 1;
    double *d_local_shellf = GPUManager::getInstance().allocate_comm_buffer<double>(NN * num_var);
    int *d_local_weight = GPUManager::getInstance().allocate_comm_buffer<int>(NN);
    GPUManager::getInstance().memset_device_memory(d_local_shellf, NN * num_var);
    GPUManager::getInstance().memset_device_memory(d_local_weight, NN);
    GPUManager::getInstance().synchronize_memory();

    MyList<Block> *Bp = blb;
    while (Bp) {
        Block *BP = Bp->data;
        if (myrank == BP->rank) {
            for (int i = 0; i < dim; i++) {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
                llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + (ghost_width - 0.5) * DH[i];
                uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - (ghost_width - 0.5) * DH[i];
#else
#ifdef Cell
                llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + ghost_width * DH[i];
                uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - ghost_width * DH[i];
#else
#error Neither Vertex nor Cell is defined
#endif
#endif
            }

            int shape_0 = BP->shape[0];
            int shape_1 = (dim > 1) ? BP->shape[1] : 1;
            int shape_2 = (dim > 2) ? BP->shape[2] : 1;

            int shape_pad_0 = BP->shape_pad[0];
            int shape_pad_1 = (dim > 1) ? BP->shape_pad[1] : 1;
            int shape_pad_2 = (dim > 2) ? BP->shape_pad[2] : 1;

            double llb_0 = llb[0], llb_1 = (dim > 1) ? llb[1] : 0.0, llb_2 = (dim > 2) ? llb[2] : 0.0;
            double uub_0 = uub[0], uub_1 = (dim > 1) ? uub[1] : 0.0, uub_2 = (dim > 2) ? uub[2] : 0.0;

            double DH_0 = DH[0];
            double DH_1 = (dim > 1) ? DH[1] : 0.0;
            double DH_2 = (dim > 2) ? DH[2] : 0.0;
            varl = VarList;
            int k = 0;
            while (varl) {
                gpu_global_interp_launch(
                    BP->stream,
                    NN, dim,
                    d_XX_0, d_XX_1, d_XX_2,
                    shape_0, shape_1, shape_2,
                    shape_pad_0, shape_pad_1, shape_pad_2,
                    BP->d_X[0], BP->d_X[1], BP->d_X[2],
                    BP->d_fgfs[varl->data->sgfn],
                    llb_0, llb_1, llb_2,
                    uub_0, uub_1, uub_2,
                    DH_0, DH_1, DH_2,
                    ordn, varl->data->SoA[0], varl->data->SoA[1], varl->data->SoA[2],
                    Symmetry, k, num_var, d_local_shellf, d_local_weight
                );
                varl = varl->next;
                k++;
            }
        }
        if (Bp == ble) break;
        Bp = Bp->next;
    }

    GPUManager::getInstance().synchronize_all();

    // MPI reduction for single point
#if MPI_CUDA_AWARE
    MPI_Allreduce(d_local_shellf, d_Shellf, NN * num_var, MPI_DOUBLE, MPI_SUM, Comm_here);
    int *d_global_weight = GPUManager::getInstance().allocate_comm_buffer<int>(NN, stream);
    MPI_Allreduce(d_local_weight, d_global_weight, NN, MPI_INT, MPI_SUM, Comm_here);
#else
    double *h_local_shellf = new double[NN * num_var];
    int    *h_local_weight = new int[NN];
    GPUManager::getInstance().sync_to_cpu(h_local_shellf, d_local_shellf, NN * num_var);
    GPUManager::getInstance().sync_to_cpu(h_local_weight, d_local_weight, NN);

    double *h_global_shellf = new double[NN * num_var];
    int    *h_global_weight = new int[NN];

    GPUManager::getInstance().synchronize_memory();

    MPI_Allreduce(h_local_shellf, h_global_shellf, NN * num_var, MPI_DOUBLE, MPI_SUM, Comm_here);
    MPI_Allreduce(h_local_weight, h_global_weight, NN, MPI_INT, MPI_SUM, Comm_here);

    GPUManager::getInstance().sync_to_gpu(h_global_shellf, d_Shellf, NN * num_var, stream);
    int *d_global_weight = GPUManager::getInstance().allocate_comm_buffer<int>(NN, stream);
    GPUManager::getInstance().sync_to_gpu(h_global_weight, d_global_weight, NN, stream);

    delete[] h_local_shellf; delete[] h_local_weight;
    delete[] h_global_shellf; delete[] h_global_weight;
#endif

    // Normalize (point at block boundary may get weight > 1)
    int *d_err_flag = GPUManager::getInstance().allocate_device_memory<int>(1, stream);
    GPUManager::getInstance().memset_device_memory(d_err_flag, 1, stream);
    gpu_normalize_shellf_launch(stream, NN, num_var, d_Shellf, d_global_weight, d_err_flag);

    int h_err_flag = 0;
    GPUManager::getInstance().sync_to_cpu(&h_err_flag, d_err_flag, 1, stream);
    GPUManager::synchronize_stream(stream);

    int h_global_weight_val = 0;
    GPUManager::getInstance().sync_to_cpu(&h_global_weight_val, d_global_weight, 1, stream);
    GPUManager::synchronize_stream(stream);

    if (h_err_flag > 0 || h_global_weight_val == 0) {
        if (lmyrank == 0) {
            double h_XX[3];
            GPUManager::getInstance().sync_to_cpu(&h_XX[0], d_XX_0, 1, stream);
            GPUManager::getInstance().sync_to_cpu(&h_XX[1], d_XX_1, 1, stream);
            GPUManager::getInstance().sync_to_cpu(&h_XX[2], d_XX_2, 1, stream);
            GPUManager::synchronize_stream(stream);

            cout << "ERROR: Patch::Interp_ONE_Point_GPU fails to find point (";
            for (int j = 0; j < dim; j++) {
                cout << h_XX[j];
                if (j < dim - 1) cout << ",";
                else cout << ")";
            }
            cout << " on Patch (";
            for (int j = 0; j < dim; j++) {
                cout << bbox[j] << "+" << lli[j] * getdX(j);
                if (j < dim - 1) cout << ",";
                else cout << ")--";
            }
            cout << "(";
            for (int j = 0; j < dim; j++) {
                cout << bbox[dim + j] << "-" << uui[j] * getdX(j);
                if (j < dim - 1) cout << ",";
                else cout << ")" << endl;
            }
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    delete[] DH; delete[] llb; delete[] uub;
    GPUManager::getInstance().free_comm_buffer(d_local_shellf, NN * num_var);
    GPUManager::getInstance().free_comm_buffer(d_local_weight, NN);
    GPUManager::getInstance().free_comm_buffer(d_global_weight, NN, stream);
    GPUManager::getInstance().free_device_memory(d_err_flag, 1, stream);

    return (h_global_weight_val > 0);
}

// ============================================================================
// Find_Maximum_GPU (MPI_COMM_WORLD) -- delegates to MPI_Comm version
// ============================================================================

void Patch::Find_Maximum_GPU(
    cudaStream_t stream,
    MyList<var> *VarList,
    double *d_XX, double *d_Shellf)
{
    Find_Maximum_GPU(stream, VarList, d_XX, d_Shellf, MPI_COMM_WORLD);
}

// ============================================================================
// Find_Maximum_GPU (MPI_Comm version)
// ============================================================================

void Patch::Find_Maximum_GPU(
    cudaStream_t stream,
    MyList<var> *VarList,
    double *d_XX, double *d_Shellf, MPI_Comm Comm_here)
{
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    MyList<var> *varl = VarList;
    int num_var = 0;
    while (varl) {
        num_var++;
        varl = varl->next;
    }

    // Host arrays for reduction
    double *shellf = new double[num_var];
    double *xx = new double[dim * num_var];
    memset(shellf, 0, sizeof(double) * num_var);
    memset(xx, 0, sizeof(double) * dim * num_var);

    double *DH = new double[dim];
    for (int i = 0; i < dim; i++) DH[i] = getdX(i);

    // Device buffers (reused across blocks/variables)
    double *d_max_val = GPUManager::getInstance().allocate_device_memory<double>(1);
    double *d_max_pos_x = GPUManager::getInstance().allocate_device_memory<double>(1);
    double *d_max_pos_y = GPUManager::getInstance().allocate_device_memory<double>(1);
    double *d_max_pos_z = GPUManager::getInstance().allocate_device_memory<double>(1);

    MyList<Block> *Bp = blb;
    while (Bp) {
        Block *BP = Bp->data;

        if (myrank == BP->rank) {
            // Compute integer llb/uub for find_maximum (ghost zone counts)
            int illb[dim], iuub[dim];
            for (int i = 0; i < dim; i++) {
                illb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? lli[i] : ghost_width;
                iuub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? uui[i] : ghost_width;
            }

            int shape_0 = BP->shape[0];
            int shape_1 = (dim > 1) ? BP->shape[1] : 1;
            int shape_2 = (dim > 2) ? BP->shape[2] : 1;

            int pshape_0 = BP->shape_pad[0];
            int pshape_1 = (dim > 1) ? BP->shape_pad[1] : 1;
            int pshape_2 = (dim > 2) ? BP->shape_pad[2] : 1;

            int llb_0 = illb[0], llb_1 = (dim > 1) ? illb[1] : 0, llb_2 = (dim > 2) ? illb[2] : 0;
            int uub_0 = iuub[0], uub_1 = (dim > 1) ? iuub[1] : 0, uub_2 = (dim > 2) ? iuub[2] : 0;

            varl = VarList;
            int k = 0;
            while (varl) {
                gpu_find_maximum_launch(
                    BP->stream,
                    shape_0, shape_1, shape_2,
                    pshape_0, pshape_1, pshape_2,
                    BP->d_X[0], BP->d_X[1], BP->d_X[2],
                    BP->d_fgfs[varl->data->sgfn],
                    llb_0, llb_1, llb_2,
                    uub_0, uub_1, uub_2,
                    d_max_val, d_max_pos_x, d_max_pos_y, d_max_pos_z
                );

                // Sync and copy results to host
                GPUManager::synchronize_stream(BP->stream);

                double h_max_val = 0.0;
                double h_max_pos[3] = {0.0, 0.0, 0.0};
                GPUManager::getInstance().sync_to_cpu(&h_max_val, d_max_val, 1, BP->stream);
                GPUManager::getInstance().sync_to_cpu(&h_max_pos[0], d_max_pos_x, 1, BP->stream);
                GPUManager::getInstance().sync_to_cpu(&h_max_pos[1], d_max_pos_y, 1, BP->stream);
                GPUManager::getInstance().sync_to_cpu(&h_max_pos[2], d_max_pos_z, 1, BP->stream);
                GPUManager::synchronize_stream(BP->stream);

                if (h_max_val > shellf[k]) {
                    shellf[k] = h_max_val;
                    for (int i = 0; i < dim; i++)
                        xx[dim * k + i] = h_max_pos[i];
                }

                varl = varl->next;
                k++;
            }
        }

        if (Bp == ble) break;
        Bp = Bp->next;
    }

    // MPI_Allreduce with MPI_DOUBLE_INT MPI_MAXLOC
    struct mloc {
        double val;
        int rank;
    };

    mloc *IN = new mloc[num_var];
    mloc *OUT = new mloc[num_var];
    for (int i = 0; i < num_var; i++) {
        IN[i].val = shellf[i];
        IN[i].rank = myrank;
    }

    MPI_Allreduce(IN, OUT, num_var, MPI_DOUBLE_INT, MPI_MAXLOC, Comm_here);

    // Zero out positions for non-winning ranks
    for (int i = 0; i < num_var; i++) {
        shellf[i] = OUT[i].val;
        if (myrank != OUT[i].rank) {
            for (int d = 0; d < dim; d++)
                xx[dim * i + d] = 0;
        }
    }

    // Sum positions (only the winning rank has non-zero)
    double *global_xx = new double[dim * num_var];
    MPI_Allreduce(xx, global_xx, dim * num_var, MPI_DOUBLE, MPI_SUM, Comm_here);

    // Copy results to device
    GPUManager::getInstance().sync_to_gpu(shellf, d_Shellf, num_var, stream);
    GPUManager::getInstance().sync_to_gpu(global_xx, d_XX, dim * num_var, stream);
    GPUManager::synchronize_stream(stream);

    // Free all memory
    delete[] IN;
    delete[] OUT;
    delete[] shellf;
    delete[] xx;
    delete[] global_xx;
    delete[] DH;

    GPUManager::getInstance().free_device_memory(d_max_val, 1);
    GPUManager::getInstance().free_device_memory(d_max_pos_x, 1);
    GPUManager::getInstance().free_device_memory(d_max_pos_y, 1);
    GPUManager::getInstance().free_device_memory(d_max_pos_z, 1);
}
