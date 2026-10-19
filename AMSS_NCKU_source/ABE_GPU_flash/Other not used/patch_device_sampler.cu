#include "MPatch.h"

#include "mpi_cuda_aware_runtime.h"

#include <cuda_runtime.h>
#include "fmisc.h"
#include "macrodef.h"

using namespace std;

namespace {

void synchronize_patch_streams(Patch* patch, int myrank) {
    MyList<Block>* block_node = patch->blb;
    while (block_node) {
        Block* block = block_node->data;
        if (myrank == block->rank) {
            CUDA_CHECK(cudaStreamSynchronize(block->stream));
        }
        if (block_node == patch->ble) break;
        block_node = block_node->next;
    }
}

double interpolation_block_trim(double spacing) {
#ifdef Vertex
    return (static_cast<double>(ghost_width) - 0.5) * spacing;
#else
    return static_cast<double>(ghost_width) * spacing;
#endif
}

}

__global__ void shellf_weight_normalizer_kernel(
    int NN, int num_var, 
    double* d_Shellf, 
    const int* d_Weight, 
    int* d_err_flag
) {
    int point_index = blockIdx.x * blockDim.x + threadIdx.x;
    if (point_index >= NN) return;

    int weight = d_Weight[point_index];
    if (weight > 1) {
        for (int var_index = 0; var_index < num_var; var_index++) {
            d_Shellf[var_index + point_index * num_var] /= (double)weight;
        }
    } else if (weight == 0) {
        *d_err_flag = 1;
    }
}

void launch_shellf_weight_normalizer(
    cudaStream_t stream,
    int NN, int num_var, 
    double* d_Shellf, 
    const int* d_Weight, 
    int* d_err_flag
) {
    int threads_per_block = 256;
    int block_count = (NN + threads_per_block - 1) / threads_per_block;
    shellf_weight_normalizer_kernel<<<block_count, threads_per_block, 0, stream>>>(
        NN, num_var, d_Shellf, d_Weight, d_err_flag
    );
}

bool Patch::Sample_Points_Device_Local(
    MyList<var> *VarList, 
    int NN, double *d_XX_0, double *d_XX_1, double *d_XX_2,
    double *d_shellf, int *d_weight, int Symmetry
) {
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int ordn = 2 * ghost_width;
    int num_var = 0;
    MyList<var> *variable_cursor = VarList;
    while (variable_cursor) {
        num_var++;
        variable_cursor = variable_cursor->next;
    }

    double cell_spacing[3];
    for (int axis = 0; axis < dim; axis++) {
        cell_spacing[axis] = getdX(axis);
    }

    MyList<Block> *block_node = blb;
    while (block_node) {
        Block *block = block_node->data;
        if (myrank == block->rank) {
            double lower_bound[3], upper_bound[3];
            for (int axis = 0; axis < dim; axis++) {
                const double block_trim = interpolation_block_trim(cell_spacing[axis]);
                lower_bound[axis] = (feq(block->bbox[axis], bbox[axis], cell_spacing[axis] / 2))
                    ? block->bbox[axis] + lli[axis] * cell_spacing[axis]
                    : block->bbox[axis] + block_trim;
                upper_bound[axis] = (feq(block->bbox[dim + axis], bbox[dim + axis], cell_spacing[axis] / 2))
                    ? block->bbox[dim + axis] - uui[axis] * cell_spacing[axis]
                    : block->bbox[dim + axis] - block_trim;
            }

            variable_cursor = VarList;
            int variable_slot = 0;
            while (variable_cursor) {
                gpu_global_interp_launch(
                    block->stream, NN, dim,
                    d_XX_0, d_XX_1, d_XX_2,
                    block->shape[0], block->shape[1], block->shape[2],
                    block->d_X[0], block->d_X[1], block->d_X[2],
                    block->d_fgfs[variable_cursor->data->sgfn],
                    lower_bound[0], lower_bound[1], lower_bound[2],
                    upper_bound[0], upper_bound[1], upper_bound[2],
                    ordn,
                    variable_cursor->data->SoA[0],
                    variable_cursor->data->SoA[1],
                    variable_cursor->data->SoA[2],
                    Symmetry, variable_slot, num_var, d_shellf, d_weight
                );
                variable_cursor = variable_cursor->next;
                variable_slot++;
            }
        }
        block_node = block_node->next;
    }

    synchronize_patch_streams(this, myrank);
    
    double *h_shellf_local = new double[NN * num_var];
    double *h_shellf_global = new double[NN * num_var];
    GpuRuntimeCoordinator::instance().copy_to_host(h_shellf_local, d_shellf, NN * num_var);
    GpuRuntimeCoordinator::instance().wait_for_memory_ops();
    MPI_Allreduce(h_shellf_local, h_shellf_global, NN * num_var, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    GpuRuntimeCoordinator::instance().copy_to_device(h_shellf_global, d_shellf, NN * num_var);

    delete[] h_shellf_local;
    delete[] h_shellf_global;
    GpuRuntimeCoordinator::instance().wait_for_memory_ops();
    return true;
}

void Patch::Sample_Points_Device(
    cudaStream_t stream,
    MyList<var> *VarList,
    int NN, double **d_XX,
    double *d_Shellf, int Symmetry
) {
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int ordn = 2 * ghost_width;
    int num_var = 0;
    MyList<var> *variable_cursor = VarList;
    while (variable_cursor) {
        num_var++;
        variable_cursor = variable_cursor->next;
    }

    double *cell_spacing = new double[dim];
    for (int axis = 0; axis < dim; axis++) {
        cell_spacing[axis] = getdX(axis);
    }
    double *lower_bound = new double[dim];
    double *upper_bound = new double[dim];

    double *d_local_shellf = GpuRuntimeCoordinator::instance().allocate_device_buffer<double>(NN * num_var);
    int *d_local_weight = GpuRuntimeCoordinator::instance().allocate_device_buffer<int>(NN);
    GpuRuntimeCoordinator::instance().wait_for_memory_ops();

    MyList<Block> *block_node = blb;
    while (block_node) {
        Block *block = block_node->data;
        if (myrank == block->rank) {
            for (int axis = 0; axis < dim; axis++) {
                const double block_trim = interpolation_block_trim(cell_spacing[axis]);
                lower_bound[axis] = (feq(block->bbox[axis], bbox[axis], cell_spacing[axis] / 2))
                    ? block->bbox[axis] + lli[axis] * cell_spacing[axis]
                    : block->bbox[axis] + block_trim;
                upper_bound[axis] = (feq(block->bbox[dim + axis], bbox[dim + axis], cell_spacing[axis] / 2))
                    ? block->bbox[dim + axis] - uui[axis] * cell_spacing[axis]
                    : block->bbox[dim + axis] - block_trim;
            }

            int shape_x = block->shape[0];
            int shape_y = (dim > 1) ? block->shape[1] : 1;
            int shape_z = (dim > 2) ? block->shape[2] : 1;

            double lower_x = lower_bound[0];
            double lower_y = (dim > 1) ? lower_bound[1] : 0.0;
            double lower_z = (dim > 2) ? lower_bound[2] : 0.0;
            double upper_x = upper_bound[0];
            double upper_y = (dim > 1) ? upper_bound[1] : 0.0;
            double upper_z = (dim > 2) ? upper_bound[2] : 0.0;

            variable_cursor = VarList;
            int variable_slot = 0;
            while (variable_cursor) {
                gpu_global_interp_launch(
                    block->stream,
                    NN, dim,
                    d_XX[0], d_XX[1], d_XX[2],
                    shape_x, shape_y, shape_z,
                    block->d_X[0], block->d_X[1], block->d_X[2],
                    block->d_fgfs[variable_cursor->data->sgfn],
                    lower_x, lower_y, lower_z,
                    upper_x, upper_y, upper_z,
                    ordn,
                    variable_cursor->data->SoA[0],
                    variable_cursor->data->SoA[1],
                    variable_cursor->data->SoA[2],
                    Symmetry, variable_slot, num_var, d_local_shellf, d_local_weight
                );
                variable_cursor = variable_cursor->next;
                variable_slot++;
            }
        }
        if (block_node == ble) break;
        block_node = block_node->next;
    }

    synchronize_patch_streams(this, myrank);

    double *h_local_shellf = new double[NN * num_var];
    int *h_local_weight = new int[NN];
    int *h_global_weight = h_local_weight;
    int *d_global_weight = GpuRuntimeCoordinator::instance().allocate_device_buffer<int>(NN, stream, false);

    if (amss_mpi_cuda_aware_runtime_enabled()) {
        MPI_Allreduce(d_local_shellf, d_Shellf, NN * num_var, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(d_local_weight, d_global_weight, NN, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        GpuRuntimeCoordinator::instance().copy_to_host(h_global_weight, d_global_weight, NN, stream);
    }
    else {
        GpuRuntimeCoordinator::instance().copy_to_host(h_local_shellf, d_local_shellf, NN * num_var, stream);
        GpuRuntimeCoordinator::instance().copy_to_host(h_local_weight, d_local_weight, NN, stream);

        GpuRuntimeCoordinator::wait_for_stream(stream);

        MPI_Allreduce(MPI_IN_PLACE, h_local_shellf, NN * num_var, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, h_local_weight, NN, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        GpuRuntimeCoordinator::instance().copy_to_device(h_local_shellf, d_Shellf, NN * num_var, stream);
        GpuRuntimeCoordinator::instance().copy_to_device(h_local_weight, d_global_weight, NN, stream);
    }

    int *d_err_flag = GpuRuntimeCoordinator::instance().allocate_device_buffer<int>(1, stream);

    launch_shellf_weight_normalizer(stream, NN, num_var, d_Shellf, d_global_weight, d_err_flag);

    int h_err_flag = 0;
    GpuRuntimeCoordinator::instance().copy_to_host(&h_err_flag, d_err_flag, 1, stream);

    GpuRuntimeCoordinator::wait_for_stream(stream);

    if (h_err_flag > 0) {
        if (myrank == 0) {
            double *h_coords = new double[dim * NN];
            double *h_XX[3] = {0, 0, 0};
            for (int axis = 0; axis < dim; axis++) {
                h_XX[axis] = h_coords + axis * NN;
                GpuRuntimeCoordinator::instance().copy_to_host(h_XX[axis], d_XX[axis], NN, stream);
            }
            GpuRuntimeCoordinator::wait_for_stream(stream);

            for (int point_index = 0; point_index < NN; point_index++) {
                if (h_global_weight[point_index] > 1) {
                    cout << "WARNING: Patch::Interp_Points meets multiple weight at point " << point_index << endl;
                } else if (h_global_weight[point_index] == 0) {
                    cout << "ERROR: Patch::Interp_Points fails to find point (";
                    for (int axis = 0; axis < dim; axis++) {
                        cout << h_XX[axis][point_index];
                        if (axis < dim - 1) cout << ",";
                        else cout << ") ! Out of domain bounds." << endl;
                    }
                }
            }
            delete[] h_coords;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    delete[] cell_spacing;
    delete[] lower_bound;
    delete[] upper_bound;
    delete[] h_local_shellf;
    delete[] h_local_weight;
    GpuRuntimeCoordinator::instance().release_device_buffer(d_local_shellf, NN * num_var, stream);
    GpuRuntimeCoordinator::instance().release_device_buffer(d_local_weight, NN, stream);
    GpuRuntimeCoordinator::instance().release_device_buffer(d_global_weight, NN, stream);
    GpuRuntimeCoordinator::instance().release_device_buffer(d_err_flag, 1, stream);
    GpuRuntimeCoordinator::wait_for_stream(stream);
}
