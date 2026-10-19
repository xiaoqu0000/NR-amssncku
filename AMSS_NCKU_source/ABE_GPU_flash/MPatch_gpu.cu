#include "MPatch.h"

#include "mpi_cuda_aware_runtime.h"

#include <cuda_runtime.h>
#include "fmisc.h"
#include "macrodef.h"

using namespace std;

namespace {

constexpr int kNormalizeBlockSize = 256;

void synchronize_patch_streams(Patch* patch, int mpi_rank) {
    for (MyList<Block>* block_node = patch->blb; block_node; block_node = block_node->next) {
        Block* block = block_node->data;
        if (mpi_rank == block->rank) {
            CUDA_CHECK(cudaStreamSynchronize(block->stream));
        }
        if (block_node == patch->ble) {
            break;
        }
    }
}

struct BlockInterpolationRegion {
    int extent[3];
    double lower[3];
    double upper[3];
};

__global__ void rebalance_shellf_kernel(
    int point_count,
    int variable_count,
    double* shell_values,
    const int* hit_count,
    int* error_flag
) {
    const int point_index = blockIdx.x * blockDim.x + threadIdx.x;
    if (point_index >= point_count) {
        return;
    }

    const int local_hits = hit_count[point_index];
    if (local_hits > 1) {
        for (int variable_index = 0; variable_index < variable_count; ++variable_index) {
            shell_values[variable_index + point_index * variable_count] /= static_cast<double>(local_hits);
        }
        return;
    }

    if (local_hits == 0) {
        *error_flag = 1;
    }
}

int launch_blocks_for(int item_count) {
    return (item_count + kNormalizeBlockSize - 1) / kNormalizeBlockSize;
}

int count_interp_variables(MyList<var>* variable_list) {
    int variable_count = 0;
    for (MyList<var>* variable_node = variable_list; variable_node; variable_node = variable_node->next) {
        ++variable_count;
    }
    return variable_count;
}

void fill_patch_spacing(Patch* patch, double spacing[3]) {
    for (int axis = 0; axis < 3; ++axis) {
        spacing[axis] = (axis < dim) ? patch->getdX(axis) : 0.0;
    }
}

BlockInterpolationRegion build_interp_bounds(Patch* patch, Block* block, const double spacing[3]) {
    BlockInterpolationRegion bounds{};
    bounds.extent[0] = block->shape[0];
    bounds.extent[1] = (dim > 1) ? block->shape[1] : 1;
    bounds.extent[2] = (dim > 2) ? block->shape[2] : 1;

    for (int axis = 0; axis < dim; ++axis) {
        double block_trim = ghost_width * spacing[axis];
#ifdef Vertex
        block_trim = (static_cast<double>(ghost_width) - 0.5) * spacing[axis];
#endif
        bounds.lower[axis] = feq(block->bbox[axis], patch->bbox[axis], spacing[axis] / 2.0)
            ? block->bbox[axis] + patch->lli[axis] * spacing[axis]
            : block->bbox[axis] + block_trim;
        bounds.upper[axis] = feq(block->bbox[dim + axis], patch->bbox[dim + axis], spacing[axis] / 2.0)
            ? block->bbox[dim + axis] - patch->uui[axis] * spacing[axis]
            : block->bbox[dim + axis] - block_trim;
    }

    return bounds;
}

void interpolate_block_points(
    Block* block,
    MyList<var>* variable_list,
    int point_count,
    int dimension_count,
    double** device_points,
    const BlockInterpolationRegion& bounds,
    int interpolation_order,
    int symmetry_mode,
    int variable_count,
    double* device_shellf,
    int* device_weight
) {
    const int extent_x = bounds.extent[0];
    const int extent_y = bounds.extent[1];
    const int extent_z = bounds.extent[2];
    const double lower_x = bounds.lower[0];
    const double lower_y = (dimension_count > 1) ? bounds.lower[1] : 0.0;
    const double lower_z = (dimension_count > 2) ? bounds.lower[2] : 0.0;
    const double upper_x = bounds.upper[0];
    const double upper_y = (dimension_count > 1) ? bounds.upper[1] : 0.0;
    const double upper_z = (dimension_count > 2) ? bounds.upper[2] : 0.0;

    int variable_slot = 0;
    for (MyList<var>* variable_node = variable_list; variable_node; variable_node = variable_node->next, ++variable_slot) {
        gpu_global_interp_launch(
            block->stream,
            point_count,
            dimension_count,
            device_points[0],
            device_points[1],
            device_points[2],
            extent_x,
            extent_y,
            extent_z,
            block->d_X[0],
            block->d_X[1],
            block->d_X[2],
            block->d_fgfs[variable_node->data->sgfn],
            lower_x,
            lower_y,
            lower_z,
            upper_x,
            upper_y,
            upper_z,
            interpolation_order,
            variable_node->data->SoA[0],
            variable_node->data->SoA[1],
            variable_node->data->SoA[2],
            symmetry_mode,
            variable_slot,
            variable_count,
            device_shellf,
            device_weight
        );
    }
}

void interpolate_block_coordinates(
    Block* block,
    MyList<var>* variable_list,
    int point_count,
    int dimension_count,
    double* device_x,
    double* device_y,
    double* device_z,
    const BlockInterpolationRegion& bounds,
    int interpolation_order,
    int symmetry_mode,
    int variable_count,
    double* device_shellf,
    int* device_weight
) {
    double* device_points[3] = {device_x, device_y, device_z};
    interpolate_block_points(
        block,
        variable_list,
        point_count,
        dimension_count,
        device_points,
        bounds,
        interpolation_order,
        symmetry_mode,
        variable_count,
        device_shellf,
        device_weight
    );
}

int* reduce_interp_results(
    Patch* patch,
    cudaStream_t stream,
    int point_count,
    int variable_count,
    double* device_local_shellf,
    int* device_local_weight,
    double* device_global_shellf
) {
    int myrank = 0;
    int cpusize = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &cpusize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int* device_global_weight = device_local_weight;

    if (cpusize == 1) {
        if (device_global_shellf != device_local_shellf) {
            CUDA_CHECK(cudaMemcpyAsync(
                device_global_shellf,
                device_local_shellf,
                point_count * variable_count * sizeof(double),
                cudaMemcpyDeviceToDevice,
                stream
            ));
        }
    }
    else if (amss_mpi_cuda_aware_runtime_enabled()) {
        MPI_Allreduce(device_local_shellf, device_global_shellf, point_count * variable_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(device_local_weight, patch->d_interp_global_weight, point_count, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        device_global_weight = patch->d_interp_global_weight;
    }
    else {
        GpuRuntimeCoordinator::instance().copy_to_host(patch->h_interp_shellf, device_local_shellf, point_count * variable_count, stream);
        GpuRuntimeCoordinator::instance().copy_to_host(patch->h_interp_weight, device_local_weight, point_count, stream);
        GpuRuntimeCoordinator::wait_for_stream(stream);

        MPI_Allreduce(MPI_IN_PLACE, patch->h_interp_shellf, point_count * variable_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, patch->h_interp_weight, point_count, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        GpuRuntimeCoordinator::instance().copy_to_device(patch->h_interp_shellf, device_global_shellf, point_count * variable_count, stream);
        GpuRuntimeCoordinator::instance().copy_to_device(patch->h_interp_weight, patch->d_interp_global_weight, point_count, stream);
        device_global_weight = patch->d_interp_global_weight;
    }

    return device_global_weight;
}

void print_interp_failures(
    int mpi_rank,
    int point_count,
    int dimension_count,
    double** device_points,
    cudaStream_t stream,
    const int* device_global_weight
) {
    if (mpi_rank != 0) {
        return;
    }

    int* host_global_weight = new int[point_count];
    GpuRuntimeCoordinator::instance().copy_to_host(host_global_weight, device_global_weight, point_count, stream);

    double* host_points[3] = {nullptr, nullptr, nullptr};
    for (int axis = 0; axis < dimension_count; ++axis) {
        host_points[axis] = new double[point_count];
        GpuRuntimeCoordinator::instance().copy_to_host(host_points[axis], device_points[axis], point_count, stream);
    }
    GpuRuntimeCoordinator::wait_for_stream(stream);

    for (int point_index = 0; point_index < point_count; ++point_index) {
        if (host_global_weight[point_index] > 1) {
            cout << "WARNING: Patch::Interp_Points meets multiple weight at point " << point_index << endl;
        } else if (host_global_weight[point_index] == 0) {
            cout << "ERROR: Patch::Interp_Points fails to find point (";
            for (int axis = 0; axis < dimension_count; ++axis) {
                cout << host_points[axis][point_index];
                if (axis + 1 < dimension_count) {
                    cout << ",";
                } else {
                    cout << ") ! Out of domain bounds." << endl;
                }
            }
        }
    }

    for (int axis = 0; axis < dimension_count; ++axis) {
        delete[] host_points[axis];
    }
    delete[] host_global_weight;
}

}  // namespace

void gpu_normalize_shellf_launch(
    cudaStream_t stream,
    int NN, int num_var,
    double* d_Shellf,
    const int* d_Weight,
    int* d_err_flag
) {
    rebalance_shellf_kernel<<<launch_blocks_for(NN), kNormalizeBlockSize, 0, stream>>>(
        NN,
        num_var,
        d_Shellf,
        d_Weight,
        d_err_flag
    );
}

void Patch::Interp_Points_GPU(
    cudaStream_t stream,
    MyList<var> *VarList,
    int NN, double **d_XX,
    double *d_Shellf, int Symmetry
) {
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    const int interpolation_order = 2 * ghost_width;
    const int variable_count = count_interp_variables(VarList);

    double spacing[3] = {0.0, 0.0, 0.0};
    fill_patch_spacing(this, spacing);

    ensure_interp_device_cache(NN, variable_count, stream);
    ensure_interp_host_cache(NN, variable_count);

    double *d_local_shellf = d_interp_local_shellf;
    int *d_local_weight = d_interp_local_weight;

    CUDA_CHECK(cudaMemsetAsync(d_local_shellf, 0, NN * variable_count * sizeof(double), stream));
    CUDA_CHECK(cudaMemsetAsync(d_local_weight, 0, NN * sizeof(int), stream));
    CUDA_CHECK(cudaMemsetAsync(d_interp_err_flag, 0, sizeof(int), stream));
    GpuRuntimeCoordinator::wait_for_stream(stream);

    for (MyList<Block> *block_node = blb; block_node; block_node = block_node->next) {
        Block *block = block_node->data;
        if (myrank == block->rank) {
            const BlockInterpolationRegion bounds = build_interp_bounds(this, block, spacing);
            interpolate_block_points(
                block,
                VarList,
                NN,
                dim,
                d_XX,
                bounds,
                interpolation_order,
                Symmetry,
                variable_count,
                d_local_shellf,
                d_local_weight
            );
        }
        if (block_node == ble) {
            break;
        }
    }

    synchronize_patch_streams(this, myrank);

    int *d_global_weight = reduce_interp_results(
        this,
        stream,
        NN,
        variable_count,
        d_local_shellf,
        d_local_weight,
        d_Shellf
    );

    int *d_err_flag = d_interp_err_flag;
    gpu_normalize_shellf_launch(stream, NN, variable_count, d_Shellf, d_global_weight, d_err_flag);

    int h_err_flag = 0;
    GpuRuntimeCoordinator::instance().copy_to_host(&h_err_flag, d_err_flag, 1, stream);
    GpuRuntimeCoordinator::wait_for_stream(stream);

    if (h_err_flag > 0) {
        print_interp_failures(myrank, NN, dim, d_XX, stream, d_global_weight);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    GpuRuntimeCoordinator::wait_for_stream(stream);
}

bool Patch::Interp_N_Points_GPU(
    MyList<var> *VarList,
    int NN, double *d_XX_0, double *d_XX_1, double *d_XX_2,
    double *d_shellf, int *d_weight, int Symmetry
) {
    int myrank;
    int cpusize = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &cpusize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    const int interpolation_order = 2 * ghost_width;
    const int variable_count = count_interp_variables(VarList);

    double spacing[3] = {0.0, 0.0, 0.0};
    fill_patch_spacing(this, spacing);

    for (MyList<Block> *block_node = blb; block_node; block_node = block_node->next) {
        Block *block = block_node->data;
        if (myrank == block->rank) {
            const BlockInterpolationRegion bounds = build_interp_bounds(this, block, spacing);
            interpolate_block_coordinates(
                block,
                VarList,
                NN,
                dim,
                d_XX_0,
                d_XX_1,
                d_XX_2,
                bounds,
                interpolation_order,
                Symmetry,
                variable_count,
                d_shellf,
                d_weight
            );
        }
    }

    ensure_interp_host_cache(NN, variable_count);
    synchronize_patch_streams(this, myrank);

    GpuRuntimeCoordinator::instance().copy_to_host(h_interp_shellf, d_shellf, NN * variable_count);
    GpuRuntimeCoordinator::instance().wait_for_memory_ops();
    if (cpusize > 1) {
        MPI_Allreduce(MPI_IN_PLACE, h_interp_shellf, NN * variable_count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    GpuRuntimeCoordinator::instance().copy_to_device(h_interp_shellf, d_shellf, NN * variable_count);

    GpuRuntimeCoordinator::instance().wait_for_memory_ops();
    return true;
}
