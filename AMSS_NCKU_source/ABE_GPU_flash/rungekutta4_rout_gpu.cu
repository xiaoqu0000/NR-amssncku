#include "rungekutta4_rout.h"

#include <cuda_runtime.h>
#include <iostream>
#include <math.h>
#include <vector>

#include "gpu_runtime_coordinator.h"

constexpr int kRungeKuttaBlockSize = 256;
constexpr double kHalfStepWeight = 0.5;
constexpr double kAccumulationFactor = 2.0;
constexpr double kFinalBlendWeight = 1.0 / 6.0;

namespace {

struct RungeKuttaBatchScratch {
    cudaStream_t stream;
    const double** device_state_origins;
    double** device_state_stages;
    double** device_rhs_accumulators;
    int capacity;
};

std::vector<RungeKuttaBatchScratch>& rungekutta_batch_scratches() {
    static std::vector<RungeKuttaBatchScratch> scratches;
    return scratches;
}

RungeKuttaBatchScratch& rungekutta_batch_scratch_for_stream(
    std::vector<RungeKuttaBatchScratch>& scratches,
    cudaStream_t stream
) {
    for (std::vector<RungeKuttaBatchScratch>::iterator scratch = scratches.begin();
         scratch != scratches.end();
         ++scratch) {
        if (scratch->stream == stream) {
            return *scratch;
        }
    }

    RungeKuttaBatchScratch scratch = {stream, 0, 0, 0, 0};
    scratches.push_back(scratch);
    return scratches.back();
}

void ensure_rungekutta_batch_capacity(RungeKuttaBatchScratch& scratch, int field_count) {
    if (field_count <= scratch.capacity) {
        return;
    }

    if (scratch.device_state_origins) {
        CUDA_CHECK(cudaFree(scratch.device_state_origins));
        scratch.device_state_origins = 0;
    }
    if (scratch.device_state_stages) {
        CUDA_CHECK(cudaFree(scratch.device_state_stages));
        scratch.device_state_stages = 0;
    }
    if (scratch.device_rhs_accumulators) {
        CUDA_CHECK(cudaFree(scratch.device_rhs_accumulators));
        scratch.device_rhs_accumulators = 0;
    }

    CUDA_CHECK(cudaMalloc(&scratch.device_state_origins, field_count * sizeof(double*)));
    CUDA_CHECK(cudaMalloc(&scratch.device_state_stages, field_count * sizeof(double*)));
    CUDA_CHECK(cudaMalloc(&scratch.device_rhs_accumulators, field_count * sizeof(double*)));
    scratch.capacity = field_count;
}

}  // namespace

__global__ void advance_rungekutta_stage_kernel(
    int cell_count,
    const double* state_origin, double* state_stage, double* rhs_accumulator, double time_step, int rk_stage
);

__global__ void advance_rungekutta_stage_batch_kernel(
    int cell_count,
    int field_count,
    const double* const* state_origins,
    double* const* state_stages,
    double* const* rhs_accumulators,
    double time_step,
    int rk_stage
);

inline int flatten_extent(const int extent[3]) {
    return extent[0] * extent[1] * extent[2];
}

inline int launch_blocks_for(int item_count) {
    return (item_count + kRungeKuttaBlockSize - 1) / kRungeKuttaBlockSize;
}

void gpu_rungekutta4_rout_launch(
    cudaStream_t &stream,
    int ex[3], double dT,
    double* d_f0, double* d_f1, double* d_f_rhs, int RK4
) {
    const int cell_count = flatten_extent(ex);
    const int grid_span = launch_blocks_for(cell_count);

    advance_rungekutta_stage_kernel<<<grid_span, kRungeKuttaBlockSize, 0, stream>>>(
        cell_count,
        d_f0, d_f1, d_f_rhs, dT, RK4
    );
}

void gpu_rungekutta4_rout_batch_launch(
    cudaStream_t &stream,
    int ex[3], double dT,
    const std::vector<double*>& d_f0_fields,
    const std::vector<double*>& d_f1_fields,
    const std::vector<double*>& d_f_rhs_fields,
    int RK4
) {
    const int field_count = static_cast<int>(d_f0_fields.size());
    if (field_count == 0) {
        return;
    }

    if (d_f1_fields.size() != d_f0_fields.size() || d_f_rhs_fields.size() != d_f0_fields.size()) {
        std::cerr << "gpu_rungekutta4_rout_batch_launch: field list sizes do not match" << std::endl;
        return;
    }

    const int cell_count = flatten_extent(ex);
    const int grid_span = launch_blocks_for(cell_count);
    RungeKuttaBatchScratch& scratch = rungekutta_batch_scratch_for_stream(
        rungekutta_batch_scratches(),
        stream
    );
    ensure_rungekutta_batch_capacity(scratch, field_count);

    CUDA_CHECK(cudaMemcpyAsync(
        scratch.device_state_origins,
        d_f0_fields.data(),
        field_count * sizeof(double*),
        cudaMemcpyHostToDevice,
        stream
    ));
    CUDA_CHECK(cudaMemcpyAsync(
        scratch.device_state_stages,
        d_f1_fields.data(),
        field_count * sizeof(double*),
        cudaMemcpyHostToDevice,
        stream
    ));
    CUDA_CHECK(cudaMemcpyAsync(
        scratch.device_rhs_accumulators,
        d_f_rhs_fields.data(),
        field_count * sizeof(double*),
        cudaMemcpyHostToDevice,
        stream
    ));

    const dim3 grid(grid_span, field_count, 1);
    advance_rungekutta_stage_batch_kernel<<<grid, kRungeKuttaBlockSize, 0, stream>>>(
        cell_count,
        field_count,
        scratch.device_state_origins,
        scratch.device_state_stages,
        scratch.device_rhs_accumulators,
        dT,
        RK4
    );
}

__global__ void advance_rungekutta_stage_kernel(
    int cell_count,
    const double* state_origin, double* state_stage, double* rhs_accumulator, double time_step, int rk_stage
) {
    const int point_id = blockIdx.x * blockDim.x + threadIdx.x;
    if (point_id >= cell_count) {
        return;
    }

    const double base_value = state_origin[point_id];
    const double stage_value = state_stage[point_id];
    const double rhs_value = rhs_accumulator[point_id];

    switch (rk_stage) {
        case 0:
            state_stage[point_id] = base_value + kHalfStepWeight * time_step * rhs_value;
            break;
        case 1:
            rhs_accumulator[point_id] = rhs_value + kAccumulationFactor * stage_value;
            state_stage[point_id] = base_value + kHalfStepWeight * time_step * stage_value;
            break;
        case 2:
            rhs_accumulator[point_id] = rhs_value + kAccumulationFactor * stage_value;
            state_stage[point_id] = base_value + time_step * stage_value;
            break;
        case 3:
            state_stage[point_id] = base_value + kFinalBlendWeight * time_step * (stage_value + rhs_value);
            break;
        default:
            break;
    }
}

__global__ void advance_rungekutta_stage_batch_kernel(
    int cell_count,
    int field_count,
    const double* const* state_origins,
    double* const* state_stages,
    double* const* rhs_accumulators,
    double time_step,
    int rk_stage
) {
    const int field_index = blockIdx.y;
    if (field_index >= field_count) {
        return;
    }

    const int point_id = blockIdx.x * blockDim.x + threadIdx.x;
    if (point_id >= cell_count) {
        return;
    }

    const double* state_origin = state_origins[field_index];
    double* state_stage = state_stages[field_index];
    double* rhs_accumulator = rhs_accumulators[field_index];

    const double base_value = state_origin[point_id];
    const double stage_value = state_stage[point_id];
    const double rhs_value = rhs_accumulator[point_id];

    switch (rk_stage) {
        case 0:
            state_stage[point_id] = base_value + kHalfStepWeight * time_step * rhs_value;
            break;
        case 1:
            rhs_accumulator[point_id] = rhs_value + kAccumulationFactor * stage_value;
            state_stage[point_id] = base_value + kHalfStepWeight * time_step * stage_value;
            break;
        case 2:
            rhs_accumulator[point_id] = rhs_value + kAccumulationFactor * stage_value;
            state_stage[point_id] = base_value + time_step * stage_value;
            break;
        case 3:
            state_stage[point_id] = base_value + kFinalBlendWeight * time_step * (stage_value + rhs_value);
            break;
        default:
            break;
    }
}
