#ifndef SHELL_PATCH_GPU_H
#define SHELL_PATCH_GPU_H

#include <cuda_runtime.h>

#include "ShellPatch.h"

void shell_patch_intertransfer_gpu_launch(
    cudaStream_t stream,
    const ShellPatchDeviceTransferList& transfer,
    double* const* d_source_fields,
    double* const* d_target_fields,
    int source_block_count,
    int target_block_count,
    double SoA_0,
    double SoA_1,
    double SoA_2,
    int Symmetry,
    int* d_error_flag);

#endif /* SHELL_PATCH_GPU_H */
