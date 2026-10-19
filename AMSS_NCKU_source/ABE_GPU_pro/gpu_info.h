#ifndef GPU_INFO_H
#define GPU_INFO_H

#define GPU_DIM_X 8
#define GPU_DIM_Y 8
#define GPU_DIM_Z 4

constexpr int GPU_ALIGNMENT_ASYNC_B = 128; // 128 bytes alignment for best performance of cuda::memcpy_async
constexpr int SIZE_DOUBLE_B = sizeof(double);
constexpr int GPU_ALIGNMENT_ASYNC_DOUBLE = GPU_ALIGNMENT_ASYNC_B / SIZE_DOUBLE_B; // 16 doubles alignment

constexpr int GPU_ALIGNMENT_B = 32; // 32 bytes alignment
constexpr int GPU_ALIGNMENT_DOUBLE = GPU_ALIGNMENT_B / SIZE_DOUBLE_B; // 4 doubles alignment

// x-direction pad offset: ghost cells aligned to GPU_ALIGNMENT_DOUBLE for efficient memory access
constexpr int GPU_PAD_X_LEFT = (ghost_width + GPU_ALIGNMENT_DOUBLE - 1) / GPU_ALIGNMENT_DOUBLE * GPU_ALIGNMENT_DOUBLE;

#endif   /* GPU_INFO_H */