#ifndef GPU_RUNTIME_COORDINATOR_H
#define GPU_RUNTIME_COORDINATOR_H

#include <iostream>
#include <cuda_runtime.h>

#define CUDA_CHECK(call) \
    do { \
        cudaError_t err = (call); \
        if (err != cudaSuccess) { \
            std::cerr << "CUDA Error: " << cudaGetErrorString(err) \
                      << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
            exit(EXIT_FAILURE); \
        } \
    } while (0)

class GpuRuntimeCoordinator {
private:
    struct RuntimeFacility;
    RuntimeFacility* facility_;

    GpuRuntimeCoordinator();
    ~GpuRuntimeCoordinator();

public:
    template<typename T> T* allocate_device_buffer(size_t element_count, cudaStream_t stream = 0, bool zero_initialize = true);
    template<typename T> void release_device_buffer(T* device_ptr, size_t element_count = 0, cudaStream_t stream = 0);

    template<typename T> void copy_to_device(const T* host_ptr, T* device_ptr, size_t element_count, cudaStream_t stream = 0);
    template<typename T> void copy_to_host(T* host_ptr, const T* device_ptr, size_t element_count, cudaStream_t stream = 0);

    static GpuRuntimeCoordinator& instance();
    cudaStream_t acquire_stream();
    void wait_for_memory_ops();
    void wait_for_device();
    static void wait_for_stream(cudaStream_t stream);
    void clear_pool();
};

extern template double* GpuRuntimeCoordinator::allocate_device_buffer<double>(size_t element_count, cudaStream_t stream, bool zero_initialize);
extern template int* GpuRuntimeCoordinator::allocate_device_buffer<int>(size_t element_count, cudaStream_t stream, bool zero_initialize);
extern template void GpuRuntimeCoordinator::release_device_buffer<double>(double* device_ptr, size_t element_count, cudaStream_t stream);
extern template void GpuRuntimeCoordinator::release_device_buffer<int>(int* device_ptr, size_t element_count, cudaStream_t stream);
extern template void GpuRuntimeCoordinator::copy_to_device<double>(const double* host_ptr, double* device_ptr, size_t element_count, cudaStream_t stream);
extern template void GpuRuntimeCoordinator::copy_to_device<int>(const int* host_ptr, int* device_ptr, size_t element_count, cudaStream_t stream);
extern template void GpuRuntimeCoordinator::copy_to_host<double>(double* host_ptr, const double* device_ptr, size_t element_count, cudaStream_t stream);
extern template void GpuRuntimeCoordinator::copy_to_host<int>(int* host_ptr, const int* device_ptr, size_t element_count, cudaStream_t stream);

#endif /* GPU_RUNTIME_COORDINATOR_H */
