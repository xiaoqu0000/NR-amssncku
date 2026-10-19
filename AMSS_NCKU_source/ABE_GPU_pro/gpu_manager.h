#ifndef GPU_MANAGER_H
#define GPU_MANAGER_H

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
    } while(0)

class GPUManager {
private:
    // 声明一个不透明的结构体指针，具体实现在 .cu 中
    struct Impl;
    Impl* pimpl;

    GPUManager();
    ~GPUManager();

public:
    static GPUManager& getInstance();

    template<typename T> T* allocate_device_memory(size_t num_elements, cudaStream_t stream = 0, bool pooled = true);
    template<typename T> void free_device_memory(T* d_ptr, size_t num_elements = 0, cudaStream_t stream = 0, bool pooled = true);

    // Allocate/free GPU buffers intended to be handed to MPI (CUDA-aware) calls.
    // Backed by synchronous cudaMalloc/cudaFree so the virtual-to-physical mapping is
    // stable for IB memory-region registration; pooled in a separate free-list to amortize.
    template<typename T> T* allocate_comm_buffer(size_t num_elements, cudaStream_t stream = 0);
    template<typename T> void free_comm_buffer(T* d_ptr, size_t num_elements = 0, cudaStream_t stream = 0);

    template<typename T> void memset_device_memory(T* d_ptr, size_t num_elements, cudaStream_t stream = 0);

    template<typename T> void sync_to_gpu(const T* h_ptr, T* d_ptr, size_t num_elements, cudaStream_t stream = 0);
    template<typename T> void sync_to_cpu(T* h_ptr, const T* d_ptr, size_t num_elements, cudaStream_t stream = 0);

    void unpack_to_pad(const double* src, double* dst, int shape[3], int shape_pad[3], cudaStream_t stream = 0);
    void pack_to_pad(const double* src, double* dst, int shape[3], int shape_pad[3], cudaStream_t stream = 0);

    cudaEvent_t event_acquire();
    void event_release(cudaEvent_t e);
    void event_record(cudaEvent_t e, cudaStream_t stream);
    void event_host_wait(cudaEvent_t e);

    cudaStream_t get_stream();
    // Internal memory-stream — used by allocate/memset and as a join point for
    // cross-stream rendezvous in gpu_transfer. Returned to callers so they can
    // post cudaStreamWaitEvent / cudaEventRecord against a non-default stream
    // (avoids the legacy default-stream's implicit cross-stream sync).
    cudaStream_t memory_stream() const;
    void synchronize_all();
    void synchronize_pool();
    void synchronize_memory();
    // Device-side fence: record an event on memory_stream and have every pool
    // stream cudaStreamWaitEvent it. Equivalent to synchronize_memory() in terms
    // of ordering (subsequent kernels on pool streams see all enqueued memory_stream
    // work) but does NOT block the host, so the host can race ahead and queue more
    // work. Use this where the host doesn't actually need to know the work is done.
    void fence_pool_after_memory();
    static void synchronize_stream(cudaStream_t stream);
    void drain_buffer_pool();
    void set_buffer_pool_cap(size_t bytes);
    void drain_comm_pool();
    void set_comm_pool_cap(size_t bytes);
};

extern template double* GPUManager::allocate_device_memory<double>(size_t num_elements, cudaStream_t stream, bool pooled);
extern template int* GPUManager::allocate_device_memory<int>(size_t num_elements, cudaStream_t stream, bool pooled);
extern template void GPUManager::free_device_memory<double>(double* d_ptr, size_t num_elements, cudaStream_t stream, bool pooled);
extern template void GPUManager::free_device_memory<int>(int* d_ptr, size_t num_elements, cudaStream_t stream, bool pooled);
extern template double* GPUManager::allocate_comm_buffer<double>(size_t num_elements, cudaStream_t stream);
extern template int* GPUManager::allocate_comm_buffer<int>(size_t num_elements, cudaStream_t stream);
extern template void GPUManager::free_comm_buffer<double>(double* d_ptr, size_t num_elements, cudaStream_t stream);
extern template void GPUManager::free_comm_buffer<int>(int* d_ptr, size_t num_elements, cudaStream_t stream);
extern template void GPUManager::memset_device_memory<double>(double* d_ptr, size_t num_elements, cudaStream_t stream);
extern template void GPUManager::memset_device_memory<int>(int* d_ptr, size_t num_elements, cudaStream_t stream);
extern template void GPUManager::sync_to_gpu<double>(const double* h_ptr, double* d_ptr, size_t num_elements, cudaStream_t stream);
extern template void GPUManager::sync_to_gpu<int>(const int* h_ptr, int* d_ptr, size_t num_elements, cudaStream_t stream);
extern template void GPUManager::sync_to_cpu<double>(double* h_ptr, const double* d_ptr, size_t num_elements, cudaStream_t stream);
extern template void GPUManager::sync_to_cpu<int>(int* h_ptr, const int* d_ptr, size_t num_elements, cudaStream_t stream);

#endif /* GPU_MANAGER_H */