#include "gpu_runtime_coordinator.h"
#include "ShellPatch.h"

#include <atomic>
#include <cstdint>
#include <cstdlib>

namespace {

cudaStream_t choose_operation_lane(cudaStream_t preferred_lane, cudaStream_t fallback_lane) {
    if (preferred_lane != 0) {
        return preferred_lane;
    }

    return fallback_lane;
}

unsigned int reserve_lane_index(std::atomic<unsigned int>& lane_cursor) {
    return lane_cursor.fetch_add(1);
}

void report_allocation_failure(
    cudaError_t allocation_error,
    size_t element_count,
    size_t element_size,
    cudaStream_t stream,
    const char* file,
    int line
) {
    size_t free_bytes = 0;
    size_t total_bytes = 0;
    cudaError_t mem_info_error = cudaMemGetInfo(&free_bytes, &total_bytes);

    std::cerr << "CUDA Error: " << cudaGetErrorString(allocation_error)
              << " at " << file << ":" << line
              << " while allocating device buffer"
              << " elements=" << element_count
              << " element_size=" << element_size
              << " bytes=" << (element_count * element_size)
              << " stream=" << stream;
    if (mem_info_error == cudaSuccess) {
        std::cerr << " cuda_free=" << free_bytes
                  << " cuda_total=" << total_bytes;
    } else {
        std::cerr << " cudaMemGetInfo_error=" << cudaGetErrorString(mem_info_error);
    }
    std::cerr << std::endl;
    exit(EXIT_FAILURE);
}

}  // namespace

struct GpuRuntimeCoordinator::RuntimeFacility {
    static constexpr int kWorkerLaneCount = 16;
    cudaStream_t work_lanes[kWorkerLaneCount];
    std::atomic<unsigned int> lane_cursor{0};
    cudaStream_t staging_lane;
};

cudaStream_t GpuRuntimeCoordinator::acquire_stream() {
    const unsigned int lane_index = reserve_lane_index(facility_->lane_cursor);
    return facility_->work_lanes[lane_index % RuntimeFacility::kWorkerLaneCount];
}

GpuRuntimeCoordinator& GpuRuntimeCoordinator::instance() {
    static GpuRuntimeCoordinator runtime_coordinator;
    return runtime_coordinator;
}

void GpuRuntimeCoordinator::wait_for_memory_ops() {
    CUDA_CHECK(cudaStreamSynchronize(facility_->staging_lane));
}

GpuRuntimeCoordinator::GpuRuntimeCoordinator() : facility_(new RuntimeFacility()) {
    for (int lane_index = 0; lane_index < RuntimeFacility::kWorkerLaneCount; ++lane_index) {
        CUDA_CHECK(cudaStreamCreate(&facility_->work_lanes[lane_index]));
    }
    CUDA_CHECK(cudaStreamCreate(&facility_->staging_lane));
}

template<typename T>
T* GpuRuntimeCoordinator::allocate_device_buffer(size_t element_count, cudaStream_t stream, bool zero_initialize) {
    T* device_ptr = nullptr;
    const cudaStream_t allocation_lane = choose_operation_lane(stream, facility_->staging_lane);
    if (element_count > SIZE_MAX / sizeof(T)) {
        std::cerr << "CUDA allocation size overflow at " << __FILE__ << ":" << __LINE__
                  << " elements=" << element_count
                  << " element_size=" << sizeof(T)
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    cudaError_t allocation_error = cudaMallocAsync((void**)&device_ptr, element_count * sizeof(T), allocation_lane);
    if (allocation_error != cudaSuccess) {
        report_allocation_failure(allocation_error, element_count, sizeof(T), allocation_lane, __FILE__, __LINE__);
    }
    if (zero_initialize && element_count > 0) {
        CUDA_CHECK(cudaMemsetAsync(device_ptr, 0, element_count * sizeof(T), allocation_lane));
    }
    return device_ptr;
}

template<typename T>
void GpuRuntimeCoordinator::copy_to_host(T* host_ptr, const T* device_ptr, size_t element_count, cudaStream_t stream) {
    const cudaStream_t transfer_lane = choose_operation_lane(stream, instance().facility_->staging_lane);
    CUDA_CHECK(cudaMemcpyAsync(host_ptr, device_ptr, element_count * sizeof(T), cudaMemcpyDeviceToHost, transfer_lane));
}

void GpuRuntimeCoordinator::wait_for_stream(cudaStream_t stream) {
    CUDA_CHECK(cudaStreamSynchronize(stream));
}

GpuRuntimeCoordinator::~GpuRuntimeCoordinator() {
    for (int lane_index = 0; lane_index < RuntimeFacility::kWorkerLaneCount; ++lane_index) {
        cudaStreamDestroy(facility_->work_lanes[lane_index]);
    }
    cudaStreamDestroy(facility_->staging_lane);
    delete facility_;
}

template<typename T>
void GpuRuntimeCoordinator::copy_to_device(const T* host_ptr, T* device_ptr, size_t element_count, cudaStream_t stream) {
    const cudaStream_t transfer_lane = choose_operation_lane(stream, instance().facility_->staging_lane);
    CUDA_CHECK(cudaMemcpyAsync(device_ptr, host_ptr, element_count * sizeof(T), cudaMemcpyHostToDevice, transfer_lane));
}

void GpuRuntimeCoordinator::wait_for_device() {
    CUDA_CHECK(cudaDeviceSynchronize());
}

template<typename T>
void GpuRuntimeCoordinator::release_device_buffer(T* device_ptr, size_t element_count, cudaStream_t stream) {
    (void)element_count;
    if (device_ptr == nullptr) {
        return;
    }

    const cudaStream_t release_lane = choose_operation_lane(stream, facility_->staging_lane);
    CUDA_CHECK(cudaFreeAsync(device_ptr, release_lane));
}

template double* GpuRuntimeCoordinator::allocate_device_buffer<double>(size_t element_count, cudaStream_t stream, bool zero_initialize);
template int* GpuRuntimeCoordinator::allocate_device_buffer<int>(size_t element_count, cudaStream_t stream, bool zero_initialize);
template void GpuRuntimeCoordinator::release_device_buffer<double>(double* device_ptr, size_t element_count, cudaStream_t stream);
template void GpuRuntimeCoordinator::release_device_buffer<int>(int* device_ptr, size_t element_count, cudaStream_t stream);
template void GpuRuntimeCoordinator::copy_to_device<double>(const double* host_ptr, double* device_ptr, size_t element_count, cudaStream_t stream);
template void GpuRuntimeCoordinator::copy_to_device<int>(const int* host_ptr, int* device_ptr, size_t element_count, cudaStream_t stream);
template void GpuRuntimeCoordinator::copy_to_host<double>(double* host_ptr, const double* device_ptr, size_t element_count, cudaStream_t stream);
template void GpuRuntimeCoordinator::copy_to_host<int>(int* host_ptr, const int* device_ptr, size_t element_count, cudaStream_t stream);
template ShellPatchDevicePoint* GpuRuntimeCoordinator::allocate_device_buffer<ShellPatchDevicePoint>(size_t element_count, cudaStream_t stream, bool zero_initialize);
template void GpuRuntimeCoordinator::release_device_buffer<ShellPatchDevicePoint>(ShellPatchDevicePoint* device_ptr, size_t element_count, cudaStream_t stream);
template void GpuRuntimeCoordinator::copy_to_device<ShellPatchDevicePoint>(const ShellPatchDevicePoint* host_ptr, ShellPatchDevicePoint* device_ptr, size_t element_count, cudaStream_t stream);
template void GpuRuntimeCoordinator::copy_to_host<ShellPatchDevicePoint>(ShellPatchDevicePoint* host_ptr, const ShellPatchDevicePoint* device_ptr, size_t element_count, cudaStream_t stream);
