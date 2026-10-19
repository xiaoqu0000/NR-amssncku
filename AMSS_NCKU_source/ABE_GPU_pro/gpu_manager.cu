#include "gpu_manager.h"
#include "gpu_common.cuh"
#include <vector>
#include <unordered_map>
#include <mutex>
#include <atomic>
#include <cstdlib>

#include <cassert>

struct GPUManager::Impl {
    static constexpr int NUM_STREAMS = 16;
    cudaStream_t stream_pool[NUM_STREAMS];
    std::atomic<unsigned int> stream_idx{0};
    cudaStream_t memory_stream;

    // --- Buffer Pool ---
    struct PoolEntry {
        void* ptr;
        size_t capacity;    // bytes (size-class aligned)
        PoolEntry* next;
        // For comm pool: event recorded on releaser's stream so the next
        // acquirer can stream-wait on it instead of busy-syncing. nullptr
        // means the releaser didn't pass a stream (caller's responsibility
        // to sync) — in that case no wait is needed.
        cudaEvent_t release_event = nullptr;
        // For comm pool: the stream the buffer was last used on. Acquire
        // prefers same-stream entries (same-stream queue is auto-serialized
        // so no cudaStreamWaitEvent is needed).
        cudaStream_t last_stream = nullptr;
    };

    // Free lists keyed by rounded size class
    std::unordered_map<size_t, PoolEntry*> free_lists;

    // In-use tracking: ptr -> capacity in bytes
    std::unordered_map<void*, size_t> in_use;

    // Total bytes cached in free lists
    size_t total_pooled_bytes = 0;

    // Pool memory cap (default 4GB)
    size_t pool_cap_bytes = 4ULL * 1024 * 1024 * 1024;

    // --- Communication Buffer Pool (cudaMalloc-backed for stable IB MR) ---
    // Stream-ordered cudaMallocAsync would remap physical pages and invalidate
    // IB memory-region cache entries; comm buffers use plain cudaMalloc so the
    // VA->PA mapping is stable until cudaFree.
    std::unordered_map<size_t, PoolEntry*> comm_free_lists;
    std::unordered_map<void*, size_t> comm_in_use;
    size_t comm_total_pooled_bytes = 0;
    size_t comm_pool_cap_bytes = 1ULL * 1024 * 1024 * 1024;  // 1 GiB
    bool comm_pool_cap_warned = false;

    static size_t size_class(size_t bytes) {
        if (bytes == 0) return 128;
        if (bytes <= 4096) {
            return ((bytes + 127) / 128) * 128;
        }
        return ((bytes + 4095) / 4096) * 4096;
    }

    void* pool_acquire(size_t bytes, cudaStream_t stream) {
        size_t sc = size_class(bytes);
        // CUDA's stream-ordered allocator handles cross-stream VA reuse with
        // implicit dependencies: a cudaMallocAsync that returns a VA freed
        // earlier on a different stream auto-inserts a wait on the freeing
        // stream's tail. Drop our custom free-list (which had a race because
        // it bypassed those implicit deps) and trust the driver.
        void* ptr = nullptr;
        CUDA_CHECK(cudaMallocAsync(&ptr, sc, stream));
        CUDA_CHECK(cudaMemsetAsync(ptr, 0, sc, stream));
        in_use[ptr] = sc;
        return ptr;
    }

    void pool_release(void* ptr, cudaStream_t stream) {
        if (ptr == nullptr) return;
        auto it = in_use.find(ptr);
        if (it == in_use.end()) return;
        in_use.erase(it);
        // Caller's `stream` must be the last stream that queued work on ptr;
        // cudaFreeAsync defers the actual free to that stream's tail so any
        // subsequent cudaMallocAsync returning the same VA picks up the
        // implicit dependency. free_device_memory remaps stream==0 to
        // memory_stream upstream of this call.
        CUDA_CHECK(cudaFreeAsync(ptr, stream));
    }

    // Pool drain/drain_all are kept as no-ops: the custom free-list is gone
    // and cudaMallocAsync owns the underlying memory pool. The setters/
    // env-knob (AMSS_POOL_CAP_MB) are preserved as configuration hints in
    // case we ever add an internal cache again.
    void pool_drain() {}

    void pool_drain_all() {
        free_lists.clear();
        total_pooled_bytes = 0;
    }

    // --- Comm pool methods ---
    // Comm buffers stay on plain cudaMalloc/cudaFree (not stream-ordered)
    // to keep the VA->PA mapping stable for IB MR cache. Recycling via
    // free_lists is correct as long as cross-stream reuses honor the
    // releaser's pending ops via per-entry event sync.
    void* comm_pool_acquire(size_t bytes, cudaStream_t stream) {
        size_t sc = size_class(bytes);
        auto it = comm_free_lists.find(sc);
        if (it != comm_free_lists.end() && it->second != nullptr) {
            // First pass: same-stream match — no waitEvent needed (the
            // releaser's tail and the new use are already serialized on
            // this stream).
            PoolEntry** pp = &it->second;
            PoolEntry* picked = nullptr;
            while (*pp) {
                if ((*pp)->last_stream == stream) {
                    picked = *pp;
                    *pp = picked->next;
                    break;
                }
                pp = &(*pp)->next;
            }
            // Fallback: any entry — cross-stream reuse, do waitEvent.
            if (!picked) {
                picked = it->second;
                it->second = picked->next;
            }
            void* ptr = picked->ptr;
            cudaEvent_t e = picked->release_event;
            cudaStream_t prev = picked->last_stream;
            delete picked;
            if (it->second == nullptr) comm_free_lists.erase(it);

            comm_in_use[ptr] = sc;
            comm_total_pooled_bytes -= sc;

            if (e) {
                if (prev != stream && prev != nullptr) {
                    // Cross-stream: ensure releaser's pending ops complete
                    // before the new owner's first op on this VA.
                    CUDA_CHECK(cudaStreamWaitEvent(stream, e, 0));
                }
                event_release_impl(e);
            }
            CUDA_CHECK(cudaMemsetAsync(ptr, 0, sc, stream));
            return ptr;
        }

        // Pool miss: synchronous cudaMalloc keeps VA stable for IB MR.
        void* ptr = nullptr;
        CUDA_CHECK(cudaMalloc(&ptr, sc));
        CUDA_CHECK(cudaMemsetAsync(ptr, 0, sc, stream));
        comm_in_use[ptr] = sc;
        return ptr;
    }

    void comm_pool_release(void* ptr, cudaStream_t stream) {
        if (ptr == nullptr) return;
        auto it = comm_in_use.find(ptr);
        if (it == comm_in_use.end()) return;

        size_t sc = it->second;
        comm_in_use.erase(it);

        // Record an event on the releaser's stream so the next acquirer can
        // cudaStreamWaitEvent on it. stream==0 means the caller already
        // synced (free_comm_buffer doesn't remap; this is the convention).
        // Map to memory_stream so the event is still recordable; same-stream
        // optimization still works for callers consistently passing 0.
        if (stream == 0) stream = memory_stream;
        cudaEvent_t e = event_acquire_impl();
        CUDA_CHECK(cudaEventRecord(e, stream));

        PoolEntry* entry = new PoolEntry{ptr, sc, nullptr, e, stream};
        auto& head = comm_free_lists[sc];
        entry->next = head;
        head = entry;
        comm_total_pooled_bytes += sc;

        if (comm_total_pooled_bytes > comm_pool_cap_bytes) {
            if (!comm_pool_cap_warned) {
                std::cerr << "[GPUManager] comm pool exceeded cap ("
                          << comm_total_pooled_bytes << " > " << comm_pool_cap_bytes
                          << " bytes); evicting. Set AMSS_COMM_POOL_CAP_MB to tune."
                          << std::endl;
                comm_pool_cap_warned = true;
            }
            comm_pool_drain();
        }
    }

    void comm_pool_drain() {
        size_t target = comm_pool_cap_bytes / 2;
        if (comm_total_pooled_bytes <= target || comm_free_lists.empty()) return;

        // Released buffers carry a release_event; their owner's stream may
        // still be running. cudaFree is synchronous and triggers UCX rcache
        // invalidation, so freeing an in-flight buffer corrupts the next
        // cudaStreamSynchronize. Sync the device once before freeing.
        CUDA_CHECK(cudaDeviceSynchronize());

        while (comm_total_pooled_bytes > target && !comm_free_lists.empty()) {
            size_t max_sc = 0;
            for (auto& kv : comm_free_lists) {
                if (kv.second && kv.first > max_sc) max_sc = kv.first;
            }
            auto it = comm_free_lists.find(max_sc);
            if (it == comm_free_lists.end() || it->second == nullptr) break;

            PoolEntry* e = it->second;
            it->second = e->next;
            comm_total_pooled_bytes -= e->capacity;
            if (e->release_event) event_release_impl(e->release_event);
            CUDA_CHECK(cudaFree(e->ptr));
            delete e;

            if (it->second == nullptr) comm_free_lists.erase(it);
        }
    }

    void comm_pool_drain_all() {
        if (comm_free_lists.empty()) {
            comm_total_pooled_bytes = 0;
            return;
        }
        CUDA_CHECK(cudaDeviceSynchronize());
        for (auto& kv : comm_free_lists) {
            PoolEntry* e = kv.second;
            while (e) {
                PoolEntry* next = e->next;
                if (e->release_event) event_release_impl(e->release_event);
                cudaFree(e->ptr);
                delete e;
                e = next;
            }
        }
        comm_free_lists.clear();
        comm_total_pooled_bytes = 0;
    }

    // --- Event pool ---
    std::vector<cudaEvent_t> event_pool;

    cudaEvent_t event_acquire_impl() {
        if (!event_pool.empty()) {
            cudaEvent_t e = event_pool.back();
            event_pool.pop_back();
            return e;
        }
        cudaEvent_t e;
        CUDA_CHECK(cudaEventCreateWithFlags(&e, cudaEventDisableTiming));
        return e;
    }

    void event_release_impl(cudaEvent_t e) {
        if (e) event_pool.push_back(e);
    }

    void event_pool_drain_all() {
        for (auto e : event_pool) cudaEventDestroy(e);
        event_pool.clear();
    }
};

// 单例获取
GPUManager& GPUManager::getInstance() {
    static GPUManager instance;
    return instance;
}

// 构造与析构
GPUManager::GPUManager() : pimpl(new Impl()) {
    for (int i = 0; i < Impl::NUM_STREAMS; ++i) {
        CUDA_CHECK(cudaStreamCreate(&pimpl->stream_pool[i]));
    }
    CUDA_CHECK(cudaStreamCreate(&pimpl->memory_stream));

    if (const char* env = std::getenv("AMSS_COMM_POOL_CAP_MB")) {
        size_t mb = static_cast<size_t>(std::strtoull(env, nullptr, 10));
        if (mb > 0) pimpl->comm_pool_cap_bytes = mb * 1024ULL * 1024ULL;
    }
}

GPUManager::~GPUManager() {
    pimpl->pool_drain_all();
    pimpl->comm_pool_drain_all();
    pimpl->event_pool_drain_all();
    // Free any buffers still in use (all GPU work should be done at this point)
    cudaDeviceSynchronize();
    for (auto& kv : pimpl->in_use) {
        cudaFree(kv.first);
    }
    pimpl->in_use.clear();
    for (auto& kv : pimpl->comm_in_use) {
        cudaFree(kv.first);
    }
    pimpl->comm_in_use.clear();
    for (int i = 0; i < Impl::NUM_STREAMS; ++i) {
        cudaStreamDestroy(pimpl->stream_pool[i]);
    }
    cudaStreamDestroy(pimpl->memory_stream);
    delete pimpl;
}

template<typename T>
T* GPUManager::allocate_device_memory(size_t num_elements, cudaStream_t stream, bool pooled) {
    if (stream == 0) {
        stream = pimpl->memory_stream;
    }
    size_t bytes = num_elements * sizeof(T);
    if (!pooled) {
        void* raw = nullptr;
        CUDA_CHECK(cudaMallocAsync(&raw, bytes, stream));
        CUDA_CHECK(cudaMemsetAsync(raw, 0, bytes, stream));
        return static_cast<T*>(raw);
    }
    void* raw = pimpl->pool_acquire(bytes, stream);
    return static_cast<T*>(raw);
}

template<typename T>
void GPUManager::free_device_memory(T* d_ptr, size_t num_elements, cudaStream_t stream, bool pooled) {
    if (d_ptr == nullptr) return;
    // Mirror allocate_device_memory's stream==0 -> memory_stream remap so the
    // event we record in pool_release lives on the same stream the alloc used.
    if (stream == 0) {
        stream = pimpl->memory_stream;
    }
    if (!pooled) {
        CUDA_CHECK(cudaFreeAsync(static_cast<void*>(d_ptr), stream));
        return;
    }
    pimpl->pool_release(static_cast<void*>(d_ptr), stream);
}

template<typename T>
T* GPUManager::allocate_comm_buffer(size_t num_elements, cudaStream_t stream) {
    if (stream == 0) {
        stream = pimpl->memory_stream;
    }
    size_t bytes = num_elements * sizeof(T);
    if (bytes == 0) bytes = sizeof(T);  // pool always returns a valid pointer
    void* raw = pimpl->comm_pool_acquire(bytes, stream);
    return static_cast<T*>(raw);
}

template<typename T>
void GPUManager::free_comm_buffer(T* d_ptr, size_t /*num_elements*/, cudaStream_t stream) {
    if (d_ptr == nullptr) return;
    // Pass the releaser's stream through so comm_pool_release can record an
    // event for cross-stream reuse sync. stream == 0 means the caller has
    // already synced (no pending ops); comm_pool_release maps it to
    // memory_stream so the recorded event is still well-defined.
    pimpl->comm_pool_release(static_cast<void*>(d_ptr), stream);
}

template<typename T>
void GPUManager::memset_device_memory(T* d_ptr, size_t num_elements, cudaStream_t stream) {
    if (d_ptr == nullptr) return;
    if (stream == 0) {
        stream = pimpl->memory_stream;
    }
    CUDA_CHECK(cudaMemsetAsync(d_ptr, 0, num_elements * sizeof(T), stream));
}

cudaStream_t GPUManager::get_stream() {
    unsigned int cur = pimpl->stream_idx.fetch_add(1);
    return pimpl->stream_pool[cur % Impl::NUM_STREAMS];
}

template<typename T>
void GPUManager::sync_to_gpu(const T* h_ptr, T* d_ptr, size_t num_elements, cudaStream_t stream) {
    if (stream == 0) {
        stream = getInstance().pimpl->memory_stream;
    }
    CUDA_CHECK(cudaMemcpyAsync(d_ptr, h_ptr, num_elements * sizeof(T), cudaMemcpyHostToDevice, stream));
}

template<typename T>
void GPUManager::sync_to_cpu(T* h_ptr, const T* d_ptr, size_t num_elements, cudaStream_t stream) {
    if (stream == 0) {
        stream = getInstance().pimpl->memory_stream;
    }
    CUDA_CHECK(cudaMemcpyAsync(h_ptr, d_ptr, num_elements * sizeof(T), cudaMemcpyDeviceToHost, stream));
}

void GPUManager::synchronize_memory() {
    CUDA_CHECK(cudaStreamSynchronize(pimpl->memory_stream));
}

cudaStream_t GPUManager::memory_stream() const {
    return pimpl->memory_stream;
}

void GPUManager::fence_pool_after_memory() {
    // Record one event capturing memory_stream's current queue, then have every
    // pool stream cudaStreamWaitEvent it. No host block.
    cudaEvent_t e = pimpl->event_acquire_impl();
    CUDA_CHECK(cudaEventRecord(e, pimpl->memory_stream));
    for (int i = 0; i < Impl::NUM_STREAMS; ++i)
        CUDA_CHECK(cudaStreamWaitEvent(pimpl->stream_pool[i], e, 0));
    pimpl->event_release_impl(e);  // waitEvent ops already captured e's state
}

void GPUManager::synchronize_all() {
    CUDA_CHECK(cudaDeviceSynchronize());
}

void GPUManager::synchronize_pool() {
    for (int i = 0; i < Impl::NUM_STREAMS; ++i)
        CUDA_CHECK(cudaStreamSynchronize(pimpl->stream_pool[i]));
    CUDA_CHECK(cudaStreamSynchronize(pimpl->memory_stream));
}

void GPUManager::synchronize_stream(cudaStream_t stream) {
    CUDA_CHECK(cudaStreamSynchronize(stream));
}

void GPUManager::drain_buffer_pool() {
    pimpl->pool_drain_all();
}

void GPUManager::set_buffer_pool_cap(size_t bytes) {
    pimpl->pool_cap_bytes = bytes;
}

void GPUManager::drain_comm_pool() {
    pimpl->comm_pool_drain_all();
}

void GPUManager::set_comm_pool_cap(size_t bytes) {
    pimpl->comm_pool_cap_bytes = bytes;
}

template double* GPUManager::allocate_device_memory<double>(size_t num_elements, cudaStream_t stream, bool pooled);
template int* GPUManager::allocate_device_memory<int>(size_t num_elements, cudaStream_t stream, bool pooled);
template void GPUManager::free_device_memory<double>(double* d_ptr, size_t num_elements, cudaStream_t stream, bool pooled);
template void GPUManager::free_device_memory<int>(int* d_ptr, size_t num_elements, cudaStream_t stream, bool pooled);
template double* GPUManager::allocate_comm_buffer<double>(size_t num_elements, cudaStream_t stream);
template int* GPUManager::allocate_comm_buffer<int>(size_t num_elements, cudaStream_t stream);
template void GPUManager::free_comm_buffer<double>(double* d_ptr, size_t num_elements, cudaStream_t stream);
template void GPUManager::free_comm_buffer<int>(int* d_ptr, size_t num_elements, cudaStream_t stream);
template void GPUManager::memset_device_memory<double>(double* d_ptr, size_t num_elements, cudaStream_t stream);
template void GPUManager::memset_device_memory<int>(int* d_ptr, size_t num_elements, cudaStream_t stream);
template void GPUManager::sync_to_gpu<double>(const double* h_ptr, double* d_ptr, size_t num_elements, cudaStream_t stream);
template void GPUManager::sync_to_gpu<int>(const int* h_ptr, int* d_ptr, size_t num_elements, cudaStream_t stream);
template void GPUManager::sync_to_cpu<double>(double* h_ptr, const double* d_ptr, size_t num_elements, cudaStream_t stream);
template void GPUManager::sync_to_cpu<int>(int* h_ptr, const int* d_ptr, size_t num_elements, cudaStream_t stream);
template unsigned long long* GPUManager::allocate_device_memory<unsigned long long>(size_t num_elements, cudaStream_t stream, bool pooled);
template void GPUManager::free_device_memory<unsigned long long>(unsigned long long* d_ptr, size_t num_elements, cudaStream_t stream, bool pooled);
template void GPUManager::sync_to_gpu<unsigned long long>(const unsigned long long* h_ptr, unsigned long long* d_ptr, size_t num_elements, cudaStream_t stream);
template void GPUManager::sync_to_cpu<unsigned long long>(unsigned long long* h_ptr, const unsigned long long* d_ptr, size_t num_elements, cudaStream_t stream);
template unsigned long* GPUManager::allocate_device_memory<unsigned long>(size_t num_elements, cudaStream_t stream, bool pooled);
template void GPUManager::free_device_memory<unsigned long>(unsigned long* d_ptr, size_t num_elements, cudaStream_t stream, bool pooled);
template void GPUManager::sync_to_gpu<unsigned long>(const unsigned long * h_ptr, unsigned long* d_ptr, size_t num_elements, cudaStream_t stream);
template void GPUManager::sync_to_cpu<unsigned long>(unsigned long* h_ptr, const unsigned long* d_ptr, size_t num_elements, cudaStream_t stream);

__global__ void unpack_padding_kernel(
    const double* __restrict__ src_unpadded, double* __restrict__ dst_shifted,
    int nx, int ny, int nz, int px, int py
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= nx * ny * nz) return;
    int i, j, k;
    linear_to_ijk(tid, nx, ny, i, j, k);

    dst_shifted[(k * py + j) * px + i] = src_unpadded[(k * ny + j) * nx + i];
}

__global__ void pack_padding_kernel(
    const double* __restrict__ src_shifted, double* __restrict__ dst_unpadded,
    int nx, int ny, int nz, int px, int py
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= nx * ny * nz) return;
    int i, j, k;
    linear_to_ijk(tid, nx, ny, i, j, k);

    dst_unpadded[(k * ny + j) * nx + i] = src_shifted[(k * py + j) * px + i];
}

void GPUManager::unpack_to_pad(const double* src, double* dst, int shape[3], int shape_pad[3], cudaStream_t stream) {
    if (stream == 0) {
        stream = pimpl->memory_stream;
    }

    int total = shape[0] * shape[1] * shape[2];
    dim3 block(256);
    dim3 grid((total + 255) / 256);

    unpack_padding_kernel<<<grid, block, 0, stream>>>(
        src, dst,
        shape[0], shape[1], shape[2],
        shape_pad[0], shape_pad[1]
    );
}

void GPUManager::pack_to_pad(const double* src, double* dst, int shape[3], int shape_pad[3], cudaStream_t stream) {
    if (stream == 0) {
        stream = pimpl->memory_stream;
    }

    int total = shape[0] * shape[1] * shape[2];
    dim3 block(256);
    dim3 grid((total + 255) / 256);

    pack_padding_kernel<<<grid, block, 0, stream>>>(
        src, dst,
        shape[0], shape[1], shape[2],
        shape_pad[0], shape_pad[1]
    );
}

cudaEvent_t GPUManager::event_acquire() {
    return pimpl->event_acquire_impl();
}

void GPUManager::event_release(cudaEvent_t e) {
    pimpl->event_release_impl(e);
}

void GPUManager::event_record(cudaEvent_t e, cudaStream_t stream) {
    CUDA_CHECK(cudaEventRecord(e, stream));
}

void GPUManager::event_host_wait(cudaEvent_t e) {
    CUDA_CHECK(cudaEventSynchronize(e));
}
