#include "Parallel.h"

#include "gpu_runtime_coordinator.h"
#include "mpi_cuda_aware_runtime.h"
#include "prolongrestrict.h"
#include "misc.h"
#include "parameters.h"
#include "fmisc.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <vector>

#include "MPatch.h"
#include "data_motion_helper.h"

namespace {

struct GridOffset {
    int x;
    int y;
    int z;
};

struct StencilTransferBatch {
    cudaStream_t stream;
    std::vector<GpuStencilTransfer> transfers;
};

struct UnpackTransferBatch {
    cudaStream_t stream;
    int max_blocks_x;
    std::vector<GpuUnpackWindowTransfer> transfers;
};

struct InterpScratch {
    double* d_local_shellf;
    int* d_local_weight;
    int* d_global_weight;
    int* d_err_flag;
    double* h_local_shellf;
    int* h_local_weight;
    int shellf_capacity;
    int weight_capacity;
};

InterpScratch& interp_scratch() {
    static InterpScratch scratch = {0, 0, 0, 0, 0, 0, 0, 0};
    return scratch;
}

void warn_gpu_cpu_fallback_once(const char* context) {
    static bool warned = false;
    if (warned) return;
    warned = true;

    int myrank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank == 0) {
        cout << "WARNING: GPU path is using a CPU/host fallback: " << context
             << ". Verify host/device synchronization and performance before production use." << endl;
    }
}

void report_cuda_ipc_transfer_once() {
    static bool reported = false;
    if (reported) return;
    reported = true;

    int myrank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank == 0) {
        cout << " MPI CUDA-aware status: runtime CUDA buffers unavailable; using same-node CUDA IPC/P2P grid transfer" << endl;
    }
}

void report_cuda_ipc_disabled_once(const char* context) {
    static bool reported = false;
    if (reported) return;
    reported = true;

    int myrank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank == 0) {
        cout << " MPI CUDA IPC/P2P grid transfer unavailable: " << context
             << ". Falling back to host-staged MPI grid transfer." << endl;
    }
}

bool cuda_ipc_disabled_by_environment() {
    const char* value = std::getenv("AMSS_NCKU_DISABLE_CUDA_IPC");
    if (!value) {
        return false;
    }

    return std::strcmp(value, "0") != 0
        && std::strcmp(value, "false") != 0
        && std::strcmp(value, "FALSE") != 0
        && std::strcmp(value, "off") != 0
        && std::strcmp(value, "OFF") != 0;
}

bool all_world_ranks_share_node(int cpusize) {
    static int cached_cpusize = 0;
    static bool cached_result = false;
    if (cached_cpusize == cpusize) {
        return cached_result;
    }

    MPI_Comm shared_comm = MPI_COMM_NULL;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shared_comm);
    int shared_size = 0;
    MPI_Comm_size(shared_comm, &shared_size);
    MPI_Comm_free(&shared_comm);
    cached_cpusize = cpusize;
    cached_result = shared_size == cpusize;
    return cached_result;
}

void ensure_interp_scratch_capacity(InterpScratch& scratch, int point_count, int variable_count, cudaStream_t stream) {
    const int required_shellf_capacity = point_count * variable_count;
    if (required_shellf_capacity > scratch.shellf_capacity) {
        GpuRuntimeCoordinator::instance().release_device_buffer(scratch.d_local_shellf, scratch.shellf_capacity, stream);
        scratch.d_local_shellf = GpuRuntimeCoordinator::instance().allocate_device_buffer<double>(required_shellf_capacity, stream, false);
        delete[] scratch.h_local_shellf;
        scratch.h_local_shellf = new double[required_shellf_capacity];
        scratch.shellf_capacity = required_shellf_capacity;
    }

    if (point_count > scratch.weight_capacity) {
        GpuRuntimeCoordinator::instance().release_device_buffer(scratch.d_local_weight, scratch.weight_capacity, stream);
        GpuRuntimeCoordinator::instance().release_device_buffer(scratch.d_global_weight, scratch.weight_capacity, stream);
        scratch.d_local_weight = GpuRuntimeCoordinator::instance().allocate_device_buffer<int>(point_count, stream, false);
        scratch.d_global_weight = GpuRuntimeCoordinator::instance().allocate_device_buffer<int>(point_count, stream, false);
        delete[] scratch.h_local_weight;
        scratch.h_local_weight = new int[point_count];
        scratch.weight_capacity = point_count;
    }

    if (!scratch.d_err_flag) {
        scratch.d_err_flag = GpuRuntimeCoordinator::instance().allocate_device_buffer<int>(1, stream, false);
    }
}

bool variable_lists_match(MyList<var>* source_variables, MyList<var>* target_variables) {
    while (source_variables && target_variables) {
        source_variables = source_variables->next;
        target_variables = target_variables->next;
    }
    return source_variables == nullptr && target_variables == nullptr;
}

int classify_transfer_mode(const Parallel::gridseg* source_segment, const Parallel::gridseg* target_segment) {
    if (source_segment->Bg->lev == target_segment->Bg->lev) {
        return 1;
    }
    return (source_segment->Bg->lev > target_segment->Bg->lev) ? 2 : 3;
}

int segment_volume(const Parallel::gridseg* segment) {
    return segment->shape[0] * segment->shape[1] * segment->shape[2];
}

GridOffset grid_offset_from_bounds(const Parallel::gridseg* segment, Block* block) {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
    const double dx = block->shape[0] > 1 ? (block->bbox[3] - block->bbox[0]) / (block->shape[0] - 1) : 1.0;
    const double dy = block->shape[1] > 1 ? (block->bbox[4] - block->bbox[1]) / (block->shape[1] - 1) : 1.0;
    const double dz = block->shape[2] > 1 ? (block->bbox[5] - block->bbox[2]) / (block->shape[2] - 1) : 1.0;
#else
#ifdef Cell
    const double dx = (block->bbox[3] - block->bbox[0]) / block->shape[0];
    const double dy = (block->bbox[4] - block->bbox[1]) / block->shape[1];
    const double dz = (block->bbox[5] - block->bbox[2]) / block->shape[2];
#else
#error Not define Vertex nor Cell
#endif
#endif
    return GridOffset{
        static_cast<int>(std::trunc((segment->llb[0] - block->bbox[0]) / dx + 0.4)),
        static_cast<int>(std::trunc((segment->llb[1] - block->bbox[1]) / dy + 0.4)),
        static_cast<int>(std::trunc((segment->llb[2] - block->bbox[2]) / dz + 0.4))
    };
}

void pack_same_level_segment(const Parallel::gridseg* source_segment, const Parallel::gridseg* target_segment, double* device_source, double* device_target) {
    const GridOffset offset = grid_offset_from_bounds(target_segment, source_segment->Bg);
    gpu_pack_launch(
        source_segment->Bg->stream,
        device_source,
        device_target,
        source_segment->Bg->shape[0],
        source_segment->Bg->shape[1],
        target_segment->shape[0],
        target_segment->shape[1],
        target_segment->shape[2],
        offset.x,
        offset.y,
        offset.z
    );
}

void unpack_same_level_segment(const Parallel::gridseg* source_segment, const Parallel::gridseg* target_segment, double* device_source, double* device_target) {
    const GridOffset offset = grid_offset_from_bounds(target_segment, target_segment->Bg);
    gpu_unpack_launch(
        target_segment->Bg->stream,
        device_source,
        device_target,
        target_segment->Bg->shape[0],
        target_segment->Bg->shape[1],
        target_segment->shape[0],
        target_segment->shape[1],
        target_segment->shape[2],
        offset.x,
        offset.y,
        offset.z
    );
}

void pack_same_level_fields(
    const Parallel::gridseg* source_segment,
    const Parallel::gridseg* target_segment,
    const std::vector<const double*>& source_fields,
    const std::vector<double*>& target_fields
) {
    const GridOffset offset = grid_offset_from_bounds(target_segment, source_segment->Bg);
    gpu_pack_fields_launch(
        source_segment->Bg->stream,
        source_fields.data(),
        target_fields.data(),
        static_cast<int>(source_fields.size()),
        source_segment->Bg->shape[0],
        source_segment->Bg->shape[1],
        target_segment->shape[0],
        target_segment->shape[1],
        target_segment->shape[2],
        offset.x,
        offset.y,
        offset.z
    );
}

void unpack_same_level_fields(
    const Parallel::gridseg* source_segment,
    const Parallel::gridseg* target_segment,
    const std::vector<const double*>& source_fields,
    const std::vector<double*>& target_fields
) {
    const GridOffset offset = grid_offset_from_bounds(target_segment, target_segment->Bg);
    gpu_unpack_fields_launch(
        target_segment->Bg->stream,
        source_fields.data(),
        target_fields.data(),
        static_cast<int>(source_fields.size()),
        target_segment->Bg->shape[0],
        target_segment->Bg->shape[1],
        target_segment->shape[0],
        target_segment->shape[1],
        target_segment->shape[2],
        offset.x,
        offset.y,
        offset.z
    );
}

bool same_stencil_launch_family(const GpuStencilTransfer& left, const GpuStencilTransfer& right) {
    return left.ni == right.ni && left.nj == right.nj && left.nk == right.nk
        && left.i_start == right.i_start && left.j_start == right.j_start && left.k_start == right.k_start
        && left.lbc0 == right.lbc0 && left.lbc1 == right.lbc1 && left.lbc2 == right.lbc2
        && left.lbf0 == right.lbf0 && left.lbf1 == right.lbf1 && left.lbf2 == right.lbf2
        && left.extc0 == right.extc0 && left.extc1 == right.extc1 && left.extc2 == right.extc2
        && left.extf0 == right.extf0 && left.extf1 == right.extf1 && left.extf2 == right.extf2;
}

StencilTransferBatch& stencil_transfer_batch(
    std::vector<StencilTransferBatch>& batches,
    cudaStream_t stream,
    const GpuStencilTransfer& transfer
) {
    for (std::vector<StencilTransferBatch>::iterator batch = batches.begin(); batch != batches.end(); ++batch) {
        if (batch->stream == stream && !batch->transfers.empty() && same_stencil_launch_family(batch->transfers.front(), transfer)) {
            return *batch;
        }
    }

    StencilTransferBatch batch = {stream, std::vector<GpuStencilTransfer>()};
    batches.push_back(batch);
    return batches.back();
}

void append_stencil_transfer_batch(
    std::vector<StencilTransferBatch>& batches,
    cudaStream_t stream,
    const GpuStencilTransfer& transfer
) {
    StencilTransferBatch& batch = stencil_transfer_batch(batches, stream, transfer);
    batch.transfers.push_back(transfer);
}

UnpackTransferBatch& unpack_transfer_batch(
    std::vector<UnpackTransferBatch>& batches,
    cudaStream_t stream
) {
    for (std::vector<UnpackTransferBatch>::iterator batch = batches.begin(); batch != batches.end(); ++batch) {
        if (batch->stream == stream) {
            return *batch;
        }
    }

    UnpackTransferBatch batch = {stream, 0, std::vector<GpuUnpackWindowTransfer>()};
    batches.push_back(batch);
    return batches.back();
}

void append_unpack_transfer_batch(
    std::vector<UnpackTransferBatch>& batches,
    cudaStream_t stream,
    const GpuUnpackWindowTransfer& transfer,
    int blocks_x
) {
    UnpackTransferBatch& batch = unpack_transfer_batch(batches, stream);
    if (blocks_x > batch.max_blocks_x) {
        batch.max_blocks_x = blocks_x;
    }
    batch.transfers.push_back(transfer);
}

void append_unique_stream(std::vector<cudaStream_t>& streams, cudaStream_t stream) {
    for (std::vector<cudaStream_t>::const_iterator existing = streams.begin(); existing != streams.end(); ++existing) {
        if (*existing == stream) {
            return;
        }
    }

    streams.push_back(stream);
}

void zero_transfer_segment(double* device_target, int element_count, cudaStream_t stream) {
    if (element_count <= 0) {
        return;
    }
    CUDA_CHECK(cudaMemsetAsync(device_target, 0, element_count * sizeof(double), stream));
}

void gpu_transfer_single_rank_local(
    MyList<Parallel::gridseg> **src,
    MyList<Parallel::gridseg> **dst,
    MyList<var> *VarList1,
    MyList<var> *VarList2,
    int Symmetry,
    int myrank
) {
    const int length = Parallel::gpu_data_packer(0, src[myrank], dst[myrank], myrank, PACK, VarList1, VarList2, Symmetry);
    if (!length) {
        return;
    }

    double *local_data = GpuRuntimeCoordinator::instance().allocate_device_buffer<double>(length, 0, false);
    GpuRuntimeCoordinator::instance().wait_for_memory_ops();
    Parallel::gpu_data_packer(local_data, src[myrank], dst[myrank], myrank, PACK, VarList1, VarList2, Symmetry);
    Parallel::gpu_data_packer(local_data, src[myrank], dst[myrank], myrank, UNPACK, VarList1, VarList2, Symmetry);
    GpuRuntimeCoordinator::instance().release_device_buffer(local_data, 0);
}

void gpu_transfer_host_staged(
    MyList<Parallel::gridseg> **src, MyList<Parallel::gridseg> **dst,
    MyList<var> *VarList1,
    MyList<var> *VarList2,
    int Symmetry
) {
    warn_gpu_cpu_fallback_once("host-staged MPI grid transfer because CUDA-aware MPI is unavailable");

    int myrank, cpusize;
    MPI_Comm_size(MPI_COMM_WORLD, &cpusize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (cpusize == 1) {
        gpu_transfer_single_rank_local(src, dst, VarList1, VarList2, Symmetry, myrank);
        return;
    }

    int node;
    MPI_Request *reqs = new MPI_Request[2 * cpusize];
    MPI_Status *stats = new MPI_Status[2 * cpusize];
    int req_no = 0;

    double **send_data_d = new double *[cpusize];
    double **rec_data_d  = new double *[cpusize];
    double **send_data_h = new double *[cpusize];
    double **rec_data_h  = new double *[cpusize];
    int length;

    for (node = 0; node < cpusize; node++) {
        send_data_d[node] = rec_data_d[node] = nullptr;
        send_data_h[node] = rec_data_h[node] = nullptr;

        if (node == myrank) {
            length = gpu_data_packer(0, src[myrank], dst[myrank], node, PACK, VarList1, VarList2, Symmetry);
            if (length) {
                rec_data_d[node] = GpuRuntimeCoordinator::instance().allocate_device_buffer<double>(length, 0, false);
                GpuRuntimeCoordinator::instance().wait_for_memory_ops();
                gpu_data_packer(rec_data_d[node], src[myrank], dst[myrank], node, PACK, VarList1, VarList2, Symmetry);
            }
        }
        else {
            length = gpu_data_packer(0, src[myrank], dst[myrank], node, PACK, VarList1, VarList2, Symmetry);
            if (length) {
                send_data_d[node] = GpuRuntimeCoordinator::instance().allocate_device_buffer<double>(length, 0, false);
                send_data_h[node] = new double[length];
                GpuRuntimeCoordinator::instance().wait_for_memory_ops();
                gpu_data_packer(send_data_d[node], src[myrank], dst[myrank], node, PACK, VarList1, VarList2, Symmetry);
                GpuRuntimeCoordinator::instance().copy_to_host(send_data_h[node], send_data_d[node], length);
                GpuRuntimeCoordinator::instance().wait_for_memory_ops();
                MPI_Isend((void *)send_data_h[node], length, MPI_DOUBLE, node, 1, MPI_COMM_WORLD, reqs + req_no++);
            }
            length = gpu_data_packer(0, src[node], dst[node], node, UNPACK, VarList1, VarList2, Symmetry);
            if (length) {
                rec_data_d[node] = GpuRuntimeCoordinator::instance().allocate_device_buffer<double>(length, 0, false);
                rec_data_h[node] = new double[length];
                MPI_Irecv((void *)rec_data_h[node], length, MPI_DOUBLE, node, 1, MPI_COMM_WORLD, reqs + req_no++);
            }
        }
    }

    MPI_Waitall(req_no, reqs, stats);

    for (node = 0; node < cpusize; node++) {
        if (node == myrank) {
            if (rec_data_d[node]) {
                gpu_data_packer(rec_data_d[node], src[node], dst[node], node, UNPACK, VarList1, VarList2, Symmetry);
            }
        }
        else {
            if (rec_data_h[node]) {
                length = gpu_data_packer(0, src[node], dst[node], node, UNPACK, VarList1, VarList2, Symmetry);
                GpuRuntimeCoordinator::instance().copy_to_device(rec_data_h[node], rec_data_d[node], length);
                GpuRuntimeCoordinator::instance().wait_for_memory_ops();
                gpu_data_packer(rec_data_d[node], src[node], dst[node], node, UNPACK, VarList1, VarList2, Symmetry);
            }
        }
    }

    for (node = 0; node < cpusize; node++) {
        if (send_data_d[node]) GpuRuntimeCoordinator::instance().release_device_buffer(send_data_d[node], 0);
        if (rec_data_d[node]) GpuRuntimeCoordinator::instance().release_device_buffer(rec_data_d[node], 0);
        if (send_data_h[node]) delete[] send_data_h[node];
        if (rec_data_h[node]) delete[] rec_data_h[node];
    }

    delete[] reqs;
    delete[] stats;
    delete[] send_data_d;
    delete[] rec_data_d;
    delete[] send_data_h;
    delete[] rec_data_h;
}

struct CudaIpcTransferDescriptor {
    int length;
    int valid;
    int generation;
    cudaIpcMemHandle_t handle;
};

struct CudaIpcPeerBuffer {
    double* send_data;
    int send_capacity;
    int send_generation;
    cudaIpcMemHandle_t send_handle;
    double* rec_data;
    int rec_capacity;
    void* remote_ptr;
    int remote_generation;
    std::vector<double*> retired_send_buffers;

    CudaIpcPeerBuffer()
        : send_data(nullptr),
          send_capacity(0),
          send_generation(0),
          rec_data(nullptr),
          rec_capacity(0),
          remote_ptr(nullptr),
          remote_generation(0) {
        std::memset(&send_handle, 0, sizeof(send_handle));
    }
};

struct CudaIpcGridTransferCache {
    int cpusize;
    std::vector<CudaIpcPeerBuffer> peers;

    CudaIpcGridTransferCache() : cpusize(0), peers() {}
};

CudaIpcGridTransferCache& cuda_ipc_grid_transfer_cache() {
    static CudaIpcGridTransferCache cache;
    return cache;
}

int grow_capacity(int current_capacity, int required_capacity) {
    int new_capacity = current_capacity > 0 ? current_capacity : 1024;
    while (new_capacity < required_capacity) {
        if (new_capacity > 1073741823) {
            return required_capacity;
        }
        new_capacity *= 2;
    }
    return new_capacity;
}

bool ensure_cuda_ipc_send_capacity(CudaIpcPeerBuffer& peer, int required_capacity, const char** failure_stage) {
    if (required_capacity <= peer.send_capacity) {
        return true;
    }

    const int new_capacity = grow_capacity(peer.send_capacity, required_capacity);
    double* new_buffer = nullptr;
    cudaError_t allocation_error = cudaMalloc((void**)&new_buffer, new_capacity * sizeof(double));
    if (allocation_error != cudaSuccess) {
        *failure_stage = "cudaMalloc failed while growing CUDA IPC send buffer";
        return false;
    }

    cudaIpcMemHandle_t new_handle;
    cudaError_t handle_error = cudaIpcGetMemHandle(&new_handle, new_buffer);
    if (handle_error != cudaSuccess) {
        cudaFree(new_buffer);
        *failure_stage = "cudaIpcGetMemHandle failed while growing CUDA IPC send buffer";
        return false;
    }

    if (peer.send_data) {
        peer.retired_send_buffers.push_back(peer.send_data);
    }
    peer.send_data = new_buffer;
    peer.send_capacity = new_capacity;
    peer.send_generation += 1;
    peer.send_handle = new_handle;
    return true;
}

bool ensure_cuda_ipc_receive_capacity(CudaIpcPeerBuffer& peer, int required_capacity, const char** failure_stage) {
    if (required_capacity <= peer.rec_capacity) {
        return true;
    }

    const int new_capacity = grow_capacity(peer.rec_capacity, required_capacity);
    double* new_buffer = nullptr;
    cudaError_t allocation_error = cudaMalloc((void**)&new_buffer, new_capacity * sizeof(double));
    if (allocation_error != cudaSuccess) {
        *failure_stage = "cudaMalloc failed while growing CUDA IPC receive buffer";
        return false;
    }

    if (peer.rec_data) {
        cudaFree(peer.rec_data);
    }
    peer.rec_data = new_buffer;
    peer.rec_capacity = new_capacity;
    return true;
}

bool gpu_transfer_cuda_ipc(
    MyList<Parallel::gridseg> **src, MyList<Parallel::gridseg> **dst,
    MyList<var> *VarList1,
    MyList<var> *VarList2,
    int Symmetry
) {
    static bool permanently_disabled = false;
    static const char* permanent_disable_reason = nullptr;
    if (permanently_disabled) {
        report_cuda_ipc_disabled_once(
            permanent_disable_reason ? permanent_disable_reason : "previous CUDA IPC setup failed"
        );
        return false;
    }

    if (cuda_ipc_disabled_by_environment()) {
        report_cuda_ipc_disabled_once("disabled by AMSS_NCKU_DISABLE_CUDA_IPC");
        return false;
    }

    int myrank = 0;
    int cpusize = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &cpusize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (cpusize == 1) {
        gpu_transfer_single_rank_local(src, dst, VarList1, VarList2, Symmetry, myrank);
        return true;
    }

    if (!all_world_ranks_share_node(cpusize)) {
        permanently_disabled = true;
        permanent_disable_reason = "MPI ranks are not all on the same node";
        report_cuda_ipc_disabled_once(permanent_disable_reason);
        return false;
    }

    CudaIpcGridTransferCache& cache = cuda_ipc_grid_transfer_cache();
    if (cache.cpusize == 0) {
        cache.cpusize = cpusize;
        cache.peers.resize(cpusize);
    } else if (cache.cpusize != cpusize) {
        permanently_disabled = true;
        permanent_disable_reason = "MPI world size changed after CUDA IPC cache initialization";
        report_cuda_ipc_disabled_once(permanent_disable_reason);
        return false;
    }

    std::vector<int> recv_lengths(cpusize, 0);
    std::vector<CudaIpcTransferDescriptor> local_descriptors(cpusize);
    std::vector<CudaIpcTransferDescriptor> all_descriptors(cpusize * cpusize);

    std::memset(local_descriptors.data(), 0, local_descriptors.size() * sizeof(CudaIpcTransferDescriptor));

    int local_ok = 1;
    const char* failure_stage = nullptr;

    bool copy_launched = false;
    for (int node = 0; node < cpusize; ++node) {
        if (node == myrank) {
            const int length = Parallel::gpu_data_packer(0, src[myrank], dst[myrank], node, PACK, VarList1, VarList2, Symmetry);
            recv_lengths[node] = length;
            if (length) {
                CudaIpcPeerBuffer& peer = cache.peers[node];
                if (!ensure_cuda_ipc_receive_capacity(peer, length, &failure_stage)) {
                    local_ok = 0;
                    break;
                }
                Parallel::gpu_data_packer(peer.rec_data, src[myrank], dst[myrank], node, PACK, VarList1, VarList2, Symmetry);
            }
            continue;
        }

        const int send_length = Parallel::gpu_data_packer(0, src[myrank], dst[myrank], node, PACK, VarList1, VarList2, Symmetry);
        local_descriptors[node].length = send_length;
        local_descriptors[node].valid = 0;
        if (send_length) {
            CudaIpcPeerBuffer& peer = cache.peers[node];
            if (!ensure_cuda_ipc_send_capacity(peer, send_length, &failure_stage)) {
                local_ok = 0;
                break;
            }

            Parallel::gpu_data_packer(peer.send_data, src[myrank], dst[myrank], node, PACK, VarList1, VarList2, Symmetry);

            local_descriptors[node].handle = peer.send_handle;
            local_descriptors[node].generation = peer.send_generation;
            local_descriptors[node].valid = 1;
        }

        const int recv_length = Parallel::gpu_data_packer(0, src[node], dst[node], node, UNPACK, VarList1, VarList2, Symmetry);
        recv_lengths[node] = recv_length;
        if (recv_length) {
            CudaIpcPeerBuffer& peer = cache.peers[node];
            if (!ensure_cuda_ipc_receive_capacity(peer, recv_length, &failure_stage)) {
                local_ok = 0;
                break;
            }
        }
    }

    int global_ok = 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (!global_ok) {
        permanently_disabled = true;
        permanent_disable_reason = failure_stage ? failure_stage : "local CUDA IPC setup failed";
        report_cuda_ipc_disabled_once(permanent_disable_reason);
        return false;
    }

    MPI_Allgather(
        local_descriptors.data(),
        static_cast<int>(local_descriptors.size() * sizeof(CudaIpcTransferDescriptor)),
        MPI_BYTE,
        all_descriptors.data(),
        static_cast<int>(local_descriptors.size() * sizeof(CudaIpcTransferDescriptor)),
        MPI_BYTE,
        MPI_COMM_WORLD
    );

    for (int node = 0; node < cpusize; ++node) {
        if (node == myrank || recv_lengths[node] == 0) {
            continue;
        }

        const CudaIpcTransferDescriptor& descriptor = all_descriptors[node * cpusize + myrank];
        if (!descriptor.valid || descriptor.length != recv_lengths[node]) {
            local_ok = 0;
            failure_stage = "CUDA IPC descriptor length/validity mismatch";
            continue;
        }

        CudaIpcPeerBuffer& peer = cache.peers[node];
        if (peer.remote_generation != descriptor.generation || !peer.remote_ptr) {
            if (peer.remote_ptr) {
                cudaIpcCloseMemHandle(peer.remote_ptr);
                peer.remote_ptr = nullptr;
                peer.remote_generation = 0;
            }

            cudaError_t open_error = cudaIpcOpenMemHandle(
                &peer.remote_ptr,
                descriptor.handle,
                cudaIpcMemLazyEnablePeerAccess
            );
            if (open_error != cudaSuccess) {
                local_ok = 0;
                failure_stage = "cudaIpcOpenMemHandle failed for peer send buffer";
                continue;
            }
            peer.remote_generation = descriptor.generation;
        }

        cudaError_t copy_error = cudaMemcpyAsync(
            peer.rec_data,
            peer.remote_ptr,
            recv_lengths[node] * sizeof(double),
            cudaMemcpyDeviceToDevice,
            0
        );
        if (copy_error != cudaSuccess) {
            local_ok = 0;
            failure_stage = "cudaMemcpyAsync failed for CUDA IPC peer copy";
        } else {
            copy_launched = true;
        }
    }

    if (copy_launched) {
        cudaError_t sync_error = cudaStreamSynchronize(0);
        if (sync_error != cudaSuccess) {
            local_ok = 0;
            failure_stage = "cudaStreamSynchronize failed after CUDA IPC peer copies";
        }
    }

    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (!global_ok) {
        permanently_disabled = true;
        permanent_disable_reason = failure_stage ? failure_stage : "CUDA IPC peer copy failed";
        report_cuda_ipc_disabled_once(permanent_disable_reason);
        return false;
    }

    report_cuda_ipc_transfer_once();

    for (int node = 0; node < cpusize; ++node) {
        if (recv_lengths[node]) {
            Parallel::gpu_data_packer(cache.peers[node].rec_data, src[node], dst[node], node, UNPACK, VarList1, VarList2, Symmetry);
        }
    }

    return true;
}

void launch_stencil_transfer_batches(
    const std::vector<StencilTransferBatch>& batches,
    int transfer_mode
) {
    for (std::vector<StencilTransferBatch>::const_iterator batch = batches.begin(); batch != batches.end(); ++batch) {
        if (batch->transfers.empty()) {
            continue;
        }

        if (transfer_mode == 2) {
            gpu_restrict3_batch_launch(batch->stream, batch->transfers.data(), static_cast<int>(batch->transfers.size()));
        } else {
            gpu_prolong3_batch_launch(batch->stream, batch->transfers.data(), static_cast<int>(batch->transfers.size()));
        }
    }
}

void launch_unpack_transfer_batches(const std::vector<UnpackTransferBatch>& batches) {
    for (std::vector<UnpackTransferBatch>::const_iterator batch = batches.begin(); batch != batches.end(); ++batch) {
        if (batch->transfers.empty()) {
            continue;
        }

        gpu_unpack_batch_launch(
            batch->stream,
            batch->transfers.data(),
            static_cast<int>(batch->transfers.size()),
            batch->max_blocks_x
        );
    }
}

void launch_transfer_kernel(
    int transfer_mode,
    int direction,
    Parallel::gridseg* source_segment,
    Parallel::gridseg* target_segment,
    double* device_source,
    double* device_target,
    const double soa[3],
    int symmetry_mode
) {
    if (direction == PACK) {
        switch (transfer_mode) {
            case 1:
                pack_same_level_segment(source_segment, target_segment, device_source, device_target);
                break;
            case 2:
                gpu_restrict3_launch(
                    source_segment->Bg->stream,
                    device_source,
                    device_target,
                    target_segment->llb,
                    target_segment->uub,
                    target_segment->shape,
                    source_segment->Bg->bbox,
                    source_segment->Bg->bbox + dim,
                    source_segment->Bg->shape,
                    target_segment->llb,
                    target_segment->uub,
                    soa,
                    symmetry_mode
                );
                break;
            default:
                gpu_prolong3_launch(
                    source_segment->Bg->stream,
                    device_source,
                    device_target,
                    source_segment->Bg->bbox,
                    source_segment->Bg->bbox + dim,
                    source_segment->Bg->shape,
                    target_segment->llb,
                    target_segment->uub,
                    target_segment->shape,
                    target_segment->llb,
                    target_segment->uub,
                    soa,
                    symmetry_mode
                );
                break;
        }
        return;
    }

    unpack_same_level_segment(source_segment, target_segment, device_source, device_target);
}

}  // namespace

int Parallel::gpu_data_packer(
    double *d_data, MyList<Parallel::gridseg> *src, MyList<Parallel::gridseg> *dst, int rank_in, int dir,
    MyList<var> *VarLists, MyList<var> *VarListd, int Symmetry
) {
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (dir != PACK && dir != UNPACK) {
        cout << "error dir " << dir << " for data_packer " << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int size_out = 0;

    if (!src || !dst) return size_out;

    if (!variable_lists_match(VarLists, VarListd)) {
        cout << "error in short data packer, var lists does not match." << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    const int transfer_mode = classify_transfer_mode(src->data, dst->data);
    std::vector<StencilTransferBatch> stencil_batches;
    std::vector<UnpackTransferBatch> unpack_batches;
    std::vector<cudaStream_t> used_streams;
    
    while (src && dst) {
        if (
            (dir == PACK && dst->data->Bg->rank == rank_in && src->data->Bg->rank == myrank) ||
            (dir == UNPACK && src->data->Bg->rank == rank_in && dst->data->Bg->rank == myrank)
        ) {
            MyList<var> *varls = VarLists;
            MyList<var> *varld = VarListd;
            std::vector<const double*> same_level_source_fields;
            std::vector<double*> same_level_target_fields;
            while (varls && varld) {
                if (d_data) {
                    double* device_source = (dir == PACK)
                        ? src->data->Bg->d_fgfs[varls->data->sgfn]
                        : d_data + size_out;
                    double* device_target = (dir == PACK)
                        ? d_data + size_out
                        : dst->data->Bg->d_fgfs[varld->data->sgfn];
                    if (dir == UNPACK) {
                        same_level_source_fields.push_back(device_source);
                        same_level_target_fields.push_back(device_target);
                        dst->data->Bg->mark_gpu_modified(varld->data->sgfn);
                    } else if (transfer_mode == 1) {
                        same_level_source_fields.push_back(device_source);
                        same_level_target_fields.push_back(device_target);
                    } else if (dir == PACK && transfer_mode == 2) {
                        GpuStencilTransfer transfer;
                        if (gpu_build_restrict3_transfer(
                                device_source,
                                device_target,
                                dst->data->llb,
                                dst->data->uub,
                                dst->data->shape,
                                src->data->Bg->bbox,
                                src->data->Bg->bbox + dim,
                                src->data->Bg->shape,
                                dst->data->llb,
                                dst->data->uub,
                                varls->data->SoA,
                                Symmetry,
                                &transfer
                            )) {
                            append_stencil_transfer_batch(stencil_batches, src->data->Bg->stream, transfer);
                        } else {
                            zero_transfer_segment(device_target, segment_volume(dst->data), src->data->Bg->stream);
                            append_unique_stream(used_streams, src->data->Bg->stream);
                        }
                    } else if (dir == PACK && transfer_mode == 3) {
                        GpuStencilTransfer transfer;
                        if (gpu_build_prolong3_transfer(
                                device_source,
                                device_target,
                                src->data->Bg->bbox,
                                src->data->Bg->bbox + dim,
                                src->data->Bg->shape,
                                dst->data->llb,
                                dst->data->uub,
                                dst->data->shape,
                                dst->data->llb,
                                dst->data->uub,
                                varls->data->SoA,
                                Symmetry,
                                &transfer
                            )) {
                            append_stencil_transfer_batch(stencil_batches, src->data->Bg->stream, transfer);
                        } else {
                            zero_transfer_segment(device_target, segment_volume(dst->data), src->data->Bg->stream);
                            append_unique_stream(used_streams, src->data->Bg->stream);
                        }
                    } else if (transfer_mode != 1) {
                        launch_transfer_kernel(
                            transfer_mode,
                            dir,
                            src->data,
                            dst->data,
                            device_source,
                            device_target,
                            varls->data->SoA,
                            Symmetry
                        );
                    }
                }
                size_out += segment_volume(dst->data);
                varls = varls->next;
                varld = varld->next;
            }

            if (d_data && dir == PACK && transfer_mode == 1 && !same_level_source_fields.empty()) {
                if (dir == PACK) {
                    pack_same_level_fields(src->data, dst->data, same_level_source_fields, same_level_target_fields);
                    append_unique_stream(used_streams, src->data->Bg->stream);
                }
            }
            if (d_data && dir == UNPACK && !same_level_source_fields.empty()) {
                unpack_same_level_fields(src->data, dst->data, same_level_source_fields, same_level_target_fields);
                append_unique_stream(used_streams, dst->data->Bg->stream);
            }
        }
        dst = dst->next;
        src = src->next;
    }

    if (d_data && dir == PACK && (transfer_mode == 2 || transfer_mode == 3) && !stencil_batches.empty()) {
        launch_stencil_transfer_batches(stencil_batches, transfer_mode);
        for (std::vector<StencilTransferBatch>::const_iterator batch = stencil_batches.begin(); batch != stencil_batches.end(); ++batch) {
            if (!batch->transfers.empty()) {
                append_unique_stream(used_streams, batch->stream);
            }
        }
    }

    if (d_data && dir == UNPACK && !unpack_batches.empty()) {
        launch_unpack_transfer_batches(unpack_batches);
        for (std::vector<UnpackTransferBatch>::const_iterator batch = unpack_batches.begin(); batch != unpack_batches.end(); ++batch) {
            if (!batch->transfers.empty()) {
                append_unique_stream(used_streams, batch->stream);
            }
        }
    }

    for (std::vector<cudaStream_t>::const_iterator stream = used_streams.begin(); stream != used_streams.end(); ++stream) {
        GpuRuntimeCoordinator::wait_for_stream(*stream);
    }

    return size_out;
}

#if MPI_CUDA_AWARE
void Parallel::gpu_transfer(
    MyList<gridseg> **src, MyList<gridseg> **dst,
    MyList<var> *VarList1 /* source */, MyList<var> *VarList2 /*target */,
    int Symmetry
) {
    int myrank = 0;
    int cpusize = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &cpusize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (cpusize == 1) {
        gpu_transfer_single_rank_local(src, dst, VarList1, VarList2, Symmetry, myrank);
        return;
    }

    if (!amss_mpi_cuda_aware_runtime_enabled()) {
        if (gpu_transfer_cuda_ipc(src, dst, VarList1, VarList2, Symmetry)) {
            return;
        }
        gpu_transfer_host_staged(src, dst, VarList1, VarList2, Symmetry);
        return;
    }

    int node;

    MPI_Request *reqs;
    MPI_Status *stats;
    reqs = new MPI_Request[2 * cpusize];
    stats = new MPI_Status[2 * cpusize];
    int req_no = 0;

    double **send_data, **rec_data;
    send_data = new double *[cpusize];
    rec_data = new double *[cpusize];
    int length;
    for (node = 0; node < cpusize; node++) {
        send_data[node] = rec_data[node] = 0;
        if (node == myrank) {
            length = gpu_data_packer(0, src[myrank], dst[myrank], node, PACK, VarList1, VarList2, Symmetry);
            if (length) {
                rec_data[node] = GpuRuntimeCoordinator::instance().allocate_device_buffer<double>(length, 0, false);
                gpu_data_packer(rec_data[node], src[myrank], dst[myrank], node, PACK, VarList1, VarList2, Symmetry);
            }
        }
        else {
            length = gpu_data_packer(0, src[myrank], dst[myrank], node, PACK, VarList1, VarList2, Symmetry);
            if (length) {
                send_data[node] = GpuRuntimeCoordinator::instance().allocate_device_buffer<double>(length, 0, false);
                gpu_data_packer(send_data[node], src[myrank], dst[myrank], node, PACK, VarList1, VarList2, Symmetry);
                MPI_Isend((void *)send_data[node], length, MPI_DOUBLE, node, 1, MPI_COMM_WORLD, reqs + req_no ++);
            }
            length = gpu_data_packer(0, src[node], dst[node], node, UNPACK, VarList1, VarList2, Symmetry);
            if (length) {
                rec_data[node] = GpuRuntimeCoordinator::instance().allocate_device_buffer<double>(length, 0, false);
                MPI_Irecv((void *)rec_data[node], length, MPI_DOUBLE, node, 1, MPI_COMM_WORLD, reqs + req_no ++);
            }
        }
    }
    MPI_Waitall(req_no, reqs, stats);
    for (node = 0; node < cpusize; node++) {
        if (rec_data[node]) gpu_data_packer(rec_data[node], src[node], dst[node], node, UNPACK, VarList1, VarList2, Symmetry);
    }
        
    for (node = 0; node < cpusize; node++) {
        if (send_data[node]) GpuRuntimeCoordinator::instance().release_device_buffer(send_data[node], 0);
        if (rec_data[node]) GpuRuntimeCoordinator::instance().release_device_buffer(rec_data[node], 0);
    }

    delete[] reqs;
    delete[] stats;
    delete[] send_data;
    delete[] rec_data;
}
#else
void Parallel::gpu_transfer(
    MyList<gridseg> **src, MyList<gridseg> **dst,
    MyList<var> *VarList1 /* source */, MyList<var> *VarList2 /*target */,
    int Symmetry
) {
    int myrank = 0;
    int cpusize = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &cpusize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (cpusize == 1) {
        gpu_transfer_single_rank_local(src, dst, VarList1, VarList2, Symmetry, myrank);
        return;
    }

    gpu_transfer_host_staged(src, dst, VarList1, VarList2, Symmetry);
}
#endif

void Parallel::Sync_GPU(Patch *Pat, MyList<var> *VarList, int Symmetry) {
    int cpusize;
    MPI_Comm_size(MPI_COMM_WORLD, &cpusize);

    MyList<Parallel::gridseg> *dst;
    MyList<Parallel::gridseg> **src, **transfer_src, **transfer_dst;
    src = new MyList<Parallel::gridseg> *[cpusize];
    transfer_src = new MyList<Parallel::gridseg> *[cpusize];
    transfer_dst = new MyList<Parallel::gridseg> *[cpusize];

    dst = build_ghost_gsl(Pat); // ghost region only
    for (int node = 0; node < cpusize; node++) {
        src[node] = build_owned_gsl0(Pat, node);                              // for the part without ghost points and do not extend
        build_gstl(src[node], dst, &transfer_src[node], &transfer_dst[node]); // for transfer_src[node], data locate on cpu#node;                                                                                                                                            // but for transfer_dst[node] the data may locate on any node
    }

    gpu_transfer(transfer_src, transfer_dst, VarList, VarList, Symmetry);

    if (dst) dst->destroyList();
    for (int node = 0; node < cpusize; node++) {
        if (src[node]) src[node]->destroyList();
        if (transfer_src[node]) transfer_src[node]->destroyList();
        if (transfer_dst[node]) transfer_dst[node]->destroyList();
    }

    delete[] src;
    delete[] transfer_src;
    delete[] transfer_dst;

}

void Parallel::Sync_GPU(MyList<Patch> *PatL, MyList<var> *VarList, int Symmetry) {
    // Patch inner Synch
    MyList<Patch> *Pp = PatL;
    while (Pp) {
        Sync_GPU(Pp->data, VarList, Symmetry);
        Pp = Pp->next;
    }

    // Patch inter Synch
    int cpusize;
    MPI_Comm_size(MPI_COMM_WORLD, &cpusize);

    MyList<Parallel::gridseg> *dst;
    MyList<Parallel::gridseg> **src, **transfer_src, **transfer_dst;
    src = new MyList<Parallel::gridseg> *[cpusize];
    transfer_src = new MyList<Parallel::gridseg> *[cpusize];
    transfer_dst = new MyList<Parallel::gridseg> *[cpusize];

    dst = build_buffer_gsl(PatL); // buffer region only
    for (int node = 0; node < cpusize; node++) {
        src[node] = build_owned_gsl(PatL, node, 5, Symmetry);                 // for the part without ghost nor buffer points and do not extend
        build_gstl(src[node], dst, &transfer_src[node], &transfer_dst[node]); // for transfer[node], data locate on cpu#node
    }

    gpu_transfer(transfer_src, transfer_dst, VarList, VarList, Symmetry);

    if (dst) dst->destroyList();
    for (int node = 0; node < cpusize; node++) {
        if (src[node]) src[node]->destroyList();
        if (transfer_src[node]) transfer_src[node]->destroyList();
        if (transfer_dst[node]) transfer_dst[node]->destroyList();
    }

    delete[] src;
    delete[] transfer_src;
    delete[] transfer_dst;
}

void Parallel::Restrict_GPU(
    MyList<Patch> *PatcL, MyList<Patch> *PatfL,
    MyList<var> *VarList1 /* source */, MyList<var> *VarList2 /* target */,
    int Symmetry
) {
    if (PatcL->data->lev >= PatfL->data->lev)
    {
        cout << "Parallel::Restrict: meet requst of Restrict from lev#" << PatfL->data->lev << " to lev#" << PatcL->data->lev << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    int cpusize;
    MPI_Comm_size(MPI_COMM_WORLD, &cpusize);

    MyList<Parallel::gridseg> *dst;
    MyList<Parallel::gridseg> **src, **transfer_src, **transfer_dst;
    src = new MyList<Parallel::gridseg> *[cpusize];
    transfer_src = new MyList<Parallel::gridseg> *[cpusize];
    transfer_dst = new MyList<Parallel::gridseg> *[cpusize];

    dst = build_complete_gsl(PatcL); // including ghost
    for (int node = 0; node < cpusize; node++)
    {
        // it seems bam always use this
        src[node] = build_owned_gsl(PatfL, node, 2, Symmetry); // - buffer - ghost
        build_gstl(src[node], dst, &transfer_src[node], &transfer_dst[node]); // for transfer[node], data locate on cpu#node
    }

    gpu_transfer(transfer_src, transfer_dst, VarList1, VarList2, Symmetry);

    if (dst)
        dst->destroyList();
    for (int node = 0; node < cpusize; node++)
    {
        if (src[node])
            src[node]->destroyList();
        if (transfer_src[node])
            transfer_src[node]->destroyList();
        if (transfer_dst[node])
            transfer_dst[node]->destroyList();
    }

    delete[] src;
    delete[] transfer_src;
    delete[] transfer_dst;
}

void Parallel::OutBdLow2Hi_GPU(
    Patch *Patc, Patch *Patf,
    MyList<var> *VarList1 /* source */, MyList<var> *VarList2 /* target */,
    int Symmetry
) {
    if (Patc->lev >= Patf->lev)
    {
        cout << "Parallel::OutBdLow2Hi: meet requst of Prolong from lev#" << Patc->lev << " to lev#" << Patf->lev << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    int cpusize;
    MPI_Comm_size(MPI_COMM_WORLD, &cpusize);

    MyList<Parallel::gridseg> *dst;
    MyList<Parallel::gridseg> **src, **transfer_src, **transfer_dst;
    src = new MyList<Parallel::gridseg> *[cpusize];
    transfer_src = new MyList<Parallel::gridseg> *[cpusize];
    transfer_dst = new MyList<Parallel::gridseg> *[cpusize];

    dst = build_buffer_gsl(Patf); // buffer region only

    for (int node = 0; node < cpusize; node++)
    {
        src[node] = build_owned_gsl4(Patc, node, Symmetry);                   // - buffer - ghost - BD ghost
        build_gstl(src[node], dst, &transfer_src[node], &transfer_dst[node]); // for transfer[node], data locate on cpu#node
    }

    gpu_transfer(transfer_src, transfer_dst, VarList1, VarList2, Symmetry);

    if (dst)
        dst->destroyList();
    for (int node = 0; node < cpusize; node++)
    {
        if (src[node])
            src[node]->destroyList();
        if (transfer_src[node])
            transfer_src[node]->destroyList();
        if (transfer_dst[node])
            transfer_dst[node]->destroyList();
    }

    delete[] src;
    delete[] transfer_src;
    delete[] transfer_dst;
}

void Parallel::gpu_prepare_inter_time_level(
    Patch *Pat,
    MyList<var> *VarList1 /* source (t+dt) */, MyList<var> *VarList2 /* source (t) */,
    MyList<var> *VarList3 /* target (t+a*dt) */, int tindex
) {
    int myrank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    MyList<var> *varl1;
    MyList<var> *varl2;
    MyList<var> *varl3;

    MyList<Block> *BP = Pat->blb;
    while (BP)
    {
        Block *cg = BP->data;
        if (myrank == cg->rank)
        {
            varl1 = VarList1;
            varl2 = VarList2;
            varl3 = VarList3;
	            while (varl1)
	            {
	                if (tindex == 0) gpu_average_launch(cg->stream, cg->shape, cg->d_fgfs[varl1->data->sgfn], cg->d_fgfs[varl2->data->sgfn], cg->d_fgfs[varl3->data->sgfn]);
	                else if (tindex == 1) gpu_average3_launch(cg->stream, cg->shape, cg->d_fgfs[varl1->data->sgfn], cg->d_fgfs[varl2->data->sgfn], cg->d_fgfs[varl3->data->sgfn]);
	                else if (tindex == -1) gpu_average3_launch(cg->stream, cg->shape, cg->d_fgfs[varl2->data->sgfn], cg->d_fgfs[varl1->data->sgfn], cg->d_fgfs[varl3->data->sgfn]);
                else {
	                    cout << "error tindex in Parallel::prepare_inter_time_level" << endl;
	                    MPI_Abort(MPI_COMM_WORLD, 1);
	                }
	                cg->mark_gpu_modified(varl3->data->sgfn);
	                varl1 = varl1->next;
	                varl2 = varl2->next;
	                varl3 = varl3->next;
	            }
        }
        if (BP == Pat->ble)
            break;
        BP = BP->next;
    }
    GpuRuntimeCoordinator::instance().wait_for_device();
}

void Parallel::gpu_prepare_inter_time_level(
    Patch *Pat,
    MyList<var> *VarList1 /* source (t+dt) */, MyList<var> *VarList2 /* source (t) */,
    MyList<var> *VarList3 /* source (t-dt) */, MyList<var> *VarList4 /* target (t+a*dt) */, int tindex
) {
    int myrank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    MyList<var> *varl1;
    MyList<var> *varl2;
    MyList<var> *varl3;
    MyList<var> *varl4;

    MyList<Block> *BP = Pat->blb;
    while (BP)
    {
        Block *cg = BP->data;
        if (myrank == cg->rank)
        {
            varl1 = VarList1;
            varl2 = VarList2;
            varl3 = VarList3;
            varl4 = VarList4;
            while (varl1)
            {
                if (tindex == 0) gpu_average2_launch(
                    cg->stream, cg->shape, cg->d_fgfs[varl1->data->sgfn], cg->d_fgfs[varl2->data->sgfn],
                    cg->d_fgfs[varl3->data->sgfn], cg->d_fgfs[varl4->data->sgfn]
                );
                else if (tindex == 1) gpu_average2p_launch(
                    cg->stream, cg->shape, cg->d_fgfs[varl1->data->sgfn], cg->d_fgfs[varl2->data->sgfn],
                    cg->d_fgfs[varl3->data->sgfn], cg->d_fgfs[varl4->data->sgfn]
                );
                else if (tindex == -1) gpu_average2m_launch(
                    cg->stream, cg->shape, cg->d_fgfs[varl1->data->sgfn], cg->d_fgfs[varl2->data->sgfn],
                    cg->d_fgfs[varl3->data->sgfn], cg->d_fgfs[varl4->data->sgfn]
                );
	                else {
	                    cout << "error tindex in long cgh::prepare_inter_time_level" << endl;
	                    MPI_Abort(MPI_COMM_WORLD, 1);
	                }
	                cg->mark_gpu_modified(varl4->data->sgfn);
	                varl1 = varl1->next;
	                varl2 = varl2->next;
	                varl3 = varl3->next;
                varl4 = varl4->next;
            }
        }
        if (BP == Pat->ble)
            break;
        BP = BP->next;
    }
    GpuRuntimeCoordinator::instance().wait_for_device();
}

bool Parallel::PatList_Interp_Points_GPU(
    cudaStream_t stream,
    MyList<Patch> *PatL, MyList<var> *VarList,
    int NN, double *d_XX[3], 
    double *d_Shellf,
    int Symmetry
) {
    return Parallel::PatList_Interp_Points_GPU(stream, PatL, VarList, NN, d_XX, d_Shellf, 0, 0, Symmetry);
}

bool Parallel::PatList_Interp_Points_GPU(
    cudaStream_t stream,
    MyList<Patch> *PatL, MyList<var> *VarList,
    int NN, double *d_XX[3],
    double *d_Shellf,
    double *h_Shellf,
    int *h_PointWeight,
    int Symmetry
) {
    int myrank;
    int cpusize = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &cpusize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    MyList<var> *varl = VarList;
    int num_var = 0;
    while (varl) {
        num_var++;
        varl = varl->next;
    }

    InterpScratch& scratch = interp_scratch();
    ensure_interp_scratch_capacity(scratch, NN, num_var, stream);

    double *d_local_shellf = scratch.d_local_shellf;
    int *d_local_weight = scratch.d_local_weight;

    CUDA_CHECK(cudaMemsetAsync(d_local_shellf, 0, NN * num_var * sizeof(double), stream));
    CUDA_CHECK(cudaMemsetAsync(d_local_weight, 0, NN * sizeof(int), stream));

    int* d_owner = nullptr;
    CUDA_CHECK(cudaMalloc(&d_owner, NN * sizeof(int)));
    CUDA_CHECK(cudaMemsetAsync(d_owner, 0xff, NN * sizeof(int), stream));

    int ordn = 2 * ghost_width;
    int block_tag = 0;

    if (PatL && PatL->data) {
        MyList<Patch> *PL = PatL;
        while (PL) {
            Patch *patch = PL->data;

            double DH[3] = {0.0, 0.0, 0.0};
            for (int i = 0; i < dim; i++) DH[i] = patch->getdX(i);

            MyList<Block> *Bp = patch->blb;
            while (Bp) {
                Block *BP = Bp->data;
                if (myrank == BP->rank) {
                    double llb[3] = {0}, uub[3] = {0};
                    for (int i = 0; i < dim; i++) {
#ifdef Vertex
                        const double block_trim = (static_cast<double>(ghost_width) - 0.5) * DH[i];
#else
                        const double block_trim = static_cast<double>(ghost_width) * DH[i];
#endif
                        const double block_lower = (feq(BP->bbox[i], patch->bbox[i], DH[i] / 2)) ? BP->bbox[i] + patch->lli[i] * DH[i] : BP->bbox[i] + block_trim;
                        const double block_upper = (feq(BP->bbox[dim + i], patch->bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - patch->uui[i] * DH[i] : BP->bbox[dim + i] - block_trim;
                        const double patch_lower = patch->bbox[i] + patch->lli[i] * DH[i];
                        const double patch_upper = patch->bbox[dim + i] - patch->uui[i] * DH[i];
                        llb[i] = std::max(block_lower - 0.5 * DH[i], patch_lower);
                        uub[i] = std::min(block_upper + 0.5 * DH[i], patch_upper);
                    }

                    int shape_0 = BP->shape[0];
                    int shape_1 = BP->shape[1];
                    int shape_2 = BP->shape[2];

                    int my_tag = ++block_tag;

                    varl = VarList;
                    int k = 0;
                    while (varl) {
                        gpu_global_interp_first_match_launch(
                            stream,
                            NN, dim,
                            d_XX[0], d_XX[1], d_XX[2],
                            shape_0, shape_1, shape_2,
                            BP->d_X[0], BP->d_X[1], BP->d_X[2],
                            BP->d_fgfs[varl->data->sgfn],
                            llb[0], llb[1], llb[2],
                            uub[0], uub[1], uub[2],
                            ordn, varl->data->SoA[0], varl->data->SoA[1], varl->data->SoA[2],
                            Symmetry, k, num_var, d_local_shellf, d_local_weight,
                            my_tag, d_owner
                        );
                        varl = varl->next;
                        k ++;
                    }
                }
                Bp = Bp->next;
            }
            PL = PL->next;
        }
    }

    CUDA_CHECK(cudaFree(d_owner));

    double *h_local_shellf = scratch.h_local_shellf;
    int    *h_local_weight = scratch.h_local_weight;
    GpuRuntimeCoordinator::instance().copy_to_host(h_local_shellf, d_local_shellf, NN * num_var, stream);
    GpuRuntimeCoordinator::instance().copy_to_host(h_local_weight, d_local_weight, NN, stream);

    GpuRuntimeCoordinator::wait_for_stream(stream);
    if (cpusize > 1) {
        MPI_Allreduce(MPI_IN_PLACE, h_local_shellf, NN * num_var, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, h_local_weight, NN, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    }

    int h_global_err_flag = 0;
    for (int point_index = 0; point_index < NN; point_index++) {
        const int local_hits = h_local_weight[point_index];
        if (local_hits > 1) {
            const double inv_hits = 1.0 / static_cast<double>(local_hits);
            for (int variable_index = 0; variable_index < num_var; variable_index++) {
                h_local_shellf[variable_index + point_index * num_var] *= inv_hits;
            }
        }
        else if (local_hits == 0) {
            h_global_err_flag = 1;
        }
        bool finite_shellf = (local_hits > 0);
        if (finite_shellf) {
            for (int variable_index = 0; variable_index < num_var; variable_index++) {
                const int slot = variable_index + point_index * num_var;
                if (!std::isfinite(h_local_shellf[slot])) {
                    finite_shellf = false;
                    break;
                }
            }
        }
        if (local_hits == 0 || !finite_shellf) {
            h_global_err_flag = 1;
            h_local_weight[point_index] = 0;
            for (int variable_index = 0; variable_index < num_var; variable_index++) {
                h_local_shellf[variable_index + point_index * num_var] = 0.0;
            }
        } else {
            h_local_weight[point_index] = 1;
        }
    }

    if (h_Shellf) {
        std::memcpy(h_Shellf, h_local_shellf, NN * num_var * sizeof(double));
    }

    if (h_PointWeight) {
        std::memcpy(h_PointWeight, h_local_weight, NN * sizeof(int));
    }

    GpuRuntimeCoordinator::instance().copy_to_device(h_local_shellf, d_Shellf, NN * num_var, stream);

    bool success = (h_global_err_flag == 0);

    if (!success) {
        checkpatchlist(PatL, false);
    }

    GpuRuntimeCoordinator::wait_for_stream(stream);
    return success;
}

void Parallel::gpu_prepare_inter_time_level(
    MyList<Patch> *PatL,
    MyList<var> *VarList1 /* source (t+dt) */, MyList<var> *VarList2 /* source (t) */,
    MyList<var> *VarList3 /* source (t-dt) */, MyList<var> *VarList4 /* target (t+a*dt) */, int tindex
) {
    while (PatL) {
        gpu_prepare_inter_time_level(PatL->data, VarList1, VarList2, VarList3, VarList4, tindex);
        PatL = PatL->next;
    }
}

void Parallel::gpu_prepare_inter_time_level(
    MyList<Patch> *PatL,
    MyList<var> *VarList1 /* source (t+dt) */, MyList<var> *VarList2 /* source (t) */,
    MyList<var> *VarList3 /* target (t+a*dt) */, int tindex
) {
    while (PatL) {
        gpu_prepare_inter_time_level(PatL->data, VarList1, VarList2, VarList3, tindex);
        PatL = PatL->next;
    }
}

void Parallel::gpu_fill_level_data(
    MyList<Patch> *PatLd, MyList<Patch> *PatLs, MyList<Patch> *PatcL,
    MyList<var> *OldList, MyList<var> *StateList, MyList<var> *FutureList,
    MyList<var> *tmList, int Symmetry, bool BB, bool CC
) {
    if (PatLd->data->lev != PatLs->data->lev) {
        cout << "Parallel::fill_level_data: meet requst from lev#" << PatLs->data->lev << " to lev#" << PatLd->data->lev << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (PatLd->data->lev <= PatcL->data->lev) {
        cout << "Parallel::fill_level_data: meet prolong requst from lev#" << PatcL->data->lev << " to lev#" << PatLd->data->lev << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int cpusize;
    MPI_Comm_size(MPI_COMM_WORLD, &cpusize);

    MyList<var> *VarList = 0;
    MyList<var> *p;
    p = StateList;
    while (p) {
        if (VarList) VarList->insert(p->data);
        else VarList = new MyList<var>(p->data);
        p = p->next;
    }
    p = FutureList;
    while (p) {
        if (VarList) VarList->insert(p->data);
        else VarList = new MyList<var>(p->data);
        p = p->next;
    }

    MyList<Parallel::gridseg> *dst;
    MyList<Parallel::gridseg> **src, **transfer_src, **transfer_dst;
    src = new MyList<Parallel::gridseg> *[cpusize];
    transfer_src = new MyList<Parallel::gridseg> *[cpusize];
    transfer_dst = new MyList<Parallel::gridseg> *[cpusize];

    dst = build_complete_gsl(PatLd); // including ghost
    // copy part
    for (int node = 0; node < cpusize; node++) {
        src[node] = build_owned_gsl(PatLs, node, 0, Symmetry);                // similar to Sync
        build_gstl(src[node], dst, &transfer_src[node], &transfer_dst[node]); // for transfer[node], data locate on cpu#node
    }

    gpu_transfer(transfer_src, transfer_dst, VarList, VarList, Symmetry);

    for (int node = 0; node < cpusize; node++) {
        if (src[node]) src[node]->destroyList();
        if (transfer_src[node]) transfer_src[node]->destroyList();
        if (transfer_dst[node]) transfer_dst[node]->destroyList();
    }

    MyList<Parallel::gridseg> *dsts, *dstd;
    dsts = build_complete_gsl_virtual(PatLs);
    dstd = dst;
    dst = gsl_subtract(dstd, dsts);
    if (dstd) dstd->destroyList();
    if (dsts) dsts->destroyList();

    if (dst) {
        // prolongation part
        for (int node = 0; node < cpusize; node++) {
            src[node] = build_owned_gsl(PatcL, node, 4, Symmetry);                // - buffer - ghost - BD ghost
            build_gstl(src[node], dst, &transfer_src[node], &transfer_dst[node]); // for transfer[node], data locate on cpu#node
        }

        if (CC) {
            // for FutureList
            // restrict first~~~>
            {
                Restrict_GPU(PatcL, PatLs, FutureList, FutureList, Symmetry);
                Sync_GPU(PatcL, FutureList, Symmetry);
            }
            //<~~~prolong then
            gpu_transfer(transfer_src, transfer_dst, FutureList, FutureList, Symmetry);

            // for StateList
            // time interpolation part
            if (BB) gpu_prepare_inter_time_level(PatcL, FutureList, StateList, OldList, tmList, 0); // use SynchList_pre as temporal storage space
            else gpu_prepare_inter_time_level(PatcL, FutureList, StateList, tmList, 0); // use SynchList_pre as temporal storage space
            
            {
                Restrict_GPU(PatcL, PatLs, StateList, tmList, Symmetry);
                Sync_GPU(PatcL, tmList, Symmetry);
            }
            //<~~~prolong then
            gpu_transfer(transfer_src, transfer_dst, tmList, StateList, Symmetry);
        }
        else {
            // for both FutureList and StateList
            // restrict first~~~>
            {
                Restrict_GPU(PatcL, PatLs, VarList, VarList, Symmetry);
                Sync_GPU(PatcL, VarList, Symmetry);
            }
            //<~~~prolong then
            gpu_transfer(transfer_src, transfer_dst, VarList, VarList, Symmetry);
        }

        for (int node = 0; node < cpusize; node++) {
            if (src[node])
                src[node]->destroyList();
            if (transfer_src[node])
                transfer_src[node]->destroyList();
            if (transfer_dst[node])
                transfer_dst[node]->destroyList();
        }

        dst->destroyList();
    }

    delete[] src;
    delete[] transfer_src;
    delete[] transfer_dst;

    VarList->clearList();
}

double Parallel::L2Norm_GPU(Patch *Pat, var *vf) {
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    double tvf, dtvf = 0;
    int BDW = ghost_width;

    MyList<Block> *BP = Pat->blb;
    while (BP) {
        Block *cg = BP->data;
        if (myrank == cg->rank) {
            gpu_l2normhelper_launch(
                cg->stream, 
                cg->shape, cg->X[0], cg->X[1], cg->X[2],
                Pat->bbox[0], Pat->bbox[1], Pat->bbox[2],
                Pat->bbox[3], Pat->bbox[4], Pat->bbox[5],
                cg->d_fgfs[vf->sgfn], tvf, BDW
            );
            dtvf += tvf;
        }
        if (BP == Pat->ble)
            break;
        BP = BP->next;
    }

    MPI_Allreduce(&dtvf, &tvf, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    tvf = sqrt(tvf);

    return tvf;
}

void Parallel::Dump_Data_GPU(MyList<Patch> *PL, MyList<var> *DumpList, char *tag, double time, double dT) {
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    DataMotionHelper::pull_all_data_to_host(PL, myrank, DumpList);
    MyList<Patch> *Pp;
    Pp = PL;
    int grd = 0;
    while (Pp)
    {
        Patch *PP = Pp->data;
        Dump_Data(PP, DumpList, tag, time, dT, grd);
        grd++;
        Pp = Pp->next;
    }
}

void Parallel::d2Dump_Data_GPU(MyList<Patch> *PL, MyList<var> *DumpList, char *tag, double time, double dT) {
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    DataMotionHelper::pull_all_data_to_host(PL, myrank, DumpList);
    MyList<Patch> *Pp;
    Pp = PL;
    int grd = 0;
    while (Pp)
    {
        Patch *PP = Pp->data;
        d2Dump_Data(PP, DumpList, tag, time, dT, grd);
        grd++;
        Pp = Pp->next;
    }
}
