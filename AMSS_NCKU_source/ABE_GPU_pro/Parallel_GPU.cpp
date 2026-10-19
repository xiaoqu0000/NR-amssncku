#include "Parallel.h"

#include "gpu_manager.h"
#include "prolongrestrict.h"
#include "misc.h"
#include "parameters.h"
#include "fmisc.h"

#include "gpu_manager.h"
#include "MPatch.h"
#include "helper.h"
#include "gpu_info.h"
#include "reduce_helper.h"

#include <vector>
#include <cstdlib>

int Parallel::gpu_data_packer(
    double *d_data, MyList<Parallel::gridseg> *src, MyList<Parallel::gridseg> *dst, int rank_in, int dir,
    MyList<var> *VarLists, MyList<var> *VarListd, int Symmetry,
    std::vector<cudaStream_t>* out_streams
) { // Fused version
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int DIM = dim;

    if (dir != PACK && dir != UNPACK) {
        cout << "error dir " << dir << " for data_packer " << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int size_out = 0;

    if (!src || !dst) return size_out;

    // If the caller passed out_streams, accumulate the Block streams we launch
    // kernels on (with simple dedup). Caller then host_waits each stream once,
    // avoiding cudaEventRecord on the legacy default stream (which would
    // implicitly synchronize all streams on every call).
    auto track_stream = [&](cudaStream_t kstream) {
        if (!out_streams) return;
        for (auto s : *out_streams) if (s == kstream) return;
        out_streams->push_back(kstream);
    };

    MyList<var> *varls, *varld;

    varls = VarLists;
    varld = VarListd;
    while (varls && varld) {
        varls = varls->next;
        varld = varld->next;
    }

    if (varls || varld) {
        cout << "error in short data packer, var lists does not match." << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int type; /* 1 copy, 2 restrict, 3 prolong */
    if (src->data->Bg->lev == dst->data->Bg->lev) type = 1;
    else if (src->data->Bg->lev > dst->data->Bg->lev) type = 2;
    else type = 3;
    
    while (src && dst) {
        if (
            (dir == PACK && dst->data->Bg->rank == rank_in && src->data->Bg->rank == myrank) ||
            (dir == UNPACK && src->data->Bg->rank == rank_in && dst->data->Bg->rank == myrank)
        ) {
            varls = VarLists;
            varld = VarListd;
            
            // 增加内层分批循环，防止变量总数超过 MAX_FUSED_VARS
            while (varls && varld) {
                FusedVars fused_vars;
                fused_vars.num_vars = 0;
                int fused_dst_sgfn[MAX_FUSED_VARS];

                while (varls && varld && fused_vars.num_vars < MAX_FUSED_VARS) {
                    if (dir == PACK) {
                        fused_vars.ptrs[fused_vars.num_vars] = src->data->Bg->d_fgfs[varls->data->sgfn];
                    } else if (dir == UNPACK) {
                        fused_vars.ptrs[fused_vars.num_vars] = dst->data->Bg->d_fgfs[varld->data->sgfn];
                        fused_dst_sgfn[fused_vars.num_vars] = varld->data->sgfn;
                    }

                    fused_vars.SoA0[fused_vars.num_vars] = varls->data->SoA[0];
                    fused_vars.SoA1[fused_vars.num_vars] = varls->data->SoA[1];
                    fused_vars.SoA2[fused_vars.num_vars] = varls->data->SoA[2];

                    fused_vars.num_vars++;
                    varls = varls->next;
                    varld = varld->next;
                }

                if (d_data) {
                    if (dir == PACK) {
                        double* d_dst_ptr = d_data + size_out; 

                        switch (type) {
                        case 1: { // COPY (PACK)
#ifdef Vertex
                            double dx = (src->data->Bg->bbox[3] - src->data->Bg->bbox[0]) / (src->data->Bg->shape[0] - 1);
                            double dy = (src->data->Bg->bbox[4] - src->data->Bg->bbox[1]) / (src->data->Bg->shape[1] - 1);
                            double dz = (src->data->Bg->bbox[5] - src->data->Bg->bbox[2]) / (src->data->Bg->shape[2] - 1);
#else
                            double dx = (src->data->Bg->bbox[3] - src->data->Bg->bbox[0]) / src->data->Bg->shape[0];
                            double dy = (src->data->Bg->bbox[4] - src->data->Bg->bbox[1]) / src->data->Bg->shape[1];
                            double dz = (src->data->Bg->bbox[5] - src->data->Bg->bbox[2]) / src->data->Bg->shape[2];
#endif

                            int off_x = (int)std::trunc((dst->data->llb[0] - src->data->Bg->bbox[0]) / dx + 0.4);
                            int off_y = (int)std::trunc((dst->data->llb[1] - src->data->Bg->bbox[1]) / dy + 0.4);
                            int off_z = (int)std::trunc((dst->data->llb[2] - src->data->Bg->bbox[2]) / dz + 0.4);

                            gpu_fused_pack_launch(
                                src->data->Bg->stream,
                                fused_vars, d_dst_ptr,
                                src->data->Bg->shape_pad[0], src->data->Bg->shape_pad[1],
                                dst->data->shape[0], dst->data->shape[1], dst->data->shape[2],
                                off_x, off_y, off_z
                            );
                            break;
                        }

                        case 2: { // RESTRICT
                            gpu_apply_symmetry_bd_vars_launch(
                                src->data->Bg->stream,
                                src->data->Bg->shape[0], src->data->Bg->shape[1], src->data->Bg->shape[2],
                                src->data->Bg->shape_pad[0], src->data->Bg->shape_pad[1], src->data->Bg->shape_pad[2],
                                Symmetry,
                                fused_vars
                            );
                            gpu_fused_restrict3_launch(
                                src->data->Bg->stream,
                                fused_vars, d_dst_ptr,
                                dst->data->llb, dst->data->uub, dst->data->shape,
                                src->data->Bg->bbox, src->data->Bg->bbox + dim, src->data->Bg->shape, src->data->Bg->shape_pad,
                                dst->data->llb, dst->data->uub,
                                Symmetry
                            );
                            break;
                        }

                        case 3: { // PROLONG
                            gpu_apply_symmetry_bd_vars_launch(
                                src->data->Bg->stream,
                                src->data->Bg->shape[0], src->data->Bg->shape[1], src->data->Bg->shape[2],
                                src->data->Bg->shape_pad[0], src->data->Bg->shape_pad[1], src->data->Bg->shape_pad[2],
                                Symmetry,
                                fused_vars
                            );
                            gpu_fused_prolong3_launch(
                                src->data->Bg->stream,
                                fused_vars, d_dst_ptr,
                                src->data->Bg->bbox, src->data->Bg->bbox + dim, src->data->Bg->shape, src->data->Bg->shape_pad,
                                dst->data->llb, dst->data->uub, dst->data->shape,
                                dst->data->llb, dst->data->uub,
                                Symmetry
                            );
                            break;
                        }
                        }
                    }
                    
                    if (dir == UNPACK) {
                        double* d_src_ptr = d_data + size_out;

#ifdef Vertex
                        double dx = (dst->data->Bg->bbox[3] - dst->data->Bg->bbox[0]) / (dst->data->Bg->shape[0] - 1);
                        double dy = (dst->data->Bg->bbox[4] - dst->data->Bg->bbox[1]) / (dst->data->Bg->shape[1] - 1);
                        double dz = (dst->data->Bg->bbox[5] - dst->data->Bg->bbox[2]) / (dst->data->Bg->shape[2] - 1);
#else
                        double dx = (dst->data->Bg->bbox[3] - dst->data->Bg->bbox[0]) / dst->data->Bg->shape[0];
                        double dy = (dst->data->Bg->bbox[4] - dst->data->Bg->bbox[1]) / dst->data->Bg->shape[1];
                        double dz = (dst->data->Bg->bbox[5] - dst->data->Bg->bbox[2]) / dst->data->Bg->shape[2];
#endif

                        int off_x = (int)std::trunc((dst->data->llb[0] - dst->data->Bg->bbox[0]) / dx + 0.4);
                        int off_y = (int)std::trunc((dst->data->llb[1] - dst->data->Bg->bbox[1]) / dy + 0.4);
                        int off_z = (int)std::trunc((dst->data->llb[2] - dst->data->Bg->bbox[2]) / dz + 0.4);

                        gpu_fused_unpack_launch(
                            dst->data->Bg->stream,
                            d_src_ptr, fused_vars,
                            dst->data->Bg->shape_pad[0], dst->data->Bg->shape_pad[1],
                            dst->data->shape[0], dst->data->shape[1], dst->data->shape[2],
                            off_x, off_y, off_z
                        );
                        for (int v = 0; v < fused_vars.num_vars; v++) {
                            dst->data->Bg->mark_gpu_modified(fused_dst_sgfn[v]);
                        }
                    }

                    // Record which Block stream this fused batch was launched on so
                    // the caller can host_wait it directly (instead of going through
                    // a default-stream record that would implicitly sync all streams).
                    track_stream((dir == PACK) ? src->data->Bg->stream
                                               : dst->data->Bg->stream);
                }

                // 统一累加 size_out，保证下一次循环的起始指针完全正确
                int single_var_size = dst->data->shape[0] * dst->data->shape[1] * dst->data->shape[2];
                size_out += single_var_size * fused_vars.num_vars;
            }
        }
        dst = dst->next;
        src = src->next;
    }
    return size_out;
}

#if MPI_CUDA_AWARE
void Parallel::gpu_transfer(
    MyList<gridseg> **src, MyList<gridseg> **dst,
    MyList<var> *VarList1 /* source */, MyList<var> *VarList2 /*target */,
    int Symmetry
) {
    int myrank, cpusize;
    MPI_Comm_size(MPI_COMM_WORLD, &cpusize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    double **send_data_d = new double *[cpusize];
    double **rec_data_d  = new double *[cpusize];
    int *send_len = new int[cpusize];
    int *recv_len = new int[cpusize];

    // Phase 1: Size calc + allocate
    for (int node = 0; node < cpusize; node++) {
        send_data_d[node] = rec_data_d[node] = nullptr;
        send_len[node] = recv_len[node] = 0;

        if (node == myrank) {
            recv_len[node] = gpu_data_packer(0, src[myrank], dst[myrank], node, PACK, VarList1, VarList2, Symmetry);
            if (recv_len[node] > 0)
                rec_data_d[node] = GPUManager::getInstance().allocate_comm_buffer<double>(recv_len[node]);
        } else {
            send_len[node] = gpu_data_packer(0, src[myrank], dst[myrank], node, PACK, VarList1, VarList2, Symmetry);
            if (send_len[node] > 0)
                send_data_d[node] = GPUManager::getInstance().allocate_comm_buffer<double>(send_len[node]);

            recv_len[node] = gpu_data_packer(0, src[node], dst[node], node, UNPACK, VarList1, VarList2, Symmetry);
            if (recv_len[node] > 0)
                rec_data_d[node] = GPUManager::getInstance().allocate_comm_buffer<double>(recv_len[node]);
        }
    }

    // Phase 1.5: Ensure all async memsets complete on memory_stream before PACK
    // kernels start. The device-side fence (fence_pool_after_memory) introduced
    // by be0769a caused an RMS regression — likely event-pool re-use letting a
    // subsequent cudaEventRecord overwrite the wait edge before it was honoured.
    // Host-blocking sync is only ~1us and the perf delta vs the fence variant
    // is in the noise floor.
    GPUManager::getInstance().synchronize_memory();

    // Phase 2: PACK on each src Block stream. Track which Block streams were
    // touched so we can host_wait only those streams below — not the whole device.
    std::vector<cudaStream_t> pack_streams;
    pack_streams.reserve(8);
    for (int node = 0; node < cpusize; node++) {
        if (node == myrank) {
            if (recv_len[node] > 0)
                gpu_data_packer(rec_data_d[node], src[myrank], dst[myrank], node, PACK,
                                VarList1, VarList2, Symmetry, &pack_streams);
        } else {
            if (send_len[node] > 0)
                gpu_data_packer(send_data_d[node], src[myrank], dst[myrank], node, PACK,
                                VarList1, VarList2, Symmetry, &pack_streams);
        }
    }

    // Phase 3: host-wait every Block stream actually used by PACK. This replaces
    // the old cudaEventRecord(.,0) per-peer pattern (which silently triggered a
    // full implicit device sync on every record) with a targeted per-stream sync.
    for (auto s : pack_streams) cudaStreamSynchronize(s);

    // Phase 4: MPI communication — Irecv first to let UCX rendezvous early.
    MPI_Request *reqs = new MPI_Request[2 * cpusize];
    MPI_Status *stats = new MPI_Status[2 * cpusize];
    int req_no = 0;

    for (int node = 0; node < cpusize; node++) {
        if (node != myrank && recv_len[node] > 0)
            MPI_Irecv((void *)rec_data_d[node], recv_len[node], MPI_DOUBLE, node, 1, MPI_COMM_WORLD, reqs + req_no++);
    }
    for (int node = 0; node < cpusize; node++) {
        if (node != myrank && send_len[node] > 0)
            MPI_Isend((void *)send_data_d[node], send_len[node], MPI_DOUBLE, node, 1, MPI_COMM_WORLD, reqs + req_no++);
    }
    MPI_Waitall(req_no, reqs, stats);

    // Phase 5: All UNPACK on dst Block streams. Track them for the targeted sync.
    std::vector<cudaStream_t> unpack_streams;
    unpack_streams.reserve(8);
    for (int node = 0; node < cpusize; node++) {
        if (rec_data_d[node])
            gpu_data_packer(rec_data_d[node], src[node], dst[node], node, UNPACK,
                            VarList1, VarList2, Symmetry, &unpack_streams);
    }

    // Phase 6: host-wait every Block stream UNPACK was launched on, then free
    // buffers. Targeted per-stream sync — no full-device barrier.
    for (auto s : unpack_streams) cudaStreamSynchronize(s);

    for (int node = 0; node < cpusize; node++) {
        if (send_data_d[node]) GPUManager::getInstance().free_comm_buffer(send_data_d[node]);
        if (rec_data_d[node]) GPUManager::getInstance().free_comm_buffer(rec_data_d[node]);
    }

    delete[] reqs;
    delete[] stats;
    delete[] send_data_d;
    delete[] rec_data_d;
    delete[] send_len;
    delete[] recv_len;
}

// Pipelined two-phase transfer. See header for rationale.
void Parallel::gpu_transfer_pipelined(
    MyList<gridseg> **src_a, MyList<gridseg> **dst_a,
    MyList<gridseg> **src_b, MyList<gridseg> **dst_b,
    MyList<var> *VarList, int Symmetry
) {
    int myrank, cpusize;
    MPI_Comm_size(MPI_COMM_WORLD, &cpusize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    double **send_d_a = new double *[cpusize];
    double **recv_d_a = new double *[cpusize];
    double **send_d_b = new double *[cpusize];
    double **recv_d_b = new double *[cpusize];
    int *send_len_a = new int[cpusize];
    int *recv_len_a = new int[cpusize];
    int *send_len_b = new int[cpusize];
    int *recv_len_b = new int[cpusize];

    auto size_and_alloc = [&](MyList<gridseg> **src, MyList<gridseg> **dst,
                              double **send_d, double **recv_d,
                              int *send_len, int *recv_len) {
        for (int node = 0; node < cpusize; node++) {
            send_d[node] = recv_d[node] = nullptr;
            send_len[node] = recv_len[node] = 0;
            if (node == myrank) {
                recv_len[node] = gpu_data_packer(0, src[myrank], dst[myrank], node, PACK, VarList, VarList, Symmetry);
                if (recv_len[node] > 0)
                    recv_d[node] = GPUManager::getInstance().allocate_comm_buffer<double>(recv_len[node]);
            } else {
                send_len[node] = gpu_data_packer(0, src[myrank], dst[myrank], node, PACK, VarList, VarList, Symmetry);
                if (send_len[node] > 0)
                    send_d[node] = GPUManager::getInstance().allocate_comm_buffer<double>(send_len[node]);
                recv_len[node] = gpu_data_packer(0, src[node], dst[node], node, UNPACK, VarList, VarList, Symmetry);
                if (recv_len[node] > 0)
                    recv_d[node] = GPUManager::getInstance().allocate_comm_buffer<double>(recv_len[node]);
            }
        }
    };

    // Phase 1: size calc + comm-pool alloc for both A and B (host-only)
    size_and_alloc(src_a, dst_a, send_d_a, recv_d_a, send_len_a, recv_len_a);
    size_and_alloc(src_b, dst_b, send_d_b, recv_d_b, send_len_b, recv_len_b);

    // Phase 1.5: host-wait the pool's cudaMemsetAsync ops on memory_stream.
    // Without this, the early Irecv posted in Phase 2 races with the pending
    // memset: a non-GDR transport (cuda_copy) will cudaMemcpyAsync the incoming
    // payload into recv_d_*, and the trailing memset can wipe it. Legacy path
    // does not race because its Irecv happens after the PACK host_wait, which
    // transitively waits for the memset via the Block-stream wait event.
    GPUManager::getInstance().synchronize_memory();

    // Phase 2: post Irecv A (tag=1) and Irecv B (tag=101) EARLY — UCX rendezvous
    // for B starts during A's PACK/transfer. All recvs must be posted before any
    // sends to avoid eager/rendezvous deadlock when both phases are in flight.
    MPI_Request *reqs = new MPI_Request[4 * cpusize];
    MPI_Status *stats = new MPI_Status[4 * cpusize];
    int req_no = 0;
    for (int node = 0; node < cpusize; node++) {
        if (node != myrank && recv_len_a[node] > 0)
            MPI_Irecv((void*)recv_d_a[node], recv_len_a[node], MPI_DOUBLE, node, 1, MPI_COMM_WORLD, reqs + req_no++);
    }
    for (int node = 0; node < cpusize; node++) {
        if (node != myrank && recv_len_b[node] > 0)
            MPI_Irecv((void*)recv_d_b[node], recv_len_b[node], MPI_DOUBLE, node, 101, MPI_COMM_WORLD, reqs + req_no++);
    }

    // Phase 3: host-block until memory_stream is drained, ensuring the pool's
    // cudaMemsetAsync ops are visible to subsequent PACK kernels on Block streams.
    // Replaces the device-side fence_pool_after_memory variant, which races via
    // event-pool re-use (see be0769a fix in non-pipelined gpu_transfer).
    GPUManager::getInstance().synchronize_memory();

    // Phase 4: PACK A then PACK B; sharing pack_streams gives free cross-phase dedup.
    std::vector<cudaStream_t> pack_streams;
    pack_streams.reserve(16);
    for (int node = 0; node < cpusize; node++) {
        if (node == myrank) {
            if (recv_len_a[node] > 0)
                gpu_data_packer(recv_d_a[node], src_a[myrank], dst_a[myrank], node, PACK,
                                VarList, VarList, Symmetry, &pack_streams);
        } else {
            if (send_len_a[node] > 0)
                gpu_data_packer(send_d_a[node], src_a[myrank], dst_a[myrank], node, PACK,
                                VarList, VarList, Symmetry, &pack_streams);
        }
    }
    for (int node = 0; node < cpusize; node++) {
        if (node == myrank) {
            if (recv_len_b[node] > 0)
                gpu_data_packer(recv_d_b[node], src_b[myrank], dst_b[myrank], node, PACK,
                                VarList, VarList, Symmetry, &pack_streams);
        } else {
            if (send_len_b[node] > 0)
                gpu_data_packer(send_d_b[node], src_b[myrank], dst_b[myrank], node, PACK,
                                VarList, VarList, Symmetry, &pack_streams);
        }
    }

    // Phase 5: host_wait merged pack stream set, then post Isend A + Isend B
    for (auto s : pack_streams) cudaStreamSynchronize(s);

    for (int node = 0; node < cpusize; node++) {
        if (node != myrank && send_len_a[node] > 0)
            MPI_Isend((void*)send_d_a[node], send_len_a[node], MPI_DOUBLE, node, 1, MPI_COMM_WORLD, reqs + req_no++);
    }
    for (int node = 0; node < cpusize; node++) {
        if (node != myrank && send_len_b[node] > 0)
            MPI_Isend((void*)send_d_b[node], send_len_b[node], MPI_DOUBLE, node, 101, MPI_COMM_WORLD, reqs + req_no++);
    }

    // Phase 6: single Waitall for all A + B requests
    MPI_Waitall(req_no, reqs, stats);

    // Phase 7: UNPACK A then UNPACK B
    std::vector<cudaStream_t> unpack_streams;
    unpack_streams.reserve(16);
    for (int node = 0; node < cpusize; node++) {
        if (recv_d_a[node])
            gpu_data_packer(recv_d_a[node], src_a[node], dst_a[node], node, UNPACK,
                            VarList, VarList, Symmetry, &unpack_streams);
    }
    for (int node = 0; node < cpusize; node++) {
        if (recv_d_b[node])
            gpu_data_packer(recv_d_b[node], src_b[node], dst_b[node], node, UNPACK,
                            VarList, VarList, Symmetry, &unpack_streams);
    }

    // Phase 8: host_wait merged unpack streams, free all comm buffers
    for (auto s : unpack_streams) cudaStreamSynchronize(s);

    for (int node = 0; node < cpusize; node++) {
        if (send_d_a[node]) GPUManager::getInstance().free_comm_buffer(send_d_a[node]);
        if (recv_d_a[node]) GPUManager::getInstance().free_comm_buffer(recv_d_a[node]);
        if (send_d_b[node]) GPUManager::getInstance().free_comm_buffer(send_d_b[node]);
        if (recv_d_b[node]) GPUManager::getInstance().free_comm_buffer(recv_d_b[node]);
    }

    delete[] reqs;
    delete[] stats;
    delete[] send_d_a;
    delete[] recv_d_a;
    delete[] send_d_b;
    delete[] recv_d_b;
    delete[] send_len_a;
    delete[] recv_len_a;
    delete[] send_len_b;
    delete[] recv_len_b;
}
#else
void Parallel::gpu_transfer(
    MyList<gridseg> **src, MyList<gridseg> **dst,
    MyList<var> *VarList1 /* source */, MyList<var> *VarList2 /*target */,
    int Symmetry
) {
    int myrank, cpusize;
    MPI_Comm_size(MPI_COMM_WORLD, &cpusize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    double **send_data_d = new double *[cpusize];
    double **rec_data_d  = new double *[cpusize];
    double **send_data_h = new double *[cpusize];
    double **rec_data_h  = new double *[cpusize];
    int *send_len = new int[cpusize];
    int *recv_len = new int[cpusize];

    // Phase 1: Size calc + allocate (no sync)
    for (int node = 0; node < cpusize; node++) {
        send_data_d[node] = rec_data_d[node] = nullptr;
        send_data_h[node] = rec_data_h[node] = nullptr;
        send_len[node] = recv_len[node] = 0;

        if (node == myrank) {
            recv_len[node] = gpu_data_packer(0, src[myrank], dst[myrank], node, PACK, VarList1, VarList2, Symmetry);
            if (recv_len[node] > 0)
                rec_data_d[node] = GPUManager::getInstance().allocate_device_memory<double>(recv_len[node]);
        } else {
            send_len[node] = gpu_data_packer(0, src[myrank], dst[myrank], node, PACK, VarList1, VarList2, Symmetry);
            if (send_len[node] > 0) {
                send_data_d[node] = GPUManager::getInstance().allocate_device_memory<double>(send_len[node]);
                send_data_h[node] = new double[send_len[node]];
            }
            recv_len[node] = gpu_data_packer(0, src[node], dst[node], node, UNPACK, VarList1, VarList2, Symmetry);
            if (recv_len[node] > 0) {
                rec_data_d[node] = GPUManager::getInstance().allocate_device_memory<double>(recv_len[node]);
                rec_data_h[node] = new double[recv_len[node]];
            }
        }
    }

    // Phase 1.5: Ensure all async allocations complete on memory_stream before PACK kernels start
    GPUManager::getInstance().synchronize_memory();

    // Phase 2: PACK on each src Block stream, tracking which streams were touched.
    std::vector<cudaStream_t> pack_streams;
    pack_streams.reserve(8);
    for (int node = 0; node < cpusize; node++) {
        if (node == myrank) {
            if (recv_len[node] > 0)
                gpu_data_packer(rec_data_d[node], src[myrank], dst[myrank], node, PACK,
                                VarList1, VarList2, Symmetry, &pack_streams);
        } else {
            if (send_len[node] > 0)
                gpu_data_packer(send_data_d[node], src[myrank], dst[myrank], node, PACK,
                                VarList1, VarList2, Symmetry, &pack_streams);
        }
    }

    // Phase 3: targeted per-stream sync, then D2H copies.
    for (auto s : pack_streams) cudaStreamSynchronize(s);
    for (int node = 0; node < cpusize; node++) {
        if (node != myrank && send_len[node] > 0)
            GPUManager::getInstance().sync_to_cpu(send_data_h[node], send_data_d[node], send_len[node]);
    }
    GPUManager::getInstance().synchronize_memory();

    // Phase 4: MPI — post Irecv first, then Isend
    MPI_Request *reqs = new MPI_Request[2 * cpusize];
    MPI_Status *stats = new MPI_Status[2 * cpusize];
    int req_no = 0;
    for (int node = 0; node < cpusize; node++) {
        if (node != myrank && recv_len[node] > 0)
            MPI_Irecv((void *)rec_data_h[node], recv_len[node], MPI_DOUBLE, node, 1, MPI_COMM_WORLD, reqs + req_no++);
    }
    for (int node = 0; node < cpusize; node++) {
        if (node != myrank && send_len[node] > 0)
            MPI_Isend((void *)send_data_h[node], send_len[node], MPI_DOUBLE, node, 1, MPI_COMM_WORLD, reqs + req_no++);
    }
    MPI_Waitall(req_no, reqs, stats);

    // Phase 5: H2D copies for UNPACK inputs
    for (int node = 0; node < cpusize; node++) {
        if (node != myrank && recv_len[node] > 0)
            GPUManager::getInstance().sync_to_gpu(rec_data_h[node], rec_data_d[node], recv_len[node]);
    }
    GPUManager::getInstance().synchronize_memory();

    // Phase 6: All UNPACK kernels on dst Block streams, tracking streams used.
    std::vector<cudaStream_t> unpack_streams;
    unpack_streams.reserve(8);
    for (int node = 0; node < cpusize; node++) {
        if (node == myrank) {
            if (rec_data_d[node])
                gpu_data_packer(rec_data_d[node], src[node], dst[node], node, UNPACK,
                                VarList1, VarList2, Symmetry, &unpack_streams);
        } else {
            if (rec_data_h[node])
                gpu_data_packer(rec_data_d[node], src[node], dst[node], node, UNPACK,
                                VarList1, VarList2, Symmetry, &unpack_streams);
        }
    }

    // Phase 7: targeted per-stream sync for UNPACK, then free buffers.
    for (auto s : unpack_streams) cudaStreamSynchronize(s);

    for (int node = 0; node < cpusize; node++) {
        if (send_data_d[node]) GPUManager::getInstance().free_device_memory(send_data_d[node]);
        if (rec_data_d[node]) GPUManager::getInstance().free_device_memory(rec_data_d[node]);
        if (send_data_h[node]) delete[] send_data_h[node];
        if (rec_data_h[node]) delete[] rec_data_h[node];
    }

    delete[] reqs;
    delete[] stats;
    delete[] send_data_d;
    delete[] rec_data_d;
    delete[] send_data_h;
    delete[] rec_data_h;
    delete[] send_len;
    delete[] recv_len;
}

// Pipelined two-phase transfer (host-staged variant). See header for rationale.
void Parallel::gpu_transfer_pipelined(
    MyList<gridseg> **src_a, MyList<gridseg> **dst_a,
    MyList<gridseg> **src_b, MyList<gridseg> **dst_b,
    MyList<var> *VarList, int Symmetry
) {
    int myrank, cpusize;
    MPI_Comm_size(MPI_COMM_WORLD, &cpusize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    double **send_d_a = new double *[cpusize];
    double **recv_d_a = new double *[cpusize];
    double **send_h_a = new double *[cpusize];
    double **recv_h_a = new double *[cpusize];
    double **send_d_b = new double *[cpusize];
    double **recv_d_b = new double *[cpusize];
    double **send_h_b = new double *[cpusize];
    double **recv_h_b = new double *[cpusize];
    int *send_len_a = new int[cpusize];
    int *recv_len_a = new int[cpusize];
    int *send_len_b = new int[cpusize];
    int *recv_len_b = new int[cpusize];

    auto size_and_alloc = [&](MyList<gridseg> **src, MyList<gridseg> **dst,
                              double **send_d, double **recv_d,
                              double **send_h, double **recv_h,
                              int *send_len, int *recv_len) {
        for (int node = 0; node < cpusize; node++) {
            send_d[node] = recv_d[node] = nullptr;
            send_h[node] = recv_h[node] = nullptr;
            send_len[node] = recv_len[node] = 0;
            if (node == myrank) {
                recv_len[node] = gpu_data_packer(0, src[myrank], dst[myrank], node, PACK, VarList, VarList, Symmetry);
                if (recv_len[node] > 0)
                    recv_d[node] = GPUManager::getInstance().allocate_device_memory<double>(recv_len[node]);
            } else {
                send_len[node] = gpu_data_packer(0, src[myrank], dst[myrank], node, PACK, VarList, VarList, Symmetry);
                if (send_len[node] > 0) {
                    send_d[node] = GPUManager::getInstance().allocate_device_memory<double>(send_len[node]);
                    send_h[node] = new double[send_len[node]];
                }
                recv_len[node] = gpu_data_packer(0, src[node], dst[node], node, UNPACK, VarList, VarList, Symmetry);
                if (recv_len[node] > 0) {
                    recv_d[node] = GPUManager::getInstance().allocate_device_memory<double>(recv_len[node]);
                    recv_h[node] = new double[recv_len[node]];
                }
            }
        }
    };

    // Phase 1: size calc + alloc for both A and B
    size_and_alloc(src_a, dst_a, send_d_a, recv_d_a, send_h_a, recv_h_a, send_len_a, recv_len_a);
    size_and_alloc(src_b, dst_b, send_d_b, recv_d_b, send_h_b, recv_h_b, send_len_b, recv_len_b);

    // Phase 2: post Irecv A (tag=1) + Irecv B (tag=101) early into host buffers.
    MPI_Request *reqs = new MPI_Request[4 * cpusize];
    MPI_Status *stats = new MPI_Status[4 * cpusize];
    int req_no = 0;
    for (int node = 0; node < cpusize; node++) {
        if (node != myrank && recv_len_a[node] > 0)
            MPI_Irecv((void*)recv_h_a[node], recv_len_a[node], MPI_DOUBLE, node, 1, MPI_COMM_WORLD, reqs + req_no++);
    }
    for (int node = 0; node < cpusize; node++) {
        if (node != myrank && recv_len_b[node] > 0)
            MPI_Irecv((void*)recv_h_b[node], recv_len_b[node], MPI_DOUBLE, node, 101, MPI_COMM_WORLD, reqs + req_no++);
    }

    // Phase 3: ensure async device allocations complete on memory_stream
    GPUManager::getInstance().synchronize_memory();

    // Phase 4: PACK A then PACK B (shared pack_streams set dedups across phases)
    std::vector<cudaStream_t> pack_streams;
    pack_streams.reserve(16);
    for (int node = 0; node < cpusize; node++) {
        if (node == myrank) {
            if (recv_len_a[node] > 0)
                gpu_data_packer(recv_d_a[node], src_a[myrank], dst_a[myrank], node, PACK,
                                VarList, VarList, Symmetry, &pack_streams);
        } else {
            if (send_len_a[node] > 0)
                gpu_data_packer(send_d_a[node], src_a[myrank], dst_a[myrank], node, PACK,
                                VarList, VarList, Symmetry, &pack_streams);
        }
    }
    for (int node = 0; node < cpusize; node++) {
        if (node == myrank) {
            if (recv_len_b[node] > 0)
                gpu_data_packer(recv_d_b[node], src_b[myrank], dst_b[myrank], node, PACK,
                                VarList, VarList, Symmetry, &pack_streams);
        } else {
            if (send_len_b[node] > 0)
                gpu_data_packer(send_d_b[node], src_b[myrank], dst_b[myrank], node, PACK,
                                VarList, VarList, Symmetry, &pack_streams);
        }
    }

    // Phase 5: targeted per-stream sync, then D2H copies for both phases.
    for (auto s : pack_streams) cudaStreamSynchronize(s);
    for (int node = 0; node < cpusize; node++) {
        if (node != myrank && send_len_a[node] > 0)
            GPUManager::getInstance().sync_to_cpu(send_h_a[node], send_d_a[node], send_len_a[node]);
        if (node != myrank && send_len_b[node] > 0)
            GPUManager::getInstance().sync_to_cpu(send_h_b[node], send_d_b[node], send_len_b[node]);
    }
    GPUManager::getInstance().synchronize_memory();

    // Phase 6: post Isend A + Isend B
    for (int node = 0; node < cpusize; node++) {
        if (node != myrank && send_len_a[node] > 0)
            MPI_Isend((void*)send_h_a[node], send_len_a[node], MPI_DOUBLE, node, 1, MPI_COMM_WORLD, reqs + req_no++);
    }
    for (int node = 0; node < cpusize; node++) {
        if (node != myrank && send_len_b[node] > 0)
            MPI_Isend((void*)send_h_b[node], send_len_b[node], MPI_DOUBLE, node, 101, MPI_COMM_WORLD, reqs + req_no++);
    }
    MPI_Waitall(req_no, reqs, stats);

    // Phase 7: H2D copies for both phases' recv data
    for (int node = 0; node < cpusize; node++) {
        if (node != myrank && recv_len_a[node] > 0)
            GPUManager::getInstance().sync_to_gpu(recv_h_a[node], recv_d_a[node], recv_len_a[node]);
        if (node != myrank && recv_len_b[node] > 0)
            GPUManager::getInstance().sync_to_gpu(recv_h_b[node], recv_d_b[node], recv_len_b[node]);
    }
    GPUManager::getInstance().synchronize_memory();

    // Phase 8: UNPACK A then UNPACK B
    std::vector<cudaStream_t> unpack_streams;
    unpack_streams.reserve(16);
    for (int node = 0; node < cpusize; node++) {
        if (node == myrank) {
            if (recv_d_a[node])
                gpu_data_packer(recv_d_a[node], src_a[node], dst_a[node], node, UNPACK,
                                VarList, VarList, Symmetry, &unpack_streams);
        } else {
            if (recv_h_a[node])
                gpu_data_packer(recv_d_a[node], src_a[node], dst_a[node], node, UNPACK,
                                VarList, VarList, Symmetry, &unpack_streams);
        }
    }
    for (int node = 0; node < cpusize; node++) {
        if (node == myrank) {
            if (recv_d_b[node])
                gpu_data_packer(recv_d_b[node], src_b[node], dst_b[node], node, UNPACK,
                                VarList, VarList, Symmetry, &unpack_streams);
        } else {
            if (recv_h_b[node])
                gpu_data_packer(recv_d_b[node], src_b[node], dst_b[node], node, UNPACK,
                                VarList, VarList, Symmetry, &unpack_streams);
        }
    }

    // Phase 9: host_wait merged unpack streams, free all buffers
    for (auto s : unpack_streams) cudaStreamSynchronize(s);

    for (int node = 0; node < cpusize; node++) {
        if (send_d_a[node]) GPUManager::getInstance().free_device_memory(send_d_a[node]);
        if (recv_d_a[node]) GPUManager::getInstance().free_device_memory(recv_d_a[node]);
        if (send_d_b[node]) GPUManager::getInstance().free_device_memory(send_d_b[node]);
        if (recv_d_b[node]) GPUManager::getInstance().free_device_memory(recv_d_b[node]);
        if (send_h_a[node]) delete[] send_h_a[node];
        if (recv_h_a[node]) delete[] recv_h_a[node];
        if (send_h_b[node]) delete[] send_h_b[node];
        if (recv_h_b[node]) delete[] recv_h_b[node];
    }

    delete[] reqs;
    delete[] stats;
    delete[] send_d_a;
    delete[] recv_d_a;
    delete[] send_d_b;
    delete[] recv_d_b;
    delete[] send_h_a;
    delete[] recv_h_a;
    delete[] send_h_b;
    delete[] recv_h_b;
    delete[] send_len_a;
    delete[] recv_len_a;
    delete[] send_len_b;
    delete[] recv_len_b;
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
    int cpusize;
    MPI_Comm_size(MPI_COMM_WORLD, &cpusize);

    // Kill switch: AMSS_SYNC_PIPELINED=0 forces the legacy two-call path.
    // Cached after first check so the env lookup doesn't hit every Sync.
    static int pipelined_enabled = -1;
    if (pipelined_enabled < 0) {
        const char *env = std::getenv("AMSS_SYNC_PIPELINED");
        pipelined_enabled = (env && env[0] == '0') ? 0 : 1;
    }

    // ---- Build Phase A (ghost concat across patches) ----
    MyList<Parallel::gridseg> **all_xfer_src = new MyList<Parallel::gridseg>*[cpusize];
    MyList<Parallel::gridseg> **all_xfer_dst = new MyList<Parallel::gridseg>*[cpusize];
    for (int node = 0; node < cpusize; node++) {
        all_xfer_src[node] = nullptr;
        all_xfer_dst[node] = nullptr;
    }

    MyList<Patch> *Pp = PatL;
    while (Pp) {
        MyList<Parallel::gridseg> *dst_pp = build_ghost_gsl(Pp->data);
        for (int node = 0; node < cpusize; node++) {
            MyList<Parallel::gridseg> *src_pp = build_owned_gsl0(Pp->data, node);
            MyList<Parallel::gridseg> *xs = nullptr, *xd = nullptr;
            build_gstl(src_pp, dst_pp, &xs, &xd);
            if (all_xfer_src[node]) all_xfer_src[node]->catList(xs);
            else                     all_xfer_src[node] = xs;
            if (all_xfer_dst[node]) all_xfer_dst[node]->catList(xd);
            else                     all_xfer_dst[node] = xd;
            if (src_pp) src_pp->destroyList();
        }
        if (dst_pp) dst_pp->destroyList();
        Pp = Pp->next;
    }

    if (!pipelined_enabled) {
        // Legacy fallback: two independent gpu_transfer calls
        gpu_transfer(all_xfer_src, all_xfer_dst, VarList, VarList, Symmetry);

        for (int node = 0; node < cpusize; node++) {
            if (all_xfer_src[node]) all_xfer_src[node]->destroyList();
            if (all_xfer_dst[node]) all_xfer_dst[node]->destroyList();
        }
        delete[] all_xfer_src;
        delete[] all_xfer_dst;

        MyList<Parallel::gridseg> *dst;
        MyList<Parallel::gridseg> **src, **transfer_src, **transfer_dst;
        src = new MyList<Parallel::gridseg> *[cpusize];
        transfer_src = new MyList<Parallel::gridseg> *[cpusize];
        transfer_dst = new MyList<Parallel::gridseg> *[cpusize];

        dst = build_buffer_gsl(PatL);
        for (int node = 0; node < cpusize; node++) {
            src[node] = build_owned_gsl(PatL, node, 5, Symmetry);
            build_gstl(src[node], dst, &transfer_src[node], &transfer_dst[node]);
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
        return;
    }

    // ---- Build Phase B (inter-patch buffer) ----
    MyList<Parallel::gridseg> *dst_buf = build_buffer_gsl(PatL);
    MyList<Parallel::gridseg> **src_buf       = new MyList<Parallel::gridseg> *[cpusize];
    MyList<Parallel::gridseg> **xfer_src_b    = new MyList<Parallel::gridseg> *[cpusize];
    MyList<Parallel::gridseg> **xfer_dst_b    = new MyList<Parallel::gridseg> *[cpusize];
    for (int node = 0; node < cpusize; node++) {
        src_buf[node] = build_owned_gsl(PatL, node, 5, Symmetry);
        build_gstl(src_buf[node], dst_buf, &xfer_src_b[node], &xfer_dst_b[node]);
    }

    // Single pipelined transfer: A's Irecv + B's Irecv both posted before any PACK
    gpu_transfer_pipelined(all_xfer_src, all_xfer_dst,
                           xfer_src_b, xfer_dst_b, VarList, Symmetry);

    // Cleanup phase A
    for (int node = 0; node < cpusize; node++) {
        if (all_xfer_src[node]) all_xfer_src[node]->destroyList();
        if (all_xfer_dst[node]) all_xfer_dst[node]->destroyList();
    }
    delete[] all_xfer_src;
    delete[] all_xfer_dst;

    // Cleanup phase B
    if (dst_buf) dst_buf->destroyList();
    for (int node = 0; node < cpusize; node++) {
        if (src_buf[node]) src_buf[node]->destroyList();
        if (xfer_src_b[node]) xfer_src_b[node]->destroyList();
        if (xfer_dst_b[node]) xfer_dst_b[node]->destroyList();
    }
    delete[] src_buf;
    delete[] xfer_src_b;
    delete[] xfer_dst_b;
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

// List-level OutBdLow2Hi_GPU: fuse N_c × N_f single-Patch fanout into a
// single gpu_transfer call by concatenating the per-(Patc, Patf) transfer plans.
//
// Per-pair build_owned_gsl4(Patc, node, Symmetry) / build_buffer_gsl(Patf) /
// build_gstl produce the same transfer pairs as the original double loop;
// only the MPI exchange is collapsed. Different fine Patches' buffer regions
// don't overlap, so this is geometrically equivalent to the original.
void Parallel::OutBdLow2Hi_GPU(
    MyList<Patch> *PatcL, MyList<Patch> *PatfL,
    MyList<var> *VarList1 /* source */, MyList<var> *VarList2 /* target */,
    int Symmetry
) {
    if (!PatcL || !PatfL) return;
    if (PatcL->data->lev >= PatfL->data->lev) {
        cout << "Parallel::OutBdLow2Hi: meet requst of Prolong from lev#" << PatfL->data->lev
             << " to lev#" << PatcL->data->lev << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int cpusize;
    MPI_Comm_size(MPI_COMM_WORLD, &cpusize);

    MyList<Parallel::gridseg> **all_xfer_src = new MyList<Parallel::gridseg>*[cpusize];
    MyList<Parallel::gridseg> **all_xfer_dst = new MyList<Parallel::gridseg>*[cpusize];
    for (int node = 0; node < cpusize; node++) {
        all_xfer_src[node] = nullptr;
        all_xfer_dst[node] = nullptr;
    }

    MyList<Patch> *Ppc = PatcL;
    while (Ppc) {
        MyList<Patch> *Pp = PatfL;
        while (Pp) {
            MyList<Parallel::gridseg> *dst_pp = build_buffer_gsl(Pp->data); // buffer region of fine Patch
            for (int node = 0; node < cpusize; node++) {
                MyList<Parallel::gridseg> *src_pp = build_owned_gsl4(Ppc->data, node, Symmetry);
                MyList<Parallel::gridseg> *xs = nullptr, *xd = nullptr;
                build_gstl(src_pp, dst_pp, &xs, &xd);
                if (all_xfer_src[node]) all_xfer_src[node]->catList(xs);
                else                     all_xfer_src[node] = xs;
                if (all_xfer_dst[node]) all_xfer_dst[node]->catList(xd);
                else                     all_xfer_dst[node] = xd;
                if (src_pp) src_pp->destroyList();
            }
            if (dst_pp) dst_pp->destroyList();
            Pp = Pp->next;
        }
        Ppc = Ppc->next;
    }

    gpu_transfer(all_xfer_src, all_xfer_dst, VarList1, VarList2, Symmetry);

    for (int node = 0; node < cpusize; node++) {
        if (all_xfer_src[node]) all_xfer_src[node]->destroyList();
        if (all_xfer_dst[node]) all_xfer_dst[node]->destroyList();
    }
    delete[] all_xfer_src;
    delete[] all_xfer_dst;
}

void Parallel::gpu_prepare_inter_time_level(
    Patch *Pat,
    MyList<var> *VarList1 /* source (t+dt) */, MyList<var> *VarList2 /* source (t) */,
    MyList<var> *VarList3 /* target (t+a*dt) */, int tindex
) { // Fused version
    int myrank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (tindex != 0 && tindex != 1 && tindex != -1) {
        cout << "error tindex in Parallel::prepare_inter_time_level" << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MyList<Block> *BP = Pat->blb;
    while (BP) {
        Block *cg = BP->data;
        if (myrank == cg->rank) {
            MyList<var> *varl1 = VarList1;
            MyList<var> *varl2 = VarList2;
            MyList<var> *varl3 = VarList3;

            // 分批打包发送给 GPU
            while (varl1 && varl2 && varl3) {
                FusedTimeVars3 fused_vars;
                fused_vars.num_vars = 0;
                
                while (varl1 && varl2 && varl3 && fused_vars.num_vars < MAX_FUSED_VARS) {
                    fused_vars.ptrs1[fused_vars.num_vars] = cg->d_fgfs[varl1->data->sgfn];
                    fused_vars.ptrs2[fused_vars.num_vars] = cg->d_fgfs[varl2->data->sgfn];
                    fused_vars.ptrs3[fused_vars.num_vars] = cg->d_fgfs[varl3->data->sgfn];
                    fused_vars.num_vars++;
                    
                    varl1 = varl1->next;
                    varl2 = varl2->next;
                    varl3 = varl3->next;
                }
                
                gpu_fused_time_interp_3_launch(cg->stream, cg->shape, cg->shape_pad, fused_vars, tindex);
            }
        }
        if (BP == Pat->ble)
            break;
        BP = BP->next;
    }
}

void Parallel::gpu_prepare_inter_time_level(
    Patch *Pat,
    MyList<var> *VarList1 /* source (t+dt) */, MyList<var> *VarList2 /* source (t) */,
    MyList<var> *VarList3 /* source (t-dt) */, MyList<var> *VarList4 /* target (t+a*dt) */, int tindex
) { // Fused version
    int myrank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (tindex != 0 && tindex != 1 && tindex != -1) {
        cout << "error tindex in long cgh::prepare_inter_time_level" << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MyList<Block> *BP = Pat->blb;
    while (BP) {
        Block *cg = BP->data;
        if (myrank == cg->rank) {
            MyList<var> *varl1 = VarList1;
            MyList<var> *varl2 = VarList2;
            MyList<var> *varl3 = VarList3;
            MyList<var> *varl4 = VarList4;

            while (varl1 && varl2 && varl3 && varl4) {
                FusedTimeVars4 fused_vars;
                fused_vars.num_vars = 0;

                while (varl1 && varl2 && varl3 && varl4 && fused_vars.num_vars < MAX_FUSED_VARS) {
                    fused_vars.ptrs1[fused_vars.num_vars] = cg->d_fgfs[varl1->data->sgfn];
                    fused_vars.ptrs2[fused_vars.num_vars] = cg->d_fgfs[varl2->data->sgfn];
                    fused_vars.ptrs3[fused_vars.num_vars] = cg->d_fgfs[varl3->data->sgfn];
                    fused_vars.ptrs4[fused_vars.num_vars] = cg->d_fgfs[varl4->data->sgfn];
                    fused_vars.num_vars++;

                    varl1 = varl1->next;
                    varl2 = varl2->next;
                    varl3 = varl3->next;
                    varl4 = varl4->next;
                }

                gpu_fused_time_interp_4_launch(cg->stream, cg->shape, cg->shape_pad, fused_vars, tindex);
            }
        }
        if (BP == Pat->ble)
            break;
        BP = BP->next;
    }
}

bool Parallel::PatList_Interp_Points_GPU(
    cudaStream_t stream,
    MyList<Patch> *PatL, MyList<var> *VarList,
    int NN, double *d_XX[3], 
    double *d_Shellf,
    int Symmetry
) {
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    MyList<var> *varl = VarList;
    int num_var = 0;
    while (varl) {
        num_var++;
        varl = varl->next;
    }

    if (!PatL || !PatL->data) {
        return false;
    }

    int ordn = 2 * ghost_width;

    Patch *selected_patch = nullptr;
    Block *selected_block = nullptr;
    if (NN == 1) {
        double h_XX[3] = {0.0, 0.0, 0.0};
        for (int i = 0; i < dim; i++) {
            GPUManager::getInstance().sync_to_cpu(&h_XX[i], d_XX[i], 1, stream);
        }
        GPUManager::synchronize_stream(stream);

        MyList<Patch> *PL_select = PatL;
        while (PL_select && !selected_block) {
            Patch *patch = PL_select->data;
            double DH[3] = {0.0, 0.0, 0.0};
            bool patch_hit = true;
            for (int i = 0; i < dim; i++) {
                DH[i] = patch->getdX(i);
                double lld = patch->lli[i] * DH[i];
                double uud = patch->uui[i] * DH[i];
                if (h_XX[i] < patch->bbox[i] + lld || h_XX[i] > patch->bbox[dim + i] - uud) {
                    patch_hit = false;
                    break;
                }
            }
            if (patch_hit) {
                MyList<Block> *Bp_select = patch->blb;
                while (Bp_select) {
                    Block *BP = Bp_select->data;
                    bool block_hit = true;
                    for (int i = 0; i < dim; i++) {
                        double llb = (feq(BP->bbox[i], patch->bbox[i], DH[i] / 2)) ? BP->bbox[i] + patch->lli[i] * DH[i] : BP->bbox[i] + ghost_width * DH[i];
                        double uub = (feq(BP->bbox[dim + i], patch->bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - patch->uui[i] * DH[i] : BP->bbox[dim + i] - ghost_width * DH[i];
                        if (h_XX[i] - llb < -DH[i] / 2 || h_XX[i] - uub > DH[i] / 2) {
                            block_hit = false;
                            break;
                        }
                    }
                    if (block_hit) {
                        selected_patch = patch;
                        selected_block = BP;
                        break;
                    }
                    if (Bp_select == patch->ble) break;
                    Bp_select = Bp_select->next;
                }
            }
            PL_select = PL_select->next;
        }

        // NN==1 hot path: redirect *all* subsequent GPU work for this call
        // onto the selected block's own stream. The PACK kernel below reads
        // selected_block->d_fgfs[...] which was just written by the RHS/RK4
        // chain on that same stream — by launching PACK there too we get
        // automatic same-stream serialization (no event-wait, no host stall).
        // The comm-pool buffers we allocate just below also follow this
        // stream, so PACK reads the freshly memset'd zeros in order.
        // Previously this function ran PACK on the caller's pool stream,
        // which required either a cudaStreamSynchronize or cudaStreamWaitEvent
        // against selected_block->stream to avoid a stale-data race
        // amplified into ~1e-3 BH drift over ~10 RK4 steps (chaotic BBH).
        if (selected_block && selected_block->rank == myrank) {
            stream = selected_block->stream;
        }
    }

    double *d_local_shellf = GPUManager::getInstance().allocate_comm_buffer<double>(NN * num_var, stream);
    int *d_local_weight = GPUManager::getInstance().allocate_comm_buffer<int>(NN, stream);
    GPUManager::getInstance().memset_device_memory(d_local_shellf, NN * num_var, stream);
    GPUManager::getInstance().memset_device_memory(d_local_weight, NN, stream);

    MyList<Patch> *PL = PatL;
    while (PL) {
        Patch *patch = PL->data;
        if (NN == 1 && selected_patch && patch != selected_patch) {
            PL = PL->next;
            continue;
        }

        double *DH = new double[dim];
        for (int i = 0; i < dim; i++) DH[i] = patch->getdX(i);

        MyList<Block> *Bp = patch->blb;
        while (Bp) {
            Block *BP = Bp->data;
            if (NN == 1 && selected_block && BP != selected_block) {
                if (Bp == patch->ble) break;
                Bp = Bp->next;
                continue;
            }
            if (myrank == BP->rank) {
                double llb[3] = {0}, uub[3] = {0};
                for (int i = 0; i < dim; i++) {
                    llb[i] = (feq(BP->bbox[i], patch->bbox[i], DH[i] / 2)) ? BP->bbox[i] + patch->lli[i] * DH[i] : BP->bbox[i] + ghost_width * DH[i];
                    uub[i] = (feq(BP->bbox[dim + i], patch->bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - patch->uui[i] * DH[i] : BP->bbox[dim + i] - ghost_width * DH[i];
                }

                int shape_0 = BP->shape[0];
                int shape_1 = BP->shape[1];
                int shape_2 = BP->shape[2];

                varl = VarList;
                int k = 0;

                while (varl) {
                    FusedVars fused_vars;
                    fused_vars.num_vars = 0;

                    while (varl && fused_vars.num_vars < MAX_FUSED_VARS) {
                        fused_vars.ptrs[fused_vars.num_vars] = BP->d_fgfs[varl->data->sgfn];
                        fused_vars.SoA0[fused_vars.num_vars] = varl->data->SoA[0];
                        fused_vars.SoA1[fused_vars.num_vars] = varl->data->SoA[1];
                        fused_vars.SoA2[fused_vars.num_vars] = varl->data->SoA[2];

                        fused_vars.num_vars++;
                        varl = varl->next;
                    }

                    gpu_fused_global_interp_launch(
                        stream,
                        NN, dim,
                        d_XX[0], d_XX[1], d_XX[2],
                        shape_0, shape_1, shape_2,
                        BP->shape_pad[0], BP->shape_pad[1], BP->shape_pad[2],
                        BP->d_X[0], BP->d_X[1], BP->d_X[2],
                        fused_vars,
                        llb[0], llb[1], llb[2],
                        uub[0], uub[1], uub[2],
                        DH[0], DH[1], DH[2],
                        ordn, Symmetry,
                        k, num_var, d_local_shellf, d_local_weight
                    );

                    k += fused_vars.num_vars;
                }
            }
            if (Bp == patch->ble) break;
            Bp = Bp->next;
        }
        delete[] DH;
        PL = PL->next;
    }

    int *d_global_weight = GPUManager::getInstance().allocate_comm_buffer<int>(NN, stream);
    GPUManager::getInstance().memset_device_memory(d_global_weight, NN, stream);

    // BH-tracking hot path: this function fires per RK4 substage × BH × level
    // (compute_Porg_rhs descent), and the message is tiny — NN=1, num_var=3 →
    // 3 doubles + 1 int per Allreduce, sub-µs of wire time vs ~100 µs of MPI
    // protocol overhead. Host-stage + manual gather-then-bcast wins regardless
    // of MPI_CUDA_AWARE: cuda-aware MPI adds ~30 µs of comm-buffer alloc and
    // device rendezvous per call on top of the protocol cost.
    //
    // Pack shellf + weight into one double buffer for a single gather+bcast:
    // doubles the wire payload (weight goes from int to double) but halves the
    // protocol-overhead count, which dominates at this message size. Layout:
    //   [0 .. NN*num_var)              : shellf
    //   [NN*num_var .. NN*num_var + NN) : weight (cast back to int on unpack)
    int packed_n = NN * num_var + NN;
    std::vector<double> h_packed(packed_n);
    std::vector<int>    h_weight(NN);
    GPUManager::getInstance().sync_to_cpu(h_packed.data(), d_local_shellf, NN * num_var, stream);
    GPUManager::getInstance().sync_to_cpu(h_weight.data(), d_local_weight, NN, stream);
    GPUManager::synchronize_stream(stream);
    for (int i = 0; i < NN; i++) h_packed[NN * num_var + i] = static_cast<double>(h_weight[i]);

    amss_reduce::allreduce_sum_inplace(h_packed.data(), packed_n, MPI_COMM_WORLD);

    // Host-side normalize + found-check: replaces the gpu_normalize_shellf_launch
    // + d_err_flag D2H + synchronize_stream that previously cost one host stall
    // per BH-tracking call. Data is already on host after the allreduce, NN is 1
    // for the BH-tracking caller, so dividing on the host is essentially free.
    bool success = true;
    for (int i = 0; i < NN; i++) {
        int w = static_cast<int>(h_packed[NN * num_var + i]);
        h_weight[i] = w;
        if (w <= 0) { success = false; continue; }
        for (int k = 0; k < num_var; k++)
            h_packed[k + i * num_var] /= static_cast<double>(w);
    }
    if (!success) checkpatchlist(PatL, false);

    // H2D back to device on the same pool stream the caller passed in. Caller's
    // post-call read of d_Shellf is on memory_stream (sync_to_cpu default) and
    // is preceded by its own synchronize_memory(), so we record an event on
    // stream and have memory_stream wait — no host stall, but the cross-stream
    // dependency is honored.
    GPUManager::getInstance().sync_to_gpu(h_packed.data(), d_Shellf, NN * num_var, stream);
    GPUManager::getInstance().sync_to_gpu(h_weight.data(), d_global_weight, NN, stream);
    {
        cudaEvent_t e = GPUManager::getInstance().event_acquire();
        GPUManager::getInstance().event_record(e, stream);
        cudaStreamWaitEvent(GPUManager::getInstance().memory_stream(), e, 0);
        GPUManager::getInstance().event_release(e);
    }

    GPUManager::getInstance().free_comm_buffer(d_local_shellf, NN * num_var, stream);
    GPUManager::getInstance().free_comm_buffer(d_local_weight, NN, stream);
    GPUManager::getInstance().free_comm_buffer(d_global_weight, NN, stream);

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
        prepare_inter_time_level(PatL->data, VarList1, VarList2, VarList3, tindex);
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

    double tvf = 0.0;
    int BDW = ghost_width;

    // 1. 在显存上分配当前进程的本地求和标量，并清零
    double* d_local_sum = GPUManager::getInstance().allocate_comm_buffer<double>(1);
    GPUManager::getInstance().memset_device_memory(d_local_sum, 1);
    GPUManager::getInstance().synchronize_memory();

    // 2. 异步发射所有 Block 的 Kernel，它们会自动并发地归约到 d_local_sum

    // 2. 异步发射所有 Block 的 Kernel，它们会自动并发地归约到 d_local_sum
    MyList<Block> *BP = Pat->blb;
    while (BP) {
        Block *cg = BP->data;
        if (myrank == cg->rank) {
            gpu_l2normhelper_launch(
                cg->stream,
                cg->shape, cg->shape_pad,
                cg->X[0], cg->X[1], cg->X[2],
                Pat->bbox[0], Pat->bbox[1], Pat->bbox[2],
                Pat->bbox[3], Pat->bbox[4], Pat->bbox[5],
                cg->d_fgfs[vf->sgfn], d_local_sum, BDW
            );
        }
        if (BP == Pat->ble)
            break;
        BP = BP->next;
    }

    GPUManager::getInstance().synchronize_all();

    // 1-double Allreduce: cuda-aware path adds ~30 µs of comm-buffer alloc +
    // device-side rendezvous setup for a single double — strictly slower than
    // staging through host. Always host-stage tiny reductions.
    //
    // Only rank 0 ever reads the returned ConV[] (writefile gates on I_Print),
    // so use sum-to-root: non-root ranks return their local partial — value is
    // junk on them, but no caller looks. Saves the Allreduce all-ranks-result
    // overhead at small N.
    double h_local_sum = 0.0;
    GPUManager::getInstance().sync_to_cpu(&h_local_sum, d_local_sum, 1);
    GPUManager::getInstance().synchronize_memory();

    tvf = h_local_sum;
    amss_reduce::reduce_sum_to_root_inplace(&tvf, 1, 0, MPI_COMM_WORLD);

    GPUManager::getInstance().free_comm_buffer(d_local_sum, 1);
    return sqrt(tvf);
}

void Parallel::Dump_Data_GPU(MyList<Patch> *PL, MyList<var> *DumpList, char *tag, double time, double dT) {
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    Helper::move_to_cpu_whole(PL, myrank, DumpList);
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
    Helper::move_to_cpu_whole(PL, myrank, DumpList);
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