// ShellPatch_GPU.cpp
// GPU-native Shell patch synchronization.
// Mirrors the CPU logic in ShellPatch.C (Sync, intertransfer, Synch, CS_Inter)
// but operates on d_fgfs arrays using gpu_transfer and GPU interpolation kernels.

#include "ShellPatch.h"
#include "Parallel.h"
#include "helper.h"
#include "gpu_manager.h"
#include "fmisc.h"
#include <cstdint>
#include <algorithm>
#include <set>
#include <numeric>

#ifdef USE_GPU

// GPU kernel declarations (implemented in ShellPatch_gpu.cu)
#include "ShellPatch_gpu.cuh"
#include "macrodef.h"
#include <numeric>

#define PI M_PI

// ---------------------------------------------------------------------------
// Helper: ensure coef/sind are computed for a pointstru (lazy init)
// This mirrors the lazy-computation block in ShellPatch::interdata_packer.
// Must be called before using coef/sind in GPU kernels.
// ---------------------------------------------------------------------------
static void ensure_interp_coeffs(ShellPatch::pointstru *ps) {
    if (ps->coef != 0) return; // already computed

    int ordn = 2 * ghost_width;
    int DIMh = (ps->dumyd == -1) ? dim : 1;

    if (DIMh == 3) {
        ps->coef = new double[ordn * DIMh];
        ps->sind = new int[dim];
        for (int i = 0; i < DIMh; i++) {
            double dd = ps->Bg->getdX(i);
            ps->sind[i] = int((ps->lpox[i] - ps->Bg->X[i][0]) / dd) - ordn / 2 + 1;
            double h1;
            for (int j = 0; j < ordn; j++) {
                h1 = ps->Bg->X[i][0] + (ps->sind[i] + j) * dd;
                ps->coef[i * ordn + j] = 1;
                for (int k = 0; k < j; k++) {
                    double h2 = ps->Bg->X[i][0] + (ps->sind[i] + k) * dd;
                    ps->coef[i * ordn + j] *= (ps->lpox[i] - h2) / (h1 - h2);
                }
                for (int k = j + 1; k < ordn; k++) {
                    double h2 = ps->Bg->X[i][0] + (ps->sind[i] + k) * dd;
                    ps->coef[i * ordn + j] *= (ps->lpox[i] - h2) / (h1 - h2);
                }
            }
        }
    } else {
        // DIMh == 1
        int actd = 1 - ps->dumyd;
        ps->sind = new int[dim];
        ps->coef = new double[ordn];
        double dd = ps->Bg->getdX(actd);
        ps->sind[0] = int((ps->lpox[actd] - ps->Bg->X[actd][0]) / dd) - ordn / 2 + 1;
        double h1;
        for (int j = 0; j < ordn; j++) {
            h1 = ps->Bg->X[actd][0] + (ps->sind[0] + j) * dd;
            ps->coef[j] = 1;
            for (int k = 0; k < j; k++) {
                double h2 = ps->Bg->X[actd][0] + (ps->sind[0] + k) * dd;
                ps->coef[j] *= (ps->lpox[actd] - h2) / (h1 - h2);
            }
            for (int k = j + 1; k < ordn; k++) {
                double h2 = ps->Bg->X[actd][0] + (ps->sind[0] + k) * dd;
                ps->coef[j] *= (ps->lpox[actd] - h2) / (h1 - h2);
            }
        }
        ps->sind[2] = int((ps->lpox[2] - ps->Bg->X[2][0]) / ps->Bg->getdX(2) + 0.001);
        ps->sind[1] = int((ps->lpox[ps->dumyd] - ps->Bg->X[ps->dumyd][0]) / ps->Bg->getdX(ps->dumyd) + 0.001);
    }
}

// ---------------------------------------------------------------------------
// ss_patch::Sync_GPU
// Ghost zone exchange within a single ss_patch, using gpu_transfer.
// Directly mirrors ss_patch::Sync (ShellPatch.C line 177).
// ---------------------------------------------------------------------------
void ss_patch::Sync_GPU(MyList<var> *VarList, int Symmetry) {
    int cpusize;
    MPI_Comm_size(MPI_COMM_WORLD, &cpusize);

    MyList<Parallel::gridseg> *dst;
    MyList<Parallel::gridseg> **src, **transfer_src, **transfer_dst;
    src = new MyList<Parallel::gridseg> *[cpusize];
    transfer_src = new MyList<Parallel::gridseg> *[cpusize];
    transfer_dst = new MyList<Parallel::gridseg> *[cpusize];

    dst = build_ghost_gsl();
    for (int node = 0; node < cpusize; node++) {
        src[node] = build_owned_gsl0(node);
        Parallel::build_gstl(src[node], dst, &transfer_src[node], &transfer_dst[node]);
    }

    Parallel::gpu_transfer(transfer_src, transfer_dst, VarList, VarList, Symmetry);

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

// ---------------------------------------------------------------------------
// ShellPatch::interdata_packer_GPU
// GPU version of ShellPatch::interdata_packer.
//
// PACK: reads from source d_fgfs using GPU interpolation kernels, writes
//       interpolated values into a device buffer.
//       Optimized: allocates sind/coef device buffers ONCE, uploads per
//       pointstru, launches ALL var kernels without intermediate sync.
//
// UNPACK: writes from device buffer to destination d_fgfs.
//         Optimized: bulk D2H of entire buffer at start, then uses host
//         values with gpu_pointcopy_cell_launch (avoids per-point D2H).
//
// size calc pass (d_data == nullptr): just counts points, no GPU ops.
//
// Parameters:
//   d_data   - device pointer for PACK/UNPACK, nullptr for size calc
//   buf_size - total number of doubles in d_data (0 for size calc).
//              Used for bulk D2H in UNPACK path.
// ---------------------------------------------------------------------------
int ShellPatch::interdata_packer_GPU(
    gpu_ps_cache &cache,
    double *d_data,
    int buf_size,
    MyList<pointstru> *src, MyList<pointstru> *dst,
    int list_rank, int rank_in, int dir,
    MyList<var> *VarLists, MyList<var> *VarListd,
    int Symmetry
) {
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int ordn = 2 * ghost_width;

    if (dir != PACK && dir != UNPACK) {
        cout << "error dir " << dir << " for interdata_packer_GPU " << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int size_out = 0;
    if (!src || !dst || cache.num_points == 0) return size_out;

    // Verify VarList lengths match
    MyList<var> *varls = VarLists, *varld = VarListd;
    int num_vars = 0;
    while (varls && varld) { num_vars++; varls = varls->next; varld = varld->next; }
    if (varls || varld) {
        cout << "error in interdata_packer_GPU, var lists does not match." << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (num_vars == 0) return 0;

    // Rank base offset from cache
    int rank_base = cache.h_rank_offsets ? cache.h_rank_offsets[list_rank] : 0;

    // -----------------------------------------------------------------------
    // Size-calc pass (d_data == nullptr): just count points, no GPU work.
    // -----------------------------------------------------------------------
    if (!d_data) {
        while (src && dst) {
            if ((dir == PACK && dst->data->Bg->rank == rank_in && src->data->Bg->rank == myrank) ||
                (dir == UNPACK && src->data->Bg->rank == rank_in && dst->data->Bg->rank == myrank))
            {
                size_out += num_vars;
            }
            dst = dst->next;
            src = src->next;
        }
        return size_out;
    }

    // -----------------------------------------------------------------------
    // Phase A: Build active_indices list and variable metadata on host
    // -----------------------------------------------------------------------

    // Count matching points first
    int num_active = 0;
    {
        MyList<pointstru> *s = src, *d = dst;
        while (s && d) {
            if ((dir == PACK && d->data->Bg->rank == rank_in && s->data->Bg->rank == myrank) ||
                (dir == UNPACK && s->data->Bg->rank == rank_in && d->data->Bg->rank == myrank))
                num_active++;
            s = s->next; d = d->next;
        }
    }
    if (num_active == 0) return 0;

    // Build active_indices: cache_idx = rank_base + position_in_sublist
    int *h_active = new int[num_active];
    {
        int idx = 0;
        int pos = 0;
        MyList<pointstru> *s = src, *d = dst;
        while (s && d) {
            if ((dir == PACK && d->data->Bg->rank == rank_in && s->data->Bg->rank == myrank) ||
                (dir == UNPACK && s->data->Bg->rank == rank_in && d->data->Bg->rank == myrank))
            {
                h_active[idx++] = rank_base + pos;
            }
            pos++;
            s = s->next; d = d->next;
        }
    }

    // Create explicit stream for consistent H2D + kernel launch ordering
    cudaStream_t local_stream = 0;
    CUDA_CHECK(cudaStreamCreate(&local_stream));

    // Upload active indices
    int *d_active = GPUManager::getInstance().allocate_device_memory<int>(num_active, local_stream);
    GPUManager::getInstance().sync_to_gpu(h_active, d_active, num_active, local_stream);

    // Build variable metadata arrays
    double *h_SoA = new double[num_vars * 3];
    int    *h_var_fgfs = new int[num_vars];
    {
        varls = VarLists; varld = VarListd;
        for (int v = 0; v < num_vars; v++) {
            h_SoA[v * 3 + 0] = varls->data->SoA[0];
            h_SoA[v * 3 + 1] = varls->data->SoA[1];
            h_SoA[v * 3 + 2] = varls->data->SoA[2];
            h_var_fgfs[v] = varls->data->sgfn;
            varls = varls->next; varld = varld->next;
        }
    }

    double *d_SoA = GPUManager::getInstance().allocate_device_memory<double>(num_vars * 3, local_stream);
    int    *d_var_fgfs = GPUManager::getInstance().allocate_device_memory<int>(num_vars, local_stream);
    GPUManager::getInstance().sync_to_gpu(h_SoA, d_SoA, num_vars * 3, local_stream);
    GPUManager::getInstance().sync_to_gpu(h_var_fgfs, d_var_fgfs, num_vars, local_stream);

    // -----------------------------------------------------------------------
    // Phase B: Launch batched kernel on same stream as H2D
    // -----------------------------------------------------------------------
    if (dir == PACK) {
        // PACK must wait for the selected source writers.
        std::set<Block *> waited;
        auto *s = src;
        auto *d = dst;
        for (; s && d; s = s->next, d = d->next) {
            Block *source = s->data->Bg;
            if (d->data->Bg->rank != rank_in || source->rank != myrank ||
                !waited.insert(source).second) continue;
            cudaEvent_t ready;
            CUDA_CHECK(cudaEventCreateWithFlags(&ready, cudaEventDisableTiming));
            CUDA_CHECK(cudaEventRecord(ready, source->stream));
            CUDA_CHECK(cudaStreamWaitEvent(local_stream, ready, 0));
            CUDA_CHECK(cudaEventDestroy(ready));
        }
        gpu_batched_interpind_pack_launch(
            local_stream,
            num_active, num_vars,
            d_active,
            cache.d_coef, cache.d_sind,
            cache.d_dumyd, cache.d_ssst,
            cache.d_src_block_idx,
            cache.max_coef,
            d_gpu_block_fgfs,
            d_gpu_block_shape, d_gpu_block_shape_pad,
            ordn, fngfs,
            d_SoA, d_var_fgfs,
            d_data + size_out
        );
    } else {
        // UNPACK: write directly from GPU buffer to d_fgfs — no D2H!
        gpu_batched_pointcopy_unpack_launch(
            local_stream,
            num_active, num_vars,
            d_active,
            cache.d_dst_lpox,
            cache.d_dst_block_idx,
            d_gpu_block_fgfs,
            d_gpu_block_shape, d_gpu_block_shape_pad,
            d_gpu_block_bbox,
            fngfs,
            d_var_fgfs,
            d_data + size_out
        );
    }

    size_out = num_active * num_vars;

    // Keep allocations, uploads, kernel, and frees ordered on one stream.
    GPUManager::getInstance().free_device_memory(d_active, 0, local_stream);
    GPUManager::getInstance().free_device_memory(d_SoA, 0, local_stream);
    GPUManager::getInstance().free_device_memory(d_var_fgfs, 0, local_stream);
    // H2D copies may still read the host staging arrays until this completes.
    CUDA_CHECK(cudaStreamSynchronize(local_stream));
    delete[] h_active;
    delete[] h_SoA;
    delete[] h_var_fgfs;
    CUDA_CHECK(cudaStreamDestroy(local_stream));

    return size_out;
}

// ---------------------------------------------------------------------------
// ShellPatch::intertransfer_GPU
// GPU version of ShellPatch::intertransfer.
// Uses interdata_packer_GPU for PACK/UNPACK and device buffers for MPI.
// Directly mirrors the structure of Parallel::gpu_transfer.
// ---------------------------------------------------------------------------
void ShellPatch::intertransfer_GPU(
    MyList<pointstru> **src, MyList<pointstru> **dst,
    MyList<var> *VarList1, MyList<var> *VarList2,
    int Symmetry
) {
    int myrank, cpusize;
    MPI_Comm_size(MPI_COMM_WORLD, &cpusize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // Determine which GPU cache to use based on src array
    gpu_ps_cache &cache = (src == ss_src) ? gpu_ss_cache :
                          (src == csatc_src) ? gpu_csatc_cache : gpu_csats_cache;

    MPI_Request *reqs = new MPI_Request[2 * cpusize];
    MPI_Status *stats = new MPI_Status[2 * cpusize];
    int req_no = 0;

#if MPI_CUDA_AWARE
    double **send_data_d = new double *[cpusize];
    double **rec_data_d  = new double *[cpusize];
#else
    double **send_data_d = new double *[cpusize];
    double **rec_data_d  = new double *[cpusize];
    double **send_data_h = new double *[cpusize];
    double **rec_data_h  = new double *[cpusize];
#endif

    int *send_len = new int[cpusize];
    int *recv_len = new int[cpusize];

    // Phase 1: Size calc + allocate
    for (int node = 0; node < cpusize; node++) {
        send_data_d[node] = rec_data_d[node] = nullptr;
        send_len[node] = recv_len[node] = 0;
#if !MPI_CUDA_AWARE
        send_data_h[node] = rec_data_h[node] = nullptr;
#endif

        if (node == myrank) {
            // Self-copy: PACK reads from src[myrank], UNPACK writes to dst[myrank]
            recv_len[node] = interdata_packer_GPU(cache, nullptr, 0, src[myrank], dst[myrank], myrank, node, PACK, VarList1, VarList2, Symmetry);
            if (recv_len[node] > 0)
                rec_data_d[node] = GPUManager::getInstance().allocate_comm_buffer<double>(recv_len[node]);
        } else {
            // Send: PACK from our src to send buffer
            send_len[node] = interdata_packer_GPU(cache, nullptr, 0, src[myrank], dst[myrank], myrank, node, PACK, VarList1, VarList2, Symmetry);
            if (send_len[node] > 0) {
                send_data_d[node] = GPUManager::getInstance().allocate_comm_buffer<double>(send_len[node]);
#if !MPI_CUDA_AWARE
                send_data_h[node] = new double[send_len[node]];
#endif
            }

            // Receive: UNPACK from rec buffer to our dst
            recv_len[node] = interdata_packer_GPU(cache, nullptr, 0, src[node], dst[node], node, node, UNPACK, VarList1, VarList2, Symmetry);
            if (recv_len[node] > 0) {
                rec_data_d[node] = GPUManager::getInstance().allocate_comm_buffer<double>(recv_len[node]);
#if !MPI_CUDA_AWARE
                rec_data_h[node] = new double[recv_len[node]];
#endif
            }
        }
    }

    // Sync recv counts via MPI_Alltoall to avoid MPI_ERR_TRUNCATE
    {
        int *sync_recv = new int[cpusize];
        MPI_Alltoall(send_len, 1, MPI_INT, sync_recv, 1, MPI_INT, MPI_COMM_WORLD);
        for (int node = 0; node < cpusize; node++) {
            if (node != myrank && sync_recv[node] != recv_len[node]) {
                if (rec_data_d[node])
                    GPUManager::getInstance().free_device_memory(rec_data_d[node]);
#if !MPI_CUDA_AWARE
                delete[] rec_data_h[node];
                rec_data_h[node] = nullptr;
#endif
                recv_len[node] = sync_recv[node];
                rec_data_d[node] = nullptr;
                if (recv_len[node] > 0) {
                    rec_data_d[node] = GPUManager::getInstance().allocate_device_memory<double>(recv_len[node]);
#if !MPI_CUDA_AWARE
                    rec_data_h[node] = new double[recv_len[node]];
#endif
                }
            }
        }
        delete[] sync_recv;
    }

    // Phase 1.5: Ensure all async allocations complete
    GPUManager::getInstance().synchronize_memory();

    // Field slots rotate during RK4.
    for (int b = 0; b < gpu_num_blocks; b++)
        for (int v = 0; v < fngfs; v++)
            h_gpu_block_fgfs[b * fngfs + v] = reinterpret_cast<unsigned long long>(gpu_blocks[b]->d_fgfs[v]);
    GPUManager::getInstance().sync_to_gpu(h_gpu_block_fgfs.data(), d_gpu_block_fgfs, gpu_num_blocks * fngfs);
    GPUManager::getInstance().synchronize_memory();

    // Phase 2: All PACK
    for (int node = 0; node < cpusize; node++) {
        if (node == myrank) {
            if (recv_len[node] > 0)
                interdata_packer_GPU(cache, rec_data_d[node], recv_len[node], src[myrank], dst[myrank], myrank, node, PACK, VarList1, VarList2, Symmetry);
        } else {
            if (send_len[node] > 0)
                interdata_packer_GPU(cache, send_data_d[node], send_len[node], src[myrank], dst[myrank], myrank, node, PACK, VarList1, VarList2, Symmetry);
        }
    }

    // Phase 3: Sync all GPU work before MPI
    GPUManager::getInstance().synchronize_all();

#if MPI_CUDA_AWARE
    // Phase 4: MPI communication with device buffers
    for (int node = 0; node < cpusize; node++) {
        if (node != myrank) {
            if (send_len[node] > 0)
                MPI_Isend((void *)send_data_d[node], send_len[node], MPI_DOUBLE, node, 1, MPI_COMM_WORLD, reqs + req_no++);
            if (recv_len[node] > 0)
                MPI_Irecv((void *)rec_data_d[node], recv_len[node], MPI_DOUBLE, node, 1, MPI_COMM_WORLD, reqs + req_no++);
        }
    }
    MPI_Waitall(req_no, reqs, stats);
#else
    // Phase 3.5: D2H copies for send data
    for (int node = 0; node < cpusize; node++) {
        if (node != myrank && send_len[node] > 0)
            GPUManager::getInstance().sync_to_cpu(send_data_h[node], send_data_d[node], send_len[node]);
    }
    GPUManager::getInstance().synchronize_memory();

    // Phase 4: MPI communication with host buffers
    for (int node = 0; node < cpusize; node++) {
        if (node != myrank) {
            if (send_len[node] > 0)
                MPI_Isend((void *)send_data_h[node], send_len[node], MPI_DOUBLE, node, 1, MPI_COMM_WORLD, reqs + req_no++);
            if (recv_len[node] > 0)
                MPI_Irecv((void *)rec_data_h[node], recv_len[node], MPI_DOUBLE, node, 1, MPI_COMM_WORLD, reqs + req_no++);
        }
    }
    MPI_Waitall(req_no, reqs, stats);

    // Phase 5: H2D copies for receive data
    for (int node = 0; node < cpusize; node++) {
        if (node != myrank && recv_len[node] > 0)
            GPUManager::getInstance().sync_to_gpu(rec_data_h[node], rec_data_d[node], recv_len[node]);
    }
    GPUManager::getInstance().synchronize_memory();
#endif

    // Phase 6: All UNPACK
    for (int node = 0; node < cpusize; node++) {
        if (rec_data_d[node])
            interdata_packer_GPU(cache, rec_data_d[node], recv_len[node], src[node], dst[node], node, node, UNPACK, VarList1, VarList2, Symmetry);
    }

    // Phase 7: Sync + cleanup
    GPUManager::getInstance().synchronize_all();

    for (int node = 0; node < cpusize; node++) {
        if (send_data_d[node]) GPUManager::getInstance().free_comm_buffer(send_data_d[node]);
        if (rec_data_d[node])  GPUManager::getInstance().free_comm_buffer(rec_data_d[node]);
#if !MPI_CUDA_AWARE
        if (send_data_h[node]) delete[] send_data_h[node];
        if (rec_data_h[node])  delete[] rec_data_h[node];
#endif
    }

    delete[] reqs;
    delete[] stats;
    delete[] send_data_d;
    delete[] rec_data_d;
    delete[] send_len;
    delete[] recv_len;
#if !MPI_CUDA_AWARE
    delete[] send_data_h;
    delete[] rec_data_h;
#endif
}

// ---------------------------------------------------------------------------
// ShellPatch::Synch_GPU
// GPU version of ShellPatch::Synch.
// Part 1: Ghost zone exchange within each ss_patch using gpu_transfer.
// Part 2: Shell-Shell inter-process interpolation transfer on GPU.
// ---------------------------------------------------------------------------
void ShellPatch::Synch_GPU(MyList<var> *VarList, int Symmetry) {
    // Part 1: Shell patch ghost zone exchange
    MyList<ss_patch> *Pp = PatL;
    while (Pp) {
        Pp->data->Sync_GPU(VarList, Symmetry);
        Pp = Pp->next;
    }

    // Part 2: Shell-Shell inter-process transfer
    intertransfer_GPU(ss_src, ss_dst, VarList, VarList, Symmetry);
}

// ---------------------------------------------------------------------------
// ShellPatch::CS_Inter_GPU
// GPU version of ShellPatch::CS_Inter.
// Cart-Shell inter-process interpolation transfer on GPU.
// ---------------------------------------------------------------------------
void ShellPatch::CS_Inter_GPU(MyList<var> *VarList, int Symmetry) {
    // Fill shell first
    intertransfer_GPU(csats_src, csats_dst, VarList, VarList, Symmetry);
    // Fill box then
    intertransfer_GPU(csatc_src, csatc_dst, VarList, VarList, Symmetry);
}

// ---------------------------------------------------------------------------
// ShellPatch::build_gpu_block_table
// Builds a GPU-resident block metadata table that maps block indices to
// device-side d_fgfs pointers, shape, shape_pad, and bbox.
// Called once during init. Shared by all pointstru caches.
// ---------------------------------------------------------------------------
void ShellPatch::build_gpu_block_table(MyList<Patch> *CPatL) {
    // Count local blocks (those with GPU memory allocated on this rank)
    int max_blocks = 0;
    MyList<Patch> *Pp = CPatL;
    while (Pp) { MyList<Block> *b = Pp->data->blb; while (b) { max_blocks++; b = b->next; } Pp = Pp->next; }
    MyList<ss_patch> *sPp = PatL;
    while (sPp) { MyList<Block> *b = sPp->data->blb; while (b) { max_blocks++; b = b->next; } sPp = sPp->next; }

    Block **h_blocks = new Block *[max_blocks];
    int *h_shape = new int[max_blocks * 3];
    int *h_shape_pad = new int[max_blocks * 3];
    double *h_bbox = new double[max_blocks * 6];
    int num_blocks = 0;

    auto add_block = [&](Block *b) {
        // Only include blocks that belong to this MPI rank (have GPU memory)
        if (b->rank != myrank) return;
        h_blocks[num_blocks] = b;
        for (int d = 0; d < 3; d++) {
            h_shape[num_blocks * 3 + d] = b->shape[d];
            h_shape_pad[num_blocks * 3 + d] = b->shape_pad[d];
        }
        for (int d = 0; d < 6; d++)
            h_bbox[num_blocks * 6 + d] = b->bbox[d];
        num_blocks++;
    };

    Pp = CPatL;
    while (Pp) { MyList<Block> *b = Pp->data->blb; while (b) { add_block(b->data); b = b->next; } Pp = Pp->next; }
    sPp = PatL;
    while (sPp) { MyList<Block> *b = sPp->data->blb; while (b) { add_block(b->data); b = b->next; } sPp = sPp->next; }

    gpu_num_blocks = num_blocks;
    gpu_blocks.assign(h_blocks, h_blocks + num_blocks);
    h_gpu_block_fgfs.resize(num_blocks * fngfs);

    d_gpu_block_shape = GPUManager::getInstance().allocate_device_memory<int>(num_blocks * 3);
    d_gpu_block_shape_pad = GPUManager::getInstance().allocate_device_memory<int>(num_blocks * 3);
    d_gpu_block_bbox = GPUManager::getInstance().allocate_device_memory<double>(num_blocks * 6);
    d_gpu_block_fgfs = GPUManager::getInstance().allocate_device_memory<unsigned long long>(num_blocks * fngfs);

    GPUManager::getInstance().sync_to_gpu(h_shape, d_gpu_block_shape, num_blocks * 3);
    GPUManager::getInstance().sync_to_gpu(h_shape_pad, d_gpu_block_shape_pad, num_blocks * 3);
    GPUManager::getInstance().sync_to_gpu(h_bbox, d_gpu_block_bbox, num_blocks * 6);

    for (int b = 0; b < num_blocks; b++)
        for (int v = 0; v < fngfs; v++)
            h_gpu_block_fgfs[b * fngfs + v] = reinterpret_cast<unsigned long long>(gpu_blocks[b]->d_fgfs[v]);
    GPUManager::getInstance().sync_to_gpu(h_gpu_block_fgfs.data(), d_gpu_block_fgfs, num_blocks * fngfs);
    GPUManager::getInstance().synchronize_memory();

    delete[] h_shape;
    delete[] h_shape_pad;
    delete[] h_bbox;
    delete[] h_blocks;

    if (myrank == 0)
        cout << " ShellPatch:: GPU block table built: " << num_blocks << " blocks" << endl;
}

void ShellPatch::free_gpu_block_table() {
    if (d_gpu_block_shape)     { GPUManager::getInstance().free_device_memory(d_gpu_block_shape); d_gpu_block_shape = nullptr; }
    if (d_gpu_block_shape_pad) { GPUManager::getInstance().free_device_memory(d_gpu_block_shape_pad); d_gpu_block_shape_pad = nullptr; }
    if (d_gpu_block_bbox)      { GPUManager::getInstance().free_device_memory(d_gpu_block_bbox); d_gpu_block_bbox = nullptr; }
    if (d_gpu_block_fgfs)      { GPUManager::getInstance().free_device_memory(d_gpu_block_fgfs); d_gpu_block_fgfs = nullptr; }
    gpu_blocks.clear();
    h_gpu_block_fgfs.clear();
    gpu_num_blocks = 0;
}

// Helper: find block index in the block list (linear search — blocks are few)
static int find_block_index(Block *target, Block **blocks, int num_blocks) {
    for (int i = 0; i < num_blocks; i++)
        if (blocks[i] == target) return i;
    return -1; // should not happen
}

// ---------------------------------------------------------------------------
// ShellPatch::build_gpu_pointstru_cache
// Build persistent GPU cache for a set of pointstru lists (ss, csatc, csats).
// Pre-computes all coef/sind, uploads them once to GPU.
// Called once during init, after setupintintstuff.
// ---------------------------------------------------------------------------
void ShellPatch::build_gpu_pointstru_cache(
    gpu_ps_cache &cache,
    MyList<pointstru> **src_lists,
    MyList<pointstru> **dst_lists,
    int cpusize,
    MyList<Patch> *CPatL)
{
    int ordn = 2 * ghost_width;
    cache.ordn = ordn;
    cache.max_coef = ordn * dim;

    // Build block array for index lookup (same order as build_gpu_block_table)
    // Only include blocks with GPU memory allocated on this rank
    int max_blocks = 0;
    MyList<Patch> *Pp = CPatL;
    while (Pp) { MyList<Block> *b = Pp->data->blb; while (b) { max_blocks++; b = b->next; } Pp = Pp->next; }
    MyList<ss_patch> *sPp = PatL;
    while (sPp) { MyList<Block> *b = sPp->data->blb; while (b) { max_blocks++; b = b->next; } sPp = sPp->next; }

    Block **block_list = new Block *[max_blocks];
    int block_count = 0;
    Pp = CPatL;
    while (Pp) { MyList<Block> *b = Pp->data->blb; while (b) { Block *blk = b->data; if (blk->rank == myrank) { block_list[block_count++] = blk; } b = b->next; } Pp = Pp->next; }
    sPp = PatL;
    while (sPp) { MyList<Block> *b = sPp->data->blb; while (b) { Block *blk = b->data; if (blk->rank == myrank) { block_list[block_count++] = blk; } b = b->next; } sPp = sPp->next; }

    // Phase 1: Count total points per rank (all points, including non-local)
    int *h_rank_counts = new int[cpusize + 1];
    h_rank_counts[0] = 0;
    for (int r = 0; r < cpusize; r++) {
        int cnt = 0;
        MyList<pointstru> *ps = src_lists[r];
        while (ps) { cnt++; ps = ps->next; }
        h_rank_counts[r + 1] = cnt;
    }
    for (int r = 0; r < cpusize; r++)
        h_rank_counts[r + 1] += h_rank_counts[r];
    int total_points = h_rank_counts[cpusize];
    cache.num_points = total_points;

    if (total_points == 0) {
        delete[] h_rank_counts;
        delete[] block_list;
        return;
    }

    // Phase 2: Allocate host arrays
    int coef_total = total_points * cache.max_coef;
    double *h_coef  = new double[coef_total]();
    int    *h_sind  = new int[total_points * 3];
    int    *h_dumyd = new int[total_points];
    int    *h_ssst  = new int[total_points];
    double *h_dst_lpox = new double[total_points * 3];
    int    *h_src_block_idx = new int[total_points];
    int    *h_dst_block_idx = new int[total_points];

    // Phase 3: Fill host arrays (keep original ordering, use sentinels for non-local)
    int global_idx = 0;
    for (int r = 0; r < cpusize; r++) {
        MyList<pointstru> *ps_src = src_lists[r];
        MyList<pointstru> *ps_dst = dst_lists[r];
        while (ps_src && ps_dst) {
            pointstru *src = ps_src->data;
            pointstru *dst = ps_dst->data;

            bool src_local = (src->Bg->rank == myrank);
            bool dst_local = (dst->Bg->rank == myrank);

            // [FIXED] 分离发送端 (Source) 的插值系数计算逻辑
            if (src_local) {
                // Compute coef/sind directly into flat arrays
                int DIMh = (src->dumyd == -1) ? dim : 1;
                double *my_coef = &h_coef[global_idx * cache.max_coef];
                int    *my_sind = &h_sind[global_idx * 3];

                if (DIMh == 3) {
                    for (int d = 0; d < dim; d++) {
                        double dd = src->Bg->getdX(d);
                        my_sind[d] = int((src->lpox[d] - src->Bg->X[d][0]) / dd) - ordn / 2 + 1;
                        for (int j = 0; j < ordn; j++) {
                            double h1 = src->Bg->X[d][0] + (my_sind[d] + j) * dd;
                            my_coef[d * ordn + j] = 1.0;
                            for (int k = 0; k < j; k++) {
                                double h2 = src->Bg->X[d][0] + (my_sind[d] + k) * dd;
                                my_coef[d * ordn + j] *= (src->lpox[d] - h2) / (h1 - h2);
                            }
                            for (int k = j + 1; k < ordn; k++) {
                                double h2 = src->Bg->X[d][0] + (my_sind[d] + k) * dd;
                                my_coef[d * ordn + j] *= (src->lpox[d] - h2) / (h1 - h2);
                            }
                        }
                    }
                } else {
                    // DIMh == 1
                    int actd = 1 - src->dumyd;
                    double dd = src->Bg->getdX(actd);
                    my_sind[0] = int((src->lpox[actd] - src->Bg->X[actd][0]) / dd) - ordn / 2 + 1;
                    for (int j = 0; j < ordn; j++) {
                        double h1 = src->Bg->X[actd][0] + (my_sind[0] + j) * dd;
                        my_coef[j] = 1.0;
                        for (int k = 0; k < j; k++) {
                            double h2 = src->Bg->X[actd][0] + (my_sind[0] + k) * dd;
                            my_coef[j] *= (src->lpox[actd] - h2) / (h1 - h2);
                        }
                        for (int k = j + 1; k < ordn; k++) {
                            double h2 = src->Bg->X[actd][0] + (my_sind[0] + k) * dd;
                            my_coef[j] *= (src->lpox[actd] - h2) / (h1 - h2);
                        }
                    }
                    my_sind[2] = int((src->lpox[2] - src->Bg->X[2][0]) / src->Bg->getdX(2) + 0.001);
                    my_sind[1] = int((src->lpox[src->dumyd] - src->Bg->X[src->dumyd][0]) / src->Bg->getdX(src->dumyd) + 0.001);
                }
                // Resolve Source block index
                h_src_block_idx[global_idx] = find_block_index(src->Bg, block_list, block_count);
            } else {
                h_src_block_idx[global_idx] = -1;
            }

            // [FIXED] 分离接收端 (Destination) 的坐标赋值逻辑
            if (dst_local) {
                h_dst_lpox[global_idx * 3 + 0] = dst->lpox[0];
                h_dst_lpox[global_idx * 3 + 1] = dst->lpox[1];
                h_dst_lpox[global_idx * 3 + 2] = dst->lpox[2];
                // Resolve Destination block index
                h_dst_block_idx[global_idx] = find_block_index(dst->Bg, block_list, block_count);
            } else {
                h_dst_lpox[global_idx * 3 + 0] = 0;
                h_dst_lpox[global_idx * 3 + 1] = 0;
                h_dst_lpox[global_idx * 3 + 2] = 0;
                h_dst_block_idx[global_idx] = -1;
            }

            // [FIXED] 共享属性：发包 (PACK) 和收包 (UNPACK) 的内核都需要使用
            h_dumyd[global_idx] = src->dumyd;
            h_ssst[global_idx]  = src->ssst;

            global_idx++;
            ps_src = ps_src->next;
            ps_dst = ps_dst->next;
        }
    }

    // Phase 4: Upload to GPU
    cache.d_coef = GPUManager::getInstance().allocate_device_memory<double>(coef_total);
    cache.d_sind = GPUManager::getInstance().allocate_device_memory<int>(total_points * 3);
    cache.d_dumyd = GPUManager::getInstance().allocate_device_memory<int>(total_points);
    cache.d_ssst = GPUManager::getInstance().allocate_device_memory<int>(total_points);
    cache.d_dst_lpox = GPUManager::getInstance().allocate_device_memory<double>(total_points * 3);
    cache.d_src_block_idx = GPUManager::getInstance().allocate_device_memory<int>(total_points);
    cache.d_dst_block_idx = GPUManager::getInstance().allocate_device_memory<int>(total_points);
    cache.d_rank_offsets = GPUManager::getInstance().allocate_device_memory<int>(cpusize + 1);
    cache.h_rank_offsets = h_rank_counts; // keep host copy (h_rank_counts already has the data)

    GPUManager::getInstance().sync_to_gpu(h_coef, cache.d_coef, coef_total);
    GPUManager::getInstance().sync_to_gpu(h_sind, cache.d_sind, total_points * 3);
    GPUManager::getInstance().sync_to_gpu(h_dumyd, cache.d_dumyd, total_points);
    GPUManager::getInstance().sync_to_gpu(h_ssst, cache.d_ssst, total_points);
    GPUManager::getInstance().sync_to_gpu(h_dst_lpox, cache.d_dst_lpox, total_points * 3);
    GPUManager::getInstance().sync_to_gpu(h_src_block_idx, cache.d_src_block_idx, total_points);
    GPUManager::getInstance().sync_to_gpu(h_dst_block_idx, cache.d_dst_block_idx, total_points);
    GPUManager::getInstance().sync_to_gpu(h_rank_counts, cache.d_rank_offsets, cpusize + 1);
    GPUManager::getInstance().synchronize_memory();

    delete[] block_list;
    // h_rank_counts is stored as cache.h_rank_offsets — do NOT delete
    delete[] h_coef;
    delete[] h_sind;
    delete[] h_dumyd;
    delete[] h_ssst;
    delete[] h_dst_lpox;
    delete[] h_src_block_idx;
    delete[] h_dst_block_idx;

    if (myrank == 0)
        cout << " ShellPatch:: GPU pointstru cache built: " << total_points
             << " points across " << cpusize << " ranks" << endl;
}

void ShellPatch::free_gpu_pointstru_cache(gpu_ps_cache &cache) {
    if (cache.d_coef)          { GPUManager::getInstance().free_device_memory(cache.d_coef); cache.d_coef = nullptr; }
    if (cache.d_sind)          { GPUManager::getInstance().free_device_memory(cache.d_sind); cache.d_sind = nullptr; }
    if (cache.d_dumyd)         { GPUManager::getInstance().free_device_memory(cache.d_dumyd); cache.d_dumyd = nullptr; }
    if (cache.d_ssst)          { GPUManager::getInstance().free_device_memory(cache.d_ssst); cache.d_ssst = nullptr; }
    if (cache.d_dst_lpox)      { GPUManager::getInstance().free_device_memory(cache.d_dst_lpox); cache.d_dst_lpox = nullptr; }
    if (cache.d_src_block_idx) { GPUManager::getInstance().free_device_memory(cache.d_src_block_idx); cache.d_src_block_idx = nullptr; }
    if (cache.d_dst_block_idx) { GPUManager::getInstance().free_device_memory(cache.d_dst_block_idx); cache.d_dst_block_idx = nullptr; }
    if (cache.d_rank_offsets)  { GPUManager::getInstance().free_device_memory(cache.d_rank_offsets); cache.d_rank_offsets = nullptr; }
    if (cache.h_rank_offsets)  { delete[] cache.h_rank_offsets; cache.h_rank_offsets = nullptr; }
    cache.num_points = 0;
}

void ShellPatch::gpu_setup_phase2(int cpusize, MyList<Patch> *CPatL, int Symmetry,
                                  double Rrange[2], double DH_radial) {
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // --- Phase A: 收集 Cartesian block 元数据 ---
    int num_blocks = 0;
    MyList<Patch> *Pp = CPatL;
    while (Pp) { MyList<Block> *b = Pp->data->blb; while (b) { num_blocks++; b = b->next; } Pp = Pp->next; }
    if (num_blocks == 0) return;

    int *h_shape = new int[num_blocks * 3];
    double *h_bbox = new double[num_blocks * 6];
    double *h_CDH  = new double[num_blocks * 3];
    int *h_offsets = new int[num_blocks + 1];
    int total_pts = 0;
    h_offsets[0] = 0;
    int bi = 0;
    Pp = CPatL;
    while (Pp) {
        MyList<Block> *b = Pp->data->blb;
        while (b) {
            Block *Bg = b->data;
            for (int d = 0; d < 3; d++) {
                h_shape[bi * 3 + d] = Bg->shape[d];
                h_bbox[bi * 6 + d] = Bg->bbox[d];
                h_bbox[bi * 6 + 3 + d] = Bg->bbox[dim + d];
            }
            h_CDH[bi * 3 + 0] = Pp->data->getdX(0);
            h_CDH[bi * 3 + 1] = Pp->data->getdX(1);
            h_CDH[bi * 3 + 2] = Pp->data->getdX(2);
            total_pts += Bg->shape[0] * Bg->shape[1] * Bg->shape[2];
            bi++;
            h_offsets[bi] = total_pts;
            b = b->next;
        }
        Pp = Pp->next;
    }

    // --- Phase B: 上传元数据到 GPU ---
    cudaStream_t local_stream;
    cudaStreamCreate(&local_stream);

    // --- Phase B: 分配内存并使用 Local Stream 异步上传 ---
    int *d_offsets = GPUManager::getInstance().allocate_device_memory<int>(num_blocks + 1, local_stream);
    int *d_shape   = GPUManager::getInstance().allocate_device_memory<int>(num_blocks * 3, local_stream);
    double *d_bbox = GPUManager::getInstance().allocate_device_memory<double>(num_blocks * 6, local_stream);
    double *d_CDH  = GPUManager::getInstance().allocate_device_memory<double>(num_blocks * 3, local_stream);
    
    // 【修改】调用带 stream 参数的异步 API，确保与 Kernel 在同一流水线
    GPUManager::getInstance().sync_to_gpu(h_offsets, d_offsets, num_blocks + 1, local_stream);
    GPUManager::getInstance().sync_to_gpu(h_shape, d_shape, num_blocks * 3, local_stream);
    GPUManager::getInstance().sync_to_gpu(h_bbox, d_bbox, num_blocks * 6, local_stream);
    GPUManager::getInstance().sync_to_gpu(h_CDH, d_CDH, num_blocks * 3, local_stream);

    int max_output = total_pts; 
    int *d_count   = GPUManager::getInstance().allocate_device_memory<int>(1, local_stream);
    double *d_gpox = GPUManager::getInstance().allocate_device_memory<double>(max_output * 3, local_stream);
    double *d_lpox = GPUManager::getInstance().allocate_device_memory<double>(max_output * 3, local_stream);
    int    *d_sst  = GPUManager::getInstance().allocate_device_memory<int>(max_output, local_stream);
    uint64_t *d_sort_key = GPUManager::getInstance().allocate_device_memory<uint64_t>(max_output, local_stream);

    // --- Phase C: 启动 GPU 内核，传入 Local Stream ---
    gpu_setup_phase2_launch(
        local_stream, // 【修复】替换掉导致隐患的 0 (默认流)
        total_pts, d_offsets, d_shape, d_bbox, d_CDH,
        CS_width, Rrange[0], Rrange[1], DH_radial,
        d_count, d_gpox, d_lpox, d_sst, d_sort_key, max_output 
    );

    // --- Phase D: 异步下载结果并执行局部同步 ---
    int h_count = 0;
    // 注意：d_count 需要优先同步下来，因为紧接着 CPU 需要它来判断后续拷贝量
    // 这里对 count 单独做一次同步是必要的，或者安全起见直接在此处同步整个流
    GPUManager::getInstance().sync_to_cpu(&h_count, d_count, 1, local_stream);
    cudaStreamSynchronize(local_stream); // 【修改】细粒度流同步，替代暴力的 synchronize_all()

    if (h_count > max_output) h_count = max_output; 

    double *h_gpox = new double[h_count * 3];
    double *h_lpox = new double[h_count * 3];
    int    *h_sst  = new int[h_count];
    uint64_t *h_sort_key = new uint64_t[h_count]; 
    
    // 【修改】继续使用 local_stream 异步取回主要数据数组
    GPUManager::getInstance().sync_to_cpu(h_gpox, d_gpox, h_count * 3, local_stream);
    GPUManager::getInstance().sync_to_cpu(h_lpox, d_lpox, h_count * 3, local_stream);
    GPUManager::getInstance().sync_to_cpu(h_sst, d_sst, h_count, local_stream);
    GPUManager::getInstance().sync_to_cpu(h_sort_key, d_sort_key, h_count, local_stream);

    // 必须等待所有数据回传完毕再进行 CPU 排序
    cudaStreamSynchronize(local_stream);

    // 销毁局部流，释放资源
    cudaStreamDestroy(local_stream);

    // --- 强制拓扑排序 ---
    int *sort_indices = new int[h_count];
    std::iota(sort_indices, sort_indices + h_count, 0); 
    std::sort(sort_indices, sort_indices + h_count, [&](int a, int b) {
        return h_sort_key[a] < h_sort_key[b];
    });

    double *sorted_gpox = new double[h_count * 3];
    double *sorted_lpox = new double[h_count * 3];
    int    *sorted_sst  = new int[h_count];
    
    for (int i = 0; i < h_count; i++) {
        int src_idx = sort_indices[i];
        sorted_gpox[i * 3 + 0] = h_gpox[src_idx * 3 + 0];
        sorted_gpox[i * 3 + 1] = h_gpox[src_idx * 3 + 1];
        sorted_gpox[i * 3 + 2] = h_gpox[src_idx * 3 + 2];
        sorted_lpox[i * 3 + 0] = h_lpox[src_idx * 3 + 0];
        sorted_lpox[i * 3 + 1] = h_lpox[src_idx * 3 + 1];
        sorted_lpox[i * 3 + 2] = h_lpox[src_idx * 3 + 2];
        sorted_sst[i]          = h_sst[src_idx];
    }

    delete[] h_gpox; h_gpox = sorted_gpox;
    delete[] h_lpox; h_lpox = sorted_lpox;
    delete[] h_sst;  h_sst  = sorted_sst;
    delete[] sort_indices;
    delete[] h_sort_key;

    // --- Phase E: 重建链表并进行 Shell Block 搜索 ---
    // [Step 2] 恢复 `ble` 边界截断，防止脏块进入
    int num_sh_blocks = 0;
    MyList<ss_patch> *sPt = PatL;
    while (sPt) { 
        MyList<Block> *b = sPt->data->blb; 
        while (b) { 
            num_sh_blocks++; 
            if (b == sPt->data->ble) break; 
            b = b->next; 
        } 
        sPt = sPt->next; 
    }
    
    Block **sh_blks = new Block *[num_sh_blocks];
    int    *sh_sst  = new int[num_sh_blocks];
    double *sh_patch_bbox = new double[num_sh_blocks * 6]; 
    { 
        int si = 0; sPt = PatL;
        while (sPt) { 
            MyList<Block> *b = sPt->data->blb; 
            while (b) { 
                sh_blks[si] = b->data; 
                sh_sst[si] = sPt->data->sst; 
                for(int d=0; d<6; d++) sh_patch_bbox[si*6 + d] = sPt->data->bbox[d];
                si++; 
                if (b == sPt->data->ble) break; 
                b = b->next; 
            } 
            sPt = sPt->next; 
        } 
    }

    double DH[dim];
    for (int d = 0; d < dim; d++) DH[d] = getdX(d);

    MyList<pointstru> **csatc_heads = new MyList<pointstru>*[cpusize];
    MyList<pointstru> **csatc_tails = new MyList<pointstru>*[cpusize];
    for (int r = 0; r < cpusize; r++) {
        csatc_heads[r] = csatc_src[r];
        csatc_tails[r] = csatc_heads[r];
        while (csatc_tails[r] && csatc_tails[r]->next) csatc_tails[r] = csatc_tails[r]->next;
    }

    for (int p = 0; p < h_count; p++) {
        double gx = h_gpox[p * 3 + 0], gy = h_gpox[p * 3 + 1], gz = h_gpox[p * 3 + 2];
        double lx = h_lpox[p * 3 + 0], ly = h_lpox[p * 3 + 1], lz = h_lpox[p * 3 + 2];
        int sst = h_sst[p];

        // [Step 3] 修复 MPI 归属权优先级，寻找最小的 rank 块
        int min_rank = cpusize + 1; // 哨兵上限值
        Block* best_bg = nullptr;

        for (int sb = 0; sb < num_sh_blocks; sb++) {
            if (sh_sst[sb] == sst) {
                Block *Bg = sh_blks[sb];
                double *p_bbox = &sh_patch_bbox[sb * 6];
                double llb[3], uub[3];

                for (int d = 0; d < 2; d++) {
                    if (feq(Bg->bbox[d], -PI / 4 - overghost * DH[d], DH[d] / 2))
                        llb[d] = -PI / 4;
                    else if (feq(Bg->bbox[d], p_bbox[d], DH[d] / 2))
                        llb[d] = Bg->bbox[d];
#ifdef Vertex
                    else llb[d] = Bg->bbox[d] + (ghost_width - 1) * DH[d];
#else
                    else llb[d] = Bg->bbox[d] + ghost_width * DH[d];
#endif

                    if (feq(Bg->bbox[dim + d], PI / 4 + overghost * DH[d], DH[d] / 2))
                        uub[d] = PI / 4;
                    else if (feq(Bg->bbox[dim + d], p_bbox[dim + d], DH[d] / 2))
                        uub[d] = Bg->bbox[dim + d];
                    else uub[d] = Bg->bbox[dim + d] - ghost_width * DH[d];
                }

                if (feq(Bg->bbox[2], p_bbox[2], DH[2] / 2))
                    llb[2] = Bg->bbox[2];
#ifdef Vertex
                else llb[2] = Bg->bbox[2] + (ghost_width - 1) * DH[2];
#else
                else llb[2] = Bg->bbox[2] + ghost_width * DH[2];
#endif
                if (feq(Bg->bbox[5], p_bbox[5], DH[2] / 2))
                    uub[2] = Bg->bbox[5];
                else uub[2] = Bg->bbox[5] - ghost_width * DH[2];

                if (lx > llb[0] - 0.0001 * DH[0] && lx < uub[0] + 0.0001 * DH[0] &&
                    ly > llb[1] - 0.0001 * DH[1] && ly < uub[1] + 0.0001 * DH[1] &&
                    lz > llb[2] - 0.0001 * DH[2] && lz < uub[2] + 0.0001 * DH[2]) {
                    
                    // [Step 3] 记录当前能覆盖该点的最小 rank，不提前退出
                    if (Bg->rank < min_rank) {
                        min_rank = Bg->rank;
                        best_bg = Bg;
                    }
                }
            }
        }

        // 遍历完毕，把点指派给最小 rank 的进程（复刻 CPU 的行为）
        if (best_bg != nullptr) {
            MyList<pointstru> *ps = new MyList<pointstru>;
            ps->data = new pointstru;
            ps->data->Bg = best_bg;
            ps->data->gpox[0] = gx; ps->data->gpox[1] = gy; ps->data->gpox[2] = gz;
            ps->data->lpox[0] = lx; ps->data->lpox[1] = ly; ps->data->lpox[2] = lz;
            ps->data->ssst = sst;
            ps->data->tsst = -1;
            ps->data->dumyd = getdumydimension(ps->data->ssst, ps->data->tsst);
            ps->data->coef = 0; ps->data->sind = 0; ps->next = 0;

            if (csatc_tails[min_rank]) {
                csatc_tails[min_rank]->next = ps;
                csatc_tails[min_rank] = ps;
            } else {
                csatc_heads[min_rank] = csatc_tails[min_rank] = ps;
            }
        } else {
            // [补充] 补齐原 CPU 版本中因为找不到有效块产生的异常终止机制
            if (myrank == 0) {
                std::cout << "ShellPatch::prolongpointstru (GPU Phase) fail to find shell source point for" << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
    }
    for (int r = 0; r < cpusize; r++) csatc_src[r] = csatc_heads[r];

    // --- 清理 ---
    delete[] csatc_heads; delete[] csatc_tails;
    delete[] sh_blks; delete[] sh_sst; delete[] sh_patch_bbox;
    delete[] h_gpox; delete[] h_lpox; delete[] h_sst;
    delete[] h_shape; delete[] h_bbox; delete[] h_CDH; delete[] h_offsets;
    GPUManager::getInstance().free_device_memory(d_offsets);
    GPUManager::getInstance().free_device_memory(d_shape);
    GPUManager::getInstance().free_device_memory(d_bbox);
    GPUManager::getInstance().free_device_memory(d_CDH);
    GPUManager::getInstance().free_device_memory(d_count);
    GPUManager::getInstance().free_device_memory(d_gpox);
    GPUManager::getInstance().free_device_memory(d_lpox);
    GPUManager::getInstance().free_device_memory(d_sst);
    GPUManager::getInstance().free_device_memory(d_sort_key); 
}

#endif // USE_GPU
