#include "interp_shm.h"

#include <hwloc.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cfloat>
#include <algorithm>

// ---------------------------------------------------------------------------
// Constructor / Destructor
// ---------------------------------------------------------------------------

InterpSHM::InterpSHM()
    : numa_comm(MPI_COMM_NULL),
      inter_numa_comm(MPI_COMM_NULL),
      numa_rank(-1), numa_size(0),
      numa_id(-1),
      world_rank(-1), world_size(0),
      is_leader(false),
      win_fgfs(MPI_WIN_NULL),
      win_shellf(MPI_WIN_NULL),
      win_weight(MPI_WIN_NULL),
      fgfs_shm(nullptr),
      shellf_shm(nullptr),
      weight_shm(nullptr),
      fgfs_shm_size(0),
      shellf_shm_max(0),
      weight_shm_max(0),
      initialized(false)
{
}

InterpSHM::~InterpSHM()
{
    if (initialized) finalize();
}

// ---------------------------------------------------------------------------
// detect_numa  —  use hwloc to find which NUMA node this process sits on
// ---------------------------------------------------------------------------

void InterpSHM::detect_numa()
{
    hwloc_topology_t topo;
    hwloc_topology_init(&topo);
    hwloc_topology_load(topo);

    // Get cpuset of calling thread
    hwloc_cpuset_t cpuset = hwloc_bitmap_alloc();
    hwloc_get_cpubind(topo, cpuset, HWLOC_CPUBIND_PROCESS);

    // Walk NUMA nodes, find the one whose cpuset intersects ours
    int n_numa = hwloc_get_nbobjs_by_type(topo, HWLOC_OBJ_NUMANODE);
    numa_id = 0; // fallback
    for (int i = 0; i < n_numa; i++) {
        hwloc_obj_t obj = hwloc_get_obj_by_type(topo, HWLOC_OBJ_NUMANODE, i);
        if (obj && obj->cpuset && hwloc_bitmap_intersects(cpuset, obj->cpuset)) {
            numa_id = i;
            break;
        }
    }

    hwloc_bitmap_free(cpuset);
    hwloc_topology_destroy(topo);
}

// ---------------------------------------------------------------------------
// init  —  create communicators and allocate shared-memory windows
// ---------------------------------------------------------------------------

void InterpSHM::init(int max_nn_total, int max_num_var)
{
    if (initialized) return;

    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // 1. Detect NUMA node
    detect_numa();

    // 2. Split WORLD by numa_id → numa_comm
    MPI_Comm_split(MPI_COMM_WORLD, numa_id, world_rank, &numa_comm);
    MPI_Comm_rank(numa_comm, &numa_rank);
    MPI_Comm_size(numa_comm, &numa_size);

    // 3. Create inter-NUMA communicator (leaders only = rank 0 of each numa_comm)
    is_leader = (numa_rank == 0);
    int color_inter = is_leader ? 0 : MPI_UNDEFINED;
    MPI_Comm_split(MPI_COMM_WORLD, color_inter, world_rank, &inter_numa_comm);
    // inter_numa_comm is MPI_COMM_NULL on non-leaders

    // 4. Shared-memory window for fgfs (field + coordinate data)
    //    Default 32 MB on rank 0, others query the pointer.
    // fgfs_shm_size = 128UL * 1024 * 1024; // 128 MB
    // 修改为更大的值 xiaoqu 2026/07/08
    fgfs_shm_size = 1024UL * 1024 * 1024; // 1024 MB
    {
        MPI_Aint local_sz = (numa_rank == 0)
                                ? (MPI_Aint)fgfs_shm_size
                                : 0;
        int disp_unit = sizeof(double);
        MPI_Win_allocate_shared(local_sz, disp_unit, MPI_INFO_NULL,
                                numa_comm, &fgfs_shm, &win_fgfs);
        if (numa_rank != 0) {
            // Query the base pointer from rank 0
            MPI_Aint sz_query;
            int disp_query;
            MPI_Win_shared_query(win_fgfs, 0, &sz_query, &disp_query,
                                 &fgfs_shm);
        }
    }

    // 5. Shared-memory window for shellf (interpolation output accumulator)
    shellf_shm_max = (size_t)max_nn_total * max_num_var;
    {
        MPI_Aint local_sz = (numa_rank == 0)
                                ? (MPI_Aint)(shellf_shm_max * sizeof(double))
                                : 0;
        int disp_unit = sizeof(double);
        MPI_Win_allocate_shared(local_sz, disp_unit, MPI_INFO_NULL,
                                numa_comm, &shellf_shm, &win_shellf);
        if (numa_rank != 0) {
            MPI_Aint sz_query;
            int disp_query;
            MPI_Win_shared_query(win_shellf, 0, &sz_query, &disp_query,
                                 &shellf_shm);
        }
    }

    // 6. Shared-memory window for weight (overlap counter)
    weight_shm_max = max_nn_total;
    {
        MPI_Aint local_sz = (numa_rank == 0)
                                ? (MPI_Aint)(weight_shm_max * (int)sizeof(int))
                                : 0;
        int disp_unit = sizeof(int);
        MPI_Win_allocate_shared(local_sz, disp_unit, MPI_INFO_NULL,
                                numa_comm, &weight_shm, &win_weight);
        if (numa_rank != 0) {
            MPI_Aint sz_query;
            int disp_query;
            MPI_Win_shared_query(win_weight, 0, &sz_query, &disp_query,
                                 &weight_shm);
        }
    }

    initialized = true;
}

// ---------------------------------------------------------------------------
// gather  —  collect block data into shared memory for NUMA-local access
//
// Returns: number of blocks gathered (size of block_table)
// ---------------------------------------------------------------------------

int InterpSHM::gather(MyList<Block> *blb, MyList<Block> *ble,
                      int *var_sgfns, int num_vars,
                      double **XX, int NN)
{
    block_table.clear();

    // Bounds check: shellf and weight windows
    if ((size_t)NN * num_vars > shellf_shm_max) {
        if (world_rank == 0)
            fprintf(stderr, "InterpSHM::gather: shellf overflow! Need %zu, have %zu\n",
                    (size_t)NN * num_vars, shellf_shm_max);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (NN > weight_shm_max) {
        if (world_rank == 0)
            fprintf(stderr, "InterpSHM::gather: weight overflow! Need %d, have %d\n",
                    NN, weight_shm_max);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // ------------------------------------------------------------------
    // 1. Compute bounding box of all interpolation points
    // ------------------------------------------------------------------
    double pt_lo[dim], pt_hi[dim];
    for (int d = 0; d < dim; d++) {
        pt_lo[d] =  DBL_MAX;
        pt_hi[d] = -DBL_MAX;
    }
    for (int j = 0; j < NN; j++) {
        for (int d = 0; d < dim; d++) {
            if (XX[d][j] < pt_lo[d]) pt_lo[d] = XX[d][j];
            if (XX[d][j] > pt_hi[d]) pt_hi[d] = XX[d][j];
        }
    }

    // ------------------------------------------------------------------
    // 2. Build world_rank → numa_id mapping via MPI_Allgather
    //    (16 ints, negligible cost)
    // ------------------------------------------------------------------
    std::vector<int> rank_numa(world_size);
    MPI_Allgather(&numa_id, 1, MPI_INT,
                  rank_numa.data(), 1, MPI_INT, MPI_COMM_WORLD);

    // ------------------------------------------------------------------
    // 3. Scan block list, find overlapping blocks in our NUMA group
    // ------------------------------------------------------------------
    int block_idx = 0;
    MyList<Block> *Bp = blb;
    while (Bp) {
        Block *B = Bp->data;

        // Check if this block's owner is in our NUMA group
        if (rank_numa[B->rank] == numa_id) {
            // Check bbox overlap with point bounding box
            bool overlap = true;
            for (int d = 0; d < dim; d++) {
                if (B->bbox[d] > pt_hi[d] || B->bbox[dim + d] < pt_lo[d]) {
                    overlap = false;
                    break;
                }
            }
            if (overlap) {
                ShmBlockInfo info;
                info.block_id   = block_idx;
                info.owner_rank = B->rank;
                info.nn = 1;
                for (int d = 0; d < dim; d++) {
                    info.shape[d] = B->shape[d];
                    info.bbox[d]  = B->bbox[d];
                    info.bbox[dim + d] = B->bbox[dim + d];
                    info.nn *= B->shape[d];
                }
                // Offsets will be computed below
                for (int d = 0; d < dim; d++) info.X_offset[d] = 0;
                info.fgfs_offset = 0;
                block_table.push_back(info);
            }
        }

        if (Bp == ble) break;
        Bp = Bp->next;
        block_idx++;
    }

    // ------------------------------------------------------------------
    // 4. Compute shared buffer layout (cache-line aligned to 64 bytes)
    // ------------------------------------------------------------------
    size_t offset = 0;
    for (size_t i = 0; i < block_table.size(); i++) {
        ShmBlockInfo &info = block_table[i];
        size_t nn_bytes = (size_t)info.nn * sizeof(double);

        // X[0], X[1], X[2] — each has shape[d] elements (NOT nn)
        for (int d = 0; d < dim; d++) {
            info.X_offset[d] = offset;
            offset += (size_t)info.shape[d] * sizeof(double);
            offset = (offset + 63) & ~(size_t)63; // align to 64 bytes
        }

        // fgfs: num_vars fields, each nn doubles
        info.fgfs_offset = offset;
        offset += (size_t)num_vars * nn_bytes;
        offset = (offset + 63) & ~(size_t)63; // align to 64 bytes
    }

    // ------------------------------------------------------------------
    // 5. Check buffer size
    // ------------------------------------------------------------------
    if (offset > fgfs_shm_size) {
        if (world_rank == 0) {
            fprintf(stderr,
                    "InterpSHM::gather: shared buffer too small! "
                    "Need %zu bytes, have %zu bytes.\n",
                    offset, fgfs_shm_size);
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // ------------------------------------------------------------------
    // 6. Leader zeros the buffer, then barrier
    // ------------------------------------------------------------------
    if (numa_rank == 0) {
        memset(fgfs_shm, 0, offset);
    }
    MPI_Barrier(numa_comm);

    // ------------------------------------------------------------------
    // 7. Each owning rank copies its X[] and fgfs[var_sgfns[v]] into SHM
    // ------------------------------------------------------------------
    char *base = (char *)fgfs_shm;

    for (size_t i = 0; i < block_table.size(); i++) {
        const ShmBlockInfo &info = block_table[i];

        // Only the owner copies data
        if (info.owner_rank != world_rank) continue;

        // Find the actual Block pointer by walking the list
        MyList<Block> *Bq = blb;
        int idx = 0;
        while (Bq) {
            if (idx == info.block_id) break;
            if (Bq == ble) { Bq = nullptr; break; }
            Bq = Bq->next;
            idx++;
        }
        if (!Bq) continue;
        Block *B = Bq->data;

        size_t nn_bytes = (size_t)info.nn * sizeof(double);

        // Copy coordinate arrays X[d] — shape[d] elements each
        for (int d = 0; d < dim; d++) {
            memcpy(base + info.X_offset[d], B->X[d],
                   (size_t)info.shape[d] * sizeof(double));
        }

        // Copy field data for each requested variable
        for (int v = 0; v < num_vars; v++) {
            size_t dst_off = info.fgfs_offset + (size_t)v * info.nn * sizeof(double);
            memcpy(base + dst_off, B->fgfs[var_sgfns[v]], nn_bytes);
        }
    }

    // ------------------------------------------------------------------
    // 8. Barrier — ensure all copies complete before interpolation begins
    // ------------------------------------------------------------------
    MPI_Barrier(numa_comm);

    // ------------------------------------------------------------------
    // 9. Each rank zeros its slice of shellf_shm and weight_shm
    //    Distribute the zeroing across NUMA ranks for speed.
    // ------------------------------------------------------------------
    {
        size_t total_shellf = (size_t)NN * num_vars;
        size_t chunk = (total_shellf + numa_size - 1) / numa_size;
        size_t my_start = (size_t)numa_rank * chunk;
        size_t my_end   = my_start + chunk;
        if (my_end > total_shellf) my_end = total_shellf;
        if (my_start < total_shellf) {
            memset(shellf_shm + my_start, 0,
                   (my_end - my_start) * sizeof(double));
        }

        int w_chunk = (NN + numa_size - 1) / numa_size;
        int w_start = numa_rank * w_chunk;
        int w_end   = w_start + w_chunk;
        if (w_end > NN) w_end = NN;
        if (w_start < NN) {
            memset(weight_shm + w_start, 0,
                   (w_end - w_start) * sizeof(int));
        }
    }

    // ------------------------------------------------------------------
    // 10. Barrier — zeroing complete
    // ------------------------------------------------------------------
    MPI_Barrier(numa_comm);

    return (int)block_table.size();
}

// ---------------------------------------------------------------------------
// reduce  —  reduce shellf/weight across NUMA groups, then broadcast within
// ---------------------------------------------------------------------------

void InterpSHM::reduce(double *Shellf, int *Weight, int NN, int num_var)
{
    // 1. Barrier — ensure all ranks have finished writing to shellf_shm/weight_shm
    MPI_Barrier(numa_comm);

    // 2. Leaders: reduce across NUMA groups
    if (is_leader) {
        MPI_Allreduce(shellf_shm, Shellf,
                      NN * num_var, MPI_DOUBLE, MPI_SUM,
                      inter_numa_comm);
        MPI_Allreduce(weight_shm, Weight,
                      NN, MPI_INT, MPI_SUM,
                      inter_numa_comm);
    }

    // 3. Broadcast reduced results from leader to all ranks in NUMA group
    MPI_Bcast(Shellf, NN * num_var, MPI_DOUBLE, 0, numa_comm);
    MPI_Bcast(Weight, NN, MPI_INT, 0, numa_comm);
}

// ---------------------------------------------------------------------------
// finalize  —  free windows and communicators
// ---------------------------------------------------------------------------

void InterpSHM::finalize()
{
    if (!initialized) return;

    if (win_weight != MPI_WIN_NULL) MPI_Win_free(&win_weight);
    if (win_shellf != MPI_WIN_NULL) MPI_Win_free(&win_shellf);
    if (win_fgfs   != MPI_WIN_NULL) MPI_Win_free(&win_fgfs);

    if (inter_numa_comm != MPI_COMM_NULL) MPI_Comm_free(&inter_numa_comm);
    if (numa_comm       != MPI_COMM_NULL) MPI_Comm_free(&numa_comm);

    fgfs_shm   = nullptr;
    shellf_shm = nullptr;
    weight_shm = nullptr;
    win_fgfs   = MPI_WIN_NULL;
    win_shellf = MPI_WIN_NULL;
    win_weight = MPI_WIN_NULL;
    numa_comm  = MPI_COMM_NULL;
    inter_numa_comm = MPI_COMM_NULL;

    initialized = false;
}
