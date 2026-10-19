#ifndef AMSS_REDUCE_HELPER_H
#define AMSS_REDUCE_HELPER_H

#include <mpi.h>
#include <vector>
#include <cstdlib>

// Reserved MPI tag for manual reduce. Picked to avoid collision with the
// data-transfer tag (1) used by gpu_transfer / intertransfer_GPU.
#ifndef AMSS_REDUCE_TAG
#define AMSS_REDUCE_TAG 0xA12D
#endif

namespace amss_reduce {

template<typename T> struct MPIType;
template<> struct MPIType<double> { static MPI_Datatype get() { return MPI_DOUBLE; } };
template<> struct MPIType<int>    { static MPI_Datatype get() { return MPI_INT; } };
template<> struct MPIType<float>  { static MPI_Datatype get() { return MPI_FLOAT; } };

// Gather + local sum to root. After return, root's `buf` holds the global sum;
// non-root processes' `buf` is unchanged (still their local partial). Use when
// only root reads the result (monitor file writes, etc.).
//
// Why this exists: on small N (<= 4 GPU processes) and small messages
// (NN ~ 7-21 doubles), MPI_Allreduce / MPI_Reduce protocol overhead dominates
// the actual data transfer. Hand-rolled point-to-point gather + local add is
// measurably faster — the user observed 2x-plus speedup at this scale.
template<typename T>
void gather_sum_to_root(T* buf, int count, int root, MPI_Comm comm) {
    int myrank, cpusize;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &cpusize);
    if (cpusize <= 1 || count <= 0) return;

    if (myrank == root) {
        std::vector<T> recv(static_cast<size_t>(cpusize - 1) * count);
        std::vector<MPI_Request> reqs;
        reqs.reserve(cpusize - 1);
        int idx = 0;
        for (int r = 0; r < cpusize; r++) {
            if (r == root) continue;
            MPI_Request req;
            MPI_Irecv(recv.data() + static_cast<size_t>(idx) * count,
                      count, MPIType<T>::get(), r, AMSS_REDUCE_TAG, comm, &req);
            reqs.push_back(req);
            idx++;
        }
        MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);
        // local accumulation
        for (int r = 0; r < cpusize - 1; r++) {
            const T* src = recv.data() + static_cast<size_t>(r) * count;
            for (int i = 0; i < count; i++) buf[i] += src[i];
        }
    } else {
        MPI_Send(buf, count, MPIType<T>::get(), root, AMSS_REDUCE_TAG, comm);
    }
}

// AMSS_MANUAL_REDUCE: env switch (default ON). Set to "0" to fall back to
// MPI_Allreduce — useful for A/B testing or for backends where the MPI
// implementation has a fast small-msg Reduce.
inline bool manual_reduce_enabled() {
    static int cached = -1;
    if (cached < 0) {
        const char* env = std::getenv("AMSS_MANUAL_REDUCE");
        cached = (env && env[0] == '0') ? 0 : 1;
    }
    return cached != 0;
}

// Drop-in replacement for `MPI_Allreduce(MPI_IN_PLACE, buf, count, type, MPI_SUM, comm)`
// for the case where only `root` reads the result. Falls back to MPI_Allreduce when
// the env switch disables the manual path, so the caller's later read on rank 0
// always sees a valid sum.
template<typename T>
void reduce_sum_to_root_inplace(T* buf, int count, int root, MPI_Comm comm) {
    if (manual_reduce_enabled()) {
        gather_sum_to_root(buf, count, root, comm);
    } else {
        // legacy path — keeps the result on every rank
        MPI_Allreduce(MPI_IN_PLACE, buf, count, MPIType<T>::get(), MPI_SUM, comm);
    }
}

// Manual point-to-point broadcast: root sends to every other rank, others receive.
// Mirrors MPI_Bcast's contract but uses one Send per peer. Worth it at small N for
// the same reason gather_sum_to_root beats Allreduce — protocol overhead per op
// dominates the wire cost when message is small and N <= ~4.
template<typename T>
void bcast_from_root(T* buf, int count, int root, MPI_Comm comm) {
    int myrank, cpusize;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &cpusize);
    if (cpusize <= 1 || count <= 0) return;

    if (myrank == root) {
        std::vector<MPI_Request> reqs;
        reqs.reserve(cpusize - 1);
        for (int r = 0; r < cpusize; r++) {
            if (r == root) continue;
            MPI_Request req;
            MPI_Isend(buf, count, MPIType<T>::get(), r,
                      AMSS_REDUCE_TAG, comm, &req);
            reqs.push_back(req);
        }
        MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(),
                    MPI_STATUSES_IGNORE);
    } else {
        MPI_Recv(buf, count, MPIType<T>::get(), root,
                 AMSS_REDUCE_TAG, comm, MPI_STATUS_IGNORE);
    }
}

// Drop-in replacement for `MPI_Allreduce(MPI_IN_PLACE, buf, count, type, MPI_SUM, comm)`
// when ALL ranks need the result. Runs gather-sum-to-root then bcast-from-root.
// Two passes through the same protocol overhead is still cheaper than the
// log-N tree of a real Allreduce when N and message size are small.
template<typename T>
void allreduce_sum_inplace(T* buf, int count, MPI_Comm comm) {
    if (manual_reduce_enabled()) {
        gather_sum_to_root(buf, count, 0, comm);
        bcast_from_root(buf, count, 0, comm);
    } else {
        MPI_Allreduce(MPI_IN_PLACE, buf, count, MPIType<T>::get(), MPI_SUM, comm);
    }
}

} // namespace amss_reduce

#endif // AMSS_REDUCE_HELPER_H
