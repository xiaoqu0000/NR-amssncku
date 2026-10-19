#ifndef MPI_ERROR_BATCH_H
#define MPI_ERROR_BATCH_H

#include <mpi.h>
#include <cstdlib>
#include <cstring>
#include <iostream>

// NaN / fatal-error guard. By default `mark(...)` does a local check and
// immediately calls MPI_Abort if it sees a non-zero value — MPI_Abort on
// MPI_COMM_WORLD cascades to every rank, so the previous batched Allreduce
// pattern is not actually needed for correctness.
//
// Set AMSS_NAN_ABORT_IMMEDIATE=0 to fall back to the old behavior: mark()
// only ORs into a local flag, flush() does an MPI_Allreduce every N steps
// (N = AMSS_NAN_CHECK_EVERY, default 8), and the caller decides whether to
// abort. Use this if some site needs all ranks to observe an error tag
// before aborting (none do today, but it's there as an escape hatch).
class MpiErrorBatch {
public:
    enum Tag {
        TAG_BSSN_NAN,
        TAG_Z4C_NAN,
        TAG_SCALAR_NAN,
        TAG_EM_NAN,
        TAG_ESCALAR_NAN,
        TAG_NEWTON_FAIL,
        TAG_TVF,
        TAG_CONSTRAINT,
        TAG_BOUNDARY,
        N_TAGS
    };

    void mark(Tag t, int val) {
        if (t < 0 || t >= N_TAGS) return;
        if (immediate_abort()) {
            if (val != 0) {
                int rank = 0;
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                std::cerr << "[MpiErrorBatch] rank " << rank
                          << " hit tag " << static_cast<int>(t)
                          << " (val=" << val << "); aborting MPI_COMM_WORLD\n";
                std::cerr.flush();
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            return;
        }
        local_flags_[t] |= val;
    }

    void flush(MPI_Comm comm = MPI_COMM_WORLD) {
        if (immediate_abort()) return;  // already aborted on the spot if needed
        MPI_Allreduce(local_flags_, global_flags_, N_TAGS, MPI_INT,
                      MPI_LOR, comm);
    }

    bool any() const {
        if (immediate_abort()) return false;
        for (int i = 0; i < N_TAGS; i++)
            if (global_flags_[i]) return true;
        return false;
    }

    bool any(Tag t) const {
        if (immediate_abort()) return false;
        if (t >= 0 && t < N_TAGS)
            return global_flags_[t] != 0;
        return false;
    }

    void reset() {
        memset(local_flags_, 0, sizeof(local_flags_));
        memset(global_flags_, 0, sizeof(global_flags_));
    }

    static MpiErrorBatch& instance() {
        static MpiErrorBatch inst;
        return inst;
    }

    static int check_every() {
        static int n = -1;
        if (n < 0) {
            const char* env = std::getenv("AMSS_NAN_CHECK_EVERY");
            n = env ? std::atoi(env) : 8;
            if (n <= 0) n = 8;
        }
        return n;
    }

private:
    MpiErrorBatch() { reset(); }

    // Default ON: any rank seeing a non-zero error immediately MPI_Aborts —
    // mpirun then kills every other rank, so there's no need for a collective.
    static bool immediate_abort() {
        static int cached = -1;
        if (cached < 0) {
            const char* env = std::getenv("AMSS_NAN_ABORT_IMMEDIATE");
            cached = (env && env[0] == '0') ? 0 : 1;
        }
        return cached != 0;
    }

    int local_flags_[N_TAGS];
    int global_flags_[N_TAGS];
};

#endif
