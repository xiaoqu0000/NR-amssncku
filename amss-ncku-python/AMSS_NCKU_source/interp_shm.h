#ifndef INTERP_SHM_H
#define INTERP_SHM_H

#include <mpi.h>
#include <cstring>
#include <vector>
#include "macrodef.h"
#include "MyList.h"
#include "var.h"
#include "Block.h"

struct ShmBlockInfo {
    int block_id;           // index in block list
    int owner_rank;         // MPI_COMM_WORLD rank
    int shape[dim];
    double bbox[2 * dim];
    int nn;                 // shape[0]*shape[1]*shape[2]
    size_t X_offset[dim];   // byte offset into fgfs_shm for X[0], X[1], X[2]
    size_t fgfs_offset;     // byte offset into fgfs_shm for fgfs data start
    // fgfs for var k at: fgfs_offset + k * nn * sizeof(double)
};

class InterpSHM {
public:
    MPI_Comm numa_comm;
    MPI_Comm inter_numa_comm;   // MPI_COMM_NULL on non-leaders
    int numa_rank, numa_size;
    int numa_id;
    int world_rank, world_size;
    bool is_leader;

    MPI_Win win_fgfs;
    MPI_Win win_shellf;
    MPI_Win win_weight;

    double *fgfs_shm;
    double *shellf_shm;
    int    *weight_shm;

    size_t fgfs_shm_size;
    size_t shellf_shm_max;
    int    weight_shm_max;

    std::vector<ShmBlockInfo> block_table;

    InterpSHM();
    ~InterpSHM();
    void init(int max_nn_total, int max_num_var);
    void finalize();

    int gather(MyList<Block> *blb, MyList<Block> *ble,
               int *var_sgfns, int num_vars,
               double **XX, int NN);

    void reduce(double *Shellf, int *Weight, int NN, int num_var);

private:
    bool initialized;
    void detect_numa();
};

#endif
