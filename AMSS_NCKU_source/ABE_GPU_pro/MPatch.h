
#ifndef PATCH_H
#define PATCH_H

#include <mpi.h>
#include "MyList.h"
#include "Block.h"
#include "var.h"
#include "macrodef.h" //need dim here; Vertex or Cell; ghost_width

#ifdef USE_GPU
#include <cuda_runtime.h>
#endif

class Patch
{

public:
    int lev;
    int shape[dim];
    double bbox[2 * dim]; // this bbox includes buffer points
    MyList<Block> *blb, *ble;
    int lli[dim], uui[dim]; // denote the buffer points on each boundary

public:
    Patch() {};
    Patch(int DIM, int *shapei, double *bboxi, int levi, bool buflog, int Symmetry);

    ~Patch();

    void checkPatch(bool buflog);
    void checkPatch(bool buflog, const int out_rank);
    void checkBlock();
    void Interp_Points(MyList<var> *VarList,
                             int NN, double **XX,
                             double *Shellf, int Symmetry);
    bool Interp_ONE_Point(MyList<var> *VarList, double *XX,
                                 double *Shellf, int Symmetry);
    double getdX(int dir);

    void Find_Maximum(MyList<var> *VarList, double *XX,
                            double *Shellf);

    bool Find_Point(double *XX);

    void Interp_Points(MyList<var> *VarList,
                             int NN, double **XX,
                             double *Shellf, int Symmetry, MPI_Comm Comm_here);
    bool Interp_ONE_Point(MyList<var> *VarList, double *XX,
                                 double *Shellf, int Symmetry, MPI_Comm Comm_here);
    void Find_Maximum(MyList<var> *VarList, double *XX,
                            double *Shellf, MPI_Comm Comm_here);
#ifdef USE_GPU
    void Interp_Points_GPU(
        cudaStream_t stream,
        MyList<var> *VarList,
        int NN, double **d_XX,
        double *d_Shellf, int Symmetry
    );
    void Interp_Points_GPU(
        cudaStream_t stream,
        MyList<var> *VarList,
        int NN, double **d_XX,
        double *d_Shellf, int Symmetry, MPI_Comm Comm_here);
    bool Interp_ONE_Point_GPU(
        cudaStream_t stream,
        MyList<var> *VarList,
        double *d_XX_0, double *d_XX_1, double *d_XX_2,
        double *d_Shellf, int Symmetry);
    bool Interp_ONE_Point_GPU(
        cudaStream_t stream,
        MyList<var> *VarList,
        double *d_XX_0, double *d_XX_1, double *d_XX_2,
        double *d_Shellf, int Symmetry, MPI_Comm Comm_here);
    void Find_Maximum_GPU(
        cudaStream_t stream,
        MyList<var> *VarList,
        double *d_XX, double *d_Shellf);
    void Find_Maximum_GPU(
        cudaStream_t stream,
        MyList<var> *VarList,
        double *d_XX, double *d_Shellf, MPI_Comm Comm_here);
    bool Interp_N_Points_GPU(
        MyList<var> *VarList,
        int NN, double *d_XX_0, double *d_XX_1, double *d_XX_2,
        double *d_shellf, int *d_weight, int Symmetry
    );
#endif /* USE_GPU */
};

#ifdef USE_GPU
void gpu_normalize_shellf_launch(
    cudaStream_t stream,
    int NN, int num_var,
    double* d_Shellf,
    const int* d_Weight,
    int* d_err_flag
);

void gpu_find_maximum_launch(
    cudaStream_t stream,
    int shape_0, int shape_1, int shape_2,
    int pshape_0, int pshape_1, int pshape_2,
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_fun,
    int llb_0, int llb_1, int llb_2,
    int uub_0, int uub_1, int uub_2,
    double* d_max_val, double* d_max_pos_x, double* d_max_pos_y, double* d_max_pos_z
);

#endif /* USE_GPU */

#endif /* PATCH_H */
