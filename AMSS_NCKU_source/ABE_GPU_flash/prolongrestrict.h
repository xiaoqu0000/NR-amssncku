
#ifndef PROLONGRESTRICT_H
#define PROLONGRESTRICT_H

#include "macrodef.fh"

#ifdef fortran1
#define f_prolong3 prolong3
#define f_prolongmix3 prolongmix3
#define f_prolongcopy3 prolongcopy3
#define f_restrict3 restrict3
#endif

#ifdef fortran2
#define f_prolong3 PROLONG3
#define f_prolongmix3 PROLONGMIX3
#define f_prolongcopy3 PROLONGCOPY3
#define f_restrict3 RESTRICT3
#endif

#ifdef fortran3
#define f_prolong3 prolong3_
#define f_prolongmix3 prolongmix3_
#define f_prolongcopy3 prolongcopy3_
#define f_restrict3 restrict3_
#endif

#include <cuda_runtime.h>

extern "C"
{
	int f_prolong3(int &, double *, double *, int *, double *,
				   double *, double *, int *, double *,
				   double *, double *, double *, int &);
}

extern "C"
{
	void f_restrict3(int &, double *, double *, int *, double *,
					 double *, double *, int *, double *,
					 double *, double *, double *, int &);
}

extern "C"
{
	int f_prolongmix3(int &, double *, double *, int *, double *,
					  double *, double *, int *, double *,
					  double *, double *, double *, int &,
					  double *, double *);
}

extern "C"
{
	int f_prolongcopy3(int &, double *, double *, int *, double *,
					   double *, double *, int *, double *,
					   double *, double *, double *, int &);
}

__device__ int d_idint(double a);

#ifdef Cell
__device__ void d_prolong3_device(
    int i, int j, int k,
    const double* llbc, const double* uubc, const int* extc, const double* func,
    const double* llbf, const double* uubf, const int* extf, double* funf, 
    const double* llbp, const double* uubp,
    const double* SoA, int Symmetry
);

__device__ void d_restrict3_device(
    int i, int j, int k, 
    const double* llbc, const double* uubc, const int* extc, double* func, // func is output
    const double* llbf, const double* uubf, const int* extf, const double* funf, // funf is input
    const double* llbr, const double* uubr,
    const double* SoA, int Symmetry
);
#endif  // Cell

void gpu_prolong3_launch(
    cudaStream_t stream,
    const double* d_src_c, double* d_dst_f,
    const double* llbc, const double* uubc, const int* extc,
    const double* llbf, const double* uubf, const int* extf,
    const double* llbt, const double* uubt,
    const double* SoA, int Symmetry
);

struct GpuStencilTransfer {
    const double* source_field;
    double* target_field;
    int ni;
    int nj;
    int nk;
    int i_start;
    int j_start;
    int k_start;
    int lbc0;
    int lbc1;
    int lbc2;
    int lbf0;
    int lbf1;
    int lbf2;
    int extc0;
    int extc1;
    int extc2;
    int extf0;
    int extf1;
    int extf2;
    double SoA0;
    double SoA1;
    double SoA2;
    int Symmetry;
};

bool gpu_build_prolong3_transfer(
    const double* d_src_c, double* d_dst_f,
    const double* llbc, const double* uubc, const int* extc,
    const double* llbf, const double* uubf, const int* extf,
    const double* llbt, const double* uubt,
    const double* SoA, int Symmetry,
    GpuStencilTransfer* transfer
);

bool gpu_build_restrict3_transfer(
    const double* d_src_f, double* d_dst_c,
    const double* llbc, const double* uubc, const int* extc,
    const double* llbf, const double* uubf, const int* extf,
    const double* llbt, const double* uubt,
    const double* SoA, int Symmetry,
    GpuStencilTransfer* transfer
);

void gpu_prolong3_batch_launch(
    cudaStream_t stream,
    const GpuStencilTransfer* transfers,
    int transfer_count
);

void gpu_restrict3_launch(
    cudaStream_t stream,
    const double* d_src_f, double* d_dst_c,
    const double* llbc, const double* uubc, const int* extc,
    const double* llbf, const double* uubf, const int* extf,
    const double* llbt, const double* uubt,
    const double* SoA, int Symmetry
);

void gpu_restrict3_batch_launch(
    cudaStream_t stream,
    const GpuStencilTransfer* transfers,
    int transfer_count
);

#endif /* PROLONGRESTRICT_H */
