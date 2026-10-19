
#ifndef GETNPEM2_H
#define GETNPEM2_H

#ifdef fortran1
#define f_getnpem2 getnpem2
#define f_getnpem2_point getnpem2_point
#define f_getnpem1_point getnpem1_point
#define f_getnpem2_ss getnpem2_ss
#define f_getnpem1 getnpem1
#define f_getnpem1_ss getnpem1_ss
#endif
#ifdef fortran2
#define f_getnpem2 GETNPEM2
#define f_getnpem2_point GETNPEM2_POINT
#define f_getnpem1_point GETNPEM1_POINT
#define f_getnpem2_ss GETNPEM2_SS
#define f_getnpem1 GETNPEM1
#define f_getnpem1_ss GETNPEM1_SS
#endif
#ifdef fortran3
#define f_getnpem2 getnpem2_
#define f_getnpem2_point getnpem2_point_
#define f_getnpem1_point getnpem1_point_
#define f_getnpem2_ss getnpem2_ss_
#define f_getnpem1 getnpem1_
#define f_getnpem1_ss getnpem1_ss_
#endif

extern "C"
{
        void f_getnpem2(int *, double *, double *, double *,
                        double *, double *, double *, double *, double *, double *, double *,
                        double *, double *, double *, double *, double *, double *,
                        double *, double *, int &);
}

extern "C"
{
        void f_getnpem2_point(double &, double &, double &,
                              double &, double &, double &, double &, double &, double &, double &,
                              double &, double &, double &, double &, double &, double &,
                              double &, double &);
}

extern "C"
{
        void f_getnpem2_ss(int *, double *, double *, double *, double *, double *, double *,
                           double *, double *, double *,
                           double *, double *, double *,
                           double *, double *, double *,
                           double *, double *, double *, double *, double *, double *,
                           double *, double *, double *, double *, double *, double *,
                           double *, double *, double *, double *, double *, double *,
                           double *, double *, double *, double *, double *, double *, double *,
                           double *, double *, double *, double *, double *, double *,
                           double *, double *, int &, int &);
}

extern "C"
{
        void f_getnpem1(int *, double *, double *, double *,
                        double *, double *, double *, double *, double *, double *, double *,
                        double *, double *, double *, double *, double *, double *,
                        double *, double *, int &);
}

extern "C"
{
        void f_getnpem1_ss(int *, double *, double *, double *, double *, double *, double *,
                           double *, double *, double *,
                           double *, double *, double *,
                           double *, double *, double *,
                           double *, double *, double *, double *, double *, double *,
                           double *, double *, double *, double *, double *, double *,
                           double *, double *, double *, double *, double *, double *,
                           double *, double *, double *, double *, double *, double *, double *,
                           double *, double *, double *, double *, double *, double *,
                           double *, double *, int &, int &);
}

extern "C"
{
        void f_getnpem1_point(double &, double &, double &,
                              double &, double &, double &, double &, double &, double &, double &,
                              double &, double &, double &, double &, double &, double &,
                              double &, double &);
}

#ifdef USE_GPU
// GPU launch wrappers for EM Newman-Penrose scalars
void gpu_getnpem2_launch(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_chi,
    const double* d_dxx, const double* d_gxy, const double* d_gxz,
    const double* d_dyy, const double* d_gyz, const double* d_dzz,
    const double* d_Ex, const double* d_Ey, const double* d_Ez,
    const double* d_Bx, const double* d_By, const double* d_Bz,
    double* d_Rphi2, double* d_Iphi2,
    int symmetry
);

void gpu_getnpem1_launch(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_chi,
    const double* d_dxx, const double* d_gxy, const double* d_gxz,
    const double* d_dyy, const double* d_gyz, const double* d_dzz,
    const double* d_Ex, const double* d_Ey, const double* d_Ez,
    const double* d_Bx, const double* d_By, const double* d_Bz,
    double* d_Rphi1, double* d_Iphi1,
    int symmetry
);

void gpu_getnpem2_ss_launch(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_chi,
    const double* d_dxx, const double* d_gxy, const double* d_gxz,
    const double* d_dyy, const double* d_gyz, const double* d_dzz,
    const double* d_Ex, const double* d_Ey, const double* d_Ez,
    const double* d_Bx, const double* d_By, const double* d_Bz,
    double* d_Rphi2, double* d_Iphi2,
    int symmetry
);

void gpu_getnpem1_ss_launch(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_chi,
    const double* d_dxx, const double* d_gxy, const double* d_gxz,
    const double* d_dyy, const double* d_gyz, const double* d_dzz,
    const double* d_Ex, const double* d_Ey, const double* d_Ez,
    const double* d_Bx, const double* d_By, const double* d_Bz,
    double* d_Rphi1, double* d_Iphi1,
    int symmetry
);
#endif /* USE_GPU */

#endif /* GETNPEM2_H */
