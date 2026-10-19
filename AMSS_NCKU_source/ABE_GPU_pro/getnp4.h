
#ifndef GETNP4_H
#define GETNP4_H

#ifdef fortran1
#define f_getnp4old getnp4old
#define f_getnp4oldscalar getnp4oldscalar
#define f_getnp4oldscalar_ss getnp4oldscalar_ss
#define f_getnp4 getnp4
#define f_getnp4_point getnp4_point
#define f_getnp4_ss getnp4_ss
#define f_getnp4old_ss getnp4old_ss
#define f_getnp4scalar getnp4scalar
#define f_getnp4scalar_ss getnp4scalar_ss
#endif
#ifdef fortran2
#define f_getnp4 GETNP4
#define f_getnp4_point GETNP4_POINT
#define f_getnp4 GETNP4OLD
#define f_getnp4scalar GETNP4OLDSCALAR
#define f_getnp4_ss GETNP4_SS
#define f_getnp4old_ss GETNP4OLD_SS
#define f_getnp4oldscalar_ss GETNP4OLDSCALAR_SS
#define f_getnp4scalar GETNP4SCALAR
#define f_getnp4scalar_ss GETNP4SCALAR_SS
#endif
#ifdef fortran3
#define f_getnp4old getnp4old_
#define f_getnp4_point getnp4_point_
#define f_getnp4oldscalar getnp4oldscalar_
#define f_getnp4oldscalar_ss getnp4oldscalar_ss_
#define f_getnp4 getnp4_
#define f_getnp4_ss getnp4_ss_
#define f_getnp4old_ss getnp4old_ss_
#define f_getnp4scalar getnp4scalar_
#define f_getnp4scalar_ss getnp4scalar_ss_
#endif

extern "C"
{
        void f_getnp4old(int *, double *, double *, double *,
                         double *, double *,
                         double *, double *, double *, double *, double *, double *,
                         double *, double *, double *, double *, double *, double *,
                         double *, double *, double *,
                         double *, double *, double *, double *,
                         double *, double *, int &);
}

extern "C"
{
        void f_getnp4old_ss(int *, double *, double *, double *, double *, double *, double *,
                            double *, double *, double *,
                            double *, double *, double *,
                            double *, double *, double *,
                            double *, double *, double *, double *, double *, double *,
                            double *, double *, double *, double *, double *, double *,
                            double *, double *, double *, double *, double *, double *,
                            double *, double *,
                            double *, double *, double *, double *, double *, double *,
                            double *, double *, double *, double *, double *, double *,
                            double *, double *, double *,
                            double *, double *, double *, double *,
                            double *, double *, int &, int &);
}

extern "C"
{
        void f_getnp4oldscalar(int *, double *, double *, double *,
                               double *, double *, double *,
                               double *, double *, double *, double *, double *, double *,
                               double *, double *, double *, double *, double *, double *,
                               double *, double *, double *,
                               double *, double *, double *, double *,
                               double *, double *, int &);
}

extern "C"
{
        void f_getnp4oldscalar_ss(int *, double *, double *, double *, double *, double *, double *,
                                  double *, double *, double *,
                                  double *, double *, double *,
                                  double *, double *, double *,
                                  double *, double *, double *, double *, double *, double *,
                                  double *, double *, double *, double *, double *, double *,
                                  double *, double *, double *, double *, double *, double *,
                                  double *, double *, double *,
                                  double *, double *, double *, double *, double *, double *,
                                  double *, double *, double *, double *, double *, double *,
                                  double *, double *, double *,
                                  double *, double *, double *, double *,
                                  double *, double *, int &, int &);
}

extern "C"
{
        void f_getnp4(int *, double *, double *, double *,
                      double *, double *,
                      double *, double *, double *, double *, double *, double *,
                      double *, double *, double *, double *, double *, double *,
                      double *, double *, double *, double *, double *, double *,
                      double *, double *, double *, double *, double *, double *,
                      double *, double *, double *, double *, double *, double *,
                      double *, double *, double *, double *, double *, double *,
                      double *, double *, int &);
}

extern "C"
{
        void f_getnp4_point(double &, double &, double &,                               // XYZ
                            double &, double &,                                         // chi,trK
                            double &, double &, double &, double &, double &, double &, // gamma_ij
                            double &, double &, double &, double &, double &, double &, // A_ij
                            double &, double &, double &,                               // chi_i
                            double &, double &, double &,                               // trK_i
                            double &, double &, double &,                               // A_ijk
                            double &, double &, double &,
                            double &, double &, double &,
                            double &, double &, double &,
                            double &, double &, double &,
                            double &, double &, double &,
                            double &, double &, double &, double &, double &, double &, // Gam_ijk
                            double &, double &, double &, double &, double &, double &,
                            double &, double &, double &, double &, double &, double &,
                            double &, double &, double &, double &, double &, double &, // R_ij
                            double &, double &);
}

extern "C"
{
        void f_getnp4_ss(int *, double *, double *, double *, double *, double *, double *,
                         double *, double *, double *,
                         double *, double *, double *,
                         double *, double *, double *,
                         double *, double *, double *, double *, double *, double *,
                         double *, double *, double *, double *, double *, double *,
                         double *, double *, double *, double *, double *, double *,
                         double *, double *,
                         double *, double *, double *, double *, double *, double *,
                         double *, double *, double *, double *, double *, double *,
                         double *, double *, double *, double *, double *, double *,
                         double *, double *, double *, double *, double *, double *,
                         double *, double *, double *, double *, double *, double *,
                         double *, double *, double *, double *, double *, double *,
                         double *, double *, int &, int &);
}

extern "C"
{
        void f_getnp4scalar(int *, double *, double *, double *,
                            double *, double *, double *,
                            double *, double *, double *, double *, double *, double *,
                            double *, double *, double *, double *, double *, double *,
                            double *, double *, double *, double *, double *, double *,
                            double *, double *, double *, double *, double *, double *,
                            double *, double *, double *, double *, double *, double *,
                            double *, double *, double *, double *, double *, double *,
                            double *, double *, int &);
}

extern "C"
{
        void f_getnp4scalar_ss(int *, double *, double *, double *, double *, double *, double *,
                               double *, double *, double *,
                               double *, double *, double *,
                               double *, double *, double *,
                               double *, double *, double *, double *, double *, double *,
                               double *, double *, double *, double *, double *, double *,
                               double *, double *, double *, double *, double *, double *,
                               double *, double *, double *,
                               double *, double *, double *, double *, double *, double *,
                               double *, double *, double *, double *, double *, double *,
                               double *, double *, double *, double *, double *, double *,
                               double *, double *, double *, double *, double *, double *,
                               double *, double *, double *, double *, double *, double *,
                               double *, double *, double *, double *, double *, double *,
                               double *, double *, int &, int &);
}

#ifdef USE_GPU
// GPU launch function (backward-compatible: GW=3, tetradtype=2)
void gpu_getnp4_launch(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_chi, const double* d_trK,
    const double* d_dxx, const double* d_gxy, const double* d_gxz,
    const double* d_dyy, const double* d_gyz, const double* d_dzz,
    const double* d_Axx, const double* d_Axy, const double* d_Axz,
    const double* d_Ayy, const double* d_Ayz, const double* d_Azz,
    const double* d_Gamxxx, const double* d_Gamxxy, const double* d_Gamxxz,
    const double* d_Gamxyy, const double* d_Gamxyz, const double* d_Gamxzz,
    const double* d_Gamyxx, const double* d_Gamyxy, const double* d_Gamyxz,
    const double* d_Gamyyy, const double* d_Gamyyz, const double* d_Gamyzz,
    const double* d_Gamzxx, const double* d_Gamzxy, const double* d_Gamzxz,
    const double* d_Gamzyy, const double* d_Gamzyz, const double* d_Gamzzz,
    const double* d_Rxx, const double* d_Rxy, const double* d_Rxz,
    const double* d_Ryy, const double* d_Ryz, const double* d_Rzz,
    double* d_Rpsi4, double* d_Ipsi4,
    int symmetry
);

// GPU launch function for Shell patches (backward-compatible: GW=3, tetradtype=2)
void gpu_getnp4_ss_launch(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_phyx, const double* d_phyy, const double* d_phyz,
    const double* d_drhodx, const double* d_drhody, const double* d_drhodz,
    const double* d_dsigmadx, const double* d_dsigmady, const double* d_dsigmadz,
    const double* d_dRdx, const double* d_dRdy, const double* d_dRdz,
    const double* d_chi, const double* d_trK,
    const double* d_dxx, const double* d_gxy, const double* d_gxz,
    const double* d_dyy, const double* d_gyz, const double* d_dzz,
    const double* d_Axx, const double* d_Axy, const double* d_Axz,
    const double* d_Ayy, const double* d_Ayz, const double* d_Azz,
    const double* d_Gamxxx, const double* d_Gamxxy, const double* d_Gamxxz,
    const double* d_Gamxyy, const double* d_Gamxyz, const double* d_Gamxzz,
    const double* d_Gamyxx, const double* d_Gamyxy, const double* d_Gamyxz,
    const double* d_Gamyyy, const double* d_Gamyyz, const double* d_Gamyzz,
    const double* d_Gamzxx, const double* d_Gamzxy, const double* d_Gamzxz,
    const double* d_Gamzyy, const double* d_Gamzyz, const double* d_Gamzzz,
    const double* d_Rxx, const double* d_Rxy, const double* d_Rxz,
    const double* d_Ryy, const double* d_Ryz, const double* d_Rzz,
    double* d_Rpsi4, double* d_Ipsi4,
    int symmetry, int sst
);

// Template GPU launch function for Shell patches: configurable GW and tetradtype
template<int GW, int TETRAD>
void gpu_getnp4_ss_launch_impl(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_phyx, const double* d_phyy, const double* d_phyz,
    const double* d_drhodx, const double* d_drhody, const double* d_drhodz,
    const double* d_dsigmadx, const double* d_dsigmady, const double* d_dsigmadz,
    const double* d_dRdx, const double* d_dRdy, const double* d_dRdz,
    const double* d_chi, const double* d_trK,
    const double* d_dxx, const double* d_gxy, const double* d_gxz,
    const double* d_dyy, const double* d_gyz, const double* d_dzz,
    const double* d_Axx, const double* d_Axy, const double* d_Axz,
    const double* d_Ayy, const double* d_Ayz, const double* d_Azz,
    const double* d_Gamxxx, const double* d_Gamxxy, const double* d_Gamxxz,
    const double* d_Gamxyy, const double* d_Gamxyz, const double* d_Gamxzz,
    const double* d_Gamyxx, const double* d_Gamyxy, const double* d_Gamyxz,
    const double* d_Gamyyy, const double* d_Gamyyz, const double* d_Gamyzz,
    const double* d_Gamzxx, const double* d_Gamzxy, const double* d_Gamzxz,
    const double* d_Gamzyy, const double* d_Gamzyz, const double* d_Gamzzz,
    const double* d_Rxx, const double* d_Rxy, const double* d_Rxz,
    const double* d_Ryy, const double* d_Ryz, const double* d_Rzz,
    double* d_Rpsi4, double* d_Ipsi4,
    int symmetry, int sst
);

// Template GPU launch function: configurable GW and tetradtype
template<int GW, int TETRAD>
void gpu_getnp4_launch_impl(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_chi, const double* d_trK,
    const double* d_dxx, const double* d_gxy, const double* d_gxz,
    const double* d_dyy, const double* d_gyz, const double* d_dzz,
    const double* d_Axx, const double* d_Axy, const double* d_Axz,
    const double* d_Ayy, const double* d_Ayz, const double* d_Azz,
    const double* d_Gamxxx, const double* d_Gamxxy, const double* d_Gamxxz,
    const double* d_Gamxyy, const double* d_Gamxyz, const double* d_Gamxzz,
    const double* d_Gamyxx, const double* d_Gamyxy, const double* d_Gamyxz,
    const double* d_Gamyyy, const double* d_Gamyyz, const double* d_Gamyzz,
    const double* d_Gamzxx, const double* d_Gamzxy, const double* d_Gamzxz,
    const double* d_Gamzyy, const double* d_Gamzyz, const double* d_Gamzzz,
    const double* d_Rxx, const double* d_Rxy, const double* d_Rxz,
    const double* d_Ryy, const double* d_Ryz, const double* d_Rzz,
    double* d_Rpsi4, double* d_Ipsi4,
    int symmetry
);

// Scalar-field Psi4 GPU launch (wraps getnp4 + scalar rescaling)
template<int GW, int TETRAD>
void gpu_getnp4scalar_launch_impl(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_chi, const double* d_trK, const double* d_Sphi,
    const double* d_dxx, const double* d_gxy, const double* d_gxz,
    const double* d_dyy, const double* d_gyz, const double* d_dzz,
    const double* d_Axx, const double* d_Axy, const double* d_Axz,
    const double* d_Ayy, const double* d_Ayz, const double* d_Azz,
    const double* d_Gamxxx, const double* d_Gamxxy, const double* d_Gamxxz,
    const double* d_Gamxyy, const double* d_Gamxyz, const double* d_Gamxzz,
    const double* d_Gamyxx, const double* d_Gamyxy, const double* d_Gamyxz,
    const double* d_Gamyyy, const double* d_Gamyyz, const double* d_Gamyzz,
    const double* d_Gamzxx, const double* d_Gamzxy, const double* d_Gamzxz,
    const double* d_Gamzyy, const double* d_Gamzyz, const double* d_Gamzzz,
    const double* d_Rxx, const double* d_Rxy, const double* d_Rxz,
    const double* d_Ryy, const double* d_Ryz, const double* d_Rzz,
    double* d_Rpsi4, double* d_Ipsi4,
    int symmetry
);

// Backward-compatible wrapper (GW=ghost_width, TETRAD=tetradtype via macros)
void gpu_getnp4scalar_launch(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_chi, const double* d_trK, const double* d_Sphi,
    const double* d_dxx, const double* d_gxy, const double* d_gxz,
    const double* d_dyy, const double* d_gyz, const double* d_dzz,
    const double* d_Axx, const double* d_Axy, const double* d_Axz,
    const double* d_Ayy, const double* d_Ayz, const double* d_Azz,
    const double* d_Gamxxx, const double* d_Gamxxy, const double* d_Gamxxz,
    const double* d_Gamxyy, const double* d_Gamxyz, const double* d_Gamxzz,
    const double* d_Gamyxx, const double* d_Gamyxy, const double* d_Gamyxz,
    const double* d_Gamyyy, const double* d_Gamyyz, const double* d_Gamyzz,
    const double* d_Gamzxx, const double* d_Gamzxy, const double* d_Gamzxz,
    const double* d_Gamzyy, const double* d_Gamzyz, const double* d_Gamzzz,
    const double* d_Rxx, const double* d_Rxy, const double* d_Rxz,
    const double* d_Ryy, const double* d_Ryz, const double* d_Rzz,
    double* d_Rpsi4, double* d_Ipsi4,
    int symmetry
);

// ABV preprocessing: BSSN -> ADM -> Christoffel + Ricci
template<int GW>
void gpu_bssn2adm_ricci_launch(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_chi, const double* d_trK,
    const double* d_dxx, const double* d_gxy, const double* d_gxz,
    const double* d_dyy, const double* d_gyz, const double* d_dzz,
    const double* d_Axx, const double* d_Axy, const double* d_Axz,
    const double* d_Ayy, const double* d_Ayz, const double* d_Azz,
    double* d_Gamxxx, double* d_Gamxxy, double* d_Gamxxz,
    double* d_Gamxyy, double* d_Gamxyz, double* d_Gamxzz,
    double* d_Gamyxx, double* d_Gamyxy, double* d_Gamyxz,
    double* d_Gamyyy, double* d_Gamyyz, double* d_Gamyzz,
    double* d_Gamzxx, double* d_Gamzxy, double* d_Gamzxz,
    double* d_Gamzyy, double* d_Gamzyz, double* d_Gamzzz,
    double* d_Rxx, double* d_Rxy, double* d_Rxz,
    double* d_Ryy, double* d_Ryz, double* d_Rzz,
    double* d_temp_adm_dxx, double* d_temp_adm_gxy, double* d_temp_adm_gxz,
    double* d_temp_adm_dyy, double* d_temp_adm_gyz, double* d_temp_adm_dzz,
    int symmetry
);

// Non-template wrapper for calling from C++ code
void gpu_bssn2adm_ricci_launch_wrapper(
    cudaStream_t stream,
    int ex[3], int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_chi, const double* d_trK,
    const double* d_dxx, const double* d_gxy, const double* d_gxz,
    const double* d_dyy, const double* d_gyz, const double* d_dzz,
    const double* d_Axx, const double* d_Axy, const double* d_Axz,
    const double* d_Ayy, const double* d_Ayz, const double* d_Azz,
    double* d_Gamxxx, double* d_Gamxxy, double* d_Gamxxz,
    double* d_Gamxyy, double* d_Gamxyz, double* d_Gamxzz,
    double* d_Gamyxx, double* d_Gamyxy, double* d_Gamyxz,
    double* d_Gamyyy, double* d_Gamyyz, double* d_Gamyzz,
    double* d_Gamzxx, double* d_Gamzxy, double* d_Gamzxz,
    double* d_Gamzyy, double* d_Gamzyz, double* d_Gamzzz,
    double* d_Rxx, double* d_Rxy, double* d_Rxz,
    double* d_Ryy, double* d_Ryz, double* d_Rzz,
    double* d_temp_adm_dxx, double* d_temp_adm_gxy, double* d_temp_adm_gxz,
    double* d_temp_adm_dyy, double* d_temp_adm_gyz, double* d_temp_adm_dzz,
    int symmetry
);
#endif /* USE_GPU */

#endif /* GETNP4_H */
