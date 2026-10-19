
#ifndef BSSN_H
#define BSSN_H

#ifdef fortran1
#define f_compute_rhs_bssn compute_rhs_bssn
#define f_compute_rhs_bssn_ss compute_rhs_bssn_ss
#define f_compute_rhs_bssn_escalar compute_rhs_bssn_escalar
#define f_compute_rhs_bssn_escalar_ss compute_rhs_bssn_escalar_ss
#define f_compute_rhs_Z4c compute_rhs_z4c
#define f_compute_rhs_Z4cnot compute_rhs_z4cnot
#define f_compute_rhs_Z4c_ss compute_rhs_z4c_ss
#define f_compute_constraint_fr compute_constraint_fr
#endif
#ifdef fortran2
#define f_compute_rhs_bssn COMPUTE_RHS_BSSN
#define f_compute_rhs_bssn_ss COMPUTE_RHS_BSSN_SS
#define f_compute_rhs_bssn_escalar COMPUTE_RHS_BSSN_ESCALAR
#define f_compute_rhs_bssn_escalar_ss COMPUTE_RHS_BSSN_ESCALAR_SS
#define f_compute_rhs_Z4c COMPUTE_RHS_Z4C
#define f_compute_rhs_Z4cnot COMPUTE_RHS_Z4CNOT
#define f_compute_rhs_Z4c_ss COMPUTE_RHS_Z4C_SS
#define f_compute_constraint_fr COMPUTE_CONSTRAINT_FR
#endif
#ifdef fortran3
#define f_compute_rhs_bssn compute_rhs_bssn_
#define f_compute_rhs_bssn_ss compute_rhs_bssn_ss_
#define f_compute_rhs_bssn_escalar compute_rhs_bssn_escalar_
#define f_compute_rhs_bssn_escalar_ss compute_rhs_bssn_escalar_ss_
#define f_compute_rhs_Z4c compute_rhs_z4c_
#define f_compute_rhs_Z4cnot compute_rhs_z4cnot_
#define f_compute_rhs_Z4c_ss compute_rhs_z4c_ss_
#define f_compute_constraint_fr compute_constraint_fr_
#endif
extern "C"
{
        int f_compute_rhs_bssn(int *, double &, double *, double *, double *,                                                      // ex,T,X,Y,Z
                               double *, double *,                                                                                 // chi, trK
                               double *, double *, double *, double *, double *, double *,                                         // gij
                               double *, double *, double *, double *, double *, double *,                                         // Aij
                               double *, double *, double *,                                                                       // Gam
                               double *, double *, double *, double *, double *, double *, double *,                               // Gauge
                               double *, double *,                                                                                 // chi, trK
                               double *, double *, double *, double *, double *, double *,                                         // gij
                               double *, double *, double *, double *, double *, double *,                                         // Aij
                               double *, double *, double *,                                                                       // Gam
                               double *, double *, double *, double *, double *, double *, double *,                               // Gauge
                               double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, // stress-energy
                               double *, double *, double *, double *, double *, double *,                                         // Christoffel
                               double *, double *, double *, double *, double *, double *,                                         // Christoffel
                               double *, double *, double *, double *, double *, double *,                                         // Christoffel
                               double *, double *, double *, double *, double *, double *,                                         // Ricci
                               double *, double *, double *, double *, double *, double *, double *,                               // constraint violation
                               int &, int &, double &, int &);
}

extern "C"
{
        int f_compute_rhs_bssn_ss(int *, double &, double *, double *, double *,                                                      // ex,T,rho,sigma,R
                                  double *, double *, double *,                                                                       // X,Y,Z
                                  double *, double *, double *,                                                                       // drhodx,drhody,drhodz
                                  double *, double *, double *,                                                                       // dsigmadx,dsigmady,dsigmadz
                                  double *, double *, double *,                                                                       // dRdx,dRdy,dRdz
                                  double *, double *, double *, double *, double *, double *,                                         // drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz
                                  double *, double *, double *, double *, double *, double *,                                         // dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz
                                  double *, double *, double *, double *, double *, double *,                                         // dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz
                                  double *, double *,                                                                                 // chi, trK
                                  double *, double *, double *, double *, double *, double *,                                         // gij
                                  double *, double *, double *, double *, double *, double *,                                         // Aij
                                  double *, double *, double *,                                                                       // Gam
                                  double *, double *, double *, double *, double *, double *, double *,                               // Gauge
                                  double *, double *,                                                                                 // chi, trK
                                  double *, double *, double *, double *, double *, double *,                                         // gij
                                  double *, double *, double *, double *, double *, double *,                                         // Aij
                                  double *, double *, double *,                                                                       // Gam
                                  double *, double *, double *, double *, double *, double *, double *,                               // Gauge
                                  double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, // stress-energy
                                  double *, double *, double *, double *, double *, double *,                                         // Christoffel
                                  double *, double *, double *, double *, double *, double *,                                         // Christoffel
                                  double *, double *, double *, double *, double *, double *,                                         // Christoffel
                                  double *, double *, double *, double *, double *, double *,                                         // Ricci
                                  double *, double *, double *, double *, double *, double *, double *,                               // constraint violation
                                  int &, int &, double &, int &, int &);
}

extern "C"
{
        int f_compute_rhs_bssn_escalar(int *, double &, double *, double *, double *,                                                      // ex,T,X,Y,Z
                                       double *, double *,                                                                                 // chi, trK
                                       double *, double *, double *, double *, double *, double *,                                         // gij
                                       double *, double *, double *, double *, double *, double *,                                         // Aij
                                       double *, double *, double *,                                                                       // Gam
                                       double *, double *, double *, double *, double *, double *, double *,                               // Gauge
                                       double *, double *,                                                                                 // Sphi, Spi
                                       double *, double *,                                                                                 // chi, trK
                                       double *, double *, double *, double *, double *, double *,                                         // gij
                                       double *, double *, double *, double *, double *, double *,                                         // Aij
                                       double *, double *, double *,                                                                       // Gam
                                       double *, double *, double *, double *, double *, double *, double *,                               // Gauge
                                       double *, double *,                                                                                 // Sphi, Spi
                                       double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, // stress-energy
                                       double *, double *, double *, double *, double *, double *,                                         // Christoffel
                                       double *, double *, double *, double *, double *, double *,                                         // Christoffel
                                       double *, double *, double *, double *, double *, double *,                                         // Christoffel
                                       double *, double *, double *, double *, double *, double *,                                         // Ricci
                                       double *, double *, double *, double *, double *, double *, double *,                               // constraint violation
                                       int &, int &, double &, int &);
}

extern "C"
{
        int f_compute_rhs_bssn_escalar_ss(int *, double &, double *, double *, double *,                                                      // ex,T,rho,sigma,R
                                          double *, double *, double *,                                                                       // X,Y,Z
                                          double *, double *, double *,                                                                       // drhodx,drhody,drhodz
                                          double *, double *, double *,                                                                       // dsigmadx,dsigmady,dsigmadz
                                          double *, double *, double *,                                                                       // dRdx,dRdy,dRdz
                                          double *, double *, double *, double *, double *, double *,                                         // drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz
                                          double *, double *, double *, double *, double *, double *,                                         // dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz
                                          double *, double *, double *, double *, double *, double *,                                         // dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz
                                          double *, double *,                                                                                 // chi, trK
                                          double *, double *, double *, double *, double *, double *,                                         // gij
                                          double *, double *, double *, double *, double *, double *,                                         // Aij
                                          double *, double *, double *,                                                                       // Gam
                                          double *, double *, double *, double *, double *, double *, double *,                               // Gauge
                                          double *, double *,                                                                                 // Sphi,Spi
                                          double *, double *,                                                                                 // chi, trK
                                          double *, double *, double *, double *, double *, double *,                                         // gij
                                          double *, double *, double *, double *, double *, double *,                                         // Aij
                                          double *, double *, double *,                                                                       // Gam
                                          double *, double *, double *, double *, double *, double *, double *,                               // Gauge
                                          double *, double *,                                                                                 // Sphi,Spi
                                          double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, // stress-energy
                                          double *, double *, double *, double *, double *, double *,                                         // Christoffel
                                          double *, double *, double *, double *, double *, double *,                                         // Christoffel
                                          double *, double *, double *, double *, double *, double *,                                         // Christoffel
                                          double *, double *, double *, double *, double *, double *,                                         // Ricci
                                          double *, double *, double *, double *, double *, double *, double *,                               // constraint violation
                                          int &, int &, double &, int &, int &);
}

extern "C"
{
        int f_compute_rhs_Z4c(int *, double &, double *, double *, double *,                        // ex,T,X,Y,Z
                              double *, double *,                                                   // chi, trK
                              double *, double *, double *, double *, double *, double *,           // gij
                              double *, double *, double *, double *, double *, double *,           // Aij
                              double *, double *, double *,                                         // Gam
                              double *, double *, double *, double *, double *, double *, double *, // Gauge
                              double *,                                                             // Z4
                              double *, double *,                                                   // chi, trK
                              double *, double *, double *, double *, double *, double *,           // gij
                              double *, double *, double *, double *, double *, double *,           // Aij
                              double *, double *, double *,                                         // Gam
                              double *, double *, double *, double *, double *, double *, double *, // Gauge
                              double *,                                                             // Z4
                              double *, double *, double *, double *, double *, double *, double *, double *, double *, double *,
                              double *, double *, double *, double *, double *, double *,
                              double *, double *, double *, double *, double *, double *,
                              double *, double *, double *, double *, double *, double *,
                              double *, double *, double *, double *, double *, double *,
                              double *, double *, double *, double *, double *, double *, double *,
                              int &, int &, double &, int &);
}

extern "C"
{
        int f_compute_rhs_Z4c_ss(int *, double &, double *, double *, double *,                                                      // ex,T,rho,sigma,R
                                 double *, double *, double *,                                                                       // X,Y,Z
                                 double *, double *, double *,                                                                       // drhodx,drhody,drhodz
                                 double *, double *, double *,                                                                       // dsigmadx,dsigmady,dsigmadz
                                 double *, double *, double *,                                                                       // dRdx,dRdy,dRdz
                                 double *, double *, double *, double *, double *, double *,                                         // drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz
                                 double *, double *, double *, double *, double *, double *,                                         // dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz
                                 double *, double *, double *, double *, double *, double *,                                         // dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz
                                 double *, double *,                                                                                 // chi, trK
                                 double *, double *, double *, double *, double *, double *,                                         // gij
                                 double *, double *, double *, double *, double *, double *,                                         // Aij
                                 double *, double *, double *,                                                                       // Gam
                                 double *, double *, double *, double *, double *, double *, double *,                               // Gauge
                                 double *,                                                                                           // TZ
                                 double *, double *,                                                                                 // chi, trK
                                 double *, double *, double *, double *, double *, double *,                                         // gij
                                 double *, double *, double *, double *, double *, double *,                                         // Aij
                                 double *, double *, double *,                                                                       // Gam
                                 double *, double *, double *, double *, double *, double *, double *,                               // Gauge
                                 double *,                                                                                           // TZ
                                 double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, // stress-energy
                                 double *, double *, double *, double *, double *, double *,                                         // Christoffel
                                 double *, double *, double *, double *, double *, double *,                                         // Christoffel
                                 double *, double *, double *, double *, double *, double *,                                         // Christoffel
                                 double *, double *, double *, double *, double *, double *,                                         // Ricci
                                 double *, double *, double *, double *, double *, double *, double *,                               // constraint violation
                                 int &, int &, double &, int &, int &);
}

extern "C"
{
        int f_compute_rhs_Z4cnot(int *, double &, double *, double *, double *,                        // ex,T,X,Y,Z
                                 double *, double *,                                                   // chi, trK
                                 double *, double *, double *, double *, double *, double *,           // gij
                                 double *, double *, double *, double *, double *, double *,           // Aij
                                 double *, double *, double *,                                         // Gam
                                 double *, double *, double *, double *, double *, double *, double *, // Gauge
                                 double *,                                                             // Z4
                                 double *, double *,                                                   // chi, trK
                                 double *, double *, double *, double *, double *, double *,           // gij
                                 double *, double *, double *, double *, double *, double *,           // Aij
                                 double *, double *, double *,                                         // Gam
                                 double *, double *, double *, double *, double *, double *, double *, // Gauge
                                 double *,                                                             // Z4
                                 double *, double *, double *, double *, double *, double *, double *, double *, double *, double *,
                                 double *, double *, double *, double *, double *, double *,
                                 double *, double *, double *, double *, double *, double *,
                                 double *, double *, double *, double *, double *, double *,
                                 double *, double *, double *, double *, double *, double *,
                                 double *, double *, double *, double *, double *, double *, double *,
                                 int &, int &, double &, int &, double &);
}

extern "C"
{
        void f_compute_constraint_fr(int *, double *, double *, double *,                        // ex,X,Y,Z
                                     double *, double *, double *, double *,                     // chi, trK,rho,Sphi
                                     double *, double *, double *, double *, double *, double *, // gij
                                     double *, double *, double *, double *, double *, double *, // Aij
                                     double *, double *, double *, double *, double *, double *, // Rij
                                     double *, double *, double *, double *, double *, double *, // Sij
                                     double *);
} // FR_cons

void gpu_compute_rhs_bssn_launch( // launch kernel with device pointers
    cudaStream_t &stream,
    int* ex, double T, double* d_X, double* d_Y, double* d_Z,
    double* d_chi, double* d_trK,
    double* d_dxx, double* d_gxy, double* d_gxz,
    double* d_dyy, double* d_gyz, double* d_dzz,
    double* d_Axx, double* d_Axy, double* d_Axz,
    double* d_Ayy, double* d_Ayz, double* d_Azz,
    double* d_Gamx, double* d_Gamy, double* d_Gamz,
    double* d_Lap,
    double* d_betax, double* d_betay, double* d_betaz,
    double* d_dtSfx, double* d_dtSfy, double* d_dtSfz,
    double* d_chi_rhs, double* d_trK_rhs,
    double* d_gxx_rhs, double* d_gxy_rhs, double* d_gxz_rhs,
    double* d_gyy_rhs, double* d_gyz_rhs, double* d_gzz_rhs,
    double* d_Axx_rhs, double* d_Axy_rhs, double* d_Axz_rhs,
    double* d_Ayy_rhs, double* d_Ayz_rhs, double* d_Azz_rhs,
    double* d_Gamx_rhs, double* d_Gamy_rhs, double* d_Gamz_rhs,
    double* d_Lap_rhs,
    double* d_betax_rhs, double* d_betay_rhs, double* d_betaz_rhs,
    double* d_dtSfx_rhs, double* d_dtSfy_rhs, double* d_dtSfz_rhs,
    double* d_rho, double* d_Sx, double* d_Sy, double* d_Sz,
    double* d_Sxx, double* d_Sxy, double* d_Sxz,
    double* d_Syy, double* d_Syz, double* d_Szz,
    double* d_Gamxxx, double* d_Gamxxy, double* d_Gamxxz,
    double* d_Gamxyy, double* d_Gamxyz, double* d_Gamxzz,
    double* d_Gamyxx, double* d_Gamyxy, double* d_Gamyxz,
    double* d_Gamyyy, double* d_Gamyyz, double* d_Gamyzz,
    double* d_Gamzxx, double* d_Gamzxy, double* d_Gamzxz,
    double* d_Gamzyy, double* d_Gamzyz, double* d_Gamzzz,
    double* d_Rxx, double* d_Rxy, double* d_Rxz,
    double* d_Ryy, double* d_Ryz, double* d_Rzz,
    double* d_ham_Res, double* d_movx_Res, double* d_movy_Res, double* d_movz_Res,
    double* d_Gmx_Res, double* d_Gmy_Res, double* d_Gmz_Res,
    int symmetry, int lev, double eps, int co
);

void gpu_compute_rhs_bssn_escalar_launch( // launch BSSN-EScalar RHS kernels with device pointers
    cudaStream_t &stream,
    int* ex, double T, double* d_X, double* d_Y, double* d_Z,
    double* d_chi, double* d_trK,
    double* d_dxx, double* d_gxy, double* d_gxz,
    double* d_dyy, double* d_gyz, double* d_dzz,
    double* d_Axx, double* d_Axy, double* d_Axz,
    double* d_Ayy, double* d_Ayz, double* d_Azz,
    double* d_Gamx, double* d_Gamy, double* d_Gamz,
    double* d_Lap,
    double* d_betax, double* d_betay, double* d_betaz,
    double* d_dtSfx, double* d_dtSfy, double* d_dtSfz,
    double* d_Sphi, double* d_Spi,
    double* d_chi_rhs, double* d_trK_rhs,
    double* d_gxx_rhs, double* d_gxy_rhs, double* d_gxz_rhs,
    double* d_gyy_rhs, double* d_gyz_rhs, double* d_gzz_rhs,
    double* d_Axx_rhs, double* d_Axy_rhs, double* d_Axz_rhs,
    double* d_Ayy_rhs, double* d_Ayz_rhs, double* d_Azz_rhs,
    double* d_Gamx_rhs, double* d_Gamy_rhs, double* d_Gamz_rhs,
    double* d_Lap_rhs,
    double* d_betax_rhs, double* d_betay_rhs, double* d_betaz_rhs,
    double* d_dtSfx_rhs, double* d_dtSfy_rhs, double* d_dtSfz_rhs,
    double* d_Sphi_rhs, double* d_Spi_rhs,
    double* d_rho, double* d_Sx, double* d_Sy, double* d_Sz,
    double* d_Sxx, double* d_Sxy, double* d_Sxz,
    double* d_Syy, double* d_Syz, double* d_Szz,
    double* d_Gamxxx, double* d_Gamxxy, double* d_Gamxxz,
    double* d_Gamxyy, double* d_Gamxyz, double* d_Gamxzz,
    double* d_Gamyxx, double* d_Gamyxy, double* d_Gamyxz,
    double* d_Gamyyy, double* d_Gamyyz, double* d_Gamyzz,
    double* d_Gamzxx, double* d_Gamzxy, double* d_Gamzxz,
    double* d_Gamzyy, double* d_Gamzyz, double* d_Gamzzz,
    double* d_Rxx, double* d_Rxy, double* d_Rxz,
    double* d_Ryy, double* d_Ryz, double* d_Rzz,
    double* d_ham_Res, double* d_movx_Res, double* d_movy_Res, double* d_movz_Res,
    double* d_Gmx_Res, double* d_Gmy_Res, double* d_Gmz_Res,
    int symmetry, int lev, double eps, int co, double a2
);

void gpu_compute_constraint_fr_launch(
    cudaStream_t &stream,
    int* ex,
    double* d_chi,
    double* d_trK,
    double* d_rho,
    double* d_Sphi,
    double* d_dxx,
    double* d_gxy,
    double* d_gxz,
    double* d_dyy,
    double* d_gyz,
    double* d_dzz,
    double* d_Axx,
    double* d_Axy,
    double* d_Axz,
    double* d_Ayy,
    double* d_Ayz,
    double* d_Azz,
    double* d_Rxx,
    double* d_Rxy,
    double* d_Rxz,
    double* d_Ryy,
    double* d_Ryz,
    double* d_Rzz,
    double* d_Sxx,
    double* d_Sxy,
    double* d_Sxz,
    double* d_Syy,
    double* d_Syz,
    double* d_Szz,
    double* d_Cons_fR,
    double a2
);

void gpu_compute_rhs_bssn_em_launch( // launch BSSN-EM RHS kernels with device pointers
    cudaStream_t &stream,
    int* ex, double T, double* d_X, double* d_Y, double* d_Z,
    double* d_chi, double* d_trK,
    double* d_dxx, double* d_gxy, double* d_gxz,
    double* d_dyy, double* d_gyz, double* d_dzz,
    double* d_Axx, double* d_Axy, double* d_Axz,
    double* d_Ayy, double* d_Ayz, double* d_Azz,
    double* d_Gamx, double* d_Gamy, double* d_Gamz,
    double* d_Lap,
    double* d_betax, double* d_betay, double* d_betaz,
    double* d_dtSfx, double* d_dtSfy, double* d_dtSfz,
    double* d_Ex, double* d_Ey, double* d_Ez,
    double* d_Bx, double* d_By, double* d_Bz,
    double* d_Kpsi, double* d_Kphi,
    double* d_Jx, double* d_Jy, double* d_Jz, double* d_qchar,
    double* d_chi_rhs, double* d_trK_rhs,
    double* d_gxx_rhs, double* d_gxy_rhs, double* d_gxz_rhs,
    double* d_gyy_rhs, double* d_gyz_rhs, double* d_gzz_rhs,
    double* d_Axx_rhs, double* d_Axy_rhs, double* d_Axz_rhs,
    double* d_Ayy_rhs, double* d_Ayz_rhs, double* d_Azz_rhs,
    double* d_Gamx_rhs, double* d_Gamy_rhs, double* d_Gamz_rhs,
    double* d_Lap_rhs,
    double* d_betax_rhs, double* d_betay_rhs, double* d_betaz_rhs,
    double* d_dtSfx_rhs, double* d_dtSfy_rhs, double* d_dtSfz_rhs,
    double* d_Ex_rhs, double* d_Ey_rhs, double* d_Ez_rhs,
    double* d_Bx_rhs, double* d_By_rhs, double* d_Bz_rhs,
    double* d_Kpsi_rhs, double* d_Kphi_rhs,
    double* d_rho, double* d_Sx, double* d_Sy, double* d_Sz,
    double* d_Sxx, double* d_Sxy, double* d_Sxz,
    double* d_Syy, double* d_Syz, double* d_Szz,
    double* d_Gamxxx, double* d_Gamxxy, double* d_Gamxxz,
    double* d_Gamxyy, double* d_Gamxyz, double* d_Gamxzz,
    double* d_Gamyxx, double* d_Gamyxy, double* d_Gamyxz,
    double* d_Gamyyy, double* d_Gamyyz, double* d_Gamyzz,
    double* d_Gamzxx, double* d_Gamzxy, double* d_Gamzxz,
    double* d_Gamzyy, double* d_Gamzyz, double* d_Gamzzz,
    double* d_Rxx, double* d_Rxy, double* d_Rxz,
    double* d_Ryy, double* d_Ryz, double* d_Rzz,
    double* d_ham_Res, double* d_movx_Res, double* d_movy_Res, double* d_movz_Res,
    double* d_Gmx_Res, double* d_Gmy_Res, double* d_Gmz_Res,
    int symmetry, int lev, double eps, int co
);

void gpu_compute_rhs_z4c_launch( // launch Z4C RHS kernels with device pointers
    cudaStream_t &stream,
    int* ex, double T, double* d_X, double* d_Y, double* d_Z,
    double* d_chi, double* d_trK,
    double* d_dxx, double* d_gxy, double* d_gxz,
    double* d_dyy, double* d_gyz, double* d_dzz,
    double* d_Axx, double* d_Axy, double* d_Axz,
    double* d_Ayy, double* d_Ayz, double* d_Azz,
    double* d_Gamx, double* d_Gamy, double* d_Gamz,
    double* d_Lap,
    double* d_betax, double* d_betay, double* d_betaz,
    double* d_dtSfx, double* d_dtSfy, double* d_dtSfz,
    double* d_TZ,
    double* d_chi_rhs, double* d_trK_rhs,
    double* d_gxx_rhs, double* d_gxy_rhs, double* d_gxz_rhs,
    double* d_gyy_rhs, double* d_gyz_rhs, double* d_gzz_rhs,
    double* d_Axx_rhs, double* d_Axy_rhs, double* d_Axz_rhs,
    double* d_Ayy_rhs, double* d_Ayz_rhs, double* d_Azz_rhs,
    double* d_Gamx_rhs, double* d_Gamy_rhs, double* d_Gamz_rhs,
    double* d_Lap_rhs,
    double* d_betax_rhs, double* d_betay_rhs, double* d_betaz_rhs,
    double* d_dtSfx_rhs, double* d_dtSfy_rhs, double* d_dtSfz_rhs,
    double* d_TZ_rhs,
    double* d_rho, double* d_Sx, double* d_Sy, double* d_Sz,
    double* d_Sxx, double* d_Sxy, double* d_Sxz,
    double* d_Syy, double* d_Syz, double* d_Szz,
    double* d_Gamxxx, double* d_Gamxxy, double* d_Gamxxz,
    double* d_Gamxyy, double* d_Gamxyz, double* d_Gamxzz,
    double* d_Gamyxx, double* d_Gamyxy, double* d_Gamyxz,
    double* d_Gamyyy, double* d_Gamyyz, double* d_Gamyzz,
    double* d_Gamzxx, double* d_Gamzxy, double* d_Gamzxz,
    double* d_Gamzyy, double* d_Gamzyz, double* d_Gamzzz,
    double* d_Rxx, double* d_Rxy, double* d_Rxz,
    double* d_Ryy, double* d_Ryz, double* d_Rzz,
    double* d_ham_Res, double* d_movx_Res, double* d_movy_Res, double* d_movz_Res,
    double* d_Gmx_Res, double* d_Gmy_Res, double* d_Gmz_Res,
    int symmetry, int lev, double eps, int co
);

#endif /* BSSN_H */
