
#ifndef BSSN_H
#define BSSN_H

#ifdef USE_GPU
#include <cuda_runtime_api.h>
#endif

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

#ifdef USE_GPU
void gpu_compute_rhs_bssn_escalar_launch(
    cudaStream_t &stream, int* ex, int* pex, double T, double* X, double* Y, double* Z,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Gamx, double* Gamy, double* Gamz, double* Lap,
    double* betax, double* betay, double* betaz,
    double* dtSfx, double* dtSfy, double* dtSfz,
    double* Sphi, double* Spi,
    double* chi_rhs, double* trK_rhs,
    double* gxx_rhs, double* gxy_rhs, double* gxz_rhs, double* gyy_rhs, double* gyz_rhs, double* gzz_rhs,
    double* Axx_rhs, double* Axy_rhs, double* Axz_rhs, double* Ayy_rhs, double* Ayz_rhs, double* Azz_rhs,
    double* Gamx_rhs, double* Gamy_rhs, double* Gamz_rhs, double* Lap_rhs,
    double* betax_rhs, double* betay_rhs, double* betaz_rhs,
    double* dtSfx_rhs, double* dtSfy_rhs, double* dtSfz_rhs,
    double* Sphi_rhs, double* Spi_rhs,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Sxx, double* Sxy, double* Sxz, double* Syy, double* Syz, double* Szz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz, double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz, double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz, double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* ham_Res, double* movx_Res, double* movy_Res, double* movz_Res,
    double* Gmx_Res, double* Gmy_Res, double* Gmz_Res,
    int symmetry, int lev, double eps, int co);

void gpu_compute_rhs_bssn_escalar_ss_launch(
    cudaStream_t &stream, int* ex, int* pex, double T,
    double* crho, double* sigma, double* R,
    double* X, double* Y, double* Z,
    double* drhodx, double* drhody, double* drhodz,
    double* dsigmadx, double* dsigmady, double* dsigmadz,
    double* dRdx, double* dRdy, double* dRdz,
    double* drhodxx, double* drhodxy, double* drhodxz,
    double* drhodyy, double* drhodyz, double* drhodzz,
    double* dsigmadxx, double* dsigmadxy, double* dsigmadxz,
    double* dsigmadyy, double* dsigmadyz, double* dsigmadzz,
    double* dRdxx, double* dRdxy, double* dRdxz,
    double* dRdyy, double* dRdyz, double* dRdzz,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Gamx, double* Gamy, double* Gamz, double* Lap,
    double* betax, double* betay, double* betaz,
    double* dtSfx, double* dtSfy, double* dtSfz,
    double* Sphi, double* Spi,
    double* chi_rhs, double* trK_rhs,
    double* gxx_rhs, double* gxy_rhs, double* gxz_rhs, double* gyy_rhs, double* gyz_rhs, double* gzz_rhs,
    double* Axx_rhs, double* Axy_rhs, double* Axz_rhs, double* Ayy_rhs, double* Ayz_rhs, double* Azz_rhs,
    double* Gamx_rhs, double* Gamy_rhs, double* Gamz_rhs, double* Lap_rhs,
    double* betax_rhs, double* betay_rhs, double* betaz_rhs,
    double* dtSfx_rhs, double* dtSfy_rhs, double* dtSfz_rhs,
    double* Sphi_rhs, double* Spi_rhs,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Sxx, double* Sxy, double* Sxz, double* Syy, double* Syz, double* Szz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz, double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz, double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz, double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* ham_Res, double* movx_Res, double* movy_Res, double* movz_Res,
    double* Gmx_Res, double* Gmy_Res, double* Gmz_Res,
    int symmetry, int lev, double eps, int sst, int co);
#endif /* USE_GPU */

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

#ifdef USE_GPU
void gpu_compute_rhs_bssn_launch(
    cudaStream_t &stream, int* ex, int* pex, double T, double* X, double* Y, double* Z,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Gamx, double* Gamy, double* Gamz, double* Lap,
    double* betax, double* betay, double* betaz,
    double* dtSfx, double* dtSfy, double* dtSfz,
    double* chi_rhs, double* trK_rhs,
    double* gxx_rhs, double* gxy_rhs, double* gxz_rhs, double* gyy_rhs, double* gyz_rhs, double* gzz_rhs,
    double* Axx_rhs, double* Axy_rhs, double* Axz_rhs, double* Ayy_rhs, double* Ayz_rhs, double* Azz_rhs,
    double* Gamx_rhs, double* Gamy_rhs, double* Gamz_rhs, double* Lap_rhs,
    double* betax_rhs, double* betay_rhs, double* betaz_rhs,
    double* dtSfx_rhs, double* dtSfy_rhs, double* dtSfz_rhs,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Sxx, double* Sxy, double* Sxz, double* Syy, double* Syz, double* Szz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz, double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz, double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz, double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* ham_Res, double* movx_Res, double* movy_Res, double* movz_Res,
    double* Gmx_Res, double* Gmy_Res, double* Gmz_Res,
    int symmetry, int lev, double eps, int co
);

void gpu_compute_rhs_bssn_ss_launch(
    cudaStream_t &stream, int* ex, int* pex, double T,
    double* crho, double* sigma, double* R,
    double* X, double* Y, double* Z,
    double* drhodx, double* drhody, double* drhodz,
    double* dsigmadx, double* dsigmady, double* dsigmadz,
    double* dRdx, double* dRdy, double* dRdz,
    double* drhodxx, double* drhodxy, double* drhodxz,
    double* drhodyy, double* drhodyz, double* drhodzz,
    double* dsigmadxx, double* dsigmadxy, double* dsigmadxz,
    double* dsigmadyy, double* dsigmadyz, double* dsigmadzz,
    double* dRdxx, double* dRdxy, double* dRdxz,
    double* dRdyy, double* dRdyz, double* dRdzz,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Gamx, double* Gamy, double* Gamz, double* Lap,
    double* betax, double* betay, double* betaz,
    double* dtSfx, double* dtSfy, double* dtSfz,
    double* chi_rhs, double* trK_rhs,
    double* gxx_rhs, double* gxy_rhs, double* gxz_rhs, double* gyy_rhs, double* gyz_rhs, double* gzz_rhs,
    double* Axx_rhs, double* Axy_rhs, double* Axz_rhs, double* Ayy_rhs, double* Ayz_rhs, double* Azz_rhs,
    double* Gamx_rhs, double* Gamy_rhs, double* Gamz_rhs, double* Lap_rhs,
    double* betax_rhs, double* betay_rhs, double* betaz_rhs,
    double* dtSfx_rhs, double* dtSfy_rhs, double* dtSfz_rhs,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Sxx, double* Sxy, double* Sxz, double* Syy, double* Syz, double* Szz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz, double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz, double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz, double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* ham_Res, double* movx_Res, double* movy_Res, double* movz_Res,
    double* Gmx_Res, double* Gmy_Res, double* Gmz_Res,
    int symmetry, int lev, double eps, int co, int sst
);

void gpu_ricci_only_launch(
    cudaStream_t &stream, int* ex, int* pex, double T, double* X, double* Y, double* Z,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Gamx, double* Gamy, double* Gamz, double* Lap,
    double* betax, double* betay, double* betaz,
    double* dtSfx, double* dtSfy, double* dtSfz,
    double* chi_rhs, double* trK_rhs,
    double* gxx_rhs, double* gxy_rhs, double* gxz_rhs, double* gyy_rhs, double* gyz_rhs, double* gzz_rhs,
    double* Axx_rhs, double* Axy_rhs, double* Axz_rhs, double* Ayy_rhs, double* Ayz_rhs, double* Azz_rhs,
    double* Gamx_rhs, double* Gamy_rhs, double* Gamz_rhs, double* Lap_rhs,
    double* betax_rhs, double* betay_rhs, double* betaz_rhs,
    double* dtSfx_rhs, double* dtSfy_rhs, double* dtSfz_rhs,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Sxx, double* Sxy, double* Sxz, double* Syy, double* Syz, double* Szz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz, double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz, double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz, double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* ham_Res, double* movx_Res, double* movy_Res, double* movz_Res,
    double* Gmx_Res, double* Gmy_Res, double* Gmz_Res,
    int symmetry, int lev, double eps, int co
);

void gpu_rhs_rest_launch(
    cudaStream_t &stream, int* ex, int* pex, double T, double* X, double* Y, double* Z,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Gamx, double* Gamy, double* Gamz, double* Lap,
    double* betax, double* betay, double* betaz,
    double* dtSfx, double* dtSfy, double* dtSfz,
    double* chi_rhs, double* trK_rhs,
    double* gxx_rhs, double* gxy_rhs, double* gxz_rhs, double* gyy_rhs, double* gyz_rhs, double* gzz_rhs,
    double* Axx_rhs, double* Axy_rhs, double* Axz_rhs, double* Ayy_rhs, double* Ayz_rhs, double* Azz_rhs,
    double* Gamx_rhs, double* Gamy_rhs, double* Gamz_rhs, double* Lap_rhs,
    double* betax_rhs, double* betay_rhs, double* betaz_rhs,
    double* dtSfx_rhs, double* dtSfy_rhs, double* dtSfz_rhs,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Sxx, double* Sxy, double* Sxz, double* Syy, double* Syz, double* Szz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz, double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz, double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz, double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* ham_Res, double* movx_Res, double* movy_Res, double* movz_Res,
    double* Gmx_Res, double* Gmy_Res, double* Gmz_Res,
    int symmetry, int lev, double eps, int co
);

void gpu_compute_constraints_only_launch(
    cudaStream_t &stream, int* ex, int* pex, double T, double* X, double* Y, double* Z,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Gamx, double* Gamy, double* Gamz, double* Lap,
    double* betax, double* betay, double* betaz,
    double* dtSfx, double* dtSfy, double* dtSfz,
    double* chi_rhs, double* trK_rhs,
    double* gxx_rhs, double* gxy_rhs, double* gxz_rhs, double* gyy_rhs, double* gyz_rhs, double* gzz_rhs,
    double* Axx_rhs, double* Axy_rhs, double* Axz_rhs, double* Ayy_rhs, double* Ayz_rhs, double* Azz_rhs,
    double* Gamx_rhs, double* Gamy_rhs, double* Gamz_rhs, double* Lap_rhs,
    double* betax_rhs, double* betay_rhs, double* betaz_rhs,
    double* dtSfx_rhs, double* dtSfy_rhs, double* dtSfz_rhs,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Sxx, double* Sxy, double* Sxz, double* Syy, double* Syz, double* Szz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz, double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz, double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz, double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* ham_Res, double* movx_Res, double* movy_Res, double* movz_Res,
    double* Gmx_Res, double* Gmy_Res, double* Gmz_Res,
    int symmetry, int lev, double eps, int co
);
// EScalar f(R)-constraint pipeline: get4ricciscalar -> frfprim -> finalize.
// Caller provides a scratch device buffer of the block's padded extent for RR.
// Cons_fR ends up holding `Sphi - sqrt(3/pi)/4 * log(f'(R))`.
extern "C" void gpu_compute_constraint_fr_launch(
    cudaStream_t stream, int ex[3], int pex[3],
    double* chi, double* trK, double* rho, double* Sphi,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* Sxx, double* Sxy, double* Sxz, double* Syy, double* Syz, double* Szz,
    double* Cons_fR, double* d_scratch_RR);
#endif /* USE_GPU */

#endif /* BSSN_H */
