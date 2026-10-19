
#ifndef CPBC_H
#define CPBC_H

#ifdef USE_GPU
#include <cuda_runtime.h>
#endif

#ifdef fortran1
#define f_david_milton_extroplate_ss david_milton_extroplate_ss
#define f_david_milton_cpbc_ss david_milton_cpbc_ss
#endif
#ifdef fortran2
#define f_david_milton_extroplate_ss DAVID_MILTON_EXTROPLATE_SS
#define f_david_milton_cpbc_ss DAVID_MILTON_CPBC_SS
#endif
#ifdef fortran3
#define f_david_milton_extroplate_ss david_milton_extroplate_ss_
#define f_david_milton_cpbc_ss david_milton_cpbc_ss_
#endif
extern "C"
{
	int f_david_milton_extroplate_ss(int *, double *, double *, double *,								   // ex,crho,sigma,R
									 double *, double *, double *,										   // TZ, chi, trK
									 double *, double *, double *, double *, double *, double *,		   // gij
									 double *, double *, double *, double *, double *, double *,		   // Aij
									 double *, double *, double *,										   // Gam
									 double *, double *, double *, double *, double *, double *, double *, // Gauge
									 double &, double &);
} // zmin,zmax

#ifdef USE_GPU
// GPU launch functions
void gpu_david_milton_cpbc_ss_launch(
    cudaStream_t &stream,
    int ex[3], int pex[3],
    double *crho, double *sigma, double *R,
    double *x, double *y, double *z,
    double *drhodx, double *drhody, double *drhodz,
    double *dsigmadx, double *dsigmady, double *dsigmadz,
    double *dRdx, double *dRdy, double *dRdz,
    double *drhodxx, double *drhodxy, double *drhodxz, double *drhodyy, double *drhodyz, double *drhodzz,
    double *dsigmadxx, double *dsigmadxy, double *dsigmadxz, double *dsigmadyy, double *dsigmadyz, double *dsigmadzz,
    double *dRdxx, double *dRdxy, double *dRdxz, double *dRdyy, double *dRdyz, double *dRdzz,
    double xmin, double ymin, double zmin, double xmax, double ymax, double zmax,
    double *TZ, double *chi, double *trK,
    double *dxx, double *gxy, double *gxz, double *dyy, double *gyz, double *dzz,
    double *Axx, double *Axy, double *Axz, double *Ayy, double *Ayz, double *Azz,
    double *Gamx, double *Gamy, double *Gamz,
    double *Lap, double *betax, double *betay, double *betaz,
    double *dtSfx, double *dtSfy, double *dtSfz,
    double *TZ_rhs, double *chi_rhs, double *trK_rhs,
    double *gxx_rhs, double *gxy_rhs, double *gxz_rhs,
    double *gyy_rhs, double *gyz_rhs, double *gzz_rhs,
    double *Axx_rhs, double *Axy_rhs, double *Axz_rhs,
    double *Ayy_rhs, double *Ayz_rhs, double *Azz_rhs,
    double *Gamx_rhs, double *Gamy_rhs, double *Gamz_rhs,
    double *Lap_rhs, double *betax_rhs, double *betay_rhs, double *betaz_rhs,
    double *dtSfx_rhs, double *dtSfy_rhs, double *dtSfz_rhs,
    double *pGamxxx, double *pGamxxy, double *pGamxxz,
    double *pGamxyy, double *pGamxyz, double *pGamxzz,
    double *pGamyxx, double *pGamyxy, double *pGamyxz,
    double *pGamyyy, double *pGamyyz, double *pGamyzz,
    double *pGamzxx, double *pGamzxy, double *pGamzxz,
    double *pGamzyy, double *pGamzyz, double *pGamzzz,
    double *Rxx, double *Rxy, double *Rxz,
    double *Ryy, double *Ryz, double *Rzz,
    double *Gmxcon, double *Gmycon, double *Gmzcon,
    int Symmetry, double eps, int sst,
    bool has_outer);

void gpu_david_milton_extroplate_ss_launch(
    cudaStream_t &stream,
    int ex[3], int pex[3],
    double *crho, double *sigma, double *R,
    double *TZ, double *chi, double *trK,
    double *dxx, double *gxy, double *gxz,
    double *dyy, double *gyz, double *dzz,
    double *Axx, double *Axy, double *Axz,
    double *Ayy, double *Ayz, double *Azz,
    double *Gmx, double *Gmy, double *Gmz,
    double *Lap, double *betax, double *betay, double *betaz,
    double *dtSfx, double *dtSfy, double *dtSfz,
    double zmin, double zmax,
    bool has_outer);
#endif /* USE_GPU */

extern "C"
{
	int f_david_milton_cpbc_ss(int *, double *, double *, double *,									 // ex,crho,sigma,R
							   double *, double *, double *,										 // x,y,z
							   double *, double *, double *,										 // drhodx,drhody,drhodz
							   double *, double *, double *,										 // dsigmadx,dsigmady,dsigmadz
							   double *, double *, double *,										 // dRdx,dRdy,dRdz
							   double *, double *, double *, double *, double *, double *,			 // drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz
							   double *, double *, double *, double *, double *, double *,			 // dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz
							   double *, double *, double *, double *, double *, double *,			 // dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz
							   double &, double &, double &, double &, double &, double &,			 // xmin,ymin,zmin,xmax,ymax,zmax
							   double *, double *, double *,										 // TZ,chi, trK
							   double *, double *, double *, double *, double *, double *,			 // gij
							   double *, double *, double *, double *, double *, double *,			 // Aij
							   double *, double *, double *,										 // Gam
							   double *, double *, double *, double *, double *, double *, double *, // Gauge
							   double *, double *, double *,										 // TZ, chi, trK
							   double *, double *, double *, double *, double *, double *,			 // gij
							   double *, double *, double *, double *, double *, double *,			 // Aij
							   double *, double *, double *,										 // Gam
							   double *, double *, double *, double *, double *, double *, double *, // Gauge
							   double *, double *, double *, double *, double *, double *,			 // Christoffel
							   double *, double *, double *, double *, double *, double *,			 // Christoffel
							   double *, double *, double *, double *, double *, double *,			 // Christoffel
							   double *, double *, double *, double *, double *, double *,			 // Ricci
							   double *, double *, double *,										 // Gama constraint
							   int &, double &, int &);
} // Symmetry, eps, sst
#endif /* CPBC_H */
