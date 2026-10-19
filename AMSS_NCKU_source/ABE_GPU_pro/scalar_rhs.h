
#ifndef SCALAR_RHS_H
#define SCALAR_RHS_H

#ifdef fortran1
#define f_compute_rhs_scalar compute_rhs_scalar
#define f_compute_rhs_scalar_ss compute_rhs_scalar_ss
#endif
#ifdef fortran2
#define f_compute_rhs_scalar COMPUTE_RHS_SCALAR
#define f_compute_rhs_scalar_ss COMPUTE_RHS_SCALAR_SS
#endif
#ifdef fortran3
#define f_compute_rhs_scalar compute_rhs_scalar_
#define f_compute_rhs_scalar_ss compute_rhs_scalar_ss_
#endif
extern "C"
{
	int f_compute_rhs_scalar(int *, double &, double *, double *, double *, // ex,T,X,Y,Z
							 double *, double *,							// Sphi,Spi
							 double *, double *,							// Sphi_rhs,Spi_rhs
							 int &, int &, double &);
}

extern "C"
{
	int f_compute_rhs_scalar_ss(int *, double &, double *, double *, double *,				// ex,T,rho,sigma,R
								double *, double *, double *,								// X,Y,Z
								double *, double *, double *,								// drhodx,drhody,drhodz
								double *, double *, double *,								// dsigmadx,dsigmady,dsigmadz
								double *, double *, double *,								// dRdx,dRdy,dRdz
								double *, double *, double *, double *, double *, double *, // drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz
								double *, double *, double *, double *, double *, double *, // dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz
								double *, double *, double *, double *, double *, double *, // dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz
								double *, double *,											// Sphi,Spi
								double *, double *,											// Sphi_rhs,Spi_rhs
								int &, int &, double &, int &);
}

#ifdef USE_GPU
#include <cuda_runtime_api.h>

void gpu_compute_rhs_scalar_launch(cudaStream_t& stream, int* ex, int* pex, double T,
								   double* X, double* Y, double* Z,
								   double* Sphi, double* Spi,
								   double* Sphi_rhs, double* Spi_rhs,
								   int symmetry, int lev, double eps);

void gpu_compute_rhs_scalar_ss_launch(cudaStream_t& stream, int* ex, int* pex, double T,
									  double* crho, double* sigma, double* R,
									  double* x, double* y, double* z,
									  double* drhodx, double* drhody, double* drhodz,
									  double* dsigmadx, double* dsigmady, double* dsigmadz,
									  double* dRdx, double* dRdy, double* dRdz,
									  double* drhodxx, double* drhodxy, double* drhodxz,
									  double* drhodyy, double* drhodyz, double* drhodzz,
									  double* dsigmadxx, double* dsigmadxy, double* dsigmadxz,
									  double* dsigmadyy, double* dsigmadyz, double* dsigmadzz,
									  double* dRdxx, double* dRdxy, double* dRdxz,
									  double* dRdyy, double* dRdyz, double* dRdzz,
									  double* Sphi, double* Spi,
									  double* Sphi_rhs, double* Spi_rhs,
									  int symmetry, int lev, double eps, int sst);
#endif /* USE_GPU */

#endif /* SCALAR_RHS_H */
