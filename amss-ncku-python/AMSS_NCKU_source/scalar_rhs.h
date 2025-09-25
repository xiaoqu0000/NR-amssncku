
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
#endif /* SCALAR_RHS_H */
