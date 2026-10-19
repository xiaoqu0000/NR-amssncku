
#ifndef CPBC_H
#define CPBC_H

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
