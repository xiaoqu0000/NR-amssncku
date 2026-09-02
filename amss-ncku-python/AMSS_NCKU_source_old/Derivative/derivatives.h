
#ifndef DERIVATIVES
#define DERIVATIVES

#ifdef fortran1
#define f_fderivs fderivs
#define f_fderivs_sh fderivs_sh
#define f_fderivs_shc fderivs_shc
#define f_fdderivs_shc fdderivs_shc
#define f_fdderivs fdderivs
#endif
#ifdef fortran2
#define f_fderivs FDERIVS
#define f_fderivs_sh FDERIVS_SH
#define f_fderivs_shc FDERIVS_SHC
#define f_fdderivs_shc FDDERIVS_SHC
#define f_fdderivs FDDERIVS
#endif
#ifdef fortran3
#define f_fderivs fderivs_
#define f_fderivs_sh fderivs_sh_
#define f_fderivs_shc fderivs_shc_
#define f_fdderivs_shc fdderivs_shc_
#define f_fdderivs fdderivs_
#endif

extern "C"
{
	void f_fderivs(int *, double *,
				   double *, double *, double *,
				   double *, double *, double *,
				   double &, double &, double &, int &, int &);
}

extern "C"
{
	void f_fderivs_sh(int *, double *,
					  double *, double *, double *,
					  double *, double *, double *,
					  double &, double &, double &, int &, int &, int &);
}

extern "C"
{
	void f_fderivs_shc(int *, double *,
					   double *, double *, double *,
					   double *, double *, double *,
					   double &, double &, double &, int &, int &, int &,
					   double *, double *, double *,
					   double *, double *, double *,
					   double *, double *, double *);
}

extern "C"
{
	void f_fdderivs_shc(int *, double *,
						double *, double *, double *, double *, double *, double *,
						double *, double *, double *,
						double &, double &, double &, int &, int &, int &,
						double *, double *, double *,
						double *, double *, double *,
						double *, double *, double *,
						double *, double *, double *, double *, double *, double *,
						double *, double *, double *, double *, double *, double *,
						double *, double *, double *, double *, double *, double *);
}

extern "C"
{
	void f_fdderivs(int *, double *,
					double *, double *, double *, double *, double *, double *,
					double *, double *, double *,
					double &, double &, double &, int &, int &);
}

#endif /* DERIVATIVES */
