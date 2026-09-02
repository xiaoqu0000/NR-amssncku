
#ifndef FADMQUANTITES_H
#define FADMQUANTITES_H

#ifdef fortran1
#define f_admmass_bssn admmass_bssn
#define f_admmass_bssn_ss admmass_bssn_ss
#define f_admmomentum_bssn admmomentum_bssn
#endif
#ifdef fortran2
#define f_admmass_bssn ADMMASS_BSSN
#define f_admmass_bssn_ss ADMMASS_BSSN_SS
#define f_admmomentum_bssn ADMMOMENTUM_BSSN
#endif
#ifdef fortran3
#define f_admmass_bssn admmass_bssn_
#define f_admmass_bssn_ss admmass_bssn_ss_
#define f_admmomentum_bssn admmomentum_bssn_
#endif

extern "C"
{
	void f_admmass_bssn(int *, double *, double *, double *,
						double *, double *,
						double *, double *, double *, double *, double *, double *,
						double *, double *, double *, double *, double *, double *,
						double *, double *, double *,
						double *, double *, double *,
						int &);
}

extern "C"
{
	void f_admmass_bssn_ss(int *, double *, double *, double *,
						   double *, double *, double *,
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
						   double *, double *, double *,
						   int &, int &);
}

extern "C"
{
	void f_admmomentum_bssn(int *, double *, double *, double *,
							double *, double *,
							double *, double *, double *, double *, double *, double *,
							double *, double *, double *, double *, double *, double *,
							double *, double *, double *,
							double *, double *, double *,
							double *, double *, double *, double *, double *, double *);
}
#endif /* FADMQUANTITES_H */
