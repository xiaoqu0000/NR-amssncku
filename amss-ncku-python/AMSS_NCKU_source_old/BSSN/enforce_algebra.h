
#ifndef ENFORCE_ALGEBRA_H
#define ENFORCE_ALGEBRA_H

#ifdef fortran1
#define f_enforce_ag enforce_ag
#define f_enforce_ga enforce_ga
#endif
#ifdef fortran2
#define f_enforce_ag ENFORCE_AG
#define f_enforce_ga ENFORCE_GA
#endif
#ifdef fortran3
#define f_enforce_ag enforce_ag_
#define f_enforce_ga enforce_ga_
#endif

extern "C"
{
	void f_enforce_ag(int *,
					  double *, double *, double *, double *, double *, double *,
					  double *, double *, double *, double *, double *, double *);
}
extern "C"
{
	void f_enforce_ga(int *,
					  double *, double *, double *, double *, double *, double *,
					  double *, double *, double *, double *, double *, double *);
}
#endif /* ENFORCE_ALGEBRA_H */
