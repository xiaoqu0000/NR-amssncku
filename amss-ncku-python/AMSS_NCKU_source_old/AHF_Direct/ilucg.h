
#ifndef ILUCG_H
#define ILUCG_H

#ifdef fortran1
#define f_ilucg ilucg
#endif
#ifdef fortran2
#define f_ilucg ILUCG
#endif
#ifdef fortran3
#define f_ilucg ilucg_
#endif

extern "C"
{
	void f_ilucg(const int &N,
				 const int *IA, const int *JA, const double *A,
				 const double *B, double *X,
				 int *ITEMP, double *RTEMP,
				 const double &EPS, const int &ITER, int &ISTATUS);
}

#endif /* ILUCG_H */
