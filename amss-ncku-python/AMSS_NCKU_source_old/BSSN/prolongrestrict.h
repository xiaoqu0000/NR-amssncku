
#ifndef PROLONGRESTRICT_H
#define PROLONGRESTRICT_H

#ifdef fortran1
#define f_prolong3 prolong3
#define f_prolongmix3 prolongmix3
#define f_prolongcopy3 prolongcopy3
#define f_restrict3 restrict3
#endif

#ifdef fortran2
#define f_prolong3 PROLONG3
#define f_prolongmix3 PROLONGMIX3
#define f_prolongcopy3 PROLONGCOPY3
#define f_restrict3 RESTRICT3
#endif

#ifdef fortran3
#define f_prolong3 prolong3_
#define f_prolongmix3 prolongmix3_
#define f_prolongcopy3 prolongcopy3_
#define f_restrict3 restrict3_
#endif

extern "C"
{
	int f_prolong3(int &, double *, double *, int *, double *,
				   double *, double *, int *, double *,
				   double *, double *, double *, int &);
}

extern "C"
{
	void f_restrict3(int &, double *, double *, int *, double *,
					 double *, double *, int *, double *,
					 double *, double *, double *, int &);
}

extern "C"
{
	int f_prolongmix3(int &, double *, double *, int *, double *,
					  double *, double *, int *, double *,
					  double *, double *, double *, int &,
					  double *, double *);
}

extern "C"
{
	int f_prolongcopy3(int &, double *, double *, int *, double *,
					   double *, double *, int *, double *,
					   double *, double *, double *, int &);
}

#endif /* PROLONGRESTRICT_H */
