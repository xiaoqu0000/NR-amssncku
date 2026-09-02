
#ifndef INITIAL_NULL2_H
#define INITIAL_NULL2_H

#ifdef fortran1
#define f_get_initial_null2 get_initial_null2
#define f_get_initial_null3 get_initial_null3
#define f_get_gauge_g00 get_gauge_g00
#define f_get_gauge_g00_K get_gauge_g00_k
#define f_get_gauge_g00_real get_gauge_g00_real
#define f_get_null_boundary2 get_null_boundary2
#define f_get_null_boundary3 get_null_boundary3
#define f_get_g00_with_t get_g00_with_t
#endif
#ifdef fortran2
#define f_get_initial_null2 GET_INITIAL_NULL2
#define f_get_initial_null3 GET_INITIAL_NULL3
#define f_get_gauge_g00 GET_GAUGE_G00
#define f_get_gauge_g00_K GET_GAUGE_G00_K
#define f_get_gauge_g00_real GET_GAUGE_G00_REAL
#define f_get_null_boundary2 GET_NULL_BOUNDARY2
#define f_get_null_boundary3 GET_NULL_BOUNDARY3
#define f_get_g00_with_t GET_G00_WITH_T
#endif
#ifdef fortran3
#define f_get_initial_null2 get_initial_null2_
#define f_get_initial_null3 get_initial_null3_
#define f_get_gauge_g00 get_gauge_g00_
#define f_get_gauge_g00_K get_gauge_g00_k_
#define f_get_gauge_g00_real get_gauge_g00_real_
#define f_get_null_boundary2 get_null_boundary2_
#define f_get_null_boundary3 get_null_boundary3_
#define f_get_g00_with_t get_g00_with_t_
#endif

extern "C"
{
	void f_get_initial_null2(int *, double *, double *, double *,
							 double *, double *, double *,
							 int &, double &);
}

extern "C"
{
	void f_get_gauge_g00(int *, double *, double *, double *,
						 double *, double *, double *,
						 double *, double *, double *,
						 double *, double &, int &);
}

extern "C"
{
	void f_get_gauge_g00_K(int *, double *, double *, double *,
						   double *, double *, double *,
						   double *, double *, double *,
						   double *, double &);
}

extern "C"
{
	void f_get_gauge_g00_real(int *, double *, double *, double *,
							  double *, double *, double *,
							  double *, double *, double *,
							  double *, double &);
}

extern "C"
{
	void f_get_null_boundary2(int *, double *, double *, double *,
							  double *, double *, double *,
							  double *, double *, double *, double *, double *,
							  double *, double *, double *,
							  double &);
}

extern "C"
{
	void f_get_g00_with_t(double &, int *, double *, double *, double *,
						  double *, double &, int &);
}

extern "C"
{
	void f_get_null_boundary3(double &, int *, double *, double *, double *,
							  double *, double *, double *,
							  double *, double *, double *, double *, double *,
							  double *, double *, double *,
							  double &, int &);
}

extern "C"
{
	void f_get_initial_null3(int *, double *, double *, double *,
							 double *, double *, double *,
							 int &, double &);
}

#endif /* INITIAL_NULL2_H */
