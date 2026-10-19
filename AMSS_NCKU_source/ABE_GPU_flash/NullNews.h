
#ifndef NULLNEWS_H
#define NULLNEWS_H

#ifdef fortran1
#define f_drive_null_news drive_null_news
#define f_get_null_news2 get_null_news2
#define f_drive_null_news_diff drive_null_news_diff
#define f_omega_rhs omega_rhs
#define f_get_exact_omega get_exact_omega
#define f_get_omega_and_dtomega_pre get_omega_and_dtomega_pre
#define f_get_omega_and_dtomega_LN get_omega_and_dtomega_ln
#define f_get_dtomega get_dtomega
#endif
#ifdef fortran2
#define f_drive_null_news DRIVE_NULL_NEWS
#define f_get_null_news2 GET_NULL_NEWS2
#define f_drive_null_news_diff DRIVE_NULL_NEWS_DIFF
#define f_omega_rhs OMEGA_RHS
#define f_get_exact_omega GET_EXACT_OMEGA
#define f_get_omega_and_dtomega_pre GET_OMEGA_AND_DTOMEGA_PRE
#define f_get_omega_and_dtomega_LN GET_OMEGA_AND_DTOMEGA_LN
#define f_get_dtomega GET_DTOMEGA
#endif
#ifdef fortran3
#define f_drive_null_news drive_null_news_
#define f_get_null_news2 get_null_news2_
#define f_drive_null_news_diff drive_null_news_diff_
#define f_omega_rhs omega_rhs_
#define f_get_exact_omega get_exact_omega_
#define f_get_omega_and_dtomega_pre get_omega_and_dtomega_pre_
#define f_get_omega_and_dtomega_LN get_omega_and_dtomega_ln_
#define f_get_dtomega get_dtomega_
#endif

extern "C"
{
	void f_drive_null_news(int *, double *, double *, double *,
						   double *, double *, double *, double *, double *, double *, double *, double *,
						   double *, double *, double *, double *,
						   double *, double *, double *, double *,
						   double *, double *,
						   double *, double *, double *, double *,
						   double *, double *, double *, double *,
						   double *, double *, double *, double *,
						   double *, double *, double &, int &);
}

extern "C"
{
	void f_drive_null_news_diff(int *, double *, double *, double *,
								double *, double *, double *, double *, double *, double *, double *, double *,
								double *, double *, double *, double *,
								double *, double *, double *, double *,
								double *, double *,
								double *, double *, double *, double *,
								double *, double *, double *, double *,
								double *, double *, double *, double *,
								double *, double *, double &, int &, double &);
}

extern "C"
{
	void f_omega_rhs(int *, double *, double *, double *,
					 double *, double *, double *, double *,
					 double *, double *, double *, double *,
					 double *, double *);
}

extern "C"
{
	void f_get_exact_omega(int *, double *, double *, double *,
						   double *,
						   int &, double &, double &);
}

extern "C"
{
	void f_get_null_news2(int *, double *, double *, double *,
						  double *, double *,
						  double *, double *, double *, double *,
						  double *, double *, double *,
						  double *, double *, double *,
						  double *, double *, double &, int &);
}

extern "C"
{
	void f_get_omega_and_dtomega_pre(int *, double *, double *, double *,
									 double *, double *, double *,
									 double *, double *, double &);
}

extern "C"
{
	void f_get_dtomega(int *, double *, double *, double *,
					   double *, double *, double *,
					   double *, double *, double &);
}

extern "C"
{
	void f_get_omega_and_dtomega_LN(double &, int *, double *, double *, double *,
									double *, double *, double &, int &);
}
#endif /* NULLNEWS_H */
