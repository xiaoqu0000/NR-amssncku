
#ifndef INITIAL_NULL_H
#define INITIAL_NULL_H

#ifdef fortran1
#define f_get_initial_nbhs_null get_initial_nbhs_null
#define f_get_initial_null get_initial_null
#define f_get_exact_null get_exact_null
#define f_get_exact_null_theta get_exact_null_theta
#define f_get_null_boundary get_null_boundary
#define f_get_null_boundary_c get_null_boundary_c
#define f_get_exact_omegau get_exact_omegau
#endif
#ifdef fortran2
#define f_get_initial_nbhs_null GET_INITIAL_NBHS_NULL
#define f_get_initial_null GET_INITIAL_NULL
#define f_get_exact_null GET_EXACT_NULL
#define f_get_exact_null_theta GET_EXACT_NULL_THETA
#define f_get_null_boundary GET_NULL_BOUNDARY
#define f_get_null_boundary_c GET_NULL_BOUNDARY_C
#define f_get_exact_omegau GET_EXACT_OMEGAU
#endif
#ifdef fortran3
#define f_get_initial_nbhs_null get_initial_nbhs_null_
#define f_get_initial_null get_initial_null_
#define f_get_exact_null get_exact_null_
#define f_get_exact_null_theta get_exact_null_theta_
#define f_get_null_boundary get_null_boundary_
#define f_get_null_boundary_c get_null_boundary_c_
#define f_get_exact_omegau get_exact_omegau_
#endif

extern "C"
{
      void f_get_initial_nbhs_null(int *, double *, double *, double *,
                                   double *, double *, double *,
                                   int &, double &);
}

extern "C"
{
      void f_get_initial_null(int *, double *, double *, double *,
                              double *, double *,
                              int &, double &);
}

extern "C"
{
      void f_get_null_boundary(int *, double *, double *, double *,
                               double *, double *, double *, double *, double *, double *, double *, double *,
                               double *, double *, double *, double *, double *, double *, double *, double *,
                               double *, double *,
                               double *, double *, double *, double *, double *, double *, double *, double *,
                               double *, double *, double *, double *,
                               double &, double &, int &);
}

extern "C"
{
      void f_get_null_boundary_c(int *, double *, double *, double *,
                                 double *, double *, double *, double *, double *, double *, double *, double *,
                                 double *, double *, double *, double *, double *, double *, double *, double *,
                                 double *, double *,
                                 double *, double *, double *, double *, double *, double *, double *, double *,
                                 double *, double *, double *, double *,
                                 double &, double &, int &);
}

extern "C"
{
      void f_get_exact_null(int *, double *, double *, double *,
                            double *, double *, int &, double &, double &,
                            double *, double *, double *, double *, double *, double *, double *, double *,
                            double *, double *,
                            double *, double *, double *, double *, double *, double *, double *, double *,
                            double *, double *, double *, double *);
}

extern "C"
{
      void f_get_exact_null_theta(int *, double *, double *, double *,
                                  double *, double *, int &, double &, double &,
                                  double *, double *, double *, double *, double *, double *, double *, double *,
                                  double *, double *,
                                  double *, double *, double *, double *, double *, double *, double *, double *,
                                  double *, double *, double *, double *);
}

extern "C"
{
      void f_get_exact_omegau(int *, double *, double *, double *,
                              double *,
                              double *, double *, double *, double *, double *, double *, double *, double *,
                              double *, double *,
                              double *, double *, double *, double *, double *, double *, double *, double *,
                              double *, double *, double *, double *,
                              double &, double &, int &);
}

#endif /* INITIAL_NULL_H */
