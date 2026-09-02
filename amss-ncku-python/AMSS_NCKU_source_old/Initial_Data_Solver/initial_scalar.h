
#ifndef GET_INITIAL_SCALAR_H
#define GET_INITIAL_SCALAR_H

#ifdef fortran1
#define f_get_initial_scalar get_initial_scalar
#define f_get_initial_scalar_sh get_initial_scalar_sh
#endif
#ifdef fortran2
#define f_get_initial_scalar GET_INITIAL_SCALAR
#define f_get_initial_scalar_sh GET_INITIAL_SCALAR_SH
#endif
#ifdef fortran3
#define f_get_initial_scalar get_initial_scalar_
#define f_get_initial_scalar_sh get_initial_scalar_sh_
#endif

extern "C"
{
    void f_get_initial_scalar(int *, double *, double *, double *,
                              double *, double *,
                              double &, double &, double &);
}

extern "C"
{
    void f_get_initial_scalar_sh(int *, double *, double *, double *,
                                 double *, double *,
                                 double &, double &, double &);
}
#endif /* GET_INITIAL_SCALAR_H */
