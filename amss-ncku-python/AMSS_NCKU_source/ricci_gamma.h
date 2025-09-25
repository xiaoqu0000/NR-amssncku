
#ifndef RICCI_GAMMA_H
#define RICCI_GAMMA_H

#ifdef fortran1
#define f_ricci_gamma ricci_gamma
#define f_ricci_gamma_ss ricci_gamma_ss
#endif
#ifdef fortran2
#define f_ricci_gamma RICCI_GAMMA
#define f_ricci_gamma_ss RICCI_GAMMA_SS
#endif
#ifdef fortran3
#define f_ricci_gamma ricci_gamma_
#define f_ricci_gamma_ss ricci_gamma_ss_
#endif
extern "C"
{
        void f_ricci_gamma(int *, double *, double *, double *,
                           double *,
                           double *, double *, double *, double *, double *, double *,
                           double *, double *, double *,
                           double *, double *, double *, double *, double *, double *,
                           double *, double *, double *, double *, double *, double *,
                           double *, double *, double *, double *, double *, double *,
                           double *, double *, double *, double *, double *, double *,
                           int &);
}

extern "C"
{
        void f_ricci_gamma_ss(int *, double *, double *, double *, double *, double *, double *,
                              double *, double *, double *,
                              double *, double *, double *,
                              double *, double *, double *,
                              double *, double *, double *, double *, double *, double *,
                              double *, double *, double *, double *, double *, double *,
                              double *, double *, double *, double *, double *, double *,
                              double *,
                              double *, double *, double *, double *, double *, double *,
                              double *, double *, double *,
                              double *, double *, double *, double *, double *, double *,
                              double *, double *, double *, double *, double *, double *,
                              double *, double *, double *, double *, double *, double *,
                              double *, double *, double *, double *, double *, double *,
                              int &, int &, int &);
}
#endif /* RICCI_GAMMA_H */
