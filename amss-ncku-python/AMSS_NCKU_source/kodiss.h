
#ifndef KODISS_H
#define KODISS_H

#ifdef fortran1
#define f_kodis_sh kodis_sh
#define f_kodis_shcr kodis_shcr
#define f_kodis_shor kodis_shor
#endif
#ifdef fortran2
#define f_kodis_sh KODIS_SH
#define f_kodis_shcr KODIS_SHCR
#define f_kodis_shor KODIS_SHOR
#endif
#ifdef fortran3
#define f_kodis_sh kodis_sh_
#define f_kodis_shcr kodis_shcr_
#define f_kodis_shor kodis_shor_
#endif

extern "C"
{
    void f_kodis_sh(int *, double *, double *, double *,
                    double *, double *,
                    double *, int &, double &, int &);
}

extern "C"
{
    void f_kodis_shcr(int *, double *, double *, double *,
                      double *, double *,
                      double *, int &, double &, int &);
}

extern "C"
{
    void f_kodis_shor(int *, double *, double *, double *,
                      double *, double *,
                      double *, int &, double &, int &);
}

#endif /* KODISS_H */
