
#ifndef EMPART_H
#define EMPART_H

#ifdef fortran1
#define f_compute_rhs_empart compute_rhs_empart
#define f_compute_rhs_empart_ss compute_rhs_empart_ss
#endif
#ifdef fortran2
#define f_compute_rhs_empart COMPUTE_RHS_EMPART
#define f_compute_rhs_empart_ss COMPUTE_RHS_EMPART_SS
#endif
#ifdef fortran3
#define f_compute_rhs_empart compute_rhs_empart_
#define f_compute_rhs_empart_ss compute_rhs_empart_ss_
#endif

extern "C"
{
    int f_compute_rhs_empart(int *, double *, double *, double *,
                             double *, double *, double *, double *, double *, double *, double *,
                             double *, double *, double *, double *, double *,
                             double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *,
                             double *, double *, double *, double *, double *, double *, double *, double *,
                             double *, double *, double *, double *, double *, double *, double *, double *, double *, double *,
                             int &, int &, double &);
}

extern "C"
{
    int f_compute_rhs_empart_ss(int *, double *, double *, double *, double *, double *, double *,
                                double *, double *, double *,
                                double *, double *, double *,
                                double *, double *, double *,
                                double *, double *, double *, double *, double *, double *,
                                double *, double *, double *, double *, double *, double *,
                                double *, double *, double *, double *, double *, double *,
                                double *, double *, double *, double *, double *, double *, double *,
                                double *, double *, double *, double *, double *,
                                double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *,
                                double *, double *, double *, double *, double *, double *, double *, double *,
                                double *, double *, double *, double *, double *, double *, double *, double *, double *, double *,
                                int &, int &, double &, int &);
}
#endif /* EMPART_H */
