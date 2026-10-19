
#ifndef GET_INITIAL_MAXWELL_H
#define GET_INITIAL_MAXWELL_H

#ifdef fortran1
#define f_get_initial_nbhsem get_initial_nbhsem
#define f_get_initial_nbhsem_ss get_initial_nbhsem_ss
#define f_get_ansorg_nbhs_em get_ansorg_nbhs_em
#define f_get_ansorg_nbhs_ss_em get_ansorg_nbhs_ss_em
#endif
#ifdef fortran2
#define f_get_initial_nbhsem GET_INITIAL_NBHSEM
#define f_get_initial_nbhsem_ss GET_INITIAL_NBHSEM_SS
#define f_get_ansorg_nbhs_em GET_ANSORG_NBHS_EM
#define f_get_ansorg_nbhs_ss_em GET_ANSORG_NBHS_SS_EM
#endif
#ifdef fortran3
#define f_get_initial_nbhsem get_initial_nbhsem_
#define f_get_initial_nbhsem_ss get_initial_nbhsem_ss_
#define f_get_ansorg_nbhs_em get_ansorg_nbhs_em_
#define f_get_ansorg_nbhs_ss_em get_ansorg_nbhs_ss_em_
#endif

extern "C"
{
    void f_get_initial_nbhsem(int *, double *, double *, double *,
                              double *, double *,
                              double *, double *, double *, double *, double *, double *,
                              double *, double *, double *, double *, double *, double *,
                              double *, double *, double *,
                              double *, double *, double *, double *,
                              double *, double *, double *,
                              double *, double *, double *, double *, double *, double *, double *, double *,
                              double *, double *, double *, double *, double *, int &);
}

extern "C"
{
    void f_get_initial_nbhsem_ss(int *, double *, double *, double *,
                                 double *, double *,
                                 double *, double *, double *, double *, double *, double *,
                                 double *, double *, double *, double *, double *, double *,
                                 double *, double *, double *,
                                 double *, double *, double *, double *,
                                 double *, double *, double *,
                                 double *, double *, double *, double *, double *, double *, double *, double *,
                                 double *, double *, double *, double *, double *, int &);
}

extern "C"
{
    void f_get_ansorg_nbhs_em(int *, double *, double *, double *,
                              double *, double *,
                              double *, double *, double *, double *, double *, double *,
                              double *, double *, double *, double *, double *, double *,
                              double *, double *, double *,
                              double *, double *, double *, double *,
                              double *, double *, double *,
                              double *, double *, double *, double *, double *, double *, double *, double *,
                              double *, double *, double *, double *, double *, int &);
}

extern "C"
{
    void f_get_ansorg_nbhs_ss_em(int *, double *, double *, double *,
                                 double *, double *,
                                 double *, double *, double *, double *, double *, double *,
                                 double *, double *, double *, double *, double *, double *,
                                 double *, double *, double *,
                                 double *, double *, double *, double *,
                                 double *, double *, double *,
                                 double *, double *, double *, double *, double *, double *, double *, double *,
                                 double *, double *, double *, double *, double *, int &);
}

#endif /* GET_INITIAL_MAXWELL_H */
