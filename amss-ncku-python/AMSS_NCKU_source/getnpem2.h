
#ifndef GETNPEM2_H
#define GETNPEM2_H

#ifdef fortran1
#define f_getnpem2 getnpem2
#define f_getnpem2_point getnpem2_point
#define f_getnpem1_point getnpem1_point
#define f_getnpem2_ss getnpem2_ss
#define f_getnpem1 getnpem1
#define f_getnpem1_ss getnpem1_ss
#endif
#ifdef fortran2
#define f_getnpem2 GETNPEM2
#define f_getnpem2_point GETNPEM2_POINT
#define f_getnpem1_point GETNPEM1_POINT
#define f_getnpem2_ss GETNPEM2_SS
#define f_getnpem1 GETNPEM1
#define f_getnpem1_ss GETNPEM1_SS
#endif
#ifdef fortran3
#define f_getnpem2 getnpem2_
#define f_getnpem2_point getnpem2_point_
#define f_getnpem1_point getnpem1_point_
#define f_getnpem2_ss getnpem2_ss_
#define f_getnpem1 getnpem1_
#define f_getnpem1_ss getnpem1_ss_
#endif

extern "C"
{
        void f_getnpem2(int *, double *, double *, double *,
                        double *, double *, double *, double *, double *, double *, double *,
                        double *, double *, double *, double *, double *, double *,
                        double *, double *, int &);
}

extern "C"
{
        void f_getnpem2_point(double &, double &, double &,
                              double &, double &, double &, double &, double &, double &, double &,
                              double &, double &, double &, double &, double &, double &,
                              double &, double &);
}

extern "C"
{
        void f_getnpem2_ss(int *, double *, double *, double *, double *, double *, double *,
                           double *, double *, double *,
                           double *, double *, double *,
                           double *, double *, double *,
                           double *, double *, double *, double *, double *, double *,
                           double *, double *, double *, double *, double *, double *,
                           double *, double *, double *, double *, double *, double *,
                           double *, double *, double *, double *, double *, double *, double *,
                           double *, double *, double *, double *, double *, double *,
                           double *, double *, int &, int &);
}

extern "C"
{
        void f_getnpem1(int *, double *, double *, double *,
                        double *, double *, double *, double *, double *, double *, double *,
                        double *, double *, double *, double *, double *, double *,
                        double *, double *, int &);
}

extern "C"
{
        void f_getnpem1_ss(int *, double *, double *, double *, double *, double *, double *,
                           double *, double *, double *,
                           double *, double *, double *,
                           double *, double *, double *,
                           double *, double *, double *, double *, double *, double *,
                           double *, double *, double *, double *, double *, double *,
                           double *, double *, double *, double *, double *, double *,
                           double *, double *, double *, double *, double *, double *, double *,
                           double *, double *, double *, double *, double *, double *,
                           double *, double *, int &, int &);
}

extern "C"
{
        void f_getnpem1_point(double &, double &, double &,
                              double &, double &, double &, double &, double &, double &, double &,
                              double &, double &, double &, double &, double &, double &,
                              double &, double &);
}

#endif /* GETNPEM2_H */
