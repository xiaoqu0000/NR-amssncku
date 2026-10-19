
#ifndef GETNP4_H
#define GETNP4_H

#ifdef fortran1
#define f_getnp4old getnp4old
#define f_getnp4oldscalar getnp4oldscalar
#define f_getnp4oldscalar_ss getnp4oldscalar_ss
#define f_getnp4 getnp4
#define f_getnp4_point getnp4_point
#define f_getnp4_ss getnp4_ss
#define f_getnp4old_ss getnp4old_ss
#define f_getnp4scalar getnp4scalar
#define f_getnp4scalar_ss getnp4scalar_ss
#endif
#ifdef fortran2
#define f_getnp4 GETNP4
#define f_getnp4_point GETNP4_POINT
#define f_getnp4 GETNP4OLD
#define f_getnp4scalar GETNP4OLDSCALAR
#define f_getnp4_ss GETNP4_SS
#define f_getnp4old_ss GETNP4OLD_SS
#define f_getnp4oldscalar_ss GETNP4OLDSCALAR_SS
#define f_getnp4scalar GETNP4SCALAR
#define f_getnp4scalar_ss GETNP4SCALAR_SS
#endif
#ifdef fortran3
#define f_getnp4old getnp4old_
#define f_getnp4_point getnp4_point_
#define f_getnp4oldscalar getnp4oldscalar_
#define f_getnp4oldscalar_ss getnp4oldscalar_ss_
#define f_getnp4 getnp4_
#define f_getnp4_ss getnp4_ss_
#define f_getnp4old_ss getnp4old_ss_
#define f_getnp4scalar getnp4scalar_
#define f_getnp4scalar_ss getnp4scalar_ss_
#endif

extern "C"
{
        void f_getnp4old(int *, double *, double *, double *,
                         double *, double *,
                         double *, double *, double *, double *, double *, double *,
                         double *, double *, double *, double *, double *, double *,
                         double *, double *, double *,
                         double *, double *, double *, double *,
                         double *, double *, int &);
}

extern "C"
{
        void f_getnp4old_ss(int *, double *, double *, double *, double *, double *, double *,
                            double *, double *, double *,
                            double *, double *, double *,
                            double *, double *, double *,
                            double *, double *, double *, double *, double *, double *,
                            double *, double *, double *, double *, double *, double *,
                            double *, double *, double *, double *, double *, double *,
                            double *, double *,
                            double *, double *, double *, double *, double *, double *,
                            double *, double *, double *, double *, double *, double *,
                            double *, double *, double *,
                            double *, double *, double *, double *,
                            double *, double *, int &, int &);
}

extern "C"
{
        void f_getnp4oldscalar(int *, double *, double *, double *,
                               double *, double *, double *,
                               double *, double *, double *, double *, double *, double *,
                               double *, double *, double *, double *, double *, double *,
                               double *, double *, double *,
                               double *, double *, double *, double *,
                               double *, double *, int &);
}

extern "C"
{
        void f_getnp4oldscalar_ss(int *, double *, double *, double *, double *, double *, double *,
                                  double *, double *, double *,
                                  double *, double *, double *,
                                  double *, double *, double *,
                                  double *, double *, double *, double *, double *, double *,
                                  double *, double *, double *, double *, double *, double *,
                                  double *, double *, double *, double *, double *, double *,
                                  double *, double *, double *,
                                  double *, double *, double *, double *, double *, double *,
                                  double *, double *, double *, double *, double *, double *,
                                  double *, double *, double *,
                                  double *, double *, double *, double *,
                                  double *, double *, int &, int &);
}

extern "C"
{
        void f_getnp4(int *, double *, double *, double *,
                      double *, double *,
                      double *, double *, double *, double *, double *, double *,
                      double *, double *, double *, double *, double *, double *,
                      double *, double *, double *, double *, double *, double *,
                      double *, double *, double *, double *, double *, double *,
                      double *, double *, double *, double *, double *, double *,
                      double *, double *, double *, double *, double *, double *,
                      double *, double *, int &);
}

extern "C"
{
        void f_getnp4_point(double &, double &, double &,                               // XYZ
                            double &, double &,                                         // chi,trK
                            double &, double &, double &, double &, double &, double &, // gamma_ij
                            double &, double &, double &, double &, double &, double &, // A_ij
                            double &, double &, double &,                               // chi_i
                            double &, double &, double &,                               // trK_i
                            double &, double &, double &,                               // A_ijk
                            double &, double &, double &,
                            double &, double &, double &,
                            double &, double &, double &,
                            double &, double &, double &,
                            double &, double &, double &,
                            double &, double &, double &, double &, double &, double &, // Gam_ijk
                            double &, double &, double &, double &, double &, double &,
                            double &, double &, double &, double &, double &, double &,
                            double &, double &, double &, double &, double &, double &, // R_ij
                            double &, double &);
}

extern "C"
{
        void f_getnp4_ss(int *, double *, double *, double *, double *, double *, double *,
                         double *, double *, double *,
                         double *, double *, double *,
                         double *, double *, double *,
                         double *, double *, double *, double *, double *, double *,
                         double *, double *, double *, double *, double *, double *,
                         double *, double *, double *, double *, double *, double *,
                         double *, double *,
                         double *, double *, double *, double *, double *, double *,
                         double *, double *, double *, double *, double *, double *,
                         double *, double *, double *, double *, double *, double *,
                         double *, double *, double *, double *, double *, double *,
                         double *, double *, double *, double *, double *, double *,
                         double *, double *, double *, double *, double *, double *,
                         double *, double *, int &, int &);
}

extern "C"
{
        void f_getnp4scalar(int *, double *, double *, double *,
                            double *, double *, double *,
                            double *, double *, double *, double *, double *, double *,
                            double *, double *, double *, double *, double *, double *,
                            double *, double *, double *, double *, double *, double *,
                            double *, double *, double *, double *, double *, double *,
                            double *, double *, double *, double *, double *, double *,
                            double *, double *, double *, double *, double *, double *,
                            double *, double *, int &);
}

extern "C"
{
        void f_getnp4scalar_ss(int *, double *, double *, double *, double *, double *, double *,
                               double *, double *, double *,
                               double *, double *, double *,
                               double *, double *, double *,
                               double *, double *, double *, double *, double *, double *,
                               double *, double *, double *, double *, double *, double *,
                               double *, double *, double *, double *, double *, double *,
                               double *, double *, double *,
                               double *, double *, double *, double *, double *, double *,
                               double *, double *, double *, double *, double *, double *,
                               double *, double *, double *, double *, double *, double *,
                               double *, double *, double *, double *, double *, double *,
                               double *, double *, double *, double *, double *, double *,
                               double *, double *, double *, double *, double *, double *,
                               double *, double *, int &, int &);
}

#endif /* GETNP4_H */
