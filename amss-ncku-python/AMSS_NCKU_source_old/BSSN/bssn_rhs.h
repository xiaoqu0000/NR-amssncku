
#ifndef BSSN_H
#define BSSN_H

#ifdef fortran1
#define f_compute_rhs_bssn compute_rhs_bssn
#define f_compute_rhs_bssn_ss compute_rhs_bssn_ss
#define f_compute_rhs_bssn_escalar compute_rhs_bssn_escalar
#define f_compute_rhs_bssn_escalar_ss compute_rhs_bssn_escalar_ss
#define f_compute_rhs_Z4c compute_rhs_z4c
#define f_compute_rhs_Z4cnot compute_rhs_z4cnot
#define f_compute_rhs_Z4c_ss compute_rhs_z4c_ss
#define f_compute_constraint_fr compute_constraint_fr
#endif
#ifdef fortran2
#define f_compute_rhs_bssn COMPUTE_RHS_BSSN
#define f_compute_rhs_bssn_ss COMPUTE_RHS_BSSN_SS
#define f_compute_rhs_bssn_escalar COMPUTE_RHS_BSSN_ESCALAR
#define f_compute_rhs_bssn_escalar_ss COMPUTE_RHS_BSSN_ESCALAR_SS
#define f_compute_rhs_Z4c COMPUTE_RHS_Z4C
#define f_compute_rhs_Z4cnot COMPUTE_RHS_Z4CNOT
#define f_compute_rhs_Z4c_ss COMPUTE_RHS_Z4C_SS
#define f_compute_constraint_fr COMPUTE_CONSTRAINT_FR
#endif
#ifdef fortran3
#define f_compute_rhs_bssn compute_rhs_bssn_
#define f_compute_rhs_bssn_ss compute_rhs_bssn_ss_
#define f_compute_rhs_bssn_escalar compute_rhs_bssn_escalar_
#define f_compute_rhs_bssn_escalar_ss compute_rhs_bssn_escalar_ss_
#define f_compute_rhs_Z4c compute_rhs_z4c_
#define f_compute_rhs_Z4cnot compute_rhs_z4cnot_
#define f_compute_rhs_Z4c_ss compute_rhs_z4c_ss_
#define f_compute_constraint_fr compute_constraint_fr_
#endif
extern "C"
{
        int f_compute_rhs_bssn(int *, double &, double *, double *, double *,                                                      // ex,T,X,Y,Z
                               double *, double *,                                                                                 // chi, trK
                               double *, double *, double *, double *, double *, double *,                                         // gij
                               double *, double *, double *, double *, double *, double *,                                         // Aij
                               double *, double *, double *,                                                                       // Gam
                               double *, double *, double *, double *, double *, double *, double *,                               // Gauge
                               double *, double *,                                                                                 // chi, trK
                               double *, double *, double *, double *, double *, double *,                                         // gij
                               double *, double *, double *, double *, double *, double *,                                         // Aij
                               double *, double *, double *,                                                                       // Gam
                               double *, double *, double *, double *, double *, double *, double *,                               // Gauge
                               double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, // stress-energy
                               double *, double *, double *, double *, double *, double *,                                         // Christoffel
                               double *, double *, double *, double *, double *, double *,                                         // Christoffel
                               double *, double *, double *, double *, double *, double *,                                         // Christoffel
                               double *, double *, double *, double *, double *, double *,                                         // Ricci
                               double *, double *, double *, double *, double *, double *, double *,                               // constraint violation
                               int &, int &, double &, int &);
}

extern "C"
{
        int f_compute_rhs_bssn_ss(int *, double &, double *, double *, double *,                                                      // ex,T,rho,sigma,R
                                  double *, double *, double *,                                                                       // X,Y,Z
                                  double *, double *, double *,                                                                       // drhodx,drhody,drhodz
                                  double *, double *, double *,                                                                       // dsigmadx,dsigmady,dsigmadz
                                  double *, double *, double *,                                                                       // dRdx,dRdy,dRdz
                                  double *, double *, double *, double *, double *, double *,                                         // drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz
                                  double *, double *, double *, double *, double *, double *,                                         // dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz
                                  double *, double *, double *, double *, double *, double *,                                         // dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz
                                  double *, double *,                                                                                 // chi, trK
                                  double *, double *, double *, double *, double *, double *,                                         // gij
                                  double *, double *, double *, double *, double *, double *,                                         // Aij
                                  double *, double *, double *,                                                                       // Gam
                                  double *, double *, double *, double *, double *, double *, double *,                               // Gauge
                                  double *, double *,                                                                                 // chi, trK
                                  double *, double *, double *, double *, double *, double *,                                         // gij
                                  double *, double *, double *, double *, double *, double *,                                         // Aij
                                  double *, double *, double *,                                                                       // Gam
                                  double *, double *, double *, double *, double *, double *, double *,                               // Gauge
                                  double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, // stress-energy
                                  double *, double *, double *, double *, double *, double *,                                         // Christoffel
                                  double *, double *, double *, double *, double *, double *,                                         // Christoffel
                                  double *, double *, double *, double *, double *, double *,                                         // Christoffel
                                  double *, double *, double *, double *, double *, double *,                                         // Ricci
                                  double *, double *, double *, double *, double *, double *, double *,                               // constraint violation
                                  int &, int &, double &, int &, int &);
}

extern "C"
{
        int f_compute_rhs_bssn_escalar(int *, double &, double *, double *, double *,                                                      // ex,T,X,Y,Z
                                       double *, double *,                                                                                 // chi, trK
                                       double *, double *, double *, double *, double *, double *,                                         // gij
                                       double *, double *, double *, double *, double *, double *,                                         // Aij
                                       double *, double *, double *,                                                                       // Gam
                                       double *, double *, double *, double *, double *, double *, double *,                               // Gauge
                                       double *, double *,                                                                                 // Sphi, Spi
                                       double *, double *,                                                                                 // chi, trK
                                       double *, double *, double *, double *, double *, double *,                                         // gij
                                       double *, double *, double *, double *, double *, double *,                                         // Aij
                                       double *, double *, double *,                                                                       // Gam
                                       double *, double *, double *, double *, double *, double *, double *,                               // Gauge
                                       double *, double *,                                                                                 // Sphi, Spi
                                       double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, // stress-energy
                                       double *, double *, double *, double *, double *, double *,                                         // Christoffel
                                       double *, double *, double *, double *, double *, double *,                                         // Christoffel
                                       double *, double *, double *, double *, double *, double *,                                         // Christoffel
                                       double *, double *, double *, double *, double *, double *,                                         // Ricci
                                       double *, double *, double *, double *, double *, double *, double *,                               // constraint violation
                                       int &, int &, double &, int &);
}

extern "C"
{
        int f_compute_rhs_bssn_escalar_ss(int *, double &, double *, double *, double *,                                                      // ex,T,rho,sigma,R
                                          double *, double *, double *,                                                                       // X,Y,Z
                                          double *, double *, double *,                                                                       // drhodx,drhody,drhodz
                                          double *, double *, double *,                                                                       // dsigmadx,dsigmady,dsigmadz
                                          double *, double *, double *,                                                                       // dRdx,dRdy,dRdz
                                          double *, double *, double *, double *, double *, double *,                                         // drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz
                                          double *, double *, double *, double *, double *, double *,                                         // dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz
                                          double *, double *, double *, double *, double *, double *,                                         // dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz
                                          double *, double *,                                                                                 // chi, trK
                                          double *, double *, double *, double *, double *, double *,                                         // gij
                                          double *, double *, double *, double *, double *, double *,                                         // Aij
                                          double *, double *, double *,                                                                       // Gam
                                          double *, double *, double *, double *, double *, double *, double *,                               // Gauge
                                          double *, double *,                                                                                 // Sphi,Spi
                                          double *, double *,                                                                                 // chi, trK
                                          double *, double *, double *, double *, double *, double *,                                         // gij
                                          double *, double *, double *, double *, double *, double *,                                         // Aij
                                          double *, double *, double *,                                                                       // Gam
                                          double *, double *, double *, double *, double *, double *, double *,                               // Gauge
                                          double *, double *,                                                                                 // Sphi,Spi
                                          double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, // stress-energy
                                          double *, double *, double *, double *, double *, double *,                                         // Christoffel
                                          double *, double *, double *, double *, double *, double *,                                         // Christoffel
                                          double *, double *, double *, double *, double *, double *,                                         // Christoffel
                                          double *, double *, double *, double *, double *, double *,                                         // Ricci
                                          double *, double *, double *, double *, double *, double *, double *,                               // constraint violation
                                          int &, int &, double &, int &, int &);
}

extern "C"
{
        int f_compute_rhs_Z4c(int *, double &, double *, double *, double *,                        // ex,T,X,Y,Z
                              double *, double *,                                                   // chi, trK
                              double *, double *, double *, double *, double *, double *,           // gij
                              double *, double *, double *, double *, double *, double *,           // Aij
                              double *, double *, double *,                                         // Gam
                              double *, double *, double *, double *, double *, double *, double *, // Gauge
                              double *,                                                             // Z4
                              double *, double *,                                                   // chi, trK
                              double *, double *, double *, double *, double *, double *,           // gij
                              double *, double *, double *, double *, double *, double *,           // Aij
                              double *, double *, double *,                                         // Gam
                              double *, double *, double *, double *, double *, double *, double *, // Gauge
                              double *,                                                             // Z4
                              double *, double *, double *, double *, double *, double *, double *, double *, double *, double *,
                              double *, double *, double *, double *, double *, double *,
                              double *, double *, double *, double *, double *, double *,
                              double *, double *, double *, double *, double *, double *,
                              double *, double *, double *, double *, double *, double *,
                              double *, double *, double *, double *, double *, double *, double *,
                              int &, int &, double &, int &);
}

extern "C"
{
        int f_compute_rhs_Z4c_ss(int *, double &, double *, double *, double *,                                                      // ex,T,rho,sigma,R
                                 double *, double *, double *,                                                                       // X,Y,Z
                                 double *, double *, double *,                                                                       // drhodx,drhody,drhodz
                                 double *, double *, double *,                                                                       // dsigmadx,dsigmady,dsigmadz
                                 double *, double *, double *,                                                                       // dRdx,dRdy,dRdz
                                 double *, double *, double *, double *, double *, double *,                                         // drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz
                                 double *, double *, double *, double *, double *, double *,                                         // dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz
                                 double *, double *, double *, double *, double *, double *,                                         // dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz
                                 double *, double *,                                                                                 // chi, trK
                                 double *, double *, double *, double *, double *, double *,                                         // gij
                                 double *, double *, double *, double *, double *, double *,                                         // Aij
                                 double *, double *, double *,                                                                       // Gam
                                 double *, double *, double *, double *, double *, double *, double *,                               // Gauge
                                 double *,                                                                                           // TZ
                                 double *, double *,                                                                                 // chi, trK
                                 double *, double *, double *, double *, double *, double *,                                         // gij
                                 double *, double *, double *, double *, double *, double *,                                         // Aij
                                 double *, double *, double *,                                                                       // Gam
                                 double *, double *, double *, double *, double *, double *, double *,                               // Gauge
                                 double *,                                                                                           // TZ
                                 double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, // stress-energy
                                 double *, double *, double *, double *, double *, double *,                                         // Christoffel
                                 double *, double *, double *, double *, double *, double *,                                         // Christoffel
                                 double *, double *, double *, double *, double *, double *,                                         // Christoffel
                                 double *, double *, double *, double *, double *, double *,                                         // Ricci
                                 double *, double *, double *, double *, double *, double *, double *,                               // constraint violation
                                 int &, int &, double &, int &, int &);
}

extern "C"
{
        int f_compute_rhs_Z4cnot(int *, double &, double *, double *, double *,                        // ex,T,X,Y,Z
                                 double *, double *,                                                   // chi, trK
                                 double *, double *, double *, double *, double *, double *,           // gij
                                 double *, double *, double *, double *, double *, double *,           // Aij
                                 double *, double *, double *,                                         // Gam
                                 double *, double *, double *, double *, double *, double *, double *, // Gauge
                                 double *,                                                             // Z4
                                 double *, double *,                                                   // chi, trK
                                 double *, double *, double *, double *, double *, double *,           // gij
                                 double *, double *, double *, double *, double *, double *,           // Aij
                                 double *, double *, double *,                                         // Gam
                                 double *, double *, double *, double *, double *, double *, double *, // Gauge
                                 double *,                                                             // Z4
                                 double *, double *, double *, double *, double *, double *, double *, double *, double *, double *,
                                 double *, double *, double *, double *, double *, double *,
                                 double *, double *, double *, double *, double *, double *,
                                 double *, double *, double *, double *, double *, double *,
                                 double *, double *, double *, double *, double *, double *,
                                 double *, double *, double *, double *, double *, double *, double *,
                                 int &, int &, double &, int &, double &);
}

extern "C"
{
        void f_compute_constraint_fr(int *, double *, double *, double *,                        // ex,X,Y,Z
                                     double *, double *, double *, double *,                     // chi, trK,rho,Sphi
                                     double *, double *, double *, double *, double *, double *, // gij
                                     double *, double *, double *, double *, double *, double *, // Aij
                                     double *, double *, double *, double *, double *, double *, // Rij
                                     double *, double *, double *, double *, double *, double *, // Sij
                                     double *);
} // FR_cons

#endif /* BSSN_H */
