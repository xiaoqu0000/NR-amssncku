#ifndef SHELLFUNCTIONS_GPU_CUH
#define SHELLFUNCTIONS_GPU_CUH

#include <cmath>

// ============================================================
// Shell coordinate parameters (Cartesian identity)
// ============================================================

__device__ __forceinline__
void d_shellcordpar(double &A, double &B, double &r0, double &eps)
{
    A   = 1.0;
    B   = 0.0;
    r0  = 0.0;
    eps = 1.0;
}

// ============================================================
// Coordinate mapping: R -> Cartesian r
// ============================================================

__device__ __forceinline__
double d_getcartr_scalar(double R_val)
{
    double A, B, r0, eps;
    d_shellcordpar(A, B, r0, eps);
    double f = R_val + B;
    double AAsq = A * A;
    double BBsq = B * B;
    double AAe  = AAsq - BBsq / eps;
    return r0 + (A * f - B * sqrt(AAsq + (f * f - BBsq) / eps)) / AAe;
}

// ============================================================
// dR / d(cartesian r)
// ============================================================

__device__ __forceinline__
double d_getdRdcartr_scalar(double R_val, double cartr_val)
{
    (void)R_val;
    double A, B, r0, eps;
    d_shellcordpar(A, B, r0, eps);
    double cmr = cartr_val - r0;
    return A + B * cmr / sqrt(eps * eps + eps * cmr * cmr);
}

// ============================================================
// d^2 R / d(cartesian r)^2
// ============================================================

__device__ __forceinline__
double d_getdRdcartrcartr_scalar(double R_val, double cartr_val)
{
    (void)R_val;
    double A, B, r0, eps;
    d_shellcordpar(A, B, r0, eps);
    double cmr = cartr_val - r0;
    double den = sqrt(eps + cmr * cmr);
    return B * sqrt(eps) / (den * den * den);
}

// ============================================================
// Face coordinate getters: (rho, sigma, R) -> (x, y, z)
// ============================================================

__device__ __forceinline__
void d_zp_getxyz(double rho_val, double sigma_val, double R_val,
                 double &x, double &y, double &z)
{
    double cartr   = d_getcartr_scalar(R_val);
    double tgrho   = tan(rho_val);
    double tgsigma = tan(sigma_val);
    double invden  = 1.0 / sqrt(1.0 + tgrho * tgrho + tgsigma * tgsigma);
    z = cartr * invden;
    x = z * tgrho;
    y = z * tgsigma;
}

__device__ __forceinline__
void d_zm_getxyz(double rho_val, double sigma_val, double R_val,
                 double &x, double &y, double &z)
{
    double cartr   = d_getcartr_scalar(R_val);
    double tgrho   = tan(rho_val);
    double tgsigma = tan(sigma_val);
    double invden  = 1.0 / sqrt(1.0 + tgrho * tgrho + tgsigma * tgsigma);
    z = -cartr * invden;
    x = z * tgrho;
    y = z * tgsigma;
}

__device__ __forceinline__
void d_yp_getxyz(double rho_val, double sigma_val, double R_val,
                 double &x, double &y, double &z)
{
    double cartr   = d_getcartr_scalar(R_val);
    double tgrho   = tan(rho_val);
    double tgsigma = tan(sigma_val);
    double invden  = 1.0 / sqrt(1.0 + tgrho * tgrho + tgsigma * tgsigma);
    y = cartr * invden;
    x = y * tgrho;
    z = y * tgsigma;
}

__device__ __forceinline__
void d_ym_getxyz(double rho_val, double sigma_val, double R_val,
                 double &x, double &y, double &z)
{
    double cartr   = d_getcartr_scalar(R_val);
    double tgrho   = tan(rho_val);
    double tgsigma = tan(sigma_val);
    double invden  = 1.0 / sqrt(1.0 + tgrho * tgrho + tgsigma * tgsigma);
    y = -cartr * invden;
    x = y * tgrho;
    z = y * tgsigma;
}

__device__ __forceinline__
void d_xp_getxyz(double rho_val, double sigma_val, double R_val,
                 double &x, double &y, double &z)
{
    double cartr   = d_getcartr_scalar(R_val);
    double tgrho   = tan(rho_val);
    double tgsigma = tan(sigma_val);
    double invden  = 1.0 / sqrt(1.0 + tgrho * tgrho + tgsigma * tgsigma);
    x = cartr * invden;
    y = x * tgrho;
    z = x * tgsigma;
}

__device__ __forceinline__
void d_xm_getxyz(double rho_val, double sigma_val, double R_val,
                 double &x, double &y, double &z)
{
    double cartr   = d_getcartr_scalar(R_val);
    double tgrho   = tan(rho_val);
    double tgsigma = tan(sigma_val);
    double invden  = 1.0 / sqrt(1.0 + tgrho * tgrho + tgsigma * tgsigma);
    x = -cartr * invden;
    y = x * tgrho;
    z = x * tgsigma;
}

// ============================================================
// Common helper: second derivatives of R w.r.t. Cartesian coords
// ============================================================

__device__ __forceinline__
void d_compute_dR_second_derivs(
    double x, double y, double z,
    double srt, double dRdcartr, double dRdcartrcartr,
    double &dRdxx, double &dRdxy, double &dRdxz,
    double &dRdyy, double &dRdyz, double &dRdzz)
{
    double srt2 = srt * srt;
    double srt3 = srt2 * srt;
    dRdxx = dRdcartrcartr * x * x / srt2 + dRdcartr * (y * y + z * z) / srt3;
    dRdxy = dRdcartrcartr * x * y / srt2 - dRdcartr * (x * y) / srt3;
    dRdxz = dRdcartrcartr * x * z / srt2 - dRdcartr * (x * z) / srt3;
    dRdyy = dRdcartrcartr * y * y / srt2 + dRdcartr * (x * x + z * z) / srt3;
    dRdyz = dRdcartrcartr * y * z / srt2 - dRdcartr * (y * z) / srt3;
    dRdzz = dRdcartrcartr * z * z / srt2 + dRdcartr * (x * x + y * y) / srt3;
}

// ============================================================
// Jacobian for xp / xm faces
// ============================================================

__device__ __forceinline__
void d_xpm_getjacobian(
    double x, double y, double z, double R_val,
    double &drhodx,  double &drhody,  double &drhodz,
    double &dsigmadx, double &dsigmady, double &dsigmadz,
    double &dRdx,     double &dRdy,     double &dRdz,
    double &drhodxx,  double &drhodxy,  double &drhodxz,
    double &drhodyy,  double &drhodyz,  double &drhodzz,
    double &dsigmadxx, double &dsigmadxy, double &dsigmadxz,
    double &dsigmadyy, double &dsigmadyz, double &dsigmadzz,
    double &dRdxx,    double &dRdxy,    double &dRdxz,
    double &dRdyy,    double &dRdyz,    double &dRdzz)
{
    double cartr_val = d_getcartr_scalar(R_val);
    double dRdcartr     = d_getdRdcartr_scalar(R_val, cartr_val);
    double dRdcartrcartr = d_getdRdcartrcartr_scalar(R_val, cartr_val);

    double xxyy = x * x + y * y;
    double xxzz = x * x + z * z;
    double srt  = sqrt(xxyy + z * z);

    dRdx = x / srt * dRdcartr;
    dRdy = y / srt * dRdcartr;
    dRdz = z / srt * dRdcartr;

    drhodx  = -y / xxyy;
    drhody  =  x / xxyy;
    drhodz  = 0.0;
    dsigmadx = -z / xxzz;
    dsigmady = 0.0;
    dsigmadz =  x / xxzz;

    d_compute_dR_second_derivs(x, y, z, srt, dRdcartr, dRdcartrcartr,
                               dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz);

    double xxyy2 = xxyy * xxyy;
    double xxzz2 = xxzz * xxzz;

    drhodxx  =  2.0 * x * y / xxyy2;
    drhodxy  = (-x * x + y * y) / xxyy2;
    drhodxz  = 0.0;
    drhodyy  = -drhodxx;
    drhodyz  = 0.0;
    drhodzz  = 0.0;

    dsigmadxx =  2.0 * x * z / xxzz2;
    dsigmadxy = 0.0;
    dsigmadxz = (-x * x + z * z) / xxzz2;
    dsigmadyy = 0.0;
    dsigmadyz = 0.0;
    dsigmadzz = -dsigmadxx;
}

// ============================================================
// Jacobian for yp / ym faces
// ============================================================

__device__ __forceinline__
void d_ypm_getjacobian(
    double x, double y, double z, double R_val,
    double &drhodx,  double &drhody,  double &drhodz,
    double &dsigmadx, double &dsigmady, double &dsigmadz,
    double &dRdx,     double &dRdy,     double &dRdz,
    double &drhodxx,  double &drhodxy,  double &drhodxz,
    double &drhodyy,  double &drhodyz,  double &drhodzz,
    double &dsigmadxx, double &dsigmadxy, double &dsigmadxz,
    double &dsigmadyy, double &dsigmadyz, double &dsigmadzz,
    double &dRdxx,    double &dRdxy,    double &dRdxz,
    double &dRdyy,    double &dRdyz,    double &dRdzz)
{
    double cartr_val = d_getcartr_scalar(R_val);
    double dRdcartr     = d_getdRdcartr_scalar(R_val, cartr_val);
    double dRdcartrcartr = d_getdRdcartrcartr_scalar(R_val, cartr_val);

    double xxyy = x * x + y * y;
    double yyzz = y * y + z * z;
    double srt  = sqrt(xxyy + z * z);

    dRdx = x / srt * dRdcartr;
    dRdy = y / srt * dRdcartr;
    dRdz = z / srt * dRdcartr;

    drhodx  =  y / xxyy;
    drhody  = -x / xxyy;
    drhodz  = 0.0;
    dsigmadx = 0.0;
    dsigmady = -z / yyzz;
    dsigmadz =  y / yyzz;

    d_compute_dR_second_derivs(x, y, z, srt, dRdcartr, dRdcartrcartr,
                               dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz);

    double xxyy2 = xxyy * xxyy;
    double yyzz2 = yyzz * yyzz;

    drhodxx  = -2.0 * x * y / xxyy2;
    drhodxy  = (x * x - y * y) / xxyy2;
    drhodxz  = 0.0;
    drhodyy  = -drhodxx;
    drhodyz  = 0.0;
    drhodzz  = 0.0;

    dsigmadxx = 0.0;
    dsigmadxy = 0.0;
    dsigmadxz = 0.0;
    dsigmadyy =  2.0 * y * z / yyzz2;
    dsigmadyz = (-y * y + z * z) / yyzz2;
    dsigmadzz = -dsigmadyy;
}

// ============================================================
// Jacobian for zp / zm faces
// ============================================================

__device__ __forceinline__
void d_zpm_getjacobian(
    double x, double y, double z, double R_val,
    double &drhodx,  double &drhody,  double &drhodz,
    double &dsigmadx, double &dsigmady, double &dsigmadz,
    double &dRdx,     double &dRdy,     double &dRdz,
    double &drhodxx,  double &drhodxy,  double &drhodxz,
    double &drhodyy,  double &drhodyz,  double &drhodzz,
    double &dsigmadxx, double &dsigmadxy, double &dsigmadxz,
    double &dsigmadyy, double &dsigmadyz, double &dsigmadzz,
    double &dRdxx,    double &dRdxy,    double &dRdxz,
    double &dRdyy,    double &dRdyz,    double &dRdzz)
{
    double cartr_val = d_getcartr_scalar(R_val);
    double dRdcartr     = d_getdRdcartr_scalar(R_val, cartr_val);
    double dRdcartrcartr = d_getdRdcartrcartr_scalar(R_val, cartr_val);

    double xxzz = x * x + z * z;
    double yyzz = y * y + z * z;
    double srt  = sqrt(xxzz + y * y);

    dRdx = x / srt * dRdcartr;
    dRdy = y / srt * dRdcartr;
    dRdz = z / srt * dRdcartr;

    drhodx  =  z / xxzz;
    drhody  = 0.0;
    drhodz  = -x / xxzz;
    dsigmadx = 0.0;
    dsigmady =  z / yyzz;
    dsigmadz = -y / yyzz;

    d_compute_dR_second_derivs(x, y, z, srt, dRdcartr, dRdcartrcartr,
                               dRdxx, dRdxy, dRdxz, dRdyy, dRdyz, dRdzz);

    double xxzz2 = xxzz * xxzz;
    double yyzz2 = yyzz * yyzz;

    drhodxx  = -2.0 * x * z / xxzz2;
    drhodxy  = 0.0;
    drhodxz  = (x * x - z * z) / xxzz2;
    drhodyy  = 0.0;
    drhodyz  = 0.0;
    drhodzz  = -drhodxx;

    dsigmadxx = 0.0;
    dsigmadxy = 0.0;
    dsigmadxz = 0.0;
    dsigmadyy = -2.0 * y * z / yyzz2;
    dsigmadyz = (y * y - z * z) / yyzz2;
    dsigmadzz = -dsigmadyy;
}

#endif // SHELLFUNCTIONS_GPU_CUH
