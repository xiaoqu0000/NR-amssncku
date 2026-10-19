#ifndef GETNPEM2_GPU_CUH
#define GETNPEM2_GPU_CUH

#include <cmath>

// ---------------------------------------------------------
// Device helper: metric inversion, tetrad initialization,
// and Gram-Schmidt orthonormalization.
// Shared by both Phi2 and Phi1 computation.
// ---------------------------------------------------------
template<int TETRAD>
__device__ __forceinline__ void d_setup_em_tetrad(
    double X, double Y, double Z,
    double dxx, double gxy, double gxz, double dyy, double gyz, double dzz,
    double chi,
    double& vx, double& vy, double& vz,
    double& ux, double& uy, double& uz,
    double& wx, double& wy, double& wz,
    double& gupxx, double& gupxy, double& gupxz,
    double& gupyy, double& gupyz, double& gupzz
)
{
    const double ZEO = 0.0;
    const double ONE = 1.0;
    const double TWO = 2.0;
    const double TINYRR = 1e-14;

    // Physical metric
    double l_gxx = dxx + ONE;
    double l_gyy = dyy + ONE;
    double l_gzz = dzz + ONE;
    double l_gxy = gxy;
    double l_gxz = gxz;
    double l_gyz = gyz;

    double chipn1 = chi + ONE;
    double f = ONE / chipn1;

    // Metric inversion
    double detG = l_gxx * l_gyy * l_gzz + l_gxy * l_gyz * l_gxz + l_gxz * l_gxy * l_gyz -
                  l_gxz * l_gyy * l_gxz - l_gxy * l_gxy * l_gzz - l_gxx * l_gyz * l_gyz;

    gupxx =  (l_gyy * l_gzz - l_gyz * l_gyz) / detG;
    gupxy = -(l_gxy * l_gzz - l_gyz * l_gxz) / detG;
    gupxz =  (l_gxy * l_gyz - l_gyy * l_gxz) / detG;
    gupyy =  (l_gxx * l_gzz - l_gxz * l_gxz) / detG;
    gupyz = -(l_gxx * l_gyz - l_gxy * l_gxz) / detG;
    gupzz =  (l_gxx * l_gyy - l_gxy * l_gxy) / detG;

    // Initialize tetrad vectors
    // v: radial direction, u: phi direction, w: theta direction
    if (fabs(X) < TINYRR && fabs(Y) < TINYRR && fabs(Z) < TINYRR) {
        vx = TINYRR; vy = TINYRR; vz = TINYRR;
    } else {
        vx = X; vy = Y; vz = Z;
    }

    if (fabs(X) < TINYRR && fabs(Y) < TINYRR) {
        ux = -TINYRR; uy = TINYRR; uz = ZEO;
        wx = TINYRR * Z; wy = TINYRR * Z; wz = -TWO * TINYRR * TINYRR;
    } else {
        ux = -Y; uy = X; uz = ZEO;
        wx = X * Z; wy = Y * Z; wz = -(X * X + Y * Y);
    }

    // Gram-Schmidt orthonormalization (order depends on tetradtype)
    if constexpr (TETRAD == 2) {
        // Raise v with inverse metric, then same GS as TETRAD=0
        double fxt = vx, fyt = vy, fzt = vz;
        vx = gupxx * fxt + gupxy * fyt + gupxz * fzt;
        vy = gupxy * fxt + gupyy * fyt + gupyz * fzt;
        vz = gupxz * fxt + gupyz * fyt + gupzz * fzt;

        double temp;
        temp = l_gxx*vx*vx + l_gyy*vy*vy + l_gzz*vz*vz + (l_gxy*vx*vy + l_gxz*vx*vz + l_gyz*vy*vz)*TWO;
        temp = sqrt(temp*f);
        vx /= temp; vy /= temp; vz /= temp;

        temp = l_gxx*vx*ux + l_gxy*vx*uy + l_gxz*vx*uz +
               l_gxy*vy*ux + l_gyy*vy*uy + l_gyz*vy*uz +
               l_gxz*vz*ux + l_gyz*vz*uy + l_gzz*vz*uz;
        temp *= f;
        ux -= temp*vx; uy -= temp*vy; uz -= temp*vz;

        temp = l_gxx*ux*ux + l_gyy*uy*uy + l_gzz*uz*uz + (l_gxy*ux*uy + l_gxz*ux*uz + l_gyz*uy*uz)*TWO;
        temp = sqrt(temp*f);
        ux /= temp; uy /= temp; uz /= temp;

        temp = l_gxx*vx*wx + l_gxy*vx*wy + l_gxz*vx*wz +
               l_gxy*vy*wx + l_gyy*vy*wy + l_gyz*vy*wz +
               l_gxz*vz*wx + l_gyz*vz*wy + l_gzz*vz*wz;
        temp *= f;
        wx -= temp*vx; wy -= temp*vy; wz -= temp*vz;

        temp = l_gxx*ux*wx + l_gxy*ux*wy + l_gxz*ux*wz +
               l_gxy*uy*wx + l_gyy*uy*wy + l_gyz*uy*wz +
               l_gxz*uz*wx + l_gyz*uz*wy + l_gzz*uz*wz;
        temp *= f;
        wx -= temp*ux; wy -= temp*uy; wz -= temp*uz;

        temp = l_gxx*wx*wx + l_gyy*wy*wy + l_gzz*wz*wz + (l_gxy*wx*wy + l_gxz*wx*wz + l_gyz*wy*wz)*TWO;
        temp = sqrt(temp*f);
        wx /= temp; wy /= temp; wz /= temp;

    } else if constexpr (TETRAD == 0) {
        // Normalize v, project out v from u, normalize u,
        // project out v/u from w, normalize w
        double temp;
        temp = l_gxx*vx*vx + l_gyy*vy*vy + l_gzz*vz*vz + (l_gxy*vx*vy + l_gxz*vx*vz + l_gyz*vy*vz)*TWO;
        temp = sqrt(temp*f);
        vx /= temp; vy /= temp; vz /= temp;

        temp = l_gxx*vx*ux + l_gxy*vx*uy + l_gxz*vx*uz +
               l_gxy*vy*ux + l_gyy*vy*uy + l_gyz*vy*uz +
               l_gxz*vz*ux + l_gyz*vz*uy + l_gzz*vz*uz;
        temp *= f;
        ux -= temp*vx; uy -= temp*vy; uz -= temp*vz;

        temp = l_gxx*ux*ux + l_gyy*uy*uy + l_gzz*uz*uz + (l_gxy*ux*uy + l_gxz*ux*uz + l_gyz*uy*uz)*TWO;
        temp = sqrt(temp*f);
        ux /= temp; uy /= temp; uz /= temp;

        temp = l_gxx*vx*wx + l_gxy*vx*wy + l_gxz*vx*wz +
               l_gxy*vy*wx + l_gyy*vy*wy + l_gyz*vy*wz +
               l_gxz*vz*wx + l_gyz*vz*wy + l_gzz*vz*wz;
        temp *= f;
        wx -= temp*vx; wy -= temp*vy; wz -= temp*vz;

        temp = l_gxx*ux*wx + l_gxy*ux*wy + l_gxz*ux*wz +
               l_gxy*uy*wx + l_gyy*uy*wy + l_gyz*uy*wz +
               l_gxz*uz*wx + l_gyz*uz*wy + l_gzz*uz*wz;
        temp *= f;
        wx -= temp*ux; wy -= temp*uy; wz -= temp*uz;

        temp = l_gxx*wx*wx + l_gyy*wy*wy + l_gzz*wz*wz + (l_gxy*wx*wy + l_gxz*wx*wz + l_gyz*wy*wz)*TWO;
        temp = sqrt(temp*f);
        wx /= temp; wy /= temp; wz /= temp;

    } else if constexpr (TETRAD == 1) {
        // Normalize w, project out w from u, normalize u,
        // project out w/u from v, normalize v
        double temp;

        temp = l_gxx*wx*wx + l_gyy*wy*wy + l_gzz*wz*wz + (l_gxy*wx*wy + l_gxz*wx*wz + l_gyz*wy*wz)*TWO;
        temp = sqrt(temp*f);
        wx /= temp; wy /= temp; wz /= temp;

        temp = l_gxx*wx*ux + l_gxy*wx*uy + l_gxz*wx*uz +
               l_gxy*wy*ux + l_gyy*wy*uy + l_gyz*wy*uz +
               l_gxz*wz*ux + l_gyz*wz*uy + l_gzz*wz*uz;
        temp *= f;
        ux -= temp*wx; uy -= temp*wy; uz -= temp*wz;

        temp = l_gxx*ux*ux + l_gyy*uy*uy + l_gzz*uz*uz + (l_gxy*ux*uy + l_gxz*ux*uz + l_gyz*uy*uz)*TWO;
        temp = sqrt(temp*f);
        ux /= temp; uy /= temp; uz /= temp;

        temp = l_gxx*vx*wx + l_gxy*vx*wy + l_gxz*vx*wz +
               l_gxy*vy*wx + l_gyy*vy*wy + l_gyz*vy*wz +
               l_gxz*vz*wx + l_gyz*vz*wy + l_gzz*vz*wz;
        temp *= f;
        vx -= temp*wx; vy -= temp*wy; vz -= temp*wz;

        temp = l_gxx*ux*vx + l_gxy*ux*vy + l_gxz*ux*vz +
               l_gxy*uy*vx + l_gyy*uy*vy + l_gyz*uy*vz +
               l_gxz*uz*vx + l_gyz*uz*vy + l_gzz*uz*vz;
        temp *= f;
        vx -= temp*ux; vy -= temp*uy; vz -= temp*uz;

        temp = l_gxx*vx*vx + l_gyy*vy*vy + l_gzz*vz*vz + (l_gxy*vx*vy + l_gxz*vx*vz + l_gyz*vy*vz)*TWO;
        temp = sqrt(temp*f);
        vx /= temp; vy /= temp; vz /= temp;
    }
}

// ---------------------------------------------------------
// Device point function for Phi1 computation
// ---------------------------------------------------------
template<int TETRAD>
__device__ __forceinline__ void d_getnpem1_point(
    double X, double Y, double Z,
    double chi,
    double dxx, double gxy, double gxz, double dyy, double gyz, double dzz,
    double Ex, double Ey, double Ez,
    double Bx, double By, double Bz,
    double& Rphi1, double& Iphi1
)
{
    const double ONE = 1.0;

    double vx, vy, vz, ux, uy, uz, wx, wy, wz;
    double gupxx, gupxy, gupxz, gupyy, gupyz, gupzz;

    d_setup_em_tetrad<TETRAD>(
        X, Y, Z,
        dxx, gxy, gxz, dyy, gyz, dzz,
        chi,
        vx, vy, vz, ux, uy, uz, wx, wy, wz,
        gupxx, gupxy, gupxz, gupyy, gupyz, gupzz
    );

    double chipn1 = chi + ONE;

    // Physical metric (for lowering E)
    double l_gxx = dxx + ONE;
    double l_gyy = dyy + ONE;
    double l_gzz = dzz + ONE;
    double l_gxy = gxy;
    double l_gxz = gxz;
    double l_gyz = gyz;

    // Lower E with physical metric, scale by chipn1
    double HEx = (l_gxx*Ex + l_gxy*Ey + l_gxz*Ez) * chipn1;
    double HEy = (l_gxy*Ex + l_gyy*Ey + l_gyz*Ez) * chipn1;
    double HEz = (l_gxz*Ex + l_gyz*Ey + l_gzz*Ez) * chipn1;

    // Real part: HE dot v
    Rphi1 = HEx*vx + HEy*vy + HEz*vz;

    // f = chipn1^(-3/2) for EM field scaling
    double t = sqrt(ONE / chipn1);
    double f_em = t * t * t;

    // HB = f * (u x w)
    double HBx = (uy*wz - uz*wy) * f_em;
    double HBy = (uz*wx - ux*wz) * f_em;
    double HBz = (ux*wy - uy*wx) * f_em;

    // Imaginary part: HB dot B
    Iphi1 = HBx*Bx + HBy*By + HBz*Bz;

    Rphi1 /= 2.0;
    Iphi1 /= 2.0;
}

__device__ inline void getnpem2_point_device(
    double X, double Y, double Z,
    double chi, double dxx, double gxy, double gxz, double dyy, double gyz, double dzz,
    double Ex, double Ey, double Ez, double Bx, double By, double Bz,
    double &Rphi2, double &Iphi2
) {
    const double TINYRR = 1.0e-14;

    double gxx = dxx + 1.0;
    double gyy = dyy + 1.0;
    double gzz = dzz + 1.0;
    double chipn1 = chi + 1.0;

    // invert tilted metric
    double gupzz_denom = gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz -
                         gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz;
    double gupxx =  (gyy * gzz - gyz * gyz) / gupzz_denom;
    double gupxy = -(gxy * gzz - gyz * gxz) / gupzz_denom;
    double gupxz =  (gxy * gyz - gyy * gxz) / gupzz_denom;
    double gupyy =  (gxx * gzz - gxz * gxz) / gupzz_denom;
    double gupyz = -(gxx * gyz - gxy * gxz) / gupzz_denom;
    double gupzz =  (gxx * gyy - gxy * gxy) / gupzz_denom;

    // initialize U, V, W vectors (tetradtype == 2)
    double vx, vy, vz, ux, uy, uz, wx, wy, wz;

    if (fabs(X) < TINYRR && fabs(Y) < TINYRR && fabs(Z) < TINYRR) {
        vx = TINYRR; vy = TINYRR; vz = TINYRR;
    } else {
        vx = X; vy = Y; vz = Z;
    }

    if (fabs(X) < TINYRR && fabs(Y) < TINYRR) {
        ux = -TINYRR; uy = TINYRR; uz = 0.0;
        wx = TINYRR * Z; wy = TINYRR * Z; wz = -2.0 * TINYRR * TINYRR;
    } else {
        ux = -Y; uy = X; uz = 0.0;
        wx = X * Z; wy = Y * Z; wz = -(X * X + Y * Y);
    }

    // tetradtype 2: raise v first
    double fx = vx, fy = vy, fz = vz;
    vx = gupxx * fx + gupxy * fy + gupxz * fz;
    vy = gupxy * fx + gupyy * fy + gupyz * fz;
    vz = gupxz * fx + gupyz * fy + gupzz * fz;

    double f = 1.0 / chipn1;

    // normalize v
    fx = gxx * vx * vx + gyy * vy * vy + gzz * vz * vz +
         (gxy * vx * vy + gxz * vx * vz + gyz * vy * vz) * 2.0;
    fx = sqrt(fx * f);
    vx /= fx; vy /= fx; vz /= fx;

    // orthogonalize u w.r.t. v
    fx = gxx * vx * ux + gxy * vx * uy + gxz * vx * uz +
         gxy * vy * ux + gyy * vy * uy + gyz * vy * uz +
         gxz * vz * ux + gyz * vz * uy + gzz * vz * uz;
    fx *= f;
    ux -= fx * vx; uy -= fx * vy; uz -= fx * vz;

    // normalize u
    fx = gxx * ux * ux + gyy * uy * uy + gzz * uz * uz +
         (gxy * ux * uy + gxz * ux * uz + gyz * uy * uz) * 2.0;
    fx = sqrt(fx * f);
    ux /= fx; uy /= fx; uz /= fx;

    // orthogonalize w w.r.t. v
    fx = gxx * vx * wx + gxy * vx * wy + gxz * vx * wz +
         gxy * vy * wx + gyy * vy * wy + gyz * vy * wz +
         gxz * vz * wx + gyz * vz * wy + gzz * vz * wz;
    fx *= f;
    wx -= fx * vx; wy -= fx * vy; wz -= fx * vz;

    // orthogonalize w w.r.t. u
    fx = gxx * ux * wx + gxy * ux * wy + gxz * ux * wz +
         gxy * uy * wx + gyy * uy * wy + gyz * uy * wz +
         gxz * uz * wx + gyz * uz * wy + gzz * uz * wz;
    fx *= f;
    wx -= fx * ux; wy -= fx * uy; wz -= fx * uz;

    // normalize w
    fx = gxx * wx * wx + gyy * wy * wy + gzz * wz * wz +
         (gxy * wx * wy + gxz * wx * wz + gyz * wy * wz) * 2.0;
    fx = sqrt(fx * f);
    wx /= fx; wy /= fx; wz /= fx;

    // E_i (lower index with metric, raised by chipn1)
    double HEx = (gxx * Ex + gxy * Ey + gxz * Ez) * chipn1;
    double HEy = (gxy * Ex + gyy * Ey + gyz * Ez) * chipn1;
    double HEz = (gxz * Ex + gyz * Ey + gzz * Ez) * chipn1;

    // sqrt(gamma) * r x B
    f = sqrt(f);
    double f3 = f * f * f; // f^3

    // sqrt(gamma) * r x B
    double HBx = (vy * Bz - vz * By) * f3;
    double HBy = (vz * Bx - vx * Bz) * f3;
    double HBz = (vx * By - vy * Bx) * f3;

    // m = (phi - i theta)/sqrt(2) (tetradtype != 1)
    // real part: (HEx-HBx)*ux + (HEy-HBy)*uy + (HEz-HBz)*uz
    Rphi2 =  (HEx - HBx) * ux + (HEy - HBy) * uy + (HEz - HBz) * uz;
    // imaginary part: (HEx-HBx)*wx + (HEy-HBy)*wy + (HEz-HBz)*wz
    Iphi2 =  (HEx - HBx) * wx + (HEy - HBy) * wy + (HEz - HBz) * wz;

    Rphi2 *= 0.5;
    Iphi2 *= 0.5;
}

__device__ inline void getnpem1_point_device(
    double X, double Y, double Z, double chi,
    double gxx, double gxy, double gxz, double gyy, double gyz, double gzz,
    double Ex, double Ey, double Ez, double Bx, double By, double Bz,
    double &Rphi1, double &Iphi1
) {
#if (tetradtype == 0)
    d_getnpem1_point<0>(X, Y, Z, chi, gxx, gxy, gxz, gyy, gyz, gzz, Ex, Ey, Ez, Bx, By, Bz, Rphi1, Iphi1);
#elif (tetradtype == 1)
    d_getnpem1_point<1>(X, Y, Z, chi, gxx, gxy, gxz, gyy, gyz, gzz, Ex, Ey, Ez, Bx, By, Bz, Rphi1, Iphi1);
#else
    d_getnpem1_point<2>(X, Y, Z, chi, gxx, gxy, gxz, gyy, gyz, gzz, Ex, Ey, Ez, Bx, By, Bz, Rphi1, Iphi1);
#endif
}

#endif /* GETNPEM2_GPU_CUH */
