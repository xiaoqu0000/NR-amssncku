#ifndef GETNP4_GPU_CUH
#define GETNP4_GPU_CUH

#include <cmath>

// Single-point Newman-Penrose Weyl scalar Psi4 computation.
// Ported from getnp4_point in getnp4.f90.
// Pre-computed derivatives are passed as input (no fderivs calls).
// Tetradtype: 0 = V directly (no raise), order V->U->W
//              1 = W first, order W->U->V
//              2 = V raised to contravariant, order V->U->W
template<int TETRAD>
__device__ __forceinline__ void d_getnp4_point(
    double X, double Y, double Z,
    double chi, double trK,
    double dxx, double gxy, double gxz, double dyy, double gyz, double dzz,
    double Axx, double Axy, double Axz, double Ayy, double Ayz, double Azz,
    double chix, double chiy, double chiz,
    double trKx, double trKy, double trKz,
    double& Axxx, double& Axxy, double& Axxz,
    double& Axyx, double& Axyy, double& Axyz,
    double& Axzx, double& Axzy, double& Axzz,
    double& Ayyx, double& Ayyy, double& Ayyz,
    double& Ayzx, double& Ayzy, double& Ayzz,
    double& Azzx, double& Azzy, double& Azzz,
    double Gamxxx, double Gamxxy, double Gamxxz,
    double Gamxyy, double Gamxyz, double Gamxzz,
    double Gamyxx, double Gamyxy, double Gamyxz,
    double Gamyyy, double Gamyyz, double Gamyzz,
    double Gamzxx, double Gamzxy, double Gamzxz,
    double Gamzyy, double Gamzyz, double Gamzzz,
    double Rxx, double Rxy, double Rxz,
    double Ryy, double Ryz, double Rzz,
    double& Rpsi4, double& Ipsi4
) {
    const double ZEO = 0.0;
    const double ONE = 1.0;
    const double TWO = 2.0;
    const double F1o3 = 1.0 / 3.0;
    const double TINYRR = 1e-14;

    double gxx = dxx + ONE;
    double gyy = dyy + ONE;
    double gzz = dzz + ONE;
    double chipn1 = chi + ONE;

    // Invert metric
    double gupzz = gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz
                 - gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz;
    double gupxx =  (gyy * gzz - gyz * gyz) / gupzz;
    double gupxy = -(gxy * gzz - gyz * gxz) / gupzz;
    double gupxz =  (gxy * gyz - gyy * gxz) / gupzz;
    double gupyy =  (gxx * gzz - gxz * gxz) / gupzz;
    double gupyz = -(gxx * gyz - gxy * gxz) / gupzz;
    double gupzz_ = (gxx * gyy - gxy * gxy) / gupzz;

    // Initialize tetrad vectors
    double vx, vy, vz, ux, uy, uz, wx, wy, wz;

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

    double f = ONE / chipn1;

    // Gram-Schmidt orthogonalization
    if constexpr (TETRAD == 2) {
        double fxt = vx, fyt = vy, fzt = vz;
        vx = gupxx * fxt + gupxy * fyt + gupxz * fzt;
        vy = gupxy * fxt + gupyy * fyt + gupyz * fzt;
        vz = gupxz * fxt + gupyz * fyt + gupzz_ * fzt;
    }

    if constexpr (TETRAD == 1) {
        // W first
        double temp;
        temp = gxx*wx*wx + gyy*wy*wy + gzz*wz*wz + (gxy*wx*wy + gxz*wx*wz + gyz*wy*wz)*TWO;
        temp = sqrt(temp*f);
        wx /= temp; wy /= temp; wz /= temp;

        temp = gxx*wx*ux + gxy*wx*uy + gxz*wx*uz +
               gxy*wy*ux + gyy*wy*uy + gyz*wy*uz +
               gxz*wz*ux + gyz*wz*uy + gzz*wz*uz;
        temp *= f;
        ux -= temp*wx; uy -= temp*wy; uz -= temp*wz;

        temp = gxx*ux*ux + gyy*uy*uy + gzz*uz*uz + (gxy*ux*uy + gxz*ux*uz + gyz*uy*uz)*TWO;
        temp = sqrt(temp*f);
        ux /= temp; uy /= temp; uz /= temp;

        temp = gxx*vx*wx + gxy*vx*wy + gxz*vx*wz +
               gxy*vy*wx + gyy*vy*wy + gyz*vy*wz +
               gxz*vz*wx + gyz*vz*wy + gzz*vz*wz;
        temp *= f;
        vx -= temp*wx; vy -= temp*wy; vz -= temp*wz;

        temp = gxx*ux*vx + gxy*ux*vy + gxz*ux*vz +
               gxy*uy*vx + gyy*uy*vy + gyz*uy*vz +
               gxz*uz*vx + gyz*uz*vy + gzz*uz*vz;
        temp *= f;
        vx -= temp*ux; vy -= temp*uy; vz -= temp*uz;

        temp = gxx*vx*vx + gyy*vy*vy + gzz*vz*vz + (gxy*vx*vy + gxz*vx*vz + gyz*vy*vz)*TWO;
        temp = sqrt(temp*f);
        vx /= temp; vy /= temp; vz /= temp;

    } else {
        // tetradtype == 0 or 2: V -> U -> W order
        double temp;
        temp = gxx*vx*vx + gyy*vy*vy + gzz*vz*vz + (gxy*vx*vy + gxz*vx*vz + gyz*vy*vz)*TWO;
        temp = sqrt(temp*f);
        vx /= temp; vy /= temp; vz /= temp;

        temp = gxx*vx*ux + gxy*vx*uy + gxz*vx*uz +
               gxy*vy*ux + gyy*vy*uy + gyz*vy*uz +
               gxz*vz*ux + gyz*vz*uy + gzz*vz*uz;
        temp *= f;
        ux -= temp*vx; uy -= temp*vy; uz -= temp*vz;

        temp = gxx*ux*ux + gyy*uy*uy + gzz*uz*uz + (gxy*ux*uy + gxz*ux*uz + gyz*uy*uz)*TWO;
        temp = sqrt(temp*f);
        ux /= temp; uy /= temp; uz /= temp;

        temp = gxx*vx*wx + gxy*vx*wy + gxz*vx*wz +
               gxy*vy*wx + gyy*vy*wy + gyz*vy*wz +
               gxz*vz*wx + gyz*vz*wy + gzz*vz*wz;
        temp *= f;
        wx -= temp*vx; wy -= temp*vy; wz -= temp*vz;

        temp = gxx*ux*wx + gxy*ux*wy + gxz*ux*wz +
               gxy*uy*wx + gyy*uy*wy + gyz*uy*wz +
               gxz*uz*wx + gyz*uz*wy + gzz*uz*wz;
        temp *= f;
        wx -= temp*ux; wy -= temp*uy; wz -= temp*uz;

        temp = gxx*wx*wx + gyy*wy*wy + gzz*wz*wz + (gxy*wx*wy + gxz*wx*wz + gyz*wy*wz)*TWO;
        temp = sqrt(temp*f);
        wx /= temp; wy /= temp; wz /= temp;
    }

    // Compute D_k K_ij
    Axxx = Axxx - (Gamxxx*Axx + Gamyxx*Axy + Gamzxx*Axz)*TWO - chix/chipn1*Axx + F1o3*gxx*trKx;
    Axxy = Axxy - (Gamxxy*Axx + Gamyxy*Axy + Gamzxy*Axz)*TWO - chiy/chipn1*Axx + F1o3*gxx*trKy;
    Axxz = Axxz - (Gamxxz*Axx + Gamyxz*Axy + Gamzxz*Axz)*TWO - chiz/chipn1*Axx + F1o3*gxx*trKz;
    Ayyx = Ayyx - (Gamxxy*Axy + Gamyxy*Ayy + Gamzxy*Ayz)*TWO - chix/chipn1*Ayy + F1o3*gyy*trKx;
    Ayyy = Ayyy - (Gamxyy*Axy + Gamyyy*Ayy + Gamzyy*Ayz)*TWO - chiy/chipn1*Ayy + F1o3*gyy*trKy;
    Ayyz = Ayyz - (Gamxyz*Axy + Gamyyz*Ayy + Gamzyz*Ayz)*TWO - chiz/chipn1*Ayy + F1o3*gyy*trKz;
    Azzx = Azzx - (Gamxxz*Axz + Gamyxz*Ayz + Gamzxz*Azz)*TWO - chix/chipn1*Azz + F1o3*gzz*trKx;
    Azzy = Azzy - (Gamxyz*Axz + Gamyyz*Ayz + Gamzyz*Azz)*TWO - chiy/chipn1*Azz + F1o3*gzz*trKy;
    Azzz = Azzz - (Gamxzz*Axz + Gamyzz*Ayz + Gamzzz*Azz)*TWO - chiz/chipn1*Azz + F1o3*gzz*trKz;
    Axyx = Axyx - (Gamxxy*Axx + Gamyxy*Axy + Gamzxy*Axz + Gamxxx*Axy + Gamyxx*Ayy + Gamzxx*Ayz) - chix/chipn1*Axy + F1o3*gxy*trKx;
    Axyy = Axyy - (Gamxyy*Axx + Gamyyy*Axy + Gamzyy*Axz + Gamxxy*Axy + Gamyxy*Ayy + Gamzxy*Ayz) - chiy/chipn1*Axy + F1o3*gxy*trKy;
    Axyz = Axyz - (Gamxyz*Axx + Gamyyz*Axy + Gamzyz*Axz + Gamxxz*Axy + Gamyxz*Ayy + Gamzxz*Ayz) - chiz/chipn1*Axy + F1o3*gxy*trKz;
    Axzx = Axzx - (Gamxxz*Axx + Gamyxz*Axy + Gamzxz*Axz + Gamxxx*Axz + Gamyxx*Ayz + Gamzxx*Azz) - chix/chipn1*Axz + F1o3*gxz*trKx;
    Axzy = Axzy - (Gamxyz*Axx + Gamyyz*Axy + Gamzyz*Axz + Gamxxy*Axz + Gamyxy*Ayz + Gamzxy*Azz) - chiy/chipn1*Axz + F1o3*gxz*trKy;
    Axzz = Axzz - (Gamxzz*Axx + Gamyzz*Axy + Gamzzz*Axz + Gamxxz*Axz + Gamyxz*Ayz + Gamzxz*Azz) - chiz/chipn1*Axz + F1o3*gxz*trKz;
    Ayzx = Ayzx - (Gamxxz*Axy + Gamyxz*Ayy + Gamzxz*Ayz + Gamxxy*Axz + Gamyxy*Ayz + Gamzxy*Azz) - chix/chipn1*Ayz + F1o3*gyz*trKx;
    Ayzy = Ayzy - (Gamxyz*Axy + Gamyyz*Ayy + Gamzyz*Ayz + Gamxyy*Axz + Gamyyy*Ayz + Gamzyy*Azz) - chiy/chipn1*Ayz + F1o3*gyz*trKy;
    Ayzz = Ayzz - (Gamxzz*Axy + Gamyzz*Ayy + Gamzzz*Ayz + Gamxyz*Axz + Gamyyz*Ayz + Gamzyz*Azz) - chiz/chipn1*Ayz + F1o3*gyz*trKz;

    // Symmetrize B_ij
    double Bxx = (vy*(Axxy - Axyx) + vz*(Axxz - Axzx)) * f;
    double Byy = (vx*(Ayyx - Axyy) + vz*(Ayyz - Ayzy)) * f;
    double Bzz = (vx*(Azzx - Axzz) + vy*(Azzy - Ayzz)) * f;
    double Bxy = (vx*(Axyx - (Axxy+Axyx)/TWO) + vy*(Axyy-Ayyx)/TWO + vz*(Axyz - (Axzy+Ayzx)/TWO)) * f;
    double Bxz = (vx*(Axzx - (Axxz+Axzx)/TWO) + vy*(Axzy - (Axyz+Ayzx)/TWO) + vz*(Axzz-Azzx)/TWO) * f;
    double Byz = (vx*(Ayzx - (Axyz+Axzy)/TWO) + vy*(Ayzy - (Ayyz+Ayzy)/TWO) + vz*(Ayzz-Azzy)/TWO) * f;

    // Physical K_ij = A_ij + (1/3)*trK*g_ij
    double Kxx_c = Axx + F1o3 * trK * gxx;
    double Kxy_c = Axy + F1o3 * trK * gxy;
    double Kxz_c = Axz + F1o3 * trK * gxz;
    double Kyy_c = Ayy + F1o3 * trK * gyy;
    double Kyz_c = Ayz + F1o3 * trK * gyz;
    double Kzz_c = Azz + F1o3 * trK * gzz;

    // Compute E_ij
    double Exx = gupxx*Kxx_c*Kxx_c + gupyy*Kxy_c*Kxy_c + gupzz_*Kxz_c*Kxz_c +
                 TWO*(gupxy*Kxx_c*Kxy_c + gupxz*Kxx_c*Kxz_c + gupyz*Kxy_c*Kxz_c);
    double Eyy = gupxx*Kxy_c*Kxy_c + gupyy*Kyy_c*Kyy_c + gupzz_*Kyz_c*Kyz_c +
                 TWO*(gupxy*Kxy_c*Kyy_c + gupxz*Kxy_c*Kyz_c + gupyz*Kyy_c*Kyz_c);
    double Ezz = gupxx*Kxz_c*Kxz_c + gupyy*Kyz_c*Kyz_c + gupzz_*Kzz_c*Kzz_c +
                 TWO*(gupxy*Kxz_c*Kyz_c + gupxz*Kxz_c*Kzz_c + gupyz*Kyz_c*Kzz_c);
    double Exy = gupxx*Kxx_c*Kxy_c + gupyy*Kxy_c*Kyy_c + gupzz_*Kxz_c*Kyz_c +
                 gupxy*(Kxx_c*Kyy_c + Kxy_c*Kxy_c) +
                 gupxz*(Kxx_c*Kyz_c + Kxz_c*Kxy_c) +
                 gupyz*(Kxy_c*Kyz_c + Kxz_c*Kyy_c);
    double Exz = gupxx*Kxx_c*Kxz_c + gupyy*Kxy_c*Kyz_c + gupzz_*Kxz_c*Kzz_c +
                 gupxy*(Kxx_c*Kyz_c + Kxy_c*Kxz_c) +
                 gupxz*(Kxx_c*Kzz_c + Kxz_c*Kxz_c) +
                 gupyz*(Kxy_c*Kzz_c + Kxz_c*Kyz_c);
    double Eyz = gupxx*Kxy_c*Kxz_c + gupyy*Kyy_c*Kyz_c + gupzz_*Kyz_c*Kzz_c +
                 gupxy*(Kxy_c*Kyz_c + Kyy_c*Kxz_c) +
                 gupxz*(Kxy_c*Kzz_c + Kyz_c*Kxz_c) +
                 gupyz*(Kyy_c*Kzz_c + Kyz_c*Kyz_c);

    Exx = Rxx - (Exx - Kxx_c*trK)*f - Bxx;
    Exy = Rxy - (Exy - Kxy_c*trK)*f - Bxy;
    Exz = Rxz - (Exz - Kxz_c*trK)*f - Bxz;
    Eyy = Ryy - (Eyy - Kyy_c*trK)*f - Byy;
    Eyz = Ryz - (Eyz - Kyz_c*trK)*f - Byz;
    Ezz = Rzz - (Ezz - Kzz_c*trK)*f - Bzz;

    // Projection tensors
    double uuwwxx = ux * ux - wx * wx;
    double uuwwxy = ux * uy - wx * wy;
    double uuwwxz = ux * uz - wx * wz;
    double uuwwyy = uy * uy - wy * wy;
    double uuwwyz = uy * uz - wy * wz;
    double uuwwzz = uz * uz - wz * wz;

    double uwxx = ux * wx + wx * ux;
    double uwxy = ux * wy + wx * uy;
    double uwxz = ux * wz + wx * uz;
    double uwyy = uy * wy + wy * uy;
    double uwyz = uy * wz + wy * uz;
    double uwzz = uz * wz + wz * uz;

    // Combine Psi4
    Rpsi4 = -(Exx*uuwwxx + Eyy*uuwwyy + Ezz*uuwwzz + (Exy*uuwwxy + Exz*uuwwxz + Eyz*uuwwyz)*TWO) / TWO;
    Ipsi4 = -(Exx*uwxx + Eyy*uwyy + Ezz*uwzz + (Exy*uwxy + Exz*uwxz + Eyz*uwyz)*TWO) / TWO;
}

#endif // GETNP4_GPU_CUH
