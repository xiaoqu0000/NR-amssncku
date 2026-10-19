#ifndef GETNPEM2_GPU_CUH
#define GETNPEM2_GPU_CUH

#include <math.h>

#include "macrodef.h"

namespace amss_gpu_em {

struct NpValue {
    double real;
    double imag;
};

__device__ __forceinline__ void em_inverse_metric(
    double gxx,
    double gxy,
    double gxz,
    double gyy,
    double gyz,
    double gzz,
    double& gupxx,
    double& gupxy,
    double& gupxz,
    double& gupyy,
    double& gupyz,
    double& gupzz
) {
    double determinant = gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz
        - gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz;
    gupxx = (gyy * gzz - gyz * gyz) / determinant;
    gupxy = -(gxy * gzz - gyz * gxz) / determinant;
    gupxz = (gxy * gyz - gyy * gxz) / determinant;
    gupyy = (gxx * gzz - gxz * gxz) / determinant;
    gupyz = -(gxx * gyz - gxy * gxz) / determinant;
    gupzz = (gxx * gyy - gxy * gxy) / determinant;
}

__device__ __forceinline__ double em_metric_dot(
    double ax,
    double ay,
    double az,
    double bx,
    double by,
    double bz,
    double gxx,
    double gxy,
    double gxz,
    double gyy,
    double gyz,
    double gzz
) {
    return gxx * ax * bx + gxy * ax * by + gxz * ax * bz
        + gxy * ay * bx + gyy * ay * by + gyz * ay * bz
        + gxz * az * bx + gyz * az * by + gzz * az * bz;
}

__device__ __forceinline__ void em_normalize_vector(
    double& ax,
    double& ay,
    double& az,
    double conformal_factor,
    double gxx,
    double gxy,
    double gxz,
    double gyy,
    double gyz,
    double gzz
) {
    const double norm = sqrt(em_metric_dot(ax, ay, az, ax, ay, az, gxx, gxy, gxz, gyy, gyz, gzz) * conformal_factor);
    ax /= norm;
    ay /= norm;
    az /= norm;
}

__device__ __forceinline__ void em_project_out(
    double basis_x,
    double basis_y,
    double basis_z,
    double conformal_factor,
    double gxx,
    double gxy,
    double gxz,
    double gyy,
    double gyz,
    double gzz,
    double& ax,
    double& ay,
    double& az
) {
    const double projection = em_metric_dot(
        basis_x,
        basis_y,
        basis_z,
        ax,
        ay,
        az,
        gxx,
        gxy,
        gxz,
        gyy,
        gyz,
        gzz
    ) * conformal_factor;
    ax -= projection * basis_x;
    ay -= projection * basis_y;
    az -= projection * basis_z;
}

__device__ __forceinline__ void em_initialize_tetrad(
    double x,
    double y,
    double z,
    double chi,
    double dxx,
    double gxy,
    double gxz,
    double dyy,
    double gyz,
    double dzz,
    double& gxx,
    double& gyy,
    double& gzz,
    double& chipn1,
    double& vx,
    double& vy,
    double& vz,
    double& ux,
    double& uy,
    double& uz,
    double& wx,
    double& wy,
    double& wz
) {
    constexpr double kZero = 0.0;
    constexpr double kOne = 1.0;
    constexpr double kTwo = 2.0;
    constexpr double kTinyRadius = 1.0e-14;

    gxx = dxx + kOne;
    gyy = dyy + kOne;
    gzz = dzz + kOne;
    chipn1 = chi + kOne;

    double gupxx = 0.0;
    double gupxy = 0.0;
    double gupxz = 0.0;
    double gupyy = 0.0;
    double gupyz = 0.0;
    double gupzz = 0.0;
    em_inverse_metric(gxx, gxy, gxz, gyy, gyz, gzz, gupxx, gupxy, gupxz, gupyy, gupyz, gupzz);

    if (fabs(x) < kTinyRadius && fabs(y) < kTinyRadius && fabs(z) < kTinyRadius) {
        vx = kTinyRadius;
        vy = kTinyRadius;
        vz = kTinyRadius;
    } else {
        vx = x;
        vy = y;
        vz = z;
    }

    if (fabs(x) < kTinyRadius && fabs(y) < kTinyRadius) {
        ux = -kTinyRadius;
        uy = kTinyRadius;
        uz = kZero;
        wx = kTinyRadius * z;
        wy = kTinyRadius * z;
        wz = -kTwo * kTinyRadius * kTinyRadius;
    } else {
        ux = -y;
        uy = x;
        uz = kZero;
        wx = x * z;
        wy = y * z;
        wz = -(x * x + y * y);
    }

    const double conformal_factor = kOne / chipn1;

#if (tetradtype == 0)
    em_normalize_vector(vx, vy, vz, conformal_factor, gxx, gxy, gxz, gyy, gyz, gzz);
    em_project_out(vx, vy, vz, conformal_factor, gxx, gxy, gxz, gyy, gyz, gzz, ux, uy, uz);
    em_normalize_vector(ux, uy, uz, conformal_factor, gxx, gxy, gxz, gyy, gyz, gzz);
    em_project_out(vx, vy, vz, conformal_factor, gxx, gxy, gxz, gyy, gyz, gzz, wx, wy, wz);
    em_project_out(ux, uy, uz, conformal_factor, gxx, gxy, gxz, gyy, gyz, gzz, wx, wy, wz);
    em_normalize_vector(wx, wy, wz, conformal_factor, gxx, gxy, gxz, gyy, gyz, gzz);
#elif (tetradtype == 1)
    em_normalize_vector(wx, wy, wz, conformal_factor, gxx, gxy, gxz, gyy, gyz, gzz);
    em_project_out(wx, wy, wz, conformal_factor, gxx, gxy, gxz, gyy, gyz, gzz, ux, uy, uz);
    em_normalize_vector(ux, uy, uz, conformal_factor, gxx, gxy, gxz, gyy, gyz, gzz);
    em_project_out(wx, wy, wz, conformal_factor, gxx, gxy, gxz, gyy, gyz, gzz, vx, vy, vz);
    em_project_out(ux, uy, uz, conformal_factor, gxx, gxy, gxz, gyy, gyz, gzz, vx, vy, vz);
    em_normalize_vector(vx, vy, vz, conformal_factor, gxx, gxy, gxz, gyy, gyz, gzz);
#elif (tetradtype == 2)
    const double raw_vx = vx;
    const double raw_vy = vy;
    const double raw_vz = vz;
    vx = gupxx * raw_vx + gupxy * raw_vy + gupxz * raw_vz;
    vy = gupxy * raw_vx + gupyy * raw_vy + gupyz * raw_vz;
    vz = gupxz * raw_vx + gupyz * raw_vy + gupzz * raw_vz;
    em_normalize_vector(vx, vy, vz, conformal_factor, gxx, gxy, gxz, gyy, gyz, gzz);
    em_project_out(vx, vy, vz, conformal_factor, gxx, gxy, gxz, gyy, gyz, gzz, ux, uy, uz);
    em_normalize_vector(ux, uy, uz, conformal_factor, gxx, gxy, gxz, gyy, gyz, gzz);
    em_project_out(vx, vy, vz, conformal_factor, gxx, gxy, gxz, gyy, gyz, gzz, wx, wy, wz);
    em_project_out(ux, uy, uz, conformal_factor, gxx, gxy, gxz, gyy, gyz, gzz, wx, wy, wz);
    em_normalize_vector(wx, wy, wz, conformal_factor, gxx, gxy, gxz, gyy, gyz, gzz);
#else
#error "not recognized tetradtype"
#endif
}

__device__ __forceinline__ NpValue compute_phi2_point(
    double x,
    double y,
    double z,
    double chi,
    double dxx,
    double gxy,
    double gxz,
    double dyy,
    double gyz,
    double dzz,
    double ex,
    double ey,
    double ez,
    double bx,
    double by,
    double bz
) {
    double gxx = 0.0;
    double gyy = 0.0;
    double gzz = 0.0;
    double chipn1 = 0.0;
    double vx = 0.0;
    double vy = 0.0;
    double vz = 0.0;
    double ux = 0.0;
    double uy = 0.0;
    double uz = 0.0;
    double wx = 0.0;
    double wy = 0.0;
    double wz = 0.0;

    em_initialize_tetrad(x, y, z, chi, dxx, gxy, gxz, dyy, gyz, dzz,
        gxx, gyy, gzz, chipn1, vx, vy, vz, ux, uy, uz, wx, wy, wz);

    const double h_ex = (gxx * ex + gxy * ey + gxz * ez) * chipn1;
    const double h_ey = (gxy * ex + gyy * ey + gyz * ez) * chipn1;
    const double h_ez = (gxz * ex + gyz * ey + gzz * ez) * chipn1;

    const double sqrt_chi_inverse = sqrt(1.0 / chipn1);
    const double volume_scale = sqrt_chi_inverse * sqrt_chi_inverse * sqrt_chi_inverse;

    const double h_bx = (vy * bz - vz * by) * volume_scale;
    const double h_by = (vz * bx - vx * bz) * volume_scale;
    const double h_bz = (vx * by - vy * bx) * volume_scale;

    const double field_x = h_ex - h_bx;
    const double field_y = h_ey - h_by;
    const double field_z = h_ez - h_bz;

    NpValue value;
#if (tetradtype == 1)
    value.real = 0.5 * (field_x * wx + field_y * wy + field_z * wz);
    value.imag = -0.5 * (field_x * ux + field_y * uy + field_z * uz);
#else
    value.real = 0.5 * (field_x * ux + field_y * uy + field_z * uz);
    value.imag = 0.5 * (field_x * wx + field_y * wy + field_z * wz);
#endif
    return value;
}

__device__ __forceinline__ NpValue compute_phi1_point(
    double x,
    double y,
    double z,
    double chi,
    double dxx,
    double gxy,
    double gxz,
    double dyy,
    double gyz,
    double dzz,
    double ex,
    double ey,
    double ez,
    double bx,
    double by,
    double bz
) {
    double gxx = 0.0;
    double gyy = 0.0;
    double gzz = 0.0;
    double chipn1 = 0.0;
    double vx = 0.0;
    double vy = 0.0;
    double vz = 0.0;
    double ux = 0.0;
    double uy = 0.0;
    double uz = 0.0;
    double wx = 0.0;
    double wy = 0.0;
    double wz = 0.0;

    em_initialize_tetrad(x, y, z, chi, dxx, gxy, gxz, dyy, gyz, dzz,
        gxx, gyy, gzz, chipn1, vx, vy, vz, ux, uy, uz, wx, wy, wz);

    const double h_ex = (gxx * ex + gxy * ey + gxz * ez) * chipn1;
    const double h_ey = (gxy * ex + gyy * ey + gyz * ez) * chipn1;
    const double h_ez = (gxz * ex + gyz * ey + gzz * ez) * chipn1;

    const double sqrt_chi_inverse = sqrt(1.0 / chipn1);
    const double volume_scale = sqrt_chi_inverse * sqrt_chi_inverse * sqrt_chi_inverse;

    const double h_bx = (uy * wz - uz * wy) * volume_scale;
    const double h_by = (uz * wx - ux * wz) * volume_scale;
    const double h_bz = (ux * wy - uy * wx) * volume_scale;

    NpValue value;
    value.real = 0.5 * (h_ex * vx + h_ey * vy + h_ez * vz);
    value.imag = 0.5 * (h_bx * bx + h_by * by + h_bz * bz);
    return value;
}

}  // namespace amss_gpu_em

#endif /* GETNPEM2_GPU_CUH */
