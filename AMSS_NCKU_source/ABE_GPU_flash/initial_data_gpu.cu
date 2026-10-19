#include "initial_data_gpu.h"

#include <math.h>

#include "gpu_runtime_coordinator.h"

namespace {

__device__ __forceinline__ double fallback_radius(const double* x, int nx) {
    if (nx > 1) return 0.5 * fabs(x[1] - x[0]);
    return 1.0e-14;
}

__device__ __forceinline__ double conformal_correction(
    int method,
    double mass,
    double px,
    double py,
    double pz,
    double sx,
    double sy,
    double sz,
    double nx,
    double ny,
    double nz,
    double rr) {
    const double one = 1.0;
    const double pp = sqrt(px * px + py * py + pz * pz);
    const double ss = sqrt(sx * sx + sy * sy + sz * sz);
    const double mup = (pp > 0.0) ? (px * nx + py * ny + pz * nz) / pp : 0.0;
    const double mus = (ss > 0.0) ? (sx * nx + sy * ny + sz * nz) / ss : 0.0;
    const double b = 2.0 * rr / mass;
    const double ell = 1.0 / (1.0 + b);
    const double ell2 = ell * ell;
    const double ell3 = ell2 * ell;
    const double ell4 = ell2 * ell2;
    const double ell5 = ell4 * ell;
    const double tmp = (py * sz - pz * sy) * nx +
                       (pz * sx - px * sz) * ny +
                       (px * sy - py * sx) * nz;

    const double u1 = 5.0 / 8.0 * ell *
                      (1.0 - 2.0 * ell + 2.0 * ell2 - ell3 + ell4 / 5.0);
    const double u2 = (15.0 + 117.0 * ell - 79.0 * ell2 + 43.0 * ell3 -
                       14.0 * ell4 + 2.0 * ell5 + 84.0 * log(ell) / b) /
                      (40.0 * b * b);

    if (method == AMSS_GPU_INIT_CAO) {
        const double u3 = ell / 20.0 * (1.0 + ell + ell2 - 4.0 * ell3 + 2.0 * ell4);
        const double u4 = ell2 / 10.0 * (10.0 - 25.0 * ell + 21.0 * ell2 - 6.0 * ell3);
        return pp * pp / (mass * mass) * (u1 + u2 * (3.0 * mup * mup - one)) +
               6.0 * u3 / (mass * mass * mass * mass) * ss * ss * (1.0 + mus * mus) +
               u4 / (mass * mass * mass) * tmp;
    }

    const double u3 = ell + ell2 + ell3 - 4.0 * ell4 + 2.0 * ell5;
    const double u4 = -b * b * ell5;
    const double u5 = b * (1.0 + 5.0 * b + 10.0 * b * b) * ell5 / 80.0;
    return pp * pp / (mass * mass) * (u1 + u2 * (3.0 * mup * mup - one)) +
           2.0 * ss * ss / (5.0 * mass * mass * mass * mass) *
               (u3 + u4 * (3.0 * mus * mus - one)) +
           u5 * tmp;
}

__device__ __forceinline__ void add_puncture_extrinsic_curvature(
    double px,
    double py,
    double pz,
    double sx,
    double sy,
    double sz,
    double nx,
    double ny,
    double nz,
    double rr,
    double& axx,
    double& axy,
    double& axz,
    double& ayy,
    double& ayz,
    double& azz) {
    const double one = 1.0;
    const double half = 0.5;
    const double three = 3.0;
    const double tmp = px * nx + py * ny + pz * nz;
    const double scale = three / (rr * rr);

    axx += (half * (px * nx + nx * px - (one - nx * nx) * tmp) +
            (nx * sy * nz - nx * sz * ny + nx * sy * nz - nx * sz * ny) / rr) * scale;
    ayy += (half * (py * ny + ny * py - (one - ny * ny) * tmp) +
            (ny * sz * nx - ny * sx * nz + ny * sz * nx - ny * sx * nz) / rr) * scale;
    azz += (half * (pz * nz + nz * pz - (one - nz * nz) * tmp) +
            (nz * sx * ny - nz * sy * nx + nz * sx * ny - nz * sy * nx) / rr) * scale;
    axy += (half * (px * ny + nx * py + nx * ny * tmp) +
            (nx * sz * nx - nx * sx * nz + ny * sy * nz - ny * sz * ny) / rr) * scale;
    axz += (half * (px * nz + nx * pz + nx * nz * tmp) +
            (nx * sx * ny - nx * sy * nx + nz * sy * nz - nz * sz * ny) / rr) * scale;
    ayz += (half * (py * nz + ny * pz + ny * nz * tmp) +
            (ny * sx * ny - ny * sy * nx + nz * sz * nx - nz * sx * nz) / rr) * scale;
}

__global__ void initial_data_puncture_kernel(
    int method,
    int nx_count,
    int ny_count,
    int nz_count,
    const double* x,
    const double* y,
    const double* z,
    double* chi,
    double* trK,
    double* gxx,
    double* gxy,
    double* gxz,
    double* gyy,
    double* gyz,
    double* gzz,
    double* Axx,
    double* Axy,
    double* Axz,
    double* Ayy,
    double* Ayz,
    double* Azz,
    double* Gmx,
    double* Gmy,
    double* Gmz,
    double* Lap,
    double* Sfx,
    double* Sfy,
    double* Sfz,
    double* dtSfx,
    double* dtSfy,
    double* dtSfz,
    const double* mass,
    const double* porg,
    const double* pmom,
    const double* spin,
    int bh_count) {
    const int point_count = nx_count * ny_count * nz_count;
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= point_count) return;

    const int i = index % nx_count;
    const int j = (index / nx_count) % ny_count;
    const int k = index / (nx_count * ny_count);
    const double xx = x[i];
    const double yy = y[j];
    const double zz = z[k];
    const double half = 0.5;
    const double tiny_rr = 1.0e-14;
    const double radius_floor = fallback_radius(x, nx_count);

    double psi = 1.0;
    double axx = 0.0;
    double axy = 0.0;
    double axz = 0.0;
    double ayy = 0.0;
    double ayz = 0.0;
    double azz = 0.0;

    for (int bh = 0; bh < bh_count; ++bh) {
        const int bh3 = 3 * bh;
        const double m = mass[bh];
        const double px = pmom[bh3 + 0];
        const double py = pmom[bh3 + 1];
        const double pz = pmom[bh3 + 2];
        const double sx = spin[bh3 + 0];
        const double sy = spin[bh3 + 1];
        const double sz = spin[bh3 + 2];
        double dx = xx - porg[bh3 + 0];
        double dy = yy - porg[bh3 + 1];
        double dz = zz - porg[bh3 + 2];
        double rr = sqrt(dx * dx + dy * dy + dz * dz);
        rr = (rr < tiny_rr) ? radius_floor : rr;
        const double inv_rr = 1.0 / rr;
        const double ex = dx * inv_rr;
        const double ey = dy * inv_rr;
        const double ez = dz * inv_rr;

        psi += conformal_correction(method, m, px, py, pz, sx, sy, sz, ex, ey, ez, rr) +
               half * m * inv_rr;
        add_puncture_extrinsic_curvature(px, py, pz, sx, sy, sz, ex, ey, ez, rr,
                                         axx, axy, axz, ayy, ayz, azz);
    }

    const double psi2 = psi * psi;
    const double psi4 = psi2 * psi2;
    const double psi6 = psi4 * psi2;
    chi[index] = 1.0 / psi4 - 1.0;
    Lap[index] = 1.0 / psi2 - 1.0;
    Axx[index] = axx / psi6;
    Axy[index] = axy / psi6;
    Axz[index] = axz / psi6;
    Ayy[index] = ayy / psi6;
    Ayz[index] = ayz / psi6;
    Azz[index] = azz / psi6;

    trK[index] = 0.0;
    gxx[index] = 0.0;
    gxy[index] = 0.0;
    gxz[index] = 0.0;
    gyy[index] = 0.0;
    gyz[index] = 0.0;
    gzz[index] = 0.0;
    Gmx[index] = 0.0;
    Gmy[index] = 0.0;
    Gmz[index] = 0.0;
    Sfx[index] = 0.0;
    Sfy[index] = 0.0;
    Sfz[index] = 0.0;
    dtSfx[index] = 0.0;
    dtSfy[index] = 0.0;
    dtSfz[index] = 0.0;
}

}  // namespace

void gpu_initial_data_puncture_launch(
    cudaStream_t stream,
    int method,
    const int shape[3],
    const double* x,
    const double* y,
    const double* z,
    double* chi,
    double* trK,
    double* gxx,
    double* gxy,
    double* gxz,
    double* gyy,
    double* gyz,
    double* gzz,
    double* Axx,
    double* Axy,
    double* Axz,
    double* Ayy,
    double* Ayz,
    double* Azz,
    double* Gmx,
    double* Gmy,
    double* Gmz,
    double* Lap,
    double* Sfx,
    double* Sfy,
    double* Sfz,
    double* dtSfx,
    double* dtSfy,
    double* dtSfz,
    const double* mass,
    const double* porg,
    const double* pmom,
    const double* spin,
    int bh_count) {
    const int point_count = shape[0] * shape[1] * shape[2];
    const int threads_per_block = 256;
    const int block_count = (point_count + threads_per_block - 1) / threads_per_block;
    initial_data_puncture_kernel<<<block_count, threads_per_block, 0, stream>>>(
        method, shape[0], shape[1], shape[2], x, y, z, chi, trK, gxx, gxy, gxz, gyy, gyz, gzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz, Gmx, Gmy, Gmz, Lap, Sfx, Sfy, Sfz, dtSfx, dtSfy, dtSfz,
        mass, porg, pmom, spin, bh_count);
    CUDA_CHECK(cudaGetLastError());
}
