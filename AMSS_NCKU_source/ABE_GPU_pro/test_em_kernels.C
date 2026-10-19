// ===========================================================================
// test_em_kernels.C — Standalone accuracy test for EM GPU kernels
//
// Compares GPU kernel output against CPU Fortran equivalents.
// Must be linked against compiled .o files from the GPU build.
//
// IMPORTANT: GW must match the build's ghost_width (defined in macrodef.fh):
//   ghost_width = 2  (2nd-order)
//   ghost_width = 3  (4th-order)
//   ghost_width = 4  (6th-order, CURRENT DEFAULT)
//   ghost_width = 5  (8th-order)
//
// Compilation example:
//   mpicxx -DUSE_GPU -I/path/to/AMSS_NCKU_source \
//          -o test_em_kernels test_em_kernels.C \
//          empart_gpu.o getnpem2_gpu.o bssnEM_gpu_class.o \
//          empart.o getnpem2.o ... \
//          -L/usr/local/cuda/lib64 -lcudart -lstdc++ -lgfortran
// ===========================================================================

#include <cstdio>
#include <cmath>
#include <cstdlib>

#include <cuda_runtime.h>

#ifndef USE_GPU
#define USE_GPU
#endif

#include "empart.h"
#include "getnpem2.h"
#include "bssn_rhs.h"
#include "enforce_algebra.h"
#include "fmisc.h"

// ===========================================================================
// CUDA error-checking macro
// ===========================================================================
#define CUDA_CHECK(call) do {                                       \
    cudaError_t err = call;                                         \
    if (err != cudaSuccess) {                                       \
        fprintf(stderr, "CUDA error at %s:%d: %s\n", __FILE__, __LINE__, \
                cudaGetErrorString(err));                           \
        exit(1);                                                    \
    }                                                               \
} while (0)

// ===========================================================================
// Grid parameters — MUST match the build's ghost_width from macrodef.fh.
//   ghost_width = 2  (2nd-order)
//   ghost_width = 3  (4th-order)
//   ghost_width = 4  (6th-order) <-- current build default
//   ghost_width = 5  (8th-order)
//
// The ghost_width is generated from Finite_Diffenence_Method in AMSS_NCKU_Input.py
// and written to macrodef.fh by generate_macrodef.py at build time.
// Check the generated GW/macrodef.fh for the active value.
// ===========================================================================
static const int GW = 4;

// Tolerance for floating-point comparison (double precision)
static const double TOL = 1e-12;

// ===========================================================================
// fill_test_data — fill all grid arrays with smooth analytic functions
// (sin/cos products) over the full padded domain including ghost zones.
// ===========================================================================
static void fill_test_data(int ex[3],
                           double *X, double *Y, double *Z,
                           double *chi,
                           double *gxx0, double *gxy0, double *gxz0,
                           double *gyy0, double *gyz0, double *gzz0,
                           double *Lap0,
                           double *Sfx0, double *Sfy0, double *Sfz0,
                           double *trK0,
                           double *Ex0, double *Ey0, double *Ez0,
                           double *Bx0, double *By0, double *Bz0,
                           double *Kpsi0, double *Kphi0,
                           double *Jx, double *Jy, double *Jz, double *qchar)
{
    double dx = 0.1;
    for (int k = 0; k < ex[2]; k++) {
        for (int j = 0; j < ex[1]; j++) {
            for (int i = 0; i < ex[0]; i++) {
                int idx = i + j * ex[0] + k * ex[0] * ex[1];
                double x = (i - GW) * dx;
                double y = (j - GW) * dx;
                double z = (k - GW) * dx;

                X[idx] = x;
                Y[idx] = y;
                Z[idx] = z;

                chi[idx] = 1.0 + 0.01 * sin(x) * cos(y);

                gxx0[idx] = 0.0;  gxy0[idx] = 0.0;  gxz0[idx] = 0.0;
                gyy0[idx] = 0.0;  gyz0[idx] = 0.0;  gzz0[idx] = 0.0;

                Lap0[idx] = 1.0;

                Sfx0[idx] = 0.0;  Sfy0[idx] = 0.0;  Sfz0[idx] = 0.0;
                trK0[idx] = 0.0;

                Ex0[idx] = sin(y) * cos(z);
                Ey0[idx] = cos(x) * sin(z);
                Ez0[idx] = sin(x) * cos(y);

                Bx0[idx] = cos(y) * sin(z);
                By0[idx] = sin(x) * cos(z);
                Bz0[idx] = cos(x) * sin(y);

                Kpsi0[idx] = 0.0;
                Kphi0[idx] = 0.0;

                Jx[idx]   = 0.0;  Jy[idx]   = 0.0;
                Jz[idx]   = 0.0;  qchar[idx] = 0.0;
            }
        }
    }
}

// ===========================================================================
// compare_field — compare GPU vs CPU for a single field (interior only)
// Returns 1 if max absolute difference exceeds tolerance, 0 otherwise.
// ===========================================================================
static int compare_field(const char *name,
                         const double *gpu, const double *cpu,
                         int ex[3], double tol)
{
    double maxdiff = 0.0;
    int interior_count = 0;
    for (int k = GW; k < ex[2] - GW; k++) {
        for (int j = GW; j < ex[1] - GW; j++) {
            for (int i = GW; i < ex[0] - GW; i++) {
                int idx = i + j * ex[0] + k * ex[0] * ex[1];
                double d = fabs(gpu[idx] - cpu[idx]);
                if (d > maxdiff) maxdiff = d;
                interior_count++;
            }
        }
    }
    int fail = (maxdiff > tol) ? 1 : 0;
    printf("  max|D%-14s|= %#.2e  [%s]  (interior=%d)\n",
           name, maxdiff, fail ? "FAIL" : "PASS", interior_count);
    return fail;
}

// ===========================================================================
// Test 1: EM RHS (compute_rhs_empart)
// Compares gpu_compute_rhs_empart_launch vs f_compute_rhs_empart.
// ===========================================================================
static int test_em_rhs()
{
    printf("[1/4] EM RHS (compute_rhs_empart)\n");

    int ex[3]  = {8 + 2*GW, 8 + 2*GW, 8 + 2*GW};
    int pex[3] = {ex[0], ex[1], ex[2]};
    int n = ex[0] * ex[1] * ex[2];

    // ---- allocate host input arrays ----
    double *X    = new double[n];
    double *Y    = new double[n];
    double *Z    = new double[n];
    double *chi  = new double[n];
    double *gxx0 = new double[n];  double *gxy0 = new double[n];  double *gxz0 = new double[n];
    double *gyy0 = new double[n];  double *gyz0 = new double[n];  double *gzz0 = new double[n];
    double *Lap0 = new double[n];
    double *Sfx0 = new double[n];  double *Sfy0 = new double[n];  double *Sfz0 = new double[n];
    double *trK0 = new double[n];
    double *Ex0  = new double[n];  double *Ey0  = new double[n];  double *Ez0  = new double[n];
    double *Bx0  = new double[n];  double *By0  = new double[n];  double *Bz0  = new double[n];
    double *Kpsi0 = new double[n]; double *Kphi0 = new double[n];
    double *Jx   = new double[n];  double *Jy   = new double[n];
    double *Jz   = new double[n];  double *qchar = new double[n];

    fill_test_data(ex, X, Y, Z, chi,
                   gxx0, gxy0, gxz0, gyy0, gyz0, gzz0,
                   Lap0, Sfx0, Sfy0, Sfz0, trK0,
                   Ex0, Ey0, Ez0, Bx0, By0, Bz0, Kpsi0, Kphi0,
                   Jx, Jy, Jz, qchar);

    // ---- allocate GPU memory ----
    size_t nb = n * sizeof(double);
    double *d_X, *d_Y, *d_Z;
    double *d_chi, *d_Jx, *d_Jy, *d_Jz, *d_qchar;
    double *d_gxx0, *d_gxy0, *d_gxz0, *d_gyy0, *d_gyz0, *d_gzz0;
    double *d_Lap0, *d_Sfx0, *d_Sfy0, *d_Sfz0, *d_trK0;
    double *d_Ex0, *d_Ey0, *d_Ez0, *d_Bx0, *d_By0, *d_Bz0;
    double *d_Kpsi0, *d_Kphi0;
    double *d_Ex_rhs, *d_Ey_rhs, *d_Ez_rhs;
    double *d_Bx_rhs, *d_By_rhs, *d_Bz_rhs;
    double *d_Kpsi_rhs, *d_Kphi_rhs;
    double *d_rho, *d_Sx, *d_Sy, *d_Sz;
    double *d_Sxx, *d_Sxy, *d_Sxz, *d_Syy, *d_Syz, *d_Szz;

    CUDA_CHECK(cudaMalloc(&d_X,    nb));
    CUDA_CHECK(cudaMalloc(&d_Y,    nb));
    CUDA_CHECK(cudaMalloc(&d_Z,    nb));
    CUDA_CHECK(cudaMalloc(&d_chi,  nb));
    CUDA_CHECK(cudaMalloc(&d_Jx,   nb));
    CUDA_CHECK(cudaMalloc(&d_Jy,   nb));
    CUDA_CHECK(cudaMalloc(&d_Jz,   nb));
    CUDA_CHECK(cudaMalloc(&d_qchar, nb));
    CUDA_CHECK(cudaMalloc(&d_gxx0, nb));
    CUDA_CHECK(cudaMalloc(&d_gxy0, nb));
    CUDA_CHECK(cudaMalloc(&d_gxz0, nb));
    CUDA_CHECK(cudaMalloc(&d_gyy0, nb));
    CUDA_CHECK(cudaMalloc(&d_gyz0, nb));
    CUDA_CHECK(cudaMalloc(&d_gzz0, nb));
    CUDA_CHECK(cudaMalloc(&d_Lap0, nb));
    CUDA_CHECK(cudaMalloc(&d_Sfx0, nb));
    CUDA_CHECK(cudaMalloc(&d_Sfy0, nb));
    CUDA_CHECK(cudaMalloc(&d_Sfz0, nb));
    CUDA_CHECK(cudaMalloc(&d_trK0, nb));
    CUDA_CHECK(cudaMalloc(&d_Ex0,  nb));
    CUDA_CHECK(cudaMalloc(&d_Ey0,  nb));
    CUDA_CHECK(cudaMalloc(&d_Ez0,  nb));
    CUDA_CHECK(cudaMalloc(&d_Bx0,  nb));
    CUDA_CHECK(cudaMalloc(&d_By0,  nb));
    CUDA_CHECK(cudaMalloc(&d_Bz0,  nb));
    CUDA_CHECK(cudaMalloc(&d_Kpsi0, nb));
    CUDA_CHECK(cudaMalloc(&d_Kphi0, nb));
    CUDA_CHECK(cudaMalloc(&d_Ex_rhs,   nb));
    CUDA_CHECK(cudaMalloc(&d_Ey_rhs,   nb));
    CUDA_CHECK(cudaMalloc(&d_Ez_rhs,   nb));
    CUDA_CHECK(cudaMalloc(&d_Bx_rhs,   nb));
    CUDA_CHECK(cudaMalloc(&d_By_rhs,   nb));
    CUDA_CHECK(cudaMalloc(&d_Bz_rhs,   nb));
    CUDA_CHECK(cudaMalloc(&d_Kpsi_rhs, nb));
    CUDA_CHECK(cudaMalloc(&d_Kphi_rhs, nb));
    CUDA_CHECK(cudaMalloc(&d_rho,  nb));
    CUDA_CHECK(cudaMalloc(&d_Sx,   nb));
    CUDA_CHECK(cudaMalloc(&d_Sy,   nb));
    CUDA_CHECK(cudaMalloc(&d_Sz,   nb));
    CUDA_CHECK(cudaMalloc(&d_Sxx,  nb));
    CUDA_CHECK(cudaMalloc(&d_Sxy,  nb));
    CUDA_CHECK(cudaMalloc(&d_Sxz,  nb));
    CUDA_CHECK(cudaMalloc(&d_Syy,  nb));
    CUDA_CHECK(cudaMalloc(&d_Syz,  nb));
    CUDA_CHECK(cudaMalloc(&d_Szz,  nb));

    // ---- copy inputs to GPU ----
    cudaStream_t stream = 0;
    CUDA_CHECK(cudaMemcpy(d_X,    X,    nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Y,    Y,    nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Z,    Z,    nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_chi,  chi,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Jx,   Jx,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Jy,   Jy,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Jz,   Jz,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_qchar, qchar, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gxx0, gxx0, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gxy0, gxy0, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gxz0, gxz0, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gyy0, gyy0, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gyz0, gyz0, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gzz0, gzz0, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Lap0, Lap0, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Sfx0, Sfx0, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Sfy0, Sfy0, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Sfz0, Sfz0, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_trK0, trK0, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Ex0,  Ex0,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Ey0,  Ey0,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Ez0,  Ez0,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Bx0,  Bx0,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_By0,  By0,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Bz0,  Bz0,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Kpsi0, Kpsi0, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Kphi0, Kphi0, nb, cudaMemcpyHostToDevice));

    // ---- run GPU kernel ----
    // GPU parameter order (from empart.h):
    //   stream, ex, pex, X,Y,Z, chi, Jx,Jy,Jz,qchar,
    //   gxx,gxy,gxz,gyy,gyz,gzz, Lap, Sfx,Sfy,Sfz,trK,
    //   Ex,Ey,Ez,Bx,By,Bz, Kpsi,Kphi,
    //   Ex_rhs,Ey_rhs,Ez_rhs, Bx_rhs,By_rhs,Bz_rhs, Kpsi_rhs,Kphi_rhs,
    //   rho,Sx,Sy,Sz, Sxx,Sxy,Sxz,Syy,Syz,Szz,
    //   symmetry, lev, eps
    gpu_compute_rhs_empart_launch(
        stream, ex, pex,
        d_X, d_Y, d_Z,
        d_chi, d_Jx, d_Jy, d_Jz, d_qchar,
        d_gxx0, d_gxy0, d_gxz0, d_gyy0, d_gyz0, d_gzz0,
        d_Lap0, d_Sfx0, d_Sfy0, d_Sfz0, d_trK0,
        d_Ex0, d_Ey0, d_Ez0, d_Bx0, d_By0, d_Bz0,
        d_Kpsi0, d_Kphi0,
        d_Ex_rhs, d_Ey_rhs, d_Ez_rhs,
        d_Bx_rhs, d_By_rhs, d_Bz_rhs,
        d_Kpsi_rhs, d_Kphi_rhs,
        d_rho, d_Sx, d_Sy, d_Sz,
        d_Sxx, d_Sxy, d_Sxz,
        d_Syy, d_Syz, d_Szz,
        0 /*symmetry*/, 0 /*lev*/, 0.0 /*eps*/);
    CUDA_CHECK(cudaDeviceSynchronize());

    // ---- copy GPU outputs back to host ----
    double *gpu_Ex_rhs   = new double[n];
    double *gpu_Ey_rhs   = new double[n];
    double *gpu_Ez_rhs   = new double[n];
    double *gpu_Bx_rhs   = new double[n];
    double *gpu_By_rhs   = new double[n];
    double *gpu_Bz_rhs   = new double[n];
    double *gpu_Kpsi_rhs = new double[n];
    double *gpu_Kphi_rhs = new double[n];
    double *gpu_rho = new double[n];
    double *gpu_Sx  = new double[n];
    double *gpu_Sy  = new double[n];
    double *gpu_Sz  = new double[n];
    double *gpu_Sxx = new double[n];
    double *gpu_Sxy = new double[n];
    double *gpu_Sxz = new double[n];
    double *gpu_Syy = new double[n];
    double *gpu_Syz = new double[n];
    double *gpu_Szz = new double[n];

    CUDA_CHECK(cudaMemcpy(gpu_Ex_rhs,   d_Ex_rhs,   nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Ey_rhs,   d_Ey_rhs,   nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Ez_rhs,   d_Ez_rhs,   nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Bx_rhs,   d_Bx_rhs,   nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_By_rhs,   d_By_rhs,   nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Bz_rhs,   d_Bz_rhs,   nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Kpsi_rhs, d_Kpsi_rhs, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Kphi_rhs, d_Kphi_rhs, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_rho, d_rho, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Sx,  d_Sx,  nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Sy,  d_Sy,  nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Sz,  d_Sz,  nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Sxx, d_Sxx, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Sxy, d_Sxy, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Sxz, d_Sxz, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Syy, d_Syy, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Syz, d_Syz, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Szz, d_Szz, nb, cudaMemcpyDeviceToHost));

    // ---- run CPU Fortran function ----
    // CPU parameter order (from empart.h / bssnEM_class.C:857-876):
    //   ex, X,Y,Z, chi, gxx,gxy,gxz,gyy,gyz,gzz,
    //   Lap, Sfx,Sfy,Sfz,trK,
    //   Ex,Ey,Ez,Bx,By,Bz, Kpsi,Kphi, Jx,Jy,Jz,qchar,
    //   Ex_rhs,Ey_rhs,Ez_rhs, Bx_rhs,By_rhs,Bz_rhs, Kpsi_rhs,Kphi_rhs,
    //   rho,Sx,Sy,Sz, Sxx,Sxy,Sxz,Syy,Syz,Szz,
    //   symmetry, lev, ndeps (references)
    double *cpu_Ex_rhs   = new double[n];
    double *cpu_Ey_rhs   = new double[n];
    double *cpu_Ez_rhs   = new double[n];
    double *cpu_Bx_rhs   = new double[n];
    double *cpu_By_rhs   = new double[n];
    double *cpu_Bz_rhs   = new double[n];
    double *cpu_Kpsi_rhs = new double[n];
    double *cpu_Kphi_rhs = new double[n];
    double *cpu_rho = new double[n];
    double *cpu_Sx  = new double[n];
    double *cpu_Sy  = new double[n];
    double *cpu_Sz  = new double[n];
    double *cpu_Sxx = new double[n];
    double *cpu_Sxy = new double[n];
    double *cpu_Sxz = new double[n];
    double *cpu_Syy = new double[n];
    double *cpu_Syz = new double[n];
    double *cpu_Szz = new double[n];

    int    sym_ref = 0;
    int    lev_ref = 0;
    double ndeps_ref = 0.0;

    int ret = f_compute_rhs_empart(
        ex, X, Y, Z,
        chi,
        gxx0, gxy0, gxz0, gyy0, gyz0, gzz0,
        Lap0,
        Sfx0, Sfy0, Sfz0,
        trK0,
        Ex0, Ey0, Ez0, Bx0, By0, Bz0, Kpsi0, Kphi0,
        Jx, Jy, Jz, qchar,
        cpu_Ex_rhs, cpu_Ey_rhs, cpu_Ez_rhs,
        cpu_Bx_rhs, cpu_By_rhs, cpu_Bz_rhs,
        cpu_Kpsi_rhs, cpu_Kphi_rhs,
        cpu_rho, cpu_Sx, cpu_Sy, cpu_Sz,
        cpu_Sxx, cpu_Sxy, cpu_Sxz,
        cpu_Syy, cpu_Syz, cpu_Szz,
        sym_ref, lev_ref, ndeps_ref);

    printf("  CPU f_compute_rhs_empart ret=%d\n", ret);

    // ---- compare (interior only) ----
    int fails = 0;
    fails += compare_field("Ex_rhs",   gpu_Ex_rhs,   cpu_Ex_rhs,   ex, TOL);
    fails += compare_field("Ey_rhs",   gpu_Ey_rhs,   cpu_Ey_rhs,   ex, TOL);
    fails += compare_field("Ez_rhs",   gpu_Ez_rhs,   cpu_Ez_rhs,   ex, TOL);
    fails += compare_field("Bx_rhs",   gpu_Bx_rhs,   cpu_Bx_rhs,   ex, TOL);
    fails += compare_field("By_rhs",   gpu_By_rhs,   cpu_By_rhs,   ex, TOL);
    fails += compare_field("Bz_rhs",   gpu_Bz_rhs,   cpu_Bz_rhs,   ex, TOL);
    fails += compare_field("Kpsi_rhs", gpu_Kpsi_rhs, cpu_Kpsi_rhs, ex, TOL);
    fails += compare_field("Kphi_rhs", gpu_Kphi_rhs, cpu_Kphi_rhs, ex, TOL);
    fails += compare_field("rho",      gpu_rho,      cpu_rho,      ex, TOL);
    fails += compare_field("Sx",       gpu_Sx,       cpu_Sx,       ex, TOL);
    fails += compare_field("Sy",       gpu_Sy,       cpu_Sy,       ex, TOL);
    fails += compare_field("Sz",       gpu_Sz,       cpu_Sz,       ex, TOL);
    fails += compare_field("Sxx",      gpu_Sxx,      cpu_Sxx,      ex, TOL);
    fails += compare_field("Sxy",      gpu_Sxy,      cpu_Sxy,      ex, TOL);
    fails += compare_field("Sxz",      gpu_Sxz,      cpu_Sxz,      ex, TOL);
    fails += compare_field("Syy",      gpu_Syy,      cpu_Syy,      ex, TOL);
    fails += compare_field("Syz",      gpu_Syz,      cpu_Syz,      ex, TOL);
    fails += compare_field("Szz",      gpu_Szz,      cpu_Szz,      ex, TOL);

    // ---- cleanup ----
    delete[] X;    delete[] Y;    delete[] Z;
    delete[] chi;
    delete[] gxx0; delete[] gxy0; delete[] gxz0;
    delete[] gyy0; delete[] gyz0; delete[] gzz0;
    delete[] Lap0;
    delete[] Sfx0; delete[] Sfy0; delete[] Sfz0; delete[] trK0;
    delete[] Ex0;  delete[] Ey0;  delete[] Ez0;
    delete[] Bx0;  delete[] By0;  delete[] Bz0;
    delete[] Kpsi0; delete[] Kphi0;
    delete[] Jx;   delete[] Jy;   delete[] Jz;   delete[] qchar;

    delete[] gpu_Ex_rhs;   delete[] gpu_Ey_rhs;   delete[] gpu_Ez_rhs;
    delete[] gpu_Bx_rhs;   delete[] gpu_By_rhs;   delete[] gpu_Bz_rhs;
    delete[] gpu_Kpsi_rhs; delete[] gpu_Kphi_rhs;
    delete[] gpu_rho;
    delete[] gpu_Sx;  delete[] gpu_Sy;  delete[] gpu_Sz;
    delete[] gpu_Sxx; delete[] gpu_Sxy; delete[] gpu_Sxz;
    delete[] gpu_Syy; delete[] gpu_Syz; delete[] gpu_Szz;

    delete[] cpu_Ex_rhs;   delete[] cpu_Ey_rhs;   delete[] cpu_Ez_rhs;
    delete[] cpu_Bx_rhs;   delete[] cpu_By_rhs;   delete[] cpu_Bz_rhs;
    delete[] cpu_Kpsi_rhs; delete[] cpu_Kphi_rhs;
    delete[] cpu_rho;
    delete[] cpu_Sx;  delete[] cpu_Sy;  delete[] cpu_Sz;
    delete[] cpu_Sxx; delete[] cpu_Sxy; delete[] cpu_Sxz;
    delete[] cpu_Syy; delete[] cpu_Syz; delete[] cpu_Szz;

    CUDA_CHECK(cudaFree(d_X));    CUDA_CHECK(cudaFree(d_Y));    CUDA_CHECK(cudaFree(d_Z));
    CUDA_CHECK(cudaFree(d_chi));  CUDA_CHECK(cudaFree(d_Jx));   CUDA_CHECK(cudaFree(d_Jy));
    CUDA_CHECK(cudaFree(d_Jz));   CUDA_CHECK(cudaFree(d_qchar));
    CUDA_CHECK(cudaFree(d_gxx0)); CUDA_CHECK(cudaFree(d_gxy0)); CUDA_CHECK(cudaFree(d_gxz0));
    CUDA_CHECK(cudaFree(d_gyy0)); CUDA_CHECK(cudaFree(d_gyz0)); CUDA_CHECK(cudaFree(d_gzz0));
    CUDA_CHECK(cudaFree(d_Lap0));
    CUDA_CHECK(cudaFree(d_Sfx0)); CUDA_CHECK(cudaFree(d_Sfy0)); CUDA_CHECK(cudaFree(d_Sfz0));
    CUDA_CHECK(cudaFree(d_trK0));
    CUDA_CHECK(cudaFree(d_Ex0));  CUDA_CHECK(cudaFree(d_Ey0));  CUDA_CHECK(cudaFree(d_Ez0));
    CUDA_CHECK(cudaFree(d_Bx0));  CUDA_CHECK(cudaFree(d_By0));  CUDA_CHECK(cudaFree(d_Bz0));
    CUDA_CHECK(cudaFree(d_Kpsi0)); CUDA_CHECK(cudaFree(d_Kphi0));
    CUDA_CHECK(cudaFree(d_Ex_rhs));   CUDA_CHECK(cudaFree(d_Ey_rhs));
    CUDA_CHECK(cudaFree(d_Ez_rhs));   CUDA_CHECK(cudaFree(d_Bx_rhs));
    CUDA_CHECK(cudaFree(d_By_rhs));   CUDA_CHECK(cudaFree(d_Bz_rhs));
    CUDA_CHECK(cudaFree(d_Kpsi_rhs)); CUDA_CHECK(cudaFree(d_Kphi_rhs));
    CUDA_CHECK(cudaFree(d_rho));
    CUDA_CHECK(cudaFree(d_Sx));   CUDA_CHECK(cudaFree(d_Sy));   CUDA_CHECK(cudaFree(d_Sz));
    CUDA_CHECK(cudaFree(d_Sxx));  CUDA_CHECK(cudaFree(d_Sxy));  CUDA_CHECK(cudaFree(d_Sxz));
    CUDA_CHECK(cudaFree(d_Syy));  CUDA_CHECK(cudaFree(d_Syz));  CUDA_CHECK(cudaFree(d_Szz));

    return fails;
}

// ===========================================================================
// Test 2: EM Symmetry Boundary
// The GPU symmetry BD kernel applies parity reflection to ghost zones.
// There is no CPU Fortran equivalent for this operation, so we verify the
// parity rules directly by comparing ghost-zone values against reflected
// interior values.
//
// Parity rules (from empart_gpu.cu:1200-1202):
//   Ex:   ANTI in x, SYM  in y, SYM  in z
//   Ey:   SYM  in x, ANTI in y, SYM  in z
//   Ez:   SYM  in x, SYM  in y, ANTI in z
//   Bx:   SYM  in x, ANTI in y, ANTI in z
//   By:   ANTI in x, SYM  in y, ANTI in z
//   Bz:   ANTI in x, ANTI in y, SYM  in z
//   Kpsi: SYM  in x, SYM  in y, SYM  in z
//   Kphi: SYM  in x, SYM  in y, SYM  in z
//
// Reflection formula (without Vertex): ri = -i-1 for i<0
// ===========================================================================
static int test_em_symmetry_bd()
{
    printf("[2/4] EM Symmetry BD (apply_symmetry_bd_empart)\n");

    int ex[3]  = {8 + 2*GW, 8 + 2*GW, 8 + 2*GW};
    int pex[3] = {ex[0], ex[1], ex[2]};
    int n = ex[0] * ex[1] * ex[2];

    // ---- allocate and fill host arrays ----
    double *Ex   = new double[n];
    double *Ey   = new double[n];
    double *Ez   = new double[n];
    double *Bx   = new double[n];
    double *By   = new double[n];
    double *Bz   = new double[n];
    double *Kpsi = new double[n];
    double *Kphi = new double[n];

    double dx = 0.1;
    for (int k = 0; k < ex[2]; k++) {
        for (int j = 0; j < ex[1]; j++) {
            for (int i = 0; i < ex[0]; i++) {
                int idx = i + j * ex[0] + k * ex[0] * ex[1];
                double x = (i - GW) * dx;
                double y = (j - GW) * dx;
                double z = (k - GW) * dx;

                Ex[idx]   = sin(y) * cos(z) + 1.0;
                Ey[idx]   = cos(x) * sin(z) + 1.0;
                Ez[idx]   = sin(x) * cos(y) + 1.0;
                Bx[idx]   = cos(y) * sin(z) + 2.0;
                By[idx]   = sin(x) * cos(z) + 2.0;
                Bz[idx]   = cos(x) * sin(y) + 2.0;
                Kpsi[idx] = sin(x) * cos(y) + 3.0;
                Kphi[idx] = cos(x) * sin(y) + 3.0;
            }
        }
    }

    // ---- allocate GPU memory and copy ----
    size_t nb = n * sizeof(double);
    double *d_Ex, *d_Ey, *d_Ez, *d_Bx, *d_By, *d_Bz, *d_Kpsi, *d_Kphi;
    CUDA_CHECK(cudaMalloc(&d_Ex,   nb)); CUDA_CHECK(cudaMalloc(&d_Ey,   nb));
    CUDA_CHECK(cudaMalloc(&d_Ez,   nb)); CUDA_CHECK(cudaMalloc(&d_Bx,   nb));
    CUDA_CHECK(cudaMalloc(&d_By,   nb)); CUDA_CHECK(cudaMalloc(&d_Bz,   nb));
    CUDA_CHECK(cudaMalloc(&d_Kpsi, nb)); CUDA_CHECK(cudaMalloc(&d_Kphi, nb));

    cudaStream_t stream = 0;
    CUDA_CHECK(cudaMemcpy(d_Ex,   Ex,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Ey,   Ey,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Ez,   Ez,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Bx,   Bx,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_By,   By,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Bz,   Bz,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Kpsi, Kpsi, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Kphi, Kphi, nb, cudaMemcpyHostToDevice));

    // ---- run GPU symmetry BD ----
    gpu_apply_symmetry_bd_empart_launch(
        stream,
        ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
        d_Ex, d_Ey, d_Ez,
        d_Bx, d_By, d_Bz,
        d_Kpsi, d_Kphi);
    CUDA_CHECK(cudaDeviceSynchronize());

    // ---- copy results back ----
    CUDA_CHECK(cudaMemcpy(Ex,   d_Ex,   nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(Ey,   d_Ey,   nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(Ez,   d_Ez,   nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(Bx,   d_Bx,   nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(By,   d_By,   nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(Bz,   d_Bz,   nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(Kpsi, d_Kpsi, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(Kphi, d_Kphi, nb, cudaMemcpyDeviceToHost));

    // ---- verify parity rules for ghost zones ----
    // Check only the lower-z face (k < GW in padded index) as a representative
    // sample. The kernel applies symmetry to all 6 faces; checking one face
    // provides a strong correctness signal.
    //
    // For each ghost point at logical k<0 (padded k < GW), the reflected
    // source is at logical rk = -k-1 (padded = GW + (-k-1) = GW - k - 1).
    // Skipping the Vertex case (standard in this codebase).
    //
    // The parity factor depends on which dimensions are negative:
    //   fx = (i < GW) ? -1.0 : 1.0,  etc.

    struct {
        const char *name;
        double *data;
        // parity multipliers: x_anti, y_anti, z_anti
        // true = anti-symmetric (flip sign), false = symmetric
        bool x_anti, y_anti, z_anti;
    } fields[] = {
        {"Ex",   Ex,   true,  false, false},
        {"Ey",   Ey,   false, true,  false},
        {"Ez",   Ez,   false, false, true },
        {"Bx",   Bx,   false, true,  true },
        {"By",   By,   true,  false, true },
        {"Bz",   Bz,   true,  true,  false},
        {"Kpsi", Kpsi, false, false, false},
        {"Kphi", Kphi, false, false, false},
    };
    const int nfields = sizeof(fields) / sizeof(fields[0]);

    int fails = 0;
    for (int f = 0; f < nfields; f++) {
        double maxerr = 0.0;
        int checked = 0;

        // Check ghost zone at lower z-face: k < GW, i,j in [GW, ex-GW)
        for (int k = 0; k < GW; k++) {
            int rk_raw = GW - k - 1;  // reflected raw index (Cell: ri=-i-1): padded_source = GW + rk_raw
            for (int j = GW; j < ex[1] - GW; j++) {
                for (int i = GW; i < ex[0] - GW; i++) {
                    int gidx = i + j * ex[0] + k * ex[0] * ex[1];
                    int sidx = i + j * ex[0] + (rk_raw + GW) * ex[0] * ex[1];

                    // Only k is negative (k < GW), so only z_anti matters
                    double factor = fields[f].z_anti ? -1.0 : 1.0;
                    double expected = fields[f].data[sidx] * factor;
                    double diff = fabs(fields[f].data[gidx] - expected);
                    if (diff > maxerr) maxerr = diff;
                    checked++;
                }
            }
        }

        printf("  max|D%-14s|= %#.2e  [%s]  (checked=%d, z-face ghost)\n",
               fields[f].name, maxerr, (maxerr > TOL) ? "FAIL" : "PASS", checked);
        if (maxerr > TOL) fails++;
    }

    // ---- cleanup ----
    delete[] Ex; delete[] Ey; delete[] Ez;
    delete[] Bx; delete[] By; delete[] Bz;
    delete[] Kpsi; delete[] Kphi;
    CUDA_CHECK(cudaFree(d_Ex));   CUDA_CHECK(cudaFree(d_Ey));
    CUDA_CHECK(cudaFree(d_Ez));   CUDA_CHECK(cudaFree(d_Bx));
    CUDA_CHECK(cudaFree(d_By));   CUDA_CHECK(cudaFree(d_Bz));
    CUDA_CHECK(cudaFree(d_Kpsi)); CUDA_CHECK(cudaFree(d_Kphi));

    return fails;
}

// ===========================================================================
// fill_np_test_data — fill arrays needed for NP scalar tests.
// Uses the same analytic functions as fill_test_data, omitting fields
// not needed for NP scalars.
// ===========================================================================
static void fill_np_test_data(int ex[3],
                              double *X, double *Y, double *Z,
                              double *chi,
                              double *gxx0, double *gxy0, double *gxz0,
                              double *gyy0, double *gyz0, double *gzz0,
                              double *Ex0, double *Ey0, double *Ez0,
                              double *Bx0, double *By0, double *Bz0)
{
    double dx = 0.1;
    for (int k = 0; k < ex[2]; k++) {
        for (int j = 0; j < ex[1]; j++) {
            for (int i = 0; i < ex[0]; i++) {
                int idx = i + j * ex[0] + k * ex[0] * ex[1];
                double x = (i - GW) * dx;
                double y = (j - GW) * dx;
                double z = (k - GW) * dx;

                X[idx] = x;
                Y[idx] = y;
                Z[idx] = z;

                chi[idx] = 1.0 + 0.01 * sin(x) * cos(y);

                gxx0[idx] = 0.0;  gxy0[idx] = 0.0;  gxz0[idx] = 0.0;
                gyy0[idx] = 0.0;  gyz0[idx] = 0.0;  gzz0[idx] = 0.0;

                Ex0[idx] = sin(y) * cos(z);
                Ey0[idx] = cos(x) * sin(z);
                Ez0[idx] = sin(x) * cos(y);

                Bx0[idx] = cos(y) * sin(z);
                By0[idx] = sin(x) * cos(z);
                Bz0[idx] = cos(x) * sin(y);
            }
        }
    }
}

// ===========================================================================
// Test 3: NP scalar Phi2 (gpu_getnpem2_launch vs f_getnpem2)
// ===========================================================================
static int test_np_phi2()
{
    printf("[3/4] NP scalar Phi2 (getnpem2)\n");

    int ex[3]  = {8 + 2*GW, 8 + 2*GW, 8 + 2*GW};
    int pex[3] = {ex[0], ex[1], ex[2]};
    int n = ex[0] * ex[1] * ex[2];

    // ---- allocate host arrays ----
    double *X    = new double[n];
    double *Y    = new double[n];
    double *Z    = new double[n];
    double *chi  = new double[n];
    double *gxx0 = new double[n]; double *gxy0 = new double[n]; double *gxz0 = new double[n];
    double *gyy0 = new double[n]; double *gyz0 = new double[n]; double *gzz0 = new double[n];
    double *Ex0  = new double[n]; double *Ey0  = new double[n]; double *Ez0  = new double[n];
    double *Bx0  = new double[n]; double *By0  = new double[n]; double *Bz0  = new double[n];

    fill_np_test_data(ex, X, Y, Z, chi,
                      gxx0, gxy0, gxz0, gyy0, gyz0, gzz0,
                      Ex0, Ey0, Ez0, Bx0, By0, Bz0);

    // ---- allocate GPU memory and copy inputs ----
    size_t nb = n * sizeof(double);
    double *d_X, *d_Y, *d_Z, *d_chi;
    double *d_gxx0, *d_gxy0, *d_gxz0, *d_gyy0, *d_gyz0, *d_gzz0;
    double *d_Ex0, *d_Ey0, *d_Ez0, *d_Bx0, *d_By0, *d_Bz0;
    double *d_Rphi2, *d_Iphi2;

    CUDA_CHECK(cudaMalloc(&d_X,  nb)); CUDA_CHECK(cudaMalloc(&d_Y,  nb));
    CUDA_CHECK(cudaMalloc(&d_Z,  nb)); CUDA_CHECK(cudaMalloc(&d_chi, nb));
    CUDA_CHECK(cudaMalloc(&d_gxx0, nb)); CUDA_CHECK(cudaMalloc(&d_gxy0, nb));
    CUDA_CHECK(cudaMalloc(&d_gxz0, nb)); CUDA_CHECK(cudaMalloc(&d_gyy0, nb));
    CUDA_CHECK(cudaMalloc(&d_gyz0, nb)); CUDA_CHECK(cudaMalloc(&d_gzz0, nb));
    CUDA_CHECK(cudaMalloc(&d_Ex0, nb)); CUDA_CHECK(cudaMalloc(&d_Ey0, nb));
    CUDA_CHECK(cudaMalloc(&d_Ez0, nb)); CUDA_CHECK(cudaMalloc(&d_Bx0, nb));
    CUDA_CHECK(cudaMalloc(&d_By0, nb)); CUDA_CHECK(cudaMalloc(&d_Bz0, nb));
    CUDA_CHECK(cudaMalloc(&d_Rphi2, nb)); CUDA_CHECK(cudaMalloc(&d_Iphi2, nb));

    cudaStream_t stream = 0;
    CUDA_CHECK(cudaMemcpy(d_X,    X,    nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Y,    Y,    nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Z,    Z,    nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_chi,  chi,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gxx0, gxx0, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gxy0, gxy0, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gxz0, gxz0, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gyy0, gyy0, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gyz0, gyz0, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gzz0, gzz0, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Ex0,  Ex0,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Ey0,  Ey0,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Ez0,  Ez0,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Bx0,  Bx0,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_By0,  By0,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Bz0,  Bz0,  nb, cudaMemcpyHostToDevice));

    // ---- run GPU kernel ----
    // GPU parameter order (from getnpem2.h):
    //   stream, ex, pex, X,Y,Z, chi,
    //   gxx,gxy,gxz,gyy,gyz,gzz,
    //   Ex,Ey,Ez,Bx,By,Bz,
    //   Rphi2, Iphi2, symmetry
    gpu_getnpem2_launch(
        stream, ex, pex,
        d_X, d_Y, d_Z, d_chi,
        d_gxx0, d_gxy0, d_gxz0, d_gyy0, d_gyz0, d_gzz0,
        d_Ex0, d_Ey0, d_Ez0, d_Bx0, d_By0, d_Bz0,
        d_Rphi2, d_Iphi2,
        0 /*symmetry*/);
    CUDA_CHECK(cudaDeviceSynchronize());

    // ---- copy GPU outputs back ----
    double *gpu_Rphi2 = new double[n];
    double *gpu_Iphi2 = new double[n];
    CUDA_CHECK(cudaMemcpy(gpu_Rphi2, d_Rphi2, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Iphi2, d_Iphi2, nb, cudaMemcpyDeviceToHost));

    // ---- run CPU Fortran function ----
    // CPU parameter order (from getnpem2.h / bssnEM_class.C:1884-1891):
    //   ex, X,Y,Z, chi, gxx,gxy,gxz,gyy,gyz,gzz,
    //   Ex,Ey,Ez,Bx,By,Bz,
    //   Rphi2, Iphi2, symmetry (reference)
    double *cpu_Rphi2 = new double[n];
    double *cpu_Iphi2 = new double[n];
    int sym_ref = 0;

    f_getnpem2(ex, X, Y, Z,
               chi,
               gxx0, gxy0, gxz0, gyy0, gyz0, gzz0,
               Ex0, Ey0, Ez0, Bx0, By0, Bz0,
               cpu_Rphi2, cpu_Iphi2,
               sym_ref);

    // ---- compare (interior only) ----
    int fails = 0;
    fails += compare_field("Rphi2", gpu_Rphi2, cpu_Rphi2, ex, TOL);
    fails += compare_field("Iphi2", gpu_Iphi2, cpu_Iphi2, ex, TOL);

    // ---- cleanup ----
    delete[] X; delete[] Y; delete[] Z; delete[] chi;
    delete[] gxx0; delete[] gxy0; delete[] gxz0;
    delete[] gyy0; delete[] gyz0; delete[] gzz0;
    delete[] Ex0;  delete[] Ey0;  delete[] Ez0;
    delete[] Bx0;  delete[] By0;  delete[] Bz0;
    delete[] gpu_Rphi2; delete[] gpu_Iphi2;
    delete[] cpu_Rphi2; delete[] cpu_Iphi2;

    CUDA_CHECK(cudaFree(d_X));  CUDA_CHECK(cudaFree(d_Y));  CUDA_CHECK(cudaFree(d_Z));
    CUDA_CHECK(cudaFree(d_chi));
    CUDA_CHECK(cudaFree(d_gxx0)); CUDA_CHECK(cudaFree(d_gxy0)); CUDA_CHECK(cudaFree(d_gxz0));
    CUDA_CHECK(cudaFree(d_gyy0)); CUDA_CHECK(cudaFree(d_gyz0)); CUDA_CHECK(cudaFree(d_gzz0));
    CUDA_CHECK(cudaFree(d_Ex0)); CUDA_CHECK(cudaFree(d_Ey0)); CUDA_CHECK(cudaFree(d_Ez0));
    CUDA_CHECK(cudaFree(d_Bx0)); CUDA_CHECK(cudaFree(d_By0)); CUDA_CHECK(cudaFree(d_Bz0));
    CUDA_CHECK(cudaFree(d_Rphi2)); CUDA_CHECK(cudaFree(d_Iphi2));

    return fails;
}

// ===========================================================================
// Test 4: NP scalar Phi1 (gpu_getnpem1_launch vs f_getnpem1)
// ===========================================================================
static int test_np_phi1()
{
    printf("[4/4] NP scalar Phi1 (getnpem1)\n");

    int ex[3]  = {8 + 2*GW, 8 + 2*GW, 8 + 2*GW};
    int pex[3] = {ex[0], ex[1], ex[2]};
    int n = ex[0] * ex[1] * ex[2];

    // ---- allocate host arrays ----
    double *X    = new double[n];
    double *Y    = new double[n];
    double *Z    = new double[n];
    double *chi  = new double[n];
    double *gxx0 = new double[n]; double *gxy0 = new double[n]; double *gxz0 = new double[n];
    double *gyy0 = new double[n]; double *gyz0 = new double[n]; double *gzz0 = new double[n];
    double *Ex0  = new double[n]; double *Ey0  = new double[n]; double *Ez0  = new double[n];
    double *Bx0  = new double[n]; double *By0  = new double[n]; double *Bz0  = new double[n];

    fill_np_test_data(ex, X, Y, Z, chi,
                      gxx0, gxy0, gxz0, gyy0, gyz0, gzz0,
                      Ex0, Ey0, Ez0, Bx0, By0, Bz0);

    // ---- allocate GPU memory and copy inputs ----
    size_t nb = n * sizeof(double);
    double *d_X, *d_Y, *d_Z, *d_chi;
    double *d_gxx0, *d_gxy0, *d_gxz0, *d_gyy0, *d_gyz0, *d_gzz0;
    double *d_Ex0, *d_Ey0, *d_Ez0, *d_Bx0, *d_By0, *d_Bz0;
    double *d_Rphi1, *d_Iphi1;

    CUDA_CHECK(cudaMalloc(&d_X,  nb)); CUDA_CHECK(cudaMalloc(&d_Y,  nb));
    CUDA_CHECK(cudaMalloc(&d_Z,  nb)); CUDA_CHECK(cudaMalloc(&d_chi, nb));
    CUDA_CHECK(cudaMalloc(&d_gxx0, nb)); CUDA_CHECK(cudaMalloc(&d_gxy0, nb));
    CUDA_CHECK(cudaMalloc(&d_gxz0, nb)); CUDA_CHECK(cudaMalloc(&d_gyy0, nb));
    CUDA_CHECK(cudaMalloc(&d_gyz0, nb)); CUDA_CHECK(cudaMalloc(&d_gzz0, nb));
    CUDA_CHECK(cudaMalloc(&d_Ex0, nb)); CUDA_CHECK(cudaMalloc(&d_Ey0, nb));
    CUDA_CHECK(cudaMalloc(&d_Ez0, nb)); CUDA_CHECK(cudaMalloc(&d_Bx0, nb));
    CUDA_CHECK(cudaMalloc(&d_By0, nb)); CUDA_CHECK(cudaMalloc(&d_Bz0, nb));
    CUDA_CHECK(cudaMalloc(&d_Rphi1, nb)); CUDA_CHECK(cudaMalloc(&d_Iphi1, nb));

    cudaStream_t stream = 0;
    CUDA_CHECK(cudaMemcpy(d_X,    X,    nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Y,    Y,    nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Z,    Z,    nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_chi,  chi,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gxx0, gxx0, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gxy0, gxy0, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gxz0, gxz0, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gyy0, gyy0, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gyz0, gyz0, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gzz0, gzz0, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Ex0,  Ex0,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Ey0,  Ey0,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Ez0,  Ez0,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Bx0,  Bx0,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_By0,  By0,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Bz0,  Bz0,  nb, cudaMemcpyHostToDevice));

    // ---- run GPU kernel ----
    gpu_getnpem1_launch(
        stream, ex, pex,
        d_X, d_Y, d_Z, d_chi,
        d_gxx0, d_gxy0, d_gxz0, d_gyy0, d_gyz0, d_gzz0,
        d_Ex0, d_Ey0, d_Ez0, d_Bx0, d_By0, d_Bz0,
        d_Rphi1, d_Iphi1,
        0 /*symmetry*/);
    CUDA_CHECK(cudaDeviceSynchronize());

    // ---- copy GPU outputs back ----
    double *gpu_Rphi1 = new double[n];
    double *gpu_Iphi1 = new double[n];
    CUDA_CHECK(cudaMemcpy(gpu_Rphi1, d_Rphi1, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Iphi1, d_Iphi1, nb, cudaMemcpyDeviceToHost));

    // ---- run CPU Fortran function ----
    double *cpu_Rphi1 = new double[n];
    double *cpu_Iphi1 = new double[n];
    int sym_ref = 0;

    f_getnpem1(ex, X, Y, Z,
               chi,
               gxx0, gxy0, gxz0, gyy0, gyz0, gzz0,
               Ex0, Ey0, Ez0, Bx0, By0, Bz0,
               cpu_Rphi1, cpu_Iphi1,
               sym_ref);

    // ---- compare (interior only) ----
    int fails = 0;
    fails += compare_field("Rphi1", gpu_Rphi1, cpu_Rphi1, ex, TOL);
    fails += compare_field("Iphi1", gpu_Iphi1, cpu_Iphi1, ex, TOL);

    // ---- cleanup ----
    delete[] X; delete[] Y; delete[] Z; delete[] chi;
    delete[] gxx0; delete[] gxy0; delete[] gxz0;
    delete[] gyy0; delete[] gyz0; delete[] gzz0;
    delete[] Ex0;  delete[] Ey0;  delete[] Ez0;
    delete[] Bx0;  delete[] By0;  delete[] Bz0;
    delete[] gpu_Rphi1; delete[] gpu_Iphi1;
    delete[] cpu_Rphi1; delete[] cpu_Iphi1;

    CUDA_CHECK(cudaFree(d_X));  CUDA_CHECK(cudaFree(d_Y));  CUDA_CHECK(cudaFree(d_Z));
    CUDA_CHECK(cudaFree(d_chi));
    CUDA_CHECK(cudaFree(d_gxx0)); CUDA_CHECK(cudaFree(d_gxy0)); CUDA_CHECK(cudaFree(d_gxz0));
    CUDA_CHECK(cudaFree(d_gyy0)); CUDA_CHECK(cudaFree(d_gyz0)); CUDA_CHECK(cudaFree(d_gzz0));
    CUDA_CHECK(cudaFree(d_Ex0)); CUDA_CHECK(cudaFree(d_Ey0)); CUDA_CHECK(cudaFree(d_Ez0));
    CUDA_CHECK(cudaFree(d_Bx0)); CUDA_CHECK(cudaFree(d_By0)); CUDA_CHECK(cudaFree(d_Bz0));
    CUDA_CHECK(cudaFree(d_Rphi1)); CUDA_CHECK(cudaFree(d_Iphi1));

    return fails;
}

// ===========================================================================
// fill_bssn_test_data — fill BSSN field arrays with smooth analytic functions
// ===========================================================================
static void fill_bssn_test_data(int ex[3],
                                double *X, double *Y, double *Z,
                                double *chi, double *trK,
                                double *gxx, double *gxy, double *gxz,
                                double *gyy, double *gyz, double *gzz,
                                double *Axx, double *Axy, double *Axz,
                                double *Ayy, double *Ayz, double *Azz,
                                double *Gamx, double *Gamy, double *Gamz,
                                double *Lap, double *betax, double *betay, double *betaz,
                                double *dtSfx, double *dtSfy, double *dtSfz)
{
    double dx = 0.1;
    for (int k = 0; k < ex[2]; k++) {
        for (int j = 0; j < ex[1]; j++) {
            for (int i = 0; i < ex[0]; i++) {
                int idx = i + j * ex[0] + k * ex[0] * ex[1];
                double x = (i - GW) * dx;
                double y = (j - GW) * dx;
                double z = (k - GW) * dx;

                X[idx] = x;
                Y[idx] = y;
                Z[idx] = z;

                // Small perturbations from flat Minkowski
                chi[idx]   = 1.0 + 0.01 * sin(x) * cos(y);
                trK[idx]   = 0.001 * sin(z);

                gxx[idx]   = 0.0;
                gxy[idx]   = 0.01 * sin(x) * cos(y) * sin(z);
                gxz[idx]   = 0.005 * cos(x) * sin(y) * cos(z);
                gyy[idx]   = 0.0;
                gyz[idx]   = 0.008 * sin(x) * sin(y) * cos(z);
                gzz[idx]   = 0.0;

                Axx[idx]   = 0.001 * cos(x) * sin(y) * sin(z);
                Axy[idx]   = 0.002 * sin(x) * sin(y) * cos(z);
                Axz[idx]   = 0.0015 * cos(x) * cos(y) * sin(z);
                Ayy[idx]   = -0.001 * cos(x) * sin(y) * sin(z);  // trace-free (diag sum zero)
                Ayz[idx]   = 0.003 * sin(x) * cos(y) * cos(z);
                Azz[idx]   = 0.0;

                Gamx[idx]  = 0.001 * sin(x + y) * cos(z);
                Gamy[idx]  = 0.001 * cos(x) * sin(y + z);
                Gamz[idx]  = 0.001 * cos(x + z) * sin(y);

                Lap[idx]   = 1.0;
                betax[idx] = 0.001 * sin(y) * cos(z);
                betay[idx] = 0.001 * cos(x) * sin(z);
                betaz[idx] = 0.001 * sin(x) * cos(y);
                dtSfx[idx] = 0.0;
                dtSfy[idx] = 0.0;
                dtSfz[idx] = 0.0;
            }
        }
    }
}

// ===========================================================================
// fill_puncture_test_data — fill BSSN fields with Brill-Lindquist puncture data
//
// Conformal factor: chi = psi^(-4) - 1, where psi = 1 + M/(2r)
// Near the puncture (r -> 0): psi is large, chi -> -1
// Far from the puncture (r -> inf): psi -> 1, chi -> 0
// ===========================================================================
static void fill_puncture_test_data(int ex[3], double M, double h,
                                     double *X, double *Y, double *Z,
                                     double *chi, double *trK,
                                     double *gxx, double *gxy, double *gxz,
                                     double *gyy, double *gyz, double *gzz,
                                     double *Axx, double *Axy, double *Axz,
                                     double *Ayy, double *Ayz, double *Azz,
                                     double *Gamx, double *Gamy, double *Gamz,
                                     double *Lap, double *betax, double *betay, double *betaz,
                                     double *dtSfx, double *dtSfy, double *dtSfz)
{
    int cx = ex[0] / 2;
    int cy = ex[1] / 2;
    int cz = ex[2] / 2;

    for (int k = 0; k < ex[2]; k++) {
        for (int j = 0; j < ex[1]; j++) {
            for (int i = 0; i < ex[0]; i++) {
                int idx = i + j * ex[0] + k * ex[0] * ex[1];
                double dx = (i - cx) * h;
                double dy = (j - cy) * h;
                double dz = (k - cz) * h;
                double r = sqrt(dx*dx + dy*dy + dz*dz) + 1e-6;

                // Brill-Lindquist puncture: psi = 1 + M/(2r), chi = psi^(-4) - 1
                double psi = 1.0 + M / (2.0 * r);
                double chi_val = 1.0 / (psi * psi * psi * psi) - 1.0;

                X[idx] = dx;
                Y[idx] = dy;
                Z[idx] = dz;

                chi[idx]   = chi_val;
                trK[idx]   = 0.001 * exp(-r);

                gxx[idx]   = 0.0;
                gxy[idx]   = 0.0;
                gxz[idx]   = 0.0;
                gyy[idx]   = 0.0;
                gyz[idx]   = 0.0;
                gzz[idx]   = 0.0;

                // Trace-free perturbation on Aij (diag sum = 0)
                double amp = 0.001 * exp(-r);
                double n2 = dx*dx + dy*dy + dz*dz + 1e-12;
                Axx[idx]   = amp * (dx*dx - dy*dy) / n2;
                Axy[idx]   = amp * dx*dy / n2;
                Axz[idx]   = amp * dx*dz / n2;
                Ayy[idx]   = amp * (dy*dy - dz*dz) / n2;
                Ayz[idx]   = amp * dy*dz / n2;
                Azz[idx]   = amp * (dz*dz - dx*dx) / n2;

                Gamx[idx]  = 0.001 * exp(-r) * dx / (r + 1e-6);
                Gamy[idx]  = 0.001 * exp(-r) * dy / (r + 1e-6);
                Gamz[idx]  = 0.001 * exp(-r) * dz / (r + 1e-6);

                Lap[idx]   = 1.0;
                betax[idx] = 0.001 * exp(-r) * dx / (r + 1e-6);
                betay[idx] = 0.001 * exp(-r) * dy / (r + 1e-6);
                betaz[idx] = 0.001 * exp(-r) * dz / (r + 1e-6);
                dtSfx[idx] = 0.0;
                dtSfy[idx] = 0.0;
                dtSfz[idx] = 0.0;
            }
        }
    }
}

// ===========================================================================
// fill_stress_energy_test_data — fill stress-energy arrays
// ===========================================================================
static void fill_stress_energy_test_data(int ex[3],
                                          double *rho, double *Sx, double *Sy, double *Sz,
                                          double *Sxx, double *Sxy, double *Sxz,
                                          double *Syy, double *Syz, double *Szz)
{
    double dx = 0.1;
    for (int k = 0; k < ex[2]; k++) {
        for (int j = 0; j < ex[1]; j++) {
            for (int i = 0; i < ex[0]; i++) {
                int idx = i + j * ex[0] + k * ex[0] * ex[1];
                double x = (i - GW) * dx;
                double y = (j - GW) * dx;
                double z = (k - GW) * dx;

                rho[idx] = 0.001 * sin(x) * sin(y) * sin(z);
                Sx[idx]  = 0.001 * cos(x) * sin(y) * sin(z);
                Sy[idx]  = 0.001 * sin(x) * cos(y) * sin(z);
                Sz[idx]  = 0.001 * sin(x) * sin(y) * cos(z);

                Sxx[idx] = 0.0005 * (1.0 + sin(x));
                Sxy[idx] = 0.0005 * cos(x) * cos(y);
                Sxz[idx] = 0.0005 * cos(x) * cos(z);
                Syy[idx] = 0.0005 * (1.0 + sin(y));
                Syz[idx] = 0.0005 * cos(y) * cos(z);
                Szz[idx] = 0.0005 * (1.0 + sin(z));
            }
        }
    }
}

// ===========================================================================
// fill_christoffel_test_data — fill Christoffel symbol arrays (18 components)
// ===========================================================================
static void fill_christoffel_test_data(int ex[3],
                                        double *Gamxxx, double *Gamxxy, double *Gamxxz,
                                        double *Gamxyy, double *Gamxyz, double *Gamxzz,
                                        double *Gamyxx, double *Gamyxy, double *Gamyxz,
                                        double *Gamyyy, double *Gamyyz, double *Gamyzz,
                                        double *Gamzxx, double *Gamzxy, double *Gamzxz,
                                        double *Gamzyy, double *Gamzyz, double *Gamzzz)
{
    double dx = 0.1;
    for (int k = 0; k < ex[2]; k++) {
        for (int j = 0; j < ex[1]; j++) {
            for (int i = 0; i < ex[0]; i++) {
                int idx = i + j * ex[0] + k * ex[0] * ex[1];
                double x = (i - GW) * dx;
                double y = (j - GW) * dx;
                double z = (k - GW) * dx;

                // Fill all 18 with small perturbations (for Minkowski they'd be zero)
                Gamxxx[idx] = 0.0001 * cos(x) * sin(y) * cos(z);
                Gamxxy[idx] = 0.0001 * sin(x) * cos(y) * sin(z);
                Gamxxz[idx] = 0.0001 * cos(x) * cos(y) * sin(z);
                Gamxyy[idx] = 0.0001 * sin(x) * sin(y) * cos(z);
                Gamxyz[idx] = 0.0001 * cos(x) * sin(y) * sin(z);
                Gamxzz[idx] = 0.0001 * sin(x) * cos(y) * cos(z);

                Gamyxx[idx] = 0.0001 * sin(x) * sin(y) * cos(z);
                Gamyxy[idx] = 0.0001 * cos(x) * cos(y) * sin(z);
                Gamyxz[idx] = 0.0001 * sin(x) * cos(y) * sin(z);
                Gamyyy[idx] = 0.0001 * cos(x) * sin(y) * cos(z);
                Gamyyz[idx] = 0.0001 * sin(x) * sin(y) * sin(z);
                Gamyzz[idx] = 0.0001 * cos(x) * cos(y) * cos(z);

                Gamzxx[idx] = 0.0001 * cos(x) * sin(y) * sin(z);
                Gamzxy[idx] = 0.0001 * sin(x) * cos(y) * cos(z);
                Gamzxz[idx] = 0.0001 * cos(x) * cos(y) * sin(z);
                Gamzyy[idx] = 0.0001 * sin(x) * sin(y) * cos(z);
                Gamzyz[idx] = 0.0001 * cos(x) * sin(y) * sin(z);
                Gamzzz[idx] = 0.0001 * sin(x) * cos(y) * cos(z);
            }
        }
    }
}

// ===========================================================================
// fill_ricci_test_data — fill Ricci tensor arrays (6 components)
// ===========================================================================
static void fill_ricci_test_data(int ex[3],
                                  double *Rxx, double *Rxy, double *Rxz,
                                  double *Ryy, double *Ryz, double *Rzz)
{
    double dx = 0.1;
    for (int k = 0; k < ex[2]; k++) {
        for (int j = 0; j < ex[1]; j++) {
            for (int i = 0; i < ex[0]; i++) {
                int idx = i + j * ex[0] + k * ex[0] * ex[1];
                double x = (i - GW) * dx;
                double y = (j - GW) * dx;
                double z = (k - GW) * dx;

                Rxx[idx] = 0.0001 * sin(x) * cos(y) * cos(z);
                Rxy[idx] = 0.0001 * cos(x) * sin(y) * cos(z);
                Rxz[idx] = 0.0001 * cos(x) * cos(y) * sin(z);
                Ryy[idx] = 0.0001 * sin(x) * sin(y) * cos(z);
                Ryz[idx] = 0.0001 * sin(x) * cos(y) * sin(z);
                Rzz[idx] = 0.0001 * cos(x) * sin(y) * sin(z);
            }
        }
    }
}

// ===========================================================================
// Test 5: BSSN RHS (gpu_compute_rhs_bssn_launch vs f_compute_rhs_bssn)
// ===========================================================================
static int test_bssn_rhs()
{
    printf("[5/7] BSSN RHS (compute_rhs_bssn)\n");

    int ex[3]  = {8 + 2*GW, 8 + 2*GW, 8 + 2*GW};
    int pex[3] = {ex[0], ex[1], ex[2]};
    int n = ex[0] * ex[1] * ex[2];

    // ---- allocate host arrays ----
    double *X    = new double[n];
    double *Y    = new double[n];
    double *Z    = new double[n];

    // 23 BSSN state arrays
    double *chi   = new double[n];
    double *trK   = new double[n];
    double *gxx   = new double[n]; double *gxy   = new double[n]; double *gxz   = new double[n];
    double *gyy   = new double[n]; double *gyz   = new double[n]; double *gzz   = new double[n];
    double *Axx   = new double[n]; double *Axy   = new double[n]; double *Axz   = new double[n];
    double *Ayy   = new double[n]; double *Ayz   = new double[n]; double *Azz   = new double[n];
    double *Gamx  = new double[n]; double *Gamy  = new double[n]; double *Gamz  = new double[n];
    double *Lap   = new double[n];
    double *betax = new double[n]; double *betay = new double[n]; double *betaz = new double[n];
    double *dtSfx = new double[n]; double *dtSfy = new double[n]; double *dtSfz = new double[n];

    // 10 stress-energy arrays
    double *rho = new double[n];
    double *Sx  = new double[n]; double *Sy  = new double[n]; double *Sz  = new double[n];
    double *Sxx = new double[n]; double *Sxy = new double[n]; double *Sxz = new double[n];
    double *Syy = new double[n]; double *Syz = new double[n]; double *Szz = new double[n];

    // 18 Christoffel arrays
    double *Gamxxx = new double[n]; double *Gamxxy = new double[n]; double *Gamxxz = new double[n];
    double *Gamxyy = new double[n]; double *Gamxyz = new double[n]; double *Gamxzz = new double[n];
    double *Gamyxx = new double[n]; double *Gamyxy = new double[n]; double *Gamyxz = new double[n];
    double *Gamyyy = new double[n]; double *Gamyyz = new double[n]; double *Gamyzz = new double[n];
    double *Gamzxx = new double[n]; double *Gamzxy = new double[n]; double *Gamzxz = new double[n];
    double *Gamzyy = new double[n]; double *Gamzyz = new double[n]; double *Gamzzz = new double[n];

    // 6 Ricci arrays
    double *Ricxx = new double[n]; double *Ricxy = new double[n]; double *Ricxz = new double[n];
    double *Ricyy = new double[n]; double *Ricyz = new double[n]; double *Riczz = new double[n];

    // Fill all input arrays
    fill_bssn_test_data(ex, X, Y, Z,
                        chi, trK,
                        gxx, gxy, gxz, gyy, gyz, gzz,
                        Axx, Axy, Axz, Ayy, Ayz, Azz,
                        Gamx, Gamy, Gamz,
                        Lap, betax, betay, betaz,
                        dtSfx, dtSfy, dtSfz);
    fill_stress_energy_test_data(ex, rho, Sx, Sy, Sz,
                                  Sxx, Sxy, Sxz, Syy, Syz, Szz);
    fill_christoffel_test_data(ex,
                                Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
                                Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
                                Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz);
    fill_ricci_test_data(ex, Ricxx, Ricxy, Ricxz, Ricyy, Ricyz, Riczz);

    // ---- allocate GPU memory ----
    size_t nb = n * sizeof(double);
    double *d_X, *d_Y, *d_Z;
    double *d_chi, *d_trK;
    double *d_gxx, *d_gxy, *d_gxz, *d_gyy, *d_gyz, *d_gzz;
    double *d_Axx, *d_Axy, *d_Axz, *d_Ayy, *d_Ayz, *d_Azz;
    double *d_Gamx, *d_Gamy, *d_Gamz, *d_Lap;
    double *d_betax, *d_betay, *d_betaz;
    double *d_dtSfx, *d_dtSfy, *d_dtSfz;
    double *d_chi_rhs, *d_trK_rhs;
    double *d_gxx_rhs, *d_gxy_rhs, *d_gxz_rhs, *d_gyy_rhs, *d_gyz_rhs, *d_gzz_rhs;
    double *d_Axx_rhs, *d_Axy_rhs, *d_Axz_rhs, *d_Ayy_rhs, *d_Ayz_rhs, *d_Azz_rhs;
    double *d_Gamx_rhs, *d_Gamy_rhs, *d_Gamz_rhs, *d_Lap_rhs;
    double *d_betax_rhs, *d_betay_rhs, *d_betaz_rhs;
    double *d_dtSfx_rhs, *d_dtSfy_rhs, *d_dtSfz_rhs;
    double *d_rho, *d_Sx, *d_Sy, *d_Sz;
    double *d_Sxx, *d_Sxy, *d_Sxz, *d_Syy, *d_Syz, *d_Szz;
    double *d_Gamxxx, *d_Gamxxy, *d_Gamxxz, *d_Gamxyy, *d_Gamxyz, *d_Gamxzz;
    double *d_Gamyxx, *d_Gamyxy, *d_Gamyxz, *d_Gamyyy, *d_Gamyyz, *d_Gamyzz;
    double *d_Gamzxx, *d_Gamzxy, *d_Gamzxz, *d_Gamzyy, *d_Gamzyz, *d_Gamzzz;
    double *d_Rxx, *d_Rxy, *d_Rxz, *d_Ryy, *d_Ryz, *d_Rzz;
    double *d_ham, *d_movx, *d_movy, *d_movz, *d_Gmx, *d_Gmy, *d_Gmz;

    CUDA_CHECK(cudaMalloc(&d_X,    nb)); CUDA_CHECK(cudaMalloc(&d_Y,    nb)); CUDA_CHECK(cudaMalloc(&d_Z,    nb));
    CUDA_CHECK(cudaMalloc(&d_chi,  nb)); CUDA_CHECK(cudaMalloc(&d_trK,  nb));
    CUDA_CHECK(cudaMalloc(&d_gxx,  nb)); CUDA_CHECK(cudaMalloc(&d_gxy,  nb)); CUDA_CHECK(cudaMalloc(&d_gxz,  nb));
    CUDA_CHECK(cudaMalloc(&d_gyy,  nb)); CUDA_CHECK(cudaMalloc(&d_gyz,  nb)); CUDA_CHECK(cudaMalloc(&d_gzz,  nb));
    CUDA_CHECK(cudaMalloc(&d_Axx,  nb)); CUDA_CHECK(cudaMalloc(&d_Axy,  nb)); CUDA_CHECK(cudaMalloc(&d_Axz,  nb));
    CUDA_CHECK(cudaMalloc(&d_Ayy,  nb)); CUDA_CHECK(cudaMalloc(&d_Ayz,  nb)); CUDA_CHECK(cudaMalloc(&d_Azz,  nb));
    CUDA_CHECK(cudaMalloc(&d_Gamx, nb)); CUDA_CHECK(cudaMalloc(&d_Gamy, nb)); CUDA_CHECK(cudaMalloc(&d_Gamz, nb));
    CUDA_CHECK(cudaMalloc(&d_Lap,  nb));
    CUDA_CHECK(cudaMalloc(&d_betax, nb)); CUDA_CHECK(cudaMalloc(&d_betay, nb)); CUDA_CHECK(cudaMalloc(&d_betaz, nb));
    CUDA_CHECK(cudaMalloc(&d_dtSfx, nb)); CUDA_CHECK(cudaMalloc(&d_dtSfy, nb)); CUDA_CHECK(cudaMalloc(&d_dtSfz, nb));
    CUDA_CHECK(cudaMalloc(&d_chi_rhs,  nb)); CUDA_CHECK(cudaMalloc(&d_trK_rhs,  nb));
    CUDA_CHECK(cudaMalloc(&d_gxx_rhs,  nb)); CUDA_CHECK(cudaMalloc(&d_gxy_rhs,  nb)); CUDA_CHECK(cudaMalloc(&d_gxz_rhs,  nb));
    CUDA_CHECK(cudaMalloc(&d_gyy_rhs,  nb)); CUDA_CHECK(cudaMalloc(&d_gyz_rhs,  nb)); CUDA_CHECK(cudaMalloc(&d_gzz_rhs,  nb));
    CUDA_CHECK(cudaMalloc(&d_Axx_rhs,  nb)); CUDA_CHECK(cudaMalloc(&d_Axy_rhs,  nb)); CUDA_CHECK(cudaMalloc(&d_Axz_rhs,  nb));
    CUDA_CHECK(cudaMalloc(&d_Ayy_rhs,  nb)); CUDA_CHECK(cudaMalloc(&d_Ayz_rhs,  nb)); CUDA_CHECK(cudaMalloc(&d_Azz_rhs,  nb));
    CUDA_CHECK(cudaMalloc(&d_Gamx_rhs, nb)); CUDA_CHECK(cudaMalloc(&d_Gamy_rhs, nb)); CUDA_CHECK(cudaMalloc(&d_Gamz_rhs, nb));
    CUDA_CHECK(cudaMalloc(&d_Lap_rhs,  nb));
    CUDA_CHECK(cudaMalloc(&d_betax_rhs, nb)); CUDA_CHECK(cudaMalloc(&d_betay_rhs, nb)); CUDA_CHECK(cudaMalloc(&d_betaz_rhs, nb));
    CUDA_CHECK(cudaMalloc(&d_dtSfx_rhs, nb)); CUDA_CHECK(cudaMalloc(&d_dtSfy_rhs, nb)); CUDA_CHECK(cudaMalloc(&d_dtSfz_rhs, nb));
    CUDA_CHECK(cudaMalloc(&d_rho,   nb)); CUDA_CHECK(cudaMalloc(&d_Sx,  nb)); CUDA_CHECK(cudaMalloc(&d_Sy,  nb)); CUDA_CHECK(cudaMalloc(&d_Sz, nb));
    CUDA_CHECK(cudaMalloc(&d_Sxx,   nb)); CUDA_CHECK(cudaMalloc(&d_Sxy, nb)); CUDA_CHECK(cudaMalloc(&d_Sxz, nb));
    CUDA_CHECK(cudaMalloc(&d_Syy,   nb)); CUDA_CHECK(cudaMalloc(&d_Syz, nb)); CUDA_CHECK(cudaMalloc(&d_Szz, nb));
    CUDA_CHECK(cudaMalloc(&d_Gamxxx, nb)); CUDA_CHECK(cudaMalloc(&d_Gamxxy, nb)); CUDA_CHECK(cudaMalloc(&d_Gamxxz, nb));
    CUDA_CHECK(cudaMalloc(&d_Gamxyy, nb)); CUDA_CHECK(cudaMalloc(&d_Gamxyz, nb)); CUDA_CHECK(cudaMalloc(&d_Gamxzz, nb));
    CUDA_CHECK(cudaMalloc(&d_Gamyxx, nb)); CUDA_CHECK(cudaMalloc(&d_Gamyxy, nb)); CUDA_CHECK(cudaMalloc(&d_Gamyxz, nb));
    CUDA_CHECK(cudaMalloc(&d_Gamyyy, nb)); CUDA_CHECK(cudaMalloc(&d_Gamyyz, nb)); CUDA_CHECK(cudaMalloc(&d_Gamyzz, nb));
    CUDA_CHECK(cudaMalloc(&d_Gamzxx, nb)); CUDA_CHECK(cudaMalloc(&d_Gamzxy, nb)); CUDA_CHECK(cudaMalloc(&d_Gamzxz, nb));
    CUDA_CHECK(cudaMalloc(&d_Gamzyy, nb)); CUDA_CHECK(cudaMalloc(&d_Gamzyz, nb)); CUDA_CHECK(cudaMalloc(&d_Gamzzz, nb));
    CUDA_CHECK(cudaMalloc(&d_Rxx,   nb)); CUDA_CHECK(cudaMalloc(&d_Rxy,   nb)); CUDA_CHECK(cudaMalloc(&d_Rxz,   nb));
    CUDA_CHECK(cudaMalloc(&d_Ryy,   nb)); CUDA_CHECK(cudaMalloc(&d_Ryz,   nb)); CUDA_CHECK(cudaMalloc(&d_Rzz,   nb));
    CUDA_CHECK(cudaMalloc(&d_ham,   nb)); CUDA_CHECK(cudaMalloc(&d_movx,  nb));
    CUDA_CHECK(cudaMalloc(&d_movy,  nb)); CUDA_CHECK(cudaMalloc(&d_movz,  nb));
    CUDA_CHECK(cudaMalloc(&d_Gmx,   nb)); CUDA_CHECK(cudaMalloc(&d_Gmy,   nb)); CUDA_CHECK(cudaMalloc(&d_Gmz,  nb));

    // ---- copy inputs to GPU ----
    cudaStream_t stream = 0;
    CUDA_CHECK(cudaMemcpy(d_X,    X,    nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Y,    Y,    nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Z,    Z,    nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_chi,  chi,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_trK,  trK,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gxx,  gxx,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gxy,  gxy,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gxz,  gxz,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gyy,  gyy,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gyz,  gyz,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gzz,  gzz,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Axx,  Axx,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Axy,  Axy,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Axz,  Axz,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Ayy,  Ayy,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Ayz,  Ayz,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Azz,  Azz,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamx, Gamx, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamy, Gamy, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamz, Gamz, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Lap,  Lap,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_betax, betax, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_betay, betay, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_betaz, betaz, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_dtSfx, dtSfx, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_dtSfy, dtSfy, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_dtSfz, dtSfz, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_rho,   rho,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Sx,    Sx,    nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Sy,    Sy,    nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Sz,    Sz,    nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Sxx,   Sxx,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Sxy,   Sxy,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Sxz,   Sxz,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Syy,   Syy,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Syz,   Syz,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Szz,   Szz,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamxxx, Gamxxx, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamxxy, Gamxxy, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamxxz, Gamxxz, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamxyy, Gamxyy, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamxyz, Gamxyz, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamxzz, Gamxzz, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamyxx, Gamyxx, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamyxy, Gamyxy, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamyxz, Gamyxz, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamyyy, Gamyyy, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamyyz, Gamyyz, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamyzz, Gamyzz, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamzxx, Gamzxx, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamzxy, Gamzxy, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamzxz, Gamzxz, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamzyy, Gamzyy, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamzyz, Gamzyz, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamzzz, Gamzzz, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Rxx,   Ricxx, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Rxy,   Ricxy, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Rxz,   Ricxz, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Ryy,   Ricyy, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Ryz,   Ricyz, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Rzz,   Riczz, nb, cudaMemcpyHostToDevice));

    // ---- run GPU kernel ----
    // GPU parameter order (from bssn_rhs.h:301-324):
    //   stream, ex, pex, T, X,Y,Z,
    //   chi,trK, gxx..gzz, Axx..Azz, Gamx..Gamz, Lap,betax..dtSfz,
    //   chi_rhs..dtSfz_rhs (23 RHS outputs),
    //   rho..Szz (10 stress-energy),
    //   Gamxxx..Gamzzz (18 Christoffel),
    //   Rxx..Rzz (6 Ricci),
    //   ham_Res..Gmz_Res (7 constraints),
    //   symmetry, lev, eps, co
    double T = 0.0;
    gpu_compute_rhs_bssn_launch(
        stream, ex, pex, T,
        d_X, d_Y, d_Z,
        d_chi, d_trK,
        d_gxx, d_gxy, d_gxz, d_gyy, d_gyz, d_gzz,
        d_Axx, d_Axy, d_Axz, d_Ayy, d_Ayz, d_Azz,
        d_Gamx, d_Gamy, d_Gamz, d_Lap,
        d_betax, d_betay, d_betaz,
        d_dtSfx, d_dtSfy, d_dtSfz,
        d_chi_rhs, d_trK_rhs,
        d_gxx_rhs, d_gxy_rhs, d_gxz_rhs, d_gyy_rhs, d_gyz_rhs, d_gzz_rhs,
        d_Axx_rhs, d_Axy_rhs, d_Axz_rhs, d_Ayy_rhs, d_Ayz_rhs, d_Azz_rhs,
        d_Gamx_rhs, d_Gamy_rhs, d_Gamz_rhs, d_Lap_rhs,
        d_betax_rhs, d_betay_rhs, d_betaz_rhs,
        d_dtSfx_rhs, d_dtSfy_rhs, d_dtSfz_rhs,
        d_rho, d_Sx, d_Sy, d_Sz,
        d_Sxx, d_Sxy, d_Sxz, d_Syy, d_Syz, d_Szz,
        d_Gamxxx, d_Gamxxy, d_Gamxxz, d_Gamxyy, d_Gamxyz, d_Gamxzz,
        d_Gamyxx, d_Gamyxy, d_Gamyxz, d_Gamyyy, d_Gamyyz, d_Gamyzz,
        d_Gamzxx, d_Gamzxy, d_Gamzxz, d_Gamzyy, d_Gamzyz, d_Gamzzz,
        d_Rxx, d_Rxy, d_Rxz, d_Ryy, d_Ryz, d_Rzz,
        d_ham, d_movx, d_movy, d_movz, d_Gmx, d_Gmy, d_Gmz,
        0 /*symmetry*/, 0 /*lev*/, 0.0 /*eps*/, 0 /*co*/);
    CUDA_CHECK(cudaDeviceSynchronize());

    // ---- copy GPU outputs back to host ----
    double *gpu_chi_rhs  = new double[n]; double *gpu_trK_rhs  = new double[n];
    double *gpu_gxx_rhs  = new double[n]; double *gpu_gxy_rhs  = new double[n]; double *gpu_gxz_rhs  = new double[n];
    double *gpu_gyy_rhs  = new double[n]; double *gpu_gyz_rhs  = new double[n]; double *gpu_gzz_rhs  = new double[n];
    double *gpu_Axx_rhs  = new double[n]; double *gpu_Axy_rhs  = new double[n]; double *gpu_Axz_rhs  = new double[n];
    double *gpu_Ayy_rhs  = new double[n]; double *gpu_Ayz_rhs  = new double[n]; double *gpu_Azz_rhs  = new double[n];
    double *gpu_Gamx_rhs = new double[n]; double *gpu_Gamy_rhs = new double[n]; double *gpu_Gamz_rhs = new double[n];
    double *gpu_Lap_rhs  = new double[n];
    double *gpu_betax_rhs = new double[n]; double *gpu_betay_rhs = new double[n]; double *gpu_betaz_rhs = new double[n];
    double *gpu_dtSfx_rhs = new double[n]; double *gpu_dtSfy_rhs = new double[n]; double *gpu_dtSfz_rhs = new double[n];
    double *gpu_ham  = new double[n]; double *gpu_movx = new double[n]; double *gpu_movy = new double[n];
    double *gpu_movz = new double[n]; double *gpu_Gmx  = new double[n]; double *gpu_Gmy  = new double[n]; double *gpu_Gmz  = new double[n];
    // Also copy back stress-energy to verify it's unchanged
    double *gpu_rho = new double[n]; double *gpu_Sx  = new double[n]; double *gpu_Sy  = new double[n]; double *gpu_Sz  = new double[n];
    double *gpu_Sxx = new double[n]; double *gpu_Sxy = new double[n]; double *gpu_Sxz = new double[n];
    double *gpu_Syy = new double[n]; double *gpu_Syz = new double[n]; double *gpu_Szz = new double[n];

    CUDA_CHECK(cudaMemcpy(gpu_chi_rhs,  d_chi_rhs,  nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_trK_rhs,  d_trK_rhs,  nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_gxx_rhs,  d_gxx_rhs,  nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_gxy_rhs,  d_gxy_rhs,  nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_gxz_rhs,  d_gxz_rhs,  nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_gyy_rhs,  d_gyy_rhs,  nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_gyz_rhs,  d_gyz_rhs,  nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_gzz_rhs,  d_gzz_rhs,  nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Axx_rhs,  d_Axx_rhs,  nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Axy_rhs,  d_Axy_rhs,  nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Axz_rhs,  d_Axz_rhs,  nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Ayy_rhs,  d_Ayy_rhs,  nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Ayz_rhs,  d_Ayz_rhs,  nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Azz_rhs,  d_Azz_rhs,  nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Gamx_rhs, d_Gamx_rhs, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Gamy_rhs, d_Gamy_rhs, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Gamz_rhs, d_Gamz_rhs, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Lap_rhs,  d_Lap_rhs,  nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_betax_rhs, d_betax_rhs, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_betay_rhs, d_betay_rhs, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_betaz_rhs, d_betaz_rhs, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_dtSfx_rhs, d_dtSfx_rhs, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_dtSfy_rhs, d_dtSfy_rhs, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_dtSfz_rhs, d_dtSfz_rhs, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_ham,  d_ham,  nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_movx, d_movx, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_movy, d_movy, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_movz, d_movz, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Gmx,  d_Gmx,  nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Gmy,  d_Gmy,  nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Gmz,  d_Gmz,  nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_rho,  d_rho,  nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Sx,   d_Sx,   nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Sy,   d_Sy,   nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Sz,   d_Sz,   nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Sxx,  d_Sxx,  nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Sxy,  d_Sxy,  nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Sxz,  d_Sxz,  nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Syy,  d_Syy,  nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Syz,  d_Syz,  nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Szz,  d_Szz,  nb, cudaMemcpyDeviceToHost));

    // ---- run CPU Fortran function ----
    // CPU parameter order (from bssn_rhs.h:41-58):
    //   ex, T, X, Y, Z,
    //   chi, trK, gij(6), Aij(6), Gam(3), Gauge(7),
    //   chi_rhs, trK_rhs, gij_rhs(6), Aij_rhs(6), Gam_rhs(3), Gauge_rhs(7),
    //   stress-energy(10), Christoffel(18), Ricci(6), constraints(7),
    //   symmetry&, lev&, eps&, co&
    double *cpu_chi_rhs   = new double[n]; double *cpu_trK_rhs   = new double[n];
    double *cpu_gxx_rhs   = new double[n]; double *cpu_gxy_rhs   = new double[n]; double *cpu_gxz_rhs   = new double[n];
    double *cpu_gyy_rhs   = new double[n]; double *cpu_gyz_rhs   = new double[n]; double *cpu_gzz_rhs   = new double[n];
    double *cpu_Axx_rhs   = new double[n]; double *cpu_Axy_rhs   = new double[n]; double *cpu_Axz_rhs   = new double[n];
    double *cpu_Ayy_rhs   = new double[n]; double *cpu_Ayz_rhs   = new double[n]; double *cpu_Azz_rhs   = new double[n];
    double *cpu_Gamx_rhs  = new double[n]; double *cpu_Gamy_rhs  = new double[n]; double *cpu_Gamz_rhs  = new double[n];
    double *cpu_Lap_rhs   = new double[n];
    double *cpu_betax_rhs = new double[n]; double *cpu_betay_rhs = new double[n]; double *cpu_betaz_rhs = new double[n];
    double *cpu_dtSfx_rhs = new double[n]; double *cpu_dtSfy_rhs = new double[n]; double *cpu_dtSfz_rhs = new double[n];
    double *cpu_ham  = new double[n]; double *cpu_movx = new double[n]; double *cpu_movy = new double[n];
    double *cpu_movz = new double[n]; double *cpu_Gmx  = new double[n]; double *cpu_Gmy  = new double[n]; double *cpu_Gmz  = new double[n];

    int    sym_ref = 0;
    int    lev_ref = 0;
    double eps_ref = 0.0;
    int    co_ref  = 0;

    int ret = f_compute_rhs_bssn(
        ex, T,
        X, Y, Z,
        chi, trK,
        gxx, gxy, gxz, gyy, gyz, gzz,
        Axx, Axy, Axz, Ayy, Ayz, Azz,
        Gamx, Gamy, Gamz,
        Lap, betax, betay, betaz, dtSfx, dtSfy, dtSfz,
        cpu_chi_rhs, cpu_trK_rhs,
        cpu_gxx_rhs, cpu_gxy_rhs, cpu_gxz_rhs, cpu_gyy_rhs, cpu_gyz_rhs, cpu_gzz_rhs,
        cpu_Axx_rhs, cpu_Axy_rhs, cpu_Axz_rhs, cpu_Ayy_rhs, cpu_Ayz_rhs, cpu_Azz_rhs,
        cpu_Gamx_rhs, cpu_Gamy_rhs, cpu_Gamz_rhs,
        cpu_Lap_rhs, cpu_betax_rhs, cpu_betay_rhs, cpu_betaz_rhs,
        cpu_dtSfx_rhs, cpu_dtSfy_rhs, cpu_dtSfz_rhs,
        rho, Sx, Sy, Sz,
        Sxx, Sxy, Sxz, Syy, Syz, Szz,
        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
        Ricxx, Ricxy, Ricxz, Ricyy, Ricyz, Riczz,
        cpu_ham, cpu_movx, cpu_movy, cpu_movz, cpu_Gmx, cpu_Gmy, cpu_Gmz,
        sym_ref, lev_ref, eps_ref, co_ref);

    printf("  CPU f_compute_rhs_bssn ret=%d\n", ret);

    // ---- compare RHS outputs (interior only) ----
    int fails = 0;
    fails += compare_field("chi_rhs",    gpu_chi_rhs,    cpu_chi_rhs,    ex, TOL);
    fails += compare_field("trK_rhs",    gpu_trK_rhs,    cpu_trK_rhs,    ex, TOL);
    fails += compare_field("gxx_rhs",    gpu_gxx_rhs,    cpu_gxx_rhs,    ex, TOL);
    fails += compare_field("gxy_rhs",    gpu_gxy_rhs,    cpu_gxy_rhs,    ex, TOL);
    fails += compare_field("gxz_rhs",    gpu_gxz_rhs,    cpu_gxz_rhs,    ex, TOL);
    fails += compare_field("gyy_rhs",    gpu_gyy_rhs,    cpu_gyy_rhs,    ex, TOL);
    fails += compare_field("gyz_rhs",    gpu_gyz_rhs,    cpu_gyz_rhs,    ex, TOL);
    fails += compare_field("gzz_rhs",    gpu_gzz_rhs,    cpu_gzz_rhs,    ex, TOL);
    fails += compare_field("Axx_rhs",    gpu_Axx_rhs,    cpu_Axx_rhs,    ex, TOL);
    fails += compare_field("Axy_rhs",    gpu_Axy_rhs,    cpu_Axy_rhs,    ex, TOL);
    fails += compare_field("Axz_rhs",    gpu_Axz_rhs,    cpu_Axz_rhs,    ex, TOL);
    fails += compare_field("Ayy_rhs",    gpu_Ayy_rhs,    cpu_Ayy_rhs,    ex, TOL);
    fails += compare_field("Ayz_rhs",    gpu_Ayz_rhs,    cpu_Ayz_rhs,    ex, TOL);
    fails += compare_field("Azz_rhs",    gpu_Azz_rhs,    cpu_Azz_rhs,    ex, TOL);
    fails += compare_field("Gamx_rhs",   gpu_Gamx_rhs,   cpu_Gamx_rhs,   ex, TOL);
    fails += compare_field("Gamy_rhs",   gpu_Gamy_rhs,   cpu_Gamy_rhs,   ex, TOL);
    fails += compare_field("Gamz_rhs",   gpu_Gamz_rhs,   cpu_Gamz_rhs,   ex, TOL);
    fails += compare_field("Lap_rhs",    gpu_Lap_rhs,    cpu_Lap_rhs,    ex, TOL);
    fails += compare_field("betax_rhs",  gpu_betax_rhs,  cpu_betax_rhs,  ex, TOL);
    fails += compare_field("betay_rhs",  gpu_betay_rhs,  cpu_betay_rhs,  ex, TOL);
    fails += compare_field("betaz_rhs",  gpu_betaz_rhs,  cpu_betaz_rhs,  ex, TOL);
    fails += compare_field("dtSfx_rhs",  gpu_dtSfx_rhs,  cpu_dtSfx_rhs,  ex, TOL);
    fails += compare_field("dtSfy_rhs",  gpu_dtSfy_rhs,  cpu_dtSfy_rhs,  ex, TOL);
    fails += compare_field("dtSfz_rhs",  gpu_dtSfz_rhs,  cpu_dtSfz_rhs,  ex, TOL);
    fails += compare_field("ham_Res",    gpu_ham,        cpu_ham,        ex, TOL);
    fails += compare_field("movx_Res",   gpu_movx,       cpu_movx,       ex, TOL);
    fails += compare_field("movy_Res",   gpu_movy,       cpu_movy,       ex, TOL);
    fails += compare_field("movz_Res",   gpu_movz,       cpu_movz,       ex, TOL);
    fails += compare_field("Gmx_Res",    gpu_Gmx,        cpu_Gmx,        ex, TOL);
    fails += compare_field("Gmy_Res",    gpu_Gmy,        cpu_Gmy,        ex, TOL);
    fails += compare_field("Gmz_Res",    gpu_Gmz,        cpu_Gmz,        ex, TOL);

    // ---- verify stress-energy is unchanged (BSSN RHS reads them, doesn't modify) ----
    printf("  --- stress-energy unchanged check ---\n");
    fails += compare_field("rho(unch)",  gpu_rho, rho, ex, TOL);
    fails += compare_field("Sx(unch)",   gpu_Sx,  Sx,  ex, TOL);
    fails += compare_field("Sy(unch)",   gpu_Sy,  Sy,  ex, TOL);
    fails += compare_field("Sz(unch)",   gpu_Sz,  Sz,  ex, TOL);
    fails += compare_field("Sxx(unch)",  gpu_Sxx, Sxx, ex, TOL);
    fails += compare_field("Sxy(unch)",  gpu_Sxy, Sxy, ex, TOL);
    fails += compare_field("Sxz(unch)",  gpu_Sxz, Sxz, ex, TOL);
    fails += compare_field("Syy(unch)",  gpu_Syy, Syy, ex, TOL);
    fails += compare_field("Syz(unch)",  gpu_Syz, Syz, ex, TOL);
    fails += compare_field("Szz(unch)",  gpu_Szz, Szz, ex, TOL);

    // ---- cleanup host ----
    delete[] X;    delete[] Y;    delete[] Z;
    delete[] chi;  delete[] trK;
    delete[] gxx;  delete[] gxy;  delete[] gxz;
    delete[] gyy;  delete[] gyz;  delete[] gzz;
    delete[] Axx;  delete[] Axy;  delete[] Axz;
    delete[] Ayy;  delete[] Ayz;  delete[] Azz;
    delete[] Gamx; delete[] Gamy; delete[] Gamz;
    delete[] Lap;  delete[] betax; delete[] betay; delete[] betaz;
    delete[] dtSfx; delete[] dtSfy; delete[] dtSfz;
    delete[] rho; delete[] Sx;  delete[] Sy;  delete[] Sz;
    delete[] Sxx; delete[] Sxy; delete[] Sxz;
    delete[] Syy; delete[] Syz; delete[] Szz;
    delete[] Gamxxx; delete[] Gamxxy; delete[] Gamxxz;
    delete[] Gamxyy; delete[] Gamxyz; delete[] Gamxzz;
    delete[] Gamyxx; delete[] Gamyxy; delete[] Gamyxz;
    delete[] Gamyyy; delete[] Gamyyz; delete[] Gamyzz;
    delete[] Gamzxx; delete[] Gamzxy; delete[] Gamzxz;
    delete[] Gamzyy; delete[] Gamzyz; delete[] Gamzzz;
    delete[] Ricxx; delete[] Ricxy; delete[] Ricxz;
    delete[] Ricyy; delete[] Ricyz; delete[] Riczz;

    delete[] gpu_chi_rhs;  delete[] gpu_trK_rhs;
    delete[] gpu_gxx_rhs;  delete[] gpu_gxy_rhs;  delete[] gpu_gxz_rhs;
    delete[] gpu_gyy_rhs;  delete[] gpu_gyz_rhs;  delete[] gpu_gzz_rhs;
    delete[] gpu_Axx_rhs;  delete[] gpu_Axy_rhs;  delete[] gpu_Axz_rhs;
    delete[] gpu_Ayy_rhs;  delete[] gpu_Ayz_rhs;  delete[] gpu_Azz_rhs;
    delete[] gpu_Gamx_rhs; delete[] gpu_Gamy_rhs; delete[] gpu_Gamz_rhs;
    delete[] gpu_Lap_rhs;
    delete[] gpu_betax_rhs; delete[] gpu_betay_rhs; delete[] gpu_betaz_rhs;
    delete[] gpu_dtSfx_rhs; delete[] gpu_dtSfy_rhs; delete[] gpu_dtSfz_rhs;
    delete[] gpu_ham;  delete[] gpu_movx; delete[] gpu_movy;
    delete[] gpu_movz; delete[] gpu_Gmx;  delete[] gpu_Gmy; delete[] gpu_Gmz;
    delete[] gpu_rho; delete[] gpu_Sx;  delete[] gpu_Sy;  delete[] gpu_Sz;
    delete[] gpu_Sxx; delete[] gpu_Sxy; delete[] gpu_Sxz;
    delete[] gpu_Syy; delete[] gpu_Syz; delete[] gpu_Szz;

    delete[] cpu_chi_rhs;  delete[] cpu_trK_rhs;
    delete[] cpu_gxx_rhs;  delete[] cpu_gxy_rhs;  delete[] cpu_gxz_rhs;
    delete[] cpu_gyy_rhs;  delete[] cpu_gyz_rhs;  delete[] cpu_gzz_rhs;
    delete[] cpu_Axx_rhs;  delete[] cpu_Axy_rhs;  delete[] cpu_Axz_rhs;
    delete[] cpu_Ayy_rhs;  delete[] cpu_Ayz_rhs;  delete[] cpu_Azz_rhs;
    delete[] cpu_Gamx_rhs; delete[] cpu_Gamy_rhs; delete[] cpu_Gamz_rhs;
    delete[] cpu_Lap_rhs;
    delete[] cpu_betax_rhs; delete[] cpu_betay_rhs; delete[] cpu_betaz_rhs;
    delete[] cpu_dtSfx_rhs; delete[] cpu_dtSfy_rhs; delete[] cpu_dtSfz_rhs;
    delete[] cpu_ham;  delete[] cpu_movx; delete[] cpu_movy;
    delete[] cpu_movz; delete[] cpu_Gmx;  delete[] cpu_Gmy; delete[] cpu_Gmz;

    // ---- cleanup GPU ----
    CUDA_CHECK(cudaFree(d_X));    CUDA_CHECK(cudaFree(d_Y));    CUDA_CHECK(cudaFree(d_Z));
    CUDA_CHECK(cudaFree(d_chi));  CUDA_CHECK(cudaFree(d_trK));
    CUDA_CHECK(cudaFree(d_gxx));  CUDA_CHECK(cudaFree(d_gxy));  CUDA_CHECK(cudaFree(d_gxz));
    CUDA_CHECK(cudaFree(d_gyy));  CUDA_CHECK(cudaFree(d_gyz));  CUDA_CHECK(cudaFree(d_gzz));
    CUDA_CHECK(cudaFree(d_Axx));  CUDA_CHECK(cudaFree(d_Axy));  CUDA_CHECK(cudaFree(d_Axz));
    CUDA_CHECK(cudaFree(d_Ayy));  CUDA_CHECK(cudaFree(d_Ayz));  CUDA_CHECK(cudaFree(d_Azz));
    CUDA_CHECK(cudaFree(d_Gamx)); CUDA_CHECK(cudaFree(d_Gamy)); CUDA_CHECK(cudaFree(d_Gamz));
    CUDA_CHECK(cudaFree(d_Lap));
    CUDA_CHECK(cudaFree(d_betax)); CUDA_CHECK(cudaFree(d_betay)); CUDA_CHECK(cudaFree(d_betaz));
    CUDA_CHECK(cudaFree(d_dtSfx)); CUDA_CHECK(cudaFree(d_dtSfy)); CUDA_CHECK(cudaFree(d_dtSfz));
    CUDA_CHECK(cudaFree(d_chi_rhs));  CUDA_CHECK(cudaFree(d_trK_rhs));
    CUDA_CHECK(cudaFree(d_gxx_rhs));  CUDA_CHECK(cudaFree(d_gxy_rhs));  CUDA_CHECK(cudaFree(d_gxz_rhs));
    CUDA_CHECK(cudaFree(d_gyy_rhs));  CUDA_CHECK(cudaFree(d_gyz_rhs));  CUDA_CHECK(cudaFree(d_gzz_rhs));
    CUDA_CHECK(cudaFree(d_Axx_rhs));  CUDA_CHECK(cudaFree(d_Axy_rhs));  CUDA_CHECK(cudaFree(d_Axz_rhs));
    CUDA_CHECK(cudaFree(d_Ayy_rhs));  CUDA_CHECK(cudaFree(d_Ayz_rhs));  CUDA_CHECK(cudaFree(d_Azz_rhs));
    CUDA_CHECK(cudaFree(d_Gamx_rhs)); CUDA_CHECK(cudaFree(d_Gamy_rhs)); CUDA_CHECK(cudaFree(d_Gamz_rhs));
    CUDA_CHECK(cudaFree(d_Lap_rhs));
    CUDA_CHECK(cudaFree(d_betax_rhs)); CUDA_CHECK(cudaFree(d_betay_rhs)); CUDA_CHECK(cudaFree(d_betaz_rhs));
    CUDA_CHECK(cudaFree(d_dtSfx_rhs)); CUDA_CHECK(cudaFree(d_dtSfy_rhs)); CUDA_CHECK(cudaFree(d_dtSfz_rhs));
    CUDA_CHECK(cudaFree(d_rho));   CUDA_CHECK(cudaFree(d_Sx));  CUDA_CHECK(cudaFree(d_Sy));  CUDA_CHECK(cudaFree(d_Sz));
    CUDA_CHECK(cudaFree(d_Sxx));   CUDA_CHECK(cudaFree(d_Sxy)); CUDA_CHECK(cudaFree(d_Sxz));
    CUDA_CHECK(cudaFree(d_Syy));   CUDA_CHECK(cudaFree(d_Syz)); CUDA_CHECK(cudaFree(d_Szz));
    CUDA_CHECK(cudaFree(d_Gamxxx)); CUDA_CHECK(cudaFree(d_Gamxxy)); CUDA_CHECK(cudaFree(d_Gamxxz));
    CUDA_CHECK(cudaFree(d_Gamxyy)); CUDA_CHECK(cudaFree(d_Gamxyz)); CUDA_CHECK(cudaFree(d_Gamxzz));
    CUDA_CHECK(cudaFree(d_Gamyxx)); CUDA_CHECK(cudaFree(d_Gamyxy)); CUDA_CHECK(cudaFree(d_Gamyxz));
    CUDA_CHECK(cudaFree(d_Gamyyy)); CUDA_CHECK(cudaFree(d_Gamyyz)); CUDA_CHECK(cudaFree(d_Gamyzz));
    CUDA_CHECK(cudaFree(d_Gamzxx)); CUDA_CHECK(cudaFree(d_Gamzxy)); CUDA_CHECK(cudaFree(d_Gamzxz));
    CUDA_CHECK(cudaFree(d_Gamzyy)); CUDA_CHECK(cudaFree(d_Gamzyz)); CUDA_CHECK(cudaFree(d_Gamzzz));
    CUDA_CHECK(cudaFree(d_Rxx));   CUDA_CHECK(cudaFree(d_Rxy));   CUDA_CHECK(cudaFree(d_Rxz));
    CUDA_CHECK(cudaFree(d_Ryy));   CUDA_CHECK(cudaFree(d_Ryz));   CUDA_CHECK(cudaFree(d_Rzz));
    CUDA_CHECK(cudaFree(d_ham));   CUDA_CHECK(cudaFree(d_movx));
    CUDA_CHECK(cudaFree(d_movy));  CUDA_CHECK(cudaFree(d_movz));
    CUDA_CHECK(cudaFree(d_Gmx));   CUDA_CHECK(cudaFree(d_Gmy));   CUDA_CHECK(cudaFree(d_Gmz));

    return fails;
}

// ===========================================================================
// Test: BSSN RHS with puncture-like data (chi -> -1 at center)
// Runs same GPU-vs-CPU comparison as test_bssn_rhs() but with Brill-Lindquist
// puncture initial data for various puncture masses M.
// ===========================================================================
static int test_bssn_rhs_puncture()
{
    double M_values[] = {0.01, 0.1, 0.5, 1.0};
    const int nM = 4;
    double h = 0.1;
    int total_fails = 0;

    int ex[3]  = {8 + 2*GW, 8 + 2*GW, 8 + 2*GW};
    int pex[3] = {ex[0], ex[1], ex[2]};
    int n = ex[0] * ex[1] * ex[2];
    size_t nb = n * sizeof(double);

    // ---- allocate host arrays (reused for all M) ----
    double *X    = new double[n];
    double *Y    = new double[n];
    double *Z    = new double[n];

    double *chi   = new double[n];
    double *trK   = new double[n];
    double *gxx   = new double[n]; double *gxy   = new double[n]; double *gxz   = new double[n];
    double *gyy   = new double[n]; double *gyz   = new double[n]; double *gzz   = new double[n];
    double *Axx   = new double[n]; double *Axy   = new double[n]; double *Axz   = new double[n];
    double *Ayy   = new double[n]; double *Ayz   = new double[n]; double *Azz   = new double[n];
    double *Gamx  = new double[n]; double *Gamy  = new double[n]; double *Gamz  = new double[n];
    double *Lap   = new double[n];
    double *betax = new double[n]; double *betay = new double[n]; double *betaz = new double[n];
    double *dtSfx = new double[n]; double *dtSfy = new double[n]; double *dtSfz = new double[n];

    double *rho = new double[n];
    double *Sx  = new double[n]; double *Sy  = new double[n]; double *Sz  = new double[n];
    double *Sxx = new double[n]; double *Sxy = new double[n]; double *Sxz = new double[n];
    double *Syy = new double[n]; double *Syz = new double[n]; double *Szz = new double[n];

    double *Gamxxx = new double[n]; double *Gamxxy = new double[n]; double *Gamxxz = new double[n];
    double *Gamxyy = new double[n]; double *Gamxyz = new double[n]; double *Gamxzz = new double[n];
    double *Gamyxx = new double[n]; double *Gamyxy = new double[n]; double *Gamyxz = new double[n];
    double *Gamyyy = new double[n]; double *Gamyyz = new double[n]; double *Gamyzz = new double[n];
    double *Gamzxx = new double[n]; double *Gamzxy = new double[n]; double *Gamzxz = new double[n];
    double *Gamzyy = new double[n]; double *Gamzyz = new double[n]; double *Gamzzz = new double[n];

    double *Ricxx = new double[n]; double *Ricxy = new double[n]; double *Ricxz = new double[n];
    double *Ricyy = new double[n]; double *Ricyz = new double[n]; double *Riczz = new double[n];

    // Fill stress-energy, Christoffel, Ricci once (independent of M)
    fill_stress_energy_test_data(ex, rho, Sx, Sy, Sz,
                                  Sxx, Sxy, Sxz, Syy, Syz, Szz);
    fill_christoffel_test_data(ex,
                                Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
                                Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
                                Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz);
    fill_ricci_test_data(ex, Ricxx, Ricxy, Ricxz, Ricyy, Ricyz, Riczz);

    // ---- allocate GPU memory (reused for all M) ----
    double *d_X, *d_Y, *d_Z;
    double *d_chi, *d_trK;
    double *d_gxx, *d_gxy, *d_gxz, *d_gyy, *d_gyz, *d_gzz;
    double *d_Axx, *d_Axy, *d_Axz, *d_Ayy, *d_Ayz, *d_Azz;
    double *d_Gamx, *d_Gamy, *d_Gamz, *d_Lap;
    double *d_betax, *d_betay, *d_betaz;
    double *d_dtSfx, *d_dtSfy, *d_dtSfz;
    double *d_chi_rhs, *d_trK_rhs;
    double *d_gxx_rhs, *d_gxy_rhs, *d_gxz_rhs, *d_gyy_rhs, *d_gyz_rhs, *d_gzz_rhs;
    double *d_Axx_rhs, *d_Axy_rhs, *d_Axz_rhs, *d_Ayy_rhs, *d_Ayz_rhs, *d_Azz_rhs;
    double *d_Gamx_rhs, *d_Gamy_rhs, *d_Gamz_rhs, *d_Lap_rhs;
    double *d_betax_rhs, *d_betay_rhs, *d_betaz_rhs;
    double *d_dtSfx_rhs, *d_dtSfy_rhs, *d_dtSfz_rhs;
    double *d_rho, *d_Sx, *d_Sy, *d_Sz;
    double *d_Sxx, *d_Sxy, *d_Sxz, *d_Syy, *d_Syz, *d_Szz;
    double *d_Gamxxx, *d_Gamxxy, *d_Gamxxz, *d_Gamxyy, *d_Gamxyz, *d_Gamxzz;
    double *d_Gamyxx, *d_Gamyxy, *d_Gamyxz, *d_Gamyyy, *d_Gamyyz, *d_Gamyzz;
    double *d_Gamzxx, *d_Gamzxy, *d_Gamzxz, *d_Gamzyy, *d_Gamzyz, *d_Gamzzz;
    double *d_Rxx, *d_Rxy, *d_Rxz, *d_Ryy, *d_Ryz, *d_Rzz;
    double *d_ham, *d_movx, *d_movy, *d_movz, *d_Gmx, *d_Gmy, *d_Gmz;

    CUDA_CHECK(cudaMalloc(&d_X,    nb)); CUDA_CHECK(cudaMalloc(&d_Y,    nb)); CUDA_CHECK(cudaMalloc(&d_Z,    nb));
    CUDA_CHECK(cudaMalloc(&d_chi,  nb)); CUDA_CHECK(cudaMalloc(&d_trK,  nb));
    CUDA_CHECK(cudaMalloc(&d_gxx,  nb)); CUDA_CHECK(cudaMalloc(&d_gxy,  nb)); CUDA_CHECK(cudaMalloc(&d_gxz,  nb));
    CUDA_CHECK(cudaMalloc(&d_gyy,  nb)); CUDA_CHECK(cudaMalloc(&d_gyz,  nb)); CUDA_CHECK(cudaMalloc(&d_gzz,  nb));
    CUDA_CHECK(cudaMalloc(&d_Axx,  nb)); CUDA_CHECK(cudaMalloc(&d_Axy,  nb)); CUDA_CHECK(cudaMalloc(&d_Axz,  nb));
    CUDA_CHECK(cudaMalloc(&d_Ayy,  nb)); CUDA_CHECK(cudaMalloc(&d_Ayz,  nb)); CUDA_CHECK(cudaMalloc(&d_Azz,  nb));
    CUDA_CHECK(cudaMalloc(&d_Gamx, nb)); CUDA_CHECK(cudaMalloc(&d_Gamy, nb)); CUDA_CHECK(cudaMalloc(&d_Gamz, nb));
    CUDA_CHECK(cudaMalloc(&d_Lap,  nb));
    CUDA_CHECK(cudaMalloc(&d_betax, nb)); CUDA_CHECK(cudaMalloc(&d_betay, nb)); CUDA_CHECK(cudaMalloc(&d_betaz, nb));
    CUDA_CHECK(cudaMalloc(&d_dtSfx, nb)); CUDA_CHECK(cudaMalloc(&d_dtSfy, nb)); CUDA_CHECK(cudaMalloc(&d_dtSfz, nb));
    CUDA_CHECK(cudaMalloc(&d_chi_rhs,  nb)); CUDA_CHECK(cudaMalloc(&d_trK_rhs,  nb));
    CUDA_CHECK(cudaMalloc(&d_gxx_rhs,  nb)); CUDA_CHECK(cudaMalloc(&d_gxy_rhs,  nb)); CUDA_CHECK(cudaMalloc(&d_gxz_rhs,  nb));
    CUDA_CHECK(cudaMalloc(&d_gyy_rhs,  nb)); CUDA_CHECK(cudaMalloc(&d_gyz_rhs,  nb)); CUDA_CHECK(cudaMalloc(&d_gzz_rhs,  nb));
    CUDA_CHECK(cudaMalloc(&d_Axx_rhs,  nb)); CUDA_CHECK(cudaMalloc(&d_Axy_rhs,  nb)); CUDA_CHECK(cudaMalloc(&d_Axz_rhs,  nb));
    CUDA_CHECK(cudaMalloc(&d_Ayy_rhs,  nb)); CUDA_CHECK(cudaMalloc(&d_Ayz_rhs,  nb)); CUDA_CHECK(cudaMalloc(&d_Azz_rhs,  nb));
    CUDA_CHECK(cudaMalloc(&d_Gamx_rhs, nb)); CUDA_CHECK(cudaMalloc(&d_Gamy_rhs, nb)); CUDA_CHECK(cudaMalloc(&d_Gamz_rhs, nb));
    CUDA_CHECK(cudaMalloc(&d_Lap_rhs,  nb));
    CUDA_CHECK(cudaMalloc(&d_betax_rhs, nb)); CUDA_CHECK(cudaMalloc(&d_betay_rhs, nb)); CUDA_CHECK(cudaMalloc(&d_betaz_rhs, nb));
    CUDA_CHECK(cudaMalloc(&d_dtSfx_rhs, nb)); CUDA_CHECK(cudaMalloc(&d_dtSfy_rhs, nb)); CUDA_CHECK(cudaMalloc(&d_dtSfz_rhs, nb));
    CUDA_CHECK(cudaMalloc(&d_rho,   nb)); CUDA_CHECK(cudaMalloc(&d_Sx,  nb)); CUDA_CHECK(cudaMalloc(&d_Sy,  nb)); CUDA_CHECK(cudaMalloc(&d_Sz, nb));
    CUDA_CHECK(cudaMalloc(&d_Sxx,   nb)); CUDA_CHECK(cudaMalloc(&d_Sxy, nb)); CUDA_CHECK(cudaMalloc(&d_Sxz, nb));
    CUDA_CHECK(cudaMalloc(&d_Syy,   nb)); CUDA_CHECK(cudaMalloc(&d_Syz, nb)); CUDA_CHECK(cudaMalloc(&d_Szz, nb));
    CUDA_CHECK(cudaMalloc(&d_Gamxxx, nb)); CUDA_CHECK(cudaMalloc(&d_Gamxxy, nb)); CUDA_CHECK(cudaMalloc(&d_Gamxxz, nb));
    CUDA_CHECK(cudaMalloc(&d_Gamxyy, nb)); CUDA_CHECK(cudaMalloc(&d_Gamxyz, nb)); CUDA_CHECK(cudaMalloc(&d_Gamxzz, nb));
    CUDA_CHECK(cudaMalloc(&d_Gamyxx, nb)); CUDA_CHECK(cudaMalloc(&d_Gamyxy, nb)); CUDA_CHECK(cudaMalloc(&d_Gamyxz, nb));
    CUDA_CHECK(cudaMalloc(&d_Gamyyy, nb)); CUDA_CHECK(cudaMalloc(&d_Gamyyz, nb)); CUDA_CHECK(cudaMalloc(&d_Gamyzz, nb));
    CUDA_CHECK(cudaMalloc(&d_Gamzxx, nb)); CUDA_CHECK(cudaMalloc(&d_Gamzxy, nb)); CUDA_CHECK(cudaMalloc(&d_Gamzxz, nb));
    CUDA_CHECK(cudaMalloc(&d_Gamzyy, nb)); CUDA_CHECK(cudaMalloc(&d_Gamzyz, nb)); CUDA_CHECK(cudaMalloc(&d_Gamzzz, nb));
    CUDA_CHECK(cudaMalloc(&d_Rxx,   nb)); CUDA_CHECK(cudaMalloc(&d_Rxy,   nb)); CUDA_CHECK(cudaMalloc(&d_Rxz,   nb));
    CUDA_CHECK(cudaMalloc(&d_Ryy,   nb)); CUDA_CHECK(cudaMalloc(&d_Ryz,   nb)); CUDA_CHECK(cudaMalloc(&d_Rzz,   nb));
    CUDA_CHECK(cudaMalloc(&d_ham,   nb)); CUDA_CHECK(cudaMalloc(&d_movx,  nb));
    CUDA_CHECK(cudaMalloc(&d_movy,  nb)); CUDA_CHECK(cudaMalloc(&d_movz,  nb));
    CUDA_CHECK(cudaMalloc(&d_Gmx,   nb)); CUDA_CHECK(cudaMalloc(&d_Gmy,   nb)); CUDA_CHECK(cudaMalloc(&d_Gmz,  nb));

    // Copy stress-energy/Christoffel/Ricci to GPU once (same for all M)
    cudaStream_t stream = 0;
    CUDA_CHECK(cudaMemcpy(d_rho,   rho,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Sx,    Sx,    nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Sy,    Sy,    nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Sz,    Sz,    nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Sxx,   Sxx,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Sxy,   Sxy,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Sxz,   Sxz,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Syy,   Syy,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Syz,   Syz,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Szz,   Szz,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamxxx, Gamxxx, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamxxy, Gamxxy, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamxxz, Gamxxz, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamxyy, Gamxyy, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamxyz, Gamxyz, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamxzz, Gamxzz, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamyxx, Gamyxx, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamyxy, Gamyxy, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamyxz, Gamyxz, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamyyy, Gamyyy, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamyyz, Gamyyz, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamyzz, Gamyzz, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamzxx, Gamzxx, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamzxy, Gamzxy, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamzxz, Gamzxz, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamzyy, Gamzyy, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamzyz, Gamzyz, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamzzz, Gamzzz, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Rxx,   Ricxx, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Rxy,   Ricxy, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Rxz,   Ricxz, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Ryy,   Ricyy, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Ryz,   Ricyz, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Rzz,   Riczz, nb, cudaMemcpyHostToDevice));

    // ---- test loop over M values ----
    for (int m_idx = 0; m_idx < nM; m_idx++) {
        double M = M_values[m_idx];
        printf("[BSSN RHS puncture test M=%.2g]\n", M);

        // Fill BSSN state with puncture data for this M
        fill_puncture_test_data(ex, M, h, X, Y, Z,
                                chi, trK,
                                gxx, gxy, gxz, gyy, gyz, gzz,
                                Axx, Axy, Axz, Ayy, Ayz, Azz,
                                Gamx, Gamy, Gamz,
                                Lap, betax, betay, betaz,
                                dtSfx, dtSfy, dtSfz);

        // ---- copy BSSN inputs to GPU ----
        CUDA_CHECK(cudaMemcpy(d_X,    X,    nb, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_Y,    Y,    nb, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_Z,    Z,    nb, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_chi,  chi,  nb, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_trK,  trK,  nb, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_gxx,  gxx,  nb, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_gxy,  gxy,  nb, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_gxz,  gxz,  nb, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_gyy,  gyy,  nb, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_gyz,  gyz,  nb, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_gzz,  gzz,  nb, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_Axx,  Axx,  nb, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_Axy,  Axy,  nb, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_Axz,  Axz,  nb, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_Ayy,  Ayy,  nb, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_Ayz,  Ayz,  nb, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_Azz,  Azz,  nb, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_Gamx, Gamx, nb, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_Gamy, Gamy, nb, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_Gamz, Gamz, nb, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_Lap,  Lap,  nb, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_betax, betax, nb, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_betay, betay, nb, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_betaz, betaz, nb, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_dtSfx, dtSfx, nb, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_dtSfy, dtSfy, nb, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_dtSfz, dtSfz, nb, cudaMemcpyHostToDevice));

        // ---- run GPU kernel ----
        double T = 0.0;
        gpu_compute_rhs_bssn_launch(
            stream, ex, pex, T,
            d_X, d_Y, d_Z,
            d_chi, d_trK,
            d_gxx, d_gxy, d_gxz, d_gyy, d_gyz, d_gzz,
            d_Axx, d_Axy, d_Axz, d_Ayy, d_Ayz, d_Azz,
            d_Gamx, d_Gamy, d_Gamz, d_Lap,
            d_betax, d_betay, d_betaz,
            d_dtSfx, d_dtSfy, d_dtSfz,
            d_chi_rhs, d_trK_rhs,
            d_gxx_rhs, d_gxy_rhs, d_gxz_rhs, d_gyy_rhs, d_gyz_rhs, d_gzz_rhs,
            d_Axx_rhs, d_Axy_rhs, d_Axz_rhs, d_Ayy_rhs, d_Ayz_rhs, d_Azz_rhs,
            d_Gamx_rhs, d_Gamy_rhs, d_Gamz_rhs, d_Lap_rhs,
            d_betax_rhs, d_betay_rhs, d_betaz_rhs,
            d_dtSfx_rhs, d_dtSfy_rhs, d_dtSfz_rhs,
            d_rho, d_Sx, d_Sy, d_Sz,
            d_Sxx, d_Sxy, d_Sxz, d_Syy, d_Syz, d_Szz,
            d_Gamxxx, d_Gamxxy, d_Gamxxz, d_Gamxyy, d_Gamxyz, d_Gamxzz,
            d_Gamyxx, d_Gamyxy, d_Gamyxz, d_Gamyyy, d_Gamyyz, d_Gamyzz,
            d_Gamzxx, d_Gamzxy, d_Gamzxz, d_Gamzyy, d_Gamzyz, d_Gamzzz,
            d_Rxx, d_Rxy, d_Rxz, d_Ryy, d_Ryz, d_Rzz,
            d_ham, d_movx, d_movy, d_movz, d_Gmx, d_Gmy, d_Gmz,
            0 /*symmetry*/, 0 /*lev*/, 0.0 /*eps*/, 0 /*co*/);
        CUDA_CHECK(cudaDeviceSynchronize());

        // ---- copy GPU outputs back to host ----
        double *gpu_chi_rhs  = new double[n]; double *gpu_trK_rhs  = new double[n];
        double *gpu_gxx_rhs  = new double[n]; double *gpu_gxy_rhs  = new double[n]; double *gpu_gxz_rhs  = new double[n];
        double *gpu_gyy_rhs  = new double[n]; double *gpu_gyz_rhs  = new double[n]; double *gpu_gzz_rhs  = new double[n];
        double *gpu_Axx_rhs  = new double[n]; double *gpu_Axy_rhs  = new double[n]; double *gpu_Axz_rhs  = new double[n];
        double *gpu_Ayy_rhs  = new double[n]; double *gpu_Ayz_rhs  = new double[n]; double *gpu_Azz_rhs  = new double[n];
        double *gpu_Gamx_rhs = new double[n]; double *gpu_Gamy_rhs = new double[n]; double *gpu_Gamz_rhs = new double[n];
        double *gpu_Lap_rhs  = new double[n];
        double *gpu_betax_rhs = new double[n]; double *gpu_betay_rhs = new double[n]; double *gpu_betaz_rhs = new double[n];
        double *gpu_dtSfx_rhs = new double[n]; double *gpu_dtSfy_rhs = new double[n]; double *gpu_dtSfz_rhs = new double[n];
        double *gpu_ham  = new double[n]; double *gpu_movx = new double[n]; double *gpu_movy = new double[n];
        double *gpu_movz = new double[n]; double *gpu_Gmx  = new double[n]; double *gpu_Gmy  = new double[n]; double *gpu_Gmz  = new double[n];

        CUDA_CHECK(cudaMemcpy(gpu_chi_rhs,  d_chi_rhs,  nb, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(gpu_trK_rhs,  d_trK_rhs,  nb, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(gpu_gxx_rhs,  d_gxx_rhs,  nb, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(gpu_gxy_rhs,  d_gxy_rhs,  nb, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(gpu_gxz_rhs,  d_gxz_rhs,  nb, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(gpu_gyy_rhs,  d_gyy_rhs,  nb, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(gpu_gyz_rhs,  d_gyz_rhs,  nb, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(gpu_gzz_rhs,  d_gzz_rhs,  nb, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(gpu_Axx_rhs,  d_Axx_rhs,  nb, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(gpu_Axy_rhs,  d_Axy_rhs,  nb, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(gpu_Axz_rhs,  d_Axz_rhs,  nb, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(gpu_Ayy_rhs,  d_Ayy_rhs,  nb, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(gpu_Ayz_rhs,  d_Ayz_rhs,  nb, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(gpu_Azz_rhs,  d_Azz_rhs,  nb, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(gpu_Gamx_rhs, d_Gamx_rhs, nb, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(gpu_Gamy_rhs, d_Gamy_rhs, nb, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(gpu_Gamz_rhs, d_Gamz_rhs, nb, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(gpu_Lap_rhs,  d_Lap_rhs,  nb, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(gpu_betax_rhs, d_betax_rhs, nb, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(gpu_betay_rhs, d_betay_rhs, nb, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(gpu_betaz_rhs, d_betaz_rhs, nb, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(gpu_dtSfx_rhs, d_dtSfx_rhs, nb, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(gpu_dtSfy_rhs, d_dtSfy_rhs, nb, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(gpu_dtSfz_rhs, d_dtSfz_rhs, nb, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(gpu_ham,  d_ham,  nb, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(gpu_movx, d_movx, nb, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(gpu_movy, d_movy, nb, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(gpu_movz, d_movz, nb, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(gpu_Gmx,  d_Gmx,  nb, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(gpu_Gmy,  d_Gmy,  nb, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(gpu_Gmz,  d_Gmz,  nb, cudaMemcpyDeviceToHost));

        // ---- run CPU Fortran function ----
        double *cpu_chi_rhs   = new double[n]; double *cpu_trK_rhs   = new double[n];
        double *cpu_gxx_rhs   = new double[n]; double *cpu_gxy_rhs   = new double[n]; double *cpu_gxz_rhs   = new double[n];
        double *cpu_gyy_rhs   = new double[n]; double *cpu_gyz_rhs   = new double[n]; double *cpu_gzz_rhs   = new double[n];
        double *cpu_Axx_rhs   = new double[n]; double *cpu_Axy_rhs   = new double[n]; double *cpu_Axz_rhs   = new double[n];
        double *cpu_Ayy_rhs   = new double[n]; double *cpu_Ayz_rhs   = new double[n]; double *cpu_Azz_rhs   = new double[n];
        double *cpu_Gamx_rhs  = new double[n]; double *cpu_Gamy_rhs  = new double[n]; double *cpu_Gamz_rhs  = new double[n];
        double *cpu_Lap_rhs   = new double[n];
        double *cpu_betax_rhs = new double[n]; double *cpu_betay_rhs = new double[n]; double *cpu_betaz_rhs = new double[n];
        double *cpu_dtSfx_rhs = new double[n]; double *cpu_dtSfy_rhs = new double[n]; double *cpu_dtSfz_rhs = new double[n];
        double *cpu_ham  = new double[n]; double *cpu_movx = new double[n]; double *cpu_movy = new double[n];
        double *cpu_movz = new double[n]; double *cpu_Gmx  = new double[n]; double *cpu_Gmy  = new double[n]; double *cpu_Gmz  = new double[n];

        int    sym_ref = 0;
        int    lev_ref = 0;
        double eps_ref = 0.0;
        int    co_ref  = 0;

        int ret = f_compute_rhs_bssn(
            ex, T,
            X, Y, Z,
            chi, trK,
            gxx, gxy, gxz, gyy, gyz, gzz,
            Axx, Axy, Axz, Ayy, Ayz, Azz,
            Gamx, Gamy, Gamz,
            Lap, betax, betay, betaz, dtSfx, dtSfy, dtSfz,
            cpu_chi_rhs, cpu_trK_rhs,
            cpu_gxx_rhs, cpu_gxy_rhs, cpu_gxz_rhs, cpu_gyy_rhs, cpu_gyz_rhs, cpu_gzz_rhs,
            cpu_Axx_rhs, cpu_Axy_rhs, cpu_Axz_rhs, cpu_Ayy_rhs, cpu_Ayz_rhs, cpu_Azz_rhs,
            cpu_Gamx_rhs, cpu_Gamy_rhs, cpu_Gamz_rhs,
            cpu_Lap_rhs, cpu_betax_rhs, cpu_betay_rhs, cpu_betaz_rhs,
            cpu_dtSfx_rhs, cpu_dtSfy_rhs, cpu_dtSfz_rhs,
            rho, Sx, Sy, Sz,
            Sxx, Sxy, Sxz, Syy, Syz, Szz,
            Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
            Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
            Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
            Ricxx, Ricxy, Ricxz, Ricyy, Ricyz, Riczz,
            cpu_ham, cpu_movx, cpu_movy, cpu_movz, cpu_Gmx, cpu_Gmy, cpu_Gmz,
            sym_ref, lev_ref, eps_ref, co_ref);

        printf("  CPU f_compute_rhs_bssn ret=%d\n", ret);

        // ---- compare RHS outputs (interior only) ----
        int fails = 0;
        fails += compare_field("chi_rhs",    gpu_chi_rhs,    cpu_chi_rhs,    ex, TOL);
        fails += compare_field("trK_rhs",    gpu_trK_rhs,    cpu_trK_rhs,    ex, TOL);
        fails += compare_field("gxx_rhs",    gpu_gxx_rhs,    cpu_gxx_rhs,    ex, TOL);
        fails += compare_field("gxy_rhs",    gpu_gxy_rhs,    cpu_gxy_rhs,    ex, TOL);
        fails += compare_field("gxz_rhs",    gpu_gxz_rhs,    cpu_gxz_rhs,    ex, TOL);
        fails += compare_field("gyy_rhs",    gpu_gyy_rhs,    cpu_gyy_rhs,    ex, TOL);
        fails += compare_field("gyz_rhs",    gpu_gyz_rhs,    cpu_gyz_rhs,    ex, TOL);
        fails += compare_field("gzz_rhs",    gpu_gzz_rhs,    cpu_gzz_rhs,    ex, TOL);
        fails += compare_field("Axx_rhs",    gpu_Axx_rhs,    cpu_Axx_rhs,    ex, TOL);
        fails += compare_field("Axy_rhs",    gpu_Axy_rhs,    cpu_Axy_rhs,    ex, TOL);
        fails += compare_field("Axz_rhs",    gpu_Axz_rhs,    cpu_Axz_rhs,    ex, TOL);
        fails += compare_field("Ayy_rhs",    gpu_Ayy_rhs,    cpu_Ayy_rhs,    ex, TOL);
        fails += compare_field("Ayz_rhs",    gpu_Ayz_rhs,    cpu_Ayz_rhs,    ex, TOL);
        fails += compare_field("Azz_rhs",    gpu_Azz_rhs,    cpu_Azz_rhs,    ex, TOL);
        fails += compare_field("Gamx_rhs",   gpu_Gamx_rhs,   cpu_Gamx_rhs,   ex, TOL);
        fails += compare_field("Gamy_rhs",   gpu_Gamy_rhs,   cpu_Gamy_rhs,   ex, TOL);
        fails += compare_field("Gamz_rhs",   gpu_Gamz_rhs,   cpu_Gamz_rhs,   ex, TOL);
        fails += compare_field("Lap_rhs",    gpu_Lap_rhs,    cpu_Lap_rhs,    ex, TOL);
        fails += compare_field("betax_rhs",  gpu_betax_rhs,  cpu_betax_rhs,  ex, TOL);
        fails += compare_field("betay_rhs",  gpu_betay_rhs,  cpu_betay_rhs,  ex, TOL);
        fails += compare_field("betaz_rhs",  gpu_betaz_rhs,  cpu_betaz_rhs,  ex, TOL);
        fails += compare_field("dtSfx_rhs",  gpu_dtSfx_rhs,  cpu_dtSfx_rhs,  ex, TOL);
        fails += compare_field("dtSfy_rhs",  gpu_dtSfy_rhs,  cpu_dtSfy_rhs,  ex, TOL);
        fails += compare_field("dtSfz_rhs",  gpu_dtSfz_rhs,  cpu_dtSfz_rhs,  ex, TOL);
        fails += compare_field("ham_Res",    gpu_ham,        cpu_ham,        ex, TOL);
        fails += compare_field("movx_Res",   gpu_movx,       cpu_movx,       ex, TOL);
        fails += compare_field("movy_Res",   gpu_movy,       cpu_movy,       ex, TOL);
        fails += compare_field("movz_Res",   gpu_movz,       cpu_movz,       ex, TOL);
        fails += compare_field("Gmx_Res",    gpu_Gmx,        cpu_Gmx,        ex, TOL);
        fails += compare_field("Gmy_Res",    gpu_Gmy,        cpu_Gmy,        ex, TOL);
        fails += compare_field("Gmz_Res",    gpu_Gmz,        cpu_Gmz,        ex, TOL);

        total_fails += fails;

        // ---- cleanup per-M host arrays ----
        delete[] gpu_chi_rhs;  delete[] gpu_trK_rhs;
        delete[] gpu_gxx_rhs;  delete[] gpu_gxy_rhs;  delete[] gpu_gxz_rhs;
        delete[] gpu_gyy_rhs;  delete[] gpu_gyz_rhs;  delete[] gpu_gzz_rhs;
        delete[] gpu_Axx_rhs;  delete[] gpu_Axy_rhs;  delete[] gpu_Axz_rhs;
        delete[] gpu_Ayy_rhs;  delete[] gpu_Ayz_rhs;  delete[] gpu_Azz_rhs;
        delete[] gpu_Gamx_rhs; delete[] gpu_Gamy_rhs; delete[] gpu_Gamz_rhs;
        delete[] gpu_Lap_rhs;
        delete[] gpu_betax_rhs; delete[] gpu_betay_rhs; delete[] gpu_betaz_rhs;
        delete[] gpu_dtSfx_rhs; delete[] gpu_dtSfy_rhs; delete[] gpu_dtSfz_rhs;
        delete[] gpu_ham;  delete[] gpu_movx; delete[] gpu_movy;
        delete[] gpu_movz; delete[] gpu_Gmx;  delete[] gpu_Gmy; delete[] gpu_Gmz;

        delete[] cpu_chi_rhs;  delete[] cpu_trK_rhs;
        delete[] cpu_gxx_rhs;  delete[] cpu_gxy_rhs;  delete[] cpu_gxz_rhs;
        delete[] cpu_gyy_rhs;  delete[] cpu_gyz_rhs;  delete[] cpu_gzz_rhs;
        delete[] cpu_Axx_rhs;  delete[] cpu_Axy_rhs;  delete[] cpu_Axz_rhs;
        delete[] cpu_Ayy_rhs;  delete[] cpu_Ayz_rhs;  delete[] cpu_Azz_rhs;
        delete[] cpu_Gamx_rhs; delete[] cpu_Gamy_rhs; delete[] cpu_Gamz_rhs;
        delete[] cpu_Lap_rhs;
        delete[] cpu_betax_rhs; delete[] cpu_betay_rhs; delete[] cpu_betaz_rhs;
        delete[] cpu_dtSfx_rhs; delete[] cpu_dtSfy_rhs; delete[] cpu_dtSfz_rhs;
        delete[] cpu_ham;  delete[] cpu_movx; delete[] cpu_movy;
        delete[] cpu_movz; delete[] cpu_Gmx;  delete[] cpu_Gmy; delete[] cpu_Gmz;
    }

    // ---- cleanup host ----
    delete[] X;    delete[] Y;    delete[] Z;
    delete[] chi;  delete[] trK;
    delete[] gxx;  delete[] gxy;  delete[] gxz;
    delete[] gyy;  delete[] gyz;  delete[] gzz;
    delete[] Axx;  delete[] Axy;  delete[] Axz;
    delete[] Ayy;  delete[] Ayz;  delete[] Azz;
    delete[] Gamx; delete[] Gamy; delete[] Gamz;
    delete[] Lap;  delete[] betax; delete[] betay; delete[] betaz;
    delete[] dtSfx; delete[] dtSfy; delete[] dtSfz;
    delete[] rho; delete[] Sx;  delete[] Sy;  delete[] Sz;
    delete[] Sxx; delete[] Sxy; delete[] Sxz;
    delete[] Syy; delete[] Syz; delete[] Szz;
    delete[] Gamxxx; delete[] Gamxxy; delete[] Gamxxz;
    delete[] Gamxyy; delete[] Gamxyz; delete[] Gamxzz;
    delete[] Gamyxx; delete[] Gamyxy; delete[] Gamyxz;
    delete[] Gamyyy; delete[] Gamyyz; delete[] Gamyzz;
    delete[] Gamzxx; delete[] Gamzxy; delete[] Gamzxz;
    delete[] Gamzyy; delete[] Gamzyz; delete[] Gamzzz;
    delete[] Ricxx; delete[] Ricxy; delete[] Ricxz;
    delete[] Ricyy; delete[] Ricyz; delete[] Riczz;

    // ---- cleanup GPU ----
    CUDA_CHECK(cudaFree(d_X));    CUDA_CHECK(cudaFree(d_Y));    CUDA_CHECK(cudaFree(d_Z));
    CUDA_CHECK(cudaFree(d_chi));  CUDA_CHECK(cudaFree(d_trK));
    CUDA_CHECK(cudaFree(d_gxx));  CUDA_CHECK(cudaFree(d_gxy));  CUDA_CHECK(cudaFree(d_gxz));
    CUDA_CHECK(cudaFree(d_gyy));  CUDA_CHECK(cudaFree(d_gyz));  CUDA_CHECK(cudaFree(d_gzz));
    CUDA_CHECK(cudaFree(d_Axx));  CUDA_CHECK(cudaFree(d_Axy));  CUDA_CHECK(cudaFree(d_Axz));
    CUDA_CHECK(cudaFree(d_Ayy));  CUDA_CHECK(cudaFree(d_Ayz));  CUDA_CHECK(cudaFree(d_Azz));
    CUDA_CHECK(cudaFree(d_Gamx)); CUDA_CHECK(cudaFree(d_Gamy)); CUDA_CHECK(cudaFree(d_Gamz));
    CUDA_CHECK(cudaFree(d_Lap));
    CUDA_CHECK(cudaFree(d_betax)); CUDA_CHECK(cudaFree(d_betay)); CUDA_CHECK(cudaFree(d_betaz));
    CUDA_CHECK(cudaFree(d_dtSfx)); CUDA_CHECK(cudaFree(d_dtSfy)); CUDA_CHECK(cudaFree(d_dtSfz));
    CUDA_CHECK(cudaFree(d_chi_rhs));  CUDA_CHECK(cudaFree(d_trK_rhs));
    CUDA_CHECK(cudaFree(d_gxx_rhs));  CUDA_CHECK(cudaFree(d_gxy_rhs));  CUDA_CHECK(cudaFree(d_gxz_rhs));
    CUDA_CHECK(cudaFree(d_gyy_rhs));  CUDA_CHECK(cudaFree(d_gyz_rhs));  CUDA_CHECK(cudaFree(d_gzz_rhs));
    CUDA_CHECK(cudaFree(d_Axx_rhs));  CUDA_CHECK(cudaFree(d_Axy_rhs));  CUDA_CHECK(cudaFree(d_Axz_rhs));
    CUDA_CHECK(cudaFree(d_Ayy_rhs));  CUDA_CHECK(cudaFree(d_Ayz_rhs));  CUDA_CHECK(cudaFree(d_Azz_rhs));
    CUDA_CHECK(cudaFree(d_Gamx_rhs)); CUDA_CHECK(cudaFree(d_Gamy_rhs)); CUDA_CHECK(cudaFree(d_Gamz_rhs));
    CUDA_CHECK(cudaFree(d_Lap_rhs));
    CUDA_CHECK(cudaFree(d_betax_rhs)); CUDA_CHECK(cudaFree(d_betay_rhs)); CUDA_CHECK(cudaFree(d_betaz_rhs));
    CUDA_CHECK(cudaFree(d_dtSfx_rhs)); CUDA_CHECK(cudaFree(d_dtSfy_rhs)); CUDA_CHECK(cudaFree(d_dtSfz_rhs));
    CUDA_CHECK(cudaFree(d_rho));   CUDA_CHECK(cudaFree(d_Sx));  CUDA_CHECK(cudaFree(d_Sy));  CUDA_CHECK(cudaFree(d_Sz));
    CUDA_CHECK(cudaFree(d_Sxx));   CUDA_CHECK(cudaFree(d_Sxy)); CUDA_CHECK(cudaFree(d_Sxz));
    CUDA_CHECK(cudaFree(d_Syy));   CUDA_CHECK(cudaFree(d_Syz)); CUDA_CHECK(cudaFree(d_Szz));
    CUDA_CHECK(cudaFree(d_Gamxxx)); CUDA_CHECK(cudaFree(d_Gamxxy)); CUDA_CHECK(cudaFree(d_Gamxxz));
    CUDA_CHECK(cudaFree(d_Gamxyy)); CUDA_CHECK(cudaFree(d_Gamxyz)); CUDA_CHECK(cudaFree(d_Gamxzz));
    CUDA_CHECK(cudaFree(d_Gamyxx)); CUDA_CHECK(cudaFree(d_Gamyxy)); CUDA_CHECK(cudaFree(d_Gamyxz));
    CUDA_CHECK(cudaFree(d_Gamyyy)); CUDA_CHECK(cudaFree(d_Gamyyz)); CUDA_CHECK(cudaFree(d_Gamyzz));
    CUDA_CHECK(cudaFree(d_Gamzxx)); CUDA_CHECK(cudaFree(d_Gamzxy)); CUDA_CHECK(cudaFree(d_Gamzxz));
    CUDA_CHECK(cudaFree(d_Gamzyy)); CUDA_CHECK(cudaFree(d_Gamzyz)); CUDA_CHECK(cudaFree(d_Gamzzz));
    CUDA_CHECK(cudaFree(d_Rxx));   CUDA_CHECK(cudaFree(d_Rxy));   CUDA_CHECK(cudaFree(d_Rxz));
    CUDA_CHECK(cudaFree(d_Ryy));   CUDA_CHECK(cudaFree(d_Ryz));   CUDA_CHECK(cudaFree(d_Rzz));
    CUDA_CHECK(cudaFree(d_ham));   CUDA_CHECK(cudaFree(d_movx));
    CUDA_CHECK(cudaFree(d_movy));  CUDA_CHECK(cudaFree(d_movz));
    CUDA_CHECK(cudaFree(d_Gmx));   CUDA_CHECK(cudaFree(d_Gmy));   CUDA_CHECK(cudaFree(d_Gmz));

    return total_fails;
}

// ===========================================================================
// Test 6: BSSN Symmetry BD (gpu_apply_symmetry_bd_fused_launch)
// All 19 BSSN fields have SYM,SYM,SYM parity -> no sign flip.
// Verify ghost zone = reflected interior.
// ===========================================================================
static int test_bssn_symmetry_bd()
{
    printf("[6/7] BSSN Symmetry BD (apply_symmetry_bd_fused)\n");

    int ex[3]  = {8 + 2*GW, 8 + 2*GW, 8 + 2*GW};
    int pex[3] = {ex[0], ex[1], ex[2]};
    int n = ex[0] * ex[1] * ex[2];

    // ---- allocate and fill host arrays (19 BSSN fields) ----
    double *chi   = new double[n]; double *trK   = new double[n];
    double *gxx   = new double[n]; double *gxy   = new double[n]; double *gxz   = new double[n];
    double *gyy   = new double[n]; double *gyz   = new double[n]; double *gzz   = new double[n];
    double *Axx   = new double[n]; double *Axy   = new double[n]; double *Axz   = new double[n];
    double *Ayy   = new double[n]; double *Ayz   = new double[n]; double *Azz   = new double[n];
    double *Gamx  = new double[n]; double *Gamy  = new double[n]; double *Gamz  = new double[n];
    double *Lap   = new double[n];
    double *betax = new double[n]; double *betay = new double[n]; double *betaz = new double[n];
    double *dtSfx = new double[n]; double *dtSfy = new double[n]; double *dtSfz = new double[n];

    double dx = 0.1;
    for (int k = 0; k < ex[2]; k++) {
        for (int j = 0; j < ex[1]; j++) {
            for (int i = 0; i < ex[0]; i++) {
                int idx = i + j * ex[0] + k * ex[0] * ex[1];
                double x = (i - GW) * dx;
                double y = (j - GW) * dx;
                double z = (k - GW) * dx;

                chi[idx]   =  sin(x) * cos(y) * sin(z) + 1.0;
                trK[idx]   =  cos(x) * sin(y) * cos(z) + 0.1;
                gxx[idx]   =  cos(x) * cos(y) * sin(z) + 0.2;
                gxy[idx]   =  sin(x) * sin(y) * cos(z) + 0.3;
                gxz[idx]   =  cos(x) * sin(y) * sin(z) + 0.4;
                gyy[idx]   =  sin(x) * cos(y) * cos(z) + 0.5;
                gyz[idx]   =  cos(x) * cos(y) * cos(z) + 0.6;
                gzz[idx]   =  sin(x) * sin(y) * sin(z) + 0.7;
                Axx[idx]   =  sin(x) * sin(y) * cos(z) + 1.0;
                Axy[idx]   =  cos(x) * cos(y) * sin(z) + 1.1;
                Axz[idx]   =  sin(x) * cos(y) * cos(z) + 1.2;
                Ayy[idx]   =  cos(x) * sin(y) * sin(z) + 1.3;
                Ayz[idx]   =  sin(x) * sin(y) * sin(z) + 1.4;
                Azz[idx]   =  cos(x) * cos(y) * cos(z) + 1.5;
                Gamx[idx]  =  sin(x + y) * cos(z)     + 2.0;
                Gamy[idx]  =  cos(x) * sin(y + z)     + 2.1;
                Gamz[idx]  =  sin(x + y + z)          + 2.2;
                Lap[idx]   =  sin(x) * sin(y) * sin(z) + 1.0;
                betax[idx] =  cos(y) * sin(z)          + 0.1;
                betay[idx] =  sin(x) * cos(z)          + 0.2;
                betaz[idx] =  cos(x) * sin(y)          + 0.3;
                dtSfx[idx] =  sin(x) * cos(y) * cos(z) + 0.1;
                dtSfy[idx] =  cos(x) * sin(y) * cos(z) + 0.2;
                dtSfz[idx] =  cos(x) * cos(y) * sin(z) + 0.3;
            }
        }
    }

    // ---- allocate GPU memory and copy ----
    size_t nb = n * sizeof(double);
    double *d_chi, *d_trK;
    double *d_gxx, *d_gxy, *d_gxz, *d_gyy, *d_gyz, *d_gzz;
    double *d_Axx, *d_Axy, *d_Axz, *d_Ayy, *d_Ayz, *d_Azz;
    double *d_Gamx, *d_Gamy, *d_Gamz, *d_Lap;
    double *d_betax, *d_betay, *d_betaz;
    double *d_dtSfx, *d_dtSfy, *d_dtSfz;

    CUDA_CHECK(cudaMalloc(&d_chi,   nb)); CUDA_CHECK(cudaMalloc(&d_trK,  nb));
    CUDA_CHECK(cudaMalloc(&d_gxx,   nb)); CUDA_CHECK(cudaMalloc(&d_gxy,   nb)); CUDA_CHECK(cudaMalloc(&d_gxz,   nb));
    CUDA_CHECK(cudaMalloc(&d_gyy,   nb)); CUDA_CHECK(cudaMalloc(&d_gyz,   nb)); CUDA_CHECK(cudaMalloc(&d_gzz,   nb));
    CUDA_CHECK(cudaMalloc(&d_Axx,   nb)); CUDA_CHECK(cudaMalloc(&d_Axy,   nb)); CUDA_CHECK(cudaMalloc(&d_Axz,   nb));
    CUDA_CHECK(cudaMalloc(&d_Ayy,   nb)); CUDA_CHECK(cudaMalloc(&d_Ayz,   nb)); CUDA_CHECK(cudaMalloc(&d_Azz,   nb));
    CUDA_CHECK(cudaMalloc(&d_Gamx,  nb)); CUDA_CHECK(cudaMalloc(&d_Gamy,  nb)); CUDA_CHECK(cudaMalloc(&d_Gamz,  nb));
    CUDA_CHECK(cudaMalloc(&d_Lap,   nb));
    CUDA_CHECK(cudaMalloc(&d_betax, nb)); CUDA_CHECK(cudaMalloc(&d_betay, nb)); CUDA_CHECK(cudaMalloc(&d_betaz, nb));
    CUDA_CHECK(cudaMalloc(&d_dtSfx, nb)); CUDA_CHECK(cudaMalloc(&d_dtSfy, nb)); CUDA_CHECK(cudaMalloc(&d_dtSfz, nb));

    cudaStream_t stream = 0;
    CUDA_CHECK(cudaMemcpy(d_chi,   chi,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_trK,   trK,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gxx,   gxx,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gxy,   gxy,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gxz,   gxz,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gyy,   gyy,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gyz,   gyz,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gzz,   gzz,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Axx,   Axx,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Axy,   Axy,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Axz,   Axz,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Ayy,   Ayy,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Ayz,   Ayz,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Azz,   Azz,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamx,  Gamx,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamy,  Gamy,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Gamz,  Gamz,  nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Lap,   Lap,   nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_betax, betax, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_betay, betay, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_betaz, betaz, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_dtSfx, dtSfx, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_dtSfy, dtSfy, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_dtSfz, dtSfz, nb, cudaMemcpyHostToDevice));

    // ---- run GPU symmetry BD ----
    gpu_apply_symmetry_bd_fused_launch(
        stream,
        ex[0], ex[1], ex[2],
        pex[0], pex[1], pex[2],
        d_chi, d_trK,
        d_gxx, d_gxy, d_gxz, d_gyy, d_gyz, d_gzz,
        d_Axx, d_Axy, d_Axz, d_Ayy, d_Ayz, d_Azz,
        d_Gamx, d_Gamy, d_Gamz, d_Lap,
        d_betax, d_betay, d_betaz,
        d_dtSfx, d_dtSfy, d_dtSfz);
    CUDA_CHECK(cudaDeviceSynchronize());

    // ---- copy results back ----
    CUDA_CHECK(cudaMemcpy(chi,   d_chi,   nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(trK,   d_trK,   nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gxx,   d_gxx,   nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gxy,   d_gxy,   nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gxz,   d_gxz,   nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gyy,   d_gyy,   nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gyz,   d_gyz,   nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gzz,   d_gzz,   nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(Axx,   d_Axx,   nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(Axy,   d_Axy,   nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(Axz,   d_Axz,   nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(Ayy,   d_Ayy,   nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(Ayz,   d_Ayz,   nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(Azz,   d_Azz,   nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(Gamx,  d_Gamx,  nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(Gamy,  d_Gamy,  nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(Gamz,  d_Gamz,  nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(Lap,   d_Lap,   nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(betax, d_betax, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(betay, d_betay, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(betaz, d_betaz, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(dtSfx, d_dtSfx, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(dtSfy, d_dtSfy, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(dtSfz, d_dtSfz, nb, cudaMemcpyDeviceToHost));

    // ---- verify: all 19 fields have SYM,SYM,SYM parity -> ghost = reflected interior ----
    // Check lower z-face ghost zone: k < GW, i,j in [GW, ex-GW)
    struct {
        const char *name;
        double *data;
    } fields[] = {
        {"chi",   chi},
        {"trK",   trK},
        {"gxx",   gxx},   {"gxy",   gxy},   {"gxz",   gxz},
        {"gyy",   gyy},   {"gyz",   gyz},   {"gzz",   gzz},
        {"Axx",   Axx},   {"Axy",   Axy},   {"Axz",   Axz},
        {"Ayy",   Ayy},   {"Ayz",   Ayz},   {"Azz",   Azz},
        {"Gamx",  Gamx},  {"Gamy",  Gamy},  {"Gamz",  Gamz},
        {"Lap",   Lap},
        {"betax", betax}, {"betay", betay}, {"betaz", betaz},
        {"dtSfx", dtSfx}, {"dtSfy", dtSfy}, {"dtSfz", dtSfz},
    };
    const int nfields = sizeof(fields) / sizeof(fields[0]);

    int fails = 0;
    for (int f = 0; f < nfields; f++) {
        double maxerr = 0.0;
        int checked = 0;

        for (int k = 0; k < GW; k++) {
            int rk_raw = GW - k - 1;  // reflected raw index: padded_source = GW + rk_raw
            for (int j = GW; j < ex[1] - GW; j++) {
                for (int i = GW; i < ex[0] - GW; i++) {
                    int gidx = i + j * ex[0] + k * ex[0] * ex[1];
                    int sidx = i + j * ex[0] + (rk_raw + GW) * ex[0] * ex[1];

                    double expected = fields[f].data[sidx];  // no sign flip: all SYM
                    double diff = fabs(fields[f].data[gidx] - expected);
                    if (diff > maxerr) maxerr = diff;
                    checked++;
                }
            }
        }

        printf("  max|D%-14s|= %#.2e  [%s]  (checked=%d, z-face ghost)\n",
               fields[f].name, maxerr, (maxerr > TOL) ? "FAIL" : "PASS", checked);
        if (maxerr > TOL) fails++;
    }

    // ---- cleanup ----
    delete[] chi;   delete[] trK;
    delete[] gxx;   delete[] gxy;   delete[] gxz;
    delete[] gyy;   delete[] gyz;   delete[] gzz;
    delete[] Axx;   delete[] Axy;   delete[] Axz;
    delete[] Ayy;   delete[] Ayz;   delete[] Azz;
    delete[] Gamx;  delete[] Gamy;  delete[] Gamz;
    delete[] Lap;   delete[] betax; delete[] betay; delete[] betaz;
    delete[] dtSfx; delete[] dtSfy; delete[] dtSfz;

    CUDA_CHECK(cudaFree(d_chi));   CUDA_CHECK(cudaFree(d_trK));
    CUDA_CHECK(cudaFree(d_gxx));   CUDA_CHECK(cudaFree(d_gxy));   CUDA_CHECK(cudaFree(d_gxz));
    CUDA_CHECK(cudaFree(d_gyy));   CUDA_CHECK(cudaFree(d_gyz));   CUDA_CHECK(cudaFree(d_gzz));
    CUDA_CHECK(cudaFree(d_Axx));   CUDA_CHECK(cudaFree(d_Axy));   CUDA_CHECK(cudaFree(d_Axz));
    CUDA_CHECK(cudaFree(d_Ayy));   CUDA_CHECK(cudaFree(d_Ayz));   CUDA_CHECK(cudaFree(d_Azz));
    CUDA_CHECK(cudaFree(d_Gamx));  CUDA_CHECK(cudaFree(d_Gamy));  CUDA_CHECK(cudaFree(d_Gamz));
    CUDA_CHECK(cudaFree(d_Lap));
    CUDA_CHECK(cudaFree(d_betax)); CUDA_CHECK(cudaFree(d_betay)); CUDA_CHECK(cudaFree(d_betaz));
    CUDA_CHECK(cudaFree(d_dtSfx)); CUDA_CHECK(cudaFree(d_dtSfy)); CUDA_CHECK(cudaFree(d_dtSfz));

    return fails;
}

// ===========================================================================
// compute_determinant_3x3 — compute determinant of a 3x3 symmetric matrix
// ===========================================================================
static double compute_determinant_3x3(double gxx, double gxy, double gxz,
                                       double gyy, double gyz, double gzz)
{
    return gxx * (gyy * gzz - gyz * gyz)
         - gxy * (gxy * gzz - gyz * gxz)
         + gxz * (gxy * gyz - gyy * gxz);
}

// ===========================================================================
// Test 7: enforce_ga (gpu_enforce_ga_launch)
// This enforces det(g) = 1 and trace-free A_ij.
// GPU vs CPU comparison + property verification.
// ===========================================================================
static int test_enforce_ga()
{
    printf("[7/7] Enforce GA (gpu_enforce_ga_launch vs f_enforce_ga)\n");

    int ex[3]  = {8 + 2*GW, 8 + 2*GW, 8 + 2*GW};
    int pex[3] = {ex[0], ex[1], ex[2]};
    int n = ex[0] * ex[1] * ex[2];

    // ---- allocate two sets: one for GPU, one for CPU ----
    double *gpu_gxx = new double[n]; double *gpu_gxy = new double[n]; double *gpu_gxz = new double[n];
    double *gpu_gyy = new double[n]; double *gpu_gyz = new double[n]; double *gpu_gzz = new double[n];
    double *gpu_Axx = new double[n]; double *gpu_Axy = new double[n]; double *gpu_Axz = new double[n];
    double *gpu_Ayy = new double[n]; double *gpu_Ayz = new double[n]; double *gpu_Azz = new double[n];

    double *cpu_gxx = new double[n]; double *cpu_gxy = new double[n]; double *cpu_gxz = new double[n];
    double *cpu_gyy = new double[n]; double *cpu_gyz = new double[n]; double *cpu_gzz = new double[n];
    double *cpu_Axx = new double[n]; double *cpu_Axy = new double[n]; double *cpu_Axz = new double[n];
    double *cpu_Ayy = new double[n]; double *cpu_Ayz = new double[n]; double *cpu_Azz = new double[n];

    // ---- fill with conformal metric that has det != 1 and A_ij with tr != 0 ----
    double dx = 0.1;
    for (int k = 0; k < ex[2]; k++) {
        for (int j = 0; j < ex[1]; j++) {
            for (int i = 0; i < ex[0]; i++) {
                int idx = i + j * ex[0] + k * ex[0] * ex[1];
                double x = (i - GW) * dx;
                double y = (j - GW) * dx;
                double z = (k - GW) * dx;

                // Construct a conformal metric with det != 1
                // g_ij = delta_ij + epsilon * h_ij where epsilon ~ 0.1
                // gives det ~ 1 + epsilon*tr(h) + O(epsilon^2)
                double hxx = sin(x) * cos(y);
                double hxy = 0.1 * cos(x) * sin(y);
                double hxz = 0.05 * sin(x) * sin(z);
                double hyy = cos(x) * sin(z);
                double hyz = 0.08 * cos(y) * cos(z);
                double hzz = sin(y) * cos(z);

                double eps_metric = 0.1;
                gpu_gxx[idx] = 1.0 + eps_metric * hxx;
                gpu_gxy[idx] = eps_metric * hxy;
                gpu_gxz[idx] = eps_metric * hxz;
                gpu_gyy[idx] = 1.0 + eps_metric * hyy;
                gpu_gyz[idx] = eps_metric * hyz;
                gpu_gzz[idx] = 1.0 + eps_metric * hzz;

                // A_ij with trace != 0
                gpu_Axx[idx] = sin(x) * cos(y);
                gpu_Axy[idx] = 0.5 * cos(x) * sin(z);
                gpu_Axz[idx] = 0.3 * sin(y) * cos(z);
                gpu_Ayy[idx] = cos(x) * sin(y);
                gpu_Ayz[idx] = 0.4 * sin(x) * sin(z);
                gpu_Azz[idx] = sin(y) * cos(z);

                // Same for CPU
                cpu_gxx[idx] = gpu_gxx[idx];
                cpu_gxy[idx] = gpu_gxy[idx];
                cpu_gxz[idx] = gpu_gxz[idx];
                cpu_gyy[idx] = gpu_gyy[idx];
                cpu_gyz[idx] = gpu_gyz[idx];
                cpu_gzz[idx] = gpu_gzz[idx];
                cpu_Axx[idx] = gpu_Axx[idx];
                cpu_Axy[idx] = gpu_Axy[idx];
                cpu_Axz[idx] = gpu_Axz[idx];
                cpu_Ayy[idx] = gpu_Ayy[idx];
                cpu_Ayz[idx] = gpu_Ayz[idx];
                cpu_Azz[idx] = gpu_Azz[idx];
            }
        }
    }

    // ---- run CPU Fortran function first ----
    int ret_int = ex[0] * ex[1] * ex[2];
    f_enforce_ga(&ret_int,
                 cpu_gxx, cpu_gxy, cpu_gxz, cpu_gyy, cpu_gyz, cpu_gzz,
                 cpu_Axx, cpu_Axy, cpu_Axz, cpu_Ayy, cpu_Ayz, cpu_Azz);

    // ---- allocate GPU memory, copy, run, copy back ----
    size_t nb = n * sizeof(double);
    double *d_gxx, *d_gxy, *d_gxz, *d_gyy, *d_gyz, *d_gzz;
    double *d_Axx, *d_Axy, *d_Axz, *d_Ayy, *d_Ayz, *d_Azz;

    CUDA_CHECK(cudaMalloc(&d_gxx, nb)); CUDA_CHECK(cudaMalloc(&d_gxy, nb)); CUDA_CHECK(cudaMalloc(&d_gxz, nb));
    CUDA_CHECK(cudaMalloc(&d_gyy, nb)); CUDA_CHECK(cudaMalloc(&d_gyz, nb)); CUDA_CHECK(cudaMalloc(&d_gzz, nb));
    CUDA_CHECK(cudaMalloc(&d_Axx, nb)); CUDA_CHECK(cudaMalloc(&d_Axy, nb)); CUDA_CHECK(cudaMalloc(&d_Axz, nb));
    CUDA_CHECK(cudaMalloc(&d_Ayy, nb)); CUDA_CHECK(cudaMalloc(&d_Ayz, nb)); CUDA_CHECK(cudaMalloc(&d_Azz, nb));

    cudaStream_t stream = 0;
    CUDA_CHECK(cudaMemcpy(d_gxx, gpu_gxx, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gxy, gpu_gxy, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gxz, gpu_gxz, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gyy, gpu_gyy, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gyz, gpu_gyz, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_gzz, gpu_gzz, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Axx, gpu_Axx, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Axy, gpu_Axy, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Axz, gpu_Axz, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Ayy, gpu_Ayy, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Ayz, gpu_Ayz, nb, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_Azz, gpu_Azz, nb, cudaMemcpyHostToDevice));

    gpu_enforce_ga_launch(stream, ex, pex,
                          d_gxx, d_gxy, d_gxz, d_gyy, d_gyz, d_gzz,
                          d_Axx, d_Axy, d_Axz, d_Ayy, d_Ayz, d_Azz);
    CUDA_CHECK(cudaDeviceSynchronize());

    CUDA_CHECK(cudaMemcpy(gpu_gxx, d_gxx, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_gxy, d_gxy, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_gxz, d_gxz, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_gyy, d_gyy, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_gyz, d_gyz, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_gzz, d_gzz, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Axx, d_Axx, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Axy, d_Axy, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Axz, d_Axz, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Ayy, d_Ayy, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Ayz, d_Ayz, nb, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(gpu_Azz, d_Azz, nb, cudaMemcpyDeviceToHost));

    // ---- compare GPU vs CPU outputs (interior only) ----
    int fails = 0;
    fails += compare_field("gxx(ga)", gpu_gxx, cpu_gxx, ex, TOL);
    fails += compare_field("gxy(ga)", gpu_gxy, cpu_gxy, ex, TOL);
    fails += compare_field("gxz(ga)", gpu_gxz, cpu_gxz, ex, TOL);
    fails += compare_field("gyy(ga)", gpu_gyy, cpu_gyy, ex, TOL);
    fails += compare_field("gyz(ga)", gpu_gyz, cpu_gyz, ex, TOL);
    fails += compare_field("gzz(ga)", gpu_gzz, cpu_gzz, ex, TOL);
    fails += compare_field("Axx(ga)", gpu_Axx, cpu_Axx, ex, TOL);
    fails += compare_field("Axy(ga)", gpu_Axy, cpu_Axy, ex, TOL);
    fails += compare_field("Axz(ga)", gpu_Axz, cpu_Axz, ex, TOL);
    fails += compare_field("Ayy(ga)", gpu_Ayy, cpu_Ayy, ex, TOL);
    fails += compare_field("Ayz(ga)", gpu_Ayz, cpu_Ayz, ex, TOL);
    fails += compare_field("Azz(ga)", gpu_Azz, cpu_Azz, ex, TOL);

    // ---- verify mathematical properties on GPU output (interior only) ----
    double max_det_err = 0.0;
    double max_trA_err = 0.0;
    int interior_count = 0;
    for (int k = GW; k < ex[2] - GW; k++) {
        for (int j = GW; j < ex[1] - GW; j++) {
            for (int i = GW; i < ex[0] - GW; i++) {
                int idx = i + j * ex[0] + k * ex[0] * ex[1];
                double det = compute_determinant_3x3(gpu_gxx[idx], gpu_gxy[idx], gpu_gxz[idx],
                                                      gpu_gyy[idx], gpu_gyz[idx], gpu_gzz[idx]);
                double det_err = fabs(det - 1.0);
                if (det_err > max_det_err) max_det_err = det_err;

                double trA = gpu_Axx[idx] + gpu_Ayy[idx] + gpu_Azz[idx];
                double trA_err = fabs(trA);
                if (trA_err > max_trA_err) max_trA_err = trA_err;
                interior_count++;
            }
        }
    }
    printf("  max|det(g)-1|      = %#.2e  [%s]  (interior=%d)\n",
           max_det_err, (max_det_err > TOL) ? "FAIL" : "PASS", interior_count);
    printf("  max|tr(A)|         = %#.2e  [%s]  (interior=%d)\n",
           max_trA_err, (max_trA_err > TOL) ? "FAIL" : "PASS", interior_count);
    if (max_det_err > TOL) fails++;
    if (max_trA_err > TOL) fails++;

    // ---- cleanup ----
    delete[] gpu_gxx; delete[] gpu_gxy; delete[] gpu_gxz;
    delete[] gpu_gyy; delete[] gpu_gyz; delete[] gpu_gzz;
    delete[] gpu_Axx; delete[] gpu_Axy; delete[] gpu_Axz;
    delete[] gpu_Ayy; delete[] gpu_Ayz; delete[] gpu_Azz;
    delete[] cpu_gxx; delete[] cpu_gxy; delete[] cpu_gxz;
    delete[] cpu_gyy; delete[] cpu_gyz; delete[] cpu_gzz;
    delete[] cpu_Axx; delete[] cpu_Axy; delete[] cpu_Axz;
    delete[] cpu_Ayy; delete[] cpu_Ayz; delete[] cpu_Azz;

    CUDA_CHECK(cudaFree(d_gxx)); CUDA_CHECK(cudaFree(d_gxy)); CUDA_CHECK(cudaFree(d_gxz));
    CUDA_CHECK(cudaFree(d_gyy)); CUDA_CHECK(cudaFree(d_gyz)); CUDA_CHECK(cudaFree(d_gzz));
    CUDA_CHECK(cudaFree(d_Axx)); CUDA_CHECK(cudaFree(d_Axy)); CUDA_CHECK(cudaFree(d_Axz));
    CUDA_CHECK(cudaFree(d_Ayy)); CUDA_CHECK(cudaFree(d_Ayz)); CUDA_CHECK(cudaFree(d_Azz));

    return fails;
}

// ===========================================================================
// main
// ===========================================================================
int main()
{
    int ex[3] = {8 + 2*GW, 8 + 2*GW, 8 + 2*GW};
    int interior_size = (ex[0] - 2*GW) * (ex[1] - 2*GW) * (ex[2] - 2*GW);

    printf("=== EM + BSSN GPU Kernel Accuracy Test ===\n");
    printf("Grid: %dx%dx%d, interior: %dx%dx%d, GW=%d, tol=%.0e\n\n",
           ex[0], ex[1], ex[2],
           ex[0]-2*GW, ex[1]-2*GW, ex[2]-2*GW,
           GW, TOL);

    int total_fails = 0;
    total_fails += test_em_rhs();
    printf("\n");
    total_fails += test_em_symmetry_bd();
    printf("\n");
    total_fails += test_np_phi2();
    printf("\n");
    total_fails += test_np_phi1();
    printf("\n");
    total_fails += test_bssn_rhs();
    printf("\n");
    total_fails += test_bssn_rhs_puncture();
    printf("\n");
    total_fails += test_bssn_symmetry_bd();
    printf("\n");
    total_fails += test_enforce_ga();
    printf("\n");

    if (total_fails == 0) {
        printf("RESULT: ALL PASS (tolerance=%.0e)\n", TOL);
    } else {
        printf("RESULT: %d FAILURES (tolerance=%.0e)\n", total_fails, TOL);
    }

    return (total_fails == 0) ? 0 : 1;
}
