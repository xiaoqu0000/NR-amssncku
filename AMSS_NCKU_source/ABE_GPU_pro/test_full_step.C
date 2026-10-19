// ===========================================================================
// test_full_step.C — Standalone full RK4 step test: GPU vs CPU comparison
//
// Compares the ENTIRE single Step_GPU output (predictor + 3 correctors,
// including EM + BSSN RHS, RK4, Sommerfeld, enforce_ga, symmetry BD)
// against an equivalent CPU computation using Fortran functions.
//
// Must be compiled as part of the GPU build (see Makefile).
//
// Grid:  GW=4 (6th-order), internal=8, dx=0.1, dT=0.05
// ===========================================================================

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <cuda_runtime.h>

#ifndef USE_GPU
#define USE_GPU
#endif

#include "empart.h"
#include "bssn_rhs.h"
#include "enforce_algebra.h"
#include "fmisc.h"
#include "rungekutta4_rout.h"
#include "sommerfeld_rout.h"
#include "gpu_fusion.h"

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
// Grid and tolerance constants
// ===========================================================================
static const int GW = 4;
static const double TOL = 1e-10;
static double dT = 0.05;
static double chitiny_val = -0.999;
static const double dx = 0.1;
static double chitiny = 1e-12;  // non-const for Fortran double& params

// ===========================================================================
// Variable indices for grouped arrays
// ===========================================================================
// BSSN state (23 vars): indices 0..22
enum {
    I_CHI = 0, I_TRK,
    I_GXX, I_GXY, I_GXZ, I_GYY, I_GYZ, I_GZZ,
    I_AXX, I_AXY, I_AXZ, I_AYY, I_AYZ, I_AZZ,
    I_GAMX, I_GAMY, I_GAMZ,
    I_LAP,
    I_BETAX, I_BETAY, I_BETAZ,
    I_DTSFX, I_DTSFY, I_DTSFZ
};
#define NVAR_BSSN 26

// EM state (8 vars): indices 26..33 (=NVAR_BSSN .. NVAR_BSSN+7)
enum {
    I_Ex = NVAR_BSSN, I_Ey, I_Ez, I_Bx, I_By, I_Bz, I_KPSI, I_KPHI
};
#define NVAR_EM 8
#define NVAR_TOTAL (NVAR_BSSN + NVAR_EM)  // 34

// Extra variable indices for grouped per-variable ops
#define NVAR_EXTRA_BSSN 0
#define NVAR_EXTRA_EM 0

// ===========================================================================
// SoA (symmetry-on-axis) lookup tables
// ===========================================================================
// BSSN SoA values matching var initialization in bssn_class.C
static double soa_bssn[NVAR_BSSN][3] = {
    // chi, trK
    { 1, 1, 1}, { 1, 1, 1},
    // gij
    { 1, 1, 1}, {-1,-1, 1}, {-1, 1,-1},
    { 1, 1, 1}, { 1,-1,-1}, { 1, 1, 1},
    // Aij
    { 1, 1, 1}, {-1,-1, 1}, {-1, 1,-1},
    { 1, 1, 1}, { 1,-1,-1}, { 1, 1, 1},
    // Gam_i
    {-1, 1, 1}, { 1,-1, 1}, { 1, 1,-1},
    // Lap
    { 1, 1, 1},
    // betai
    {-1, 1, 1}, { 1,-1, 1}, { 1, 1,-1},
    // dtSfi
    {-1, 1, 1}, { 1,-1, 1}, { 1, 1,-1}
};

// EM SoA values matching var initialization in bssnEM_class.C
static double soa_em[NVAR_EM][3] = {
    // Ex, Ey, Ez
    {-1, 1, 1}, { 1,-1, 1}, { 1, 1,-1},
    // Bx, By, Bz
    { 1,-1,-1}, {-1, 1,-1}, {-1,-1, 1},
    // Kpsi, Kphi
    { 1, 1, 1}, { 1, 1, 1}
};

// propspeed for Sommerfeld — matches bssn_class.C initialization
// Only chi, trK, Lap have a non-zero propspeed. All others use 0 (freeze BC).
// The actual value is set to v (light speed = 1.0 typically for BSSN gauge vars).
static const double propspeed_bssn[NVAR_BSSN] = {
    1.0, 1.0,    // chi, trK
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  // gij
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  // Aij
    0.0, 0.0, 0.0,                  // Gam_i
    1.0,                             // Lap
    0.0, 0.0, 0.0,                  // betai
    0.0, 0.0, 0.0                   // dtSfi
};

static const double propspeed_em[NVAR_EM] = {
    1.0, 1.0, 1.0,   // Ex, Ey, Ez
    1.0, 1.0, 1.0,   // Bx, By, Bz
    0.0, 0.0          // Kpsi, Kphi
};

// ===========================================================================
// fill_test_data — fill all grid arrays with smooth analytic test functions
// Physical domain [GW, ex[i]-GW) has data; ghost zones get 0.
// ===========================================================================
static void fill_test_data(int ex[3],
                           double *X, double *Y, double *Z,
                           // BSSN state (23)
                           double *chi, double *trK,
                           double *gxx, double *gxy, double *gxz,
                           double *gyy, double *gyz, double *gzz,
                           double *Axx, double *Axy, double *Axz,
                           double *Ayy, double *Ayz, double *Azz,
                           double *Gamx, double *Gamy, double *Gamz,
                           double *Lap,
                           double *betax, double *betay, double *betaz,
                           double *dtSfx, double *dtSfy, double *dtSfz,
                           // EM state (8)
                           double *Ex, double *Ey, double *Ez,
                           double *Bx, double *By, double *Bz,
                           double *Kpsi, double *Kphi,
                           // EM source (4) — vacuum
                           double *Jx, double *Jy, double *Jz, double *qchar,
                           // stress-energy (10)
                           double *rho, double *Sx, double *Sy, double *Sz,
                           double *Sxx, double *Sxy, double *Sxz,
                           double *Syy, double *Syz, double *Szz,
                           // Christoffel (18)
                           double *Gamxxx, double *Gamxxy, double *Gamxxz,
                           double *Gamxyy, double *Gamxyz, double *Gamxzz,
                           double *Gamyxx, double *Gamyxy, double *Gamyxz,
                           double *Gamyyy, double *Gamyyz, double *Gamyzz,
                           double *Gamzxx, double *Gamzxy, double *Gamzxz,
                           double *Gamzyy, double *Gamzyz, double *Gamzzz,
                           // Ricci (6)
                           double *Ricxx, double *Ricxy, double *Ricxz,
                           double *Ricyy, double *Ricyz, double *Riczz)
{
    for (int k = 0; k < ex[2]; k++) {
        for (int j = 0; j < ex[1]; j++) {
            for (int i = 0; i < ex[0]; i++) {
                int idx = i + j * ex[0] + k * ex[0] * ex[1];

                // Coordinates (global for ghost-aware filling)
                double x = (i - GW) * dx;
                double y = (j - GW) * dx;
                double z = (k - GW) * dx;
                X[idx] = x;
                Y[idx] = y;
                Z[idx] = z;

                // Only fill physical domain with non-zero data; ghost = 0
                bool interior = (i >= GW && i < ex[0] - GW &&
                                 j >= GW && j < ex[1] - GW &&
                                 k >= GW && k < ex[2] - GW);

                // BSSN state
                if (interior) {
                    chi[idx]   = 1.0 + 0.01 * sin(x) * cos(y);
                    trK[idx]   = 0.001 * sin(z);

                    gxx[idx]   = 0.0;
                    gxy[idx]   = 0.01 * sin(x) * cos(y) * sin(z);
                    gxz[idx]   = 0.005 * cos(x) * sin(y) * cos(z);
                    gyy[idx]   = 0.0;
                    gyz[idx]   = 0.008 * sin(x) * sin(y) * cos(z);
                    gzz[idx]   = 0.0;

                    Axx[idx]   =  0.001 * cos(x) * sin(y) * sin(z);
                    Axy[idx]   =  0.002 * sin(x) * sin(y) * cos(z);
                    Axz[idx]   =  0.0015 * cos(x) * cos(y) * sin(z);
                    Ayy[idx]   = -0.001 * cos(x) * sin(y) * sin(z);
                    Ayz[idx]   =  0.003 * sin(x) * cos(y) * cos(z);
                    Azz[idx]   =  0.0;

                    Gamx[idx]  = 0.001 * sin(x + y) * cos(z);
                    Gamy[idx]  = 0.001 * cos(x) * sin(y + z);
                    Gamz[idx]  = 0.001 * cos(x + z) * sin(y);

                    Lap[idx]   = 1.0;
                    betax[idx] = 0.0;
                    betay[idx] = 0.0;
                    betaz[idx] = 0.0;
                    dtSfx[idx] = 0.0;
                    dtSfy[idx] = 0.0;
                    dtSfz[idx] = 0.0;

                    // EM state
                    Ex[idx]   = sin(y) * cos(z);
                    Ey[idx]   = cos(x) * sin(z);
                    Ez[idx]   = sin(x) * cos(y);
                    Bx[idx]   = cos(y) * sin(z);
                    By[idx]   = sin(x) * cos(z);
                    Bz[idx]   = cos(x) * sin(y);
                    Kpsi[idx] = 0.0;
                    Kphi[idx] = 0.0;

                    // EM source (vacuum)
                    Jx[idx]   = 0.0;
                    Jy[idx]   = 0.0;
                    Jz[idx]   = 0.0;
                    qchar[idx] = 0.0;

                    // stress-energy
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

                    // Christoffel
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

                    // Ricci
                    Ricxx[idx] = 0.0001 * sin(x) * cos(y) * cos(z);
                    Ricxy[idx] = 0.0001 * cos(x) * sin(y) * cos(z);
                    Ricxz[idx] = 0.0001 * cos(x) * cos(y) * sin(z);
                    Ricyy[idx] = 0.0001 * sin(x) * sin(y) * cos(z);
                    Ricyz[idx] = 0.0001 * sin(x) * cos(y) * sin(z);
                    Riczz[idx] = 0.0001 * cos(x) * sin(y) * sin(z);
                } else {
                    // Ghost zones: zero everything
                    chi[idx]   = 0.0;  trK[idx]   = 0.0;
                    gxx[idx]   = 0.0;  gxy[idx]   = 0.0;  gxz[idx]   = 0.0;
                    gyy[idx]   = 0.0;  gyz[idx]   = 0.0;  gzz[idx]   = 0.0;
                    Axx[idx]   = 0.0;  Axy[idx]   = 0.0;  Axz[idx]   = 0.0;
                    Ayy[idx]   = 0.0;  Ayz[idx]   = 0.0;  Azz[idx]   = 0.0;
                    Gamx[idx]  = 0.0;  Gamy[idx]  = 0.0;  Gamz[idx]  = 0.0;
                    Lap[idx]   = 0.0;
                    betax[idx] = 0.0;  betay[idx] = 0.0;  betaz[idx] = 0.0;
                    dtSfx[idx] = 0.0;  dtSfy[idx] = 0.0;  dtSfz[idx] = 0.0;
                    Ex[idx]   = 0.0;  Ey[idx]   = 0.0;  Ez[idx]   = 0.0;
                    Bx[idx]   = 0.0;  By[idx]   = 0.0;  Bz[idx]   = 0.0;
                    Kpsi[idx] = 0.0;  Kphi[idx] = 0.0;
                    Jx[idx]   = 0.0;  Jy[idx]   = 0.0;  Jz[idx]   = 0.0;  qchar[idx] = 0.0;
                    rho[idx]  = 0.0;  Sx[idx]   = 0.0;  Sy[idx]   = 0.0;  Sz[idx]   = 0.0;
                    Sxx[idx]  = 0.0;  Sxy[idx]  = 0.0;  Sxz[idx]  = 0.0;
                    Syy[idx]  = 0.0;  Syz[idx]  = 0.0;  Szz[idx]  = 0.0;
                    Gamxxx[idx] = 0.0; Gamxxy[idx] = 0.0; Gamxxz[idx] = 0.0;
                    Gamxyy[idx] = 0.0; Gamxyz[idx] = 0.0; Gamxzz[idx] = 0.0;
                    Gamyxx[idx] = 0.0; Gamyxy[idx] = 0.0; Gamyxz[idx] = 0.0;
                    Gamyyy[idx] = 0.0; Gamyyz[idx] = 0.0; Gamyzz[idx] = 0.0;
                    Gamzxx[idx] = 0.0; Gamzxy[idx] = 0.0; Gamzxz[idx] = 0.0;
                    Gamzyy[idx] = 0.0; Gamzyz[idx] = 0.0; Gamzzz[idx] = 0.0;
                    Ricxx[idx] = 0.0; Ricxy[idx] = 0.0; Ricxz[idx] = 0.0;
                    Ricyy[idx] = 0.0; Ricyz[idx] = 0.0; Riczz[idx] = 0.0;
                }
            }
        }
    }
}

// ===========================================================================
// CPU symmetry BD helpers — replicate GPU symmetry boundary conditions
// Uses Cell-centered reflection: ghost i=-1 reflects to i=0, etc.
// ===========================================================================

// Apply symmetry BD to a single EM field with given parity factors
static void cpu_apply_symmetry_bd_empart(int ex[3],
                                          double *Ex, double *Ey, double *Ez,
                                          double *Bx, double *By, double *Bz,
                                          double *Kpsi, double *Kphi)
{
    int ntot = ex[0] * ex[1] * ex[2];

    // Process all indices — only ghost zones where ANY dimension < GW
    for (int k_pad = 0; k_pad < ex[2]; k_pad++) {
        for (int j_pad = 0; j_pad < ex[1]; j_pad++) {
            for (int i_pad = 0; i_pad < ex[0]; i_pad++) {
                int i = i_pad - GW;
                int j = j_pad - GW;
                int k = k_pad - GW;

                // Skip interior points (all coords >= 0)
                if (i >= 0 && j >= 0 && k >= 0) continue;
                // Skip upper ghost corners beyond ext
                if (i >= ex[0] - 2*GW || j >= ex[1] - 2*GW || k >= ex[2] - 2*GW) continue;
                // Skip if ALL coords are out of the lower ghost range
                if (i + (ex[0] - 2*GW) <= 0 || j + (ex[1] - 2*GW) <= 0 || k + (ex[2] - 2*GW) <= 0) continue;

                int out_idx = i_pad + j_pad * ex[0] + k_pad * ex[0] * ex[1];

                bool i_neg = (i < 0);
                bool j_neg = (j < 0);
                bool k_neg = (k < 0);

                // Cell-centered reflection
                int ri = i_neg ? (-i - 1) : i;
                int rj = j_neg ? (-j - 1) : j;
                int rk = k_neg ? (-k - 1) : k;
                int src_i = ri + GW;
                int src_j = rj + GW;
                int src_k = rk + GW;
                int src_idx = src_i + src_j * ex[0] + src_k * ex[0] * ex[1];

                double fi = i_neg ? -1.0 : 1.0;
                double fj = j_neg ? -1.0 : 1.0;
                double fk = k_neg ? -1.0 : 1.0;

                // EM parity: Ex=ANTI,SYM,SYM; Ey=SYM,ANTI,SYM; Ez=SYM,SYM,ANTI
                //             Bx=SYM,ANTI,ANTI; By=ANTI,SYM,ANTI; Bz=ANTI,ANTI,SYM
                //             Kpsi=SYM,SYM,SYM; Kphi=SYM,SYM,SYM
                Ex[out_idx]   = Ex[src_idx] * fi;
                Ey[out_idx]   = Ey[src_idx] * fj;
                Ez[out_idx]   = Ez[src_idx] * fk;
                Bx[out_idx]   = Bx[src_idx] * fj * fk;
                By[out_idx]   = By[src_idx] * fi * fk;
                Bz[out_idx]   = Bz[src_idx] * fi * fj;
                Kpsi[out_idx] = Kpsi[src_idx];
                Kphi[out_idx] = Kphi[src_idx];
            }
        }
    }
}

// Apply symmetry BD to a single scalar field using its SoA
static void cpu_apply_symmetry_bd_scalar(int ex[3], double *f,
                                          double SoA0, double SoA1, double SoA2)
{
    for (int k_pad = 0; k_pad < ex[2]; k_pad++) {
        for (int j_pad = 0; j_pad < ex[1]; j_pad++) {
            for (int i_pad = 0; i_pad < ex[0]; i_pad++) {
                int i = i_pad - GW;
                int j = j_pad - GW;
                int k = k_pad - GW;

                if (i >= 0 && j >= 0 && k >= 0) continue;
                if (i >= ex[0] - 2*GW || j >= ex[1] - 2*GW || k >= ex[2] - 2*GW) continue;
                if (i + (ex[0] - 2*GW) <= 0 || j + (ex[1] - 2*GW) <= 0 || k + (ex[2] - 2*GW) <= 0) continue;

                int out_idx = i_pad + j_pad * ex[0] + k_pad * ex[0] * ex[1];

                bool i_neg = (i < 0);
                bool j_neg = (j < 0);
                bool k_neg = (k < 0);

                int ri = i_neg ? (-i - 1) : i;
                int rj = j_neg ? (-j - 1) : j;
                int rk = k_neg ? (-k - 1) : k;
                int src_i = ri + GW;
                int src_j = rj + GW;
                int src_k = rk + GW;
                int src_idx = src_i + src_j * ex[0] + src_k * ex[0] * ex[1];

                double factor = (i_neg ? SoA0 : 1.0)
                              * (j_neg ? SoA1 : 1.0)
                              * (k_neg ? SoA2 : 1.0);
                f[out_idx] = f[src_idx] * factor;
            }
        }
    }
}

// Apply symmetry BD to all BSSN variables
static void cpu_apply_symmetry_bd_bssn(int ex[3],
                                        double *bssn_state[NVAR_BSSN])
{
    for (int v = 0; v < NVAR_BSSN; v++) {
        cpu_apply_symmetry_bd_scalar(ex, bssn_state[v],
                                      soa_bssn[v][0], soa_bssn[v][1], soa_bssn[v][2]);
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
    for (int k = GW; k < ex[2] - GW; k++) {
        for (int j = GW; j < ex[1] - GW; j++) {
            for (int i = GW; i < ex[0] - GW; i++) {
                int idx = i + j * ex[0] + k * ex[0] * ex[1];
                double d = fabs(gpu[idx] - cpu[idx]);
                if (d > maxdiff) maxdiff = d;
            }
        }
    }
    int fail = (maxdiff > tol) ? 1 : 0;
    printf("  max|D%-14s|= %#.2e  [%s]\n",
           name, maxdiff, fail ? "FAIL" : "PASS");
    return fail;
}

// ===========================================================================
// GPU alloc/free helpers
// ===========================================================================
static double* gpu_alloc(int n) {
    double *p;
    CUDA_CHECK(cudaMalloc(&p, n * sizeof(double)));
    return p;
}

static void gpu_copy_h2d(double *d, const double *h, int n) {
    CUDA_CHECK(cudaMemcpy(d, h, n * sizeof(double), cudaMemcpyHostToDevice));
}

static void gpu_copy_d2h(double *h, const double *d, int n) {
    CUDA_CHECK(cudaMemcpy(h, d, n * sizeof(double), cudaMemcpyDeviceToHost));
}

// ===========================================================================
// GPU path: run full Step_GPU logic on a single block
// ===========================================================================
static void run_gpu_step(int ex[3], int pex[3], int n, cudaStream_t stream,
                         // Coordinates
                         double *d_X, double *d_Y, double *d_Z,
                         // StateList pointers
                         double *d_SL_bssn[NVAR_BSSN], double *d_SL_em[NVAR_EM],
                         // SynchList_pre pointers
                         double *d_SP_bssn[NVAR_BSSN], double *d_SP_em[NVAR_EM],
                         // SynchList_cor pointers
                         double *d_SC_bssn[NVAR_BSSN], double *d_SC_em[NVAR_EM],
                         // RHSList pointers
                         double *d_RL_bssn[NVAR_BSSN], double *d_RL_em[NVAR_EM],
                         // Misc: EM source
                         double *d_Jx, double *d_Jy, double *d_Jz, double *d_qchar,
                         // Misc: stress-energy
                         double *d_rho, double *d_Sx, double *d_Sy, double *d_Sz,
                         double *d_Sxx, double *d_Sxy, double *d_Sxz,
                         double *d_Syy, double *d_Syz, double *d_Szz,
                         // Extra: Christoffel
                         double *d_Gamxxx, double *d_Gamxxy, double *d_Gamxxz,
                         double *d_Gamxyy, double *d_Gamxyz, double *d_Gamxzz,
                         double *d_Gamyxx, double *d_Gamyxy, double *d_Gamyxz,
                         double *d_Gamyyy, double *d_Gamyyz, double *d_Gamyzz,
                         double *d_Gamzxx, double *d_Gamzxy, double *d_Gamzxz,
                         double *d_Gamzyy, double *d_Gamzyz, double *d_Gamzzz,
                         // Extra: Ricci
                         double *d_Rxx, double *d_Rxy, double *d_Rxz,
                         double *d_Ryy, double *d_Ryz, double *d_Rzz,
                         // Extra: constraints
                         double *d_ham, double *d_movx, double *d_movy,
                         double *d_movz, double *d_Gmx, double *d_Gmy, double *d_Gmz)
{
    // Bbox for Sommerfeld (not used in interior, but needed by kernel)
    double bbox[6] = {
        double(GW) * dx, double(GW) * dx, double(GW) * dx,
        double(ex[0] - GW) * dx, double(ex[1] - GW) * dx, double(ex[2] - GW) * dx
    };

    double TRK4 = 0.0;
    int pre = 0, cor = 1;

    // ========================
    // Predictor (iter_count=0)
    // ========================
    {
        int iter_count = 0;

        // enforce_ga on StateList
        gpu_enforce_ga_launch(stream, ex, pex,
            d_SL_bssn[I_GXX], d_SL_bssn[I_GXY], d_SL_bssn[I_GXZ],
            d_SL_bssn[I_GYY], d_SL_bssn[I_GYZ], d_SL_bssn[I_GZZ],
            d_SL_bssn[I_AXX], d_SL_bssn[I_AXY], d_SL_bssn[I_AXZ],
            d_SL_bssn[I_AYY], d_SL_bssn[I_AYZ], d_SL_bssn[I_AZZ]);
        CUDA_CHECK(cudaDeviceSynchronize());

        // EM symmetry BD
        gpu_apply_symmetry_bd_empart_launch(stream,
            ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
            d_SL_em[I_Ex - NVAR_BSSN], d_SL_em[I_Ey - NVAR_BSSN], d_SL_em[I_Ez - NVAR_BSSN],
            d_SL_em[I_Bx - NVAR_BSSN], d_SL_em[I_By - NVAR_BSSN], d_SL_em[I_Bz - NVAR_BSSN],
            d_SL_em[I_KPSI - NVAR_BSSN], d_SL_em[I_KPHI - NVAR_BSSN]);
        CUDA_CHECK(cudaDeviceSynchronize());

        // EM RHS (predictor: state=StateList, RHS->RHSList)
        gpu_compute_rhs_empart_launch(stream, ex, pex,
            d_X, d_Y, d_Z,
            d_SL_bssn[I_CHI],
            d_Jx, d_Jy, d_Jz, d_qchar,
            d_SL_bssn[I_GXX], d_SL_bssn[I_GXY], d_SL_bssn[I_GXZ],
            d_SL_bssn[I_GYY], d_SL_bssn[I_GYZ], d_SL_bssn[I_GZZ],
            d_SL_bssn[I_LAP],
            d_SL_bssn[I_BETAX], d_SL_bssn[I_BETAY], d_SL_bssn[I_BETAZ],
            d_SL_bssn[I_TRK],
            d_SL_em[I_Ex - NVAR_BSSN], d_SL_em[I_Ey - NVAR_BSSN], d_SL_em[I_Ez - NVAR_BSSN],
            d_SL_em[I_Bx - NVAR_BSSN], d_SL_em[I_By - NVAR_BSSN], d_SL_em[I_Bz - NVAR_BSSN],
            d_SL_em[I_KPSI - NVAR_BSSN], d_SL_em[I_KPHI - NVAR_BSSN],
            d_RL_em[I_Ex - NVAR_BSSN], d_RL_em[I_Ey - NVAR_BSSN], d_RL_em[I_Ez - NVAR_BSSN],
            d_RL_em[I_Bx - NVAR_BSSN], d_RL_em[I_By - NVAR_BSSN], d_RL_em[I_Bz - NVAR_BSSN],
            d_RL_em[I_KPSI - NVAR_BSSN], d_RL_em[I_KPHI - NVAR_BSSN],
            d_rho, d_Sx, d_Sy, d_Sz,
            d_Sxx, d_Sxy, d_Sxz,
            d_Syy, d_Syz, d_Szz,
            0 /*symmetry*/, 0 /*lev*/, 0.0 /*eps*/);
        CUDA_CHECK(cudaDeviceSynchronize());

        // BSSN symmetry BD
        gpu_apply_symmetry_bd_fused_launch(stream,
            ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
            d_SL_bssn[I_CHI], d_SL_bssn[I_TRK],
            d_SL_bssn[I_GXX], d_SL_bssn[I_GXY], d_SL_bssn[I_GXZ],
            d_SL_bssn[I_GYY], d_SL_bssn[I_GYZ], d_SL_bssn[I_GZZ],
            d_SL_bssn[I_AXX], d_SL_bssn[I_AXY], d_SL_bssn[I_AXZ],
            d_SL_bssn[I_AYY], d_SL_bssn[I_AYZ], d_SL_bssn[I_AZZ],
            d_SL_bssn[I_GAMX], d_SL_bssn[I_GAMY], d_SL_bssn[I_GAMZ],
            d_SL_bssn[I_LAP],
            d_SL_bssn[I_BETAX], d_SL_bssn[I_BETAY], d_SL_bssn[I_BETAZ],
            d_SL_bssn[I_DTSFX], d_SL_bssn[I_DTSFY], d_SL_bssn[I_DTSFZ]);
        CUDA_CHECK(cudaDeviceSynchronize());

        // BSSN RHS (predictor: state=StateList, RHS->RHSList)
        gpu_compute_rhs_bssn_launch(stream, ex, pex, TRK4,
            d_X, d_Y, d_Z,
            d_SL_bssn[I_CHI], d_SL_bssn[I_TRK],
            d_SL_bssn[I_GXX], d_SL_bssn[I_GXY], d_SL_bssn[I_GXZ],
            d_SL_bssn[I_GYY], d_SL_bssn[I_GYZ], d_SL_bssn[I_GZZ],
            d_SL_bssn[I_AXX], d_SL_bssn[I_AXY], d_SL_bssn[I_AXZ],
            d_SL_bssn[I_AYY], d_SL_bssn[I_AYZ], d_SL_bssn[I_AZZ],
            d_SL_bssn[I_GAMX], d_SL_bssn[I_GAMY], d_SL_bssn[I_GAMZ],
            d_SL_bssn[I_LAP],
            d_SL_bssn[I_BETAX], d_SL_bssn[I_BETAY], d_SL_bssn[I_BETAZ],
            d_SL_bssn[I_DTSFX], d_SL_bssn[I_DTSFY], d_SL_bssn[I_DTSFZ],
            d_RL_bssn[I_CHI], d_RL_bssn[I_TRK],
            d_RL_bssn[I_GXX], d_RL_bssn[I_GXY], d_RL_bssn[I_GXZ],
            d_RL_bssn[I_GYY], d_RL_bssn[I_GYZ], d_RL_bssn[I_GZZ],
            d_RL_bssn[I_AXX], d_RL_bssn[I_AXY], d_RL_bssn[I_AXZ],
            d_RL_bssn[I_AYY], d_RL_bssn[I_AYZ], d_RL_bssn[I_AZZ],
            d_RL_bssn[I_GAMX], d_RL_bssn[I_GAMY], d_RL_bssn[I_GAMZ],
            d_RL_bssn[I_LAP],
            d_RL_bssn[I_BETAX], d_RL_bssn[I_BETAY], d_RL_bssn[I_BETAZ],
            d_RL_bssn[I_DTSFX], d_RL_bssn[I_DTSFY], d_RL_bssn[I_DTSFZ],
            d_rho, d_Sx, d_Sy, d_Sz,
            d_Sxx, d_Sxy, d_Sxz,
            d_Syy, d_Syz, d_Szz,
            d_Gamxxx, d_Gamxxy, d_Gamxxz, d_Gamxyy, d_Gamxyz, d_Gamxzz,
            d_Gamyxx, d_Gamyxy, d_Gamyxz, d_Gamyyy, d_Gamyyz, d_Gamyzz,
            d_Gamzxx, d_Gamzxy, d_Gamzxz, d_Gamzyy, d_Gamzyz, d_Gamzzz,
            d_Rxx, d_Rxy, d_Rxz, d_Ryy, d_Ryz, d_Rzz,
            d_ham, d_movx, d_movy, d_movz, d_Gmx, d_Gmy, d_Gmz,
            0 /*symmetry*/, 0 /*lev*/, 0.0 /*eps*/, pre);
        CUDA_CHECK(cudaDeviceSynchronize());

        // Fused RK4 (predictor): f0=StateList, f1=SynchList_pre, f_rhs=RHSList
        // Fused RK4 (predictor): batch 1 — BSSN vars
        {
            FusedRK4Vars rk4_vars;
            rk4_vars.num_vars = 0;
            for (int v = 0; v < NVAR_BSSN; v++) {
                rk4_vars.f0[rk4_vars.num_vars] = d_SL_bssn[v];
                rk4_vars.f1[rk4_vars.num_vars] = d_SP_bssn[v];
                rk4_vars.f_rhs[rk4_vars.num_vars] = d_RL_bssn[v];
                rk4_vars.num_vars++;
            }
            gpu_fused_rungekutta4_rout_launch(stream, ex, pex, dT, rk4_vars, iter_count);
        }
        // Fused RK4 (predictor): batch 2 — EM vars
        {
            FusedRK4Vars rk4_vars;
            rk4_vars.num_vars = 0;
            for (int v = 0; v < NVAR_EM; v++) {
                rk4_vars.f0[rk4_vars.num_vars] = d_SL_em[v];
                rk4_vars.f1[rk4_vars.num_vars] = d_SP_em[v];
                rk4_vars.f_rhs[rk4_vars.num_vars] = d_RL_em[v];
                rk4_vars.num_vars++;
            }
            gpu_fused_rungekutta4_rout_launch(stream, ex, pex, dT, rk4_vars, iter_count);
            CUDA_CHECK(cudaDeviceSynchronize());
        }

        // Fused Sommerfeld (predictor, lev>0): f0=StateList, f_out=SynchList_pre
        {
            FusedSommerfeldVars somm_vars;
            somm_vars.num_vars = 0;
            for (int v = 0; v < NVAR_BSSN; v++) {
                somm_vars.f0[somm_vars.num_vars] = d_SL_bssn[v];
                somm_vars.f_out[somm_vars.num_vars] = d_SP_bssn[v];
                somm_vars.propspeed[somm_vars.num_vars] = propspeed_bssn[v];
                somm_vars.SoA0[somm_vars.num_vars] = soa_bssn[v][0];
                somm_vars.SoA1[somm_vars.num_vars] = soa_bssn[v][1];
                somm_vars.SoA2[somm_vars.num_vars] = soa_bssn[v][2];
                somm_vars.num_vars++;
            }
            gpu_fused_sommerfeld_rout_launch(stream, ex, pex, d_X, d_Y, d_Z,
                bbox[0], bbox[1], bbox[2], bbox[3], bbox[4], bbox[5],
                dT, d_SL_bssn[I_CHI], d_SL_bssn[I_LAP], somm_vars, 0, 1);
        }
        {
            FusedSommerfeldVars somm_vars;
            somm_vars.num_vars = 0;
            for (int v = 0; v < NVAR_EM; v++) {
                somm_vars.f0[somm_vars.num_vars] = d_SL_em[v];
                somm_vars.f_out[somm_vars.num_vars] = d_SP_em[v];
                somm_vars.propspeed[somm_vars.num_vars] = propspeed_em[v];
                somm_vars.SoA0[somm_vars.num_vars] = soa_em[v][0];
                somm_vars.SoA1[somm_vars.num_vars] = soa_em[v][1];
                somm_vars.SoA2[somm_vars.num_vars] = soa_em[v][2];
                somm_vars.num_vars++;
            }
            gpu_fused_sommerfeld_rout_launch(stream, ex, pex,
                d_X, d_Y, d_Z,
                bbox[0], bbox[1], bbox[2], bbox[3], bbox[4], bbox[5],
                dT, d_SL_bssn[I_CHI], d_SL_bssn[I_LAP],
                somm_vars, 0 /*Symmetry*/, cor);
            CUDA_CHECK(cudaDeviceSynchronize());
        }

        // lowerboundset on SynchList_pre phi
        gpu_lowerboundset_launch(stream, ex, pex, d_SP_bssn[I_CHI], chitiny);
        CUDA_CHECK(cudaDeviceSynchronize());
    }

    // ============================================================
    // Correctors (iter_count = 1, 2, 3)
    // ============================================================
    for (int iter_count = 1; iter_count < 4; iter_count++) {
        if (iter_count == 1 || iter_count == 3) TRK4 += dT / 2.0;

        // enforce_ga on SynchList_pre
        gpu_enforce_ga_launch(stream, ex, pex,
            d_SP_bssn[I_GXX], d_SP_bssn[I_GXY], d_SP_bssn[I_GXZ],
            d_SP_bssn[I_GYY], d_SP_bssn[I_GYZ], d_SP_bssn[I_GZZ],
            d_SP_bssn[I_AXX], d_SP_bssn[I_AXY], d_SP_bssn[I_AXZ],
            d_SP_bssn[I_AYY], d_SP_bssn[I_AYZ], d_SP_bssn[I_AZZ]);
        CUDA_CHECK(cudaDeviceSynchronize());

        // EM symmetry BD on SynchList_pre EM vars
        gpu_apply_symmetry_bd_empart_launch(stream,
            ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
            d_SP_em[I_Ex - NVAR_BSSN], d_SP_em[I_Ey - NVAR_BSSN], d_SP_em[I_Ez - NVAR_BSSN],
            d_SP_em[I_Bx - NVAR_BSSN], d_SP_em[I_By - NVAR_BSSN], d_SP_em[I_Bz - NVAR_BSSN],
            d_SP_em[I_KPSI - NVAR_BSSN], d_SP_em[I_KPHI - NVAR_BSSN]);
        CUDA_CHECK(cudaDeviceSynchronize());

        // EM RHS (corrector: state=SynchList_pre, RHS->SynchList_cor)
        gpu_compute_rhs_empart_launch(stream, ex, pex,
            d_X, d_Y, d_Z,
            d_SP_bssn[I_CHI],
            d_Jx, d_Jy, d_Jz, d_qchar,
            d_SP_bssn[I_GXX], d_SP_bssn[I_GXY], d_SP_bssn[I_GXZ],
            d_SP_bssn[I_GYY], d_SP_bssn[I_GYZ], d_SP_bssn[I_GZZ],
            d_SP_bssn[I_LAP],
            d_SP_bssn[I_BETAX], d_SP_bssn[I_BETAY], d_SP_bssn[I_BETAZ],
            d_SP_bssn[I_TRK],
            d_SP_em[I_Ex - NVAR_BSSN], d_SP_em[I_Ey - NVAR_BSSN], d_SP_em[I_Ez - NVAR_BSSN],
            d_SP_em[I_Bx - NVAR_BSSN], d_SP_em[I_By - NVAR_BSSN], d_SP_em[I_Bz - NVAR_BSSN],
            d_SP_em[I_KPSI - NVAR_BSSN], d_SP_em[I_KPHI - NVAR_BSSN],
            d_SC_em[I_Ex - NVAR_BSSN], d_SC_em[I_Ey - NVAR_BSSN], d_SC_em[I_Ez - NVAR_BSSN],
            d_SC_em[I_Bx - NVAR_BSSN], d_SC_em[I_By - NVAR_BSSN], d_SC_em[I_Bz - NVAR_BSSN],
            d_SC_em[I_KPSI - NVAR_BSSN], d_SC_em[I_KPHI - NVAR_BSSN],
            d_rho, d_Sx, d_Sy, d_Sz,
            d_Sxx, d_Sxy, d_Sxz,
            d_Syy, d_Syz, d_Szz,
            0 /*symmetry*/, 0 /*lev*/, 0.0 /*eps*/);
        CUDA_CHECK(cudaDeviceSynchronize());

        // BSSN symmetry BD on SynchList_pre BSSN vars
        gpu_apply_symmetry_bd_fused_launch(stream,
            ex[0], ex[1], ex[2], pex[0], pex[1], pex[2],
            d_SP_bssn[I_CHI], d_SP_bssn[I_TRK],
            d_SP_bssn[I_GXX], d_SP_bssn[I_GXY], d_SP_bssn[I_GXZ],
            d_SP_bssn[I_GYY], d_SP_bssn[I_GYZ], d_SP_bssn[I_GZZ],
            d_SP_bssn[I_AXX], d_SP_bssn[I_AXY], d_SP_bssn[I_AXZ],
            d_SP_bssn[I_AYY], d_SP_bssn[I_AYZ], d_SP_bssn[I_AZZ],
            d_SP_bssn[I_GAMX], d_SP_bssn[I_GAMY], d_SP_bssn[I_GAMZ],
            d_SP_bssn[I_LAP],
            d_SP_bssn[I_BETAX], d_SP_bssn[I_BETAY], d_SP_bssn[I_BETAZ],
            d_SP_bssn[I_DTSFX], d_SP_bssn[I_DTSFY], d_SP_bssn[I_DTSFZ]);
        CUDA_CHECK(cudaDeviceSynchronize());

        // BSSN RHS (corrector: state=SynchList_pre, RHS->SynchList_cor)
        gpu_compute_rhs_bssn_launch(stream, ex, pex, TRK4,
            d_X, d_Y, d_Z,
            d_SP_bssn[I_CHI], d_SP_bssn[I_TRK],
            d_SP_bssn[I_GXX], d_SP_bssn[I_GXY], d_SP_bssn[I_GXZ],
            d_SP_bssn[I_GYY], d_SP_bssn[I_GYZ], d_SP_bssn[I_GZZ],
            d_SP_bssn[I_AXX], d_SP_bssn[I_AXY], d_SP_bssn[I_AXZ],
            d_SP_bssn[I_AYY], d_SP_bssn[I_AYZ], d_SP_bssn[I_AZZ],
            d_SP_bssn[I_GAMX], d_SP_bssn[I_GAMY], d_SP_bssn[I_GAMZ],
            d_SP_bssn[I_LAP],
            d_SP_bssn[I_BETAX], d_SP_bssn[I_BETAY], d_SP_bssn[I_BETAZ],
            d_SP_bssn[I_DTSFX], d_SP_bssn[I_DTSFY], d_SP_bssn[I_DTSFZ],
            d_SC_bssn[I_CHI], d_SC_bssn[I_TRK],
            d_SC_bssn[I_GXX], d_SC_bssn[I_GXY], d_SC_bssn[I_GXZ],
            d_SC_bssn[I_GYY], d_SC_bssn[I_GYZ], d_SC_bssn[I_GZZ],
            d_SC_bssn[I_AXX], d_SC_bssn[I_AXY], d_SC_bssn[I_AXZ],
            d_SC_bssn[I_AYY], d_SC_bssn[I_AYZ], d_SC_bssn[I_AZZ],
            d_SC_bssn[I_GAMX], d_SC_bssn[I_GAMY], d_SC_bssn[I_GAMZ],
            d_SC_bssn[I_LAP],
            d_SC_bssn[I_BETAX], d_SC_bssn[I_BETAY], d_SC_bssn[I_BETAZ],
            d_SC_bssn[I_DTSFX], d_SC_bssn[I_DTSFY], d_SC_bssn[I_DTSFZ],
            d_rho, d_Sx, d_Sy, d_Sz,
            d_Sxx, d_Sxy, d_Sxz,
            d_Syy, d_Syz, d_Szz,
            d_Gamxxx, d_Gamxxy, d_Gamxxz, d_Gamxyy, d_Gamxyz, d_Gamxzz,
            d_Gamyxx, d_Gamyxy, d_Gamyxz, d_Gamyyy, d_Gamyyz, d_Gamyzz,
            d_Gamzxx, d_Gamzxy, d_Gamzxz, d_Gamzyy, d_Gamzyz, d_Gamzzz,
            d_Rxx, d_Rxy, d_Rxz, d_Ryy, d_Ryz, d_Rzz,
            d_ham, d_movx, d_movy, d_movz, d_Gmx, d_Gmy, d_Gmz,
            0 /*symmetry*/, 0 /*lev*/, 0.0 /*eps*/, cor);
        CUDA_CHECK(cudaDeviceSynchronize());

        // Fused RK4 (corrector): f0=StateList, f1=SynchList_cor, f_rhs=RHSList
        {
            FusedRK4Vars rk4_vars;
            rk4_vars.num_vars = 0;
            for (int v = 0; v < NVAR_BSSN; v++) {
                rk4_vars.f0[rk4_vars.num_vars] = d_SL_bssn[v];
                rk4_vars.f1[rk4_vars.num_vars] = d_SC_bssn[v];
                rk4_vars.f_rhs[rk4_vars.num_vars] = d_RL_bssn[v];
                rk4_vars.num_vars++;
            }
            gpu_fused_rungekutta4_rout_launch(stream, ex, pex, dT, rk4_vars, iter_count);
        }
        {
            FusedRK4Vars rk4_vars;
            rk4_vars.num_vars = 0;
            for (int v = 0; v < NVAR_EM; v++) {
                rk4_vars.f0[rk4_vars.num_vars] = d_SL_em[v];
                rk4_vars.f1[rk4_vars.num_vars] = d_SC_em[v];
                rk4_vars.f_rhs[rk4_vars.num_vars] = d_RL_em[v];
                rk4_vars.num_vars++;
            }
            gpu_fused_rungekutta4_rout_launch(stream, ex, pex, dT, rk4_vars, iter_count);
            CUDA_CHECK(cudaDeviceSynchronize());
        }

        // Fused Sommerfeld (corrector, lev>0): f0=StateList, f_out=SynchList_cor
        {
            FusedSommerfeldVars somm_vars;
            somm_vars.num_vars = 0;
            for (int v = 0; v < NVAR_BSSN; v++) {
                somm_vars.f0[somm_vars.num_vars] = d_SL_bssn[v];
                somm_vars.f_out[somm_vars.num_vars] = d_SC_bssn[v];
                somm_vars.propspeed[somm_vars.num_vars] = propspeed_bssn[v];
                somm_vars.SoA0[somm_vars.num_vars] = soa_bssn[v][0];
                somm_vars.SoA1[somm_vars.num_vars] = soa_bssn[v][1];
                somm_vars.SoA2[somm_vars.num_vars] = soa_bssn[v][2];
                somm_vars.num_vars++;
            }
            for (int v = 0; v < NVAR_EM; v++) {
                somm_vars.f0[somm_vars.num_vars] = d_SL_em[v];
                somm_vars.f_out[somm_vars.num_vars] = d_SC_em[v];
                somm_vars.propspeed[somm_vars.num_vars] = propspeed_em[v];
                somm_vars.SoA0[somm_vars.num_vars] = soa_em[v][0];
                somm_vars.SoA1[somm_vars.num_vars] = soa_em[v][1];
                somm_vars.SoA2[somm_vars.num_vars] = soa_em[v][2];
                somm_vars.num_vars++;
            }
            gpu_fused_sommerfeld_rout_launch(stream, ex, pex,
                d_X, d_Y, d_Z,
                bbox[0], bbox[1], bbox[2], bbox[3], bbox[4], bbox[5],
                dT, d_SL_bssn[I_CHI], d_SL_bssn[I_LAP],
                somm_vars, 0 /*Symmetry*/, cor);
            CUDA_CHECK(cudaDeviceSynchronize());
        }

        // lowerboundset on SynchList_cor phi
        gpu_lowerboundset_launch(stream, ex, pex, d_SC_bssn[I_CHI], chitiny);
        CUDA_CHECK(cudaDeviceSynchronize());

        // swapList: swap SynchList_pre <-> SynchList_cor
        // On GPU, this means swapping the d_xxx pointers in fgfs slots.
        // For our flat arrays, we swap the data by copying.
        // Actually, in the real code swapList swaps stored pointers.
        // Here we swap the array contents.
        {
            double *tmp;
            for (int v = 0; v < NVAR_BSSN; v++) {
                tmp = d_SP_bssn[v]; d_SP_bssn[v] = d_SC_bssn[v]; d_SC_bssn[v] = tmp;
            }
            for (int v = 0; v < NVAR_EM; v++) {
                tmp = d_SP_em[v]; d_SP_em[v] = d_SC_em[v]; d_SC_em[v] = tmp;
            }
        }
    }

    // Final swap: StateList <-> SynchList_cor (after last swapList, SynchList_cor
    // holds the final result; in the real code it's StateList<-SynchList_cor and
    // OldStateList<-SynchList_cor)
    // The final result after corrector 3 swapList is in SynchList_pre.
    // So we swap SL <-> SP.
    {
        double *tmp;
        for (int v = 0; v < NVAR_BSSN; v++) {
            tmp = d_SL_bssn[v]; d_SL_bssn[v] = d_SP_bssn[v]; d_SP_bssn[v] = tmp;
        }
        for (int v = 0; v < NVAR_EM; v++) {
            tmp = d_SL_em[v]; d_SL_em[v] = d_SP_em[v]; d_SP_em[v] = tmp;
        }
    }
}

// ===========================================================================
// CPU path: run equivalent full step using Fortran functions
// ===========================================================================
static void run_cpu_step(int ex[3], int pex[3],
                         double *X, double *Y, double *Z,
                         // StateList
                         double *SL_bssn[NVAR_BSSN], double *SL_em[NVAR_EM],
                         // SynchList_pre
                         double *SP_bssn[NVAR_BSSN], double *SP_em[NVAR_EM],
                         // SynchList_cor
                         double *SC_bssn[NVAR_BSSN], double *SC_em[NVAR_EM],
                         // RHSList
                         double *RL_bssn[NVAR_BSSN], double *RL_em[NVAR_EM],
                         // EM source
                         double *Jx, double *Jy, double *Jz, double *qchar,
                         // stress-energy
                         double *rho, double *Sx, double *Sy, double *Sz,
                         double *Sxx, double *Sxy, double *Sxz,
                         double *Syy, double *Syz, double *Szz,
                         // Christoffel
                         double *Gamxxx, double *Gamxxy, double *Gamxxz,
                         double *Gamxyy, double *Gamxyz, double *Gamxzz,
                         double *Gamyxx, double *Gamyxy, double *Gamyxz,
                         double *Gamyyy, double *Gamyyz, double *Gamyzz,
                         double *Gamzxx, double *Gamzxy, double *Gamzxz,
                         double *Gamzyy, double *Gamzyz, double *Gamzzz,
                         // Ricci
                         double *Ricxx, double *Ricxy, double *Ricxz,
                         double *Ricyy, double *Ricyz, double *Riczz,
                         // constraints
                         double *ham, double *movx, double *movy, double *movz,
                         double *Gmx, double *Gmy, double *Gmz)
{
    double bbox[6] = {
        double(GW) * dx, double(GW) * dx, double(GW) * dx,
        double(ex[0] - GW) * dx, double(ex[1] - GW) * dx, double(ex[2] - GW) * dx
    };

    double TRK4 = 0.0;
    int pre = 0, cor = 1;

    // ========================
    // Predictor (iter_count=0)
    // ========================
    {
        int iter_count = 0;

        // enforce_ga on StateList
        f_enforce_ga(ex,
            SL_bssn[I_GXX], SL_bssn[I_GXY], SL_bssn[I_GXZ],
            SL_bssn[I_GYY], SL_bssn[I_GYZ], SL_bssn[I_GZZ],
            SL_bssn[I_AXX], SL_bssn[I_AXY], SL_bssn[I_AXZ],
            SL_bssn[I_AYY], SL_bssn[I_AYZ], SL_bssn[I_AZZ]);

        // EM symmetry BD
        cpu_apply_symmetry_bd_empart(ex,
            SL_em[I_Ex - NVAR_BSSN], SL_em[I_Ey - NVAR_BSSN], SL_em[I_Ez - NVAR_BSSN],
            SL_em[I_Bx - NVAR_BSSN], SL_em[I_By - NVAR_BSSN], SL_em[I_Bz - NVAR_BSSN],
            SL_em[I_KPSI - NVAR_BSSN], SL_em[I_KPHI - NVAR_BSSN]);

        // EM RHS (predictor): state=StateList, RHS->RHSList
        {
            int sym_ref = 0, lev_ref = 0;
            double eps_ref = 0.0;
            f_compute_rhs_empart(ex, X, Y, Z,
                SL_bssn[I_CHI],
                Jx, Jy, Jz, qchar,
                SL_bssn[I_GXX], SL_bssn[I_GXY], SL_bssn[I_GXZ],
                SL_bssn[I_GYY], SL_bssn[I_GYZ], SL_bssn[I_GZZ],
                SL_bssn[I_LAP],
                SL_bssn[I_BETAX], SL_bssn[I_BETAY], SL_bssn[I_BETAZ],
                SL_bssn[I_TRK],
                SL_em[I_Ex - NVAR_BSSN], SL_em[I_Ey - NVAR_BSSN], SL_em[I_Ez - NVAR_BSSN],
                SL_em[I_Bx - NVAR_BSSN], SL_em[I_By - NVAR_BSSN], SL_em[I_Bz - NVAR_BSSN],
                SL_em[I_KPSI - NVAR_BSSN], SL_em[I_KPHI - NVAR_BSSN],
                RL_em[I_Ex - NVAR_BSSN], RL_em[I_Ey - NVAR_BSSN], RL_em[I_Ez - NVAR_BSSN],
                RL_em[I_Bx - NVAR_BSSN], RL_em[I_By - NVAR_BSSN], RL_em[I_Bz - NVAR_BSSN],
                RL_em[I_KPSI - NVAR_BSSN], RL_em[I_KPHI - NVAR_BSSN],
                rho, Sx, Sy, Sz,
                Sxx, Sxy, Sxz,
                Syy, Syz, Szz,
                sym_ref, lev_ref, eps_ref);
        }

        // BSSN symmetry BD
        cpu_apply_symmetry_bd_bssn(ex, SL_bssn);

        // BSSN RHS (predictor): state=StateList, RHS->RHSList
        {
            int sym_ref = 0, lev_ref = 0;
            double eps_ref = 0.0;
            f_compute_rhs_bssn(ex, TRK4, X, Y, Z,
                SL_bssn[I_CHI], SL_bssn[I_TRK],
                SL_bssn[I_GXX], SL_bssn[I_GXY], SL_bssn[I_GXZ],
                SL_bssn[I_GYY], SL_bssn[I_GYZ], SL_bssn[I_GZZ],
                SL_bssn[I_AXX], SL_bssn[I_AXY], SL_bssn[I_AXZ],
                SL_bssn[I_AYY], SL_bssn[I_AYZ], SL_bssn[I_AZZ],
                SL_bssn[I_GAMX], SL_bssn[I_GAMY], SL_bssn[I_GAMZ],
                SL_bssn[I_LAP],
                SL_bssn[I_BETAX], SL_bssn[I_BETAY], SL_bssn[I_BETAZ],
                SL_bssn[I_DTSFX], SL_bssn[I_DTSFY], SL_bssn[I_DTSFZ],
                RL_bssn[I_CHI], RL_bssn[I_TRK],
                RL_bssn[I_GXX], RL_bssn[I_GXY], RL_bssn[I_GXZ],
                RL_bssn[I_GYY], RL_bssn[I_GYZ], RL_bssn[I_GZZ],
                RL_bssn[I_AXX], RL_bssn[I_AXY], RL_bssn[I_AXZ],
                RL_bssn[I_AYY], RL_bssn[I_AYZ], RL_bssn[I_AZZ],
                RL_bssn[I_GAMX], RL_bssn[I_GAMY], RL_bssn[I_GAMZ],
                RL_bssn[I_LAP],
                RL_bssn[I_BETAX], RL_bssn[I_BETAY], RL_bssn[I_BETAZ],
                RL_bssn[I_DTSFX], RL_bssn[I_DTSFY], RL_bssn[I_DTSFZ],
                rho, Sx, Sy, Sz,
                Sxx, Sxy, Sxz,
                Syy, Syz, Szz,
                Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
                Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
                Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
                Ricxx, Ricxy, Ricxz, Ricyy, Ricyz, Riczz,
                ham, movx, movy, movz, Gmx, Gmy, Gmz,
                sym_ref, lev_ref, eps_ref, pre);
        }

        // Per-variable RK4 (predictor): f0=StateList, f1=SynchList_pre, f_rhs=RHSList
        for (int v = 0; v < NVAR_BSSN; v++) {
            f_rungekutta4_rout(ex, dT, SL_bssn[v], SP_bssn[v], RL_bssn[v], iter_count);
        }
        for (int v = 0; v < NVAR_EM; v++) {
            f_rungekutta4_rout(ex, dT, SL_em[v], SP_em[v], RL_em[v], iter_count);
        }

        // Per-variable Sommerfeld (predictor, lev>0): f0=StateList, f=SynchList_pre
        {
            int Symmetry = 0;
            for (int v = 0; v < NVAR_BSSN; v++) {
                f_sommerfeld_rout(ex, X, Y, Z,
                    bbox[0], bbox[1], bbox[2], bbox[3], bbox[4], bbox[5],
                    dT, SL_bssn[I_CHI], SL_bssn[I_LAP],
                    SL_bssn[v], SP_bssn[v],
                    soa_bssn[v],
                    Symmetry, cor);
            }
            for (int v = 0; v < NVAR_EM; v++) {
                f_sommerfeld_rout(ex, X, Y, Z,
                    bbox[0], bbox[1], bbox[2], bbox[3], bbox[4], bbox[5],
                    dT, SL_bssn[I_CHI], SL_bssn[I_LAP],
                    SL_em[v], SP_em[v],
                    soa_em[v],
                    Symmetry, cor);
            }
        }

        // lowerboundset on SynchList_pre phi
        f_lowerboundset(ex, SP_bssn[I_CHI], chitiny_val);
    }

    // ============================================================
    // Correctors (iter_count = 1, 2, 3)
    // ============================================================
    for (int iter_count = 1; iter_count < 4; iter_count++) {
        if (iter_count == 1 || iter_count == 3) TRK4 += dT / 2.0;

        // enforce_ga on SynchList_pre
        f_enforce_ga(ex,
            SP_bssn[I_GXX], SP_bssn[I_GXY], SP_bssn[I_GXZ],
            SP_bssn[I_GYY], SP_bssn[I_GYZ], SP_bssn[I_GZZ],
            SP_bssn[I_AXX], SP_bssn[I_AXY], SP_bssn[I_AXZ],
            SP_bssn[I_AYY], SP_bssn[I_AYZ], SP_bssn[I_AZZ]);

        // EM symmetry BD on SynchList_pre EM vars
        cpu_apply_symmetry_bd_empart(ex,
            SP_em[I_Ex - NVAR_BSSN], SP_em[I_Ey - NVAR_BSSN], SP_em[I_Ez - NVAR_BSSN],
            SP_em[I_Bx - NVAR_BSSN], SP_em[I_By - NVAR_BSSN], SP_em[I_Bz - NVAR_BSSN],
            SP_em[I_KPSI - NVAR_BSSN], SP_em[I_KPHI - NVAR_BSSN]);

        // EM RHS (corrector): state=SynchList_pre, RHS->SynchList_cor
        {
            int sym_ref = 0, lev_ref = 0;
            double eps_ref = 0.0;
            f_compute_rhs_empart(ex, X, Y, Z,
                SP_bssn[I_CHI],
                Jx, Jy, Jz, qchar,
                SP_bssn[I_GXX], SP_bssn[I_GXY], SP_bssn[I_GXZ],
                SP_bssn[I_GYY], SP_bssn[I_GYZ], SP_bssn[I_GZZ],
                SP_bssn[I_LAP],
                SP_bssn[I_BETAX], SP_bssn[I_BETAY], SP_bssn[I_BETAZ],
                SP_bssn[I_TRK],
                SP_em[I_Ex - NVAR_BSSN], SP_em[I_Ey - NVAR_BSSN], SP_em[I_Ez - NVAR_BSSN],
                SP_em[I_Bx - NVAR_BSSN], SP_em[I_By - NVAR_BSSN], SP_em[I_Bz - NVAR_BSSN],
                SP_em[I_KPSI - NVAR_BSSN], SP_em[I_KPHI - NVAR_BSSN],
                SC_em[I_Ex - NVAR_BSSN], SC_em[I_Ey - NVAR_BSSN], SC_em[I_Ez - NVAR_BSSN],
                SC_em[I_Bx - NVAR_BSSN], SC_em[I_By - NVAR_BSSN], SC_em[I_Bz - NVAR_BSSN],
                SC_em[I_KPSI - NVAR_BSSN], SC_em[I_KPHI - NVAR_BSSN],
                rho, Sx, Sy, Sz,
                Sxx, Sxy, Sxz,
                Syy, Syz, Szz,
                sym_ref, lev_ref, eps_ref);
        }

        // BSSN symmetry BD on SynchList_pre BSSN vars
        cpu_apply_symmetry_bd_bssn(ex, SP_bssn);

        // BSSN RHS (corrector): state=SynchList_pre, RHS->SynchList_cor
        {
            int sym_ref = 0, lev_ref = 0;
            double eps_ref = 0.0;
            f_compute_rhs_bssn(ex, TRK4, X, Y, Z,
                SP_bssn[I_CHI], SP_bssn[I_TRK],
                SP_bssn[I_GXX], SP_bssn[I_GXY], SP_bssn[I_GXZ],
                SP_bssn[I_GYY], SP_bssn[I_GYZ], SP_bssn[I_GZZ],
                SP_bssn[I_AXX], SP_bssn[I_AXY], SP_bssn[I_AXZ],
                SP_bssn[I_AYY], SP_bssn[I_AYZ], SP_bssn[I_AZZ],
                SP_bssn[I_GAMX], SP_bssn[I_GAMY], SP_bssn[I_GAMZ],
                SP_bssn[I_LAP],
                SP_bssn[I_BETAX], SP_bssn[I_BETAY], SP_bssn[I_BETAZ],
                SP_bssn[I_DTSFX], SP_bssn[I_DTSFY], SP_bssn[I_DTSFZ],
                SC_bssn[I_CHI], SC_bssn[I_TRK],
                SC_bssn[I_GXX], SC_bssn[I_GXY], SC_bssn[I_GXZ],
                SC_bssn[I_GYY], SC_bssn[I_GYZ], SC_bssn[I_GZZ],
                SC_bssn[I_AXX], SC_bssn[I_AXY], SC_bssn[I_AXZ],
                SC_bssn[I_AYY], SC_bssn[I_AYZ], SC_bssn[I_AZZ],
                SC_bssn[I_GAMX], SC_bssn[I_GAMY], SC_bssn[I_GAMZ],
                SC_bssn[I_LAP],
                SC_bssn[I_BETAX], SC_bssn[I_BETAY], SC_bssn[I_BETAZ],
                SC_bssn[I_DTSFX], SC_bssn[I_DTSFY], SC_bssn[I_DTSFZ],
                rho, Sx, Sy, Sz,
                Sxx, Sxy, Sxz,
                Syy, Syz, Szz,
                Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
                Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
                Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
                Ricxx, Ricxy, Ricxz, Ricyy, Ricyz, Riczz,
                ham, movx, movy, movz, Gmx, Gmy, Gmz,
                sym_ref, lev_ref, eps_ref, cor);
        }

        // Per-variable RK4 (corrector): f0=StateList, f1=SynchList_cor, f_rhs=RHSList
        for (int v = 0; v < NVAR_BSSN; v++) {
            f_rungekutta4_rout(ex, dT, SL_bssn[v], SC_bssn[v], RL_bssn[v], iter_count);
        }
        for (int v = 0; v < NVAR_EM; v++) {
            f_rungekutta4_rout(ex, dT, SL_em[v], SC_em[v], RL_em[v], iter_count);
        }

        // Per-variable Sommerfeld (corrector, lev>0): f0=StateList, f=SynchList_cor
        {
            int Symmetry = 0;
            for (int v = 0; v < NVAR_BSSN; v++) {
                f_sommerfeld_rout(ex, X, Y, Z,
                    bbox[0], bbox[1], bbox[2], bbox[3], bbox[4], bbox[5],
                    dT, SL_bssn[I_CHI], SL_bssn[I_LAP],
                    SL_bssn[v], SC_bssn[v],
                    soa_bssn[v],
                    Symmetry, cor);
            }
            for (int v = 0; v < NVAR_EM; v++) {
                f_sommerfeld_rout(ex, X, Y, Z,
                    bbox[0], bbox[1], bbox[2], bbox[3], bbox[4], bbox[5],
                    dT, SL_bssn[I_CHI], SL_bssn[I_LAP],
                    SL_em[v], SC_em[v],
                    soa_em[v],
                    Symmetry, cor);
            }
        }

        // lowerboundset on SynchList_cor phi
        f_lowerboundset(ex, SC_bssn[I_CHI], chitiny_val);

        // swapList: swap SynchList_pre <-> SynchList_cor
        {
            double *tmp;
            for (int v = 0; v < NVAR_BSSN; v++) {
                tmp = SP_bssn[v]; SP_bssn[v] = SC_bssn[v]; SC_bssn[v] = tmp;
            }
            for (int v = 0; v < NVAR_EM; v++) {
                tmp = SP_em[v]; SP_em[v] = SC_em[v]; SC_em[v] = tmp;
            }
        }
    }

    // Final swap: StateList <-> SynchList_pre
    {
        double *tmp;
        for (int v = 0; v < NVAR_BSSN; v++) {
            tmp = SL_bssn[v]; SL_bssn[v] = SP_bssn[v]; SP_bssn[v] = tmp;
        }
        for (int v = 0; v < NVAR_EM; v++) {
            tmp = SL_em[v]; SL_em[v] = SP_em[v]; SP_em[v] = tmp;
        }
    }
}

// ===========================================================================
// main — run GPU and CPU paths, compare results
// ===========================================================================
int main()
{
    int ex[3]  = {GW + 8 + GW, GW + 8 + GW, GW + 8 + GW};
    int pex[3] = {ex[0], ex[1], ex[2]};
    int n = ex[0] * ex[1] * ex[2];
    int interior_count = (ex[0] - 2*GW) * (ex[1] - 2*GW) * (ex[2] - 2*GW);

    printf("=== Full Step GPU vs CPU Accuracy Test ===\n");
    printf("Grid: %dx%dx%d, interior: %dx%dx%d, GW=%d, dT=%g, tol=%.0e\n\n",
           ex[0], ex[1], ex[2],
           ex[0] - 2*GW, ex[1] - 2*GW, ex[2] - 2*GW,
           GW, dT, TOL);

    size_t nb = n * sizeof(double);

    // ================================================================
    // Allocate CPU input arrays
    // ================================================================
    double *h_X = new double[n], *h_Y = new double[n], *h_Z = new double[n];

    // StateList (initial data for both paths)
    double *h_SL_bssn[NVAR_BSSN], *h_SL_em[NVAR_EM];
    for (int v = 0; v < NVAR_BSSN; v++) h_SL_bssn[v] = new double[n];
    for (int v = 0; v < NVAR_EM;   v++) h_SL_em[v]   = new double[n];

    // EM source
    double *h_Jx = new double[n], *h_Jy = new double[n];
    double *h_Jz = new double[n], *h_qchar = new double[n];

    // stress-energy
    double *h_rho = new double[n], *h_Sx_cpu = new double[n];
    double *h_Sy = new double[n], *h_Sz = new double[n];
    double *h_Sxx = new double[n], *h_Sxy = new double[n], *h_Sxz = new double[n];
    double *h_Syy = new double[n], *h_Syz = new double[n], *h_Szz = new double[n];

    // Christoffel
    double *h_Gamxxx = new double[n], *h_Gamxxy = new double[n], *h_Gamxxz = new double[n];
    double *h_Gamxyy = new double[n], *h_Gamxyz = new double[n], *h_Gamxzz = new double[n];
    double *h_Gamyxx = new double[n], *h_Gamyxy = new double[n], *h_Gamyxz = new double[n];
    double *h_Gamyyy = new double[n], *h_Gamyyz = new double[n], *h_Gamyzz = new double[n];
    double *h_Gamzxx = new double[n], *h_Gamzxy = new double[n], *h_Gamzxz = new double[n];
    double *h_Gamzyy = new double[n], *h_Gamzyz = new double[n], *h_Gamzzz = new double[n];

    // Ricci
    double *h_Rxx = new double[n], *h_Rxy = new double[n], *h_Rxz = new double[n];
    double *h_Ryy = new double[n], *h_Ryz = new double[n], *h_Rzz = new double[n];

    // Constraints (scratch)
    double *h_ham = new double[n], *h_movx = new double[n], *h_movy = new double[n];
    double *h_movz = new double[n], *h_Gmx_c = new double[n], *h_Gmy = new double[n], *h_Gmz = new double[n];

    // Fill initial data
    fill_test_data(ex, h_X, h_Y, h_Z,
        h_SL_bssn[I_CHI], h_SL_bssn[I_TRK],
        h_SL_bssn[I_GXX], h_SL_bssn[I_GXY], h_SL_bssn[I_GXZ],
        h_SL_bssn[I_GYY], h_SL_bssn[I_GYZ], h_SL_bssn[I_GZZ],
        h_SL_bssn[I_AXX], h_SL_bssn[I_AXY], h_SL_bssn[I_AXZ],
        h_SL_bssn[I_AYY], h_SL_bssn[I_AYZ], h_SL_bssn[I_AZZ],
        h_SL_bssn[I_GAMX], h_SL_bssn[I_GAMY], h_SL_bssn[I_GAMZ],
        h_SL_bssn[I_LAP],
        h_SL_bssn[I_BETAX], h_SL_bssn[I_BETAY], h_SL_bssn[I_BETAZ],
        h_SL_bssn[I_DTSFX], h_SL_bssn[I_DTSFY], h_SL_bssn[I_DTSFZ],
        h_SL_em[I_Ex - NVAR_BSSN], h_SL_em[I_Ey - NVAR_BSSN], h_SL_em[I_Ez - NVAR_BSSN],
        h_SL_em[I_Bx - NVAR_BSSN], h_SL_em[I_By - NVAR_BSSN], h_SL_em[I_Bz - NVAR_BSSN],
        h_SL_em[I_KPSI - NVAR_BSSN], h_SL_em[I_KPHI - NVAR_BSSN],
        h_Jx, h_Jy, h_Jz, h_qchar,
        h_rho, h_Sx_cpu, h_Sy, h_Sz,
        h_Sxx, h_Sxy, h_Sxz,
        h_Syy, h_Syz, h_Szz,
        h_Gamxxx, h_Gamxxy, h_Gamxxz, h_Gamxyy, h_Gamxyz, h_Gamxzz,
        h_Gamyxx, h_Gamyxy, h_Gamyxz, h_Gamyyy, h_Gamyyz, h_Gamyzz,
        h_Gamzxx, h_Gamzxy, h_Gamzxz, h_Gamzyy, h_Gamzyz, h_Gamzzz,
        h_Rxx, h_Rxy, h_Rxz, h_Ryy, h_Ryz, h_Rzz);

    // ================================================================
    // GPU PATH -- allocate GPU arrays and run
    // ================================================================
    printf("=== GPU Path ===\n");

    cudaStream_t stream = 0;

    // GPU copies of all input arrays
    double *d_X = gpu_alloc(n), *d_Y = gpu_alloc(n), *d_Z = gpu_alloc(n);

    double *d_SL_bssn[NVAR_BSSN], *d_SL_em[NVAR_EM];
    double *d_SP_bssn[NVAR_BSSN], *d_SP_em[NVAR_EM];
    double *d_SC_bssn[NVAR_BSSN], *d_SC_em[NVAR_EM];
    double *d_RL_bssn[NVAR_BSSN], *d_RL_em[NVAR_EM];
    for (int v = 0; v < NVAR_BSSN; v++) {
        d_SL_bssn[v] = gpu_alloc(n);
        d_SP_bssn[v] = gpu_alloc(n);
        d_SC_bssn[v] = gpu_alloc(n);
        d_RL_bssn[v] = gpu_alloc(n);
    }
    for (int v = 0; v < NVAR_EM; v++) {
        d_SL_em[v] = gpu_alloc(n);
        d_SP_em[v] = gpu_alloc(n);
        d_SC_em[v] = gpu_alloc(n);
        d_RL_em[v] = gpu_alloc(n);
    }

    double *d_Jx = gpu_alloc(n), *d_Jy = gpu_alloc(n), *d_Jz = gpu_alloc(n), *d_qchar = gpu_alloc(n);
    double *d_rho = gpu_alloc(n), *d_Sx = gpu_alloc(n), *d_Sy = gpu_alloc(n), *d_Sz = gpu_alloc(n);
    double *d_Sxx = gpu_alloc(n), *d_Sxy = gpu_alloc(n), *d_Sxz = gpu_alloc(n);
    double *d_Syy = gpu_alloc(n), *d_Syz = gpu_alloc(n), *d_Szz = gpu_alloc(n);

    double *d_Gamxxx = gpu_alloc(n), *d_Gamxxy = gpu_alloc(n), *d_Gamxxz = gpu_alloc(n);
    double *d_Gamxyy = gpu_alloc(n), *d_Gamxyz = gpu_alloc(n), *d_Gamxzz = gpu_alloc(n);
    double *d_Gamyxx = gpu_alloc(n), *d_Gamyxy = gpu_alloc(n), *d_Gamyxz = gpu_alloc(n);
    double *d_Gamyyy = gpu_alloc(n), *d_Gamyyz = gpu_alloc(n), *d_Gamyzz = gpu_alloc(n);
    double *d_Gamzxx = gpu_alloc(n), *d_Gamzxy = gpu_alloc(n), *d_Gamzxz = gpu_alloc(n);
    double *d_Gamzyy = gpu_alloc(n), *d_Gamzyz = gpu_alloc(n), *d_Gamzzz = gpu_alloc(n);

    double *d_Rxx = gpu_alloc(n), *d_Rxy = gpu_alloc(n), *d_Rxz = gpu_alloc(n);
    double *d_Ryy = gpu_alloc(n), *d_Ryz = gpu_alloc(n), *d_Rzz = gpu_alloc(n);

    double *d_ham = gpu_alloc(n), *d_movx = gpu_alloc(n), *d_movy = gpu_alloc(n);
    double *d_movz = gpu_alloc(n), *d_Gmx = gpu_alloc(n), *d_Gmy = gpu_alloc(n), *d_Gmz = gpu_alloc(n);

    // Copy inputs to GPU
    gpu_copy_h2d(d_X, h_X, n); gpu_copy_h2d(d_Y, h_Y, n); gpu_copy_h2d(d_Z, h_Z, n);
    for (int v = 0; v < NVAR_BSSN; v++) {
        gpu_copy_h2d(d_SL_bssn[v], h_SL_bssn[v], n);
        CUDA_CHECK(cudaMemset(d_SP_bssn[v], 0, nb));
        CUDA_CHECK(cudaMemset(d_SC_bssn[v], 0, nb));
        CUDA_CHECK(cudaMemset(d_RL_bssn[v], 0, nb));
    }
    for (int v = 0; v < NVAR_EM; v++) {
        gpu_copy_h2d(d_SL_em[v], h_SL_em[v], n);
        CUDA_CHECK(cudaMemset(d_SP_em[v], 0, nb));
        CUDA_CHECK(cudaMemset(d_SC_em[v], 0, nb));
        CUDA_CHECK(cudaMemset(d_RL_em[v], 0, nb));
    }
    gpu_copy_h2d(d_Jx, h_Jx, n); gpu_copy_h2d(d_Jy, h_Jy, n);
    gpu_copy_h2d(d_Jz, h_Jz, n); gpu_copy_h2d(d_qchar, h_qchar, n);
    gpu_copy_h2d(d_rho, h_rho, n); gpu_copy_h2d(d_Sx, h_Sx_cpu, n);
    gpu_copy_h2d(d_Sy, h_Sy, n); gpu_copy_h2d(d_Sz, h_Sz, n);
    gpu_copy_h2d(d_Sxx, h_Sxx, n); gpu_copy_h2d(d_Sxy, h_Sxy, n);
    gpu_copy_h2d(d_Sxz, h_Sxz, n);
    gpu_copy_h2d(d_Syy, h_Syy, n); gpu_copy_h2d(d_Syz, h_Syz, n);
    gpu_copy_h2d(d_Szz, h_Szz, n);
    gpu_copy_h2d(d_Gamxxx, h_Gamxxx, n); gpu_copy_h2d(d_Gamxxy, h_Gamxxy, n);
    gpu_copy_h2d(d_Gamxxz, h_Gamxxz, n); gpu_copy_h2d(d_Gamxyy, h_Gamxyy, n);
    gpu_copy_h2d(d_Gamxyz, h_Gamxyz, n); gpu_copy_h2d(d_Gamxzz, h_Gamxzz, n);
    gpu_copy_h2d(d_Gamyxx, h_Gamyxx, n); gpu_copy_h2d(d_Gamyxy, h_Gamyxy, n);
    gpu_copy_h2d(d_Gamyxz, h_Gamyxz, n); gpu_copy_h2d(d_Gamyyy, h_Gamyyy, n);
    gpu_copy_h2d(d_Gamyyz, h_Gamyyz, n); gpu_copy_h2d(d_Gamyzz, h_Gamyzz, n);
    gpu_copy_h2d(d_Gamzxx, h_Gamzxx, n); gpu_copy_h2d(d_Gamzxy, h_Gamzxy, n);
    gpu_copy_h2d(d_Gamzxz, h_Gamzxz, n); gpu_copy_h2d(d_Gamzyy, h_Gamzyy, n);
    gpu_copy_h2d(d_Gamzyz, h_Gamzyz, n); gpu_copy_h2d(d_Gamzzz, h_Gamzzz, n);
    gpu_copy_h2d(d_Rxx, h_Rxx, n); gpu_copy_h2d(d_Rxy, h_Rxy, n);
    gpu_copy_h2d(d_Rxz, h_Rxz, n); gpu_copy_h2d(d_Ryy, h_Ryy, n);
    gpu_copy_h2d(d_Ryz, h_Ryz, n); gpu_copy_h2d(d_Rzz, h_Rzz, n);

    // Run GPU step
    run_gpu_step(ex, pex, n, stream,
        d_X, d_Y, d_Z,
        d_SL_bssn, d_SL_em,
        d_SP_bssn, d_SP_em,
        d_SC_bssn, d_SC_em,
        d_RL_bssn, d_RL_em,
        d_Jx, d_Jy, d_Jz, d_qchar,
        d_rho, d_Sx, d_Sy, d_Sz,
        d_Sxx, d_Sxy, d_Sxz,
        d_Syy, d_Syz, d_Szz,
        d_Gamxxx, d_Gamxxy, d_Gamxxz, d_Gamxyy, d_Gamxyz, d_Gamxzz,
        d_Gamyxx, d_Gamyxy, d_Gamyxz, d_Gamyyy, d_Gamyyz, d_Gamyzz,
        d_Gamzxx, d_Gamzxy, d_Gamzxz, d_Gamzyy, d_Gamzyz, d_Gamzzz,
        d_Rxx, d_Rxy, d_Rxz, d_Ryy, d_Ryz, d_Rzz,
        d_ham, d_movx, d_movy, d_movz, d_Gmx, d_Gmy, d_Gmz);

    // Copy GPU results back for comparison
    printf("  Copying GPU results back...\n");
    double *g_SL_bssn[NVAR_BSSN], *g_SL_em[NVAR_EM];
    for (int v = 0; v < NVAR_BSSN; v++) {
        g_SL_bssn[v] = new double[n];
        gpu_copy_d2h(g_SL_bssn[v], d_SL_bssn[v], n);
    }
    for (int v = 0; v < NVAR_EM; v++) {
        g_SL_em[v] = new double[n];
        gpu_copy_d2h(g_SL_em[v], d_SL_em[v], n);
    }

    printf("  GPU step done.\n\n");

    // ================================================================
    // CPU PATH -- use second copy of initial data, run Fortran step
    // ================================================================
    printf("=== CPU Path ===\n");

    // Allocate CPU working arrays
    double *c_SL_bssn[NVAR_BSSN], *c_SL_em[NVAR_EM];
    double *c_SP_bssn[NVAR_BSSN], *c_SP_em[NVAR_EM];
    double *c_SC_bssn[NVAR_BSSN], *c_SC_em[NVAR_EM];
    double *c_RL_bssn[NVAR_BSSN], *c_RL_em[NVAR_EM];
    for (int v = 0; v < NVAR_BSSN; v++) {
        c_SL_bssn[v] = new double[n]; memcpy(c_SL_bssn[v], h_SL_bssn[v], nb);
        c_SP_bssn[v] = new double[n]; memset(c_SP_bssn[v], 0, nb);
        c_SC_bssn[v] = new double[n]; memset(c_SC_bssn[v], 0, nb);
        c_RL_bssn[v] = new double[n]; memset(c_RL_bssn[v], 0, nb);
    }
    for (int v = 0; v < NVAR_EM; v++) {
        c_SL_em[v] = new double[n]; memcpy(c_SL_em[v], h_SL_em[v], nb);
        c_SP_em[v] = new double[n]; memset(c_SP_em[v], 0, nb);
        c_SC_em[v] = new double[n]; memset(c_SC_em[v], 0, nb);
        c_RL_em[v] = new double[n]; memset(c_RL_em[v], 0, nb);
    }

    // CPU working copies of fixed input arrays
    double *c_X = new double[n]; memcpy(c_X, h_X, nb);
    double *c_Y = new double[n]; memcpy(c_Y, h_Y, nb);
    double *c_Z = new double[n]; memcpy(c_Z, h_Z, nb);
    double *c_Jx = new double[n]; memcpy(c_Jx, h_Jx, nb);
    double *c_Jy = new double[n]; memcpy(c_Jy, h_Jy, nb);
    double *c_Jz = new double[n]; memcpy(c_Jz, h_Jz, nb);
    double *c_qchar = new double[n]; memcpy(c_qchar, h_qchar, nb);
    double *c_rho = new double[n]; memcpy(c_rho, h_rho, nb);
    double *c_Sx = new double[n]; memcpy(c_Sx, h_Sx_cpu, nb);
    double *c_Sy = new double[n]; memcpy(c_Sy, h_Sy, nb);
    double *c_Sz = new double[n]; memcpy(c_Sz, h_Sz, nb);
    double *c_Sxx = new double[n]; memcpy(c_Sxx, h_Sxx, nb);
    double *c_Sxy = new double[n]; memcpy(c_Sxy, h_Sxy, nb);
    double *c_Sxz = new double[n]; memcpy(c_Sxz, h_Sxz, nb);
    double *c_Syy = new double[n]; memcpy(c_Syy, h_Syy, nb);
    double *c_Syz = new double[n]; memcpy(c_Syz, h_Syz, nb);
    double *c_Szz = new double[n]; memcpy(c_Szz, h_Szz, nb);
    double *c_Gamxxx = new double[n]; memcpy(c_Gamxxx, h_Gamxxx, nb);
    double *c_Gamxxy = new double[n]; memcpy(c_Gamxxy, h_Gamxxy, nb);
    double *c_Gamxxz = new double[n]; memcpy(c_Gamxxz, h_Gamxxz, nb);
    double *c_Gamxyy = new double[n]; memcpy(c_Gamxyy, h_Gamxyy, nb);
    double *c_Gamxyz = new double[n]; memcpy(c_Gamxyz, h_Gamxyz, nb);
    double *c_Gamxzz = new double[n]; memcpy(c_Gamxzz, h_Gamxzz, nb);
    double *c_Gamyxx = new double[n]; memcpy(c_Gamyxx, h_Gamyxx, nb);
    double *c_Gamyxy = new double[n]; memcpy(c_Gamyxy, h_Gamyxy, nb);
    double *c_Gamyxz = new double[n]; memcpy(c_Gamyxz, h_Gamyxz, nb);
    double *c_Gamyyy = new double[n]; memcpy(c_Gamyyy, h_Gamyyy, nb);
    double *c_Gamyyz = new double[n]; memcpy(c_Gamyyz, h_Gamyyz, nb);
    double *c_Gamyzz = new double[n]; memcpy(c_Gamyzz, h_Gamyzz, nb);
    double *c_Gamzxx = new double[n]; memcpy(c_Gamzxx, h_Gamzxx, nb);
    double *c_Gamzxy = new double[n]; memcpy(c_Gamzxy, h_Gamzxy, nb);
    double *c_Gamzxz = new double[n]; memcpy(c_Gamzxz, h_Gamzxz, nb);
    double *c_Gamzyy = new double[n]; memcpy(c_Gamzyy, h_Gamzyy, nb);
    double *c_Gamzyz = new double[n]; memcpy(c_Gamzyz, h_Gamzyz, nb);
    double *c_Gamzzz = new double[n]; memcpy(c_Gamzzz, h_Gamzzz, nb);
    double *c_Rxx = new double[n]; memcpy(c_Rxx, h_Rxx, nb);
    double *c_Rxy = new double[n]; memcpy(c_Rxy, h_Rxy, nb);
    double *c_Rxz = new double[n]; memcpy(c_Rxz, h_Rxz, nb);
    double *c_Ryy = new double[n]; memcpy(c_Ryy, h_Ryy, nb);
    double *c_Ryz = new double[n]; memcpy(c_Ryz, h_Ryz, nb);
    double *c_Rzz = new double[n]; memcpy(c_Rzz, h_Rzz, nb);
    double *c_ham = new double[n]; double *c_movx = new double[n];
    double *c_movy = new double[n]; double *c_movz = new double[n];
    double *c_Gmx = new double[n]; double *c_Gmy_c = new double[n]; double *c_Gmz = new double[n];

    // Run CPU step
    run_cpu_step(ex, pex,
        c_X, c_Y, c_Z,
        c_SL_bssn, c_SL_em,
        c_SP_bssn, c_SP_em,
        c_SC_bssn, c_SC_em,
        c_RL_bssn, c_RL_em,
        c_Jx, c_Jy, c_Jz, c_qchar,
        c_rho, c_Sx, c_Sy, c_Sz,
        c_Sxx, c_Sxy, c_Sxz,
        c_Syy, c_Syz, c_Szz,
        c_Gamxxx, c_Gamxxy, c_Gamxxz, c_Gamxyy, c_Gamxyz, c_Gamxzz,
        c_Gamyxx, c_Gamyxy, c_Gamyxz, c_Gamyyy, c_Gamyyz, c_Gamyzz,
        c_Gamzxx, c_Gamzxy, c_Gamzxz, c_Gamzyy, c_Gamzyz, c_Gamzzz,
        c_Rxx, c_Rxy, c_Rxz, c_Ryy, c_Ryz, c_Rzz,
        c_ham, c_movx, c_movy, c_movz, c_Gmx, c_Gmy_c, c_Gmz);

    printf("  CPU step done.\n\n");

    // ================================================================
    // Compare results (StateList only — final evolved fields)
    // ================================================================
    printf("=== Comparison (interior [%d:%d) only) ===\n", GW, ex[0] - GW);

    const char *bssn_names[NVAR_BSSN] = {
        "chi",  "trK",
        "gxx",  "gxy",  "gxz",  "gyy",  "gyz",  "gzz",
        "Axx",  "Axy",  "Axz",  "Ayy",  "Ayz",  "Azz",
        "Gamx", "Gamy", "Gamz",
        "Lap",
        "betax","betay","betaz",
        "dtSfx","dtSfy","dtSfz"
    };
    const char *em_names[NVAR_EM] = {
        "Ex","Ey","Ez","Bx","By","Bz","Kpsi","Kphi"
    };

    int total_fails = 0;
    for (int v = 0; v < NVAR_BSSN; v++) {
        total_fails += compare_field(bssn_names[v], g_SL_bssn[v], c_SL_bssn[v], ex, TOL);
    }
    for (int v = 0; v < NVAR_EM; v++) {
        total_fails += compare_field(em_names[v], g_SL_em[v], c_SL_em[v], ex, TOL);
    }

    printf("\n");
    if (total_fails == 0) {
        printf("RESULT: ALL PASS (tolerance=%.0e, interior=%d points)\n", TOL, interior_count);
    } else {
        printf("RESULT: %d FAILURES (tolerance=%.0e, interior=%d points)\n",
               total_fails, TOL, interior_count);
    }

    // ================================================================
    // Cleanup
    // ================================================================
    delete[] h_X; delete[] h_Y; delete[] h_Z;
    for (int v = 0; v < NVAR_BSSN; v++) delete[] h_SL_bssn[v];
    for (int v = 0; v < NVAR_EM;   v++) delete[] h_SL_em[v];
    delete[] h_Jx; delete[] h_Jy; delete[] h_Jz; delete[] h_qchar;
    delete[] h_rho; delete[] h_Sx_cpu; delete[] h_Sy; delete[] h_Sz;
    delete[] h_Sxx; delete[] h_Sxy; delete[] h_Sxz;
    delete[] h_Syy; delete[] h_Syz; delete[] h_Szz;
    delete[] h_Gamxxx; delete[] h_Gamxxy; delete[] h_Gamxxz;
    delete[] h_Gamxyy; delete[] h_Gamxyz; delete[] h_Gamxzz;
    delete[] h_Gamyxx; delete[] h_Gamyxy; delete[] h_Gamyxz;
    delete[] h_Gamyyy; delete[] h_Gamyyz; delete[] h_Gamyzz;
    delete[] h_Gamzxx; delete[] h_Gamzxy; delete[] h_Gamzxz;
    delete[] h_Gamzyy; delete[] h_Gamzyz; delete[] h_Gamzzz;
    delete[] h_Rxx; delete[] h_Rxy; delete[] h_Rxz;
    delete[] h_Ryy; delete[] h_Ryz; delete[] h_Rzz;
    delete[] h_ham; delete[] h_movx; delete[] h_movy;
    delete[] h_movz; delete[] h_Gmx_c; delete[] h_Gmy; delete[] h_Gmz;

    for (int v = 0; v < NVAR_BSSN; v++) delete[] g_SL_bssn[v];
    for (int v = 0; v < NVAR_EM;   v++) delete[] g_SL_em[v];

    for (int v = 0; v < NVAR_BSSN; v++) {
        delete[] c_SL_bssn[v]; delete[] c_SP_bssn[v];
        delete[] c_SC_bssn[v]; delete[] c_RL_bssn[v];
    }
    for (int v = 0; v < NVAR_EM; v++) {
        delete[] c_SL_em[v]; delete[] c_SP_em[v];
        delete[] c_SC_em[v]; delete[] c_RL_em[v];
    }
    delete[] c_X; delete[] c_Y; delete[] c_Z;
    delete[] c_Jx; delete[] c_Jy; delete[] c_Jz; delete[] c_qchar;
    delete[] c_rho; delete[] c_Sx; delete[] c_Sy; delete[] c_Sz;
    delete[] c_Sxx; delete[] c_Sxy; delete[] c_Sxz;
    delete[] c_Syy; delete[] c_Syz; delete[] c_Szz;
    delete[] c_Gamxxx; delete[] c_Gamxxy; delete[] c_Gamxxz;
    delete[] c_Gamxyy; delete[] c_Gamxyz; delete[] c_Gamxzz;
    delete[] c_Gamyxx; delete[] c_Gamyxy; delete[] c_Gamyxz;
    delete[] c_Gamyyy; delete[] c_Gamyyz; delete[] c_Gamyzz;
    delete[] c_Gamzxx; delete[] c_Gamzxy; delete[] c_Gamzxz;
    delete[] c_Gamzyy; delete[] c_Gamzyz; delete[] c_Gamzzz;
    delete[] c_Rxx; delete[] c_Rxy; delete[] c_Rxz;
    delete[] c_Ryy; delete[] c_Ryz; delete[] c_Rzz;
    delete[] c_ham; delete[] c_movx; delete[] c_movy;
    delete[] c_movz; delete[] c_Gmx; delete[] c_Gmy_c; delete[] c_Gmz;

    // GPU free
    CUDA_CHECK(cudaFree(d_X)); CUDA_CHECK(cudaFree(d_Y)); CUDA_CHECK(cudaFree(d_Z));
    for (int v = 0; v < NVAR_BSSN; v++) {
        CUDA_CHECK(cudaFree(d_SL_bssn[v])); CUDA_CHECK(cudaFree(d_SP_bssn[v]));
        CUDA_CHECK(cudaFree(d_SC_bssn[v])); CUDA_CHECK(cudaFree(d_RL_bssn[v]));
    }
    for (int v = 0; v < NVAR_EM; v++) {
        CUDA_CHECK(cudaFree(d_SL_em[v])); CUDA_CHECK(cudaFree(d_SP_em[v]));
        CUDA_CHECK(cudaFree(d_SC_em[v])); CUDA_CHECK(cudaFree(d_RL_em[v]));
    }
    CUDA_CHECK(cudaFree(d_Jx)); CUDA_CHECK(cudaFree(d_Jy));
    CUDA_CHECK(cudaFree(d_Jz)); CUDA_CHECK(cudaFree(d_qchar));
    CUDA_CHECK(cudaFree(d_rho)); CUDA_CHECK(cudaFree(d_Sx));
    CUDA_CHECK(cudaFree(d_Sy));  CUDA_CHECK(cudaFree(d_Sz));
    CUDA_CHECK(cudaFree(d_Sxx)); CUDA_CHECK(cudaFree(d_Sxy));
    CUDA_CHECK(cudaFree(d_Sxz)); CUDA_CHECK(cudaFree(d_Syy));
    CUDA_CHECK(cudaFree(d_Syz)); CUDA_CHECK(cudaFree(d_Szz));
    CUDA_CHECK(cudaFree(d_Gamxxx)); CUDA_CHECK(cudaFree(d_Gamxxy));
    CUDA_CHECK(cudaFree(d_Gamxxz)); CUDA_CHECK(cudaFree(d_Gamxyy));
    CUDA_CHECK(cudaFree(d_Gamxyz)); CUDA_CHECK(cudaFree(d_Gamxzz));
    CUDA_CHECK(cudaFree(d_Gamyxx)); CUDA_CHECK(cudaFree(d_Gamyxy));
    CUDA_CHECK(cudaFree(d_Gamyxz)); CUDA_CHECK(cudaFree(d_Gamyyy));
    CUDA_CHECK(cudaFree(d_Gamyyz)); CUDA_CHECK(cudaFree(d_Gamyzz));
    CUDA_CHECK(cudaFree(d_Gamzxx)); CUDA_CHECK(cudaFree(d_Gamzxy));
    CUDA_CHECK(cudaFree(d_Gamzxz)); CUDA_CHECK(cudaFree(d_Gamzyy));
    CUDA_CHECK(cudaFree(d_Gamzyz)); CUDA_CHECK(cudaFree(d_Gamzzz));
    CUDA_CHECK(cudaFree(d_Rxx)); CUDA_CHECK(cudaFree(d_Rxy));
    CUDA_CHECK(cudaFree(d_Rxz)); CUDA_CHECK(cudaFree(d_Ryy));
    CUDA_CHECK(cudaFree(d_Ryz)); CUDA_CHECK(cudaFree(d_Rzz));
    CUDA_CHECK(cudaFree(d_ham));  CUDA_CHECK(cudaFree(d_movx));
    CUDA_CHECK(cudaFree(d_movy)); CUDA_CHECK(cudaFree(d_movz));
    CUDA_CHECK(cudaFree(d_Gmx));  CUDA_CHECK(cudaFree(d_Gmy));
    CUDA_CHECK(cudaFree(d_Gmz));

    return (total_fails == 0) ? 0 : 1;
}
