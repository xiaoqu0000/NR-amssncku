// EXTERNAL AHF-ONLY candidate. Fortran REAL stencil literal, original symmetry/fetches.
// No CPU coordinate injection. Caller passes device coordinates, runtime 2*ghost_width.
#ifndef AHF_NEVILLE_EXACTXYZ_CUH
#define AHF_NEVILLE_EXACTXYZ_CUH
#include "macrodef.h"
#include "fmisc_gpu.cuh"
__device__ __forceinline__ void ahf_neville_interp_device(
    const int* ex, const int* pex,
    const double* X, const double* Y, const double* Z,
    const double* f, double* f_int,
    double x1, double y1, double z1,
    int ORDN, const double* SoA, int symmetry
) {
    if (ORDN <= 0 || ORDN > MAX_ORDN) { f_int[0] = NAN; gpu_stop(); return; }

    const int NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2;
    (void)EQUATORIAL;

    int imin = 1, jmin = 1, kmin = 1;

    double dX = X_at_1b(X, imin + 1) - X_at_1b(X, imin);
    double dY = X_at_1b(Y, jmin + 1) - X_at_1b(Y, jmin);
    double dZ = X_at_1b(Z, kmin + 1) - X_at_1b(Z, kmin);

    // Match Fortran's default-REAL literal before promotion to double.
    constexpr double fortran_shift = static_cast<double>(0.4f);
    int cxI[3];
    cxI[0] = (int)((x1 - X_at_1b(X, 1)) / dX + fortran_shift) + 1;
    cxI[1] = (int)((y1 - X_at_1b(Y, 1)) / dY + fortran_shift) + 1;
    cxI[2] = (int)((z1 - X_at_1b(Z, 1)) / dZ + fortran_shift) + 1;

    int cxB[3], cxT[3], cmin[3], cmax[3];
    for (int m = 0; m < 3; ++m) {
        cxB[m] = cxI[m] - ORDN / 2 + 1;
        cxT[m] = cxB[m] + ORDN - 1;
        cmin[m] = 1;
        cmax[m] = ex[m];
    }

#ifdef Vertex
    int sym_shift = -ORDN / 2 + 2;
#else
    int sym_shift = -ORDN / 2 + 1;
#endif

    if (symmetry == OCTANT && fabs(X_at_1b(X, 1)) < dX) cmin[0] = sym_shift;
    if (symmetry == OCTANT && fabs(X_at_1b(Y, 1)) < dY) cmin[1] = sym_shift;
    if (symmetry != NO_SYMM && fabs(X_at_1b(Z, 1)) < dZ) cmin[2] = sym_shift;

    for (int m = 0; m < 3; ++m) {
        if (cxB[m] < cmin[m]) { cxB[m] = cmin[m]; cxT[m] = cxB[m] + ORDN - 1; }
        if (cxT[m] > cmax[m]) { cxT[m] = cmax[m]; cxB[m] = cxT[m] + 1 - ORDN; }
    }

#ifdef Vertex
    int ref_shift = 2;
#else
    int ref_shift = 1;
#endif

    double cx[3];
    cx[0] = (cxB[0] > 0) ? (x1 - X_at_1b(X, cxB[0])) / dX : (x1 + X_at_1b(X, ref_shift - cxB[0])) / dX;
    cx[1] = (cxB[1] > 0) ? (y1 - X_at_1b(Y, cxB[1])) / dY : (y1 + X_at_1b(Y, ref_shift - cxB[1])) / dY;
    cx[2] = (cxB[2] > 0) ? (z1 - X_at_1b(Z, cxB[2])) / dZ : (z1 + X_at_1b(Z, ref_shift - cxB[2])) / dZ;

    // polin3_device consumes ((i*ORDN+j)*ORDN+k), NOT Fortran x-fastest.
    double nodes[MAX_ORDN], ya[MAX_ORDN * MAX_ORDN * MAX_ORDN], dy;
    for (int i=0;i<ORDN;++i) nodes[i]=static_cast<double>(i);
    for (int i=0;i<ORDN;++i) for (int j=0;j<ORDN;++j) for (int k=0;k<ORDN;++k) {
        double value = fetch_f_sym_pad(f, ex, pex, cxB[0]+i, cxB[1]+j, cxB[2]+k, SoA);
        if (isnan(value)) { f_int[0]=NAN; gpu_stop(); return; }
        ya[(i*ORDN+j)*ORDN+k]=value;
    }
    f_int[0]=polin3_device(nodes,nodes,nodes,ya,cx[0],cx[1],cx[2],&dy,ORDN);
}
__device__ __forceinline__ void ahf_neville_interp_ss_device(
    const int* ex, const int* pex,
    const double* X, const double* Y, const double* Z,
    const double* f, double* f_int,
    double x1, double y1, double z1,
    int ORDN, const double* SoA, int symmetry, int sst
) {
    if (ORDN <= 0 || ORDN > MAX_ORDN) { f_int[0] = NAN; gpu_stop(); return; }

    double dX = X_at_1b(X, 2) - X_at_1b(X, 1);
    double dY = X_at_1b(Y, 2) - X_at_1b(Y, 1);
    double dZ = X_at_1b(Z, 2) - X_at_1b(Z, 1);

    // Fortran's default-REAL 0.4 is rounded to float before double promotion.
    constexpr double fortran_shift = static_cast<double>(0.4f);
    int cxI[3];
    cxI[0] = (int)((x1 - X_at_1b(X, 1)) / dX + fortran_shift) + 1;
    cxI[1] = (int)((y1 - X_at_1b(Y, 1)) / dY + fortran_shift) + 1;
    cxI[2] = (int)((z1 - X_at_1b(Z, 1)) / dZ + fortran_shift) + 1;

    int cxB[3], cxT[3], cmin[3], cmax[3];
    for (int m = 0; m < 3; ++m) {
        cxB[m] = cxI[m] - ORDN / 2 + 1;
        cxT[m] = cxB[m] + ORDN - 1;
        cmin[m] = 1;
        cmax[m] = ex[m];
    }

    // [修复点 1]：必须根据 sst 生成 SoAh
    double SoAh[3];
    shell_soa3(sst, SoA, SoAh);

    if (sst == -1) {
        cmin[0] = -ORDN + 1; cmin[1] = -ORDN + 1; cmin[2] = -ORDN + 1;
        cmax[0] = ex[0]; cmax[1] = ex[1]; cmax[2] = ex[2];
    } else if (sst == 0 || sst == 1 || sst == 2 || sst == 3 || sst == 4 || sst == 5) {
        cmin[0] = -ORDN + 1; cmin[1] = -ORDN + 1;
        cmax[0] = ex[0] + ORDN; cmax[1] = ex[1] + ORDN;
    }

    for (int m = 0; m < 3; ++m) {
        if (cxB[m] < cmin[m]) { cxB[m] = cmin[m]; cxT[m] = cxB[m] + ORDN - 1; }
        if (cxT[m] > cmax[m]) { cxT[m] = cmax[m]; cxB[m] = cxT[m] + 1 - ORDN; }
    }

    // [修复点 2]：忠于 Fortran，直接代数平移计算 cx，不能带 if 条件
    double cx[3];
    cx[0] = (x1 - X_at_1b(X, 1)) / dX - cxB[0] + 1.0;
    cx[1] = (y1 - X_at_1b(Y, 1)) / dY - cxB[1] + 1.0;
    cx[2] = (z1 - X_at_1b(Z, 1)) / dZ - cxB[2] + 1.0;

    // polin3_device consumes ((i*ORDN+j)*ORDN+k), NOT Fortran x-fastest.
    double nodes[MAX_ORDN], ya[MAX_ORDN * MAX_ORDN * MAX_ORDN], dy;
    for (int i=0;i<ORDN;++i) nodes[i]=static_cast<double>(i);
    for (int i=0;i<ORDN;++i) for (int j=0;j<ORDN;++j) for (int k=0;k<ORDN;++k) {
        double value = fetch_f_tbd_0b_pad(f, ex, pex, cxB[0]+i-1, cxB[1]+j-1, cxB[2]+k-1, SoAh);
        if (isnan(value)) { f_int[0]=NAN; gpu_stop(); return; }
        ya[(i*ORDN+j)*ORDN+k]=value;
    }
    f_int[0]=polin3_device(nodes,nodes,nodes,ya,cx[0],cx[1],cx[2],&dy,ORDN);
}
#endif
