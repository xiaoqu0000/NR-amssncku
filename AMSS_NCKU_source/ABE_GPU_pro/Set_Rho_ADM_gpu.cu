#include "macrodef.h"
// =============================================================================
// Set_Rho_ADM_gpu.cu
//
// CUDA GPU kernels for ADM density and scalar field potential functions.
// Port of Set_Rho_ADM.f90 for f(R) scalar-tensor theory.
//
// 5 scalar field configurations (CC=1..5):
//   CC=1: Case C of 1112.3928, V=0
//         phi = phi0 * tanh((r-r0)/sigma)
//   CC=2: shell with a2^2*phi0/(1+a2^2), f(R)=R+a2*R^2 induced V
//         phi = phi0_scaled * exp(-(r-r0)^2/sigma)
//   CC=3: ground state of Schrodinger-Newton system, f(R)=R+a2*R^2 induced V
//         phi = 5-Gaussian sum / l2
//   CC=4: a2 = oo, phi = phi0*0.5*(tanh((r+r0)/sigma)-tanh((r-r0)/sigma)), V=0
//   CC=5: shell with phi0*exp(-(r-r0)^2/sigma), V=0
//
// Model: f(R) = R + a2 * R^2
// =============================================================================

#include "gpu_common.cuh"
#include "macrodef.fh"

#include <cuda_runtime.h>
#include <cmath>

// =============================================================================
// Pointwise device helper functions
// =============================================================================

// ---------------------------------------------------------------------------
// d_phi<CC> -- scalar field value phi(X,Y,Z) at a point
// ---------------------------------------------------------------------------
template<int CC>
__device__ __forceinline__
double d_phi(double X, double Y, double Z,
             double a2, double r0, double phi0, double sigma, double l2)
{
    const double r = sqrt(X * X + Y * Y + Z * Z);

    if constexpr (CC == 1) {
        // Configuration 1: phi0 * tanh((r-r0)/sigma)
        return phi0 * tanh((r - r0) / sigma);
    }
    else if constexpr (CC == 2) {
        // Configuration 2: scaled phi0 * exp(-(r-r0)^2 / sigma)
        const double phi0_s = a2 * a2 * phi0 / (1.0 + a2 * a2);
        const double dr = r - r0;
        return phi0_s * exp(-dr * dr / sigma);
    }
    else if constexpr (CC == 3) {
        // Configuration 3: ground state of Schrodinger-Newton
        // Sum of 5 Gaussians divided by l2
        return (0.0481646  * exp(-0.0581545  * (r - 1.8039e-8)  * (r - 1.8039e-8)  / l2) +
                0.298408   * exp(-0.111412   * (r + 9.6741e-9)  * (r + 9.6741e-9)  / l2) +
                0.42755    * exp(-0.207156   * (r - 1.09822e-8) * (r - 1.09822e-8) / l2) +
                0.204229   * exp(-0.37742    * (r + 2.13778e-8) * (r + 2.13778e-8) / l2) +
                0.021649   * exp(-0.68406    * (r - 8.78608e-8) * (r - 8.78608e-8) / l2)) / l2;
    }
    else if constexpr (CC == 4) {
        // Configuration 4: a2 = oo, two tanh walls
        const double phi0_h = 0.5 * phi0;
        return phi0_h * (tanh((r + r0) / sigma) - tanh((r - r0) / sigma));
    }
    else if constexpr (CC == 5) {
        // Configuration 5: phi0 * exp(-(r-r0)^2 / sigma)
        const double dr = r - r0;
        return phi0 * exp(-dr * dr / sigma);
    }
    return 0.0;
}

// ---------------------------------------------------------------------------
// d_dphi<CC> -- radial derivative d(phi)/dr at a point
// ---------------------------------------------------------------------------
template<int CC>
__device__ __forceinline__
double d_dphi(double X, double Y, double Z,
              double a2, double r0, double phi0, double sigma, double l2)
{
    const double r = sqrt(X * X + Y * Y + Z * Z);

    if constexpr (CC == 1) {
        // d/dr: phi0 * tanh((r-r0)/sigma)
        // = phi0/sigma * (1 - tanh((r-r0)/sigma)^2)
        const double t = tanh((r - r0) / sigma);
        return phi0 / sigma * (1.0 - t * t);
    }
    else if constexpr (CC == 2) {
        // d/dr: phi0_scaled * exp(-(r-r0)^2/sigma)
        const double phi0_s = a2 * a2 * phi0 / (1.0 + a2 * a2);
        const double dr = r - r0;
        return -2.0 * phi0_s * dr / sigma * exp(-dr * dr / sigma);
    }
    else if constexpr (CC == 3) {
        // d/dr of the 5-Gaussian sum / l2
        return (-0.005601976461  * (r - 0.18039e-7)  / l2 * exp(-0.0581545  * (r - 0.18039e-7)  * (r - 0.18039e-7)  / l2)
                -0.06649246419  * (r + 0.96741e-8)  / l2 * exp(-0.111412   * (r + 0.96741e-8)  * (r + 0.96741e-8)  / l2)
                -0.1771390956   * (r - 0.109822e-7) / l2 * exp(-0.207156   * (r - 0.109822e-7) * (r - 0.109822e-7) / l2)
                -0.1541602184   * (r + 0.213778e-7) / l2 * exp(-0.37742    * (r + 0.213778e-7) * (r + 0.213778e-7) / l2)
                -0.02961842988  * (r - 0.878608e-7) / l2 * exp(-0.68406    * (r - 0.878608e-7) * (r - 0.878608e-7) / l2)) / l2;
    }
    else if constexpr (CC == 4) {
        // d/dr of two tanh walls
        const double phi0_h = 0.5 * phi0;
        const double tp = tanh((r + r0) / sigma);
        const double tm = tanh((r - r0) / sigma);
        return phi0_h * ((1.0 - tp * tp) / sigma - (1.0 - tm * tm) / sigma);
    }
    else if constexpr (CC == 5) {
        // d/dr: phi0 * exp(-(r-r0)^2/sigma)
        const double dr = r - r0;
        return -2.0 * phi0 * dr / sigma * exp(-dr * dr / sigma);
    }
    return 0.0;
}

// ---------------------------------------------------------------------------
// d_potential<CC> -- f(R) scalar field potential V(phi) at a point
// ---------------------------------------------------------------------------
template<int CC>
__device__ __forceinline__
double d_potential(double X, double Y, double Z,
                   double a2, double r0, double phi0, double sigma, double l2)
{
    if constexpr (CC == 1 || CC == 4 || CC == 5) {
        // Zero potential for these configurations
        return 0.0;
    }
    else if constexpr (CC == 2 || CC == 3) {
        const double PI = acos(-1.0);
        const double sqrt_PI_over_3 = sqrt(PI / 3.0);

        // Get phi at this point (includes CC-specific scaling internally)
        const double v = d_phi<CC>(X, Y, Z, a2, r0, phi0, sigma, l2);

        // V = exp(-8 sqrt(PI/3) phi) * (1 - exp(4 sqrt(PI/3) phi))^2 / (32 PI a2)
        const double e_term = exp(4.0 * sqrt_PI_over_3 * v);
        return exp(-8.0 * sqrt_PI_over_3 * v)
             * (1.0 - e_term) * (1.0 - e_term)
             / (32.0 * PI * a2);
    }
    return 0.0;
}

// =============================================================================
// 3D GPU Kernels (evolution-time grid operations)
// =============================================================================

// ---------------------------------------------------------------------------
// frpotential_kernel<CC>
//   Compute f(R) potential V and its derivative dV/dSphi on a 3D grid.
//   For CC=1,4,5: V=0, dVdSphi=0.
//   For CC=2,3: V = exp(-8 sqrt(PI/3) Sphi) * (1-exp(4 sqrt(PI/3) Sphi))^2 / (32 PI a2)
//               dVdSphi = 1/(12 a2) sqrt(3/PI) exp(-8 sqrt(PI/3) Sphi) * (-1+exp(4 sqrt(PI/3) Sphi))
// ---------------------------------------------------------------------------
template<int CC>
__global__ void frpotential_kernel(
    int ex0, int ex1, int ex2,
    const double* __restrict__ Sphi,
    double* __restrict__ V,
    double* __restrict__ dVdSphi,
    double a2, double r0, double phi0, double sigma, double l2,
    int pex0, int pex1)
{
    const int n = ex0 * ex1 * ex2;
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n) return;

    int i, j, k;
    linear_to_ijk(tid, ex0, ex1, i, j, k);
    const int idx = i + pex0 * (j + pex1 * k);

    if constexpr (CC == 1 || CC == 4 || CC == 5) {
        // Zero potential and derivative
        V[idx] = 0.0;
        dVdSphi[idx] = 0.0;
    }
    else if constexpr (CC == 2 || CC == 3) {
        const double PI = acos(-1.0);
        const double sqrt_PI_over_3 = sqrt(PI / 3.0);

        const double sphi = Sphi[idx];
        const double e4 = exp(4.0 * sqrt_PI_over_3 * sphi);
        const double e8_neg = exp(-8.0 * sqrt_PI_over_3 * sphi);

        // V = exp(-8 sqrt(PI/3) Sphi) * (1 - exp(4 sqrt(PI/3) Sphi))^2 / (32 PI a2)
        V[idx] = e8_neg * (1.0 - e4) * (1.0 - e4) / (32.0 * PI * a2);

        // dVdSphi = 1/(12 a2) * sqrt(3/PI) * exp(-8 sqrt(PI/3) Sphi) * (-1 + exp(4 sqrt(PI/3) Sphi))
        dVdSphi[idx] = (1.0 / (a2 * 12.0)) * sqrt(3.0 / PI) * e8_neg * (-1.0 + e4);
    }
}

// ---------------------------------------------------------------------------
// frfprim_kernel<CC>
//   Compute f'(R) = 1 + 2*a2*R for f(R) = R + a2*R^2 gravity.
//   For CC=1,4,5: fprim = 1.
//   For CC=2,3: fprim = 1 + 2*a2*RR.
// ---------------------------------------------------------------------------
template<int CC>
__global__ void frfprim_kernel(
    int ex0, int ex1, int ex2,
    const double* __restrict__ RR,
    double* __restrict__ fprim,
    double a2, double r0, double phi0, double sigma, double l2,
    int pex0, int pex1)
{
    const int n = ex0 * ex1 * ex2;
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n) return;

    int i, j, k;
    linear_to_ijk(tid, ex0, ex1, i, j, k);
    const int idx = i + pex0 * (j + pex1 * k);

    if constexpr (CC == 1 || CC == 4 || CC == 5) {
        fprim[idx] = 1.0;
    }
    else if constexpr (CC == 2 || CC == 3) {
        fprim[idx] = 1.0 + 2.0 * a2 * RR[idx];
    }
}

// =============================================================================
// 1D GPU Kernels (for ADM initial data setup)
// =============================================================================

// ---------------------------------------------------------------------------
// set_rho_adm1_kernel<CC>
//   rho(i) = potential(X(i), Y(i), Z(i))   -- per-point potential
// ---------------------------------------------------------------------------
template<int CC>
__global__ void set_rho_adm1_kernel(
    int ex0,
    const double* __restrict__ X,
    const double* __restrict__ Y,
    const double* __restrict__ Z,
    double* __restrict__ rho,
    double a2, double r0, double phi0, double sigma, double l2)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= ex0) return;

    rho[tid] = d_potential<CC>(X[tid], Y[tid], Z[tid],
                               a2, r0, phi0, sigma, l2);
}

// ---------------------------------------------------------------------------
// set_rho_adm2_kernel<CC>
//   rho(i) = dphi(X(i), Y(i), Z(i))^2   -- squared radial derivative
// ---------------------------------------------------------------------------
template<int CC>
__global__ void set_rho_adm2_kernel(
    int ex0,
    const double* __restrict__ X,
    const double* __restrict__ Y,
    const double* __restrict__ Z,
    double* __restrict__ rho,
    double a2, double r0, double phi0, double sigma, double l2)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= ex0) return;

    const double dph = d_dphi<CC>(X[tid], Y[tid], Z[tid],
                                  a2, r0, phi0, sigma, l2);
    rho[tid] = dph * dph;
}

// ---------------------------------------------------------------------------
// set_rho_adm_kernel
//   rho (INOUT):  rho -> rho^4 -> rho^8*rho1 + rho^4*rho2
//   Equivalent to Fortran:
//     rho = rho**4
//     rho = rho**2 * rho1 + rho * rho2
// ---------------------------------------------------------------------------
__global__ void set_rho_adm_kernel(
    int ex0,
    double* __restrict__ rho,
    const double* __restrict__ rho1,
    const double* __restrict__ rho2)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= ex0) return;

    // rho -> psi^4
    double psi4 = rho[tid];
    psi4 = psi4 * psi4;   // psi^2
    psi4 = psi4 * psi4;   // psi^4

    // rho = psi^8 * rho1 + psi^4 * rho2
    rho[tid] = psi4 * psi4 * rho1[tid] + psi4 * rho2[tid];
}

// =============================================================================
// Host C++ launch functions (template on CC, used by extern "C" wrappers)
// =============================================================================

template<int CC>
void launch_frpotential(
    cudaStream_t stream,
    int ex[3], int pex[3],
    double* Sphi, double* V, double* dVdSphi,
    double a2, double r0, double phi0, double sigma, double l2)
{
    const int total = ex[0] * ex[1] * ex[2];
    if (total <= 0) return;
    const dim3 block(256);
    const dim3 grid((total + 255) / 256);
    frpotential_kernel<CC><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2],
        Sphi, V, dVdSphi,
        a2, r0, phi0, sigma, l2,
        pex[0], pex[1]);
}

template<int CC>
void launch_frfprim(
    cudaStream_t stream,
    int ex[3], int pex[3],
    double* RR, double* fprim,
    double a2, double r0, double phi0, double sigma, double l2)
{
    const int total = ex[0] * ex[1] * ex[2];
    if (total <= 0) return;
    const dim3 block(256);
    const dim3 grid((total + 255) / 256);
    frfprim_kernel<CC><<<grid, block, 0, stream>>>(
        ex[0], ex[1], ex[2],
        RR, fprim,
        a2, r0, phi0, sigma, l2,
        pex[0], pex[1]);
}

template<int CC>
void launch_set_rho_adm1(
    cudaStream_t stream,
    int ex0,
    const double* X, const double* Y, const double* Z,
    double* rho,
    double a2, double r0, double phi0, double sigma, double l2)
{
    if (ex0 <= 0) return;
    const dim3 block(256);
    const dim3 grid((ex0 + 255) / 256);
    set_rho_adm1_kernel<CC><<<grid, block, 0, stream>>>(
        ex0, X, Y, Z, rho,
        a2, r0, phi0, sigma, l2);
}

template<int CC>
void launch_set_rho_adm2(
    cudaStream_t stream,
    int ex0,
    const double* X, const double* Y, const double* Z,
    double* rho,
    double a2, double r0, double phi0, double sigma, double l2)
{
    if (ex0 <= 0) return;
    const dim3 block(256);
    const dim3 grid((ex0 + 255) / 256);
    set_rho_adm2_kernel<CC><<<grid, block, 0, stream>>>(
        ex0, X, Y, Z, rho,
        a2, r0, phi0, sigma, l2);
}

void launch_set_rho_adm(
    cudaStream_t stream,
    int ex0,
    double* rho,
    const double* rho1, const double* rho2)
{
    if (ex0 <= 0) return;
    const dim3 block(256);
    const dim3 grid((ex0 + 255) / 256);
    set_rho_adm_kernel<<<grid, block, 0, stream>>>(
        ex0, rho, rho1, rho2);
}

// =============================================================================
// Explicit template instantiation macro (CC = 1..5)
//   Uses a helper macro to avoid repeating the long signatures.
//   Macro parameter name (SCALAR_CC) does not clash with C++ identifiers.
// =============================================================================
#define INSTANTIATE_SET_RHO_ADM_FOR_CC(SCALAR_CC) \
template __global__ void frpotential_kernel<SCALAR_CC>( \
    int, int, int, \
    const double* __restrict__, double* __restrict__, double* __restrict__, \
    double, double, double, double, double, \
    int, int); \
template __global__ void frfprim_kernel<SCALAR_CC>( \
    int, int, int, \
    const double* __restrict__, double* __restrict__, \
    double, double, double, double, double, \
    int, int); \
template __global__ void set_rho_adm1_kernel<SCALAR_CC>( \
    int, \
    const double* __restrict__, const double* __restrict__, const double* __restrict__, \
    double* __restrict__, \
    double, double, double, double, double); \
template __global__ void set_rho_adm2_kernel<SCALAR_CC>( \
    int, \
    const double* __restrict__, const double* __restrict__, const double* __restrict__, \
    double* __restrict__, \
    double, double, double, double, double); \
template void launch_frpotential<SCALAR_CC>( \
    cudaStream_t, int[3], int[3], \
    double*, double*, double*, \
    double, double, double, double, double); \
template void launch_frfprim<SCALAR_CC>( \
    cudaStream_t, int[3], int[3], \
    double*, double*, \
    double, double, double, double, double); \
template void launch_set_rho_adm1<SCALAR_CC>( \
    cudaStream_t, int, \
    const double*, const double*, const double*, \
    double*, \
    double, double, double, double, double); \
template void launch_set_rho_adm2<SCALAR_CC>( \
    cudaStream_t, int, \
    const double*, const double*, const double*, \
    double*, \
    double, double, double, double, double)

INSTANTIATE_SET_RHO_ADM_FOR_CC(1);
INSTANTIATE_SET_RHO_ADM_FOR_CC(2);
INSTANTIATE_SET_RHO_ADM_FOR_CC(3);
INSTANTIATE_SET_RHO_ADM_FOR_CC(4);
INSTANTIATE_SET_RHO_ADM_FOR_CC(5);

#undef INSTANTIATE_SET_RHO_ADM_FOR_CC

// =============================================================================
// extern "C" linkage wrappers for Fortran / C interop
//
// These read EScalar_CC from macrodef.fh (preprocessor constant) and dispatch
// to the appropriate template instantiation at compile time.
//
// Host-side: the caller is responsible for calling set_escalar_parameter
// to obtain the 5 parameters (a2, phi0, r0, sigma, l2) and passing them in.
// =============================================================================

extern "C" void gpu_frpotential_launch(
    cudaStream_t stream,
    int ex[3], int pex[3],
    double* Sphi, double* V, double* dVdSphi,
    double a2, double r0, double phi0, double sigma, double l2)
{
    // EScalar_CC is a preprocessor macro from macrodef.fh.
    // Use if constexpr-compatible dispatch. Since all instantiations exist,
    // a plain switch with literal cases is safe and explicit.
    switch (EScalar_CC) {
        case 1: launch_frpotential<1>(stream, ex, pex, Sphi, V, dVdSphi, a2, r0, phi0, sigma, l2); break;
        case 2: launch_frpotential<2>(stream, ex, pex, Sphi, V, dVdSphi, a2, r0, phi0, sigma, l2); break;
        case 3: launch_frpotential<3>(stream, ex, pex, Sphi, V, dVdSphi, a2, r0, phi0, sigma, l2); break;
        case 4: launch_frpotential<4>(stream, ex, pex, Sphi, V, dVdSphi, a2, r0, phi0, sigma, l2); break;
        case 5: launch_frpotential<5>(stream, ex, pex, Sphi, V, dVdSphi, a2, r0, phi0, sigma, l2); break;
    }
}

extern "C" void gpu_frfprim_launch(
    cudaStream_t stream,
    int ex[3], int pex[3],
    double* RR, double* fprim,
    double a2, double r0, double phi0, double sigma, double l2)
{
    switch (EScalar_CC) {
        case 1: launch_frfprim<1>(stream, ex, pex, RR, fprim, a2, r0, phi0, sigma, l2); break;
        case 2: launch_frfprim<2>(stream, ex, pex, RR, fprim, a2, r0, phi0, sigma, l2); break;
        case 3: launch_frfprim<3>(stream, ex, pex, RR, fprim, a2, r0, phi0, sigma, l2); break;
        case 4: launch_frfprim<4>(stream, ex, pex, RR, fprim, a2, r0, phi0, sigma, l2); break;
        case 5: launch_frfprim<5>(stream, ex, pex, RR, fprim, a2, r0, phi0, sigma, l2); break;
    }
}

extern "C" void gpu_set_rho_adm1_launch(
    cudaStream_t stream,
    int ex0,
    const double* X, const double* Y, const double* Z,
    double* rho,
    double a2, double r0, double phi0, double sigma, double l2)
{
    switch (EScalar_CC) {
        case 1: launch_set_rho_adm1<1>(stream, ex0, X, Y, Z, rho, a2, r0, phi0, sigma, l2); break;
        case 2: launch_set_rho_adm1<2>(stream, ex0, X, Y, Z, rho, a2, r0, phi0, sigma, l2); break;
        case 3: launch_set_rho_adm1<3>(stream, ex0, X, Y, Z, rho, a2, r0, phi0, sigma, l2); break;
        case 4: launch_set_rho_adm1<4>(stream, ex0, X, Y, Z, rho, a2, r0, phi0, sigma, l2); break;
        case 5: launch_set_rho_adm1<5>(stream, ex0, X, Y, Z, rho, a2, r0, phi0, sigma, l2); break;
    }
}

extern "C" void gpu_set_rho_adm2_launch(
    cudaStream_t stream,
    int ex0,
    const double* X, const double* Y, const double* Z,
    double* rho,
    double a2, double r0, double phi0, double sigma, double l2)
{
    switch (EScalar_CC) {
        case 1: launch_set_rho_adm2<1>(stream, ex0, X, Y, Z, rho, a2, r0, phi0, sigma, l2); break;
        case 2: launch_set_rho_adm2<2>(stream, ex0, X, Y, Z, rho, a2, r0, phi0, sigma, l2); break;
        case 3: launch_set_rho_adm2<3>(stream, ex0, X, Y, Z, rho, a2, r0, phi0, sigma, l2); break;
        case 4: launch_set_rho_adm2<4>(stream, ex0, X, Y, Z, rho, a2, r0, phi0, sigma, l2); break;
        case 5: launch_set_rho_adm2<5>(stream, ex0, X, Y, Z, rho, a2, r0, phi0, sigma, l2); break;
    }
}

extern "C" void gpu_set_rho_adm_launch(
    cudaStream_t stream,
    int ex0,
    double* rho,
    const double* rho1, const double* rho2)
{
    launch_set_rho_adm(stream, ex0, rho, rho1, rho2);
}
