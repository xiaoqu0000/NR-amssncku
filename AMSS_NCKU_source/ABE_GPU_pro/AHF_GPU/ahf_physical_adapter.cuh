#ifndef AHF_PHYSICAL_ADAPTER_CUH
#define AHF_PHYSICAL_ADAPTER_CUH
#if defined(__CUDACC__)
#include "ahf_pow_compat.cuh"
#endif
__device__ inline bool ahf_physical_point(const double* raw,double* physical) {
#if defined(__CUDA_ARCH__)
    // Preserve original input values even when raw and physical alias.
    double r[35];
    for(int f=0;f<35;++f)r[f]=raw[f];
    for(int f=0;f<35;++f)physical[f]=r[f];
    const double chi=r[24]+1.0;
    if(!(chi>0.0)||!isfinite(chi))return false;
    const auto psi_result=ahf_power::pow_compat(chi,-0.25);
    const auto p125_result=ahf_power::pow_compat(chi,-1.25);
    if(psi_result.status!=ahf_power::Status::ok||p125_result.status!=ahf_power::Status::ok)return false;
    const double psi=psi_result.value,p125=p125_result.value;
    double dp[3];
    for(int a=0;a<3;++a)dp[a]=-p125/4*r[25+a];
    const double psi3=(psi*psi)*psi,psi4=(psi*psi)*(psi*psi),tr3=r[34]/3.0;
    physical[24]=psi;
    for(int a=0;a<3;++a)physical[25+a]=dp[a];
    for(int j=0;j<6;++j){
        const double g=(j==0||j==3||j==5)?r[4*j]+1.0:r[4*j];
        physical[4*j]=g*psi4;
        for(int a=0;a<3;++a)physical[4*j+1+a]=4.0*psi3*dp[a]*g+psi4*r[4*j+1+a];
        physical[28+j]=psi4*(r[28+j]+g*tr3);
    }
    for(int f=0;f<35;++f)if(!isfinite(physical[f]))return false;
    return true;
#else
    // nvcc host parsing only: never substitute host or CUDA-native powers.
    return false;
#endif
}
#endif
