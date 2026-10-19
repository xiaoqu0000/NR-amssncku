#include "macrodef.h"
#if defined(USE_GPU) && defined(With_AHF)
#include "ahf_gpu.h"
#include "../diff_new_gpu.cuh"
#include "../point_diff_new_sh_gpu.cuh"
#include <stdexcept>
#include <climits>

namespace AHFinderDirect { namespace gpu {
namespace {
struct DerivativeInputs {
    int shape[3],padded[3],symmetry,shell_type;
    const double *coordinates[3],*fields[7],*transform[9];
    double *outputs[21];
};
__global__ void derivative_kernel(DerivativeInputs in){
    int linear=blockIdx.x*blockDim.x+threadIdx.x;
    int n=in.shape[0]*in.shape[1]*in.shape[2];if(linear>=n)return;
    int i=linear%in.shape[0],j=(linear/in.shape[0])%in.shape[1],k=linear/(in.shape[0]*in.shape[1]);
    int index=(k*in.padded[1]+j)*in.padded[0]+i;
    // phi, gxx, gxy, gxz, gyy, gyz, gzz in precisely the CPU order.
    const double parity[7][3]={{1,1,1},{1,1,1},{-1,-1,1},{-1,1,-1},{1,1,1},{1,-1,-1},{1,1,1}};
    for(int v=0;v<7;++v){
        double dx,dy,dz;
        if(in.shell_type<0)d_fderivs_point(in.shape,in.padded,in.fields[v],&dx,&dy,&dz,
            in.coordinates[0],in.coordinates[1],in.coordinates[2],parity[v][0],parity[v][1],parity[v][2],in.symmetry,0,i,j,k);
        else d_point_fderivs_shc_gw<ghost_width>(in.shape,in.padded,in.fields[v],&dx,&dy,&dz,
            in.coordinates[0],in.coordinates[1],in.coordinates[2],parity[v][0],parity[v][1],parity[v][2],in.symmetry,0,in.shell_type,
            in.transform[0],in.transform[1],in.transform[2],in.transform[3],in.transform[4],in.transform[5],in.transform[6],in.transform[7],in.transform[8],i,j,k);
        in.outputs[3*v][index]=dx;in.outputs[3*v+1][index]=dy;in.outputs[3*v+2][index]=dz;
    }
}
}
void derivatives(cudaStream_t stream,const int shape[3],const int padded[3],const double *const coordinates[3],
    const double *const fields[7],double *const outputs[21],int symmetry,int shell_type,const double *const transform[9]){
    if(!shape||!padded||!coordinates||!fields||!outputs||symmetry<0||symmetry>2)
        throw std::runtime_error("AHF invalid derivative inputs");
    DerivativeInputs in{};in.symmetry=symmetry;in.shell_type=shell_type;
    int n=1,np=1;
    for(int a=0;a<3;++a){
        if(shape[a]<2||padded[a]<shape[a]||!coordinates[a]||n>INT_MAX/shape[a]||np>INT_MAX/padded[a])
            throw std::runtime_error("AHF invalid or overflowing derivative shape");
        n*=shape[a];np*=padded[a];
        in.shape[a]=shape[a];in.padded[a]=padded[a];in.coordinates[a]=coordinates[a];
    }
    for(int v=0;v<7;++v){if(!fields[v])throw std::runtime_error("AHF missing derivative field");in.fields[v]=fields[v];}
    for(int v=0;v<21;++v){if(!outputs[v])throw std::runtime_error("AHF missing derivative output");in.outputs[v]=outputs[v];}
    if(shell_type>=0){
        if(!transform)throw std::runtime_error("AHF shell derivatives require Cartesian coordinate Jacobian");
        for(int a=0;a<9;++a){if(!transform[a])throw std::runtime_error("AHF missing shell Jacobian field");in.transform[a]=transform[a];}
    }
    if(n>INT_MAX-127)throw std::runtime_error("AHF derivative launch size overflow");
    derivative_kernel<<<(n+127)/128,128,0,stream>>>(in);check(cudaGetLastError(),"AHF padded point derivatives");
}
}}
#endif
