#include "macrodef.h"
#if defined(USE_GPU) && defined(With_AHF)
#include "ahf_gpu.h"
#include "ahf_ilucg_gpu.h"
#include "../AHF/driver.h"
#include "ahf_native_geometry.cuh"
#include "ahf_physical_adapter.cuh"
#include <cmath>
#include <cstdio>
#include <stdexcept>

namespace AHFinderDirect { namespace gpu {
namespace {
// Share the validated ORDER4 expression with the diagnostic implementation.
using ahf_diagnostics::DevicePlan;
__device__ void angular(DeviceSurface s,int q,double a[5]) {
    ahf_diagnostics::DeviceSurface view;
    view.stencils=s.diagnostic_plan.stencils;view.points=s.diagnostic_plan.points;view.ghost_h=s.ghost_h;
    // Defined below in a separate namespace to retain the exact accepted source.
    extern __device__ void original_angular(ahf_diagnostics::DeviceSurface,int,double*);
    original_angular(view,q,a);
}
namespace stencil {
using DeviceSurface=ahf_diagnostics::DeviceSurface;
#include "ahf_order4_stencil.cuh"
}
__device__ void original_angular(ahf_diagnostics::DeviceSurface s,int q,double *a){stencil::angular(s,q,a);}
__global__ void sync_kernel(DeviceSurface s){
    if(blockIdx.x||threadIdx.x)return;
    for(int q=0;q<s.n;++q)s.ghost_h[s.nominal_ghost[q]]=s.h[q];
    const auto p=s.diagnostic_plan;
    for(int i=0;i<p.no;++i){auto op=p.operations[i];double v=0;
        for(int j=op.begin;j<op.end;++j)v=__dadd_rn(v,__dmul_rn(p.weights[j],s.ghost_h[p.sources[j]]));
        s.ghost_h[op.destination]=v;
    }
}
__global__ void coordinate_kernel(DeviceSurface s){
    int q=blockIdx.x*blockDim.x+threadIdx.x;if(q>=s.n)return;
    auto m=s.metadata[q];for(int a=0;a<3;++a)m.origin[a]=s.origin[a];
    double xyz[3];
    if(!ahf_native::native_geometry(m,s.h[q],s.runtime_geometry[q],xyz)){
        for(int a=0;a<3;++a)s.xyz[a*s.n+q]=NAN;return;
    }
    for(int a=0;a<3;++a)s.xyz[a*s.n+q]=xyz[a];
}
__global__ void init_kernel(DeviceSurface s,double cx,double cy,double cz,double a,double b,double c){
    int q=blockIdx.x*blockDim.x+threadIdx.x;if(q>=s.n)return;
    const double xcos=s.initial_direction[3*q],ycos=s.initial_direction[3*q+1],zcos=s.initial_direction[3*q+2];
    const double AU=cx-s.origin[0],BV=cy-s.origin[1],CW=cz-s.origin[2];
    double r_plus,r_minus;
    // Original initial_guess.C generated ellipsoid expression, unchanged order.
    double t1,t2,t3,t5,t6,t7,t9,t10,t12,t28,t30,t33,t35,t36,t40,t42,t43,t48,t49,t52,t55;
    t1=a*a;t2=b*b;t3=t1*t2;t5=t3*zcos*CW;t6=c*c;t7=t1*t6;
    t9=t7*ycos*BV;t10=t2*t6;t12=t10*xcos*AU;t28=xcos*xcos;t30=CW*CW;
    t33=BV*BV;t35=t10*t28;t36=ycos*ycos;t40=AU*AU;t42=t7*t36;t43=zcos*zcos;t48=t3*t43;
    t49=-2.0*t1*zcos*CW*ycos*BV-2.0*t2*zcos*CW*xcos*AU-2.0*t6*ycos*BV*xcos*AU+t2*t28*t30+t6*t28*t33-t35+t1*t36*t30+t6*t36*t40-t42+t1*t43*t33+t2*t43*t40-t48;
    t52=sqrt(-t3*t6*t49);t55=1/(t35+t42+t48);
    r_plus=(t5+t9+t12+t52)*t55;r_minus=(t5+t9+t12-t52)*t55;
    s.h[q]=r_plus>0&&r_minus<0?r_plus:r_plus<0&&r_minus>0?r_minus:NAN;
}
__global__ void transform_kernel(int n,const double *raw,double *out){
    int q=blockIdx.x*blockDim.x+threadIdx.x;if(q>=n)return;
    // A rejected power/domain must not leave a partly valid physical point.
    if(!ahf_physical_point(raw+35*q,out+35*q))for(int j=0;j<35;++j)out[35*q+j]=NAN;
}
__global__ void theta_kernel(DeviceSurface s,double add,bool coeff){
    int q=blockIdx.x*blockDim.x+threadIdx.x;if(q>=s.n)return;
    double a[5],trace[51],cf[5];angular(s,q,a);
    auto m=s.metadata[q];for(int j=0;j<3;++j)m.origin[j]=s.origin[j];
    double xyz[3];
    bool valid=ahf_native::native_geometry(m,s.h[q],s.runtime_geometry[q],xyz);
    if(valid){
        const double x=xyz[0],y=xyz[1],z=xyz[2];
        double *scratch=s.scratch6+6*q;
        scratch[0]=x*x;scratch[1]=x*y;scratch[2]=x*z;scratch[3]=y*y;scratch[4]=y*z;scratch[5]=z*z;
        valid=ahf_native::original_scalar(s.runtime_geometry[q],a,s.geometry+35*q,add,trace,cf);
    }
    if(!valid){s.theta[q]=NAN;atomicExch(s.status,1);return;}
    s.theta[q]=trace[50];if(coeff)for(int j=0;j<5;++j)s.coefficients[5*q+j]=cf[j];
}
constexpr int reduction_threads=256;
__global__ void norm_kernel(DeviceSurface s){
    __shared__ double maxima[reduction_threads],sums[reduction_threads];
    __shared__ int invalid[reduction_threads];
    const int t=threadIdx.x;
    double inf=0,mean=0;int bad=0;
    for(int q=t;q<s.n;q+=reduction_threads){
        if(!isfinite(s.theta[q])||!isfinite(s.h[q])||s.h[q]<=0)bad=1;
        inf=fmax(inf,fabs(s.theta[q]));mean=__dadd_rn(mean,s.h[q]);
    }
    maxima[t]=inf;sums[t]=mean;invalid[t]=bad;
    __syncthreads();
    for(int stride=reduction_threads/2;stride;stride/=2){
        if(t<stride){
            maxima[t]=fmax(maxima[t],maxima[t+stride]);
            sums[t]=__dadd_rn(sums[t],sums[t+stride]);
            invalid[t]|=invalid[t+stride];
        }
        __syncthreads();
    }
    if(t==0){
        s.scalars[0]=maxima[0];s.scalars[1]=__ddiv_rn(sums[0],double(s.n));
        s.scalars[2]=invalid[0]?1:*s.status;
    }
}
__global__ void radial_kernel(DeviceSurface s){
    int q=blockIdx.x*blockDim.x+threadIdx.x;if(q>=s.n)return;
    for(int j=s.angular.rows[q];j<s.angular.rows[q+1];++j)if(s.angular.columns[j]==q)s.values[j]+=(s.theta[q]-s.saved_theta[q])/1e-6;
    s.theta[q]=s.saved_theta[q];
}
__global__ void signature_scratch_kernel(DeviceSurface s){
    int q=blockIdx.x*blockDim.x+threadIdx.x;if(q>=s.n)return;
    for(int j=0;j<4;++j)s.scratch6[6*q+j]=s.gauge[4*q+j];
    // yz retains the last expansion moment; zz is overwritten by signature.
    s.scratch6[6*q+5]=s.signature[q];
}
__global__ void signature_classify_kernel(DeviceSurface s){
    __shared__ int pole_counts[reduction_threads],other_counts[reduction_threads];
    const int t=threadIdx.x;
    int poles=0,other=0;
    for(int q=t;q<s.n;q+=reduction_threads)if(!isfinite(s.signature[q])){
        const auto g=s.runtime_geometry[q];
        if(isnan(s.signature[q])&&g.xyz[0]==0.0&&g.xyz[1]==0.0)++poles;else ++other;
    }
    pole_counts[t]=poles;other_counts[t]=other;
    __syncthreads();
    for(int stride=reduction_threads/2;stride;stride/=2){
        if(t<stride){
            pole_counts[t]+=pole_counts[t+stride];other_counts[t]+=other_counts[t+stride];
        }
        __syncthreads();
    }
    if(t==0){s.scalars[0]=pole_counts[0];s.scalars[1]=other_counts[0];}
}
__global__ void recenter_kernel(DeviceSurface s){
    __shared__ bool centered;
    const int t=threadIdx.x;
    if(t==0){
        centered=true;
        for(int j=0;j<3;++j)centered=centered&&fabs(s.diagnostics[13+2*j]+s.diagnostics[12+2*j]-2.0*s.origin[j])<.45*(s.diagnostics[13+2*j]-s.diagnostics[12+2*j]);
    }
    __syncthreads();
    if(!centered){
        for(int q=t;q<s.n;q+=reduction_threads){
            double x=s.xyz[q]-s.diagnostics[0],y=s.xyz[s.n+q]-s.diagnostics[1],z=s.xyz[2*s.n+q]-s.diagnostics[2];
            s.h[q]=sqrt(x*x+y*y+z*z);
        }
    }
    __syncthreads();
    if(t==0&&!centered)for(int j=0;j<3;++j)s.origin[j]=s.diagnostics[j];
}
}
void initialize(Surface &s,const double c[3],const double r[3]){
    init_kernel<<<(s.d.n+127)/128,128,0,s.stream>>>(s.d,c[0],c[1],c[2],r[0],r[1],r[2]);check(cudaGetLastError(),"AHF original initial guess");
}
void save_old_shape(Surface &s){check(cudaMemcpyAsync(s.d.oldh,s.d.h,s.d.n*sizeof(double),cudaMemcpyDeviceToDevice,s.stream),"AHF oldh");}
void synchronize(Surface &s){sync_kernel<<<1,1,0,s.stream>>>(s.d);check(cudaGetLastError(),"AHF ordered ghost sync");}
void coordinates(Surface &s){coordinate_kernel<<<(s.d.n+127)/128,128,0,s.stream>>>(s.d);check(cudaGetLastError(),"AHF native coordinates");}
void transform_geometry(cudaStream_t st,int n,const double *raw,double *out){transform_kernel<<<(n+63)/64,64,0,st>>>(n,raw,out);check(cudaGetLastError(),"AHF original physical geometry");}
void evaluate(Surface &s,double add,bool coeff){
    s.exact_pole_signature_nans=0;
    check(cudaMemsetAsync(s.d.status,0,sizeof(int),s.stream),"AHF status reset");synchronize(s);
    theta_kernel<<<(s.d.n+63)/64,64,0,s.stream>>>(s.d,add,coeff);check(cudaGetLastError(),"AHF original scalar expansion");
}
void norms(Surface &s,double result[3]){
    norm_kernel<<<1,reduction_threads,0,s.stream>>>(s.d);check(cudaGetLastError(),"AHF norms kernel");
    check(cudaMemcpyAsync(result,s.d.scalars,3*sizeof(double),cudaMemcpyDeviceToHost,s.stream),"AHF norms");
    check(cudaStreamSynchronize(s.stream),"AHF norms completion");transfer_stats().scalar_d2h+=3*sizeof(double);
    if(result[2])std::fprintf(stderr,"AHF unsupported/nonfinite native geometry, physical conversion, or scalar domain; no CPU fallback\n");
}
void assemble(Surface &s){
    check(ahf_native_baseline_mean(s.stream,s.d.n,s.d.ng,s.d.nominal_ghost,s.d.ghost_h,s.d.newton),"AHF baseline mean");
    assemble_ordered(s.stream,{s.d.n,s.d.instruction_rows,s.d.instructions},s.d.coefficients,s.d.values);
}
void perturb(Surface &s){
    check(cudaMemcpyAsync(s.d.saved_theta,s.d.theta,s.d.n*sizeof(double),cudaMemcpyDeviceToDevice,s.stream),"AHF save Theta");
    check(ahf_native_radial_plus(s.stream,s.d.n,s.d.ng,s.d.nominal_ghost,s.d.h,s.d.ghost_h,s.d.newton),"AHF original radial plus");
}
void restore_and_radial(Surface &s,bool ok){
    // Original failed radial expansion returns without restoring the perturbed h.
    if(!ok)return;
    radial_kernel<<<(s.d.n+127)/128,128,0,s.stream>>>(s.d);check(cudaGetLastError(),"AHF radial diagonal");
    check(ahf_native_radial_restore(s.stream,s.d.n,s.d.ng,s.d.nominal_ghost,s.d.h,s.d.ghost_h,s.d.newton),"AHF original current-ghost restore");
    // No synchronize, coordinate rebuild, saved-h bitcopy or base-field restore here.
}
int linear_step(Surface &s){
    AhfNewtonState state;
    if(cudaMemcpyAsync(&state,s.d.newton,sizeof(state),cudaMemcpyDeviceToHost,s.stream)!=cudaSuccess)
        return GPU_AHF_LINEAR_CUDA_FAILURE;
    if(cudaStreamSynchronize(s.stream)!=cudaSuccess)return GPU_AHF_LINEAR_CUDA_FAILURE;
    transfer_stats().scalar_d2h+=sizeof(state);
    if(state.status||state.phase!=3)return GPU_AHF_LINEAR_NEWTON_STATE;
    return gpu_ahf_linear_solve(s.linear_solver,s.d.values,s.d.theta,s.d.delta);
}
void update(Surface &s){check(ahf_native_newton_update(s.stream,s.d.n,s.d.ng,s.d.nominal_ghost,s.d.h,s.d.ghost_h,s.d.delta,s.d.newton),"AHF original Newton update");}
void diagnostics(Surface &s){
    ahf_diagnostics::compute_from_scratch(s.d.diagnostic_plan,s.d.h,s.d.origin,s.d.geometry,s.d.scratch6,s.d.diagnostics,s.stream);
}
void signature(Surface &s,double dt){
    ahf_diagnostics::signature(s.d.diagnostic_plan,s.d.h,s.d.geometry,s.d.oldh,s.d.gauge,dt,s.d.signature,s.stream);
    signature_scratch_kernel<<<(s.d.n+127)/128,128,0,s.stream>>>(s.d);
    signature_classify_kernel<<<1,reduction_threads,0,s.stream>>>(s.d);check(cudaGetLastError(),"AHF original signature scratch/classification");
    double counts[2];check(cudaMemcpyAsync(counts,s.d.scalars,sizeof(counts),cudaMemcpyDeviceToHost,s.stream),"AHF signature classification");
    check(cudaStreamSynchronize(s.stream),"AHF signature classification completion");transfer_stats().scalar_d2h+=sizeof(counts);
    s.exact_pole_signature_nans=static_cast<int>(counts[0]);
    if(counts[0])std::fprintf(stderr,"AHF reference defect: %.0f exact-pole signature NaNs retained; not a finite-parity pass\n",counts[0]);
    if(counts[1])std::fprintf(stderr,"AHF unclassified nonfinite signature: %.0f points; not a finite-parity pass\n",counts[1]);
    // Classification must not suppress the original writer or diagnostics sequence.
}
void recenter(Surface &s){
    recenter_kernel<<<1,reduction_threads,0,s.stream>>>(s.d);check(cudaGetLastError(),"AHF recenter");
    double origin[3];check(cudaMemcpyAsync(origin,s.d.origin,sizeof(origin),cudaMemcpyDeviceToHost,s.stream),"AHF origin metadata output");
    check(cudaStreamSynchronize(s.stream),"AHF recenter completion");transfer_stats().scalar_d2h+=sizeof(origin);
    s.reference->recentering(origin[0],origin[1],origin[2]);
}
}}
#endif
