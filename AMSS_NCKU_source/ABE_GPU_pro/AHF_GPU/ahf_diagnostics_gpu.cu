#include "ahf_diagnostics_gpu.h"
#include <cmath>
#include <stdexcept>
#include <string>
namespace ahf_diagnostics {
namespace {
#include "ahf_order4_stencil.cuh"
#include "ahf_original_arithmetic.cuh"
constexpr double pi=3.14159265358979323846;
void check(cudaError_t e,const char *op){if(e!=cudaSuccess)throw std::runtime_error(std::string(op)+": "+cudaGetErrorString(e));}
template<class T> void alloc(T *&p,size_t n){check(cudaMalloc(reinterpret_cast<void**>(&p),n*sizeof(T)),"allocate diagnostics scratch");}
template<class T> void upload(T *&p,const std::vector<T> &v){if(v.empty())return;alloc(p,v.size());check(cudaMemcpy(p,v.data(),v.size()*sizeof(T),cudaMemcpyHostToDevice),"upload static diagnostics plan");}
__device__ void local_xyz(const Point &m,double r,double &x,double &y,double &z){
 if(m.family==0){z=m.sign*r*m.normalizer;x=m.second*z;y=m.first*z;}
 else if(m.family==1){x=m.sign*r*m.normalizer;z=m.first*x;y=m.second*x;}
 else{y=m.sign*r*m.normalizer;z=m.first*y;x=m.second*y;}
}
__device__ void angular_ud(const Point &m,double x,double y,double z,double ud[2][3]){
 for(int i=0;i<6;++i)ud[i/3][i%3]=0.;
 if(m.family==1){ud[0][0]=z/(x*x+z*z);ud[0][2]=-x/(x*x+z*z);}
 else{ud[0][1]=z/(y*y+z*z);ud[0][2]=-y/(y*y+z*z);}
 if(m.family==0){ud[1][0]=z/(x*x+z*z);ud[1][2]=-x/(x*x+z*z);}
 else{ud[1][0]=-y/(x*x+y*y);ud[1][1]=x/(x*x+y*y);}
}
__global__ void seed(DevicePlan p,const double *h){int q=blockIdx.x*blockDim.x+threadIdx.x;if(q<p.n)p.ghost_h[p.points[q].ghost]=h[q];}
__global__ void ghosts(DevicePlan p){if(blockIdx.x||threadIdx.x)return;for(int i=0;i<p.no;++i){auto op=p.operations[i];double v=0.;for(int j=op.begin;j<op.end;++j)v=__dadd_rn(v,__dmul_rn(p.weights[j],p.ghost_h[p.sources[j]]));p.ghost_h[op.destination]=v;}}
// work per point: global xyz, local xyz, a_rho,a_sigma, sqrt(abs(J)),g_rr,g_ss.
__global__ void points(DevicePlan p,const double *h,const double *origin,const double *physical){
 int q=blockIdx.x*blockDim.x+threadIdx.x;if(q>=p.n)return;
 const Point m=p.points[q];double *w=p.work+11*q;double x,y,z;local_xyz(m,h[q],x,y,z);
 w[0]=origin[0]+x;w[1]=origin[1]+y;w[2]=origin[2]+z;w[3]=x;w[4]=y;w[5]=z;
 DeviceSurface s;s.stencils=p.stencils;s.points=p.points;s.ghost_h=p.ghost_h;double a[5];angular(s,q,a);w[6]=a[0];w[7]=a[1];
 double v[9];if(m.family==0)tangent_0(m,h[q],x,y,z,v);else if(m.family==1)tangent_1(m,h[q],x,y,z,v);else tangent_2(m,h[q],x,y,z,v);
 const double *g=physical+35*q;double grr,grs,gss;
 double jac=induced_metric(v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],a[0],a[1],g[0],g[4],g[8],g[12],g[16],g[20],grr,grs,gss);
 w[8]=sqrt(fabs(jac));w[9]=grr;w[10]=gss;
}
template<bool FromScratch=false>
__device__ double fn_value(DevicePlan p,const double *h,int q,int field,const double *scratch=nullptr){
 const double *w=p.work+11*q;
 if(field==0)return 1.;if(field==1)return h[q];if(field<5)return w[field-2];
 if constexpr(FromScratch)return scratch[6*q+field-5];
 switch(field){case 5:return w[0]*w[0];case 6:return w[0]*w[1];case 7:return w[0]*w[2];case 8:return w[1]*w[1];case 9:return w[1]*w[2];default:return w[2]*w[2];}
}
__device__ double symmetry(double integral,int sym,bool z,bool y,bool x){
 if(sym==1)integral*=z?2.:0.;
 if(sym==2){integral*=z?2.:0.;integral*=y?2.:0.;integral*=x?2.:0.;}
 return integral;
}
template<bool FromScratch=false>
__device__ double integrate(DevicePlan p,const double *h,int field,const double *scratch=nullptr){
 double integral=0.;for(int pn=0;pn<p.np;++pn){const Patch patch=p.patches[pn];double sum=0.;
  for(int q=patch.begin;q<patch.end;++q){const Point m=p.points[q];sum+=m.cr*m.cs*fn_value<FromScratch>(p,h,q,field,scratch)*p.work[11*q+8];}
  integral+=patch.dr*patch.ds*sum;
 }
 const bool even[11][3]={{1,1,1},{1,1,1},{1,1,0},{1,0,1},{0,1,1},{1,1,1},{1,0,0},{0,1,0},{1,1,1},{0,0,1},{1,1,1}};
 return symmetry(integral,p.symmetry,even[field][0],even[field][1],even[field][2]);
}
template<bool FromScratch=false>
__global__ void moments(DevicePlan p,const double *h,const double *origin,double *d,const double *scratch=nullptr){
 if(blockIdx.x||threadIdx.x)return;
 double hrmin=INFINITY,hrmax=-INFINITY,lo[3]={INFINITY,INFINITY,INFINITY},hi[3]={-INFINITY,-INFINITY,-INFINITY};
 for(int q=0;q<p.n;++q){double r=fabs(h[q]);if(r<hrmin)hrmin=r;if(r>hrmax)hrmax=r;for(int a=0;a<3;++a){double v=p.work[11*q+a];if(v<lo[a])lo[a]=v;if(v>hi[a])hi[a]=v;}}
 d[9]=hrmin;d[10]=hrmax;
 if(p.symmetry==1)lo[2]=origin[2]-(hi[2]-origin[2]);
 if(p.symmetry==2)for(int a=0;a<3;++a)lo[a]=origin[a]-(hi[a]-origin[a]);
 for(int a=0;a<3;++a){d[12+2*a]=lo[a];d[13+2*a]=hi[a];}
 double integral[11];for(int i=0;i<11;++i)integral[i]=integrate<FromScratch>(p,h,i,scratch);
 for(int a=0;a<3;++a)d[a]=integral[2+a]/integral[0];
 d[3]=integral[5]/integral[0]-d[0]*d[0];d[4]=integral[6]/integral[0]-d[0]*d[1];d[5]=integral[7]/integral[0]-d[0]*d[2];
 d[6]=integral[8]/integral[0]-d[1]*d[1];d[7]=integral[9]/integral[0]-d[1]*d[2];d[8]=integral[10]/integral[0]-d[2]*d[2];
 d[11]=integral[1]/integral[0];d[21]=integral[0];d[22]=sqrt(d[21]/(16.*pi));d[23]=sqrt(d[21]/(4.*pi));
 double arc_total[3]={0.,0.,0.};
 for(int i=0;i<p.na;++i){const Arc arc=p.arcs[i];double sum=0.;for(int j=arc.begin;j<arc.end;++j){int q=p.arc_points[j];sum+=p.arc_coefficients[j]*sqrt(p.work[11*q+9+arc.axis]);}arc_total[arc.plane]+=arc.delta*sum;}
 for(int a=0;a<3;++a){if(p.symmetry==1)arc_total[a]*=a==0?1.:2.;if(p.symmetry==2)arc_total[a]*=4.;d[18+a]=arc_total[a];}
}
__global__ void ps_points(DevicePlan p,const double *h,const double *physical,const double *d){
 int q=blockIdx.x*blockDim.x+threadIdx.x;if(q>=p.n)return;double *w=p.work+11*q;double ud[2][3];angular_ud(p.points[q],w[3],w[4],w[5],ud);double a[2]={w[6],w[7]};
 original_ps(h[q],w[3],w[4],w[5],w,d,a,physical+35*q,ud,p.integrals+6*q);
}
__global__ void ps_reduce(DevicePlan p,double *d){
 if(blockIdx.x||threadIdx.x)return;const bool even[6][3]={{1,1,0},{1,0,1},{0,1,1},{0,0,1},{0,1,0},{1,0,0}};
 for(int f=0;f<6;++f){double integral=0.;for(int pn=0;pn<p.np;++pn){const Patch patch=p.patches[pn];double sum=0.;for(int q=patch.begin;q<patch.end;++q){const Point m=p.points[q];sum+=m.cr*m.cs*p.integrals[6*q+f]*p.work[11*q+8];}integral+=patch.dr*patch.ds*sum;}integral=symmetry(integral,p.symmetry,even[f][0],even[f][1],even[f][2]);const double F1o8pi=1.0/8/pi;d[24+f]=integral*F1o8pi;}
}
__global__ void signature_kernel(DevicePlan p,const double *h,const double *physical,const double *oldh,const double *gauge,double dt,double *out){
 int q=blockIdx.x*blockDim.x+threadIdx.x;if(q>=p.n)return;double x,y,z;local_xyz(p.points[q],h[q],x,y,z);out[q]=original_signature(h[q],x,y,z,physical+35*q,gauge+4*q,oldh[q],dt);
}
}
DevicePlan upload_plan(const HostPlan &h){
 DevicePlan p;p.n=h.n;p.ng=h.ng;p.symmetry=h.symmetry;p.np=h.patches.size();p.no=h.operations.size();p.na=h.arcs.size();
 try{upload(p.points,h.points);upload(p.stencils,h.stencils);upload(p.sources,h.sources);upload(p.weights,h.weights);upload(p.operations,h.operations);upload(p.patches,h.patches);upload(p.arcs,h.arcs);upload(p.arc_points,h.arc_points);upload(p.arc_coefficients,h.arc_coefficients);alloc(p.ghost_h,p.ng);alloc(p.work,11*p.n);alloc(p.integrals,6*p.n);}catch(...){release_plan(p);throw;}return p;
}
void release_plan(DevicePlan &p) noexcept{cudaFree(p.points);cudaFree(p.stencils);cudaFree(p.sources);cudaFree(p.weights);cudaFree(p.operations);cudaFree(p.patches);cudaFree(p.arcs);cudaFree(p.arc_points);cudaFree(p.arc_coefficients);cudaFree(p.ghost_h);cudaFree(p.work);cudaFree(p.integrals);p=DevicePlan{};}
void compute(DevicePlan &p,const double *h,const double *origin,const double *physical,double *slots,cudaStream_t stream){
 if(p.n<=0||!h||!origin||!physical||!slots)throw std::invalid_argument("Missing diagnostics device inputs");int blocks=(p.n+127)/128;
 check(cudaMemsetAsync(p.ghost_h,0,p.ng*sizeof(double),stream),"initialize ghost scratch");seed<<<blocks,128,0,stream>>>(p,h);ghosts<<<1,1,0,stream>>>(p);points<<<blocks,128,0,stream>>>(p,h,origin,physical);moments<<<1,1,0,stream>>>(p,h,origin,slots);ps_points<<<blocks,128,0,stream>>>(p,h,physical,slots);ps_reduce<<<1,1,0,stream>>>(p,slots);check(cudaGetLastError(),"launch diagnostics kernels");
}
void compute_from_scratch(DevicePlan &p,const double *h,const double *origin,const double *physical,const double *scratch,double *slots,cudaStream_t stream){
 if(p.n<=0||!h||!origin||!physical||!scratch||!slots)throw std::invalid_argument("Missing actualsequence device inputs");int blocks=(p.n+127)/128;
 check(cudaMemsetAsync(p.ghost_h,0,p.ng*sizeof(double),stream),"initialize sequence ghost scratch");seed<<<blocks,128,0,stream>>>(p,h);ghosts<<<1,1,0,stream>>>(p);points<<<blocks,128,0,stream>>>(p,h,origin,physical);moments<true><<<1,1,0,stream>>>(p,h,origin,slots,scratch);ps_points<<<blocks,128,0,stream>>>(p,h,physical,slots);ps_reduce<<<1,1,0,stream>>>(p,slots);check(cudaGetLastError(),"launch actualsequence diagnostics");
}
void signature(DevicePlan &p,const double *h,const double *physical,const double *oldh,const double *gauge,double dt,double *out,cudaStream_t stream){
 if(p.n<=0||!h||!physical||!oldh||!gauge||!out||!std::isfinite(dt)||dt<=0)throw std::invalid_argument("Missing signature inputs or invalid dt");signature_kernel<<<(p.n+127)/128,128,0,stream>>>(p,h,physical,oldh,gauge,dt,out);check(cudaGetLastError(),"launch original signature");
}
}
