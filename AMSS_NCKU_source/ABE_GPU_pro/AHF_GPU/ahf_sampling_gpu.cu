#include "macrodef.h"
#if defined(USE_GPU) && defined(With_AHF)
#include "ahf_gpu.h"
#include "../bssn_gpu_class.h"
#include "../MPatch.h"
#include "../gpu_fusion.h"
#include "../fmisc_gpu.cuh"
#include "../shellfunctions.h"
#include "ahf_neville_interp.cuh"
#include "ahf_device_atan.cuh"
#include <mpi.h>
#include <vector>
#include <algorithm>
#include <climits>
#include <cfloat>
#include <cstdlib>
#include <cstdio>
#include <stdexcept>

namespace AHFinderDirect { namespace gpu {
namespace {
constexpr int max_fields = 64;   // 度量 35 场 + 规范 4 场
struct GridBlock {
    int rank, shell_type, shape[3], padded[3];
    double patch_lower[3], patch_upper[3], lower[3], upper[3], spacing[3];
    double *coordinates[3];
};
struct Buffers {
    cudaStream_t stream=nullptr;
    int capacity=0, block_capacity=0;
    double *xyz=nullptr,*local_xyz=nullptr,*sample=nullptr,*soa=nullptr;
    int *selection=nullptr,*bad=nullptr;
    GridBlock *blocks=nullptr;
    double **field_ptrs=nullptr;
    cudaEvent_t *events=nullptr;
    ~Buffers(){release();}
    void release(){
        if(stream)cudaStreamSynchronize(stream);
        cudaFree(xyz);cudaFree(local_xyz);cudaFree(sample);cudaFree(selection);
        cudaFree(bad);cudaFree(blocks);cudaFree(field_ptrs);cudaFree(soa);
        if(events){for(int b=0;b<block_capacity;++b)if(events[b])cudaEventDestroy(events[b]);std::free(events);}
        if(stream)cudaStreamDestroy(stream);
        stream=nullptr;capacity=block_capacity=0;
        xyz=local_xyz=sample=soa=nullptr;selection=bad=nullptr;blocks=nullptr;field_ptrs=nullptr;events=nullptr;
    }
    template<class T> void allocate(T *&p,std::size_t n){check(cudaMalloc(reinterpret_cast<void**>(&p),n*sizeof(T)),"AHF collective allocation");}
    void reserve(int n,int nb){
        if(n<=capacity&&nb<=block_capacity)return;
        release();check(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking),"AHF sampling stream");
        capacity=std::max(n,1);block_capacity=std::max(nb,1);
        allocate(xyz,3*capacity);allocate(local_xyz,3*capacity);allocate(sample,35*capacity);
        allocate(selection,capacity);allocate(bad,1);allocate(blocks,block_capacity);
        allocate(field_ptrs,std::size_t(block_capacity)*max_fields);
        allocate(soa,3*max_fields);
        events=static_cast<cudaEvent_t*>(std::calloc(block_capacity,sizeof(cudaEvent_t)));
        if(!events)throw std::runtime_error("AHF sampling event allocation");
        for(int b=0;b<block_capacity;++b)check(cudaEventCreateWithFlags(events+b,cudaEventDisableTiming),"AHF sampling event");
    }
};
Buffers buffers;

struct RadialMap { double A=1,B=0,r0=0,eps=1; };
__device__ void ahf_shell_map(double x,double y,double z,int &face,double p[3],RadialMap m){
    const double r=sqrt(x*x+y*y+z*z);
    const double f=m.A*(r-m.r0)+m.B*sqrt(1+(r-m.r0)*(r-m.r0)/m.eps);
    p[2]=f+m.A*m.r0-m.B*sqrt(1+m.r0*m.r0/m.eps);
    if(fabs(x)<=z&&fabs(y)<=z){face=0;p[0]=ahf_device_atan(x/z);p[1]=ahf_device_atan(y/z);}
    else if(fabs(x)<=-z&&fabs(y)<=-z){face=1;p[0]=ahf_device_atan(x/z);p[1]=ahf_device_atan(y/z);}
    else if(fabs(y)<=x&&fabs(z)<=x){face=2;p[0]=ahf_device_atan(y/x);p[1]=ahf_device_atan(z/x);}
    else if(fabs(y)<=-x&&fabs(z)<=-x){face=3;p[0]=ahf_device_atan(y/x);p[1]=ahf_device_atan(z/x);}
    else if(fabs(x)<=y&&fabs(z)<=y){face=4;p[0]=ahf_device_atan(x/y);p[1]=ahf_device_atan(z/y);}
    else if(fabs(x)<=-y&&fabs(z)<=-y){face=5;p[0]=ahf_device_atan(x/y);p[1]=ahf_device_atan(z/y);}
    else{face=-1;p[0]=p[1]=NAN;}
}
__global__ void select_kernel(int n,int nb,const GridBlock *blocks,const double *xyz,
    double *local,int *selected,int *bad,double rmin,double rmax,RadialMap radial){
    int q=blockIdx.x*blockDim.x+threadIdx.x;if(q>=n)return;
    double p[3]={xyz[q],xyz[n+q],xyz[2*n+q]};
    selected[q]=-1;
    bool finite=isfinite(p[0])&&isfinite(p[1])&&isfinite(p[2]);
    if(!finite){atomicExch(bad,1);return;}
    int sst=-1;double sp[3]={0,0,0};bool mapped=false;
    // Metadata is ordered finest level, first patch, first block; every rank
    // traverses identical metadata and selects exactly one globally owned block.
    for(int b=0;b<nb;++b){
        GridBlock g=blocks[b];bool inside=true;const double *point=p;
        if(g.shell_type>=0){
            if(!mapped){ahf_shell_map(p[0],p[1],p[2],sst,sp,radial);mapped=true;}
            if(!isfinite(sp[0])||!isfinite(sp[1])||!isfinite(sp[2])){atomicExch(bad,1);return;}
            if(sst!=g.shell_type||sp[2]<rmin||sp[2]>rmax)continue;
            point=sp;
        }
        for(int a=0;a<3;++a)inside=inside&&point[a]>=g.patch_lower[a]&&point[a]<=g.patch_upper[a]
            &&!(point[a]-g.lower[a]<-g.spacing[a]/2||point[a]-g.upper[a]>g.spacing[a]/2);
        if(!inside)continue;
        selected[q]=b;
        for(int a=0;a<3;++a)local[a*n+q]=point[a];
        return;
    }
    atomicExch(bad,1);
}
// 模板区间与三向权重每点只算一次，全部场共享；不用 ya[ORDN^3] 立方缓冲（避免溢出 local memory）。
__device__ __forceinline__ void fused_cartesian_point(
    const GridBlock &g,const double *soa,const double *const *field_ptrs,
    std::size_t base,int nf,double x,double y,double z,int symmetry,double *out)
{
    const int ORDN=2*ghost_width;
    const int NO_SYMM=0,EQUATORIAL=1,OCTANT=2;(void)EQUATORIAL;
    const double *X=g.coordinates[0],*Y=g.coordinates[1],*Z=g.coordinates[2];
    const int *ex=g.shape,*pex=g.padded;
    double dX=X_at_1b(X,2)-X_at_1b(X,1);
    double dY=X_at_1b(Y,2)-X_at_1b(Y,1);
    double dZ=X_at_1b(Z,2)-X_at_1b(Z,1);
    constexpr double fshift=static_cast<double>(0.4f);
    int cxI[3];
    cxI[0]=(int)((x-X_at_1b(X,1))/dX+fshift)+1;
    cxI[1]=(int)((y-X_at_1b(Y,1))/dY+fshift)+1;
    cxI[2]=(int)((z-X_at_1b(Z,1))/dZ+fshift)+1;
    int cxB[3],cxT[3],cmin[3],cmax[3];
    for(int m=0;m<3;++m){cxB[m]=cxI[m]-ORDN/2+1;cxT[m]=cxB[m]+ORDN-1;cmin[m]=1;cmax[m]=ex[m];}
#ifdef Vertex
    const int sym_shift=-ORDN/2+2;const int ref_shift=2;
#else
    const int sym_shift=-ORDN/2+1;const int ref_shift=1;
#endif
    if(symmetry==OCTANT&&fabs(X_at_1b(X,1))<dX)cmin[0]=sym_shift;
    if(symmetry==OCTANT&&fabs(X_at_1b(Y,1))<dY)cmin[1]=sym_shift;
    if(symmetry!=NO_SYMM&&fabs(X_at_1b(Z,1))<dZ)cmin[2]=sym_shift;
    for(int m=0;m<3;++m){
        if(cxB[m]<cmin[m]){cxB[m]=cmin[m];cxT[m]=cxB[m]+ORDN-1;}
        if(cxT[m]>cmax[m]){cxT[m]=cmax[m];cxB[m]=cxT[m]+1-ORDN;}
    }
    double cx[3];
    cx[0]=(cxB[0]>0)?(x-X_at_1b(X,cxB[0]))/dX:(x+X_at_1b(X,ref_shift-cxB[0]))/dX;
    cx[1]=(cxB[1]>0)?(y-X_at_1b(Y,cxB[1]))/dY:(y+X_at_1b(Y,ref_shift-cxB[1]))/dY;
    cx[2]=(cxB[2]>0)?(z-X_at_1b(Z,cxB[2]))/dZ:(z+X_at_1b(Z,ref_shift-cxB[2]))/dZ;
    double wx[MAX_ORDN],wy[MAX_ORDN],wz[MAX_ORDN];
    compute_lagrange_weights(cx[0],ORDN,wx);
    compute_lagrange_weights(cx[1],ORDN,wy);
    compute_lagrange_weights(cx[2],ORDN,wz);
    for(int v=0;v<nf;++v){
        const double s3[3]={soa[3*v],soa[3*v+1],soa[3*v+2]};
        const double *f=field_ptrs[base+(size_t)v];
        double val=0.0;bool valid=true;
        for(int k=0;k<ORDN&&valid;++k){
            const int gz=cxB[2]+k;
            for(int j=0;j<ORDN&&valid;++j){
                const int gy=cxB[1]+j;
                for(int i=0;i<ORDN;++i){
                    const int gx=cxB[0]+i;
                    const double fv=fetch_f_sym_pad(f,ex,pex,gx,gy,gz,s3);
                    if(isnan(fv)){valid=false;break;}
                    val+=wx[i]*wy[j]*wz[k]*fv;
                }
            }
        }
        out[v]=valid?val:NAN;
    }
}
__global__ void sample_fused_kernel(int n,int nf,int my_rank,const GridBlock *blocks,
    const int *selection,const double *local,const double *const *field_ptrs,
    const double *soa,double *sample,int symmetry){
    int q=blockIdx.x*blockDim.x+threadIdx.x;if(q>=n)return;
    const int b=selection[q];if(b<0)return;
    const GridBlock g=blocks[b];if(g.rank!=my_rank)return;
    const double x=local[q],y=local[n+q],z=local[2*n+q];
    double out[max_fields];
    if(g.shell_type<0)
        fused_cartesian_point(g,soa,field_ptrs,(size_t)b*(size_t)nf,nf,x,y,z,symmetry,out);
    else
        for(int v=0;v<nf;++v){
            const double s3[3]={soa[3*v],soa[3*v+1],soa[3*v+2]};double value;
            ahf_neville_interp_ss_device(g.shape,g.padded,g.coordinates[0],g.coordinates[1],g.coordinates[2],
                field_ptrs[(size_t)b*nf+v],&value,x,y,z,2*ghost_width,s3,symmetry,g.shell_type);
            out[v]=value;
        }
    for(int v=0;v<nf;++v)sample[(size_t)q*nf+v]=out[v];
}
__global__ void finite_kernel(int count,const double *values,int *bad){
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<count&&!isfinite(values[i]))atomicExch(bad,1);
}

void append_block(std::vector<GridBlock> &metadata,std::vector<Block*> &addresses,
    Block &block,const double *patch_bbox,const double *spacing,const int *lli,const int *uui,int sst,int rank){
    GridBlock g{};g.rank=block.rank;g.shell_type=sst;
#ifdef Vertex
    const double inner=ghost_width-.5;
#else
    const double inner=ghost_width;
#endif
    for(int a=0;a<3;++a){
        double dh=spacing[a];g.spacing[a]=dh;
        g.shape[a]=block.shape[a];g.padded[a]=block.rank==rank?block.shape_pad[a]:block.shape[a];
        bool low=feq(block.bbox[a],patch_bbox[a],dh/2);
        bool high=feq(block.bbox[a+3],patch_bbox[a+3],dh/2);
        g.patch_lower[a]=sst<0?patch_bbox[a]+lli[a]*dh-dh/100:-DBL_MAX;
        g.patch_upper[a]=sst<0?patch_bbox[a+3]-uui[a]*dh+dh/100:DBL_MAX;
        g.lower[a]=block.bbox[a]+(low?(sst<0?lli[a]*dh:0):inner*dh);
        g.upper[a]=block.bbox[a+3]-(high?(sst<0?uui[a]*dh:0):inner*dh);
        // No remote device address enters metadata.
        g.coordinates[a]=block.rank==rank?block.d_X[a]:nullptr;
    }
    metadata.push_back(g);addresses.push_back(block.rank==rank?&block:nullptr);
}
}

void cleanup_sampling(){buffers.release();}

int sample_collective(Surface *surface,int owner,bssn_class &adm,MyList<var> *fields,int symmetry,bool gauge){
    int rank=0,np=0;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&np);
    // Reject inconsistent/invalid roots and missing owner storage before Bcast.
    int collective_owner=0,all_ok=0;
    MPI_Allreduce(&owner,&collective_owner,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    int local_ok=owner==collective_owner&&owner>=0&&owner<np&&adm.GH
        &&(rank==owner?surface&&surface->d.xyz&&(gauge?surface->d.gauge:surface->d.raw):!surface);
    MPI_Allreduce(&local_ok,&all_ok,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);if(!all_ok)return 0;
    int n=surface?surface->d.n:0;MPI_Bcast(&n,1,MPI_INT,owner,MPI_COMM_WORLD);
    int nf=0;std::vector<var*> variables;
    for(auto *v=fields;v;v=v->next){variables.push_back(v->data);++nf;}
    local_ok=n>0&&nf==(gauge?4:35)&&n<=INT_MAX/35;
    MPI_Allreduce(&local_ok,&all_ok,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);if(!all_ok)return 0;
    std::vector<GridBlock> metadata;std::vector<Block*> addresses;
    for(int lev=adm.GH->levels-1;lev>=0;--lev)for(auto *p=adm.GH->PatL[lev];p;p=p->next){
        double dh[3];for(int a=0;a<3;++a)dh[a]=p->data->getdX(a);
        for(auto *b=p->data->blb;b;b=b->next){
            append_block(metadata,addresses,*b->data,p->data->bbox,dh,p->data->lli,p->data->uui,-1,rank);
            if(b==p->data->ble)break;
        }
    }
    double rmin=0,rmax=0;RadialMap radial;
#ifdef WithShell
    if(adm.SH){
        f_shellcordpar(radial.A,radial.B,radial.r0,radial.eps);
        rmin=adm.SH->Rrange[0];rmax=adm.SH->Rrange[1];
        double dh[3];for(int a=0;a<3;++a)dh[a]=adm.SH->getdX(a);
        const int zero[3]={0,0,0};
        for(auto *p=adm.SH->PatL;p;p=p->next)for(auto *b=p->data->blb;b;b=b->next){
            append_block(metadata,addresses,*b->data,p->data->bbox,dh,zero,zero,p->data->sst,rank);
            if(b==p->data->ble)break;
        }
    }
#endif
    std::vector<double> coordinates_host,send,receive;
    int bad=0;
    try{
        if(np>1){coordinates_host.resize(3*n);send.resize(n*nf);receive.resize(n*nf);}
        buffers.reserve(n,metadata.size());
        check(cudaMemcpyAsync(buffers.blocks,metadata.data(),metadata.size()*sizeof(GridBlock),cudaMemcpyHostToDevice,buffers.stream),"AHF refreshed ownership metadata");
        transfer_stats().metadata_h2d+=metadata.size()*sizeof(GridBlock);
        if(np==1){
            // Nonblocking streams have no implicit dependency on the writer.
            check(cudaStreamSynchronize(surface->stream),"AHF coordinates producer completion");
            check(cudaMemcpyAsync(buffers.xyz,surface->d.xyz,3*n*sizeof(double),cudaMemcpyDeviceToDevice,buffers.stream),"AHF local coordinates copy");
        }else if(surface){check(cudaMemcpyAsync(coordinates_host.data(),surface->d.xyz,3*n*sizeof(double),cudaMemcpyDeviceToHost,surface->stream),"AHF MPI coordinates staging");check(cudaStreamSynchronize(surface->stream),"AHF coordinates producer completion");transfer_stats().mpi_d2h+=3*n*sizeof(double);}
    }catch(const std::exception &e){std::fprintf(stderr,"AHF rank %d sampling preparation: %s\n",rank,e.what());bad=1;}
    MPI_Allreduce(&bad,&all_ok,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);if(all_ok)return 0;
    // Portable staging is deliberate for MPI communication. CUDA-aware buffers
    // must not be selected merely from an unverified environment variable.
    if(np>1)MPI_Bcast(coordinates_host.data(),3*n,MPI_DOUBLE,owner,MPI_COMM_WORLD);
    try{
        if(np>1){check(cudaMemcpyAsync(buffers.xyz,coordinates_host.data(),3*n*sizeof(double),cudaMemcpyHostToDevice,buffers.stream),"AHF MPI coordinates upload");transfer_stats().mpi_h2d+=3*n*sizeof(double);}
        check(cudaMemsetAsync(buffers.sample,0,n*nf*sizeof(double),buffers.stream),"AHF sample clear");
        check(cudaMemsetAsync(buffers.bad,0,sizeof(int),buffers.stream),"AHF sampling status clear");
        select_kernel<<<(n+127)/128,128,0,buffers.stream>>>(n,metadata.size(),buffers.blocks,buffers.xyz,buffers.local_xyz,buffers.selection,buffers.bad,rmin,rmax,radial);
        if(nf>max_fields)throw std::runtime_error("AHF sampling field count exceeds fused table capacity");
        std::vector<double*> ptr_table(metadata.size()*(size_t)nf,nullptr);
        std::vector<double> soa_table((size_t)3*nf);
        for(int v=0;v<nf;++v){soa_table[3*(size_t)v]=variables[v]->SoA[0];soa_table[3*(size_t)v+1]=variables[v]->SoA[1];soa_table[3*(size_t)v+2]=variables[v]->SoA[2];}
        for(std::size_t b=0;b<metadata.size();++b){
            if(metadata[b].rank!=rank)continue;
            Block &block=*addresses[b];
            // 主机权威字段在此自动同步到设备（require_on_gpu 不会覆盖更新的 GPU 数据）
            bool uploaded=false;
            for(var *v:variables)if(!block.gpu_valid[v->sgfn]){block.require_on_gpu(v->sgfn);uploaded=true;
                if(!block.gpu_valid[v->sgfn])throw std::runtime_error("AHF sampled field is neither cpu- nor gpu-authoritative");}
            if(uploaded)GPUManager::getInstance().fence_pool_after_memory();
            for(int v=0;v<nf;++v)ptr_table[b*(size_t)nf+v]=block.d_fgfs[variables[v]->sgfn];
            // 采样在另一条流上执行，需先等待该块流上的先前写入（设备侧事件）
            check(cudaEventRecord(buffers.events[b],block.stream),"AHF sampled block writer event");
            check(cudaStreamWaitEvent(buffers.stream,buffers.events[b],0),"AHF sampled block writer wait");
        }
        check(cudaMemcpyAsync(buffers.field_ptrs,ptr_table.data(),ptr_table.size()*sizeof(double*),cudaMemcpyHostToDevice,buffers.stream),"AHF field pointer table");
        check(cudaMemcpyAsync(buffers.soa,soa_table.data(),soa_table.size()*sizeof(double),cudaMemcpyHostToDevice,buffers.stream),"AHF variable parity table");
        sample_fused_kernel<<<(n+127)/128,128,0,buffers.stream>>>(n,nf,rank,buffers.blocks,buffers.selection,buffers.local_xyz,buffers.field_ptrs,buffers.soa,buffers.sample,symmetry);
        finite_kernel<<<(n*nf+127)/128,128,0,buffers.stream>>>(n*nf,buffers.sample,buffers.bad);
        check(cudaGetLastError(),"AHF sample kernels");
        check(cudaMemcpyAsync(&bad,buffers.bad,sizeof(int),cudaMemcpyDeviceToHost,buffers.stream),"AHF sampling status scalar");
        if(np>1){check(cudaMemcpyAsync(send.data(),buffers.sample,n*nf*sizeof(double),cudaMemcpyDeviceToHost,buffers.stream),"AHF MPI samples staging");transfer_stats().mpi_d2h+=n*nf*sizeof(double);}
        check(cudaStreamSynchronize(buffers.stream),"AHF sampling completion");transfer_stats().scalar_d2h+=sizeof(int);
    }catch(const std::exception &e){std::fprintf(stderr,"AHF rank %d sampling: %s\n",rank,e.what());bad=1;}
    MPI_Allreduce(&bad,&all_ok,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);if(all_ok)return 0;
    if(np>1)MPI_Reduce(send.data(),receive.data(),n*nf,MPI_DOUBLE,MPI_SUM,owner,MPI_COMM_WORLD);
    // Only one rank contributes each point, so reduction is collection, not
    // numerical averaging. Host arrays are communication staging only.
    try{
        if(np==1){
            // Sampling has completed above. Queue on the consumer stream and
            // finish the copy before a later call can reuse sampling scratch.
            check(cudaMemcpyAsync(gauge?surface->d.gauge:surface->d.raw,buffers.sample,n*nf*sizeof(double),cudaMemcpyDeviceToDevice,surface->stream),"AHF local samples copy");
            check(cudaStreamSynchronize(surface->stream),"AHF local samples completion");
        }else if(surface){check(cudaMemcpyAsync(gauge?surface->d.gauge:surface->d.raw,receive.data(),n*nf*sizeof(double),cudaMemcpyHostToDevice,surface->stream),"AHF collected samples upload");check(cudaStreamSynchronize(surface->stream),"AHF staging lifetime completion");transfer_stats().mpi_h2d+=n*nf*sizeof(double);}
    }catch(const std::exception &e){std::fprintf(stderr,"AHF rank %d collected samples: %s\n",rank,e.what());bad=1;}
    MPI_Allreduce(&bad,&all_ok,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);return !all_ok;
}
}}
#endif
