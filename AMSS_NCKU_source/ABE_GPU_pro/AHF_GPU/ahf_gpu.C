#include "macrodef.h"
#if defined(USE_GPU) && defined(With_AHF)
#include "ahf_gpu.h"
#include "ahf_ilucg_gpu.h"
#include "../Parallel.h"
#include "../gpu_manager.h"
#include <mpi.h>
#include "../AHF/driver.h"
#include "../AHF/myglobal.h"
#include <vector>
#include "ahf_static_metadata.h"
#include <cmath>
#include <stdexcept>
#include <string>
#include <memory>
#include <cstdio>

namespace AHFinderDirect { extern struct state state; namespace gpu {
namespace {
TransferStats stats;
std::vector<std::unique_ptr<Surface> > surfaces;
template<class T> void allocate(T *&p,std::size_t n){check(cudaMalloc(reinterpret_cast<void **>(&p),sizeof(T)*n),"AHF allocation");}
template<class T> void upload(T *&p,const std::vector<T> &v){allocate(p,v.size());check(cudaMemcpy(p,v.data(),v.size()*sizeof(T),cudaMemcpyHostToDevice),"AHF static metadata");stats.metadata_h2d+=v.size()*sizeof(T);}
void copy_output(double *host,const double *device,std::size_t n,cudaStream_t stream){check(cudaMemcpyAsync(host,device,n*sizeof(double),cudaMemcpyDeviceToHost,stream),"AHF output array");stats.output_d2h+=n*sizeof(double);}
void release_surface(DeviceSurface &d,cudaStream_t stream) noexcept {
    if(stream)cudaStreamSynchronize(stream);
    ahf_diagnostics::release_plan(d.diagnostic_plan);
    cudaFree(d.metadata);cudaFree(d.runtime_geometry);cudaFree(d.initial_direction);
    cudaFree(d.nominal_ghost);cudaFree(d.newton);
    cudaFree(d.instruction_rows);cudaFree(d.instructions);
    cudaFree(d.angular.rows);cudaFree(d.angular.columns);
    double *ptrs[]={d.h,d.oldh,d.ghost_h,d.xyz,d.raw,d.geometry,d.gauge,d.theta,d.coefficients,d.values,d.delta,d.saved_theta,d.scratch6,d.scalars,d.diagnostics,d.signature,d.origin};
    for(double *p:ptrs)cudaFree(p);
    cudaFree(d.status);
    if(stream)cudaStreamDestroy(stream);
    d=DeviceSurface{};
}
}
TransferStats &transfer_stats(){return stats;}
void check(cudaError_t e,const char *op){if(e!=cudaSuccess)throw std::runtime_error(std::string(op)+": "+cudaGetErrorString(e));}
Surface::Surface(patch_system &ps,int symmetry):reference(&ps){
    try {
        check(cudaStreamCreateWithFlags(&stream,cudaStreamNonBlocking),"AHF stream");
        d.n=ps.N_grid_points();d.ng=ps.ghosted_N_grid_points();d.symmetry=symmetry;
        // This setup seeds a constant scratch h solely to initialize static ghost tables.
        const auto topology=ahf_diagnostics::make_plan(ps,symmetry);
        const auto ordered=make_ordered_plan(ps);
        d.diagnostic_plan=ahf_diagnostics::upload_plan(topology);
        d.nnz=ordered.columns.size();
        upload(d.angular.rows,ordered.rows);upload(d.angular.columns,ordered.columns);
        upload(d.instruction_rows,ordered.instruction_rows);upload(d.instructions,ordered.instructions);
        std::vector<ahf_native::StaticMetadata> metadata(d.n);
        std::vector<int> nominal(d.n);std::vector<double> directions(3*d.n);
        const double origin[3]={ps.origin_x(),ps.origin_y(),ps.origin_z()};
        for(int q=0;q<d.n;++q){
            int r,t;const patch &p=ps.patch_irho_isigma_of_gpn(q,r,t);
            const double rho=p.rho_of_irho(r),sigma=p.sigma_of_isigma(t);
            int family=p.ctype()=='z'?0:p.ctype()=='x'?1:2;
            if(!ahf_native::make_static_metadata(family,rho,sigma,origin,metadata[q]))
                throw std::runtime_error("AHF unsupported static angles");
            nominal[q]=topology.points[q].ghost;
            // Static direction cosines are used only by original ellipsoid initialization.
            p.xyzcos_of_rho_sigma(rho,sigma,directions[3*q],directions[3*q+1],directions[3*q+2]);
        }
        upload(d.metadata,metadata);upload(d.nominal_ghost,nominal);upload(d.initial_direction,directions);
        allocate(d.runtime_geometry,d.n);allocate(d.newton,1);
        allocate(d.h,d.n);allocate(d.oldh,d.n);allocate(d.ghost_h,d.ng);allocate(d.xyz,3*d.n);
        allocate(d.raw,35*d.n);allocate(d.geometry,35*d.n);allocate(d.gauge,4*d.n);
        allocate(d.theta,d.n);allocate(d.coefficients,5*d.n);allocate(d.values,d.nnz);
        allocate(d.delta,d.n);allocate(d.saved_theta,d.n);allocate(d.scratch6,6*d.n);
        allocate(d.scalars,3);allocate(d.diagnostics,30);allocate(d.signature,d.n);
        allocate(d.origin,3);allocate(d.status,1);
        linear_solver=gpu_ahf_linear_create(stream,d.n,d.nnz,d.angular.rows,d.angular.columns);
        check(cudaMemsetAsync(d.ghost_h,0,d.ng*sizeof(double),stream),"AHF initialize unused ghosts");
        check(cudaMemsetAsync(d.status,0,sizeof(int),stream),"AHF initialize status");
        check(cudaMemcpyAsync(d.origin,origin,sizeof(origin),cudaMemcpyHostToDevice,stream),"AHF origin metadata");
        check(cudaStreamSynchronize(stream),"AHF metadata lifetime");stats.metadata_h2d+=sizeof(origin);
    } catch (...) {
        gpu_ahf_linear_destroy(linear_solver);linear_solver=nullptr;
        release_surface(d,stream);stream=nullptr;throw;
    }
}
Surface::~Surface(){gpu_ahf_linear_destroy(linear_solver);release_surface(d,stream);}
void download_diagnostics(Surface &s,BH_diagnostics &v){
    double a[30];copy_output(a,s.d.diagnostics,30,s.stream);check(cudaStreamSynchronize(s.stream),"AHF diagnostics output");
    if(std::isnan(a[8]))std::fprintf(stderr,s.exact_pole_signature_nans?
        "AHF reference defect: scratch quadrupole_zz is NaN after exact-pole signature; not a finite-parity pass\n":
        "AHF unclassified NaN in scratch quadrupole_zz; not a finite-parity pass\n");
    for(int j=0;j<30;++j)if(!std::isfinite(a[j])&&!(j==8&&std::isnan(a[j])&&s.exact_pole_signature_nans))
        std::fprintf(stderr,"AHF unclassified nonfinite diagnostic slot %d; not a finite-parity pass\n",j);
    v.copy_from_buffer(a);v.Px=a[24];v.Py=a[25];v.Pz=a[26];v.Sx=a[27];v.Sy=a[28];v.Sz=a[29];
}
void download_output(Surface &s,bool sig){
    synchronize(s);copy_output(s.reference->ghosted_gridfn_data(gfns::gfn__h),s.d.ghost_h,s.d.ng,s.stream);
    coordinates(s);
    for(int a=0;a<3;++a)copy_output(s.reference->gridfn_data(gfns::gfn__global_x+a),s.d.xyz+a*s.d.n,s.d.n,s.stream);
    if(sig)copy_output(s.reference->gridfn_data(gfns::gfn__global_zz),s.d.signature,s.d.n,s.stream);
    check(cudaStreamSynchronize(s.stream),"AHF surface output completion");
}
void setup(){
    if(state.N_horizons<0||state.N_active_procs<0||state.N_active_procs>state.N_procs||
       (state.N_horizons>0&&state.N_active_procs==0))
        throw std::runtime_error("AHF invalid horizon ownership metadata");
    surfaces.clear();surfaces.resize(state.N_horizons+1);
    for(int hn=1;hn<=state.N_horizons;++hn)if((hn-1)%state.N_active_procs==state.my_proc){
        if(!state.AH_data_array[hn]||!state.AH_data_array[hn]->ps_ptr)
            throw std::runtime_error("AHF missing owner angular patch system");
        surfaces[hn].reset(new Surface(*state.AH_data_array[hn]->ps_ptr,state.Symmetry));
    }
}
void cleanup(){surfaces.clear();cleanup_sampling();}
void prepare_derivatives(bssn_class &adm,MyList<var> *output_list){
    GPUManager::getInstance().synchronize_all();
    var *input[7]={adm.phi0,adm.gxx0,adm.gxy0,adm.gxz0,adm.gyy0,adm.gyz0,adm.gzz0};
    var *output[21]={adm.dtSfx_rhs,adm.dtSfy_rhs,adm.dtSfz_rhs,adm.Gamxxx,adm.Gamyxx,adm.Gamzxx,
        adm.Gamxxy,adm.Gamyxy,adm.Gamzxy,adm.Gamxxz,adm.Gamyxz,adm.Gamzxz,adm.Gamxyy,adm.Gamyyy,adm.Gamzyy,
        adm.Gamxyz,adm.Gamyyz,adm.Gamzyz,adm.Gamxzz,adm.Gamyzz,adm.Gamzzz};
    auto launch=[&](Block &block,int sst){
        if(block.rank!=adm.myrank)return;
        // 导数内核与先前写入同在块流上，流内天然有序，无需逐块同步
        const double *fields[7],*coords[3],*transform[9];double *out[21];
        bool uploaded=false;
        for(int a=0;a<3;++a)coords[a]=block.d_X[a];
        for(int v=0;v<7;++v){
            if(!block.gpu_valid[input[v]->sgfn]){block.require_on_gpu(input[v]->sgfn);uploaded=true;
                if(!block.gpu_valid[input[v]->sgfn])throw std::runtime_error("AHF derivative input is neither cpu- nor gpu-authoritative");}
            fields[v]=block.d_fgfs[input[v]->sgfn];
        }
        if(sst>=0)for(int a=0;a<9;++a){int index=adm.SH->fngfs+ShellPatch::drhodx+a;
            if(!block.gpu_valid[index]){block.require_on_gpu(index);uploaded=true;
                if(!block.gpu_valid[index])throw std::runtime_error("AHF shell coordinate Jacobian is neither cpu- nor gpu-authoritative");}
            transform[a]=block.d_fgfs[index];}
        if(uploaded)GPUManager::getInstance().fence_pool_after_memory();
        for(int v=0;v<21;++v)out[v]=block.d_fgfs[output[v]->sgfn];
        derivatives(block.stream,block.shape,block.shape_pad,coords,fields,out,adm.Symmetry,sst,sst>=0?transform:nullptr);
        for(int v=0;v<21;++v)block.mark_gpu_modified(output[v]->sgfn);
    };
    for(int lev=0;lev<adm.GH->levels;++lev){
        for(auto *p=adm.GH->PatL[lev];p;p=p->next)for(auto *b=p->data->blb;b;b=b->next){launch(*b->data,-1);if(b==p->data->ble)break;}
        GPUManager::getInstance().synchronize_all();
        Parallel::Sync_GPU(adm.GH->PatL[lev],output_list,adm.Symmetry);
    }
#ifdef WithShell
    if(adm.SH){
        for(auto *p=adm.SH->PatL;p;p=p->next)for(auto *b=p->data->blb;b;b=b->next){launch(*b->data,p->data->sst);if(b==p->data->ble)break;}
        GPUManager::getInstance().synchronize_all();adm.SH->Synch_GPU(output_list,adm.Symmetry);
#ifndef Vertex
        // Cell sampling needs derivatives in Cartesian owner cells at the Shell seam.
        GPUManager::getInstance().synchronize_all();adm.SH->CS_Inter_GPU(output_list,adm.Symmetry);
#endif
        // Propagate Cell seam values to MPI ghosts; preserve the validated
        // post-Shell Cartesian refresh in both grid configurations.
        for(int lev=0;lev<adm.GH->levels;++lev)
            Parallel::Sync_GPU(adm.GH->PatL[lev],output_list,adm.Symmetry);
    }
#endif
    GPUManager::getInstance().synchronize_all();
}
void find_horizons(int count,int *dumpid,double *xc,double *yc,double *zc,double *xr,double *yr,double *zr,bool *trigger,double *dt){
    // Serial horizon scheduling deliberately keeps every MPI rank in precisely the
    // same collectives, including ranks with no horizons. Ownership remains round-robin.
    for(int hn=1;hn<=count;++hn){
        int owner=(hn-1)%state.N_active_procs;Surface *s=state.my_proc==owner?surfaces[hn].get():nullptr;
        AH_data *a=s?state.AH_data_array[hn]:nullptr;
        int enabled=s?(trigger[hn-1]&&!a->stop_finding):0;MPI_Bcast(&enabled,1,MPI_INT,owner,MPI_COMM_WORLD);if(!enabled)continue;
        if(s){
            if(a->found_flag)a->initial_find_flag=false;
            else if(!a->recentering_flag){double c[3]={xc[hn-1],yc[hn-1],zc[hn-1]},r[3]={xr[hn-1],yr[hn-1],zr[hn-1]};
                // Match Newton.C's initial recenter test using the supplied
                // bounding-box metadata, not an unconditional center change.
                double origin[3]={s->reference->origin_x(),s->reference->origin_y(),s->reference->origin_z()};
                bool centered=true;for(int axis=0;axis<3;++axis){
                    const double hi=c[axis]+r[axis]/2,lo=c[axis]-r[axis]/2;
                    centered=centered&&std::abs(hi+lo-2.0*origin[axis])<.45*(hi-lo);
                }
                if(!centered){check(cudaMemcpyAsync(s->d.origin,c,sizeof(c),cudaMemcpyHostToDevice,s->stream),"AHF initial origin");
                    check(cudaStreamSynchronize(s->stream),"AHF initial origin lifetime");stats.metadata_h2d+=sizeof(c);s->reference->recentering(c[0],c[1],c[2]);}
                initialize(*s,c,r);a->initial_find_flag=true;}
            save_old_shape(*s);
        }
        int maxit=s?(a->initial_find_flag?80:20):0;MPI_Bcast(&maxit,1,MPI_INT,owner,MPI_COMM_WORLD);
        int found=0;
        for(int it=1;it<=maxit;++it){
            if(s)coordinates(*s);
            int sampled=sample_collective(s,owner,*state.ADM,state.AHList,state.Symmetry);
            double n[3]={0,0,1};if(s&&sampled){transform_geometry(s->stream,s->d.n,s->d.raw,s->d.geometry);evaluate(*s,-a->surface_expansion);norms(*s,n);}
            int action=0;if(s&&sampled&&!n[2]&&n[0]<=1e10&&n[1]<=1e10)action=n[0]<=1e-11?2:(it<maxit?1:0);
            MPI_Bcast(&action,1,MPI_INT,owner,MPI_COMM_WORLD);
            // if(action==2){found=1;if(s)std::printf(" found horizon #%d with %d steps (device Theta infinity %.17g)\n",hn,it,n[0]);break;}
            if(action==2){found=1;if(s)std::printf(" found horizon #%d with %d steps \n",hn,it,n[0]);break;}
            if(action==0)break;
            if(s){assemble(*s);perturb(*s);coordinates(*s);}
            int radial=sample_collective(s,owner,*state.ADM,state.AHList,state.Symmetry);
            int success=radial;if(s&&radial){transform_geometry(s->stream,s->d.n,s->d.raw,s->d.geometry);evaluate(*s,-a->surface_expansion,false);norms(*s,n);success=!n[2];}
            MPI_Bcast(&success,1,MPI_INT,owner,MPI_COMM_WORLD);
            if(s)restore_and_radial(*s,success);
            if(!success)break;
            int linear=s?linear_step(*s):0;MPI_Bcast(&linear,1,MPI_INT,owner,MPI_COMM_WORLD);if(linear<=0)break;
            if(s)update(*s);
        }
        MPI_Bcast(&found,1,MPI_INT,owner,MPI_COMM_WORLD);if(a)a->found_flag=found;
        if(found){
            int want_signature=dumpid[hn-1]>0&&dt[hn-1]>0;
            if(want_signature&&!sample_collective(s,owner,*state.ADM,state.GaugeList,state.Symmetry,true))throw std::runtime_error("AHF gauge interpolation failed after convergence");
            if(s){
                if(dumpid[hn-1]>0){if(want_signature)signature(*s,dt[hn-1]);download_output(*s,want_signature);char name[100];std::snprintf(name,sizeof(name),"ah%02d_%05d.dat",hn,dumpid[hn-1]);
                    if(want_signature)s->reference->print_gridfn_with_xyz(gfns::gfn__global_zz,true,gfns::gfn__h,name);
                    else s->reference->print_ghosted_gridfn_with_xyz(gfns::gfn__h,true,gfns::gfn__h,name,false);}
                diagnostics(*s);download_diagnostics(*s,a->BH_diagnostics);
                if(!a->BH_diagnostics_fileptr){
                    char name[50];std::snprintf(name,sizeof(name),"infoah%02d.dat",hn);
                    a->BH_diagnostics_fileptr=std::fopen(name,"w");
                    if(!a->BH_diagnostics_fileptr)throw std::runtime_error("AHF diagnostics output file could not be opened");
                    std::fprintf(a->BH_diagnostics_fileptr,"#time Mass x y z Px Py Pz Sx Sy Sz\n");
                }
                a->BH_diagnostics.output(a->BH_diagnostics_fileptr,*state.PhysTime);recenter(*s);a->recentering_flag=true;
            }
        }else if(state.my_proc==owner)std::fprintf(stderr,"AHF: device solve failed for horizon #%d (no CPU fallback)\n",hn);
    }
}
}}
#endif
