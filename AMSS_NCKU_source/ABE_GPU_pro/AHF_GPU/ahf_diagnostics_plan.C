#include "macrodef.h"
#include "AHF/driver.h"
#include "ahf_diagnostics_gpu.h"
#include <stdexcept>
#include <climits>
namespace ahf_diagnostics {
#include "ahf_static_angles.h"
using namespace AHFinderDirect;
HostPlan make_plan(patch_system &ps,int symmetry){
 if(symmetry<0||symmetry>2||fd_grid::molecule_radius()!=2)throw std::invalid_argument("Expected symmetry 0/1/2 and AHF ORDER4 radius2");
 const auto expected=symmetry==0?patch_system::patch_system__full_sphere:symmetry==1?patch_system::patch_system__plus_z_hemisphere:patch_system::patch_system__plus_xyz_octant_mirrored;
 if(ps.type()!=expected)throw std::invalid_argument("Symmetry does not match patch system");
 HostPlan plan;plan.n=ps.N_grid_points();plan.ng=ps.ghosted_N_grid_points();plan.symmetry=symmetry;
 if(plan.n<=0||plan.ng<=0)throw std::invalid_argument("Empty topology");
 for(int pn=0;pn<ps.N_patches();++pn){patch &p=ps.ith_patch(pn);for(int r=p.min_irho();r<=p.max_irho();++r)for(int t=p.min_isigma();t<=p.max_isigma();++t)p.ghosted_gridfn(gfns::gfn__h,r,t)=1.;}
 ps.synchronize();ps.compute_synchronize_Jacobian();
 // Accepted ordered ghost schedule: symmetry, interpatch, symmetry corners.
 for(int phase=0;phase<3;++phase)for(int pn=0;pn<ps.N_patches();++pn)for(int mn=0;mn<=1;++mn)for(int rho=0;rho<=1;++rho){
  const ghost_zone &gz=ps.ith_patch(pn).minmax_ang_ghost_zone(mn,rho);
  if((phase==1)!=gz.is_interpatch())continue;
  const patch_edge &e=gz.my_edge(),&src=gz.other_edge();
  for(int perp=gz.min_iperp();perp<=gz.max_iperp();++perp)for(int par=gz.min_ipar(perp);par<=gz.max_ipar(perp);++par){
   if(phase!=1&&!e.ipar_is_in_selected_part(phase==2,phase==0,par))continue;
   GhostOperation op{ps.ghosted_gpn_of_patch_irho_isigma(e.my_patch(),e.irho_of_iperp_ipar(perp,par),e.isigma_of_iperp_ipar(perp,par)),static_cast<int>(plan.sources.size()),0};
   if(gz.is_symmetry()){
    const auto &sg=gz.cast_to_symmetry_ghost_zone();int pp=sg.iperp_map_of_iperp(perp),pa=sg.ipar_map_of_ipar(par);
    plan.sources.push_back(ps.ghosted_gpn_of_patch_irho_isigma(src.my_patch(),src.irho_of_iperp_ipar(pp,pa),src.isigma_of_iperp_ipar(pp,pa)));plan.weights.push_back(1.);
   }else{
    int pp=gz.Jacobian_y_iperp(perp),pos=gz.Jacobian_y_ipar_posn(perp,par);
    for(int m=gz.Jacobian_min_y_ipar_m();m<=gz.Jacobian_max_y_ipar_m();++m){plan.sources.push_back(ps.ghosted_gpn_of_patch_irho_isigma(src.my_patch(),src.irho_of_iperp_ipar(pp,pos+m),src.isigma_of_iperp_ipar(pp,pos+m)));plan.weights.push_back(gz.Jacobian(perp,par,m));}
   }
   op.end=static_cast<int>(plan.sources.size());plan.operations.push_back(op);
  }
 }
 const double unused_origin[3]={0,0,0};
 for(int pn=0;pn<ps.N_patches();++pn){
  const patch &p=ps.ith_patch(pn);Patch range{static_cast<int>(plan.points.size()),0,p.delta_rho(),p.delta_sigma()};
  for(int r=p.min_irho();r<=p.max_irho();++r)for(int t=p.min_isigma();t<=p.max_isigma();++t){
   if(ps.gpn_of_patch_irho_isigma(p,r,t)!=static_cast<int>(plan.points.size()))throw std::runtime_error("Original point order mismatch");
   Point m;int family=p.ctype()=='z'?0:p.ctype()=='x'?1:2;
   if(!make_static_metadata(family,p.rho_of_irho(r),p.sigma_of_isigma(t),unused_origin,m))throw std::runtime_error("Illegal static angles");
   m.ghost=ps.ghosted_gpn_of_patch_irho_isigma(p,r,t);m.dr=p.inverse_delta_rho();m.ds=p.inverse_delta_sigma();
   m.cr=patch::integration_coeff(patch::integration_method__automatic_choice,p.max_irho()-p.min_irho(),r-p.min_irho());
   m.cs=patch::integration_coeff(patch::integration_method__automatic_choice,p.max_isigma()-p.min_isigma(),t-p.min_isigma());
   plan.points.push_back(m);
   for(int mr=-2;mr<=2;++mr)for(int mt=-2;mt<=2;++mt)plan.stencils.push_back(ps.ghosted_gpn_of_patch_irho_isigma(p,r+mr,t+mt));
  }
  range.end=static_cast<int>(plan.points.size());plan.patches.push_back(range);
  const char *planes[3]={"xy","xz","yz"};
  for(int a=0;a<3;++a){
   if(p.ctype()!=planes[a][0]&&p.ctype()!=planes[a][1])continue;
   bool rho=(a==1&&p.ctype()=='x')||a==2;int fixed=0;bool have=false;
   for(double ang:{0.,90.,180.,-90.})if(rho?p.is_valid_dsigma(ang):p.is_valid_drho(ang)){fixed=rho?p.isigma_of_sigma(p.sigma_of_dsigma(ang)):p.irho_of_rho(p.rho_of_drho(ang));have=true;break;}
   if(!have)throw std::runtime_error("No circumference plane");
   Arc arc{a,rho?0:1,static_cast<int>(plan.arc_points.size()),0,rho?p.delta_rho():p.delta_sigma()};
   int lo=rho?p.min_irho():p.min_isigma(),hi=rho?p.max_irho():p.max_isigma();
   for(int i=lo;i<=hi;++i){plan.arc_points.push_back(ps.gpn_of_patch_irho_isigma(p,rho?i:fixed,rho?fixed:i));plan.arc_coefficients.push_back(patch::integration_coeff(patch::integration_method__automatic_choice,hi-lo,i-lo));}
   arc.end=static_cast<int>(plan.arc_points.size());plan.arcs.push_back(arc);
  }
 }
 return plan;
}
}
