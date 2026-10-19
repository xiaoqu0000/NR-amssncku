

#include "macrodef.fh"

  function compute_rhs_z4cnot(ex, T,X, Y, Z,                                      &
               chi    ,   trK    ,                                             &
               dxx    ,   gxy    ,   gxz    ,   dyy    ,   gyz    ,   dzz,     &
               Axx    ,   Axy    ,   Axz    ,   Ayy    ,   Ayz    ,   Azz,     &
               Gamx   ,  Gamy    ,  Gamz    ,                                  &
               Lap    ,  betax   ,  betay   ,  betaz   ,                       &
               dtSfx  ,  dtSfy   ,  dtSfz   ,                                  &
               TZ     ,                                                        &
               chi_rhs,   trK_rhs,                                             &
               gxx_rhs,   gxy_rhs,   gxz_rhs,   gyy_rhs,   gyz_rhs,   gzz_rhs, &
               Axx_rhs,   Axy_rhs,   Axz_rhs,   Ayy_rhs,   Ayz_rhs,   Azz_rhs, &
               Gamx_rhs,  Gamy_rhs,  Gamz_rhs,                                 &
               Lap_rhs,  betax_rhs,  betay_rhs,  betaz_rhs,                    &
               dtSfx_rhs,  dtSfy_rhs,  dtSfz_rhs,                              &
               TZ_rhs   ,                                                      &
               rho,Sx,Sy,Sz,Sxx,Sxy,Sxz,Syy,Syz,Szz,                           &
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,                      &
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,                      &
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,                      &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,                                        &
               Hcon,Mxcon,Mycon,Mzcon,Gmxcon,Gmycon,Gmzcon,                    &
! co is not used here, we always compute constraint               
               Symmetry,Lev,eps,co,chitiny)  result(gont)
  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3), Symmetry,Lev,co
  real*8, intent(in ):: T,chitiny
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: chi,dxx,dyy,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: trK
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: gxy,gxz,gyz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamx,Gamy,Gamz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Lap, betax, betay, betaz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: dtSfx,  dtSfy,  dtSfz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: TZ
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: chi_rhs,trK_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gxx_rhs,gxy_rhs,gxz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gyy_rhs,gyz_rhs,gzz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Axx_rhs,Axy_rhs,Axz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Ayy_rhs,Ayz_rhs,Azz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamx_rhs,Gamy_rhs,Gamz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Lap_rhs, betax_rhs, betay_rhs, betaz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: dtSfx_rhs,dtSfy_rhs,dtSfz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: TZ_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: rho,Sx,Sy,Sz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Sxx,Sxy,Sxz,Syy,Syz,Szz
! when out, physical second kind of connection  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamzyy, Gamzyz, Gamzzz
! when out, physical Ricci tensor  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Rxx,Rxy,Rxz,Ryy,Ryz,Rzz
! when out, constraint violation  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Hcon,Mxcon,Mycon,Mzcon,Gmxcon,Gmycon,Gmzcon
  real*8,intent(in) :: eps
!  gont = 0: success; gont = 1: something wrong
  integer::gont,compute_rhs_z4c

  real*8, dimension(ex(1),ex(2),ex(3)) :: chihere

  chihere = chi
  call lowerboundset(ex,chihere,chitiny)

  gont = compute_rhs_z4c(ex, T,X, Y, Z,                                      &
               chihere,   trK    ,                                             &
               dxx    ,   gxy    ,   gxz    ,   dyy    ,   gyz    ,   dzz,     &
               Axx    ,   Axy    ,   Axz    ,   Ayy    ,   Ayz    ,   Azz,     &
               Gamx   ,  Gamy    ,  Gamz    ,                                  &
               Lap    ,  betax   ,  betay   ,  betaz   ,                       &
               dtSfx  ,  dtSfy   ,  dtSfz   ,                                  &
               TZ     ,                                                        &
               chi_rhs,   trK_rhs,                                             &
               gxx_rhs,   gxy_rhs,   gxz_rhs,   gyy_rhs,   gyz_rhs,   gzz_rhs, &
               Axx_rhs,   Axy_rhs,   Axz_rhs,   Ayy_rhs,   Ayz_rhs,   Azz_rhs, &
               Gamx_rhs,  Gamy_rhs,  Gamz_rhs,                                 &
               Lap_rhs,  betax_rhs,  betay_rhs,  betaz_rhs,                    &
               dtSfx_rhs,  dtSfy_rhs,  dtSfz_rhs,                              &
               TZ_rhs   ,                                                      &
               rho,Sx,Sy,Sz,Sxx,Sxy,Sxz,Syy,Syz,Szz,                           &
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,                      &
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,                      &
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,                      &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,                                        &
               Hcon,Mxcon,Mycon,Mzcon,Gmxcon,Gmycon,Gmzcon,                    &            
               Symmetry,Lev,eps,co)

#if (ABV == 0)  
  call ricci_gamma(ex, X, Y, Z,                                      &
               chi,                                                  &
               dxx    ,   gxy    ,   gxz    ,   dyy    ,   gyz    ,   dzz,&
               Gamx   ,  Gamy    ,  Gamz    , &
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,&
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,&
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,&
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,&
               Symmetry)
#endif
  call constraint_bssn(ex, X, Y, Z,&
               chi,trK, &
               dxx,gxy,gxz,dyy,gyz,dzz, &
               Axx,Axy,Axz,Ayy,Ayz,Azz, &
               Gamx,Gamy,Gamz,&
               Lap,betax,betay,betaz,rho,Sx,Sy,Sz,&
               Gamxxx, Gamxxy, Gamxxz,Gamxyy, Gamxyz, Gamxzz, &
               Gamyxx, Gamyxy, Gamyxz,Gamyyy, Gamyyz, Gamyzz, &
               Gamzxx, Gamzxy, Gamzxz,Gamzyy, Gamzyz, Gamzzz, &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz, &
               Hcon,Mxcon,Mycon,Mzcon,Gmxcon,Gmycon,Gmzcon, &
               Symmetry)

  return

  end function compute_rhs_Z4cnot

#if 1
  function compute_rhs_z4c(ex, T,X, Y, Z,                                      &
               chi    ,   trK    ,                                             &
               dxx    ,   gxy    ,   gxz    ,   dyy    ,   gyz    ,   dzz,     &
               Axx    ,   Axy    ,   Axz    ,   Ayy    ,   Ayz    ,   Azz,     &
               Gamx   ,  Gamy    ,  Gamz    ,                                  &
               Lap    ,  betax   ,  betay   ,  betaz   ,                       &
               dtSfx  ,  dtSfy   ,  dtSfz   ,                                  &
               TZ     ,                                                        &
               chi_rhs,   trK_rhs,                                             &
               gxx_rhs,   gxy_rhs,   gxz_rhs,   gyy_rhs,   gyz_rhs,   gzz_rhs, &
               Axx_rhs,   Axy_rhs,   Axz_rhs,   Ayy_rhs,   Ayz_rhs,   Azz_rhs, &
               Gamx_rhs,  Gamy_rhs,  Gamz_rhs,                                 &
               Lap_rhs,  betax_rhs,  betay_rhs,  betaz_rhs,                    &
               dtSfx_rhs,  dtSfy_rhs,  dtSfz_rhs,                              &
               TZ_rhs   ,                                                      &
               rho,Sx,Sy,Sz,Sxx,Sxy,Sxz,Syy,Syz,Szz,                           &
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,                      &
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,                      &
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,                      &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,                                        &
               Hcon,Mxcon,Mycon,Mzcon,Gmxcon,Gmycon,Gmzcon,                    &
               Symmetry,Lev,eps,co)  result(gont)
! Wrapper function: dispatches to compute_rhs_z4c_opt.
! Mirrors compute_rhs_bssn / compute_rhs_bssn_opt in bssn_rhs.f90.
! A serial _cor reference can be added later for opt-vs-cor comparison.
  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3), Symmetry,Lev,co
  real*8, intent(in ):: T
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: chi,dxx,dyy,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: trK
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: gxy,gxz,gyz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamx,Gamy,Gamz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Lap, betax, betay, betaz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: dtSfx,  dtSfy,  dtSfz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: TZ
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: chi_rhs,trK_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gxx_rhs,gxy_rhs,gxz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gyy_rhs,gyz_rhs,gzz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Axx_rhs,Axy_rhs,Axz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Ayy_rhs,Ayz_rhs,Azz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamx_rhs,Gamy_rhs,Gamz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Lap_rhs, betax_rhs, betay_rhs, betaz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: dtSfx_rhs,dtSfy_rhs,dtSfz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: TZ_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: rho,Sx,Sy,Sz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Sxx,Sxy,Sxz,Syy,Syz,Szz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamzyy, Gamzyz, Gamzzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Rxx,Rxy,Rxz,Ryy,Ryz,Rzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Hcon,Mxcon,Mycon,Mzcon,Gmxcon,Gmycon,Gmzcon
  real*8,intent(in) :: eps
  integer :: gont, gont_opt

!~~~~~~> Dispatch to the optimized version:

  gont_opt = compute_rhs_z4c_opt(ex, T,X, Y, Z,                                &
               chi    ,   trK    ,                                             &
               dxx    ,   gxy    ,   gxz    ,   dyy    ,   gyz    ,   dzz,     &
               Axx    ,   Axy    ,   Axz    ,   Ayy    ,   Ayz    ,   Azz,     &
               Gamx   ,  Gamy    ,  Gamz    ,                                  &
               Lap    ,  betax   ,  betay   ,  betaz   ,                       &
               dtSfx  ,  dtSfy   ,  dtSfz   ,                                  &
               TZ     ,                                                        &
               chi_rhs,   trK_rhs,                                             &
               gxx_rhs,   gxy_rhs,   gxz_rhs,   gyy_rhs,   gyz_rhs,   gzz_rhs, &
               Axx_rhs,   Axy_rhs,   Axz_rhs,   Ayy_rhs,   Ayz_rhs,   Azz_rhs, &
               Gamx_rhs,  Gamy_rhs,  Gamz_rhs,                                 &
               Lap_rhs,  betax_rhs,  betay_rhs,  betaz_rhs,                    &
               dtSfx_rhs,  dtSfy_rhs,  dtSfz_rhs,                              &
               TZ_rhs   ,                                                      &
               rho,Sx,Sy,Sz,Sxx,Sxy,Sxz,Syy,Syz,Szz,                           &
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,                      &
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,                      &
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,                      &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,                                        &
               Hcon,Mxcon,Mycon,Mzcon,Gmxcon,Gmycon,Gmzcon,                    &
               Symmetry,Lev,eps,co)
  gont = gont_opt
  return

contains

  function compute_rhs_z4c_opt(ex, T,X, Y, Z,                                  &
               chi    ,   trK    ,                                             &
               dxx    ,   gxy    ,   gxz    ,   dyy    ,   gyz    ,   dzz,     &
               Axx    ,   Axy    ,   Axz    ,   Ayy    ,   Ayz    ,   Azz,     &
               Gamx   ,  Gamy    ,  Gamz    ,                                  &
               Lap    ,  betax   ,  betay   ,  betaz   ,                       &
               dtSfx  ,  dtSfy   ,  dtSfz   ,                                  &
               TZ     ,                                                        &
               chi_rhs,   trK_rhs,                                             &
               gxx_rhs,   gxy_rhs,   gxz_rhs,   gyy_rhs,   gyz_rhs,   gzz_rhs, &
               Axx_rhs,   Axy_rhs,   Axz_rhs,   Ayy_rhs,   Ayz_rhs,   Azz_rhs, &
               Gamx_rhs,  Gamy_rhs,  Gamz_rhs,                                 &
               Lap_rhs,  betax_rhs,  betay_rhs,  betaz_rhs,                    &
               dtSfx_rhs,  dtSfy_rhs,  dtSfz_rhs,                              &
               TZ_rhs   ,                                                      &
               rho,Sx,Sy,Sz,Sxx,Sxy,Sxz,Syy,Syz,Szz,                           &
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,                      &
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,                      &
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,                      &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,                                        &
               Hcon,Mxcon,Mycon,Mzcon,Gmxcon,Gmycon,Gmzcon,                    &
! co is not used here, we always compute constraint
               Symmetry,Lev,eps,co)  result(gont)
  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3), Symmetry,Lev,co
  real*8, intent(in ):: T
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: chi,dxx,dyy,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: trK
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: gxy,gxz,gyz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamx,Gamy,Gamz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Lap, betax, betay, betaz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: dtSfx,  dtSfy,  dtSfz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: TZ
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: chi_rhs,trK_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gxx_rhs,gxy_rhs,gxz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gyy_rhs,gyz_rhs,gzz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Axx_rhs,Axy_rhs,Axz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Ayy_rhs,Ayz_rhs,Azz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamx_rhs,Gamy_rhs,Gamz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Lap_rhs, betax_rhs, betay_rhs, betaz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: dtSfx_rhs,dtSfy_rhs,dtSfz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: TZ_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: rho,Sx,Sy,Sz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Sxx,Sxy,Sxz,Syy,Syz,Szz
! when out, physical second kind of connection  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamzyy, Gamzyz, Gamzzz
! when out, physical Ricci tensor  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Rxx,Rxy,Rxz,Ryy,Ryz,Rzz
! when out, constraint violation  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Hcon,Mxcon,Mycon,Mzcon,Gmxcon,Gmycon,Gmzcon
  real*8,intent(in) :: eps
!  gont = 0: success; gont = 1: something wrong
  integer::gont

!~~~~~~> Other variables:

  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gyy,gzz,alpn1,chin1
  real*8, dimension(ex(1),ex(2),ex(3)) :: trKd
  real*8, dimension(ex(1),ex(2),ex(3)) :: chix,chiy,chiz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxx,gxyx,gxzx,gyyx,gyzx,gzzx
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxy,gxyy,gxzy,gyyy,gyzy,gzzy
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxz,gxyz,gxzz,gyyz,gyzz,gzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Lapx,Lapy,Lapz
  real*8, dimension(ex(1),ex(2),ex(3)) :: betaxx,betaxy,betaxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: betayx,betayy,betayz
  real*8, dimension(ex(1),ex(2),ex(3)) :: betazx,betazy,betazz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamxx,Gamxy,Gamxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamyx,Gamyy,Gamyz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamzx,Gamzy,Gamzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Kx,Ky,Kz,div_beta,S
  real*8, dimension(ex(1),ex(2),ex(3)) :: f,fxx,fxy,fxz,fyy,fyz,fzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxxx,gxxxy,gxxxz,gxxyy,gxxyz,gxxzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxyxx,gxyxy,gxyxz,gxyyy,gxyyz,gxyzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxzxx,gxzxy,gxzxz,gxzyy,gxzyz,gxzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyyxx,gyyxy,gyyxz,gyyyy,gyyyz,gyyzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyzxx,gyzxy,gyzxz,gyzyy,gyzyz,gyzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gzzxx,gzzxy,gzzxz,gzzyy,gzzyz,gzzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamxa,Gamya,Gamza
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupyy,gupyz,gupzz

  real*8,dimension(3) ::SSS,AAS,ASA,SAA,ASS,SAS,SSA
  real*8            :: dX, dY, dZ, PI
  real*8, parameter :: ZEO=0.d0,ONE = 1.D0, TWO = 2.D0, FOUR = 4.D0,F16=1.6d1
  real*8, parameter :: EIGHT = 8.D0, HALF = 0.5D0, THR = 3.d0,F8=8.d0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  real*8, parameter :: F1o3 = 1.D0/3.D0, F2o3 = 2.D0/3.D0,F3o2=1.5d0, F1o6 = 1.D0/6.D0

! constraint damping terms stuffs PRD 81, 084003 (2010)
  real*8 :: kappa1,kappa2,kappa3,FF,eta
  real*8 :: tz_acc

  integer :: i, j, k

  call get_Z4cparameters(kappa1,kappa2,kappa3,FF,eta)

!!! sanity check (disabled for production)
!  dX = sum(chi)+sum(trK)+sum(dxx)+sum(gxy)+sum(gxz)+sum(dyy)+sum(gyz)+sum(dzz) &
!      +sum(Axx)+sum(Axy)+sum(Axz)+sum(Ayy)+sum(Ayz)+sum(Azz)                   &
!      +sum(Gamx)+sum(Gamy)+sum(Gamz)                                           &
!      +sum(Lap)+sum(betax)+sum(betay)+sum(betaz)+sum(dtSfx)+sum(dtSfy)+sum(dtSfz) &
!      +sum(TZ)
!  if(dX.ne.dX) then
!     if(sum(chi).ne.sum(chi))write(*,*)"Z4c_rhs.f90: find NaN in chi"
!     if(sum(trK).ne.sum(trK))write(*,*)"Z4c_rhs.f90: find NaN in trk"
!     if(sum(dxx).ne.sum(dxx))write(*,*)"Z4c_rhs.f90: find NaN in gxx"
!     if(sum(gxy).ne.sum(gxy))write(*,*)"Z4c_rhs.f90: find NaN in gxy"
!     if(sum(gxz).ne.sum(gxz))write(*,*)"Z4c_rhs.f90: find NaN in gxz"
!     if(sum(dyy).ne.sum(dyy))write(*,*)"Z4c_rhs.f90: find NaN in gyy"
!     if(sum(gyz).ne.sum(gyz))write(*,*)"Z4c_rhs.f90: find NaN in gyz"
!     if(sum(dzz).ne.sum(dzz))write(*,*)"Z4c_rhs.f90: find NaN in gzz"
!     if(sum(Axx).ne.sum(Axx))write(*,*)"Z4c_rhs.f90: find NaN in Axx"
!     if(sum(Axy).ne.sum(Axy))write(*,*)"Z4c_rhs.f90: find NaN in Axy"
!     if(sum(Axz).ne.sum(Axz))write(*,*)"Z4c_rhs.f90: find NaN in Axz"
!     if(sum(Ayy).ne.sum(Ayy))write(*,*)"Z4c_rhs.f90: find NaN in Ayy"
!     if(sum(Ayz).ne.sum(Ayz))write(*,*)"Z4c_rhs.f90: find NaN in Ayz"
!     if(sum(Azz).ne.sum(Azz))write(*,*)"Z4c_rhs.f90: find NaN in Azz"
!     if(sum(Gamx).ne.sum(Gamx))write(*,*)"Z4c_rhs.f90: find NaN in Gamx"
!     if(sum(Gamy).ne.sum(Gamy))write(*,*)"Z4c_rhs.f90: find NaN in Gamy"
!     if(sum(Gamz).ne.sum(Gamz))write(*,*)"Z4c_rhs.f90: find NaN in Gamz"
!     if(sum(Lap).ne.sum(Lap))write(*,*)"Z4c_rhs.f90: find NaN in Lap"
!     if(sum(betax).ne.sum(betax))write(*,*)"Z4c_rhs.f90: find NaN in betax"
!     if(sum(betay).ne.sum(betay))write(*,*)"Z4c_rhs.f90: find NaN in betay"
!     if(sum(betaz).ne.sum(betaz))write(*,*)"Z4c_rhs.f90: find NaN in betaz"
!     if(sum(dtSfx).ne.sum(dtSfx))write(*,*)"Z4c_rhs.f90: find NaN in dtSfx"
!     if(sum(dtSfy).ne.sum(dtSfy))write(*,*)"Z4c_rhs.f90: find NaN in dtSfy"
!     if(sum(dtSfz).ne.sum(dtSfz))write(*,*)"Z4c_rhs.f90: find NaN in dtSfz"
!     if(sum(TZ).ne.sum(Tz))write(*,*)"Z4c_rhs.f90: find NaN in TZ"
!     gont = 1
!     return
!  endif

  PI = dacos(-ONE)

  dX = X(2) - X(1)
  dY = Y(2) - Y(1)
  dZ = Z(2) - Z(1)

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        alpn1(i,j,k) = Lap(i,j,k) + ONE
        chin1(i,j,k) = chi(i,j,k) + ONE
        gxx(i,j,k)   = dxx(i,j,k) + ONE
        gyy(i,j,k)   = dyy(i,j,k) + ONE
        gzz(i,j,k)   = dzz(i,j,k) + ONE
        trKd(i,j,k)  = trK(i,j,k) + TWO * TZ(i,j,k)
      end do
    end do
  end do
!$OMP END PARALLEL DO

!this beta^i_,j will be kept till the end of this routine
  call fderivs(ex,betax,betaxx,betaxy,betaxz,X,Y,Z,ANTI, SYM, SYM,Symmetry,Lev)
  call fderivs(ex,betay,betayx,betayy,betayz,X,Y,Z, SYM,ANTI, SYM,Symmetry,Lev)
  call fderivs(ex,betaz,betazx,betazy,betazz,X,Y,Z, SYM, SYM,ANTI,Symmetry,Lev)

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        div_beta(i,j,k) = betaxx(i,j,k) + betayy(i,j,k) + betazz(i,j,k)
      end do
    end do
  end do
!$OMP END PARALLEL DO

  call fderivs(ex,chi,chix,chiy,chiz,X,Y,Z,SYM,SYM,SYM,symmetry,Lev)

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        chi_rhs(i,j,k) = F2o3 * chin1(i,j,k) * ( alpn1(i,j,k) * trKd(i,j,k) - div_beta(i,j,k) )
      end do
    end do
  end do
!$OMP END PARALLEL DO

  call fderivs(ex,dxx,gxxx,gxxy,gxxz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)
  call fderivs(ex,gxy,gxyx,gxyy,gxyz,X,Y,Z,ANTI,ANTI,SYM ,Symmetry,Lev)
  call fderivs(ex,gxz,gxzx,gxzy,gxzz,X,Y,Z,ANTI,SYM ,ANTI,Symmetry,Lev)
  call fderivs(ex,dyy,gyyx,gyyy,gyyz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)
  call fderivs(ex,gyz,gyzx,gyzy,gyzz,X,Y,Z,SYM ,ANTI,ANTI,Symmetry,Lev)
  call fderivs(ex,dzz,gzzx,gzzy,gzzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        gxx_rhs(i,j,k) = - TWO * alpn1(i,j,k) * Axx(i,j,k) - F2o3 * gxx(i,j,k) * div_beta(i,j,k) + &
                           TWO * ( gxx(i,j,k) * betaxx(i,j,k) + gxy(i,j,k) * betayx(i,j,k) + gxz(i,j,k) * betazx(i,j,k) )

        gyy_rhs(i,j,k) = - TWO * alpn1(i,j,k) * Ayy(i,j,k) - F2o3 * gyy(i,j,k) * div_beta(i,j,k) + &
                           TWO * ( gxy(i,j,k) * betaxy(i,j,k) + gyy(i,j,k) * betayy(i,j,k) + gyz(i,j,k) * betazy(i,j,k) )

        gzz_rhs(i,j,k) = - TWO * alpn1(i,j,k) * Azz(i,j,k) - F2o3 * gzz(i,j,k) * div_beta(i,j,k) + &
                           TWO * ( gxz(i,j,k) * betaxz(i,j,k) + gyz(i,j,k) * betayz(i,j,k) + gzz(i,j,k) * betazz(i,j,k) )

        gxy_rhs(i,j,k) = - TWO * alpn1(i,j,k) * Axy(i,j,k) + F1o3 * gxy(i,j,k) * div_beta(i,j,k) + &
                           gxx(i,j,k) * betaxy(i,j,k) + gxz(i,j,k) * betazy(i,j,k) + &
                           gyy(i,j,k) * betayx(i,j,k) + gyz(i,j,k) * betazx(i,j,k) - &
                           gxy(i,j,k) * betazz(i,j,k)

        gyz_rhs(i,j,k) = - TWO * alpn1(i,j,k) * Ayz(i,j,k) + F1o3 * gyz(i,j,k) * div_beta(i,j,k) + &
                           gxy(i,j,k) * betaxz(i,j,k) + gyy(i,j,k) * betayz(i,j,k) + &
                           gxz(i,j,k) * betaxy(i,j,k) + gzz(i,j,k) * betazy(i,j,k) - &
                           gyz(i,j,k) * betaxx(i,j,k)

        gxz_rhs(i,j,k) = - TWO * alpn1(i,j,k) * Axz(i,j,k) + F1o3 * gxz(i,j,k) * div_beta(i,j,k) + &
                           gxx(i,j,k) * betaxz(i,j,k) + gxy(i,j,k) * betayz(i,j,k) + &
                           gyz(i,j,k) * betayx(i,j,k) + gzz(i,j,k) * betazx(i,j,k) - &
                           gxz(i,j,k) * betayy(i,j,k)
      end do
    end do
  end do
!$OMP END PARALLEL DO

! invert tilted metric
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        gupzz(i,j,k) =   gxx(i,j,k) * gyy(i,j,k) * gzz(i,j,k) &
                       + gxy(i,j,k) * gyz(i,j,k) * gxz(i,j,k) &
                       + gxz(i,j,k) * gxy(i,j,k) * gyz(i,j,k) &
                       - gxz(i,j,k) * gyy(i,j,k) * gxz(i,j,k) &
                       - gxy(i,j,k) * gxy(i,j,k) * gzz(i,j,k) &
                       - gxx(i,j,k) * gyz(i,j,k) * gyz(i,j,k)
        gupxx(i,j,k) =   ( gyy(i,j,k) * gzz(i,j,k) - gyz(i,j,k) * gyz(i,j,k) ) / gupzz(i,j,k)
        gupxy(i,j,k) = - ( gxy(i,j,k) * gzz(i,j,k) - gyz(i,j,k) * gxz(i,j,k) ) / gupzz(i,j,k)
        gupxz(i,j,k) =   ( gxy(i,j,k) * gyz(i,j,k) - gyy(i,j,k) * gxz(i,j,k) ) / gupzz(i,j,k)
        gupyy(i,j,k) =   ( gxx(i,j,k) * gzz(i,j,k) - gxz(i,j,k) * gxz(i,j,k) ) / gupzz(i,j,k)
        gupyz(i,j,k) = - ( gxx(i,j,k) * gyz(i,j,k) - gxy(i,j,k) * gxz(i,j,k) ) / gupzz(i,j,k)
        gupzz(i,j,k) =   ( gxx(i,j,k) * gyy(i,j,k) - gxy(i,j,k) * gxy(i,j,k) ) / gupzz(i,j,k)
      end do
    end do
  end do
!$OMP END PARALLEL DO
! gij_,kl will be stored till end of this routine
  call fdderivs(ex,dxx,gxxxx,gxxxy,gxxxz,gxxyy,gxxyz,gxxzz,X,Y,Z,SYM ,SYM ,SYM ,symmetry,Lev)
  call fdderivs(ex,dyy,gyyxx,gyyxy,gyyxz,gyyyy,gyyyz,gyyzz,X,Y,Z,SYM ,SYM ,SYM ,symmetry,Lev)
  call fdderivs(ex,dzz,gzzxx,gzzxy,gzzxz,gzzyy,gzzyz,gzzzz,X,Y,Z,SYM ,SYM ,SYM ,symmetry,Lev)
  call fdderivs(ex,gxy,gxyxx,gxyxy,gxyxz,gxyyy,gxyyz,gxyzz,X,Y,Z,ANTI,ANTI,SYM ,symmetry,Lev)
  call fdderivs(ex,gxz,gxzxx,gxzxy,gxzxz,gxzyy,gxzyz,gxzzz,X,Y,Z,ANTI,SYM ,ANTI,symmetry,Lev)
  call fdderivs(ex,gyz,gyzxx,gyzxy,gyzxz,gyzyy,gyzyz,gyzzz,X,Y,Z,SYM ,ANTI,ANTI,symmetry,Lev)
! second kind of connection
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        Gamxxx(i,j,k) = HALF*( gupxx(i,j,k)*gxxx(i,j,k) + gupxy(i,j,k)*(TWO*gxyx(i,j,k) - gxxy(i,j,k)) + gupxz(i,j,k)*(TWO*gxzx(i,j,k) - gxxz(i,j,k)) )
        Gamyxx(i,j,k) = HALF*( gupxy(i,j,k)*gxxx(i,j,k) + gupyy(i,j,k)*(TWO*gxyx(i,j,k) - gxxy(i,j,k)) + gupyz(i,j,k)*(TWO*gxzx(i,j,k) - gxxz(i,j,k)) )
        Gamzxx(i,j,k) = HALF*( gupxz(i,j,k)*gxxx(i,j,k) + gupyz(i,j,k)*(TWO*gxyx(i,j,k) - gxxy(i,j,k)) + gupzz(i,j,k)*(TWO*gxzx(i,j,k) - gxxz(i,j,k)) )

        Gamxyy(i,j,k) = HALF*( gupxx(i,j,k)*(TWO*gxyy(i,j,k) - gyyx(i,j,k)) + gupxy(i,j,k)*gyyy(i,j,k) + gupxz(i,j,k)*(TWO*gyzy(i,j,k) - gyyz(i,j,k)) )
        Gamyyy(i,j,k) = HALF*( gupxy(i,j,k)*(TWO*gxyy(i,j,k) - gyyx(i,j,k)) + gupyy(i,j,k)*gyyy(i,j,k) + gupyz(i,j,k)*(TWO*gyzy(i,j,k) - gyyz(i,j,k)) )
        Gamzyy(i,j,k) = HALF*( gupxz(i,j,k)*(TWO*gxyy(i,j,k) - gyyx(i,j,k)) + gupyz(i,j,k)*gyyy(i,j,k) + gupzz(i,j,k)*(TWO*gyzy(i,j,k) - gyyz(i,j,k)) )

        Gamxzz(i,j,k) = HALF*( gupxx(i,j,k)*(TWO*gxzz(i,j,k) - gzzx(i,j,k)) + gupxy(i,j,k)*(TWO*gyzz(i,j,k) - gzzy(i,j,k)) + gupxz(i,j,k)*gzzz(i,j,k) )
        Gamyzz(i,j,k) = HALF*( gupxy(i,j,k)*(TWO*gxzz(i,j,k) - gzzx(i,j,k)) + gupyy(i,j,k)*(TWO*gyzz(i,j,k) - gzzy(i,j,k)) + gupyz(i,j,k)*gzzz(i,j,k) )
        Gamzzz(i,j,k) = HALF*( gupxz(i,j,k)*(TWO*gxzz(i,j,k) - gzzx(i,j,k)) + gupyz(i,j,k)*(TWO*gyzz(i,j,k) - gzzy(i,j,k)) + gupzz(i,j,k)*gzzz(i,j,k) )

        Gamxxy(i,j,k) = HALF*( gupxx(i,j,k)*gxxy(i,j,k) + gupxy(i,j,k)*gyyx(i,j,k) + gupxz(i,j,k)*( gxzy(i,j,k) + gyzx(i,j,k) - gxyz(i,j,k) ) )
        Gamyxy(i,j,k) = HALF*( gupxy(i,j,k)*gxxy(i,j,k) + gupyy(i,j,k)*gyyx(i,j,k) + gupyz(i,j,k)*( gxzy(i,j,k) + gyzx(i,j,k) - gxyz(i,j,k) ) )
        Gamzxy(i,j,k) = HALF*( gupxz(i,j,k)*gxxy(i,j,k) + gupyz(i,j,k)*gyyx(i,j,k) + gupzz(i,j,k)*( gxzy(i,j,k) + gyzx(i,j,k) - gxyz(i,j,k) ) )

        Gamxxz(i,j,k) = HALF*( gupxx(i,j,k)*gxxz(i,j,k) + gupxy(i,j,k)*( gxyz(i,j,k) + gyzx(i,j,k) - gxzy(i,j,k) ) + gupxz(i,j,k)*gzzx(i,j,k) )
        Gamyxz(i,j,k) = HALF*( gupxy(i,j,k)*gxxz(i,j,k) + gupyy(i,j,k)*( gxyz(i,j,k) + gyzx(i,j,k) - gxzy(i,j,k) ) + gupyz(i,j,k)*gzzx(i,j,k) )
        Gamzxz(i,j,k) = HALF*( gupxz(i,j,k)*gxxz(i,j,k) + gupyz(i,j,k)*( gxyz(i,j,k) + gyzx(i,j,k) - gxzy(i,j,k) ) + gupzz(i,j,k)*gzzx(i,j,k) )

        Gamxyz(i,j,k) = HALF*( gupxx(i,j,k)*( gxyz(i,j,k) + gxzy(i,j,k) - gyzx(i,j,k) ) + gupxy(i,j,k)*gyyz(i,j,k) + gupxz(i,j,k)*gzzy(i,j,k) )
        Gamyyz(i,j,k) = HALF*( gupxy(i,j,k)*( gxyz(i,j,k) + gxzy(i,j,k) - gyzx(i,j,k) ) + gupyy(i,j,k)*gyyz(i,j,k) + gupyz(i,j,k)*gzzy(i,j,k) )
        Gamzyz(i,j,k) = HALF*( gupxz(i,j,k)*( gxyz(i,j,k) + gxzy(i,j,k) - gyzx(i,j,k) ) + gupyz(i,j,k)*gyyz(i,j,k) + gupzz(i,j,k)*gzzy(i,j,k) )

        ! the so called Gamma_d
        Gamxa(i,j,k) =        gupxx(i,j,k) * Gamxxx(i,j,k) + gupyy(i,j,k) * Gamxyy(i,j,k) + gupzz(i,j,k) * Gamxzz(i,j,k) + &
                        TWO*( gupxy(i,j,k) * Gamxxy(i,j,k) + gupxz(i,j,k) * Gamxxz(i,j,k) + gupyz(i,j,k) * Gamxyz(i,j,k) )
        Gamya(i,j,k) =        gupxx(i,j,k) * Gamyxx(i,j,k) + gupyy(i,j,k) * Gamyyy(i,j,k) + gupzz(i,j,k) * Gamyzz(i,j,k) + &
                        TWO*( gupxy(i,j,k) * Gamyxy(i,j,k) + gupxz(i,j,k) * Gamyxz(i,j,k) + gupyz(i,j,k) * Gamyyz(i,j,k) )
        Gamza(i,j,k) =        gupxx(i,j,k) * Gamzxx(i,j,k) + gupyy(i,j,k) * Gamzyy(i,j,k) + gupzz(i,j,k) * Gamzzz(i,j,k) + &
                        TWO*( gupxy(i,j,k) * Gamzxy(i,j,k) + gupxz(i,j,k) * Gamzxz(i,j,k) + gupyz(i,j,k) * Gamzyz(i,j,k) )

        !!!!!!!!!!!! because gij_,k will be overwrite later, we calculate TWO*d_k Z^k here
        ! use Gamma^i as more as possible
        Gmxcon(i,j,k) = Gamx(i,j,k) - Gamxa(i,j,k)
        Gmycon(i,j,k) = Gamy(i,j,k) - Gamya(i,j,k)
        Gmzcon(i,j,k) = Gamz(i,j,k) - Gamza(i,j,k)
      end do
    end do
  end do
!$OMP END PARALLEL DO

!Maple generated code for g^ki*g^jm*g^ln*g_mn,k*g_ij,l
! Inlined as triple-nested loop with scalar accumulator
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k,tz_acc)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
      tz_acc = 3*gupxz(i,j,k)**2*gupzz(i,j,k)*gxzz(i,j,k)**2+gupxx(i,j,k)*gupxz(i,j,k)**2*gxxz(i,j,k)**2+2*gxyx(i,j,k)*gupxy(i,j,k)**3*gxyy(i,j,k)+ &
           2*gxyx(i,j,k)*gupxy(i,j,k)**3*gyyx(i,j,k)+gupxx(i,j,k)**2*gupzz(i,j,k)*gxzx(i,j,k)**2+3*gupxx(i,j,k)*gupxy(i,j,k)**2*gxyx(i,j,k)**2+ &
           6*gxyx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gyzy(i,j,k)+gupxx(i,j,k)**2*gupyy(i,j,k)*gxyx(i,j,k)**2+ &
           2*gxyz(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)**2*gyyz(i,j,k)+2*gxxz(i,j,k)*gupxx(i,j,k)**2*gupyz(i,j,k)*gxyx(i,j,k)+ &
           gupxz(i,j,k)**2*gupyy(i,j,k)*gyzx(i,j,k)**2+2*gxxy(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxxz(i,j,k)+ &
           2*gyzx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*gzzx(i,j,k)+3*gupyy(i,j,k)*gupyz(i,j,k)**2*gyzy(i,j,k)**2+ &
           2*gyyy(i,j,k)*gupyz(i,j,k)**3*gzzz(i,j,k)+2*gxxz(i,j,k)*gupxz(i,j,k)**3*gxzz(i,j,k)+ &
           4*gxzy(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gxyx(i,j,k)+gupyy(i,j,k)*gupyz(i,j,k)**2*gyyz(i,j,k)**2
      tz_acc = tz_acc +2*gxxz(i,j,k)*gupxy(i,j,k)**2*gupzz(i,j,k)*gyzy(i,j,k)+4*gxyz(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxxx(i,j,k)+ &
           6*gxzz(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gyzy(i,j,k)+2*gxxy(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gxzz(i,j,k)+ &
           3*gupxy(i,j,k)**2*gupyy(i,j,k)*gxyy(i,j,k)**2+2*gxyz(i,j,k)*gupxx(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*gyzx(i,j,k)+ &
           4*gxyy(i,j,k)*gupxx(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gyzx(i,j,k)+6*gxyy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gxzz(i,j,k)+ &
           4*gxzz(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gyzx(i,j,k)+3*gupxx(i,j,k)*gupxz(i,j,k)**2*gxzx(i,j,k)**2+ &
           4*gxyz(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gxyx(i,j,k)+2*gxxz(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gxyz(i,j,k)+ &
           2*gxxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gyzz(i,j,k)+2*gxzx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gyyz(i,j,k)+ &
           gupyz(i,j,k)**2*gupzz(i,j,k)*gzzy(i,j,k)**2+gupxz(i,j,k)**2*gupzz(i,j,k)*gzzx(i,j,k)**2+ &
           gupyy(i,j,k)*gupzz(i,j,k)**2*gyzz(i,j,k)**2+2*gyzy(i,j,k)*gupyz(i,j,k)**3*gzzy(i,j,k)+gupxx(i,j,k)*gupzz(i,j,k)**2*gxzz(i,j,k)**2
      tz_acc = tz_acc +gupxx(i,j,k)*gupyz(i,j,k)**2*gxzy(i,j,k)**2+2*gxzx(i,j,k)*gupxz(i,j,k)**3*gzzx(i,j,k)+ &
           3*gupyz(i,j,k)**2*gupzz(i,j,k)*gyzz(i,j,k)**2+2*gyzy(i,j,k)*gupyz(i,j,k)**3*gyzz(i,j,k)+gupyy(i,j,k)**2*gupzz(i,j,k)*gyzy(i,j,k)**2+ &
           gupxy(i,j,k)**2*gupzz(i,j,k)*gyzx(i,j,k)**2+2*gyyz(i,j,k)*gupyz(i,j,k)**3*gyzz(i,j,k)+gupxy(i,j,k)**2*gupyy(i,j,k)*gyyx(i,j,k)**2+ &
           gupxx(i,j,k)*gupyz(i,j,k)**2*gxyz(i,j,k)**2+gupxx(i,j,k)*gupyy(i,j,k)**2*gxyy(i,j,k)**2+ &
           gupxy(i,j,k)**2*gupzz(i,j,k)*gxzy(i,j,k)**2+2*gxzx(i,j,k)*gupxz(i,j,k)**3*gxzz(i,j,k)+ &
           2*gyyx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gyzx(i,j,k)+gupxx(i,j,k)*gupxy(i,j,k)**2*gxxy(i,j,k)**2+ &
           2*gxxx(i,j,k)*gupxz(i,j,k)**3*gzzz(i,j,k)+2*gxxx(i,j,k)*gupxy(i,j,k)**3*gyyy(i,j,k)+gupxz(i,j,k)**2*gupyy(i,j,k)*gxyz(i,j,k)**2+ &
           2*gxyy(i,j,k)*gupxy(i,j,k)**3*gxxy(i,j,k)
      tz_acc = tz_acc +2*gxyy(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)**2*gyzy(i,j,k)+6*gxyy(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gxzx(i,j,k)+ &
           4*gxyy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gxyz(i,j,k)+2*gyzx(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*gzzy(i,j,k)+ &
           2*gxzy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gxyy(i,j,k)+4*gxzy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*gxzz(i,j,k)+ &
           2*gyyx(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gyzz(i,j,k)+6*gxyx(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gxzz(i,j,k)+ &
           2*gxyz(i,j,k)*gupxy(i,j,k)**2*gupzz(i,j,k)*gxzy(i,j,k)+2*gxyz(i,j,k)*gupxy(i,j,k)**2*gupyz(i,j,k)*gxyy(i,j,k)+ &
           2*gxyz(i,j,k)*gupxy(i,j,k)**2*gupxz(i,j,k)*gxxy(i,j,k)+2*gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gxyz(i,j,k)**2+ &
           4*gxyy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)**2*gzzz(i,j,k)+2*gxyy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)**2*gzzy(i,j,k)+ &
           4*gxyy(i,j,k)*gupxy(i,j,k)**2*gupyz(i,j,k)*gxzy(i,j,k)+2*gxyy(i,j,k)*gupxy(i,j,k)**2*gupxz(i,j,k)*gxxz(i,j,k)+ &
           4*gxyy(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)**2*gxxx(i,j,k)+2*gxyx(i,j,k)*gupxy(i,j,k)**2*gupxz(i,j,k)*gxzy(i,j,k)+ &
           2*gxyx(i,j,k)*gupxy(i,j,k)**2*gupyz(i,j,k)*gyzy(i,j,k)
      tz_acc = tz_acc +2*gxyx(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)**2*gxxy(i,j,k)+4*gyzz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)**2*gzzz(i,j,k)+ &
           4*gxzy(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gxzx(i,j,k)+2*gxzy(i,j,k)*gupxx(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*gyzx(i,j,k)+ &
           4*gxxx(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gyzx(i,j,k)+2*gxyx(i,j,k)*gupxx(i,j,k)**2*gupyz(i,j,k)*gxzx(i,j,k)+ &
           2*gxyx(i,j,k)*gupxy(i,j,k)**2*gupxz(i,j,k)*gxyz(i,j,k)+2*gxzy(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gyyz(i,j,k)+ &
           4*gxzy(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gyyy(i,j,k)+2*gxzy(i,j,k)*gupxx(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gyyx(i,j,k)+ &
           2*gyyx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gxzy(i,j,k)+2*gyyx(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gyyz(i,j,k)+ &
           2*gyyx(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gyzy(i,j,k)+4*gxzy(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*gyzz(i,j,k)+ &
           2*gyyx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gxzz(i,j,k)+2*gxyz(i,j,k)*gupxx(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*gxzy(i,j,k)+ &
           2*gxyy(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gzzy(i,j,k)
      tz_acc = tz_acc +2*gxyy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gzzx(i,j,k)+2*gxyy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gyzx(i,j,k)+ &
           2*gxyy(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gyyz(i,j,k)+2*gxyy(i,j,k)*gupxx(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gxzy(i,j,k)+ &
           2*gxxy(i,j,k)*gupxy(i,j,k)**2*gupxz(i,j,k)*gxzy(i,j,k)+2*gxxy(i,j,k)*gupxy(i,j,k)**2*gupyz(i,j,k)*gyzy(i,j,k)+ &
           2*gxxy(i,j,k)*gupxy(i,j,k)**2*gupyy(i,j,k)*gyyy(i,j,k)+2*gxxy(i,j,k)*gupxx(i,j,k)**2*gupyz(i,j,k)*gxzx(i,j,k)+ &
           2*gxxy(i,j,k)*gupxx(i,j,k)**2*gupyy(i,j,k)*gxyx(i,j,k)+2*gxxx(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)**2*gzzx(i,j,k)+ &
           4*gxxx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)**2*gyzz(i,j,k)+4*gxxx(i,j,k)*gupxy(i,j,k)**2*gupxz(i,j,k)*gyzy(i,j,k)+ &
           2*gxxx(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)**2*gyyx(i,j,k)+4*gxxx(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)**2*gxzz(i,j,k)+ &
           4*gxxx(i,j,k)*gupxx(i,j,k)**2*gupxz(i,j,k)*gxzx(i,j,k)+2*gxxx(i,j,k)*gupxx(i,j,k)**2*gupxz(i,j,k)*gxxz(i,j,k)+ &
           4*gxyz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)**2*gyzz(i,j,k)+2*gxyz(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)**2*gyzy(i,j,k)+ &
           2*gxzy(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*gyzy(i,j,k)
      tz_acc = tz_acc +2*gxyy(i,j,k)*gupxx(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gxyz(i,j,k)+6*gxzz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gyzz(i,j,k)+ &
           4*gxzy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gzzz(i,j,k)+gupyy(i,j,k)**3*gyyy(i,j,k)**2+ &
           2*gxzy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gzzy(i,j,k)+2*gxzy(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gzzx(i,j,k)+ &
           2*gxyz(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gxzz(i,j,k)+2*gxzy(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gxzz(i,j,k)+ &
           2*gyzy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gzzx(i,j,k)+2*gyzy(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gxzy(i,j,k)+ &
           6*gyzy(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gyzz(i,j,k)+4*gyzx(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gyzy(i,j,k)+ &
           4*gyzx(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gyzz(i,j,k)+2*gxxy(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*gxyy(i,j,k)+ &
           4*gyzx(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gzzz(i,j,k)+2*gyzx(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*gyzy(i,j,k)+ &
           2*gyyz(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gzzy(i,j,k)+2*gyyz(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gzzx(i,j,k)
      tz_acc = tz_acc +2*gyyz(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gyzz(i,j,k)+2*gyyz(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*gyzx(i,j,k)+ &
           2*gyyz(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gxzz(i,j,k)+2*gxxy(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gyzx(i,j,k)+ &
           4*gyyy(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gyzx(i,j,k)+2*gyyx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gzzx(i,j,k)+ &
           2*gxyz(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gyzz(i,j,k)+2*gxxz(i,j,k)*gupxz(i,j,k)**2*gupzz(i,j,k)*gzzz(i,j,k)+ &
           2*gxxz(i,j,k)*gupxz(i,j,k)**2*gupyz(i,j,k)*gyzz(i,j,k)+2*gxxz(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)**2*gxzy(i,j,k)+ &
           2*gxxz(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)**2*gxzx(i,j,k)+2*gxxz(i,j,k)*gupxy(i,j,k)**2*gupyz(i,j,k)*gyyy(i,j,k)+ &
           2*gxxz(i,j,k)*gupxx(i,j,k)**2*gupzz(i,j,k)*gxzx(i,j,k)+2*gxxy(i,j,k)*gupxz(i,j,k)**2*gupyz(i,j,k)*gzzz(i,j,k)+ &
           2*gxxy(i,j,k)*gupxz(i,j,k)**2*gupyy(i,j,k)*gyzz(i,j,k)+2*gxxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)**2*gxzz(i,j,k)+ &
           2*gzzx(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gzzy(i,j,k)+2*gyzz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gzzx(i,j,k)+ &
           2*gxzx(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*gzzx(i,j,k)+2*gyzx(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*gyzz(i,j,k)
      tz_acc = tz_acc +gupzz(i,j,k)**3*gzzz(i,j,k)**2+2*gxzz(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*gyzx(i,j,k)+ &
           6*gxzx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gyzy(i,j,k)+2*gxxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gzzy(i,j,k)+ &
           4*gxzz(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)**2*gyyy(i,j,k)+2*gxzy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)**2*gyzz(i,j,k)+ &
           2*gxzy(i,j,k)*gupxz(i,j,k)**2*gupyz(i,j,k)*gxzz(i,j,k)+2*gxzy(i,j,k)*gupxz(i,j,k)**2*gupyy(i,j,k)*gxyz(i,j,k)+ &
           2*gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gxzy(i,j,k)**2+4*gxzx(i,j,k)*gupxz(i,j,k)**2*gupzz(i,j,k)*gzzz(i,j,k)+ &
           2*gxzx(i,j,k)*gupxz(i,j,k)**2*gupyz(i,j,k)*gyzz(i,j,k)+2*gxyz(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*gzzx(i,j,k)+ &
           2*gxyz(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*gzzy(i,j,k)+2*gxyx(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gxyz(i,j,k)+ &
           2*gxzz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)**2*gyyz(i,j,k)+2*gxxy(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxzx(i,j,k)+ &
           2*gyyx(i,j,k)*gupxy(i,j,k)**2*gupxz(i,j,k)*gxzx(i,j,k)
      tz_acc = tz_acc +2*gxyx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gzzy(i,j,k)+2*gyzy(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gzzy(i,j,k)+ &
           2*gxyx(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gyzx(i,j,k)+2*gyyx(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)**2*gyzz(i,j,k)+ &
           2*gyyx(i,j,k)*gupxy(i,j,k)**2*gupyz(i,j,k)*gyzx(i,j,k)+2*gyyx(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)**2*gzzz(i,j,k)+ &
           2*gyyx(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)**2*gyyy(i,j,k)+2*gxyz(i,j,k)*gupxy(i,j,k)**2*gupzz(i,j,k)*gyzx(i,j,k)+ &
           2*gxyz(i,j,k)*gupxy(i,j,k)**2*gupyz(i,j,k)*gyyx(i,j,k)+2*gxyy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)**2*gyzz(i,j,k)+ &
           2*gxyy(i,j,k)*gupxy(i,j,k)**2*gupyz(i,j,k)*gyzx(i,j,k)+2*gxyy(i,j,k)*gupxy(i,j,k)**2*gupyy(i,j,k)*gyyx(i,j,k)+ &
           2*gxyx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)**2*gzzx(i,j,k)+2*gxyx(i,j,k)*gupxy(i,j,k)**2*gupyz(i,j,k)*gyyz(i,j,k)+ &
           4*gxzz(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)**2*gzzz(i,j,k)+2*gxzz(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)**2*gzzy(i,j,k)+ &
           2*gxzz(i,j,k)*gupxx(i,j,k)*gupzz(i,j,k)**2*gzzx(i,j,k)+6*gxyx(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxzx(i,j,k)+ &
           2*gxyz(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gyzx(i,j,k)
      tz_acc = tz_acc +2*gyyx(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)**2*gyzy(i,j,k)+2*gyyx(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gzzy(i,j,k)+ &
           2*gxxz(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gxyy(i,j,k)+2*gyzx(i,j,k)*gupxz(i,j,k)**2*gupyy(i,j,k)*gxzy(i,j,k)+ &
           4*gyzx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)**2*gxzx(i,j,k)+2*gyzx(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)**2*gyzz(i,j,k)+ &
           2*gyzx(i,j,k)*gupxz(i,j,k)**2*gupyz(i,j,k)*gxzz(i,j,k)+2*gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gyzx(i,j,k)**2+ &
           2*gyyz(i,j,k)*gupyz(i,j,k)**2*gupzz(i,j,k)*gzzz(i,j,k)+2*gyyz(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)**2*gyzy(i,j,k)+ &
           2*gyyz(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)**2*gyzx(i,j,k)+2*gyyz(i,j,k)*gupyy(i,j,k)**2*gupzz(i,j,k)*gyzy(i,j,k)+ &
           2*gyyz(i,j,k)*gupxy(i,j,k)**2*gupzz(i,j,k)*gxzx(i,j,k)+2*gyyy(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)**2*gzzy(i,j,k)+ &
           2*gyyy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)**2*gzzx(i,j,k)+4*gyyy(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)**2*gyzz(i,j,k)+ &
           4*gyyy(i,j,k)*gupyy(i,j,k)**2*gupyz(i,j,k)*gyzy(i,j,k)+2*gyyy(i,j,k)*gupyy(i,j,k)**2*gupyz(i,j,k)*gyyz(i,j,k)
      tz_acc = tz_acc +2*gxyz(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gyzy(i,j,k)+2*gxyz(i,j,k)*gupxx(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gyyx(i,j,k)+ &
           2*gzzx(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)**2*gzzz(i,j,k)+2*gxzy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gyzx(i,j,k)+ &
           2*gyzz(i,j,k)*gupyz(i,j,k)**2*gupzz(i,j,k)*gzzy(i,j,k)+2*gyzy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)**2*gzzx(i,j,k)+ &
           2*gyzx(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)**2*gzzy(i,j,k)+2*gyzx(i,j,k)*gupxz(i,j,k)**2*gupyz(i,j,k)*gzzx(i,j,k)+ &
           2*gxzz(i,j,k)*gupxz(i,j,k)**2*gupzz(i,j,k)*gzzx(i,j,k)+2*gxzz(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)**2*gyzz(i,j,k)+ &
           2*gxzy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)**2*gzzy(i,j,k)+2*gxzy(i,j,k)*gupxz(i,j,k)**2*gupyz(i,j,k)*gzzx(i,j,k)+ &
           2*gxzx(i,j,k)*gupxz(i,j,k)**2*gupyz(i,j,k)*gzzy(i,j,k)+2*gyzz(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)**2*gzzy(i,j,k)+ &
           2*gyzz(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)**2*gzzx(i,j,k)+4*gyzy(i,j,k)*gupyz(i,j,k)**2*gupzz(i,j,k)*gzzz(i,j,k)+ &
           2*gyzy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)**2*gyzx(i,j,k)+2*gyzy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)**2*gxzz(i,j,k)+ &
           2*gxzy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gyzz(i,j,k)+2*gxyx(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gxzy(i,j,k)
      tz_acc = tz_acc +gupxx(i,j,k)**3*gxxx(i,j,k)**2+2*gzzy(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)**2*gzzz(i,j,k)+ &
           6*gxyx(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*gxyy(i,j,k)+2*gxzz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)* gupzz(i,j,k)*gzzy(i,j,k)+ &
           6*gxyx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gyzz(i,j,k)+2*gxzx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)**2*gxzy(i,j,k)+ &
           2*gxyx(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*gyyx(i,j,k)+2*gxyx(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gzzx(i,j,k)+ &
           2*gxyx(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxxz(i,j,k)+4*gxyx(i,j,k)*gupxx(i,j,k)**2*gupxy(i,j,k)*gxxx(i,j,k)+ &
           2*gxyx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gyyz(i,j,k)+6*gxyy(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gyzy(i,j,k)+ &
           2*gxyx(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gyzx(i,j,k)+6*gxyy(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gyzz(i,j,k)+ &
           4*gxyz(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*gxzx(i,j,k)+2*gxyz(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*gxzz(i,j,k)+ &
           4*gxyx(i,j,k)*gupxy(i,j,k)**2*gupyy(i,j,k)*gyyy(i,j,k)+2*gxyz(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gyyz(i,j,k)
      tz_acc = tz_acc +4*gxyz(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gyyy(i,j,k)+2*gxyx(i,j,k)*gupxz(i,j,k)**2*gupyy(i,j,k)*gyzz(i,j,k)+ &
           2*gxyz(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*gyzz(i,j,k)+2*gxyx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)**2*gxzz(i,j,k)+ &
           4*gxyz(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*gyzy(i,j,k)+2*gxzx(i,j,k)*gupxy(i,j,k)**2*gupzz(i,j,k)*gyzy(i,j,k)+ &
           2*gxyz(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gxzx(i,j,k)+4*gxyx(i,j,k)*gupxz(i,j,k)**2*gupyz(i,j,k)*gzzz(i,j,k)+ &
           4*gxzx(i,j,k)*gupxy(i,j,k)**2*gupyz(i,j,k)*gyyy(i,j,k)+2*gyyz(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*gxzy(i,j,k)+ &
           2*gxyz(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gxzy(i,j,k)+2*gxyz(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gzzx(i,j,k)+ &
           4*gxyy(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)**2*gyyy(i,j,k)+2*gxyy(i,j,k)*gupxx(i,j,k)*gupyy(i,j,k)**2*gyyx(i,j,k)+ &
           2*gxyy(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)**2*gzzx(i,j,k)+2*gxyz(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gzzy(i,j,k)+ &
           2*gxyy(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)**2*gyyz(i,j,k)+4*gxyz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gzzz(i,j,k)+ &
           2*gxxy(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gxyz(i,j,k)
      tz_acc = tz_acc +2*gxzx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)**2*gxyz(i,j,k)+2*gxxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gyzy(i,j,k)+ &
           4*gxxx(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxzy(i,j,k)+2*gxxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gyyz(i,j,k)+ &
           2*gxxy(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gyzx(i,j,k)+2*gxxy(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gzzx(i,j,k)+ &
           2*gxzx(i,j,k)*gupxy(i,j,k)**2*gupxz(i,j,k)*gxyy(i,j,k)+2*gxxy(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gxzy(i,j,k)+ &
           2*gxyz(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)**2*gxxz(i,j,k)+2*gxxy(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*gyyx(i,j,k)+ &
           2*gxyz(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)**2*gyzx(i,j,k)+4*gxyz(i,j,k)*gupxz(i,j,k)**2*gupyz(i,j,k)*gxzz(i,j,k)+ &
           2*gxxz(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*gxzy(i,j,k)+2*gxxz(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*gxzz(i,j,k)+ &
           2*gxxx(i,j,k)*gupxx(i,j,k)**2*gupxy(i,j,k)*gxxy(i,j,k)+2*gxxz(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gyyx(i,j,k)+ &
           2*gxxz(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*gyzz(i,j,k)+2*gxxz(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*gyzx(i,j,k)
      TZ_rhs(i,j,k) = tz_acc +2*gxxz(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gyyz(i,j,k)+2*gxxz(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gyzx(i,j,k)+ &
           2*gxxz(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*gzzx(i,j,k)+2*gxxz(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gyzy(i,j,k)+ &
           2*gxzx(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*gxzy(i,j,k)+2*gxxz(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*gzzy(i,j,k)+ &
           6*gxzx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*gyzz(i,j,k)+2*gxzx(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*gyzx(i,j,k)+ &
           2*gxzx(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gyyx(i,j,k)+6*gxzx(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*gxzz(i,j,k)+ &
           2*gxxx(i,j,k)*gupxy(i,j,k)**2*gupxz(i,j,k)*gyyz(i,j,k)+2*gxzx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*gzzy(i,j,k)+ &
           2*gxzx(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gyzx(i,j,k)+2*gxxx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)**2*gzzy(i,j,k)+ &
           4*gxzy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)**2*gyzy(i,j,k)+2*gxzy(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)**2*gyzx(i,j,k)+ &
           2*gxzz(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)**2*gyyx(i,j,k)+4*gxyx(i,j,k)*gupxy(i,j,k)**2*gupxz(i,j,k)*gyzx(i,j,k)+ &
           2*gxyx(i,j,k)*gupxz(i,j,k)**2*gupyy(i,j,k)*gzzy(i,j,k)+2*gxyy(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)**2*gxzz(i,j,k)
      end do
    end do
  end do
!$OMP END PARALLEL DO

! Gami_,j will be kept till the end of this routine
  call fderivs(ex,Gamx,Gamxx,Gamxy,Gamxz,X,Y,Z,ANTI,SYM ,SYM ,Symmetry,Lev)
  call fderivs(ex,Gamy,Gamyx,Gamyy,Gamyz,X,Y,Z,SYM ,ANTI,SYM ,Symmetry,Lev)
  call fderivs(ex,Gamz,Gamzx,Gamzy,Gamzz,X,Y,Z,SYM ,SYM ,ANTI,Symmetry,Lev)

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        TZ_rhs(i,j,k) = chix(i,j,k)*Gmxcon(i,j,k)+chiy(i,j,k)*Gmycon(i,j,k)+chiz(i,j,k)*Gmzcon(i,j,k) &
          +chin1(i,j,k)*(Gamxx(i,j,k)+Gamyy(i,j,k)+Gamzz(i,j,k) -          &
          (TWO*(gupxz(i,j,k)*gupyz(i,j,k)*gyzxz(i,j,k)+gupxx(i,j,k)*gupyy(i,j,k)*gxyxy(i,j,k)+gupxy(i,j,k)*gupyz(i,j,k)*gxzyy(i,j,k)+ &
                gupxx(i,j,k)*gupxy(i,j,k)*gxxxy(i,j,k)+gupxx(i,j,k)*gupxz(i,j,k)*gxxxz(i,j,k)+gupxx(i,j,k)*gupxy(i,j,k)*gxyxx(i,j,k)+ &
                gupxx(i,j,k)*gupyz(i,j,k)*gxyxz(i,j,k)+gupxx(i,j,k)*gupxz(i,j,k)*gxzxx(i,j,k)+gupxx(i,j,k)*gupyz(i,j,k)*gxzxy(i,j,k)+ &
                gupxx(i,j,k)*gupzz(i,j,k)*gxzxz(i,j,k)+gupxy(i,j,k)*gupxz(i,j,k)*gxxyz(i,j,k)+gupxy(i,j,k)*gupyy(i,j,k)*gxyyy(i,j,k)+ &
                gupxy(i,j,k)*gupyz(i,j,k)*gxyyz(i,j,k)+gupxy(i,j,k)*gupxz(i,j,k)*gxzxy(i,j,k)+gupxy(i,j,k)*gupzz(i,j,k)*gxzyz(i,j,k)+ &
                gupxy(i,j,k)*gupxz(i,j,k)*gxyxz(i,j,k)+gupxz(i,j,k)*gupyy(i,j,k)*gxyyz(i,j,k)+gupxz(i,j,k)*gupyz(i,j,k)*gxyzz(i,j,k)+ &
                gupxz(i,j,k)*gupyz(i,j,k)*gxzyz(i,j,k)+gupxz(i,j,k)*gupzz(i,j,k)*gxzzz(i,j,k)+gupxy(i,j,k)*gupyy(i,j,k)*gyyxy(i,j,k)+ &
                gupxy(i,j,k)*gupyz(i,j,k)*gyyxz(i,j,k)+gupxy(i,j,k)*gupxz(i,j,k)*gyzxx(i,j,k)+gupxy(i,j,k)*gupyz(i,j,k)*gyzxy(i,j,k)+ &
                gupxy(i,j,k)*gupzz(i,j,k)*gyzxz(i,j,k)+gupyy(i,j,k)*gupyz(i,j,k)*gyyyz(i,j,k)+gupxz(i,j,k)*gupyy(i,j,k)*gyzxy(i,j,k)+ &
                gupyy(i,j,k)*gupyz(i,j,k)*gyzyy(i,j,k)+gupyy(i,j,k)*gupzz(i,j,k)*gyzyz(i,j,k)+gupyz(i,j,k)*gupzz(i,j,k)*gyzzz(i,j,k)+ &
                gupxz(i,j,k)*gupyz(i,j,k)*gzzxy(i,j,k)+gupxz(i,j,k)*gupzz(i,j,k)*gzzxz(i,j,k)+gupyz(i,j,k)*gupzz(i,j,k)*gzzyz(i,j,k)+ &
                gupxy(i,j,k)*gupxy(i,j,k)*gxyxy(i,j,k)+gupxz(i,j,k)*gupxz(i,j,k)*gxzxz(i,j,k)+gupyz(i,j,k)*gupyz(i,j,k)*gyzyz(i,j,k)) &
               +gupxx(i,j,k)*gupxx(i,j,k)*gxxxx(i,j,k)+gupxy(i,j,k)*gupxy(i,j,k)*gxxyy(i,j,k)+gupxz(i,j,k)*gupxz(i,j,k)*gxxzz(i,j,k)+ &
                gupxy(i,j,k)*gupxy(i,j,k)*gyyxx(i,j,k)+gupyy(i,j,k)*gupyy(i,j,k)*gyyyy(i,j,k)+gupyz(i,j,k)*gupyz(i,j,k)*gyyzz(i,j,k)+ &
                gupxz(i,j,k)*gupxz(i,j,k)*gzzxx(i,j,k)+gupyz(i,j,k)*gupyz(i,j,k)*gzzyy(i,j,k)+gupzz(i,j,k)*gupzz(i,j,k)*gzzzz(i,j,k))+&
               (gxx(i,j,k)*Gamxa(i,j,k)*Gamxa(i,j,k)+gyy(i,j,k)*Gamya(i,j,k)*Gamya(i,j,k)+gzz(i,j,k)*Gamza(i,j,k)*Gamza(i,j,k)       +&
           TWO*(gxy(i,j,k)*Gamxa(i,j,k)*Gamya(i,j,k)+gxz(i,j,k)*Gamxa(i,j,k)*Gamza(i,j,k)+gyz(i,j,k)*Gamya(i,j,k)*Gamza(i,j,k))) + TZ_rhs(i,j,k))
      end do
    end do
  end do
!$OMP END PARALLEL DO

! Raise indices of \tilde A_{ij} and store in R_ij

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        Rxx(i,j,k) =    gupxx(i,j,k) * gupxx(i,j,k) * Axx(i,j,k) + gupxy(i,j,k) * gupxy(i,j,k) * Ayy(i,j,k) + gupxz(i,j,k) * gupxz(i,j,k) * Azz(i,j,k) + &
            TWO*(gupxx(i,j,k) * gupxy(i,j,k) * Axy(i,j,k) + gupxx(i,j,k) * gupxz(i,j,k) * Axz(i,j,k) + gupxy(i,j,k) * gupxz(i,j,k) * Ayz(i,j,k))

        Ryy(i,j,k) =    gupxy(i,j,k) * gupxy(i,j,k) * Axx(i,j,k) + gupyy(i,j,k) * gupyy(i,j,k) * Ayy(i,j,k) + gupyz(i,j,k) * gupyz(i,j,k) * Azz(i,j,k) + &
            TWO*(gupxy(i,j,k) * gupyy(i,j,k) * Axy(i,j,k) + gupxy(i,j,k) * gupyz(i,j,k) * Axz(i,j,k) + gupyy(i,j,k) * gupyz(i,j,k) * Ayz(i,j,k))

        Rzz(i,j,k) =    gupxz(i,j,k) * gupxz(i,j,k) * Axx(i,j,k) + gupyz(i,j,k) * gupyz(i,j,k) * Ayy(i,j,k) + gupzz(i,j,k) * gupzz(i,j,k) * Azz(i,j,k) + &
            TWO*(gupxz(i,j,k) * gupyz(i,j,k) * Axy(i,j,k) + gupxz(i,j,k) * gupzz(i,j,k) * Axz(i,j,k) + gupyz(i,j,k) * gupzz(i,j,k) * Ayz(i,j,k))

        Rxy(i,j,k) =    gupxx(i,j,k) * gupxy(i,j,k) * Axx(i,j,k) + gupxy(i,j,k) * gupyy(i,j,k) * Ayy(i,j,k) + gupxz(i,j,k) * gupyz(i,j,k) * Azz(i,j,k) + &
                (gupxx(i,j,k) * gupyy(i,j,k)       + gupxy(i,j,k) * gupxy(i,j,k))* Axy(i,j,k)                       + &
                (gupxx(i,j,k) * gupyz(i,j,k)       + gupxz(i,j,k) * gupxy(i,j,k))* Axz(i,j,k)                       + &
                (gupxy(i,j,k) * gupyz(i,j,k)       + gupxz(i,j,k) * gupyy(i,j,k))* Ayz(i,j,k)

        Rxz(i,j,k) =    gupxx(i,j,k) * gupxz(i,j,k) * Axx(i,j,k) + gupxy(i,j,k) * gupyz(i,j,k) * Ayy(i,j,k) + gupxz(i,j,k) * gupzz(i,j,k) * Azz(i,j,k) + &
                (gupxx(i,j,k) * gupyz(i,j,k)       + gupxy(i,j,k) * gupxz(i,j,k))* Axy(i,j,k)                       + &
                (gupxx(i,j,k) * gupzz(i,j,k)       + gupxz(i,j,k) * gupxz(i,j,k))* Axz(i,j,k)                       + &
                (gupxy(i,j,k) * gupzz(i,j,k)       + gupxz(i,j,k) * gupyz(i,j,k))* Ayz(i,j,k)

        Ryz(i,j,k) =    gupxy(i,j,k) * gupxz(i,j,k) * Axx(i,j,k) + gupyy(i,j,k) * gupyz(i,j,k) * Ayy(i,j,k) + gupyz(i,j,k) * gupzz(i,j,k) * Azz(i,j,k) + &
                (gupxy(i,j,k) * gupyz(i,j,k)       + gupyy(i,j,k) * gupxz(i,j,k))* Axy(i,j,k)                       + &
                (gupxy(i,j,k) * gupzz(i,j,k)       + gupyz(i,j,k) * gupxz(i,j,k))* Axz(i,j,k)                       + &
                (gupyy(i,j,k) * gupzz(i,j,k)       + gupyz(i,j,k) * gupyz(i,j,k))* Ayz(i,j,k)
      end do
    end do
  end do
!$OMP END PARALLEL DO

! Right hand side for Gam^i without shift terms...
! Lap_,i will be kept till the end of this routine
  call fderivs(ex,Lap,Lapx,Lapy,Lapz,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)
! K_,i stored K_,i+TZ_,i/2 indeed, will be kept till the end of this routine  
  call fderivs(ex,trK,Kx,Ky,Kz,X,Y,Z,SYM,SYM,SYM,symmetry,Lev)
  call fderivs(ex,TZ,fxx,fxy,fxz,X,Y,Z,SYM,SYM,SYM,symmetry,Lev)

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        Kx(i,j,k) = Kx(i,j,k) + fxx(i,j,k)/TWO
        Ky(i,j,k) = Ky(i,j,k) + fxy(i,j,k)/TWO
        Kz(i,j,k) = Kz(i,j,k) + fxz(i,j,k)/TWO

        Gamx_rhs(i,j,k) = - TWO * (   Lapx(i,j,k) * Rxx(i,j,k) +   Lapy(i,j,k) * Rxy(i,j,k) +   Lapz(i,j,k) * Rxz(i,j,k) ) + &
             TWO * alpn1(i,j,k) * (                                                &
             -F3o2/chin1(i,j,k) * (   chix(i,j,k) * Rxx(i,j,k) +   chiy(i,j,k) * Rxy(i,j,k) +   chiz(i,j,k) * Rxz(i,j,k) ) - &
                   gupxx(i,j,k) * (   F2o3 * Kx(i,j,k)  +  EIGHT * PI * Sx(i,j,k)            ) - &
                   gupxy(i,j,k) * (   F2o3 * Ky(i,j,k)  +  EIGHT * PI * Sy(i,j,k)            ) - &
                   gupxz(i,j,k) * (   F2o3 * Kz(i,j,k)  +  EIGHT * PI * Sz(i,j,k)            ) + &
                             Gamxxx(i,j,k) * Rxx(i,j,k) + Gamxyy(i,j,k) * Ryy(i,j,k) + Gamxzz(i,j,k) * Rzz(i,j,k)   + &
                     TWO * ( Gamxxy(i,j,k) * Rxy(i,j,k) + Gamxxz(i,j,k) * Rxz(i,j,k) + Gamxyz(i,j,k) * Ryz(i,j,k) ) )

        Gamy_rhs(i,j,k) = - TWO * (   Lapx(i,j,k) * Rxy(i,j,k) +   Lapy(i,j,k) * Ryy(i,j,k) +   Lapz(i,j,k) * Ryz(i,j,k) ) + &
             TWO * alpn1(i,j,k) * (                                                &
             -F3o2/chin1(i,j,k) * (   chix(i,j,k) * Rxy(i,j,k) +  chiy(i,j,k) * Ryy(i,j,k) +    chiz(i,j,k) * Ryz(i,j,k) ) - &
                   gupxy(i,j,k) * (   F2o3 * Kx(i,j,k)  +  EIGHT * PI * Sx(i,j,k)            ) - &
                   gupyy(i,j,k) * (   F2o3 * Ky(i,j,k)  +  EIGHT * PI * Sy(i,j,k)            ) - &
                   gupyz(i,j,k) * (   F2o3 * Kz(i,j,k)  +  EIGHT * PI * Sz(i,j,k)            ) + &
                             Gamyxx(i,j,k) * Rxx(i,j,k) + Gamyyy(i,j,k) * Ryy(i,j,k) + Gamyzz(i,j,k) * Rzz(i,j,k)   + &
                     TWO * ( Gamyxy(i,j,k) * Rxy(i,j,k) + Gamyxz(i,j,k) * Rxz(i,j,k) + Gamyyz(i,j,k) * Ryz(i,j,k) ) )

        Gamz_rhs(i,j,k) = - TWO * (   Lapx(i,j,k) * Rxz(i,j,k) +   Lapy(i,j,k) * Ryz(i,j,k) +   Lapz(i,j,k) * Rzz(i,j,k) ) + &
             TWO * alpn1(i,j,k) * (                                                &
             -F3o2/chin1(i,j,k) * (   chix(i,j,k) * Rxz(i,j,k) +  chiy(i,j,k) * Ryz(i,j,k) +    chiz(i,j,k) * Rzz(i,j,k) ) - &
                   gupxz(i,j,k) * (   F2o3 * Kx(i,j,k)  +  EIGHT * PI * Sx(i,j,k)            ) - &
                   gupyz(i,j,k) * (   F2o3 * Ky(i,j,k)  +  EIGHT * PI * Sy(i,j,k)            ) - &
                   gupzz(i,j,k) * (   F2o3 * Kz(i,j,k)  +  EIGHT * PI * Sz(i,j,k)            ) + &
                             Gamzxx(i,j,k) * Rxx(i,j,k) + Gamzyy(i,j,k) * Ryy(i,j,k) + Gamzzz(i,j,k) * Rzz(i,j,k)   + &
                     TWO * ( Gamzxy(i,j,k) * Rxy(i,j,k) + Gamzxz(i,j,k) * Rxz(i,j,k) + Gamzyz(i,j,k) * Ryz(i,j,k) ) )
      end do
    end do
  end do
!$OMP END PARALLEL DO
         
    call fdderivs(ex,betax,gxxx,gxyx,gxzx,gyyx,gyzx,gzzx,&
                X,Y,Z,ANTI,SYM, SYM ,Symmetry,Lev)
    call fdderivs(ex,betay,gxxy,gxyy,gxzy,gyyy,gyzy,gzzy,&
                X,Y,Z,SYM ,ANTI,SYM ,Symmetry,Lev)
    call fdderivs(ex,betaz,gxxz,gxyz,gxzz,gyyz,gyzz,gzzz,&
                X,Y,Z,SYM ,SYM, ANTI,Symmetry,Lev)

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        fxx(i,j,k) = gxxx(i,j,k) + gxyy(i,j,k) + gxzz(i,j,k)
        fxy(i,j,k) = gxyx(i,j,k) + gyyy(i,j,k) + gyzz(i,j,k)
        fxz(i,j,k) = gxzx(i,j,k) + gyzy(i,j,k) + gzzz(i,j,k)

        Gamx_rhs(i,j,k) =               Gamx_rhs(i,j,k) +  F2o3 *  Gamxa(i,j,k) * div_beta(i,j,k)        - &
                           Gamxa(i,j,k) * betaxx(i,j,k) - Gamya(i,j,k) * betaxy(i,j,k) - Gamza(i,j,k) * betaxz(i,j,k)  + &
                   F1o3 * (gupxx(i,j,k) * fxx(i,j,k)    + gupxy(i,j,k) * fxy(i,j,k)    + gupxz(i,j,k) * fxz(i,j,k)    ) + &
                           gupxx(i,j,k) * gxxx(i,j,k)   + gupyy(i,j,k) * gyyx(i,j,k)   + gupzz(i,j,k) * gzzx(i,j,k)    + &
                    TWO * (gupxy(i,j,k) * gxyx(i,j,k)   + gupxz(i,j,k) * gxzx(i,j,k)   + gupyz(i,j,k) * gyzx(i,j,k)  )

        Gamy_rhs(i,j,k) =               Gamy_rhs(i,j,k) +  F2o3 *  Gamya(i,j,k) * div_beta(i,j,k)        - &
                           Gamxa(i,j,k) * betayx(i,j,k) - Gamya(i,j,k) * betayy(i,j,k) - Gamza(i,j,k) * betayz(i,j,k)  + &
                   F1o3 * (gupxy(i,j,k) * fxx(i,j,k)    + gupyy(i,j,k) * fxy(i,j,k)    + gupyz(i,j,k) * fxz(i,j,k)    ) + &
                           gupxx(i,j,k) * gxxy(i,j,k)   + gupyy(i,j,k) * gyyy(i,j,k)   + gupzz(i,j,k) * gzzy(i,j,k)    + &
                    TWO * (gupxy(i,j,k) * gxyy(i,j,k)   + gupxz(i,j,k) * gxzy(i,j,k)   + gupyz(i,j,k) * gyzy(i,j,k)  )

        Gamz_rhs(i,j,k) =               Gamz_rhs(i,j,k) +  F2o3 *  Gamza(i,j,k) * div_beta(i,j,k)        - &
                           Gamxa(i,j,k) * betazx(i,j,k) - Gamya(i,j,k) * betazy(i,j,k) - Gamza(i,j,k) * betazz(i,j,k)  + &
                   F1o3 * (gupxz(i,j,k) * fxx(i,j,k)    + gupyz(i,j,k) * fxy(i,j,k)    + gupzz(i,j,k) * fxz(i,j,k)    ) + &
                           gupxx(i,j,k) * gxxz(i,j,k)   + gupyy(i,j,k) * gyyz(i,j,k)   + gupzz(i,j,k) * gzzz(i,j,k)    + &
                    TWO * (gupxy(i,j,k) * gxyz(i,j,k)   + gupxz(i,j,k) * gxzz(i,j,k)   + gupyz(i,j,k) * gyzz(i,j,k)  )    !rhs for Gam^i
      end do
    end do
  end do
!$OMP END PARALLEL DO

!first kind of connection stored in gij,k
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        gxxx(i,j,k) = gxx(i,j,k) * Gamxxx(i,j,k) + gxy(i,j,k) * Gamyxx(i,j,k) + gxz(i,j,k) * Gamzxx(i,j,k)
        gxyx(i,j,k) = gxx(i,j,k) * Gamxxy(i,j,k) + gxy(i,j,k) * Gamyxy(i,j,k) + gxz(i,j,k) * Gamzxy(i,j,k)
        gxzx(i,j,k) = gxx(i,j,k) * Gamxxz(i,j,k) + gxy(i,j,k) * Gamyxz(i,j,k) + gxz(i,j,k) * Gamzxz(i,j,k)
        gyyx(i,j,k) = gxx(i,j,k) * Gamxyy(i,j,k) + gxy(i,j,k) * Gamyyy(i,j,k) + gxz(i,j,k) * Gamzyy(i,j,k)
        gyzx(i,j,k) = gxx(i,j,k) * Gamxyz(i,j,k) + gxy(i,j,k) * Gamyyz(i,j,k) + gxz(i,j,k) * Gamzyz(i,j,k)
        gzzx(i,j,k) = gxx(i,j,k) * Gamxzz(i,j,k) + gxy(i,j,k) * Gamyzz(i,j,k) + gxz(i,j,k) * Gamzzz(i,j,k)

        gxxy(i,j,k) = gxy(i,j,k) * Gamxxx(i,j,k) + gyy(i,j,k) * Gamyxx(i,j,k) + gyz(i,j,k) * Gamzxx(i,j,k)
        gxyy(i,j,k) = gxy(i,j,k) * Gamxxy(i,j,k) + gyy(i,j,k) * Gamyxy(i,j,k) + gyz(i,j,k) * Gamzxy(i,j,k)
        gxzy(i,j,k) = gxy(i,j,k) * Gamxxz(i,j,k) + gyy(i,j,k) * Gamyxz(i,j,k) + gyz(i,j,k) * Gamzxz(i,j,k)
        gyyy(i,j,k) = gxy(i,j,k) * Gamxyy(i,j,k) + gyy(i,j,k) * Gamyyy(i,j,k) + gyz(i,j,k) * Gamzyy(i,j,k)
        gyzy(i,j,k) = gxy(i,j,k) * Gamxyz(i,j,k) + gyy(i,j,k) * Gamyyz(i,j,k) + gyz(i,j,k) * Gamzyz(i,j,k)
        gzzy(i,j,k) = gxy(i,j,k) * Gamxzz(i,j,k) + gyy(i,j,k) * Gamyzz(i,j,k) + gyz(i,j,k) * Gamzzz(i,j,k)

        gxxz(i,j,k) = gxz(i,j,k) * Gamxxx(i,j,k) + gyz(i,j,k) * Gamyxx(i,j,k) + gzz(i,j,k) * Gamzxx(i,j,k)
        gxyz(i,j,k) = gxz(i,j,k) * Gamxxy(i,j,k) + gyz(i,j,k) * Gamyxy(i,j,k) + gzz(i,j,k) * Gamzxy(i,j,k)
        gxzz(i,j,k) = gxz(i,j,k) * Gamxxz(i,j,k) + gyz(i,j,k) * Gamyxz(i,j,k) + gzz(i,j,k) * Gamzxz(i,j,k)
        gyyz(i,j,k) = gxz(i,j,k) * Gamxyy(i,j,k) + gyz(i,j,k) * Gamyyy(i,j,k) + gzz(i,j,k) * Gamzyy(i,j,k)
        gyzz(i,j,k) = gxz(i,j,k) * Gamxyz(i,j,k) + gyz(i,j,k) * Gamyyz(i,j,k) + gzz(i,j,k) * Gamzyz(i,j,k)
        gzzz(i,j,k) = gxz(i,j,k) * Gamxzz(i,j,k) + gyz(i,j,k) * Gamyzz(i,j,k) + gzz(i,j,k) * Gamzzz(i,j,k)
      end do
    end do
  end do
!$OMP END PARALLEL DO

!compute Ricci tensor for tilted metric
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        Rxx(i,j,k) =   gupxx(i,j,k) * gxxxx(i,j,k) + gupyy(i,j,k) * gxxyy(i,j,k) + gupzz(i,j,k) * gxxzz(i,j,k) + &
              ( gupxy(i,j,k) * gxxxy(i,j,k) + gupxz(i,j,k) * gxxxz(i,j,k) + gupyz(i,j,k) * gxxyz(i,j,k) ) * TWO

        Ryy(i,j,k) =   gupxx(i,j,k) * gyyxx(i,j,k) + gupyy(i,j,k) * gyyyy(i,j,k) + gupzz(i,j,k) * gyyzz(i,j,k) + &
              ( gupxy(i,j,k) * gyyxy(i,j,k) + gupxz(i,j,k) * gyyxz(i,j,k) + gupyz(i,j,k) * gyyyz(i,j,k) ) * TWO

        Rzz(i,j,k) =   gupxx(i,j,k) * gzzxx(i,j,k) + gupyy(i,j,k) * gzzyy(i,j,k) + gupzz(i,j,k) * gzzzz(i,j,k) + &
              ( gupxy(i,j,k) * gzzxy(i,j,k) + gupxz(i,j,k) * gzzxz(i,j,k) + gupyz(i,j,k) * gzzyz(i,j,k) ) * TWO

        Rxy(i,j,k) =   gupxx(i,j,k) * gxyxx(i,j,k) + gupyy(i,j,k) * gxyyy(i,j,k) + gupzz(i,j,k) * gxyzz(i,j,k) + &
              ( gupxy(i,j,k) * gxyxy(i,j,k) + gupxz(i,j,k) * gxyxz(i,j,k) + gupyz(i,j,k) * gxyyz(i,j,k) ) * TWO

        Rxz(i,j,k) =   gupxx(i,j,k) * gxzxx(i,j,k) + gupyy(i,j,k) * gxzyy(i,j,k) + gupzz(i,j,k) * gxzzz(i,j,k) + &
              ( gupxy(i,j,k) * gxzxy(i,j,k) + gupxz(i,j,k) * gxzxz(i,j,k) + gupyz(i,j,k) * gxzyz(i,j,k) ) * TWO

        Ryz(i,j,k) =   gupxx(i,j,k) * gyzxx(i,j,k) + gupyy(i,j,k) * gyzyy(i,j,k) + gupzz(i,j,k) * gyzzz(i,j,k) + &
              ( gupxy(i,j,k) * gyzxy(i,j,k) + gupxz(i,j,k) * gyzxz(i,j,k) + gupyz(i,j,k) * gyzyz(i,j,k) ) * TWO
      end do
    end do
  end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        Rxx(i,j,k) =     - HALF * Rxx(i,j,k)                                   + &
               gxx(i,j,k) * Gamxx(i,j,k)+ gxy(i,j,k) * Gamyx(i,j,k)   +    gxz(i,j,k) * Gamzx(i,j,k) + &
             Gamxa(i,j,k) * gxxx(i,j,k) +  Gamya(i,j,k) * gxyx(i,j,k) +  Gamza(i,j,k) * gxzx(i,j,k)  + &
   gupxx(i,j,k) *(                                                  &
       TWO*(Gamxxx(i,j,k) * gxxx(i,j,k) + Gamyxx(i,j,k) * gxyx(i,j,k) + Gamzxx(i,j,k) * gxzx(i,j,k)) + &
            Gamxxx(i,j,k) * gxxx(i,j,k) + Gamyxx(i,j,k) * gxxy(i,j,k) + Gamzxx(i,j,k) * gxxz(i,j,k) )+ &
   gupxy(i,j,k) *(                                                  &
       TWO*(Gamxxx(i,j,k) * gxyx(i,j,k) + Gamyxx(i,j,k) * gyyx(i,j,k) + Gamzxx(i,j,k) * gyzx(i,j,k)  + &
            Gamxxy(i,j,k) * gxxx(i,j,k) + Gamyxy(i,j,k) * gxyx(i,j,k) + Gamzxy(i,j,k) * gxzx(i,j,k)) + &
            Gamxxy(i,j,k) * gxxx(i,j,k) + Gamyxy(i,j,k) * gxxy(i,j,k) + Gamzxy(i,j,k) * gxxz(i,j,k)  + &
            Gamxxx(i,j,k) * gxyx(i,j,k) + Gamyxx(i,j,k) * gxyy(i,j,k) + Gamzxx(i,j,k) * gxyz(i,j,k) )+ &
   gupxz(i,j,k) *(                                                  &
       TWO*(Gamxxx(i,j,k) * gxzx(i,j,k) + Gamyxx(i,j,k) * gyzx(i,j,k) + Gamzxx(i,j,k) * gzzx(i,j,k)  + &
            Gamxxz(i,j,k) * gxxx(i,j,k) + Gamyxz(i,j,k) * gxyx(i,j,k) + Gamzxz(i,j,k) * gxzx(i,j,k)) + &
            Gamxxz(i,j,k) * gxxx(i,j,k) + Gamyxz(i,j,k) * gxxy(i,j,k) + Gamzxz(i,j,k) * gxxz(i,j,k)  + &
            Gamxxx(i,j,k) * gxzx(i,j,k) + Gamyxx(i,j,k) * gxzy(i,j,k) + Gamzxx(i,j,k) * gxzz(i,j,k) )+ &
   gupyy(i,j,k) *(                                                  &
       TWO*(Gamxxy(i,j,k) * gxyx(i,j,k) + Gamyxy(i,j,k) * gyyx(i,j,k) + Gamzxy(i,j,k) * gyzx(i,j,k)) + &
            Gamxxy(i,j,k) * gxyx(i,j,k) + Gamyxy(i,j,k) * gxyy(i,j,k) + Gamzxy(i,j,k) * gxyz(i,j,k) )+ &
   gupyz(i,j,k) *(                                                  &
       TWO*(Gamxxy(i,j,k) * gxzx(i,j,k) + Gamyxy(i,j,k) * gyzx(i,j,k) + Gamzxy(i,j,k) * gzzx(i,j,k)  + &
            Gamxxz(i,j,k) * gxyx(i,j,k) + Gamyxz(i,j,k) * gyyx(i,j,k) + Gamzxz(i,j,k) * gyzx(i,j,k)) + &
            Gamxxz(i,j,k) * gxyx(i,j,k) + Gamyxz(i,j,k) * gxyy(i,j,k) + Gamzxz(i,j,k) * gxyz(i,j,k)  + &
            Gamxxy(i,j,k) * gxzx(i,j,k) + Gamyxy(i,j,k) * gxzy(i,j,k) + Gamzxy(i,j,k) * gxzz(i,j,k) )+ &
   gupzz(i,j,k) *(                                                  &
       TWO*(Gamxxz(i,j,k) * gxzx(i,j,k) + Gamyxz(i,j,k) * gyzx(i,j,k) + Gamzxz(i,j,k) * gzzx(i,j,k)) + &
            Gamxxz(i,j,k) * gxzx(i,j,k) + Gamyxz(i,j,k) * gxzy(i,j,k) + Gamzxz(i,j,k) * gxzz(i,j,k) )

        Ryy(i,j,k) =     - HALF * Ryy(i,j,k)                                   + &
               gxy(i,j,k) * Gamxy(i,j,k)+  gyy(i,j,k) * Gamyy(i,j,k)  +  gyz(i,j,k) * Gamzy(i,j,k)   + &
             Gamxa(i,j,k) * gxyy(i,j,k) +  Gamya(i,j,k) * gyyy(i,j,k) +  Gamza(i,j,k) * gyzy(i,j,k)  + &
   gupxx(i,j,k) *(                                                  &
       TWO*(Gamxxy(i,j,k) * gxxy(i,j,k) + Gamyxy(i,j,k) * gxyy(i,j,k) + Gamzxy(i,j,k) * gxzy(i,j,k)) + &
            Gamxxy(i,j,k) * gxyx(i,j,k) + Gamyxy(i,j,k) * gxyy(i,j,k) + Gamzxy(i,j,k) * gxyz(i,j,k) )+ &
   gupxy(i,j,k) *(                                                  &
       TWO*(Gamxxy(i,j,k) * gxyy(i,j,k) + Gamyxy(i,j,k) * gyyy(i,j,k) + Gamzxy(i,j,k) * gyzy(i,j,k)  + &
            Gamxyy(i,j,k) * gxxy(i,j,k) + Gamyyy(i,j,k) * gxyy(i,j,k) + Gamzyy(i,j,k) * gxzy(i,j,k)) + &
            Gamxyy(i,j,k) * gxyx(i,j,k) + Gamyyy(i,j,k) * gxyy(i,j,k) + Gamzyy(i,j,k) * gxyz(i,j,k)  + &
            Gamxxy(i,j,k) * gyyx(i,j,k) + Gamyxy(i,j,k) * gyyy(i,j,k) + Gamzxy(i,j,k) * gyyz(i,j,k) )+ &
   gupxz(i,j,k) *(                                                  &
       TWO*(Gamxxy(i,j,k) * gxzy(i,j,k) + Gamyxy(i,j,k) * gyzy(i,j,k) + Gamzxy(i,j,k) * gzzy(i,j,k)  + &
            Gamxyz(i,j,k) * gxxy(i,j,k) + Gamyyz(i,j,k) * gxyy(i,j,k) + Gamzyz(i,j,k) * gxzy(i,j,k)) + &
            Gamxyz(i,j,k) * gxyx(i,j,k) + Gamyyz(i,j,k) * gxyy(i,j,k) + Gamzyz(i,j,k) * gxyz(i,j,k)  + &
            Gamxxy(i,j,k) * gyzx(i,j,k) + Gamyxy(i,j,k) * gyzy(i,j,k) + Gamzxy(i,j,k) * gyzz(i,j,k) )+ &
   gupyy(i,j,k) *(                                                  &
       TWO*(Gamxyy(i,j,k) * gxyy(i,j,k) + Gamyyy(i,j,k) * gyyy(i,j,k) + Gamzyy(i,j,k) * gyzy(i,j,k)) + &
            Gamxyy(i,j,k) * gyyx(i,j,k) + Gamyyy(i,j,k) * gyyy(i,j,k) + Gamzyy(i,j,k) * gyyz(i,j,k) )+ &
   gupyz(i,j,k) *(                                                  &
       TWO*(Gamxyy(i,j,k) * gxzy(i,j,k) + Gamyyy(i,j,k) * gyzy(i,j,k) + Gamzyy(i,j,k) * gzzy(i,j,k)  + &
            Gamxyz(i,j,k) * gxyy(i,j,k) + Gamyyz(i,j,k) * gyyy(i,j,k) + Gamzyz(i,j,k) * gyzy(i,j,k)) + &
            Gamxyz(i,j,k) * gyyx(i,j,k) + Gamyyz(i,j,k) * gyyy(i,j,k) + Gamzyz(i,j,k) * gyyz(i,j,k)  + &
            Gamxyy(i,j,k) * gyzx(i,j,k) + Gamyyy(i,j,k) * gyzy(i,j,k) + Gamzyy(i,j,k) * gyzz(i,j,k) )+ &
   gupzz(i,j,k) *(                                                  &
       TWO*(Gamxyz(i,j,k) * gxzy(i,j,k) + Gamyyz(i,j,k) * gyzy(i,j,k) + Gamzyz(i,j,k) * gzzy(i,j,k)) + &
            Gamxyz(i,j,k) * gyzx(i,j,k) + Gamyyz(i,j,k) * gyzy(i,j,k) + Gamzyz(i,j,k) * gyzz(i,j,k) )

        Rzz(i,j,k) =     - HALF * Rzz(i,j,k)                                   + &
               gxz(i,j,k) * Gamxz(i,j,k)+ gyz(i,j,k) * Gamyz(i,j,k)  +    gzz(i,j,k) * Gamzz(i,j,k)  + &
             Gamxa(i,j,k) * gxzz(i,j,k) +  Gamya(i,j,k) * gyzz(i,j,k) +  Gamza(i,j,k) * gzzz(i,j,k)  + &
   gupxx(i,j,k) *(                                                  &
       TWO*(Gamxxz(i,j,k) * gxxz(i,j,k) + Gamyxz(i,j,k) * gxyz(i,j,k) + Gamzxz(i,j,k) * gxzz(i,j,k)) + &
            Gamxxz(i,j,k) * gxzx(i,j,k) + Gamyxz(i,j,k) * gxzy(i,j,k) + Gamzxz(i,j,k) * gxzz(i,j,k) )+ &
   gupxy(i,j,k) *(                                                  &
       TWO*(Gamxxz(i,j,k) * gxyz(i,j,k) + Gamyxz(i,j,k) * gyyz(i,j,k) + Gamzxz(i,j,k) * gyzz(i,j,k)  + &
            Gamxyz(i,j,k) * gxxz(i,j,k) + Gamyyz(i,j,k) * gxyz(i,j,k) + Gamzyz(i,j,k) * gxzz(i,j,k)) + &
            Gamxyz(i,j,k) * gxzx(i,j,k) + Gamyyz(i,j,k) * gxzy(i,j,k) + Gamzyz(i,j,k) * gxzz(i,j,k)  + &
            Gamxxz(i,j,k) * gyzx(i,j,k) + Gamyxz(i,j,k) * gyzy(i,j,k) + Gamzxz(i,j,k) * gyzz(i,j,k) )+ &
   gupxz(i,j,k) *(                                                  &
       TWO*(Gamxxz(i,j,k) * gxzz(i,j,k) + Gamyxz(i,j,k) * gyzz(i,j,k) + Gamzxz(i,j,k) * gzzz(i,j,k)  + &
            Gamxzz(i,j,k) * gxxz(i,j,k) + Gamyzz(i,j,k) * gxyz(i,j,k) + Gamzzz(i,j,k) * gxzz(i,j,k)) + &
            Gamxzz(i,j,k) * gxzx(i,j,k) + Gamyzz(i,j,k) * gxzy(i,j,k) + Gamzzz(i,j,k) * gxzz(i,j,k)  + &
            Gamxxz(i,j,k) * gzzx(i,j,k) + Gamyxz(i,j,k) * gzzy(i,j,k) + Gamzxz(i,j,k) * gzzz(i,j,k) )+ &
   gupyy(i,j,k) *(                                                  &
       TWO*(Gamxyz(i,j,k) * gxyz(i,j,k) + Gamyyz(i,j,k) * gyyz(i,j,k) + Gamzyz(i,j,k) * gyzz(i,j,k)) + &
            Gamxyz(i,j,k) * gyzx(i,j,k) + Gamyyz(i,j,k) * gyzy(i,j,k) + Gamzyz(i,j,k) * gyzz(i,j,k) )+ &
   gupyz(i,j,k) *(                                                  &
       TWO*(Gamxyz(i,j,k) * gxzz(i,j,k) + Gamyyz(i,j,k) * gyzz(i,j,k) + Gamzyz(i,j,k) * gzzz(i,j,k)  + &
            Gamxzz(i,j,k) * gxyz(i,j,k) + Gamyzz(i,j,k) * gyyz(i,j,k) + Gamzzz(i,j,k) * gyzz(i,j,k)) + &
            Gamxzz(i,j,k) * gyzx(i,j,k) + Gamyzz(i,j,k) * gyzy(i,j,k) + Gamzzz(i,j,k) * gyzz(i,j,k)  + &
            Gamxyz(i,j,k) * gzzx(i,j,k) + Gamyyz(i,j,k) * gzzy(i,j,k) + Gamzyz(i,j,k) * gzzz(i,j,k) )+ &
   gupzz(i,j,k) *(                                                  &
       TWO*(Gamxzz(i,j,k) * gxzz(i,j,k) + Gamyzz(i,j,k) * gyzz(i,j,k) + Gamzzz(i,j,k) * gzzz(i,j,k)) + &
            Gamxzz(i,j,k) * gzzx(i,j,k) + Gamyzz(i,j,k) * gzzy(i,j,k) + Gamzzz(i,j,k) * gzzz(i,j,k) )

        Rxy(i,j,k) = HALF*(     - Rxy(i,j,k)                                   + &
               gxx(i,j,k) * Gamxy(i,j,k) +    gxy(i,j,k) * Gamyy(i,j,k) + gxz(i,j,k) * Gamzy(i,j,k)  + &
               gxy(i,j,k) * Gamxx(i,j,k) +    gyy(i,j,k) * Gamyx(i,j,k) + gyz(i,j,k) * Gamzx(i,j,k)  + &
             Gamxa(i,j,k) * gxyx(i,j,k) +  Gamya(i,j,k) * gyyx(i,j,k) +  Gamza(i,j,k) * gyzx(i,j,k)  + &
             Gamxa(i,j,k) * gxxy(i,j,k) +  Gamya(i,j,k) * gxyy(i,j,k) +  Gamza(i,j,k) * gxzy(i,j,k) )+ &
   gupxx(i,j,k) *(                                                  &
            Gamxxx(i,j,k) * gxxy(i,j,k) + Gamyxx(i,j,k) * gxyy(i,j,k) + Gamzxx(i,j,k) * gxzy(i,j,k)  + &
            Gamxxy(i,j,k) * gxxx(i,j,k) + Gamyxy(i,j,k) * gxyx(i,j,k) + Gamzxy(i,j,k) * gxzx(i,j,k)  + &
            Gamxxx(i,j,k) * gxyx(i,j,k) + Gamyxx(i,j,k) * gxyy(i,j,k) + Gamzxx(i,j,k) * gxyz(i,j,k) )+ &
   gupxy(i,j,k) *(                                                  &
            Gamxxx(i,j,k) * gxyy(i,j,k) + Gamyxx(i,j,k) * gyyy(i,j,k) + Gamzxx(i,j,k) * gyzy(i,j,k)  + &
            Gamxxy(i,j,k) * gxyx(i,j,k) + Gamyxy(i,j,k) * gyyx(i,j,k) + Gamzxy(i,j,k) * gyzx(i,j,k)  + &
            Gamxxy(i,j,k) * gxyx(i,j,k) + Gamyxy(i,j,k) * gxyy(i,j,k) + Gamzxy(i,j,k) * gxyz(i,j,k)  + &
            Gamxxy(i,j,k) * gxxy(i,j,k) + Gamyxy(i,j,k) * gxyy(i,j,k) + Gamzxy(i,j,k) * gxzy(i,j,k)  + &
            Gamxyy(i,j,k) * gxxx(i,j,k) + Gamyyy(i,j,k) * gxyx(i,j,k) + Gamzyy(i,j,k) * gxzx(i,j,k)  + &
            Gamxxx(i,j,k) * gyyx(i,j,k) + Gamyxx(i,j,k) * gyyy(i,j,k) + Gamzxx(i,j,k) * gyyz(i,j,k) )+ &
   gupxz(i,j,k) *(                                                  &
            Gamxxx(i,j,k) * gxzy(i,j,k) + Gamyxx(i,j,k) * gyzy(i,j,k) + Gamzxx(i,j,k) * gzzy(i,j,k)  + &
            Gamxxy(i,j,k) * gxzx(i,j,k) + Gamyxy(i,j,k) * gyzx(i,j,k) + Gamzxy(i,j,k) * gzzx(i,j,k)  + &
            Gamxxz(i,j,k) * gxyx(i,j,k) + Gamyxz(i,j,k) * gxyy(i,j,k) + Gamzxz(i,j,k) * gxyz(i,j,k)  + &
            Gamxxz(i,j,k) * gxxy(i,j,k) + Gamyxz(i,j,k) * gxyy(i,j,k) + Gamzxz(i,j,k) * gxzy(i,j,k)  + &
            Gamxyz(i,j,k) * gxxx(i,j,k) + Gamyyz(i,j,k) * gxyx(i,j,k) + Gamzyz(i,j,k) * gxzx(i,j,k)  + &
            Gamxxx(i,j,k) * gyzx(i,j,k) + Gamyxx(i,j,k) * gyzy(i,j,k) + Gamzxx(i,j,k) * gyzz(i,j,k) )+ &
   gupyy(i,j,k) *(                                                  &
            Gamxxy(i,j,k) * gxyy(i,j,k) + Gamyxy(i,j,k) * gyyy(i,j,k) + Gamzxy(i,j,k) * gyzy(i,j,k)  + &
            Gamxyy(i,j,k) * gxyx(i,j,k) + Gamyyy(i,j,k) * gyyx(i,j,k) + Gamzyy(i,j,k) * gyzx(i,j,k)  + &
            Gamxxy(i,j,k) * gyyx(i,j,k) + Gamyxy(i,j,k) * gyyy(i,j,k) + Gamzxy(i,j,k) * gyyz(i,j,k) )+ &
   gupyz(i,j,k) *(                                                  &
            Gamxxy(i,j,k) * gxzy(i,j,k) + Gamyxy(i,j,k) * gyzy(i,j,k) + Gamzxy(i,j,k) * gzzy(i,j,k)  + &
            Gamxyy(i,j,k) * gxzx(i,j,k) + Gamyyy(i,j,k) * gyzx(i,j,k) + Gamzyy(i,j,k) * gzzx(i,j,k)  + &
            Gamxxz(i,j,k) * gyyx(i,j,k) + Gamyxz(i,j,k) * gyyy(i,j,k) + Gamzxz(i,j,k) * gyyz(i,j,k)  + &
            Gamxxz(i,j,k) * gxyy(i,j,k) + Gamyxz(i,j,k) * gyyy(i,j,k) + Gamzxz(i,j,k) * gyzy(i,j,k)  + &
            Gamxyz(i,j,k) * gxyx(i,j,k) + Gamyyz(i,j,k) * gyyx(i,j,k) + Gamzyz(i,j,k) * gyzx(i,j,k)  + &
            Gamxxy(i,j,k) * gyzx(i,j,k) + Gamyxy(i,j,k) * gyzy(i,j,k) + Gamzxy(i,j,k) * gyzz(i,j,k) )+ &
   gupzz(i,j,k) *(                                                  &
            Gamxxz(i,j,k) * gxzy(i,j,k) + Gamyxz(i,j,k) * gyzy(i,j,k) + Gamzxz(i,j,k) * gzzy(i,j,k)  + &
            Gamxyz(i,j,k) * gxzx(i,j,k) + Gamyyz(i,j,k) * gyzx(i,j,k) + Gamzyz(i,j,k) * gzzx(i,j,k)  + &
            Gamxxz(i,j,k) * gyzx(i,j,k) + Gamyxz(i,j,k) * gyzy(i,j,k) + Gamzxz(i,j,k) * gyzz(i,j,k) )

        Rxz(i,j,k) = HALF*(     - Rxz(i,j,k)                                   + &
               gxx(i,j,k) * Gamxz(i,j,k) +  gxy(i,j,k) * Gamyz(i,j,k) + gxz(i,j,k) * Gamzz(i,j,k)    + &
               gxz(i,j,k) * Gamxx(i,j,k) +  gyz(i,j,k) * Gamyx(i,j,k) + gzz(i,j,k) * Gamzx(i,j,k)    + &
             Gamxa(i,j,k) * gxzx(i,j,k) +  Gamya(i,j,k) * gyzx(i,j,k) +  Gamza(i,j,k) * gzzx(i,j,k)  + &
             Gamxa(i,j,k) * gxxz(i,j,k) +  Gamya(i,j,k) * gxyz(i,j,k) +  Gamza(i,j,k) * gxzz(i,j,k) )+ &
   gupxx(i,j,k) *(                                                  &
            Gamxxx(i,j,k) * gxxz(i,j,k) + Gamyxx(i,j,k) * gxyz(i,j,k) + Gamzxx(i,j,k) * gxzz(i,j,k)  + &
            Gamxxz(i,j,k) * gxxx(i,j,k) + Gamyxz(i,j,k) * gxyx(i,j,k) + Gamzxz(i,j,k) * gxzx(i,j,k)  + &
            Gamxxx(i,j,k) * gxzx(i,j,k) + Gamyxx(i,j,k) * gxzy(i,j,k) + Gamzxx(i,j,k) * gxzz(i,j,k) )+ &
   gupxy(i,j,k) *(                                                  &
            Gamxxx(i,j,k) * gxyz(i,j,k) + Gamyxx(i,j,k) * gyyz(i,j,k) + Gamzxx(i,j,k) * gyzz(i,j,k)  + &
            Gamxxz(i,j,k) * gxyx(i,j,k) + Gamyxz(i,j,k) * gyyx(i,j,k) + Gamzxz(i,j,k) * gyzx(i,j,k)  + &
            Gamxxy(i,j,k) * gxzx(i,j,k) + Gamyxy(i,j,k) * gxzy(i,j,k) + Gamzxy(i,j,k) * gxzz(i,j,k)  + &
            Gamxxy(i,j,k) * gxxz(i,j,k) + Gamyxy(i,j,k) * gxyz(i,j,k) + Gamzxy(i,j,k) * gxzz(i,j,k)  + &
            Gamxyz(i,j,k) * gxxx(i,j,k) + Gamyyz(i,j,k) * gxyx(i,j,k) + Gamzyz(i,j,k) * gxzx(i,j,k)  + &
            Gamxxx(i,j,k) * gyzx(i,j,k) + Gamyxx(i,j,k) * gyzy(i,j,k) + Gamzxx(i,j,k) * gyzz(i,j,k) )+ &
   gupxz(i,j,k) *(                                                  &
            Gamxxx(i,j,k) * gxzz(i,j,k) + Gamyxx(i,j,k) * gyzz(i,j,k) + Gamzxx(i,j,k) * gzzz(i,j,k)  + &
            Gamxxz(i,j,k) * gxzx(i,j,k) + Gamyxz(i,j,k) * gyzx(i,j,k) + Gamzxz(i,j,k) * gzzx(i,j,k)  + &
            Gamxxz(i,j,k) * gxzx(i,j,k) + Gamyxz(i,j,k) * gxzy(i,j,k) + Gamzxz(i,j,k) * gxzz(i,j,k)  + &
            Gamxxz(i,j,k) * gxxz(i,j,k) + Gamyxz(i,j,k) * gxyz(i,j,k) + Gamzxz(i,j,k) * gxzz(i,j,k)  + &
            Gamxzz(i,j,k) * gxxx(i,j,k) + Gamyzz(i,j,k) * gxyx(i,j,k) + Gamzzz(i,j,k) * gxzx(i,j,k)  + &
            Gamxxx(i,j,k) * gzzx(i,j,k) + Gamyxx(i,j,k) * gzzy(i,j,k) + Gamzxx(i,j,k) * gzzz(i,j,k) )+ &
   gupyy(i,j,k) *(                                                  &
            Gamxxy(i,j,k) * gxyz(i,j,k) + Gamyxy(i,j,k) * gyyz(i,j,k) + Gamzxy(i,j,k) * gyzz(i,j,k)  + &
            Gamxyz(i,j,k) * gxyx(i,j,k) + Gamyyz(i,j,k) * gyyx(i,j,k) + Gamzyz(i,j,k) * gyzx(i,j,k)  + &
            Gamxxy(i,j,k) * gyzx(i,j,k) + Gamyxy(i,j,k) * gyzy(i,j,k) + Gamzxy(i,j,k) * gyzz(i,j,k) )+ &
   gupyz(i,j,k) *(                                                  &
            Gamxxy(i,j,k) * gxzz(i,j,k) + Gamyxy(i,j,k) * gyzz(i,j,k) + Gamzxy(i,j,k) * gzzz(i,j,k)  + &
            Gamxyz(i,j,k) * gxzx(i,j,k) + Gamyyz(i,j,k) * gyzx(i,j,k) + Gamzyz(i,j,k) * gzzx(i,j,k)  + &
            Gamxxz(i,j,k) * gyzx(i,j,k) + Gamyxz(i,j,k) * gyzy(i,j,k) + Gamzxz(i,j,k) * gyzz(i,j,k)  + &
            Gamxxz(i,j,k) * gxyz(i,j,k) + Gamyxz(i,j,k) * gyyz(i,j,k) + Gamzxz(i,j,k) * gyzz(i,j,k)  + &
            Gamxzz(i,j,k) * gxyx(i,j,k) + Gamyzz(i,j,k) * gyyx(i,j,k) + Gamzzz(i,j,k) * gyzx(i,j,k)  + &
            Gamxxy(i,j,k) * gzzx(i,j,k) + Gamyxy(i,j,k) * gzzy(i,j,k) + Gamzxy(i,j,k) * gzzz(i,j,k) )+ &
   gupzz(i,j,k) *(                                                  &
            Gamxxz(i,j,k) * gxzz(i,j,k) + Gamyxz(i,j,k) * gyzz(i,j,k) + Gamzxz(i,j,k) * gzzz(i,j,k)  + &
            Gamxzz(i,j,k) * gxzx(i,j,k) + Gamyzz(i,j,k) * gyzx(i,j,k) + Gamzzz(i,j,k) * gzzx(i,j,k)  + &
            Gamxxz(i,j,k) * gzzx(i,j,k) + Gamyxz(i,j,k) * gzzy(i,j,k) + Gamzxz(i,j,k) * gzzz(i,j,k) )

        Ryz(i,j,k) = HALF*(     - Ryz(i,j,k)                                   + &
               gxy(i,j,k) * Gamxz(i,j,k) + gyy(i,j,k) * Gamyz(i,j,k) + gyz(i,j,k) * Gamzz(i,j,k)     + &
               gxz(i,j,k) * Gamxy(i,j,k) + gyz(i,j,k) * Gamyy(i,j,k) + gzz(i,j,k) * Gamzy(i,j,k)     + &
             Gamxa(i,j,k) * gxzy(i,j,k) +  Gamya(i,j,k) * gyzy(i,j,k) +  Gamza(i,j,k) * gzzy(i,j,k)  + &
             Gamxa(i,j,k) * gxyz(i,j,k) +  Gamya(i,j,k) * gyyz(i,j,k) +  Gamza(i,j,k) * gyzz(i,j,k) )+ &
   gupxx(i,j,k) *(                                                  &
            Gamxxy(i,j,k) * gxxz(i,j,k) + Gamyxy(i,j,k) * gxyz(i,j,k) + Gamzxy(i,j,k) * gxzz(i,j,k)  + &
            Gamxxz(i,j,k) * gxxy(i,j,k) + Gamyxz(i,j,k) * gxyy(i,j,k) + Gamzxz(i,j,k) * gxzy(i,j,k)  + &
            Gamxxy(i,j,k) * gxzx(i,j,k) + Gamyxy(i,j,k) * gxzy(i,j,k) + Gamzxy(i,j,k) * gxzz(i,j,k) )+ &
   gupxy(i,j,k) *(                                                  &
            Gamxxy(i,j,k) * gxyz(i,j,k) + Gamyxy(i,j,k) * gyyz(i,j,k) + Gamzxy(i,j,k) * gyzz(i,j,k)  + &
            Gamxxz(i,j,k) * gxyy(i,j,k) + Gamyxz(i,j,k) * gyyy(i,j,k) + Gamzxz(i,j,k) * gyzy(i,j,k)  + &
            Gamxyy(i,j,k) * gxzx(i,j,k) + Gamyyy(i,j,k) * gxzy(i,j,k) + Gamzyy(i,j,k) * gxzz(i,j,k)  + &
            Gamxyy(i,j,k) * gxxz(i,j,k) + Gamyyy(i,j,k) * gxyz(i,j,k) + Gamzyy(i,j,k) * gxzz(i,j,k)  + &
            Gamxyz(i,j,k) * gxxy(i,j,k) + Gamyyz(i,j,k) * gxyy(i,j,k) + Gamzyz(i,j,k) * gxzy(i,j,k)  + &
            Gamxxy(i,j,k) * gyzx(i,j,k) + Gamyxy(i,j,k) * gyzy(i,j,k) + Gamzxy(i,j,k) * gyzz(i,j,k) )+ &
   gupxz(i,j,k) *(                                                  &
            Gamxxy(i,j,k) * gxzz(i,j,k) + Gamyxy(i,j,k) * gyzz(i,j,k) + Gamzxy(i,j,k) * gzzz(i,j,k)  + &
            Gamxxz(i,j,k) * gxzy(i,j,k) + Gamyxz(i,j,k) * gyzy(i,j,k) + Gamzxz(i,j,k) * gzzy(i,j,k)  + &
            Gamxyz(i,j,k) * gxzx(i,j,k) + Gamyyz(i,j,k) * gxzy(i,j,k) + Gamzyz(i,j,k) * gxzz(i,j,k)  + &
            Gamxyz(i,j,k) * gxxz(i,j,k) + Gamyyz(i,j,k) * gxyz(i,j,k) + Gamzyz(i,j,k) * gxzz(i,j,k)  + &
            Gamxzz(i,j,k) * gxxy(i,j,k) + Gamyzz(i,j,k) * gxyy(i,j,k) + Gamzzz(i,j,k) * gxzy(i,j,k)  + &
            Gamxxy(i,j,k) * gzzx(i,j,k) + Gamyxy(i,j,k) * gzzy(i,j,k) + Gamzxy(i,j,k) * gzzz(i,j,k) )+ &
   gupyy(i,j,k) *(                                                  &
            Gamxyy(i,j,k) * gxyz(i,j,k) + Gamyyy(i,j,k) * gyyz(i,j,k) + Gamzyy(i,j,k) * gyzz(i,j,k)  + &
            Gamxyz(i,j,k) * gxyy(i,j,k) + Gamyyz(i,j,k) * gyyy(i,j,k) + Gamzyz(i,j,k) * gyzy(i,j,k)  + &
            Gamxyy(i,j,k) * gyzx(i,j,k) + Gamyyy(i,j,k) * gyzy(i,j,k) + Gamzyy(i,j,k) * gyzz(i,j,k) )+ &
   gupyz(i,j,k) *(                                                  &
            Gamxyy(i,j,k) * gxzz(i,j,k) + Gamyyy(i,j,k) * gyzz(i,j,k) + Gamzyy(i,j,k) * gzzz(i,j,k)  + &
            Gamxyz(i,j,k) * gxzy(i,j,k) + Gamyyz(i,j,k) * gyzy(i,j,k) + Gamzyz(i,j,k) * gzzy(i,j,k)  + &
            Gamxyz(i,j,k) * gyzx(i,j,k) + Gamyyz(i,j,k) * gyzy(i,j,k) + Gamzyz(i,j,k) * gyzz(i,j,k)  + &
            Gamxyz(i,j,k) * gxyz(i,j,k) + Gamyyz(i,j,k) * gyyz(i,j,k) + Gamzyz(i,j,k) * gyzz(i,j,k)  + &
            Gamxzz(i,j,k) * gxyy(i,j,k) + Gamyzz(i,j,k) * gyyy(i,j,k) + Gamzzz(i,j,k) * gyzy(i,j,k)  + &
            Gamxyy(i,j,k) * gzzx(i,j,k) + Gamyyy(i,j,k) * gzzy(i,j,k) + Gamzyy(i,j,k) * gzzz(i,j,k) )+ &
   gupzz(i,j,k) *(                                                  &
            Gamxyz(i,j,k) * gxzz(i,j,k) + Gamyyz(i,j,k) * gyzz(i,j,k) + Gamzyz(i,j,k) * gzzz(i,j,k)  + &
            Gamxzz(i,j,k) * gxzy(i,j,k) + Gamyzz(i,j,k) * gyzy(i,j,k) + Gamzzz(i,j,k) * gzzy(i,j,k)  + &
            Gamxyz(i,j,k) * gzzx(i,j,k) + Gamyyz(i,j,k) * gzzy(i,j,k) + Gamzyz(i,j,k) * gzzz(i,j,k) )
      end do
    end do
  end do
!$OMP END PARALLEL DO
!covariant second derivative of chi respect to tilted metric

! Store D^l D_l chi - 3/(2*chi) D^l chi D_l chi in f

  call fdderivs(ex,chi,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        fxx(i,j,k) = fxx(i,j,k) - Gamxxx(i,j,k) * chix(i,j,k) - Gamyxx(i,j,k) * chiy(i,j,k) - Gamzxx(i,j,k) * chiz(i,j,k)
        fxy(i,j,k) = fxy(i,j,k) - Gamxxy(i,j,k) * chix(i,j,k) - Gamyxy(i,j,k) * chiy(i,j,k) - Gamzxy(i,j,k) * chiz(i,j,k)
        fxz(i,j,k) = fxz(i,j,k) - Gamxxz(i,j,k) * chix(i,j,k) - Gamyxz(i,j,k) * chiy(i,j,k) - Gamzxz(i,j,k) * chiz(i,j,k)
        fyy(i,j,k) = fyy(i,j,k) - Gamxyy(i,j,k) * chix(i,j,k) - Gamyyy(i,j,k) * chiy(i,j,k) - Gamzyy(i,j,k) * chiz(i,j,k)
        fyz(i,j,k) = fyz(i,j,k) - Gamxyz(i,j,k) * chix(i,j,k) - Gamyyz(i,j,k) * chiy(i,j,k) - Gamzyz(i,j,k) * chiz(i,j,k)
        fzz(i,j,k) = fzz(i,j,k) - Gamxzz(i,j,k) * chix(i,j,k) - Gamyzz(i,j,k) * chiy(i,j,k) - Gamzzz(i,j,k) * chiz(i,j,k)

        f(i,j,k) =        gupxx(i,j,k) * ( fxx(i,j,k) - F3o2/chin1(i,j,k) * chix(i,j,k) * chix(i,j,k) ) + &
                          gupyy(i,j,k) * ( fyy(i,j,k) - F3o2/chin1(i,j,k) * chiy(i,j,k) * chiy(i,j,k) ) + &
                          gupzz(i,j,k) * ( fzz(i,j,k) - F3o2/chin1(i,j,k) * chiz(i,j,k) * chiz(i,j,k) ) + &
                    TWO * gupxy(i,j,k) * ( fxy(i,j,k) - F3o2/chin1(i,j,k) * chix(i,j,k) * chiy(i,j,k) ) + &
                    TWO * gupxz(i,j,k) * ( fxz(i,j,k) - F3o2/chin1(i,j,k) * chix(i,j,k) * chiz(i,j,k) ) + &
                    TWO * gupyz(i,j,k) * ( fyz(i,j,k) - F3o2/chin1(i,j,k) * chiy(i,j,k) * chiz(i,j,k) )

! Add chi part to Ricci tensor:

        fxx(i,j,k) = Rxx(i,j,k) + (fxx(i,j,k) - chix(i,j,k)*chix(i,j,k)/chin1(i,j,k)/TWO + gxx(i,j,k) * f(i,j,k))/chin1(i,j,k)/TWO
        fyy(i,j,k) = Ryy(i,j,k) + (fyy(i,j,k) - chiy(i,j,k)*chiy(i,j,k)/chin1(i,j,k)/TWO + gyy(i,j,k) * f(i,j,k))/chin1(i,j,k)/TWO
        fzz(i,j,k) = Rzz(i,j,k) + (fzz(i,j,k) - chiz(i,j,k)*chiz(i,j,k)/chin1(i,j,k)/TWO + gzz(i,j,k) * f(i,j,k))/chin1(i,j,k)/TWO
        fxy(i,j,k) = Rxy(i,j,k) + (fxy(i,j,k) - chix(i,j,k)*chiy(i,j,k)/chin1(i,j,k)/TWO + gxy(i,j,k) * f(i,j,k))/chin1(i,j,k)/TWO
        fxz(i,j,k) = Rxz(i,j,k) + (fxz(i,j,k) - chix(i,j,k)*chiz(i,j,k)/chin1(i,j,k)/TWO + gxz(i,j,k) * f(i,j,k))/chin1(i,j,k)/TWO
        fyz(i,j,k) = Ryz(i,j,k) + (fyz(i,j,k) - chiy(i,j,k)*chiz(i,j,k)/chin1(i,j,k)/TWO + gyz(i,j,k) * f(i,j,k))/chin1(i,j,k)/TWO
! store R/chi in Hcon
        Hcon(i,j,k) =   gupxx(i,j,k) * fxx(i,j,k) + gupyy(i,j,k) * fyy(i,j,k) + gupzz(i,j,k) * fzz(i,j,k) + &
              TWO* ( gupxy(i,j,k) * fxy(i,j,k) + gupxz(i,j,k) * fxz(i,j,k) + gupyz(i,j,k) * fyz(i,j,k) )

        Rxx(i,j,k) = fxx(i,j,k)
        Ryy(i,j,k) = fyy(i,j,k)
        Rzz(i,j,k) = fzz(i,j,k)
        Rxy(i,j,k) = fxy(i,j,k)
        Rxz(i,j,k) = fxz(i,j,k)
        Ryz(i,j,k) = fyz(i,j,k)

        gxxx(i,j,k) = (gupxx(i,j,k) * chix(i,j,k) + gupxy(i,j,k) * chiy(i,j,k) + gupxz(i,j,k) * chiz(i,j,k))/chin1(i,j,k)
        gxxy(i,j,k) = (gupxy(i,j,k) * chix(i,j,k) + gupyy(i,j,k) * chiy(i,j,k) + gupyz(i,j,k) * chiz(i,j,k))/chin1(i,j,k)
        gxxz(i,j,k) = (gupxz(i,j,k) * chix(i,j,k) + gupyz(i,j,k) * chiy(i,j,k) + gupzz(i,j,k) * chiz(i,j,k))/chin1(i,j,k)
! now get physical second kind of connection
        Gamxxx(i,j,k) = Gamxxx(i,j,k) - ( (chix(i,j,k) + chix(i,j,k))/chin1(i,j,k) - gxx(i,j,k) * gxxx(i,j,k) )*HALF
        Gamyxx(i,j,k) = Gamyxx(i,j,k) - (                     - gxx(i,j,k) * gxxy(i,j,k) )*HALF
        Gamzxx(i,j,k) = Gamzxx(i,j,k) - (                     - gxx(i,j,k) * gxxz(i,j,k) )*HALF
        Gamxyy(i,j,k) = Gamxyy(i,j,k) - (                     - gyy(i,j,k) * gxxx(i,j,k) )*HALF
        Gamyyy(i,j,k) = Gamyyy(i,j,k) - ( (chiy(i,j,k) + chiy(i,j,k))/chin1(i,j,k) - gyy(i,j,k) * gxxy(i,j,k) )*HALF
        Gamzyy(i,j,k) = Gamzyy(i,j,k) - (                     - gyy(i,j,k) * gxxz(i,j,k) )*HALF
        Gamxzz(i,j,k) = Gamxzz(i,j,k) - (                     - gzz(i,j,k) * gxxx(i,j,k) )*HALF
        Gamyzz(i,j,k) = Gamyzz(i,j,k) - (                     - gzz(i,j,k) * gxxy(i,j,k) )*HALF
        Gamzzz(i,j,k) = Gamzzz(i,j,k) - ( (chiz(i,j,k) + chiz(i,j,k))/chin1(i,j,k) - gzz(i,j,k) * gxxz(i,j,k) )*HALF
        Gamxxy(i,j,k) = Gamxxy(i,j,k) - (  chiy(i,j,k)        /chin1(i,j,k) - gxy(i,j,k) * gxxx(i,j,k) )*HALF
        Gamyxy(i,j,k) = Gamyxy(i,j,k) - (         chix(i,j,k) /chin1(i,j,k) - gxy(i,j,k) * gxxy(i,j,k) )*HALF
        Gamzxy(i,j,k) = Gamzxy(i,j,k) - (                     - gxy(i,j,k) * gxxz(i,j,k) )*HALF
        Gamxxz(i,j,k) = Gamxxz(i,j,k) - (  chiz(i,j,k)        /chin1(i,j,k) - gxz(i,j,k) * gxxx(i,j,k) )*HALF
        Gamyxz(i,j,k) = Gamyxz(i,j,k) - (                     - gxz(i,j,k) * gxxy(i,j,k) )*HALF
        Gamzxz(i,j,k) = Gamzxz(i,j,k) - (         chix(i,j,k) /chin1(i,j,k) - gxz(i,j,k) * gxxz(i,j,k) )*HALF
        Gamxyz(i,j,k) = Gamxyz(i,j,k) - (                     - gyz(i,j,k) * gxxx(i,j,k) )*HALF
        Gamyyz(i,j,k) = Gamyyz(i,j,k) - (  chiz(i,j,k)        /chin1(i,j,k) - gyz(i,j,k) * gxxy(i,j,k) )*HALF
        Gamzyz(i,j,k) = Gamzyz(i,j,k) - (         chiy(i,j,k) /chin1(i,j,k) - gyz(i,j,k) * gxxz(i,j,k) )*HALF
      end do
    end do
  end do
!$OMP END PARALLEL DO

! covariant second derivatives of the lapse respect to physical metric

   call fdderivs(ex,Lap,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z, &
                SYM,SYM,SYM,symmetry,Lev)

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        fxx(i,j,k) = fxx(i,j,k) - Gamxxx(i,j,k)*Lapx(i,j,k) - Gamyxx(i,j,k)*Lapy(i,j,k) - Gamzxx(i,j,k)*Lapz(i,j,k)
        fyy(i,j,k) = fyy(i,j,k) - Gamxyy(i,j,k)*Lapx(i,j,k) - Gamyyy(i,j,k)*Lapy(i,j,k) - Gamzyy(i,j,k)*Lapz(i,j,k)
        fzz(i,j,k) = fzz(i,j,k) - Gamxzz(i,j,k)*Lapx(i,j,k) - Gamyzz(i,j,k)*Lapy(i,j,k) - Gamzzz(i,j,k)*Lapz(i,j,k)
        fxy(i,j,k) = fxy(i,j,k) - Gamxxy(i,j,k)*Lapx(i,j,k) - Gamyxy(i,j,k)*Lapy(i,j,k) - Gamzxy(i,j,k)*Lapz(i,j,k)
        fxz(i,j,k) = fxz(i,j,k) - Gamxxz(i,j,k)*Lapx(i,j,k) - Gamyxz(i,j,k)*Lapy(i,j,k) - Gamzxz(i,j,k)*Lapz(i,j,k)
        fyz(i,j,k) = fyz(i,j,k) - Gamxyz(i,j,k)*Lapx(i,j,k) - Gamyyz(i,j,k)*Lapy(i,j,k) - Gamzyz(i,j,k)*Lapz(i,j,k)

! store D^i D_i Lap in trK_rhs upto chi
        trK_rhs(i,j,k) =    gupxx(i,j,k) * fxx(i,j,k) + gupyy(i,j,k) * fyy(i,j,k) + gupzz(i,j,k) * fzz(i,j,k) + &
              TWO* ( gupxy(i,j,k) * fxy(i,j,k) + gupxz(i,j,k) * fxz(i,j,k) + gupyz(i,j,k) * fyz(i,j,k) )
! Add lapse and S_ij parts to Ricci tensor:

        fxx(i,j,k) = EIGHT * PI * alpn1(i,j,k) * Sxx(i,j,k) + fxx(i,j,k)
        fxy(i,j,k) = EIGHT * PI * alpn1(i,j,k) * Sxy(i,j,k) + fxy(i,j,k)
        fxz(i,j,k) = EIGHT * PI * alpn1(i,j,k) * Sxz(i,j,k) + fxz(i,j,k)
        fyy(i,j,k) = EIGHT * PI * alpn1(i,j,k) * Syy(i,j,k) + fyy(i,j,k)
        fyz(i,j,k) = EIGHT * PI * alpn1(i,j,k) * Syz(i,j,k) + fyz(i,j,k)
        fzz(i,j,k) = EIGHT * PI * alpn1(i,j,k) * Szz(i,j,k) + fzz(i,j,k)

! Compute trace-free part (note: chi^-1 and chi cancel!):
        f(i,j,k) =        gupxx(i,j,k) * fxx(i,j,k) + gupyy(i,j,k) * fyy(i,j,k) + gupzz(i,j,k) * fzz(i,j,k) + &
            TWO* ( gupxy(i,j,k) * fxy(i,j,k) + gupxz(i,j,k) * fxz(i,j,k) + gupyz(i,j,k) * fyz(i,j,k) )

        f(i,j,k) = F1o3 * (Hcon(i,j,k)*alpn1(i,j,k) - f(i,j,k))

        fxx(i,j,k) = alpn1(i,j,k) * Rxx(i,j,k) - fxx(i,j,k)
        fxy(i,j,k) = alpn1(i,j,k) * Rxy(i,j,k) - fxy(i,j,k)
        fxz(i,j,k) = alpn1(i,j,k) * Rxz(i,j,k) - fxz(i,j,k)
        fyy(i,j,k) = alpn1(i,j,k) * Ryy(i,j,k) - fyy(i,j,k)
        fyz(i,j,k) = alpn1(i,j,k) * Ryz(i,j,k) - fyz(i,j,k)
        fzz(i,j,k) = alpn1(i,j,k) * Rzz(i,j,k) - fzz(i,j,k)

        Axx_rhs(i,j,k) = fxx(i,j,k) - gxx(i,j,k) * f(i,j,k)
        Ayy_rhs(i,j,k) = fyy(i,j,k) - gyy(i,j,k) * f(i,j,k)
        Azz_rhs(i,j,k) = fzz(i,j,k) - gzz(i,j,k) * f(i,j,k)
        Axy_rhs(i,j,k) = fxy(i,j,k) - gxy(i,j,k) * f(i,j,k)
        Axz_rhs(i,j,k) = fxz(i,j,k) - gxz(i,j,k) * f(i,j,k)
        Ayz_rhs(i,j,k) = fyz(i,j,k) - gyz(i,j,k) * f(i,j,k)

! Now: store A_il A^l_j into fij:

        fxx(i,j,k) =       gupxx(i,j,k) * Axx(i,j,k) * Axx(i,j,k) + gupyy(i,j,k) * Axy(i,j,k) * Axy(i,j,k) + gupzz(i,j,k) * Axz(i,j,k) * Axz(i,j,k) + &
             TWO * (gupxy(i,j,k) * Axx(i,j,k) * Axy(i,j,k) + gupxz(i,j,k) * Axx(i,j,k) * Axz(i,j,k) + gupyz(i,j,k) * Axy(i,j,k) * Axz(i,j,k))
        fyy(i,j,k) =       gupxx(i,j,k) * Axy(i,j,k) * Axy(i,j,k) + gupyy(i,j,k) * Ayy(i,j,k) * Ayy(i,j,k) + gupzz(i,j,k) * Ayz(i,j,k) * Ayz(i,j,k) + &
             TWO * (gupxy(i,j,k) * Axy(i,j,k) * Ayy(i,j,k) + gupxz(i,j,k) * Axy(i,j,k) * Ayz(i,j,k) + gupyz(i,j,k) * Ayy(i,j,k) * Ayz(i,j,k))
        fzz(i,j,k) =       gupxx(i,j,k) * Axz(i,j,k) * Axz(i,j,k) + gupyy(i,j,k) * Ayz(i,j,k) * Ayz(i,j,k) + gupzz(i,j,k) * Azz(i,j,k) * Azz(i,j,k) + &
             TWO * (gupxy(i,j,k) * Axz(i,j,k) * Ayz(i,j,k) + gupxz(i,j,k) * Axz(i,j,k) * Azz(i,j,k) + gupyz(i,j,k) * Ayz(i,j,k) * Azz(i,j,k))
        fxy(i,j,k) =       gupxx(i,j,k) * Axx(i,j,k) * Axy(i,j,k) + gupyy(i,j,k) * Axy(i,j,k) * Ayy(i,j,k) + gupzz(i,j,k) * Axz(i,j,k) * Ayz(i,j,k) + &
                    gupxy(i,j,k) *(Axx(i,j,k) * Ayy(i,j,k) + Axy(i,j,k) * Axy(i,j,k))                            + &
                    gupxz(i,j,k) *(Axx(i,j,k) * Ayz(i,j,k) + Axz(i,j,k) * Axy(i,j,k))                            + &
                    gupyz(i,j,k) *(Axy(i,j,k) * Ayz(i,j,k) + Axz(i,j,k) * Ayy(i,j,k))
        fxz(i,j,k) =       gupxx(i,j,k) * Axx(i,j,k) * Axz(i,j,k) + gupyy(i,j,k) * Axy(i,j,k) * Ayz(i,j,k) + gupzz(i,j,k) * Axz(i,j,k) * Azz(i,j,k) + &
                    gupxy(i,j,k) *(Axx(i,j,k) * Ayz(i,j,k) + Axy(i,j,k) * Axz(i,j,k))                            + &
                    gupxz(i,j,k) *(Axx(i,j,k) * Azz(i,j,k) + Axz(i,j,k) * Axz(i,j,k))                            + &
                    gupyz(i,j,k) *(Axy(i,j,k) * Azz(i,j,k) + Axz(i,j,k) * Ayz(i,j,k))
        fyz(i,j,k) =       gupxx(i,j,k) * Axy(i,j,k) * Axz(i,j,k) + gupyy(i,j,k) * Ayy(i,j,k) * Ayz(i,j,k) + gupzz(i,j,k) * Ayz(i,j,k) * Azz(i,j,k) + &
                    gupxy(i,j,k) *(Axy(i,j,k) * Ayz(i,j,k) + Ayy(i,j,k) * Axz(i,j,k))                            + &
                    gupxz(i,j,k) *(Axy(i,j,k) * Azz(i,j,k) + Ayz(i,j,k) * Axz(i,j,k))                            + &
                    gupyz(i,j,k) *(Ayy(i,j,k) * Azz(i,j,k) + Ayz(i,j,k) * Ayz(i,j,k))

        f(i,j,k) = chin1(i,j,k)
! store D^i D_i Lap in trK_rhs
        trK_rhs(i,j,k) = f(i,j,k)*trK_rhs(i,j,k)

        Axx_rhs(i,j,k) =           f(i,j,k) * Axx_rhs(i,j,k)+ alpn1(i,j,k) * (trKd(i,j,k) * Axx(i,j,k) - TWO * fxx(i,j,k))  + &
                 TWO * (  Axx(i,j,k) * betaxx(i,j,k) +   Axy(i,j,k) * betayx(i,j,k) +   Axz(i,j,k) * betazx(i,j,k) ) - &
                   F2o3 * Axx(i,j,k) * div_beta(i,j,k)

        Ayy_rhs(i,j,k) =           f(i,j,k) * Ayy_rhs(i,j,k)+ alpn1(i,j,k) * (trKd(i,j,k) * Ayy(i,j,k) - TWO * fyy(i,j,k))  + &
                 TWO * (  Axy(i,j,k) * betaxy(i,j,k) +   Ayy(i,j,k) * betayy(i,j,k) +   Ayz(i,j,k) * betazy(i,j,k) ) - &
                   F2o3 * Ayy(i,j,k) * div_beta(i,j,k)

        Azz_rhs(i,j,k) =           f(i,j,k) * Azz_rhs(i,j,k)+ alpn1(i,j,k) * (trKd(i,j,k) * Azz(i,j,k) - TWO * fzz(i,j,k))  + &
                 TWO * (  Axz(i,j,k) * betaxz(i,j,k) +   Ayz(i,j,k) * betayz(i,j,k) +   Azz(i,j,k) * betazz(i,j,k) ) - &
                   F2o3 * Azz(i,j,k) * div_beta(i,j,k)

        Axy_rhs(i,j,k) =           f(i,j,k) * Axy_rhs(i,j,k)+ alpn1(i,j,k) *( trKd(i,j,k) * Axy(i,j,k)  - TWO * fxy(i,j,k) )+ &
                          Axx(i,j,k) * betaxy(i,j,k)                  +   Axz(i,j,k) * betazy(i,j,k)   + &
                                           Ayy(i,j,k) * betayx(i,j,k) +   Ayz(i,j,k) * betazx(i,j,k)   + &
                   F1o3 * Axy(i,j,k) * div_beta(i,j,k)                -   Axy(i,j,k) * betazz(i,j,k)

        Ayz_rhs(i,j,k) =           f(i,j,k) * Ayz_rhs(i,j,k)+ alpn1(i,j,k) *( trKd(i,j,k) * Ayz(i,j,k)  - TWO * fyz(i,j,k) )+ &
                          Axy(i,j,k) * betaxz(i,j,k) +   Ayy(i,j,k) * betayz(i,j,k)                    + &
                          Axz(i,j,k) * betaxy(i,j,k)                  +   Azz(i,j,k) * betazy(i,j,k)   + &
                   F1o3 * Ayz(i,j,k) * div_beta(i,j,k)                -   Ayz(i,j,k) * betaxx(i,j,k)

        Axz_rhs(i,j,k) =           f(i,j,k) * Axz_rhs(i,j,k)+ alpn1(i,j,k) *( trKd(i,j,k) * Axz(i,j,k)  - TWO * fxz(i,j,k) )+ &
                          Axx(i,j,k) * betaxz(i,j,k) +   Axy(i,j,k) * betayz(i,j,k)                    + &
                                           Ayz(i,j,k) * betayx(i,j,k) +   Azz(i,j,k) * betazx(i,j,k)   + &
                   F1o3 * Axz(i,j,k) * div_beta(i,j,k)                -   Axz(i,j,k) * betayy(i,j,k)      !rhs for Aij

! Compute trace of S_ij

        S(i,j,k) =  f(i,j,k) * ( gupxx(i,j,k) * Sxx(i,j,k) + gupyy(i,j,k) * Syy(i,j,k) + gupzz(i,j,k) * Szz(i,j,k) + &
           TWO * ( gupxy(i,j,k) * Sxy(i,j,k) + gupxz(i,j,k) * Sxz(i,j,k) + gupyz(i,j,k) * Syz(i,j,k) ) )

        trK_rhs(i,j,k) = - trK_rhs(i,j,k) + alpn1(i,j,k) *( F1o3 * trKd(i,j,k) * trKd(i,j,k) + &
                      gupxx(i,j,k) * fxx(i,j,k) + gupyy(i,j,k) * fyy(i,j,k) + gupzz(i,j,k) * fzz(i,j,k)   + &
              TWO * ( gupxy(i,j,k) * fxy(i,j,k) + gupxz(i,j,k) * fxz(i,j,k) + gupyz(i,j,k) * fyz(i,j,k) ) + &
             FOUR * PI * ( rho(i,j,k) + S(i,j,k) ))                                !rhs for trK

!!!!!gauge variable part
        Lap_rhs(i,j,k) = -TWO*alpn1(i,j,k)*trK(i,j,k)
      end do
    end do
  end do
!$OMP END PARALLEL DO

#if (GAUGE == 0)
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        betax_rhs(i,j,k) = FF*dtSfx(i,j,k)
        betay_rhs(i,j,k) = FF*dtSfy(i,j,k)
        betaz_rhs(i,j,k) = FF*dtSfz(i,j,k)

        dtSfx_rhs(i,j,k) = Gamx_rhs(i,j,k) - eta*dtSfx(i,j,k)
        dtSfy_rhs(i,j,k) = Gamy_rhs(i,j,k) - eta*dtSfy(i,j,k)
        dtSfz_rhs(i,j,k) = Gamz_rhs(i,j,k) - eta*dtSfz(i,j,k)
      end do
    end do
  end do
!$OMP END PARALLEL DO
#elif  (GAUGE == 1)
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        betax_rhs(i,j,k) = Gamx(i,j,k) - eta*betax(i,j,k)
        betay_rhs(i,j,k) = Gamy(i,j,k) - eta*betay(i,j,k)
        betaz_rhs(i,j,k) = Gamz(i,j,k) - eta*betaz(i,j,k)

        dtSfx_rhs(i,j,k) = ZEO
        dtSfy_rhs(i,j,k) = ZEO
        dtSfz_rhs(i,j,k) = ZEO
      end do
    end do
  end do
!$OMP END PARALLEL DO
#endif
!!!!!Z4 part
! H = trR + 2/3 * trKd^2 - A_ij * A^ij - 16 * PI * rho
! here trR is respect to physical metric

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        Hcon(i,j,k) = chin1(i,j,k)*Hcon(i,j,k) + F2o3 * trKd(i,j,k) * trKd(i,j,k) -(&
             gupxx(i,j,k) * ( &
             gupxx(i,j,k) * Axx(i,j,k) * Axx(i,j,k) + gupyy(i,j,k) * Axy(i,j,k) * Axy(i,j,k) + gupzz(i,j,k) * Axz(i,j,k) * Axz(i,j,k) + &
             TWO * (gupxy(i,j,k) * Axx(i,j,k) * Axy(i,j,k) + gupxz(i,j,k) * Axx(i,j,k) * Axz(i,j,k) + gupyz(i,j,k) * Axy(i,j,k) * Axz(i,j,k)) ) + &
             gupyy(i,j,k) * ( &
             gupxx(i,j,k) * Axy(i,j,k) * Axy(i,j,k) + gupyy(i,j,k) * Ayy(i,j,k) * Ayy(i,j,k) + gupzz(i,j,k) * Ayz(i,j,k) * Ayz(i,j,k) + &
             TWO * (gupxy(i,j,k) * Axy(i,j,k) * Ayy(i,j,k) + gupxz(i,j,k) * Axy(i,j,k) * Ayz(i,j,k) + gupyz(i,j,k) * Ayy(i,j,k) * Ayz(i,j,k)) ) + &
             gupzz(i,j,k) * ( &
             gupxx(i,j,k) * Axz(i,j,k) * Axz(i,j,k) + gupyy(i,j,k) * Ayz(i,j,k) * Ayz(i,j,k) + gupzz(i,j,k) * Azz(i,j,k) * Azz(i,j,k) + &
             TWO * (gupxy(i,j,k) * Axz(i,j,k) * Ayz(i,j,k) + gupxz(i,j,k) * Axz(i,j,k) * Azz(i,j,k) + gupyz(i,j,k) * Ayz(i,j,k) * Azz(i,j,k)) ) + &
             TWO * ( &
             gupxy(i,j,k) * ( &
             gupxx(i,j,k) * Axx(i,j,k) * Axy(i,j,k) + gupyy(i,j,k) * Axy(i,j,k) * Ayy(i,j,k) + gupzz(i,j,k) * Axz(i,j,k) * Ayz(i,j,k) + &
             gupxy(i,j,k) * (Axx(i,j,k) * Ayy(i,j,k) + Axy(i,j,k) * Axy(i,j,k)) + &
             gupxz(i,j,k) * (Axx(i,j,k) * Ayz(i,j,k) + Axz(i,j,k) * Axy(i,j,k)) + &
             gupyz(i,j,k) * (Axy(i,j,k) * Ayz(i,j,k) + Axz(i,j,k) * Ayy(i,j,k)) ) + &
             gupxz(i,j,k) * ( &
             gupxx(i,j,k) * Axx(i,j,k) * Axz(i,j,k) + gupyy(i,j,k) * Axy(i,j,k) * Ayz(i,j,k) + gupzz(i,j,k) * Axz(i,j,k) * Azz(i,j,k) + &
             gupxy(i,j,k) * (Axx(i,j,k) * Ayz(i,j,k) + Axy(i,j,k) * Axz(i,j,k)) + &
             gupxz(i,j,k) * (Axx(i,j,k) * Azz(i,j,k) + Axz(i,j,k) * Axz(i,j,k)) + &
             gupyz(i,j,k) * (Axy(i,j,k) * Azz(i,j,k) + Axz(i,j,k) * Ayz(i,j,k)) ) + &
             gupyz(i,j,k) * ( &
             gupxx(i,j,k) * Axy(i,j,k) * Axz(i,j,k) + gupyy(i,j,k) * Ayy(i,j,k) * Ayz(i,j,k) + gupzz(i,j,k) * Ayz(i,j,k) * Azz(i,j,k) + &
             gupxy(i,j,k) * (Axy(i,j,k) * Ayz(i,j,k) + Ayy(i,j,k) * Axz(i,j,k)) + &
             gupxz(i,j,k) * (Axy(i,j,k) * Azz(i,j,k) + Ayz(i,j,k) * Axz(i,j,k)) + &
             gupyz(i,j,k) * (Ayy(i,j,k) * Azz(i,j,k) + Ayz(i,j,k) * Ayz(i,j,k)) ) ))- F16 * PI * rho(i,j,k)
      end do
    end do
  end do
!$OMP END PARALLEL DO
! M_j = gupki*(-1/chi d_k chi*A_ij + D_k A_ij) - 2/3 d_j trK - 8 PI s_j where D respect to physical metric
! store D_i A_jk - 1/chi d_i chi*A_jk in gjk_i
  call fderivs(ex,Axx,gxxx,gxxy,gxxz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,lev)
  call fderivs(ex,Axy,gxyx,gxyy,gxyz,X,Y,Z,ANTI,ANTI,SYM ,Symmetry,lev)
  call fderivs(ex,Axz,gxzx,gxzy,gxzz,X,Y,Z,ANTI,SYM ,ANTI,Symmetry,lev)
  call fderivs(ex,Ayy,gyyx,gyyy,gyyz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,lev)
  call fderivs(ex,Ayz,gyzx,gyzy,gyzz,X,Y,Z,SYM ,ANTI,ANTI,Symmetry,lev)
  call fderivs(ex,Azz,gzzx,gzzy,gzzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,lev)

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        gxxx(i,j,k) = gxxx(i,j,k) - (  Gamxxx(i,j,k) * Axx(i,j,k) + Gamyxx(i,j,k) * Axy(i,j,k) + Gamzxx(i,j,k) * Axz(i,j,k) &
                       + Gamxxx(i,j,k) * Axx(i,j,k) + Gamyxx(i,j,k) * Axy(i,j,k) + Gamzxx(i,j,k) * Axz(i,j,k)) - chix(i,j,k)*Axx(i,j,k)/chin1(i,j,k)
        gxyx(i,j,k) = gxyx(i,j,k) - (  Gamxxy(i,j,k) * Axx(i,j,k) + Gamyxy(i,j,k) * Axy(i,j,k) + Gamzxy(i,j,k) * Axz(i,j,k) &
                       + Gamxxx(i,j,k) * Axy(i,j,k) + Gamyxx(i,j,k) * Ayy(i,j,k) + Gamzxx(i,j,k) * Ayz(i,j,k)) - chix(i,j,k)*Axy(i,j,k)/chin1(i,j,k)
        gxzx(i,j,k) = gxzx(i,j,k) - (  Gamxxz(i,j,k) * Axx(i,j,k) + Gamyxz(i,j,k) * Axy(i,j,k) + Gamzxz(i,j,k) * Axz(i,j,k) &
                       + Gamxxx(i,j,k) * Axz(i,j,k) + Gamyxx(i,j,k) * Ayz(i,j,k) + Gamzxx(i,j,k) * Azz(i,j,k)) - chix(i,j,k)*Axz(i,j,k)/chin1(i,j,k)
        gyyx(i,j,k) = gyyx(i,j,k) - (  Gamxxy(i,j,k) * Axy(i,j,k) + Gamyxy(i,j,k) * Ayy(i,j,k) + Gamzxy(i,j,k) * Ayz(i,j,k) &
                       + Gamxxy(i,j,k) * Axy(i,j,k) + Gamyxy(i,j,k) * Ayy(i,j,k) + Gamzxy(i,j,k) * Ayz(i,j,k)) - chix(i,j,k)*Ayy(i,j,k)/chin1(i,j,k)
        gyzx(i,j,k) = gyzx(i,j,k) - (  Gamxxz(i,j,k) * Axy(i,j,k) + Gamyxz(i,j,k) * Ayy(i,j,k) + Gamzxz(i,j,k) * Ayz(i,j,k) &
                       + Gamxxy(i,j,k) * Axz(i,j,k) + Gamyxy(i,j,k) * Ayz(i,j,k) + Gamzxy(i,j,k) * Azz(i,j,k)) - chix(i,j,k)*Ayz(i,j,k)/chin1(i,j,k)
        gzzx(i,j,k) = gzzx(i,j,k) - (  Gamxxz(i,j,k) * Axz(i,j,k) + Gamyxz(i,j,k) * Ayz(i,j,k) + Gamzxz(i,j,k) * Azz(i,j,k) &
                       + Gamxxz(i,j,k) * Axz(i,j,k) + Gamyxz(i,j,k) * Ayz(i,j,k) + Gamzxz(i,j,k) * Azz(i,j,k)) - chix(i,j,k)*Azz(i,j,k)/chin1(i,j,k)
        gxxy(i,j,k) = gxxy(i,j,k) - (  Gamxxy(i,j,k) * Axx(i,j,k) + Gamyxy(i,j,k) * Axy(i,j,k) + Gamzxy(i,j,k) * Axz(i,j,k) &
                       + Gamxxy(i,j,k) * Axx(i,j,k) + Gamyxy(i,j,k) * Axy(i,j,k) + Gamzxy(i,j,k) * Axz(i,j,k)) - chiy(i,j,k)*Axx(i,j,k)/chin1(i,j,k)
        gxyy(i,j,k) = gxyy(i,j,k) - (  Gamxyy(i,j,k) * Axx(i,j,k) + Gamyyy(i,j,k) * Axy(i,j,k) + Gamzyy(i,j,k) * Axz(i,j,k) &
                       + Gamxxy(i,j,k) * Axy(i,j,k) + Gamyxy(i,j,k) * Ayy(i,j,k) + Gamzxy(i,j,k) * Ayz(i,j,k)) - chiy(i,j,k)*Axy(i,j,k)/chin1(i,j,k)
        gxzy(i,j,k) = gxzy(i,j,k) - (  Gamxyz(i,j,k) * Axx(i,j,k) + Gamyyz(i,j,k) * Axy(i,j,k) + Gamzyz(i,j,k) * Axz(i,j,k) &
                       + Gamxxy(i,j,k) * Axz(i,j,k) + Gamyxy(i,j,k) * Ayz(i,j,k) + Gamzxy(i,j,k) * Azz(i,j,k)) - chiy(i,j,k)*Axz(i,j,k)/chin1(i,j,k)
        gyyy(i,j,k) = gyyy(i,j,k) - (  Gamxyy(i,j,k) * Axy(i,j,k) + Gamyyy(i,j,k) * Ayy(i,j,k) + Gamzyy(i,j,k) * Ayz(i,j,k) &
                       + Gamxyy(i,j,k) * Axy(i,j,k) + Gamyyy(i,j,k) * Ayy(i,j,k) + Gamzyy(i,j,k) * Ayz(i,j,k)) - chiy(i,j,k)*Ayy(i,j,k)/chin1(i,j,k)
        gyzy(i,j,k) = gyzy(i,j,k) - (  Gamxyz(i,j,k) * Axy(i,j,k) + Gamyyz(i,j,k) * Ayy(i,j,k) + Gamzyz(i,j,k) * Ayz(i,j,k) &
                       + Gamxyy(i,j,k) * Axz(i,j,k) + Gamyyy(i,j,k) * Ayz(i,j,k) + Gamzyy(i,j,k) * Azz(i,j,k)) - chiy(i,j,k)*Ayz(i,j,k)/chin1(i,j,k)
        gzzy(i,j,k) = gzzy(i,j,k) - (  Gamxyz(i,j,k) * Axz(i,j,k) + Gamyyz(i,j,k) * Ayz(i,j,k) + Gamzyz(i,j,k) * Azz(i,j,k) &
                       + Gamxyz(i,j,k) * Axz(i,j,k) + Gamyyz(i,j,k) * Ayz(i,j,k) + Gamzyz(i,j,k) * Azz(i,j,k)) - chiy(i,j,k)*Azz(i,j,k)/chin1(i,j,k)
        gxxz(i,j,k) = gxxz(i,j,k) - (  Gamxxz(i,j,k) * Axx(i,j,k) + Gamyxz(i,j,k) * Axy(i,j,k) + Gamzxz(i,j,k) * Axz(i,j,k) &
                       + Gamxxz(i,j,k) * Axx(i,j,k) + Gamyxz(i,j,k) * Axy(i,j,k) + Gamzxz(i,j,k) * Axz(i,j,k)) - chiz(i,j,k)*Axx(i,j,k)/chin1(i,j,k)
        gxyz(i,j,k) = gxyz(i,j,k) - (  Gamxyz(i,j,k) * Axx(i,j,k) + Gamyyz(i,j,k) * Axy(i,j,k) + Gamzyz(i,j,k) * Axz(i,j,k) &
                       + Gamxxz(i,j,k) * Axy(i,j,k) + Gamyxz(i,j,k) * Ayy(i,j,k) + Gamzxz(i,j,k) * Ayz(i,j,k)) - chiz(i,j,k)*Axy(i,j,k)/chin1(i,j,k)
        gxzz(i,j,k) = gxzz(i,j,k) - (  Gamxzz(i,j,k) * Axx(i,j,k) + Gamyzz(i,j,k) * Axy(i,j,k) + Gamzzz(i,j,k) * Axz(i,j,k) &
                       + Gamxxz(i,j,k) * Axz(i,j,k) + Gamyxz(i,j,k) * Ayz(i,j,k) + Gamzxz(i,j,k) * Azz(i,j,k)) - chiz(i,j,k)*Axz(i,j,k)/chin1(i,j,k)
        gyyz(i,j,k) = gyyz(i,j,k) - (  Gamxyz(i,j,k) * Axy(i,j,k) + Gamyyz(i,j,k) * Ayy(i,j,k) + Gamzyz(i,j,k) * Ayz(i,j,k) &
                       + Gamxyz(i,j,k) * Axy(i,j,k) + Gamyyz(i,j,k) * Ayy(i,j,k) + Gamzyz(i,j,k) * Ayz(i,j,k)) - chiz(i,j,k)*Ayy(i,j,k)/chin1(i,j,k)
        gyzz(i,j,k) = gyzz(i,j,k) - (  Gamxzz(i,j,k) * Axy(i,j,k) + Gamyzz(i,j,k) * Ayy(i,j,k) + Gamzzz(i,j,k) * Ayz(i,j,k) &
                       + Gamxyz(i,j,k) * Axz(i,j,k) + Gamyyz(i,j,k) * Ayz(i,j,k) + Gamzyz(i,j,k) * Azz(i,j,k)) - chiz(i,j,k)*Ayz(i,j,k)/chin1(i,j,k)
        gzzz(i,j,k) = gzzz(i,j,k) - (  Gamxzz(i,j,k) * Axz(i,j,k) + Gamyzz(i,j,k) * Ayz(i,j,k) + Gamzzz(i,j,k) * Azz(i,j,k) &
                       + Gamxzz(i,j,k) * Axz(i,j,k) + Gamyzz(i,j,k) * Ayz(i,j,k) + Gamzzz(i,j,k) * Azz(i,j,k)) - chiz(i,j,k)*Azz(i,j,k)/chin1(i,j,k)
        Mxcon(i,j,k)  = gupxx(i,j,k)*gxxx(i,j,k) + gupyy(i,j,k)*gxyy(i,j,k) + gupzz(i,j,k)*gxzz(i,j,k) &
                +gupxy(i,j,k)*gxyx(i,j,k) + gupxz(i,j,k)*gxzx(i,j,k) + gupyz(i,j,k)*gxzy(i,j,k) &
                +gupxy(i,j,k)*gxxy(i,j,k) + gupxz(i,j,k)*gxxz(i,j,k) + gupyz(i,j,k)*gxyz(i,j,k)
        Mycon(i,j,k)  = gupxx(i,j,k)*gxyx(i,j,k) + gupyy(i,j,k)*gyyy(i,j,k) + gupzz(i,j,k)*gyzz(i,j,k) &
                +gupxy(i,j,k)*gyyx(i,j,k) + gupxz(i,j,k)*gyzx(i,j,k) + gupyz(i,j,k)*gyzy(i,j,k) &
                +gupxy(i,j,k)*gxyy(i,j,k) + gupxz(i,j,k)*gxyz(i,j,k) + gupyz(i,j,k)*gyyz(i,j,k)
        Mzcon(i,j,k)  = gupxx(i,j,k)*gxzx(i,j,k) + gupyy(i,j,k)*gyzy(i,j,k) + gupzz(i,j,k)*gzzz(i,j,k) &
                +gupxy(i,j,k)*gyzx(i,j,k) + gupxz(i,j,k)*gzzx(i,j,k) + gupyz(i,j,k)*gzzy(i,j,k) &
                +gupxy(i,j,k)*gxzy(i,j,k) + gupxz(i,j,k)*gxzz(i,j,k) + gupyz(i,j,k)*gyzz(i,j,k)
! we have already considered TZ_,i in K_,i here, or to say here Micon =
! Micon+TZ_,i indeed
        Mxcon(i,j,k)  = Mxcon(i,j,k) - F2o3*Kx(i,j,k) - F8*PI*sx(i,j,k)
        Mycon(i,j,k)  = Mycon(i,j,k) - F2o3*Ky(i,j,k) - F8*PI*sy(i,j,k)
        Mzcon(i,j,k)  = Mzcon(i,j,k) - F2o3*Kz(i,j,k) - F8*PI*sz(i,j,k)

        f(i,j,k) = TZ_rhs(i,j,k)

        TZ_rhs(i,j,k) = alpn1(i,j,k)*Hcon(i,j,k)/TWO
      end do
    end do
  end do
!$OMP END PARALLEL DO
! delete TWO*Z^i_,i From Hcon' to get Hcon, this is wrong
!  Hcon = Hcon - f

  SSS(1)=SYM
  SSS(2)=SYM
  SSS(3)=SYM

  AAS(1)=ANTI
  AAS(2)=ANTI
  AAS(3)=SYM

  ASA(1)=ANTI
  ASA(2)=SYM
  ASA(3)=ANTI

  SAA(1)=SYM
  SAA(2)=ANTI
  SAA(3)=ANTI

  ASS(1)=ANTI
  ASS(2)=SYM
  ASS(3)=SYM

  SAS(1)=SYM
  SAS(2)=ANTI
  SAS(3)=SYM

  SSA(1)=SYM
  SSA(2)=SYM
  SSA(3)=ANTI
!!!!!!!!!advection term part
  call lopsided(ex,X,Y,Z,gxx,gxx_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,gxy,gxy_rhs,betax,betay,betaz,Symmetry,AAS)
  call lopsided(ex,X,Y,Z,gxz,gxz_rhs,betax,betay,betaz,Symmetry,ASA)
  call lopsided(ex,X,Y,Z,gyy,gyy_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,gyz,gyz_rhs,betax,betay,betaz,Symmetry,SAA)
  call lopsided(ex,X,Y,Z,gzz,gzz_rhs,betax,betay,betaz,Symmetry,SSS)

  call lopsided(ex,X,Y,Z,Axx,Axx_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,Axy,Axy_rhs,betax,betay,betaz,Symmetry,AAS)
  call lopsided(ex,X,Y,Z,Axz,Axz_rhs,betax,betay,betaz,Symmetry,ASA)
  call lopsided(ex,X,Y,Z,Ayy,Ayy_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,Ayz,Ayz_rhs,betax,betay,betaz,Symmetry,SAA)
  call lopsided(ex,X,Y,Z,Azz,Azz_rhs,betax,betay,betaz,Symmetry,SSS)

  call lopsided(ex,X,Y,Z,chi,chi_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,trK,trK_rhs,betax,betay,betaz,Symmetry,SSS)

  call lopsided(ex,X,Y,Z,Gamx,Gamx_rhs,betax,betay,betaz,Symmetry,ASS)
  call lopsided(ex,X,Y,Z,Gamy,Gamy_rhs,betax,betay,betaz,Symmetry,SAS)
  call lopsided(ex,X,Y,Z,Gamz,Gamz_rhs,betax,betay,betaz,Symmetry,SSA)

  call lopsided(ex,X,Y,Z,Lap,Lap_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,betax,betax_rhs,betax,betay,betaz,Symmetry,ASS)
  call lopsided(ex,X,Y,Z,betay,betay_rhs,betax,betay,betaz,Symmetry,SAS)
  call lopsided(ex,X,Y,Z,betaz,betaz_rhs,betax,betay,betaz,Symmetry,SSA)
  call lopsided(ex,X,Y,Z,dtSfx,dtSfx_rhs,betax,betay,betaz,Symmetry,ASS)
  call lopsided(ex,X,Y,Z,dtSfy,dtSfy_rhs,betax,betay,betaz,Symmetry,SAS)
  call lopsided(ex,X,Y,Z,dtSfz,dtSfz_rhs,betax,betay,betaz,Symmetry,SSA)

  call lopsided(ex,X,Y,Z,TZ,TZ_rhs,betax,betay,betaz,Symmetry,SSS)

! constraint damping terms
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
    do k = 1, ex(3)
      do j = 1, ex(2)
        do i = 1, ex(1)
          TZ_rhs(i,j,k) = TZ_rhs(i,j,k) - alpn1(i,j,k)*(TWO+kappa2)*kappa1*TZ(i,j,k)
          trK_rhs(i,j,k) = trK_rhs(i,j,k) + alpn1(i,j,k)*kappa1*(ONE-kappa2)*TZ(i,j,k)
          Gamx_rhs(i,j,k) = Gamx_rhs(i,j,k) - TWO*alpn1(i,j,k)*kappa1*(Gamx(i,j,k)-Gamxa(i,j,k))
          Gamy_rhs(i,j,k) = Gamy_rhs(i,j,k) - TWO*alpn1(i,j,k)*kappa1*(Gamy(i,j,k)-Gamya(i,j,k))
          Gamz_rhs(i,j,k) = Gamz_rhs(i,j,k) - TWO*alpn1(i,j,k)*kappa1*(Gamz(i,j,k)-Gamza(i,j,k))
        end do
      end do
    end do
!$OMP END PARALLEL DO

! numerical dissipation part  
  if(eps>0)then 
! usual Kreiss-Oliger dissipation    
  call kodis(ex,X,Y,Z,chi,chi_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,trK,trK_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,gxx,gxx_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,gxy,gxy_rhs,AAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,gxz,gxz_rhs,ASA,Symmetry,eps)
  call kodis(ex,X,Y,Z,gyy,gyy_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,gyz,gyz_rhs,SAA,Symmetry,eps)
  call kodis(ex,X,Y,Z,gzz,gzz_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Axx,Axx_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Axy,Axy_rhs,AAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Axz,Axz_rhs,ASA,Symmetry,eps)
  call kodis(ex,X,Y,Z,Ayy,Ayy_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Ayz,Ayz_rhs,SAA,Symmetry,eps)
  call kodis(ex,X,Y,Z,Azz,Azz_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Gamx,Gamx_rhs,ASS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Gamy,Gamy_rhs,SAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Gamz,Gamz_rhs,SSA,Symmetry,eps)
  call kodis(ex,X,Y,Z,Lap,Lap_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,betax,betax_rhs,ASS,Symmetry,eps)
  call kodis(ex,X,Y,Z,betay,betay_rhs,SAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,betaz,betaz_rhs,SSA,Symmetry,eps)
  call kodis(ex,X,Y,Z,dtSfx,dtSfx_rhs,ASS,Symmetry,eps)
  call kodis(ex,X,Y,Z,dtSfy,dtSfy_rhs,SAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,dtSfz,dtSfz_rhs,SSA,Symmetry,eps)
  call kodis(ex,X,Y,Z,TZ,TZ_rhs,SSS,Symmetry,eps)

  endif

#if (ABV == 0)  
  call ricci_gamma(ex, X, Y, Z,                                      &
               chi,                                                  &
               dxx    ,   gxy    ,   gxz    ,   dyy    ,   gyz    ,   dzz,&
               Gamx   ,  Gamy    ,  Gamz    , &
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,&
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,&
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,&
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,&
               Symmetry)
#endif

  call constraint_bssn(ex, X, Y, Z,&
               chi,trK, &
               dxx,gxy,gxz,dyy,gyz,dzz, &
               Axx,Axy,Axz,Ayy,Ayz,Azz, &
               Gamx,Gamy,Gamz,&
               Lap,betax,betay,betaz,rho,Sx,Sy,Sz,&
               Gamxxx, Gamxxy, Gamxxz,Gamxyy, Gamxyz, Gamxzz, &
               Gamyxx, Gamyxy, Gamyxz,Gamyyy, Gamyyz, Gamyzz, &
               Gamzxx, Gamzxy, Gamzxz,Gamzyy, Gamzyz, Gamzzz, &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz, &
               Hcon,Mxcon,Mycon,Mzcon,Gmxcon,Gmycon,Gmzcon, &
               Symmetry)

  gont = 0

  return

  end function compute_rhs_z4c_opt

  end function compute_rhs_z4c
#endif


!! using David Z4c-rhs code
#if 0
  function compute_rhs_z4c(ex, T,X, Y, Z,                      &
               chi    ,   trK    ,                                             &
               dxx    ,   gxy    ,   gxz    ,   dyy    ,   gyz    ,   dzz,     &
               Axx    ,   Axy    ,   Axz    ,   Ayy    ,   Ayz    ,   Azz,     &
               Gamx   ,  Gamy    ,  Gamz    ,                                  &
               Lap    ,  betax   ,  betay   ,  betaz   ,                       &
               dtSfx  ,  dtSfy   ,  dtSfz   ,                                  &
               TZ     ,                                                        &
               chi_rhs,   trK_rhs,                                             &
               gxx_rhs,   gxy_rhs,   gxz_rhs,   gyy_rhs,   gyz_rhs,   gzz_rhs, &
               Axx_rhs,   Axy_rhs,   Axz_rhs,   Ayy_rhs,   Ayz_rhs,   Azz_rhs, &
               Gamx_rhs,  Gamy_rhs,  Gamz_rhs,                                 &
               Lap_rhs,  betax_rhs,  betay_rhs,  betaz_rhs,                    &
               dtSfx_rhs,  dtSfy_rhs,  dtSfz_rhs,                              &
               TZ_rhs   ,                                                      &
               rho,Sx,Sy,Sz,Sxx,Sxy,Sxz,Syy,Syz,Szz,                           &
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,                      &
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,                      &
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,                      &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,                                        &
               Hcon,Mxcon,Mycon,Mzcon,Gmxcon,Gmycon,Gmzcon,                    &
! co is not used here, we always compute constraint               
               Symmetry,Lev,eps,co)  result(gont)
  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3), Symmetry,Lev,co
  real*8, intent(in ):: T
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: chi,dxx,dyy,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: trK
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: gxy,gxz,gyz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamx,Gamy,Gamz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Lap, betax, betay, betaz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: dtSfx,  dtSfy,  dtSfz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: TZ
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: chi_rhs,trK_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gxx_rhs,gxy_rhs,gxz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gyy_rhs,gyz_rhs,gzz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Axx_rhs,Axy_rhs,Axz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Ayy_rhs,Ayz_rhs,Azz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamx_rhs,Gamy_rhs,Gamz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Lap_rhs, betax_rhs, betay_rhs, betaz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: dtSfx_rhs,dtSfy_rhs,dtSfz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: TZ_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: rho,Sx,Sy,Sz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Sxx,Sxy,Sxz,Syy,Syz,Szz
! when out, physical second kind of connection  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Gamzyy, Gamzyz, Gamzzz
! when out, physical Ricci tensor  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Rxx,Rxy,Rxz,Ryy,Ryz,Rzz
! when out, constraint violation  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Hcon,Mxcon,Mycon,Mzcon,Gmxcon,Gmycon,Gmzcon
  real*8,intent(in) :: eps
!  gont = 0: success; gont = 1: something wrong
  integer::gont

!~~~~~~> Other variables:

  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gyy,gzz,alpn1,chin1
  real*8, dimension(ex(1),ex(2),ex(3)) :: chix,chiy,chiz,chixx,chixy,chixz,chiyy,chiyz,chizz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxx,gxyx,gxzx,gyyx,gyzx,gzzx
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxy,gxyy,gxzy,gyyy,gyzy,gzzy
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxz,gxyz,gxzz,gyyz,gyzz,gzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Lapx,Lapy,Lapz,Lapxx,Lapxy,Lapxz,Lapyy,Lapyz,Lapzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: betaxx,betaxy,betaxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: betayx,betayy,betayz
  real*8, dimension(ex(1),ex(2),ex(3)) :: betazx,betazy,betazz
  real*8, dimension(ex(1),ex(2),ex(3)) :: dBxx,dBxy,dBxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: dByx,dByy,dByz
  real*8, dimension(ex(1),ex(2),ex(3)) :: dBzx,dBzy,dBzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: sfxxx,sfxxy,sfxxz,sfxyy,sfxyz,sfxzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: sfyxx,sfyxy,sfyxz,sfyyy,sfyyz,sfyzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: sfzxx,sfzxy,sfzxz,sfzyy,sfzyz,sfzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamxx,Gamxy,Gamxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamyx,Gamyy,Gamyz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamzx,Gamzy,Gamzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Kx,Ky,Kz,TZx,TZy,TZz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxxx,gxxxy,gxxxz,gxxyy,gxxyz,gxxzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxyxx,gxyxy,gxyxz,gxyyy,gxyyz,gxyzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxzxx,gxzxy,gxzxz,gxzyy,gxzyz,gxzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyyxx,gyyxy,gyyxz,gyyyy,gyyyz,gyyzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyzxx,gyzxy,gyzxz,gyzyy,gyzyz,gyzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gzzxx,gzzxy,gzzxz,gzzyy,gzzyz,gzzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Axxx,Axyx,Axzx,Ayyx,Ayzx,Azzx
  real*8, dimension(ex(1),ex(2),ex(3)) :: Axxy,Axyy,Axzy,Ayyy,Ayzy,Azzy
  real*8, dimension(ex(1),ex(2),ex(3)) :: Axxz,Axyz,Axzz,Ayyz,Ayzz,Azzz

  real*8,dimension(3) ::SSS,AAS,ASA,SAA,ASS,SAS,SSA
  real*8            :: dX, dY, dZ, PI
  real*8, parameter :: ZEO=0.d0,ONE = 1.D0, TWO = 2.D0, FOUR = 4.D0,F16=1.6d1
  real*8, parameter :: EIGHT = 8.D0, HALF = 0.5D0, THR = 3.d0,F8=8.d0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  real*8, parameter :: F1o3 = 1.D0/3.D0, F2o3 = 2.D0/3.D0,F3o2=1.5d0, F1o6 = 1.D0/6.D0
  integer :: i,j,k

! constraint damping terms stuffs PRD 81, 084003 (2010)
  real*8 :: kappa1,kappa2,kappa3,FF,eta

  real*8,parameter :: chiDivfloor=1.d-5

  call get_Z4cparameters(kappa1,kappa2,kappa3,FF,eta)

!!! sanity check
  dX = sum(chi)+sum(trK)+sum(dxx)+sum(gxy)+sum(gxz)+sum(dyy)+sum(gyz)+sum(dzz) &
      +sum(Axx)+sum(Axy)+sum(Axz)+sum(Ayy)+sum(Ayz)+sum(Azz)                   &
      +sum(Gamx)+sum(Gamy)+sum(Gamz)                                           &
      +sum(Lap)+sum(betax)+sum(betay)+sum(betaz)+sum(dtSfx)+sum(dtSfy)+sum(dtSfz) &
      +sum(TZ)
  if(dX.ne.dX) then
     if(sum(chi).ne.sum(chi))write(*,*)"Z4c_rhs.f90: find NaN in chi"
     if(sum(trK).ne.sum(trK))write(*,*)"Z4c_rhs.f90: find NaN in trk"
     if(sum(dxx).ne.sum(dxx))write(*,*)"Z4c_rhs.f90: find NaN in gxx"
     if(sum(gxy).ne.sum(gxy))write(*,*)"Z4c_rhs.f90: find NaN in gxy"
     if(sum(gxz).ne.sum(gxz))write(*,*)"Z4c_rhs.f90: find NaN in gxz"
     if(sum(dyy).ne.sum(dyy))write(*,*)"Z4c_rhs.f90: find NaN in gyy"
     if(sum(gyz).ne.sum(gyz))write(*,*)"Z4c_rhs.f90: find NaN in gyz"
     if(sum(dzz).ne.sum(dzz))write(*,*)"Z4c_rhs.f90: find NaN in gzz"
     if(sum(Axx).ne.sum(Axx))write(*,*)"Z4c_rhs.f90: find NaN in Axx"
     if(sum(Axy).ne.sum(Axy))write(*,*)"Z4c_rhs.f90: find NaN in Axy"
     if(sum(Axz).ne.sum(Axz))write(*,*)"Z4c_rhs.f90: find NaN in Axz"
     if(sum(Ayy).ne.sum(Ayy))write(*,*)"Z4c_rhs.f90: find NaN in Ayy"
     if(sum(Ayz).ne.sum(Ayz))write(*,*)"Z4c_rhs.f90: find NaN in Ayz"
     if(sum(Azz).ne.sum(Azz))write(*,*)"Z4c_rhs.f90: find NaN in Azz"
     if(sum(Gamx).ne.sum(Gamx))write(*,*)"Z4c_rhs.f90: find NaN in Gamx"
     if(sum(Gamy).ne.sum(Gamy))write(*,*)"Z4c_rhs.f90: find NaN in Gamy"
     if(sum(Gamz).ne.sum(Gamz))write(*,*)"Z4c_rhs.f90: find NaN in Gamz"
     if(sum(Lap).ne.sum(Lap))write(*,*)"Z4c_rhs.f90: find NaN in Lap"
     if(sum(betax).ne.sum(betax))write(*,*)"Z4c_rhs.f90: find NaN in betax"
     if(sum(betay).ne.sum(betay))write(*,*)"Z4c_rhs.f90: find NaN in betay"
     if(sum(betaz).ne.sum(betaz))write(*,*)"Z4c_rhs.f90: find NaN in betaz"
     if(sum(dtSfx).ne.sum(dtSfx))write(*,*)"Z4c_rhs.f90: find NaN in dtSfx"
     if(sum(dtSfy).ne.sum(dtSfy))write(*,*)"Z4c_rhs.f90: find NaN in dtSfy"
     if(sum(dtSfz).ne.sum(dtSfz))write(*,*)"Z4c_rhs.f90: find NaN in dtSfz"
     if(sum(TZ).ne.sum(Tz))write(*,*)"Z4c_rhs.f90: find NaN in TZ"
     gont = 1
     return
  endif

  PI = dacos(-ONE)

  dX = X(2) - X(1)
  dY = Y(2) - Y(1)
  dZ = Z(2) - Z(1)

  alpn1 = Lap + ONE
  chin1 = chi + ONE
  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE

  call fderivs(ex,betax,betaxx,betaxy,betaxz,X,Y,Z,ANTI, SYM, SYM,Symmetry,Lev)
  call fderivs(ex,betay,betayx,betayy,betayz,X,Y,Z, SYM,ANTI, SYM,Symmetry,Lev)
  call fderivs(ex,betaz,betazx,betazy,betazz,X,Y,Z, SYM, SYM,ANTI,Symmetry,Lev)
  call fderivs(ex,dtSfx,dBxx,dBxy,dBxz,X,Y,Z,ANTI, SYM, SYM,Symmetry,Lev)
  call fderivs(ex,dtSfy,dByx,dByy,dByz,X,Y,Z, SYM,ANTI, SYM,Symmetry,Lev)
  call fderivs(ex,dtSfz,dBzx,dBzy,dBzz,X,Y,Z, SYM, SYM,ANTI,Symmetry,Lev)
  call fderivs(ex,chi,chix,chiy,chiz,X,Y,Z, SYM, SYM,SYM,Symmetry,Lev)
  call fderivs(ex,dxx,gxxx,gxxy,gxxz,X,Y,Z, SYM, SYM,SYM,Symmetry,Lev)
  call fderivs(ex,gxy,gxyx,gxyy,gxyz,X,Y,Z,ANTI,ANTI,SYM,Symmetry,Lev)
  call fderivs(ex,gxz,gxzx,gxzy,gxzz,X,Y,Z,ANTI,SYM,ANTI,Symmetry,Lev)
  call fderivs(ex,dyy,gyyx,gyyy,gyyz,X,Y,Z, SYM, SYM,SYM,Symmetry,Lev)
  call fderivs(ex,gyz,gyzx,gyzy,gyzz,X,Y,Z,SYM,ANTI,ANTI,Symmetry,Lev)
  call fderivs(ex,dzz,gzzx,gzzy,gzzz,X,Y,Z, SYM, SYM,SYM,Symmetry,Lev)

  call fdderivs(ex,dxx,gxxxx,gxxxy,gxxxz,gxxyy,gxxyz,gxxzz,X,Y,Z, SYM, SYM,SYM ,Symmetry,Lev)
  call fdderivs(ex,dyy,gyyxx,gyyxy,gyyxz,gyyyy,gyyyz,gyyzz,X,Y,Z, SYM, SYM,SYM ,Symmetry,Lev)
  call fdderivs(ex,dzz,gzzxx,gzzxy,gzzxz,gzzyy,gzzyz,gzzzz,X,Y,Z, SYM, SYM,SYM ,Symmetry,Lev)
  call fdderivs(ex,gxy,gxyxx,gxyxy,gxyxz,gxyyy,gxyyz,gxyzz,X,Y,Z,ANTI,ANTI,SYM ,Symmetry,Lev)
  call fdderivs(ex,gxz,gxzxx,gxzxy,gxzxz,gxzyy,gxzyz,gxzzz,X,Y,Z,ANTI,SYM ,ANTI,Symmetry,Lev)
  call fdderivs(ex,gyz,gyzxx,gyzxy,gyzxz,gyzyy,gyzyz,gyzzz,X,Y,Z,SYM ,ANTI,ANTI,Symmetry,Lev)

  call fderivs(ex,Gamx,Gamxx,Gamxy,Gamxz,X,Y,Z,ANTI,SYM ,SYM,Symmetry,Lev)
  call fderivs(ex,Gamy,Gamyx,Gamyy,Gamyz,X,Y,Z,SYM ,ANTI,SYM,Symmetry,Lev)
  call fderivs(ex,Gamz,Gamzx,Gamzy,Gamzz,X,Y,Z,SYM ,SYM,ANTI,Symmetry,Lev)

  call fderivs(ex,Lap,Lapx,Lapy,Lapz,X,Y,Z, SYM, SYM,SYM,Symmetry,Lev)
  call fderivs(ex,trK,Kx,Ky,Kz,X,Y,Z, SYM, SYM,SYM,Symmetry,Lev)

  call fderivs(ex,TZ,TZx,TZy,TZz,X,Y,Z, SYM, SYM,SYM,Symmetry,Lev)
         
  call fdderivs(ex,betax,sfxxx,sfxxy,sfxxz,sfxyy,sfxyz,sfxzz,X,Y,Z,ANTI, SYM, SYM,Symmetry,Lev)
  call fdderivs(ex,betay,sfyxx,sfyxy,sfyxz,sfyyy,sfyyz,sfyzz,X,Y,Z, SYM,ANTI, SYM,Symmetry,Lev)
  call fdderivs(ex,betaz,sfzxx,sfzxy,sfzxz,sfzyy,sfzyz,sfzzz,X,Y,Z, SYM, SYM,ANTI,Symmetry,Lev)

  call fdderivs(ex,chi,chixx,chixy,chixz,chiyy,chiyz,chizz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)

  call fdderivs(ex,Lap,Lapxx,Lapxy,Lapxz,Lapyy,Lapyz,Lapzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)

  call fderivs(ex,Axx,Axxx,Axxy,Axxz,X,Y,Z, SYM, SYM,SYM,Symmetry,Lev)
  call fderivs(ex,Axy,Axyx,Axyy,Axyz,X,Y,Z,ANTI,ANTI,SYM,Symmetry,Lev)
  call fderivs(ex,Axz,Axzx,Axzy,Axzz,X,Y,Z,ANTI,SYM,ANTI,Symmetry,Lev)
  call fderivs(ex,Ayy,Ayyx,Ayyy,Ayyz,X,Y,Z, SYM, SYM,SYM,Symmetry,Lev)
  call fderivs(ex,Ayz,Ayzx,Ayzy,Ayzz,X,Y,Z,SYM,ANTI,ANTI,Symmetry,Lev)
  call fderivs(ex,Azz,Azzx,Azzy,Azzz,X,Y,Z, SYM, SYM,SYM,Symmetry,Lev)

  do i=1,ex(1)
  do j=1,ex(2)
  do k=1,ex(3)
     call z4c_rhs_point(Axx(i,j,k),Axy(i,j,k),Axz(i,j,k),Ayy(i,j,k),Ayz(i,j,k),Azz(i,j,k), &
                        alpn1(i,j,k),dtSfx(i,j,k),dtSfy(i,j,k),dtSfz(i,j,k), &
                        betax(i,j,k),betay(i,j,k),betaz(i,j,k), &
                        chin1(i,j,k),chiDivfloor, &
                        Lapx(i,j,k), &
                        Axxx(i,j,k),Axyx(i,j,k),Axzx(i,j,k),Ayyx(i,j,k),Ayzx(i,j,k),Azzx(i,j,k), &
                        Lapy(i,j,k), &
                        Axxy(i,j,k),Axyy(i,j,k),Axzy(i,j,k),Ayyy(i,j,k),Ayzy(i,j,k),Azzy(i,j,k), &
                        Lapz(i,j,k), &
                        Axxz(i,j,k),Axyz(i,j,k),Axzz(i,j,k),Ayyz(i,j,k),Ayzz(i,j,k),Azzz(i,j,k), &
                        betaxx(i,j,k),dBxx(i,j,k),betayx(i,j,k),dByx(i,j,k),betazx(i,j,k),dBzx(i,j,k), &
                        betaxy(i,j,k),dBxy(i,j,k),betayy(i,j,k),dByy(i,j,k),betazy(i,j,k),dBzy(i,j,k), &
                        betaxz(i,j,k),dBxz(i,j,k),betayz(i,j,k),dByz(i,j,k),betazz(i,j,k),dBzz(i,j,k), &
                        chix(i,j,k),chiy(i,j,k),chiz(i,j,k), &
                        Lapxx(i,j,k),Lapxy(i,j,k),Lapxz(i,j,k),Lapyy(i,j,k),Lapyz(i,j,k),Lapzz(i,j,k), &
                        sfxxx(i,j,k),sfyxx(i,j,k),sfzxx(i,j,k), &
                        sfxxy(i,j,k),sfyxy(i,j,k),sfzxy(i,j,k), &
                        sfxxz(i,j,k),sfyxz(i,j,k),sfzxz(i,j,k), &
                        sfxyy(i,j,k),sfyyy(i,j,k),sfzyy(i,j,k), &
                        sfxyz(i,j,k),sfyyz(i,j,k),sfzyz(i,j,k), &
                        sfxzz(i,j,k),sfyzz(i,j,k),sfzzz(i,j,k), &
                        chixx(i,j,k),chixy(i,j,k),chixz(i,j,k),chiyy(i,j,k),chiyz(i,j,k),chizz(i,j,k), &
                        gxxxx(i,j,k),gxyxx(i,j,k),gxzxx(i,j,k),gyyxx(i,j,k),gyzxx(i,j,k),gzzxx(i,j,k), &
                        gxxxy(i,j,k),gxyxy(i,j,k),gxzxy(i,j,k),gyyxy(i,j,k),gyzxy(i,j,k),gzzxy(i,j,k), &
                        gxxxz(i,j,k),gxyxz(i,j,k),gxzxz(i,j,k),gyyxz(i,j,k),gyzxz(i,j,k),gzzxz(i,j,k), &
                        gxxyy(i,j,k),gxyyy(i,j,k),gxzyy(i,j,k),gyyyy(i,j,k),gyzyy(i,j,k),gzzyy(i,j,k), &
                        gxxyz(i,j,k),gxyyz(i,j,k),gxzyz(i,j,k),gyyyz(i,j,k),gyzyz(i,j,k),gzzyz(i,j,k), &
                        gxxzz(i,j,k),gxyzz(i,j,k),gxzzz(i,j,k),gyyzz(i,j,k),gyzzz(i,j,k),gzzzz(i,j,k), &
                        Gamxx(i,j,k),gxxx(i,j,k),gxyx(i,j,k),gxzx(i,j,k), &
                        Gamyx(i,j,k),gyyx(i,j,k),gyzx(i,j,k), &
                        Gamzx(i,j,k),gzzx(i,j,k), &
                        Gamxy(i,j,k),gxxy(i,j,k),gxyy(i,j,k),gxzy(i,j,k), &
                        Gamyy(i,j,k),gyyy(i,j,k),gyzy(i,j,k), &
                        Gamzy(i,j,k),gzzy(i,j,k), &
                        Gamxz(i,j,k),gxxz(i,j,k),gxyz(i,j,k),gxzz(i,j,k), &
                        Gamyz(i,j,k),gyyz(i,j,k),gyzz(i,j,k), &
                        Gamzz(i,j,k),gzzz(i,j,k), &
                        Kx(i,j,k),Ky(i,j,k),Kz(i,j,k), &
                        TZx(i,j,k),TZy(i,j,k),TZz(i,j,k), &
                        Gamx(i,j,k),gxx(i,j,k),gxy(i,j,k),gxz(i,j,k), &
                        Gamy(i,j,k),gyy(i,j,k),gyz(i,j,k), &
                        Gamz(i,j,k),gzz(i,j,k), &
                        kappa1,kappa2, &
                        trK(i,j,k), &
                        Axx_rhs(i,j,k),Axy_rhs(i,j,k),Axz_rhs(i,j,k),Ayy_rhs(i,j,k),Ayz_rhs(i,j,k),Azz_rhs(i,j,k), &
                        chi_rhs(i,j,k), &
                        Gamx_rhs(i,j,k),gxx_rhs(i,j,k),gxy_rhs(i,j,k),gxz_rhs(i,j,k), &
                        Gamy_rhs(i,j,k),gyy_rhs(i,j,k),gyz_rhs(i,j,k), &
                        Gamz_rhs(i,j,k),gzz_rhs(i,j,k),trK_rhs(i,j,k),TZ_rhs(i,j,k),TZ(i,j,k))
  enddo
  enddo
  enddo

!!!!!gauge variable part
  Lap_rhs = -TWO*alpn1*trK
#if (GAUGE == 0)
  betax_rhs = FF*dtSfx
  betay_rhs = FF*dtSfy
  betaz_rhs = FF*dtSfz

  dtSfx_rhs = Gamx_rhs - eta*dtSfx
  dtSfy_rhs = Gamy_rhs - eta*dtSfy
  dtSfz_rhs = Gamz_rhs - eta*dtSfz
#elif (GAUGE == 1)
  betax_rhs = Gamx - eta*betax
  betay_rhs = Gamy - eta*betay
  betaz_rhs = Gamz - eta*betaz

  dtSfx_rhs = ZEO
  dtSfy_rhs = ZEO
  dtSfz_rhs = ZEO
#endif  

  SSS(1)=SYM
  SSS(2)=SYM
  SSS(3)=SYM

  AAS(1)=ANTI
  AAS(2)=ANTI
  AAS(3)=SYM

  ASA(1)=ANTI
  ASA(2)=SYM
  ASA(3)=ANTI

  SAA(1)=SYM
  SAA(2)=ANTI
  SAA(3)=ANTI

  ASS(1)=ANTI
  ASS(2)=SYM
  ASS(3)=SYM

  SAS(1)=SYM
  SAS(2)=ANTI
  SAS(3)=SYM

  SSA(1)=SYM
  SSA(2)=SYM
  SSA(3)=ANTI
!!!!!!!!!advection term part
  call lopsided(ex,X,Y,Z,gxx,gxx_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,gxy,gxy_rhs,betax,betay,betaz,Symmetry,AAS)
  call lopsided(ex,X,Y,Z,gxz,gxz_rhs,betax,betay,betaz,Symmetry,ASA)
  call lopsided(ex,X,Y,Z,gyy,gyy_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,gyz,gyz_rhs,betax,betay,betaz,Symmetry,SAA)
  call lopsided(ex,X,Y,Z,gzz,gzz_rhs,betax,betay,betaz,Symmetry,SSS)

  call lopsided(ex,X,Y,Z,Axx,Axx_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,Axy,Axy_rhs,betax,betay,betaz,Symmetry,AAS)
  call lopsided(ex,X,Y,Z,Axz,Axz_rhs,betax,betay,betaz,Symmetry,ASA)
  call lopsided(ex,X,Y,Z,Ayy,Ayy_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,Ayz,Ayz_rhs,betax,betay,betaz,Symmetry,SAA)
  call lopsided(ex,X,Y,Z,Azz,Azz_rhs,betax,betay,betaz,Symmetry,SSS)

  call lopsided(ex,X,Y,Z,chi,chi_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,trK,trK_rhs,betax,betay,betaz,Symmetry,SSS)

  call lopsided(ex,X,Y,Z,Gamx,Gamx_rhs,betax,betay,betaz,Symmetry,ASS)
  call lopsided(ex,X,Y,Z,Gamy,Gamy_rhs,betax,betay,betaz,Symmetry,SAS)
  call lopsided(ex,X,Y,Z,Gamz,Gamz_rhs,betax,betay,betaz,Symmetry,SSA)

  call lopsided(ex,X,Y,Z,Lap,Lap_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,betax,betax_rhs,betax,betay,betaz,Symmetry,ASS)
  call lopsided(ex,X,Y,Z,betay,betay_rhs,betax,betay,betaz,Symmetry,SAS)
  call lopsided(ex,X,Y,Z,betaz,betaz_rhs,betax,betay,betaz,Symmetry,SSA)

#if (GAUGE == 0)
  call lopsided(ex,X,Y,Z,dtSfx,dtSfx_rhs,betax,betay,betaz,Symmetry,ASS)
  call lopsided(ex,X,Y,Z,dtSfy,dtSfy_rhs,betax,betay,betaz,Symmetry,SAS)
  call lopsided(ex,X,Y,Z,dtSfz,dtSfz_rhs,betax,betay,betaz,Symmetry,SSA)
#endif

  call lopsided(ex,X,Y,Z,TZ,TZ_rhs,betax,betay,betaz,Symmetry,SSS)
! numerical dissipation part  
  if(eps>0)then 
! usual Kreiss-Oliger dissipation   

  call kodis(ex,X,Y,Z,chi,chi_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,trK,trK_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,dxx,gxx_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,gxy,gxy_rhs,AAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,gxz,gxz_rhs,ASA,Symmetry,eps)
  call kodis(ex,X,Y,Z,dyy,gyy_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,gyz,gyz_rhs,SAA,Symmetry,eps)
  call kodis(ex,X,Y,Z,dzz,gzz_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Axx,Axx_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Axy,Axy_rhs,AAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Axz,Axz_rhs,ASA,Symmetry,eps)
  call kodis(ex,X,Y,Z,Ayy,Ayy_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Ayz,Ayz_rhs,SAA,Symmetry,eps)
  call kodis(ex,X,Y,Z,Azz,Azz_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Gamx,Gamx_rhs,ASS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Gamy,Gamy_rhs,SAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Gamz,Gamz_rhs,SSA,Symmetry,eps)
  call kodis(ex,X,Y,Z,Lap,Lap_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,betax,betax_rhs,ASS,Symmetry,eps)
  call kodis(ex,X,Y,Z,betay,betay_rhs,SAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,betaz,betaz_rhs,SSA,Symmetry,eps)
#if (GAUGE == 0)
  call kodis(ex,X,Y,Z,dtSfx,dtSfx_rhs,ASS,Symmetry,eps)
  call kodis(ex,X,Y,Z,dtSfy,dtSfy_rhs,SAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,dtSfz,dtSfz_rhs,SSA,Symmetry,eps)
#endif
  call kodis(ex,X,Y,Z,TZ,TZ_rhs,SSS,Symmetry,eps)

  endif

#if (ABV == 0)  
  call ricci_gamma(ex, X, Y, Z,                                      &
               chi,                                                  &
               dxx    ,   gxy    ,   gxz    ,   dyy    ,   gyz    ,   dzz,&
               Gamx   ,  Gamy    ,  Gamz    , &
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,&
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,&
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,&
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,&
               Symmetry)
#endif

  call constraint_bssn(ex, X, Y, Z,&
               chi,trK, &
               dxx,gxy,gxz,dyy,gyz,dzz, &
               Axx,Axy,Axz,Ayy,Ayz,Azz, &
               Gamx,Gamy,Gamz,&
               Lap,betax,betay,betaz,rho,Sx,Sy,Sz,&
               Gamxxx, Gamxxy, Gamxxz,Gamxyy, Gamxyz, Gamxzz, &
               Gamyxx, Gamyxy, Gamyxz,Gamyyy, Gamyyz, Gamyzz, &
               Gamzxx, Gamzxy, Gamzxz,Gamzyy, Gamzyz, Gamzzz, &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz, &
               Hcon,Mxcon,Mycon,Mzcon,Gmxcon,Gmycon,Gmzcon, &
               Symmetry)

  gont = 0

  return

  end function compute_rhs_Z4c
#endif
