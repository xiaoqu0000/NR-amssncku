

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

  trKd = trK+TWO*TZ
!this beta^i_,j will be kept till the end of this routine
  call fderivs(ex,betax,betaxx,betaxy,betaxz,X,Y,Z,ANTI, SYM, SYM,Symmetry,Lev)
  call fderivs(ex,betay,betayx,betayy,betayz,X,Y,Z, SYM,ANTI, SYM,Symmetry,Lev)
  call fderivs(ex,betaz,betazx,betazy,betazz,X,Y,Z, SYM, SYM,ANTI,Symmetry,Lev)

  div_beta = betaxx + betayy + betazz
 
  call fderivs(ex,chi,chix,chiy,chiz,X,Y,Z,SYM,SYM,SYM,symmetry,Lev)

  chi_rhs = F2o3 *chin1*( alpn1 * trKd - div_beta ) !rhs for chi

  call fderivs(ex,dxx,gxxx,gxxy,gxxz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)
  call fderivs(ex,gxy,gxyx,gxyy,gxyz,X,Y,Z,ANTI,ANTI,SYM ,Symmetry,Lev)
  call fderivs(ex,gxz,gxzx,gxzy,gxzz,X,Y,Z,ANTI,SYM ,ANTI,Symmetry,Lev)
  call fderivs(ex,dyy,gyyx,gyyy,gyyz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)
  call fderivs(ex,gyz,gyzx,gyzy,gyzz,X,Y,Z,SYM ,ANTI,ANTI,Symmetry,Lev)
  call fderivs(ex,dzz,gzzx,gzzy,gzzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)

  gxx_rhs = - TWO * alpn1 * Axx    -  F2o3 * gxx * div_beta          + &
              TWO *(  gxx * betaxx +   gxy * betayx +   gxz * betazx)

  gyy_rhs = - TWO * alpn1 * Ayy    -  F2o3 * gyy * div_beta          + &
              TWO *(  gxy * betaxy +   gyy * betayy +   gyz * betazy)

  gzz_rhs = - TWO * alpn1 * Azz    -  F2o3 * gzz * div_beta          + &
              TWO *(  gxz * betaxz +   gyz * betayz +   gzz * betazz)

  gxy_rhs = - TWO * alpn1 * Axy    +  F1o3 * gxy    * div_beta       + &
                      gxx * betaxy                  +   gxz * betazy + &
                                       gyy * betayx +   gyz * betazx   &
                                                    -   gxy * betazz

  gyz_rhs = - TWO * alpn1 * Ayz    +  F1o3 * gyz    * div_beta       + &
                      gxy * betaxz +   gyy * betayz                  + &
                      gxz * betaxy                  +   gzz * betazy   &
                                                    -   gyz * betaxx
 
  gxz_rhs = - TWO * alpn1 * Axz    +  F1o3 * gxz    * div_beta       + &
                      gxx * betaxz +   gxy * betayz                  + &
                                       gyz * betayx +   gzz * betazx   &
                                                    -   gxz * betayy     !rhs for gij

! invert tilted metric
  gupzz =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz - &
           gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz
  gupxx =   ( gyy * gzz - gyz * gyz ) / gupzz
  gupxy = - ( gxy * gzz - gyz * gxz ) / gupzz
  gupxz =   ( gxy * gyz - gyy * gxz ) / gupzz
  gupyy =   ( gxx * gzz - gxz * gxz ) / gupzz
  gupyz = - ( gxx * gyz - gxy * gxz ) / gupzz
  gupzz =   ( gxx * gyy - gxy * gxy ) / gupzz
! gij_,kl will be stored till end of this routine
  call fdderivs(ex,dxx,gxxxx,gxxxy,gxxxz,gxxyy,gxxyz,gxxzz,X,Y,Z,SYM ,SYM ,SYM ,symmetry,Lev)
  call fdderivs(ex,dyy,gyyxx,gyyxy,gyyxz,gyyyy,gyyyz,gyyzz,X,Y,Z,SYM ,SYM ,SYM ,symmetry,Lev)
  call fdderivs(ex,dzz,gzzxx,gzzxy,gzzxz,gzzyy,gzzyz,gzzzz,X,Y,Z,SYM ,SYM ,SYM ,symmetry,Lev)
  call fdderivs(ex,gxy,gxyxx,gxyxy,gxyxz,gxyyy,gxyyz,gxyzz,X,Y,Z,ANTI,ANTI,SYM ,symmetry,Lev)
  call fdderivs(ex,gxz,gxzxx,gxzxy,gxzxz,gxzyy,gxzyz,gxzzz,X,Y,Z,ANTI,SYM ,ANTI,symmetry,Lev)
  call fdderivs(ex,gyz,gyzxx,gyzxy,gyzxz,gyzyy,gyzyz,gyzzz,X,Y,Z,SYM ,ANTI,ANTI,symmetry,Lev)
! second kind of connection
  Gamxxx =HALF*( gupxx*gxxx + gupxy*(TWO*gxyx - gxxy ) + gupxz*(TWO*gxzx - gxxz ))
  Gamyxx =HALF*( gupxy*gxxx + gupyy*(TWO*gxyx - gxxy ) + gupyz*(TWO*gxzx - gxxz ))
  Gamzxx =HALF*( gupxz*gxxx + gupyz*(TWO*gxyx - gxxy ) + gupzz*(TWO*gxzx - gxxz ))
 
  Gamxyy =HALF*( gupxx*(TWO*gxyy - gyyx ) + gupxy*gyyy + gupxz*(TWO*gyzy - gyyz ))
  Gamyyy =HALF*( gupxy*(TWO*gxyy - gyyx ) + gupyy*gyyy + gupyz*(TWO*gyzy - gyyz ))
  Gamzyy =HALF*( gupxz*(TWO*gxyy - gyyx ) + gupyz*gyyy + gupzz*(TWO*gyzy - gyyz ))

  Gamxzz =HALF*( gupxx*(TWO*gxzz - gzzx ) + gupxy*(TWO*gyzz - gzzy ) + gupxz*gzzz)
  Gamyzz =HALF*( gupxy*(TWO*gxzz - gzzx ) + gupyy*(TWO*gyzz - gzzy ) + gupyz*gzzz)
  Gamzzz =HALF*( gupxz*(TWO*gxzz - gzzx ) + gupyz*(TWO*gyzz - gzzy ) + gupzz*gzzz)

  Gamxxy =HALF*( gupxx*gxxy + gupxy*gyyx + gupxz*( gxzy + gyzx - gxyz ) )
  Gamyxy =HALF*( gupxy*gxxy + gupyy*gyyx + gupyz*( gxzy + gyzx - gxyz ) )
  Gamzxy =HALF*( gupxz*gxxy + gupyz*gyyx + gupzz*( gxzy + gyzx - gxyz ) )

  Gamxxz =HALF*( gupxx*gxxz + gupxy*( gxyz + gyzx - gxzy ) + gupxz*gzzx )
  Gamyxz =HALF*( gupxy*gxxz + gupyy*( gxyz + gyzx - gxzy ) + gupyz*gzzx )
  Gamzxz =HALF*( gupxz*gxxz + gupyz*( gxyz + gyzx - gxzy ) + gupzz*gzzx )

  Gamxyz =HALF*( gupxx*( gxyz + gxzy - gyzx ) + gupxy*gyyz + gupxz*gzzy )
  Gamyyz =HALF*( gupxy*( gxyz + gxzy - gyzx ) + gupyy*gyyz + gupyz*gzzy )
  Gamzyz =HALF*( gupxz*( gxyz + gxzy - gyzx ) + gupyz*gyyz + gupzz*gzzy )
! the so called Gamma_d
  Gamxa =       gupxx * Gamxxx + gupyy * Gamxyy + gupzz * Gamxzz + &
          TWO*( gupxy * Gamxxy + gupxz * Gamxxz + gupyz * Gamxyz )
  Gamya =       gupxx * Gamyxx + gupyy * Gamyyy + gupzz * Gamyzz + &
          TWO*( gupxy * Gamyxy + gupxz * Gamyxz + gupyz * Gamyyz )
  Gamza =       gupxx * Gamzxx + gupyy * Gamzyy + gupzz * Gamzzz + &
          TWO*( gupxy * Gamzxy + gupxz * Gamzxz + gupyz * Gamzyz )

!!!!!!!!!!!!because gij_,k will be overwrite later, we calculate TWO*d_k Z^k here
! use Gamma^i as more as possible
  Gmxcon = Gamx - Gamxa
  Gmycon = Gamy - Gamya
  Gmzcon = Gamz - Gamza

!Maple generated code for g^ki*g^jm*g^ln*g_mn,k*g_ij,l
! Gami_,j are used as maple temp variables
      Gamyy = 3*gupxz**2*gupzz*gxzz**2+gupxx*gupxz**2*gxxz**2+2*gxyx*gupxy**3*gxyy+ &
           2*gxyx*gupxy**3*gyyx+gupxx**2*gupzz*gxzx**2+3*gupxx*gupxy**2*gxyx**2+ &
           6*gxyx*gupxy*gupxz*gupyy*gyzy+gupxx**2*gupyy*gxyx**2+ &
           2*gxyz*gupxy*gupyz**2*gyyz+2*gxxz*gupxx**2*gupyz*gxyx+ &
           gupxz**2*gupyy*gyzx**2+2*gxxy*gupxx*gupxy*gupxz*gxxz+ &
           2*gyzx*gupxy*gupxz*gupzz*gzzx+3*gupyy*gupyz**2*gyzy**2+ &
           2*gyyy*gupyz**3*gzzz+2*gxxz*gupxz**3*gxzz+ &
           4*gxzy*gupxx*gupxz*gupyy*gxyx+gupyy*gupyz**2*gyyz**2
      Gamxz = Gamyy+2*gxxz*gupxy**2*gupzz*gyzy+4*gxyz*gupxx*gupxy*gupxz*gxxx+ &
           6*gxzz*gupxy*gupyz*gupzz*gyzy+2*gxxy*gupxx*gupxz*gupyz*gxzz+ &
           3*gupxy**2*gupyy*gxyy**2+2*gxyz*gupxx*gupyy*gupzz*gyzx+ &
           4*gxyy*gupxx*gupyy*gupyz*gyzx+6*gxyy*gupxy*gupxz*gupyz*gxzz+ &
           4*gxzz*gupxx*gupyz*gupzz*gyzx+3*gupxx*gupxz**2*gxzx**2+ &
           4*gxyz*gupxx*gupxy*gupyz*gxyx+2*gxxz*gupxx*gupxz*gupyz*gxyz+ &
           2*gxxy*gupxy*gupxz*gupyz*gyzz+2*gxzx*gupxy*gupxz*gupyz*gyyz+ &
           gupyz**2*gupzz*gzzy**2+gupxz**2*gupzz*gzzx**2+ &
           gupyy*gupzz**2*gyzz**2+2*gyzy*gupyz**3*gzzy+gupxx*gupzz**2*gxzz**2
      Gamyy = Gamxz+gupxx*gupyz**2*gxzy**2+2*gxzx*gupxz**3*gzzx+ &
           3*gupyz**2*gupzz*gyzz**2+2*gyzy*gupyz**3*gyzz+gupyy**2*gupzz*gyzy**2+ &
           gupxy**2*gupzz*gyzx**2+2*gyyz*gupyz**3*gyzz+gupxy**2*gupyy*gyyx**2+ &
           gupxx*gupyz**2*gxyz**2+gupxx*gupyy**2*gxyy**2+ &
           gupxy**2*gupzz*gxzy**2+2*gxzx*gupxz**3*gxzz+ &
           2*gyyx*gupxy*gupxz*gupyy*gyzx+gupxx*gupxy**2*gxxy**2+ &
           2*gxxx*gupxz**3*gzzz+2*gxxx*gupxy**3*gyyy+gupxz**2*gupyy*gxyz**2+ &
           2*gxyy*gupxy**3*gxxy
      Gamxy = Gamyy+2*gxyy*gupxz*gupyy**2*gyzy+6*gxyy*gupxx*gupxy*gupyz*gxzx+ &
           4*gxyy*gupxy*gupxz*gupyy*gxyz+2*gyzx*gupxz*gupyy*gupzz*gzzy+ &
           2*gxzy*gupxy*gupxz*gupyy*gxyy+4*gxzy*gupxy*gupxz*gupzz*gxzz+ &
           2*gyyx*gupxz*gupyy*gupyz*gyzz+6*gxyx*gupxx*gupxz*gupyz*gxzz+ &
           2*gxyz*gupxy**2*gupzz*gxzy+2*gxyz*gupxy**2*gupyz*gxyy+ &
           2*gxyz*gupxy**2*gupxz*gxxy+2*gupxy*gupxz*gupyz*gxyz**2+ &
           4*gxyy*gupxz*gupyz**2*gzzz+2*gxyy*gupxy*gupyz**2*gzzy+ &
           4*gxyy*gupxy**2*gupyz*gxzy+2*gxyy*gupxy**2*gupxz*gxxz+ &
           4*gxyy*gupxx*gupxy**2*gxxx+2*gxyx*gupxy**2*gupxz*gxzy+ &
           2*gxyx*gupxy**2*gupyz*gyzy
      Gamyy = Gamxy+2*gxyx*gupxx*gupxy**2*gxxy+4*gyzz*gupyz*gupzz**2*gzzz+ &
           4*gxzy*gupxx*gupxz*gupyz*gxzx+2*gxzy*gupxx*gupyy*gupzz*gyzx+ &
           4*gxxx*gupxx*gupxy*gupxz*gyzx+2*gxyx*gupxx**2*gupyz*gxzx+ &
           2*gxyx*gupxy**2*gupxz*gxyz+2*gxzy*gupxz*gupyy*gupyz*gyyz+ &
           4*gxzy*gupxy*gupyy*gupyz*gyyy+2*gxzy*gupxx*gupyy*gupyz*gyyx+ &
           2*gyyx*gupxy*gupxz*gupyy*gxzy+2*gyyx*gupxy*gupyy*gupyz*gyyz+ &
           2*gyyx*gupxy*gupyy*gupyz*gyzy+4*gxzy*gupxz*gupyy*gupzz*gyzz+ &
           2*gyyx*gupxy*gupxz*gupyz*gxzz+2*gxyz*gupxx*gupyy*gupzz*gxzy+ &
           2*gxyy*gupxz*gupyy*gupyz*gzzy
      Gamxz = Gamyy+2*gxyy*gupxy*gupxz*gupyz*gzzx+2*gxyy*gupxy*gupxz*gupyy*gyzx+ &
           2*gxyy*gupxy*gupyy*gupyz*gyyz+2*gxyy*gupxx*gupyy*gupyz*gxzy+ &
           2*gxxy*gupxy**2*gupxz*gxzy+2*gxxy*gupxy**2*gupyz*gyzy+ &
           2*gxxy*gupxy**2*gupyy*gyyy+2*gxxy*gupxx**2*gupyz*gxzx+ &
           2*gxxy*gupxx**2*gupyy*gxyx+2*gxxx*gupxx*gupxz**2*gzzx+ &
           4*gxxx*gupxy*gupxz**2*gyzz+4*gxxx*gupxy**2*gupxz*gyzy+ &
           2*gxxx*gupxx*gupxy**2*gyyx+4*gxxx*gupxx*gupxz**2*gxzz+ &
           4*gxxx*gupxx**2*gupxz*gxzx+2*gxxx*gupxx**2*gupxz*gxxz+ &
           4*gxyz*gupxz*gupyz**2*gyzz+2*gxyz*gupxy*gupyz**2*gyzy+ &
           2*gxzy*gupxy*gupyy*gupzz*gyzy
      Gamyy = Gamxz+2*gxyy*gupxx*gupyy*gupyz*gxyz+6*gxzz*gupxz*gupyz*gupzz*gyzz+ &
           4*gxzy*gupxz*gupyz*gupzz*gzzz+gupyy**3*gyyy**2+ &
           2*gxzy*gupxy*gupyz*gupzz*gzzy+2*gxzy*gupxx*gupyz*gupzz*gzzx+ &
           2*gxyz*gupxx*gupyz*gupzz*gxzz+2*gxzy*gupxx*gupyz*gupzz*gxzz+ &
           2*gyzy*gupxy*gupyz*gupzz*gzzx+2*gyzy*gupxz*gupyy*gupyz*gxzy+ &
           6*gyzy*gupyy*gupyz*gupzz*gyzz+4*gyzx*gupxz*gupyy*gupyz*gyzy+ &
           4*gyzx*gupxy*gupyz*gupzz*gyzz+2*gxxy*gupxx*gupxy*gupyy*gxyy+ &
           4*gyzx*gupxz*gupyz*gupzz*gzzz+2*gyzx*gupxy*gupyy*gupzz*gyzy+ &
           2*gyyz*gupyy*gupyz*gupzz*gzzy+2*gyyz*gupxy*gupyz*gupzz*gzzx
      Gamxx = Gamyy+2*gyyz*gupyy*gupyz*gupzz*gyzz+2*gyyz*gupxy*gupyy*gupzz*gyzx+ &
           2*gyyz*gupxy*gupyz*gupzz*gxzz+2*gxxy*gupxx*gupxy*gupyz*gyzx+ &
           4*gyyy*gupxy*gupyy*gupyz*gyzx+2*gyyx*gupxy*gupxz*gupyz*gzzx+ &
           2*gxyz*gupxy*gupyz*gupzz*gyzz+2*gxxz*gupxz**2*gupzz*gzzz+ &
           2*gxxz*gupxz**2*gupyz*gyzz+2*gxxz*gupxy*gupxz**2*gxzy+ &
           2*gxxz*gupxx*gupxz**2*gxzx+2*gxxz*gupxy**2*gupyz*gyyy+ &
           2*gxxz*gupxx**2*gupzz*gxzx+2*gxxy*gupxz**2*gupyz*gzzz+ &
           2*gxxy*gupxz**2*gupyy*gyzz+2*gxxy*gupxy*gupxz**2*gxzz+ &
           2*gzzx*gupxz*gupyz*gupzz*gzzy+2*gyzz*gupxz*gupyz*gupzz*gzzx+ &
           2*gxzx*gupxx*gupxz*gupzz*gzzx+2*gyzx*gupxz*gupyy*gupzz*gyzz
      Gamyy = Gamxx+gupzz**3*gzzz**2+2*gxzz*gupxy*gupxz*gupzz*gyzx+ &
           6*gxzx*gupxy*gupxz*gupyz*gyzy+2*gxxy*gupxy*gupxz*gupyz*gzzy+ &
           4*gxzz*gupxy*gupyz**2*gyyy+2*gxzy*gupxz*gupyz**2*gyzz+ &
           2*gxzy*gupxz**2*gupyz*gxzz+2*gxzy*gupxz**2*gupyy*gxyz+ &
           2*gupxy*gupxz*gupyz*gxzy**2+4*gxzx*gupxz**2*gupzz*gzzz+ &
           2*gxzx*gupxz**2*gupyz*gyzz+2*gxyz*gupxy*gupxz*gupzz*gzzx+ &
           2*gxyz*gupxz*gupyy*gupzz*gzzy+2*gxyx*gupxx*gupxz*gupyy*gxyz+ &
           2*gxzz*gupxz*gupyz**2*gyyz+2*gxxy*gupxx*gupxy*gupxz*gxzx+ &
           2*gyyx*gupxy**2*gupxz*gxzx
      Gamxz = Gamyy+2*gxyx*gupxy*gupxz*gupyz*gzzy+2*gyzy*gupyy*gupyz*gupzz*gzzy+ &
           2*gxyx*gupxx*gupxz*gupyy*gyzx+2*gyyx*gupxy*gupyz**2*gyzz+ &
           2*gyyx*gupxy**2*gupyz*gyzx+2*gyyx*gupxz*gupyz**2*gzzz+ &
           2*gyyx*gupxy*gupyy**2*gyyy+2*gxyz*gupxy**2*gupzz*gyzx+ &
           2*gxyz*gupxy**2*gupyz*gyyx+2*gxyy*gupxy*gupyz**2*gyzz+ &
           2*gxyy*gupxy**2*gupyz*gyzx+2*gxyy*gupxy**2*gupyy*gyyx+ &
           2*gxyx*gupxy*gupxz**2*gzzx+2*gxyx*gupxy**2*gupyz*gyyz+ &
           4*gxzz*gupxz*gupzz**2*gzzz+2*gxzz*gupxy*gupzz**2*gzzy+ &
           2*gxzz*gupxx*gupzz**2*gzzx+6*gxyx*gupxx*gupxy*gupxz*gxzx+ &
           2*gxyz*gupxy*gupxz*gupyz*gyzx
      Gamyy = Gamxz+2*gyyx*gupxz*gupyy**2*gyzy+2*gyyx*gupxz*gupyy*gupyz*gzzy+ &
           2*gxxz*gupxx*gupxy*gupyz*gxyy+2*gyzx*gupxz**2*gupyy*gxzy+ &
           4*gyzx*gupxy*gupxz**2*gxzx+2*gyzx*gupxz*gupyz**2*gyzz+ &
           2*gyzx*gupxz**2*gupyz*gxzz+2*gupxy*gupxz*gupyz*gyzx**2+ &
           2*gyyz*gupyz**2*gupzz*gzzz+2*gyyz*gupyy*gupyz**2*gyzy+ &
           2*gyyz*gupxy*gupyz**2*gyzx+2*gyyz*gupyy**2*gupzz*gyzy+ &
           2*gyyz*gupxy**2*gupzz*gxzx+2*gyyy*gupyy*gupyz**2*gzzy+ &
           2*gyyy*gupxy*gupyz**2*gzzx+4*gyyy*gupyy*gupyz**2*gyzz+ &
           4*gyyy*gupyy**2*gupyz*gyzy+2*gyyy*gupyy**2*gupyz*gyyz
      Gamxy = Gamyy+2*gxyz*gupxz*gupyy*gupyz*gyzy+2*gxyz*gupxx*gupyy*gupyz*gyyx+ &
           2*gzzx*gupxz*gupzz**2*gzzz+2*gxzy*gupxy*gupxz*gupyz*gyzx+ &
           2*gyzz*gupyz**2*gupzz*gzzy+2*gyzy*gupxz*gupyz**2*gzzx+ &
           2*gyzx*gupxz*gupyz**2*gzzy+2*gyzx*gupxz**2*gupyz*gzzx+ &
           2*gxzz*gupxz**2*gupzz*gzzx+2*gxzz*gupxy*gupzz**2*gyzz+ &
           2*gxzy*gupxz*gupyz**2*gzzy+2*gxzy*gupxz**2*gupyz*gzzx+ &
           2*gxzx*gupxz**2*gupyz*gzzy+2*gyzz*gupyy*gupzz**2*gzzy+ &
           2*gyzz*gupxy*gupzz**2*gzzx+4*gyzy*gupyz**2*gupzz*gzzz+ &
           2*gyzy*gupxy*gupyz**2*gyzx+2*gyzy*gupxz*gupyz**2*gxzz+ &
           2*gxzy*gupxy*gupyz*gupzz*gyzz+2*gxyx*gupxx*gupxy*gupyz*gxzy
      Gamyy = Gamxy+gupxx**3*gxxx**2+2*gzzy*gupyz*gupzz**2*gzzz+ &
           6*gxyx*gupxx*gupxy*gupyy*gxyy+2*gxzz*gupxz*gupyz* gupzz*gzzy+ &
           6*gxyx*gupxy*gupxz*gupyz*gyzz+2*gxzx*gupxy*gupxz**2*gxzy+ &
           2*gxyx*gupxx*gupxy*gupyy*gyyx+2*gxyx*gupxx*gupxz*gupyz*gzzx+ &
           2*gxyx*gupxx*gupxy*gupxz*gxxz+4*gxyx*gupxx**2*gupxy*gxxx+ &
           2*gxyx*gupxy*gupxz*gupyy*gyyz+6*gxyy*gupxy*gupyy*gupyz*gyzy+ &
           2*gxyx*gupxx*gupxy*gupyz*gyzx+6*gxyy*gupxz*gupyy*gupyz*gyzz+ &
           4*gxyz*gupxx*gupxy*gupzz*gxzx+2*gxyz*gupxy*gupxz*gupzz*gxzz+ &
           4*gxyx*gupxy**2*gupyy*gyyy+2*gxyz*gupxz*gupyy*gupyz*gyyz
      Gamxz = Gamyy+4*gxyz*gupxy*gupyy*gupyz*gyyy+2*gxyx*gupxz**2*gupyy*gyzz+ &
           2*gxyz*gupxz*gupyy*gupzz*gyzz+2*gxyx*gupxy*gupxz**2*gxzz+ &
           4*gxyz*gupxy*gupyy*gupzz*gyzy+2*gxzx*gupxy**2*gupzz*gyzy+ &
           2*gxyz*gupxx*gupxz*gupyz*gxzx+4*gxyx*gupxz**2*gupyz*gzzz+ &
           4*gxzx*gupxy**2*gupyz*gyyy+2*gyyz*gupxy*gupyy*gupzz*gxzy+ &
           2*gxyz*gupxy*gupxz*gupyz*gxzy+2*gxyz*gupxx*gupyz*gupzz*gzzx+ &
           4*gxyy*gupxy*gupyy**2*gyyy+2*gxyy*gupxx*gupyy**2*gyyx+ &
           2*gxyy*gupxx*gupyz**2*gzzx+2*gxyz*gupxy*gupyz*gupzz*gzzy+ &
           2*gxyy*gupxz*gupyy**2*gyyz+4*gxyz*gupxz*gupyz*gupzz*gzzz+ &
           2*gxxy*gupxx*gupxz*gupyy*gxyz
      Gamyy = Gamxz+2*gxzx*gupxy*gupxz**2*gxyz+2*gxxy*gupxy*gupxz*gupyy*gyzy+ &
           4*gxxx*gupxx*gupxy*gupxz*gxzy+2*gxxy*gupxy*gupxz*gupyy*gyyz+ &
           2*gxxy*gupxx*gupxz*gupyy*gyzx+2*gxxy*gupxx*gupxz*gupyz*gzzx+ &
           2*gxzx*gupxy**2*gupxz*gxyy+2*gxxy*gupxx*gupxy*gupyz*gxzy+ &
           2*gxyz*gupxy*gupxz**2*gxxz+2*gxxy*gupxx*gupxy*gupyy*gyyx+ &
           2*gxyz*gupxx*gupyz**2*gyzx+4*gxyz*gupxz**2*gupyz*gxzz+ &
           2*gxxz*gupxx*gupxy*gupzz*gxzy+2*gxxz*gupxx*gupxz*gupzz*gxzz+ &
           2*gxxx*gupxx**2*gupxy*gxxy+2*gxxz*gupxx*gupxy*gupyz*gyyx+ &
           2*gxxz*gupxy*gupxz*gupzz*gyzz+2*gxxz*gupxx*gupxy*gupzz*gyzx
      TZ_rhs = Gamyy+2*gxxz*gupxy*gupxz*gupyz*gyyz+2*gxxz*gupxx*gupxz*gupyz*gyzx+ &
           2*gxxz*gupxx*gupxz*gupzz*gzzx+2*gxxz*gupxy*gupxz*gupyz*gyzy+ &
           2*gxzx*gupxx*gupxy*gupzz*gxzy+2*gxxz*gupxy*gupxz*gupzz*gzzy+ &
           6*gxzx*gupxy*gupxz*gupzz*gyzz+2*gxzx*gupxx*gupxy*gupzz*gyzx+ &
           2*gxzx*gupxx*gupxy*gupyz*gyyx+6*gxzx*gupxx*gupxz*gupzz*gxzz+ &
           2*gxxx*gupxy**2*gupxz*gyyz+2*gxzx*gupxy*gupxz*gupzz*gzzy+ &
           2*gxzx*gupxx*gupxz*gupyz*gyzx+2*gxxx*gupxy*gupxz**2*gzzy+ &
           4*gxzy*gupxy*gupyz**2*gyzy+2*gxzy*gupxx*gupyz**2*gyzx+ &
           2*gxzz*gupxx*gupyz**2*gyyx+4*gxyx*gupxy**2*gupxz*gyzx+ &
           2*gxyx*gupxz**2*gupyy*gzzy+2*gxyy*gupxx*gupyz**2*gxzz

! Gami_,j will be kept till the end of this routine
  call fderivs(ex,Gamx,Gamxx,Gamxy,Gamxz,X,Y,Z,ANTI,SYM ,SYM ,Symmetry,Lev)
  call fderivs(ex,Gamy,Gamyx,Gamyy,Gamyz,X,Y,Z,SYM ,ANTI,SYM ,Symmetry,Lev)
  call fderivs(ex,Gamz,Gamzx,Gamzy,Gamzz,X,Y,Z,SYM ,SYM ,ANTI,Symmetry,Lev)

  TZ_rhs = chix*Gmxcon+chiy*Gmycon+chiz*Gmzcon &
          +chin1*(Gamxx+Gamyy+Gamzz -          &
          (TWO*(gupxz*gupyz*gyzxz+gupxx*gupyy*gxyxy+gupxy*gupyz*gxzyy+ &
                gupxx*gupxy*gxxxy+gupxx*gupxz*gxxxz+gupxx*gupxy*gxyxx+ &
                gupxx*gupyz*gxyxz+gupxx*gupxz*gxzxx+gupxx*gupyz*gxzxy+ &
                gupxx*gupzz*gxzxz+gupxy*gupxz*gxxyz+gupxy*gupyy*gxyyy+ &
                gupxy*gupyz*gxyyz+gupxy*gupxz*gxzxy+gupxy*gupzz*gxzyz+ &
                gupxy*gupxz*gxyxz+gupxz*gupyy*gxyyz+gupxz*gupyz*gxyzz+ &
                gupxz*gupyz*gxzyz+gupxz*gupzz*gxzzz+gupxy*gupyy*gyyxy+ &
                gupxy*gupyz*gyyxz+gupxy*gupxz*gyzxx+gupxy*gupyz*gyzxy+ &
                gupxy*gupzz*gyzxz+gupyy*gupyz*gyyyz+gupxz*gupyy*gyzxy+ &
                gupyy*gupyz*gyzyy+gupyy*gupzz*gyzyz+gupyz*gupzz*gyzzz+ &
                gupxz*gupyz*gzzxy+gupxz*gupzz*gzzxz+gupyz*gupzz*gzzyz+ &
                gupxy*gupxy*gxyxy+gupxz*gupxz*gxzxz+gupyz*gupyz*gyzyz) &
               +gupxx*gupxx*gxxxx+gupxy*gupxy*gxxyy+gupxz*gupxz*gxxzz+ &
                gupxy*gupxy*gyyxx+gupyy*gupyy*gyyyy+gupyz*gupyz*gyyzz+ &
                gupxz*gupxz*gzzxx+gupyz*gupyz*gzzyy+gupzz*gupzz*gzzzz)+&
               (gxx*Gamxa*Gamxa+gyy*Gamya*Gamya+gzz*Gamza*Gamza       +&
           TWO*(gxy*Gamxa*Gamya+gxz*Gamxa*Gamza+gyz*Gamya*Gamza)) + TZ_rhs)

! Raise indices of \tilde A_{ij} and store in R_ij

  Rxx =    gupxx * gupxx * Axx + gupxy * gupxy * Ayy + gupxz * gupxz * Azz + &
      TWO*(gupxx * gupxy * Axy + gupxx * gupxz * Axz + gupxy * gupxz * Ayz)

  Ryy =    gupxy * gupxy * Axx + gupyy * gupyy * Ayy + gupyz * gupyz * Azz + &
      TWO*(gupxy * gupyy * Axy + gupxy * gupyz * Axz + gupyy * gupyz * Ayz)

  Rzz =    gupxz * gupxz * Axx + gupyz * gupyz * Ayy + gupzz * gupzz * Azz + &
      TWO*(gupxz * gupyz * Axy + gupxz * gupzz * Axz + gupyz * gupzz * Ayz)

  Rxy =    gupxx * gupxy * Axx + gupxy * gupyy * Ayy + gupxz * gupyz * Azz + &
          (gupxx * gupyy       + gupxy * gupxy)* Axy                       + &
          (gupxx * gupyz       + gupxz * gupxy)* Axz                       + &
          (gupxy * gupyz       + gupxz * gupyy)* Ayz

  Rxz =    gupxx * gupxz * Axx + gupxy * gupyz * Ayy + gupxz * gupzz * Azz + &
          (gupxx * gupyz       + gupxy * gupxz)* Axy                       + &
          (gupxx * gupzz       + gupxz * gupxz)* Axz                       + &
          (gupxy * gupzz       + gupxz * gupyz)* Ayz

  Ryz =    gupxy * gupxz * Axx + gupyy * gupyz * Ayy + gupyz * gupzz * Azz + &
          (gupxy * gupyz       + gupyy * gupxz)* Axy                       + &
          (gupxy * gupzz       + gupyz * gupxz)* Axz                       + &
          (gupyy * gupzz       + gupyz * gupyz)* Ayz

! Right hand side for Gam^i without shift terms...
! Lap_,i will be kept till the end of this routine
  call fderivs(ex,Lap,Lapx,Lapy,Lapz,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)
! K_,i stored K_,i+TZ_,i/2 indeed, will be kept till the end of this routine  
  call fderivs(ex,trK,Kx,Ky,Kz,X,Y,Z,SYM,SYM,SYM,symmetry,Lev)
  call fderivs(ex,TZ,fxx,fxy,fxz,X,Y,Z,SYM,SYM,SYM,symmetry,Lev)

  Kx = Kx + fxx/TWO
  Ky = Ky + fxy/TWO
  Kz = Kz + fxz/TWO

   Gamx_rhs = - TWO * (   Lapx * Rxx +   Lapy * Rxy +   Lapz * Rxz ) + &
        TWO * alpn1 * (                                                &
        -F3o2/chin1 * (   chix * Rxx +   chiy * Rxy +   chiz * Rxz ) - &
              gupxx * (   F2o3 * Kx  +  EIGHT * PI * Sx            ) - &
              gupxy * (   F2o3 * Ky  +  EIGHT * PI * Sy            ) - &
              gupxz * (   F2o3 * Kz  +  EIGHT * PI * Sz            ) + &
                        Gamxxx * Rxx + Gamxyy * Ryy + Gamxzz * Rzz   + &
                TWO * ( Gamxxy * Rxy + Gamxxz * Rxz + Gamxyz * Ryz ) )

   Gamy_rhs = - TWO * (   Lapx * Rxy +   Lapy * Ryy +   Lapz * Ryz ) + &
        TWO * alpn1 * (                                                &
        -F3o2/chin1 * (   chix * Rxy +  chiy * Ryy +    chiz * Ryz ) - &
              gupxy * (   F2o3 * Kx  +  EIGHT * PI * Sx            ) - &
              gupyy * (   F2o3 * Ky  +  EIGHT * PI * Sy            ) - &
              gupyz * (   F2o3 * Kz  +  EIGHT * PI * Sz            ) + &
                        Gamyxx * Rxx + Gamyyy * Ryy + Gamyzz * Rzz   + &
                TWO * ( Gamyxy * Rxy + Gamyxz * Rxz + Gamyyz * Ryz ) )

   Gamz_rhs = - TWO * (   Lapx * Rxz +   Lapy * Ryz +   Lapz * Rzz ) + &
        TWO * alpn1 * (                                                &
        -F3o2/chin1 * (   chix * Rxz +  chiy * Ryz +    chiz * Rzz ) - &
              gupxz * (   F2o3 * Kx  +  EIGHT * PI * Sx            ) - &
              gupyz * (   F2o3 * Ky  +  EIGHT * PI * Sy            ) - &
              gupzz * (   F2o3 * Kz  +  EIGHT * PI * Sz            ) + &
                        Gamzxx * Rxx + Gamzyy * Ryy + Gamzzz * Rzz   + &
                TWO * ( Gamzxy * Rxy + Gamzxz * Rxz + Gamzyz * Ryz ) )
         
    call fdderivs(ex,betax,gxxx,gxyx,gxzx,gyyx,gyzx,gzzx,&
                X,Y,Z,ANTI,SYM, SYM ,Symmetry,Lev)
    call fdderivs(ex,betay,gxxy,gxyy,gxzy,gyyy,gyzy,gzzy,&
                X,Y,Z,SYM ,ANTI,SYM ,Symmetry,Lev)
    call fdderivs(ex,betaz,gxxz,gxyz,gxzz,gyyz,gyzz,gzzz,&
                X,Y,Z,SYM ,SYM, ANTI,Symmetry,Lev)

  fxx = gxxx + gxyy + gxzz
  fxy = gxyx + gyyy + gyzz
  fxz = gxzx + gyzy + gzzz

  Gamx_rhs =               Gamx_rhs +  F2o3 *  Gamxa * div_beta        - &
                     Gamxa * betaxx - Gamya * betaxy - Gamza * betaxz  + &
             F1o3 * (gupxx * fxx    + gupxy * fxy    + gupxz * fxz    ) + &
                     gupxx * gxxx   + gupyy * gyyx   + gupzz * gzzx    + &
              TWO * (gupxy * gxyx   + gupxz * gxzx   + gupyz * gyzx  )

  Gamy_rhs =               Gamy_rhs +  F2o3 *  Gamya * div_beta        - &
                     Gamxa * betayx - Gamya * betayy - Gamza * betayz  + &
             F1o3 * (gupxy * fxx    + gupyy * fxy    + gupyz * fxz    ) + &
                     gupxx * gxxy   + gupyy * gyyy   + gupzz * gzzy    + &
              TWO * (gupxy * gxyy   + gupxz * gxzy   + gupyz * gyzy  )

  Gamz_rhs =               Gamz_rhs +  F2o3 *  Gamza * div_beta        - &
                     Gamxa * betazx - Gamya * betazy - Gamza * betazz  + &
             F1o3 * (gupxz * fxx    + gupyz * fxy    + gupzz * fxz    ) + &
                     gupxx * gxxz   + gupyy * gyyz   + gupzz * gzzz    + &
              TWO * (gupxy * gxyz   + gupxz * gxzz   + gupyz * gyzz  )    !rhs for Gam^i

!first kind of connection stored in gij,k
  gxxx = gxx * Gamxxx + gxy * Gamyxx + gxz * Gamzxx
  gxyx = gxx * Gamxxy + gxy * Gamyxy + gxz * Gamzxy
  gxzx = gxx * Gamxxz + gxy * Gamyxz + gxz * Gamzxz
  gyyx = gxx * Gamxyy + gxy * Gamyyy + gxz * Gamzyy
  gyzx = gxx * Gamxyz + gxy * Gamyyz + gxz * Gamzyz
  gzzx = gxx * Gamxzz + gxy * Gamyzz + gxz * Gamzzz

  gxxy = gxy * Gamxxx + gyy * Gamyxx + gyz * Gamzxx
  gxyy = gxy * Gamxxy + gyy * Gamyxy + gyz * Gamzxy
  gxzy = gxy * Gamxxz + gyy * Gamyxz + gyz * Gamzxz
  gyyy = gxy * Gamxyy + gyy * Gamyyy + gyz * Gamzyy
  gyzy = gxy * Gamxyz + gyy * Gamyyz + gyz * Gamzyz
  gzzy = gxy * Gamxzz + gyy * Gamyzz + gyz * Gamzzz

  gxxz = gxz * Gamxxx + gyz * Gamyxx + gzz * Gamzxx
  gxyz = gxz * Gamxxy + gyz * Gamyxy + gzz * Gamzxy
  gxzz = gxz * Gamxxz + gyz * Gamyxz + gzz * Gamzxz
  gyyz = gxz * Gamxyy + gyz * Gamyyy + gzz * Gamzyy
  gyzz = gxz * Gamxyz + gyz * Gamyyz + gzz * Gamzyz
  gzzz = gxz * Gamxzz + gyz * Gamyzz + gzz * Gamzzz

!compute Ricci tensor for tilted metric
   Rxx =   gupxx * gxxxx + gupyy * gxxyy + gupzz * gxxzz + &
         ( gupxy * gxxxy + gupxz * gxxxz + gupyz * gxxyz ) * TWO

   Ryy =   gupxx * gyyxx + gupyy * gyyyy + gupzz * gyyzz + &
         ( gupxy * gyyxy + gupxz * gyyxz + gupyz * gyyyz ) * TWO

   Rzz =   gupxx * gzzxx + gupyy * gzzyy + gupzz * gzzzz + &
         ( gupxy * gzzxy + gupxz * gzzxz + gupyz * gzzyz ) * TWO

   Rxy =   gupxx * gxyxx + gupyy * gxyyy + gupzz * gxyzz + &
         ( gupxy * gxyxy + gupxz * gxyxz + gupyz * gxyyz ) * TWO

   Rxz =   gupxx * gxzxx + gupyy * gxzyy + gupzz * gxzzz + &
         ( gupxy * gxzxy + gupxz * gxzxz + gupyz * gxzyz ) * TWO

   Ryz =   gupxx * gyzxx + gupyy * gyzyy + gupzz * gyzzz + &
         ( gupxy * gyzxy + gupxz * gyzxz + gupyz * gyzyz ) * TWO

  Rxx =     - HALF * Rxx                                   + &
               gxx * Gamxx+ gxy * Gamyx   +    gxz * Gamzx + &
             Gamxa * gxxx +  Gamya * gxyx +  Gamza * gxzx  + &
   gupxx *(                                                  &
       TWO*(Gamxxx * gxxx + Gamyxx * gxyx + Gamzxx * gxzx) + &
            Gamxxx * gxxx + Gamyxx * gxxy + Gamzxx * gxxz )+ &
   gupxy *(                                                  &
       TWO*(Gamxxx * gxyx + Gamyxx * gyyx + Gamzxx * gyzx  + &
            Gamxxy * gxxx + Gamyxy * gxyx + Gamzxy * gxzx) + &
            Gamxxy * gxxx + Gamyxy * gxxy + Gamzxy * gxxz  + &
            Gamxxx * gxyx + Gamyxx * gxyy + Gamzxx * gxyz )+ &
   gupxz *(                                                  &
       TWO*(Gamxxx * gxzx + Gamyxx * gyzx + Gamzxx * gzzx  + &
            Gamxxz * gxxx + Gamyxz * gxyx + Gamzxz * gxzx) + &
            Gamxxz * gxxx + Gamyxz * gxxy + Gamzxz * gxxz  + &
            Gamxxx * gxzx + Gamyxx * gxzy + Gamzxx * gxzz )+ &
   gupyy *(                                                  &
       TWO*(Gamxxy * gxyx + Gamyxy * gyyx + Gamzxy * gyzx) + &
            Gamxxy * gxyx + Gamyxy * gxyy + Gamzxy * gxyz )+ &
   gupyz *(                                                  &
       TWO*(Gamxxy * gxzx + Gamyxy * gyzx + Gamzxy * gzzx  + &
            Gamxxz * gxyx + Gamyxz * gyyx + Gamzxz * gyzx) + &
            Gamxxz * gxyx + Gamyxz * gxyy + Gamzxz * gxyz  + &
            Gamxxy * gxzx + Gamyxy * gxzy + Gamzxy * gxzz )+ &
   gupzz *(                                                  &
       TWO*(Gamxxz * gxzx + Gamyxz * gyzx + Gamzxz * gzzx) + &
            Gamxxz * gxzx + Gamyxz * gxzy + Gamzxz * gxzz )

  Ryy =     - HALF * Ryy                                   + &
               gxy * Gamxy+  gyy * Gamyy  +  gyz * Gamzy   + &
             Gamxa * gxyy +  Gamya * gyyy +  Gamza * gyzy  + &
   gupxx *(                                                  &
       TWO*(Gamxxy * gxxy + Gamyxy * gxyy + Gamzxy * gxzy) + &
            Gamxxy * gxyx + Gamyxy * gxyy + Gamzxy * gxyz )+ &
   gupxy *(                                                  &
       TWO*(Gamxxy * gxyy + Gamyxy * gyyy + Gamzxy * gyzy  + &
            Gamxyy * gxxy + Gamyyy * gxyy + Gamzyy * gxzy) + &
            Gamxyy * gxyx + Gamyyy * gxyy + Gamzyy * gxyz  + &
            Gamxxy * gyyx + Gamyxy * gyyy + Gamzxy * gyyz )+ &
   gupxz *(                                                  &
       TWO*(Gamxxy * gxzy + Gamyxy * gyzy + Gamzxy * gzzy  + &
            Gamxyz * gxxy + Gamyyz * gxyy + Gamzyz * gxzy) + &
            Gamxyz * gxyx + Gamyyz * gxyy + Gamzyz * gxyz  + &
            Gamxxy * gyzx + Gamyxy * gyzy + Gamzxy * gyzz )+ &
   gupyy *(                                                  &
       TWO*(Gamxyy * gxyy + Gamyyy * gyyy + Gamzyy * gyzy) + &
            Gamxyy * gyyx + Gamyyy * gyyy + Gamzyy * gyyz )+ &
   gupyz *(                                                  &
       TWO*(Gamxyy * gxzy + Gamyyy * gyzy + Gamzyy * gzzy  + &
            Gamxyz * gxyy + Gamyyz * gyyy + Gamzyz * gyzy) + &
            Gamxyz * gyyx + Gamyyz * gyyy + Gamzyz * gyyz  + &
            Gamxyy * gyzx + Gamyyy * gyzy + Gamzyy * gyzz )+ &
   gupzz *(                                                  &
       TWO*(Gamxyz * gxzy + Gamyyz * gyzy + Gamzyz * gzzy) + &
            Gamxyz * gyzx + Gamyyz * gyzy + Gamzyz * gyzz )

  Rzz =     - HALF * Rzz                                   + &
               gxz * Gamxz+ gyz * Gamyz  +    gzz * Gamzz  + &
             Gamxa * gxzz +  Gamya * gyzz +  Gamza * gzzz  + &
   gupxx *(                                                  &
       TWO*(Gamxxz * gxxz + Gamyxz * gxyz + Gamzxz * gxzz) + &
            Gamxxz * gxzx + Gamyxz * gxzy + Gamzxz * gxzz )+ &
   gupxy *(                                                  &
       TWO*(Gamxxz * gxyz + Gamyxz * gyyz + Gamzxz * gyzz  + &
            Gamxyz * gxxz + Gamyyz * gxyz + Gamzyz * gxzz) + &
            Gamxyz * gxzx + Gamyyz * gxzy + Gamzyz * gxzz  + &
            Gamxxz * gyzx + Gamyxz * gyzy + Gamzxz * gyzz )+ &
   gupxz *(                                                  &
       TWO*(Gamxxz * gxzz + Gamyxz * gyzz + Gamzxz * gzzz  + &
            Gamxzz * gxxz + Gamyzz * gxyz + Gamzzz * gxzz) + &
            Gamxzz * gxzx + Gamyzz * gxzy + Gamzzz * gxzz  + &
            Gamxxz * gzzx + Gamyxz * gzzy + Gamzxz * gzzz )+ &
   gupyy *(                                                  &
       TWO*(Gamxyz * gxyz + Gamyyz * gyyz + Gamzyz * gyzz) + &
            Gamxyz * gyzx + Gamyyz * gyzy + Gamzyz * gyzz )+ &
   gupyz *(                                                  &
       TWO*(Gamxyz * gxzz + Gamyyz * gyzz + Gamzyz * gzzz  + &
            Gamxzz * gxyz + Gamyzz * gyyz + Gamzzz * gyzz) + &
            Gamxzz * gyzx + Gamyzz * gyzy + Gamzzz * gyzz  + &
            Gamxyz * gzzx + Gamyyz * gzzy + Gamzyz * gzzz )+ &
   gupzz *(                                                  &
       TWO*(Gamxzz * gxzz + Gamyzz * gyzz + Gamzzz * gzzz) + &
            Gamxzz * gzzx + Gamyzz * gzzy + Gamzzz * gzzz )

  Rxy = HALF*(     - Rxy                                   + &
               gxx * Gamxy +    gxy * Gamyy + gxz * Gamzy  + &
               gxy * Gamxx +    gyy * Gamyx + gyz * Gamzx  + &
             Gamxa * gxyx +  Gamya * gyyx +  Gamza * gyzx  + &
             Gamxa * gxxy +  Gamya * gxyy +  Gamza * gxzy )+ &
   gupxx *(                                                  &
            Gamxxx * gxxy + Gamyxx * gxyy + Gamzxx * gxzy  + &
            Gamxxy * gxxx + Gamyxy * gxyx + Gamzxy * gxzx  + &
            Gamxxx * gxyx + Gamyxx * gxyy + Gamzxx * gxyz )+ &
   gupxy *(                                                  &
            Gamxxx * gxyy + Gamyxx * gyyy + Gamzxx * gyzy  + &
            Gamxxy * gxyx + Gamyxy * gyyx + Gamzxy * gyzx  + &
            Gamxxy * gxyx + Gamyxy * gxyy + Gamzxy * gxyz  + &
            Gamxxy * gxxy + Gamyxy * gxyy + Gamzxy * gxzy  + &
            Gamxyy * gxxx + Gamyyy * gxyx + Gamzyy * gxzx  + &
            Gamxxx * gyyx + Gamyxx * gyyy + Gamzxx * gyyz )+ &
   gupxz *(                                                  &
            Gamxxx * gxzy + Gamyxx * gyzy + Gamzxx * gzzy  + &
            Gamxxy * gxzx + Gamyxy * gyzx + Gamzxy * gzzx  + &
            Gamxxz * gxyx + Gamyxz * gxyy + Gamzxz * gxyz  + &
            Gamxxz * gxxy + Gamyxz * gxyy + Gamzxz * gxzy  + &
            Gamxyz * gxxx + Gamyyz * gxyx + Gamzyz * gxzx  + &
            Gamxxx * gyzx + Gamyxx * gyzy + Gamzxx * gyzz )+ &
   gupyy *(                                                  &
            Gamxxy * gxyy + Gamyxy * gyyy + Gamzxy * gyzy  + &
            Gamxyy * gxyx + Gamyyy * gyyx + Gamzyy * gyzx  + &
            Gamxxy * gyyx + Gamyxy * gyyy + Gamzxy * gyyz )+ &
   gupyz *(                                                  &
            Gamxxy * gxzy + Gamyxy * gyzy + Gamzxy * gzzy  + &
            Gamxyy * gxzx + Gamyyy * gyzx + Gamzyy * gzzx  + &
            Gamxxz * gyyx + Gamyxz * gyyy + Gamzxz * gyyz  + &
            Gamxxz * gxyy + Gamyxz * gyyy + Gamzxz * gyzy  + &
            Gamxyz * gxyx + Gamyyz * gyyx + Gamzyz * gyzx  + &
            Gamxxy * gyzx + Gamyxy * gyzy + Gamzxy * gyzz )+ &
   gupzz *(                                                  &
            Gamxxz * gxzy + Gamyxz * gyzy + Gamzxz * gzzy  + &
            Gamxyz * gxzx + Gamyyz * gyzx + Gamzyz * gzzx  + &
            Gamxxz * gyzx + Gamyxz * gyzy + Gamzxz * gyzz )

  Rxz = HALF*(     - Rxz                                   + &
               gxx * Gamxz +  gxy * Gamyz + gxz * Gamzz    + &
               gxz * Gamxx +  gyz * Gamyx + gzz * Gamzx    + &
             Gamxa * gxzx +  Gamya * gyzx +  Gamza * gzzx  + &
             Gamxa * gxxz +  Gamya * gxyz +  Gamza * gxzz )+ &
   gupxx *(                                                  &
            Gamxxx * gxxz + Gamyxx * gxyz + Gamzxx * gxzz  + &
            Gamxxz * gxxx + Gamyxz * gxyx + Gamzxz * gxzx  + &
            Gamxxx * gxzx + Gamyxx * gxzy + Gamzxx * gxzz )+ &
   gupxy *(                                                  &
            Gamxxx * gxyz + Gamyxx * gyyz + Gamzxx * gyzz  + &
            Gamxxz * gxyx + Gamyxz * gyyx + Gamzxz * gyzx  + &
            Gamxxy * gxzx + Gamyxy * gxzy + Gamzxy * gxzz  + &
            Gamxxy * gxxz + Gamyxy * gxyz + Gamzxy * gxzz  + &
            Gamxyz * gxxx + Gamyyz * gxyx + Gamzyz * gxzx  + &
            Gamxxx * gyzx + Gamyxx * gyzy + Gamzxx * gyzz )+ &
   gupxz *(                                                  &
            Gamxxx * gxzz + Gamyxx * gyzz + Gamzxx * gzzz  + &
            Gamxxz * gxzx + Gamyxz * gyzx + Gamzxz * gzzx  + &
            Gamxxz * gxzx + Gamyxz * gxzy + Gamzxz * gxzz  + &
            Gamxxz * gxxz + Gamyxz * gxyz + Gamzxz * gxzz  + &
            Gamxzz * gxxx + Gamyzz * gxyx + Gamzzz * gxzx  + &
            Gamxxx * gzzx + Gamyxx * gzzy + Gamzxx * gzzz )+ &
   gupyy *(                                                  &
            Gamxxy * gxyz + Gamyxy * gyyz + Gamzxy * gyzz  + &
            Gamxyz * gxyx + Gamyyz * gyyx + Gamzyz * gyzx  + &
            Gamxxy * gyzx + Gamyxy * gyzy + Gamzxy * gyzz )+ &
   gupyz *(                                                  &
            Gamxxy * gxzz + Gamyxy * gyzz + Gamzxy * gzzz  + &
            Gamxyz * gxzx + Gamyyz * gyzx + Gamzyz * gzzx  + &
            Gamxxz * gyzx + Gamyxz * gyzy + Gamzxz * gyzz  + &
            Gamxxz * gxyz + Gamyxz * gyyz + Gamzxz * gyzz  + &
            Gamxzz * gxyx + Gamyzz * gyyx + Gamzzz * gyzx  + &
            Gamxxy * gzzx + Gamyxy * gzzy + Gamzxy * gzzz )+ &
   gupzz *(                                                  &
            Gamxxz * gxzz + Gamyxz * gyzz + Gamzxz * gzzz  + &
            Gamxzz * gxzx + Gamyzz * gyzx + Gamzzz * gzzx  + &
            Gamxxz * gzzx + Gamyxz * gzzy + Gamzxz * gzzz )

  Ryz = HALF*(     - Ryz                                   + &
               gxy * Gamxz + gyy * Gamyz + gyz * Gamzz     + &
               gxz * Gamxy + gyz * Gamyy + gzz * Gamzy     + &
             Gamxa * gxzy +  Gamya * gyzy +  Gamza * gzzy  + &
             Gamxa * gxyz +  Gamya * gyyz +  Gamza * gyzz )+ &
   gupxx *(                                                  &
            Gamxxy * gxxz + Gamyxy * gxyz + Gamzxy * gxzz  + &
            Gamxxz * gxxy + Gamyxz * gxyy + Gamzxz * gxzy  + &
            Gamxxy * gxzx + Gamyxy * gxzy + Gamzxy * gxzz )+ &
   gupxy *(                                                  &
            Gamxxy * gxyz + Gamyxy * gyyz + Gamzxy * gyzz  + &
            Gamxxz * gxyy + Gamyxz * gyyy + Gamzxz * gyzy  + &
            Gamxyy * gxzx + Gamyyy * gxzy + Gamzyy * gxzz  + &
            Gamxyy * gxxz + Gamyyy * gxyz + Gamzyy * gxzz  + &
            Gamxyz * gxxy + Gamyyz * gxyy + Gamzyz * gxzy  + &
            Gamxxy * gyzx + Gamyxy * gyzy + Gamzxy * gyzz )+ &
   gupxz *(                                                  &
            Gamxxy * gxzz + Gamyxy * gyzz + Gamzxy * gzzz  + &
            Gamxxz * gxzy + Gamyxz * gyzy + Gamzxz * gzzy  + &
            Gamxyz * gxzx + Gamyyz * gxzy + Gamzyz * gxzz  + &
            Gamxyz * gxxz + Gamyyz * gxyz + Gamzyz * gxzz  + &
            Gamxzz * gxxy + Gamyzz * gxyy + Gamzzz * gxzy  + &
            Gamxxy * gzzx + Gamyxy * gzzy + Gamzxy * gzzz )+ &
   gupyy *(                                                  &
            Gamxyy * gxyz + Gamyyy * gyyz + Gamzyy * gyzz  + &
            Gamxyz * gxyy + Gamyyz * gyyy + Gamzyz * gyzy  + &
            Gamxyy * gyzx + Gamyyy * gyzy + Gamzyy * gyzz )+ &
   gupyz *(                                                  &
            Gamxyy * gxzz + Gamyyy * gyzz + Gamzyy * gzzz  + &
            Gamxyz * gxzy + Gamyyz * gyzy + Gamzyz * gzzy  + &
            Gamxyz * gyzx + Gamyyz * gyzy + Gamzyz * gyzz  + &
            Gamxyz * gxyz + Gamyyz * gyyz + Gamzyz * gyzz  + &
            Gamxzz * gxyy + Gamyzz * gyyy + Gamzzz * gyzy  + &
            Gamxyy * gzzx + Gamyyy * gzzy + Gamzyy * gzzz )+ &
   gupzz *(                                                  &
            Gamxyz * gxzz + Gamyyz * gyzz + Gamzyz * gzzz  + &
            Gamxzz * gxzy + Gamyzz * gyzy + Gamzzz * gzzy  + &
            Gamxyz * gzzx + Gamyyz * gzzy + Gamzyz * gzzz )
!covariant second derivative of chi respect to tilted metric

! Store D^l D_l chi - 3/(2*chi) D^l chi D_l chi in f

  call fdderivs(ex,chi,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)

  fxx = fxx - Gamxxx * chix - Gamyxx * chiy - Gamzxx * chiz
  fxy = fxy - Gamxxy * chix - Gamyxy * chiy - Gamzxy * chiz
  fxz = fxz - Gamxxz * chix - Gamyxz * chiy - Gamzxz * chiz
  fyy = fyy - Gamxyy * chix - Gamyyy * chiy - Gamzyy * chiz
  fyz = fyz - Gamxyz * chix - Gamyyz * chiy - Gamzyz * chiz
  fzz = fzz - Gamxzz * chix - Gamyzz * chiy - Gamzzz * chiz

  f =        gupxx * ( fxx - F3o2/chin1 * chix * chix ) + &
             gupyy * ( fyy - F3o2/chin1 * chiy * chiy ) + &
             gupzz * ( fzz - F3o2/chin1 * chiz * chiz ) + &
       TWO * gupxy * ( fxy - F3o2/chin1 * chix * chiy ) + &
       TWO * gupxz * ( fxz - F3o2/chin1 * chix * chiz ) + &
       TWO * gupyz * ( fyz - F3o2/chin1 * chiy * chiz ) 

! Add chi part to Ricci tensor:

  fxx = Rxx + (fxx - chix*chix/chin1/TWO + gxx * f)/chin1/TWO
  fyy = Ryy + (fyy - chiy*chiy/chin1/TWO + gyy * f)/chin1/TWO
  fzz = Rzz + (fzz - chiz*chiz/chin1/TWO + gzz * f)/chin1/TWO
  fxy = Rxy + (fxy - chix*chiy/chin1/TWO + gxy * f)/chin1/TWO
  fxz = Rxz + (fxz - chix*chiz/chin1/TWO + gxz * f)/chin1/TWO
  fyz = Ryz + (fyz - chiy*chiz/chin1/TWO + gyz * f)/chin1/TWO  
! store R/chi in Hcon
  Hcon =   gupxx * fxx + gupyy * fyy + gupzz * fzz + &
        TWO* ( gupxy * fxy + gupxz * fxz + gupyz * fyz )

  Rxx = fxx
  Ryy = fyy
  Rzz = fzz
  Rxy = fxy
  Rxz = fxz
  Ryz = fyz

  gxxx = (gupxx * chix + gupxy * chiy + gupxz * chiz)/chin1
  gxxy = (gupxy * chix + gupyy * chiy + gupyz * chiz)/chin1
  gxxz = (gupxz * chix + gupyz * chiy + gupzz * chiz)/chin1
! now get physical second kind of connection
  Gamxxx = Gamxxx - ( (chix + chix)/chin1 - gxx * gxxx )*HALF
  Gamyxx = Gamyxx - (                     - gxx * gxxy )*HALF
  Gamzxx = Gamzxx - (                     - gxx * gxxz )*HALF
  Gamxyy = Gamxyy - (                     - gyy * gxxx )*HALF
  Gamyyy = Gamyyy - ( (chiy + chiy)/chin1 - gyy * gxxy )*HALF
  Gamzyy = Gamzyy - (                     - gyy * gxxz )*HALF
  Gamxzz = Gamxzz - (                     - gzz * gxxx )*HALF
  Gamyzz = Gamyzz - (                     - gzz * gxxy )*HALF
  Gamzzz = Gamzzz - ( (chiz + chiz)/chin1 - gzz * gxxz )*HALF
  Gamxxy = Gamxxy - (  chiy        /chin1 - gxy * gxxx )*HALF
  Gamyxy = Gamyxy - (         chix /chin1 - gxy * gxxy )*HALF
  Gamzxy = Gamzxy - (                     - gxy * gxxz )*HALF
  Gamxxz = Gamxxz - (  chiz        /chin1 - gxz * gxxx )*HALF
  Gamyxz = Gamyxz - (                     - gxz * gxxy )*HALF
  Gamzxz = Gamzxz - (         chix /chin1 - gxz * gxxz )*HALF
  Gamxyz = Gamxyz - (                     - gyz * gxxx )*HALF
  Gamyyz = Gamyyz - (  chiz        /chin1 - gyz * gxxy )*HALF
  Gamzyz = Gamzyz - (         chiy /chin1 - gyz * gxxz )*HALF

! covariant second derivatives of the lapse respect to physical metric

   call fdderivs(ex,Lap,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z, &
                SYM,SYM,SYM,symmetry,Lev)

  fxx = fxx - Gamxxx*Lapx - Gamyxx*Lapy - Gamzxx*Lapz
  fyy = fyy - Gamxyy*Lapx - Gamyyy*Lapy - Gamzyy*Lapz
  fzz = fzz - Gamxzz*Lapx - Gamyzz*Lapy - Gamzzz*Lapz
  fxy = fxy - Gamxxy*Lapx - Gamyxy*Lapy - Gamzxy*Lapz
  fxz = fxz - Gamxxz*Lapx - Gamyxz*Lapy - Gamzxz*Lapz
  fyz = fyz - Gamxyz*Lapx - Gamyyz*Lapy - Gamzyz*Lapz

! store D^i D_i Lap in trK_rhs upto chi
  trK_rhs =    gupxx * fxx + gupyy * fyy + gupzz * fzz + &
        TWO* ( gupxy * fxy + gupxz * fxz + gupyz * fyz )
! Add lapse and S_ij parts to Ricci tensor:

  fxx = EIGHT * PI * alpn1 * Sxx + fxx
  fxy = EIGHT * PI * alpn1 * Sxy + fxy
  fxz = EIGHT * PI * alpn1 * Sxz + fxz
  fyy = EIGHT * PI * alpn1 * Syy + fyy
  fyz = EIGHT * PI * alpn1 * Syz + fyz
  fzz = EIGHT * PI * alpn1 * Szz + fzz

! Compute trace-free part (note: chi^-1 and chi cancel!):
  f =        gupxx * fxx + gupyy * fyy + gupzz * fzz + &
      TWO* ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) 

  f = F1o3 * (Hcon*alpn1 - f)

  fxx = alpn1 * Rxx - fxx
  fxy = alpn1 * Rxy - fxy
  fxz = alpn1 * Rxz - fxz
  fyy = alpn1 * Ryy - fyy
  fyz = alpn1 * Ryz - fyz
  fzz = alpn1 * Rzz - fzz

  Axx_rhs = fxx - gxx * f
  Ayy_rhs = fyy - gyy * f
  Azz_rhs = fzz - gzz * f
  Axy_rhs = fxy - gxy * f
  Axz_rhs = fxz - gxz * f
  Ayz_rhs = fyz - gyz * f

! Now: store A_il A^l_j into fij:

  fxx =       gupxx * Axx * Axx + gupyy * Axy * Axy + gupzz * Axz * Axz + &
       TWO * (gupxy * Axx * Axy + gupxz * Axx * Axz + gupyz * Axy * Axz)
  fyy =       gupxx * Axy * Axy + gupyy * Ayy * Ayy + gupzz * Ayz * Ayz + &
       TWO * (gupxy * Axy * Ayy + gupxz * Axy * Ayz + gupyz * Ayy * Ayz)
  fzz =       gupxx * Axz * Axz + gupyy * Ayz * Ayz + gupzz * Azz * Azz + &
       TWO * (gupxy * Axz * Ayz + gupxz * Axz * Azz + gupyz * Ayz * Azz)
  fxy =       gupxx * Axx * Axy + gupyy * Axy * Ayy + gupzz * Axz * Ayz + &
              gupxy *(Axx * Ayy + Axy * Axy)                            + &
              gupxz *(Axx * Ayz + Axz * Axy)                            + &
              gupyz *(Axy * Ayz + Axz * Ayy)
  fxz =       gupxx * Axx * Axz + gupyy * Axy * Ayz + gupzz * Axz * Azz + &
              gupxy *(Axx * Ayz + Axy * Axz)                            + &
              gupxz *(Axx * Azz + Axz * Axz)                            + &
              gupyz *(Axy * Azz + Axz * Ayz)
  fyz =       gupxx * Axy * Axz + gupyy * Ayy * Ayz + gupzz * Ayz * Azz + &
              gupxy *(Axy * Ayz + Ayy * Axz)                            + &
              gupxz *(Axy * Azz + Ayz * Axz)                            + &
              gupyz *(Ayy * Azz + Ayz * Ayz)

  f = chin1
! store D^i D_i Lap in trK_rhs
  trK_rhs = f*trK_rhs
          
  Axx_rhs =           f * Axx_rhs+ alpn1 * (trKd * Axx - TWO * fxx)  + &
           TWO * (  Axx * betaxx +   Axy * betayx +   Axz * betazx ) - &
             F2o3 * Axx * div_beta

  Ayy_rhs =           f * Ayy_rhs+ alpn1 * (trKd * Ayy - TWO * fyy)  + &
           TWO * (  Axy * betaxy +   Ayy * betayy +   Ayz * betazy ) - &
             F2o3 * Ayy * div_beta

  Azz_rhs =           f * Azz_rhs+ alpn1 * (trKd * Azz - TWO * fzz)  + &
           TWO * (  Axz * betaxz +   Ayz * betayz +   Azz * betazz ) - &
             F2o3 * Azz * div_beta

  Axy_rhs =           f * Axy_rhs+ alpn1 *( trKd * Axy  - TWO * fxy )+ &
                    Axx * betaxy                  +   Axz * betazy   + &
                                     Ayy * betayx +   Ayz * betazx   + &
             F1o3 * Axy * div_beta                -   Axy * betazz

  Ayz_rhs =           f * Ayz_rhs+ alpn1 *( trKd * Ayz  - TWO * fyz )+ &
                    Axy * betaxz +   Ayy * betayz                    + &
                    Axz * betaxy                  +   Azz * betazy   + &
             F1o3 * Ayz * div_beta                -   Ayz * betaxx
 
  Axz_rhs =           f * Axz_rhs+ alpn1 *( trKd * Axz  - TWO * fxz )+ &
                    Axx * betaxz +   Axy * betayz                    + &
                                     Ayz * betayx +   Azz * betazx   + &
             F1o3 * Axz * div_beta                -   Axz * betayy      !rhs for Aij

! Compute trace of S_ij

  S =  f * ( gupxx * Sxx + gupyy * Syy + gupzz * Szz + &
     TWO * ( gupxy * Sxy + gupxz * Sxz + gupyz * Syz ) )

  trK_rhs = - trK_rhs + alpn1 *( F1o3 * trKd * trKd + &
                gupxx * fxx + gupyy * fyy + gupzz * fzz   + &
        TWO * ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) + &
       FOUR * PI * ( rho + S ))                                !rhs for trK

!!!!!gauge variable part
  Lap_rhs = -TWO*alpn1*trK

#if (GAUGE == 0)
  betax_rhs = FF*dtSfx
  betay_rhs = FF*dtSfy
  betaz_rhs = FF*dtSfz

  dtSfx_rhs = Gamx_rhs - eta*dtSfx
  dtSfy_rhs = Gamy_rhs - eta*dtSfy
  dtSfz_rhs = Gamz_rhs - eta*dtSfz
#elif  (GAUGE == 1)
  betax_rhs = Gamx - eta*betax
  betay_rhs = Gamy - eta*betay
  betaz_rhs = Gamz - eta*betaz

  dtSfx_rhs = ZEO
  dtSfy_rhs = ZEO
  dtSfz_rhs = ZEO
#endif  
!!!!!Z4 part
! H = trR + 2/3 * trKd^2 - A_ij * A^ij - 16 * PI * rho
! here trR is respect to physical metric

  Hcon = chin1*Hcon + F2o3 * trKd * trKd -(&
       gupxx * ( &
       gupxx * Axx * Axx + gupyy * Axy * Axy + gupzz * Axz * Axz + &
       TWO * (gupxy * Axx * Axy + gupxz * Axx * Axz + gupyz * Axy * Axz) ) + &
       gupyy * ( &
       gupxx * Axy * Axy + gupyy * Ayy * Ayy + gupzz * Ayz * Ayz + &
       TWO * (gupxy * Axy * Ayy + gupxz * Axy * Ayz + gupyz * Ayy * Ayz) ) + &
       gupzz * ( &
       gupxx * Axz * Axz + gupyy * Ayz * Ayz + gupzz * Azz * Azz + &
       TWO * (gupxy * Axz * Ayz + gupxz * Axz * Azz + gupyz * Ayz * Azz) ) + &
       TWO * ( &
       gupxy * ( &
       gupxx * Axx * Axy + gupyy * Axy * Ayy + gupzz * Axz * Ayz + &
       gupxy * (Axx * Ayy + Axy * Axy) + &
       gupxz * (Axx * Ayz + Axz * Axy) + &
       gupyz * (Axy * Ayz + Axz * Ayy) ) + &
       gupxz * ( &
       gupxx * Axx * Axz + gupyy * Axy * Ayz + gupzz * Axz * Azz + &
       gupxy * (Axx * Ayz + Axy * Axz) + &
       gupxz * (Axx * Azz + Axz * Axz) + &
       gupyz * (Axy * Azz + Axz * Ayz) ) + &
       gupyz * ( &
       gupxx * Axy * Axz + gupyy * Ayy * Ayz + gupzz * Ayz * Azz + &
       gupxy * (Axy * Ayz + Ayy * Axz) + &
       gupxz * (Axy * Azz + Ayz * Axz) + &
       gupyz * (Ayy * Azz + Ayz * Ayz) ) ))- F16 * PI * rho
! M_j = gupki*(-1/chi d_k chi*A_ij + D_k A_ij) - 2/3 d_j trK - 8 PI s_j where D respect to physical metric
! store D_i A_jk - 1/chi d_i chi*A_jk in gjk_i
  call fderivs(ex,Axx,gxxx,gxxy,gxxz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,lev)
  call fderivs(ex,Axy,gxyx,gxyy,gxyz,X,Y,Z,ANTI,ANTI,SYM ,Symmetry,lev)
  call fderivs(ex,Axz,gxzx,gxzy,gxzz,X,Y,Z,ANTI,SYM ,ANTI,Symmetry,lev)
  call fderivs(ex,Ayy,gyyx,gyyy,gyyz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,lev)
  call fderivs(ex,Ayz,gyzx,gyzy,gyzz,X,Y,Z,SYM ,ANTI,ANTI,Symmetry,lev)
  call fderivs(ex,Azz,gzzx,gzzy,gzzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,lev)

  gxxx = gxxx - (  Gamxxx * Axx + Gamyxx * Axy + Gamzxx * Axz &
                 + Gamxxx * Axx + Gamyxx * Axy + Gamzxx * Axz) - chix*Axx/chin1
  gxyx = gxyx - (  Gamxxy * Axx + Gamyxy * Axy + Gamzxy * Axz &
                 + Gamxxx * Axy + Gamyxx * Ayy + Gamzxx * Ayz) - chix*Axy/chin1
  gxzx = gxzx - (  Gamxxz * Axx + Gamyxz * Axy + Gamzxz * Axz &
                 + Gamxxx * Axz + Gamyxx * Ayz + Gamzxx * Azz) - chix*Axz/chin1
  gyyx = gyyx - (  Gamxxy * Axy + Gamyxy * Ayy + Gamzxy * Ayz &
                 + Gamxxy * Axy + Gamyxy * Ayy + Gamzxy * Ayz) - chix*Ayy/chin1
  gyzx = gyzx - (  Gamxxz * Axy + Gamyxz * Ayy + Gamzxz * Ayz &
                 + Gamxxy * Axz + Gamyxy * Ayz + Gamzxy * Azz) - chix*Ayz/chin1
  gzzx = gzzx - (  Gamxxz * Axz + Gamyxz * Ayz + Gamzxz * Azz &
                 + Gamxxz * Axz + Gamyxz * Ayz + Gamzxz * Azz) - chix*Azz/chin1
  gxxy = gxxy - (  Gamxxy * Axx + Gamyxy * Axy + Gamzxy * Axz &
                 + Gamxxy * Axx + Gamyxy * Axy + Gamzxy * Axz) - chiy*Axx/chin1
  gxyy = gxyy - (  Gamxyy * Axx + Gamyyy * Axy + Gamzyy * Axz &
                 + Gamxxy * Axy + Gamyxy * Ayy + Gamzxy * Ayz) - chiy*Axy/chin1
  gxzy = gxzy - (  Gamxyz * Axx + Gamyyz * Axy + Gamzyz * Axz &
                 + Gamxxy * Axz + Gamyxy * Ayz + Gamzxy * Azz) - chiy*Axz/chin1
  gyyy = gyyy - (  Gamxyy * Axy + Gamyyy * Ayy + Gamzyy * Ayz &
                 + Gamxyy * Axy + Gamyyy * Ayy + Gamzyy * Ayz) - chiy*Ayy/chin1
  gyzy = gyzy - (  Gamxyz * Axy + Gamyyz * Ayy + Gamzyz * Ayz &
                 + Gamxyy * Axz + Gamyyy * Ayz + Gamzyy * Azz) - chiy*Ayz/chin1
  gzzy = gzzy - (  Gamxyz * Axz + Gamyyz * Ayz + Gamzyz * Azz &
                 + Gamxyz * Axz + Gamyyz * Ayz + Gamzyz * Azz) - chiy*Azz/chin1
  gxxz = gxxz - (  Gamxxz * Axx + Gamyxz * Axy + Gamzxz * Axz &
                 + Gamxxz * Axx + Gamyxz * Axy + Gamzxz * Axz) - chiz*Axx/chin1
  gxyz = gxyz - (  Gamxyz * Axx + Gamyyz * Axy + Gamzyz * Axz &
                 + Gamxxz * Axy + Gamyxz * Ayy + Gamzxz * Ayz) - chiz*Axy/chin1
  gxzz = gxzz - (  Gamxzz * Axx + Gamyzz * Axy + Gamzzz * Axz &
                 + Gamxxz * Axz + Gamyxz * Ayz + Gamzxz * Azz) - chiz*Axz/chin1
  gyyz = gyyz - (  Gamxyz * Axy + Gamyyz * Ayy + Gamzyz * Ayz &
                 + Gamxyz * Axy + Gamyyz * Ayy + Gamzyz * Ayz) - chiz*Ayy/chin1
  gyzz = gyzz - (  Gamxzz * Axy + Gamyzz * Ayy + Gamzzz * Ayz &
                 + Gamxyz * Axz + Gamyyz * Ayz + Gamzyz * Azz) - chiz*Ayz/chin1
  gzzz = gzzz - (  Gamxzz * Axz + Gamyzz * Ayz + Gamzzz * Azz &
                 + Gamxzz * Axz + Gamyzz * Ayz + Gamzzz * Azz) - chiz*Azz/chin1
  Mxcon  = gupxx*gxxx + gupyy*gxyy + gupzz*gxzz &
          +gupxy*gxyx + gupxz*gxzx + gupyz*gxzy &
          +gupxy*gxxy + gupxz*gxxz + gupyz*gxyz
  Mycon  = gupxx*gxyx + gupyy*gyyy + gupzz*gyzz &
          +gupxy*gyyx + gupxz*gyzx + gupyz*gyzy &
          +gupxy*gxyy + gupxz*gxyz + gupyz*gyyz
  Mzcon  = gupxx*gxzx + gupyy*gyzy + gupzz*gzzz &
          +gupxy*gyzx + gupxz*gzzx + gupyz*gzzy &
          +gupxy*gxzy + gupxz*gxzz + gupyz*gyzz
! we have already considered TZ_,i in K_,i here, or to say here Micon =
! Micon+TZ_,i indeed
  Mxcon  = Mxcon - F2o3*Kx - F8*PI*sx
  Mycon  = Mycon - F2o3*Ky - F8*PI*sy
  Mzcon  = Mzcon - F2o3*Kz - F8*PI*sz

  f = TZ_rhs

  TZ_rhs = alpn1*Hcon/TWO
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
    TZ_rhs = TZ_rhs - alpn1*(TWO+kappa2)*kappa1*TZ
    trK_rhs = trK_rhs + alpn1*kappa1*(ONE-kappa2)*TZ
    Gamx_rhs = Gamx_rhs - TWO*alpn1*kappa1*(Gamx-Gamxa)
    Gamy_rhs = Gamy_rhs - TWO*alpn1*kappa1*(Gamy-Gamya)
    Gamz_rhs = Gamz_rhs - TWO*alpn1*kappa1*(Gamz-Gamza)

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

  end function compute_rhs_Z4c
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
