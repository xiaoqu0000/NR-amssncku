

!! note that the potential for scalar field in F(R) gravity
!!       is defined in the file Set_Rho_ADM.f90

#include "macrodef.fh"

! rhs for scalar and GR variables
! here we consider vacuum spacetime only
  function compute_rhs_bssn_escalar(ex, T,X, Y, Z,                             &
               chi    ,   trK    ,                                             &
               dxx    ,   gxy    ,   gxz    ,   dyy    ,   gyz    ,   dzz,     &
               Axx    ,   Axy    ,   Axz    ,   Ayy    ,   Ayz    ,   Azz,     &
               Gamx   ,  Gamy    ,  Gamz    ,                                  &
               Lap    ,  betax   ,  betay   ,  betaz   ,                       &
               dtSfx  ,  dtSfy   ,  dtSfz   ,                                  &
               Sphi   ,   Spi    ,                                             &
               chi_rhs,   trK_rhs,                                             &
               gxx_rhs,   gxy_rhs,   gxz_rhs,   gyy_rhs,   gyz_rhs,   gzz_rhs, &
               Axx_rhs,   Axy_rhs,   Axz_rhs,   Ayy_rhs,   Ayz_rhs,   Azz_rhs, &
               Gamx_rhs,  Gamy_rhs,  Gamz_rhs,                                 &
               Lap_rhs,  betax_rhs,  betay_rhs,  betaz_rhs,                    &
               dtSfx_rhs,  dtSfy_rhs,  dtSfz_rhs,                              &
               Sphi_rhs ,   Spi_rhs ,                                          &
               rho,Sx,Sy,Sz,Sxx,Sxy,Sxz,Syy,Syz,Szz,                           &
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,                      &
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,                      &
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,                      &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,                                        &
               ham_Res, movx_Res, movy_Res, movz_Res,                          &
                        Gmx_Res, Gmy_Res, Gmz_Res,                             &
               Symmetry,Lev,eps,co)  result(gont)
! calculate constraint violation when co=0       
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
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Sphi,Spi
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: chi_rhs,trK_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gxx_rhs,gxy_rhs,gxz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gyy_rhs,gyz_rhs,gzz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Axx_rhs,Axy_rhs,Axz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Ayy_rhs,Ayz_rhs,Azz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamx_rhs,Gamy_rhs,Gamz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Lap_rhs, betax_rhs, betay_rhs, betaz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: dtSfx_rhs,dtSfy_rhs,dtSfz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Sphi_rhs,Spi_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: rho,Sx,Sy,Sz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Sxx,Sxy,Sxz,Syy,Syz,Szz
! when out, physical second kind of connection  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamzyy, Gamzyz, Gamzzz
! when out, physical Ricci tensor  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Rxx,Rxy,Rxz,Ryy,Ryz,Rzz
  real*8,intent(in) :: eps
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: ham_Res, movx_Res, movy_Res, movz_Res
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gmx_Res, Gmy_Res, Gmz_Res
!  gont = 0: success; gont = 1: something wrong
  integer::gont

!~~~~~~> Other variables:

  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gyy,gzz
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
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamxa,Gamya,Gamza,alpn1,chin1
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupyy,gupyz,gupzz

  real*8,dimension(3) ::SSS,AAS,ASA,SAA,ASS,SAS,SSA
  real*8            :: dX, dY, dZ, PI
  real*8, parameter :: ZEO=0.d0, ONE = 1.D0, TWO = 2.D0, FOUR = 4.D0
  real*8, parameter :: EIGHT = 8.D0, HALF = 0.5D0, THR = 3.d0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  double precision,parameter::FF = 0.75d0,eta=2.d0
  real*8, parameter :: F1o3 = 1.D0/3.D0, F2o3 = 2.D0/3.D0,F3o2=1.5d0, F1o6 = 1.D0/6.D0
  real*8, parameter :: F16=1.6d1,F8=8.d0

  integer :: i,j,k

!!! sanity check
  dX = sum(chi)+sum(trK)+sum(dxx)+sum(gxy)+sum(gxz)+sum(dyy)+sum(gyz)+sum(dzz) &
      +sum(Axx)+sum(Axy)+sum(Axz)+sum(Ayy)+sum(Ayz)+sum(Azz)                   &
      +sum(Gamx)+sum(Gamy)+sum(Gamz)                                           &
      +sum(Lap)+sum(betax)+sum(betay)+sum(betaz)                               &
      +sum(Sphi)+sum(Spi)
  if(dX.ne.dX) then
     if(sum(chi).ne.sum(chi))write(*,*)"bssn.f90: find NaN in chi"
     if(sum(trK).ne.sum(trK))write(*,*)"bssn.f90: find NaN in trk"
     if(sum(dxx).ne.sum(dxx))write(*,*)"bssn.f90: find NaN in dxx"
     if(sum(gxy).ne.sum(gxy))write(*,*)"bssn.f90: find NaN in gxy"
     if(sum(gxz).ne.sum(gxz))write(*,*)"bssn.f90: find NaN in gxz"
     if(sum(dyy).ne.sum(dyy))write(*,*)"bssn.f90: find NaN in dyy"
     if(sum(gyz).ne.sum(gyz))write(*,*)"bssn.f90: find NaN in gyz"
     if(sum(dzz).ne.sum(dzz))write(*,*)"bssn.f90: find NaN in dzz"
     if(sum(Axx).ne.sum(Axx))write(*,*)"bssn.f90: find NaN in Axx"
     if(sum(Axy).ne.sum(Axy))write(*,*)"bssn.f90: find NaN in Axy"
     if(sum(Axz).ne.sum(Axz))write(*,*)"bssn.f90: find NaN in Axz"
     if(sum(Ayy).ne.sum(Ayy))write(*,*)"bssn.f90: find NaN in Ayy"
     if(sum(Ayz).ne.sum(Ayz))write(*,*)"bssn.f90: find NaN in Ayz"
     if(sum(Azz).ne.sum(Azz))write(*,*)"bssn.f90: find NaN in Azz"
     if(sum(Gamx).ne.sum(Gamx))write(*,*)"bssn.f90: find NaN in Gamx"
     if(sum(Gamy).ne.sum(Gamy))write(*,*)"bssn.f90: find NaN in Gamy"
     if(sum(Gamz).ne.sum(Gamz))write(*,*)"bssn.f90: find NaN in Gamz"
     if(sum(Lap).ne.sum(Lap))write(*,*)"bssn.f90: find NaN in Lap"
     if(sum(betax).ne.sum(betax))write(*,*)"bssn.f90: find NaN in betax"
     if(sum(betay).ne.sum(betay))write(*,*)"bssn.f90: find NaN in betay"
     if(sum(betaz).ne.sum(betaz))write(*,*)"bssn.f90: find NaN in betaz"
     if(sum(Sphi).ne.sum(Sphi))write(*,*)"bssn.f90: find NaN in Sphi"
     if(sum(Spi).ne.sum(Spi))write(*,*)"bssn.f90: find NaN in Spi"
     gont = 1
     return
  endif

  PI = dacos(-ONE)

  dX = X(2) - X(1)
  dY = Y(2) - Y(1)
  dZ = Z(2) - Z(1)

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
  alpn1(i,j,k) = Lap(i,j,k) + ONE
  chin1(i,j,k) = chi(i,j,k) + ONE
  gxx(i,j,k) = dxx(i,j,k) + ONE
  gyy(i,j,k) = dyy(i,j,k) + ONE
  gzz(i,j,k) = dzz(i,j,k) + ONE
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

  call fderivs(ex,betax,betaxx,betaxy,betaxz,X,Y,Z,ANTI, SYM, SYM,Symmetry,Lev)
  call fderivs(ex,betay,betayx,betayy,betayz,X,Y,Z, SYM,ANTI, SYM,Symmetry,Lev)
  call fderivs(ex,betaz,betazx,betazy,betazz,X,Y,Z, SYM, SYM,ANTI,Symmetry,Lev)
  
  div_beta = betaxx + betayy + betazz
 
  call fderivs(ex,chi,chix,chiy,chiz,X,Y,Z,SYM,SYM,SYM,symmetry,Lev)

  chi_rhs = F2o3 *chin1*( alpn1 * trK - div_beta ) !rhs for chi

  call fderivs(ex,Lap,Lapx,Lapy,Lapz,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)
! invert tilted metric
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
  gupzz(i,j,k) =  gxx(i,j,k) * gyy(i,j,k) * gzz(i,j,k) + gxy(i,j,k) * gyz(i,j,k) * gxz(i,j,k) + gxz(i,j,k) * gxy(i,j,k) * gyz(i,j,k) - &
           gxz(i,j,k) * gyy(i,j,k) * gxz(i,j,k) - gxy(i,j,k) * gxy(i,j,k) * gzz(i,j,k) - gxx(i,j,k) * gyz(i,j,k) * gyz(i,j,k)
  gupxx(i,j,k) =   ( gyy(i,j,k) * gzz(i,j,k) - gyz(i,j,k) * gyz(i,j,k) ) / gupzz(i,j,k)
  gupxy(i,j,k) = - ( gxy(i,j,k) * gzz(i,j,k) - gyz(i,j,k) * gxz(i,j,k) ) / gupzz(i,j,k)
  gupxz(i,j,k) =   ( gxy(i,j,k) * gyz(i,j,k) - gyy(i,j,k) * gxz(i,j,k) ) / gupzz(i,j,k)
  gupyy(i,j,k) =   ( gxx(i,j,k) * gzz(i,j,k) - gxz(i,j,k) * gxz(i,j,k) ) / gupzz(i,j,k)
  gupyz(i,j,k) = - ( gxx(i,j,k) * gyz(i,j,k) - gxy(i,j,k) * gxz(i,j,k) ) / gupzz(i,j,k)
  gupzz(i,j,k) =   ( gxx(i,j,k) * gyy(i,j,k) - gxy(i,j,k) * gxy(i,j,k) ) / gupzz(i,j,k)
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

#if 1 
  Sphi_rhs = alpn1 * Spi !rhs for Scalar phi
!rhs for Spi
   call fderivs(ex,Sphi,Kx,Ky,Kz,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)
   call fdderivs(ex,Sphi,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
   Spi_rhs(i,j,k) =   gupxx(i,j,k) * fxx(i,j,k) + gupyy(i,j,k) * fyy(i,j,k) + gupzz(i,j,k) * fzz(i,j,k) +           &
             ( gupxy(i,j,k) * fxy(i,j,k) + gupxz(i,j,k) * fxz(i,j,k) + gupyz(i,j,k) * fyz(i,j,k) ) * TWO -   &
             ((Gamx(i,j,k)+(gupxx(i,j,k)*chix(i,j,k)+gupxy(i,j,k)*chiy(i,j,k)+gupxz(i,j,k)*chiz(i,j,k))/TWO/chin1(i,j,k))*Kx(i,j,k)  &
            + (Gamy(i,j,k)+(gupxy(i,j,k)*chix(i,j,k)+gupyy(i,j,k)*chiy(i,j,k)+gupyz(i,j,k)*chiz(i,j,k))/TWO/chin1(i,j,k))*Ky(i,j,k)  &
            + (Gamz(i,j,k)+(gupxz(i,j,k)*chix(i,j,k)+gupyz(i,j,k)*chiy(i,j,k)+gupzz(i,j,k)*chiz(i,j,k))/TWO/chin1(i,j,k))*Kz(i,j,k))
   Spi_rhs(i,j,k) = Spi_rhs(i,j,k)*alpn1(i,j,k) +                                       &
            (gupxx(i,j,k)*Lapx(i,j,k)*Kx(i,j,k) + gupxy(i,j,k)*Lapx(i,j,k)*Ky(i,j,k) + gupxz(i,j,k)*Lapx(i,j,k)*Kz(i,j,k)&
            +gupxy(i,j,k)*Lapy(i,j,k)*Kx(i,j,k) + gupyy(i,j,k)*Lapy(i,j,k)*Ky(i,j,k) + gupyz(i,j,k)*Lapy(i,j,k)*Kz(i,j,k)&
            +gupxz(i,j,k)*Lapz(i,j,k)*Kx(i,j,k) + gupyz(i,j,k)*Lapz(i,j,k)*Ky(i,j,k) + gupzz(i,j,k)*Lapz(i,j,k)*Kz(i,j,k))
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO
 
   call frpotential(ex,Sphi,f,S)
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
   Spi_rhs(i,j,k) = Spi_rhs(i,j,k)*chin1(i,j,k) + alpn1(i,j,k)*(trK(i,j,k)*Spi(i,j,k) - S(i,j,k))
! matter content of scalar field T_ab 
  rho(i,j,k) = chin1(i,j,k)*((gupxx(i,j,k) * Kx(i,j,k) * Kx(i,j,k) + gupyy(i,j,k) * Ky(i,j,k) * Ky(i,j,k) + gupzz(i,j,k) * Kz(i,j,k) * Kz(i,j,k))/TWO +  &
                gupxy(i,j,k) * Kx(i,j,k) * Ky(i,j,k) + gupxz(i,j,k) * Kx(i,j,k) * Kz(i,j,k) + gupyz(i,j,k) * Ky(i,j,k) * Kz(i,j,k) ) &
        + Spi(i,j,k)*Spi(i,j,k)/TWO+f(i,j,k)  
  Sx(i,j,k) = -Spi(i,j,k)*Kx(i,j,k)
  Sy(i,j,k) = -Spi(i,j,k)*Ky(i,j,k)
  Sz(i,j,k) = -Spi(i,j,k)*Kz(i,j,k)
  f(i,j,k) = (rho(i,j,k) - Spi(i,j,k)*Spi(i,j,k))/chin1(i,j,k)
  Sxx(i,j,k) = Kx(i,j,k)*Kx(i,j,k)-f(i,j,k)*gxx(i,j,k)
  Sxy(i,j,k) = Kx(i,j,k)*Ky(i,j,k)-f(i,j,k)*gxy(i,j,k)
  Sxz(i,j,k) = Kx(i,j,k)*Kz(i,j,k)-f(i,j,k)*gxz(i,j,k)
  Syy(i,j,k) = Ky(i,j,k)*Ky(i,j,k)-f(i,j,k)*gyy(i,j,k)
  Syz(i,j,k) = Ky(i,j,k)*Kz(i,j,k)-f(i,j,k)*gyz(i,j,k)
  Szz(i,j,k) = Kz(i,j,k)*Kz(i,j,k)-f(i,j,k)*gzz(i,j,k)
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO
#else
  rho = ZEO
  Sx = ZEO
  Sy = ZEO
  Sz = ZEO
  Sxx = ZEO
  Sxy = ZEO
  Sxz = ZEO
  Syy = ZEO
  Syz = ZEO
  Szz = ZEO
#endif  

  call fderivs(ex,dxx,gxxx,gxxy,gxxz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)
  call fderivs(ex,gxy,gxyx,gxyy,gxyz,X,Y,Z,ANTI,ANTI,SYM ,Symmetry,Lev)
  call fderivs(ex,gxz,gxzx,gxzy,gxzz,X,Y,Z,ANTI,SYM ,ANTI,Symmetry,Lev)
  call fderivs(ex,dyy,gyyx,gyyy,gyyz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)
  call fderivs(ex,gyz,gyzx,gyzy,gyzz,X,Y,Z,SYM ,ANTI,ANTI,Symmetry,Lev)
  call fderivs(ex,dzz,gzzx,gzzy,gzzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)

  if(co == 0)then
! Gam^i_Res = Gam^i + gup^ij_,j
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
  Gmx_Res(i,j,k) = Gamx(i,j,k) - (gupxx(i,j,k)*(gupxx(i,j,k)*gxxx(i,j,k)+gupxy(i,j,k)*gxyx(i,j,k)+gupxz(i,j,k)*gxzx(i,j,k))&
                   +gupxy(i,j,k)*(gupxx(i,j,k)*gxyx(i,j,k)+gupxy(i,j,k)*gyyx(i,j,k)+gupxz(i,j,k)*gyzx(i,j,k))&
                   +gupxz(i,j,k)*(gupxx(i,j,k)*gxzx(i,j,k)+gupxy(i,j,k)*gyzx(i,j,k)+gupxz(i,j,k)*gzzx(i,j,k))&
                   +gupxx(i,j,k)*(gupxy(i,j,k)*gxxy(i,j,k)+gupyy(i,j,k)*gxyy(i,j,k)+gupyz(i,j,k)*gxzy(i,j,k))&
                   +gupxy(i,j,k)*(gupxy(i,j,k)*gxyy(i,j,k)+gupyy(i,j,k)*gyyy(i,j,k)+gupyz(i,j,k)*gyzy(i,j,k))&
                   +gupxz(i,j,k)*(gupxy(i,j,k)*gxzy(i,j,k)+gupyy(i,j,k)*gyzy(i,j,k)+gupyz(i,j,k)*gzzy(i,j,k))&
                   +gupxx(i,j,k)*(gupxz(i,j,k)*gxxz(i,j,k)+gupyz(i,j,k)*gxyz(i,j,k)+gupzz(i,j,k)*gxzz(i,j,k))&
                   +gupxy(i,j,k)*(gupxz(i,j,k)*gxyz(i,j,k)+gupyz(i,j,k)*gyyz(i,j,k)+gupzz(i,j,k)*gyzz(i,j,k))&
                   +gupxz(i,j,k)*(gupxz(i,j,k)*gxzz(i,j,k)+gupyz(i,j,k)*gyzz(i,j,k)+gupzz(i,j,k)*gzzz(i,j,k)))
  Gmy_Res(i,j,k) = Gamy(i,j,k) - (gupxx(i,j,k)*(gupxy(i,j,k)*gxxx(i,j,k)+gupyy(i,j,k)*gxyx(i,j,k)+gupyz(i,j,k)*gxzx(i,j,k))&
                   +gupxy(i,j,k)*(gupxy(i,j,k)*gxyx(i,j,k)+gupyy(i,j,k)*gyyx(i,j,k)+gupyz(i,j,k)*gyzx(i,j,k))&
                   +gupxz(i,j,k)*(gupxy(i,j,k)*gxzx(i,j,k)+gupyy(i,j,k)*gyzx(i,j,k)+gupyz(i,j,k)*gzzx(i,j,k))&
                   +gupxy(i,j,k)*(gupxy(i,j,k)*gxxy(i,j,k)+gupyy(i,j,k)*gxyy(i,j,k)+gupyz(i,j,k)*gxzy(i,j,k))&
                   +gupyy(i,j,k)*(gupxy(i,j,k)*gxyy(i,j,k)+gupyy(i,j,k)*gyyy(i,j,k)+gupyz(i,j,k)*gyzy(i,j,k))&
                   +gupyz(i,j,k)*(gupxy(i,j,k)*gxzy(i,j,k)+gupyy(i,j,k)*gyzy(i,j,k)+gupyz(i,j,k)*gzzy(i,j,k))&
                   +gupxy(i,j,k)*(gupxz(i,j,k)*gxxz(i,j,k)+gupyz(i,j,k)*gxyz(i,j,k)+gupzz(i,j,k)*gxzz(i,j,k))&
                   +gupyy(i,j,k)*(gupxz(i,j,k)*gxyz(i,j,k)+gupyz(i,j,k)*gyyz(i,j,k)+gupzz(i,j,k)*gyzz(i,j,k))&
                   +gupyz(i,j,k)*(gupxz(i,j,k)*gxzz(i,j,k)+gupyz(i,j,k)*gyzz(i,j,k)+gupzz(i,j,k)*gzzz(i,j,k)))
  Gmz_Res(i,j,k) = Gamz(i,j,k) - (gupxx(i,j,k)*(gupxz(i,j,k)*gxxx(i,j,k)+gupyz(i,j,k)*gxyx(i,j,k)+gupzz(i,j,k)*gxzx(i,j,k))&
                   +gupxy(i,j,k)*(gupxz(i,j,k)*gxyx(i,j,k)+gupyz(i,j,k)*gyyx(i,j,k)+gupzz(i,j,k)*gyzx(i,j,k))&
                   +gupxz(i,j,k)*(gupxz(i,j,k)*gxzx(i,j,k)+gupyz(i,j,k)*gyzx(i,j,k)+gupzz(i,j,k)*gzzx(i,j,k))&
                   +gupxy(i,j,k)*(gupxz(i,j,k)*gxxy(i,j,k)+gupyz(i,j,k)*gxyy(i,j,k)+gupzz(i,j,k)*gxzy(i,j,k))&
                   +gupyy(i,j,k)*(gupxz(i,j,k)*gxyy(i,j,k)+gupyz(i,j,k)*gyyy(i,j,k)+gupzz(i,j,k)*gyzy(i,j,k))&
                   +gupyz(i,j,k)*(gupxz(i,j,k)*gxzy(i,j,k)+gupyz(i,j,k)*gyzy(i,j,k)+gupzz(i,j,k)*gzzy(i,j,k))&
                   +gupxz(i,j,k)*(gupxz(i,j,k)*gxxz(i,j,k)+gupyz(i,j,k)*gxyz(i,j,k)+gupzz(i,j,k)*gxzz(i,j,k))&
                   +gupyz(i,j,k)*(gupxz(i,j,k)*gxyz(i,j,k)+gupyz(i,j,k)*gyyz(i,j,k)+gupzz(i,j,k)*gyzz(i,j,k))&
                   +gupzz(i,j,k)*(gupxz(i,j,k)*gxzz(i,j,k)+gupyz(i,j,k)*gyzz(i,j,k)+gupzz(i,j,k)*gzzz(i,j,k)))
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO
  endif

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
  gxx_rhs(i,j,k) = - TWO * alpn1(i,j,k) * Axx(i,j,k)    -  F2o3 * gxx(i,j,k) * div_beta(i,j,k)          + &
              TWO *(  gxx(i,j,k) * betaxx(i,j,k) +   gxy(i,j,k) * betayx(i,j,k) +   gxz(i,j,k) * betazx(i,j,k))

  gyy_rhs(i,j,k) = - TWO * alpn1(i,j,k) * Ayy(i,j,k)    -  F2o3 * gyy(i,j,k) * div_beta(i,j,k)          + &
              TWO *(  gxy(i,j,k) * betaxy(i,j,k) +   gyy(i,j,k) * betayy(i,j,k) +   gyz(i,j,k) * betazy(i,j,k))

  gzz_rhs(i,j,k) = - TWO * alpn1(i,j,k) * Azz(i,j,k)    -  F2o3 * gzz(i,j,k) * div_beta(i,j,k)          + &
              TWO *(  gxz(i,j,k) * betaxz(i,j,k) +   gyz(i,j,k) * betayz(i,j,k) +   gzz(i,j,k) * betazz(i,j,k))

  gxy_rhs(i,j,k) = - TWO * alpn1(i,j,k) * Axy(i,j,k)    +  F1o3 * gxy(i,j,k)    * div_beta(i,j,k)       + &
                      gxx(i,j,k) * betaxy(i,j,k)                  +   gxz(i,j,k) * betazy(i,j,k) + &
                                       gyy(i,j,k) * betayx(i,j,k) +   gyz(i,j,k) * betazx(i,j,k)   &
                                                    -   gxy(i,j,k) * betazz(i,j,k)

  gyz_rhs(i,j,k) = - TWO * alpn1(i,j,k) * Ayz(i,j,k)    +  F1o3 * gyz(i,j,k)    * div_beta(i,j,k)       + &
                      gxy(i,j,k) * betaxz(i,j,k) +   gyy(i,j,k) * betayz(i,j,k)                  + &
                      gxz(i,j,k) * betaxy(i,j,k)                  +   gzz(i,j,k) * betazy(i,j,k)   &
                                                    -   gyz(i,j,k) * betaxx(i,j,k)
 
  gxz_rhs(i,j,k) = - TWO * alpn1(i,j,k) * Axz(i,j,k)    +  F1o3 * gxz(i,j,k)    * div_beta(i,j,k)       + &
                      gxx(i,j,k) * betaxz(i,j,k) +   gxy(i,j,k) * betayz(i,j,k)                  + &
                                       gyz(i,j,k) * betayx(i,j,k) +   gzz(i,j,k) * betazx(i,j,k)   &
                                                    -   gxz(i,j,k) * betayy(i,j,k)     !rhs for gij

! second kind of connection
  Gamxxx(i,j,k) =HALF*( gupxx(i,j,k)*gxxx(i,j,k) + gupxy(i,j,k)*(TWO*gxyx(i,j,k) - gxxy(i,j,k) ) + gupxz(i,j,k)*(TWO*gxzx(i,j,k) - gxxz(i,j,k) ))
  Gamyxx(i,j,k) =HALF*( gupxy(i,j,k)*gxxx(i,j,k) + gupyy(i,j,k)*(TWO*gxyx(i,j,k) - gxxy(i,j,k) ) + gupyz(i,j,k)*(TWO*gxzx(i,j,k) - gxxz(i,j,k) ))
  Gamzxx(i,j,k) =HALF*( gupxz(i,j,k)*gxxx(i,j,k) + gupyz(i,j,k)*(TWO*gxyx(i,j,k) - gxxy(i,j,k) ) + gupzz(i,j,k)*(TWO*gxzx(i,j,k) - gxxz(i,j,k) ))
 
  Gamxyy(i,j,k) =HALF*( gupxx(i,j,k)*(TWO*gxyy(i,j,k) - gyyx(i,j,k) ) + gupxy(i,j,k)*gyyy(i,j,k) + gupxz(i,j,k)*(TWO*gyzy(i,j,k) - gyyz(i,j,k) ))
  Gamyyy(i,j,k) =HALF*( gupxy(i,j,k)*(TWO*gxyy(i,j,k) - gyyx(i,j,k) ) + gupyy(i,j,k)*gyyy(i,j,k) + gupyz(i,j,k)*(TWO*gyzy(i,j,k) - gyyz(i,j,k) ))
  Gamzyy(i,j,k) =HALF*( gupxz(i,j,k)*(TWO*gxyy(i,j,k) - gyyx(i,j,k) ) + gupyz(i,j,k)*gyyy(i,j,k) + gupzz(i,j,k)*(TWO*gyzy(i,j,k) - gyyz(i,j,k) ))

  Gamxzz(i,j,k) =HALF*( gupxx(i,j,k)*(TWO*gxzz(i,j,k) - gzzx(i,j,k) ) + gupxy(i,j,k)*(TWO*gyzz(i,j,k) - gzzy(i,j,k) ) + gupxz(i,j,k)*gzzz(i,j,k))
  Gamyzz(i,j,k) =HALF*( gupxy(i,j,k)*(TWO*gxzz(i,j,k) - gzzx(i,j,k) ) + gupyy(i,j,k)*(TWO*gyzz(i,j,k) - gzzy(i,j,k) ) + gupyz(i,j,k)*gzzz(i,j,k))
  Gamzzz(i,j,k) =HALF*( gupxz(i,j,k)*(TWO*gxzz(i,j,k) - gzzx(i,j,k) ) + gupyz(i,j,k)*(TWO*gyzz(i,j,k) - gzzy(i,j,k) ) + gupzz(i,j,k)*gzzz(i,j,k))

  Gamxxy(i,j,k) =HALF*( gupxx(i,j,k)*gxxy(i,j,k) + gupxy(i,j,k)*gyyx(i,j,k) + gupxz(i,j,k)*( gxzy(i,j,k) + gyzx(i,j,k) - gxyz(i,j,k) ) )
  Gamyxy(i,j,k) =HALF*( gupxy(i,j,k)*gxxy(i,j,k) + gupyy(i,j,k)*gyyx(i,j,k) + gupyz(i,j,k)*( gxzy(i,j,k) + gyzx(i,j,k) - gxyz(i,j,k) ) )
  Gamzxy(i,j,k) =HALF*( gupxz(i,j,k)*gxxy(i,j,k) + gupyz(i,j,k)*gyyx(i,j,k) + gupzz(i,j,k)*( gxzy(i,j,k) + gyzx(i,j,k) - gxyz(i,j,k) ) )

  Gamxxz(i,j,k) =HALF*( gupxx(i,j,k)*gxxz(i,j,k) + gupxy(i,j,k)*( gxyz(i,j,k) + gyzx(i,j,k) - gxzy(i,j,k) ) + gupxz(i,j,k)*gzzx(i,j,k) )
  Gamyxz(i,j,k) =HALF*( gupxy(i,j,k)*gxxz(i,j,k) + gupyy(i,j,k)*( gxyz(i,j,k) + gyzx(i,j,k) - gxzy(i,j,k) ) + gupyz(i,j,k)*gzzx(i,j,k) )
  Gamzxz(i,j,k) =HALF*( gupxz(i,j,k)*gxxz(i,j,k) + gupyz(i,j,k)*( gxyz(i,j,k) + gyzx(i,j,k) - gxzy(i,j,k) ) + gupzz(i,j,k)*gzzx(i,j,k) )

  Gamxyz(i,j,k) =HALF*( gupxx(i,j,k)*( gxyz(i,j,k) + gxzy(i,j,k) - gyzx(i,j,k) ) + gupxy(i,j,k)*gyyz(i,j,k) + gupxz(i,j,k)*gzzy(i,j,k) )
  Gamyyz(i,j,k) =HALF*( gupxy(i,j,k)*( gxyz(i,j,k) + gxzy(i,j,k) - gyzx(i,j,k) ) + gupyy(i,j,k)*gyyz(i,j,k) + gupyz(i,j,k)*gzzy(i,j,k) )
  Gamzyz(i,j,k) =HALF*( gupxz(i,j,k)*( gxyz(i,j,k) + gxzy(i,j,k) - gyzx(i,j,k) ) + gupyz(i,j,k)*gyyz(i,j,k) + gupzz(i,j,k)*gzzy(i,j,k) )
! Raise indices of \tilde A_{ij} and store in R_ij

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
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

! Right hand side for Gam^i without shift terms...
  call fderivs(ex,trK,Kx,Ky,Kz,X,Y,Z,SYM,SYM,SYM,symmetry,Lev)

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
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
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

  call fdderivs(ex,betax,gxxx,gxyx,gxzx,gyyx,gyzx,gzzx,&
                X,Y,Z,ANTI,SYM, SYM ,Symmetry,Lev)
  call fdderivs(ex,betay,gxxy,gxyy,gxzy,gyyy,gyzy,gzzy,&
                X,Y,Z,SYM ,ANTI,SYM ,Symmetry,Lev)
  call fdderivs(ex,betaz,gxxz,gxyz,gxzz,gyyz,gyzz,gzzz,&
                X,Y,Z,SYM ,SYM, ANTI,Symmetry,Lev)

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
  fxx(i,j,k) = gxxx(i,j,k) + gxyy(i,j,k) + gxzz(i,j,k)
  fxy(i,j,k) = gxyx(i,j,k) + gyyy(i,j,k) + gyzz(i,j,k)
  fxz(i,j,k) = gxzx(i,j,k) + gyzy(i,j,k) + gzzz(i,j,k)

  Gamxa(i,j,k) =       gupxx(i,j,k) * Gamxxx(i,j,k) + gupyy(i,j,k) * Gamxyy(i,j,k) + gupzz(i,j,k) * Gamxzz(i,j,k) + &
          TWO*( gupxy(i,j,k) * Gamxxy(i,j,k) + gupxz(i,j,k) * Gamxxz(i,j,k) + gupyz(i,j,k) * Gamxyz(i,j,k) )
  Gamya(i,j,k) =       gupxx(i,j,k) * Gamyxx(i,j,k) + gupyy(i,j,k) * Gamyyy(i,j,k) + gupzz(i,j,k) * Gamyzz(i,j,k) + &
          TWO*( gupxy(i,j,k) * Gamyxy(i,j,k) + gupxz(i,j,k) * Gamyxz(i,j,k) + gupyz(i,j,k) * Gamyyz(i,j,k) )
  Gamza(i,j,k) =       gupxx(i,j,k) * Gamzxx(i,j,k) + gupyy(i,j,k) * Gamzyy(i,j,k) + gupzz(i,j,k) * Gamzzz(i,j,k) + &
          TWO*( gupxy(i,j,k) * Gamzxy(i,j,k) + gupxz(i,j,k) * Gamzxz(i,j,k) + gupyz(i,j,k) * Gamzyz(i,j,k) )
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

  call fderivs(ex,Gamx,Gamxx,Gamxy,Gamxz,X,Y,Z,ANTI,SYM ,SYM ,Symmetry,Lev)
  call fderivs(ex,Gamy,Gamyx,Gamyy,Gamyz,X,Y,Z,SYM ,ANTI,SYM ,Symmetry,Lev)
  call fderivs(ex,Gamz,Gamzx,Gamzy,Gamzz,X,Y,Z,SYM ,SYM ,ANTI,Symmetry,Lev)

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
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

!first kind of connection stored in gij,k
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
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

!compute Ricci tensor for tilted metric
   call fdderivs(ex,dxx,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM ,SYM ,SYM ,symmetry,Lev)
   Rxx =   gupxx * fxx + gupyy * fyy + gupzz * fzz + &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) * TWO

   call fdderivs(ex,dyy,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM ,SYM ,SYM ,symmetry,Lev)
   Ryy =   gupxx * fxx + gupyy * fyy + gupzz * fzz + &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) * TWO

   call fdderivs(ex,dzz,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM ,SYM ,SYM ,symmetry,Lev)
   Rzz =   gupxx * fxx + gupyy * fyy + gupzz * fzz + &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) * TWO

   call fdderivs(ex,gxy,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,ANTI, ANTI,SYM ,symmetry,Lev)
   Rxy =   gupxx * fxx + gupyy * fyy + gupzz * fzz + &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) * TWO

   call fdderivs(ex,gxz,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,ANTI ,SYM ,ANTI,symmetry,Lev)
   Rxz =   gupxx * fxx + gupyy * fyy + gupzz * fzz + &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) * TWO

   call fdderivs(ex,gyz,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM ,ANTI ,ANTI,symmetry,Lev)
   Ryz =   gupxx * fxx + gupyy * fyy + gupzz * fzz + &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) * TWO

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
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
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO
!covariant second derivative of chi respect to tilted metric
  call fdderivs(ex,chi,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
  fxx(i,j,k) = fxx(i,j,k) - Gamxxx(i,j,k) * chix(i,j,k) - Gamyxx(i,j,k) * chiy(i,j,k) - Gamzxx(i,j,k) * chiz(i,j,k)
  fxy(i,j,k) = fxy(i,j,k) - Gamxxy(i,j,k) * chix(i,j,k) - Gamyxy(i,j,k) * chiy(i,j,k) - Gamzxy(i,j,k) * chiz(i,j,k)
  fxz(i,j,k) = fxz(i,j,k) - Gamxxz(i,j,k) * chix(i,j,k) - Gamyxz(i,j,k) * chiy(i,j,k) - Gamzxz(i,j,k) * chiz(i,j,k)
  fyy(i,j,k) = fyy(i,j,k) - Gamxyy(i,j,k) * chix(i,j,k) - Gamyyy(i,j,k) * chiy(i,j,k) - Gamzyy(i,j,k) * chiz(i,j,k)
  fyz(i,j,k) = fyz(i,j,k) - Gamxyz(i,j,k) * chix(i,j,k) - Gamyyz(i,j,k) * chiy(i,j,k) - Gamzyz(i,j,k) * chiz(i,j,k)
  fzz(i,j,k) = fzz(i,j,k) - Gamxzz(i,j,k) * chix(i,j,k) - Gamyzz(i,j,k) * chiy(i,j,k) - Gamzzz(i,j,k) * chiz(i,j,k)
! Store D^l D_l chi(i,j,k) - 3/(2*chi(i,j,k)) D^l chi(i,j,k) D_l chi(i,j,k) in f(i,j,k)

  f(i,j,k) =        gupxx(i,j,k) * ( fxx(i,j,k) - F3o2/chin1(i,j,k) * chix(i,j,k) * chix(i,j,k) ) + &
             gupyy(i,j,k) * ( fyy(i,j,k) - F3o2/chin1(i,j,k) * chiy(i,j,k) * chiy(i,j,k) ) + &
             gupzz(i,j,k) * ( fzz(i,j,k) - F3o2/chin1(i,j,k) * chiz(i,j,k) * chiz(i,j,k) ) + &
       TWO * gupxy(i,j,k) * ( fxy(i,j,k) - F3o2/chin1(i,j,k) * chix(i,j,k) * chiy(i,j,k) ) + &
       TWO * gupxz(i,j,k) * ( fxz(i,j,k) - F3o2/chin1(i,j,k) * chix(i,j,k) * chiz(i,j,k) ) + &
       TWO * gupyz(i,j,k) * ( fyz(i,j,k) - F3o2/chin1(i,j,k) * chiy(i,j,k) * chiz(i,j,k) ) 
! Add chi(i,j,k) part to Ricci tensor:

  Rxx(i,j,k) = Rxx(i,j,k) + (fxx(i,j,k) - chix(i,j,k)*chix(i,j,k)/chin1(i,j,k)/TWO + gxx(i,j,k) * f(i,j,k))/chin1(i,j,k)/TWO
  Ryy(i,j,k) = Ryy(i,j,k) + (fyy(i,j,k) - chiy(i,j,k)*chiy(i,j,k)/chin1(i,j,k)/TWO + gyy(i,j,k) * f(i,j,k))/chin1(i,j,k)/TWO
  Rzz(i,j,k) = Rzz(i,j,k) + (fzz(i,j,k) - chiz(i,j,k)*chiz(i,j,k)/chin1(i,j,k)/TWO + gzz(i,j,k) * f(i,j,k))/chin1(i,j,k)/TWO
  Rxy(i,j,k) = Rxy(i,j,k) + (fxy(i,j,k) - chix(i,j,k)*chiy(i,j,k)/chin1(i,j,k)/TWO + gxy(i,j,k) * f(i,j,k))/chin1(i,j,k)/TWO
  Rxz(i,j,k) = Rxz(i,j,k) + (fxz(i,j,k) - chix(i,j,k)*chiz(i,j,k)/chin1(i,j,k)/TWO + gxz(i,j,k) * f(i,j,k))/chin1(i,j,k)/TWO
  Ryz(i,j,k) = Ryz(i,j,k) + (fyz(i,j,k) - chiy(i,j,k)*chiz(i,j,k)/chin1(i,j,k)/TWO + gyz(i,j,k) * f(i,j,k))/chin1(i,j,k)/TWO

! now prepare to get physical second kind of connection
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
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

! covariant second derivatives of the lapse respect to physical metric
  call fdderivs(ex,Lap,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z, &
                SYM,SYM,SYM,symmetry,Lev)

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
  fxx(i,j,k) = fxx(i,j,k) - Gamxxx(i,j,k)*Lapx(i,j,k) - Gamyxx(i,j,k)*Lapy(i,j,k) - Gamzxx(i,j,k)*Lapz(i,j,k)
  fyy(i,j,k) = fyy(i,j,k) - Gamxyy(i,j,k)*Lapx(i,j,k) - Gamyyy(i,j,k)*Lapy(i,j,k) - Gamzyy(i,j,k)*Lapz(i,j,k)
  fzz(i,j,k) = fzz(i,j,k) - Gamxzz(i,j,k)*Lapx(i,j,k) - Gamyzz(i,j,k)*Lapy(i,j,k) - Gamzzz(i,j,k)*Lapz(i,j,k)
  fxy(i,j,k) = fxy(i,j,k) - Gamxxy(i,j,k)*Lapx(i,j,k) - Gamyxy(i,j,k)*Lapy(i,j,k) - Gamzxy(i,j,k)*Lapz(i,j,k)
  fxz(i,j,k) = fxz(i,j,k) - Gamxxz(i,j,k)*Lapx(i,j,k) - Gamyxz(i,j,k)*Lapy(i,j,k) - Gamzxz(i,j,k)*Lapz(i,j,k)
  fyz(i,j,k) = fyz(i,j,k) - Gamxyz(i,j,k)*Lapx(i,j,k) - Gamyyz(i,j,k)*Lapy(i,j,k) - Gamzyz(i,j,k)*Lapz(i,j,k)

! mater content only appears in the right hand side of trK(i,j,k) and A_ij which are
! comming now
! store D^i D_i Lap(i,j,k) in trK_rhs(i,j,k) upto chi(i,j,k)
  trK_rhs(i,j,k) =    gupxx(i,j,k) * fxx(i,j,k) + gupyy(i,j,k) * fyy(i,j,k) + gupzz(i,j,k) * fzz(i,j,k) + &
        TWO* ( gupxy(i,j,k) * fxy(i,j,k) + gupxz(i,j,k) * fxz(i,j,k) + gupyz(i,j,k) * fyz(i,j,k) )
! Add lapse and S_ij parts to Ricci tensor:

  fxx(i,j,k) = alpn1(i,j,k) * (Rxx(i,j,k) - EIGHT * PI * Sxx(i,j,k)) - fxx(i,j,k)
  fxy(i,j,k) = alpn1(i,j,k) * (Rxy(i,j,k) - EIGHT * PI * Sxy(i,j,k)) - fxy(i,j,k)
  fxz(i,j,k) = alpn1(i,j,k) * (Rxz(i,j,k) - EIGHT * PI * Sxz(i,j,k)) - fxz(i,j,k)
  fyy(i,j,k) = alpn1(i,j,k) * (Ryy(i,j,k) - EIGHT * PI * Syy(i,j,k)) - fyy(i,j,k)
  fyz(i,j,k) = alpn1(i,j,k) * (Ryz(i,j,k) - EIGHT * PI * Syz(i,j,k)) - fyz(i,j,k)
  fzz(i,j,k) = alpn1(i,j,k) * (Rzz(i,j,k) - EIGHT * PI * Szz(i,j,k)) - fzz(i,j,k)

! Compute trace-free part (note: chi(i,j,k)^-1 and chi(i,j,k) cancel!):

  f(i,j,k) = F1o3 *(  gupxx(i,j,k) * fxx(i,j,k) + gupyy(i,j,k) * fyy(i,j,k) + gupzz(i,j,k) * fzz(i,j,k) + &
        TWO* ( gupxy(i,j,k) * fxy(i,j,k) + gupxz(i,j,k) * fxz(i,j,k) + gupyz(i,j,k) * fyz(i,j,k) ) )

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
! store D^i D_i Lap(i,j,k) in trK_rhs(i,j,k)
  trK_rhs(i,j,k) = f(i,j,k)*trK_rhs(i,j,k)
          
  Axx_rhs(i,j,k) =           f(i,j,k) * Axx_rhs(i,j,k)+ alpn1(i,j,k) * (trK(i,j,k) * Axx(i,j,k) - TWO * fxx(i,j,k))  + &
           TWO * (  Axx(i,j,k) * betaxx(i,j,k) +   Axy(i,j,k) * betayx(i,j,k) +   Axz(i,j,k) * betazx(i,j,k) )- &
             F2o3 * Axx(i,j,k) * div_beta(i,j,k)

  Ayy_rhs(i,j,k) =           f(i,j,k) * Ayy_rhs(i,j,k)+ alpn1(i,j,k) * (trK(i,j,k) * Ayy(i,j,k) - TWO * fyy(i,j,k))  + &
           TWO * (  Axy(i,j,k) * betaxy(i,j,k) +   Ayy(i,j,k) * betayy(i,j,k) +   Ayz(i,j,k) * betazy(i,j,k) )- &
             F2o3 * Ayy(i,j,k) * div_beta(i,j,k)

  Azz_rhs(i,j,k) =           f(i,j,k) * Azz_rhs(i,j,k)+ alpn1(i,j,k) * (trK(i,j,k) * Azz(i,j,k) - TWO * fzz(i,j,k))  + &
           TWO * (  Axz(i,j,k) * betaxz(i,j,k) +   Ayz(i,j,k) * betayz(i,j,k) +   Azz(i,j,k) * betazz(i,j,k) )- &
             F2o3 * Azz(i,j,k) * div_beta(i,j,k)

  Axy_rhs(i,j,k) =           f(i,j,k) * Axy_rhs(i,j,k)+ alpn1(i,j,k) *( trK(i,j,k) * Axy(i,j,k)  - TWO * fxy(i,j,k) )+ &
                    Axx(i,j,k) * betaxy(i,j,k)                  +   Axz(i,j,k) * betazy(i,j,k)  + &
                                     Ayy(i,j,k) * betayx(i,j,k) +   Ayz(i,j,k) * betazx(i,j,k)  + &
             F1o3 * Axy(i,j,k) * div_beta(i,j,k)                -   Axy(i,j,k) * betazz(i,j,k)

  Ayz_rhs(i,j,k) =           f(i,j,k) * Ayz_rhs(i,j,k)+ alpn1(i,j,k) *( trK(i,j,k) * Ayz(i,j,k)  - TWO * fyz(i,j,k) )+ &
                    Axy(i,j,k) * betaxz(i,j,k) +   Ayy(i,j,k) * betayz(i,j,k)                   + &
                    Axz(i,j,k) * betaxy(i,j,k)                  +   Azz(i,j,k) * betazy(i,j,k)  + &
             F1o3 * Ayz(i,j,k) * div_beta(i,j,k)                -   Ayz(i,j,k) * betaxx(i,j,k)
 
  Axz_rhs(i,j,k) =           f(i,j,k) * Axz_rhs(i,j,k)+ alpn1(i,j,k) *( trK(i,j,k) * Axz(i,j,k)  - TWO * fxz(i,j,k) )+ &
                    Axx(i,j,k) * betaxz(i,j,k) +   Axy(i,j,k) * betayz(i,j,k)                   + &
                                     Ayz(i,j,k) * betayx(i,j,k) +   Azz(i,j,k) * betazx(i,j,k)  + &
             F1o3 * Axz(i,j,k) * div_beta(i,j,k)                -   Axz(i,j,k) * betayy(i,j,k)      !rhs for Aij

! Compute trace of S_ij

  S(i,j,k) =  f(i,j,k) * ( gupxx(i,j,k) * Sxx(i,j,k) + gupyy(i,j,k) * Syy(i,j,k) + gupzz(i,j,k) * Szz(i,j,k) + &
     TWO * ( gupxy(i,j,k) * Sxy(i,j,k) + gupxz(i,j,k) * Sxz(i,j,k) + gupyz(i,j,k) * Syz(i,j,k) ) )

  trK_rhs(i,j,k) = - trK_rhs(i,j,k) + alpn1(i,j,k) *( F1o3 * trK(i,j,k) * trK(i,j,k)         + &
                gupxx(i,j,k) * fxx(i,j,k) + gupyy(i,j,k) * fyy(i,j,k) + gupzz(i,j,k) * fzz(i,j,k)   + &
        TWO * ( gupxy(i,j,k) * fxy(i,j,k) + gupxz(i,j,k) * fxz(i,j,k) + gupyz(i,j,k) * fyz(i,j,k) ) + &
       FOUR * PI * ( rho(i,j,k) + S(i,j,k) ))                                !rhs for trK(i,j,k)
  
!!!! gauge variable part

  Lap_rhs(i,j,k) = -TWO*alpn1(i,j,k)*trK(i,j,k)

  betax_rhs(i,j,k) = FF*dtSfx(i,j,k)
  betay_rhs(i,j,k) = FF*dtSfy(i,j,k)
  betaz_rhs(i,j,k) = FF*dtSfz(i,j,k)

  dtSfx_rhs(i,j,k) = Gamx_rhs(i,j,k) - eta*dtSfx(i,j,k)
  dtSfy_rhs(i,j,k) = Gamy_rhs(i,j,k) - eta*dtSfy(i,j,k)
  dtSfz_rhs(i,j,k) = Gamz_rhs(i,j,k) - eta*dtSfz(i,j,k)
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

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
!!
  call lopsided(ex,X,Y,Z,Lap,Lap_rhs,betax,betay,betaz,Symmetry,SSS)

  call lopsided(ex,X,Y,Z,betax,betax_rhs,betax,betay,betaz,Symmetry,ASS)
  call lopsided(ex,X,Y,Z,betay,betay_rhs,betax,betay,betaz,Symmetry,SAS)
  call lopsided(ex,X,Y,Z,betaz,betaz_rhs,betax,betay,betaz,Symmetry,SSA)

  call lopsided(ex,X,Y,Z,dtSfx,dtSfx_rhs,betax,betay,betaz,Symmetry,ASS)
  call lopsided(ex,X,Y,Z,dtSfy,dtSfy_rhs,betax,betay,betaz,Symmetry,SAS)
  call lopsided(ex,X,Y,Z,dtSfz,dtSfz_rhs,betax,betay,betaz,Symmetry,SSA)
#if 1 
!!
  call lopsided(ex,X,Y,Z,Sphi,Sphi_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,Spi , Spi_rhs,betax,betay,betaz,Symmetry,SSS)
#endif
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
  call kodis(ex,X,Y,Z,dtSfx,dtSfx_rhs,ASS,Symmetry,eps)
  call kodis(ex,X,Y,Z,dtSfy,dtSfy_rhs,SAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,dtSfz,dtSfz_rhs,SSA,Symmetry,eps)

#if 1
  call kodis(ex,X,Y,Z,Sphi,Sphi_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Spi ,Spi_rhs ,SSS,Symmetry,eps)
#endif
  endif

#if 0
   Sphi_rhs = ZEO
   Spi_rhs = ZEO
#endif   

  if(co == 0)then
! ham_Res = trR + 2/3 * K^2 - A_ij * A^ij - 16 * PI * rho
! here trR is respect to physical metric
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
  ham_Res(i,j,k) =   gupxx(i,j,k) * Rxx(i,j,k) + gupyy(i,j,k) * Ryy(i,j,k) + gupzz(i,j,k) * Rzz(i,j,k) + &
        TWO* ( gupxy(i,j,k) * Rxy(i,j,k) + gupxz(i,j,k) * Rxz(i,j,k) + gupyz(i,j,k) * Ryz(i,j,k) )

  ham_Res(i,j,k) = chin1(i,j,k)*ham_Res(i,j,k) + F2o3 * trK(i,j,k) * trK(i,j,k) -(&
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
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

! mov_Res_j = gupkj*(-1/chi d_k chi*A_ij + D_k A_ij) - 2/3 d_j trK - 8 PI s_j where D respect to physical metric
! store D_i A_jk - 1/chi d_i chi*A_jk in gjk_i
  call fderivs(ex,Axx,gxxx,gxxy,gxxz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)
  call fderivs(ex,Axy,gxyx,gxyy,gxyz,X,Y,Z,ANTI,ANTI,SYM ,Symmetry,0)
  call fderivs(ex,Axz,gxzx,gxzy,gxzz,X,Y,Z,ANTI,SYM ,ANTI,Symmetry,0)
  call fderivs(ex,Ayy,gyyx,gyyy,gyyz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)
  call fderivs(ex,Ayz,gyzx,gyzy,gyzz,X,Y,Z,SYM ,ANTI,ANTI,Symmetry,0)
  call fderivs(ex,Azz,gzzx,gzzy,gzzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
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
movx_Res(i,j,k) = gupxx(i,j,k)*gxxx(i,j,k) + gupyy(i,j,k)*gxyy(i,j,k) + gupzz(i,j,k)*gxzz(i,j,k) &
          +gupxy(i,j,k)*gxyx(i,j,k) + gupxz(i,j,k)*gxzx(i,j,k) + gupyz(i,j,k)*gxzy(i,j,k) &
          +gupxy(i,j,k)*gxxy(i,j,k) + gupxz(i,j,k)*gxxz(i,j,k) + gupyz(i,j,k)*gxyz(i,j,k)
movy_Res(i,j,k) = gupxx(i,j,k)*gxyx(i,j,k) + gupyy(i,j,k)*gyyy(i,j,k) + gupzz(i,j,k)*gyzz(i,j,k) &
          +gupxy(i,j,k)*gyyx(i,j,k) + gupxz(i,j,k)*gyzx(i,j,k) + gupyz(i,j,k)*gyzy(i,j,k) &
          +gupxy(i,j,k)*gxyy(i,j,k) + gupxz(i,j,k)*gxyz(i,j,k) + gupyz(i,j,k)*gyyz(i,j,k)
movz_Res(i,j,k) = gupxx(i,j,k)*gxzx(i,j,k) + gupyy(i,j,k)*gyzy(i,j,k) + gupzz(i,j,k)*gzzz(i,j,k) &
          +gupxy(i,j,k)*gyzx(i,j,k) + gupxz(i,j,k)*gzzx(i,j,k) + gupyz(i,j,k)*gzzy(i,j,k) &
          +gupxy(i,j,k)*gxzy(i,j,k) + gupxz(i,j,k)*gxzz(i,j,k) + gupyz(i,j,k)*gyzz(i,j,k)

movx_Res(i,j,k) = movx_Res(i,j,k) - F2o3*Kx(i,j,k) - F8*PI*sx(i,j,k)
movy_Res(i,j,k) = movy_Res(i,j,k) - F2o3*Ky(i,j,k) - F8*PI*sy(i,j,k)
movz_Res(i,j,k) = movz_Res(i,j,k) - F2o3*Kz(i,j,k) - F8*PI*sz(i,j,k)
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO
  endif

  gont = 0

  return

  end function compute_rhs_bssn_escalar
!! for shell part
  function compute_rhs_bssn_escalar_ss(ex, T,crho,sigma,R,x,y,z,               &
               drhodx, drhody, drhodz,                                         &
               dsigmadx,dsigmady,dsigmadz,                                     &
               dRdx,dRdy,dRdz,                                                 &
               drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                &
               dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,    &
               dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,                            &
               chi    ,   trK    ,                                             &
               dxx    ,   gxy    ,   gxz    ,   dyy    ,   gyz    ,   dzz,     &
               Axx    ,   Axy    ,   Axz    ,   Ayy    ,   Ayz    ,   Azz,     &
               Gamx   ,  Gamy    ,  Gamz    ,                                  &
               Lap    ,  betax   ,  betay   ,  betaz   ,                       &
               dtSfx  ,  dtSfy   ,  dtSfz   ,                                  &
               Sphi   ,   Spi    ,                                             &
               chi_rhs,   trK_rhs,                                             &
               gxx_rhs,   gxy_rhs,   gxz_rhs,   gyy_rhs,   gyz_rhs,   gzz_rhs, &
               Axx_rhs,   Axy_rhs,   Axz_rhs,   Ayy_rhs,   Ayz_rhs,   Azz_rhs, &
               Gamx_rhs,  Gamy_rhs,  Gamz_rhs,                                 &
               Lap_rhs,  betax_rhs,  betay_rhs,  betaz_rhs,                    &
               dtSfx_rhs,  dtSfy_rhs,  dtSfz_rhs,                              &
               Sphi_rhs ,   Spi_rhs ,                                          &
               rho,Sx,Sy,Sz,Sxx,Sxy,Sxz,Syy,Syz,Szz,                           &
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,                      &
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,                      &
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,                      &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,                                        &
               ham_Res, movx_Res, movy_Res, movz_Res,                          &
                        Gmx_Res, Gmy_Res, Gmz_Res,                             &
               Symmetry,Lev,eps,sst,co)  result(gont)
! calculate constraint violation when co=0       
  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3), Symmetry,Lev,sst,co
  real*8, intent(in ):: T
  double precision,intent(in),dimension(ex(1))::crho
  double precision,intent(in),dimension(ex(2))::sigma
  double precision,intent(in),dimension(ex(3))::R
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::x,y,z
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::drhodx, drhody, drhodz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::dsigmadx,dsigmady,dsigmadz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::dRdx,dRdy,dRdz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: chi,dxx,dyy,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: trK
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: gxy,gxz,gyz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamx,Gamy,Gamz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Lap, betax, betay, betaz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: dtSfx,  dtSfy,  dtSfz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Sphi,Spi
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: chi_rhs,trK_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gxx_rhs,gxy_rhs,gxz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gyy_rhs,gyz_rhs,gzz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Axx_rhs,Axy_rhs,Axz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Ayy_rhs,Ayz_rhs,Azz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamx_rhs,Gamy_rhs,Gamz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Lap_rhs, betax_rhs, betay_rhs, betaz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: dtSfx_rhs,dtSfy_rhs,dtSfz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Sphi_rhs,Spi_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: rho,Sx,Sy,Sz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Sxx,Sxy,Sxz,Syy,Syz,Szz
! when out, physical second kind of connection  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamzyy, Gamzyz, Gamzzz
! when out, physical Ricci tensor  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Rxx,Rxy,Rxz,Ryy,Ryz,Rzz
  real*8,intent(in) :: eps
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: ham_Res, movx_Res, movy_Res, movz_Res
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gmx_Res, Gmy_Res, Gmz_Res
!  gont = 0: success; gont = 1: something wrong
  integer::gont

!~~~~~~> Other variables:

  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gyy,gzz
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
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamxa,Gamya,Gamza,alpn1,chin1
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupyy,gupyz,gupzz

  real*8,dimension(3) ::SSS,AAS,ASA,SAA,ASS,SAS,SSA
  real*8            :: dX, dY, dZ, PI
  real*8, parameter :: ZEO=0.d0, ONE = 1.D0, TWO = 2.D0, FOUR = 4.D0
  real*8, parameter :: EIGHT = 8.D0, HALF = 0.5D0, THR = 3.d0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  double precision,parameter::FF = 0.75d0,eta=2.d0
  real*8, parameter :: F1o3 = 1.D0/3.D0, F2o3 = 2.D0/3.D0,F3o2=1.5d0, F1o6 = 1.D0/6.D0
  real*8, parameter :: F16=1.6d1,F8=8.d0

  integer :: i,j,k

!!! sanity check
  dX = sum(chi)+sum(trK)+sum(dxx)+sum(gxy)+sum(gxz)+sum(dyy)+sum(gyz)+sum(dzz) &
      +sum(Axx)+sum(Axy)+sum(Axz)+sum(Ayy)+sum(Ayz)+sum(Azz)                   &
      +sum(Gamx)+sum(Gamy)+sum(Gamz)                                           &
      +sum(Lap)+sum(betax)+sum(betay)+sum(betaz)                               &
      +sum(Sphi)+sum(Spi)
  if(dX.ne.dX) then
     if(sum(chi).ne.sum(chi))write(*,*)"bssn.f90: find NaN in chi"
     if(sum(trK).ne.sum(trK))write(*,*)"bssn.f90: find NaN in trk"
     if(sum(dxx).ne.sum(dxx))write(*,*)"bssn.f90: find NaN in dxx"
     if(sum(gxy).ne.sum(gxy))write(*,*)"bssn.f90: find NaN in gxy"
     if(sum(gxz).ne.sum(gxz))write(*,*)"bssn.f90: find NaN in gxz"
     if(sum(dyy).ne.sum(dyy))write(*,*)"bssn.f90: find NaN in dyy"
     if(sum(gyz).ne.sum(gyz))write(*,*)"bssn.f90: find NaN in gyz"
     if(sum(dzz).ne.sum(dzz))write(*,*)"bssn.f90: find NaN in dzz"
     if(sum(Axx).ne.sum(Axx))write(*,*)"bssn.f90: find NaN in Axx"
     if(sum(Axy).ne.sum(Axy))write(*,*)"bssn.f90: find NaN in Axy"
     if(sum(Axz).ne.sum(Axz))write(*,*)"bssn.f90: find NaN in Axz"
     if(sum(Ayy).ne.sum(Ayy))write(*,*)"bssn.f90: find NaN in Ayy"
     if(sum(Ayz).ne.sum(Ayz))write(*,*)"bssn.f90: find NaN in Ayz"
     if(sum(Azz).ne.sum(Azz))write(*,*)"bssn.f90: find NaN in Azz"
     if(sum(Gamx).ne.sum(Gamx))write(*,*)"bssn.f90: find NaN in Gamx"
     if(sum(Gamy).ne.sum(Gamy))write(*,*)"bssn.f90: find NaN in Gamy"
     if(sum(Gamz).ne.sum(Gamz))write(*,*)"bssn.f90: find NaN in Gamz"
     if(sum(Lap).ne.sum(Lap))write(*,*)"bssn.f90: find NaN in Lap"
     if(sum(betax).ne.sum(betax))write(*,*)"bssn.f90: find NaN in betax"
     if(sum(betay).ne.sum(betay))write(*,*)"bssn.f90: find NaN in betay"
     if(sum(betaz).ne.sum(betaz))write(*,*)"bssn.f90: find NaN in betaz"
     if(sum(Sphi).ne.sum(Sphi))write(*,*)"bssn.f90: find NaN in Sphi"
     if(sum(Spi).ne.sum(Spi))write(*,*)"bssn.f90: find NaN in Spi"
     gont = 1
     return
  endif

  PI = dacos(-ONE)

  dX = crho(2) - crho(1)
  dY = sigma(2) - sigma(1)
  dZ = R(2) - R(1)

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
  alpn1(i,j,k) = Lap(i,j,k) + ONE
  chin1(i,j,k) = chi(i,j,k) + ONE
  gxx(i,j,k) = dxx(i,j,k) + ONE
  gyy(i,j,k) = dyy(i,j,k) + ONE
  gzz(i,j,k) = dzz(i,j,k) + ONE
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

  call fderivs_shc(ex,betax,betaxx,betaxy,betaxz,crho,sigma,R,ANTI, SYM, SYM,Symmetry,Lev,sst, &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,betay,betayx,betayy,betayz,crho,sigma,R, SYM,ANTI, SYM,Symmetry,Lev,sst, &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,betaz,betazx,betazy,betazz,crho,sigma,R, SYM, SYM,ANTI,Symmetry,Lev,sst, &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  
  div_beta = betaxx + betayy + betazz
 
  call fderivs_shc(ex,chi,chix,chiy,chiz,crho,sigma,R, SYM, SYM,SYM,Symmetry,Lev,sst,          &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)

  chi_rhs = F2o3 *chin1*( alpn1 * trK - div_beta ) !rhs for chi

  call fderivs_shc(ex,Lap,Lapx,Lapy,Lapz,crho,sigma,R, SYM, SYM,SYM,Symmetry,Lev,sst,          &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
! invert tilted metric
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
  gupzz(i,j,k) =  gxx(i,j,k) * gyy(i,j,k) * gzz(i,j,k) + gxy(i,j,k) * gyz(i,j,k) * gxz(i,j,k) + gxz(i,j,k) * gxy(i,j,k) * gyz(i,j,k) - &
           gxz(i,j,k) * gyy(i,j,k) * gxz(i,j,k) - gxy(i,j,k) * gxy(i,j,k) * gzz(i,j,k) - gxx(i,j,k) * gyz(i,j,k) * gyz(i,j,k)
  gupxx(i,j,k) =   ( gyy(i,j,k) * gzz(i,j,k) - gyz(i,j,k) * gyz(i,j,k) ) / gupzz(i,j,k)
  gupxy(i,j,k) = - ( gxy(i,j,k) * gzz(i,j,k) - gyz(i,j,k) * gxz(i,j,k) ) / gupzz(i,j,k)
  gupxz(i,j,k) =   ( gxy(i,j,k) * gyz(i,j,k) - gyy(i,j,k) * gxz(i,j,k) ) / gupzz(i,j,k)
  gupyy(i,j,k) =   ( gxx(i,j,k) * gzz(i,j,k) - gxz(i,j,k) * gxz(i,j,k) ) / gupzz(i,j,k)
  gupyz(i,j,k) = - ( gxx(i,j,k) * gyz(i,j,k) - gxy(i,j,k) * gxz(i,j,k) ) / gupzz(i,j,k)
  gupzz(i,j,k) =   ( gxx(i,j,k) * gyy(i,j,k) - gxy(i,j,k) * gxy(i,j,k) ) / gupzz(i,j,k)
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

#if 1 
  Sphi_rhs = alpn1 * Spi !rhs for Scalar phi
!rhs for Spi
  call fderivs_shc(ex,Sphi,Kx,Ky,Kz,crho,sigma,R, SYM, SYM,SYM,Symmetry,Lev,sst,          &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fdderivs_shc(ex,Sphi,fxx,fxy,fxz,fyy,fyz,fzz,crho,sigma,R, SYM, SYM,SYM ,Symmetry,Lev,sst,  &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz)
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
   Spi_rhs(i,j,k) =   gupxx(i,j,k) * fxx(i,j,k) + gupyy(i,j,k) * fyy(i,j,k) + gupzz(i,j,k) * fzz(i,j,k) +           &
             ( gupxy(i,j,k) * fxy(i,j,k) + gupxz(i,j,k) * fxz(i,j,k) + gupyz(i,j,k) * fyz(i,j,k) ) * TWO -   &
             ((Gamx(i,j,k)+(gupxx(i,j,k)*chix(i,j,k)+gupxy(i,j,k)*chiy(i,j,k)+gupxz(i,j,k)*chiz(i,j,k))/TWO/chin1(i,j,k))*Kx(i,j,k)  &
            + (Gamy(i,j,k)+(gupxy(i,j,k)*chix(i,j,k)+gupyy(i,j,k)*chiy(i,j,k)+gupyz(i,j,k)*chiz(i,j,k))/TWO/chin1(i,j,k))*Ky(i,j,k)  &
            + (Gamz(i,j,k)+(gupxz(i,j,k)*chix(i,j,k)+gupyz(i,j,k)*chiy(i,j,k)+gupzz(i,j,k)*chiz(i,j,k))/TWO/chin1(i,j,k))*Kz(i,j,k))
   Spi_rhs(i,j,k) = Spi_rhs(i,j,k)*alpn1(i,j,k) +                                       &
            (gupxx(i,j,k)*Lapx(i,j,k)*Kx(i,j,k) + gupxy(i,j,k)*Lapx(i,j,k)*Ky(i,j,k) + gupxz(i,j,k)*Lapx(i,j,k)*Kz(i,j,k)&
            +gupxy(i,j,k)*Lapy(i,j,k)*Kx(i,j,k) + gupyy(i,j,k)*Lapy(i,j,k)*Ky(i,j,k) + gupyz(i,j,k)*Lapy(i,j,k)*Kz(i,j,k)&
            +gupxz(i,j,k)*Lapz(i,j,k)*Kx(i,j,k) + gupyz(i,j,k)*Lapz(i,j,k)*Ky(i,j,k) + gupzz(i,j,k)*Lapz(i,j,k)*Kz(i,j,k))
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO
 
   call frpotential(ex,Sphi,f,S)
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
   Spi_rhs(i,j,k) = Spi_rhs(i,j,k)*chin1(i,j,k) + alpn1(i,j,k)*(trK(i,j,k)*Spi(i,j,k) - S(i,j,k))
! matter content of scalar field T_ab 
  rho(i,j,k) = chin1(i,j,k)*((gupxx(i,j,k) * Kx(i,j,k) * Kx(i,j,k) + gupyy(i,j,k) * Ky(i,j,k) * Ky(i,j,k) + gupzz(i,j,k) * Kz(i,j,k) * Kz(i,j,k))/TWO +  &
                gupxy(i,j,k) * Kx(i,j,k) * Ky(i,j,k) + gupxz(i,j,k) * Kx(i,j,k) * Kz(i,j,k) + gupyz(i,j,k) * Ky(i,j,k) * Kz(i,j,k) ) &
        + Spi(i,j,k)*Spi(i,j,k)/TWO+f(i,j,k)  
  Sx(i,j,k) = -Spi(i,j,k)*Kx(i,j,k)
  Sy(i,j,k) = -Spi(i,j,k)*Ky(i,j,k)
  Sz(i,j,k) = -Spi(i,j,k)*Kz(i,j,k)
  f(i,j,k) = (rho(i,j,k) - Spi(i,j,k)*Spi(i,j,k))/chin1(i,j,k)
  Sxx(i,j,k) = Kx(i,j,k)*Kx(i,j,k)-f(i,j,k)*gxx(i,j,k)
  Sxy(i,j,k) = Kx(i,j,k)*Ky(i,j,k)-f(i,j,k)*gxy(i,j,k)
  Sxz(i,j,k) = Kx(i,j,k)*Kz(i,j,k)-f(i,j,k)*gxz(i,j,k)
  Syy(i,j,k) = Ky(i,j,k)*Ky(i,j,k)-f(i,j,k)*gyy(i,j,k)
  Syz(i,j,k) = Ky(i,j,k)*Kz(i,j,k)-f(i,j,k)*gyz(i,j,k)
  Szz(i,j,k) = Kz(i,j,k)*Kz(i,j,k)-f(i,j,k)*gzz(i,j,k)
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO
#else
  rho = ZEO
  Sx = ZEO
  Sy = ZEO
  Sz = ZEO
  Sxx = ZEO
  Sxy = ZEO
  Sxz = ZEO
  Syy = ZEO
  Syz = ZEO
  Szz = ZEO
#endif  

  call fderivs_shc(ex,dxx,gxxx,gxxy,gxxz,crho,sigma,R, SYM, SYM,SYM,Symmetry,Lev,sst,          &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,gxy,gxyx,gxyy,gxyz,crho,sigma,R,ANTI,ANTI,SYM,Symmetry,Lev,sst,          &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,gxz,gxzx,gxzy,gxzz,crho,sigma,R,ANTI,SYM ,ANTI,Symmetry,Lev,sst,         &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,dyy,gyyx,gyyy,gyyz,crho,sigma,R, SYM, SYM,SYM,Symmetry,Lev,sst,          &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,gyz,gyzx,gyzy,gyzz,crho,sigma,R,SYM ,ANTI,ANTI,Symmetry,Lev,sst,         &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,dzz,gzzx,gzzy,gzzz,crho,sigma,R, SYM, SYM,SYM,Symmetry,Lev,sst,          &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)

  if(co == 0)then
! Gam^i_Res = Gam^i + gup^ij_,j
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
  Gmx_Res(i,j,k) = Gamx(i,j,k) - (gupxx(i,j,k)*(gupxx(i,j,k)*gxxx(i,j,k)+gupxy(i,j,k)*gxyx(i,j,k)+gupxz(i,j,k)*gxzx(i,j,k))&
                   +gupxy(i,j,k)*(gupxx(i,j,k)*gxyx(i,j,k)+gupxy(i,j,k)*gyyx(i,j,k)+gupxz(i,j,k)*gyzx(i,j,k))&
                   +gupxz(i,j,k)*(gupxx(i,j,k)*gxzx(i,j,k)+gupxy(i,j,k)*gyzx(i,j,k)+gupxz(i,j,k)*gzzx(i,j,k))&
                   +gupxx(i,j,k)*(gupxy(i,j,k)*gxxy(i,j,k)+gupyy(i,j,k)*gxyy(i,j,k)+gupyz(i,j,k)*gxzy(i,j,k))&
                   +gupxy(i,j,k)*(gupxy(i,j,k)*gxyy(i,j,k)+gupyy(i,j,k)*gyyy(i,j,k)+gupyz(i,j,k)*gyzy(i,j,k))&
                   +gupxz(i,j,k)*(gupxy(i,j,k)*gxzy(i,j,k)+gupyy(i,j,k)*gyzy(i,j,k)+gupyz(i,j,k)*gzzy(i,j,k))&
                   +gupxx(i,j,k)*(gupxz(i,j,k)*gxxz(i,j,k)+gupyz(i,j,k)*gxyz(i,j,k)+gupzz(i,j,k)*gxzz(i,j,k))&
                   +gupxy(i,j,k)*(gupxz(i,j,k)*gxyz(i,j,k)+gupyz(i,j,k)*gyyz(i,j,k)+gupzz(i,j,k)*gyzz(i,j,k))&
                   +gupxz(i,j,k)*(gupxz(i,j,k)*gxzz(i,j,k)+gupyz(i,j,k)*gyzz(i,j,k)+gupzz(i,j,k)*gzzz(i,j,k)))
  Gmy_Res(i,j,k) = Gamy(i,j,k) - (gupxx(i,j,k)*(gupxy(i,j,k)*gxxx(i,j,k)+gupyy(i,j,k)*gxyx(i,j,k)+gupyz(i,j,k)*gxzx(i,j,k))&
                   +gupxy(i,j,k)*(gupxy(i,j,k)*gxyx(i,j,k)+gupyy(i,j,k)*gyyx(i,j,k)+gupyz(i,j,k)*gyzx(i,j,k))&
                   +gupxz(i,j,k)*(gupxy(i,j,k)*gxzx(i,j,k)+gupyy(i,j,k)*gyzx(i,j,k)+gupyz(i,j,k)*gzzx(i,j,k))&
                   +gupxy(i,j,k)*(gupxy(i,j,k)*gxxy(i,j,k)+gupyy(i,j,k)*gxyy(i,j,k)+gupyz(i,j,k)*gxzy(i,j,k))&
                   +gupyy(i,j,k)*(gupxy(i,j,k)*gxyy(i,j,k)+gupyy(i,j,k)*gyyy(i,j,k)+gupyz(i,j,k)*gyzy(i,j,k))&
                   +gupyz(i,j,k)*(gupxy(i,j,k)*gxzy(i,j,k)+gupyy(i,j,k)*gyzy(i,j,k)+gupyz(i,j,k)*gzzy(i,j,k))&
                   +gupxy(i,j,k)*(gupxz(i,j,k)*gxxz(i,j,k)+gupyz(i,j,k)*gxyz(i,j,k)+gupzz(i,j,k)*gxzz(i,j,k))&
                   +gupyy(i,j,k)*(gupxz(i,j,k)*gxyz(i,j,k)+gupyz(i,j,k)*gyyz(i,j,k)+gupzz(i,j,k)*gyzz(i,j,k))&
                   +gupyz(i,j,k)*(gupxz(i,j,k)*gxzz(i,j,k)+gupyz(i,j,k)*gyzz(i,j,k)+gupzz(i,j,k)*gzzz(i,j,k)))
  Gmz_Res(i,j,k) = Gamz(i,j,k) - (gupxx(i,j,k)*(gupxz(i,j,k)*gxxx(i,j,k)+gupyz(i,j,k)*gxyx(i,j,k)+gupzz(i,j,k)*gxzx(i,j,k))&
                   +gupxy(i,j,k)*(gupxz(i,j,k)*gxyx(i,j,k)+gupyz(i,j,k)*gyyx(i,j,k)+gupzz(i,j,k)*gyzx(i,j,k))&
                   +gupxz(i,j,k)*(gupxz(i,j,k)*gxzx(i,j,k)+gupyz(i,j,k)*gyzx(i,j,k)+gupzz(i,j,k)*gzzx(i,j,k))&
                   +gupxy(i,j,k)*(gupxz(i,j,k)*gxxy(i,j,k)+gupyz(i,j,k)*gxyy(i,j,k)+gupzz(i,j,k)*gxzy(i,j,k))&
                   +gupyy(i,j,k)*(gupxz(i,j,k)*gxyy(i,j,k)+gupyz(i,j,k)*gyyy(i,j,k)+gupzz(i,j,k)*gyzy(i,j,k))&
                   +gupyz(i,j,k)*(gupxz(i,j,k)*gxzy(i,j,k)+gupyz(i,j,k)*gyzy(i,j,k)+gupzz(i,j,k)*gzzy(i,j,k))&
                   +gupxz(i,j,k)*(gupxz(i,j,k)*gxxz(i,j,k)+gupyz(i,j,k)*gxyz(i,j,k)+gupzz(i,j,k)*gxzz(i,j,k))&
                   +gupyz(i,j,k)*(gupxz(i,j,k)*gxyz(i,j,k)+gupyz(i,j,k)*gyyz(i,j,k)+gupzz(i,j,k)*gyzz(i,j,k))&
                   +gupzz(i,j,k)*(gupxz(i,j,k)*gxzz(i,j,k)+gupyz(i,j,k)*gyzz(i,j,k)+gupzz(i,j,k)*gzzz(i,j,k)))
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO
  endif

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
  gxx_rhs(i,j,k) = - TWO * alpn1(i,j,k) * Axx(i,j,k)    -  F2o3 * gxx(i,j,k) * div_beta(i,j,k)          + &
              TWO *(  gxx(i,j,k) * betaxx(i,j,k) +   gxy(i,j,k) * betayx(i,j,k) +   gxz(i,j,k) * betazx(i,j,k))

  gyy_rhs(i,j,k) = - TWO * alpn1(i,j,k) * Ayy(i,j,k)    -  F2o3 * gyy(i,j,k) * div_beta(i,j,k)          + &
              TWO *(  gxy(i,j,k) * betaxy(i,j,k) +   gyy(i,j,k) * betayy(i,j,k) +   gyz(i,j,k) * betazy(i,j,k))

  gzz_rhs(i,j,k) = - TWO * alpn1(i,j,k) * Azz(i,j,k)    -  F2o3 * gzz(i,j,k) * div_beta(i,j,k)          + &
              TWO *(  gxz(i,j,k) * betaxz(i,j,k) +   gyz(i,j,k) * betayz(i,j,k) +   gzz(i,j,k) * betazz(i,j,k))

  gxy_rhs(i,j,k) = - TWO * alpn1(i,j,k) * Axy(i,j,k)    +  F1o3 * gxy(i,j,k)    * div_beta(i,j,k)       + &
                      gxx(i,j,k) * betaxy(i,j,k)                  +   gxz(i,j,k) * betazy(i,j,k) + &
                                       gyy(i,j,k) * betayx(i,j,k) +   gyz(i,j,k) * betazx(i,j,k)   &
                                                    -   gxy(i,j,k) * betazz(i,j,k)

  gyz_rhs(i,j,k) = - TWO * alpn1(i,j,k) * Ayz(i,j,k)    +  F1o3 * gyz(i,j,k)    * div_beta(i,j,k)       + &
                      gxy(i,j,k) * betaxz(i,j,k) +   gyy(i,j,k) * betayz(i,j,k)                  + &
                      gxz(i,j,k) * betaxy(i,j,k)                  +   gzz(i,j,k) * betazy(i,j,k)   &
                                                    -   gyz(i,j,k) * betaxx(i,j,k)
 
  gxz_rhs(i,j,k) = - TWO * alpn1(i,j,k) * Axz(i,j,k)    +  F1o3 * gxz(i,j,k)    * div_beta(i,j,k)       + &
                      gxx(i,j,k) * betaxz(i,j,k) +   gxy(i,j,k) * betayz(i,j,k)                  + &
                                       gyz(i,j,k) * betayx(i,j,k) +   gzz(i,j,k) * betazx(i,j,k)   &
                                                    -   gxz(i,j,k) * betayy(i,j,k)     !rhs for gij

! second kind of connection
  Gamxxx(i,j,k) =HALF*( gupxx(i,j,k)*gxxx(i,j,k) + gupxy(i,j,k)*(TWO*gxyx(i,j,k) - gxxy(i,j,k) ) + gupxz(i,j,k)*(TWO*gxzx(i,j,k) - gxxz(i,j,k) ))
  Gamyxx(i,j,k) =HALF*( gupxy(i,j,k)*gxxx(i,j,k) + gupyy(i,j,k)*(TWO*gxyx(i,j,k) - gxxy(i,j,k) ) + gupyz(i,j,k)*(TWO*gxzx(i,j,k) - gxxz(i,j,k) ))
  Gamzxx(i,j,k) =HALF*( gupxz(i,j,k)*gxxx(i,j,k) + gupyz(i,j,k)*(TWO*gxyx(i,j,k) - gxxy(i,j,k) ) + gupzz(i,j,k)*(TWO*gxzx(i,j,k) - gxxz(i,j,k) ))
 
  Gamxyy(i,j,k) =HALF*( gupxx(i,j,k)*(TWO*gxyy(i,j,k) - gyyx(i,j,k) ) + gupxy(i,j,k)*gyyy(i,j,k) + gupxz(i,j,k)*(TWO*gyzy(i,j,k) - gyyz(i,j,k) ))
  Gamyyy(i,j,k) =HALF*( gupxy(i,j,k)*(TWO*gxyy(i,j,k) - gyyx(i,j,k) ) + gupyy(i,j,k)*gyyy(i,j,k) + gupyz(i,j,k)*(TWO*gyzy(i,j,k) - gyyz(i,j,k) ))
  Gamzyy(i,j,k) =HALF*( gupxz(i,j,k)*(TWO*gxyy(i,j,k) - gyyx(i,j,k) ) + gupyz(i,j,k)*gyyy(i,j,k) + gupzz(i,j,k)*(TWO*gyzy(i,j,k) - gyyz(i,j,k) ))

  Gamxzz(i,j,k) =HALF*( gupxx(i,j,k)*(TWO*gxzz(i,j,k) - gzzx(i,j,k) ) + gupxy(i,j,k)*(TWO*gyzz(i,j,k) - gzzy(i,j,k) ) + gupxz(i,j,k)*gzzz(i,j,k))
  Gamyzz(i,j,k) =HALF*( gupxy(i,j,k)*(TWO*gxzz(i,j,k) - gzzx(i,j,k) ) + gupyy(i,j,k)*(TWO*gyzz(i,j,k) - gzzy(i,j,k) ) + gupyz(i,j,k)*gzzz(i,j,k))
  Gamzzz(i,j,k) =HALF*( gupxz(i,j,k)*(TWO*gxzz(i,j,k) - gzzx(i,j,k) ) + gupyz(i,j,k)*(TWO*gyzz(i,j,k) - gzzy(i,j,k) ) + gupzz(i,j,k)*gzzz(i,j,k))

  Gamxxy(i,j,k) =HALF*( gupxx(i,j,k)*gxxy(i,j,k) + gupxy(i,j,k)*gyyx(i,j,k) + gupxz(i,j,k)*( gxzy(i,j,k) + gyzx(i,j,k) - gxyz(i,j,k) ) )
  Gamyxy(i,j,k) =HALF*( gupxy(i,j,k)*gxxy(i,j,k) + gupyy(i,j,k)*gyyx(i,j,k) + gupyz(i,j,k)*( gxzy(i,j,k) + gyzx(i,j,k) - gxyz(i,j,k) ) )
  Gamzxy(i,j,k) =HALF*( gupxz(i,j,k)*gxxy(i,j,k) + gupyz(i,j,k)*gyyx(i,j,k) + gupzz(i,j,k)*( gxzy(i,j,k) + gyzx(i,j,k) - gxyz(i,j,k) ) )

  Gamxxz(i,j,k) =HALF*( gupxx(i,j,k)*gxxz(i,j,k) + gupxy(i,j,k)*( gxyz(i,j,k) + gyzx(i,j,k) - gxzy(i,j,k) ) + gupxz(i,j,k)*gzzx(i,j,k) )
  Gamyxz(i,j,k) =HALF*( gupxy(i,j,k)*gxxz(i,j,k) + gupyy(i,j,k)*( gxyz(i,j,k) + gyzx(i,j,k) - gxzy(i,j,k) ) + gupyz(i,j,k)*gzzx(i,j,k) )
  Gamzxz(i,j,k) =HALF*( gupxz(i,j,k)*gxxz(i,j,k) + gupyz(i,j,k)*( gxyz(i,j,k) + gyzx(i,j,k) - gxzy(i,j,k) ) + gupzz(i,j,k)*gzzx(i,j,k) )

  Gamxyz(i,j,k) =HALF*( gupxx(i,j,k)*( gxyz(i,j,k) + gxzy(i,j,k) - gyzx(i,j,k) ) + gupxy(i,j,k)*gyyz(i,j,k) + gupxz(i,j,k)*gzzy(i,j,k) )
  Gamyyz(i,j,k) =HALF*( gupxy(i,j,k)*( gxyz(i,j,k) + gxzy(i,j,k) - gyzx(i,j,k) ) + gupyy(i,j,k)*gyyz(i,j,k) + gupyz(i,j,k)*gzzy(i,j,k) )
  Gamzyz(i,j,k) =HALF*( gupxz(i,j,k)*( gxyz(i,j,k) + gxzy(i,j,k) - gyzx(i,j,k) ) + gupyz(i,j,k)*gyyz(i,j,k) + gupzz(i,j,k)*gzzy(i,j,k) )
! Raise indices of \tilde A_{ij} and store in R_ij

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
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

! Right hand side for Gam^i without shift terms...
  call fderivs_shc(ex,trK,Kx,Ky,Kz,crho,sigma,R, SYM, SYM,SYM,Symmetry,Lev,sst,                &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
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
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

  call fdderivs_shc(ex,betax,gxxx,gxyx,gxzx,gyyx,gyzx,gzzx,crho,sigma,R,ANTI, SYM, SYM,Symmetry,Lev,sst,&
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz)
  call fdderivs_shc(ex,betay,gxxy,gxyy,gxzy,gyyy,gyzy,gzzy,crho,sigma,R, SYM,ANTI, SYM,Symmetry,Lev,sst,&
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz)
  call fdderivs_shc(ex,betaz,gxxz,gxyz,gxzz,gyyz,gyzz,gzzz,crho,sigma,R, SYM, SYM,ANTI,Symmetry,Lev,sst,&
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz)

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
  fxx(i,j,k) = gxxx(i,j,k) + gxyy(i,j,k) + gxzz(i,j,k)
  fxy(i,j,k) = gxyx(i,j,k) + gyyy(i,j,k) + gyzz(i,j,k)
  fxz(i,j,k) = gxzx(i,j,k) + gyzy(i,j,k) + gzzz(i,j,k)

  Gamxa(i,j,k) =       gupxx(i,j,k) * Gamxxx(i,j,k) + gupyy(i,j,k) * Gamxyy(i,j,k) + gupzz(i,j,k) * Gamxzz(i,j,k) + &
          TWO*( gupxy(i,j,k) * Gamxxy(i,j,k) + gupxz(i,j,k) * Gamxxz(i,j,k) + gupyz(i,j,k) * Gamxyz(i,j,k) )
  Gamya(i,j,k) =       gupxx(i,j,k) * Gamyxx(i,j,k) + gupyy(i,j,k) * Gamyyy(i,j,k) + gupzz(i,j,k) * Gamyzz(i,j,k) + &
          TWO*( gupxy(i,j,k) * Gamyxy(i,j,k) + gupxz(i,j,k) * Gamyxz(i,j,k) + gupyz(i,j,k) * Gamyyz(i,j,k) )
  Gamza(i,j,k) =       gupxx(i,j,k) * Gamzxx(i,j,k) + gupyy(i,j,k) * Gamzyy(i,j,k) + gupzz(i,j,k) * Gamzzz(i,j,k) + &
          TWO*( gupxy(i,j,k) * Gamzxy(i,j,k) + gupxz(i,j,k) * Gamzxz(i,j,k) + gupyz(i,j,k) * Gamzyz(i,j,k) )
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

  call fderivs_shc(ex,Gamx,Gamxx,Gamxy,Gamxz,crho,sigma,R,ANTI,SYM ,SYM,Symmetry,Lev,sst,      &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,Gamy,Gamyx,Gamyy,Gamyz,crho,sigma,R,SYM ,ANTI,SYM,Symmetry,Lev,sst,      &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,Gamz,Gamzx,Gamzy,Gamzz,crho,sigma,R,SYM ,SYM ,ANTI,Symmetry,Lev,sst,     &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
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

!first kind of connection stored in gij,k
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
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

!compute Ricci tensor for tilted metric
  call fdderivs_shc(ex,dxx,fxx,fxy,fxz,fyy,fyz,fzz,crho,sigma,R, SYM, SYM,SYM ,Symmetry,Lev,sst,  &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz)
   Rxx =   gupxx * fxx + gupyy * fyy + gupzz * fzz + &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) * TWO

  call fdderivs_shc(ex,dyy,fxx,fxy,fxz,fyy,fyz,fzz,crho,sigma,R, SYM, SYM,SYM ,Symmetry,Lev,sst,  &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz)
   Ryy =   gupxx * fxx + gupyy * fyy + gupzz * fzz + &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) * TWO

  call fdderivs_shc(ex,dzz,fxx,fxy,fxz,fyy,fyz,fzz,crho,sigma,R, SYM, SYM,SYM ,Symmetry,Lev,sst,  &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz)
   Rzz =   gupxx * fxx + gupyy * fyy + gupzz * fzz + &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) * TWO

  call fdderivs_shc(ex,gxy,fxx,fxy,fxz,fyy,fyz,fzz,crho,sigma,R,ANTI,ANTI,SYM ,Symmetry,Lev,sst,  &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz)
   Rxy =   gupxx * fxx + gupyy * fyy + gupzz * fzz + &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) * TWO

  call fdderivs_shc(ex,gxz,fxx,fxy,fxz,fyy,fyz,fzz,crho,sigma,R,ANTI,SYM ,ANTI,Symmetry,Lev,sst,  &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz)
   Rxz =   gupxx * fxx + gupyy * fyy + gupzz * fzz + &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) * TWO

  call fdderivs_shc(ex,gyz,fxx,fxy,fxz,fyy,fyz,fzz,crho,sigma,R,SYM ,ANTI,ANTI,Symmetry,Lev,sst,  &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz)
   Ryz =   gupxx * fxx + gupyy * fyy + gupzz * fzz + &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) * TWO

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
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
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO
!covariant second derivative of chi respect to tilted metric
  call fdderivs_shc(ex,chi,fxx,fxy,fxz,fyy,fyz,fzz,crho,sigma,R,SYM ,SYM ,SYM ,Symmetry,Lev,sst,  &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz)

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
  fxx(i,j,k) = fxx(i,j,k) - Gamxxx(i,j,k) * chix(i,j,k) - Gamyxx(i,j,k) * chiy(i,j,k) - Gamzxx(i,j,k) * chiz(i,j,k)
  fxy(i,j,k) = fxy(i,j,k) - Gamxxy(i,j,k) * chix(i,j,k) - Gamyxy(i,j,k) * chiy(i,j,k) - Gamzxy(i,j,k) * chiz(i,j,k)
  fxz(i,j,k) = fxz(i,j,k) - Gamxxz(i,j,k) * chix(i,j,k) - Gamyxz(i,j,k) * chiy(i,j,k) - Gamzxz(i,j,k) * chiz(i,j,k)
  fyy(i,j,k) = fyy(i,j,k) - Gamxyy(i,j,k) * chix(i,j,k) - Gamyyy(i,j,k) * chiy(i,j,k) - Gamzyy(i,j,k) * chiz(i,j,k)
  fyz(i,j,k) = fyz(i,j,k) - Gamxyz(i,j,k) * chix(i,j,k) - Gamyyz(i,j,k) * chiy(i,j,k) - Gamzyz(i,j,k) * chiz(i,j,k)
  fzz(i,j,k) = fzz(i,j,k) - Gamxzz(i,j,k) * chix(i,j,k) - Gamyzz(i,j,k) * chiy(i,j,k) - Gamzzz(i,j,k) * chiz(i,j,k)
! Store D^l D_l chi(i,j,k) - 3/(2*chi(i,j,k)) D^l chi(i,j,k) D_l chi(i,j,k) in f(i,j,k)

  f(i,j,k) =        gupxx(i,j,k) * ( fxx(i,j,k) - F3o2/chin1(i,j,k) * chix(i,j,k) * chix(i,j,k) ) + &
             gupyy(i,j,k) * ( fyy(i,j,k) - F3o2/chin1(i,j,k) * chiy(i,j,k) * chiy(i,j,k) ) + &
             gupzz(i,j,k) * ( fzz(i,j,k) - F3o2/chin1(i,j,k) * chiz(i,j,k) * chiz(i,j,k) ) + &
       TWO * gupxy(i,j,k) * ( fxy(i,j,k) - F3o2/chin1(i,j,k) * chix(i,j,k) * chiy(i,j,k) ) + &
       TWO * gupxz(i,j,k) * ( fxz(i,j,k) - F3o2/chin1(i,j,k) * chix(i,j,k) * chiz(i,j,k) ) + &
       TWO * gupyz(i,j,k) * ( fyz(i,j,k) - F3o2/chin1(i,j,k) * chiy(i,j,k) * chiz(i,j,k) ) 
! Add chi(i,j,k) part to Ricci tensor:

  Rxx(i,j,k) = Rxx(i,j,k) + (fxx(i,j,k) - chix(i,j,k)*chix(i,j,k)/chin1(i,j,k)/TWO + gxx(i,j,k) * f(i,j,k))/chin1(i,j,k)/TWO
  Ryy(i,j,k) = Ryy(i,j,k) + (fyy(i,j,k) - chiy(i,j,k)*chiy(i,j,k)/chin1(i,j,k)/TWO + gyy(i,j,k) * f(i,j,k))/chin1(i,j,k)/TWO
  Rzz(i,j,k) = Rzz(i,j,k) + (fzz(i,j,k) - chiz(i,j,k)*chiz(i,j,k)/chin1(i,j,k)/TWO + gzz(i,j,k) * f(i,j,k))/chin1(i,j,k)/TWO
  Rxy(i,j,k) = Rxy(i,j,k) + (fxy(i,j,k) - chix(i,j,k)*chiy(i,j,k)/chin1(i,j,k)/TWO + gxy(i,j,k) * f(i,j,k))/chin1(i,j,k)/TWO
  Rxz(i,j,k) = Rxz(i,j,k) + (fxz(i,j,k) - chix(i,j,k)*chiz(i,j,k)/chin1(i,j,k)/TWO + gxz(i,j,k) * f(i,j,k))/chin1(i,j,k)/TWO
  Ryz(i,j,k) = Ryz(i,j,k) + (fyz(i,j,k) - chiy(i,j,k)*chiz(i,j,k)/chin1(i,j,k)/TWO + gyz(i,j,k) * f(i,j,k))/chin1(i,j,k)/TWO

! now prepare to get physical second kind of connection
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
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

! covariant second derivatives of the lapse respect to physical metric
  call fdderivs_shc(ex,Lap,fxx,fxy,fxz,fyy,fyz,fzz,crho,sigma,R,SYM ,SYM ,SYM ,Symmetry,Lev,sst,  &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz)

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
  fxx(i,j,k) = fxx(i,j,k) - Gamxxx(i,j,k)*Lapx(i,j,k) - Gamyxx(i,j,k)*Lapy(i,j,k) - Gamzxx(i,j,k)*Lapz(i,j,k)
  fyy(i,j,k) = fyy(i,j,k) - Gamxyy(i,j,k)*Lapx(i,j,k) - Gamyyy(i,j,k)*Lapy(i,j,k) - Gamzyy(i,j,k)*Lapz(i,j,k)
  fzz(i,j,k) = fzz(i,j,k) - Gamxzz(i,j,k)*Lapx(i,j,k) - Gamyzz(i,j,k)*Lapy(i,j,k) - Gamzzz(i,j,k)*Lapz(i,j,k)
  fxy(i,j,k) = fxy(i,j,k) - Gamxxy(i,j,k)*Lapx(i,j,k) - Gamyxy(i,j,k)*Lapy(i,j,k) - Gamzxy(i,j,k)*Lapz(i,j,k)
  fxz(i,j,k) = fxz(i,j,k) - Gamxxz(i,j,k)*Lapx(i,j,k) - Gamyxz(i,j,k)*Lapy(i,j,k) - Gamzxz(i,j,k)*Lapz(i,j,k)
  fyz(i,j,k) = fyz(i,j,k) - Gamxyz(i,j,k)*Lapx(i,j,k) - Gamyyz(i,j,k)*Lapy(i,j,k) - Gamzyz(i,j,k)*Lapz(i,j,k)

! mater content only appears in the right hand side of trK(i,j,k) and A_ij which are
! comming now
! store D^i D_i Lap(i,j,k) in trK_rhs(i,j,k) upto chi(i,j,k)
  trK_rhs(i,j,k) =    gupxx(i,j,k) * fxx(i,j,k) + gupyy(i,j,k) * fyy(i,j,k) + gupzz(i,j,k) * fzz(i,j,k) + &
        TWO* ( gupxy(i,j,k) * fxy(i,j,k) + gupxz(i,j,k) * fxz(i,j,k) + gupyz(i,j,k) * fyz(i,j,k) )
! Add lapse and S_ij parts to Ricci tensor:

  fxx(i,j,k) = alpn1(i,j,k) * (Rxx(i,j,k) - EIGHT * PI * Sxx(i,j,k)) - fxx(i,j,k)
  fxy(i,j,k) = alpn1(i,j,k) * (Rxy(i,j,k) - EIGHT * PI * Sxy(i,j,k)) - fxy(i,j,k)
  fxz(i,j,k) = alpn1(i,j,k) * (Rxz(i,j,k) - EIGHT * PI * Sxz(i,j,k)) - fxz(i,j,k)
  fyy(i,j,k) = alpn1(i,j,k) * (Ryy(i,j,k) - EIGHT * PI * Syy(i,j,k)) - fyy(i,j,k)
  fyz(i,j,k) = alpn1(i,j,k) * (Ryz(i,j,k) - EIGHT * PI * Syz(i,j,k)) - fyz(i,j,k)
  fzz(i,j,k) = alpn1(i,j,k) * (Rzz(i,j,k) - EIGHT * PI * Szz(i,j,k)) - fzz(i,j,k)

! Compute trace-free part (note: chi(i,j,k)^-1 and chi(i,j,k) cancel!):

  f(i,j,k) = F1o3 *(  gupxx(i,j,k) * fxx(i,j,k) + gupyy(i,j,k) * fyy(i,j,k) + gupzz(i,j,k) * fzz(i,j,k) + &
        TWO* ( gupxy(i,j,k) * fxy(i,j,k) + gupxz(i,j,k) * fxz(i,j,k) + gupyz(i,j,k) * fyz(i,j,k) ) )

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
! store D^i D_i Lap(i,j,k) in trK_rhs(i,j,k)
  trK_rhs(i,j,k) = f(i,j,k)*trK_rhs(i,j,k)
          
  Axx_rhs(i,j,k) =           f(i,j,k) * Axx_rhs(i,j,k)+ alpn1(i,j,k) * (trK(i,j,k) * Axx(i,j,k) - TWO * fxx(i,j,k))  + &
           TWO * (  Axx(i,j,k) * betaxx(i,j,k) +   Axy(i,j,k) * betayx(i,j,k) +   Axz(i,j,k) * betazx(i,j,k) )- &
             F2o3 * Axx(i,j,k) * div_beta(i,j,k)

  Ayy_rhs(i,j,k) =           f(i,j,k) * Ayy_rhs(i,j,k)+ alpn1(i,j,k) * (trK(i,j,k) * Ayy(i,j,k) - TWO * fyy(i,j,k))  + &
           TWO * (  Axy(i,j,k) * betaxy(i,j,k) +   Ayy(i,j,k) * betayy(i,j,k) +   Ayz(i,j,k) * betazy(i,j,k) )- &
             F2o3 * Ayy(i,j,k) * div_beta(i,j,k)

  Azz_rhs(i,j,k) =           f(i,j,k) * Azz_rhs(i,j,k)+ alpn1(i,j,k) * (trK(i,j,k) * Azz(i,j,k) - TWO * fzz(i,j,k))  + &
           TWO * (  Axz(i,j,k) * betaxz(i,j,k) +   Ayz(i,j,k) * betayz(i,j,k) +   Azz(i,j,k) * betazz(i,j,k) )- &
             F2o3 * Azz(i,j,k) * div_beta(i,j,k)

  Axy_rhs(i,j,k) =           f(i,j,k) * Axy_rhs(i,j,k)+ alpn1(i,j,k) *( trK(i,j,k) * Axy(i,j,k)  - TWO * fxy(i,j,k) )+ &
                    Axx(i,j,k) * betaxy(i,j,k)                  +   Axz(i,j,k) * betazy(i,j,k)  + &
                                     Ayy(i,j,k) * betayx(i,j,k) +   Ayz(i,j,k) * betazx(i,j,k)  + &
             F1o3 * Axy(i,j,k) * div_beta(i,j,k)                -   Axy(i,j,k) * betazz(i,j,k)

  Ayz_rhs(i,j,k) =           f(i,j,k) * Ayz_rhs(i,j,k)+ alpn1(i,j,k) *( trK(i,j,k) * Ayz(i,j,k)  - TWO * fyz(i,j,k) )+ &
                    Axy(i,j,k) * betaxz(i,j,k) +   Ayy(i,j,k) * betayz(i,j,k)                   + &
                    Axz(i,j,k) * betaxy(i,j,k)                  +   Azz(i,j,k) * betazy(i,j,k)  + &
             F1o3 * Ayz(i,j,k) * div_beta(i,j,k)                -   Ayz(i,j,k) * betaxx(i,j,k)
 
  Axz_rhs(i,j,k) =           f(i,j,k) * Axz_rhs(i,j,k)+ alpn1(i,j,k) *( trK(i,j,k) * Axz(i,j,k)  - TWO * fxz(i,j,k) )+ &
                    Axx(i,j,k) * betaxz(i,j,k) +   Axy(i,j,k) * betayz(i,j,k)                   + &
                                     Ayz(i,j,k) * betayx(i,j,k) +   Azz(i,j,k) * betazx(i,j,k)  + &
             F1o3 * Axz(i,j,k) * div_beta(i,j,k)                -   Axz(i,j,k) * betayy(i,j,k)      !rhs for Aij

! Compute trace of S_ij

  S(i,j,k) =  f(i,j,k) * ( gupxx(i,j,k) * Sxx(i,j,k) + gupyy(i,j,k) * Syy(i,j,k) + gupzz(i,j,k) * Szz(i,j,k) + &
     TWO * ( gupxy(i,j,k) * Sxy(i,j,k) + gupxz(i,j,k) * Sxz(i,j,k) + gupyz(i,j,k) * Syz(i,j,k) ) )

  trK_rhs(i,j,k) = - trK_rhs(i,j,k) + alpn1(i,j,k) *( F1o3 * trK(i,j,k) * trK(i,j,k)         + &
                gupxx(i,j,k) * fxx(i,j,k) + gupyy(i,j,k) * fyy(i,j,k) + gupzz(i,j,k) * fzz(i,j,k)   + &
        TWO * ( gupxy(i,j,k) * fxy(i,j,k) + gupxz(i,j,k) * fxz(i,j,k) + gupyz(i,j,k) * fyz(i,j,k) ) + &
       FOUR * PI * ( rho(i,j,k) + S(i,j,k) ))                                !rhs for trK(i,j,k)
  
!!!! gauge variable part

  Lap_rhs(i,j,k) = -TWO*alpn1(i,j,k)*trK(i,j,k)

  betax_rhs(i,j,k) = FF*dtSfx(i,j,k)
  betay_rhs(i,j,k) = FF*dtSfy(i,j,k)
  betaz_rhs(i,j,k) = FF*dtSfz(i,j,k)

  dtSfx_rhs(i,j,k) = Gamx_rhs(i,j,k) - eta*dtSfx(i,j,k)
  dtSfy_rhs(i,j,k) = Gamy_rhs(i,j,k) - eta*dtSfy(i,j,k)
  dtSfz_rhs(i,j,k) = Gamz_rhs(i,j,k) - eta*dtSfz(i,j,k)
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

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
!g_ij
  call fderivs_shc(ex,dxx,fxx,fxy,fxz,crho,sigma,R,SYM,SYM ,SYM,Symmetry,Lev,sst,              &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  gxx_rhs = gxx_rhs + betax*fxx+betay*fxy+betaz*fxz
  call fderivs_shc(ex,gxy,fxx,fxy,fxz,crho,sigma,R,ANTI,ANTI,SYM,Symmetry,Lev,sst,             &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  gxy_rhs = gxy_rhs + betax*fxx+betay*fxy+betaz*fxz
  call fderivs_shc(ex,gxz,fxx,fxy,fxz,crho,sigma,R,ANTI,SYM ,ANTI,Symmetry,Lev,sst,            &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  gxz_rhs = gxz_rhs + betax*fxx+betay*fxy+betaz*fxz
  call fderivs_shc(ex,dyy,fxx,fxy,fxz,crho,sigma,R,SYM ,SYM ,SYM ,Symmetry,Lev,sst,            &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  gyy_rhs = gyy_rhs + betax*fxx+betay*fxy+betaz*fxz
  call fderivs_shc(ex,gyz,fxx,fxy,fxz,crho,sigma,R,SYM ,ANTI,ANTI,Symmetry,Lev,sst,            &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  gyz_rhs = gyz_rhs + betax*fxx+betay*fxy+betaz*fxz
  call fderivs_shc(ex,dzz,fxx,fxy,fxz,crho,sigma,R,SYM ,SYM ,SYM ,Symmetry,Lev,sst,            &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  gzz_rhs = gzz_rhs + betax*fxx+betay*fxy+betaz*fxz
!A_ij
  call fderivs_shc(ex,Axx,fxx,fxy,fxz,crho,sigma,R,SYM,SYM ,SYM,Symmetry,Lev,sst,              &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  Axx_rhs = Axx_rhs + betax*fxx+betay*fxy+betaz*fxz
  call fderivs_shc(ex,Axy,fxx,fxy,fxz,crho,sigma,R,ANTI,ANTI,SYM,Symmetry,Lev,sst,             &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  Axy_rhs = Axy_rhs + betax*fxx+betay*fxy+betaz*fxz
  call fderivs_shc(ex,Axz,fxx,fxy,fxz,crho,sigma,R,ANTI,SYM ,ANTI,Symmetry,Lev,sst,            &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  Axz_rhs = Axz_rhs + betax*fxx+betay*fxy+betaz*fxz
  call fderivs_shc(ex,Ayy,fxx,fxy,fxz,crho,sigma,R,SYM ,SYM ,SYM ,Symmetry,Lev,sst,            &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  Ayy_rhs = Ayy_rhs + betax*fxx+betay*fxy+betaz*fxz
  call fderivs_shc(ex,Ayz,fxx,fxy,fxz,crho,sigma,R,SYM ,ANTI,ANTI,Symmetry,Lev,sst,            &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  Ayz_rhs = Ayz_rhs + betax*fxx+betay*fxy+betaz*fxz
  call fderivs_shc(ex,Azz,fxx,fxy,fxz,crho,sigma,R,SYM ,SYM ,SYM ,Symmetry,Lev,sst,            &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  Azz_rhs = Azz_rhs + betax*fxx+betay*fxy+betaz*fxz
!chi and trK
  call fderivs_shc(ex,chi,fxx,fxy,fxz,crho,sigma,R,SYM,SYM ,SYM,Symmetry,Lev,sst,              &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  chi_rhs = chi_rhs + betax*fxx+betay*fxy+betaz*fxz
  call fderivs_shc(ex,trK,fxx,fxy,fxz,crho,sigma,R,SYM,SYM ,SYM,Symmetry,Lev,sst,              &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  trK_rhs = trK_rhs + betax*fxx+betay*fxy+betaz*fxz
!Gam^i  
  call fderivs_shc(ex,Gamx,fxx,fxy,fxz,crho,sigma,R,ANTI,SYM ,SYM,Symmetry,Lev,sst,            &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  Gamx_rhs = Gamx_rhs + betax*fxx+betay*fxy+betaz*fxz
  call fderivs_shc(ex,Gamy,fxx,fxy,fxz,crho,sigma,R,SYM ,ANTI,SYM,Symmetry,Lev,sst,            &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  Gamy_rhs = Gamy_rhs + betax*fxx+betay*fxy+betaz*fxz
  call fderivs_shc(ex,Gamz,fxx,fxy,fxz,crho,sigma,R,SYM ,SYM,ANTI,Symmetry,Lev,sst,            &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  Gamz_rhs = Gamz_rhs + betax*fxx+betay*fxy+betaz*fxz
!gauge variables  
  call fderivs_shc(ex,Lap,fxx,fxy,fxz,crho,sigma,R,SYM,SYM ,SYM,Symmetry,Lev,sst,              &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  Lap_rhs = Lap_rhs + betax*fxx+betay*fxy+betaz*fxz
  call fderivs_shc(ex,betax,fxx,fxy,fxz,crho,sigma,R,ANTI,SYM ,SYM,Symmetry,Lev,sst,           &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  betax_rhs = betax_rhs + betax*fxx+betay*fxy+betaz*fxz
  call fderivs_shc(ex,betay,fxx,fxy,fxz,crho,sigma,R,SYM ,ANTI,SYM,Symmetry,Lev,sst,           &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  betay_rhs = betay_rhs + betax*fxx+betay*fxy+betaz*fxz
  call fderivs_shc(ex,betaz,fxx,fxy,fxz,crho,sigma,R,SYM ,SYM,ANTI,Symmetry,Lev,sst,           &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  betaz_rhs = betaz_rhs + betax*fxx+betay*fxy+betaz*fxz
  call fderivs_shc(ex,dtSfx,fxx,fxy,fxz,crho,sigma,R,ANTI,SYM ,SYM,Symmetry,Lev,sst,           &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  dtSfx_rhs = dtSfx_rhs + betax*fxx+betay*fxy+betaz*fxz
  call fderivs_shc(ex,dtSfy,fxx,fxy,fxz,crho,sigma,R,SYM ,ANTI,SYM,Symmetry,Lev,sst,           &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  dtSfy_rhs = dtSfy_rhs + betax*fxx+betay*fxy+betaz*fxz
  call fderivs_shc(ex,dtSfz,fxx,fxy,fxz,crho,sigma,R,SYM ,SYM,ANTI,Symmetry,Lev,sst,           &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  dtSfz_rhs = dtSfz_rhs + betax*fxx+betay*fxy+betaz*fxz
#if 1 
!!
  call fderivs_shc(ex,Sphi,fxx,fxy,fxz,crho,sigma,R,SYM,SYM ,SYM,Symmetry,Lev,sst,              &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  Sphi_rhs = Sphi_rhs + betax*fxx+betay*fxy+betaz*fxz
  call fderivs_shc(ex,Spi,fxx,fxy,fxz,crho,sigma,R,SYM,SYM ,SYM,Symmetry,Lev,sst,              &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  Spi_rhs = Spi_rhs + betax*fxx+betay*fxy+betaz*fxz
#endif
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
  call kodis(ex,X,Y,Z,dtSfx,dtSfx_rhs,ASS,Symmetry,eps)
  call kodis(ex,X,Y,Z,dtSfy,dtSfy_rhs,SAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,dtSfz,dtSfz_rhs,SSA,Symmetry,eps)

#if 1
  call kodis(ex,X,Y,Z,Sphi,Sphi_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Spi ,Spi_rhs ,SSS,Symmetry,eps)
#endif
  endif

#if 0
   Sphi_rhs = ZEO
   Spi_rhs = ZEO
#endif   

  if(co == 0)then
! ham_Res = trR + 2/3 * K^2 - A_ij * A^ij - 16 * PI * rho
! here trR is respect to physical metric
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
  ham_Res(i,j,k) =   gupxx(i,j,k) * Rxx(i,j,k) + gupyy(i,j,k) * Ryy(i,j,k) + gupzz(i,j,k) * Rzz(i,j,k) + &
        TWO* ( gupxy(i,j,k) * Rxy(i,j,k) + gupxz(i,j,k) * Rxz(i,j,k) + gupyz(i,j,k) * Ryz(i,j,k) )

  ham_Res(i,j,k) = chin1(i,j,k)*ham_Res(i,j,k) + F2o3 * trK(i,j,k) * trK(i,j,k) -(&
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
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

! mov_Res_j = gupkj*(-1/chi d_k chi*A_ij + D_k A_ij) - 2/3 d_j trK - 8 PI s_j where D respect to physical metric
! store D_i A_jk - 1/chi d_i chi*A_jk in gjk_i

  call fderivs_shc(ex,Axx,gxxx,gxxy,gxxz,crho,sigma,R, SYM, SYM,SYM,Symmetry,Lev,sst,          &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,Axy,gxyx,gxyy,gxyz,crho,sigma,R,ANTI,ANTI,SYM,Symmetry,Lev,sst,          &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,Axz,gxzx,gxzy,gxzz,crho,sigma,R,ANTI,SYM ,ANTI,Symmetry,Lev,sst,         &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,Ayy,gyyx,gyyy,gyyz,crho,sigma,R, SYM, SYM,SYM,Symmetry,Lev,sst,          &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,Ayz,gyzx,gyzy,gyzz,crho,sigma,R,SYM ,ANTI,ANTI,Symmetry,Lev,sst,         &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,Azz,gzzx,gzzy,gzzz,crho,sigma,R, SYM, SYM,SYM,Symmetry,Lev,sst,          &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
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
movx_Res(i,j,k) = gupxx(i,j,k)*gxxx(i,j,k) + gupyy(i,j,k)*gxyy(i,j,k) + gupzz(i,j,k)*gxzz(i,j,k) &
          +gupxy(i,j,k)*gxyx(i,j,k) + gupxz(i,j,k)*gxzx(i,j,k) + gupyz(i,j,k)*gxzy(i,j,k) &
          +gupxy(i,j,k)*gxxy(i,j,k) + gupxz(i,j,k)*gxxz(i,j,k) + gupyz(i,j,k)*gxyz(i,j,k)
movy_Res(i,j,k) = gupxx(i,j,k)*gxyx(i,j,k) + gupyy(i,j,k)*gyyy(i,j,k) + gupzz(i,j,k)*gyzz(i,j,k) &
          +gupxy(i,j,k)*gyyx(i,j,k) + gupxz(i,j,k)*gyzx(i,j,k) + gupyz(i,j,k)*gyzy(i,j,k) &
          +gupxy(i,j,k)*gxyy(i,j,k) + gupxz(i,j,k)*gxyz(i,j,k) + gupyz(i,j,k)*gyyz(i,j,k)
movz_Res(i,j,k) = gupxx(i,j,k)*gxzx(i,j,k) + gupyy(i,j,k)*gyzy(i,j,k) + gupzz(i,j,k)*gzzz(i,j,k) &
          +gupxy(i,j,k)*gyzx(i,j,k) + gupxz(i,j,k)*gzzx(i,j,k) + gupyz(i,j,k)*gzzy(i,j,k) &
          +gupxy(i,j,k)*gxzy(i,j,k) + gupxz(i,j,k)*gxzz(i,j,k) + gupyz(i,j,k)*gyzz(i,j,k)

movx_Res(i,j,k) = movx_Res(i,j,k) - F2o3*Kx(i,j,k) - F8*PI*sx(i,j,k)
movy_Res(i,j,k) = movy_Res(i,j,k) - F2o3*Ky(i,j,k) - F8*PI*sy(i,j,k)
movz_Res(i,j,k) = movz_Res(i,j,k) - F2o3*Kz(i,j,k) - F8*PI*sz(i,j,k)
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO
  endif

  gont = 0

  return

  end function compute_rhs_bssn_escalar_ss
!-----------------------------------------------------------------------------
!
! compute constraint introduced by dynamical equation reduction of fR
! this routine is valid for both box and shell
!
!-----------------------------------------------------------------------------

  subroutine compute_constraint_fr(ex, X, Y, Z, &
               chi, trK, rho, Sphi,&
               dxx,gxy,gxz,dyy,gyz,dzz, &
               Axx,Axy,Axz,Ayy,Ayz,Azz, &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,&
               Sxx,Sxy,Sxz,Syy,Syz,Szz,&
               Cons_fr)

  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3)
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: dxx,gxy,gxz,dyy,gyz,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: chi,trK,rho,Sphi
! physical Ricci tensor  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: Rxx,Rxy,Rxz,Ryy,Ryz,Rzz
! matter 
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: Sxx,Sxy,Sxz,Syy,Syz,Szz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out):: Cons_fr

! 4D Ricci scalar
  real*8, dimension(ex(1),ex(2),ex(3)) :: RR
  real*8,parameter :: ONE=1.d0,THR=3.d0,FOU=4.d0
  real*8 :: PI

  PI = dacos(-ONE)

  call get4ricciscalar(ex, X, Y, Z, &
               chi, trK, rho, &
               dxx,gxy,gxz,dyy,gyz,dzz, &
               Axx,Axy,Axz,Ayy,Ayz,Azz, &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,&
               Sxx,Sxy,Sxz,Syy,Syz,Szz,&
               RR)
  call frfprim(ex,RR,Cons_fr)

  Cons_fr = Sphi-dsqrt(THR/PI)/FOU*dlog(Cons_fr)

  return

  end subroutine compute_constraint_fr
