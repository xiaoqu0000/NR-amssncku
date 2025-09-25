

#include "macrodef.fh"

  subroutine get_Z4cparameters(kappa1,kappa2,kappa3,FF,eta)

  implicit none

  real*8,intent(out) :: kappa1,kappa2,kappa3,FF,eta

  kappa1 = 2.d-2
  kappa2 = 0.d0
  kappa3 = 0.d0

  FF = 0.75d0
  eta=2.0d0

  return

  end subroutine get_Z4cparameters
#if 1 
! need CPBC_ghost_width
!PRD 83, 024025 (2011)
  subroutine david_milton_cpbc_ss(ex,crho,sigma,R,x,y,z,                       &
               drhodx, drhody, drhodz,                                         &
               dsigmadx,dsigmady,dsigmadz,                                     &
               dRdx,dRdy,dRdz,                                                 &
               drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                &
               dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,    &
               dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,                            &
               xmin,ymin,zmin,xmax,ymax,zmax,                                  &
               TZ,chi,trK,                                                     &
               dxx,gxy,gxz,dyy,gyz,dzz,                                        &
               Axx,Axy,Axz,Ayy,Ayz,Azz,                                        &
               Gamx,Gamy,Gamz,                                                 &
               Lap    ,  betax   ,  betay   ,  betaz   ,                       &
               dtSfx  ,  dtSfy   ,  dtSfz   ,                                  &
               TZ_rhs,chi_rhs,trK_rhs,                                         &
               gxx_rhs,gxy_rhs,gxz_rhs,gyy_rhs,gyz_rhs,gzz_rhs,                &
               Axx_rhs,Axy_rhs,Axz_rhs,Ayy_rhs,Ayz_rhs,Azz_rhs,                &
               Gamx_rhs,Gamy_rhs,Gamz_rhs,                                     &
               Lap_rhs,  betax_rhs,  betay_rhs,  betaz_rhs,                    &
               dtSfx_rhs,  dtSfy_rhs,  dtSfz_rhs,                              &
               pGamxxx,pGamxxy,pGamxxz,pGamxyy,pGamxyz,pGamxzz,                &
               pGamyxx,pGamyxy,pGamyxz,pGamyyy,pGamyyz,pGamyzz,                &
               pGamzxx,pGamzxy,pGamzxz,pGamzyy,pGamzyz,pGamzzz,                &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,                                        &
               Gmxcon,Gmycon,Gmzcon,                                           &
               Symmetry,eps,sst) 

! NOTE: we need Kreiss-Oliger dissipation here instead of rhs calculation routine
  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3), Symmetry,sst
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
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: TZ,chi,dxx,dyy,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: trK
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: gxy,gxz,gyz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamx,Gamy,Gamz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Lap, betax, betay, betaz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: dtSfx,  dtSfy,  dtSfz,Gmxcon,Gmycon,Gmzcon
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: TZ_rhs,chi_rhs,trK_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: gxx_rhs,gxy_rhs,gxz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: gyy_rhs,gyz_rhs,gzz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Axx_rhs,Axy_rhs,Axz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Ayy_rhs,Ayz_rhs,Azz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Gamx_rhs,Gamy_rhs,Gamz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Lap_rhs, betax_rhs, betay_rhs, betaz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: dtSfx_rhs,dtSfy_rhs,dtSfz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Rxx,Rxy,Rxz,Ryy,Ryz,Rzz
  real*8,  intent(in):: xmin,ymin,zmin,xmax,ymax,zmax
  real*8,intent(in) :: eps
! physical second kind of connection  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: pGamxxx, pGamxxy, pGamxxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: pGamxyy, pGamxyz, pGamxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: pGamyxx, pGamyxy, pGamyxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: pGamyyy, pGamyyz, pGamyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: pGamzxx, pGamzxy, pGamzxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: pGamzyy, pGamzyz, pGamzzz

!~~~~~~~~~~~> local variables
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gyy,gzz,alpn1,chin1
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupyy,gupyz,gupzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: qxx,qxy,qxz,qyy,qyz,qzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: qupxx,qupxy,qupxz,qupyy,qupyz,qupzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: qulxx,qulxy,qulxz,qulyx,qulyy,qulyz,qulzx,qulzy,qulzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: slx,sly,slz,ulx,uly,ulz,wlx,wly,wlz
  real*8, dimension(ex(1),ex(2),ex(3)) :: vx,vy,vz,ux,uy,uz,wx,wy,wz
  real*8, dimension(ex(1),ex(2),ex(3)) :: fx,fy,fz
  logical :: gont
  real*8 :: dR
  integer :: i, j, k
  integer :: layer(1:6,1:6),gp
! index of layer, first one: i,j,k; second one: front back etc. boundary
  integer :: kmin,kmax
! derivatives
  real*8 :: Lapx,Lapy,Lapz,Lapxx,Lapxy,Lapxz,Lapyy,Lapyz,Lapzz
  real*8 :: sfxx,sfxy,sfxz,sfyx,sfyy,sfyz,sfzx,sfzy,sfzz
  real*8 :: sfxxx,sfxxy,sfxxz,sfxyy,sfxyz,sfxzz
  real*8 :: sfyxx,sfyxy,sfyxz,sfyyy,sfyyz,sfyzz
  real*8 :: sfzxx,sfzxy,sfzxz,sfzyy,sfzyz,sfzzz
  real*8 :: TZx,TZy,TZz
  real*8 :: chix,chiy,chiz,Kx,Ky,Kz
  real*8 :: chixx,chixy,chixz,chiyy,chiyz,chizz
  real*8 :: Axxx,Axxy,Axxz
  real*8 :: Axyx,Axyy,Axyz
  real*8 :: Axzx,Axzy,Axzz
  real*8 :: Ayyx,Ayyy,Ayyz
  real*8 :: Ayzx,Ayzy,Ayzz
  real*8 :: Azzx,Azzy,Azzz
  real*8 :: gxxx,gxxy,gxxz
  real*8 :: gxyx,gxyy,gxyz
  real*8 :: gxzx,gxzy,gxzz
  real*8 :: gyyx,gyyy,gyyz
  real*8 :: gyzx,gyzy,gyzz
  real*8 :: gzzx,gzzy,gzzz
  real*8 :: gxxxx,gxxxy,gxxxz,gxxyy,gxxyz,gxxzz
  real*8 :: gxyxx,gxyxy,gxyxz,gxyyy,gxyyz,gxyzz
  real*8 :: gxzxx,gxzxy,gxzxz,gxzyy,gxzyz,gxzzz
  real*8 :: gyyxx,gyyxy,gyyxz,gyyyy,gyyyz,gyyzz
  real*8 :: gyzxx,gyzxy,gyzxz,gyzyy,gyzyz,gyzzz
  real*8 :: gzzxx,gzzxy,gzzxz,gzzyy,gzzyz,gzzzz
  real*8 :: Gamxx,Gamxy,Gamxz
  real*8 :: Gamyx,Gamyy,Gamyz
  real*8 :: Gamzx,Gamzy,Gamzz
  
  real*8,dimension(3) ::SSS,AAS,ASA,SAA,ASS,SAS,SSA
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  real*8, parameter :: ZEO = 0.d0, ONE = 1.d0, TWO=2.d0,HALF=0.5d0
  real*8,parameter::TINYRR=1.d-14
! in order to synchronize the following parameters with Z4c_rhs calculation, we
! call a routine
  real*8 :: kappa1,kappa2,kappa3,FF,eta

!  real*8,parameter :: ha=0.d0,thbs=0.d0,hu=0.d0,hw=0.d0,Rhpsi0=0.d0,Ihpsi0=0.d0

  call get_Z4cparameters(kappa1,kappa2,kappa3,FF,eta)

  dR = R(2) - R(1)

  kmax = ex(3)

  kmin = 1

layer(1:3,:) = 1
layer(4:6,:) =-1
 
if(dabs(R(ex(3))-zmax) < dR)then
  layer(1,3) = 1
  layer(2,3) = 1
  layer(3,3) = ex(3) - CPBC_ghost_width
  layer(4,3) = ex(1)
  layer(5,3) = ex(2)
  layer(6,3) = ex(3) - CPBC_ghost_width
endif

if(dabs(R(1)-zmin) < dR)then
   layer(1,6) = 1
   layer(2,6) = 1
   layer(3,6) = 1
   layer(4,6) = ex(1)
   layer(5,6) = ex(2)
   layer(6,6) = 1
endif
! fix BD
  gp = 6

  gont = any( layer(:,gp) == - 1 )
 
   if( .not. gont ) then

    do k = layer(3,gp), layer(6,gp)
     do j = layer(2,gp), layer(5,gp)
      do i = layer(1,gp), layer(4,gp)
! z direction   
        TZ_rhs(i,j,k) = ZEO
        chi_rhs(i,j,k) = ZEO
        trK_rhs(i,j,k) = ZEO
        gxx_rhs(i,j,k) = ZEO
        gxy_rhs(i,j,k) = ZEO
        gxz_rhs(i,j,k) = ZEO
        gyy_rhs(i,j,k) = ZEO
        gyz_rhs(i,j,k) = ZEO
        gzz_rhs(i,j,k) = ZEO
        Axx_rhs(i,j,k) = ZEO
        Axy_rhs(i,j,k) = ZEO
        Axz_rhs(i,j,k) = ZEO
        Ayy_rhs(i,j,k) = ZEO
        Ayz_rhs(i,j,k) = ZEO
        Azz_rhs(i,j,k) = ZEO
        Gamx_rhs(i,j,k) = ZEO
        Gamy_rhs(i,j,k) = ZEO
        Gamz_rhs(i,j,k) = ZEO
        Lap_rhs(i,j,k) = ZEO
        betax_rhs(i,j,k) = ZEO
        betay_rhs(i,j,k) = ZEO
        betaz_rhs(i,j,k) = ZEO
        dtSfx_rhs(i,j,k) = ZEO
        dtSfy_rhs(i,j,k) = ZEO
        dtSfz_rhs(i,j,k) = ZEO
      enddo
     enddo
    enddo
   endif

! constraint preserving BD
  gp = 3

  gont = any( layer(:,gp) == - 1 )
 
  if( .not. gont ) then

! cpbc real starts

  alpn1 = Lap + ONE
  chin1 = chi + ONE
  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE  

    do k = layer(3,gp), layer(6,gp)
     do j = layer(2,gp), layer(5,gp)
      do i = layer(1,gp), layer(4,gp)
!calculate the involved derivatives     
#if 0
         Kx = 0.d0
         Ky = 0.d0
         Kz = 0.d0
         chix = 0.d0
         chiy = 0.d0
         chiz = 0.d0
         Lapx = 0.d0
         Lapy = 0.d0
         Lapz = 0.d0
         TZx = 0.d0
         TZy = 0.d0
         TZz = 0.d0
         Gamxx = 0.d0
         Gamxy = 0.d0
         Gamxz = 0.d0
         Gamyx = 0.d0
         Gamyy = 0.d0
         Gamyz = 0.d0
         Gamzx = 0.d0
         Gamzy = 0.d0
         Gamzz = 0.d0
         sfxx = 0.d0
         sfxy = 0.d0
         sfxz = 0.d0
         sfyx = 0.d0
         sfyy = 0.d0
         sfyz = 0.d0
         sfzx = 0.d0
         sfzy = 0.d0
         sfzz = 0.d0
         Axxx = 0.d0
         Axxy = 0.d0
         Axxz = 0.d0
         Axyx = 0.d0
         Axyy = 0.d0
         Axyz = 0.d0
         Axzx = 0.d0
         Axzy = 0.d0
         Axzz = 0.d0
         Ayyx = 0.d0
         Ayyy = 0.d0
         Ayyz = 0.d0
         Ayzx = 0.d0
         Ayzy = 0.d0
         Ayzz = 0.d0
         Azzx = 0.d0
         Azzy = 0.d0
         Azzz = 0.d0
         gxxx = 0.d0
         gxxy = 0.d0
         gxxz = 0.d0
         gxyx = 0.d0
         gxyy = 0.d0
         gxyz = 0.d0
         gxzx = 0.d0
         gxzy = 0.d0
         gxzz = 0.d0
         gyyx = 0.d0
         gyyy = 0.d0
         gyyz = 0.d0
         gyzx = 0.d0
         gyzy = 0.d0
         gyzz = 0.d0
         gzzx = 0.d0
         gzzy = 0.d0
         gzzz = 0.d0
#else
         call point_fderivs_shc(ex,trK,Kx,Ky,Kz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,chi,chix,chiy,chiz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Lap,Lapx,Lapy,Lapz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,TZ,TZx,TZy,TZz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Gamx,Gamxx,Gamxy,Gamxz,crho,sigma,R,ANTI,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Gamy,Gamyx,Gamyy,Gamyz,crho,sigma,R,SYM,ANTI,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Gamz,Gamzx,Gamzy,Gamzz,crho,sigma,R,SYM,SYM,ANTI,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
#if 0
         sfxx = 0.d0
         sfxy = 0.d0
         sfxz = 0.d0
         sfyx = 0.d0
         sfyy = 0.d0
         sfyz = 0.d0
         sfzx = 0.d0
         sfzy = 0.d0
         sfzz = 0.d0
#else
         call point_fderivs_shc(ex,betax,sfxx,sfxy,sfxz,crho,sigma,R,ANTI,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,betay,sfyx,sfyy,sfyz,crho,sigma,R,SYM,ANTI,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,betaz,sfzx,sfzy,sfzz,crho,sigma,R,SYM,SYM,ANTI,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
#endif                                
         call point_fderivs_shc(ex,Axx,Axxx,Axxy,Axxz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Axy,Axyx,Axyy,Axyz,crho,sigma,R,ANTI,ANTI,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Axz,Axzx,Axzy,Axzz,crho,sigma,R,ANTI,SYM,ANTI,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Ayy,Ayyx,Ayyy,Ayyz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Ayz,Ayzx,Ayzy,Ayzz,crho,sigma,R,SYM,ANTI,ANTI,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Azz,Azzx,Azzy,Azzz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,dxx,gxxx,gxxy,gxxz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,gxy,gxyx,gxyy,gxyz,crho,sigma,R,ANTI,ANTI,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,gxz,gxzx,gxzy,gxzz,crho,sigma,R,ANTI,SYM,ANTI,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,dyy,gyyx,gyyy,gyyz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,gyz,gyzx,gyzy,gyzz,crho,sigma,R,SYM,ANTI,ANTI,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,dzz,gzzx,gzzy,gzzz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
#endif                                

#if 0
         Lapxx = 0.d0
         Lapxy = 0.d0
         Lapxz = 0.d0
         Lapyy = 0.d0
         Lapyz = 0.d0
         Lapzz = 0.d0
         chixx = 0.d0
         chixy = 0.d0
         chixz = 0.d0
         chiyy = 0.d0
         chiyz = 0.d0
         chizz = 0.d0
         gxxxx = 0.d0
         gxxxy = 0.d0
         gxxxz = 0.d0
         gxxyy = 0.d0
         gxxyz = 0.d0
         gxxzz = 0.d0
         gyyxx = 0.d0
         gyyxy = 0.d0
         gyyxz = 0.d0
         gyyyy = 0.d0
         gyyyz = 0.d0
         gyyzz = 0.d0
         gzzxx = 0.d0
         gzzxy = 0.d0
         gzzxz = 0.d0
         gzzyy = 0.d0
         gzzyz = 0.d0
         gzzzz = 0.d0
         gxyxx = 0.d0
         gxyxy = 0.d0
         gxyxz = 0.d0
         gxyyy = 0.d0
         gxyyz = 0.d0
         gxyzz = 0.d0
         gxzxx = 0.d0
         gxzxy = 0.d0
         gxzxz = 0.d0
         gxzyy = 0.d0
         gxzyz = 0.d0
         gxzzz = 0.d0
         gyzxx = 0.d0
         gyzxy = 0.d0
         gyzxz = 0.d0
         gyzyy = 0.d0
         gyzyz = 0.d0
         gyzzz = 0.d0
         sfxxx = 0.d0
         sfxxy = 0.d0
         sfxxz = 0.d0
         sfxyy = 0.d0
         sfxyz = 0.d0
         sfxzz = 0.d0
         sfyxx = 0.d0
         sfyxy = 0.d0
         sfyxz = 0.d0
         sfyyy = 0.d0
         sfyyz = 0.d0
         sfyzz = 0.d0
         sfzxx = 0.d0
         sfzxy = 0.d0
         sfzxz = 0.d0
         sfzyy = 0.d0
         sfzyz = 0.d0
         sfzzz = 0.d0
#else         
         call point_fdderivs_shc(ex,Lap,Lapxx,Lapxy,Lapxz,Lapyy,Lapyz,Lapzz,crho,sigma,R,SYM ,SYM ,SYM ,Symmetry,0,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)
         call point_fdderivs_shc(ex,chi,chixx,chixy,chixz,chiyy,chiyz,chizz,crho,sigma,R,SYM ,SYM ,SYM ,Symmetry,0,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)
         call point_fdderivs_shc(ex,dxx,gxxxx,gxxxy,gxxxz,gxxyy,gxxyz,gxxzz,crho,sigma,R,SYM ,SYM ,SYM ,Symmetry,0,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)
         call point_fdderivs_shc(ex,dyy,gyyxx,gyyxy,gyyxz,gyyyy,gyyyz,gyyzz,crho,sigma,R,SYM ,SYM ,SYM ,Symmetry,0,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)
         call point_fdderivs_shc(ex,dzz,gzzxx,gzzxy,gzzxz,gzzyy,gzzyz,gzzzz,crho,sigma,R,SYM ,SYM ,SYM ,Symmetry,0,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)
         call point_fdderivs_shc(ex,gxy,gxyxx,gxyxy,gxyxz,gxyyy,gxyyz,gxyzz,crho,sigma,R,ANTI,ANTI,SYM ,Symmetry,0,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)
         call point_fdderivs_shc(ex,gxz,gxzxx,gxzxy,gxzxz,gxzyy,gxzyz,gxzzz,crho,sigma,R,ANTI,SYM ,ANTI,Symmetry,0,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)
         call point_fdderivs_shc(ex,gyz,gyzxx,gyzxy,gyzxz,gyzyy,gyzyz,gyzzz,crho,sigma,R,SYM ,ANTI,ANTI ,Symmetry,0,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)
         call point_fdderivs_shc(ex,betax,sfxxx,sfxxy,sfxxz,sfxyy,sfxyz,sfxzz,crho,sigma,R,ANTI,SYM ,SYM ,Symmetry,0,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)
         call point_fdderivs_shc(ex,betay,sfyxx,sfyxy,sfyxz,sfyyy,sfyyz,sfyzz,crho,sigma,R,SYM ,ANTI,SYM ,Symmetry,0,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)
         call point_fdderivs_shc(ex,betaz,sfzxx,sfzxy,sfzxz,sfzyy,sfzyz,sfzzz,crho,sigma,R,SYM ,SYM ,ANTI,Symmetry,0,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)
#endif

        call cpbc_point(R(k),x(i,j,k),y(i,j,k),z(i,j,k),TZ(i,j,k),chin1(i,j,k),trK(i,j,k), &
                         gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k), &
                         Axx(i,j,k),Axy(i,j,k),Axz(i,j,k),Ayy(i,j,k),Ayz(i,j,k),Azz(i,j,k), &
                         Gamx(i,j,k),Gamy(i,j,k),Gamz(i,j,k),                               &
                         alpn1(i,j,k),betax(i,j,k),betay(i,j,k),betaz(i,j,k),               &
                         Lapx,Lapy,Lapz,Lapxx,Lapxy,Lapxz,Lapyy,Lapyz,Lapzz,                & 
                         sfxx,sfxy,sfxz,                                                    &
                         sfyx,sfyy,sfyz,                                                    &
                         sfzx,sfzy,sfzz,                                                    &
                         sfxxx,sfxxy,sfxxz,sfxyy,sfxyz,sfxzz,                               &
                         sfyxx,sfyxy,sfyxz,sfyyy,sfyyz,sfyzz,                               &
                         sfzxx,sfzxy,sfzxz,sfzyy,sfzyz,sfzzz,                               &
                         chix,chiy,chiz,chixx,chixy,chixz,chiyy,chiyz,chizz,                &
                         gxxx,gxyx,gxzx,gyyx,gyzx,gzzx,                                     &
                         gxxy,gxyy,gxzy,gyyy,gyzy,gzzy,                                     &
                         gxxz,gxyz,gxzz,gyyz,gyzz,gzzz,                                     &
                         gxxxx,gxxxy,gxxxz,gxxyy,gxxyz,gxxzz,                               &
                         gxyxx,gxyxy,gxyxz,gxyyy,gxyyz,gxyzz,                               &
                         gxzxx,gxzxy,gxzxz,gxzyy,gxzyz,gxzzz,                               &
                         gyyxx,gyyxy,gyyxz,gyyyy,gyyyz,gyyzz,                               &
                         gyzxx,gyzxy,gyzxz,gyzyy,gyzyz,gyzzz,                               &
                         gzzxx,gzzxy,gzzxz,gzzyy,gzzyz,gzzzz,                               &
                         Kx,Ky,Kz,                                                          &
                         Axxx,Axyx,Axzx,Ayyx,Ayzx,Azzx,                                     &
                         Axxy,Axyy,Axzy,Ayyy,Ayzy,Azzy,                                     &
                         Axxz,Axyz,Axzz,Ayyz,Ayzz,Azzz,                                     &
                         Gamxx,Gamxy,Gamxz,                                                 &
                         Gamyx,Gamyy,Gamyz,                                                 &
                         Gamzx,Gamzy,Gamzz,                                                 &
                         TZx,TZy,TZz,                                                       &
                         trK_rhs(i,j,k),TZ_rhs(i,j,k),                                      &
                         Axx_rhs(i,j,k),Axy_rhs(i,j,k),Axz_rhs(i,j,k),Ayy_rhs(i,j,k),Ayz_rhs(i,j,k),Azz_rhs(i,j,k), &
                         Gamx_rhs(i,j,k),Gamy_rhs(i,j,k),Gamz_rhs(i,j,k),kappa1,kappa2,eta)  
      enddo
     enddo
    enddo

  endif

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

! NOTE: we need Kreiss-Oliger dissipation here instead of rhs calculation routine
  if(eps>0)then 
! usual Kreiss-Oliger dissipation      
  call kodis_sh(ex,crho,sigma,R,chi,chi_rhs,SSS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,trK,trK_rhs,SSS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,dxx,gxx_rhs,SSS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,gxy,gxy_rhs,AAS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,gxz,gxz_rhs,ASA,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,dyy,gyy_rhs,SSS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,gyz,gyz_rhs,SAA,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,dzz,gzz_rhs,SSS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,Axx,Axx_rhs,SSS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,Axy,Axy_rhs,AAS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,Axz,Axz_rhs,ASA,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,Ayy,Ayy_rhs,SSS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,Ayz,Ayz_rhs,SAA,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,Azz,Azz_rhs,SSS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,Gamx,Gamx_rhs,ASS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,Gamy,Gamy_rhs,SAS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,Gamz,Gamz_rhs,SSA,Symmetry,eps,sst)

  call kodis_sh(ex,crho,sigma,R,Lap,Lap_rhs,SSS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,betax,betax_rhs,ASS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,betay,betay_rhs,SAS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,betaz,betaz_rhs,SSA,Symmetry,eps,sst)

#if 0
  call kodis_sh(ex,crho,sigma,R,dtSfx,dtSfx_rhs,ASS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,dtSfy,dtSfy_rhs,SAS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,dtSfz,dtSfz_rhs,SSA,Symmetry,eps,sst)
#endif

  call kodis_sh(ex,crho,sigma,R,TZ,TZ_rhs,SSS,Symmetry,eps,sst)
  endif

  return

  end subroutine david_milton_cpbc_ss
#elif 1
#error "did you change sommerfeld routine for buffer points considering?"
!!! CV == 0: Sommerfeld on everything after decomposing
!!! CV == 1: Sommerfeld on only the CPBC vars after decomposing
!!! CV == 1 and replace Sommerfeld to CPBC one by one
#define CV 1
! Sommefeld after 2+1 decomposation
  subroutine david_milton_cpbc_ss(ex,crho,sigma,R,x,y,z,                       &
               drhodx, drhody, drhodz,                                         &
               dsigmadx,dsigmady,dsigmadz,                                     &
               dRdx,dRdy,dRdz,                                                 &
               drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                &
               dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,    &
               dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,                            &
               xmin,ymin,zmin,xmax,ymax,zmax,                                  &
               TZ,chi,trK,                                                     &
               dxx,gxy,gxz,dyy,gyz,dzz,                                        &
               Axx,Axy,Axz,Ayy,Ayz,Azz,                                        &
               Gamx,Gamy,Gamz,                                                 &
               Lap    ,  betax   ,  betay   ,  betaz   ,                       &
               dtSfx  ,  dtSfy   ,  dtSfz   ,                                  &
               TZ_rhs,chi_rhs,trK_rhs,                                         &
               gxx_rhs,gxy_rhs,gxz_rhs,gyy_rhs,gyz_rhs,gzz_rhs,                &
               Axx_rhs,Axy_rhs,Axz_rhs,Ayy_rhs,Ayz_rhs,Azz_rhs,                &
               Gamx_rhs,Gamy_rhs,Gamz_rhs,                                     &
               Lap_rhs,  betax_rhs,  betay_rhs,  betaz_rhs,                    &
               dtSfx_rhs,  dtSfy_rhs,  dtSfz_rhs,                              &
               pGamxxx,pGamxxy,pGamxxz,pGamxyy,pGamxyz,pGamxzz,                &
               pGamyxx,pGamyxy,pGamyxz,pGamyyy,pGamyyz,pGamyzz,                &
               pGamzxx,pGamzxy,pGamzxz,pGamzyy,pGamzyz,pGamzzz,                &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,                                        &
               Gmxcon,Gmycon,Gmzcon,                                           &
               Symmetry,eps,sst) 

! NOTE: we need Kreiss-Oliger dissipation here instead of rhs calculation routine
  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3), Symmetry,sst
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
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: TZ,chi,dxx,dyy,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: trK
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: gxy,gxz,gyz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamx,Gamy,Gamz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Lap, betax, betay, betaz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: dtSfx,  dtSfy,  dtSfz,Gmxcon,Gmycon,Gmzcon
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: TZ_rhs,chi_rhs,trK_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: gxx_rhs,gxy_rhs,gxz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: gyy_rhs,gyz_rhs,gzz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Axx_rhs,Axy_rhs,Axz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Ayy_rhs,Ayz_rhs,Azz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Gamx_rhs,Gamy_rhs,Gamz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Lap_rhs, betax_rhs, betay_rhs, betaz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: dtSfx_rhs,dtSfy_rhs,dtSfz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Rxx,Rxy,Rxz,Ryy,Ryz,Rzz
  real*8,  intent(in):: xmin,ymin,zmin,xmax,ymax,zmax
  real*8,intent(in) :: eps
! physical second kind of connection  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: pGamxxx, pGamxxy, pGamxxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: pGamxyy, pGamxyz, pGamxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: pGamyxx, pGamyxy, pGamyxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: pGamyyy, pGamyyz, pGamyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: pGamzxx, pGamzxy, pGamzxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: pGamzyy, pGamzyz, pGamzzz

!~~~~~~~~~~~> local variables

  real*8, dimension(ex(1),ex(2),ex(3)) :: chin1,gxx,gyy,gzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: toAqq,toAss,toAsx,toAsy,toAsz
  real*8, dimension(ex(1),ex(2),ex(3)) :: toAxx,toAxy,toAxz,toAyy,toAyz,toAzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: toAqq_rhs,toAss_rhs,toAsx_rhs,toAsy_rhs,toAsz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)) :: toAxx_rhs,toAxy_rhs,toAxz_rhs,toAyy_rhs,toAyz_rhs,toAzz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)) :: toGams,toGamx,toGamy,toGamz
  real*8, dimension(ex(1),ex(2),ex(3)) :: toGams_rhs,toGamx_rhs,toGamy_rhs,toGamz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)) :: tobetas,tobetax,tobetay,tobetaz
  real*8, dimension(ex(1),ex(2),ex(3)) :: tobetas_rhs,tobetax_rhs,tobetay_rhs,tobetaz_rhs

  real*8,dimension(3) ::SSS,AAS,ASA,SAA,ASS,SAS,SSA
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  real*8, parameter :: ZEO = 0.d0

  logical :: gont
  real*8 :: dR
  integer :: i, j, k
  integer :: layer(1:6,1:6),gp

! in order to synchronize the following parameters with Z4c_rhs calculation, we
! call a routine
  real*8 :: kappa1,kappa2,kappa3,FF,eta

!  real*8,parameter :: ha=0.d0,thbs=0.d0,hu=0.d0,hw=0.d0,Rhpsi0=0.d0,Ihpsi0=0.d0

  call get_Z4cparameters(kappa1,kappa2,kappa3,FF,eta)

  dR = R(2) - R(1)

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

#if 1
  chin1 = chi+1.d0
  gxx = dxx+1.d0
  gyy = dyy+1.d0
  gzz = dzz+1.d0

! decompose
    do k = 1, ex(3)
     do j = 1, ex(2)
      do i = 1, ex(1)
#if (CV == 0)  
      call decompose2p1_1(R(k),x(i,j,k),y(i,j,k),z(i,j,k), chin1(i,j,k), &
                        gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k), &
                        betax(i,j,k),betay(i,j,k),betaz(i,j,k), &
                        tobetas(i,j,k),tobetax(i,j,k),tobetay(i,j,k),tobetaz(i,j,k))
#endif                        
      call decompose2p1_1(R(k),x(i,j,k),y(i,j,k),z(i,j,k), chin1(i,j,k), &
                        gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k), &
                        Gamx(i,j,k),Gamy(i,j,k),Gamz(i,j,k), &
                        toGams(i,j,k),toGamx(i,j,k),toGamy(i,j,k),toGamz(i,j,k))
      call decompose2p1_2(R(k),x(i,j,k),y(i,j,k),z(i,j,k), chin1(i,j,k), &
                        gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k), &
                        Axx(i,j,k),Axy(i,j,k),Axz(i,j,k),Ayy(i,j,k),Ayz(i,j,k),Azz(i,j,k), &
                        toAqq(i,j,k),toAss(i,j,k),toAsx(i,j,k),toAsy(i,j,k),toAsz(i,j,k), &
                        toAxx(i,j,k),toAxy(i,j,k),toAxz(i,j,k),toAyy(i,j,k),toAyz(i,j,k),toAzz(i,j,k))

      enddo
     enddo
    enddo

! sommerfeld boundary
! cpbc variables       
#if 0
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,trK_rhs,trK,1.d0,SSS,Symmetry)    
#else
  call cpbcrtrK(ex,crho,sigma,R,x,y,z,                                  &
               drhodx, drhody, drhodz,                                         &
               dsigmadx,dsigmady,dsigmadz,                                     &
               dRdx,dRdy,dRdz,                                                 &
               drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                &
               dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,    &
               dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,                            &
               xmin,ymin,zmin,xmax,ymax,zmax,trK_rhs,&
               chi,trK,dxx,gxy,gxz,dyy,gyz,dzz, &
               Lap,betax,betay,betaz,TZ,Symmetry,sst,kappa1,kappa2)  
#endif
#if 0
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,TZ_rhs,TZ,1.d0,SSS,Symmetry)  
#else
  call cpbcrtheta(ex,crho,sigma,R,x,y,z,                                  &
               drhodx, drhody, drhodz,                                         &
               dsigmadx,dsigmady,dsigmadz,                                     &
               dRdx,dRdy,dRdz,                                                 &
               drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                &
               dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,    &
               dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,                            &
               xmin,ymin,zmin,xmax,ymax,zmax,TZ_rhs,&
               chi,dxx,gxy,gxz,dyy,gyz,dzz, &
               Lap,betax,betay,betaz,TZ,Symmetry,sst,kappa1,kappa2) 
#endif
#if 1
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,toGams_rhs,toGams,1.d0,SSS,Symmetry)    
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,toGamx_rhs,toGamx,1.d0,ASS,Symmetry)         
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,toGamy_rhs,toGamy,1.d0,SAS,Symmetry)         
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,toGamz_rhs,toGamz,1.d0,SSA,Symmetry)
#else  
  call cpbcrgam(ex,crho,sigma,R,x,y,z,                                  &
               drhodx, drhody, drhodz,                                         &
               dsigmadx,dsigmady,dsigmadz,                                     &
               dRdx,dRdy,dRdz,                                                 &
               drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                &
               dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,    &
               dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,                            &
               xmin,ymin,zmin,xmax,ymax,zmax,toGamx_rhs,toGamy_rhs,toGamz_rhs,toGams_rhs,&
               chi,trK,dxx,gxy,gxz,dyy,gyz,dzz, &
               Lap,betax,betay,betaz,TZ,Gamx,Gamy,Gamz,Symmetry,sst,eta) 
#endif
#if 1
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,toAss_rhs,toAss,1.d0,SSS,Symmetry)    
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,toAsx_rhs,toAsx,1.d0,ASS,Symmetry)    
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,toAsy_rhs,toAsy,1.d0,SAS,Symmetry)    
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,toAsz_rhs,toAsz,1.d0,SSA,Symmetry)    
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,toAxx_rhs,toAxx,1.d0,SSS,Symmetry)    
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,toAxy_rhs,toAxy,1.d0,AAS,Symmetry)    
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,toAxz_rhs,toAxz,1.d0,ASA,Symmetry)    
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,toAyy_rhs,toAyy,1.d0,SSS,Symmetry)    
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,toAyz_rhs,toAyz,1.d0,SAA,Symmetry)   
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,toAzz_rhs,toAzz,1.d0,SSS,Symmetry)  
#else
  call cpbcra(ex,crho,sigma,R,x,y,z,                                     &
               drhodx, drhody, drhodz,                                         &
               dsigmadx,dsigmady,dsigmadz,                                     &
               dRdx,dRdy,dRdz,                                                 &
               drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                &
               dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,    &
               dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,                            &
               xmin,ymin,zmin,xmax,ymax,zmax, &
               toAxx_rhs,toAxy_rhs,toAxz_rhs,toAyy_rhs,toAyz_rhs,toAzz_rhs,&
               toAsx_rhs,toAsy_rhs,toAsz_rhs,toAss_rhs, &
               chi,trK,dxx,gxy,gxz,dyy,gyz,dzz, &
               Axx,Axy,Axz,Ayy,Ayz,Azz, &
               Lap,betax,betay,betaz,TZ,Gamx,Gamy,Gamz,Symmetry,sst,kappa1)
#endif
! non-cpbc variables
#if (CV == 0)
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,toAqq_rhs,toAqq,1.d0,SSS,Symmetry) 
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,chi_rhs,chi,1.d0,SSS,Symmetry)        
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,gxx_rhs,dxx,1.d0,SSS,Symmetry)    
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,gxy_rhs,gxy,1.d0,AAS,Symmetry)    
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,gxz_rhs,gxz,1.d0,ASA,Symmetry)    
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,gyy_rhs,dyy,1.d0,SSS,Symmetry)    
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,gyz_rhs,gyz,1.d0,SAA,Symmetry)   
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,gzz_rhs,dzz,1.d0,SSS,Symmetry)           
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,Lap_rhs,Lap,1.d0,SSS,Symmetry)    
#if 1  
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,tobetas_rhs,tobetas,1.d0,SSS,Symmetry)    
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,tobetax_rhs,tobetax,1.d0,ASS,Symmetry)         
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,tobetay_rhs,tobetay,1.d0,SAS,Symmetry)         
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,tobetaz_rhs,tobetaz,1.d0,SSA,Symmetry)     
#else  
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,betax_rhs,betax,1.d0,ASS,Symmetry)          
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,betay_rhs,betay,1.d0,SAS,Symmetry)          
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,betaz_rhs,betaz,1.d0,SSA,Symmetry)     
#endif  
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,dtSfx_rhs,dtSfx,1.d0,ASS,Symmetry)          
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,dtSfy_rhs,dtSfy,1.d0,SAS,Symmetry)          
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,dtSfz_rhs,dtSfz,1.d0,SSA,Symmetry)   

#else  
  call cpbcrACqq(ex,crho,sigma,R,x,y,z,                                  &
               drhodx, drhody, drhodz,                                         &
               dsigmadx,dsigmady,dsigmadz,                                     &
               dRdx,dRdy,dRdz,                                                 &
               drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                &
               dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,    &
               dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,                            &
               xmin,ymin,zmin,xmax,ymax,zmax,toAqq_rhs,&
               chi,dxx,gxy,gxz,dyy,gyz,dzz, &
               Lap,betax,betay,betaz,Axx,Axy,Axz,Ayy,Ayz,Azz,toAss_rhs,Symmetry,sst)
#endif
! reconstruct
layer(1:3,:) = 1
layer(4:6,:) =-1
 
if(dabs(R(ex(3))-zmax) < dR)then
  layer(1,3) = 1
  layer(2,3) = 1
  layer(3,3) = ex(3) - CPBC_ghost_width
  layer(4,3) = ex(1)
  layer(5,3) = ex(2)
  layer(6,3) = ex(3) - CPBC_ghost_width
endif

if(dabs(R(1)-zmin) < dR)then
   layer(1,6) = 1
   layer(2,6) = 1
   layer(3,6) = 1
   layer(4,6) = ex(1)
   layer(5,6) = ex(2)
   layer(6,6) = 1
endif
! fix BD
  gp = 6

  gont = any( layer(:,gp) == - 1 )
 
   if( .not. gont ) then

    do k = layer(3,gp), layer(6,gp)
     do j = layer(2,gp), layer(5,gp)
      do i = layer(1,gp), layer(4,gp)
! z direction   
        TZ_rhs(i,j,k) = ZEO
        chi_rhs(i,j,k) = ZEO
        trK_rhs(i,j,k) = ZEO
        gxx_rhs(i,j,k) = ZEO
        gxy_rhs(i,j,k) = ZEO
        gxz_rhs(i,j,k) = ZEO
        gyy_rhs(i,j,k) = ZEO
        gyz_rhs(i,j,k) = ZEO
        gzz_rhs(i,j,k) = ZEO
        Axx_rhs(i,j,k) = ZEO
        Axy_rhs(i,j,k) = ZEO
        Axz_rhs(i,j,k) = ZEO
        Ayy_rhs(i,j,k) = ZEO
        Ayz_rhs(i,j,k) = ZEO
        Azz_rhs(i,j,k) = ZEO
        Gamx_rhs(i,j,k) = ZEO
        Gamy_rhs(i,j,k) = ZEO
        Gamz_rhs(i,j,k) = ZEO
        Lap_rhs(i,j,k) = ZEO
        betax_rhs(i,j,k) = ZEO
        betay_rhs(i,j,k) = ZEO
        betaz_rhs(i,j,k) = ZEO
        dtSfx_rhs(i,j,k) = ZEO
        dtSfy_rhs(i,j,k) = ZEO
        dtSfz_rhs(i,j,k) = ZEO
      enddo
     enddo
    enddo
   endif

! constraint preserving BD
  gp = 3

  gont = any( layer(:,gp) == - 1 )
 
  if( .not. gont ) then

    do k = layer(3,gp), layer(6,gp)
     do j = layer(2,gp), layer(5,gp)
      do i = layer(1,gp), layer(4,gp)
#if (CV == 0)  
      call compose2p1_1(R(k),x(i,j,k),y(i,j,k),z(i,j,k), chin1(i,j,k), &
                        gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k), &
                        betax_rhs(i,j,k),betay_rhs(i,j,k),betaz_rhs(i,j,k), &
                        tobetas_rhs(i,j,k),tobetax_rhs(i,j,k),tobetay_rhs(i,j,k),tobetaz_rhs(i,j,k))
#endif                        
      call compose2p1_1(R(k),x(i,j,k),y(i,j,k),z(i,j,k), chin1(i,j,k), &
                        gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k), &
                        Gamx_rhs(i,j,k),Gamy_rhs(i,j,k),Gamz_rhs(i,j,k), &
                        toGams_rhs(i,j,k),toGamx_rhs(i,j,k),toGamy_rhs(i,j,k),toGamz_rhs(i,j,k))
      call compose2p1_2(R(k),x(i,j,k),y(i,j,k),z(i,j,k), chin1(i,j,k), &
                        gxx(i,j,k),gxy(i,j,k),gxz(i,j,k),gyy(i,j,k),gyz(i,j,k),gzz(i,j,k), &
                        Axx_rhs(i,j,k),Axy_rhs(i,j,k),Axz_rhs(i,j,k),Ayy_rhs(i,j,k),Ayz_rhs(i,j,k),Azz_rhs(i,j,k), &
                        toAqq_rhs(i,j,k),toAss_rhs(i,j,k),toAsx_rhs(i,j,k),toAsy_rhs(i,j,k),toAsz_rhs(i,j,k), &
                        toAxx_rhs(i,j,k),toAxy_rhs(i,j,k),toAxz_rhs(i,j,k),toAyy_rhs(i,j,k),toAyz_rhs(i,j,k),toAzz_rhs(i,j,k))

      enddo
     enddo
    enddo

  endif

! check direct Sommerfeld BD
#else
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,trK_rhs,trK,1.d0,SSS,Symmetry)       
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,TZ_rhs,TZ,1.d0,SSS,Symmetry)      
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,Gamx_rhs,Gamx,1.d0,ASS,Symmetry)         
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,Gamy_rhs,Gamy,1.d0,SAS,Symmetry)         
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,Gamz_rhs,Gamz,1.d0,SSA,Symmetry)   
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,Axx_rhs,Axx,1.d0,SSS,Symmetry)    
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,Axy_rhs,Axy,1.d0,AAS,Symmetry)    
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,Axz_rhs,Axz,1.d0,ASA,Symmetry)    
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,Ayy_rhs,Ayy,1.d0,SSS,Symmetry)    
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,Ayz_rhs,Ayz,1.d0,SAA,Symmetry)   
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,Azz_rhs,Azz,1.d0,SSS,Symmetry)  
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,chi_rhs,chi,1.d0,SSS,Symmetry)        
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,gxx_rhs,dxx,1.d0,SSS,Symmetry)    
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,gxy_rhs,gxy,1.d0,AAS,Symmetry)    
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,gxz_rhs,gxz,1.d0,ASA,Symmetry)    
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,gyy_rhs,dyy,1.d0,SSS,Symmetry)    
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,gyz_rhs,gyz,1.d0,SAA,Symmetry)   
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,gzz_rhs,dzz,1.d0,SSS,Symmetry)           
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,Lap_rhs,Lap,1.d0,SSS,Symmetry)   
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,betax_rhs,betax,1.d0,ASS,Symmetry)          
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,betay_rhs,betay,1.d0,SAS,Symmetry)          
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,betaz_rhs,betaz,1.d0,SSA,Symmetry)    
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,dtSfx_rhs,dtSfx,1.d0,ASS,Symmetry)          
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,dtSfy_rhs,dtSfy,1.d0,SAS,Symmetry)          
  call sommerfeld_routbam_ss(ex,crho,sigma,R,xmin,ymin,zmin,xmax,ymax,zmax,dtSfz_rhs,dtSfz,1.d0,SSA,Symmetry)  
#endif  

! NOTE: we need Kreiss-Oliger dissipation here instead of rhs calculation routine
  if(eps>0)then 
! usual Kreiss-Oliger dissipation      
  call kodis_sh(ex,crho,sigma,R,chi,chi_rhs,SSS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,trK,trK_rhs,SSS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,dxx,gxx_rhs,SSS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,gxy,gxy_rhs,AAS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,gxz,gxz_rhs,ASA,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,dyy,gyy_rhs,SSS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,gyz,gyz_rhs,SAA,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,dzz,gzz_rhs,SSS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,Axx,Axx_rhs,SSS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,Axy,Axy_rhs,AAS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,Axz,Axz_rhs,ASA,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,Ayy,Ayy_rhs,SSS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,Ayz,Ayz_rhs,SAA,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,Azz,Azz_rhs,SSS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,Gamx,Gamx_rhs,ASS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,Gamy,Gamy_rhs,SAS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,Gamz,Gamz_rhs,SSA,Symmetry,eps,sst)

  call kodis_sh(ex,crho,sigma,R,Lap,Lap_rhs,SSS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,betax,betax_rhs,ASS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,betay,betay_rhs,SAS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,betaz,betaz_rhs,SSA,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,dtSfx,dtSfx_rhs,ASS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,dtSfy,dtSfy_rhs,SAS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,dtSfz,dtSfz_rhs,SSA,Symmetry,eps,sst)

  call kodis_sh(ex,crho,sigma,R,TZ,TZ_rhs,SSS,Symmetry,eps,sst)
  endif

  return

  end subroutine david_milton_cpbc_ss
#undef CV  
#else
!out of time code, never debuged
! need CPBC_ghost_width
!PRD 83, 024025 (2011)
  subroutine david_milton_cpbc_ss(ex,crho,sigma,R,x,y,z,                       &
               drhodx, drhody, drhodz,                                         &
               dsigmadx,dsigmady,dsigmadz,                                     &
               dRdx,dRdy,dRdz,                                                 &
               drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                &
               dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,    &
               dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,                            &
               xmin,ymin,zmin,xmax,ymax,zmax,                                  &
               TZ,chi,trK,                                                     &
               dxx,gxy,gxz,dyy,gyz,dzz,                                        &
               Axx,Axy,Axz,Ayy,Ayz,Azz,                                        &
               Gamx,Gamy,Gamz,                                                 &
               Lap    ,  betax   ,  betay   ,  betaz   ,                       &
               dtSfx  ,  dtSfy   ,  dtSfz   ,                                  &
               TZ_rhs,chi_rhs,trK_rhs,                                         &
               gxx_rhs,gxy_rhs,gxz_rhs,gyy_rhs,gyz_rhs,gzz_rhs,                &
               Axx_rhs,Axy_rhs,Axz_rhs,Ayy_rhs,Ayz_rhs,Azz_rhs,                &
               Gamx_rhs,Gamy_rhs,Gamz_rhs,                                     &
               Lap_rhs,  betax_rhs,  betay_rhs,  betaz_rhs,                    &
               dtSfx_rhs,  dtSfy_rhs,  dtSfz_rhs,                              &
               pGamxxx,pGamxxy,pGamxxz,pGamxyy,pGamxyz,pGamxzz,                &
               pGamyxx,pGamyxy,pGamyxz,pGamyyy,pGamyyz,pGamyzz,                &
               pGamzxx,pGamzxy,pGamzxz,pGamzyy,pGamzyz,pGamzzz,                &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,                                        &
               Gmxcon,Gmycon,Gmzcon,                                           &
               Symmetry,eps,sst) 

! NOTE: we need Kreiss-Oliger dissipation here instead of rhs calculation routine
  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3), Symmetry,sst
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
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: TZ,chi,dxx,dyy,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: trK
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: gxy,gxz,gyz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamx,Gamy,Gamz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Lap, betax, betay, betaz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: dtSfx,  dtSfy,  dtSfz,Gmxcon,Gmycon,Gmzcon
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: TZ_rhs,chi_rhs,trK_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: gxx_rhs,gxy_rhs,gxz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: gyy_rhs,gyz_rhs,gzz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Axx_rhs,Axy_rhs,Axz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Ayy_rhs,Ayz_rhs,Azz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Gamx_rhs,Gamy_rhs,Gamz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Lap_rhs, betax_rhs, betay_rhs, betaz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: dtSfx_rhs,dtSfy_rhs,dtSfz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Rxx,Rxy,Rxz,Ryy,Ryz,Rzz
  real*8,  intent(in):: xmin,ymin,zmin,xmax,ymax,zmax
  real*8,intent(in) :: eps
! physical second kind of connection  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: pGamxxx, pGamxxy, pGamxxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: pGamxyy, pGamxyz, pGamxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: pGamyxx, pGamyxy, pGamyxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: pGamyyy, pGamyyz, pGamyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: pGamzxx, pGamzxy, pGamzxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: pGamzyy, pGamzyz, pGamzzz

!~~~~~~~~~~~> local variables
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gyy,gzz,alpn1,chin1
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupyy,gupyz,gupzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: qxx,qxy,qxz,qyy,qyz,qzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: qupxx,qupxy,qupxz,qupyy,qupyz,qupzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: qulxx,qulxy,qulxz,qulyx,qulyy,qulyz,qulzx,qulzy,qulzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: slx,sly,slz,ulx,uly,ulz,wlx,wly,wlz
  real*8, dimension(ex(1),ex(2),ex(3)) :: vx,vy,vz,ux,uy,uz,wx,wy,wz
  real*8, dimension(ex(1),ex(2),ex(3)) :: fx,fy,fz
  logical :: gont
  real*8 :: dR
  integer :: i, j, k
  integer :: layer(1:6,1:6),gp
! index of layer, first one: i,j,k; second one: front back etc. boundary
  integer :: kmin,kmax
  real*8 :: toAss_rhs,toAqq_rhs,toAs1_rhs,toAs2_rhs,toA11_rhs,toA12_rhs,toA22_rhs
  real*8 :: toGams_rhs,toGam1_rhs,toGam2_rhs
  real*8 :: totrK_rhs,toTZ_rhs
! derivatives
  real*8 :: Lapx,Lapy,Lapz,Lapxx,Lapxy,Lapxz,Lapyy,Lapyz,Lapzz
  real*8 :: sfxx,sfxy,sfxz,sfyx,sfyy,sfyz,sfzx,sfzy,sfzz
  real*8 :: sfxxx,sfxxy,sfxxz,sfxyy,sfxyz,sfxzz
  real*8 :: sfyxx,sfyxy,sfyxz,sfyyy,sfyyz,sfyzz
  real*8 :: sfzxx,sfzxy,sfzxz,sfzyy,sfzyz,sfzzz
  real*8 :: TZx,TZy,TZz
  real*8 :: chix,chiy,chiz,Kx,Ky,Kz
  real*8 :: Axxx,Axxy,Axxz
  real*8 :: Axyx,Axyy,Axyz
  real*8 :: Axzx,Axzy,Axzz
  real*8 :: Ayyx,Ayyy,Ayyz
  real*8 :: Ayzx,Ayzy,Ayzz
  real*8 :: Azzx,Azzy,Azzz
  real*8 :: gxxx,gxxy,gxxz
  real*8 :: gxyx,gxyy,gxyz
  real*8 :: gxzx,gxzy,gxzz
  real*8 :: gyyx,gyyy,gyyz
  real*8 :: gyzx,gyzy,gyzz
  real*8 :: gzzx,gzzy,gzzz
  real*8 :: gxxxx,gxxxy,gxxxz,gxxyy,gxxyz,gxxzz
  real*8 :: gxyxx,gxyxy,gxyxz,gxyyy,gxyyz,gxyzz
  real*8 :: gxzxx,gxzxy,gxzxz,gxzyy,gxzyz,gxzzz
  real*8 :: gyyxx,gyyxy,gyyxz,gyyyy,gyyyz,gyyzz
  real*8 :: gyzxx,gyzxy,gyzxz,gyzyy,gyzyz,gyzzz
  real*8 :: gzzxx,gzzxy,gzzxz,gzzyy,gzzyz,gzzzz
  real*8 :: Gamxx,Gamxy,Gamxz
  real*8 :: Gamyx,Gamyy,Gamyz
  real*8 :: Gamzx,Gamzy,Gamzz
  real*8 :: Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz
  real*8 :: Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz
  real*8 :: Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz
  real*8 :: Gamxa,Gamya,Gamza
  real*8 :: CAZxx,CAZxy,CAZxz
  real*8 :: CAZyx,CAZyy,CAZyz
  real*8 :: CAZzx,CAZzy,CAZzz
! tilted A^k_iA_kj
  real*8 :: AAxx,AAxy,AAxz,AAyy,AAyz,AAzz
  real*8 :: Ainvxx,Ainvxy,Ainvxz,Ainvyy,Ainvyz,Ainvzz
  real*8 :: liegxx,liegxy,liegxz,liegyy,liegyz,liegzz
  real*8 :: fxx,fxy,fxz,fyy,fyz,fzz
  real*8 :: TFxx,TFxy,TFxz,TFyy,TFyz,TFzz

  real*8 :: MapleGenVar1,MapleGenVar2,MapleGenVar3,MapleGenVar4
  real*8 :: f,betas
  
  real*8,dimension(3) ::SSS,AAS,ASA,SAA,ASS,SAS,SSA
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  real*8, parameter :: ZEO = 0.d0, ONE = 1.d0, TWO=2.d0,HALF=0.5d0
  real*8,parameter::TINYRR=1.d-14
! in order to synchronize the following parameters with Z4c_rhs calculation, we
! call a routine
  real*8 :: muL,tmuSL,tmuST
  real*8 :: kappa1,kappa2,kappa3,FF,eta

  real*8,parameter :: ha=0.d0,thbs=0.d0,hu=0.d0,hw=0.d0,Rhpsi0=0.d0,Ihpsi0=0.d0

  call get_Z4cparameters(kappa1,kappa2,kappa3,FF,eta)

  dR = R(2) - R(1)

  kmax = ex(3)

  kmin = 1

layer(1:3,:) = 1
layer(4:6,:) =-1
 
if(dabs(R(ex(3))-zmax) < dR)then
  layer(1,3) = 1
  layer(2,3) = 1
  layer(3,3) = ex(3) - CPBC_ghost_width
  layer(4,3) = ex(1)
  layer(5,3) = ex(2)
  layer(6,3) = ex(3) - CPBC_ghost_width
endif

if(dabs(R(1)-zmin) < dR)then
   layer(1,6) = 1
   layer(2,6) = 1
   layer(3,6) = 1
   layer(4,6) = ex(1)
   layer(5,6) = ex(2)
   layer(6,6) = 1
endif
! fix BD
  gp = 6

  gont = any( layer(:,gp) == - 1 )
 
   if( .not. gont ) then

    do k = layer(3,gp), layer(6,gp)
     do j = layer(2,gp), layer(5,gp)
      do i = layer(1,gp), layer(4,gp)
! z direction   
        TZ_rhs(i,j,k) = ZEO
        chi_rhs(i,j,k) = ZEO
        trK_rhs(i,j,k) = ZEO
        gxx_rhs(i,j,k) = ZEO
        gxy_rhs(i,j,k) = ZEO
        gxz_rhs(i,j,k) = ZEO
        gyy_rhs(i,j,k) = ZEO
        gyz_rhs(i,j,k) = ZEO
        gzz_rhs(i,j,k) = ZEO
        Axx_rhs(i,j,k) = ZEO
        Axy_rhs(i,j,k) = ZEO
        Axz_rhs(i,j,k) = ZEO
        Ayy_rhs(i,j,k) = ZEO
        Ayz_rhs(i,j,k) = ZEO
        Azz_rhs(i,j,k) = ZEO
        Gamx_rhs(i,j,k) = ZEO
        Gamy_rhs(i,j,k) = ZEO
        Gamz_rhs(i,j,k) = ZEO
        Lap_rhs(i,j,k) = ZEO
        betax_rhs(i,j,k) = ZEO
        betay_rhs(i,j,k) = ZEO
        betaz_rhs(i,j,k) = ZEO
        dtSfx_rhs(i,j,k) = ZEO
        dtSfy_rhs(i,j,k) = ZEO
        dtSfz_rhs(i,j,k) = ZEO
      enddo
     enddo
    enddo
   endif

! constraint preserving BD
  gp = 3

  gont = any( layer(:,gp) == - 1 )
 
  if( .not. gont ) then

! cpbc real starts

  alpn1 = Lap + ONE
  chin1 = chi + ONE
  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE  
! invert tilted metric
  gupzz =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz - &
           gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz
  gupxx =   ( gyy * gzz - gyz * gyz ) / gupzz
  gupxy = - ( gxy * gzz - gyz * gxz ) / gupzz
  gupxz =   ( gxy * gyz - gyy * gxz ) / gupzz
  gupyy =   ( gxx * gzz - gxz * gxz ) / gupzz
  gupyz = - ( gxx * gyz - gxy * gxz ) / gupzz
  gupzz =   ( gxx * gyy - gxy * gxy ) / gupzz
! tetrad for 2+1 decomposation
  do i=1,ex(1)
  do j=1,ex(2)
  do k=1,ex(3)
     if(abs(X(i,j,k)) < TINYRR .and. abs(Y(i,j,k)) < TINYRR .and. abs(Z(i,j,k)) < TINYRR)then
        vx(i,j,k) = TINYRR
        vy(i,j,k) = TINYRR
        vz(i,j,k) = TINYRR
     else
        vx(i,j,k) = X(i,j,k)
        vy(i,j,k) = Y(i,j,k)
        vz(i,j,k) = Z(i,j,k)
     endif
     if(abs(X(i,j,k)) < TINYRR .and. abs(Y(i,j,k)) < TINYRR)then
        ux(i,j,k) = - TINYRR
        uy(i,j,k) = TINYRR
        uz(i,j,k) = ZEO
        wx(i,j,k) = TINYRR*Z(i,j,k)
        wy(i,j,k) = TINYRR*Z(i,j,k)
        wz(i,j,k) = -2*TINYRR*TINYRR
     else
        ux(i,j,k) = - Y(i,j,k)
        uy(i,j,k) = X(i,j,k)
        uz(i,j,k) = ZEO
        wx(i,j,k) = X(i,j,k)*Z(i,j,k)
        wy(i,j,k) = Y(i,j,k)*Z(i,j,k)
        wz(i,j,k) = -(X(i,j,k)*X(i,j,k) + Y(i,j,k)*Y(i,j,k))
     endif
  enddo
  enddo
  enddo

! v^i corresponds to s^i  
  fx = vx
  fy = vy
  fz = vz
  slx = vx
  sly = vy
  slz = vz
  vx = gupxx*fx + gupxy*fy + gupxz*fz
  vy = gupxy*fx + gupyy*fy + gupyz*fz
  vz = gupxz*fx + gupyz*fy + gupzz*fz

  fx = gxx*vx*vx + gyy*vy*vy + gzz*vz*vz &
     +(gxy*vx*vy + gxz*vx*vz + gyz*vy*vz)*TWO
  fx = dsqrt(fx*chin1)
  vx = vx/fx
  vy = vy/fx
  vz = vz/fx
  slx = slx/fx
  sly = sly/fx
  slz = slz/fx
! 2+1: 1->u, 2->w
  fx = gxx*vx*ux + gxy*vx*uy + gxz*vx*uz + &
       gxy*vy*ux + gyy*vy*uy + gyz*vy*uz + &
       gxz*vz*ux + gyz*vz*uy + gzz*vz*uz
  fx = fx/chin1
  ux = ux - fx*vx
  uy = uy - fx*vy
  uz = uz - fx*vz
  fx = gxx*ux*ux + gyy*uy*uy + gzz*uz*uz &
     +(gxy*ux*uy + gxz*ux*uz + gyz*uy*uz)*TWO
  fx = dsqrt(fx/chin1) 
  ux = ux/fx
  uy = uy/fx
  uz = uz/fx
  ulx = (gxx*ux+gxy*uy+gxz*uz)/chin1
  uly = (gxy*ux+gyy*uy+gyz*uz)/chin1
  ulz = (gxz*ux+gyz*uy+gzz*uz)/chin1

  fx = gxx*vx*wx + gxy*vx*wy + gxz*vx*wz + &
       gxy*vy*wx + gyy*vy*wy + gyz*vy*wz + &
       gxz*vz*wx + gyz*vz*wy + gzz*vz*wz
  fx = fx/chin1
  wx = wx - fx*vx
  wy = wy - fx*vy
  wz = wz - fx*vz
  fx = gxx*ux*wx + gxy*ux*wy + gxz*ux*wz + &
       gxy*uy*wx + gyy*uy*wy + gyz*uy*wz + &
       gxz*uz*wx + gyz*uz*wy + gzz*uz*wz
  fx = fx/chin1
  wx = wx - fx*ux
  wy = wy - fx*uy
  wz = wz - fx*uz
  fx = gxx*wx*wx + gyy*wy*wy + gzz*wz*wz &
     +(gxy*wx*wy + gxz*wx*wz + gyz*wy*wz)*TWO
  fx = dsqrt(fx/chin1)
  wx = wx/fx
  wy = wy/fx
  wz = wz/fx
  wlx = (gxx*wx+gxy*wy+gxz*wz)/chin1
  wly = (gxy*wx+gyy*wy+gyz*wz)/chin1
  wlz = (gxz*wx+gyz*wy+gzz*wz)/chin1
!~ end tetrad  

  qupxx = gupxx*chin1 - vx*vx
  qupxy = gupxy*chin1 - vx*vy
  qupxz = gupxz*chin1 - vx*vz
  qupyy = gupyy*chin1 - vy*vy
  qupyz = gupyz*chin1 - vy*vz
  qupzz = gupzz*chin1 - vz*vz

  qxx = gxx/chin1 - slx*slx
  qxy = gxy/chin1 - slx*sly
  qxz = gxz/chin1 - slx*slz
  qyy = gyy/chin1 - sly*sly
  qyz = gyz/chin1 - sly*slz
  qzz = gzz/chin1 - slz*slz

  qulxx = ONE - vx*slx
  qulyy = ONE - vy*sly
  qulzz = ONE - vz*slz
  qulxy =     - vx*sly
  qulyx =     - vy*slx
  qulxz =     - vx*slz
  qulzx =     - vz*slx
  qulyz =     - vy*slz
  qulzy =     - vz*sly

    do k = layer(3,gp), layer(6,gp)
     do j = layer(2,gp), layer(5,gp)
      do i = layer(1,gp), layer(4,gp)
!calculate the involved derivatives      
         call point_fderivs_shc(ex,trK,Kx,Ky,Kz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,chi,chix,chiy,chiz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Lap,Lapx,Lapy,Lapz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,TZ,TZx,TZy,TZz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Gamx,Gamxx,Gamxy,Gamxz,crho,sigma,R,ANTI,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Gamy,Gamyx,Gamyy,Gamyz,crho,sigma,R,SYM,ANTI,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Gamz,Gamzx,Gamzy,Gamzz,crho,sigma,R,SYM,SYM,ANTI,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,betax,sfxx,sfxy,sfxz,crho,sigma,R,ANTI,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,betay,sfyx,sfyy,sfyz,crho,sigma,R,SYM,ANTI,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,betaz,sfzx,sfzy,sfzz,crho,sigma,R,SYM,SYM,ANTI,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Axx,Axxx,Axxy,Axxz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Axy,Axyx,Axyy,Axyz,crho,sigma,R,ANTI,ANTI,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Axz,Axzx,Axzy,Axzz,crho,sigma,R,ANTI,SYM,ANTI,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Ayy,Ayyx,Ayyy,Ayyz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Ayz,Ayzx,Ayzy,Ayzz,crho,sigma,R,SYM,ANTI,ANTI,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Azz,Azzx,Azzy,Azzz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,dxx,gxxx,gxxy,gxxz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,gxy,gxyx,gxyy,gxyz,crho,sigma,R,ANTI,ANTI,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,gxz,gxzx,gxzy,gxzz,crho,sigma,R,ANTI,SYM,ANTI,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,dyy,gyyx,gyyy,gyyz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,gyz,gyzx,gyzy,gyzz,crho,sigma,R,SYM,ANTI,ANTI,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,dzz,gzzx,gzzy,gzzz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)

  liegxx = betax(i,j,k)*gxxx+gxx(i,j,k)*sfxx+betay(i,j,k)*gxxy-gxx(i,j,k)*sfyy+2.0*sfyx*gxy(i,j,k)+betaz(i,j,k)*gxxz-gxx(i,j,k)*sfzz+2.0*sfzx*gxz(i,j,k)
  liegxy = betax(i,j,k)*gxyx+sfxy*gxx(i,j,k)+betay(i,j,k)*gxyy+sfyx*gyy(i,j,k)+betaz(i,j,k)*gxyz-gxy(i,j,k)*sfzz+sfzx*gyz(i,j,k)+sfzy*gxz(i,j,k)
  liegxz = betax(i,j,k)*gxzx+sfxz*gxx(i,j,k)+betay(i,j,k)*gxzy-gxz(i,j,k)*sfyy+sfyx*gyz(i,j,k)+sfyz*gxy(i,j,k)+betaz(i,j,k)*gxzz+sfzx*gzz(i,j,k)
  liegyy = betax(i,j,k)*gyyx-gyy(i,j,k)*sfxx+2.0*sfxy*gxy(i,j,k)+betay(i,j,k)*gyyy+gyy(i,j,k)*sfyy+betaz(i,j,k)*gyyz-gyy(i,j,k)*sfzz+2.0*sfzy*gyz(i,j,k)
  liegyz = betax(i,j,k)*gyzx-gyz(i,j,k)*sfxx+sfxy*gxz(i,j,k)+sfxz*gxy(i,j,k)+betay(i,j,k)*gyzy+sfyz*gyy(i,j,k)+betaz(i,j,k)*gyzz+sfzy*gzz(i,j,k)
  liegzz = betax(i,j,k)*gzzx-gzz(i,j,k)*sfxx+2.0*sfxz*gxz(i,j,k)+betay(i,j,k)*gzzy-gzz(i,j,k)*sfyy+2.0*sfyz*gyz(i,j,k)+betaz(i,j,k)*gzzz+gzz(i,j,k)*sfzz

         call point_fdderivs_shc(ex,dxx,gxxxx,gxxxy,gxxxz,gxxyy,gxxyz,gxxzz,crho,sigma,R,SYM ,SYM ,SYM ,Symmetry,0,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)
         call point_fdderivs_shc(ex,dyy,gyyxx,gyyxy,gyyxz,gyyyy,gyyyz,gyyzz,crho,sigma,R,SYM ,SYM ,SYM ,Symmetry,0,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)
         call point_fdderivs_shc(ex,dzz,gzzxx,gzzxy,gzzxz,gzzyy,gzzyz,gzzzz,crho,sigma,R,SYM ,SYM ,SYM ,Symmetry,0,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)
         call point_fdderivs_shc(ex,gxy,gxyxx,gxyxy,gxyxz,gxyyy,gxyyz,gxyzz,crho,sigma,R,ANTI,ANTI,SYM ,Symmetry,0,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)
         call point_fdderivs_shc(ex,gxz,gxzxx,gxzxy,gxzxz,gxzyy,gxzyz,gxzzz,crho,sigma,R,ANTI,SYM ,ANTI,Symmetry,0,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)
         call point_fdderivs_shc(ex,gyz,gyzxx,gyzxy,gyzxz,gyzyy,gyzyz,gyzzz,crho,sigma,R,SYM ,SYM ,SYM ,Symmetry,0,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)

      MapleGenVar3 = gupxy(i,j,k)*gupxy(i,j,k)*gxxyy-4.0*gupxz(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*gxzx*gxzz+gupxy(i,j,k)*&
gupxy(i,j,k)*gxyxy-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*gxyy*gxyy+gupxx(i,j,k)*gupyz(i,j,k)*gxyxz+gupxz(i,j,k)*gupxz(i,j,k)*gxxzz&
-2.0*gupxx(i,j,k)*gupxx(i,j,k)*gupzz(i,j,k)*gxxz*gxxz-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*gxzz*gxzz-2.0*gupxx(i,j,k)*&
gupxx(i,j,k)*gupyy(i,j,k)*gxxy*gxxy+gupxx(i,j,k)*gupyy(i,j,k)*gxyxy+gupxz(i,j,k)*gupzz(i,j,k)*gxzzz-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupxx(i,j,k)&
*gxzx*gxzx+2.0*gupxx(i,j,k)*gupxz(i,j,k)*gxxxz+gupxz(i,j,k)*gupyy(i,j,k)*gxyyz+gupxx(i,j,k)*gupyz(i,j,k)*gxzxy+gupxy(i,j,k)*&
gupzz(i,j,k)*gxzyz+2.0*gupxy(i,j,k)*gupxz(i,j,k)*gxxyz+2.0*gupxx(i,j,k)*gupxy(i,j,k)*gxxxy+gupxz(i,j,k)*gupxy(i,j,k)*gxyxz+gupxz(i,j,k)&
*gupyz(i,j,k)*gxyzz-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxx(i,j,k)*gxyx*gxyx-4.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*gxyx*gxyy+&
gupxx(i,j,k)*gupxx(i,j,k)*gxxxx+gupxx(i,j,k)*gupxz(i,j,k)*gxzxx-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gxzy*gxzy+gupxx(i,j,k)*gupzz(i,j,k)&
*gxzxz
      MapleGenVar4 = gupxz(i,j,k)*gupyz(i,j,k)*gxzyz+gupxy(i,j,k)*gupyz(i,j,k)*gxyyz+gupxy(i,j,k)*gupyz(i,j,k)*gxzyy+&
gupxx(i,j,k)*gupxy(i,j,k)*gxyxx+gupxy(i,j,k)*gupxz(i,j,k)*gxzxy-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*gxyz*gxyz+gupxy(i,j,k)*gupyy(i,j,k)&
*gxyyy+gupxz(i,j,k)*gupxz(i,j,k)*gxzxz-2.0*gupxx(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*gxyx*gxxz-2.0*gupxx(i,j,k)*gupxx(i,j,k)*gupxx(i,j,k)*&
gxxx*gxxx-6.0*gupxx(i,j,k)*gupzz(i,j,k)*gupxy(i,j,k)*gxxz*gxyz-6.0*gupxx(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gxxx*gxyx-6.0*&
gupxx(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gxxx*gxzx
      MapleGenVar2 = MapleGenVar4-gupxx(i,j,k)*gupxx(i,j,k)*gupyy(i,j,k)*(gxxx*gyyx+gxyx*gxyx)-2.0*&
gupxx(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*(gxxx*gyzx+gxyx*gxzx)-gupxx(i,j,k)*gupxx(i,j,k)*gupzz(i,j,k)*(gxxx*gzzx+gxzx*gxzx)&
-4.0*gupxx(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gxxx*gxxy-2.0*gupxx(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*gxxx*gxyy-4.0*gupxx(i,j,k)*&
gupxy(i,j,k)*gupxy(i,j,k)*(gxxx*gxyy+gxyx*gxxy)-gupxx(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*(gxxx*gyyy+gxyx*gxyy)-gupxx(i,j,k)&
*gupxy(i,j,k)*gupyz(i,j,k)*(gxxx*gyzy+gxyx*gxzy)-4.0*gupxx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxxx*gxzy+gxzx*gxxy)-&
gupxx(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxxx*gyzy+gxzx*gxyy)-gupxx(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gxxx*gzzy+gxzx*gxzy)&
-4.0*gupxx(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gxxx*gxxz-2.0*gupxx(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*gxxx*gxzz+MapleGenVar3
      MapleGenVar4 = -4.0*gupxx(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*(gxxx*gxyz+gxyx*gxxz)-gupxx(i,j,k)*gupxz(i,j,k)*&
gupyy(i,j,k)*(gxxx*gyyz+gxyx*gxyz)-gupxx(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxxx*gyzz+gxyx*gxzz)-4.0*gupxx(i,j,k)*&
gupxz(i,j,k)*gupxz(i,j,k)*(gxxx*gxzz+gxzx*gxxz)-gupxx(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxxx*gyzz+gxzx*gxyz)-gupxx(i,j,k)&
*gupxz(i,j,k)*gupzz(i,j,k)*(gxxx*gzzz+gxzx*gxzz)-2.0*gupxx(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*gxyx*gxxy-gupxx(i,j,k)*gupxy(i,j,k)*&
gupyy(i,j,k)*(gxxy*gyyx+gxyx*gxyy)-gupxx(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxxy*gyzx+gxzx*gxyy)-gupxx(i,j,k)*gupxy(i,j,k)&
*gupyz(i,j,k)*(gxxy*gyzx+gxyx*gxzy)-gupxx(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gxxy*gzzx+gxzx*gxzy)-gupxx(i,j,k)*&
gupyy(i,j,k)*gupyy(i,j,k)*(gxxy*gyyy+gxyy*gxyy)-2.0*gupxx(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*(gxxy*gyzy+gxyy*gxzy)
      MapleGenVar3 = MapleGenVar4-gupxx(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*(gxxy*gzzy+gxzy*gxzy)-4.0*&
gupxx(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*gxxy*gxxz-4.0*gupxx(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*(gxxy*gxyz+gxyy*gxxz)-gupxx(i,j,k)*&
gupyz(i,j,k)*gupyy(i,j,k)*(gxxy*gyyz+gxyy*gxyz)-gupxx(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxxy*gyzz+gxyy*gxzz)-4.0*&
gupxx(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxxy*gxzz+gxzy*gxxz)-gupxx(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxxy*gyzz+gxzy*gxyz)&
-gupxx(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(gxxy*gzzz+gxzy*gxzz)-2.0*gupxx(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*gxzx*gxxz-gupxx(i,j,k)*&
gupxz(i,j,k)*gupyy(i,j,k)*(gxxz*gyyx+gxyx*gxyz)-gupxx(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxxz*gyzx+gxzx*gxyz)-gupxx(i,j,k)&
*gupxz(i,j,k)*gupyz(i,j,k)*(gxxz*gyzx+gxyx*gxzz)-gupxx(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gxxz*gzzx+gxzx*gxzz)-&
gupxx(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*(gxxz*gyyy+gxyy*gxyz)
      MapleGenVar4 = -gupxx(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxxz*gyzy+gxzy*gxyz)-gupxx(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)&
*(gxxz*gyzy+gxyy*gxzz)-gupxx(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(gxxz*gzzy+gxzy*gxzz)-gupxx(i,j,k)*gupzz(i,j,k)*&
gupyy(i,j,k)*(gxxz*gyyz+gxyz*gxyz)-2.0*gupxx(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*(gxxz*gyzz+gxyz*gxzz)-gupxx(i,j,k)*&
gupzz(i,j,k)*gupzz(i,j,k)*(gxxz*gzzz+gxzz*gxzz)-gupxy(i,j,k)*gupxy(i,j,k)*gupxx(i,j,k)*(gxxx*gyyx+gxyx*gxyx)-2.0*&
gupxy(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*(gxxx*gyzx+gxyx*gxzx)-gupxy(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*(gxyx*gyzx+gxzx*gyyx)&
-gupxy(i,j,k)*gupxx(i,j,k)*gupzz(i,j,k)*(gxyx*gzzx+gxzx*gyzx)-gupxy(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*(gxxx*gyyy+gxyx*gxyy&
)-gupxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxxx*gyzy+gxyx*gxzy)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*gxyx*gyyy-2.0*&
gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gxyx*gyzy
      MapleGenVar1 = MapleGenVar4-4.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxyx*gxzy+gxzx*gxyy)-&
gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxyx*gyzy+gxzx*gyyy)-gupxy(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gxyx*gzzy+gxzx*gyzy)&
-gupxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxxx*gyyz+gxyx*gxyz)-gupxy(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxxx*gyzz+gxyx*gxzz&
)-4.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxyx*gxyz-2.0*gupxx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxxx*gxzy-2.0*gupxx(i,j,k)*&
gupxz(i,j,k)*gupxy(i,j,k)*gxxx*gxyz-2.0*gupxx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxzx*gxxy-6.0*gupxx(i,j,k)*gupyy(i,j,k)*gupxy(i,j,k)*&
gxxy*gxyy-6.0*gupxx(i,j,k)*gupyy(i,j,k)*gupxz(i,j,k)*gxxy*gxzy-2.0*gupxx(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*gxxy*gxyz+&
MapleGenVar3+MapleGenVar2
      MapleGenVar4 = MapleGenVar1-2.0*gupxx(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*gxxy*gxzz-2.0*gupxx(i,j,k)*&
gupyz(i,j,k)*gupxy(i,j,k)*gxyy*gxxz-2.0*gupxx(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*gxzy*gxxz-6.0*gupxx(i,j,k)*gupzz(i,j,k)*gupxz(i,j,k)*&
gxxz*gxzz-2.0*gupxy(i,j,k)*gupxx(i,j,k)*gupyy(i,j,k)*gxyx*gyyx-2.0*gupxy(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*gxyx*gyzx-4.0*&
gupxy(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gxyx*gxzx-2.0*gupxy(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gxyx*gyyz-2.0*gupxy(i,j,k)*gupxz(i,j,k)*&
gupyz(i,j,k)*gxyx*gyzz-4.0*gupxy(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxyx*gxzz+gxzx*gxyz)-gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(&
gxyx*gyzz+gxzx*gyyz)-gupxy(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gxyx*gzzz+gxzx*gyzz)
      MapleGenVar3 = MapleGenVar4-gupxy(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*(gxxy*gyyx+gxyx*gxyy)-gupxy(i,j,k)&
*gupxy(i,j,k)*gupxz(i,j,k)*(gxxy*gyzx+gxzx*gxyy)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*gxyy*gyyx-2.0*gupxy(i,j,k)*&
gupxy(i,j,k)*gupyz(i,j,k)*gxyy*gyzx-gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxyy*gyzx+gxzy*gyyx)-gupxy(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)&
*(gxyy*gzzx+gxzy*gyzx)-gupxy(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*(gxxy*gyyy+gxyy*gxyy)-2.0*gupxy(i,j,k)*gupyy(i,j,k)*&
gupxz(i,j,k)*(gxxy*gyzy+gxyy*gxzy)-2.0*gupxy(i,j,k)*gupyy(i,j,k)*gupyy(i,j,k)*gxyy*gyyy-gupxy(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*(&
gxyy*gyzy+gxzy*gyyy)-gupxy(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*(gxyy*gzzy+gxzy*gyzy)-gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*&
(gxxy*gyyz+gxyy*gxyz)-gupxy(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxxy*gyzz+gxyy*gxzz)
      MapleGenVar4 = MapleGenVar3-4.0*gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gxyy*gxyz-2.0*gupxy(i,j,k)*&
gupyz(i,j,k)*gupyz(i,j,k)*gxyy*gyzz-4.0*gupxy(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxyy*gxzz+gxzy*gxyz)-gupxy(i,j,k)*gupyz(i,j,k)*&
gupyz(i,j,k)*(gxyy*gyzz+gxzy*gyyz)-gupxy(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(gxyy*gzzz+gxzy*gyzz)-gupxy(i,j,k)*gupxy(i,j,k)&
*gupxz(i,j,k)*(gxxz*gyyx+gxyx*gxyz)-gupxy(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxxz*gyzx+gxzx*gxyz)-gupxy(i,j,k)*&
gupxz(i,j,k)*gupyz(i,j,k)*(gxyz*gyzx+gxzz*gyyx)-gupxy(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gxyz*gzzx+gxzz*gyzx)-gupxy(i,j,k)&
*gupxy(i,j,k)*gupyz(i,j,k)*(gxxz*gyyy+gxyy*gxyz)-gupxy(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxxz*gyzy+gxzy*gxyz)-2.0*&
gupxy(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*gxyz*gyzy-gupxy(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxyz*gyzy+gxzz*gyyy)
      MapleGenVar2 = MapleGenVar4-gupxy(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(gxyz*gzzy+gxzz*gyzy)-gupxy(i,j,k)&
*gupxy(i,j,k)*gupzz(i,j,k)*(gxxz*gyyz+gxyz*gxyz)-2.0*gupxy(i,j,k)*gupzz(i,j,k)*gupxz(i,j,k)*(gxxz*gyzz+gxyz*gxzz)-&
gupxy(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*(gxyz*gyzz+gxzz*gyyz)-gupxy(i,j,k)*gupzz(i,j,k)*gupzz(i,j,k)*(gxyz*gzzz+gxzz*gyzz)&
-gupxz(i,j,k)*gupxz(i,j,k)*gupxx(i,j,k)*(gxxx*gzzx+gxzx*gxzx)-gupxz(i,j,k)*gupxx(i,j,k)*gupyy(i,j,k)*(gxyx*gyzx+gxzx*gyyx&
)-gupxz(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*(gxyx*gzzx+gxzx*gyzx)-gupxz(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*(gxxx*gyzy+gxzx*&
gxyy)-gupxz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*(gxxx*gzzy+gxzx*gxzy)-gupxz(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*(gxyx*gyzy+gxzx&
*gyyy)-gupxz(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxyx*gzzy+gxzx*gyzy)-4.0*gupxz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*gxzx*gxzy-&
gupxz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*(gxxx*gyzz+gxzx*gxyz)
      MapleGenVar4 = MapleGenVar2-gupxz(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxxx*gzzz+gxzx*gxzz)-gupxz(i,j,k)&
*gupxz(i,j,k)*gupyy(i,j,k)*(gxyx*gyzz+gxzx*gyyz)-gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxyx*gzzz+gxzx*gyzz)-2.0*&
gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gxzx*gyzz-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*gxzx*gzzz-gupxz(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*(&
gxxy*gyzx+gxyx*gxzy)-gupxz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*(gxxy*gzzx+gxzx*gxzy)-gupxz(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*&
(gxyy*gyzx+gxzy*gyyx)-gupxz(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxyy*gzzx+gxzy*gyzx)-gupxz(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)&
*(gxxy*gzzy+gxzy*gxzy)-gupxz(i,j,k)*gupyy(i,j,k)*gupyy(i,j,k)*(gxyy*gyzy+gxzy*gyyy)-gupxz(i,j,k)*gupyy(i,j,k)*&
gupyz(i,j,k)*(gxyy*gzzy+gxzy*gyzy)
      MapleGenVar3 = MapleGenVar4-gupxz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*(gxxy*gyzz+gxzy*gxyz)-gupxz(i,j,k)&
*gupxz(i,j,k)*gupyz(i,j,k)*(gxxy*gzzz+gxzy*gxzz)-gupxz(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*(gxyy*gyzz+gxzy*gyyz)-&
gupxz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxyy*gzzz+gxzy*gyzz)-4.0*gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gxzy*gxzz-2.0*&
gupxz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*gxzy*gyzz-gupxz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*(gxxz*gyzx+gxyx*gxzz)-gupxz(i,j,k)*gupxz(i,j,k)&
*gupxz(i,j,k)*(gxxz*gzzx+gxzx*gxzz)-gupxz(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*(gxyz*gyzx+gxzz*gyyx)-gupxz(i,j,k)*&
gupxz(i,j,k)*gupyz(i,j,k)*(gxyz*gzzx+gxzz*gyzx)-2.0*gupxy(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gxyy*gyzy-4.0*gupxy(i,j,k)*&
gupyy(i,j,k)*gupxz(i,j,k)*gxyy*gxzy-2.0*gupxy(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*gxyy*gyyz-2.0*gupxy(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*&
gxyz*gyyx
      MapleGenVar4 = MapleGenVar3-2.0*gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gxyz*gyzx-2.0*gupxy(i,j,k)*&
gupyz(i,j,k)*gupyy(i,j,k)*gxyz*gyyy-2.0*gupxy(i,j,k)*gupzz(i,j,k)*gupyy(i,j,k)*gxyz*gyyz-2.0*gupxy(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*&
gxyz*gyzz-4.0*gupxy(i,j,k)*gupzz(i,j,k)*gupxz(i,j,k)*gxyz*gxzz-2.0*gupxz(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*gxzx*gyzx-2.0*&
gupxz(i,j,k)*gupxx(i,j,k)*gupzz(i,j,k)*gxzx*gzzx-2.0*gupxz(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gxzx*gyzy-2.0*gupxz(i,j,k)*gupxy(i,j,k)*&
gupzz(i,j,k)*gxzx*gzzy-2.0*gupxz(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gxzy*gyzx-2.0*gupxz(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*gxzy*gzzx&
-2.0*gupxz(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gxzy*gyzy-2.0*gupxz(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*gxzy*gzzy
      CAZxx = Gamxx - (MapleGenVar4-2.0*gupxz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gxzy*gzzz-2.0*gupxz(i,j,k)*&
gupxz(i,j,k)*gupyz(i,j,k)*gxzz*gyzx-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*gxzz*gzzx-gupxz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*(gxxz*&
gyzy+gxyy*gxzz)-gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxxz*gzzy+gxzy*gxzz)-gupxz(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*(gxyz&
*gyzy+gxzz*gyyy)-gupxz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxyz*gzzy+gxzz*gyzy)-2.0*gupxz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*&
gxzz*gyzy-gupxz(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gxxz*gzzz+gxzz*gxzz)-gupxz(i,j,k)*gupzz(i,j,k)*gupyy(i,j,k)*(gxyz*gyzz+&
gxzz*gyyz)-gupxz(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*(gxyz*gzzz+gxzz*gyzz)-2.0*gupxz(i,j,k)*gupzz(i,j,k)*gupzz(i,j,k)*gxzz*&
gzzz-2.0*gupxz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gxzz*gzzy-2.0*gupxz(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*gxzz*gyzz)
      MapleGenVar3 = -2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*gxzy*gyzz+gupxy(i,j,k)*gupyz(i,j,k)*gxzxy+2.0*&
gupxy(i,j,k)*gupyy(i,j,k)*gxyxy+gupxy(i,j,k)*gupxz(i,j,k)*gxzxx+gupxy(i,j,k)*gupxy(i,j,k)*gxyxx-2.0*gupxy(i,j,k)*gupxx(i,j,k)*gupxx(i,j,k)*&
gxxx*gxxx+2.0*gupxy(i,j,k)*gupyz(i,j,k)*gxyxz-2.0*gupyy(i,j,k)*gupyy(i,j,k)*gupyy(i,j,k)*gxyy*gyyy+gupyz(i,j,k)*gupxz(i,j,k)*&
gxxzz+gupxy(i,j,k)*gupxy(i,j,k)*gxxxy+gupxy(i,j,k)*gupxx(i,j,k)*gxxxx+gupyy(i,j,k)*gupxx(i,j,k)*gxxxy+gupxy(i,j,k)*gupxz(i,j,k)*gxxxz+&
gupxy(i,j,k)*gupzz(i,j,k)*gxzxz-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*gxzz*gyzy+gupyy(i,j,k)*gupxy(i,j,k)*gxxyy+gupyz(i,j,k)*gupxz(i,j,k)&
*gxzxz+2.0*gupyy(i,j,k)*gupyz(i,j,k)*gxyyz+gupyz(i,j,k)*gupzz(i,j,k)*gxzzz+gupyy(i,j,k)*gupzz(i,j,k)*gxzyz-2.0*gupxy(i,j,k)*&
gupxy(i,j,k)*gupxy(i,j,k)*gxyx*gxxy-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*gxxx*gxyy+gupyy(i,j,k)*gupxz(i,j,k)*gxxyz+gupyz(i,j,k)*&
gupyz(i,j,k)*gxzyz-2.0*gupyy(i,j,k)*gupyy(i,j,k)*gupxy(i,j,k)*gxyy*gxyy
      MapleGenVar4 = gupyz(i,j,k)*gupxx(i,j,k)*gxxxz+gupyz(i,j,k)*gupxy(i,j,k)*gxxyz+gupyy(i,j,k)*gupxz(i,j,k)*gxzxy+&
gupyz(i,j,k)*gupyz(i,j,k)*gxyzz+gupyy(i,j,k)*gupyz(i,j,k)*gxzyy+gupyy(i,j,k)*gupyy(i,j,k)*gxyyy-gupyz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*(gxyy*&
gzzx+gxzy*gyzx)-4.0*gupxy(i,j,k)*gupxz(i,j,k)*gupxx(i,j,k)*gxxx*gxxz-4.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxx(i,j,k)*gxxx*gxyx&
-2.0*gupxy(i,j,k)*gupxx(i,j,k)*gupyy(i,j,k)*(gxxx*gyyx+gxyx*gxyx)-3.0*gupxy(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*(gxxx*gyzx+&
gxyx*gxzx)-gupxy(i,j,k)*gupxx(i,j,k)*gupzz(i,j,k)*(gxxx*gzzx+gxzx*gxzx)-4.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxx(i,j,k)*gxxx*&
gxxy
      MapleGenVar2 = MapleGenVar4-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxxx*gxzy-2.0*gupxy(i,j,k)*&
gupxy(i,j,k)*gupxy(i,j,k)*(gxxx*gxyy+gxyx*gxxy)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*(gxxx*gyyy+gxyx*gxyy)-&
gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxxx*gyzy+gxyx*gxzy)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxxx*gxzy+gxzx*&
gxxy)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxxx*gyzy+gxzx*gxyy)-gupxy(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gxxx*gzzy+&
gxzx*gxzy)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxxx*gxyz-2.0*gupxy(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*gxxx*gxzz-2.0*&
gupxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxxx*gxyz+gxyx*gxxz)-2.0*gupxy(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*(gxxx*gyyz+gxyx*&
gxyz)-gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxxx*gyzz+gxyx*gxzz)-2.0*gupxy(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxxx*gxzz+&
gxzx*gxxz)+MapleGenVar3
      MapleGenVar4 = -2.0*gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxxx*gyzz+gxzx*gxyz)-gupxy(i,j,k)*gupxz(i,j,k)*&
gupzz(i,j,k)*(gxxx*gzzz+gxzx*gxzz)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxzx*gxxy-2.0*gupxy(i,j,k)*gupxy(i,j,k)*&
gupyy(i,j,k)*(gxxy*gyyx+gxyx*gxyy)-gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxxy*gyzx+gxzx*gxyy)-2.0*gupxy(i,j,k)*&
gupxy(i,j,k)*gupyz(i,j,k)*(gxxy*gyzx+gxyx*gxzy)-gupxy(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gxxy*gzzx+gxzx*gxzy)-2.0*&
gupxy(i,j,k)*gupyy(i,j,k)*gupxx(i,j,k)*gxxy*gxxy-4.0*gupxy(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*gxxy*gxyy-2.0*gupxy(i,j,k)*gupyy(i,j,k)*&
gupyy(i,j,k)*(gxxy*gyyy+gxyy*gxyy)-3.0*gupxy(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*(gxxy*gyzy+gxyy*gxzy)-gupxy(i,j,k)*&
gupyy(i,j,k)*gupzz(i,j,k)*(gxxy*gzzy+gxzy*gxzy)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gxxy*gxyz
      MapleGenVar3 = MapleGenVar4-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxxy*gxyz+gxyy*gxxz)&
-2.0*gupxy(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*(gxxy*gyyz+gxyy*gxyz)-gupxy(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxxy*gyzz+gxyy*&
gxzz)-2.0*gupxy(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxxy*gxzz+gxzy*gxxz)-2.0*gupxy(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxxy*&
gyzz+gxzy*gxyz)-gupxy(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(gxxy*gzzz+gxzy*gxzz)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*&
gxyx*gxxz-2.0*gupxy(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*gxzx*gxxz-2.0*gupxy(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*(gxxz*gyyx+gxyx&
*gxyz)-gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxxz*gyzx+gxzx*gxyz)-2.0*gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxxz*gyzx&
+gxyx*gxzz)-gupxy(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gxxz*gzzx+gxzx*gxzz)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gxyy*&
gxxz
      MapleGenVar4 = MapleGenVar3-2.0*gupxy(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*(gxxz*gyyy+gxyy*gxyz)-&
gupxy(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxxz*gyzy+gxzy*gxyz)-2.0*gupxy(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxxz*gyzy+gxyy*&
gxzz)-gupxy(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(gxxz*gzzy+gxzy*gxzz)-2.0*gupxy(i,j,k)*gupzz(i,j,k)*gupxx(i,j,k)*gxxz*gxxz&
-4.0*gupxy(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*gxxz*gxyz-2.0*gupxy(i,j,k)*gupzz(i,j,k)*gupyy(i,j,k)*(gxxz*gyyz+gxyz*gxyz)&
-3.0*gupxy(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*(gxxz*gyzz+gxyz*gxzz)-gupxy(i,j,k)*gupzz(i,j,k)*gupzz(i,j,k)*(gxxz*gzzz+gxzz*&
gxzz)-2.0*gupyy(i,j,k)*gupxx(i,j,k)*gupxx(i,j,k)*gxxx*gxyx-gupyy(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*(gxxx*gyzx+gxyx*gxzx)&
-2.0*gupyy(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gxyx*gxyx
      MapleGenVar1 = MapleGenVar4-4.0*gupxy(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gxxx*gxzx-4.0*gupxy(i,j,k)*&
gupyy(i,j,k)*gupxz(i,j,k)*gxxy*gxzy-4.0*gupxy(i,j,k)*gupyz(i,j,k)*gupxx(i,j,k)*gxxy*gxxz-2.0*gupxy(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*&
gxxy*gxzz-2.0*gupxy(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*gxzy*gxxz-4.0*gupxy(i,j,k)*gupzz(i,j,k)*gupxz(i,j,k)*gxxz*gxzz-2.0*&
gupyy(i,j,k)*gupyy(i,j,k)*gupxx(i,j,k)*gxyx*gyyx-2.0*gupyy(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*(gxyx*gyzx+gxzx*gyyx)-gupyy(i,j,k)*&
gupxx(i,j,k)*gupzz(i,j,k)*(gxyx*gzzx+gxzx*gyzx)-2.0*gupyy(i,j,k)*gupxy(i,j,k)*gupxx(i,j,k)*(gxxx*gxyy+gxyx*gxxy)-&
gupyy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxxx*gyzy+gxyx*gxzy)-4.0*gupyy(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*gxyx*gxyy-2.0*&
gupyy(i,j,k)*gupyy(i,j,k)*gupxy(i,j,k)*gxyx*gyyy+MapleGenVar2
      MapleGenVar4 = -2.0*gupyy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxyx*gxzy+gxzx*gxyy)-2.0*gupyy(i,j,k)*&
gupxy(i,j,k)*gupyz(i,j,k)*(gxyx*gyzy+gxzx*gyyy)-gupyy(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gxyx*gzzy+gxzx*gyzy)-2.0*&
gupyy(i,j,k)*gupxz(i,j,k)*gupxx(i,j,k)*(gxxx*gxyz+gxyx*gxxz)-gupyy(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxxx*gyzz+gxyx*gxzz)&
-2.0*gupyy(i,j,k)*gupyy(i,j,k)*gupxz(i,j,k)*gxyx*gyyz-2.0*gupyy(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxyx*gxzz+gxzx*gxyz)&
-2.0*gupyy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxyx*gyzz+gxzx*gyyz)-gupyy(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gxyx*gzzz+gxzx*&
gyzz)-gupyy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxxy*gyzx+gxzx*gxyy)-2.0*gupyy(i,j,k)*gupyy(i,j,k)*gupxy(i,j,k)*gxyy*gyyx&
-2.0*gupyy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxyy*gyzx+gxzy*gyyx)
      MapleGenVar3 = MapleGenVar4-gupyy(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gxyy*gzzx+gxzy*gyzx)-2.0*&
gupyy(i,j,k)*gupyy(i,j,k)*gupxx(i,j,k)*gxxy*gxyy-gupyy(i,j,k)*gupyy(i,j,k)*gupxz(i,j,k)*(gxxy*gyzy+gxyy*gxzy)-2.0*gupyy(i,j,k)*&
gupyy(i,j,k)*gupyz(i,j,k)*gxyy*gyzy-2.0*gupyy(i,j,k)*gupyy(i,j,k)*gupxz(i,j,k)*gxyy*gxzy-2.0*gupyy(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*(&
gxyy*gyzy+gxzy*gyyy)-gupyy(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*(gxyy*gzzy+gxzy*gyzy)-2.0*gupyy(i,j,k)*gupyz(i,j,k)*&
gupxx(i,j,k)*(gxxy*gxyz+gxyy*gxxz)-gupyy(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxxy*gyzz+gxyy*gxzz)-2.0*gupyy(i,j,k)*&
gupyy(i,j,k)*gupyz(i,j,k)*gxyy*gyyz-2.0*gupyy(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*gxyy*gyzz-2.0*gupyy(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(&
gxyy*gxzz+gxzy*gxyz)-2.0*gupyy(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxyy*gyzz+gxzy*gyyz)
      MapleGenVar4 = MapleGenVar3-gupyy(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(gxyy*gzzz+gxzy*gyzz)-gupyy(i,j,k)&
*gupxz(i,j,k)*gupxz(i,j,k)*(gxxz*gyzx+gxzx*gxyz)-2.0*gupyy(i,j,k)*gupyy(i,j,k)*gupxz(i,j,k)*gxyz*gyyx-2.0*gupyy(i,j,k)*&
gupxz(i,j,k)*gupyz(i,j,k)*(gxyz*gyzx+gxzz*gyyx)-gupyy(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gxyz*gzzx+gxzz*gyzx)-gupyy(i,j,k)&
*gupyz(i,j,k)*gupxz(i,j,k)*(gxxz*gyzy+gxzy*gxyz)-2.0*gupyy(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gxyz*gyyy-2.0*gupyy(i,j,k)*&
gupyz(i,j,k)*gupyz(i,j,k)*gxyz*gyzy-2.0*gupyy(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxyz*gyzy+gxzz*gyyy)-gupyy(i,j,k)*gupyz(i,j,k)*&
gupzz(i,j,k)*(gxyz*gzzy+gxzz*gyzy)-gupyy(i,j,k)*gupzz(i,j,k)*gupxz(i,j,k)*(gxxz*gyzz+gxyz*gxzz)-2.0*gupyy(i,j,k)*&
gupzz(i,j,k)*gupxy(i,j,k)*gxyz*gxyz
      MapleGenVar2 = MapleGenVar4-2.0*gupyy(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*gxyz*gyyz-2.0*gupyy(i,j,k)*&
gupzz(i,j,k)*gupyz(i,j,k)*(gxyz*gyzz+gxzz*gyyz)-gupyy(i,j,k)*gupzz(i,j,k)*gupzz(i,j,k)*(gxyz*gzzz+gxzz*gyzz)-2.0*&
gupyz(i,j,k)*gupxx(i,j,k)*gupxx(i,j,k)*gxxx*gxzx-gupyz(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*(gxxx*gzzx+gxzx*gxzx)-gupyz(i,j,k)*gupyz(i,j,k)&
*gupxx(i,j,k)*(gxyx*gzzx+gxzx*gyzx)-2.0*gupyz(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gxzx*gxzx-2.0*gupyz(i,j,k)*gupyz(i,j,k)*&
gupxx(i,j,k)*gxzx*gyzx-2.0*gupyz(i,j,k)*gupxy(i,j,k)*gupxx(i,j,k)*(gxxx*gxzy+gxzx*gxxy)-gupyz(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(&
gxxx*gzzy+gxzx*gxzy)-2.0*gupyz(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*(gxyx*gxzy+gxzx*gxyy)-gupyz(i,j,k)*gupyz(i,j,k)*&
gupxy(i,j,k)*(gxyx*gzzy+gxzx*gyzy)-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*gxzx*gyzy-2.0*gupyz(i,j,k)*gupxz(i,j,k)*&
gupxx(i,j,k)*(gxxx*gxzz+gxzx*gxxz)
      MapleGenVar4 = MapleGenVar2-gupyz(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxxx*gzzz+gxzx*gxzz)-2.0*&
gupyz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*(gxyx*gxzz+gxzx*gxyz)-gupyz(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxyx*gzzz+gxzx*gyzz)&
-4.0*gupyz(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*gxzx*gxzz-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*gxzx*gyzz-gupyz(i,j,k)*gupxy(i,j,k)*&
gupxz(i,j,k)*(gxxy*gzzx+gxzx*gxzy)-2.0*gupyy(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*gxyx*gyzx-2.0*gupyy(i,j,k)*gupxx(i,j,k)*&
gupxz(i,j,k)*gxyx*gxzx-2.0*gupyy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gxyx*gyzy-4.0*gupyy(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*gxyx*gxyz&
-2.0*gupyy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gxyx*gyzz-2.0*gupyy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gxyy*gyzx
      MapleGenVar3 = MapleGenVar4-4.0*gupyy(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*gxyy*gxyz-2.0*gupyy(i,j,k)*&
gupxz(i,j,k)*gupyz(i,j,k)*gxyz*gyzx-2.0*gupyy(i,j,k)*gupzz(i,j,k)*gupxx(i,j,k)*gxxz*gxyz-2.0*gupyy(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*&
gxyz*gyzz-2.0*gupyy(i,j,k)*gupzz(i,j,k)*gupxz(i,j,k)*gxyz*gxzz-2.0*gupyz(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gxyx*gxzx-2.0*&
gupyz(i,j,k)*gupxx(i,j,k)*gupzz(i,j,k)*gxzx*gzzx-4.0*gupyz(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxzx*gxzy-2.0*gupyz(i,j,k)*gupxy(i,j,k)*&
gupzz(i,j,k)*gxzx*gzzy-2.0*gupyz(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*gxzx*gzzz-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*gxzy*gyzx&
-gupyz(i,j,k)*gupyy(i,j,k)*gupxz(i,j,k)*(gxxy*gzzy+gxzy*gxzy)-gupyz(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*(gxyy*gzzy+gxzy*gyzy&
)
      MapleGenVar4 = MapleGenVar3-2.0*gupyz(i,j,k)*gupyy(i,j,k)*gupxz(i,j,k)*gxzy*gxzy-2.0*gupyz(i,j,k)*&
gupyz(i,j,k)*gupyy(i,j,k)*gxzy*gyzy-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxx(i,j,k)*(gxxy*gxzz+gxzy*gxxz)-gupyz(i,j,k)*gupyz(i,j,k)*&
gupxz(i,j,k)*(gxxy*gzzz+gxzy*gxzz)-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*(gxyy*gxzz+gxzy*gxyz)-gupyz(i,j,k)*&
gupyz(i,j,k)*gupyz(i,j,k)*(gxyy*gzzz+gxzy*gyzz)-4.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*gxzy*gxzz-2.0*gupyz(i,j,k)*&
gupyz(i,j,k)*gupzz(i,j,k)*gxzy*gzzz-gupyz(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxxz*gzzx+gxzx*gxzz)-gupyz(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)&
*(gxyz*gzzx+gxzz*gyzx)-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*gxzz*gyzx-gupyz(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxxz*&
gzzy+gxzy*gxzz)-gupyz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxyz*gzzy+gxzz*gyzy)
      CAZyx = Gamyx - (MapleGenVar4-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gxzz*gzzy-gupyz(i,j,k)*gupzz(i,j,k)*&
gupxz(i,j,k)*(gxxz*gzzz+gxzz*gxzz)-gupyz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(gxyz*gzzz+gxzz*gyzz)-2.0*gupyz(i,j,k)*&
gupzz(i,j,k)*gupxz(i,j,k)*gxzz*gxzz-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gxzz*gyzz-2.0*gupyz(i,j,k)*gupzz(i,j,k)*gupzz(i,j,k)*&
gxzz*gzzz-2.0*gupyz(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*gxzy*gzzx-2.0*gupyz(i,j,k)*gupyy(i,j,k)*gupxx(i,j,k)*gxxy*gxzy-2.0*&
gupyz(i,j,k)*gupyy(i,j,k)*gupxy(i,j,k)*gxyy*gxzy-2.0*gupyz(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*gxzy*gzzy-2.0*gupyz(i,j,k)*gupxz(i,j,k)*&
gupzz(i,j,k)*gxzz*gzzx-2.0*gupyz(i,j,k)*gupzz(i,j,k)*gupxx(i,j,k)*gxxz*gxzz-2.0*gupyz(i,j,k)*gupzz(i,j,k)*gupxy(i,j,k)*gxyz*gxzz&
+MapleGenVar1)
      MapleGenVar3 = gupzz(i,j,k)*gupxz(i,j,k)*gxxzz+gupyz(i,j,k)*gupxz(i,j,k)*gxxyz+gupyz(i,j,k)*gupyz(i,j,k)*gxyyz+2.0*&
gupyz(i,j,k)*gupzz(i,j,k)*gxzyz+gupyz(i,j,k)*gupxx(i,j,k)*gxxxy-2.0*gupzz(i,j,k)*gupzz(i,j,k)*gupxz(i,j,k)*gxzz*gxzz-2.0*gupyz(i,j,k)*&
gupyz(i,j,k)*gupyz(i,j,k)*gxyz*gyzy+gupzz(i,j,k)*gupxy(i,j,k)*gxyxz+gupyz(i,j,k)*gupxy(i,j,k)*gxxyy+gupzz(i,j,k)*gupyz(i,j,k)*gxyzz+&
gupxz(i,j,k)*gupxz(i,j,k)*gxxxz+2.0*gupxz(i,j,k)*gupzz(i,j,k)*gxzxz-2.0*gupxz(i,j,k)*gupxx(i,j,k)*gupxx(i,j,k)*gxxx*gxxx+gupzz(i,j,k)*&
gupxy(i,j,k)*gxxyz+gupxz(i,j,k)*gupyz(i,j,k)*gxyxz+gupxz(i,j,k)*gupxy(i,j,k)*gxxxy+gupzz(i,j,k)*gupzz(i,j,k)*gxzzz+gupyz(i,j,k)*gupyy(i,j,k)*&
gxyyy+gupxz(i,j,k)*gupyy(i,j,k)*gxyxy+gupyz(i,j,k)*gupyz(i,j,k)*gxzyy+gupxz(i,j,k)*gupxz(i,j,k)*gxzxx+gupzz(i,j,k)*gupxx(i,j,k)*gxxxz+&
2.0*gupxz(i,j,k)*gupyz(i,j,k)*gxzxy-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*gxzx*gxxz-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*&
gxyy*gyzz
      MapleGenVar4 = MapleGenVar3+gupzz(i,j,k)*gupyy(i,j,k)*gxyyz+gupxz(i,j,k)*gupxy(i,j,k)*gxyxx-2.0*gupzz(i,j,k)&
*gupzz(i,j,k)*gupzz(i,j,k)*gxzz*gzzz+gupxz(i,j,k)*gupxx(i,j,k)*gxxxx+gupyz(i,j,k)*gupxy(i,j,k)*gxyxy-2.0*gupxz(i,j,k)*gupxz(i,j,k)*&
gupxz(i,j,k)*gxxx*gxzz-2.0*gupyz(i,j,k)*gupyy(i,j,k)*gupxx(i,j,k)*gxxy*gxyy-4.0*gupxz(i,j,k)*gupxz(i,j,k)*gupxx(i,j,k)*gxxx*gxzx&
-gupxz(i,j,k)*gupxx(i,j,k)*gupyy(i,j,k)*(gxxx*gyyx+gxyx*gxyx)-3.0*gupxz(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*(gxxx*gyzx+gxyx*&
gxzx)-2.0*gupxz(i,j,k)*gupxx(i,j,k)*gupzz(i,j,k)*(gxxx*gzzx+gxzx*gxzx)-2.0*gupxz(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*gxxx*&
gxyy
      MapleGenVar2 = MapleGenVar4-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*gxxx*gxzy-2.0*gupxz(i,j,k)*&
gupxy(i,j,k)*gupxy(i,j,k)*(gxxx*gxyy+gxyx*gxxy)-gupxz(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*(gxxx*gyyy+gxyx*gxyy)-2.0*&
gupxz(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxxx*gyzy+gxyx*gxzy)-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*(gxxx*gxzy+gxzx*&
gxxy)-gupxz(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxxx*gyzy+gxzx*gxyy)-2.0*gupxz(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gxxx*gzzy+&
gxzx*gxzy)-4.0*gupxz(i,j,k)*gupxz(i,j,k)*gupxx(i,j,k)*gxxx*gxxz-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*gxxx*gxyz-2.0*&
gupxz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*(gxxx*gxyz+gxyx*gxxz)-gupxz(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*(gxxx*gyyz+gxyx*gxyz)&
-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxxx*gyzz+gxyx*gxzz)-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxxx*gxzz+&
gxzx*gxxz)-gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxxx*gyzz+gxzx*gxyz)
      MapleGenVar4 = MapleGenVar2-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gxxx*gzzz+gxzx*gxzz)&
-2.0*gupxz(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*gxyx*gxxy-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*gxzx*gxxy-gupxz(i,j,k)*gupxy(i,j,k)*&
gupyy(i,j,k)*(gxxy*gyyx+gxyx*gxyy)-2.0*gupxz(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxxy*gyzx+gxzx*gxyy)-gupxz(i,j,k)*&
gupxy(i,j,k)*gupyz(i,j,k)*(gxxy*gyzx+gxyx*gxzy)-2.0*gupxz(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gxxy*gzzx+gxzx*gxzy)&
-2.0*gupxz(i,j,k)*gupyy(i,j,k)*gupxx(i,j,k)*gxxy*gxxy-4.0*gupxz(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gxxy*gxzy-gupxz(i,j,k)*gupyy(i,j,k)*&
gupyy(i,j,k)*(gxxy*gyyy+gxyy*gxyy)-3.0*gupxz(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*(gxxy*gyzy+gxyy*gxzy)-2.0*&
gupxz(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*(gxxy*gzzy+gxzy*gxzy)
      MapleGenVar3 = MapleGenVar4-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gxxy*gxzz-2.0*gupxz(i,j,k)*&
gupyz(i,j,k)*gupxy(i,j,k)*(gxxy*gxyz+gxyy*gxxz)-gupxz(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*(gxxy*gyyz+gxyy*gxyz)-2.0*&
gupxz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxxy*gyzz+gxyy*gxzz)-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxxy*gxzz+gxzy*&
gxxz)-gupxz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxxy*gyzz+gxzy*gxyz)-2.0*gupxz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(gxxy*gzzz+&
gxzy*gxzz)-4.0*gupxz(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gxxx*gxyx-4.0*gupxz(i,j,k)*gupxy(i,j,k)*gupxx(i,j,k)*gxxx*gxxy-4.0*&
gupxz(i,j,k)*gupyy(i,j,k)*gupxy(i,j,k)*gxxy*gxyy-4.0*gupxz(i,j,k)*gupyz(i,j,k)*gupxx(i,j,k)*gxxy*gxxz-2.0*gupxz(i,j,k)*gupyz(i,j,k)*&
gupxy(i,j,k)*gxxy*gxyz-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*gxyx*gxxz
      MapleGenVar4 = MapleGenVar3-gupxz(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*(gxxz*gyyx+gxyx*gxyz)-2.0*&
gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxxz*gyzx+gxzx*gxyz)-gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxxz*gyzx+gxyx*gxzz)&
-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gxxz*gzzx+gxzx*gxzz)-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gxzy*gxxz-&
gupxz(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*(gxxz*gyyy+gxyy*gxyz)-2.0*gupxz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxxz*gyzy+gxzy*&
gxyz)-gupxz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxxz*gyzy+gxyy*gxzz)-2.0*gupxz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(gxxz*gzzy+&
gxzy*gxzz)-2.0*gupxz(i,j,k)*gupzz(i,j,k)*gupxx(i,j,k)*gxxz*gxxz-4.0*gupxz(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*gxxz*gxzz-&
gupxz(i,j,k)*gupzz(i,j,k)*gupyy(i,j,k)*(gxxz*gyyz+gxyz*gxyz)
      MapleGenVar1 = MapleGenVar4-3.0*gupxz(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*(gxxz*gyzz+gxyz*gxzz)&
-2.0*gupxz(i,j,k)*gupzz(i,j,k)*gupzz(i,j,k)*(gxxz*gzzz+gxzz*gxzz)-2.0*gupyz(i,j,k)*gupxx(i,j,k)*gupxx(i,j,k)*gxxx*gxyx-&
gupyz(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*(gxxx*gyyx+gxyx*gxyx)-2.0*gupyz(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gxyx*gxyx-2.0*&
gupyz(i,j,k)*gupyz(i,j,k)*gupxx(i,j,k)*gxyx*gyzx-gupyz(i,j,k)*gupyz(i,j,k)*gupxx(i,j,k)*(gxyx*gyzx+gxzx*gyyx)-2.0*gupyz(i,j,k)*&
gupxx(i,j,k)*gupzz(i,j,k)*(gxyx*gzzx+gxzx*gyzx)-2.0*gupyz(i,j,k)*gupxy(i,j,k)*gupxx(i,j,k)*(gxxx*gxyy+gxyx*gxxy)-&
gupyz(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*(gxxx*gyyy+gxyx*gxyy)-4.0*gupyz(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*gxyx*gxyy-2.0*&
gupyz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*gxyx*gyzy-2.0*gupyz(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxyx*gxzy+gxzx*gxyy)-gupyz(i,j,k)*&
gupyz(i,j,k)*gupxy(i,j,k)*(gxyx*gyzy+gxzx*gyyy)
      MapleGenVar4 = MapleGenVar1-2.0*gupyz(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gxyx*gzzy+gxzx*gyzy)&
-2.0*gupyz(i,j,k)*gupxz(i,j,k)*gupxx(i,j,k)*(gxxx*gxyz+gxyx*gxxz)-gupyz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*(gxxx*gyyz+gxyx*&
gxyz)-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*gxyx*gyzz-2.0*gupyz(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxyx*gxzz+gxzx*&
gxyz)-gupyz(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxyx*gyzz+gxzx*gyyz)-2.0*gupyz(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gxyx*gzzz+&
gxzx*gyzz)-gupyz(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*(gxxy*gyyx+gxyx*gxyy)-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*gxyy*&
gyzx-gupyz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*(gxyy*gyzx+gxzy*gyyx)-2.0*gupyz(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gxyy*gzzx+&
gxzy*gyzx)
      MapleGenVar3 = MapleGenVar4-gupyz(i,j,k)*gupyy(i,j,k)*gupxy(i,j,k)*(gxxy*gyyy+gxyy*gxyy)-2.0*&
gupyz(i,j,k)*gupyy(i,j,k)*gupxy(i,j,k)*gxyy*gxyy-2.0*gupyz(i,j,k)*gupyy(i,j,k)*gupyy(i,j,k)*gxyy*gyyy-2.0*gupyz(i,j,k)*gupyz(i,j,k)*&
gupyy(i,j,k)*gxyy*gyzy-gupyz(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*(gxyy*gyzy+gxzy*gyyy)-2.0*gupyz(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*(&
gxyy*gzzy+gxzy*gyzy)-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxx(i,j,k)*(gxxy*gxyz+gxyy*gxxz)-gupyz(i,j,k)*gupyz(i,j,k)*&
gupxy(i,j,k)*(gxxy*gyyz+gxyy*gxyz)-4.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*gxyy*gxyz-2.0*gupyz(i,j,k)*gupyz(i,j,k)*&
gupyy(i,j,k)*gxyy*gyyz-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxyy*gxzz+gxzy*gxyz)-gupyz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(&
gxyy*gyzz+gxzy*gyyz)-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(gxyy*gzzz+gxzy*gyzz)
      MapleGenVar4 = MapleGenVar3-gupyz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*(gxxz*gyyx+gxyx*gxyz)-2.0*&
gupyz(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*gxyz*gyzx-gupyz(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxyz*gyzx+gxzz*gyyx)-2.0*gupyz(i,j,k)*&
gupxz(i,j,k)*gupzz(i,j,k)*(gxyz*gzzx+gxzz*gyzx)-gupyz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*(gxxz*gyyy+gxyy*gxyz)-2.0*&
gupyz(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*gxyz*gyyy-gupyz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxyz*gyzy+gxzz*gyyy)-2.0*gupyz(i,j,k)*&
gupyz(i,j,k)*gupzz(i,j,k)*(gxyz*gzzy+gxzz*gyzy)-gupyz(i,j,k)*gupzz(i,j,k)*gupxy(i,j,k)*(gxxz*gyyz+gxyz*gxyz)-2.0*&
gupyz(i,j,k)*gupzz(i,j,k)*gupxy(i,j,k)*gxyz*gxyz-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gxyz*gyzz-gupyz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(&
gxyz*gyzz+gxzz*gyyz)
      MapleGenVar2 = MapleGenVar4-2.0*gupyz(i,j,k)*gupzz(i,j,k)*gupzz(i,j,k)*(gxyz*gzzz+gxzz*gyzz)&
-2.0*gupxz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*gxyy*gxxz-4.0*gupxz(i,j,k)*gupzz(i,j,k)*gupxy(i,j,k)*gxxz*gxyz-2.0*gupyz(i,j,k)*&
gupxx(i,j,k)*gupyy(i,j,k)*gxyx*gyyx-2.0*gupyz(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gxyx*gxzx-2.0*gupyz(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*&
gxyx*gyyy-4.0*gupyz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*gxyx*gxyz-2.0*gupyz(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gxyx*gyyz-2.0*&
gupyz(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*gxyy*gyyx-2.0*gupyz(i,j,k)*gupyy(i,j,k)*gupxz(i,j,k)*gxyy*gxzy-2.0*gupyz(i,j,k)*gupxz(i,j,k)*&
gupyy(i,j,k)*gxyz*gyyx-2.0*gupyz(i,j,k)*gupzz(i,j,k)*gupxx(i,j,k)*gxxz*gxyz-2.0*gupyz(i,j,k)*gupzz(i,j,k)*gupyy(i,j,k)*gxyz*gyyz&
-2.0*gupyz(i,j,k)*gupzz(i,j,k)*gupxz(i,j,k)*gxyz*gxzz
      MapleGenVar4 = MapleGenVar2-2.0*gupzz(i,j,k)*gupxx(i,j,k)*gupxx(i,j,k)*gxxx*gxzx-gupzz(i,j,k)*gupxx(i,j,k)*&
gupxy(i,j,k)*(gxxx*gyzx+gxyx*gxzx)-gupzz(i,j,k)*gupxx(i,j,k)*gupyy(i,j,k)*(gxyx*gyzx+gxzx*gyyx)-2.0*gupzz(i,j,k)*&
gupxx(i,j,k)*gupxz(i,j,k)*gxzx*gxzx-2.0*gupzz(i,j,k)*gupzz(i,j,k)*gupxx(i,j,k)*gxzx*gzzx-2.0*gupzz(i,j,k)*gupxy(i,j,k)*gupxx(i,j,k)*(&
gxxx*gxzy+gxzx*gxxy)-gupzz(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*(gxxx*gyzy+gxzx*gxyy)-2.0*gupzz(i,j,k)*gupxy(i,j,k)*&
gupxy(i,j,k)*(gxyx*gxzy+gxzx*gxyy)-gupzz(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*(gxyx*gyzy+gxzx*gyyy)-2.0*gupzz(i,j,k)*&
gupzz(i,j,k)*gupxy(i,j,k)*gxzx*gzzy-2.0*gupzz(i,j,k)*gupxz(i,j,k)*gupxx(i,j,k)*(gxxx*gxzz+gxzx*gxxz)-gupzz(i,j,k)*gupxz(i,j,k)*&
gupxy(i,j,k)*(gxxx*gyzz+gxzx*gxyz)
      MapleGenVar3 = MapleGenVar4-2.0*gupzz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*(gxyx*gxzz+gxzx*gxyz)-&
gupzz(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*(gxyx*gyzz+gxzx*gyyz)-4.0*gupzz(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*gxzx*gxzz-2.0*&
gupzz(i,j,k)*gupzz(i,j,k)*gupxz(i,j,k)*gxzx*gzzz-gupzz(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*(gxxy*gyzx+gxyx*gxzy)-gupzz(i,j,k)*gupxy(i,j,k)&
*gupyy(i,j,k)*(gxyy*gyzx+gxzy*gyyx)-2.0*gupzz(i,j,k)*gupzz(i,j,k)*gupxy(i,j,k)*gxzy*gzzx-gupzz(i,j,k)*gupyy(i,j,k)*gupxy(i,j,k)*&
(gxxy*gyzy+gxyy*gxzy)-gupzz(i,j,k)*gupyy(i,j,k)*gupyy(i,j,k)*(gxyy*gyzy+gxzy*gyyy)-2.0*gupzz(i,j,k)*gupyy(i,j,k)*&
gupxz(i,j,k)*gxzy*gxzy-2.0*gupzz(i,j,k)*gupzz(i,j,k)*gupyy(i,j,k)*gxzy*gzzy-2.0*gupzz(i,j,k)*gupyz(i,j,k)*gupxx(i,j,k)*(gxxy*&
gxzz+gxzy*gxxz)-gupzz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*(gxxy*gyzz+gxzy*gxyz)
      MapleGenVar4 = MapleGenVar3-2.0*gupzz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*(gxyy*gxzz+gxzy*gxyz)-&
gupzz(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*(gxyy*gyzz+gxzy*gyyz)-2.0*gupzz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*gxzy*gyzz-2.0*&
gupzz(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*gxzy*gzzz-gupzz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*(gxxz*gyzx+gxyx*gxzz)-gupzz(i,j,k)*gupxz(i,j,k)&
*gupyy(i,j,k)*(gxyz*gyzx+gxzz*gyyx)-2.0*gupzz(i,j,k)*gupzz(i,j,k)*gupxz(i,j,k)*gxzz*gzzx-gupzz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*&
(gxxz*gyzy+gxyy*gxzz)-gupzz(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*(gxyz*gyzy+gxzz*gyyy)-2.0*gupzz(i,j,k)*gupyz(i,j,k)*&
gupyz(i,j,k)*gxzz*gyzy-2.0*gupzz(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*gxzz*gzzy-2.0*gupzz(i,j,k)*gupzz(i,j,k)*gupxx(i,j,k)*gxxz*gxzz&
-gupzz(i,j,k)*gupzz(i,j,k)*gupxy(i,j,k)*(gxxz*gyzz+gxyz*gxzz)
      CAZzx = Gamzx -(MapleGenVar4-2.0*gupzz(i,j,k)*gupzz(i,j,k)*gupxy(i,j,k)*gxyz*gxzz-gupzz(i,j,k)*gupzz(i,j,k)*&
gupyy(i,j,k)*(gxyz*gyzz+gxzz*gyyz)-2.0*gupzz(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*gxzz*gyzz-2.0*gupzz(i,j,k)*gupxx(i,j,k)*&
gupxy(i,j,k)*gxyx*gxzx-2.0*gupzz(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*gxzx*gyzx-4.0*gupzz(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxzx*gxzy&
-2.0*gupzz(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gxzx*gyzy-2.0*gupzz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gxzx*gyzz-2.0*gupzz(i,j,k)*&
gupxy(i,j,k)*gupyz(i,j,k)*gxzy*gyzx-2.0*gupzz(i,j,k)*gupyy(i,j,k)*gupxx(i,j,k)*gxxy*gxzy-2.0*gupzz(i,j,k)*gupyy(i,j,k)*gupxy(i,j,k)*&
gxyy*gxzy-2.0*gupzz(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gxzy*gyzy-4.0*gupzz(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*gxzy*gxzz-2.0*&
gupzz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gxzz*gyzx)
      MapleGenVar3 = gupxz(i,j,k)*gupxz(i,j,k)*gyzxz+2.0*gupxy(i,j,k)*gupxz(i,j,k)*gxyyz+gupxx(i,j,k)*gupyz(i,j,k)*gyzxy&
-2.0*gupxy(i,j,k)*gupyy(i,j,k)*gupyy(i,j,k)*gyyy*gyyy+gupxx(i,j,k)*gupxx(i,j,k)*gxyxx+gupxy(i,j,k)*gupyz(i,j,k)*gyzyy-2.0*gupxz(i,j,k)*&
gupxz(i,j,k)*gupxz(i,j,k)*gxzx*gyzz-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*gxyx*gyyy+gupxy(i,j,k)*gupyy(i,j,k)*gyyyy+gupxz(i,j,k)*&
gupyz(i,j,k)*gyyzz+gupxy(i,j,k)*gupzz(i,j,k)*gyzyz+gupxx(i,j,k)*gupyz(i,j,k)*gyyxz+gupxy(i,j,k)*gupyz(i,j,k)*gyyyz-2.0*gupxx(i,j,k)*&
gupxx(i,j,k)*gupxx(i,j,k)*gxxx*gxyx+gupxy(i,j,k)*gupxy(i,j,k)*gxyyy+gupxz(i,j,k)*gupyy(i,j,k)*gyyyz+gupxz(i,j,k)*gupxy(i,j,k)*gyyxz+&
gupxx(i,j,k)*gupxz(i,j,k)*gyzxx+gupxz(i,j,k)*gupxz(i,j,k)*gxyzz+2.0*gupxx(i,j,k)*gupxy(i,j,k)*gxyxy+2.0*gupxx(i,j,k)*gupxz(i,j,k)*gxyxz&
+gupxz(i,j,k)*gupzz(i,j,k)*gyzzz-2.0*gupxx(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gxyx*gxyx+gupxy(i,j,k)*gupxy(i,j,k)*gyyxy+gupxz(i,j,k)*&
gupyz(i,j,k)*gyzyz
      MapleGenVar4 = MapleGenVar3-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*gxzz*gyzx-2.0*gupxy(i,j,k)*&
gupxy(i,j,k)*gupxy(i,j,k)*gxyy*gyyx+gupxx(i,j,k)*gupzz(i,j,k)*gyzxz-2.0*gupxx(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxyx*gyzz+gxyz*&
gyzx)+gupxx(i,j,k)*gupxy(i,j,k)*gyyxx-2.0*gupxy(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gyyy*gyzz+gyzy*gyyz)-4.0*gupxy(i,j,k)*&
gupxz(i,j,k)*gupyy(i,j,k)*gyyx*gyyz-2.0*gupxx(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxxz*gyzx+gxyx*gxzz)+gupxx(i,j,k)*gupyy(i,j,k)*&
gyyxy-gupxy(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(gyyy*gzzz+gyzy*gyzz)+gupxy(i,j,k)*gupxz(i,j,k)*gyzxy-4.0*gupxx(i,j,k)*gupxz(i,j,k)&
*gupxy(i,j,k)*gxyx*gxyz
      MapleGenVar2 = MapleGenVar4-2.0*gupxx(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gxyx*gxzx-2.0*gupxx(i,j,k)*&
gupxx(i,j,k)*gupxy(i,j,k)*(gxxx*gyyx+gxyx*gxyx)-2.0*gupxx(i,j,k)*gupxx(i,j,k)*gupyy(i,j,k)*gxyx*gyyx-gupxx(i,j,k)*gupxx(i,j,k)*&
gupyz(i,j,k)*(gxyx*gyzx+gxzx*gyyx)-2.0*gupxx(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*(gxxx*gyzx+gxyx*gxzx)-2.0*&
gupxx(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*gxyx*gyzx-gupxx(i,j,k)*gupxx(i,j,k)*gupzz(i,j,k)*(gxyx*gzzx+gxzx*gyzx)-2.0*gupxx(i,j,k)*&
gupxx(i,j,k)*gupxy(i,j,k)*gxyx*gxxy-4.0*gupxx(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*gxyx*gxyy-2.0*gupxx(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*(&
gxxy*gyyx+gxyx*gxyy)-2.0*gupxx(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*(gxyx*gyyy+gxyy*gyyx)-gupxx(i,j,k)*gupxy(i,j,k)*&
gupyz(i,j,k)*(gxyx*gyzy+gxzy*gyyx)-2.0*gupxx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxxy*gyzx+gxyx*gxzy)-2.0*&
gupxx(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxyx*gyzy+gxyy*gyzx)
      MapleGenVar4 = MapleGenVar2-gupxx(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gxyx*gzzy+gxzy*gyzx)-2.0*&
gupxx(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gxyx*gxxz-2.0*gupxx(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*gxyx*gxzz-2.0*gupxx(i,j,k)*gupxz(i,j,k)*&
gupxy(i,j,k)*(gxxz*gyyx+gxyx*gxyz)-2.0*gupxx(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*(gxyx*gyyz+gxyz*gyyx)-gupxx(i,j,k)*&
gupxz(i,j,k)*gupyz(i,j,k)*(gxyx*gyzz+gxzz*gyyx)-2.0*gupxx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxyx*gxzy-gupxx(i,j,k)*gupxz(i,j,k)*&
gupzz(i,j,k)*(gxyx*gzzz+gxzz*gyzx)-2.0*gupxx(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gxxx*gxyy-2.0*gupxx(i,j,k)*gupxy(i,j,k)*&
gupxy(i,j,k)*(gxxx*gyyy+gxyx*gxyy)-gupxx(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxyy*gyzx+gxzx*gyyy)-2.0*gupxx(i,j,k)*&
gupxy(i,j,k)*gupxz(i,j,k)*(gxxx*gyzy+gxzx*gxyy)
      MapleGenVar3 = MapleGenVar4-gupxx(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gxyy*gzzx+gxzx*gyzy)-2.0*&
gupxx(i,j,k)*gupxx(i,j,k)*gupyy(i,j,k)*gxxy*gxyy-2.0*gupxx(i,j,k)*gupyy(i,j,k)*gupxy(i,j,k)*gxyy*gxyy-2.0*gupxx(i,j,k)*gupyy(i,j,k)*&
gupxy(i,j,k)*(gxxy*gyyy+gxyy*gxyy)-2.0*gupxx(i,j,k)*gupyy(i,j,k)*gupyy(i,j,k)*gxyy*gyyy-gupxx(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*(&
gxyy*gyzy+gxzy*gyyy)-2.0*gupxx(i,j,k)*gupyy(i,j,k)*gupxz(i,j,k)*(gxxy*gyzy+gxyy*gxzy)-gupxx(i,j,k)*gupyy(i,j,k)*&
gupzz(i,j,k)*(gxyy*gzzy+gxzy*gyzy)-2.0*gupxx(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*gxyy*gxxz-2.0*gupxx(i,j,k)*gupyz(i,j,k)*&
gupxy(i,j,k)*(gxxz*gyyy+gxyy*gxyz)-2.0*gupxx(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*(gxyy*gyyz+gxyz*gyyy)-gupxx(i,j,k)*&
gupyz(i,j,k)*gupyz(i,j,k)*(gxyy*gyzz+gxzz*gyyy)-2.0*gupxx(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxxz*gyzy+gxyy*gxzz)
      MapleGenVar4 = MapleGenVar3-2.0*gupxx(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxyy*gyzz+gxyz*gyzy)-&
gupxx(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(gxyy*gzzz+gxzz*gyzy)-2.0*gupxx(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gxxx*gxyz-2.0*&
gupxx(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*gxzx*gxyz-2.0*gupxx(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*(gxxx*gyyz+gxyx*gxyz)-gupxx(i,j,k)*&
gupxz(i,j,k)*gupyz(i,j,k)*(gxyz*gyzx+gxzx*gyyz)-2.0*gupxx(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxxx*gyzz+gxzx*gxyz)-&
gupxx(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gxyz*gzzx+gxzx*gyzz)-2.0*gupxx(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*gxxy*gxyz-2.0*&
gupxx(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*(gxxy*gyyz+gxyy*gxyz)-gupxx(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxyz*gyzy+gxzy*gyyz)&
-2.0*gupxx(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxxy*gyzz+gxzy*gxyz)
      MapleGenVar1 = MapleGenVar4-gupxx(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(gxyz*gzzy+gxzy*gyzz)-2.0*&
gupxx(i,j,k)*gupxx(i,j,k)*gupzz(i,j,k)*gxxz*gxyz-2.0*gupxx(i,j,k)*gupzz(i,j,k)*gupxy(i,j,k)*gxyz*gxyz-2.0*gupxx(i,j,k)*gupzz(i,j,k)*&
gupxy(i,j,k)*(gxxz*gyyz+gxyz*gxyz)-gupxx(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*(gxyz*gyzz+gxzz*gyyz)-2.0*gupxx(i,j,k)*&
gupzz(i,j,k)*gupxz(i,j,k)*(gxxz*gyzz+gxyz*gxzz)-gupxx(i,j,k)*gupzz(i,j,k)*gupzz(i,j,k)*(gxyz*gzzz+gxzz*gyzz)-4.0*&
gupxy(i,j,k)*gupxy(i,j,k)*gupxx(i,j,k)*gxyx*gyyx-3.0*gupxy(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*(gxyx*gyzx+gxzx*gyyx)-2.0*&
gupxy(i,j,k)*gupxx(i,j,k)*gupyy(i,j,k)*gyyx*gyyx-gupxy(i,j,k)*gupxx(i,j,k)*gupzz(i,j,k)*(gyyx*gzzx+gyzx*gyzx)-2.0*gupxy(i,j,k)*&
gupxy(i,j,k)*gupxy(i,j,k)*(gxyx*gyyy+gxyy*gyyx)-gupxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxyx*gyzy+gxzy*gyyx)-4.0*&
gupxy(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*gyyx*gyyy
      MapleGenVar4 = MapleGenVar1-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gyyx*gyzy-2.0*gupxy(i,j,k)*&
gupxy(i,j,k)*gupxz(i,j,k)*(gxyy*gyzx+gxzy*gyyx)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gyyx*gyzy+gyzx*gyyy)-&
gupxy(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gyyx*gzzy+gyzx*gyzy)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxyx*gyyz+gxyz*&
gyyx)-gupxy(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxyx*gyzz+gxzz*gyyx)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxyz*gyyx&
-2.0*gupxy(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxyz*gyzx+gxzz*gyyx)-2.0*gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gyyx*gyzz+&
gyzx*gyyz)-gupxy(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gyyx*gzzz+gyzx*gyzz)-gupxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxyy*gyzx&
+gxzx*gyyy)
      MapleGenVar3 = MapleGenVar4-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gyzx*gyyy-2.0*gupxy(i,j,k)*&
gupxy(i,j,k)*gupxz(i,j,k)*(gxyx*gyzy+gxzx*gyyy)-gupxy(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gyyy*gzzx+gyzx*gyzy)-4.0*&
gupxy(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*gxyy*gyyy-3.0*gupxy(i,j,k)*gupyy(i,j,k)*gupxz(i,j,k)*(gxyy*gyzy+gxzy*gyyy)-gupxy(i,j,k)*&
gupyy(i,j,k)*gupzz(i,j,k)*(gyyy*gzzy+gyzy*gyzy)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxyy*gyyz+gxyz*gyyy)-&
gupxy(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxyy*gyzz+gxzz*gyyy)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gxyz*gyyy-2.0*&
gupxy(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*gyyy*gyzz-2.0*gupxy(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxyz*gyzy+gxzz*gyyy)-2.0*&
gupxx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxzx*gxyy-2.0*gupxx(i,j,k)*gupyy(i,j,k)*gupxz(i,j,k)*gxyy*gxzy
      MapleGenVar4 = MapleGenVar3-2.0*gupxx(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gxyy*gyzy-4.0*gupxx(i,j,k)*&
gupyz(i,j,k)*gupxy(i,j,k)*gxyy*gxyz-2.0*gupxx(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*gxyy*gxzz-2.0*gupxx(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*&
gxzy*gxyz-2.0*gupxx(i,j,k)*gupzz(i,j,k)*gupxz(i,j,k)*gxyz*gxzz-2.0*gupxx(i,j,k)*gupzz(i,j,k)*gupyy(i,j,k)*gxyz*gyyz-2.0*&
gupxx(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*gxyz*gyzz-4.0*gupxy(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*gyyx*gyzx-2.0*gupxy(i,j,k)*gupxz(i,j,k)*&
gupyz(i,j,k)*gyyx*gyzz-4.0*gupxy(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gyyy*gyzy-4.0*gupxy(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*gyyy*gyyz&
-gupxy(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxyz*gyzx+gxzx*gyyz)
      MapleGenVar2 = MapleGenVar4-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxyx*gyyz-2.0*gupxy(i,j,k)*&
gupxz(i,j,k)*gupxz(i,j,k)*(gxyx*gyzz+gxzx*gyyz)-gupxy(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gyyz*gzzx+gyzx*gyzz)-gupxy(i,j,k)&
*gupyz(i,j,k)*gupxz(i,j,k)*(gxyz*gyzy+gxzy*gyyz)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gxyy*gyyz-2.0*gupxy(i,j,k)*&
gupyz(i,j,k)*gupyz(i,j,k)*gyzy*gyyz-2.0*gupxy(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxyy*gyzz+gxzy*gyyz)-gupxy(i,j,k)*gupyz(i,j,k)*&
gupzz(i,j,k)*(gyyz*gzzy+gyzy*gyzz)-4.0*gupxy(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*gxyz*gyyz-3.0*gupxy(i,j,k)*gupzz(i,j,k)*&
gupxz(i,j,k)*(gxyz*gyzz+gxzz*gyyz)-2.0*gupxy(i,j,k)*gupzz(i,j,k)*gupyy(i,j,k)*gyyz*gyyz-gupxy(i,j,k)*gupzz(i,j,k)*gupzz(i,j,k)*(&
gyyz*gzzz+gyzz*gyzz)-gupxz(i,j,k)*gupxz(i,j,k)*gupxx(i,j,k)*(gxyx*gzzx+gxzx*gyzx)-gupxz(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*&
(gyyx*gzzx+gyzx*gyzx)
      MapleGenVar4 = MapleGenVar2-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupxx(i,j,k)*gxzx*gyzx-2.0*gupxz(i,j,k)*&
gupxx(i,j,k)*gupyz(i,j,k)*gyzx*gyzx-2.0*gupxz(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*(gxyx*gyzy+gxyy*gyzx)-gupxz(i,j,k)*gupxz(i,j,k)*&
gupxy(i,j,k)*(gxyx*gzzy+gxzy*gyzx)-2.0*gupxz(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*(gyyx*gyzy+gyzx*gyyy)-gupxz(i,j,k)*&
gupxy(i,j,k)*gupyz(i,j,k)*(gyyx*gzzy+gyzx*gyzy)-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*gxzy*gyzx-2.0*gupxz(i,j,k)*&
gupxz(i,j,k)*gupxy(i,j,k)*(gxyx*gyzz+gxyz*gyzx)-gupxz(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxyx*gzzz+gxzz*gyzx)-2.0*&
gupxz(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*(gyyx*gyzz+gyzx*gyyz)-gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gyyx*gzzz+gyzx*gyzz)&
-4.0*gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gyzx*gyzz
      MapleGenVar3 = MapleGenVar4-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*gyzx*gzzz-gupxz(i,j,k)*gupxz(i,j,k)*&
gupxy(i,j,k)*(gxyy*gzzx+gxzx*gyzy)-gupxz(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gyyy*gzzx+gyzx*gyzy)-2.0*gupxz(i,j,k)*&
gupxz(i,j,k)*gupxy(i,j,k)*gxzx*gyzy-gupxz(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*(gxyy*gzzy+gxzy*gyzy)-2.0*gupxz(i,j,k)*gupyy(i,j,k)*&
gupyy(i,j,k)*gyyy*gyzy-gupxz(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*(gyyy*gzzy+gyzy*gyzy)-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*&
gxzy*gyzy-2.0*gupxz(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gyzy*gyzy-2.0*gupxz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*(gxyy*gyzz+gxyz&
*gyzy)-gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxyy*gzzz+gxzz*gyzy)-2.0*gupxz(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*(gyyy*gyzz&
+gyzy*gyyz)-gupxz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gyyy*gzzz+gyzy*gyzz)
      MapleGenVar4 = MapleGenVar3-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gxzz*gyzy-4.0*gupxz(i,j,k)*&
gupyz(i,j,k)*gupyz(i,j,k)*gyzy*gyzz-gupxz(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxyz*gzzx+gxzx*gyzz)-gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)&
*(gyyz*gzzx+gyzx*gyzz)-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*gyzz*gzzx-gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxyz*&
gzzy+gxzy*gyzz)-gupxz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gyyz*gzzy+gyzy*gyzz)-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*&
gxzy*gyzz-gupxz(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gxyz*gzzz+gxzz*gyzz)-gupxz(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*(gyyz*gzzz+&
gyzz*gyzz)-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*gxzz*gyzz-2.0*gupxz(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*gyzz*gyzz-2.0*&
gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gyzx*gyyz
      CAZxy = Gamxy - (MapleGenVar4-4.0*gupxy(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*gyyz*gyzz-2.0*gupxz(i,j,k)*&
gupxx(i,j,k)*gupxy(i,j,k)*gxyx*gyzx-2.0*gupxz(i,j,k)*gupxx(i,j,k)*gupyy(i,j,k)*gyyx*gyzx-2.0*gupxz(i,j,k)*gupxx(i,j,k)*gupzz(i,j,k)*&
gyzx*gzzx-4.0*gupxz(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gyzx*gyzy-2.0*gupxz(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*gyzx*gzzy-2.0*&
gupxz(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*gyzy*gzzx-2.0*gupxz(i,j,k)*gupyy(i,j,k)*gupxy(i,j,k)*gxyy*gyzy-2.0*gupxz(i,j,k)*gupyy(i,j,k)*&
gupzz(i,j,k)*gyzy*gzzy-2.0*gupxz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gyzy*gzzz-2.0*gupxz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gyzz*gzzy&
-2.0*gupxz(i,j,k)*gupzz(i,j,k)*gupxy(i,j,k)*gxyz*gyzz-2.0*gupxz(i,j,k)*gupzz(i,j,k)*gupyy(i,j,k)*gyyz*gyzz-2.0*gupxz(i,j,k)*&
gupzz(i,j,k)*gupzz(i,j,k)*gyzz*gzzz)
      MapleGenVar3 = gupxy(i,j,k)*gupxy(i,j,k)*gxyxy-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*gxyy*gxyy+gupxx(i,j,k)*&
gupyz(i,j,k)*gxyxz-2.0*gupyy(i,j,k)*gupyy(i,j,k)*gupxx(i,j,k)*gyyx*gyyx-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gyzz*gyzz+&
gupxx(i,j,k)*gupyy(i,j,k)*gxyxy-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*gyzy*gyzy+gupxy(i,j,k)*gupyz(i,j,k)*gyzxy+gupxz(i,j,k)*gupyy(i,j,k)&
*gxyyz+gupyz(i,j,k)*gupyz(i,j,k)*gyzyz+2.0*gupyy(i,j,k)*gupyz(i,j,k)*gyyyz+gupyy(i,j,k)*gupzz(i,j,k)*gyzyz+2.0*gupxy(i,j,k)*&
gupyy(i,j,k)*gyyxy+gupxz(i,j,k)*gupxy(i,j,k)*gxyxz+gupxz(i,j,k)*gupyz(i,j,k)*gxyzz-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxx(i,j,k)*gxyx*gxyx&
+gupxy(i,j,k)*gupxz(i,j,k)*gyzxx+gupyy(i,j,k)*gupxz(i,j,k)*gyzxy-4.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*gxyx*gxyy+gupyz(i,j,k)*&
gupxz(i,j,k)*gyzxz+gupyy(i,j,k)*gupyz(i,j,k)*gyzyy-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxx(i,j,k)*gyzx*gyzx+gupyz(i,j,k)*gupyz(i,j,k)*gyyzz&
-4.0*gupyz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*gyzy*gyzz+gupyz(i,j,k)*gupzz(i,j,k)*gyzzz+gupxy(i,j,k)*gupyz(i,j,k)*gxyyz
      MapleGenVar4 = MapleGenVar3+gupxy(i,j,k)*gupzz(i,j,k)*gyzxz+gupxx(i,j,k)*gupxy(i,j,k)*gxyxx-2.0*gupxy(i,j,k)&
*gupxy(i,j,k)*gupzz(i,j,k)*gxyz*gxyz+gupxy(i,j,k)*gupyy(i,j,k)*gxyyy-2.0*gupyy(i,j,k)*gupyy(i,j,k)*gupyy(i,j,k)*gyyy*gyyy+gupxy(i,j,k)*&
gupxy(i,j,k)*gyyxx+gupyy(i,j,k)*gupyy(i,j,k)*gyyyy+2.0*gupxy(i,j,k)*gupyz(i,j,k)*gyyxz-2.0*gupyy(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*gyyz*&
gyyz-2.0*gupxx(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*gxyx*gxxz-2.0*gupxx(i,j,k)*gupzz(i,j,k)*gupxy(i,j,k)*gxxz*gxyz-2.0*gupxx(i,j,k)*&
gupxx(i,j,k)*gupxy(i,j,k)*gxxx*gxyx
      MapleGenVar2 = MapleGenVar4-gupxx(i,j,k)*gupxx(i,j,k)*gupyy(i,j,k)*(gxxx*gyyx+gxyx*gxyx)-gupxx(i,j,k)&
*gupxx(i,j,k)*gupyz(i,j,k)*(gxxx*gyzx+gxyx*gxzx)-2.0*gupxx(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*gxxx*gxyy-gupxx(i,j,k)*gupxy(i,j,k)*&
gupyy(i,j,k)*(gxxx*gyyy+gxyx*gxyy)-gupxx(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxxx*gyzy+gxzx*gxyy)-gupxx(i,j,k)*gupxz(i,j,k)&
*gupyy(i,j,k)*(gxxx*gyyz+gxyx*gxyz)-gupxx(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxxx*gyzz+gxzx*gxyz)-2.0*gupxx(i,j,k)*&
gupxy(i,j,k)*gupxy(i,j,k)*gxyx*gxxy-gupxx(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*(gxxy*gyyx+gxyx*gxyy)-gupxx(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)&
*(gxxy*gyzx+gxyx*gxzy)-gupxx(i,j,k)*gupyy(i,j,k)*gupyy(i,j,k)*(gxxy*gyyy+gxyy*gxyy)-gupxx(i,j,k)*gupyy(i,j,k)*&
gupyz(i,j,k)*(gxxy*gyzy+gxyy*gxzy)-gupxx(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*(gxxy*gyyz+gxyy*gxyz)-gupxx(i,j,k)*gupyz(i,j,k)&
*gupyz(i,j,k)*(gxxy*gyzz+gxzy*gxyz)
      MapleGenVar4 = MapleGenVar2-gupxx(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*(gxxz*gyyx+gxyx*gxyz)-gupxx(i,j,k)&
*gupxz(i,j,k)*gupyz(i,j,k)*(gxxz*gyzx+gxyx*gxzz)-gupxx(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*(gxxz*gyyy+gxyy*gxyz)-&
gupxx(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxxz*gyzy+gxyy*gxzz)-gupxx(i,j,k)*gupzz(i,j,k)*gupyy(i,j,k)*(gxxz*gyyz+gxyz*gxyz)&
-gupxx(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*(gxxz*gyzz+gxyz*gxzz)-gupxy(i,j,k)*gupxy(i,j,k)*gupxx(i,j,k)*(gxxx*gyyx+gxyx*gxyx&
)-gupxy(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*(gxxx*gyzx+gxyx*gxzx)-2.0*gupxy(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*(gxyx*gyzx+gxzx&
*gyyx)-gupxy(i,j,k)*gupxx(i,j,k)*gupzz(i,j,k)*(gxyx*gzzx+gxzx*gyzx)-gupxy(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*(gxxx*gyyy+&
gxyx*gxyy)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*gxyx*gyyy
      MapleGenVar3 = MapleGenVar4-gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxyx*gyzy+gxzx*gyyy)-gupxy(i,j,k)&
*gupxy(i,j,k)*gupxz(i,j,k)*(gxxx*gyyz+gxyx*gxyz)-4.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxyx*gxyz-2.0*gupxx(i,j,k)*&
gupxz(i,j,k)*gupxy(i,j,k)*gxxx*gxyz-2.0*gupxx(i,j,k)*gupyy(i,j,k)*gupxy(i,j,k)*gxxy*gxyy-2.0*gupxx(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*&
gxxy*gxyz-2.0*gupxx(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*gxyy*gxxz-6.0*gupxy(i,j,k)*gupxx(i,j,k)*gupyy(i,j,k)*gxyx*gyyx-4.0*&
gupxy(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*gxyx*gyzx-2.0*gupxy(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gxyx*gxzx-2.0*gupxy(i,j,k)*gupxz(i,j,k)*&
gupyy(i,j,k)*gxyx*gyyz-gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxyx*gyzz+gxzx*gyyz)-gupxy(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*(gxxy&
*gyyx+gxyx*gxyy)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*gxyy*gyyx
      MapleGenVar4 = MapleGenVar3-gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxyy*gyzx+gxzy*gyyx)-gupxy(i,j,k)&
*gupxy(i,j,k)*gupyy(i,j,k)*(gxxy*gyyy+gxyy*gxyy)-gupxy(i,j,k)*gupyy(i,j,k)*gupxz(i,j,k)*(gxxy*gyzy+gxyy*gxzy)-6.0*&
gupxy(i,j,k)*gupyy(i,j,k)*gupyy(i,j,k)*gxyy*gyyy-2.0*gupxy(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*(gxyy*gyzy+gxzy*gyyy)-gupxy(i,j,k)*&
gupyy(i,j,k)*gupzz(i,j,k)*(gxyy*gzzy+gxzy*gyzy)-gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxxy*gyyz+gxyy*gxyz)-4.0*&
gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gxyy*gxyz-gupxy(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxyy*gyzz+gxzy*gyyz)-gupxy(i,j,k)*gupxy(i,j,k)&
*gupxz(i,j,k)*(gxxz*gyyx+gxyx*gxyz)-gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxyz*gyzx+gxzz*gyyx)-gupxy(i,j,k)*&
gupxy(i,j,k)*gupyz(i,j,k)*(gxxz*gyyy+gxyy*gxyz)-gupxy(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxyz*gyzy+gxzz*gyyy)
      MapleGenVar1 = MapleGenVar4-gupxy(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gxxz*gyyz+gxyz*gxyz)-gupxy(i,j,k)&
*gupzz(i,j,k)*gupxz(i,j,k)*(gxxz*gyzz+gxyz*gxzz)-2.0*gupxy(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*(gxyz*gyzz+gxzz*gyyz)-&
gupxy(i,j,k)*gupzz(i,j,k)*gupzz(i,j,k)*(gxyz*gzzz+gxzz*gyzz)-2.0*gupxz(i,j,k)*gupxx(i,j,k)*gupyy(i,j,k)*(gxyx*gyzx+gxzx*&
gyyx)-gupxz(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*(gxyx*gzzx+gxzx*gyzx)-gupxz(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*(gxxx*gyzy+gxzx&
*gxyy)-gupxz(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*(gxyx*gyzy+gxzx*gyyy)-gupxz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*(gxxx*gyzz+&
gxzx*gxyz)-gupxz(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*(gxyx*gyzz+gxzx*gyyz)-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gxzx*&
gyzz-gupxz(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*(gxxy*gyzx+gxyx*gxzy)-gupxz(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*(gxyy*gyzx+gxzy*&
gyyx)-2.0*gupxz(i,j,k)*gupyy(i,j,k)*gupyy(i,j,k)*(gxyy*gyzy+gxzy*gyyy)
      MapleGenVar4 = MapleGenVar1-gupxz(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*(gxyy*gzzy+gxzy*gyzy)-gupxz(i,j,k)&
*gupyz(i,j,k)*gupxy(i,j,k)*(gxxy*gyzz+gxzy*gxyz)-gupxz(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*(gxyy*gyzz+gxzy*gyyz)-2.0*&
gupxz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*gxzy*gyzz-gupxz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*(gxxz*gyzx+gxyx*gxzz)-gupxz(i,j,k)*gupxz(i,j,k)&
*gupyy(i,j,k)*(gxyz*gyzx+gxzz*gyyx)-4.0*gupxy(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gxyy*gyzy-2.0*gupxy(i,j,k)*gupyy(i,j,k)*&
gupxz(i,j,k)*gxyy*gxzy-2.0*gupxy(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*gxyy*gyyz-2.0*gupxy(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gxyz*gyyx&
-2.0*gupxy(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*gxyz*gyyy-6.0*gupxy(i,j,k)*gupzz(i,j,k)*gupyy(i,j,k)*gxyz*gyyz
      MapleGenVar3 = MapleGenVar4-4.0*gupxy(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*gxyz*gyzz-2.0*gupxy(i,j,k)*&
gupzz(i,j,k)*gupxz(i,j,k)*gxyz*gxzz-2.0*gupxz(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*gxzx*gyzx-2.0*gupxz(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*&
gxzx*gyzy-2.0*gupxz(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gxzy*gyzx-2.0*gupxz(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gxzy*gyzy-2.0*&
gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gxzz*gyzx-gupxz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*(gxxz*gyzy+gxyy*gxzz)-gupxz(i,j,k)*gupyz(i,j,k)&
*gupyy(i,j,k)*(gxyz*gyzy+gxzz*gyyy)-2.0*gupxz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*gxzz*gyzy-2.0*gupxz(i,j,k)*gupzz(i,j,k)*&
gupyy(i,j,k)*(gxyz*gyzz+gxzz*gyyz)-gupxz(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*(gxyz*gzzz+gxzz*gyzz)-2.0*gupxz(i,j,k)*&
gupzz(i,j,k)*gupyz(i,j,k)*gxzz*gyzz
      MapleGenVar4 = MapleGenVar3-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxyx*gxzy-4.0*gupxy(i,j,k)*&
gupxy(i,j,k)*gupyy(i,j,k)*(gxyx*gyyy+gxyy*gyyx)-gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxyx*gyzy+gxzy*gyyx)-4.0*&
gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxyx*gyzy+gxyy*gyzx)-gupxy(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gxyx*gzzy+gxzy*gyzx)&
-2.0*gupxy(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*gxyx*gxzz-4.0*gupxy(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*(gxyx*gyyz+gxyz*gyyx)-&
gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxyx*gyzz+gxzz*gyyx)-4.0*gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxyx*gyzz+gxyz*&
gyzx)-gupxy(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gxyx*gzzz+gxzz*gyzx)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxzx*gxyy-&
gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxyy*gyzx+gxzx*gyyy)-gupxy(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gxyy*gzzx+gxzx*gyzy)

      MapleGenVar2 = MapleGenVar4-4.0*gupxy(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*(gxyy*gyyz+gxyz*gyyy)-&
gupxy(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxyy*gyzz+gxzz*gyyy)-4.0*gupxy(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxyy*gyzz+gxyz*&
gyzy)-gupxy(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(gxyy*gzzz+gxzz*gyzy)-2.0*gupxy(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*gxzx*gxyz-&
gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxyz*gyzx+gxzx*gyyz)-gupxy(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gxyz*gzzx+gxzx*gyzz)&
-gupxy(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxyz*gyzy+gxzy*gyyz)-gupxy(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(gxyz*gzzy+gxzy*gyzz&
)-gupyy(i,j,k)*gupxx(i,j,k)*gupzz(i,j,k)*(gyyx*gzzx+gyzx*gyzx)-gupyy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxyx*gyzy+gxzy*&
gyyx)-4.0*gupyy(i,j,k)*gupyy(i,j,k)*gupxy(i,j,k)*gyyx*gyyy-4.0*gupyy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gyyx*gyzy+gyzx*&
gyyy)-gupyy(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gyyx*gzzy+gyzx*gyzy)
      MapleGenVar4 = MapleGenVar2-gupyy(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxyx*gyzz+gxzz*gyyx)-4.0*&
gupyy(i,j,k)*gupyy(i,j,k)*gupxz(i,j,k)*gyyx*gyyz-4.0*gupyy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gyyx*gyzz+gyzx*gyyz)-gupyy(i,j,k)*&
gupxz(i,j,k)*gupzz(i,j,k)*(gyyx*gzzz+gyzx*gyzz)-gupyy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxyy*gyzx+gxzx*gyyy)-2.0*&
gupxy(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*gxyy*gxzz-2.0*gupxy(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*gxzy*gxyz-6.0*gupyy(i,j,k)*gupxx(i,j,k)*&
gupyz(i,j,k)*gyyx*gyzx-2.0*gupyy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gyyx*gyzy-2.0*gupyy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gyyx*gyzz&
-gupyy(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gyyy*gzzx+gyzx*gyzy)-6.0*gupyy(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gyyy*gyzy
      MapleGenVar3 = MapleGenVar4-gupyy(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*(gyyy*gzzy+gyzy*gyzy)-gupyy(i,j,k)&
*gupyz(i,j,k)*gupxz(i,j,k)*(gxyy*gyzz+gxzz*gyyy)-4.0*gupyy(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gyyy*gyyz-2.0*gupyy(i,j,k)*&
gupyz(i,j,k)*gupyz(i,j,k)*gyyy*gyzz-4.0*gupyy(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gyyy*gyzz+gyzy*gyyz)-gupyy(i,j,k)*gupyz(i,j,k)*&
gupzz(i,j,k)*(gyyy*gzzz+gyzy*gyzz)-gupyy(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxyz*gyzx+gxzx*gyyz)-gupyy(i,j,k)*gupxz(i,j,k)&
*gupzz(i,j,k)*(gyyz*gzzx+gyzx*gyzz)-gupyy(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxyz*gyzy+gxzy*gyyz)-2.0*gupyy(i,j,k)*&
gupyz(i,j,k)*gupyz(i,j,k)*gyzy*gyyz-gupyy(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(gyyz*gzzy+gyzy*gyzz)-gupyy(i,j,k)*gupzz(i,j,k)*gupzz(i,j,k)&
*(gyyz*gzzz+gyzz*gyzz)-gupyz(i,j,k)*gupyz(i,j,k)*gupxx(i,j,k)*(gyyx*gzzx+gyzx*gyzx)-gupyz(i,j,k)*gupxy(i,j,k)*&
gupxz(i,j,k)*(gxyx*gzzy+gxzy*gyzx)
      MapleGenVar4 = MapleGenVar3-gupyz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*(gyyx*gzzy+gyzx*gyzy)-4.0*&
gupyz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*gyzx*gyzy-gupyz(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxyx*gzzz+gxzz*gyzx)-gupyz(i,j,k)*gupyz(i,j,k)&
*gupxz(i,j,k)*(gyyx*gzzz+gyzx*gyzz)-4.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*gyzx*gyzz-gupyz(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*&
(gxyy*gzzx+gxzx*gyzy)-gupyz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*(gyyy*gzzx+gyzx*gyzy)-gupyz(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)&
*(gyyy*gzzy+gyzy*gyzy)-gupyz(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxyy*gzzz+gxzz*gyzy)-gupyz(i,j,k)*gupyz(i,j,k)*&
gupyz(i,j,k)*(gyyy*gzzz+gyzy*gyzz)-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gyzy*gzzz-gupyz(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(&
gxyz*gzzx+gxzx*gyzz)-gupyz(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gyyz*gzzx+gyzx*gyzz)
      CAZyy = Gamyy - (MapleGenVar4-2.0*gupyy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gyzx*gyyy-2.0*gupyy(i,j,k)*&
gupxz(i,j,k)*gupyz(i,j,k)*gyzx*gyyz-6.0*gupyy(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*gyyz*gyzz-2.0*gupyz(i,j,k)*gupxx(i,j,k)*gupzz(i,j,k)*&
gyzx*gzzx-2.0*gupyz(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*gyzx*gzzy-2.0*gupyz(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*gyzx*gzzz-2.0*&
gupyz(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*gyzy*gzzx-2.0*gupyz(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*gyzy*gzzy-2.0*gupyz(i,j,k)*gupxz(i,j,k)*&
gupzz(i,j,k)*gyzz*gzzx-gupyz(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxyz*gzzy+gxzy*gyzz)-gupyz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gyyz&
*gzzy+gyzy*gyzz)-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gyzz*gzzy-gupyz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(gyyz*gzzz+&
gyzz*gyzz)-2.0*gupyz(i,j,k)*gupzz(i,j,k)*gupzz(i,j,k)*gyzz*gzzz)
      MapleGenVar3 = gupyz(i,j,k)*gupxx(i,j,k)*gxyxy+gupxz(i,j,k)*gupyy(i,j,k)*gyyxy+gupzz(i,j,k)*gupxy(i,j,k)*gyyxz+2.0*&
gupxz(i,j,k)*gupzz(i,j,k)*gyzxz+gupzz(i,j,k)*gupxy(i,j,k)*gxyyz+gupyz(i,j,k)*gupxy(i,j,k)*gxyyy+gupyz(i,j,k)*gupxy(i,j,k)*gyyxy+gupyz(i,j,k)*&
gupyy(i,j,k)*gyyyy+gupxz(i,j,k)*gupxz(i,j,k)*gyzxx+gupxz(i,j,k)*gupyz(i,j,k)*gyyxz-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*gxzx*gxyz&
-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*gyyy*gyzz+gupxz(i,j,k)*gupxy(i,j,k)*gxyxy+gupzz(i,j,k)*gupyy(i,j,k)*gyyyz+gupxz(i,j,k)*&
gupxy(i,j,k)*gyyxx-2.0*gupzz(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*gyzz*gyzz+gupyz(i,j,k)*gupyz(i,j,k)*gyzyy+2.0*gupyz(i,j,k)*gupzz(i,j,k)*&
gyzyz+gupzz(i,j,k)*gupxx(i,j,k)*gxyxz+gupyz(i,j,k)*gupxz(i,j,k)*gxyyz+gupxz(i,j,k)*gupxx(i,j,k)*gxyxx+2.0*gupxz(i,j,k)*gupyz(i,j,k)*&
gyzxy-2.0*gupyz(i,j,k)*gupyy(i,j,k)*gupyy(i,j,k)*gyyy*gyyy-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*gyzy*gyyz+gupzz(i,j,k)*&
gupzz(i,j,k)*gyzzz
      MapleGenVar4 = MapleGenVar3+gupzz(i,j,k)*gupxz(i,j,k)*gxyzz-2.0*gupzz(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxzy*&
gyzx+gupyz(i,j,k)*gupyz(i,j,k)*gyyyz+gupxz(i,j,k)*gupxz(i,j,k)*gxyxz-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*gxyx*gxzz-2.0*&
gupzz(i,j,k)*gupzz(i,j,k)*gupzz(i,j,k)*gyzz*gzzz+gupzz(i,j,k)*gupyz(i,j,k)*gyyzz-2.0*gupyz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*gxyz*gyyx&
-2.0*gupxz(i,j,k)*gupxx(i,j,k)*gupxx(i,j,k)*gxxx*gxyx-2.0*gupxz(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gxyx*gxyx-2.0*gupxz(i,j,k)*&
gupxz(i,j,k)*gupxx(i,j,k)*gxyx*gxzx-gupxz(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*(gxxx*gyyx+gxyx*gxyx)
      MapleGenVar2 = MapleGenVar4-3.0*gupxz(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*(gxyx*gyzx+gxzx*gyyx)-&
gupxz(i,j,k)*gupxz(i,j,k)*gupxx(i,j,k)*(gxxx*gyzx+gxyx*gxzx)-2.0*gupxz(i,j,k)*gupxx(i,j,k)*gupzz(i,j,k)*(gxyx*gzzx+gxzx*&
gyzx)-4.0*gupxz(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*gxyx*gxyy-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*gxyx*gxzy-gupxz(i,j,k)*&
gupxy(i,j,k)*gupxy(i,j,k)*(gxxy*gyyx+gxyx*gxyy)-2.0*gupxz(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*(gxyx*gyyy+gxyy*gyyx)&
-2.0*gupxz(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxyx*gyzy+gxzy*gyyx)-gupxz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*(gxxy*gyzx+gxyx*&
gxzy)-2.0*gupxz(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxyx*gyzy+gxyy*gyzx)-2.0*gupxz(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gxyx*&
gzzy+gxzy*gyzx)-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupxx(i,j,k)*gxyx*gxxz-4.0*gupxz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*gxyx*gxyz&
-gupxz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*(gxxz*gyyx+gxyx*gxyz)
      MapleGenVar4 = MapleGenVar2-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*(gxyx*gyyz+gxyz*gyyx)&
-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxyx*gyzz+gxzz*gyyx)-gupxz(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxxz*gyzx+gxyx*&
gxzz)-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxyx*gyzz+gxyz*gyzx)-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gxyx*&
gzzz+gxzz*gyzx)-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*gxzx*gxyy-gupxz(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*(gxxx*gyyy+&
gxyx*gxyy)-2.0*gupxz(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxyy*gyzx+gxzx*gyyy)-gupxz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*(gxxx*&
gyzy+gxzx*gxyy)-2.0*gupxz(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gxyy*gzzx+gxzx*gyzy)-2.0*gupxz(i,j,k)*gupyy(i,j,k)*&
gupxy(i,j,k)*gxyy*gxyy-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gxyy*gxzy
      MapleGenVar3 = MapleGenVar4-gupxz(i,j,k)*gupyy(i,j,k)*gupxy(i,j,k)*(gxxy*gyyy+gxyy*gxyy)-2.0*&
gupxz(i,j,k)*gupyy(i,j,k)*gupyy(i,j,k)*gxyy*gyyy-3.0*gupxz(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*(gxyy*gyzy+gxzy*gyyy)-gupxz(i,j,k)*&
gupxz(i,j,k)*gupyy(i,j,k)*(gxxy*gyzy+gxyy*gxzy)-2.0*gupxz(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*(gxyy*gzzy+gxzy*gyzy)&
-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gxyy*gxzz-gupxz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*(gxxz*gyyy+gxyy*gxyz)-2.0*&
gupxz(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*(gxyy*gyyz+gxyz*gyyy)-2.0*gupxz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxyy*gyzz+gxzz*&
gyyy)-gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxxz*gyzy+gxyy*gxzz)-2.0*gupxz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxyy*gyzz+&
gxyz*gyzy)-2.0*gupxz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(gxyy*gzzz+gxzz*gyzy)-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupxx(i,j,k)*&
gxxx*gxyz
      MapleGenVar4 = MapleGenVar3-gupxz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*(gxxx*gyyz+gxyx*gxyz)-2.0*&
gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxyz*gyzx+gxzx*gyyz)-gupxz(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxxx*gyzz+gxzx*gxyz)&
-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gxyz*gzzx+gxzx*gyzz)-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gxzy*gxyz-&
gupxz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*(gxxy*gyyz+gxyy*gxyz)-2.0*gupxz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxyz*gyzy+gxzy*&
gyyz)-gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxxy*gyzz+gxzy*gxyz)-2.0*gupxz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(gxyz*gzzy+&
gxzy*gyzz)-2.0*gupxz(i,j,k)*gupzz(i,j,k)*gupxy(i,j,k)*gxyz*gxyz-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*gxyz*gxzz-&
gupxz(i,j,k)*gupzz(i,j,k)*gupxy(i,j,k)*(gxxz*gyyz+gxyz*gxyz)
      MapleGenVar1 = MapleGenVar4-3.0*gupxz(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*(gxyz*gyzz+gxzz*gyyz)-&
gupxz(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gxxz*gyzz+gxyz*gxzz)-2.0*gupxz(i,j,k)*gupzz(i,j,k)*gupzz(i,j,k)*(gxyz*gzzz+gxzz*&
gyzz)-gupyz(i,j,k)*gupxx(i,j,k)*gupxx(i,j,k)*(gxxx*gyyx+gxyx*gxyx)-2.0*gupyz(i,j,k)*gupxx(i,j,k)*gupyy(i,j,k)*gyyx*gyyx&
-4.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxx(i,j,k)*gyyx*gyzx-2.0*gupxz(i,j,k)*gupxx(i,j,k)*gupyy(i,j,k)*gxyx*gyyx-2.0*gupxz(i,j,k)*&
gupxx(i,j,k)*gupyz(i,j,k)*gxyx*gyzx-2.0*gupxz(i,j,k)*gupxy(i,j,k)*gupxx(i,j,k)*gxyx*gxxy-2.0*gupxz(i,j,k)*gupxy(i,j,k)*gupxx(i,j,k)*&
gxxx*gxyy-2.0*gupxz(i,j,k)*gupyy(i,j,k)*gupxx(i,j,k)*gxxy*gxyy-2.0*gupxz(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gxyy*gyzy-2.0*&
gupxz(i,j,k)*gupyz(i,j,k)*gupxx(i,j,k)*gxyy*gxxz-4.0*gupxz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*gxyy*gxyz
      MapleGenVar4 = MapleGenVar1-2.0*gupxz(i,j,k)*gupyz(i,j,k)*gupxx(i,j,k)*gxxy*gxyz-2.0*gupxz(i,j,k)*&
gupzz(i,j,k)*gupxx(i,j,k)*gxxz*gxyz-2.0*gupxz(i,j,k)*gupzz(i,j,k)*gupyy(i,j,k)*gxyz*gyyz-2.0*gupxz(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*&
gxyz*gyzz-4.0*gupyz(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gxyx*gyyx-2.0*gupyz(i,j,k)*gupxx(i,j,k)*gupzz(i,j,k)*(gyyx*gzzx+gyzx&
*gyzx)-gupyz(i,j,k)*gupxy(i,j,k)*gupxx(i,j,k)*(gxxy*gyyx+gxyx*gxyy)-2.0*gupyz(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*(gxyx*gyyy&
+gxyy*gyyx)-2.0*gupyz(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*gxyy*gyyx-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*gyyx*gyzy-&
gupyz(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxyy*gyzx+gxzy*gyyx)
      MapleGenVar3 = MapleGenVar4-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*(gyyx*gyzy+gyzx*gyyy)&
-2.0*gupyz(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gyyx*gzzy+gyzx*gyzy)-gupyz(i,j,k)*gupxz(i,j,k)*gupxx(i,j,k)*(gxxz*gyyx+gxyx*&
gxyz)-2.0*gupyz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*(gxyx*gyyz+gxyz*gyyx)-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*gyyx*&
gyzz-gupyz(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxyz*gyzx+gxzz*gyyx)-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gyyx*gyzz+&
gyzx*gyyz)-2.0*gupyz(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gyyx*gzzz+gyzx*gyzz)-gupyz(i,j,k)*gupxy(i,j,k)*gupxx(i,j,k)*(gxxx*&
gyyy+gxyx*gxyy)-2.0*gupyz(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*gxyx*gyyy-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*gyzx*gyyy&
-gupyz(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxyx*gyzy+gxzx*gyyy)-2.0*gupyz(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gyyy*gzzx+gyzx*&
gyzy)
      MapleGenVar4 = MapleGenVar3-gupyz(i,j,k)*gupyy(i,j,k)*gupxx(i,j,k)*(gxxy*gyyy+gxyy*gxyy)-4.0*&
gupyz(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*gyyy*gyzy-2.0*gupyz(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*(gyyy*gzzy+gyzy*gyzy)-gupyz(i,j,k)*&
gupyz(i,j,k)*gupxx(i,j,k)*(gxxz*gyyy+gxyy*gxyz)-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*(gxyy*gyyz+gxyz*gyyy)&
-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*gxyz*gyyy-4.0*gupyz(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*gyyy*gyyz-gupyz(i,j,k)*gupyz(i,j,k)*&
gupxz(i,j,k)*(gxyz*gyzy+gxzz*gyyy)-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gyyy*gyzz+gyzy*gyyz)-2.0*&
gupyz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(gyyy*gzzz+gyzy*gyzz)-gupyz(i,j,k)*gupxz(i,j,k)*gupxx(i,j,k)*(gxxx*gyyz+gxyx*gxyz)&
-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*gyzx*gyyz
      MapleGenVar2 = MapleGenVar4-gupyz(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxyx*gyzz+gxzx*gyyz)-2.0*&
gupyz(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gyyz*gzzx+gyzx*gyzz)-gupyz(i,j,k)*gupyz(i,j,k)*gupxx(i,j,k)*(gxxy*gyyz+gxyy*gxyz)&
-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*gxyy*gyyz-gupyz(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxyy*gyzz+gxzy*gyyz)-2.0*&
gupyz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(gyyz*gzzy+gyzy*gyzz)-gupyz(i,j,k)*gupzz(i,j,k)*gupxx(i,j,k)*(gxxz*gyyz+gxyz*gxyz)&
-2.0*gupyz(i,j,k)*gupzz(i,j,k)*gupyy(i,j,k)*gyyz*gyyz-4.0*gupyz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gyyz*gyzz-2.0*gupyz(i,j,k)*&
gupzz(i,j,k)*gupzz(i,j,k)*(gyyz*gzzz+gyzz*gyzz)-gupzz(i,j,k)*gupxx(i,j,k)*gupxx(i,j,k)*(gxxx*gyzx+gxyx*gxzx)-gupzz(i,j,k)&
*gupxx(i,j,k)*gupxy(i,j,k)*(gxyx*gyzx+gxzx*gyyx)-2.0*gupzz(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*gyzx*gyzx-2.0*gupzz(i,j,k)*&
gupzz(i,j,k)*gupxx(i,j,k)*gyzx*gzzx
      MapleGenVar4 = MapleGenVar2-gupzz(i,j,k)*gupxy(i,j,k)*gupxx(i,j,k)*(gxxy*gyzx+gxyx*gxzy)-2.0*&
gupzz(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*(gxyx*gyzy+gxyy*gyzx)-gupzz(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*(gxyy*gyzx+gxzy*gyyx)&
-2.0*gupzz(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*(gyyx*gyzy+gyzx*gyyy)-2.0*gupzz(i,j,k)*gupzz(i,j,k)*gupxy(i,j,k)*gyzx*gzzy-&
gupzz(i,j,k)*gupxz(i,j,k)*gupxx(i,j,k)*(gxxz*gyzx+gxyx*gxzz)-2.0*gupzz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*(gxyx*gyzz+gxyz*&
gyzx)-gupzz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*(gxyz*gyzx+gxzz*gyyx)-2.0*gupzz(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*(gyyx*gyzz+&
gyzx*gyyz)-2.0*gupzz(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*gxzz*gyzx-2.0*gupzz(i,j,k)*gupzz(i,j,k)*gupxz(i,j,k)*gyzx*gzzz-&
gupzz(i,j,k)*gupxy(i,j,k)*gupxx(i,j,k)*(gxxx*gyzy+gxzx*gxyy)
      MapleGenVar3 = MapleGenVar4-gupzz(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*(gxyx*gyzy+gxzx*gyyy)-4.0*&
gupyz(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*gyyx*gyyy-4.0*gupyz(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gyyx*gyyz-4.0*gupyz(i,j,k)*gupyy(i,j,k)*&
gupxy(i,j,k)*gxyy*gyyy-2.0*gupyz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*gxyx*gyyz-4.0*gupyz(i,j,k)*gupzz(i,j,k)*gupxy(i,j,k)*gxyz*gyyz&
-2.0*gupzz(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gxyx*gyzx-2.0*gupzz(i,j,k)*gupxx(i,j,k)*gupyy(i,j,k)*gyyx*gyzx-2.0*gupzz(i,j,k)*&
gupxx(i,j,k)*gupxz(i,j,k)*gxzx*gyzx-4.0*gupzz(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gyzx*gyzy-4.0*gupzz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*&
gyzx*gyzz-2.0*gupzz(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxzx*gyzy-2.0*gupzz(i,j,k)*gupzz(i,j,k)*gupxy(i,j,k)*gyzy*gzzx
      MapleGenVar4 = MapleGenVar3-gupzz(i,j,k)*gupyy(i,j,k)*gupxx(i,j,k)*(gxxy*gyzy+gxyy*gxzy)-gupzz(i,j,k)&
*gupyy(i,j,k)*gupxy(i,j,k)*(gxyy*gyzy+gxzy*gyyy)-2.0*gupzz(i,j,k)*gupyy(i,j,k)*gupyy(i,j,k)*gyyy*gyzy-2.0*gupzz(i,j,k)*&
gupyy(i,j,k)*gupyz(i,j,k)*gyzy*gyzy-2.0*gupzz(i,j,k)*gupzz(i,j,k)*gupyy(i,j,k)*gyzy*gzzy-gupzz(i,j,k)*gupyz(i,j,k)*gupxx(i,j,k)*(gxxz*&
gyzy+gxyy*gxzz)-2.0*gupzz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*(gxyy*gyzz+gxyz*gyzy)-gupzz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*(&
gxyz*gyzy+gxzz*gyyy)-2.0*gupzz(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*(gyyy*gyzz+gyzy*gyyz)-4.0*gupzz(i,j,k)*&
gupyz(i,j,k)*gupyz(i,j,k)*gyzy*gyzz-2.0*gupzz(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*gyzy*gzzz-gupzz(i,j,k)*gupxz(i,j,k)*gupxx(i,j,k)*(gxxx*&
gyzz+gxzx*gxyz)-gupzz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*(gxyx*gyzz+gxzx*gyyz)
      CAZzy  = Gamzy -(MapleGenVar4-2.0*gupzz(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*gxzx*gyzz-2.0*gupzz(i,j,k)*&
gupzz(i,j,k)*gupxz(i,j,k)*gyzz*gzzx-gupzz(i,j,k)*gupyz(i,j,k)*gupxx(i,j,k)*(gxxy*gyzz+gxzy*gxyz)-gupzz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)&
*(gxyy*gyzz+gxzy*gyyz)-2.0*gupzz(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*gyzz*gzzy-gupzz(i,j,k)*gupzz(i,j,k)*gupxx(i,j,k)*(gxxz*&
gyzz+gxyz*gxzz)-2.0*gupzz(i,j,k)*gupzz(i,j,k)*gupxy(i,j,k)*gxyz*gyzz-gupzz(i,j,k)*gupzz(i,j,k)*gupxy(i,j,k)*(gxyz*gyzz+&
gxzz*gyyz)-2.0*gupzz(i,j,k)*gupzz(i,j,k)*gupyy(i,j,k)*gyyz*gyzz-2.0*gupzz(i,j,k)*gupzz(i,j,k)*gupxz(i,j,k)*gxzz*gyzz-2.0*&
gupzz(i,j,k)*gupyy(i,j,k)*gupxy(i,j,k)*gxyy*gyzy-2.0*gupzz(i,j,k)*gupyy(i,j,k)*gupxz(i,j,k)*gxzy*gyzy-2.0*gupzz(i,j,k)*gupyz(i,j,k)*&
gupxz(i,j,k)*gxzz*gyzy-2.0*gupzz(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*gxzy*gyzz)
      MapleGenVar3 = gupxx(i,j,k)*gupyy(i,j,k)*gyzxy-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*gxzx*gzzz+gupxz(i,j,k)*&
gupzz(i,j,k)*gzzzz+gupxy(i,j,k)*gupxy(i,j,k)*gxzyy+2.0*gupxx(i,j,k)*gupxy(i,j,k)*gxzxy+2.0*gupxy(i,j,k)*gupxz(i,j,k)*gxzyz-2.0*&
gupxx(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gxzx*gxzx+gupxz(i,j,k)*gupxy(i,j,k)*gyzxz+gupxx(i,j,k)*gupyz(i,j,k)*gzzxy+gupxy(i,j,k)*gupxy(i,j,k)*&
gyzxy+gupxz(i,j,k)*gupyz(i,j,k)*gzzyz+gupxy(i,j,k)*gupyy(i,j,k)*gyzyy+gupxz(i,j,k)*gupyy(i,j,k)*gyzyz+gupxz(i,j,k)*gupxz(i,j,k)*gzzxz&
-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*gxyy*gyzx+2.0*gupxx(i,j,k)*gupxz(i,j,k)*gxzxz+gupxx(i,j,k)*gupxx(i,j,k)*gxzxx-2.0*&
gupxz(i,j,k)*gupzz(i,j,k)*gupzz(i,j,k)*gzzz*gzzz+gupxx(i,j,k)*gupyz(i,j,k)*gyzxz+gupxz(i,j,k)*gupxz(i,j,k)*gxzzz-2.0*gupxx(i,j,k)*gupxx(i,j,k)&
*gupxx(i,j,k)*gxxx*gxzx+gupxy(i,j,k)*gupxz(i,j,k)*gzzxy+gupxx(i,j,k)*gupxy(i,j,k)*gyzxx+gupxy(i,j,k)*gupyz(i,j,k)*gyzyz-2.0*&
gupxy(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*gxyx*gyzy
      MapleGenVar4 = MapleGenVar3+gupxx(i,j,k)*gupxz(i,j,k)*gzzxx+gupxx(i,j,k)*gupzz(i,j,k)*gzzxz+gupxy(i,j,k)*&
gupyz(i,j,k)*gzzyy+gupxy(i,j,k)*gupzz(i,j,k)*gzzyz+gupxz(i,j,k)*gupyz(i,j,k)*gyzzz-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*gxzz*gzzx&
-2.0*gupxx(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gxyx*gxzx-2.0*gupxx(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*(gxxx*gyzx+gxyx*gxzx)-&
gupxx(i,j,k)*gupxx(i,j,k)*gupyy(i,j,k)*(gxyx*gyzx+gxzx*gyyx)-2.0*gupxx(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*gxzx*gyzx-2.0*&
gupxx(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*(gxxx*gzzx+gxzx*gxzx)-gupxx(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*(gxyx*gzzx+gxzx*gyzx)

      MapleGenVar2 = MapleGenVar4-2.0*gupxx(i,j,k)*gupxx(i,j,k)*gupzz(i,j,k)*gxzx*gzzx-2.0*gupxx(i,j,k)*&
gupxx(i,j,k)*gupxy(i,j,k)*gxzx*gxxy-2.0*gupxx(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*gxzx*gxyy-2.0*gupxx(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*(&
gxxy*gyzx+gxzx*gxyy)-gupxx(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*(gxyy*gyzx+gxzx*gyyy)-2.0*gupxx(i,j,k)*gupxy(i,j,k)*&
gupyz(i,j,k)*(gxzx*gyzy+gxzy*gyzx)-2.0*gupxx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxxy*gzzx+gxzx*gxzy)-gupxx(i,j,k)*&
gupxy(i,j,k)*gupyz(i,j,k)*(gxyy*gzzx+gxzx*gyzy)-2.0*gupxx(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gxzx*gzzy+gxzy*gzzx)&
-2.0*gupxx(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gxzx*gxxz-4.0*gupxx(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*gxzx*gxzz-2.0*gupxx(i,j,k)*&
gupxz(i,j,k)*gupxy(i,j,k)*(gxxz*gyzx+gxzx*gxyz)-gupxx(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*(gxyz*gyzx+gxzx*gyyz)-2.0*&
gupxx(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxzx*gyzz+gxzz*gyzx)
      MapleGenVar4 = MapleGenVar2-2.0*gupxx(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxxz*gzzx+gxzx*gxzz)-&
gupxx(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxyz*gzzx+gxzx*gyzz)-2.0*gupxx(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gxzx*gzzz+gxzz*&
gzzx)-2.0*gupxx(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gxxx*gxzy-2.0*gupxx(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*gxyx*gxzy-2.0*gupxx(i,j,k)&
*gupxy(i,j,k)*gupxy(i,j,k)*(gxxx*gyzy+gxyx*gxzy)-gupxx(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*(gxyx*gyzy+gxzy*gyyx)-2.0*&
gupxx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxxx*gzzy+gxzx*gxzy)-gupxx(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxyx*gzzy+gxzy*gyzx)&
-2.0*gupxx(i,j,k)*gupxx(i,j,k)*gupyy(i,j,k)*gxxy*gxzy-2.0*gupxx(i,j,k)*gupyy(i,j,k)*gupxz(i,j,k)*gxzy*gxzy-2.0*gupxx(i,j,k)*&
gupyy(i,j,k)*gupxy(i,j,k)*(gxxy*gyzy+gxyy*gxzy)
      MapleGenVar3 = MapleGenVar4-gupxx(i,j,k)*gupyy(i,j,k)*gupyy(i,j,k)*(gxyy*gyzy+gxzy*gyyy)-2.0*&
gupxx(i,j,k)*gupyy(i,j,k)*gupxz(i,j,k)*(gxxy*gzzy+gxzy*gxzy)-gupxx(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*(gxyy*gzzy+gxzy*gyzy)&
-2.0*gupxx(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*gxzy*gxxz-2.0*gupxx(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*(gxxz*gyzy+gxzy*gxyz)-&
gupxx(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*(gxyz*gyzy+gxzy*gyyz)-2.0*gupxx(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxzy*gyzz+gxzz*&
gyzy)-2.0*gupxx(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxxz*gzzy+gxzy*gxzz)-gupxx(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxyz*gzzy+&
gxzy*gyzz)-2.0*gupxx(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(gxzy*gzzz+gxzz*gzzy)-2.0*gupxx(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*&
gxxx*gxzz-4.0*gupxx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxzx*gxzy-2.0*gupxx(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*gxzx*gxyz
      MapleGenVar4 = MapleGenVar3-2.0*gupxx(i,j,k)*gupyy(i,j,k)*gupxy(i,j,k)*gxyy*gxzy-2.0*gupxx(i,j,k)*&
gupyy(i,j,k)*gupyz(i,j,k)*gxzy*gyzy-2.0*gupxx(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*gxzy*gzzy-2.0*gupxx(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*&
gxzy*gxyz-4.0*gupxx(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*gxzy*gxzz-2.0*gupxx(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*gxyx*gxzz-2.0*&
gupxx(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*(gxxx*gyzz+gxyx*gxzz)-gupxx(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*(gxyx*gyzz+gxzz*gyyx)&
-2.0*gupxx(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxxx*gzzz+gxzx*gxzz)-gupxx(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxyx*gzzz+gxzz*&
gyzx)-2.0*gupxx(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*gxxy*gxzz-2.0*gupxx(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*(gxxy*gyzz+gxyy*&
gxzz)
      MapleGenVar1 = MapleGenVar4-gupxx(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*(gxyy*gyzz+gxzz*gyyy)-2.0*&
gupxx(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxxy*gzzz+gxzy*gxzz)-gupxx(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxyy*gzzz+gxzz*gyzy)&
-2.0*gupxx(i,j,k)*gupxx(i,j,k)*gupzz(i,j,k)*gxxz*gxzz-2.0*gupxx(i,j,k)*gupzz(i,j,k)*gupxz(i,j,k)*gxzz*gxzz-2.0*gupxx(i,j,k)*&
gupzz(i,j,k)*gupxy(i,j,k)*(gxxz*gyzz+gxyz*gxzz)-gupxx(i,j,k)*gupzz(i,j,k)*gupyy(i,j,k)*(gxyz*gyzz+gxzz*gyyz)-2.0*&
gupxx(i,j,k)*gupzz(i,j,k)*gupxz(i,j,k)*(gxxz*gzzz+gxzz*gxzz)-gupxx(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*(gxyz*gzzz+gxzz*gyzz)&
-2.0*gupxx(i,j,k)*gupzz(i,j,k)*gupzz(i,j,k)*gxzz*gzzz-gupxy(i,j,k)*gupxy(i,j,k)*gupxx(i,j,k)*(gxyx*gyzx+gxzx*gyyx)-2.0*&
gupxy(i,j,k)*gupxy(i,j,k)*gupxx(i,j,k)*gxyx*gyzx-2.0*gupxy(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*gyzx*gyzx-3.0*gupxy(i,j,k)*gupxx(i,j,k)*&
gupxz(i,j,k)*(gxyx*gzzx+gxzx*gyzx)
      MapleGenVar4 = MapleGenVar1-gupxy(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*(gyyx*gzzx+gyzx*gyzx)-gupxy(i,j,k)&
*gupxy(i,j,k)*gupxy(i,j,k)*(gxyy*gyzx+gxzx*gyyy)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxzx*gyzy+gxzy*gyzx)&
-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*gyzx*gyyy-4.0*gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gyzx*gyzy-2.0*gupxy(i,j,k)*&
gupxy(i,j,k)*gupxz(i,j,k)*(gxyy*gzzx+gxzy*gyzx)-gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gyyy*gzzx+gyzx*gyzy)-2.0*&
gupxy(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gyzx*gzzy+gyzy*gzzx)-gupxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxyz*gyzx+gxzx*gyyz)&
-2.0*gupxy(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxzx*gyzz+gxzz*gyzx)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxyz*gyzx
      MapleGenVar3 = MapleGenVar4-2.0*gupxy(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxyz*gzzx+gxzz*gyzx)-&
gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gyyz*gzzx+gyzx*gyzz)-2.0*gupxy(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gyzx*gzzz+gyzz*&
gzzx)-gupxy(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*(gxyx*gyzy+gxzy*gyyx)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*gyyx*gyzy&
-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxyx*gzzy+gxzx*gyzy)-gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gyyx*gzzy+gyzx*&
gyzy)-gupxy(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*(gxyy*gyzy+gxzy*gyyy)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*gxyy*gyzy&
-2.0*gupxy(i,j,k)*gupyy(i,j,k)*gupyy(i,j,k)*gyyy*gyzy-2.0*gupxy(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gyzy*gyzy-3.0*gupxy(i,j,k)*&
gupyy(i,j,k)*gupxz(i,j,k)*(gxyy*gzzy+gxzy*gyzy)-gupxy(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*(gyyy*gzzy+gyzy*gyzy)
      MapleGenVar4 = MapleGenVar3-gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxyz*gyzy+gxzy*gyyz)-2.0*&
gupxy(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxzy*gyzz+gxzz*gyzy)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gxyz*gyzy-4.0*&
gupxy(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*gyzy*gyzz-2.0*gupxy(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxyz*gzzy+gxzz*gyzy)-gupxy(i,j,k)*&
gupyz(i,j,k)*gupyz(i,j,k)*(gyyz*gzzy+gyzy*gyzz)-2.0*gupxy(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(gyzy*gzzz+gyzz*gzzy)-&
gupxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxyx*gyzz+gxzz*gyyx)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxyx*gyzz-2.0*&
gupxy(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxyx*gzzz+gxzx*gyzz)-gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gyyx*gzzz+gyzx*gyzz)&
-gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxyy*gyzz+gxzz*gyyy)
      MapleGenVar2 = MapleGenVar4-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gxyy*gyzz-2.0*gupxy(i,j,k)*&
gupyz(i,j,k)*gupxz(i,j,k)*(gxyy*gzzz+gxzy*gyzz)-gupxy(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gyyy*gzzz+gyzy*gyzz)-gupxy(i,j,k)&
*gupxy(i,j,k)*gupzz(i,j,k)*(gxyz*gyzz+gxzz*gyyz)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*gxyz*gyzz-2.0*gupxy(i,j,k)*&
gupzz(i,j,k)*gupyz(i,j,k)*gyzz*gyzz-3.0*gupxy(i,j,k)*gupzz(i,j,k)*gupxz(i,j,k)*(gxyz*gzzz+gxzz*gyzz)-gupxy(i,j,k)*gupzz(i,j,k)*&
gupyz(i,j,k)*(gyyz*gzzz+gyzz*gyzz)-2.0*gupxy(i,j,k)*gupzz(i,j,k)*gupzz(i,j,k)*gyzz*gzzz-2.0*gupxx(i,j,k)*gupyz(i,j,k)*&
gupxy(i,j,k)*gxyy*gxzz-2.0*gupxx(i,j,k)*gupzz(i,j,k)*gupxy(i,j,k)*gxyz*gxzz-2.0*gupxx(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*gxzz*gyzz&
-2.0*gupxy(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gxzx*gyzx-2.0*gupxy(i,j,k)*gupxx(i,j,k)*gupyy(i,j,k)*gyyx*gyzx
      MapleGenVar4 = MapleGenVar2-2.0*gupxy(i,j,k)*gupxx(i,j,k)*gupzz(i,j,k)*gyzx*gzzx-2.0*gupxy(i,j,k)*&
gupxz(i,j,k)*gupyy(i,j,k)*gyzx*gyyz-4.0*gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gyzx*gyzz-2.0*gupxy(i,j,k)*gupyy(i,j,k)*gupxz(i,j,k)*&
gxzy*gyzy-2.0*gupxy(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*gyzy*gzzy-2.0*gupxy(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*gyzy*gyyz-2.0*&
gupxy(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gyyx*gyzz-2.0*gupxy(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*gyyy*gyzz-2.0*gupxy(i,j,k)*gupzz(i,j,k)*&
gupxz(i,j,k)*gxzz*gyzz-2.0*gupxy(i,j,k)*gupzz(i,j,k)*gupyy(i,j,k)*gyyz*gyzz-4.0*gupxz(i,j,k)*gupxz(i,j,k)*gupxx(i,j,k)*gxzx*gzzx&
-gupxz(i,j,k)*gupxx(i,j,k)*gupyy(i,j,k)*(gyyx*gzzx+gyzx*gyzx)
      MapleGenVar3 = MapleGenVar4-2.0*gupxz(i,j,k)*gupxx(i,j,k)*gupzz(i,j,k)*gzzx*gzzx-gupxz(i,j,k)*gupxy(i,j,k)*&
gupxy(i,j,k)*(gxyy*gzzx+gxzx*gyzy)-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*(gxzx*gzzy+gxzy*gzzx)-gupxz(i,j,k)*&
gupxy(i,j,k)*gupyy(i,j,k)*(gyyy*gzzx+gyzx*gyzy)-2.0*gupxz(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gyzx*gzzy+gyzy*gzzx)&
-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*gxzy*gzzx-gupxz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*(gxyz*gzzx+gxzx*gyzz)-2.0*&
gupxz(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxzx*gzzz+gxzz*gzzx)-gupxz(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*(gyyz*gzzx+gyzx*gyzz)&
-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gyzx*gzzz+gyzz*gzzx)-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gyzz*gzzx&
-4.0*gupxz(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*gzzx*gzzz-gupxz(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*(gxyx*gzzy+gxzy*gyzx)
      MapleGenVar4 = MapleGenVar3-gupxz(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*(gyyx*gzzy+gyzx*gyzy)-2.0*&
gupxz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*gxzx*gzzy-4.0*gupxz(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gxzy*gzzy-gupxz(i,j,k)*gupyy(i,j,k)*gupyy(i,j,k)*(&
gyyy*gzzy+gyzy*gyzy)-2.0*gupxz(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*gzzy*gzzy-gupxz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*(gxyz*&
gzzy+gxzy*gyzz)-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxzy*gzzz+gxzz*gzzy)-gupxz(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*(&
gyyz*gzzy+gyzy*gyzz)-2.0*gupxz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gyzy*gzzz+gyzz*gzzy)-2.0*gupxz(i,j,k)*&
gupxz(i,j,k)*gupyz(i,j,k)*gxzz*gzzy-2.0*gupxz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*gyzz*gzzy-gupxz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*(gxyx*&
gzzz+gxzz*gyzx)-gupxz(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*(gyyx*gzzz+gyzx*gyzz)
      CAZxz = Gamxz - (MapleGenVar4-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gyzx*gzzz-gupxz(i,j,k)*gupyz(i,j,k)*&
gupxy(i,j,k)*(gxyy*gzzz+gxzz*gyzy)-gupxz(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*(gyyy*gzzz+gyzy*gyzz)-2.0*gupxz(i,j,k)*&
gupxz(i,j,k)*gupyz(i,j,k)*gxzy*gzzz-2.0*gupxz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*gyzy*gzzz-4.0*gupxz(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*&
gxzz*gzzz-gupxz(i,j,k)*gupzz(i,j,k)*gupyy(i,j,k)*(gyyz*gzzz+gyzz*gyzz)-4.0*gupxz(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*gyzx*&
gzzx-2.0*gupxz(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gyzy*gzzx-4.0*gupxz(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*gzzx*gzzy-2.0*gupxz(i,j,k)*&
gupxy(i,j,k)*gupyz(i,j,k)*gyzx*gzzy-4.0*gupxz(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gyzy*gzzy-4.0*gupxz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*&
gzzy*gzzz-4.0*gupxz(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*gyzz*gzzz)
      MapleGenVar3 = gupyz(i,j,k)*gupyz(i,j,k)*gzzyz+gupyy(i,j,k)*gupyy(i,j,k)*gyzyy+gupyz(i,j,k)*gupzz(i,j,k)*gzzzz-2.0*&
gupyy(i,j,k)*gupyy(i,j,k)*gupyy(i,j,k)*gyyy*gyzy+2.0*gupyy(i,j,k)*gupyz(i,j,k)*gyzyz+gupxy(i,j,k)*gupyz(i,j,k)*gzzxy-2.0*gupyy(i,j,k)*&
gupyy(i,j,k)*gupyz(i,j,k)*gyzy*gyzy+gupyy(i,j,k)*gupxx(i,j,k)*gxzxy+gupyz(i,j,k)*gupyz(i,j,k)*gyzzz+gupxy(i,j,k)*gupxz(i,j,k)*gxzxz+2.0&
*gupxy(i,j,k)*gupyz(i,j,k)*gyzxz+gupxy(i,j,k)*gupxx(i,j,k)*gxzxx+gupyy(i,j,k)*gupxz(i,j,k)*gzzxy+gupyy(i,j,k)*gupxz(i,j,k)*gxzyz+gupxy(i,j,k)*&
gupzz(i,j,k)*gzzxz-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*gxzx*gxyy+gupyz(i,j,k)*gupxy(i,j,k)*gxzyz-2.0*gupyz(i,j,k)*gupzz(i,j,k)*&
gupzz(i,j,k)*gzzz*gzzz+gupxy(i,j,k)*gupxz(i,j,k)*gzzxx-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*gyzy*gzzz-2.0*gupyz(i,j,k)*&
gupyz(i,j,k)*gupyz(i,j,k)*gyzz*gzzy-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxx(i,j,k)*gxxx*gxzy+2.0*gupxy(i,j,k)*gupyy(i,j,k)*gyzxy+&
gupyz(i,j,k)*gupxz(i,j,k)*gxzzz-2.0*gupxy(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gxzx*gzzz+gxzz*gzzx)
      MapleGenVar4 = MapleGenVar3+gupxy(i,j,k)*gupxy(i,j,k)*gyzxx+gupxy(i,j,k)*gupxy(i,j,k)*gxzxy+gupyy(i,j,k)*&
gupyz(i,j,k)*gzzyy+gupyz(i,j,k)*gupxz(i,j,k)*gzzxz+gupyz(i,j,k)*gupxx(i,j,k)*gxzxz+gupyy(i,j,k)*gupzz(i,j,k)*gzzyz-2.0*gupxy(i,j,k)*&
gupxy(i,j,k)*gupxy(i,j,k)*gxyx*gxzy-2.0*gupyy(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*gxyy*gyzz+gupyy(i,j,k)*gupxy(i,j,k)*gxzyy-2.0*&
gupxy(i,j,k)*gupxx(i,j,k)*gupxx(i,j,k)*gxxx*gxzx-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxx(i,j,k)*gxyx*gxzx-2.0*gupxy(i,j,k)*gupxx(i,j,k)*&
gupxz(i,j,k)*gxzx*gxzx
      MapleGenVar2 = MapleGenVar4-gupxy(i,j,k)*gupxy(i,j,k)*gupxx(i,j,k)*(gxxx*gyzx+gxyx*gxzx)-2.0*&
gupxy(i,j,k)*gupxx(i,j,k)*gupyy(i,j,k)*(gxyx*gyzx+gxzx*gyyx)-gupxy(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*(gxxx*gzzx+gxzx*gxzx)&
-3.0*gupxy(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*(gxyx*gzzx+gxzx*gyzx)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxx(i,j,k)*gxzx*gxxy&
-4.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxzx*gxzy-gupxy(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*(gxxy*gyzx+gxzx*gxyy)-2.0*&
gupxy(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*(gxyy*gyzx+gxzx*gyyy)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxzx*gyzy+gxzy*&
gyzx)-gupxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxxy*gzzx+gxzx*gxzy)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxyy*gzzx+&
gxzx*gyzy)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gxzx*gzzy+gxzy*gzzx)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*&
gxzx*gxyz-4.0*gupxy(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*gxzx*gxzz
      MapleGenVar4 = MapleGenVar2-gupxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxxz*gyzx+gxzx*gxyz)-2.0*&
gupxy(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*(gxyz*gyzx+gxzx*gyyz)-2.0*gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxzx*gyzz+gxzz*&
gyzx)-gupxy(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxxz*gzzx+gxzx*gxzz)-2.0*gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxyz*gzzx+&
gxzx*gyzz)-2.0*gupxy(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*gxzx*gyzx-2.0*gupxy(i,j,k)*gupxx(i,j,k)*gupzz(i,j,k)*gxzx*gzzx-2.0*&
gupxy(i,j,k)*gupxz(i,j,k)*gupxx(i,j,k)*gxzx*gxxz-gupxy(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*(gxxx*gyzy+gxyx*gxzy)-2.0*gupxy(i,j,k)*&
gupxy(i,j,k)*gupyy(i,j,k)*(gxyx*gyzy+gxzy*gyyx)-gupxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxxx*gzzy+gxzx*gxzy)-2.0*&
gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxyx*gzzy+gxzy*gyzx)
      MapleGenVar3 = MapleGenVar4-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*gxyy*gxzy-2.0*gupxy(i,j,k)*&
gupyy(i,j,k)*gupxz(i,j,k)*gxzy*gxzy-gupxy(i,j,k)*gupxy(i,j,k)*gupyy(i,j,k)*(gxxy*gyzy+gxyy*gxzy)-2.0*gupxy(i,j,k)*gupyy(i,j,k)*&
gupyy(i,j,k)*(gxyy*gyzy+gxzy*gyyy)-gupxy(i,j,k)*gupyy(i,j,k)*gupxz(i,j,k)*(gxxy*gzzy+gxzy*gxzy)-3.0*gupxy(i,j,k)*&
gupyy(i,j,k)*gupyz(i,j,k)*(gxyy*gzzy+gxzy*gyzy)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gxzy*gxyz-gupxy(i,j,k)*gupxy(i,j,k)*&
gupyz(i,j,k)*(gxxz*gyzy+gxzy*gxyz)-2.0*gupxy(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*(gxyz*gyzy+gxzy*gyyz)-2.0*&
gupxy(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxzy*gyzz+gxzz*gyzy)-gupxy(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxxz*gzzy+gxzy*gxzz)&
-2.0*gupxy(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxyz*gzzy+gxzy*gyzz)-2.0*gupxy(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(gxzy*gzzz+&
gxzz*gzzy)
      MapleGenVar4 = MapleGenVar3-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxyx*gxzz-gupxy(i,j,k)*gupxy(i,j,k)*&
gupxz(i,j,k)*(gxxx*gyzz+gxyx*gxzz)-2.0*gupxy(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*(gxyx*gyzz+gxzz*gyyx)-gupxy(i,j,k)*&
gupxz(i,j,k)*gupxz(i,j,k)*(gxxx*gzzz+gxzx*gxzz)-2.0*gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxyx*gzzz+gxzz*gyzx)&
-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gxyy*gxzz-gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxxy*gyzz+gxyy*gxzz)-2.0*&
gupxy(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*(gxyy*gyzz+gxzz*gyyy)-gupxy(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxxy*gzzz+gxzy*gxzz)&
-2.0*gupxy(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxyy*gzzz+gxzz*gyzy)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*gxyz*gxzz&
-2.0*gupxy(i,j,k)*gupzz(i,j,k)*gupxz(i,j,k)*gxzz*gxzz
      MapleGenVar1 = MapleGenVar4-gupxy(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gxxz*gyzz+gxyz*gxzz)-2.0*&
gupxy(i,j,k)*gupzz(i,j,k)*gupyy(i,j,k)*(gxyz*gyzz+gxzz*gyyz)-gupxy(i,j,k)*gupzz(i,j,k)*gupxz(i,j,k)*(gxxz*gzzz+gxzz*gxzz)&
-3.0*gupxy(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*(gxyz*gzzz+gxzz*gyzz)-2.0*gupxy(i,j,k)*gupzz(i,j,k)*gupzz(i,j,k)*gxzz*gzzz-&
gupyy(i,j,k)*gupxx(i,j,k)*gupxx(i,j,k)*(gxxx*gyzx+gxyx*gxzx)-2.0*gupyy(i,j,k)*gupyy(i,j,k)*gupxx(i,j,k)*gyyx*gyzx-2.0*&
gupyy(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*gyzx*gyzx-gupyy(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*(gxyx*gzzx+gxzx*gyzx)-2.0*gupyy(i,j,k)*&
gupxx(i,j,k)*gupyz(i,j,k)*(gyyx*gzzx+gyzx*gyzx)-gupyy(i,j,k)*gupxy(i,j,k)*gupxx(i,j,k)*(gxxy*gyzx+gxzx*gxyy)-2.0*&
gupyy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxzx*gyzy+gxzy*gyzx)-2.0*gupyy(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*gxyy*gyzx-2.0*&
gupyy(i,j,k)*gupyy(i,j,k)*gupxy(i,j,k)*gyzx*gyyy
      MapleGenVar4 = MapleGenVar1-gupyy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxyy*gzzx+gxzy*gyzx)-2.0*&
gupyy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gyyy*gzzx+gyzx*gyzy)-2.0*gupyy(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gyzx*gzzy+gyzy*&
gzzx)-gupyy(i,j,k)*gupxz(i,j,k)*gupxx(i,j,k)*(gxxz*gyzx+gxzx*gxyz)-2.0*gupyy(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxzx*gyzz+&
gxzz*gyzx)-2.0*gupyy(i,j,k)*gupyy(i,j,k)*gupxz(i,j,k)*gyzx*gyyz-gupyy(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxyz*gzzx+gxzz*&
gyzx)-2.0*gupyy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gyyz*gzzx+gyzx*gyzz)-2.0*gupyy(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gyzx*&
gzzz+gyzz*gzzx)-gupyy(i,j,k)*gupxy(i,j,k)*gupxx(i,j,k)*(gxxx*gyzy+gxyx*gxzy)-2.0*gupyy(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*&
gxyx*gyzy
      MapleGenVar3 = MapleGenVar4-2.0*gupyy(i,j,k)*gupyy(i,j,k)*gupxy(i,j,k)*gyyx*gyzy-gupyy(i,j,k)*gupxy(i,j,k)*&
gupxz(i,j,k)*(gxyx*gzzy+gxzx*gyzy)-2.0*gupyy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gyyx*gzzy+gyzx*gyzy)-gupyy(i,j,k)*&
gupyy(i,j,k)*gupxx(i,j,k)*(gxxy*gyzy+gxyy*gxzy)-2.0*gupyy(i,j,k)*gupyy(i,j,k)*gupxz(i,j,k)*gxzy*gyzy-2.0*gupyy(i,j,k)*&
gupyy(i,j,k)*gupxy(i,j,k)*gxyy*gyzy-gupyy(i,j,k)*gupyy(i,j,k)*gupxz(i,j,k)*(gxyy*gzzy+gxzy*gyzy)-2.0*gupyy(i,j,k)*gupyy(i,j,k)*&
gupyz(i,j,k)*(gyyy*gzzy+gyzy*gyzy)-2.0*gupyy(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*gyzy*gzzy-gupyy(i,j,k)*gupyz(i,j,k)*gupxx(i,j,k)*(&
gxxz*gyzy+gxzy*gxyz)-2.0*gupyy(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxzy*gyzz+gxzz*gyzy)-2.0*gupyy(i,j,k)*&
gupyy(i,j,k)*gupyz(i,j,k)*gyzy*gyyz-4.0*gupyy(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*gyzy*gyzz
      MapleGenVar4 = MapleGenVar3-gupyy(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxyz*gzzy+gxzz*gyzy)-2.0*&
gupyy(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gyyz*gzzy+gyzy*gyzz)-2.0*gupyy(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(gyzy*gzzz+gyzz*&
gzzy)-gupyy(i,j,k)*gupxz(i,j,k)*gupxx(i,j,k)*(gxxx*gyzz+gxyx*gxzz)-2.0*gupxy(i,j,k)*gupyy(i,j,k)*gupxx(i,j,k)*gxxy*gxzy&
-2.0*gupxy(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gxzy*gyzy-2.0*gupxy(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*gxzy*gzzy-2.0*gupxy(i,j,k)*&
gupyz(i,j,k)*gupxx(i,j,k)*gxzy*gxxz-4.0*gupxy(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*gxzy*gxzz-2.0*gupxy(i,j,k)*gupxz(i,j,k)*gupxx(i,j,k)*&
gxxx*gxzz-2.0*gupxy(i,j,k)*gupyz(i,j,k)*gupxx(i,j,k)*gxxy*gxzz-2.0*gupxy(i,j,k)*gupzz(i,j,k)*gupxx(i,j,k)*gxxz*gxzz
      MapleGenVar2 = MapleGenVar4-2.0*gupxy(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*gxzz*gyzz-2.0*gupyy(i,j,k)*&
gupxx(i,j,k)*gupxz(i,j,k)*gxzx*gyzx-2.0*gupyy(i,j,k)*gupxx(i,j,k)*gupxy(i,j,k)*gxyx*gyzx-2.0*gupyy(i,j,k)*gupxx(i,j,k)*gupzz(i,j,k)*&
gyzx*gzzx-4.0*gupyy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gyzx*gyzy-2.0*gupyy(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*gxyz*gyzx-4.0*&
gupyy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gyzx*gyzz-2.0*gupyy(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*gxyz*gyzy-2.0*gupyy(i,j,k)*gupxz(i,j,k)*&
gupxy(i,j,k)*gxyx*gyzz-2.0*gupyy(i,j,k)*gupyy(i,j,k)*gupxz(i,j,k)*gyyx*gyzz-gupyy(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxyx*gzzz+&
gxzx*gyzz)-2.0*gupyy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gyyx*gzzz+gyzx*gyzz)-gupyy(i,j,k)*gupyz(i,j,k)*gupxx(i,j,k)*(gxxy*&
gyzz+gxyy*gxzz)-2.0*gupyy(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gyyy*gyzz
      MapleGenVar4 = MapleGenVar2-gupyy(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxyy*gzzz+gxzy*gyzz)-2.0*&
gupyy(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gyyy*gzzz+gyzy*gyzz)-gupyy(i,j,k)*gupzz(i,j,k)*gupxx(i,j,k)*(gxxz*gyzz+gxyz*gxzz)&
-2.0*gupyy(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*gyyz*gyzz-2.0*gupyy(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*gyzz*gyzz-gupyy(i,j,k)*gupzz(i,j,k)*&
gupxz(i,j,k)*(gxyz*gzzz+gxzz*gyzz)-2.0*gupyy(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*(gyyz*gzzz+gyzz*gyzz)-2.0*&
gupyy(i,j,k)*gupzz(i,j,k)*gupzz(i,j,k)*gyzz*gzzz-gupyz(i,j,k)*gupxx(i,j,k)*gupxx(i,j,k)*(gxxx*gzzx+gxzx*gxzx)-4.0*gupyz(i,j,k)*&
gupyz(i,j,k)*gupxx(i,j,k)*gyzx*gzzx-2.0*gupyz(i,j,k)*gupxx(i,j,k)*gupzz(i,j,k)*gzzx*gzzx-gupyz(i,j,k)*gupxy(i,j,k)*gupxx(i,j,k)*(gxxy*&
gzzx+gxzx*gxzy)
      MapleGenVar3 = MapleGenVar4-2.0*gupyz(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxzx*gzzy+gxzy*gzzx)-&
gupyz(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*(gxyy*gzzx+gxzy*gyzx)-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*(gyzx*gzzy+gyzy*&
gzzx)-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*gyzy*gzzx-gupyz(i,j,k)*gupxz(i,j,k)*gupxx(i,j,k)*(gxxz*gzzx+gxzx*gxzz)&
-2.0*gupyz(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxzx*gzzz+gxzz*gzzx)-gupyz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*(gxyz*gzzx+gxzz*&
gyzx)-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gyzx*gzzz+gyzz*gzzx)-2.0*gupyz(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*gxzz*&
gzzx-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*gyzz*gzzx-gupyz(i,j,k)*gupxy(i,j,k)*gupxx(i,j,k)*(gxxx*gzzy+gxzx*gxzy)-&
gupyz(i,j,k)*gupxy(i,j,k)*gupxy(i,j,k)*(gxyx*gzzy+gxzx*gyzy)-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*gyzx*gzzy
      MapleGenVar4 = MapleGenVar3-gupyz(i,j,k)*gupyy(i,j,k)*gupxx(i,j,k)*(gxxy*gzzy+gxzy*gxzy)-4.0*&
gupyz(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*gyzy*gzzy-2.0*gupyz(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*gzzy*gzzy-gupyz(i,j,k)*gupyz(i,j,k)*gupxx(i,j,k)*(&
gxxz*gzzy+gxzy*gxzz)-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxzy*gzzz+gxzz*gzzy)-gupyz(i,j,k)*gupyz(i,j,k)*&
gupxy(i,j,k)*(gxyz*gzzy+gxzz*gyzy)-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gyzy*gzzz+gyzz*gzzy)-2.0*&
gupyz(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*gxzz*gzzy-4.0*gupyz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gzzy*gzzz-gupyz(i,j,k)*gupxz(i,j,k)*gupxx(i,j,k)*(&
gxxx*gzzz+gxzx*gxzz)-gupyz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*(gxyx*gzzz+gxzx*gyzz)-2.0*gupyz(i,j,k)*gupxz(i,j,k)*&
gupxz(i,j,k)*gxzx*gzzz-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*gyzx*gzzz
      CAZyz = Gamyz - (MapleGenVar4-gupyz(i,j,k)*gupyz(i,j,k)*gupxx(i,j,k)*(gxxy*gzzz+gxzy*gxzz)-gupyz(i,j,k)*&
gupyz(i,j,k)*gupxy(i,j,k)*(gxyy*gzzz+gxzy*gyzz)-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*gxzy*gzzz-gupyz(i,j,k)*gupzz(i,j,k)*&
gupxx(i,j,k)*(gxxz*gzzz+gxzz*gxzz)-4.0*gupyz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gyzz*gzzz-2.0*gupyy(i,j,k)*gupzz(i,j,k)*&
gupxz(i,j,k)*gxzz*gyzz-2.0*gupyy(i,j,k)*gupzz(i,j,k)*gupxy(i,j,k)*gxyz*gyzz-4.0*gupyz(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gxzx*gzzx&
-2.0*gupyz(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxzy*gzzx-4.0*gupyz(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*gzzx*gzzy-4.0*gupyz(i,j,k)*&
gupxz(i,j,k)*gupzz(i,j,k)*gzzx*gzzz-2.0*gupyz(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxzx*gzzy-4.0*gupyz(i,j,k)*gupyy(i,j,k)*gupxz(i,j,k)*&
gxzy*gzzy-4.0*gupyz(i,j,k)*gupzz(i,j,k)*gupxz(i,j,k)*gxzz*gzzz)
      MapleGenVar3 = -4.0*gupxz(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*gxzx*gxzz-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*&
gxzz*gxzz-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gyzz*gyzz-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*gyzy*gyzy+gupxz(i,j,k)&
*gupzz(i,j,k)*gxzzz-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupxx(i,j,k)*gxzx*gxzx+gupxy(i,j,k)*gupyz(i,j,k)*gyzxy+gupyz(i,j,k)*gupyz(i,j,k)*&
gyzyz+gupxx(i,j,k)*gupyz(i,j,k)*gxzxy+gupxy(i,j,k)*gupzz(i,j,k)*gxzyz+gupyy(i,j,k)*gupzz(i,j,k)*gyzyz+2.0*gupxz(i,j,k)*gupzz(i,j,k)*&
gzzxz+gupxy(i,j,k)*gupxz(i,j,k)*gyzxx+gupyy(i,j,k)*gupxz(i,j,k)*gyzxy+gupyz(i,j,k)*gupxz(i,j,k)*gyzxz+gupyy(i,j,k)*gupyz(i,j,k)*gyzyy&
-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxx(i,j,k)*gyzx*gyzx+gupxx(i,j,k)*gupxz(i,j,k)*gxzxx-2.0*gupzz(i,j,k)*gupzz(i,j,k)*gupxx(i,j,k)*gzzx*&
gzzx-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupyy(i,j,k)*gxzy*gxzy-4.0*gupyz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*gyzy*gyzz+gupxx(i,j,k)*&
gupzz(i,j,k)*gxzxz+gupxz(i,j,k)*gupyz(i,j,k)*gxzyz+gupyz(i,j,k)*gupzz(i,j,k)*gyzzz+gupxz(i,j,k)*gupxz(i,j,k)*gzzxx+gupxy(i,j,k)*gupyz(i,j,k)*&
gxzyy
      MapleGenVar4 = MapleGenVar3+gupzz(i,j,k)*gupzz(i,j,k)*gzzzz+gupxy(i,j,k)*gupzz(i,j,k)*gyzxz+2.0*gupyz(i,j,k)&
*gupzz(i,j,k)*gzzyz+gupyz(i,j,k)*gupyz(i,j,k)*gzzyy+gupxy(i,j,k)*gupxz(i,j,k)*gxzxy+gupxz(i,j,k)*gupxz(i,j,k)*gxzxz-2.0*gupzz(i,j,k)*&
gupzz(i,j,k)*gupyy(i,j,k)*gzzy*gzzy+2.0*gupxz(i,j,k)*gupyz(i,j,k)*gzzxy-2.0*gupxx(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gxxx*gxzx-&
gupxx(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*(gxxx*gyzx+gxyx*gxzx)-gupxx(i,j,k)*gupxx(i,j,k)*gupzz(i,j,k)*(gxxx*gzzx+gxzx*gxzx)&
-gupxx(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxxx*gyzy+gxyx*gxzy)
      MapleGenVar2 = MapleGenVar4-gupxx(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gxxx*gzzy+gxzx*gxzy)-2.0*&
gupxx(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*gxxx*gxzz-gupxx(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxxx*gyzz+gxyx*gxzz)-gupxx(i,j,k)*gupxz(i,j,k)&
*gupzz(i,j,k)*(gxxx*gzzz+gxzx*gxzz)-gupxx(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxxy*gyzx+gxzx*gxyy)-gupxx(i,j,k)*&
gupxy(i,j,k)*gupzz(i,j,k)*(gxxy*gzzx+gxzx*gxzy)-gupxx(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*(gxxy*gyzy+gxyy*gxzy)-gupxx(i,j,k)&
*gupyy(i,j,k)*gupzz(i,j,k)*(gxxy*gzzy+gxzy*gxzy)-gupxx(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxxy*gyzz+gxyy*gxzz)-&
gupxx(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(gxxy*gzzz+gxzy*gxzz)-2.0*gupxx(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*gxzx*gxxz-gupxx(i,j,k)*&
gupxz(i,j,k)*gupyz(i,j,k)*(gxxz*gyzx+gxzx*gxyz)-gupxx(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gxxz*gzzx+gxzx*gxzz)-gupxx(i,j,k)&
*gupyz(i,j,k)*gupyz(i,j,k)*(gxxz*gyzy+gxzy*gxyz)
      MapleGenVar4 = MapleGenVar2-gupxx(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(gxxz*gzzy+gxzy*gxzz)-gupxx(i,j,k)&
*gupzz(i,j,k)*gupyz(i,j,k)*(gxxz*gyzz+gxyz*gxzz)-gupxx(i,j,k)*gupzz(i,j,k)*gupzz(i,j,k)*(gxxz*gzzz+gxzz*gxzz)-&
gupxy(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*(gxxx*gyzx+gxyx*gxzx)-gupxy(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*(gxyx*gyzx+gxzx*gyyx)&
-2.0*gupxy(i,j,k)*gupxx(i,j,k)*gupzz(i,j,k)*(gxyx*gzzx+gxzx*gyzx)-gupxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxxx*gyzy+gxyx*&
gxzy)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gxyx*gyzy-gupxy(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gxyx*gzzy+gxzx*gyzy)-&
gupxy(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxxx*gyzz+gxyx*gxzz)-2.0*gupxx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxxx*gxzy-2.0*&
gupxx(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxzx*gxxy
      MapleGenVar3 = MapleGenVar4-2.0*gupxx(i,j,k)*gupyy(i,j,k)*gupxz(i,j,k)*gxxy*gxzy-2.0*gupxx(i,j,k)*&
gupyz(i,j,k)*gupxz(i,j,k)*gxxy*gxzz-2.0*gupxx(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*gxzy*gxxz-2.0*gupxx(i,j,k)*gupzz(i,j,k)*gupxz(i,j,k)*&
gxxz*gxzz-2.0*gupxy(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*gxyx*gyzx-2.0*gupxy(i,j,k)*gupxx(i,j,k)*gupxz(i,j,k)*gxyx*gxzx-2.0*&
gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gxyx*gyzz-gupxy(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gxyx*gzzz+gxzx*gyzz)-gupxy(i,j,k)*gupxy(i,j,k)&
*gupxz(i,j,k)*(gxxy*gyzx+gxzx*gxyy)-2.0*gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gxyy*gyzx-gupxy(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*&
(gxyy*gzzx+gxzy*gyzx)-gupxy(i,j,k)*gupyy(i,j,k)*gupxz(i,j,k)*(gxxy*gyzy+gxyy*gxzy)-gupxy(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)&
*(gxyy*gyzy+gxzy*gyyy)-2.0*gupxy(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*(gxyy*gzzy+gxzy*gyzy)
      MapleGenVar4 = MapleGenVar3-gupxy(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxxy*gyzz+gxyy*gxzz)-2.0*&
gupxy(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*gxyy*gyzz-gupxy(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(gxyy*gzzz+gxzy*gyzz)-gupxy(i,j,k)*gupxz(i,j,k)&
*gupxz(i,j,k)*(gxxz*gyzx+gxzx*gxyz)-gupxy(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gxyz*gzzx+gxzz*gyzx)-gupxy(i,j,k)*&
gupyz(i,j,k)*gupxz(i,j,k)*(gxxz*gyzy+gxzy*gxyz)-2.0*gupxy(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*gxyz*gyzy-gupxy(i,j,k)*gupyz(i,j,k)*&
gupzz(i,j,k)*(gxyz*gzzy+gxzz*gyzy)-gupxy(i,j,k)*gupzz(i,j,k)*gupxz(i,j,k)*(gxxz*gyzz+gxyz*gxzz)-gupxy(i,j,k)*gupzz(i,j,k)&
*gupyz(i,j,k)*(gxyz*gyzz+gxzz*gyyz)-2.0*gupxy(i,j,k)*gupzz(i,j,k)*gupzz(i,j,k)*(gxyz*gzzz+gxzz*gyzz)-gupxz(i,j,k)*&
gupxz(i,j,k)*gupxx(i,j,k)*(gxxx*gzzx+gxzx*gxzx)-gupxz(i,j,k)*gupxx(i,j,k)*gupyy(i,j,k)*(gxyx*gyzx+gxzx*gyyx)
      MapleGenVar1 = MapleGenVar4-2.0*gupxz(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*(gxyx*gzzx+gxzx*gyzx)-&
gupxz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*(gxxx*gzzy+gxzx*gxzy)-gupxz(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxyx*gzzy+gxzx*gyzy)&
-4.0*gupxz(i,j,k)*gupxz(i,j,k)*gupxy(i,j,k)*gxzx*gxzy-gupxz(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxxx*gzzz+gxzx*gxzz)-gupxz(i,j,k)*&
gupxz(i,j,k)*gupyz(i,j,k)*(gxyx*gzzz+gxzx*gyzz)-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*gxzx*gzzz-gupxz(i,j,k)*gupxz(i,j,k)*&
gupxy(i,j,k)*(gxxy*gzzx+gxzx*gxzy)-gupxz(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxyy*gzzx+gxzy*gyzx)-gupxz(i,j,k)*gupxz(i,j,k)&
*gupyy(i,j,k)*(gxxy*gzzy+gxzy*gxzy)-gupxz(i,j,k)*gupyy(i,j,k)*gupyy(i,j,k)*(gxyy*gyzy+gxzy*gyyy)-2.0*gupxz(i,j,k)*&
gupyy(i,j,k)*gupyz(i,j,k)*(gxyy*gzzy+gxzy*gyzy)-gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxxy*gzzz+gxzy*gxzz)-gupxz(i,j,k)&
*gupyz(i,j,k)*gupyz(i,j,k)*(gxyy*gzzz+gxzy*gyzz)
      MapleGenVar4 = MapleGenVar1-4.0*gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gxzy*gxzz-gupxz(i,j,k)*gupxz(i,j,k)*&
gupxz(i,j,k)*(gxxz*gzzx+gxzx*gxzz)-gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxyz*gzzx+gxzz*gyzx)-2.0*gupxy(i,j,k)*&
gupyy(i,j,k)*gupyz(i,j,k)*gxyy*gyzy-2.0*gupxy(i,j,k)*gupyy(i,j,k)*gupxz(i,j,k)*gxyy*gxzy-2.0*gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*&
gxyz*gyzx-2.0*gupxy(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*gxyz*gyzz-2.0*gupxy(i,j,k)*gupzz(i,j,k)*gupxz(i,j,k)*gxyz*gxzz-4.0*&
gupxz(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*gxzx*gyzx-6.0*gupxz(i,j,k)*gupxx(i,j,k)*gupzz(i,j,k)*gxzx*gzzx-2.0*gupxz(i,j,k)*gupxy(i,j,k)*&
gupzz(i,j,k)*gxzx*gzzy-2.0*gupxz(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*gxzy*gzzx
      MapleGenVar3 = MapleGenVar4-4.0*gupxz(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gxzy*gyzy-6.0*gupxz(i,j,k)*&
gupyy(i,j,k)*gupzz(i,j,k)*gxzy*gzzy-2.0*gupxz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gxzy*gzzz-2.0*gupxz(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*&
gxzz*gzzx-gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxxz*gzzy+gxzy*gxzz)-gupxz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxyz*gzzy+&
gxzz*gyzy)-gupxz(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gxxz*gzzz+gxzz*gxzz)-gupxz(i,j,k)*gupzz(i,j,k)*gupyy(i,j,k)*(gxyz*gyzz&
+gxzz*gyyz)-2.0*gupxz(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*(gxyz*gzzz+gxzz*gyzz)-6.0*gupxz(i,j,k)*gupzz(i,j,k)*gupzz(i,j,k)*&
gxzz*gzzz-2.0*gupxz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gxzz*gzzy-4.0*gupxz(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*gxzz*gyzz-2.0*&
gupxy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*gxyx*gxzy
      MapleGenVar4 = MapleGenVar3-gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxyx*gyzy+gxzy*gyyx)-gupxy(i,j,k)&
*gupxy(i,j,k)*gupzz(i,j,k)*(gxyx*gzzy+gxzy*gyzx)-2.0*gupxy(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*gxyx*gxzz-gupxy(i,j,k)*gupxz(i,j,k)*&
gupyz(i,j,k)*(gxyx*gyzz+gxzz*gyyx)-gupxy(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gxyx*gzzz+gxzz*gyzx)-2.0*gupxy(i,j,k)*&
gupxy(i,j,k)*gupxz(i,j,k)*gxzx*gxyy-gupxy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*(gxyy*gyzx+gxzx*gyyy)-gupxy(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)&
*(gxyy*gzzx+gxzx*gyzy)-gupxy(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxyy*gyzz+gxzz*gyyy)-gupxy(i,j,k)*gupyz(i,j,k)*&
gupzz(i,j,k)*(gxyy*gzzz+gxzz*gyzy)-2.0*gupxy(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*gxzx*gxyz-gupxy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(&
gxyz*gyzx+gxzx*gyyz)-gupxy(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gxyz*gzzx+gxzx*gyzz)
      MapleGenVar2 = MapleGenVar4-gupxy(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxyz*gyzy+gxzy*gyyz)-gupxy(i,j,k)&
*gupyz(i,j,k)*gupzz(i,j,k)*(gxyz*gzzy+gxzy*gyzz)-gupyy(i,j,k)*gupxx(i,j,k)*gupzz(i,j,k)*(gyyx*gzzx+gyzx*gyzx)-&
gupyy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxyx*gyzy+gxzy*gyyx)-gupyy(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gyyx*gzzy+gyzx*gyzy)&
-gupyy(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxyx*gyzz+gxzz*gyyx)-gupyy(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gyyx*gzzz+gyzx*gyzz&
)-gupyy(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxyy*gyzx+gxzx*gyyy)-2.0*gupxy(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*gxyy*gxzz-2.0*&
gupxy(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*gxzy*gxyz-2.0*gupyy(i,j,k)*gupxx(i,j,k)*gupyz(i,j,k)*gyyx*gyzx-2.0*gupyy(i,j,k)*gupxy(i,j,k)*&
gupyz(i,j,k)*gyyx*gyzy-2.0*gupyy(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*gyyx*gyzz-gupyy(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gyyy*gzzx+&
gyzx*gyzy)
      MapleGenVar4 = MapleGenVar2-2.0*gupyy(i,j,k)*gupyy(i,j,k)*gupyz(i,j,k)*gyyy*gyzy-gupyy(i,j,k)*gupyy(i,j,k)*&
gupzz(i,j,k)*(gyyy*gzzy+gyzy*gyzy)-gupyy(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxyy*gyzz+gxzz*gyyy)-2.0*gupyy(i,j,k)*&
gupyz(i,j,k)*gupyz(i,j,k)*gyyy*gyzz-gupyy(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(gyyy*gzzz+gyzy*gyzz)-gupyy(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)&
*(gxyz*gyzx+gxzx*gyyz)-gupyy(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gyyz*gzzx+gyzx*gyzz)-gupyy(i,j,k)*gupyz(i,j,k)*&
gupxz(i,j,k)*(gxyz*gyzy+gxzy*gyyz)-2.0*gupyy(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*gyzy*gyyz-gupyy(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(&
gyyz*gzzy+gyzy*gyzz)-gupyy(i,j,k)*gupzz(i,j,k)*gupzz(i,j,k)*(gyyz*gzzz+gyzz*gyzz)-gupyz(i,j,k)*gupyz(i,j,k)*gupxx(i,j,k)*&
(gyyx*gzzx+gyzx*gyzx)
      MapleGenVar3 = MapleGenVar4-gupyz(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxyx*gzzy+gxzy*gyzx)-gupyz(i,j,k)&
*gupyz(i,j,k)*gupxy(i,j,k)*(gyyx*gzzy+gyzx*gyzy)-4.0*gupyz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)*gyzx*gyzy-gupyz(i,j,k)*gupxz(i,j,k)*&
gupxz(i,j,k)*(gxyx*gzzz+gxzz*gyzx)-gupyz(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gyyx*gzzz+gyzx*gyzz)-4.0*gupyz(i,j,k)*&
gupyz(i,j,k)*gupxz(i,j,k)*gyzx*gyzz-gupyz(i,j,k)*gupxy(i,j,k)*gupxz(i,j,k)*(gxyy*gzzx+gxzx*gyzy)-gupyz(i,j,k)*gupyz(i,j,k)*gupxy(i,j,k)&
*(gyyy*gzzx+gyzx*gyzy)-gupyz(i,j,k)*gupyz(i,j,k)*gupyy(i,j,k)*(gyyy*gzzy+gyzy*gyzy)-gupyz(i,j,k)*gupyz(i,j,k)*&
gupxz(i,j,k)*(gxyy*gzzz+gxzz*gyzy)-gupyz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gyyy*gzzz+gyzy*gyzz)-2.0*gupyz(i,j,k)*&
gupyz(i,j,k)*gupzz(i,j,k)*gyzy*gzzz-gupyz(i,j,k)*gupxz(i,j,k)*gupxz(i,j,k)*(gxyz*gzzx+gxzx*gyzz)-gupyz(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)&
*(gyyz*gzzx+gyzx*gyzz)
      MapleGenVar4 = MapleGenVar3-2.0*gupyy(i,j,k)*gupxy(i,j,k)*gupyz(i,j,k)*gyzx*gyyy-2.0*gupyy(i,j,k)*&
gupxz(i,j,k)*gupyz(i,j,k)*gyzx*gyyz-2.0*gupyy(i,j,k)*gupzz(i,j,k)*gupyz(i,j,k)*gyyz*gyzz-6.0*gupyz(i,j,k)*gupxx(i,j,k)*gupzz(i,j,k)*&
gyzx*gzzx-2.0*gupyz(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*gyzx*gzzy-2.0*gupyz(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*gyzx*gzzz-2.0*&
gupyz(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*gyzy*gzzx-6.0*gupyz(i,j,k)*gupyy(i,j,k)*gupzz(i,j,k)*gyzy*gzzy-2.0*gupyz(i,j,k)*gupxz(i,j,k)*&
gupzz(i,j,k)*gyzz*gzzx-gupyz(i,j,k)*gupyz(i,j,k)*gupxz(i,j,k)*(gxyz*gzzy+gxzy*gyzz)-gupyz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gyyz&
*gzzy+gyzy*gyzz)-2.0*gupyz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*gyzz*gzzy-gupyz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(gyyz*gzzz+&
gyzz*gyzz)
      CAZzz = Gamzz - (MapleGenVar4-6.0*gupyz(i,j,k)*gupzz(i,j,k)*gupzz(i,j,k)*gyzz*gzzz-4.0*gupxz(i,j,k)*&
gupxy(i,j,k)*gupyz(i,j,k)*(gxzx*gyzy+gxzy*gyzx)-4.0*gupxz(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gxzx*gzzy+gxzy*gzzx)&
-4.0*gupxz(i,j,k)*gupxz(i,j,k)*gupyz(i,j,k)*(gxzx*gyzz+gxzz*gyzx)-4.0*gupxz(i,j,k)*gupxz(i,j,k)*gupzz(i,j,k)*(gxzx*gzzz+&
gxzz*gzzx)-4.0*gupxz(i,j,k)*gupyz(i,j,k)*gupyz(i,j,k)*(gxzy*gyzz+gxzz*gyzy)-4.0*gupxz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(&
gxzy*gzzz+gxzz*gzzy)-4.0*gupyz(i,j,k)*gupxy(i,j,k)*gupzz(i,j,k)*(gyzx*gzzy+gyzy*gzzx)-4.0*gupyz(i,j,k)*&
gupxz(i,j,k)*gupzz(i,j,k)*(gyzx*gzzz+gyzz*gzzx)-4.0*gupyz(i,j,k)*gupyz(i,j,k)*gupzz(i,j,k)*(gyzy*gzzz+gyzz*gzzy)&
-4.0*gupzz(i,j,k)*gupzz(i,j,k)*gupxy(i,j,k)*gzzx*gzzy-4.0*gupzz(i,j,k)*gupzz(i,j,k)*gupxz(i,j,k)*gzzx*gzzz-4.0*gupzz(i,j,k)*&
gupzz(i,j,k)*gupyz(i,j,k)*gzzy*gzzz-2.0*gupzz(i,j,k)*gupzz(i,j,k)*gupzz(i,j,k)*gzzz*gzzz)                                

! second kind of connection
  Gamxxx =HALF*( gupxx(i,j,k)*gxxx + gupxy(i,j,k)*(TWO*gxyx - gxxy ) + gupxz(i,j,k)*(TWO*gxzx - gxxz ))
  Gamyxx =HALF*( gupxy(i,j,k)*gxxx + gupyy(i,j,k)*(TWO*gxyx - gxxy ) + gupyz(i,j,k)*(TWO*gxzx - gxxz ))
  Gamzxx =HALF*( gupxz(i,j,k)*gxxx + gupyz(i,j,k)*(TWO*gxyx - gxxy ) + gupzz(i,j,k)*(TWO*gxzx - gxxz ))
 
  Gamxyy =HALF*( gupxx(i,j,k)*(TWO*gxyy - gyyx ) + gupxy(i,j,k)*gyyy + gupxz(i,j,k)*(TWO*gyzy - gyyz ))
  Gamyyy =HALF*( gupxy(i,j,k)*(TWO*gxyy - gyyx ) + gupyy(i,j,k)*gyyy + gupyz(i,j,k)*(TWO*gyzy - gyyz ))
  Gamzyy =HALF*( gupxz(i,j,k)*(TWO*gxyy - gyyx ) + gupyz(i,j,k)*gyyy + gupzz(i,j,k)*(TWO*gyzy - gyyz ))

  Gamxzz =HALF*( gupxx(i,j,k)*(TWO*gxzz - gzzx ) + gupxy(i,j,k)*(TWO*gyzz - gzzy ) + gupxz(i,j,k)*gzzz)
  Gamyzz =HALF*( gupxy(i,j,k)*(TWO*gxzz - gzzx ) + gupyy(i,j,k)*(TWO*gyzz - gzzy ) + gupyz(i,j,k)*gzzz)
  Gamzzz =HALF*( gupxz(i,j,k)*(TWO*gxzz - gzzx ) + gupyz(i,j,k)*(TWO*gyzz - gzzy ) + gupzz(i,j,k)*gzzz)

  Gamxxy =HALF*( gupxx(i,j,k)*gxxy + gupxy(i,j,k)*gyyx + gupxz(i,j,k)*( gxzy + gyzx - gxyz ) )
  Gamyxy =HALF*( gupxy(i,j,k)*gxxy + gupyy(i,j,k)*gyyx + gupyz(i,j,k)*( gxzy + gyzx - gxyz ) )
  Gamzxy =HALF*( gupxz(i,j,k)*gxxy + gupyz(i,j,k)*gyyx + gupzz(i,j,k)*( gxzy + gyzx - gxyz ) )

  Gamxxz =HALF*( gupxx(i,j,k)*gxxz + gupxy(i,j,k)*( gxyz + gyzx - gxzy ) + gupxz(i,j,k)*gzzx )
  Gamyxz =HALF*( gupxy(i,j,k)*gxxz + gupyy(i,j,k)*( gxyz + gyzx - gxzy ) + gupyz(i,j,k)*gzzx )
  Gamzxz =HALF*( gupxz(i,j,k)*gxxz + gupyz(i,j,k)*( gxyz + gyzx - gxzy ) + gupzz(i,j,k)*gzzx )

  Gamxyz =HALF*( gupxx(i,j,k)*( gxyz + gxzy - gyzx ) + gupxy(i,j,k)*gyyz + gupxz(i,j,k)*gzzy )
  Gamyyz =HALF*( gupxy(i,j,k)*( gxyz + gxzy - gyzx ) + gupyy(i,j,k)*gyyz + gupyz(i,j,k)*gzzy )
  Gamzyz =HALF*( gupxz(i,j,k)*( gxyz + gxzy - gyzx ) + gupyz(i,j,k)*gyyz + gupzz(i,j,k)*gzzy )

  Gamxa =       gupxx(i,j,k) * Gamxxx + gupyy(i,j,k) * Gamxyy + gupzz(i,j,k) * Gamxzz + &
          TWO*( gupxy(i,j,k) * Gamxxy + gupxz(i,j,k) * Gamxxz + gupyz(i,j,k) * Gamxyz )
  Gamya =       gupxx(i,j,k) * Gamyxx + gupyy(i,j,k) * Gamyyy + gupzz(i,j,k) * Gamyzz + &
          TWO*( gupxy(i,j,k) * Gamyxy + gupxz(i,j,k) * Gamyxz + gupyz(i,j,k) * Gamyyz )
  Gamza =       gupxx(i,j,k) * Gamzxx + gupyy(i,j,k) * Gamzyy + gupzz(i,j,k) * Gamzzz + &
          TWO*( gupxy(i,j,k) * Gamzxy + gupxz(i,j,k) * Gamzxz + gupyz(i,j,k) * Gamzyz )

         call point_fdderivs_shc(ex,Lap,Lapxx,Lapxy,Lapxz,Lapyy,Lapyz,Lapzz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst,&
                                 drhodx, drhody, drhodz,                                                    &
                                 dsigmadx,dsigmady,dsigmadz,                                                &
                                 dRdx,dRdy,dRdz,                                                            &
                                 drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                                 dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                                 dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)     
         call point_fdderivs_shc(ex,betax,sfxxx,sfxxy,sfxxz,sfxyy,sfxyz,sfxzz,crho,sigma,R,ANTI,SYM,SYM,Symmetry,0,sst,&
                                 drhodx, drhody, drhodz,                                                    &
                                 dsigmadx,dsigmady,dsigmadz,                                                &
                                 dRdx,dRdy,dRdz,                                                            &
                                 drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                                 dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                                 dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)        
         call point_fdderivs_shc(ex,betay,sfyxx,sfyxy,sfyxz,sfyyy,sfyyz,sfyzz,crho,sigma,R,SYM,ANTI,SYM,Symmetry,0,sst,&
                                 drhodx, drhody, drhodz,                                                    &
                                 dsigmadx,dsigmady,dsigmadz,                                                &
                                 dRdx,dRdy,dRdz,                                                            &
                                 drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                                 dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                                 dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)     
         call point_fdderivs_shc(ex,betaz,sfzxx,sfzxy,sfzxz,sfzyy,sfzyz,sfzzz,crho,sigma,R,SYM,SYM,ANTI,Symmetry,0,sst,&
                                 drhodx, drhody, drhodz,                                                    &
                                 dsigmadx,dsigmady,dsigmadz,                                                &
                                 dRdx,dRdy,dRdz,                                                            &
                                 drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                                 dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                                 dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)

         AAxx =       gupxx(i,j,k) * Axx(i,j,k) * Axx(i,j,k) + gupyy(i,j,k) * Axy(i,j,k) * Axy(i,j,k) + gupzz(i,j,k) * Axz(i,j,k) * Axz(i,j,k) + &
               TWO * (gupxy(i,j,k) * Axx(i,j,k) * Axy(i,j,k) + gupxz(i,j,k) * Axx(i,j,k) * Axz(i,j,k) + gupyz(i,j,k) * Axy(i,j,k) * Axz(i,j,k))
         AAyy =       gupxx(i,j,k) * Axy(i,j,k) * Axy(i,j,k) + gupyy(i,j,k) * Ayy(i,j,k) * Ayy(i,j,k) + gupzz(i,j,k) * Ayz(i,j,k) * Ayz(i,j,k) + &
               TWO * (gupxy(i,j,k) * Axy(i,j,k) * Ayy(i,j,k) + gupxz(i,j,k) * Axy(i,j,k) * Ayz(i,j,k) + gupyz(i,j,k) * Ayy(i,j,k) * Ayz(i,j,k))
         AAzz =       gupxx(i,j,k) * Axz(i,j,k) * Axz(i,j,k) + gupyy(i,j,k) * Ayz(i,j,k) * Ayz(i,j,k) + gupzz(i,j,k) * Azz(i,j,k) * Azz(i,j,k) + &
               TWO * (gupxy(i,j,k) * Axz(i,j,k) * Ayz(i,j,k) + gupxz(i,j,k) * Axz(i,j,k) * Azz(i,j,k) + gupyz(i,j,k) * Ayz(i,j,k) * Azz(i,j,k))
         AAxy =       gupxx(i,j,k) * Axx(i,j,k) * Axy(i,j,k) + gupyy(i,j,k) * Axy(i,j,k) * Ayy(i,j,k) + gupzz(i,j,k) * Axz(i,j,k) * Ayz(i,j,k) + &
                      gupxy(i,j,k) *(Axx(i,j,k) * Ayy(i,j,k) + Axy(i,j,k) * Axy(i,j,k))                            + &
                      gupxz(i,j,k) *(Axx(i,j,k) * Ayz(i,j,k) + Axz(i,j,k) * Axy(i,j,k))                            + &
                      gupyz(i,j,k) *(Axy(i,j,k) * Ayz(i,j,k) + Axz(i,j,k) * Ayy(i,j,k))
         AAxz =       gupxx(i,j,k) * Axx(i,j,k) * Axz(i,j,k) + gupyy(i,j,k) * Axy(i,j,k) * Ayz(i,j,k) + gupzz(i,j,k) * Axz(i,j,k) * Azz(i,j,k) + &
                      gupxy(i,j,k) *(Axx(i,j,k) * Ayz(i,j,k) + Axy(i,j,k) * Axz(i,j,k))                            + &
                      gupxz(i,j,k) *(Axx(i,j,k) * Azz(i,j,k) + Axz(i,j,k) * Axz(i,j,k))                            + &
                      gupyz(i,j,k) *(Axy(i,j,k) * Azz(i,j,k) + Axz(i,j,k) * Ayz(i,j,k))
         AAyz =       gupxx(i,j,k) * Axy(i,j,k) * Axz(i,j,k) + gupyy(i,j,k) * Ayy(i,j,k) * Ayz(i,j,k) + gupzz(i,j,k) * Ayz(i,j,k) * Azz(i,j,k) + &
                      gupxy(i,j,k) *(Axy(i,j,k) * Ayz(i,j,k) + Ayy(i,j,k) * Axz(i,j,k))                            + &
                      gupxz(i,j,k) *(Axy(i,j,k) * Azz(i,j,k) + Ayz(i,j,k) * Axz(i,j,k))                            + &
                      gupyz(i,j,k) *(Ayy(i,j,k) * Azz(i,j,k) + Ayz(i,j,k) * Ayz(i,j,k))                                 

         betas = betax(i,j,k)*slx(i,j,k)+betay(i,j,k)*sly(i,j,k)+betaz(i,j,k)*slz(i,j,k)
         fxx = trK(i,j,k)+TWO*TZ(i,j,k)
         fxy = fxx*Axy(i,j,k)-TWO*AAxy
         fxz = fxx*Axz(i,j,k)-TWO*AAxz
         fyy = fxx*Ayy(i,j,k)-TWO*AAyy
         fyz = fxx*Ayz(i,j,k)-TWO*AAyz
         fzz = fxx*Azz(i,j,k)-TWO*AAzz
         fxx = fxx*Axx(i,j,k)-TWO*AAxx

         muL = 2.d0/alpn1(i,j,k)
         tmuSL = chin1(i,j,k)*2.d0/dsqrt(3.d0)/alpn1(i,j,k)**2
         tmuST = chin1(i,j,k)/alpn1(i,j,k)**2
! Eq.(17)
         totrK_rhs =  (betax(i,j,k)*Kx+betay(i,j,k)*Ky+betaz(i,j,k)*Kz) &
                     -dsqrt(muL)*alpn1(i,j,k)*(vx(i,j,k)*Kx+vy(i,j,k)*Ky+vz(i,j,k)*Kz+trK(i,j,k)/R(k))  
#if 0                     
                     -0.5d0*(qupxx(i,j,k)*Lapxx+qupyy(i,j,k)*Lapyy+qupzz(i,j,k)*Lapzz+ &
                     TWO*(qupxy(i,j,k)*Lapxy+qupxz(i,j,k)*Lapxz+qupyz(i,j,k)*Lapyz)) &
                     -trK(i,j,k)/R(k)*betas                                                             &
                     -0.5d0*alpn1(i,j,k)*(gupxx(i,j,k)*AAxx+gupyy(i,j,k)*AAyy+gupzz(i,j,k)*AAzz+ &
                     TWO*(gupxy(i,j,k)*AAxy+gupxz(i,j,k)*AAxz+gupyz(i,j,k)*AAyz)+(trK(i,j,k)+TWO*TZ(i,j,k))**2/3.d0 &
                     +kappa1*(ONE-kappa2)*TZ(i,j,k))+(ONE+betas/dsqrt(muL)/alpn1(i,j,k))/R(k)*(vx(i,j,k)*Lapx+vy(i,j,k)*Lapy+vz(i,j,k)*Lapz) &
                     +ha/R(k)**4-kappa3*alpn1(i,j,k)*trK(i,j,k)
#endif                     

! Eq.(18)                     
         toGams_rhs = -alpn1(i,j,k)*dsqrt(tmuSL)*(Gamxx+Gamyy+Gamzz) + ( &
                     slx(i,j,k)*(qupxx(i,j,k)*sfxxx+qupyy(i,j,k)*sfxyy+qupzz(i,j,k)*sfxzz+TWO*(qupxy(i,j,k)*sfxxy+qupxz(i,j,k)*sfxxz+qupyz(i,j,k)*sfxyz)) &
                    +sly(i,j,k)*(qupxx(i,j,k)*sfyxx+qupyy(i,j,k)*sfyyy+qupzz(i,j,k)*sfyzz+TWO*(qupxy(i,j,k)*sfyxy+qupxz(i,j,k)*sfyxz+qupyz(i,j,k)*sfyyz)) &
                    +slz(i,j,k)*(qupxx(i,j,k)*sfzxx+qupyy(i,j,k)*sfzyy+qupzz(i,j,k)*sfzzz+TWO*(qupxy(i,j,k)*sfzxy+qupxz(i,j,k)*sfzxz+qupyz(i,j,k)*sfzyz)) ) &
                    /chin1(i,j,k)                                    - ( &
                        vx(i,j,k)*(qulxx(i,j,k)*sfxxx+qulxy(i,j,k)*sfxxy+qulxz(i,j,k)*sfxxz+ &
                                   qulxy(i,j,k)*sfyxx+qulyy(i,j,k)*sfyxy+qulyz(i,j,k)*sfyxz+ &
                                   qulxz(i,j,k)*sfzxx+qulyz(i,j,k)*sfzxy+qulzz(i,j,k)*sfzxz) &
                       +vy(i,j,k)*(qulxx(i,j,k)*sfxxy+qulxy(i,j,k)*sfxyy+qulxz(i,j,k)*sfxyz+ &
                                   qulxy(i,j,k)*sfyxy+qulyy(i,j,k)*sfyyy+qulyz(i,j,k)*sfyyz+ &
                                   qulxz(i,j,k)*sfzxy+qulyz(i,j,k)*sfzyy+qulzz(i,j,k)*sfzyz) &
                       +vz(i,j,k)*(qulxx(i,j,k)*sfxxz+qulxy(i,j,k)*sfxyz+qulxz(i,j,k)*sfxzz+ &
                                   qulxy(i,j,k)*sfyxz+qulyy(i,j,k)*sfyyz+qulyz(i,j,k)*sfyzz+ &
                                   qulxz(i,j,k)*sfzxz+qulyz(i,j,k)*sfzyz+qulzz(i,j,k)*sfzzz) )/chin1(i,j,k) &
                      -4.d0*alpn1(i,j,k)*dsqrt(muL)/3.d0/(dsqrt(tmuSL)+dsqrt(muL))/chin1(i,j,k)*(vx(i,j,k)*Kx+vy(i,j,k)*Ky+vz(i,j,k)*Kz) &
                      -2.d0*alpn1(i,j,k)/3.d0/(dsqrt(tmuSL)+ONE)/chin1(i,j,k)*(vx(i,j,k)*TZx+vy(i,j,k)*TZy+vz(i,j,k)*TZz) &
                      +thbs-kappa3*alpn1(i,j,k)*(slx(i,j,k)*Gamx(i,j,k)+sly(i,j,k)*Gamy(i,j,k)+slz(i,j,k)*Gamz(i,j,k)) &
                     +(slx(i,j,k)*(betax(i,j,k)*Gamxx+betay(i,j,k)*Gamxy+betaz(i,j,k)*Gamxz)      &
                      +sly(i,j,k)*(betax(i,j,k)*Gamyx+betay(i,j,k)*Gamyy+betaz(i,j,k)*Gamyz)      &
                      +slz(i,j,k)*(betax(i,j,k)*Gamzx+betay(i,j,k)*Gamzy+betaz(i,j,k)*Gamzz))

         toTZ_rhs = -alpn1(i,j,k)*(vx(i,j,k)*TZx+vy(i,j,k)*TZy+vz(i,j,k)*TZz)+(betax(i,j,k)*TZx+betay(i,j,k)*TZy+betaz(i,j,k)*TZz)

         toAss_rhs = -alpn1(i,j,k)*chin1(i,j,k)*(                                &
                    TWO*((gupxx(i,j,k)*(Axxx-(Gamxxx*Axx(i,j,k)+Gamyxx*Axy(i,j,k)+Gamzxx*Axz(i,j,k)))  &
                       +  gupxy(i,j,k)*(Axxy-(Gamxxx*Axy(i,j,k)+Gamyxx*Ayy(i,j,k)+Gamzxx*Ayz(i,j,k)))  &
                       +  gupxz(i,j,k)*(Axxz-(Gamxxx*Axz(i,j,k)+Gamyxx*Ayz(i,j,k)+Gamzxx*Azz(i,j,k)))  & 
                       +  gupxy(i,j,k)*(Axxy-(Gamxxy*Axx(i,j,k)+Gamyxy*Axy(i,j,k)+Gamzxy*Axz(i,j,k)))  &
                       +  gupyy(i,j,k)*(Axyy-(Gamxxy*Axy(i,j,k)+Gamyxy*Ayy(i,j,k)+Gamzxy*Ayz(i,j,k)))  &
                       +  gupyz(i,j,k)*(Axyz-(Gamxxy*Axz(i,j,k)+Gamyxy*Ayz(i,j,k)+Gamzxy*Azz(i,j,k)))  &
                       +  gupxz(i,j,k)*(Axxz-(Gamxxz*Axx(i,j,k)+Gamyxz*Axy(i,j,k)+Gamzxz*Axz(i,j,k)))  &
                       +  gupyz(i,j,k)*(Axyz-(Gamxxz*Axy(i,j,k)+Gamyxz*Ayy(i,j,k)+Gamzxz*Ayz(i,j,k)))  &
                       +  gupzz(i,j,k)*(Axzz-(Gamxxz*Axz(i,j,k)+Gamyxz*Ayz(i,j,k)+Gamzxz*Azz(i,j,k)))  &
                       -  (Gamxa*Axx(i,j,k)+Gamya*Axy(i,j,k)+Gamza*Axz(i,j,k)) )*vx(i,j,k)             &
                       + (gupxx(i,j,k)*(Axyx-(Gamxxy*Axx(i,j,k)+Gamyxy*Axy(i,j,k)+Gamzxy*Axz(i,j,k)))  &
                       +  gupxy(i,j,k)*(Axyy-(Gamxxy*Axy(i,j,k)+Gamyxy*Ayy(i,j,k)+Gamzxy*Ayz(i,j,k)))  &
                       +  gupxz(i,j,k)*(Axyz-(Gamxxy*Axz(i,j,k)+Gamyxy*Ayz(i,j,k)+Gamzxy*Azz(i,j,k)))  & 
                       +  gupxy(i,j,k)*(Axyy-(Gamxyy*Axx(i,j,k)+Gamyyy*Axy(i,j,k)+Gamzyy*Axz(i,j,k)))  &
                       +  gupyy(i,j,k)*(Ayyy-(Gamxyy*Axy(i,j,k)+Gamyyy*Ayy(i,j,k)+Gamzyy*Ayz(i,j,k)))  &
                       +  gupyz(i,j,k)*(Ayyz-(Gamxyy*Axz(i,j,k)+Gamyyy*Ayz(i,j,k)+Gamzyy*Azz(i,j,k)))  &
                       +  gupxz(i,j,k)*(Axyz-(Gamxyz*Axx(i,j,k)+Gamyyz*Axy(i,j,k)+Gamzyz*Axz(i,j,k)))  &
                       +  gupyz(i,j,k)*(Ayyz-(Gamxyz*Axy(i,j,k)+Gamyyz*Ayy(i,j,k)+Gamzyz*Ayz(i,j,k)))  &
                       +  gupzz(i,j,k)*(Ayzz-(Gamxyz*Axz(i,j,k)+Gamyyz*Ayz(i,j,k)+Gamzyz*Azz(i,j,k)))  &
                       -  (Gamxa*Axy(i,j,k)+Gamya*Ayy(i,j,k)+Gamza*Ayz(i,j,k)) )*vy(i,j,k)             &
                       + (gupxx(i,j,k)*(Axzx-(Gamxxz*Axx(i,j,k)+Gamyxz*Axy(i,j,k)+Gamzxz*Axz(i,j,k)))  &
                       +  gupxy(i,j,k)*(Axzy-(Gamxxz*Axy(i,j,k)+Gamyxz*Ayy(i,j,k)+Gamzxz*Ayz(i,j,k)))  &
                       +  gupxz(i,j,k)*(Axzz-(Gamxxz*Axz(i,j,k)+Gamyxz*Ayz(i,j,k)+Gamzxz*Azz(i,j,k)))  & 
                       +  gupxy(i,j,k)*(Axzy-(Gamxyz*Axx(i,j,k)+Gamyyz*Axy(i,j,k)+Gamzyz*Axz(i,j,k)))  &
                       +  gupyy(i,j,k)*(Ayzy-(Gamxyz*Axy(i,j,k)+Gamyyz*Ayy(i,j,k)+Gamzyz*Ayz(i,j,k)))  &
                       +  gupyz(i,j,k)*(Ayzz-(Gamxyz*Axz(i,j,k)+Gamyyz*Ayz(i,j,k)+Gamzyz*Azz(i,j,k)))  &
                       +  gupxz(i,j,k)*(Axzz-(Gamxzz*Axx(i,j,k)+Gamyzz*Axy(i,j,k)+Gamzzz*Axz(i,j,k)))  &
                       +  gupyz(i,j,k)*(Ayzz-(Gamxzz*Axy(i,j,k)+Gamyzz*Ayy(i,j,k)+Gamzzz*Ayz(i,j,k)))  &
                       +  gupzz(i,j,k)*(Azzz-(Gamxzz*Axz(i,j,k)+Gamyzz*Ayz(i,j,k)+Gamzzz*Azz(i,j,k)))  &
                       -  (Gamxa*Axz(i,j,k)+Gamya*Ayz(i,j,k)+Gamza*Azz(i,j,k)) )*vz(i,j,k)          )  &
                   -2.d0/3.d0*chin1(i,j,k)*(vx(i,j,k)*(TWO*Kx+TZx)+vy(i,j,k)*(TWO*Ky+TZy)+vz(i,j,k)*(TWO*Kz+TZz)) &
                   -2.d0/3.d0*(Rxx(i,j,k)*vx(i,j,k)*vx(i,j,k)+Ryy(i,j,k)*vy(i,j,k)*vy(i,j,k)+Rzz(i,j,k)*vz(i,j,k)*vz(i,j,k) &
                         +TWO*(Rxy(i,j,k)*vx(i,j,k)*vy(i,j,k)+Rxz(i,j,k)*vx(i,j,k)*vz(i,j,k)+Ryz(i,j,k)*vy(i,j,k)*vz(i,j,k))) &
                   + ONE/3.d0*(Rxx(i,j,k)*qupxx(i,j,k)+Ryy(i,j,k)*qupyy(i,j,k)+Rzz(i,j,k)*qupzz(i,j,k) &
                         +TWO*(Rxy(i,j,k)*qupxy(i,j,k)+Rxz(i,j,k)*qupxz(i,j,k)+Ryz(i,j,k)*qupyz(i,j,k))) &
                   +2.d0/3.d0*chin1(i,j,k)*(slx(i,j,k)*vx(i,j,k)*CAZxx+slx(i,j,k)*vy(i,j,k)*CAZxy+slx(i,j,k)*vz(i,j,k)*CAZxz &
                                           +sly(i,j,k)*vx(i,j,k)*CAZyx+sly(i,j,k)*vy(i,j,k)*CAZyy+sly(i,j,k)*vz(i,j,k)*CAZyz &
                                           +slz(i,j,k)*vx(i,j,k)*CAZzx+slz(i,j,k)*vy(i,j,k)*CAZzy+slz(i,j,k)*vz(i,j,k)*CAZzz) &
                   -ONE/3.d0*chin1(i,j,k)*(qulxx(i,j,k)*CAZxx+qulyx(i,j,k)*CAZxy+qulzx(i,j,k)*CAZxz &
                                          +qulxy(i,j,k)*CAZyx+qulyy(i,j,k)*CAZyy+qulzy(i,j,k)*CAZyz &
                                          +qulxz(i,j,k)*CAZzx+qulyz(i,j,k)*CAZzy+qulzz(i,j,k)*CAZzz) &
                   -3.d0/chin1(i,j,k)*(vx(i,j,k)*(gupxx(i,j,k)*chix*Axx(i,j,k)+gupxy(i,j,k)*chix*Axy(i,j,k)+gupxz(i,j,k)*chix*Axz(i,j,k) &
                                           +gupxy(i,j,k)*chiy*Axx(i,j,k)+gupyy(i,j,k)*chiy*Axy(i,j,k)+gupyz(i,j,k)*chiy*Axz(i,j,k) &
                                           +gupxz(i,j,k)*chiz*Axx(i,j,k)+gupyz(i,j,k)*chiz*Axy(i,j,k)+gupzz(i,j,k)*chiz*Axz(i,j,k)) &
                                       +vy(i,j,k)*(gupxx(i,j,k)*chix*Axy(i,j,k)+gupxy(i,j,k)*chix*Ayy(i,j,k)+gupxz(i,j,k)*chix*Ayz(i,j,k) &
                                           +gupxy(i,j,k)*chiy*Axy(i,j,k)+gupyy(i,j,k)*chiy*Ayy(i,j,k)+gupyz(i,j,k)*chiy*Ayz(i,j,k) &
                                           +gupxz(i,j,k)*chiz*Axy(i,j,k)+gupyz(i,j,k)*chiz*Ayy(i,j,k)+gupzz(i,j,k)*chiz*Ayz(i,j,k)) &
                                       +vz(i,j,k)*(gupxx(i,j,k)*chix*Axz(i,j,k)+gupxy(i,j,k)*chix*Ayz(i,j,k)+gupxz(i,j,k)*chix*Azz(i,j,k) &
                                           +gupxy(i,j,k)*chiy*Axz(i,j,k)+gupyy(i,j,k)*chiy*Ayz(i,j,k)+gupyz(i,j,k)*chiy*Azz(i,j,k) &
                                           +gupxz(i,j,k)*chiz*Axz(i,j,k)+gupyz(i,j,k)*chiz*Ayz(i,j,k)+gupzz(i,j,k)*chiz*Azz(i,j,k)) ) &
                   -kappa1*(vx(i,j,k)*Gmxcon(i,j,k)+vy(i,j,k)*Gmycon(i,j,k)+vz(i,j,k)*Gmzcon(i,j,k))  )          &
                   +alpn1(i,j,k)*(fxx*vx(i,j,k)*vx(i,j,k)+fyy*vy(i,j,k)*vy(i,j,k)+fzz*vz(i,j,k)*vz(i,j,k) &
                            +TWO*(fxy*vx(i,j,k)*vy(i,j,k)+fxz*vx(i,j,k)*vz(i,j,k)+fyz*vy(i,j,k)*vz(i,j,k)))

! Eq.(22)                            
  toAs1_rhs = alpn1(i,j,k)*(fxx*vx(i,j,k)*ux(i,j,k)+fxy*vy(i,j,k)*ux(i,j,k)+fxz*vz(i,j,k)*ux(i,j,k) &
                           +fxy*vx(i,j,k)*uy(i,j,k)+fyy*vy(i,j,k)*uy(i,j,k)+fyz*vz(i,j,k)*uy(i,j,k) &
                           +fxz*vx(i,j,k)*uz(i,j,k)+fyz*vy(i,j,k)*uz(i,j,k)+fzz*vz(i,j,k)*uz(i,j,k))

  toAs2_rhs = alpn1(i,j,k)*(fxx*vx(i,j,k)*wx(i,j,k)+fxy*vy(i,j,k)*wx(i,j,k)+fxz*vz(i,j,k)*wx(i,j,k) &
                           +fxy*vx(i,j,k)*wy(i,j,k)+fyy*vy(i,j,k)*wy(i,j,k)+fyz*vz(i,j,k)*wy(i,j,k) &
                           +fxz*vx(i,j,k)*wz(i,j,k)+fyz*vy(i,j,k)*wz(i,j,k)+fzz*vz(i,j,k)*wz(i,j,k))

  fxx = Lapxx - (Gamxxx-((chix+chix)/chin1(i,j,k)-gxx(i,j,k)*gxxx)*HALF)*Lapx - (Gamyxx+gxx(i,j,k)*gxxy*HALF)*Lapy - (Gamzxx+gxx(i,j,k)*gxxz*HALF)*Lapz
  fyy = Lapyy - (Gamxyy+gyy(i,j,k)*gxxx*HALF)*Lapx - (Gamyyy-((chiy+chiy)/chin1(i,j,k)-gyy(i,j,k)*gxxy)*HALF)*Lapy - (Gamzyy+gyy(i,j,k)*gxxz*HALF)*Lapz
  fzz = Lapzz - (Gamxzz+gzz(i,j,k)*gxxx*HALF)*Lapx - (Gamyzz+gzz(i,j,k)*gxxy*HALF)*Lapy - (Gamzzz-((chiz+chiz)/chin1(i,j,k)-gzz(i,j,k)*gxxz)*HALF)*Lapz
  fxy = Lapxy - (Gamxxy-(chiy/chin1(i,j,k)-gxy(i,j,k)*gxxx)*HALF)*Lapx - (Gamyxy-(chix/chin1(i,j,k)-gxy(i,j,k)*gxxy)*HALF)*Lapy&
              - (Gamzxy+gxy(i,j,k)*gxxz*HALF)*Lapz
  fxz = Lapxz - (Gamxxz-(chiz/chin1(i,j,k)-gxz(i,j,k)*gxxx)*HALF)*Lapx - (Gamyxz+gxz(i,j,k)*gxxy*HALF)*Lapy&
              - (Gamzxz-(chix/chin1(i,j,k)-gxz(i,j,k)*gxxz)*HALF)*Lapz
  fyz = Lapyz - (Gamxyz+gyz(i,j,k)*gxxx*HALF)*Lapx - (Gamyyz-(chiz/chin1(i,j,k)-gyz(i,j,k)*gxxy)*HALF)*Lapy&
              - (Gamzyz-(chiy/chin1(i,j,k)-gyz(i,j,k)*gxxz)*HALF)*Lapz

  TFxx = -chin1(i,j,k)*fxx
  TFxy = -chin1(i,j,k)*fxy
  TFxz = -chin1(i,j,k)*fxz
  TFyy = -chin1(i,j,k)*fyy
  TFyz = -chin1(i,j,k)*fyz
  TFzz = -chin1(i,j,k)*fzz
  toAss_rhs = toAss_rhs -2.d0/3.d0*chin1(i,j,k)*(fxx*vx(i,j,k)*vx(i,j,k)+fyy*vy(i,j,k)*vy(i,j,k)+fzz*vz(i,j,k)*vz(i,j,k) &
                                            +TWO*(fxy*vx(i,j,k)*vy(i,j,k)+fxz*vx(i,j,k)*vz(i,j,k)+fyz*vy(i,j,k)*vz(i,j,k))) &
              +ONE/3.d0*chin1(i,j,k)*(fxx*qupxx(i,j,k)+fyy*qupyy(i,j,k)+fzz*qupzz(i,j,k) &
                                 +TWO*(fxy*qupxy(i,j,k)+fxz*qupxz(i,j,k)+fyz*qupyz(i,j,k)))
  toAs1_rhs =   toAs1_rhs -chin1(i,j,k)*(fxx*vx(i,j,k)*ux(i,j,k)+fxy*vy(i,j,k)*ux(i,j,k)+fxz*vz(i,j,k)*ux(i,j,k) &
                                        +fxy*vx(i,j,k)*uy(i,j,k)+fyy*vy(i,j,k)*uy(i,j,k)+fyz*vz(i,j,k)*uy(i,j,k) &
                                        +fxz*vx(i,j,k)*uz(i,j,k)+fyz*vy(i,j,k)*uz(i,j,k)+fzz*vz(i,j,k)*uz(i,j,k))
  toAs2_rhs =   toAs2_rhs -chin1(i,j,k)*(fxx*vx(i,j,k)*wx(i,j,k)+fxy*vy(i,j,k)*wx(i,j,k)+fxz*vz(i,j,k)*wx(i,j,k) &
                                        +fxy*vx(i,j,k)*wy(i,j,k)+fyy*vy(i,j,k)*wy(i,j,k)+fyz*vz(i,j,k)*wy(i,j,k) &
                                        +fxz*vx(i,j,k)*wz(i,j,k)+fyz*vy(i,j,k)*wz(i,j,k)+fzz*vz(i,j,k)*wz(i,j,k))

  fxx = (betax(i,j,k)*Axxx+betay(i,j,k)*Axxy+betaz(i,j,k)*Axxz)-TWO*(Axx(i,j,k)*sfxx+Axy(i,j,k)*sfyx+Axz(i,j,k)*sfzx)
  fxy = (betax(i,j,k)*Axyx+betay(i,j,k)*Axyy+betaz(i,j,k)*Axyz)- &
        (Axx(i,j,k)*sfxy+Axy(i,j,k)*sfyy+Axz(i,j,k)*sfzy)-(Axy(i,j,k)*sfxx+Ayy(i,j,k)*sfyx+Ayz(i,j,k)*sfzx)
  fxz = (betax(i,j,k)*Axzx+betay(i,j,k)*Axzy+betaz(i,j,k)*Axzz)- &
        (Axx(i,j,k)*sfxz+Axy(i,j,k)*sfyz+Axz(i,j,k)*sfzz)-(Axz(i,j,k)*sfxx+Ayz(i,j,k)*sfyx+Azz(i,j,k)*sfzx)
  fyy = (betax(i,j,k)*Ayyx+betay(i,j,k)*Ayyy+betaz(i,j,k)*Ayyz)-TWO*(Axy(i,j,k)*sfxy+Ayy(i,j,k)*sfyy+Ayz(i,j,k)*sfzy)
  fyz = (betax(i,j,k)*Ayzx+betay(i,j,k)*Ayzy+betaz(i,j,k)*Ayzz)- &
        (Axy(i,j,k)*sfxz+Ayy(i,j,k)*sfyz+Ayz(i,j,k)*sfzz)-(Axz(i,j,k)*sfxy+Ayz(i,j,k)*sfyy+Azz(i,j,k)*sfzy)
  fzz = (betax(i,j,k)*Azzx+betay(i,j,k)*Azzy+betaz(i,j,k)*Azzz)-TWO*(Axz(i,j,k)*sfxz+Ayz(i,j,k)*sfyz+Azz(i,j,k)*sfzz)
  TFxx = TFxx+fxx
  TFxy = TFxy+fxy
  TFxz = TFxz+fxz
  TFyy = TFyy+fyy
  TFyz = TFyz+fyz
  TFzz = TFzz+fzz

  toAss_rhs = toAss_rhs + (fxx*vx(i,j,k)*vx(i,j,k)+fyy*vy(i,j,k)*vy(i,j,k)+fzz*vz(i,j,k)*vz(i,j,k) &
                     +TWO*(fxy*vx(i,j,k)*vy(i,j,k)+fxz*vx(i,j,k)*vz(i,j,k)+fyz*vy(i,j,k)*vz(i,j,k)))
  toAs1_rhs =   toAs1_rhs               +(fxx*vx(i,j,k)*ux(i,j,k)+fxy*vy(i,j,k)*ux(i,j,k)+fxz*vz(i,j,k)*ux(i,j,k) &
                                        + fxy*vx(i,j,k)*uy(i,j,k)+fyy*vy(i,j,k)*uy(i,j,k)+fyz*vz(i,j,k)*uy(i,j,k) &
                                        + fxz*vx(i,j,k)*uz(i,j,k)+fyz*vy(i,j,k)*uz(i,j,k)+fzz*vz(i,j,k)*uz(i,j,k))
  toAs2_rhs =   toAs2_rhs               +(fxx*vx(i,j,k)*wx(i,j,k)+fxy*vy(i,j,k)*wx(i,j,k)+fxz*vz(i,j,k)*wx(i,j,k) &
                                        + fxy*vx(i,j,k)*wy(i,j,k)+fyy*vy(i,j,k)*wy(i,j,k)+fyz*vz(i,j,k)*wy(i,j,k) &
                                        + fxz*vx(i,j,k)*wz(i,j,k)+fyz*vy(i,j,k)*wz(i,j,k)+fzz*vz(i,j,k)*wz(i,j,k))
  toAs1_rhs = toAs1_rhs-alpn1(i,j,k)*chin1(i,j,k)*(                              &
                         (gupxx(i,j,k)*(Axxx-(Gamxxx*Axx(i,j,k)+Gamyxx*Axy(i,j,k)+Gamzxx*Axz(i,j,k)))  &
                       +  gupxy(i,j,k)*(Axxy-(Gamxxx*Axy(i,j,k)+Gamyxx*Ayy(i,j,k)+Gamzxx*Ayz(i,j,k)))  &
                       +  gupxz(i,j,k)*(Axxz-(Gamxxx*Axz(i,j,k)+Gamyxx*Ayz(i,j,k)+Gamzxx*Azz(i,j,k)))  & 
                       +  gupxy(i,j,k)*(Axxy-(Gamxxy*Axx(i,j,k)+Gamyxy*Axy(i,j,k)+Gamzxy*Axz(i,j,k)))  &
                       +  gupyy(i,j,k)*(Axyy-(Gamxxy*Axy(i,j,k)+Gamyxy*Ayy(i,j,k)+Gamzxy*Ayz(i,j,k)))  &
                       +  gupyz(i,j,k)*(Axyz-(Gamxxy*Axz(i,j,k)+Gamyxy*Ayz(i,j,k)+Gamzxy*Azz(i,j,k)))  &
                       +  gupxz(i,j,k)*(Axxz-(Gamxxz*Axx(i,j,k)+Gamyxz*Axy(i,j,k)+Gamzxz*Axz(i,j,k)))  &
                       +  gupyz(i,j,k)*(Axyz-(Gamxxz*Axy(i,j,k)+Gamyxz*Ayy(i,j,k)+Gamzxz*Ayz(i,j,k)))  &
                       +  gupzz(i,j,k)*(Axzz-(Gamxxz*Axz(i,j,k)+Gamyxz*Ayz(i,j,k)+Gamzxz*Azz(i,j,k)))  &
                       -  (Gamxa*Axx(i,j,k)+Gamya*Axy(i,j,k)+Gamza*Axz(i,j,k)) )*ux(i,j,k)             &
                       + (gupxx(i,j,k)*(Axyx-(Gamxxy*Axx(i,j,k)+Gamyxy*Axy(i,j,k)+Gamzxy*Axz(i,j,k)))  &
                       +  gupxy(i,j,k)*(Axyy-(Gamxxy*Axy(i,j,k)+Gamyxy*Ayy(i,j,k)+Gamzxy*Ayz(i,j,k)))  &
                       +  gupxz(i,j,k)*(Axyz-(Gamxxy*Axz(i,j,k)+Gamyxy*Ayz(i,j,k)+Gamzxy*Azz(i,j,k)))  & 
                       +  gupxy(i,j,k)*(Axyy-(Gamxyy*Axx(i,j,k)+Gamyyy*Axy(i,j,k)+Gamzyy*Axz(i,j,k)))  &
                       +  gupyy(i,j,k)*(Ayyy-(Gamxyy*Axy(i,j,k)+Gamyyy*Ayy(i,j,k)+Gamzyy*Ayz(i,j,k)))  &
                       +  gupyz(i,j,k)*(Ayyz-(Gamxyy*Axz(i,j,k)+Gamyyy*Ayz(i,j,k)+Gamzyy*Azz(i,j,k)))  &
                       +  gupxz(i,j,k)*(Axyz-(Gamxyz*Axx(i,j,k)+Gamyyz*Axy(i,j,k)+Gamzyz*Axz(i,j,k)))  &
                       +  gupyz(i,j,k)*(Ayyz-(Gamxyz*Axy(i,j,k)+Gamyyz*Ayy(i,j,k)+Gamzyz*Ayz(i,j,k)))  &
                       +  gupzz(i,j,k)*(Ayzz-(Gamxyz*Axz(i,j,k)+Gamyyz*Ayz(i,j,k)+Gamzyz*Azz(i,j,k)))  &
                       -  (Gamxa*Axy(i,j,k)+Gamya*Ayy(i,j,k)+Gamza*Ayz(i,j,k)) )*uy(i,j,k)             &
                       + (gupxx(i,j,k)*(Axzx-(Gamxxz*Axx(i,j,k)+Gamyxz*Axy(i,j,k)+Gamzxz*Axz(i,j,k)))  &
                       +  gupxy(i,j,k)*(Axzy-(Gamxxz*Axy(i,j,k)+Gamyxz*Ayy(i,j,k)+Gamzxz*Ayz(i,j,k)))  &
                       +  gupxz(i,j,k)*(Axzz-(Gamxxz*Axz(i,j,k)+Gamyxz*Ayz(i,j,k)+Gamzxz*Azz(i,j,k)))  & 
                       +  gupxy(i,j,k)*(Axzy-(Gamxyz*Axx(i,j,k)+Gamyyz*Axy(i,j,k)+Gamzyz*Axz(i,j,k)))  &
                       +  gupyy(i,j,k)*(Ayzy-(Gamxyz*Axy(i,j,k)+Gamyyz*Ayy(i,j,k)+Gamzyz*Ayz(i,j,k)))  &
                       +  gupyz(i,j,k)*(Ayzz-(Gamxyz*Axz(i,j,k)+Gamyyz*Ayz(i,j,k)+Gamzyz*Azz(i,j,k)))  &
                       +  gupxz(i,j,k)*(Axzz-(Gamxzz*Axx(i,j,k)+Gamyzz*Axy(i,j,k)+Gamzzz*Axz(i,j,k)))  &
                       +  gupyz(i,j,k)*(Ayzz-(Gamxzz*Axy(i,j,k)+Gamyzz*Ayy(i,j,k)+Gamzzz*Ayz(i,j,k)))  &
                       +  gupzz(i,j,k)*(Azzz-(Gamxzz*Axz(i,j,k)+Gamyzz*Ayz(i,j,k)+Gamzzz*Azz(i,j,k)))  &
                       -  (Gamxa*Axz(i,j,k)+Gamya*Ayz(i,j,k)+Gamza*Azz(i,j,k)) )*uz(i,j,k)             &                 
                       -2.d0/3.d0*(Kx*ux(i,j,k)+Ky*uy(i,j,k)+Kz*uz(i,j,k))        &
                       -ONE/3.d0* (TZx*ux(i,j,k)+TZy*uy(i,j,k)+TZz*uz(i,j,k))     &
                       -1.5d0/chin1(i,j,k)*                                      &
                                (ux(i,j,k)*(gupxx(i,j,k)*chix*Axx(i,j,k)+gupxy(i,j,k)*chix*Axy(i,j,k)+gupxz(i,j,k)*chix*Axz(i,j,k) &
                                           +gupxy(i,j,k)*chiy*Axx(i,j,k)+gupyy(i,j,k)*chiy*Axy(i,j,k)+gupyz(i,j,k)*chiy*Axz(i,j,k) &
                                           +gupxz(i,j,k)*chiz*Axx(i,j,k)+gupyz(i,j,k)*chiz*Axy(i,j,k)+gupzz(i,j,k)*chiz*Axz(i,j,k)) &
                                +uy(i,j,k)*(gupxx(i,j,k)*chix*Axy(i,j,k)+gupxy(i,j,k)*chix*Ayy(i,j,k)+gupxz(i,j,k)*chix*Ayz(i,j,k) &
                                           +gupxy(i,j,k)*chiy*Axy(i,j,k)+gupyy(i,j,k)*chiy*Ayy(i,j,k)+gupyz(i,j,k)*chiy*Ayz(i,j,k) &
                                           +gupxz(i,j,k)*chiz*Axy(i,j,k)+gupyz(i,j,k)*chiz*Ayy(i,j,k)+gupzz(i,j,k)*chiz*Ayz(i,j,k)) &
                                +uz(i,j,k)*(gupxx(i,j,k)*chix*Axz(i,j,k)+gupxy(i,j,k)*chix*Ayz(i,j,k)+gupxz(i,j,k)*chix*Azz(i,j,k) &
                                           +gupxy(i,j,k)*chiy*Axz(i,j,k)+gupyy(i,j,k)*chiy*Ayz(i,j,k)+gupyz(i,j,k)*chiy*Azz(i,j,k) &
                                           +gupxz(i,j,k)*chiz*Axz(i,j,k)+gupyz(i,j,k)*chiz*Ayz(i,j,k)+gupzz(i,j,k)*chiz*Azz(i,j,k)) ) &
                       -0.5d0*kappa1*(Gmxcon(i,j,k)*ulx(i,j,k)+Gmycon(i,j,k)*uly(i,j,k)+Gmzcon(i,j,k)*ulz(i,j,k))  &
                       -(Rxx(i,j,k)*vx(i,j,k)*ux(i,j,k)+Rxy(i,j,k)*vy(i,j,k)*ux(i,j,k)+Rxz(i,j,k)*vz(i,j,k)*ux(i,j,k) &
                        +Rxy(i,j,k)*vx(i,j,k)*uy(i,j,k)+Ryy(i,j,k)*vy(i,j,k)*uy(i,j,k)+Ryz(i,j,k)*vz(i,j,k)*uy(i,j,k) &
                        +Rxz(i,j,k)*vx(i,j,k)*uz(i,j,k)+Ryz(i,j,k)*vy(i,j,k)*uz(i,j,k)+Rzz(i,j,k)*vz(i,j,k)*uz(i,j,k)) &
                       +0.5d0*chin1(i,j,k)*(ulx(i,j,k)*vx(i,j,k)*CAZxx+ulx(i,j,k)*vy(i,j,k)*CAZxy+ulx(i,j,k)*vz(i,j,k)*CAZxz    &
                                           +uly(i,j,k)*vx(i,j,k)*CAZyx+uly(i,j,k)*vy(i,j,k)*CAZyy+uly(i,j,k)*vz(i,j,k)*CAZyz    &
                                           +ulz(i,j,k)*vx(i,j,k)*CAZzx+ulz(i,j,k)*vy(i,j,k)*CAZzy+ulz(i,j,k)*vz(i,j,k)*CAZzz))
  toAs2_rhs = toAs2_rhs-alpn1(i,j,k)*chin1(i,j,k)*(                              &
                         (gupxx(i,j,k)*(Axxx-(Gamxxx*Axx(i,j,k)+Gamyxx*Axy(i,j,k)+Gamzxx*Axz(i,j,k)))  &
                       +  gupxy(i,j,k)*(Axxy-(Gamxxx*Axy(i,j,k)+Gamyxx*Ayy(i,j,k)+Gamzxx*Ayz(i,j,k)))  &
                       +  gupxz(i,j,k)*(Axxz-(Gamxxx*Axz(i,j,k)+Gamyxx*Ayz(i,j,k)+Gamzxx*Azz(i,j,k)))  & 
                       +  gupxy(i,j,k)*(Axxy-(Gamxxy*Axx(i,j,k)+Gamyxy*Axy(i,j,k)+Gamzxy*Axz(i,j,k)))  &
                       +  gupyy(i,j,k)*(Axyy-(Gamxxy*Axy(i,j,k)+Gamyxy*Ayy(i,j,k)+Gamzxy*Ayz(i,j,k)))  &
                       +  gupyz(i,j,k)*(Axyz-(Gamxxy*Axz(i,j,k)+Gamyxy*Ayz(i,j,k)+Gamzxy*Azz(i,j,k)))  &
                       +  gupxz(i,j,k)*(Axxz-(Gamxxz*Axx(i,j,k)+Gamyxz*Axy(i,j,k)+Gamzxz*Axz(i,j,k)))  &
                       +  gupyz(i,j,k)*(Axyz-(Gamxxz*Axy(i,j,k)+Gamyxz*Ayy(i,j,k)+Gamzxz*Ayz(i,j,k)))  &
                       +  gupzz(i,j,k)*(Axzz-(Gamxxz*Axz(i,j,k)+Gamyxz*Ayz(i,j,k)+Gamzxz*Azz(i,j,k)))  &
                       -  (Gamxa*Axx(i,j,k)+Gamya*Axy(i,j,k)+Gamza*Axz(i,j,k)) )*wx(i,j,k)             &
                       + (gupxx(i,j,k)*(Axyx-(Gamxxy*Axx(i,j,k)+Gamyxy*Axy(i,j,k)+Gamzxy*Axz(i,j,k)))  &
                       +  gupxy(i,j,k)*(Axyy-(Gamxxy*Axy(i,j,k)+Gamyxy*Ayy(i,j,k)+Gamzxy*Ayz(i,j,k)))  &
                       +  gupxz(i,j,k)*(Axyz-(Gamxxy*Axz(i,j,k)+Gamyxy*Ayz(i,j,k)+Gamzxy*Azz(i,j,k)))  & 
                       +  gupxy(i,j,k)*(Axyy-(Gamxyy*Axx(i,j,k)+Gamyyy*Axy(i,j,k)+Gamzyy*Axz(i,j,k)))  &
                       +  gupyy(i,j,k)*(Ayyy-(Gamxyy*Axy(i,j,k)+Gamyyy*Ayy(i,j,k)+Gamzyy*Ayz(i,j,k)))  &
                       +  gupyz(i,j,k)*(Ayyz-(Gamxyy*Axz(i,j,k)+Gamyyy*Ayz(i,j,k)+Gamzyy*Azz(i,j,k)))  &
                       +  gupxz(i,j,k)*(Axyz-(Gamxyz*Axx(i,j,k)+Gamyyz*Axy(i,j,k)+Gamzyz*Axz(i,j,k)))  &
                       +  gupyz(i,j,k)*(Ayyz-(Gamxyz*Axy(i,j,k)+Gamyyz*Ayy(i,j,k)+Gamzyz*Ayz(i,j,k)))  &
                       +  gupzz(i,j,k)*(Ayzz-(Gamxyz*Axz(i,j,k)+Gamyyz*Ayz(i,j,k)+Gamzyz*Azz(i,j,k)))  &
                       -  (Gamxa*Axy(i,j,k)+Gamya*Ayy(i,j,k)+Gamza*Ayz(i,j,k)) )*wy(i,j,k)             &
                       + (gupxx(i,j,k)*(Axzx-(Gamxxz*Axx(i,j,k)+Gamyxz*Axy(i,j,k)+Gamzxz*Axz(i,j,k)))  &
                       +  gupxy(i,j,k)*(Axzy-(Gamxxz*Axy(i,j,k)+Gamyxz*Ayy(i,j,k)+Gamzxz*Ayz(i,j,k)))  &
                       +  gupxz(i,j,k)*(Axzz-(Gamxxz*Axz(i,j,k)+Gamyxz*Ayz(i,j,k)+Gamzxz*Azz(i,j,k)))  & 
                       +  gupxy(i,j,k)*(Axzy-(Gamxyz*Axx(i,j,k)+Gamyyz*Axy(i,j,k)+Gamzyz*Axz(i,j,k)))  &
                       +  gupyy(i,j,k)*(Ayzy-(Gamxyz*Axy(i,j,k)+Gamyyz*Ayy(i,j,k)+Gamzyz*Ayz(i,j,k)))  &
                       +  gupyz(i,j,k)*(Ayzz-(Gamxyz*Axz(i,j,k)+Gamyyz*Ayz(i,j,k)+Gamzyz*Azz(i,j,k)))  &
                       +  gupxz(i,j,k)*(Axzz-(Gamxzz*Axx(i,j,k)+Gamyzz*Axy(i,j,k)+Gamzzz*Axz(i,j,k)))  &
                       +  gupyz(i,j,k)*(Ayzz-(Gamxzz*Axy(i,j,k)+Gamyzz*Ayy(i,j,k)+Gamzzz*Ayz(i,j,k)))  &
                       +  gupzz(i,j,k)*(Azzz-(Gamxzz*Axz(i,j,k)+Gamyzz*Ayz(i,j,k)+Gamzzz*Azz(i,j,k)))  &
                       -  (Gamxa*Axz(i,j,k)+Gamya*Ayz(i,j,k)+Gamza*Azz(i,j,k)) )*wz(i,j,k)             &                 
                       -2.d0/3.d0*(Kx*wx(i,j,k)+ky*wy(i,j,k)+Kz*wz(i,j,k))        &
                       -ONE/3.d0* (TZx*wx(i,j,k)+TZy*wy(i,j,k)+TZz*wz(i,j,k))     &
                       -1.5d0/chin1(i,j,k)*                                      &
                                (wx(i,j,k)*(gupxx(i,j,k)*chix*Axx(i,j,k)+gupxy(i,j,k)*chix*Axy(i,j,k)+gupxz(i,j,k)*chix*Axz(i,j,k) &
                                           +gupxy(i,j,k)*chiy*Axx(i,j,k)+gupyy(i,j,k)*chiy*Axy(i,j,k)+gupyz(i,j,k)*chiy*Axz(i,j,k) &
                                           +gupxz(i,j,k)*chiz*Axx(i,j,k)+gupyz(i,j,k)*chiz*Axy(i,j,k)+gupzz(i,j,k)*chiz*Axz(i,j,k)) &
                                +wy(i,j,k)*(gupxx(i,j,k)*chix*Axy(i,j,k)+gupxy(i,j,k)*chix*Ayy(i,j,k)+gupxz(i,j,k)*chix*Ayz(i,j,k) &
                                           +gupxy(i,j,k)*chiy*Axy(i,j,k)+gupyy(i,j,k)*chiy*Ayy(i,j,k)+gupyz(i,j,k)*chiy*Ayz(i,j,k) &
                                           +gupxz(i,j,k)*chiz*Axy(i,j,k)+gupyz(i,j,k)*chiz*Ayy(i,j,k)+gupzz(i,j,k)*chiz*Ayz(i,j,k)) &
                                +wz(i,j,k)*(gupxx(i,j,k)*chix*Axz(i,j,k)+gupxy(i,j,k)*chix*Ayz(i,j,k)+gupxz(i,j,k)*chix*Azz(i,j,k) &
                                           +gupxy(i,j,k)*chiy*Axz(i,j,k)+gupyy(i,j,k)*chiy*Ayz(i,j,k)+gupyz(i,j,k)*chiy*Azz(i,j,k) &
                                           +gupxz(i,j,k)*chiz*Axz(i,j,k)+gupyz(i,j,k)*chiz*Ayz(i,j,k)+gupzz(i,j,k)*chiz*Azz(i,j,k)) ) &
                       -0.5d0*kappa1*(Gmxcon(i,j,k)*wlx(i,j,k)+Gmycon(i,j,k)*wly(i,j,k)+Gmzcon(i,j,k)*wlz(i,j,k))  &
                       -(Rxx(i,j,k)*vx(i,j,k)*wx(i,j,k)+Rxy(i,j,k)*vy(i,j,k)*wx(i,j,k)+Rxz(i,j,k)*vz(i,j,k)*wx(i,j,k) &
                        +Rxy(i,j,k)*vx(i,j,k)*wy(i,j,k)+Ryy(i,j,k)*vy(i,j,k)*wy(i,j,k)+Ryz(i,j,k)*vz(i,j,k)*wy(i,j,k) &
                        +Rxz(i,j,k)*vx(i,j,k)*wz(i,j,k)+Ryz(i,j,k)*vy(i,j,k)*wz(i,j,k)+Rzz(i,j,k)*vz(i,j,k)*wz(i,j,k)) &
                       +0.5d0*chin1(i,j,k)*(wlx(i,j,k)*vx(i,j,k)*CAZxx+wlx(i,j,k)*vy(i,j,k)*CAZxy+wlx(i,j,k)*vz(i,j,k)*CAZxz    &
                                           +wly(i,j,k)*vx(i,j,k)*CAZyx+wly(i,j,k)*vy(i,j,k)*CAZyy+wly(i,j,k)*vz(i,j,k)*CAZyz    &
                                           +wlz(i,j,k)*vx(i,j,k)*CAZzx+wlz(i,j,k)*vy(i,j,k)*CAZzy+wlz(i,j,k)*vz(i,j,k)*CAZzz))

  toGam1_rhs = -alpn1(i,j,k)*dsqrt(tmuST)*((Gamxx*vx(i,j,k)*ulx(i,j,k)+Gamxy*vy(i,j,k)*ulx(i,j,k)+Gamxz*vz(i,j,k)*ulx(i,j,k) &
                                           +Gamyx*vx(i,j,k)*uly(i,j,k)+Gamyy*vy(i,j,k)*uly(i,j,k)+Gamyz*vz(i,j,k)*uly(i,j,k) &
                                           +Gamzx*vx(i,j,k)*ulz(i,j,k)+Gamzy*vy(i,j,k)*ulz(i,j,k)+Gamzz*vz(i,j,k)*ulz(i,j,k)) &
                                          -(Gamxx*ux(i,j,k)*slx(i,j,k)+Gamxy*uy(i,j,k)*slx(i,j,k)+Gamxz*uz(i,j,k)*slx(i,j,k) &
                                           +Gamyx*ux(i,j,k)*sly(i,j,k)+Gamyy*uy(i,j,k)*sly(i,j,k)+Gamyz*uz(i,j,k)*sly(i,j,k) &
                                           +Gamzx*ux(i,j,k)*slz(i,j,k)+Gamzy*uy(i,j,k)*slz(i,j,k)+Gamzz*uz(i,j,k)*slz(i,j,k))/chin1(i,j,k) ) &
               +((qupxx(i,j,k)*sfxxx+qupxx(i,j,k)*sfxxx+qupxx(i,j,k)*sfxxx                  &
            +TWO*(qupxy(i,j,k)*sfxxy+qupxz(i,j,k)*sfxxz+qupyz(i,j,k)*sfxyz))*ulx(i,j,k)     &
                +(qupxx(i,j,k)*sfyxx+qupxx(i,j,k)*sfyxx+qupxx(i,j,k)*sfyxx                  &
            +TWO*(qupxy(i,j,k)*sfyxy+qupxz(i,j,k)*sfyxz+qupyz(i,j,k)*sfyyz))*uly(i,j,k)     &
                +(qupxx(i,j,k)*sfzxx+qupxx(i,j,k)*sfzxx+qupxx(i,j,k)*sfzxx                  &
            +TWO*(qupxy(i,j,k)*sfzxy+qupxz(i,j,k)*sfzxz+qupyz(i,j,k)*sfzyz))*ulz(i,j,k)     &
                )/chin1(i,j,k)                                                              &
            +4.d0/3.d0/chin1(i,j,k)*(ux(i,j,k)*(vx(i,j,k)*slx(i,j,k)*sfxxx+vx(i,j,k)*sly(i,j,k)*sfyxx+vx(i,j,k)*slz(i,j,k)*sfzxx  &
                                               +vy(i,j,k)*slx(i,j,k)*sfxxy+vy(i,j,k)*sly(i,j,k)*sfyxy+vy(i,j,k)*slz(i,j,k)*sfzxy  &
                                               +vz(i,j,k)*slx(i,j,k)*sfxxz+vz(i,j,k)*sly(i,j,k)*sfyxz+vz(i,j,k)*slz(i,j,k)*sfzxz) &
                                    +uy(i,j,k)*(vx(i,j,k)*slx(i,j,k)*sfxxy+vx(i,j,k)*sly(i,j,k)*sfyxy+vx(i,j,k)*slz(i,j,k)*sfzxy  &
                                               +vy(i,j,k)*slx(i,j,k)*sfxyy+vy(i,j,k)*sly(i,j,k)*sfyyy+vy(i,j,k)*slz(i,j,k)*sfzyy  &
                                               +vz(i,j,k)*slx(i,j,k)*sfxyz+vz(i,j,k)*sly(i,j,k)*sfyyz+vz(i,j,k)*slz(i,j,k)*sfzyz) &
                                    +uz(i,j,k)*(vx(i,j,k)*slx(i,j,k)*sfxxz+vx(i,j,k)*sly(i,j,k)*sfyxz+vx(i,j,k)*slz(i,j,k)*sfzxz  &
                                               +vy(i,j,k)*slx(i,j,k)*sfxyz+vy(i,j,k)*sly(i,j,k)*sfyyz+vy(i,j,k)*slz(i,j,k)*sfzyz  &
                                               +vz(i,j,k)*slx(i,j,k)*sfxzz+vz(i,j,k)*sly(i,j,k)*sfyzz+vz(i,j,k)*slz(i,j,k)*sfzzz)) &
            +ONE/3.d0/chin1(i,j,k)* (ux(i,j,k)*(qulxx(i,j,k)*sfxxx+qulxy(i,j,k)*sfyxx+qulxz(i,j,k)*sfzxx  &
                                               +qulyx(i,j,k)*sfxxy+qulyy(i,j,k)*sfyxy+qulyz(i,j,k)*sfzxy  &
                                               +qulzx(i,j,k)*sfxxz+qulzy(i,j,k)*sfyxz+qulzz(i,j,k)*sfzxz) &
                                    +uy(i,j,k)*(qulxx(i,j,k)*sfxxy+qulxy(i,j,k)*sfyxy+qulxz(i,j,k)*sfzxy  &
                                               +qulyx(i,j,k)*sfxyy+qulyy(i,j,k)*sfyyy+qulyz(i,j,k)*sfzyy  &
                                               +qulzx(i,j,k)*sfxyz+qulzy(i,j,k)*sfyyz+qulzz(i,j,k)*sfzyz) &
                                    +uz(i,j,k)*(qulxx(i,j,k)*sfxxz+qulxy(i,j,k)*sfyxz+qulxz(i,j,k)*sfzxz  &
                                               +qulyx(i,j,k)*sfxyz+qulyy(i,j,k)*sfyyz+qulyz(i,j,k)*sfzyz  &
                                               +qulzx(i,j,k)*sfxzz+qulzy(i,j,k)*sfyzz+qulzz(i,j,k)*sfzzz)) &
            -2.d0/3.d0*alpn1(i,j,k)/chin1(i,j,k)*(ux(i,j,k)*(TWO*Kx+TZx)+uy(i,j,k)*(TWO*Ky+TZy)+uz(i,j,k)*(TWO*Kz+TZz))  &
            +hu-kappa3*alpn1(i,j,k)*(Gamx(i,j,k)*ulx(i,j,k)+Gamx(i,j,k)*uly(i,j,k)+Gamz(i,j,k)*ulz(i,j,k))               &
            +(betax(i,j,k)*Gamxx+betay(i,j,k)*Gamxy+betaz(i,j,k)*Gamxz)*ulx(i,j,k)                                       &
            +(betax(i,j,k)*Gamyx+betay(i,j,k)*Gamyy+betaz(i,j,k)*Gamyz)*uly(i,j,k)                                       &
            +(betax(i,j,k)*Gamzx+betay(i,j,k)*Gamzy+betaz(i,j,k)*Gamzz)*ulz(i,j,k)

  toGam2_rhs = -alpn1(i,j,k)*dsqrt(tmuST)*((Gamxx*vx(i,j,k)*wlx(i,j,k)+Gamxy*vy(i,j,k)*wlx(i,j,k)+Gamxz*vz(i,j,k)*wlx(i,j,k) &
                                           +Gamyx*vx(i,j,k)*wly(i,j,k)+Gamyy*vy(i,j,k)*wly(i,j,k)+Gamyz*vz(i,j,k)*wly(i,j,k) &
                                           +Gamzx*vx(i,j,k)*wlz(i,j,k)+Gamzy*vy(i,j,k)*wlz(i,j,k)+Gamzz*vz(i,j,k)*wlz(i,j,k)) &
                                          -(Gamxx*wx(i,j,k)*slx(i,j,k)+Gamxy*wy(i,j,k)*slx(i,j,k)+Gamxz*wz(i,j,k)*slx(i,j,k) &
                                           +Gamyx*wx(i,j,k)*sly(i,j,k)+Gamyy*wy(i,j,k)*sly(i,j,k)+Gamyz*wz(i,j,k)*sly(i,j,k) &
                                           +Gamzx*wx(i,j,k)*slz(i,j,k)+Gamzy*wy(i,j,k)*slz(i,j,k)+Gamzz*wz(i,j,k)*slz(i,j,k))/chin1(i,j,k) ) &
               +((qupxx(i,j,k)*sfxxx+qupxx(i,j,k)*sfxxx+qupxx(i,j,k)*sfxxx                  &
            +TWO*(qupxy(i,j,k)*sfxxy+qupxz(i,j,k)*sfxxz+qupyz(i,j,k)*sfxyz))*wlx(i,j,k)     &
                +(qupxx(i,j,k)*sfyxx+qupxx(i,j,k)*sfyxx+qupxx(i,j,k)*sfyxx                  &
            +TWO*(qupxy(i,j,k)*sfyxy+qupxz(i,j,k)*sfyxz+qupyz(i,j,k)*sfyyz))*wly(i,j,k)     &
                +(qupxx(i,j,k)*sfzxx+qupxx(i,j,k)*sfzxx+qupxx(i,j,k)*sfzxx                  &
            +TWO*(qupxy(i,j,k)*sfzxy+qupxz(i,j,k)*sfzxz+qupyz(i,j,k)*sfzyz))*wlz(i,j,k)     &
                )/chin1(i,j,k)                                                              &
            +4.d0/3.d0/chin1(i,j,k)*(wx(i,j,k)*(vx(i,j,k)*slx(i,j,k)*sfxxx+vx(i,j,k)*sly(i,j,k)*sfyxx+vx(i,j,k)*slz(i,j,k)*sfzxx  &
                                                +vy(i,j,k)*slx(i,j,k)*sfxxy+vy(i,j,k)*sly(i,j,k)*sfyxy+vy(i,j,k)*slz(i,j,k)*sfzxy  &
                                                +vz(i,j,k)*slx(i,j,k)*sfxxz+vz(i,j,k)*sly(i,j,k)*sfyxz+vz(i,j,k)*slz(i,j,k)*sfzxz) &
                                     +wy(i,j,k)*(vx(i,j,k)*slx(i,j,k)*sfxxy+vx(i,j,k)*sly(i,j,k)*sfyxy+vx(i,j,k)*slz(i,j,k)*sfzxy  &
                                                +vy(i,j,k)*slx(i,j,k)*sfxyy+vy(i,j,k)*sly(i,j,k)*sfyyy+vy(i,j,k)*slz(i,j,k)*sfzyy  &
                                                +vz(i,j,k)*slx(i,j,k)*sfxyz+vz(i,j,k)*sly(i,j,k)*sfyyz+vz(i,j,k)*slz(i,j,k)*sfzyz) &
                                     +wz(i,j,k)*(vx(i,j,k)*slx(i,j,k)*sfxxz+vx(i,j,k)*sly(i,j,k)*sfyxz+vx(i,j,k)*slz(i,j,k)*sfzxz  &
                                                +vy(i,j,k)*slx(i,j,k)*sfxyz+vy(i,j,k)*sly(i,j,k)*sfyyz+vy(i,j,k)*slz(i,j,k)*sfzyz  &
                                                +vz(i,j,k)*slx(i,j,k)*sfxzz+vz(i,j,k)*sly(i,j,k)*sfyzz+vz(i,j,k)*slz(i,j,k)*sfzzz)) &
            +ONE/3.d0/chin1(i,j,k)* (wx(i,j,k)*(qulxx(i,j,k)*sfxxx+qulxy(i,j,k)*sfyxx+qulxz(i,j,k)*sfzxx  &
                                                +qulyx(i,j,k)*sfxxy+qulyy(i,j,k)*sfyxy+qulyz(i,j,k)*sfzxy  &
                                                +qulzx(i,j,k)*sfxxz+qulzy(i,j,k)*sfyxz+qulzz(i,j,k)*sfzxz) &
                                     +wy(i,j,k)*(qulxx(i,j,k)*sfxxy+qulxy(i,j,k)*sfyxy+qulxz(i,j,k)*sfzxy  &
                                                +qulyx(i,j,k)*sfxyy+qulyy(i,j,k)*sfyyy+qulyz(i,j,k)*sfzyy  &
                                                +qulzx(i,j,k)*sfxyz+qulzy(i,j,k)*sfyyz+qulzz(i,j,k)*sfzyz) &
                                     +wz(i,j,k)*(qulxx(i,j,k)*sfxxz+qulxy(i,j,k)*sfyxz+qulxz(i,j,k)*sfzxz  &
                                                +qulyx(i,j,k)*sfxyz+qulyy(i,j,k)*sfyyz+qulyz(i,j,k)*sfzyz  &
                                                +qulzx(i,j,k)*sfxzz+qulzy(i,j,k)*sfyzz+qulzz(i,j,k)*sfzzz)) &
            -2.d0/3.d0*alpn1(i,j,k)/chin1(i,j,k)*(wx(i,j,k)*(TWO*Kx+TZx)+wy(i,j,k)*(TWO*Ky+TZy)+wz(i,j,k)*(TWO*Kz+TZz))  &
            +hw-kappa3*alpn1(i,j,k)*(Gamx(i,j,k)*wlx(i,j,k)+Gamx(i,j,k)*wly(i,j,k)+Gamz(i,j,k)*wlz(i,j,k))               &
            +(betax(i,j,k)*Gamxx+betay(i,j,k)*Gamxy+betaz(i,j,k)*Gamxz)*wlx(i,j,k)                                       &
            +(betax(i,j,k)*Gamyx+betay(i,j,k)*Gamyy+betaz(i,j,k)*Gamyz)*wly(i,j,k)                                       &
            +(betax(i,j,k)*Gamzx+betay(i,j,k)*Gamzy+betaz(i,j,k)*Gamzz)*wlz(i,j,k)

! \tilde{D} A_ij            
  gxxx = Axxx-TWO*(Gamxxx*Axx(i,j,k)+Gamyxx*Axy(i,j,k)+Gamzxx*Axz(i,j,k))
  gxxy = Axxy-TWO*(Gamxxy*Axx(i,j,k)+Gamyxy*Axy(i,j,k)+Gamzxy*Axz(i,j,k))
  gxxz = Axxz-TWO*(Gamxxz*Axx(i,j,k)+Gamyxz*Axy(i,j,k)+Gamzxz*Axz(i,j,k))
  gyyx = Ayyx-TWO*(Gamxxy*Axy(i,j,k)+Gamyxy*Ayy(i,j,k)+Gamzxy*Ayz(i,j,k))
  gyyy = Ayyy-TWO*(Gamxyy*Axy(i,j,k)+Gamyyy*Ayy(i,j,k)+Gamzyy*Ayz(i,j,k))
  gyyz = Ayyz-TWO*(Gamxyz*Axy(i,j,k)+Gamyyz*Ayy(i,j,k)+Gamzyz*Ayz(i,j,k))
  gzzx = Azzx-TWO*(Gamxxz*Axz(i,j,k)+Gamyxz*Ayz(i,j,k)+Gamzxz*Azz(i,j,k))
  gzzy = Azzy-TWO*(Gamxyz*Axz(i,j,k)+Gamyyz*Ayz(i,j,k)+Gamzyz*Azz(i,j,k))
  gzzz = Azzz-TWO*(Gamxzz*Axz(i,j,k)+Gamyzz*Ayz(i,j,k)+Gamzzz*Azz(i,j,k))
  gxyx = Axyx-(Gamxxy*Axx(i,j,k)+Gamyxy*Axy(i,j,k)+Gamzxy*Axz(i,j,k)+Gamxxx*Axy(i,j,k)+Gamyxx*Ayy(i,j,k)+Gamzxx*Ayz(i,j,k))
  gxyy = Axyy-(Gamxyy*Axx(i,j,k)+Gamyyy*Axy(i,j,k)+Gamzyy*Axz(i,j,k)+Gamxxy*Axy(i,j,k)+Gamyxy*Ayy(i,j,k)+Gamzxy*Ayz(i,j,k))
  gxyz = Axyz-(Gamxyz*Axx(i,j,k)+Gamyyz*Axy(i,j,k)+Gamzyz*Axz(i,j,k)+Gamxxz*Axy(i,j,k)+Gamyxz*Ayy(i,j,k)+Gamzxz*Ayz(i,j,k))
  gxzx = Axzx-(Gamxxz*Axx(i,j,k)+Gamyxz*Axy(i,j,k)+Gamzxz*Axz(i,j,k)+Gamxxx*Axz(i,j,k)+Gamyxx*Ayz(i,j,k)+Gamzxx*Azz(i,j,k))
  gxzy = Axzy-(Gamxyz*Axx(i,j,k)+Gamyyz*Axy(i,j,k)+Gamzyz*Axz(i,j,k)+Gamxxy*Axz(i,j,k)+Gamyxy*Ayz(i,j,k)+Gamzxy*Azz(i,j,k))
  gxzz = Axzz-(Gamxzz*Axx(i,j,k)+Gamyzz*Axy(i,j,k)+Gamzzz*Axz(i,j,k)+Gamxxz*Axz(i,j,k)+Gamyxz*Ayz(i,j,k)+Gamzxz*Azz(i,j,k))
  gyzx = Ayzx-(Gamxxz*Axy(i,j,k)+Gamyxz*Ayy(i,j,k)+Gamzxz*Ayz(i,j,k)+Gamxxy*Axz(i,j,k)+Gamyxy*Ayz(i,j,k)+Gamzxy*Azz(i,j,k))
  gyzy = Ayzy-(Gamxyz*Axy(i,j,k)+Gamyyz*Ayy(i,j,k)+Gamzyz*Ayz(i,j,k)+Gamxyy*Axz(i,j,k)+Gamyyy*Ayz(i,j,k)+Gamzyy*Azz(i,j,k))
  gyzz = Ayzz-(Gamxzz*Axy(i,j,k)+Gamyzz*Ayy(i,j,k)+Gamzzz*Ayz(i,j,k)+Gamxyz*Axz(i,j,k)+Gamyyz*Ayz(i,j,k)+Gamzyz*Azz(i,j,k))

  f = (trK(i,j,k)+TWO*TZ(i,j,k))*TWO/3.d0
  fxx =  (vx(i,j,k)*gxxx + vy(i,j,k)*gxxy + vz(i,j,k)*gxxz &
        -(vx(i,j,k)*gxxx + vy(i,j,k)*gxyx + vz(i,j,k)*gxzx)) &
        +AAxx-f*Axx(i,j,k)
  fyy =  (vx(i,j,k)*gyyx + vy(i,j,k)*gyyy + vz(i,j,k)*gyyz &
        -(vx(i,j,k)*gxyy + vy(i,j,k)*gyyy + vz(i,j,k)*gyzy)) &
        +AAyy-f*Ayy(i,j,k)
  fzz =  (vx(i,j,k)*gzzx + vy(i,j,k)*gzzy + vz(i,j,k)*gzzz &
        -(vx(i,j,k)*gxzz + vy(i,j,k)*gyzz + vz(i,j,k)*gzzz)) &
        +AAzz-f*Azz(i,j,k)
  fxy =  (vx(i,j,k)*gxyx + vy(i,j,k)*gxyy + vz(i,j,k)*gxyz & 
        -(vx(i,j,k)*gxxy + vy(i,j,k)*gxyy + vz(i,j,k)*gxzy + vx(i,j,k)*gxyx + vy(i,j,k)*gyyx + vz(i,j,k)*gyzx)/TWO) &
        +AAxy-f*Axy(i,j,k)
  fxz =  (vx(i,j,k)*gxzx + vy(i,j,k)*gxzy + vz(i,j,k)*gxzz & 
        -(vx(i,j,k)*gxxz + vy(i,j,k)*gxyz + vz(i,j,k)*gxzz + vx(i,j,k)*gyzx + vy(i,j,k)*gyzx + vz(i,j,k)*gzzx)/TWO) &
        +AAxz-f*Axz(i,j,k)
  fyz =  (vx(i,j,k)*gyzx + vy(i,j,k)*gyzy + vz(i,j,k)*gyzz & 
        -(vx(i,j,k)*gxyz + vy(i,j,k)*gyyz + vz(i,j,k)*gyzz + vx(i,j,k)*gyzy + vy(i,j,k)*gyzy + vz(i,j,k)*gzzy)/TWO) &
        +AAyz-f*Ayz(i,j,k)

! 1/2 A_ij D_k(ln chi)        
  gxxx = Axx(i,j,k)*chix/TWO/chin1(i,j,k)
  gxxy = Axx(i,j,k)*chiy/TWO/chin1(i,j,k)
  gxxz = Axx(i,j,k)*chiz/TWO/chin1(i,j,k)
  gxyx = Axy(i,j,k)*chix/TWO/chin1(i,j,k)
  gxyy = Axy(i,j,k)*chiy/TWO/chin1(i,j,k)
  gxyz = Axy(i,j,k)*chiz/TWO/chin1(i,j,k)
  gxzx = Axz(i,j,k)*chix/TWO/chin1(i,j,k)
  gxzy = Axz(i,j,k)*chiy/TWO/chin1(i,j,k)
  gxzz = Axz(i,j,k)*chiz/TWO/chin1(i,j,k)
  gyyx = Ayy(i,j,k)*chix/TWO/chin1(i,j,k)
  gyyy = Ayy(i,j,k)*chiy/TWO/chin1(i,j,k)
  gyyz = Ayy(i,j,k)*chiz/TWO/chin1(i,j,k)
  gyzx = Ayz(i,j,k)*chix/TWO/chin1(i,j,k)
  gyzy = Ayz(i,j,k)*chiy/TWO/chin1(i,j,k)
  gyzz = Ayz(i,j,k)*chiz/TWO/chin1(i,j,k)
  gzzx = Azz(i,j,k)*chix/TWO/chin1(i,j,k)
  gzzy = Azz(i,j,k)*chiy/TWO/chin1(i,j,k)
  gzzz = Azz(i,j,k)*chiz/TWO/chin1(i,j,k)

  fxx = fxx - (vx(i,j,k)*gxxx + vy(i,j,k)*gxxy + vz(i,j,k)*gxxz &
             -(vx(i,j,k)*gxxx + vy(i,j,k)*gxyx + vz(i,j,k)*gxzx))
  fyy = fyy - (vx(i,j,k)*gyyx + vy(i,j,k)*gyyy + vz(i,j,k)*gyyz &
             -(vx(i,j,k)*gxyy + vy(i,j,k)*gyyy + vz(i,j,k)*gyzy))
  fzz = fzz - (vx(i,j,k)*gzzx + vy(i,j,k)*gzzy + vz(i,j,k)*gzzz &
             -(vx(i,j,k)*gxzz + vy(i,j,k)*gyzz + vz(i,j,k)*gzzz))
  fxy = fxy - (vx(i,j,k)*gxyx + vy(i,j,k)*gxyy + vz(i,j,k)*gxyz & 
             -(vx(i,j,k)*gxxy + vy(i,j,k)*gxyy + vz(i,j,k)*gxzy + vx(i,j,k)*gxyx + vy(i,j,k)*gyyx + vz(i,j,k)*gyzx)/TWO)
  fxz = fxz - (vx(i,j,k)*gxzx + vy(i,j,k)*gxzy + vz(i,j,k)*gxzz & 
             -(vx(i,j,k)*gxxz + vy(i,j,k)*gxyz + vz(i,j,k)*gxzz + vx(i,j,k)*gyzx + vy(i,j,k)*gyzx + vz(i,j,k)*gzzx)/TWO)
  fyz = fyz - (vx(i,j,k)*gyzx + vy(i,j,k)*gyzy + vz(i,j,k)*gyzz & 
             -(vx(i,j,k)*gxyz + vy(i,j,k)*gyyz + vz(i,j,k)*gyzz + vx(i,j,k)*gyzy + vy(i,j,k)*gyzy + vz(i,j,k)*gzzy)/TWO)

  TFxx = TFxx-alpn1(i,j,k)*fxx
  TFxy = TFxy-alpn1(i,j,k)*fxy
  TFxz = TFxz-alpn1(i,j,k)*fxz
  TFyy = TFyy-alpn1(i,j,k)*fyy
  TFyz = TFyz-alpn1(i,j,k)*fyz
  TFzz = TFzz-alpn1(i,j,k)*fzz

  f = 0.5d0*(qupxx(i,j,k)*TFxx+qupyy(i,j,k)*TFyy+qupzz(i,j,k)*TFzz &
       +TWO*(qupxy(i,j,k)*TFxy+qupxz(i,j,k)*TFxz+qupyz(i,j,k)*TFyz))

  toA11_rhs =  ux(i,j,k)*ux(i,j,k)*TFxx+uy(i,j,k)*uy(i,j,k)*TFyy+uz(i,j,k)*uz(i,j,k)*TFzz+ &
          TWO*(ux(i,j,k)*uy(i,j,k)*TFxy+ux(i,j,k)*uz(i,j,k)*TFxz+uy(i,j,k)*uz(i,j,k)*TFyz)-f
  toA22_rhs =  wx(i,j,k)*wx(i,j,k)*TFxx+wy(i,j,k)*wy(i,j,k)*TFyy+wz(i,j,k)*wz(i,j,k)*TFzz+ &
          TWO*(wx(i,j,k)*wy(i,j,k)*TFxy+wx(i,j,k)*wz(i,j,k)*TFxz+wy(i,j,k)*wz(i,j,k)*TFyz)-f
  toA12_rhs =  ux(i,j,k)*wx(i,j,k)*TFxx+ux(i,j,k)*wy(i,j,k)*TFxy+ux(i,j,k)*wz(i,j,k)*TFxz &
              +uy(i,j,k)*wx(i,j,k)*TFxy+uy(i,j,k)*wy(i,j,k)*TFyy+uy(i,j,k)*wz(i,j,k)*TFyz &
              +uz(i,j,k)*wx(i,j,k)*TFxz+uz(i,j,k)*wy(i,j,k)*TFyz+uz(i,j,k)*wz(i,j,k)*TFzz

  toA11_rhs = toA11_rhs +alpn1(i,j,k)*chin1(i,j,k)*Rhpsi0
  toA22_rhs = toA22_rhs -alpn1(i,j,k)*chin1(i,j,k)*Rhpsi0
  toA12_rhs = toA12_rhs +alpn1(i,j,k)*chin1(i,j,k)*Ihpsi0

#if 0  
  toAqq_rhs = qupxx(i,j,k)*Axx_rhs(i,j,k)+qupyy(i,j,k)*Ayy_rhs(i,j,k)+qupzz(i,j,k)*Azz_rhs(i,j,k) &
             +TWO*(qupxy(i,j,k)*Axy_rhs(i,j,k)+qupxz(i,j,k)*Axz_rhs(i,j,k)+qupyz(i,j,k)*Ayz_rhs(i,j,k))
#else
  Ainvxx = gupxx(i,j,k)*gupxx(i,j,k)*Axx(i,j,k)+2.0*gupxx(i,j,k)*gupxy(i,j,k)*Axy(i,j,k)+ &
       2.0*gupxx(i,j,k)*gupxz(i,j,k)*Axz(i,j,k)+gupxy(i,j,k)*gupxy(i,j,k)*Ayy(i,j,k)+     &
       2.0*gupxy(i,j,k)*gupxz(i,j,k)*Ayz(i,j,k)+gupxz(i,j,k)*gupxz(i,j,k)*Azz(i,j,k)

  Ainvxy = gupxx(i,j,k)*gupxy(i,j,k)*Axx(i,j,k)+gupxx(i,j,k)*gupyy(i,j,k)*Axy(i,j,k)+     &
           gupxx(i,j,k)*gupyz(i,j,k)*Axz(i,j,k)+gupxy(i,j,k)*gupxy(i,j,k)*Axy(i,j,k)+     &
           gupxy(i,j,k)*gupyy(i,j,k)*Ayy(i,j,k)+gupxy(i,j,k)*gupyz(i,j,k)*Ayz(i,j,k)+     &
           gupxz(i,j,k)*gupxy(i,j,k)*Axz(i,j,k)+gupxz(i,j,k)*gupyy(i,j,k)*Ayz(i,j,k)+gupxz(i,j,k)*gupyz(i,j,k)*Azz(i,j,k)

  Ainvxz = gupxx(i,j,k)*gupxz(i,j,k)*Axx(i,j,k)+gupxx(i,j,k)*gupyz(i,j,k)*Axy(i,j,k)+     &
           gupxx(i,j,k)*gupzz(i,j,k)*Axz(i,j,k)+gupxy(i,j,k)*gupxz(i,j,k)*Axy(i,j,k)+     &
           gupxy(i,j,k)*gupyz(i,j,k)*Ayy(i,j,k)+gupxy(i,j,k)*gupzz(i,j,k)*Ayz(i,j,k)+     &
           gupxz(i,j,k)*gupxz(i,j,k)*Axz(i,j,k)+gupxz(i,j,k)*gupyz(i,j,k)*Ayz(i,j,k)+gupxz(i,j,k)*gupzz(i,j,k)*Azz(i,j,k)
  Ainvyy = gupxy(i,j,k)*gupxy(i,j,k)*Axx(i,j,k)+2.0*gupxy(i,j,k)*gupyy(i,j,k)*Axy(i,j,k)+ &
       2.0*gupxy(i,j,k)*gupyz(i,j,k)*Axz(i,j,k)+gupyy(i,j,k)*gupyy(i,j,k)*Ayy(i,j,k)+     &
       2.0*gupyy(i,j,k)*gupyz(i,j,k)*Ayz(i,j,k)+gupyz(i,j,k)*gupyz(i,j,k)*Azz(i,j,k)

  Ainvyz = gupxy(i,j,k)*gupxz(i,j,k)*Axx(i,j,k)+gupxy(i,j,k)*gupyz(i,j,k)*Axy(i,j,k)+     &
           gupxy(i,j,k)*gupzz(i,j,k)*Axz(i,j,k)+gupyy(i,j,k)*gupxz(i,j,k)*Axy(i,j,k)+     &
           gupyy(i,j,k)*gupyz(i,j,k)*Ayy(i,j,k)+gupyy(i,j,k)*gupzz(i,j,k)*Ayz(i,j,k)+     &
           gupyz(i,j,k)*gupxz(i,j,k)*Axz(i,j,k)+gupyz(i,j,k)*gupyz(i,j,k)*Ayz(i,j,k)+gupyz(i,j,k)*gupzz(i,j,k)*Azz(i,j,k)
  Ainvzz = gupxz(i,j,k)*gupxz(i,j,k)*Axx(i,j,k)+2.0*gupxz(i,j,k)*gupyz(i,j,k)*Axy(i,j,k)+ &
       2.0*gupxz(i,j,k)*gupzz(i,j,k)*Axz(i,j,k)+gupyz(i,j,k)*gupyz(i,j,k)*Ayy(i,j,k)+     &
       2.0*gupyz(i,j,k)*gupzz(i,j,k)*Ayz(i,j,k)+gupzz(i,j,k)*gupzz(i,j,k)*Azz(i,j,k)

  toAqq_rhs = -TWO*alpn1(i,j,k)*chin1(i,j,k)*(gupxx(i,j,k)*AAxx+gupyy(i,j,k)*AAyy+gupzz(i,j,k)*AAzz &
             +TWO*(gupxy(i,j,k)*AAxy+gupxz(i,j,k)*AAxz+gupyz(i,j,k)*AAyz))+chin1(i,j,k)*(Ainvxx+liegxx+Ainvyy*liegyy+Ainvzz*liegzz &
             +TWO*(Ainvxy*liegxy+Ainvxz*liegxz+Ainvyz*liegyz))-toAss_rhs
#endif
! reconstruct rhs for dynamical variables
  trK_rhs(i,j,k) = totrK_rhs
  TZ_rhs(i,j,k)  = toTZ_rhs
  Gamx_rhs(i,j,k) = toGams_rhs*vx(i,j,k)+toGam1_rhs*ux(i,j,k)+toGam2_rhs*wx(i,j,k)
  Gamy_rhs(i,j,k) = toGams_rhs*vy(i,j,k)+toGam1_rhs*uy(i,j,k)+toGam2_rhs*wy(i,j,k)
  Gamz_rhs(i,j,k) = toGams_rhs*vz(i,j,k)+toGam1_rhs*uz(i,j,k)+toGam2_rhs*wz(i,j,k)
  Axx_rhs(i,j,k) = (ulx(i,j,k)*ulx(i,j,k)-0.5d0*qxx(i,j,k))*toA11_rhs+(wlx(i,j,k)*wlx(i,j,k)-0.5d0*qxx(i,j,k))*toA22_rhs+ulx(i,j,k)*wlx(i,j,k)*toA12_rhs &
                  + ulx(i,j,k)*slx(i,j,k)*toAs1_rhs+wlx(i,j,k)*slx(i,j,k)*toAs2_rhs+slx(i,j,k)*slx(i,j,k)*toAss_rhs &
                  + 0.5d0*qxx(i,j,k)*toAqq_rhs
  Ayy_rhs(i,j,k) = (uly(i,j,k)*uly(i,j,k)-0.5d0*qyy(i,j,k))*toA11_rhs+(wly(i,j,k)*wly(i,j,k)-0.5d0*qyy(i,j,k))*toA22_rhs+uly(i,j,k)*wly(i,j,k)*toA12_rhs &
                  + uly(i,j,k)*sly(i,j,k)*toAs1_rhs+wly(i,j,k)*sly(i,j,k)*toAs2_rhs+sly(i,j,k)*sly(i,j,k)*toAss_rhs &
                  + 0.5d0*qyy(i,j,k)*toAqq_rhs
  Azz_rhs(i,j,k) = (ulz(i,j,k)*ulz(i,j,k)-0.5d0*qzz(i,j,k))*toA11_rhs+(wlz(i,j,k)*wlz(i,j,k)-0.5d0*qzz(i,j,k))*toA22_rhs+ulz(i,j,k)*wlz(i,j,k)*toA12_rhs &
                  + ulz(i,j,k)*slz(i,j,k)*toAs1_rhs+wlz(i,j,k)*slz(i,j,k)*toAs2_rhs+slz(i,j,k)*slz(i,j,k)*toAss_rhs &
                  + 0.5d0*qzz(i,j,k)*toAqq_rhs
  Axy_rhs(i,j,k) = (ulx(i,j,k)*uly(i,j,k)-0.5d0*qxy(i,j,k))*toA11_rhs+(wlx(i,j,k)*wly(i,j,k)-0.5d0*qxy(i,j,k))*toA22_rhs+ &
                   (ulx(i,j,k)*wly(i,j,k)+uly(i,j,k)*wlx(i,j,k))/TWO*toA12_rhs &
                  +(ulx(i,j,k)*sly(i,j,k)+uly(i,j,k)*slx(i,j,k))/TWO*toAs1_rhs &
                  +(wlx(i,j,k)*sly(i,j,k)+wly(i,j,k)*slx(i,j,k))/TWO*toAs2_rhs &
                  +(slx(i,j,k)*sly(i,j,k)+sly(i,j,k)*slx(i,j,k))/TWO*toAss_rhs &
                  + 0.5d0*qxy(i,j,k)*toAqq_rhs
  Axz_rhs(i,j,k) = (ulx(i,j,k)*ulz(i,j,k)-0.5d0*qxz(i,j,k))*toA11_rhs+(wlx(i,j,k)*wlz(i,j,k)-0.5d0*qxz(i,j,k))*toA22_rhs+ &
                   (ulx(i,j,k)*wlz(i,j,k)+ulz(i,j,k)*wlx(i,j,k))/TWO*toA12_rhs &
                  +(ulx(i,j,k)*slz(i,j,k)+ulz(i,j,k)*slx(i,j,k))/TWO*toAs1_rhs &
                  +(wlx(i,j,k)*slz(i,j,k)+wlz(i,j,k)*slx(i,j,k))/TWO*toAs2_rhs &
                  +(slx(i,j,k)*slz(i,j,k)+slz(i,j,k)*slx(i,j,k))/TWO*toAss_rhs &
                  + 0.5d0*qxz(i,j,k)*toAqq_rhs
  Ayz_rhs(i,j,k) = (uly(i,j,k)*ulz(i,j,k)-0.5d0*qyz(i,j,k))*toA11_rhs+(wlz(i,j,k)*wlz(i,j,k)-0.5d0*qyz(i,j,k))*toA22_rhs+ &
                   (uly(i,j,k)*wlz(i,j,k)+ulz(i,j,k)*wly(i,j,k))/TWO*toA12_rhs &
                  +(uly(i,j,k)*slz(i,j,k)+ulz(i,j,k)*sly(i,j,k))/TWO*toAs1_rhs &
                  +(wly(i,j,k)*slz(i,j,k)+wlz(i,j,k)*sly(i,j,k))/TWO*toAs2_rhs &
                  +(sly(i,j,k)*slz(i,j,k)+slz(i,j,k)*sly(i,j,k))/TWO*toAss_rhs &
                  + 0.5d0*qyz(i,j,k)*toAqq_rhs
      enddo
     enddo
    enddo

  endif

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

! NOTE: we need Kreiss-Oliger dissipation here instead of rhs calculation routine
  if(eps>0)then 
! usual Kreiss-Oliger dissipation      
  call kodis_sh(ex,crho,sigma,R,chi,chi_rhs,SSS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,trK,trK_rhs,SSS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,dxx,gxx_rhs,SSS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,gxy,gxy_rhs,AAS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,gxz,gxz_rhs,ASA,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,dyy,gyy_rhs,SSS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,gyz,gyz_rhs,SAA,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,dzz,gzz_rhs,SSS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,Axx,Axx_rhs,SSS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,Axy,Axy_rhs,AAS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,Axz,Axz_rhs,ASA,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,Ayy,Ayy_rhs,SSS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,Ayz,Ayz_rhs,SAA,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,Azz,Azz_rhs,SSS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,Gamx,Gamx_rhs,ASS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,Gamy,Gamy_rhs,SAS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,Gamz,Gamz_rhs,SSA,Symmetry,eps,sst)

  call kodis_sh(ex,crho,sigma,R,Lap,Lap_rhs,SSS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,betax,betax_rhs,ASS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,betay,betay_rhs,SAS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,betaz,betaz_rhs,SSA,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,dtSfx,dtSfx_rhs,ASS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,dtSfy,dtSfy_rhs,SAS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,dtSfz,dtSfz_rhs,SSA,Symmetry,eps,sst)

  call kodis_sh(ex,crho,sigma,R,TZ,TZ_rhs,SSS,Symmetry,eps,sst)
  endif

  return

  end subroutine david_milton_cpbc_ss
#endif  
! repopulate the buffer points of outer boundary through extroplation
! need CPBC_ghost_width
  subroutine repo_extro_ss(ex,x,y,z,f,zmin,zmax,tpp) 
  implicit none
  integer,intent(in ):: ex(1:3)
  double precision,intent(in),dimension(ex(1))::x
  double precision,intent(in),dimension(ex(2))::y
  double precision,intent(in),dimension(ex(3))::z
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: f
  real*8,  intent(in):: zmin,zmax
! extraplate type  
! 0: Lagange polynomial; 1: D+^n f = 0
  integer,intent(in) :: tpp
!~~~~~~~~~~~> local variables
  logical :: gont
  real*8 :: dZ
  integer :: i, j, k
  integer :: layer(1:6,1:6),gp
  real*8 :: extroplate_lag,extroplate_cg

  integer :: NP

!sanity check
  if(ex(3) .le. CPBC_ghost_width +(ghost_width*2+1))then
    write(*,*) "repo_extro_ss has assumed ex(3) > CPBC_ghost_width +(ghost_width*2+1) but ex(3) = ",ex(3),"CPBC_ghost_width = ",CPBC_ghost_width
    stop
  endif

  dZ = Z(2) - Z(1)

layer(1:3,:) = 1
layer(4:6,:) =-1
 
if(dabs(Z(ex(3))-zmax) < dZ)then
  layer(1,3) = 1
  layer(2,3) = 1
  layer(3,3) = ex(3) - CPBC_ghost_width
  layer(4,3) = ex(1)
  layer(5,3) = ex(2)
  layer(6,3) = ex(3) - CPBC_ghost_width
endif
! extroplate point by point
  gp = 3

  gont = any( layer(:,gp) == - 1 )
 
  if( .not. gont ) then

!!! fixme: note the assumption points requirement is enough or not          
     select case (tpp)
     case (0)             
       NP = ghost_width*2+1
!      NP = ghost_width*2-1

    do k = layer(3,gp) + 1,ex(3)
     do j = layer(2,gp), layer(5,gp)
      do i = layer(1,gp), layer(4,gp)
         f(i,j,k) = extroplate_lag(NP,f(i,j,k-NP:k-1))
      enddo
     enddo
    enddo

     case (1)
!       NP = (ghost_width-1)*2
       NP = ghost_width*2

    do k = layer(3,gp) + 1,ex(3)
     do j = layer(2,gp), layer(5,gp)
      do i = layer(1,gp), layer(4,gp)
         f(i,j,k) = extroplate_cg(NP,f(i,j,k-NP:k-1))
      enddo
     enddo
    enddo

     case (2)
       NP = ghost_width*2+1
!      NP = ghost_width*2-1

       NP = NP + CPBC_ghost_width
     do j = layer(2,gp), layer(5,gp)
      do i = layer(1,gp), layer(4,gp)
         call extroplate_lag2(NP,f(i,j,ex(3)-NP+1:ex(3)))
      enddo
     enddo

     case default
             write(*,*) "repo_extro_ss: not recognized extraplation type = ",tpp
      return
     end select


  endif

  return

  end subroutine repo_extro_ss
! extroplate for unigrid with Lagange polynomial
  function extroplate_lag(N,f)  result(gont)             
  implicit none
  integer,intent(in ) :: N
  real*8,dimension(N),intent(in) :: f

  real*8 :: gont

  real*8,parameter :: THR=3.d0
  real*8,parameter :: FIV=5.d0,TEN=1.d1,NIN=9.d0
  real*8,parameter :: SEV=7.d0,TYO=2.1d1,F35=3.5d1
  real*8,parameter :: F36=3.6d1,F84=8.4d1,F126=1.26d2
  real*8,parameter :: F11=1.1d1,F55=5.5d1,F165=1.65d2,F330=3.3d2,F462=4.62d2

! Lagange polynomial
  select case (N)
! for 2nd order code
  case (3)
     gont = THR*f(3)-THR*f(2)+f(1)
! for 2nd order code
  case (5)
     gont = FIV*f(5)-TEN*f(4)+TEN*f(3)-FIV*f(2)+f(1)
! for 4th order code
  case (7)
     gont = SEV*f(7)-TYO*f(6)+F35*f(5)-F35*f(4)+TYO*f(3)-SEV*f(2)+f(1)
! for 6th order code
  case (9)
     gont = NIN*f(9)-F36*f(8)+F84*f(7)-F126*f(6)+F126*f(5)-F84*f(4)+F36*f(3)-NIN*f(2)+f(1)
! for 8th order code
  case (11)
     gont = F11*f(11)-F55*f(10)+F165*f(9)-F330*f(8)+F462*f(7)-F462*f(6)+F330*f(5)-F165*f(4)+F55*f(3)-F11*f(2)+f(1)
  end select

  return

  end function extroplate_lag
! extroplate for unigrid with Lagange polynomial
! but using inner N-ghost_width points for all of the outer ghost_width points
  subroutine extroplate_lag2(N,f)
  implicit none
  integer,intent(in ) :: N
  real*8,dimension(N),intent(inout) :: f

  integer :: NI,i
  real*8 :: s1,s2

  NI = N - CPBC_ghost_width

  do i=1,CPBC_ghost_width

! Lagange polynomial
    select case (NI)
! for 2nd order code
    case (3)
        f(NI+i) = i**2*f(1)/2+i*f(1)/2-i**2*f(2)-2*i*f(2)+f(3)*i**2/2+3.D0/2.D0*f(3)*i+f(3)
! for 2nd order code
    case (5)
        f(NI+i) = i**4*f(1)/24+i**3*f(1)/4+11.D0/24.D0*i**2*f(1)+i*f(1)/4-i**4*f(2)/6  &
                 -7.D0/6.D0*i**3*f(2)-7.D0/3.D0*i**2*f(2)-4.D0/3.D0*i*f(2)+f(3)*i**4/4 &
                 +2*f(3)*i**3+19.D0/4.D0*f(3)*i**2+3*f(3)*i-i**4*f(4)/6                &
                 -3.D0/2.D0*i**3*f(4)-13.D0/3.D0*i**2*f(4)-4*i*f(4)                    &
                 +f(5)*i**4/24+5.D0/12.D0*f(5)*i**3+35.D0/24.D0*f(5)*i**2              &
                 +25.D0/12.D0*f(5)*i+f(5)       
! for 4th order code
    case (7)
       s1 = 33.D0/4.D0*f(3)*i**2+15.D0/2.D0*f(5)*i+117.D0/8.D0*f(5)*i**2         &
           -121.D0/36.D0*i**4*f(4)+i**5*f(1)/48-i**6*f(2)/120-20.D0/3.D0*i*f(4)  &
           +35.D0/144.D0*f(7)*i**4+137.D0/48.D0*f(5)*i**4+107.D0/48.D0*f(3)*i**4 &
           +203.D0/90.D0*f(7)*i**2-31.D0/3.D0*i**3*f(4)-27.D0/10.D0*i**2*f(2)    &
           +f(5)*i**6/48+15.D0/4.D0*f(3)*i+17.D0/144.D0*i**4*f(1)-i**5*f(4)/2    &
           -i**5*f(6)/6-13.D0/6.D0*i**3*f(2)+137.D0/360.D0*i**2*f(1)             &
           +49.D0/48.D0*f(7)*i**3
       f(NI+i) = s1-19.D0/24.D0*i**4*f(2)+7.D0/240.D0*f(7)*i**5+f(7)*i**6/720    &
                +461.D0/48.D0*f(5)*i**3-6*i*f(6)-87.D0/10.D0*i**2*f(6)           &
                +i**6*f(1)/720-127.D0/9.D0*i**2*f(4)+19.D0/48.D0*f(5)*i**5       &
                +5.D0/16.D0*i**3*f(1)-i**6*f(6)/120+49.D0/20.D0*f(7)*i           &
                -29.D0/6.D0*i**3*f(6)-31.D0/24.D0*i**4*f(6)-6.D0/5.D0*i*f(2)     &
                +i*f(1)/6-2.D0/15.D0*i**5*f(2)-i**6*f(4)/36                      &
                +307.D0/48.D0*f(3)*i**3+17.D0/48.D0*f(3)*i**5+f(3)*i**6/48+f(7)
! for 6th order code
    case (9)
       s1 = -8*i*f(8)-527.D0/180.D0*i**3*f(2)+2803.D0/480.D0*f(3)*i**4           &
            -i**8*f(8)/5040+18353.D0/720.D0*f(7)*i**3-391.D0/720.D0*i**6*f(4)    &
            +1457.D0/36.D0*f(5)*i**3+9.D0/80.D0*f(9)*i**5+23.D0/2880.D0*i**6*f(1) &
            +17.D0/720.D0*f(7)*i**7-56.D0/3.D0*i*f(6)+761.D0/280.D0*f(9)*i        &
            -67.D0/45.D0*i**4*f(2)+14*f(7)*i+f(3)*i**7/48-268.D0/15.D0*i**4*f(6)  &
            -2003.D0/45.D0*i**2*f(6)+13.D0/960.D0*f(9)*i**6-73.D0/720.D0*i**6*f(8) &
            +f(7)*i**8/1440+363.D0/1120.D0*i**2*f(1)-797.D0/20.D0*i**3*f(6)       &
            -11.D0/240.D0*i**7*f(6)-329.D0/90.D0*i**4*f(8)+179.D0/36.D0*f(5)*i**5 &
            +967.D0/5760.D0*i**4*f(1)-103.D0/35.D0*i**2*f(2)-481.D0/35.D0*i**2*f(8) &
            +179.D0/72.D0*f(7)*i**5-349.D0/36.D0*i**3*f(8)+61.D0/240.D0*f(3)*i**6  &
            -115.D0/144.D0*i**5*f(8)+f(5)*i**8/576-56.D0/5.D0*i*f(4)               &
            +187.D0/16.D0*f(3)*i**3-149.D0/240.D0*i**6*f(6)
       f(NI+i) = s1+i**8*f(1)/40320+2143.D0/180.D0*f(3)*i**2+469.D0/1440.D0*i**3*f(1) &
               +f(5)*i**7/18+621.D0/20.D0*f(7)*i**2+f(9)+267.D0/160.D0*f(9)*i**3     &
               +7.D0/144.D0*i**5*f(1)-i**7*f(8)/144+1069.D0/1920.D0*f(9)*i**4        &
               -141.D0/5.D0*i**2*f(4)-29.D0/5040.D0*i**7*f(2)-i**8*f(4)/720          &
               +691.D0/16.D0*f(5)*i**2+f(9)*i**7/1120-2581.D0/720.D0*i**5*f(4)       &
               -i**8*f(2)/5040+10993.D0/576.D0*f(5)*i**4+14.D0/3.D0*f(3)*i           &
               -i**8*f(6)/720-4891.D0/180.D0*i**3*f(4)+13.D0/8.D0*f(3)*i**5          &
               +239.D0/720.D0*f(7)*i**6+15289.D0/1440.D0*f(7)*i**4+f(3)*i**8/1440    &
               -49.D0/720.D0*i**6*f(2)+35.D0/2.D0*f(5)*i-71.D0/16.D0*i**5*f(6)       &
               +f(9)*i**8/40320+29531.D0/10080.D0*f(9)*i**2+209.D0/288.D0*f(5)*i**6  &
               -1193.D0/90.D0*i**4*f(4)+i*f(1)/8-61.D0/144.D0*i**5*f(2)+i**7*f(1)/1440 &
               -31.D0/720.D0*i**7*f(4)-8.D0/7.D0*i*f(2)
! for 8th order code
    case (11)
       s2 = -433739.D0/7560.D0*i**4*f(8)+7129.D0/25200.D0*i**2*f(1)         &
            -6947.D0/8640.D0*i**5*f(2)+3013.D0/172800.D0*i**6*f(1)          &
            -107.D0/1440.D0*i**8*f(6)-119.D0/4320.D0*i**7*f(2)+i*f(1)/10    &
            +f(7)*i**10/17280+i**10*f(1)/3628800+59.D0/5040.D0*f(3)*i**8    &
            -67.D0/1440.D0*i**7*f(10)+105.D0/2.D0*f(7)*i                    &
            +84095.D0/36288.D0*f(11)*i**3-62549.D0/720.D0*i**4*f(6)         &
            -263.D0/84.D0*i**2*f(2)-5419.D0/1440.D0*i**6*f(8)               &
            +11.D0/30240.D0*f(11)*i**8+757.D0/5760.D0*f(3)*i**7             &
            +39867.D0/2240.D0*f(3)*i**3-i**10*f(8)/30240-8.D0/9.D0*i**7*f(6) &
            -10*i*f(10)+6961.D0/72.D0*f(5)*i**2-i**10*f(10)/362880          &
            +i**9*f(1)/80640+728587.D0/8640.D0*f(7)*i**4-41.D0/1260.D0*i**8*f(4)
      s1 = s2-6709.D0/17280.D0*i**6*f(10)+47.D0/80640.D0*f(3)*i**9          &
          +49.D0/17280.D0*f(5)*i**9-1253.D0/480.D0*i**6*f(4)                &
          +6751.D0/48.D0*f(7)*i**2+10427.D0/11520.D0*f(3)*i**6-10.D0/9.D0*i*f(2) &
          -23.D0/181440.D0*i**9*f(2)-i**9*f(10)/6720+45449.D0/11520.D0*f(3)*i**5 &
          +2281.D0/2880.D0*f(7)*i**7-252.D0/5.D0*i*f(6)+f(9)*i**10/80640         &
          +29.D0/120960.D0*i**8*f(1)-6541.D0/63.D0*i**2*f(8)                     &
          +461789.D0/4320.D0*f(5)*i**3-161353.D0/45360.D0*i**3*f(2)              &
          +435893.D0/40320.D0*f(3)*i**4-97.D0/2520.D0*i**8*f(8)                  &
          +1123.D0/5760.D0*f(9)*i**7-i**10*f(4)/30240-1003.D0/21.D0*i**2*f(4)    &
          -4861.D0/252.D0*i**2*f(10)-40*i*f(8)-i**9*f(6)/288-i**10*f(6)/14400    &
          -13.D0/7560.D0*i**9*f(8)+1303.D0/4032.D0*i**3*f(1)
      s2 = s1-151.D0/60480.D0*i**8*f(2)+19.D0/256.D0*i**5*f(1)                   &
          +71689.D0/480.D0*f(7)*i**3-6877.D0/50.D0*i**2*f(6)+34343.D0/5760.D0*f(7)*i**6 &
          -211.D0/60480.D0*i**8*f(10)+129067.D0/5760.D0*f(5)*i**5+35*f(5)*i     &
          -8321.D0/720.D0*i**5*f(4)+i**7*f(1)/384-1197.D0/8.D0*i**3*f(6)        &
          +3533.D0/224.D0*f(3)*i**2+18047.D0/11520.D0*f(9)*i**6                 &
          +163313.D0/5760.D0*f(7)*i**5+f(5)*i**10/17280+f(11)*i**10/3628800     &
          +11.D0/725760.D0*f(11)*i**9-120.D0/7.D0*i*f(4)-i**9*f(4)/630          &
          +177133.D0/50400.D0*f(11)*i**2-93773.D0/14400.D0*i**6*f(6)            &
          +121.D0/24192.D0*f(11)*i**7+28603.D0/5760.D0*f(5)*i**6                &
          -197741.D0/90720.D0*i**4*f(2)+f(3)*i**10/80640+7381.D0/2520.D0*f(11)*i-3229.D0/17280.D0*i**6*f(2)
     f(NI+i) = s2+17.D0/5760.D0*f(7)*i**9-22439.D0/420.D0*i**3*f(4)             &
              -1877.D0/5040.D0*i**7*f(4)-i**10*f(2)/362880-43319.D0/1440.D0*i**5*f(6) &
              +31.D0/480.D0*f(7)*i**8+1999.D0/2880.D0*f(5)*i**7+45.D0/8.D0*f(3)*i     &
              +19.D0/320.D0*f(5)*i**8-349.D0/720.D0*i**7*f(8)                         &
              +273431.D0/4320.D0*f(5)*i**4-242639.D0/7560.D0*i**4*f(4)                &
              +4523.D0/22680.D0*i**4*f(1)+607.D0/40320.D0*f(9)*i**8                   &
              +92771.D0/11520.D0*f(9)*i**5+264767.D0/10080.D0*f(9)*i**4               &
              +115923.D0/2240.D0*f(9)*i**3+6121.D0/112.D0*f(9)*i**2+45.D0/2.D0*f(9)*i &
              +53.D0/80640.D0*f(9)*i**9-6041.D0/2880.D0*i**5*f(10)                    &
              -663941.D0/90720.D0*i**4*f(10)-79913.D0/5040.D0*i**3*f(10)              &
              +7513.D0/172800.D0*f(11)*i**6+8591.D0/34560.D0*f(11)*i**5               &
              +341693.D0/362880.D0*f(11)*i**4-13349.D0/720.D0*i**5*f(8)               &
              -400579.D0/3780.D0*i**3*f(8)+f(11)

    end select

  enddo

  return

  end subroutine extroplate_lag2
! extroplate for unigrid with Calabrese Gundlach type, Eq.(16) of CQG 23 S343 (2006)
  function extroplate_cg(N,f)  result(gont)             
  implicit none
  integer,intent(in ) :: N
  real*8,dimension(N),intent(in) :: f

  real*8 :: gont

! Eq.(16) of CQG 23 S343 (2006)
  select case (N)
! for 2nd order code
  case (2)
     gont = 2.d0*f(2)-f(1)
! for 4th order code
  case (4)
     gont = 4.d0*f(4)-6.d0*f(3)+4.d0*f(2)-f(1)
! for 6th order code
  case (6)
! Eq.(C7) of PRD 83, 024025          
     gont = 6.d0*f(6)-1.5d1*f(5)+2.d1*f(4)-1.5d1*f(3)+6.d0*f(2)-f(1)
! for 8th order code
  case (8)
     gont = 8.d0*f(8)-2.8d1*f(7)+5.6d1*f(6)-7.d1*f(5)+5.6d1*f(4)-2.8d1*f(3)+8.d0*f(2)-f(1)
  end select

  return

  end function extroplate_cg
! need CPBC_ghost_width
  subroutine david_milton_extroplate_ss(ex,crho,sigma,R,                       &
               TZ,chi,trK,                                                     &
               dxx,gxy,gxz,dyy,gyz,dzz,                                        &
               Axx,Axy,Axz,Ayy,Ayz,Azz,                                        &
               Gmx,Gmy,Gmz,                                                    &
               Lap    ,  betax   ,  betay   ,  betaz   ,                       &
               dtSfx  ,  dtSfy   ,  dtSfz,zmin,zmax) 

! NOTE: we need Kreiss-Oliger dissipation here instead of rhs calculation routine
  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3)
  double precision,intent(in),dimension(ex(1))::crho
  double precision,intent(in),dimension(ex(2))::sigma
  double precision,intent(in),dimension(ex(3))::R
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: TZ,chi,dxx,dyy,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: trK
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: gxy,gxz,gyz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Gmx,Gmy,Gmz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Lap, betax, betay, betaz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: dtSfx,  dtSfy,  dtSfz
  real*8,  intent(in):: zmin,zmax

#define tptype 1
#if (tptype == 0)
! default we always use hp (tpp=0)

  call repo_extro_ss(ex,crho,sigma,R,TZ,zmin,zmax,0) 
  call repo_extro_ss(ex,crho,sigma,R,chi,zmin,zmax,0) 
  call repo_extro_ss(ex,crho,sigma,R,dxx,zmin,zmax,0) 
  call repo_extro_ss(ex,crho,sigma,R,dyy,zmin,zmax,0) 
  call repo_extro_ss(ex,crho,sigma,R,dzz,zmin,zmax,0) 
  call repo_extro_ss(ex,crho,sigma,R,gxy,zmin,zmax,0) 
  call repo_extro_ss(ex,crho,sigma,R,gxz,zmin,zmax,0) 
  call repo_extro_ss(ex,crho,sigma,R,gyz,zmin,zmax,0) 
  call repo_extro_ss(ex,crho,sigma,R,trK,zmin,zmax,0) 
  call repo_extro_ss(ex,crho,sigma,R,Axx,zmin,zmax,0) 
  call repo_extro_ss(ex,crho,sigma,R,Ayy,zmin,zmax,0) 
  call repo_extro_ss(ex,crho,sigma,R,Azz,zmin,zmax,0) 
  call repo_extro_ss(ex,crho,sigma,R,Axy,zmin,zmax,0) 
  call repo_extro_ss(ex,crho,sigma,R,Axz,zmin,zmax,0) 
  call repo_extro_ss(ex,crho,sigma,R,Ayz,zmin,zmax,0) 
  call repo_extro_ss(ex,crho,sigma,R,Gmx,zmin,zmax,0) 
  call repo_extro_ss(ex,crho,sigma,R,Gmy,zmin,zmax,0) 
  call repo_extro_ss(ex,crho,sigma,R,Gmz,zmin,zmax,0) 
  call repo_extro_ss(ex,crho,sigma,R,Lap,zmin,zmax,0) 
  call repo_extro_ss(ex,crho,sigma,R,betax,zmin,zmax,0) 
  call repo_extro_ss(ex,crho,sigma,R,betay,zmin,zmax,0) 
  call repo_extro_ss(ex,crho,sigma,R,betaz,zmin,zmax,0) 
  call repo_extro_ss(ex,crho,sigma,R,dtSfx,zmin,zmax,0) 
  call repo_extro_ss(ex,crho,sigma,R,dtSfy,zmin,zmax,0) 
  call repo_extro_ss(ex,crho,sigma,R,dtSfz,zmin,zmax,0) 
#elif (tptype == 1)
! all D+ f = 0 (tpp=1)

  call repo_extro_ss(ex,crho,sigma,R,TZ,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,chi,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,dxx,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,dyy,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,dzz,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,gxy,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,gxz,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,gyz,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,trK,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,Axx,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,Ayy,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,Azz,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,Axy,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,Axz,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,Ayz,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,Gmx,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,Gmy,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,Gmz,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,Lap,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,betax,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,betay,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,betaz,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,dtSfx,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,dtSfy,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,dtSfz,zmin,zmax,1) 
#elif (tptype == 2)
! Lagange polynomial but all used inner points (tpp=2)

  call repo_extro_ss(ex,crho,sigma,R,TZ,zmin,zmax,2) 
  call repo_extro_ss(ex,crho,sigma,R,chi,zmin,zmax,2) 
  call repo_extro_ss(ex,crho,sigma,R,dxx,zmin,zmax,2) 
  call repo_extro_ss(ex,crho,sigma,R,dyy,zmin,zmax,2) 
  call repo_extro_ss(ex,crho,sigma,R,dzz,zmin,zmax,2) 
  call repo_extro_ss(ex,crho,sigma,R,gxy,zmin,zmax,2) 
  call repo_extro_ss(ex,crho,sigma,R,gxz,zmin,zmax,2) 
  call repo_extro_ss(ex,crho,sigma,R,gyz,zmin,zmax,2) 
  call repo_extro_ss(ex,crho,sigma,R,trK,zmin,zmax,2) 
  call repo_extro_ss(ex,crho,sigma,R,Axx,zmin,zmax,2) 
  call repo_extro_ss(ex,crho,sigma,R,Ayy,zmin,zmax,2) 
  call repo_extro_ss(ex,crho,sigma,R,Azz,zmin,zmax,2) 
  call repo_extro_ss(ex,crho,sigma,R,Axy,zmin,zmax,2) 
  call repo_extro_ss(ex,crho,sigma,R,Axz,zmin,zmax,2) 
  call repo_extro_ss(ex,crho,sigma,R,Ayz,zmin,zmax,2) 
  call repo_extro_ss(ex,crho,sigma,R,Gmx,zmin,zmax,2) 
  call repo_extro_ss(ex,crho,sigma,R,Gmy,zmin,zmax,2) 
  call repo_extro_ss(ex,crho,sigma,R,Gmz,zmin,zmax,2) 
  call repo_extro_ss(ex,crho,sigma,R,Lap,zmin,zmax,2) 
  call repo_extro_ss(ex,crho,sigma,R,betax,zmin,zmax,2) 
  call repo_extro_ss(ex,crho,sigma,R,betay,zmin,zmax,2) 
  call repo_extro_ss(ex,crho,sigma,R,betaz,zmin,zmax,2) 
  call repo_extro_ss(ex,crho,sigma,R,dtSfx,zmin,zmax,2) 
  call repo_extro_ss(ex,crho,sigma,R,dtSfy,zmin,zmax,2) 
  call repo_extro_ss(ex,crho,sigma,R,dtSfz,zmin,zmax,2) 

#elif (tptype == 3)
! thumb of rule: D+ f = 0 (tpp=1) for outgoing ones; hp (tpp=0) for ingoing ones

  call repo_extro_ss(ex,crho,sigma,R,TZ,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,chi,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,dxx,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,dyy,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,dzz,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,gxy,zmin,zmax,0) 
  call repo_extro_ss(ex,crho,sigma,R,gxz,zmin,zmax,0) 
  call repo_extro_ss(ex,crho,sigma,R,gyz,zmin,zmax,0) 
  call repo_extro_ss(ex,crho,sigma,R,trK,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,Axx,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,Ayy,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,Azz,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,Axy,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,Axz,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,Ayz,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,Gmx,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,Gmy,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,Gmz,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,Lap,zmin,zmax,0) 
  call repo_extro_ss(ex,crho,sigma,R,betax,zmin,zmax,0) 
  call repo_extro_ss(ex,crho,sigma,R,betay,zmin,zmax,0) 
  call repo_extro_ss(ex,crho,sigma,R,betaz,zmin,zmax,0) 
  call repo_extro_ss(ex,crho,sigma,R,dtSfx,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,dtSfy,zmin,zmax,1) 
  call repo_extro_ss(ex,crho,sigma,R,dtSfz,zmin,zmax,1) 

#else  
#error "not recognized tptype"
#endif

#undef tptype

  return

  end subroutine david_milton_extroplate_ss
!construct rACqq rhs
  subroutine cpbcrACqq(ex,crho,sigma,R,x,y,z,                                  &
               drhodx, drhody, drhodz,                                         &
               dsigmadx,dsigmady,dsigmadz,                                     &
               dRdx,dRdy,dRdz,                                                 &
               drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                &
               dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,    &
               dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,                            &
               xmin,ymin,zmin,xmax,ymax,zmax,rACqq,&
               chi,dxx,gxy,gxz,dyy,gyz,dzz, &
               Lap,Sfx,Sfy,Sfz,Axx,Axy,Axz,Ayy,Ayz,Azz,rACss,Symmetry,sst)

  implicit none
 
!~~~~~~> Input parameters:
  integer, intent(in):: ex(1:3),Symmetry,sst
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
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: chi,dxx,gxy,gxz,dyy,gyz,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Lap,Sfx,Sfy,Sfz,Axx,Axy,Axz,Ayy,Ayz,Azz,rACss
  real*8,  intent(in):: xmin,ymin,zmin,xmax,ymax,zmax
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout)::rACqq
!~~~~~~> Other variables:
  real*8 :: chin1,alpha,gxx,gyy,gzz
  real*8 :: sfxx,sfxy,sfxz,sfyx,sfyy,sfyz,sfzx,sfzy,sfzz
  real*8 :: gxxx,gxxy,gxxz
  real*8 :: gxyx,gxyy,gxyz
  real*8 :: gxzx,gxzy,gxzz
  real*8 :: gyyx,gyyy,gyyz
  real*8 :: gyzx,gyzy,gyzz
  real*8 :: gzzx,gzzy,gzzz
  logical :: gont
  integer :: i, j, k
  integer :: layer(1:6,1:6),gp
! index of layer, first one: i,j,k; second one: front back etc. boundary
  integer :: kmin,kmax
  real*8, parameter :: ZEO = 0.d0, ONE = 1.d0, TWO=2.d0
  real*8 :: dR
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0

  dR = R(2)-R(1)
  
  kmax = ex(3)

  kmin = 1

layer(1:3,:) = 1
layer(4:6,:) =-1
 
if(dabs(R(ex(3))-zmax) < dR)then
  layer(1,3) = 1
  layer(2,3) = 1
  layer(4,3) = ex(1)
  layer(5,3) = ex(2)
! consider buffer points near boundary 
  layer(3,3) = ex(3) - CPBC_ghost_width
  layer(6,3) = ex(3) - CPBC_ghost_width
endif

  gp = 3

  gont = any( layer(:,gp) == - 1 )
 
   if( .not. gont ) then

    do k = layer(3,gp), layer(6,gp)
     do j = layer(2,gp), layer(5,gp)
      do i = layer(1,gp), layer(4,gp)
         alpha = Lap(i,j,k)+1.d0
         chin1 = chi(i,j,k)+1.d0
         gxx = dxx(i,j,k)+1.d0
         gyy = dyy(i,j,k)+1.d0
         gzz = dzz(i,j,k)+1.d0
         call point_fderivs_shc(ex,Sfx,sfxx,sfxy,sfxz,crho,sigma,R,ANTI,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Sfy,sfyx,sfyy,sfyz,crho,sigma,R,SYM,ANTI,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Sfz,sfzx,sfzy,sfzz,crho,sigma,R,SYM,SYM,ANTI,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,dxx,gxxx,gxxy,gxxz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,gxy,gxyx,gxyy,gxyz,crho,sigma,R,ANTI,ANTI,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,gxz,gxzx,gxzy,gxzz,crho,sigma,R,ANTI,SYM,ANTI,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,dyy,gyyx,gyyy,gyyz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,gyz,gyzx,gyzy,gyzz,crho,sigma,R,SYM,ANTI,ANTI,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,dzz,gzzx,gzzy,gzzz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call racqq_point(Axx(i,j,k),Axy(i,j,k),Axz(i,j,k),Ayy(i,j,k),Ayz(i,j,k),Azz(i,j,k), &
                    alpha,Sfx(i,j,k),Sfy(i,j,k),Sfz(i,j,k),chin1, &
                    sfxx,sfyx,sfzx,sfxy,sfyy,sfzy,sfxz,sfyz,sfzz, &
                    gxxx,gxyx,gxzx,gyyx,gyzx,gzzx, &
                    gxxy,gxyy,gxzy,gyyy,gyzy,gzzy, &
                    gxxz,gxyz,gxzz,gyyz,gyzz,gzzz, &
                    gxx,gxy(i,j,k),gxz(i,j,k),gyy,gyz(i,j,k),gzz, &
                    rACqq(i,j,k),rACss(i,j,k))
      enddo
     enddo
    enddo
   endif

  return

  end subroutine cpbcrACqq
!construct rtrK rhs
  subroutine cpbcrtrK(ex,crho,sigma,R,x,y,z,                                  &
               drhodx, drhody, drhodz,                                         &
               dsigmadx,dsigmady,dsigmadz,                                     &
               dRdx,dRdy,dRdz,                                                 &
               drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                &
               dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,    &
               dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,                            &
               xmin,ymin,zmin,xmax,ymax,zmax,rtrK,&
               chi,trK,dxx,gxy,gxz,dyy,gyz,dzz, &
               Lap,Sfx,Sfy,Sfz,TZ,Symmetry,sst,kappa1,kappa2)

  implicit none
 
!~~~~~~> Input parameters:
  integer, intent(in):: ex(1:3),Symmetry,sst
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
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: chi,trK,dxx,gxy,gxz,dyy,gyz,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Lap,Sfx,Sfy,Sfz,TZ
  real*8,  intent(in):: xmin,ymin,zmin,xmax,ymax,zmax,kappa1,kappa2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout)::rtrK
!~~~~~~> Other variables:
  real*8 :: chin1,alpha,gxx,gyy,gzz
  real*8 :: Kx,Ky,Kz,TZx,TZy,TZz
  logical :: gont
  integer :: i, j, k
  integer :: layer(1:6,1:6),gp
! index of layer, first one: i,j,k; second one: front back etc. boundary
  integer :: kmin,kmax
  real*8, parameter :: ZEO = 0.d0, ONE = 1.d0, TWO=2.d0
  real*8 :: dR
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0

  dR = R(2)-R(1)
  
  kmax = ex(3)

  kmin = 1

layer(1:3,:) = 1
layer(4:6,:) =-1
 
if(dabs(R(ex(3))-zmax) < dR)then
  layer(1,3) = 1
  layer(2,3) = 1
  layer(4,3) = ex(1)
  layer(5,3) = ex(2)
! consider buffer points near boundary 
  layer(3,3) = ex(3) - CPBC_ghost_width
  layer(6,3) = ex(3) - CPBC_ghost_width
endif

  gp = 3

  gont = any( layer(:,gp) == - 1 )
 
   if( .not. gont ) then

    do k = layer(3,gp), layer(6,gp)
     do j = layer(2,gp), layer(5,gp)
      do i = layer(1,gp), layer(4,gp)
         alpha = Lap(i,j,k)+1.d0
         chin1 = chi(i,j,k)+1.d0
         gxx = dxx(i,j,k)+1.d0
         gyy = dyy(i,j,k)+1.d0
         gzz = dzz(i,j,k)+1.d0
         call point_fderivs_shc(ex,trK,Kx,Ky,Kz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,TZ,TZx,TZy,TZz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call rkhat_point(alpha,Sfx(i,j,k),Sfy(i,j,k),Sfz(i,j,k),chin1, &
                    Kx,Ky,Kz,TZx,TZy,TZz, &
                    gxx,gxy(i,j,k),gxz(i,j,k),gyy,gyz(i,j,k),gzz,kappa1,kappa2, &
                    trK(i,j,k),R(k),rtrK(i,j,k),TZ(i,j,k),x(i,j,k),y(i,j,k),z(i,j,k))
      enddo
     enddo
    enddo
   endif

  return

  end subroutine cpbcrtrK
!construct rTZ rhs
  subroutine cpbcrtheta(ex,crho,sigma,R,x,y,z,                                  &
               drhodx, drhody, drhodz,                                         &
               dsigmadx,dsigmady,dsigmadz,                                     &
               dRdx,dRdy,dRdz,                                                 &
               drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                &
               dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,    &
               dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,                            &
               xmin,ymin,zmin,xmax,ymax,zmax,rTheta,&
               chi,dxx,gxy,gxz,dyy,gyz,dzz, &
               Lap,Sfx,Sfy,Sfz,TZ,Symmetry,sst,kappa1,kappa2)

  implicit none
 
!~~~~~~> Input parameters:
  integer, intent(in):: ex(1:3),Symmetry,sst
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
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: chi,dxx,gxy,gxz,dyy,gyz,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Lap,Sfx,Sfy,Sfz,TZ
  real*8,  intent(in):: xmin,ymin,zmin,xmax,ymax,zmax,kappa1,kappa2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout)::rTheta
!~~~~~~> Other variables:
  real*8 :: alpha,chin1,gxx,gyy,gzz
  real*8 :: TZx,TZy,TZz
  logical :: gont
  integer :: i, j, k
  integer :: layer(1:6,1:6),gp
! index of layer, first one: i,j,k; second one: front back etc. boundary
  integer :: kmin,kmax
  real*8, parameter :: ZEO = 0.d0, ONE = 1.d0, TWO=2.d0
  real*8 :: dR
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0

  dR = R(2)-R(1)
  
  kmax = ex(3)

  kmin = 1

layer(1:3,:) = 1
layer(4:6,:) =-1
 
if(dabs(R(ex(3))-zmax) < dR)then
  layer(1,3) = 1
  layer(2,3) = 1
  layer(4,3) = ex(1)
  layer(5,3) = ex(2)
! consider buffer points near boundary 
  layer(3,3) = ex(3) - CPBC_ghost_width
  layer(6,3) = ex(3) - CPBC_ghost_width
endif

  gp = 3

  gont = any( layer(:,gp) == - 1 )
 
   if( .not. gont ) then

    do k = layer(3,gp), layer(6,gp)
     do j = layer(2,gp), layer(5,gp)
      do i = layer(1,gp), layer(4,gp)
         alpha = Lap(i,j,k)+1.d0
         chin1 = chi(i,j,k)+1.d0
         gxx = dxx(i,j,k)+1.d0
         gyy = dyy(i,j,k)+1.d0
         gzz = dzz(i,j,k)+1.d0
         call point_fderivs_shc(ex,TZ,TZx,TZy,TZz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call rtheta_point(alpha,Sfx(i,j,k),Sfy(i,j,k),Sfz(i,j,k),chin1, &
                    TZx,TZy,TZz, &
                    gxx,gxy(i,j,k),gxz(i,j,k),gyy,gyz(i,j,k),gzz,kappa1,kappa2, &
                    R(k),rTheta(i,j,k),TZ(i,j,k),x(i,j,k),y(i,j,k),z(i,j,k))
      enddo
     enddo
    enddo
   endif

  return

  end subroutine cpbcrtheta
!construct rGam rhs
  subroutine cpbcrgam(ex,crho,sigma,R,x,y,z,                                  &
               drhodx, drhody, drhodz,                                         &
               dsigmadx,dsigmady,dsigmadz,                                     &
               dRdx,dRdy,dRdz,                                                 &
               drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                &
               dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,    &
               dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,                            &
               xmin,ymin,zmin,xmax,ymax,zmax,rGamAx,rGamAy,rGamAz,rGams,&
               chi,trK,dxx,gxy,gxz,dyy,gyz,dzz, &
               Lap,Sfx,Sfy,Sfz,TZ,Gamx,Gamy,Gamz,Symmetry,sst,eta)

  implicit none
 
!~~~~~~> Input parameters:
  integer, intent(in):: ex(1:3),Symmetry,sst
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
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: trK,chi,dxx,gxy,gxz,dyy,gyz,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Lap,Sfx,Sfy,Sfz,TZ,Gamx,Gamy,Gamz
  real*8,  intent(in):: xmin,ymin,zmin,xmax,ymax,zmax,eta
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout)::rGamAx,rGamAy,rGamAz,rGams
!~~~~~~> Other variables:
  real*8 :: alpha,chin1,gxx,gyy,gzz
  real*8 :: sfxx,sfyx,sfzx
  real*8 :: sfxy,sfyy,sfzy
  real*8 :: sfxz,sfyz,sfzz
  real*8 :: sfxxx,sfyxx,sfzxx
  real*8 :: sfxxy,sfyxy,sfzxy
  real*8 :: sfxxz,sfyxz,sfzxz
  real*8 :: sfxyy,sfyyy,sfzyy
  real*8 :: sfxyz,sfyyz,sfzyz
  real*8 :: sfxzz,sfyzz,sfzzz
  real*8 :: Gamxx,Gamyx,Gamzx
  real*8 :: Gamxy,Gamyy,Gamzy
  real*8 :: Gamxz,Gamyz,Gamzz
  real*8 :: Kx,Ky,Kz,TZx,TZy,TZz
  logical :: gont
  integer :: i, j, k
  integer :: layer(1:6,1:6),gp
! index of layer, first one: i,j,k; second one: front back etc. boundary
  integer :: kmin,kmax
  real*8, parameter :: ZEO = 0.d0, ONE = 1.d0, TWO=2.d0
  real*8 :: dR
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0

  dR = R(2)-R(1)
  
  kmax = ex(3)

  kmin = 1

layer(1:3,:) = 1
layer(4:6,:) =-1
 
if(dabs(R(ex(3))-zmax) < dR)then
  layer(1,3) = 1
  layer(2,3) = 1
  layer(4,3) = ex(1)
  layer(5,3) = ex(2)
! consider buffer points near boundary 
  layer(3,3) = ex(3) - CPBC_ghost_width
  layer(6,3) = ex(3) - CPBC_ghost_width
endif

  gp = 3

  gont = any( layer(:,gp) == - 1 )
 
   if( .not. gont ) then

    do k = layer(3,gp), layer(6,gp)
     do j = layer(2,gp), layer(5,gp)
      do i = layer(1,gp), layer(4,gp)
         alpha = Lap(i,j,k)+1.d0
         chin1 = chi(i,j,k)+1.d0
         gxx = dxx(i,j,k)+1.d0
         gyy = dyy(i,j,k)+1.d0
         gzz = dzz(i,j,k)+1.d0
         call point_fderivs_shc(ex,trK,Kx,Ky,Kz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,TZ,TZx,TZy,TZz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Gamx,Gamxx,Gamxy,Gamxz,crho,sigma,R,ANTI,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Gamy,Gamyx,Gamyy,Gamyz,crho,sigma,R,SYM,ANTI,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Gamz,Gamzx,Gamzy,Gamzz,crho,sigma,R,SYM,SYM,ANTI,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Sfx,sfxx,sfxy,sfxz,crho,sigma,R,ANTI,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Sfy,sfyx,sfyy,sfyz,crho,sigma,R,SYM,ANTI,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Sfz,sfzx,sfzy,sfzz,crho,sigma,R,SYM,SYM,ANTI,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fdderivs_shc(ex,Sfx,sfxxx,sfxxy,sfxxz,sfxyy,sfxyz,sfxzz,crho,sigma,R,ANTI,SYM ,SYM ,Symmetry,0,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)
         call point_fdderivs_shc(ex,Sfy,sfyxx,sfyxy,sfyxz,sfyyy,sfyyz,sfyzz,crho,sigma,R,SYM ,ANTI,SYM ,Symmetry,0,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)
         call point_fdderivs_shc(ex,Sfz,sfzxx,sfzxy,sfzxz,sfzyy,sfzyz,sfzzz,crho,sigma,R,SYM ,SYM ,ANTI,Symmetry,0,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)
         call rgam_point(alpha,Sfx(i,j,k),Sfy(i,j,k),Sfz(i,j,k),chin1, &
                    sfxx,sfyx,sfzx, &
                    sfxy,sfyy,sfzy, &
                    sfxz,sfyz,sfzz, &
                    sfxxx,sfyxx,sfzxx, &
                    sfxxy,sfyxy,sfzxy, &
                    sfxxz,sfyxz,sfzxz, &
                    sfxyy,sfyyy,sfzyy, &
                    sfxyz,sfyyz,sfzyz, &
                    sfxzz,sfyzz,sfzzz, &
                    Gamxx,Gamyx,Gamzx, &
                    Gamxy,Gamyy,Gamzy, &
                    Gamxz,Gamyz,Gamzz, &
                    Kx,Ky,Kz,TZx,TZy,TZz, &
                    gxx,gxy(i,j,k),gxz(i,j,k),gyy,gyz(i,j,k),gzz,&
                    R(k),rGamAx(i,j,k),rGamAy(i,j,k),rGamAz(i,j,k),rGams(i,j,k), &
                    eta,x(i,j,k),y(i,j,k),z(i,j,k))
      enddo
     enddo
    enddo
   endif

  return

  end subroutine cpbcrgam
!construct rA rhs
  subroutine cpbcra(ex,crho,sigma,R,x,y,z,                                     &
               drhodx, drhody, drhodz,                                         &
               dsigmadx,dsigmady,dsigmadz,                                     &
               dRdx,dRdy,dRdz,                                                 &
               drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                &
               dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,    &
               dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,                            &
               xmin,ymin,zmin,xmax,ymax,zmax, &
               rACABTFxx,rACABTFxy,rACABTFxz,rACABTFyy,rACABTFyz,rACABTFzz,&
               rACsAx,rACsAy,rACsAz,rACss, &
               chi,trK,dxx,gxy,gxz,dyy,gyz,dzz, &
               Axx,Axy,Axz,Ayy,Ayz,Azz, &
               Lap,Sfx,Sfy,Sfz,TZ,Gamx,Gamy,Gamz,Symmetry,sst,kappa1)

  implicit none
 
!~~~~~~> Input parameters:
  integer, intent(in):: ex(1:3),Symmetry,sst
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
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: trK,chi,dxx,gxy,gxz,dyy,gyz,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Lap,Sfx,Sfy,Sfz,TZ,Gamx,Gamy,Gamz,Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8,  intent(in):: xmin,ymin,zmin,xmax,ymax,zmax,kappa1
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout)::rACABTFxx,rACABTFxy,rACABTFxz,rACABTFyy,rACABTFyz,rACABTFzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout)::rACsAx,rACsAy,rACsAz,rACss
!~~~~~~> Other variables:
  real*8 :: alpha,chin1,gxx,gyy,gzz
  real*8 :: sfxx,sfyx,sfzx
  real*8 :: sfxy,sfyy,sfzy
  real*8 :: sfxz,sfyz,sfzz
  real*8 :: sfxxx,sfyxx,sfzxx
  real*8 :: sfxxy,sfyxy,sfzxy
  real*8 :: sfxxz,sfyxz,sfzxz
  real*8 :: sfxyy,sfyyy,sfzyy
  real*8 :: sfxyz,sfyyz,sfzyz
  real*8 :: sfxzz,sfyzz,sfzzz
  real*8 :: Gamxx,Gamyx,Gamzx
  real*8 :: Gamxy,Gamyy,Gamzy
  real*8 :: Gamxz,Gamyz,Gamzz
  real*8 :: Kx,Ky,Kz,TZx,TZy,TZz
  real*8 :: Lapx,Lapy,Lapz,Lapxx,Lapxy,Lapxz,Lapyy,Lapyz,Lapzz
  real*8 :: chix,chiy,chiz
  real*8 :: chixx,chixy,chixz,chiyy,chiyz,chizz
  real*8 :: Axxx,Axxy,Axxz
  real*8 :: Axyx,Axyy,Axyz
  real*8 :: Axzx,Axzy,Axzz
  real*8 :: Ayyx,Ayyy,Ayyz
  real*8 :: Ayzx,Ayzy,Ayzz
  real*8 :: Azzx,Azzy,Azzz
  real*8 :: gxxx,gxxy,gxxz
  real*8 :: gxyx,gxyy,gxyz
  real*8 :: gxzx,gxzy,gxzz
  real*8 :: gyyx,gyyy,gyyz
  real*8 :: gyzx,gyzy,gyzz
  real*8 :: gzzx,gzzy,gzzz
  real*8 :: gxxxx,gxxxy,gxxxz,gxxyy,gxxyz,gxxzz
  real*8 :: gxyxx,gxyxy,gxyxz,gxyyy,gxyyz,gxyzz
  real*8 :: gxzxx,gxzxy,gxzxz,gxzyy,gxzyz,gxzzz
  real*8 :: gyyxx,gyyxy,gyyxz,gyyyy,gyyyz,gyyzz
  real*8 :: gyzxx,gyzxy,gyzxz,gyzyy,gyzyz,gyzzz
  real*8 :: gzzxx,gzzxy,gzzxz,gzzyy,gzzyz,gzzzz
  logical :: gont
  integer :: i, j, k
  integer :: layer(1:6,1:6),gp
! index of layer, first one: i,j,k; second one: front back etc. boundary
  integer :: kmin,kmax
  real*8, parameter :: ZEO = 0.d0, ONE = 1.d0, TWO=2.d0
  real*8 :: dR
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0

  dR = R(2)-R(1)
  
  kmax = ex(3)

  kmin = 1

layer(1:3,:) = 1
layer(4:6,:) =-1
 
if(dabs(R(ex(3))-zmax) < dR)then
  layer(1,3) = 1
  layer(2,3) = 1
  layer(4,3) = ex(1)
  layer(5,3) = ex(2)
! consider buffer points near boundary 
  layer(3,3) = ex(3) - CPBC_ghost_width
  layer(6,3) = ex(3) - CPBC_ghost_width
endif

  gp = 3

  gont = any( layer(:,gp) == - 1 )
 
   if( .not. gont ) then

    do k = layer(3,gp), layer(6,gp)
     do j = layer(2,gp), layer(5,gp)
      do i = layer(1,gp), layer(4,gp)
         alpha = Lap(i,j,k)+1.d0
         chin1 = chi(i,j,k)+1.d0
         gxx = dxx(i,j,k)+1.d0
         gyy = dyy(i,j,k)+1.d0
         gzz = dzz(i,j,k)+1.d0
         call point_fderivs_shc(ex,trK,Kx,Ky,Kz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,TZ,TZx,TZy,TZz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Axx,Axxx,Axxy,Axxz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Axy,Axyx,Axyy,Axyz,crho,sigma,R,ANTI,ANTI,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Axz,Axzx,Axzy,Axzz,crho,sigma,R,ANTI,SYM,ANTI,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Ayy,Ayyx,Ayyy,Ayyz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Ayz,Ayzx,Ayzy,Ayzz,crho,sigma,R,SYM,ANTI,ANTI,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Azz,Azzx,Azzy,Azzz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,dxx,gxxx,gxxy,gxxz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,gxy,gxyx,gxyy,gxyz,crho,sigma,R,ANTI,ANTI,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,gxz,gxzx,gxzy,gxzz,crho,sigma,R,ANTI,SYM,ANTI,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,dyy,gyyx,gyyy,gyyz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,gyz,gyzx,gyzy,gyzz,crho,sigma,R,SYM,ANTI,ANTI,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,dzz,gzzx,gzzy,gzzz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Gamx,Gamxx,Gamxy,Gamxz,crho,sigma,R,ANTI,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Gamy,Gamyx,Gamyy,Gamyz,crho,sigma,R,SYM,ANTI,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Gamz,Gamzx,Gamzy,Gamzz,crho,sigma,R,SYM,SYM,ANTI,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Sfx,sfxx,sfxy,sfxz,crho,sigma,R,ANTI,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Sfy,sfyx,sfyy,sfyz,crho,sigma,R,SYM,ANTI,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Sfz,sfzx,sfzy,sfzz,crho,sigma,R,SYM,SYM,ANTI,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fdderivs_shc(ex,Sfx,sfxxx,sfxxy,sfxxz,sfxyy,sfxyz,sfxzz,crho,sigma,R,ANTI,SYM ,SYM ,Symmetry,0,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)
         call point_fdderivs_shc(ex,Sfy,sfyxx,sfyxy,sfyxz,sfyyy,sfyyz,sfyzz,crho,sigma,R,SYM ,ANTI,SYM ,Symmetry,0,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)
         call point_fdderivs_shc(ex,Sfz,sfzxx,sfzxy,sfzxz,sfzyy,sfzyz,sfzzz,crho,sigma,R,SYM ,SYM ,ANTI,Symmetry,0,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)

         call point_fderivs_shc(ex,chi,chix,chiy,chiz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fderivs_shc(ex,Lap,Lapx,Lapy,Lapz,crho,sigma,R,SYM,SYM,SYM,Symmetry,0,sst, &
                                drhodx, drhody, drhodz,                                  &
                                dsigmadx,dsigmady,dsigmadz,                              &
                                dRdx,dRdy,dRdz,i,j,k)
         call point_fdderivs_shc(ex,Lap,Lapxx,Lapxy,Lapxz,Lapyy,Lapyz,Lapzz,crho,sigma,R,SYM ,SYM ,SYM ,Symmetry,0,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)
         call point_fdderivs_shc(ex,chi,chixx,chixy,chixz,chiyy,chiyz,chizz,crho,sigma,R,SYM ,SYM ,SYM ,Symmetry,0,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)
         call point_fdderivs_shc(ex,dxx,gxxxx,gxxxy,gxxxz,gxxyy,gxxyz,gxxzz,crho,sigma,R,SYM ,SYM ,SYM ,Symmetry,0,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)
         call point_fdderivs_shc(ex,dyy,gyyxx,gyyxy,gyyxz,gyyyy,gyyyz,gyyzz,crho,sigma,R,SYM ,SYM ,SYM ,Symmetry,0,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)
         call point_fdderivs_shc(ex,dzz,gzzxx,gzzxy,gzzxz,gzzyy,gzzyz,gzzzz,crho,sigma,R,SYM ,SYM ,SYM ,Symmetry,0,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)
         call point_fdderivs_shc(ex,gxy,gxyxx,gxyxy,gxyxz,gxyyy,gxyyz,gxyzz,crho,sigma,R,ANTI,ANTI,SYM ,Symmetry,0,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)
         call point_fdderivs_shc(ex,gxz,gxzxx,gxzxy,gxzxz,gxzyy,gxzyz,gxzzz,crho,sigma,R,ANTI,SYM ,ANTI,Symmetry,0,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)
         call point_fdderivs_shc(ex,gyz,gyzxx,gyzxy,gyzxz,gyzyy,gyzyz,gyzzz,crho,sigma,R,SYM ,ANTI,ANTI ,Symmetry,0,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)
         call ra_point(Axx(i,j,k),Axy(i,j,k),Axz(i,j,k),Ayy(i,j,k),Ayz(i,j,k),Azz(i,j,k), &
                    alpha,Sfx(i,j,k),Sfy(i,j,k),Sfz(i,j,k),chin1, &
                    Lapx,Axxx,Axyx,Axzx,Ayyx,Ayzx,Azzx, &
                    Lapy,Axxy,Axyy,Axzy,Ayyy,Ayzy,Azzy, &
                    Lapz,Axxz,Axyz,Axzz,Ayyz,Ayzz,Azzz, &
                    sfxx,sfyx,sfzx, &
                    sfxy,sfyy,sfzy, &
                    sfxz,sfyz,sfzz, &
                    chix,chiy,chiz, &
                    Lapxx,Lapxy,Lapxz,Lapyy,Lapyz,Lapzz, &
                    sfxxx,sfyxx,sfzxx, &
                    sfxxy,sfyxy,sfzxy, &
                    sfxxz,sfyxz,sfzxz, &
                    sfxyy,sfyyy,sfzyy, &
                    sfxyz,sfyyz,sfzyz, &
                    sfxzz,sfyzz,sfzzz, &
                    chixx,chixy,chixz,chiyy,chiyz,chizz, &
                    gxxxx,gxyxx,gxzxx,gyyxx,gyzxx,gzzxx, &
                    gxxxy,gxyxy,gxzxy,gyyxy,gyzxy,gzzxy, &
                    gxxxz,gxyxz,gxzxz,gyyxz,gyzxz,gzzxz, &
                    gxxyy,gxyyy,gxzyy,gyyyy,gyzyy,gzzyy, &
                    gxxyz,gxyyz,gxzyz,gyyyz,gyzyz,gzzyz, &
                    gxxzz,gxyzz,gxzzz,gyyzz,gyzzz,gzzzz, &
                    Gamxx,gxxx,gxyx,gxzx, &
                    Gamyx,gyyx,gyzx, &
                    Gamzx,gzzx, &
                    Gamxy,gxxy,gxyy,gxzy, &
                    Gamyy,gyyy,gyzy, &
                    Gamzy,gzzy, &
                    Gamxz,gxxz,gxyz,gxzz, &
                    Gamyz,gyyz,gyzz, &
                    Gamzz,gzzz, &
                    Kx,Ky,Kz,TZx,TZy,TZz, &
                    Gamx(i,j,k),gxx,gxy(i,j,k),gxz(i,j,k), &
                    Gamy(i,j,k),gyy,gyz(i,j,k), &
                    Gamz(i,j,k),gzz, &
                    kappa1,trK(i,j,k), &
                    R(k),rACABTFxx(i,j,k),rACABTFxy(i,j,k),rACABTFxz(i,j,k), &
                    rACABTFyy(i,j,k),rACABTFyz(i,j,k),rACABTFzz(i,j,k), &
                    rACsAx(i,j,k),rACsAy(i,j,k),rACsAz(i,j,k), &
                    rACss(i,j,k),TZ(i,j,k), &
                    x(i,j,k),y(i,j,k),z(i,j,k))
      enddo
     enddo
    enddo
   endif

  return

  end subroutine cpbcra
