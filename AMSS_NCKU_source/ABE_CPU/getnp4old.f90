

#include "macrodef.fh"

!-----------------------------------------------------------------------------
!
! compute rhw Newman-Penrose Weyl scalar Psi4
! for BSSN dynamical variables
!
!-----------------------------------------------------------------------------

  subroutine getnp4old(ex, X, Y, Z,  chi,  trK,             &
                      dxx,  gxy,  gxz,  dyy,  gyz,  dzz, &
                      Axx,  Axy,  Axz,  Ayy,  Ayz,  Azz, &
                      Gmx,  Gmy,  Gmz,                   & 
                      Lap,  Sfx,  Sfy,  Sfz,Rpsi4,Ipsi4, symmetry)

  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3),symmetry
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: chi,trK
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: dxx,gxy,gxz,dyy,gyz,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gmx,Gmy,Gmz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Lap,Sfx,Sfy,Sfz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Rpsi4,Ipsi4

!~~~~~~> Other variables:

  real*8, dimension(ex(1),ex(2),ex(3)) :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupyy,gupyz,gupzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: ep4phi,alpn1 
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gyy,gzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: phi,phix,phiy,phiz
  real*8, dimension(ex(1),ex(2),ex(3)) :: phixx,phixy,phixz,phiyy,phiyz,phizz
  real*8, dimension(ex(1),ex(2),ex(3)) :: tRxyxy, tRxyxz, tRxyyz, tRxzxz, tRxzyz, tRyzyz
  real*8, dimension(ex(1),ex(2),ex(3)) ::  Rxyxy,  Rxyxz,  Rxyyz,  Rxzxz,  Rxzyz,  Ryzyz
  real*8, dimension(ex(1),ex(2),ex(3)) :: vx,vy,vz,ux,uy,uz,wx,wy,wz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamzyy, Gamzyz, Gamzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: tRxx,tRxy,tRxz,tRyy,tRyz,tRzz
  real*8, dimension(ex(1),ex(2),ex(3)) ::  Rxx, Rxy, Rxz, Ryy, Ryz, Rzz

  real*8, dimension(ex(1),ex(2),ex(3)) :: Kxx,Kxy,Kxz,Kyy,Kyz,Kzz
!D_i K_jk ---> DKijk
  real*8, dimension(ex(1),ex(2),ex(3)) :: DKxxx,DKxxy,DKxxz,DKxyy,DKxyz,DKxzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: DKyxx,DKyxy,DKyxz,DKyyy,DKyyz,DKyzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: DKzxx,DKzxy,DKzxz,DKzyy,DKzyz,DKzzz
! Aij,k --> stored as Aijk
  real*8, dimension(ex(1),ex(2),ex(3))::Axxx,Axxy,Axxz
  real*8, dimension(ex(1),ex(2),ex(3))::Axyx,Axyy,Axyz
  real*8, dimension(ex(1),ex(2),ex(3))::Axzx,Axzy,Axzz
  real*8, dimension(ex(1),ex(2),ex(3))::Ayyx,Ayyy,Ayyz
  real*8, dimension(ex(1),ex(2),ex(3))::Ayzx,Ayzy,Ayzz
  real*8, dimension(ex(1),ex(2),ex(3))::Azzx,Azzy,Azzz
! trK,i
  real*8, dimension(ex(1),ex(2),ex(3))::Kx,Ky,Kz

  real*8, dimension(ex(1),ex(2),ex(3)) :: ass_Gamxxx,ass_Gamxxy,ass_Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: ass_Gamxyy,ass_Gamxyz,ass_Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: ass_Gamyxx,ass_Gamyxy,ass_Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: ass_Gamyyy,ass_Gamyyz,ass_Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: ass_Gamzxx,ass_Gamzxy,ass_Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: ass_Gamzyy,ass_Gamzyz,ass_Gamzzz
! first order partial derivative of metric
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxx,gxxy,gxxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxyx,gxyy,gxyz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxzx,gxzy,gxzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyyx,gyyy,gyyz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyzx,gyzy,gyzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gzzx,gzzy,gzzz
! second order partial derivative of metric
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxxx,gxxxy,gxxxz,gxxyy,gxxyz,gxxzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxyxx,gxyxy,gxyxz,gxyyy,gxyyz,gxyzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxzxx,gxzxy,gxzxz,gxzyy,gxzyz,gxzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyyxx,gyyxy,gyyxz,gyyyy,gyyyz,gyyzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyzxx,gyzxy,gyzxz,gyzyy,gyzyz,gyzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gzzxx,gzzxy,gzzxz,gzzyy,gzzyz,gzzzz

  real*8, parameter :: F1o4=2.5d-1,ONE=1.d0,TWO=2.d0,FOUR=4.d0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0

  phi = -0.25d0*dlog(chi+ONE)
!~~~~~~

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

  alpn1 = Lap + ONE

  ep4phi = dexp( FOUR * Phi )

!~~~~~~>

  call d1metric(ex,X,Y,Z,                   &
                dxx ,gxy ,gxz ,dyy ,gyz ,dzz , &
                gxxx,gxyx,gxzx,gyyx,gyzx,gzzx, &
                gxxy,gxyy,gxzy,gyyy,gyzy,gzzy, &
                gxxz,gxyz,gxzz,gyyz,gyzz,gzzz, symmetry)

  call d2metric(ex,X,Y,Z,                         &
                  dxx,  gxy,  gxz,  dyy,  gyz,  dzz, &
                gxxxx,gxxxy,gxxxz,gxxyy,gxxyz,gxxzz, &
                gxyxx,gxyxy,gxyxz,gxyyy,gxyyz,gxyzz, &
                gxzxx,gxzxy,gxzxz,gxzyy,gxzyz,gxzzz, &
                gyyxx,gyyxy,gyyxz,gyyyy,gyyyz,gyyzz, &
                gyzxx,gyzxy,gyzxz,gyzyy,gyzyz,gyzzz, &
                gzzxx,gzzxy,gzzxz,gzzyy,gzzyz,gzzzz, symmetry)

  call kind1_connection(ex, gxxx,gxyx,gxzx,gyyx,gyzx,gzzx, &
                            gxxy,gxyy,gxzy,gyyy,gyzy,gzzy, &
                            gxxz,gxyz,gxzz,gyyz,gyzz,gzzz, &
                        ass_Gamxxx, ass_Gamxxy, ass_Gamxxz,&
                        ass_Gamxyy, ass_Gamxyz, ass_Gamxzz,&
                        ass_Gamyxx, ass_Gamyxy, ass_Gamyxz,&
                        ass_Gamyyy, ass_Gamyyz, ass_Gamyzz,&
                        ass_Gamzxx, ass_Gamzxy, ass_Gamzxz,&
                        ass_Gamzyy, ass_Gamzyz, ass_Gamzzz)

  call kind2_connection(ex, gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
                        ass_Gamxxx, ass_Gamxxy, ass_Gamxxz,      &
                        ass_Gamxyy, ass_Gamxyz, ass_Gamxzz,      &
                        ass_Gamyxx, ass_Gamyxy, ass_Gamyxz,      & 
                        ass_Gamyyy, ass_Gamyyz, ass_Gamyzz,      &
                        ass_Gamzxx, ass_Gamzxy, ass_Gamzxz,      &
                        ass_Gamzyy, ass_Gamzyz, ass_Gamzzz,      &
                        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz, &
                        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz, &
                        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz)

!~~~~~~> derivs of conformal factor

  call fderivs(ex,phi,phix,phiy,phiz,X,Y,Z,SYM,SYM,SYM,Symmetry,0)

  call fdderivs(ex,phi,phixx,phixy,phixz,phiyy,phiyz,phizz,X,Y,Z, &
                 SYM,SYM,SYM,symmetry,0)

  call xcov_deriv(ex, phix, phiy, phiz,                           &
                   phixx,  phixy,  phixz,  phiyy,  phiyz,  phizz, &
                  Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz, &
                  Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz, &
                  Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz)

!~~~~~~> get spatial Riemann curvature

  call adm_riemann(ex, gxxxx,gxxxy,gxxxz,gxxyy,gxxyz,gxxzz, &
                       gxyxx,gxyxy,gxyxz,gxyyy,gxyyz,gxyzz, &
                       gxzxx,gxzxy,gxzxz,gxzyy,gxzyz,gxzzz, &
                       gyyxx,gyyxy,gyyxz,gyyyy,gyyyz,gyyzz, &
                       gyzxx,gyzxy,gyzxz,gyzyy,gyzyz,gyzzz, &
                       gzzxx,gzzxy,gzzxz,gzzyy,gzzyz,gzzzz, &
                       Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz, &
                       Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz, &
                       Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz, &
                       ass_Gamxxx,ass_Gamxxy,ass_Gamxxz, &
                       ass_Gamxyy,ass_Gamxyz,ass_Gamxzz, &
                       ass_Gamyxx,ass_Gamyxy,ass_Gamyxz, &
                       ass_Gamyyy,ass_Gamyyz,ass_Gamyzz, &
                       ass_Gamzxx,ass_Gamzxy,ass_Gamzxz, &
                       ass_Gamzyy,ass_Gamzyz,ass_Gamzzz, &
                       tRxyxy, tRxyxz, tRxyyz, tRxzxz, tRxzyz, tRyzyz)

  call get_physical_riemann(ex, ep4phi,                             &
                      dxx ,   gxy ,   gxz ,   dyy ,   gyz ,   dzz , &
                    gupxx , gupxy , gupxz , gupyy , gupyz , gupzz , &
                     phix ,  phiy ,  phiz ,                         &
                    phixx , phixy , phixz , phiyy , phiyz , phizz , &
                    tRxyxy, tRxyxz, tRxyyz, tRxzxz, tRxzyz, tRyzyz, &
                     Rxyxy,  Rxyxz,  Rxyyz,  Rxzxz,  Rxzyz,  Ryzyz)

!~~~~~~> get spatial Ricci tensor

   call adm_ricci(ex, gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, &
                     tRxyxy,tRxyxz,tRxyyz,tRxzxz,tRxzyz,tRyzyz, &
                       tRxx,  tRxy,  tRxz,  tRyy,  tRyz,  tRzz)

  call get_physical_ricci(ex,dxx,gxy,gxz,dyy,gyz,dzz,phix,phiy,phiz, &
                          phixx,phixy,phixz,phiyy,phiyz,phizz, &
                          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
                           tRxx, tRxy, tRxz, tRyy, tRyz, tRzz, &
                            Rxx,  Rxy,  Rxz,  Ryy,  Ryz,  Rzz)

!~~~~~~> get the real spatial extrinsic curvature

  call get_physical_k(ex, phi, trK, dxx, gxy, gxz, dyy, gyz, dzz, &
                                    Axx, Axy, Axz, Ayy, Ayz, Azz, &
                                    Kxx, Kxy, Kxz, Kyy, Kyz, Kzz)

!~~~~~~> derivs of trace of extrinsic curvature

  call fderivs(ex,trK, Kx, Ky, Kz,X,Y,Z,SYM,SYM,SYM,Symmetry,0)

!~~~~~~> derivs of tilde extrinsic curvature

  call fderivs(ex,Axx,Axxx,Axxy,Axxz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)
  call fderivs(ex,Axy,Axyx,Axyy,Axyz,X,Y,Z,ANTI,ANTI,SYM ,Symmetry,0)
  call fderivs(ex,Axz,Axzx,Axzy,Axzz,X,Y,Z,ANTI,SYM ,ANTI,Symmetry,0)
  call fderivs(ex,Ayy,Ayyx,Ayyy,Ayyz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)
  call fderivs(ex,Ayz,Ayzx,Ayzy,Ayzz,X,Y,Z,SYM ,ANTI,ANTI,Symmetry,0)
  call fderivs(ex,Azz,Azzx,Azzy,Azzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)

!~~~~~~> derivs of extrinsic curvature, Kij

  call get_diff_physical_k(ex, phi, trK, Kx, Ky, Kz, phix, phiy, phiz, &
                             dxx,  gxy,  gxz,  dyy,  gyz,  dzz, &
                           gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
                             Axx,  Axy,  Axz,  Ayy,  Ayz,  Azz, &
                            Axxx, Axxy, Axxz, Axyx, Axyy, Axyz, &
                            Axzx, Axzy, Axzz, Ayyx, Ayyy, Ayyz, &
                            Ayzx, Ayzy, Ayzz, Azzx, Azzy, Azzz, &
                            Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz, &
                            Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz, &
                            Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz, &
                               Kxx,    Kxy,    Kxz,    Kyy,    Kyz,    Kzz, &
                             DKxxx,  DKxxy,  DKxxz,  DKxyy,  DKxyz,  DKxzz, &
                             DKyxx,  DKyxy,  DKyxz,  DKyyy,  DKyyz,  DKyzz, &
                             DKzxx,  DKzxy,  DKzxz,  DKzyy,  DKzyz,  DKzzz)

!~~~~~~> get the Gram-Schmidt orthonormalize triad coordinate
#if (tetradtype == 0)
  call get_triad0(ex, X, Y, Z, ep4phi, gxx, gxy, gxz, gyy, gyz, gzz, &
                 vx,vy,vz,ux,uy,uz,wx,wy,wz)
#elif (tetradtype == 1)
  call get_triad1(ex, X, Y, Z, ep4phi, gxx, gxy, gxz, gyy, gyz, gzz, &
                 vx,vy,vz,ux,uy,uz,wx,wy,wz)
#elif (tetradtype == 2)  
  call get_triad2(ex, X, Y, Z, ep4phi, gxx, gxy, gxz, gyy, gyz, gzz, &
                 vx,vy,vz,ux,uy,uz,wx,wy,wz)
#endif

!~~~~~~> compute the Newnamm-Penrose psi4 which split real and image part

   ep4phi = ONE / ep4phi

   call bssn_compute_psi4(ex,ep4phi, alpn1,   Sfx,  Sfy,  Sfz, &
                          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
                          vx,vy,vz,ux,uy,uz,wx,wy,wz,          &
                          trK,Kxx,Kxy,Kxz,Kyy,Kyz,Kzz,         &
                          Rxyxy,Rxyxz,Rxyyz,Rxzxz,Rxzyz,Ryzyz, &
                          Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,        &
                          DKxxx,DKxxy,DKxxz,DKxyy,DKxyz,DKxzz, &
                          DKyxx,DKyxy,DKyxz,DKyyy,DKyyz,DKyzz, &
                          DKzxx,DKzxy,DKzxz,DKzyy,DKzyz,DKzzz, Rpsi4, Ipsi4)

  return

  end subroutine getnp4old
!-----------------------------------------------------------------------------------
! for shell
!

  subroutine getnp4old_ss(ex,crho,sigma,R, X, Y, Z,                               &
               drhodx, drhody, drhodz,                                         &
               dsigmadx,dsigmady,dsigmadz,                                     &
               dRdx,dRdy,dRdz,                                                 &
               drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                &
               dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,    &
               dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,                            &
                      chi,  trK,             &
                      dxx,  gxy,  gxz,  dyy,  gyz,  dzz, &
                      Axx,  Axy,  Axz,  Ayy,  Ayz,  Azz, &
                      Gmx,  Gmy,  Gmz,                   & 
                      Lap,  Sfx,  Sfy,  Sfz,Rpsi4,Ipsi4, symmetry,sst)

  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3),symmetry,sst
  double precision,intent(in),dimension(ex(1))::crho
  double precision,intent(in),dimension(ex(2))::sigma
  double precision,intent(in),dimension(ex(3))::R
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: X,Y,Z
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::drhodx, drhody, drhodz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::dsigmadx,dsigmady,dsigmadz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::dRdx,dRdy,dRdz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: chi,trK
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: dxx,gxy,gxz,dyy,gyz,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gmx,Gmy,Gmz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Lap,Sfx,Sfy,Sfz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Rpsi4,Ipsi4

!~~~~~~> Other variables:

  real*8, dimension(ex(1),ex(2),ex(3)) :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupyy,gupyz,gupzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: ep4phi,alpn1 
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gyy,gzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: phi,phix,phiy,phiz
  real*8, dimension(ex(1),ex(2),ex(3)) :: phixx,phixy,phixz,phiyy,phiyz,phizz
  real*8, dimension(ex(1),ex(2),ex(3)) :: tRxyxy, tRxyxz, tRxyyz, tRxzxz, tRxzyz, tRyzyz
  real*8, dimension(ex(1),ex(2),ex(3)) ::  Rxyxy,  Rxyxz,  Rxyyz,  Rxzxz,  Rxzyz,  Ryzyz
  real*8, dimension(ex(1),ex(2),ex(3)) :: vx,vy,vz,ux,uy,uz,wx,wy,wz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamzyy, Gamzyz, Gamzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: tRxx,tRxy,tRxz,tRyy,tRyz,tRzz
  real*8, dimension(ex(1),ex(2),ex(3)) ::  Rxx, Rxy, Rxz, Ryy, Ryz, Rzz

  real*8, dimension(ex(1),ex(2),ex(3)) :: Kxx,Kxy,Kxz,Kyy,Kyz,Kzz
!D_i K_jk ---> DKijk
  real*8, dimension(ex(1),ex(2),ex(3)) :: DKxxx,DKxxy,DKxxz,DKxyy,DKxyz,DKxzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: DKyxx,DKyxy,DKyxz,DKyyy,DKyyz,DKyzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: DKzxx,DKzxy,DKzxz,DKzyy,DKzyz,DKzzz
! Aij,k --> stored as Aijk
  real*8, dimension(ex(1),ex(2),ex(3))::Axxx,Axxy,Axxz
  real*8, dimension(ex(1),ex(2),ex(3))::Axyx,Axyy,Axyz
  real*8, dimension(ex(1),ex(2),ex(3))::Axzx,Axzy,Axzz
  real*8, dimension(ex(1),ex(2),ex(3))::Ayyx,Ayyy,Ayyz
  real*8, dimension(ex(1),ex(2),ex(3))::Ayzx,Ayzy,Ayzz
  real*8, dimension(ex(1),ex(2),ex(3))::Azzx,Azzy,Azzz
! trK,i
  real*8, dimension(ex(1),ex(2),ex(3))::Kx,Ky,Kz

  real*8, dimension(ex(1),ex(2),ex(3)) :: ass_Gamxxx,ass_Gamxxy,ass_Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: ass_Gamxyy,ass_Gamxyz,ass_Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: ass_Gamyxx,ass_Gamyxy,ass_Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: ass_Gamyyy,ass_Gamyyz,ass_Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: ass_Gamzxx,ass_Gamzxy,ass_Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: ass_Gamzyy,ass_Gamzyz,ass_Gamzzz
! first order partial derivative of metric
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxx,gxxy,gxxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxyx,gxyy,gxyz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxzx,gxzy,gxzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyyx,gyyy,gyyz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyzx,gyzy,gyzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gzzx,gzzy,gzzz
! second order partial derivative of metric
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxxx,gxxxy,gxxxz,gxxyy,gxxyz,gxxzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxyxx,gxyxy,gxyxz,gxyyy,gxyyz,gxyzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxzxx,gxzxy,gxzxz,gxzyy,gxzyz,gxzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyyxx,gyyxy,gyyxz,gyyyy,gyyyz,gyyzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyzxx,gyzxy,gyzxz,gyzyy,gyzyz,gyzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gzzxx,gzzxy,gzzxz,gzzyy,gzzyz,gzzzz

  real*8, parameter :: F1o4=2.5d-1,ONE=1.d0,TWO=2.d0,FOUR=4.d0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  real*8,parameter::TINYRR=1.d-14

  phi = -0.25d0*dlog(chi+ONE)
!~~~~~~

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

  alpn1 = Lap + ONE

  ep4phi = dexp( FOUR * Phi )

!~~~~~~>

  call fderivs_shc(ex,dxx,gxxx,gxxy,gxxz,crho,sigma,R, SYM, SYM,SYM,Symmetry,0,sst,          &
                       drhodx, drhody, drhodz,                                               &
                       dsigmadx,dsigmady,dsigmadz,                                           &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,gxy,gxyx,gxyy,gxyz,crho,sigma,R,ANTI,ANTI,SYM,Symmetry,0,sst,          &
                       drhodx, drhody, drhodz,                                               &
                       dsigmadx,dsigmady,dsigmadz,                                           &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,gxz,gxzx,gxzy,gxzz,crho,sigma,R,ANTI,SYM ,ANTI,Symmetry,0,sst,         &
                       drhodx, drhody, drhodz,                                               &
                       dsigmadx,dsigmady,dsigmadz,                                           &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,dyy,gyyx,gyyy,gyyz,crho,sigma,R, SYM, SYM,SYM,Symmetry,0,sst,          &
                       drhodx, drhody, drhodz,                                               &
                       dsigmadx,dsigmady,dsigmadz,                                           &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,gyz,gyzx,gyzy,gyzz,crho,sigma,R,SYM ,ANTI,ANTI,Symmetry,0,sst,         &
                       drhodx, drhody, drhodz,                                               &
                       dsigmadx,dsigmady,dsigmadz,                                           &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,dzz,gzzx,gzzy,gzzz,crho,sigma,R, SYM, SYM,SYM,Symmetry,0,sst,          &
                       drhodx, drhody, drhodz,                                               &
                       dsigmadx,dsigmady,dsigmadz,                                           &
                       dRdx,dRdy,dRdz)

  call fdderivs_shc(ex,dxx,gxxxx,gxxxy,gxxxz,gxxyy,gxxyz,gxxzz,crho,sigma,R, SYM, SYM,SYM ,Symmetry,0,sst,  &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz)
  call fdderivs_shc(ex,dyy,gyyxx,gyyxy,gyyxz,gyyyy,gyyyz,gyyzz,crho,sigma,R, SYM, SYM,SYM ,Symmetry,0,sst,  &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz)
  call fdderivs_shc(ex,dzz,gzzxx,gzzxy,gzzxz,gzzyy,gzzyz,gzzzz,crho,sigma,R, SYM, SYM,SYM ,Symmetry,0,sst,  &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz)
  call fdderivs_shc(ex,gxy,gxyxx,gxyxy,gxyxz,gxyyy,gxyyz,gxyzz,crho,sigma,R,ANTI,ANTI,SYM ,Symmetry,0,sst,  &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz)
  call fdderivs_shc(ex,gxz,gxzxx,gxzxy,gxzxz,gxzyy,gxzyz,gxzzz,crho,sigma,R,ANTI,SYM ,ANTI,Symmetry,0,sst,  &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz)
  call fdderivs_shc(ex,gyz,gyzxx,gyzxy,gyzxz,gyzyy,gyzyz,gyzzz,crho,sigma,R,SYM ,ANTI,ANTI,Symmetry,0,sst,  &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz)

  call kind1_connection(ex, gxxx,gxyx,gxzx,gyyx,gyzx,gzzx, &
                            gxxy,gxyy,gxzy,gyyy,gyzy,gzzy, &
                            gxxz,gxyz,gxzz,gyyz,gyzz,gzzz, &
                        ass_Gamxxx, ass_Gamxxy, ass_Gamxxz,&
                        ass_Gamxyy, ass_Gamxyz, ass_Gamxzz,&
                        ass_Gamyxx, ass_Gamyxy, ass_Gamyxz,&
                        ass_Gamyyy, ass_Gamyyz, ass_Gamyzz,&
                        ass_Gamzxx, ass_Gamzxy, ass_Gamzxz,&
                        ass_Gamzyy, ass_Gamzyz, ass_Gamzzz)

  call kind2_connection(ex, gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
                        ass_Gamxxx, ass_Gamxxy, ass_Gamxxz,      &
                        ass_Gamxyy, ass_Gamxyz, ass_Gamxzz,      &
                        ass_Gamyxx, ass_Gamyxy, ass_Gamyxz,      & 
                        ass_Gamyyy, ass_Gamyyz, ass_Gamyzz,      &
                        ass_Gamzxx, ass_Gamzxy, ass_Gamzxz,      &
                        ass_Gamzyy, ass_Gamzyz, ass_Gamzzz,      &
                        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz, &
                        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz, &
                        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz)

!~~~~~~> derivs of conformal factor
  call fderivs_shc(ex,phi,phix,phiy,phiz,crho,sigma,R, SYM, SYM,SYM,Symmetry,0,sst,          &
                       drhodx, drhody, drhodz,                                               &
                       dsigmadx,dsigmady,dsigmadz,                                           &
                       dRdx,dRdy,dRdz)

  call fdderivs_shc(ex,phi,phixx,phixy,phixz,phiyy,phiyz,phizz,crho,sigma,R,SYM ,SYM ,SYM ,Symmetry,0,sst,  &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz)

  call xcov_deriv(ex, phix, phiy, phiz,                           &
                   phixx,  phixy,  phixz,  phiyy,  phiyz,  phizz, &
                  Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz, &
                  Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz, &
                  Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz)

!~~~~~~> get spatial Riemann curvature

  call adm_riemann(ex, gxxxx,gxxxy,gxxxz,gxxyy,gxxyz,gxxzz, &
                       gxyxx,gxyxy,gxyxz,gxyyy,gxyyz,gxyzz, &
                       gxzxx,gxzxy,gxzxz,gxzyy,gxzyz,gxzzz, &
                       gyyxx,gyyxy,gyyxz,gyyyy,gyyyz,gyyzz, &
                       gyzxx,gyzxy,gyzxz,gyzyy,gyzyz,gyzzz, &
                       gzzxx,gzzxy,gzzxz,gzzyy,gzzyz,gzzzz, &
                       Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz, &
                       Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz, &
                       Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz, &
                       ass_Gamxxx,ass_Gamxxy,ass_Gamxxz, &
                       ass_Gamxyy,ass_Gamxyz,ass_Gamxzz, &
                       ass_Gamyxx,ass_Gamyxy,ass_Gamyxz, &
                       ass_Gamyyy,ass_Gamyyz,ass_Gamyzz, &
                       ass_Gamzxx,ass_Gamzxy,ass_Gamzxz, &
                       ass_Gamzyy,ass_Gamzyz,ass_Gamzzz, &
                       tRxyxy, tRxyxz, tRxyyz, tRxzxz, tRxzyz, tRyzyz)

  call get_physical_riemann(ex, ep4phi,                             &
                      dxx ,   gxy ,   gxz ,   dyy ,   gyz ,   dzz , &
                    gupxx , gupxy , gupxz , gupyy , gupyz , gupzz , &
                     phix ,  phiy ,  phiz ,                         &
                    phixx , phixy , phixz , phiyy , phiyz , phizz , &
                    tRxyxy, tRxyxz, tRxyyz, tRxzxz, tRxzyz, tRyzyz, &
                     Rxyxy,  Rxyxz,  Rxyyz,  Rxzxz,  Rxzyz,  Ryzyz)

!~~~~~~> get spatial Ricci tensor

   call adm_ricci(ex, gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, &
                     tRxyxy,tRxyxz,tRxyyz,tRxzxz,tRxzyz,tRyzyz, &
                       tRxx,  tRxy,  tRxz,  tRyy,  tRyz,  tRzz)

  call get_physical_ricci(ex,dxx,gxy,gxz,dyy,gyz,dzz,phix,phiy,phiz, &
                          phixx,phixy,phixz,phiyy,phiyz,phizz, &
                          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
                           tRxx, tRxy, tRxz, tRyy, tRyz, tRzz, &
                            Rxx,  Rxy,  Rxz,  Ryy,  Ryz,  Rzz)

!~~~~~~> get the real spatial extrinsic curvature

  call get_physical_k(ex, phi, trK, dxx, gxy, gxz, dyy, gyz, dzz, &
                                    Axx, Axy, Axz, Ayy, Ayz, Azz, &
                                    Kxx, Kxy, Kxz, Kyy, Kyz, Kzz)

!~~~~~~> derivs of trace of extrinsic curvature
  call fderivs_shc(ex,trK,Kx,Ky,Kz,crho,sigma,R, SYM, SYM,SYM,Symmetry,0,sst,                &
                       drhodx, drhody, drhodz,                                               &
                       dsigmadx,dsigmady,dsigmadz,                                           &
                       dRdx,dRdy,dRdz)

!~~~~~~> derivs of tilde extrinsic curvature

  call fderivs_shc(ex,Axx,Axxx,Axxy,Axxz,crho,sigma,R, SYM, SYM,SYM,Symmetry,0,sst,         &
                       drhodx, drhody, drhodz,                                              &
                       dsigmadx,dsigmady,dsigmadz,                                          &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,Axy,Axyx,Axyy,Axyz,crho,sigma,R,ANTI,ANTI,SYM,Symmetry,0,sst,         &
                       drhodx, drhody, drhodz,                                              &
                       dsigmadx,dsigmady,dsigmadz,                                          &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,Axz,Axzx,Axzy,Axzz,crho,sigma,R,ANTI,SYM ,ANTI,Symmetry,0,sst,        &
                       drhodx, drhody, drhodz,                                              &
                       dsigmadx,dsigmady,dsigmadz,                                          &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,Ayy,Ayyx,Ayyy,Ayyz,crho,sigma,R, SYM, SYM,SYM,Symmetry,0,sst,         &
                       drhodx, drhody, drhodz,                                              &
                       dsigmadx,dsigmady,dsigmadz,                                          &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,Ayz,Ayzx,Ayzy,Ayzz,crho,sigma,R,SYM ,ANTI,ANTI,Symmetry,0,sst,        &
                       drhodx, drhody, drhodz,                                              &
                       dsigmadx,dsigmady,dsigmadz,                                          &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,Azz,Azzx,Azzy,Azzz,crho,sigma,R, SYM, SYM,SYM,Symmetry,0,sst,         &
                       drhodx, drhody, drhodz,                                              &
                       dsigmadx,dsigmady,dsigmadz,                                          &
                       dRdx,dRdy,dRdz)

!~~~~~~> derivs of extrinsic curvature, Kij

  call get_diff_physical_k(ex, phi, trK, Kx, Ky, Kz, phix, phiy, phiz, &
                             dxx,  gxy,  gxz,  dyy,  gyz,  dzz, &
                           gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
                             Axx,  Axy,  Axz,  Ayy,  Ayz,  Azz, &
                            Axxx, Axxy, Axxz, Axyx, Axyy, Axyz, &
                            Axzx, Axzy, Axzz, Ayyx, Ayyy, Ayyz, &
                            Ayzx, Ayzy, Ayzz, Azzx, Azzy, Azzz, &
                            Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz, &
                            Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz, &
                            Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz, &
                               Kxx,    Kxy,    Kxz,    Kyy,    Kyz,    Kzz, &
                             DKxxx,  DKxxy,  DKxxz,  DKxyy,  DKxyz,  DKxzz, &
                             DKyxx,  DKyxy,  DKyxz,  DKyyy,  DKyyz,  DKyzz, &
                             DKzxx,  DKzxy,  DKzxz,  DKzyy,  DKzyz,  DKzzz)

!~~~~~~> get the Gram-Schmidt orthonormalize triad coordinate

#if (tetradtype == 0)
  call get_triad0_ss(ex, X, Y, Z, ep4phi, gxx, gxy, gxz, gyy, gyz, gzz, &
                 vx,vy,vz,ux,uy,uz,wx,wy,wz)
#elif (tetradtype == 1)
  call get_triad1_ss(ex, X, Y, Z, ep4phi, gxx, gxy, gxz, gyy, gyz, gzz, &
                 vx,vy,vz,ux,uy,uz,wx,wy,wz)
#elif (tetradtype == 2)  
  call get_triad2_ss(ex, X, Y, Z, ep4phi, gxx, gxy, gxz, gyy, gyz, gzz, &
                 vx,vy,vz,ux,uy,uz,wx,wy,wz)
#endif

!~~~~~~> compute the Newnamm-Penrose psi4 which split real and image part

   ep4phi = ONE / ep4phi

   call bssn_compute_psi4(ex,ep4phi, alpn1,   Sfx,  Sfy,  Sfz, &
                          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
                          vx,vy,vz,ux,uy,uz,wx,wy,wz,          &
                          trK,Kxx,Kxy,Kxz,Kyy,Kyz,Kzz,         &
                          Rxyxy,Rxyxz,Rxyyz,Rxzxz,Rxzyz,Ryzyz, &
                          Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,        &
                          DKxxx,DKxxy,DKxxz,DKxyy,DKxyz,DKxzz, &
                          DKyxx,DKyxy,DKyxz,DKyyy,DKyyz,DKyzz, &
                          DKzxx,DKzxy,DKzxz,DKzyy,DKzyz,DKzzz, Rpsi4, Ipsi4)

  return

  end subroutine getnp4old_ss
!----------------------------------------------------------!
!                                                          !
!  derivatives related to 3-dimensional Riemann slice      !
!                                                          !
!----------------------------------------------------------!

!-----------------------------------------------------------------------------
!  Interface to compute the first order derivative of metric
!-----------------------------------------------------------------------------

  subroutine d1metric(ex,X,Y,Z,                   &
                      dxx ,gxy ,gxz ,dyy ,gyz ,dzz , &
                      gxxx,gxyx,gxzx,gyyx,gyzx,gzzx, &
                      gxxy,gxyy,gxzy,gyyy,gyzy,gzzy, &
                      gxxz,gxyz,gxzz,gyyz,gyzz,gzzz, symmetry)

  implicit none

!~~~~~~ Input parameters:

  integer,                              intent(in ) :: ex(1:3),symmetry
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: dxx,gxy,gxz,dyy,gyz,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gxxx,gxyx,gxzx,gyyx,gyzx,gzzx
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gxxy,gxyy,gxzy,gyyy,gyzy,gzzy
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gxxz,gxyz,gxzz,gyyz,gyzz,gzzz

!~~~~~~ local variables

  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0

!~~~~~~ 1st derivs of matric

   call fderivs(ex,dxx,gxxx,gxxy,gxxz,X,Y,Z,SYM ,SYM ,SYM ,symmetry,0)
   call fderivs(ex,gxy,gxyx,gxyy,gxyz,X,Y,Z,ANTI,ANTI,SYM ,symmetry,0)
   call fderivs(ex,gxz,gxzx,gxzy,gxzz,X,Y,Z,ANTI,SYM ,ANTI,symmetry,0)
   call fderivs(ex,dyy,gyyx,gyyy,gyyz,X,Y,Z,SYM ,SYM ,SYM ,symmetry,0)
   call fderivs(ex,gyz,gyzx,gyzy,gyzz,X,Y,Z,SYM ,ANTI,ANTI,symmetry,0)
   call fderivs(ex,dzz,gzzx,gzzy,gzzz,X,Y,Z,SYM ,SYM ,SYM ,symmetry,0)

   return

   end subroutine d1metric

!-----------------------------------------------------------------------------
!  Interface to compute the second order derivative of metric
!-----------------------------------------------------------------------------

  subroutine d2metric(ex,X,Y,Z,                         &
                        dxx,  gxy,  gxz,  dyy,  gyz,  dzz, &
                      gxxxx,gxxxy,gxxxz,gxxyy,gxxyz,gxxzz, &
                      gxyxx,gxyxy,gxyxz,gxyyy,gxyyz,gxyzz, &
                      gxzxx,gxzxy,gxzxz,gxzyy,gxzyz,gxzzz, &
                      gyyxx,gyyxy,gyyxz,gyyyy,gyyyz,gyyzz, &
                      gyzxx,gyzxy,gyzxz,gyzyy,gyzyz,gyzzz, &
                      gzzxx,gzzxy,gzzxz,gzzyy,gzzyz,gzzzz, symmetry)

  implicit none

!~~~~~~ Input parameters:

  integer,                              intent(in ) :: ex(1:3),symmetry
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: dxx,gxy,gxz,dyy,gyz,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gxxxx,gxxxy,gxxxz,gxxyy,gxxyz,gxxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gxyxx,gxyxy,gxyxz,gxyyy,gxyyz,gxyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gxzxx,gxzxy,gxzxz,gxzyy,gxzyz,gxzzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gyyxx,gyyxy,gyyxz,gyyyy,gyyyz,gyyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gyzxx,gyzxy,gyzxz,gyzyy,gyzyz,gyzzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gzzxx,gzzxy,gzzxz,gzzyy,gzzyz,gzzzz

!~~~~~~ local variables

  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0

!~~~~~~ 2nd derivs of matric

   call fdderivs(ex,dxx,gxxxx,gxxxy,gxxxz,gxxyy,gxxyz,gxxzz,X,Y,Z, &
                    SYM ,SYM ,SYM ,symmetry,0)
   call fdderivs(ex,gxy,gxyxx,gxyxy,gxyxz,gxyyy,gxyyz,gxyzz,X,Y,Z, &
                    ANTI,ANTI,SYM ,symmetry,0)
   call fdderivs(ex,gxz,gxzxx,gxzxy,gxzxz,gxzyy,gxzyz,gxzzz,X,Y,Z, &
                    ANTI,SYM ,ANTI,symmetry,0)
   call fdderivs(ex,dyy,gyyxx,gyyxy,gyyxz,gyyyy,gyyyz,gyyzz,X,Y,Z, &
                    SYM, SYM ,SYM ,symmetry,0)
   call fdderivs(ex,gyz,gyzxx,gyzxy,gyzxz,gyzyy,gyzyz,gyzzz,X,Y,Z, &
                    SYM ,ANTI,ANTI,symmetry,0)
   call fdderivs(ex,dzz,gzzxx,gzzxy,gzzxz,gzzyy,gzzyz,gzzzz,X,Y,Z, &
                    SYM ,SYM ,SYM ,symmetry,0)

   return

   end subroutine d2metric
!----------------------------------------------------------!
!                                                          !
!  algebraic computation based on geometric quantites      !
!  and their partial derivatives related to 3-dimensional  !
!  Riemann slice                                           !
!                                                          !
!----------------------------------------------------------!

!-----------------------------------------------------------------------------
! Get first kind of connection coefficients
! based on first order derivative of metric
! ass_Gam_ijk = 1/2 *(g_ij,k + g_ki,j - g_jk,i)
!-----------------------------------------------------------------------------

  subroutine kind1_connection(ex,gxxx,gxyx,gxzx,gyyx,gyzx,gzzx, &
                                 gxxy,gxyy,gxzy,gyyy,gyzy,gzzy, &
                                 gxxz,gxyz,gxzz,gyyz,gyzz,gzzz, &
                            ass_Gamxxx, ass_Gamxxy, ass_Gamxxz, &
                            ass_Gamxyy, ass_Gamxyz, ass_Gamxzz, &
                            ass_Gamyxx, ass_Gamyxy, ass_Gamyxz, &
                            ass_Gamyyy, ass_Gamyyz, ass_Gamyzz, &
                            ass_Gamzxx, ass_Gamzxy, ass_Gamzxz, &
                            ass_Gamzyy, ass_Gamzyz, ass_Gamzzz)
  implicit none

!~~~~~~> Input parameters:
 
  integer, intent(in) :: ex(1:3)
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: gxxx,gxyx,gxzx
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: gyyx,gyzx,gzzx
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: gxxy,gxyy,gxzy
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: gyyy,gyzy,gzzy
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: gxxz,gxyz,gxzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: gyyz,gyzz,gzzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out) :: ass_Gamxxx,ass_Gamxxy,ass_Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out) :: ass_Gamxyy,ass_Gamxyz,ass_Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out) :: ass_Gamyxx,ass_Gamyxy,ass_Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out) :: ass_Gamyyy,ass_Gamyyz,ass_Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out) :: ass_Gamzxx,ass_Gamzxy,ass_Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out) :: ass_Gamzyy,ass_Gamzyz,ass_Gamzzz
 
!~~~~~~> Other variables:
 
  real*8, parameter :: HLF=0.5d0

!~~~~~~=  Get Connection coefficients
! ass_Gam_ijk = 1/2 *(g_ij,k + g_ki,j - g_jk,i)

  ass_Gamxxx = HLF * ( gxxx               )
  ass_Gamyxx = HLF * ( gxyx + gxyx - gxxy )
  ass_Gamzxx = HLF * ( gxzx + gxzx - gxxz )
  ass_Gamxyy = HLF * ( gxyy + gxyy - gyyx )
  ass_Gamyyy = HLF * ( gyyy               )
  ass_Gamzyy = HLF * ( gyzy + gyzy - gyyz )
  ass_Gamxzz = HLF * ( gxzz + gxzz - gzzx )
  ass_Gamyzz = HLF * ( gyzz + gyzz - gzzy )
  ass_Gamzzz = HLF * ( gzzz               )
  ass_Gamxxy = HLF * ( gxxy + gxyx - gxyx )
  ass_Gamyxy = HLF * ( gxyy + gyyx - gxyy )
  ass_Gamzxy = HLF * ( gxzy + gyzx - gxyz )
  ass_Gamxxz = HLF * ( gxxz + gxzx - gxzx )
  ass_Gamyxz = HLF * ( gxyz + gyzx - gxzy )
  ass_Gamzxz = HLF * ( gxzz + gzzx - gxzz )
  ass_Gamxyz = HLF * ( gxyz + gxzy - gyzx )
  ass_Gamyyz = HLF * ( gyyz + gyzy - gyzy )
  ass_Gamzyz = HLF * ( gyzz + gzzy - gyzz )

  return

  end subroutine kind1_connection

!-----------------------------------------------------------------------------
! Get second kind of connection coefficients
! based on first kind of connection coefficients
! and gup
!-----------------------------------------------------------------------------

  subroutine kind2_connection(ex,gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
                                  ass_Gamxxx, ass_Gamxxy, ass_Gamxxz, &
                                  ass_Gamxyy, ass_Gamxyz, ass_Gamxzz, &
                                  ass_Gamyxx, ass_Gamyxy, ass_Gamyxz, &
                                  ass_Gamyyy, ass_Gamyyz, ass_Gamyzz, &
                                  ass_Gamzxx, ass_Gamzxy, ass_Gamzxz, &
                                  ass_Gamzyy, ass_Gamzyz, ass_Gamzzz, &
                      Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz, &
                      Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz, &
                      Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz)
  implicit none

!~~~~~~> Input parameters:
 
  integer, intent(in) :: ex(1:3)
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in) :: gupxx,gupxy,gupxz,gupyy,gupyz,gupzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in) :: ass_Gamxxx,ass_Gamxxy,ass_Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in) :: ass_Gamxyy,ass_Gamxyz,ass_Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in) :: ass_Gamyxx,ass_Gamyxy,ass_Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in) :: ass_Gamyyy,ass_Gamyyz,ass_Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in) :: ass_Gamzxx,ass_Gamzxy,ass_Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in) :: ass_Gamzyy,ass_Gamzyz,ass_Gamzzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out) :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out) :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out) :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out) :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out) :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out) :: Gamzyy, Gamzyz, Gamzzz
 
!~~~~~~> Other variables:

  Gamxxx = gupxx * ass_Gamxxx + gupxy * ass_Gamyxx + gupxz * ass_Gamzxx
  Gamxxy = gupxx * ass_Gamxxy + gupxy * ass_Gamyxy + gupxz * ass_Gamzxy
  Gamxxz = gupxx * ass_Gamxxz + gupxy * ass_Gamyxz + gupxz * ass_Gamzxz
  Gamxyy = gupxx * ass_Gamxyy + gupxy * ass_Gamyyy + gupxz * ass_Gamzyy
  Gamxyz = gupxx * ass_Gamxyz + gupxy * ass_Gamyyz + gupxz * ass_Gamzyz
  Gamxzz = gupxx * ass_Gamxzz + gupxy * ass_Gamyzz + gupxz * ass_Gamzzz

  Gamyxx = gupxy * ass_Gamxxx + gupyy * ass_Gamyxx + gupyz * ass_Gamzxx
  Gamyxy = gupxy * ass_Gamxxy + gupyy * ass_Gamyxy + gupyz * ass_Gamzxy
  Gamyxz = gupxy * ass_Gamxxz + gupyy * ass_Gamyxz + gupyz * ass_Gamzxz
  Gamyyy = gupxy * ass_Gamxyy + gupyy * ass_Gamyyy + gupyz * ass_Gamzyy
  Gamyyz = gupxy * ass_Gamxyz + gupyy * ass_Gamyyz + gupyz * ass_Gamzyz
  Gamyzz = gupxy * ass_Gamxzz + gupyy * ass_Gamyzz + gupyz * ass_Gamzzz

  Gamzxx = gupxz * ass_Gamxxx + gupyz * ass_Gamyxx + gupzz * ass_Gamzxx
  Gamzxy = gupxz * ass_Gamxxy + gupyz * ass_Gamyxy + gupzz * ass_Gamzxy
  Gamzxz = gupxz * ass_Gamxxz + gupyz * ass_Gamyxz + gupzz * ass_Gamzxz
  Gamzyy = gupxz * ass_Gamxyy + gupyz * ass_Gamyyy + gupzz * ass_Gamzyy
  Gamzyz = gupxz * ass_Gamxyz + gupyz * ass_Gamyyz + gupzz * ass_Gamzyz
  Gamzzz = gupxz * ass_Gamxzz + gupyz * ass_Gamyzz + gupzz * ass_Gamzzz

  return

  end subroutine kind2_connection

!----------------------------------------------------------------------
! compute Riemann tensor for three dimensional space
! based on second derivatives of metric 
!      and first knid and second kind of connection
!----------------------------------------------------------------------

  subroutine adm_riemann(ex,gxxxx,gxxxy,gxxxz,gxxyy,gxxyz,gxxzz, &
                            gxyxx,gxyxy,gxyxz,gxyyy,gxyyz,gxyzz, &
                            gxzxx,gxzxy,gxzxz,gxzyy,gxzyz,gxzzz, &
                            gyyxx,gyyxy,gyyxz,gyyyy,gyyyz,gyyzz, &
                            gyzxx,gyzxy,gyzxz,gyzyy,gyzyz,gyzzz, &
                            gzzxx,gzzxy,gzzxz,gzzyy,gzzyz,gzzzz, &
                 Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz, &
                 Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz, &
                 Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz, &
                               ass_Gamxxx,ass_Gamxxy,ass_Gamxxz, &
                               ass_Gamxyy,ass_Gamxyz,ass_Gamxzz, &
                               ass_Gamyxx,ass_Gamyxy,ass_Gamyxz, &
                               ass_Gamyyy,ass_Gamyyz,ass_Gamyzz, &
                               ass_Gamzxx,ass_Gamzxy,ass_Gamzxz, &
                               ass_Gamzyy,ass_Gamzyz,ass_Gamzzz, &
                       Rxyxy, Rxyxz, Rxyyz, Rxzxz, Rxzyz, Ryzyz)

  implicit none

!~~~~~~ argument variables

  integer,intent(in ) :: ex(1:3)
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: gxxxx,gxxxy,gxxxz,gxxyy,gxxyz,gxxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: gxyxx,gxyxy,gxyxz,gxyyy,gxyyz,gxyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: gxzxx,gxzxy,gxzxz,gxzyy,gxzyz,gxzzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: gyyxx,gyyxy,gyyxz,gyyyy,gyyyz,gyyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: gyzxx,gyzxy,gyzxz,gyzyy,gyzyz,gyzzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: gzzxx,gzzxy,gzzxz,gzzyy,gzzyz,gzzzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamzyy, Gamzyz, Gamzzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: ass_Gamxxx,ass_Gamxxy,ass_Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: ass_Gamxyy,ass_Gamxyz,ass_Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: ass_Gamyxx,ass_Gamyxy,ass_Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: ass_Gamyyy,ass_Gamyyz,ass_Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: ass_Gamzxx,ass_Gamzxy,ass_Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: ass_Gamzyy,ass_Gamzyz,ass_Gamzzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) ::  Rxyxy,  Rxyxz,  Rxyyz,  Rxzxz,  Rxzyz,  Ryzyz

!~~~~~~local variables

  real*8, parameter :: HLF=0.5d0

!R_ijkl = HLF *(@_jk g_il + @_il g_jk - @_jl g_ik - @_ik g_jl)
!               + Gam_rjk Gam^r_il - Gam_rjl Gam^r_ik

  Rxyxy = HLF *( gxyxy + gxyxy - gxxyy - gyyxx )                            + &
          (ass_Gamxxy * Gamxxy + ass_Gamyxy * Gamyxy + ass_Gamzxy * Gamzxy) - &
          (ass_Gamxyy * Gamxxx + ass_Gamyyy * Gamyxx + ass_Gamzyy * Gamzxx) 

  Rxyxz = HLF *( gxzxy + gxyxz - gxxyz - gyzxx )                            + &
          (ass_Gamxxy * Gamxxz + ass_Gamyxy * Gamyxz + ass_Gamzxy * Gamzxz) - &
          (ass_Gamxyz * Gamxxx + ass_Gamyyz * Gamyxx + ass_Gamzyz * Gamzxx) 

  Rxyyz = HLF *( gxzyy + gyyxz - gxyyz - gyzxy )                            + &
          (ass_Gamxyy * Gamxxz + ass_Gamyyy * Gamyxz + ass_Gamzyy * Gamzxz) - &
          (ass_Gamxyz * Gamxxy + ass_Gamyyz * Gamyxy + ass_Gamzyz * Gamzxy) 

  Rxzxz = HLF *( gxzxz + gxzxz - gxxzz - gzzxx )                            + &
          (ass_Gamxxz * Gamxxz + ass_Gamyxz * Gamyxz + ass_Gamzxz * Gamzxz) - &
          (ass_Gamxzz * Gamxxx + ass_Gamyzz * Gamyxx + ass_Gamzzz * Gamzxx) 

  Rxzyz = HLF *( gxzyz + gyzxz - gxyzz - gzzxy )                            + &
          (ass_Gamxyz * Gamxxz + ass_Gamyyz * Gamyxz + ass_Gamzyz * Gamzxz) - &
          (ass_Gamxzz * Gamxxy + ass_Gamyzz * Gamyxy + ass_Gamzzz * Gamzxy) 

  Ryzyz = HLF *( gyzyz + gyzyz - gyyzz - gzzyy )                            + &
          (ass_Gamxyz * Gamxyz + ass_Gamyyz * Gamyyz + ass_Gamzyz * Gamzyz) - &
          (ass_Gamxzz * Gamxyy + ass_Gamyzz * Gamyyy + ass_Gamzzz * Gamzyy) 

  return

  end subroutine adm_riemann

!-----------------------------------------------------------------------------
! Get Ricci tensor of metric g from Riemann tensor
! for adm form
! R_ij = gup^kl * R_ikjl
!-----------------------------------------------------------------------------

  subroutine adm_ricci(ex, gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, &
                           Rxyxy, Rxyxz, Rxyyz, Rxzxz, Rxzyz, Ryzyz, &
                             Rxx,   Rxy,   Rxz,   Ryy,   Ryz,   Rzz)
  implicit none

!~~~~~~> Input parameters:
 
  integer,                              intent(in ) :: ex(1:3)
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: gupyy,gupyz,gupzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: Rxyxy, Rxyxz, Rxyyz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: Rxzxz, Rxzyz, Ryzyz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out) :: Rxx,Rxy,Rxz,Ryy,Ryz,Rzz

  Rxx =   gupyy * Rxyxy + gupyz * Rxyxz + gupyz * Rxyxz + gupzz * Rxzxz             
  Rxy = - gupxy * Rxyxy + gupyz * Rxyyz - gupxz * Rxyxz + gupzz * Rxzyz
  Rxz = - gupxy * Rxyxz - gupyy * Rxyyz - gupxz * Rxzxz - gupyz * Rxzyz
  Ryy =   gupxx * Rxyxy - gupxz * Rxyyz - gupxz * Rxyyz + gupzz * Ryzyz 
  Ryz =   gupxx * Rxyxz + gupxy * Rxyyz - gupxz * Rxzyz - gupyz * Ryzyz
  Rzz =   gupxx * Rxzxz + gupxy * Rxzyz + gupxy * Rxzyz + gupyy * Ryzyz 

  return

  end subroutine adm_ricci

!-----------------------------------------------------------------------------
! raise index
!-----------------------------------------------------------------------------

  subroutine raise(ex,fx,fy,fz,fupx,fupy,fupz, &
                   gupxx,gupxy,gupxz,gupyy,gupyz,gupzz)
  implicit none

!~~~~~~ Input parameters:

  integer,                              intent(in ) :: ex(1:3)
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: fx,fy,fz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: gupxx, gupxy, gupxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: gupyy, gupyz, gupzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out) :: fupx,fupy,fupz

  fupx = gupxx * fx + gupxy * fy + gupxz * fz
  fupy = gupxy * fx + gupyy * fy + gupyz * fz
  fupz = gupxz * fx + gupyz * fy + gupzz * fz

  return

  end subroutine raise

!-----------------------------------------------------------------------------
! lower index
!-----------------------------------------------------------------------------

  subroutine lower(ex,fx,fy,fz,Lfx,Lfy,Lfz,gxx,gxy,gxz,gyy,gyz,gzz)
  implicit none

!~~~~~~ Input parameters:

  integer,                              intent(in ) :: ex(1:3)
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: fx,fy,fz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out) :: Lfx,Lfy,Lfz

  Lfx = gxx * fx + gxy * fy + gxz * fz
  Lfy = gxy * fx + gyy * fy + gyz * fz
  Lfz = gxz * fx + gyz * fy + gzz * fz

  return

  end subroutine lower

!----------------------------------------------------------------------------------
!  inner product of two three dimensional vectors with metric g_ij
! metric here do not upto ONE
!----------------------------------------------------------------------------------

  subroutine InnerProd(ex,norm,ux,uy,uz,vx,vy,vz,gxx,gxy,gxz,gyy,gyz,gzz)
  implicit none

!~~~~~~ argument variables

  integer,intent(in ):: ex(1:3)
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in)::ux,uy,uz,vx,vy,vz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in)::gxx,gxy,gxz,gyy,gyz,gzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out)::norm

  norm =   gxx * ux * vx + gxy * ux * vy + gxz * ux * vz &
         + gxy * uy * vx + gyy * uy * vy + gyz * uy * vz &
         + gxz * uz * vx + gyz * uz * vy + gzz * uz * vz

  return

  end subroutine InnerProd
!----------------------------------------------------------!
!                                                          !
!  algebraic computation based on geometric quantites      !
!  and their partial derivatives related to 3-dimensional  !
!  Riemann slice                                           !
!                                                          !
!      * for BSSN form  *                                  !
!----------------------------------------------------------!

!-----------------------------------------------------------------------------
! second order covariant derivatives w.r.t. *untilded* (i.e. physical) metric 
! of *symmetric* variable of scalar field
!-----------------------------------------------------------------------------

  subroutine fnt_cov_s_dderiv(ex,   fx,   fy,   fz,             &
                             fxx,  fxy,  fxz,  fyy,  fyz,  fzz, &
                            phix, phiy, phiz,                   &
                             dxx,  gxy,  gxz,  dyy,  gyz,  dzz, &
                           gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
                           Gmxxx,Gmxxy,Gmxxz,Gmxyy,Gmxyz,Gmxzz, &
                           Gmyxx,Gmyxy,Gmyxz,Gmyyy,Gmyyz,Gmyzz, &
                           Gmzxx,Gmzxy,Gmzxz,Gmzyy,Gmzyz,Gmzzz)
  implicit none

!~~~~~~ Input arguments

  integer,                              intent(in ) :: ex(1:3)
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: fx,fy,fz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: phix,phiy,phiz
! tilted Christofel symble
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: Gmxxx, Gmxxy, Gmxxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: Gmxyy, Gmxyz, Gmxzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: Gmyxx, Gmyxy, Gmyxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: Gmyyy, Gmyyz, Gmyzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: Gmzxx, Gmzxy, Gmzxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: Gmzyy, Gmzyz, Gmzzz 
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: dxx,gxy,gxz,dyy,gyz,dzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: gupyy,gupyz,gupzz
! input partial derivatives, output covariant derivative respect to physical metric
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: fxx,fxy,fxz,fyy,fyz,fzz

!~~~~~~ Other variables:

  real*8, dimension(ex(1),ex(2),ex(3)) :: phiupx,phiupy,phiupz
  real*8,parameter :: TWO = 2.d0

!~~~~~~ Make untilded Gamma's out of tilded ones - first raise index on phi_i...

  phiupx = gupxx * phix + gupxy * phiy + gupxz * phiz
  phiupy = gupxy * phix + gupyy * phiy + gupyz * phiz
  phiupz = gupxz * phix + gupyz * phiy + gupzz * phiz

!~~~~~~ ... and then add reconstructed *untilded* Christofels...

  fxx = fxx - ( Gmxxx + TWO * ( phix + phix - dxx * phiupx - phiupx ))* fx - &
              ( Gmyxx + TWO * (             - dxx * phiupy - phiupy ))* fy - &
              ( Gmzxx + TWO * (             - dxx * phiupz - phiupz ))* fz

  fyy = fyy - ( Gmxyy + TWO * (             - dyy * phiupx - phiupx ))* fx - &
              ( Gmyyy + TWO * ( phiy + phiy - dyy * phiupy - phiupy ))* fy - &
              ( Gmzyy + TWO * (             - dyy * phiupz - phiupz ))* fz

  fzz = fzz - ( Gmxzz + TWO * (             - dzz * phiupx - phiupx ))* fx - &
              ( Gmyzz + TWO * (             - dzz * phiupy - phiupy ))* fy - &
              ( Gmzzz + TWO * ( phiz + phiz - dzz * phiupz - phiupz ))* fz
 
  fxy = fxy - ( Gmxxy + TWO * ( phiy        - gxy * phiupx          ))* fx - &
              ( Gmyxy + TWO * (        phix - gxy * phiupy          ))* fy - &
              ( Gmzxy + TWO * (             - gxy * phiupz          ))* fz

  fxz = fxz - ( Gmxxz + TWO * ( phiz        - gxz * phiupx          ))* fx - &
              ( Gmyxz + TWO * (             - gxz * phiupy          ))* fy - &
              ( Gmzxz + TWO * (        phix - gxz * phiupz          ))* fz 

  fyz = fyz - ( Gmxyz + TWO * (             - gyz * phiupx          ))* fx - &
              ( Gmyyz + TWO * ( phiz        - gyz * phiupy          ))* fy - &
              ( Gmzyz + TWO * (        phiy - gyz * phiupz          ))* fz

  return

  end subroutine fnt_cov_s_dderiv

!-----------------------------------------------------------------------------
!
! Get physical riemann tensor 
!
!-----------------------------------------------------------------------------

  subroutine get_physical_riemann(ex, ep4phi,                                 &
                                 dxx,    gxy,    gxz,    dyy,    gyz,    dzz, &
                               gupxx,  gupxy,  gupxz,  gupyy,  gupyz,  gupzz, &
                                phix,   phiy,   phiz,                         &
                               phixx,  phixy,  phixz,  phiyy,  phiyz,  phizz, &
                              tRxyxy, tRxyxz, tRxyyz, tRxzxz, tRxzyz, tRyzyz, &
                               Rxyxy,  Rxyxz,  Rxyyz,  Rxzxz,  Rxzyz,  Ryzyz)
  implicit none

!~~~~~~> Input parameters:
 
  integer,                              intent(in ):: ex(1:3)
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ):: ep4phi
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ):: dxx,gxy,gxz,dyy,gyz,dzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ):: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ):: gupyy,gupyz,gupzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ):: phix,phiy,phiz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ):: phixx,phixy,phixz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ):: phiyy,phiyz,phizz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ):: tRxyxy,tRxyxz,tRxyyz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ):: tRxzxz,tRxzyz,tRyzyz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out)::  Rxyxy, Rxyxz, Rxyyz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out)::  Rxzxz, Rxzyz, Ryzyz
 
!~~~~~~> Other variables:

  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gyy,gzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: tmp
  real*8,parameter::ONE = 1.d0, TWO = 2.d0, FOUR = 4.d0

  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE

!~~~~~~> R_ijkl = tilde R_ijkl + TWO *( gli * D_j D_k phi - glj * D_i D_k phi -
!                                       gki * D_j D_l phi + gkj * D_i D_l phi )
!                              + FOUR*( gjl * D_i phi * D_k phi - gil * D_j phi * D_k phi -
!                                       gjk * D_i phi * D_l phi + gik * D_j phi * D_l phi )
!                              + FOUR*( gjk * gil - gik * gjl )* g^mn * D_m phi * D_n phi

    tmp =        gupxx * phix * phix + gupyy * phiy * phiy + gupzz * phiz * phiz + &
          TWO *( gupxy * phix * phiy + gupxz * phix * phiz + gupyz * phiy * phiz )

!~~~~~~> R_ijkl = tilde R_ijkl + TWO *( gli * phi_jk - glj * phi_ik - 
!                                       gki * phi_jl + gkj * phi_il )
!                              + FOUR*( gjl * phi_i * phi_k - gil * phi_j * phi_k - 
!                                       gjk * phi_i * phi_l + gik * phi_j * phi_l )
!                              + FOUR*( gjk * gil - gik * gjl )* tmp

  Rxyxy = tRxyxy + TWO *( gxy * phixy - gyy * phixx - gxx * phiyy + gxy * phixy ) &
                 + FOUR*( gyy * phix * phix - gxy * phiy * phix - &
                          gxy * phix * phiy + gxx * phiy * phiy ) &
                 + FOUR*( gxy * gxy - gxx * gyy )* tmp

  Rxyxz = tRxyxz + TWO *( gxz * phixy - gyz * phixx - gxx * phiyz + gxy * phixz ) &
                 + FOUR*( gyz * phix * phix - gxz * phiy * phix - &
                          gxy * phix * phiz + gxx * phiy * phiz ) &
                 + FOUR*( gxy * gxz - gxx * gyz )* tmp

  Rxyyz = tRxyyz + TWO *( gxz * phiyy - gyz * phixy - gxy * phiyz + gyy * phixz ) &
                 + FOUR*( gyz * phix * phiy - gxz * phiy * phiy - &
                          gyy * phix * phiz + gxy * phiy * phiz ) &
                 + FOUR*( gyy * gxz - gxy * gyz )* tmp

  Rxzxz = tRxzxz + TWO *( gxz * phixz - gzz * phixx - gxx * phizz + gxz * phixz ) &
                 + FOUR*( gzz * phix * phix - gxz * phiz * phix - &
                          gxz * phix * phiz + gxx * phiz * phiz ) &
                 + FOUR*( gxz * gxz - gxx * gzz )* tmp

  Rxzyz = tRxzyz + TWO *( gxz * phiyz - gzz * phixy - gxy * phizz + gyz * phixz ) &
                 + FOUR*( gzz * phix * phiy - gxz * phiz * phiy - &
                          gyz * phix * phiz + gxy * phiz * phiz ) &
                 + FOUR*( gyz * gxz - gxy * gzz )* tmp

  Ryzyz = tRyzyz + TWO *( gyz * phiyz - gzz * phiyy - gyy * phizz + gyz * phiyz ) &
                 + FOUR*( gzz * phiy * phiy - gyz * phiz * phiy - &
                          gyz * phiy * phiz + gyy * phiz * phiz ) &
                 + FOUR*( gyz * gyz - gyy * gzz )* tmp

!multipli with factor exp( 4 * phi)

  Rxyxy = Rxyxy * ep4phi
  Rxyxz = Rxyxz * ep4phi
  Rxyyz = Rxyyz * ep4phi
  Rxzxz = Rxzxz * ep4phi
  Rxzyz = Rxzyz * ep4phi
  Ryzyz = Ryzyz * ep4phi

  return

  end subroutine get_physical_riemann

!-----------------------------------------------------------------------------
!
! Get physical Ricci tensor 
!
!-----------------------------------------------------------------------------

  subroutine get_physical_ricci(ex,dxx,gxy,gxz,dyy,gyz,dzz,phix,phiy,phiz, &
                                phixx,phixy,phixz,phiyy,phiyz,phizz, &
                                gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
                                 tRxx, tRxy, tRxz, tRyy, tRyz, tRzz, &
                                  Rxx,  Rxy,  Rxz,  Ryy,  Ryz,  Rzz)

  implicit none

!~~~~~~ argument variables

  integer, intent(in) :: ex(1:3)
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in):: dxx,gxy,gxz,dyy,gyz,dzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in):: phix,phiy,phiz
! covariant derivative respect to tilted metric
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in):: phixx,phixy,phixz,phiyy,phiyz,phizz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in):: gupxx,gupxy,gupxz,gupyy,gupyz,gupzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in) :: tRxx,tRxy,tRxz,tRyy,tRyz,tRzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out) :: Rxx, Rxy, Rxz, Ryy, Ryz, Rzz

!~~~~~~ local variables

  real*8, dimension(ex(1),ex(2),ex(3)) :: tempf
  real*8,parameter::TWO = 2.d0, FOUR = 4.d0

!~~~~~~

  tempf = TWO * (gupxx * ( phixx + TWO * phix * phix ) + &
                 gupyy * ( phiyy + TWO * phiy * phiy ) + &
                 gupzz * ( phizz + TWO * phiz * phiz ) + &
           TWO * gupxy * ( phixy + TWO * phix * phiy ) + &
           TWO * gupxz * ( phixz + TWO * phix * phiz ) + &
           TWO * gupyz * ( phiyz + TWO * phiy * phiz ) )

! Add phi part to Ricci tensor:

  Rxx = tRxx - TWO * phixx + FOUR * phix * phix - dxx * tempf - tempf
  Ryy = tRyy - TWO * phiyy + FOUR * phiy * phiy - dyy * tempf - tempf
  Rzz = tRzz - TWO * phizz + FOUR * phiz * phiz - dzz * tempf - tempf
  Rxy = tRxy - TWO * phixy + FOUR * phix * phiy - gxy * tempf
  Rxz = tRxz - TWO * phixz + FOUR * phix * phiz - gxz * tempf
  Ryz = tRyz - TWO * phiyz + FOUR * phiy * phiz - gyz * tempf

  return

  end subroutine get_physical_ricci

!-----------------------------------------------------------------------------
!
! compute physical extrinic curver:
! Kij = exp( 4 * phi ) ( tilde Aij + F1o3 * tilde gij * trK )
!
!-----------------------------------------------------------------------------

  subroutine get_physical_k(ex, phi, trK, dxx, gxy, gxz, dyy, gyz, dzz, &
                                          Axx, Axy, Axz, Ayy, Ayz, Azz, &
                                          Kxx, Kxy, Kxz, Kyy, Kyz, Kzz)
  implicit none

!~~~~~~> Input parameters:

  integer,dimension(3)                , intent(in) :: ex
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in) :: phi, trK
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in) :: dxx,gxy,gxz,dyy,gyz,dzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out):: Kxx,Kxy,Kxz,Kyy,Kyz,Kzz

!~~~~~~> Other variables: 

  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx, gyy, gzz
  real*8, parameter :: F1o3 = 1.d0 / 3.d0, ONE = 1.d0, FOUR = 4.d0

!~~~~~~>

  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE

  Kzz = exp( FOUR * phi )

!~~~~~~>

  Kxx = ( Axx + F1o3 * gxx * trK )* Kzz
  Kxy = ( Axy + F1o3 * gxy * trK )* Kzz
  Kxz = ( Axz + F1o3 * gxz * trK )* Kzz
  Kyy = ( Ayy + F1o3 * gyy * trK )* Kzz
  Kyz = ( Ayz + F1o3 * gyz * trK )* Kzz
  Kzz = ( Azz + F1o3 * gzz * trK )* Kzz

  return

  end subroutine get_physical_k

!-------------------------------------------------------------------------------------------------------
!
! compute covariant derivatives of extrinic curver
!
!D_i K_jk stored as DKijk
!
!   DKijk =   e^(4 phi) (A_jk,i - Gam^l_ij A_lk - Gam^l_ik A_jl + 1/3 g_jk trK,i)
!           - 2 K_ik phi,j + 2 g_ij g^lm phi,m K_lk
!           - 2 K_ij phi,k + 2 g_ik g^lm phi,m K_lj
!-------------------------------------------------------------------------------------------------------

  subroutine get_diff_physical_k(ex, phi, trK, Kx, Ky, Kz, phix, phiy, phiz, &
                                 dxx,  gxy,  gxz,  dyy,  gyz,  dzz, &
                                 gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
                                 Axx,  Axy,  Axz,  Ayy,  Ayz,  Azz, &
                                 Axxx, Axxy, Axxz, Axyx, Axyy, Axyz, &
                                 Axzx, Axzy, Axzz, Ayyx, Ayyy, Ayyz, &
                                 Ayzx, Ayzy, Ayzz, Azzx, Azzy, Azzz, &
                                 Gmxxx,Gmxxy,Gmxxz,Gmxyy,Gmxyz,Gmxzz, &
                                 Gmyxx,Gmyxy,Gmyxz,Gmyyy,Gmyyz,Gmyzz, &
                                 Gmzxx,Gmzxy,Gmzxz,Gmzyy,Gmzyz,Gmzzz, &
                                 Kxx,  Kxy,  Kxz,  Kyy,  Kyz,  Kzz, &
                                 DKxxx,DKxxy,DKxxz,DKxyy,DKxyz,DKxzz, &
                                 DKyxx,DKyxy,DKyxz,DKyyy,DKyyz,DKyzz, &
                                 DKzxx,DKzxy,DKzxz,DKzyy,DKzyz,DKzzz)

  implicit none

!~~~~~~> Input parameters:

  integer,dimension(3), intent(in) :: ex
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ):: phi,trK
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ):: Kx,Ky,Kz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ):: phix,phiy,phiz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ):: dxx,gxy,gxz,dyy,gyz,dzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ):: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ):: gupyy,gupyz,gupzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ):: Axx,Axy,Axz,Ayy,Ayz,Azz
! Aij,k --> stored as Aijk
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ):: Axxx,Axxy,Axxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ):: Axyx,Axyy,Axyz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ):: Axzx,Axzy,Axzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ):: Ayyx,Ayyy,Ayyz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ):: Ayzx,Ayzy,Ayzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ):: Azzx,Azzy,Azzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ):: Gmxxx,Gmxxy,Gmxxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ):: Gmxyy,Gmxyz,Gmxzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ):: Gmyxx,Gmyxy,Gmyxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ):: Gmyyy,Gmyyz,Gmyzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ):: Gmzxx,Gmzxy,Gmzxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ):: Gmzyy,Gmzyz,Gmzzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ):: Kxx,Kxy,Kxz,Kyy,Kyz,Kzz
! D_i K_jk --> stored as DKijk
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out):: DKxxx,DKxxy,DKxxz,DKxyy,DKxyz,DKxzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out):: DKyxx,DKyxy,DKyxz,DKyyy,DKyyz,DKyzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out):: DKzxx,DKzxy,DKzxz,DKzyy,DKzyz,DKzzz

!~~~~~~> Other variables:

  real*8, dimension(ex(1),ex(2),ex(3)):: phiupx,phiupy,phiupz
  real*8, dimension(ex(1),ex(2),ex(3)):: phiupKx,phiupKy,phiupKz
  real*8, dimension(ex(1),ex(2),ex(3)):: e4phi
  real*8, dimension(ex(1),ex(2),ex(3)):: gxx,gyy,gzz

  real*8,parameter::ONE = 1.d0, TWO = 2.d0, FOUR = 4.d0
  real*8,parameter::F1o3 = 1.d0/3.d0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0

!~~~~~~> Input translation

   gxx = dxx + ONE
   gyy = dyy + ONE
   gzz = dzz + ONE

   e4phi = dexp(FOUR * phi)

!~~~~~~>

   phiupx = gupxx * phix + gupxy * phiy + gupxz * phiz
   phiupy = gupxy * phix + gupyy * phiy + gupyz * phiz
   phiupz = gupxz * phix + gupyz * phiy + gupzz * phiz

   phiupKx = phiupx * Kxx + phiupy * Kxy + phiupz * Kxz
   phiupKy = phiupx * Kxy + phiupy * Kyy + phiupz * Kyz
   phiupKz = phiupx * Kxz + phiupy * Kyz + phiupz * Kzz

!~~~~~~> tmp = - Gam^l_ij A_lk - Gam^l_ik A_jl

  DKxxx = - Gmxxx * Axx - Gmyxx * Axy - Gmzxx * Axz &
          - Gmxxx * Axx - Gmyxx * Axy - Gmzxx * Axz 

  DKxxy = - Gmxxx * Axy - Gmyxx * Ayy - Gmzxx * Ayz &
          - Gmxxy * Axx - Gmyxy * Axy - Gmzxy * Axz 

  DKxxz = - Gmxxx * Axz - Gmyxx * Ayz - Gmzxx * Azz &
          - Gmxxz * Axx - Gmyxz * Axy - Gmzxz * Axz 

  DKxyy = - Gmxxy * Axy - Gmyxy * Ayy - Gmzxy * Ayz &
          - Gmxxy * Axy - Gmyxy * Ayy - Gmzxy * Ayz 

  DKxyz = - Gmxxy * Axz - Gmyxy * Ayz - Gmzxy * Azz &
          - Gmxxz * Axy - Gmyxz * Ayy - Gmzxz * Ayz 

  DKxzz = - Gmxxz * Axz - Gmyxz * Ayz - Gmzxz * Azz &
          - Gmxxz * Axz - Gmyxz * Ayz - Gmzxz * Azz 

  DKyxx = - Gmxxy * Axx - Gmyxy * Axy - Gmzxy * Axz &
          - Gmxxy * Axx - Gmyxy * Axy - Gmzxy * Axz 

  DKyxy = - Gmxxy * Axy - Gmyxy * Ayy - Gmzxy * Ayz &
          - Gmxyy * Axx - Gmyyy * Axy - Gmzyy * Axz 

  DKyxz = - Gmxxy * Axz - Gmyxy * Ayz - Gmzxy * Azz &
          - Gmxyz * Axx - Gmyyz * Axy - Gmzyz * Axz 

  DKyyy = - Gmxyy * Axy - Gmyyy * Ayy - Gmzyy * Ayz &
          - Gmxyy * Axy - Gmyyy * Ayy - Gmzyy * Ayz 

  DKyyz = - Gmxyy * Axz - Gmyyy * Ayz - Gmzyy * Azz &
          - Gmxyz * Axy - Gmyyz * Ayy - Gmzyz * Ayz 

  DKyzz = - Gmxyz * Axz - Gmyyz * Ayz - Gmzyz * Azz &
          - Gmxyz * Axz - Gmyyz * Ayz - Gmzyz * Azz 

  DKzxx = - Gmxxz * Axx - Gmyxz * Axy - Gmzxz * Axz &
          - Gmxxz * Axx - Gmyxz * Axy - Gmzxz * Axz 

  DKzxy = - Gmxxz * Axy - Gmyxz * Ayy - Gmzxz * Ayz &
          - Gmxyz * Axx - Gmyyz * Axy - Gmzyz * Axz 

  DKzxz = - Gmxxz * Axz - Gmyxz * Ayz - Gmzxz * Azz &
          - Gmxzz * Axx - Gmyzz * Axy - Gmzzz * Axz 

  DKzyy = - Gmxyz * Axy - Gmyyz * Ayy - Gmzyz * Ayz &
          - Gmxyz * Axy - Gmyyz * Ayy - Gmzyz * Ayz 

  DKzyz = - Gmxyz * Axz - Gmyyz * Ayz - Gmzyz * Azz &
          - Gmxzz * Axy - Gmyzz * Ayy - Gmzzz * Ayz 

  DKzzz = - Gmxzz * Axz - Gmyzz * Ayz - Gmzzz * Azz &
          - Gmxzz * Axz - Gmyzz * Ayz - Gmzzz * Azz 

!~~~~~~>   DKijk =   e^(4 phi) (A_jk,i + tmp + 1/3 g_jk K_i)
!                  - 2 K_ik phi,j + 2 g_ij phiupK_k
!                  - 2 K_ij phi,k + 2 g_ik phiupK_j

  DKxxx =    e4phi * (Axxx + DKxxx  + F1o3 * gxx * Kx) &
           - TWO * Kxx * phix + TWO * gxx * phiupKx &
           - TWO * Kxx * phix + TWO * gxx * phiupKx

  DKxxy =    e4phi * (Axyx + DKxxy  + F1o3 * gxy * Kx) &
           - TWO * Kxy * phix + TWO * gxx * phiupKy &
           - TWO * Kxx * phiy + TWO * gxy * phiupKx

  DKxxz =   e4phi * (Axzx + DKxxz  + F1o3 * gxz * Kx) &
           - TWO * Kxz * phix + TWO * gxx * phiupKz &
           - TWO * Kxx * phiz + TWO * gxz * phiupKx

  DKxyy =    e4phi * (Ayyx + DKxyy  + F1o3 * gyy * Kx) &
           - TWO * Kxy * phiy + TWO * gxy * phiupKy &
           - TWO * Kxy * phiy + TWO * gxy * phiupKy

  DKxyz =    e4phi * (Ayzx + DKxyz  + F1o3 * gyz * Kx) &
           - TWO * Kxz * phiy + TWO * gxy * phiupKz &
           - TWO * Kxy * phiz + TWO * gxz * phiupKy

  DKxzz =    e4phi * (Azzx + DKxzz  + F1o3 * gzz * Kx) &
           - TWO * Kxz * phiz + TWO * gxz * phiupKz &
           - TWO * Kxz * phiz + TWO * gxz * phiupKz

!~~~~~~>

  DKyxx =    e4phi * (Axxy + DKyxx  + F1o3 * gxx * Ky) &
           - TWO * Kxy * phix + TWO * gxy * phiupKx &
           - TWO * Kxy * phix + TWO * gxy * phiupKx

  DKyxy =    e4phi * (Axyy + DKyxy  + F1o3 * gxy * Ky) &
           - TWO * Kyy * phix + TWO * gxy * phiupKy &
           - TWO * Kxy * phiy + TWO * gyy * phiupKx

  DKyxz =    e4phi * (Axzy + DKyxz  + F1o3 * gxz * Ky) &
           - TWO * Kyz * phix + TWO * gxy * phiupKz &
           - TWO * Kxy * phiz + TWO * gyz * phiupKx

  DKyyy =    e4phi * (Ayyy + DKyyy  + F1o3 * gyy * Ky) &
           - TWO * Kyy * phiy + TWO * gyy * phiupKy &
           - TWO * Kyy * phiy + TWO * gyy * phiupKy

  DKyyz =    e4phi * (Ayzy + DKyyz  + F1o3 * gyz * Ky) &
           - TWO * Kyz * phiy + TWO * gyy * phiupKz &
           - TWO * Kyy * phiz + TWO * gyz * phiupKy

  DKyzz =    e4phi * (Azzy + DKyzz  + F1o3 * gzz * Ky) &
           - TWO * Kyz * phiz + TWO * gyz * phiupKz &
           - TWO * Kyz * phiz + TWO * gyz * phiupKz

!~~~~~~>

  DKzxx =    e4phi * (Axxz + DKzxx  + F1o3 * gxx * Kz) &
           - TWO * Kxz * phix + TWO * gxz * phiupKx &
           - TWO * Kxz * phix + TWO * gxz * phiupKx

  DKzxy =    e4phi * (Axyz + DKzxy  + F1o3 * gxy * Kz) &
           - TWO * Kyz * phix + TWO * gxz * phiupKy &
           - TWO * Kxz * phiy + TWO * gyz * phiupKx

  DKzxz =    e4phi * (Axzz + DKzxz  + F1o3 * gxz * Kz) &
           - TWO * Kzz * phix + TWO * gxz * phiupKz &
           - TWO * Kxz * phiz + TWO * gzz * phiupKx

  DKzyy =    e4phi * (Ayyz + DKzyy  + F1o3 * gyy * Kz) &
           - TWO * Kyz * phiy + TWO * gyz * phiupKy &
           - TWO * Kyz * phiy + TWO * gyz * phiupKy

  DKzyz =    e4phi * (Ayzz + DKzyz  + F1o3 * gyz * Kz) &
           - TWO * Kzz * phiy + TWO * gyz * phiupKz &
           - TWO * Kyz * phiz + TWO * gzz * phiupKy

  DKzzz =    e4phi * (Azzz + DKzzz  + F1o3 * gzz * Kz) &
           - TWO * Kzz * phiz + TWO * gzz * phiupKz &
           - TWO * Kzz * phiz + TWO * gzz * phiupKz

  return

  end subroutine get_diff_physical_k

!----------------------------------------------------------------------
!------>Begin to compute Psi4
!------>based on quantites:
!------>triad v^i, u^i, w^i
!------>lapse and shift vector beta^i
!------>extrinsic curvature K_ij and trK
!------>covariant derivative of extrinsic curvature D_i K_jk
!------>Ricci tensor: R_ij
!------>gup^ij
!------>Riemann tensor R_ijkl
!----------------------------------------------------------------------

  subroutine bssn_compute_psi4(ex, em4phi,lapse, betax,betay,betaz, &
                          gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,      &
                          vx,vy,vz,ux,uy,uz,wx,wy,wz,               &
                          trK,Kxx,Kxy,Kxz,Kyy,Kyz,Kzz,              &
                          Rxyxy, Rxyxz, Rxyyz, Rxzxz, Rxzyz, Ryzyz, &
                          Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,             &
                          DKxxx,DKxxy,DKxxz,DKxyy,DKxyz,DKxzz,      &
                          DKyxx,DKyxy,DKyxz,DKyyy,DKyyz,DKyzz,      &
                          DKzxx,DKzxy,DKzxz,DKzyy,DKzyz,DKzzz, Rpsi4, Ipsi4)

  implicit none

!~~~~~~ argument variables

  integer,intent(in ):: ex(1:3)
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: em4phi
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: lapse
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: betax, betay, betaz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: gupyy,gupyz,gupzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: vx,vy,vz,ux,uy,uz,wx,wy,wz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: trK,Kxx,Kxy,Kxz,Kyy,Kyz,Kzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Rxyxy, Rxyxz, Rxyyz, Rxzxz, Rxzyz, Ryzyz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Rxx, Rxy, Rxz, Ryy, Ryz, Rzz
!D_i K_jk ---> DKijk
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: DKxxx,DKxxy,DKxxz,DKxyy,DKxyz,DKxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: DKyxx,DKyxy,DKyxz,DKyyy,DKyyz,DKyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: DKzxx,DKzxy,DKzxz,DKzyy,DKzyz,DKzzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Rpsi4,Ipsi4

!~~~~~~ local variables

!n^i upto 1/sqrt(2)
  real*8, dimension(ex(1),ex(2),ex(3)) :: nx,ny,nz
!n^i * n^k upto 1/2
  real*8, dimension(ex(1),ex(2),ex(3)) :: nnxx,nnxy,nnxz,nnyy,nnyz,nnzz
!u^j * u^l - w^j * w^l
  real*8, dimension(ex(1),ex(2),ex(3)) :: uuwwxx,uuwwxy,uuwwxz,uuwwyy,uuwwyz,uuwwzz
!- u^j * w^l - w^j * u^l
  real*8, dimension(ex(1),ex(2),ex(3)) :: uwxx,uwxy,uwxz,uwyy,uwyz,uwzz
! temp variables
  real*8, dimension(ex(1),ex(2),ex(3)) ::temRxx, temRxy, temRxz, temRyy, temRyz, temRzz
  real*8, dimension(ex(1),ex(2),ex(3)) ::temRxyxy,temRxyxz,temRxyyz,temRxzxz,temRxzyz,temRyzyz
  real*8, dimension(ex(1),ex(2),ex(3)) ::lapse2
! K^i_j
  real*8, dimension(ex(1),ex(2),ex(3)) ::Kupxx,Kupxy,Kupxz,Kupyy,Kupyz,Kupzz

  real*8, parameter :: TWO = 2.d0, F1o4 = 1.d0/4.d0

!~~~~~~

! compute n^i = - beta^i/lapse - v^i
  nx = - betax/lapse - vx
  ny = - betay/lapse - vy
  nz = - betaz/lapse - vz

! compute nn^ij = n^i * n^j
  nnxx = nx * nx
  nnxy = nx * ny
  nnxz = nx * nz
  nnyy = ny * ny
  nnyz = ny * nz
  nnzz = nz * nz

! compute uuww^ij = u^i * u^j - w^i * w^j
  uuwwxx = ux * ux - wx * wx
  uuwwxy = ux * uy - wx * wy
  uuwwxz = ux * uz - wx * wz
  uuwwyy = uy * uy - wy * wy
  uuwwyz = uy * uz - wy * wz
  uuwwzz = uz * uz - wz * wz

! compute uw^ij = - u^i * w^j - w^i * u^j
  uwxx = ux * wx + wx * ux
  uwxy = ux * wy + wx * uy
  uwxz = ux * wz + wx * uz
  uwyy = uy * wy + wy * uy
  uwyz = uy * wz + wy * uz
  uwzz = uz * wz + wz * uz

!Commonterm_jl = -1/4 * (  (R_ijkl + K_ik * K_jl - K_il * K_jk) * nn^ik
!                        - 2 * (D_l K_jk - D_k K_jl) * n^0 * n^k
!                        + (R_jl - K_jm * K^m_l + K * K_jl) * n^0 * n^0
!                       )

!add trK * K_jl to R_jl
  temRxx = Rxx + trK * Kxx
  temRxy = Rxy + trK * Kxy
  temRxz = Rxz + trK * Kxz
  temRyy = Ryy + trK * Kyy
  temRyz = Ryz + trK * Kyz
  temRzz = Rzz + trK * Kzz

!add - K_jm * K^m_l to R_jl

! compute K^m_l
  call raise(ex,Kxx,Kxy,Kxz,Kupxx,Kupxy,Kupxz, &
             gupxx,gupxy,gupxz,gupyy,gupyz,gupzz)

  call raise(ex,Kxy,Kyy,Kyz,Kupxy,Kupyy,Kupyz, &
             gupxx,gupxy,gupxz,gupyy,gupyz,gupzz)  

  call raise(ex,Kxz,Kyz,Kzz,Kupxz,Kupyz,Kupzz, &
             gupxx,gupxy,gupxz,gupyy,gupyz,gupzz)

  temRxx = temRxx - em4phi * ( Kupxx * Kxx + Kupxy * Kxy + Kupxz * Kxz  )

  temRxy = temRxy - em4phi * ( Kupxx * Kxy + Kupxy * Kyy + Kupxz * Kyz  )

  temRxz = temRxz - em4phi * ( Kupxx * Kxz + Kupxy * Kyz + Kupxz * Kzz  )

  temRyy = temRyy - em4phi * ( Kupxy * Kxy + Kupyy * Kyy + Kupyz * Kyz  )

  temRyz = temRyz - em4phi * ( Kupxy * Kxz + Kupyy * Kyz + Kupyz * Kzz  )

  temRzz = temRzz - em4phi * ( Kupxz * Kxz + Kupyz * Kyz + Kupzz * Kzz  )

! multiply with n^0 * n^0 upto 1/2
! n^0 = 1/(sqrt(2) * lapse)
  lapse2 = lapse * lapse

  temRxx = temRxx/lapse2
  temRxy = temRxy/lapse2
  temRxz = temRxz/lapse2
  temRyy = temRyy/lapse2
  temRyz = temRyz/lapse2
  temRzz = temRzz/lapse2

!add (K_ik * K_jl - K_il * K_jk) to R_ijkl, note they have the same symmetric index

  temRxyxy = Rxyxy + Kxx * Kyy - Kxy * Kxy
  temRxyxz = Rxyxz + Kxx * Kyz - Kxz * Kxy
  temRxyyz = Rxyyz + Kxy * Kyz - Kxz * Kyy
  temRxzxz = Rxzxz + Kxx * Kzz - Kxz * Kxz
  temRxzyz = Rxzyz + Kxy * Kzz - Kxz * Kyz
  temRyzyz = Ryzyz + Kyy * Kzz - Kyz * Kyz

!add (R_ijkl + K_ik * K_jl - K_il * K_jk) * nn^ik to R_jl, upto 1/2
! note they have the same symmetric index
  temRxx = temRxx + temRxyxy * nnyy + temRxyxz * nnyz + temRxyxz * nnyz + temRxzxz * nnzz
  temRxy = temRxy - temRxyxy * nnxy + temRxyyz * nnyz - temRxyxz * nnxz + temRxzyz * nnzz
  temRxz = temRxz - temRxyxz * nnxy - temRxyyz * nnyy - temRxzxz * nnxz - temRxzyz * nnyz
  temRyy = temRyy + temRxyxy * nnxx - temRxyyz * nnxz - temRxyyz * nnxz + temRyzyz * nnzz
  temRyz = temRyz + temRxyxz * nnxx + temRxyyz * nnxy - temRxzyz * nnxz - temRyzyz * nnyz
  temRzz = temRzz + temRxzxz * nnxx + temRxzyz * nnxy + temRxzyz * nnxy + temRyzyz * nnyy

!add 2 * (D_k K_jl  * n^0 * n^k) to R_jl, upto 1/2
  temRxx = temRxx + TWO * ( DKxxx * nx + DKyxx * ny + DKzxx * nz)/lapse
  temRxy = temRxy + TWO * ( DKxxy * nx + DKyxy * ny + DKzxy * nz)/lapse
  temRxz = temRxz + TWO * ( DKxxz * nx + DKyxz * ny + DKzxz * nz)/lapse
  temRyy = temRyy + TWO * ( DKxyy * nx + DKyyy * ny + DKzyy * nz)/lapse
  temRyz = temRyz + TWO * ( DKxyz * nx + DKyyz * ny + DKzyz * nz)/lapse
  temRzz = temRzz + TWO * ( DKxzz * nx + DKyzz * ny + DKzzz * nz)/lapse

!add - (D_l K_jk + D_j K_lk) * n^0 * ^k to R_jl, upto 1/2
! note we symmetrize the index here
  temRxx = temRxx - ((DKxxx + DKxxx) * nx + (DKxxy + DKxxy) * ny + (DKxxz + DKxxz) * nz)/lapse
  temRxy = temRxy - ((DKyxx + DKxxy) * nx + (DKyxy + DKxyy) * ny + (DKyxz + DKxyz) * nz)/lapse
  temRxz = temRxz - ((DKzxx + DKxxz) * nx + (DKzxy + DKxyz) * ny + (DKzxz + DKxzz) * nz)/lapse
  temRyy = temRyy - ((DKyxy + DKyxy) * nx + (DKyyy + DKyyy) * ny + (DKyyz + DKyyz) * nz)/lapse
  temRyz = temRyz - ((DKzxy + DKyxz) * nx + (DKzyy + DKyyz) * ny + (DKzyz + DKyzz) * nz)/lapse
  temRzz = temRzz - ((DKzxz + DKzxz) * nx + (DKzyz + DKzyz) * ny + (DKzzz + DKzzz) * nz)/lapse

!the real part of Psi4
  Rpsi4 =   temRxx * uuwwxx + temRyy * uuwwyy + temRzz * uuwwzz &
          + (temRxy * uuwwxy + temRxz * uuwwxz + temRyz * uuwwyz) * TWO

!the imaginary part of Psi4
  Ipsi4 =   temRxx * uwxx + temRyy * uwyy + temRzz * uwzz &
          + (temRxy * uwxy + temRxz * uwxz + temRyz * uwyz) * TWO

!multiply with -1/4 
  Rpsi4 = - F1o4 * Rpsi4
  Ipsi4 = - F1o4 * Ipsi4

  return

  end subroutine bssn_compute_psi4
!-----------------------------------------------------------------------------
! covariant derivatives w.r.t *tilded metric* of *symmetric* variable
!-----------------------------------------------------------------------------

  subroutine xcov_deriv(ex,fx,fy,fz,fxx,fxy,fxz,fyy,fyz,fzz,            &
                        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz, &
                        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz, &
                        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz)
  implicit none

!~~~~~~ Input arguments

  integer,                              intent(in ) :: ex(1:3)
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: fx, fy, fz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out) :: fxx,fxy,fxz,fyy,fyz,fzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in ) :: Gamzyy, Gamzyz, Gamzzz 

!~~~~~~ Add Connection terms

  fxx = fxx - Gamxxx * fx - Gamyxx * fy - Gamzxx * fz
  fxy = fxy - Gamxxy * fx - Gamyxy * fy - Gamzxy * fz
  fxz = fxz - Gamxxz * fx - Gamyxz * fy - Gamzxz * fz
  fyy = fyy - Gamxyy * fx - Gamyyy * fy - Gamzyy * fz
  fyz = fyz - Gamxyz * fx - Gamyyz * fy - Gamzyz * fz
  fzz = fzz - Gamxzz * fx - Gamyzz * fy - Gamzzz * fz

  return

  end subroutine xcov_deriv

!--------------------------------------------------------------------
!    Gram-Schmidt orthonormal in Cartesin coordinate
!    V1 = [ x, y, z  ]
!	 V2 = [-y, x, ZEO]
!	 V3 = [xz,yz,-(x^2+y^2)]
!	 V1 -> V1 / sqrt(W11)
!	 V2 -> ( V2 - V1 * W12 ) / sqrt(W22)
!	 V3 -> ( V3 - V1 * W13 - V2 * W23 ) / sqrt(W33)
!	 W_ij = g_ab * Vi^a * Vj^b
!	 it is metric, not tilde metric
! V1 first
!--------------------------------------------------------------------

  subroutine get_triad0(ex, X, Y, Z, ep4phi, &
                       gxxi,gxyi,gxzi,gyyi,gyzi,gzzi, &
                       vx,vy,vz,ux,uy,uz,wx,wy,wz)
  
  implicit none

!~~~~~~ argument variables

  integer,intent(in ):: ex(1:3)
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: ep4phi  !exp(4 * phi)
! tilted metric
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: gxxi,gxyi,gxzi,gyyi,gyzi,gzzi
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: vx,vy,vz,ux,uy,uz,wx,wy,wz

!~~~~~~ local variables

  real*8, dimension(ex(1),ex(2),ex(3)) :: norm
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gxy,gxz,gyy,gyz,gzz

  real*8,parameter:: ZEO = 0.d0, ONE = 1.d0
  integer::i,j,k

!~~~~~~

  gxx = gxxi * ep4phi
  gxy = gxyi * ep4phi
  gxz = gxzi * ep4phi
  gyy = gyyi * ep4phi
  gyz = gyzi * ep4phi
  gzz = gzzi * ep4phi

! initialize U, V, W vetors	
  do i=1,ex(1)
   do j=1,ex(2)
    do k=1,ex(3)

     vx(i,j,k) = X(i)
     vy(i,j,k) = Y(j)
     vz(i,j,k) = Z(k)
     ux(i,j,k) = - Y(j)
     uy(i,j,k) = X(i)
     uz(i,j,k) = ZEO
     wx(i,j,k) = X(i)*Z(k)
     wy(i,j,k) = Y(j)*Z(k)
     wz(i,j,k) = -(X(i)*X(i) + Y(j)*Y(j))

    enddo
   enddo
  enddo

! Gram-Schmidt orthonormalization
  call InnerProd(ex,norm,vx,vy,vz,vx,vy,vz,gxx,gxy,gxz,gyy,gyz,gzz)
  vx = vx/dsqrt(norm)
  vy = vy/dsqrt(norm)
  vz = vz/dsqrt(norm)

  call InnerProd(ex,norm,ux,uy,uz,vx,vy,vz,gxx,gxy,gxz,gyy,gyz,gzz)
  ux = ux - norm*vx
  uy = uy - norm*vy
  uz = uz - norm*vz

  call InnerProd(ex,norm,ux,uy,uz,ux,uy,uz,gxx,gxy,gxz,gyy,gyz,gzz)
  ux = ux/dsqrt(norm)
  uy = uy/dsqrt(norm)
  uz = uz/dsqrt(norm)

  call InnerProd(ex,norm,wx,wy,wz,vx,vy,vz,gxx,gxy,gxz,gyy,gyz,gzz)
  wx = wx - norm*vx
  wy = wy - norm*vy
  wz = wz - norm*vz

  call InnerProd(ex,norm,wx,wy,wz,ux,uy,uz,gxx,gxy,gxz,gyy,gyz,gzz)
  wx = wx - norm*ux
  wy = wy - norm*uy
  wz = wz - norm*uz

  call InnerProd(ex,norm,wx,wy,wz,wx,wy,wz,gxx,gxy,gxz,gyy,gyz,gzz)
  wx = wx/dsqrt(norm)
  wy = wy/dsqrt(norm)
  wz = wz/dsqrt(norm)

  return

  end subroutine get_triad0
!--------------------------------------------------------------------
!    Gram-Schmidt orthonormal in Cartesin coordinate
!    V1 = [ x, y, z  ]
!	 V2 = [-y, x, ZEO]
!	 V3 = [xz,yz,-(x^2+y^2)]
!	 V1 -> V1 / sqrt(W11)
!	 V2 -> ( V2 - V1 * W12 ) / sqrt(W22)
!	 V3 -> ( V3 - V1 * W13 - V2 * W23 ) / sqrt(W33)
!	 W_ij = g_ab * Vi^a * Vj^b
!	 it is metric, not tilde metric
! V2 first
!--------------------------------------------------------------------

  subroutine get_triad1(ex, X, Y, Z, ep4phi, &
                       gxxi,gxyi,gxzi,gyyi,gyzi,gzzi, &
                       vx,vy,vz,ux,uy,uz,wx,wy,wz)
  
  implicit none

!~~~~~~ argument variables

  integer,intent(in ):: ex(1:3)
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: ep4phi  !exp(4 * phi)
! tilted metric
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: gxxi,gxyi,gxzi,gyyi,gyzi,gzzi
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: vx,vy,vz,ux,uy,uz,wx,wy,wz

!~~~~~~ local variables

  real*8, dimension(ex(1),ex(2),ex(3)) :: norm
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gxy,gxz,gyy,gyz,gzz

  real*8,parameter:: ZEO = 0.d0, ONE = 1.d0
  integer::i,j,k

!~~~~~~

  gxx = gxxi * ep4phi
  gxy = gxyi * ep4phi
  gxz = gxzi * ep4phi
  gyy = gyyi * ep4phi
  gyz = gyzi * ep4phi
  gzz = gzzi * ep4phi

! initialize U, V, W vetors	
  do i=1,ex(1)
   do j=1,ex(2)
    do k=1,ex(3)

     vx(i,j,k) = X(i)
     vy(i,j,k) = Y(j)
     vz(i,j,k) = Z(k)
     ux(i,j,k) = - Y(j)
     uy(i,j,k) = X(i)
     uz(i,j,k) = ZEO
     wx(i,j,k) = X(i)*Z(k)
     wy(i,j,k) = Y(j)*Z(k)
     wz(i,j,k) = -(X(i)*X(i) + Y(j)*Y(j))

    enddo
   enddo
  enddo

! Gram-Schmidt orthonormalization
  call InnerProd(ex,norm,ux,uy,uz,ux,uy,uz,gxx,gxy,gxz,gyy,gyz,gzz)
  ux = ux/dsqrt(norm)
  uy = uy/dsqrt(norm)
  uz = uz/dsqrt(norm)

  call InnerProd(ex,norm,vx,vy,vz,ux,uy,uz,gxx,gxy,gxz,gyy,gyz,gzz)
  vx = vx - norm*ux
  vy = vy - norm*uy
  vz = vz - norm*uz

  call InnerProd(ex,norm,vx,vy,vz,vx,vy,vz,gxx,gxy,gxz,gyy,gyz,gzz)
  vx = vx/dsqrt(norm)
  vy = vy/dsqrt(norm)
  vz = vz/dsqrt(norm)

  call InnerProd(ex,norm,wx,wy,wz,ux,uy,uz,gxx,gxy,gxz,gyy,gyz,gzz)
  wx = wx - norm*ux
  wy = wy - norm*uy
  wz = wz - norm*uz

  call InnerProd(ex,norm,wx,wy,wz,vx,vy,vz,gxx,gxy,gxz,gyy,gyz,gzz)
  wx = wx - norm*vx
  wy = wy - norm*vy
  wz = wz - norm*vz

  call InnerProd(ex,norm,wx,wy,wz,wx,wy,wz,gxx,gxy,gxz,gyy,gyz,gzz)
  wx = wx/dsqrt(norm)
  wy = wy/dsqrt(norm)
  wz = wz/dsqrt(norm)

  return

  end subroutine get_triad1
!--------------------------------------------------------------------
!    Gram-Schmidt orthonormal in Cartesin coordinate
!    V1 = [ x, y, z  ]
!	 V2 = [-y, x, ZEO]
!	 V3 = [xz,yz,-(x^2+y^2)]
!	 V1 -> V1 / sqrt(W11)
!	 V2 -> ( V2 - V1 * W12 ) / sqrt(W22)
!	 V3 -> ( V3 - V1 * W13 - V2 * W23 ) / sqrt(W33)
!	 W_ij = g_ab * Vi^a * Vj^b
!	 it is metric, not tilde metric
! raise V1, then V1 first
!--------------------------------------------------------------------

  subroutine get_triad2(ex, X, Y, Z, ep4phi, &
                       gxxi,gxyi,gxzi,gyyi,gyzi,gzzi, &
                       vx,vy,vz,ux,uy,uz,wx,wy,wz)
  
  implicit none

!~~~~~~ argument variables

  integer,intent(in ):: ex(1:3)
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: ep4phi  !exp(4 * phi)
! tilted metric
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: gxxi,gxyi,gxzi,gyyi,gyzi,gzzi
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: vx,vy,vz,ux,uy,uz,wx,wy,wz

!~~~~~~ local variables

  real*8, dimension(ex(1),ex(2),ex(3)) :: norm,fx,fy,fz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupxx,gupxy,gupxz,gupyy,gupyz,gupzz

  real*8,parameter:: ZEO = 0.d0, ONE = 1.d0
  integer::i,j,k

!~~~~~~

  gxx = gxxi * ep4phi
  gxy = gxyi * ep4phi
  gxz = gxzi * ep4phi
  gyy = gyyi * ep4phi
  gyz = gyzi * ep4phi
  gzz = gzzi * ep4phi
! invert metric
  gupzz =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz - &
           gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz
  gupxx =   ( gyy * gzz - gyz * gyz ) / gupzz
  gupxy = - ( gxy * gzz - gyz * gxz ) / gupzz
  gupxz =   ( gxy * gyz - gyy * gxz ) / gupzz
  gupyy =   ( gxx * gzz - gxz * gxz ) / gupzz
  gupyz = - ( gxx * gyz - gxy * gxz ) / gupzz
  gupzz =   ( gxx * gyy - gxy * gxy ) / gupzz
! initialize U, V, W vetors	
  do i=1,ex(1)
   do j=1,ex(2)
    do k=1,ex(3)

     vx(i,j,k) = X(i)
     vy(i,j,k) = Y(j)
     vz(i,j,k) = Z(k)
     ux(i,j,k) = - Y(j)
     uy(i,j,k) = X(i)
     uz(i,j,k) = ZEO
     wx(i,j,k) = X(i)*Z(k)
     wy(i,j,k) = Y(j)*Z(k)
     wz(i,j,k) = -(X(i)*X(i) + Y(j)*Y(j))

    enddo
   enddo
  enddo

  fx = vx
  fy = vy
  fz = vz
  call raise(ex,fx,fy,fz,vx,vy,vz, &
             gupxx,gupxy,gupxz,gupyy,gupyz,gupzz)
! Gram-Schmidt orthonormalization
  call InnerProd(ex,norm,vx,vy,vz,vx,vy,vz,gxx,gxy,gxz,gyy,gyz,gzz)
  vx = vx/dsqrt(norm)
  vy = vy/dsqrt(norm)
  vz = vz/dsqrt(norm)

  call InnerProd(ex,norm,ux,uy,uz,vx,vy,vz,gxx,gxy,gxz,gyy,gyz,gzz)
  ux = ux - norm*vx
  uy = uy - norm*vy
  uz = uz - norm*vz

  call InnerProd(ex,norm,ux,uy,uz,ux,uy,uz,gxx,gxy,gxz,gyy,gyz,gzz)
  ux = ux/dsqrt(norm)
  uy = uy/dsqrt(norm)
  uz = uz/dsqrt(norm)

  call InnerProd(ex,norm,wx,wy,wz,vx,vy,vz,gxx,gxy,gxz,gyy,gyz,gzz)
  wx = wx - norm*vx
  wy = wy - norm*vy
  wz = wz - norm*vz

  call InnerProd(ex,norm,wx,wy,wz,ux,uy,uz,gxx,gxy,gxz,gyy,gyz,gzz)
  wx = wx - norm*ux
  wy = wy - norm*uy
  wz = wz - norm*uz

  call InnerProd(ex,norm,wx,wy,wz,wx,wy,wz,gxx,gxy,gxz,gyy,gyz,gzz)
  wx = wx/dsqrt(norm)
  wy = wy/dsqrt(norm)
  wz = wz/dsqrt(norm)

  return

  end subroutine get_triad2
!***********for shell*********************
!--------------------------------------------------------------------
!    Gram-Schmidt orthonormal in Cartesin coordinate
!    V1 = [ x, y, z  ]
!	 V2 = [-y, x, ZEO]
!	 V3 = [xz,yz,-(x^2+y^2)]
!	 V1 -> V1 / sqrt(W11)
!	 V2 -> ( V2 - V1 * W12 ) / sqrt(W22)
!	 V3 -> ( V3 - V1 * W13 - V2 * W23 ) / sqrt(W33)
!	 W_ij = g_ab * Vi^a * Vj^b
!	 it is metric, not tilde metric
! V1 first
!--------------------------------------------------------------------

  subroutine get_triad0_ss(ex, X, Y, Z, ep4phi, &
                       gxxi,gxyi,gxzi,gyyi,gyzi,gzzi, &
                       vx,vy,vz,ux,uy,uz,wx,wy,wz)
  
  implicit none

!~~~~~~ argument variables

  integer,intent(in ):: ex(1:3)
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: ep4phi  !exp(4 * phi)
! tilted metric
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: gxxi,gxyi,gxzi,gyyi,gyzi,gzzi
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: vx,vy,vz,ux,uy,uz,wx,wy,wz

!~~~~~~ local variables

  real*8, dimension(ex(1),ex(2),ex(3)) :: norm
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gxy,gxz,gyy,gyz,gzz

  real*8,parameter:: ZEO = 0.d0, ONE = 1.d0
  integer::i,j,k
  real*8,parameter::TINYRR=1.d-14

!~~~~~~

  gxx = gxxi * ep4phi
  gxy = gxyi * ep4phi
  gxz = gxzi * ep4phi
  gyy = gyyi * ep4phi
  gyz = gyzi * ep4phi
  gzz = gzzi * ep4phi

! initialize U, V, W vetors	
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

! Gram-Schmidt orthonormalization
  call InnerProd(ex,norm,vx,vy,vz,vx,vy,vz,gxx,gxy,gxz,gyy,gyz,gzz)
  vx = vx/dsqrt(norm)
  vy = vy/dsqrt(norm)
  vz = vz/dsqrt(norm)

  call InnerProd(ex,norm,ux,uy,uz,vx,vy,vz,gxx,gxy,gxz,gyy,gyz,gzz)
  ux = ux - norm*vx
  uy = uy - norm*vy
  uz = uz - norm*vz

  call InnerProd(ex,norm,ux,uy,uz,ux,uy,uz,gxx,gxy,gxz,gyy,gyz,gzz)
  ux = ux/dsqrt(norm)
  uy = uy/dsqrt(norm)
  uz = uz/dsqrt(norm)

  call InnerProd(ex,norm,wx,wy,wz,vx,vy,vz,gxx,gxy,gxz,gyy,gyz,gzz)
  wx = wx - norm*vx
  wy = wy - norm*vy
  wz = wz - norm*vz

  call InnerProd(ex,norm,wx,wy,wz,ux,uy,uz,gxx,gxy,gxz,gyy,gyz,gzz)
  wx = wx - norm*ux
  wy = wy - norm*uy
  wz = wz - norm*uz

  call InnerProd(ex,norm,wx,wy,wz,wx,wy,wz,gxx,gxy,gxz,gyy,gyz,gzz)
  wx = wx/dsqrt(norm)
  wy = wy/dsqrt(norm)
  wz = wz/dsqrt(norm)

  return

  end subroutine get_triad0_ss
!--------------------------------------------------------------------
!    Gram-Schmidt orthonormal in Cartesin coordinate
!    V1 = [ x, y, z  ]
!	 V2 = [-y, x, ZEO]
!	 V3 = [xz,yz,-(x^2+y^2)]
!	 V1 -> V1 / sqrt(W11)
!	 V2 -> ( V2 - V1 * W12 ) / sqrt(W22)
!	 V3 -> ( V3 - V1 * W13 - V2 * W23 ) / sqrt(W33)
!	 W_ij = g_ab * Vi^a * Vj^b
!	 it is metric, not tilde metric
! V2 first
!--------------------------------------------------------------------

  subroutine get_triad1_ss(ex, X, Y, Z, ep4phi, &
                       gxxi,gxyi,gxzi,gyyi,gyzi,gzzi, &
                       vx,vy,vz,ux,uy,uz,wx,wy,wz)
  
  implicit none

!~~~~~~ argument variables

  integer,intent(in ):: ex(1:3)
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: ep4phi  !exp(4 * phi)
! tilted metric
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: gxxi,gxyi,gxzi,gyyi,gyzi,gzzi
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: vx,vy,vz,ux,uy,uz,wx,wy,wz

!~~~~~~ local variables

  real*8, dimension(ex(1),ex(2),ex(3)) :: norm
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gxy,gxz,gyy,gyz,gzz

  real*8,parameter:: ZEO = 0.d0, ONE = 1.d0
  integer::i,j,k
  real*8,parameter::TINYRR=1.d-14

!~~~~~~

  gxx = gxxi * ep4phi
  gxy = gxyi * ep4phi
  gxz = gxzi * ep4phi
  gyy = gyyi * ep4phi
  gyz = gyzi * ep4phi
  gzz = gzzi * ep4phi

! initialize U, V, W vetors	
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

! Gram-Schmidt orthonormalization
  call InnerProd(ex,norm,ux,uy,uz,ux,uy,uz,gxx,gxy,gxz,gyy,gyz,gzz)
  ux = ux/dsqrt(norm)
  uy = uy/dsqrt(norm)
  uz = uz/dsqrt(norm)

  call InnerProd(ex,norm,vx,vy,vz,ux,uy,uz,gxx,gxy,gxz,gyy,gyz,gzz)
  vx = vx - norm*ux
  vy = vy - norm*uy
  vz = vz - norm*uz

  call InnerProd(ex,norm,vx,vy,vz,vx,vy,vz,gxx,gxy,gxz,gyy,gyz,gzz)
  vx = vx/dsqrt(norm)
  vy = vy/dsqrt(norm)
  vz = vz/dsqrt(norm)

  call InnerProd(ex,norm,wx,wy,wz,ux,uy,uz,gxx,gxy,gxz,gyy,gyz,gzz)
  wx = wx - norm*ux
  wy = wy - norm*uy
  wz = wz - norm*uz

  call InnerProd(ex,norm,wx,wy,wz,vx,vy,vz,gxx,gxy,gxz,gyy,gyz,gzz)
  wx = wx - norm*vx
  wy = wy - norm*vy
  wz = wz - norm*vz

  call InnerProd(ex,norm,wx,wy,wz,wx,wy,wz,gxx,gxy,gxz,gyy,gyz,gzz)
  wx = wx/dsqrt(norm)
  wy = wy/dsqrt(norm)
  wz = wz/dsqrt(norm)

  return

  end subroutine get_triad1_ss
!--------------------------------------------------------------------
!    Gram-Schmidt orthonormal in Cartesin coordinate
!    V1 = [ x, y, z  ]
!	 V2 = [-y, x, ZEO]
!	 V3 = [xz,yz,-(x^2+y^2)]
!	 V1 -> V1 / sqrt(W11)
!	 V2 -> ( V2 - V1 * W12 ) / sqrt(W22)
!	 V3 -> ( V3 - V1 * W13 - V2 * W23 ) / sqrt(W33)
!	 W_ij = g_ab * Vi^a * Vj^b
!	 it is metric, not tilde metric
! raise V1, then V1 first
!--------------------------------------------------------------------

  subroutine get_triad2_ss(ex, X, Y, Z, ep4phi, &
                       gxxi,gxyi,gxzi,gyyi,gyzi,gzzi, &
                       vx,vy,vz,ux,uy,uz,wx,wy,wz)
  
  implicit none

!~~~~~~ argument variables

  integer,intent(in ):: ex(1:3)
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: X,Y,Z
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: ep4phi  !exp(4 * phi)
! tilted metric
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: gxxi,gxyi,gxzi,gyyi,gyzi,gzzi
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: vx,vy,vz,ux,uy,uz,wx,wy,wz

!~~~~~~ local variables

  real*8, dimension(ex(1),ex(2),ex(3)) :: norm,fx,fy,fz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupxx,gupxy,gupxz,gupyy,gupyz,gupzz

  real*8,parameter:: ZEO = 0.d0, ONE = 1.d0
  integer::i,j,k
  real*8,parameter::TINYRR=1.d-14

!~~~~~~

  gxx = gxxi * ep4phi
  gxy = gxyi * ep4phi
  gxz = gxzi * ep4phi
  gyy = gyyi * ep4phi
  gyz = gyzi * ep4phi
  gzz = gzzi * ep4phi
! invert metric
  gupzz =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz - &
           gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz
  gupxx =   ( gyy * gzz - gyz * gyz ) / gupzz
  gupxy = - ( gxy * gzz - gyz * gxz ) / gupzz
  gupxz =   ( gxy * gyz - gyy * gxz ) / gupzz
  gupyy =   ( gxx * gzz - gxz * gxz ) / gupzz
  gupyz = - ( gxx * gyz - gxy * gxz ) / gupzz
  gupzz =   ( gxx * gyy - gxy * gxy ) / gupzz
! initialize U, V, W vetors	
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

  fx = vx
  fy = vy
  fz = vz
  call raise(ex,fx,fy,fz,vx,vy,vz, &
             gupxx,gupxy,gupxz,gupyy,gupyz,gupzz)
! Gram-Schmidt orthonormalization
  call InnerProd(ex,norm,vx,vy,vz,vx,vy,vz,gxx,gxy,gxz,gyy,gyz,gzz)
  vx = vx/dsqrt(norm)
  vy = vy/dsqrt(norm)
  vz = vz/dsqrt(norm)

  call InnerProd(ex,norm,ux,uy,uz,vx,vy,vz,gxx,gxy,gxz,gyy,gyz,gzz)
  ux = ux - norm*vx
  uy = uy - norm*vy
  uz = uz - norm*vz

  call InnerProd(ex,norm,ux,uy,uz,ux,uy,uz,gxx,gxy,gxz,gyy,gyz,gzz)
  ux = ux/dsqrt(norm)
  uy = uy/dsqrt(norm)
  uz = uz/dsqrt(norm)

  call InnerProd(ex,norm,wx,wy,wz,vx,vy,vz,gxx,gxy,gxz,gyy,gyz,gzz)
  wx = wx - norm*vx
  wy = wy - norm*vy
  wz = wz - norm*vz

  call InnerProd(ex,norm,wx,wy,wz,ux,uy,uz,gxx,gxy,gxz,gyy,gyz,gzz)
  wx = wx - norm*ux
  wy = wy - norm*uy
  wz = wz - norm*uz

  call InnerProd(ex,norm,wx,wy,wz,wx,wy,wz,gxx,gxy,gxz,gyy,gyz,gzz)
  wx = wx/dsqrt(norm)
  wy = wy/dsqrt(norm)
  wz = wz/dsqrt(norm)

  return

  end subroutine get_triad2_ss
