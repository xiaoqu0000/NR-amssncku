

#include "macrodef.fh"

!-----------------------------------------------------------------------------
!
! compute the Newman-Penrose Weyl scalar Psi4
! for BSSN dynamical variables
!
!-----------------------------------------------------------------------------

  subroutine getnp4scalar(ex, X, Y, Z,                            &
               chi, trK, Sphi,&
               dxx,gxy,gxz,dyy,gyz,dzz, &
               Axx,Axy,Axz,Ayy,Ayz,Azz, &
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,&
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,&
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,&
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,&
               Rpsi4, Ipsi4, &
               symmetry)

  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3),symmetry
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: dxx,gxy,gxz,dyy,gyz,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: chi,trK,Sphi
! physical second kind of connection  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Gamzyy, Gamzyz, Gamzzz
! physical Ricci tensor  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Rxx,Rxy,Rxz,Ryy,Ryz,Rzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out):: Rpsi4,Ipsi4

!~~~~~~> Other variables:

  real*8, dimension(ex(1),ex(2),ex(3)) :: f,fx,fy,fz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gyy,gzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: chix,chiy,chiz,chipn1
  real*8, dimension(ex(1),ex(2),ex(3)) :: vx,vy,vz,ux,uy,uz,wx,wy,wz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Exx,Exy,Exz,Eyy,Eyz,Ezz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Bxx,Bxy,Bxz,Byy,Byz,Bzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Axxx,Axxy,Axxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Axyx,Axyy,Axyz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Axzx,Axzy,Axzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Ayyx,Ayyy,Ayyz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Ayzx,Ayzy,Ayzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Azzx,Azzy,Azzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupyy,gupyz,gupzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: uuwwxx,uuwwxy,uuwwxz,uuwwyy,uuwwyz,uuwwzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: uwxx,uwxy,uwxz,uwyy,uwyz,uwzz

  real*8, parameter :: ZEO = 0.d0, ONE = 1.d0, TWO = 2.d0
  real*8, parameter :: F1o3 = 1.d0/3.d0, FOUR = 4.d0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  real*8            :: dX, dY, dZ
  integer::i,j,k
  real*8,parameter::TINYRR=1.d-14
  real*8 :: PI

  PI = dacos(-ONE)

  call getnp4(ex, X, Y, Z,                            &
               chi, trK, &
               dxx,gxy,gxz,dyy,gyz,dzz, &
               Axx,Axy,Axz,Ayy,Ayz,Azz, &
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,&
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,&
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,&
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,&
               Rpsi4, Ipsi4, &
               symmetry)

  Rpsi4 = dexp(-FOUR*dsqrt(PI/3)*Sphi)*Rpsi4
  Ipsi4 = dexp(-FOUR*dsqrt(PI/3)*Sphi)*Ipsi4

  return

  end subroutine getnp4scalar
! 4D method  
  subroutine getnp4oldscalar(ex, X, Y, Z,  chi,  trK,Sphi,        &
                      dxx,  gxy,  gxz,  dyy,  gyz,  dzz, &
                      Axx,  Axy,  Axz,  Ayy,  Ayz,  Azz, &
                      Gmx,  Gmy,  Gmz,                   & 
                      Lap,  Sfx,  Sfy,  Sfz,Rpsi4,Ipsi4, symmetry)

  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3),symmetry
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: chi,trK,Sphi
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: dxx,gxy,gxz,dyy,gyz,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gmx,Gmy,Gmz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Lap,Sfx,Sfy,Sfz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Rpsi4,Ipsi4

  real*8 :: PI

  real*8, parameter :: ZEO = 0.d0, ONE = 1.d0, TWO = 2.d0
  real*8, parameter :: F1o3 = 1.d0/3.d0, FOUR = 4.d0

  PI = dacos(-ONE)

  call getnp4old(ex, X, Y, Z,  chi,  trK,             &
                      dxx,  gxy,  gxz,  dyy,  gyz,  dzz, &
                      Axx,  Axy,  Axz,  Ayy,  Ayz,  Azz, &
                      Gmx,  Gmy,  Gmz,                   & 
                      Lap,  Sfx,  Sfy,  Sfz,Rpsi4,Ipsi4, symmetry)

  Rpsi4 = dexp(-FOUR*dsqrt(PI/3)*Sphi)*Rpsi4
  Ipsi4 = dexp(-FOUR*dsqrt(PI/3)*Sphi)*Ipsi4

  return

  end subroutine getnp4oldscalar
!-----------------------------------------------------------------------------
! for shell  
!-----------------------------------------------------------------------------

  subroutine getnp4scalar_ss(ex,crho,sigma,R, X, Y, Z,                               &
               drhodx, drhody, drhodz,                                         &
               dsigmadx,dsigmady,dsigmadz,                                     &
               dRdx,dRdy,dRdz,                                                 &
               drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                &
               dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,    &
               dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,                            &
               chi, trK, Sphi,&
               dxx,gxy,gxz,dyy,gyz,dzz, &
               Axx,Axy,Axz,Ayy,Ayz,Azz, &
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,&
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,&
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,&
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,&
               Rpsi4, Ipsi4, &
               symmetry,sst)

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
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: dxx,gxy,gxz,dyy,gyz,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: chi,trK,Sphi
! physical second kind of connection  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Gamzyy, Gamzyz, Gamzzz
! physical Ricci tensor  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Rxx,Rxy,Rxz,Ryy,Ryz,Rzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out):: Rpsi4,Ipsi4

!~~~~~~> Other variables:

  real*8, dimension(ex(1),ex(2),ex(3)) :: f,fx,fy,fz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gyy,gzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: chix,chiy,chiz,chipn1
  real*8, dimension(ex(1),ex(2),ex(3)) :: vx,vy,vz,ux,uy,uz,wx,wy,wz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Exx,Exy,Exz,Eyy,Eyz,Ezz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Bxx,Bxy,Bxz,Byy,Byz,Bzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Axxx,Axxy,Axxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Axyx,Axyy,Axyz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Axzx,Axzy,Axzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Ayyx,Ayyy,Ayyz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Ayzx,Ayzy,Ayzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Azzx,Azzy,Azzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupyy,gupyz,gupzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: uuwwxx,uuwwxy,uuwwxz,uuwwyy,uuwwyz,uuwwzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: uwxx,uwxy,uwxz,uwyy,uwyz,uwzz

  real*8, parameter :: ZEO = 0.d0, ONE = 1.d0, TWO = 2.d0
  real*8, parameter :: F1o3 = 1.d0/3.d0, FOUR = 4.d0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  integer::i,j,k
  real*8,parameter::TINYRR=1.d-14
  real*8 :: PI

  PI = dacos(-ONE)

  call getnp4_ss(ex,crho,sigma,R, X, Y, Z,                               &
               drhodx, drhody, drhodz,                                         &
               dsigmadx,dsigmady,dsigmadz,                                     &
               dRdx,dRdy,dRdz,                                                 &
               drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                &
               dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,    &
               dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,                            &
               chi, trK, &
               dxx,gxy,gxz,dyy,gyz,dzz, &
               Axx,Axy,Axz,Ayy,Ayz,Azz, &
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,&
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,&
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,&
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,&
               Rpsi4, Ipsi4, &
               symmetry,sst)

  Rpsi4 = dexp(-FOUR*dsqrt(PI/3)*Sphi)*Rpsi4
  Ipsi4 = dexp(-FOUR*dsqrt(PI/3)*Sphi)*Ipsi4

  return

  end subroutine getnp4scalar_ss
! 4D method  
  subroutine getnp4oldscalar_ss(ex,crho,sigma,R, X, Y, Z,                               &
               drhodx, drhody, drhodz,                                         &
               dsigmadx,dsigmady,dsigmadz,                                     &
               dRdx,dRdy,dRdz,                                                 &
               drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                &
               dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,    &
               dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,                            &
                      chi,  trK, Sphi,             &
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
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: chi,trK,Sphi
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: dxx,gxy,gxz,dyy,gyz,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gmx,Gmy,Gmz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Lap,Sfx,Sfy,Sfz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Rpsi4,Ipsi4

  real*8 :: PI

  real*8, parameter :: ZEO = 0.d0, ONE = 1.d0, TWO = 2.d0
  real*8, parameter :: F1o3 = 1.d0/3.d0, FOUR = 4.d0

  PI = dacos(-ONE)

  call getnp4old_ss(ex,crho,sigma,R, X, Y, Z,                               &
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

  Rpsi4 = dexp(-FOUR*dsqrt(PI/3)*Sphi)*Rpsi4
  Ipsi4 = dexp(-FOUR*dsqrt(PI/3)*Sphi)*Ipsi4

  return

  end subroutine getnp4oldscalar_ss
