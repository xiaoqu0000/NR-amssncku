

#include "macrodef.fh"

!-----------------------------------------------------------------------------
!
! compute 4 dimensional Ricci scalar
! this routine is valid for both box and shell
!
!-----------------------------------------------------------------------------

  subroutine get4ricciscalar(ex, X, Y, Z, &
               chi, trK, rho, &
               dxx,gxy,gxz,dyy,gyz,dzz, &
               Axx,Axy,Axz,Ayy,Ayz,Azz, &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,&
               Sxx,Sxy,Sxz,Syy,Syz,Szz,&
               RR)

  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3)
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: dxx,gxy,gxz,dyy,gyz,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: chi,trK,rho
! physical Ricci tensor  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: Rxx,Rxy,Rxz,Ryy,Ryz,Rzz
! matter 
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: Sxx,Sxy,Sxz,Syy,Syz,Szz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out):: RR

!~~~~~~> Other variables:

  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gyy,gzz,chipn1
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupyy,gupyz,gupzz
  real*8, parameter :: ONE = 1.d0, TWO = 2.d0, THR = 3.d0, F8 = 8.d0, F2o3 = 2.d0/3.d0
  real*8 :: PI

  PI = dacos(-ONE)

  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE
  chipn1= chi + ONE

! invert tilted metric
  gupzz =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz - &
           gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz
  gupxx =   ( gyy * gzz - gyz * gyz ) / gupzz
  gupxy = - ( gxy * gzz - gyz * gxz ) / gupzz
  gupxz =   ( gxy * gyz - gyy * gxz ) / gupzz
  gupyy =   ( gxx * gzz - gxz * gxz ) / gupzz
  gupyz = - ( gxx * gyz - gxy * gxz ) / gupzz
  gupzz =   ( gxx * gyy - gxy * gxy ) / gupzz

  RR =(gupxx * ( &
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
       gupyz * (Ayy * Azz + Ayz * Ayz) ) )) - F2o3*trK*trK &
       -(gupxx*Rxx+gupyy*Ryy+gupzz*Rzz+TWO*(gupxy*Rxy+gupxz*Rxz+gupyz*Ryz))*chipn1 &
       -F8*PI*(THR*rho- &
       (gupxx*Sxx+gupyy*Syy+gupzz*Szz+TWO*(gupxy*Sxy+gupxz*Sxz+gupyz*Syz))*chipn1)

  return

  end subroutine get4ricciscalar  
