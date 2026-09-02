

#include "macrodef.fh"

#if (ABV == 0)  
!! using BSSN variables
!-------------------------------------------------------------------------------!
! computed constraint for bssn formalism                                        !
!-------------------------------------------------------------------------------!
  subroutine constraint_bssn(ex, X, Y, Z,&
               chi,trK, &
               dxx,gxy,gxz,dyy,gyz,dzz, &
               Axx,Axy,Axz,Ayy,Ayz,Azz, &
               Gmx,Gmy,Gmz,&
               Lap,Sfx,Sfy,Sfz,rho,Sx,Sy,Sz,&
               Gamxxx, Gamxxy, Gamxxz,Gamxyy, Gamxyz, Gamxzz, &
               Gamyxx, Gamyxy, Gamyxz,Gamyyy, Gamyyz, Gamyzz, &
               Gamzxx, Gamzxy, Gamzxz,Gamzyy, Gamzyz, Gamzzz, &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz, &
               ham_Res, movx_Res, movy_Res, movz_Res, Gmx_Res, Gmy_Res, Gmz_Res, &
               Symmetry)

  implicit none
!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3),symmetry
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: chi,trK
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: dxx,gxy,gxz,dyy,gyz,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gmx,Gmy,Gmz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Lap,Sfx,Sfy,Sfz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: rho,Sx,Sy,Sz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Rxx,Rxy,Rxz,Ryy,Ryz,Rzz
! second kind of Christofel symble Gamma^i_jk respect to physical metric
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamzyy, Gamzyz, Gamzzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: ham_Res, movx_Res, movy_Res, movz_Res
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gmx_Res, Gmy_Res, Gmz_Res
!~~~~~~> Other variables:
!  inverse metric
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupyy,gupyz,gupzz
! first order derivative of metric, @_k g_ij
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxx,gxyx,gxzx
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyyx,gyzx,gzzx
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxy,gxyy,gxzy
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyyy,gyzy,gzzy
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxz,gxyz,gxzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyyz,gyzz,gzzz
!  partial derivative of chi, chi_i
  real*8, dimension(ex(1),ex(2),ex(3)) :: chin1,chix,chiy,chiz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gyy,gzz

  integer, parameter :: NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2
  real*8, parameter :: ZERO = 0.D0, HALF = 0.5d0, ONE = 1.d0, TWO = 2.d0, FOUR = 4.d0
  real*8, parameter :: F2o3 = 2.d0/3.d0, F8 = 8.d0, F16 = 1.6d1, SIX = 6.d0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  real*8            :: PI
  integer :: i, j, k

  PI = dacos(-ONE)

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
  gxx(i,j,k) = dxx(i,j,k) + ONE
  gyy(i,j,k) = dyy(i,j,k) + ONE
  gzz(i,j,k) = dzz(i,j,k) + ONE
  chin1(i,j,k) = chi(i,j,k)+ONE
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO
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

  call fderivs(ex,dxx,gxxx,gxxy,gxxz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)
  call fderivs(ex,gxy,gxyx,gxyy,gxyz,X,Y,Z,ANTI,ANTI,SYM ,Symmetry,0)
  call fderivs(ex,gxz,gxzx,gxzy,gxzz,X,Y,Z,ANTI,SYM ,ANTI,Symmetry,0)
  call fderivs(ex,dyy,gyyx,gyyy,gyyz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)
  call fderivs(ex,gyz,gyzx,gyzy,gyzz,X,Y,Z,SYM ,ANTI,ANTI,Symmetry,0)
  call fderivs(ex,dzz,gzzx,gzzy,gzzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)

! Gam^i_Res = Gam^i + gup^ij_,j
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
  Gmx_Res(i,j,k) = Gmx(i,j,k) - (gupxx(i,j,k)*(gupxx(i,j,k)*gxxx(i,j,k)+gupxy(i,j,k)*gxyx(i,j,k)+gupxz(i,j,k)*gxzx(i,j,k))&
                  +gupxy(i,j,k)*(gupxx(i,j,k)*gxyx(i,j,k)+gupxy(i,j,k)*gyyx(i,j,k)+gupxz(i,j,k)*gyzx(i,j,k))&
                  +gupxz(i,j,k)*(gupxx(i,j,k)*gxzx(i,j,k)+gupxy(i,j,k)*gyzx(i,j,k)+gupxz(i,j,k)*gzzx(i,j,k))&
                  +gupxx(i,j,k)*(gupxy(i,j,k)*gxxy(i,j,k)+gupyy(i,j,k)*gxyy(i,j,k)+gupyz(i,j,k)*gxzy(i,j,k))&
                  +gupxy(i,j,k)*(gupxy(i,j,k)*gxyy(i,j,k)+gupyy(i,j,k)*gyyy(i,j,k)+gupyz(i,j,k)*gyzy(i,j,k))&
                  +gupxz(i,j,k)*(gupxy(i,j,k)*gxzy(i,j,k)+gupyy(i,j,k)*gyzy(i,j,k)+gupyz(i,j,k)*gzzy(i,j,k))&
                  +gupxx(i,j,k)*(gupxz(i,j,k)*gxxz(i,j,k)+gupyz(i,j,k)*gxyz(i,j,k)+gupzz(i,j,k)*gxzz(i,j,k))&
                  +gupxy(i,j,k)*(gupxz(i,j,k)*gxyz(i,j,k)+gupyz(i,j,k)*gyyz(i,j,k)+gupzz(i,j,k)*gyzz(i,j,k))&
                  +gupxz(i,j,k)*(gupxz(i,j,k)*gxzz(i,j,k)+gupyz(i,j,k)*gyzz(i,j,k)+gupzz(i,j,k)*gzzz(i,j,k)))
  Gmy_Res(i,j,k) = Gmy(i,j,k) - (gupxx(i,j,k)*(gupxy(i,j,k)*gxxx(i,j,k)+gupyy(i,j,k)*gxyx(i,j,k)+gupyz(i,j,k)*gxzx(i,j,k))&
                  +gupxy(i,j,k)*(gupxy(i,j,k)*gxyx(i,j,k)+gupyy(i,j,k)*gyyx(i,j,k)+gupyz(i,j,k)*gyzx(i,j,k))&
                  +gupxz(i,j,k)*(gupxy(i,j,k)*gxzx(i,j,k)+gupyy(i,j,k)*gyzx(i,j,k)+gupyz(i,j,k)*gzzx(i,j,k))&
                  +gupxy(i,j,k)*(gupxy(i,j,k)*gxxy(i,j,k)+gupyy(i,j,k)*gxyy(i,j,k)+gupyz(i,j,k)*gxzy(i,j,k))&
                  +gupyy(i,j,k)*(gupxy(i,j,k)*gxyy(i,j,k)+gupyy(i,j,k)*gyyy(i,j,k)+gupyz(i,j,k)*gyzy(i,j,k))&
                  +gupyz(i,j,k)*(gupxy(i,j,k)*gxzy(i,j,k)+gupyy(i,j,k)*gyzy(i,j,k)+gupyz(i,j,k)*gzzy(i,j,k))&
                  +gupxy(i,j,k)*(gupxz(i,j,k)*gxxz(i,j,k)+gupyz(i,j,k)*gxyz(i,j,k)+gupzz(i,j,k)*gxzz(i,j,k))&
                  +gupyy(i,j,k)*(gupxz(i,j,k)*gxyz(i,j,k)+gupyz(i,j,k)*gyyz(i,j,k)+gupzz(i,j,k)*gyzz(i,j,k))&
                  +gupyz(i,j,k)*(gupxz(i,j,k)*gxzz(i,j,k)+gupyz(i,j,k)*gyzz(i,j,k)+gupzz(i,j,k)*gzzz(i,j,k)))
  Gmz_Res(i,j,k) = Gmz(i,j,k) - (gupxx(i,j,k)*(gupxz(i,j,k)*gxxx(i,j,k)+gupyz(i,j,k)*gxyx(i,j,k)+gupzz(i,j,k)*gxzx(i,j,k))&
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

! M_j = gupki*(-1/chi d_k chi*A_ij + D_k A_ij) - 2/3 d_j trK - 8 PI s_j where D respect to physical metric
! store D_i A_jk - 1/chi d_i chi*A_jk in gjk_i
  call fderivs(ex,chi,chix,chiy,chiz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)
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
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

!store K,i in chi,i
  call fderivs(ex,trK,chix,chiy,chiz,X,Y,Z,SYM,SYM,SYM,Symmetry,0)

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
movx_Res(i,j,k) = movx_Res(i,j,k) - F2o3*chix(i,j,k) - F8*PI*sx(i,j,k)
movy_Res(i,j,k) = movy_Res(i,j,k) - F2o3*chiy(i,j,k) - F8*PI*sy(i,j,k)
movz_Res(i,j,k) = movz_Res(i,j,k) - F2o3*chiz(i,j,k) - F8*PI*sz(i,j,k)
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

  return

  end subroutine constraint_bssn
!-------------------------------------------------------------------------------!
! computed constraint for bssn formalism for shell                              !
!-------------------------------------------------------------------------------!
  subroutine constraint_bssn_ss(ex,crho,sigma,R,X, Y, Z,  &
               drhodx, drhody, drhodz,                                         &
               dsigmadx,dsigmady,dsigmadz,                                     &
               dRdx,dRdy,dRdz,                                                 &
               drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                &
               dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,    &
               dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,                            &
               chi,trK, &
               dxx,gxy,gxz,dyy,gyz,dzz, &
               Axx,Axy,Axz,Ayy,Ayz,Azz, &
               Gmx,Gmy,Gmz,&
               Lap,Sfx,Sfy,Sfz,rho,Sx,Sy,Sz,&
               Gamxxx, Gamxxy, Gamxxz,Gamxyy, Gamxyz, Gamxzz, &
               Gamyxx, Gamyxy, Gamyxz,Gamyyy, Gamyyz, Gamyzz, &
               Gamzxx, Gamzxy, Gamzxz,Gamzyy, Gamzyz, Gamzzz, &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz, &
               ham_Res, movx_Res, movy_Res, movz_Res, Gmx_Res, Gmy_Res, Gmz_Res, &
               Symmetry,Lev,sst)

  implicit none
!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3),symmetry,Lev,sst
  double precision,intent(in),dimension(ex(1))::crho
  double precision,intent(in),dimension(ex(2))::sigma
  double precision,intent(in),dimension(ex(3))::R
  real*8, intent(in ),dimension(ex(1),ex(2),ex(3)):: X,Y,Z
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
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: rho,Sx,Sy,Sz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Rxx,Rxy,Rxz,Ryy,Ryz,Rzz
! second kind of Christofel symble Gamma^i_jk respect to physical metric
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamzyy, Gamzyz, Gamzzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: ham_Res, movx_Res, movy_Res, movz_Res
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gmx_Res, Gmy_Res, Gmz_Res
!~~~~~~> Other variables:
!  inverse metric
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupyy,gupyz,gupzz
! first order derivative of metric, @_k g_ij
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxx,gxyx,gxzx
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyyx,gyzx,gzzx
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxy,gxyy,gxzy
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyyy,gyzy,gzzy
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxz,gxyz,gxzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyyz,gyzz,gzzz
!  partial derivative of chi, chi_i
  real*8, dimension(ex(1),ex(2),ex(3)) :: chin1,chix,chiy,chiz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gyy,gzz

  integer, parameter :: NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2
  real*8, parameter :: ZERO = 0.D0, HALF = 0.5d0, ONE = 1.d0, TWO = 2.d0, FOUR = 4.d0
  real*8, parameter :: F2o3 = 2.d0/3.d0, F8 = 8.d0, F16 = 1.6d1, SIX = 6.d0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  real*8            :: PI
  integer :: i, j, k

  PI = dacos(-ONE)

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
  gxx(i,j,k) = dxx(i,j,k) + ONE
  gyy(i,j,k) = dyy(i,j,k) + ONE
  gzz(i,j,k) = dzz(i,j,k) + ONE
  chin1(i,j,k) = chi(i,j,k)+ONE
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO
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

! Gam^i_Res = Gam^i + gup^ij_,j
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
  Gmx_Res(i,j,k) = Gmx(i,j,k) - (gupxx(i,j,k)*(gupxx(i,j,k)*gxxx(i,j,k)+gupxy(i,j,k)*gxyx(i,j,k)+gupxz(i,j,k)*gxzx(i,j,k))&
                  +gupxy(i,j,k)*(gupxx(i,j,k)*gxyx(i,j,k)+gupxy(i,j,k)*gyyx(i,j,k)+gupxz(i,j,k)*gyzx(i,j,k))&
                  +gupxz(i,j,k)*(gupxx(i,j,k)*gxzx(i,j,k)+gupxy(i,j,k)*gyzx(i,j,k)+gupxz(i,j,k)*gzzx(i,j,k))&
                  +gupxx(i,j,k)*(gupxy(i,j,k)*gxxy(i,j,k)+gupyy(i,j,k)*gxyy(i,j,k)+gupyz(i,j,k)*gxzy(i,j,k))&
                  +gupxy(i,j,k)*(gupxy(i,j,k)*gxyy(i,j,k)+gupyy(i,j,k)*gyyy(i,j,k)+gupyz(i,j,k)*gyzy(i,j,k))&
                  +gupxz(i,j,k)*(gupxy(i,j,k)*gxzy(i,j,k)+gupyy(i,j,k)*gyzy(i,j,k)+gupyz(i,j,k)*gzzy(i,j,k))&
                  +gupxx(i,j,k)*(gupxz(i,j,k)*gxxz(i,j,k)+gupyz(i,j,k)*gxyz(i,j,k)+gupzz(i,j,k)*gxzz(i,j,k))&
                  +gupxy(i,j,k)*(gupxz(i,j,k)*gxyz(i,j,k)+gupyz(i,j,k)*gyyz(i,j,k)+gupzz(i,j,k)*gyzz(i,j,k))&
                  +gupxz(i,j,k)*(gupxz(i,j,k)*gxzz(i,j,k)+gupyz(i,j,k)*gyzz(i,j,k)+gupzz(i,j,k)*gzzz(i,j,k)))
  Gmy_Res(i,j,k) = Gmy(i,j,k) - (gupxx(i,j,k)*(gupxy(i,j,k)*gxxx(i,j,k)+gupyy(i,j,k)*gxyx(i,j,k)+gupyz(i,j,k)*gxzx(i,j,k))&
                  +gupxy(i,j,k)*(gupxy(i,j,k)*gxyx(i,j,k)+gupyy(i,j,k)*gyyx(i,j,k)+gupyz(i,j,k)*gyzx(i,j,k))&
                  +gupxz(i,j,k)*(gupxy(i,j,k)*gxzx(i,j,k)+gupyy(i,j,k)*gyzx(i,j,k)+gupyz(i,j,k)*gzzx(i,j,k))&
                  +gupxy(i,j,k)*(gupxy(i,j,k)*gxxy(i,j,k)+gupyy(i,j,k)*gxyy(i,j,k)+gupyz(i,j,k)*gxzy(i,j,k))&
                  +gupyy(i,j,k)*(gupxy(i,j,k)*gxyy(i,j,k)+gupyy(i,j,k)*gyyy(i,j,k)+gupyz(i,j,k)*gyzy(i,j,k))&
                  +gupyz(i,j,k)*(gupxy(i,j,k)*gxzy(i,j,k)+gupyy(i,j,k)*gyzy(i,j,k)+gupyz(i,j,k)*gzzy(i,j,k))&
                  +gupxy(i,j,k)*(gupxz(i,j,k)*gxxz(i,j,k)+gupyz(i,j,k)*gxyz(i,j,k)+gupzz(i,j,k)*gxzz(i,j,k))&
                  +gupyy(i,j,k)*(gupxz(i,j,k)*gxyz(i,j,k)+gupyz(i,j,k)*gyyz(i,j,k)+gupzz(i,j,k)*gyzz(i,j,k))&
                  +gupyz(i,j,k)*(gupxz(i,j,k)*gxzz(i,j,k)+gupyz(i,j,k)*gyzz(i,j,k)+gupzz(i,j,k)*gzzz(i,j,k)))
  Gmz_Res(i,j,k) = Gmz(i,j,k) - (gupxx(i,j,k)*(gupxz(i,j,k)*gxxx(i,j,k)+gupyz(i,j,k)*gxyx(i,j,k)+gupzz(i,j,k)*gxzx(i,j,k))&
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

! M_j = gupki*(-1/chi d_k chi*A_ij + D_k A_ij) - 2/3 d_j trK - 8 PI s_j where D respect to physical metric
! store D_i A_jk - 1/chi d_i chi*A_jk in gjk_i
  call fderivs_shc(ex,chi,chix,chiy,chiz,crho,sigma,R, SYM, SYM,SYM,Symmetry,Lev,sst,          &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
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
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

!store K,i in chi,i 
  call fderivs_shc(ex,trK,chix,chiy,chiz,crho,sigma,R, SYM, SYM,SYM,Symmetry,Lev,sst,                &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
movx_Res(i,j,k) = movx_Res(i,j,k) - F2o3*chix(i,j,k) - F8*PI*sx(i,j,k)
movy_Res(i,j,k) = movy_Res(i,j,k) - F2o3*chiy(i,j,k) - F8*PI*sy(i,j,k)
movz_Res(i,j,k) = movz_Res(i,j,k) - F2o3*chiz(i,j,k) - F8*PI*sz(i,j,k)
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

  return

  end subroutine constraint_bssn_ss
#elif (ABV == 1)  
!! using ADM variables
!-------------------------------------------------------------------------------!
! computed constraint for bssn formalism                                        !
!-------------------------------------------------------------------------------!
  subroutine constraint_bssn(ex, X, Y, Z,&
               chi,trK, &
               dxx,gxy,gxz,dyy,gyz,dzz, &
               Axx,Axy,Axz,Ayy,Ayz,Azz, &
               Gmx,Gmy,Gmz,&
               Lap,Sfx,Sfy,Sfz,rho,Sx,Sy,Sz,&
               Gamxxx, Gamxxy, Gamxxz,Gamxyy, Gamxyz, Gamxzz, &
               Gamyxx, Gamyxy, Gamyxz,Gamyyy, Gamyyz, Gamyzz, &
               Gamzxx, Gamzxy, Gamzxz,Gamzyy, Gamzyz, Gamzzz, &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz, &
               ham_Res, movx_Res, movy_Res, movz_Res, Gmx_Res, Gmy_Res, Gmz_Res, &
               Symmetry)

  implicit none
!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3),symmetry
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: chi,trK
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: dxx,gxy,gxz,dyy,gyz,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gmx,Gmy,Gmz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Lap,Sfx,Sfy,Sfz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: rho,Sx,Sy,Sz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Rxx,Rxy,Rxz,Ryy,Ryz,Rzz
! second kind of Christofel symble Gamma^i_jk respect to physical metric
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamzyy, Gamzyz, Gamzzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: ham_Res, movx_Res, movy_Res, movz_Res
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gmx_Res, Gmy_Res, Gmz_Res
!~~~~~~> Other variables:
!  inverse metric
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupyy,gupyz,gupzz
! first order derivative of metric, @_k g_ij
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxx,gxyx,gxzx
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyyx,gyzx,gzzx
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxy,gxyy,gxzy
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyyy,gyzy,gzzy
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxz,gxyz,gxzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyyz,gyzz,gzzz
!  partial derivative of chi, chi_i
  real*8, dimension(ex(1),ex(2),ex(3)) :: chin1
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gyy,gzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: adm_dxx,adm_dyy,adm_dzz,adm_gxy,adm_gxz,adm_gyz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Kxx,Kyy,Kzz,Kxy,Kxz,Kyz

  integer, parameter :: NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2
  real*8, parameter :: ZERO = 0.D0, HALF = 0.5d0, ONE = 1.d0, TWO = 2.d0, FOUR = 4.d0
  real*8, parameter :: F2o3 = 2.d0/3.d0, F8 = 8.d0, F16 = 1.6d1, SIX = 6.d0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  real*8            :: PI
  integer :: i, j, k

  PI = dacos(-ONE)

  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE
  chin1 = chi+ONE
! invert tilted metric
  gupzz =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz - &
           gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz
  gupxx =   ( gyy * gzz - gyz * gyz ) / gupzz
  gupxy = - ( gxy * gzz - gyz * gxz ) / gupzz
  gupxz =   ( gxy * gyz - gyy * gxz ) / gupzz
  gupyy =   ( gxx * gzz - gxz * gxz ) / gupzz
  gupyz = - ( gxx * gyz - gxy * gxz ) / gupzz
  gupzz =   ( gxx * gyy - gxy * gxy ) / gupzz

  call fderivs(ex,dxx,gxxx,gxxy,gxxz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)
  call fderivs(ex,gxy,gxyx,gxyy,gxyz,X,Y,Z,ANTI,ANTI,SYM ,Symmetry,0)
  call fderivs(ex,gxz,gxzx,gxzy,gxzz,X,Y,Z,ANTI,SYM ,ANTI,Symmetry,0)
  call fderivs(ex,dyy,gyyx,gyyy,gyyz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)
  call fderivs(ex,gyz,gyzx,gyzy,gyzz,X,Y,Z,SYM ,ANTI,ANTI,Symmetry,0)
  call fderivs(ex,dzz,gzzx,gzzy,gzzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)

! Gam^i_Res = Gam^i + gup^ij_,j
  Gmx_Res = Gmx - (gupxx*(gupxx*gxxx+gupxy*gxyx+gupxz*gxzx)&
                  +gupxy*(gupxx*gxyx+gupxy*gyyx+gupxz*gyzx)&
                  +gupxz*(gupxx*gxzx+gupxy*gyzx+gupxz*gzzx)&
                  +gupxx*(gupxy*gxxy+gupyy*gxyy+gupyz*gxzy)&
                  +gupxy*(gupxy*gxyy+gupyy*gyyy+gupyz*gyzy)&
                  +gupxz*(gupxy*gxzy+gupyy*gyzy+gupyz*gzzy)&
                  +gupxx*(gupxz*gxxz+gupyz*gxyz+gupzz*gxzz)&
                  +gupxy*(gupxz*gxyz+gupyz*gyyz+gupzz*gyzz)&
                  +gupxz*(gupxz*gxzz+gupyz*gyzz+gupzz*gzzz))
  Gmy_Res = Gmy - (gupxx*(gupxy*gxxx+gupyy*gxyx+gupyz*gxzx)&
                  +gupxy*(gupxy*gxyx+gupyy*gyyx+gupyz*gyzx)&
                  +gupxz*(gupxy*gxzx+gupyy*gyzx+gupyz*gzzx)&
                  +gupxy*(gupxy*gxxy+gupyy*gxyy+gupyz*gxzy)&
                  +gupyy*(gupxy*gxyy+gupyy*gyyy+gupyz*gyzy)&
                  +gupyz*(gupxy*gxzy+gupyy*gyzy+gupyz*gzzy)&
                  +gupxy*(gupxz*gxxz+gupyz*gxyz+gupzz*gxzz)&
                  +gupyy*(gupxz*gxyz+gupyz*gyyz+gupzz*gyzz)&
                  +gupyz*(gupxz*gxzz+gupyz*gyzz+gupzz*gzzz))
  Gmz_Res = Gmz - (gupxx*(gupxz*gxxx+gupyz*gxyx+gupzz*gxzx)&
                  +gupxy*(gupxz*gxyx+gupyz*gyyx+gupzz*gyzx)&
                  +gupxz*(gupxz*gxzx+gupyz*gyzx+gupzz*gzzx)&
                  +gupxy*(gupxz*gxxy+gupyz*gxyy+gupzz*gxzy)&
                  +gupyy*(gupxz*gxyy+gupyz*gyyy+gupzz*gyzy)&
                  +gupyz*(gupxz*gxzy+gupyz*gyzy+gupzz*gzzy)&
                  +gupxz*(gupxz*gxxz+gupyz*gxyz+gupzz*gxzz)&
                  +gupyz*(gupxz*gxyz+gupyz*gyyz+gupzz*gyzz)&
                  +gupzz*(gupxz*gxzz+gupyz*gyzz+gupzz*gzzz))

  call bssn2adm(ex,chin1,trK,gxx,gxy,gxz,gyy,gyz,gzz, &
                             Axx,Axy,Axz,Ayy,Ayz,Azz, &
              adm_dxx,adm_gxy,adm_gxz,adm_dyy,adm_gyz,adm_dzz, &
              Kxx,Kxy,Kxz,Kyy,Kyz,Kzz)    
  adm_dxx = adm_dxx - ONE  
  adm_dyy = adm_dyy - ONE  
  adm_dzz = adm_dzz - ONE  

  call constraint_adm(ex, X, Y, Z,&
               adm_dxx,adm_gxy,adm_gxz,adm_dyy,adm_gyz,adm_dzz, &
               Kxx,Kxy,Kxz,Kyy,Kyz,Kzz, &
               Lap,Sfx,Sfy,Sfz,rho,Sx,Sy,Sz,&
               ham_Res, movx_Res, movy_Res, movz_Res, &
               Symmetry)

  return

  end subroutine constraint_bssn
!-------------------------------------------------------------------------------!
! computed constraint for bssn formalism for shell                              !
!-------------------------------------------------------------------------------!
  subroutine constraint_bssn_ss(ex,crho,sigma,R,X, Y, Z,  &
               drhodx, drhody, drhodz,                                         &
               dsigmadx,dsigmady,dsigmadz,                                     &
               dRdx,dRdy,dRdz,                                                 &
               drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                &
               dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,    &
               dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,                            &
               chi,trK, &
               dxx,gxy,gxz,dyy,gyz,dzz, &
               Axx,Axy,Axz,Ayy,Ayz,Azz, &
               Gmx,Gmy,Gmz,&
               Lap,Sfx,Sfy,Sfz,rho,Sx,Sy,Sz,&
               Gamxxx, Gamxxy, Gamxxz,Gamxyy, Gamxyz, Gamxzz, &
               Gamyxx, Gamyxy, Gamyxz,Gamyyy, Gamyyz, Gamyzz, &
               Gamzxx, Gamzxy, Gamzxz,Gamzyy, Gamzyz, Gamzzz, &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz, &
               ham_Res, movx_Res, movy_Res, movz_Res, Gmx_Res, Gmy_Res, Gmz_Res, &
               Symmetry,Lev,sst)

  implicit none
!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3),symmetry,Lev,sst
  double precision,intent(in),dimension(ex(1))::crho
  double precision,intent(in),dimension(ex(2))::sigma
  double precision,intent(in),dimension(ex(3))::R
  real*8, intent(in ),dimension(ex(1),ex(2),ex(3)):: X,Y,Z
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
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: rho,Sx,Sy,Sz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Rxx,Rxy,Rxz,Ryy,Ryz,Rzz
! second kind of Christofel symble Gamma^i_jk respect to physical metric
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamzyy, Gamzyz, Gamzzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: ham_Res, movx_Res, movy_Res, movz_Res
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gmx_Res, Gmy_Res, Gmz_Res
!~~~~~~> Other variables:
!  inverse metric
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupyy,gupyz,gupzz
! first order derivative of metric, @_k g_ij
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxx,gxyx,gxzx
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyyx,gyzx,gzzx
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxy,gxyy,gxzy
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyyy,gyzy,gzzy
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxz,gxyz,gxzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyyz,gyzz,gzzz
!  partial derivative of chi, chi_i
  real*8, dimension(ex(1),ex(2),ex(3)) :: chin1
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gyy,gzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: adm_dxx,adm_dyy,adm_dzz,adm_gxy,adm_gxz,adm_gyz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Kxx,Kyy,Kzz,Kxy,Kxz,Kyz

  integer, parameter :: NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2
  real*8, parameter :: ZERO = 0.D0, HALF = 0.5d0, ONE = 1.d0, TWO = 2.d0, FOUR = 4.d0
  real*8, parameter :: F2o3 = 2.d0/3.d0, F8 = 8.d0, F16 = 1.6d1, SIX = 6.d0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  real*8            :: PI
  integer :: i, j, k

  PI = dacos(-ONE)

  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE
  chin1 = chi+ONE
! invert tilted metric
  gupzz =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz - &
           gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz
  gupxx =   ( gyy * gzz - gyz * gyz ) / gupzz
  gupxy = - ( gxy * gzz - gyz * gxz ) / gupzz
  gupxz =   ( gxy * gyz - gyy * gxz ) / gupzz
  gupyy =   ( gxx * gzz - gxz * gxz ) / gupzz
  gupyz = - ( gxx * gyz - gxy * gxz ) / gupzz
  gupzz =   ( gxx * gyy - gxy * gxy ) / gupzz

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

! Gam^i_Res = Gam^i + gup^ij_,j
  Gmx_Res = Gmx - (gupxx*(gupxx*gxxx+gupxy*gxyx+gupxz*gxzx)&
                  +gupxy*(gupxx*gxyx+gupxy*gyyx+gupxz*gyzx)&
                  +gupxz*(gupxx*gxzx+gupxy*gyzx+gupxz*gzzx)&
                  +gupxx*(gupxy*gxxy+gupyy*gxyy+gupyz*gxzy)&
                  +gupxy*(gupxy*gxyy+gupyy*gyyy+gupyz*gyzy)&
                  +gupxz*(gupxy*gxzy+gupyy*gyzy+gupyz*gzzy)&
                  +gupxx*(gupxz*gxxz+gupyz*gxyz+gupzz*gxzz)&
                  +gupxy*(gupxz*gxyz+gupyz*gyyz+gupzz*gyzz)&
                  +gupxz*(gupxz*gxzz+gupyz*gyzz+gupzz*gzzz))
  Gmy_Res = Gmy - (gupxx*(gupxy*gxxx+gupyy*gxyx+gupyz*gxzx)&
                  +gupxy*(gupxy*gxyx+gupyy*gyyx+gupyz*gyzx)&
                  +gupxz*(gupxy*gxzx+gupyy*gyzx+gupyz*gzzx)&
                  +gupxy*(gupxy*gxxy+gupyy*gxyy+gupyz*gxzy)&
                  +gupyy*(gupxy*gxyy+gupyy*gyyy+gupyz*gyzy)&
                  +gupyz*(gupxy*gxzy+gupyy*gyzy+gupyz*gzzy)&
                  +gupxy*(gupxz*gxxz+gupyz*gxyz+gupzz*gxzz)&
                  +gupyy*(gupxz*gxyz+gupyz*gyyz+gupzz*gyzz)&
                  +gupyz*(gupxz*gxzz+gupyz*gyzz+gupzz*gzzz))
  Gmz_Res = Gmz - (gupxx*(gupxz*gxxx+gupyz*gxyx+gupzz*gxzx)&
                  +gupxy*(gupxz*gxyx+gupyz*gyyx+gupzz*gyzx)&
                  +gupxz*(gupxz*gxzx+gupyz*gyzx+gupzz*gzzx)&
                  +gupxy*(gupxz*gxxy+gupyz*gxyy+gupzz*gxzy)&
                  +gupyy*(gupxz*gxyy+gupyz*gyyy+gupzz*gyzy)&
                  +gupyz*(gupxz*gxzy+gupyz*gyzy+gupzz*gzzy)&
                  +gupxz*(gupxz*gxxz+gupyz*gxyz+gupzz*gxzz)&
                  +gupyz*(gupxz*gxyz+gupyz*gyyz+gupzz*gyzz)&
                  +gupzz*(gupxz*gxzz+gupyz*gyzz+gupzz*gzzz))

  call bssn2adm(ex,chin1,trK,gxx,gxy,gxz,gyy,gyz,gzz, &
                             Axx,Axy,Axz,Ayy,Ayz,Azz, &
              adm_dxx,adm_gxy,adm_gxz,adm_dyy,adm_gyz,adm_dzz, &
              Kxx,Kxy,Kxz,Kyy,Kyz,Kzz)    
  adm_dxx = adm_dxx - ONE  
  adm_dyy = adm_dyy - ONE  
  adm_dzz = adm_dzz - ONE  

  call constraint_adm_ss(ex,crho,sigma,R, X, Y, Z,&
               drhodx, drhody, drhodz,                                         &
               dsigmadx,dsigmady,dsigmadz,                                     &
               dRdx,dRdy,dRdz,                                                 &
               drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                &
               dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,    &
               dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,                            &
               adm_dxx,adm_gxy,adm_gxz,adm_dyy,adm_gyz,adm_dzz, &
               Kxx,Kxy,Kxz,Kyy,Kyz,Kzz, &
               Lap,Sfx,Sfy,Sfz,rho,Sx,Sy,Sz,&
               Gamxxx, Gamxxy, Gamxxz,Gamxyy, Gamxyz, Gamxzz, &
               Gamyxx, Gamyxy, Gamyxz,Gamyyy, Gamyyz, Gamyzz, &
               Gamzxx, Gamzxy, Gamzxz,Gamzyy, Gamzyz, Gamzzz, &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz, &
               ham_Res, movx_Res, movy_Res, movz_Res, &
               Symmetry,Lev,sst)

  return

  end subroutine constraint_bssn_ss
#else
#error "not recognized ABV"
#endif
