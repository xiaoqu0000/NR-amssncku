

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

! ham_Res = trR + 2/3 * K^2 - A_ij * A^ij - 16 * PI * rho
! here trR is respect to physical metric
  ham_Res =   gupxx * Rxx + gupyy * Ryy + gupzz * Rzz + &
        TWO* ( gupxy * Rxy + gupxz * Rxz + gupyz * Ryz )

  ham_Res = chin1*ham_Res + F2o3 * trK * trK -(&
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
  call fderivs(ex,chi,chix,chiy,chiz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)
  call fderivs(ex,Axx,gxxx,gxxy,gxxz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)
  call fderivs(ex,Axy,gxyx,gxyy,gxyz,X,Y,Z,ANTI,ANTI,SYM ,Symmetry,0)
  call fderivs(ex,Axz,gxzx,gxzy,gxzz,X,Y,Z,ANTI,SYM ,ANTI,Symmetry,0)
  call fderivs(ex,Ayy,gyyx,gyyy,gyyz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)
  call fderivs(ex,Ayz,gyzx,gyzy,gyzz,X,Y,Z,SYM ,ANTI,ANTI,Symmetry,0)
  call fderivs(ex,Azz,gzzx,gzzy,gzzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)

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
movx_Res = gupxx*gxxx + gupyy*gxyy + gupzz*gxzz &
          +gupxy*gxyx + gupxz*gxzx + gupyz*gxzy &
          +gupxy*gxxy + gupxz*gxxz + gupyz*gxyz
movy_Res = gupxx*gxyx + gupyy*gyyy + gupzz*gyzz &
          +gupxy*gyyx + gupxz*gyzx + gupyz*gyzy &
          +gupxy*gxyy + gupxz*gxyz + gupyz*gyyz
movz_Res = gupxx*gxzx + gupyy*gyzy + gupzz*gzzz &
          +gupxy*gyzx + gupxz*gzzx + gupyz*gzzy &
          +gupxy*gxzy + gupxz*gxzz + gupyz*gyzz

!store K,i in chi,i
  call fderivs(ex,trK,chix,chiy,chiz,X,Y,Z,SYM,SYM,SYM,Symmetry,0)

movx_Res = movx_Res - F2o3*chix - F8*PI*sx
movy_Res = movy_Res - F2o3*chiy - F8*PI*sy
movz_Res = movz_Res - F2o3*chiz - F8*PI*sz

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

! ham_Res = trR + 2/3 * K^2 - A_ij * A^ij - 16 * PI * rho
! here trR is respect to physical metric
  ham_Res =   gupxx * Rxx + gupyy * Ryy + gupzz * Rzz + &
        TWO* ( gupxy * Rxy + gupxz * Rxz + gupyz * Ryz )

  ham_Res = chin1*ham_Res + F2o3 * trK * trK -(&
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
movx_Res = gupxx*gxxx + gupyy*gxyy + gupzz*gxzz &
          +gupxy*gxyx + gupxz*gxzx + gupyz*gxzy &
          +gupxy*gxxy + gupxz*gxxz + gupyz*gxyz
movy_Res = gupxx*gxyx + gupyy*gyyy + gupzz*gyzz &
          +gupxy*gyyx + gupxz*gyzx + gupyz*gyzy &
          +gupxy*gxyy + gupxz*gxyz + gupyz*gyyz
movz_Res = gupxx*gxzx + gupyy*gyzy + gupzz*gzzz &
          +gupxy*gyzx + gupxz*gzzx + gupyz*gzzy &
          +gupxy*gxzy + gupxz*gxzz + gupyz*gyzz

!store K,i in chi,i 
  call fderivs_shc(ex,trK,chix,chiy,chiz,crho,sigma,R, SYM, SYM,SYM,Symmetry,Lev,sst,                &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)

movx_Res = movx_Res - F2o3*chix - F8*PI*sx
movy_Res = movy_Res - F2o3*chiy - F8*PI*sy
movz_Res = movz_Res - F2o3*chiz - F8*PI*sz

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
