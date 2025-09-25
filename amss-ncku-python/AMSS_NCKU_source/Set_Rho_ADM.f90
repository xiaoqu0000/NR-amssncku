
! 定义与 F(R) 标量张量理论相关的标量场分布与相互作用势
! 1: Case C of 1112.3928, V=0
! 2: shell with a2^2*phi0/(1+a2^2), f(R) = R+a2*R^2 induced V
! 3: ground state of Schrodinger-Newton system, f(R) = R+a2*R^2 induced V
! 4: a2 = oo and \phi = \phi_0*0.5*(tanh((r+r_0)/\sigma)-tanh((r-r_0)/\sigma))
! 5: shell with phi0*dexp(-(r-r0)**2/sigma), V = 0

! 原来的方式，手动定义预处理宏
! #define CC 2
! 新的方式，根据 "macrodef.fh" 文件中的预处理宏来定义
#include "macrodef.fh"
#define CC EScalar_CC

subroutine setparameters(a2,r0,phi0,sigma,l2)
implicit none
real*8,intent(out) :: a2,r0,phi0,sigma,l2

! 原来的操作：一个一个读入参数
! call seta2(a2)
! call setphi0(phi0)

! 现在的操作：一次性读入所有参数
call set_escalar_parameter(a2, phi0, r0, sigma, l2)

! r0=120.d0
! sigma=8.d0
! l2=1.d4

! write(*,*)
! write(*,*) " Set_Rho_ADM.f90 a2     = ", a2
! write(*,*) " Set_Rho_ADM.f90 phi0   = ", phi0
! write(*,*) " Set_Rho_ADM.f90 r0     = ", r0
! write(*,*) " Set_Rho_ADM.f90 sigma0 = ", sigma
! write(*,*) " Set_Rho_ADM.f90 l2     = ", l2
! write(*,*)

return

end subroutine setparameters
!===================================================================
function phi(X,Y,Z) result(gont)
implicit none

double precision,intent(in)::X
double precision,intent(in)::Y
double precision,intent(in)::Z
real*8 :: gont

real*8 ::r
real*8 :: a2,r0,phi0,sigma,l2

  call setparameters(a2,r0,phi0,sigma,l2)
  r=dsqrt(X*X+Y*Y+Z*Z)
#if  ( CC == 1)  
! configuration 1  
  gont = phi0*dtanh((r-r0)/sigma)
#elif ( CC == 2)
! configuration 2  
  phi0 = a2**2*phi0/(1+a2**2)
  gont = phi0*dexp(-(r-r0)**2/sigma)
#elif ( CC == 3)
  gont =  (0.0481646d0*dexp(-0.0581545d0*(r-1.8039d-8)*(r-1.8039d-8)/l2) &
           +0.298408d0*dexp(-0.111412d0*(r+9.6741d-9)*(r+9.6741d-9)/l2)+   &
             0.42755d0*dexp(-0.207156d0*(r-1.09822d-8)*(r-1.09822d-8)/l2)+   &
            0.204229d0*dexp(-0.37742d0*(r+2.13778d-8)*(r+2.13778d-8)/l2)+   &
            0.021649d0*dexp(-0.68406d0*(r-8.78608d-8)*(r-8.78608d-8)/l2))/l2
#elif ( CC == 4)
! configuration 4, a2 = oo
  phi0 = 0.5d0*phi0
  gont = phi0*(dtanh((r+r0)/sigma)-dtanh((r-r0)/sigma))
#elif ( CC == 5)
! configuration 5  
  gont = phi0*dexp(-(r-r0)**2/sigma)
#endif

return

end function phi

! d phi/dr
function dphi(X,Y,Z) result(gont)
implicit none

double precision,intent(in)::X
double precision,intent(in)::Y
double precision,intent(in)::Z
real*8 :: gont

real*8 ::r
real*8 :: a2,r0,phi0,sigma,l2

  call setparameters(a2,r0,phi0,sigma,l2)
  r=dsqrt(X*X+Y*Y+Z*Z)
#if  ( CC == 1)  
! configuration 1  
  gont = phi0/sigma*(1-(dtanh((r-r0)/sigma))**2)
#elif ( CC == 2)
! configuration 2
  phi0 = a2**2*phi0/(1+a2**2)
  gont = -2.d0*phi0*(r-r0)/sigma*exp(-(r-r0)**2/sigma)
#elif ( CC == 3)
  gont = (-0.5601976461d-2*(r-0.18039d-7)/l2*dexp(-0.581545d-1*(r-0.18039d-7)**2/l2) &
          -0.6649246419d-1*(r+0.96741d-8)/l2*dexp(-0.111412d0*(r+.96741e-8)**2/l2)    &
          -0.1771390956d0*(r-0.109822d-7)/l2*dexp(-0.207156d0*(r-0.109822d-7)**2/l2)  &
          -0.1541602184d0*(r+0.213778d-7)/l2*dexp(-0.37742d0*(r+0.213778d-7)**2/l2)   &
          -0.2961842988d-1*(r-0.878608d-7)/l2*dexp(-0.68406*(r-0.878608d-7)**2/l2))/l2
#elif ( CC == 4)
! configuration 4, a2 = oo
  phi0 = 0.5d0*phi0
  gont = phi0*((1-dtanh((r+r0)/sigma)**2)/sigma- &
               (1-dtanh((r-r0)/sigma)**2)/sigma)
#elif ( CC == 5)
! configuration 5
  gont = -2.d0*phi0*(r-r0)/sigma*exp(-(r-r0)**2/sigma)
#endif

return

end function dphi
!==================================================================
function potential(X,Y,Z) result(gont)
implicit none

double precision,intent(in)::X
double precision,intent(in)::Y
double precision,intent(in)::Z
real*8 :: gont

real*8 :: phi
real*8 :: PI,v

real*8 :: a2,r0,phi0,sigma,l2

#if  ( CC == 1 || CC == 4 || CC == 5)  
  gont = 0.d0

#elif ( CC == 2 || CC == 3)
  call setparameters(a2,r0,phi0,sigma,l2)
  PI = dacos(-1.d0)

  v = phi(X,Y,Z)

  gont = dexp(-8.d0*dsqrt(PI/3)*v)*(1-dexp(4*dsqrt(PI/3)*v))**2/32/PI/a2
#endif

return

end function potential
!==================================================================
!Note this part is for evolution 
!not just for initial configuration

!f(R) potential F=R+a_2R^2
subroutine frpotential(ex,Sphi,V,dVdSphi)

implicit none

integer,intent(in ):: ex(1:3)
real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Sphi
real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: V,dVdSphi

real*8 :: a2,r0,phi0,sigma,l2
real*8, parameter :: Four = 4.d0, TWO = 2.d0,ONE = 1.d0,ZEO=0.d0
real*8 :: PI

  PI = dacos(-ONE)

#if  ( CC == 1 || CC == 4 || CC == 5)  
  V = ZEO
  dVdSphi = ZEO
#elif ( CC == 2 || CC == 3)
  call setparameters(a2,r0,phi0,sigma,l2)
  V = dexp(-8.d0*dsqrt(PI/3)*Sphi)*(1-dexp(4*dsqrt(PI/3)*Sphi))**2/32/PI/a2
  dVdSphi = 1.d0/a2/1.2d1*dsqrt(3.d0/PI)*dexp(-8.d0*dsqrt(PI/3.d0)*Sphi)*(-1+dexp(4*dsqrt(Pi/3)*Sphi))
#endif

return

end subroutine frpotential
!==================================================================
!f(R) potential F=R+a_2R^2
!fprim(R) = 1+2*a_2*R
subroutine frfprim(ex,RR,fprim)

implicit none

integer,intent(in ):: ex(1:3)
real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: RR
real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: fprim

real*8 :: a2,r0,phi0,sigma,l2
real*8, parameter :: ONE=1.d0, TWO = 2.d0

#if  ( CC == 1 || CC == 4 || CC == 5)  
  fprim = ONE
#elif ( CC == 2 || CC == 3)
  call setparameters(a2,r0,phi0,sigma,l2)
  fprim = ONE+TWO*a2*RR
#endif

return

end subroutine frfprim
!==================================================================
subroutine set_rho_adm2(ex,rho,X,Y,Z)

implicit none
! argument variables
integer,intent(in)::ex
double precision,intent(in),dimension(ex)::X
double precision,intent(in),dimension(ex)::Y
double precision,intent(in),dimension(ex)::Z
double precision,intent(out),dimension(ex)::rho

integer :: i
real*8 :: dphi

  do i=1,ex
    ! rho(i) = dphi(X,Y,Z)
    rho(i) = dphi(X(i),Y(i),Z(i))
    rho(i) = rho(i)*rho(i)
  enddo

  return

end subroutine set_rho_adm2

subroutine set_rho_adm1(ex,rho,X,Y,Z)

implicit none
! argument variables
integer,intent(in)::ex
double precision,intent(in),dimension(ex)::X
double precision,intent(in),dimension(ex)::Y
double precision,intent(in),dimension(ex)::Z
double precision,intent(out),dimension(ex)::rho

real*8 :: potential
integer :: i

  do i=1,ex
    rho(i) = potential(X(i),Y(i),Z(i))
  enddo

  return

end subroutine set_rho_adm1

subroutine set_rho_adm(ex,rho,X,Y,Z)

implicit none
! argument variables
integer,intent(in)::ex
double precision,intent(in),dimension(ex)::X
double precision,intent(in),dimension(ex)::Y
double precision,intent(in),dimension(ex)::Z
! in psivac, out rho_adm
double precision,intent(inout),dimension(ex)::rho

double precision,dimension(ex)::rho1,rho2

  call set_rho_adm1(ex,rho1,X,Y,Z)
  call set_rho_adm2(ex,rho2,X,Y,Z)

  rho = rho**4
  rho = rho**2*rho1+rho*rho2

  return

end subroutine set_rho_adm
