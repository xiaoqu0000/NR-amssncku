

#include "macrodef.fh"

  subroutine get_initial_nbhs_null(ex,crho,sigma,x,RJ,IJ,omega,sst,Rmin)
  
  implicit none
!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3),sst
  real*8,intent(in ) :: Rmin
  double precision,intent(in),dimension(ex(1))::crho
  double precision,intent(in),dimension(ex(2))::sigma
  double precision,intent(in),dimension(ex(3))::x
  double precision,intent(inout),dimension(ex(1),ex(2),ex(3))::RJ,IJ,omega

!~~~~~~> Other variables:
  real*8 :: xe
  real*8,dimension(ex(1),ex(2)) :: RJe,IJe
  integer :: k
  
  xe = x(1)
  RJe = RJ(:,:,1)
  IJe = IJ(:,:,1)

  do k=1,ex(3)
    RJ(:,:,k) = RJe*(1.d0-x(k))*xe/(1-xe)/x(k)
    IJ(:,:,k) = IJe*(1.d0-x(k))*xe/(1-xe)/x(k)
  enddo

  omega = 1.d0

  return

  end subroutine get_initial_nbhs_null
!-----------------------------------
!Eq.(10) of CQG 24, S327 (2007)
!----------------------------------
  function Zslm(s,l,m,the,phi) result(gont)
  implicit none
  integer,intent(in) :: s,l,m
  real*8,intent(in) :: the,phi

  double complex :: Yslm,gont,II

  II=dcmplx(0.d0,1.d0)

  if(m>0)then
    gont = Yslm(s,l,m,the,phi)
    if(m/2*2==m)then
      gont = gont+Yslm(s,l,-m,the,phi)
    else
      gont = gont-Yslm(s,l,-m,the,phi)
    endif
    gont = gont/dsqrt(2.d0)
  elseif(m<0)then
    gont = -Yslm(s,l,-m,the,phi)
    if(m/2*2==m)then
      gont = gont+Yslm(s,l,m,the,phi)
    else
      gont = gont-Yslm(s,l,m,the,phi)
    endif
    gont = II*gont/dsqrt(2.d0)
  else
    gont = Yslm(s,l,m,the,phi)
  endif

  return

  end function Zslm

!#define SCH

#ifdef SCH

subroutine get_initial_null(ex,crho,sigma,R,RJ,IJ,sst,Rmin)
implicit none
! argument variables
integer, intent(in ):: ex(1:3),sst
real*8,intent(in) :: Rmin
double precision,intent(in),dimension(ex(1))::crho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::R
real*8,dimension(ex(1),ex(2),ex(3)),intent(out)::RJ,IJ

RJ = 0.d0
IJ = 0.d0

return

end subroutine get_initial_null
!-------------------------
  subroutine get_null_boundary(ex,crho,sigma,R,beta,RQ,IQ,RU,IU,W,RTheta,ITheta, &
                        quR1,quR2,quI1,quI2,qlR1,qlR2,qlI1,qlI2,          &
                        gR,gI,                                            &
                        dquR1,dquR2,dquI1,dquI2,                          &
                        bdquR1,bdquR2,bdquI1,bdquI2,                      &
                        dgR,dgI,bdgR,bdgI,T,Rmin,sst)

  implicit none

!~~~~~~% Input parameters:
  integer,intent(in) :: ex(3),sst
  real*8,intent(in) :: T,Rmin
  double precision,intent(in),dimension(ex(1))::crho
  double precision,intent(in),dimension(ex(2))::sigma
  double precision,intent(in),dimension(ex(3))::R
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout) :: beta,RQ,IQ,RU,IU
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout) :: W,RTheta,ITheta
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: quR1,quR2,quI1,quI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: qlR1,qlR2,qlI1,qlI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: gR,gI
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dquR1,dquR2,dquI1,dquI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: bdquR1,bdquR2,bdquI1,bdquI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dgR,dgI,bdgR,bdgI

  integer :: k
     k=1

     beta(:,:,k) = 0.d0
     W(:,:,k) =-2.d0/R(k)**2/Rmin**2*(1.d0-R(k))**2

     RQ(:,:,k) = 0.d0
     IQ(:,:,k) = 0.d0

     RTheta(:,:,k) = 0.d0
     ITheta(:,:,k) = 0.d0

     RU(:,:,k) = 0.d0
     IU(:,:,k) = 0.d0

  return

  end subroutine get_null_boundary
!-------------------------------------------------------------  
subroutine get_exact_null(ex,crho,sigma,R,RJ,IJ,sst,Rmin,T,               &
                        quR1,quR2,quI1,quI2,qlR1,qlR2,qlI1,qlI2,          &
                        gR,gI,                                            &
                        dquR1,dquR2,dquI1,dquI2,                          &
                        bdquR1,bdquR2,bdquI1,bdquI2,                      &
                        dgR,dgI,bdgR,bdgI)
implicit none
! argument variables
integer, intent(in ):: ex(1:3),sst
real*8,intent(in) :: Rmin,T
double precision,intent(in),dimension(ex(1))::crho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::R
real*8,dimension(ex(1),ex(2),ex(3)),intent(out)::RJ,IJ
real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: quR1,quR2,quI1,quI2
real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: qlR1,qlR2,qlI1,qlI2
real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: gR,gI
real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dquR1,dquR2,dquI1,dquI2
real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: bdquR1,bdquR2,bdquI1,bdquI2
real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dgR,dgI,bdgR,bdgI

RJ = 0.d0
IJ = 0.d0

return

end subroutine get_exact_null
!-------------------------------------------------------------------------------------------
  subroutine get_null_boundary_c(ex,crho,sigma,R,beta,RQ,IQ,RU,IU,W,RTheta,ITheta, &
                        quR1,quR2,quI1,quI2,qlR1,qlR2,qlI1,qlI2,          &
                        gR,gI,                                            &
                        dquR1,dquR2,dquI1,dquI2,                          &
                        bdquR1,bdquR2,bdquI1,bdquI2,                      &
                        dgR,dgI,bdgR,bdgI,T,Rmin,sst)

  implicit none

!~~~~~~% Input parameters:
  integer,intent(in) :: ex(3),sst
  real*8,intent(in) :: T,Rmin
  double precision,intent(in),dimension(ex(1))::crho
  double precision,intent(in),dimension(ex(2))::sigma
  double precision,intent(in),dimension(ex(3))::R
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout) :: beta,RQ,IQ,RU,IU
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout) :: W,RTheta,ITheta
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: quR1,quR2,quI1,quI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: qlR1,qlR2,qlI1,qlI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: gR,gI
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dquR1,dquR2,dquI1,dquI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: bdquR1,bdquR2,bdquI1,bdquI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dgR,dgI,bdgR,bdgI

  integer :: k

  do k=1,ex(3)

     beta(:,:,k) = 0.d0
     W(:,:,k) =-2.d0/R(k)**2/Rmin**2*(1.d0-R(k))**2

     RQ(:,:,k) = 0.d0
     IQ(:,:,k) = 0.d0

     RTheta(:,:,k) = 0.d0
     ITheta(:,:,k) = 0.d0

     RU(:,:,k) = 0.d0
     IU(:,:,k) = 0.d0
 enddo

  return

  end subroutine get_null_boundary_c

#else

#if 0
! for some trival check
#if 1
!-------------------------------------------------------------  
! Linear wave given in CQG 24S327
!-------------------------------------------------------------
subroutine get_initial_null(ex,crho,sigma,R,RJ,IJ,sst,Rmin)
implicit none
! argument variables
integer, intent(in ):: ex(1:3),sst
real*8,intent(in) :: Rmin
double precision,intent(in),dimension(ex(1))::crho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::R
real*8,dimension(ex(1),ex(2),ex(3)),intent(out)::RJ,IJ

integer :: i,j,k
real*8 ::x,y,z,gr,gt,gp,tgrho,tgsigma
double complex :: Yslm,II,Jr

double complex :: beta0,C1,C2
integer :: nu,m

call initial_null_paramter(beta0,C1,C2,nu,m)

  II = dcmplx(0.d0,1.d0)

  do i=1,ex(1)
  do j=1,ex(2)
  do k=1,ex(3)
! here fake global coordinate is enough  
    gr = 1.d0
    tgrho = dtan(crho(i))
    tgsigma = dtan(sigma(j))
    select case (sst)
    case (0)
      z = gr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = z*tgrho
      y = z*tgsigma
    case (1)
      z = -gr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = z*tgrho
      y = z*tgsigma
    case (2)
      x = gr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      y = x*tgrho
      z = x*tgsigma
    case (3)
      x = -gr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      y = x*tgrho
      z = x*tgsigma
    case (4)
      y = gr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = y*tgrho
      z = y*tgsigma
    case (5)
      y = -gr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = y*tgrho
      z = y*tgsigma
    case default
      write(*,*) "get_initial_null: not recognized sst = ",sst
      return
    end select
    gt = dacos(z/gr)
    gp = datan2(y,x)

    gr = (1.d0-R(k))/R(k)/Rmin
    Jr = (2.4d1*beta0+3.d0*II*nu*C1-II*nu**3*C2)/3.6d1+C1/4.d0*gr-C2/1.2d1*gr**3
    gr = dreal(Jr)
    Jr = Yslm(0,2,m,gt,gp)
    RJ(i,j,k) = gr*dreal(Jr)
    IJ(i,j,k) = gr*dimag(Jr)

#if 0 
    RJ(i,j,k) = 0.25d0*dsqrt(5.d0/3.1415926)*(3/(1.d0+tgrho*tgrho+tgsigma*tgsigma)-1.d0)
    IJ(i,j,k) = 0.d0
#endif
  enddo
  enddo
  enddo

return

end subroutine get_initial_null
#else
! for check usage
subroutine get_initial_null(ex,crho,sigma,R,RJ,IJ,sst,Rmin)
implicit none
! argument variables
integer, intent(in ):: ex(1:3),sst
real*8,intent(in) :: Rmin
double precision,intent(in),dimension(ex(1))::crho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::R
real*8,dimension(ex(1),ex(2),ex(3)),intent(out)::RJ,IJ

real*8 :: thetac,thetas,sr,ss,cr,cs,srss,crcs,tcts,tcts2
real*8 :: sr2,ss2,cr2,cs2,tc2,ts2
integer :: i,j,k
real*8 :: ggr,tgrho,tgsigma

real*8 :: PI

PI = dacos(-1.d0)

  do i=1,ex(1)
  do j=1,ex(2)
  do k=1,ex(3)
     sr = dsin(crho(i))
     ss = dsin(sigma(j))
     cr = dcos(crho(i))
     cs = dcos(sigma(j))
     srss = sr*ss
     crcs = cr*cs
     sr2 = sr*sr
     ss2 = ss*ss
     cr2 = cr*cr
     cs2 = cs*cs
     thetac = dsqrt((1.d0-srss)/2.d0) 
     thetas = dsqrt((1.d0+srss)/2.d0) 
     tc2 = thetac*thetac
     ts2 = thetas*thetas
     tcts = thetac*thetas
     tcts2 = tcts*tcts
! q^Aq^B@_A@_B Y20
     RJ(i,j,k) =-1.5d0*dsqrt(5.d0/PI)*sr*ss*(4.d0*cr2*cs2+cs2+cr2)
     IJ(i,j,k) = 3.d0*dsqrt(5.d0/PI)*thetac*thetas*(cs2-cr2)

! @_rho@_rho Y20
     RJ(i,j,k) = 1.5d0*dsqrt(5.d0/PI)*cs2*cs2*(-cr2*cs2+2*cr2*cr2*cs2-2*cr2*cr2-cs2+3*cr2) &
                 /(3*cr2*cs2*cs2*cs2-3*cr2*cr2*cr2*cs2*cs2-3*cr2*cr2*cs2*cs2*cs2+          &
                 cr2*cr2*cr2*cs2*cs2*cs2-3*cr2*cr2*cs2-cs2*cs2*cs2-cr2*cr2*cr2+            &
                 3*cs2*cr2*cr2*cr2+6*cr2*cr2*cs2*cs2-3*cs2*cs2*cr2)
     IJ(i,j,k) = 0.d0
! q^Aq^B h_AB
     RJ(i,j,k) = 0.d0
     IJ(i,j,k) = 0.d0
  enddo
  enddo
  enddo

return

end subroutine get_initial_null
#endif
#endif
!======================================================================================
!-------------------------------------------------------------  
! Linear wave given in CQG 24S327
!-------------------------------------------------------------
subroutine get_initial_null(ex,crho,sigma,R,RJ,IJ,sst,Rmin)
implicit none
! argument variables
integer, intent(in ):: ex(1:3),sst
real*8,intent(in) :: Rmin
double precision,intent(in),dimension(ex(1))::crho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::R
real*8,dimension(ex(1),ex(2),ex(3)),intent(out)::RJ,IJ

integer :: i,j,k
real*8 ::x,y,z,gr,gt,gp,tgrho,tgsigma,tc,ts
double complex :: Yslm,II,Jr

double complex :: beta0,C1,C2
integer :: nu,m

double complex :: swtf,ff

call initial_null_paramter(beta0,C1,C2,nu,m)

  II = dcmplx(0.d0,1.d0)

  do i=1,ex(1)
  do j=1,ex(2)
  do k=1,ex(3)
!fake global coordinate is enough here  
    gr = 1.d0
    tgrho = dtan(crho(i))
    tgsigma = dtan(sigma(j))
    tc = dsqrt((1.d0-dsin(crho(i))*dsin(sigma(j)))/2.d0)
    ts = dsqrt((1.d0+dsin(crho(i))*dsin(sigma(j)))/2.d0)
    select case (sst)
    case (0)
      z = gr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = z*tgrho
      y = z*tgsigma
    case (1)
      z = -gr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = z*tgrho
      y = z*tgsigma
    case (2)
      x = gr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      y = x*tgrho
      z = x*tgsigma
    case (3)
      x = -gr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      y = x*tgrho
      z = x*tgsigma
    case (4)
      y = gr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = y*tgrho
      z = y*tgsigma
    case (5)
      y = -gr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = y*tgrho
      z = y*tgsigma
    case default
      write(*,*) "get_initial_null: not recognized sst = ",sst
      return
    end select
    gt = dacos(z/gr)
    gp = datan2(y,x)
    swtf = 2.d0*tc*ts*(ts+II*tc)/dcos(sigma(j))
    if(sst==1 .or. sst==3 .or. sst==4) swtf = dconjg(swtf)
    select case (sst)
    case (0,1)
      swtf = swtf/(dcos(gp)+II*dcos(gt)*dsin(gp))*(dcos(gt)**2+dsin(gt)**2*dcos(gp)**2)
    case (2,3)
      swtf = II*swtf*dsin(gt)
    case (4,5)
      swtf = -II*swtf*dsin(gt)
    end select

    gr = (1.d0-R(k))/R(k)/Rmin
    Jr = (2.4d1*beta0+3.d0*II*nu*C1-II*nu**3*C2)/3.6d1+C1/4.d0*gr-C2/1.2d1*gr**3
    gr = dreal(Jr)
    Jr = Yslm(2,2,m,gt,gp)
    ff = dsqrt(dble((2-1)*2*(2+1)*(2+2)))*gr*Jr*swtf**2
    RJ(i,j,k) = dreal(ff)
    IJ(i,j,k) = dimag(ff)

  enddo
  enddo
  enddo

return

end subroutine get_initial_null
!==============================================================================================

#if 0
! for checking derivs_eth and dderivs_eth
!-------------------------
  subroutine get_null_boundary(ex,crho,sigma,R,beta,RQ,IQ,RU,IU,W,RTheta,ITheta, &
                        quR1,quR2,quI1,quI2,qlR1,qlR2,qlI1,qlI2,          &
                        gR,gI,                                            &
                        dquR1,dquR2,dquI1,dquI2,                          &
                        bdquR1,bdquR2,bdquI1,bdquI2,                      &
                        dgR,dgI,bdgR,bdgI,T,Rmin,sst)

  implicit none

!~~~~~~% Input parameters:
  integer,intent(in) :: ex(3),sst
  real*8,intent(in) :: T,Rmin
  double precision,intent(in),dimension(ex(1))::crho
  double precision,intent(in),dimension(ex(2))::sigma
  double precision,intent(in),dimension(ex(3))::R
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout) :: beta,RQ,IQ,RU,IU
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout) :: W,RTheta,ITheta
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: quR1,quR2,quI1,quI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: qlR1,qlR2,qlI1,qlI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: gR,gI
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dquR1,dquR2,dquI1,dquI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: bdquR1,bdquR2,bdquI1,bdquI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dgR,dgI,bdgR,bdgI

  double complex,dimension(ex(1),ex(2)) :: Y20,dY20,ddY20,f
  integer :: i,j,k
  real*8 ::x,y,z,hgr,gt,gp,tgrho,tgsigma
  double complex :: Yslm,II,Jr

double complex :: beta0,C1,C2
integer :: nu,m

call initial_null_paramter(beta0,C1,C2,nu,m)

  II = dcmplx(0.d0,1.d0)

  k=1
  do i=1,ex(1)
  do j=1,ex(2)
! fake global coordinate is enough  
    hgr = 1.d0
    tgrho = dtan(crho(i))
    tgsigma = dtan(sigma(j))
    select case (sst)
    case (0)
      z = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = z*tgrho
      y = z*tgsigma
    case (1)
      z = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = z*tgrho
      y = z*tgsigma
    case (2)
      x = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      y = x*tgrho
      z = x*tgsigma
    case (3)
      x = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      y = x*tgrho
      z = x*tgsigma
    case (4)
      y = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = y*tgrho
      z = y*tgsigma
    case (5)
      y = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = y*tgrho
      z = y*tgsigma
    case default
      write(*,*) "get_initial_null: not recognized sst = ",sst
      return
    end select
    gt = dacos(z/hgr)
    gp = datan2(y,x)

    hgr = (1.d0-R(k))/R(k)/Rmin
    Jr = (2.4d1*beta0+3.d0*II*nu*C1-II*nu**3*C2)/3.6d1+C1/4.d0*hgr-C2/1.2d1*hgr**3
    RTheta(i,j,k) = dreal(Jr*nu*(II*dcos(nu*T)-dsin(nu*T)))
    f(i,j) = dreal(Jr*(dcos(nu*T)+II*dsin(nu*T)))

    Jr = (-2.4d1*II*nu*beta0+3.d0*nu*nu*C1-nu**4*C2)/36.d0+2.d0*beta0*hgr&
        +C1/2.d0*hgr*hgr+II*nu*C2/3.d0*hgr**3+C2/4.d0*hgr**4
    RU(i,j,k) = dreal(Jr*(dcos(nu*T)+II*dsin(nu*T)))

! Re(r^2*Ul_,r*exp(i nu T)) of CQG 24S327, (12) indeed    
    Jr = -2.d0*beta0-C1*hgr-II*nu*C2*hgr**2-C2*hgr**3
    RQ(i,j,k) = dreal(Jr*(dcos(nu*T)+II*dsin(nu*T)))

    Jr = (2.4d1*II*nu*beta0-3.d0*nu*nu*C1+nu**4*C2)/6.d0+(3.d0*II*nu*C1-6.d0*beta0-II*nu**3*C2)/3.d0*hgr&
        -nu*nu*C2*hgr*hgr+II*nu*C2*hgr**3+C2/2.d0*hgr**4
    W(i,j,k) = dreal(Jr*(dcos(nu*T)+II*dsin(nu*T)))

    Y20(i,j) = Yslm(0,2,m,gt,gp)

  enddo
  enddo

     call derivs_eth(ex(1:2),crho,sigma,Y20,dY20,0,1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call dderivs_eth(ex(1:2),crho,sigma,Y20,ddY20,0,1,1,  &
                      quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k),&
                      dquR1(:,:,k),dquR2(:,:,k),dquI1(:,:,k),dquI2(:,:,k),       &
                      bdquR1(:,:,k),bdquR2(:,:,k),bdquI1(:,:,k),bdquI2(:,:,k),   &
                      dgR(:,:,k),dgI(:,:,k),bdgR(:,:,k),bdgI(:,:,k))

     beta(:,:,k) = dreal(beta0*(dcos(nu*T)+II*dsin(nu*T)))*dreal(Y20)
     W(:,:,k) = W(:,:,k)*dreal(Y20)

     f = dexp(-2.d0*beta(:,:,k))*(f*ddY20*RQ(:,:,k)*dconjg(dY20)+dsqrt(1.d0+abs(f*ddY20))*RQ(:,:,k)*dY20)
     RQ(:,:,k) = dreal(f)
     IQ(:,:,k) = dimag(f)

     f = ddY20*RTheta(:,:,k)
     RTheta(:,:,k) = dreal(f)
     ITheta(:,:,k) = dimag(f)

     f = dY20*RU(:,:,k)
     RU(:,:,k) = dreal(f)
     IU(:,:,k) = dimag(f)

  return

  end subroutine get_null_boundary
#else
!-------------------------
  subroutine get_null_boundary(ex,crho,sigma,R,beta,RQ,IQ,RU,IU,W,RTheta,ITheta, &
                        quR1,quR2,quI1,quI2,qlR1,qlR2,qlI1,qlI2,          &
                        gR,gI,                                            &
                        dquR1,dquR2,dquI1,dquI2,                          &
                        bdquR1,bdquR2,bdquI1,bdquI2,                      &
                        dgR,dgI,bdgR,bdgI,T,Rmin,sst)

  implicit none

!~~~~~~% Input parameters:
  integer,intent(in) :: ex(3),sst
  real*8,intent(in) :: T,Rmin
  double precision,intent(in),dimension(ex(1))::crho
  double precision,intent(in),dimension(ex(2))::sigma
  double precision,intent(in),dimension(ex(3))::R
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout) :: beta,RQ,IQ,RU,IU
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout) :: W,RTheta,ITheta
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: quR1,quR2,quI1,quI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: qlR1,qlR2,qlI1,qlI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: gR,gI
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dquR1,dquR2,dquI1,dquI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: bdquR1,bdquR2,bdquI1,bdquI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dgR,dgI,bdgR,bdgI

integer :: i,j,k
real*8 ::x,y,z,hgr,gt,gp,tgrho,tgsigma,tc,ts,rf
double complex :: Yslm,II,Jr,swtf,ff

double complex :: beta0,C1,C2
integer :: nu,m

call initial_null_paramter(beta0,C1,C2,nu,m)

  II = dcmplx(0.d0,1.d0)

  k=1
  do i=1,ex(1)
  do j=1,ex(2)
! fake global coordinate is enough here  
    hgr = 1.d0
    tgrho = dtan(crho(i))
    tgsigma = dtan(sigma(j))
    tc = dsqrt((1.d0-dsin(crho(i))*dsin(sigma(j)))/2.d0)
    ts = dsqrt((1.d0+dsin(crho(i))*dsin(sigma(j)))/2.d0)
    select case (sst)
    case (0)
      z = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = z*tgrho
      y = z*tgsigma
    case (1)
      z = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = z*tgrho
      y = z*tgsigma
    case (2)
      x = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      y = x*tgrho
      z = x*tgsigma
    case (3)
      x = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      y = x*tgrho
      z = x*tgsigma
    case (4)
      y = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = y*tgrho
      z = y*tgsigma
    case (5)
      y = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = y*tgrho
      z = y*tgsigma
    case default
      write(*,*) "get_null_boundary: not recognized sst = ",sst
      return
    end select
    gt = dacos(z/hgr)
    gp = datan2(y,x)
    swtf = 2.d0*tc*ts*(ts+II*tc)/dcos(sigma(j))
    if(sst==1 .or. sst==3 .or. sst==4) swtf = dconjg(swtf)
    select case (sst)
    case (0,1)
      swtf = swtf/(dcos(gp)+II*dcos(gt)*dsin(gp))*(dcos(gt)**2+dsin(gt)**2*dcos(gp)**2)
    case (2,3)
      swtf = II*swtf*dsin(gt)
    case (4,5)
      swtf = -II*swtf*dsin(gt)
    end select

    hgr = (1.d0-R(k))/R(k)/Rmin
    Jr = (2.4d1*beta0+3.d0*II*nu*C1-II*nu**3*C2)/3.6d1+C1/4.d0*hgr-C2/1.2d1*hgr**3
    rf = dreal(Jr*nu*(II*dcos(nu*T)-dsin(nu*T)))
    Jr = Yslm(2,2,m,gt,gp)
    ff = dsqrt(dble((2-1)*2*(2+1)*(2+2)))*rf*Jr*swtf**2
    RTheta(i,j,k) = dreal(ff)
    ITheta(i,j,k) = dimag(ff)

    rf = dreal(Yslm(0,2,m,gt,gp))
    beta(i,j,k) = rf*dreal(beta0*(dcos(nu*T)+II*dsin(nu*T)))
    Jr = (2.4d1*II*nu*beta0-3.d0*nu*nu*C1+nu**4*C2)/6.d0+(3.d0*II*nu*C1-6.d0*beta0-II*nu**3*C2)/3.d0*hgr&
        -nu*nu*C2*hgr*hgr+II*nu*C2*hgr**3+C2/2.d0*hgr**4
    W(i,j,k) = rf*dreal(Jr*(dcos(nu*T)+II*dsin(nu*T)))

    Jr = (-2.4d1*II*nu*beta0+3.d0*nu*nu*C1-nu**4*C2)/36.d0+2.d0*beta0*hgr&
        +C1/2.d0*hgr*hgr+II*nu*C2/3.d0*hgr**3+C2/4.d0*hgr**4
    rf = dreal(Jr*(dcos(nu*T)+II*dsin(nu*T)))
    Jr = Yslm(1,2,m,gt,gp)
    ff = dsqrt(dble(2*(2+1)))*rf*Jr*swtf
    RU(i,j,k) = dreal(ff)
    IU(i,j,k) = dimag(ff)
  
    Jr = -2.d0*beta0-C1*hgr-II*nu*C2*hgr**2-C2*hgr**3
    rf = dreal(Jr*(dcos(nu*T)+II*dsin(nu*T)))
    Jr = Yslm(1,2,m,gt,gp)
    ff = dsqrt(dble(2*(2+1)))*rf*Jr*swtf  !! U_,r
    Jr = (2.4d1*beta0+3.d0*II*nu*C1-II*nu**3*C2)/3.6d1+C1/4.d0*hgr-C2/1.2d1*hgr**3
    rf = dreal(Jr*(dcos(nu*T)+II*dsin(nu*T)))
    Jr = Yslm(2,2,m,gt,gp)
    Jr = dsqrt(dble((2-1)*2*(2+1)*(2+2)))*rf*Jr*swtf**2  !! J
    rf = dsqrt(1.d0+abs(Jr)**2)  !! K
    ff = dexp(-2.d0*beta(i,j,k))*(Jr*dconjg(ff)+rf*ff)
    RQ(i,j,k) = dreal(ff)
    IQ(i,j,k) = dimag(ff)

  enddo
  enddo

  return

  end subroutine get_null_boundary

#endif

#if 0
! for checking dderivs_eth
!-------------------------------------------------------------  
! Linear wave given in CQG 24S327
!-------------------------------------------------------------
subroutine get_exact_null(ex,crho,sigma,R,RJ,IJ,sst,Rmin,T,               &
                        quR1,quR2,quI1,quI2,qlR1,qlR2,qlI1,qlI2,          &
                        gR,gI,                                            &
                        dquR1,dquR2,dquI1,dquI2,                          &
                        bdquR1,bdquR2,bdquI1,bdquI2,                      &
                        dgR,dgI,bdgR,bdgI)
implicit none
! argument variables
integer, intent(in ):: ex(1:3),sst
real*8,intent(in) :: Rmin,T
double precision,intent(in),dimension(ex(1))::crho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::R
real*8,dimension(ex(1),ex(2),ex(3)),intent(out)::RJ,IJ
real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: quR1,quR2,quI1,quI2
real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: qlR1,qlR2,qlI1,qlI2
real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: gR,gI
real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dquR1,dquR2,dquI1,dquI2
real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: bdquR1,bdquR2,bdquI1,bdquI2
real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dgR,dgI,bdgR,bdgI

double complex,dimension(ex(1),ex(2)) :: Y20,ddY20,f
integer :: i,j,k
real*8 ::x,y,z,hgr,gt,gp,tgrho,tgsigma
double complex :: Yslm,II,Jr

double complex :: beta0,C1,C2
integer :: nu,m

call initial_null_paramter(beta0,C1,C2,nu,m)

  II = dcmplx(0.d0,1.d0)

  do i=1,ex(1)
  do j=1,ex(2)
  do k=1,ex(3)
! fake global coordinate is enough here  
    hgr = 1.d0
    tgrho = dtan(crho(i))
    tgsigma = dtan(sigma(j))
    select case (sst)
    case (0)
      z = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = z*tgrho
      y = z*tgsigma
    case (1)
      z = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = z*tgrho
      y = z*tgsigma
    case (2)
      x = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      y = x*tgrho
      z = x*tgsigma
    case (3)
      x = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      y = x*tgrho
      z = x*tgsigma
    case (4)
      y = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = y*tgrho
      z = y*tgsigma
    case (5)
      y = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = y*tgrho
      z = y*tgsigma
    case default
      write(*,*) "get_exact_null: not recognized sst = ",sst
      return
    end select
    gt = dacos(z/hgr)
    gp = datan2(y,x)

    hgr = (1.d0-R(k))/R(k)/Rmin
    Jr = (2.4d1*beta0+3.d0*II*nu*C1-II*nu**3*C2)/3.6d1+C1/4.d0*hgr-C2/1.2d1*hgr**3
    Y20(i,j) = Yslm(0,2,m,gt,gp)
    RJ(i,j,k) = dreal(Jr*(dcos(nu*T)+II*dsin(nu*T)))
    IJ(i,j,k) = 0.d0

  enddo
  enddo
  enddo

  k=1
  call dderivs_eth(ex(1:2),crho,sigma,Y20,ddY20,0,1,1,  &
                      quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k),&
                      dquR1(:,:,k),dquR2(:,:,k),dquI1(:,:,k),dquI2(:,:,k),       &
                      bdquR1(:,:,k),bdquR2(:,:,k),bdquI1(:,:,k),bdquI2(:,:,k),   &
                      dgR(:,:,k),dgI(:,:,k),bdgR(:,:,k),bdgI(:,:,k))

  do k=1,ex(3)
     f = ddY20*RJ(:,:,k)
     RJ(:,:,k) = dreal(f)
     IJ(:,:,k) = dimag(f)
  enddo

return

end subroutine get_exact_null

#else
!-------------------------------------------------------------  
! Linear wave given in CQG 24S327
!-------------------------------------------------------------
subroutine get_exact_null(ex,crho,sigma,R,RJ,IJ,sst,Rmin,T,               &
                        quR1,quR2,quI1,quI2,qlR1,qlR2,qlI1,qlI2,          &
                        gR,gI,                                            &
                        dquR1,dquR2,dquI1,dquI2,                          &
                        bdquR1,bdquR2,bdquI1,bdquI2,                      &
                        dgR,dgI,bdgR,bdgI)
implicit none
! argument variables
integer, intent(in ):: ex(1:3),sst
real*8,intent(in) :: Rmin,T
double precision,intent(in),dimension(ex(1))::crho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::R
real*8,dimension(ex(1),ex(2),ex(3)),intent(out)::RJ,IJ
real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: quR1,quR2,quI1,quI2
real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: qlR1,qlR2,qlI1,qlI2
real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: gR,gI
real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dquR1,dquR2,dquI1,dquI2
real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: bdquR1,bdquR2,bdquI1,bdquI2
real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dgR,dgI,bdgR,bdgI

integer :: i,j,k
real*8 ::x,y,z,hgr,gt,gp,tgrho,tgsigma,tc,ts
double complex :: Yslm,II,Jr,swtf,ff

double complex :: beta0,C1,C2
integer :: nu,m

call initial_null_paramter(beta0,C1,C2,nu,m)

  II = dcmplx(0.d0,1.d0)

  do i=1,ex(1) 
  do j=1,ex(2)
  do k=1,ex(3)
! fake global coordinate is enough here
    hgr = 1.d0
    tgrho = dtan(crho(i))
    tgsigma = dtan(sigma(j))
    tc = dsqrt((1.d0-dsin(crho(i))*dsin(sigma(j)))/2.d0)
    ts = dsqrt((1.d0+dsin(crho(i))*dsin(sigma(j)))/2.d0)
    select case (sst)
    case (0)
      z = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = z*tgrho
      y = z*tgsigma
    case (1)
      z = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = z*tgrho
      y = z*tgsigma
    case (2)
      x = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      y = x*tgrho
      z = x*tgsigma
    case (3)
      x = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      y = x*tgrho
      z = x*tgsigma
    case (4)
      y = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = y*tgrho
      z = y*tgsigma
    case (5)
      y = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = y*tgrho
      z = y*tgsigma
    case default
      write(*,*) "get_exact_null: not recognized sst = ",sst
      return
    end select
    gt = dacos(z/hgr)
    gp = datan2(y,x)
    swtf = 2.d0*tc*ts*(ts+II*tc)/dcos(sigma(j))
    if(sst==1 .or. sst==3 .or. sst==4) swtf = dconjg(swtf)
    select case (sst)
    case (0,1)
      swtf = swtf/(dcos(gp)+II*dcos(gt)*dsin(gp))*(dcos(gt)**2+dsin(gt)**2*dcos(gp)**2)
    case (2,3)
      swtf = II*swtf*dsin(gt)
    case (4,5)
      swtf = -II*swtf*dsin(gt)
    end select

    hgr = (1.d0-R(k))/R(k)/Rmin
    Jr = (2.4d1*beta0+3.d0*II*nu*C1-II*nu**3*C2)/3.6d1+C1/4.d0*hgr-C2/1.2d1*hgr**3
    hgr = dreal(Jr*(dcos(nu*T)+II*dsin(nu*T)))
    Jr = Yslm(2,2,m,gt,gp)
    ff = dsqrt(dble((2-1)*2*(2+1)*(2+2)))*hgr*Jr*swtf**2
    RJ(i,j,k) = dreal(ff)
    IJ(i,j,k) = dimag(ff)
  enddo
  enddo
  enddo

return

end subroutine get_exact_null

#endif

#if 0
! for checking derivs_eth and dderivs_eth
!-------------------------
  subroutine get_null_boundary_c(ex,crho,sigma,R,beta,RQ,IQ,RU,IU,W,RTheta,ITheta, &
                        quR1,quR2,quI1,quI2,qlR1,qlR2,qlI1,qlI2,          &
                        gR,gI,                                            &
                        dquR1,dquR2,dquI1,dquI2,                          &
                        bdquR1,bdquR2,bdquI1,bdquI2,                      &
                        dgR,dgI,bdgR,bdgI,T,Rmin,sst)

  implicit none

!~~~~~~% Input parameters:
  integer,intent(in) :: ex(3),sst
  real*8,intent(in) :: T,Rmin
  double precision,intent(in),dimension(ex(1))::crho
  double precision,intent(in),dimension(ex(2))::sigma
  double precision,intent(in),dimension(ex(3))::R
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout) :: beta,RQ,IQ,RU,IU
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout) :: W,RTheta,ITheta
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: quR1,quR2,quI1,quI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: qlR1,qlR2,qlI1,qlI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: gR,gI
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dquR1,dquR2,dquI1,dquI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: bdquR1,bdquR2,bdquI1,bdquI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dgR,dgI,bdgR,bdgI

  double complex,dimension(ex(1),ex(2)) :: Y20,dY20,ddY20,f
  integer :: i,j,k
  real*8 ::x,y,z,hgr,gt,gp,tgrho,tgsigma
  double complex :: Yslm,II,Jr

double complex :: beta0,C1,C2
integer :: nu,m

call initial_null_paramter(beta0,C1,C2,nu,m)

  II = dcmplx(0.d0,1.d0)

!  write(*,*) abs(II)   confirms abs == cabs

  do k=1,ex(3)
  do i=1,ex(1)
  do j=1,ex(2)
! fake global coordinate is enough here  
    hgr = 1.d0
    tgrho = dtan(crho(i))
    tgsigma = dtan(sigma(j))
    select case (sst)
    case (0)
      z = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = z*tgrho
      y = z*tgsigma
    case (1)
      z = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = z*tgrho
      y = z*tgsigma
    case (2)
      x = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      y = x*tgrho
      z = x*tgsigma
    case (3)
      x = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      y = x*tgrho
      z = x*tgsigma
    case (4)
      y = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = y*tgrho
      z = y*tgsigma
    case (5)
      y = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = y*tgrho
      z = y*tgsigma
    case default
      write(*,*) "get_null_boundary_c: not recognized sst = ",sst
      return
    end select
    gt = dacos(z/hgr)
    gp = datan2(y,x)

    hgr = (1.d0-R(k))/R(k)/Rmin
    Jr = (2.4d1*beta0+3.d0*II*nu*C1-II*nu**3*C2)/3.6d1+C1/4.d0*hgr-C2/1.2d1*hgr**3
    RTheta(i,j,k) = dreal(Jr*nu*(II*dcos(nu*T)-dsin(nu*T)))
    f(i,j) = dreal(Jr*(dcos(nu*T)+II*dsin(nu*T)))

    Jr = (-2.4d1*II*nu*beta0+3.d0*nu*nu*C1-nu**4*C2)/36.d0+2.d0*beta0*hgr&
        +C1/2.d0*hgr*hgr+II*nu*C2/3.d0*hgr**3+C2/4.d0*hgr**4
    RU(i,j,k) = dreal(Jr*(dcos(nu*T)+II*dsin(nu*T)))

! Re(r^2*Ul_,r*exp(i nu T)) of CQG 24S327, (12) indeed    
    Jr = -2.d0*beta0-C1*hgr-II*nu*C2*hgr**2-C2*hgr**3
    RQ(i,j,k) = dreal(Jr*(dcos(nu*T)+II*dsin(nu*T)))

    Jr = (2.4d1*II*nu*beta0-3.d0*nu*nu*C1+nu**4*C2)/6.d0+(3.d0*II*nu*C1-6.d0*beta0-II*nu**3*C2)/3.d0*hgr&
        -nu*nu*C2*hgr*hgr+II*nu*C2*hgr**3+C2/2.d0*hgr**4
    W(i,j,k) = dreal(Jr*(dcos(nu*T)+II*dsin(nu*T)))

    Y20(i,j) = Yslm(0,2,m,gt,gp)

  enddo
  enddo

     call derivs_eth(ex(1:2),crho,sigma,Y20,dY20,0,1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call dderivs_eth(ex(1:2),crho,sigma,Y20,ddY20,0,1,1,  &
                      quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k),&
                      dquR1(:,:,k),dquR2(:,:,k),dquI1(:,:,k),dquI2(:,:,k),       &
                      bdquR1(:,:,k),bdquR2(:,:,k),bdquI1(:,:,k),bdquI2(:,:,k),   &
                      dgR(:,:,k),dgI(:,:,k),bdgR(:,:,k),bdgI(:,:,k))

     beta(:,:,k) = dreal(beta0*cdexp(II*nu*T))*dreal(Y20)
     W(:,:,k) = W(:,:,k)*dreal(Y20)

     f = dexp(-2.d0*beta(:,:,k))*(f*ddY20*RQ(:,:,k)*dconjg(dY20)+dsqrt(1.d0+abs(f*ddY20))*RQ(:,:,k)*dY20)
     RQ(:,:,k) = dreal(f)
     IQ(:,:,k) = dimag(f)

     f = ddY20*RTheta(:,:,k)
     RTheta(:,:,k) = dreal(f)
     ITheta(:,:,k) = dimag(f)

     f = dY20*RU(:,:,k)
     RU(:,:,k) = dreal(f)
     IU(:,:,k) = dimag(f)
   enddo

  return

  end subroutine get_null_boundary_c

#else

!-------------------------
  subroutine get_null_boundary_c(ex,crho,sigma,R,beta,RQ,IQ,RU,IU,W,RTheta,ITheta, &
                        quR1,quR2,quI1,quI2,qlR1,qlR2,qlI1,qlI2,          &
                        gR,gI,                                            &
                        dquR1,dquR2,dquI1,dquI2,                          &
                        bdquR1,bdquR2,bdquI1,bdquI2,                      &
                        dgR,dgI,bdgR,bdgI,T,Rmin,sst)

  implicit none

!~~~~~~% Input parameters:
  integer,intent(in) :: ex(3),sst
  real*8,intent(in) :: T,Rmin
  double precision,intent(in),dimension(ex(1))::crho
  double precision,intent(in),dimension(ex(2))::sigma
  double precision,intent(in),dimension(ex(3))::R
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout) :: beta,RQ,IQ,RU,IU
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout) :: W,RTheta,ITheta
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: quR1,quR2,quI1,quI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: qlR1,qlR2,qlI1,qlI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: gR,gI
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dquR1,dquR2,dquI1,dquI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: bdquR1,bdquR2,bdquI1,bdquI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dgR,dgI,bdgR,bdgI

integer :: i,j,k
real*8 ::x,y,z,hgr,gt,gp,tgrho,tgsigma,tc,ts,rf
double complex :: Yslm,II,Jr,swtf,ff

double complex :: beta0,C1,C2
integer :: nu,m

#if 0
real*8 :: betax,KK,KKx,Wx
double complex :: CJ,DCJ,CJx,CJxx,DCJx,CU,CUx,DCU,DCUx,bDCU,bDCUx,CB,DCB,bDCB,CBx
double complex :: Cnu,Cnux,Ck,fCTheta,fCThetax,Theta_rhs


T=0.25d0
  i=1
  j=1
  k=1
  hgr = 1.d0
  beta(i,j,k) = dreal(beta0*cdexp(II*nu*T))
  CB = beta(i,j,k)
  DCB = CB
  bDCB = CB
  betax = 0.d0
  Jr = (2.4d1*II*nu*beta0-3.d0*nu*nu*C1+nu**4*C2)/6.d0+(3.d0*II*nu*C1-6.d0*beta0-II*nu**3*C2)/3.d0/hgr&
        -nu*nu*C2/hgr/hgr+II*nu*C2/hgr**3+C2/2.d0/hgr**4
  W(i,j,k) = dreal(Jr*cdexp(II*nu*T))
  Jr = -(3.d0*II*nu*C1-6.d0*beta0-II*nu**3*C2)/3.d0/hgr**2&
        +2.d0*nu*nu*C2/hgr**3-3.d0*II*nu*C2/hgr**4-2.d0*C2/hgr**5
  Wx = dreal(Jr*cdexp(II*nu*T))*(Rmin+hgr)**2/Rmin
  Jr = (2.4d1*beta0+3.d0*II*nu*C1-II*nu**3*C2)/3.6d1+C1/4.d0/hgr-C2/1.2d1/hgr**3
  CJ = dreal(Jr*cdexp(II*nu*T))
  fCTheta = dreal(Jr*II*nu*cdexp(II*nu*T))
  DCJ=CJ
  KK = dsqrt(1.d0+cdabs(CJ)**2)
  Jr = -C1/4.d0/hgr**2+C2/4.d0/hgr**4
  rf = dreal(Jr*cdexp(II*nu*T))
  CJx = rf*(Rmin+hgr)**2/Rmin
  fCThetax = dreal(Jr*II*nu*cdexp(II*nu*T))*(Rmin+hgr)**2/Rmin
  Jr = C1/2.d0/hgr**3-C2/hgr**5
  rf = dreal(Jr*cdexp(II*nu*T))
  CJxx = rf*(Rmin+hgr)**4/Rmin**2+2.d0*(Rmin+hgr)/Rmin*CJx
  DCJx = CJx
  KKx = dreal(CJ*dconjg(CJx))/KK
  Jr = (-2.4d1*II*nu*beta0+3.d0*nu*nu*C1-nu**4*C2)/36.d0+2.d0*beta0/hgr&
        +C1/2.d0/hgr/hgr+II*nu*C2/3.d0/hgr**3+C2/4.d0/hgr**4
  CU = dreal(Jr*cdexp(II*nu*T))
  bDCU = CU
  DCU=CU
  Jr = -2.d0*beta0/hgr/hgr-C1/hgr**3-II*nu*C2/hgr**4-C2/hgr**5
  rf = dreal(Jr*cdexp(II*nu*T))
  CUx = rf*(Rmin+hgr)**2/Rmin
  DCUx=CUx
  bDCUx=CUx

  Cnu = CJ
  Cnux = CJx
  Ck = 0.d0
  hgr = hgr/(Rmin+hgr)


        call getndxs(T,crho(i),sigma(j),hgr,beta(i,j,k),KK,CU,bDCU,DCU, &
                     CB,DCB,W(i,j,k),CJ,DCJ,bDCB,Cnu,Ck,fCTheta,sst,Rmin)
        call getdxs(T,crho(i),sigma(j),hgr,betax,KKx,CUx,DCUx,bDCUx, &
                    Wx,CJx,CJxx,DCJx,Cnux,fCThetax,sst,Rmin)
!  write(*,*) 2.d0*hgr*(1.d0-hgr)*fCThetax-(-(hgr*(1-hgr)*DCUx+2.d0*DCU)+2.d0/hgr/Rmin*(1.d0-hgr)*DCB &
!         +(1.d0-hgr)**3/Rmin*(2.d0*CJx+hgr*CJxx)-2.d0*fCTheta)
!         stop
  write(*,*) fCThetax-Theta_rhs(hgr,Rmin,beta(i,j,k),betax,KK,KKx,CU,CUx,DCUx,bDCU,bDCUx, &
                       DCU,CB,DCB,W(i,j,k),Wx,CJ,DCJ,    &
                       CJx,CJxx,DCJx,bDCB,Cnu,Cnux,Ck,fCTheta)
  stop              
#endif
  call initial_null_paramter(beta0,C1,C2,nu,m)

  II = dcmplx(0.d0,1.d0)

  do i=1,ex(1)
  do j=1,ex(2)
  do k=1,ex(3)
!fake global coordinate is enough here  
    hgr = 1.d0
    tgrho = dtan(crho(i))
    tgsigma = dtan(sigma(j))
    tc = dsqrt((1.d0-dsin(crho(i))*dsin(sigma(j)))/2.d0)
    ts = dsqrt((1.d0+dsin(crho(i))*dsin(sigma(j)))/2.d0)
    select case (sst)
    case (0)
      z = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = z*tgrho
      y = z*tgsigma
    case (1)
      z = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = z*tgrho
      y = z*tgsigma
    case (2)
      x = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      y = x*tgrho
      z = x*tgsigma
    case (3)
      x = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      y = x*tgrho
      z = x*tgsigma
    case (4)
      y = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = y*tgrho
      z = y*tgsigma
    case (5)
      y = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = y*tgrho
      z = y*tgsigma
    case default
      write(*,*) "get_null_boundary: not recognized sst = ",sst
      return
    end select
    gt = dacos(z/hgr)
    gp = datan2(y,x)
    swtf = 2.d0*tc*ts*(ts+II*tc)/dcos(sigma(j))
    if(sst==1 .or. sst==3 .or. sst==4) swtf = dconjg(swtf)
    select case (sst)
    case (0,1)
      swtf = swtf/(dcos(gp)+II*dcos(gt)*dsin(gp))*(dcos(gt)**2+dsin(gt)**2*dcos(gp)**2)
    case (2,3)
      swtf = II*swtf*dsin(gt)
    case (4,5)
      swtf = -II*swtf*dsin(gt)
    end select
 
    hgr = (1.d0-R(k))/R(k)/Rmin
    Jr = (2.4d1*beta0+3.d0*II*nu*C1-II*nu**3*C2)/3.6d1+C1/4.d0*hgr-C2/1.2d1*hgr**3
    rf = dreal(Jr*nu*II*cdexp(II*nu*T))
    Jr = Yslm(2,2,m,gt,gp)*swtf**2
    ff = dsqrt(dble((2-1)*2*(2+1)*(2+2)))*rf*Jr
    RTheta(i,j,k) = dreal(ff)
    ITheta(i,j,k) = dimag(ff)

    rf = dreal(Yslm(0,2,m,gt,gp))
    beta(i,j,k) = rf*dreal(beta0*cdexp(II*nu*T))
    Jr = (2.4d1*II*nu*beta0-3.d0*nu*nu*C1+nu**4*C2)/6.d0+(3.d0*II*nu*C1-6.d0*beta0-II*nu**3*C2)/3.d0*hgr&
        -nu*nu*C2*hgr*hgr+II*nu*C2*hgr**3+C2/2.d0*hgr**4
    W(i,j,k) = rf*dreal(Jr*cdexp(II*nu*T))

    Jr = (-2.4d1*II*nu*beta0+3.d0*nu*nu*C1-nu**4*C2)/36.d0+2.d0*beta0*hgr&
        +C1/2.d0*hgr*hgr+II*nu*C2/3.d0*hgr**3+C2/4.d0*hgr**4
    rf = dreal(Jr*cdexp(II*nu*T))
    Jr = Yslm(1,2,m,gt,gp)*swtf
    ff = dsqrt(dble(2*(2+1)))*rf*Jr
    RU(i,j,k) = dreal(ff)
    IU(i,j,k) = dimag(ff)
  
    Jr = -2.d0*beta0-C1*hgr-II*nu*C2*hgr**2-C2*hgr**3
    rf = dreal(Jr*cdexp(II*nu*T))
    Jr = Yslm(1,2,m,gt,gp)*swtf
    ff = dsqrt(dble(2*(2+1)))*rf*Jr  !! U_,r
    Jr = (2.4d1*beta0+3.d0*II*nu*C1-II*nu**3*C2)/3.6d1+C1/4.d0*hgr-C2/1.2d1*hgr**3
    rf = dreal(Jr*cdexp(II*nu*T))
    Jr = Yslm(2,2,m,gt,gp)*swtf**2
    Jr = dsqrt(dble((2-1)*2*(2+1)*(2+2)))*rf*Jr  !! J
    rf = dsqrt(1.d0+cdabs(Jr)**2)  !! K
    ff = dexp(-2.d0*beta(i,j,k))*(Jr*dconjg(ff)+rf*ff)
    RQ(i,j,k) = dreal(ff)
    IQ(i,j,k) = dimag(ff)

  enddo
  enddo
  enddo

  return

  end subroutine get_null_boundary_c
#endif  

!==========================================================
subroutine initial_null_paramter(beta0,C1,C2,nu,m)

implicit none

double complex,intent(out) :: beta0,C1,C2
integer,intent(out) :: nu,m

nu=1
m=0
beta0 = dcmplx(0.d0,1.d-6)
C1    = dcmplx(3.d-6,0.d0)
C2    = dcmplx(1.d-6,0.d0)

end subroutine initial_null_paramter

#if 1
subroutine get_exact_null_theta(ex,crho,sigma,R,RTheta,ITheta,sst,Rmin,T,               &
                        quR1,quR2,quI1,quI2,qlR1,qlR2,qlI1,qlI2,          &
                        gR,gI,                                            &
                        dquR1,dquR2,dquI1,dquI2,                          &
                        bdquR1,bdquR2,bdquI1,bdquI2,                      &
                        dgR,dgI,bdgR,bdgI)
implicit none
! argument variables
integer, intent(in ):: ex(1:3),sst
real*8,intent(in) :: Rmin,T
double precision,intent(in),dimension(ex(1))::crho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::R
real*8,dimension(ex(1),ex(2),ex(3)),intent(out)::RTheta,ITheta
real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: quR1,quR2,quI1,quI2
real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: qlR1,qlR2,qlI1,qlI2
real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: gR,gI
real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dquR1,dquR2,dquI1,dquI2
real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: bdquR1,bdquR2,bdquI1,bdquI2
real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dgR,dgI,bdgR,bdgI

integer :: i,j,k
real*8 ::x,y,z,hgr,gt,gp,tgrho,tgsigma,tc,ts
double complex :: Yslm,II,Jr,swtf,ff

double complex :: beta0,C1,C2
integer :: nu,m

call initial_null_paramter(beta0,C1,C2,nu,m)

  II = dcmplx(0.d0,1.d0)

  do i=1,ex(1) 
  do j=1,ex(2)
  do k=1,ex(3)
! fake global coordinate is enough here
    hgr = 1.d0
    tgrho = dtan(crho(i))
    tgsigma = dtan(sigma(j))
    tc = dsqrt((1.d0-dsin(crho(i))*dsin(sigma(j)))/2.d0)
    ts = dsqrt((1.d0+dsin(crho(i))*dsin(sigma(j)))/2.d0)
    select case (sst)
    case (0)
      z = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = z*tgrho
      y = z*tgsigma
    case (1)
      z = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = z*tgrho
      y = z*tgsigma
    case (2)
      x = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      y = x*tgrho
      z = x*tgsigma
    case (3)
      x = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      y = x*tgrho
      z = x*tgsigma
    case (4)
      y = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = y*tgrho
      z = y*tgsigma
    case (5)
      y = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = y*tgrho
      z = y*tgsigma
    case default
      write(*,*) "get_exact_null_theta: not recognized sst = ",sst
      return
    end select
    gt = dacos(z/hgr)
    gp = datan2(y,x)
    swtf = 2.d0*tc*ts*(ts+II*tc)/dcos(sigma(j))
    if(sst==1 .or. sst==3 .or. sst==4) swtf = dconjg(swtf)
    select case (sst)
    case (0,1)
      swtf = swtf/(dcos(gp)+II*dcos(gt)*dsin(gp))*(dcos(gt)**2+dsin(gt)**2*dcos(gp)**2)
    case (2,3)
      swtf = II*swtf*dsin(gt)
    case (4,5)
      swtf = -II*swtf*dsin(gt)
    end select

    hgr = (1.d0-R(k))/R(k)/Rmin
    Jr = (2.4d1*beta0+3.d0*II*nu*C1-II*nu**3*C2)/3.6d1+C1/4.d0*hgr-C2/1.2d1*hgr**3
    hgr = dreal(Jr*nu*(-dsin(nu*T)+II*dcos(nu*T)))
    Jr = Yslm(2,2,m,gt,gp)
    ff = dsqrt(dble((2-1)*2*(2+1)*(2+2)))*hgr*Jr*swtf**2
    RTheta(i,j,k) = dreal(ff)
    ITheta(i,j,k) = dimag(ff)

  enddo
  enddo
  enddo

return

end subroutine get_exact_null_theta
!------------------------------------------------------------------------------------------------
subroutine get_exact_null_theta_x(ex,crho,sigma,R,RThetax,IThetax,sst,Rmin,T,               &
                        quR1,quR2,quI1,quI2,qlR1,qlR2,qlI1,qlI2,          &
                        gR,gI,                                            &
                        dquR1,dquR2,dquI1,dquI2,                          &
                        bdquR1,bdquR2,bdquI1,bdquI2,                      &
                        dgR,dgI,bdgR,bdgI)
implicit none
! argument variables
integer, intent(in ):: ex(1:3),sst
real*8,intent(in) :: Rmin,T
double precision,intent(in),dimension(ex(1))::crho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::R
real*8,dimension(ex(1),ex(2),ex(3)),intent(out)::RThetax,IThetax
real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: quR1,quR2,quI1,quI2
real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: qlR1,qlR2,qlI1,qlI2
real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: gR,gI
real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dquR1,dquR2,dquI1,dquI2
real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: bdquR1,bdquR2,bdquI1,bdquI2
real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dgR,dgI,bdgR,bdgI

integer :: i,j,k
real*8 ::x,y,z,hgr,gt,gp,tgrho,tgsigma,tc,ts
double complex :: Yslm,II,Jr,swtf,ff

double complex :: beta0,C1,C2
integer :: nu,m

call initial_null_paramter(beta0,C1,C2,nu,m)

  II = dcmplx(0.d0,1.d0)

  do i=1,ex(1) 
  do j=1,ex(2)
  do k=1,ex(3)
! fake global coordinate is enough here
    hgr = 1.d0
    tgrho = dtan(crho(i))
    tgsigma = dtan(sigma(j))
    tc = dsqrt((1.d0-dsin(crho(i))*dsin(sigma(j)))/2.d0)
    ts = dsqrt((1.d0+dsin(crho(i))*dsin(sigma(j)))/2.d0)
    select case (sst)
    case (0)
      z = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = z*tgrho
      y = z*tgsigma
    case (1)
      z = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = z*tgrho
      y = z*tgsigma
    case (2)
      x = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      y = x*tgrho
      z = x*tgsigma
    case (3)
      x = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      y = x*tgrho
      z = x*tgsigma
    case (4)
      y = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = y*tgrho
      z = y*tgsigma
    case (5)
      y = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = y*tgrho
      z = y*tgsigma
    case default
      write(*,*) "get_exact_null_theta_x: not recognized sst = ",sst
      return
    end select
    gt = dacos(z/hgr)
    gp = datan2(y,x)
    swtf = 2.d0*tc*ts*(ts+II*tc)/dcos(sigma(j))
    if(sst==1 .or. sst==3 .or. sst==4) swtf = dconjg(swtf)
    select case (sst)
    case (0,1)
      swtf = swtf/(dcos(gp)+II*dcos(gt)*dsin(gp))*(dcos(gt)**2+dsin(gt)**2*dcos(gp)**2)
    case (2,3)
      swtf = II*swtf*dsin(gt)
    case (4,5)
      swtf = -II*swtf*dsin(gt)
    end select

    hgr = (1.d0-R(k))/R(k)/Rmin
    Jr = C1/4.d0-C2/1.2d1*3*hgr**2
    Jr = -Jr/Rmin/R(k)/R(k)
    hgr = dreal(Jr*nu*(-dsin(nu*T)+II*dcos(nu*T)))
    Jr = Yslm(2,2,m,gt,gp)
    ff = dsqrt(dble((2-1)*2*(2+1)*(2+2)))*hgr*Jr*swtf**2
    RThetax(i,j,k) = dreal(ff)
    IThetax(i,j,k) = dimag(ff)

  enddo
  enddo
  enddo

return

end subroutine get_exact_null_theta_x
!-------------------------
  subroutine get_exact_Jul(ex,crho,sigma,R,RJul,IJul, &
                        quR1,quR2,quI1,quI2,qlR1,qlR2,qlI1,qlI2,          &
                        gR,gI,                                            &
                        dquR1,dquR2,dquI1,dquI2,                          &
                        bdquR1,bdquR2,bdquI1,bdquI2,                      &
                        dgR,dgI,bdgR,bdgI,T,Rmin,sst)

  implicit none

!~~~~~~% Input parameters:
  integer,intent(in) :: ex(3),sst
  real*8,intent(in) :: T,Rmin
  double precision,intent(in),dimension(ex(1))::crho
  double precision,intent(in),dimension(ex(2))::sigma
  double precision,intent(in),dimension(ex(3))::R
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout) :: RJul,IJul
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: quR1,quR2,quI1,quI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: qlR1,qlR2,qlI1,qlI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: gR,gI
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dquR1,dquR2,dquI1,dquI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: bdquR1,bdquR2,bdquI1,bdquI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dgR,dgI,bdgR,bdgI

integer :: i,j,k
real*8 ::x,y,z,hgr,gt,gp,tgrho,tgsigma,tc,ts,rf
double complex :: Yslm,II,Jr,swtf,ff

double complex :: beta0,C1,C2
integer :: nu,m

  call initial_null_paramter(beta0,C1,C2,nu,m)

  II = dcmplx(0.d0,1.d0)

  do i=1,ex(1)
  do j=1,ex(2)
  do k=1,ex(3)
!fake global coordinate is enough here  
    hgr = 1.d0
    tgrho = dtan(crho(i))
    tgsigma = dtan(sigma(j))
    tc = dsqrt((1.d0-dsin(crho(i))*dsin(sigma(j)))/2.d0)
    ts = dsqrt((1.d0+dsin(crho(i))*dsin(sigma(j)))/2.d0)
    select case (sst)
    case (0)
      z = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = z*tgrho
      y = z*tgsigma
    case (1)
      z = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = z*tgrho
      y = z*tgsigma
    case (2)
      x = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      y = x*tgrho
      z = x*tgsigma
    case (3)
      x = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      y = x*tgrho
      z = x*tgsigma
    case (4)
      y = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = y*tgrho
      z = y*tgsigma
    case (5)
      y = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = y*tgrho
      z = y*tgsigma
    case default
      write(*,*) "get_exact_Jul: not recognized sst = ",sst
      return
    end select
    gt = dacos(z/hgr)
    gp = datan2(y,x)
    swtf = 2.d0*tc*ts*(ts+II*tc)/dcos(sigma(j))
    if(sst==1 .or. sst==3 .or. sst==4) swtf = dconjg(swtf)
    select case (sst)
    case (0,1)
      swtf = swtf/(dcos(gp)+II*dcos(gt)*dsin(gp))*(dcos(gt)**2+dsin(gt)**2*dcos(gp)**2)
    case (2,3)
      swtf = II*swtf*dsin(gt)
    case (4,5)
      swtf = -II*swtf*dsin(gt)
    end select
 
    hgr = (1.d0-R(k))/R(k)/Rmin
    Jr = C1/4.d0-C2/4.d0*hgr**2
    rf = dreal(Jr*nu*II*cdexp(II*nu*T))
    Jr = Yslm(2,2,m,gt,gp)*swtf**2
    ff = dsqrt(dble((2-1)*2*(2+1)*(2+2)))*rf*Jr
    RJul(i,j,k) = dreal(ff)
    IJul(i,j,k) = dimag(ff)

  enddo
  enddo
  enddo

  return

  end subroutine get_exact_Jul
!-------------------------
  subroutine get_fake_Ju(ex,crho,sigma,R,RJul,IJul, &
                        quR1,quR2,quI1,quI2,qlR1,qlR2,qlI1,qlI2,          &
                        gR,gI,                                            &
                        dquR1,dquR2,dquI1,dquI2,                          &
                        bdquR1,bdquR2,bdquI1,bdquI2,                      &
                        dgR,dgI,bdgR,bdgI,T,Rmin,sst)

  implicit none

!~~~~~~% Input parameters:
  integer,intent(in) :: ex(3),sst
  real*8,intent(in) :: T,Rmin
  double precision,intent(in),dimension(ex(1))::crho
  double precision,intent(in),dimension(ex(2))::sigma
  double precision,intent(in),dimension(ex(3))::R
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout) :: RJul,IJul
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: quR1,quR2,quI1,quI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: qlR1,qlR2,qlI1,qlI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: gR,gI
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dquR1,dquR2,dquI1,dquI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: bdquR1,bdquR2,bdquI1,bdquI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dgR,dgI,bdgR,bdgI

integer :: i,j,k
real*8 ::x,y,z,hgr,gt,gp,tgrho,tgsigma,tc,ts,rf
double complex :: Yslm,II,Jr,swtf,ff

double complex :: beta0,C1,C2
integer :: nu,m

  call initial_null_paramter(beta0,C1,C2,nu,m)

  II = dcmplx(0.d0,1.d0)

  do i=1,ex(1)
  do j=1,ex(2)
  do k=1,ex(3)
!fake global coordinate is enough here  
    hgr = 1.d0
    tgrho = dtan(crho(i))
    tgsigma = dtan(sigma(j))
    tc = dsqrt((1.d0-dsin(crho(i))*dsin(sigma(j)))/2.d0)
    ts = dsqrt((1.d0+dsin(crho(i))*dsin(sigma(j)))/2.d0)
    select case (sst)
    case (0)
      z = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = z*tgrho
      y = z*tgsigma
    case (1)
      z = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = z*tgrho
      y = z*tgsigma
    case (2)
      x = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      y = x*tgrho
      z = x*tgsigma
    case (3)
      x = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      y = x*tgrho
      z = x*tgsigma
    case (4)
      y = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = y*tgrho
      z = y*tgsigma
    case (5)
      y = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = y*tgrho
      z = y*tgsigma
    case default
      write(*,*) "get_fake_Ju: not recognized sst = ",sst
      return
    end select
    gt = dacos(z/hgr)
    gp = datan2(y,x)
    swtf = 2.d0*tc*ts*(ts+II*tc)/dcos(sigma(j))
    if(sst==1 .or. sst==3 .or. sst==4) swtf = dconjg(swtf)
    select case (sst)
    case (0,1)
      swtf = swtf/(dcos(gp)+II*dcos(gt)*dsin(gp))*(dcos(gt)**2+dsin(gt)**2*dcos(gp)**2)
    case (2,3)
      swtf = II*swtf*dsin(gt)
    case (4,5)
      swtf = -II*swtf*dsin(gt)
    end select
 
    hgr = (1.d0-R(k))/R(k)/Rmin
    Jr = (2.4d1*beta0+3.d0*II*nu*C1-II*nu**3*C2)/3.6d1+C1/4.d0*hgr-C2/1.2d1*hgr**3
    rf = dreal(Jr*nu*II*cdexp(II*nu*T))
    ff = dcmplx(rf,0.d0)
    RJul(i,j,k) = dreal(ff)
    IJul(i,j,k) = dimag(ff)

  enddo
  enddo
  enddo

if(sst == 0 .and. crho(1) < -dacos(-1.d0)/4 .and. sigma(1) < -dacos(-1.d0)/4)then  
  write(*,*)"T=",T,"exp(i T)=",cdexp(II*nu*T)
  write(*,*)RJul(ex(1)/2,ex(2)/2,ex(3)),RJul(ex(1)/2,ex(2)/2,ex(3)-1),R(2)-R(1)
endif

  return

  end subroutine get_fake_Ju
!-------------------------
  subroutine get_exact_omegau(ex,crho,sigma,R,omegau, &
                        quR1,quR2,quI1,quI2,qlR1,qlR2,qlI1,qlI2,          &
                        gR,gI,                                            &
                        dquR1,dquR2,dquI1,dquI2,                          &
                        bdquR1,bdquR2,bdquI1,bdquI2,                      &
                        dgR,dgI,bdgR,bdgI,T,Rmin,sst)

  implicit none

!~~~~~~% Input parameters:
  integer,intent(in) :: ex(3),sst
  real*8,intent(in) :: T,Rmin
  double precision,intent(in),dimension(ex(1))::crho
  double precision,intent(in),dimension(ex(2))::sigma
  double precision,intent(in),dimension(ex(3))::R
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout) :: omegau
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: quR1,quR2,quI1,quI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: qlR1,qlR2,qlI1,qlI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: gR,gI
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dquR1,dquR2,dquI1,dquI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: bdquR1,bdquR2,bdquI1,bdquI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dgR,dgI,bdgR,bdgI

integer :: i,j,k
real*8 ::x,y,z,hgr,gt,gp,tgrho,tgsigma,tc,ts,rf
double complex :: Yslm,II,Jr,swtf,ff

double complex :: beta0,C1,C2
integer :: nu,m

  call initial_null_paramter(beta0,C1,C2,nu,m)

  II = dcmplx(0.d0,1.d0)

  do i=1,ex(1)
  do j=1,ex(2)
  do k=1,ex(3)
!fake global coordinate is enough here  
    hgr = 1.d0
    tgrho = dtan(crho(i))
    tgsigma = dtan(sigma(j))
    tc = dsqrt((1.d0-dsin(crho(i))*dsin(sigma(j)))/2.d0)
    ts = dsqrt((1.d0+dsin(crho(i))*dsin(sigma(j)))/2.d0)
    select case (sst)
    case (0)
      z = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = z*tgrho
      y = z*tgsigma
    case (1)
      z = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = z*tgrho
      y = z*tgsigma
    case (2)
      x = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      y = x*tgrho
      z = x*tgsigma
    case (3)
      x = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      y = x*tgrho
      z = x*tgsigma
    case (4)
      y = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = y*tgrho
      z = y*tgsigma
    case (5)
      y = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = y*tgrho
      z = y*tgsigma
    case default
      write(*,*) "get_exact_omegau: not recognized sst = ",sst
      return
    end select
    gt = dacos(z/hgr)
    gp = datan2(y,x)
    swtf = 2.d0*tc*ts*(ts+II*tc)/dcos(sigma(j))
    if(sst==1 .or. sst==3 .or. sst==4) swtf = dconjg(swtf)
    select case (sst)
    case (0,1)
      swtf = swtf/(dcos(gp)+II*dcos(gt)*dsin(gp))*(dcos(gt)**2+dsin(gt)**2*dcos(gp)**2)
    case (2,3)
      swtf = II*swtf*dsin(gt)
    case (4,5)
      swtf = -II*swtf*dsin(gt)
    end select
 
    hgr = (1.d0-R(k))/R(k)/Rmin
    Jr = (2.4d1*beta0+3.d0*II*nu*C1-II*nu**3*C2)/3.6d1+C1/4.d0*hgr-C2/1.2d1*hgr**3
    rf = dreal(Jr*nu*II*cdexp(II*nu*T))
    Jr = Yslm(0,2,m,gt,gp)
    ff = -dble(2*(2+1))/2.d0*rf*Jr
    omegau(i,j,k) = dreal(ff)

  enddo
  enddo
  enddo

  return

  end subroutine get_exact_omegau
!-------------------------
  subroutine get_exact_eth2omega(ex,crho,sigma,R,Reth2omega,Ieth2omega, &
                        quR1,quR2,quI1,quI2,qlR1,qlR2,qlI1,qlI2,          &
                        gR,gI,                                            &
                        dquR1,dquR2,dquI1,dquI2,                          &
                        bdquR1,bdquR2,bdquI1,bdquI2,                      &
                        dgR,dgI,bdgR,bdgI,T,Rmin,sst)

  implicit none

!~~~~~~% Input parameters:
  integer,intent(in) :: ex(3),sst
  real*8,intent(in) :: T,Rmin
  double precision,intent(in),dimension(ex(1))::crho
  double precision,intent(in),dimension(ex(2))::sigma
  double precision,intent(in),dimension(ex(3))::R
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout) :: Reth2omega,Ieth2omega
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: quR1,quR2,quI1,quI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: qlR1,qlR2,qlI1,qlI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: gR,gI
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dquR1,dquR2,dquI1,dquI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: bdquR1,bdquR2,bdquI1,bdquI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dgR,dgI,bdgR,bdgI

integer :: i,j,k
real*8 ::x,y,z,hgr,gt,gp,tgrho,tgsigma,tc,ts,rf
double complex :: Yslm,II,Jr,swtf,ff

double complex :: beta0,C1,C2
integer :: nu,m

  call initial_null_paramter(beta0,C1,C2,nu,m)

  II = dcmplx(0.d0,1.d0)

  do i=1,ex(1)
  do j=1,ex(2)
  do k=1,ex(3)
!fake global coordinate is enough here  
    hgr = 1.d0
    tgrho = dtan(crho(i))
    tgsigma = dtan(sigma(j))
    tc = dsqrt((1.d0-dsin(crho(i))*dsin(sigma(j)))/2.d0)
    ts = dsqrt((1.d0+dsin(crho(i))*dsin(sigma(j)))/2.d0)
    select case (sst)
    case (0)
      z = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = z*tgrho
      y = z*tgsigma
    case (1)
      z = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = z*tgrho
      y = z*tgsigma
    case (2)
      x = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      y = x*tgrho
      z = x*tgsigma
    case (3)
      x = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      y = x*tgrho
      z = x*tgsigma
    case (4)
      y = hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = y*tgrho
      z = y*tgsigma
    case (5)
      y = -hgr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      x = y*tgrho
      z = y*tgsigma
    case default
      write(*,*) "get_exact_eth2omega: not recognized sst = ",sst
      return
    end select
    gt = dacos(z/hgr)
    gp = datan2(y,x)
    swtf = 2.d0*tc*ts*(ts+II*tc)/dcos(sigma(j))
    if(sst==1 .or. sst==3 .or. sst==4) swtf = dconjg(swtf)
    select case (sst)
    case (0,1)
      swtf = swtf/(dcos(gp)+II*dcos(gt)*dsin(gp))*(dcos(gt)**2+dsin(gt)**2*dcos(gp)**2)
    case (2,3)
      swtf = II*swtf*dsin(gt)
    case (4,5)
      swtf = -II*swtf*dsin(gt)
    end select
 
    hgr = (1.d0-R(k))/R(k)/Rmin
    Jr = (2.4d1*beta0+3.d0*II*nu*C1-II*nu**3*C2)/3.6d1+C1/4.d0*hgr-C2/1.2d1*hgr**3
    rf = dreal(Jr*cdexp(II*nu*T))
    Jr = Yslm(2,2,m,gt,gp)
    ff = -dble(2*(2+1))/2.d0*dsqrt(dble((2-1)*2*(2+1)*(2+2)))*rf*Jr
    Reth2omega(i,j,k) = dreal(ff)
    Ieth2omega(i,j,k) = dimag(ff)

  enddo
  enddo
  enddo

  return

  end subroutine get_exact_eth2omega
#endif  
#endif  
