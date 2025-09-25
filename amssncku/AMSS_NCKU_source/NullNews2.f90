

#include "macrodef.fh"

!------------------------------------------------------------------------------
! input R is X indeed
! input g00 is g00/r^2 indeed
! input g0A is g0A/r^2 indeed
! input gAB is gAB/r^2 indeed
! output Gamma is Gamma of omega^2 g_{munu}/r^2 at r = infinity or to say X = 1
! ** in coordinate (u,X,x,y) **
subroutine get_christoffel(Rmin,g00,g01,g02,g03, &
                           g22,g23,g33, &
                           dgt22,dgt23,dgt33,&
                           dg22,dg23,dg33,&
                           dgx02,dgx03,dgx22,dgx23,dgx33,&
                           dgy02,dgy03,dgy22,dgy23,dgy33,&
                           omega,dtomega,dxomega,dyomega,&
                           Gamuxx,Gamuxy,Gamuyy, &
                           Gamrxx,Gamrxy,Gamryy, &
                           Gamxxx,Gamxxy,Gamxyy, &
                           Gamyxx,Gamyxy,Gamyyy)

  implicit none

  real*8,intent(in)::Rmin
  real*8,intent(in)::g00,g01,g02,g03,g22,g23,g33
  real*8,intent(in)::dgt22,dgt23,dgt33
  real*8,intent(in)::dg22,dg23,dg33
  real*8,intent(in)::dgx02,dgx03,dgx22,dgx23,dgx33
  real*8,intent(in)::dgy02,dgy03,dgy22,dgy23,dgy33
  real*8,intent(in) :: omega,dtomega,dxomega,dyomega
  real*8,intent(out) :: Gamuxx,Gamuxy,Gamuyy
  real*8,intent(out) :: Gamrxx,Gamrxy,Gamryy
  real*8,intent(out) :: Gamxxx,Gamxxy,Gamxyy
  real*8,intent(out) :: Gamyxx,Gamyxy,Gamyyy

  real*8 :: t1;
  real*8 :: t10;
  real*8 :: t11;
  real*8 :: t117;
  real*8 :: t12;
  real*8 :: t121;
  real*8 :: t138;
  real*8 :: t142;
  real*8 :: t147;
  real*8 :: t18;
  real*8 :: t184;
  real*8 :: t190;
  real*8 :: t194;
  real*8 :: t198;
  real*8 :: t2;
  real*8 :: t204;
  real*8 :: t206;
  real*8 :: t208;
  real*8 :: t214;
  real*8 :: t216;
  real*8 :: t220;
  real*8 :: t222;
  real*8 :: t227;
  real*8 :: t230;
  real*8 :: t233;
  real*8 :: t239;
  real*8 :: t24;
  real*8 :: t241;
  real*8 :: t242;
  real*8 :: t244;
  real*8 :: t249;
  real*8 :: t25;
  real*8 :: t252;
  real*8 :: t28;
  real*8 :: t29;
  real*8 :: t32;
  real*8 :: t37;
  real*8 :: t47;
  real*8 :: t53;
  real*8 :: t54;
  real*8 :: t58;
  real*8 :: t64;
  real*8 :: t65;
  real*8 :: t66;
  real*8 :: t68;
  real*8 :: t71;
  real*8 :: t72;
  real*8 :: t73;
  real*8 :: t75;
  real*8 :: t76;
  real*8 :: t77;
  real*8 :: t80;
  real*8 :: t82;
  real*8 :: t84;
  real*8 :: t85;
  real*8 :: t88;
  real*8 :: t9;
  real*8 :: t91;

    t1 = 1/g01;
    t2 = Rmin*t1;
    t9 = 1/omega;
    t10 = Rmin*t9;
    t11 = g01*omega;
    t12 = g22*g03;
    t18 = g23*g02;
    t24 = g01*g22;
    t25 = t18*dyomega;
    t28 = g23*g03;
    t29 = t28*dxomega;
    t32 = g33*g02;
    t37 = g22*g33;
    t47 = g23*g23;
    t53 = g22*g22;
    t54 = g01*t53;
    t58 = t47*dtomega;
    t64 = Rmin*dg22;
    t65 = t64*omega;
    t66 = t37*g00;
    t68 = t18*g03;
    t71 = omega*g22;
    t72 = g03*g03;
    t73 = t71*t72;
    t75 = omega*g33;
    t76 = g02*g02;
    t77 = t75*t76;
    t80 = omega*t47*g00;
    t82 = 2.0*t24*t32*dxomega-2.0*t11*t47*dgx02+t11*t47*dgt22-2.0*t54*g33*dtomega &
          +2.0*t24*t58+2.0*t54*g03*dyomega+t65*t66+2.0*t65*t68-t64*t73-t64*t77-t64*t80;
    t84 = g01*g01;
    t85 = 1/t84;
    t88 = 1/(t37-t47);
    t91 = Rmin*dg23;
    t117 = g01*g33;
    t121 = g01*t47;
    t138 = t91*omega;
    t142 = -t11*t12*dgx33+t11*t18*dgx33+2.0*t117*t18*dxomega-2.0*t121*g03*dxomega &
           -2.0*t121*g02*dyomega+t11*t47*dgt23-t11*t47*dgx03-t11*t47*dgy02+2.0*g01*t47*g23*dtomega+t138*t66+2.0*t138*t68;
    t147 = Rmin*dg33;
    t184 = g33*g33;
    t190 = g01*t184;
    t194 = t147*omega;
    t198 = -2.0*t117*t25-2.0*t117*t29-t11*t12*dgy33+t11*t18*dgy33-2.0*t11*t47*dgy03+t11*t47*dgt33-2.0*t24*t184*dtomega &
           +2.0*t117*t58+2.0*t190*g02*dxomega+t194*t66+2.0*t194*t68;
    t204 = g02*dg22*Rmin;
    t206 = omega*g23;
    t208 = g03*dg22*Rmin;
    t214 = 2.0*t24*g33*dxomega;
    t216 = t11*g23*dgy22;
    t220 = g23*dyomega;
    t222 = 2.0*t24*t220;
    t227 = t1*t88;
    t230 = g02*dg23*Rmin;
    t233 = g03*dg23*Rmin;
    t239 = 2.0*t24*g33*dyomega;
    t241 = t11*g23*dgx33;
    t242 = g23*dxomega;
    t244 = 2.0*t117*t242;
    t249 = g02*dg33*Rmin;
    t252 = g03*dg33*Rmin;
    Gamuxx = -t2*dg22/2.0;
    Gamuxy = -t2*dg23/2.0;
    Gamuyy = -t2*dg33/2.0;
    Gamrxx = t10*(-2.0*t11*t12*dgx23+t11*t12*dgy22+2.0*t11*t18*dgx23-t11*t18*dgy22+t11*t28*dgx22-t11*t32*dgx22 &
             -t11*t37*dgt22+2.0*t11*t37*dgx02-2.0*t24*t25-2.0*t24*t29+t82)*t85*t88/2.0;
    Gamrxy = t10*(-t91*t73-t91*t77-t91*t80-2.0*t24*g33*g23*dtomega-t11*t37*dgt23+t11*t37*dgx03+t11*t37*dgy02 &
             -t11*t32*dgy22+t11*t28*dgy22+2.0*t24*t28*dyomega+t142)*t85*t88/2.0;
    Gamryy = t10*(-t147*t73-t147*t77-t147*t80+2.0*t11*t37*dgy03-t11*t37*dgt33+2.0*t24*g33*g03*dyomega &
             -2.0*t11*t32*dgy23+t11*t32*dgx33+2.0*t11*t28*dgy23-t11*t28*dgx33+t198)*t85*t88/2.0;
    Gamxxx = t9*(-2.0*t11*g23*dgx23+t11*g33*dgx22+t75*t204-4.0*t121*dxomega-t206*t208+t214+t216+t222)*t227/2.0;
    Gamxxy = t9*(t11*g33*dgy22+t75*t230-t206*t233+t239-t241-t244)*t227/2.0;
    Gamxyy = t9*(-t11*g23*dgy33-t11*g33*dgx33+2.0*t11*g33*dgy23+t75*t249-2.0*t190*dxomega+2.0*t117*t220-t206*t252)*t227/2.0;
    Gamyxx = -t9*(-2.0*t11*g22*dgx23+t11*g22*dgy22+t11*g23*dgx22-2.0*t24*t242+2.0*t54*dyomega-t71*t208+t206*t204)*t227/2.0;
    Gamyxy = -(-t11*g22*dgx33-t71*t233+t206*t230-t214+t216+t222)*t9*t227/2.0;
    Gamyyy = t9*(t11*g22*dgy33-2.0*t11*g23*dgy23+t71*t252-4.0*t121*dyomega-t206*t249+t239+t241+t244)*t227/2.0;

  return

end subroutine get_christoffel
!!----------------------------------------------------------------------------------------
subroutine get_News(crho,sigma,&
                           dxxomega,dxyomega,dyyomega,&
                           omega,dtomega,dxomega,dyomega,&
                           Gamuxx,Gamuxy,Gamuyy, &
                           Gamrxx,Gamrxy,Gamryy, &
                           Gamxxx,Gamxxy,Gamxyy, &
                           Gamyxx,Gamyxy,Gamyyy,RNew,INew,sst)

  implicit none

  integer,intent(in) :: sst
  real*8,intent(in)::crho,sigma
  real*8,intent(in) :: dxxomega,dxyomega,dyyomega
  real*8,intent(in) :: omega,dtomega,dxomega,dyomega
  real*8,intent(in) :: Gamuxx,Gamuxy,Gamuyy
  real*8,intent(in) :: Gamrxx,Gamrxy,Gamryy
  real*8,intent(in) :: Gamxxx,Gamxxy,Gamxyy
  real*8,intent(in) :: Gamyxx,Gamyxy,Gamyyy

  real*8,intent(out) :: RNew,INew


  real*8 :: cs,cr,ss,sr,tc,ts
  real*8 :: WWxx,WWxy,WWyy
  real*8 :: Rmmxx,Rmmxy,Rmmyy
  real*8 :: Immxx,Immxy,Immyy

  real*8 :: gr,tgrho,tgsigma,x,y,z,gt,gp

  double complex :: swtf,II
write(*,*) Gamrxx,Gamrxy,Gamryy
  WWxx = (dxxomega-(Gamuxx*dtomega+Gamxxx*dxomega+Gamyxx*dyomega))/omega/2
  WWxy = (dxyomega-(Gamuxy*dtomega+Gamxxy*dxomega+Gamyxy*dyomega))/omega/2
  WWyy = (dyyomega-(Gamuyy*dtomega+Gamxyy*dxomega+Gamyyy*dyomega))/omega/2

  cs = dcos(sigma)
  cr = dcos(crho)
  ss = dsin(sigma)
  sr = dsin(crho)
  tc = dsqrt((1-sr*ss)/2)
  ts = dsqrt((1+sr*ss)/2)
  Rmmxx = 4*tc*tc*ts*ts*(ts*ts-tc*tc)/cs/cs
  Rmmxy = 4*tc*tc*ts*ts*(ts*ts+tc*tc)/cs/cr
  Rmmyy = 4*tc*tc*ts*ts*(ts*ts-tc*tc)/cr/cr
  Immxx = 8*tc*tc*ts*ts*ts*tc/cs/cs
  Immxy = 0
  Immyy = -8*tc*tc*ts*ts*ts*tc/cr/cr
    
  if(sst==1 .or. sst==3 .or. sst==4)then
    Immxx = -Immxx
    Immxy = -Immxy
    Immyy = -Immyy
  endif

  RNew = Rmmxx*WWxx+2*Rmmxy*WWxy+Rmmyy*WWyy
  INew = Immxx*WWxx+2*Immxy*WWxy+Immyy*WWyy
!! change to tetrad theta phi
!fake global coordinate is enough here  

    II = dcmplx(0.d0,1.d0)
    gr = 1.d0
    tgrho = dtan(crho)
    tgsigma = dtan(sigma)
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
      write(*,*) "get_News: not recognized sst = ",sst
      return
    end select
    gt = dacos(z/gr)
    gp = datan2(y,x)
    swtf = 2.d0*tc*ts*(ts+II*tc)/dcos(sigma)
    if(sst==1 .or. sst==3 .or. sst==4) swtf = dconjg(swtf)
    select case (sst)
    case (0,1)
      swtf = swtf/(dcos(gp)+II*dcos(gt)*dsin(gp))*(dcos(gt)**2+dsin(gt)**2*dcos(gp)**2)
    case (2,3)
      swtf = II*swtf*dsin(gt)
    case (4,5)
      swtf = -II*swtf*dsin(gt)
    end select

    swtf = (RNew+II*INew)/swtf**2

    RNew = dreal(swtf)
    INew = dimag(swtf)

  return

end subroutine get_News  
!------------------------------------------------------------------------------------------------------------
subroutine get_null_news2(ex,crho,sigma,R,omega,dtomega, &
                         g00,g01,g02,g03,g22,g23,g33, &
                         dtg22,dtg23,dtg33, &
                         RNews,INews,Rmin,sst)

implicit none

integer,intent(in) :: ex(3),sst
real*8,intent(in) :: Rmin
real*8,intent(in),dimension(ex(1))::crho
real*8,intent(in),dimension(ex(2))::sigma
real*8,intent(in),dimension(ex(3))::R
real*8,dimension(ex(1),ex(2),ex(3)),intent(in ) :: omega,dtomega
real*8,dimension(ex(1),ex(2),ex(3)),intent(in ) :: g00,g01,g02,g03,g22,g23,g33
real*8,dimension(ex(1),ex(2),ex(3)),intent(in ) :: dtg22,dtg23,dtg33
real*8,dimension(ex(1),ex(2),ex(3)),intent(out) :: RNews,INews

real*8 :: Gamuxx,Gamuxy,Gamuyy
real*8 :: Gamrxx,Gamrxy,Gamryy
real*8 :: Gamxxx,Gamxxy,Gamxyy
real*8 :: Gamyxx,Gamyxy,Gamyyy
real*8 :: dg22,dg23,dg33
real*8 :: dgx22,dgx23,dgx33
real*8 :: dgx02,dgx03
real*8 :: dgy22,dgy23,dgy33
real*8 :: dgy02,dgy03
real*8 :: dxomega,dyomega
real*8 :: dxxomega,dxyomega,dyyomega

integer :: i,j,k

k = ex(3)
do i=1,ex(1)
do j=1,ex(2)
    call rderivs_x_point(ex(3),R,g22(i,j,:),dg22,k)
    call rderivs_x_point(ex(3),R,g23(i,j,:),dg23,k)
    call rderivs_x_point(ex(3),R,g33(i,j,:),dg33,k)

    call rderivs_x_point(ex(1),crho,g02(:,j,k),dgx02,i)
    call rderivs_x_point(ex(1),crho,g03(:,j,k),dgx03,i)
    call rderivs_x_point(ex(1),crho,g22(:,j,k),dgx22,i)
    call rderivs_x_point(ex(1),crho,g23(:,j,k),dgx23,i)
    call rderivs_x_point(ex(1),crho,g33(:,j,k),dgx33,i)
    call rderivs_x_point(ex(1),crho,omega(:,j,k),dxomega,i)

    call rderivs_x_point(ex(2),sigma,g02(i,:,k),dgy02,j)
    call rderivs_x_point(ex(2),sigma,g03(i,:,k),dgy03,j)
    call rderivs_x_point(ex(2),sigma,g22(i,:,k),dgy22,j)
    call rderivs_x_point(ex(2),sigma,g23(i,:,k),dgy23,j)
    call rderivs_x_point(ex(2),sigma,g33(i,:,k),dgy33,j)
    call rderivs_x_point(ex(2),sigma,omega(i,:,k),dyomega,j)

    call get_christoffel(Rmin,g00(i,j,k),g01(i,j,k),g02(i,j,k),g03(i,j,k), &
                           g22(i,j,k),g23(i,j,k),g33(i,j,k), &
                           dtg22(i,j,k),dtg23(i,j,k),dtg33(i,j,k),&
                           dg22,dg23,dg33,&
                           dgx02,dgx03,dgx22,dgx23,dgx33,&
                           dgy02,dgy03,dgy22,dgy23,dgy33,&
                           omega(i,j,k),dtomega(i,j,k),dxomega,dyomega,&
                           Gamuxx,Gamuxy,Gamuyy, &
                           Gamrxx,Gamrxy,Gamryy, &
                           Gamxxx,Gamxxy,Gamxyy, &
                           Gamyxx,Gamyxy,Gamyyy)

    call rdderivs_x_point(ex(1),crho,omega(:,j,k),dxxomega,i)
    call rdderivs_x_point(ex(2),crho,omega(i,:,k),dyyomega,j)
    call rdderivs_xy_point(ex(1),ex(2),crho,sigma,omega(:,:,k),dxyomega,i,j)

    call get_News(crho(i),sigma(j),&
                           dxxomega,dxyomega,dyyomega,&
                           omega(i,j,k),dtomega(i,j,k),dxomega,dyomega,&
                           Gamuxx,Gamuxy,Gamuyy, &
                           Gamrxx,Gamrxy,Gamryy, &
                           Gamxxx,Gamxxy,Gamxyy, &
                           Gamyxx,Gamyxy,Gamyyy,RNews(i,j,k),INews(i,j,k),sst)
enddo
enddo

    return
        
end subroutine get_null_news2
!!------------------------------------------------------------------------------------------------------------
!! input g_AB and Theta_AB are divided by r^2 indeed
!! input g_00 is also divided by r^2 indeed
! the output g00 is K
subroutine get_omega_and_dtomega_pre(ex,crho,sigma,X,g22,g23,g33, &
                         omega,dtomega, Rmin)
implicit none
! argument variables
integer, intent(in ):: ex(1:3)
real*8,intent(in) :: Rmin
double precision,intent(in),dimension(ex(1))::crho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::X
real*8,dimension(ex(1),ex(2),ex(3)),intent(in)::g22,g23,g33
real*8,dimension(ex(1),ex(2),ex(3)),intent(out)::omega,dtomega


double precision,dimension(ex(3))::R
real*8,dimension(ex(1),ex(2),ex(3))::det,gup22,gup23,gup33,KK

real*8 :: sr,ss,cr,cs,sr2,ss2,cr2,cs2,tg22,tg23,tg33
real*8 :: fr,fs,frr,fss,frs,covf

integer :: i,j,k

real*8 :: m0,Pp0,Pm0,ap,am,bp,bm,cp,cm,gam

call get_RT_parameters(m0,Pp0,Pm0,ap,am,bp,bm,cp,cm,gam)

R = X*Rmin/(1-X)
det = g22*g33-g23*g23
gup22 = g33/det
gup23 = -g23/det
gup33 = g22/det

  do i=1,ex(1)
  do j=1,ex(2)
  do k=1,ex(3)
     sr = dsin(crho(i))
     ss = dsin(sigma(j))
     cr = dcos(crho(i))
     cs = dcos(sigma(j))
     sr2 = sr*sr
     ss2 = ss*ss
     cr2 = cr*cr
     cs2 = cs*cs

     tg22 = 1-sr2*ss2
     tg22 = 1/tg22/tg22

     tg23 = -sr*cr*ss*cs*tg22
     tg33 = cr2*tg22
     tg22 = cs2*tg22

! ghat/(g/r^4) indeed
     det(i,j,k) = (tg22*tg33-tg23*tg23)/det(i,j,k)
   enddo
   enddo
   enddo

   omega = dsqrt(det)
   k = ex(3)
  do i=1,ex(1)
  do j=1,ex(2)
    
     call rderivs_x_point(ex(1),crho,det(:,j,k),fr,i)
     call rderivs_x_point(ex(2),sigma,det(i,:,k),fs,j)
     call rdderivs_xy_point(ex(1),ex(2),crho,sigma,det(:,:,k),frs,i,j)
     call rdderivs_x_point(ex(1),crho,det(:,j,k),frr,i)
     call rdderivs_x_point(ex(2),sigma,det(i,:,k),fss,j)

     call std_covdiff(crho(i),sigma(j),fs,fr,fss,frr,frs,covf)

     KK(i,j,k) = dsqrt(det(i,j,k))*(1-0.25*covf/R(k)**2)
   enddo
   enddo

   dtomega = KK

   return

end subroutine get_omega_and_dtomega_pre
!------------------------------------------------------------------------------------------------------
subroutine get_dtomega(ex,crho,sigma,X,g22,g23,g33, &
                         omega,dtomega, Rmin)
implicit none
! argument variables
integer, intent(in ):: ex(1:3)
real*8,intent(in) :: Rmin
double precision,intent(in),dimension(ex(1))::crho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::X
real*8,dimension(ex(1),ex(2),ex(3)),intent(in)::omega,g22,g23,g33
real*8,dimension(ex(1),ex(2),ex(3)),intent(inout)::dtomega


double precision,dimension(ex(3))::R
real*8,dimension(ex(1),ex(2),ex(3))::det,gup22,gup23,gup33,KK

real*8 :: sr,ss,cr,cs,sr2,ss2,cr2,cs2,tg22,tg23,tg33
real*8 :: fr,fs,frr,fss,frs,covf

integer :: i,j,k

real*8 :: m0,Pp0,Pm0,ap,am,bp,bm,cp,cm,gam

call get_RT_parameters(m0,Pp0,Pm0,ap,am,bp,bm,cp,cm,gam)

   KK = dtomega

   k = ex(3)
  do i=1,ex(1)
  do j=1,ex(2)
    
     call rderivs_x_point(ex(1),crho,KK(:,j,k),fr,i)
     call rderivs_x_point(ex(2),sigma,KK(i,:,k),fs,j)
     call rdderivs_xy_point(ex(1),ex(2),crho,sigma,KK(:,:,k),frs,i,j)
     call rdderivs_x_point(ex(1),crho,KK(:,j,k),frr,i)
     call rdderivs_x_point(ex(2),sigma,KK(i,:,k),fss,j)

     call std_covdiff(crho(i),sigma(j),fs,fr,fss,frr,frs,covf)

     dtomega(i,j,k) = -covf*omega(i,j,k)**3/6/m0/2
   enddo
   enddo

   return

end subroutine get_dtomega
!!------------------------------------------------------------------------------------------------------------
!! input g_AB and Theta_AB are divided by r^2 indeed
!! input g_00 is also divided by r^2 indeed
subroutine get_omega_and_dtomega_LN(time,ex,crho,sigma,XX, &
                         omega,dtomega, Rmin,sst)
implicit none
! argument variables
integer, intent(in ):: ex(1:3),sst
real*8,intent(in) :: time,Rmin
double precision,intent(in),dimension(ex(1))::crho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::XX
real*8,dimension(ex(1),ex(2),ex(3)),intent(out)::omega,dtomega

integer :: i,j,k
real*8 :: gr,gt,gp,tgrho,tgsigma,tc,ts,x,y,z

double complex :: II,Jr,Jrt
double complex :: Zslm,z020

double complex :: beta0,C1,C2,mx,my,mlx,mly
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
      write(*,*) "get_null_boundary3: not recognized sst = ",sst
      return
    end select
    gt = dacos(z/gr)
    gp = datan2(y,x)

    z020 = Zslm(0,2,m,gt,gp)

    Jr = (2.4d1*beta0+3.d0*II*nu*C1-II*nu**3*C2)/3.6d1
    Jr = Jr*exp(II*nu*time)
    Jrt = II*nu*Jr*exp(II*nu*time)

    Jr = dsqrt(dble((2-1)))*dreal(Jr)*z020
    Jrt = dsqrt(dble((2-1)))*dreal(Jrt)*z020

    omega(i,j,k) = 1-dreal(Jr)
    dtomega(i,j,k) = -dreal(Jrt)

  enddo
  enddo
  enddo

   return

end subroutine get_omega_and_dtomega_LN
