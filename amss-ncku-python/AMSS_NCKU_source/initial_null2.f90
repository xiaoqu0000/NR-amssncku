

#include "macrodef.fh"

subroutine get_RT_parameters(m0o,Pp0o,Pm0o,apo,amo,bpo,bmo,cpo,cmo,gamo)
implicit none
real*8,intent(out) :: m0o,Pp0o,Pm0o,apo,amo,bpo,bmo,cpo,cmo,gamo

real*8,parameter::m0=1.d0,Pp0=1.d0,Pm0=1.d0,ap=1.d0,am=1.d0
real*8,parameter::bp=0.d0,bm=0.d0,cp=0.d0,cm=0.d0
real*8,parameter::gam=0.5d0

m0o = m0
Pp0o = Pp0
Pm0o = Pm0
apo = ap
amo = am
bpo = bp
bmo = bm
cpo = cp
cmo = cm
gamo = gam
end subroutine get_RT_parameters
!!!---------------------------------------------------------------------------------------------
  function boostbhP(P0,gam,a,b,c,gt,gp)  result(gont)
  implicit none

!~~~~~~> Input parameters:

  real*8, intent(in ):: P0,gam,a,b,c,gt,gp

  real*8::gont

  gont = dcosh(gam)+a*dsinh(gam)*dcos(gt)+dsinh(gam)*dsin(gt)*(b*dcos(gp)+c*dsin(gp))

  gont = P0*gont

  end function boostbhP
!!!!-------------------------------------------------------------------------------------------
#if 1
!! RT ID
subroutine get_initial_null2(ex,crho,sigma,XX,g22,g23,g33,sst,Rmin)
implicit none
! argument variables
integer, intent(in ):: ex(1:3),sst
real*8,intent(in) :: Rmin
double precision,intent(in),dimension(ex(1))::crho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::XX
real*8,dimension(ex(1),ex(2),ex(3)),intent(out)::g22,g23,g33

double precision,dimension(ex(3))::R
real*8 :: sr,ss,cr,cs
real*8 :: sr2,ss2,cr2,cs2
integer :: i,j,k
real*8 :: ggr,tgrho,tgsigma
real*8 ::x,y,z,gr,gt,gp
real*8,dimension(ex(1),ex(2),ex(3))::P

real*8 :: PI

real*8 :: m0,Pp0,Pm0,ap,am,bp,bm,cp,cm,gam

real*8::boostbhP

call get_RT_parameters(m0,Pp0,Pm0,ap,am,bp,bm,cp,cm,gam)

R = XX*Rmin/(1-XX)

PI = dacos(-1.d0)

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

     g22(i,j,k) = 1-sr2*ss2
     g22(i,j,k) = 1/g22(i,j,k)/g22(i,j,k)

     g23(i,j,k) = -sr*cr*ss*cs*g22(i,j,k)
     g33(i,j,k) = cr2*g22(i,j,k)
     g22(i,j,k) = cs2*g22(i,j,k)

! we want g_AB/r^2 instead of g_AB     
!     g22(i,j,k) = R(k)*R(k)*g22(i,j,k)
!     g23(i,j,k) = R(k)*R(k)*g23(i,j,k)
!     g33(i,j,k) = R(k)*R(k)*g33(i,j,k)

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
      write(*,*) "get_initial_null2: not recognized sst = ",sst
      return
    end select
    gt = dacos(z/gr)
    gp = datan2(y,x)

    P(i,j,k) = 1/(1/dsqrt(boostbhP(Pp0,gam,ap,bp,cp,gt,gp))+1/dsqrt(boostbhP(Pm0,gam,am,bm,cm,gt,gp)))**2

  enddo
  enddo
  enddo

  g22 = g22/P**2
  g23 = g23/P**2
  g33 = g33/P**2

return

end subroutine get_initial_null2
#else
!! fake RT for test
subroutine get_initial_null2(ex,crho,sigma,XX,g22,g23,g33,sst,Rmin)
implicit none
! argument variables
integer, intent(in ):: ex(1:3),sst
real*8,intent(in) :: Rmin
double precision,intent(in),dimension(ex(1))::crho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::XX
real*8,dimension(ex(1),ex(2),ex(3)),intent(out)::g22,g23,g33

double precision,dimension(ex(3))::R
real*8 :: sr,ss,cr,cs
real*8 :: sr2,ss2,cr2,cs2
integer :: i,j,k
real*8 :: ggr,tgrho,tgsigma
real*8 ::x,y,z,gr,gt,gp
real*8,dimension(ex(1),ex(2),ex(3))::P

real*8 :: PI

real*8 :: m0,Pp0,Pm0,ap,am,bp,bm,cp,cm,gam

real*8::boostbhP

call get_RT_parameters(m0,Pp0,Pm0,ap,am,bp,bm,cp,cm,gam)

R = XX*Rmin/(1-XX)

PI = dacos(-1.d0)

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

     g22(i,j,k) = 1-sr2*ss2
     g22(i,j,k) = 1/g22(i,j,k)/g22(i,j,k)

     g23(i,j,k) = -sr*cr*ss*cs*g22(i,j,k)
     g33(i,j,k) = cr2*g22(i,j,k)
     g22(i,j,k) = cs2*g22(i,j,k)

! we want g_AB/r^2 instead of g_AB     
!     g22(i,j,k) = R(k)*R(k)*g22(i,j,k)
!     g23(i,j,k) = R(k)*R(k)*g23(i,j,k)
!     g33(i,j,k) = R(k)*R(k)*g33(i,j,k)

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
      write(*,*) "get_initial_null2: not recognized sst = ",sst
      return
    end select
    gt = dacos(z/gr)
    gp = datan2(y,x)

    P(i,j,k) = 1/(1/dsqrt(boostbhP(Pp0,gam,ap,bp,cp,gt,gp))+1/dsqrt(boostbhP(Pm0,gam,am,bm,cm,gt,gp)))**2

  enddo
  enddo
  enddo

  g22 = P

return

end subroutine get_initial_null2
#endif
!!------------------------------------------------------------------------------------------------------------
subroutine std_covdiff(rho,sigma,fs,fr,fss,frr,frs,covf)
implicit none
! argument variables
real*8,intent(in) :: rho,sigma,fs,fr,fss,frr,frs
real*8,intent(out):: covf

real*8 :: t1,t2,t3,t4,t5,t6,t7,t8,t11,t12,t13,t15,t16,t19,t20
real*8 :: t27,t28,t29,t32,t33,t34,t38,t39,t51,t54,t55,t58,t59
real*8 :: t62,t71,t72,t88,t90,t91,t92,t93,t94,t95,t97,t98,t99
real*8 :: t100,t104,t107,t108,t109,t112,t113,t117,t118,t121,t128,t132,t133,t136,t137,t140,t141
real*8 :: t144,t152,t153,t154,t155,t160,t166,t169,t172,t175,t178,t181,t187,t199,t204,t205,t208
real*8 :: t209,t216,t217,t223,t226,t227,t243,t250,t256,t267,t276,t284,t287,t290,t301,t303,t306
real*8 :: t307,t310,t313,t314,t316,t319,t323,t326,t329,t338,t346,t356,t359,t368,t371,t376,t377
real*8 :: t380,t385,t387,t391,t394,t398,t401,t404,t407,t412,t415,t420,t427,t450,t451,t456,t459
real*8 :: t486,t487,t511,t516,t522,t532,t537,t546,t575,t586,t591,t595,t599,t295,t298

      t1 = cos(sigma);
      t2 = t1*t1;
      t3 = t2*t2;
      t4 = t3*t2;
      t5 = t4*fss;
      t6 = 2.0*sigma;
      t7 = cos(t6);
      t8 = t7*t7;
      t11 = cos(rho);
      t12 = t11*t11;
      t13 = t12*t11;
      t15 = sin(rho);
      t16 = fr*t15;
      t19 = t12*t12;
      t20 = t19*t13;
      t27 = t19*t11;
      t28 = t27*fr;
      t29 = t15*t8;
      t32 = 2.0*rho;
      t33 = cos(t32);
      t34 = t33*t33;
      t38 = t11*fr;
      t39 = t15*t3;
      t51 = t19*frr;
      t54 = t19*t12;
      t55 = t54*frr;
      t58 = -2.0*t5*t8-8.0*t2*t13*t16-64.0*t3*t20*t16+32.0*t4*t20*t16+4.0*t28*t29 &
            +4.0*t28*t15*t34-4.0*t38*t39+8.0*t3*t13*t16+32.0*t4*t13*t16+8.0*t2*t27*t16 &
            +4.0*t51*t2-2.0*t55*t34;
      t59 = t3*fss;
      t62 = t12*frr;
      t71 = t19*t19;
      t72 = t71*frr;
      t88 = -32.0*t59*t54+2.0*t62*t3-2.0*t55*t8+64.0*t55*t4-2.0*t5*t34-32.0*t72*t2 &
            +64.0*t72*t3-32.0*t72*t4-4.0*t55*t2-4.0*t51*t3-62.0*t55*t3+60.0*t3*t27*t16;
      t90 = sin(t32);
      t91 = sin(t6);
      t92 = t90*t91;
      t93 = t92*frs;
      t94 = t3*t8;
      t95 = t94*t34;
      t97 = t3*t1;
      t98 = t97*fs;
      t99 = sin(sigma);
      t100 = t98*t99;
      t104 = t2*fss;
      t107 = t54*fr;
      t108 = t90*t33;
      t109 = t108*t2;
      t112 = t19*t8;
      t113 = t112*t34;
      t117 = t12*fr;
      t118 = t108*t3;
      t121 = t12*t3;
      t128 = t19*t2;
      t132 = t3*t3;
      t133 = t132*fss;
      t136 = t93*t95-4.0*t100-32.0*t51*t4+2.0*t104*t19+8.0*t107*t109+t93*t113-62.0*t5*t19 &
            -4.0*t117*t118+2.0*t93*t121*t8+32.0*t2*t20*t16+2.0*t93*t128*t8-32.0*t133*t12;
      t137 = t3*t19;
      t140 = t2*t8;
      t141 = t140*t34;
      t144 = t8*t34;
      t152 = t107*t91;
      t153 = t90*t99;
      t154 = t2*t1;
      t155 = t153*t154;
      t160 = t33*t3*t8;
      t166 = frs*t19;
      t169 = t19*fr;
      t172 = frs*t3;
      t175 = t107*t90;
      t178 = -t93*t137*t8-4.0*t55*t141-2.0*t93*t128*t144+2.0*t62*t95+t93*t137*t144+16.0*t152*t155 &
             +4.0*t117*t90*t160+4.0*t107*t108*t8-t92*t166*t8+8.0*t169*t118-t92*t172*t34+4.0*t175*t160;
      t181 = t169*t90;
      t187 = t33*t2*t8;
      t199 = frs*t2;
      t204 = fs*t3*t154;
      t205 = t19*t99;
      t208 = fs*t154;
      t209 = t54*t99;
      t216 = -8.0*t181*t160-4.0*t107*t118-8.0*t175*t187-t93*t137*t34+4.0*t51*t141+2.0*t93*t128*t34 &
             -4.0*t51*t95+2.0*t92*t199*t12-64.0*t204*t205+32.0*t208*t209-64.0*t98*t209+32.0*t204*t209;
      t217 = t99*t8;
      t223 = t1*fs;
      t226 = t4*fs;
      t227 = t91*t7;
      t243 = t12*t99;
      t250 = 4.0*t98*t217+4.0*t98*t99*t34-4.0*t223*t205-4.0*t226*t227-64.0*t4*t27*t16+2.0*t93*t121*t34 &
            -t92*t166*t34-2.0*t93*t121*t144+8.0*t208*t205+8.0*t98*t243+60.0*t98*t205+32.0*t204*t243;
      t256 = t2*t34;
      t267 = t3*t34;
      t276 = -8.0*t208*t243+t92*t172+t92*t166-4.0*t51*t256-4.0*t51*t140+4.0*t51*t94+2.0*t55*t144 &
             -2.0*t62*t94-2.0*t62*t267-2.0*t55*t94+4.0*t55*t140+4.0*t51*t267;
      t284 = fs*t91*t33;
      t287 = t243*t34;
      t290 = t205*t34;
      t295 = fs*t27*t15;
      t298 = t92*t4;
      t301 = t92*t3;
      t303 = fs*t13*t15;
      t306 = t208*t99;
      t307 = t144*t12;
      t310 = t227*t19;
      t313 = t3*fs;
      t314 = t313*t91;
      t316 = t7*t19*t34;
      t319 = 4.0*t55*t256-2.0*t55*t267+2.0*t55*t95-32.0*t137*t284-8.0*t98*t287+4.0*t98*t290 &
            -8.0*t92*t2*t295-8.0*t298*t295-16.0*t301*t303-8.0*t306*t307-4.0*t226*t310-8.0*t314*t316;
      t323 = t217*t12;
      t326 = t227*t12;
      t329 = t226*t91;
      t338 = t7*t12*t34;
      t346 = t2*fs;
      t356 = 8.0*t208*t323+8.0*t226*t326+4.0*t329*t316+8.0*t208*t287+4.0*t226*t227*t34+8.0*t314*t338 &
           -16.0*t92*t199*t54-8.0*t329*t338-4.0*t346*t310+8.0*t313*t310-8.0*t313*t326+4.0*t346*t91*t316;
      t359 = t205*t8;
      t368 = t153*t97;
      t371 = t169*t91;
      t376 = t13*fr;
      t377 = t29*t2;
      t380 = t376*t15;
      t385 = t4*t12;
      t387 = fr*t7*t90;
      t391 = t15*t2*t34;
      t394 = 8.0*t181*t187+4.0*t98*t359-8.0*t169*t109-8.0*t152*t153*t1-8.0*t117*t91*t368+16.0*t371*t368 &
            -8.0*t152*t368+8.0*t376*t377-8.0*t380*t141-16.0*t371*t155-16.0*t385*t387+8.0*t376*t391;
      t398 = t4*t54;
      t401 = t3*t54;
      t404 = t2*t54;
      t407 = t4*t19;
      t412 = t39*t8;
      t415 = t28*t15;
      t420 = t39*t34;
      t427 = -32.0*t137*t387-16.0*t398*t387+32.0*t401*t387-16.0*t404*t387+32.0*t407*t387-8.0*t28*t377 &
             -8.0*t376*t412-4.0*t415*t95+4.0*t28*t412+4.0*t38*t420+4.0*t28*t420-8.0*t28*t391;
      t450 = t2*t12;
      t451 = t450*t34;
      t456 = -8.0*t376*t420+8.0*t415*t141+8.0*t380*t95-4.0*t28*t29*t34+4.0*t38*t412-8.0*t208*t290 &
             +32.0*t92*t172*t54+32.0*t401*t284-4.0*t100*t113+4.0*t223*t290-2.0*t93*t451-16.0*t398*t284;
      t459 = frs*t4;
      t486 = -16.0*t92*t459*t12+32.0*t407*t284-4.0*t98*t217*t34+2.0*t55-16.0*t385*t284-16.0*t404*t284 &
             -2.0*t104*t112-8.0*t208*t359-8.0*t98*t323-t92*t172*t8+4.0*t223*t359+32.0*t92*t459*t19;
      t487 = t19*t34;
      t511 = t140*t12;
      t516 = 4.0*t59*t487+8.0*t306*t113+8.0*t100*t307-2.0*t5*t487-4.0*t223*t99*t113+16.0*t298*t303 &
            -8.0*t298*fs*t11*t15+16.0*t301*t295+2.0*t104*t113+4.0*t59*t307-2.0*t93*t511-4.0*t59*t113;
      t522 = t8*t12;
      t532 = t144*t450;
      t537 = t12*t34;
      t546 = 2.0*t5*t113-4.0*t5*t307-4.0*t59*t522+2.0*t5-31.0*t92*t172*t19-2.0*t92*t199*t19+2.0*t93*t532 &
            -2.0*t5*t112+4.0*t5*t537-4.0*t59*t537+2.0*t5*t144+4.0*t5*t522;
      t575 = 4.0*t59*t112-2.0*t104*t487-4.0*t107*t108-4.0*t38*t15*t95-16.0*t92*t459*t54-2.0*t92*t172*t12 &
            -4.0*t5*t12+4.0*t59*t12+64.0*t5*t54+64.0*t133*t19-32.0*t133*t54-4.0*t59*t19-4.0*t415;
      t586 = t34*t34;
      t591 = t8*t8;
      t595 = 256.0*t137-32.0*t450+32.0*t451+32.0*t511-32.0*t532+1.0-2.0*t34+t586-2.0*t8+4.0*t144 &
            -2.0*t8*t586+t591-2.0*t591*t34+t591*t586;
      covf = -8.0*(t58+t88+t136+t178+t216+t250+t276+t319+t356+t394+t427+t456+t486+t516+t546+t575)/t595;

return

end subroutine std_covdiff
!!------------------------------------------------------------------------------------------------------------
!! input g_AB and Theta_AB are divided by r^2 indeed
!! input g_00 is also divided by r^2 indeed
! the output g00 is K
#if 1
subroutine get_gauge_g00_K(ex,crho,sigma,X,g22,g23,g33, &
                         Theta22,Theta23,Theta33, g00, Rmin)
implicit none
! argument variables
integer, intent(in ):: ex(1:3)
real*8,intent(in) :: Rmin
double precision,intent(in),dimension(ex(1))::crho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::X
real*8,dimension(ex(1),ex(2),ex(3)),intent(in)::g22,g23,g33
real*8,dimension(ex(1),ex(2),ex(3)),intent(in)::Theta22,Theta23,Theta33
real*8,dimension(ex(1),ex(2),ex(3)),intent(out)::g00


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

  do i=1,ex(1)
  do j=1,ex(2)
  do k=1,ex(3)
    
     call rderivs_x_point(ex(1),crho,det(:,j,k),fr,i)
     call rderivs_x_point(ex(2),sigma,det(i,:,k),fs,j)
     call rdderivs_xy_point(ex(1),ex(2),crho,sigma,det(:,:,k),frs,i,j)
     call rdderivs_x_point(ex(1),crho,det(:,j,k),frr,i)
     call rdderivs_x_point(ex(2),sigma,det(i,:,k),fss,j)

     call std_covdiff(crho(i),sigma(j),fs,fr,fss,frr,frs,covf)

     KK(i,j,k) = dsqrt(det(i,j,k))*(1-0.25*covf/R(k)**2)
   enddo
   enddo
   enddo

   g00 = KK

   return

end subroutine get_gauge_g00_K
! the input g00 is K
subroutine get_gauge_g00(ex,crho,sigma,X,g22,g23,g33, &
                         Theta22,Theta23,Theta33, g00, Rmin,fp)
implicit none
! argument variables
integer, intent(in ):: ex(1:3),fp
real*8,intent(in) :: Rmin
double precision,intent(in),dimension(ex(1))::crho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::X
real*8,dimension(ex(1),ex(2),ex(3)),intent(in)::g22,g23,g33
real*8,dimension(ex(1),ex(2),ex(3)),intent(out)::Theta22,Theta23,Theta33
real*8,dimension(ex(1),ex(2),ex(3)),intent(in)::g00


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

     Theta22(i,j,k) = tg22/6/m0
     Theta23(i,j,k) = tg23/6/m0
     Theta33(i,j,k) = tg33/6/m0
   enddo
   enddo
   enddo

  KK = g00

  if(fp == 0)then
  k = 1
  do i=1,ex(1)
  do j=1,ex(2)
    
     call rderivs_x_point(ex(1),crho,KK(:,j,k),fr,i)
     call rderivs_x_point(ex(2),sigma,KK(i,:,k),fs,j)
     call rdderivs_xy_point(ex(1),ex(2),crho,sigma,KK(:,:,k),frs,i,j)
     call rdderivs_x_point(ex(1),crho,KK(:,j,k),frr,i)
     call rdderivs_x_point(ex(2),sigma,KK(i,:,k),fss,j)

     call std_covdiff(crho(i),sigma(j),fs,fr,fss,frr,frs,covf)

     Theta22(i,j,k) = covf*Theta22(i,j,k)
     Theta23(i,j,k) = covf*Theta23(i,j,k)
     Theta33(i,j,k) = covf*Theta33(i,j,k)
   enddo
   enddo
   else
  do i=1,ex(1)
  do j=1,ex(2)
  do k=1,ex(3)
    
     call rderivs_x_point(ex(1),crho,KK(:,j,k),fr,i)
     call rderivs_x_point(ex(2),sigma,KK(i,:,k),fs,j)
     call rdderivs_xy_point(ex(1),ex(2),crho,sigma,KK(:,:,k),frs,i,j)
     call rdderivs_x_point(ex(1),crho,KK(:,j,k),frr,i)
     call rdderivs_x_point(ex(2),sigma,KK(i,:,k),fss,j)

     call std_covdiff(crho(i),sigma(j),fs,fr,fss,frr,frs,covf)

     Theta22(i,j,k) = covf*Theta22(i,j,k)
     Theta23(i,j,k) = covf*Theta23(i,j,k)
     Theta33(i,j,k) = covf*Theta33(i,j,k)
   enddo
   enddo
   enddo
   endif

   return

end subroutine get_gauge_g00
#else
subroutine get_gauge_g00_K(ex,crho,sigma,X,g22,g23,g33, &
                         Theta22,Theta23,Theta33, g00, Rmin)
implicit none
! argument variables
integer, intent(in ):: ex(1:3)
real*8,intent(in) :: Rmin
double precision,intent(in),dimension(ex(1))::crho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::X
real*8,dimension(ex(1),ex(2),ex(3)),intent(in)::g22,g23,g33
real*8,dimension(ex(1),ex(2),ex(3)),intent(in)::Theta22,Theta23,Theta33
real*8,dimension(ex(1),ex(2),ex(3)),intent(out)::g00


double precision,dimension(ex(3))::R
real*8,dimension(ex(1),ex(2),ex(3))::det,gup22,gup23,gup33,KK

real*8 :: sr,ss,cr,cs,sr2,ss2,cr2,cs2,tg22,tg23,tg33
real*8 :: fr,fs,frr,fss,frs,covf

integer :: i,j,k

real*8 :: m0,Pp0,Pm0,ap,am,bp,bm,cp,cm,gam

call get_RT_parameters(m0,Pp0,Pm0,ap,am,bp,bm,cp,cm,gam)

R = X*Rmin/(1-X)
! g22 is P
det = dlog(g22**2)

  do i=1,ex(1)
  do j=1,ex(2)
  do k=1,ex(3)
    
     call rderivs_x_point(ex(1),crho,det(:,j,k),fr,i)
     call rderivs_x_point(ex(2),sigma,det(i,:,k),fs,j)
     call rdderivs_xy_point(ex(1),ex(2),crho,sigma,det(:,:,k),frs,i,j)
     call rdderivs_x_point(ex(1),crho,det(:,j,k),frr,i)
     call rdderivs_x_point(ex(2),sigma,det(i,:,k),fss,j)

     call std_covdiff(crho(i),sigma(j),fs,fr,fss,frr,frs,covf)

     KK(i,j,k) = covf
   enddo
   enddo
   enddo

   g00 = g22**2*(1+0.5*KK)

   return

end subroutine get_gauge_g00_K
! the input g00 is K
subroutine get_gauge_g00(ex,crho,sigma,X,g22,g23,g33, &
                         Theta22,Theta23,Theta33, g00, Rmin,fp)
implicit none
! argument variables
integer, intent(in ):: ex(1:3),fp
real*8,intent(in) :: Rmin
double precision,intent(in),dimension(ex(1))::crho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::X
real*8,dimension(ex(1),ex(2),ex(3)),intent(in)::g22,g23,g33
real*8,dimension(ex(1),ex(2),ex(3)),intent(out)::Theta22,Theta23,Theta33
real*8,dimension(ex(1),ex(2),ex(3)),intent(in)::g00


double precision,dimension(ex(3))::R
real*8,dimension(ex(1),ex(2),ex(3))::det,gup22,gup23,gup33,KK

real*8 :: sr,ss,cr,cs,sr2,ss2,cr2,cs2,tg22,tg23,tg33
real*8 :: fr,fs,frr,fss,frs,covf

integer :: i,j,k

real*8 :: m0,Pp0,Pm0,ap,am,bp,bm,cp,cm,gam

call get_RT_parameters(m0,Pp0,Pm0,ap,am,bp,bm,cp,cm,gam)

R = X*Rmin/(1-X)

  KK = g00

  if(fp == 0)then
  k = 1
  do i=1,ex(1)
  do j=1,ex(2)
    
     call rderivs_x_point(ex(1),crho,KK(:,j,k),fr,i)
     call rderivs_x_point(ex(2),sigma,KK(i,:,k),fs,j)
     call rdderivs_xy_point(ex(1),ex(2),crho,sigma,KK(:,:,k),frs,i,j)
     call rdderivs_x_point(ex(1),crho,KK(:,j,k),frr,i)
     call rdderivs_x_point(ex(2),sigma,KK(i,:,k),fss,j)

     call std_covdiff(crho(i),sigma(j),fs,fr,fss,frr,frs,covf)

     Theta22(i,j,k) = covf
   enddo
   enddo
   else
  do i=1,ex(1)
  do j=1,ex(2)
  do k=1,ex(3)
    
     call rderivs_x_point(ex(1),crho,KK(:,j,k),fr,i)
     call rderivs_x_point(ex(2),sigma,KK(i,:,k),fs,j)
     call rdderivs_xy_point(ex(1),ex(2),crho,sigma,KK(:,:,k),frs,i,j)
     call rdderivs_x_point(ex(1),crho,KK(:,j,k),frr,i)
     call rdderivs_x_point(ex(2),sigma,KK(i,:,k),fss,j)

     call std_covdiff(crho(i),sigma(j),fs,fr,fss,frr,frs,covf)

     Theta22(i,j,k) = covf
   enddo
   enddo
   enddo
   endif

   Theta22 = -Theta22/12/m0*g22**3
   return

end subroutine get_gauge_g00
#endif
!!---------------------------------------------------------------------------
subroutine get_gauge_g00_real(ex,crho,sigma,X,g22,g23,g33, &
                         Theta22,Theta23,Theta33, g00, Rmin)
implicit none
! argument variables
integer, intent(in ):: ex(1:3)
real*8,intent(in) :: Rmin
double precision,intent(in),dimension(ex(1))::crho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::X
real*8,dimension(ex(1),ex(2),ex(3)),intent(in)::g22,g23,g33
real*8,dimension(ex(1),ex(2),ex(3)),intent(in)::Theta22,Theta23,Theta33
real*8,dimension(ex(1),ex(2),ex(3)),intent(out)::g00


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

  do i=1,ex(1)
  do j=1,ex(2)
  do k=1,ex(3)
    
     call rderivs_x_point(ex(1),crho,det(:,j,k),fr,i)
     call rderivs_x_point(ex(2),sigma,det(i,:,k),fs,j)
     call rdderivs_xy_point(ex(1),ex(2),crho,sigma,det(:,:,k),frs,i,j)
     call rdderivs_x_point(ex(1),crho,det(:,j,k),frr,i)
     call rdderivs_x_point(ex(2),sigma,det(i,:,k),fss,j)

     call std_covdiff(crho(i),sigma(j),fs,fr,fss,frr,frs,covf)

     KK(i,j,k) = dsqrt(det(i,j,k))*(1-0.25*covf/R(k)**2)

     g00(i,j,k) = 2*m0/R(k)**3-KK(i,j,k)/R(k)**2 &
                 -(gup22(i,j,k)*Theta22(i,j,k)+2*gup23(i,j,k)*Theta23(i,j,k)+gup33(i,j,k)*Theta33(i,j,k))/2/R(k)
   enddo
   enddo
   enddo

   return

end subroutine get_gauge_g00_real
!!------------------------------------------------------------------------------------------------------------
subroutine get_null_boundary2(ex,crho,sigma,X,g22,g23,g33, &
                         g01,p02,p03,g02,g03,Theta22,Theta23,Theta33, Rmin)
implicit none
! argument variables
integer, intent(in ):: ex(1:3)
real*8,intent(in) :: Rmin
double precision,intent(in),dimension(ex(1))::crho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::X
real*8,dimension(ex(1),ex(2),ex(3)),intent(inout)::g22,g23,g33
real*8,dimension(ex(1),ex(2),ex(3)),intent(inout)::g01,p02,p03,g02,g03
real*8,dimension(ex(1),ex(2),ex(3)),intent(inout)::Theta22,Theta23,Theta33

#if 1
real*8 :: fact

!fact = X(1)/X(2)*((1-X(2))/(1-X(1)))
!fact = fact**2
! since we used gAB/r^2 instead of gAB, so fact = 1
fact = 1.d0

g22(:,:,1) = g22(:,:,2)*fact
g23(:,:,1) = g23(:,:,2)*fact
g33(:,:,1) = g33(:,:,2)*fact

g01(:,:,1) = -1.d0

p02(:,:,1) = 0.d0
p03(:,:,1) = 0.d0
g02(:,:,1) = 0.d0
g03(:,:,1) = 0.d0

! have done in get_gauge_g00 
!Theta22(:,:,1) = Theta22(:,:,2)*fact
!Theta23(:,:,1) = Theta23(:,:,2)*fact
!Theta33(:,:,1) = Theta33(:,:,2)*fact
#else
g01 = -1
g02 = 0
g03 = 0
#endif
return

end subroutine get_null_boundary2
!!!--------------------------------------------------------------------------------------------------------------
subroutine get_initial_null3(ex,crho,sigma,XX,g22,g23,g33,sst,Rmin)
implicit none
! argument variables
integer, intent(in ):: ex(1:3),sst
real*8,intent(in) :: Rmin
double precision,intent(in),dimension(ex(1))::crho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::XX
real*8,dimension(ex(1),ex(2),ex(3)),intent(out)::g22,g23,g33

double precision,dimension(ex(3))::R
real*8 :: sr,ss,cr,cs
real*8 :: sr2,ss2,cr2,cs2
integer :: i,j,k
real*8 :: ggr,tgrho,tgsigma
real*8 ::x,y,z,gr,gt,gp

real*8 :: gxx,gxy,gyy,tc,ts,PI

double complex :: Zslm,II,Jr,ctp
double complex :: swtf,z220

double complex :: beta0,C1,C2,mx,my,mlx,mly
integer :: nu,m

call initial_null_paramter(beta0,C1,C2,nu,m)

  II = dcmplx(0.d0,1.d0)

R = XX*Rmin/(1-XX)

PI = dacos(-1.d0)

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

     gxx = 1-sr2*ss2
     gxx = 1/gxx/gxx

     gxy = -sr*cr*ss*cs*gxx
     gyy = cr2*gxx
     gxx = cs2*gxx
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
      write(*,*) "get_initial_null2: not recognized sst = ",sst
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

    z220 = Zslm(2,2,m,gt,gp)*swtf**2

    if(sst==1 .or. sst==3 .or. sst==4)then
      mx = 2*tc*ts*(ts-II*tc)/dcos(sigma(j))
      my = 2*tc*ts*(ts+II*tc)/dcos(crho(i))
    else
      mx = 2*tc*ts*(ts+II*tc)/dcos(sigma(j))
      my = 2*tc*ts*(ts-II*tc)/dcos(crho(i))
    endif
    mlx = gxx*mx+gxy*my
    mly = gxy*mx+gyy*my
    Jr = (2.4d1*beta0+3.d0*II*nu*C1-II*nu**3*C2)/3.6d1+C1/4.d0/R(k)-C2/1.2d1/R(k)**3
    Jr = dsqrt(dble((2-1)*2*(2+1)*(2+2)))*dreal(Jr)*z220

    ctp = Jr*mlx*mlx+dsqrt(1+abs(Jr)**2)*dconjg(mlx)*mlx
    g22(i,j,k) = dreal(ctp)
    ctp = Jr*mlx*mly+dsqrt(1+abs(Jr)**2)*dconjg(mlx)*mly
    g23(i,j,k) = dreal(ctp)
    ctp = Jr*mly*mly+dsqrt(1+abs(Jr)**2)*dconjg(mly)*mly
    g33(i,j,k) = dreal(ctp)

  enddo
  enddo
  enddo

return

end subroutine get_initial_null3
!!!--------------------------------------------------------------------------------------------------------------
subroutine get_g00_with_t(time,ex,crho,sigma,XX,g00,Rmin,sst)
implicit none
! argument variables
integer, intent(in ):: ex(1:3),sst
real*8,intent(in) :: time,Rmin
double precision,intent(in),dimension(ex(1))::crho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::XX
real*8,dimension(ex(1),ex(2),ex(3)),intent(out)::g00

double precision,dimension(ex(3))::R
real*8 :: sr,ss,cr,cs
real*8 :: sr2,ss2,cr2,cs2
integer :: i,j,k
real*8 :: ggr,tgrho,tgsigma
real*8 ::x,y,z,gr,gt,gp

real*8 :: tc,ts,PI

double complex :: Zslm,II,Jr,Ur,Wr
double complex :: swtf,z020,z120,z220

double complex :: beta0,C1,C2
integer :: nu,m

call initial_null_paramter(beta0,C1,C2,nu,m)

  II = dcmplx(0.d0,1.d0)

R = XX*Rmin/(1-XX)

PI = dacos(-1.d0)

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
      write(*,*) "get_g00_with_t: not recognized sst = ",sst
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

    z020 = Zslm(0,2,m,gt,gp)
    z120 = Zslm(1,2,m,gt,gp)*swtf
    z220 = Zslm(2,2,m,gt,gp)*swtf**2

    Jr = (2.4d1*beta0+3.d0*II*nu*C1-II*nu**3*C2)/3.6d1+C1/4.d0/R(k)-C2/1.2d1/R(k)**3
    Ur = (-24*II*nu*beta0+3*nu*nu*C1-nu**4*C2)/36+2*beta0/R(k)+C1/2/R(k)**2+ &
         II*nu*C2/3/R(k)**3+C2/4/R(k)**4
    Wr = (24*II*nu*beta0-2*nu*C1+nu**4*C2)/6+ &
         (3*II*nu*C1-6*beta0-II*nu**3*C2)/3/R(k) - &
         nu**2*C2/R(k)**2+II*nu*C2/R(k)**3+C2/2/R(k)**4

    Jr = Jr*exp(II*nu*time)
    Ur = Ur*exp(II*nu*time)
    Wr = Wr*exp(II*nu*time)

    g00(i,j,k) = 2*(2*(2+1)*dsqrt(dble((2-1)*2*(2+1)*(2+2)))*dreal(Ur)**2* &
                    dreal(Jr)*dreal(z120**2*dconjg(z220))+ &
                    2*(2+1)*dsqrt(1+(2-1)*2*(2+1)*(2+2)*dreal(Jr)**2*abs(z220)**2)* &
                    dreal(Ur)**2*abs(z120)**2)-(1/R(k)**2+dreal(z020*Wr)/R(k))* &
                    exp(2*dreal(z020*beta0*exp(II*nu*time)))

  enddo
  enddo
  enddo

!if(sst==0 .and. crho(1) <-0.9 .and. sigma(1) <-0.9 .and. XX(1)<0.18182)write(*,*)"time = ",time,g00(1,1,1)

return

end subroutine get_g00_with_t
!!------------------------------------------------------------------------------------------------------------
subroutine get_null_boundary3(time,ex,crho,sigma,XX,g22,g23,g33, &
                         g01,p02,p03,g02,g03,Theta22,Theta23,Theta33, Rmin,sst)
implicit none
! argument variables
integer, intent(in ):: ex(1:3),sst
real*8,intent(in) ::time,Rmin
double precision,intent(in),dimension(ex(1))::crho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::XX
real*8,dimension(ex(1),ex(2),ex(3)),intent(inout)::g22,g23,g33
real*8,dimension(ex(1),ex(2),ex(3)),intent(inout)::g01,p02,p03,g02,g03
real*8,dimension(ex(1),ex(2),ex(3)),intent(inout)::Theta22,Theta23,Theta33

double precision,dimension(ex(3))::R
real*8 :: sr,ss,cr,cs
real*8 :: sr2,ss2,cr2,cs2
integer :: i,j,k
real*8 :: ggr,tgrho,tgsigma
real*8 ::x,y,z,gr,gt,gp

real*8 :: gxx,gxy,gyy,tc,ts,PI

double complex :: Zslm,II,Jr,ctp,Jrp,Jrt,Ur,Urp,Wr
double complex :: swtf,z020,z120,z220

double complex :: beta0,C1,C2,mx,my,mlx,mly
integer :: nu,m

call initial_null_paramter(beta0,C1,C2,nu,m)

  II = dcmplx(0.d0,1.d0)

R = XX*Rmin/(1-XX)

PI = dacos(-1.d0)

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

     gxx = 1-sr2*ss2
     gxx = 1/gxx/gxx

     gxy = -sr*cr*ss*cs*gxx
     gyy = cr2*gxx
     gxx = cs2*gxx
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

    z020 = Zslm(0,2,m,gt,gp)
    z120 = Zslm(1,2,m,gt,gp)*swtf
    z220 = Zslm(2,2,m,gt,gp)*swtf**2

    if(sst==1 .or. sst==3 .or. sst==4)then
      mx = 2*tc*ts*(ts-II*tc)/dcos(sigma(j))
      my = 2*tc*ts*(ts+II*tc)/dcos(crho(i))
    else
      mx = 2*tc*ts*(ts+II*tc)/dcos(sigma(j))
      my = 2*tc*ts*(ts-II*tc)/dcos(crho(i))
    endif
    mlx = gxx*mx+gxy*my
    mly = gxy*mx+gyy*my

    Jr = (2.4d1*beta0+3.d0*II*nu*C1-II*nu**3*C2)/3.6d1+C1/4.d0/R(k)-C2/1.2d1/R(k)**3
! Jrp = d Jr/d X instead of d Jr/d r    
    Jrp = -C1/4.d0/Rmin/XX(k)**2+C2/1.2d1*3/R(k)**2/Rmin/XX(k)**2
    Ur = (-24*II*nu*beta0+3*nu*nu*C1-nu**4*C2)/36+2*beta0/R(k)+C1/2/R(k)**2+ &
         II*nu*C2/3/R(k)**3+C2/4/R(k)**4
    Urp = -2*beta0/Rmin/XX(k)**2-C1/R(k)/Rmin/XX(k)**2- &
         II*nu*C2/R(k)**2/Rmin/XX(k)**2-C2/R(k)**3/Rmin/XX(k)**2
    Wr = (24*II*nu*beta0-2*nu*C1+nu**4*C2)/6+ &
         (3*II*nu*C1-6*beta0-II*nu**3*C2)/3/R(k) - &
         nu**2*C2/R(k)**2+II*nu*C2/R(k)**3+C2/2/R(k)**4

    Jr = Jr*exp(II*nu*time)
    Jrp = Jrp*exp(II*nu*time)
    Jrt = II*nu*Jr*exp(II*nu*time)
    Ur = Ur*exp(II*nu*time)
    Urp = Urp*exp(II*nu*time)
    Wr = Wr*exp(II*nu*time)

    Jr = dsqrt(dble((2-1)*2*(2+1)*(2+2)))*dreal(Jr)*z220
    Jrt = dsqrt(dble((2-1)*2*(2+1)*(2+2)))*dreal(Jrt)*z220
    Jrp = dsqrt(dble((2-1)*2*(2+1)*(2+2)))*dreal(Jrp)*z220

    g01(i,j,k) = -dexp(2*dreal(z020*beta0*exp(II*nu*time)))
#if 1
    g02(i,j,k) = -dsqrt(dble(2*(2+1)))*dreal(Ur)*dreal(mlx*(z120*dconjg(z220)* &
                  dconjg(Jr)+dsqrt(1+abs(Jr)**2)*dconjg(z120)))
    g03(i,j,k) = -dsqrt(dble(2*(2+1)))*dreal(Ur)*dreal(mly*(z120*dconjg(z220)* &
                  dconjg(Jr)+dsqrt(1+abs(Jr)**2)*dconjg(z120)))
#elif 0
    mlx = mlx/swtf
    mly = mly/swtf
    g02(i,j,k) = dreal(mlx)
    g03(i,j,k) = dreal(mly)
    !if(sst==0 .and. crho(1) <-0.9 .and. sigma(1) <-0.9 .and. XX(1)<0.18182)write(*,*)g02(i,j,k),g03(i,j,k)
    !if(crho(1) <-0.9 .and. sigma(1) <-0.9 .and. XX(1)<0.18182)write(*,*)g02(i,j,k),g03(i,j,k)
#else   
    select case (sst)
    case (0,1)
      tc =-dcos(gp)/(dcos(gt)**2*dcos(gp)**2-dcos(gt)**2-dcos(gp)**2)
      ts = dsin(gp)/(1-dsin(gt)**2*dcos(gp)**2)
    case (2,3)
      tc = 0
      ts = dcos(gp)/(dcos(gt)**2*dcos(gp)**2-dcos(gt)**2-dcos(gp)**2)
    case (4,5)
      tc = 0
      ts =-dsin(gp)/(1-dsin(gt)**2*dcos(gp)**2)
    end select
    g02(i,j,k) = gxx*tc+gxy*ts
    g03(i,j,k) = gxy*tc+gyy*ts
    !if(sst==0 .and. crho(1) <-0.9 .and. sigma(1) <-0.9 .and. XX(1)<0.18182)write(*,*)g02(i,j,k),g03(i,j,k)
    if(crho(1) <-0.9 .and. sigma(1) <-0.9 .and. XX(1)<0.18182)write(*,*)g02(i,j,k),g03(i,j,k),sst
    stop
#endif
    p02(i,j,k) = -dsqrt(dble(2*(2+1)))*dreal(Urp)*dreal(mlx*(z120*dconjg(Jr)+ &
                  dsqrt(1+abs(Jr)**2)*dconjg(z120))) &
                 -dsqrt(dble(2*(2+1)))*dreal(Ur)*dreal(mlx*(z120*dconjg(Jrp)+ &
                  abs(Jrp)*abs(Jr)/dsqrt(1+abs(Jr)**2)*dconjg(z120)))
    p03(i,j,k) = -dsqrt(dble(2*(2+1)))*dreal(Urp)*dreal(mly*(z120*dconjg(Jr)+ &
                  dsqrt(1+abs(Jr)**2)*dconjg(z120))) &
                 -dsqrt(dble(2*(2+1)))*dreal(Ur)*dreal(mly*(z120*dconjg(Jrp)+ &
                  abs(Jrp)*abs(Jr)/dsqrt(1+abs(Jr)**2)*dconjg(z120)))

    ctp = dconjg(Jr)*mlx*mlx+dsqrt(1+abs(Jr)**2)*dconjg(mlx)*mlx
    g22(i,j,k) = dreal(ctp)
    ctp = dconjg(Jr)*mlx*mly+dsqrt(1+abs(Jr)**2)*dconjg(mlx)*mly
    g23(i,j,k) = dreal(ctp)
    ctp = dconjg(Jr)*mly*mly+dsqrt(1+abs(Jr)**2)*dconjg(mly)*mly
    g33(i,j,k) = dreal(ctp)

    ctp = dconjg(Jrt)*mlx*mlx+abs(Jr)*abs(Jrt)/dsqrt(1+abs(Jr)**2)*dconjg(mlx)*mlx
    Theta22(i,j,k) = dreal(ctp)
    ctp = dconjg(Jrt)*mlx*mly+abs(Jr)*abs(Jrt)/dsqrt(1+abs(Jr)**2)*dconjg(mlx)*mly
    Theta23(i,j,k) = dreal(ctp)
    ctp = dconjg(Jrt)*mly*mly+abs(Jr)*abs(Jrt)/dsqrt(1+abs(Jr)**2)*dconjg(mly)*mly
    Theta33(i,j,k) = dreal(ctp)

  enddo
  enddo
  enddo

return

end subroutine get_null_boundary3
