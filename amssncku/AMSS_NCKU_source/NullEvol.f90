

#include "macrodef.fh"

!#define OLD

! 0: rk4, 1: Adams-Moulton

#define RKorAM 0 

function beta_rhs(xx,CJx,Kx)    result(gont)   
  implicit none
  double complex,intent(in) :: CJx
  real*8,intent(in) :: xx,Kx

  real*8 :: gont

  gont = xx*(1.d0-xx)/8.d0*(dreal(CJx*dconjg(CJx))-Kx*Kx)

  return

end function beta_rhs

function Q_rhs(xx,CJ,CJx,DCJx,KK,Ck,Ckx,Cnux,KKx,CBx,Cnu,DCJ,CB,CQ)  result(gont)   
  implicit none
  double complex,intent(in) :: CJ,CJx,DCJx,Ck,Ckx,Cnux,CBx,Cnu,DCJ,CB,CQ
  real*8,intent(in) :: xx,KK,KKx

  double complex :: gont

  gont = -KK*(Ckx+Cnux)+Cnu*KKx+CJ*dconjg(Ckx)+2.d0*CBx &
         +dconjg(Cnu)*CJx+dconjg(CJ)*DCJx-dconjg(Ck)*CJx &
         +(dconjg(Cnu)*(CJx-CJ*CJ*dconjg(CJx)) &
         +DCJ*(dconjg(CJx)-dconjg(CJ*CJ)*CJx)/2.d0/KK/KK) &
         -2.d0*(2.d0*CB+CQ)/xx/(1.d0-xx)

  return

end function Q_rhs

function U_rhs(xx,Rmin,beta,KK,CQ,CJ)  result(gont)   
  implicit none
  double complex,intent(in) :: CQ,CJ
  real*8,intent(in) :: xx,Rmin,beta,KK

  double complex :: gont

#if 1
  gont = dexp(2.d0*beta)/Rmin/xx/xx*(KK*CQ-CJ*dconjg(CQ))
#else
  gont = CQ/Rmin/xx/xx
#endif

#if 0  
  if(cdabs(gont)>1)then
          write(*,*)beta,KK,CQ,CJ
          stop
  endif
#endif

  return

end function U_rhs

function W_rhs(xx,Rmin,beta,KK,DCB,CB,CJ,Cnu,Ck,W, &
               CQ,bDCk,bDCnu,bDCB,bDCU,bDCUx,DCJ)  result(gont)   
  implicit none
  double complex,intent(in) :: DCB,CB,CJ,Cnu,Ck,CQ,bDCk,bDCnu,bDCB,bDCU,bDCUx,DCJ
  real*8,intent(in) :: xx,Rmin,beta,KK,W

  real*8 :: Ric,gont

  Ric = dreal(2.d0*KK+bDCnu-bDCk+(DCJ*dconjg(DCJ)-Cnu*dconjg(Cnu))/4.d0/KK)

  gont = dreal(dexp(2.d0*beta)*(Ric/2.d0-KK*(bDCB+CB*dconjg(CB))+dconjg(CJ)*(bDCB+CB*CB) &
        +(Cnu-Ck)*dconjg(CB))-1.d0+2.d0*Rmin*xx/(1.d0-xx)*(bDCU-W)                     &
        +Rmin*xx*xx/2.d0*bDCUx-dexp(2.d0*beta)/4.d0*                                    &
        (KK*KK-CJ*dconjg(CJ))*(KK*dconjg(CQ)-dconjg(CJ)*CQ)*CQ)

  gont = gont/Rmin/xx/xx

  return

end function W_rhs

function Theta_rhs(xx,Rmin,beta,betax,KK,KKx,CU,CUx,DCUx,bDCU,bDCUx,DCU,CB,DCB,W,Wx,CJ,DCJ,CJx,CJxx, &
                   DCJx,bDCB,Cnu,Cnux,Ck,Theta)  result(gont)   
  implicit none
  double complex,intent(in) :: CU,CUx,DCUx,bDCU,bDCUx,DCU,CB,DCB,CJ,DCJ,CJx,CJxx,DCJx
  double complex,intent(in) :: bDCB,Cnu,Cnux,Ck,Theta
  real*8,intent(in) :: xx,Rmin,beta,betax,KK,KKx,W,Wx

  double complex :: JH,II,gont
  real*8 :: V,Vx,Pu

  II = dcmplx(0.d0,1.d0)

  V = xx*Rmin/(1.d0-xx)*(1.d0+xx*Rmin/(1.d0-xx)*W)

  Vx = Rmin/(1.d0-xx)**2+2.d0*xx*Rmin*Rmin/(1.d0-xx)**3*W+xx*xx*Rmin*Rmin/(1.d0-xx)**2*Wx

  Pu = 2.d0*xx*(1.d0-xx)/KK*dreal(Theta*(dconjg(CJx)*KK-dconjg(CJ)*KKx))

  JH = (1.d0-xx)/xx/Rmin*dexp(2.d0*beta)*(-KK*DCJ*dconjg(CB)+ &
       (KK*Cnu+(KK*KK-1.d0)*DCJ-2.d0*KK*Ck)*CB+CJ* &
       ((2.d0*Ck-Cnu)*dconjg(CB)-2.d0*KK*(bDCB+CB*dconjg(CB))+ &
       2.d0*dreal((Cnu-Ck)*dconjg(CB)+dconjg(CJ)*(DCB+CB*CB)))) &
       +0.5d0*Rmin*xx**3*(1.d0-xx)*dexp(-2.d0*beta)* &
       ((KK*CUx+CJ*dconjg(CUx))**2- &
       CJ*dreal(dconjg(CUx)*(KK*CUx+CJ*dconjg(CUx)))) &
       -0.5d0*(Cnu*(xx*(1.d0-xx)*CUx+2.d0*CU)+DCJ*(xx*(1.d0-xx)*dconjg(CUx)+ &
       2.d0*dconjg(CU)))+CJ*II*dimag(xx*(1.d0-xx)*bDCUx+2.d0*bDCU) &
       -xx*(1.d0-xx)*CJx*dreal(bDCU) &
       +xx*(1.d0-xx)*(dconjg(CU)*DCJ+CU*Cnu)*II*dimag(CJ*dconjg(CJx)) &
       -xx*(1.d0-xx)*(dconjg(CU)*DCJx+CU*Cnux) &
       -2.d0*xx*(1.d0-xx)*(CJ*KKx-KK*CJx)*(dreal(dconjg(CU)*Ck)+ &
       II*dimag(KK*bDCU-dconjg(CJ)*DCU)) &
       -8.d0*CJ*((1.d0-xx)**2/Rmin+xx*(1.d0-xx)*W)*betax

  gont = -KK*(xx*(1-xx)*DCUx+2.d0*DCU)+2.d0*(1.d0-xx)/xx/Rmin*dexp(2.d0*beta)*(DCB+CB*CB) &
         -(xx*(1.d0-xx)*Wx+W)*CJ+JH+CJ*Pu-2.d0*Theta &
         -(1.d0-xx)*(1.d0-xx)/xx/xx/Rmin/Rmin*V*(CJ+xx*(1.d0-xx)*CJx) &
         +(1.d0-xx)*(1.d0-xx)*(1.d0-xx)/xx/Rmin/Rmin*Vx*(CJ+xx*(1.d0-xx)*CJx) &
         +(1.d0-xx)**4/xx/Rmin/Rmin*V*(2.d0*CJx+xx*CJxx)
#if 0 
  gont = -(xx*(1-xx)*DCUx+2.d0*DCU)+2.d0*(1.d0-xx)/xx/Rmin*DCB &
         -2.d0*Theta &
         +(1.d0-xx)**3/Rmin*(2.d0*CJx+xx*CJxx)
#endif

  gont = gont/2.d0/xx/(1.d0-xx)

  return

end function Theta_rhs
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////
subroutine fake_Theta_rhs(lx,X,rhs,Theta)
  implicit none
  integer,intent(in) :: lx
  double complex,dimension(lx),intent(in) :: Theta
  double complex,dimension(lx),intent(out) :: rhs
  real*8,dimension(lx),intent(in) :: X

  call cderivs_x(lx,X,Theta,rhs)

  return

end subroutine fake_Theta_rhs
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////
! try other guy's old method
function Theta_rhs_o(xx,Rmin,beta,betax,KK,KKx,CU,CUx,DCUx,bDCU,bDCUx,DCU,CB,DCB,W,Wx,CJ,DCJ,CJx,CJxx, &
                   DCJx,bDCB,Cnu,Cnux,Ck,Theta)  result(gont)   
  implicit none
  double complex,intent(in) :: CU,CUx,DCUx,bDCU,bDCUx,DCU,CB,DCB,CJ,DCJ,CJx,CJxx,DCJx
  double complex,intent(in) :: Cnu,Cnux,Ck,bDCB,Theta
  real*8,intent(in) :: xx,Rmin,beta,betax,KK,KKx,W,Wx

  double complex :: JH,II,gont
  real*8 :: V,Vx,Pu

  II = dcmplx(0.d0,1.d0)

  V = xx*Rmin/(1.d0-xx)*(1.d0+xx*Rmin/(1.d0-xx)*W)

  Vx = Rmin/(1.d0-xx)**2+2.d0*xx*Rmin*Rmin/(1.d0-xx)**3*W+xx*xx*Rmin*Rmin/(1.d0-xx)**2*Wx

  Pu = 2.d0*xx*(1.d0-xx)/KK*dreal(Theta*(dconjg(CJx)*KK-dconjg(CJ)*KKx))

  JH = (1.d0-xx)/xx/Rmin*dexp(2.d0*beta)*(-KK*DCJ*dconjg(CB)+ &
       (KK*Cnu+(KK*KK-1.d0)*DCJ-2.d0*KK*Ck)*CB+CJ* &
       ((2.d0*Ck-Cnu)*dconjg(CB)-2.d0*KK*(bDCB+CB*dconjg(CB))+ &
       2.d0*dreal((Cnu-Ck)*dconjg(CB)+dconjg(CJ)*(DCB+CB*CB)))) &
       +0.5d0*Rmin*xx**3*(1.d0-xx)*dexp(-2.d0*beta)* &
       ((KK*CUx+CJ*dconjg(CUx))**2- &
       CJ*dreal(dconjg(CUx)*(KK*CUx+CJ*dconjg(CUx)))) &
       -0.5d0*(Cnu*(xx*(1.d0-xx)*CUx+2.d0*CU)+DCJ*(xx*(1.d0-xx)*dconjg(CUx)+ &
       2.d0*dconjg(CU)))+CJ*II*dimag(xx*(1.d0-xx)*bDCUx+2.d0*bDCU) &
       -xx*(1.d0-xx)*CJx*dreal(bDCU) &
       +xx*(1.d0-xx)*(dconjg(CU)*DCJ+CU*Cnu)*II*dimag(CJ*dconjg(CJx)) &
       -xx*(1.d0-xx)*(dconjg(CU)*DCJx+CU*Cnux) &
       -2.d0*xx*(1.d0-xx)*(CJ*KKx-KK*CJx)*(dreal(dconjg(CU)*Ck)+ &
       II*dimag(KK*bDCU-dconjg(CJ)*DCU)) &
       -8.d0*CJ*((1.d0-xx)**2/Rmin+xx*(1.d0-xx)*W)*betax

  gont = -KK*(xx*(1-xx)*DCUx+2.d0*DCU)+2.d0*(1.d0-xx)/xx/Rmin*dexp(2.d0*beta)*(DCB+CB*CB) &
         -(xx*(1.d0-xx)*Wx+W)*CJ+JH+CJ*Pu &
         -(1.d0-xx)*(1.d0-xx)/xx/xx/Rmin/Rmin*V*(CJ+xx*(1.d0-xx)*CJx) &
         +(1.d0-xx)*(1.d0-xx)*(1.d0-xx)/xx/Rmin/Rmin*Vx*(CJ+xx*(1.d0-xx)*CJx) &
         +(1.d0-xx)**4/xx/Rmin/Rmin*V*(2.d0*CJx+xx*CJxx)

  return

end function Theta_rhs_o

#if (RKorAM == 0)

!--------------------------------------------------------------------
! this R is indeed x
function NullEvol_Theta_o(ex,crho,sigma,R,RJ,IJ,RU,IU,beta,RB,IB, &
                        Rnu,Inu,Rk,Ik,RTheta,ITheta,W,Rmin,       &
                        qlR1,qlR2,qlI1,qlI2,quR1,quR2,quI1,quI2,gR,gI) result(gont)   
  implicit none
  integer,intent(in ):: ex(1:3)
  real*8,intent(in),dimension(ex(1))::crho
  real*8,intent(in),dimension(ex(2))::sigma
  real*8,intent(in),dimension(ex(3))::R
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: beta,W
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: RJ,IJ,RU,IU,RB,IB,Rnu,Inu,Rk,Ik
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: RTheta,ITheta
  real*8,intent(in) :: Rmin
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: qlR1,qlR2,qlI1,qlI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: quR1,quR2,quI1,quI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: gR,gI
!  gont = 0: success; gont = 1: something wrong
  integer::gont

  double complex,dimension(ex(1),ex(2),ex(3)) :: CU,DCU,bDCU,CB,DCB,bDCB,CJ,DCJ
  double complex :: CTheta0,CTheta,CTheta1,RHS
  integer :: i,j,k,RK4
  double complex,dimension(ex(3)) :: Cnu,Ck,HCnu,HCk,HCU,HDCU,CUx,HCUx,DCUx,HDCUx,HbDCU,bDCUx,HbDCUx
  double complex,dimension(ex(3)) :: Cnux,HCnux,HCJ,HDCJ,CJx,HCJx,CJxx,HCJxx,DCJx,HDCJx,HCB,HDCB,HbDCB
  real*8,dimension(ex(3)) :: KK,KKx,HKK,HKKx,Hbeta,betax,Hbetax,HW,Wx,HWx
  double complex :: Theta_rhs_o
  real*8 :: dR

!!! sanity check
  dR = sum(RJ)+sum(IJ)+sum(RU)+sum(IU)+sum(beta)+sum(RB)+sum(IB) + &
       sum(Rnu)+sum(Inu)+sum(Rk)+sum(Ik)+sum(RTheta)+sum(ITheta)
  if(dR.ne.dR) then
     if(sum(RJ).ne.sum(RJ))write(*,*)"NullEvol_Theta: find NaN in RJ"
     if(sum(IJ).ne.sum(IJ))write(*,*)"NullEvol_Theta: find NaN in IJ"
     if(sum(RU).ne.sum(RU))write(*,*)"NullEvol_Theta: find NaN in RU"
     if(sum(IU).ne.sum(IU))write(*,*)"NullEvol_Theta: find NaN in IU"
     if(sum(beta).ne.sum(beta))write(*,*)"NullEvol_Theta: find NaN in beta"
     if(sum(RB).ne.sum(RB))write(*,*)"NullEvol_Theta: find NaN in RB"
     if(sum(IB).ne.sum(IB))write(*,*)"NullEvol_Theta: find NaN in IB"
     if(sum(Rnu).ne.sum(Rnu))write(*,*)"NullEvol_Theta: find NaN in Rnu"
     if(sum(Inu).ne.sum(Inu))write(*,*)"NullEvol_Theta: find NaN in Inu"
     if(sum(Rk).ne.sum(Rk))write(*,*)"NullEvol_Theta: find NaN in Rk"
     if(sum(Ik).ne.sum(Ik))write(*,*)"NullEvol_Theta: find NaN in Ik"
     if(sum(RTheta).ne.sum(RTheta))write(*,*)"NullEvol_Theta: find NaN in RTheta"
     if(sum(ITheta).ne.sum(ITheta))write(*,*)"NullEvol_Theta: find NaN in ITheta"
     gont = 1
     return
  endif

  dR = R(2) - R(1)
  
  CU = dcmplx(RU,IU)
  CB = dcmplx(RB,IB)
  CJ = dcmplx(RJ,IJ)

  do k=1,ex(3)
     call derivs_eth(ex(1:2),crho,sigma,CU(:,:,k),DCU(:,:,k),1,1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CU(:,:,k),bDCU(:,:,k),1,-1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CB(:,:,k),DCB(:,:,k),1,1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CB(:,:,k),bDCB(:,:,k),1,-1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CJ(:,:,k),DCJ(:,:,k),2,1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
  enddo

  do j=1,ex(2)
  do i=1,ex(1)
     CTheta0 = dcmplx(RTheta(i,j,1),ITheta(i,j,1))
     Cnu = dcmplx(Rnu(i,j,:),Inu(i,j,:))
     Ck = dcmplx(Rk(i,j,:),Ik(i,j,:))
     call cget_half_x(ex(3),CB(i,j,:),HCB)
     call cget_half_x(ex(3),DCB(i,j,:),HDCB)
     call cget_half_x(ex(3),bDCB(i,j,:),HbDCB)
     call cget_half_x(ex(3),Cnu,HCnu)
     call cderivs_x(ex(3),R,Cnu,Cnux)
     call cget_half_x(ex(3),Cnux,HCnux)
     call cget_half_x(ex(3),Ck,HCk)
     call rget_half_x(ex(3),beta(i,j,:),Hbeta)
     call rderivs_x(ex(3),R,beta(i,j,:),betax)
     call rget_half_x(ex(3),betax,Hbetax)
     KK = dsqrt(1.d0+RJ(i,j,:)*RJ(i,j,:)+IJ(i,j,:)*IJ(i,j,:))
     call rget_half_x(ex(3),KK,HKK)
     call rderivs_x(ex(3),R,KK,KKx)
     call rget_half_x(ex(3),KKx,HKKx)
     call rderivs_x(ex(3),R,W,Wx)
     call rget_half_x(ex(3),Wx,HWx)
     call rget_half_x(ex(3),W(i,j,:),HW)
     call cget_half_x(ex(3),CU(i,j,:),HCU)
     call cderivs_x(ex(3),R,DCU(i,j,:),DCUx)
     call cderivs_x(ex(3),R,CU(i,j,:),CUx)
     call cget_half_x(ex(3),DCUx,HDCUx)
     call cget_half_x(ex(3),CUx,HCUx)
     call cget_half_x(ex(3),DCU(i,j,:),HDCU)
     call cderivs_x(ex(3),R,bDCU(i,j,:),bDCUx)
     call cget_half_x(ex(3),bDCUx,HbDCUx)
     call cget_half_x(ex(3),bDCU(i,j,:),HbDCU)
     call cderivs_x(ex(3),R,CJ(i,j,:),CJx)
     call cdderivs_x(ex(3),R,CJ(i,j,:),CJxx)
     call cget_half_x(ex(3),CJx,HCJx)
     call cget_half_x(ex(3),CJxx,HCJxx)
     call cderivs_x(ex(3),R,DCJ(i,j,:),DCJx)
     call cget_half_x(ex(3),DCJx,HDCJx)
     do k=1,ex(3)-1
        RHS = Theta_rhs_o(R(k)+dR/2.d0,Rmin,Hbeta(k),betax(k),HKK(k),KKx(k),HCU(k),CUx(k),DCUx(k),HbDCU(k),bDCUx(k), &
                        HDCU(k),HCB(k),HDCB(k),HW(k),Wx(k),HCJ(k),HDCJ(k),    &
                        CJx(k),CJxx(k),DCJx(k),HbDCB(k),HCnu(k),Cnux(k),HCk(k),CTheta0)
        CTheta1 = RHS-(1-2.d0*(R(k)+dR/2.d0)*(1.d0-R(k)-dR/2.d0)/dR)*CTheta0
        CTheta1 = CTheta1/(1+2.d0*(R(k)+dR/2.d0)*(1.d0-R(k)-dR/2.d0)/dR)
        CTheta0 = CTheta1

        RTheta(i,j,k+1) = dreal(CTheta0)
        ITheta(i,j,k+1) = dimag(CTheta0)
     enddo
  enddo
  enddo

  gont = 0
  return

end function NullEvol_Theta_o
!------------------------------------------------------------------------------
! this R is indeed x
function NullEvol_beta(ex,crho,sigma,R,RJ,IJ,beta,KKx,HKKx)    result(gont)   
  implicit none
  integer,intent(in ):: ex(1:3)
  real*8,intent(in),dimension(ex(1))::crho
  real*8,intent(in),dimension(ex(2))::sigma
  real*8,intent(in),dimension(ex(3))::R
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: beta
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: RJ,IJ,KKx,HKKx
!  gont = 0: success; gont = 1: something wrong
  integer::gont
  real*8 :: dR

  double complex, dimension(ex(3)):: CJ,CJx,HCJx
  real*8 :: betah0,betah1,betah,rhs
  integer :: i,j,k,RK4
  real*8 :: beta_rhs

!!! sanity check
  dR = sum(RJ)+sum(IJ)+sum(beta)+sum(KKx)+sum(HKKx)
  if(dR.ne.dR) then
     if(sum(RJ).ne.sum(RJ))write(*,*)"NullEvol_beta: find NaN in RJ"
     if(sum(IJ).ne.sum(IJ))write(*,*)"NullEvol_beta: find NaN in IJ"
     if(sum(beta).ne.sum(beta))write(*,*)"NullEvol_beta: find NaN in beta"
     if(sum(KKx).ne.sum(KKx))write(*,*)"NullEvol_beta: find NaN in KKx"
     if(sum(HKKx).ne.sum(HKKx))write(*,*)"NullEvol_beta: find NaN in HKKx"
     gont = 1
     return
  endif

  dR = R(2) - R(1)

  do j=1,ex(2)
  do i=1,ex(1)
     betah0 = beta(i,j,1)
     CJ = dcmplx(RJ(i,j,:),IJ(i,j,:))
     call cderivs_x(ex(3),R,CJ,CJx)
     call cget_half_x(ex(3),CJx,HCJx)
#ifdef OLD     
     do k = 1,ex(3)-1
! note our CJx(ex(3)) = (CJ(ex(3))-CJ(ex(3)-1))/dR
! note our KKx(ex(3)) = (KK(ex(3))-KK(ex(3)-1))/dR
       rhs = beta_rhs(R(k)+dR/2.d0,CJx(k+1),KKx(i,j,k+1)) 
       beta(i,j,k+1) = beta(i,j,k) + rhs*dR
     enddo
#else
     do k=1,ex(3)-1
        RK4 = 0
        rhs = beta_rhs(R(k),CJx(k),KKx(i,j,k)) 
        call rungekutta4_scalar(dR,betah0,betah,rhs,RK4)

        RK4 = 1
        betah1 = beta_rhs(R(k)+dR/2.d0,HCJx(k),HKKx(i,j,k)) 
        call rungekutta4_scalar(dR,betah0,betah1,rhs,RK4)
        call rswap(betah,betah1)

        RK4 = 2
        betah1 = beta_rhs(R(k)+dR/2.d0,HCJx(k),HKKx(i,j,k)) 
        call rungekutta4_scalar(dR,betah0,betah1,rhs,RK4)
        call rswap(betah,betah1)

        RK4 = 3
        betah1 = beta_rhs(R(k+1),CJx(k+1),KKx(i,j,k+1)) 
        call rungekutta4_scalar(dR,betah0,betah1,rhs,RK4)
        call rswap(betah0,betah1)

        beta(i,j,k+1) = betah0
     enddo
! above k takes ex(3)-1 then do not need this closing step
#if 1
! closing step
     k = ex(3)-1
! note our CJx(ex(3)) = (CJ(ex(3))-CJ(ex(3)-1))/dR
! note our KKx(ex(3)) = (KK(ex(3))-KK(ex(3)-1))/dR
     rhs = beta_rhs(R(k)+dR/2.d0,CJx(k+1),KKx(i,j,k+1)) 
     beta(i,j,k+1) = beta(i,j,k) + rhs*dR
#endif

#endif
  enddo
  enddo

  gont = 0

  return

end function NullEvol_beta
!------------------------------------------------------------------------------
! this R is indeed x
function NullEvol_Q(ex,crho,sigma,R,RJ,IJ,Rk,Ik,Rnu,Inu,RB,IB,RQ,IQ,KK,Hkk,KKx,HKKx, &
                    qlR1,qlR2,qlI1,qlI2,quR1,quR2,quI1,quI2,gR,gI)    result(gont)   
  implicit none
  integer,intent(in ):: ex(1:3)
  real*8,intent(in),dimension(ex(1))::crho
  real*8,intent(in),dimension(ex(2))::sigma
  real*8,intent(in),dimension(ex(3))::R
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: RQ,IQ
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: RJ,IJ,KK,Hkk,KKx,HKKx
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: Rk,Ik,Rnu,Inu,RB,IB
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: qlR1,qlR2,qlI1,qlI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: quR1,quR2,quI1,quI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: gR,gI
!  gont = 0: success; gont = 1: something wrong
  integer::gont
  real*8 :: xx,dR

  double complex :: CQ0,CQ,CQ1,RHS
  double complex,dimension(ex(3)) :: CJx,HCJx,DCJx,HDCJx,Ck,Ckx,HCkx,Cnu,Cnux,HCnux,CB,CBx,HCBx
  double complex,dimension(ex(3)) :: HCJ,HCk,HCnu,HCB,HDCJ
  double complex, dimension(ex(1),ex(2),ex(3)) :: CJ,DCJ
  integer :: i,j,k,RK4
  double complex :: Q_rhs

!!! sanity check
  dR = sum(RJ)+sum(IJ)+sum(RQ)+sum(IQ) &
      +sum(RK)+sum(IK)+sum(Rnu)+sum(Inu)+sum(RB)+sum(IB) &
      +sum(KK)+sum(HKK)+sum(KKx)+sum(HKKx)
  if(dR.ne.dR) then
     if(sum(RJ).ne.sum(RJ))write(*,*)"NullEvol_Q: find NaN in RJ"
     if(sum(IJ).ne.sum(IJ))write(*,*)"NullEvol_Q: find NaN in IJ"
     if(sum(RQ).ne.sum(RQ))write(*,*)"NullEvol_Q: find NaN in RQ"
     if(sum(IQ).ne.sum(IQ))write(*,*)"NullEvol_Q: find NaN in IQ"
     if(sum(RK).ne.sum(RK))write(*,*)"NullEvol_Q: find NaN in RK"
     if(sum(IK).ne.sum(IK))write(*,*)"NullEvol_Q: find NaN in IK"
     if(sum(Rnu).ne.sum(Rnu))write(*,*)"NullEvol_Q: find NaN in Rnu"
     if(sum(Inu).ne.sum(Inu))write(*,*)"NullEvol_Q: find NaN in Inu"
     if(sum(RB).ne.sum(RB))write(*,*)"NullEvol_Q: find NaN in RB"
     if(sum(IB).ne.sum(IB))write(*,*)"NullEvol_Q: find NaN in IB"
     if(sum(KK).ne.sum(KK))write(*,*)"NullEvol_Q: find NaN in KK"
     if(sum(HKK).ne.sum(HKK))write(*,*)"NullEvol_Q: find NaN in HKK"
     if(sum(KKx).ne.sum(KKx))write(*,*)"NullEvol_Q: find NaN in KKx"
     if(sum(HKKx).ne.sum(HKKx))write(*,*)"NullEvol_Q: find NaN in HKKx"
     gont = 1
     return
  endif

  dR = R(2) - R(1)

  CJ = dcmplx(RJ,IJ)
  do k=1,ex(3)
     call derivs_eth(ex(1:2),crho,sigma,CJ(:,:,k),DCJ(:,:,k),2,1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
  enddo
  do j=1,ex(2)
  do i=1,ex(1)
     CQ0 = dcmplx(RQ(i,j,1),IQ(i,j,1))
     call cderivs_x(ex(3),R,CJ(i,j,:),CJx)
     call cget_half_x(ex(3),CJx,HCJx)
     call cget_half_x(ex(3),CJ,HCJ)
     call cderivs_x(ex(3),R,DCJ(i,j,:),DCJx)
     call cget_half_x(ex(3),DCJx,HDCJx)
     call cget_half_x(ex(3),DCJ,HDCJ)
     Ck = dcmplx(Rk(i,j,:),Ik(i,j,:))
     call cderivs_x(ex(3),R,Ck,Ckx)
     call cget_half_x(ex(3),Ckx,HCkx)
     call cget_half_x(ex(3),Ck,HCk)
     Cnu = dcmplx(Rnu(i,j,:),Inu(i,j,:))
     call cderivs_x(ex(3),R,Cnu,Cnux)
     call cget_half_x(ex(3),Cnux,HCnux)
     call cget_half_x(ex(3),Cnu,HCnu)
     CB = dcmplx(RB(i,j,:),IB(i,j,:))
     call cderivs_x(ex(3),R,CB,CBx)
     call cget_half_x(ex(3),CBx,HCBx)
     call cget_half_x(ex(3),CB,HCB)
#ifdef OLD    
     do k = 1,ex(3)-1
       xx = R(k)+dR/2.d0
! note our HCJ(ex(3)-1) = (CJ(ex(3))+CJ(ex(3)-1))/2
! note our KKx(ex(3)) = (KK(ex(3))-KK(ex(3)-1))/dR
       RHS = Q_rhs(xx,HCJ(k),CJx(k+1),DCJx(k+1),HKK(i,j,k),HCk(k),Ckx(k+1),Cnux(k+1),KKx(i,j,k+1),CBx(k+1),HCnu(k),HDCJ(k),HCB(k),0)
       RHS = RHS+CQ0*(1.d0/dR-1.d0/xx/(1.d0-xx))
       CQ0 = RHS/(1.d0/dR+1.d0/xx/(1.d0-xx))
       RQ(i,j,k+1) = dreal(CQ0)
       IQ(i,j,k+1) = dimag(CQ0)
       enddo
#else
     do k=1,ex(3)-2
        RK4 = 0
        RHS = Q_rhs(R(k),CJ(i,j,k),CJx(k),DCJx(k),KK(i,j,k),Ck(k),Ckx(k),Cnux(k),KKx(i,j,k),CBx(k),Cnu(k),DCJ(i,j,k),CB(k),CQ0)
        call rungekutta4_cplxscalar(dR,CQ0,CQ,RHS,RK4)

        RK4 = 1
        CQ1 = Q_rhs(R(k)+dR/2.d0,HCJ(k),HCJx(k),HDCJx(k),HKK(i,j,k),HCk(k),HCkx(k),HCnux(k),HKKx(i,j,k), &
                    HCBx(k),HCnu(k),HDCJ(k),HCB(k),CQ)
        call rungekutta4_cplxscalar(dR,CQ0,CQ1,RHS,RK4)
        call cswap(CQ,CQ1)

        RK4 = 2
        CQ1 = Q_rhs(R(k)+dR/2.d0,HCJ(k),HCJx(k),HDCJx(k),HKK(i,j,k),HCk(k),HCkx(k),HCnux(k),HKKx(i,j,k), &
                    HCBx(k),HCnu(k),HDCJ(k),HCB(k),CQ)
        call rungekutta4_cplxscalar(dR,CQ0,CQ1,RHS,RK4)
        call cswap(CQ,CQ1)
        
        RK4 = 3
        CQ1 = Q_rhs(R(k+1),CJ(i,j,k+1),CJx(k+1),DCJx(k+1),KK(i,j,k+1),Ck(k+1),Ckx(k+1),Cnux(k+1),KKx(i,j,k+1), &
                    CBx(k+1),Cnu(k+1),DCJ(i,j,k+1),CB(k+1),CQ)
        call rungekutta4_cplxscalar(dR,CQ0,CQ1,RHS,RK4)
        call cswap(CQ0,CQ1)

        RQ(i,j,k+1) = dreal(CQ0)
        IQ(i,j,k+1) = dimag(CQ0)
     enddo
#if 0
     k = ex(3)
     CQ0 = -2*CB(k)
     RQ(i,j,k+1) = dreal(CQ0)
     IQ(i,j,k+1) = dimag(CQ0)
#else
! closing step
     k = ex(3)-1
     CQ0 = dcmplx(RQ(i,j,k),IQ(i,j,k))
     xx = R(k)+dR/2.d0
! note our HCJ(ex(3)-1) = (CJ(ex(3))+CJ(ex(3)-1))/2
! note our KKx(ex(3)) = (KK(ex(3))-KK(ex(3)-1))/dR
     RHS = Q_rhs(xx,HCJ(k),CJx(k+1),DCJx(k+1),HKK(i,j,k), &
                 HCk(k),Ckx(k+1),Cnux(k+1),KKx(i,j,k+1),CBx(k+1),HCnu(k),HDCJ(k),HCB(k),dcmplx(0.d0,0.d0))
     RHS = RHS+CQ0*(1.d0/dR-1.d0/xx/(1.d0-xx))
     CQ0 = RHS/(1.d0/dR+1.d0/xx/(1.d0-xx))
     RQ(i,j,k+1) = dreal(CQ0)
     IQ(i,j,k+1) = dimag(CQ0)
#endif

#endif     
  enddo
  enddo

  gont = 0
  return

end function NullEvol_Q
!--------------------------------------------------------------------
! this R is indeed x
function NullEvol_U(ex,crho,sigma,R,RJ,IJ,RQ,IQ,KK,HKK,beta,RU,IU, &
                    Rmin) result(gont)   
  implicit none
  integer,intent(in ):: ex(1:3)
  real*8,intent(in),dimension(ex(1))::crho
  real*8,intent(in),dimension(ex(2))::sigma
  real*8,intent(in),dimension(ex(3))::R
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: RJ,IJ,RQ,IQ,beta,KK,HKK
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: RU,IU
  real*8,intent(in) :: Rmin
!  gont = 0: success; gont = 1: something wrong
  integer::gont
  real*8 :: dR

  double complex :: CU0,CU,CU1,RHS
  integer :: i,j,k,RK4
  double complex,dimension(ex(3)) :: CJ,CQ,HCJ,HCQ
  real*8,dimension(ex(3)) :: Hbeta
  double complex :: U_rhs

!!! sanity check
  dR = sum(RJ)+sum(IJ)+sum(RQ)+sum(IQ)+sum(beta)+sum(RU)+sum(IU)+sum(KK)+sum(HKK)
  if(dR.ne.dR) then
     if(sum(RJ).ne.sum(RJ))write(*,*)"NullEvol_U: find NaN in RJ"
     if(sum(IJ).ne.sum(IJ))write(*,*)"NullEvol_U: find NaN in IJ"
     if(sum(RQ).ne.sum(RQ))write(*,*)"NullEvol_U: find NaN in RQ"
     if(sum(IQ).ne.sum(IQ))write(*,*)"NullEvol_U: find NaN in IQ"
     if(sum(beta).ne.sum(beta))write(*,*)"NullEvol_U: find NaN in beta"
     if(sum(RU).ne.sum(RU))write(*,*)"NullEvol_U: find NaN in RU"
     if(sum(IU).ne.sum(IU))write(*,*)"NullEvol_U: find NaN in IU"
     if(sum(KK).ne.sum(KK))write(*,*)"NullEvol_U: find NaN in KK"
     if(sum(HKK).ne.sum(HKK))write(*,*)"NullEvol_U: find NaN in HKK"
     gont = 1
     return
  endif

  dR = R(2) - R(1)
  
  do j=1,ex(2)
  do i=1,ex(1)
     CU0 = dcmplx(RU(i,j,1),IU(i,j,1))
     CJ = dcmplx(RJ(i,j,:),IJ(i,j,:))
     CQ = dcmplx(RQ(i,j,:),IQ(i,j,:))
     call cget_half_x(ex(3),CJ,HCJ)
     call cget_half_x(ex(3),CQ,HCQ)
     call rget_half_x(ex(3),beta(i,j,:),Hbeta)
#ifdef OLD     
     do k = 1,ex(3)-1
! note our HCJ(ex(3)-1) = (CJ(ex(3))+CJ(ex(3)-1))/2
! note our KKx(ex(3)) = (KK(ex(3))-KK(ex(3)-1))/dR
       RHS = U_rhs(R(k)+dR/2,Rmin,Hbeta(k),HKK(i,j,k),HCQ(k),HCJ(k))
       CU0 = CU0+RHS*dR
       RU(i,j,k+1) = dreal(CU0)
       IU(i,j,k+1) = dimag(CU0)
     enddo
#else

     do k=1,ex(3)-2
 
        RK4 = 0
        RHS = U_rhs(R(k),Rmin,beta(i,j,k),KK(i,j,k),CQ(k),CJ(k))
        call rungekutta4_cplxscalar(dR,CU0,CU,RHS,RK4)

        RK4 = 1
        CU1 = U_rhs(R(k)+dR/2.d0,Rmin,Hbeta(k),HKK(i,j,k),HCQ(k),HCJ(k))
        call rungekutta4_cplxscalar(dR,CU0,CU1,RHS,RK4)
        call cswap(CU,CU1)

        RK4 = 2
        CU1 = U_rhs(R(k)+dR/2.d0,Rmin,Hbeta(k),HKK(i,j,k),HCQ(k),HCJ(k))
        call rungekutta4_cplxscalar(dR,CU0,CU1,RHS,RK4)
        call cswap(CU,CU1)
        
        RK4 = 3
        CU1 = U_rhs(R(k+1),Rmin,beta(i,j,k+1),KK(i,j,k+1),CQ(k+1),CJ(k+1))
        call rungekutta4_cplxscalar(dR,CU0,CU1,RHS,RK4)
        call cswap(CU0,CU1)

        RU(i,j,k+1) = dreal(CU0)
        IU(i,j,k+1) = dimag(CU0)

     enddo
! above k takes ex(3)-1 then do not need closing step     
#if 1
! closing step
     k = ex(3)-1
     CU0 = dcmplx(RU(i,j,k),IU(i,j,k))
! note our HCJ(ex(3)-1) = (CJ(ex(3))+CJ(ex(3)-1))/2
! note our KKx(ex(3)) = (KK(ex(3))-KK(ex(3)-1))/dR
     RHS = U_rhs(R(k)+dR/2,Rmin,Hbeta(k),HKK(i,j,k),HCQ(k),HCJ(k))
     CU0 = CU0+RHS*dR
     RU(i,j,k+1) = dreal(CU0)
     IU(i,j,k+1) = dimag(CU0)
#endif

#endif     
  enddo
  enddo

  gont = 0
  return

end function NullEvol_U
!----------------------------------------------------------------------------------------
! this R is indeed x
function NullEvol_W(ex,crho,sigma,R,RJ,IJ,RB,IB,Rnu,Inu,Rk,Ik, &
                    RU,IU,RQ,IQ,W,beta,KK,HKK,Rmin, &
                    qlR1,qlR2,qlI1,qlI2,quR1,quR2,quI1,quI2,gR,gI)    result(gont)   
  implicit none
  integer,intent(in ):: ex(1:3)
  real*8,intent(in),dimension(ex(1))::crho
  real*8,intent(in),dimension(ex(2))::sigma
  real*8,intent(in),dimension(ex(3))::R
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: W
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: RJ,IJ,RB,IB
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: Rnu,Inu,Rk,Ik
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: RU,IU,RQ,IQ,beta,KK,HKK
  real*8,intent(in ) :: Rmin
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: qlR1,qlR2,qlI1,qlI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: quR1,quR2,quI1,quI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: gR,gI
!  gont = 0: success; gont = 1: something wrong
  integer::gont
  real*8 :: dR

  real*8, dimension(ex(3)) :: Hbeta
  double complex, dimension(ex(1),ex(2),ex(3)) :: CU,DCU,bDCU
  double complex, dimension(ex(1),ex(2),ex(3)) :: CB,DCB,bDCB,CJ,DCJ,Cnu,bDCnu,Ck,bDCk
  double complex, dimension(ex(3)) :: HCB,HDCB,HbDCB,HCJ,HDCJ,HCnu,HbDCnu,HCk,HbDCk
  double complex, dimension(ex(3)) :: HbDCU,bDCUx,HbDCUx,CQ,HCQ
  real*8 :: Wh0,Wh1,Wh,rhs
  integer :: i,j,k,RK4
  real*8 :: xx,W_rhs

!!! sanity check
  dR = sum(RJ)+sum(IJ)+sum(beta)+sum(RB)+sum(IB)+sum(Rnu)+sum(Inu) &
      +sum(Rk)+sum(Ik)+sum(W)+sum(RU)+sum(IU)+sum(RQ)+sum(IQ)&
      +sum(KK)+sum(HKK)
  if(dR.ne.dR) then
     if(sum(RJ).ne.sum(RJ))write(*,*)"NullEvol_W: find NaN in RJ"
     if(sum(IJ).ne.sum(IJ))write(*,*)"NullEvol_W: find NaN in IJ"
     if(sum(beta).ne.sum(beta))write(*,*)"NullEvol_W: find NaN in beta"
     if(sum(RB).ne.sum(RB))write(*,*)"NullEvol_W: find NaN in RB"
     if(sum(IB).ne.sum(IB))write(*,*)"NullEvol_W: find NaN in IB"
     if(sum(Rnu).ne.sum(Rnu))write(*,*)"NullEvol_W: find NaN in Rnu"
     if(sum(Inu).ne.sum(Inu))write(*,*)"NullEvol_W: find NaN in Inu"
     if(sum(Rk).ne.sum(Rk))write(*,*)"NullEvol_W: find NaN in Rk"
     if(sum(Ik).ne.sum(Ik))write(*,*)"NullEvol_W: find NaN in Ik"
     if(sum(RU).ne.sum(RU))write(*,*)"NullEvol_W: find NaN in RU"
     if(sum(IU).ne.sum(IU))write(*,*)"NullEvol_W: find NaN in IU"
     if(sum(RQ).ne.sum(RQ))write(*,*)"NullEvol_W: find NaN in RQ"
     if(sum(IQ).ne.sum(IQ))write(*,*)"NullEvol_W: find NaN in IQ"
     if(sum(W).ne.sum(W))write(*,*)"NullEvol_W: find NaN in W"
     if(sum(KK).ne.sum(KK))write(*,*)"NullEvol_W: find NaN in KK"
     if(sum(HKK).ne.sum(HKK))write(*,*)"NullEvol_W: find NaN in HKK"
     gont = 1
     return
  endif

  dR = R(2) - R(1)
  
  CB = dcmplx(RB,IB)
  CU = dcmplx(RU,IU)
  Ck = dcmplx(Rk,Ik)
  Cnu = dcmplx(Rnu,Inu)
  CJ = dcmplx(RJ,IJ)
  do k=1,ex(3)
     call derivs_eth(ex(1:2),crho,sigma,CJ(:,:,k),DCJ(:,:,k),2,1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CB(:,:,k),DCB(:,:,k),1,1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CB(:,:,k),bDCB(:,:,k),1,-1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CU(:,:,k),bDCU(:,:,k),1,-1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,Ck(:,:,k),bDCk(:,:,k),1,-1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,Cnu(:,:,k),bDCnu(:,:,k),1,-1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
  enddo

  do j=1,ex(2)
  do i=1,ex(1)
     Wh0 = W(i,j,1)
     call rget_half_x(ex(3),beta(i,j,:),Hbeta)
     call cderivs_x(ex(3),R,bDCU,bDCUx)
     call cget_half_x(ex(3),bDCUx,HbDCUx)
     call cget_half_x(ex(3),bDCU,HbDCU)
     call cget_half_x(ex(3),DCJ,HDCJ)
     call cget_half_x(ex(3),DCB,HDCB)
     call cget_half_x(ex(3),bDCB,HbDCB)
     call cget_half_x(ex(3),CB,HCB)
     call cget_half_x(ex(3),CJ,HCJ)
     call cget_half_x(ex(3),Cnu,HCnu)
     call cget_half_x(ex(3),Ck,HCk)
     CQ = dcmplx(RQ(i,j,:),IQ(i,j,:))
     call cget_half_x(ex(3),CQ,HCQ)
     call cget_half_x(ex(3),bDCk,HbDCk)
     call cget_half_x(ex(3),bDCnu,HbDCnu)
#ifdef OLD   
     do k = 1,ex(3)-1
       xx = R(k)+dR/2
! note our HCJ(ex(3)-1) = (CJ(ex(3))+CJ(ex(3)-1))/2
! note our KKx(ex(3)) = (KK(ex(3))-KK(ex(3)-1))/dR
       rhs = W_rhs(xx,Rmin,Hbeta(k),HKK(i,j,k),HDCB(k),HCB(k),HCJ(k),HCnu(k),HCk(k),0, &
                 HCQ(k),HbDCk(k),HbDCnu(k),HbDCB(k),HbDCU(k),bDCUx(k+1),HDCJ(k))
       rhs = rhs+Wh0*(1.d0/dR-1.d0/xx/(1.d0-xx))                 
       W(i,j,k+1) = rhs/(1.d0/dR+1.d0/xx/(1.d0-xx)) 
     enddo
#else
     do k=1,ex(3)-2
        RK4 = 0
        rhs = W_rhs(R(k),Rmin,beta(i,j,k),KK(i,j,k),DCB(i,j,k),CB(i,j,k),CJ(i,j,k),Cnu(i,j,k),Ck(i,j,k),Wh0, &
                    CQ(k),bDCk(i,j,k),bDCnu(i,j,k),bDCB(i,j,k),bDCU(i,j,k),bDCUx(k),DCJ(i,j,k))
        call rungekutta4_scalar(dR,Wh0,Wh,rhs,RK4)

        RK4 = 1
        Wh1 = W_rhs(R(k)+dR/2.d0,Rmin,Hbeta(k),HKK(i,j,k),HDCB(k),HCB(k),HCJ(k),HCnu(k),HCk(k),Wh, &
                    HCQ(k),HbDCk(k),HbDCnu(k),HbDCB(k),HbDCU(k),HbDCUx(k),HDCJ(k))
        call rungekutta4_scalar(dR,Wh0,Wh1,rhs,RK4)
        call rswap(Wh,Wh1)

        RK4 = 2       
        Wh1 = W_rhs(R(k)+dR/2.d0,Rmin,Hbeta(k),HKK(i,j,k),HDCB(k),HCB(k),HCJ(k),HCnu(k),HCk(k),Wh, &
                    HCQ(k),HbDCk(k),HbDCnu(k),HbDCB(k),HbDCU(k),HbDCUx(k),HDCJ(k))
        call rungekutta4_scalar(dR,Wh0,Wh1,rhs,RK4)
        call rswap(Wh,Wh1)

        RK4 = 3
        Wh1 = W_rhs(R(k+1),Rmin,beta(i,j,k+1),KK(i,j,k+1),DCB(i,j,k+1),CB(i,j,k+1),CJ(i,j,k+1),Cnu(i,j,k+1),Ck(i,j,k+1),Wh, &
                    CQ(k+1),bDCk(i,j,k+1),bDCnu(i,j,k+1),bDCB(i,j,k+1),bDCU(i,j,k+1),bDCUx(k+1),DCJ(i,j,k+1))
        call rungekutta4_scalar(dR,Wh0,Wh1,rhs,RK4)
        call rswap(Wh0,Wh1)

        W(i,j,k+1) = Wh0
     enddo
#if 0
     k = ex(3)
     W(i,j,k) = dreal(bDCU(i,j,k))
#else
! closing step
     k = ex(3)-1
     Wh0 = W(i,j,k)
     xx = R(k)+dR/2
! note our HCJ(ex(3)-1) = (CJ(ex(3))+CJ(ex(3)-1))/2
! note our KKx(ex(3)) = (KK(ex(3))-KK(ex(3)-1))/dR
     rhs = W_rhs(xx,Rmin,Hbeta(k),HKK(i,j,k),HDCB(k),HCB(k),HCJ(k),HCnu(k),HCk(k),0.d0, &
                 HCQ(k),HbDCk(k),HbDCnu(k),HbDCB(k),HbDCU(k),bDCUx(k+1),HDCJ(k))
     rhs = rhs+Wh0*(1.d0/dR-1.d0/xx/(1.d0-xx))
     W(i,j,k+1) = rhs/(1.d0/dR+1.d0/xx/(1.d0-xx))
#endif

#endif     
  enddo
  enddo

  gont = 0
  return

end function NullEvol_W
!-----------------------------------------------------------------------------------------------
! given exact Theta_x
! this R is indeed x
function NullEvol_Theta_givenx(ex,crho,sigma,R,RJ,IJ,RU,IU,beta,RB,IB, &
                        Rnu,Inu,Rk,Ik,RTheta,ITheta,W,KK,HKK,KKx,HKKx,Rmin,       &
                        qlR1,qlR2,qlI1,qlI2,quR1,quR2,quI1,quI2,gR,gI, &
                        dquR1,dquR2,dquI1,dquI2,                          &
                        bdquR1,bdquR2,bdquI1,bdquI2,                      &
                        dgR,dgI,bdgR,bdgI,T,sst) result(gont)   
  implicit none
  integer,intent(in ):: ex(1:3),sst
  real*8,intent(in),dimension(ex(1))::crho
  real*8,intent(in),dimension(ex(2))::sigma
  real*8,intent(in),dimension(ex(3))::R
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: beta,W,KK,HKK,KKx,HKKx
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: RJ,IJ,RU,IU,RB,IB,Rnu,Inu,Rk,Ik
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: RTheta,ITheta
  real*8,intent(in) :: Rmin,T
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: qlR1,qlR2,qlI1,qlI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: quR1,quR2,quI1,quI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: gR,gI
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dquR1,dquR2,dquI1,dquI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: bdquR1,bdquR2,bdquI1,bdquI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dgR,dgI,bdgR,bdgI
!  gont = 0: success; gont = 1: something wrong
  integer::gont

  real*8,dimension(ex(3))::HR
  real*8,dimension(ex(1),ex(2),ex(3)) :: RThetax,IThetax,HRThetax,HIThetax
  double complex,dimension(ex(3)) :: fRHS,HfRHS
  real*8 :: xx,dR
  integer :: i,j,k,RK4
  double complex :: CTheta0,CTheta,CTheta1,RHS
  integer,parameter :: ks=1

!!! sanity check
  dR = sum(RJ)+sum(IJ)+sum(RU)+sum(IU)+sum(beta)+sum(RB)+sum(IB) + &
       sum(Rnu)+sum(Inu)+sum(Rk)+sum(Ik)+sum(RTheta)+sum(ITheta) + &
       sum(KK)+sum(HKK)+sum(KKx)+sum(HKKx)
  if(dR.ne.dR) then
     if(sum(RJ).ne.sum(RJ))write(*,*)"NullEvol_Theta: find NaN in RJ"
     if(sum(IJ).ne.sum(IJ))write(*,*)"NullEvol_Theta: find NaN in IJ"
     if(sum(RU).ne.sum(RU))write(*,*)"NullEvol_Theta: find NaN in RU"
     if(sum(IU).ne.sum(IU))write(*,*)"NullEvol_Theta: find NaN in IU"
     if(sum(beta).ne.sum(beta))write(*,*)"NullEvol_Theta: find NaN in beta"
     if(sum(RB).ne.sum(RB))write(*,*)"NullEvol_Theta: find NaN in RB"
     if(sum(IB).ne.sum(IB))write(*,*)"NullEvol_Theta: find NaN in IB"
     if(sum(Rnu).ne.sum(Rnu))write(*,*)"NullEvol_Theta: find NaN in Rnu"
     if(sum(Inu).ne.sum(Inu))write(*,*)"NullEvol_Theta: find NaN in Inu"
     if(sum(Rk).ne.sum(Rk))write(*,*)"NullEvol_Theta: find NaN in Rk"
     if(sum(Ik).ne.sum(Ik))write(*,*)"NullEvol_Theta: find NaN in Ik"
     if(sum(RTheta).ne.sum(RTheta))write(*,*)"NullEvol_Theta: find NaN in RTheta"
     if(sum(ITheta).ne.sum(ITheta))write(*,*)"NullEvol_Theta: find NaN in ITheta"
     if(sum(KK).ne.sum(KK))write(*,*)"NullEvol_Theta: find NaN in KK"
     if(sum(HKK).ne.sum(HKK))write(*,*)"NullEvol_Theta: find NaN in HKK"
     if(sum(KKx).ne.sum(KKx))write(*,*)"NullEvol_Theta: find NaN in KKx"
     if(sum(HKKx).ne.sum(HKKx))write(*,*)"NullEvol_Theta: find NaN in HKKx"
     gont = 1
     return
  endif

  dR = R(2) - R(1)
  HR = R+dR/2
  
  call get_exact_null_theta_x(ex,crho,sigma,R,RThetax,IThetax,sst,Rmin,T,               &
                        quR1,quR2,quI1,quI2,qlR1,qlR2,qlI1,qlI2,          &
                        gR,gI,                                            &
                        dquR1,dquR2,dquI1,dquI2,                          &
                        bdquR1,bdquR2,bdquI1,bdquI2,                      &
                        dgR,dgI,bdgR,bdgI)
  call get_exact_null_theta_x(ex,crho,sigma,HR,HRThetax,HIThetax,sst,Rmin,T,               &
                        quR1,quR2,quI1,quI2,qlR1,qlR2,qlI1,qlI2,          &
                        gR,gI,                                            &
                        dquR1,dquR2,dquI1,dquI2,                          &
                        bdquR1,bdquR2,bdquI1,bdquI2,                      &
                        dgR,dgI,bdgR,bdgI)
  do j=1,ex(2)
  do i=1,ex(1)
     CTheta0 = dcmplx(RTheta(i,j,ks),ITheta(i,j,ks))
     fRHS = dcmplx(RThetax(i,j,:),IThetax(i,j,:))
     HfRHS = dcmplx(HRThetax(i,j,:),HIThetax(i,j,:))
     ! call cget_half_x(ex(3),fRHS,HfRHS)

     do k=ks,ex(3)-1
        RK4 = 0
        RHS = fRHS(k)
        call rungekutta4_cplxscalar(dR,CTheta0,CTheta,RHS,RK4)

        RK4 = 1
        CTheta1 = HfRHS(k)
        call rungekutta4_cplxscalar(dR,CTheta0,CTheta1,RHS,RK4)
        call cswap(CTheta,CTheta1)

        RK4 = 2
        CTheta1 = HfRHS(k)
        call rungekutta4_cplxscalar(dR,CTheta0,CTheta1,RHS,RK4)
        call cswap(CTheta,CTheta1)
        
        RK4 = 3
        CTheta1 = fRHS(k+1)
        call rungekutta4_cplxscalar(dR,CTheta0,CTheta1,RHS,RK4)
        call cswap(CTheta0,CTheta1)

        RTheta(i,j,k+1) = dreal(CTheta0)
        ITheta(i,j,k+1) = dimag(CTheta0)
     enddo

#if 0 
! closing step
     k = ex(3)-1
     RHS = fRHS(k)               
     CTheta0 = dcmplx(RTheta(i,j,k),ITheta(i,j,k))+RHS*dR
     RTheta(i,j,k+1) = dreal(CTheta0)
     ITheta(i,j,k+1) = dimag(CTheta0)
#endif

  enddo
  enddo

  gont = 0
  return

end function NullEvol_Theta_givenx
!-----------------------------------------------------------------------------------------------
#if 1
! real evolve
! for eth_x, eth first, _x later
! this R is indeed x
function NullEvol_Theta(ex,crho,sigma,R,RJ,IJ,RU,IU,beta,RB,IB, &
                        Rnu,Inu,Rk,Ik,RTheta,ITheta,W,KK,HKK,KKx,HKKx,Rmin,       &
                        qlR1,qlR2,qlI1,qlI2,quR1,quR2,quI1,quI2,gR,gI) result(gont)   
  implicit none
  integer,intent(in ):: ex(1:3)
  real*8,intent(in),dimension(ex(1))::crho
  real*8,intent(in),dimension(ex(2))::sigma
  real*8,intent(in),dimension(ex(3))::R
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: beta,W,KK,HKK,KKx,HKKx
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: RJ,IJ,RU,IU,RB,IB,Rnu,Inu,Rk,Ik
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: RTheta,ITheta
  real*8,intent(in) :: Rmin
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: qlR1,qlR2,qlI1,qlI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: quR1,quR2,quI1,quI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: gR,gI
!  gont = 0: success; gont = 1: something wrong
  integer::gont

  double complex,dimension(ex(1),ex(2),ex(3)) :: CU,DCU,bDCU,CB,DCB,bDCB,CJ,DCJ
  double complex :: CTheta0,CTheta,CTheta1,RHS
  integer :: i,j,k,RK4
  double complex,dimension(ex(3)) :: Cnu,Ck,HCnu,HCk,HCU,HDCU,CUx,HCUx,DCUx,HDCUx,HbDCU,bDCUx,HbDCUx
  double complex,dimension(ex(3)) :: Cnux,HCnux,HCJ,HDCJ,CJx,HCJx,CJxx,HCJxx,DCJx,HDCJx,HCB,HDCB,HbDCB
  real*8,dimension(ex(3)) :: Hbeta,betax,Hbetax,HW,Wx,HWx
  double complex :: Theta_rhs,Theta_rhs_o
  real*8 :: xx,dR

!!! sanity check
  dR = sum(RJ)+sum(IJ)+sum(RU)+sum(IU)+sum(beta)+sum(RB)+sum(IB) + &
       sum(Rnu)+sum(Inu)+sum(Rk)+sum(Ik)+sum(RTheta)+sum(ITheta) + &
       sum(KK)+sum(HKK)+sum(KKx)+sum(HKKx)
  if(dR.ne.dR) then
     if(sum(RJ).ne.sum(RJ))write(*,*)"NullEvol_Theta: find NaN in RJ"
     if(sum(IJ).ne.sum(IJ))write(*,*)"NullEvol_Theta: find NaN in IJ"
     if(sum(RU).ne.sum(RU))write(*,*)"NullEvol_Theta: find NaN in RU"
     if(sum(IU).ne.sum(IU))write(*,*)"NullEvol_Theta: find NaN in IU"
     if(sum(beta).ne.sum(beta))write(*,*)"NullEvol_Theta: find NaN in beta"
     if(sum(RB).ne.sum(RB))write(*,*)"NullEvol_Theta: find NaN in RB"
     if(sum(IB).ne.sum(IB))write(*,*)"NullEvol_Theta: find NaN in IB"
     if(sum(Rnu).ne.sum(Rnu))write(*,*)"NullEvol_Theta: find NaN in Rnu"
     if(sum(Inu).ne.sum(Inu))write(*,*)"NullEvol_Theta: find NaN in Inu"
     if(sum(Rk).ne.sum(Rk))write(*,*)"NullEvol_Theta: find NaN in Rk"
     if(sum(Ik).ne.sum(Ik))write(*,*)"NullEvol_Theta: find NaN in Ik"
     if(sum(RTheta).ne.sum(RTheta))write(*,*)"NullEvol_Theta: find NaN in RTheta"
     if(sum(ITheta).ne.sum(ITheta))write(*,*)"NullEvol_Theta: find NaN in ITheta"
     if(sum(KK).ne.sum(KK))write(*,*)"NullEvol_Theta: find NaN in KK"
     if(sum(HKK).ne.sum(HKK))write(*,*)"NullEvol_Theta: find NaN in HKK"
     if(sum(KKx).ne.sum(KKx))write(*,*)"NullEvol_Theta: find NaN in KKx"
     if(sum(HKKx).ne.sum(HKKx))write(*,*)"NullEvol_Theta: find NaN in HKKx"
     gont = 1
     return
  endif

  dR = R(2) - R(1)
  
  CU = dcmplx(RU,IU)
  CB = dcmplx(RB,IB)
  CJ = dcmplx(RJ,IJ)

  do k=1,ex(3)
     call derivs_eth(ex(1:2),crho,sigma,CU(:,:,k),DCU(:,:,k),1,1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CU(:,:,k),bDCU(:,:,k),1,-1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CB(:,:,k),DCB(:,:,k),1,1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CB(:,:,k),bDCB(:,:,k),1,-1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CJ(:,:,k),DCJ(:,:,k),2,1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
  enddo

  do j=1,ex(2)
  do i=1,ex(1)
     CTheta0 = dcmplx(RTheta(i,j,1),ITheta(i,j,1))
     Cnu = dcmplx(Rnu(i,j,:),Inu(i,j,:))
     Ck = dcmplx(Rk(i,j,:),Ik(i,j,:))
     call cget_half_x(ex(3),CB(i,j,:),HCB)
     call cget_half_x(ex(3),DCB(i,j,:),HDCB)
     call cget_half_x(ex(3),bDCB(i,j,:),HbDCB)
     call cget_half_x(ex(3),Cnu,HCnu)
     call cderivs_x(ex(3),R,Cnu,Cnux)
     call cget_half_x(ex(3),Cnux,HCnux)
     call cget_half_x(ex(3),Ck,HCk)
     call rget_half_x(ex(3),beta(i,j,:),Hbeta)
     call rderivs_x(ex(3),R,beta(i,j,:),betax)
     call rget_half_x(ex(3),betax,Hbetax)
     call rderivs_x(ex(3),R,W,Wx)
     call rget_half_x(ex(3),Wx,HWx)
     call rget_half_x(ex(3),W(i,j,:),HW)
     call cget_half_x(ex(3),CU(i,j,:),HCU)
     call cderivs_x(ex(3),R,DCU(i,j,:),DCUx)
     call cderivs_x(ex(3),R,CU(i,j,:),CUx)
     call cget_half_x(ex(3),DCUx,HDCUx)
     call cget_half_x(ex(3),CUx,HCUx)
     call cget_half_x(ex(3),DCU(i,j,:),HDCU)
     call cderivs_x(ex(3),R,bDCU(i,j,:),bDCUx)
     call cget_half_x(ex(3),bDCUx,HbDCUx)
     call cget_half_x(ex(3),bDCU(i,j,:),HbDCU)
     call cderivs_x(ex(3),R,CJ(i,j,:),CJx)
     call cdderivs_x(ex(3),R,CJ(i,j,:),CJxx)
     call cget_half_x(ex(3),CJx,HCJx)
     call cget_half_x(ex(3),CJxx,HCJxx)
     call cderivs_x(ex(3),R,DCJ(i,j,:),DCJx)
     call cget_half_x(ex(3),DCJx,HDCJx)
! old type code: PRD 54, 6153, Eq.(32) etc.
#if 0
! start up part
     k = 1
     RHS = Theta_rhs_o(R(k)+dR/2.d0,Rmin,Hbeta(k),betax(k),HKK(i,j,k),KKx(i,j,k),HCU(k),CUx(k),DCUx(k),HbDCU(k),bDCUx(k), &
                        HDCU(k),HCB(k),HDCB(k),HW(k),Wx(k),HCJ(k),HDCJ(k),    &
                        CJx(k),CJxx(k),DCJx(k),HbDCB(k),HCnu(k),Cnux(k),HCk(k),CTheta0)
     CTheta1 = RHS-(1-2.d0*(R(k)+dR/2.d0)*(1.d0-R(k)-dR/2.d0)/dR)*CTheta0
     CTheta0 = CTheta1/(1+2.d0*(R(k)+dR/2.d0)*(1.d0-R(k)-dR/2.d0)/dR)

     RTheta(i,j,k+1) = dreal(CTheta0)
     ITheta(i,j,k+1) = dimag(CTheta0)

     do k=1,ex(3)-2
        RHS = Theta_rhs_o(R(k+1),Rmin,beta(i,j,k+1),betax(k+1),KK(i,j,k+1),KKx(i,j,k+1),CU(i,j,k+1),CUx(k+1),DCUx(k+1),bDCU(i,j,k+1),bDCUx(k+1), &
                        DCU(i,j,k+1),CB(i,j,k+1),DCB(i,j,k+1),W(i,j,k+1),Wx(k+1),CJ(i,j,k+1),DCJ(i,j,k+1),    &
                        CJx(k+1),CJxx(k+1),DCJx(k+1),bDCB(i,j,k+1),Cnu(k+1),Cnux(k+1),Ck(k+1),CTheta0)
        CTheta1 = RHS-(1-R(k+1)*(1.d0-R(k+1))/dR)*(dcmplx(RTheta(i,j,k),ITheta(i,j,k)))
        CTheta0 = CTheta1/(1+R(k+1)*(1.d0-R(k+1))/dR)

        RTheta(i,j,k+2) = dreal(CTheta0)
        ITheta(i,j,k+2) = dimag(CTheta0)
     enddo
#endif     

#ifdef OLD
     do k = 1,ex(3)-1
       xx = R(k)+dR/2
! note our HCJ(ex(3)-1) = (CJ(ex(3))+CJ(ex(3)-1))/2
! note our KKx(ex(3)) = (KK(ex(3))-KK(ex(3)-1))/dR
! note our fxx(ex(3)) = (f(ex(3))-2.d0*f(ex(3)-1)+f(ex(3)-2))/dR
       RHS = Theta_rhs(xx,Rmin,Hbeta(k),betax(k+1),HKK(i,j,k),KKx(i,j,k+1),HCU(k),CUx(k+1),DCUx(k+1),HbDCU(k),bDCUx(k+1), &
                     HDCU(k),HCB(k),HDCB(k),HW(k),Wx(k+1),HCJ(k),HDCJ(k),    &
                     CJx(k+1),CJxx(k+1),DCJx(k+1),HbDCB(k),HCnu(k),Cnux(k+1),HCk(k),0)
       RHS = RHS+CTheta0*(1.d0/dR-0.5d0/xx/(1.d0-xx))                 
       CTheta0 = RHS/(1.d0/dR+0.5d0/xx/(1.d0-xx)) 
       RTheta(i,j,k+1) = dreal(CTheta0)
       ITheta(i,j,k+1) = dimag(CTheta0)
     enddo
#else
     do k=1,ex(3)-2
        RK4 = 0
        RHS = Theta_rhs(R(k),Rmin,beta(i,j,k),betax(k),KK(i,j,k),KKx(i,j,k),CU(i,j,k),CUx(k),DCUx(k),bDCU(i,j,k),bDCUx(k), &
                        DCU(i,j,k),CB(i,j,k),DCB(i,j,k),W(i,j,k),Wx(k),CJ(i,j,k),DCJ(i,j,k),    &
                        CJx(k),CJxx(k),DCJx(k),bDCB(i,j,k),Cnu(k),Cnux(k),Ck(k),CTheta0)
        call rungekutta4_cplxscalar(dR,CTheta0,CTheta,RHS,RK4)

        RK4 = 1
        CTheta1 = Theta_rhs(R(k)+dR/2.d0,Rmin,Hbeta(k),Hbetax(k),HKK(i,j,k),HKKx(i,j,k), &
                        HCU(k),HCUx(k),HDCUx(k),HbDCU(k),HbDCUx(k), &
                        HDCU(k),HCB(k),HDCB(k),HW(k),HWx(k),HCJ(k),HDCJ(k),    &
                        HCJx(k),HCJxx(k),HDCJx(k),HbDCB(k),HCnu(k),HCnux(k),HCk(k),CTheta)
        call rungekutta4_cplxscalar(dR,CTheta0,CTheta1,RHS,RK4)
        call cswap(CTheta,CTheta1)

        RK4 = 2
        CTheta1 = Theta_rhs(R(k)+dR/2.d0,Rmin,Hbeta(k),Hbetax(k),HKK(i,j,k),HKKx(i,j,k), &
                        HCU(k),HCUx(k),HDCUx(k),HbDCU(k),HbDCUx(k), &
                        HDCU(k),HCB(k),HDCB(k),HW(k),HWx(k),HCJ(k),HDCJ(k),    &
                        HCJx(k),HCJxx(k),HDCJx(k),HbDCB(k),HCnu(k),HCnux(k),HCk(k),CTheta)
        call rungekutta4_cplxscalar(dR,CTheta0,CTheta1,RHS,RK4)
        call cswap(CTheta,CTheta1)
        
        RK4 = 3
        CTheta1 = Theta_rhs(R(k+1),Rmin,beta(i,j,k+1),betax(k+1),KK(i,j,k+1),KKx(i,j,k+1), &
                        CU(i,j,k+1),CUx(k+1),DCUx(k+1),bDCU(i,j,k+1),bDCUx(k+1), &
                        DCU(i,j,k+1),CB(i,j,k+1),DCB(i,j,k+1),W(i,j,k+1),Wx(k+1),CJ(i,j,k+1),DCJ(i,j,k+1),    &
                        CJx(k+1),CJxx(k+1),DCJx(k+1),bDCB(i,j,k+1),Cnu(k+1),Cnux(k+1),Ck(k+1),CTheta)
        call rungekutta4_cplxscalar(dR,CTheta0,CTheta1,RHS,RK4)
        call cswap(CTheta0,CTheta1)

        RTheta(i,j,k+1) = dreal(CTheta0)
        ITheta(i,j,k+1) = dimag(CTheta0)
     enddo
#if 0
     k = ex(3)
     CTheta0 = -KK(i,j,k)*DCU(i,j,k)-(CU(i,j,k)*Cnu(k)+dconjg(CU(i,j,k))*DCJ(i,j,k))/2 &
                  +CJ(i,j,k)*(bDCU(i,j,k)-dconjg(bDCU(i,j,k)))/2 - W(i,j,k)*CJ(i,j,k)/2

     RTheta(i,j,k) = dreal(CTheta0)
     ITheta(i,j,k) = dimag(CTheta0)
#else
! closing step
     k = ex(3)-1
     CTheta0 = dcmplx(RTheta(i,j,k),ITheta(i,j,k))
     xx = R(k)+dR/2
! note our HCJ(ex(3)-1) = (CJ(ex(3))+CJ(ex(3)-1))/2
! note our KKx(ex(3)) = (KK(ex(3))-KK(ex(3)-1))/dR
! note our fxx(ex(3)) = (f(ex(3))-2.d0*f(ex(3)-1)+f(ex(3)-2))/dR
     RHS = Theta_rhs(xx,Rmin,Hbeta(k),betax(k+1),HKK(i,j,k),KKx(i,j,k+1),HCU(k),CUx(k+1),DCUx(k+1),HbDCU(k),bDCUx(k+1), &
                     HDCU(k),HCB(k),HDCB(k),HW(k),Wx(k+1),HCJ(k),HDCJ(k),    &
                     CJx(k+1),CJxx(k+1),DCJx(k+1),HbDCB(k),HCnu(k),Cnux(k+1),HCk(k),dcmplx(0.d0,0.d0))
     RHS = RHS+CTheta0*(1.d0/dR-0.5d0/xx/(1.d0-xx))                 
     CTheta0 = RHS/(1.d0/dR+0.5d0/xx/(1.d0-xx)) 
     RTheta(i,j,k+1) = dreal(CTheta0)
     ITheta(i,j,k+1) = dimag(CTheta0)
#endif

#endif
  enddo
  enddo

  gont = 0
  return

end function NullEvol_Theta
!--------------------------------------------------------------------
! check with fake_Theta_rhs
#elif 0
! this R is indeed x
function NullEvol_Theta(ex,crho,sigma,R,RJ,IJ,RU,IU,beta,RB,IB, &
                        Rnu,Inu,Rk,Ik,RTheta,ITheta,W,KK,HKK,KKx,HKKx,Rmin,       &
                        qlR1,qlR2,qlI1,qlI2,quR1,quR2,quI1,quI2,gR,gI) result(gont)   
  implicit none
  integer,intent(in ):: ex(1:3)
  real*8,intent(in),dimension(ex(1))::crho
  real*8,intent(in),dimension(ex(2))::sigma
  real*8,intent(in),dimension(ex(3))::R
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: beta,W,KK,HKK,KKx,HKKx
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: RJ,IJ,RU,IU,RB,IB,Rnu,Inu,Rk,Ik
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: RTheta,ITheta
  real*8,intent(in) :: Rmin
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: qlR1,qlR2,qlI1,qlI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: quR1,quR2,quI1,quI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: gR,gI
!  gont = 0: success; gont = 1: something wrong
  integer::gont

  double complex,dimension(ex(1),ex(2),ex(3)) :: CU,DCU,bDCU,CB,DCB,bDCB,CJ,DCJ
  double complex :: CTheta0,CTheta,CTheta1,RHS
  integer :: i,j,k,RK4
  integer,parameter :: ks=1
  double complex,dimension(ex(3)) :: Cnu,Ck,HCnu,HCk,HCU,HDCU,CUx,HCUx,DCUx,HDCUx,HbDCU,bDCUx,HbDCUx
  double complex,dimension(ex(3)) :: Cnux,HCnux,HCJ,HDCJ,CJx,HCJx,CJxx,HCJxx,DCJx,HDCJx,HCB,HDCB,HbDCB
  real*8,dimension(ex(3)) :: Hbeta,betax,Hbetax,HW,Wx,HWx
  double complex :: Theta_rhs,Theta_rhs_o
  real*8 :: xx,dR

!!! sanity check
  dR = sum(RJ)+sum(IJ)+sum(RU)+sum(IU)+sum(beta)+sum(RB)+sum(IB) + &
       sum(Rnu)+sum(Inu)+sum(Rk)+sum(Ik)+sum(RTheta)+sum(ITheta) + &
       sum(KK)+sum(HKK)+sum(KKx)+sum(HKKx)
  if(dR.ne.dR) then
     if(sum(RJ).ne.sum(RJ))write(*,*)"NullEvol_Theta: find NaN in RJ"
     if(sum(IJ).ne.sum(IJ))write(*,*)"NullEvol_Theta: find NaN in IJ"
     if(sum(RU).ne.sum(RU))write(*,*)"NullEvol_Theta: find NaN in RU"
     if(sum(IU).ne.sum(IU))write(*,*)"NullEvol_Theta: find NaN in IU"
     if(sum(beta).ne.sum(beta))write(*,*)"NullEvol_Theta: find NaN in beta"
     if(sum(RB).ne.sum(RB))write(*,*)"NullEvol_Theta: find NaN in RB"
     if(sum(IB).ne.sum(IB))write(*,*)"NullEvol_Theta: find NaN in IB"
     if(sum(Rnu).ne.sum(Rnu))write(*,*)"NullEvol_Theta: find NaN in Rnu"
     if(sum(Inu).ne.sum(Inu))write(*,*)"NullEvol_Theta: find NaN in Inu"
     if(sum(Rk).ne.sum(Rk))write(*,*)"NullEvol_Theta: find NaN in Rk"
     if(sum(Ik).ne.sum(Ik))write(*,*)"NullEvol_Theta: find NaN in Ik"
     if(sum(RTheta).ne.sum(RTheta))write(*,*)"NullEvol_Theta: find NaN in RTheta"
     if(sum(ITheta).ne.sum(ITheta))write(*,*)"NullEvol_Theta: find NaN in ITheta"
     if(sum(KK).ne.sum(KK))write(*,*)"NullEvol_Theta: find NaN in KK"
     if(sum(HKK).ne.sum(HKK))write(*,*)"NullEvol_Theta: find NaN in HKK"
     if(sum(KKx).ne.sum(KKx))write(*,*)"NullEvol_Theta: find NaN in KKx"
     if(sum(HKKx).ne.sum(HKKx))write(*,*)"NullEvol_Theta: find NaN in HKKx"
     gont = 1
     return
  endif

  dR = R(2) - R(1)
  
  CU = dcmplx(RU,IU)
  CB = dcmplx(RB,IB)
  CJ = dcmplx(RJ,IJ)

  do k=1,ex(3)
     call derivs_eth(ex(1:2),crho,sigma,CU(:,:,k),DCU(:,:,k),1,1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CU(:,:,k),bDCU(:,:,k),1,-1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CB(:,:,k),DCB(:,:,k),1,1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CB(:,:,k),bDCB(:,:,k),1,-1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CJ(:,:,k),DCJ(:,:,k),2,1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
  enddo

  do j=1,ex(2)
  do i=1,ex(1)
     CTheta0 = dcmplx(RTheta(i,j,ks),ITheta(i,j,ks))
     Cnu = dcmplx(RTheta(i,j,:),ITheta(i,j,:))
     call fake_Theta_rhs(ex(3),R,Ck,Cnu)
     call cget_half_x(ex(3),Ck,HCk)

     do k=ks,ex(3)-1
        RK4 = 0
        RHS = Ck(k)
        call rungekutta4_cplxscalar(dR,CTheta0,CTheta,RHS,RK4)

        RK4 = 1
        CTheta1 = HCk(k)
        call rungekutta4_cplxscalar(dR,CTheta0,CTheta1,RHS,RK4)
        call cswap(CTheta,CTheta1)

        RK4 = 2
        CTheta1 = HCk(k)
        call rungekutta4_cplxscalar(dR,CTheta0,CTheta1,RHS,RK4)
        call cswap(CTheta,CTheta1)
        
        RK4 = 3
        CTheta1 = Ck(k+1)
        call rungekutta4_cplxscalar(dR,CTheta0,CTheta1,RHS,RK4)
        call cswap(CTheta0,CTheta1)

        RTheta(i,j,k+1) = dreal(CTheta0)
        ITheta(i,j,k+1) = dimag(CTheta0)
     enddo

  enddo
  enddo

  gont = 0
  return

end function NullEvol_Theta

#else
! for eth_x, _x first, eth second
! this R is indeed x
function NullEvol_Theta(ex,crho,sigma,R,RJ,IJ,RU,IU,beta,RB,IB, &
                        Rnu,Inu,Rk,Ik,RTheta,ITheta,W,KK,HKK,KKx,HKKx,Rmin,       &
                        qlR1,qlR2,qlI1,qlI2,quR1,quR2,quI1,quI2,gR,gI) result(gont)   
  implicit none
  integer,intent(in ):: ex(1:3)
  real*8,intent(in),dimension(ex(1))::crho
  real*8,intent(in),dimension(ex(2))::sigma
  real*8,intent(in),dimension(ex(3))::R
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: beta,W,KK,HKK,KKx,HKKx
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: RJ,IJ,RU,IU,RB,IB,Rnu,Inu,Rk,Ik
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: RTheta,ITheta
  real*8,intent(in) :: Rmin
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: qlR1,qlR2,qlI1,qlI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: quR1,quR2,quI1,quI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: gR,gI
!  gont = 0: success; gont = 1: something wrong
  integer::gont

  double complex,dimension(ex(1),ex(2),ex(3)) :: CU,DCU,bDCU,CB,DCB,bDCB,CJ,DCJ
  double complex :: CTheta0,CTheta,CTheta1,RHS
  integer :: i,j,k,RK4
  double complex,dimension(ex(1),ex(2),ex(3)) :: Cnu,Ck,HCnu,HCk,HCU,HDCU,CUx,HCUx,DCUx,HDCUx,HbDCU,bDCUx,HbDCUx
  double complex,dimension(ex(1),ex(2),ex(3)) :: Cnux,HCnux,HCJ,HDCJ,CJx,HCJx,CJxx,HCJxx,DCJx,HDCJx,HCB,HDCB,HbDCB
  real*8,dimension(ex(1),ex(2),ex(3)) :: Hbeta,betax,Hbetax,HW,Wx,HWx
  double complex :: Theta_rhs
  real*8 :: xx,dR

!!! sanity check
  dR = sum(RJ)+sum(IJ)+sum(RU)+sum(IU)+sum(beta)+sum(RB)+sum(IB) + &
       sum(Rnu)+sum(Inu)+sum(Rk)+sum(Ik)+sum(RTheta)+sum(ITheta) + &
       sum(KK)+sum(HKK)+sum(KKx)+sum(HKKx)
  if(dR.ne.dR) then
     if(sum(RJ).ne.sum(RJ))write(*,*)"NullEvol_Theta: find NaN in RJ"
     if(sum(IJ).ne.sum(IJ))write(*,*)"NullEvol_Theta: find NaN in IJ"
     if(sum(RU).ne.sum(RU))write(*,*)"NullEvol_Theta: find NaN in RU"
     if(sum(IU).ne.sum(IU))write(*,*)"NullEvol_Theta: find NaN in IU"
     if(sum(beta).ne.sum(beta))write(*,*)"NullEvol_Theta: find NaN in beta"
     if(sum(RB).ne.sum(RB))write(*,*)"NullEvol_Theta: find NaN in RB"
     if(sum(IB).ne.sum(IB))write(*,*)"NullEvol_Theta: find NaN in IB"
     if(sum(Rnu).ne.sum(Rnu))write(*,*)"NullEvol_Theta: find NaN in Rnu"
     if(sum(Inu).ne.sum(Inu))write(*,*)"NullEvol_Theta: find NaN in Inu"
     if(sum(Rk).ne.sum(Rk))write(*,*)"NullEvol_Theta: find NaN in Rk"
     if(sum(Ik).ne.sum(Ik))write(*,*)"NullEvol_Theta: find NaN in Ik"
     if(sum(RTheta).ne.sum(RTheta))write(*,*)"NullEvol_Theta: find NaN in RTheta"
     if(sum(ITheta).ne.sum(ITheta))write(*,*)"NullEvol_Theta: find NaN in ITheta"
     if(sum(KK).ne.sum(KK))write(*,*)"NullEvol_Theta: find NaN in KK"
     if(sum(HKK).ne.sum(HKK))write(*,*)"NullEvol_Theta: find NaN in HKK"
     if(sum(KKx).ne.sum(KKx))write(*,*)"NullEvol_Theta: find NaN in KKx"
     if(sum(HKKx).ne.sum(HKKx))write(*,*)"NullEvol_Theta: find NaN in HKKx"
     gont = 1
     return
  endif

  dR = R(2) - R(1)
  
  CU = dcmplx(RU,IU)
  CB = dcmplx(RB,IB)
  CJ = dcmplx(RJ,IJ)
  Cnu = dcmplx(Rnu,Inu)
  Ck = dcmplx(Rk,Ik)

  do j=1,ex(2)
  do i=1,ex(1)
     call cderivs_x(ex(3),R,Cnu(i,j,:),Cnux(i,j,:))
     call rderivs_x(ex(3),R,beta(i,j,:),betax(i,j,:))
     call rderivs_x(ex(3),R,W(i,j,:),Wx(i,j,:))
     call cderivs_x(ex(3),R,CU(i,j,:),CUx(i,j,:))
     call cderivs_x(ex(3),R,CJ(i,j,:),CJx(i,j,:))
     call cdderivs_x(ex(3),R,CJ(i,j,:),CJxx(i,j,:))
  enddo
  enddo

  do k=1,ex(3)
     call derivs_eth(ex(1:2),crho,sigma,CU(:,:,k),DCU(:,:,k),1,1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CUx(:,:,k),DCUx(:,:,k),1,1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CU(:,:,k),bDCU(:,:,k),1,-1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CUx(:,:,k),bDCUx(:,:,k),1,-1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CB(:,:,k),DCB(:,:,k),1,1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CB(:,:,k),bDCB(:,:,k),1,-1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CJ(:,:,k),DCJ(:,:,k),2,1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CJx(:,:,k),DCJx(:,:,k),2,1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
  enddo

  do j=1,ex(2)
  do i=1,ex(1)
     call cget_half_x(ex(3),CB(i,j,:),HCB(i,j,:))
     call cget_half_x(ex(3),DCB(i,j,:),HDCB(i,j,:))
     call cget_half_x(ex(3),bDCB(i,j,:),HbDCB(i,j,:))
     call cget_half_x(ex(3),Cnu(i,j,:),HCnu(i,j,:))
     call cget_half_x(ex(3),Cnux(i,j,:),HCnux(i,j,:))
     call cget_half_x(ex(3),Ck(i,j,:),HCk(i,j,:))
     call rget_half_x(ex(3),beta(i,j,:),Hbeta(i,j,:))
     call rget_half_x(ex(3),betax(i,j,:),Hbetax(i,j,:))
     call rget_half_x(ex(3),Wx(i,j,:),HWx(i,j,:))
     call rget_half_x(ex(3),W(i,j,:),HW(i,j,:))
     call cget_half_x(ex(3),CU(i,j,:),HCU(i,j,:))
     call cget_half_x(ex(3),DCUx(i,j,:),HDCUx(i,j,:))
     call cget_half_x(ex(3),CUx(i,j,:),HCUx(i,j,:))
     call cget_half_x(ex(3),DCU(i,j,:),HDCU(i,j,:))
     call cget_half_x(ex(3),bDCUx(i,j,:),HbDCUx(i,j,:))
     call cget_half_x(ex(3),bDCU(i,j,:),HbDCU(i,j,:))
     call cget_half_x(ex(3),CJx(i,j,:),HCJx(i,j,:))
     call cget_half_x(ex(3),CJxx(i,j,:),HCJxx(i,j,:))
     call cget_half_x(ex(3),DCJx(i,j,:),HDCJx(i,j,:))
  enddo
  enddo

  do j=1,ex(2)
  do i=1,ex(1)
     CTheta0 = dcmplx(RTheta(i,j,1),ITheta(i,j,1))

     do k=1,ex(3)-2
        RK4 = 0
        RHS = Theta_rhs(R(k),Rmin,beta(i,j,k),betax(i,j,k),KK(i,j,k),KKx(i,j,k),CU(i,j,k),CUx(i,j,k),DCUx(i,j,k),bDCU(i,j,k),bDCUx(i,j,k), &
                        DCU(i,j,k),CB(i,j,k),DCB(i,j,k),W(i,j,k),Wx(i,j,k),CJ(i,j,k),DCJ(i,j,k),    &
                        CJx(i,j,k),CJxx(i,j,k),DCJx(i,j,k),bDCB(i,j,k),Cnu(i,j,k),Cnux(i,j,k),Ck(i,j,k),CTheta0)
        call rungekutta4_cplxscalar(dR,CTheta0,CTheta,RHS,RK4)

        RK4 = 1
        CTheta1 = Theta_rhs(R(k)+dR/2.d0,Rmin,Hbeta(i,j,k),Hbetax(i,j,k),HKK(i,j,k),HKKx(i,j,k), &
                   HCU(i,j,k),HCUx(i,j,k),HDCUx(i,j,k),HbDCU(i,j,k),HbDCUx(i,j,k), &
                        HDCU(i,j,k),HCB(i,j,k),HDCB(i,j,k),HW(i,j,k),HWx(i,j,k),HCJ(i,j,k),HDCJ(i,j,k),    &
                        HCJx(i,j,k),HCJxx(i,j,k),HDCJx(i,j,k),HbDCB(i,j,k),HCnu(i,j,k),HCnux(i,j,k),HCk(i,j,k),CTheta)
        call rungekutta4_cplxscalar(dR,CTheta0,CTheta1,RHS,RK4)
        call cswap(CTheta,CTheta1)

        RK4 = 2
        CTheta1 = Theta_rhs(R(k)+dR/2.d0,Rmin,Hbeta(i,j,k),Hbetax(i,j,k),HKK(i,j,k),HKKx(i,j,k), &
                   HCU(i,j,k),HCUx(i,j,k),HDCUx(i,j,k),HbDCU(i,j,k),HbDCUx(i,j,k), &
                        HDCU(i,j,k),HCB(i,j,k),HDCB(i,j,k),HW(i,j,k),HWx(i,j,k),HCJ(i,j,k),HDCJ(i,j,k),    &
                        HCJx(i,j,k),HCJxx(i,j,k),HDCJx(i,j,k),HbDCB(i,j,k),HCnu(i,j,k),HCnux(i,j,k),HCk(i,j,k),CTheta)
        call rungekutta4_cplxscalar(dR,CTheta0,CTheta1,RHS,RK4)
        call cswap(CTheta,CTheta1)
        
        RK4 = 3
        CTheta1 = Theta_rhs(R(k+1),Rmin,beta(i,j,k+1),betax(i,j,k+1),KK(i,j,k+1), &
                 KKx(i,j,k+1),CU(i,j,k+1),CUx(i,j,k+1),DCUx(i,j,k+1),bDCU(i,j,k+1),bDCUx(i,j,k+1), &
                        DCU(i,j,k+1),CB(i,j,k+1),DCB(i,j,k+1),W(i,j,k+1),Wx(i,j,k+1),CJ(i,j,k+1),DCJ(i,j,k+1),    &
                        CJx(i,j,k+1),CJxx(i,j,k+1),DCJx(i,j,k+1),bDCB(i,j,k+1),Cnu(i,j,k+1),Cnux(i,j,k+1),Ck(i,j,k+1),CTheta)
        call rungekutta4_cplxscalar(dR,CTheta0,CTheta1,RHS,RK4)
        call cswap(CTheta0,CTheta1)

        RTheta(i,j,k+1) = dreal(CTheta0)
        ITheta(i,j,k+1) = dimag(CTheta0)
     enddo

     k = ex(3)
     CTheta0 = -KK(i,j,k)*DCU(i,j,k)-(CU(i,j,k)*Cnu(i,j,k)+dconjg(CU(i,j,k))*DCJ(i,j,k))/2 &
                  +CJ(i,j,k)*(bDCU(i,j,k)-dconjg(bDCU(i,j,k)))/2 - W(i,j,k)*CJ(i,j,k)/2

     RTheta(i,j,k) = dreal(CTheta0)
     ITheta(i,j,k) = dimag(CTheta0)

  enddo
  enddo

  gont = 0
  return

end function NullEvol_Theta
#endif

#elif (RKorAM == 1)
!------------------------------------------------------------------------------
! this R is indeed x
function NullEvol_beta(ex,crho,sigma,R,RJ,IJ,beta,KKx,HKKx)    result(gont)   
  implicit none
  integer,intent(in ):: ex(1:3)
  real*8,intent(in),dimension(ex(1))::crho
  real*8,intent(in),dimension(ex(2))::sigma
  real*8,intent(in),dimension(ex(3))::R
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: beta
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: RJ,IJ,KKx,HKKx
!  gont = 0: success; gont = 1: something wrong
  integer::gont
  real*8 :: dR,beta_rhs

  double complex, dimension(ex(3)):: CJ,CJx
  real*8, dimension(ex(3)) :: rhs
  integer :: i,j,k

  real*8,parameter :: F5o12=2.d0/1.2d1,F2o3=2.d0/3.d0,F1o12=1.d0/1.2d1
  real*8,parameter :: F3o8=3.d0/8.d0,F19o24=1.9d1/2.4d1,F5o24=5.d0/2.4d1,F1o24=1.d0/2.4d1

!!! sanity check
  dR = sum(RJ)+sum(IJ)+sum(beta)+sum(KKx)+sum(HKKx)
  if(dR.ne.dR) then
     if(sum(RJ).ne.sum(RJ))write(*,*)"NullEvol_beta: find NaN in RJ"
     if(sum(IJ).ne.sum(IJ))write(*,*)"NullEvol_beta: find NaN in IJ"
     if(sum(beta).ne.sum(beta))write(*,*)"NullEvol_beta: find NaN in beta"
     if(sum(KKx).ne.sum(KKx))write(*,*)"NullEvol_beta: find NaN in KKx"
     if(sum(HKKx).ne.sum(HKKx))write(*,*)"NullEvol_beta: find NaN in HKKx"
     gont = 1
     return
  endif

  dR = R(2) - R(1)

  do j=1,ex(2)
  do i=1,ex(1)
     CJ = dcmplx(RJ(i,j,:),IJ(i,j,:))
#if 0
     call cderivs_sw_x(ex(3),R,CJ,CJx)
#else
     call cderivs_x(ex(3),R,CJ,CJx)
#endif

     do k=1,ex(3)
        rhs(k) = beta_rhs(R(k),CJx(k),KKx(i,j,k)) 
     enddo

     k = 1
     beta(i,j,k+1) = beta(i,j,k) + (rhs(k+1)+rhs(k))*dR/2

     k = 2
     beta(i,j,k+1) = beta(i,j,k) + (F5o12*rhs(k+1) + F2o3*rhs(k) - F1o12*rhs(k-1))*dR

     do k=3,ex(3)-1
        beta(i,j,k+1) = beta(i,j,k) + (F3o8*rhs(k+1) + F19o24*rhs(k) - F5o24*rhs(k-1) + F1o24*rhs(k-2))*dR
     enddo

  enddo
  enddo

  gont = 0

  return

end function NullEvol_beta
!------------------------------------------------------------------------------
! this R is indeed x
function NullEvol_Q(ex,crho,sigma,R,RJ,IJ,Rk,Ik,Rnu,Inu,RB,IB,RQ,IQ,KK,Hkk,KKx,HKKx, &
                    qlR1,qlR2,qlI1,qlI2,quR1,quR2,quI1,quI2,gR,gI)    result(gont)   
  implicit none
  integer,intent(in ):: ex(1:3)
  real*8,intent(in),dimension(ex(1))::crho
  real*8,intent(in),dimension(ex(2))::sigma
  real*8,intent(in),dimension(ex(3))::R
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: RQ,IQ
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: RJ,IJ,KK,KKx,HKK,HKKx
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: Rk,Ik,Rnu,Inu,RB,IB
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: qlR1,qlR2,qlI1,qlI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: quR1,quR2,quI1,quI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: gR,gI
!  gont = 0: success; gont = 1: something wrong
  integer::gont
  real*8 :: xx,dR

  double complex,dimension(ex(3)) :: CQ,RHS
  real*8, dimension(ex(3)) :: gunc
  double complex,dimension(ex(3)) :: CJx,DCJx,Ck,Ckx,Cnu,Cnux,CB,CBx
  double complex, dimension(ex(1),ex(2),ex(3)) :: CJ,DCJ
  integer :: i,j,k
  double complex :: ZEO,Q_rhs

  real*8,parameter :: F5o12=2.d0/1.2d1,F2o3=2.d0/3.d0,F1o12=1.d0/1.2d1
  real*8,parameter :: F3o8=3.d0/8.d0,F19o24=1.9d1/2.4d1,F5o24=5.d0/2.4d1,F1o24=1.d0/2.4d1

!!! sanity check
  dR = sum(RJ)+sum(IJ)+sum(RQ)+sum(IQ) &
      +sum(RK)+sum(IK)+sum(Rnu)+sum(Inu)+sum(RB)+sum(IB) &
      +sum(KK)+sum(KKx)
  if(dR.ne.dR) then
     if(sum(RJ).ne.sum(RJ))write(*,*)"NullEvol_Q: find NaN in RJ"
     if(sum(IJ).ne.sum(IJ))write(*,*)"NullEvol_Q: find NaN in IJ"
     if(sum(RQ).ne.sum(RQ))write(*,*)"NullEvol_Q: find NaN in RQ"
     if(sum(IQ).ne.sum(IQ))write(*,*)"NullEvol_Q: find NaN in IQ"
     if(sum(RK).ne.sum(RK))write(*,*)"NullEvol_Q: find NaN in RK"
     if(sum(IK).ne.sum(IK))write(*,*)"NullEvol_Q: find NaN in IK"
     if(sum(Rnu).ne.sum(Rnu))write(*,*)"NullEvol_Q: find NaN in Rnu"
     if(sum(Inu).ne.sum(Inu))write(*,*)"NullEvol_Q: find NaN in Inu"
     if(sum(RB).ne.sum(RB))write(*,*)"NullEvol_Q: find NaN in RB"
     if(sum(IB).ne.sum(IB))write(*,*)"NullEvol_Q: find NaN in IB"
     if(sum(KK).ne.sum(KK))write(*,*)"NullEvol_Q: find NaN in KK"
     if(sum(KKx).ne.sum(KKx))write(*,*)"NullEvol_Q: find NaN in KKx"
     gont = 1
     return
  endif

  dR = R(2) - R(1)
  ZEO = dcmplx(0.d0,0.d0)

  CJ = dcmplx(RJ,IJ)
  do k=1,ex(3)
     call derivs_eth(ex(1:2),crho,sigma,CJ(:,:,k),DCJ(:,:,k),2,1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
  enddo
  do j=1,ex(2)
  do i=1,ex(1)

     CQ(1) = dcmplx(RQ(i,j,1),IQ(i,j,1))
     Ck = dcmplx(Rk(i,j,:),Ik(i,j,:))
     Cnu = dcmplx(Rnu(i,j,:),Inu(i,j,:))
     CB = dcmplx(RB(i,j,:),IB(i,j,:))
#if 0     
     call cderivs_sw_x(ex(3),R,CJ(i,j,:),CJx)
     call cderivs_sw_x(ex(3),R,DCJ(i,j,:),DCJx)
     call cderivs_sw_x(ex(3),R,Ck,Ckx)
     call cderivs_sw_x(ex(3),R,Cnu,Cnux)
     call cderivs_sw_x(ex(3),R,CB,CBx)
#else
     call cderivs_x(ex(3),R,CJ(i,j,:),CJx)
     call cderivs_x(ex(3),R,DCJ(i,j,:),DCJx)
     call cderivs_x(ex(3),R,Ck,Ckx)
     call cderivs_x(ex(3),R,Cnu,Cnux)
     call cderivs_x(ex(3),R,CB,CBx)
#endif

     do k = 1,ex(3)
        RHS(k) = Q_rhs(R(k),CJ(i,j,k),CJx(k),DCJx(k),KK(i,j,k),Ck(k),Ckx(k),Cnux(k),KKx(i,j,k),CBx(k),Cnu(k),DCJ(i,j,k),CB(k),ZEO)
        gunc(k) = -2/R(k)/(1-R(k))
     enddo

     k = 1
     CQ(k+1) = CQ(k) + (RHS(k+1)+RHS(k)+CQ(k)*gunc(k))*dR/2
     CQ(k+1) = CQ(k+1)/(1-0.5*dR*gunc(k+1))

     k = 2
     CQ(k+1) = CQ(k) + (F5o12*RHS(k+1) + F2o3*(RHS(k)+CQ(k)*gunc(k)) - F1o12*(RHS(k-1)+CQ(k-1)*gunc(k-1)))*dR
     CQ(k+1) = CQ(k+1)/(1-F5o12*dR*gunc(k+1))

     do k=3,ex(3)-2
        CQ(k+1) = CQ(k) + (F3o8*RHS(k+1) + F19o24*(RHS(k)+CQ(k)*gunc(k)) - F5o24*(RHS(k-1)+CQ(k-1)*gunc(k-1)) &
                  + F1o24*(RHS(k-2)+CQ(k-2)*gunc(k-2)))*dR
        CQ(k+1) = CQ(k+1)/(1-F3o8*dR*gunc(k+1))
     enddo

     k = ex(3)
     CQ(k) = -2*CB(k)

     RQ(i,j,:) = dreal(CQ)
     IQ(i,j,:) = dimag(CQ)

  enddo
  enddo

  gont = 0

  return

end function NullEvol_Q
!--------------------------------------------------------------------
! this R is indeed x
function NullEvol_U(ex,crho,sigma,R,RJ,IJ,RQ,IQ,KK,HKK,beta,RU,IU, &
                    Rmin) result(gont)   
  implicit none
  integer,intent(in ):: ex(1:3)
  real*8,intent(in),dimension(ex(1))::crho
  real*8,intent(in),dimension(ex(2))::sigma
  real*8,intent(in),dimension(ex(3))::R
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: RJ,IJ,RQ,IQ,beta,KK,HKK
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: RU,IU
  real*8,intent(in) :: Rmin
!  gont = 0: success; gont = 1: something wrong
  integer::gont
  real*8 :: dR

  double complex,dimension(ex(3)) :: CU0,RHS
  integer :: i,j,k
  double complex :: U_rhs
  double complex,dimension(ex(3)) :: CJ,CQ

  real*8,parameter :: F5o12=2.d0/1.2d1,F2o3=2.d0/3.d0,F1o12=1.d0/1.2d1
  real*8,parameter :: F3o8=3.d0/8.d0,F19o24=1.9d1/2.4d1,F5o24=5.d0/2.4d1,F1o24=1.d0/2.4d1

!!! sanity check
  dR = sum(RJ)+sum(IJ)+sum(RQ)+sum(IQ)+sum(beta)+sum(RU)+sum(IU)+sum(KK)
  if(dR.ne.dR) then
     if(sum(RJ).ne.sum(RJ))write(*,*)"NullEvol_U: find NaN in RJ"
     if(sum(IJ).ne.sum(IJ))write(*,*)"NullEvol_U: find NaN in IJ"
     if(sum(RQ).ne.sum(RQ))write(*,*)"NullEvol_U: find NaN in RQ"
     if(sum(IQ).ne.sum(IQ))write(*,*)"NullEvol_U: find NaN in IQ"
     if(sum(beta).ne.sum(beta))write(*,*)"NullEvol_U: find NaN in beta"
     if(sum(RU).ne.sum(RU))write(*,*)"NullEvol_U: find NaN in RU"
     if(sum(IU).ne.sum(IU))write(*,*)"NullEvol_U: find NaN in IU"
     if(sum(KK).ne.sum(KK))write(*,*)"NullEvol_U: find NaN in KK"
     gont = 1
     return
  endif

  dR = R(2) - R(1)
  
  do j=1,ex(2)
  do i=1,ex(1)
     CU0(1) = dcmplx(RU(i,j,1),IU(i,j,1))
     CJ = dcmplx(RJ(i,j,:),IJ(i,j,:))
     CQ = dcmplx(RQ(i,j,:),IQ(i,j,:))

     do k = 1,ex(3)
        RHS(k) = U_rhs(R(k),Rmin,beta(i,j,k),KK(i,j,k),CQ(k),CJ(k))
     enddo

     k = 1
     CU0(k+1) = CU0(k) + (RHS(k+1)+RHS(k))*dR/2

     k = 2
     CU0(k+1) = CU0(k) + (F5o12*RHS(k+1) + F2o3*RHS(k) - F1o12*RHS(k-1))*dR

     do k=3,ex(3)-1
        CU0(k+1) = CU0(k) + (F3o8*RHS(k+1) + F19o24*RHS(k) - F5o24*RHS(k-1) &
                  + F1o24*RHS(k-2))*dR
     enddo

     RU(i,j,:) = dreal(CU0)
     IU(i,j,:) = dimag(CU0)

  enddo
  enddo

  gont = 0
  return

end function NullEvol_U
!----------------------------------------------------------------------------------------
! this R is indeed x
function NullEvol_W(ex,crho,sigma,R,RJ,IJ,RB,IB,Rnu,Inu,Rk,Ik, &
                    RU,IU,RQ,IQ,W,beta,KK,HKK,Rmin, &
                    qlR1,qlR2,qlI1,qlI2,quR1,quR2,quI1,quI2,gR,gI)    result(gont)   
  implicit none
  integer,intent(in ):: ex(1:3)
  real*8,intent(in),dimension(ex(1))::crho
  real*8,intent(in),dimension(ex(2))::sigma
  real*8,intent(in),dimension(ex(3))::R
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: W
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: RJ,IJ,RB,IB
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: Rnu,Inu,Rk,Ik
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: RU,IU,RQ,IQ,beta,KK,HKK
  real*8,intent(in ) :: Rmin
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: qlR1,qlR2,qlI1,qlI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: quR1,quR2,quI1,quI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: gR,gI
!  gont = 0: success; gont = 1: something wrong
  integer::gont
  real*8 :: dR

  double complex, dimension(ex(1),ex(2),ex(3)) :: CU,DCU,bDCU
  double complex, dimension(ex(1),ex(2),ex(3)) :: CB,DCB,bDCB,CJ,DCJ,Cnu,bDCnu,Ck,bDCk
  double complex, dimension(ex(3)) :: bDCUx,CQ
  integer :: i,j,k
  real*8, dimension(ex(3)) :: rhs,gunc
  real*8 :: zeo,W_rhs

  real*8,parameter :: F5o12=2.d0/1.2d1,F2o3=2.d0/3.d0,F1o12=1.d0/1.2d1
  real*8,parameter :: F3o8=3.d0/8.d0,F19o24=1.9d1/2.4d1,F5o24=5.d0/2.4d1,F1o24=1.d0/2.4d1

!!! sanity check
  dR = sum(RJ)+sum(IJ)+sum(beta)+sum(RB)+sum(IB)+sum(Rnu)+sum(Inu) &
      +sum(Rk)+sum(Ik)+sum(W)+sum(RU)+sum(IU)+sum(RQ)+sum(IQ)&
      +sum(KK)
  if(dR.ne.dR) then
     if(sum(RJ).ne.sum(RJ))write(*,*)"NullEvol_W: find NaN in RJ"
     if(sum(IJ).ne.sum(IJ))write(*,*)"NullEvol_W: find NaN in IJ"
     if(sum(beta).ne.sum(beta))write(*,*)"NullEvol_W: find NaN in beta"
     if(sum(RB).ne.sum(RB))write(*,*)"NullEvol_W: find NaN in RB"
     if(sum(IB).ne.sum(IB))write(*,*)"NullEvol_W: find NaN in IB"
     if(sum(Rnu).ne.sum(Rnu))write(*,*)"NullEvol_W: find NaN in Rnu"
     if(sum(Inu).ne.sum(Inu))write(*,*)"NullEvol_W: find NaN in Inu"
     if(sum(Rk).ne.sum(Rk))write(*,*)"NullEvol_W: find NaN in Rk"
     if(sum(Ik).ne.sum(Ik))write(*,*)"NullEvol_W: find NaN in Ik"
     if(sum(RU).ne.sum(RU))write(*,*)"NullEvol_W: find NaN in RU"
     if(sum(IU).ne.sum(IU))write(*,*)"NullEvol_W: find NaN in IU"
     if(sum(RQ).ne.sum(RQ))write(*,*)"NullEvol_W: find NaN in RQ"
     if(sum(IQ).ne.sum(IQ))write(*,*)"NullEvol_W: find NaN in IQ"
     if(sum(W).ne.sum(W))write(*,*)"NullEvol_W: find NaN in W"
     if(sum(KK).ne.sum(KK))write(*,*)"NullEvol_W: find NaN in KK"
     gont = 1
     return
  endif

  dR = R(2) - R(1)
  zeo = 0.d0
  
  CB = dcmplx(RB,IB)
  CU = dcmplx(RU,IU)
  Ck = dcmplx(Rk,Ik)
  Cnu = dcmplx(Rnu,Inu)
  CJ = dcmplx(RJ,IJ)
  do k=1,ex(3)
     call derivs_eth(ex(1:2),crho,sigma,CJ(:,:,k),DCJ(:,:,k),2,1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CB(:,:,k),DCB(:,:,k),1,1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CB(:,:,k),bDCB(:,:,k),1,-1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CU(:,:,k),bDCU(:,:,k),1,-1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,Ck(:,:,k),bDCk(:,:,k),1,-1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,Cnu(:,:,k),bDCnu(:,:,k),1,-1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
  enddo

  do j=1,ex(2)
  do i=1,ex(1)
#if 0
     call cderivs_sw_x(ex(3),R,bDCU,bDCUx)
#else
     call cderivs_x(ex(3),R,bDCU,bDCUx)
#endif

     CQ = dcmplx(RQ(i,j,:),IQ(i,j,:))

     do k = 1,ex(3)
        rhs(k) = W_rhs(R(k),Rmin,beta(i,j,k),KK(i,j,k),DCB(i,j,k),CB(i,j,k),CJ(i,j,k),Cnu(i,j,k),Ck(i,j,k),zeo, &
                    CQ(k),bDCk(i,j,k),bDCnu(i,j,k),bDCB(i,j,k),bDCU(i,j,k),bDCUx(k),DCJ(i,j,k))
        gunc(k) = -2/R(k)/(1-R(k))
     enddo

     k = 1
     W(i,j,k+1) = W(i,j,k) + (rhs(k+1)+rhs(k)+W(i,j,k)*gunc(k))*dR/2
     W(i,j,k+1) = W(i,j,k+1)/(1-0.5*dR*gunc(k+1))

     k = 2
     W(i,j,k+1) = W(i,j,k) + (F5o12*rhs(k+1) + F2o3*(rhs(k)+W(i,j,k)*gunc(k)) - F1o12*(rhs(k-1)+W(i,j,k-1)*gunc(k-1)))*dR
     W(i,j,k+1) = W(i,j,k+1)/(1-F5o12*dR*gunc(k+1))

     do k=3,ex(3)-2
        W(i,j,k+1) = W(i,j,k) + (F3o8*rhs(k+1) + F19o24*(rhs(k)+W(i,j,k)*gunc(k)) - F5o24*(rhs(k-1)+W(i,j,k-1)*gunc(k-1)) &
                  + F1o24*(rhs(k-2)+W(i,j,k-2)*gunc(k-2)))*dR
        W(i,j,k+1) = W(i,j,k+1)/(1-F3o8*dR*gunc(k+1))
     enddo

     k = ex(3)
     W(i,j,k) = dreal(bDCU(i,j,k))

  enddo
  enddo

  gont = 0
  return

end function NullEvol_W
!--------------------------------------------------------------------
! this R is indeed x
function NullEvol_Theta(ex,crho,sigma,R,RJ,IJ,RU,IU,beta,RB,IB, &
                        Rnu,Inu,Rk,Ik,RTheta,ITheta,W,KK,HKK,KKx,HKKx,Rmin,       &
                        qlR1,qlR2,qlI1,qlI2,quR1,quR2,quI1,quI2,gR,gI) result(gont)   
  implicit none
  integer,intent(in ):: ex(1:3)
  real*8,intent(in),dimension(ex(1))::crho
  real*8,intent(in),dimension(ex(2))::sigma
  real*8,intent(in),dimension(ex(3))::R
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: beta,W,KK,KKx,HKK,HKKx
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: RJ,IJ,RU,IU,RB,IB,Rnu,Inu,Rk,Ik
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: RTheta,ITheta
  real*8,intent(in) :: Rmin
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: qlR1,qlR2,qlI1,qlI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: quR1,quR2,quI1,quI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: gR,gI
!  gont = 0: success; gont = 1: something wrong
  integer::gont

  double complex,dimension(ex(1),ex(2),ex(3)) :: CU,DCU,bDCU,CB,DCB,bDCB,CJ,DCJ
  double complex,dimension(ex(3)) :: CTheta0,RHS
  integer :: i,j,k,RK4
  double complex,dimension(ex(3)) :: Cnu,Ck,CUx,DCUx,bDCUx
  double complex,dimension(ex(3)) :: Cnux,CJx,CJxx,DCJx
  real*8,dimension(ex(3)) :: betax,Wx,gunc
  double complex :: Theta_rhs,ZEO
  real*8 :: dR

  real*8,parameter :: F5o12=2.d0/1.2d1,F2o3=2.d0/3.d0,F1o12=1.d0/1.2d1
  real*8,parameter :: F3o8=3.d0/8.d0,F19o24=1.9d1/2.4d1,F5o24=5.d0/2.4d1,F1o24=1.d0/2.4d1

!!! sanity check
  dR = sum(RJ)+sum(IJ)+sum(RU)+sum(IU)+sum(beta)+sum(RB)+sum(IB) + &
       sum(Rnu)+sum(Inu)+sum(Rk)+sum(Ik)+sum(RTheta)+sum(ITheta) + &
       sum(KK)+sum(KKx)+sum(W)
  if(dR.ne.dR) then
     if(sum(RJ).ne.sum(RJ))write(*,*)"NullEvol_Theta: find NaN in RJ"
     if(sum(IJ).ne.sum(IJ))write(*,*)"NullEvol_Theta: find NaN in IJ"
     if(sum(RU).ne.sum(RU))write(*,*)"NullEvol_Theta: find NaN in RU"
     if(sum(IU).ne.sum(IU))write(*,*)"NullEvol_Theta: find NaN in IU"
     if(sum(beta).ne.sum(beta))write(*,*)"NullEvol_Theta: find NaN in beta"
     if(sum(RB).ne.sum(RB))write(*,*)"NullEvol_Theta: find NaN in RB"
     if(sum(IB).ne.sum(IB))write(*,*)"NullEvol_Theta: find NaN in IB"
     if(sum(Rnu).ne.sum(Rnu))write(*,*)"NullEvol_Theta: find NaN in Rnu"
     if(sum(Inu).ne.sum(Inu))write(*,*)"NullEvol_Theta: find NaN in Inu"
     if(sum(Rk).ne.sum(Rk))write(*,*)"NullEvol_Theta: find NaN in Rk"
     if(sum(Ik).ne.sum(Ik))write(*,*)"NullEvol_Theta: find NaN in Ik"
     if(sum(RTheta).ne.sum(RTheta))write(*,*)"NullEvol_Theta: find NaN in RTheta"
     if(sum(ITheta).ne.sum(ITheta))write(*,*)"NullEvol_Theta: find NaN in ITheta"
     if(sum(KK).ne.sum(KK))write(*,*)"NullEvol_Theta: find NaN in KK"
     if(sum(KKx).ne.sum(KKx))write(*,*)"NullEvol_Theta: find NaN in KKx"
     if(sum(W).ne.sum(W))write(*,*)"NullEvol_Theta: find NaN in W"
     gont = 1
     return
  endif

  dR = R(2) - R(1)
  ZEO = dcmplx(0.d0,0.d0)
  
  CU = dcmplx(RU,IU)
  CB = dcmplx(RB,IB)
  CJ = dcmplx(RJ,IJ)

  do k=1,ex(3)
     call derivs_eth(ex(1:2),crho,sigma,CU(:,:,k),DCU(:,:,k),1,1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CU(:,:,k),bDCU(:,:,k),1,-1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CB(:,:,k),DCB(:,:,k),1,1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CB(:,:,k),bDCB(:,:,k),1,-1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CJ(:,:,k),DCJ(:,:,k),2,1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
  enddo

  do j=1,ex(2)
  do i=1,ex(1)
     CTheta0(1) = dcmplx(RTheta(i,j,1),ITheta(i,j,1))
     Cnu = dcmplx(Rnu(i,j,:),Inu(i,j,:))
     Ck = dcmplx(Rk(i,j,:),Ik(i,j,:))
#if 0     
     call cderivs_sw_x(ex(3),R,Cnu,Cnux)
     call rderivs_sw_x(ex(3),R,beta(i,j,:),betax)
     call rderivs_sw_x(ex(3),R,W,Wx)
     call cderivs_sw_x(ex(3),R,DCU(i,j,:),DCUx)
     call cderivs_sw_x(ex(3),R,CU(i,j,:),CUx)
     call cderivs_sw_x(ex(3),R,bDCU(i,j,:),bDCUx)
     call cderivs_sw_x(ex(3),R,CJ(i,j,:),CJx)
     call cdderivs_sw_x(ex(3),R,CJ(i,j,:),CJxx)
     call cderivs_sw_x(ex(3),R,DCJ(i,j,:),DCJx)
#else
     call cderivs_x(ex(3),R,Cnu,Cnux)
     call rderivs_x(ex(3),R,beta(i,j,:),betax)
     call rderivs_x(ex(3),R,W,Wx)
     call cderivs_x(ex(3),R,DCU(i,j,:),DCUx)
     call cderivs_x(ex(3),R,CU(i,j,:),CUx)
     call cderivs_x(ex(3),R,bDCU(i,j,:),bDCUx)
     call cderivs_x(ex(3),R,CJ(i,j,:),CJx)
     call cdderivs_x(ex(3),R,CJ(i,j,:),CJxx)
     call cderivs_x(ex(3),R,DCJ(i,j,:),DCJx)
#endif
     do k = 1,ex(3)
        rhs(k) = Theta_rhs(R(k),Rmin,beta(i,j,k),betax(k),KK(i,j,k),KKx(i,j,k),CU(i,j,k),CUx(k),DCUx(k),bDCU(i,j,k),bDCUx(k), &
                        DCU(i,j,k),CB(i,j,k),DCB(i,j,k),W(i,j,k),Wx(k),CJ(i,j,k),DCJ(i,j,k),    &
                        CJx(k),CJxx(k),DCJx(k),bDCB(i,j,k),Cnu(k),Cnux(k),Ck(k),ZEO)
        gunc(k) = -1/R(k)/(1-R(k))
     enddo

     k = 1
     CTheta0(k+1) = CTheta0(k) + (RHS(k+1)+RHS(k)+CTheta0(k)*gunc(k))*dR/2
     CTheta0(k+1) = CTheta0(k+1)/(1-0.5*dR*gunc(k+1))

     k = 2
     CTheta0(k+1) = CTheta0(k) + (F5o12*RHS(k+1) + F2o3*(RHS(k)+CTheta0(k)*gunc(k)) - F1o12*(RHS(k-1)+CTheta0(k-1)*gunc(k-1)))*dR
     CTheta0(k+1) = CTheta0(k+1)/(1-F5o12*dR*gunc(k+1))

     do k=3,ex(3)-2
        CTheta0(k+1) = CTheta0(k) + (F3o8*RHS(k+1) + F19o24*(RHS(k)+CTheta0(k)*gunc(k)) - F5o24*(RHS(k-1)+CTheta0(k-1)*gunc(k-1)) &
                  + F1o24*(RHS(k-2)+CTheta0(k-2)*gunc(k-2)))*dR
        CTheta0(k+1) = CTheta0(k+1)/(1-F3o8*dR*gunc(k+1))
     enddo

     k = ex(3)
     CTheta0(k) = -KK(i,j,k)*DCU(i,j,k)-(CU(i,j,k)*Cnu(k)+dconjg(CU(i,j,k))*DCJ(i,j,k))/2 &
                  +CJ(i,j,k)*(bDCU(i,j,k)-dconjg(bDCU(i,j,k)))/2 - W(i,j,k)*CJ(i,j,k)/2

     RTheta(i,j,:) = dreal(CTheta0)
     ITheta(i,j,:) = dimag(CTheta0)
  enddo
  enddo

  gont = 0
  return

end function NullEvol_Theta

#else
#error "not recognized RKorAM"
#endif

!=====================================================================================================================================
! basic tool routines
  subroutine rswap(r1,r2)

  implicit none

!~~~~~~% Input parameters:

  real*8,intent(inout) :: r1,r2
  
  real*8 :: r

  r = r1
  r1= r2
  r2= r

  return

  end subroutine rswap
!----
  subroutine cswap(r1,r2)

  implicit none

!~~~~~~% Input parameters:

  double complex,intent(inout) :: r1,r2
  
  double complex :: r

  r = r1
  r1= r2
  r2= r

  return

  end subroutine cswap

! center type finite difference
!====================================================================================
!----
  subroutine rderivs_x(lx,X,f,fx)

  implicit none

!~~~~~~% Input parameters:
  integer,intent(in) :: lx
  real*8,intent(in),dimension(lx) :: X
  real*8,intent(in),dimension(lx) :: f
  real*8,intent(out),dimension(lx) :: fx

  real*8 :: dX

  dX = X(2)-X(1)

#ifdef OLD
  fx(1:lx-1) = (f(2:lx)-f(1:lx-1))/dX
  fx(lx)     = (f(lx)-f(lx-1))/dX
#else

#if (ghost_width == 2)
  fx(2:lx-1) = (f(3:lx)-f(1:lx-2))/2.d0/dX
  fx(1)      = (f(2)-f(1))/dX
  fx(lx)     = (f(lx)-f(lx-1))/dX
#elif (ghost_width == 3)
  fx(3:lx-2) = (f(1:lx-4)-8.d0*f(2:lx-3)+8.d0*f(4:lx-1)-f(5:lx))/1.2d1/dX
  fx(2)      = (f(3)-f(1))/2.d0/dX
  fx(lx-1)   = (f(lx)-f(lx-2))/2.d0/dX
  fx(1)      = (f(2)-f(1))/dX
  fx(lx)     = (f(lx)-f(lx-1))/dX
!  fx(1)    =-(2.5d1*f(1)-4.8d1*f(2)+3.6d1*f(3)-1.6d1*f(4)+3.d0*f(5))/1.2d1/dX
!  fx(2)    =-(3.d0*f(1)+1.d1*f(2)-1.8d1*f(3)+6.d0*f(4)-f(5))/1.2d1/dX
#elif (ghost_width == 4)
  fx(4:lx-3) = (-f(1:lx-6)+9.d0*f(2:lx-5)-4.5d1*f(3:lx-4)+4.5d1*f(5:lx-2)-9.d0*f(6:lx-1)+f(7:lx))/6.d1/dX
  fx(3)      = (f(1)-8.d0*f(2)+8.d0*f(4)-f(5))/1.2d1/dX
  fx(lx-2)   = (f(lx-4)-8.d0*f(lx-3)+8.d0*f(lx-1)-f(lx))/1.2d1/dX
  fx(2)      = (f(3)-f(1))/2.d0/dX
  fx(lx-1)   = (f(lx)-f(lx-2))/2.d0/dX
  fx(1)      = (f(2)-f(1))/dX
  fx(lx)     = (f(lx)-f(lx-1))/dX
#elif (ghost_width == 5)
  fx(5:lx-4) = (3.d0*f(1:lx-8)-3.2d1*f(2:lx-7)+1.68d2*f(3:lx-6)-6.72d2*f(4:lx-5)+ &
                6.72d2*f(6:lx-3)-1.68d2*f(7:lx-2)+3.2d1*f(8:lx-1)-3.d0*f(9:lx))/8.4d2/dX
  fx(4)      = (-f(1)+9.d0*f(2)-4.5d1*f(3)+4.5d1*f(5)-9.d0*f(6)+f(7))/6.d1/dX
  fx(lx-3)   = (-f(lx-6)+9.d0*f(lx-5)-4.5d1*f(lx-4)+4.5d1*f(lx-2)-9.d0*f(lx-1)+f(lx))/6.d1/dX
  fx(3)      = (f(1)-8.d0*f(2)+8.d0*f(4)-f(5))/1.2d1/dX
  fx(lx-2)   = (f(lx-4)-8.d0*f(lx-3)+8.d0*f(lx-1)-f(lx))/1.2d1/dX
  fx(2)      = (f(3)-f(1))/2.d0/dX
  fx(lx-1)   = (f(lx)-f(lx-2))/2.d0/dX
  fx(1)      = (f(2)-f(1))/dX
  fx(lx)     = (f(lx)-f(lx-1))/dX
#endif  

#endif
  return

  end subroutine rderivs_x
!----
  subroutine rderivs_x_point(lx,X,f,fx,k)

  implicit none

!~~~~~~% Input parameters:
  integer,intent(in) :: lx,k
  real*8,intent(in),dimension(lx) :: X
  real*8,intent(in),dimension(lx) :: f
  real*8,intent(out) :: fx

  real*8 :: dX

  dX = X(2)-X(1)

#ifdef OLD
  if(k .eq. lx)then
    fx = (f(lx)-f(lx-1))/dX
  else
    fx = (f(k+1)-f(k))/dX
  endif 
#else

#if (ghost_width == 2)
  if(k .gt. 1 .and. k .lt. lx) then
    fx = (f(k+1)-f(k-1))/2.d0/dX
  elseif(k.eq.1) then
    fx = (f(2)-f(1))/dX
  elseif(k.eq.lx) then
    fx = (f(lx)-f(lx-1))/dX
  endif
#elif (ghost_width == 3)
  if(k .gt. 2 .and. k .lt. lx-1) then
    fx = (f(k-2)-8.d0*f(k-1)+8.d0*f(k+1)-f(k+2))/1.2d1/dX
  elseif(k.eq.1) then
    fx = (f(2)-f(1))/dX
  elseif(k.eq.lx) then
    fx = (f(lx)-f(lx-1))/dX
  elseif(k.eq.2) then
    fx = (f(3)-f(1))/2.d0/dX
  elseif(k.eq.lx-1) then
    fx = (f(lx)-f(lx-2))/2.d0/dX
  endif
#elif (ghost_width == 4)
  if(k .gt. 3 .and. k .lt. lx-2) then
    fx = (-f(k-3)+9.d0*f(k-2)-4.5d1*f(k-1)+4.5d1*f(k+1)-9.d0*f(k+2)+f(k+3))/6.d1/dX
  elseif(k.eq.1) then
    fx = (f(2)-f(1))/dX
  elseif(k.eq.lx) then
    fx = (f(lx)-f(lx-1))/dX
  elseif(k.eq.2) then
    fx = (f(3)-f(1))/2.d0/dX
  elseif(k.eq.lx-1) then
    fx = (f(lx)-f(lx-2))/2.d0/dX
  elseif(k.eq.3) then
    fx = (f(1)-8.d0*f(2)+8.d0*f(4)-f(5))/1.2d1/dX
  elseif(k.eq.lx-2) then
    fx = (f(lx-4)-8.d0*f(lx-3)+8.d0*f(lx-1)-f(lx))/1.2d1/dX
  endif
#elif (ghost_width == 5)
  if(k .gt. 4 .and. k .lt. lx-3) then
    fx = (3.d0*f(k-4)-3.2d1*f(k-3)+1.68d2*f(k-2)-6.72d2*f(k-1)+ &
                6.72d2*f(k+1)-1.68d2*f(k+2)+3.2d1*f(k+3)-3.d0*f(k+4))/8.4d2/dX
  elseif(k.eq.1) then
    fx = (f(2)-f(1))/dX
  elseif(k.eq.lx) then
    fx = (f(lx)-f(lx-1))/dX
  elseif(k.eq.2) then
    fx = (f(3)-f(1))/2.d0/dX
  elseif(k.eq.lx-1) then
    fx = (f(lx)-f(lx-2))/2.d0/dX
  elseif(k.eq.3) then
    fx = (f(1)-8.d0*f(2)+8.d0*f(4)-f(5))/1.2d1/dX
  elseif(k.eq.lx-2) then
    fx = (f(lx-4)-8.d0*f(lx-3)+8.d0*f(lx-1)-f(lx))/1.2d1/dX
  elseif(k.eq.4) then
    fx = (-f(1)+9.d0*f(2)-4.5d1*f(3)+4.5d1*f(5)-9.d0*f(6)+f(7))/6.d1/dX
  elseif(k.eq.lx-3) then
    fx = (-f(lx-6)+9.d0*f(lx-5)-4.5d1*f(lx-4)+4.5d1*f(lx-2)-9.d0*f(lx-1)+f(lx))/6.d1/dX
  endif
#endif  

#endif
  return

  end subroutine rderivs_x_point
!----
  subroutine cderivs_x(lx,X,f,fx)

  implicit none

!~~~~~~% Input parameters:
  integer,intent(in) :: lx
  real*8,intent(in),dimension(lx) :: X
  double complex,intent(in),dimension(lx) :: f
  double complex,intent(out),dimension(lx) :: fx

  real*8 :: dX

  dX = X(2)-X(1)

#ifdef OLD
  fx(1:lx-1) = (f(2:lx)-f(1:lx-1))/dX
  fx(lx)     = (f(lx)-f(lx-1))/dX
#else

#if (ghost_width == 2)
  fx(2:lx-1) = (f(3:lx)-f(1:lx-2))/2.d0/dX
  fx(1)      = (f(2)-f(1))/dX
  fx(lx)     = (f(lx)-f(lx-1))/dX
#elif (ghost_width == 3)
  fx(3:lx-2) = (f(1:lx-4)-8.d0*f(2:lx-3)+8.d0*f(4:lx-1)-f(5:lx))/1.2d1/dX
  fx(2)      = (f(3)-f(1))/2.d0/dX
  fx(lx-1)   = (f(lx)-f(lx-2))/2.d0/dX
  fx(1)      = (f(2)-f(1))/dX
  fx(lx)     = (f(lx)-f(lx-1))/dX
!  fx(1)    =-(2.5d1*f(1)-4.8d1*f(2)+3.6d1*f(3)-1.6d1*f(4)+3.d0*f(5))/1.2d1/dX
!  fx(2)    =-(3.d0*f(1)+1.d1*f(2)-1.8d1*f(3)+6.d0*f(4)-f(5))/1.2d1/dX
#elif (ghost_width == 4)
  fx(4:lx-3) = (-f(1:lx-6)+9.d0*f(2:lx-5)-4.5d1*f(3:lx-4)+4.5d1*f(5:lx-2)-9.d0*f(6:lx-1)+f(7:lx))/6.d1/dX
  fx(3)      = (f(1)-8.d0*f(2)+8.d0*f(4)-f(5))/1.2d1/dX
  fx(lx-2)   = (f(lx-4)-8.d0*f(lx-3)+8.d0*f(lx-1)-f(lx))/1.2d1/dX
  fx(2)      = (f(3)-f(1))/2.d0/dX
  fx(lx-1)   = (f(lx)-f(lx-2))/2.d0/dX
  fx(1)      = (f(2)-f(1))/dX
  fx(lx)     = (f(lx)-f(lx-1))/dX
#elif (ghost_width == 5)
  fx(5:lx-4) = (3.d0*f(1:lx-8)-3.2d1*f(2:lx-7)+1.68d2*f(3:lx-6)-6.72d2*f(4:lx-5)+ &
                6.72d2*f(6:lx-3)-1.68d2*f(7:lx-2)+3.2d1*f(8:lx-1)-3.d0*f(9:lx))/8.4d2/dX
  fx(4)      = (-f(1)+9.d0*f(2)-4.5d1*f(3)+4.5d1*f(5)-9.d0*f(6)+f(7))/6.d1/dX
  fx(lx-3)   = (-f(lx-6)+9.d0*f(lx-5)-4.5d1*f(lx-4)+4.5d1*f(lx-2)-9.d0*f(lx-1)+f(lx))/6.d1/dX
  fx(3)      = (f(1)-8.d0*f(2)+8.d0*f(4)-f(5))/1.2d1/dX
  fx(lx-2)   = (f(lx-4)-8.d0*f(lx-3)+8.d0*f(lx-1)-f(lx))/1.2d1/dX
  fx(2)      = (f(3)-f(1))/2.d0/dX
  fx(lx-1)   = (f(lx)-f(lx-2))/2.d0/dX
  fx(1)      = (f(2)-f(1))/dX
  fx(lx)     = (f(lx)-f(lx-1))/dX
#endif  

#endif

  return

  end subroutine cderivs_x
!----
  subroutine cdderivs_x(lx,X,f,fxx)

  implicit none

!~~~~~~% Input parameters:
  integer,intent(in) :: lx
  real*8,intent(in),dimension(lx) :: X
  double complex,intent(in),dimension(lx) :: f
  double complex,intent(out),dimension(lx) :: fxx

  real*8 :: dX

  dX = X(2)-X(1)
  dX = dX*dX

#ifdef OLD
  fxx(1:lx-2) = (f(3:lx)-2.0*f(2:lx-1)+f(1:lx-2))/dX
  fxx(lx-1)   = (f(lx)-2.0*f(lx-1)+f(lx-2))/dX
  fxx(lx  )   = (f(lx)-2.0*f(lx-1)+f(lx-2))/dX
#else

#if (ghost_width == 2)
  fxx(2:lx-1) = (f(3:lx)-2.d0*f(2:lx-1)+f(1:lx-2))/dX
  fxx(1)      = (f(3)-2.d0*f(2)+f(1))/dX
  fxx(lx)     = (f(lx)-2.d0*f(lx-1)+f(lx-2))/dX
#elif (ghost_width == 3)
  fxx(3:lx-2) = (-f(1:lx-4)+1.6d1*f(2:lx-3)-3.d1*f(3:lx-2)+1.6d1*f(4:lx-1)-f(5:lx))/1.2d1/dX
  fxx(2)      = (f(3)-2.d0*f(2)+f(1))/dX
  fxx(lx-1)   = (f(lx)-2.d0*f(lx-1)+f(lx-2))/dX
  fxx(1)      = (f(3)-2.d0*f(2)+f(1))/dX
  fxx(lx)     = (f(lx)-2.d0*f(lx-1)+f(lx-2))/dX
#elif (ghost_width == 4)
  fxx(4:lx-3) = (2.d0*f(1:lx-6)-2.7d1*f(2:lx-5)+2.7d2*f(3:lx-4)-4.9d2*f(4:lx-3) &
                +2.7d2*f(5:lx-2)-2.7d1*f(6:lx-1)+2.d0*f(7:lx))/1.8d2/dX
  fxx(3)      = (-f(1)+1.6d1*f(2)-3.d1*f(3)+1.6d1*f(4)-f(5))/1.2d1/dX
  fxx(lx-2)   = (-f(lx-4)+1.6d1*f(lx-3)-3.d1*f(lx-2)+1.6d1*f(lx-1)-f(lx))/1.2d1/dX
  fxx(2)      = (f(3)-2.d0*f(2)+f(1))/dX
  fxx(lx-1)   = (f(lx)-2.d0*f(lx-1)+f(lx-2))/dX
  fxx(1)      = (f(3)-2.d0*f(2)+f(1))/dX
  fxx(lx)     = (f(lx)-2.d0*f(lx-1)+f(lx-2))/dX
#elif (ghost_width == 5)
  fxx(5:lx-4) = (-9.d0*f(1:lx-8)+1.28d2*f(2:lx-7)-1.008d3*f(3:lx-6)+8.064d3*f(4:lx-5)-1.435d4*f(5:lx-4) &
                +8.064d3*f(6:lx-3)-1.008d3*f(7:lx-2)+1.28d2*f(8:lx-1)-9.d0*f(9:lx))/5.04d3/dX
  fxx(4) = (2.d0*f(1)-2.7d1*f(2)+2.7d2*f(3)-4.9d2*f(4) &
           +2.7d2*f(5)-2.7d1*f(6)+2.d0*f(7))/1.8d2/dX
  fxx(lx-3) = (2.d0*f(lx-6)-2.7d1*f(lx-5)+2.7d2*f(lx-4)-4.9d2*f(lx-3) &
              +2.7d2*f(lx-2)-2.7d1*f(lx-1)+2.d0*f(lx))/1.8d2/dX
  fxx(3)    = (-f(1)+1.6d1*f(2)-3.d1*f(3)+1.6d1*f(4)-f(5))/1.2d1/dX
  fxx(lx-2) = (-f(lx-4)+1.6d1*f(lx-3)-3.d1*f(lx-2)+1.6d1*f(lx-1)-f(lx))/1.2d1/dX
  fxx(2)    = (f(3)-2.d0*f(2)+f(1))/dX
  fxx(lx-1) = (f(lx)-2.d0*f(lx-1)+f(lx-2))/dX
  fxx(1)    = (f(3)-2.d0*f(2)+f(1))/dX
  fxx(lx)   = (f(lx)-2.d0*f(lx-1)+f(lx-2))/dX
#endif  

#endif

  return

  end subroutine cdderivs_x
!----
  subroutine rdderivs_x(lx,X,f,fxx)

  implicit none

!~~~~~~% Input parameters:
  integer,intent(in) :: lx
  real*8,intent(in),dimension(lx) :: X
  real*8,intent(in),dimension(lx) :: f
  real*8,intent(out),dimension(lx) :: fxx

  real*8 :: dX

  dX = X(2)-X(1)
  dX = dX*dX

#ifdef OLD
  fxx(1:lx-2) = (f(3:lx)-2.0*f(2:lx-1)+f(1:lx-2))/dX
  fxx(lx-1)   = (f(lx)-2.0*f(lx-1)+f(lx-2))/dX
  fxx(lx  )   = (f(lx)-2.0*f(lx-1)+f(lx-2))/dX
#else

#if (ghost_width == 2)
  fxx(2:lx-1) = (f(3:lx)-2.d0*f(2:lx-1)+f(1:lx-2))/dX
  fxx(1)      = (f(3)-2.d0*f(2)+f(1))/dX
  fxx(lx)     = (f(lx)-2.d0*f(lx-1)+f(lx-2))/dX
#elif (ghost_width == 3)
  fxx(3:lx-2) = (-f(1:lx-4)+1.6d1*f(2:lx-3)-3.d1*f(3:lx-2)+1.6d1*f(4:lx-1)-f(5:lx))/1.2d1/dX
  fxx(2)      = (f(3)-2.d0*f(2)+f(1))/dX
  fxx(lx-1)   = (f(lx)-2.d0*f(lx-1)+f(lx-2))/dX
  fxx(1)      = (f(3)-2.d0*f(2)+f(1))/dX
  fxx(lx)     = (f(lx)-2.d0*f(lx-1)+f(lx-2))/dX
#elif (ghost_width == 4)
  fxx(4:lx-3) = (2.d0*f(1:lx-6)-2.7d1*f(2:lx-5)+2.7d2*f(3:lx-4)-4.9d2*f(4:lx-3) &
                +2.7d2*f(5:lx-2)-2.7d1*f(6:lx-1)+2.d0*f(7:lx))/1.8d2/dX
  fxx(3)    = (-f(1)+1.6d1*f(2)-3.d1*f(3)+1.6d1*f(4)-f(5))/1.2d1/dX
  fxx(lx-2) = (-f(lx-4)+1.6d1*f(lx-3)-3.d1*f(lx-2)+1.6d1*f(lx-1)-f(lx))/1.2d1/dX
  fxx(2)    = (f(3)-2.d0*f(2)+f(1))/dX
  fxx(lx-1) = (f(lx)-2.d0*f(lx-1)+f(lx-2))/dX
  fxx(1)    = (f(3)-2.d0*f(2)+f(1))/dX
  fxx(lx)   = (f(lx)-2.d0*f(lx-1)+f(lx-2))/dX
#elif (ghost_width == 5)
  fxx(5:lx-4) = (-9.d0*f(1:lx-8)+1.28d2*f(2:lx-7)-1.008d3*f(3:lx-6)+8.064d3*f(4:lx-5)-1.435d4*f(5:lx-4) &
                +8.064d3*f(6:lx-3)-1.008d3*f(7:lx-2)+1.28d2*f(8:lx-1)-9.d0*f(9:lx))/5.04d3/dX
  fxx(4) = (2.d0*f(1)-2.7d1*f(2)+2.7d2*f(3)-4.9d2*f(4) &
           +2.7d2*f(5)-2.7d1*f(6)+2.d0*f(7))/1.8d2/dX
  fxx(lx-3) = (2.d0*f(lx-6)-2.7d1*f(lx-5)+2.7d2*f(lx-4)-4.9d2*f(lx-3) &
              +2.7d2*f(lx-2)-2.7d1*f(lx-1)+2.d0*f(lx))/1.8d2/dX
  fxx(3)    = (-f(1)+1.6d1*f(2)-3.d1*f(3)+1.6d1*f(4)-f(5))/1.2d1/dX
  fxx(lx-2) = (-f(lx-4)+1.6d1*f(lx-3)-3.d1*f(lx-2)+1.6d1*f(lx-1)-f(lx))/1.2d1/dX
  fxx(2)    = (f(3)-2.d0*f(2)+f(1))/dX
  fxx(lx-1) = (f(lx)-2.d0*f(lx-1)+f(lx-2))/dX
  fxx(1)    = (f(3)-2.d0*f(2)+f(1))/dX
  fxx(lx)   = (f(lx)-2.d0*f(lx-1)+f(lx-2))/dX
#endif  

#endif

  return

  end subroutine rdderivs_x
!----
  subroutine rdderivs_x_point(lx,X,f,fxx,k)

  implicit none

!~~~~~~% Input parameters:
  integer,intent(in) :: lx,k
  real*8,intent(in),dimension(lx) :: X
  real*8,intent(in),dimension(lx) :: f
  real*8,intent(out) :: fxx

  real*8 :: dX

  dX = X(2)-X(1)
  dX = dX*dX

#ifdef OLD
  if(k.lt.lx-1) then
    fxx = (f(k+2)-2.0*f(k+1)+f(k))/dX
  elseif(k.eq.lx-1) then
    fxx = (f(lx)-2.0*f(lx-1)+f(lx-2))/dX
  elseif(k.eq.lx) then
    fxx = (f(lx)-2.0*f(lx-1)+f(lx-2))/dX
  endif
#else

#if (ghost_width == 2)
  if(k.gt.1 .and. k.lt.lx) then
    fxx = (f(k+1)-2.d0*f(k)+f(k-1))/dX
  elseif(k.eq.1) then
    fxx = (f(3)-2.d0*f(2)+f(1))/dX
  elseif(k.eq.lx) then
    fxx = (f(lx)-2.d0*f(lx-1)+f(lx-2))/dX
  endif
#elif (ghost_width == 3)
  if(k.gt.2 .and. k.lt.lx-1) then
    fxx = (-f(k-2)+1.6d1*f(k-1)-3.d1*f(k)+1.6d1*f(k+1)-f(k+2))/1.2d1/dX
  elseif(k.eq.1) then
    fxx = (f(3)-2.d0*f(2)+f(1))/dX
  elseif(k.eq.lx) then
    fxx = (f(lx)-2.d0*f(lx-1)+f(lx-2))/dX
  elseif(k.eq.2) then
    fxx = (f(3)-2.d0*f(2)+f(1))/dX
  elseif(k.eq.lx-1) then
    fxx = (f(lx)-2.d0*f(lx-1)+f(lx-2))/dX
  endif
#elif (ghost_width == 4)
  if(k.gt.3 .and. k.lt.lx-2)then
    fxx = (2.d0*f(k-3)-2.7d1*f(k-2)+2.7d2*f(k-1)-4.9d2*f(k) &
                +2.7d2*f(k+1)-2.7d1*f(k+2)+2.d0*f(k+3))/1.8d2/dX
  elseif(k.eq.1) then
    fxx = (f(3)-2.d0*f(2)+f(1))/dX
  elseif(k.eq.lx) then
    fxx = (f(lx)-2.d0*f(lx-1)+f(lx-2))/dX
  elseif(k.eq.2) then
    fxx = (f(3)-2.d0*f(2)+f(1))/dX
  elseif(k.eq.lx-1) then
    fxx = (f(lx)-2.d0*f(lx-1)+f(lx-2))/dX
  elseif(k.eq.3) then
    fxx = (-f(1)+1.6d1*f(2)-3.d1*f(3)+1.6d1*f(4)-f(5))/1.2d1/dX
  elseif(k.eq.lx-2) then
    fxx = (-f(lx-4)+1.6d1*f(lx-3)-3.d1*f(lx-2)+1.6d1*f(lx-1)-f(lx))/1.2d1/dX
  endif
#elif (ghost_width == 5)
  if(k.gt.4 .and. k.lt.lx-3) then
    fxx = (-9.d0*f(k-4)+1.28d2*f(k-3)-1.008d3*f(k-2)+8.064d3*f(k-1)-1.435d4*f(k) &
                +8.064d3*f(k+1)-1.008d3*f(k+2)+1.28d2*f(k+3)-9.d0*f(k+4))/5.04d3/dX
  elseif(k.eq.1) then
    fxx = (f(3)-2.d0*f(2)+f(1))/dX
  elseif(k.eq.lx) then
    fxx = (f(lx)-2.d0*f(lx-1)+f(lx-2))/dX
  elseif(k.eq.2) then
    fxx = (f(3)-2.d0*f(2)+f(1))/dX
  elseif(k.eq.lx-1) then
    fxx = (f(lx)-2.d0*f(lx-1)+f(lx-2))/dX
  elseif(k.eq.3) then
    fxx = (-f(1)+1.6d1*f(2)-3.d1*f(3)+1.6d1*f(4)-f(5))/1.2d1/dX
  elseif(k.eq.lx-2) then
    fxx = (-f(lx-4)+1.6d1*f(lx-3)-3.d1*f(lx-2)+1.6d1*f(lx-1)-f(lx))/1.2d1/dX
  elseif(k.eq.4) then
    fxx = (2.d0*f(1)-2.7d1*f(2)+2.7d2*f(3)-4.9d2*f(4) &
              +2.7d2*f(5)-2.7d1*f(6)+2.d0*f(7))/1.8d2/dX
  elseif(k.eq.lx-3) then
    fxx = (2.d0*f(lx-6)-2.7d1*f(lx-5)+2.7d2*f(lx-4)-4.9d2*f(lx-3) &
              +2.7d2*f(lx-2)-2.7d1*f(lx-1)+2.d0*f(lx))/1.8d2/dX
  endif
#endif  

#endif

  return

  end subroutine rdderivs_x_point
!----
  subroutine rdderivs_xy_point(lx,ly,X,Y,f,fxy,i,j)

  implicit none

!~~~~~~% Input parameters:
  integer,intent(in) :: lx,ly,i,j
  real*8,intent(in),dimension(lx) :: X
  real*8,intent(in),dimension(ly) :: Y
  real*8,intent(in),dimension(lx,ly) :: f
  real*8,intent(out) :: fxy

  real*8 :: dX,dY

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
!! we only consider inner points
#if (ghost_width == 2)
  if(i>1 .and. j>1.and.i<lx .and.j<ly)then
   fxy = (f(i-1,j-1)-f(i+1,j-1)-f(i-1,j+1)+f(i+1,j+1))/(4*dX*dY)
  else
   fxy = 0.d0
  endif
#elif (ghost_width == 3)
  if(i>2 .and. j>2.and.i<lx-1 .and.j<ly-1)then
   fxy = (     (f(i-2,j-2)-8*f(i-1,j-2)+8*f(i+1,j-2)-f(i+2,j-2))  &
                       -8 *(f(i-2,j-1)-8*f(i-1,j-1)+8*f(i+1,j-1)-f(i+2,j-1))  &
                       +8 *(f(i-2,j+1)-8*f(i-1,j+1)+8*f(i+1,j+1)-f(i+2,j+1))  &
                       -   (f(i-2,j+2)-8*f(i-1,j+2)+8*f(i+1,j+2)-f(i+2,j+2)))/(144*dX*dY)
  else
   fxy = 0.d0
  endif
#elif (ghost_width == 4)
  if(i>3 .and. j>3.and.i<lx-2 .and.j<ly-2)then
   fxy = (-    (-f(i-3,j-3)+9*f(i-2,j-3)-45*f(i-1,j-3)+45*f(i+1,j-3)-9*f(i+2,j-3)+f(i+3,j-3))  &
                       +9 *(-f(i-3,j-2)+9*f(i-2,j-2)-45*f(i-1,j-2)+45*f(i+1,j-2)-9*f(i+2,j-2)+f(i+3,j-2))  &
                       -45*(-f(i-3,j-1)+9*f(i-2,j-1)-45*f(i-1,j-1)+45*f(i+1,j-1)-9*f(i+2,j-1)+f(i+3,j-1))  &
                       +45*(-f(i-3,j+1)+9*f(i-2,j+1)-45*f(i-1,j+1)+45*f(i+1,j+1)-9*f(i+2,j+1)+f(i+3,j+1))  &
                       -9 *(-f(i-3,j+2)+9*f(i-2,j+2)-45*f(i-1,j+2)+45*f(i+1,j+2)-9*f(i+2,j+2)+f(i+3,j+2))  &
                       +    (-f(i-3,j+3)+9*f(i-2,j+3)-45*f(i-1,j+3)+45*f(i+1,j+3)-9*f(i+2,j+3)+f(i+3,j+3)))/(3.6d3*dX*dY)
  else
   fxy = 0.d0
  endif
#elif (ghost_width == 5)
  if(i>4 .and. j>4.and.i<lx-3 .and.j<ly-3)then
   fxy = ( 3 *( 3*f(i-4,j-4)-32*f(i-3,j-4)+168*f(i-2,j-4)-672*f(i-1,j-4)        &
                              -3*f(i+4,j-4)+32*f(i+3,j-4)-168*f(i+2,j-4)+672*f(i+1,j-4))       &
                       -32 *( 3*f(i-4,j-3)-32*f(i-3,j-3)+168*f(i-2,j-3)-672*f(i-1,j-3)        &
                              -3*f(i+4,j-3)+32*f(i+3,j-3)-168*f(i+2,j-3)+672*f(i+1,j-3))       &
                       +168*( 3*f(i-4,j-2)-32*f(i-3,j-2)+168*f(i-2,j-2)-672*f(i-1,j-2)        &
                              -3*f(i+4,j-2)+32*f(i+3,j-2)-168*f(i+2,j-2)+672*f(i+1,j-2))       &
                       -672*( 3*f(i-4,j-1)-32*f(i-3,j-1)+168*f(i-2,j-1)-672*f(i-1,j-1)        &
                              -3*f(i+4,j-1)+32*f(i+3,j-1)-168*f(i+2,j-1)+672*f(i+1,j-1))       &
                       +672*( 3*f(i-4,j+1)-32*f(i-3,j+1)+168*f(i-2,j+1)-672*f(i-1,j+1)        &
                              -3*f(i+4,j+1)+32*f(i+3,j+1)-168*f(i+2,j+1)+672*f(i+1,j+1))       &
                       -168*( 3*f(i-4,j+2)-32*f(i-3,j+2)+168*f(i-2,j+2)-672*f(i-1,j+2)        &
                              -3*f(i+4,j+2)+32*f(i+3,j+2)-168*f(i+2,j+2)+672*f(i+1,j+2))       &
                       +32 *( 3*f(i-4,j+3)-32*f(i-3,j+3)+168*f(i-2,j+3)-672*f(i-1,j+3)        &
                              -3*f(i+4,j+3)+32*f(i+3,j+3)-168*f(i+2,j+3)+672*f(i+1,j+3))       &
                       -3 *( 3*f(i-4,j+4)-32*f(i-3,j+4)+168*f(i-2,j+4)-672*f(i-1,j+4)        &
                              -3*f(i+4,j+4)+32*f(i+3,j+4)-168*f(i+2,j+4)+672*f(i+1,j+4)) )/(7.056d5*dX*dY)
  else
   fxy = 0.d0
  endif
#endif

  return

  end subroutine rdderivs_xy_point
! side-winded type finite difference (arXiv:1208.3891)  
!====================================================================================

  subroutine rderivs_sw_x(lx,X,f,fx)

  implicit none

!~~~~~~% Input parameters:
  integer,intent(in) :: lx
  real*8,intent(in),dimension(lx) :: X
  real*8,intent(in),dimension(lx) :: f
  real*8,intent(out),dimension(lx) :: fx

  real*8 :: dX
  integer :: i

  dX = X(2)-X(1)

#if (ghost_width == 2)
!#error
  do i=1,lx-3
    fx(i) = (-3.d0*f(i)+4.d0*f(i+1)-f(i+2))/2.d0/dX      ! 向前差分，2价精度
  enddo
  do i=lx-2,lx-1
    fx(i) = (f(i+1)-f(i-1))/2.d0/dX                      ! 中心差分，2价精度
  enddo
  i=lx
    fx(i) = (3.d0*f(i)-4.d0*f(i-1)+f(i-2))/2.d0/dX       ! 向后差分，2价精度
    
#elif (ghost_width == 3)
  do i=1,lx-5
    fx(i) = (-2.5d1/1.2d1*f(i)+4*f(i+1)-3*f(i+2)+4.d0/3.d0*f(i+3)-0.25d0*f(i+4))/dX  ! 向前差分，4价精度
  enddo
  do i=lx-4,lx-2
    fx(i) = (f(i-2)-8.d0*f(i-1)+8.d0*f(i+1)-f(i+2))/1.2d1/dX                         ! 中心差分，4价精度
  enddo
  do i=lx-1,lx
    fx(i) = (2.5d1/1.2d1*f(i)-4*f(i-1)+3*f(i-2)-4.d0/3.d0*f(i-3)+0.25d0*f(i-4))/dX   ! 向后差分，4价精度
  enddo

#elif (ghost_width == 4)
!#error
  do i=1,lx-7
    fx(i) = (-147.d0*f(i)+360.d0*f(i+1)-450.d0*f(i+2)+400.d0*f(i+3)-225.d0*f(i+4)+72.d0*f(i+5)-10.d0*f(i+6)) /60.d0/dX  ! 向前差分，6价精度
  enddo
  do i=lx-6,lx-3  
    fx(i) = (f(i-3)-9.d0*f(i-2)+45.d0*f(i-1)-45.d0*f(i+1)+9.0*f(i+2)-f(i+3)) /60.d0/dX                                  ! 中心差分，6价精度
  enddo
  do i=lx-2,lx
    fx(i) = (147.d0*f(i)-360.d0*f(i-1)+450.d0*f(i-2)-400.d0*f(i-3)+225.d0*f(i-4)-72.d0*f(i-5)+10.d0*f(i-6)) /60.d0/dX   ! 向后差分，6价精度
  enddo
  
#elif (ghost_width == 5)
!#error
  do i=1,lx-9    
    fx(i) = (-(761.d0/280.d0)*f(i)+8.d0*f(i+1)-14.d0*f(i+2)+(56.d0/3.d0)*f(i+3)-(35.d0/2.d0)*f(i+4)  &
             +(56.d0/5.d0)*f(i+5)-(14.d0/3.d0)*f(i+6)+(8.d0/7.d0)*f(i+7)-(1.d0/8.d0)*f(i+8)) /dX           
            ! 向前差分，8价精度
  enddo
  do i=lx-8,lx-4          
    fx(i) = (3.d0*f(i-4)-32.d0*f(i-3)+168.d0*f(i-2)-672.d0*f(i-1)+672.d0*f(i+1)-168.d0*f(i+2)+32.d0*f(i+3)-3.d0*f(i+4)) /840.d0/dX    
            ! 中心差分，8价精度
  enddo
  do i=lx-3,lx
    fx(i) = ((761.d0/280.d0)*f(i)-8.d0*f(i-1)+14.d0*f(i-2)-(56.d0/3.d0)*f(i-3)+(35.d0/2.d0)*f(i-4)  &
             -(56.d0/5.d0)*f(i-5)+(14.d0/3.d0)*f(i-6)-(8.d0/7.d0)*f(i-7)+(1.d0/8.d0)*f(i-8)) /dX           
            ! 向后差分，8价精度
  enddo

#endif  

  return

  end subroutine rderivs_sw_x
!----
  subroutine cderivs_sw_x(lx,X,f,fx)

  implicit none

!~~~~~~% Input parameters:
  integer,intent(in) :: lx
  real*8,intent(in),dimension(lx) :: X
  double complex,intent(in),dimension(lx) :: f
  double complex,intent(out),dimension(lx) :: fx

  real*8 :: dX
  integer :: i

  dX = X(2)-X(1)

#if (ghost_width == 2)
!#error
  do i=1,lx-3
    fx(i) = (-3.d0*f(i)+4.d0*f(i+1)-f(i+2))/2.d0/dX  ! 向前差分，2价精度
  enddo
  do i=lx-2,lx-1                                     ! 宽度为 2
    fx(i) = (f(i+1)-f(i-1))/2.d0/dX                  ! 中心差分，2价精度
  enddo
  i=lx
    fx(i) = (3.d0*f(i)-4.d0*f(i-1)+f(i-2))/2.d0/dX   ! 向后差分，2价精度
    
#elif (ghost_width == 3)
  do i=1,lx-5
    fx(i) = (-2.5d1/1.2d1*f(i)+4*f(i+1)-3*f(i+2)+4.d0/3.d0*f(i+3)-0.25d0*f(i+4))/dX  ! 向前差分，4价精度 
  enddo
  do i=lx-4,lx-2                                                                     ! 宽度为 3
    fx(i) = (f(i-2)-8.d0*f(i-1)+8.d0*f(i+1)-f(i+2))/1.2d1/dX                         ! 中心差分，4价精度 
  enddo
  do i=lx-1,lx
    fx(i) = (2.5d1/1.2d1*f(i)-4*f(i-1)+3*f(i-2)-4.d0/3.d0*f(i-3)+0.25d0*f(i-4))/dX   ! 向后差分，4价精度
  enddo
#elif (ghost_width == 4)
!#error
  do i=1,lx-7
    fx(i) = (-147.d0*f(i)+360.d0*f(i+1)-450.d0*f(i+2)+400.d0*f(i+3)-225.d0*f(i+4)+72.d0*f(i+5)-10.d0*f(i+6)) /60.d0/dX 
            ! 向前差分，6价精度
  enddo
  do i=lx-6,lx-3                                                                        ! 宽度为 4
    fx(i) = (f(i-3)-9.d0*f(i-2)+45.d0*f(i-1)-45.d0*f(i+1)+9.0*f(i+2)-f(i+3)) /60.d0/dX  ! 中心差分，6价精度
  enddo
  do i=lx-2,lx
    fx(i) = (147.d0*f(i)-360.d0*f(i-1)+450.d0*f(i-2)-400.d0*f(i-3)+225.d0*f(i-4)-72.d0*f(i-5)+10.d0*f(i-6)) /60.d0/dX  
            ! 向后差分，6价精度
  enddo
  
#elif (ghost_width == 5)
!#error
  do i=1,lx-9    
    fx(i) = (-(761.d0/280.d0)*f(i)+8.d0*f(i+1)-14.d0*f(i+2)+(56.d0/3.d0)*f(i+3)-(35.d0/2.d0)*f(i+4)  &
             +(56.d0/5.d0)*f(i+5)-(14.d0/3.d0)*f(i+6)+(8.d0/7.d0)*f(i+7)-(1.d0/8.d0)*f(i+8)) /dX         
            ! 向前差分，8价精度
  enddo
  do i=lx-8,lx-4        ! 宽度为 5
    fx(i) = (3.d0*f(i-4)-32.d0*f(i-3)+168.d0*f(i-2)-672.d0*f(i-1)+672.d0*f(i+1)-168.d0*f(i+2)+32.d0*f(i+3)-3.d0*f(i+4)) /840.d0/dX  
            ! 中心差分，8价精度
  enddo
  do i=lx-3,lx
    fx(i) = ((761.d0/280.d0)*f(i)-8.d0*f(i-1)+14.d0*f(i-2)-(56.d0/3.d0)*f(i-3)+(35.d0/2.d0)*f(i-4)  &
             -(56.d0/5.d0)*f(i-5)+(14.d0/3.d0)*f(i-6)-(8.d0/7.d0)*f(i-7)+(1.d0/8.d0)*f(i-8)) /dX         
            ! 向后差分，8价精度
  enddo
  
#endif  

  return

  end subroutine cderivs_sw_x
!----
  subroutine cdderivs_sw_x(lx,X,f,fxx)

  implicit none

!~~~~~~% Input parameters:
  integer,intent(in) :: lx
  real*8,intent(in),dimension(lx) :: X
  double complex,intent(in),dimension(lx) :: f
  double complex,intent(out),dimension(lx) :: fxx

  real*8 :: dX
  integer :: i

  dX = X(2)-X(1)
  dX = dX*dX

#if (ghost_width == 2)
!#error
  do i=1,lx-3
    fxx(i) = (f(i)-2.d0*f(i+1)+f(i+2))/dX  ! 向前差分，2价精度
  enddo
  do i=lx-2,lx-1
    fxx(i) = (f(i-1)-2.d0*f(i)+f(i+1))/dX  ! 中心差分，2价精度
  enddo
  i=lx
    fxx(i) = (f(i)-2.d0*f(i-1)+f(i-2))/dX  ! 向后差分，2价精度
    
#elif (ghost_width == 3)
  do i=1,lx-5
    fxx(i) = ((15.d0/4)*f(i)-(77.d0/6)*f(i+1)+(107.d0/6)*f(i+2)-13.d0*f(i+3)+(61.d0/12.d0)*f(i+4)-(5.d0/6.d0)*f(i+5))/dX ! 向前差分，4价精度
  enddo
  do i=lx-4,lx-2         ! 宽度为 3
    fxx(i) = (-f(i-2)+16.d0*f(i-1)-30.d0*f(i)+16.d0*f(i+1)-f(i+2))/12.d0/dX                                              ! 中心差分，4价精度
  enddo
  do i=lx-1,lx
    fxx(i) = ((15.d0/4)*f(i)-(77.d0/6)*f(i-1)+(107.d0/6)*f(i-2)-13.d0*f(i-3)+(61.d0/12.d0)*f(i-4)-(5.d0/6.d0)*f(i-5))/dX ! 向后差分，4价精度
  enddo
#elif (ghost_width == 4)
  do i=1,lx-7
    fxx(i) = (812.d0*f(i)-3132.d0*f(i+1)+5265.d0*f(i+2)-5080.d0*f(i+3)+2970.d0*f(i+4)-972.d0*f(i+5)+137.d0*f(i+6)) /180.d0/dX ! 向前差分，6价精度
  enddo
  do i=lx-6,lx-3         ! 宽度为 4
    fxx(i) = (2.d0*f(i-3)-27.d0*f(i-2)+270.d0*f(i-1)-490.d0*f(i)+270.d0*f(i+1)-27.d0*f(i+2)+2.d0*f(i+3)) /180.d0/dX           ! 中心差分，6价精度
  enddo
  do i=lx-2,lx
    fxx(i) = (812.d0*f(i)-3132.d0*f(i-1)+5265.d0*f(i-2)-5080.d0*f(i-3)+2970.d0*f(i-4)-972.d0*f(i-5)+137.d0*f(i-6)) /180.d0/dX ! 向后差分，6价精度
  enddo
  
#elif (ghost_width == 5)
  do i=1,lx-9
    fxx(i) = (   (29531.d0/5040.d0)*f(i) - (962.d0/35.d0)*f(i+1) + (621.d0/10.d0)*f(i+2)         &
               - (4006.d0/45.d0)*f(i+3)  + (691.d0/8.d0)*f(i+4)  - (282.d0/5.d0)*f(i+5)          &
               + (2143.d0/90.d0)*f(i+6)  - (206.d0/35.d0)*f(i+7) + (363.d0/560.d0)*f(i+8) ) /dX            
             ! 向前差分，8阶精度
  enddo
  do i=lx-8,lx-4            ! 宽度为 5
    fxx(i) = ( - 9.d0*f(i-4) + 128.d0*f(i-3) - 1008.d0*f(i-2) + 8064.d0*f(i-1)                   &
               - 14350.d0*f(i) + 8064.d0*f(i+1) - 1008.d0*f(i+2) + 128.d0*f(i+3) - 9.d0*f(i+4) ) / 5040.d0/dX   
             ! 中心差分，8阶精度
  enddo
  do i=lx-3,lx
    fxx(i) = (   (29531.d0/5040.d0)*f(i) - (962.d0/35.d0)*f(i-1) + (621.d0/10.d0)*f(i-2)         &
               - (4006.d0/45.d0)*f(i-3)  + (691.d0/8.d0)*f(i-4)  - (282.d0/5.d0)*f(i-5)          &
               + (2143.d0/90.d0)*f(i-6)  - (206.d0/35.d0)*f(i-7) + (363.d0/560.d0)*f(i-8) ) /dX            
             ! 向后差分，8阶精度
  enddo
#endif  

  return

  end subroutine cdderivs_sw_x

!--------------------
  subroutine rget_half_x(lx,f,hf)

  implicit none

!~~~~~~% Input parameters:

  integer,intent(in) :: lx
  real*8,intent(in),dimension(lx) :: f
  real*8,intent(out),dimension(lx) :: hf

#ifdef OLD
  hf(1:lx-1) = (f(2:lx)+f(1:lx-1))/2.d0
  hf(lx)     = hf(lx-1)

#else
  hf(lx) = 0.d0
#if (ghost_width == 2)
  hf(1:lx-1) = (f(1:lx-1)+f(2:lx))/2.d0
#elif (ghost_width == 3)
  hf(2:lx-2) = (-f(1:lx-3)+9.d0*f(2:lx-2)+9.d0*f(3:lx-1)-f(4:lx))/1.6d1
  hf(1)      = (f(2)+f(1))/2.d0
  hf(lx-1)   = (f(lx)+f(lx-1))/2.d0
#elif (ghost_width == 4)
  hf(3:lx-3) = 3.d0/2.56d2*(f(1:lx-5)+f(6:lx))-2.5d1/2.56d2*(f(2:lx-4)+f(5:lx-1))+7.5d1/1.28d2*(f(3:lx-3)+f(4:lx-2))
  hf(2)      = (-f(1)+9.d0*f(2)+9.d0*f(3)-f(4))/1.6d1
  hf(lx-2)   = (-f(lx-1)+9.d0*f(lx-2)+9.d0*f(lx-3)-f(lx-4))/1.6d1
  hf(1)      = (f(2)+f(1))/2.d0
  hf(lx-1)   = (f(lx)+f(lx-1))/2.d0
#elif (ghost_width == 5)
  hf(4:lx-4) = (-5.d0*(f(1:lx-7)+f(8:lx))+4.9d1*(f(2:lx-6)+f(7:lx-1))-2.45d2*(f(3:lx-5)+f(6:lx-2))+1.225d3*(f(4:lx-4)+f(5:lx-3))) &
               /2.048d3
  hf(3)      = 3.d0/2.56d2*(f(1)+f(6))-2.5d1/2.56d2*(f(2)+f(5))+7.5d1/1.28d2*(f(3)+f(4))
  hf(lx-3)   = 3.d0/2.56d2*(f(lx)+f(lx-5))-2.5d1/2.56d2*(f(lx-1)+f(lx-4))+7.5d1/1.28d2*(f(lx-2)+f(lx-3))
  hf(2)      = (-f(1)+9.d0*f(2)+9.d0*f(3)-f(4))/1.6d1
  hf(lx-2)   = (-f(lx-1)+9.d0*f(lx-2)+9.d0*f(lx-3)-f(lx-4))/1.6d1
  hf(1)      = (f(2)+f(1))/2.d0
  hf(lx-1)   = (f(lx)+f(lx-1))/2.d0
#endif  

#endif

  return

  end subroutine rget_half_x
!--------------------
  subroutine rget_half_x_point(lx,f,hf,k)

  implicit none

!~~~~~~% Input parameters:

  integer,intent(in) :: lx,k
  real*8,intent(in),dimension(lx) :: f
  real*8,intent(out) :: hf

#ifdef OLD
  if(k.eq.lx .or. k.eq.lx-1)then
    hf = (f(lx)+f(lx-1))/2.d0
  else
    hf = (f(k+1)+f(k))/2.d0
  endif

#else
  if(k .eq. lx) hf = 0.d0
#if (ghost_width == 2)
  if(k .lt. lx) hf = (f(k+1)+f(k))/2.d0
#elif (ghost_width == 3)
  if(k .eq. 1)  hf = (f(2)+f(1))/2.d0
  if(k .eq. lx-1) hf = (f(lx)+f(lx-1))/2.d0
  if(k .gt. 1 .and. k .lt. lx-1) hf = (-f(k-1)+9.d0*f(k)+9.d0*f(k+1)-f(k+2))/1.6d1
#elif (ghost_width == 4)
  if(k .eq. 1) hf = (f(2)+f(1))/2.d0
  if(k .eq. lx-1) hf = (f(lx)+f(lx-1))/2.d0
  if(k .eq. 2) hf = (-f(1)+9.d0*f(2)+9.d0*f(3)-f(4))/1.6d1
  if(k .eq. lx-2) hf =(-f(lx-1)+9.d0*f(lx-2)+9.d0*f(lx-3)-f(lx-4))/1.6d1
  if(k .gt. 2 .and. k .lt. lx-2) hf = 3.d0/2.56d2*(f(k-2)+f(k+3))-2.5d1/2.56d2*(f(k-1)+f(k+2))+7.5d1/1.28d2*(f(k)+f(k+1))
#elif (ghost_width == 5)
  if(k .eq. 1) hf = (f(2)+f(1))/2.d0
  if(k .eq. lx-1) hf = (f(lx)+f(lx-1))/2.d0
  if(k .eq. 2) hf = (-f(1)+9.d0*f(2)+9.d0*f(3)-f(4))/1.6d1
  if(k .eq. lx-2) hf =(-f(lx-1)+9.d0*f(lx-2)+9.d0*f(lx-3)-f(lx-4))/1.6d1
  if(k .eq. 3) hf = 3.d0/2.56d2*(f(1)+f(6))-2.5d1/2.56d2*(f(2)+f(5))+7.5d1/1.28d2*(f(3)+f(4))
  if(k .eq. lx-3) hf = 3.d0/2.56d2*(f(lx)+f(lx-5))-2.5d1/2.56d2*(f(lx-1)+f(lx-4))+7.5d1/1.28d2*(f(lx-2)+f(lx-3))
  if(k .gt. 3 .and. k .lt. lx-3) then
    hf = (-5.d0*(f(k-3)+f(k+4))+4.9d1*(f(k-2)+f(k+3))-2.45d2*(f(k-1)+f(k+2))+1.225d3*(f(k)+f(k+1)))/2.048d3
  endif
#endif  

#endif

  return

  end subroutine rget_half_x_point
!--------------------
  subroutine cget_half_x(lx,f,hf)

  implicit none

!~~~~~~% Input parameters:

  integer,intent(in) :: lx
  double complex,intent(in),dimension(lx) :: f
  double complex,intent(out),dimension(lx) :: hf

#ifdef OLD
  hf(1:lx-1) = (f(2:lx)+f(1:lx-1))/2.d0
  hf(lx)     = hf(lx-1)

#else
  hf(lx) = 0.d0
#if (ghost_width == 2)
  hf(1:lx-1) = (f(1:lx-1)+f(2:lx))/2.d0
#elif (ghost_width == 3)
  hf(2:lx-2) = (-f(1:lx-3)+9.d0*f(2:lx-2)+9.d0*f(3:lx-1)-f(4:lx))/1.6d1
  hf(1)      = (f(2)+f(1))/2.d0
  hf(lx-1)   = (f(lx)+f(lx-1))/2.d0
#elif (ghost_width == 4)
  hf(3:lx-3) = 3.d0/2.56d2*(f(1:lx-5)+f(6:lx))-2.5d1/2.56d2*(f(2:lx-4)+f(5:lx-1))+7.5d1/1.28d2*(f(3:lx-3)+f(4:lx-2))
  hf(2)      = (-f(1)+9.d0*f(2)+9.d0*f(3)-f(4))/1.6d1
  hf(lx-2)   = (-f(lx-1)+9.d0*f(lx-2)+9.d0*f(lx-3)-f(lx-4))/1.6d1
  hf(1)      = (f(2)+f(1))/2.d0
  hf(lx-1)   = (f(lx)+f(lx-1))/2.d0
#elif (ghost_width == 5)
  hf(4:lx-4) = (-5.d0*(f(1:lx-7)+f(8:lx))+4.9d1*(f(2:lx-6)+f(7:lx-1))-2.45d2*(f(3:lx-5)+f(6:lx-2))+1.225d3*(f(4:lx-4)+f(5:lx-3))) &
               /2.048d3
  hf(3)      = 3.d0/2.56d2*(f(1)+f(6))-2.5d1/2.56d2*(f(2)+f(5))+7.5d1/1.28d2*(f(3)+f(4))
  hf(lx-3)   = 3.d0/2.56d2*(f(lx)+f(lx-5))-2.5d1/2.56d2*(f(lx-1)+f(lx-4))+7.5d1/1.28d2*(f(lx-2)+f(lx-3))
  hf(2)      = (-f(1)+9.d0*f(2)+9.d0*f(3)-f(4))/1.6d1
  hf(lx-2)   = (-f(lx-1)+9.d0*f(lx-2)+9.d0*f(lx-3)-f(lx-4))/1.6d1
  hf(1)      = (f(2)+f(1))/2.d0
  hf(lx-1)   = (f(lx)+f(lx-1))/2.d0
#endif  

#endif
  return

  end subroutine cget_half_x
!---------
  subroutine derivs_eth(ex,X,Y,f,eth_f,spin,e,   &
                        quR1,quR2,quI1,quI2,gR,gI)

  implicit none

!~~~~~~% Input parameters:
! spin: spin weight of f; e: eth (1) or eth bar (-1)
  integer,intent(in) :: spin,e
  integer,intent(in),dimension(2) :: ex
  real*8,intent(in),dimension(ex(1)) :: X
  real*8,intent(in),dimension(ex(2)) :: Y
  double complex,intent(in),dimension(ex(1),ex(2)) :: f
  double complex,intent(out),dimension(ex(1),ex(2)) :: eth_f
  real*8,intent(in),dimension(ex(1),ex(2)) :: quR1,quR2,quI1,quI2,gR,gI
    
  double complex,dimension(ex(1),ex(2)) :: fx,fy
  double complex,dimension(ex(1),ex(2)) :: qu1,qu2,gama
  integer :: i,j

!sanity check  
  if(e.ne.1 .and.e.ne.-1)then
     write(*,*) "derivs_eth: bad input e = ",e
     return
  endif

  qu1 = dcmplx(quR1,quI1)
  qu2 = dcmplx(quR2,quI2)
  gama = dcmplx(gR,gI)

  do j=1,ex(2)
    call cderivs_x(ex(1),X,f(:,j),fx(:,j))
  enddo

  do i=1,ex(1)
    call cderivs_x(ex(2),Y,f(i,:),fy(i,:))
  enddo

  if(e.eq.1)then
    eth_f = qu1*fx+qu2*fy+spin*gama*f
  else
    eth_f = dconjg(qu1)*fx+dconjg(qu2)*fy-spin*dconjg(gama)*f
  endif

  return

  end subroutine derivs_eth
!---------
! dqu: q^B @_B q^A
! bdqu: bar{q}^B @_B q^A
! dg: q^B @_B Gamma
! bdg: q^B @_B bar{Gamma}
  subroutine dderivs_eth(ex,X,Y,f,eth_f,spin,e1,e2,   &
                         quR1,quR2,quI1,quI2,gR,gI, &
                         dquR1,dquR2,dquI1,dquI2,   &
                         bdquR1,bdquR2,bdquI1,bdquI2,   &
                         dgR,dgI,bdgR,bdgI)

  implicit none

!~~~~~~% Input parameters:
! spin: spin weight of f; e: eth (1) or eth bar (-1)
  integer,intent(in) :: spin,e1,e2
  integer,intent(in),dimension(2) :: ex
  real*8,intent(in),dimension(ex(1)) :: X
  real*8,intent(in),dimension(ex(2)) :: Y
  double complex,intent(in),dimension(ex(1),ex(2)) :: f
  double complex,intent(out),dimension(ex(1),ex(2)) :: eth_f
  real*8,intent(in),dimension(ex(1),ex(2)) :: quR1,quR2,quI1,quI2,gR,gI
  real*8,intent(in),dimension(ex(1),ex(2)) :: dquR1,dquR2,dquI1,dquI2
  real*8,intent(in),dimension(ex(1),ex(2)) :: bdquR1,bdquR2,bdquI1,bdquI2
  real*8,intent(in),dimension(ex(1),ex(2)) :: dgR,dgI,bdgR,bdgI
    
  double complex,dimension(ex(1),ex(2)) :: fx,fy,fxx,fxy,fyy
  double complex,dimension(ex(1),ex(2)) :: qu1,qu2,dqu1,dqu2,bdqu1,bdqu2,gama,dgama,bdgama
  integer :: i,j

!sanity check  
  if((e1.ne.1 .and. e1.ne.-1).or.(e2.ne.1 .and. e2.ne.-1))then
     write(*,*) "dderivs_eth: bad input e1 = ",e1,"e2 = ",e2
     return
  endif

  qu1 = dcmplx(quR1,quI1)
  qu2 = dcmplx(quR2,quI2)
  gama = dcmplx(gR,gI)
  dqu1 = dcmplx(dquR1,dquI1)
  dqu2 = dcmplx(dquR2,dquI2)
  bdqu1 = dcmplx(bdquR1,bdquI1)
  bdqu2 = dcmplx(bdquR2,bdquI2)
  dgama = dcmplx(dgR,dgI)
  bdgama = dcmplx(bdgR,bdgI)

  do j=1,ex(2)
    call cderivs_x(ex(1),X,f(:,j),fx(:,j))
    call cdderivs_x(ex(1),X,f(:,j),fxx(:,j))
  enddo

  do i=1,ex(1)
    call cderivs_x(ex(2),Y,f(i,:),fy(i,:))
    call cdderivs_x(ex(2),Y,f(i,:),fyy(i,:))
  enddo

  do j=1,ex(2)
    call cderivs_x(ex(1),X,fy(:,j),fxy(:,j))
  enddo

  if(e1.eq.1.and.e2.eq.1)then
! eth_eth          
    eth_f = qu1*qu1*fxx+2.d0*qu1*qu2*fxy+qu2*qu2*fyy &
           +(dqu1+(2*spin+1)*gama*qu1)*fx+(dqu2+(2*spin+1)*gama*qu2)*fy &
           +spin*(dgama+(spin+1)*gama*gama)*f
  elseif(e1.eq.1.and.e2.eq.-1)then
! eth_ethb          
    eth_f = qu1*dconjg(qu1)*fxx+qu1*dconjg(qu2)*fxy+qu2*dconjg(qu1)*fxy+qu2*dconjg(qu2)*fyy &
           +(dconjg(bdqu1)-spin*dconjg(gama)*qu1+(spin-1)*gama*dconjg(qu1))*fx &
           +(dconjg(bdqu2)-spin*dconjg(gama)*qu2+(spin-1)*gama*dconjg(qu2))*fy &
           -spin*(bdgama+(spin-1)*gama*dconjg(gama))*f
  elseif(e1.eq.-1.and.e2.eq.1)then
! ethb_eth
    eth_f = qu1*dconjg(qu1)*fxx+qu1*dconjg(qu2)*fxy+qu2*dconjg(qu1)*fxy+qu2*dconjg(qu2)*fyy &
           +(bdqu1+spin*gama*dconjg(qu1)-(spin+1)*dconjg(gama)*qu1)*fx &
           +(bdqu2+spin*gama*dconjg(qu2)-(spin+1)*dconjg(gama)*qu2)*fy &
           +spin*(dconjg(bdgama)-(spin+1)*gama*dconjg(gama))*f
  else
! ethb_ethb          
    eth_f = dconjg(qu1*qu1)*fxx+2.d0*dconjg(qu1*qu2)*fxy+dconjg(qu2*qu2)*fyy &
           +(dconjg(dqu1)-(2*spin-1)*dconjg(gama*qu1))*fx+(dconjg(dqu2)-(2*spin-1)*dconjg(gama*qu2))*fy &
           -spin*(dconjg(dgama)-(spin-1)*dconjg(gama*gama))*f
  endif

  return

  end subroutine dderivs_eth
!------------------------------------------------------------------------------------------
! CQG 24S327, Eq.(19)--(22)
!------------------------------------------------------------------------------------------
  subroutine setup_dyad(ex,crho,sigma,R,      &
                        quR1,quR2,quI1,quI2,  &
                        qlR1,qlR2,qlI1,qlI2,  &
                        gR,gI,                &
                        dquR1,dquR2,dquI1,dquI2,   &
                        bdquR1,bdquR2,bdquI1,bdquI2,   &
                        dgR,dgI,bdgR,bdgI,    &
                        gx,gy,gz,sst,Rmin)

  implicit none

!~~~~~~% Input parameters:
  integer,intent(in ):: ex(1:3),sst
  real*8,intent(in),dimension(ex(1))::crho
  real*8,intent(in),dimension(ex(2))::sigma
  real*8,intent(in),dimension(ex(3))::R
  real*8,intent(in) :: Rmin
  real*8,intent(out),dimension(ex(1),ex(2),ex(3)) :: quR1,quR2,quI1,quI2
  real*8,intent(out),dimension(ex(1),ex(2),ex(3)) :: qlR1,qlR2,qlI1,qlI2
  real*8,intent(out),dimension(ex(1),ex(2),ex(3)) :: gR,gI
  real*8,intent(out),dimension(ex(1),ex(2),ex(3)) :: dquR1,dquR2,dquI1,dquI2
  real*8,intent(out),dimension(ex(1),ex(2),ex(3)) :: bdquR1,bdquR2,bdquI1,bdquI2
  real*8,intent(out),dimension(ex(1),ex(2),ex(3)) :: dgR,dgI,bdgR,bdgI
  real*8,intent(out),dimension(ex(1),ex(2),ex(3)) :: gx,gy,gz
    
  real*8 :: thetac,thetas,sr,ss,cr,cs,srss,crcs,tcts,tcts2
  real*8 :: sr2,ss2,cr2,cs2,tc2,ts2
  integer :: i,j,k
  real*8 :: ggr,tgrho,tgsigma
  real*8,dimension(ex(1),ex(2),ex(3)) :: tp

  do j=1,ex(2)
  do i=1,ex(1)
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

     qlR1(i,j,:) = thetac*cs/4.d0/tcts2
     qlI1(i,j,:) = thetas*cs/4.d0/tcts2
     qlR2(i,j,:) = thetac*cr/4.d0/tcts2
     qlI2(i,j,:) =-thetas*cr/4.d0/tcts2

     quR1(i,j,:) = 2.d0*tcts*thetas/cs
     quI1(i,j,:) = 2.d0*tcts*thetac/cs
     quR2(i,j,:) = 2.d0*tcts*thetas/cr
     quI2(i,j,:) =-2.d0*tcts*thetac/cr

     gR(i,j,:) = (crcs*crcs*(sr+ss)+(cr*cr-cs*cs)*(ss-sr))/4.d0/thetac/crcs
     gI(i,j,:) = (crcs*crcs*(sr-ss)+(cs*cs-cr*cr)*(ss+sr))/4.d0/thetas/crcs

     dquR1(i,j,:) = 3.d0*cr2*ss*(cs2+cr2*ss2)/2.d0/cr/cs2
     dquI1(i,j,:) =-thetac*thetas*sr*(cs2+3.d0*cr2*ss2)/cr/cs2
     dquR2(i,j,:) = 3.d0*cs2*sr*(cr2+cs2*sr2)/2.d0/cs/cr2
     dquI2(i,j,:) = thetac*thetas*ss*(cr2+3.d0*cs2*sr2)/cs/cr2

     bdquR1(i,j,:) = sr*cs2*(cs2+cr2*ss2)/2.d0/cr/cs2
     bdquI1(i,j,:) = thetac*thetas*ss*(8.d0*tc2*ts2-3.d0*cs2*sr2-cr2)/cr/cs2
     bdquR2(i,j,:) = ss*cr2*(cr2+cs2*sr2)/2.d0/cs/cr2
     bdquI2(i,j,:) =-thetac*thetas*sr*(8.d0*tc2*ts2-3.d0*cr2*ss2-cs2)/cs/cr2

     ggr = tc2*ts2/cr2/cs2
     dgR(i,j,:) =-3.d0*ggr*sr*ss*(cs2+cr2)
     dgI(i,j,:) = 6.d0*ggr*thetac*thetas*(cs2-cr2)
     bdgR(i,j,:) = ggr*(2.d0*cr2*cs2+cr2+cs2)
     bdgI(i,j,:) = 0.d0

    do k=1,ex(3)
     ggr = R(k)*Rmin/(1.d0-R(k))
     tgrho = dtan(crho(i))
     tgsigma = dtan(sigma(j))
     select case (sst)
     case (0)
      gz(i,j,k) = ggr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      gx(i,j,k) = gz(i,j,k)*tgrho
      gy(i,j,k) = gz(i,j,k)*tgsigma
    case (1)
      gz(i,j,k) = -ggr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      gx(i,j,k) = gz(i,j,k)*tgrho
      gy(i,j,k) = gz(i,j,k)*tgsigma
    case (2)
      gx(i,j,k) = ggr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      gy(i,j,k) = gx(i,j,k)*tgrho
      gz(i,j,k) = gx(i,j,k)*tgsigma
    case (3)
      gx(i,j,k) = -ggr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      gy(i,j,k) = gx(i,j,k)*tgrho
      gz(i,j,k) = gx(i,j,k)*tgsigma
    case (4)
      gy(i,j,k) = ggr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      gx(i,j,k) = gy(i,j,k)*tgrho
      gz(i,j,k) = gy(i,j,k)*tgsigma
    case (5)
      gy(i,j,k) = -ggr/dsqrt(1+tgrho*tgrho+tgsigma*tgsigma)
      gx(i,j,k) = gy(i,j,k)*tgrho
      gz(i,j,k) = gy(i,j,k)*tgsigma
     case default
      write(*,*) "setup_dyad: not recognized sst = ",sst
      return
     end select
    enddo
  enddo
  enddo

! note the difference between right hand coordinate and left hand coordinate
  if(sst==1 .or. sst==3 .or. sst==4)then
     quI1 = -quI1
     quI2 = -quI2
     qlI1 = -qlI1
     qlI2 = -qlI2
     gI = -gI

     dquI1 = -dquI1
     dquI2 = -dquI2

     bdquI1 = -bdquI1
     bdquI2 = -bdquI2

     dgI = -dgI
     bdgI = -bdgI
  endif

  return

  end subroutine setup_dyad
!---------------------------------------------------------------------------------
! interface to c
!---------------------------------------------------------------------------------
subroutine eth_derivs(ex,X,Y,Rfi,Ifi,Reth_f,Ieth_f,spin,e,   &
                        quR1,quR2,quI1,quI2,gR,gI)

  implicit none

!~~~~~~% Input parameters:
! spin: spin weight of f; e: eth (1) or eth bar (-1)
  integer,intent(in) :: spin,e
  integer,intent(in),dimension(3) :: ex
  real*8,intent(in),dimension(ex(1)) :: X
  real*8,intent(in),dimension(ex(2)) :: Y
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) ::  Rfi,Ifi
  real*8,intent(out),dimension(ex(1),ex(2),ex(3)) :: Reth_f,Ieth_f
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: quR1,quR2,quI1,quI2,gR,gI
    
  double complex,dimension(ex(1),ex(2)) :: f,eth_f,fx,fy
  double complex,dimension(ex(1),ex(2)) :: qu1,qu2,gama
  integer :: k

  do k=1,ex(3)
     f = dcmplx(Rfi(:,:,k),Ifi(:,:,k))
     call derivs_eth(ex(1:2),X,Y,f,eth_f,spin,e,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     Reth_f(:,:,k) = dreal(eth_f)
     Ieth_f(:,:,k) = dimag(eth_f)
  enddo

  return

  end subroutine eth_derivs
!---------------------------------------------------------------------------------
subroutine eth_dderivs(ex,X,Y,Rfi,Ifi,Reth_f,Ieth_f,spin,e1,e2,   &
                       quR1,quR2,quI1,quI2,gR,gI, &
                       dquR1,dquR2,dquI1,dquI2,   &
                       bdquR1,bdquR2,bdquI1,bdquI2,   &
                       dgR,dgI,bdgR,bdgI)

  implicit none

!~~~~~~% Input parameters:
! spin: spin weight of f; e: eth (1) or eth bar (-1)
  integer,intent(in) :: spin,e1,e2
  integer,intent(in),dimension(3) :: ex
  real*8,intent(in),dimension(ex(1)) :: X
  real*8,intent(in),dimension(ex(2)) :: Y
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) ::  Rfi,Ifi
  real*8,intent(out),dimension(ex(1),ex(2),ex(3)) :: Reth_f,Ieth_f
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: quR1,quR2,quI1,quI2,gR,gI
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dquR1,dquR2,dquI1,dquI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: bdquR1,bdquR2,bdquI1,bdquI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: dgR,dgI,bdgR,bdgI
    
  double complex,dimension(ex(1),ex(2)) :: f,eth_f,fx,fy
  double complex,dimension(ex(1),ex(2)) :: qu1,qu2,gama
  integer :: k

  do k=1,ex(3)
     f = dcmplx(Rfi(:,:,k),Ifi(:,:,k))
     call dderivs_eth(ex(1:2),X,Y,f,eth_f,spin,e1,e2,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k), &
                     dquR1(:,:,k),dquR2(:,:,k),dquI1(:,:,k),dquI2(:,:,k), &
                     bdquR1(:,:,k),bdquR2(:,:,k),bdquI1(:,:,k),bdquI2(:,:,k), &
                     dgR(:,:,k),dgI(:,:,k),bdgR(:,:,k),bdgI(:,:,k))
     Reth_f(:,:,k) = dreal(eth_f)
     Ieth_f(:,:,k) = dimag(eth_f)
  enddo

  return

  end subroutine eth_dderivs
!---------------------------------------------------------------------------------
! fill symmetric boundary buffer points
!---------------------------------------------------------------------------------
subroutine fill_symmetric_boundarybuffer(ex,crho,sigma,R,drho,dsigma,   &
                        quR1,quR2,quI1,quI2,qlR1,qlR2,qlI1,qlI2,        &
                        Rv,Iv,Symmetry,sst,spin)

  implicit none

!~~~~~~% Input parameters:
  integer,intent(in),dimension(3) :: ex
  integer,intent(in) :: Symmetry,sst,spin
  real*8,intent(in),dimension(ex(1))::crho
  real*8,intent(in),dimension(ex(2))::sigma
  real*8,intent(in),dimension(ex(3))::R
  real*8,intent(in) :: drho,dsigma
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: quR1,quR2,quI1,quI2
  real*8,intent(in),dimension(ex(1),ex(2),ex(3)) :: qlR1,qlR2,qlI1,qlI2
  real*8,intent(inout),dimension(ex(1),ex(2),ex(3)) :: Rv,Iv
    
  integer :: i,j,k,t

  select case (Symmetry)
  case (0)
     return
  case (1)
     if((sst==2.or.sst==4).and.dabs(sigma(1)+ghost_width*dsigma) < dsigma/2.d0)then
       do k=1,ex(3)
       do j=1,ghost_width
       do i=1,ex(1)
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif 
          t = 2*ghost_width+2-j
#endif 
#ifdef Cell
#ifdef Vertex
#error Both Cell and Vertex are defined
#endif 
          t = 2*ghost_width+1-j
#endif      
          call symmetrymap(quR1(i,j,k),quR2(i,j,k),quI1(i,j,k),quI2(i,j,k),qlR1(i,j,k),qlR2(i,j,k),qlI1(i,j,k),qlI2(i,j,k), &
                      quR1(i,t,k),quR2(i,t,k),quI1(i,t,k),quI2(i,t,k),qlR1(i,t,k),qlR2(i,t,k),qlI1(i,t,k),qlI2(i,t,k), &
                      Rv(i,j,k),Iv(i,j,k),Rv(i,t,k),Iv(i,t,k),2,spin)                
       enddo
       enddo
       enddo
     endif
     if((sst==3.or.sst==5).and.dabs(sigma(ex(2))-ghost_width*dsigma) < dsigma/2.d0)then
       do k=1,ex(3)
       do j=ex(2)-ghost_width+1,ex(2)
       do i=1,ex(1)
          t = ex(2)-j+1
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif 
          t = ex(2)-2*ghost_width-1+t
#endif 
#ifdef Cell
#ifdef Vertex
#error Both Cell and Vertex are defined
#endif   
          t = ex(2)-2*ghost_width+t
#endif        
          call symmetrymap(quR1(i,j,k),quR2(i,j,k),quI1(i,j,k),quI2(i,j,k),qlR1(i,j,k),qlR2(i,j,k),qlI1(i,j,k),qlI2(i,j,k), &
                      quR1(i,t,k),quR2(i,t,k),quI1(i,t,k),quI2(i,t,k),qlR1(i,t,k),qlR2(i,t,k),qlI1(i,t,k),qlI2(i,t,k), &
                      Rv(i,j,k),Iv(i,j,k),Rv(i,t,k),Iv(i,t,k),2,spin) 
       enddo
       enddo
       enddo
     endif
  case (2)
    if(dabs(crho(1)+ghost_width*drho) < drho/2.d0)then
     if(dabs(sigma(1)+ghost_width*dsigma) < dsigma/2.d0)then
       do k=1,ex(3)
       do j=1,ghost_width
       do i=ghost_width+1,ex(1)
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif 
          t = 2*ghost_width+2-j
#endif 
#ifdef Cell
#ifdef Vertex
#error Both Cell and Vertex are defined
#endif 
          t = 2*ghost_width+1-j
#endif      
          call symmetrymap(quR1(i,j,k),quR2(i,j,k),quI1(i,j,k),quI2(i,j,k),qlR1(i,j,k),qlR2(i,j,k),qlI1(i,j,k),qlI2(i,j,k), &
                      quR1(i,t,k),quR2(i,t,k),quI1(i,t,k),quI2(i,t,k),qlR1(i,t,k),qlR2(i,t,k),qlI1(i,t,k),qlI2(i,t,k), &
                      Rv(i,j,k),Iv(i,j,k),Rv(i,t,k),Iv(i,t,k),2,spin)              
       enddo
       enddo
       enddo
      endif
       do k=1,ex(3)
       do j=1,ex(2)
       do i=1,ghost_width
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif 
          t = 2*ghost_width+2-i
#endif 
#ifdef Cell
#ifdef Vertex
#error Both Cell and Vertex are defined
#endif 
          t = 2*ghost_width+1-i
#endif      
          call symmetrymap(quR1(i,j,k),quR2(i,j,k),quI1(i,j,k),quI2(i,j,k),qlR1(i,j,k),qlR2(i,j,k),qlI1(i,j,k),qlI2(i,j,k), &
                      quR1(t,j,k),quR2(t,j,k),quI1(t,j,k),quI2(t,j,k),qlR1(t,j,k),qlR2(t,j,k),qlI1(t,j,k),qlI2(t,j,k), &
                      Rv(i,j,k),Iv(i,j,k),Rv(t,j,k),Iv(t,j,k),1,spin)              
       enddo
       enddo
       enddo
    else
     if(dabs(sigma(1)+ghost_width*dsigma) < dsigma/2.d0)then
       do k=1,ex(3)
       do j=1,ghost_width
       do i=1,ex(1)
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif 
          t = 2*ghost_width+2-j
#endif 
#ifdef Cell
#ifdef Vertex
#error Both Cell and Vertex are defined
#endif 
          t = 2*ghost_width+1-j
#endif      
          call symmetrymap(quR1(i,j,k),quR2(i,j,k),quI1(i,j,k),quI2(i,j,k),qlR1(i,j,k),qlR2(i,j,k),qlI1(i,j,k),qlI2(i,j,k), &
                      quR1(i,t,k),quR2(i,t,k),quI1(i,t,k),quI2(i,t,k),qlR1(i,t,k),qlR2(i,t,k),qlI1(i,t,k),qlI2(i,t,k), &
                      Rv(i,j,k),Iv(i,j,k),Rv(i,t,k),Iv(i,t,k),2,spin)               
       enddo
       enddo
       enddo
     endif
    endif
  end select

  return

  end subroutine fill_symmetric_boundarybuffer
!-------------------------------------------------------------------------------------------
! note Theta has spin weight 2, but the its determinant does not satisfy our
! assumption, so in principle this routine can not be applied to Theta
! but since we never calculate eth derivative on Theta, it does not matter
!-------------------------------------------------------------------------------------------
  subroutine symmetrymap(oquR1,oquR2,oquI1,oquI2,oqlR1,oqlR2,oqlI1,oqlI2,       &
                         iquR1,iquR2,iquI1,iquI2,iqlR1,iqlR2,iqlI1,iqlI2,       &
                         oRv,oIv,iRv,iIv,ind,spin)

  implicit none

!~~~~~~% Input parameters:
  integer,intent(in) :: ind,spin
  real*8,intent(in) :: iquR1,iquR2,iquI1,iquI2,iqlR1,iqlR2,iqlI1,iqlI2
  real*8,intent(in) :: iRv,iIv
  real*8,intent(in) :: oquR1,oquR2,oquI1,oquI2,oqlR1,oqlR2,oqlI1,oqlI2
  real*8,intent(out) :: oRv,oIv

  double complex :: iqu1,iqu2,iql1,iql2,oqu1,oqu2,oql1,oql2
  real*8 :: h11,h12,h22,KK
  double complex :: r

  iqu1 = dcmplx(iquR1,iquI1)
  iqu2 = dcmplx(iquR2,iquI2)
  iql1 = dcmplx(iqlR1,iqlI1)
  iql2 = dcmplx(iqlR2,iqlI2)

  oqu1 = dcmplx(oquR1,oquI1)
  oqu2 = dcmplx(oquR2,oquI2)
  oql1 = dcmplx(oqlR1,oqlI1)
  oql2 = dcmplx(oqlR2,oqlI2)

! spin weight
  select case (spin)
  case (-2)
    r = dcmplx(iRv,iIv)
#if 1
    KK=dsqrt(iRv*iRv+iIv*iIv+1.d0)
    h11 = dreal(r*iql1*iql1+KK*dconjg(iql1)*iql1)
    h12 =-dreal(r*iql1*iql2+KK*dconjg(iql1)*iql2)
    h22 = dreal(r*iql2*iql2+KK*dconjg(iql2)*iql2)
#else
    h11 = dreal(r*iql1*iql1)
    h12 =-dreal(r*iql1*iql2)
    h22 = dreal(r*iql2*iql2)
#endif
    r = h11*dconjg(oqu1)*dconjg(oqu1)+2.d0*h12*dconjg(oqu1)*dconjg(oqu2)+h22*dconjg(oqu2)*dconjg(oqu2)
    oRv = dreal(r)
    oIv = dimag(r)
  case (-1)
    r = dcmplx(iRv,iIv)
    h11 = dreal(r*iql1)
    h22 = dreal(r*iql2)
    if(ind == 1) h11 = -h11
    if(ind == 2) h22 = -h22
    r = h11*dconjg(oqu1)+h22*dconjg(oqu2)
    oRv = dreal(r)
    oIv = dimag(r)
  case (0)
    oRv = iRv
    oIv = iIv
  case(1) 
    r = dcmplx(iRv,iIv)
    h11 = dreal(r*dconjg(iql1))
    h22 = dreal(r*dconjg(iql2))
    if(ind == 1) h11 = -h11
    if(ind == 2) h22 = -h22
    r = h11*oqu1+h22*oqu2
    oRv = dreal(r)
    oIv = dimag(r)
  case (2)
    r = dcmplx(iRv,iIv)
#if 1
    KK=dsqrt(iRv*iRv+iIv*iIv+1.d0)
    h11 = dreal(r*dconjg(iql1)*dconjg(iql1)+KK*dconjg(iql1)*iql1)
    h12 =-dreal(r*dconjg(iql1)*dconjg(iql2)+KK*dconjg(iql1)*iql2)
    h22 = dreal(r*dconjg(iql2)*dconjg(iql2)+KK*dconjg(iql2)*iql2)
#else    
    h11 = dreal(r*dconjg(iql1)*dconjg(iql1))
    h12 =-dreal(r*dconjg(iql1)*dconjg(iql2))
    h22 = dreal(r*dconjg(iql2)*dconjg(iql2))
#endif
    r = h11*oqu1*oqu1+2.d0*h12*oqu1*oqu2+h22*oqu2*oqu2
    oRv = dreal(r)
    oIv = dimag(r)
  case default
    write(*,*)"symmetrymap does not yet support spin weight ",spin
  end select

  return

  end subroutine symmetrymap
!-------------------------
  subroutine calculate_K(ex,X,Y,Z,RJ,IJ,KK,HKK,KKx,HKKx)

  implicit none

!~~~~~~% Input parameters:
  integer,intent(in) :: ex(3)
  real*8,intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8,dimension(ex(1),ex(2),ex(3)),intent(in) :: RJ,IJ
  real*8,dimension(ex(1),ex(2),ex(3)),intent(out) :: KK,HKK,KKx,HKKx

  integer :: i,j

  KK = dsqrt(1.d0+RJ*RJ+IJ*IJ)

  do j=1,ex(2)
  do i=1,ex(1)
     call rget_half_x(ex(3),KK(i,j,:),HKK(i,j,:))
     call rderivs_x(ex(3),Z,KK(i,j,:),KKx(i,j,:))
     call rget_half_x(ex(3),KKx(i,j,:),HKKx(i,j,:))
  enddo
  enddo

  return

  end subroutine calculate_K  
!--------------------------------------------------------------------
! this R is indeed x
function Eq_Theta(ex,crho,sigma,R,RJ,IJ,RU,IU,beta,RB,IB, &
                        Rnu,Inu,Rk,Ik,RTheta,ITheta,W,Rmin,       &
                        qlR1,qlR2,qlI1,qlI2,quR1,quR2,quI1,quI2,gR,gI) result(gont)   
  implicit none
  integer,intent(in ):: ex(1:3)
  real*8,intent(in),dimension(ex(1))::crho
  real*8,intent(in),dimension(ex(2))::sigma
  real*8,intent(in),dimension(ex(3))::R
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: beta,W
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: RJ,IJ,RU,IU,RB,IB,Rnu,Inu,Rk,Ik
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: RTheta,ITheta
  real*8,intent(in) :: Rmin
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: qlR1,qlR2,qlI1,qlI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: quR1,quR2,quI1,quI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: gR,gI
!  gont = 0: success; gont = 1: something wrong
  integer::gont

  double complex,dimension(ex(1),ex(2),ex(3)) :: CU,DCU,bDCU,CB,DCB,bDCB,CJ,DCJ
  double complex :: CTheta0,CTheta,CTheta1,RHS
  integer :: i,j,k,RK4
  double complex,dimension(ex(3)) :: Cnu,Ck,CUx,DCUx,bDCUx
  double complex,dimension(ex(3)) :: Cnux,CJx,CJxx,DCJx
  double complex,dimension(ex(3)) :: fCTheta,CThetax
  real*8,dimension(ex(3)) :: KK,KKx,betax,Wx
  double complex :: Theta_rhs,Theta_rhs_o
  real*8 :: dR

!!! sanity check
  dR = sum(RJ)+sum(IJ)+sum(RU)+sum(IU)+sum(beta)+sum(RB)+sum(IB) + &
       sum(Rnu)+sum(Inu)+sum(Rk)+sum(Ik)+sum(RTheta)+sum(ITheta)
  if(dR.ne.dR) then
     if(sum(RJ).ne.sum(RJ))write(*,*)"NullEvol_Theta: find NaN in RJ"
     if(sum(IJ).ne.sum(IJ))write(*,*)"NullEvol_Theta: find NaN in IJ"
     if(sum(RU).ne.sum(RU))write(*,*)"NullEvol_Theta: find NaN in RU"
     if(sum(IU).ne.sum(IU))write(*,*)"NullEvol_Theta: find NaN in IU"
     if(sum(beta).ne.sum(beta))write(*,*)"NullEvol_Theta: find NaN in beta"
     if(sum(RB).ne.sum(RB))write(*,*)"NullEvol_Theta: find NaN in RB"
     if(sum(IB).ne.sum(IB))write(*,*)"NullEvol_Theta: find NaN in IB"
     if(sum(Rnu).ne.sum(Rnu))write(*,*)"NullEvol_Theta: find NaN in Rnu"
     if(sum(Inu).ne.sum(Inu))write(*,*)"NullEvol_Theta: find NaN in Inu"
     if(sum(Rk).ne.sum(Rk))write(*,*)"NullEvol_Theta: find NaN in Rk"
     if(sum(Ik).ne.sum(Ik))write(*,*)"NullEvol_Theta: find NaN in Ik"
     if(sum(RTheta).ne.sum(RTheta))write(*,*)"NullEvol_Theta: find NaN in RTheta"
     if(sum(ITheta).ne.sum(ITheta))write(*,*)"NullEvol_Theta: find NaN in ITheta"
     gont = 1
     return
  endif

  dR = R(2) - R(1)
  
  CU = dcmplx(RU,IU)
  CB = dcmplx(RB,IB)
  CJ = dcmplx(RJ,IJ)

  do k=1,ex(3)
     call derivs_eth(ex(1:2),crho,sigma,CU(:,:,k),DCU(:,:,k),1,1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CU(:,:,k),bDCU(:,:,k),1,-1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CB(:,:,k),DCB(:,:,k),1,1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CB(:,:,k),bDCB(:,:,k),1,-1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CJ(:,:,k),DCJ(:,:,k),2,1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
  enddo

  do j=1,ex(2)
  do i=1,ex(1)
     CTheta0 = dcmplx(RTheta(i,j,1),ITheta(i,j,1))
     fCTheta = dcmplx(RTheta(i,j,:),ITheta(i,j,:))
     call cderivs_x(ex(3),R,fCTheta,CThetax)
     Cnu = dcmplx(Rnu(i,j,:),Inu(i,j,:))
     Ck = dcmplx(Rk(i,j,:),Ik(i,j,:))
     call cderivs_x(ex(3),R,Cnu,Cnux)
     call rderivs_x(ex(3),R,beta(i,j,:),betax)
     KK = dsqrt(1.d0+RJ(i,j,:)*RJ(i,j,:)+IJ(i,j,:)*IJ(i,j,:))
     call rderivs_x(ex(3),R,KK,KKx)
     call rderivs_x(ex(3),R,W,Wx)
     call cderivs_x(ex(3),R,DCU(i,j,:),DCUx)
     call cderivs_x(ex(3),R,CU(i,j,:),CUx)
     call cderivs_x(ex(3),R,bDCU(i,j,:),bDCUx)
     call cderivs_x(ex(3),R,CJ(i,j,:),CJx)
     call cdderivs_x(ex(3),R,CJ(i,j,:),CJxx)
     call cderivs_x(ex(3),R,DCJ(i,j,:),DCJx)

     RTheta(i,j,1) = 0.d0
     ITheta(i,j,1) = 0.d0
     do k=1,ex(3)-1
#if 0        
        call check_daxiao(beta(i,j,k))
        call check_daxiao(CU(i,j,k))
        call check_daxiao(bDCU(i,j,k))
        call check_daxiao(DCU(i,j,k))
        call check_daxiao(CB(i,j,k))
        call check_daxiao(DCB(i,j,k))
        call check_daxiao(W(i,j,k))
        call check_daxiao(CJ(i,j,k))
        call check_daxiao(DCJ(i,j,k))
        call check_daxiao(bDCB(i,j,k))
        call check_daxiao(Cnu(k))
        call check_daxiao(Ck(k))
        call check_daxiao(fCTheta(k))

        call check_daxiao(betax(k))
        call check_daxiao(KKx(k))
        call check_daxiao(CUx(k))
        call check_daxiao(DCUx(k))
        call check_daxiao(bDCUx(k))
        call check_daxiao(Wx(k))
        call check_daxiao(CJx(k))
        call check_daxiao(CJxx(k))
        call check_daxiao(DCJx(k))
        call check_daxiao(Cnux(k))
#endif
        RHS = Theta_rhs(R(k),Rmin,beta(i,j,k),betax(k),KK(k),KKx(k),CU(i,j,k),CUx(k),DCUx(k),bDCU(i,j,k),bDCUx(k), &
                        DCU(i,j,k),CB(i,j,k),DCB(i,j,k),W(i,j,k),Wx(k),CJ(i,j,k),DCJ(i,j,k),    &
                        CJx(k),CJxx(k),DCJx(k),bDCB(i,j,k),Cnu(k),Cnux(k),Ck(k),fCTheta(k))
        RHS = RHS - CThetax(k)
        
        RTheta(i,j,k+1) = dreal(RHS)
        ITheta(i,j,k+1) = dimag(RHS)
     enddo
     write(*,*)RTheta(i,j,:)
     stop
  enddo
  enddo

  gont = 0
  return

end function Eq_Theta

subroutine check_daxiao(f)
implicit none
double complex,intent(in) :: f

real*8,parameter :: eps=1.d-5

if(cdabs(f)>eps) write(*,*) f

return

end subroutine check_daxiao
subroutine check_factor(T,crho,sigma,R,sst,Rmin)
implicit none
integer,intent(in) :: sst
real*8,intent(in) :: T,crho,sigma,R,Rmin

real*8 ::x,y,z,hgr,gt,gp,tgrho,tgsigma,tc,ts,rf
double complex :: Yslm,II,Jr,swtf,ff

    hgr = R*Rmin/(1.d0-R)
    tgrho = dtan(crho)
    tgsigma = dtan(sigma)
    tc = dsqrt((1.d0-dsin(crho)*dsin(sigma))/2.d0)
    ts = dsqrt((1.d0+dsin(crho)*dsin(sigma))/2.d0)
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

  write(*,*) dsqrt(dble((2-1)*2*(2+1)*(2+2)))*Yslm(2,2,0,gt,gp)*swtf**2

  return

  end subroutine check_factor

subroutine getdxs(T,crho,sigma,R,betax,KKx,CUx,DCUx,bDCUx,Wx,CJx,CJxx,DCJx,Cnux,CThetax,sst,Rmin)
implicit none
integer,intent(in) :: sst
real*8,intent(in) :: T,crho,sigma,R,Rmin
real*8,intent(out) :: betax,KKx,Wx
double complex,intent(out) :: CUx,DCUx,bDCUx,CJx,CJxx,DCJx,Cnux,CThetax

real*8 ::x,y,z,hgr,gt,gp,tgrho,tgsigma,tc,ts,rf
double complex :: Yslm,II,Jr,swtf,ff
double complex :: beta0,C1,C2
integer :: nu,m

  call initial_null_paramter(beta0,C1,C2,nu,m)

  II = dcmplx(0.d0,1.d0)

    hgr = R*Rmin/(1.d0-R)
    tgrho = dtan(crho)
    tgsigma = dtan(sigma)
    tc = dsqrt((1.d0-dsin(crho)*dsin(sigma))/2.d0)
    ts = dsqrt((1.d0+dsin(crho)*dsin(sigma))/2.d0)
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

  betax = 0.d0

  Jr = -(3.d0*II*nu*C1-6.d0*beta0-II*nu**3*C2)/3.d0/hgr**2&
        +2.d0*nu*nu*C2/hgr**3-3.d0*II*nu*C2/hgr**4-2.d0*C2/hgr**5
  Wx = dreal(Yslm(0,2,m,gt,gp))*dreal(Jr*cdexp(II*nu*T))*(Rmin+hgr)**2/Rmin
!  Wx = dreal(Jr*cdexp(II*nu*T))*(Rmin+hgr)**2/Rmin
  KKx = 0.d0

  Jr = -2.d0*beta0/hgr/hgr-C1/hgr**3-II*nu*C2/hgr**4-C2/hgr**5
  rf = dreal(Jr*cdexp(II*nu*T))
  CUx = dsqrt(dble(2*(2+1)))*Yslm(1,2,m,gt,gp)*swtf*rf*(Rmin+hgr)**2/Rmin
!  CUx = rf*(Rmin+hgr)**2/Rmin
  DCUx = dsqrt(dble((2-1)*2*(2+1)*(2+2)))*Yslm(2,2,m,gt,gp)*swtf**2*rf*(Rmin+hgr)**2/Rmin
!  DCUx = rf*(Rmin+hgr)**2/Rmin
  bDCUx =-dble(2*(2+1))*Yslm(0,2,m,gt,gp)*rf*(Rmin+hgr)**2/Rmin
!  bDCUx = rf*(Rmin+hgr)**2/Rmin

  Jr = -C1/4.d0/hgr**2+C2/4.d0/hgr**4
  rf = dreal(Jr*cdexp(II*nu*T))
  CJx = dsqrt(dble((2-1)*2*(2+1)*(2+2)))*Yslm(2,2,m,gt,gp)*swtf**2*rf*(Rmin+hgr)**2/Rmin
!  CJx = rf*(Rmin+hgr)**2/Rmin
  Cnux =-dble((2-1)*(2+2))*dsqrt(dble(2*(2+1)))*Yslm(1,2,m,gt,gp)*swtf*rf*(Rmin+hgr)**2/Rmin
!  Cnux = rf*(Rmin+hgr)**2/Rmin
  DCJx = 0.d0
  rf = dreal(Jr*II*nu*cdexp(II*nu*T))
  CThetax = dsqrt(dble((2-1)*2*(2+1)*(2+2)))*Yslm(2,2,m,gt,gp)*swtf**2*rf*(Rmin+hgr)**2/Rmin
!  CThetax = rf*(Rmin+hgr)**2/Rmin
  Jr = C1/2.d0/hgr**3-C2/hgr**5
  rf = dreal(Jr*cdexp(II*nu*T))
  CJxx = dsqrt(dble((2-1)*2*(2+1)*(2+2)))*Yslm(2,2,m,gt,gp)*swtf**2*rf*(Rmin+hgr)**4/Rmin**2+2.d0*(Rmin+hgr)/Rmin*CJx
!  CJxx = rf*(Rmin+hgr)**4/Rmin**2+2.d0*(Rmin+hgr)/Rmin*CJx

#if 0
  DCUx = DCUx*dsqrt(dble((2-1)*2*(2+1)*(2+2)))*Yslm(2,2,m,gt,gp)*swtf**2
  CJx = CJx*dsqrt(dble((2-1)*2*(2+1)*(2+2)))*Yslm(2,2,m,gt,gp)*swtf**2
  CJxx = CJxx*dsqrt(dble((2-1)*2*(2+1)*(2+2)))*Yslm(2,2,m,gt,gp)*swtf**2
  CThetax = CThetax*dsqrt(dble((2-1)*2*(2+1)*(2+2)))*Yslm(2,2,m,gt,gp)*swtf**2
#endif
  return

  end subroutine getdxs

subroutine getndxs(T,crho,sigma,R,beta,KK,CU,bDCU,DCU,CB,DCB,W,CJ,DCJ,bDCB,Cnu,Ck,CTheta,sst,Rmin)
implicit none
integer,intent(in) :: sst
real*8,intent(in) :: T,crho,sigma,R,Rmin
real*8,intent(out) :: beta,KK,W
double complex,intent(out) :: CU,bDCU,DCU,CB,DCB,CJ,DCJ,bDCB,Cnu,Ck,CTheta

real*8 ::x,y,z,hgr,gt,gp,tgrho,tgsigma,tc,ts,rf
double complex :: Yslm,II,Jr,swtf,ff
double complex :: beta0,C1,C2
integer :: nu,m

  call initial_null_paramter(beta0,C1,C2,nu,m)

  II = dcmplx(0.d0,1.d0)

    hgr = R*Rmin/(1.d0-R)
    tgrho = dtan(crho)
    tgsigma = dtan(sigma)
    tc = dsqrt((1.d0-dsin(crho)*dsin(sigma))/2.d0)
    ts = dsqrt((1.d0+dsin(crho)*dsin(sigma))/2.d0)
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

  beta =  dreal(Yslm(0,2,m,gt,gp))*dreal(beta0*cdexp(II*nu*T))
!  beta =  dreal(beta0*cdexp(II*nu*T))
  CB = dsqrt(dble(2*(2+1)))*Yslm(1,2,m,gt,gp)*swtf*dreal(beta0*cdexp(II*nu*T))
!  CB = dreal(beta0*cdexp(II*nu*T))
  DCB = dsqrt(dble((2-1)*2*(2+1)*(2+2)))*Yslm(2,2,m,gt,gp)*swtf**2*dreal(beta0*cdexp(II*nu*T))
!  DCB = dreal(beta0*cdexp(II*nu*T))
  bDCB =-dble(2*(2+1))*Yslm(0,2,m,gt,gp)*dreal(beta0*cdexp(II*nu*T))
!  bDCB = dreal(beta0*cdexp(II*nu*T))

  Jr = (2.4d1*II*nu*beta0-3.d0*nu*nu*C1+nu**4*C2)/6.d0+(3.d0*II*nu*C1-6.d0*beta0-II*nu**3*C2)/3.d0/hgr&
        -nu*nu*C2/hgr/hgr+II*nu*C2/hgr**3+C2/2.d0/hgr**4
  W = dreal(Yslm(0,2,m,gt,gp))*dreal(Jr*cdexp(II*nu*T))
!  W = dreal(Jr*cdexp(II*nu*T))

  Jr = (2.4d1*beta0+3.d0*II*nu*C1-II*nu**3*C2)/3.6d1+C1/4.d0/hgr-C2/1.2d1/hgr**3
  rf = dreal(Jr*cdexp(II*nu*T))
  CJ = dsqrt(dble((2-1)*2*(2+1)*(2+2)))*Yslm(2,2,m,gt,gp)*swtf**2*rf
!  CJ = rf
  DCJ = 0.d0
  Cnu =-dsqrt(dble((2+2)*(2-2+1)*(2-1)*2*(2+1)*(2+2)))*Yslm(1,2,m,gt,gp)*swtf*rf
!  Cnu = rf
  KK = dsqrt(1.d0+cdabs(CJ)**2)
  Ck = 0.d0
  rf = dreal(Jr*II*nu*cdexp(II*nu*T))
  CTheta = dsqrt(dble((2-1)*2*(2+1)*(2+2)))*Yslm(2,2,m,gt,gp)*swtf**2*rf
!  CTheta = rf

  Jr = (-2.4d1*II*nu*beta0+3.d0*nu*nu*C1-nu**4*C2)/36.d0+2.d0*beta0/hgr&
        +C1/2.d0/hgr/hgr+II*nu*C2/3.d0/hgr**3+C2/4.d0/hgr**4
  rf = dreal(Jr*cdexp(II*nu*T))
  CU = dsqrt(dble(2*(2+1)))*Yslm(1,2,m,gt,gp)*swtf*rf
!  CU = rf
  DCU = dsqrt(dble((2-1)*2*(2+1)*(2+2)))*Yslm(2,2,m,gt,gp)*swtf**2*rf
!  DCU = rf
  bDCU =-dble(2*(2+1))*Yslm(0,2,m,gt,gp)*rf
!  bDCU = rf

#if 0
  DCU = DCU*dsqrt(dble((2-1)*2*(2+1)*(2+2)))*Yslm(2,2,m,gt,gp)*swtf**2
  DCB = DCB*dsqrt(dble((2-1)*2*(2+1)*(2+2)))*Yslm(2,2,m,gt,gp)*swtf**2
  CTheta = CTheta*dsqrt(dble((2-1)*2*(2+1)*(2+2)))*Yslm(2,2,m,gt,gp)*swtf**2
#endif
  return

  end subroutine getndxs
!--------------------------------------------------------------------
! this R is indeed x
function Eq_Theta_2(ex,crho,sigma,R,RJ,IJ,RU,IU,beta,RB,IB, &
                        Rnu,Inu,Rk,Ik,RTheta,ITheta,W,Rmin,       &
                        qlR1,qlR2,qlI1,qlI2,quR1,quR2,quI1,quI2,gR,gI, &
                        T,sst) result(gont)   
  implicit none
  integer,intent(in ):: ex(1:3),sst
  real*8,intent(in),dimension(ex(1))::crho
  real*8,intent(in),dimension(ex(2))::sigma
  real*8,intent(in),dimension(ex(3))::R
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: beta,W
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: RJ,IJ,RU,IU,RB,IB,Rnu,Inu,Rk,Ik
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: RTheta,ITheta
  real*8,intent(in) :: Rmin,T
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: qlR1,qlR2,qlI1,qlI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: quR1,quR2,quI1,quI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: gR,gI
!  gont = 0: success; gont = 1: something wrong
  integer::gont

  double complex,dimension(ex(1),ex(2),ex(3)) :: CU,DCU,bDCU,CB,DCB,bDCB,CJ,DCJ
  double complex :: CTheta0,CTheta,CTheta1,RHS
  integer :: i,j,k,RK4
  double complex,dimension(ex(3)) :: Cnu,Ck,HCnu,HCk,HCU,HDCU,CUx,HCUx,DCUx,HDCUx,HbDCU,bDCUx,HbDCUx
  double complex,dimension(ex(3)) :: Cnux,HCnux,HCJ,HDCJ,CJx,HCJx,CJxx,HCJxx,DCJx,HDCJx,HCB,HDCB,HbDCB
  double complex,dimension(ex(3)) :: fCTheta,CThetax
  real*8,dimension(ex(3)) :: KK,KKx,HKK,HKKx,Hbeta,betax,Hbetax,HW,Wx,HWx
  double complex :: Theta_rhs,Theta_rhs_o
  real*8 :: dR

!!! sanity check
  dR = sum(RJ)+sum(IJ)+sum(RU)+sum(IU)+sum(beta)+sum(RB)+sum(IB) + &
       sum(Rnu)+sum(Inu)+sum(Rk)+sum(Ik)+sum(RTheta)+sum(ITheta)
  if(dR.ne.dR) then
     if(sum(RJ).ne.sum(RJ))write(*,*)"NullEvol_Theta: find NaN in RJ"
     if(sum(IJ).ne.sum(IJ))write(*,*)"NullEvol_Theta: find NaN in IJ"
     if(sum(RU).ne.sum(RU))write(*,*)"NullEvol_Theta: find NaN in RU"
     if(sum(IU).ne.sum(IU))write(*,*)"NullEvol_Theta: find NaN in IU"
     if(sum(beta).ne.sum(beta))write(*,*)"NullEvol_Theta: find NaN in beta"
     if(sum(RB).ne.sum(RB))write(*,*)"NullEvol_Theta: find NaN in RB"
     if(sum(IB).ne.sum(IB))write(*,*)"NullEvol_Theta: find NaN in IB"
     if(sum(Rnu).ne.sum(Rnu))write(*,*)"NullEvol_Theta: find NaN in Rnu"
     if(sum(Inu).ne.sum(Inu))write(*,*)"NullEvol_Theta: find NaN in Inu"
     if(sum(Rk).ne.sum(Rk))write(*,*)"NullEvol_Theta: find NaN in Rk"
     if(sum(Ik).ne.sum(Ik))write(*,*)"NullEvol_Theta: find NaN in Ik"
     if(sum(RTheta).ne.sum(RTheta))write(*,*)"NullEvol_Theta: find NaN in RTheta"
     if(sum(ITheta).ne.sum(ITheta))write(*,*)"NullEvol_Theta: find NaN in ITheta"
     gont = 1
     return
  endif

  dR = R(2) - R(1)
  
  CU = dcmplx(RU,IU)
  CB = dcmplx(RB,IB)
  CJ = dcmplx(RJ,IJ)

  do k=1,ex(3)
     call derivs_eth(ex(1:2),crho,sigma,CU(:,:,k),DCU(:,:,k),1,1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CU(:,:,k),bDCU(:,:,k),1,-1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CB(:,:,k),DCB(:,:,k),1,1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CB(:,:,k),bDCB(:,:,k),1,-1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,CJ(:,:,k),DCJ(:,:,k),2,1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
  enddo

  do j=ghost_width+1,ex(2)-ghost_width
  do i=ghost_width+1,ex(1)-ghost_width
     CTheta0 = dcmplx(RTheta(i,j,1),ITheta(i,j,1))
     fCTheta = dcmplx(RTheta(i,j,:),ITheta(i,j,:))
     call cderivs_x(ex(3),R,fCTheta,CThetax)
     Cnu = dcmplx(Rnu(i,j,:),Inu(i,j,:))
     Ck = dcmplx(Rk(i,j,:),Ik(i,j,:))
     call cget_half_x(ex(3),CB(i,j,:),HCB)
     call cget_half_x(ex(3),DCB(i,j,:),HDCB)
     call cget_half_x(ex(3),bDCB(i,j,:),HbDCB)
     call cget_half_x(ex(3),Cnu,HCnu)
     call cderivs_x(ex(3),R,Cnu,Cnux)
     call cget_half_x(ex(3),Cnux,HCnux)
     call cget_half_x(ex(3),Ck,HCk)
     call rget_half_x(ex(3),beta(i,j,:),Hbeta)
     call rderivs_x(ex(3),R,beta(i,j,:),betax)
     call rget_half_x(ex(3),betax,Hbetax)
     KK = dsqrt(1.d0+RJ(i,j,:)*RJ(i,j,:)+IJ(i,j,:)*IJ(i,j,:))
     call rget_half_x(ex(3),KK,HKK)
     call rderivs_x(ex(3),R,KK,KKx)
     call rget_half_x(ex(3),KKx,HKKx)
     call rderivs_x(ex(3),R,W,Wx)
     call rget_half_x(ex(3),Wx,HWx)
     call rget_half_x(ex(3),W(i,j,:),HW)
     call cget_half_x(ex(3),CU(i,j,:),HCU)
     call cderivs_x(ex(3),R,DCU(i,j,:),DCUx)
     call cderivs_x(ex(3),R,CU(i,j,:),CUx)
     call cget_half_x(ex(3),DCUx,HDCUx)
     call cget_half_x(ex(3),CUx,HCUx)
     call cget_half_x(ex(3),DCU(i,j,:),HDCU)
     call cderivs_x(ex(3),R,bDCU(i,j,:),bDCUx)
     call cget_half_x(ex(3),bDCUx,HbDCUx)
     call cget_half_x(ex(3),bDCU(i,j,:),HbDCU)
     call cderivs_x(ex(3),R,CJ(i,j,:),CJx)
     call cdderivs_x(ex(3),R,CJ(i,j,:),CJxx)
     call cget_half_x(ex(3),CJx,HCJx)
     call cget_half_x(ex(3),CJxx,HCJxx)
     call cderivs_x(ex(3),R,DCJ(i,j,:),DCJx)
     call cget_half_x(ex(3),DCJx,HDCJx)

     RTheta(i,j,1) = 0.d0
     ITheta(i,j,1) = 0.d0
     do k=1,ex(3)-1
!        call getndxs(T,crho(i),sigma(j),R(k),beta(i,j,k),KK(k),CU(i,j,k),bDCU(i,j,k),DCU(i,j,k), &
!                     CB(i,j,k),DCB(i,j,k),W(i,j,k),CJ(i,j,k),DCJ(i,j,k),bDCB(i,j,k),Cnu(k),Ck(k),fCTheta(k),sst,Rmin)
!        call getdxs(T,crho(i),sigma(j),R(k),betax(k),KKx(k),CUx(k),DCUx(k),bDCUx(k), &
!                    Wx(k),CJx(k),CJxx(k),DCJx(k),Cnux(k),CThetax(k),sst,Rmin)
        RHS = Theta_rhs(R(k),Rmin,beta(i,j,k),betax(k),KK(k),KKx(k),CU(i,j,k),CUx(k),DCUx(k),bDCU(i,j,k),bDCUx(k), &
                        DCU(i,j,k),CB(i,j,k),DCB(i,j,k),W(i,j,k),Wx(k),CJ(i,j,k),DCJ(i,j,k),    &
                        CJx(k),CJxx(k),DCJx(k),bDCB(i,j,k),Cnu(k),Cnux(k),Ck(k),fCTheta(k))
        RHS = RHS - CThetax(k)
#if 0        
        if(cdabs(RHS)>1.d-9)then
#if 0                
                write(*,*)beta(i,j,k),KK(k),CU(i,j,k),bDCU(i,j,k),DCU(i,j,k),CB(i,j,k),DCB(i,j,k)
                write(*,*)W(i,j,k),CJ(i,j,k),DCJ(i,j,k),bDCB(i,j,k),Cnu(k),Ck(k),fCTheta(k)
        call getndxs(T,crho(i),sigma(j),R(k),beta(i,j,k),KK(k),CU(i,j,k),bDCU(i,j,k),DCU(i,j,k), &
                     CB(i,j,k),DCB(i,j,k),W(i,j,k),CJ(i,j,k),DCJ(i,j,k),bDCB(i,j,k),Cnu(k),Ck(k),fCTheta(k),sst,Rmin)
                write(*,*)"VS"
                write(*,*)beta(i,j,k),KK(k),CU(i,j,k),bDCU(i,j,k),DCU(i,j,k),CB(i,j,k),DCB(i,j,k)
                write(*,*)W(i,j,k),CJ(i,j,k),DCJ(i,j,k),bDCB(i,j,k),Cnu(k),Ck(k),fCTheta(k)
#endif               
                write(*,*)betax(k),KKx(k),CUx(k),DCUx(k),bDCUx(k)
                write(*,*)Wx(k),CJx(k),CJxx(k),DCJx(k),Cnux(k),CThetax(k)
        call getdxs(T,crho(i),sigma(j),R(k),betax(k),KKx(k),CUx(k),DCUx(k),bDCUx(k), &
                    Wx(k),CJx(k),CJxx(k),DCJx(k),Cnux(k),CThetax(k),sst,Rmin)
                write(*,*)"VS"
                write(*,*)betax(k),KKx(k),CUx(k),DCUx(k),bDCUx(k)
                write(*,*)Wx(k),CJx(k),CJxx(k),DCJx(k),Cnux(k),CThetax(k)
!                write(*,*)RHS
!                call check_factor(T,crho(i),sigma(j),R(k),sst,Rmin)
                stop
        endif
#endif        
        RTheta(i,j,k+1) = dreal(RHS)
        ITheta(i,j,k+1) = dimag(RHS)
     enddo
  enddo
  enddo

  gont = 0
  return

end function Eq_Theta_2
