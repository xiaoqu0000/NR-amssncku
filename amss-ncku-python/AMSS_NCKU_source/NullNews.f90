

#include "macrodef.fh"

!------------------------------------------------------------------------------
function omega_rhs(ex,crho,sigma,R,omega,RU,IU,omegarhs, &
                   quR1,quR2,quI1,quI2,gR,gI)    result(gont)   

  implicit none

  integer,intent(in) :: ex(3)
  real*8,intent(in),dimension(ex(1))::crho
  real*8,intent(in),dimension(ex(2))::sigma
  real*8,intent(in),dimension(ex(3))::R
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: omega,RU,IU
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: omegarhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: quR1,quR2,quI1,quI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: gR,gI
!  gont = 0: success; gont = 1: something wrong
  integer::gont

  double complex, dimension(ex(1),ex(2),ex(3)) :: comega,eth_omega,U,eth_Ub
  real*8 :: dR
  integer :: k

!!! sanity check
  dR = sum(omega)+sum(RU)+sum(IU)
  if(dR.ne.dR) then
     if(sum(omega).ne.sum(omega))write(*,*)"NullEvol_beta: find NaN in omega"
     if(sum(RU).ne.sum(RU))write(*,*)"NullEvol_beta: find NaN in RU"
     if(sum(IU).ne.sum(IU))write(*,*)"NullEvol_beta: find NaN in IU"
     gont = 1
     return
  endif

  comega = dcmplx(omega,0.d0)
  U = dcmplx(RU,IU)

  do k=1,ex(3)
     call derivs_eth(ex(1:2),crho,sigma,comega(:,:,k),eth_omega(:,:,k),0,1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
     call derivs_eth(ex(1:2),crho,sigma,U(:,:,k),eth_Ub(:,:,k),1,-1,   &
                     quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k),gR(:,:,k),gI(:,:,k))
  enddo

!!! The term * e^{-2beta} has been added so as to be consistent with HPN. Nigel
    !omega_u = - dble(eth_omega * conjg(U) + 0.5d0 * omega * eth_Ub * exp(-2*beta))

!!! - update .. I thought this may have been wrong so I removed the
!!! e^{-2beta} for testing. Yosef
!  omegarhs = - dreal(eth_omega * dconjg(U) + 0.5d0 * omega * eth_Ub)

  omegarhs = - 0.5d0*dreal(eth_Ub)

  gont = 0

  return

end function omega_rhs
!---------------------------------------------------------------------------------------------------------
subroutine drive_null_news(ex,crho,sigma,R,RJ,IJ,RU,IU,RTheta,ITheta,omega,beta, &
                         qlR1,qlR2,qlI1,qlI2, &
                         quR1,quR2,quI1,quI2, &
                         gR,gI,               &
                         dquR1,dquR2,dquI1,dquI2, &
                         bdquR1,bdquR2,bdquI1,bdquI2, &
                         dgR,dgI,bdgR,bdgI,RNews,INews,Rmin,sst)

  implicit none

  integer,intent(in) :: ex(3),sst
  real*8,intent(in) :: Rmin
  real*8,intent(in),dimension(ex(1))::crho
  real*8,intent(in),dimension(ex(2))::sigma
  real*8,intent(in),dimension(ex(3))::R
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: RJ,IJ,RU,IU,RTheta,ITheta,omega,beta
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: RNews,INews
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: qlR1,qlR2,qlI1,qlI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: quR1,quR2,quI1,quI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: gR,gI
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: dquR1,dquR2,dquI1,dquI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: bdquR1,bdquR2,bdquI1,bdquI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: dgR,dgI,bdgR,bdgI

  integer :: i,j,k
  double complex, dimension(ex(1),ex(2),ex(3)) :: CJ,U,J_u,J_l,J_l_u,News
#if 0
  call get_fake_Ju(ex,crho,sigma,R,RTheta,ITheta,                         &
                        quR1,quR2,quI1,quI2,qlR1,qlR2,qlI1,qlI2,          &
                        gR,gI,                                            &
                        dquR1,dquR2,dquI1,dquI2,                          &
                        bdquR1,bdquR2,bdquI1,bdquI2,                      &
                        dgR,dgI,bdgR,bdgI,dacos(-1.d0)/2,Rmin,sst)
#endif                        

  CJ = dcmplx(RJ,IJ)
  U = dcmplx(RU,IU)
  J_u = dcmplx(RTheta,ITheta)

  do j=1,ex(2)
  do i=1,ex(1)
     call cderivs_x(ex(3),R,CJ(i,j,:),J_l(i,j,:))
     call cderivs_x(ex(3),R,J_u(i,j,:),J_l_u(i,j,:))
     J_l(i,j,:) = -J_l(i,j,:)*Rmin*R**2
     J_l_u(i,j,:) = -J_l_u(i,j,:)*Rmin*R**2
  enddo
  enddo

#if 0  
if(sst == 0 .and. crho(1) < -dacos(-1.d0)/4 .and. sigma(1) < -dacos(-1.d0)/4)then  
  call get_exact_Jul(ex,crho,sigma,R,RNews,INews,                         &
                        quR1,quR2,quI1,quI2,qlR1,qlR2,qlI1,qlI2,          &
                        gR,gI,                                            &
                        dquR1,dquR2,dquI1,dquI2,                          &
                        bdquR1,bdquR2,bdquI1,bdquI2,                      &
                        dgR,dgI,bdgR,bdgI,dacos(-1.d0)/2,Rmin,sst)
write(*,*) J_u(ex(1)/2,ex(2)/2,ex(3)-1),J_u(ex(1)/2,ex(2)/2,ex(3))                        
write(*,*) RNews(ex(1)/2,ex(2)/2,ex(3)),INews(ex(1)/2,ex(2)/2,ex(3))     
write(*,*) J_l_u(ex(1)/2,ex(2)/2,ex(3))                   
write(*,*)dcmplx(RNews(ex(1)/2,ex(2)/2,ex(3)),INews(ex(1)/2,ex(2)/2,ex(3)))/J_l_u(ex(1)/2,ex(2)/2,ex(3))
endif
stop
#endif

  do k=1,ex(3)
    call get_null_news(ex(1:2),crho,sigma,CJ(:,:,k),U(:,:,k),J_u(:,:,k),J_l(:,:,k),J_l_u(:,:,k),omega(:,:,k),beta(:,:,k), &
                         qlR1(:,:,k),qlR2(:,:,k),qlI1(:,:,k),qlI2(:,:,k), &
                         quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k), &
                         gR(:,:,k),gI(:,:,k),               &
                         dquR1(:,:,k),dquR2(:,:,k),dquI1(:,:,k),dquI2(:,:,k), &
                         bdquR1(:,:,k),bdquR2(:,:,k),bdquI1(:,:,k),bdquI2(:,:,k), &
                         dgR(:,:,k),dgI(:,:,k),bdgR(:,:,k),bdgI(:,:,k),News(:,:,k))
  enddo

  RNews = dreal(News)
  INews = dimag(News)

#if 0        
if(sst ==0 .and. crho(1) < -dacos(-1.d0)/4 .and. sigma(1) < -dacos(-1.d0)/4)then  
  call get_exact_eth2omega(ex,crho,sigma,R,RNews,INews,                   &
                        quR1,quR2,quI1,quI2,qlR1,qlR2,qlI1,qlI2,          &
                        gR,gI,                                            &
                        dquR1,dquR2,dquI1,dquI2,                          &
                        bdquR1,bdquR2,bdquI1,bdquI2,                      &
                        dgR,dgI,bdgR,bdgI,dacos(-1.d0)/2,Rmin,sst)    
write(*,*) RNews(ex(1)/2,ex(2)/2,ex(3)),INews(ex(1)/2,ex(2)/2,ex(3))
endif
stop
#endif                        

#if 0
! check orthornormality
  RNews = RJ
  INews = IJ

  RNews = 0.5d0*dreal(J_l_u)
  INews = 0.5d0*dimag(J_l_u)
#endif

  call six2spher(ex,crho,sigma,R,RNews,INews,2,Rmin,sst)

  return

end subroutine drive_null_news
!---------------------------------------------------------------------------------------------------------
subroutine drive_null_news_diff(ex,crho,sigma,R,RJ,IJ,RU,IU,RTheta,ITheta,omega,beta, &
                         qlR1,qlR2,qlI1,qlI2, &
                         quR1,quR2,quI1,quI2, &
                         gR,gI,               &
                         dquR1,dquR2,dquI1,dquI2, &
                         bdquR1,bdquR2,bdquI1,bdquI2, &
                         dgR,dgI,bdgR,bdgI,RNews,INews,Rmin,sst,Time)

  implicit none

  integer,intent(in) :: ex(3),sst
  real*8,intent(in) :: Rmin,Time
  real*8,intent(in),dimension(ex(1))::crho
  real*8,intent(in),dimension(ex(2))::sigma
  real*8,intent(in),dimension(ex(3))::R
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: RJ,IJ,RU,IU,RTheta,ITheta,omega,beta
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: RNews,INews
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: qlR1,qlR2,qlI1,qlI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: quR1,quR2,quI1,quI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: gR,gI
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: dquR1,dquR2,dquI1,dquI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: bdquR1,bdquR2,bdquI1,bdquI2
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: dgR,dgI,bdgR,bdgI

  integer :: i,j,k
  double complex, dimension(ex(1),ex(2),ex(3)) :: CJ,U,J_u,J_l,J_l_u,News
#if 0
  call get_fake_Ju(ex,crho,sigma,R,RTheta,ITheta,                         &
                        quR1,quR2,quI1,quI2,qlR1,qlR2,qlI1,qlI2,          &
                        gR,gI,                                            &
                        dquR1,dquR2,dquI1,dquI2,                          &
                        bdquR1,bdquR2,bdquI1,bdquI2,                      &
                        dgR,dgI,bdgR,bdgI,dacos(-1.d0)/2,Rmin,sst)
#endif                        

  CJ = dcmplx(RJ,IJ)
  U = dcmplx(RU,IU)
  J_u = dcmplx(RTheta,ITheta)

  do j=1,ex(2)
  do i=1,ex(1)
     call cderivs_x(ex(3),R,CJ(i,j,:),J_l(i,j,:))
     call cderivs_x(ex(3),R,J_u(i,j,:),J_l_u(i,j,:))
     J_l(i,j,:) = -J_l(i,j,:)*Rmin*R**2
     J_l_u(i,j,:) = -J_l_u(i,j,:)*Rmin*R**2
  enddo
  enddo

#if 0  
if(sst == 0 .and. crho(1) < -dacos(-1.d0)/4 .and. sigma(1) < -dacos(-1.d0)/4)then  
  call get_exact_Jul(ex,crho,sigma,R,RNews,INews,                         &
                        quR1,quR2,quI1,quI2,qlR1,qlR2,qlI1,qlI2,          &
                        gR,gI,                                            &
                        dquR1,dquR2,dquI1,dquI2,                          &
                        bdquR1,bdquR2,bdquI1,bdquI2,                      &
                        dgR,dgI,bdgR,bdgI,dacos(-1.d0)/2,Rmin,sst)
write(*,*) J_u(ex(1)/2,ex(2)/2,ex(3)-1),J_u(ex(1)/2,ex(2)/2,ex(3))                        
write(*,*) RNews(ex(1)/2,ex(2)/2,ex(3)),INews(ex(1)/2,ex(2)/2,ex(3))     
write(*,*) J_l_u(ex(1)/2,ex(2)/2,ex(3))                   
write(*,*)dcmplx(RNews(ex(1)/2,ex(2)/2,ex(3)),INews(ex(1)/2,ex(2)/2,ex(3)))/J_l_u(ex(1)/2,ex(2)/2,ex(3))
endif
stop
#endif

  do k=1,ex(3)
    call get_null_news(ex(1:2),crho,sigma,CJ(:,:,k),U(:,:,k),J_u(:,:,k),J_l(:,:,k),J_l_u(:,:,k),omega(:,:,k),beta(:,:,k), &
                         qlR1(:,:,k),qlR2(:,:,k),qlI1(:,:,k),qlI2(:,:,k), &
                         quR1(:,:,k),quR2(:,:,k),quI1(:,:,k),quI2(:,:,k), &
                         gR(:,:,k),gI(:,:,k),               &
                         dquR1(:,:,k),dquR2(:,:,k),dquI1(:,:,k),dquI2(:,:,k), &
                         bdquR1(:,:,k),bdquR2(:,:,k),bdquI1(:,:,k),bdquI2(:,:,k), &
                         dgR(:,:,k),dgI(:,:,k),bdgR(:,:,k),bdgI(:,:,k),News(:,:,k))
  enddo

  call get_exact_news(ex,crho,sigma,R,RNews,INews,sst,Rmin,Time)

  RNews = dreal(News) - Rnews
  INews = dimag(News) - INews

!this part is nonsence  
  RNews(:,:,1:ex(3)-1) = 0.d0
  INews(:,:,1:ex(3)-1) = 0.d0
#if 0        
if(sst ==0 .and. crho(1) < -dacos(-1.d0)/4 .and. sigma(1) < -dacos(-1.d0)/4)then  
  call get_exact_eth2omega(ex,crho,sigma,R,RNews,INews,                   &
                        quR1,quR2,quI1,quI2,qlR1,qlR2,qlI1,qlI2,          &
                        gR,gI,                                            &
                        dquR1,dquR2,dquI1,dquI2,                          &
                        bdquR1,bdquR2,bdquI1,bdquI2,                      &
                        dgR,dgI,bdgR,bdgI,dacos(-1.d0)/2,Rmin,sst)    
write(*,*) RNews(ex(1)/2,ex(2)/2,ex(3)),INews(ex(1)/2,ex(2)/2,ex(3))
endif
stop
#endif                        

#if 0
! check orthornormality
  RNews = RJ
  INews = IJ

  RNews = 0.5d0*dreal(J_l_u)
  INews = 0.5d0*dimag(J_l_u)
#endif

  call six2spher(ex,crho,sigma,R,RNews,INews,2,Rmin,sst)

  return

end subroutine drive_null_news_diff
!------------------------------------------------------------------------------------------------------------
subroutine get_null_news(ex,crho,sigma,J,U,J_u,J_l,J_l_u,omega,beta, &
                         qlR1,qlR2,qlI1,qlI2, &
                         quR1,quR2,quI1,quI2, &
                         gR,gI,               &
                         dquR1,dquR2,dquI1,dquI2, &
                         bdquR1,bdquR2,bdquI1,bdquI2, &
                         dgR,dgI,bdgR,bdgI,News)

implicit none

integer,intent(in) :: ex(2)
real*8,intent(in),dimension(ex(1))::crho
real*8,intent(in),dimension(ex(2))::sigma
double complex,dimension(ex(1),ex(2)),intent(in) :: J,U
double complex,dimension(ex(1),ex(2)),intent(in) :: J_u,J_l,J_l_u
real*8,dimension(ex(1),ex(2)),intent(in) :: omega,beta
real*8,dimension(ex(1),ex(2)),intent(in) :: qlR1,qlR2,qlI1,qlI2
real*8,dimension(ex(1),ex(2)),intent(in) :: quR1,quR2,quI1,quI2
real*8,dimension(ex(1),ex(2)),intent(in) :: gR,gI
real*8,dimension(ex(1),ex(2)),intent(in) :: dquR1,dquR2,dquI1,dquI2
real*8,dimension(ex(1),ex(2)),intent(in) :: bdquR1,bdquR2,bdquI1,bdquI2
real*8,dimension(ex(1),ex(2)),intent(in) :: dgR,dgI,bdgR,bdgI
double complex,dimension(ex(1),ex(2)),intent(out) :: News

! local variables
real*8,dimension(ex(1),ex(2)) :: K,K_u,K_l,K_l_u
real*8,dimension(ex(1),ex(2)) :: a
double complex,dimension(ex(1),ex(2)) :: Comega,Cbeta
double complex,dimension(ex(1),ex(2)) :: Jb,Ub
double complex,dimension(ex(1),ex(2)) :: eth_a,eth2_a,eth_ethb_a
double complex,dimension(ex(1),ex(2)) :: s1,s2,s3,s4,s5
double complex,dimension(ex(1),ex(2)) :: eth_U,ethb_U,eth_J,ethb_J
double complex,dimension(ex(1),ex(2)) :: eth_J_l,ethb_J_l,eth_K_l,eth_K
double complex,dimension(ex(1),ex(2)) :: eth_omega,eth_beta
double complex,dimension(ex(1),ex(2)) :: eth2_omega,eth2_beta
double complex,dimension(ex(1),ex(2)) :: eth_ethb_omega,eth_ethb_beta

    Comega = dcmplx(omega,0.d0)
    Cbeta = dcmplx(beta,0.d0)
    call derivs_eth(ex,crho,sigma,Comega,eth_omega,0,1,quR1,quR2,quI1,quI2,gR,gI)
    call derivs_eth(ex,crho,sigma,Cbeta,eth_beta,0,1,quR1,quR2,quI1,quI2,gR,gI)
    call dderivs_eth(ex,crho,sigma,Comega,eth2_omega,0,1,1,   &
                         quR1,quR2,quI1,quI2,gR,gI, &
                         dquR1,dquR2,dquI1,dquI2,   &
                         bdquR1,bdquR2,bdquI1,bdquI2,   &
                         dgR,dgI,bdgR,bdgI)
    call dderivs_eth(ex,crho,sigma,Cbeta,eth2_beta,0,1,1,   &
                         quR1,quR2,quI1,quI2,gR,gI, &
                         dquR1,dquR2,dquI1,dquI2,   &
                         bdquR1,bdquR2,bdquI1,bdquI2,   &
                         dgR,dgI,bdgR,bdgI)
    call dderivs_eth(ex,crho,sigma,Comega,eth_ethb_omega,0,-1,1,   &
                         quR1,quR2,quI1,quI2,gR,gI, &
                         dquR1,dquR2,dquI1,dquI2,   &
                         bdquR1,bdquR2,bdquI1,bdquI2,   &
                         dgR,dgI,bdgR,bdgI)
    call dderivs_eth(ex,crho,sigma,Cbeta,eth_ethb_beta,0,-1,1,   &
                         quR1,quR2,quI1,quI2,gR,gI, &
                         dquR1,dquR2,dquI1,dquI2,   &
                         bdquR1,bdquR2,bdquI1,bdquI2,   &
                         dgR,dgI,bdgR,bdgI)
    call derivs_eth(ex,crho,sigma,U,eth_U,1,1,quR1,quR2,quI1,quI2,gR,gI)
    call derivs_eth(ex,crho,sigma,U,ethb_U,1,-1,quR1,quR2,quI1,quI2,gR,gI)
    call derivs_eth(ex,crho,sigma,J,eth_J,2,1,quR1,quR2,quI1,quI2,gR,gI)
    call derivs_eth(ex,crho,sigma,J,ethb_J,2,-1,quR1,quR2,quI1,quI2,gR,gI)
    call derivs_eth(ex,crho,sigma,J_l,eth_J_l,2,1,quR1,quR2,quI1,quI2,gR,gI)
    call derivs_eth(ex,crho,sigma,J_l,ethb_J_l,2,-1,quR1,quR2,quI1,quI2,gR,gI)

    Jb = dconjg(J)
    Ub = dconjg(U)
    K = dsqrt(1.0d0 + cdabs(J)**2)
! temp storage
    Comega=dcmplx(K,0.d0)
    call derivs_eth(ex,crho,sigma,Comega,eth_K,0,1,quR1,quR2,quI1,quI2,gR,gI)

    K_u = dreal( J_u * Jb ) / K
    K_l = dreal( J_l * Jb ) / K
! temp storage
    Comega=dcmplx(K_l,0.d0)
    call derivs_eth(ex,crho,sigma,Comega,eth_K_l,0,1,quR1,quR2,quI1,quI2,gR,gI)
    K_l_u = dreal( J_u * dconjg(J_l) + J_l_u * Jb )/ K - K_l * K_u / K 

    a = omega * dexp(2.0d0 * beta)

    eth_a = dexp(2.0d0 * beta) * ( eth_omega + 2.0d0 * omega * eth_beta )

    eth2_a = dexp(2.0d0 * beta) * ( 4.0d0 * eth_beta * eth_omega &
         + 4.0d0 * omega * eth_beta**2 &
         + eth2_omega + 2.0d0 * omega * eth2_beta )

    eth_ethb_a = dexp(2.0d0 * beta) * ( 4.0d0 * dreal(eth_beta * dconjg(eth_omega)) &
         + 4.0d0 * omega * eth_beta * dconjg(eth_beta) &
         + eth_ethb_omega + 2.0d0 * omega * eth_ethb_beta )

    s1 = ( -2.0d0 * K_l_u * J * (K + 1.0d0) + J_l_u * (K + 1.0d0)**2 &
         + dconjg(J_l_u) * J**2 ) / (K + 1.0d0)

    s2 = 0.5d0 / ( K + 1.0d0) * ( &
         (K + 1.0d0)* (eth_J_l *Ub * (K+1.0d0) - 2.0d0* eth_K_l * J *Ub ) &
         + eth_U * (K+1.0d0)* ( -2.0d0 * J * dconjg(J_l) + K_l * 2.0d0 * (K+1.0d0) ) &
         + dconjg(ethb_U) * (K+1.0d0) * ( -2.0d0* J * K_l + J_l * 2.0d0 * (K+1.0d0) ) &
         + ethb_J_l * U * (K+1.0d0)**2 - dconjg(eth_K_l) * 2.0d0 * U * J * (K+1.0d0) &
         + ethb_U * 2.0d0 * J * ( J * dconjg(J_l) - (K+1.0d0) * K_l) &
         + J**2 * ( U * dconjg(eth_J_l) + dconjg(ethb_J_l * U) ) &
         + J * 2.0d0 * dconjg(eth_U) * ( J * K_l - J_l * (K+1.0d0) ) )

    s3 = ( J_l * (K + 1.0d0)**2 -2.0d0 * K_l * J * (K + 1.0d0) &
         + dconjg(J_l) * J**2) / (K + 1.0d0)

    s4 = 0.5d0 / ( K + 1.0d0) * ( eth_a * eth_omega * (K + 1.0d0)**2 &
         - (K+1.0d0) * J * 2.0d0* dreal( eth_a * dconjg(eth_omega) ) &
         + J**2 * dconjg(eth_a * eth_omega) )

    s5 = 0.25d0 / ( K + 1.0d0) * ( 2.0d0 * eth2_a * (K+1.0d0)**2 &
         + 2.0d0 * J**2 * dconjg(eth2_a) &
         - 4.0d0 * eth_ethb_a * J * (K+1.0d0) &
         + Jb * eth_a * eth_J* (K+1.0d0)**2 &
         + J * eth_a * dconjg(ethb_J) * (K+1.0d0)**2 &
         - eth_a * eth_K * 2.0d0 * (K+1.0d0) * ( J*Jb + (K+1.0d0) ) & 
         + eth_a * ethb_J * (K+1.0d0) * ( -J*Jb + (K+1.0d0) ) &
         - J**2 * eth_a * dconjg(eth_J) * K &
         +  J**2 * Jb * 2.0d0* eth_a * dconjg(eth_K) &
         - dconjg(eth_a) * eth_J * (K+1.0d0) * ( J*Jb + K+1.0d0 ) &
         - dconjg(ethb_J) * dconjg(eth_a) * J**2 * ( K + 2.0d0) &
         + J * 2.0d0 * (K+1.0d0)**2 * eth_K * dconjg(eth_a)  &
         + J**2 * Jb * ethb_J * dconjg(eth_a) &
         + J**3 * dconjg(eth_a * eth_J) &
         - 2.0d0* J**2 *K*dconjg(eth_K * eth_a) )

    !   News = 0.25d0 * ( s1 + s2 + 0.5d0 * dble(ethb_U) * s3 &
    !   - 4.0d0 * s4 / omega**2 + 2.0d0 * s5 / omega ) / ( omega**2 * exp(2.0d0 * beta) )

    ! change sign of s3 to compensate for a bug in Eqs. 30, 37, and 38 of
    ! HPN
#if 1 
    News = 0.25d0 * ( s1 + s2 - 0.5d0 * dreal(ethb_U) * s3 &
         - 4.0d0 * s4 / omega**2 + 2.0d0 * s5 / omega ) / ( omega**2 * dexp(2.0d0 * beta) )
#else   
#if 0
if(crho(1) < -dacos(-1.d0)/4 .and. sigma(1) < -dacos(-1.d0)/4)then 
write(*,*) eth2_omega(ex(1)/2,ex(2)/2)                   
endif
#endif
     News = 0.5d0*J_l_u+eth2_beta+0.5d0*eth2_omega  ! if given omega error is about 6e-9
!    News = 0.5d0*J_l_u+eth2_beta-1.5d0*J    ! error is about 6e-9
#endif
    return
        
end subroutine get_null_news        
!--------------------------------------------------------------------------------------------------
! change spin weighted function from 6 patches to spherical coordinate
  subroutine six2spher(ex,crho,sigma,R,RU,IU,spin,Rmin,sst)

  implicit none

!~~~~~~% Input parameters:
  integer,intent(in) :: ex(3),sst,spin
  real*8,intent(in) :: Rmin
  double precision,intent(in),dimension(ex(1))::crho
  double precision,intent(in),dimension(ex(2))::sigma
  double precision,intent(in),dimension(ex(3))::R
  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout) :: RU,IU

integer :: i,j,k
real*8 ::x,y,z,hgr,gt,gp,tgrho,tgsigma,tc,ts,rf
double complex :: II,swtf,ff

  II = dcmplx(0.d0,1.d0)
  hgr = 1.d0

  do i=1,ex(1)
  do j=1,ex(2)
  do k=1,ex(3)
!    hgr = R(k)*Rmin/(1.d0-R(k))  R is not invovled indeed, to avoid NaN, we set
!    it to 1 above
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
      write(*,*) "six2spher: not recognized sst = ",sst
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

    ff=dcmplx(RU(i,j,k),IU(i,j,k))/swtf**spin

    RU(i,j,k) = dreal(ff)
    IU(i,j,k) = dimag(ff)
  enddo
  enddo
  enddo

  return

  end subroutine six2spher
!-------------------------------------------------------------  
! Linear wave given in Eq.(27) of CQG 22, 2393 (2005)
!-------------------------------------------------------------
subroutine get_exact_omega(ex,crho,sigma,R,omega,sst,Rmin,T)
implicit none
! argument variables
integer, intent(in ):: ex(1:3),sst
real*8,intent(in) :: Rmin,T
double precision,intent(in),dimension(ex(1))::crho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::R
real*8,dimension(ex(1),ex(2),ex(3)),intent(out)::omega

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
      write(*,*) "get_exact_omega: not recognized sst = ",sst
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
    gr = dreal(Jr*cdexp(II*nu*T))
    Jr = Yslm(0,2,m,gt,gp)
    omega(i,j,k) = 1.d0-2.d0*(2+1)/2.d0*gr*dreal(Jr)

  enddo
  enddo
  enddo

return

end subroutine get_exact_omega
!-------------------------------------------------------------  
! Linear wave given in Eq.(16) of CQG 24S327
!-------------------------------------------------------------
subroutine get_exact_news(ex,crho,sigma,R,RNews,INews,sst,Rmin,Time)
implicit none
! argument variables
integer, intent(in ):: ex(1:3),sst
real*8,intent(in) :: Rmin,Time
double precision,intent(in),dimension(ex(1))::crho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::R
real*8,dimension(ex(1),ex(2),ex(3)),intent(out)::RNews,INews

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

    Jr = II*nu**3*C2/dsqrt(2.4d1)
    gr = dreal(Jr)
    Jr = Yslm(2,2,m,gt,gp)
    ff = gr*Jr*swtf**2
    RNews(i,j,k) = dreal(ff)
    INews(i,j,k) = dimag(ff)

  enddo
  enddo
  enddo

return

end subroutine get_exact_news
