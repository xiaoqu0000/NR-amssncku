

#include "macrodef.fh"

!---------------------------------------------------------------------------------
! fill symmetric boundary buffer points
!---------------------------------------------------------------------------------
subroutine fill_symmetric_boundarybuffer2(ex,crho,sigma,R,drho,dsigma,   &
                                          var,Symmetry,sst,AoS)

  implicit none

!~~~~~~% Input parameters:
  integer,intent(in),dimension(3) :: ex
  integer,intent(in) :: Symmetry,sst
  real*8,dimension(3) :: AoS
  real*8,intent(in),dimension(ex(1))::crho
  real*8,intent(in),dimension(ex(2))::sigma
  real*8,intent(in),dimension(ex(3))::R
  real*8,intent(in) :: drho,dsigma
  real*8,intent(inout),dimension(ex(1),ex(2),ex(3)) :: var
    
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
      var(i,j,k) = AoS(2)*var(i,t,k)               
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
      var(i,j,k) = AoS(2)*var(i,t,k)    
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
      var(i,j,k) = AoS(2)*var(i,t,k)    
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
      var(i,j,k) = AoS(1)*var(t,j,k)                 
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
      var(i,j,k) = AoS(2)*var(i,t,k)                 
       enddo
       enddo
       enddo
     endif
    endif
  end select

  return

  end subroutine fill_symmetric_boundarybuffer2
!---------------------------------------------------------------------------------
!!!! using r^2g_AB instead of g_AB
!!!! using r^2g_0A instead of g_0A
!!!! using r^2g_00 instead of g_00
!!!! using x in the metric form directly instead of r
!---------------------------------------------------------------------------------
! this R is indeed x
function NullEvol_g01(ex,crho,sigma,R, &
                      g22,g23,g33,g01,Rmin)    result(gont)   
  implicit none
  integer,intent(in ):: ex(1:3)
  real*8,intent(in ):: Rmin
  real*8,intent(in),dimension(ex(1))::crho
  real*8,intent(in),dimension(ex(2))::sigma
  real*8,intent(in),dimension(ex(3))::R
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: g01
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: g22,g23,g33
!  gont = 0: success; gont = 1: something wrong
  integer::gont
  real*8 :: dR

  real*8, dimension(ex(3)):: dg22,dg23,dg33,ddg22,ddg23,ddg33
  real*8, dimension(ex(3)):: Hg22,Hg23,Hg33
  real*8, dimension(ex(3)):: Hdg22,Hdg23,Hdg33,Hddg22,Hddg23,Hddg33
  real*8 :: g010,g011,g01h,rhs
  integer :: i,j,k,RK4

!!! sanity check
  dR = sum(g22)+sum(g23)+sum(g33)+sum(g01)
  if(dR.ne.dR) then
     if(sum(g22).ne.sum(g22))write(*,*)"NullEvol_g01: find NaN in g22"
     if(sum(g23).ne.sum(g23))write(*,*)"NullEvol_g01: find NaN in g23"
     if(sum(g33).ne.sum(g33))write(*,*)"NullEvol_g01: find NaN in g33"
     if(sum(g01).ne.sum(g01))write(*,*)"NullEvol_g01: find NaN in g01"
     gont = 1
     return
  endif

  dR = R(2) - R(1)

  do j=1,ex(2)
  do i=1,ex(1)
     g010 = g01(i,j,1)

     call rderivs_x(ex(3),R,g22(i,j,:),dg22)
     call rderivs_x(ex(3),R,g23(i,j,:),dg23)
     call rderivs_x(ex(3),R,g33(i,j,:),dg33)
     call rdderivs_x(ex(3),R,g22(i,j,:),ddg22)
     call rdderivs_x(ex(3),R,g23(i,j,:),ddg23)
     call rdderivs_x(ex(3),R,g33(i,j,:),ddg33)

     call rget_half_x(ex(3),g22(i,j,:),Hg22)
     call rget_half_x(ex(3),g23(i,j,:),Hg23)
     call rget_half_x(ex(3),g33(i,j,:),Hg33)

     call rget_half_x(ex(3),dg22,Hdg22)
     call rget_half_x(ex(3),dg23,Hdg23)
     call rget_half_x(ex(3),dg33,Hdg33)

     call rget_half_x(ex(3),ddg22,Hddg22)
     call rget_half_x(ex(3),ddg23,Hddg23)
     call rget_half_x(ex(3),ddg33,Hddg33)

     do k=1,ex(3)-2
        RK4 = 0
        call get_g01_rhs(R(k),g22(i,j,k),g23(i,j,k),g33(i,j,k),dg22(k),dg23(k), &
                    dg33(k),ddg22(k),ddg23(k),ddg33(k),g01(i,j,k),rhs)
        call rungekutta4_scalar(dR,g010,g01h,rhs,RK4)

        RK4 = 1
        call get_g01_rhs(R(k)+dR/2,Hg22(k),Hg23(k),Hg33(k),Hdg22(k),Hdg23(k), &
                    Hdg33(k),Hddg22(k),Hddg23(k),Hddg33(k),g01h,g011)
        call rungekutta4_scalar(dR,g010,g011,rhs,RK4)
        call rswap(g01h,g011)

        RK4 = 2
        call get_g01_rhs(R(k)+dR/2,Hg22(k),Hg23(k),Hg33(k),Hdg22(k),Hdg23(k), &
                    Hdg33(k),Hddg22(k),Hddg23(k),Hddg33(k),g01h,g011)
        call rungekutta4_scalar(dR,g010,g011,rhs,RK4)
        call rswap(g01h,g011)

        RK4 = 3
        call get_g01_rhs(R(k+1),g22(i,j,k+1),g23(i,j,k+1),g33(i,j,k+1),dg22(k+1),dg23(k+1), &
                    dg33(k+1),ddg22(k+1),ddg23(k+1),ddg33(k+1),g01h,g011)
        call rungekutta4_scalar(dR,g010,g011,rhs,RK4)
        call rswap(g010,g011)

        g01(i,j,k+1) = g010
     enddo
! closing step
     k = ex(3)-1
     call get_g01_rhs(R(k),g22(i,j,k),g23(i,j,k),g33(i,j,k),dg22(k),dg23(k), &
                 dg33(k),ddg22(k),ddg23(k),ddg33(k),g01(i,j,k),rhs)
     g01(i,j,k+1) = g01(i,j,k) + rhs*dR

  enddo
  enddo

  gont = 0

  return

end function NullEvol_g01
!------------------------------------------------------------------------------
! this R is indeed x
function NullEvol_pg0a(ex,crho,sigma,R, &
                      g22,g23,g33,g01,p02,p03,g02,g03,Rmin)    result(gont)   
  implicit none
  integer,intent(in ):: ex(1:3)
  real*8,intent(in ):: Rmin
  real*8,intent(in),dimension(ex(1))::crho
  real*8,intent(in),dimension(ex(2))::sigma
  real*8,intent(in),dimension(ex(3))::R
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: p02,p03,g02,g03
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: g22,g23,g33,g01
!  gont = 0: success; gont = 1: something wrong
  integer::gont
  real*8 :: dR

  real*8, dimension(ex(3)) :: Hg01
  real*8, dimension(ex(3)) :: Hg22,Hg23,Hg33

  real*8, dimension(ex(3)) :: dg01,dg02,dg03
  real*8, dimension(ex(3)) :: dgx01,dgx22,dgx23,dgx33
  real*8, dimension(ex(3)) :: dgy01,dgy22,dgy23,dgy33
  real*8, dimension(ex(3)) :: ddgxr01,ddgxr22,ddgxr23,ddgxr33
  real*8, dimension(ex(3)) :: ddgyr01,ddgyr22,ddgyr23,ddgyr33
  real*8, dimension(ex(3)) :: dg22,dg23,dg33,ddg22,ddg23,ddg33
  real*8, dimension(ex(3)) :: Hdg01,Hdg02,Hdg03
  real*8, dimension(ex(3)) :: Hdgx01,Hdgx22,Hdgx23,Hdgx33
  real*8, dimension(ex(3)) :: Hdgy01,Hdgy22,Hdgy23,Hdgy33
  real*8, dimension(ex(3)) :: Hddgxr01,Hddgxr22,Hddgxr23,Hddgxr33
  real*8, dimension(ex(3)) :: Hddgyr01,Hddgyr22,Hddgyr23,Hddgyr33
  real*8, dimension(ex(3)) :: Hdg22,Hdg23,Hdg33,Hddg22,Hddg23,Hddg33

  real*8 :: p020,p021,p02h,p02_rhs
  real*8 :: p030,p031,p03h,p03_rhs
  real*8 :: g020,g021,g02h,g02_rhs
  real*8 :: g030,g031,g03h,g03_rhs
  integer :: i,j,k,RK4

!!! sanity check
  dR = sum(g22)+sum(g23)+sum(g33)+sum(g01) &
      +sum(p02)+sum(p03)+sum(g02)+sum(g03)
  if(dR.ne.dR) then
     if(sum(g22).ne.sum(g22))write(*,*)"NullEvol_pg0a: find NaN in g22"
     if(sum(g23).ne.sum(g23))write(*,*)"NullEvol_pg0a: find NaN in g23"
     if(sum(g33).ne.sum(g33))write(*,*)"NullEvol_pg0a: find NaN in g33"
     if(sum(g01).ne.sum(g01))write(*,*)"NullEvol_pg0a: find NaN in g01"
     if(sum(p02).ne.sum(p02))write(*,*)"NullEvol_pg0a: find NaN in p02"
     if(sum(p03).ne.sum(p03))write(*,*)"NullEvol_pg0a: find NaN in p03"
     if(sum(g02).ne.sum(g02))write(*,*)"NullEvol_pg0a: find NaN in g02"
     if(sum(g03).ne.sum(g03))write(*,*)"NullEvol_pg0a: find NaN in g03"
     gont = 1
     return
  endif

  dR = R(2) - R(1)

  do j=1,ex(2)
  do i=1,ex(1)

     call rderivs_x(ex(3),R,g01(i,j,:),dg01)
     dg02 = p02(i,j,:)
     dg03 = p03(i,j,:)
     call rderivs_x(ex(3),R,g22(i,j,:),dg22)
     call rderivs_x(ex(3),R,g23(i,j,:),dg23)
     call rderivs_x(ex(3),R,g33(i,j,:),dg33)
     call rdderivs_x(ex(3),R,g22(i,j,:),ddg22)
     call rdderivs_x(ex(3),R,g23(i,j,:),ddg23)
     call rdderivs_x(ex(3),R,g33(i,j,:),ddg33)

     do k=1,ex(3)
       call rderivs_x_point(ex(1),crho,g01(:,j,k),dgx01(k),i)
       call rderivs_x_point(ex(1),crho,g22(:,j,k),dgx22(k),i)
       call rderivs_x_point(ex(1),crho,g23(:,j,k),dgx23(k),i)
       call rderivs_x_point(ex(1),crho,g33(:,j,k),dgx33(k),i)

       call rderivs_x_point(ex(2),sigma,g01(i,:,k),dgy01(k),j)
       call rderivs_x_point(ex(2),sigma,g22(i,:,k),dgy22(k),j)
       call rderivs_x_point(ex(2),sigma,g23(i,:,k),dgy23(k),j)
       call rderivs_x_point(ex(2),sigma,g33(i,:,k),dgy33(k),j)

       call rdderivs_xy_point(ex(1),ex(3),crho,R,g01(:,j,:),ddgxr01(k),i,k)
       call rdderivs_xy_point(ex(1),ex(3),crho,R,g22(:,j,:),ddgxr22(k),i,k)
       call rdderivs_xy_point(ex(1),ex(3),crho,R,g23(:,j,:),ddgxr23(k),i,k)
       call rdderivs_xy_point(ex(1),ex(3),crho,R,g33(:,j,:),ddgxr33(k),i,k)

       call rdderivs_xy_point(ex(2),ex(3),sigma,R,g01(i,:,:),ddgyr01(k),j,k)
       call rdderivs_xy_point(ex(2),ex(3),sigma,R,g22(i,:,:),ddgyr22(k),j,k)
       call rdderivs_xy_point(ex(2),ex(3),sigma,R,g23(i,:,:),ddgyr23(k),j,k)
       call rdderivs_xy_point(ex(2),ex(3),sigma,R,g33(i,:,:),ddgyr33(k),j,k)
     enddo

     call rget_half_x(ex(3),g01(i,j,:),Hg01)
     call rget_half_x(ex(3),g22(i,j,:),Hg22)
     call rget_half_x(ex(3),g23(i,j,:),Hg23)
     call rget_half_x(ex(3),g33(i,j,:),Hg33)

     call rget_half_x(ex(3),dg01,Hdg01)
     call rget_half_x(ex(3),dg02,Hdg02)
     call rget_half_x(ex(3),dg03,Hdg03)

     call rget_half_x(ex(3),dgx01,Hdgx01)
     call rget_half_x(ex(3),dgy01,Hdgy01)

     call rget_half_x(ex(3),dgx22,Hdgx22)
     call rget_half_x(ex(3),dgx23,Hdgx23)
     call rget_half_x(ex(3),dgx33,Hdgx33)
     call rget_half_x(ex(3),dgy22,Hdgy22)
     call rget_half_x(ex(3),dgy23,Hdgy23)
     call rget_half_x(ex(3),dgy33,Hdgy33)

     call rget_half_x(ex(3),ddgxr01,Hddgxr01)
     call rget_half_x(ex(3),ddgyr01,Hddgyr01)

     call rget_half_x(ex(3),ddgxr22,Hddgxr22)
     call rget_half_x(ex(3),ddgxr23,Hddgxr23)
     call rget_half_x(ex(3),ddgxr33,Hddgxr33)
     call rget_half_x(ex(3),ddgyr22,Hddgyr22)
     call rget_half_x(ex(3),ddgyr23,Hddgyr23)
     call rget_half_x(ex(3),ddgyr33,Hddgyr33)

     call rget_half_x(ex(3),dg22,Hdg22)
     call rget_half_x(ex(3),dg23,Hdg23)
     call rget_half_x(ex(3),dg33,Hdg33)
     call rget_half_x(ex(3),ddg22,Hddg22)
     call rget_half_x(ex(3),ddg23,Hddg23)
     call rget_half_x(ex(3),ddg33,Hddg33)

#if 0
     g020 = g02(i,j,1)
     g030 = g03(i,j,1)
     p020 = p02(i,j,1)
     p030 = p03(i,j,1)

  do k=1,ex(3)-2
        RK4 = 0
        call pg0a_rhs(Rmin,R(k),p020,p030,g020,g030,g22(i,j,k),g23(i,j,k),g33(i,j,k),dg22(k),dg23(k), &
                     dg33(k),ddg22(k),ddg23(k),ddg33(k),g01(i,j,k), &
                     dg01(k),dg02(k),dg03(k),         &
                    dgx01(k),dgx22(k),dgx23(k),dgx33(k), &
                    dgy01(k),dgy22(k),dgy23(k),dgy33(k), &
                    ddgxr01(k),ddgxr22(k),ddgxr23(k),ddgxr33(k), &
                    ddgyr01(k),ddgyr22(k),ddgyr23(k),ddgyr33(k), &
                    g02_rhs,g03_rhs,p02_rhs,p03_rhs)
        call rungekutta4_scalar(dR,g020,g02h,g02_rhs,RK4)
        call rungekutta4_scalar(dR,g030,g03h,g03_rhs,RK4)
        call rungekutta4_scalar(dR,p020,p02h,p02_rhs,RK4)
        call rungekutta4_scalar(dR,p030,p03h,p03_rhs,RK4)

        RK4 = 1
        call pg0a_rhs(Rmin,R(k)+dR/2,p02h,p03h,g02h,g03h,Hg22(k),Hg23(k),Hg33(k),Hdg22(k),Hdg23(k), &
                     Hdg33(k),Hddg22(k),Hddg23(k),Hddg33(k),Hg01(k), &
                     Hdg01(k),Hdg02(k),Hdg03(k),         &
                    Hdgx01(k),Hdgx22(k),Hdgx23(k),Hdgx33(k), &
                    Hdgy01(k),Hdgy22(k),Hdgy23(k),Hdgy33(k), &
                    Hddgxr01(k),Hddgxr22(k),Hddgxr23(k),Hddgxr33(k), &
                    Hddgyr01(k),Hddgyr22(k),Hddgyr23(k),Hddgyr33(k), &
                    g021,g031,p021,p031)
        call rungekutta4_scalar(dR,g020,g021,g02_rhs,RK4)
        call rungekutta4_scalar(dR,g030,g031,g03_rhs,RK4)
        call rungekutta4_scalar(dR,p020,p021,p02_rhs,RK4)
        call rungekutta4_scalar(dR,p030,p031,p03_rhs,RK4)
        call rswap(g02h,g021)
        call rswap(g03h,g031)
        call rswap(p02h,p021)
        call rswap(p03h,p031)

        RK4 = 2
        call pg0a_rhs(Rmin,R(k)+dR/2,p02h,p03h,g02h,g03h,Hg22(k),Hg23(k),Hg33(k),Hdg22(k),Hdg23(k), &
                     Hdg33(k),Hddg22(k),Hddg23(k),Hddg33(k),Hg01(k), &
                     Hdg01(k),Hdg02(k),Hdg03(k),         &
                    Hdgx01(k),Hdgx22(k),Hdgx23(k),Hdgx33(k), &
                    Hdgy01(k),Hdgy22(k),Hdgy23(k),Hdgy33(k), &
                    Hddgxr01(k),Hddgxr22(k),Hddgxr23(k),Hddgxr33(k), &
                    Hddgyr01(k),Hddgyr22(k),Hddgyr23(k),Hddgyr33(k), &
                    g021,g031,p021,p031)
        call rungekutta4_scalar(dR,g020,g021,g02_rhs,RK4)
        call rungekutta4_scalar(dR,g030,g031,g03_rhs,RK4)
        call rungekutta4_scalar(dR,p020,p021,p02_rhs,RK4)
        call rungekutta4_scalar(dR,p030,p031,p03_rhs,RK4)
        call rswap(g02h,g021)
        call rswap(g03h,g031)
        call rswap(p02h,p021)
        call rswap(p03h,p031)

        RK4 = 3
        call pg0a_rhs(Rmin,R(k+1),p02h,p03h,g02h,g03h,Hg22(k+1),Hg23(k+1),Hg33(k+1),Hdg22(k+1),Hdg23(k+1), &
                     Hdg33(k+1),Hddg22(k+1),Hddg23(k+1),Hddg33(k+1),Hg01(k+1), &
                     Hdg01(k+1),Hdg02(k+1),Hdg03(k+1),         &
                    Hdgx01(k+1),Hdgx22(k+1),Hdgx23(k+1),Hdgx33(k+1), &
                    Hdgy01(k+1),Hdgy22(k+1),Hdgy23(k+1),Hdgy33(k+1), &
                    Hddgxr01(k+1),Hddgxr22(k+1),Hddgxr23(k+1),Hddgxr33(k+1), &
                    Hddgyr01(k+1),Hddgyr22(k+1),Hddgyr23(k+1),Hddgyr33(k+1), &
                    g021,g031,p021,p031)
        call rungekutta4_scalar(dR,g020,g021,g02_rhs,RK4)
        call rungekutta4_scalar(dR,g030,g031,g03_rhs,RK4)
        call rungekutta4_scalar(dR,p020,p021,p02_rhs,RK4)
        call rungekutta4_scalar(dR,p030,p031,p03_rhs,RK4)
        call rswap(g020,g021)
        call rswap(g030,g031)
        call rswap(p020,p021)
        call rswap(p030,p031)

        g02(i,j,k+1) = g020
        g03(i,j,k+1) = g030
        p02(i,j,k+1) = p020
        p03(i,j,k+1) = p030

  enddo
       k=ex(3)-1
!    closing step
        call pg0a_rhs(Rmin,R(k),p020,p030,g020,g030,g22(i,j,k),g23(i,j,k),g33(i,j,k),dg22(k),dg23(k), &
                     dg33(k),ddg22(k),ddg23(k),ddg33(k),g01(i,j,k), &
                     dg01(k),dg02(k),dg03(k),         &
                    dgx01(k),dgx22(k),dgx23(k),dgx33(k), &
                    dgy01(k),dgy22(k),dgy23(k),dgy33(k), &
                    ddgxr01(k),ddgxr22(k),ddgxr23(k),ddgxr33(k), &
                    ddgyr01(k),ddgyr22(k),ddgyr23(k),ddgyr33(k), &
                    g02_rhs,g03_rhs,p02_rhs,p03_rhs)
     g02(i,j,k+1) = g02(i,j,k) + g02_rhs*dR
     g03(i,j,k+1) = g03(i,j,k) + g03_rhs*dR
     p02(i,j,k+1) = p02(i,j,k) + p02_rhs*dR
     p03(i,j,k+1) = p03(i,j,k) + p03_rhs*dR
#endif

  enddo
  enddo

  gont = 0

  return

end function NullEvol_pg0a
!------------------------------------------------------------------------------
! this R is indeed x
function NullEvol_Theta2(ex,crho,sigma,R, &
                      g22,g23,g33,g00,g01,g02,g03,p02,p03, &
                      Theta22,Theta23,Theta33,Rmin)    result(gont)   
  implicit none
  integer,intent(in ):: ex(1:3)
  real*8,intent(in ):: Rmin
  real*8,intent(in),dimension(ex(1))::crho
  real*8,intent(in),dimension(ex(2))::sigma
  real*8,intent(in),dimension(ex(3))::R
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: g00
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: g02,g03,p02,p03
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in   ) :: g22,g23,g33,g01
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Theta22,Theta23,Theta33
!  gont = 0: success; gont = 1: something wrong
  integer::gont
  real*8 :: dR

  real*8,dimension(ex(3)) :: dg22,dg23,dg33,ddg22,ddg23,ddg33                   
  real*8,dimension(ex(3)) :: dg00,dg01,dg02,dg03
  real*8,dimension(ex(3)) :: dgx01,dgx02,dgx03
  real*8,dimension(ex(3)) :: dgy01,dgy02,dgy03
  real*8,dimension(ex(3)) :: dgx22,dgx23,dgx33
  real*8,dimension(ex(3)) :: dgy22,dgy23,dgy33
  real*8,dimension(ex(3)) :: ddgxx01,ddgxx33,ddgyy01,ddgyy22,ddgxy23
  real*8,dimension(ex(3)) :: ddgxy01,ddgxr02,ddgxr03,ddgyr02,ddgyr03
  real*8,dimension(ex(3)) :: ddgxr22,ddgxr23,ddgxr33,ddgyr22,ddgyr23,ddgyr33

  real*8,dimension(ex(3)) :: Hdg22,Hdg23,Hdg33,Hddg22,Hddg23,Hddg33                
  real*8,dimension(ex(3)) :: Hdg00,Hdg01,Hdg02,Hdg03
  real*8,dimension(ex(3)) :: Hdgx01,Hdgx02,Hdgx03
  real*8,dimension(ex(3)) :: Hdgy01,Hdgy02,Hdgy03
  real*8,dimension(ex(3)) :: Hdgx22,Hdgx23,Hdgx33
  real*8,dimension(ex(3)) :: Hdgy22,Hdgy23,Hdgy33
  real*8,dimension(ex(3)) :: Hddgxx01,Hddgxx33,Hddgyy01,Hddgyy22,Hddgxy23
  real*8,dimension(ex(3)) :: Hddgxy01,Hddgxr02,Hddgxr03,Hddgyr02,Hddgyr03
  real*8,dimension(ex(3)) :: Hddgxr22,Hddgxr23,Hddgxr33,Hddgyr22,Hddgyr23,Hddgyr33

  real*8,dimension(ex(3)) :: Hg00,Hg01,Hg02,Hg03,Hg22,Hg23,Hg33
  real*8,dimension(ex(3)) :: HTheta22,HTheta23,HTheta33

  real*8 :: Theta220,Theta221,Theta22h,Theta22_rhs
  real*8 :: Theta230,Theta231,Theta23h,Theta23_rhs
  real*8 :: Theta330,Theta331,Theta33h,Theta33_rhs
  integer :: i,j,k,RK4

!!! sanity check
  dR = sum(g22)+sum(g23)+sum(g33)+sum(g01) &
      +sum(g00)+sum(g02)+sum(g03) &
      +sum(Theta22)+sum(Theta23)+sum(Theta33)
  if(dR.ne.dR) then
     if(sum(g22).ne.sum(g22))write(*,*)"NullEvol_Theta: find NaN in g22"
     if(sum(g23).ne.sum(g23))write(*,*)"NullEvol_Theta: find NaN in g23"
     if(sum(g33).ne.sum(g33))write(*,*)"NullEvol_Theta: find NaN in g33"
     if(sum(g01).ne.sum(g01))write(*,*)"NullEvol_Theta: find NaN in g01"
     if(sum(g00).ne.sum(g00))write(*,*)"NullEvol_Theta: find NaN in g00"
     if(sum(g02).ne.sum(g02))write(*,*)"NullEvol_Theta: find NaN in g02"
     if(sum(g03).ne.sum(g03))write(*,*)"NullEvol_Theta: find NaN in g03"
     if(sum(Theta22).ne.sum(Theta22))write(*,*)"NullEvol_Theta: find NaN in Theta22"
     if(sum(Theta23).ne.sum(Theta23))write(*,*)"NullEvol_Theta: find NaN in Theta23"
     if(sum(Theta33).ne.sum(Theta33))write(*,*)"NullEvol_Theta: find NaN in Theta33"
     gont = 1
     return
  endif

  dR = R(2) - R(1)

  do j=1,ex(2)
  do i=1,ex(1)
     call rderivs_x(ex(3),R,g00(i,j,:),dg00)
     call rderivs_x(ex(3),R,g01(i,j,:),dg01)
     dg02 = p02(i,j,:)
     dg03 = p03(i,j,:)

     call rderivs_x(ex(3),R,g22(i,j,:),dg22)
     call rderivs_x(ex(3),R,g23(i,j,:),dg23)
     call rderivs_x(ex(3),R,g33(i,j,:),dg33)
     call rdderivs_x(ex(3),R,g22(i,j,:),ddg22)
     call rdderivs_x(ex(3),R,g23(i,j,:),ddg23)
     call rdderivs_x(ex(3),R,g33(i,j,:),ddg33)

   do k=1,ex(3)
     call rderivs_x_point(ex(1),crho,g01(:,j,k),dgx01(k),i)
     call rderivs_x_point(ex(1),crho,g02(:,j,k),dgx02(k),i)
     call rderivs_x_point(ex(1),crho,g03(:,j,k),dgx03(k),i)

     call rderivs_x_point(ex(2),sigma,g01(i,:,k),dgy01(k),j)
     call rderivs_x_point(ex(2),sigma,g02(i,:,k),dgy02(k),j)
     call rderivs_x_point(ex(2),sigma,g03(i,:,k),dgy03(k),j)

     call rderivs_x_point(ex(1),crho,g22(:,j,k),dgx22(k),i)
     call rderivs_x_point(ex(1),crho,g23(:,j,k),dgx23(k),i)
     call rderivs_x_point(ex(1),crho,g33(:,j,k),dgx33(k),i)

     call rderivs_x_point(ex(2),sigma,g22(i,:,k),dgy22(k),j)
     call rderivs_x_point(ex(2),sigma,g23(i,:,k),dgy23(k),j)
     call rderivs_x_point(ex(2),sigma,g33(i,:,k),dgy33(k),j)

     call rdderivs_x_point(ex(1),crho,g01(:,j,k),ddgxx01(k),i)
     call rdderivs_x_point(ex(1),crho,g33(:,j,k),ddgxx33(k),i)

     call rdderivs_x_point(ex(2),sigma,g01(i,:,k),ddgyy01(k),j)
     call rdderivs_x_point(ex(2),sigma,g22(i,:,k),ddgyy22(k),j)

     call rderivs_x_point(ex(1),crho,p02(:,j,k),ddgxr02(k),i)
     call rderivs_x_point(ex(1),crho,p03(:,j,k),ddgxr03(k),i)

     call rderivs_x_point(ex(2),sigma,p02(i,:,k),ddgyr02(k),j)
     call rderivs_x_point(ex(2),sigma,p03(i,:,k),ddgyr03(k),j)

     call rdderivs_xy_point(ex(1),ex(2),crho,sigma,g01(:,:,k),ddgxy01(k),i,j)
     call rdderivs_xy_point(ex(1),ex(2),crho,sigma,g23(:,:,k),ddgxy23(k),i,j)

     call rdderivs_xy_point(ex(1),ex(3),crho,R,g22(:,j,:),ddgxr22(k),i,k)
     call rdderivs_xy_point(ex(1),ex(3),crho,R,g23(:,j,:),ddgxr23(k),i,k)
     call rdderivs_xy_point(ex(1),ex(3),crho,R,g33(:,j,:),ddgxr33(k),i,k)

     call rdderivs_xy_point(ex(2),ex(3),sigma,R,g22(i,:,:),ddgyr22(k),j,k)
     call rdderivs_xy_point(ex(2),ex(3),sigma,R,g23(i,:,:),ddgyr23(k),j,k)
     call rdderivs_xy_point(ex(2),ex(3),sigma,R,g33(i,:,:),ddgyr33(k),j,k)
   enddo

     call rget_half_x(ex(3),g00(i,j,:),Hg00)
     call rget_half_x(ex(3),g01(i,j,:),Hg01)
     call rget_half_x(ex(3),g02(i,j,:),Hg02)
     call rget_half_x(ex(3),g03(i,j,:),Hg03)
     call rget_half_x(ex(3),g22(i,j,:),Hg22)
     call rget_half_x(ex(3),g23(i,j,:),Hg23)
     call rget_half_x(ex(3),g33(i,j,:),Hg33)
     call rget_half_x(ex(3),Theta22(i,j,:),HTheta22)
     call rget_half_x(ex(3),Theta23(i,j,:),HTheta23)
     call rget_half_x(ex(3),Theta33(i,j,:),HTheta33)

     call rget_half_x(ex(3),dg22,Hdg22)
     call rget_half_x(ex(3),dg23,Hdg23)
     call rget_half_x(ex(3),dg33,Hdg33)
     call rget_half_x(ex(3),ddg22,Hddg22)
     call rget_half_x(ex(3),ddg23,Hddg23)
     call rget_half_x(ex(3),ddg33,Hddg33)
     call rget_half_x(ex(3),dg00,Hdg00)
     call rget_half_x(ex(3),dg01,Hdg01)
     call rget_half_x(ex(3),dg02,Hdg02)
     call rget_half_x(ex(3),dg03,Hdg03)
     call rget_half_x(ex(3),dgx01,Hdgx01)
     call rget_half_x(ex(3),dgx02,Hdgx02)
     call rget_half_x(ex(3),dgx03,Hdgx03)
     call rget_half_x(ex(3),dgy01,Hdgy01)
     call rget_half_x(ex(3),dgy02,Hdgy02)
     call rget_half_x(ex(3),dgy03,Hdgy03)
     call rget_half_x(ex(3),dgx22,Hdgx22)
     call rget_half_x(ex(3),dgx23,Hdgx23)
     call rget_half_x(ex(3),dgx33,Hdgx33)
     call rget_half_x(ex(3),dgy22,Hdgy22)
     call rget_half_x(ex(3),dgy23,Hdgy23)
     call rget_half_x(ex(3),dgy33,Hdgy33)
     call rget_half_x(ex(3),ddgxx01,Hddgxx01)
     call rget_half_x(ex(3),ddgxx33,Hddgxx33)
     call rget_half_x(ex(3),ddgyy01,Hddgyy01)
     call rget_half_x(ex(3),ddgyy22,Hddgyy22)
     call rget_half_x(ex(3),ddgxy23,Hddgxy23)
     call rget_half_x(ex(3),ddgxy01,Hddgxy01)
     call rget_half_x(ex(3),ddgxr02,Hddgxr02)
     call rget_half_x(ex(3),ddgxr03,Hddgxr03)
     call rget_half_x(ex(3),ddgyr02,Hddgyr02)
     call rget_half_x(ex(3),ddgyr03,Hddgyr03)
     call rget_half_x(ex(3),ddgxr22,Hddgxr22)
     call rget_half_x(ex(3),ddgxr23,Hddgxr23)
     call rget_half_x(ex(3),ddgxr33,Hddgxr33)
     call rget_half_x(ex(3),ddgyr22,Hddgyr22)
     call rget_half_x(ex(3),ddgyr23,Hddgyr23)
     call rget_half_x(ex(3),ddgyr33,Hddgyr33)

#if 0
     Theta220 = Theta22(i,j,1)
     Theta230 = Theta23(i,j,1)
     Theta330 = Theta33(i,j,1)

  do k=1,ex(3)-2
        RK4 = 0
        call Theta_rhs2(Rmin,R(k),g00(i,j,k),g02(i,j,k),g03(i,j,k),g22(i,j,k),g23(i,j,k),g33(i,j,k), &
                     dg22(k),dg23(k),dg33(k),ddg22(k),ddg23(k),ddg33(k),g01(i,j,k), &
                     Theta220,Theta230,Theta330, &
                    dg01(k),dg02(k),dg03(k),          &
                    dgx01(k),dgx02(k),dgx03(k),       &
                    dgy01(k),dgy02(k),dgy03(k),       &
                    dgx22(k),dgx23(k),dgx33(k),       &
                    dgy22(k),dgy23(k),dgy33(k),       &
                    dg00(k), &
                    ddgxx01(k), &
                    ddgxx33(k), &
                    ddgyy01(k), &
                    ddgyy22(k), &
                    ddgxy23(k), &
                    ddgxy01(k), &
                    ddgxr02(k),ddgxr03(k), &
                    ddgyr02(k),ddgyr03(k), &
                    ddgxr22(k),ddgxr23(k),ddgxr33(k), &
                    ddgyr22(k),ddgyr23(k),ddgyr33(k), &
                    Theta22_rhs,Theta23_rhs,Theta33_rhs)
        call rungekutta4_scalar(dR,Theta220,Theta22h,Theta22_rhs,RK4)
        call rungekutta4_scalar(dR,Theta230,Theta23h,Theta23_rhs,RK4)
        call rungekutta4_scalar(dR,Theta330,Theta33h,Theta33_rhs,RK4)

        RK4 = 1

        call Theta_rhs2(Rmin,R(k)+dR/2,Hg00(k),Hg02(k),Hg03(k),Hg22(k),Hg23(k),Hg33(k), &
                     Hdg22(k),Hdg23(k),Hdg33(k),Hddg22(k),Hddg23(k),Hddg33(k),Hg01(k), &
                     Theta22h,Theta23h,Theta33h, &
                    Hdg01(k),Hdg02(k),Hdg03(k),          &
                    Hdgx01(k),Hdgx02(k),Hdgx03(k),       &
                    Hdgy01(k),Hdgy02(k),Hdgy03(k),       &
                    Hdgx22(k),Hdgx23(k),Hdgx33(k),       &
                    Hdgy22(k),Hdgy23(k),Hdgy33(k),       &
                    Hdg00(k), &
                    Hddgxx01(k), &
                    Hddgxx33(k), &
                    Hddgyy01(k), &
                    Hddgyy22(k), &
                    Hddgxy23(k), &
                    Hddgxy01(k), &
                    Hddgxr02(k),Hddgxr03(k), &
                    Hddgyr02(k),Hddgyr03(k), &
                    Hddgxr22(k),Hddgxr23(k),Hddgxr33(k), &
                    Hddgyr22(k),Hddgyr23(k),Hddgyr33(k), &
                    Theta221,Theta231,Theta331)

        call rungekutta4_scalar(dR,Theta220,Theta221,Theta22_rhs,RK4)
        call rungekutta4_scalar(dR,Theta230,Theta231,Theta23_rhs,RK4)
        call rungekutta4_scalar(dR,Theta330,Theta331,Theta33_rhs,RK4)
        call rswap(Theta22h,Theta221)
        call rswap(Theta23h,Theta231)
        call rswap(Theta33h,Theta331)

        RK4 = 2
        call Theta_rhs2(Rmin,R(k)+dR/2,Hg00(k),Hg02(k),Hg03(k),Hg22(k),Hg23(k),Hg33(k), &
                     Hdg22(k),Hdg23(k),Hdg33(k),Hddg22(k),Hddg23(k),Hddg33(k),Hg01(k), &
                     Theta22h,Theta23h,Theta33h, &
                    Hdg01(k),Hdg02(k),Hdg03(k),          &
                    Hdgx01(k),Hdgx02(k),Hdgx03(k),       &
                    Hdgy01(k),Hdgy02(k),Hdgy03(k),       &
                    Hdgx22(k),Hdgx23(k),Hdgx33(k),       &
                    Hdgy22(k),Hdgy23(k),Hdgy33(k),       &
                    Hdg00(k), &
                    Hddgxx01(k), &
                    Hddgxx33(k), &
                    Hddgyy01(k), &
                    Hddgyy22(k), &
                    Hddgxy23(k), &
                    Hddgxy01(k), &
                    Hddgxr02(k),Hddgxr03(k), &
                    Hddgyr02(k),Hddgyr03(k), &
                    Hddgxr22(k),Hddgxr23(k),Hddgxr33(k), &
                    Hddgyr22(k),Hddgyr23(k),Hddgyr33(k), &
                    Theta221,Theta231,Theta331)

        call rungekutta4_scalar(dR,Theta220,Theta221,Theta22_rhs,RK4)
        call rungekutta4_scalar(dR,Theta230,Theta231,Theta23_rhs,RK4)
        call rungekutta4_scalar(dR,Theta330,Theta331,Theta33_rhs,RK4)
        call rswap(Theta22h,Theta221)
        call rswap(Theta23h,Theta231)
        call rswap(Theta33h,Theta331)

        RK4 = 3
        call Theta_rhs2(Rmin,R(k+1),g00(i,j,k+1),g02(i,j,k+1),g03(i,j,k+1),g22(i,j,k+1),g23(i,j,k+1),g33(i,j,k+1), &
                     dg22(k+1),dg23(k+1),dg33(k+1),ddg22(k+1),ddg23(k+1),ddg33(k+1),g01(i,j,k+1), &
                     Theta22h,Theta23h,Theta33h, &
                    dg01(k+1),dg02(k+1),dg03(k+1),          &
                    dgx01(k+1),dgx02(k+1),dgx03(k+1),       &
                    dgy01(k+1),dgy02(k+1),dgy03(k+1),       &
                    dgx22(k+1),dgx23(k+1),dgx33(k+1),       &
                    dgy22(k+1),dgy23(k+1),dgy33(k+1),       &
                    dg00(k+1), &
                    ddgxx01(k+1), &
                    ddgxx33(k+1), &
                    ddgyy01(k+1), &
                    ddgyy22(k+1), &
                    ddgxy23(k+1), &
                    ddgxy01(k+1), &
                    ddgxr02(k+1),ddgxr03(k+1), &
                    ddgyr02(k+1),ddgyr03(k+1), &
                    ddgxr22(k+1),ddgxr23(k+1),ddgxr33(k+1), &
                    ddgyr22(k+1),ddgyr23(k+1),ddgyr33(k+1), &
                    Theta221,Theta231,Theta331)

        call rungekutta4_scalar(dR,Theta220,Theta221,Theta22_rhs,RK4)
        call rungekutta4_scalar(dR,Theta230,Theta231,Theta23_rhs,RK4)
        call rungekutta4_scalar(dR,Theta330,Theta331,Theta33_rhs,RK4)
        call rswap(Theta220,Theta221)
        call rswap(Theta230,Theta231)
        call rswap(Theta330,Theta331)

        Theta22(i,j,k+1) = Theta220
        Theta23(i,j,k+1) = Theta230
        Theta33(i,j,k+1) = Theta330
      enddo

        k=ex(3)-1
!    closing step

        call Theta_rhs2(Rmin,R(k),g00(i,j,k),g02(i,j,k),g03(i,j,k),g22(i,j,k),g23(i,j,k),g33(i,j,k), &
                     dg22(k),dg23(k),dg33(k),ddg22(k),ddg23(k),ddg33(k),g01(i,j,k), &
                     Theta22(i,j,k),Theta23(i,j,k),Theta33(i,j,k), &
                    dg01(k),dg02(k),dg03(k),          &
                    dgx01(k),dgx02(k),dgx03(k),       &
                    dgy01(k),dgy02(k),dgy03(k),       &
                    dgx22(k),dgx23(k),dgx33(k),       &
                    dgy22(k),dgy23(k),dgy33(k),       &
                    dg00(k), &
                    ddgxx01(k), &
                    ddgxx33(k), &
                    ddgyy01(k), &
                    ddgyy22(k), &
                    ddgxy23(k), &
                    ddgxy01(k), &
                    ddgxr02(k),ddgxr03(k), &
                    ddgyr02(k),ddgyr03(k), &
                    ddgxr22(k),ddgxr23(k),ddgxr33(k), &
                    ddgyr22(k),ddgyr23(k),ddgyr33(k), &
                    Theta22_rhs,Theta23_rhs,Theta33_rhs)

     Theta22(i,j,k+1) = Theta22(i,j,k) + Theta22_rhs*dR
     Theta23(i,j,k+1) = Theta23(i,j,k) + Theta23_rhs*dR
     Theta33(i,j,k+1) = Theta33(i,j,k) + Theta33_rhs*dR

#endif
  enddo
  enddo

  gont = 0

  return

end function NullEvol_Theta2
!---------------------------------------------------------------------------------
subroutine Theta_rhs2(Rmin,r,g00,g02,g03,g22,g23,g33,dg22,dg23,dg33,ddg22,ddg23,ddg33,g01, &
                     Theta22,Theta23,Theta33, &
                    dg01,dg02,dg03,          &
                    dgx01,dgx02,dgx03,       &
                    dgy01,dgy02,dgy03,       &
                    dgx22,dgx23,dgx33,       &
                    dgy22,dgy23,dgy33,       &
                    dg00, &
                    ddgxx01, &
                    ddgxx33, &
                    ddgyy01, &
                    ddgyy22, &
                    ddgxy23, &
                    ddgxy01, &
                    ddgxr02,ddgxr03, &
                    ddgyr02,ddgyr03, &
                    ddgxr22,ddgxr23,ddgxr33, &
                    ddgyr22,ddgyr23,ddgyr33, &
                    Theta22_rhs,Theta23_rhs,Theta33_rhs)

  implicit none

!~~~~~~% Input parameters:
  real*8,intent(in) :: Rmin,r,g00,g02,g03,g22,g23,g33,dg22,dg23,dg33,ddg22,ddg23,ddg33,g01
  real*8,intent(in) :: Theta22,Theta23,Theta33,dg01,dg02,dg03
  real*8,intent(in) :: dgx01,dgx02,dgx03,dgx22,dgx23,dgx33
  real*8,intent(in) :: dgy01,dgy02,dgy03,dgy22,dgy23,dgy33
  real*8,intent(in) :: dg00
  real*8,intent(out) :: Theta22_rhs,Theta23_rhs,Theta33_rhs
  real*8,intent(in) ::  ddgxx01
  real*8,intent(in) ::  ddgxx33
  real*8,intent(in) ::  ddgyy01
  real*8,intent(in) ::  ddgyy22
  real*8,intent(in) ::  ddgxy23
  real*8,intent(in) ::  ddgxy01
  real*8,intent(in) ::  ddgxr02,ddgxr03
  real*8,intent(in) ::  ddgyr02,ddgyr03
  real*8,intent(in) ::  ddgxr22,ddgxr23,ddgxr33,ddgyr22,ddgyr23,ddgyr33

  real*8 ::  t1;
  real*8 ::  t100;
  real*8 ::  t1001;
  real*8 ::  t1009;
  real*8 ::  t1010;
  real*8 ::  t1011;
  real*8 ::  t1015;
  real*8 ::  t1019;
  real*8 ::  t1023;
  real*8 ::  t1037;
  real*8 ::  t104;
  real*8 ::  t1041;
  real*8 ::  t1042;
  real*8 ::  t1049;
  real*8 ::  t1065;
  real*8 ::  t1070;
  real*8 ::  t1090;
  real*8 ::  t1094;
  real*8 ::  t1099;
  real*8 ::  t11;
  real*8 ::  t111;
  real*8 ::  t1113;
  real*8 ::  t112;
  real*8 ::  t1123;
  real*8 ::  t1126;
  real*8 ::  t1130;
  real*8 ::  t1134;
  real*8 ::  t1160;
  real*8 ::  t1173;
  real*8 ::  t1174;
  real*8 ::  t1180;
  real*8 ::  t12;
  real*8 ::  t1207;
  real*8 ::  t1211;
  real*8 ::  t1218;
  real*8 ::  t1222;
  real*8 ::  t1223;
  real*8 ::  t1226;
  real*8 ::  t1227;
  real*8 ::  t1230;
  real*8 ::  t1231;
  real*8 ::  t1234;
  real*8 ::  t1240;
  real*8 ::  t1242;
  real*8 ::  t1245;
  real*8 ::  t1248;
  real*8 ::  t125;
  real*8 ::  t1250;
  real*8 ::  t1254;
  real*8 ::  t1265;
  real*8 ::  t1272;
  real*8 ::  t1277;
  real*8 ::  t1281;
  real*8 ::  t1282;
  real*8 ::  t1287;
  real*8 ::  t1296;
  real*8 ::  t13;
  real*8 ::  t1301;
  real*8 ::  t1308;
  real*8 ::  t1311;
  real*8 ::  t1325;
  real*8 ::  t1326;
  real*8 ::  t1330;
  real*8 ::  t1334;
  real*8 ::  t1335;
  real*8 ::  t1338;
  real*8 ::  t1348;
  real*8 ::  t1351;
  real*8 ::  t1354;
  real*8 ::  t1386;
  real*8 ::  t1398;
  real*8 ::  t1411;
  real*8 ::  t142;
  real*8 ::  t1426;
  real*8 ::  t143;
  real*8 ::  t1432;
  real*8 ::  t1437;
  real*8 ::  t144;
  real*8 ::  t1441;
  real*8 ::  t1449;
  real*8 ::  t1475;
  real*8 ::  t148;
  real*8 ::  t1483;
  real*8 ::  t1496;
  real*8 ::  t1506;
  real*8 ::  t152;
  real*8 ::  t1522;
  real*8 ::  t1523;
  real*8 ::  t1526;
  real*8 ::  t1529;
  real*8 ::  t1532;
  real*8 ::  t1535;
  real*8 ::  t1536;
  real*8 ::  t1539;
  real*8 ::  t1540;
  real*8 ::  t1543;
  real*8 ::  t1547;
  real*8 ::  t1556;
  real*8 ::  t1592;
  real*8 ::  t1598;
  real*8 ::  t1601;
  real*8 ::  t1604;
  real*8 ::  t162;
  real*8 ::  t1629;
  real*8 ::  t1636;
  real*8 ::  t1641;
  real*8 ::  t1646;
  real*8 ::  t1647;
  real*8 ::  t1652;
  real*8 ::  t1653;
  real*8 ::  t1654;
  real*8 ::  t1668;
  real*8 ::  t1673;
  real*8 ::  t1674;
  real*8 ::  t1678;
  real*8 ::  t1682;
  real*8 ::  t1686;
  real*8 ::  t1691;
  real*8 ::  t1694;
  real*8 ::  t1695;
  real*8 ::  t1697;
  real*8 ::  t17;
  real*8 ::  t170;
  real*8 ::  t1700;
  real*8 ::  t1701;
  real*8 ::  t1703;
  real*8 ::  t1706;
  real*8 ::  t1707;
  real*8 ::  t1710;
  real*8 ::  t1711;
  real*8 ::  t1712;
  real*8 ::  t1716;
  real*8 ::  t1717;
  real*8 ::  t1718;
  real*8 ::  t1720;
  real*8 ::  t1727;
  real*8 ::  t1728;
  real*8 ::  t1731;
  real*8 ::  t1733;
  real*8 ::  t1737;
  real*8 ::  t1740;
  real*8 ::  t1744;
  real*8 ::  t1747;
  real*8 ::  t1760;
  real*8 ::  t1764;
  real*8 ::  t1768;
  real*8 ::  t177;
  real*8 ::  t1787;
  real*8 ::  t18;
  real*8 ::  t1813;
  real*8 ::  t1817;
  real*8 ::  t1820;
  real*8 ::  t1822;
  real*8 ::  t1825;
  real*8 ::  t1828;
  real*8 ::  t1833;
  real*8 ::  t1847;
  real*8 ::  t185;
  real*8 ::  t1873;
  real*8 ::  t1876;
  real*8 ::  t1882;
  real*8 ::  t1884;
  real*8 ::  t1887;
  real*8 ::  t1891;
  real*8 ::  t1896;
  real*8 ::  t1897;
  real*8 ::  t19;
  real*8 ::  t1901;
  real*8 ::  t1904;
  real*8 ::  t1906;
  real*8 ::  t1909;
  real*8 ::  t1910;
  real*8 ::  t1914;
  real*8 ::  t192;
  real*8 ::  t1932;
  real*8 ::  t1934;
  real*8 ::  t1935;
  real*8 ::  t1936;
  real*8 ::  t1939;
  real*8 ::  t1942;
  real*8 ::  t1943;
  real*8 ::  t1946;
  real*8 ::  t1949;
  real*8 ::  t197;
  real*8 ::  t1973;
  real*8 ::  t198;
  real*8 ::  t1982;
  real*8 ::  t199;
  real*8 ::  t1995;
  real*8 ::  t1998;
  real*8 ::  t20;
  real*8 ::  t201;
  real*8 ::  t202;
  real*8 ::  t2035;
  real*8 ::  t205;
  real*8 ::  t207;
  real*8 ::  t211;
  real*8 ::  t22;
  real*8 ::  t234;
  real*8 ::  t249;
  real*8 ::  t25;
  real*8 ::  t265;
  real*8 ::  t266;
  real*8 ::  t267;
  real*8 ::  t27;
  real*8 ::  t270;
  real*8 ::  t273;
  real*8 ::  t274;
  real*8 ::  t277;
  real*8 ::  t278;
  real*8 ::  t279;
  real*8 ::  t285;
  real*8 ::  t3;
  real*8 ::  t301;
  real*8 ::  t304;
  real*8 ::  t305;
  real*8 ::  t306;
  real*8 ::  t31;
  real*8 ::  t315;
  real*8 ::  t320;
  real*8 ::  t321;
  real*8 ::  t325;
  real*8 ::  t326;
  real*8 ::  t327;
  real*8 ::  t329;
  real*8 ::  t333;
  real*8 ::  t336;
  real*8 ::  t337;
  real*8 ::  t338;
  real*8 ::  t339;
  real*8 ::  t341;
  real*8 ::  t348;
  real*8 ::  t35;
  real*8 ::  t355;
  real*8 ::  t364;
  real*8 ::  t365;
  real*8 ::  t366;
  real*8 ::  t367;
  real*8 ::  t368;
  real*8 ::  t371;
  real*8 ::  t372;
  real*8 ::  t373;
  real*8 ::  t377;
  real*8 ::  t378;
  real*8 ::  t382;
  real*8 ::  t385;
  real*8 ::  t386;
  real*8 ::  t387;
  real*8 ::  t388;
  real*8 ::  t39;
  real*8 ::  t392;
  real*8 ::  t393;
  real*8 ::  t397;
  real*8 ::  t4;
  real*8 ::  t401;
  real*8 ::  t402;
  real*8 ::  t406;
  real*8 ::  t407;
  real*8 ::  t408;
  real*8 ::  t411;
  real*8 ::  t412;
  real*8 ::  t415;
  real*8 ::  t416;
  real*8 ::  t417;
  real*8 ::  t42;
  real*8 ::  t420;
  real*8 ::  t421;
  real*8 ::  t422;
  real*8 ::  t426;
  real*8 ::  t427;
  real*8 ::  t43;
  real*8 ::  t430;
  real*8 ::  t431;
  real*8 ::  t432;
  real*8 ::  t435;
  real*8 ::  t436;
  real*8 ::  t437;
  real*8 ::  t440;
  real*8 ::  t441;
  real*8 ::  t444;
  real*8 ::  t448;
  real*8 ::  t449;
  real*8 ::  t453;
  real*8 ::  t454;
  real*8 ::  t455;
  real*8 ::  t458;
  real*8 ::  t461;
  real*8 ::  t462;
  real*8 ::  t465;
  real*8 ::  t466;
  real*8 ::  t469;
  real*8 ::  t470;
  real*8 ::  t473;
  real*8 ::  t474;
  real*8 ::  t477;
  real*8 ::  t479;
  real*8 ::  t48;
  real*8 ::  t480;
  real*8 ::  t483;
  real*8 ::  t484;
  real*8 ::  t487;
  real*8 ::  t488;
  real*8 ::  t491;
  real*8 ::  t495;
  real*8 ::  t496;
  real*8 ::  t5;
  real*8 ::  t500;
  real*8 ::  t501;
  real*8 ::  t504;
  real*8 ::  t505;
  real*8 ::  t508;
  real*8 ::  t509;
  real*8 ::  t510;
  real*8 ::  t516;
  real*8 ::  t519;
  real*8 ::  t52;
  real*8 ::  t522;
  real*8 ::  t523;
  real*8 ::  t524;
  real*8 ::  t525;
  real*8 ::  t528;
  real*8 ::  t529;
  real*8 ::  t532;
  real*8 ::  t535;
  real*8 ::  t541;
  real*8 ::  t549;
  real*8 ::  t55;
  real*8 ::  t552;
  real*8 ::  t553;
  real*8 ::  t56;
  real*8 ::  t561;
  real*8 ::  t564;
  real*8 ::  t569;
  real*8 ::  t57;
  real*8 ::  t572;
  real*8 ::  t575;
  real*8 ::  t576;
  real*8 ::  t577;
  real*8 ::  t579;
  real*8 ::  t582;
  real*8 ::  t586;
  real*8 ::  t589;
  real*8 ::  t590;
  real*8 ::  t591;
  real*8 ::  t594;
  real*8 ::  t595;
  real*8 ::  t6;
  real*8 ::  t605;
  real*8 ::  t61;
  real*8 ::  t610;
  real*8 ::  t611;
  real*8 ::  t618;
  real*8 ::  t622;
  real*8 ::  t623;
  real*8 ::  t624;
  real*8 ::  t627;
  real*8 ::  t631;
  real*8 ::  t634;
  real*8 ::  t638;
  real*8 ::  t639;
  real*8 ::  t640;
  real*8 ::  t643;
  real*8 ::  t644;
  real*8 ::  t645;
  real*8 ::  t648;
  real*8 ::  t649;
  real*8 ::  t658;
  real*8 ::  t659;
  real*8 ::  t660;
  real*8 ::  t663;
  real*8 ::  t664;
  real*8 ::  t668;
  real*8 ::  t671;
  real*8 ::  t686;
  real*8 ::  t7;
  real*8 ::  t70;
  real*8 ::  t706;
  real*8 ::  t710;
  real*8 ::  t713;
  real*8 ::  t717;
  real*8 ::  t723;
  real*8 ::  t725;
  real*8 ::  t728;
  real*8 ::  t731;
  real*8 ::  t733;
  real*8 ::  t738;
  real*8 ::  t741;
  real*8 ::  t742;
  real*8 ::  t746;
  real*8 ::  t749;
  real*8 ::  t750;
  real*8 ::  t751;
  real*8 ::  t754;
  real*8 ::  t755;
  real*8 ::  t758;
  real*8 ::  t77;
  real*8 ::  t775;
  real*8 ::  t780;
  real*8 ::  t782;
  real*8 ::  t783;
  real*8 ::  t786;
  real*8 ::  t787;
  real*8 ::  t788;
  real*8 ::  t792;
  real*8 ::  t796;
  real*8 ::  t799;
  real*8 ::  t800;
  real*8 ::  t804;
  real*8 ::  t811;
  real*8 ::  t812;
  real*8 ::  t822;
  real*8 ::  t831;
  real*8 ::  t832;
  real*8 ::  t835;
  real*8 ::  t836;
  real*8 ::  t837;
  real*8 ::  t84;
  real*8 ::  t850;
  real*8 ::  t855;
  real*8 ::  t856;
  real*8 ::  t857;
  real*8 ::  t860;
  real*8 ::  t862;
  real*8 ::  t865;
  real*8 ::  t871;
  real*8 ::  t876;
  real*8 ::  t88;
  real*8 ::  t880;
  real*8 ::  t884;
  real*8 ::  t888;
  real*8 ::  t889;
  real*8 ::  t892;
  real*8 ::  t895;
  real*8 ::  t898;
  real*8 ::  t901;
  real*8 ::  t904;
  real*8 ::  t92;
  real*8 ::  t922;
  real*8 ::  t925;
  real*8 ::  t928;
  real*8 ::  t929;
  real*8 ::  t93;
  real*8 ::  t932;
  real*8 ::  t935;
  real*8 ::  t938;
  real*8 ::  t956;
  real*8 ::  t959;
  real*8 ::  t960;
  real*8 ::  t963;
  real*8 ::  t970;
  real*8 ::  t975;
  real*8 ::  t979;
  real*8 ::  t980;
  real*8 ::  t983;
  real*8 ::  t985;
  real*8 ::  t991;
  real*8 ::  t996;

  real*8 ::  t10;
  real*8 ::  t1006;
  real*8 ::  t1007;
  real*8 ::  t1012;
  real*8 ::  t1030;
  real*8 ::  t1039;
  real*8 ::  t1044;
  real*8 ::  t1067;
  real*8 ::  t1084;
  real*8 ::  t1092;
  real*8 ::  t1100;
  real*8 ::  t1112;
  real*8 ::  t1117;
  real*8 ::  t1121;
  real*8 ::  t1122;
  real*8 ::  t1127;
  real*8 ::  t1131;
  real*8 ::  t1133;
  real*8 ::  t1138;
  real*8 ::  t1141;
  real*8 ::  t1142;
  real*8 ::  t1143;
  real*8 ::  t1144;
  real*8 ::  t1148;
  real*8 ::  t1166;
  real*8 ::  t1177;
  real*8 ::  t1181;
  real*8 ::  t1191;
  real*8 ::  t120;
  real*8 ::  t1203;
  real*8 ::  t1204;
  real*8 ::  t121;
  real*8 ::  t1212;
  real*8 ::  t1235;
  real*8 ::  t1239;
  real*8 ::  t124;
  real*8 ::  t1249;
  real*8 ::  t1252;
  real*8 ::  t1253;
  real*8 ::  t1256;
  real*8 ::  t128;
  real*8 ::  t1289;
  real*8 ::  t129;
  real*8 ::  t1291;
  real*8 ::  t1293;
  real*8 ::  t130;
  real*8 ::  t131;
  real*8 ::  t1313;
  real*8 ::  t1314;
  real*8 ::  t1317;
  real*8 ::  t132;
  real*8 ::  t1320;
  real*8 ::  t1322;
  real*8 ::  t133;
  real*8 ::  t1331;
  real*8 ::  t1332;
  real*8 ::  t1342;
  real*8 ::  t1355;
  real*8 ::  t1357;
  real*8 ::  t1359;
  real*8 ::  t136;
  real*8 ::  t1362;
  real*8 ::  t1366;
  real*8 ::  t137;
  real*8 ::  t1374;
  real*8 ::  t1379;
  real*8 ::  t138;
  real*8 ::  t1380;
  real*8 ::  t1381;
  real*8 ::  t1384;
  real*8 ::  t1385;
  real*8 ::  t1388;
  real*8 ::  t139;
  real*8 ::  t1391;
  real*8 ::  t1392;
  real*8 ::  t140;
  real*8 ::  t1405;
  real*8 ::  t1406;
  real*8 ::  t1409;
  real*8 ::  t1410;
  real*8 ::  t1413;
  real*8 ::  t1416;
  real*8 ::  t1417;
  real*8 ::  t1419;
  real*8 ::  t1421;
  real*8 ::  t1428;
  real*8 ::  t1434;
  real*8 ::  t1440;
  real*8 ::  t1444;
  real*8 ::  t145;
  real*8 ::  t1450;
  real*8 ::  t1454;
  real*8 ::  t1457;
  real*8 ::  t1473;
  real*8 ::  t1476;
  real*8 ::  t1488;
  real*8 ::  t1490;
  real*8 ::  t1501;
  real*8 ::  t1505;
  real*8 ::  t1510;
  real*8 ::  t1516;
  real*8 ::  t1577;
  real*8 ::  t1612;
  real*8 ::  t1615;
  real*8 ::  t1619;
  real*8 ::  t1624;
  real*8 ::  t1625;
  real*8 ::  t163;
  real*8 ::  t1634;
  real*8 ::  t1640;
  real*8 ::  t1644;
  real*8 ::  t1648;
  real*8 ::  t1651;
  real*8 ::  t1655;
  real*8 ::  t1660;
  real*8 ::  t1663;
  real*8 ::  t1664;
  real*8 ::  t168;
  real*8 ::  t1689;
  real*8 ::  t169;
  real*8 ::  t1690;
  real*8 ::  t1693;
  real*8 ::  t1696;
  real*8 ::  t1708;
  real*8 ::  t171;
  real*8 ::  t1724;
  real*8 ::  t174;
  real*8 ::  t1741;
  real*8 ::  t175;
  real*8 ::  t1752;
  real*8 ::  t176;
  real*8 ::  t1775;
  real*8 ::  t1783;
  real*8 ::  t1788;
  real*8 ::  t1791;
  real*8 ::  t1795;
  real*8 ::  t181;
  real*8 ::  t1823;
  real*8 ::  t1824;
  real*8 ::  t183;
  real*8 ::  t1836;
  real*8 ::  t1842;
  real*8 ::  t1852;
  real*8 ::  t1856;
  real*8 ::  t1859;
  real*8 ::  t186;
  real*8 ::  t1863;
  real*8 ::  t187;
  real*8 ::  t1875;
  real*8 ::  t1878;
  real*8 ::  t1883;
  real*8 ::  t189;
  real*8 ::  t1890;
  real*8 ::  t191;
  real*8 ::  t1918;
  real*8 ::  t1921;
  real*8 ::  t1927;
  real*8 ::  t1931;
  real*8 ::  t194;
  real*8 ::  t1952;
  real*8 ::  t196;
  real*8 ::  t1970;
  real*8 ::  t200;
  real*8 ::  t2003;
  real*8 ::  t2004;
  real*8 ::  t2008;
  real*8 ::  t2017;
  real*8 ::  t2024;
  real*8 ::  t2032;
  real*8 ::  t204;
  real*8 ::  t206;
  real*8 ::  t2065;
  real*8 ::  t208;
  real*8 ::  t2085;
  real*8 ::  t209;
  real*8 ::  t2091;
  real*8 ::  t2093;
  real*8 ::  t21;
  real*8 ::  t212;
  real*8 ::  t2122;
  real*8 ::  t213;
  real*8 ::  t2133;
  real*8 ::  t2138;
  real*8 ::  t214;
  real*8 ::  t215;
  real*8 ::  t2166;
  real*8 ::  t219;
  real*8 ::  t2192;
  real*8 ::  t2201;
  real*8 ::  t222;
  real*8 ::  t226;
  real*8 ::  t23;
  real*8 ::  t233;
  real*8 ::  t236;
  real*8 ::  t237;
  real*8 ::  t238;
  real*8 ::  t239;
  real*8 ::  t24;
  real*8 ::  t247;
  real*8 ::  t248;
  real*8 ::  t251;
  real*8 ::  t252;
  real*8 ::  t255;
  real*8 ::  t258;
  real*8 ::  t259;
  real*8 ::  t268;
  real*8 ::  t28;
  real*8 ::  t282;
  real*8 ::  t283;
  real*8 ::  t287;
  real*8 ::  t29;
  real*8 ::  t290;
  real*8 ::  t293;
  real*8 ::  t296;
  real*8 ::  t297;
  real*8 ::  t298;
  real*8 ::  t302;
  real*8 ::  t310;
  real*8 ::  t311;
  real*8 ::  t316;
  real*8 ::  t317;
  real*8 ::  t32;
  real*8 ::  t322;
  real*8 ::  t323;
  real*8 ::  t324;
  real*8 ::  t33;
  real*8 ::  t330;
  real*8 ::  t34;
  real*8 ::  t344;
  real*8 ::  t376;
  real*8 ::  t384;
  real*8 ::  t389;
  real*8 ::  t394;
  real*8 ::  t399;
  real*8 ::  t40;
  real*8 ::  t404;
  real*8 ::  t41;
  real*8 ::  t419;
  real*8 ::  t438;
  real*8 ::  t439;
  real*8 ::  t443;
  real*8 ::  t445;
  real*8 ::  t450;
  real*8 ::  t451;
  real*8 ::  t459;
  real*8 ::  t46;
  real*8 ::  t460;
  real*8 ::  t463;
  real*8 ::  t464;
  real*8 ::  t47;
  real*8 ::  t478;
  real*8 ::  t482;
  real*8 ::  t49;
  real*8 ::  t492;
  real*8 ::  t50;
  real*8 ::  t503;
  real*8 ::  t507;
  real*8 ::  t511;
  real*8 ::  t514;
  real*8 ::  t515;
  real*8 ::  t517;
  real*8 ::  t518;
  real*8 ::  t520;
  real*8 ::  t533;
  real*8 ::  t537;
  real*8 ::  t538;
  real*8 ::  t54;
  real*8 ::  t544;
  real*8 ::  t545;
  real*8 ::  t548;
  real*8 ::  t551;
  real*8 ::  t555;
  real*8 ::  t556;
  real*8 ::  t560;
  real*8 ::  t563;
  real*8 ::  t566;
  real*8 ::  t574;
  real*8 ::  t59;
  real*8 ::  t597;
  real*8 ::  t599;
  real*8 ::  t600;
  real*8 ::  t603;
  real*8 ::  t606;
  real*8 ::  t609;
  real*8 ::  t614;
  real*8 ::  t616;
  real*8 ::  t617;
  real*8 ::  t620;
  real*8 ::  t621;
  real*8 ::  t625;
  real*8 ::  t63;
  real*8 ::  t64;
  real*8 ::  t646;
  real*8 ::  t655;
  real*8 ::  t662;
  real*8 ::  t667;
  real*8 ::  t672;
  real*8 ::  t676;
  real*8 ::  t677;
  real*8 ::  t68;
  real*8 ::  t680;
  real*8 ::  t683;
  real*8 ::  t696;
  real*8 ::  t72;
  real*8 ::  t739;
  real*8 ::  t740;
  real*8 ::  t743;
  real*8 ::  t744;
  real*8 ::  t745;
  real*8 ::  t748;
  real*8 ::  t752;
  real*8 ::  t756;
  real*8 ::  t76;
  real*8 ::  t766;
  real*8 ::  t769;
  real*8 ::  t770;
  real*8 ::  t771;
  real*8 ::  t774;
  real*8 ::  t778;
  real*8 ::  t785;
  real*8 ::  t789;
  real*8 ::  t793;
  real*8 ::  t798;
  real*8 ::  t8;
  real*8 ::  t801;
  real*8 ::  t803;
  real*8 ::  t806;
  real*8 ::  t808;
  real*8 ::  t813;
  real*8 ::  t816;
  real*8 ::  t817;
  real*8 ::  t818;
  real*8 ::  t823;
  real*8 ::  t824;
  real*8 ::  t838;
  real*8 ::  t842;
  real*8 ::  t843;
  real*8 ::  t844;
  real*8 ::  t849;
  real*8 ::  t868;
  real*8 ::  t873;
  real*8 ::  t874;
  real*8 ::  t877;
  real*8 ::  t878;
  real*8 ::  t881;
  real*8 ::  t882;
  real*8 ::  t885;
  real*8 ::  t89;
  real*8 ::  t9;
  real*8 ::  t900;
  real*8 ::  t902;
  real*8 ::  t913;
  real*8 ::  t915;
  real*8 ::  t919;
  real*8 ::  t921;
  real*8 ::  t923;
  real*8 ::  t941;
  real*8 ::  t944;
  real*8 ::  t946;
  real*8 ::  t949;
  real*8 ::  t953;
  real*8 ::  t958;
  real*8 ::  t96;
  real*8 ::  t968;
  real*8 ::  t969;
  real*8 ::  t973;
  real*8 ::  t978;
  real*8 ::  t990;
  real*8 ::  t995;
  real*8 ::  t998;

  real*8 ::  t1004;
  real*8 ::  t1005;
  real*8 ::  t1008;
  real*8 ::  t1020;
  real*8 ::  t103;
  real*8 ::  t105;
  real*8 ::  t1051;
  real*8 ::  t1057;
  real*8 ::  t1064;
  real*8 ::  t107;
  real*8 ::  t1072;
  real*8 ::  t1077;
  real*8 ::  t108;
  real*8 ::  t1086;
  real*8 ::  t1111;
  real*8 ::  t1115;
  real*8 ::  t1119;
  real*8 ::  t1125;
  real*8 ::  t1128;
  real*8 ::  t113;
  real*8 ::  t1136;
  real*8 ::  t1139;
  real*8 ::  t1146;
  real*8 ::  t1149;
  real*8 ::  t1151;
  real*8 ::  t1154;
  real*8 ::  t1158;
  real*8 ::  t1164;
  real*8 ::  t1168;
  real*8 ::  t117;
  real*8 ::  t1175;
  real*8 ::  t1185;
  real*8 ::  t1186;
  real*8 ::  t1208;
  real*8 ::  t1216;
  real*8 ::  t1232;
  real*8 ::  t1276;
  real*8 ::  t1284;
  real*8 ::  t1286;
  real*8 ::  t1290;
  real*8 ::  t1295;
  real*8 ::  t1298;
  real*8 ::  t1299;
  real*8 ::  t1302;
  real*8 ::  t1305;
  real*8 ::  t1309;
  real*8 ::  t1323;
  real*8 ::  t1324;
  real*8 ::  t1328;
  real*8 ::  t1333;
  real*8 ::  t1336;
  real*8 ::  t1337;
  real*8 ::  t1340;
  real*8 ::  t1341;
  real*8 ::  t1344;
  real*8 ::  t1365;
  real*8 ::  t1382;
  real*8 ::  t1389;
  real*8 ::  t1390;
  real*8 ::  t1393;
  real*8 ::  t1396;
  real*8 ::  t1397;
  real*8 ::  t14;
  real*8 ::  t1400;
  real*8 ::  t1401;
  real*8 ::  t1404;
  real*8 ::  t1408;
  real*8 ::  t1412;
  real*8 ::  t1415;
  real*8 ::  t1418;
  real*8 ::  t1424;
  real*8 ::  t1427;
  real*8 ::  t1430;
  real*8 ::  t1436;
  real*8 ::  t1439;
  real*8 ::  t1442;
  real*8 ::  t1445;
  real*8 ::  t1448;
  real*8 ::  t1453;
  real*8 ::  t1456;
  real*8 ::  t1471;
  real*8 ::  t1477;
  real*8 ::  t1481;
  real*8 ::  t1485;
  real*8 ::  t1487;
  real*8 ::  t1493;
  real*8 ::  t15;
  real*8 ::  t1500;
  real*8 ::  t1504;
  real*8 ::  t1507;
  real*8 ::  t1509;
  real*8 ::  t151;
  real*8 ::  t1513;
  real*8 ::  t1517;
  real*8 ::  t1521;
  real*8 ::  t1527;
  real*8 ::  t1541;
  real*8 ::  t1550;
  real*8 ::  t1553;
  real*8 ::  t1557;
  real*8 ::  t1558;
  real*8 ::  t1578;
  real*8 ::  t1579;
  real*8 ::  t158;
  real*8 ::  t1582;
  real*8 ::  t1583;
  real*8 ::  t1587;
  real*8 ::  t1595;
  real*8 ::  t1600;
  real*8 ::  t1603;
  real*8 ::  t1608;
  real*8 ::  t161;
  real*8 ::  t1620;
  real*8 ::  t1626;
  real*8 ::  t165;
  real*8 ::  t1656;
  real*8 ::  t1661;
  real*8 ::  t1662;
  real*8 ::  t1665;
  real*8 ::  t1666;
  real*8 ::  t1677;
  real*8 ::  t1681;
  real*8 ::  t1685;
  real*8 ::  t1692;
  real*8 ::  t1721;
  real*8 ::  t1722;
  real*8 ::  t1726;
  real*8 ::  t173;
  real*8 ::  t1730;
  real*8 ::  t1743;
  real*8 ::  t1745;
  real*8 ::  t1756;
  real*8 ::  t1761;
  real*8 ::  t1780;
  real*8 ::  t1807;
  real*8 ::  t1812;
  real*8 ::  t1845;
  real*8 ::  t1846;
  real*8 ::  t1854;
  real*8 ::  t1855;
  real*8 ::  t1860;
  real*8 ::  t1864;
  real*8 ::  t1869;
  real*8 ::  t188;
  real*8 ::  t1888;
  real*8 ::  t1894;
  real*8 ::  t1944;
  real*8 ::  t195;
  real*8 ::  t1956;
  real*8 ::  t1988;
  real*8 ::  t1997;
  real*8 ::  t2038;
  real*8 ::  t225;
  real*8 ::  t227;
  real*8 ::  t228;
  real*8 ::  t230;
  real*8 ::  t235;
  real*8 ::  t240;
  real*8 ::  t241;
  real*8 ::  t243;
  real*8 ::  t246;
  real*8 ::  t250;
  real*8 ::  t253;
  real*8 ::  t254;
  real*8 ::  t256;
  real*8 ::  t257;
  real*8 ::  t260;
  real*8 ::  t261;
  real*8 ::  t262;
  real*8 ::  t263;
  real*8 ::  t276;
  real*8 ::  t280;
  real*8 ::  t284;
  real*8 ::  t288;
  real*8 ::  t291;
  real*8 ::  t292;
  real*8 ::  t294;
  real*8 ::  t295;
  real*8 ::  t299;
  real*8 ::  t307;
  real*8 ::  t308;
  real*8 ::  t312;
  real*8 ::  t314;
  real*8 ::  t319;
  real*8 ::  t328;
  real*8 ::  t331;
  real*8 ::  t334;
  real*8 ::  t345;
  real*8 ::  t347;
  real*8 ::  t349;
  real*8 ::  t352;
  real*8 ::  t353;
  real*8 ::  t356;
  real*8 ::  t357;
  real*8 ::  t358;
  real*8 ::  t36;
  real*8 ::  t363;
  real*8 ::  t369;
  real*8 ::  t381;
  real*8 ::  t403;
  real*8 ::  t405;
  real*8 ::  t409;
  real*8 ::  t410;
  real*8 ::  t428;
  real*8 ::  t433;
  real*8 ::  t442;
  real*8 ::  t45;
  real*8 ::  t456;
  real*8 ::  t481;
  real*8 ::  t539;
  real*8 ::  t540;
  real*8 ::  t570;
  real*8 ::  t58;
  real*8 ::  t613;
  real*8 ::  t615;
  real*8 ::  t619;
  real*8 ::  t629;
  real*8 ::  t630;
  real*8 ::  t633;
  real*8 ::  t656;
  real*8 ::  t669;
  real*8 ::  t682;
  real*8 ::  t685;
  real*8 ::  t732;
  real*8 ::  t737;
  real*8 ::  t747;
  real*8 ::  t764;
  real*8 ::  t765;
  real*8 ::  t772;
  real*8 ::  t776;
  real*8 ::  t78;
  real*8 ::  t781;
  real*8 ::  t79;
  real*8 ::  t790;
  real*8 ::  t807;
  real*8 ::  t809;
  real*8 ::  t819;
  real*8 ::  t820;
  real*8 ::  t821;
  real*8 ::  t825;
  real*8 ::  t826;
  real*8 ::  t829;
  real*8 ::  t83;
  real*8 ::  t830;
  real*8 ::  t833;
  real*8 ::  t851;
  real*8 ::  t853;
  real*8 ::  t859;
  real*8 ::  t864;
  real*8 ::  t869;
  real*8 ::  t87;
  real*8 ::  t872;
  real*8 ::  t875;
  real*8 ::  t883;
  real*8 ::  t886;
  real*8 ::  t887;
  real*8 ::  t891;
  real*8 ::  t896;
  real*8 ::  t910;
  real*8 ::  t911;
  real*8 ::  t924;
  real*8 ::  t937;
  real*8 ::  t940;
  real*8 ::  t95;
  real*8 ::  t952;
  real*8 ::  t961;
  real*8 ::  t971;
  real*8 ::  t981;
  real*8 ::  t988;
  real*8 ::  t994;
  real*8 ::  t997;
  real*8 ::  t999;

    t1 = g01*g01;
    t3 = g22*g22;
    t4 = g01*t3;
    t5 = g33*g33;
    t6 = dgx01*dgx01;
    t7 = t5*t6;
    t11 = g01*g22;
    t12 = g23*g23;
    t13 = g33*t12;
    t17 = t12*t12;
    t18 = t1*t17;
    t19 = r*r;
    t20 = t19*t19;
    t22 = ddgxr02*t20*Rmin;
    t25 = t19*r;
    t27 = ddgxr02*t25*Rmin;
    t31 = ddgxr02*t19*Rmin;
    t35 = Theta22*t19*Rmin;
    t39 = dgx02*t19*Rmin;
    t42 = t1*t3;
    t43 = t5*ddgxx01;
    t48 = Theta22*r*Rmin;
    t52 = dgx02*r*Rmin;
    t55 = t1*g01;
    t56 = t3*t55;
    t57 = dgy33*dgx23;
    t61 = dgy33*dgy22;
    t70 = g33*ddgxy23;
    t77 = g33*ddgyy22;
    t84 = g33*ddgxx33;
    t88 = -4.0*t4*t7*r-4.0*t11*t13*t6+4.0*t18*t22-8.0*t18*t27+4.0*t18*t31+4.0*t18*t35 &
          -8.0*t18*t39-4.0*t42*t43*t19-4.0*t18*t48+8.0*t18*t52-2.0*t56*t57*t19        &
          +t56*t61*t19+4.0*t56*t57*r-2.0*t56*t61*r+4.0*t56*t70*t19-8.0*t56*t70*r      &
          -2.0*t56*t77*t19+4.0*t56*t77*r-2.0*t56*t84*t19;
    t92 = g22*t55;
    t93 = t12*ddgxy23;
    t100 = t12*ddgyy22;
    t104 = t12*ddgxx33;
    t111 = dgy22*dgy22;
    t112 = g33*t111;
    t125 = g23*dgx33;
    t142 = t12*g23;
    t143 = t1*t142;
    t144 = dgx01*dgx23;
    t148 = dgx01*dgy22;
    t152 = dgy01*dgx22;
    t162 = 4.0*t56*t84*r-4.0*t92*t93*t19+8.0*t92*t93*r-4.0*t92*t100*r+2.0*t92*t104*t19             &
           -4.0*t92*t104*r+t92*t112*t19-2.0*t92*t112*r+t92*g33*dgx33*dgx22-2.0*t92*g33*dgy23*dgx22 &
           -t92*t125*dgy22+t92*g23*dgy33*dgx22-2.0*t92*t125*dgx23+4.0*t92*g23*dgx23*dgy23          &
           -2.0*t92*g23*dgy23*dgy22+4.0*t143*t144*t19-2.0*t143*t148*t19+2.0*t143*t152*t19+8.0*t42*t43*r-8.0*t143*t144*r;
    t170 = g33*dgy01;
    t177 = t1*g22;
    t185 = t12*dgy01;
    t192 = t1*g33;
    t197 = g01*t17;
    t198 = dg02*dg02;
    t199 = t20*t19;
    t201 = Rmin*Rmin;
    t202 = t198*t199*t201;
    t205 = t20*r;
    t207 = t198*t205*t201;
    t211 = t198*t20*t201;
    t234 = 4.0*t143*t148*r-4.0*t143*t152*r+4.0*t42*t170*dgx23-2.0*t42*t170*dgy22                 &
           +2.0*t177*t5*dgx01*dgx22+8.0*t177*t13*ddgxx01-4.0*t177*t185*dgx23+2.0*t177*t185*dgy22 &
           -2.0*t192*t12*dgx01*dgx22-2.0*t197*t202+4.0*t197*t207-2.0*t197*t211+2.0*t4*t7*t19     &
           +2.0*t92*t100*t19-4.0*t18*ddgxx01*t19+8.0*t18*ddgxx01*r                               &
           -4.0*t42*t43+4.0*t143*t144-2.0*t143*t148+2.0*t143*t152;
    t249 = dgx33*dgx33;
    t265 = t201*t20;
    t266 = t265*t11;
    t267 = g03*dg33;
    t270 = t267*dg22*g23*g02;
    t273 = t201*t199;
    t274 = t273*t11;
    t277 = g33*g23;
    t278 = t11*t277;
    t279 = g02*dg02;
    t285 = g02*dg22;
    t301 = 2.0*t197*t6*t19-4.0*t197*t6*r+2.0*t4*t7+4.0*t56*t70-2.0*t56*t77-2.0*t56*t84+t56*t249*t19         &
           -2.0*t56*t249*r-2.0*t56*t57+t56*t61-4.0*t92*t93+2.0*t92*t100+2.0*t92*t104+t92*t112-2.0*t266*t270 &
           -2.0*t274*t270-4.0*t278*t279*dg23*t199*t201+4.0*t278*t285*dg03*t199*t201                         &
           +8.0*t278*t279*dg23*t205*t201-8.0*t278*t285*dg03*t205*t201;
    t304 = g02*g03;
    t305 = dg22*t20;
    t306 = t305*t201;
    t315 = dg03*t20;
    t320 = dg22*t25;
    t321 = t320*t201;
    t325 = t205*t201;
    t326 = t325*t11;
    t327 = g00*dg23;
    t329 = t277*t327*dg22;
    t333 = t277*t304*ddg22;
    t336 = g22*g33;
    t337 = t325*t336;
    t338 = g23*g02;
    t339 = g03*dg01;
    t341 = t338*t339*dg22;
    t348 = t265*t336;
    t355 = t273*t336;
    t364 = t273*g22;
    t365 = g02*g02;
    t366 = t5*t365;
    t367 = dg01*dg22;
    t368 = t366*t367;
    t371 = g03*g03;
    t372 = t12*t371;
    t373 = t372*t367;
    t377 = t12*t365;
    t378 = t377*t367;
    t382 = t304*t367;
    t385 = -20.0*t278*t304*t306-4.0*t278*t279*dg23*t20*t201+4.0*t278*t285*t315*t201+12.0*t278*t304*t321                     &
           -4.0*t326*t329-8.0*t326*t333+8.0*t337*t341+2.0*t266*t329+4.0*t266*t333-4.0*t348*t341+2.0*t274*t329+4.0*t274*t333 &
           -4.0*t355*t341+8.0*t326*t277*t304*dg22+4.0*t326*t270+2.0*t364*t368-2.0*t364*t373-2.0*t273*g33*t378+4.0*t273*t142*t382;
    t386 = t177*g33;
    t387 = t12*dgx02;
    t388 = r*Rmin;
    t392 = g23*dgx01;
    t393 = dgx23*t19;
    t397 = dgy22*t19;
    t401 = g23*dgy01;
    t402 = dgx22*t19;
    t406 = t177*t12;
    t407 = g03*dgx23;
    t408 = t407*t388;
    t411 = g03*dgy22;
    t412 = t411*t388;
    t415 = t192*t12;
    t416 = g02*dgx22;
    t417 = t416*t388;
    t420 = t4*g33;
    t421 = t12*dg00;
    t422 = t25*t201;
    t426 = t371*dg22;
    t427 = t426*t422;
    t430 = t4*t12;
    t431 = g03*dg03;
    t432 = t431*t422;
    t435 = t11*t142;
    t436 = g02*dg03;
    t437 = t436*t422;
    t440 = g03*dg02;
    t441 = t440*t422;
    t444 = t11*t12;
    t448 = g01*t142*g02;
    t449 = g03*dg22;
    t453 = t42*g33;
    t454 = t19*Rmin;
    t455 = t407*t454;
    t458 = t411*t454;
    t461 = Theta23*dg23;
    t462 = t461*t454;
    t465 = dg23*dgy02;
    t466 = t465*t454;
    t469 = dgx23*dg03;
    t470 = t469*t454;
    t473 = dgy22*dg03;
    t474 = t473*t454;
    t477 = -20.0*t386*t387*t388-4.0*t386*t392*t393+2.0*t386*t392*t397-2.0*t386*t401*t402+8.0*t406*t408 &
           -8.0*t406*t412+4.0*t415*t417-8.0*t420*t421*t422-10.0*t420*t427                              &
           +8.0*t430*t432-8.0*t435*t437-8.0*t435*t441+12.0*t444*t427                                   &
           -16.0*t448*t449*t422+12.0*t453*t455-6.0*t453*t458+4.0*t453*t462-4.0*t453*t466-4.0*t453*t470+2.0*t453*t474;
    t479 = t177*t5;
    t480 = t416*t454;
    t483 = dg02*dgx22;
    t484 = t483*t454;
    t487 = Theta22*dg22;
    t488 = t487*t454;
    t491 = t12*ddgxr02;
    t495 = t11*g33;
    t496 = t12*t198;
    t500 = dg02*dg23;
    t501 = t500*t273;
    t504 = dg22*dg03;
    t505 = t504*t273;
    t508 = t4*t5;
    t509 = dg00*dg22;
    t510 = t509*t325;
    t516 = t500*t325;
    t519 = t504*t325;
    t522 = t3*g22;
    t523 = g01*t522;
    t524 = t523*g33;
    t525 = t431*t265;
    t528 = g00*dg22;
    t529 = t528*t265;
    t532 = t279*t265;
    t535 = t509*t265;
    t541 = t426*t265;
    t549 = t436*t265;
    t552 = 6.0*t479*t480-2.0*t479*t484+2.0*t479*t488-8.0*t386*t491*t454+4.0*t495*t496*t273+4.0*t448*t501               &
           -4.0*t448*t505-4.0*t508*t510-8.0*t495*t496*t325-8.0*t448*t516+8.0*t448*t519+8.0*t524*t525-14.0*t508*t529    &
           +8.0*t508*t532+2.0*t508*t535+8.0*t420*t421*t265+14.0*t420*t541-8.0*t430*t525+4.0*t495*t496*t265+8.0*t435*t549;
    t553 = t440*t265;
    t561 = t500*t265;
    t564 = t504*t265;
    t569 = t528*t422;
    t572 = t279*t422;
    t575 = Rmin*t25;
    t576 = t575*t1;
    t577 = t3*g33;
    t579 = t577*ddgxr23*g03;
    t582 = t454*t1;
    t586 = t12*dgx33*t285;
    t589 = t12*dg23;
    t590 = dgx23*g02;
    t591 = t589*t590;
    t594 = dgy22*g02;
    t595 = t589*t594;
    t605 = t12*dg33*t416;
    t610 = t325*g01;
    t611 = t3*t5;
    t618 = 8.0*t435*t553-16.0*t444*t541+24.0*t448*t449*t265+4.0*t448*t561-4.0*t448*t564-8.0*t524*t432 &
           +10.0*t508*t569-8.0*t508*t572+8.0*t576*t579-4.0*t582*t579+4.0*t576*t586                    &
           -8.0*t576*t591+8.0*t576*t595-2.0*t582*t586+4.0*t582*t591                                   &
           -4.0*t582*t595-4.0*t576*t605+2.0*t582*t605+4.0*t610*t611*t528-4.0*t610*t577*t426;
    t622 = g22*t5;
    t623 = t365*dg22;
    t624 = t622*t623;
    t627 = g22*t12;
    t631 = t13*t623;
    t634 = t142*g02;
    t638 = t273*g01;
    t639 = g00*ddg22;
    t640 = t611*t639;
    t643 = dg23*dg23;
    t644 = g00*t643;
    t645 = t577*t644;
    t648 = t371*ddg22;
    t649 = t577*t648;
    t658 = dg22*dg22;
    t659 = g00*t658;
    t660 = t622*t659;
    t663 = t365*ddg22;
    t664 = t622*t663;
    t668 = t336*t365*t643;
    t671 = t12*Theta22;
    t686 = -4.0*t610*t624+4.0*t610*t627*t426+4.0*t610*t631-8.0*t610*t634*t449+2.0*t638*t640-2.0*t638*t645-2.0*t638*t649 &
           -4.0*t610*t640+4.0*t610*t645+4.0*t610*t649+2.0*t610*t660+4.0*t610*t664                                       &
           -4.0*t610*t668-10.0*t386*t671*t454+20.0*t386*t387*t454                                                       &
           -8.0*t406*t455+8.0*t406*t458-4.0*t406*t462+4.0*t406*t466;
    t706 = t19*t201;
    t710 = t304*t706;
    t713 = t12*g00;
    t717 = t265*g01;
    t723 = t336*t371*t658;
    t725 = t627*t644;
    t728 = t627*t648;
    t731 = t13*t659;
    t733 = t13*t663;
    t738 = 4.0*t406*t470-2.0*t406*t474-4.0*t415*t480+2.0*t415*t484-2.0*t415*t488-12.0*t453*t408+6.0*t453*t412 &
           -6.0*t479*t417+10.0*t386*t671*t388+12.0*t495*t377*t706                                             &
           -24.0*t435*t710-24.0*t420*t713*t706+10.0*t717*t624+2.0*t717*t668                                   &
           +t717*t723+2.0*t717*t725+2.0*t717*t728+t717*t731+2.0*t717*t733-12.0*t717*t631;
    t741 = dg23*dg22;
    t742 = t142*g00*t741;
    t746 = t634*g03*ddg22;
    t749 = t265*t3;
    t750 = t5*g00;
    t751 = t750*t367;
    t754 = g33*t371;
    t755 = t754*t367;
    t758 = t265*g22;
    t775 = t273*t3;
    t780 = t4*t25;
    t782 = g00*dg33;
    t783 = t201*t12*t782;
    t786 = t177*Rmin;
    t787 = t19*t12;
    t788 = g02*dgx33;
    t792 = r*t12;
    t796 = dg33*Theta22;
    t799 = t20*t12;
    t800 = dg33*dgx02;
    t804 = t25*t12;
    t811 = -2.0*t717*t742-4.0*t717*t746-2.0*t749*t751+2.0*t749*t755+2.0*t758*t368-2.0*t758*t373   &
           -2.0*t265*g33*t378+4.0*t265*t142*t382+2.0*t638*t733-2.0*t638*t742                      &
           -4.0*t638*t746-2.0*t775*t751+2.0*t775*t755-2.0*t780*t783+4.0*t786*t787*t788            &
           -4.0*t786*t792*t788+t786*t787*t796-2.0*t786*t799*t800+4.0*t786*t804*t800-2.0*t786*t787*t800;
    t812 = dg22*dgy03;
    t822 = Theta33*dg22;
    t831 = t201*g33;
    t832 = t831*t782;
    t835 = t4*t20;
    t836 = t365*dg33;
    t837 = t831*t836;
    t850 = t422*g01;
    t855 = t422*t3;
    t856 = g00*dg01;
    t857 = t13*t856;
    t860 = t422*g22;
    t862 = t13*t365*dg01;
    t865 = t634*t339;
    t871 = -2.0*t786*t799*t812+4.0*t786*t804*t812-2.0*t786*t787*t812+t786*t799*t822-2.0*t786*t804*t822+t786*t787*t822           &
           -2.0*t523*t20*t832+2.0*t835*t837+2.0*t835*t783-2.0*t610*t723-4.0*t610*t725-4.0*t610*t728-2.0*t610*t731-6.0*t850*t624 &
           +8.0*t850*t631+8.0*t855*t857-4.0*t860*t862+8.0*t860*t865-t638*t660-2.0*t638*t664;
    t876 = dgx23*r;
    t880 = dgy22*r;
    t884 = dgx22*r;
    t888 = t20*Rmin;
    t889 = t461*t888;
    t892 = t465*t888;
    t895 = t469*t888;
    t898 = t473*t888;
    t901 = t483*t888;
    t904 = t487*t888;
    t922 = t461*t575;
    t925 = t465*t575;
    t928 = 2.0*t638*t668+8.0*t386*t392*t876-4.0*t386*t392*t880+4.0*t386*t401*t884+4.0*t453*t889-4.0*t453*t892 &
           -4.0*t453*t895+2.0*t453*t898-2.0*t479*t901+2.0*t479*t904                                           &
           -8.0*t386*t491*t888-4.0*t406*t889+4.0*t406*t892+4.0*t406*t895                                      &
           -2.0*t406*t898+2.0*t415*t901-2.0*t415*t904-8.0*t453*t922+8.0*t453*t925;
    t929 = t469*t575;
    t932 = t473*t575;
    t935 = t483*t575;
    t938 = t487*t575;
    t956 = t509*t273;
    t959 = t3*dg33;
    t960 = t959*t407;
    t963 = t959*t411;
    t970 = t577*ddgyr22*g03;
    t975 = t888*t1;
    t979 = dg33*dg22;
    t980 = t3*t371*t979;
    t983 = 8.0*t453*t929-4.0*t453*t932+4.0*t479*t935-4.0*t479*t938+16.0*t386*t491*t575+8.0*t406*t922-8.0*t406*t925          &
           -8.0*t406*t929+4.0*t406*t932-4.0*t415*t935+4.0*t415*t938+2.0*t508*t956-4.0*t576*t960+2.0*t576*t963+2.0*t582*t960 &
           -t582*t963-8.0*t576*t970+4.0*t582*t970+2.0*t975*t605-2.0*t610*t980;
    t985 = t377*t979;
    t991 = t627*t836;
    t996 = t3*g23*t371*dg23;
    t1001 = t201*t142*t327;
    t1009 = t42*Rmin;
    t1010 = r*g23;
    t1011 = g02*dgy33;
    t1015 = t20*g33;
    t1019 = t25*g33;
    t1023 = t19*g33;
    t1037 = g02*dgy23;
    t1041 = t19*g23;
    t1042 = g03*dgx33;
    t1049 = -4.0*t610*t985+t717*t980+2.0*t717*t985+4.0*t850*t991-4.0*t850*t996+4.0*t11*t25*t1001+t786*t799*t796 &
            -2.0*t786*t804*t796+2.0*t1009*t1010*t1011+2.0*t1009*t1015*t812                                      &
            -4.0*t1009*t1019*t812+2.0*t1009*t1023*t812-t1009*t1015*t822                                         &
            +2.0*t1009*t1019*t822-t1009*t1023*t822-2.0*t1009*t1023*t788                                         &
            +4.0*t1009*t1023*t1037-2.0*t1009*t1041*t1042+t638*t723+2.0*t638*t725;
    t1065 = t325*t3;
    t1070 = t325*g22;
    t1090 = g03*dgy23;
    t1094 = 2.0*t638*t728+t638*t731-8.0*t749*t857+4.0*t758*t862-8.0*t758*t865-4.0*t610*t733 &
            +4.0*t610*t742+8.0*t610*t746+4.0*t1065*t751                                     &
            -4.0*t1065*t755-4.0*t1070*t368+4.0*t1070*t373+4.0*t325*g33*t378                 &
            -8.0*t325*t142*t382+2.0*t717*t640-2.0*t717*t645-2.0*t717*t649                   &
            -t717*t660-2.0*t717*t664-4.0*t1009*t1041*t1090;
    t1099 = r*g33;
    t1113 = dgy33*dg22;
    t1123 = t19*g03;
    t1126 = g23*Theta23;
    t1130 = g23*dgx03;
    t1134 = g23*dgy02;
    t1160 = 2.0*t1009*t1099*t788-4.0*t1009*t1099*t1037+2.0*t1009*t1010*t1042                                           &
            +4.0*t1009*t1010*t1090-t1009*t20*g03*t1113+2.0*t1009*t25*g03*t1113                                         &
            -2.0*t1009*t1041*t1011-t1009*t1123*t1113-4.0*t1009*t1023*t1126+4.0*t1009*t1023*t1130+4.0*t1009*t1023*t1134 &    
            +2.0*t1009*t1019*t796-t1009*t1023*t796+2.0*t1009*t1015*t800                                                &
            -4.0*t1009*t1019*t800+2.0*t1009*t1023*t800-t1009*t1015*t796+4.0*t1009*t1099*t1126-4.0*t1009*t1099*t1130;
    t1173 = g22*g03;
    t1174 = t12*ddgyr22*t1173;
    t1180 = t12*ddgxr23*t1173;
    t1207 = -4.0*t1009*t1099*t1134+2.0*t523*t25*t832-2.0*t780*t837-4.0*t11*t20*t1001                             &
            +8.0*t576*t1174-4.0*t582*t1174-8.0*t576*t1180+4.0*t582*t1180-4.0*t975*t1174+4.0*t975*t1180+t638*t980 &
            +2.0*t638*t985-4.0*t717*t991+4.0*t717*t996+2.0*t975*t960-t975*t963                                   &
            +4.0*t975*t970-4.0*t975*t579-2.0*t975*t586+4.0*t975*t591;
    t1211 = t143*g02;
    t1218 = t143*g03;
    t1222 = dgy01*dgx23;
    t1223 = t1222*r;
    t1226 = dgy01*dgy22;
    t1227 = t1226*r;
    t1230 = dgx01*dgx22;
    t1231 = t1230*r;
    t1234 = t12*ddgxx01;
    t1240 = t11*t17;
    t1242 = dg00*t25*t201;
    t1245 = t197*g00;
    t1248 = t523*t5;
    t1250 = g00*t19*t201;
    t1254 = t371*t19*t201;
    t1265 = t12*t6;
    t1272 = t143*dg02;
    t1277 = t143*Theta23;
    t1281 = -4.0*t975*t595-8.0*t1211*t876*Rmin+4.0*t1211*t880*Rmin-4.0*t1218*t884*Rmin                      &
            -8.0*t453*t1223+4.0*t453*t1227-4.0*t479*t1231-16.0*t386*t1234*r+8.0*t406*t1223+4.0*t1240*t1242  &
            +8.0*t1245*t321+12.0*t1248*t1250-12.0*t524*t1254-12.0*t508*t365*t19*t201                        &
            +12.0*t430*t1254+12.0*t1240*t1250-4.0*t495*t1265*t19+8.0*t495*t1265*r-4.0*t1272*dgy22*t25*Rmin-4.0*t1277*t320*Rmin;
    t1282 = t143*dg23;
    t1287 = t143*dg22;
    t1296 = t143*dgx22;
    t1301 = t42*t5;
    t1308 = t393*Rmin;
    t1311 = t397*Rmin;
    t1325 = t92*g23;
    t1326 = dgx23*dgy23;
    t1330 = dgy23*dgy22;
    t1334 = t92*g33;
    t1335 = dgx33*dgx22;
    t1338 = dgy23*dgx22;
    t1348 = dgx33*dgy22;
    t1351 = -4.0*t1282*Theta22*t25*Rmin-4.0*t1287*dgx03*t25*Rmin+4.0*t1287*dgy02*t25*Rmin                &
            +4.0*t1296*dg03*t25*Rmin+4.0*t1301*t31+6.0*t1301*t35-12.0*t1301*t39+8.0*t1211*t1308          &
            -4.0*t1211*t1311+4.0*t1218*t402*Rmin-4.0*t1272*t1308                                         &
            +2.0*t1272*t1311+2.0*t1277*dg22*t19*Rmin-8.0*t1325*t1326*r+4.0*t1325*t1330*r+t1334*t1335*t19 &
            -2.0*t1334*t1338*t19-2.0*t1334*t1335*r+4.0*t1334*t1338*r-t1325*t1348*t19;
    t1354 = dgy33*dgx22;
    t1386 = t197*dg00;
    t1398 = dg00*t20*t201;
    t1411 = t1325*t1354*t19+2.0*t1325*t1348*r-2.0*t1325*t1354*r+2.0*t1282*Theta22*t20*Rmin                          &
            +2.0*t1287*dgx03*t20*Rmin-2.0*t1287*dgy02*t20*Rmin-2.0*t1296*t315*Rmin-8.0*t1301*t27                    &
            +8.0*t1272*dgx23*t25*Rmin-2.0*t508*t202+2.0*t1386*dg22*t199*t201+4.0*t508*t207-4.0*t1386*dg22*t205*t201 &
            -4.0*t1248*t1398-2.0*t508*t211-4.0*t1240*t1398-12.0*t1245*t306+2.0*t1386*t306+4.0*t1248*t1242;
    t1426 = t142*ddgyr22*g02;
    t1432 = t142*ddgxr23*g02;
    t1437 = t377*t643;
    t1441 = t522*t371*dg33;
    t1449 = t522*t1*Rmin;
    t1475 = 4.0*t1301*t22-4.0*t1272*dgx23*t20*Rmin+2.0*t1272*dgy22*t20*Rmin+2.0*t1277*t305*Rmin         &
            -8.0*t576*t1426+4.0*t582*t1426+8.0*t576*t1432-4.0*t582*t1432-2.0*t638*t1437-2.0*t717*t1441  &
            +4.0*t975*t1426-4.0*t975*t1432+2.0*t1449*t1123*dgy33                                        &
            -2.0*t1449*r*g03*dgy33-4.0*t1449*t1023*dgy03+4.0*t1009*t787*dgy03                           &
            +2.0*t1449*t1023*Theta33-2.0*t1009*t787*Theta33+4.0*t1449*t1099*dgy03-4.0*t1009*t792*dgy03;
    t1483 = dgx33*dgx23;
    t1496 = r*t142;
    t1506 = t19*t142;
    t1522 = t422*t522;
    t1523 = t750*dg01;
    t1526 = t754*dg01;
    t1529 = t366*dg01;
    t1532 = t372*dg01;
    t1535 = t17*g00;
    t1536 = t1535*dg01;
    t1539 = -2.0*t1449*t1099*Theta33+2.0*t1009*t792*Theta33+4.0*t1325*t1483*r-2.0*t1325*t1483*t19+4.0*t1325*t1326*t19    &
            -2.0*t1325*t1330*t19-4.0*t786*t1496*Theta23+4.0*t786*t1496*dgx03+4.0*t786*t1496*dgy02+4.0*t786*t1506*Theta23 &
            -4.0*t786*t1506*dgx03-4.0*t786*t1506*dgy02+4.0*t610*t1437-2.0*t717*t1437+2.0*t850*t1441                      &
            -4.0*t1522*t1523+4.0*t1522*t1526+4.0*t855*t1529-4.0*t855*t1532-4.0*t860*t1536;
    t1540 = t366*t658;
    t1543 = t1535*ddg22;
    t1547 = t856*dg22;
    t1556 = t265*t522;
    t1592 = -2.0*t610*t1540-4.0*t610*t1543+4.0*t325*t17*t1547+t717*t1540+2.0*t717*t1543        &
            -2.0*t265*t17*t1547+4.0*t1556*t1523-4.0*t1556*t1526-4.0*t749*t1529+4.0*t749*t1532  &
            +4.0*t758*t1536+t638*t1540+2.0*t638*t1543-2.0*t273*t17*t1547                       &
            +4.0*t610*t1535*dg22+2.0*t1282*t35+2.0*t1287*dgx03*t19*Rmin                        &
            -2.0*t1287*dgy02*t19*Rmin-2.0*t1296*dg03*t19*Rmin-6.0*t1301*t48;
    t1598 = t1222*t19;
    t1601 = t1226*t19;
    t1604 = t1230*t19;
    t1629 = t11*t13;
    t1636 = t365*dg23;
    t1641 = t11*t12*g03;
    t1646 = 12.0*t1301*t52+4.0*t453*t1598-2.0*t453*t1601+2.0*t479*t1604+8.0*t386*t1234*t19-4.0*t406*t1598  &
            +2.0*t406*t1601-2.0*t415*t1604-4.0*t406*t1227+4.0*t415*t1231-4.0*t386*t392*dgx23               &
            +2.0*t386*t392*dgy22-2.0*t386*t401*dgx22                                                       &
            +26.0*t1629*t529-8.0*t1629*t532-4.0*t1629*t535+4.0*t278*t1636*t265-4.0*t1641*t561+4.0*t1641*t564;
    t1647 = t4*t277;
    t1652 = g33*g02;
    t1653 = t4*t1652;
    t1654 = g03*dg23;
    t1668 = t4*g33*g03;
    t1673 = t177*t277;
    t1674 = g03*dgx22;
    t1678 = dg02*dgx23;
    t1682 = dg02*dgy22;
    t1686 = Theta23*dg22;
    t1691 = g23*t371*t741;
    t1694 = g01*g33;
    t1695 = t325*t1694;
    t1697 = g23*t365*t741;
    t1700 = t1*g23;
    t1701 = t454*t1700;
    t1703 = dg23*dgx23*t1173;
    t1706 = dg23*dgy22;
    t1707 = t1706*t1173;
    t1710 = t575*t1700;
    t1711 = dg33*dgx22;
    t1712 = t1711*t1173;
    t1716 = t888*t1700;
    t1717 = dgx33*dg22;
    t1718 = t1717*t1173;
    t1720 = 8.0*t1647*t437+8.0*t1647*t441+8.0*t1653*t1654*t422-18.0*t1629*t569+8.0*t1629*t572  &
            -4.0*t278*t1636*t422+24.0*t1647*t710+4.0*t1668*t501                                &
            -4.0*t1668*t505+6.0*t1673*t1674*t388+4.0*t1673*t1678*t888-2.0*t1673*t1682*t888     &
            -2.0*t1673*t1686*t888+4.0*t326*t1691+4.0*t1695*t1697-4.0*t1701*t1703               &
            +2.0*t1701*t1707+2.0*t1710*t1712-t1701*t1712+t1716*t1718;
    t1727 = t265*t4;
    t1728 = t338*t267;
    t1731 = t888*t177;
    t1733 = dg33*dgx23*t338;
    t1737 = dg33*dgy22*t338;
    t1740 = g33*ddgyr22*t338;
    t1744 = g33*ddgxr23*t338;
    t1747 = dgx22*dg03;
    t1760 = dg23*Theta22;
    t1764 = dg22*dgx03;
    t1768 = dg22*dgy02;
    t1787 = -4.0*t1716*t1703+2.0*t1716*t1707-t1716*t1712+4.0*t1727*t1728-2.0*t1731*t1733+t1731*t1737 &
            -4.0*t1731*t1740+4.0*t1731*t1744+2.0*t1673*t1747*t888-8.0*t1673*t1678*t575               &
            +4.0*t1673*t1682*t575+4.0*t1673*t1686*t575+4.0*t1673*t1760*t575                          &
            +4.0*t1673*t1764*t575-4.0*t1673*t1768*t575-4.0*t1673*t1747*t575-12.0*t1673*t590*t454     &
            +2.0*t1673*t594*t454-6.0*t1673*t1674*t454+4.0*t1673*t1678*t454;
    t1813 = t338*t339;
    t1817 = t338*g03*t658;
    t1820 = g01*t12;
    t1822 = t304*t741;
    t1825 = t713*t367;
    t1828 = t13*t639;
    t1833 = t265*t1694;
    t1847 = -2.0*t1673*t1682*t454-2.0*t1673*t1686*t454-2.0*t1673*t1760*t454-2.0*t1673*t1764*t454 &
            +2.0*t1673*t1768*t454+2.0*t1673*t1747*t454+12.0*t1673*t590*t388-2.0*t1673*t594*t388  &
            -8.0*t422*t577*t1813+4.0*t1695*t1817-8.0*t325*t1820*t1822                            &
            -8.0*t337*t1825-4.0*t266*t1828-2.0*t266*t1691-2.0*t1833*t1697-2.0*t1833*t1817        &
            +4.0*t265*t1820*t1822+4.0*t348*t1825-4.0*t1629*t956-4.0*t1641*t501;
    t1873 = t713*t979;
    t1876 = g33*t365*t979;
    t1882 = t1652*t1711;
    t1884 = t177*t575;
    t1887 = t177*t454;
    t1891 = g33*g00*t979;
    t1896 = 4.0*t1641*t505-8.0*t1668*t516+8.0*t1668*t519+8.0*t1629*t510+8.0*t1641*t516-8.0*t1641*t519-8.0*t1647*t549 &
            -8.0*t1647*t553-8.0*t1653*t1654*t265+4.0*t1668*t561-4.0*t1668*t564-t274*t1873                            &
            +2.0*t326*t1876+2.0*t326*t1873-t266*t1876                                                                &
            -t1731*t1882+2.0*t1884*t1882-t1887*t1882+t4*t273*t1891-2.0*t4*t325*t1891;
    t1897 = t277*t327;
    t1901 = t4*t422;
    t1904 = t1652*t1717;
    t1906 = t1652*t1706;
    t1909 = dgy23*dg22;
    t1910 = t1652*t1909;
    t1914 = g23*g03*t1909;
    t1932 = t338*t1113;
    t1934 = t888*t192;
    t1935 = dg23*dgx22;
    t1936 = t1935*t1173;
    t1939 = t1935*t338;
    t1942 = dg22*dgy22;
    t1943 = t1942*t1173;
    t1946 = t1942*t338;
    t1949 = 4.0*t1727*t1897+t1727*t1891-4.0*t1901*t1897+t1731*t1904+2.0*t1731*t1906                      &
            -2.0*t1731*t1910+2.0*t1731*t1914-2.0*t1884*t1904-4.0*t1884*t1906                             &
            +4.0*t1884*t1910-4.0*t1884*t1914+t1887*t1904+2.0*t1887*t1906-2.0*t1887*t1910                 &
            +2.0*t1887*t1914+t1731*t1932+2.0*t1934*t1936-2.0*t1934*t1939-2.0*t1934*t1943+2.0*t1934*t1946;
    t1973 = t575*t192;
    t1982 = t454*t192;
    t1995 = -4.0*t1901*t1728+4.0*t1884*t1733-2.0*t1884*t1737-2.0*t1887*t1733+t1887*t1737             &
            +8.0*t1884*t1740-4.0*t1887*t1740-2.0*t1673*t1760*t888-2.0*t1673*t1764*t888               &
            +2.0*t1673*t1768*t888-4.0*t1973*t1936+4.0*t1973*t1939+4.0*t1973*t1943-4.0*t1973*t1946    &
            +2.0*t1982*t1936-2.0*t1982*t1939-2.0*t1982*t1943-t266*t1873+8.0*t265*t577*t1813-4.0*t274*t1828;
    t1998 = t273*t1694;
    t2035 = -2.0*t274*t1691-2.0*t1998*t1697+2.0*t1982*t1946-2.0*t1710*t1718+8.0*t1710*t1703       &
            -4.0*t1710*t1707+t1701*t1718-8.0*t1884*t1744+4.0*t1887*t1744                          &
            -2.0*t1998*t1817+4.0*t273*t1820*t1822+4.0*t355*t1825-8.0*t326*t13*t528+8.0*t326*t1828 &
            -2.0*t1884*t1932+t1887*t1932-t274*t1876+2.0*t197*t6-4.0*t18*ddgxx01+t56*t249;
    Theta22_rhs = 1/t1*(t811+t738+t686+t618+t552+t477+t1720+t1646+t1592+t1539+t1475+t2035          &
                        +t385+t1995+t1411+t1351+t1281+t301+t1207+t1160+t1949+t1094+t1049+t983+t928 &
                        +t871+t1896+t1847+t1787+t234+t162+t88)                                     &
                  /(-2.0*t627*g33-2.0*t17*r-2.0*t611*r+t17*t19+t611*t19-2.0*t627*t1023+4.0*t627*t1099+t17+t611)/Rmin/t19/4.0;

    t1 = g01*g01;
    t3 = g23*g23;
    t4 = t3*t3;
    t5 = t1*t4;
    t8 = t3*g23;
    t9 = t1*g01;
    t10 = t8*t9;
    t17 = g01*g22;
    t18 = t3*g33;
    t19 = t17*t18;
    t20 = dg00*dg23;
    t21 = r*r;
    t22 = t21*t21;
    t23 = Rmin*Rmin;
    t24 = t22*t23;
    t25 = t20*t24;
    t28 = dg02*dg03;
    t29 = t28*t24;
    t32 = g33*g23;
    t33 = t17*t32;
    t34 = g03*g03;
    t35 = t34*dg22;
    t39 = g01*g33;
    t40 = t3*g02;
    t41 = t39*t40;
    t42 = g03*dg22;
    t46 = t1*g22;
    t47 = t46*t32;
    t48 = dg23*dgx03;
    t49 = t21*r;
    t50 = t49*Rmin;
    t54 = dg23*dgy02;
    t59 = dgy22*dg03;
    t63 = g02*dgx33;
    t64 = t21*Rmin;
    t68 = g03*dgy22;
    t72 = dg02*dgx33;
    t76 = Theta33*dg22;
    t89 = r*Rmin;
    t96 = t22*Rmin;
    t100 = -2.0*t47*t48*t64-4.0*t47*t59*t50-2.0*t47*t54*t64+2.0*t47*t59*t64   &
           -6.0*t47*t63*t64+6.0*t47*t63*t89-6.0*t47*t68*t64+2.0*t47*t72*t64-2.0*t47*t76*t64+6.0*t47*t68*t89+2.0*t47*t72*t96;
    t120 = dg02*dg23;
    t121 = t120*t24;
    t124 = dg22*dg03;
    t125 = t124*t24;
    t128 = g22*g22;
    t129 = g01*t128;
    t130 = t129*t32;
    t131 = g03*dg03;
    t132 = t49*t23;
    t133 = t131*t132;
    t136 = g33*g33;
    t137 = t136*g23;
    t138 = t17*t137;
    t139 = g02*dg02;
    t140 = t139*t132;
    t144 = g00*dg23;
    t145 = t144*t132;
    t148 = g02*dg03;
    t152 = g03*dg02;
    t162 = g02*g03;
    t163 = t21*t23;
    t168 = t17*t136*g02;
    t169 = t22*t21;
    t170 = t169*t23;
    t171 = t120*t170;
    t174 = t24*t17;
    t175 = g33*ddg23;
    t176 = g00*t3;
    t177 = t175*t176;
    t181 = dg33*dg23;
    t183 = t181*g00*g33;
    t185 = g02*g02;
    t186 = t185*g33;
    t187 = t181*t186;
    t189 = t181*t176;
    t191 = -20.0*t19*t145+8.0*t19*t148*t132+8.0*t19*t152*t132-2.0*t33*t35*t132-4.0*t41*t42*t132 &
           +24.0*t19*t162*t163-2.0*t168*t171-4.0*t174*t177-t24*t129*t183+t174*t187+t174*t189;
    t194 = g01*t136;
    t196 = dg23*dg22;
    t197 = g00*g22;
    t198 = t196*t197;
    t200 = t24*t39;
    t201 = t34*g22;
    t202 = t196*t201;
    t204 = t196*t176;
    t206 = g01*t3;
    t208 = dg33*dg22;
    t209 = t208*t162;
    t212 = g22*g33;
    t213 = t24*t212;
    t214 = dg01*dg23;
    t215 = t214*t176;
    t219 = t3*dg01*t162;
    t222 = t170*t17;
    t226 = t170*t39;
    t233 = t170*t212;
    t236 = t22*r;
    t237 = t236*t23;
    t238 = t237*t17;
    t239 = g33*dg23;
    t247 = t1*t128;
    t248 = t136*ddgxy01;
    t251 = t1*t8;
    t252 = dgx01*dgx33;
    t255 = dgy01*dgy22;
    t258 = g01*t4;
    t259 = dgx01*dgy01;
    t268 = 2.0*t170*t206*t209+4.0*t233*t215-8.0*t238*t239*t176-4.0*t222*t177 &
           -t170*t129*t183-4.0*t247*t248+2.0*t251*t252+2.0*t251*t255+2.0*t258*t259-4.0*t10*ddgxy23*t21+8.0*t10*ddgxy23*r;
    t282 = t3*t9;
    t283 = dgx33*dgy22;
    t285 = dgy33*dgx22;
    t287 = dgx33*dgx23;
    t290 = dgx23*dgy23;
    t293 = dgy23*dgy22;
    t296 = g23*t9;
    t297 = dgy22*dgy22;
    t298 = g33*t297;
    t301 = dgx33*dgx33;
    t302 = g22*t301;
    t310 = dg33*t22;
    t311 = t310*t23;
    t315 = g02*dg23;
    t316 = dg03*t22;
    t317 = t316*t23;
    t321 = t237*t212;
    t322 = g23*g02;
    t323 = t322*g03;
    t324 = t214*t323;
    t327 = t175*t323;
    t330 = t181*t323;
    t333 = t196*t323;
    t336 = g01*g23;
    t337 = t24*t336;
    t339 = t208*t197*g33;
    t344 = t296*t302-4.0*t5*ddgxy01*t21+8.0*t5*ddgxy01*r-2.0*t33*t139*t311+2.0*t33*t315*t317    &
           +8.0*t321*t324+4.0*t174*t327-2.0*t174*t330-2.0*t200*t333+2.0*t337*t339-4.0*t213*t324;
    t366 = t237*t39;
    t373 = g03*dg33;
    t376 = t373*dg22*g33*g02;
    t384 = dg23*t169*t23;
    t389 = dg03*t169*t23;
    t394 = dg23*t236*t23;
    t399 = dg03*t236*t23;
    t404 = dg23*t22*t23;
    t415 = dg23*t49*t23;
    t419 = -2.0*t33*t42*t317-2.0*t33*t42*t389+4.0*t33*t42*t399+2.0*t33*t152*t384-4.0*t33*t152*t394             &
           +2.0*t33*t152*t404-32.0*t33*t162*t404+24.0*t33*t162*t415-2.0*t174*t376-2.0*t222*t376+4.0*t238*t376;
    t437 = g01*t8;
    t438 = t437*g02;
    t439 = dg23*dg03;
    t440 = t439*t237;
    t443 = t129*g33;
    t444 = t34*dg23;
    t445 = t444*t24;
    t448 = t1*g33;
    t449 = t448*t3;
    t450 = dg22*dgy02;
    t451 = t450*t64;
    t453 = t247*g33;
    t454 = g03*dgx33;
    t455 = t454*t89;
    t458 = t46*t136;
    t459 = dgy22*g02;
    t460 = t459*t89;
    t463 = t46*g33;
    t464 = t3*Theta23;
    t469 = t3*dgx03;
    t473 = t3*dgy02;
    t477 = g23*dgx01;
    t478 = dgx33*t21;
    t482 = g23*dgy01;
    t483 = dgy22*t21;
    t487 = t46*t3;
    t492 = dgx33*r;
    t496 = dgy22*r;
    t500 = Theta33*dg23;
    t501 = t500*t96;
    t503 = dgx33*dg03;
    t504 = t503*t96;
    t507 = dg02*dgy22;
    t508 = t507*t96;
    t511 = -12.0*t463*t469*t89-12.0*t463*t473*t89-2.0*t463*t477*t478-2.0*t463*t482*t483    &
           +6.0*t487*t455+6.0*t449*t460+4.0*t463*t477*t492+4.0*t463*t482*t496+t453*t501-2.0*t453*t504-2.0*t458*t508;
    t514 = Theta23*dg22;
    t515 = t514*t96;
    t517 = dg22*dgx03;
    t518 = t517*t96;
    t520 = t450*t96;
    t522 = t3*ddgxr03;
    t533 = t454*t64;
    t537 = t129*t136;
    t538 = t20*t237;
    t541 = t28*t237;
    t544 = t437*g03;
    t545 = t120*t237;
    t548 = t124*t237;
    t551 = g23*dg00;
    t555 = dg33*Theta23;
    t556 = t555*t50;
    t560 = t500*t50;
    t563 = t503*t50;
    t566 = t507*t50;
    t569 = t514*t50;
    t572 = -4.0*t537*t538+4.0*t537*t541+4.0*t544*t545-4.0*t544*t548-4.0*t537*t551*t24-2.0*t453*t556 &
           -t449*t520-2.0*t453*t560+4.0*t453*t563+4.0*t458*t566-2.0*t458*t569;
    t574 = t517*t50;
    t577 = t450*t50;
    t597 = t500*t64;
    t599 = 8.0*t463*t522*t50-4.0*t449*t566+2.0*t449*t569-2.0*t449*t574+2.0*t449*t577+4.0*t453*t533+t453*t597 &
           +2.0*t458*t574-2.0*t458*t577+2.0*t487*t560-4.0*t487*t563;
    t600 = t503*t64;
    t603 = t459*t64;
    t606 = t507*t64;
    t609 = t514*t64;
    t611 = t517*t64;
    t614 = t24*g01;
    t616 = g00*dg22;
    t617 = t8*dg33*t616;
    t620 = dg23*dg23;
    t621 = g23*t620;
    t622 = t621*t201;
    t625 = t621*t186;
    t634 = -2.0*t453*t600+4.0*t458*t603-2.0*t458*t606+t458*t609-t458*t611+t458*t451        &
           -2.0*t614*t617-2.0*t614*t622-2.0*t614*t625-4.0*t463*t522*t64-12.0*t463*t464*t64;
    t645 = t17*g33;
    t646 = t3*dgx01;
    t655 = t144*t24;
    t662 = t17*t136;
    t663 = t185*dg23;
    t664 = t663*t24;
    t667 = t8*dg00;
    t671 = t17*t8;
    t672 = t131*t24;
    t676 = t39*t8;
    t677 = t139*t24;
    t680 = t39*t3;
    t683 = g03*dg23;
    t696 = t663*t132;
    t706 = -8.0*t676*t677-8.0*t680*t664+24.0*t438*t683*t24-2.0*t544*t121+2.0*t544*t125+4.0*t537*t551*t132 &
           +8.0*t537*t145-8.0*t662*t696-8.0*t645*t667*t132+8.0*t671*t133+8.0*t676*t140;
    t738 = t237*g01;
    t739 = t128*dg33;
    t740 = t739*t444;
    t743 = t136*dg23;
    t744 = dg22*t185;
    t745 = t743*t744;
    t748 = t3*ddg23;
    t749 = t748*t201;
    t752 = t748*t186;
    t756 = t8*ddg23*t162;
    t766 = t3*t620*t162;
    t769 = t237*t128;
    t770 = t136*dg01;
    t771 = t770*t144;
    t774 = -t449*t609+4.0*t738*t617+4.0*t738*t622+4.0*t738*t625-2.0*t738*t740-2.0*t738*t745 &
           -4.0*t738*t749-4.0*t738*t752+8.0*t738*t756-4.0*t738*t766+4.0*t769*t771;
    t778 = g33*dg01*t444;
    t782 = t770*t663;
    t785 = t237*t3;
    t786 = t214*t201;
    t789 = t214*t186;
    t793 = t214*t162;
    t796 = t128*t136;
    t798 = t796*ddg23*g00;
    t801 = t128*g33;
    t803 = t801*ddg23*t34;
    t806 = g22*t136;
    t808 = t806*ddg23*t185;
    t811 = t24*t128;
    t812 = dg01*t34;
    t813 = t32*t812;
    t816 = t24*g22;
    t817 = dg01*t185;
    t818 = t137*t817;
    t823 = g00*dg01;
    t824 = g33*t8*t823;
    t831 = dg02*dg33;
    t832 = t831*t170;
    t835 = t439*t170;
    t838 = t831*t237;
    t842 = t3*t1*Rmin;
    t843 = r*g33;
    t844 = dgx23*g02;
    t849 = dgy23*dg22;
    t857 = t17*t3;
    t860 = -8.0*t816*t824+t614*t740+t614*t745+2.0*t614*t749+2.0*t438*t832-2.0*t438*t835                 &
           -4.0*t438*t838+4.0*t842*t843*t844+2.0*t842*t22*g03*t849-4.0*t842*t49*g03*t849-8.0*t857*t445;
    t862 = t831*t24;
    t865 = t439*t24;
    t868 = t444*t132;
    t873 = dg33*dgx03;
    t874 = t873*t50;
    t877 = dg33*dgy02;
    t878 = t877*t50;
    t881 = dg23*Theta22;
    t882 = t881*t50;
    t885 = t3*ddgyr02;
    t895 = 2.0*t438*t862-2.0*t438*t865-8.0*t443*t868+4.0*t857*t868-2.0*t453*t874+2.0*t453*t878-2.0*t458*t882 &
           +8.0*t463*t885*t50+2.0*t487*t556+2.0*t487*t874-2.0*t487*t878;
    t898 = t555*t64;
    t900 = t873*t64;
    t902 = t877*t64;
    t904 = t881*t64;
    t913 = t555*t96;
    t915 = -4.0*t463*t885*t64+2.0*t449*t882-t449*t904+t453*t898+t453*t900-t453*t902 &
           +t453*t913+t458*t904-t487*t898-t487*t900+t487*t902;
    t919 = t873*t96;
    t921 = t877*t96;
    t923 = t881*t96;
    t932 = t20*t170;
    t935 = t28*t170;
    t941 = t124*t170;
    t944 = t50*t1;
    t946 = t806*ddgyr22*g02;
    t949 = t64*t1;
    t953 = t806*ddgxr23*g02;
    t958 = g02*dgx22;
    t959 = t743*t958;
    t963 = t136*dg22*t459;
    t968 = t3*dg23;
    t969 = g03*dgx23;
    t970 = t968*t969;
    t973 = -2.0*t544*t171+2.0*t544*t941+4.0*t944*t946-4.0*t944*t953+2.0*t944*t959-2.0*t944*t963+4.0*t944*t970 &
           -2.0*t949*t946+2.0*t949*t953-t949*t959+t949*t963;
    t978 = g33*g02;
    t979 = t3*ddgyr22*t978;
    t985 = t3*ddgxr23*t978;
    t990 = t137*t744;
    t995 = t128*t34*dg33*g23;
    t998 = t96*t1;
    t1006 = t132*t128;
    t1007 = t137*t823;
    t1012 = t132*g22;
    t1023 = t132*g01;
    t1030 = t998*t963-4.0*t1006*t1007+4.0*t1006*t813+4.0*t1012*t818+8.0*t1012*t824+2.0*t614*t752-4.0*t614*t756 &
            -2.0*t998*t985+2.0*t1023*t990+2.0*t1023*t995+2.0*t614*t766;
    t1039 = t24*t3;
    t1044 = t24*t8;
    t1067 = t170*g01;
    t1084 = 4.0*t738*t968*t186-8.0*t738*t8*dg23*t162+2.0*t1067*t798-2.0*t1067*t803-2.0*t1067*t808+t1067*t740+t1067*t745 &
            +2.0*t1067*t749+2.0*t1067*t752-4.0*t1067*t756-2.0*t1067*t617;
    t1092 = t170*t128;
    t1100 = t170*t3;
    t1112 = -2.0*t1067*t622-2.0*t1067*t625+2.0*t1067*t766-2.0*t1092*t771+2.0*t1092*t778+2.0*t170*g22*t782 &
            -2.0*t1100*t786-2.0*t1100*t789+4.0*t170*t8*t793+4.0*t811*t1007-4.0*t738*t798;
    t1117 = t21*g03;
    t1121 = r*g22;
    t1122 = g03*dgy23;
    t1126 = t21*g22;
    t1127 = g02*dgy33;
    t1131 = t437*t49;
    t1133 = g00*dg33;
    t1134 = t23*g22*t1133;
    t1138 = t23*g33*t616;
    t1141 = t1*g23;
    t1142 = t1141*Rmin;
    t1143 = t21*t128;
    t1144 = g03*dgy33;
    t1148 = r*t128;
    t1160 = 4.0*t738*t803+4.0*t738*t808+2.0*t842*t1117*t849+4.0*t842*t1121*t1122        &
            -2.0*t842*t1126*t1127-2.0*t1131*t1134-2.0*t1131*t1138+2.0*t1142*t1143*t1144 &
            -2.0*t1142*t1148*t1144+2.0*t1142*t21*t136*t958-2.0*t1142*r*t136*t958;
    t1166 = t136*dgx02;
    t1173 = t136*Theta22;
    t1177 = g33*dgy03;
    t1181 = g33*Theta33;
    t1191 = t437*t22;
    t1203 = t21*g33;
    t1204 = g03*dgx22;
    t1212 = dg33*dgx23;
    t1231 = t128*dgy33*t683;
    t1234 = g22*g03;
    t1235 = t3*ddgyr23*t1234;
    t1239 = t3*ddgxr33*t1234;
    t1242 = 2.0*t842*t1121*t1127-2.0*t842*t1203*t1204+2.0*t842*t843*t1204+2.0*t842*t22*g02*t1212 &
            -4.0*t842*t49*g02*t1212-4.0*t842*t1203*t844+2.0*t842*t21*g02*t1212                   &
            -4.0*t842*t1126*t1122-t998*t1231-2.0*t998*t1235+2.0*t998*t1239;
    t1245 = t801*ddgyr23*g03;
    t1249 = t801*ddgxr33*g03;
    t1252 = g02*dgy23;
    t1253 = t968*t1252;
    t1256 = t739*t454;
    t1289 = t96*t46;
    t1291 = dg33*dgx33*t322;
    t1293 = 2.0*t944*t1231-t949*t1231+4.0*t944*t1235-2.0*t949*t1235-4.0*t944*t1239             &
            +2.0*t949*t1239-2.0*t944*t1256+t949*t1256-t1289*t1291-2.0*t998*t970+2.0*t998*t979;
    t1313 = t1141*t64;
    t1314 = t212*t1252;
    t1317 = t1141*t89;
    t1320 = t1141*t96;
    t1322 = t212*dg33*dgx02;
    t1326 = t1141*t50;
    t1331 = t336*t132;
    t1332 = t801*t1133;
    t1335 = t806*t616;
    t1342 = t212*t969;
    t1348 = t212*dg22*dgy03;
    t1355 = 2.0*t1313*t1322+4.0*t1313*t1342+2.0*t1313*t1348-4.0*t1317*t1342+2.0*t1320*t1348                &
            -4.0*t1326*t1322-4.0*t1326*t1348+2.0*t1331*t1332+2.0*t1331*t1335-2.0*t337*t1332-2.0*t337*t1335;
    t1357 = t50*t46;
    t1359 = g33*ddgxr33*t322;
    t1362 = t64*t46;
    t1366 = dg23*dgy23*t1234;
    t1374 = dgy33*dg23*t322;
    t1379 = t50*t448;
    t1380 = dgx33*dg22;
    t1381 = t1380*t1234;
    t1384 = dg23*dgy22;
    t1385 = t1384*t1234;
    t1388 = t849*t1234;
    t1391 = 2.0*t1357*t1291-t1362*t1291+2.0*t1313*t1366-4.0*t1326*t1366-4.0*t1357*t1359-2.0*t1357*t1374 &
            +2.0*t1362*t1359+t1362*t1374-2.0*t1379*t1381-2.0*t1379*t1385+4.0*t1379*t1388;
    t1392 = t64*t448;
    t1398 = g33*ddgyr23*t322;
    t1405 = dgx33*dg23;
    t1406 = t1405*t978;
    t1409 = g23*g03;
    t1410 = t1405*t1409;
    t1413 = t1212*t978;
    t1416 = dg33*dgy22;
    t1417 = t1416*t978;
    t1419 = t1416*t1409;
    t1421 = t1392*t1381+t1392*t1385-2.0*t1392*t1388-2.0*t1289*t1398+2.0*t1289*t1359                    &
            +2.0*t1320*t1366-2.0*t1357*t1406+2.0*t1357*t1410-2.0*t1362*t1413+t1362*t1417-t1362*t1419;
    t1428 = g33*ddgyr22*t1409;
    t1434 = g33*ddgxr23*t1409;
    t1440 = dg23*dgx22*t1409;
    t1444 = dg22*dgy22*t1409;
    t1450 = t1380*t978;
    t1454 = dg23*dgx23*t978;
    t1457 = t1384*t978;
    t1473 = g03*t620*t978;
    t1476 = 2.0*t237*t129*t183-t1313*t1450+2.0*t1313*t1454-t1313*t1457+2.0*t1326*t1450-4.0*t1326*t1454 &
            +2.0*t1326*t1457+2.0*t222*t1473+8.0*t238*t177+t222*t187-2.0*t238*t187;
    t1488 = t96*t448;
    t1490 = dg33*Theta22;
    t1501 = t185*dg33;
    t1505 = t17*t40;
    t1510 = t129*g33*g03;
    t1516 = t17*t3*g03;
    t1529 = -2.0*t33*t1501*t132-2.0*t47*t1490*t96-4.0*t1505*t373*t132+4.0*t1510*t440        &
            +2.0*t1510*t832-2.0*t1510*t835-4.0*t1510*t838-4.0*t1516*t440-2.0*t1516*t832+2.0*t1516*t835+4.0*t1516*t838;
    t1577 = 8.0*t19*t538-8.0*t19*t541-4.0*t41*t545+4.0*t41*t548-2.0*t41*t941-2.0*t168*t121+2.0*t168*t125 &
            +8.0*t130*t672+8.0*t138*t677+4.0*t168*t545-4.0*t168*t548;
    t1598 = 28.0*t19*t655-8.0*t19*t148*t24-8.0*t19*t152*t24+t1289*t1374+t1488*t1381+t1488*t1385 &
            -2.0*t1488*t1388+4.0*t1357*t1398-2.0*t1362*t1398-t1488*t1444-t1320*t1450;
    t1612 = t437*t185;
    t1615 = dg33*t49;
    t1619 = t251*dg23;
    t1624 = t251*dgy22;
    t1625 = dg03*t49;
    t1629 = 2.0*t1320*t1454-t1320*t1457-4.0*t238*t1473+2.0*t174*t1473+4.0*t1357*t1413-2.0*t1357*t1417         &
            +2.0*t1357*t1419-4.0*t1612*t311+4.0*t1612*t1615*t23-4.0*t1619*dgy02*t49*Rmin+4.0*t1624*t1625*Rmin;
    t1634 = t247*t136;
    t1636 = ddgxr03*t21*Rmin;
    t1640 = Theta23*t21*Rmin;
    t1644 = dgx03*t21*Rmin;
    t1648 = dgy02*t21*Rmin;
    t1651 = t251*g02;
    t1652 = t478*Rmin;
    t1655 = t251*g03;
    t1660 = t251*dg02;
    t1663 = t251*Theta33;
    t1664 = dg22*t21;
    t1678 = Theta23*r*Rmin;
    t1682 = dgx03*r*Rmin;
    t1686 = dgy02*r*Rmin;
    t1689 = dgy01*dgx33;
    t1690 = t1689*t21;
    t1693 = t258*dg00;
    t1696 = t258*dg02;
    t1703 = t258*g00;
    t1706 = 2.0*t1619*t1648-2.0*t1624*dg03*t21*Rmin-4.0*t1634*t1678+4.0*t1634*t1682+4.0*t1634*t1686      &
            +2.0*t453*t1690+2.0*t1693*t384-2.0*t1696*t389-4.0*t1693*t394+4.0*t1696*t399-16.0*t1703*t404;
    t1708 = t258*g02;
    t1711 = t258*g03;
    t1720 = dg01*g02*g03;
    t1724 = ddgyr02*t21*Rmin;
    t1727 = t251*dg33;
    t1733 = ddgyr02*t22*Rmin;
    t1741 = ddgyr02*t49*Rmin;
    t1752 = dg22*t22;
    t1768 = ddgxr03*t49*Rmin;
    t1775 = dg22*t49;
    t1783 = t251*Rmin;
    t1788 = t8*ddgyr23*g02;
    t1791 = 4.0*t282*t287*r+2.0*t1663*t1752*Rmin+2.0*t1619*dgx03*t22*Rmin+2.0*t1619*dgy02*t22*Rmin &
            -2.0*t1624*t316*Rmin-4.0*t1634*t1768+4.0*t1660*dgx33*t49*Rmin                          &
            -4.0*t1663*t1775*Rmin-4.0*t1619*dgx03*t49*Rmin-2.0*t1783*t1203*Theta22+2.0*t998*t1788;
    t1795 = t8*ddgxr33*g02;
    t1823 = t296*g33;
    t1824 = dgy23*dgx22;
    t1828 = dgx33*dgx22;
    t1836 = t8*ddgyr22*g03;
    t1842 = t8*ddgxr23*g03;
    t1852 = t4*ddg23*g00;
    t1856 = t214*g00;
    t1859 = 4.0*t1823*t1824*r-2.0*t1823*t1828*r-2.0*t1823*t1824*t21+4.0*t237*t4*t1856+4.0*t944*t1836    &
            -2.0*t949*t1836-2.0*t998*t1836-4.0*t944*t1842+2.0*t949*t1842+2.0*t998*t1842-4.0*t738*t1852;
    t1863 = t24*t4;
    t1875 = t812*g22;
    t1878 = t817*g33;
    t1883 = t132*t8;
    t1890 = 2.0*t614*t1852-2.0*t1863*t1856+2.0*t1067*t1852-2.0*t170*t4*t1856 &
            +4.0*t738*t4*dg23*g00+4.0*t1044*t1875+4.0*t1044*t1878            &
            -8.0*t1863*t1720-4.0*t1883*t1875-4.0*t1883*t1878-2.0*t1696*t317;
    t1891 = t437*t34;
    t1918 = t259*t21;
    t1921 = t259*r;
    t1927 = -4.0*t1891*t1752*t23+12.0*t1703*t415-8.0*t1708*t1625*t23-8.0*t1711*dg02*t49*t23 &
            +4.0*t1891*t1775*t23+12.0*t671*t34*t21*t23+12.0*t676*t185*t21*t23               &
            -24.0*t1708*t1117*t23+2.0*t537*t1918-4.0*t537*t1921-4.0*t645*t646*dgy01;
    t1931 = dgx01*dgy22;
    t1932 = t1931*t21;
    t1935 = t3*ddgxy01;
    t1949 = t1689*r;
    t1952 = t1931*r;
    t1970 = ddgxr03*t22*Rmin;
    t1995 = 4.0*t449*t1952-2.0*t463*t477*dgx33-2.0*t463*t482*dgy22+2.0*t1634*t1970-2.0*t1660*dgx33*t22*Rmin &
            +2.0*t1783*t843*Theta22-2.0*t1783*t1752*dgy03                                                   &
            +4.0*t1783*t1126*dgy03+4.0*t1783*t1775*dgy03-2.0*t1783*t1664*dgy03-2.0*t1783*t1126*Theta33;
    t2003 = t296*g22;
    t2004 = dgy33*dgx23;
    t2008 = dgy33*dgy22;
    t2017 = g33*ddgxy23;
    t2024 = g33*ddgyy22;
    t2032 = g33*ddgxx33;
    t2065 = -2.0*t2003*t2032*t21+4.0*t2003*t2032*r+t1823*t1828*t21+2.0*t282*t283*r &
            -2.0*t282*t285*r-2.0*t282*t287*t21+4.0*t282*t290*t21-2.0*t282*t293*t21 &
            +4.0*t296*t212*ddgxy23-2.0*t296*t212*ddgyy22-2.0*t296*t212*ddgxx33;
    t2085 = g22*dgy33;
    t2091 = t4*g23;
    t2093 = t2091*dg01*g00;
    t2122 = t136*dgx01;
    t2133 = 8.0*t5*t1682+8.0*t5*t1686+2.0*t251*t252*t21+2.0*t251*t255*t21+8.0*t247*t248*r-4.0*t251*t252*r-4.0*t251*t255*r &
            +2.0*t247*g33*dgy01*dgx33+2.0*t46*t2122*dgy22+8.0*t46*t18*ddgxy01-2.0*t46*t3*dgy01*dgx33;
    t2138 = g01*t2091;
    t2166 = -2.0*t448*t646*dgy22                                                   &
            -4.0*t2138*dg00*t22*t23+4.0*t2138*dg00*t49*t23+12.0*t2138*g00*t21*t23  &
            +2.0*t258*t1918-4.0*t258*t1921+2.0*t129*t2122*dgy01+2.0*t5*t1733-4.0*t5*t1741+2.0*t5*t1724+2.0*t5*t1970;
    t2192 = -4.0*t5*t1768+2.0*t5*t1636+8.0*t5*t1640-8.0*t5*t1644-8.0*t5*t1648-4.0*t247*t248*t21 &
            -8.0*t5*t1678-8.0*t282*t290*r+4.0*t282*t293*r-t282*t283*t21+t282*t285*t21;
    t2201 = t3*g22;
    Theta23_rhs = 1/t1&
                  *(-12.0*t443*g23*t34*t163-12.0*t662*g23*t185*t163           &
                    -24.0*t645*t8*g00*t163+8.0*t1711*dg02*t22*t23             &
                    +2.0*t1727*Theta22*t21*Rmin+2.0*t1727*Theta22*t22*Rmin    &
                    -4.0*t1727*Theta22*t49*Rmin-2.0*t1783*t21*dg33*dgx02      &
                    +t296*g33*dgx33*dgx22-2.0*t296*g33*dgy23*dgx22            &
                    -4.0*t645*t646*dgy01*t21+8.0*t645*t646*dgy01*r            &
                    +12.0*t537*g23*g00*t163-2.0*t33*t139*dg33*t169*t23        &
                    +4.0*t33*t139*dg33*t236*t23                               &
                    +t100+t344+t268+t191+t2133+t2065+t1927                    &
                    +t706+t599+t634+t511+t572+t419+t1706+t1476                &
                    +t1629+t1421+t1598+t1391+t1355+t1529-4.0*t453*t455        &
                    +t1859+t1995+t1791+t2166+t2192+t1293+t1242                &
                    +t1890+t1112+t1160+t1084                                  &
                    -2.0*t1783*t310*dgx02+4.0*t1783*t1203*dgx02               &
                    +4.0*t1783*t1615*dgx02-4.0*t1783*t843*dgx02               &
                    +8.0*t463*t1935*t21-8.0*t1651*t492*Rmin                   &
                    -8.0*t1655*t496*Rmin-16.0*t463*t1935*r                    &
                    -4.0*t1783*t1121*dgy03+2.0*t1783*t1121*Theta33            &
                    -2.0*t2003*t2004*t21+t2003*t2008*t21+4.0*t2003*t2004*r    &
                    -2.0*t2003*t2008*r+4.0*t2003*t2017*t21                    &
                    -8.0*t2003*t2017*r-2.0*t2003*t2024*t21+4.0*t2003*t2024*r  &
                    +t296*t298*t21-2.0*t296*t298*r+t296*t302*t21              &
                    -2.0*t296*t302*r-2.0*t296*t2085*dgx23+t296*t2085*dgy22    &
                    +4.0*t738*t968*t201-4.0*t1142*t1126*t1166                 &
                    +4.0*t1142*t1121*t1166-2.0*t1142*t1121*t1173              &
                    -4.0*t1142*t1143*t1177+2.0*t1142*t1143*t1181              &
                    +4.0*t1142*t1148*t1177-2.0*t1142*t1148*t1181              &
                    +2.0*t1142*t1126*t1173-8.0*t132*t212*t219                 &
                    +2.0*t237*t194*t198-4.0*t237*t206*t209+4.0*t47*t1490*t50  &
                    -2.0*t47*t1490*t64+2.0*t33*t1501*t24+4.0*t1505*t373*t24   &
                    +8.0*t1655*t483*Rmin+2.0*t1663*t1664*Rmin                 &
                    +8.0*t132*t4*t1720+2.0*t33*t315*t389-4.0*t33*t315*t399    &
                    +12.0*t463*t464*t89-4.0*t463*t522*t96+12.0*t463*t469*t64  &
                    +12.0*t463*t473*t64+8.0*t645*t667*t24-16.0*t438*t683*t132 &
                    -4.0*t237*g22*t782-8.0*t237*t8*t793-4.0*t463*t885*t96     &
                    +4.0*t738*t796*t144-4.0*t738*t801*t444-4.0*t738*t806*t663 &
                    -t24*t194*t198+2.0*t24*t206*t209-t170*t194*t198           &
                    +2.0*t10*ddgyy22*t21-4.0*t10*ddgyy22*r+2.0*t10*ddgxx33*t21  &
                    -4.0*t10*ddgxx33*r+8.0*t238*t239*t323+2.0*t170*t336*t339    &
                    -4.0*t237*t336*t339+2.0*t33*t35*t24                         &
                    +4.0*t41*t42*t24+4.0*t47*t48*t50+4.0*t47*t54*t50            &
                    -2.0*t47*t76*t96-2.0*t47*t48*t96-2.0*t47*t54*t96            &
                    +2.0*t47*t59*t96-4.0*t47*t72*t50+4.0*t47*t76*t50+t1030+t973 &
                    +t915+t895+t774+t860+t1577+4.0*t785*t786+2.0*t1619*t1644    &
                    -4.0*t1634*t1648+8.0*t1651*t1652+t449*t611-2.0*t1660*t1652  &
                    +2.0*t1634*t1636+4.0*t1634*t1640-4.0*t1634*t1644            &
                    +2.0*t1516*t865+2.0*t168*t941-4.0*t19*t932+4.0*t19*t935     &
                    +2.0*t41*t171+2.0*t1510*t862-2.0*t1510*t865                 &
                    -2.0*t1516*t862+4.0*t222*t327-2.0*t222*t330-2.0*t226*t333   &
                    -8.0*t238*t327+4.0*t238*t330+4.0*t366*t333-4.0*t233*t324    &
                    +t282*t285-2.0*t282*t287+4.0*t282*t290-2.0*t282*t293        &
                    +t296*t298-t282*t283+t222*t189+t226*t202+t226*t204          &
                    +4.0*t213*t215+8.0*t213*t219+t200*t204+t200*t202            &
                    +4.0*t24*t2093-4.0*t132*t2093-4.0*t453*t1949                &
                    -4.0*t458*t1952+4.0*t487*t1949+2.0*t458*t1932               &
                    -2.0*t487*t1690-2.0*t449*t1932-t487*t597                    &
                    +2.0*t487*t600-6.0*t449*t603+2.0*t449*t606+4.0*t680*t696    &
                    -8.0*t671*t672-12.0*t537*t655+2.0*t537*t25-2.0*t537*t29     &
                    +12.0*t662*t664                                             &
                    +2.0*t449*t508-t449*t515+t449*t518-6.0*t487*t533            &
                    +t458*t515-t458*t518+t458*t520-t487*t501                    &
                    +2.0*t487*t504-4.0*t458*t460+4.0*t438*t440                  &
                    +12.0*t443*t445-t449*t451                                   &
                    -8.0*t138*t140-t1289*t1419+t1289*t1406-t1289*t1410          &
                    +2.0*t1289*t1428-2.0*t1289*t1434+t1488*t1440                &
                    -2.0*t1289*t1413+t1289*t1417+4.0*t1357*t1434                &
                    -2.0*t1362*t1434-2.0*t1379*t1440+2.0*t1379*t1444            &
                    +t1392*t1440-t1392*t1444+t1362*t1406-t1362*t1410            &
                    -4.0*t1357*t1428+2.0*t1362*t1428+2.0*t1320*t1322            &
                    -8.0*t321*t215+4.0*t1313*t1314                              &
                    -4.0*t1317*t1314+2.0*t41*t121-2.0*t41*t125-8.0*t130*t133    &
                    -2.0*t998*t1795-4.0*t944*t1788                              &
                    +2.0*t949*t1788+4.0*t944*t1795                              &
                    -2.0*t949*t1795-2.0*t238*t189-2.0*t366*t202                 &
                    -2.0*t366*t204-4.0*t944*t1245+2.0*t949*t1245+4.0*t944*t1249 &
                    -2.0*t949*t1249+4.0*t944*t1253-2.0*t949*t1253               &
                    +2.0*t998*t1245-2.0*t998*t1249-2.0*t998*t1253               &
                    +t998*t1256+2.0*t1191*t1134+2.0*t1191*t1138                 &
                    -4.0*t5*ddgxy01-4.0*t10*ddgxy23                             &
                    +2.0*t10*ddgyy22+2.0*t10*ddgxx33                            &
                    -4.0*t19*t25+4.0*t19*t29+2.0*t816*t782-2.0*t1039*t786       &
                     -2.0*t1039*t789+4.0*t1044*t793+2.0*t1634*t1733             &
                    -4.0*t1634*t1741+8.0*t1708*t317+2.0*t1693*t404              &
                    +2.0*t1634*t1724-2.0*t811*t771+2.0*t811*t778-t998*t959      &
                    -4.0*t944*t979+2.0*t949*t979+4.0*t944*t985-2.0*t949*t985    &
                    -2.0*t614*t990-2.0*t614*t995-2.0*t998*t946+2.0*t998*t953    &
                    -2.0*t949*t970-t449*t923+2.0*t537*t932                      &
                    -2.0*t537*t935+t453*t919-t453*t921+t458*t923-t487*t913      &
                    -t487*t919+t487*t921-4.0*t811*t813-4.0*t816*t818            &
                    +4.0*t785*t789+2.0*t614*t798-2.0*t614*t803-2.0*t614*t808    &
                    -4.0*t769*t778) &
                   /(-2.0*t2201*g33-2.0*t4*r-2.0*t796*r+t4*t21+t796*t21-2.0*t2201*t1203+4.0*t2201*t843+t4+t796)/Rmin/t21/4.0;

    t1 = g01*g01;
    t3 = g23*g23;
    t4 = t3*t3;
    t5 = g01*t4;
    t6 = dgy01*dgy01;
    t9 = t1*t4;
    t12 = g33*g33;
    t13 = t1*g01;
    t14 = t12*t13;
    t15 = dgy22*dgy22;
    t17 = g22*ddgxy23;
    t18 = r*r;
    t23 = ddgyr03*t18*Rmin;
    t27 = Theta33*t18*Rmin;
    t31 = dgy03*t18*Rmin;
    t34 = g22*g22;
    t35 = t1*t34;
    t36 = t12*ddgyy01;
    t41 = Theta33*r*Rmin;
    t45 = dgy03*r*Rmin;
    t48 = t3*g23;
    t49 = t1*t48;
    t50 = dgx01*dgy33;
    t54 = dgy01*dgx33;
    t58 = dgy01*dgy23;
    t78 = t1*g22;
    t79 = t12*dgx01;
    t83 = 2.0*t5*t6-4.0*t9*ddgyy01+t14*t15+4.0*t14*t17*t18+4.0*t9*t23+4.0*t9*t27  &
          -8.0*t9*t31-4.0*t35*t36*t18-4.0*t9*t41+8.0*t9*t45+2.0*t49*t50*t18       &
          -2.0*t49*t54*t18+4.0*t49*t58*t18+8.0*t35*t36*r-4.0*t49*t50*r            &
          +4.0*t49*t54*r-8.0*t49*t58*r+2.0*t35*g33*dgy01*dgy33-2.0*t78*t79*dgx33;
    t87 = t3*g33;
    t95 = t1*g33;
    t96 = t3*dgx01;
    t103 = dg03*dg03;
    t104 = t18*t18;
    t105 = t104*t18;
    t107 = Rmin*Rmin;
    t108 = t103*t105*t107;
    t111 = t104*r;
    t113 = t103*t111*t107;
    t117 = t103*t104*t107;
    t120 = g01*t34;
    t121 = t12*t6;
    t128 = g01*g22;
    t133 = ddgyr03*t104*Rmin;
    t136 = t18*r;
    t138 = ddgyr03*t136*Rmin;
    t144 = g22*ddgyy22;
    t151 = g22*ddgxx33;
    t158 = dgx33*dgx22;
    t161 = dgy23*dgx22;
    t165 = 4.0*t78*t79*dgy23+8.0*t78*t87*ddgyy01-2.0*t78*t3*dgy01*dgy33+2.0*t95*t96*dgx33 &
           -4.0*t95*t96*dgy23-2.0*t5*t108+4.0*t5*t113-2.0*t5*t117+2.0*t120*t121*t18       &
           -4.0*t120*t121*r-4.0*t128*t87*t6+4.0*t9*t133-8.0*t9*t138-8.0*t14*t17*r         &
           -2.0*t14*t144*t18+4.0*t14*t144*r-2.0*t14*t151*t18+4.0*t14*t151*r+t14*t158*t18-2.0*t14*t161*t18;
    t173 = g33*t13;
    t174 = t3*ddgxy23;
    t181 = t3*ddgyy22;
    t188 = t3*ddgxx33;
    t195 = g23*dgx33;
    t212 = dgx33*dgx33;
    t213 = g22*t212;
    t219 = g22*dgy33;
    t225 = t111*t107;
    t226 = t225*t128;
    t227 = g23*g02;
    t228 = dg33*dg33;
    t230 = t227*g03*t228;
    t233 = g03*g03;
    t235 = dg33*dg23;
    t236 = g23*t233*t235;
    t239 = g01*g33;
    t240 = t225*t239;
    t241 = g02*g02;
    t243 = g23*t241*t235;
    t246 = -2.0*t14*t158*r+4.0*t14*t161*r-4.0*t173*t174*t18+8.0*t173*t174*r+2.0*t173*t181*t18                        &
           -4.0*t173*t181*r+2.0*t173*t188*t18-4.0*t173*t188*r-t173*t195*dgy22                                        &
           +t173*g23*dgy33*dgx22-2.0*t173*t195*dgx23+4.0*t173*g23*dgx23*dgy23-2.0*t173*g23*dgy23*dgy22+t173*t213*t18 &
           -2.0*t173*t213*r-2.0*t173*t219*dgx23+t173*t219*dgy22+4.0*t226*t230+4.0*t226*t236+4.0*t240*t243;
    t247 = g01*t3;
    t249 = g02*g03;
    t250 = t249*t235;
    t253 = g22*g33;
    t254 = t225*t253;
    t255 = g00*t3;
    t256 = dg01*dg33;
    t257 = t255*t256;
    t260 = t104*t107;
    t261 = t260*t128;
    t262 = g00*ddg33;
    t263 = t87*t262;
    t270 = t260*t239;
    t276 = t260*t253;
    t279 = t105*t107;
    t280 = t279*t128;
    t283 = t12*g01;
    t284 = t136*t107;
    t285 = t283*t284;
    t287 = g00*dg23;
    t288 = g22*g23*t287;
    t291 = Rmin*t104;
    t292 = t95*t291;
    t293 = g22*g03;
    t294 = dgx33*dg23;
    t295 = t293*t294;
    t298 = dgy33*dg22;
    t299 = t293*t298;
    t301 = Rmin*t136;
    t302 = t95*t301;
    t307 = Rmin*t18;
    t308 = t95*t307;
    t312 = t291*t78;
    t314 = g23*g03;
    t315 = g33*ddgyr23*t314;
    t319 = g33*ddgxr33*t314;
    t322 = t1*g23;
    t323 = t291*t322;
    t324 = g33*g02;
    t325 = t294*t324;
    t328 = dg33*dg22;
    t329 = t255*t328;
    t331 = -8.0*t225*t247*t250-8.0*t254*t257-4.0*t261*t263-2.0*t261*t230-2.0*t261*t236-2.0*t270*t243 &
           +4.0*t260*t247*t250+4.0*t276*t257-4.0*t280*t263-4.0*t285*t288 &
           +2.0*t292*t295-t292*t299-4.0*t302*t295+2.0*t302*t299+2.0*t308*t295-t308*t299 &
           +4.0*t312*t315-4.0*t312*t319+2.0*t323*t325-t270*t329;
    t334 = t239*t279;
    t336 = g22*t233*t328;
    t344 = dg33*dgx22;
    t345 = t314*t344;
    t347 = t128*t87;
    t348 = g03*dg03;
    t349 = t348*t260;
    t352 = dg00*dg33;
    t353 = t352*t260;
    t356 = g33*g23;
    t357 = t128*t356;
    t358 = t233*dg23;
    t363 = t239*t3*g02;
    t364 = dg02*dg33;
    t365 = t364*t260;
    t368 = dg23*dg03;
    t369 = t368*t260;
    t372 = t78*t356;
    t373 = dg02*dgy33;
    t377 = Theta33*dg23;
    t381 = dg33*Theta23;
    t385 = dg33*dgx03;
    t389 = dg33*dgy02;
    t393 = dgx33*dg03;
    t397 = dgy23*dg03;
    t403 = -t334*t336-t334*t329+2.0*t240*t336+2.0*t240*t329-t270*t336+t292*t345         &
           -8.0*t347*t349-4.0*t347*t353+4.0*t357*t358*t260+4.0*t363*t365-4.0*t363*t369  &
           -4.0*t372*t373*t301+4.0*t372*t377*t301+4.0*t372*t381*t301-4.0*t372*t385*t301 &
           +4.0*t372*t389*t301+4.0*t372*t393*t301-8.0*t372*t397*t301-2.0*t302*t345;
    t405 = dg33*dgx23;
    t406 = t293*t405;
    t409 = dg33*dgy22;
    t410 = t293*t409;
    t412 = t227*t405;
    t428 = g22*g00*t328;
    t433 = t283*t260;
    t437 = g02*dg22;
    t438 = t437*t314;
    t442 = dgx33*dg22*t314;
    t445 = dgy23*dg22*t314;
    t448 = t301*t78;
    t451 = t307*t78;
    t456 = t308*t345-2.0*t292*t406+t292*t410+2.0*t292*t412+4.0*t302*t406-2.0*t302*t410-4.0*t302*t412 &
           -2.0*t308*t406+t308*t410+2.0*t308*t412+t283*t279*t428                                     &
           -2.0*t283*t225*t428+4.0*t433*t288+t433*t428+4.0*t433*t438+t292*t442                       &
           -2.0*t292*t445-8.0*t448*t315+4.0*t451*t315+8.0*t448*t319;
    t460 = t301*t322;
    t463 = t298*t324;
    t466 = t307*t322;
    t470 = t409*t324;
    t474 = dg23*dgy23*t324;
    t480 = dg33*dgx33;
    t481 = t480*t324;
    t484 = t480*t314;
    t487 = dgy33*dg23;
    t488 = t487*t324;
    t491 = t487*t314;
    t508 = -4.0*t451*t319-4.0*t460*t325+2.0*t460*t463+2.0*t466*t325-t466*t463               &
           -2.0*t460*t470+8.0*t460*t474+t466*t470-4.0*t466*t474+4.0*t448*t481-4.0*t448*t484 &
           -4.0*t448*t488+4.0*t448*t491-2.0*t451*t481+2.0*t451*t484+2.0*t451*t488           &
           -2.0*t451*t491-4.0*t285*t438-2.0*t302*t442+4.0*t302*t445;
    t514 = g03*dgy23;
    t539 = g02*dgy33;
    t540 = r*Rmin;
    t544 = g03*dgx33;
    t563 = t279*t253;
    t566 = t308*t442-2.0*t308*t445+8.0*t226*t263-12.0*t372*t514*t307+2.0*t372*t373*t307  &
           -2.0*t372*t377*t307-2.0*t372*t381*t307+2.0*t372*t385*t307-2.0*t372*t389*t307  &
           -2.0*t372*t393*t307+4.0*t372*t397*t307+6.0*t372*t539*t540-2.0*t372*t544*t540  &
           +12.0*t372*t514*t540+2.0*t372*t373*t291-2.0*t280*t230-2.0*t280*t236           &
           -2.0*t334*t243+4.0*t279*t247*t250+4.0*t563*t257;
    t570 = g00*dg33;
    t574 = g22*t12;
    t576 = g03*dg01;
    t577 = t227*t576;
    t613 = t12*g23;
    t614 = t128*t613;
    t615 = g02*dg03;
    t616 = t615*t284;
    t619 = g03*dg02;
    t620 = t619*t284;
    t624 = t128*t12*g02;
    t625 = g03*dg23;
    t629 = -8.0*t226*t87*t570+8.0*t260*t574*t577-8.0*t284*t574*t577-t323*t463+t323*t470           &
           -4.0*t323*t474-2.0*t312*t481+2.0*t312*t484+2.0*t312*t488-2.0*t312*t491                 &
           -2.0*t372*t377*t291-2.0*t372*t381*t291+2.0*t372*t385*t291-2.0*t372*t389*t291           &
           -2.0*t372*t393*t291+4.0*t372*t397*t291+8.0*t614*t616+8.0*t614*t620+8.0*t624*t625*t284;
    t630 = t570*t284;
    t633 = t348*t284;
    t639 = t18*t107;
    t640 = t249*t639;
    t643 = t364*t279;
    t646 = t368*t279;
    t649 = t352*t279;
    t656 = t364*t225;
    t659 = t368*t225;
    t662 = t352*t225;
    t669 = t615*t260;
    t672 = t619*t260;
    t682 = t570*t260;
    t685 = -18.0*t347*t630+8.0*t347*t633-4.0*t357*t358*t284+24.0*t614*t640                      &
           -4.0*t624*t643+4.0*t624*t646-4.0*t347*t649+4.0*t363*t643-4.0*t363*t646+8.0*t624*t656 &
           -8.0*t624*t659+8.0*t347*t662-8.0*t363*t656+8.0*t363*t659-8.0*t614*t669-8.0*t614*t672 &
           -8.0*t624*t625*t260-4.0*t624*t365+4.0*t624*t369+26.0*t347*t682;
    t732 = -6.0*t372*t539*t307+2.0*t372*t544*t307+8.0*t9*ddgyy01*r                                  &
           -4.0*t35*t36+2.0*t49*t50-2.0*t49*t54+4.0*t49*t58+2.0*t5*t6*t18-4.0*t5*t6*r+2.0*t120*t121 &
           +4.0*t14*t17-2.0*t14*t144-2.0*t14*t151+t14*t15*t18-2.0*t14*t15*r+t14*t158                &
           -2.0*t14*t161-4.0*t173*t174+2.0*t173*t181+2.0*t173*t188;
    t737 = t78*t3;
    t738 = dgy01*dgy33;
    t739 = t738*r;
    t742 = t95*t3;
    t743 = dgx01*dgx33;
    t744 = t743*r;
    t747 = dgx01*dgy23;
    t748 = t747*r;
    t751 = t78*g33;
    t752 = g23*dgx01;
    t756 = g23*dgy01;
    t764 = t12*t1*Rmin;
    t765 = t18*t3;
    t769 = t12*g33;
    t771 = t769*t1*Rmin;
    t772 = r*g22;
    t776 = r*t3;
    t781 = g00*dg01;
    t782 = t781*dg33;
    t785 = t260*g01;
    t787 = t34*t233*t228;
    t789 = t4*g00;
    t790 = t789*ddg33;
    t796 = t279*g01;
    t803 = t225*g01;
    t807 = t260*t34;
    t809 = t769*g00*dg01;
    t812 = t173*t213-4.0*t9*ddgyy01*t18+4.0*t737*t739-4.0*t742*t744+8.0*t742*t748           &
           -2.0*t751*t752*dgy33+2.0*t751*t756*dgx33-4.0*t751*t756*dgy23+4.0*t764*t765*dgx02 &
           +4.0*t771*t772*dgx02-4.0*t764*t776*dgx02+4.0*t225*t4*t782+t785*t787              &
           +2.0*t785*t790-2.0*t260*t4*t782+t796*t787+2.0*t796*t790-2.0*t279*t4*t782+4.0*t803*t789*dg33+4.0*t807*t809;
    t816 = t12*t233*dg01;
    t819 = t260*g22;
    t820 = t769*t241;
    t821 = t820*dg01;
    t825 = t3*t241;
    t826 = t825*dg01;
    t829 = t260*g33;
    t830 = t789*dg01;
    t833 = t284*t34;
    t838 = dg33*t49;
    t851 = t49*dgx33;
    t853 = dg03*t18*Rmin;
    t856 = t49*dgy23;
    t859 = t35*t12;
    t864 = t35*g33;
    t865 = t738*t18;
    t868 = t78*t12;
    t869 = t743*t18;
    t872 = t747*t18;
    t875 = t3*ddgyy01;
    t883 = -4.0*t807*t816-4.0*t819*t821+4.0*t260*t12*t826+4.0*t829*t830-4.0*t833*t809+4.0*t833*t816      &
           +2.0*t838*Theta23*t18*Rmin-2.0*t838*dgx03*t18*Rmin                                            &
           +2.0*t838*dgy02*t18*Rmin+2.0*t851*t853-4.0*t856*t853-6.0*t859*t41+12.0*t859*t45+2.0*t864*t865 &
           -2.0*t868*t869+4.0*t868*t872+8.0*t751*t875*t18-2.0*t737*t865+2.0*t742*t869;
    t886 = t49*g02;
    t887 = dgy33*r;
    t891 = t49*g03;
    t892 = dgx33*r;
    t896 = dgy23*r;
    t910 = t173*g23;
    t911 = dgx33*dgx23;
    t915 = dgx23*dgy23;
    t919 = dgy23*dgy22;
    t923 = t173*g22;
    t924 = dgy33*dgx23;
    t928 = dgy33*dgy22;
    t937 = dgx33*dgy22;
    t940 = dgy33*dgx22;
    t952 = -4.0*t742*t872-4.0*t886*t887*Rmin+4.0*t891*t892*Rmin-8.0*t891*t896*Rmin          &
           -4.0*t864*t739+4.0*t868*t744-2.0*t771*t772*Theta22+2.0*t764*t776*Theta22         &
           +4.0*t910*t911*r-8.0*t910*t915*r+4.0*t910*t919*r-2.0*t923*t924*t18+t923*t928*t18 &
           +4.0*t923*t924*r-2.0*t923*t928*r-t910*t937*t18+t910*t940*t18+2.0*t910*t937*r     &
           -2.0*t910*t940*r-2.0*t910*t911*t18;
    t960 = t95*Rmin;
    t961 = r*t48;
    t971 = t18*t48;
    t981 = t18*g22;
    t988 = t291*t1;
    t990 = t48*ddgyr23*g03;
    t994 = t48*ddgxr33*g03;
    t997 = t3*t233;
    t998 = dg23*dg23;
    t999 = t997*t998;
    t1004 = t284*g01;
    t1005 = t820*dg22;
    t1008 = t301*t1;
    t1011 = t307*t1;
    t1020 = 4.0*t910*t915*t18-2.0*t910*t919*t18-4.0*t960*t961*Theta23                               &
            +4.0*t960*t961*dgx03+4.0*t960*t961*dgy02+4.0*t960*t971*Theta23-4.0*t960*t971*dgx03      &
            -4.0*t960*t971*dgy02+2.0*t771*t981*Theta22-2.0*t764*t765*Theta22                        &
            -4.0*t988*t990+4.0*t988*t994+4.0*t803*t999-2.0*t785*t999+2.0*t1004*t1005+8.0*t1008*t990 &
            -4.0*t1011*t990-8.0*t1008*t994+4.0*t1011*t994-2.0*t796*t999;
    t1023 = t18*g02;
    t1051 = t284*g22;
    t1057 = t284*g33;
    t1064 = dg03*t104;
    t1065 = t1064*Rmin;
    t1072 = t49*dg02;
    t1077 = t49*Theta33;
    t1086 = -2.0*t785*t1005+2.0*t771*t1023*dgx22-2.0*t771*r*g02*dgx22-4.0*t771*t981*dgx02-8.0*t868*t748 &
            -16.0*t751*t875*r+2.0*t838*Theta23*t104*Rmin                                                &
            -2.0*t838*dgx03*t104*Rmin+2.0*t838*dgy02*t104*Rmin                                          &
            +4.0*t1051*t821-4.0*t284*t12*t826-4.0*t1057*t830-2.0*t803*t787                              &
            -4.0*t803*t790+2.0*t851*t1065-4.0*t856*t1065                                                &
            -8.0*t859*t138+4.0*t1072*dgy33*t136*Rmin-4.0*t1077*dg23*t136*Rmin-4.0*t838*Theta23*t136*Rmin;
    t1100 = dg03*t136*Rmin;
    t1111 = dgy33*t18;
    t1112 = t1111*Rmin;
    t1115 = dgx33*t18;
    t1119 = dgy23*t18;
    t1125 = t120*t12;
    t1128 = t5*dg00;
    t1130 = dg33*t105*t107;
    t1136 = dg33*t111*t107;
    t1139 = t120*t769;
    t1141 = dg00*t104*t107;
    t1146 = t239*t4;
    t1149 = t5*g00;
    t1151 = dg33*t104*t107;
    t1154 = 4.0*t838*dgx03*t136*Rmin-4.0*t838*dgy02*t136*Rmin-4.0*t851*t1100+8.0*t856*t1100  &
            +4.0*t859*t23+6.0*t859*t27-12.0*t859*t31+4.0*t886*t1112                          &
            -4.0*t891*t1115*Rmin+8.0*t891*t1119*Rmin-2.0*t1072*t1112-2.0*t1125*t108          &
            +2.0*t1128*t1130+4.0*t1125*t113-4.0*t1128*t1136-4.0*t1139*t1141                  &
            -2.0*t1125*t117-4.0*t1146*t1141-12.0*t1149*t1151;
    t1158 = dg00*t136*t107;
    t1164 = dg33*t136*t107;
    t1168 = g00*t18*t107;
    t1175 = t128*t769;
    t1177 = t241*t18*t107;
    t1180 = t283*t3;
    t1185 = t128*g33;
    t1186 = t3*t6;
    t1208 = t227*t576*dg33;
    t1212 = t356*t570*dg23;
    t1216 = t356*t249*ddg33;
    t1223 = 2.0*t1128*t1151+4.0*t1139*t1158+4.0*t1146*t1158+8.0*t1149*t1164                 &
            +12.0*t1139*t1168-12.0*t1125*t233*t18*t107-12.0*t1175*t1177+12.0*t1180*t1177    &
            +12.0*t1146*t1168-4.0*t1185*t1186*t18+8.0*t1185*t1186*r+2.0*t1077*dg23*t18*Rmin &
            +4.0*t859*t133-2.0*t1072*dgy33*t104*Rmin+2.0*t1077*dg23*t104*Rmin               &
            +8.0*t254*t1208+2.0*t261*t1212+4.0*t261*t1216-4.0*t276*t1208+2.0*t280*t1212;
    t1232 = g02*dg33*dg22*g23*g03;
    t1276 = t241*ddg33;
    t1277 = t87*t1276;
    t1281 = t48*g00*t235;
    t1284 = t48*g02;
    t1286 = t1284*g03*ddg33;
    t1290 = t12*g00*t256;
    t1293 = 4.0*t280*t1216-4.0*t563*t1208+4.0*t240*t1232-2.0*t270*t1232-2.0*t334*t1232            &
            +4.0*t357*t619*t1130-4.0*t357*t625*dg03*t105*t107-8.0*t357*t619*t1136                 &
            +8.0*t357*t625*dg03*t111*t107-20.0*t357*t249*t1151+4.0*t357*t619*t1151                &
            -4.0*t357*t625*t1064*t107+12.0*t357*t249*t1164+8.0*t226*t356*t249*dg33-4.0*t226*t1212 &
            -8.0*t226*t1216+2.0*t785*t1277-2.0*t785*t1281-4.0*t785*t1286-2.0*t807*t1290;
    t1295 = g33*t233*t256;
    t1298 = t12*t241;
    t1299 = t1298*t256;
    t1302 = t997*t256;
    t1305 = t825*t256;
    t1309 = t249*t256;
    t1313 = t12*t3*t781;
    t1317 = t87*t233*dg01;
    t1320 = t1284*t576;
    t1323 = t3*g22;
    t1324 = t233*dg33;
    t1325 = t1323*t1324;
    t1328 = t241*dg33;
    t1332 = g03*dgy33;
    t1333 = t1332*t540;
    t1336 = g02*dgx33;
    t1337 = t1336*t540;
    t1340 = g02*dgy23;
    t1341 = t1340*t540;
    t1344 = t3*Theta33;
    t1348 = t3*dgy03;
    t1365 = 2.0*t807*t1295+2.0*t819*t1299-2.0*t819*t1302-2.0*t829*t1305+4.0*t260*t48*t1309          &
            -8.0*t819*t1313+4.0*t819*t1317-8.0*t829*t1320+4.0*t803*t1325                            &
             +4.0*t803*t87*t1328-6.0*t864*t1333+6.0*t868*t1337-12.0*t868*t1341+10.0*t751*t1344*t540 &
            -20.0*t751*t1348*t540-2.0*t751*t752*t1111+2.0*t751*t756*t1115                           &
            -4.0*t751*t756*t1119+4.0*t737*t1333-8.0*t742*t1337;
    t1381 = dgy33*dg03;
    t1382 = t1381*t307;
    t1385 = t3*ddgyr03;
    t1389 = Theta33*dg33;
    t1390 = t1389*t291;
    t1393 = t1381*t291;
    t1396 = dg02*dgx33;
    t1397 = t1396*t291;
    t1400 = dg02*dgy23;
    t1401 = t1400*t291;
    t1404 = Theta23*dg23;
    t1405 = t1404*t291;
    t1408 = dg23*dgx03;
    t1409 = t1408*t291;
    t1412 = t1389*t301;
    t1415 = t1381*t301;
    t1418 = t1396*t301;
    t1421 = t1400*t301;
    t1424 = t1404*t301;
    t1427 = t1408*t301;
    t1430 = 8.0*t742*t1341+4.0*t751*t752*t887-4.0*t751*t756*t892+8.0*t751*t756*t896           &
            +2.0*t1125*t649+2.0*t737*t1382-8.0*t751*t1385*t291-2.0*t737*t1390+2.0*t737*t1393  &
            -2.0*t742*t1397+4.0*t742*t1401-4.0*t742*t1405+4.0*t742*t1409-4.0*t864*t1412 &
            +4.0*t864*t1415-4.0*t868*t1418+8.0*t868*t1421-8.0*t868*t1424+8.0*t868*t1427;
    t1436 = t1336*t307;
    t1439 = t1340*t307;
    t1442 = t1396*t307;
    t1445 = t1400*t307;
    t1448 = t1404*t307;
    t1453 = t1332*t307;
    t1456 = t1389*t307;
    t1471 = t1408*t307;
    t1477 = g03*dgx23;
    t1481 = g03*dgy22;
    t1485 = 16.0*t751*t1385*t301+4.0*t737*t1412+8.0*t742*t1436-8.0*t742*t1439-2.0*t742*t1442 &
            +4.0*t742*t1445-4.0*t742*t1448-8.0*t742*t1427+6.0*t864*t1453                     &
            +2.0*t864*t1456-2.0*t864*t1382-6.0*t868*t1436+12.0*t868*t1439                    &
            +2.0*t868*t1442-4.0*t868*t1445+4.0*t868*t1448-4.0*t868*t1471                     &
            -8.0*t751*t1385*t307-4.0*t764*t772*t1477+2.0*t764*t772*t1481;
    t1487 = r*g23;
    t1488 = dgx23*g02;
    t1493 = t12*dgx33*t437;
    t1496 = t12*dgy23*t437;
    t1500 = t3*ddgyr23*t324;
    t1504 = t3*ddgxr33*t324;
    t1507 = t1298*t328;
    t1509 = t997*t328;
    t1513 = t613*t241*dg23;
    t1516 = t233*dg22;
    t1517 = t87*t1516;
    t1521 = t574*ddgyr23*g02;
    t1527 = t574*ddgxr33*g02;
    t1540 = t34*g33;
    t1541 = t1540*t1324;
    t1550 = 4.0*t764*t1487*t1488-t988*t1493+2.0*t988*t1496+4.0*t988*t1500-4.0*t988*t1504+t796*t1507+2.0*t796*t1509 &
            +4.0*t785*t1513-4.0*t785*t1517+8.0*t1008*t1521                                                         &
            -4.0*t1011*t1521-8.0*t1008*t1527-10.0*t751*t1344*t307+20.0*t751*t1348*t307-4.0*t737*t1453              &
            -2.0*t737*t1456-6.0*t1004*t1541+8.0*t1004*t1325+8.0*t1051*t1313-4.0*t1051*t1317;
    t1553 = t3*t103;
    t1557 = g01*t48;
    t1558 = t1557*g03;
    t1578 = g02*dg02;
    t1579 = t1578*t260;
    t1582 = t128*t12;
    t1583 = t3*dg00;
    t1587 = t1328*t260;
    t1595 = t239*t48;
    t1600 = t239*t3;
    t1603 = t1557*g02;
    t1604 = g03*dg33;
    t1608 = 8.0*t1057*t1320+4.0*t1185*t1553*t279-4.0*t1558*t643+4.0*t1558*t646-4.0*t1125*t662     &
            -8.0*t1185*t1553*t225+8.0*t1558*t656-8.0*t1558*t659-14.0*t1125*t682                   &
            +8.0*t1125*t349+2.0*t1125*t353+8.0*t1175*t1579+8.0*t1582*t1583*t260+14.0*t1582*t1587  &
            +4.0*t1185*t1553*t260-8.0*t1180*t1579+8.0*t1595*t669+8.0*t1595*t672-16.0*t1600*t1587+24.0*t1603*t1604*t260;
    t1620 = t1578*t284;
    t1626 = t1328*t284;
    t1656 = -4.0*t1558*t365+4.0*t1558*t369+10.0*t1125*t630-8.0*t1125*t633-8.0*t1175*t1620-8.0*t1582*t1583*t284 &
            -10.0*t1582*t1626+8.0*t1180*t1620-8.0*t1595*t616                                                   &
            -8.0*t1595*t620+12.0*t1600*t1626-16.0*t1603*t1604*t284-24.0*t1582*t255*t639                        &
            +12.0*t1185*t997*t639-24.0*t1595*t640+2.0*t864*t1390-2.0*t864*t1393+2.0*t868*t1397-4.0*t868*t1401;
    t1661 = g00*t998;
    t1662 = t87*t1661;
    t1665 = t18*g23;
    t1666 = dgy22*g02;
    t1673 = g23*Theta23;
    t1677 = g23*dgx03;
    t1681 = g23*dgy02;
    t1685 = g22*t136;
    t1686 = dg33*Theta22;
    t1692 = t104*g22;
    t1693 = dg33*dgx02;
    t1717 = t34*t12;
    t1718 = t1717*t262;
    t1721 = g00*t228;
    t1722 = t1540*t1721;
    t1724 = 4.0*t868*t1405-4.0*t868*t1409+2.0*t785*t1662-2.0*t764*t1665*t1666+2.0*t764*t1487*t1666 &
            -4.0*t764*t981*t1673+4.0*t764*t981*t1677+4.0*t764*t981*t1681                           &
            +2.0*t764*t1685*t1686-t764*t981*t1686+2.0*t764*t1692*t1693-4.0*t764*t1685*t1693        &
            +2.0*t764*t981*t1693-t764*t1692*t1686+4.0*t764*t772*t1673-4.0*t764*t772*t1677          &
            -4.0*t764*t772*t1681-8.0*t803*t1284*t1604+2.0*t796*t1718-t796*t1722;
    t1726 = t233*ddg33;
    t1727 = t1540*t1726;
    t1730 = t574*t1661;
    t1733 = t574*t1276;
    t1737 = t253*t241*t228;
    t1740 = t253*t233*t998;
    t1743 = t1323*t1721;
    t1745 = t1323*t1726;
    t1756 = t279*t34;
    t1761 = t279*g22;
    t1780 = -2.0*t796*t1727-2.0*t796*t1730-2.0*t796*t1733+t796*t1737+2.0*t796*t1740+t796*t1743 &
            +2.0*t796*t1745+2.0*t796*t1662+2.0*t796*t1277-2.0*t796*t1281                       &
            -4.0*t796*t1286-2.0*t1756*t1290+2.0*t1756*t1295+2.0*t1761*t1299-2.0*t1761*t1302    &
            -2.0*t279*g33*t1305+4.0*t279*t48*t1309+4.0*t803*t1717*t570-4.0*t803*t1541-4.0*t803*t574*t1328;
    t1807 = t225*t34;
    t1812 = t225*g22;
    t1825 = -4.0*t803*t1718+2.0*t803*t1722+4.0*t803*t1727+4.0*t803*t1730+4.0*t803*t1733   &
            -2.0*t803*t1737-4.0*t803*t1740-2.0*t803*t1743-4.0*t803*t1745-4.0*t803*t1662   &
            -4.0*t803*t1277+4.0*t803*t1281+8.0*t803*t1286+4.0*t1807*t1290-4.0*t1807*t1295 &
            -4.0*t1812*t1299+4.0*t1812*t1302+4.0*t225*g33*t1305-8.0*t225*t48*t1309+2.0*t785*t1718;
    t1845 = t136*t3;
    t1846 = Theta33*dg22;
    t1852 = t769*g01;
    t1854 = t107*g22;
    t1855 = g00*dg22;
    t1856 = t1854*t1855;
    t1859 = t283*t104;
    t1860 = t1854*t1516;
    t1864 = t107*t3*t1855;
    t1869 = t107*t48*t287;
    t1875 = t104*t3;
    t1888 = g03*dgx22;
    t1894 = -4.0*t737*t1415+4.0*t742*t1418-8.0*t742*t1421+8.0*t742*t1424+4.0*t764*t981*t1477  &
            -2.0*t764*t981*t1481-4.0*t764*t1665*t1488-2.0*t960*t1845*t1846                    &
            +t960*t765*t1846-2.0*t1852*t104*t1856+2.0*t1859*t1860+2.0*t1859*t1864             &
            -4.0*t239*t104*t1869+4.0*t239*t136*t1869+t960*t1875*t1686-2.0*t960*t1845*t1686    &
            -t764*t104*g02*t344+2.0*t764*t136*g02*t344-2.0*t764*t1665*t1888-t764*t1023*t344;
    t1901 = t3*dgx33*t625;
    t1906 = t3*dgy33*dg22*g03;
    t1914 = t3*dg33*t1481;
    t1918 = t3*dg23*t514;
    t1944 = 2.0*t764*t1487*t1888+4.0*t1011*t1527+8.0*t1008*t1901-4.0*t1008*t1906             &
           -4.0*t1011*t1901+2.0*t1011*t1906+4.0*t1008*t1914-8.0*t1008*t1918-2.0*t1011*t1914  &
           +4.0*t1011*t1918+2.0*t1008*t1493-4.0*t1008*t1496-t1011*t1493+2.0*t1011*t1496 &
           -8.0*t1008*t1500+4.0*t1011*t1500+8.0*t1008*t1504-4.0*t1011*t1504-2.0*t803*t1507-4.0*t803*t1509;
    t1956 = t283*t136;
    t1988 = t785*t1507+2.0*t785*t1509-4.0*t1004*t1513+4.0*t1004*t1517+2.0*t1852*t136*t1856  &
            -2.0*t1956*t1860-2.0*t1956*t1864+4.0*t960*t765*t1481-4.0*t960*t776*t1481        &
            -4.0*t988*t1521+4.0*t988*t1527-4.0*t988*t1901+2.0*t988*t1906-2.0*t988*t1914     &
            +4.0*t988*t1918-t785*t1722-2.0*t785*t1727+10.0*t785*t1541-2.0*t785*t1730-2.0*t785*t1733;
    t1997 = dg22*dgy03;
    t2038 = t785*t1737+2.0*t785*t1740+t785*t1743+2.0*t785*t1745-12.0*t785*t1325+2.0*t764*t1692*t1997 &
            -4.0*t764*t1685*t1997+2.0*t764*t981*t1997-t764*t1692*t1846                               &
            +2.0*t764*t1685*t1846-t764*t981*t1846+t960*t765*t1686-2.0*t960*t1875*t1693               &
            +4.0*t960*t1845*t1693-2.0*t960*t765*t1693-2.0*t960*t1875*t1997+4.0*t960*t1845*t1997      &
           -2.0*t960*t765*t1997+t960*t1875*t1846+4.0*t742*t1471;
    Theta33_rhs = 1/t1 &
                  *(t2038+t83+t1988+t1944+t1894+t1825+t1430+t1365+t812+t566+t508+t1550+t685+t629 & 
                  +t456+t403+t1086+t1020+t1223+t1154+t1485+t1293+t331+t246+t165 &
                  +t1780+t1608+t952+t883+t1724+t1656+t732) /                    &
               (4.0*t1323*g33*r-2.0*t1323*g33*t18-2.0*t1323*g33-2.0*t4*r-2.0*t1717*r+t4*t18+t1717*t18+t4+t1717)/Rmin/t18/4.0;

  return

end subroutine Theta_rhs2
!---------------------------------------------------------------------------------
subroutine pg0a_rhs(Rmin,r,p02,p03,g02,g03,g22,g23,g33,dg22,dg23,dg33,ddg22,ddg23,ddg33,g01, &
                    dg01,dg02,dg03,          &
                    dgx01,dgx22,dgx23,dgx33, &
                    dgy01,dgy22,dgy23,dgy33, &
                    ddgxr01,ddgxr22,ddgxr23,ddgxr33, &
                    ddgyr01,ddgyr22,ddgyr23,ddgyr33, &
                    g02_rhs,g03_rhs,p02_rhs,p03_rhs)

  implicit none

!~~~~~~% Input parameters:
  real*8,intent(in) :: Rmin,r,p02,p03,g02,g03,g22,g23,g33,dg22,dg23,dg33,ddg22,ddg23,ddg33,g01
  real*8,intent(in) :: dg01,dg02,dg03
  real*8,intent(in) :: dgx01,dgx22,dgx23,dgx33
  real*8,intent(in) :: dgy01,dgy22,dgy23,dgy33
  real*8,intent(in) :: ddgxr01,ddgxr22,ddgxr23,ddgxr33
  real*8,intent(in) :: ddgyr01,ddgyr22,ddgyr23,ddgyr33
  real*8,intent(out) :: g02_rhs,g03_rhs,p02_rhs,p03_rhs

  real*8 :: t1;
  real*8 :: t10;
  real*8 :: t100;
  real*8 :: t101;
  real*8 :: t104;
  real*8 :: t105;
  real*8 :: t108;
  real*8 :: t11;
  real*8 :: t110;
  real*8 :: t112;
  real*8 :: t117;
  real*8 :: t118;
  real*8 :: t123;
  real*8 :: t125;
  real*8 :: t126;
  real*8 :: t129;
  real*8 :: t130;
  real*8 :: t132;
  real*8 :: t136;
  real*8 :: t14;
  real*8 :: t141;
  real*8 :: t142;
  real*8 :: t143;
  real*8 :: t150;
  real*8 :: t151;
  real*8 :: t155;
  real*8 :: t158;
  real*8 :: t16;
  real*8 :: t161;
  real*8 :: t162;
  real*8 :: t167;
  real*8 :: t168;
  real*8 :: t17;
  real*8 :: t172;
  real*8 :: t177;
  real*8 :: t179;
  real*8 :: t185;
  real*8 :: t186;
  real*8 :: t190;
  real*8 :: t191;
  real*8 :: t192;
  real*8 :: t196;
  real*8 :: t2;
  real*8 :: t20;
  real*8 :: t202;
  real*8 :: t203;
  real*8 :: t21;
  real*8 :: t213;
  real*8 :: t216;
  real*8 :: t22;
  real*8 :: t220;
  real*8 :: t221;
  real*8 :: t224;
  real*8 :: t225;
  real*8 :: t23;
  real*8 :: t230;
  real*8 :: t231;
  real*8 :: t234;
  real*8 :: t235;
  real*8 :: t243;
  real*8 :: t250;
  real*8 :: t255;
  real*8 :: t256;
  real*8 :: t26;
  real*8 :: t260;
  real*8 :: t261;
  real*8 :: t264;
  real*8 :: t265;
  real*8 :: t269;
  real*8 :: t27;
  real*8 :: t272;
  real*8 :: t275;
  real*8 :: t279;
  real*8 :: t283;
  real*8 :: t292;
  real*8 :: t30;
  real*8 :: t302;
  real*8 :: t304;
  real*8 :: t309;
  real*8 :: t31;
  real*8 :: t310;
  real*8 :: t315;
  real*8 :: t316;
  real*8 :: t317;
  real*8 :: t32;
  real*8 :: t320;
  real*8 :: t327;
  real*8 :: t328;
  real*8 :: t329;
  real*8 :: t339;
  real*8 :: t342;
  real*8 :: t345;
  real*8 :: t353;
  real*8 :: t36;
  real*8 :: t363;
  real*8 :: t374;
  real*8 :: t378;
  real*8 :: t381;
  real*8 :: t39;
  real*8 :: t392;
  real*8 :: t394;
  real*8 :: t398;
  real*8 :: t401;
  real*8 :: t403;
  real*8 :: t406;
  real*8 :: t410;
  real*8 :: t415;
  real*8 :: t43;
  real*8 :: t430;
  real*8 :: t431;
  real*8 :: t433;
  real*8 :: t434;
  real*8 :: t436;
  real*8 :: t44;
  real*8 :: t442;
  real*8 :: t445;
  real*8 :: t448;
  real*8 :: t45;
  real*8 :: t451;
  real*8 :: t453;
  real*8 :: t459;
  real*8 :: t462;
  real*8 :: t464;
  real*8 :: t466;
  real*8 :: t469;
  real*8 :: t475;
  real*8 :: t48;
  real*8 :: t483;
  real*8 :: t487;
  real*8 :: t49;
  real*8 :: t492;
  real*8 :: t496;
  real*8 :: t499;
  real*8 :: t514;
  real*8 :: t518;
  real*8 :: t530;
  real*8 :: t541;
  real*8 :: t544;
  real*8 :: t56;
  real*8 :: t568;
  real*8 :: t58;
  real*8 :: t594;
  real*8 :: t6;
  real*8 :: t67;
  real*8 :: t69;
  real*8 :: t7;
  real*8 :: t71;
  real*8 :: t73;
  real*8 :: t77;
  real*8 :: t8;
  real*8 :: t80;
  real*8 :: t81;
  real*8 :: t82;
  real*8 :: t86;
  real*8 :: t87;
  real*8 :: t89;
  real*8 :: t9;
  real*8 :: t93;
  real*8 :: t94;
  real*8 :: t97;
  real*8 :: t98;
  real*8 :: t99;

  real*8 :: t111;
  real*8 :: t115;
  real*8 :: t12;
  real*8 :: t121;
  real*8 :: t13;
  real*8 :: t133;
  real*8 :: t134;
  real*8 :: t137;
  real*8 :: t139;
  real*8 :: t140;
  real*8 :: t144;
  real*8 :: t147;
  real*8 :: t148;
  real*8 :: t15;
  real*8 :: t153;
  real*8 :: t164;
  real*8 :: t170;
  real*8 :: t18;
  real*8 :: t182;
  real*8 :: t188;
  real*8 :: t19;
  real*8 :: t193;
  real*8 :: t197;
  real*8 :: t206;
  real*8 :: t215;
  real*8 :: t222;
  real*8 :: t227;
  real*8 :: t238;
  real*8 :: t239;
  real*8 :: t24;
  real*8 :: t240;
  real*8 :: t241;
  real*8 :: t244;
  real*8 :: t245;
  real*8 :: t249;
  real*8 :: t25;
  real*8 :: t252;
  real*8 :: t257;
  real*8 :: t259;
  real*8 :: t263;
  real*8 :: t266;
  real*8 :: t270;
  real*8 :: t274;
  real*8 :: t288;
  real*8 :: t29;
  real*8 :: t293;
  real*8 :: t294;
  real*8 :: t301;
  real*8 :: t323;
  real*8 :: t326;
  real*8 :: t330;
  real*8 :: t331;
  real*8 :: t334;
  real*8 :: t335;
  real*8 :: t338;
  real*8 :: t343;
  real*8 :: t35;
  real*8 :: t350;
  real*8 :: t351;
  real*8 :: t356;
  real*8 :: t357;
  real*8 :: t361;
  real*8 :: t375;
  real*8 :: t38;
  real*8 :: t385;
  real*8 :: t388;
  real*8 :: t389;
  real*8 :: t40;
  real*8 :: t407;
  real*8 :: t41;
  real*8 :: t411;
  real*8 :: t419;
  real*8 :: t422;
  real*8 :: t428;
  real*8 :: t443;
  real*8 :: t450;
  real*8 :: t456;
  real*8 :: t46;
  real*8 :: t465;
  real*8 :: t471;
  real*8 :: t481;
  real*8 :: t486;
  real*8 :: t50;
  real*8 :: t504;
  real*8 :: t51;
  real*8 :: t534;
  real*8 :: t547;
  real*8 :: t55;
  real*8 :: t562;
  real*8 :: t592;
  real*8 :: t62;
  real*8 :: t63;
  real*8 :: t66;
  real*8 :: t70;
  real*8 :: t72;
  real*8 :: t76;
  real*8 :: t84;
  real*8 :: t92;

    t1 = r*r;
    t2 = t1*r;
    t6 = g22*g01;
    t7 = t6*Rmin;
    t8 = t1*t1;
    t9 = t8*dg02;
    t10 = g23*g23;
    t11 = dg33*t10;
    t14 = Rmin*t2;
    t16 = t10*t10;
    t17 = g01*dg02*t16;
    t20 = g33*g33;
    t21 = t20*g01;
    t22 = t21*dgx01;
    t23 = r*dg22;
    t26 = g01*t1;
    t27 = ddgxr01*t16;
    t30 = t10*g23;
    t31 = g01*g01;
    t32 = t30*t31;
    t36 = g01*r;
    t39 = dgx01*t1;
    t43 = dgx01*g01;
    t44 = g22*g22;
    t45 = t44*t20;
    t48 = dgx01*dg01;
    t49 = r*t16;
    t56 = t14*g02;
    t58 = t21*ddg22*g22;
    t67 = t14*g03;
    t69 = dg22*dg22;
    t71 = g23*g33*g01*t69;
    t73 = t30*g01;
    t77 = t73*ddg22;
    t80 = g23*t31;
    t81 = t80*g22;
    t82 = r*dg33;
    t86 = g03*g23;
    t87 = t14*t86;
    t89 = g22*g33;
    t93 = t10*t31;
    t94 = dgx22*r;
    t97 = Rmin*t8;
    t98 = t97*g03;
    t99 = t10*g01;
    t100 = dg22*dg23;
    t101 = t99*t100;
    t104 = t30*dg01;
    t105 = t104*dg22;
    t108 = -2.0*t56*t58+2.0*t32*r*ddgyr22-2.0*t32*t1*ddgyr22-t67*t71+4.0*t67*t73*dg22-2.0*t67*t77 &
           -2.0*t81*t82*dgx23-4.0*t87*g01*dg22*t89+t93*t94*dg33-2.0*t98*t101+2.0*t67*t105;
    t110 = t80*r;
    t112 = ddgyr22*g22*g33;
    t117 = dg33*dg22;
    t118 = Rmin*g02*t117;
    t123 = g02*g33;
    t125 = g23*g01;
    t126 = t125*t100;
    t129 = t97*g02;
    t130 = g33*g01;
    t132 = t130*ddg22*t10;
    t136 = g01*ddg22*t89;
    t141 = dg23*dg23;
    t142 = g01*t141;
    t143 = t89*t142;
    t150 = dg22*g22;
    t151 = t20*dg01*t150;
    t155 = dg01*dg22*t89;
    t158 = -2.0*t110*t112+2.0*t99*t2*t118+4.0*t43*t16+3.0*t97*t123*t126-2.0*t129*t132+2.0*t87*t136 &
           +2.0*t67*t101+2.0*t56*t143-3.0*t14*t123*t126+2.0*t56*t151-2.0*t87*t155;
    t161 = g22*t31;
    t162 = t161*g33;
    t167 = dg22*t10;
    t168 = g33*dg01*t167;
    t172 = g23*g22*t142;
    t177 = t21*t69;
    t179 = t97*t86;
    t185 = t1*g03;
    t186 = dg23*t10;
    t190 = t44*g01;
    t191 = t190*Rmin;
    t192 = dg23*g33;
    t196 = -2.0*t129*t143-2.0*t162*t23*dgy23+2.0*t129*t168-2.0*t67*t172+2.0*t98*t172-t129*t177 &
           +2.0*t179*t155+t98*t71+2.0*t98*t77+8.0*t7*t185*t186-8.0*t191*t185*t192;
    t202 = t130*Rmin;
    t203 = t1*g02;
    t213 = dgy33*t1;
    t216 = t21*Rmin;
    t220 = t44*t31;
    t221 = dgy33*r;
    t224 = r*dg23;
    t225 = t224*dgy22;
    t230 = t6*dgx01;
    t231 = t1*dg33;
    t234 = Rmin*t1;
    t235 = t125*t234;
    t243 = t2*dg02;
    t250 = dg02*dg01*t16;
    t255 = t125*dgx01;
    t256 = t1*dg23;
    t260 = dg01*t44;
    t261 = t260*t20;
    t264 = t73*Rmin;
    t265 = t2*dg03;
    t269 = t97*dg02;
    t272 = t230*t231*t10+8.0*t235*g03*dg22*t89-2.0*t99*t8*t118-t216*t243*t150-2.0*t73*t39*dg23 &
          +2.0*t97*t250+t202*t243*t167+2.0*t255*t256*t89-2.0*t39*t261-2.0*t264*t265*dg22+2.0*t269*t261;
    t275 = r*t44*t20;
    t279 = t89*t10;
    t283 = ddgxr01*t44*t20;
    t292 = t1*dg22;
    t302 = t14*g01;
    t304 = dg02*t44*t20;
    t309 = 2.0*t48*t275+4.0*t36*ddgxr01*t279+2.0*t26*t283-8.0*t264*t203*dg23-8.0*t264*t185*dg22 &
          -t22*t292*g22-2.0*t36*t283+4.0*t39*dg01*t279+8.0*t234*t17-4.0*t302*t304-t202*t9*t167;
    t310 = t234*g01;
    t315 = dg01*g22;
    t316 = g33*t10;
    t317 = t315*t316;
    t320 = t8*dg03;
    t327 = t125*t14;
    t328 = g03*g22;
    t329 = t328*t117;
    t339 = t14*dg02;
    t342 = t256*dgy22;
    t345 = 8.0*t310*t304+t216*t9*t150-4.0*t269*t317+2.0*t191*t320*t192-2.0*t7*t320*t186-t327*t329 &
          +2.0*t81*t231*dgx23-2.0*t14*t250-4.0*t26*ddgxr01*t279+4.0*t339*t317+2.0*t93*t342;
    t353 = dg03*dg22*t89;
    t363 = dg33*g33;
    t374 = dgx33*t1;
    t378 = dgx33*r;
    t381 = r*ddgyr23;
    t392 = t378*dg22;
    t394 = t80*g33;
    t398 = dgx22*t1;
    t401 = t80*t1;
    t403 = ddgxr23*g22*g33;
    t406 = t130*dgx01;
    t410 = dg02*g22*t316;
    t415 = t220*t378*dg33-2.0*t161*t381*t10-2.0*t93*t256*dgx23+2.0*t73*dgx01*r*dg23-t93*t392 &
          -t394*dgy22*t1*dg22-t93*t398*dg33-2.0*t401*t403-t406*t23*t10+8.0*t302*t410-16.0*t310*t410;
    t430 = dg33*dg23;
    t431 = t190*t430;
    t433 = g02*g23;
    t434 = t14*t433;
    t436 = g01*ddg23*t89;
    t442 = t104*dg23;
    t445 = t73*ddg23;
    t448 = t315*t186;
    t451 = -t394*t94*dg23+t81*t374*dg23+8.0*t235*g02*dg23*t89+t81*t221*dg22+4.0*t67*t190*t192+t67*t431 &
           +2.0*t434*t436-4.0*t67*t6*t186+2.0*t56*t442-2.0*t56*t445+2.0*t98*t448;
    t453 = t6*ddg23*t10;
    t459 = t260*t192;
    t462 = t374*dg22;
    t464 = t97*t433;
    t466 = dg01*dg23*t89;
    t469 = t1*ddgxr33;
    t475 = r*ddgxr33;
    t483 = t1*ddgyr23;
    t487 = -2.0*t98*t453+4.0*t56*t73*dg23-2.0*t98*t459-t162*t462+2.0*t464*t466-2.0*t161*t469*t10 &
           -2.0*t434*t466-2.0*t220*t475*g33-2.0*t129*t151-t191*t9*t363+2.0*t161*t483*t10;
    t492 = t190*dgx01;
    t496 = t190*ddg23*g33;
    t499 = t6*t430;
    t514 = t130*t100;
    t518 = -t7*t243*t11+t492*t82*g33+2.0*t98*t496+t464*t499+2.0*t129*t445-4.0*t434*g01*dg23*t89 &
           -t434*t499+2.0*t220*t469*g33-2.0*t67*t496+t14*t328*t514-t97*t328*t514;
    t530 = t123*t117;
    t541 = t125*t97;
    t544 = 2.0*t67*t453+2.0*t67*t459-2.0*t464*t436-t98*t431+2.0*t220*t381*g33+t6*t97*t530 &
          -2.0*t67*t448+t162*t392+2.0*t264*t320*dg22+2.0*t7*t265*t186-2.0*t541*t353;
    t568 = -t81*t378*dg23-2.0*t93*t225+t93*t462+2.0*t110*t403-2.0*t81*t256*dgy23+2.0*t401*t112-t162*t342 &
           +4.0*t56*t21*t150+t541*t329-2.0*t191*t265*t192-t6*t14*t530;
    t594 = t56*t177-2.0*t129*t442-t492*t231*g33-t230*t82*t10+t394*dgy22*r*dg22-2.0*t220*t483*g33 &
          +2.0*t161*t475*t10+2.0*t93*t224*dgx23+2.0*t129*t58+t406*t292*t10+t394*t398*dg23;
    p02_rhs = 1/t2/g01*(2.0*t81*t224*dgy23-4.0*t48*r*t279+2.0*t32*t1*ddgxr23+8.0*t202*t203*t167-t81*t213*dg22 &
             -8.0*t216*t203*t150-t220*t221*dg23+t7*t9*t11+t220*t213*dg23-t220*t374*dg33-2.0*t255*t224*t89 &
             -2.0*t32*r*ddgxr23+t22*t23*g22+2.0*t162*t292*dgy23+t191*t243*t363-4.0*t56*t130*t167-2.0*t39*dg01*t16 &
             +t108-2.0*t179*t136+2.0*t48*t49-2.0*t339*t261-2.0*t36*t27-4.0*t14*t17+4.0*t43*t45-2.0*t56*t168 &
             -2.0*t98*t105+t162*t225+2.0*t56*t132+2.0*t26*t27-8.0*t43*t279+2.0*t327*t353+t345+t158+t196+t272 &
             +t309+t415+t451+t487+t544+t518+t594+t568)/(-2.0*r*t10*t89-t16+t49+2.0*t279-t45+t275)/Rmin/2.0
!!!
    t1 = r*r;
    t2 = t1*r;
    t6 = dgy01*g01;
    t7 = g23*g23;
    t8 = t7*t7;
    t11 = Rmin*t2;
    t12 = g03*g23;
    t13 = t11*t12;
    t14 = g33*g01;
    t15 = dg22*dg23;
    t16 = t14*t15;
    t18 = g22*g01;
    t19 = dg33*dg23;
    t20 = t18*t19;
    t23 = g33*g33;
    t24 = g01*g01;
    t25 = t23*t24;
    t29 = g33*t24;
    t30 = t29*g22;
    t31 = r*dg33;
    t35 = dgx22*t1;
    t38 = t1*t1;
    t39 = Rmin*t38;
    t40 = g02*g33;
    t41 = t39*t40;
    t43 = t11*g03;
    t44 = t7*g23;
    t45 = t44*g01;
    t46 = t45*ddg23;
    t49 = t11*g02;
    t50 = t23*g01;
    t51 = dg23*g22;
    t55 = t11*t40;
    t58 = dg23*t7;
    t62 = g01*r;
    t63 = ddgyr01*t8;
    t66 = t29*g23;
    t67 = dgx22*r;
    t70 = t31*dgy22;
    t72 = g22*g22;
    t73 = t72*t23;
    t76 = dgy01*t1;
    t80 = t44*t24;
    t84 = t50*t15;
    t92 = t39*g02;
    t94 = t50*ddg23*g22;
    t97 = -4.0*t49*t14*t58-2.0*t62*t63+t66*t67*dg33+t30*t70+4.0*t6*t73-2.0*t76*dg01*t8+2.0*t80*t1*ddgyr23+t49*t84 &
          -2.0*t80*t1*ddgxr33+2.0*t80*r*ddgxr33+2.0*t92*t94;
    t99 = g23*g01;
    t100 = dg33*dg22;
    t101 = t99*t100;
    t104 = g33*dg01*t58;
    t108 = g01*t1;
    t111 = dgy01*dg01;
    t112 = r*t8;
    t115 = t39*t12;
    t117 = t1*ddgxr23;
    t121 = t1*ddgyr22;
    t126 = t23*dg01*t51;
    t129 = dg01*t44;
    t130 = t129*dg23;
    t133 = t7*g01;
    t134 = t133*t100;
    t137 = -t55*t101+2.0*t92*t104+t41*t101+2.0*t108*t63+2.0*t111*t112+t115*t16+2.0*t29*t117*t7+2.0*t25*t121*g22 &
           -2.0*t92*t126+2.0*t43*t130+2.0*t43*t134;
    t139 = g22*g33;
    t140 = g01*ddg23*t139;
    t144 = dg01*dg23*t139;
    t147 = g23*t24;
    t148 = t147*r;
    t150 = ddgyr23*g22*g33;
    t153 = t147*t1;
    t155 = ddgxr33*g22*g33;
    t162 = dg23*dg23;
    t164 = Rmin*g03*g22*t162;
    t170 = t39*g03;
    t179 = 2.0*t13*t140-2.0*t13*t144+2.0*t148*t150+2.0*t153*t155-2.0*t148*t155+2.0*t14*t2*t164+4.0*t43*t45*dg23 &
          -2.0*t170*t134-2.0*t115*t140+3.0*t115*t20-2.0*t170*t130;
    t182 = t133*t19;
    t186 = dg33*dg33;
    t188 = g23*g22*g01*t186;
    t190 = g02*g23;
    t191 = t39*t190;
    t193 = g01*ddg33*t139;
    t197 = dg01*dg33*t139;
    t202 = g33*g23*g01*t162;
    t206 = t14*ddg23*t7;
    t215 = t147*g22;
    t216 = dgy33*r;
    t222 = t7*t24;
    t224 = dgx33*t1;
    t227 = dgy33*t1;
    t230 = dgx33*r;
    t231 = t230*dg23;
    t234 = r*dg22;
    t238 = t39*dg03;
    t239 = dg01*g22;
    t240 = g33*t7;
    t241 = t239*t240;
    t244 = dg01*t72;
    t245 = t244*t23;
    t249 = dg03*dg01*t8;
    t252 = t11*dg03;
    t257 = -2.0*t153*t150-t222*t70-t215*t224*dg33-t222*t227*dg22-2.0*t222*t231-2.0*t66*t234*dgy23-4.0*t238*t241 &
           +2.0*t238*t245-2.0*t11*t249+4.0*t252*t241-2.0*t252*t245;
    t259 = dg33*t7;
    t260 = t239*t259;
    t263 = t45*ddg33;
    t266 = t129*dg33;
    t269 = t1*dg33;
    t270 = t269*dgy22;
    t274 = t224*dg23;
    t279 = t1*dg23;
    t283 = r*dg23;
    t288 = t1*dg22;
    t292 = 2.0*t170*t260-2.0*t49*t263+2.0*t49*t266+t222*t270+t215*t230*dg33+2.0*t222*t274+t215*t227*dg23 &
          -2.0*t222*t279*dgy23+2.0*t222*t283*dgy23-t92*t84+2.0*t66*t288*dgy23;
    t293 = dg33*g33;
    t294 = t244*t293;
    t301 = dg02*dg33*t139;
    t316 = t139*t7;
    t323 = ddgyr01*t72*t23;
    t326 = -2.0*t170*t294+t222*t216*dg22+2.0*t99*t11*t301+2.0*t170*t46-2.0*t49*t94+2.0*t49*t126-2.0*t49*t202 &
           -4.0*t13*g01*dg23*t139-8.0*t6*t316+4.0*t62*ddgyr01*t316-2.0*t62*t323;
    t330 = t18*Rmin;
    t331 = t38*dg03;
    t334 = t72*g01;
    t335 = t334*Rmin;
    t338 = t99*dgy01;
    t342 = t45*Rmin;
    t343 = t2*dg02;
    t350 = Rmin*t1;
    t351 = t350*g01;
    t353 = dg03*g22*t240;
    t356 = t14*Rmin;
    t357 = t1*g02;
    t361 = t50*Rmin;
    t375 = r*t72*t23;
    t381 = dg03*t72*t23;
    t385 = g01*dg03*t8;
    t388 = t2*dg03;
    t389 = dg22*t7;
    t392 = t14*dgy01;
    t398 = r*ddgyr22;
    t407 = -4.0*t111*r*t316+2.0*t111*t375+2.0*t39*t249+8.0*t351*t381-4.0*t11*t385-t356*t388*t389-t392*t234*t7 &
           -2.0*t29*t121*t7+2.0*t29*t398*t7-t25*t67*dg23+2.0*t66*t283*dgx23;
    t411 = g03*g22*t100;
    t415 = t1*g03;
    t419 = t18*dgy01;
    t422 = t334*dgy01;
    t428 = t99*t350;
    t443 = t334*t186;
    t445 = -t14*t11*t411+8.0*t350*t385-8.0*t335*t415*t293-t419*t31*t7+t422*t31*g33-8.0*t342*t415*dg23 &
           +8.0*t428*g03*dg23*t139-8.0*t342*t357*dg33+8.0*t428*g02*dg33*t139-4.0*t108*ddgyr01*t316+t43*t443;
    t450 = t334*ddg33*g33;
    t456 = t18*ddg33*t7;
    t465 = t11*g01;
    t471 = t38*dg02;
    t475 = -4.0*t43*t18*t259-2.0*t43*t450+2.0*t43*t294+2.0*t43*t456-2.0*t43*t260-t170*t443+2.0*t170*t450-t30*t274 &
           -4.0*t465*t381-2.0*t45*t76*dg23-2.0*t356*t471*t58;
    t481 = dg22*g22;
    t486 = t50*dgy01;
    t504 = 2.0*t361*t471*t51+t361*t388*t481+t356*t331*t389+t486*t234*g22+t14*t39*t411+t392*t288*t7+t30*t231 &
          +4.0*t43*t334*t293+2.0*t108*t323-2.0*t80*r*ddgyr23-t66*t283*dgy22;
    t530 = r*ddgxr23;
    t534 = t419*t269*t7-t422*t269*g33+t330*t388*t259-t335*t388*t293+2.0*t45*dgy01*r*dg23-2.0*t338*t283*t139 &
          -2.0*t92*t266+2.0*t92*t263-2.0*t25*t398*g22-2.0*t25*t117*g22+2.0*t25*t530*g22;
    t547 = t11*t190;
    t562 = -t361*t331*t481+8.0*t465*t353-t66*t35*dg33+t66*t279*dgy22-2.0*t29*t530*t7-4.0*t547*g01*dg33*t139 &
           -t49*t188+2.0*t49*t182+4.0*t49*t45*dg33+2.0*t547*t193-2.0*t547*t197;
    t592 = -2.0*t66*t279*dgx23-t486*t288*g22-2.0*t14*t38*t164-2.0*t170*t456+8.0*t330*t415*t259+2.0*t342*t471*dg33 &
           +2.0*t356*t343*t58-2.0*t361*t343*t51+2.0*t30*t269*dgx23-t25*dgy22*t1*dg22-t30*t270;
    p03_rhs = 1/t2/g01*(t504+t534+t562+t592+t137+t97+t257+t179+t326-2.0*t43*t46-2.0*t92*t206-t41*t20+2.0*t115*t144 &
             -2.0*t49*t104+2.0*t49*t206-3.0*t13*t20-t13*t16+4.0*t6*t8-16.0*t351*t353+t55*t20-2.0*t76*t245+t292 &
             +2.0*t191*t197-t330*t331*t259-2.0*t99*t39*t301+4.0*t76*dg01*t316-8.0*t361*t357*t51+8.0*t356*t357*t58 &
             +4.0*t49*t50*t51+t25*t35*dg23+2.0*t338*t279*t139-t215*t216*dg23+t335*t331*t293-2.0*t342*t343*dg33 &
             -2.0*t30*t31*dgx23+t92*t188-2.0*t92*t182+t25*dgy22*r*dg22-2.0*t191*t193+t475+t445+t407 &
             +2.0*t92*t202)/(-2.0*r*t7*t139-t8+t112+2.0*t316-t73+t375)/Rmin/2.0

      g02_rhs = p02
      g03_rhs = p03

      return

end subroutine pg0a_rhs
!------------------------------------------------------------------------------
subroutine get_g01_rhs(r,g22,g23,g33,dg22,dg23,dg33,ddg22,ddg23,ddg33,g01,g01_rhs)

  implicit none

!~~~~~~% Input parameters:
  real*8,intent(in) :: r,g22,g23,g33,dg22,dg23,dg33,ddg22,ddg23,ddg33,g01
  real*8,intent(out) :: g01_rhs

  real*8 :: t107;
  real*8 :: t11;
  real*8 :: t110;
  real*8 :: t14;
  real*8 :: t19;
  real*8 :: t2;
  real*8 :: t23;
  real*8 :: t25;
  real*8 :: t28;
  real*8 :: t3;
  real*8 :: t33;
  real*8 :: t34;
  real*8 :: t40;
  real*8 :: t45;
  real*8 :: t49;
  real*8 :: t54;
  real*8 :: t6;
  real*8 :: t7;
  real*8 :: t73;
  real*8 :: t76;
  real*8 :: t81;
  real*8 :: t89;
  real*8 :: t98;

    t2 = g23*g23;
    t3 = t2*g23;
    t6 = g22*g33;
    t7 = dg23*dg23;
    t11 = t2*r;
    t14 = g23*r;
    t19 = g22*r;
    t23 = g33*g33;
    t25 = dg22*dg22;
    t28 = r*dg22;
    t33 = g22*g22;
    t34 = t33*g33;
    t40 = g33*r;
    t45 = g22*t23;
    t49 = r*dg33;
    t54 = dg33*dg33;
    t73 = 4.0*r*ddg23*t3-2.0*t6*r*t7-2.0*t11*t7-4.0*t14*ddg23*g22*g33-2.0*t19*ddg33*t2-t23*r*t25 &
         +4.0*g33*g23*t28*dg23+2.0*r*ddg33*t34-2.0*t11*dg33*dg22-2.0*t40*ddg22*t2+2.0*r*ddg22*t45 &
         +4.0*g23*g22*t49*dg23-t33*r*t54-4.0*g33*dg22*t2-4.0*g22*dg33*t2+4.0*dg33*t33*g33 &
         +4.0*dg22*g22*t23-8.0*g23*dg23*t6+8.0*dg23*t3;
    t76 = t2*t2;
    t81 = r*r;
    t89 = dg23*g22*g33;
    t98 = dg33*t2;
    t107 = dg22*t2;
    t110 = -4.0*t76-2.0*r*dg23*t3+2.0*t81*dg23*t3+8.0*t6*t2-t28*t45+2.0*t14*t89-2.0*g23*t81*t89 &
           +t81*dg33*t34-g22*t81*t98-t49*t34+t19*t98+t81*dg22*t45-4.0*t33*t23-g33*t81*t107+t40*t107;
    g01_rhs = g01*t73*(-1.0+r)/t110/2.0

  return

end subroutine get_g01_rhs
!------------------------------------------------------------------------------
