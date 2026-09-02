

#include "macrodef.fh"

! we need only distinguish different finite difference order
! Vertex or Cell is distinguished in routine symmetry_bd which locates in
! file "fmisc.f90"

#if (ghost_width == 2)
! second order code

!-----------------------------------------------------------------------------------------------------------------
!
! General first derivatives of 2_nd oder accurate
!
!              f(i+1) - f(i-1)
!  fx(i) = -----------------------
!                    2 dx
!
!-----------------------------------------------------------------------------------------------------------------

  subroutine fderivs_sh(ex,f,fx,fy,fz,X,Y,Z,SYM1,SYM2,SYM3,symmetry,onoff,sst)
  implicit none

  integer,                               intent(in ):: ex(1:3),symmetry,onoff,sst
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in ):: f
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out):: fx,fy,fz
  real*8,                                intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8,                                intent(in ):: SYM1,SYM2,SYM3

!~~~~~~ other variables

  real*8 :: dX,dY,dZ
  real*8,dimension(0:ex(1)+1,0:ex(2)+1,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8 :: d2dx,d2dy,d2dz
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8,  parameter :: ZEO=0.d0,ONE=1.d0, F60=6.d1
  real*8,  parameter :: TWO=2.d0,EIT=8.d0
  real*8,  parameter ::  F9=9.d0,F45=4.5d1,F12=1.2d1

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = 0
     if(dabs(Y(1)) < dY) jmin = 0
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = 0
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+1
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(1,ex,f,fh,SoA)

  d2dx = ONE/TWO/dX
  d2dy = ONE/TWO/dY
  d2dz = ONE/TWO/dZ

  fx = ZEO
  fy = ZEO
  fz = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
! x direction   
        if(i+1 <= imax .and. i-1 >= imin)then
!
!              - f(i-1) + f(i+1)
!  fx(i) = --------------------------------
!                     2 dx
      fx(i,j,k)=d2dx*(-fh(i-1,j,k)+fh(i+1,j,k))

! set imax and imin 0
    endif
! y direction   
        if(j+1 <= jmax .and. j-1 >= jmin)then

     fy(i,j,k)=d2dy*(-fh(i,j-1,k)+fh(i,j+1,k))

! set jmax and jmin 0
    endif
! z direction   
        if(k+1 <= kmax .and. k-1 >= kmin)then

      fz(i,j,k)=d2dz*(-fh(i,j,k-1)+fh(i,j,k+1))

! set kmax and kmin 0
    endif

  enddo
  enddo
  enddo

  return

  end subroutine fderivs_sh
!-----------------------------------------------------------------------------
!
! single derivatives dx
!
!-----------------------------------------------------------------------------
  subroutine fdx_sh(ex,f,fx,X,Y,Z,SYM1,SYM2,SYM3,symmetry,onoff,sst)
  implicit none

  integer,                               intent(in ):: ex(1:3),symmetry,onoff,sst
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in ):: f
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out):: fx
  real*8,                                intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8,                                intent(in ):: SYM1,SYM2,SYM3

!~~~~~~ other variables

  real*8 :: dX,dY,dZ
  real*8,dimension(0:ex(1)+1,0:ex(2)+1,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8 :: d2dx
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8,  parameter :: ZEO=0.d0,ONE=1.d0, F60=6.d1
  real*8,  parameter :: TWO=2.d0,EIT=8.d0
  real*8,  parameter ::  F9=9.d0,F45=4.5d1,F12=1.2d1

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)
  
  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = 0
     if(dabs(Y(1)) < dY) jmin = 0
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = 0
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+1
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(1,ex,f,fh,SoA)

  d2dx = ONE/TWO/dX

  fx = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
! x direction   
        if(i+1 <= imax .and. i-1 >= imin)then
!
!              - f(i-1) + f(i+1)
!  fx(i) = --------------------------------
!                     2 dx
      fx(i,j,k)=d2dx*(-fh(i-1,j,k)+fh(i+1,j,k))

! set imax and imin 0
    endif

  enddo
  enddo
  enddo

  return

  end subroutine fdx_sh
!-----------------------------------------------------------------------------
!
! single derivatives dy
!
!-----------------------------------------------------------------------------
  subroutine fdy_sh(ex,f,fy,X,Y,Z,SYM1,SYM2,SYM3,symmetry,onoff,sst)
  implicit none

  integer,                               intent(in ):: ex(1:3),symmetry,onoff,sst
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in ):: f
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out):: fy
  real*8,                                intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8,                                intent(in ):: SYM1,SYM2,SYM3

!~~~~~~ other variables

  real*8 :: dX,dY,dZ
  real*8,dimension(0:ex(1)+1,0:ex(2)+1,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8 :: d2dy
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8,  parameter :: ZEO=0.d0,ONE=1.d0, F60=6.d1
  real*8,  parameter :: TWO=2.d0,EIT=8.d0
  real*8,  parameter ::  F9=9.d0,F45=4.5d1,F12=1.2d1

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)
  
  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1
  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = 0
     if(dabs(Y(1)) < dY) jmin = 0
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = 0
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+1
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(1,ex,f,fh,SoA)

  d2dy = ONE/TWO/dY

  fy = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
! y direction   
        if(j+1 <= jmax .and. j-1 >= jmin)then

     fy(i,j,k)=d2dy*(-fh(i,j-1,k)+fh(i,j+1,k))

! set jmax and jmin 0
    endif

  enddo
  enddo
  enddo

  return

  end subroutine fdy_sh
!-----------------------------------------------------------------------------
!
! single derivatives dz
!
!-----------------------------------------------------------------------------
  subroutine fdz_sh(ex,f,fz,X,Y,Z,SYM1,SYM2,SYM3,symmetry,onoff,sst)
  implicit none

  integer,                               intent(in ):: ex(1:3),symmetry,onoff,sst
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in ):: f
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out):: fz
  real*8,                                intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8,                                intent(in ):: SYM1,SYM2,SYM3

!~~~~~~ other variables
  
  real*8 :: dX,dY,dZ
  real*8,dimension(0:ex(1)+1,0:ex(2)+1,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8 :: d2dz
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8,  parameter :: ZEO=0.d0,ONE=1.d0, F60=6.d1
  real*8,  parameter :: TWO=2.d0,EIT=8.d0
  real*8,  parameter ::  F9=9.d0,F45=4.5d1,F12=1.2d1

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)
  
  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1
  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = 0
     if(dabs(Y(1)) < dY) jmin = 0
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = 0
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+1
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(1,ex,f,fh,SoA)

  d2dz = ONE/TWO/dZ

  fz = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
! z direction   
        if(k+1 <= kmax .and. k-1 >= kmin)then

      fz(i,j,k)=d2dz*(-fh(i,j,k-1)+fh(i,j,k+1))

! set kmax and kmin 0
    endif

  enddo
  enddo
  enddo

  return

  end subroutine fdz_sh
!-----------------------------------------------------------------------------------------------------------------
!
! General second derivatives of 2_nd oder accurate
!
!               f(i-1) - 2 f(i) + f(i+1)
!  fxx(i) = --------------------------------
!                        dx^2 
!
!                 f(i-1,j-1) - f(i+1,j-1) - f(i-1,j+1) + f(i+1,j+1) 
!  fxy(i,j) = -----------------------------------------------------------
!                                      4 dx dy
!
!-----------------------------------------------------------------------------------------------------------------
  subroutine fdderivs_sh(ex,f,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z, &
                      SYM1,SYM2,SYM3,symmetry,onoff,sst)
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry,onoff,sst
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fxx,fxy,fxz,fyy,fyz,fzz
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8,dimension(0:ex(1)+1,0:ex(2)+1,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdxdx,Sdydy,Sdzdz
  real*8  :: Sdxdy,Sdxdz,Sdydz
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = 0
     if(dabs(Y(1)) < dY) jmin = 0
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = 0
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+1
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(1,ex,f,fh,SoA)

  Sdxdx =  ONE /( dX * dX )
  Sdydy =  ONE /( dY * dY )
  Sdzdz =  ONE /( dZ * dZ )

  Sdxdy = F1o4 /( dX * dY )
  Sdxdz = F1o4 /( dX * dZ )
  Sdydz = F1o4 /( dY * dZ )

  fxx = ZEO
  fyy = ZEO
  fzz = ZEO
  fxy = ZEO
  fxz = ZEO
  fyz = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
!~~~~~~ fxx
        if(i+1 <= imax .and. i-1 >= imin)then
!
!               f(i-1) - 2 f(i) + f(i+1)
!  fxx(i) = --------------------------------
!                         dx^2 
   fxx(i,j,k) = Sdxdx*(fh(i-1,j,k)-TWO*fh(i,j,k) &
                      +fh(i+1,j,k)              )
   endif


!~~~~~~ fyy
       if(j+1 <= jmax .and. j-1 >= jmin)then

   fyy(i,j,k) = Sdydy*(fh(i,j-1,k)-TWO*fh(i,j,k) &
                      +fh(i,j+1,k)              )
   endif

!~~~~~~ fzz
       if(k+1 <= kmax .and. k-1 >= kmin)then

   fzz(i,j,k) = Sdzdz*(fh(i,j,k-1)-TWO*fh(i,j,k) &
                      +fh(i,j,k+1)              )
   endif
!~~~~~~ fxy
      if(i+1 <= imax .and. i-1 >= imin .and. j+1 <= jmax .and. j-1 >= jmin)then
!                 f(i-1,j-1) - f(i+1,j-1) - f(i-1,j+1) + f(i+1,j+1) 
!  fxy(i,j) = -----------------------------------------------------------
!                                      4 dx dy
   fxy(i,j,k) = Sdxdy*(fh(i-1,j-1,k)-fh(i+1,j-1,k)-fh(i-1,j+1,k)+fh(i+1,j+1,k))
   endif
!~~~~~~ fxz
      if(i+1 <= imax .and. i-1 >= imin .and. k+1 <= kmax .and. k-1 >= kmin)then
   fxz(i,j,k) = Sdxdz*(fh(i-1,j,k-1)-fh(i+1,j,k-1)-fh(i-1,j,k+1)+fh(i+1,j,k+1))
   endif
!~~~~~~ fyz
      if(j+1 <= jmax .and. j-1 >= jmin .and. k+1 <= kmax .and. k-1 >= kmin)then
   fyz(i,j,k) = Sdydz*(fh(i,j-1,k-1)-fh(i,j+1,k-1)-fh(i,j-1,k+1)+fh(i,j+1,k+1))
   endif 

   enddo
   enddo
   enddo

  return

  end subroutine fdderivs_sh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! only for compute_ricci.f90 usage
!-----------------------------------------------------------------------------
  subroutine fddxx_sh(ex,f,fxx,X,Y,Z,SYM1,SYM2,SYM3,symmetry,sst)
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry,sst
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fxx
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8,dimension(0:ex(1)+1,0:ex(2)+1,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdxdx
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = 0
     if(dabs(Y(1)) < dY) jmin = 0
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = 0
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+1
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(1,ex,f,fh,SoA)

  Sdxdx =  ONE /( dX * dX )

  fxx = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
!~~~~~~ fxx
       if(i+1 <= imax .and. i-1 >= imin)then
   fxx(i,j,k) = Sdxdx*(fh(i-1,j,k)-TWO*fh(i,j,k) &
                      +fh(i+1,j,k)              )
   endif

   enddo
   enddo
   enddo

  return

  end subroutine fddxx_sh

  subroutine fddyy_sh(ex,f,fyy,X,Y,Z,SYM1,SYM2,SYM3,symmetry,sst)
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry,sst
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fyy
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8,dimension(0:ex(1)+1,0:ex(2)+1,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdydy
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = 0
     if(dabs(Y(1)) < dY) jmin = 0
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = 0
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+1
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(1,ex,f,fh,SoA)

  Sdydy =  ONE /( dY * dY )

  fyy = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
!~~~~~~ fyy
       if(j+1 <= jmax .and. j-1 >= jmin)then

   fyy(i,j,k) = Sdydy*(fh(i,j-1,k)-TWO*fh(i,j,k) &
                      +fh(i,j+1,k)              )
   endif

   enddo
   enddo
   enddo

  return

  end subroutine fddyy_sh

  subroutine fddzz_sh(ex,f,fzz,X,Y,Z,SYM1,SYM2,SYM3,symmetry,sst)
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry,sst
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fzz
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8,dimension(0:ex(1)+1,0:ex(2)+1,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdzdz
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = 0
     if(dabs(Y(1)) < dY) jmin = 0
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = 0
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+1
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(1,ex,f,fh,SoA)

  Sdzdz =  ONE /( dZ * dZ )

  fzz = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
!~~~~~~ fzz
       if(k+1 <= kmax .and. k-1 >= kmin)then

   fzz(i,j,k) = Sdzdz*(fh(i,j,k-1)-TWO*fh(i,j,k) &
                      +fh(i,j,k+1)              )
   endif

   enddo
   enddo
   enddo

  return

  end subroutine fddzz_sh

  subroutine fddxy_sh(ex,f,fxy,X,Y,Z,SYM1,SYM2,SYM3,symmetry,sst)
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry,sst
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fxy
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8,dimension(0:ex(1)+1,0:ex(2)+1,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdxdy
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = 0
     if(dabs(Y(1)) < dY) jmin = 0
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = 0
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+1
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(1,ex,f,fh,SoA)

  Sdxdy = F1o4 /( dX * dY )

  fxy = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
!~~~~~~ fxy
      if(i+1 <= imax .and. i-1 >= imin .and. j+1 <= jmax .and. j-1 >= jmin)then

   fxy(i,j,k) = Sdxdy*(fh(i-1,j-1,k)-fh(i+1,j-1,k)-fh(i-1,j+1,k)+fh(i+1,j+1,k))
   endif

   enddo
   enddo
   enddo

  return

  end subroutine fddxy_sh

  subroutine fddxz_sh(ex,f,fxz,X,Y,Z,SYM1,SYM2,SYM3,symmetry,sst)
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry,sst
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fxz
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8,dimension(0:ex(1)+1,0:ex(2)+1,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdxdz
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = 0
     if(dabs(Y(1)) < dY) jmin = 0
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = 0
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+1
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(1,ex,f,fh,SoA)

  Sdxdz = F1o4 /( dX * dZ )

  fxz = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
!~~~~~~ fxz
      if(i+1 <= imax .and. i-1 >= imin .and. k+1 <= kmax .and. k-1 >= kmin)then
   fxz(i,j,k) = Sdxdz*(fh(i-1,j,k-1)-fh(i+1,j,k-1)-fh(i-1,j,k+1)+fh(i+1,j,k+1))
   endif

   enddo
   enddo
   enddo

  return

  end subroutine fddxz_sh

  subroutine fddyz_sh(ex,f,fyz,X,Y,Z,SYM1,SYM2,SYM3,symmetry,sst)
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry,sst
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fyz
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8,dimension(0:ex(1)+1,0:ex(2)+1,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdydz
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = 0
     if(dabs(Y(1)) < dY) jmin = 0
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = 0
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+1
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(1,ex,f,fh,SoA)

  Sdydz = F1o4 /( dY * dZ )

  fyz = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
!~~~~~~ fyz
      if(j+1 <= jmax .and. j-1 >= jmin .and. k+1 <= kmax .and. k-1 >= kmin)then
   fyz(i,j,k) = Sdydz*(fh(i,j-1,k-1)-fh(i,j+1,k-1)-fh(i,j-1,k+1)+fh(i,j+1,k+1))
   endif 

   enddo
   enddo
   enddo

  return

  end subroutine fddyz_sh

#elif (ghost_width == 3)
! fourth order code

!-----------------------------------------------------------------------------------------------------------------
!
! General first derivatives of 4_th oder accurate
!
!              f(i-2) - 8 f(i-1) + 8 f(i+1) - f(i+2)
!  fx(i) = ---------------------------------------------
!                             12 dx
!
!-----------------------------------------------------------------------------------------------------------------

  subroutine fderivs_sh(ex,f,fx,fy,fz,X,Y,Z,SYM1,SYM2,SYM3,symmetry,onoff,sst)
  implicit none

  integer,                               intent(in ):: ex(1:3),symmetry,onoff,sst
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in ):: f
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out):: fx,fy,fz
  real*8,                                intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8,                                intent(in ):: SYM1,SYM2,SYM3

!~~~~~~ other variables

  real*8 :: dX,dY,dZ
  real*8,dimension(-1:ex(1)+2,-1:ex(2)+2,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8 :: d12dx,d12dy,d12dz,d2dx,d2dy,d2dz
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8,  parameter :: ZEO=0.d0,ONE=1.d0, F60=6.d1
  real*8,  parameter :: TWO=2.d0,EIT=8.d0
  real*8,  parameter ::  F9=9.d0,F45=4.5d1,F12=1.2d1

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -1
     if(dabs(Y(1)) < dY) jmin = -1
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -1
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+2
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(2,ex,f,fh,SoA)

  d12dx = ONE/F12/dX
  d12dy = ONE/F12/dY
  d12dz = ONE/F12/dZ

  d2dx = ONE/TWO/dX
  d2dy = ONE/TWO/dY
  d2dz = ONE/TWO/dZ

  fx = ZEO
  fy = ZEO
  fz = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
! x direction   
        if(i+2 <= imax .and. i-2 >= imin)then
!
!              f(i-2) - 8 f(i-1) + 8 f(i+1) - f(i+2)
!  fx(i) = ---------------------------------------------
!                             12 dx
      fx(i,j,k)=d12dx*(fh(i-2,j,k)-EIT*fh(i-1,j,k)+EIT*fh(i+1,j,k)-fh(i+2,j,k))

    elseif(i+1 <= imax .and. i-1 >= imin)then
!
!              - f(i-1) + f(i+1)
!  fx(i) = --------------------------------
!                     2 dx
      fx(i,j,k)=d2dx*(-fh(i-1,j,k)+fh(i+1,j,k))

! set imax and imin 0
    endif
! y direction   
        if(j+2 <= jmax .and. j-2 >= jmin)then

      fy(i,j,k)=d12dy*(fh(i,j-2,k)-EIT*fh(i,j-1,k)+EIT*fh(i,j+1,k)-fh(i,j+2,k))

    elseif(j+1 <= jmax .and. j-1 >= jmin)then

     fy(i,j,k)=d2dy*(-fh(i,j-1,k)+fh(i,j+1,k))

! set jmax and jmin 0
    endif
! z direction   
        if(k+2 <= kmax .and. k-2 >= kmin)then

      fz(i,j,k)=d12dz*(fh(i,j,k-2)-EIT*fh(i,j,k-1)+EIT*fh(i,j,k+1)-fh(i,j,k+2))

    elseif(k+1 <= kmax .and. k-1 >= kmin)then

      fz(i,j,k)=d2dz*(-fh(i,j,k-1)+fh(i,j,k+1))

! set kmax and kmin 0
    endif

  enddo
  enddo
  enddo

  return

  end subroutine fderivs_sh
!-----------------------------------------------------------------------------
!
! single derivatives dx
!
!-----------------------------------------------------------------------------
  subroutine fdx_sh(ex,f,fx,X,Y,Z,SYM1,SYM2,SYM3,symmetry,onoff,sst)
  implicit none

  integer,                               intent(in ):: ex(1:3),symmetry,onoff,sst
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in ):: f
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out):: fx
  real*8,                                intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8,                                intent(in ):: SYM1,SYM2,SYM3

!~~~~~~ other variables

  real*8 :: dX,dY,dZ
  real*8,dimension(-1:ex(1)+2,-1:ex(2)+2,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8 :: d12dx,d2dx
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8,  parameter :: ZEO=0.d0,ONE=1.d0, F60=6.d1
  real*8,  parameter :: TWO=2.d0,EIT=8.d0
  real*8,  parameter ::  F9=9.d0,F45=4.5d1,F12=1.2d1

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)
  
  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -1
     if(dabs(Y(1)) < dY) jmin = -1
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -1
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+2
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(2,ex,f,fh,SoA)

  d12dx = ONE/F12/dX

  d2dx = ONE/TWO/dX

  fx = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
! x direction   
        if(i+2 <= imax .and. i-2 >= imin)then
!
!              f(i-2) - 8 f(i-1) + 8 f(i+1) - f(i+2)
!  fx(i) = ---------------------------------------------
!                             12 dx
      fx(i,j,k)=d12dx*(fh(i-2,j,k)-EIT*fh(i-1,j,k)+EIT*fh(i+1,j,k)-fh(i+2,j,k))

    elseif(i+1 <= imax .and. i-1 >= imin)then
!
!              - f(i-1) + f(i+1)
!  fx(i) = --------------------------------
!                     2 dx
      fx(i,j,k)=d2dx*(-fh(i-1,j,k)+fh(i+1,j,k))

! set imax and imin 0
    endif

  enddo
  enddo
  enddo

  return

  end subroutine fdx_sh
!-----------------------------------------------------------------------------
!
! single derivatives dy
!
!-----------------------------------------------------------------------------
  subroutine fdy_sh(ex,f,fy,X,Y,Z,SYM1,SYM2,SYM3,symmetry,onoff,sst)
  implicit none

  integer,                               intent(in ):: ex(1:3),symmetry,onoff,sst
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in ):: f
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out):: fy
  real*8,                                intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8,                                intent(in ):: SYM1,SYM2,SYM3

!~~~~~~ other variables

  real*8 :: dX,dY,dZ
  real*8,dimension(-1:ex(1)+2,-1:ex(2)+2,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8 :: d12dy,d2dy
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8,  parameter :: ZEO=0.d0,ONE=1.d0, F60=6.d1
  real*8,  parameter :: TWO=2.d0,EIT=8.d0
  real*8,  parameter ::  F9=9.d0,F45=4.5d1,F12=1.2d1

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)
  
  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -1
     if(dabs(Y(1)) < dY) jmin = -1
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -1
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+2
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(2,ex,f,fh,SoA)

  d12dy = ONE/F12/dY

  d2dy = ONE/TWO/dY

  fy = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
! y direction   
        if(j+2 <= jmax .and. j-2 >= jmin)then

      fy(i,j,k)=d12dy*(fh(i,j-2,k)-EIT*fh(i,j-1,k)+EIT*fh(i,j+1,k)-fh(i,j+2,k))

    elseif(j+1 <= jmax .and. j-1 >= jmin)then

     fy(i,j,k)=d2dy*(-fh(i,j-1,k)+fh(i,j+1,k))

! set jmax and jmin 0
    endif

  enddo
  enddo
  enddo

  return

  end subroutine fdy_sh
!-----------------------------------------------------------------------------
!
! single derivatives dz
!
!-----------------------------------------------------------------------------
  subroutine fdz_sh(ex,f,fz,X,Y,Z,SYM1,SYM2,SYM3,symmetry,onoff,sst)
  implicit none

  integer,                               intent(in ):: ex(1:3),symmetry,onoff,sst
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in ):: f
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out):: fz
  real*8,                                intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8,                                intent(in ):: SYM1,SYM2,SYM3

!~~~~~~ other variables
  
  real*8 :: dX,dY,dZ
  real*8,dimension(-1:ex(1)+2,-1:ex(2)+2,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8 :: d12dz,d2dz
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8,  parameter :: ZEO=0.d0,ONE=1.d0, F60=6.d1
  real*8,  parameter :: TWO=2.d0,EIT=8.d0
  real*8,  parameter ::  F9=9.d0,F45=4.5d1,F12=1.2d1

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)
  
  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -1
     if(dabs(Y(1)) < dY) jmin = -1
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -1
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+2
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(2,ex,f,fh,SoA)

  d12dz = ONE/F12/dZ

  d2dz = ONE/TWO/dZ

  fz = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
! z direction   
        if(k+2 <= kmax .and. k-2 >= kmin)then

      fz(i,j,k)=d12dz*(fh(i,j,k-2)-EIT*fh(i,j,k-1)+EIT*fh(i,j,k+1)-fh(i,j,k+2))

    elseif(k+1 <= kmax .and. k-1 >= kmin)then

      fz(i,j,k)=d2dz*(-fh(i,j,k-1)+fh(i,j,k+1))

! set kmax and kmin 0
    endif

  enddo
  enddo
  enddo

  return

  end subroutine fdz_sh
!-----------------------------------------------------------------------------------------------------------------
!
! General second derivatives of 4_th oder accurate
!
!               - f(i-2) + 16 f(i-1) - 30 f(i) + 16 f(i+1) - f(i+2)
!  fxx(i) = ----------------------------------------------------------
!                                  12 dx^2 
!
!             -   ( - f(i+2,j+2) + 8 f(i+1,j+2) - 8 f(i-1,j+2) + f(i-2,j+2) )
!             + 8 ( - f(i+2,j+1) + 8 f(i+1,j+1) - 8 f(i-1,j+1) + f(i-2,j+1) )
!             - 8 ( - f(i+2,j-1) + 8 f(i+1,j-1) - 8 f(i-1,j-1) + f(i-2,j-1) )
!             +   ( - f(i+2,j-2) + 8 f(i+1,j-2) - 8 f(i-1,j-2) + f(i-2,j-2) )
!  fxy(i,j) = ----------------------------------------------------------------
!                                  144 dx dy
!
!-----------------------------------------------------------------------------------------------------------------
  subroutine fdderivs_sh(ex,f,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z, &
                      SYM1,SYM2,SYM3,symmetry,onoff,sst)
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry,onoff,sst
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fxx,fxy,fxz,fyy,fyz,fzz
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8,dimension(-1:ex(1)+2,-1:ex(2)+2,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdxdx,Sdydy,Sdzdz,Fdxdx,Fdydy,Fdzdz
  real*8  :: Sdxdy,Sdxdz,Sdydz,Fdxdy,Fdxdz,Fdydz
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -1
     if(dabs(Y(1)) < dY) jmin = -1
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -1
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+2
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(2,ex,f,fh,SoA)

  Sdxdx =  ONE /( dX * dX )
  Sdydy =  ONE /( dY * dY )
  Sdzdz =  ONE /( dZ * dZ )

  Fdxdx = F1o12 /( dX * dX )
  Fdydy = F1o12 /( dY * dY )
  Fdzdz = F1o12 /( dZ * dZ )

  Sdxdy = F1o4 /( dX * dY )
  Sdxdz = F1o4 /( dX * dZ )
  Sdydz = F1o4 /( dY * dZ )

  Fdxdy = F1o144 /( dX * dY )
  Fdxdz = F1o144 /( dX * dZ )
  Fdydz = F1o144 /( dY * dZ )

  fxx = ZEO
  fyy = ZEO
  fzz = ZEO
  fxy = ZEO
  fxz = ZEO
  fyz = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
!~~~~~~ fxx
        if(i+2 <= imax .and. i-2 >= imin)then
!
!               - f(i-2) + 16 f(i-1) - 30 f(i) + 16 f(i+1) - f(i+2)
!  fxx(i) = ----------------------------------------------------------
!                                  12 dx^2 
   fxx(i,j,k) = Fdxdx*(-fh(i-2,j,k)+F16*fh(i-1,j,k)-F30*fh(i,j,k) &
                       -fh(i+2,j,k)+F16*fh(i+1,j,k)              )
   elseif(i+1 <= imax .and. i-1 >= imin)then
!
!               f(i-1) - 2 f(i) + f(i+1)
!  fxx(i) = --------------------------------
!                         dx^2 
   fxx(i,j,k) = Sdxdx*(fh(i-1,j,k)-TWO*fh(i,j,k) &
                      +fh(i+1,j,k)              )
   endif


!~~~~~~ fyy
        if(j+2 <= jmax .and. j-2 >= jmin)then

   fyy(i,j,k) = Fdydy*(-fh(i,j-2,k)+F16*fh(i,j-1,k)-F30*fh(i,j,k) &
                       -fh(i,j+2,k)+F16*fh(i,j+1,k)              )
   elseif(j+1 <= jmax .and. j-1 >= jmin)then

   fyy(i,j,k) = Sdydy*(fh(i,j-1,k)-TWO*fh(i,j,k) &
                      +fh(i,j+1,k)              )
   endif

!~~~~~~ fzz
        if(k+2 <= kmax .and. k-2 >= kmin)then

   fzz(i,j,k) = Fdzdz*(-fh(i,j,k-2)+F16*fh(i,j,k-1)-F30*fh(i,j,k) &
                       -fh(i,j,k+2)+F16*fh(i,j,k+1)              )
   elseif(k+1 <= kmax .and. k-1 >= kmin)then

   fzz(i,j,k) = Sdzdz*(fh(i,j,k-1)-TWO*fh(i,j,k) &
                      +fh(i,j,k+1)              )
   endif
!~~~~~~ fxy
       if(i+2 <= imax .and. i-2 >= imin .and. j+2 <= jmax .and. j-2 >= jmin)then
!
!                 ( f(i-2,j-2) - 8 f(i-1,j-2) + 8 f(i+1,j-2) - f(i+2,j-2) )
!             - 8 ( f(i-2,j-1) - 8 f(i-1,j-1) + 8 f(i+1,j-1) - f(i+2,j-1) )
!             + 8 ( f(i-2,j+1) - 8 f(i-1,j+1) + 8 f(i+1,j+1) - f(i+2,j+1) )
!             -   ( f(i-2,j+2) - 8 f(i-1,j+2) + 8 f(i+1,j+2) - f(i+2,j+2) )
!  fxy(i,j) = ----------------------------------------------------------------
!                                  144 dx dy
   fxy(i,j,k) = Fdxdy*(     (fh(i-2,j-2,k)-F8*fh(i-1,j-2,k)+F8*fh(i+1,j-2,k)-fh(i+2,j-2,k))  &
                       -F8 *(fh(i-2,j-1,k)-F8*fh(i-1,j-1,k)+F8*fh(i+1,j-1,k)-fh(i+2,j-1,k))  &
                       +F8 *(fh(i-2,j+1,k)-F8*fh(i-1,j+1,k)+F8*fh(i+1,j+1,k)-fh(i+2,j+1,k))  &
                       -    (fh(i-2,j+2,k)-F8*fh(i-1,j+2,k)+F8*fh(i+1,j+2,k)-fh(i+2,j+2,k)))

   elseif(i+1 <= imax .and. i-1 >= imin .and. j+1 <= jmax .and. j-1 >= jmin)then
!                 f(i-1,j-1) - f(i+1,j-1) - f(i-1,j+1) + f(i+1,j+1) 
!  fxy(i,j) = -----------------------------------------------------------
!                                      4 dx dy
   fxy(i,j,k) = Sdxdy*(fh(i-1,j-1,k)-fh(i+1,j-1,k)-fh(i-1,j+1,k)+fh(i+1,j+1,k))
   endif
!~~~~~~ fxz
       if(i+2 <= imax .and. i-2 >= imin .and. k+2 <= kmax .and. k-2 >= kmin)then
   fxz(i,j,k) = Fdxdz*(     (fh(i-2,j,k-2)-F8*fh(i-1,j,k-2)+F8*fh(i+1,j,k-2)-fh(i+2,j,k-2))  &
                       -F8 *(fh(i-2,j,k-1)-F8*fh(i-1,j,k-1)+F8*fh(i+1,j,k-1)-fh(i+2,j,k-1))  &
                       +F8 *(fh(i-2,j,k+1)-F8*fh(i-1,j,k+1)+F8*fh(i+1,j,k+1)-fh(i+2,j,k+1))  &
                       -    (fh(i-2,j,k+2)-F8*fh(i-1,j,k+2)+F8*fh(i+1,j,k+2)-fh(i+2,j,k+2)))
   elseif(i+1 <= imax .and. i-1 >= imin .and. k+1 <= kmax .and. k-1 >= kmin)then
   fxz(i,j,k) = Sdxdz*(fh(i-1,j,k-1)-fh(i+1,j,k-1)-fh(i-1,j,k+1)+fh(i+1,j,k+1))
   endif
!~~~~~~ fyz
       if(j+2 <= jmax .and. j-2 >= jmin .and. k+2 <= kmax .and. k-2 >= kmin)then
   fyz(i,j,k) = Fdydz*(     (fh(i,j-2,k-2)-F8*fh(i,j-1,k-2)+F8*fh(i,j+1,k-2)-fh(i,j+2,k-2))  &
                       -F8 *(fh(i,j-2,k-1)-F8*fh(i,j-1,k-1)+F8*fh(i,j+1,k-1)-fh(i,j+2,k-1))  &
                       +F8 *(fh(i,j-2,k+1)-F8*fh(i,j-1,k+1)+F8*fh(i,j+1,k+1)-fh(i,j+2,k+1))  &
                       -    (fh(i,j-2,k+2)-F8*fh(i,j-1,k+2)+F8*fh(i,j+1,k+2)-fh(i,j+2,k+2)))
   elseif(j+1 <= jmax .and. j-1 >= jmin .and. k+1 <= kmax .and. k-1 >= kmin)then
   fyz(i,j,k) = Sdydz*(fh(i,j-1,k-1)-fh(i,j+1,k-1)-fh(i,j-1,k+1)+fh(i,j+1,k+1))
   endif 

   enddo
   enddo
   enddo

  return

  end subroutine fdderivs_sh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! only for compute_ricci.f90 usage
!-----------------------------------------------------------------------------
  subroutine fddxx_sh(ex,f,fxx,X,Y,Z,SYM1,SYM2,SYM3,symmetry,sst)
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry,sst
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fxx
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8,dimension(-1:ex(1)+2,-1:ex(2)+2,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdxdx,Fdxdx
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -1
     if(dabs(Y(1)) < dY) jmin = -1
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -1
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+2
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(2,ex,f,fh,SoA)

  Sdxdx =  ONE /( dX * dX )

  Fdxdx = F1o12 /( dX * dX )

  fxx = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
!~~~~~~ fxx
        if(i+2 <= imax .and. i-2 >= imin)then
   fxx(i,j,k) = Fdxdx*(-fh(i-2,j,k)+F16*fh(i-1,j,k)-F30*fh(i,j,k) &
                       -fh(i+2,j,k)+F16*fh(i+1,j,k)              )
   elseif(i+1 <= imax .and. i-1 >= imin)then
   fxx(i,j,k) = Sdxdx*(fh(i-1,j,k)-TWO*fh(i,j,k) &
                      +fh(i+1,j,k)              )
   endif

   enddo
   enddo
   enddo

  return

  end subroutine fddxx_sh

  subroutine fddyy_sh(ex,f,fyy,X,Y,Z,SYM1,SYM2,SYM3,symmetry,sst)
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry,sst
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fyy
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8,dimension(-1:ex(1)+2,-1:ex(2)+2,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdydy,Fdydy
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -1
     if(dabs(Y(1)) < dY) jmin = -1
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -1
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+2
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(2,ex,f,fh,SoA)

  Sdydy =  ONE /( dY * dY )

  Fdydy = F1o12 /( dY * dY )

  fyy = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
!~~~~~~ fyy
        if(j+2 <= jmax .and. j-2 >= jmin)then

   fyy(i,j,k) = Fdydy*(-fh(i,j-2,k)+F16*fh(i,j-1,k)-F30*fh(i,j,k) &
                       -fh(i,j+2,k)+F16*fh(i,j+1,k)              )
   elseif(j+1 <= jmax .and. j-1 >= jmin)then

   fyy(i,j,k) = Sdydy*(fh(i,j-1,k)-TWO*fh(i,j,k) &
                      +fh(i,j+1,k)              )
   endif

   enddo
   enddo
   enddo

  return

  end subroutine fddyy_sh

  subroutine fddzz_sh(ex,f,fzz,X,Y,Z,SYM1,SYM2,SYM3,symmetry,sst)
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry,sst
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fzz
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8,dimension(-1:ex(1)+2,-1:ex(2)+2,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdzdz,Fdzdz
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -1
     if(dabs(Y(1)) < dY) jmin = -1
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -1
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+2
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(2,ex,f,fh,SoA)

  Sdzdz =  ONE /( dZ * dZ )

  Fdzdz = F1o12 /( dZ * dZ )

  fzz = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
!~~~~~~ fzz
        if(k+2 <= kmax .and. k-2 >= kmin)then

   fzz(i,j,k) = Fdzdz*(-fh(i,j,k-2)+F16*fh(i,j,k-1)-F30*fh(i,j,k) &
                       -fh(i,j,k+2)+F16*fh(i,j,k+1)              )
   elseif(k+1 <= kmax .and. k-1 >= kmin)then

   fzz(i,j,k) = Sdzdz*(fh(i,j,k-1)-TWO*fh(i,j,k) &
                      +fh(i,j,k+1)              )
   endif

   enddo
   enddo
   enddo

  return

  end subroutine fddzz_sh

  subroutine fddxy_sh(ex,f,fxy,X,Y,Z,SYM1,SYM2,SYM3,symmetry,sst)
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry,sst
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fxy
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8,dimension(-1:ex(1)+2,-1:ex(2)+2,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdxdy,Fdxdy
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -1
     if(dabs(Y(1)) < dY) jmin = -1
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -1
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+2
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(2,ex,f,fh,SoA)

  Sdxdy = F1o4 /( dX * dY )

  Fdxdy = F1o144 /( dX * dY )

  fxy = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
!~~~~~~ fxy
       if(i+2 <= imax .and. i-2 >= imin .and. j+2 <= jmax .and. j-2 >= jmin)then

   fxy(i,j,k) = Fdxdy*(     (fh(i-2,j-2,k)-F8*fh(i-1,j-2,k)+F8*fh(i+1,j-2,k)-fh(i+2,j-2,k))  &
                       -F8 *(fh(i-2,j-1,k)-F8*fh(i-1,j-1,k)+F8*fh(i+1,j-1,k)-fh(i+2,j-1,k))  &
                       +F8 *(fh(i-2,j+1,k)-F8*fh(i-1,j+1,k)+F8*fh(i+1,j+1,k)-fh(i+2,j+1,k))  &
                       -    (fh(i-2,j+2,k)-F8*fh(i-1,j+2,k)+F8*fh(i+1,j+2,k)-fh(i+2,j+2,k)))
   elseif(i+1 <= imax .and. i-1 >= imin .and. j+1 <= jmax .and. j-1 >= jmin)then

   fxy(i,j,k) = Sdxdy*(fh(i-1,j-1,k)-fh(i+1,j-1,k)-fh(i-1,j+1,k)+fh(i+1,j+1,k))
   endif

   enddo
   enddo
   enddo

  return

  end subroutine fddxy_sh

  subroutine fddxz_sh(ex,f,fxz,X,Y,Z,SYM1,SYM2,SYM3,symmetry,sst)
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry,sst
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fxz
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8,dimension(-1:ex(1)+2,-1:ex(2)+2,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdxdz,Fdxdz
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -1
     if(dabs(Y(1)) < dY) jmin = -1
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -1
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+2
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(2,ex,f,fh,SoA)

  Sdxdz = F1o4 /( dX * dZ )

  Fdxdz = F1o144 /( dX * dZ )

  fxz = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
!~~~~~~ fxz
       if(i+2 <= imax .and. i-2 >= imin .and. k+2 <= kmax .and. k-2 >= kmin)then
   fxz(i,j,k) = Fdxdz*(     (fh(i-2,j,k-2)-F8*fh(i-1,j,k-2)+F8*fh(i+1,j,k-2)-fh(i+2,j,k-2))  &
                       -F8 *(fh(i-2,j,k-1)-F8*fh(i-1,j,k-1)+F8*fh(i+1,j,k-1)-fh(i+2,j,k-1))  &
                       +F8 *(fh(i-2,j,k+1)-F8*fh(i-1,j,k+1)+F8*fh(i+1,j,k+1)-fh(i+2,j,k+1))  &
                       -    (fh(i-2,j,k+2)-F8*fh(i-1,j,k+2)+F8*fh(i+1,j,k+2)-fh(i+2,j,k+2)))
   elseif(i+1 <= imax .and. i-1 >= imin .and. k+1 <= kmax .and. k-1 >= kmin)then
   fxz(i,j,k) = Sdxdz*(fh(i-1,j,k-1)-fh(i+1,j,k-1)-fh(i-1,j,k+1)+fh(i+1,j,k+1))
   endif

   enddo
   enddo
   enddo

  return

  end subroutine fddxz_sh

  subroutine fddyz_sh(ex,f,fyz,X,Y,Z,SYM1,SYM2,SYM3,symmetry,sst)
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry,sst
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fyz
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8,dimension(-1:ex(1)+2,-1:ex(2)+2,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdydz,Fdydz
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -1
     if(dabs(Y(1)) < dY) jmin = -1
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -1
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+2
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(2,ex,f,fh,SoA)

  Sdydz = F1o4 /( dY * dZ )

  Fdydz = F1o144 /( dY * dZ )

  fyz = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
!~~~~~~ fyz
       if(j+2 <= jmax .and. j-2 >= jmin .and. k+2 <= kmax .and. k-2 >= kmin)then
   fyz(i,j,k) = Fdydz*(     (fh(i,j-2,k-2)-F8*fh(i,j-1,k-2)+F8*fh(i,j+1,k-2)-fh(i,j+2,k-2))  &
                       -F8 *(fh(i,j-2,k-1)-F8*fh(i,j-1,k-1)+F8*fh(i,j+1,k-1)-fh(i,j+2,k-1))  &
                       +F8 *(fh(i,j-2,k+1)-F8*fh(i,j-1,k+1)+F8*fh(i,j+1,k+1)-fh(i,j+2,k+1))  &
                       -    (fh(i,j-2,k+2)-F8*fh(i,j-1,k+2)+F8*fh(i,j+1,k+2)-fh(i,j+2,k+2)))
   elseif(j+1 <= jmax .and. j-1 >= jmin .and. k+1 <= kmax .and. k-1 >= kmin)then
   fyz(i,j,k) = Sdydz*(fh(i,j-1,k-1)-fh(i,j+1,k-1)-fh(i,j-1,k+1)+fh(i,j+1,k+1))
   endif 

   enddo
   enddo
   enddo

  return

  end subroutine fddyz_sh

#elif (ghost_width == 4)
! sixth order code

!-----------------------------------------------------------------------------------------------------------------
!
! General first derivatives of 6_th oder accurate
!
!           - f(i-3) + 9 f(i-2) - 45 f(i-1) + 45 f(i+1) - 9 f(i+2) + f(i+3)
!  fx(i) = -----------------------------------------------------------------
!                                        60 dx
!
!-----------------------------------------------------------------------------------------------------------------

  subroutine fderivs_sh(ex,f,fx,fy,fz,X,Y,Z,SYM1,SYM2,SYM3,symmetry,onoff,sst)
  implicit none

  integer,                               intent(in ):: ex(1:3),symmetry,onoff,sst
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in ):: f
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out):: fx,fy,fz
  real*8,                                intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8,                                intent(in ):: SYM1,SYM2,SYM3

!~~~~~~ other variables

  real*8 :: dX,dY,dZ
  real*8,dimension(-2:ex(1)+3,-2:ex(2)+3,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8 :: d60dx,d60dy,d60dz,d12dx,d12dy,d12dz,d2dx,d2dy,d2dz
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8,  parameter :: ZEO=0.d0,ONE=1.d0, F60=6.d1
  real*8,  parameter :: TWO=2.d0,EIT=8.d0
  real*8,  parameter ::  F9=9.d0,F45=4.5d1,F12=1.2d1

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -2
     if(dabs(Y(1)) < dY) jmin = -2
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -2
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+3
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(3,ex,f,fh,SoA)

  d60dx = ONE/F60/dX
  d60dy = ONE/F60/dY
  d60dz = ONE/F60/dZ

  d12dx = ONE/F12/dX
  d12dy = ONE/F12/dY
  d12dz = ONE/F12/dZ

  d2dx = ONE/TWO/dX
  d2dy = ONE/TWO/dY
  d2dz = ONE/TWO/dZ

  fx = ZEO
  fy = ZEO
  fz = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
! x direction   
        if(i+3 <= imax .and. i-3 >= imin)then
!                
!           - f(i-3) + 9 f(i-2) - 45 f(i-1) + 45 f(i+1) - 9 f(i+2) + f(i+3)
!  fx(i) = -----------------------------------------------------------------
!                                        60 dx
      fx(i,j,k)=d60dx*(-fh(i-3,j,k)+F9*fh(i-2,j,k)-F45*fh(i-1,j,k)+F45*fh(i+1,j,k)-F9*fh(i+2,j,k)+fh(i+3,j,k))

    elseif(i+2 <= imax .and. i-2 >= imin)then
!
!              f(i-2) - 8 f(i-1) + 8 f(i+1) - f(i+2)
!  fx(i) = ---------------------------------------------
!                             12 dx
      fx(i,j,k)=d12dx*(fh(i-2,j,k)-EIT*fh(i-1,j,k)+EIT*fh(i+1,j,k)-fh(i+2,j,k))

    elseif(i+1 <= imax .and. i-1 >= imin)then
!
!              - f(i-1) + f(i+1)
!  fx(i) = --------------------------------
!                     2 dx
      fx(i,j,k)=d2dx*(-fh(i-1,j,k)+fh(i+1,j,k))

! set imax and imin 0
    endif
! y direction   
        if(j+3 <= jmax .and. j-3 >= jmin)then

      fy(i,j,k)=d60dy*(-fh(i,j-3,k)+F9*fh(i,j-2,k)-F45*fh(i,j-1,k)+F45*fh(i,j+1,k)-F9*fh(i,j+2,k)+fh(i,j+3,k))

    elseif(j+2 <= jmax .and. j-2 >= jmin)then

      fy(i,j,k)=d12dy*(fh(i,j-2,k)-EIT*fh(i,j-1,k)+EIT*fh(i,j+1,k)-fh(i,j+2,k))

    elseif(j+1 <= jmax .and. j-1 >= jmin)then

     fy(i,j,k)=d2dy*(-fh(i,j-1,k)+fh(i,j+1,k))

! set jmax and jmin 0
    endif
! z direction   
        if(k+3 <= kmax .and. k-3 >= kmin)then

      fz(i,j,k)=d60dz*(-fh(i,j,k-3)+F9*fh(i,j,k-2)-F45*fh(i,j,k-1)+F45*fh(i,j,k+1)-F9*fh(i,j,k+2)+fh(i,j,k+3))

    elseif(k+2 <= kmax .and. k-2 >= kmin)then

      fz(i,j,k)=d12dz*(fh(i,j,k-2)-EIT*fh(i,j,k-1)+EIT*fh(i,j,k+1)-fh(i,j,k+2))

    elseif(k+1 <= kmax .and. k-1 >= kmin)then

      fz(i,j,k)=d2dz*(-fh(i,j,k-1)+fh(i,j,k+1))

! set kmax and kmin 0
    endif

  enddo
  enddo
  enddo

  return

  end subroutine fderivs_sh
!-----------------------------------------------------------------------------
!
! single derivatives dx
!
!-----------------------------------------------------------------------------
  subroutine fdx_sh(ex,f,fx,X,Y,Z,SYM1,SYM2,SYM3,symmetry,onoff,sst)
  implicit none

  integer,                               intent(in ):: ex(1:3),symmetry,onoff,sst
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in ):: f
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out):: fx
  real*8,                                intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8,                                intent(in ):: SYM1,SYM2,SYM3

!~~~~~~ other variables

  real*8 :: dX,dY,dZ
  real*8,dimension(-2:ex(1)+3,-2:ex(2)+3,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8 :: d60dx,d12dx,d2dx
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8,  parameter :: ZEO=0.d0,ONE=1.d0, F60=6.d1
  real*8,  parameter :: TWO=2.d0,EIT=8.d0
  real*8,  parameter ::  F9=9.d0,F45=4.5d1,F12=1.2d1

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)
  
  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -2
     if(dabs(Y(1)) < dY) jmin = -2
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -2
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+3
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(3,ex,f,fh,SoA)

  d60dx = ONE/F60/dX

  d12dx = ONE/F12/dX

  d2dx = ONE/TWO/dX

  fx = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
! x direction   
        if(i+3 <= imax .and. i-3 >= imin)then
!                
!           - f(i-3) + 9 f(i-2) - 45 f(i-1) + 45 f(i+1) - 9 f(i+2) + f(i+3)
!  fx(i) = -----------------------------------------------------------------
!                                        60 dx
      fx(i,j,k)=d60dx*(-fh(i-3,j,k)+F9*fh(i-2,j,k)-F45*fh(i-1,j,k)+F45*fh(i+1,j,k)-F9*fh(i+2,j,k)+fh(i+3,j,k))

    elseif(i+2 <= imax .and. i-2 >= imin)then
!
!              f(i-2) - 8 f(i-1) + 8 f(i+1) - f(i+2)
!  fx(i) = ---------------------------------------------
!                             12 dx
      fx(i,j,k)=d12dx*(fh(i-2,j,k)-EIT*fh(i-1,j,k)+EIT*fh(i+1,j,k)-fh(i+2,j,k))

    elseif(i+1 <= imax .and. i-1 >= imin)then
!
!              - f(i-1) + f(i+1)
!  fx(i) = --------------------------------
!                     2 dx
      fx(i,j,k)=d2dx*(-fh(i-1,j,k)+fh(i+1,j,k))

! set imax and imin 0
    endif

  enddo
  enddo
  enddo

  return

  end subroutine fdx_sh
!-----------------------------------------------------------------------------
!
! single derivatives dy
!
!-----------------------------------------------------------------------------
  subroutine fdy_sh(ex,f,fy,X,Y,Z,SYM1,SYM2,SYM3,symmetry,onoff,sst)
  implicit none

  integer,                               intent(in ):: ex(1:3),symmetry,onoff,sst
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in ):: f
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out):: fy
  real*8,                                intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8,                                intent(in ):: SYM1,SYM2,SYM3

!~~~~~~ other variables

  real*8 :: dX,dY,dZ
  real*8,dimension(-2:ex(1)+3,-2:ex(2)+3,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8 :: d60dy,d12dy,d2dy
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8,  parameter :: ZEO=0.d0,ONE=1.d0, F60=6.d1
  real*8,  parameter :: TWO=2.d0,EIT=8.d0
  real*8,  parameter ::  F9=9.d0,F45=4.5d1,F12=1.2d1

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)
  
  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -2
     if(dabs(Y(1)) < dY) jmin = -2
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -2
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+3
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(3,ex,f,fh,SoA)

  d60dy = ONE/F60/dY

  d12dy = ONE/F12/dY

  d2dy = ONE/TWO/dY

  fy = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
! y direction   
        if(j+3 <= jmax .and. j-3 >= jmin)then

      fy(i,j,k)=d60dy*(-fh(i,j-3,k)+F9*fh(i,j-2,k)-F45*fh(i,j-1,k)+F45*fh(i,j+1,k)-F9*fh(i,j+2,k)+fh(i,j+3,k))

    elseif(j+2 <= jmax .and. j-2 >= jmin)then

      fy(i,j,k)=d12dy*(fh(i,j-2,k)-EIT*fh(i,j-1,k)+EIT*fh(i,j+1,k)-fh(i,j+2,k))

    elseif(j+1 <= jmax .and. j-1 >= jmin)then

     fy(i,j,k)=d2dy*(-fh(i,j-1,k)+fh(i,j+1,k))

! set jmax and jmin 0
    endif

  enddo
  enddo
  enddo

  return

  end subroutine fdy_sh
!-----------------------------------------------------------------------------
!
! single derivatives dz
!
!-----------------------------------------------------------------------------
  subroutine fdz_sh(ex,f,fz,X,Y,Z,SYM1,SYM2,SYM3,symmetry,onoff,sst)
  implicit none

  integer,                               intent(in ):: ex(1:3),symmetry,onoff,sst
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in ):: f
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out):: fz
  real*8,                                intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8,                                intent(in ):: SYM1,SYM2,SYM3

!~~~~~~ other variables
  
  real*8 :: dX,dY,dZ
  real*8,dimension(-2:ex(1)+3,-2:ex(2)+3,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8 :: d60dz,d12dz,d2dz
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8,  parameter :: ZEO=0.d0,ONE=1.d0, F60=6.d1
  real*8,  parameter :: TWO=2.d0,EIT=8.d0
  real*8,  parameter ::  F9=9.d0,F45=4.5d1,F12=1.2d1

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)
  
  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -2
     if(dabs(Y(1)) < dY) jmin = -2
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -2
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+3
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(3,ex,f,fh,SoA)

  d60dz = ONE/F60/dZ

  d12dz = ONE/F12/dZ

  d2dz = ONE/TWO/dZ

  fz = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
! z direction   
        if(k+3 <= kmax .and. k-3 >= kmin)then

      fz(i,j,k)=d60dz*(-fh(i,j,k-3)+F9*fh(i,j,k-2)-F45*fh(i,j,k-1)+F45*fh(i,j,k+1)-F9*fh(i,j,k+2)+fh(i,j,k+3))

    elseif(k+2 <= kmax .and. k-2 >= kmin)then

      fz(i,j,k)=d12dz*(fh(i,j,k-2)-EIT*fh(i,j,k-1)+EIT*fh(i,j,k+1)-fh(i,j,k+2))

    elseif(k+1 <= kmax .and. k-1 >= kmin)then

      fz(i,j,k)=d2dz*(-fh(i,j,k-1)+fh(i,j,k+1))

! set kmax and kmin 0
    endif

  enddo
  enddo
  enddo

  return

  end subroutine fdz_sh
!-----------------------------------------------------------------------------------------------------------------
!
! General second derivatives of 6_th oder accurate
!
!             2 f(i-3) - 27 f(i-2) + 270 f(i-1) - 490 f(i) + 270 f(i+1) - 27 f(i+2) + 2 f(i+3)
!  fxx(i) = -----------------------------------------------------------------------------------
!                                                180 dx^2 
!
!             -    ( - f(i-3,j-3) + 9 f(i-2,j-3) - 45 f(i-1,j-3) + 45 f(i+1,j-3) - 9 f(i+2,j-3) + f(i+3,j-3) )
!             + 9  ( - f(i-3,j-2) + 9 f(i-2,j-2) - 45 f(i-1,j-2) + 45 f(i+1,j-2) - 9 f(i+2,j-2) + f(i+3,j-2) )
!             - 45 ( - f(i-3,j-1) + 9 f(i-2,j-1) - 45 f(i-1,j-1) + 45 f(i+1,j-1) - 9 f(i+2,j-1) + f(i+3,j-1) )
!             + 45 ( - f(i-3,j+1) + 9 f(i-2,j+1) - 45 f(i-1,j+1) + 45 f(i+1,j+1) - 9 f(i+2,j+1) + f(i+3,j+1) )
!             - 9  ( - f(i-3,j+2) + 9 f(i-2,j+2) - 45 f(i-1,j+2) + 45 f(i+1,j+2) - 9 f(i+2,j+2) + f(i+3,j+2) )
!             +    ( - f(i-3,j+3) + 9 f(i-2,j+3) - 45 f(i-1,j+3) + 45 f(i+1,j+3) - 9 f(i+2,j+3) + f(i+3,j+3) )  
!  fxy(i,j) = ------------------------------------------------------------------------------------------------
!                                                          3600 dx dy
!
!-----------------------------------------------------------------------------------------------------------------
  subroutine fdderivs_sh(ex,f,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z, &
                      SYM1,SYM2,SYM3,symmetry,onoff,sst)
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry,onoff,sst
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fxx,fxy,fxz,fyy,fyz,fzz
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8,dimension(-2:ex(1)+3,-2:ex(2)+3,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdxdx,Sdydy,Sdzdz,Fdxdx,Fdydy,Fdzdz,Xdxdx,Xdydy,Xdzdz
  real*8  :: Sdxdy,Sdxdz,Sdydz,Fdxdy,Fdxdz,Fdydz,Xdxdy,Xdxdz,Xdydz
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -2
     if(dabs(Y(1)) < dY) jmin = -2
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -2
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+3
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(3,ex,f,fh,SoA)

  Sdxdx =  ONE /( dX * dX )
  Sdydy =  ONE /( dY * dY )
  Sdzdz =  ONE /( dZ * dZ )

  Fdxdx = F1o12 /( dX * dX )
  Fdydy = F1o12 /( dY * dY )
  Fdzdz = F1o12 /( dZ * dZ )

  Xdxdx = F1o180 /( dX * dX )
  Xdydy = F1o180 /( dY * dY )
  Xdzdz = F1o180 /( dZ * dZ )

  Sdxdy = F1o4 /( dX * dY )
  Sdxdz = F1o4 /( dX * dZ )
  Sdydz = F1o4 /( dY * dZ )

  Fdxdy = F1o144 /( dX * dY )
  Fdxdz = F1o144 /( dX * dZ )
  Fdydz = F1o144 /( dY * dZ )

  Xdxdy = F1o3600 /( dX * dY )
  Xdxdz = F1o3600 /( dX * dZ )
  Xdydz = F1o3600 /( dY * dZ )

  fxx = ZEO
  fyy = ZEO
  fzz = ZEO
  fxy = ZEO
  fxz = ZEO
  fyz = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
!~~~~~~ fxx
        if(i+3 <= imax .and. i-3 >= imin)then
!
!             2 f(i-3) - 27 f(i-2) + 270 f(i-1) - 490 f(i) + 270 f(i+1) - 27 f(i+2) + 2 f(i+3)
!  fxx(i) = -----------------------------------------------------------------------------------
!                                                180 dx^2 
   fxx(i,j,k) = Xdxdx*(TWO*fh(i-3,j,k)-F27*fh(i-2,j,k)+F270*fh(i-1,j,k)-F490*fh(i,j,k) &
                      +TWO*fh(i+3,j,k)-F27*fh(i+2,j,k)+F270*fh(i+1,j,k)               )
   elseif(i+2 <= imax .and. i-2 >= imin)then
!
!               - f(i-2) + 16 f(i-1) - 30 f(i) + 16 f(i+1) - f(i+2)
!  fxx(i) = ----------------------------------------------------------
!                                  12 dx^2 
   fxx(i,j,k) = Fdxdx*(-fh(i-2,j,k)+F16*fh(i-1,j,k)-F30*fh(i,j,k) &
                       -fh(i+2,j,k)+F16*fh(i+1,j,k)              )
   elseif(i+1 <= imax .and. i-1 >= imin)then
!
!               f(i-1) - 2 f(i) + f(i+1)
!  fxx(i) = --------------------------------
!                         dx^2 
   fxx(i,j,k) = Sdxdx*(fh(i-1,j,k)-TWO*fh(i,j,k) &
                      +fh(i+1,j,k)              )
   endif


!~~~~~~ fyy
        if(j+3 <= jmax .and. j-3 >= jmin)then

   fyy(i,j,k) = Xdydy*(TWO*fh(i,j-3,k)-F27*fh(i,j-2,k)+F270*fh(i,j-1,k)-F490*fh(i,j,k) &
                      +TWO*fh(i,j+3,k)-F27*fh(i,j+2,k)+F270*fh(i,j+1,k)               )
   elseif(j+2 <= jmax .and. j-2 >= jmin)then

   fyy(i,j,k) = Fdydy*(-fh(i,j-2,k)+F16*fh(i,j-1,k)-F30*fh(i,j,k) &
                       -fh(i,j+2,k)+F16*fh(i,j+1,k)              )
   elseif(j+1 <= jmax .and. j-1 >= jmin)then

   fyy(i,j,k) = Sdydy*(fh(i,j-1,k)-TWO*fh(i,j,k) &
                      +fh(i,j+1,k)              )
   endif

!~~~~~~ fzz
        if(k+3 <= kmax .and. k-3 >= kmin)then

   fzz(i,j,k) = Xdzdz*(TWO*fh(i,j,k-3)-F27*fh(i,j,k-2)+F270*fh(i,j,k-1)-F490*fh(i,j,k) &
                      +TWO*fh(i,j,k+3)-F27*fh(i,j,k+2)+F270*fh(i,j,k+1)               )
   elseif(k+2 <= kmax .and. k-2 >= kmin)then

   fzz(i,j,k) = Fdzdz*(-fh(i,j,k-2)+F16*fh(i,j,k-1)-F30*fh(i,j,k) &
                       -fh(i,j,k+2)+F16*fh(i,j,k+1)              )
   elseif(k+1 <= kmax .and. k-1 >= kmin)then

   fzz(i,j,k) = Sdzdz*(fh(i,j,k-1)-TWO*fh(i,j,k) &
                      +fh(i,j,k+1)              )
   endif
!~~~~~~ fxy
      if(i+3 <= imax .and. i-3 >= imin .and. j+3 <= jmax .and. j-3 >= jmin)then
!
!             -    ( - f(i-3,j-3) + 9 f(i-2,j-3) - 45 f(i-1,j-3) + 45 f(i+1,j-3) - 9 f(i+2,j-3) + f(i+3,j-3) )
!             + 9  ( - f(i-3,j-2) + 9 f(i-2,j-2) - 45 f(i-1,j-2) + 45 f(i+1,j-2) - 9 f(i+2,j-2) + f(i+3,j-2) )
!             - 45 ( - f(i-3,j-1) + 9 f(i-2,j-1) - 45 f(i-1,j-1) + 45 f(i+1,j-1) - 9 f(i+2,j-1) + f(i+3,j-1) )
!             + 45 ( - f(i-3,j+1) + 9 f(i-2,j+1) - 45 f(i-1,j+1) + 45 f(i+1,j+1) - 9 f(i+2,j+1) + f(i+3,j+1) )
!             - 9  ( - f(i-3,j+2) + 9 f(i-2,j+2) - 45 f(i-1,j+2) + 45 f(i+1,j+2) - 9 f(i+2,j+2) + f(i+3,j+2) )
!             +    ( - f(i-3,j+3) + 9 f(i-2,j+3) - 45 f(i-1,j+3) + 45 f(i+1,j+3) - 9 f(i+2,j+3) + f(i+3,j+3) )  
!  fxy(i,j) = ------------------------------------------------------------------------------------------------
!                                                          3600 dx dy
   fxy(i,j,k) = Xdxdy*(-    (-fh(i-3,j-3,k)+F9*fh(i-2,j-3,k)-F45*fh(i-1,j-3,k)+F45*fh(i+1,j-3,k)-F9*fh(i+2,j-3,k)+fh(i+3,j-3,k))  &
                       +F9 *(-fh(i-3,j-2,k)+F9*fh(i-2,j-2,k)-F45*fh(i-1,j-2,k)+F45*fh(i+1,j-2,k)-F9*fh(i+2,j-2,k)+fh(i+3,j-2,k))  &
                       -F45*(-fh(i-3,j-1,k)+F9*fh(i-2,j-1,k)-F45*fh(i-1,j-1,k)+F45*fh(i+1,j-1,k)-F9*fh(i+2,j-1,k)+fh(i+3,j-1,k))  &
                       +F45*(-fh(i-3,j+1,k)+F9*fh(i-2,j+1,k)-F45*fh(i-1,j+1,k)+F45*fh(i+1,j+1,k)-F9*fh(i+2,j+1,k)+fh(i+3,j+1,k))  &
                       -F9 *(-fh(i-3,j+2,k)+F9*fh(i-2,j+2,k)-F45*fh(i-1,j+2,k)+F45*fh(i+1,j+2,k)-F9*fh(i+2,j+2,k)+fh(i+3,j+2,k))  &
                       +    (-fh(i-3,j+3,k)+F9*fh(i-2,j+3,k)-F45*fh(i-1,j+3,k)+F45*fh(i+1,j+3,k)-F9*fh(i+2,j+3,k)+fh(i+3,j+3,k)))
   elseif(i+2 <= imax .and. i-2 >= imin .and. j+2 <= jmax .and. j-2 >= jmin)then
!
!                 ( f(i-2,j-2) - 8 f(i-1,j-2) + 8 f(i+1,j-2) - f(i+2,j-2) )
!             - 8 ( f(i-2,j-1) - 8 f(i-1,j-1) + 8 f(i+1,j-1) - f(i+2,j-1) )
!             + 8 ( f(i-2,j+1) - 8 f(i-1,j+1) + 8 f(i+1,j+1) - f(i+2,j+1) )
!             -   ( f(i-2,j+2) - 8 f(i-1,j+2) + 8 f(i+1,j+2) - f(i+2,j+2) )
!  fxy(i,j) = ----------------------------------------------------------------
!                                  144 dx dy
   fxy(i,j,k) = Fdxdy*(     (fh(i-2,j-2,k)-F8*fh(i-1,j-2,k)+F8*fh(i+1,j-2,k)-fh(i+2,j-2,k))  &
                       -F8 *(fh(i-2,j-1,k)-F8*fh(i-1,j-1,k)+F8*fh(i+1,j-1,k)-fh(i+2,j-1,k))  &
                       +F8 *(fh(i-2,j+1,k)-F8*fh(i-1,j+1,k)+F8*fh(i+1,j+1,k)-fh(i+2,j+1,k))  &
                       -    (fh(i-2,j+2,k)-F8*fh(i-1,j+2,k)+F8*fh(i+1,j+2,k)-fh(i+2,j+2,k)))

   elseif(i+1 <= imax .and. i-1 >= imin .and. j+1 <= jmax .and. j-1 >= jmin)then
!                 f(i-1,j-1) - f(i+1,j-1) - f(i-1,j+1) + f(i+1,j+1) 
!  fxy(i,j) = -----------------------------------------------------------
!                                      4 dx dy
   fxy(i,j,k) = Sdxdy*(fh(i-1,j-1,k)-fh(i+1,j-1,k)-fh(i-1,j+1,k)+fh(i+1,j+1,k))
   endif
!~~~~~~ fxz
      if(i+3 <= imax .and. i-3 >= imin .and. k+3 <= kmax .and. k-3 >= kmin)then

   fxz(i,j,k) = Xdxdz*(-    (-fh(i-3,j,k-3)+F9*fh(i-2,j,k-3)-F45*fh(i-1,j,k-3)+F45*fh(i+1,j,k-3)-F9*fh(i+2,j,k-3)+fh(i+3,j,k-3))  &
                       +F9 *(-fh(i-3,j,k-2)+F9*fh(i-2,j,k-2)-F45*fh(i-1,j,k-2)+F45*fh(i+1,j,k-2)-F9*fh(i+2,j,k-2)+fh(i+3,j,k-2))  &
                       -F45*(-fh(i-3,j,k-1)+F9*fh(i-2,j,k-1)-F45*fh(i-1,j,k-1)+F45*fh(i+1,j,k-1)-F9*fh(i+2,j,k-1)+fh(i+3,j,k-1))  &
                       +F45*(-fh(i-3,j,k+1)+F9*fh(i-2,j,k+1)-F45*fh(i-1,j,k+1)+F45*fh(i+1,j,k+1)-F9*fh(i+2,j,k+1)+fh(i+3,j,k+1))  &
                       -F9 *(-fh(i-3,j,k+2)+F9*fh(i-2,j,k+2)-F45*fh(i-1,j,k+2)+F45*fh(i+1,j,k+2)-F9*fh(i+2,j,k+2)+fh(i+3,j,k+2))  &
                       +    (-fh(i-3,j,k+3)+F9*fh(i-2,j,k+3)-F45*fh(i-1,j,k+3)+F45*fh(i+1,j,k+3)-F9*fh(i+2,j,k+3)+fh(i+3,j,k+3)))
   elseif(i+2 <= imax .and. i-2 >= imin .and. k+2 <= kmax .and. k-2 >= kmin)then
   fxz(i,j,k) = Fdxdz*(     (fh(i-2,j,k-2)-F8*fh(i-1,j,k-2)+F8*fh(i+1,j,k-2)-fh(i+2,j,k-2))  &
                       -F8 *(fh(i-2,j,k-1)-F8*fh(i-1,j,k-1)+F8*fh(i+1,j,k-1)-fh(i+2,j,k-1))  &
                       +F8 *(fh(i-2,j,k+1)-F8*fh(i-1,j,k+1)+F8*fh(i+1,j,k+1)-fh(i+2,j,k+1))  &
                       -    (fh(i-2,j,k+2)-F8*fh(i-1,j,k+2)+F8*fh(i+1,j,k+2)-fh(i+2,j,k+2)))
   elseif(i+1 <= imax .and. i-1 >= imin .and. k+1 <= kmax .and. k-1 >= kmin)then
   fxz(i,j,k) = Sdxdz*(fh(i-1,j,k-1)-fh(i+1,j,k-1)-fh(i-1,j,k+1)+fh(i+1,j,k+1))
   endif
!~~~~~~ fyz
      if(j+3 <= jmax .and. j-3 >= jmin .and. k+3 <= kmax .and. k-3 >= kmin)then

   fyz(i,j,k) = Xdydz*(-    (-fh(i,j-3,k-3)+F9*fh(i,j-2,k-3)-F45*fh(i,j-1,k-3)+F45*fh(i,j+1,k-3)-F9*fh(i,j+2,k-3)+fh(i,j+3,k-3))  &
                       +F9 *(-fh(i,j-3,k-2)+F9*fh(i,j-2,k-2)-F45*fh(i,j-1,k-2)+F45*fh(i,j+1,k-2)-F9*fh(i,j+2,k-2)+fh(i,j+3,k-2))  &
                       -F45*(-fh(i,j-3,k-1)+F9*fh(i,j-2,k-1)-F45*fh(i,j-1,k-1)+F45*fh(i,j+1,k-1)-F9*fh(i,j+2,k-1)+fh(i,j+3,k-1))  &
                       +F45*(-fh(i,j-3,k+1)+F9*fh(i,j-2,k+1)-F45*fh(i,j-1,k+1)+F45*fh(i,j+1,k+1)-F9*fh(i,j+2,k+1)+fh(i,j+3,k+1))  &
                       -F9 *(-fh(i,j-3,k+2)+F9*fh(i,j-2,k+2)-F45*fh(i,j-1,k+2)+F45*fh(i,j+1,k+2)-F9*fh(i,j+2,k+2)+fh(i,j+3,k+2))  &
                       +    (-fh(i,j-3,k+3)+F9*fh(i,j-2,k+3)-F45*fh(i,j-1,k+3)+F45*fh(i,j+1,k+3)-F9*fh(i,j+2,k+3)+fh(i,j+3,k+3)))
   elseif(j+2 <= jmax .and. j-2 >= jmin .and. k+2 <= kmax .and. k-2 >= kmin)then
   fyz(i,j,k) = Fdydz*(     (fh(i,j-2,k-2)-F8*fh(i,j-1,k-2)+F8*fh(i,j+1,k-2)-fh(i,j+2,k-2))  &
                       -F8 *(fh(i,j-2,k-1)-F8*fh(i,j-1,k-1)+F8*fh(i,j+1,k-1)-fh(i,j+2,k-1))  &
                       +F8 *(fh(i,j-2,k+1)-F8*fh(i,j-1,k+1)+F8*fh(i,j+1,k+1)-fh(i,j+2,k+1))  &
                       -    (fh(i,j-2,k+2)-F8*fh(i,j-1,k+2)+F8*fh(i,j+1,k+2)-fh(i,j+2,k+2)))
   elseif(j+1 <= jmax .and. j-1 >= jmin .and. k+1 <= kmax .and. k-1 >= kmin)then
   fyz(i,j,k) = Sdydz*(fh(i,j-1,k-1)-fh(i,j+1,k-1)-fh(i,j-1,k+1)+fh(i,j+1,k+1))
   endif 

   enddo
   enddo
   enddo

  return

  end subroutine fdderivs_sh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! only for compute_ricci.f90 usage
!-----------------------------------------------------------------------------
  subroutine fddxx_sh(ex,f,fxx,X,Y,Z,SYM1,SYM2,SYM3,symmetry,sst)
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry,sst
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fxx
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8,dimension(-2:ex(1)+3,-2:ex(2)+3,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdxdx,Fdxdx,Xdxdx
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -2
     if(dabs(Y(1)) < dY) jmin = -2
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -2
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+3
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(3,ex,f,fh,SoA)

  Sdxdx =  ONE /( dX * dX )

  Fdxdx = F1o12 /( dX * dX )

  Xdxdx = F1o180 /( dX * dX )

  fxx = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
!~~~~~~ fxx
        if(i+3 <= imax .and. i-3 >= imin)then
   fxx(i,j,k) = Xdxdx*(TWO*fh(i-3,j,k)-F27*fh(i-2,j,k)+F270*fh(i-1,j,k)-F490*fh(i,j,k) &
                      +TWO*fh(i+3,j,k)-F27*fh(i+2,j,k)+F270*fh(i+1,j,k)               )
   elseif(i+2 <= imax .and. i-2 >= imin)then
   fxx(i,j,k) = Fdxdx*(-fh(i-2,j,k)+F16*fh(i-1,j,k)-F30*fh(i,j,k) &
                       -fh(i+2,j,k)+F16*fh(i+1,j,k)              )
   elseif(i+1 <= imax .and. i-1 >= imin)then
   fxx(i,j,k) = Sdxdx*(fh(i-1,j,k)-TWO*fh(i,j,k) &
                      +fh(i+1,j,k)              )
   endif

   enddo
   enddo
   enddo

  return

  end subroutine fddxx_sh

  subroutine fddyy_sh(ex,f,fyy,X,Y,Z,SYM1,SYM2,SYM3,symmetry,sst)
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry,sst
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fyy
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8,dimension(-2:ex(1)+3,-2:ex(2)+3,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdydy,Fdydy,Xdydy
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -2
     if(dabs(Y(1)) < dY) jmin = -2
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -2
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+3
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(3,ex,f,fh,SoA)

  Sdydy =  ONE /( dY * dY )

  Fdydy = F1o12 /( dY * dY )

  Xdydy = F1o180 /( dY * dY )

  fyy = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
!~~~~~~ fyy
        if(j+3 <= jmax .and. j-3 >= jmin)then

   fyy(i,j,k) = Xdydy*(TWO*fh(i,j-3,k)-F27*fh(i,j-2,k)+F270*fh(i,j-1,k)-F490*fh(i,j,k) &
                      +TWO*fh(i,j+3,k)-F27*fh(i,j+2,k)+F270*fh(i,j+1,k)               )
   elseif(j+2 <= jmax .and. j-2 >= jmin)then

   fyy(i,j,k) = Fdydy*(-fh(i,j-2,k)+F16*fh(i,j-1,k)-F30*fh(i,j,k) &
                       -fh(i,j+2,k)+F16*fh(i,j+1,k)              )
   elseif(j+1 <= jmax .and. j-1 >= jmin)then

   fyy(i,j,k) = Sdydy*(fh(i,j-1,k)-TWO*fh(i,j,k) &
                      +fh(i,j+1,k)              )
   endif

   enddo
   enddo
   enddo

  return

  end subroutine fddyy_sh

  subroutine fddzz_sh(ex,f,fzz,X,Y,Z,SYM1,SYM2,SYM3,symmetry,sst)
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry,sst
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fzz
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8,dimension(-2:ex(1)+3,-2:ex(2)+3,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdzdz,Fdzdz,Xdzdz
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -2
     if(dabs(Y(1)) < dY) jmin = -2
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -2
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+3
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(3,ex,f,fh,SoA)

  Sdzdz =  ONE /( dZ * dZ )

  Fdzdz = F1o12 /( dZ * dZ )

  Xdzdz = F1o180 /( dZ * dZ )

  fzz = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
!~~~~~~ fzz
        if(k+3 <= kmax .and. k-3 >= kmin)then

   fzz(i,j,k) = Xdzdz*(TWO*fh(i,j,k-3)-F27*fh(i,j,k-2)+F270*fh(i,j,k-1)-F490*fh(i,j,k) &
                      +TWO*fh(i,j,k+3)-F27*fh(i,j,k+2)+F270*fh(i,j,k+1)               )
   elseif(k+2 <= kmax .and. k-2 >= kmin)then

   fzz(i,j,k) = Fdzdz*(-fh(i,j,k-2)+F16*fh(i,j,k-1)-F30*fh(i,j,k) &
                       -fh(i,j,k+2)+F16*fh(i,j,k+1)              )
   elseif(k+1 <= kmax .and. k-1 >= kmin)then

   fzz(i,j,k) = Sdzdz*(fh(i,j,k-1)-TWO*fh(i,j,k) &
                      +fh(i,j,k+1)              )
   endif

   enddo
   enddo
   enddo

  return

  end subroutine fddzz_sh

  subroutine fddxy_sh(ex,f,fxy,X,Y,Z,SYM1,SYM2,SYM3,symmetry,sst)
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry,sst
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fxy
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8,dimension(-2:ex(1)+3,-2:ex(2)+3,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdxdy,Fdxdy,Xdxdy
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -2
     if(dabs(Y(1)) < dY) jmin = -2
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -2
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+3
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(3,ex,f,fh,SoA)

  Sdxdy = F1o4 /( dX * dY )

  Fdxdy = F1o144 /( dX * dY )

  Xdxdy = F1o3600 /( dX * dY )

  fxy = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
!~~~~~~ fxy
      if(i+3 <= imax .and. i-3 >= imin .and. j+3 <= jmax .and. j-3 >= jmin)then

   fxy(i,j,k) = Xdxdy*(-    (-fh(i-3,j-3,k)+F9*fh(i-2,j-3,k)-F45*fh(i-1,j-3,k)+F45*fh(i+1,j-3,k)-F9*fh(i+2,j-3,k)+fh(i+3,j-3,k))  &
                       +F9 *(-fh(i-3,j-2,k)+F9*fh(i-2,j-2,k)-F45*fh(i-1,j-2,k)+F45*fh(i+1,j-2,k)-F9*fh(i+2,j-2,k)+fh(i+3,j-2,k))  &
                       -F45*(-fh(i-3,j-1,k)+F9*fh(i-2,j-1,k)-F45*fh(i-1,j-1,k)+F45*fh(i+1,j-1,k)-F9*fh(i+2,j-1,k)+fh(i+3,j-1,k))  &
                       +F45*(-fh(i-3,j+1,k)+F9*fh(i-2,j+1,k)-F45*fh(i-1,j+1,k)+F45*fh(i+1,j+1,k)-F9*fh(i+2,j+1,k)+fh(i+3,j+1,k))  &
                       -F9 *(-fh(i-3,j+2,k)+F9*fh(i-2,j+2,k)-F45*fh(i-1,j+2,k)+F45*fh(i+1,j+2,k)-F9*fh(i+2,j+2,k)+fh(i+3,j+2,k))  &
                       +    (-fh(i-3,j+3,k)+F9*fh(i-2,j+3,k)-F45*fh(i-1,j+3,k)+F45*fh(i+1,j+3,k)-F9*fh(i+2,j+3,k)+fh(i+3,j+3,k)))
   elseif(i+2 <= imax .and. i-2 >= imin .and. j+2 <= jmax .and. j-2 >= jmin)then

   fxy(i,j,k) = Fdxdy*(     (fh(i-2,j-2,k)-F8*fh(i-1,j-2,k)+F8*fh(i+1,j-2,k)-fh(i+2,j-2,k))  &
                       -F8 *(fh(i-2,j-1,k)-F8*fh(i-1,j-1,k)+F8*fh(i+1,j-1,k)-fh(i+2,j-1,k))  &
                       +F8 *(fh(i-2,j+1,k)-F8*fh(i-1,j+1,k)+F8*fh(i+1,j+1,k)-fh(i+2,j+1,k))  &
                       -    (fh(i-2,j+2,k)-F8*fh(i-1,j+2,k)+F8*fh(i+1,j+2,k)-fh(i+2,j+2,k)))
   elseif(i+1 <= imax .and. i-1 >= imin .and. j+1 <= jmax .and. j-1 >= jmin)then

   fxy(i,j,k) = Sdxdy*(fh(i-1,j-1,k)-fh(i+1,j-1,k)-fh(i-1,j+1,k)+fh(i+1,j+1,k))
   endif

   enddo
   enddo
   enddo

  return

  end subroutine fddxy_sh

  subroutine fddxz_sh(ex,f,fxz,X,Y,Z,SYM1,SYM2,SYM3,symmetry,sst)
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry,sst
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fxz
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8,dimension(-2:ex(1)+3,-2:ex(2)+3,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdxdz,Fdxdz,Xdxdz
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -2
     if(dabs(Y(1)) < dY) jmin = -2
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -2
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+3
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(3,ex,f,fh,SoA)

  Sdxdz = F1o4 /( dX * dZ )

  Fdxdz = F1o144 /( dX * dZ )

  Xdxdz = F1o3600 /( dX * dZ )

  fxz = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
!~~~~~~ fxz
      if(i+3 <= imax .and. i-3 >= imin .and. k+3 <= kmax .and. k-3 >= kmin)then

   fxz(i,j,k) = Xdxdz*(-    (-fh(i-3,j,k-3)+F9*fh(i-2,j,k-3)-F45*fh(i-1,j,k-3)+F45*fh(i+1,j,k-3)-F9*fh(i+2,j,k-3)+fh(i+3,j,k-3))  &
                       +F9 *(-fh(i-3,j,k-2)+F9*fh(i-2,j,k-2)-F45*fh(i-1,j,k-2)+F45*fh(i+1,j,k-2)-F9*fh(i+2,j,k-2)+fh(i+3,j,k-2))  &
                       -F45*(-fh(i-3,j,k-1)+F9*fh(i-2,j,k-1)-F45*fh(i-1,j,k-1)+F45*fh(i+1,j,k-1)-F9*fh(i+2,j,k-1)+fh(i+3,j,k-1))  &
                       +F45*(-fh(i-3,j,k+1)+F9*fh(i-2,j,k+1)-F45*fh(i-1,j,k+1)+F45*fh(i+1,j,k+1)-F9*fh(i+2,j,k+1)+fh(i+3,j,k+1))  &
                       -F9 *(-fh(i-3,j,k+2)+F9*fh(i-2,j,k+2)-F45*fh(i-1,j,k+2)+F45*fh(i+1,j,k+2)-F9*fh(i+2,j,k+2)+fh(i+3,j,k+2))  &
                       +    (-fh(i-3,j,k+3)+F9*fh(i-2,j,k+3)-F45*fh(i-1,j,k+3)+F45*fh(i+1,j,k+3)-F9*fh(i+2,j,k+3)+fh(i+3,j,k+3)))
   elseif(i+2 <= imax .and. i-2 >= imin .and. k+2 <= kmax .and. k-2 >= kmin)then
   fxz(i,j,k) = Fdxdz*(     (fh(i-2,j,k-2)-F8*fh(i-1,j,k-2)+F8*fh(i+1,j,k-2)-fh(i+2,j,k-2))  &
                       -F8 *(fh(i-2,j,k-1)-F8*fh(i-1,j,k-1)+F8*fh(i+1,j,k-1)-fh(i+2,j,k-1))  &
                       +F8 *(fh(i-2,j,k+1)-F8*fh(i-1,j,k+1)+F8*fh(i+1,j,k+1)-fh(i+2,j,k+1))  &
                       -    (fh(i-2,j,k+2)-F8*fh(i-1,j,k+2)+F8*fh(i+1,j,k+2)-fh(i+2,j,k+2)))
   elseif(i+1 <= imax .and. i-1 >= imin .and. k+1 <= kmax .and. k-1 >= kmin)then
   fxz(i,j,k) = Sdxdz*(fh(i-1,j,k-1)-fh(i+1,j,k-1)-fh(i-1,j,k+1)+fh(i+1,j,k+1))
   endif

   enddo
   enddo
   enddo

  return

  end subroutine fddxz_sh

  subroutine fddyz_sh(ex,f,fyz,X,Y,Z,SYM1,SYM2,SYM3,symmetry,sst)
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry,sst
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fyz
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8,dimension(-2:ex(1)+3,-2:ex(2)+3,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdydz,Fdydz,Xdydz
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -2
     if(dabs(Y(1)) < dY) jmin = -2
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -2
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+3
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(3,ex,f,fh,SoA)

  Sdydz = F1o4 /( dY * dZ )

  Fdydz = F1o144 /( dY * dZ )

  Xdydz = F1o3600 /( dY * dZ )

  fyz = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
!~~~~~~ fyz
      if(j+3 <= jmax .and. j-3 >= jmin .and. k+3 <= kmax .and. k-3 >= kmin)then

   fyz(i,j,k) = Xdydz*(-    (-fh(i,j-3,k-3)+F9*fh(i,j-2,k-3)-F45*fh(i,j-1,k-3)+F45*fh(i,j+1,k-3)-F9*fh(i,j+2,k-3)+fh(i,j+3,k-3))  &
                       +F9 *(-fh(i,j-3,k-2)+F9*fh(i,j-2,k-2)-F45*fh(i,j-1,k-2)+F45*fh(i,j+1,k-2)-F9*fh(i,j+2,k-2)+fh(i,j+3,k-2))  &
                       -F45*(-fh(i,j-3,k-1)+F9*fh(i,j-2,k-1)-F45*fh(i,j-1,k-1)+F45*fh(i,j+1,k-1)-F9*fh(i,j+2,k-1)+fh(i,j+3,k-1))  &
                       +F45*(-fh(i,j-3,k+1)+F9*fh(i,j-2,k+1)-F45*fh(i,j-1,k+1)+F45*fh(i,j+1,k+1)-F9*fh(i,j+2,k+1)+fh(i,j+3,k+1))  &
                       -F9 *(-fh(i,j-3,k+2)+F9*fh(i,j-2,k+2)-F45*fh(i,j-1,k+2)+F45*fh(i,j+1,k+2)-F9*fh(i,j+2,k+2)+fh(i,j+3,k+2))  &
                       +    (-fh(i,j-3,k+3)+F9*fh(i,j-2,k+3)-F45*fh(i,j-1,k+3)+F45*fh(i,j+1,k+3)-F9*fh(i,j+2,k+3)+fh(i,j+3,k+3)))
   elseif(j+2 <= jmax .and. j-2 >= jmin .and. k+2 <= kmax .and. k-2 >= kmin)then
   fyz(i,j,k) = Fdydz*(     (fh(i,j-2,k-2)-F8*fh(i,j-1,k-2)+F8*fh(i,j+1,k-2)-fh(i,j+2,k-2))  &
                       -F8 *(fh(i,j-2,k-1)-F8*fh(i,j-1,k-1)+F8*fh(i,j+1,k-1)-fh(i,j+2,k-1))  &
                       +F8 *(fh(i,j-2,k+1)-F8*fh(i,j-1,k+1)+F8*fh(i,j+1,k+1)-fh(i,j+2,k+1))  &
                       -    (fh(i,j-2,k+2)-F8*fh(i,j-1,k+2)+F8*fh(i,j+1,k+2)-fh(i,j+2,k+2)))
   elseif(j+1 <= jmax .and. j-1 >= jmin .and. k+1 <= kmax .and. k-1 >= kmin)then
   fyz(i,j,k) = Sdydz*(fh(i,j-1,k-1)-fh(i,j+1,k-1)-fh(i,j-1,k+1)+fh(i,j+1,k+1))
   endif 

   enddo
   enddo
   enddo

  return

  end subroutine fddyz_sh

#elif (ghost_width == 5)
! eighth order code

! PRD 77, 024034 (2008)
!-----------------------------------------------------------------------------------------------------------------
!
! General first derivatives of 8_th oder accurate
!
!           3 f(i-4) - 32 f(i-3) + 168 f(i-2) - 672 f(i-1) + 672 f(i+1) - 168 f(i+2) + 32 f(i+3) - 3 f(i+4)
!  fx(i) = -------------------------------------------------------------------------------------------------
!                                                        840 dx
!
!-----------------------------------------------------------------------------------------------------------------

  subroutine fderivs_sh(ex,f,fx,fy,fz,X,Y,Z,SYM1,SYM2,SYM3,symmetry,onoff,sst)
  implicit none

  integer,                               intent(in ):: ex(1:3),symmetry,onoff,sst
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in ):: f
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out):: fx,fy,fz
  real*8,                                intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8,                                intent(in ):: SYM1,SYM2,SYM3

!~~~~~~ other variables

  real*8 :: dX,dY,dZ
  real*8,dimension(-3:ex(1)+4,-3:ex(2)+4,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8 :: d840dx,d840dy,d840dz
  real*8 :: d60dx,d60dy,d60dz,d12dx,d12dy,d12dz,d2dx,d2dy,d2dz
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8,  parameter :: ZEO=0.d0,ONE=1.d0, F60=6.d1, F32 = 3.2d1
  real*8,  parameter :: TWO=2.d0,THR=3.d0, EIT=8.d0, F168=1.68d2
  real*8,  parameter ::  F9=9.d0,F45=4.5d1,F12=1.2d1,F672=6.72d2
  real*8,  parameter :: F840=8.4d2

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -3
     if(dabs(Y(1)) < dY) jmin = -3
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -3
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+4
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(4,ex,f,fh,SoA)
  
  d840dx = ONE/F840/dX
  d840dy = ONE/F840/dY
  d840dz = ONE/F840/dZ

  d60dx = ONE/F60/dX
  d60dy = ONE/F60/dY
  d60dz = ONE/F60/dZ

  d12dx = ONE/F12/dX
  d12dy = ONE/F12/dY
  d12dz = ONE/F12/dZ

  d2dx = ONE/TWO/dX
  d2dy = ONE/TWO/dY
  d2dz = ONE/TWO/dZ

  fx = ZEO
  fy = ZEO
  fz = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
! x direction
        if(i+4 <= imax .and. i-4 >= imin)then
!           3 f(i-4) - 32 f(i-3) + 168 f(i-2) - 672 f(i-1) + 672 f(i+1) - 168 f(i+2) + 32 f(i+3) - 3 f(i+4)
!  fx(i) = -------------------------------------------------------------------------------------------------
!                                                        840 dx
      fx(i,j,k)=d840dx*( THR*fh(i-4,j,k)-F32 *fh(i-3,j,k)+F168*fh(i-2,j,k)-F672*fh(i-1,j,k)+ &
                        F672*fh(i+1,j,k)-F168*fh(i+2,j,k)+F32 *fh(i+3,j,k)-THR *fh(i+4,j,k))

    elseif(i+3 <= imax .and. i-3 >= imin)then
!                
!           - f(i-3) + 9 f(i-2) - 45 f(i-1) + 45 f(i+1) - 9 f(i+2) + f(i+3)
!  fx(i) = -----------------------------------------------------------------
!                                        60 dx
      fx(i,j,k)=d60dx*(-fh(i-3,j,k)+F9*fh(i-2,j,k)-F45*fh(i-1,j,k)+F45*fh(i+1,j,k)-F9*fh(i+2,j,k)+fh(i+3,j,k))

    elseif(i+2 <= imax .and. i-2 >= imin)then
!
!              f(i-2) - 8 f(i-1) + 8 f(i+1) - f(i+2)
!  fx(i) = ---------------------------------------------
!                             12 dx
      fx(i,j,k)=d12dx*(fh(i-2,j,k)-EIT*fh(i-1,j,k)+EIT*fh(i+1,j,k)-fh(i+2,j,k))

    elseif(i+1 <= imax .and. i-1 >= imin)then
!
!              - f(i-1) + f(i+1)
!  fx(i) = --------------------------------
!                     2 dx
      fx(i,j,k)=d2dx*(-fh(i-1,j,k)+fh(i+1,j,k))

! set imax and imin 0
    endif
! y direction   
        if(j+4 <= jmax .and. j-4 >= jmin)then

      fy(i,j,k)=d840dy*( THR*fh(i,j-4,k)-F32 *fh(i,j-3,k)+F168*fh(i,j-2,k)-F672*fh(i,j-1,k)+ &
                        F672*fh(i,j+1,k)-F168*fh(i,j+2,k)+F32 *fh(i,j+3,k)-THR *fh(i,j+4,k))

    elseif(j+3 <= jmax .and. j-3 >= jmin)then

      fy(i,j,k)=d60dy*(-fh(i,j-3,k)+F9*fh(i,j-2,k)-F45*fh(i,j-1,k)+F45*fh(i,j+1,k)-F9*fh(i,j+2,k)+fh(i,j+3,k))

    elseif(j+2 <= jmax .and. j-2 >= jmin)then

      fy(i,j,k)=d12dy*(fh(i,j-2,k)-EIT*fh(i,j-1,k)+EIT*fh(i,j+1,k)-fh(i,j+2,k))

    elseif(j+1 <= jmax .and. j-1 >= jmin)then

     fy(i,j,k)=d2dy*(-fh(i,j-1,k)+fh(i,j+1,k))

! set jmax and jmin 0
    endif
! z direction   
        if(k+4 <= kmax .and. k-4 >= kmin)then

      fz(i,j,k)=d840dz*( THR*fh(i,j,k-4)-F32 *fh(i,j,k-3)+F168*fh(i,j,k-2)-F672*fh(i,j,k-1)+ &
                        F672*fh(i,j,k+1)-F168*fh(i,j,k+2)+F32 *fh(i,j,k+3)-THR *fh(i,j,k+4))

    elseif(k+3 <= kmax .and. k-3 >= kmin)then

      fz(i,j,k)=d60dz*(-fh(i,j,k-3)+F9*fh(i,j,k-2)-F45*fh(i,j,k-1)+F45*fh(i,j,k+1)-F9*fh(i,j,k+2)+fh(i,j,k+3))

    elseif(k+2 <= kmax .and. k-2 >= kmin)then

      fz(i,j,k)=d12dz*(fh(i,j,k-2)-EIT*fh(i,j,k-1)+EIT*fh(i,j,k+1)-fh(i,j,k+2))

    elseif(k+1 <= kmax .and. k-1 >= kmin)then

      fz(i,j,k)=d2dz*(-fh(i,j,k-1)+fh(i,j,k+1))

! set kmax and kmin 0
    endif

  enddo
  enddo
  enddo

  return

  end subroutine fderivs_sh
!-----------------------------------------------------------------------------
!
! single derivatives dx
!
!-----------------------------------------------------------------------------
  subroutine fdx_sh(ex,f,fx,X,Y,Z,SYM1,SYM2,SYM3,symmetry,onoff,sst)
  implicit none

  integer,                               intent(in ):: ex(1:3),symmetry,onoff,sst
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in ):: f
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out):: fx
  real*8,                                intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8,                                intent(in ):: SYM1,SYM2,SYM3

!~~~~~~ other variables

  real*8 :: dX,dY,dZ
  real*8,dimension(-3:ex(1)+4,-3:ex(2)+4,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8 :: d840dx,d60dx,d12dx,d2dx
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8,  parameter :: ZEO=0.d0,ONE=1.d0, F60=6.d1, F32 = 3.2d1
  real*8,  parameter :: TWO=2.d0,THR=3.d0, EIT=8.d0, F168=1.68d2
  real*8,  parameter ::  F9=9.d0,F45=4.5d1,F12=1.2d1,F672=6.72d2
  real*8,  parameter :: F840=8.4d2

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)
  
  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -3
     if(dabs(Y(1)) < dY) jmin = -3
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -3
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+4
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(4,ex,f,fh,SoA)

  d840dx = ONE/F840/dX

  d60dx = ONE/F60/dX

  d12dx = ONE/F12/dX

  d2dx = ONE/TWO/dX

  fx = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
! x direction   
        if(i+4 <= imax .and. i-4 >= imin)then
!           3 f(i-4) - 32 f(i-3) + 168 f(i-2) - 672 f(i-1) + 672 f(i+1) - 168 f(i+2) + 32 f(i+3) - 3 f(i+4)
!  fx(i) = -------------------------------------------------------------------------------------------------
!                                                        840 dx
      fx(i,j,k)=d840dx*( THR*fh(i-4,j,k)-F32 *fh(i-3,j,k)+F168*fh(i-2,j,k)-F672*fh(i-1,j,k)+ &
                        F672*fh(i+1,j,k)-F168*fh(i+2,j,k)+F32 *fh(i+3,j,k)-THR *fh(i+4,j,k))

    elseif(i+3 <= imax .and. i-3 >= imin)then
!                
!           - f(i-3) + 9 f(i-2) - 45 f(i-1) + 45 f(i+1) - 9 f(i+2) + f(i+3)
!  fx(i) = -----------------------------------------------------------------
!                                        60 dx
      fx(i,j,k)=d60dx*(-fh(i-3,j,k)+F9*fh(i-2,j,k)-F45*fh(i-1,j,k)+F45*fh(i+1,j,k)-F9*fh(i+2,j,k)+fh(i+3,j,k))

    elseif(i+2 <= imax .and. i-2 >= imin)then
!
!              f(i-2) - 8 f(i-1) + 8 f(i+1) - f(i+2)
!  fx(i) = ---------------------------------------------
!                             12 dx
      fx(i,j,k)=d12dx*(fh(i-2,j,k)-EIT*fh(i-1,j,k)+EIT*fh(i+1,j,k)-fh(i+2,j,k))

    elseif(i+1 <= imax .and. i-1 >= imin)then
!
!              - f(i-1) + f(i+1)
!  fx(i) = --------------------------------
!                     2 dx
      fx(i,j,k)=d2dx*(-fh(i-1,j,k)+fh(i+1,j,k))

! set imax and imin 0
    endif

  enddo
  enddo
  enddo

  return

  end subroutine fdx_sh
!-----------------------------------------------------------------------------
!
! single derivatives dy
!
!-----------------------------------------------------------------------------
  subroutine fdy_sh(ex,f,fy,X,Y,Z,SYM1,SYM2,SYM3,symmetry,onoff,sst)
  implicit none

  integer,                               intent(in ):: ex(1:3),symmetry,onoff,sst
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in ):: f
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out):: fy
  real*8,                                intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8,                                intent(in ):: SYM1,SYM2,SYM3

!~~~~~~ other variables

  real*8 :: dX,dY,dZ
  real*8,dimension(-3:ex(1)+4,-3:ex(2)+4,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8 :: d840dy,d60dy,d12dy,d2dy
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8,  parameter :: ZEO=0.d0,ONE=1.d0, F60=6.d1, F32 = 3.2d1
  real*8,  parameter :: TWO=2.d0,THR=3.d0, EIT=8.d0, F168=1.68d2
  real*8,  parameter ::  F9=9.d0,F45=4.5d1,F12=1.2d1,F672=6.72d2
  real*8,  parameter :: F840=8.4d2

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)
  
  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -3
     if(dabs(Y(1)) < dY) jmin = -3
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -3
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+4
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(4,ex,f,fh,SoA)

  d840dy = ONE/F840/dY
  
  d60dy = ONE/F60/dY

  d12dy = ONE/F12/dY

  d2dy = ONE/TWO/dY

  fy = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
! y direction   
        if(j+4 <= jmax .and. j-4 >= jmin)then

      fy(i,j,k)=d840dy*( THR*fh(i,j-4,k)-F32 *fh(i,j-3,k)+F168*fh(i,j-2,k)-F672*fh(i,j-1,k)+ &
                        F672*fh(i,j+1,k)-F168*fh(i,j+2,k)+F32 *fh(i,j+3,k)-THR *fh(i,j+4,k))

    elseif(j+3 <= jmax .and. j-3 >= jmin)then

      fy(i,j,k)=d60dy*(-fh(i,j-3,k)+F9*fh(i,j-2,k)-F45*fh(i,j-1,k)+F45*fh(i,j+1,k)-F9*fh(i,j+2,k)+fh(i,j+3,k))

    elseif(j+2 <= jmax .and. j-2 >= jmin)then

      fy(i,j,k)=d12dy*(fh(i,j-2,k)-EIT*fh(i,j-1,k)+EIT*fh(i,j+1,k)-fh(i,j+2,k))

    elseif(j+1 <= jmax .and. j-1 >= jmin)then

     fy(i,j,k)=d2dy*(-fh(i,j-1,k)+fh(i,j+1,k))

! set jmax and jmin 0
    endif

  enddo
  enddo
  enddo

  return

  end subroutine fdy_sh
!-----------------------------------------------------------------------------
!
! single derivatives dz
!
!-----------------------------------------------------------------------------
  subroutine fdz_sh(ex,f,fz,X,Y,Z,SYM1,SYM2,SYM3,symmetry,onoff,sst)
  implicit none

  integer,                               intent(in ):: ex(1:3),symmetry,onoff,sst
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in ):: f
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out):: fz
  real*8,                                intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8,                                intent(in ):: SYM1,SYM2,SYM3

!~~~~~~ other variables
  
  real*8 :: dX,dY,dZ
  real*8,dimension(-3:ex(1)+4,-3:ex(2)+4,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8 :: d840dz,d60dz,d12dz,d2dz
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8,  parameter :: ZEO=0.d0,ONE=1.d0, F60=6.d1, F32 = 3.2d1
  real*8,  parameter :: TWO=2.d0,THR=3.d0, EIT=8.d0, F168=1.68d2
  real*8,  parameter ::  F9=9.d0,F45=4.5d1,F12=1.2d1,F672=6.72d2
  real*8,  parameter :: F840=8.4d2

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)
  
  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -3
     if(dabs(Y(1)) < dY) jmin = -3
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -3
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+4
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(4,ex,f,fh,SoA)

  d840dz = ONE/F840/dZ
  
  d60dz = ONE/F60/dZ

  d12dz = ONE/F12/dZ

  d2dz = ONE/TWO/dZ

  fz = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
! z direction   
        if(k+4 <= kmax .and. k-4 >= kmin)then

      fz(i,j,k)=d840dz*( THR*fh(i,j,k-4)-F32 *fh(i,j,k-3)+F168*fh(i,j,k-2)-F672*fh(i,j,k-1)+ &
                        F672*fh(i,j,k+1)-F168*fh(i,j,k+2)+F32 *fh(i,j,k+3)-THR *fh(i,j,k+4))

    elseif(k+3 <= kmax .and. k-3 >= kmin)then

      fz(i,j,k)=d60dz*(-fh(i,j,k-3)+F9*fh(i,j,k-2)-F45*fh(i,j,k-1)+F45*fh(i,j,k+1)-F9*fh(i,j,k+2)+fh(i,j,k+3))

    elseif(k+2 <= kmax .and. k-2 >= kmin)then

      fz(i,j,k)=d12dz*(fh(i,j,k-2)-EIT*fh(i,j,k-1)+EIT*fh(i,j,k+1)-fh(i,j,k+2))

    elseif(k+1 <= kmax .and. k-1 >= kmin)then

      fz(i,j,k)=d2dz*(-fh(i,j,k-1)+fh(i,j,k+1))

! set kmax and kmin 0
    endif

  enddo
  enddo
  enddo

  return

  end subroutine fdz_sh
!-----------------------------------------------------------------------------------------------------------------
!
! General second derivatives of 8_th oder accurate
!
!            - 9 f(i-4) + 128 f(i-3) - 1008 f(i-2) + 8064 f(i-1) - 14350 f(i) + 8064 f(i+1) - 1008 f(i+2) + 128 f(i+3) - 9 f(i+4)
!  fxx(i) = ----------------------------------------------------------------------------------------------------------------------
!                                                                   5040 dx^2 
!
!             + 3   ( 3 f(i-4,j-4) - 32 f(i-3,j-4) + 168 f(i-2,j-4) - 672 f(i-1,j-4) + 672 f(i+1,j-4) - 168 f(i+2,j-4) + 32 f(i+3,j-4) - 3 f(i+4,j-4) )
!             - 32  ( 3 f(i-4,j-3) - 32 f(i-3,j-3) + 168 f(i-2,j-3) - 672 f(i-1,j-3) + 672 f(i+1,j-3) - 168 f(i+2,j-3) + 32 f(i+3,j-3) - 3 f(i+4,j-3) )
!             + 168 ( 3 f(i-4,j-2) - 32 f(i-3,j-2) + 168 f(i-2,j-2) - 672 f(i-1,j-2) + 672 f(i+1,j-2) - 168 f(i+2,j-2) + 32 f(i+3,j-2) - 3 f(i+4,j-2) )
!             - 672 ( 3 f(i-4,j-1) - 32 f(i-3,j-1) + 168 f(i-2,j-1) - 672 f(i-1,j-1) + 672 f(i+1,j-1) - 168 f(i+2,j-1) + 32 f(i+3,j-1) - 3 f(i+4,j-1) )
!             + 672 ( 3 f(i-4,j+1) - 32 f(i-3,j+1) + 168 f(i-2,j+1) - 672 f(i-1,j+1) + 672 f(i+1,j+1) - 168 f(i+2,j+1) + 32 f(i+3,j+1) - 3 f(i+4,j+1) )
!             - 168 ( 3 f(i-4,j+2) - 32 f(i-3,j+2) + 168 f(i-2,j+2) - 672 f(i-1,j+2) + 672 f(i+1,j+2) - 168 f(i+2,j+2) + 32 f(i+3,j+2) - 3 f(i+4,j+2) )
!             + 32  ( 3 f(i-4,j+3) - 32 f(i-3,j+3) + 168 f(i-2,j+3) - 672 f(i-1,j+3) + 672 f(i+1,j+3) - 168 f(i+2,j+3) + 32 f(i+3,j+3) - 3 f(i+4,j+3) )
!             - 3   ( 3 f(i-4,j+4) - 32 f(i-3,j+4) + 168 f(i-2,j+4) - 672 f(i-1,j+4) + 672 f(i+1,j+4) - 168 f(i+2,j+4) + 32 f(i+3,j+4) - 3 f(i+4,j+4) )
!  fxy(i,j) = ------------------------------------------------------------------------------------------------------------------------------------------
!                                                                              705600 dx dy
!
!-----------------------------------------------------------------------------------------------------------------
  subroutine fdderivs_sh(ex,f,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z, &
                      SYM1,SYM2,SYM3,symmetry,onoff,sst)
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry,onoff,sst
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fxx,fxy,fxz,fyy,fyz,fzz
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8,dimension(-3:ex(1)+4,-3:ex(2)+4,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdxdx,Sdydy,Sdzdz,Fdxdx,Fdydy,Fdzdz,Xdxdx,Xdydy,Xdzdz,Edxdx,Edydy,Edzdz
  real*8  :: Sdxdy,Sdxdz,Sdydz,Fdxdy,Fdxdz,Fdydz,Xdxdy,Xdxdz,Xdydz,Edxdy,Edxdz,Edydz
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1, F128=1.28d2
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2,F1008=1.008d3
  real*8, parameter :: F8064=8.064d3,F14350=1.435d4,THR=3.d0,F32=3.2d1,F168=1.68d2,F672=6.72d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3
  real*8, parameter :: F1o5040=ONE/5.04d3,F1o705600=ONE/7.056d5

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -3
     if(dabs(Y(1)) < dY) jmin = -3
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -3
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+4
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(4,ex,f,fh,SoA)

  Sdxdx =  ONE /( dX * dX )
  Sdydy =  ONE /( dY * dY )
  Sdzdz =  ONE /( dZ * dZ )

  Fdxdx = F1o12 /( dX * dX )
  Fdydy = F1o12 /( dY * dY )
  Fdzdz = F1o12 /( dZ * dZ )

  Xdxdx = F1o180 /( dX * dX )
  Xdydy = F1o180 /( dY * dY )
  Xdzdz = F1o180 /( dZ * dZ )

  Edxdx = F1o5040 /( dX * dX )
  Edydy = F1o5040 /( dY * dY )
  Edzdz = F1o5040 /( dZ * dZ )

  Sdxdy = F1o4 /( dX * dY )
  Sdxdz = F1o4 /( dX * dZ )
  Sdydz = F1o4 /( dY * dZ )

  Fdxdy = F1o144 /( dX * dY )
  Fdxdz = F1o144 /( dX * dZ )
  Fdydz = F1o144 /( dY * dZ )

  Xdxdy = F1o3600 /( dX * dY )
  Xdxdz = F1o3600 /( dX * dZ )
  Xdydz = F1o3600 /( dY * dZ )

  Edxdy = F1o705600 /( dX * dY )
  Edxdz = F1o705600 /( dX * dZ )
  Edydz = F1o705600 /( dY * dZ )

  fxx = ZEO
  fyy = ZEO
  fzz = ZEO
  fxy = ZEO
  fxz = ZEO
  fyz = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
!~~~~~~ fxx
        if(i+4 <= imax .and. i-4 >= imin)then

!            - 9 f(i-4) + 128 f(i-3) - 1008 f(i-2) + 8064 f(i-1) - 14350 f(i) + 8064 f(i+1) - 1008 f(i+2) + 128 f(i+3) - 9 f(i+4)
!  fxx(i) = ----------------------------------------------------------------------------------------------------------------------
!                                                                   5040 dx^2 
   fxx(i,j,k) = Edxdx*(-F9*fh(i-4,j,k)+F128*fh(i-3,j,k)-F1008*fh(i-2,j,k)+F8064*fh(i-1,j,k)-F14350*fh(i,j,k) &
                       -F9*fh(i+4,j,k)+F128*fh(i+3,j,k)-F1008*fh(i+2,j,k)+F8064*fh(i+1,j,k)               )

   elseif(i+3 <= imax .and. i-3 >= imin)then

!             2 f(i-3) - 27 f(i-2) + 270 f(i-1) - 490 f(i) + 270 f(i+1) - 27 f(i+2) + 2 f(i+3)
!  fxx(i) = -----------------------------------------------------------------------------------
!                                                180 dx^2 
   fxx(i,j,k) = Xdxdx*(TWO*fh(i-3,j,k)-F27*fh(i-2,j,k)+F270*fh(i-1,j,k)-F490*fh(i,j,k) &
                      +TWO*fh(i+3,j,k)-F27*fh(i+2,j,k)+F270*fh(i+1,j,k)               )
   elseif(i+2 <= imax .and. i-2 >= imin)then
!
!               - f(i-2) + 16 f(i-1) - 30 f(i) + 16 f(i+1) - f(i+2)
!  fxx(i) = ----------------------------------------------------------
!                                  12 dx^2 
   fxx(i,j,k) = Fdxdx*(-fh(i-2,j,k)+F16*fh(i-1,j,k)-F30*fh(i,j,k) &
                       -fh(i+2,j,k)+F16*fh(i+1,j,k)              )
   elseif(i+1 <= imax .and. i-1 >= imin)then
!
!               f(i-1) - 2 f(i) + f(i+1)
!  fxx(i) = --------------------------------
!                         dx^2 
   fxx(i,j,k) = Sdxdx*(fh(i-1,j,k)-TWO*fh(i,j,k) &
                      +fh(i+1,j,k)              )
   endif

!~~~~~~ fyy
        if(j+4 <= jmax .and. j-4 >= jmin)then

   fyy(i,j,k) = Edydy*(-F9*fh(i,j-4,k)+F128*fh(i,j-3,k)-F1008*fh(i,j-2,k)+F8064*fh(i,j-1,k)-F14350*fh(i,j,k) &
                       -F9*fh(i,j+4,k)+F128*fh(i,j+3,k)-F1008*fh(i,j+2,k)+F8064*fh(i,j+1,k)               )

   elseif(j+3 <= jmax .and. j-3 >= jmin)then

   fyy(i,j,k) = Xdydy*(TWO*fh(i,j-3,k)-F27*fh(i,j-2,k)+F270*fh(i,j-1,k)-F490*fh(i,j,k) &
                      +TWO*fh(i,j+3,k)-F27*fh(i,j+2,k)+F270*fh(i,j+1,k)               )
   elseif(j+2 <= jmax .and. j-2 >= jmin)then

   fyy(i,j,k) = Fdydy*(-fh(i,j-2,k)+F16*fh(i,j-1,k)-F30*fh(i,j,k) &
                       -fh(i,j+2,k)+F16*fh(i,j+1,k)              )
   elseif(j+1 <= jmax .and. j-1 >= jmin)then

   fyy(i,j,k) = Sdydy*(fh(i,j-1,k)-TWO*fh(i,j,k) &
                      +fh(i,j+1,k)              )
   endif

!~~~~~~ fzz
        if(k+4 <= kmax .and. k-4 >= kmin)then

   fzz(i,j,k) = Edzdz*(-F9*fh(i,j,k-4)+F128*fh(i,j,k-3)-F1008*fh(i,j,k-2)+F8064*fh(i,j,k-1)-F14350*fh(i,j,k) &
                       -F9*fh(i,j,k+4)+F128*fh(i,j,k+3)-F1008*fh(i,j,k+2)+F8064*fh(i,j,k+1)               )

   elseif(k+3 <= kmax .and. k-3 >= kmin)then

   fzz(i,j,k) = Xdzdz*(TWO*fh(i,j,k-3)-F27*fh(i,j,k-2)+F270*fh(i,j,k-1)-F490*fh(i,j,k) &
                      +TWO*fh(i,j,k+3)-F27*fh(i,j,k+2)+F270*fh(i,j,k+1)               )
   elseif(k+2 <= kmax .and. k-2 >= kmin)then

   fzz(i,j,k) = Fdzdz*(-fh(i,j,k-2)+F16*fh(i,j,k-1)-F30*fh(i,j,k) &
                       -fh(i,j,k+2)+F16*fh(i,j,k+1)              )
   elseif(k+1 <= kmax .and. k-1 >= kmin)then

   fzz(i,j,k) = Sdzdz*(fh(i,j,k-1)-TWO*fh(i,j,k) &
                      +fh(i,j,k+1)              )
   endif
!~~~~~~ fxy
      if(i+4 <= imax .and. i-4 >= imin .and. j+4 <= jmax .and. j-4 >= jmin)then

!             + 3   ( 3 f(i-4,j-4) - 32 f(i-3,j-4) + 168 f(i-2,j-4) - 672 f(i-1,j-4) + 672 f(i+1,j-4) - 168 f(i+2,j-4) + 32 f(i+3,j-4) - 3 f(i+4,j-4) )
!             - 32  ( 3 f(i-4,j-3) - 32 f(i-3,j-3) + 168 f(i-2,j-3) - 672 f(i-1,j-3) + 672 f(i+1,j-3) - 168 f(i+2,j-3) + 32 f(i+3,j-3) - 3 f(i+4,j-3) )
!             + 168 ( 3 f(i-4,j-2) - 32 f(i-3,j-2) + 168 f(i-2,j-2) - 672 f(i-1,j-2) + 672 f(i+1,j-2) - 168 f(i+2,j-2) + 32 f(i+3,j-2) - 3 f(i+4,j-2) )
!             - 672 ( 3 f(i-4,j-1) - 32 f(i-3,j-1) + 168 f(i-2,j-1) - 672 f(i-1,j-1) + 672 f(i+1,j-1) - 168 f(i+2,j-1) + 32 f(i+3,j-1) - 3 f(i+4,j-1) )
!             + 672 ( 3 f(i-4,j+1) - 32 f(i-3,j+1) + 168 f(i-2,j+1) - 672 f(i-1,j+1) + 672 f(i+1,j+1) - 168 f(i+2,j+1) + 32 f(i+3,j+1) - 3 f(i+4,j+1) )
!             - 168 ( 3 f(i-4,j+2) - 32 f(i-3,j+2) + 168 f(i-2,j+2) - 672 f(i-1,j+2) + 672 f(i+1,j+2) - 168 f(i+2,j+2) + 32 f(i+3,j+2) - 3 f(i+4,j+2) )
!             + 32  ( 3 f(i-4,j+3) - 32 f(i-3,j+3) + 168 f(i-2,j+3) - 672 f(i-1,j+3) + 672 f(i+1,j+3) - 168 f(i+2,j+3) + 32 f(i+3,j+3) - 3 f(i+4,j+3) )
!             - 3   ( 3 f(i-4,j+4) - 32 f(i-3,j+4) + 168 f(i-2,j+4) - 672 f(i-1,j+4) + 672 f(i+1,j+4) - 168 f(i+2,j+4) + 32 f(i+3,j+4) - 3 f(i+4,j+4) )
!  fxy(i,j) = ------------------------------------------------------------------------------------------------------------------------------------------
!                                                                              705600 dx dy
   fxy(i,j,k) = Edxdy*( THR *( THR*fh(i-4,j-4,k)-F32*fh(i-3,j-4,k)+F168*fh(i-2,j-4,k)-F672*fh(i-1,j-4,k)        &
                              -THR*fh(i+4,j-4,k)+F32*fh(i+3,j-4,k)-F168*fh(i+2,j-4,k)+F672*fh(i+1,j-4,k))       &
                       -F32 *( THR*fh(i-4,j-3,k)-F32*fh(i-3,j-3,k)+F168*fh(i-2,j-3,k)-F672*fh(i-1,j-3,k)        &
                              -THR*fh(i+4,j-3,k)+F32*fh(i+3,j-3,k)-F168*fh(i+2,j-3,k)+F672*fh(i+1,j-3,k))       &
                       +F168*( THR*fh(i-4,j-2,k)-F32*fh(i-3,j-2,k)+F168*fh(i-2,j-2,k)-F672*fh(i-1,j-2,k)        &
                              -THR*fh(i+4,j-2,k)+F32*fh(i+3,j-2,k)-F168*fh(i+2,j-2,k)+F672*fh(i+1,j-2,k))       &
                       -F672*( THR*fh(i-4,j-1,k)-F32*fh(i-3,j-1,k)+F168*fh(i-2,j-1,k)-F672*fh(i-1,j-1,k)        &
                              -THR*fh(i+4,j-1,k)+F32*fh(i+3,j-1,k)-F168*fh(i+2,j-1,k)+F672*fh(i+1,j-1,k))       &
                       +F672*( THR*fh(i-4,j+1,k)-F32*fh(i-3,j+1,k)+F168*fh(i-2,j+1,k)-F672*fh(i-1,j+1,k)        &
                              -THR*fh(i+4,j+1,k)+F32*fh(i+3,j+1,k)-F168*fh(i+2,j+1,k)+F672*fh(i+1,j+1,k))       &
                       -F168*( THR*fh(i-4,j+2,k)-F32*fh(i-3,j+2,k)+F168*fh(i-2,j+2,k)-F672*fh(i-1,j+2,k)        &
                              -THR*fh(i+4,j+2,k)+F32*fh(i+3,j+2,k)-F168*fh(i+2,j+2,k)+F672*fh(i+1,j+2,k))       &
                       +F32 *( THR*fh(i-4,j+3,k)-F32*fh(i-3,j+3,k)+F168*fh(i-2,j+3,k)-F672*fh(i-1,j+3,k)        &
                              -THR*fh(i+4,j+3,k)+F32*fh(i+3,j+3,k)-F168*fh(i+2,j+3,k)+F672*fh(i+1,j+3,k))       &
                       -THR *( THR*fh(i-4,j+4,k)-F32*fh(i-3,j+4,k)+F168*fh(i-2,j+4,k)-F672*fh(i-1,j+4,k)        &
                              -THR*fh(i+4,j+4,k)+F32*fh(i+3,j+4,k)-F168*fh(i+2,j+4,k)+F672*fh(i+1,j+4,k)) )
   elseif(i+3 <= imax .and. i-3 >= imin .and. j+3 <= jmax .and. j-3 >= jmin)then
!
!             -    ( - f(i-3,j-3) + 9 f(i-2,j-3) - 45 f(i-1,j-3) + 45 f(i+1,j-3) - 9 f(i+2,j-3) + f(i+3,j-3) )
!             + 9  ( - f(i-3,j-2) + 9 f(i-2,j-2) - 45 f(i-1,j-2) + 45 f(i+1,j-2) - 9 f(i+2,j-2) + f(i+3,j-2) )
!             - 45 ( - f(i-3,j-1) + 9 f(i-2,j-1) - 45 f(i-1,j-1) + 45 f(i+1,j-1) - 9 f(i+2,j-1) + f(i+3,j-1) )
!             + 45 ( - f(i-3,j+1) + 9 f(i-2,j+1) - 45 f(i-1,j+1) + 45 f(i+1,j+1) - 9 f(i+2,j+1) + f(i+3,j+1) )
!             - 9  ( - f(i-3,j+2) + 9 f(i-2,j+2) - 45 f(i-1,j+2) + 45 f(i+1,j+2) - 9 f(i+2,j+2) + f(i+3,j+2) )
!             +    ( - f(i-3,j+3) + 9 f(i-2,j+3) - 45 f(i-1,j+3) + 45 f(i+1,j+3) - 9 f(i+2,j+3) + f(i+3,j+3) )  
!  fxy(i,j) = ------------------------------------------------------------------------------------------------
!                                                          3600 dx dy
   fxy(i,j,k) = Xdxdy*(-    (-fh(i-3,j-3,k)+F9*fh(i-2,j-3,k)-F45*fh(i-1,j-3,k)+F45*fh(i+1,j-3,k)-F9*fh(i+2,j-3,k)+fh(i+3,j-3,k))  &
                       +F9 *(-fh(i-3,j-2,k)+F9*fh(i-2,j-2,k)-F45*fh(i-1,j-2,k)+F45*fh(i+1,j-2,k)-F9*fh(i+2,j-2,k)+fh(i+3,j-2,k))  &
                       -F45*(-fh(i-3,j-1,k)+F9*fh(i-2,j-1,k)-F45*fh(i-1,j-1,k)+F45*fh(i+1,j-1,k)-F9*fh(i+2,j-1,k)+fh(i+3,j-1,k))  &
                       +F45*(-fh(i-3,j+1,k)+F9*fh(i-2,j+1,k)-F45*fh(i-1,j+1,k)+F45*fh(i+1,j+1,k)-F9*fh(i+2,j+1,k)+fh(i+3,j+1,k))  &
                       -F9 *(-fh(i-3,j+2,k)+F9*fh(i-2,j+2,k)-F45*fh(i-1,j+2,k)+F45*fh(i+1,j+2,k)-F9*fh(i+2,j+2,k)+fh(i+3,j+2,k))  &
                       +    (-fh(i-3,j+3,k)+F9*fh(i-2,j+3,k)-F45*fh(i-1,j+3,k)+F45*fh(i+1,j+3,k)-F9*fh(i+2,j+3,k)+fh(i+3,j+3,k)))
   elseif(i+2 <= imax .and. i-2 >= imin .and. j+2 <= jmax .and. j-2 >= jmin)then
!
!                 ( f(i-2,j-2) - 8 f(i-1,j-2) + 8 f(i+1,j-2) - f(i+2,j-2) )
!             - 8 ( f(i-2,j-1) - 8 f(i-1,j-1) + 8 f(i+1,j-1) - f(i+2,j-1) )
!             + 8 ( f(i-2,j+1) - 8 f(i-1,j+1) + 8 f(i+1,j+1) - f(i+2,j+1) )
!             -   ( f(i-2,j+2) - 8 f(i-1,j+2) + 8 f(i+1,j+2) - f(i+2,j+2) )
!  fxy(i,j) = ----------------------------------------------------------------
!                                  144 dx dy
   fxy(i,j,k) = Fdxdy*(     (fh(i-2,j-2,k)-F8*fh(i-1,j-2,k)+F8*fh(i+1,j-2,k)-fh(i+2,j-2,k))  &
                       -F8 *(fh(i-2,j-1,k)-F8*fh(i-1,j-1,k)+F8*fh(i+1,j-1,k)-fh(i+2,j-1,k))  &
                       +F8 *(fh(i-2,j+1,k)-F8*fh(i-1,j+1,k)+F8*fh(i+1,j+1,k)-fh(i+2,j+1,k))  &
                       -    (fh(i-2,j+2,k)-F8*fh(i-1,j+2,k)+F8*fh(i+1,j+2,k)-fh(i+2,j+2,k)))

   elseif(i+1 <= imax .and. i-1 >= imin .and. j+1 <= jmax .and. j-1 >= jmin)then
!                 f(i-1,j-1) - f(i+1,j-1) - f(i-1,j+1) + f(i+1,j+1) 
!  fxy(i,j) = -----------------------------------------------------------
!                                      4 dx dy
   fxy(i,j,k) = Sdxdy*(fh(i-1,j-1,k)-fh(i+1,j-1,k)-fh(i-1,j+1,k)+fh(i+1,j+1,k))
   endif
!~~~~~~ fxz
      if(i+4 <= imax .and. i-4 >= imin .and. k+4 <= kmax .and. k-4 >= kmin)then

   fxz(i,j,k) = Edxdz*( THR *( THR*fh(i-4,j,k-4)-F32*fh(i-3,j,k-4)+F168*fh(i-2,j,k-4)-F672*fh(i-1,j,k-4)        &
                              -THR*fh(i+4,j,k-4)+F32*fh(i+3,j,k-4)-F168*fh(i+2,j,k-4)+F672*fh(i+1,j,k-4))       &
                       -F32 *( THR*fh(i-4,j,k-3)-F32*fh(i-3,j,k-3)+F168*fh(i-2,j,k-3)-F672*fh(i-1,j,k-3)        &
                              -THR*fh(i+4,j,k-3)+F32*fh(i+3,j,k-3)-F168*fh(i+2,j,k-3)+F672*fh(i+1,j,k-3))       &
                       +F168*( THR*fh(i-4,j,k-2)-F32*fh(i-3,j,k-2)+F168*fh(i-2,j,k-2)-F672*fh(i-1,j,k-2)        &
                              -THR*fh(i+4,j,k-2)+F32*fh(i+3,j,k-2)-F168*fh(i+2,j,k-2)+F672*fh(i+1,j,k-2))       &
                       -F672*( THR*fh(i-4,j,k-1)-F32*fh(i-3,j,k-1)+F168*fh(i-2,j,k-1)-F672*fh(i-1,j,k-1)        &
                              -THR*fh(i+4,j,k-1)+F32*fh(i+3,j,k-1)-F168*fh(i+2,j,k-1)+F672*fh(i+1,j,k-1))       &
                       +F672*( THR*fh(i-4,j,k+1)-F32*fh(i-3,j,k+1)+F168*fh(i-2,j,k+1)-F672*fh(i-1,j,k+1)        &
                              -THR*fh(i+4,j,k+1)+F32*fh(i+3,j,k+1)-F168*fh(i+2,j,k+1)+F672*fh(i+1,j,k+1))       &
                       -F168*( THR*fh(i-4,j,k+2)-F32*fh(i-3,j,k+2)+F168*fh(i-2,j,k+2)-F672*fh(i-1,j,k+2)        &
                              -THR*fh(i+4,j,k+2)+F32*fh(i+3,j,k+2)-F168*fh(i+2,j,k+2)+F672*fh(i+1,j,k+2))       &
                       +F32 *( THR*fh(i-4,j,k+3)-F32*fh(i-3,j,k+3)+F168*fh(i-2,j,k+3)-F672*fh(i-1,j,k+3)        &
                              -THR*fh(i+4,j,k+3)+F32*fh(i+3,j,k+3)-F168*fh(i+2,j,k+3)+F672*fh(i+1,j,k+3))       &
                       -THR *( THR*fh(i-4,j,k+4)-F32*fh(i-3,j,k+4)+F168*fh(i-2,j,k+4)-F672*fh(i-1,j,k+4)        &
                              -THR*fh(i+4,j,k+4)+F32*fh(i+3,j,k+4)-F168*fh(i+2,j,k+4)+F672*fh(i+1,j,k+4)) )
   elseif(i+3 <= imax .and. i-3 >= imin .and. k+3 <= kmax .and. k-3 >= kmin)then

   fxz(i,j,k) = Xdxdz*(-    (-fh(i-3,j,k-3)+F9*fh(i-2,j,k-3)-F45*fh(i-1,j,k-3)+F45*fh(i+1,j,k-3)-F9*fh(i+2,j,k-3)+fh(i+3,j,k-3))  &
                       +F9 *(-fh(i-3,j,k-2)+F9*fh(i-2,j,k-2)-F45*fh(i-1,j,k-2)+F45*fh(i+1,j,k-2)-F9*fh(i+2,j,k-2)+fh(i+3,j,k-2))  &
                       -F45*(-fh(i-3,j,k-1)+F9*fh(i-2,j,k-1)-F45*fh(i-1,j,k-1)+F45*fh(i+1,j,k-1)-F9*fh(i+2,j,k-1)+fh(i+3,j,k-1))  &
                       +F45*(-fh(i-3,j,k+1)+F9*fh(i-2,j,k+1)-F45*fh(i-1,j,k+1)+F45*fh(i+1,j,k+1)-F9*fh(i+2,j,k+1)+fh(i+3,j,k+1))  &
                       -F9 *(-fh(i-3,j,k+2)+F9*fh(i-2,j,k+2)-F45*fh(i-1,j,k+2)+F45*fh(i+1,j,k+2)-F9*fh(i+2,j,k+2)+fh(i+3,j,k+2))  &
                       +    (-fh(i-3,j,k+3)+F9*fh(i-2,j,k+3)-F45*fh(i-1,j,k+3)+F45*fh(i+1,j,k+3)-F9*fh(i+2,j,k+3)+fh(i+3,j,k+3)))
   elseif(i+2 <= imax .and. i-2 >= imin .and. k+2 <= kmax .and. k-2 >= kmin)then
   fxz(i,j,k) = Fdxdz*(     (fh(i-2,j,k-2)-F8*fh(i-1,j,k-2)+F8*fh(i+1,j,k-2)-fh(i+2,j,k-2))  &
                       -F8 *(fh(i-2,j,k-1)-F8*fh(i-1,j,k-1)+F8*fh(i+1,j,k-1)-fh(i+2,j,k-1))  &
                       +F8 *(fh(i-2,j,k+1)-F8*fh(i-1,j,k+1)+F8*fh(i+1,j,k+1)-fh(i+2,j,k+1))  &
                       -    (fh(i-2,j,k+2)-F8*fh(i-1,j,k+2)+F8*fh(i+1,j,k+2)-fh(i+2,j,k+2)))
   elseif(i+1 <= imax .and. i-1 >= imin .and. k+1 <= kmax .and. k-1 >= kmin)then
   fxz(i,j,k) = Sdxdz*(fh(i-1,j,k-1)-fh(i+1,j,k-1)-fh(i-1,j,k+1)+fh(i+1,j,k+1))
   endif
!~~~~~~ fyz
      if(j+4 <= jmax .and. j-4 >= jmin .and. k+4 <= kmax .and. k-4 >= kmin)then

   fyz(i,j,k) = Edydz*( THR *( THR*fh(i,j-4,k-4)-F32*fh(i,j-3,k-4)+F168*fh(i,j-2,k-4)-F672*fh(i,j-1,k-4)        &
                              -THR*fh(i,j+4,k-4)+F32*fh(i,j+3,k-4)-F168*fh(i,j+2,k-4)+F672*fh(i,j+1,k-4))       &
                       -F32 *( THR*fh(i,j-4,k-3)-F32*fh(i,j-3,k-3)+F168*fh(i,j-2,k-3)-F672*fh(i,j-1,k-3)        &
                              -THR*fh(i,j+4,k-3)+F32*fh(i,j+3,k-3)-F168*fh(i,j+2,k-3)+F672*fh(i,j+1,k-3))       &
                       +F168*( THR*fh(i,j-4,k-2)-F32*fh(i,j-3,k-2)+F168*fh(i,j-2,k-2)-F672*fh(i,j-1,k-2)        &
                              -THR*fh(i,j+4,k-2)+F32*fh(i,j+3,k-2)-F168*fh(i,j+2,k-2)+F672*fh(i,j+1,k-2))       &
                       -F672*( THR*fh(i,j-4,k-1)-F32*fh(i,j-3,k-1)+F168*fh(i,j-2,k-1)-F672*fh(i,j-1,k-1)        &
                              -THR*fh(i,j+4,k-1)+F32*fh(i,j+3,k-1)-F168*fh(i,j+2,k-1)+F672*fh(i,j+1,k-1))       &
                       +F672*( THR*fh(i,j-4,k+1)-F32*fh(i,j-3,k+1)+F168*fh(i,j-2,k+1)-F672*fh(i,j-1,k+1)        &
                              -THR*fh(i,j+4,k+1)+F32*fh(i,j+3,k+1)-F168*fh(i,j+2,k+1)+F672*fh(i,j+1,k+1))       &
                       -F168*( THR*fh(i,j-4,k+2)-F32*fh(i,j-3,k+2)+F168*fh(i,j-2,k+2)-F672*fh(i,j-1,k+2)        &
                              -THR*fh(i,j+4,k+2)+F32*fh(i,j+3,k+2)-F168*fh(i,j+2,k+2)+F672*fh(i,j+1,k+2))       &
                       +F32 *( THR*fh(i,j-4,k+3)-F32*fh(i,j-3,k+3)+F168*fh(i,j-2,k+3)-F672*fh(i,j-1,k+3)        &
                              -THR*fh(i,j+4,k+3)+F32*fh(i,j+3,k+3)-F168*fh(i,j+2,k+3)+F672*fh(i,j+1,k+3))       &
                       -THR *( THR*fh(i,j-4,k+4)-F32*fh(i,j-3,k+4)+F168*fh(i,j-2,k+4)-F672*fh(i,j-1,k+4)        &
                              -THR*fh(i,j+4,k+4)+F32*fh(i,j+3,k+4)-F168*fh(i,j+2,k+4)+F672*fh(i,j+1,k+4)) )
   elseif(j+3 <= jmax .and. j-3 >= jmin .and. k+3 <= kmax .and. k-3 >= kmin)then

   fyz(i,j,k) = Xdydz*(-    (-fh(i,j-3,k-3)+F9*fh(i,j-2,k-3)-F45*fh(i,j-1,k-3)+F45*fh(i,j+1,k-3)-F9*fh(i,j+2,k-3)+fh(i,j+3,k-3))  &
                       +F9 *(-fh(i,j-3,k-2)+F9*fh(i,j-2,k-2)-F45*fh(i,j-1,k-2)+F45*fh(i,j+1,k-2)-F9*fh(i,j+2,k-2)+fh(i,j+3,k-2))  &
                       -F45*(-fh(i,j-3,k-1)+F9*fh(i,j-2,k-1)-F45*fh(i,j-1,k-1)+F45*fh(i,j+1,k-1)-F9*fh(i,j+2,k-1)+fh(i,j+3,k-1))  &
                       +F45*(-fh(i,j-3,k+1)+F9*fh(i,j-2,k+1)-F45*fh(i,j-1,k+1)+F45*fh(i,j+1,k+1)-F9*fh(i,j+2,k+1)+fh(i,j+3,k+1))  &
                       -F9 *(-fh(i,j-3,k+2)+F9*fh(i,j-2,k+2)-F45*fh(i,j-1,k+2)+F45*fh(i,j+1,k+2)-F9*fh(i,j+2,k+2)+fh(i,j+3,k+2))  &
                       +    (-fh(i,j-3,k+3)+F9*fh(i,j-2,k+3)-F45*fh(i,j-1,k+3)+F45*fh(i,j+1,k+3)-F9*fh(i,j+2,k+3)+fh(i,j+3,k+3)))
   elseif(j+2 <= jmax .and. j-2 >= jmin .and. k+2 <= kmax .and. k-2 >= kmin)then
   fyz(i,j,k) = Fdydz*(     (fh(i,j-2,k-2)-F8*fh(i,j-1,k-2)+F8*fh(i,j+1,k-2)-fh(i,j+2,k-2))  &
                       -F8 *(fh(i,j-2,k-1)-F8*fh(i,j-1,k-1)+F8*fh(i,j+1,k-1)-fh(i,j+2,k-1))  &
                       +F8 *(fh(i,j-2,k+1)-F8*fh(i,j-1,k+1)+F8*fh(i,j+1,k+1)-fh(i,j+2,k+1))  &
                       -    (fh(i,j-2,k+2)-F8*fh(i,j-1,k+2)+F8*fh(i,j+1,k+2)-fh(i,j+2,k+2)))
   elseif(j+1 <= jmax .and. j-1 >= jmin .and. k+1 <= kmax .and. k-1 >= kmin)then
   fyz(i,j,k) = Sdydz*(fh(i,j-1,k-1)-fh(i,j+1,k-1)-fh(i,j-1,k+1)+fh(i,j+1,k+1))
   endif 

   enddo
   enddo
   enddo

  return

  end subroutine fdderivs_sh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! only for compute_ricci.f90 usage
!-----------------------------------------------------------------------------
  subroutine fddxx_sh(ex,f,fxx,X,Y,Z,SYM1,SYM2,SYM3,symmetry,sst)
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry,sst
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fxx
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8,dimension(-3:ex(1)+4,-3:ex(2)+4,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdxdx,Fdxdx,Xdxdx,Edxdx
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1, F128=1.28d2
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2,F1008=1.008d3
  real*8, parameter :: F8064=8.064d3,F14350=1.435d4,THR=3.d0,F32=3.2d1,F168=1.68d2,F672=6.72d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3
  real*8, parameter :: F1o5040=ONE/5.04d3,F1o705600=ONE/7.056d5

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -3
     if(dabs(Y(1)) < dY) jmin = -3
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -3
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+4
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(4,ex,f,fh,SoA)

  Sdxdx =  ONE /( dX * dX )

  Fdxdx = F1o12 /( dX * dX )

  Xdxdx = F1o180 /( dX * dX )

  Edxdx = F1o5040 /( dX * dX )

  fxx = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
!~~~~~~ fxx
        if(i+4 <= imax .and. i-4 >= imin)then

   fxx(i,j,k) = Edxdx*(-F9*fh(i-4,j,k)+F128*fh(i-3,j,k)-F1008*fh(i-2,j,k)+F8064*fh(i-1,j,k)-F14350*fh(i,j,k) &
                       -F9*fh(i+4,j,k)+F128*fh(i+3,j,k)-F1008*fh(i+2,j,k)+F8064*fh(i+1,j,k)               )

   elseif(i+3 <= imax .and. i-3 >= imin)then
   fxx(i,j,k) = Xdxdx*(TWO*fh(i-3,j,k)-F27*fh(i-2,j,k)+F270*fh(i-1,j,k)-F490*fh(i,j,k) &
                      +TWO*fh(i+3,j,k)-F27*fh(i+2,j,k)+F270*fh(i+1,j,k)               )
   elseif(i+2 <= imax .and. i-2 >= imin)then
   fxx(i,j,k) = Fdxdx*(-fh(i-2,j,k)+F16*fh(i-1,j,k)-F30*fh(i,j,k) &
                       -fh(i+2,j,k)+F16*fh(i+1,j,k)              )
   elseif(i+1 <= imax .and. i-1 >= imin)then
   fxx(i,j,k) = Sdxdx*(fh(i-1,j,k)-TWO*fh(i,j,k) &
                      +fh(i+1,j,k)              )
   endif

   enddo
   enddo
   enddo

  return

  end subroutine fddxx_sh

  subroutine fddyy_sh(ex,f,fyy,X,Y,Z,SYM1,SYM2,SYM3,symmetry,sst)
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry,sst
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fyy
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8,dimension(-3:ex(1)+4,-3:ex(2)+4,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdydy,Fdydy,Xdydy,Edydy
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1, F128=1.28d2
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2,F1008=1.008d3
  real*8, parameter :: F8064=8.064d3,F14350=1.435d4,THR=3.d0,F32=3.2d1,F168=1.68d2,F672=6.72d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3
  real*8, parameter :: F1o5040=ONE/5.04d3,F1o705600=ONE/7.056d5

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -3
     if(dabs(Y(1)) < dY) jmin = -3
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -3
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+4
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(4,ex,f,fh,SoA)

  Sdydy =  ONE /( dY * dY )

  Fdydy = F1o12 /( dY * dY )

  Xdydy = F1o180 /( dY * dY )
  
  Edydy = F1o5040 /( dY * dY )

  fyy = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
!~~~~~~ fyy
        if(j+4 <= jmax .and. j-4 >= jmin)then

   fyy(i,j,k) = Edydy*(-F9*fh(i,j-4,k)+F128*fh(i,j-3,k)-F1008*fh(i,j-2,k)+F8064*fh(i,j-1,k)-F14350*fh(i,j,k) &
                       -F9*fh(i,j+4,k)+F128*fh(i,j+3,k)-F1008*fh(i,j+2,k)+F8064*fh(i,j+1,k)               )

   elseif(j+3 <= jmax .and. j-3 >= jmin)then

   fyy(i,j,k) = Xdydy*(TWO*fh(i,j-3,k)-F27*fh(i,j-2,k)+F270*fh(i,j-1,k)-F490*fh(i,j,k) &
                      +TWO*fh(i,j+3,k)-F27*fh(i,j+2,k)+F270*fh(i,j+1,k)               )
   elseif(j+2 <= jmax .and. j-2 >= jmin)then

   fyy(i,j,k) = Fdydy*(-fh(i,j-2,k)+F16*fh(i,j-1,k)-F30*fh(i,j,k) &
                       -fh(i,j+2,k)+F16*fh(i,j+1,k)              )
   elseif(j+1 <= jmax .and. j-1 >= jmin)then

   fyy(i,j,k) = Sdydy*(fh(i,j-1,k)-TWO*fh(i,j,k) &
                      +fh(i,j+1,k)              )
   endif

   enddo
   enddo
   enddo

  return

  end subroutine fddyy_sh

  subroutine fddzz_sh(ex,f,fzz,X,Y,Z,SYM1,SYM2,SYM3,symmetry,sst)
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry,sst
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fzz
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8,dimension(-3:ex(1)+4,-3:ex(2)+4,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdzdz,Fdzdz,Xdzdz,Edzdz
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1, F128=1.28d2
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2,F1008=1.008d3
  real*8, parameter :: F8064=8.064d3,F14350=1.435d4,THR=3.d0,F32=3.2d1,F168=1.68d2,F672=6.72d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3
  real*8, parameter :: F1o5040=ONE/5.04d3,F1o705600=ONE/7.056d5

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -3
     if(dabs(Y(1)) < dY) jmin = -3
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -3
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+4
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(4,ex,f,fh,SoA)

  Sdzdz =  ONE /( dZ * dZ )

  Fdzdz = F1o12 /( dZ * dZ )

  Xdzdz = F1o180 /( dZ * dZ )

  Edzdz = F1o5040 /( dZ * dZ )
  
  fzz = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
!~~~~~~ fzz
        if(k+4 <= kmax .and. k-4 >= kmin)then

   fzz(i,j,k) = Edzdz*(-F9*fh(i,j,k-4)+F128*fh(i,j,k-3)-F1008*fh(i,j,k-2)+F8064*fh(i,j,k-1)-F14350*fh(i,j,k) &
                       -F9*fh(i,j,k+4)+F128*fh(i,j,k+3)-F1008*fh(i,j,k+2)+F8064*fh(i,j,k+1)               )

   elseif(k+3 <= kmax .and. k-3 >= kmin)then

   fzz(i,j,k) = Xdzdz*(TWO*fh(i,j,k-3)-F27*fh(i,j,k-2)+F270*fh(i,j,k-1)-F490*fh(i,j,k) &
                      +TWO*fh(i,j,k+3)-F27*fh(i,j,k+2)+F270*fh(i,j,k+1)               )
   elseif(k+2 <= kmax .and. k-2 >= kmin)then

   fzz(i,j,k) = Fdzdz*(-fh(i,j,k-2)+F16*fh(i,j,k-1)-F30*fh(i,j,k) &
                       -fh(i,j,k+2)+F16*fh(i,j,k+1)              )
   elseif(k+1 <= kmax .and. k-1 >= kmin)then

   fzz(i,j,k) = Sdzdz*(fh(i,j,k-1)-TWO*fh(i,j,k) &
                      +fh(i,j,k+1)              )
   endif

   enddo
   enddo
   enddo

  return

  end subroutine fddzz_sh

  subroutine fddxy_sh(ex,f,fxy,X,Y,Z,SYM1,SYM2,SYM3,symmetry,sst)
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry,sst
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fxy
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8,dimension(-3:ex(1)+4,-3:ex(2)+4,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdxdy,Fdxdy,Xdxdy,Edxdy
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1, F128=1.28d2
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2,F1008=1.008d3
  real*8, parameter :: F8064=8.064d3,F14350=1.435d4,THR=3.d0,F32=3.2d1,F168=1.68d2,F672=6.72d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3
  real*8, parameter :: F1o5040=ONE/5.04d3,F1o705600=ONE/7.056d5

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -3
     if(dabs(Y(1)) < dY) jmin = -3
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -3
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+4
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(4,ex,f,fh,SoA)

  Sdxdy = F1o4 /( dX * dY )

  Fdxdy = F1o144 /( dX * dY )

  Xdxdy = F1o3600 /( dX * dY )

  Edxdy = F1o705600 /( dX * dY )

  fxy = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
!~~~~~~ fxy
      if(i+4 <= imax .and. i-4 >= imin .and. j+4 <= jmax .and. j-4 >= jmin)then

   fxy(i,j,k) = Edxdy*( THR *( THR*fh(i-4,j-4,k)-F32*fh(i-3,j-4,k)+F168*fh(i-2,j-4,k)-F672*fh(i-1,j-4,k)        &
                              -THR*fh(i+4,j-4,k)+F32*fh(i+3,j-4,k)-F168*fh(i+2,j-4,k)+F672*fh(i+1,j-4,k))       &
                       -F32 *( THR*fh(i-4,j-3,k)-F32*fh(i-3,j-3,k)+F168*fh(i-2,j-3,k)-F672*fh(i-1,j-3,k)        &
                              -THR*fh(i+4,j-3,k)+F32*fh(i+3,j-3,k)-F168*fh(i+2,j-3,k)+F672*fh(i+1,j-3,k))       &
                       +F168*( THR*fh(i-4,j-2,k)-F32*fh(i-3,j-2,k)+F168*fh(i-2,j-2,k)-F672*fh(i-1,j-2,k)        &
                              -THR*fh(i+4,j-2,k)+F32*fh(i+3,j-2,k)-F168*fh(i+2,j-2,k)+F672*fh(i+1,j-2,k))       &
                       -F672*( THR*fh(i-4,j-1,k)-F32*fh(i-3,j-1,k)+F168*fh(i-2,j-1,k)-F672*fh(i-1,j-1,k)        &
                              -THR*fh(i+4,j-1,k)+F32*fh(i+3,j-1,k)-F168*fh(i+2,j-1,k)+F672*fh(i+1,j-1,k))       &
                       +F672*( THR*fh(i-4,j+1,k)-F32*fh(i-3,j+1,k)+F168*fh(i-2,j+1,k)-F672*fh(i-1,j+1,k)        &
                              -THR*fh(i+4,j+1,k)+F32*fh(i+3,j+1,k)-F168*fh(i+2,j+1,k)+F672*fh(i+1,j+1,k))       &
                       -F168*( THR*fh(i-4,j+2,k)-F32*fh(i-3,j+2,k)+F168*fh(i-2,j+2,k)-F672*fh(i-1,j+2,k)        &
                              -THR*fh(i+4,j+2,k)+F32*fh(i+3,j+2,k)-F168*fh(i+2,j+2,k)+F672*fh(i+1,j+2,k))       &
                       +F32 *( THR*fh(i-4,j+3,k)-F32*fh(i-3,j+3,k)+F168*fh(i-2,j+3,k)-F672*fh(i-1,j+3,k)        &
                              -THR*fh(i+4,j+3,k)+F32*fh(i+3,j+3,k)-F168*fh(i+2,j+3,k)+F672*fh(i+1,j+3,k))       &
                       -THR *( THR*fh(i-4,j+4,k)-F32*fh(i-3,j+4,k)+F168*fh(i-2,j+4,k)-F672*fh(i-1,j+4,k)        &
                              -THR*fh(i+4,j+4,k)+F32*fh(i+3,j+4,k)-F168*fh(i+2,j+4,k)+F672*fh(i+1,j+4,k)) )
   elseif(i+3 <= imax .and. i-3 >= imin .and. j+3 <= jmax .and. j-3 >= jmin)then

   fxy(i,j,k) = Xdxdy*(-    (-fh(i-3,j-3,k)+F9*fh(i-2,j-3,k)-F45*fh(i-1,j-3,k)+F45*fh(i+1,j-3,k)-F9*fh(i+2,j-3,k)+fh(i+3,j-3,k))  &
                       +F9 *(-fh(i-3,j-2,k)+F9*fh(i-2,j-2,k)-F45*fh(i-1,j-2,k)+F45*fh(i+1,j-2,k)-F9*fh(i+2,j-2,k)+fh(i+3,j-2,k))  &
                       -F45*(-fh(i-3,j-1,k)+F9*fh(i-2,j-1,k)-F45*fh(i-1,j-1,k)+F45*fh(i+1,j-1,k)-F9*fh(i+2,j-1,k)+fh(i+3,j-1,k))  &
                       +F45*(-fh(i-3,j+1,k)+F9*fh(i-2,j+1,k)-F45*fh(i-1,j+1,k)+F45*fh(i+1,j+1,k)-F9*fh(i+2,j+1,k)+fh(i+3,j+1,k))  &
                       -F9 *(-fh(i-3,j+2,k)+F9*fh(i-2,j+2,k)-F45*fh(i-1,j+2,k)+F45*fh(i+1,j+2,k)-F9*fh(i+2,j+2,k)+fh(i+3,j+2,k))  &
                       +    (-fh(i-3,j+3,k)+F9*fh(i-2,j+3,k)-F45*fh(i-1,j+3,k)+F45*fh(i+1,j+3,k)-F9*fh(i+2,j+3,k)+fh(i+3,j+3,k)))
   elseif(i+2 <= imax .and. i-2 >= imin .and. j+2 <= jmax .and. j-2 >= jmin)then

   fxy(i,j,k) = Fdxdy*(     (fh(i-2,j-2,k)-F8*fh(i-1,j-2,k)+F8*fh(i+1,j-2,k)-fh(i+2,j-2,k))  &
                       -F8 *(fh(i-2,j-1,k)-F8*fh(i-1,j-1,k)+F8*fh(i+1,j-1,k)-fh(i+2,j-1,k))  &
                       +F8 *(fh(i-2,j+1,k)-F8*fh(i-1,j+1,k)+F8*fh(i+1,j+1,k)-fh(i+2,j+1,k))  &
                       -    (fh(i-2,j+2,k)-F8*fh(i-1,j+2,k)+F8*fh(i+1,j+2,k)-fh(i+2,j+2,k)))
   elseif(i+1 <= imax .and. i-1 >= imin .and. j+1 <= jmax .and. j-1 >= jmin)then

   fxy(i,j,k) = Sdxdy*(fh(i-1,j-1,k)-fh(i+1,j-1,k)-fh(i-1,j+1,k)+fh(i+1,j+1,k))
   endif

   enddo
   enddo
   enddo

  return

  end subroutine fddxy_sh

  subroutine fddxz_sh(ex,f,fxz,X,Y,Z,SYM1,SYM2,SYM3,symmetry,sst)
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry,sst
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fxz
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8,dimension(-3:ex(1)+4,-3:ex(2)+4,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdxdz,Fdxdz,Xdxdz,Edxdz
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1, F128=1.28d2
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2,F1008=1.008d3
  real*8, parameter :: F8064=8.064d3,F14350=1.435d4,THR=3.d0,F32=3.2d1,F168=1.68d2,F672=6.72d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3
  real*8, parameter :: F1o5040=ONE/5.04d3,F1o705600=ONE/7.056d5

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -3
     if(dabs(Y(1)) < dY) jmin = -3
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -3
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+4
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(4,ex,f,fh,SoA)

  Sdxdz = F1o4 /( dX * dZ )

  Fdxdz = F1o144 /( dX * dZ )

  Xdxdz = F1o3600 /( dX * dZ )

  Edxdz = F1o705600 /( dX * dZ )

  fxz = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
!~~~~~~ fxz
      if(i+4 <= imax .and. i-4 >= imin .and. k+4 <= kmax .and. k-4 >= kmin)then

   fxz(i,j,k) = Edxdz*( THR *( THR*fh(i-4,j,k-4)-F32*fh(i-3,j,k-4)+F168*fh(i-2,j,k-4)-F672*fh(i-1,j,k-4)        &
                              -THR*fh(i+4,j,k-4)+F32*fh(i+3,j,k-4)-F168*fh(i+2,j,k-4)+F672*fh(i+1,j,k-4))       &
                       -F32 *( THR*fh(i-4,j,k-3)-F32*fh(i-3,j,k-3)+F168*fh(i-2,j,k-3)-F672*fh(i-1,j,k-3)        &
                              -THR*fh(i+4,j,k-3)+F32*fh(i+3,j,k-3)-F168*fh(i+2,j,k-3)+F672*fh(i+1,j,k-3))       &
                       +F168*( THR*fh(i-4,j,k-2)-F32*fh(i-3,j,k-2)+F168*fh(i-2,j,k-2)-F672*fh(i-1,j,k-2)        &
                              -THR*fh(i+4,j,k-2)+F32*fh(i+3,j,k-2)-F168*fh(i+2,j,k-2)+F672*fh(i+1,j,k-2))       &
                       -F672*( THR*fh(i-4,j,k-1)-F32*fh(i-3,j,k-1)+F168*fh(i-2,j,k-1)-F672*fh(i-1,j,k-1)        &
                              -THR*fh(i+4,j,k-1)+F32*fh(i+3,j,k-1)-F168*fh(i+2,j,k-1)+F672*fh(i+1,j,k-1))       &
                       +F672*( THR*fh(i-4,j,k+1)-F32*fh(i-3,j,k+1)+F168*fh(i-2,j,k+1)-F672*fh(i-1,j,k+1)        &
                              -THR*fh(i+4,j,k+1)+F32*fh(i+3,j,k+1)-F168*fh(i+2,j,k+1)+F672*fh(i+1,j,k+1))       &
                       -F168*( THR*fh(i-4,j,k+2)-F32*fh(i-3,j,k+2)+F168*fh(i-2,j,k+2)-F672*fh(i-1,j,k+2)        &
                              -THR*fh(i+4,j,k+2)+F32*fh(i+3,j,k+2)-F168*fh(i+2,j,k+2)+F672*fh(i+1,j,k+2))       &
                       +F32 *( THR*fh(i-4,j,k+3)-F32*fh(i-3,j,k+3)+F168*fh(i-2,j,k+3)-F672*fh(i-1,j,k+3)        &
                              -THR*fh(i+4,j,k+3)+F32*fh(i+3,j,k+3)-F168*fh(i+2,j,k+3)+F672*fh(i+1,j,k+3))       &
                       -THR *( THR*fh(i-4,j,k+4)-F32*fh(i-3,j,k+4)+F168*fh(i-2,j,k+4)-F672*fh(i-1,j,k+4)        &
                              -THR*fh(i+4,j,k+4)+F32*fh(i+3,j,k+4)-F168*fh(i+2,j,k+4)+F672*fh(i+1,j,k+4)) )
   elseif(i+3 <= imax .and. i-3 >= imin .and. k+3 <= kmax .and. k-3 >= kmin)then

   fxz(i,j,k) = Xdxdz*(-    (-fh(i-3,j,k-3)+F9*fh(i-2,j,k-3)-F45*fh(i-1,j,k-3)+F45*fh(i+1,j,k-3)-F9*fh(i+2,j,k-3)+fh(i+3,j,k-3))  &
                       +F9 *(-fh(i-3,j,k-2)+F9*fh(i-2,j,k-2)-F45*fh(i-1,j,k-2)+F45*fh(i+1,j,k-2)-F9*fh(i+2,j,k-2)+fh(i+3,j,k-2))  &
                       -F45*(-fh(i-3,j,k-1)+F9*fh(i-2,j,k-1)-F45*fh(i-1,j,k-1)+F45*fh(i+1,j,k-1)-F9*fh(i+2,j,k-1)+fh(i+3,j,k-1))  &
                       +F45*(-fh(i-3,j,k+1)+F9*fh(i-2,j,k+1)-F45*fh(i-1,j,k+1)+F45*fh(i+1,j,k+1)-F9*fh(i+2,j,k+1)+fh(i+3,j,k+1))  &
                       -F9 *(-fh(i-3,j,k+2)+F9*fh(i-2,j,k+2)-F45*fh(i-1,j,k+2)+F45*fh(i+1,j,k+2)-F9*fh(i+2,j,k+2)+fh(i+3,j,k+2))  &
                       +    (-fh(i-3,j,k+3)+F9*fh(i-2,j,k+3)-F45*fh(i-1,j,k+3)+F45*fh(i+1,j,k+3)-F9*fh(i+2,j,k+3)+fh(i+3,j,k+3)))
   elseif(i+2 <= imax .and. i-2 >= imin .and. k+2 <= kmax .and. k-2 >= kmin)then
   fxz(i,j,k) = Fdxdz*(     (fh(i-2,j,k-2)-F8*fh(i-1,j,k-2)+F8*fh(i+1,j,k-2)-fh(i+2,j,k-2))  &
                       -F8 *(fh(i-2,j,k-1)-F8*fh(i-1,j,k-1)+F8*fh(i+1,j,k-1)-fh(i+2,j,k-1))  &
                       +F8 *(fh(i-2,j,k+1)-F8*fh(i-1,j,k+1)+F8*fh(i+1,j,k+1)-fh(i+2,j,k+1))  &
                       -    (fh(i-2,j,k+2)-F8*fh(i-1,j,k+2)+F8*fh(i+1,j,k+2)-fh(i+2,j,k+2)))
   elseif(i+1 <= imax .and. i-1 >= imin .and. k+1 <= kmax .and. k-1 >= kmin)then
   fxz(i,j,k) = Sdxdz*(fh(i-1,j,k-1)-fh(i+1,j,k-1)-fh(i-1,j,k+1)+fh(i+1,j,k+1))
   endif

   enddo
   enddo
   enddo

  return

  end subroutine fddxz_sh

  subroutine fddyz_sh(ex,f,fyz,X,Y,Z,SYM1,SYM2,SYM3,symmetry,sst)
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry,sst
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fyz
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8,dimension(-3:ex(1)+4,-3:ex(2)+4,ex(3))   :: fh
  real*8, dimension(2) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdydz,Fdydz,Xdydz,Edydz
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1, F128=1.28d2
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2,F1008=1.008d3
  real*8, parameter :: F8064=8.064d3,F14350=1.435d4,THR=3.d0,F32=3.2d1,F168=1.68d2,F672=6.72d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3
  real*8, parameter :: F1o5040=ONE/5.04d3,F1o705600=ONE/7.056d5

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry == OCTANT)then
     if(dabs(X(1)) < dX) imin = -3
     if(dabs(Y(1)) < dY) jmin = -3
  elseif(Symmetry == EQ_SYMM)then
     if((sst==2.or.sst==4).and.dabs(Y(1)) < dY) jmin = -3
     if((sst==3.or.sst==5).and.dabs(Y(ex(2))) < dY) jmax=ex(2)+4
  endif

  if(sst==0)then
     SoA(1) = SYM1
     SoA(2) = SYM2
  elseif(sst==2.or.sst==3)then
     SoA(1) = SYM2
     SoA(2) = SYM3
  elseif(sst==4.or.sst==5)then
     SoA(1) = SYM1
     SoA(2) = SYM3
  endif

  call symmetry_stbd(4,ex,f,fh,SoA)

  Sdydz = F1o4 /( dY * dZ )

  Fdydz = F1o144 /( dY * dZ )

  Xdydz = F1o3600 /( dY * dZ )

  Edydz = F1o705600 /( dY * dZ )

  fyz = ZEO

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
!~~~~~~ fyz
      if(j+4 <= jmax .and. j-4 >= jmin .and. k+4 <= kmax .and. k-4 >= kmin)then

   fyz(i,j,k) = Edydz*( THR *( THR*fh(i,j-4,k-4)-F32*fh(i,j-3,k-4)+F168*fh(i,j-2,k-4)-F672*fh(i,j-1,k-4)        &
                              -THR*fh(i,j+4,k-4)+F32*fh(i,j+3,k-4)-F168*fh(i,j+2,k-4)+F672*fh(i,j+1,k-4))       &
                       -F32 *( THR*fh(i,j-4,k-3)-F32*fh(i,j-3,k-3)+F168*fh(i,j-2,k-3)-F672*fh(i,j-1,k-3)        &
                              -THR*fh(i,j+4,k-3)+F32*fh(i,j+3,k-3)-F168*fh(i,j+2,k-3)+F672*fh(i,j+1,k-3))       &
                       +F168*( THR*fh(i,j-4,k-2)-F32*fh(i,j-3,k-2)+F168*fh(i,j-2,k-2)-F672*fh(i,j-1,k-2)        &
                              -THR*fh(i,j+4,k-2)+F32*fh(i,j+3,k-2)-F168*fh(i,j+2,k-2)+F672*fh(i,j+1,k-2))       &
                       -F672*( THR*fh(i,j-4,k-1)-F32*fh(i,j-3,k-1)+F168*fh(i,j-2,k-1)-F672*fh(i,j-1,k-1)        &
                              -THR*fh(i,j+4,k-1)+F32*fh(i,j+3,k-1)-F168*fh(i,j+2,k-1)+F672*fh(i,j+1,k-1))       &
                       +F672*( THR*fh(i,j-4,k+1)-F32*fh(i,j-3,k+1)+F168*fh(i,j-2,k+1)-F672*fh(i,j-1,k+1)        &
                              -THR*fh(i,j+4,k+1)+F32*fh(i,j+3,k+1)-F168*fh(i,j+2,k+1)+F672*fh(i,j+1,k+1))       &
                       -F168*( THR*fh(i,j-4,k+2)-F32*fh(i,j-3,k+2)+F168*fh(i,j-2,k+2)-F672*fh(i,j-1,k+2)        &
                              -THR*fh(i,j+4,k+2)+F32*fh(i,j+3,k+2)-F168*fh(i,j+2,k+2)+F672*fh(i,j+1,k+2))       &
                       +F32 *( THR*fh(i,j-4,k+3)-F32*fh(i,j-3,k+3)+F168*fh(i,j-2,k+3)-F672*fh(i,j-1,k+3)        &
                              -THR*fh(i,j+4,k+3)+F32*fh(i,j+3,k+3)-F168*fh(i,j+2,k+3)+F672*fh(i,j+1,k+3))       &
                       -THR *( THR*fh(i,j-4,k+4)-F32*fh(i,j-3,k+4)+F168*fh(i,j-2,k+4)-F672*fh(i,j-1,k+4)        &
                              -THR*fh(i,j+4,k+4)+F32*fh(i,j+3,k+4)-F168*fh(i,j+2,k+4)+F672*fh(i,j+1,k+4)) )
   elseif(j+3 <= jmax .and. j-3 >= jmin .and. k+3 <= kmax .and. k-3 >= kmin)then

   fyz(i,j,k) = Xdydz*(-    (-fh(i,j-3,k-3)+F9*fh(i,j-2,k-3)-F45*fh(i,j-1,k-3)+F45*fh(i,j+1,k-3)-F9*fh(i,j+2,k-3)+fh(i,j+3,k-3))  &
                       +F9 *(-fh(i,j-3,k-2)+F9*fh(i,j-2,k-2)-F45*fh(i,j-1,k-2)+F45*fh(i,j+1,k-2)-F9*fh(i,j+2,k-2)+fh(i,j+3,k-2))  &
                       -F45*(-fh(i,j-3,k-1)+F9*fh(i,j-2,k-1)-F45*fh(i,j-1,k-1)+F45*fh(i,j+1,k-1)-F9*fh(i,j+2,k-1)+fh(i,j+3,k-1))  &
                       +F45*(-fh(i,j-3,k+1)+F9*fh(i,j-2,k+1)-F45*fh(i,j-1,k+1)+F45*fh(i,j+1,k+1)-F9*fh(i,j+2,k+1)+fh(i,j+3,k+1))  &
                       -F9 *(-fh(i,j-3,k+2)+F9*fh(i,j-2,k+2)-F45*fh(i,j-1,k+2)+F45*fh(i,j+1,k+2)-F9*fh(i,j+2,k+2)+fh(i,j+3,k+2))  &
                       +    (-fh(i,j-3,k+3)+F9*fh(i,j-2,k+3)-F45*fh(i,j-1,k+3)+F45*fh(i,j+1,k+3)-F9*fh(i,j+2,k+3)+fh(i,j+3,k+3)))
   elseif(j+2 <= jmax .and. j-2 >= jmin .and. k+2 <= kmax .and. k-2 >= kmin)then
   fyz(i,j,k) = Fdydz*(     (fh(i,j-2,k-2)-F8*fh(i,j-1,k-2)+F8*fh(i,j+1,k-2)-fh(i,j+2,k-2))  &
                       -F8 *(fh(i,j-2,k-1)-F8*fh(i,j-1,k-1)+F8*fh(i,j+1,k-1)-fh(i,j+2,k-1))  &
                       +F8 *(fh(i,j-2,k+1)-F8*fh(i,j-1,k+1)+F8*fh(i,j+1,k+1)-fh(i,j+2,k+1))  &
                       -    (fh(i,j-2,k+2)-F8*fh(i,j-1,k+2)+F8*fh(i,j+1,k+2)-fh(i,j+2,k+2)))
   elseif(j+1 <= jmax .and. j-1 >= jmin .and. k+1 <= kmax .and. k-1 >= kmin)then
   fyz(i,j,k) = Sdydz*(fh(i,j-1,k-1)-fh(i,j+1,k-1)-fh(i,j-1,k+1)+fh(i,j+1,k+1))
   endif 

   enddo
   enddo
   enddo

  return

  end subroutine fddyz_sh

#endif 

!common code for different finite difference order
subroutine fderivs_shc(ex,f,fx,fy,fz,crho,sigma,R,SYM1,SYM2,SYM3,Symmetry,Lev,sst, &
                       drhodx, drhody, drhodz,                                     &
                       dsigmadx,dsigmady,dsigmadz,                                 &
                       dRdx,dRdy,dRdz)

  implicit none
  integer,intent(in ):: ex(1:3), Symmetry,Lev,sst
  double precision,intent(in),dimension(ex(1))::crho
  double precision,intent(in),dimension(ex(2))::sigma
  double precision,intent(in),dimension(ex(3))::R
  real*8,intent(in ):: SYM1,SYM2,SYM3
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::f
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::drhodx, drhody, drhodz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::dsigmadx,dsigmady,dsigmadz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::dRdx,dRdy,dRdz
  double precision,intent(out),dimension(ex(1),ex(2),ex(3))::fx,fy,fz

#if 0  
  integer :: i,j,k
  
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
  call point_fderivs_shc(ex,f,fx(i,j,k),fy(i,j,k),fz(i,j,k),crho,sigma,R,SYM1,SYM2,SYM3,Symmetry,Lev,sst, &
                       drhodx, drhody, drhodz,                                     &
                       dsigmadx,dsigmady,dsigmadz,                                 &
                       dRdx,dRdy,dRdz,i,j,k)
 enddo
 enddo
 enddo
#else
  double precision,dimension(ex(1),ex(2),ex(3))::gx,gy,gz

  call fderivs_sh(ex,f,gx,gy,gz,crho,sigma,R,SYM1, SYM2,SYM3,Symmetry,Lev,sst)

   fx = dRdx*gz+drhodx*gx+dsigmadx*gy
   fy = dRdy*gz+drhody*gx+dsigmady*gy
   fz = dRdz*gz+drhodz*gx+dsigmadz*gy
#endif

   return

end subroutine fderivs_shc

subroutine fdderivs_shc(ex,f,fxx,fxy,fxz,fyy,fyz,fzz,crho,sigma,R,SYM1,SYM2,SYM3,Symmetry,Lev,sst,&
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz)

  implicit none
  integer,intent(in ):: ex(1:3), Symmetry,Lev,sst
  double precision,intent(in),dimension(ex(1))::crho
  double precision,intent(in),dimension(ex(2))::sigma
  double precision,intent(in),dimension(ex(3))::R
  real*8,intent(in ):: SYM1,SYM2,SYM3
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::f
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::drhodx, drhody, drhodz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::dsigmadx,dsigmady,dsigmadz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::dRdx,dRdy,dRdz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz
  double precision,intent(out),dimension(ex(1),ex(2),ex(3))::fxx,fxy,fxz,fyy,fyz,fzz

#if 0  
  integer :: i,j,k
  
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
  call point_fdderivs_shc(ex,f,fxx(i,j,k),fxy(i,j,k),fxz(i,j,k),fyy(i,j,k),fyz(i,j,k),fzz(i,j,k),crho,sigma,R,SYM1,SYM2,SYM3,Symmetry,Lev,sst,&
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,i,j,k)
 enddo
 enddo
 enddo
#else
  double precision,dimension(ex(1),ex(2),ex(3))::gx,gy,gz,gxx,gxy,gxz,gyy,gyz,gzz
  real*8,parameter :: TWO = 2.d0

  call fderivs_sh(ex,f,gx,gy,gz,crho,sigma,R,SYM1, SYM2,SYM3,Symmetry,Lev,sst)
  call fdderivs_sh(ex,f,gxx,gxy,gxz,gyy,gyz,gzz,crho,sigma,R,SYM1,SYM2,SYM3,Symmetry,Lev,sst)

   fxx = dRdxx*gz+drhodxx*gx+dsigmadxx*gy + &
         dRdx*dRdx*gzz+drhodx*drhodx*gxx+dsigmadx*dsigmadx*gyy + &
         TWO*(dRdx*drhodx*gxz+dRdx*dsigmadx*gyz+drhodx*dsigmadx*gxy)
   fyy = dRdyy*gz+drhodyy*gx+dsigmadyy*gy + &
         dRdy*dRdy*gzz+drhody*drhody*gxx+dsigmady*dsigmady*gyy + &
         TWO*(dRdy*drhody*gxz+dRdy*dsigmady*gyz+drhody*dsigmady*gxy)
   fzz = dRdzz*gz+drhodzz*gx+dsigmadzz*gy + &
         dRdz*dRdz*gzz+drhodz*drhodz*gxx+dsigmadz*dsigmadz*gyy + &
         TWO*(dRdz*drhodz*gxz+dRdz*dsigmadz*gyz+drhodz*dsigmadz*gxy)
   fxy = dRdxy*gz+drhodxy*gx+dsigmadxy*gy + &
         dRdx*drhody*gxz+dRdx*dsigmady*gyz+drhodx*dsigmady*gxy + &
         dRdy*drhodx*gxz+dRdy*dsigmadx*gyz+drhody*dsigmadx*gxy + &
         dRdx*dRdy*gzz+drhodx*drhody*gxx+dsigmadx*dsigmady*gyy
   fxz = dRdxz*gz+drhodxz*gx+dsigmadxz*gy + &
         dRdx*drhodz*gxz+dRdx*dsigmadz*gyz+drhodx*dsigmadz*gxy + &
         dRdz*drhodx*gxz+dRdz*dsigmadx*gyz+drhodz*dsigmadx*gxy + &
         dRdx*dRdz*gzz+drhodx*drhodz*gxx+dsigmadx*dsigmadz*gyy
   fyz = dRdyz*gz+drhodyz*gx+dsigmadyz*gy + &
         dRdz*drhody*gxz+dRdz*dsigmady*gyz+drhodz*dsigmady*gxy + &
         dRdy*drhodz*gxz+dRdy*dsigmadz*gyz+drhody*dsigmadz*gxy + &
         dRdz*dRdy*gzz+drhodz*drhody*gxx+dsigmadz*dsigmady*gyy
#endif

   return

end subroutine fdderivs_shc
