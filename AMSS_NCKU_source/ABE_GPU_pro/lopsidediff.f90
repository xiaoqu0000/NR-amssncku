
! Compute advection terms in right hand sides of field equations

#include "macrodef.fh"

! we need only distinguish different finite difference order
! Vertex or Cell is distinguished in routine symmetry_bd which locates in
! file "fmisc.f90"

#if (ghost_width == 2)
! second order code

!-----------------------------------------------------------------------------
!         v
! D f = ------[ - 3 f  + 4 f   - f     ]
!  i     2dx         i      i+v   i+2v
!
! where
!
!        i
!      |B |
! v = -----
!        i
!       B
!
!-----------------------------------------------------------------------------
subroutine lopsided(ex,X,Y,Z,f,f_rhs,Sfx,Sfy,Sfz,Symmetry,SoA)
  implicit none

!~~~~~~> Input parameters:

  integer, intent(in)  :: ex(1:3),Symmetry
  real*8,  intent(in)  :: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8,dimension(ex(1),ex(2),ex(3)),intent(in)   :: f,Sfx,Sfy,Sfz

  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout):: f_rhs
  real*8,dimension(3),intent(in) ::SoA

!~~~~~~> local variables:
! note index -1,0, so we have 2 extra points
  real*8,dimension(-1:ex(1),-1:ex(2),-1:ex(3))   :: fh
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8 :: dX,dY,dZ
  real*8 :: d2dx,d2dy,d2dz
  real*8,  parameter :: ZEO=0.d0,ONE=1.d0,TWO=2.d0,THR=3.d0,FOUR=4.d0
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  d2dx = ONE/TWO/dX
  d2dy = ONE/TWO/dY
  d2dz = ONE/TWO/dZ

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1
  if(Symmetry > NO_SYMM .and. dabs(Z(1)) < dZ) kmin = -1
  if(Symmetry > EQ_SYMM .and. dabs(X(1)) < dX) imin = -1
  if(Symmetry > EQ_SYMM .and. dabs(Y(1)) < dY) jmin = -1

  call symmetry_bd(2,ex,f,fh,SoA)

! upper bound set ex-1 only for efficiency, 
! the loop body will set ex 0 also
  do k=1,ex(3)-1
  do j=1,ex(2)-1
  do i=1,ex(1)-1
! x direction   
    if(Sfx(i,j,k) >= ZEO)then
       if( i+2 <= imax .and. i >= imin)then
!         v
! D f = ------[ - 3 f  + 4 f   - f     ]
!  i     2dx         i      i+v   i+2v
     f_rhs(i,j,k)=f_rhs(i,j,k)+                           &
                  Sfx(i,j,k)*d2dx*(-THR*fh(i,j,k)+FOUR*fh(i+1,j,k)-fh(i+2,j,k))
       elseif(i+1 <= imax .and. i >= imin)then
!         v
! D f = ------[ - f  + f   ]
!  i      dx       i    i+v
     f_rhs(i,j,k)=f_rhs(i,j,k)+                           &
                  Sfx(i,j,k)*d2dx*(-fh(i,j,k)+fh(i+1,j,k))

       endif

    elseif(Sfx(i,j,k) <= ZEO)then
      if( i-2 >= imin .and. i <= imax)then
     f_rhs(i,j,k)=f_rhs(i,j,k)-                           &
                  Sfx(i,j,k)*d2dx*(-THR*fh(i,j,k)+FOUR*fh(i-1,j,k)-fh(i-2,j,k))
      elseif(i-1 >= imin .and. i <= imax)then
     f_rhs(i,j,k)=f_rhs(i,j,k)-                           &
                  Sfx(i,j,k)*d2dx*(-fh(i,j,k)+fh(i-1,j,k))
      endif

! set imax and imin 0
    endif

! y direction   
    if(Sfy(i,j,k) >= ZEO)then
       if( j+2 <= jmax .and. j >= jmin)then
!         v
! D f = ------[ - 3 f  + 4 f   - f     ]
!  i     2dx         i      i+v   i+2v
     f_rhs(i,j,k)=f_rhs(i,j,k)+                           &
                  Sfy(i,j,k)*d2dy*(-THR*fh(i,j,k)+FOUR*fh(i,j+1,k)-fh(i,j+2,k))
       elseif(j+1 <= jmax .and. j >= jmin)then
!         v
! D f = ------[ - f  + f   ]
!  i      dx       i    i+v
     f_rhs(i,j,k)=f_rhs(i,j,k)+                           &
                  Sfy(i,j,k)*d2dy*(-fh(i,j,k)+fh(i,j+1,k))
       endif

    elseif(Sfy(i,j,k) <= ZEO)then
      if( j-2 >= jmin .and. j <= jmax)then
     f_rhs(i,j,k)=f_rhs(i,j,k)-                           &
                  Sfy(i,j,k)*d2dy*(-THR*fh(i,j,k)+FOUR*fh(i,j-1,k)-fh(i,j-2,k))
      elseif(j-1 >= jmin .and. j <= jmax)then
     f_rhs(i,j,k)=f_rhs(i,j,k)-                           &
                  Sfy(i,j,k)*d2dy*(-fh(i,j,k)+fh(i,j-1,k))
      endif

! set jmin and jmax 0
     endif
!! z direction   
    if(Sfz(i,j,k) >= ZEO)then
      if( k+2 <= kmax .and. k >= kmin)then
!         v
! D f = ------[ - 3 f  + 4 f   - f     ]
!  i     2dx         i      i+v   i+2v
     f_rhs(i,j,k)=f_rhs(i,j,k)+                           &
                  Sfz(i,j,k)*d2dz*(-THR*fh(i,j,k)+FOUR*fh(i,j,k+1)-fh(i,j,k+2))
       elseif(k+1 <= kmax .and. k >= kmin)then
!         v
! D f = ------[ - f  + f   ]
!  i      dx       i    i+v
     f_rhs(i,j,k)=f_rhs(i,j,k)+                           &
                  Sfz(i,j,k)*d2dz*(-fh(i,j,k)+fh(i,j,k+1))
       endif

    elseif(Sfz(i,j,k) <= ZEO)then
      if( k-2 >= kmin .and. k <= kmax)then
     f_rhs(i,j,k)=f_rhs(i,j,k)-                           &
                  Sfz(i,j,k)*d2dz*(-THR*fh(i,j,k)+FOUR*fh(i,j,k-1)-fh(i,j,k-2))
      elseif(k-1 >= kmin .and. k <= kmax)then
     f_rhs(i,j,k)=f_rhs(i,j,k)-                           &
                  Sfz(i,j,k)*d2dz*(-fh(i,j,k)+fh(i,j,k-1))
      endif

! set kmin and kmax 0
     endif

  enddo
  enddo
  enddo

  return

  end subroutine lopsided

#elif (ghost_width == 3)
! fourth order code

!-----------------------------------------------------------------------------
!
! Compute advection terms in right hand sides of field equations
!         v
! D f = ------[ - 3f    - 10f  + 18f    - 6f     + f     ]
!  i     12dx       i-v      i      i+v     i+2v    i+3v
!
! where
!
!        i
!      |B |
! v = -----
!        i
!       B
!
!-----------------------------------------------------------------------------

subroutine lopsided(ex,X,Y,Z,f,f_rhs,Sfx,Sfy,Sfz,Symmetry,SoA)
  implicit none

!~~~~~~> Input parameters:

  integer, intent(in)  :: ex(1:3),Symmetry
  real*8,  intent(in)  :: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8,dimension(ex(1),ex(2),ex(3)),intent(in)   :: f,Sfx,Sfy,Sfz

  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout):: f_rhs
  real*8,dimension(3),intent(in) ::SoA

!~~~~~~> local variables:
! note index -2,-1,0, so we have 3 extra points
  real*8,dimension(-2:ex(1),-2:ex(2),-2:ex(3))   :: fh
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8 :: dX,dY,dZ
  real*8 :: d12dx,d12dy,d12dz,d2dx,d2dy,d2dz
  real*8,  parameter :: ZEO=0.d0,ONE=1.d0, F3=3.d0
  real*8,  parameter :: TWO=2.d0,F6=6.0d0,F18=1.8d1
  real*8,  parameter :: F12=1.2d1, F10=1.d1,EIT=8.d0
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  d12dx = ONE/F12/dX
  d12dy = ONE/F12/dY
  d12dz = ONE/F12/dZ

  d2dx = ONE/TWO/dX
  d2dy = ONE/TWO/dY
  d2dz = ONE/TWO/dZ

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1
  if(Symmetry > NO_SYMM .and. dabs(Z(1)) < dZ) kmin = -2
  if(Symmetry > EQ_SYMM .and. dabs(X(1)) < dX) imin = -2
  if(Symmetry > EQ_SYMM .and. dabs(Y(1)) < dY) jmin = -2

  call symmetry_bd(3,ex,f,fh,SoA)

! upper bound set ex-1 only for efficiency, 
! the loop body will set ex 0 also
  do k=1,ex(3)-1
  do j=1,ex(2)-1
  do i=1,ex(1)-1
#if 0  
!! old code
! x direction   
    if(Sfx(i,j,k) >= ZEO .and. i+3 <= imax .and. i-1 >= imin)then
!         v
! D f = ------[ - 3f    - 10f  + 18f    - 6f     + f     ]
!  i     12dx       i-v      i      i+v     i+2v    i+3v
     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                   &
                  Sfx(i,j,k)*d12dx*(-F3*fh(i-1,j,k)-F10*fh(i,j,k)+F18*fh(i+1,j,k) &
                                    -F6*fh(i+2,j,k)+    fh(i+3,j,k))

    elseif(Sfx(i,j,k) <= ZEO .and. i-3 >= imin .and. i+1 <= imax)then
     f_rhs(i,j,k)=f_rhs(i,j,k)-                                                   &
                  Sfx(i,j,k)*d12dx*(-F3*fh(i+1,j,k)-F10*fh(i,j,k)+F18*fh(i-1,j,k) &
                                    -F6*fh(i-2,j,k)+    fh(i-3,j,k))

     elseif(i+2 <= imax .and. i-2 >= imin)then
!
!              f(i-2) - 8 f(i-1) + 8 f(i+1) - f(i+2)
!  fx(i) = ---------------------------------------------
!                             12 dx
     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                           &
                  Sfx(i,j,k)*d12dx*(fh(i-2,j,k)-EIT*fh(i-1,j,k)+EIT*fh(i+1,j,k)-fh(i+2,j,k))

     elseif(i+1 <= imax .and. i-1 >= imin)then
!
!              - f(i-1) + f(i+1)
!  fx(i) = --------------------------------
!                     2 dx
     f_rhs(i,j,k)=f_rhs(i,j,k) + Sfx(i,j,k)*d2dx*(-fh(i-1,j,k)+fh(i+1,j,k))

! set imax and imin 0
    endif

! y direction   
    if(Sfy(i,j,k) >= ZEO .and. j+3 <= jmax .and. j-1 >= jmin)then
!         v
! D f = ------[ - 3f    - 10f  + 18f    - 6f     + f     ]
!  i     12dx       i-v      i      i+v     i+2v    i+3v
     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                   &
                  Sfy(i,j,k)*d12dy*(-F3*fh(i,j-1,k)-F10*fh(i,j,k)+F18*fh(i,j+1,k) &
                                    -F6*fh(i,j+2,k)+    fh(i,j+3,k))

    elseif(Sfy(i,j,k) <= ZEO .and. j-3 >= jmin .and. j+1 <= jmax)then
     f_rhs(i,j,k)=f_rhs(i,j,k)-                                                   &
                  Sfy(i,j,k)*d12dy*(-F3*fh(i,j+1,k)-F10*fh(i,j,k)+F18*fh(i,j-1,k) &
                                    -F6*fh(i,j-2,k)+    fh(i,j-3,k))

     elseif(j+2 <= jmax .and. j-2 >= jmin)then

     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                            &
                  Sfy(i,j,k)*d12dy*(fh(i,j-2,k)-EIT*fh(i,j-1,k)+EIT*fh(i,j+1,k)-fh(i,j+2,k))

     elseif(j+1 <= jmax .and. j-1 >= jmin)then

     f_rhs(i,j,k)=f_rhs(i,j,k) + Sfy(i,j,k)*d2dy*(-fh(i,j-1,k)+fh(i,j+1,k))
! set jmin and jmax 0
     endif
!! z direction   
    if(Sfz(i,j,k) >= ZEO .and. k+3 <= kmax .and. k-1 >= kmin)then
!         v
! D f = ------[ - 3f    - 10f  + 18f    - 6f     + f     ]
!  i     12dx       i-v      i      i+v     i+2v    i+3v
     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                   &
                  Sfz(i,j,k)*d12dz*(-F3*fh(i,j,k-1)-F10*fh(i,j,k)+F18*fh(i,j,k+1) &
                                    -F6*fh(i,j,k+2)+    fh(i,j,k+3))

    elseif(Sfz(i,j,k) <= ZEO .and. k-3 >= kmin .and. k+1 <= kmax)then
     f_rhs(i,j,k)=f_rhs(i,j,k)-                                                   &
                  Sfz(i,j,k)*d12dz*(-F3*fh(i,j,k+1)-F10*fh(i,j,k)+F18*fh(i,j,k-1) &
                                    -F6*fh(i,j,k-2)+    fh(i,j,k-3))

     elseif(k+2 <= kmax .and. k-2 >= kmin)then

     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                            &
                  Sfz(i,j,k)*d12dz*(fh(i,j,k-2)-EIT*fh(i,j,k-1)+EIT*fh(i,j,k+1)-fh(i,j,k+2))

     elseif(k+1 <= kmax .and. k-1 >= kmin)then

     f_rhs(i,j,k)=f_rhs(i,j,k)+Sfz(i,j,k)*d2dz*(-fh(i,j,k-1)+fh(i,j,k+1))
! set kmin and kmax 0
     endif
#else
!! new code, 2012dec27, based on bam
! x direction   
    if(Sfx(i,j,k) > ZEO)then
      if(i+3 <= imax)then
!         v
! D f = ------[ - 3f    - 10f  + 18f    - 6f     + f     ]
!  i     12dx       i-v      i      i+v     i+2v    i+3v
     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                   &
                  Sfx(i,j,k)*d12dx*(-F3*fh(i-1,j,k)-F10*fh(i,j,k)+F18*fh(i+1,j,k) &
                                    -F6*fh(i+2,j,k)+    fh(i+3,j,k))
     elseif(i+2 <= imax)then
!
!              f(i-2) - 8 f(i-1) + 8 f(i+1) - f(i+2)
!  fx(i) = ---------------------------------------------
!                             12 dx
     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                           &
                  Sfx(i,j,k)*d12dx*(fh(i-2,j,k)-EIT*fh(i-1,j,k)+EIT*fh(i+1,j,k)-fh(i+2,j,k))

     elseif(i+1 <= imax)then
!         v
! D f = ------[   3f    + 10f  - 18f    + 6f     - f     ]
!  i     12dx       i+v      i      i-v     i-2v    i-3v
     f_rhs(i,j,k)=f_rhs(i,j,k)-                                                   &
                  Sfx(i,j,k)*d12dx*(-F3*fh(i+1,j,k)-F10*fh(i,j,k)+F18*fh(i-1,j,k) &
                                    -F6*fh(i-2,j,k)+    fh(i-3,j,k))
! set imax and imin 0
     endif
   elseif(Sfx(i,j,k) < ZEO)then
      if(i-3 >= imin)then
!         v
! D f = ------[ - 3f    - 10f  + 18f    - 6f     + f     ]
!  i     12dx       i-v      i      i+v     i+2v    i+3v
     f_rhs(i,j,k)=f_rhs(i,j,k)-                                                   &
                  Sfx(i,j,k)*d12dx*(-F3*fh(i+1,j,k)-F10*fh(i,j,k)+F18*fh(i-1,j,k) &
                                    -F6*fh(i-2,j,k)+    fh(i-3,j,k))
     elseif(i-2 >= imin)then
!
!              f(i-2) - 8 f(i-1) + 8 f(i+1) - f(i+2)
!  fx(i) = ---------------------------------------------
!                             12 dx
     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                           &
                  Sfx(i,j,k)*d12dx*(fh(i-2,j,k)-EIT*fh(i-1,j,k)+EIT*fh(i+1,j,k)-fh(i+2,j,k))

     elseif(i-1 >= imin)then
!         v
! D f = ------[   3f    + 10f  - 18f    + 6f     - f     ]
!  i     12dx       i+v      i      i-v     i-2v    i-3v
     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                   &
                  Sfx(i,j,k)*d12dx*(-F3*fh(i-1,j,k)-F10*fh(i,j,k)+F18*fh(i+1,j,k) &
                                    -F6*fh(i+2,j,k)+    fh(i+3,j,k))
! set imax and imin 0
     endif
   endif

! y direction   
    if(Sfy(i,j,k) > ZEO)then
      if(j+3 <= jmax)then
!         v
! D f = ------[ - 3f    - 10f  + 18f    - 6f     + f     ]
!  i     12dx       i-v      i      i+v     i+2v    i+3v
     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                   &
                  Sfy(i,j,k)*d12dy*(-F3*fh(i,j-1,k)-F10*fh(i,j,k)+F18*fh(i,j+1,k) &
                                    -F6*fh(i,j+2,k)+    fh(i,j+3,k))
     elseif(j+2 <= jmax)then
!
!              f(i-2) - 8 f(i-1) + 8 f(i+1) - f(i+2)
!  fx(i) = ---------------------------------------------
!                             12 dx
     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                           &
                  Sfy(i,j,k)*d12dy*(fh(i,j-2,k)-EIT*fh(i,j-1,k)+EIT*fh(i,j+1,k)-fh(i,j+2,k))

     elseif(j+1 <= jmax)then
!         v
! D f = ------[   3f    + 10f  - 18f    + 6f     - f     ]
!  i     12dx       i+v      i      i-v     i-2v    i-3v
     f_rhs(i,j,k)=f_rhs(i,j,k)-                                                   &
                  Sfy(i,j,k)*d12dy*(-F3*fh(i,j+1,k)-F10*fh(i,j,k)+F18*fh(i,j-1,k) &
                                    -F6*fh(i,j-2,k)+    fh(i,j-3,k))
! set imax and imin 0
     endif
   elseif(Sfy(i,j,k) < ZEO)then
      if(j-3 >= jmin)then
!         v
! D f = ------[ - 3f    - 10f  + 18f    - 6f     + f     ]
!  i     12dx       i-v      i      i+v     i+2v    i+3v
     f_rhs(i,j,k)=f_rhs(i,j,k)-                                                   &
                  Sfy(i,j,k)*d12dy*(-F3*fh(i,j+1,k)-F10*fh(i,j,k)+F18*fh(i,j-1,k) &
                                    -F6*fh(i,j-2,k)+    fh(i,j-3,k))
     elseif(j-2 >= jmin)then
!
!              f(i-2) - 8 f(i-1) + 8 f(i+1) - f(i+2)
!  fx(i) = ---------------------------------------------
!                             12 dx
     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                           &
                  Sfy(i,j,k)*d12dy*(fh(i,j-2,k)-EIT*fh(i,j-1,k)+EIT*fh(i,j+1,k)-fh(i,j+2,k))

     elseif(j-1 >= jmin)then
!         v
! D f = ------[   3f    + 10f  - 18f    + 6f     - f     ]
!  i     12dx       i+v      i      i-v     i-2v    i-3v
     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                   &
                  Sfy(i,j,k)*d12dy*(-F3*fh(i,j-1,k)-F10*fh(i,j,k)+F18*fh(i,j+1,k) &
                                    -F6*fh(i,j+2,k)+    fh(i,j+3,k))
! set jmax and jmin 0
     endif
   endif

! z direction   
    if(Sfz(i,j,k) > ZEO)then
      if(k+3 <= kmax)then
!         v
! D f = ------[ - 3f    - 10f  + 18f    - 6f     + f     ]
!  i     12dx       i-v      i      i+v     i+2v    i+3v
     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                   &
                  Sfz(i,j,k)*d12dz*(-F3*fh(i,j,k-1)-F10*fh(i,j,k)+F18*fh(i,j,k+1) &
                                    -F6*fh(i,j,k+2)+    fh(i,j,k+3))
     elseif(k+2 <= kmax)then
!
!              f(i-2) - 8 f(i-1) + 8 f(i+1) - f(i+2)
!  fx(i) = ---------------------------------------------
!                             12 dx
     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                           &
                  Sfz(i,j,k)*d12dz*(fh(i,j,k-2)-EIT*fh(i,j,k-1)+EIT*fh(i,j,k+1)-fh(i,j,k+2))

     elseif(k+1 <= kmax)then
!         v
! D f = ------[   3f    + 10f  - 18f    + 6f     - f     ]
!  i     12dx       i+v      i      i-v     i-2v    i-3v
     f_rhs(i,j,k)=f_rhs(i,j,k)-                                                   &
                  Sfz(i,j,k)*d12dz*(-F3*fh(i,j,k+1)-F10*fh(i,j,k)+F18*fh(i,j,k-1) &
                                    -F6*fh(i,j,k-2)+    fh(i,j,k-3))
! set imax and imin 0
     endif
   elseif(Sfz(i,j,k) < ZEO)then
      if(k-3 >= kmin)then
!         v
! D f = ------[ - 3f    - 10f  + 18f    - 6f     + f     ]
!  i     12dx       i-v      i      i+v     i+2v    i+3v
     f_rhs(i,j,k)=f_rhs(i,j,k)-                                                   &
                  Sfz(i,j,k)*d12dz*(-F3*fh(i,j,k+1)-F10*fh(i,j,k)+F18*fh(i,j,k-1) &
                                    -F6*fh(i,j,k-2)+    fh(i,j,k-3))
     elseif(k-2 >= kmin)then
!
!              f(i-2) - 8 f(i-1) + 8 f(i+1) - f(i+2)
!  fx(i) = ---------------------------------------------
!                             12 dx
     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                           &
                  Sfz(i,j,k)*d12dz*(fh(i,j,k-2)-EIT*fh(i,j,k-1)+EIT*fh(i,j,k+1)-fh(i,j,k+2))

     elseif(k-1 >= kmin)then
!         v
! D f = ------[   3f    + 10f  - 18f    + 6f     - f     ]
!  i     12dx       i+v      i      i-v     i-2v    i-3v
     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                   &
                  Sfz(i,j,k)*d12dz*(-F3*fh(i,j,k-1)-F10*fh(i,j,k)+F18*fh(i,j,k+1) &
                                    -F6*fh(i,j,k+2)+    fh(i,j,k+3))
! set kmax and kmin 0
     endif
   endif
#endif
  enddo
  enddo
  enddo

  return

  end subroutine lopsided

#elif (ghost_width == 4)
! sixth order code
! Compute advection terms in right hand sides of field equations
!         v
! D f = ------[ 2f     - 24f    - 35f  + 80f    - 30f     + 8f     - f    ]
!  i     60dx     i-2v      i-v      i      i+v      i+2v     i+3v    i+4v
!
! where
!
!        i
!      |B |
! v = -----
!        i
!       B
!
!-----------------------------------------------------------------------------
subroutine lopsided(ex,X,Y,Z,f,f_rhs,Sfx,Sfy,Sfz,Symmetry,SoA)
  implicit none

!~~~~~~> Input parameters:

  integer, intent(in)  :: ex(1:3),Symmetry
  real*8,  intent(in)  :: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8,dimension(ex(1),ex(2),ex(3)),intent(in)   :: f,Sfx,Sfy,Sfz

  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout):: f_rhs
  real*8,dimension(3),intent(in) ::SoA

!~~~~~~> local variables:

  real*8,dimension(-3:ex(1),-3:ex(2),-3:ex(3))   :: fh
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8 :: dX,dY,dZ
  real*8 :: d60dx,d60dy,d60dz,d12dx,d12dy,d12dz,d2dx,d2dy,d2dz
  real*8,  parameter :: ZEO=0.d0,ONE=1.d0, F60=6.d1
  real*8,  parameter :: TWO=2.d0,F24=2.4d1,F35=3.5d1,F80=8.d1,F30=3.d1,EIT=8.d0
  real*8,  parameter ::  F9=9.d0,F45=4.5d1,F12=1.2d1
  real*8,  parameter ::  F10=1.d1,F77=7.7d1,F150=1.5d2,F100=1.d2,F50=5.d1,F15=1.5d1
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  d60dx = ONE/F60/dX
  d60dy = ONE/F60/dY
  d60dz = ONE/F60/dZ

  d12dx = ONE/F12/dX
  d12dy = ONE/F12/dY
  d12dz = ONE/F12/dZ

  d2dx = ONE/TWO/dX
  d2dy = ONE/TWO/dY
  d2dz = ONE/TWO/dZ

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1
  if(Symmetry > NO_SYMM .and. dabs(Z(1)) < dZ) kmin = -3
  if(Symmetry > EQ_SYMM .and. dabs(X(1)) < dX) imin = -3
  if(Symmetry > EQ_SYMM .and. dabs(Y(1)) < dY) jmin = -3

  call symmetry_bd(4,ex,f,fh,SoA)

! upper bound set ex-1 only for efficiency, 
! the loop body will set ex 0 also
  do k=1,ex(3)-1
  do j=1,ex(2)-1
  do i=1,ex(1)-1
! x direction   
    if(Sfx(i,j,k) >= ZEO .and. i+4 <= imax .and. i-2 >= imin)then
!         v
! D f = ------[ 2f     - 24f    - 35f  + 80f    - 30f     + 8f     - f    ]
!  i     60dx     i-2v      i-v      i      i+v      i+2v     i+3v    i+4v
     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                             &
                  Sfx(i,j,k)*d60dx*(TWO*fh(i-2,j,k)-F24*fh(i-1,j,k)-F35*fh(i,j,k)+F80*fh(i+1,j,k) &
                                   -F30*fh(i+2,j,k)+EIT*fh(i+3,j,k)-    fh(i+4,j,k))
    elseif(Sfx(i,j,k) >= ZEO .and. i+5 <= imax .and. i-1 >= imin)then
!         v
! D f = ------[-10f    - 77f  + 150f    - 100f     + 50f     -15f     + 2f    ]
!  i     60dx      i-v      i       i+v       i+2v      i+3v     i+4v    i+5v
     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                                        &
                  Sfx(i,j,k)*d60dx*(-F10*fh(i-1,j,k)-F77*fh(i  ,j,k)+F150*fh(i+1,j,k)-F100*fh(i+2,j,k) &
                                    +F50*fh(i+3,j,k)-F15*fh(i+4,j,k)+ TWO*fh(i+5,j,k))

    elseif(Sfx(i,j,k) <= ZEO .and. i-4 >= imin .and. i+2 <= imax)then
     f_rhs(i,j,k)=f_rhs(i,j,k)-                                                                   &
                  Sfx(i,j,k)*d60dx*(TWO*fh(i+2,j,k)-F24*fh(i+1,j,k)-F35*fh(i,j,k)+F80*fh(i-1,j,k) &
                                   -F30*fh(i-2,j,k)+EIT*fh(i-3,j,k)-    fh(i-4,j,k))
    elseif(Sfx(i,j,k) <= ZEO .and. i-5 >= imin .and. i+1 <= imax)then
     f_rhs(i,j,k)=f_rhs(i,j,k)-                                                                        &
                  Sfx(i,j,k)*d60dx*(-F10*fh(i+1,j,k)-F77*fh(i  ,j,k)+F150*fh(i-1,j,k)-F100*fh(i-2,j,k) &
                                    +F50*fh(i-3,j,k)-F15*fh(i-4,j,k)+ TWO*fh(i-5,j,k))

     elseif(i+3 <= imax .and. i-3 >= imin)then
!           - f(i-3) + 9 f(i-2) - 45 f(i-1) + 45 f(i+1) - 9 f(i+2) + f(i+3)
!  fx(i) = -----------------------------------------------------------------
!                                        60 dx
     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                                              &
                  Sfx(i,j,k)*d60dx*(-fh(i-3,j,k)+F9*fh(i-2,j,k)-F45*fh(i-1,j,k)+F45*fh(i+1,j,k)-F9*fh(i+2,j,k)+fh(i+3,j,k))

     elseif(i+2 <= imax .and. i-2 >= imin)then
!
!              f(i-2) - 8 f(i-1) + 8 f(i+1) - f(i+2)
!  fx(i) = ---------------------------------------------
!                             12 dx
     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                           &
                  Sfx(i,j,k)*d12dx*(fh(i-2,j,k)-EIT*fh(i-1,j,k)+EIT*fh(i+1,j,k)-fh(i+2,j,k))

     elseif(i+1 <= imax .and. i-1 >= imin)then
!
!              - f(i-1) + f(i+1)
!  fx(i) = --------------------------------
!                     2 dx
     f_rhs(i,j,k)=f_rhs(i,j,k) + Sfx(i,j,k)*d2dx*(-fh(i-1,j,k)+fh(i+1,j,k))

! set imax and imin 0
    endif

! y direction   
     if(Sfy(i,j,k) >= ZEO .and. j+4 <= jmax .and. j-2 >= jmin)then

     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                                   &
                  Sfy(i,j,k)*d60dy*(TWO*fh(i,j-2,k)-F24*fh(i,j-1,k)-F35*fh(i,j,k)+F80*fh(i,j+1,k) &
                                   -F30*fh(i,j+2,k)+EIT*fh(i,j+3,k)-    fh(i,j+4,k))
     elseif(Sfy(i,j,k) >= ZEO .and. j+5 <= jmax .and. j-1 >= jmin)then
     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                                        &
                  Sfy(i,j,k)*d60dy*(-F10*fh(i,j-1,k)-F77*fh(i,j  ,k)+F150*fh(i,j+1,k)-F100*fh(i,j+2,k) &
                                    +F50*fh(i,j+3,k)-F15*fh(i,j+4,k)+ TWO*fh(i,j+5,k))

     elseif(Sfy(i,j,k) <= ZEO .and. j-4 >= jmin .and. j+2 <= jmax)then

     f_rhs(i,j,k)=f_rhs(i,j,k)-                                                                   &
                  Sfy(i,j,k)*d60dy*(TWO*fh(i,j+2,k)-F24*fh(i,j+1,k)-F35*fh(i,j,k)+F80*fh(i,j-1,k) &
                                   -F30*fh(i,j-2,k)+EIT*fh(i,j-3,k)-    fh(i,j-4,k))

     elseif(Sfy(i,j,k) <= ZEO .and. j-5 >= jmin .and. j+1 <= jmax)then

     f_rhs(i,j,k)=f_rhs(i,j,k)-                                                                        &
                  Sfy(i,j,k)*d60dy*(-F10*fh(i,j+1,k)-F77*fh(i,j  ,k)+F150*fh(i,j-1,k)-F100*fh(i,j-2,k) &
                                    +F50*fh(i,j-3,k)-F15*fh(i,j-4,k)+ TWO*fh(i,j-5,k))

     elseif(j+3 <= jmax .and. j-3 >= jmin)then
          
     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                                                         &
                  Sfy(i,j,k)*d60dy*(-fh(i,j-3,k)+F9*fh(i,j-2,k)-F45*fh(i,j-1,k)+F45*fh(i,j+1,k)-F9*fh(i,j+2,k)+fh(i,j+3,k))

     elseif(j+2 <= jmax .and. j-2 >= jmin)then

     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                            &
                  Sfy(i,j,k)*d12dy*(fh(i,j-2,k)-EIT*fh(i,j-1,k)+EIT*fh(i,j+1,k)-fh(i,j+2,k))

     elseif(j+1 <= jmax .and. j-1 >= jmin)then

     f_rhs(i,j,k)=f_rhs(i,j,k) + Sfy(i,j,k)*d2dy*(-fh(i,j-1,k)+fh(i,j+1,k))
! set jmin and jmax 0
     endif
!! z direction   
     if(Sfz(i,j,k) >= ZEO .and. k+4 <= kmax .and. k-2 >= kmin)then

     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                                   &
                  Sfz(i,j,k)*d60dz*(TWO*fh(i,j,k-2)-F24*fh(i,j,k-1)-F35*fh(i,j,k)+F80*fh(i,j,k+1) &
                                   -F30*fh(i,j,k+2)+EIT*fh(i,j,k+3)-    fh(i,j,k+4))
     elseif(Sfz(i,j,k) >= ZEO .and. k+5 <= kmax .and. k-1 >= kmin)then
     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                                        &
                  Sfz(i,j,k)*d60dz*(-F10*fh(i,j,k-1)-F77*fh(i,j,k  )+F150*fh(i,j,k+1)-F100*fh(i,j,k+2) &
                                    +F50*fh(i,j,k+3)-F15*fh(i,j,k+4)+ TWO*fh(i,j,k+5))

     elseif(Sfz(i,j,k) <= ZEO .and. k-4 >= kmin .and. k+2 <= kmax)then

     f_rhs(i,j,k)=f_rhs(i,j,k)-                                                                   &
                  Sfz(i,j,k)*d60dz*(TWO*fh(i,j,k+2)-F24*fh(i,j,k+1)-F35*fh(i,j,k)+F80*fh(i,j,k-1) &
                                   -F30*fh(i,j,k-2)+EIT*fh(i,j,k-3)-    fh(i,j,k-4))

     elseif(Sfz(i,j,k) <= ZEO .and. k-5 >= kmin .and. k+1 <= kmax)then

     f_rhs(i,j,k)=f_rhs(i,j,k)-                                                                        &
                  Sfz(i,j,k)*d60dz*(-F10*fh(i,j,k+1)-F77*fh(i,j,k  )+F150*fh(i,j,k-1)-F100*fh(i,j,k-2) &
                                    +F50*fh(i,j,k-3)-F15*fh(i,j,k-4)+ TWO*fh(i,j,k-5))
     
     elseif(k+3 <= kmax .and. k-3 >= kmin)then

     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                                                         &
                  Sfz(i,j,k)*d60dz*(-fh(i,j,k-3)+F9*fh(i,j,k-2)-F45*fh(i,j,k-1)+F45*fh(i,j,k+1)-F9*fh(i,j,k+2)+fh(i,j,k+3))

     elseif(k+2 <= kmax .and. k-2 >= kmin)then

     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                            &
                  Sfz(i,j,k)*d12dz*(fh(i,j,k-2)-EIT*fh(i,j,k-1)+EIT*fh(i,j,k+1)-fh(i,j,k+2))

     elseif(k+1 <= kmax .and. k-1 >= kmin)then

     f_rhs(i,j,k)=f_rhs(i,j,k)+Sfz(i,j,k)*d2dz*(-fh(i,j,k-1)+fh(i,j,k+1))
! set kmin and kmax 0
     endif

  enddo
  enddo
  enddo

  return

  end subroutine lopsided

#elif (ghost_width == 5)
! eighth order code
!-----------------------------------------------------------------------------
! PRD 77, 024034 (2008)
! Compute advection terms in right hand sides of field equations
!        v [ - 5 f(i-3v) + 60 f(i-2v) - 420 f(i-v) - 378 f(i) + 1050 f(i+v) - 420 f(i+2v) + 140 f(i+3v) - 30 f(i+4v) + 3 f(i+5v)]
! D f = --------------------------------------------------------------------------------------------------------------------------
!  i                                                             840 dx           
!
! where
!
!        i
!      |B |
! v = -----
!        i
!       B
!
!-----------------------------------------------------------------------------
subroutine lopsided(ex,X,Y,Z,f,f_rhs,Sfx,Sfy,Sfz,Symmetry,SoA)
  implicit none

!~~~~~~> Input parameters:

  integer, intent(in)  :: ex(1:3),Symmetry
  real*8,  intent(in)  :: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8,dimension(ex(1),ex(2),ex(3)),intent(in)   :: f,Sfx,Sfy,Sfz

  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout):: f_rhs
  real*8,dimension(3),intent(in) ::SoA

!~~~~~~> local variables:

  real*8,dimension(-4:ex(1),-4:ex(2),-4:ex(3))   :: fh
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8 :: dX,dY,dZ
  real*8 :: d840dx,d840dy,d840dz,d60dx,d60dy,d60dz,d12dx,d12dy,d12dz,d2dx,d2dy,d2dz
  real*8,  parameter :: ZEO=0.d0,ONE=1.d0, F60=6.d1
  real*8,  parameter :: TWO=2.d0,F30=3.d1,EIT=8.d0
  real*8,  parameter ::  F9=9.d0,F45=4.5d1,F12=1.2d1,F140=1.4d2,THR=3.d0
  real*8,  parameter :: F840=8.4d2,F5=5.d0,F420=4.2d2,F378=3.78d2,F1050=1.05d3
  real*8,  parameter :: F32=3.2d1,F168=1.68d2,F672=6.72d2
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

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

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1
  if(Symmetry > NO_SYMM .and. dabs(Z(1)) < dZ) kmin = -4
  if(Symmetry > EQ_SYMM .and. dabs(X(1)) < dX) imin = -4
  if(Symmetry > EQ_SYMM .and. dabs(Y(1)) < dY) jmin = -4

  call symmetry_bd(5,ex,f,fh,SoA)

! upper bound set ex-1 only for efficiency, 
! the loop body will set ex 0 also
  do k=1,ex(3)-1
  do j=1,ex(2)-1
  do i=1,ex(1)-1
! x direction   
    if(Sfx(i,j,k) >= ZEO .and. i+5 <= imax .and. i-3 >= imin)then
!        v [ - 5 f(i-3v) + 60 f(i-2v) - 420 f(i-v) - 378 f(i) + 1050 f(i+v) - 420 f(i+2v) + 140 f(i+3v) - 30 f(i+4v) + 3 f(i+5v)]
! D f = --------------------------------------------------------------------------------------------------------------------------
!  i                                                             840 dx    
     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                                         &
                  Sfx(i,j,k)*d840dx*(-F5*fh(i-3,j,k)+F60 *fh(i-2,j,k)-F420*fh(i-1,j,k)-F378*fh(i  ,j,k) &
                                  +F1050*fh(i+1,j,k)-F420*fh(i+2,j,k)+F140*fh(i+3,j,k)-F30 *fh(i+4,j,k)+THR*fh(i+5,j,k))

    elseif(Sfx(i,j,k) <= ZEO .and. i-5 >= imin .and. i+3 <= imax)then
     f_rhs(i,j,k)=f_rhs(i,j,k)-                                                                          &
                  Sfx(i,j,k)*d840dx*(-F5*fh(i+3,j,k)+F60 *fh(i+2,j,k)-F420*fh(i+1,j,k)-F378*fh(i   ,j,k) &
                                  +F1050*fh(i-1,j,k)-F420*fh(i-2,j,k)+F140*fh(i-3,j,k)- F30*fh(i-4,j,k)+THR*fh(i-5,j,k))

    elseif(i+4 <= imax .and. i-4 >= imin)then
!           3 f(i-4) - 32 f(i-3) + 168 f(i-2) - 672 f(i-1) + 672 f(i+1) - 168 f(i+2) + 32 f(i+3) - 3 f(i+4)
!  fx(i) = -------------------------------------------------------------------------------------------------
!                                                        840 dx
     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                                              &
                  Sfx(i,j,k)*d840dx*( THR*fh(i-4,j,k)-F32 *fh(i-3,j,k)+F168*fh(i-2,j,k)-F672*fh(i-1,j,k)+    &
                                     F672*fh(i+1,j,k)-F168*fh(i+2,j,k)+F32 *fh(i+3,j,k)-THR *fh(i+4,j,k))

     elseif(i+3 <= imax .and. i-3 >= imin)then
!           - f(i-3) + 9 f(i-2) - 45 f(i-1) + 45 f(i+1) - 9 f(i+2) + f(i+3)
!  fx(i) = -----------------------------------------------------------------
!                                        60 dx
     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                                              &
                  Sfx(i,j,k)*d60dx*(-fh(i-3,j,k)+F9*fh(i-2,j,k)-F45*fh(i-1,j,k)+F45*fh(i+1,j,k)-F9*fh(i+2,j,k)+fh(i+3,j,k))

     elseif(i+2 <= imax .and. i-2 >= imin)then
!
!              f(i-2) - 8 f(i-1) + 8 f(i+1) - f(i+2)
!  fx(i) = ---------------------------------------------
!                             12 dx
     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                           &
                  Sfx(i,j,k)*d12dx*(fh(i-2,j,k)-EIT*fh(i-1,j,k)+EIT*fh(i+1,j,k)-fh(i+2,j,k))

     elseif(i+1 <= imax .and. i-1 >= imin)then
!
!              - f(i-1) + f(i+1)
!  fx(i) = --------------------------------
!                     2 dx
     f_rhs(i,j,k)=f_rhs(i,j,k) + Sfx(i,j,k)*d2dx*(-fh(i-1,j,k)+fh(i+1,j,k))

! set imax and imin 0
    endif

! y direction   
    if(Sfy(i,j,k) >= ZEO .and. j+5 <= jmax .and. j-3 >= jmin)then

     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                                         &
                  Sfy(i,j,k)*d840dy*(-F5*fh(i,j-3,k)+F60 *fh(i,j-2,k)-F420*fh(i,j-1,k)-F378*fh(i,j  ,k) &
                                  +F1050*fh(i,j+1,k)-F420*fh(i,j+2,k)+F140*fh(i,j+3,k)-F30 *fh(i,j+4,k)+THR*fh(i,j+5,k))

    elseif(Sfy(i,j,k) <= ZEO .and. j-5 >= jmin .and. j+3 <= jmax)then
     f_rhs(i,j,k)=f_rhs(i,j,k)-                                                                         &
                  Sfy(i,j,k)*d840dy*(-F5*fh(i,j+3,k)+F60 *fh(i,j+2,k)-F420*fh(i,j+1,k)-F378*fh(i,j  ,k) &
                                  +F1050*fh(i,j-1,k)-F420*fh(i,j-2,k)+F140*fh(i,j-3,k)- F30*fh(i,j-4,k)+THR*fh(i,j-5,k))

    elseif(j+4 <= jmax .and. j-4 >= jmin)then

     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                                              &
                  Sfy(i,j,k)*d840dy*( THR*fh(i,j-4,k)-F32 *fh(i,j-3,k)+F168*fh(i,j-2,k)-F672*fh(i,j-1,k)+    &
                                     F672*fh(i,j+1,k)-F168*fh(i,j+2,k)+F32 *fh(i,j+3,k)-THR *fh(i,j+4,k))

     elseif(j+3 <= jmax .and. j-3 >= jmin)then
          
     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                                                         &
                  Sfy(i,j,k)*d60dy*(-fh(i,j-3,k)+F9*fh(i,j-2,k)-F45*fh(i,j-1,k)+F45*fh(i,j+1,k)-F9*fh(i,j+2,k)+fh(i,j+3,k))

     elseif(j+2 <= jmax .and. j-2 >= jmin)then

     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                            &
                  Sfy(i,j,k)*d12dy*(fh(i,j-2,k)-EIT*fh(i,j-1,k)+EIT*fh(i,j+1,k)-fh(i,j+2,k))

     elseif(j+1 <= jmax .and. j-1 >= jmin)then

     f_rhs(i,j,k)=f_rhs(i,j,k) + Sfy(i,j,k)*d2dy*(-fh(i,j-1,k)+fh(i,j+1,k))
! set jmin and jmax 0
     endif
!! z direction   
    if(Sfz(i,j,k) >= ZEO .and. k+5 <= kmax .and. k-3 >= kmin)then

     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                                         &
                  Sfz(i,j,k)*d840dz*(-F5*fh(i,j,k-3)+F60 *fh(i,j,k-2)-F420*fh(i,j,k-1)-F378*fh(i,j,k  ) &
                                  +F1050*fh(i,j,k+1)-F420*fh(i,j,k+2)+F140*fh(i,j,k+3)-F30 *fh(i,j,k+4)+THR*fh(i,j,k+5))

    elseif(Sfz(i,j,k) <= ZEO .and. k-5 >= kmin .and. k+3 <= kmax)then
     f_rhs(i,j,k)=f_rhs(i,j,k)-                                                                         &
                  Sfz(i,j,k)*d840dz*(-F5*fh(i,j,k+3)+F60 *fh(i,j,k+2)-F420*fh(i,j,k+1)-F378*fh(i,j,k  ) &
                                  +F1050*fh(i,j,k-1)-F420*fh(i,j,k-2)+F140*fh(i,j,k-3)- F30*fh(i,j,k-4)+THR*fh(i,j,k-5))

    elseif(k+4 <= kmax .and. k-4 >= kmin)then

     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                                              &
                  Sfz(i,j,k)*d840dz*( THR*fh(i,j,k-4)-F32 *fh(i,j,k-3)+F168*fh(i,j,k-2)-F672*fh(i,j,k-1)+    &
                                     F672*fh(i,j,k+1)-F168*fh(i,j,k+2)+F32 *fh(i,j,k+3)-THR *fh(i,j,k+4))
     
     elseif(k+3 <= kmax .and. k-3 >= kmin)then

     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                                                         &
                  Sfz(i,j,k)*d60dz*(-fh(i,j,k-3)+F9*fh(i,j,k-2)-F45*fh(i,j,k-1)+F45*fh(i,j,k+1)-F9*fh(i,j,k+2)+fh(i,j,k+3))

     elseif(k+2 <= kmax .and. k-2 >= kmin)then

     f_rhs(i,j,k)=f_rhs(i,j,k)+                                                            &
                  Sfz(i,j,k)*d12dz*(fh(i,j,k-2)-EIT*fh(i,j,k-1)+EIT*fh(i,j,k+1)-fh(i,j,k+2))

     elseif(k+1 <= kmax .and. k-1 >= kmin)then

     f_rhs(i,j,k)=f_rhs(i,j,k)+Sfz(i,j,k)*d2dz*(-fh(i,j,k-1)+fh(i,j,k+1))
! set kmin and kmax 0
     endif

  enddo
  enddo
  enddo

  return

  end subroutine lopsided

#endif  
