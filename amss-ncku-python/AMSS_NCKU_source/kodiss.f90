

#include "macrodef.fh"

! we need only distinguish different finite difference order
! Vertex or Cell is distinguished in routine symmetry_bd which locates in
! file "fmisc.f90"

#if (ghost_width == 2)
! second order code

!------------------------------------------------------------------------------------------------------------------------------
!usual type Kreiss-Oliger type numerical dissipation
!We support cell center only
!  (D_+D_-)^2 =
!   f(i-2) - 4 f(i-1) + 6 f(i) - 4 f(i+1) + f(i+2)
! ------------------------------------------------------
!                       dx^4
!------------------------------------------------------------------------------------------------------------------------------
! do not add dissipation near boundary
subroutine kodis(ex,X,Y,Z,f,f_rhs,SoA,Symmetry,eps)

implicit none
! argument variables
integer,intent(in) :: Symmetry
integer,dimension(3),intent(in)::ex
real*8, dimension(1:3), intent(in) :: SoA
double precision,intent(in),dimension(ex(1))::X
double precision,intent(in),dimension(ex(2))::Y
double precision,intent(in),dimension(ex(3))::Z
double precision,intent(in),dimension(ex(1),ex(2),ex(3))::f
double precision,intent(inout),dimension(ex(1),ex(2),ex(3))::f_rhs
real*8,intent(in) :: eps

!~~~~~~ other variables

  real*8 :: dX,dY,dZ
  real*8,dimension(-1:ex(1),-1:ex(2),-1:ex(3))   :: fh
  integer :: imin,jmin,kmin,imax,jmax,kmax
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8,parameter   :: cof = 1.6d1 ! 2^4
  real*8,  parameter :: F4=4.d0,F6=6.d0
  integer::i,j,k

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

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

!   f(i-2) - 4 f(i-1) + 6 f(i) - 4 f(i+1) + f(i+2)
! ------------------------------------------------------
!                       dx^4

!  note the sign (-1)^r-1, now r=2
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)

  if(i-2 >= imin .and. i+2 <= imax .and. &
     j-2 >= jmin .and. j+2 <= jmax .and. &
     k-2 >= kmin .and. k+2 <= kmax) then
! x direction
   f_rhs(i,j,k)       = f_rhs(i,j,k) - eps/dX/cof * (     &
                                (fh(i-2,j,k)+fh(i+2,j,k)) &
                         - F4 * (fh(i-1,j,k)+fh(i+1,j,k)) &
                         + F6 *  fh(i,j,k) )
! y direction

   f_rhs(i,j,k)       = f_rhs(i,j,k) - eps/dY/cof * (     &
                                (fh(i,j-2,k)+fh(i,j+2,k)) &
                         - F4 * (fh(i,j-1,k)+fh(i,j+1,k)) &
                         + F6 *  fh(i,j,k) )
! z direction

   f_rhs(i,j,k)       = f_rhs(i,j,k) - eps/dZ/cof * (     &
                                (fh(i,j,k-2)+fh(i,j,k+2)) &
                         - F4 * (fh(i,j,k-1)+fh(i,j,k+1)) &
                         + F6 *  fh(i,j,k) )

  endif

  enddo
  enddo
  enddo

  return

end subroutine kodis

#elif (ghost_width == 3)
! fourth order code

!---------------------------------------------------------------------------------------------
!usual type Kreiss-Oliger type numerical dissipation
!We support cell center only
! Note the notation D_+ and D_- [P240 of B. Gustafsson, H.-O. Kreiss, and J. Oliger, Time
! Dependent Problems and Difference Methods (Wiley, New York, 1995).]
! D_+ = (f(i+1) - f(i))/h
! D_- = (f(i) - f(i-1))/h
! then we have D_+D_- = D_-D_+
!              D_+^3D_-^3 = (D_+D_-)^3 =
!    f(i-3) - 6 f(i-2) + 15 f(i-1) - 20 f(i) + 15 f(i+1) - 6 f(i+2) + f(i+3)
! -----------------------------------------------------------------------------
!                                    dx^6
! this is for 4th order accurate finite difference scheme
!---------------------------------------------------------------------------------------------
subroutine kodis(ex,X,Y,Z,f,f_rhs,SoA,Symmetry,eps)

implicit none
! argument variables
integer,intent(in) :: Symmetry
integer,dimension(3),intent(in)::ex
real*8, dimension(1:3), intent(in) :: SoA
double precision,intent(in),dimension(ex(1))::X
double precision,intent(in),dimension(ex(2))::Y
double precision,intent(in),dimension(ex(3))::Z
double precision,intent(in),dimension(ex(1),ex(2),ex(3))::f
double precision,intent(inout),dimension(ex(1),ex(2),ex(3))::f_rhs
real*8,intent(in) :: eps
! local variables
real*8,dimension(-2:ex(1),-2:ex(2),-2:ex(3))   :: fh
integer :: imin,jmin,kmin,imax,jmax,kmax
integer :: i,j,k
real*8  :: dX,dY,dZ
real*8, parameter :: ONE=1.d0,SIX=6.d0,FIT=1.5d1,TWT=2.d1
real*8,parameter::cof=6.4d1   ! 2^6
integer, parameter :: NO_SYMM=0, OCTANT=2

!rhs_i = rhs_i + eps/dx/cof*(f_i-3 - 6*f_i-2 + 15*f_i-1 - 20*f_i + 15*f_i+1 - 6*f_i+2 + f_i+3)

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)
  
  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry > NO_SYMM .and. dabs(Z(1)) < dZ) kmin = -2
  if(Symmetry == OCTANT .and. dabs(X(1)) < dX) imin = -2
  if(Symmetry == OCTANT .and. dabs(Y(1)) < dY) jmin = -2

  call symmetry_bd(3,ex,f,fh,SoA)

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)

  if(i-3 >= imin .and. i+3 <= imax .and. &
     j-3 >= jmin .and. j+3 <= jmax .and. &
     k-3 >= kmin .and. k+3 <= kmax) then
#if 0     
! x direction
   f_rhs(i,j,k)       = f_rhs(i,j,k) + eps/dX/cof * (     &
                              (fh(i-3,j,k)+fh(i+3,j,k)) - &
                          SIX*(fh(i-2,j,k)+fh(i+2,j,k)) + &
                          FIT*(fh(i-1,j,k)+fh(i+1,j,k)) - &
                          TWT* fh(i,j,k)            )
! y direction

   f_rhs(i,j,k)       = f_rhs(i,j,k) + eps/dY/cof * (     &
                              (fh(i,j-3,k)+fh(i,j+3,k)) - &
                          SIX*(fh(i,j-2,k)+fh(i,j+2,k)) + &
                          FIT*(fh(i,j-1,k)+fh(i,j+1,k)) - &
                          TWT* fh(i,j,k)            )
! z direction

   f_rhs(i,j,k)       = f_rhs(i,j,k) + eps/dZ/cof * (     &
                              (fh(i,j,k-3)+fh(i,j,k+3)) - &
                          SIX*(fh(i,j,k-2)+fh(i,j,k+2)) + &
                          FIT*(fh(i,j,k-1)+fh(i,j,k+1)) - &
                          TWT* fh(i,j,k)            )
#else
! calculation order if important ?
   f_rhs(i,j,k)       = f_rhs(i,j,k) + eps/cof *( (     &
                              (fh(i-3,j,k)+fh(i+3,j,k)) - &
                          SIX*(fh(i-2,j,k)+fh(i+2,j,k)) + &
                          FIT*(fh(i-1,j,k)+fh(i+1,j,k)) - &
                          TWT* fh(i,j,k)            )/dX + &
                                                  (     &
                              (fh(i,j-3,k)+fh(i,j+3,k)) - &
                          SIX*(fh(i,j-2,k)+fh(i,j+2,k)) + &
                          FIT*(fh(i,j-1,k)+fh(i,j+1,k)) - &
                          TWT* fh(i,j,k)            )/dY + &
                                                  (     &
                              (fh(i,j,k-3)+fh(i,j,k+3)) - &
                          SIX*(fh(i,j,k-2)+fh(i,j,k+2)) + &
                          FIT*(fh(i,j,k-1)+fh(i,j,k+1)) - &
                          TWT* fh(i,j,k)            )/dZ )
#endif
  endif

  enddo
  enddo
  enddo

  return

  end subroutine kodis

#elif (ghost_width == 4)
! sixth order code
!------------------------------------------------------------------------------------------------------------------------------
!usual type Kreiss-Oliger type numerical dissipation
!We support cell center only
!  (D_+D_-)^4 =
!   f(i-4) - 8 f(i-3) + 28 f(i-2) - 56 f(i-1) + 70 f(i) - 56 f(i+1) + 28 f(i+2) - 8 f(i+3) + f(i+4)
! ----------------------------------------------------------------------------------------------------------
!                                              dx^8
!------------------------------------------------------------------------------------------------------------------------------
! do not add dissipation near boundary
subroutine kodis(ex,X,Y,Z,f,f_rhs,SoA,Symmetry,eps)

implicit none
! argument variables
integer,intent(in) :: Symmetry
integer,dimension(3),intent(in)::ex
real*8, dimension(1:3), intent(in) :: SoA
double precision,intent(in),dimension(ex(1))::X
double precision,intent(in),dimension(ex(2))::Y
double precision,intent(in),dimension(ex(3))::Z
double precision,intent(in),dimension(ex(1),ex(2),ex(3))::f
double precision,intent(inout),dimension(ex(1),ex(2),ex(3))::f_rhs
real*8,intent(in) :: eps

!~~~~~~ other variables

  real*8 :: dX,dY,dZ
  real*8,dimension(-3:ex(1),-3:ex(2),-3:ex(3))   :: fh
  integer :: imin,jmin,kmin,imax,jmax,kmax
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8,parameter   :: cof = 2.56d2 ! 2^8
  real*8,  parameter :: F8=8.d0,F28=2.8d1,F56=5.6d1,F70=7.d1
  integer::i,j,k

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)
  
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

!   f(i-4) - 8 f(i-3) + 28 f(i-2) - 56 f(i-1) + 70 f(i) - 56 f(i+1) + 28 f(i+2) - 8 f(i+3) + f(i+4)
! ----------------------------------------------------------------------------------------------------------
!                                              dx^8

!  note the sign (-1)^r-1, now r=4
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)

  if(i>imin+3 .and. i < imax-3 .and. &
     j>jmin+3 .and. j < jmax-3 .and. &
     k>kmin+3 .and. k < kmax-3) then
! x direction
   f_rhs(i,j,k)       = f_rhs(i,j,k) - eps/dX/cof * (     &
                                (fh(i-4,j,k)+fh(i+4,j,k)) &
                         - F8 * (fh(i-3,j,k)+fh(i+3,j,k)) &
                         +F28 * (fh(i-2,j,k)+fh(i+2,j,k)) &
                         -F56 * (fh(i-1,j,k)+fh(i+1,j,k)) &
                         +F70 *  fh(i,j,k) )
! y direction

   f_rhs(i,j,k)       = f_rhs(i,j,k) - eps/dY/cof * (     &
                                (fh(i,j-4,k)+fh(i,j+4,k)) &
                         - F8 * (fh(i,j-3,k)+fh(i,j+3,k)) &
                         +F28 * (fh(i,j-2,k)+fh(i,j+2,k)) &
                         -F56 * (fh(i,j-1,k)+fh(i,j+1,k)) &
                         +F70 *  fh(i,j,k) )
! z direction

   f_rhs(i,j,k)       = f_rhs(i,j,k) - eps/dZ/cof * (     &
                                (fh(i,j,k-4)+fh(i,j,k+4)) &
                         - F8 * (fh(i,j,k-3)+fh(i,j,k+3)) &
                         +F28 * (fh(i,j,k-2)+fh(i,j,k+2)) &
                         -F56 * (fh(i,j,k-1)+fh(i,j,k+1)) &
                         +F70 *  fh(i,j,k) )

  endif

  enddo
  enddo
  enddo

  return

end subroutine kodis

#elif (ghost_width == 5)
! eighth order code
!------------------------------------------------------------------------------------------------------------------------------
!usual type Kreiss-Oliger type numerical dissipation
!We support cell center only
! Note the notation D_+ and D_- [P240 of B. Gustafsson, H.-O. Kreiss, and J. Oliger, Time
! Dependent Problems and Difference Methods (Wiley, New York, 1995).]
! D_+ = (f(i+1) - f(i))/h
! D_- = (f(i) - f(i-1))/h
! then we have D_+D_- = D_-D_+ = (f(i+1) - 2f(i) + f(i-1))/h^2
! for nth order accurate finite difference code, we need r =n/2+1
!              D_+^rD_-^r = (D_+D_-)^r 
! following the tradiation of PRD 77, 024027 (BB's calibration paper, Eq.(64),
!  correct some typo according to above book) :
! + eps*(-1)^(r-1)*h^(2r-1)/2^(2r)*(D_+D_-)^r 
!
!
! this is for 8th order accurate finite difference scheme
!  (D_+D_-)^5 =
!  f(i-5) - 10 f(i-4) + 45 f(i-3) - 120 f(i-2) + 210 f(i-1) - 252 f(i) + 210 f(i+1) - 120 f(i+2) + 45 f(i+3) - 10 f(i+4) + f(i+5)
! -------------------------------------------------------------------------------------------------------------------------------
!                                                              dx^10
!---------------------------------------------------------------------------------------------------------------------------------
! do not add dissipation near boundary
subroutine kodis(ex,X,Y,Z,f,f_rhs,SoA,Symmetry,eps)

implicit none
! argument variables
integer,intent(in) :: Symmetry
integer,dimension(3),intent(in)::ex
real*8, dimension(1:3), intent(in) :: SoA
double precision,intent(in),dimension(ex(1))::X
double precision,intent(in),dimension(ex(2))::Y
double precision,intent(in),dimension(ex(3))::Z
double precision,intent(in),dimension(ex(1),ex(2),ex(3))::f
double precision,intent(inout),dimension(ex(1),ex(2),ex(3))::f_rhs
real*8,intent(in) :: eps

!~~~~~~ other variables

  real*8 :: dX,dY,dZ
  real*8,dimension(-4:ex(1),-4:ex(2),-4:ex(3))   :: fh
  integer :: imin,jmin,kmin,imax,jmax,kmax
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8,parameter   :: cof = 1.024d3 ! 2^2r = 2^10
  real*8,  parameter :: F10=1.d1,F45=4.5d1,F120=1.2d2,F210=2.1d2,F252=2.52d2
  integer::i,j,k

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)
  
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

!  f(i-5) - 10 f(i-4) + 45 f(i-3) - 120 f(i-2) + 210 f(i-1) - 252 f(i) + 210 f(i+1) - 120 f(i+2) + 45 f(i+3) - 10 f(i+4) + f(i+5)
! -------------------------------------------------------------------------------------------------------------------------------
!                                                              dx^10

!  note the sign (-1)^r-1, now r=5
  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)

  if(i>imin+4 .and. i < imax-4 .and. &
     j>jmin+4 .and. j < jmax-4 .and. &
     k>kmin+4 .and. k < kmax-4) then
! x direction
   f_rhs(i,j,k)       = f_rhs(i,j,k) + eps/dX/cof * (      &
                                 (fh(i-5,j,k)+fh(i+5,j,k)) &
                         - F10 * (fh(i-4,j,k)+fh(i+4,j,k)) &
                         + F45 * (fh(i-3,j,k)+fh(i+3,j,k)) &
                         - F120* (fh(i-2,j,k)+fh(i+2,j,k)) &
                         + F210* (fh(i-1,j,k)+fh(i+1,j,k)) &
                         - F252 * fh(i,j,k) )
! y direction

   f_rhs(i,j,k)       = f_rhs(i,j,k) + eps/dY/cof * (      &
                                 (fh(i,j-5,k)+fh(i,j+5,k)) &
                         - F10 * (fh(i,j-4,k)+fh(i,j+4,k)) &
                         + F45 * (fh(i,j-3,k)+fh(i,j+3,k)) &
                         - F120* (fh(i,j-2,k)+fh(i,j+2,k)) &
                         + F210* (fh(i,j-1,k)+fh(i,j+1,k)) &
                         - F252 * fh(i,j,k) )
! z direction

   f_rhs(i,j,k)       = f_rhs(i,j,k) + eps/dZ/cof * (      &
                                 (fh(i,j,k-5)+fh(i,j,k+5)) &
                         - F10 * (fh(i,j,k-4)+fh(i,j,k+4)) &
                         + F45 * (fh(i,j,k-3)+fh(i,j,k+3)) &
                         - F120* (fh(i,j,k-2)+fh(i,j,k+2)) &
                         + F210* (fh(i,j,k-1)+fh(i,j,k+1)) &
                         - F252 * fh(i,j,k) )

  endif

  enddo
  enddo
  enddo

  return

end subroutine kodis

#endif  
