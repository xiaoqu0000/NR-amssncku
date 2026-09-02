! Interpolation functions extracted from fmisc.f90
! Compiled without -fp-model fast=2 (small loops, fast=2 has overhead)

#include "macrodef.fh"

#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
!-----------------------------------------------------------------------------------------------------------------
! three dimensional interpolation for vertex center grid structure  
  subroutine global_interp(ex,X,Y,Z,f,f_int,x1,y1,z1,ORDN,SoA,symmetry)
  implicit none

!~~~~~~> Input parameters:

  integer,                             intent(in) :: ex(1:3), symmetry,ORDN
  real*8,intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: f
  real*8,                              intent(out):: f_int
  real*8,                              intent(in) :: x1,y1,z1
  real*8, dimension(3),                intent(in) :: SoA

!~~~~~~> Other parameters:

  integer :: j,m,imin,jmin,kmin
  integer,dimension(3) :: cxB,cxT,cxI,cmin,cmax
  real*8,dimension(3) :: cx
  real*8, dimension(1:ORDN) :: x1a
  real*8, dimension(1:ORDN,1:ORDN,1:ORDN) :: ya
  integer, parameter :: NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2
  real*8 :: dX,dY,dZ,ddy
  real*8, parameter :: ONE=1.d0
  logical::decide3d

  imin = lbound(f,1)
  jmin = lbound(f,2)
  kmin = lbound(f,3)

  dX = X(imin+1)-X(imin)
  dY = Y(jmin+1)-Y(jmin)
  dZ = Z(kmin+1)-Z(kmin)

  forall( j = 1:ordn ) x1a(j) = ( j - 1 )* ONE

  cxI(1) = idint((x1-X(1))/dX+0.4)+1
  cxI(2) = idint((y1-Y(1))/dY+0.4)+1
  cxI(3) = idint((z1-Z(1))/dZ+0.4)+1

  cxB = cxI - ORDN/2+1
  cxT = cxB + ORDN - 1
       
  cmin = 1
  cmax = ex
  if(Symmetry == OCTANT  .and.dabs(X(1))<dX) cmin(1) = -ORDN/2+2
  if(Symmetry == OCTANT  .and.dabs(Y(1))<dY) cmin(2) = -ORDN/2+2
  if(Symmetry /= NO_SYMM .and.dabs(Z(1))<dZ) cmin(3) = -ORDN/2+2
  do m =1,3
   if(cxB(m) < cmin(m))then
      cxB(m) = cmin(m)
      cxT(m) = cxB(m) + ORDN - 1
   endif
   if(cxT(m) > cmax(m))then
      cxT(m) = cmax(m)
      cxB(m) = cxT(m) + 1 - ORDN
   endif
 enddo
 if(cxB(1)>0)then
  cx(1) = (x1 - X(cxB(1)))/dX
 else
  cx(1) = (x1 + X(2-cxB(1)))/dX
 endif
 if(cxB(2)>0)then
  cx(2) = (y1 - Y(cxB(2)))/dY
 else
  cx(2) = (y1 + Y(2-cxB(2)))/dY
 endif
 if(cxB(3)>0)then
  cx(3) = (z1 - Z(cxB(3)))/dZ
 else
  cx(3) = (z1 + Z(2-cxB(3)))/dZ
 endif

  if(decide3d(ex,f,f,cxB,cxT,SoA,ya,ORDN,Symmetry))then
     write(*,*)"global_interp position: ",x1,y1,z1
     write(*,*)"data range: ",X(1),X(ex(1)),Y(1),Y(ex(2)),Z(1),Z(ex(3))
     stop
  endif
  call polin3(x1a,x1a,x1a,ya,cx(1),cx(2),cx(3),f_int,ddy,ORDN)

  return

  end subroutine global_interp
!----------------------------------------------------------------
! decide which 3d data to be used does not surport PI-Symmetry yet 
!----------------------------------------------------------------
  function decide3d(ex,f,fpi,cxB,cxT,SoA,ya,ORDN,Symmetry)  result(gont)
  implicit none

  integer,                                 intent(in) :: ORDN,Symmetry
  integer,dimension(1:3)                 , intent(in) :: ex,cxB,cxT
  real*8, dimension(1:3)                 , intent(in) :: SoA
  real*8, dimension(ex(1),ex(2),ex(3))   , intent(in) :: f,fpi
  real*8, dimension(cxB(1):cxT(1),cxB(2):cxT(2),cxB(3):cxT(3)), intent(out):: ya
  logical::gont

  integer,dimension(1:3) :: fmin1,fmin2,fmax1,fmax2
  integer::i,j,k,m

  gont=.false.
  do m=1,3
! check cxB and cxT are NaN or not  
    if(.not.(iabs(cxB(m)).ge.0)) gont=.true.
    if(.not.(iabs(cxT(m)).ge.0)) gont=.true.
    fmin1(m) = max(1,cxB(m))
    fmax1(m) = cxT(m)
    fmin2(m) = cxB(m)
    fmax2(m) = min(0,cxT(m))
    if((fmin1(m).le.fmax1(m)).and.(  fmin1(m)<1.or.  fmax1(m)>ex(m)))gont=.true.
    if((fmin2(m).le.fmax2(m)).and.(2-fmax2(m)<1.or.2-fmin2(m)>ex(m)))gont=.true.
  enddo
!sanity check
  if(gont)then
          write(*,*)"error in decide3d"
          write(*,*)((fmin1.le.fmax1).and.(  fmin1<1.or.  fmax1>ex))
          write(*,*)((fmin2.le.fmax2).and.(2-fmax2<1.or.2-fmin2>ex))
          write(*,*)"cxB, cxT and data shape:"
          write(*,*)cxB,cxT,ex
          write(*,*)"resulted fmin1, fmax1 and fmin2, fmax2:"
          write(*,*)fmin1,fmax1,fmin2,fmax2
  else

! Optimized for Intel ifx SIMD vectorization (AVX-512)
! Use explicit loops with OMP SIMD to enable ZMM register usage
! Main region: direct copy (continuous memory access)
  do k=fmin1(3),fmax1(3)
     do j=fmin1(2),fmax1(2)
!DIR$ VECTOR ALWAYS

        do i=fmin1(1),fmax1(1)
           ya(i,j,k) = f(i,j,k)
        enddo
     enddo
  enddo

! X-direction symmetry (vertex uses 2-i)
  do k=fmin1(3),fmax1(3)
     do j=fmin1(2),fmax1(2)
!DIR$ VECTOR ALWAYS

        do i=fmin2(1),fmax2(1)
           ya(i,j,k) = f(2-i,j,k)*SoA(1)
        enddo
     enddo
  enddo

! Y-direction symmetry
  do k=fmin1(3),fmax1(3)
     do j=fmin2(2),fmax2(2)
!DIR$ VECTOR ALWAYS

        do i=fmin1(1),fmax1(1)
           ya(i,j,k) = f(i,2-j,k)*SoA(2)
        enddo
     enddo
  enddo

! XY-direction symmetry
  do k=fmin1(3),fmax1(3)
     do j=fmin2(2),fmax2(2)
!DIR$ VECTOR ALWAYS

        do i=fmin2(1),fmax2(1)
           ya(i,j,k) = f(2-i,2-j,k)*SoA(1)*SoA(2)
        enddo
     enddo
  enddo

! Z-direction symmetry
  do k=fmin2(3),fmax2(3)
     do j=fmin1(2),fmax1(2)
!DIR$ VECTOR ALWAYS

        do i=fmin1(1),fmax1(1)
           ya(i,j,k) = f(i,j,2-k)*SoA(3)
        enddo
     enddo
  enddo

! XZ-direction symmetry
  do k=fmin2(3),fmax2(3)
     do j=fmin1(2),fmax1(2)
!DIR$ VECTOR ALWAYS

        do i=fmin2(1),fmax2(1)
           ya(i,j,k) = f(2-i,j,2-k)*SoA(1)*SoA(3)
        enddo
     enddo
  enddo

! YZ-direction symmetry
  do k=fmin2(3),fmax2(3)
     do j=fmin2(2),fmax2(2)
!DIR$ VECTOR ALWAYS

        do i=fmin1(1),fmax1(1)
           ya(i,j,k) = f(i,2-j,2-k)*SoA(2)*SoA(3)
        enddo
     enddo
  enddo

! XYZ-direction symmetry
  do k=fmin2(3),fmax2(3)
     do j=fmin2(2),fmax2(2)
!DIR$ VECTOR ALWAYS

        do i=fmin2(1),fmax2(1)
           ya(i,j,k) = f(2-i,2-j,2-k)*SoA(1)*SoA(2)*SoA(3)
        enddo
     enddo
  enddo

  endif

  end function decide3d

#else
#ifdef Cell
!--------------------------------------------------------------------------
! three dimensional interpolation for cell center grid structure  
  subroutine global_interp(ex,X,Y,Z,f,f_int,x1,y1,z1,ORDN,SoA,symmetry)
  implicit none

!~~~~~~> Input parameters:

  integer,                             intent(in) :: ex(1:3), symmetry,ORDN
  real*8,intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: f
  real*8,                              intent(out):: f_int
  real*8,                              intent(in) :: x1,y1,z1
  real*8, dimension(3),                intent(in) :: SoA

!~~~~~~> Other parameters:

  integer :: j,m,imin,jmin,kmin
  integer,dimension(3) :: cxB,cxT,cxI,cmin,cmax
  real*8,dimension(3) :: cx
  real*8, dimension(1:ORDN) :: x1a
  real*8, dimension(1:ORDN,1:ORDN,1:ORDN) :: ya
  integer, parameter :: NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2
  real*8 :: dX,dY,dZ,ddy
  real*8, parameter :: ONE=1.d0
  logical::decide3d

  imin = lbound(f,1)
  jmin = lbound(f,2)
  kmin = lbound(f,3)

  dX = X(imin+1)-X(imin)
  dY = Y(jmin+1)-Y(jmin)
  dZ = Z(kmin+1)-Z(kmin)

  forall( j = 1:ordn ) x1a(j) = ( j - 1 )* ONE

  cxI(1) = idint((x1-X(1))/dX+0.4)+1
  cxI(2) = idint((y1-Y(1))/dY+0.4)+1
  cxI(3) = idint((z1-Z(1))/dZ+0.4)+1

  cxB = cxI - ORDN/2+1
  cxT = cxB + ORDN - 1
       
  cmin = 1
  cmax = ex
  if(Symmetry == OCTANT  .and.dabs(X(1))<dX) cmin(1) = -ORDN/2+1
  if(Symmetry == OCTANT  .and.dabs(Y(1))<dY) cmin(2) = -ORDN/2+1
  if(Symmetry /= NO_SYMM .and.dabs(Z(1))<dZ) cmin(3) = -ORDN/2+1

  do m =1,3
   if(cxB(m) < cmin(m))then
      cxB(m) = cmin(m)
      cxT(m) = cxB(m) + ORDN - 1
   endif
   if(cxT(m) > cmax(m))then
      cxT(m) = cmax(m)
      cxB(m) = cxT(m) + 1 - ORDN
   endif
 enddo
 if(cxB(1)>0)then
  cx(1) = (x1 - X(cxB(1)))/dX
 else
  cx(1) = (x1 + X(1-cxB(1)))/dX
 endif
 if(cxB(2)>0)then
  cx(2) = (y1 - Y(cxB(2)))/dY
 else
  cx(2) = (y1 + Y(1-cxB(2)))/dY
 endif
 if(cxB(3)>0)then
  cx(3) = (z1 - Z(cxB(3)))/dZ
 else
  cx(3) = (z1 + Z(1-cxB(3)))/dZ
 endif

  if(decide3d(ex,f,f,cxB,cxT,SoA,ya,ORDN,Symmetry))then
     write(*,*)"global_interp position: ",x1,y1,z1
     write(*,*)"data range: ",X(1),X(ex(1)),Y(1),Y(ex(2)),Z(1),Z(ex(3))
     stop
  endif
  call polin3(x1a,x1a,x1a,ya,cx(1),cx(2),cx(3),f_int,ddy,ORDN)

  return

  end subroutine global_interp
!----------------------------------------------------------------
! decide which 3d data to be used does not surport PI-Symmetry yet 
!----------------------------------------------------------------
  function decide3d(ex,f,fpi,cxB,cxT,SoA,ya,ORDN,Symmetry)  result(gont)
  implicit none

  integer,                                 intent(in) :: ORDN,Symmetry
  integer,dimension(1:3)                 , intent(in) :: ex,cxB,cxT
  real*8, dimension(1:3)                 , intent(in) :: SoA
  real*8, dimension(ex(1),ex(2),ex(3))   , intent(in) :: f,fpi
  real*8, dimension(cxB(1):cxT(1),cxB(2):cxT(2),cxB(3):cxT(3)), intent(out):: ya
  logical::gont

  integer,dimension(1:3) :: fmin1,fmin2,fmax1,fmax2
  integer::i,j,k,m

  gont=.false.
  do m=1,3
! check cxB and cxT are NaN or not  
    if(.not.(iabs(cxB(m)).ge.0)) gont=.true.
    if(.not.(iabs(cxT(m)).ge.0)) gont=.true.
    fmin1(m) = max(1,cxB(m))
    fmax1(m) = cxT(m)
    fmin2(m) = cxB(m)
    fmax2(m) = min(0,cxT(m))
    if((fmin1(m).le.fmax1(m)).and.(  fmin1(m)<1.or.  fmax1(m)>ex(m)))gont=.true.
    if((fmin2(m).le.fmax2(m)).and.(1-fmax2(m)<1.or.1-fmin2(m)>ex(m)))gont=.true.
  enddo
!sanity check
  if(gont)then
          write(*,*)"error in decide3d"
          write(*,*)((fmin1.le.fmax1).and.(  fmin1<1.or.  fmax1>ex))
          write(*,*)((fmin2.le.fmax2).and.(1-fmax2<1.or.1-fmin2>ex))
          write(*,*)"cxB, cxT and data shape:"
          write(*,*)cxB,cxT,ex
          write(*,*)"resulted fmin1, fmax1 and fmin2, fmax2:"
          write(*,*)fmin1,fmax1,fmin2,fmax2
  else

! Optimized for Intel ifx SIMD vectorization (AVX-512)
! Use explicit loops with OMP SIMD to enable ZMM register usage
! Main region: direct copy (continuous memory access)
  do k=fmin1(3),fmax1(3)
     do j=fmin1(2),fmax1(2)
!DIR$ VECTOR ALWAYS

        do i=fmin1(1),fmax1(1)
           ya(i,j,k) = f(i,j,k)
        enddo
     enddo
  enddo

! X-direction symmetry (cell uses 1-i)
  do k=fmin1(3),fmax1(3)
     do j=fmin1(2),fmax1(2)
!DIR$ VECTOR ALWAYS

        do i=fmin2(1),fmax2(1)
           ya(i,j,k) = f(1-i,j,k)*SoA(1)
        enddo
     enddo
  enddo

! Y-direction symmetry
  do k=fmin1(3),fmax1(3)
     do j=fmin2(2),fmax2(2)
!DIR$ VECTOR ALWAYS

        do i=fmin1(1),fmax1(1)
           ya(i,j,k) = f(i,1-j,k)*SoA(2)
        enddo
     enddo
  enddo

! XY-direction symmetry
  do k=fmin1(3),fmax1(3)
     do j=fmin2(2),fmax2(2)
!DIR$ VECTOR ALWAYS

        do i=fmin2(1),fmax2(1)
           ya(i,j,k) = f(1-i,1-j,k)*SoA(1)*SoA(2)
        enddo
     enddo
  enddo

! Z-direction symmetry
  do k=fmin2(3),fmax2(3)
     do j=fmin1(2),fmax1(2)
!DIR$ VECTOR ALWAYS

        do i=fmin1(1),fmax1(1)
           ya(i,j,k) = f(i,j,1-k)*SoA(3)
        enddo
     enddo
  enddo

! XZ-direction symmetry
  do k=fmin2(3),fmax2(3)
     do j=fmin1(2),fmax1(2)
!DIR$ VECTOR ALWAYS

        do i=fmin2(1),fmax2(1)
           ya(i,j,k) = f(1-i,j,1-k)*SoA(1)*SoA(3)
        enddo
     enddo
  enddo

! YZ-direction symmetry
  do k=fmin2(3),fmax2(3)
     do j=fmin2(2),fmax2(2)
!DIR$ VECTOR ALWAYS

        do i=fmin1(1),fmax1(1)
           ya(i,j,k) = f(i,1-j,1-k)*SoA(2)*SoA(3)
        enddo
     enddo
  enddo

! XYZ-direction symmetry
  do k=fmin2(3),fmax2(3)
     do j=fmin2(2),fmax2(2)
!DIR$ VECTOR ALWAYS

        do i=fmin2(1),fmax2(1)
           ya(i,j,k) = f(1-i,1-j,1-k)*SoA(1)*SoA(2)*SoA(3)
        enddo
     enddo
  enddo

  endif

  end function decide3d

#else
#error Not define Vertex nor Cell
#endif
#endif
! common code for cell and vertex
!------------------------------------------------------------------------------
! Lagrangian polynomial interpolation
!------------------------------------------------------------------------------

  subroutine polint(xa,ya,x,y,dy,ordn)

  implicit none

!~~~~~~> Input Parameter:
  integer,intent(in) :: ordn
  real*8, dimension(ordn), intent(in) :: xa,ya
  real*8, intent(in) :: x
  real*8, intent(out) :: y,dy

!~~~~~~> Other parameter:

  integer :: m,n,ns
  real*8, dimension(ordn) :: c,d,den,ho
  real*8 :: dif,dift

!~~~~~~>

  n=ordn
  m=ordn

  c=ya
  d=ya
  ho=xa-x

  ns=1
  dif=abs(x-xa(1))
  do m=1,n
   dift=abs(x-xa(m))
   if(dift < dif) then
    ns=m
    dif=dift
   end if
  end do

  y=ya(ns)
  ns=ns-1
  do m=1,n-1
    den(1:n-m)=ho(1:n-m)-ho(1+m:n)
    if (any(den(1:n-m) == 0.0))then
      write(*,*) 'failure in polint for point',x
      write(*,*) 'with input points: ',xa
      stop
    endif
    den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
    d(1:n-m)=ho(1+m:n)*den(1:n-m)
    c(1:n-m)=ho(1:n-m)*den(1:n-m)
    if (2*ns < n-m) then
      dy=c(ns+1)
    else
      dy=d(ns)
      ns=ns-1
    end if
    y=y+dy
  end do

  return

  end subroutine polint
!------------------------------------------------------------------------------
!
! interpolation in 2 dimensions, follow yx order
!
!------------------------------------------------------------------------------
  subroutine polin2(x1a,x2a,ya,x1,x2,y,dy,ordn)

  implicit none

!~~~~~~> Input parameters:
  integer,intent(in) :: ordn
  real*8, dimension(1:ordn), intent(in) :: x1a,x2a
  real*8, dimension(1:ordn,1:ordn), intent(in) :: ya
  real*8, intent(in) :: x1,x2
  real*8, intent(out) :: y,dy

!~~~~~~> Other parameters:

  integer  :: i,m
  real*8, dimension(ordn) :: ymtmp
  real*8, dimension(ordn) :: yntmp

  m=size(x1a)
  
  do i=1,m

    yntmp=ya(i,:)
    call polint(x2a,yntmp,x2,ymtmp(i),dy,ordn)

  end do

  call polint(x1a,ymtmp,x1,y,dy,ordn)

  return

  end subroutine polin2
!------------------------------------------------------------------------------
!
! interpolation in 3 dimensions, follow zyx order
!
!------------------------------------------------------------------------------
  subroutine polin3(x1a,x2a,x3a,ya,x1,x2,x3,y,dy,ordn)

  implicit none

!~~~~~~> Input parameters:
  integer,intent(in) :: ordn
  real*8, dimension(1:ordn), intent(in) :: x1a,x2a,x3a
  real*8, dimension(1:ordn,1:ordn,1:ordn), intent(in) :: ya
  real*8, intent(in) :: x1,x2,x3
  real*8, intent(out) :: y,dy

!~~~~~~> Other parameters:

  integer  :: i,j,m,n
  real*8, dimension(ordn,ordn) :: yatmp
  real*8, dimension(ordn) :: ymtmp
  real*8, dimension(ordn) :: yntmp
  real*8, dimension(ordn) :: yqtmp

  m=size(x1a)
  n=size(x2a)
  
  do i=1,m
   do j=1,n

    yqtmp=ya(i,j,:)
    call polint(x3a,yqtmp,x3,yatmp(i,j),dy,ordn)

   end do

    yntmp=yatmp(i,:)
    call polint(x2a,yntmp,x2,ymtmp(i),dy,ordn)

  end do

  call polint(x1a,ymtmp,x1,y,dy,ordn)

  return

  end subroutine polin3

  subroutine global_interpind(ex,X,Y,Z,f,f_int,x1,y1,z1,ORDN,SoA,symmetry,inds,coef,sst)
  implicit none

!~~~~~~> Input parameters:

  integer,                             intent(in) :: ex(1:3), symmetry,ORDN,sst
  real*8,intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: f
  real*8,                              intent(out):: f_int
  real*8,                              intent(in) :: x1,y1,z1
  real*8, dimension(3),                intent(in) :: SoA
  integer,dimension(3),                intent(in) :: inds
  real*8, dimension(3*ORDN),           intent(in) :: coef

!~~~~~~> Other parameters:

  real*8, dimension(-ORDN+1:ex(1)+ORDN,-ORDN+1:ex(2)+ORDN,-ORDN+1:ex(3)+ORDN) :: fh
  integer :: m
  integer,dimension(3) :: cxB,cxT
  real*8, dimension(ORDN,ORDN,ORDN) :: ya
  real*8, dimension(ORDN,ORDN) :: tmp2
  real*8, dimension(ORDN) :: tmp1
  real*8, dimension(3) :: SoAh

! +1 because c++ gives 0 for first point
  cxB = inds+1  
  cxT = cxB + ORDN - 1
      
  if(all(cxB>0).and.all(cxT<ex+1))then
     ya=f(cxB(1):cxT(1),cxB(2):cxT(2),cxB(3):cxT(3))
  elseif(any(cxB<-ORDN+1).or.any(cxT>ex+ORDN))then
     write(*,*)"error in global_interpind, cxB = ",cxB
     write(*,*)"                           cxT = ",cxT 
     write(*,*)"                           ext = ",ex
     stop
  else
     if(sst==-1)then
        SoAh = SoA
        if(any(cxT>ex)) write(*,*)"error global_interpind sst =",sst
     elseif(sst==0.or.sst==1)then
        SoAh = SoA
        SoAh(3) = 0
        if(cxB(3)<1.or.cxT(3)>ex(3)) write(*,*)"error global_interpind sst =",sst
     elseif(sst==2.or.sst==3)then
        SoAh(1) = SoA(2)
        SoAh(2) = SoA(3)
        SoAh(3) = 0
        if(cxB(3)<1.or.cxT(3)>ex(3)) write(*,*)"error global_interpind sst =",sst
     elseif(sst==4.or.sst==5)then
        SoAh(1) = SoA(1)
        SoAh(2) = SoA(3)
        SoAh(3) = 0
        if(cxB(3)<1.or.cxT(3)>ex(3)) write(*,*)"error global_interpind sst =",sst,cxB(3),cxT(3)
     endif
     call symmetry_tbd(ORDN,ex,f,fh,SoAh)
     ya=fh(cxB(1):cxT(1),cxB(2):cxT(2),cxB(3):cxT(3))
  endif 

  tmp2=0
  do m=1,ORDN
    tmp2 = tmp2 + coef(2*ORDN+m)*ya(:,:,m)
  enddo

  tmp1=0
  do m=1,ORDN
    tmp1 = tmp1 + coef(ORDN+m)*tmp2(:,m)
  enddo

  f_int=0
  do m=1,ORDN
    f_int = f_int + coef(m)*tmp1(m)
  enddo

  return

  end subroutine global_interpind
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!global interpolation with given index and coeffients
! special for shell to shell
  subroutine global_interpind2d(ex,X,Y,Z,f,f_int,x1,y1,z1,ORDN,SoA,symmetry,inds,coef,sst)
  implicit none

!~~~~~~> Input parameters:

  integer,                             intent(in) :: ex(1:3), symmetry,ORDN,sst
  real*8,intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: f
  real*8,                              intent(out):: f_int
  real*8,                              intent(in) :: x1,y1,z1
  real*8, dimension(3),                intent(in) :: SoA
  integer,dimension(3),                intent(in) :: inds
  real*8, dimension(2*ORDN),           intent(in) :: coef

!~~~~~~> Other parameters:

  real*8, dimension(-ORDN+1:ex(1)+ORDN,-ORDN+1:ex(2)+ORDN,ex(3)) :: fh
  integer :: m
  integer,dimension(2) :: cxB,cxT
  real*8, dimension(ORDN,ORDN) :: ya
  real*8, dimension(ORDN) :: tmp1
  real*8, dimension(2) :: SoAh

! +1 because c++ gives 0 for first point
  cxB = inds(1:2)+1  
  cxT = cxB + ORDN - 1
      
  if(all(cxB>0).and.all(cxT<ex(1:2)+1))then
     ya=f(cxB(1):cxT(1),cxB(2):cxT(2),inds(3))
  elseif(any(cxB<-ORDN+1).or.any(cxT>ex(1:2)+ORDN))then
     write(*,*)"error in global_interpind2d, cxB = ",cxB
     write(*,*)"                             cxT = ",cxT 
     write(*,*)"                             ext = ",ex(1:2)
     stop
  else
     if(sst==-1)then
       write(*,*)"error in global_interpind2d, sst = ",sst
       stop
     elseif(sst==0.or.sst==1)then
        SoAh = SoA(1:2)
     elseif(sst==2.or.sst==3)then
        SoAh(1) = SoA(2)
        SoAh(2) = SoA(3)
     elseif(sst==4.or.sst==5)then
        SoAh(1) = SoA(1)
        SoAh(2) = SoA(3)
     endif
     call symmetry_stbd(ORDN,ex,f,fh,SoAh)
     ya=fh(cxB(1):cxT(1),cxB(2):cxT(2),inds(3))
  endif 

  tmp1=0
  do m=1,ORDN
    tmp1 = tmp1 + coef(ORDN+m)*ya(:,m)
  enddo

  f_int=0
  do m=1,ORDN
    f_int = f_int + coef(m)*tmp1(m)
  enddo

  return

  end subroutine global_interpind2d
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!global interpolation with given index and coeffients
! special for shell to shell
! dumyd refer to source
  subroutine global_interpind1d(ex,X,Y,Z,f,f_int,x1,y1,z1,ORDN,SoA,symmetry,indsi,coef,sst,dumyd)
  implicit none

!~~~~~~> Input parameters:

  integer,                             intent(in) :: ex(1:3),symmetry,ORDN,sst,dumyd
  real*8,intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: f
  real*8,                              intent(out):: f_int
  real*8,                              intent(in) :: x1,y1,z1
  real*8, dimension(3),                intent(in) :: SoA
  integer,dimension(3),                intent(in) :: indsi
  real*8, dimension(ORDN),             intent(in) :: coef

!~~~~~~> Other parameters:

  real*8, dimension(-ORDN+1:ex(1)+ORDN,-ORDN+1:ex(2)+ORDN,ex(3)) :: fh
  integer :: m
  integer :: cxB,cxT
  real*8, dimension(ORDN) :: ya
  real*8 :: SoAh
  integer,dimension(3) :: inds

! +1 because c++ gives 0 for first point
  inds = indsi + 1
  cxB = inds(1)  
  cxT = cxB + ORDN - 1
     
! active is rho  
  if(dumyd==1)then

  if(cxB>0.and.cxT<ex(1)+1)then
     ya=f(cxB:cxT,inds(2),inds(3))
  elseif(cxB<-ORDN+1.or.cxT>ex(1)+ORDN)then
     write(*,*)"error in global_interpind1d, cxB = ",cxB
     write(*,*)"                             cxT = ",cxT 
     write(*,*)"                             ext = ",ex(1)
     stop
  else
     if(sst==-1)then
       write(*,*)"error in global_interpind1d, sst = ",sst
       stop
     elseif(sst==0.or.sst==1)then
        SoAh = SoA(1)
     elseif(sst==2.or.sst==3)then
        SoAh = SoA(2)
     elseif(sst==4.or.sst==5)then
        SoAh = SoA(1)
     endif
     call symmetry_sntbd(ORDN,ex,f,fh,SoAh,1-dumyd)
     ya=fh(cxB:cxT,inds(2),inds(3))
  endif 

! active is sigma
  elseif(dumyd==0)then

  if(cxB>0.and.cxT<ex(2)+1)then
     ya=f(inds(2),cxB:cxT,inds(3))
  elseif(cxB<-ORDN+1.or.cxT>ex(2)+ORDN)then
     write(*,*)"error in global_interpind1d, cxB = ",cxB
     write(*,*)"                             cxT = ",cxT 
     write(*,*)"                             ext = ",ex(2)
     stop
  else
     if(sst==-1)then
       write(*,*)"error in global_interpind1d, sst = ",sst
       stop
     elseif(sst==0.or.sst==1)then
        SoAh = SoA(2)
     elseif(sst==2.or.sst==3)then
        SoAh = SoA(3)
     elseif(sst==4.or.sst==5)then
        SoAh = SoA(3)
     endif
     call symmetry_sntbd(ORDN,ex,f,fh,SoAh,1-dumyd)
     ya=fh(inds(2),cxB:cxT,inds(3))
  endif 

  else
          write(*,*)"error in global_interpind1d, not recognized dumyd = ",dumyd
  endif

  f_int=0
  do m=1,ORDN
    f_int = f_int + coef(m)*ya(m)
  enddo

  return

  end subroutine global_interpind1d
!-----------------------------------------------------------------------------------------------------------------
! three dimensional interpolation for both vertex and cell center grid structure  
! for distinguishing shell and Cartesian
  subroutine global_interp_ss(ex,X,Y,Z,f,f_int,x1,y1,z1,ORDN,SoA,symmetry,sst)
  implicit none

!~~~~~~> Input parameters:

  integer,                             intent(in) :: ex(1:3), symmetry,ORDN,sst
  real*8,intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: f
  real*8,                              intent(out):: f_int
  real*8,                              intent(in) :: x1,y1,z1
  real*8, dimension(3),                intent(in) :: SoA

!~~~~~~> Other parameters:

  real*8, dimension(-ORDN+1:ex(1)+ORDN,-ORDN+1:ex(2)+ORDN,-ORDN+1:ex(3)+ORDN) :: fh
  real*8, dimension(3) :: SoAh
  integer :: j,m,imin,jmin,kmin
  integer,dimension(3) :: cxB,cxT,cxI,cmin,cmax
  real*8,dimension(3) :: cx
  real*8, dimension(1:ORDN) :: x1a
  real*8, dimension(1:ORDN,1:ORDN,1:ORDN) :: ya
  integer, parameter :: NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2
  real*8 :: dX,dY,dZ,ddy
  real*8, parameter :: ONE=1.d0

  imin = lbound(f,1)
  jmin = lbound(f,2)
  kmin = lbound(f,3)

  dX = X(imin+1)-X(imin)
  dY = Y(jmin+1)-Y(jmin)
  dZ = Z(kmin+1)-Z(kmin)

  forall( j = 1:ordn ) x1a(j) = ( j - 1 )* ONE

  cxI(1) = idint((x1-X(1))/dX+0.4)+1
  cxI(2) = idint((y1-Y(1))/dY+0.4)+1
  cxI(3) = idint((z1-Z(1))/dZ+0.4)+1

  cxB = cxI - ORDN/2+1
  cxT = cxB + ORDN - 1

  cmin = 1
  cmax = ex

  if(sst==-1)then
     SoAh = SoA
     cmin = -ORDN+1
  elseif(sst==0.or.sst==1)then
     SoAh = SoA
     SoAh(3) = 0
     cmin(1:2) = -ORDN+1
     cmax(1:2) = ex(1:2)+ORDN
  elseif(sst==2.or.sst==3)then
     SoAh(1) = SoA(2)
     SoAh(2) = SoA(3)
     SoAh(3) = 0
     cmin(1:2) = -ORDN+1
     cmax(1:2) = ex(1:2)+ORDN
  elseif(sst==4.or.sst==5)then
     SoAh(1) = SoA(1)
     SoAh(2) = SoA(3)
     SoAh(3) = 0
     cmin(1:2) = -ORDN+1
     cmax(1:2) = ex(1:2)+ORDN
  endif
  do m =1,3
   if(cxB(m) < cmin(m))then
      cxB(m) = cmin(m)
      cxT(m) = cxB(m) + ORDN - 1
   endif
   if(cxT(m) > cmax(m))then
      cxT(m) = cmax(m)
      cxB(m) = cxT(m) + 1 - ORDN
   endif
 enddo
 cx(1) = (x1 - X(1))/dX-cxB(1)+1
 cx(2) = (y1 - Y(1))/dY-cxB(2)+1
 cx(3) = (z1 - Z(1))/dZ-cxB(3)+1

  call symmetry_tbd(ORDN,ex,f,fh,SoAh)
  ya=fh(cxB(1):cxT(1),cxB(2):cxT(2),cxB(3):cxT(3))

  call polin3(x1a,x1a,x1a,ya,cx(1),cx(2),cx(3),f_int,ddy,ORDN)

  return

  end subroutine global_interp_ss
!-----------------------------------------------------------------------------------------------------------------
! two dimensional interpolation for both vertex and cell center grid structure  
! for distinguishing shell and Cartesian
  subroutine global_interp_ss_2d(ex,X,Y,indZ,f,f_int,x1,y1,ORDN,SoA,symmetry,sst)
  implicit none

!~~~~~~> Input parameters:

  integer,                             intent(in) :: ex(1:3),indZ,symmetry,ORDN,sst
  real*8,intent(in) :: X(ex(1)),Y(ex(2))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: f
  real*8,                              intent(out):: f_int
  real*8,                              intent(in) :: x1,y1
  real*8, dimension(3),                intent(in) :: SoA

!~~~~~~> Other parameters:

  real*8, dimension(-ORDN+1:ex(1)+ORDN,-ORDN+1:ex(2)+ORDN,-ORDN+1:ex(3)+ORDN) :: fh
  real*8, dimension(3) :: SoAh
  integer :: j,m,imin,jmin,kmin
  integer,dimension(2) :: cxB,cxT,cxI,cmin,cmax
  real*8,dimension(2) :: cx
  real*8, dimension(1:ORDN) :: x1a
  real*8, dimension(1:ORDN,1:ORDN) :: ya
  integer, parameter :: NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2
  real*8 :: dX,dY,ddy
  real*8, parameter :: ONE=1.d0

! sanity check
  if(indZ < 1 .or. indZ > ex(3))then
      write(*,*)"error in global_interp_ss_2d, ext = ",ex(3),"ind = ",indZ
      return
  endif

  imin = lbound(f,1)
  jmin = lbound(f,2)
  kmin = lbound(f,3)

  dX = X(imin+1)-X(imin)
  dY = Y(jmin+1)-Y(jmin)

  forall( j = 1:ordn ) x1a(j) = ( j - 1 )* ONE

  cxI(1) = idint((x1-X(1))/dX+0.4)+1
  cxI(2) = idint((y1-Y(1))/dY+0.4)+1

  cxB = cxI - ORDN/2+1
  cxT = cxB + ORDN - 1

  cmin = 1
  cmax = ex(1:2)

  if(sst==-1)then
     SoAh = SoA
     cmin = -ORDN+1
  elseif(sst==0.or.sst==1)then
     SoAh = SoA
     SoAh(3) = 0
     cmin(1:2) = -ORDN+1
     cmax(1:2) = ex(1:2)+ORDN
  elseif(sst==2.or.sst==3)then
     SoAh(1) = SoA(2)
     SoAh(2) = SoA(3)
     SoAh(3) = 0
     cmin(1:2) = -ORDN+1
     cmax(1:2) = ex(1:2)+ORDN
  elseif(sst==4.or.sst==5)then
     SoAh(1) = SoA(1)
     SoAh(2) = SoA(3)
     SoAh(3) = 0
     cmin(1:2) = -ORDN+1
     cmax(1:2) = ex(1:2)+ORDN
  endif
  do m =1,2
   if(cxB(m) < cmin(m))then
      cxB(m) = cmin(m)
      cxT(m) = cxB(m) + ORDN - 1
   endif
   if(cxT(m) > cmax(m))then
      cxT(m) = cmax(m)
      cxB(m) = cxT(m) + 1 - ORDN
   endif
 enddo
 cx(1) = (x1 - X(1))/dX-cxB(1)+1
 cx(2) = (y1 - Y(1))/dY-cxB(2)+1

  call symmetry_tbd(ORDN,ex,f,fh,SoAh)
  ya=fh(cxB(1):cxT(1),cxB(2):cxT(2),indZ)

  call polin2(x1a,x1a,ya,cx(1),cx(2),f_int,ddy,ORDN)

  return

  end subroutine global_interp_ss_2d
