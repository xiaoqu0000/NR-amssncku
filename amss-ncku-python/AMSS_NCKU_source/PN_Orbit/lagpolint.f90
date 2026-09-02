!$Id: lagpolint.f90,v 1.1.1.1 2017/09/25 02:35:56 zjcao Exp $
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
