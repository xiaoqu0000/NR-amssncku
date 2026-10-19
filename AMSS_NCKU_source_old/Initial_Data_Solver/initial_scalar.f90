
!-----------------------------------------------------------------------------
!
! Setting initial scalar with spherical Gauss profile centered at shell r=R0
! with width WD and amplitude A
!
!-----------------------------------------------------------------------------

  subroutine get_initial_scalar(ex, X, Y, Z,Sphi,Spi,R0,WD,A)
  implicit none

!~~~~~~> Input parameters

  integer,intent(in ):: ex(1:3)
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3)),R0,WD,A
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out):: Sphi,Spi

!~~~~~~> Local variables

  real*8 :: rr
  integer::i,j,k
  real*8, parameter :: ZEO = 0.d0,TWO=2.d0

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
     rr = dsqrt(X(i)*X(i)+Y(j)*Y(j)+Z(k)*Z(k))-R0
     Sphi(i,j,k) = A*dexp(-rr*rr/TWO/WD/WD)
  enddo
  enddo
  enddo

  Spi = ZEO

  return

  end subroutine get_initial_scalar
! for shell
  subroutine get_initial_scalar_sh(ex, X, Y, Z,Sphi,Spi,R0,WD,A)
  implicit none

!~~~~~~> Input parameters

  integer,intent(in ):: ex(1:3)
  real*8, intent(in ):: R0,WD,A
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in):: X, Y, Z
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out):: Sphi,Spi

!~~~~~~> Local variables

  real*8 :: rr
  integer::i,j,k
  real*8, parameter :: ZEO = 0.d0,TWO=2.d0

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
     rr = dsqrt(X(i,j,k)*X(i,j,k)+Y(i,j,k)*Y(i,j,k)+Z(i,j,k)*Z(i,j,k))-R0
     Sphi(i,j,k) = A*dexp(-rr*rr/TWO/WD/WD)
  enddo
  enddo
  enddo

  Spi = ZEO

  return

  end subroutine get_initial_scalar_sh
