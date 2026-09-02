
!-----------------------------------------------------------------------------
!
! remove the trace of Aij
! trace-free Aij and enforce the determinant of bssn metric to one
!-----------------------------------------------------------------------------

  subroutine enforce_ag(ex,  dxx,  gxy,  gxz,  dyy,  gyz,  dzz, &
                             Axx,  Axy,  Axz,  Ayy,  Ayz,  Azz)
  implicit none

!~~~~~~> Input parameters:

  integer,                              intent(in)    :: ex(1:3)
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: dxx,dyy,dzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: gxy,gxz,gyz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: Axx,Axy,Axz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: Ayy,Ayz,Azz

!~~~~~~~> Local variable:
  
  real*8, dimension(ex(1),ex(2),ex(3)) :: trA,detg
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gyy,gzz 
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupxx,gupxy,gupxz,gupyy,gupyz,gupzz
  real*8, parameter :: F1o3 = 1.D0 / 3.D0, ONE = 1.D0, TWO = 2.D0

!~~~~~~>

  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE

  detg =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz - &
          gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz
  gupxx =   ( gyy * gzz - gyz * gyz ) / detg
  gupxy = - ( gxy * gzz - gyz * gxz ) / detg
  gupxz =   ( gxy * gyz - gyy * gxz ) / detg
  gupyy =   ( gxx * gzz - gxz * gxz ) / detg
  gupyz = - ( gxx * gyz - gxy * gxz ) / detg
  gupzz =   ( gxx * gyy - gxy * gxy ) / detg

  trA =         gupxx * Axx + gupyy * Ayy + gupzz * Azz &
       + TWO * (gupxy * Axy + gupxz * Axz + gupyz * Ayz)

  Axx = Axx - F1o3 * gxx * trA
  Axy = Axy - F1o3 * gxy * trA
  Axz = Axz - F1o3 * gxz * trA
  Ayy = Ayy - F1o3 * gyy * trA
  Ayz = Ayz - F1o3 * gyz * trA
  Azz = Azz - F1o3 * gzz * trA

  detg = ONE / ( detg ** F1o3 ) 
  
  gxx = gxx * detg
  gxy = gxy * detg
  gxz = gxz * detg
  gyy = gyy * detg
  gyz = gyz * detg
  gzz = gzz * detg

  dxx = gxx - ONE
  dyy = gyy - ONE
  dzz = gzz - ONE

  return

  end subroutine enforce_ag
!----------------------------------------------------------------------------------  
! swap the turn of a and g
!----------------------------------------------------------------------------------
  subroutine enforce_ga(ex,  dxx,  gxy,  gxz,  dyy,  gyz,  dzz, &
                             Axx,  Axy,  Axz,  Ayy,  Ayz,  Azz)
  use sft_trace_mod
  use omp_lib
  implicit none
  integer(8) :: ts0_thr
  integer :: tid

!~~~~~~> Input parameters:

  integer,                              intent(in)    :: ex(1:3)
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: dxx,dyy,dzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: gxy,gxz,gyz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: Axx,Axy,Axz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: Ayy,Ayz,Azz

!~~~~~~~> Local variable:

  real*8, parameter :: F1o3 = 1.D0 / 3.D0, ONE = 1.D0, TWO = 2.D0
  real*8  :: trA, detg
  real*8  :: gxx, gyy, gzz
  real*8  :: gupxx, gupxy, gupxz, gupyy, gupyz, gupzz
  integer :: i, j, k

!~~~~~~>

!$OMP PARALLEL PRIVATE(i,j,k,gxx,gyy,gzz,gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,detg,trA,ts0_thr,tid) SHARED(dxx,dyy,dzz,gxy,gxz,gyz,Axx,Axy,Axz,Ayy,Ayz,Azz,ex)
  tid = omp_get_thread_num()
  ts0_thr = sft_get_ts()
!$OMP DO COLLAPSE(1)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)

        gxx = dxx(i,j,k) + ONE
        gyy = dyy(i,j,k) + ONE
        gzz = dzz(i,j,k) + ONE
! for g
        detg =  gxx * gyy * gzz + gxy(i,j,k) * gyz(i,j,k) * gxz(i,j,k) &
              + gxz(i,j,k) * gxy(i,j,k) * gyz(i,j,k) &
              - gxz(i,j,k) * gyy * gxz(i,j,k) &
              - gxy(i,j,k) * gxy(i,j,k) * gzz &
              - gxx * gyz(i,j,k) * gyz(i,j,k)

        detg = ONE / ( detg ** F1o3 )

        gxx           = gxx           * detg
        gxy(i,j,k)    = gxy(i,j,k)    * detg
        gxz(i,j,k)    = gxz(i,j,k)    * detg
        gyy           = gyy           * detg
        gyz(i,j,k)    = gyz(i,j,k)    * detg
        gzz           = gzz           * detg

        dxx(i,j,k) = gxx - ONE
        dyy(i,j,k) = gyy - ONE
        dzz(i,j,k) = gzz - ONE
! for A
        gupxx =   ( gyy * gzz - gyz(i,j,k) * gyz(i,j,k) )
        gupxy = - ( gxy(i,j,k) * gzz - gyz(i,j,k) * gxz(i,j,k) )
        gupxz =   ( gxy(i,j,k) * gyz(i,j,k) - gyy * gxz(i,j,k) )
        gupyy =   ( gxx * gzz - gxz(i,j,k) * gxz(i,j,k) )
        gupyz = - ( gxx * gyz(i,j,k) - gxy(i,j,k) * gxz(i,j,k) )
        gupzz =   ( gxx * gyy - gxy(i,j,k) * gxy(i,j,k) )

        trA =         gupxx * Axx(i,j,k) + gupyy * Ayy(i,j,k) + gupzz * Azz(i,j,k) &
             + TWO * (gupxy * Axy(i,j,k) + gupxz * Axz(i,j,k) + gupyz * Ayz(i,j,k))

        Axx(i,j,k) = Axx(i,j,k) - F1o3 * gxx           * trA
        Axy(i,j,k) = Axy(i,j,k) - F1o3 * gxy(i,j,k)    * trA
        Axz(i,j,k) = Axz(i,j,k) - F1o3 * gxz(i,j,k)    * trA
        Ayy(i,j,k) = Ayy(i,j,k) - F1o3 * gyy           * trA
        Ayz(i,j,k) = Ayz(i,j,k) - F1o3 * gyz(i,j,k)    * trA
        Azz(i,j,k) = Azz(i,j,k) - F1o3 * gzz           * trA

      end do
    end do
  end do
!$OMP END DO
  call sft_trace("enforce_ga", ts0_thr, sft_get_ts(), 0, tid)
!$OMP END PARALLEL

  return

  end subroutine enforce_ga
