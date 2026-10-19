
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
#if 1 
!----------------------------------------------------------------------------------  
! swap the turn of a and g
!----------------------------------------------------------------------------------
  subroutine enforce_ga(ex,  dxx,  gxy,  gxz,  dyy,  gyz,  dzz, &
                             Axx,  Axy,  Axz,  Ayy,  Ayz,  Azz)
  implicit none

!~~~~~~> Input parameters:

  integer,                              intent(in)    :: ex(1:3)
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: dxx,dyy,dzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: gxy,gxz,gyz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: Axx,Axy,Axz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: Ayy,Ayz,Azz

!~~~~~~~> Local variable:
  
  real*8, dimension(ex(1),ex(2),ex(3)) :: trA
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gyy,gzz 
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupxx,gupxy,gupxz,gupyy,gupyz,gupzz
  real*8, parameter :: F1o3 = 1.D0 / 3.D0, ONE = 1.D0, TWO = 2.D0

!~~~~~~>

  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE
! for g
  gupzz =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz - &
           gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz

  gupzz = ONE / ( gupzz ** F1o3 ) 
  
  gxx = gxx * gupzz
  gxy = gxy * gupzz
  gxz = gxz * gupzz
  gyy = gyy * gupzz
  gyz = gyz * gupzz
  gzz = gzz * gupzz

  dxx = gxx - ONE
  dyy = gyy - ONE
  dzz = gzz - ONE
! for A  

  gupxx =   ( gyy * gzz - gyz * gyz )
  gupxy = - ( gxy * gzz - gyz * gxz )
  gupxz =   ( gxy * gyz - gyy * gxz )
  gupyy =   ( gxx * gzz - gxz * gxz )
  gupyz = - ( gxx * gyz - gxy * gxz )
  gupzz =   ( gxx * gyy - gxy * gxy )

  trA =         gupxx * Axx + gupyy * Ayy + gupzz * Azz &
       + TWO * (gupxy * Axy + gupxz * Axz + gupyz * Ayz)

  Axx = Axx - F1o3 * gxx * trA
  Axy = Axy - F1o3 * gxy * trA
  Axz = Axz - F1o3 * gxz * trA
  Ayy = Ayy - F1o3 * gyy * trA
  Ayz = Ayz - F1o3 * gyz * trA
  Azz = Azz - F1o3 * gzz * trA

  return

  end subroutine enforce_ga
#else
!----------------------------------------------------------------------------------  
! duplicate bam
!----------------------------------------------------------------------------------
  subroutine enforce_ga(ex,  dxx,  gxy,  gxz,  dyy,  gyz,  dzz, &
                             Axx,  Axy,  Axz,  Ayy,  Ayz,  Azz)
  implicit none

!~~~~~~> Input parameters:

  integer,                              intent(in)    :: ex(1:3)
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: dxx,dyy,dzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: gxy,gxz,gyz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: Axx,Axy,Axz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: Ayy,Ayz,Azz

!~~~~~~~> Local variable:
  
  real*8, dimension(ex(1),ex(2),ex(3)) :: trA
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gyy,gzz 
  real*8, dimension(ex(1),ex(2),ex(3)) :: aux,detginv
  real*8, parameter :: oot = 1.D0 / 3.D0, ONE = 1.D0, TWO = 2.D0

!~~~~~~>

  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE
! for g
aux = (2.d0*gxy*gxz*gyz + gxx*gyy*gzz &
    - gzz*gxy**2 - gyy*gxz**2 - gxx*gyz**2)**(-oot)

  gxx = gxx * aux
  gxy = gxy * aux
  gxz = gxz * aux
  gyy = gyy * aux
  gyz = gyz * aux
  gzz = gzz * aux

  dxx = gxx - ONE
  dyy = gyy - ONE
  dzz = gzz - ONE
! for A  

detginv = 1/(2.d0*gxy*gxz*gyz + gxx*gyy*gzz &
    - gzz*gxy**2 - gyy*gxz**2 - gxx*gyz**2)

trA = detginv*(-2.d0*Ayz*gxx*gyz + Axx*gyy*gzz + &
    gxx*(Azz*gyy + Ayy*gzz) + 2.d0*(gxz*(Ayz*gxy - Axz*gyy + &
    Axy*gyz) + gxy*(Axz*gyz - Axy*gzz)) - Azz*gxy**2 - Ayy*gxz**2 - &
    Axx*gyz**2)

aux = -(oot*trA)

  Axx = Axx + aux * gxx
  Axy = Axy + aux * gxy
  Axz = Axz + aux * gxz
  Ayy = Ayy + aux * gyy
  Ayz = Ayz + aux * gyz
  Azz = Azz + aux * gzz

  return

  end subroutine enforce_ga
#endif
