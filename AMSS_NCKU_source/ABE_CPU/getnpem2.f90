

#include "macrodef.fh"

!-----------------------------------------------------------------------------
!
! compute the Newman-Penrose Weyl scalar Psi4
! for BSSN dynamical variables
!
!-----------------------------------------------------------------------------

  subroutine getnpem2(ext, X, Y, Z,         &
               chi,dxx,gxy,gxz,dyy,gyz,dzz, &
               Ex,Ey,Ez,Bx,By,Bz,           &
               Rphi2, Iphi2, &
               symmetry)

  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ext(1:3),symmetry
  real*8, intent(in ):: X(1:ext(1)),Y(1:ext(2)),Z(1:ext(3))
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in ) :: dxx,gxy,gxz,dyy,gyz,dzz
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in ) :: chi,Ex,Ey,Ez,Bx,By,Bz
  real*8, dimension(ext(1),ext(2),ext(3)), intent(out):: Rphi2,Iphi2

!~~~~~~> Other variables:

  real*8, dimension(ext(1),ext(2),ext(3)) :: f,fx,fy,fz
  real*8, dimension(ext(1),ext(2),ext(3)) :: gxx,gyy,gzz
  real*8, dimension(ext(1),ext(2),ext(3)) :: chipn1,chi3o2
  real*8, dimension(ext(1),ext(2),ext(3)) :: vx,vy,vz,ux,uy,uz,wx,wy,wz
  real*8, dimension(ext(1),ext(2),ext(3)) :: HEx,HEy,HEz,HBx,HBy,HBz
  real*8, dimension(ext(1),ext(2),ext(3)) :: gupxx,gupxy,gupxz
  real*8, dimension(ext(1),ext(2),ext(3)) :: gupyy,gupyz,gupzz

  real*8, parameter :: ZEO = 0.d0, ONE = 1.d0, TWO = 2.d0
  real*8, parameter :: F1o3 = 1.d0/3.d0, FOUR = 4.d0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  integer::i,j,k
  real*8,parameter::TINYRR=1.d-14

  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE
  chipn1= chi + ONE
  chi3o2  = dsqrt(chipn1)**3

! invert tilted metric
  gupzz =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz - &
           gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz
  gupxx =   ( gyy * gzz - gyz * gyz ) / gupzz
  gupxy = - ( gxy * gzz - gyz * gxz ) / gupzz
  gupxz =   ( gxy * gyz - gyy * gxz ) / gupzz
  gupyy =   ( gxx * gzz - gxz * gxz ) / gupzz
  gupyz = - ( gxx * gyz - gxy * gxz ) / gupzz
  gupzz =   ( gxx * gyy - gxy * gxy ) / gupzz

! initialize U, V, W vetors	
! v:r; u: phi; w: theta
#if (tetradtype == 0)
  do i=1,ext(1)
  do j=1,ext(2)
  do k=1,ext(3)
     if(abs(X(i)) < TINYRR .and. abs(Y(j)) < TINYRR .and. abs(Z(k)) < TINYRR)then
        vx(i,j,k) = TINYRR
        vy(i,j,k) = TINYRR
        vz(i,j,k) = TINYRR
     else
        vx(i,j,k) = X(i)
        vy(i,j,k) = Y(j)
        vz(i,j,k) = Z(k)
     endif
     if(abs(X(i)) < TINYRR .and. abs(Y(j)) < TINYRR)then
        ux(i,j,k) = - TINYRR
        uy(i,j,k) = TINYRR
        uz(i,j,k) = ZEO
        wx(i,j,k) = TINYRR*Z(k)
        wy(i,j,k) = TINYRR*Z(k)
        wz(i,j,k) = -2*TINYRR*TINYRR
     else
        ux(i,j,k) = - Y(j)
        uy(i,j,k) = X(i)
        uz(i,j,k) = ZEO
        wx(i,j,k) = X(i)*Z(k)
        wy(i,j,k) = Y(j)*Z(k)
        wz(i,j,k) = -(X(i)*X(i) + Y(j)*Y(j))
     endif
  enddo
  enddo
  enddo

  f = 1.d0/chipn1

  fx = gxx*vx*vx + gyy*vy*vy + gzz*vz*vz &
     +(gxy*vx*vy + gxz*vx*vz + gyz*vy*vz)*TWO
  fx = dsqrt(fx*f)
  vx = vx/fx
  vy = vy/fx
  vz = vz/fx

  fx = gxx*vx*ux + gxy*vx*uy + gxz*vx*uz + &
       gxy*vy*ux + gyy*vy*uy + gyz*vy*uz + &
       gxz*vz*ux + gyz*vz*uy + gzz*vz*uz
  fx = fx*f
  ux = ux - fx*vx
  uy = uy - fx*vy
  uz = uz - fx*vz
  fx = gxx*ux*ux + gyy*uy*uy + gzz*uz*uz &
     +(gxy*ux*uy + gxz*ux*uz + gyz*uy*uz)*TWO
  fx = dsqrt(fx*f) 
  ux = ux/fx
  uy = uy/fx
  uz = uz/fx

  fx = gxx*vx*wx + gxy*vx*wy + gxz*vx*wz + &
       gxy*vy*wx + gyy*vy*wy + gyz*vy*wz + &
       gxz*vz*wx + gyz*vz*wy + gzz*vz*wz
  fx = fx*f       
  wx = wx - fx*vx
  wy = wy - fx*vy
  wz = wz - fx*vz
  fx = gxx*ux*wx + gxy*ux*wy + gxz*ux*wz + &
       gxy*uy*wx + gyy*uy*wy + gyz*uy*wz + &
       gxz*uz*wx + gyz*uz*wy + gzz*uz*wz
  fx = fx*f       
  wx = wx - fx*ux
  wy = wy - fx*uy
  wz = wz - fx*uz
  fx = gxx*wx*wx + gyy*wy*wy + gzz*wz*wz &
     +(gxy*wx*wy + gxz*wx*wz + gyz*wy*wz)*TWO
  fx = dsqrt(fx*f)
  wx = wx/fx
  wy = wy/fx
  wz = wz/fx
#elif  (tetradtype == 1)
  do i=1,ext(1)
  do j=1,ext(2)
  do k=1,ext(3)
     if(abs(X(i)) < TINYRR .and. abs(Y(j)) < TINYRR .and. abs(Z(k)) < TINYRR)then
        vx(i,j,k) = TINYRR
        vy(i,j,k) = TINYRR
        vz(i,j,k) = TINYRR
     else
        vx(i,j,k) = X(i)
        vy(i,j,k) = Y(j)
        vz(i,j,k) = Z(k)
     endif
     if(abs(X(i)) < TINYRR .and. abs(Y(j)) < TINYRR)then
        ux(i,j,k) = - TINYRR
        uy(i,j,k) = TINYRR
        uz(i,j,k) = ZEO
        wx(i,j,k) = TINYRR*Z(k)
        wy(i,j,k) = TINYRR*Z(k)
        wz(i,j,k) = -2*TINYRR*TINYRR
     else
        ux(i,j,k) = - Y(j)
        uy(i,j,k) = X(i)
        uz(i,j,k) = ZEO
        wx(i,j,k) = X(i)*Z(k)
        wy(i,j,k) = Y(j)*Z(k)
        wz(i,j,k) = -(X(i)*X(i) + Y(j)*Y(j))
     endif
  enddo
  enddo
  enddo

  f = 1.d0/chipn1

  fx = gxx*wx*wx + gyy*wy*wy + gzz*wz*wz &
     +(gxy*wx*wy + gxz*wx*wz + gyz*wy*wz)*TWO
  fx = dsqrt(fx*f)
  wx = wx/fx
  wy = wy/fx
  wz = wz/fx

  fx = gxx*wx*ux + gxy*wx*uy + gxz*wx*uz + &
       gxy*wy*ux + gyy*wy*uy + gyz*wy*uz + &
       gxz*wz*ux + gyz*wz*uy + gzz*wz*uz
  fx = fx*f
  ux = ux - fx*wx
  uy = uy - fx*wy
  uz = uz - fx*wz
  fx = gxx*ux*ux + gyy*uy*uy + gzz*uz*uz &
     +(gxy*ux*uy + gxz*ux*uz + gyz*uy*uz)*TWO
  fx = dsqrt(fx*f) 
  ux = ux/fx
  uy = uy/fx
  uz = uz/fx

  fx = gxx*vx*wx + gxy*vx*wy + gxz*vx*wz + &
       gxy*vy*wx + gyy*vy*wy + gyz*vy*wz + &
       gxz*vz*wx + gyz*vz*wy + gzz*vz*wz
  fx = fx*f       
  vx = vx - fx*wx
  vy = vy - fx*wy
  vz = vz - fx*wz
  fx = gxx*ux*vx + gxy*ux*vy + gxz*ux*vz + &
       gxy*uy*vx + gyy*uy*vy + gyz*uy*vz + &
       gxz*uz*vx + gyz*uz*vy + gzz*uz*vz
  fx = fx*f       
  vx = vx - fx*ux
  vy = vy - fx*uy
  vz = vz - fx*uz
  fx = gxx*vx*vx + gyy*vy*vy + gzz*vz*vz &
     +(gxy*vx*vy + gxz*vx*vz + gyz*vy*vz)*TWO
  fx = dsqrt(fx*f)
  vx = vx/fx
  vy = vy/fx
  vz = vz/fx
#elif (tetradtype == 2)  
  do i=1,ext(1)
  do j=1,ext(2)
  do k=1,ext(3)
     if(abs(X(i)) < TINYRR .and. abs(Y(j)) < TINYRR .and. abs(Z(k)) < TINYRR)then
        vx(i,j,k) = TINYRR
        vy(i,j,k) = TINYRR
        vz(i,j,k) = TINYRR
     else
        vx(i,j,k) = X(i)
        vy(i,j,k) = Y(j)
        vz(i,j,k) = Z(k)
     endif
     if(abs(X(i)) < TINYRR .and. abs(Y(j)) < TINYRR)then
        ux(i,j,k) = - TINYRR
        uy(i,j,k) = TINYRR
        uz(i,j,k) = ZEO
        wx(i,j,k) = TINYRR*Z(k)
        wy(i,j,k) = TINYRR*Z(k)
        wz(i,j,k) = -2*TINYRR*TINYRR
     else
        ux(i,j,k) = - Y(j)
        uy(i,j,k) = X(i)
        uz(i,j,k) = ZEO
        wx(i,j,k) = X(i)*Z(k)
        wy(i,j,k) = Y(j)*Z(k)
        wz(i,j,k) = -(X(i)*X(i) + Y(j)*Y(j))
     endif
  enddo
  enddo
  enddo

  fx = vx
  fy = vy
  fz = vz
  vx = gupxx*fx + gupxy*fy + gupxz*fz
  vy = gupxy*fx + gupyy*fy + gupyz*fz
  vz = gupxz*fx + gupyz*fy + gupzz*fz

  f = 1.d0/chipn1

  fx = gxx*vx*vx + gyy*vy*vy + gzz*vz*vz &
     +(gxy*vx*vy + gxz*vx*vz + gyz*vy*vz)*TWO
  fx = dsqrt(fx*f)
  vx = vx/fx
  vy = vy/fx
  vz = vz/fx

  fx = gxx*vx*ux + gxy*vx*uy + gxz*vx*uz + &
       gxy*vy*ux + gyy*vy*uy + gyz*vy*uz + &
       gxz*vz*ux + gyz*vz*uy + gzz*vz*uz
  fx = fx*f
  ux = ux - fx*vx
  uy = uy - fx*vy
  uz = uz - fx*vz
  fx = gxx*ux*ux + gyy*uy*uy + gzz*uz*uz &
     +(gxy*ux*uy + gxz*ux*uz + gyz*uy*uz)*TWO
  fx = dsqrt(fx*f) 
  ux = ux/fx
  uy = uy/fx
  uz = uz/fx

  fx = gxx*vx*wx + gxy*vx*wy + gxz*vx*wz + &
       gxy*vy*wx + gyy*vy*wy + gyz*vy*wz + &
       gxz*vz*wx + gyz*vz*wy + gzz*vz*wz
  fx = fx*f       
  wx = wx - fx*vx
  wy = wy - fx*vy
  wz = wz - fx*vz
  fx = gxx*ux*wx + gxy*ux*wy + gxz*ux*wz + &
       gxy*uy*wx + gyy*uy*wy + gyz*uy*wz + &
       gxz*uz*wx + gyz*uz*wy + gzz*uz*wz
  fx = fx*f       
  wx = wx - fx*ux
  wy = wy - fx*uy
  wz = wz - fx*uz
  fx = gxx*wx*wx + gyy*wy*wy + gzz*wz*wz &
     +(gxy*wx*wy + gxz*wx*wz + gyz*wy*wz)*TWO
  fx = dsqrt(fx*f)
  wx = wx/fx
  wy = wy/fx
  wz = wz/fx
#endif

! E_i
  HEx = (gxx*Ex+gxy*Ey+gxz*Ez)*chipn1
  HEy = (gxy*Ex+gyy*Ey+gyz*Ez)*chipn1
  HEz = (gxz*Ex+gyz*Ey+gzz*Ez)*chipn1

  f = dsqrt(f)**3
! \sqrt(gamma)r x B
  HBx = (vy*Bz - vz*By)*f
  HBy = (vz*Bx - vx*Bz)*f
  HBz = (vx*By - vy*Bx)*f

#if (tetradtype == 1)  
!set m = (theta + i phi)/sqrt(2) following Sperhake, Eq.(3.2) of  PRD 85, 124062(2012)
!    m = (w     + i u  )/sqrt(2)
!the real part of Phi2
  Rphi2 =  (HEx-HBx)*wx+(HEy-HBy)*wy+(HEz-HBz)*wz
!the imaginary part of Phi2
  Iphi2 = -(HEx-HBx)*ux-(HEy-HBy)*uy-(HEz-HBz)*uz
         
#else  
!set m = (phi - i theta)/sqrt(2) following Frans,Eq.(8) of  PRD 75, 124018(2007)
!    m = (u   - i w    )/sqrt(2)

!the real part of Phi2
  Rphi2 =  (HEx-HBx)*ux+(HEy-HBy)*uy+(HEz-HBz)*uz
!the imaginary part of Phi2
  Iphi2 =  (HEx-HBx)*wx+(HEy-HBy)*wy+(HEz-HBz)*wz
#endif  

  Rphi2 = Rphi2/2.d0
  Iphi2 = Iphi2/2.d0

  return

  end subroutine getnpem2
!-----------------------------------------------------------------------------
!
! compute the Newman-Penrose Weyl scalar Psi4
! for BSSN dynamical variables
! for shell
!
!-----------------------------------------------------------------------------

  subroutine getnpem2_ss(ext,crho,sigma,R, X, Y, Z,                            &
               drhodx, drhody, drhodz,                                         &
               dsigmadx,dsigmady,dsigmadz,                                     &
               dRdx,dRdy,dRdz,                                                 &
               drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                &
               dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,    &
               dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,                            & 
               chi,dxx,gxy,gxz,dyy,gyz,dzz, &
               Ex,Ey,Ez,Bx,By,Bz,           &
               Rphi2, Iphi2, &
               symmetry,sst)

  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ext(1:3),symmetry,sst
  double precision,intent(in),dimension(ext(1))::crho
  double precision,intent(in),dimension(ext(2))::sigma
  double precision,intent(in),dimension(ext(3))::R
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in ) :: X,Y,Z
  double precision,intent(in),dimension(ext(1),ext(2),ext(3))::drhodx, drhody, drhodz
  double precision,intent(in),dimension(ext(1),ext(2),ext(3))::dsigmadx,dsigmady,dsigmadz
  double precision,intent(in),dimension(ext(1),ext(2),ext(3))::dRdx,dRdy,dRdz
  double precision,intent(in),dimension(ext(1),ext(2),ext(3))::drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz
  double precision,intent(in),dimension(ext(1),ext(2),ext(3))::dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz
  double precision,intent(in),dimension(ext(1),ext(2),ext(3))::dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in ) :: dxx,gxy,gxz,dyy,gyz,dzz
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in ) :: chi,Ex,Ey,Ez,Bx,By,Bz
  real*8, dimension(ext(1),ext(2),ext(3)), intent(out):: Rphi2,Iphi2

!~~~~~~> Other variables:

  real*8, dimension(ext(1),ext(2),ext(3)) :: f,fx,fy,fz
  real*8, dimension(ext(1),ext(2),ext(3)) :: gxx,gyy,gzz
  real*8, dimension(ext(1),ext(2),ext(3)) :: chipn1,chi3o2
  real*8, dimension(ext(1),ext(2),ext(3)) :: vx,vy,vz,ux,uy,uz,wx,wy,wz
  real*8, dimension(ext(1),ext(2),ext(3)) :: HEx,HEy,HEz,HBx,HBy,HBz
  real*8, dimension(ext(1),ext(2),ext(3)) :: gupxx,gupxy,gupxz
  real*8, dimension(ext(1),ext(2),ext(3)) :: gupyy,gupyz,gupzz

  real*8, parameter :: ZEO = 0.d0, ONE = 1.d0, TWO = 2.d0
  real*8, parameter :: F1o3 = 1.d0/3.d0, FOUR = 4.d0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  integer::i,j,k
  real*8,parameter::TINYRR=1.d-14

  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE
  chipn1= chi + ONE
  chi3o2  = dsqrt(chipn1)**3

! invert tilted metric
  gupzz =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz - &
           gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz
  gupxx =   ( gyy * gzz - gyz * gyz ) / gupzz
  gupxy = - ( gxy * gzz - gyz * gxz ) / gupzz
  gupxz =   ( gxy * gyz - gyy * gxz ) / gupzz
  gupyy =   ( gxx * gzz - gxz * gxz ) / gupzz
  gupyz = - ( gxx * gyz - gxy * gxz ) / gupzz
  gupzz =   ( gxx * gyy - gxy * gxy ) / gupzz

! initialize U, V, W vetors
! v:r; u: phi; w: theta
#if (tetradtype == 0)
  do i=1,ext(1)
  do j=1,ext(2)
  do k=1,ext(3)
     if(abs(X(i,j,k)) < TINYRR .and. abs(Y(i,j,k)) < TINYRR .and. abs(Z(i,j,k)) < TINYRR)then
        vx(i,j,k) = TINYRR
        vy(i,j,k) = TINYRR
        vz(i,j,k) = TINYRR
     else
        vx(i,j,k) = X(i,j,k)
        vy(i,j,k) = Y(i,j,k)
        vz(i,j,k) = Z(i,j,k)
     endif
     if(abs(X(i,j,k)) < TINYRR .and. abs(Y(i,j,k)) < TINYRR)then
        ux(i,j,k) = - TINYRR
        uy(i,j,k) = TINYRR
        uz(i,j,k) = ZEO
        wx(i,j,k) = TINYRR*Z(i,j,k)
        wy(i,j,k) = TINYRR*Z(i,j,k)
        wz(i,j,k) = -2*TINYRR*TINYRR
     else
        ux(i,j,k) = - Y(i,j,k)
        uy(i,j,k) = X(i,j,k)
        uz(i,j,k) = ZEO
        wx(i,j,k) = X(i,j,k)*Z(i,j,k)
        wy(i,j,k) = Y(i,j,k)*Z(i,j,k)
        wz(i,j,k) = -(X(i,j,k)*X(i,j,k) + Y(i,j,k)*Y(i,j,k))
     endif
  enddo
  enddo
  enddo

  f = 1.d0/chipn1

  fx = gxx*vx*vx + gyy*vy*vy + gzz*vz*vz &
     +(gxy*vx*vy + gxz*vx*vz + gyz*vy*vz)*TWO
  fx = dsqrt(fx*f)
  vx = vx/fx
  vy = vy/fx
  vz = vz/fx

  fx = gxx*vx*ux + gxy*vx*uy + gxz*vx*uz + &
       gxy*vy*ux + gyy*vy*uy + gyz*vy*uz + &
       gxz*vz*ux + gyz*vz*uy + gzz*vz*uz
  fx = fx*f
  ux = ux - fx*vx
  uy = uy - fx*vy
  uz = uz - fx*vz
  fx = gxx*ux*ux + gyy*uy*uy + gzz*uz*uz &
     +(gxy*ux*uy + gxz*ux*uz + gyz*uy*uz)*TWO
  fx = dsqrt(fx*f) 
  ux = ux/fx
  uy = uy/fx
  uz = uz/fx

  fx = gxx*vx*wx + gxy*vx*wy + gxz*vx*wz + &
       gxy*vy*wx + gyy*vy*wy + gyz*vy*wz + &
       gxz*vz*wx + gyz*vz*wy + gzz*vz*wz
  fx = fx*f       
  wx = wx - fx*vx
  wy = wy - fx*vy
  wz = wz - fx*vz
  fx = gxx*ux*wx + gxy*ux*wy + gxz*ux*wz + &
       gxy*uy*wx + gyy*uy*wy + gyz*uy*wz + &
       gxz*uz*wx + gyz*uz*wy + gzz*uz*wz
  fx = fx*f       
  wx = wx - fx*ux
  wy = wy - fx*uy
  wz = wz - fx*uz
  fx = gxx*wx*wx + gyy*wy*wy + gzz*wz*wz &
     +(gxy*wx*wy + gxz*wx*wz + gyz*wy*wz)*TWO
  fx = dsqrt(fx*f)
  wx = wx/fx
  wy = wy/fx
  wz = wz/fx
#elif  (tetradtype == 1)
  do i=1,ext(1)
  do j=1,ext(2)
  do k=1,ext(3)
     if(abs(X(i,j,k)) < TINYRR .and. abs(Y(i,j,k)) < TINYRR .and. abs(Z(i,j,k)) < TINYRR)then
        vx(i,j,k) = TINYRR
        vy(i,j,k) = TINYRR
        vz(i,j,k) = TINYRR
     else
        vx(i,j,k) = X(i,j,k)
        vy(i,j,k) = Y(i,j,k)
        vz(i,j,k) = Z(i,j,k)
     endif
     if(abs(X(i,j,k)) < TINYRR .and. abs(Y(i,j,k)) < TINYRR)then
        ux(i,j,k) = - TINYRR
        uy(i,j,k) = TINYRR
        uz(i,j,k) = ZEO
        wx(i,j,k) = TINYRR*Z(i,j,k)
        wy(i,j,k) = TINYRR*Z(i,j,k)
        wz(i,j,k) = -2*TINYRR*TINYRR
     else
        ux(i,j,k) = - Y(i,j,k)
        uy(i,j,k) = X(i,j,k)
        uz(i,j,k) = ZEO
        wx(i,j,k) = X(i,j,k)*Z(i,j,k)
        wy(i,j,k) = Y(i,j,k)*Z(i,j,k)
        wz(i,j,k) = -(X(i,j,k)*X(i,j,k) + Y(i,j,k)*Y(i,j,k))
     endif
  enddo
  enddo
  enddo

  f = 1.d0/chipn1

  fx = gxx*wx*wx + gyy*wy*wy + gzz*wz*wz &
     +(gxy*wx*wy + gxz*wx*wz + gyz*wy*wz)*TWO
  fx = dsqrt(fx*f)
  wx = wx/fx
  wy = wy/fx
  wz = wz/fx

  fx = gxx*wx*ux + gxy*wx*uy + gxz*wx*uz + &
       gxy*wy*ux + gyy*wy*uy + gyz*wy*uz + &
       gxz*wz*ux + gyz*wz*uy + gzz*wz*uz
  fx = fx*f
  ux = ux - fx*wx
  uy = uy - fx*wy
  uz = uz - fx*wz
  fx = gxx*ux*ux + gyy*uy*uy + gzz*uz*uz &
     +(gxy*ux*uy + gxz*ux*uz + gyz*uy*uz)*TWO
  fx = dsqrt(fx*f) 
  ux = ux/fx
  uy = uy/fx
  uz = uz/fx

  fx = gxx*vx*wx + gxy*vx*wy + gxz*vx*wz + &
       gxy*vy*wx + gyy*vy*wy + gyz*vy*wz + &
       gxz*vz*wx + gyz*vz*wy + gzz*vz*wz
  fx = fx*f       
  vx = vx - fx*wx
  vy = vy - fx*wy
  vz = vz - fx*wz
  fx = gxx*ux*vx + gxy*ux*vy + gxz*ux*vz + &
       gxy*uy*vx + gyy*uy*vy + gyz*uy*vz + &
       gxz*uz*vx + gyz*uz*vy + gzz*uz*vz
  fx = fx*f       
  vx = vx - fx*ux
  vy = vy - fx*uy
  vz = vz - fx*uz
  fx = gxx*vx*vx + gyy*vy*vy + gzz*vz*vz &
     +(gxy*vx*vy + gxz*vx*vz + gyz*vy*vz)*TWO
  fx = dsqrt(fx*f)
  vx = vx/fx
  vy = vy/fx
  vz = vz/fx
#elif (tetradtype == 2)  
  do i=1,ext(1)
  do j=1,ext(2)
  do k=1,ext(3)
     if(abs(X(i,j,k)) < TINYRR .and. abs(Y(i,j,k)) < TINYRR .and. abs(Z(i,j,k)) < TINYRR)then
        vx(i,j,k) = TINYRR
        vy(i,j,k) = TINYRR
        vz(i,j,k) = TINYRR
     else
        vx(i,j,k) = X(i,j,k)
        vy(i,j,k) = Y(i,j,k)
        vz(i,j,k) = Z(i,j,k)
     endif
     if(abs(X(i,j,k)) < TINYRR .and. abs(Y(i,j,k)) < TINYRR)then
        ux(i,j,k) = - TINYRR
        uy(i,j,k) = TINYRR
        uz(i,j,k) = ZEO
        wx(i,j,k) = TINYRR*Z(i,j,k)
        wy(i,j,k) = TINYRR*Z(i,j,k)
        wz(i,j,k) = -2*TINYRR*TINYRR
     else
        ux(i,j,k) = - Y(i,j,k)
        uy(i,j,k) = X(i,j,k)
        uz(i,j,k) = ZEO
        wx(i,j,k) = X(i,j,k)*Z(i,j,k)
        wy(i,j,k) = Y(i,j,k)*Z(i,j,k)
        wz(i,j,k) = -(X(i,j,k)*X(i,j,k) + Y(i,j,k)*Y(i,j,k))
     endif
  enddo
  enddo
  enddo

  fx = vx
  fy = vy
  fz = vz
  vx = gupxx*fx + gupxy*fy + gupxz*fz
  vy = gupxy*fx + gupyy*fy + gupyz*fz
  vz = gupxz*fx + gupyz*fy + gupzz*fz

  f = 1.d0/chipn1

  fx = gxx*vx*vx + gyy*vy*vy + gzz*vz*vz &
     +(gxy*vx*vy + gxz*vx*vz + gyz*vy*vz)*TWO
  fx = dsqrt(fx*f)
  vx = vx/fx
  vy = vy/fx
  vz = vz/fx

  fx = gxx*vx*ux + gxy*vx*uy + gxz*vx*uz + &
       gxy*vy*ux + gyy*vy*uy + gyz*vy*uz + &
       gxz*vz*ux + gyz*vz*uy + gzz*vz*uz
  fx = fx*f
  ux = ux - fx*vx
  uy = uy - fx*vy
  uz = uz - fx*vz
  fx = gxx*ux*ux + gyy*uy*uy + gzz*uz*uz &
     +(gxy*ux*uy + gxz*ux*uz + gyz*uy*uz)*TWO
  fx = dsqrt(fx*f) 
  ux = ux/fx
  uy = uy/fx
  uz = uz/fx

  fx = gxx*vx*wx + gxy*vx*wy + gxz*vx*wz + &
       gxy*vy*wx + gyy*vy*wy + gyz*vy*wz + &
       gxz*vz*wx + gyz*vz*wy + gzz*vz*wz
  fx = fx*f       
  wx = wx - fx*vx
  wy = wy - fx*vy
  wz = wz - fx*vz
  fx = gxx*ux*wx + gxy*ux*wy + gxz*ux*wz + &
       gxy*uy*wx + gyy*uy*wy + gyz*uy*wz + &
       gxz*uz*wx + gyz*uz*wy + gzz*uz*wz
  fx = fx*f       
  wx = wx - fx*ux
  wy = wy - fx*uy
  wz = wz - fx*uz
  fx = gxx*wx*wx + gyy*wy*wy + gzz*wz*wz &
     +(gxy*wx*wy + gxz*wx*wz + gyz*wy*wz)*TWO
  fx = dsqrt(fx*f)
  wx = wx/fx
  wy = wy/fx
  wz = wz/fx
#endif

! E_i
  HEx = (gxx*Ex+gxy*Ey+gxz*Ez)*chipn1
  HEy = (gxy*Ex+gyy*Ey+gyz*Ez)*chipn1
  HEz = (gxz*Ex+gyz*Ey+gzz*Ez)*chipn1

  f = dsqrt(f)**3
!set m = (u + iw)/sqrt(2) following Frans, PRD 75, 124018(2007)

! \sqrt(gamma)r x B
  HBx = (vy*Bz - vz*By)*f
  HBy = (vz*Bx - vx*Bz)*f
  HBz = (vx*By - vy*Bx)*f

#if (tetradtype == 1)  
!set m = (theta + i phi)/sqrt(2) following Sperhake, Eq.(3.2) of  PRD 85, 124062(2012)
!    m = (w     + i u  )/sqrt(2)
!the real part of Phi2
  Rphi2 =  (HEx-HBx)*wx+(HEy-HBy)*wy+(HEz-HBz)*wz
!the imaginary part of Phi2
  Iphi2 = -(HEx-HBx)*ux-(HEy-HBy)*uy-(HEz-HBz)*uz
         
#else  
!set m = (phi - i theta)/sqrt(2) following Frans,Eq.(8) of  PRD 75, 124018(2007)
!    m = (u   - i w    )/sqrt(2)

!the real part of Phi2
  Rphi2 =  (HEx-HBx)*ux+(HEy-HBy)*uy+(HEz-HBz)*uz
!the imaginary part of Phi2
  Iphi2 =  (HEx-HBx)*wx+(HEy-HBy)*wy+(HEz-HBz)*wz
#endif  
                  
  Rphi2 = Rphi2/2.d0
  Iphi2 = Iphi2/2.d0

  return

  end subroutine getnpem2_ss
!-----------------------------------------------------------------------------
!
! compute the EM wave phi2
! for BSSN dynamical variables
! for single point
!-----------------------------------------------------------------------------

  subroutine getnpem2_point(X, Y, Z,        &
               chi,dxx,gxy,gxz,dyy,gyz,dzz, &
               Ex,Ey,Ez,Bx,By,Bz,           &
               Rphi2, Iphi2)

  implicit none

!~~~~~~> Input parameters:
  real*8, intent(in ) :: X,Y,Z
  real*8, intent(in ) :: dxx,gxy,gxz,dyy,gyz,dzz
  real*8, intent(in ) :: chi,Ex,Ey,Ez,Bx,By,Bz
  real*8, intent(out):: Rphi2,Iphi2

!~~~~~~> Other variables:

  real*8 :: f,fx,fy,fz
  real*8 :: gxx,gyy,gzz
  real*8 :: chipn1,chi3o2
  real*8 :: vx,vy,vz,ux,uy,uz,wx,wy,wz
  real*8 :: HEx,HEy,HEz,HBx,HBy,HBz
  real*8 :: gupxx,gupxy,gupxz
  real*8 :: gupyy,gupyz,gupzz

  real*8, parameter :: ZEO = 0.d0, ONE = 1.d0, TWO = 2.d0
  real*8, parameter :: F1o3 = 1.d0/3.d0, FOUR = 4.d0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  real*8,parameter::TINYRR=1.d-14

  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE
  chipn1= chi + ONE
  chi3o2  = dsqrt(chipn1)**3

! invert tilted metric
  gupzz =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz - &
           gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz
  gupxx =   ( gyy * gzz - gyz * gyz ) / gupzz
  gupxy = - ( gxy * gzz - gyz * gxz ) / gupzz
  gupxz =   ( gxy * gyz - gyy * gxz ) / gupzz
  gupyy =   ( gxx * gzz - gxz * gxz ) / gupzz
  gupyz = - ( gxx * gyz - gxy * gxz ) / gupzz
  gupzz =   ( gxx * gyy - gxy * gxy ) / gupzz

! initialize U, V, W vetors
! v:r; u: phi; w: theta
#if (tetradtype == 0)
     if(abs(X) < TINYRR .and. abs(Y) < TINYRR .and. abs(Z) < TINYRR)then
        vx = TINYRR
        vy = TINYRR
        vz = TINYRR
     else
        vx = X
        vy = Y
        vz = Z
     endif
     if(abs(X) < TINYRR .and. abs(Y) < TINYRR)then
        ux = - TINYRR
        uy = TINYRR
        uz = ZEO
        wx = TINYRR*Z
        wy = TINYRR*Z
        wz = -2*TINYRR*TINYRR
     else
        ux = - Y
        uy = X
        uz = ZEO
        wx = X*Z
        wy = Y*Z
        wz = -(X*X + Y*Y)
     endif

  f = 1.d0/chipn1

  fx = gxx*vx*vx + gyy*vy*vy + gzz*vz*vz &
     +(gxy*vx*vy + gxz*vx*vz + gyz*vy*vz)*TWO
  fx = dsqrt(fx*f)
  vx = vx/fx
  vy = vy/fx
  vz = vz/fx

  fx = gxx*vx*ux + gxy*vx*uy + gxz*vx*uz + &
       gxy*vy*ux + gyy*vy*uy + gyz*vy*uz + &
       gxz*vz*ux + gyz*vz*uy + gzz*vz*uz
  fx = fx*f
  ux = ux - fx*vx
  uy = uy - fx*vy
  uz = uz - fx*vz
  fx = gxx*ux*ux + gyy*uy*uy + gzz*uz*uz &
     +(gxy*ux*uy + gxz*ux*uz + gyz*uy*uz)*TWO
  fx = dsqrt(fx*f) 
  ux = ux/fx
  uy = uy/fx
  uz = uz/fx

  fx = gxx*vx*wx + gxy*vx*wy + gxz*vx*wz + &
       gxy*vy*wx + gyy*vy*wy + gyz*vy*wz + &
       gxz*vz*wx + gyz*vz*wy + gzz*vz*wz
  fx = fx*f       
  wx = wx - fx*vx
  wy = wy - fx*vy
  wz = wz - fx*vz
  fx = gxx*ux*wx + gxy*ux*wy + gxz*ux*wz + &
       gxy*uy*wx + gyy*uy*wy + gyz*uy*wz + &
       gxz*uz*wx + gyz*uz*wy + gzz*uz*wz
  fx = fx*f       
  wx = wx - fx*ux
  wy = wy - fx*uy
  wz = wz - fx*uz
  fx = gxx*wx*wx + gyy*wy*wy + gzz*wz*wz &
     +(gxy*wx*wy + gxz*wx*wz + gyz*wy*wz)*TWO
  fx = dsqrt(fx*f)
  wx = wx/fx
  wy = wy/fx
  wz = wz/fx
#elif  (tetradtype == 1)
     if(abs(X) < TINYRR .and. abs(Y) < TINYRR .and. abs(Z) < TINYRR)then
        vx = TINYRR
        vy = TINYRR
        vz = TINYRR
     else
        vx = X
        vy = Y
        vz = Z
     endif
     if(abs(X) < TINYRR .and. abs(Y) < TINYRR)then
        ux = - TINYRR
        uy = TINYRR
        uz = ZEO
        wx = TINYRR*Z
        wy = TINYRR*Z
        wz = -2*TINYRR*TINYRR
     else
        ux = - Y
        uy = X
        uz = ZEO
        wx = X*Z
        wy = Y*Z
        wz = -(X*X + Y*Y)
     endif

  f = 1.d0/chipn1

  fx = gxx*wx*wx + gyy*wy*wy + gzz*wz*wz &
     +(gxy*wx*wy + gxz*wx*wz + gyz*wy*wz)*TWO
  fx = dsqrt(fx*f)
  wx = wx/fx
  wy = wy/fx
  wz = wz/fx

  fx = gxx*wx*ux + gxy*wx*uy + gxz*wx*uz + &
       gxy*wy*ux + gyy*wy*uy + gyz*wy*uz + &
       gxz*wz*ux + gyz*wz*uy + gzz*wz*uz
  fx = fx*f
  ux = ux - fx*wx
  uy = uy - fx*wy
  uz = uz - fx*wz
  fx = gxx*ux*ux + gyy*uy*uy + gzz*uz*uz &
     +(gxy*ux*uy + gxz*ux*uz + gyz*uy*uz)*TWO
  fx = dsqrt(fx*f) 
  ux = ux/fx
  uy = uy/fx
  uz = uz/fx

  fx = gxx*vx*wx + gxy*vx*wy + gxz*vx*wz + &
       gxy*vy*wx + gyy*vy*wy + gyz*vy*wz + &
       gxz*vz*wx + gyz*vz*wy + gzz*vz*wz
  fx = fx*f       
  vx = vx - fx*wx
  vy = vy - fx*wy
  vz = vz - fx*wz
  fx = gxx*ux*vx + gxy*ux*vy + gxz*ux*vz + &
       gxy*uy*vx + gyy*uy*vy + gyz*uy*vz + &
       gxz*uz*vx + gyz*uz*vy + gzz*uz*vz
  fx = fx*f       
  vx = vx - fx*ux
  vy = vy - fx*uy
  vz = vz - fx*uz
  fx = gxx*vx*vx + gyy*vy*vy + gzz*vz*vz &
     +(gxy*vx*vy + gxz*vx*vz + gyz*vy*vz)*TWO
  fx = dsqrt(fx*f)
  vx = vx/fx
  vy = vy/fx
  vz = vz/fx
#elif (tetradtype == 2)  
     if(abs(X) < TINYRR .and. abs(Y) < TINYRR .and. abs(Z) < TINYRR)then
        vx = TINYRR
        vy = TINYRR
        vz = TINYRR
     else
        vx = X
        vy = Y
        vz = Z
     endif
     if(abs(X) < TINYRR .and. abs(Y) < TINYRR)then
        ux = - TINYRR
        uy = TINYRR
        uz = ZEO
        wx = TINYRR*Z
        wy = TINYRR*Z
        wz = -2*TINYRR*TINYRR
     else
        ux = - Y
        uy = X
        uz = ZEO
        wx = X*Z
        wy = Y*Z
        wz = -(X*X + Y*Y)
     endif

  fx = vx
  fy = vy
  fz = vz
  vx = gupxx*fx + gupxy*fy + gupxz*fz
  vy = gupxy*fx + gupyy*fy + gupyz*fz
  vz = gupxz*fx + gupyz*fy + gupzz*fz

  f = 1.d0/chipn1

  fx = gxx*vx*vx + gyy*vy*vy + gzz*vz*vz &
     +(gxy*vx*vy + gxz*vx*vz + gyz*vy*vz)*TWO
  fx = dsqrt(fx*f)
  vx = vx/fx
  vy = vy/fx
  vz = vz/fx

  fx = gxx*vx*ux + gxy*vx*uy + gxz*vx*uz + &
       gxy*vy*ux + gyy*vy*uy + gyz*vy*uz + &
       gxz*vz*ux + gyz*vz*uy + gzz*vz*uz
  fx = fx*f
  ux = ux - fx*vx
  uy = uy - fx*vy
  uz = uz - fx*vz
  fx = gxx*ux*ux + gyy*uy*uy + gzz*uz*uz &
     +(gxy*ux*uy + gxz*ux*uz + gyz*uy*uz)*TWO
  fx = dsqrt(fx*f) 
  ux = ux/fx
  uy = uy/fx
  uz = uz/fx

  fx = gxx*vx*wx + gxy*vx*wy + gxz*vx*wz + &
       gxy*vy*wx + gyy*vy*wy + gyz*vy*wz + &
       gxz*vz*wx + gyz*vz*wy + gzz*vz*wz
  fx = fx*f       
  wx = wx - fx*vx
  wy = wy - fx*vy
  wz = wz - fx*vz
  fx = gxx*ux*wx + gxy*ux*wy + gxz*ux*wz + &
       gxy*uy*wx + gyy*uy*wy + gyz*uy*wz + &
       gxz*uz*wx + gyz*uz*wy + gzz*uz*wz
  fx = fx*f       
  wx = wx - fx*ux
  wy = wy - fx*uy
  wz = wz - fx*uz
  fx = gxx*wx*wx + gyy*wy*wy + gzz*wz*wz &
     +(gxy*wx*wy + gxz*wx*wz + gyz*wy*wz)*TWO
  fx = dsqrt(fx*f)
  wx = wx/fx
  wy = wy/fx
  wz = wz/fx
#endif

! E_i
  HEx = (gxx*Ex+gxy*Ey+gxz*Ez)*chipn1
  HEy = (gxy*Ex+gyy*Ey+gyz*Ez)*chipn1
  HEz = (gxz*Ex+gyz*Ey+gzz*Ez)*chipn1

  f = dsqrt(f)**3
!set m = (u + iw)/sqrt(2) following Frans, PRD 75, 124018(2007)

! \sqrt(gamma)r x B
  HBx = (vy*Bz - vz*By)*f
  HBy = (vz*Bx - vx*Bz)*f
  HBz = (vx*By - vy*Bx)*f

#if (tetradtype == 1)  
!set m = (theta + i phi)/sqrt(2) following Sperhake, Eq.(3.2) of  PRD 85, 124062(2012)
!    m = (w     + i u  )/sqrt(2)
!the real part of Phi2
  Rphi2 =  (HEx-HBx)*wx+(HEy-HBy)*wy+(HEz-HBz)*wz
!the imaginary part of Phi2
  Iphi2 = -(HEx-HBx)*ux-(HEy-HBy)*uy-(HEz-HBz)*uz
         
#else  
!set m = (phi - i theta)/sqrt(2) following Frans,Eq.(8) of  PRD 75, 124018(2007)
!    m = (u   - i w    )/sqrt(2)

!the real part of Phi2
  Rphi2 =  (HEx-HBx)*ux+(HEy-HBy)*uy+(HEz-HBz)*uz
!the imaginary part of Phi2
  Iphi2 =  (HEx-HBx)*wx+(HEy-HBy)*wy+(HEz-HBz)*wz
#endif  
         
  Rphi2 = Rphi2/2.d0
  Iphi2 = Iphi2/2.d0

  return

  end subroutine getnpem2_point
!-----------------------------------------------------------------------------
!
! compute the Newman-Penrose Weyl scalar Psi4
! for BSSN dynamical variables
!
!-----------------------------------------------------------------------------

  subroutine getnpem1(ext, X, Y, Z,         &
               chi,dxx,gxy,gxz,dyy,gyz,dzz, &
               Ex,Ey,Ez,Bx,By,Bz,           &
               Rphi1, Iphi1, &
               symmetry)

  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ext(1:3),symmetry
  real*8, intent(in ):: X(1:ext(1)),Y(1:ext(2)),Z(1:ext(3))
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in ) :: dxx,gxy,gxz,dyy,gyz,dzz
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in ) :: chi,Ex,Ey,Ez,Bx,By,Bz
  real*8, dimension(ext(1),ext(2),ext(3)), intent(out):: Rphi1,Iphi1

!~~~~~~> Other variables:

  real*8, dimension(ext(1),ext(2),ext(3)) :: f,fx,fy,fz
  real*8, dimension(ext(1),ext(2),ext(3)) :: gxx,gyy,gzz
  real*8, dimension(ext(1),ext(2),ext(3)) :: chipn1,chi3o2
  real*8, dimension(ext(1),ext(2),ext(3)) :: vx,vy,vz,ux,uy,uz,wx,wy,wz
  real*8, dimension(ext(1),ext(2),ext(3)) :: HEx,HEy,HEz,HBx,HBy,HBz
  real*8, dimension(ext(1),ext(2),ext(3)) :: gupxx,gupxy,gupxz
  real*8, dimension(ext(1),ext(2),ext(3)) :: gupyy,gupyz,gupzz

  real*8, parameter :: ZEO = 0.d0, ONE = 1.d0, TWO = 2.d0
  real*8, parameter :: F1o3 = 1.d0/3.d0, FOUR = 4.d0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  real*8            :: sqr2
  integer::i,j,k
  real*8,parameter::TINYRR=1.d-14

  sqr2 = dsqrt(2.d0)

  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE
  chipn1= chi + ONE
  chi3o2  = dsqrt(chipn1)**3

! invert tilted metric
  gupzz =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz - &
           gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz
  gupxx =   ( gyy * gzz - gyz * gyz ) / gupzz
  gupxy = - ( gxy * gzz - gyz * gxz ) / gupzz
  gupxz =   ( gxy * gyz - gyy * gxz ) / gupzz
  gupyy =   ( gxx * gzz - gxz * gxz ) / gupzz
  gupyz = - ( gxx * gyz - gxy * gxz ) / gupzz
  gupzz =   ( gxx * gyy - gxy * gxy ) / gupzz

! initialize U, V, W vetors	
#if (tetradtype == 0)
  do i=1,ext(1)
  do j=1,ext(2)
  do k=1,ext(3)
     if(abs(X(i)) < TINYRR .and. abs(Y(j)) < TINYRR .and. abs(Z(k)) < TINYRR)then
        vx(i,j,k) = TINYRR
        vy(i,j,k) = TINYRR
        vz(i,j,k) = TINYRR
     else
        vx(i,j,k) = X(i)
        vy(i,j,k) = Y(j)
        vz(i,j,k) = Z(k)
     endif
     if(abs(X(i)) < TINYRR .and. abs(Y(j)) < TINYRR)then
        ux(i,j,k) = - TINYRR
        uy(i,j,k) = TINYRR
        uz(i,j,k) = ZEO
        wx(i,j,k) = TINYRR*Z(k)
        wy(i,j,k) = TINYRR*Z(k)
        wz(i,j,k) = -2*TINYRR*TINYRR
     else
        ux(i,j,k) = - Y(j)
        uy(i,j,k) = X(i)
        uz(i,j,k) = ZEO
        wx(i,j,k) = X(i)*Z(k)
        wy(i,j,k) = Y(j)*Z(k)
        wz(i,j,k) = -(X(i)*X(i) + Y(j)*Y(j))
     endif
  enddo
  enddo
  enddo

  f = 1.d0/chipn1

  fx = gxx*vx*vx + gyy*vy*vy + gzz*vz*vz &
     +(gxy*vx*vy + gxz*vx*vz + gyz*vy*vz)*TWO
  fx = dsqrt(fx*f)
  vx = vx/fx
  vy = vy/fx
  vz = vz/fx

  fx = gxx*vx*ux + gxy*vx*uy + gxz*vx*uz + &
       gxy*vy*ux + gyy*vy*uy + gyz*vy*uz + &
       gxz*vz*ux + gyz*vz*uy + gzz*vz*uz
  fx = fx*f
  ux = ux - fx*vx
  uy = uy - fx*vy
  uz = uz - fx*vz
  fx = gxx*ux*ux + gyy*uy*uy + gzz*uz*uz &
     +(gxy*ux*uy + gxz*ux*uz + gyz*uy*uz)*TWO
  fx = dsqrt(fx*f) 
  ux = ux/fx
  uy = uy/fx
  uz = uz/fx

  fx = gxx*vx*wx + gxy*vx*wy + gxz*vx*wz + &
       gxy*vy*wx + gyy*vy*wy + gyz*vy*wz + &
       gxz*vz*wx + gyz*vz*wy + gzz*vz*wz
  fx = fx*f       
  wx = wx - fx*vx
  wy = wy - fx*vy
  wz = wz - fx*vz
  fx = gxx*ux*wx + gxy*ux*wy + gxz*ux*wz + &
       gxy*uy*wx + gyy*uy*wy + gyz*uy*wz + &
       gxz*uz*wx + gyz*uz*wy + gzz*uz*wz
  fx = fx*f       
  wx = wx - fx*ux
  wy = wy - fx*uy
  wz = wz - fx*uz
  fx = gxx*wx*wx + gyy*wy*wy + gzz*wz*wz &
     +(gxy*wx*wy + gxz*wx*wz + gyz*wy*wz)*TWO
  fx = dsqrt(fx*f)
  wx = wx/fx
  wy = wy/fx
  wz = wz/fx
#elif  (tetradtype == 1)
  do i=1,ext(1)
  do j=1,ext(2)
  do k=1,ext(3)
     if(abs(X(i)) < TINYRR .and. abs(Y(j)) < TINYRR .and. abs(Z(k)) < TINYRR)then
        vx(i,j,k) = TINYRR
        vy(i,j,k) = TINYRR
        vz(i,j,k) = TINYRR
     else
        vx(i,j,k) = X(i)
        vy(i,j,k) = Y(j)
        vz(i,j,k) = Z(k)
     endif
     if(abs(X(i)) < TINYRR .and. abs(Y(j)) < TINYRR)then
        ux(i,j,k) = - TINYRR
        uy(i,j,k) = TINYRR
        uz(i,j,k) = ZEO
        wx(i,j,k) = TINYRR*Z(k)
        wy(i,j,k) = TINYRR*Z(k)
        wz(i,j,k) = -2*TINYRR*TINYRR
     else
        ux(i,j,k) = - Y(j)
        uy(i,j,k) = X(i)
        uz(i,j,k) = ZEO
        wx(i,j,k) = X(i)*Z(k)
        wy(i,j,k) = Y(j)*Z(k)
        wz(i,j,k) = -(X(i)*X(i) + Y(j)*Y(j))
     endif
  enddo
  enddo
  enddo

  f = 1.d0/chipn1

  fx = gxx*wx*wx + gyy*wy*wy + gzz*wz*wz &
     +(gxy*wx*wy + gxz*wx*wz + gyz*wy*wz)*TWO
  fx = dsqrt(fx*f)
  wx = wx/fx
  wy = wy/fx
  wz = wz/fx

  fx = gxx*wx*ux + gxy*wx*uy + gxz*wx*uz + &
       gxy*wy*ux + gyy*wy*uy + gyz*wy*uz + &
       gxz*wz*ux + gyz*wz*uy + gzz*wz*uz
  fx = fx*f
  ux = ux - fx*wx
  uy = uy - fx*wy
  uz = uz - fx*wz
  fx = gxx*ux*ux + gyy*uy*uy + gzz*uz*uz &
     +(gxy*ux*uy + gxz*ux*uz + gyz*uy*uz)*TWO
  fx = dsqrt(fx*f) 
  ux = ux/fx
  uy = uy/fx
  uz = uz/fx

  fx = gxx*vx*wx + gxy*vx*wy + gxz*vx*wz + &
       gxy*vy*wx + gyy*vy*wy + gyz*vy*wz + &
       gxz*vz*wx + gyz*vz*wy + gzz*vz*wz
  fx = fx*f       
  vx = vx - fx*wx
  vy = vy - fx*wy
  vz = vz - fx*wz
  fx = gxx*ux*vx + gxy*ux*vy + gxz*ux*vz + &
       gxy*uy*vx + gyy*uy*vy + gyz*uy*vz + &
       gxz*uz*vx + gyz*uz*vy + gzz*uz*vz
  fx = fx*f       
  vx = vx - fx*ux
  vy = vy - fx*uy
  vz = vz - fx*uz
  fx = gxx*vx*vx + gyy*vy*vy + gzz*vz*vz &
     +(gxy*vx*vy + gxz*vx*vz + gyz*vy*vz)*TWO
  fx = dsqrt(fx*f)
  vx = vx/fx
  vy = vy/fx
  vz = vz/fx
#elif (tetradtype == 2)  
  do i=1,ext(1)
  do j=1,ext(2)
  do k=1,ext(3)
     if(abs(X(i)) < TINYRR .and. abs(Y(j)) < TINYRR .and. abs(Z(k)) < TINYRR)then
        vx(i,j,k) = TINYRR
        vy(i,j,k) = TINYRR
        vz(i,j,k) = TINYRR
     else
        vx(i,j,k) = X(i)
        vy(i,j,k) = Y(j)
        vz(i,j,k) = Z(k)
     endif
     if(abs(X(i)) < TINYRR .and. abs(Y(j)) < TINYRR)then
        ux(i,j,k) = - TINYRR
        uy(i,j,k) = TINYRR
        uz(i,j,k) = ZEO
        wx(i,j,k) = TINYRR*Z(k)
        wy(i,j,k) = TINYRR*Z(k)
        wz(i,j,k) = -2*TINYRR*TINYRR
     else
        ux(i,j,k) = - Y(j)
        uy(i,j,k) = X(i)
        uz(i,j,k) = ZEO
        wx(i,j,k) = X(i)*Z(k)
        wy(i,j,k) = Y(j)*Z(k)
        wz(i,j,k) = -(X(i)*X(i) + Y(j)*Y(j))
     endif
  enddo
  enddo
  enddo

  fx = vx
  fy = vy
  fz = vz
  vx = gupxx*fx + gupxy*fy + gupxz*fz
  vy = gupxy*fx + gupyy*fy + gupyz*fz
  vz = gupxz*fx + gupyz*fy + gupzz*fz

  f = 1.d0/chipn1

  fx = gxx*vx*vx + gyy*vy*vy + gzz*vz*vz &
     +(gxy*vx*vy + gxz*vx*vz + gyz*vy*vz)*TWO
  fx = dsqrt(fx*f)
  vx = vx/fx
  vy = vy/fx
  vz = vz/fx

  fx = gxx*vx*ux + gxy*vx*uy + gxz*vx*uz + &
       gxy*vy*ux + gyy*vy*uy + gyz*vy*uz + &
       gxz*vz*ux + gyz*vz*uy + gzz*vz*uz
  fx = fx*f
  ux = ux - fx*vx
  uy = uy - fx*vy
  uz = uz - fx*vz
  fx = gxx*ux*ux + gyy*uy*uy + gzz*uz*uz &
     +(gxy*ux*uy + gxz*ux*uz + gyz*uy*uz)*TWO
  fx = dsqrt(fx*f) 
  ux = ux/fx
  uy = uy/fx
  uz = uz/fx

  fx = gxx*vx*wx + gxy*vx*wy + gxz*vx*wz + &
       gxy*vy*wx + gyy*vy*wy + gyz*vy*wz + &
       gxz*vz*wx + gyz*vz*wy + gzz*vz*wz
  fx = fx*f       
  wx = wx - fx*vx
  wy = wy - fx*vy
  wz = wz - fx*vz
  fx = gxx*ux*wx + gxy*ux*wy + gxz*ux*wz + &
       gxy*uy*wx + gyy*uy*wy + gyz*uy*wz + &
       gxz*uz*wx + gyz*uz*wy + gzz*uz*wz
  fx = fx*f       
  wx = wx - fx*ux
  wy = wy - fx*uy
  wz = wz - fx*uz
  fx = gxx*wx*wx + gyy*wy*wy + gzz*wz*wz &
     +(gxy*wx*wy + gxz*wx*wz + gyz*wy*wz)*TWO
  fx = dsqrt(fx*f)
  wx = wx/fx
  wy = wy/fx
  wz = wz/fx
#endif

  f = dsqrt(f)**3
! E_i
  HEx = (gxx*Ex+gxy*Ey+gxz*Ez)*chipn1
  HEy = (gxy*Ex+gyy*Ey+gyz*Ez)*chipn1
  HEz = (gxz*Ex+gyz*Ey+gzz*Ez)*chipn1

  Rphi1 = HEx*vx+HEy*vy+HEz*vz

! \sqrt(gamma)u x w (theta x phi)
  HBx = (uy*wz - uz*wy)*f
  HBy = (uz*wx - ux*wz)*f
  HBz = (ux*wy - uy*wx)*f
  Iphi1 = HBx*Bx+HBy*By+HBz*Bz

  Rphi1 = Rphi1/2.d0
  Iphi1 = Iphi1/2.d0

  return

  end subroutine getnpem1
!-----------------------------------------------------------------------------
!
! compute the Newman-Penrose Weyl scalar Psi4
! for BSSN dynamical variables
! for shell
!
!-----------------------------------------------------------------------------

  subroutine getnpem1_ss(ext,crho,sigma,R, X, Y, Z,                            &
               drhodx, drhody, drhodz,                                         &
               dsigmadx,dsigmady,dsigmadz,                                     &
               dRdx,dRdy,dRdz,                                                 &
               drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                &
               dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,    &
               dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,                            & 
               chi,dxx,gxy,gxz,dyy,gyz,dzz, &
               Ex,Ey,Ez,Bx,By,Bz,           &
               Rphi1, Iphi1, &
               symmetry,sst)

  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ext(1:3),symmetry,sst
  double precision,intent(in),dimension(ext(1))::crho
  double precision,intent(in),dimension(ext(2))::sigma
  double precision,intent(in),dimension(ext(3))::R
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in ) :: X,Y,Z
  double precision,intent(in),dimension(ext(1),ext(2),ext(3))::drhodx, drhody, drhodz
  double precision,intent(in),dimension(ext(1),ext(2),ext(3))::dsigmadx,dsigmady,dsigmadz
  double precision,intent(in),dimension(ext(1),ext(2),ext(3))::dRdx,dRdy,dRdz
  double precision,intent(in),dimension(ext(1),ext(2),ext(3))::drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz
  double precision,intent(in),dimension(ext(1),ext(2),ext(3))::dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz
  double precision,intent(in),dimension(ext(1),ext(2),ext(3))::dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in ) :: dxx,gxy,gxz,dyy,gyz,dzz
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in ) :: chi,Ex,Ey,Ez,Bx,By,Bz
  real*8, dimension(ext(1),ext(2),ext(3)), intent(out):: Rphi1,Iphi1

!~~~~~~> Other variables:

  real*8, dimension(ext(1),ext(2),ext(3)) :: f,fx,fy,fz
  real*8, dimension(ext(1),ext(2),ext(3)) :: gxx,gyy,gzz
  real*8, dimension(ext(1),ext(2),ext(3)) :: chipn1,chi3o2
  real*8, dimension(ext(1),ext(2),ext(3)) :: vx,vy,vz,ux,uy,uz,wx,wy,wz
  real*8, dimension(ext(1),ext(2),ext(3)) :: HEx,HEy,HEz,HBx,HBy,HBz
  real*8, dimension(ext(1),ext(2),ext(3)) :: gupxx,gupxy,gupxz
  real*8, dimension(ext(1),ext(2),ext(3)) :: gupyy,gupyz,gupzz

  real*8, parameter :: ZEO = 0.d0, ONE = 1.d0, TWO = 2.d0
  real*8, parameter :: F1o3 = 1.d0/3.d0, FOUR = 4.d0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  real*8            :: sqr2
  integer::i,j,k
  real*8,parameter::TINYRR=1.d-14

  sqr2 = dsqrt(2.d0)

  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE
  chipn1= chi + ONE
  chi3o2  = dsqrt(chipn1)**3

! invert tilted metric
  gupzz =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz - &
           gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz
  gupxx =   ( gyy * gzz - gyz * gyz ) / gupzz
  gupxy = - ( gxy * gzz - gyz * gxz ) / gupzz
  gupxz =   ( gxy * gyz - gyy * gxz ) / gupzz
  gupyy =   ( gxx * gzz - gxz * gxz ) / gupzz
  gupyz = - ( gxx * gyz - gxy * gxz ) / gupzz
  gupzz =   ( gxx * gyy - gxy * gxy ) / gupzz

! initialize U, V, W vetors
! v:r; u: phi; w: theta
#if (tetradtype == 0)
  do i=1,ext(1)
  do j=1,ext(2)
  do k=1,ext(3)
     if(abs(X(i,j,k)) < TINYRR .and. abs(Y(i,j,k)) < TINYRR .and. abs(Z(i,j,k)) < TINYRR)then
        vx(i,j,k) = TINYRR
        vy(i,j,k) = TINYRR
        vz(i,j,k) = TINYRR
     else
        vx(i,j,k) = X(i,j,k)
        vy(i,j,k) = Y(i,j,k)
        vz(i,j,k) = Z(i,j,k)
     endif
     if(abs(X(i,j,k)) < TINYRR .and. abs(Y(i,j,k)) < TINYRR)then
        ux(i,j,k) = - TINYRR
        uy(i,j,k) = TINYRR
        uz(i,j,k) = ZEO
        wx(i,j,k) = TINYRR*Z(i,j,k)
        wy(i,j,k) = TINYRR*Z(i,j,k)
        wz(i,j,k) = -2*TINYRR*TINYRR
     else
        ux(i,j,k) = - Y(i,j,k)
        uy(i,j,k) = X(i,j,k)
        uz(i,j,k) = ZEO
        wx(i,j,k) = X(i,j,k)*Z(i,j,k)
        wy(i,j,k) = Y(i,j,k)*Z(i,j,k)
        wz(i,j,k) = -(X(i,j,k)*X(i,j,k) + Y(i,j,k)*Y(i,j,k))
     endif
  enddo
  enddo
  enddo

  f = 1.d0/chipn1

  fx = gxx*vx*vx + gyy*vy*vy + gzz*vz*vz &
     +(gxy*vx*vy + gxz*vx*vz + gyz*vy*vz)*TWO
  fx = dsqrt(fx*f)
  vx = vx/fx
  vy = vy/fx
  vz = vz/fx

  fx = gxx*vx*ux + gxy*vx*uy + gxz*vx*uz + &
       gxy*vy*ux + gyy*vy*uy + gyz*vy*uz + &
       gxz*vz*ux + gyz*vz*uy + gzz*vz*uz
  fx = fx*f
  ux = ux - fx*vx
  uy = uy - fx*vy
  uz = uz - fx*vz
  fx = gxx*ux*ux + gyy*uy*uy + gzz*uz*uz &
     +(gxy*ux*uy + gxz*ux*uz + gyz*uy*uz)*TWO
  fx = dsqrt(fx*f) 
  ux = ux/fx
  uy = uy/fx
  uz = uz/fx

  fx = gxx*vx*wx + gxy*vx*wy + gxz*vx*wz + &
       gxy*vy*wx + gyy*vy*wy + gyz*vy*wz + &
       gxz*vz*wx + gyz*vz*wy + gzz*vz*wz
  fx = fx*f       
  wx = wx - fx*vx
  wy = wy - fx*vy
  wz = wz - fx*vz
  fx = gxx*ux*wx + gxy*ux*wy + gxz*ux*wz + &
       gxy*uy*wx + gyy*uy*wy + gyz*uy*wz + &
       gxz*uz*wx + gyz*uz*wy + gzz*uz*wz
  fx = fx*f       
  wx = wx - fx*ux
  wy = wy - fx*uy
  wz = wz - fx*uz
  fx = gxx*wx*wx + gyy*wy*wy + gzz*wz*wz &
     +(gxy*wx*wy + gxz*wx*wz + gyz*wy*wz)*TWO
  fx = dsqrt(fx*f)
  wx = wx/fx
  wy = wy/fx
  wz = wz/fx
#elif  (tetradtype == 1)
  do i=1,ext(1)
  do j=1,ext(2)
  do k=1,ext(3)
     if(abs(X(i,j,k)) < TINYRR .and. abs(Y(i,j,k)) < TINYRR .and. abs(Z(i,j,k)) < TINYRR)then
        vx(i,j,k) = TINYRR
        vy(i,j,k) = TINYRR
        vz(i,j,k) = TINYRR
     else
        vx(i,j,k) = X(i,j,k)
        vy(i,j,k) = Y(i,j,k)
        vz(i,j,k) = Z(i,j,k)
     endif
     if(abs(X(i,j,k)) < TINYRR .and. abs(Y(i,j,k)) < TINYRR)then
        ux(i,j,k) = - TINYRR
        uy(i,j,k) = TINYRR
        uz(i,j,k) = ZEO
        wx(i,j,k) = TINYRR*Z(i,j,k)
        wy(i,j,k) = TINYRR*Z(i,j,k)
        wz(i,j,k) = -2*TINYRR*TINYRR
     else
        ux(i,j,k) = - Y(i,j,k)
        uy(i,j,k) = X(i,j,k)
        uz(i,j,k) = ZEO
        wx(i,j,k) = X(i,j,k)*Z(i,j,k)
        wy(i,j,k) = Y(i,j,k)*Z(i,j,k)
        wz(i,j,k) = -(X(i,j,k)*X(i,j,k) + Y(i,j,k)*Y(i,j,k))
     endif
  enddo
  enddo
  enddo

  f = 1.d0/chipn1

  fx = gxx*wx*wx + gyy*wy*wy + gzz*wz*wz &
     +(gxy*wx*wy + gxz*wx*wz + gyz*wy*wz)*TWO
  fx = dsqrt(fx*f)
  wx = wx/fx
  wy = wy/fx
  wz = wz/fx

  fx = gxx*wx*ux + gxy*wx*uy + gxz*wx*uz + &
       gxy*wy*ux + gyy*wy*uy + gyz*wy*uz + &
       gxz*wz*ux + gyz*wz*uy + gzz*wz*uz
  fx = fx*f
  ux = ux - fx*wx
  uy = uy - fx*wy
  uz = uz - fx*wz
  fx = gxx*ux*ux + gyy*uy*uy + gzz*uz*uz &
     +(gxy*ux*uy + gxz*ux*uz + gyz*uy*uz)*TWO
  fx = dsqrt(fx*f) 
  ux = ux/fx
  uy = uy/fx
  uz = uz/fx

  fx = gxx*vx*wx + gxy*vx*wy + gxz*vx*wz + &
       gxy*vy*wx + gyy*vy*wy + gyz*vy*wz + &
       gxz*vz*wx + gyz*vz*wy + gzz*vz*wz
  fx = fx*f       
  vx = vx - fx*wx
  vy = vy - fx*wy
  vz = vz - fx*wz
  fx = gxx*ux*vx + gxy*ux*vy + gxz*ux*vz + &
       gxy*uy*vx + gyy*uy*vy + gyz*uy*vz + &
       gxz*uz*vx + gyz*uz*vy + gzz*uz*vz
  fx = fx*f       
  vx = vx - fx*ux
  vy = vy - fx*uy
  vz = vz - fx*uz
  fx = gxx*vx*vx + gyy*vy*vy + gzz*vz*vz &
     +(gxy*vx*vy + gxz*vx*vz + gyz*vy*vz)*TWO
  fx = dsqrt(fx*f)
  vx = vx/fx
  vy = vy/fx
  vz = vz/fx
#elif (tetradtype == 2)  
  do i=1,ext(1)
  do j=1,ext(2)
  do k=1,ext(3)
     if(abs(X(i,j,k)) < TINYRR .and. abs(Y(i,j,k)) < TINYRR .and. abs(Z(i,j,k)) < TINYRR)then
        vx(i,j,k) = TINYRR
        vy(i,j,k) = TINYRR
        vz(i,j,k) = TINYRR
     else
        vx(i,j,k) = X(i,j,k)
        vy(i,j,k) = Y(i,j,k)
        vz(i,j,k) = Z(i,j,k)
     endif
     if(abs(X(i,j,k)) < TINYRR .and. abs(Y(i,j,k)) < TINYRR)then
        ux(i,j,k) = - TINYRR
        uy(i,j,k) = TINYRR
        uz(i,j,k) = ZEO
        wx(i,j,k) = TINYRR*Z(i,j,k)
        wy(i,j,k) = TINYRR*Z(i,j,k)
        wz(i,j,k) = -2*TINYRR*TINYRR
     else
        ux(i,j,k) = - Y(i,j,k)
        uy(i,j,k) = X(i,j,k)
        uz(i,j,k) = ZEO
        wx(i,j,k) = X(i,j,k)*Z(i,j,k)
        wy(i,j,k) = Y(i,j,k)*Z(i,j,k)
        wz(i,j,k) = -(X(i,j,k)*X(i,j,k) + Y(i,j,k)*Y(i,j,k))
     endif
  enddo
  enddo
  enddo

  fx = vx
  fy = vy
  fz = vz
  vx = gupxx*fx + gupxy*fy + gupxz*fz
  vy = gupxy*fx + gupyy*fy + gupyz*fz
  vz = gupxz*fx + gupyz*fy + gupzz*fz

  f = 1.d0/chipn1

  fx = gxx*vx*vx + gyy*vy*vy + gzz*vz*vz &
     +(gxy*vx*vy + gxz*vx*vz + gyz*vy*vz)*TWO
  fx = dsqrt(fx*f)
  vx = vx/fx
  vy = vy/fx
  vz = vz/fx

  fx = gxx*vx*ux + gxy*vx*uy + gxz*vx*uz + &
       gxy*vy*ux + gyy*vy*uy + gyz*vy*uz + &
       gxz*vz*ux + gyz*vz*uy + gzz*vz*uz
  fx = fx*f
  ux = ux - fx*vx
  uy = uy - fx*vy
  uz = uz - fx*vz
  fx = gxx*ux*ux + gyy*uy*uy + gzz*uz*uz &
     +(gxy*ux*uy + gxz*ux*uz + gyz*uy*uz)*TWO
  fx = dsqrt(fx*f) 
  ux = ux/fx
  uy = uy/fx
  uz = uz/fx

  fx = gxx*vx*wx + gxy*vx*wy + gxz*vx*wz + &
       gxy*vy*wx + gyy*vy*wy + gyz*vy*wz + &
       gxz*vz*wx + gyz*vz*wy + gzz*vz*wz
  fx = fx*f       
  wx = wx - fx*vx
  wy = wy - fx*vy
  wz = wz - fx*vz
  fx = gxx*ux*wx + gxy*ux*wy + gxz*ux*wz + &
       gxy*uy*wx + gyy*uy*wy + gyz*uy*wz + &
       gxz*uz*wx + gyz*uz*wy + gzz*uz*wz
  fx = fx*f       
  wx = wx - fx*ux
  wy = wy - fx*uy
  wz = wz - fx*uz
  fx = gxx*wx*wx + gyy*wy*wy + gzz*wz*wz &
     +(gxy*wx*wy + gxz*wx*wz + gyz*wy*wz)*TWO
  fx = dsqrt(fx*f)
  wx = wx/fx
  wy = wy/fx
  wz = wz/fx
#endif

  f = dsqrt(f)**3
! E_i
  HEx = (gxx*Ex+gxy*Ey+gxz*Ez)*chipn1
  HEy = (gxy*Ex+gyy*Ey+gyz*Ez)*chipn1
  HEz = (gxz*Ex+gyz*Ey+gzz*Ez)*chipn1

  Rphi1 = HEx*vx+HEy*vy+HEz*vz
! \sqrt(gamma)u x w (theta x phi)
  HBx = (uy*wz - uz*wy)*f
  HBy = (uz*wx - ux*wz)*f
  HBz = (ux*wy - uy*wx)*f
  Iphi1 = HBx*Bx+HBy*By+HBz*Bz

  Rphi1 = Rphi1/2.d0
  Iphi1 = Iphi1/2.d0

  return

  end subroutine getnpem1_ss
!-----------------------------------------------------------------------------
!
! compute the EM wave phi1
! for BSSN dynamical variables
! for single point
!-----------------------------------------------------------------------------

  subroutine getnpem1_point(X, Y, Z,        &
               chi,dxx,gxy,gxz,dyy,gyz,dzz, &
               Ex,Ey,Ez,Bx,By,Bz,           &
               Rphi1, Iphi1)

  implicit none

!~~~~~~> Input parameters:
  real*8, intent(in ) :: X,Y,Z
  real*8, intent(in ) :: dxx,gxy,gxz,dyy,gyz,dzz
  real*8, intent(in ) :: chi,Ex,Ey,Ez,Bx,By,Bz
  real*8, intent(out):: Rphi1,Iphi1

!~~~~~~> Other variables:

  real*8 :: f,fx,fy,fz
  real*8 :: gxx,gyy,gzz
  real*8 :: chipn1,chi3o2
  real*8 :: vx,vy,vz,ux,uy,uz,wx,wy,wz
  real*8 :: HEx,HEy,HEz,HBx,HBy,HBz
  real*8 :: gupxx,gupxy,gupxz
  real*8 :: gupyy,gupyz,gupzz

  real*8, parameter :: ZEO = 0.d0, ONE = 1.d0, TWO = 2.d0
  real*8, parameter :: F1o3 = 1.d0/3.d0, FOUR = 4.d0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  real*8,parameter::TINYRR=1.d-14

  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE
  chipn1= chi + ONE
  chi3o2  = dsqrt(chipn1)**3

! invert tilted metric
  gupzz =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz - &
           gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz
  gupxx =   ( gyy * gzz - gyz * gyz ) / gupzz
  gupxy = - ( gxy * gzz - gyz * gxz ) / gupzz
  gupxz =   ( gxy * gyz - gyy * gxz ) / gupzz
  gupyy =   ( gxx * gzz - gxz * gxz ) / gupzz
  gupyz = - ( gxx * gyz - gxy * gxz ) / gupzz
  gupzz =   ( gxx * gyy - gxy * gxy ) / gupzz

! initialize U, V, W vetors
! v:r; u: phi; w: theta
#if (tetradtype == 0)
     if(abs(X) < TINYRR .and. abs(Y) < TINYRR .and. abs(Z) < TINYRR)then
        vx = TINYRR
        vy = TINYRR
        vz = TINYRR
     else
        vx = X
        vy = Y
        vz = Z
     endif
     if(abs(X) < TINYRR .and. abs(Y) < TINYRR)then
        ux = - TINYRR
        uy = TINYRR
        uz = ZEO
        wx = TINYRR*Z
        wy = TINYRR*Z
        wz = -2*TINYRR*TINYRR
     else
        ux = - Y
        uy = X
        uz = ZEO
        wx = X*Z
        wy = Y*Z
        wz = -(X*X + Y*Y)
     endif

  f = 1.d0/chipn1

  fx = gxx*vx*vx + gyy*vy*vy + gzz*vz*vz &
     +(gxy*vx*vy + gxz*vx*vz + gyz*vy*vz)*TWO
  fx = dsqrt(fx*f)
  vx = vx/fx
  vy = vy/fx
  vz = vz/fx

  fx = gxx*vx*ux + gxy*vx*uy + gxz*vx*uz + &
       gxy*vy*ux + gyy*vy*uy + gyz*vy*uz + &
       gxz*vz*ux + gyz*vz*uy + gzz*vz*uz
  fx = fx*f
  ux = ux - fx*vx
  uy = uy - fx*vy
  uz = uz - fx*vz
  fx = gxx*ux*ux + gyy*uy*uy + gzz*uz*uz &
     +(gxy*ux*uy + gxz*ux*uz + gyz*uy*uz)*TWO
  fx = dsqrt(fx*f) 
  ux = ux/fx
  uy = uy/fx
  uz = uz/fx

  fx = gxx*vx*wx + gxy*vx*wy + gxz*vx*wz + &
       gxy*vy*wx + gyy*vy*wy + gyz*vy*wz + &
       gxz*vz*wx + gyz*vz*wy + gzz*vz*wz
  fx = fx*f       
  wx = wx - fx*vx
  wy = wy - fx*vy
  wz = wz - fx*vz
  fx = gxx*ux*wx + gxy*ux*wy + gxz*ux*wz + &
       gxy*uy*wx + gyy*uy*wy + gyz*uy*wz + &
       gxz*uz*wx + gyz*uz*wy + gzz*uz*wz
  fx = fx*f       
  wx = wx - fx*ux
  wy = wy - fx*uy
  wz = wz - fx*uz
  fx = gxx*wx*wx + gyy*wy*wy + gzz*wz*wz &
     +(gxy*wx*wy + gxz*wx*wz + gyz*wy*wz)*TWO
  fx = dsqrt(fx*f)
  wx = wx/fx
  wy = wy/fx
  wz = wz/fx
#elif  (tetradtype == 1)
     if(abs(X) < TINYRR .and. abs(Y) < TINYRR .and. abs(Z) < TINYRR)then
        vx = TINYRR
        vy = TINYRR
        vz = TINYRR
     else
        vx = X
        vy = Y
        vz = Z
     endif
     if(abs(X) < TINYRR .and. abs(Y) < TINYRR)then
        ux = - TINYRR
        uy = TINYRR
        uz = ZEO
        wx = TINYRR*Z
        wy = TINYRR*Z
        wz = -2*TINYRR*TINYRR
     else
        ux = - Y
        uy = X
        uz = ZEO
        wx = X*Z
        wy = Y*Z
        wz = -(X*X + Y*Y)
     endif

  f = 1.d0/chipn1

  fx = gxx*wx*wx + gyy*wy*wy + gzz*wz*wz &
     +(gxy*wx*wy + gxz*wx*wz + gyz*wy*wz)*TWO
  fx = dsqrt(fx*f)
  wx = wx/fx
  wy = wy/fx
  wz = wz/fx

  fx = gxx*wx*ux + gxy*wx*uy + gxz*wx*uz + &
       gxy*wy*ux + gyy*wy*uy + gyz*wy*uz + &
       gxz*wz*ux + gyz*wz*uy + gzz*wz*uz
  fx = fx*f
  ux = ux - fx*wx
  uy = uy - fx*wy
  uz = uz - fx*wz
  fx = gxx*ux*ux + gyy*uy*uy + gzz*uz*uz &
     +(gxy*ux*uy + gxz*ux*uz + gyz*uy*uz)*TWO
  fx = dsqrt(fx*f) 
  ux = ux/fx
  uy = uy/fx
  uz = uz/fx

  fx = gxx*vx*wx + gxy*vx*wy + gxz*vx*wz + &
       gxy*vy*wx + gyy*vy*wy + gyz*vy*wz + &
       gxz*vz*wx + gyz*vz*wy + gzz*vz*wz
  fx = fx*f       
  vx = vx - fx*wx
  vy = vy - fx*wy
  vz = vz - fx*wz
  fx = gxx*ux*vx + gxy*ux*vy + gxz*ux*vz + &
       gxy*uy*vx + gyy*uy*vy + gyz*uy*vz + &
       gxz*uz*vx + gyz*uz*vy + gzz*uz*vz
  fx = fx*f       
  vx = vx - fx*ux
  vy = vy - fx*uy
  vz = vz - fx*uz
  fx = gxx*vx*vx + gyy*vy*vy + gzz*vz*vz &
     +(gxy*vx*vy + gxz*vx*vz + gyz*vy*vz)*TWO
  fx = dsqrt(fx*f)
  vx = vx/fx
  vy = vy/fx
  vz = vz/fx
#elif (tetradtype == 2)  
     if(abs(X) < TINYRR .and. abs(Y) < TINYRR .and. abs(Z) < TINYRR)then
        vx = TINYRR
        vy = TINYRR
        vz = TINYRR
     else
        vx = X
        vy = Y
        vz = Z
     endif
     if(abs(X) < TINYRR .and. abs(Y) < TINYRR)then
        ux = - TINYRR
        uy = TINYRR
        uz = ZEO
        wx = TINYRR*Z
        wy = TINYRR*Z
        wz = -2*TINYRR*TINYRR
     else
        ux = - Y
        uy = X
        uz = ZEO
        wx = X*Z
        wy = Y*Z
        wz = -(X*X + Y*Y)
     endif

  fx = vx
  fy = vy
  fz = vz
  vx = gupxx*fx + gupxy*fy + gupxz*fz
  vy = gupxy*fx + gupyy*fy + gupyz*fz
  vz = gupxz*fx + gupyz*fy + gupzz*fz

  f = 1.d0/chipn1

  fx = gxx*vx*vx + gyy*vy*vy + gzz*vz*vz &
     +(gxy*vx*vy + gxz*vx*vz + gyz*vy*vz)*TWO
  fx = dsqrt(fx*f)
  vx = vx/fx
  vy = vy/fx
  vz = vz/fx

  fx = gxx*vx*ux + gxy*vx*uy + gxz*vx*uz + &
       gxy*vy*ux + gyy*vy*uy + gyz*vy*uz + &
       gxz*vz*ux + gyz*vz*uy + gzz*vz*uz
  fx = fx*f
  ux = ux - fx*vx
  uy = uy - fx*vy
  uz = uz - fx*vz
  fx = gxx*ux*ux + gyy*uy*uy + gzz*uz*uz &
     +(gxy*ux*uy + gxz*ux*uz + gyz*uy*uz)*TWO
  fx = dsqrt(fx*f) 
  ux = ux/fx
  uy = uy/fx
  uz = uz/fx

  fx = gxx*vx*wx + gxy*vx*wy + gxz*vx*wz + &
       gxy*vy*wx + gyy*vy*wy + gyz*vy*wz + &
       gxz*vz*wx + gyz*vz*wy + gzz*vz*wz
  fx = fx*f       
  wx = wx - fx*vx
  wy = wy - fx*vy
  wz = wz - fx*vz
  fx = gxx*ux*wx + gxy*ux*wy + gxz*ux*wz + &
       gxy*uy*wx + gyy*uy*wy + gyz*uy*wz + &
       gxz*uz*wx + gyz*uz*wy + gzz*uz*wz
  fx = fx*f       
  wx = wx - fx*ux
  wy = wy - fx*uy
  wz = wz - fx*uz
  fx = gxx*wx*wx + gyy*wy*wy + gzz*wz*wz &
     +(gxy*wx*wy + gxz*wx*wz + gyz*wy*wz)*TWO
  fx = dsqrt(fx*f)
  wx = wx/fx
  wy = wy/fx
  wz = wz/fx
#endif

  f = dsqrt(f)**3
! E_i
  HEx = (gxx*Ex+gxy*Ey+gxz*Ez)*chipn1
  HEy = (gxy*Ex+gyy*Ey+gyz*Ez)*chipn1
  HEz = (gxz*Ex+gyz*Ey+gzz*Ez)*chipn1

  Rphi1 = HEx*vx+HEy*vy+HEz*vz
! \sqrt(gamma)u x w (theta x phi)
  HBx = (uy*wz - uz*wy)*f
  HBy = (uz*wx - ux*wz)*f
  HBz = (ux*wy - uy*wx)*f
  Iphi1 = HBx*Bx+HBy*By+HBz*Bz

  Rphi1 = Rphi1/2.d0
  Iphi1 = Iphi1/2.d0

  return

  end subroutine getnpem1_point
