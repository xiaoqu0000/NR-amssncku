

!-----------------------------------------------------------------------------------
!
!Set up approximate puncture initial data for n black holes with lousto's
!formula PRD 77, 024034 (2008)
!
!-----------------------------------------------------------------------------------

  subroutine get_initial_nbhs_sh(ex,X,Y,Z, &
                         chi,  trK, &
                         gxx,  gxy,  gxz,  gyy,  gyz,  gzz,&
                         Axx,  Axy,  Axz,  Ayy,  Ayz,  Azz,&
                         Gmx,  Gmy,  Gmz,  &
                         Lap,  Sfx,  Sfy,  Sfz,&
                         dtSfx,dtSfy,dtSfz,Mass,Porg,Pmom,Spin,N)

  implicit none

!------= input arguments
 
  integer,intent(in) :: N
  integer, dimension(3),  intent(in) :: ex
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in) :: X,Y,Z
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: chi
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: trK,Lap,Sfx,Sfy,Sfz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: Gmx,Gmy,Gmz
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out) :: dtSfx,dtSfy,dtSfz
  real*8,  dimension(N), intent(in) :: Mass
  real*8,  dimension(3*N), intent(in) :: Porg,Pmom,Spin

!------= local variables
  real*8,dimension(ex(1),ex(2),ex(3))::psi
  integer :: i,j,k,bhi
  real*8 :: M,Px,Py,Pz,PP,Sx,Sy,Sz,SS
  real*8 :: nx,ny,nz,rr,tmp
  real*8 :: u,u1,u2,u3,u4,u5
  real*8 :: mup,mus,b,ell
  real*8,  parameter :: HLF = 5.d-1, ZEO = 0.d0, ONE = 1.d0, THR = 3.d0
  real*8,parameter::TINYRR=1.d-14

  do k = 1,ex(3)
   do j = 1,ex(2)
    do i = 1,ex(1)
! black hole 1
      M = mass(1)
      nx = x(i,j,k) - Porg(1)
      ny = y(i,j,k) - Porg(2)
      nz = z(i,j,k) - Porg(3)
      Px = Pmom(1)
      Py = Pmom(2)
      Pz = Pmom(3)
      Sx = Spin(1)
      Sy = Spin(2)
      Sz = Spin(3)

      rr = dsqrt(nx*nx+ny*ny+nz*nz)
      if(rr.lt.TINYRR) rr=(X(2,1,1)-X(1,1,1))/2.d0
      nx = nx / rr
      ny = ny / rr
      nz = nz / rr
      PP = dsqrt(Px**2 + Py**2 + Pz**2)
      if(PP .gt. 0.d0) then
        mup = (Px*nx+Py*ny+Pz*nz)/PP
      else
        mup = 0.0
      endif
      SS = dsqrt(Sx**2 + Sy**2 + Sz**2)
      if(SS .gt. 0.d0) then
        mus = (Sx*nx+Sy*ny+Sz*nz)/SS
      else
        mus = 0.0
      endif
      b = 2.d0*rr/M
      ell = 1.d0/(1.d0+b)

      u1 = 5.d0/8.d0*ell*(1.d0-2.d0*ell+2.d0*ell**2-ell**3+ell**4/5.d0)       
      u2 = (1.5d1+1.17d2*ell-7.9d1*ell**2+4.3d1*ell**3-1.4d1*ell**4+2.d0*ell**5 &
          +8.4d1*dlog(ell)/b)/4.d1/b**2
      u3 = ell+ell**2+ell**3-4.d0*ell**4+2.d0*ell**5
      u4 = -b**2*ell**5
      u5 = b*(1.d0+5.d0*b+1.d1*b**2)*ell**5/8.d1

      tmp = (Py*Sz-Pz*Sy)*nx + (Pz*Sx-Px*Sz)*ny + (Px*Sy-Py*Sx)*nz

      u = PP**2/M**2*(u1 + u2*(3.d0*mup**2-ONE)) + &
          2.d0*SS**2/5.d0/M**4*(u3 + u4*(3.d0*mus**2-ONE)) + &
          u5*tmp

     psi(i,j,k) = ONE + u + HLF*M/rr

     tmp = Px * nx + Py * ny + Pz * nz

     Axx(i,j,k) = (HLF *( Px * nx + nx * Px - ( ONE - nx * nx )* tmp ) + &
                   ( nx * Sy * nz - nx * Sz * ny + nx * Sy * nz - nx * Sz * ny ) / rr ) * &
                   THR / ( rr * rr )

     Ayy(i,j,k) = (HLF *( Py * ny + ny * Py - ( ONE - ny * ny )* tmp ) + &
                   ( ny * Sz * nx - ny * Sx * nz + ny * Sz * nx - ny * Sx * nz ) / rr ) * &
                   THR / ( rr * rr )

     Azz(i,j,k) = (HLF *( Pz * nz + nz * Pz - ( ONE - nz * nz )* tmp ) + &
                   ( nz * Sx * ny - nz * Sy * nx + nz * Sx * ny - nz * Sy * nx ) / rr ) * &
                   THR / ( rr * rr )

     Axy(i,j,k) = (HLF *( Px * ny + nx * Py + nx * ny * tmp ) + &
                   ( nx * Sz * nx - nx * Sx * nz + ny * Sy * nz - ny * Sz * ny ) / rr ) * &
                   THR / ( rr * rr )

     Axz(i,j,k) = (HLF *( Px * nz + nx * Pz + nx * nz * tmp ) + &
                   ( nx * Sx * ny - nx * Sy * nx + nz * Sy * nz - nz * Sz * ny ) / rr ) * &
                   THR / ( rr * rr )

     Ayz(i,j,k) = (HLF *( Py * nz + ny * Pz + ny * nz * tmp ) + &
                   ( ny * Sx * ny - ny * Sy * nx + nz * Sz * nx - nz * Sx * nz ) / rr ) * &
                   THR / ( rr * rr )
! black hole 2 and 3, ...
     do bhi=2,N
      M = Mass(bhi)
      nx = x(i,j,k) - Porg(3*(bhi-1)+1)
      ny = y(i,j,k) - Porg(3*(bhi-1)+2)
      nz = z(i,j,k) - Porg(3*(bhi-1)+3)
      Px = Pmom(3*(bhi-1)+1)
      Py = Pmom(3*(bhi-1)+2)
      Pz = Pmom(3*(bhi-1)+3)
      Sx = Spin(3*(bhi-1)+1)
      Sy = Spin(3*(bhi-1)+2)
      Sz = Spin(3*(bhi-1)+3)

      rr = dsqrt(nx*nx+ny*ny+nz*nz)
      if(rr.lt.TINYRR) rr=(X(2,1,1)-X(1,1,1))/2.d0
      nx = nx / rr
      ny = ny / rr
      nz = nz / rr
      PP = dsqrt(Px**2 + Py**2 + Pz**2)
      if(PP .gt. 0.d0) then
        mup = (Px*nx+Py*ny+Pz*nz)/PP
      else
        mup = 0.0
      endif
      SS = dsqrt(Sx**2 + Sy**2 + Sz**2)
      if(SS .gt. 0.d0) then
        mus = (Sx*nx+Sy*ny+Sz*nz)/SS
      else
        mus = 0.0
      endif
      b = 2.d0*rr/M
      ell = 1.d0/(1.d0+b)

      u1 = 5.d0/8.d0*ell*(1.d0-2.d0*ell+2.d0*ell**2-ell**3+ell**4/5.d0)       
      u2 = (1.5d1+1.17d2*ell-7.9d1*ell**2+4.3d1*ell**3-1.4d1*ell**4+2.d0*ell**5 &
          +8.4d1*dlog(ell)/b)/4.d1/b**2
      u3 = ell+ell**2+ell**3-4.d0*ell**4+2.d0*ell**5
      u4 = -b**2*ell**5
      u5 = b*(1.d0+5.d0*b+1.d1*b**2)*ell**5/8.d1

      tmp = (Py*Sz-Pz*Sy)*nx + (Pz*Sx-Px*Sz)*ny + (Px*Sy-Py*Sx)*nz

      u = PP**2/M**2*(u1 + u2*(3.d0*mup**2-ONE)) + &
          2.d0*SS**2/5.d0/M**4*(u3 + u4*(3.d0*mus**2-ONE)) + &
          u5*tmp

     psi(i,j,k) = psi(i,j,k) + u + HLF*M/rr

     tmp = Px * nx + Py * ny + Pz * nz

     Axx(i,j,k) = Axx(i,j,k) + &
                  (HLF *( Px * nx + nx * Px - ( ONE - nx * nx )* tmp ) + &
                  ( nx * Sy * nz - nx * Sz * ny + nx * Sy * nz - nx * Sz * ny ) / rr ) * &
                  THR / ( rr * rr )

     Ayy(i,j,k) = Ayy(i,j,k) + &
                  (HLF *( Py * ny + ny * Py - ( ONE - ny * ny )* tmp ) + &
                  ( ny * Sz * nx - ny * Sx * nz + ny * Sz * nx - ny * Sx * nz ) / rr ) * &
                  THR / ( rr * rr )

     Azz(i,j,k) = Azz(i,j,k) + &
                  (HLF *( Pz * nz + nz * Pz - ( ONE - nz * nz )* tmp ) + &
                  ( nz * Sx * ny - nz * Sy * nx + nz * Sx * ny - nz * Sy * nx ) / rr ) * &
                  THR / ( rr * rr )

     Axy(i,j,k) = Axy(i,j,k) + &
                  (HLF *( Px * ny + nx * Py + nx * ny * tmp ) + &
                  ( nx * Sz * nx - nx * Sx * nz + ny * Sy * nz - ny * Sz * ny ) / rr ) * &
                  THR / ( rr * rr )

     Axz(i,j,k) = Axz(i,j,k) + &
                  (HLF *( Px * nz + nx * Pz + nx * nz * tmp ) + &
                  ( nx * Sx * ny - nx * Sy * nx + nz * Sy * nz - nz * Sz * ny ) / rr ) * &
                  THR / ( rr * rr )

     Ayz(i,j,k) = Ayz(i,j,k) + &
                  (HLF *( Py * nz + ny * Pz + ny * nz * tmp ) + &
                  ( ny * Sx * ny - ny * Sy * nx + nz * Sz * nx - nz * Sx * nz ) / rr ) * &
                  THR / ( rr * rr )
      enddo
    enddo
   enddo
  enddo

  chi = ONE / psi **4 - ONE

  Lap = ONE / ( psi * psi ) - ONE

!~~~~~~ tilde Aij = Aij / Psi^6
  psi = psi * psi * psi * psi * psi * psi

  Axx = Axx / psi
  Ayy = Ayy / psi
  Azz = Azz / psi
  Axy = Axy / psi
  Axz = Axz / psi
  Ayz = Ayz / psi

  gxx = ZEO
  gyy = ZEO
  gzz = ZEO
  gxy = ZEO
  gxz = ZEO
  gyz = ZEO

  trK = ZEO

  Gmx = ZEO
  Gmy = ZEO
  Gmz = ZEO

  Sfx = ZEO
  Sfy = ZEO
  Sfz = ZEO

  dtSfx = ZEO
  dtSfy = ZEO
  dtSfz = ZEO
  
  return

  end subroutine get_initial_nbhs_sh
!----------------------------------------------------  
! I use this routine to unify the parameters
subroutine shellcordpar(A,B,r0,eps)
implicit none
! argument variables
double precision,intent(out)::A,B,r0,eps

A=1.d0
B=0.d0
r0=0.d0
eps=1.d0

return

end subroutine shellcordpar
!----------------------------------------------------  
! R = f(r)-f(0)
subroutine getcartr(ex,R,cartr)
implicit none
! argument variables
integer,intent(in)::ex
double precision,intent(in),dimension(ex)::R
double precision,intent(out),dimension(ex)::cartr

!~~~~~~ local variables
double precision,dimension(ex)::f
double precision :: A,B,r0,eps

call shellcordpar(A,B,r0,eps)
f = R+B
cartr = r0+(A*f-B*dsqrt(A*A+(f*f-B*B)/eps))/(A*A-B*B/eps)

return
end subroutine getcartr
! dR/dr = ...
subroutine getdRdcartr(ex,R,dRdcartr)
implicit none
! argument variables
integer,intent(in)::ex
double precision,intent(in),dimension(ex)::R
double precision,intent(out),dimension(ex)::dRdcartr

!~~~~~~ local variables
double precision,dimension(ex)::cartr
double precision :: A,B,r0,eps

  call shellcordpar(A,B,r0,eps)

  call getcartr(ex,R,cartr)
  dRdcartr = A + B*(cartr-r0)/dsqrt(eps*eps+eps*(cartr-r0)*(cartr-r0))

return
end subroutine getdRdcartr
! dR/drdr = ...
subroutine getdRdcartrcartr(ex,R,dRdcartrcartr)
implicit none
! argument variables
integer,intent(in)::ex
double precision,intent(in),dimension(ex)::R
double precision,intent(out),dimension(ex)::dRdcartrcartr

!~~~~~~ local variables
double precision,dimension(ex)::cartr

double precision :: A,B,r0,eps

  call shellcordpar(A,B,r0,eps)

  call getcartr(ex,R,cartr)
  dRdcartrcartr = B*dsqrt(eps)/(dsqrt(eps+(cartr-r0)*(cartr-r0)))**3

  return

end subroutine getdRdcartrcartr

subroutine zp_getxyz(ex,rho,sigma,R,x,y,z)

implicit none
! argument variables
integer,dimension(3),intent(in)::ex
double precision,intent(in),dimension(ex(1))::rho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::R
double precision,intent(out),dimension(ex(1),ex(2),ex(3))::x,y,z
!~~~~~~ other variables
double precision,dimension(ex(3))::cartr
double precision,dimension(ex(1))::tgrho
double precision,dimension(ex(2))::tgsigma
integer :: i,j,k

  call getcartr(ex(3),R,cartr)
  tgrho = dtan(rho)
  tgsigma = dtan(sigma)

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
    z(i,j,k) = cartr(k)/dsqrt(1+tgrho(i)*tgrho(i)+tgsigma(j)*tgsigma(j))
    x(i,j,k) = z(i,j,k)*tgrho(i)
    y(i,j,k) = z(i,j,k)*tgsigma(j)
  enddo
  enddo
  enddo

  return

end subroutine zp_getxyz

subroutine zm_getxyz(ex,rho,sigma,R,x,y,z)

implicit none
! argument variables
integer,dimension(3),intent(in)::ex
double precision,intent(in),dimension(ex(1))::rho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::R
double precision,intent(out),dimension(ex(1),ex(2),ex(3))::x,y,z
!~~~~~~ other variables
double precision,dimension(ex(3))::cartr
double precision,dimension(ex(1))::tgrho
double precision,dimension(ex(2))::tgsigma
integer :: i,j,k

  call getcartr(ex(3),R,cartr)
  tgrho = dtan(rho)
  tgsigma = dtan(sigma)

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
    z(i,j,k) = -cartr(k)/dsqrt(1+tgrho(i)*tgrho(i)+tgsigma(j)*tgsigma(j))
    x(i,j,k) = z(i,j,k)*tgrho(i)
    y(i,j,k) = z(i,j,k)*tgsigma(j)
  enddo
  enddo
  enddo

  return

end subroutine zm_getxyz

subroutine yp_getxyz(ex,rho,sigma,R,x,y,z)

implicit none
! argument variables
integer,dimension(3),intent(in)::ex
double precision,intent(in),dimension(ex(1))::rho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::R
double precision,intent(out),dimension(ex(1),ex(2),ex(3))::x,y,z
!~~~~~~ other variables
double precision,dimension(ex(3))::cartr
double precision,dimension(ex(1))::tgrho
double precision,dimension(ex(2))::tgsigma
integer :: i,j,k

  call getcartr(ex(3),R,cartr)
  tgrho = dtan(rho)
  tgsigma = dtan(sigma)

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
    y(i,j,k) = cartr(k)/dsqrt(1+tgrho(i)*tgrho(i)+tgsigma(j)*tgsigma(j))
    x(i,j,k) = y(i,j,k)*tgrho(i)
    z(i,j,k) = y(i,j,k)*tgsigma(j)
  enddo
  enddo
  enddo

  return

end subroutine yp_getxyz

subroutine ym_getxyz(ex,rho,sigma,R,x,y,z)

implicit none
! argument variables
integer,dimension(3),intent(in)::ex
double precision,intent(in),dimension(ex(1))::rho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::R
double precision,intent(out),dimension(ex(1),ex(2),ex(3))::x,y,z
!~~~~~~ other variables
double precision,dimension(ex(3))::cartr
double precision,dimension(ex(1))::tgrho
double precision,dimension(ex(2))::tgsigma
integer :: i,j,k

  call getcartr(ex(3),R,cartr)
  tgrho = dtan(rho)
  tgsigma = dtan(sigma)

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
    y(i,j,k) = -cartr(k)/dsqrt(1+tgrho(i)*tgrho(i)+tgsigma(j)*tgsigma(j))
    x(i,j,k) = y(i,j,k)*tgrho(i)
    z(i,j,k) = y(i,j,k)*tgsigma(j)
  enddo
  enddo
  enddo

  return

end subroutine ym_getxyz

subroutine xp_getxyz(ex,rho,sigma,R,x,y,z)

implicit none
! argument variables
integer,dimension(3),intent(in)::ex
double precision,intent(in),dimension(ex(1))::rho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::R
double precision,intent(out),dimension(ex(1),ex(2),ex(3))::x,y,z
!~~~~~~ other variables
double precision,dimension(ex(3))::cartr
double precision,dimension(ex(1))::tgrho
double precision,dimension(ex(2))::tgsigma
integer :: i,j,k

  call getcartr(ex(3),R,cartr)
  tgrho = dtan(rho)
  tgsigma = dtan(sigma)

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
    x(i,j,k) = cartr(k)/dsqrt(1+tgrho(i)*tgrho(i)+tgsigma(j)*tgsigma(j))
    y(i,j,k) = x(i,j,k)*tgrho(i)
    z(i,j,k) = x(i,j,k)*tgsigma(j)
  enddo
  enddo
  enddo

  return

end subroutine xp_getxyz

subroutine xm_getxyz(ex,rho,sigma,R,x,y,z)

implicit none
! argument variables
integer,dimension(3),intent(in)::ex
double precision,intent(in),dimension(ex(1))::rho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::R
double precision,intent(out),dimension(ex(1),ex(2),ex(3))::x,y,z
!~~~~~~ other variables
double precision,dimension(ex(3))::cartr
double precision,dimension(ex(1))::tgrho
double precision,dimension(ex(2))::tgsigma
integer :: i,j,k

  call getcartr(ex(3),R,cartr)
  tgrho = dtan(rho)
  tgsigma = dtan(sigma)

  do k=1,ex(3)
  do j=1,ex(2)
  do i=1,ex(1)
    x(i,j,k) = -cartr(k)/dsqrt(1+tgrho(i)*tgrho(i)+tgsigma(j)*tgsigma(j))
    y(i,j,k) = x(i,j,k)*tgrho(i)
    z(i,j,k) = x(i,j,k)*tgsigma(j)
  enddo
  enddo
  enddo

  return

end subroutine xm_getxyz
!------------------------------------------------------------------------------------------
! calculate Jacobians
subroutine xpm_getjacobian(ex,rho,sigma,R,x,y,z,                        &
          drhodx, drhody, drhodz,                                      &
          dsigmadx,dsigmady,dsigmadz,                                  &
          dRdx,dRdy,dRdz,                                              &
          drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,             &
          dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz, &
          dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz)

implicit none
! argument variables
integer,dimension(3),intent(in)::ex
double precision,intent(in),dimension(ex(1))::rho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::R
double precision,intent(in),dimension(ex(1),ex(2),ex(3))::x,y,z
double precision,intent(out),dimension(ex(1),ex(2),ex(3))::drhodx, drhody, drhodz
double precision,intent(out),dimension(ex(1),ex(2),ex(3))::dsigmadx,dsigmady,dsigmadz
double precision,intent(out),dimension(ex(1),ex(2),ex(3))::dRdx,dRdy,dRdz
double precision,intent(out),dimension(ex(1),ex(2),ex(3))::drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz
double precision,intent(out),dimension(ex(1),ex(2),ex(3))::dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz
double precision,intent(out),dimension(ex(1),ex(2),ex(3))::dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz
!~~~~~~ other variables
double precision,dimension(ex(3))::cartr
double precision,dimension(ex(1),ex(2),ex(3))::srt,xxyy,xxzz,dRdcartr,dRdcartrcartr
integer :: i,j,k
real*8,parameter :: ZEO=0.d0

  xxyy = x*x + y*y
  xxzz = x*x + z*z
  srt = dsqrt(xxyy + z*z)
  call getdRdcartr(ex(3),R,dRdcartr(1,1,:))
  call getdRdcartrcartr(ex(3),R,dRdcartrcartr(1,1,:))
  do k=1,ex(3)
    dRdcartr(:,:,k) = dRdcartr(1,1,k)
    dRdcartrcartr(:,:,k) = dRdcartrcartr(1,1,k)
  enddo

  dRdx = x/srt*dRdcartr
  dRdy = y/srt*dRdcartr
  dRdz = z/srt*dRdcartr
  drhodx = -y/xxyy
  drhody = x/xxyy
  drhodz = ZEO
  dsigmadx = -z/xxzz
  dsigmady = ZEO
  dsigmadz = x/xxzz

  dRdxx = dRdcartrcartr*x*x/srt/srt+dRdcartr*(y*y+z*z)/srt**3
  dRdxy = dRdcartrcartr*x*y/srt/srt-dRdcartr*(    x*y)/srt**3
  dRdxz = dRdcartrcartr*x*z/srt/srt-dRdcartr*(    x*z)/srt**3
  dRdyy = dRdcartrcartr*y*y/srt/srt+dRdcartr*(x*x+z*z)/srt**3
  dRdyz = dRdcartrcartr*y*z/srt/srt-dRdcartr*(    y*z)/srt**3
  dRdzz = dRdcartrcartr*z*z/srt/srt+dRdcartr*(x*x+y*y)/srt**3
  drhodxx = 2*x*y/xxyy**2
  drhodxy = (-x*x + y*y)/xxyy**2
  drhodxz = ZEO
  drhodyy = -drhodxx
  drhodyz = ZEO
  drhodzz = ZEO
  dsigmadxx = (2*x*z)/xxzz**2
  dsigmadxy = ZEO
  dsigmadxz = (-x*x + z*z)/xxzz**2
  dsigmadyy = ZEO
  dsigmadyz = ZEO
  dsigmadzz = -dsigmadxx

  return

end subroutine xpm_getjacobian
!~~~~
subroutine ypm_getjacobian(ex,rho,sigma,R,x,y,z,                        &
          drhodx, drhody, drhodz,                                      &
          dsigmadx,dsigmady,dsigmadz,                                  &
          dRdx,dRdy,dRdz,                                              &
          drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,             &
          dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz, &
          dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz)

implicit none
! argument variables
integer,dimension(3),intent(in)::ex
double precision,intent(in),dimension(ex(1))::rho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::R
double precision,intent(in),dimension(ex(1),ex(2),ex(3))::x,y,z
double precision,intent(out),dimension(ex(1),ex(2),ex(3))::drhodx, drhody, drhodz
double precision,intent(out),dimension(ex(1),ex(2),ex(3))::dsigmadx,dsigmady,dsigmadz
double precision,intent(out),dimension(ex(1),ex(2),ex(3))::dRdx,dRdy,dRdz
double precision,intent(out),dimension(ex(1),ex(2),ex(3))::drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz
double precision,intent(out),dimension(ex(1),ex(2),ex(3))::dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz
double precision,intent(out),dimension(ex(1),ex(2),ex(3))::dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz
!~~~~~~ other variables
double precision,dimension(ex(3))::cartr
double precision,dimension(ex(1),ex(2),ex(3))::srt,xxyy,yyzz,dRdcartr,dRdcartrcartr
integer :: i,j,k
real*8,parameter :: ZEO=0.d0

  xxyy = x*x + y*y
  yyzz = y*y + z*z
  srt = dsqrt(xxyy + z*z)
  call getdRdcartr(ex(3),R,dRdcartr(1,1,:))
  call getdRdcartrcartr(ex(3),R,dRdcartrcartr(1,1,:))
  do k=1,ex(3)
    dRdcartr(:,:,k) = dRdcartr(1,1,k)
    dRdcartrcartr(:,:,k) = dRdcartrcartr(1,1,k)
  enddo

  dRdx = x/srt*dRdcartr
  dRdy = y/srt*dRdcartr
  dRdz = z/srt*dRdcartr
  drhodx = y/xxyy
  drhody = -x/xxyy
  drhodz = ZEO
  dsigmadx = ZEO
  dsigmady = -z/yyzz
  dsigmadz = y/yyzz

  dRdxx = dRdcartrcartr*x*x/srt/srt+dRdcartr*(y*y+z*z)/srt**3
  dRdxy = dRdcartrcartr*x*y/srt/srt-dRdcartr*(    x*y)/srt**3
  dRdxz = dRdcartrcartr*x*z/srt/srt-dRdcartr*(    x*z)/srt**3
  dRdyy = dRdcartrcartr*y*y/srt/srt+dRdcartr*(x*x+z*z)/srt**3
  dRdyz = dRdcartrcartr*y*z/srt/srt-dRdcartr*(    y*z)/srt**3
  dRdzz = dRdcartrcartr*z*z/srt/srt+dRdcartr*(x*x+y*y)/srt**3
  drhodxx = -2*x*y/xxyy**2
  drhodxy = (x*x - y*y)/xxyy**2
  drhodxz = ZEO
  drhodyy = -drhodxx
  drhodyz = ZEO
  drhodzz = ZEO
  dsigmadxx = ZEO
  dsigmadxy = ZEO
  dsigmadxz = ZEO
  dsigmadyy = (2*y*z)/yyzz**2
  dsigmadyz = (-y*y + z*z)/yyzz**2
  dsigmadzz = -dsigmadyy

  return

end subroutine ypm_getjacobian
!~~~~
subroutine zpm_getjacobian(ex,rho,sigma,R,x,y,z,                        &
          drhodx, drhody, drhodz,                                      &
          dsigmadx,dsigmady,dsigmadz,                                  &
          dRdx,dRdy,dRdz,                                              &
          drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,             &
          dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz, &
          dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz)

implicit none
! argument variables
integer,dimension(3),intent(in)::ex
double precision,intent(in),dimension(ex(1))::rho
double precision,intent(in),dimension(ex(2))::sigma
double precision,intent(in),dimension(ex(3))::R
double precision,intent(in),dimension(ex(1),ex(2),ex(3))::x,y,z
double precision,intent(out),dimension(ex(1),ex(2),ex(3))::drhodx, drhody, drhodz
double precision,intent(out),dimension(ex(1),ex(2),ex(3))::dsigmadx,dsigmady,dsigmadz
double precision,intent(out),dimension(ex(1),ex(2),ex(3))::dRdx,dRdy,dRdz
double precision,intent(out),dimension(ex(1),ex(2),ex(3))::drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz
double precision,intent(out),dimension(ex(1),ex(2),ex(3))::dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz
double precision,intent(out),dimension(ex(1),ex(2),ex(3))::dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz
!~~~~~~ other variables
double precision,dimension(ex(3))::cartr
double precision,dimension(ex(1),ex(2),ex(3))::srt,xxzz,yyzz,dRdcartr,dRdcartrcartr
integer :: i,j,k
real*8,parameter :: ZEO=0.d0

  xxzz = x*x + z*z
  yyzz = y*y + z*z
  srt = dsqrt(xxzz + y*y)
  call getdRdcartr(ex(3),R,dRdcartr(1,1,:))
  call getdRdcartrcartr(ex(3),R,dRdcartrcartr(1,1,:))
  do k=1,ex(3)
    dRdcartr(:,:,k) = dRdcartr(1,1,k)
    dRdcartrcartr(:,:,k) = dRdcartrcartr(1,1,k)
  enddo

  dRdx = x/srt*dRdcartr
  dRdy = y/srt*dRdcartr
  dRdz = z/srt*dRdcartr
  drhodx = z/xxzz
  drhody = ZEO
  drhodz = -x/xxzz
  dsigmadx = ZEO
  dsigmady = z/yyzz
  dsigmadz = -y/yyzz

  dRdxx = dRdcartrcartr*x*x/srt/srt+dRdcartr*(y*y+z*z)/srt**3
  dRdxy = dRdcartrcartr*x*y/srt/srt-dRdcartr*(    x*y)/srt**3
  dRdxz = dRdcartrcartr*x*z/srt/srt-dRdcartr*(    x*z)/srt**3
  dRdyy = dRdcartrcartr*y*y/srt/srt+dRdcartr*(x*x+z*z)/srt**3
  dRdyz = dRdcartrcartr*y*z/srt/srt-dRdcartr*(    y*z)/srt**3
  dRdzz = dRdcartrcartr*z*z/srt/srt+dRdcartr*(x*x+y*y)/srt**3
  drhodxx = -2*x*z/xxzz**2
  drhodxy = ZEO
  drhodxz = (x*x - z*z)/xxzz**2
  drhodyy = ZEO
  drhodyz = ZEO
  drhodzz = -drhodxx
  dsigmadxx = ZEO
  dsigmadxy = ZEO
  dsigmadxz = ZEO
  dsigmadyy = -(2*y*z)/yyzz**2
  dsigmadyz = (y*y - z*z)/yyzz**2
  dsigmadzz = -dsigmadyy

  return

end subroutine zpm_getjacobian
