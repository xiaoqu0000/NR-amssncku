
!-----------------------------------------------------------------------------------
!
!Set up approximate puncture initial data for n charged black holes
!PRD 80, 104022
!-----------------------------------------------------------------------------------

  subroutine get_initial_nbhsem(ext,X,Y,Z, &
                         chi,  trK, &
                         gxx,  gxy,  gxz,  gyy,  gyz,  gzz,&
                         Axx,  Axy,  Axz,  Ayy,  Ayz,  Azz,&
                         Gmx,  Gmy,  Gmz,  &
                         Lap,  Sfx,  Sfy,  Sfz,&
                         dtSfx,dtSfy,dtSfz,&
                         Ex,Ey,Ez,Bx,By,Bz,Kpsi,Kphi,&
                         Mass,Qchar,Porg,Pmom,Spin,N)

  implicit none

!------= input arguments
 
  integer,intent(in) :: N
  integer, dimension(3),  intent(in) :: ext
  real*8,  dimension(ext(1)), intent(in) :: X
  real*8,  dimension(ext(2)), intent(in) :: Y
  real*8,  dimension(ext(3)), intent(in) :: Z
  real*8,  dimension(ext(1),ext(2),ext(3)), intent(out) :: chi
  real*8,  dimension(ext(1),ext(2),ext(3)), intent(out) :: Ex,Ey,Ez,Bx,By,Bz,Kpsi,Kphi
  real*8,  dimension(ext(1),ext(2),ext(3)), intent(out) :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8,  dimension(ext(1),ext(2),ext(3)), intent(out) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8,  dimension(ext(1),ext(2),ext(3)), intent(out) :: trK,Lap,Sfx,Sfy,Sfz
  real*8,  dimension(ext(1),ext(2),ext(3)), intent(out) :: Gmx,Gmy,Gmz
  real*8,  dimension(ext(1),ext(2),ext(3)), intent(out) :: dtSfx,dtSfy,dtSfz
  real*8,  dimension(N), intent(in) :: Mass,Qchar
  real*8,  dimension(3*N), intent(in) :: Porg,Pmom,Spin

!------= local variables
  real*8,dimension(ext(1),ext(2),ext(3))::psi,phi
  integer :: i,j,k,bhi
  real*8 :: M,Q,Px,Py,Pz,PP,Sx,Sy,Sz,SS
  real*8 :: nx,ny,nz,rr,tmp
  real*8 :: u,u1,u2,u3,u4
  real*8 :: mup,mus,b,ell
  real*8,  parameter :: HLF = 5.d-1, ZEO = 0.d0, ONE = 1.d0, THR = 3.d0,FOUR=4.d0
  real*8,parameter::TINYRR=1.d-14
!sanity check: M/Q = constant
  M = mass(1)
  Q = Qchar(1)
  u1 = M/Q
  u2 = M/Q
  do bhi=2,N
     M = mass(bhi)
     Q = Qchar(bhi)
     u1 = min(u1,M/Q)
     u2 = max(u2,M/Q)
  enddo
  if(u2-u1.gt.TINYRR)then
      write(*,*)"error in initial_punctureem.f90: get_initial_nbhsem; we need constant Mi/Qi, but"
      write(*,*)"Mass = ",mass
      write(*,*)"Qchar = ",Qchar
      stop
   endif

  do k = 1,ext(3)
   do j = 1,ext(2)
    do i = 1,ext(1)
! black hole 1
      M = mass(1)
      Q = Qchar(1)
      nx = x(i) - Porg(1)
      ny = y(j) - Porg(2)
      nz = z(k) - Porg(3)
      Px = Pmom(1)
      Py = Pmom(2)
      Pz = Pmom(3)
      Sx = Spin(1)
      Sy = Spin(2)
      Sz = Spin(3)

      rr = dsqrt(nx*nx+ny*ny+nz*nz)
      if(rr.lt.TINYRR) rr=(X(2)-X(1))/2.d0
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
     u3 = ell/2.d1*(1.d0+ell+ell**2-4.d0*ell**3+2.d0*ell**4)
     u4 = ell**2/1.d1*(1.d1-2.5d1*ell+2.1d1*ell**2-6.d0*ell**3)

     tmp = (Py*Sz-Pz*Sy)*nx + (Pz*Sx-Px*Sz)*ny + (Px*Sy-Py*Sx)*nz

     u = PP**2/M**2*(u1 + u2*(3.d0*mup**2-ONE)) + &
         6.d0*u3/M**4*SS**2*(1.d0+mus**2) + u4/M**3*tmp

     psi(i,j,k) = ONE + u + HLF*M/rr
     phi(i,j,k) = Q/rr

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

     Ex(i,j,k) = Q*nx/rr/rr                
     Ey(i,j,k) = Q*ny/rr/rr                
     Ez(i,j,k) = Q*nz/rr/rr                
! black hole 2 and 3, ...
     do bhi=2,N
      M = Mass(bhi)
      Q = Qchar(bhi)
      nx = x(i) - Porg(3*(bhi-1)+1)
      ny = y(j) - Porg(3*(bhi-1)+2)
      nz = z(k) - Porg(3*(bhi-1)+3)
      Px = Pmom(3*(bhi-1)+1)
      Py = Pmom(3*(bhi-1)+2)
      Pz = Pmom(3*(bhi-1)+3)
      Sx = Spin(3*(bhi-1)+1)
      Sy = Spin(3*(bhi-1)+2)
      Sz = Spin(3*(bhi-1)+3)

      rr = dsqrt(nx*nx+ny*ny+nz*nz)
      if(rr.lt.TINYRR) rr=(X(2)-X(1))/2.d0
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
     u3 = ell/2.d1*(1.d0+ell+ell**2-4.d0*ell**3+2.d0*ell**4)
     u4 = ell**2/1.d1*(1.d1-2.5d1*ell+2.1d1*ell**2-6.d0*ell**3)

     tmp = (Py*Sz-Pz*Sy)*nx + (Pz*Sx-Px*Sz)*ny + (Px*Sy-Py*Sx)*nz

     u = PP**2/M**2*(u1 + u2*(3.d0*mup**2-ONE)) + &
         6.d0*u3/M**4*SS**2*(1.d0+mus**2) + u4/M**3*tmp

     psi(i,j,k) = psi(i,j,k) + u + HLF*M/rr
     phi(i,j,k) = phi(i,j,k) + Q/rr

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

     Ex(i,j,k) = Ex(i,j,k) + Q*nx/rr/rr                
     Ey(i,j,k) = Ey(i,j,k) + Q*ny/rr/rr                
     Ez(i,j,k) = Ez(i,j,k) + Q*nz/rr/rr  
      enddo
    enddo
   enddo
  enddo

  psi = dsqrt(psi**2 - phi*phi/FOUR)
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

  Ex = Ex   / psi
  Ey = Ey   / psi
  Ez = Ez   / psi

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
 
  Bx = ZEO
  By = ZEO
  Bz = ZEO

  Kpsi = ZEO
  Kphi = ZEO

  return

  end subroutine get_initial_nbhsem
!-----------------------------------------------------------------------------------
!
!Set up approximate puncture initial data for n charged black holes
!PRD 80, 104022
! for shell
!-----------------------------------------------------------------------------------

  subroutine get_initial_nbhsem_ss(ext,X,Y,Z, &
                         chi,  trK, &
                         gxx,  gxy,  gxz,  gyy,  gyz,  gzz,&
                         Axx,  Axy,  Axz,  Ayy,  Ayz,  Azz,&
                         Gmx,  Gmy,  Gmz,  &
                         Lap,  Sfx,  Sfy,  Sfz,&
                         dtSfx,dtSfy,dtSfz,&
                         Ex,Ey,Ez,Bx,By,Bz,Kpsi,Kphi,&
                         Mass,Qchar,Porg,Pmom,Spin,N)

  implicit none

!------= input arguments
 
  integer,intent(in) :: N
  integer, dimension(3),  intent(in) :: ext
  real*8,  dimension(ext(1),ext(2),ext(3)), intent(in) :: X,Y,Z
  real*8,  dimension(ext(1),ext(2),ext(3)), intent(out) :: chi
  real*8,  dimension(ext(1),ext(2),ext(3)), intent(out) :: Ex,Ey,Ez,Bx,By,Bz,Kpsi,Kphi
  real*8,  dimension(ext(1),ext(2),ext(3)), intent(out) :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8,  dimension(ext(1),ext(2),ext(3)), intent(out) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8,  dimension(ext(1),ext(2),ext(3)), intent(out) :: trK,Lap,Sfx,Sfy,Sfz
  real*8,  dimension(ext(1),ext(2),ext(3)), intent(out) :: Gmx,Gmy,Gmz
  real*8,  dimension(ext(1),ext(2),ext(3)), intent(out) :: dtSfx,dtSfy,dtSfz
  real*8,  dimension(N), intent(in) :: Mass,Qchar
  real*8,  dimension(3*N), intent(in) :: Porg,Pmom,Spin

!------= local variables
  real*8,dimension(ext(1),ext(2),ext(3))::psi,phi
  integer :: i,j,k,bhi
  real*8 :: M,Q,Px,Py,Pz,PP,Sx,Sy,Sz,SS
  real*8 :: nx,ny,nz,rr,tmp
  real*8 :: u,u1,u2,u3,u4
  real*8 :: mup,mus,b,ell
  real*8,  parameter :: HLF = 5.d-1, ZEO = 0.d0, ONE = 1.d0, THR = 3.d0,FOUR=4.d0
  real*8,parameter::TINYRR=1.d-14
!sanity check: M/Q = constant
  M = mass(1)
  Q = Qchar(1)
  u1 = M/Q
  u2 = M/Q
  do bhi=2,N
     M = mass(bhi)
     Q = Qchar(bhi)
     u1 = min(u1,M/Q)
     u2 = max(u2,M/Q)
  enddo
  if(u2-u1.gt.TINYRR)then
      write(*,*)"error in initial_punctureem.f90: get_initial_nbhsem; we need constant Mi/Qi, but"
      write(*,*)"Mass = ",mass
      write(*,*)"Qchar = ",Qchar
      stop
   endif

  do k = 1,ext(3)
   do j = 1,ext(2)
    do i = 1,ext(1)
! black hole 1
      M = mass(1)
      Q = Qchar(1)
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
      if(rr.lt.TINYRR) rr=TINYRR
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
     u3 = ell/2.d1*(1.d0+ell+ell**2-4.d0*ell**3+2.d0*ell**4)
     u4 = ell**2/1.d1*(1.d1-2.5d1*ell+2.1d1*ell**2-6.d0*ell**3)

     tmp = (Py*Sz-Pz*Sy)*nx + (Pz*Sx-Px*Sz)*ny + (Px*Sy-Py*Sx)*nz

     u = PP**2/M**2*(u1 + u2*(3.d0*mup**2-ONE)) + &
         6.d0*u3/M**4*SS**2*(1.d0+mus**2) + u4/M**3*tmp

     psi(i,j,k) = ONE + u + HLF*M/rr
     phi(i,j,k) = Q/rr

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

     Ex(i,j,k) = Q*nx/rr/rr                
     Ey(i,j,k) = Q*ny/rr/rr                
     Ez(i,j,k) = Q*nz/rr/rr                
! black hole 2 and 3, ...
     do bhi=2,N
      M = Mass(bhi)
      Q = Qchar(bhi)
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
      if(rr.lt.TINYRR) rr=TINYRR
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
     u3 = ell/2.d1*(1.d0+ell+ell**2-4.d0*ell**3+2.d0*ell**4)
     u4 = ell**2/1.d1*(1.d1-2.5d1*ell+2.1d1*ell**2-6.d0*ell**3)

     tmp = (Py*Sz-Pz*Sy)*nx + (Pz*Sx-Px*Sz)*ny + (Px*Sy-Py*Sx)*nz

     u = PP**2/M**2*(u1 + u2*(3.d0*mup**2-ONE)) + &
         6.d0*u3/M**4*SS**2*(1.d0+mus**2) + u4/M**3*tmp

     psi(i,j,k) = psi(i,j,k) + u + HLF*M/rr
     phi(i,j,k) = phi(i,j,k) + Q/rr

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

     Ex(i,j,k) = Ex(i,j,k) + Q*nx/rr/rr                
     Ey(i,j,k) = Ey(i,j,k) + Q*ny/rr/rr                
     Ez(i,j,k) = Ez(i,j,k) + Q*nz/rr/rr  
      enddo
    enddo
   enddo
  enddo

  psi = dsqrt(psi**2 - phi*phi/FOUR)
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

  Ex = Ex   / psi
  Ey = Ey   / psi
  Ez = Ez   / psi

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
 
  Bx = ZEO
  By = ZEO
  Bz = ZEO

  Kpsi = ZEO
  Kphi = ZEO

  return

  end subroutine get_initial_nbhsem_ss
!-----------------------------------------------------------------------------------
!
!Set up approximate puncture initial data for n charged black holes
!aided with Ansorg's solver
!-----------------------------------------------------------------------------------

  subroutine get_ansorg_nbhs_em(ext,X,Y,Z, &
                         chi,  trK, &
                         gxx,  gxy,  gxz,  gyy,  gyz,  gzz,&
                         Axx,  Axy,  Axz,  Ayy,  Ayz,  Azz,&
                         Gmx,  Gmy,  Gmz,  &
                         Lap,  Sfx,  Sfy,  Sfz,&
                         dtSfx,dtSfy,dtSfz,&
                         Ex,Ey,Ez,Bx,By,Bz,Kpsi,Kphi,&
                         Mass,Qchar,Porg,Pmom,Spin,N)

  implicit none

!------= input arguments
 
  integer,intent(in) :: N
  integer, dimension(3),  intent(in) :: ext
  real*8,  dimension(ext(1)), intent(in) :: X
  real*8,  dimension(ext(2)), intent(in) :: Y
  real*8,  dimension(ext(3)), intent(in) :: Z
  real*8,  dimension(ext(1),ext(2),ext(3)), intent(inout) :: chi
  real*8,  dimension(ext(1),ext(2),ext(3)), intent(out) :: Ex,Ey,Ez,Bx,By,Bz,Kpsi,Kphi
  real*8,  dimension(ext(1),ext(2),ext(3)), intent(out) :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8,  dimension(ext(1),ext(2),ext(3)), intent(out) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8,  dimension(ext(1),ext(2),ext(3)), intent(out) :: trK,Lap,Sfx,Sfy,Sfz
  real*8,  dimension(ext(1),ext(2),ext(3)), intent(out) :: Gmx,Gmy,Gmz
  real*8,  dimension(ext(1),ext(2),ext(3)), intent(out) :: dtSfx,dtSfy,dtSfz
  real*8,  dimension(N), intent(in) :: Mass,Qchar
  real*8,  dimension(3*N), intent(in) :: Porg,Pmom,Spin

!------= local variables
  real*8,dimension(ext(1),ext(2),ext(3))::psi,phi
  integer :: i,j,k,bhi
  real*8 :: M,Q,Px,Py,Pz,Sx,Sy,Sz
  real*8 :: nx,ny,nz,rr,tmp
  real*8,  parameter :: HLF = 5.d-1, ZEO = 0.d0, ONE = 1.d0, THR = 3.d0,FOUR=4.d0
  real*8,parameter::TINYRR=1.d-14

  do k = 1,ext(3)
   do j = 1,ext(2)
    do i = 1,ext(1)
! black hole 1
      M = mass(1)
      Q = Qchar(1)
      nx = x(i) - Porg(1)
      ny = y(j) - Porg(2)
      nz = z(k) - Porg(3)
      Px = Pmom(1)
      Py = Pmom(2)
      Pz = Pmom(3)
      Sx = Spin(1)
      Sy = Spin(2)
      Sz = Spin(3)

      rr = dsqrt(nx*nx+ny*ny+nz*nz)
      if(rr.lt.TINYRR) rr=(X(2)-X(1))/2.d0
      nx = nx / rr
      ny = ny / rr
      nz = nz / rr

     psi(i,j,k) = ONE + chi(i,j,k) + HLF*M/rr
     phi(i,j,k) = Q/rr

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

     Ex(i,j,k) = Q*nx/rr/rr                
     Ey(i,j,k) = Q*ny/rr/rr                
     Ez(i,j,k) = Q*nz/rr/rr                
! black hole 2 and 3, ...
     do bhi=2,N
      M = Mass(bhi)
      Q = Qchar(bhi)
      nx = x(i) - Porg(3*(bhi-1)+1)
      ny = y(j) - Porg(3*(bhi-1)+2)
      nz = z(k) - Porg(3*(bhi-1)+3)
      Px = Pmom(3*(bhi-1)+1)
      Py = Pmom(3*(bhi-1)+2)
      Pz = Pmom(3*(bhi-1)+3)
      Sx = Spin(3*(bhi-1)+1)
      Sy = Spin(3*(bhi-1)+2)
      Sz = Spin(3*(bhi-1)+3)

      rr = dsqrt(nx*nx+ny*ny+nz*nz)
      if(rr.lt.TINYRR) rr=(X(2)-X(1))/2.d0
      nx = nx / rr
      ny = ny / rr
      nz = nz / rr

     psi(i,j,k) = psi(i,j,k) + HLF*M/rr
     phi(i,j,k) = phi(i,j,k) + Q/rr

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

     Ex(i,j,k) = Ex(i,j,k) + Q*nx/rr/rr                
     Ey(i,j,k) = Ey(i,j,k) + Q*ny/rr/rr                
     Ez(i,j,k) = Ez(i,j,k) + Q*nz/rr/rr  
      enddo
    enddo
   enddo
  enddo

  psi = dsqrt(psi**2 - phi*phi/FOUR)
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

  Ex = Ex   / psi
  Ey = Ey   / psi
  Ez = Ez   / psi

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
 
  Bx = ZEO
  By = ZEO
  Bz = ZEO

  Kpsi = ZEO
  Kphi = ZEO

  return

  end subroutine get_ansorg_nbhs_em
!-----------------------------------------------------------------------------------
!
!Set up approximate puncture initial data for n charged black holes
!aided with Ansorg's solver
! for shell
!-----------------------------------------------------------------------------------

  subroutine get_ansorg_nbhs_ss_em(ext,X,Y,Z, &
                         chi,  trK, &
                         gxx,  gxy,  gxz,  gyy,  gyz,  gzz,&
                         Axx,  Axy,  Axz,  Ayy,  Ayz,  Azz,&
                         Gmx,  Gmy,  Gmz,  &
                         Lap,  Sfx,  Sfy,  Sfz,&
                         dtSfx,dtSfy,dtSfz,&
                         Ex,Ey,Ez,Bx,By,Bz,Kpsi,Kphi,&
                         Mass,Qchar,Porg,Pmom,Spin,N)

  implicit none

!------= input arguments
 
  integer,intent(in) :: N
  integer, dimension(3),  intent(in) :: ext
  real*8,  dimension(ext(1),ext(2),ext(3)), intent(in) :: X,Y,Z
  real*8,  dimension(ext(1),ext(2),ext(3)), intent(inout) :: chi
  real*8,  dimension(ext(1),ext(2),ext(3)), intent(out) :: Ex,Ey,Ez,Bx,By,Bz,Kpsi,Kphi
  real*8,  dimension(ext(1),ext(2),ext(3)), intent(out) :: gxx,gxy,gxz,gyy,gyz,gzz
  real*8,  dimension(ext(1),ext(2),ext(3)), intent(out) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8,  dimension(ext(1),ext(2),ext(3)), intent(out) :: trK,Lap,Sfx,Sfy,Sfz
  real*8,  dimension(ext(1),ext(2),ext(3)), intent(out) :: Gmx,Gmy,Gmz
  real*8,  dimension(ext(1),ext(2),ext(3)), intent(out) :: dtSfx,dtSfy,dtSfz
  real*8,  dimension(N), intent(in) :: Mass,Qchar
  real*8,  dimension(3*N), intent(in) :: Porg,Pmom,Spin

!------= local variables
  real*8,dimension(ext(1),ext(2),ext(3))::psi,phi
  integer :: i,j,k,bhi
  real*8 :: M,Q,Px,Py,Pz,Sx,Sy,Sz
  real*8 :: nx,ny,nz,rr,tmp
  real*8,  parameter :: HLF = 5.d-1, ZEO = 0.d0, ONE = 1.d0, THR = 3.d0,FOUR=4.d0
  real*8,parameter::TINYRR=1.d-14

  do k = 1,ext(3)
   do j = 1,ext(2)
    do i = 1,ext(1)
! black hole 1
      M = mass(1)
      Q = Qchar(1)
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
      if(rr.lt.TINYRR) rr=TINYRR
      nx = nx / rr
      ny = ny / rr
      nz = nz / rr

     psi(i,j,k) = ONE + chi(i,j,k) + HLF*M/rr
     phi(i,j,k) = Q/rr

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

     Ex(i,j,k) = Q*nx/rr/rr                
     Ey(i,j,k) = Q*ny/rr/rr                
     Ez(i,j,k) = Q*nz/rr/rr                
! black hole 2 and 3, ...
     do bhi=2,N
      M = Mass(bhi)
      Q = Qchar(bhi)
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
      if(rr.lt.TINYRR) rr=TINYRR
      nx = nx / rr
      ny = ny / rr
      nz = nz / rr

     psi(i,j,k) = psi(i,j,k) + HLF*M/rr
     phi(i,j,k) = phi(i,j,k) + Q/rr

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

     Ex(i,j,k) = Ex(i,j,k) + Q*nx/rr/rr                
     Ey(i,j,k) = Ey(i,j,k) + Q*ny/rr/rr                
     Ez(i,j,k) = Ez(i,j,k) + Q*nz/rr/rr  
      enddo
    enddo
   enddo
  enddo

  psi = dsqrt(psi**2 - phi*phi/FOUR)
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

  Ex = Ex   / psi
  Ey = Ey   / psi
  Ez = Ez   / psi

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
 
  Bx = ZEO
  By = ZEO
  Bz = ZEO

  Kpsi = ZEO
  Kphi = ZEO

  return

  end subroutine get_ansorg_nbhs_ss_em
