
!including advection term in this routine
  function compute_rhs_empart(ext, X, Y, Z,                                          &
               chi    ,   dxx    ,   dxy    ,   dxz    ,   dyy    ,   dyz    ,   dzz,&
               Lap    ,  betax   ,  betay   ,  betaz   , trK,                        &
               Ex, Ey, Ez, Bx, By, Bz, Kpsi, Kphi,Jx,Jy,Jz,qchar,                    &
               Ex_rhs, Ey_rhs, Ez_rhs, Bx_rhs, By_rhs, Bz_rhs, Kpsi_rhs, Kphi_rhs,   &
               rho,Sx,Sy,Sz,Sxx,Sxy,Sxz,Syy,Syz,Szz,                                 &
               Symmetry,Lev,eps)  result(gont)
  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ext(1:3), Symmetry,Lev
  real*8, intent(in ):: X(1:ext(1)),Y(1:ext(2)),Z(1:ext(3))
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in ) :: chi,Jx,Jy,Jz,qchar
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in ) :: dxx,dxy,dxz,dyy,dyz,dzz
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in ) :: Lap, betax, betay, betaz, trK
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in ) :: Ex,Ey,Ez,Bx,By,Bz,Kpsi,Kphi
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out) :: Ex_rhs, Ey_rhs, Ez_rhs 
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out) :: Bx_rhs, By_rhs, Bz_rhs, Kpsi_rhs, Kphi_rhs
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out) :: rho,Sx,Sy,Sz
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out) :: Sxx,Sxy,Sxz,Syy,Syz,Szz
  real*8,intent(in) :: eps
!  gont = 0: success; gont = 1: something wrong
  integer::gont

!~~~~~~> Other variables:

  real*8, dimension(ext(1),ext(2),ext(3)) :: gxx,gyy,gzz,gxy,gxz,gyz
  real*8, dimension(ext(1),ext(2),ext(3)) :: chix,chiy,chiz,chi3o2
  real*8, dimension(ext(1),ext(2),ext(3)) :: gxxx,gxyx,gxzx,gyyx,gyzx,gzzx
  real*8, dimension(ext(1),ext(2),ext(3)) :: gxxy,gxyy,gxzy,gyyy,gyzy,gzzy
  real*8, dimension(ext(1),ext(2),ext(3)) :: gxxz,gxyz,gxzz,gyyz,gyzz,gzzz
  real*8, dimension(ext(1),ext(2),ext(3)) :: Lapx,Lapy,Lapz
  real*8, dimension(ext(1),ext(2),ext(3)) :: betaxx,betaxy,betaxz
  real*8, dimension(ext(1),ext(2),ext(3)) :: betayx,betayy,betayz
  real*8, dimension(ext(1),ext(2),ext(3)) :: betazx,betazy,betazz
  real*8, dimension(ext(1),ext(2),ext(3)) :: alpn1,chin1
  real*8, dimension(ext(1),ext(2),ext(3)) :: gupxx,gupxy,gupxz
  real*8, dimension(ext(1),ext(2),ext(3)) :: gupyy,gupyz,gupzz
  real*8, dimension(ext(1),ext(2),ext(3)) :: Exx,Exy,Exz,Eyx,Eyy,Eyz,Ezx,Ezy,Ezz
  real*8, dimension(ext(1),ext(2),ext(3)) :: Bxx,Bxy,Bxz,Byx,Byy,Byz,Bzx,Bzy,Bzz
  real*8, dimension(ext(1),ext(2),ext(3)) :: Kpsix,Kpsiy,Kpsiz
  real*8, dimension(ext(1),ext(2),ext(3)) :: Kphix,Kphiy,Kphiz
  real*8, dimension(ext(1),ext(2),ext(3)) :: lEx,lEy,lEz,lBx,lBy,lBz

  real*8,dimension(3) ::SSS,AAS,ASA,SAA,ASS,SAS,SSA
  real*8            :: dX, dY, dZ, PI
  real*8, parameter :: ONE = 1.D0, TWO = 2.D0, FOUR = 4.D0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  real*8, parameter :: F3o2=1.5d0,EIT=8.d0
  real*8,parameter  :: kappa = 1.d0
  integer :: i, j, k
!!! sanity check
  dX = sum(Ex)+sum(Ey)+sum(Ez)+sum(Bx)+sum(By)+sum(Bz)+sum(Kpsi)+sum(Kphi)
  if(dX.ne.dX) then
     if(sum(Ex).ne.sum(Ex))write(*,*)"empart.f90: find NaN in Ex"
     if(sum(Ey).ne.sum(Ey))write(*,*)"empart.f90: find NaN in Ey"
     if(sum(Ez).ne.sum(Ez))write(*,*)"empart.f90: find NaN in Ez"
     if(sum(Bx).ne.sum(Bx))write(*,*)"empart.f90: find NaN in Bx"
     if(sum(By).ne.sum(By))write(*,*)"empart.f90: find NaN in By"
     if(sum(Bz).ne.sum(Bz))write(*,*)"empart.f90: find NaN in Bz"
     if(sum(Kpsi).ne.sum(Kpsi))write(*,*)"empart.f90: find NaN in Kpsi"
     if(sum(Kphi).ne.sum(Kphi))write(*,*)"empart.f90: find NaN in Kphi"
     gont = 1
     return
  endif

  PI = dacos(-ONE)

  dX = X(2) - X(1)
  dY = Y(2) - Y(1)
  dZ = Z(2) - Z(1)

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ext(3)
  do j=1,ext(2)
  do i=1,ext(1)
  alpn1(i,j,k) = Lap(i,j,k) + ONE
  chin1(i,j,k) = chi(i,j,k) + ONE
  chi3o2(i,j,k)  = dsqrt(chin1(i,j,k))**3
  gxx(i,j,k) = dxx(i,j,k) + ONE
  gyy(i,j,k) = dyy(i,j,k) + ONE
  gzz(i,j,k) = dzz(i,j,k) + ONE
  gxy(i,j,k) = dxy(i,j,k)
  gxz(i,j,k) = dxz(i,j,k)
  gyz(i,j,k) = dyz(i,j,k)
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

  call fderivs(ext,Lap,Lapx,Lapy,Lapz,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)
  call fderivs(ext,betax,betaxx,betaxy,betaxz,X,Y,Z,ANTI, SYM, SYM,Symmetry,Lev)
  call fderivs(ext,betay,betayx,betayy,betayz,X,Y,Z, SYM,ANTI, SYM,Symmetry,Lev)
  call fderivs(ext,betaz,betazx,betazy,betazz,X,Y,Z, SYM, SYM,ANTI,Symmetry,Lev)
 
  call fderivs(ext,chi,chix,chiy,chiz,X,Y,Z,SYM,SYM,SYM,symmetry,Lev)

  call fderivs(ext,dxx,gxxx,gxxy,gxxz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)
  call fderivs(ext,gxy,gxyx,gxyy,gxyz,X,Y,Z,ANTI,ANTI,SYM ,Symmetry,Lev)
  call fderivs(ext,gxz,gxzx,gxzy,gxzz,X,Y,Z,ANTI,SYM ,ANTI,Symmetry,Lev)
  call fderivs(ext,dyy,gyyx,gyyy,gyyz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)
  call fderivs(ext,gyz,gyzx,gyzy,gyzz,X,Y,Z,SYM ,ANTI,ANTI,Symmetry,Lev)
  call fderivs(ext,dzz,gzzx,gzzy,gzzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)

  call fderivs(ext,Kpsi,Kpsix,Kpsiy,Kpsiz,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)
  call fderivs(ext,Kphi,Kphix,Kphiy,Kphiz,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)

  call fderivs(ext,Ex,Exx,Exy,Exz,X,Y,Z,ANTI,SYM,SYM ,Symmetry,Lev)
  call fderivs(ext,Ey,Eyx,Eyy,Eyz,X,Y,Z,SYM,ANTI,SYM ,Symmetry,Lev)
  call fderivs(ext,Ez,Ezx,Ezy,Ezz,X,Y,Z,SYM,SYM,ANTI ,Symmetry,Lev)

  call fderivs(ext,Bx,Bxx,Bxy,Bxz,X,Y,Z,SYM,ANTI,ANTI ,Symmetry,Lev)
  call fderivs(ext,By,Byx,Byy,Byz,X,Y,Z,ANTI,SYM,ANTI ,Symmetry,Lev)
  call fderivs(ext,Bz,Bzx,Bzy,Bzz,X,Y,Z,ANTI,ANTI,SYM ,Symmetry,Lev)

! physical gij
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ext(3)
  do j=1,ext(2)
  do i=1,ext(1)
  gxx(i,j,k) = gxx(i,j,k)/chin1(i,j,k)
  gxy(i,j,k) = gxy(i,j,k)/chin1(i,j,k)
  gxz(i,j,k) = gxz(i,j,k)/chin1(i,j,k)
  gyy(i,j,k) = gyy(i,j,k)/chin1(i,j,k)
  gyz(i,j,k) = gyz(i,j,k)/chin1(i,j,k)
  gzz(i,j,k) = gzz(i,j,k)/chin1(i,j,k)
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO
!physical gij,k  
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ext(3)
  do j=1,ext(2)
  do i=1,ext(1)
  gxxx(i,j,k) = (gxxx(i,j,k)-gxx(i,j,k)*chix(i,j,k))/chin1(i,j,k)
  gxxy(i,j,k) = (gxxy(i,j,k)-gxx(i,j,k)*chiy(i,j,k))/chin1(i,j,k)
  gxxz(i,j,k) = (gxxz(i,j,k)-gxx(i,j,k)*chiz(i,j,k))/chin1(i,j,k)
  gxyx(i,j,k) = (gxyx(i,j,k)-gxy(i,j,k)*chix(i,j,k))/chin1(i,j,k)
  gxyy(i,j,k) = (gxyy(i,j,k)-gxy(i,j,k)*chiy(i,j,k))/chin1(i,j,k)
  gxyz(i,j,k) = (gxyz(i,j,k)-gxy(i,j,k)*chiz(i,j,k))/chin1(i,j,k)
  gxzx(i,j,k) = (gxzx(i,j,k)-gxz(i,j,k)*chix(i,j,k))/chin1(i,j,k)
  gxzy(i,j,k) = (gxzy(i,j,k)-gxz(i,j,k)*chiy(i,j,k))/chin1(i,j,k)
  gxzz(i,j,k) = (gxzz(i,j,k)-gxz(i,j,k)*chiz(i,j,k))/chin1(i,j,k)
  gyyx(i,j,k) = (gyyx(i,j,k)-gyy(i,j,k)*chix(i,j,k))/chin1(i,j,k)
  gyyy(i,j,k) = (gyyy(i,j,k)-gyy(i,j,k)*chiy(i,j,k))/chin1(i,j,k)
  gyyz(i,j,k) = (gyyz(i,j,k)-gyy(i,j,k)*chiz(i,j,k))/chin1(i,j,k)
  gyzx(i,j,k) = (gyzx(i,j,k)-gyz(i,j,k)*chix(i,j,k))/chin1(i,j,k)
  gyzy(i,j,k) = (gyzy(i,j,k)-gyz(i,j,k)*chiy(i,j,k))/chin1(i,j,k)
  gyzz(i,j,k) = (gyzz(i,j,k)-gyz(i,j,k)*chiz(i,j,k))/chin1(i,j,k)
  gzzx(i,j,k) = (gzzx(i,j,k)-gzz(i,j,k)*chix(i,j,k))/chin1(i,j,k)
  gzzy(i,j,k) = (gzzy(i,j,k)-gzz(i,j,k)*chiy(i,j,k))/chin1(i,j,k)
  gzzz(i,j,k) = (gzzz(i,j,k)-gzz(i,j,k)*chiz(i,j,k))/chin1(i,j,k)
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

! physical inverse metric
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ext(3)
  do j=1,ext(2)
  do i=1,ext(1)
  gupzz(i,j,k) =  gxx(i,j,k) * gyy(i,j,k) * gzz(i,j,k) + gxy(i,j,k) * gyz(i,j,k) * gxz(i,j,k) + gxz(i,j,k) * gxy(i,j,k) * gyz(i,j,k) - &
           gxz(i,j,k) * gyy(i,j,k) * gxz(i,j,k) - gxy(i,j,k) * gxy(i,j,k) * gzz(i,j,k) - gxx(i,j,k) * gyz(i,j,k) * gyz(i,j,k)
  gupxx(i,j,k) =   ( gyy(i,j,k) * gzz(i,j,k) - gyz(i,j,k) * gyz(i,j,k) ) / gupzz(i,j,k)
  gupxy(i,j,k) = - ( gxy(i,j,k) * gzz(i,j,k) - gyz(i,j,k) * gxz(i,j,k) ) / gupzz(i,j,k)
  gupxz(i,j,k) =   ( gxy(i,j,k) * gyz(i,j,k) - gyy(i,j,k) * gxz(i,j,k) ) / gupzz(i,j,k)
  gupyy(i,j,k) =   ( gxx(i,j,k) * gzz(i,j,k) - gxz(i,j,k) * gxz(i,j,k) ) / gupzz(i,j,k)
  gupyz(i,j,k) = - ( gxx(i,j,k) * gyz(i,j,k) - gxy(i,j,k) * gxz(i,j,k) ) / gupzz(i,j,k)
  gupzz(i,j,k) =   ( gxx(i,j,k) * gyy(i,j,k) - gxy(i,j,k) * gxy(i,j,k) ) / gupzz(i,j,k)
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ext(3)
  do j=1,ext(2)
  do i=1,ext(1)
  Ex_rhs(i,j,k) = alpn1(i,j,k)*trK(i,j,k)*Ex(i,j,k)-(Ex(i,j,k)*betaxx(i,j,k)+Ey(i,j,k)*betaxy(i,j,k)+Ez(i,j,k)*betaxz(i,j,k))                  &
          -FOUR*PI*alpn1(i,j,k)*Jx(i,j,k)-alpn1(i,j,k)*(gupxx(i,j,k)*Kpsix(i,j,k)+gupxy(i,j,k)*Kpsiy(i,j,k)+gupxz(i,j,k)*Kpsiz(i,j,k))  &
          +chi3o2(i,j,k)*(                                                      &
          ((gxz(i,j,k)*Bx(i,j,k)+gyz(i,j,k)*By(i,j,k)+gzz(i,j,k)*Bz(i,j,k))*Lapy(i,j,k)+alpn1(i,j,k)*(gxz(i,j,k)*Bxy(i,j,k)+gyz(i,j,k)*Byy(i,j,k)+gzz(i,j,k)*Bzy(i,j,k))+alpn1(i,j,k)*(Bx(i,j,k)*gxzy(i,j,k)+By(i,j,k)*gyzy(i,j,k)+Bz(i,j,k)*gzzy(i,j,k)))-&
          ((gxy(i,j,k)*Bx(i,j,k)+gyy(i,j,k)*By(i,j,k)+gyz(i,j,k)*Bz(i,j,k))*Lapz(i,j,k)+alpn1(i,j,k)*(gxy(i,j,k)*Bxz(i,j,k)+gyy(i,j,k)*Byz(i,j,k)+gyz(i,j,k)*Bzz(i,j,k))+alpn1(i,j,k)*(Bx(i,j,k)*gxyz(i,j,k)+By(i,j,k)*gyyz(i,j,k)+Bz(i,j,k)*gyzz(i,j,k))))
  Ey_rhs(i,j,k) = alpn1(i,j,k)*trK(i,j,k)*Ey(i,j,k)-(Ex(i,j,k)*betayx(i,j,k)+Ey(i,j,k)*betayy(i,j,k)+Ez(i,j,k)*betayz(i,j,k))                  &
          -FOUR*PI*alpn1(i,j,k)*Jy(i,j,k)-alpn1(i,j,k)*(gupxy(i,j,k)*Kpsix(i,j,k)+gupyy(i,j,k)*Kpsiy(i,j,k)+gupyz(i,j,k)*Kpsiz(i,j,k))  &
          +chi3o2(i,j,k)*(                                                      &
          ((gxx(i,j,k)*Bx(i,j,k)+gxy(i,j,k)*By(i,j,k)+gxz(i,j,k)*Bz(i,j,k))*Lapz(i,j,k)+alpn1(i,j,k)*(gxx(i,j,k)*Bxz(i,j,k)+gxy(i,j,k)*Byz(i,j,k)+gxz(i,j,k)*Bzz(i,j,k))+alpn1(i,j,k)*(Bx(i,j,k)*gxxz(i,j,k)+By(i,j,k)*gxyz(i,j,k)+Bz(i,j,k)*gxzz(i,j,k)))-&
          ((gxz(i,j,k)*Bx(i,j,k)+gyz(i,j,k)*By(i,j,k)+gzz(i,j,k)*Bz(i,j,k))*Lapx(i,j,k)+alpn1(i,j,k)*(gxz(i,j,k)*Bxx(i,j,k)+gyz(i,j,k)*Byx(i,j,k)+gzz(i,j,k)*Bzx(i,j,k))+alpn1(i,j,k)*(Bx(i,j,k)*gxzx(i,j,k)+By(i,j,k)*gyzx(i,j,k)+Bz(i,j,k)*gzzx(i,j,k))))
  Ez_rhs(i,j,k) = alpn1(i,j,k)*trK(i,j,k)*Ez(i,j,k)-(Ex(i,j,k)*betazx(i,j,k)+Ey(i,j,k)*betazy(i,j,k)+Ez(i,j,k)*betazz(i,j,k))                  &
          -FOUR*PI*alpn1(i,j,k)*Jz(i,j,k)-alpn1(i,j,k)*(gupxz(i,j,k)*Kpsix(i,j,k)+gupyz(i,j,k)*Kpsiy(i,j,k)+gupzz(i,j,k)*Kpsiz(i,j,k))  &
          +chi3o2(i,j,k)*(                                                      &
          ((gxy(i,j,k)*Bx(i,j,k)+gyy(i,j,k)*By(i,j,k)+gyz(i,j,k)*Bz(i,j,k))*Lapx(i,j,k)+alpn1(i,j,k)*(gxy(i,j,k)*Bxx(i,j,k)+gyy(i,j,k)*Byx(i,j,k)+gyz(i,j,k)*Bzx(i,j,k))+alpn1(i,j,k)*(Bx(i,j,k)*gxyx(i,j,k)+By(i,j,k)*gyyx(i,j,k)+Bz(i,j,k)*gyzx(i,j,k)))-&
          ((gxx(i,j,k)*Bx(i,j,k)+gxy(i,j,k)*By(i,j,k)+gxz(i,j,k)*Bz(i,j,k))*Lapy(i,j,k)+alpn1(i,j,k)*(gxx(i,j,k)*Bxy(i,j,k)+gxy(i,j,k)*Byy(i,j,k)+gxz(i,j,k)*Bzy(i,j,k))+alpn1(i,j,k)*(Bx(i,j,k)*gxxy(i,j,k)+By(i,j,k)*gxyy(i,j,k)+Bz(i,j,k)*gxzy(i,j,k))))

  Bx_rhs(i,j,k) = alpn1(i,j,k)*trK(i,j,k)*Bx(i,j,k)-(Bx(i,j,k)*betaxx(i,j,k)+By(i,j,k)*betaxy(i,j,k)+Bz(i,j,k)*betaxz(i,j,k))                  &
                           -alpn1(i,j,k)*(gupxx(i,j,k)*Kphix(i,j,k)+gupxy(i,j,k)*Kphiy(i,j,k)+gupxz(i,j,k)*Kphiz(i,j,k))  &
          -chi3o2(i,j,k)*(                                                      &
          ((gxz(i,j,k)*Ex(i,j,k)+gyz(i,j,k)*Ey(i,j,k)+gzz(i,j,k)*Ez(i,j,k))*Lapy(i,j,k)+alpn1(i,j,k)*(gxz(i,j,k)*Exy(i,j,k)+gyz(i,j,k)*Eyy(i,j,k)+gzz(i,j,k)*Ezy(i,j,k))+alpn1(i,j,k)*(Ex(i,j,k)*gxzy(i,j,k)+Ey(i,j,k)*gyzy(i,j,k)+Ez(i,j,k)*gzzy(i,j,k)))-&
          ((gxy(i,j,k)*Ex(i,j,k)+gyy(i,j,k)*Ey(i,j,k)+gyz(i,j,k)*Ez(i,j,k))*Lapz(i,j,k)+alpn1(i,j,k)*(gxy(i,j,k)*Exz(i,j,k)+gyy(i,j,k)*Eyz(i,j,k)+gyz(i,j,k)*Ezz(i,j,k))+alpn1(i,j,k)*(Ex(i,j,k)*gxyz(i,j,k)+Ey(i,j,k)*gyyz(i,j,k)+Ez(i,j,k)*gyzz(i,j,k))))
  By_rhs(i,j,k) = alpn1(i,j,k)*trK(i,j,k)*By(i,j,k)-(Bx(i,j,k)*betayx(i,j,k)+By(i,j,k)*betayy(i,j,k)+Bz(i,j,k)*betayz(i,j,k))                  &
                           -alpn1(i,j,k)*(gupxy(i,j,k)*Kphix(i,j,k)+gupyy(i,j,k)*Kphiy(i,j,k)+gupyz(i,j,k)*Kphiz(i,j,k))  &
          -chi3o2(i,j,k)*(                                                      &
          ((gxx(i,j,k)*Ex(i,j,k)+gxy(i,j,k)*Ey(i,j,k)+gxz(i,j,k)*Ez(i,j,k))*Lapz(i,j,k)+alpn1(i,j,k)*(gxx(i,j,k)*Exz(i,j,k)+gxy(i,j,k)*Eyz(i,j,k)+gxz(i,j,k)*Ezz(i,j,k))+alpn1(i,j,k)*(Ex(i,j,k)*gxxz(i,j,k)+Ey(i,j,k)*gxyz(i,j,k)+Ez(i,j,k)*gxzz(i,j,k)))-&
          ((gxz(i,j,k)*Ex(i,j,k)+gyz(i,j,k)*Ey(i,j,k)+gzz(i,j,k)*Ez(i,j,k))*Lapx(i,j,k)+alpn1(i,j,k)*(gxz(i,j,k)*Exx(i,j,k)+gyz(i,j,k)*Eyx(i,j,k)+gzz(i,j,k)*Ezx(i,j,k))+alpn1(i,j,k)*(Ex(i,j,k)*gxzx(i,j,k)+Ey(i,j,k)*gyzx(i,j,k)+Ez(i,j,k)*gzzx(i,j,k))))
  Bz_rhs(i,j,k) = alpn1(i,j,k)*trK(i,j,k)*Bz(i,j,k)-(Bx(i,j,k)*betazx(i,j,k)+By(i,j,k)*betazy(i,j,k)+Bz(i,j,k)*betazz(i,j,k))                  &
                           -alpn1(i,j,k)*(gupxz(i,j,k)*Kphix(i,j,k)+gupyz(i,j,k)*Kphiy(i,j,k)+gupzz(i,j,k)*Kphiz(i,j,k))  &
          -chi3o2(i,j,k)*(                                                      &
          ((gxy(i,j,k)*Ex(i,j,k)+gyy(i,j,k)*Ey(i,j,k)+gyz(i,j,k)*Ez(i,j,k))*Lapx(i,j,k)+alpn1(i,j,k)*(gxy(i,j,k)*Exx(i,j,k)+gyy(i,j,k)*Eyx(i,j,k)+gyz(i,j,k)*Ezx(i,j,k))+alpn1(i,j,k)*(Ex(i,j,k)*gxyx(i,j,k)+Ey(i,j,k)*gyyx(i,j,k)+Ez(i,j,k)*gyzx(i,j,k)))-&
          ((gxx(i,j,k)*Ex(i,j,k)+gxy(i,j,k)*Ey(i,j,k)+gxz(i,j,k)*Ez(i,j,k))*Lapy(i,j,k)+alpn1(i,j,k)*(gxx(i,j,k)*Exy(i,j,k)+gxy(i,j,k)*Eyy(i,j,k)+gxz(i,j,k)*Ezy(i,j,k))+alpn1(i,j,k)*(Ex(i,j,k)*gxxy(i,j,k)+Ey(i,j,k)*gxyy(i,j,k)+Ez(i,j,k)*gxzy(i,j,k))))

  Kpsi_rhs(i,j,k) = FOUR*PI*alpn1(i,j,k)*qchar(i,j,k)-alpn1(i,j,k)*kappa*Kpsi(i,j,k) - &
            alpn1(i,j,k)*(Exx(i,j,k)+Eyy(i,j,k)+Ezz(i,j,k)-F3o2/chin1(i,j,k)*(chix(i,j,k)*Ex(i,j,k)+chiy(i,j,k)*Ey(i,j,k)+chiz(i,j,k)*Ez(i,j,k))) 
  Kphi_rhs(i,j,k) =                    -alpn1(i,j,k)*kappa*Kphi(i,j,k) - &
            alpn1(i,j,k)*(Bxx(i,j,k)+Byy(i,j,k)+Bzz(i,j,k)-F3o2/chin1(i,j,k)*(chix(i,j,k)*Bx(i,j,k)+chiy(i,j,k)*By(i,j,k)+chiz(i,j,k)*Bz(i,j,k))) 
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

  SSS(1)=SYM
  SSS(2)=SYM
  SSS(3)=SYM

  AAS(1)=ANTI
  AAS(2)=ANTI
  AAS(3)=SYM

  ASA(1)=ANTI
  ASA(2)=SYM
  ASA(3)=ANTI

  SAA(1)=SYM
  SAA(2)=ANTI
  SAA(3)=ANTI

  ASS(1)=ANTI
  ASS(2)=SYM
  ASS(3)=SYM

  SAS(1)=SYM
  SAS(2)=ANTI
  SAS(3)=SYM

  SSA(1)=SYM
  SSA(2)=SYM
  SSA(3)=ANTI

!!!!!!!!!advection term part  
  call lopsided(ext,X,Y,Z,KPsi,KPsi_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ext,X,Y,Z,KPhi,KPhi_rhs,betax,betay,betaz,Symmetry,SSS)

  call lopsided(ext,X,Y,Z,Ex,Ex_rhs,betax,betay,betaz,Symmetry,ASS)
  call lopsided(ext,X,Y,Z,Ey,Ey_rhs,betax,betay,betaz,Symmetry,SAS)
  call lopsided(ext,X,Y,Z,Ez,Ez_rhs,betax,betay,betaz,Symmetry,SSA)

  call lopsided(ext,X,Y,Z,Bx,Bx_rhs,betax,betay,betaz,Symmetry,SAA)
  call lopsided(ext,X,Y,Z,By,By_rhs,betax,betay,betaz,Symmetry,ASA)
  call lopsided(ext,X,Y,Z,Bz,Bz_rhs,betax,betay,betaz,Symmetry,AAS)

! numerical dissipation part
  if(eps>0)then 
! usual Kreiss-Oliger dissipation 

  call kodis(ext,X,Y,Z,Kpsi,Kpsi_rhs,SSS,Symmetry,eps)
  call kodis(ext,X,Y,Z,Kphi,Kphi_rhs,SSS,Symmetry,eps)
  call kodis(ext,X,Y,Z,Ex,Ex_rhs,ASS,Symmetry,eps)
  call kodis(ext,X,Y,Z,Ey,Ey_rhs,SAS,Symmetry,eps)
  call kodis(ext,X,Y,Z,Ez,Ez_rhs,SSA,Symmetry,eps)
  call kodis(ext,X,Y,Z,Bx,Bx_rhs,SAA,Symmetry,eps)
  call kodis(ext,X,Y,Z,By,By_rhs,ASA,Symmetry,eps)
  call kodis(ext,X,Y,Z,Bz,Bz_rhs,AAS,Symmetry,eps)

  endif
! stress-energy tensor
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ext(3)
  do j=1,ext(2)
  do i=1,ext(1)
  rho(i,j,k) = (gxx(i,j,k)*(Ex(i,j,k)*Ex(i,j,k)+Bx(i,j,k)*Bx(i,j,k))+gyy(i,j,k)*(Ey(i,j,k)*Ey(i,j,k)+By(i,j,k)*By(i,j,k))+gzz(i,j,k)*(Ez(i,j,k)*Ez(i,j,k)+Bz(i,j,k)*Bz(i,j,k)) + &
      +TWO*(gxy(i,j,k)*(Ex(i,j,k)*Ey(i,j,k)+Bx(i,j,k)*By(i,j,k))+gxz(i,j,k)*(Ex(i,j,k)*Ez(i,j,k)+Bx(i,j,k)*Bz(i,j,k))+gyz(i,j,k)*(Ey(i,j,k)*Ez(i,j,k)+By(i,j,k)*Bz(i,j,k))))/EIT/PI
  Sx(i,j,k) = (Ey(i,j,k)*Bz(i,j,k)-Ez(i,j,k)*By(i,j,k))/FOUR/PI/chi3o2(i,j,k)
  Sy(i,j,k) = (Ez(i,j,k)*Bx(i,j,k)-Ex(i,j,k)*Bz(i,j,k))/FOUR/PI/chi3o2(i,j,k)
  Sz(i,j,k) = (Ex(i,j,k)*By(i,j,k)-Ey(i,j,k)*Bx(i,j,k))/FOUR/PI/chi3o2(i,j,k)
  lEx(i,j,k) = gxx(i,j,k)*Ex(i,j,k)+gxy(i,j,k)*Ey(i,j,k)+gxz(i,j,k)*Ez(i,j,k)
  lEy(i,j,k) = gxy(i,j,k)*Ex(i,j,k)+gyy(i,j,k)*Ey(i,j,k)+gyz(i,j,k)*Ez(i,j,k)
  lEz(i,j,k) = gxz(i,j,k)*Ex(i,j,k)+gyz(i,j,k)*Ey(i,j,k)+gzz(i,j,k)*Ez(i,j,k)
  lBx(i,j,k) = gxx(i,j,k)*Bx(i,j,k)+gxy(i,j,k)*By(i,j,k)+gxz(i,j,k)*Bz(i,j,k)
  lBy(i,j,k) = gxy(i,j,k)*Bx(i,j,k)+gyy(i,j,k)*By(i,j,k)+gyz(i,j,k)*Bz(i,j,k)
  lBz(i,j,k) = gxz(i,j,k)*Bx(i,j,k)+gyz(i,j,k)*By(i,j,k)+gzz(i,j,k)*Bz(i,j,k)
  Sxx(i,j,k) = rho(i,j,k)*gxx(i,j,k)-(lEx(i,j,k)*lEx(i,j,k)+lBx(i,j,k)*lBx(i,j,k))/FOUR/PI
  Sxy(i,j,k) = rho(i,j,k)*gxy(i,j,k)-(lEx(i,j,k)*lEy(i,j,k)+lBx(i,j,k)*lBy(i,j,k))/FOUR/PI
  Sxz(i,j,k) = rho(i,j,k)*gxz(i,j,k)-(lEx(i,j,k)*lEz(i,j,k)+lBx(i,j,k)*lBz(i,j,k))/FOUR/PI
  Syy(i,j,k) = rho(i,j,k)*gyy(i,j,k)-(lEy(i,j,k)*lEy(i,j,k)+lBy(i,j,k)*lBy(i,j,k))/FOUR/PI
  Syz(i,j,k) = rho(i,j,k)*gyz(i,j,k)-(lEy(i,j,k)*lEz(i,j,k)+lBy(i,j,k)*lBz(i,j,k))/FOUR/PI
  Szz(i,j,k) = rho(i,j,k)*gzz(i,j,k)-(lEz(i,j,k)*lEz(i,j,k)+lBz(i,j,k)*lBz(i,j,k))/FOUR/PI
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

  gont = 0

  return

  end function compute_rhs_empart
!including advection term in this routine
! for shell
  function compute_rhs_empart_ss(ext,crho,sigma,R,x,y,z,                       &
               drhodx, drhody, drhodz,                                         &
               dsigmadx,dsigmady,dsigmadz,                                     &
               dRdx,dRdy,dRdz,                                                 &
               drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                &
               dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,    &
               dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,                            &
               chi    ,   dxx    ,   dxy    ,   dxz    ,   dyy    ,   dyz    ,   dzz,&
               Lap    ,  betax   ,  betay   ,  betaz   , trK,                        &
               Ex, Ey, Ez, Bx, By, Bz, Kpsi, Kphi,Jx,Jy,Jz,qchar,                    &
               Ex_rhs, Ey_rhs, Ez_rhs, Bx_rhs, By_rhs, Bz_rhs, Kpsi_rhs, Kphi_rhs,   &
               rho,Sx,Sy,Sz,Sxx,Sxy,Sxz,Syy,Syz,Szz,                                 &
               Symmetry,Lev,eps,sst)  result(gont)
  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ext(1:3), Symmetry,Lev,sst
  double precision,intent(in),dimension(ext(1))::crho
  double precision,intent(in),dimension(ext(2))::sigma
  double precision,intent(in),dimension(ext(3))::R
  double precision,intent(in),dimension(ext(1),ext(2),ext(3))::x,y,z
  double precision,intent(in),dimension(ext(1),ext(2),ext(3))::drhodx, drhody, drhodz
  double precision,intent(in),dimension(ext(1),ext(2),ext(3))::dsigmadx,dsigmady,dsigmadz
  double precision,intent(in),dimension(ext(1),ext(2),ext(3))::dRdx,dRdy,dRdz
  double precision,intent(in),dimension(ext(1),ext(2),ext(3))::drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz
  double precision,intent(in),dimension(ext(1),ext(2),ext(3))::dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz
  double precision,intent(in),dimension(ext(1),ext(2),ext(3))::dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in ) :: chi,Jx,Jy,Jz,qchar
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in ) :: dxx,dxy,dxz,dyy,dyz,dzz
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in ) :: Lap, betax, betay, betaz, trK
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in ) :: Ex,Ey,Ez,Bx,By,Bz,Kpsi,Kphi
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out) :: Ex_rhs, Ey_rhs, Ez_rhs 
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out) :: Bx_rhs, By_rhs, Bz_rhs, Kpsi_rhs, Kphi_rhs
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out) :: rho,Sx,Sy,Sz
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out) :: Sxx,Sxy,Sxz,Syy,Syz,Szz
  real*8,intent(in) :: eps
!  gont = 0: success; gont = 1: something wrong
  integer::gont

!~~~~~~> Other variables:

  real*8, dimension(ext(1),ext(2),ext(3)) :: gxx,gyy,gzz,gxy,gxz,gyz
  real*8, dimension(ext(1),ext(2),ext(3)) :: chix,chiy,chiz,chi3o2
  real*8, dimension(ext(1),ext(2),ext(3)) :: gxxx,gxyx,gxzx,gyyx,gyzx,gzzx
  real*8, dimension(ext(1),ext(2),ext(3)) :: gxxy,gxyy,gxzy,gyyy,gyzy,gzzy
  real*8, dimension(ext(1),ext(2),ext(3)) :: gxxz,gxyz,gxzz,gyyz,gyzz,gzzz
  real*8, dimension(ext(1),ext(2),ext(3)) :: Lapx,Lapy,Lapz
  real*8, dimension(ext(1),ext(2),ext(3)) :: betaxx,betaxy,betaxz
  real*8, dimension(ext(1),ext(2),ext(3)) :: betayx,betayy,betayz
  real*8, dimension(ext(1),ext(2),ext(3)) :: betazx,betazy,betazz
  real*8, dimension(ext(1),ext(2),ext(3)) :: alpn1,chin1
  real*8, dimension(ext(1),ext(2),ext(3)) :: gupxx,gupxy,gupxz
  real*8, dimension(ext(1),ext(2),ext(3)) :: gupyy,gupyz,gupzz
  real*8, dimension(ext(1),ext(2),ext(3)) :: Exx,Exy,Exz,Eyx,Eyy,Eyz,Ezx,Ezy,Ezz
  real*8, dimension(ext(1),ext(2),ext(3)) :: Bxx,Bxy,Bxz,Byx,Byy,Byz,Bzx,Bzy,Bzz
  real*8, dimension(ext(1),ext(2),ext(3)) :: Kpsix,Kpsiy,Kpsiz
  real*8, dimension(ext(1),ext(2),ext(3)) :: Kphix,Kphiy,Kphiz
  real*8, dimension(ext(1),ext(2),ext(3)) :: lEx,lEy,lEz,lBx,lBy,lBz

  real*8,dimension(3) ::SSS,AAS,ASA,SAA,ASS,SAS,SSA
  real*8            :: dX, dY, dZ, PI
  real*8, parameter :: ONE = 1.D0, TWO = 2.D0, FOUR = 4.D0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  real*8, parameter :: F3o2=1.5d0,EIT=8.d0
  real*8,parameter  :: kappa = 1.d0
  integer :: i, j, k

!!! sanity check
  dX = sum(Ex)+sum(Ey)+sum(Ez)+sum(Bx)+sum(By)+sum(Bz)+sum(Kpsi)+sum(Kphi)
  if(dX.ne.dX) then
     if(sum(Ex).ne.sum(Ex))write(*,*)"empart.f90: find NaN in Ex"
     if(sum(Ey).ne.sum(Ey))write(*,*)"empart.f90: find NaN in Ey"
     if(sum(Ez).ne.sum(Ez))write(*,*)"empart.f90: find NaN in Ez"
     if(sum(Bx).ne.sum(Bx))write(*,*)"empart.f90: find NaN in Bx"
     if(sum(By).ne.sum(By))write(*,*)"empart.f90: find NaN in By"
     if(sum(Bz).ne.sum(Bz))write(*,*)"empart.f90: find NaN in Bz"
     if(sum(Kpsi).ne.sum(Kpsi))write(*,*)"empart.f90: find NaN in Kpsi"
     if(sum(Kphi).ne.sum(Kphi))write(*,*)"empart.f90: find NaN in Kphi"
     gont = 1
     return
  endif

  PI = dacos(-ONE)

  dX = crho(2) - crho(1)
  dY = sigma(2) - sigma(1)
  dZ = R(2) - R(1)

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ext(3)
  do j=1,ext(2)
  do i=1,ext(1)
  alpn1(i,j,k) = Lap(i,j,k) + ONE
  chin1(i,j,k) = chi(i,j,k) + ONE
  chi3o2(i,j,k)  = dsqrt(chin1(i,j,k))**3
  gxx(i,j,k) = dxx(i,j,k) + ONE
  gyy(i,j,k) = dyy(i,j,k) + ONE
  gzz(i,j,k) = dzz(i,j,k) + ONE
  gxy(i,j,k) = dxy(i,j,k)
  gxz(i,j,k) = dxz(i,j,k)
  gyz(i,j,k) = dyz(i,j,k)
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

  call fderivs_shc(ext,Lap,Lapx,Lapy,Lapz,crho,sigma,R, SYM, SYM,SYM,Symmetry,Lev,sst,          &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ext,betax,betaxx,betaxy,betaxz,crho,sigma,R,ANTI, SYM, SYM,Symmetry,Lev,sst, &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ext,betay,betayx,betayy,betayz,crho,sigma,R, SYM,ANTI, SYM,Symmetry,Lev,sst, &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ext,betaz,betazx,betazy,betazz,crho,sigma,R, SYM, SYM,ANTI,Symmetry,Lev,sst, &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
 
  call fderivs_shc(ext,chi,chix,chiy,chiz,crho,sigma,R, SYM, SYM,SYM,Symmetry,Lev,sst,          &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)

  call fderivs_shc(ext,dxx,gxxx,gxxy,gxxz,crho,sigma,R, SYM, SYM,SYM,Symmetry,Lev,sst,          &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ext,gxy,gxyx,gxyy,gxyz,crho,sigma,R,ANTI,ANTI,SYM,Symmetry,Lev,sst,          &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ext,gxz,gxzx,gxzy,gxzz,crho,sigma,R,ANTI,SYM ,ANTI,Symmetry,Lev,sst,         &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ext,dyy,gyyx,gyyy,gyyz,crho,sigma,R, SYM, SYM,SYM,Symmetry,Lev,sst,          &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ext,gyz,gyzx,gyzy,gyzz,crho,sigma,R,SYM ,ANTI,ANTI,Symmetry,Lev,sst,         &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ext,dzz,gzzx,gzzy,gzzz,crho,sigma,R, SYM, SYM,SYM,Symmetry,Lev,sst,          &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)

  call fderivs_shc(ext,Kpsi,Kpsix,Kpsiy,Kpsiz,crho,sigma,R, SYM, SYM,SYM,Symmetry,Lev,sst,          &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ext,Kphi,Kphix,Kphiy,Kphiz,crho,sigma,R, SYM, SYM,SYM,Symmetry,Lev,sst,          &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)

  call fderivs_shc(ext,Ex,Exx,Exy,Exz,crho,sigma,R,ANTI, SYM, SYM,Symmetry,Lev,sst, &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ext,Ey,Eyx,Eyy,Eyz,crho,sigma,R, SYM,ANTI, SYM,Symmetry,Lev,sst, &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ext,Ez,Ezx,Ezy,Ezz,crho,sigma,R, SYM, SYM,ANTI,Symmetry,Lev,sst, &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)

#if 1    
  call fderivs_shc(ext,Bx,Bxx,Bxy,Bxz,crho,sigma,R, SYM,ANTI,ANTI,Symmetry,Lev,sst, &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ext,By,Byx,Byy,Byz,crho,sigma,R,ANTI, SYM,ANTI,Symmetry,Lev,sst, &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ext,Bz,Bzx,Bzy,Bzz,crho,sigma,R,ANTI,ANTI, SYM,Symmetry,Lev,sst, &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
#else 
  call fderivs_shc(ext,Bx,Bxx,Bxy,Bxz,crho,sigma,R,ANTI, SYM, SYM,Symmetry,Lev,sst, &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ext,By,Byx,Byy,Byz,crho,sigma,R, SYM,ANTI, SYM,Symmetry,Lev,sst, &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ext,Bz,Bzx,Bzy,Bzz,crho,sigma,R, SYM, SYM,ANTI,Symmetry,Lev,sst, &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
#endif
! check axal vector
! physical gij
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ext(3)
  do j=1,ext(2)
  do i=1,ext(1)
  gxx(i,j,k) = gxx(i,j,k)/chin1(i,j,k)
  gxy(i,j,k) = gxy(i,j,k)/chin1(i,j,k)
  gxz(i,j,k) = gxz(i,j,k)/chin1(i,j,k)
  gyy(i,j,k) = gyy(i,j,k)/chin1(i,j,k)
  gyz(i,j,k) = gyz(i,j,k)/chin1(i,j,k)
  gzz(i,j,k) = gzz(i,j,k)/chin1(i,j,k)
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO
!physical gij,k  
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ext(3)
  do j=1,ext(2)
  do i=1,ext(1)
  gxxx(i,j,k) = (gxxx(i,j,k)-gxx(i,j,k)*chix(i,j,k))/chin1(i,j,k)
  gxxy(i,j,k) = (gxxy(i,j,k)-gxx(i,j,k)*chiy(i,j,k))/chin1(i,j,k)
  gxxz(i,j,k) = (gxxz(i,j,k)-gxx(i,j,k)*chiz(i,j,k))/chin1(i,j,k)
  gxyx(i,j,k) = (gxyx(i,j,k)-gxy(i,j,k)*chix(i,j,k))/chin1(i,j,k)
  gxyy(i,j,k) = (gxyy(i,j,k)-gxy(i,j,k)*chiy(i,j,k))/chin1(i,j,k)
  gxyz(i,j,k) = (gxyz(i,j,k)-gxy(i,j,k)*chiz(i,j,k))/chin1(i,j,k)
  gxzx(i,j,k) = (gxzx(i,j,k)-gxz(i,j,k)*chix(i,j,k))/chin1(i,j,k)
  gxzy(i,j,k) = (gxzy(i,j,k)-gxz(i,j,k)*chiy(i,j,k))/chin1(i,j,k)
  gxzz(i,j,k) = (gxzz(i,j,k)-gxz(i,j,k)*chiz(i,j,k))/chin1(i,j,k)
  gyyx(i,j,k) = (gyyx(i,j,k)-gyy(i,j,k)*chix(i,j,k))/chin1(i,j,k)
  gyyy(i,j,k) = (gyyy(i,j,k)-gyy(i,j,k)*chiy(i,j,k))/chin1(i,j,k)
  gyyz(i,j,k) = (gyyz(i,j,k)-gyy(i,j,k)*chiz(i,j,k))/chin1(i,j,k)
  gyzx(i,j,k) = (gyzx(i,j,k)-gyz(i,j,k)*chix(i,j,k))/chin1(i,j,k)
  gyzy(i,j,k) = (gyzy(i,j,k)-gyz(i,j,k)*chiy(i,j,k))/chin1(i,j,k)
  gyzz(i,j,k) = (gyzz(i,j,k)-gyz(i,j,k)*chiz(i,j,k))/chin1(i,j,k)
  gzzx(i,j,k) = (gzzx(i,j,k)-gzz(i,j,k)*chix(i,j,k))/chin1(i,j,k)
  gzzy(i,j,k) = (gzzy(i,j,k)-gzz(i,j,k)*chiy(i,j,k))/chin1(i,j,k)
  gzzz(i,j,k) = (gzzz(i,j,k)-gzz(i,j,k)*chiz(i,j,k))/chin1(i,j,k)
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

! physical inverse metric
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ext(3)
  do j=1,ext(2)
  do i=1,ext(1)
  gupzz(i,j,k) =  gxx(i,j,k) * gyy(i,j,k) * gzz(i,j,k) + gxy(i,j,k) * gyz(i,j,k) * gxz(i,j,k) + gxz(i,j,k) * gxy(i,j,k) * gyz(i,j,k) - &
           gxz(i,j,k) * gyy(i,j,k) * gxz(i,j,k) - gxy(i,j,k) * gxy(i,j,k) * gzz(i,j,k) - gxx(i,j,k) * gyz(i,j,k) * gyz(i,j,k)
  gupxx(i,j,k) =   ( gyy(i,j,k) * gzz(i,j,k) - gyz(i,j,k) * gyz(i,j,k) ) / gupzz(i,j,k)
  gupxy(i,j,k) = - ( gxy(i,j,k) * gzz(i,j,k) - gyz(i,j,k) * gxz(i,j,k) ) / gupzz(i,j,k)
  gupxz(i,j,k) =   ( gxy(i,j,k) * gyz(i,j,k) - gyy(i,j,k) * gxz(i,j,k) ) / gupzz(i,j,k)
  gupyy(i,j,k) =   ( gxx(i,j,k) * gzz(i,j,k) - gxz(i,j,k) * gxz(i,j,k) ) / gupzz(i,j,k)
  gupyz(i,j,k) = - ( gxx(i,j,k) * gyz(i,j,k) - gxy(i,j,k) * gxz(i,j,k) ) / gupzz(i,j,k)
  gupzz(i,j,k) =   ( gxx(i,j,k) * gyy(i,j,k) - gxy(i,j,k) * gxy(i,j,k) ) / gupzz(i,j,k)
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ext(3)
  do j=1,ext(2)
  do i=1,ext(1)
  Ex_rhs(i,j,k) = alpn1(i,j,k)*trK(i,j,k)*Ex(i,j,k)-(Ex(i,j,k)*betaxx(i,j,k)+Ey(i,j,k)*betaxy(i,j,k)+Ez(i,j,k)*betaxz(i,j,k))                  &
          -FOUR*PI*alpn1(i,j,k)*Jx(i,j,k)-alpn1(i,j,k)*(gupxx(i,j,k)*Kpsix(i,j,k)+gupxy(i,j,k)*Kpsiy(i,j,k)+gupxz(i,j,k)*Kpsiz(i,j,k))  &
          +chi3o2(i,j,k)*(                                                      &
          ((gxz(i,j,k)*Bx(i,j,k)+gyz(i,j,k)*By(i,j,k)+gzz(i,j,k)*Bz(i,j,k))*Lapy(i,j,k)+alpn1(i,j,k)*(gxz(i,j,k)*Bxy(i,j,k)+gyz(i,j,k)*Byy(i,j,k)+gzz(i,j,k)*Bzy(i,j,k))+alpn1(i,j,k)*(Bx(i,j,k)*gxzy(i,j,k)+By(i,j,k)*gyzy(i,j,k)+Bz(i,j,k)*gzzy(i,j,k)))-&
          ((gxy(i,j,k)*Bx(i,j,k)+gyy(i,j,k)*By(i,j,k)+gyz(i,j,k)*Bz(i,j,k))*Lapz(i,j,k)+alpn1(i,j,k)*(gxy(i,j,k)*Bxz(i,j,k)+gyy(i,j,k)*Byz(i,j,k)+gyz(i,j,k)*Bzz(i,j,k))+alpn1(i,j,k)*(Bx(i,j,k)*gxyz(i,j,k)+By(i,j,k)*gyyz(i,j,k)+Bz(i,j,k)*gyzz(i,j,k))))
  Ey_rhs(i,j,k) = alpn1(i,j,k)*trK(i,j,k)*Ey(i,j,k)-(Ex(i,j,k)*betayx(i,j,k)+Ey(i,j,k)*betayy(i,j,k)+Ez(i,j,k)*betayz(i,j,k))                  &
          -FOUR*PI*alpn1(i,j,k)*Jy(i,j,k)-alpn1(i,j,k)*(gupxy(i,j,k)*Kpsix(i,j,k)+gupyy(i,j,k)*Kpsiy(i,j,k)+gupyz(i,j,k)*Kpsiz(i,j,k))  &
          +chi3o2(i,j,k)*(                                                      &
          ((gxx(i,j,k)*Bx(i,j,k)+gxy(i,j,k)*By(i,j,k)+gxz(i,j,k)*Bz(i,j,k))*Lapz(i,j,k)+alpn1(i,j,k)*(gxx(i,j,k)*Bxz(i,j,k)+gxy(i,j,k)*Byz(i,j,k)+gxz(i,j,k)*Bzz(i,j,k))+alpn1(i,j,k)*(Bx(i,j,k)*gxxz(i,j,k)+By(i,j,k)*gxyz(i,j,k)+Bz(i,j,k)*gxzz(i,j,k)))-&
          ((gxz(i,j,k)*Bx(i,j,k)+gyz(i,j,k)*By(i,j,k)+gzz(i,j,k)*Bz(i,j,k))*Lapx(i,j,k)+alpn1(i,j,k)*(gxz(i,j,k)*Bxx(i,j,k)+gyz(i,j,k)*Byx(i,j,k)+gzz(i,j,k)*Bzx(i,j,k))+alpn1(i,j,k)*(Bx(i,j,k)*gxzx(i,j,k)+By(i,j,k)*gyzx(i,j,k)+Bz(i,j,k)*gzzx(i,j,k))))
  Ez_rhs(i,j,k) = alpn1(i,j,k)*trK(i,j,k)*Ez(i,j,k)-(Ex(i,j,k)*betazx(i,j,k)+Ey(i,j,k)*betazy(i,j,k)+Ez(i,j,k)*betazz(i,j,k))                  &
          -FOUR*PI*alpn1(i,j,k)*Jz(i,j,k)-alpn1(i,j,k)*(gupxz(i,j,k)*Kpsix(i,j,k)+gupyz(i,j,k)*Kpsiy(i,j,k)+gupzz(i,j,k)*Kpsiz(i,j,k))  &
          +chi3o2(i,j,k)*(                                                      &
          ((gxy(i,j,k)*Bx(i,j,k)+gyy(i,j,k)*By(i,j,k)+gyz(i,j,k)*Bz(i,j,k))*Lapx(i,j,k)+alpn1(i,j,k)*(gxy(i,j,k)*Bxx(i,j,k)+gyy(i,j,k)*Byx(i,j,k)+gyz(i,j,k)*Bzx(i,j,k))+alpn1(i,j,k)*(Bx(i,j,k)*gxyx(i,j,k)+By(i,j,k)*gyyx(i,j,k)+Bz(i,j,k)*gyzx(i,j,k)))-&
          ((gxx(i,j,k)*Bx(i,j,k)+gxy(i,j,k)*By(i,j,k)+gxz(i,j,k)*Bz(i,j,k))*Lapy(i,j,k)+alpn1(i,j,k)*(gxx(i,j,k)*Bxy(i,j,k)+gxy(i,j,k)*Byy(i,j,k)+gxz(i,j,k)*Bzy(i,j,k))+alpn1(i,j,k)*(Bx(i,j,k)*gxxy(i,j,k)+By(i,j,k)*gxyy(i,j,k)+Bz(i,j,k)*gxzy(i,j,k))))

  Bx_rhs(i,j,k) = alpn1(i,j,k)*trK(i,j,k)*Bx(i,j,k)-(Bx(i,j,k)*betaxx(i,j,k)+By(i,j,k)*betaxy(i,j,k)+Bz(i,j,k)*betaxz(i,j,k))                  &
                           -alpn1(i,j,k)*(gupxx(i,j,k)*Kphix(i,j,k)+gupxy(i,j,k)*Kphiy(i,j,k)+gupxz(i,j,k)*Kphiz(i,j,k))  &
          -chi3o2(i,j,k)*(                                                      &
          ((gxz(i,j,k)*Ex(i,j,k)+gyz(i,j,k)*Ey(i,j,k)+gzz(i,j,k)*Ez(i,j,k))*Lapy(i,j,k)+alpn1(i,j,k)*(gxz(i,j,k)*Exy(i,j,k)+gyz(i,j,k)*Eyy(i,j,k)+gzz(i,j,k)*Ezy(i,j,k))+alpn1(i,j,k)*(Ex(i,j,k)*gxzy(i,j,k)+Ey(i,j,k)*gyzy(i,j,k)+Ez(i,j,k)*gzzy(i,j,k)))-&
          ((gxy(i,j,k)*Ex(i,j,k)+gyy(i,j,k)*Ey(i,j,k)+gyz(i,j,k)*Ez(i,j,k))*Lapz(i,j,k)+alpn1(i,j,k)*(gxy(i,j,k)*Exz(i,j,k)+gyy(i,j,k)*Eyz(i,j,k)+gyz(i,j,k)*Ezz(i,j,k))+alpn1(i,j,k)*(Ex(i,j,k)*gxyz(i,j,k)+Ey(i,j,k)*gyyz(i,j,k)+Ez(i,j,k)*gyzz(i,j,k))))
  By_rhs(i,j,k) = alpn1(i,j,k)*trK(i,j,k)*By(i,j,k)-(Bx(i,j,k)*betayx(i,j,k)+By(i,j,k)*betayy(i,j,k)+Bz(i,j,k)*betayz(i,j,k))                  &
                           -alpn1(i,j,k)*(gupxy(i,j,k)*Kphix(i,j,k)+gupyy(i,j,k)*Kphiy(i,j,k)+gupyz(i,j,k)*Kphiz(i,j,k))  &
          -chi3o2(i,j,k)*(                                                      &
          ((gxx(i,j,k)*Ex(i,j,k)+gxy(i,j,k)*Ey(i,j,k)+gxz(i,j,k)*Ez(i,j,k))*Lapz(i,j,k)+alpn1(i,j,k)*(gxx(i,j,k)*Exz(i,j,k)+gxy(i,j,k)*Eyz(i,j,k)+gxz(i,j,k)*Ezz(i,j,k))+alpn1(i,j,k)*(Ex(i,j,k)*gxxz(i,j,k)+Ey(i,j,k)*gxyz(i,j,k)+Ez(i,j,k)*gxzz(i,j,k)))-&
          ((gxz(i,j,k)*Ex(i,j,k)+gyz(i,j,k)*Ey(i,j,k)+gzz(i,j,k)*Ez(i,j,k))*Lapx(i,j,k)+alpn1(i,j,k)*(gxz(i,j,k)*Exx(i,j,k)+gyz(i,j,k)*Eyx(i,j,k)+gzz(i,j,k)*Ezx(i,j,k))+alpn1(i,j,k)*(Ex(i,j,k)*gxzx(i,j,k)+Ey(i,j,k)*gyzx(i,j,k)+Ez(i,j,k)*gzzx(i,j,k))))
  Bz_rhs(i,j,k) = alpn1(i,j,k)*trK(i,j,k)*Bz(i,j,k)-(Bx(i,j,k)*betazx(i,j,k)+By(i,j,k)*betazy(i,j,k)+Bz(i,j,k)*betazz(i,j,k))                  &
                           -alpn1(i,j,k)*(gupxz(i,j,k)*Kphix(i,j,k)+gupyz(i,j,k)*Kphiy(i,j,k)+gupzz(i,j,k)*Kphiz(i,j,k))  &
          -chi3o2(i,j,k)*(                                                      &
          ((gxy(i,j,k)*Ex(i,j,k)+gyy(i,j,k)*Ey(i,j,k)+gyz(i,j,k)*Ez(i,j,k))*Lapx(i,j,k)+alpn1(i,j,k)*(gxy(i,j,k)*Exx(i,j,k)+gyy(i,j,k)*Eyx(i,j,k)+gyz(i,j,k)*Ezx(i,j,k))+alpn1(i,j,k)*(Ex(i,j,k)*gxyx(i,j,k)+Ey(i,j,k)*gyyx(i,j,k)+Ez(i,j,k)*gyzx(i,j,k)))-&
          ((gxx(i,j,k)*Ex(i,j,k)+gxy(i,j,k)*Ey(i,j,k)+gxz(i,j,k)*Ez(i,j,k))*Lapy(i,j,k)+alpn1(i,j,k)*(gxx(i,j,k)*Exy(i,j,k)+gxy(i,j,k)*Eyy(i,j,k)+gxz(i,j,k)*Ezy(i,j,k))+alpn1(i,j,k)*(Ex(i,j,k)*gxxy(i,j,k)+Ey(i,j,k)*gxyy(i,j,k)+Ez(i,j,k)*gxzy(i,j,k))))

  Kpsi_rhs(i,j,k) = FOUR*PI*alpn1(i,j,k)*qchar(i,j,k)-alpn1(i,j,k)*kappa*Kpsi(i,j,k) - &
            alpn1(i,j,k)*(Exx(i,j,k)+Eyy(i,j,k)+Ezz(i,j,k)-F3o2/chin1(i,j,k)*(chix(i,j,k)*Ex(i,j,k)+chiy(i,j,k)*Ey(i,j,k)+chiz(i,j,k)*Ez(i,j,k))) 
  Kphi_rhs(i,j,k) =                    -alpn1(i,j,k)*kappa*Kphi(i,j,k) - &
            alpn1(i,j,k)*(Bxx(i,j,k)+Byy(i,j,k)+Bzz(i,j,k)-F3o2/chin1(i,j,k)*(chix(i,j,k)*Bx(i,j,k)+chiy(i,j,k)*By(i,j,k)+chiz(i,j,k)*Bz(i,j,k))) 
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

  SSS(1)=SYM
  SSS(2)=SYM
  SSS(3)=SYM

  AAS(1)=ANTI
  AAS(2)=ANTI
  AAS(3)=SYM

  ASA(1)=ANTI
  ASA(2)=SYM
  ASA(3)=ANTI

  SAA(1)=SYM
  SAA(2)=ANTI
  SAA(3)=ANTI

  ASS(1)=ANTI
  ASS(2)=SYM
  ASS(3)=SYM

  SAS(1)=SYM
  SAS(2)=ANTI
  SAS(3)=SYM

  SSA(1)=SYM
  SSA(2)=SYM
  SSA(3)=ANTI

!!!!!!!!!advection term part
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ext(3)
  do j=1,ext(2)
  do i=1,ext(1)
  Kpsi_rhs(i,j,k) = Kpsi_rhs(i,j,k) + betax(i,j,k)*Kpsix(i,j,k)+betay(i,j,k)*Kpsiy(i,j,k)+betaz(i,j,k)*Kpsiz(i,j,k)
  Kphi_rhs(i,j,k) = Kphi_rhs(i,j,k) + betax(i,j,k)*Kphix(i,j,k)+betay(i,j,k)*Kphiy(i,j,k)+betaz(i,j,k)*Kphiz(i,j,k)

  Ex_rhs(i,j,k) = Ex_rhs(i,j,k) + betax(i,j,k)*Exx(i,j,k)+betay(i,j,k)*Exy(i,j,k)+betaz(i,j,k)*Exz(i,j,k)
  Ey_rhs(i,j,k) = Ey_rhs(i,j,k) + betax(i,j,k)*Eyx(i,j,k)+betay(i,j,k)*Eyy(i,j,k)+betaz(i,j,k)*Eyz(i,j,k)
  Ez_rhs(i,j,k) = Ez_rhs(i,j,k) + betax(i,j,k)*Ezx(i,j,k)+betay(i,j,k)*Ezy(i,j,k)+betaz(i,j,k)*Ezz(i,j,k)

  Bx_rhs(i,j,k) = Bx_rhs(i,j,k) + betax(i,j,k)*Bxx(i,j,k)+betay(i,j,k)*Bxy(i,j,k)+betaz(i,j,k)*Bxz(i,j,k)
  By_rhs(i,j,k) = By_rhs(i,j,k) + betax(i,j,k)*Byx(i,j,k)+betay(i,j,k)*Byy(i,j,k)+betaz(i,j,k)*Byz(i,j,k)
  Bz_rhs(i,j,k) = Bz_rhs(i,j,k) + betax(i,j,k)*Bzx(i,j,k)+betay(i,j,k)*Bzy(i,j,k)+betaz(i,j,k)*Bzz(i,j,k)
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

! numerical dissipation part
  if(eps>0)then 
! usual Kreiss-Oliger dissipation 

  call kodis_sh(ext,crho,sigma,R,Kpsi,Kpsi_rhs,SSS,Symmetry,eps,sst)
  call kodis_sh(ext,crho,sigma,R,Kphi,Kphi_rhs,SSS,Symmetry,eps,sst)
  call kodis_sh(ext,crho,sigma,R,Ex,Ex_rhs,ASS,Symmetry,eps,sst)
  call kodis_sh(ext,crho,sigma,R,Ey,Ey_rhs,SAS,Symmetry,eps,sst)
  call kodis_sh(ext,crho,sigma,R,Ez,Ez_rhs,SSA,Symmetry,eps,sst)
  call kodis_sh(ext,crho,sigma,R,Bx,Bx_rhs,SAA,Symmetry,eps,sst)
  call kodis_sh(ext,crho,sigma,R,By,By_rhs,ASA,Symmetry,eps,sst)
  call kodis_sh(ext,crho,sigma,R,Bz,Bz_rhs,AAS,Symmetry,eps,sst)

  endif
! stress-energy tensor
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(SHARED) PRIVATE(i,j,k)
  do k=1,ext(3)
  do j=1,ext(2)
  do i=1,ext(1)
  rho(i,j,k) = (gxx(i,j,k)*(Ex(i,j,k)*Ex(i,j,k)+Bx(i,j,k)*Bx(i,j,k))+gyy(i,j,k)*(Ey(i,j,k)*Ey(i,j,k)+By(i,j,k)*By(i,j,k))+gzz(i,j,k)*(Ez(i,j,k)*Ez(i,j,k)+Bz(i,j,k)*Bz(i,j,k)) + &
      +TWO*(gxy(i,j,k)*(Ex(i,j,k)*Ey(i,j,k)+Bx(i,j,k)*By(i,j,k))+gxz(i,j,k)*(Ex(i,j,k)*Ez(i,j,k)+Bx(i,j,k)*Bz(i,j,k))+gyz(i,j,k)*(Ey(i,j,k)*Ez(i,j,k)+By(i,j,k)*Bz(i,j,k))))/EIT/PI
  Sx(i,j,k) = (Ey(i,j,k)*Bz(i,j,k)-Ez(i,j,k)*By(i,j,k))/FOUR/PI/chi3o2(i,j,k)
  Sy(i,j,k) = (Ez(i,j,k)*Bx(i,j,k)-Ex(i,j,k)*Bz(i,j,k))/FOUR/PI/chi3o2(i,j,k)
  Sz(i,j,k) = (Ex(i,j,k)*By(i,j,k)-Ey(i,j,k)*Bx(i,j,k))/FOUR/PI/chi3o2(i,j,k)
  lEx(i,j,k) = gxx(i,j,k)*Ex(i,j,k)+gxy(i,j,k)*Ey(i,j,k)+gxz(i,j,k)*Ez(i,j,k)
  lEy(i,j,k) = gxy(i,j,k)*Ex(i,j,k)+gyy(i,j,k)*Ey(i,j,k)+gyz(i,j,k)*Ez(i,j,k)
  lEz(i,j,k) = gxz(i,j,k)*Ex(i,j,k)+gyz(i,j,k)*Ey(i,j,k)+gzz(i,j,k)*Ez(i,j,k)
  lBx(i,j,k) = gxx(i,j,k)*Bx(i,j,k)+gxy(i,j,k)*By(i,j,k)+gxz(i,j,k)*Bz(i,j,k)
  lBy(i,j,k) = gxy(i,j,k)*Bx(i,j,k)+gyy(i,j,k)*By(i,j,k)+gyz(i,j,k)*Bz(i,j,k)
  lBz(i,j,k) = gxz(i,j,k)*Bx(i,j,k)+gyz(i,j,k)*By(i,j,k)+gzz(i,j,k)*Bz(i,j,k)
  Sxx(i,j,k) = rho(i,j,k)*gxx(i,j,k)-(lEx(i,j,k)*lEx(i,j,k)+lBx(i,j,k)*lBx(i,j,k))/FOUR/PI
  Sxy(i,j,k) = rho(i,j,k)*gxy(i,j,k)-(lEx(i,j,k)*lEy(i,j,k)+lBx(i,j,k)*lBy(i,j,k))/FOUR/PI
  Sxz(i,j,k) = rho(i,j,k)*gxz(i,j,k)-(lEx(i,j,k)*lEz(i,j,k)+lBx(i,j,k)*lBz(i,j,k))/FOUR/PI
  Syy(i,j,k) = rho(i,j,k)*gyy(i,j,k)-(lEy(i,j,k)*lEy(i,j,k)+lBy(i,j,k)*lBy(i,j,k))/FOUR/PI
  Syz(i,j,k) = rho(i,j,k)*gyz(i,j,k)-(lEy(i,j,k)*lEz(i,j,k)+lBy(i,j,k)*lBz(i,j,k))/FOUR/PI
  Szz(i,j,k) = rho(i,j,k)*gzz(i,j,k)-(lEz(i,j,k)*lEz(i,j,k)+lBz(i,j,k)*lBz(i,j,k))/FOUR/PI
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO

  gont = 0

  return

  end function compute_rhs_empart_ss
