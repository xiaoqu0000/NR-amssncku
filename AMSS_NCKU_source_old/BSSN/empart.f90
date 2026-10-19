
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

  alpn1 = Lap + ONE
  chin1 = chi + ONE
  chi3o2  = dsqrt(chin1)**3
  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE
  gxy = dxy
  gxz = dxz
  gyz = dyz

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
  gxx = gxx/chin1
  gxy = gxy/chin1
  gxz = gxz/chin1
  gyy = gyy/chin1
  gyz = gyz/chin1
  gzz = gzz/chin1
!physical gij,k  
  gxxx = (gxxx-gxx*chix)/chin1
  gxxy = (gxxy-gxx*chiy)/chin1
  gxxz = (gxxz-gxx*chiz)/chin1
  gxyx = (gxyx-gxy*chix)/chin1
  gxyy = (gxyy-gxy*chiy)/chin1
  gxyz = (gxyz-gxy*chiz)/chin1
  gxzx = (gxzx-gxz*chix)/chin1
  gxzy = (gxzy-gxz*chiy)/chin1
  gxzz = (gxzz-gxz*chiz)/chin1
  gyyx = (gyyx-gyy*chix)/chin1
  gyyy = (gyyy-gyy*chiy)/chin1
  gyyz = (gyyz-gyy*chiz)/chin1
  gyzx = (gyzx-gyz*chix)/chin1
  gyzy = (gyzy-gyz*chiy)/chin1
  gyzz = (gyzz-gyz*chiz)/chin1
  gzzx = (gzzx-gzz*chix)/chin1
  gzzy = (gzzy-gzz*chiy)/chin1
  gzzz = (gzzz-gzz*chiz)/chin1

! physical inverse metric
  gupzz =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz - &
           gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz
  gupxx =   ( gyy * gzz - gyz * gyz ) / gupzz
  gupxy = - ( gxy * gzz - gyz * gxz ) / gupzz
  gupxz =   ( gxy * gyz - gyy * gxz ) / gupzz
  gupyy =   ( gxx * gzz - gxz * gxz ) / gupzz
  gupyz = - ( gxx * gyz - gxy * gxz ) / gupzz
  gupzz =   ( gxx * gyy - gxy * gxy ) / gupzz

  Ex_rhs = alpn1*trK*Ex-(Ex*betaxx+Ey*betaxy+Ez*betaxz)                  &
          -FOUR*PI*alpn1*Jx-alpn1*(gupxx*Kpsix+gupxy*Kpsiy+gupxz*Kpsiz)  &
          +chi3o2*(                                                      &
          ((gxz*Bx+gyz*By+gzz*Bz)*Lapy+alpn1*(gxz*Bxy+gyz*Byy+gzz*Bzy)+alpn1*(Bx*gxzy+By*gyzy+Bz*gzzy))-&
          ((gxy*Bx+gyy*By+gyz*Bz)*Lapz+alpn1*(gxy*Bxz+gyy*Byz+gyz*Bzz)+alpn1*(Bx*gxyz+By*gyyz+Bz*gyzz)))
  Ey_rhs = alpn1*trK*Ey-(Ex*betayx+Ey*betayy+Ez*betayz)                  &
          -FOUR*PI*alpn1*Jy-alpn1*(gupxy*Kpsix+gupyy*Kpsiy+gupyz*Kpsiz)  &
          +chi3o2*(                                                      &
          ((gxx*Bx+gxy*By+gxz*Bz)*Lapz+alpn1*(gxx*Bxz+gxy*Byz+gxz*Bzz)+alpn1*(Bx*gxxz+By*gxyz+Bz*gxzz))-&
          ((gxz*Bx+gyz*By+gzz*Bz)*Lapx+alpn1*(gxz*Bxx+gyz*Byx+gzz*Bzx)+alpn1*(Bx*gxzx+By*gyzx+Bz*gzzx)))
  Ez_rhs = alpn1*trK*Ez-(Ex*betazx+Ey*betazy+Ez*betazz)                  &
          -FOUR*PI*alpn1*Jz-alpn1*(gupxz*Kpsix+gupyz*Kpsiy+gupzz*Kpsiz)  &
          +chi3o2*(                                                      &
          ((gxy*Bx+gyy*By+gyz*Bz)*Lapx+alpn1*(gxy*Bxx+gyy*Byx+gyz*Bzx)+alpn1*(Bx*gxyx+By*gyyx+Bz*gyzx))-&
          ((gxx*Bx+gxy*By+gxz*Bz)*Lapy+alpn1*(gxx*Bxy+gxy*Byy+gxz*Bzy)+alpn1*(Bx*gxxy+By*gxyy+Bz*gxzy)))

  Bx_rhs = alpn1*trK*Bx-(Bx*betaxx+By*betaxy+Bz*betaxz)                  &
                           -alpn1*(gupxx*Kphix+gupxy*Kphiy+gupxz*Kphiz)  &
          -chi3o2*(                                                      &
          ((gxz*Ex+gyz*Ey+gzz*Ez)*Lapy+alpn1*(gxz*Exy+gyz*Eyy+gzz*Ezy)+alpn1*(Ex*gxzy+Ey*gyzy+Ez*gzzy))-&
          ((gxy*Ex+gyy*Ey+gyz*Ez)*Lapz+alpn1*(gxy*Exz+gyy*Eyz+gyz*Ezz)+alpn1*(Ex*gxyz+Ey*gyyz+Ez*gyzz)))
  By_rhs = alpn1*trK*By-(Bx*betayx+By*betayy+Bz*betayz)                  &
                           -alpn1*(gupxy*Kphix+gupyy*Kphiy+gupyz*Kphiz)  &
          -chi3o2*(                                                      &
          ((gxx*Ex+gxy*Ey+gxz*Ez)*Lapz+alpn1*(gxx*Exz+gxy*Eyz+gxz*Ezz)+alpn1*(Ex*gxxz+Ey*gxyz+Ez*gxzz))-&
          ((gxz*Ex+gyz*Ey+gzz*Ez)*Lapx+alpn1*(gxz*Exx+gyz*Eyx+gzz*Ezx)+alpn1*(Ex*gxzx+Ey*gyzx+Ez*gzzx)))
  Bz_rhs = alpn1*trK*Bz-(Bx*betazx+By*betazy+Bz*betazz)                  &
                           -alpn1*(gupxz*Kphix+gupyz*Kphiy+gupzz*Kphiz)  &
          -chi3o2*(                                                      &
          ((gxy*Ex+gyy*Ey+gyz*Ez)*Lapx+alpn1*(gxy*Exx+gyy*Eyx+gyz*Ezx)+alpn1*(Ex*gxyx+Ey*gyyx+Ez*gyzx))-&
          ((gxx*Ex+gxy*Ey+gxz*Ez)*Lapy+alpn1*(gxx*Exy+gxy*Eyy+gxz*Ezy)+alpn1*(Ex*gxxy+Ey*gxyy+Ez*gxzy)))

  Kpsi_rhs = FOUR*PI*alpn1*qchar-alpn1*kappa*Kpsi - &
            alpn1*(Exx+Eyy+Ezz-F3o2/chin1*(chix*Ex+chiy*Ey+chiz*Ez)) 
  Kphi_rhs =                    -alpn1*kappa*Kphi - &
            alpn1*(Bxx+Byy+Bzz-F3o2/chin1*(chix*Bx+chiy*By+chiz*Bz)) 

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
  rho = (gxx*(Ex*Ex+Bx*Bx)+gyy*(Ey*Ey+By*By)+gzz*(Ez*Ez+Bz*Bz) + &
      +TWO*(gxy*(Ex*Ey+Bx*By)+gxz*(Ex*Ez+Bx*Bz)+gyz*(Ey*Ez+By*Bz)))/EIT/PI
  Sx = (Ey*Bz-Ez*By)/FOUR/PI/chi3o2
  Sy = (Ez*Bx-Ex*Bz)/FOUR/PI/chi3o2
  Sz = (Ex*By-Ey*Bx)/FOUR/PI/chi3o2
  lEx = gxx*Ex+gxy*Ey+gxz*Ez
  lEy = gxy*Ex+gyy*Ey+gyz*Ez
  lEz = gxz*Ex+gyz*Ey+gzz*Ez
  lBx = gxx*Bx+gxy*By+gxz*Bz
  lBy = gxy*Bx+gyy*By+gyz*Bz
  lBz = gxz*Bx+gyz*By+gzz*Bz
  Sxx = rho*gxx-(lEx*lEx+lBx*lBx)/FOUR/PI
  Sxy = rho*gxy-(lEx*lEy+lBx*lBy)/FOUR/PI
  Sxz = rho*gxz-(lEx*lEz+lBx*lBz)/FOUR/PI
  Syy = rho*gyy-(lEy*lEy+lBy*lBy)/FOUR/PI
  Syz = rho*gyz-(lEy*lEz+lBy*lBz)/FOUR/PI
  Szz = rho*gzz-(lEz*lEz+lBz*lBz)/FOUR/PI

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

  alpn1 = Lap + ONE
  chin1 = chi + ONE
  chi3o2  = dsqrt(chin1)**3
  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE
  gxy = dxy
  gxz = dxz
  gyz = dyz

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
  gxx = gxx/chin1
  gxy = gxy/chin1
  gxz = gxz/chin1
  gyy = gyy/chin1
  gyz = gyz/chin1
  gzz = gzz/chin1
!physical gij,k  
  gxxx = (gxxx-gxx*chix)/chin1
  gxxy = (gxxy-gxx*chiy)/chin1
  gxxz = (gxxz-gxx*chiz)/chin1
  gxyx = (gxyx-gxy*chix)/chin1
  gxyy = (gxyy-gxy*chiy)/chin1
  gxyz = (gxyz-gxy*chiz)/chin1
  gxzx = (gxzx-gxz*chix)/chin1
  gxzy = (gxzy-gxz*chiy)/chin1
  gxzz = (gxzz-gxz*chiz)/chin1
  gyyx = (gyyx-gyy*chix)/chin1
  gyyy = (gyyy-gyy*chiy)/chin1
  gyyz = (gyyz-gyy*chiz)/chin1
  gyzx = (gyzx-gyz*chix)/chin1
  gyzy = (gyzy-gyz*chiy)/chin1
  gyzz = (gyzz-gyz*chiz)/chin1
  gzzx = (gzzx-gzz*chix)/chin1
  gzzy = (gzzy-gzz*chiy)/chin1
  gzzz = (gzzz-gzz*chiz)/chin1

! physical inverse metric
  gupzz =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz - &
           gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz
  gupxx =   ( gyy * gzz - gyz * gyz ) / gupzz
  gupxy = - ( gxy * gzz - gyz * gxz ) / gupzz
  gupxz =   ( gxy * gyz - gyy * gxz ) / gupzz
  gupyy =   ( gxx * gzz - gxz * gxz ) / gupzz
  gupyz = - ( gxx * gyz - gxy * gxz ) / gupzz
  gupzz =   ( gxx * gyy - gxy * gxy ) / gupzz

  Ex_rhs = alpn1*trK*Ex-(Ex*betaxx+Ey*betaxy+Ez*betaxz)                  &
          -FOUR*PI*alpn1*Jx-alpn1*(gupxx*Kpsix+gupxy*Kpsiy+gupxz*Kpsiz)  &
          +chi3o2*(                                                      &
          ((gxz*Bx+gyz*By+gzz*Bz)*Lapy+alpn1*(gxz*Bxy+gyz*Byy+gzz*Bzy)+alpn1*(Bx*gxzy+By*gyzy+Bz*gzzy))-&
          ((gxy*Bx+gyy*By+gyz*Bz)*Lapz+alpn1*(gxy*Bxz+gyy*Byz+gyz*Bzz)+alpn1*(Bx*gxyz+By*gyyz+Bz*gyzz)))
  Ey_rhs = alpn1*trK*Ey-(Ex*betayx+Ey*betayy+Ez*betayz)                  &
          -FOUR*PI*alpn1*Jy-alpn1*(gupxy*Kpsix+gupyy*Kpsiy+gupyz*Kpsiz)  &
          +chi3o2*(                                                      &
          ((gxx*Bx+gxy*By+gxz*Bz)*Lapz+alpn1*(gxx*Bxz+gxy*Byz+gxz*Bzz)+alpn1*(Bx*gxxz+By*gxyz+Bz*gxzz))-&
          ((gxz*Bx+gyz*By+gzz*Bz)*Lapx+alpn1*(gxz*Bxx+gyz*Byx+gzz*Bzx)+alpn1*(Bx*gxzx+By*gyzx+Bz*gzzx)))
  Ez_rhs = alpn1*trK*Ez-(Ex*betazx+Ey*betazy+Ez*betazz)                  &
          -FOUR*PI*alpn1*Jz-alpn1*(gupxz*Kpsix+gupyz*Kpsiy+gupzz*Kpsiz)  &
          +chi3o2*(                                                      &
          ((gxy*Bx+gyy*By+gyz*Bz)*Lapx+alpn1*(gxy*Bxx+gyy*Byx+gyz*Bzx)+alpn1*(Bx*gxyx+By*gyyx+Bz*gyzx))-&
          ((gxx*Bx+gxy*By+gxz*Bz)*Lapy+alpn1*(gxx*Bxy+gxy*Byy+gxz*Bzy)+alpn1*(Bx*gxxy+By*gxyy+Bz*gxzy)))

  Bx_rhs = alpn1*trK*Bx-(Bx*betaxx+By*betaxy+Bz*betaxz)                  &
                           -alpn1*(gupxx*Kphix+gupxy*Kphiy+gupxz*Kphiz)  &
          -chi3o2*(                                                      &
          ((gxz*Ex+gyz*Ey+gzz*Ez)*Lapy+alpn1*(gxz*Exy+gyz*Eyy+gzz*Ezy)+alpn1*(Ex*gxzy+Ey*gyzy+Ez*gzzy))-&
          ((gxy*Ex+gyy*Ey+gyz*Ez)*Lapz+alpn1*(gxy*Exz+gyy*Eyz+gyz*Ezz)+alpn1*(Ex*gxyz+Ey*gyyz+Ez*gyzz)))
  By_rhs = alpn1*trK*By-(Bx*betayx+By*betayy+Bz*betayz)                  &
                           -alpn1*(gupxy*Kphix+gupyy*Kphiy+gupyz*Kphiz)  &
          -chi3o2*(                                                      &
          ((gxx*Ex+gxy*Ey+gxz*Ez)*Lapz+alpn1*(gxx*Exz+gxy*Eyz+gxz*Ezz)+alpn1*(Ex*gxxz+Ey*gxyz+Ez*gxzz))-&
          ((gxz*Ex+gyz*Ey+gzz*Ez)*Lapx+alpn1*(gxz*Exx+gyz*Eyx+gzz*Ezx)+alpn1*(Ex*gxzx+Ey*gyzx+Ez*gzzx)))
  Bz_rhs = alpn1*trK*Bz-(Bx*betazx+By*betazy+Bz*betazz)                  &
                           -alpn1*(gupxz*Kphix+gupyz*Kphiy+gupzz*Kphiz)  &
          -chi3o2*(                                                      &
          ((gxy*Ex+gyy*Ey+gyz*Ez)*Lapx+alpn1*(gxy*Exx+gyy*Eyx+gyz*Ezx)+alpn1*(Ex*gxyx+Ey*gyyx+Ez*gyzx))-&
          ((gxx*Ex+gxy*Ey+gxz*Ez)*Lapy+alpn1*(gxx*Exy+gxy*Eyy+gxz*Ezy)+alpn1*(Ex*gxxy+Ey*gxyy+Ez*gxzy)))

  Kpsi_rhs = FOUR*PI*alpn1*qchar-alpn1*kappa*Kpsi - &
            alpn1*(Exx+Eyy+Ezz-F3o2/chin1*(chix*Ex+chiy*Ey+chiz*Ez)) 
  Kphi_rhs =                    -alpn1*kappa*Kphi - &
            alpn1*(Bxx+Byy+Bzz-F3o2/chin1*(chix*Bx+chiy*By+chiz*Bz)) 

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
  Kpsi_rhs = Kpsi_rhs + betax*Kpsix+betay*Kpsiy+betaz*Kpsiz
  Kphi_rhs = Kphi_rhs + betax*Kphix+betay*Kphiy+betaz*Kphiz

  Ex_rhs = Ex_rhs + betax*Exx+betay*Exy+betaz*Exz
  Ey_rhs = Ey_rhs + betax*Eyx+betay*Eyy+betaz*Eyz
  Ez_rhs = Ez_rhs + betax*Ezx+betay*Ezy+betaz*Ezz

  Bx_rhs = Bx_rhs + betax*Bxx+betay*Bxy+betaz*Bxz
  By_rhs = By_rhs + betax*Byx+betay*Byy+betaz*Byz
  Bz_rhs = Bz_rhs + betax*Bzx+betay*Bzy+betaz*Bzz

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
  rho = (gxx*(Ex*Ex+Bx*Bx)+gyy*(Ey*Ey+By*By)+gzz*(Ez*Ez+Bz*Bz) + &
      +TWO*(gxy*(Ex*Ey+Bx*By)+gxz*(Ex*Ez+Bx*Bz)+gyz*(Ey*Ez+By*Bz)))/EIT/PI
  Sx = (Ey*Bz-Ez*By)/FOUR/PI/chi3o2
  Sy = (Ez*Bx-Ex*Bz)/FOUR/PI/chi3o2
  Sz = (Ex*By-Ey*Bx)/FOUR/PI/chi3o2
  lEx = gxx*Ex+gxy*Ey+gxz*Ez
  lEy = gxy*Ex+gyy*Ey+gyz*Ez
  lEz = gxz*Ex+gyz*Ey+gzz*Ez
  lBx = gxx*Bx+gxy*By+gxz*Bz
  lBy = gxy*Bx+gyy*By+gyz*Bz
  lBz = gxz*Bx+gyz*By+gzz*Bz
  Sxx = rho*gxx-(lEx*lEx+lBx*lBx)/FOUR/PI
  Sxy = rho*gxy-(lEx*lEy+lBx*lBy)/FOUR/PI
  Sxz = rho*gxz-(lEx*lEz+lBx*lBz)/FOUR/PI
  Syy = rho*gyy-(lEy*lEy+lBy*lBy)/FOUR/PI
  Syz = rho*gyz-(lEy*lEz+lBy*lBz)/FOUR/PI
  Szz = rho*gzz-(lEz*lEz+lBz*lBz)/FOUR/PI

  gont = 0

  return

  end function compute_rhs_empart_ss
