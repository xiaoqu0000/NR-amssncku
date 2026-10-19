
! PIN==0: standard scalar wave
! PIN==1: \block phi = \eta(dphi,dphi)
#define PIN 0

  function compute_rhs_scalar(ex, T, X, Y, Z,                           &
               Sphi,Spi,Sphi_rhs,Spi_rhs,                               &
               Symmetry,Lev,eps) result(gont)
  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3), Symmetry,Lev
  real*8, intent(in ):: T,X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Sphi,Spi
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Sphi_rhs,Spi_rhs
  real*8,intent(in) :: eps
!  gont = 0: success; gont = 1: something wrong
  integer::gont

!~~~~~~> Other variables:

  real*8, dimension(ex(1),ex(2),ex(3)) :: fxx,fxy,fxz,fyy,fyz,fzz
  real*8,dimension(3) ::SSS
  real*8, parameter :: HALF = 0.5d0, ONE = 1.D0, TWO = 2.D0, FOUR = 4.D0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  real*8 :: tt

!!! sanity check
  tt = sum(Sphi)+sum(Spi)
  if(tt.ne.tt) then
     if(sum(Sphi).ne.sum(Sphi))write(*,*)"scalar_rhs.f90:compute_rhs_scalar find NaN in Sphi"
     if(sum(Spi).ne.sum(Spi))write(*,*)"scalar_rhs.f90:compute_rhs_scalar find NaN in Spi"
     gont = 1
     return
  endif

  Sphi_rhs = Spi !rhs for phi

#if   (PIN == 0)  
  call fdderivs(ex,Sphi,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)
   Spi_rhs =   fxx + fyy + fzz
#elif (PIN == 1)
  call fdderivs(ex,Sphi,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)
   Spi_rhs =   Spi*Spi + fxx + fyy + fzz
  call fderivs(ex,Sphi,fxx,fyy,fzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)
   Spi_rhs = Spi_rhs - (fxx*fxx+fyy*fyy+fzz*fzz)
#endif
  if(eps>0)then 
! usual Kreiss-Oliger dissipation      
  SSS(1)=SYM
  SSS(2)=SYM
  SSS(3)=SYM

  call kodis(ex,X,Y,Z,Sphi,Sphi_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Spi,Spi_rhs,SSS,Symmetry,eps)
  endif

  gont = 0

  return

  end function compute_rhs_scalar
! for shell
  function compute_rhs_scalar_ss(ex, T,crho,sigma,R,x,y,z,                   &
               drhodx, drhody, drhodz,                                         &
               dsigmadx,dsigmady,dsigmadz,                                     &
               dRdx,dRdy,dRdz,                                                 &
               drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                &
               dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,    &
               dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,                            &
               Sphi,Spi,Sphi_rhs,Spi_rhs,                                      &
               Symmetry,Lev,eps,sst) result(gont)
  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3), Symmetry,Lev,sst
  real*8, intent(in ):: T
  double precision,intent(in),dimension(ex(1))::crho
  double precision,intent(in),dimension(ex(2))::sigma
  double precision,intent(in),dimension(ex(3))::R
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::x,y,z
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::drhodx, drhody, drhodz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::dsigmadx,dsigmady,dsigmadz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::dRdx,dRdy,dRdz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Sphi,Spi
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Sphi_rhs,Spi_rhs
  real*8,intent(in) :: eps
!  gont = 0: success; gont = 1: something wrong
  integer::gont

!~~~~~~> Other variables:

  real*8, dimension(ex(1),ex(2),ex(3)) :: fxx,fxy,fxz,fyy,fyz,fzz
  real*8,dimension(3) ::SSS
  real*8, parameter :: HALF = 0.5d0, ONE = 1.D0, TWO = 2.D0, FOUR = 4.D0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  real*8 :: tt

!!! sanity check
  tt = sum(Sphi)+sum(Spi)
  if(tt.ne.tt) then
     if(sum(Sphi).ne.sum(Sphi))write(*,*)"scalar_rhs.f90:compute_rhs_scalar_ss find NaN in Sphi"
     if(sum(Spi).ne.sum(Spi))write(*,*)"scalar_rhs.f90:compute_rhs_scalar_ss find NaN in Spi"
     gont = 1
     return
  endif

  Sphi_rhs = Spi !rhs for phi

#if   (PIN == 0)  
  call fdderivs_shc(ex,Sphi,fxx,fxy,fxz,fyy,fyz,fzz,crho,sigma,R,SYM ,SYM ,SYM ,Symmetry,Lev,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz)

  Spi_rhs = fxx+fyy+fzz
#elif (PIN == 1)
  call fdderivs_shc(ex,Sphi,fxx,fxy,fxz,fyy,fyz,fzz,crho,sigma,R,SYM ,SYM ,SYM ,Symmetry,Lev,sst, &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz)
   Spi_rhs =   Spi*Spi + fxx + fyy + fzz
  call fderivs_shc(ex,Sphi,fxx,fyy,fzz,crho,sigma,R,SYM,SYM,SYM,Symmetry,Lev,sst,  &
                       drhodx, drhody, drhodz,                                     &
                       dsigmadx,dsigmady,dsigmadz,                                 &
                       dRdx,dRdy,dRdz)
   Spi_rhs = Spi_rhs - (fxx*fxx+fyy*fyy+fzz*fzz)
#endif

  if(eps>0)then 
! usual Kreiss-Oliger dissipation      
  SSS(1)=SYM
  SSS(2)=SYM
  SSS(3)=SYM

  call kodis_sh(ex,crho,sigma,R,Sphi,Sphi_rhs,SSS,Symmetry,eps,sst)
  call kodis_sh(ex,crho,sigma,R,Spi,Spi_rhs,SSS,Symmetry,eps,sst)
  endif

  gont = 0

  return

  end function compute_rhs_scalar_ss
