
! for ADM variables
  subroutine adm_ricci_gamma(ex, X, Y, Z,                        &
               dxx    ,   gxy    ,   gxz    ,   dyy    ,   gyz    ,   dzz,&
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,&
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,&
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,&
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,&
               Symmetry)
  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3), Symmetry
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: dxx,gxy,gxz,dyy,gyz,dzz
! when out, physical second kind of connection  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamzyy, Gamzyz, Gamzzz
! when out, physical Ricci tensor  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Rxx,Rxy,Rxz,Ryy,Ryz,Rzz

!~~~~~~> Other variables:

  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gyy,gzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxx,gxyx,gxzx,gyyx,gyzx,gzzx
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxy,gxyy,gxzy,gyyy,gyzy,gzzy
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxz,gxyz,gxzz,gyyz,gyzz,gzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: ass_Gamxxx,ass_Gamxxy,ass_Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: ass_Gamxyy,ass_Gamxyz,ass_Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: ass_Gamyxx,ass_Gamyxy,ass_Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: ass_Gamyyy,ass_Gamyyz,ass_Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: ass_Gamzxx,ass_Gamzxy,ass_Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: ass_Gamzyy,ass_Gamzyz,ass_Gamzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupyy,gupyz,gupzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxxx,gxxxy,gxxxz,gxxyy,gxxyz,gxxzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxyxx,gxyxy,gxyxz,gxyyy,gxyyz,gxyzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxzxx,gxzxy,gxzxz,gxzyy,gxzyz,gxzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyyxx,gyyxy,gyyxz,gyyyy,gyyyz,gyyzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyzxx,gyzxy,gyzxz,gyzyy,gyzyz,gyzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gzzxx,gzzxy,gzzxz,gzzyy,gzzyz,gzzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Rxyxy,  Rxyxz,  Rxyyz,  Rxzxz,  Rxzyz,  Ryzyz

  real*8, parameter :: ONE = 1.D0, TWO = 2.D0, FOUR = 4.D0
  real*8, parameter :: HALF = 0.5D0, F2o3 = 2.d0/3.d0, F3o2 = 1.5d0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0

  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE
 
  call fderivs(ex,dxx,gxxx,gxxy,gxxz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)
  call fderivs(ex,gxy,gxyx,gxyy,gxyz,X,Y,Z,ANTI,ANTI,SYM ,Symmetry,0)
  call fderivs(ex,gxz,gxzx,gxzy,gxzz,X,Y,Z,ANTI,SYM ,ANTI,Symmetry,0)
  call fderivs(ex,dyy,gyyx,gyyy,gyyz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)
  call fderivs(ex,gyz,gyzx,gyzy,gyzz,X,Y,Z,SYM ,ANTI,ANTI,Symmetry,0)
  call fderivs(ex,dzz,gzzx,gzzy,gzzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)

  call kind1_connection(ex,                      &
                        gxxx,gxyx,gxzx,gyyx,gyzx,gzzx, &
                        gxxy,gxyy,gxzy,gyyy,gyzy,gzzy, &
                        gxxz,gxyz,gxzz,gyyz,gyzz,gzzz, &
                        ass_Gamxxx, ass_Gamxxy, ass_Gamxxz, ass_Gamxyy, ass_Gamxyz, ass_Gamxzz, &
                        ass_Gamyxx, ass_Gamyxy, ass_Gamyxz, ass_Gamyyy, ass_Gamyyz, ass_Gamyzz, &
                        ass_Gamzxx, ass_Gamzxy, ass_Gamzxz, ass_Gamzyy, ass_Gamzyz, ass_Gamzzz)
! invert metric
  gupzz =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz - &
           gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz
  gupxx =   ( gyy * gzz - gyz * gyz ) / gupzz
  gupxy = - ( gxy * gzz - gyz * gxz ) / gupzz
  gupxz =   ( gxy * gyz - gyy * gxz ) / gupzz
  gupyy =   ( gxx * gzz - gxz * gxz ) / gupzz
  gupyz = - ( gxx * gyz - gxy * gxz ) / gupzz
  gupzz =   ( gxx * gyy - gxy * gxy ) / gupzz

  call kind2_connection(ex,                      &
                        gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
                        ass_Gamxxx, ass_Gamxxy, ass_Gamxxz, ass_Gamxyy, ass_Gamxyz, ass_Gamxzz, &
                        ass_Gamyxx, ass_Gamyxy, ass_Gamyxz, ass_Gamyyy, ass_Gamyyz, ass_Gamyzz, &
                        ass_Gamzxx, ass_Gamzxy, ass_Gamzxz, ass_Gamzyy, ass_Gamzyz, ass_Gamzzz, &
                        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz, &
                        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz, &
                        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz)                        

  call fdderivs(ex,dxx,gxxxx,gxxxy,gxxxz,gxxyy,gxxyz,gxxzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)
  call fdderivs(ex,gxy,gxyxx,gxyxy,gxyxz,gxyyy,gxyyz,gxyzz,X,Y,Z,ANTI,ANTI,SYM ,Symmetry,0)
  call fdderivs(ex,gxz,gxzxx,gxzxy,gxzxz,gxzyy,gxzyz,gxzzz,X,Y,Z,ANTI,SYM ,ANTI,Symmetry,0)
  call fdderivs(ex,dyy,gyyxx,gyyxy,gyyxz,gyyyy,gyyyz,gyyzz,X,Y,Z,SYM, SYM ,SYM ,Symmetry,0)
  call fdderivs(ex,gyz,gyzxx,gyzxy,gyzxz,gyzyy,gyzyz,gyzzz,X,Y,Z,SYM ,ANTI,ANTI,Symmetry,0)
  call fdderivs(ex,dzz,gzzxx,gzzxy,gzzxz,gzzyy,gzzyz,gzzzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)

  call adm_riemann(ex, &
               gxxxx,gxxxy,gxxxz,gxxyy,gxxyz,gxxzz, &
               gxyxx,gxyxy,gxyxz,gxyyy,gxyyz,gxyzz, &
               gxzxx,gxzxy,gxzxz,gxzyy,gxzyz,gxzzz, &
               gyyxx,gyyxy,gyyxz,gyyyy,gyyyz,gyyzz, &
               gyzxx,gyzxy,gyzxz,gyzyy,gyzyz,gyzzz, &
               gzzxx,gzzxy,gzzxz,gzzyy,gzzyz,gzzzz, &
               Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz, &
               Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz, &
               Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz, &
               ass_Gamxxx,ass_Gamxxy,ass_Gamxxz, ass_Gamxyy,ass_Gamxyz,ass_Gamxzz, &
               ass_Gamyxx,ass_Gamyxy,ass_Gamyxz, ass_Gamyyy,ass_Gamyyz,ass_Gamyzz, &
               ass_Gamzxx,ass_Gamzxy,ass_Gamzxz, ass_Gamzyy,ass_Gamzyz,ass_Gamzzz, &
               Rxyxy, Rxyxz, Rxyyz, Rxzxz, Rxzyz, Ryzyz)

  call adm_ricci(ex,                          &
                   gupxx , gupxy , gupxz , gupyy , gupyz , gupzz , &
                   Rxyxy, Rxyxz, Rxyyz, Rxzxz, Rxzyz, Ryzyz, &
                   Rxx ,   Rxy ,   Rxz ,   Ryy ,   Ryz ,   Rzz)

  return

  end subroutine adm_ricci_gamma
!----------------------------------------------------------------------------
  subroutine adm_ricci_gamma_ss(ex,crho,sigma,R,X, Y, Z,                      &
               drhodx, drhody, drhodz,                                         &
               dsigmadx,dsigmady,dsigmadz,                                     &
               dRdx,dRdy,dRdz,                                                 &
               drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                &
               dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,    &
               dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,                            &
               dxx    ,   gxy    ,   gxz    ,   dyy    ,   gyz    ,   dzz,&
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,&
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,&
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,&
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,&
               Symmetry,Lev,sst)
  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3), Symmetry,Lev,sst
  double precision,intent(in),dimension(ex(1))::crho
  double precision,intent(in),dimension(ex(2))::sigma
  double precision,intent(in),dimension(ex(3))::R
  real*8, intent(in ),dimension(ex(1),ex(2),ex(3)):: X,Y,Z
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::drhodx, drhody, drhodz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::dsigmadx,dsigmady,dsigmadz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::dRdx,dRdy,dRdz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: dxx,gxy,gxz,dyy,gyz,dzz
! when out, physical second kind of connection  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamzyy, Gamzyz, Gamzzz
! when out, physical Ricci tensor  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Rxx,Rxy,Rxz,Ryy,Ryz,Rzz

!~~~~~~> Other variables:

  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gyy,gzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxx,gxyx,gxzx,gyyx,gyzx,gzzx
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxy,gxyy,gxzy,gyyy,gyzy,gzzy
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxz,gxyz,gxzz,gyyz,gyzz,gzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: ass_Gamxxx,ass_Gamxxy,ass_Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: ass_Gamxyy,ass_Gamxyz,ass_Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: ass_Gamyxx,ass_Gamyxy,ass_Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: ass_Gamyyy,ass_Gamyyz,ass_Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: ass_Gamzxx,ass_Gamzxy,ass_Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: ass_Gamzyy,ass_Gamzyz,ass_Gamzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupyy,gupyz,gupzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxxx,gxxxy,gxxxz,gxxyy,gxxyz,gxxzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxyxx,gxyxy,gxyxz,gxyyy,gxyyz,gxyzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxzxx,gxzxy,gxzxz,gxzyy,gxzyz,gxzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyyxx,gyyxy,gyyxz,gyyyy,gyyyz,gyyzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyzxx,gyzxy,gyzxz,gyzyy,gyzyz,gyzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gzzxx,gzzxy,gzzxz,gzzyy,gzzyz,gzzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Rxyxy,  Rxyxz,  Rxyyz,  Rxzxz,  Rxzyz,  Ryzyz

  real*8, parameter :: ONE = 1.D0, TWO = 2.D0, FOUR = 4.D0
  real*8, parameter :: HALF = 0.5D0, F2o3 = 2.d0/3.d0, F3o2 = 1.5d0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0

  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE
 
  call fderivs_shc(ex,dxx,gxxx,gxxy,gxxz,crho,sigma,R, SYM, SYM,SYM,Symmetry,Lev,sst,          &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,gxy,gxyx,gxyy,gxyz,crho,sigma,R,ANTI,ANTI,SYM,Symmetry,Lev,sst,          &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,gxz,gxzx,gxzy,gxzz,crho,sigma,R,ANTI,SYM ,ANTI,Symmetry,Lev,sst,         &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,dyy,gyyx,gyyy,gyyz,crho,sigma,R, SYM, SYM,SYM,Symmetry,Lev,sst,          &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,gyz,gyzx,gyzy,gyzz,crho,sigma,R,SYM ,ANTI,ANTI,Symmetry,Lev,sst,         &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,dzz,gzzx,gzzy,gzzz,crho,sigma,R, SYM, SYM,SYM,Symmetry,Lev,sst,          &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)

  call kind1_connection(ex,                      &
                        gxxx,gxyx,gxzx,gyyx,gyzx,gzzx, &
                        gxxy,gxyy,gxzy,gyyy,gyzy,gzzy, &
                        gxxz,gxyz,gxzz,gyyz,gyzz,gzzz, &
                        ass_Gamxxx, ass_Gamxxy, ass_Gamxxz, ass_Gamxyy, ass_Gamxyz, ass_Gamxzz, &
                        ass_Gamyxx, ass_Gamyxy, ass_Gamyxz, ass_Gamyyy, ass_Gamyyz, ass_Gamyzz, &
                        ass_Gamzxx, ass_Gamzxy, ass_Gamzxz, ass_Gamzyy, ass_Gamzyz, ass_Gamzzz)
! invert metric
  gupzz =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz - &
           gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz
  gupxx =   ( gyy * gzz - gyz * gyz ) / gupzz
  gupxy = - ( gxy * gzz - gyz * gxz ) / gupzz
  gupxz =   ( gxy * gyz - gyy * gxz ) / gupzz
  gupyy =   ( gxx * gzz - gxz * gxz ) / gupzz
  gupyz = - ( gxx * gyz - gxy * gxz ) / gupzz
  gupzz =   ( gxx * gyy - gxy * gxy ) / gupzz

  call kind2_connection(ex,                      &
                        gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
                        ass_Gamxxx, ass_Gamxxy, ass_Gamxxz, ass_Gamxyy, ass_Gamxyz, ass_Gamxzz, &
                        ass_Gamyxx, ass_Gamyxy, ass_Gamyxz, ass_Gamyyy, ass_Gamyyz, ass_Gamyzz, &
                        ass_Gamzxx, ass_Gamzxy, ass_Gamzxz, ass_Gamzyy, ass_Gamzyz, ass_Gamzzz, &
                        Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz, &
                        Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz, &
                        Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz)                        

  call fdderivs_shc(ex,dxx,gxxxx,gxxxy,gxxxz,gxxyy,gxxyz,gxxzz,crho,sigma,R, SYM, SYM,SYM ,Symmetry,Lev,sst,  &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz)
  call fdderivs_shc(ex,dyy,gyyxx,gyyxy,gyyxz,gyyyy,gyyyz,gyyzz,crho,sigma,R, SYM, SYM,SYM ,Symmetry,Lev,sst,  &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz)
  call fdderivs_shc(ex,dzz,gzzxx,gzzxy,gzzxz,gzzyy,gzzyz,gzzzz,crho,sigma,R, SYM, SYM,SYM ,Symmetry,Lev,sst,  &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz)
  call fdderivs_shc(ex,gxy,gxyxx,gxyxy,gxyxz,gxyyy,gxyyz,gxyzz,crho,sigma,R,ANTI,ANTI,SYM ,Symmetry,Lev,sst,  &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz)
  call fdderivs_shc(ex,gxz,gxzxx,gxzxy,gxzxz,gxzyy,gxzyz,gxzzz,crho,sigma,R,ANTI,SYM ,ANTI,Symmetry,Lev,sst,  &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz)
  call fdderivs_shc(ex,gyz,gyzxx,gyzxy,gyzxz,gyzyy,gyzyz,gyzzz,crho,sigma,R,SYM ,ANTI,ANTI,Symmetry,Lev,sst,  &
                       drhodx, drhody, drhodz,                                                    &
                       dsigmadx,dsigmady,dsigmadz,                                                &
                       dRdx,dRdy,dRdz,                                                            &
                       drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                           &
                       dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,               &
                       dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz)

  call adm_riemann(ex, &
               gxxxx,gxxxy,gxxxz,gxxyy,gxxyz,gxxzz, &
               gxyxx,gxyxy,gxyxz,gxyyy,gxyyz,gxyzz, &
               gxzxx,gxzxy,gxzxz,gxzyy,gxzyz,gxzzz, &
               gyyxx,gyyxy,gyyxz,gyyyy,gyyyz,gyyzz, &
               gyzxx,gyzxy,gyzxz,gyzyy,gyzyz,gyzzz, &
               gzzxx,gzzxy,gzzxz,gzzyy,gzzyz,gzzzz, &
               Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz, &
               Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz, &
               Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz, &
               ass_Gamxxx,ass_Gamxxy,ass_Gamxxz, ass_Gamxyy,ass_Gamxyz,ass_Gamxzz, &
               ass_Gamyxx,ass_Gamyxy,ass_Gamyxz, ass_Gamyyy,ass_Gamyyz,ass_Gamyzz, &
               ass_Gamzxx,ass_Gamzxy,ass_Gamzxz, ass_Gamzyy,ass_Gamzyz,ass_Gamzzz, &
               Rxyxy, Rxyxz, Rxyyz, Rxzxz, Rxzyz, Ryzyz)

  call adm_ricci(ex,                          &
                   gupxx , gupxy , gupxz , gupyy , gupyz , gupzz , &
                   Rxyxy, Rxyxz, Rxyyz, Rxzxz, Rxzyz, Ryzyz, &
                   Rxx ,   Rxy ,   Rxz ,   Ryy ,   Ryz ,   Rzz)

  return

  end subroutine adm_ricci_gamma_ss
