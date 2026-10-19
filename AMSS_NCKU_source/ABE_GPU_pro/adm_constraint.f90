
!-------------------------------------------------------------------------------!
! computed constraint for ADM formalism                                         !
!-------------------------------------------------------------------------------!
  subroutine constraint_adm(ex, X, Y, Z,&
               dxx,gxy,gxz,dyy,gyz,dzz, &
               Kxx,Kxy,Kxz,Kyy,Kyz,Kzz, &
               Lap,Sfx,Sfy,Sfz,rho,Sx,Sy,Sz,&
               ham_Res, movx_Res, movy_Res, movz_Res, &
               Symmetry)

  implicit none
!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3),symmetry
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: dxx,gxy,gxz,dyy,gyz,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Kxx,Kxy,Kxz,Kyy,Kyz,Kzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Lap,Sfx,Sfy,Sfz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: rho,Sx,Sy,Sz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: ham_Res, movx_Res, movy_Res, movz_Res
!~~~~~~> Other variables:
!  inverse metric
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupyy,gupyz,gupzz
! first order derivative of metric, @_k g_ij
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxx,gxyx,gxzx
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyyx,gyzx,gzzx
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxy,gxyy,gxzy
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyyy,gyzy,gzzy
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxz,gxyz,gxzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyyz,gyzz,gzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gyy,gzz,trK,fx,fy,fz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Rxx,Rxy,Rxz,Ryy,Ryz,Rzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamzyy, Gamzyz, Gamzzz

  integer, parameter :: NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2
  real*8, parameter :: ZERO = 0.D0, HALF = 0.5d0, ONE = 1.d0, TWO = 2.d0, FOUR = 4.d0
  real*8, parameter :: F2o3 = 2.d0/3.d0, F8 = 8.d0, F16 = 1.6d1, SIX = 6.d0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  real*8            :: PI

  call adm_ricci_gamma(ex, X, Y, Z,                        &
               dxx    ,   gxy    ,   gxz    ,   dyy    ,   gyz    ,   dzz,&
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,&
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,&
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,&
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,&
               Symmetry)

  PI = dacos(-ONE)

  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE
! invert metric
  gupzz =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz - &
           gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz
  gupxx =   ( gyy * gzz - gyz * gyz ) / gupzz
  gupxy = - ( gxy * gzz - gyz * gxz ) / gupzz
  gupxz =   ( gxy * gyz - gyy * gxz ) / gupzz
  gupyy =   ( gxx * gzz - gxz * gxz ) / gupzz
  gupyz = - ( gxx * gyz - gxy * gxz ) / gupzz
  gupzz =   ( gxx * gyy - gxy * gxy ) / gupzz

  trK =          gupxx * Kxx + gupyy * Kyy + gupzz * Kzz &
        + TWO * (gupxy * Kxy + gupxz * Kxz + gupyz * Kyz)

! ham_Res = trR + K^2 - K_ij * K^ij - 16 * PI * rho
  ham_Res =   gupxx * Rxx + gupyy * Ryy + gupzz * Rzz + &
        TWO* ( gupxy * Rxy + gupxz * Rxz + gupyz * Ryz )

  ham_Res = ham_Res + trK * trK -(&
       gupxx * ( &
       gupxx * Kxx * Kxx + gupyy * Kxy * Kxy + gupzz * Kxz * Kxz + &
       TWO * (gupxy * Kxx * Kxy + gupxz * Kxx * Kxz + gupyz * Kxy * Kxz) ) + &
       gupyy * ( &
       gupxx * Kxy * Kxy + gupyy * Kyy * Kyy + gupzz * Kyz * Kyz + &
       TWO * (gupxy * Kxy * Kyy + gupxz * Kxy * Kyz + gupyz * Kyy * Kyz) ) + &
       gupzz * ( &
       gupxx * Kxz * Kxz + gupyy * Kyz * Kyz + gupzz * Kzz * Kzz + &
       TWO * (gupxy * Kxz * Kyz + gupxz * Kxz * Kzz + gupyz * Kyz * Kzz) ) + &
       TWO * ( &
       gupxy * ( &
       gupxx * Kxx * Kxy + gupyy * Kxy * Kyy + gupzz * Kxz * Kyz + &
       gupxy * (Kxx * Kyy + Kxy * Kxy) + &
       gupxz * (Kxx * Kyz + Kxz * Kxy) + &
       gupyz * (Kxy * Kyz + Kxz * Kyy) ) + &
       gupxz * ( &
       gupxx * Kxx * Kxz + gupyy * Kxy * Kyz + gupzz * Kxz * Kzz + &
       gupxy * (Kxx * Kyz + Kxy * Kxz) + &
       gupxz * (Kxx * Kzz + Kxz * Kxz) + &
       gupyz * (Kxy * Kzz + Kxz * Kyz) ) + &
       gupyz * ( &
       gupxx * Kxy * Kxz + gupyy * Kyy * Kyz + gupzz * Kyz * Kzz + &
       gupxy * (Kxy * Kyz + Kyy * Kxz) + &
       gupxz * (Kxy * Kzz + Kyz * Kxz) + &
       gupyz * (Kyy * Kzz + Kyz * Kyz) ) ))- F16 * PI * rho

! mov_Res_j = gupkj*D_k K_ij - d_j trK - 8 PI s_j where D respect to physical metric
! store D_i K_jk
  call fderivs(ex,Kxx,gxxx,gxxy,gxxz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)
  call fderivs(ex,Kxy,gxyx,gxyy,gxyz,X,Y,Z,ANTI,ANTI,SYM ,Symmetry,0)
  call fderivs(ex,Kxz,gxzx,gxzy,gxzz,X,Y,Z,ANTI,SYM ,ANTI,Symmetry,0)
  call fderivs(ex,Kyy,gyyx,gyyy,gyyz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)
  call fderivs(ex,Kyz,gyzx,gyzy,gyzz,X,Y,Z,SYM ,ANTI,ANTI,Symmetry,0)
  call fderivs(ex,Kzz,gzzx,gzzy,gzzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)

  gxxx = gxxx - (  Gamxxx * Kxx + Gamyxx * Kxy + Gamzxx * Kxz &
                 + Gamxxx * Kxx + Gamyxx * Kxy + Gamzxx * Kxz)
  gxyx = gxyx - (  Gamxxy * Kxx + Gamyxy * Kxy + Gamzxy * Kxz &
                 + Gamxxx * Kxy + Gamyxx * Kyy + Gamzxx * Kyz)
  gxzx = gxzx - (  Gamxxz * Kxx + Gamyxz * Kxy + Gamzxz * Kxz &
                 + Gamxxx * Kxz + Gamyxx * Kyz + Gamzxx * Kzz)
  gyyx = gyyx - (  Gamxxy * Kxy + Gamyxy * Kyy + Gamzxy * Kyz &
                 + Gamxxy * Kxy + Gamyxy * Kyy + Gamzxy * Kyz)
  gyzx = gyzx - (  Gamxxz * Kxy + Gamyxz * Kyy + Gamzxz * Kyz &
                 + Gamxxy * Kxz + Gamyxy * Kyz + Gamzxy * Kzz)
  gzzx = gzzx - (  Gamxxz * Kxz + Gamyxz * Kyz + Gamzxz * Kzz &
                 + Gamxxz * Kxz + Gamyxz * Kyz + Gamzxz * Kzz)
  gxxy = gxxy - (  Gamxxy * Kxx + Gamyxy * Kxy + Gamzxy * Kxz &
                 + Gamxxy * Kxx + Gamyxy * Kxy + Gamzxy * Kxz)
  gxyy = gxyy - (  Gamxyy * Kxx + Gamyyy * Kxy + Gamzyy * Kxz &
                 + Gamxxy * Kxy + Gamyxy * Kyy + Gamzxy * Kyz)
  gxzy = gxzy - (  Gamxyz * Kxx + Gamyyz * Kxy + Gamzyz * Kxz &
                 + Gamxxy * Kxz + Gamyxy * Kyz + Gamzxy * Kzz)
  gyyy = gyyy - (  Gamxyy * Kxy + Gamyyy * Kyy + Gamzyy * Kyz &
                 + Gamxyy * Kxy + Gamyyy * Kyy + Gamzyy * Kyz)
  gyzy = gyzy - (  Gamxyz * Kxy + Gamyyz * Kyy + Gamzyz * Kyz &
                 + Gamxyy * Kxz + Gamyyy * Kyz + Gamzyy * Kzz)
  gzzy = gzzy - (  Gamxyz * Kxz + Gamyyz * Kyz + Gamzyz * Kzz &
                 + Gamxyz * Kxz + Gamyyz * Kyz + Gamzyz * Kzz)
  gxxz = gxxz - (  Gamxxz * Kxx + Gamyxz * Kxy + Gamzxz * Kxz &
                 + Gamxxz * Kxx + Gamyxz * Kxy + Gamzxz * Kxz)
  gxyz = gxyz - (  Gamxyz * Kxx + Gamyyz * Kxy + Gamzyz * Kxz &
                 + Gamxxz * Kxy + Gamyxz * Kyy + Gamzxz * Kyz)
  gxzz = gxzz - (  Gamxzz * Kxx + Gamyzz * Kxy + Gamzzz * Kxz &
                 + Gamxxz * Kxz + Gamyxz * Kyz + Gamzxz * Kzz)
  gyyz = gyyz - (  Gamxyz * Kxy + Gamyyz * Kyy + Gamzyz * Kyz &
                 + Gamxyz * Kxy + Gamyyz * Kyy + Gamzyz * Kyz)
  gyzz = gyzz - (  Gamxzz * Kxy + Gamyzz * Kyy + Gamzzz * Kyz &
                 + Gamxyz * Kxz + Gamyyz * Kyz + Gamzyz * Kzz)
  gzzz = gzzz - (  Gamxzz * Kxz + Gamyzz * Kyz + Gamzzz * Kzz &
                 + Gamxzz * Kxz + Gamyzz * Kyz + Gamzzz * Kzz)
movx_Res = gupxx*gxxx + gupyy*gxyy + gupzz*gxzz &
          +gupxy*gxyx + gupxz*gxzx + gupyz*gxzy &
          +gupxy*gxxy + gupxz*gxxz + gupyz*gxyz
movy_Res = gupxx*gxyx + gupyy*gyyy + gupzz*gyzz &
          +gupxy*gyyx + gupxz*gyzx + gupyz*gyzy &
          +gupxy*gxyy + gupxz*gxyz + gupyz*gyyz
movz_Res = gupxx*gxzx + gupyy*gyzy + gupzz*gzzz &
          +gupxy*gyzx + gupxz*gzzx + gupyz*gzzy &
          +gupxy*gxzy + gupxz*gxzz + gupyz*gyzz

  call fderivs(ex,trK,fx,fy,fz,X,Y,Z,SYM,SYM,SYM,Symmetry,0)

movx_Res = movx_Res - fx - F8*PI*sx
movy_Res = movy_Res - fy - F8*PI*sy
movz_Res = movz_Res - fz - F8*PI*sz

  return

  end subroutine constraint_adm
!-------------------------------------------------------------------------------!
! computed constraint for ADM formalism for shell                              !
!-------------------------------------------------------------------------------!
  subroutine constraint_adm_ss(ex,crho,sigma,R, X, Y, Z,&
               drhodx, drhody, drhodz,                                         &
               dsigmadx,dsigmady,dsigmadz,                                     &
               dRdx,dRdy,dRdz,                                                 &
               drhodxx,drhodxy,drhodxz,drhodyy,drhodyz,drhodzz,                &
               dsigmadxx,dsigmadxy,dsigmadxz,dsigmadyy,dsigmadyz,dsigmadzz,    &
               dRdxx,dRdxy,dRdxz,dRdyy,dRdyz,dRdzz,                            &
               dxx,gxy,gxz,dyy,gyz,dzz, &
               Kxx,Kxy,Kxz,Kyy,Kyz,Kzz, &
               Lap,Sfx,Sfy,Sfz,rho,Sx,Sy,Sz,&
               Gamxxx, Gamxxy, Gamxxz,Gamxyy, Gamxyz, Gamxzz, &
               Gamyxx, Gamyxy, Gamyxz,Gamyyy, Gamyyz, Gamyzz, &
               Gamzxx, Gamzxy, Gamzxz,Gamzyy, Gamzyz, Gamzzz, &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz, &
               ham_Res, movx_Res, movy_Res, movz_Res, &
               Symmetry,Lev,sst)

  implicit none
!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3),symmetry,Lev,sst
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
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Kxx,Kxy,Kxz,Kyy,Kyz,Kzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Lap,Sfx,Sfy,Sfz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: rho,Sx,Sy,Sz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Rxx,Rxy,Rxz,Ryy,Ryz,Rzz
! second kind of Christofel symble Gamma^i_jk respect to physical metric
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamzyy, Gamzyz, Gamzzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: ham_Res, movx_Res, movy_Res, movz_Res
!~~~~~~> Other variables:
!  inverse metric
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupyy,gupyz,gupzz
! first order derivative of metric, @_k g_ij
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxx,gxyx,gxzx
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyyx,gyzx,gzzx
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxy,gxyy,gxzy
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyyy,gyzy,gzzy
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxz,gxyz,gxzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gyyz,gyzz,gzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gyy,gzz,trK,fx,fy,fz

  integer, parameter :: NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2
  real*8, parameter :: ZERO = 0.D0, HALF = 0.5d0, ONE = 1.d0, TWO = 2.d0, FOUR = 4.d0
  real*8, parameter :: F2o3 = 2.d0/3.d0, F8 = 8.d0, F16 = 1.6d1, SIX = 6.d0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  real*8            :: PI

  call adm_ricci_gamma_ss(ex,crho,sigma,R,X, Y, Z,                      &
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

  PI = dacos(-ONE)

  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE
! invert metric
  gupzz =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz - &
           gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz
  gupxx =   ( gyy * gzz - gyz * gyz ) / gupzz
  gupxy = - ( gxy * gzz - gyz * gxz ) / gupzz
  gupxz =   ( gxy * gyz - gyy * gxz ) / gupzz
  gupyy =   ( gxx * gzz - gxz * gxz ) / gupzz
  gupyz = - ( gxx * gyz - gxy * gxz ) / gupzz
  gupzz =   ( gxx * gyy - gxy * gxy ) / gupzz

  trK =          gupxx * Kxx + gupyy * Kyy + gupzz * Kzz &
        + TWO * (gupxy * Kxy + gupxz * Kxz + gupyz * Kyz)

! ham_Res = trR + K^2 - K_ij * K^ij - 16 * PI * rho
  ham_Res =   gupxx * Rxx + gupyy * Ryy + gupzz * Rzz + &
        TWO* ( gupxy * Rxy + gupxz * Rxz + gupyz * Ryz )

  ham_Res = ham_Res + trK * trK -(&
       gupxx * ( &
       gupxx * Kxx * Kxx + gupyy * Kxy * Kxy + gupzz * Kxz * Kxz + &
       TWO * (gupxy * Kxx * Kxy + gupxz * Kxx * Kxz + gupyz * Kxy * Kxz) ) + &
       gupyy * ( &
       gupxx * Kxy * Kxy + gupyy * Kyy * Kyy + gupzz * Kyz * Kyz + &
       TWO * (gupxy * Kxy * Kyy + gupxz * Kxy * Kyz + gupyz * Kyy * Kyz) ) + &
       gupzz * ( &
       gupxx * Kxz * Kxz + gupyy * Kyz * Kyz + gupzz * Kzz * Kzz + &
       TWO * (gupxy * Kxz * Kyz + gupxz * Kxz * Kzz + gupyz * Kyz * Kzz) ) + &
       TWO * ( &
       gupxy * ( &
       gupxx * Kxx * Kxy + gupyy * Kxy * Kyy + gupzz * Kxz * Kyz + &
       gupxy * (Kxx * Kyy + Kxy * Kxy) + &
       gupxz * (Kxx * Kyz + Kxz * Kxy) + &
       gupyz * (Kxy * Kyz + Kxz * Kyy) ) + &
       gupxz * ( &
       gupxx * Kxx * Kxz + gupyy * Kxy * Kyz + gupzz * Kxz * Kzz + &
       gupxy * (Kxx * Kyz + Kxy * Kxz) + &
       gupxz * (Kxx * Kzz + Kxz * Kxz) + &
       gupyz * (Kxy * Kzz + Kxz * Kyz) ) + &
       gupyz * ( &
       gupxx * Kxy * Kxz + gupyy * Kyy * Kyz + gupzz * Kyz * Kzz + &
       gupxy * (Kxy * Kyz + Kyy * Kxz) + &
       gupxz * (Kxy * Kzz + Kyz * Kxz) + &
       gupyz * (Kyy * Kzz + Kyz * Kyz) ) ))- F16 * PI * rho

! mov_Res_j = gupkj*D_k K_ij - d_j trK - 8 PI s_j where D respect to physical metric
! store D_i K_jk
  call fderivs_shc(ex,Kxx,gxxx,gxxy,gxxz,crho,sigma,R, SYM, SYM,SYM,Symmetry,Lev,sst,          &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,Kxy,gxyx,gxyy,gxyz,crho,sigma,R,ANTI,ANTI,SYM,Symmetry,Lev,sst,          &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,Kxz,gxzx,gxzy,gxzz,crho,sigma,R,ANTI,SYM ,ANTI,Symmetry,Lev,sst,         &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,Kyy,gyyx,gyyy,gyyz,crho,sigma,R, SYM, SYM,SYM,Symmetry,Lev,sst,          &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,Kyz,gyzx,gyzy,gyzz,crho,sigma,R,SYM ,ANTI,ANTI,Symmetry,Lev,sst,         &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)
  call fderivs_shc(ex,Kzz,gzzx,gzzy,gzzz,crho,sigma,R, SYM, SYM,SYM,Symmetry,Lev,sst,          &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)

  gxxx = gxxx - (  Gamxxx * Kxx + Gamyxx * Kxy + Gamzxx * Kxz &
                 + Gamxxx * Kxx + Gamyxx * Kxy + Gamzxx * Kxz)
  gxyx = gxyx - (  Gamxxy * Kxx + Gamyxy * Kxy + Gamzxy * Kxz &
                 + Gamxxx * Kxy + Gamyxx * Kyy + Gamzxx * Kyz)
  gxzx = gxzx - (  Gamxxz * Kxx + Gamyxz * Kxy + Gamzxz * Kxz &
                 + Gamxxx * Kxz + Gamyxx * Kyz + Gamzxx * Kzz)
  gyyx = gyyx - (  Gamxxy * Kxy + Gamyxy * Kyy + Gamzxy * Kyz &
                 + Gamxxy * Kxy + Gamyxy * Kyy + Gamzxy * Kyz)
  gyzx = gyzx - (  Gamxxz * Kxy + Gamyxz * Kyy + Gamzxz * Kyz &
                 + Gamxxy * Kxz + Gamyxy * Kyz + Gamzxy * Kzz)
  gzzx = gzzx - (  Gamxxz * Kxz + Gamyxz * Kyz + Gamzxz * Kzz &
                 + Gamxxz * Kxz + Gamyxz * Kyz + Gamzxz * Kzz)
  gxxy = gxxy - (  Gamxxy * Kxx + Gamyxy * Kxy + Gamzxy * Kxz &
                 + Gamxxy * Kxx + Gamyxy * Kxy + Gamzxy * Kxz)
  gxyy = gxyy - (  Gamxyy * Kxx + Gamyyy * Kxy + Gamzyy * Kxz &
                 + Gamxxy * Kxy + Gamyxy * Kyy + Gamzxy * Kyz)
  gxzy = gxzy - (  Gamxyz * Kxx + Gamyyz * Kxy + Gamzyz * Kxz &
                 + Gamxxy * Kxz + Gamyxy * Kyz + Gamzxy * Kzz)
  gyyy = gyyy - (  Gamxyy * Kxy + Gamyyy * Kyy + Gamzyy * Kyz &
                 + Gamxyy * Kxy + Gamyyy * Kyy + Gamzyy * Kyz)
  gyzy = gyzy - (  Gamxyz * Kxy + Gamyyz * Kyy + Gamzyz * Kyz &
                 + Gamxyy * Kxz + Gamyyy * Kyz + Gamzyy * Kzz)
  gzzy = gzzy - (  Gamxyz * Kxz + Gamyyz * Kyz + Gamzyz * Kzz &
                 + Gamxyz * Kxz + Gamyyz * Kyz + Gamzyz * Kzz)
  gxxz = gxxz - (  Gamxxz * Kxx + Gamyxz * Kxy + Gamzxz * Kxz &
                 + Gamxxz * Kxx + Gamyxz * Kxy + Gamzxz * Kxz)
  gxyz = gxyz - (  Gamxyz * Kxx + Gamyyz * Kxy + Gamzyz * Kxz &
                 + Gamxxz * Kxy + Gamyxz * Kyy + Gamzxz * Kyz)
  gxzz = gxzz - (  Gamxzz * Kxx + Gamyzz * Kxy + Gamzzz * Kxz &
                 + Gamxxz * Kxz + Gamyxz * Kyz + Gamzxz * Kzz)
  gyyz = gyyz - (  Gamxyz * Kxy + Gamyyz * Kyy + Gamzyz * Kyz &
                 + Gamxyz * Kxy + Gamyyz * Kyy + Gamzyz * Kyz)
  gyzz = gyzz - (  Gamxzz * Kxy + Gamyzz * Kyy + Gamzzz * Kyz &
                 + Gamxyz * Kxz + Gamyyz * Kyz + Gamzyz * Kzz)
  gzzz = gzzz - (  Gamxzz * Kxz + Gamyzz * Kyz + Gamzzz * Kzz &
                 + Gamxzz * Kxz + Gamyzz * Kyz + Gamzzz * Kzz)
movx_Res = gupxx*gxxx + gupyy*gxyy + gupzz*gxzz &
          +gupxy*gxyx + gupxz*gxzx + gupyz*gxzy &
          +gupxy*gxxy + gupxz*gxxz + gupyz*gxyz
movy_Res = gupxx*gxyx + gupyy*gyyy + gupzz*gyzz &
          +gupxy*gyyx + gupxz*gyzx + gupyz*gyzy &
          +gupxy*gxyy + gupxz*gxyz + gupyz*gyyz
movz_Res = gupxx*gxzx + gupyy*gyzy + gupzz*gzzz &
          +gupxy*gyzx + gupxz*gzzx + gupyz*gzzy &
          +gupxy*gxzy + gupxz*gxzz + gupyz*gyzz

  call fderivs_shc(ex,trK,fx,fy,fz,crho,sigma,R, SYM, SYM,SYM,Symmetry,Lev,sst,                &
                       drhodx, drhody, drhodz,                                                 &
                       dsigmadx,dsigmady,dsigmadz,                                             &
                       dRdx,dRdy,dRdz)

movx_Res = movx_Res - fx - F8*PI*sx
movy_Res = movy_Res - fy - F8*PI*sy
movz_Res = movz_Res - fz - F8*PI*sz

  return

  end subroutine constraint_adm_ss
