

#include "macrodef.fh"

! Per-order halo / symmetry-bd ord macros for compute_rhs_bssn_opt.
! Stencil radius = ghost_width - 1, so the lower fh bound is 1-radius =
! 2 - ghost_width, and symmetry_bd_inner needs ord = ghost_width - 1.
#if (ghost_width == 2)
#define OPT_HALO_LO 0
#define OPT_SYMORD 1
#elif (ghost_width == 3)
#define OPT_HALO_LO -1
#define OPT_SYMORD 2
#elif (ghost_width == 4)
#define OPT_HALO_LO -2
#define OPT_SYMORD 3
#elif (ghost_width == 5)
#define OPT_HALO_LO -3
#define OPT_SYMORD 4
#endif

  function compute_rhs_bssn(ex, T,X, Y, Z,                                     &
               chi    ,   trK    ,                                             &
               dxx    ,   gxy    ,   gxz    ,   dyy    ,   gyz    ,   dzz,     &
               Axx    ,   Axy    ,   Axz    ,   Ayy    ,   Ayz    ,   Azz,     &
               Gamx   ,  Gamy    ,  Gamz    ,                                  &
               Lap    ,  betax   ,  betay   ,  betaz   ,                       &
               dtSfx  ,  dtSfy   ,  dtSfz   ,                                  &
               chi_rhs,   trK_rhs,                                             &
               gxx_rhs,   gxy_rhs,   gxz_rhs,   gyy_rhs,   gyz_rhs,   gzz_rhs, &
               Axx_rhs,   Axy_rhs,   Axz_rhs,   Ayy_rhs,   Ayz_rhs,   Azz_rhs, &
               Gamx_rhs,  Gamy_rhs,  Gamz_rhs,                                 &
               Lap_rhs,  betax_rhs,  betay_rhs,  betaz_rhs,                    &
               dtSfx_rhs,  dtSfy_rhs,  dtSfz_rhs,                              &
               rho,Sx,Sy,Sz,Sxx,Sxy,Sxz,Syy,Syz,Szz,                           &
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,                      &
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,                      &
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,                      &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,                                        &
               ham_Res, movx_Res, movy_Res, movz_Res,                          &
                        Gmx_Res, Gmy_Res, Gmz_Res,                             &
               Symmetry,Lev,eps,co)  result(gont)
! Wrapper function to check consistency between opt and cor versions
  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3), Symmetry,Lev,co
  real*8, intent(in ):: T
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: chi,dxx,dyy,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: trK
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: gxy,gxz,gyz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamx,Gamy,Gamz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Lap, betax, betay, betaz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: dtSfx,  dtSfy,  dtSfz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: chi_rhs,trK_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gxx_rhs,gxy_rhs,gxz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gyy_rhs,gyz_rhs,gzz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Axx_rhs,Axy_rhs,Axz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Ayy_rhs,Ayz_rhs,Azz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamx_rhs,Gamy_rhs,Gamz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Lap_rhs, betax_rhs, betay_rhs, betaz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: dtSfx_rhs,dtSfy_rhs,dtSfz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: rho,Sx,Sy,Sz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Sxx,Sxy,Sxz,Syy,Syz,Szz
! when out, physical second kind of connection  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamzyy, Gamzyz, Gamzzz
! when out, physical Ricci tensor  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Rxx,Rxy,Rxz,Ryy,Ryz,Rzz
  real*8,intent(in) :: eps
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: ham_Res, movx_Res, movy_Res, movz_Res
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Gmx_Res, Gmy_Res, Gmz_Res
!  gont = 0: success; gont = 1: something wrong
  integer::gont

!~~~~~~> Temporary arrays for cor version:

  ! real*8, dimension(ex(1),ex(2),ex(3)) :: chi_cor,dxx_cor,dyy_cor,dzz_cor
  ! real*8, dimension(ex(1),ex(2),ex(3)) :: Lap_cor,betax_cor,betay_cor,betaz_cor
  ! real*8, dimension(ex(1),ex(2),ex(3)) :: chi_rhs_cor,trK_rhs_cor
  ! real*8, dimension(ex(1),ex(2),ex(3)) :: gxx_rhs_cor,gxy_rhs_cor,gxz_rhs_cor
  ! real*8, dimension(ex(1),ex(2),ex(3)) :: gyy_rhs_cor,gyz_rhs_cor,gzz_rhs_cor
  ! real*8, dimension(ex(1),ex(2),ex(3)) :: Axx_rhs_cor,Axy_rhs_cor,Axz_rhs_cor
  ! real*8, dimension(ex(1),ex(2),ex(3)) :: Ayy_rhs_cor,Ayz_rhs_cor,Azz_rhs_cor
  ! real*8, dimension(ex(1),ex(2),ex(3)) :: Gamx_rhs_cor,Gamy_rhs_cor,Gamz_rhs_cor
  ! real*8, dimension(ex(1),ex(2),ex(3)) :: Lap_rhs_cor, betax_rhs_cor, betay_rhs_cor, betaz_rhs_cor
  ! real*8, dimension(ex(1),ex(2),ex(3)) :: dtSfx_rhs_cor,dtSfy_rhs_cor,dtSfz_rhs_cor
  ! real*8, dimension(ex(1),ex(2),ex(3)) :: Gamxxx_cor, Gamxxy_cor, Gamxxz_cor
  ! real*8, dimension(ex(1),ex(2),ex(3)) :: Gamxyy_cor, Gamxyz_cor, Gamxzz_cor
  ! real*8, dimension(ex(1),ex(2),ex(3)) :: Gamyxx_cor, Gamyxy_cor, Gamyxz_cor
  ! real*8, dimension(ex(1),ex(2),ex(3)) :: Gamyyy_cor, Gamyyz_cor, Gamyzz_cor
  ! real*8, dimension(ex(1),ex(2),ex(3)) :: Gamzxx_cor, Gamzxy_cor, Gamzxz_cor
  ! real*8, dimension(ex(1),ex(2),ex(3)) :: Gamzyy_cor, Gamzyz_cor, Gamzzz_cor
  ! real*8, dimension(ex(1),ex(2),ex(3)) :: Rxx_cor,Rxy_cor,Rxz_cor,Ryy_cor,Ryz_cor,Rzz_cor
  ! real*8, dimension(ex(1),ex(2),ex(3)) :: ham_Res_cor, movx_Res_cor, movy_Res_cor, movz_Res_cor
  ! real*8, dimension(ex(1),ex(2),ex(3)) :: Gmx_Res_cor, Gmy_Res_cor, Gmz_Res_cor
  
!~~~~~~> Variables for comparison:
  integer :: gont_opt, gont_cor
  integer :: i, j, k, diff_count
  real*8 :: max_diff, tol, curr_diff
  
  tol = 0.0d0  ! Exact comparison: any difference triggers warning
  
!~~~~~~> Copy inout parameters for cor version:
  ! chi_cor = chi
  ! dxx_cor = dxx
  ! dyy_cor = dyy
  ! dzz_cor = dzz
  ! Lap_cor = Lap
  ! betax_cor = betax
  ! betay_cor = betay
  ! betaz_cor = betaz
  ! ham_Res_cor = ham_Res
  ! movx_Res_cor = movx_Res
  ! movy_Res_cor = movy_Res
  ! movz_Res_cor = movz_Res
  ! Gmx_Res_cor = Gmx_Res
  ! Gmy_Res_cor = Gmy_Res
  ! Gmz_Res_cor = Gmz_Res

!~~~~~~> Dispatch by ghost_width. compute_rhs_bssn_opt now supports all four
! ghost_widths (2/3/4/5) via OPT_HALO_LO/OPT_SYMORD macros plus per-order
! _inner derivative implementations in diff_new.f90. The opt path uses
! symmetry_bd_inner from fmisc.f90, which is only defined under #ifdef Cell;
! Vertex grids fall back to the legacy path.
#if defined(Cell) && (ghost_width == 2 || ghost_width == 3 || ghost_width == 4 || ghost_width == 5)
  gont_opt = compute_rhs_bssn_opt(ex, T, X, Y, Z, &
               chi, trK, &
               dxx, gxy, gxz, dyy, gyz, dzz, &
               Axx, Axy, Axz, Ayy, Ayz, Azz, &
               Gamx, Gamy, Gamz, &
               Lap, betax, betay, betaz, &
               dtSfx, dtSfy, dtSfz, &
               chi_rhs, trK_rhs, &
               gxx_rhs, gxy_rhs, gxz_rhs, gyy_rhs, gyz_rhs, gzz_rhs, &
               Axx_rhs, Axy_rhs, Axz_rhs, Ayy_rhs, Ayz_rhs, Azz_rhs, &
               Gamx_rhs, Gamy_rhs, Gamz_rhs, &
               Lap_rhs, betax_rhs, betay_rhs, betaz_rhs, &
               dtSfx_rhs, dtSfy_rhs, dtSfz_rhs, &
               rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz, &
               Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz, &
               Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz, &
               Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz, &
               Rxx, Rxy, Rxz, Ryy, Ryz, Rzz, &
               ham_Res, movx_Res, movy_Res, movz_Res, &
               Gmx_Res, Gmy_Res, Gmz_Res, &
               Symmetry, Lev, eps, co)
  gont = gont_opt
#else
  gont = compute_rhs_bssn_legacy(ex, T, X, Y, Z, &
               chi, trK, &
               dxx, gxy, gxz, dyy, gyz, dzz, &
               Axx, Axy, Axz, Ayy, Ayz, Azz, &
               Gamx, Gamy, Gamz, &
               Lap, betax, betay, betaz, &
               dtSfx, dtSfy, dtSfz, &
               chi_rhs, trK_rhs, &
               gxx_rhs, gxy_rhs, gxz_rhs, gyy_rhs, gyz_rhs, gzz_rhs, &
               Axx_rhs, Axy_rhs, Axz_rhs, Ayy_rhs, Ayz_rhs, Azz_rhs, &
               Gamx_rhs, Gamy_rhs, Gamz_rhs, &
               Lap_rhs, betax_rhs, betay_rhs, betaz_rhs, &
               dtSfx_rhs, dtSfy_rhs, dtSfz_rhs, &
               rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz, &
               Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz, &
               Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz, &
               Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz, &
               Rxx, Rxy, Rxz, Ryy, Ryz, Rzz, &
               ham_Res, movx_Res, movy_Res, movz_Res, &
               Gmx_Res, Gmy_Res, Gmz_Res, &
               Symmetry, Lev, eps, co)
#endif

  return

! !~~~~~~> Call cor version:
!   gont_cor = compute_rhs_bssn_cor(ex, T, X, Y, Z, &
!                chi_cor, trK, &
!                dxx_cor, gxy, gxz, dyy_cor, gyz, dzz_cor, &
!                Axx, Axy, Axz, Ayy, Ayz, Azz, &
!                Gamx, Gamy, Gamz, &
!                Lap_cor, betax_cor, betay_cor, betaz_cor, &
!                dtSfx, dtSfy, dtSfz, &
!                chi_rhs_cor, trK_rhs_cor, &
!                gxx_rhs_cor, gxy_rhs_cor, gxz_rhs_cor, gyy_rhs_cor, gyz_rhs_cor, gzz_rhs_cor, &
!                Axx_rhs_cor, Axy_rhs_cor, Axz_rhs_cor, Ayy_rhs_cor, Ayz_rhs_cor, Azz_rhs_cor, &
!                Gamx_rhs_cor, Gamy_rhs_cor, Gamz_rhs_cor, &
!                Lap_rhs_cor, betax_rhs_cor, betay_rhs_cor, betaz_rhs_cor, &
!                dtSfx_rhs_cor, dtSfy_rhs_cor, dtSfz_rhs_cor, &
!                rho, Sx, Sy, Sz, Sxx, Sxy, Sxz, Syy, Syz, Szz, &
!                Gamxxx_cor, Gamxxy_cor, Gamxxz_cor, Gamxyy_cor, Gamxyz_cor, Gamxzz_cor, &
!                Gamyxx_cor, Gamyxy_cor, Gamyxz_cor, Gamyyy_cor, Gamyyz_cor, Gamyzz_cor, &
!                Gamzxx_cor, Gamzxy_cor, Gamzxz_cor, Gamzyy_cor, Gamzyz_cor, Gamzzz_cor, &
!                Rxx_cor, Rxy_cor, Rxz_cor, Ryy_cor, Ryz_cor, Rzz_cor, &
!                ham_Res_cor, movx_Res_cor, movy_Res_cor, movz_Res_cor, &
!                Gmx_Res_cor, Gmy_Res_cor, Gmz_Res_cor, &
!                Symmetry, Lev, eps, co)

! !~~~~~~> Compare results:
!   diff_count = 0
!   max_diff = 0.0d0

!   ! Check all RHS outputs
!   do k = 1, ex(3)
!     do j = 1, ex(2)
!       do i = 1, ex(1)
!         curr_diff = abs(chi_rhs(i,j,k) - chi_rhs_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(trK_rhs(i,j,k) - trK_rhs_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(gxx_rhs(i,j,k) - gxx_rhs_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(gxy_rhs(i,j,k) - gxy_rhs_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(gxz_rhs(i,j,k) - gxz_rhs_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(gyy_rhs(i,j,k) - gyy_rhs_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(gyz_rhs(i,j,k) - gyz_rhs_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(gzz_rhs(i,j,k) - gzz_rhs_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Axx_rhs(i,j,k) - Axx_rhs_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Axy_rhs(i,j,k) - Axy_rhs_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Axz_rhs(i,j,k) - Axz_rhs_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Ayy_rhs(i,j,k) - Ayy_rhs_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Ayz_rhs(i,j,k) - Ayz_rhs_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Azz_rhs(i,j,k) - Azz_rhs_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Gamx_rhs(i,j,k) - Gamx_rhs_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Gamy_rhs(i,j,k) - Gamy_rhs_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Gamz_rhs(i,j,k) - Gamz_rhs_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Lap_rhs(i,j,k) - Lap_rhs_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(betax_rhs(i,j,k) - betax_rhs_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(betay_rhs(i,j,k) - betay_rhs_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(betaz_rhs(i,j,k) - betaz_rhs_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(dtSfx_rhs(i,j,k) - dtSfx_rhs_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(dtSfy_rhs(i,j,k) - dtSfy_rhs_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(dtSfz_rhs(i,j,k) - dtSfz_rhs_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         ! Check Christoffel symbols
!         curr_diff = abs(Gamxxx(i,j,k) - Gamxxx_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Gamxxy(i,j,k) - Gamxxy_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Gamxxz(i,j,k) - Gamxxz_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Gamxyy(i,j,k) - Gamxyy_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Gamxyz(i,j,k) - Gamxyz_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Gamxzz(i,j,k) - Gamxzz_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Gamyxx(i,j,k) - Gamyxx_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Gamyxy(i,j,k) - Gamyxy_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Gamyxz(i,j,k) - Gamyxz_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Gamyyy(i,j,k) - Gamyyy_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Gamyyz(i,j,k) - Gamyyz_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Gamyzz(i,j,k) - Gamyzz_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Gamzxx(i,j,k) - Gamzxx_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Gamzxy(i,j,k) - Gamzxy_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Gamzxz(i,j,k) - Gamzxz_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Gamzyy(i,j,k) - Gamzyy_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Gamzyz(i,j,k) - Gamzyz_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Gamzzz(i,j,k) - Gamzzz_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         ! Check Ricci tensor
!         curr_diff = abs(Rxx(i,j,k) - Rxx_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Rxy(i,j,k) - Rxy_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Rxz(i,j,k) - Rxz_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Ryy(i,j,k) - Ryy_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Ryz(i,j,k) - Ryz_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Rzz(i,j,k) - Rzz_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         ! Check constraints
!         curr_diff = abs(ham_Res(i,j,k) - ham_Res_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(movx_Res(i,j,k) - movx_Res_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(movy_Res(i,j,k) - movy_Res_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(movz_Res(i,j,k) - movz_Res_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Gmx_Res(i,j,k) - Gmx_Res_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Gmy_Res(i,j,k) - Gmy_Res_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
        
!         curr_diff = abs(Gmz_Res(i,j,k) - Gmz_Res_cor(i,j,k))
!         if (curr_diff > tol) diff_count = diff_count + 1
!         max_diff = max(max_diff, curr_diff)
!       enddo
!     enddo
!   enddo

! !~~~~~~> Print comparison results:
!   if (diff_count > 0) then
!     write(*,'(A)') '=========================================='
!     write(*,'(A)') 'WARNING: Differences detected between opt and cor versions!'
!     write(*,'(A,I10)') 'Number of differences above tolerance: ', diff_count
!     write(*,'(A,E15.7)') 'Maximum difference: ', max_diff
!     write(*,'(A,E15.7)') 'Tolerance: ', tol
!     write(*,'(A)') '=========================================='
!   else
!     write(*,'(A)') 'compute_rhs_bssn: opt and cor versions match within tolerance.'
!   endif

! !~~~~~~> Return opt version results (already in output arrays)
!   gont = gont_opt

!   return

contains

#if defined(Cell)
  function compute_rhs_bssn_opt(ex, T,X, Y, Z,                                     &
               chi    ,   trK    ,                                             &
               dxx    ,   gxy    ,   gxz    ,   dyy    ,   gyz    ,   dzz,     &
               Axx    ,   Axy    ,   Axz    ,   Ayy    ,   Ayz    ,   Azz,     &
               Gamx   ,  Gamy    ,  Gamz    ,                                  &
               Lap    ,  betax   ,  betay   ,  betaz   ,                       &
               dtSfx  ,  dtSfy   ,  dtSfz   ,                                  &
               chi_rhs,   trK_rhs,                                             &
               gxx_rhs,   gxy_rhs,   gxz_rhs,   gyy_rhs,   gyz_rhs,   gzz_rhs, &
               Axx_rhs,   Axy_rhs,   Axz_rhs,   Ayy_rhs,   Ayz_rhs,   Azz_rhs, &
               Gamx_rhs,  Gamy_rhs,  Gamz_rhs,                                 &
               Lap_rhs,  betax_rhs,  betay_rhs,  betaz_rhs,                    &
               dtSfx_rhs,  dtSfy_rhs,  dtSfz_rhs,                              &
               rho,Sx,Sy,Sz,Sxx,Sxy,Sxz,Syy,Syz,Szz,                           &
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,                      &
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,                      &
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,                      &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,                                        &
               ham_Res, movx_Res, movy_Res, movz_Res,                          &
                        Gmx_Res, Gmy_Res, Gmz_Res,                             &
               Symmetry,Lev,eps,co)  result(gont)
! calculate constraint violation when co=0               
  use omp_lib
  use sft_trace_mod
  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3), Symmetry,Lev,co
  real*8, intent(in ):: T
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: chi,dxx,dyy,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: trK
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: gxy,gxz,gyz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamx,Gamy,Gamz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Lap, betax, betay, betaz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: dtSfx,  dtSfy,  dtSfz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: chi_rhs,trK_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gxx_rhs,gxy_rhs,gxz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gyy_rhs,gyz_rhs,gzz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Axx_rhs,Axy_rhs,Axz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Ayy_rhs,Ayz_rhs,Azz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamx_rhs,Gamy_rhs,Gamz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Lap_rhs, betax_rhs, betay_rhs, betaz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: dtSfx_rhs,dtSfy_rhs,dtSfz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: rho,Sx,Sy,Sz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Sxx,Sxy,Sxz,Syy,Syz,Szz
! when out, physical second kind of connection  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamzyy, Gamzyz, Gamzzz
! when out, physical Ricci tensor  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Rxx,Rxy,Rxz,Ryy,Ryz,Rzz
  real*8,intent(in) :: eps
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: ham_Res, movx_Res, movy_Res, movz_Res
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Gmx_Res, Gmy_Res, Gmz_Res
!  gont = 0: success; gont = 1: something wrong
  integer::gont

!~~~~~~> Other variables:

  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gyy,gzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: chix,chiy,chiz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxx,gxyx,gxzx,gyyx,gyzx,gzzx
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxy,gxyy,gxzy,gyyy,gyzy,gzzy
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxz,gxyz,gxzz,gyyz,gyzz,gzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Lapx,Lapy,Lapz
  real*8, dimension(ex(1),ex(2),ex(3)) :: betaxx,betaxy,betaxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: betayx,betayy,betayz
  real*8, dimension(ex(1),ex(2),ex(3)) :: betazx,betazy,betazz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamxx,Gamxy,Gamxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamyx,Gamyy,Gamyz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamzx,Gamzy,Gamzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Kx,Ky,Kz,div_beta,S
  real*8, dimension(ex(1),ex(2),ex(3)) :: f,fxx,fxy,fxz,fyy,fyz,fzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamxa,Gamya,Gamza,alpn1,chin1
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupyy,gupyz,gupzz
  real*8,dimension(OPT_HALO_LO:ex(1),OPT_HALO_LO:ex(2),OPT_HALO_LO:ex(3)) :: fh_d, fh_d2, fh_d3, fh_d4, fh_d5, fh_d6
  real*8, dimension(ex(1),ex(2),ex(3)) :: fxx2,fxy2,fxz2,fyy2,fyz2,fzz2
  real*8,dimension(-2:ex(1),-2:ex(2),-2:ex(3)) :: fh_l1, fh_l2, fh_l3, fh_l4
  real*8,dimension(-2:ex(1),-2:ex(2),-2:ex(3)) :: fh_l5, fh_l6, fh_l7, fh_l8
  real*8,dimension(-2:ex(1),-2:ex(2),-2:ex(3)) :: fh_l9, fh_l10,fh_l11,fh_l12
  real*8,dimension(-2:ex(1),-2:ex(2),-2:ex(3)) :: fh_l13,fh_l14,fh_l15,fh_l16
  real*8,dimension(-2:ex(1),-2:ex(2),-2:ex(3)) :: fh_l17,fh_l18,fh_l19,fh_l20
  real*8,dimension(-2:ex(1),-2:ex(2),-2:ex(3)) :: fh_l21,fh_l22,fh_l23,fh_l24
  real*8,dimension(-2:ex(1),-2:ex(2),-2:ex(3)) :: fh_k1, fh_k2, fh_k3, fh_k4
  real*8,dimension(-2:ex(1),-2:ex(2),-2:ex(3)) :: fh_k5, fh_k6, fh_k7, fh_k8
  real*8,dimension(-2:ex(1),-2:ex(2),-2:ex(3)) :: fh_k9, fh_k10,fh_k11,fh_k12
  real*8,dimension(-2:ex(1),-2:ex(2),-2:ex(3)) :: fh_k13,fh_k14,fh_k15,fh_k16
  real*8,dimension(-2:ex(1),-2:ex(2),-2:ex(3)) :: fh_k17,fh_k18,fh_k19,fh_k20
  real*8,dimension(-2:ex(1),-2:ex(2),-2:ex(3)) :: fh_k21,fh_k22,fh_k23,fh_k24

  real*8,dimension(3) ::SSS,AAS,ASA,SAA,ASS,SAS,SSA
  real*8            :: dX, dY, dZ, PI
  real*8, parameter :: ZEO = 0.d0,ONE = 1.D0, TWO = 2.D0, FOUR = 4.D0
  real*8, parameter :: EIGHT = 8.D0, HALF = 0.5D0, THR = 3.d0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  double precision,parameter::FF = 0.75d0,eta=2.d0
  real*8, parameter :: F1o3 = 1.D0/3.D0, F2o3 = 2.D0/3.D0,F3o2=1.5d0, F1o6 = 1.D0/6.D0
  real*8, parameter :: F16=1.6d1,F8=8.d0

#if (GAUGE == 2 || GAUGE == 3 || GAUGE == 4 || GAUGE == 5)
  real*8, dimension(ex(1),ex(2),ex(3)) :: reta
#endif

#if (GAUGE == 6 || GAUGE == 7)
  integer :: BHN,i,j,k
  real*8, dimension(9) :: Porg
  real*8, dimension(3) :: Mass
  real*8 :: r1,r2,M,A,w1,w2,C1,C2
  real*8, dimension(ex(1),ex(2),ex(3)) :: reta

  call getpbh(BHN,Porg,Mass)
#endif

  integer :: i, j, k, tid, nthreads
  integer(8) :: ts0_thr  ! per-thread trace timestamp

!!! sanity check (disabled for production)
!  dX = sum(chi)+sum(trK)+sum(dxx)+sum(gxy)+sum(gxz)+sum(dyy)+sum(gyz)+sum(dzz) &
!      +sum(Axx)+sum(Axy)+sum(Axz)+sum(Ayy)+sum(Ayz)+sum(Azz)                   &
!      +sum(Gamx)+sum(Gamy)+sum(Gamz)                                           &
!      +sum(Lap)+sum(betax)+sum(betay)+sum(betaz)
!  if(dX.ne.dX) then
!     if(sum(chi).ne.sum(chi))write(*,*)"bssn.f90: find NaN in chi"
!     if(sum(trK).ne.sum(trK))write(*,*)"bssn.f90: find NaN in trk"
!     if(sum(dxx).ne.sum(dxx))write(*,*)"bssn.f90: find NaN in dxx"
!     if(sum(gxy).ne.sum(gxy))write(*,*)"bssn.f90: find NaN in gxy"
!     if(sum(gxz).ne.sum(gxz))write(*,*)"bssn.f90: find NaN in gxz"
!     if(sum(dyy).ne.sum(dyy))write(*,*)"bssn.f90: find NaN in dyy"
!     if(sum(gyz).ne.sum(gyz))write(*,*)"bssn.f90: find NaN in gyz"
!     if(sum(dzz).ne.sum(dzz))write(*,*)"bssn.f90: find NaN in dzz"
!     if(sum(Axx).ne.sum(Axx))write(*,*)"bssn.f90: find NaN in Axx"
!     if(sum(Axy).ne.sum(Axy))write(*,*)"bssn.f90: find NaN in Axy"
!     if(sum(Axz).ne.sum(Axz))write(*,*)"bssn.f90: find NaN in Axz"
!     if(sum(Ayy).ne.sum(Ayy))write(*,*)"bssn.f90: find NaN in Ayy"
!     if(sum(Ayz).ne.sum(Ayz))write(*,*)"bssn.f90: find NaN in Ayz"
!     if(sum(Azz).ne.sum(Azz))write(*,*)"bssn.f90: find NaN in Azz"
!     if(sum(Gamx).ne.sum(Gamx))write(*,*)"bssn.f90: find NaN in Gamx"
!     if(sum(Gamy).ne.sum(Gamy))write(*,*)"bssn.f90: find NaN in Gamy"
!     if(sum(Gamz).ne.sum(Gamz))write(*,*)"bssn.f90: find NaN in Gamz"
!     if(sum(Lap).ne.sum(Lap))write(*,*)"bssn.f90: find NaN in Lap"
!     if(sum(betax).ne.sum(betax))write(*,*)"bssn.f90: find NaN in betax"
!     if(sum(betay).ne.sum(betay))write(*,*)"bssn.f90: find NaN in betay"
!     if(sum(betaz).ne.sum(betaz))write(*,*)"bssn.f90: find NaN in betaz"
!     gont = 1
!     return
!  endif

  PI = dacos(-ONE)

  dX = X(2) - X(1)
  dY = Y(2) - Y(1)
  dZ = Z(2) - Z(1)

  SSS(1)=SYM;  SSS(2)=SYM;  SSS(3)=SYM
  AAS(1)=ANTI; AAS(2)=ANTI; AAS(3)=SYM
  ASA(1)=ANTI; ASA(2)=SYM;  ASA(3)=ANTI
  SAA(1)=SYM;  SAA(2)=ANTI; SAA(3)=ANTI
  ASS(1)=ANTI; ASS(2)=SYM;  ASS(3)=SYM
  SAS(1)=SYM;  SAS(2)=ANTI; SAS(3)=SYM
  SSA(1)=SYM;  SSA(2)=SYM;  SSA(3)=ANTI

!$OMP PARALLEL PRIVATE(i,j,k,ts0_thr,tid) DEFAULT(SHARED)
  tid = omp_get_thread_num()
  ts0_thr = sft_get_ts()
!$OMP DO COLLAPSE(1)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        alpn1(i,j,k) = Lap(i,j,k) + ONE
        chin1(i,j,k) = chi(i,j,k) + ONE
        gxx(i,j,k) = dxx(i,j,k) + ONE
        gyy(i,j,k) = dyy(i,j,k) + ONE
        gzz(i,j,k) = dzz(i,j,k) + ONE
      end do
    end do
  end do
!$OMP END DO NOWAIT
  call sft_trace("init_alphachi", ts0_thr, sft_get_ts(), 0, tid)
  ts0_thr = sft_get_ts()

!$OMP BARRIER
    call symmetry_bd_inner(OPT_SYMORD,ex,betax,fh_d,ASS)
    call symmetry_bd_inner(OPT_SYMORD,ex,betay,fh_d2,SAS)
    call symmetry_bd_inner(OPT_SYMORD,ex,betaz,fh_d3,SSA)
    call symmetry_bd_inner(OPT_SYMORD,ex,chi,fh_d4,SSS)
!$OMP BARRIER
    call fderivs_inner(ex,fh_d,betaxx,betaxy,betaxz,X,Y,Z,Symmetry)
    call fderivs_inner(ex,fh_d2,betayx,betayy,betayz,X,Y,Z,Symmetry)
    call fderivs_inner(ex,fh_d3,betazx,betazy,betazz,X,Y,Z,Symmetry)
    call fderivs_inner(ex,fh_d4,chix,chiy,chiz,X,Y,Z,symmetry)

!$OMP BARRIER
!$OMP DO COLLAPSE(1)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        div_beta(i,j,k) = betaxx(i,j,k) + betayy(i,j,k) + betazz(i,j,k)
        chi_rhs(i,j,k) = F2o3 *chin1(i,j,k)*( alpn1(i,j,k) * trK(i,j,k) - div_beta(i,j,k) ) !rhs for chi
      end do
    end do
  end do
!$OMP END DO NOWAIT
  call sft_trace("chi_rhs", ts0_thr, sft_get_ts(), 0, tid)
  ts0_thr = sft_get_ts()

    call symmetry_bd_inner(OPT_SYMORD,ex,dxx,fh_d,SSS)
    call symmetry_bd_inner(OPT_SYMORD,ex,gxy,fh_d2,AAS)
    call symmetry_bd_inner(OPT_SYMORD,ex,gxz,fh_d3,ASA)
    call symmetry_bd_inner(OPT_SYMORD,ex,dyy,fh_d4,SSS)
    call symmetry_bd_inner(OPT_SYMORD,ex,gyz,fh_d5,SAA)
    call symmetry_bd_inner(OPT_SYMORD,ex,dzz,fh_d6,SSS)
!$OMP BARRIER
    call fderivs_inner(ex,fh_d,gxxx,gxxy,gxxz,X,Y,Z,Symmetry)
    call fderivs_inner(ex,fh_d2,gxyx,gxyy,gxyz,X,Y,Z,Symmetry)
    call fderivs_inner(ex,fh_d3,gxzx,gxzy,gxzz,X,Y,Z,Symmetry)
    call fderivs_inner(ex,fh_d4,gyyx,gyyy,gyyz,X,Y,Z,Symmetry)
    call fderivs_inner(ex,fh_d5,gyzx,gyzy,gyzz,X,Y,Z,Symmetry)
    call fderivs_inner(ex,fh_d6,gzzx,gzzy,gzzz,X,Y,Z,Symmetry)

!$OMP BARRIER
!$OMP DO COLLAPSE(1)
  do k = 1, ex(3)
    do j = 1, ex(2)

      do i = 1, ex(1)
        gxx_rhs(i,j,k) = - TWO * alpn1(i,j,k) * Axx(i,j,k)    -  F2o3 * gxx(i,j,k) * div_beta(i,j,k)          + &
                TWO *(  gxx(i,j,k) * betaxx(i,j,k) +   gxy(i,j,k) * betayx(i,j,k) +   gxz(i,j,k) * betazx(i,j,k))

        gyy_rhs(i,j,k) = - TWO * alpn1(i,j,k) * Ayy(i,j,k)    -  F2o3 * gyy(i,j,k) * div_beta(i,j,k)          + &
                TWO *(  gxy(i,j,k) * betaxy(i,j,k) +   gyy(i,j,k) * betayy(i,j,k) +   gyz(i,j,k) * betazy(i,j,k))

        gzz_rhs(i,j,k) = - TWO * alpn1(i,j,k) * Azz(i,j,k)    -  F2o3 * gzz(i,j,k) * div_beta(i,j,k)          + &
                TWO *(  gxz(i,j,k) * betaxz(i,j,k) +   gyz(i,j,k) * betayz(i,j,k) +   gzz(i,j,k) * betazz(i,j,k))

        gxy_rhs(i,j,k) = - TWO * alpn1(i,j,k) * Axy(i,j,k)    +  F1o3 * gxy(i,j,k)    * div_beta(i,j,k)       + &
                        gxx(i,j,k) * betaxy(i,j,k)                  +   gxz(i,j,k) * betazy(i,j,k) + &
                                                gyy(i,j,k) * betayx(i,j,k) +   gyz(i,j,k) * betazx(i,j,k)   &
                                                        -   gxy(i,j,k) * betazz(i,j,k)

        gyz_rhs(i,j,k) = - TWO * alpn1(i,j,k) * Ayz(i,j,k)    +  F1o3 * gyz(i,j,k)    * div_beta(i,j,k)       + &
                        gxy(i,j,k) * betaxz(i,j,k) +   gyy(i,j,k) * betayz(i,j,k)                  + &
                        gxz(i,j,k) * betaxy(i,j,k)                  +   gzz(i,j,k) * betazy(i,j,k)   &
                                                        -   gyz(i,j,k) * betaxx(i,j,k)

        gxz_rhs(i,j,k) = - TWO * alpn1(i,j,k) * Axz(i,j,k)    +  F1o3 * gxz(i,j,k)    * div_beta(i,j,k)       + &
                        gxx(i,j,k) * betaxz(i,j,k) +   gxy(i,j,k) * betayz(i,j,k)                  + &
                                                gyz(i,j,k) * betayx(i,j,k) +   gzz(i,j,k) * betazx(i,j,k)   &
                                                        -   gxz(i,j,k) * betayy(i,j,k)     !rhs for gij
! invert tilted metric
        gupzz(i,j,k) =  gxx(i,j,k) * gyy(i,j,k) * gzz(i,j,k) + gxy(i,j,k) * gyz(i,j,k) * gxz(i,j,k) + gxz(i,j,k) * gxy(i,j,k) * gyz(i,j,k) - &
             gxz(i,j,k) * gyy(i,j,k) * gxz(i,j,k) - gxy(i,j,k) * gxy(i,j,k) * gzz(i,j,k) - gxx(i,j,k) * gyz(i,j,k) * gyz(i,j,k)
        gupxx(i,j,k) =   ( gyy(i,j,k) * gzz(i,j,k) - gyz(i,j,k) * gyz(i,j,k) ) / gupzz(i,j,k)
        gupxy(i,j,k) = - ( gxy(i,j,k) * gzz(i,j,k) - gyz(i,j,k) * gxz(i,j,k) ) / gupzz(i,j,k)
        gupxz(i,j,k) =   ( gxy(i,j,k) * gyz(i,j,k) - gyy(i,j,k) * gxz(i,j,k) ) / gupzz(i,j,k)
        gupyy(i,j,k) =   ( gxx(i,j,k) * gzz(i,j,k) - gxz(i,j,k) * gxz(i,j,k) ) / gupzz(i,j,k)
        gupyz(i,j,k) = - ( gxx(i,j,k) * gyz(i,j,k) - gxy(i,j,k) * gxz(i,j,k) ) / gupzz(i,j,k)
        gupzz(i,j,k) =   ( gxx(i,j,k) * gyy(i,j,k) - gxy(i,j,k) * gxy(i,j,k) ) / gupzz(i,j,k)
      end do

    end do
  end do
!$OMP END DO NOWAIT
  call sft_trace("gxx_rhs+gup", ts0_thr, sft_get_ts(), 0, tid)
  ts0_thr = sft_get_ts()

! second kind of connection + raise A^ij + [co==0] Gmx_Res — merged
!$OMP BARRIER
!$OMP DO COLLAPSE(1)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        Gamxxx(i,j,k) =HALF*( gupxx(i,j,k)*gxxx(i,j,k) + gupxy(i,j,k)*(TWO*gxyx(i,j,k) - gxxy(i,j,k) ) + gupxz(i,j,k)*(TWO*gxzx(i,j,k) - gxxz(i,j,k) ))
        Gamyxx(i,j,k) =HALF*( gupxy(i,j,k)*gxxx(i,j,k) + gupyy(i,j,k)*(TWO*gxyx(i,j,k) - gxxy(i,j,k) ) + gupyz(i,j,k)*(TWO*gxzx(i,j,k) - gxxz(i,j,k) ))
        Gamzxx(i,j,k) =HALF*( gupxz(i,j,k)*gxxx(i,j,k) + gupyz(i,j,k)*(TWO*gxyx(i,j,k) - gxxy(i,j,k) ) + gupzz(i,j,k)*(TWO*gxzx(i,j,k) - gxxz(i,j,k) ))

        Gamxyy(i,j,k) =HALF*( gupxx(i,j,k)*(TWO*gxyy(i,j,k) - gyyx(i,j,k) ) + gupxy(i,j,k)*gyyy(i,j,k) + gupxz(i,j,k)*(TWO*gyzy(i,j,k) - gyyz(i,j,k) ))
        Gamyyy(i,j,k) =HALF*( gupxy(i,j,k)*(TWO*gxyy(i,j,k) - gyyx(i,j,k) ) + gupyy(i,j,k)*gyyy(i,j,k) + gupyz(i,j,k)*(TWO*gyzy(i,j,k) - gyyz(i,j,k) ))
        Gamzyy(i,j,k) =HALF*( gupxz(i,j,k)*(TWO*gxyy(i,j,k) - gyyx(i,j,k) ) + gupyz(i,j,k)*gyyy(i,j,k) + gupzz(i,j,k)*(TWO*gyzy(i,j,k) - gyyz(i,j,k) ))

        Gamxzz(i,j,k) =HALF*( gupxx(i,j,k)*(TWO*gxzz(i,j,k) - gzzx(i,j,k) ) + gupxy(i,j,k)*(TWO*gyzz(i,j,k) - gzzy(i,j,k) ) + gupxz(i,j,k)*gzzz(i,j,k))
        Gamyzz(i,j,k) =HALF*( gupxy(i,j,k)*(TWO*gxzz(i,j,k) - gzzx(i,j,k) ) + gupyy(i,j,k)*(TWO*gyzz(i,j,k) - gzzy(i,j,k) ) + gupyz(i,j,k)*gzzz(i,j,k))
        Gamzzz(i,j,k) =HALF*( gupxz(i,j,k)*(TWO*gxzz(i,j,k) - gzzx(i,j,k) ) + gupyz(i,j,k)*(TWO*gyzz(i,j,k) - gzzy(i,j,k) ) + gupzz(i,j,k)*gzzz(i,j,k))

        Gamxxy(i,j,k) =HALF*( gupxx(i,j,k)*gxxy(i,j,k) + gupxy(i,j,k)*gyyx(i,j,k) + gupxz(i,j,k)*( gxzy(i,j,k) + gyzx(i,j,k) - gxyz(i,j,k) ) )
        Gamyxy(i,j,k) =HALF*( gupxy(i,j,k)*gxxy(i,j,k) + gupyy(i,j,k)*gyyx(i,j,k) + gupyz(i,j,k)*( gxzy(i,j,k) + gyzx(i,j,k) - gxyz(i,j,k) ) )
        Gamzxy(i,j,k) =HALF*( gupxz(i,j,k)*gxxy(i,j,k) + gupyz(i,j,k)*gyyx(i,j,k) + gupzz(i,j,k)*( gxzy(i,j,k) + gyzx(i,j,k) - gxyz(i,j,k) ) )

        Gamxxz(i,j,k) =HALF*( gupxx(i,j,k)*gxxz(i,j,k) + gupxy(i,j,k)*( gxyz(i,j,k) + gyzx(i,j,k) - gxzy(i,j,k) ) + gupxz(i,j,k)*gzzx(i,j,k) )
        Gamyxz(i,j,k) =HALF*( gupxy(i,j,k)*gxxz(i,j,k) + gupyy(i,j,k)*( gxyz(i,j,k) + gyzx(i,j,k) - gxzy(i,j,k) ) + gupyz(i,j,k)*gzzx(i,j,k) )
        Gamzxz(i,j,k) =HALF*( gupxz(i,j,k)*gxxz(i,j,k) + gupyz(i,j,k)*( gxyz(i,j,k) + gyzx(i,j,k) - gxzy(i,j,k) ) + gupzz(i,j,k)*gzzx(i,j,k) )

        Gamxyz(i,j,k) =HALF*( gupxx(i,j,k)*( gxyz(i,j,k) + gxzy(i,j,k) - gyzx(i,j,k) ) + gupxy(i,j,k)*gyyz(i,j,k) + gupxz(i,j,k)*gzzy(i,j,k) )
        Gamyyz(i,j,k) =HALF*( gupxy(i,j,k)*( gxyz(i,j,k) + gxzy(i,j,k) - gyzx(i,j,k) ) + gupyy(i,j,k)*gyyz(i,j,k) + gupyz(i,j,k)*gzzy(i,j,k) )
        Gamzyz(i,j,k) =HALF*( gupxz(i,j,k)*( gxyz(i,j,k) + gxzy(i,j,k) - gyzx(i,j,k) ) + gupyz(i,j,k)*gyyz(i,j,k) + gupzz(i,j,k)*gzzy(i,j,k) )
! Raise indices of \tilde A_{ij} and store in R_ij
        Rxx(i,j,k) =    gupxx(i,j,k) * gupxx(i,j,k) * Axx(i,j,k) + gupxy(i,j,k) * gupxy(i,j,k) * Ayy(i,j,k) + gupxz(i,j,k) * gupxz(i,j,k) * Azz(i,j,k) + &
            TWO*(gupxx(i,j,k) * gupxy(i,j,k) * Axy(i,j,k) + gupxx(i,j,k) * gupxz(i,j,k) * Axz(i,j,k) + gupxy(i,j,k) * gupxz(i,j,k) * Ayz(i,j,k))

        Ryy(i,j,k) =    gupxy(i,j,k) * gupxy(i,j,k) * Axx(i,j,k) + gupyy(i,j,k) * gupyy(i,j,k) * Ayy(i,j,k) + gupyz(i,j,k) * gupyz(i,j,k) * Azz(i,j,k) + &
            TWO*(gupxy(i,j,k) * gupyy(i,j,k) * Axy(i,j,k) + gupxy(i,j,k) * gupyz(i,j,k) * Axz(i,j,k) + gupyy(i,j,k) * gupyz(i,j,k) * Ayz(i,j,k))

        Rzz(i,j,k) =    gupxz(i,j,k) * gupxz(i,j,k) * Axx(i,j,k) + gupyz(i,j,k) * gupyz(i,j,k) * Ayy(i,j,k) + gupzz(i,j,k) * gupzz(i,j,k) * Azz(i,j,k) + &
            TWO*(gupxz(i,j,k) * gupyz(i,j,k) * Axy(i,j,k) + gupxz(i,j,k) * gupzz(i,j,k) * Axz(i,j,k) + gupyz(i,j,k) * gupzz(i,j,k) * Ayz(i,j,k))

        Rxy(i,j,k) =    gupxx(i,j,k) * gupxy(i,j,k) * Axx(i,j,k) + gupxy(i,j,k) * gupyy(i,j,k) * Ayy(i,j,k) + gupxz(i,j,k) * gupyz(i,j,k) * Azz(i,j,k) + &
                (gupxx(i,j,k) * gupyy(i,j,k)       + gupxy(i,j,k) * gupxy(i,j,k))* Axy(i,j,k)                       + &
                (gupxx(i,j,k) * gupyz(i,j,k)       + gupxz(i,j,k) * gupxy(i,j,k))* Axz(i,j,k)                       + &
                (gupxy(i,j,k) * gupyz(i,j,k)       + gupxz(i,j,k) * gupyy(i,j,k))* Ayz(i,j,k)

        Rxz(i,j,k) =    gupxx(i,j,k) * gupxz(i,j,k) * Axx(i,j,k) + gupxy(i,j,k) * gupyz(i,j,k) * Ayy(i,j,k) + gupxz(i,j,k) * gupzz(i,j,k) * Azz(i,j,k) + &
                (gupxx(i,j,k) * gupyz(i,j,k)       + gupxy(i,j,k) * gupxz(i,j,k))* Axy(i,j,k)                       + &
                (gupxx(i,j,k) * gupzz(i,j,k)       + gupxz(i,j,k) * gupxz(i,j,k))* Axz(i,j,k)                       + &
                (gupxy(i,j,k) * gupzz(i,j,k)       + gupxz(i,j,k) * gupyz(i,j,k))* Ayz(i,j,k)

        Ryz(i,j,k) =    gupxy(i,j,k) * gupxz(i,j,k) * Axx(i,j,k) + gupyy(i,j,k) * gupyz(i,j,k) * Ayy(i,j,k) + gupyz(i,j,k) * gupzz(i,j,k) * Azz(i,j,k) + &
                (gupxy(i,j,k) * gupyz(i,j,k)       + gupyy(i,j,k) * gupxz(i,j,k))* Axy(i,j,k)                       + &
                (gupxy(i,j,k) * gupzz(i,j,k)       + gupyz(i,j,k) * gupxz(i,j,k))* Axz(i,j,k)                       + &
                (gupyy(i,j,k) * gupzz(i,j,k)       + gupyz(i,j,k) * gupyz(i,j,k))* Ayz(i,j,k)
! Gam^i_Res (only when co==0)
        if(co == 0)then
        Gmx_Res(i,j,k) = Gamx(i,j,k) - (gupxx(i,j,k)*(gupxx(i,j,k)*gxxx(i,j,k)+gupxy(i,j,k)*gxyx(i,j,k)+gupxz(i,j,k)*gxzx(i,j,k))&
                     +gupxy(i,j,k)*(gupxx(i,j,k)*gxyx(i,j,k)+gupxy(i,j,k)*gyyx(i,j,k)+gupxz(i,j,k)*gyzx(i,j,k))&
                     +gupxz(i,j,k)*(gupxx(i,j,k)*gxzx(i,j,k)+gupxy(i,j,k)*gyzx(i,j,k)+gupxz(i,j,k)*gzzx(i,j,k))&
                     +gupxx(i,j,k)*(gupxy(i,j,k)*gxxy(i,j,k)+gupyy(i,j,k)*gxyy(i,j,k)+gupyz(i,j,k)*gxzy(i,j,k))&
                     +gupxy(i,j,k)*(gupxy(i,j,k)*gxyy(i,j,k)+gupyy(i,j,k)*gyyy(i,j,k)+gupyz(i,j,k)*gyzy(i,j,k))&
                     +gupxz(i,j,k)*(gupxy(i,j,k)*gxzy(i,j,k)+gupyy(i,j,k)*gyzy(i,j,k)+gupyz(i,j,k)*gzzy(i,j,k))&
                     +gupxx(i,j,k)*(gupxz(i,j,k)*gxxz(i,j,k)+gupyz(i,j,k)*gxyz(i,j,k)+gupzz(i,j,k)*gxzz(i,j,k))&
                     +gupxy(i,j,k)*(gupxz(i,j,k)*gxyz(i,j,k)+gupyz(i,j,k)*gyyz(i,j,k)+gupzz(i,j,k)*gyzz(i,j,k))&
                     +gupxz(i,j,k)*(gupxz(i,j,k)*gxzz(i,j,k)+gupyz(i,j,k)*gyzz(i,j,k)+gupzz(i,j,k)*gzzz(i,j,k)))
        Gmy_Res(i,j,k) = Gamy(i,j,k) - (gupxx(i,j,k)*(gupxy(i,j,k)*gxxx(i,j,k)+gupyy(i,j,k)*gxyx(i,j,k)+gupyz(i,j,k)*gxzx(i,j,k))&
                     +gupxy(i,j,k)*(gupxy(i,j,k)*gxyx(i,j,k)+gupyy(i,j,k)*gyyx(i,j,k)+gupyz(i,j,k)*gyzx(i,j,k))&
                     +gupxz(i,j,k)*(gupxy(i,j,k)*gxzx(i,j,k)+gupyy(i,j,k)*gyzx(i,j,k)+gupyz(i,j,k)*gzzx(i,j,k))&
                     +gupxy(i,j,k)*(gupxy(i,j,k)*gxxy(i,j,k)+gupyy(i,j,k)*gxyy(i,j,k)+gupyz(i,j,k)*gxzy(i,j,k))&
                     +gupyy(i,j,k)*(gupxy(i,j,k)*gxyy(i,j,k)+gupyy(i,j,k)*gyyy(i,j,k)+gupyz(i,j,k)*gyzy(i,j,k))&
                     +gupyz(i,j,k)*(gupxy(i,j,k)*gxzy(i,j,k)+gupyy(i,j,k)*gyzy(i,j,k)+gupyz(i,j,k)*gzzy(i,j,k))&
                     +gupxy(i,j,k)*(gupxz(i,j,k)*gxxz(i,j,k)+gupyz(i,j,k)*gxyz(i,j,k)+gupzz(i,j,k)*gxzz(i,j,k))&
                     +gupyy(i,j,k)*(gupxz(i,j,k)*gxyz(i,j,k)+gupyz(i,j,k)*gyyz(i,j,k)+gupzz(i,j,k)*gyzz(i,j,k))&
                     +gupyz(i,j,k)*(gupxz(i,j,k)*gxzz(i,j,k)+gupyz(i,j,k)*gyzz(i,j,k)+gupzz(i,j,k)*gzzz(i,j,k)))
        Gmz_Res(i,j,k) = Gamz(i,j,k) - (gupxx(i,j,k)*(gupxz(i,j,k)*gxxx(i,j,k)+gupyz(i,j,k)*gxyx(i,j,k)+gupzz(i,j,k)*gxzx(i,j,k))&
                     +gupxy(i,j,k)*(gupxz(i,j,k)*gxyx(i,j,k)+gupyz(i,j,k)*gyyx(i,j,k)+gupzz(i,j,k)*gyzx(i,j,k))&
                     +gupxz(i,j,k)*(gupxz(i,j,k)*gxzx(i,j,k)+gupyz(i,j,k)*gyzx(i,j,k)+gupzz(i,j,k)*gzzx(i,j,k))&
                     +gupxy(i,j,k)*(gupxz(i,j,k)*gxxy(i,j,k)+gupyz(i,j,k)*gxyy(i,j,k)+gupzz(i,j,k)*gxzy(i,j,k))&
                     +gupyy(i,j,k)*(gupxz(i,j,k)*gxyy(i,j,k)+gupyz(i,j,k)*gyyy(i,j,k)+gupzz(i,j,k)*gyzy(i,j,k))&
                     +gupyz(i,j,k)*(gupxz(i,j,k)*gxzy(i,j,k)+gupyz(i,j,k)*gyzy(i,j,k)+gupzz(i,j,k)*gzzy(i,j,k))&
                     +gupxz(i,j,k)*(gupxz(i,j,k)*gxxz(i,j,k)+gupyz(i,j,k)*gxyz(i,j,k)+gupzz(i,j,k)*gxzz(i,j,k))&
                     +gupyz(i,j,k)*(gupxz(i,j,k)*gxyz(i,j,k)+gupyz(i,j,k)*gyyz(i,j,k)+gupzz(i,j,k)*gyzz(i,j,k))&
                     +gupzz(i,j,k)*(gupxz(i,j,k)*gxzz(i,j,k)+gupyz(i,j,k)*gyzz(i,j,k)+gupzz(i,j,k)*gzzz(i,j,k)))
        endif
      end do
    end do
  end do
!$OMP END DO NOWAIT
  call sft_trace("Gamxxx+Rxx+Gmx_Res", ts0_thr, sft_get_ts(), 0, tid)
  ts0_thr = sft_get_ts()

! Batch Lap + trK + Gam^i first derivatives
!$OMP BARRIER
    call symmetry_bd_inner(OPT_SYMORD,ex,Lap,fh_d,SSS)
    call symmetry_bd_inner(OPT_SYMORD,ex,trK,fh_d2,SSS)
    call symmetry_bd_inner(OPT_SYMORD,ex,Gamx,fh_d3,ASS)
    call symmetry_bd_inner(OPT_SYMORD,ex,Gamy,fh_d4,SAS)
    call symmetry_bd_inner(OPT_SYMORD,ex,Gamz,fh_d5,SSA)
!$OMP BARRIER
    call fderivs_inner(ex,fh_d,Lapx,Lapy,Lapz,X,Y,Z,Symmetry)
    call fderivs_inner(ex,fh_d2,Kx,Ky,Kz,X,Y,Z,symmetry)
    call fderivs_inner(ex,fh_d3,Gamxx,Gamxy,Gamxz,X,Y,Z,Symmetry)
    call fderivs_inner(ex,fh_d4,Gamyx,Gamyy,Gamyz,X,Y,Z,Symmetry)
    call fderivs_inner(ex,fh_d5,Gamzx,Gamzy,Gamzz,X,Y,Z,Symmetry)

! Gamx_rhs partial + Gamxa (merged) + beta 2nd deriv sym_bd overlap
!$OMP BARRIER
!$OMP DO COLLAPSE(1)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        Gamx_rhs(i,j,k) = - TWO * (   Lapx(i,j,k) * Rxx(i,j,k) +   Lapy(i,j,k) * Rxy(i,j,k) +   Lapz(i,j,k) * Rxz(i,j,k) ) + &
             TWO * alpn1(i,j,k) * (                                                &
             -F3o2/chin1(i,j,k) * (   chix(i,j,k) * Rxx(i,j,k) +   chiy(i,j,k) * Rxy(i,j,k) +   chiz(i,j,k) * Rxz(i,j,k) ) - &
                   gupxx(i,j,k) * (   F2o3 * Kx(i,j,k)  +  EIGHT * PI * Sx(i,j,k)            ) - &
                   gupxy(i,j,k) * (   F2o3 * Ky(i,j,k)  +  EIGHT * PI * Sy(i,j,k)            ) - &
                   gupxz(i,j,k) * (   F2o3 * Kz(i,j,k)  +  EIGHT * PI * Sz(i,j,k)            ) + &
                             Gamxxx(i,j,k) * Rxx(i,j,k) + Gamxyy(i,j,k) * Ryy(i,j,k) + Gamxzz(i,j,k) * Rzz(i,j,k)   + &
                     TWO * ( Gamxxy(i,j,k) * Rxy(i,j,k) + Gamxxz(i,j,k) * Rxz(i,j,k) + Gamxyz(i,j,k) * Ryz(i,j,k) ) )

        Gamy_rhs(i,j,k) = - TWO * (   Lapx(i,j,k) * Rxy(i,j,k) +   Lapy(i,j,k) * Ryy(i,j,k) +   Lapz(i,j,k) * Ryz(i,j,k) ) + &
             TWO * alpn1(i,j,k) * (                                                &
             -F3o2/chin1(i,j,k) * (   chix(i,j,k) * Rxy(i,j,k) +  chiy(i,j,k) * Ryy(i,j,k) +    chiz(i,j,k) * Ryz(i,j,k) ) - &
                   gupxy(i,j,k) * (   F2o3 * Kx(i,j,k)  +  EIGHT * PI * Sx(i,j,k)            ) - &
                   gupyy(i,j,k) * (   F2o3 * Ky(i,j,k)  +  EIGHT * PI * Sy(i,j,k)            ) - &
                   gupyz(i,j,k) * (   F2o3 * Kz(i,j,k)  +  EIGHT * PI * Sz(i,j,k)            ) + &
                             Gamyxx(i,j,k) * Rxx(i,j,k) + Gamyyy(i,j,k) * Ryy(i,j,k) + Gamyzz(i,j,k) * Rzz(i,j,k)   + &
                     TWO * ( Gamyxy(i,j,k) * Rxy(i,j,k) + Gamyxz(i,j,k) * Rxz(i,j,k) + Gamyyz(i,j,k) * Ryz(i,j,k) ) )

        Gamz_rhs(i,j,k) = - TWO * (   Lapx(i,j,k) * Rxz(i,j,k) +   Lapy(i,j,k) * Ryz(i,j,k) +   Lapz(i,j,k) * Rzz(i,j,k) ) + &
             TWO * alpn1(i,j,k) * (                                                &
             -F3o2/chin1(i,j,k) * (   chix(i,j,k) * Rxz(i,j,k) +  chiy(i,j,k) * Ryz(i,j,k) +    chiz(i,j,k) * Rzz(i,j,k) ) - &
                   gupxz(i,j,k) * (   F2o3 * Kx(i,j,k) +  EIGHT * PI * Sx(i,j,k)            ) - &
                   gupyz(i,j,k) * (   F2o3 * Ky(i,j,k)  +  EIGHT * PI * Sy(i,j,k)            ) - &
                   gupzz(i,j,k) * (   F2o3 * Kz(i,j,k)  +  EIGHT * PI * Sz(i,j,k)            ) + &
                             Gamzxx(i,j,k) * Rxx(i,j,k) + Gamzyy(i,j,k) * Ryy(i,j,k) + Gamzzz(i,j,k) * Rzz(i,j,k)   + &
                     TWO * ( Gamzxy(i,j,k) * Rxy(i,j,k) + Gamzxz(i,j,k) * Rxz(i,j,k) + Gamzyz(i,j,k) * Ryz(i,j,k) ) )

        Gamxa(i,j,k) =       gupxx(i,j,k) * Gamxxx(i,j,k) + gupyy(i,j,k) * Gamxyy(i,j,k) + gupzz(i,j,k) * Gamxzz(i,j,k) + &
                TWO*( gupxy(i,j,k) * Gamxxy(i,j,k) + gupxz(i,j,k) * Gamxxz(i,j,k) + gupyz(i,j,k) * Gamxyz(i,j,k) )
        Gamya(i,j,k) =       gupxx(i,j,k) * Gamyxx(i,j,k) + gupyy(i,j,k) * Gamyyy(i,j,k) + gupzz(i,j,k) * Gamyzz(i,j,k) + &
                TWO*( gupxy(i,j,k) * Gamyxy(i,j,k) + gupxz(i,j,k) * Gamyxz(i,j,k) + gupyz(i,j,k) * Gamyyz(i,j,k) )
        Gamza(i,j,k) =       gupxx(i,j,k) * Gamzxx(i,j,k) + gupyy(i,j,k) * Gamzyy(i,j,k) + gupzz(i,j,k) * Gamzzz(i,j,k) + &
                TWO*( gupxy(i,j,k) * Gamzxy(i,j,k) + gupxz(i,j,k) * Gamzxz(i,j,k) + gupyz(i,j,k) * Gamzyz(i,j,k) )
      end do
    end do
  end do
!$OMP END DO NOWAIT
  call sft_trace("Gamx_rhs+Gamxa", ts0_thr, sft_get_ts(), 0, tid)
  ts0_thr = sft_get_ts()

! Batch beta 2nd derivatives (overlap with above NOWAIT DO)
    call symmetry_bd_inner(OPT_SYMORD,ex,betax,fh_d,ASS)
    call symmetry_bd_inner(OPT_SYMORD,ex,betay,fh_d2,SAS)
    call symmetry_bd_inner(OPT_SYMORD,ex,betaz,fh_d3,SSA)
!$OMP BARRIER
    call fdderivs_inner(ex,fh_d,gxxx,gxyx,gxzx,gyyx,gyzx,gzzx,X,Y,Z,Symmetry)
    call fdderivs_inner(ex,fh_d2,gxxy,gxyy,gxzy,gyyy,gyzy,gzzy,X,Y,Z,Symmetry)
    call fdderivs_inner(ex,fh_d3,gxxz,gxyz,gxzz,gyyz,gyzz,gzzz,X,Y,Z,Symmetry)

! fxx + Gamx_rhs_complete (merged)
!$OMP BARRIER
!$OMP DO COLLAPSE(1)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        fxx(i,j,k) = gxxx(i,j,k) + gxyy(i,j,k) + gxzz(i,j,k)
        fxy(i,j,k) = gxyx(i,j,k) + gyyy(i,j,k) + gyzz(i,j,k)
        fxz(i,j,k) = gxzx(i,j,k) + gyzy(i,j,k) + gzzz(i,j,k)

        Gamx_rhs(i,j,k) =               Gamx_rhs(i,j,k) +  F2o3 *  Gamxa(i,j,k) * div_beta(i,j,k)        - &
                       Gamxa(i,j,k) * betaxx(i,j,k) - Gamya(i,j,k) * betaxy(i,j,k) - Gamza(i,j,k) * betaxz(i,j,k)  + &
               F1o3 * (gupxx(i,j,k) * fxx(i,j,k)    + gupxy(i,j,k) * fxy(i,j,k)    + gupxz(i,j,k) * fxz(i,j,k)    ) + &
                       gupxx(i,j,k) * gxxx(i,j,k)   + gupyy(i,j,k) * gyyx(i,j,k)   + gupzz(i,j,k) * gzzx(i,j,k)    + &
                TWO * (gupxy(i,j,k) * gxyx(i,j,k)   + gupxz(i,j,k) * gxzx(i,j,k)   + gupyz(i,j,k) * gyzx(i,j,k)  )

        Gamy_rhs(i,j,k) =               Gamy_rhs(i,j,k) +  F2o3 *  Gamya(i,j,k) * div_beta(i,j,k)        - &
                       Gamxa(i,j,k) * betayx(i,j,k) - Gamya(i,j,k) * betayy(i,j,k) - Gamza(i,j,k) * betayz(i,j,k)  + &
               F1o3 * (gupxy(i,j,k) * fxx(i,j,k)    + gupyy(i,j,k) * fxy(i,j,k)    + gupyz(i,j,k) * fxz(i,j,k)    ) + &
                       gupxx(i,j,k) * gxxy(i,j,k)   + gupyy(i,j,k) * gyyy(i,j,k)   + gupzz(i,j,k) * gzzy(i,j,k)    + &
                TWO * (gupxy(i,j,k) * gxyy(i,j,k)   + gupxz(i,j,k) * gxzy(i,j,k)   + gupyz(i,j,k) * gyzy(i,j,k)  )

        Gamz_rhs(i,j,k) =               Gamz_rhs(i,j,k) +  F2o3 *  Gamza(i,j,k) * div_beta(i,j,k)        - &
                       Gamxa(i,j,k) * betazx(i,j,k) - Gamya(i,j,k) * betazy(i,j,k) - Gamza(i,j,k) * betazz(i,j,k)  + &
               F1o3 * (gupxz(i,j,k) * fxx(i,j,k)    + gupyz(i,j,k) * fxy(i,j,k)    + gupzz(i,j,k) * fxz(i,j,k)    ) + &
                       gupxx(i,j,k) * gxxz(i,j,k)   + gupyy(i,j,k) * gyyz(i,j,k)   + gupzz(i,j,k) * gzzz(i,j,k)    + &
                TWO * (gupxy(i,j,k) * gxyz(i,j,k)   + gupxz(i,j,k) * gxzz(i,j,k)   + gupyz(i,j,k) * gyzz(i,j,k)  )    !rhs for Gam^i
      end do
    end do
  end do
!$OMP END DO NOWAIT
  call sft_trace("Gamx_rhs_complete", ts0_thr, sft_get_ts(), 0, tid)
  ts0_thr = sft_get_ts()

!first kind of connection stored in gij,k
!! Complex Ricci tensor calculations - using explicit DO + PARALLEL DO
!$OMP BARRIER
!$OMP DO COLLAPSE(1)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        gxxx(i,j,k) = gxx(i,j,k) * Gamxxx(i,j,k) + gxy(i,j,k) * Gamyxx(i,j,k) + gxz(i,j,k) * Gamzxx(i,j,k)
        gxyx(i,j,k) = gxx(i,j,k) * Gamxxy(i,j,k) + gxy(i,j,k) * Gamyxy(i,j,k) + gxz(i,j,k) * Gamzxy(i,j,k)
        gxzx(i,j,k) = gxx(i,j,k) * Gamxxz(i,j,k) + gxy(i,j,k) * Gamyxz(i,j,k) + gxz(i,j,k) * Gamzxz(i,j,k)
        gyyx(i,j,k) = gxx(i,j,k) * Gamxyy(i,j,k) + gxy(i,j,k) * Gamyyy(i,j,k) + gxz(i,j,k) * Gamzyy(i,j,k)
        gyzx(i,j,k) = gxx(i,j,k) * Gamxyz(i,j,k) + gxy(i,j,k) * Gamyyz(i,j,k) + gxz(i,j,k) * Gamzyz(i,j,k)
        gzzx(i,j,k) = gxx(i,j,k) * Gamxzz(i,j,k) + gxy(i,j,k) * Gamyzz(i,j,k) + gxz(i,j,k) * Gamzzz(i,j,k)

        gxxy(i,j,k) = gxy(i,j,k) * Gamxxx(i,j,k) + gyy(i,j,k) * Gamyxx(i,j,k) + gyz(i,j,k) * Gamzxx(i,j,k)
        gxyy(i,j,k) = gxy(i,j,k) * Gamxxy(i,j,k) + gyy(i,j,k) * Gamyxy(i,j,k) + gyz(i,j,k) * Gamzxy(i,j,k)
        gxzy(i,j,k) = gxy(i,j,k) * Gamxxz(i,j,k) + gyy(i,j,k) * Gamyxz(i,j,k) + gyz(i,j,k) * Gamzxz(i,j,k)
        gyyy(i,j,k) = gxy(i,j,k) * Gamxyy(i,j,k) + gyy(i,j,k) * Gamyyy(i,j,k) + gyz(i,j,k) * Gamzyy(i,j,k)
        gyzy(i,j,k) = gxy(i,j,k) * Gamxyz(i,j,k) + gyy(i,j,k) * Gamyyz(i,j,k) + gyz(i,j,k) * Gamzyz(i,j,k)
        gzzy(i,j,k) = gxy(i,j,k) * Gamxzz(i,j,k) + gyy(i,j,k) * Gamyzz(i,j,k) + gyz(i,j,k) * Gamzzz(i,j,k)

        gxxz(i,j,k) = gxz(i,j,k) * Gamxxx(i,j,k) + gyz(i,j,k) * Gamyxx(i,j,k) + gzz(i,j,k) * Gamzxx(i,j,k)
        gxyz(i,j,k) = gxz(i,j,k) * Gamxxy(i,j,k) + gyz(i,j,k) * Gamyxy(i,j,k) + gzz(i,j,k) * Gamzxy(i,j,k)
        gxzz(i,j,k) = gxz(i,j,k) * Gamxxz(i,j,k) + gyz(i,j,k) * Gamyxz(i,j,k) + gzz(i,j,k) * Gamzxz(i,j,k)
        gyyz(i,j,k) = gxz(i,j,k) * Gamxyy(i,j,k) + gyz(i,j,k) * Gamyyy(i,j,k) + gzz(i,j,k) * Gamzyy(i,j,k)
        gyzz(i,j,k) = gxz(i,j,k) * Gamxyz(i,j,k) + gyz(i,j,k) * Gamyyz(i,j,k) + gzz(i,j,k) * Gamzyz(i,j,k)
        gzzz(i,j,k) = gxz(i,j,k) * Gamxzz(i,j,k) + gyz(i,j,k) * Gamyzz(i,j,k) + gzz(i,j,k) * Gamzzz(i,j,k)
      end do
    end do
  end do
!$OMP END DO NOWAIT
  call sft_trace("gxxx", ts0_thr, sft_get_ts(), 0, tid)
  ts0_thr = sft_get_ts()

!compute Ricci tensor for tilted metric — batched 2 at a time
! Batch 1: dxx + dyy → Rxx, Ryy
     call symmetry_bd_inner(OPT_SYMORD,ex,dxx,fh_d,SSS)
     call symmetry_bd_inner(OPT_SYMORD,ex,dyy,fh_d2,SSS)
!$OMP BARRIER
     call fdderivs_inner(ex,fh_d,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,symmetry)
     call fdderivs_inner(ex,fh_d2,fxx2,fxy2,fxz2,fyy2,fyz2,fzz2,X,Y,Z,symmetry)
!$OMP BARRIER
!$OMP DO COLLAPSE(1)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        Rxx(i,j,k) = gupxx(i,j,k)*fxx(i,j,k) + gupyy(i,j,k)*fyy(i,j,k) + gupzz(i,j,k)*fzz(i,j,k) + &
                   ( gupxy(i,j,k)*fxy(i,j,k)  + gupxz(i,j,k)*fxz(i,j,k) + gupyz(i,j,k)*fyz(i,j,k) ) * TWO
        Ryy(i,j,k) = gupxx(i,j,k)*fxx2(i,j,k) + gupyy(i,j,k)*fyy2(i,j,k) + gupzz(i,j,k)*fzz2(i,j,k) + &
                   ( gupxy(i,j,k)*fxy2(i,j,k)  + gupxz(i,j,k)*fxz2(i,j,k) + gupyz(i,j,k)*fyz2(i,j,k) ) * TWO
      end do
    end do
  end do
!$OMP END DO NOWAIT
  call sft_trace("Rxx_Ryy", ts0_thr, sft_get_ts(), 0, tid)
  ts0_thr = sft_get_ts()

! Batch 2: dzz + gxy → Rzz, Rxy
     call symmetry_bd_inner(OPT_SYMORD,ex,dzz,fh_d,SSS)
     call symmetry_bd_inner(OPT_SYMORD,ex,gxy,fh_d2,AAS)
!$OMP BARRIER
     call fdderivs_inner(ex,fh_d,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,symmetry)
     call fdderivs_inner(ex,fh_d2,fxx2,fxy2,fxz2,fyy2,fyz2,fzz2,X,Y,Z,symmetry)
!$OMP BARRIER
!$OMP DO COLLAPSE(1)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        Rzz(i,j,k) = gupxx(i,j,k)*fxx(i,j,k) + gupyy(i,j,k)*fyy(i,j,k) + gupzz(i,j,k)*fzz(i,j,k) + &
                   ( gupxy(i,j,k)*fxy(i,j,k)  + gupxz(i,j,k)*fxz(i,j,k) + gupyz(i,j,k)*fyz(i,j,k) ) * TWO
        Rxy(i,j,k) = gupxx(i,j,k)*fxx2(i,j,k) + gupyy(i,j,k)*fyy2(i,j,k) + gupzz(i,j,k)*fzz2(i,j,k) + &
                   ( gupxy(i,j,k)*fxy2(i,j,k)  + gupxz(i,j,k)*fxz2(i,j,k) + gupyz(i,j,k)*fyz2(i,j,k) ) * TWO
      end do
    end do
  end do
!$OMP END DO NOWAIT
  call sft_trace("Rzz_Rxy", ts0_thr, sft_get_ts(), 0, tid)
  ts0_thr = sft_get_ts()

! Batch 3: gxz + gyz → Rxz, Ryz
     call symmetry_bd_inner(OPT_SYMORD,ex,gxz,fh_d,ASA)
     call symmetry_bd_inner(OPT_SYMORD,ex,gyz,fh_d2,SAA)
!$OMP BARRIER
     call fdderivs_inner(ex,fh_d,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,symmetry)
     call fdderivs_inner(ex,fh_d2,fxx2,fxy2,fxz2,fyy2,fyz2,fzz2,X,Y,Z,symmetry)
!$OMP BARRIER
!$OMP DO COLLAPSE(1)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        Rxz(i,j,k) = gupxx(i,j,k)*fxx(i,j,k) + gupyy(i,j,k)*fyy(i,j,k) + gupzz(i,j,k)*fzz(i,j,k) + &
                   ( gupxy(i,j,k)*fxy(i,j,k)  + gupxz(i,j,k)*fxz(i,j,k) + gupyz(i,j,k)*fyz(i,j,k) ) * TWO
        Ryz(i,j,k) = gupxx(i,j,k)*fxx2(i,j,k) + gupyy(i,j,k)*fyy2(i,j,k) + gupzz(i,j,k)*fzz2(i,j,k) + &
                   ( gupxy(i,j,k)*fxy2(i,j,k)  + gupxz(i,j,k)*fxz2(i,j,k) + gupyz(i,j,k)*fyz2(i,j,k) ) * TWO
      end do
    end do
  end do
!$OMP END DO NOWAIT
  call sft_trace("Rxz_Ryz", ts0_thr, sft_get_ts(), 0, tid)
  ts0_thr = sft_get_ts()

!$OMP BARRIER
!$OMP DO COLLAPSE(1)
  do k = 1, ex(3)
    do j = 1, ex(2)

      do i = 1, ex(1)
        Rxx(i,j,k) =     - HALF * Rxx(i,j,k)                                   + &
                   gxx(i,j,k) * Gamxx(i,j,k)+ gxy(i,j,k) * Gamyx(i,j,k)   +    gxz(i,j,k) * Gamzx(i,j,k) + &
                 Gamxa(i,j,k) * gxxx(i,j,k) +  Gamya(i,j,k) * gxyx(i,j,k) +  Gamza(i,j,k) * gxzx(i,j,k)  + &
         gupxx(i,j,k) *(                                                  &
             TWO*(Gamxxx(i,j,k) * gxxx(i,j,k) + Gamyxx(i,j,k) * gxyx(i,j,k) + Gamzxx(i,j,k) * gxzx(i,j,k)) + &
                  Gamxxx(i,j,k) * gxxx(i,j,k) + Gamyxx(i,j,k) * gxxy(i,j,k) + Gamzxx(i,j,k) * gxxz(i,j,k) )+ &
         gupxy(i,j,k) *(                                                  &
             TWO*(Gamxxx(i,j,k) * gxyx(i,j,k) + Gamyxx(i,j,k) * gyyx(i,j,k) + Gamzxx(i,j,k) * gyzx(i,j,k)  + &
                  Gamxxy(i,j,k) * gxxx(i,j,k) + Gamyxy(i,j,k) * gxyx(i,j,k) + Gamzxy(i,j,k) * gxzx(i,j,k)) + &
                  Gamxxy(i,j,k) * gxxx(i,j,k) + Gamyxy(i,j,k) * gxxy(i,j,k) + Gamzxy(i,j,k) * gxxz(i,j,k)  + &
                  Gamxxx(i,j,k) * gxyx(i,j,k) + Gamyxx(i,j,k) * gxyy(i,j,k) + Gamzxx(i,j,k) * gxyz(i,j,k) )+ &
         gupxz(i,j,k) *(                                                  &
             TWO*(Gamxxx(i,j,k) * gxzx(i,j,k) + Gamyxx(i,j,k) * gyzx(i,j,k) + Gamzxx(i,j,k) * gzzx(i,j,k)  + &
                  Gamxxz(i,j,k) * gxxx(i,j,k) + Gamyxz(i,j,k) * gxyx(i,j,k) + Gamzxz(i,j,k) * gxzx(i,j,k)) + &
                  Gamxxz(i,j,k) * gxxx(i,j,k) + Gamyxz(i,j,k) * gxxy(i,j,k) + Gamzxz(i,j,k) * gxxz(i,j,k)  + &
                  Gamxxx(i,j,k) * gxzx(i,j,k) + Gamyxx(i,j,k) * gxzy(i,j,k) + Gamzxx(i,j,k) * gxzz(i,j,k) )+ &
         gupyy(i,j,k) *(                                                  &
             TWO*(Gamxxy(i,j,k) * gxyx(i,j,k) + Gamyxy(i,j,k) * gyyx(i,j,k) + Gamzxy(i,j,k) * gyzx(i,j,k)) + &
                  Gamxxy(i,j,k) * gxyx(i,j,k) + Gamyxy(i,j,k) * gxyy(i,j,k) + Gamzxy(i,j,k) * gxyz(i,j,k) )+ &
         gupyz(i,j,k) *(                                                  &
             TWO*(Gamxxy(i,j,k) * gxzx(i,j,k) + Gamyxy(i,j,k) * gyzx(i,j,k) + Gamzxy(i,j,k) * gzzx(i,j,k)  + &
                  Gamxxz(i,j,k) * gxyx(i,j,k) + Gamyxz(i,j,k) * gyyx(i,j,k) + Gamzxz(i,j,k) * gyzx(i,j,k)) + &
                  Gamxxz(i,j,k) * gxyx(i,j,k) + Gamyxz(i,j,k) * gxyy(i,j,k) + Gamzxz(i,j,k) * gxyz(i,j,k)  + &
                  Gamxxy(i,j,k) * gxzx(i,j,k) + Gamyxy(i,j,k) * gxzy(i,j,k) + Gamzxy(i,j,k) * gxzz(i,j,k) )+ &
         gupzz(i,j,k) *(                                                  &
             TWO*(Gamxxz(i,j,k) * gxzx(i,j,k) + Gamyxz(i,j,k) * gyzx(i,j,k) + Gamzxz(i,j,k) * gzzx(i,j,k)) + &
                  Gamxxz(i,j,k) * gxzx(i,j,k) + Gamyxz(i,j,k) * gxzy(i,j,k) + Gamzxz(i,j,k) * gxzz(i,j,k) )

        Ryy(i,j,k) =     - HALF * Ryy(i,j,k)                                   + &
                   gxy(i,j,k) * Gamxy(i,j,k)+  gyy(i,j,k) * Gamyy(i,j,k)  +  gyz(i,j,k) * Gamzy(i,j,k)   + &
                 Gamxa(i,j,k) * gxyy(i,j,k) +  Gamya(i,j,k) * gyyy(i,j,k) +  Gamza(i,j,k) * gyzy(i,j,k)  + &
         gupxx(i,j,k) *(                                                  &
             TWO*(Gamxxy(i,j,k) * gxxy(i,j,k) + Gamyxy(i,j,k) * gxyy(i,j,k) + Gamzxy(i,j,k) * gxzy(i,j,k)) + &
                  Gamxxy(i,j,k) * gxyx(i,j,k) + Gamyxy(i,j,k) * gxyy(i,j,k) + Gamzxy(i,j,k) * gxyz(i,j,k) )+ &
         gupxy(i,j,k) *(                                                  &
             TWO*(Gamxxy(i,j,k) * gxyy(i,j,k) + Gamyxy(i,j,k) * gyyy(i,j,k) + Gamzxy(i,j,k) * gyzy(i,j,k)  + &
                  Gamxyy(i,j,k) * gxxy(i,j,k) + Gamyyy(i,j,k) * gxyy(i,j,k) + Gamzyy(i,j,k) * gxzy(i,j,k)) + &
                  Gamxyy(i,j,k) * gxyx(i,j,k) + Gamyyy(i,j,k) * gxyy(i,j,k) + Gamzyy(i,j,k) * gxyz(i,j,k)  + &
                  Gamxxy(i,j,k) * gyyx(i,j,k) + Gamyxy(i,j,k) * gyyy(i,j,k) + Gamzxy(i,j,k) * gyyz(i,j,k) )+ &
         gupxz(i,j,k) *(                                                  &
             TWO*(Gamxxy(i,j,k) * gxzy(i,j,k) + Gamyxy(i,j,k) * gyzy(i,j,k) + Gamzxy(i,j,k) * gzzy(i,j,k)  + &
                  Gamxyz(i,j,k) * gxxy(i,j,k) + Gamyyz(i,j,k) * gxyy(i,j,k) + Gamzyz(i,j,k) * gxzy(i,j,k)) + &
                  Gamxyz(i,j,k) * gxyx(i,j,k) + Gamyyz(i,j,k) * gxyy(i,j,k) + Gamzyz(i,j,k) * gxyz(i,j,k)  + &
                  Gamxxy(i,j,k) * gyzx(i,j,k) + Gamyxy(i,j,k) * gyzy(i,j,k) + Gamzxy(i,j,k) * gyzz(i,j,k) )+ &
         gupyy(i,j,k) *(                                                  &
             TWO*(Gamxyy(i,j,k) * gxyy(i,j,k) + Gamyyy(i,j,k) * gyyy(i,j,k) + Gamzyy(i,j,k) * gyzy(i,j,k)) + &
                  Gamxyy(i,j,k) * gyyx(i,j,k) + Gamyyy(i,j,k) * gyyy(i,j,k) + Gamzyy(i,j,k) * gyyz(i,j,k) )+ &
         gupyz(i,j,k) *(                                                  &
             TWO*(Gamxyy(i,j,k) * gxzy(i,j,k) + Gamyyy(i,j,k) * gyzy(i,j,k) + Gamzyy(i,j,k) * gzzy(i,j,k)  + &
                  Gamxyz(i,j,k) * gxyy(i,j,k) + Gamyyz(i,j,k) * gyyy(i,j,k) + Gamzyz(i,j,k) * gyzy(i,j,k)) + &
                  Gamxyz(i,j,k) * gyyx(i,j,k) + Gamyyz(i,j,k) * gyyy(i,j,k) + Gamzyz(i,j,k) * gyyz(i,j,k)  + &
                  Gamxyy(i,j,k) * gyzx(i,j,k) + Gamyyy(i,j,k) * gyzy(i,j,k) + Gamzyy(i,j,k) * gyzz(i,j,k) )+ &
         gupzz(i,j,k) *(                                                  &
             TWO*(Gamxyz(i,j,k) * gxzy(i,j,k) + Gamyyz(i,j,k) * gyzy(i,j,k) + Gamzyz(i,j,k) * gzzy(i,j,k)) + &
                  Gamxyz(i,j,k) * gyzx(i,j,k) + Gamyyz(i,j,k) * gyzy(i,j,k) + Gamzyz(i,j,k) * gyzz(i,j,k) )

        Rzz(i,j,k) =     - HALF * Rzz(i,j,k)                                   + &
                   gxz(i,j,k) * Gamxz(i,j,k)+ gyz(i,j,k) * Gamyz(i,j,k)  +    gzz(i,j,k) * Gamzz(i,j,k)  + &
                 Gamxa(i,j,k) * gxzz(i,j,k) +  Gamya(i,j,k) * gyzz(i,j,k) +  Gamza(i,j,k) * gzzz(i,j,k)  + &
         gupxx(i,j,k) *(                                                  &
             TWO*(Gamxxz(i,j,k) * gxxz(i,j,k) + Gamyxz(i,j,k) * gxyz(i,j,k) + Gamzxz(i,j,k) * gxzz(i,j,k)) + &
                  Gamxxz(i,j,k) * gxzx(i,j,k) + Gamyxz(i,j,k) * gxzy(i,j,k) + Gamzxz(i,j,k) * gxzz(i,j,k) )+ &
         gupxy(i,j,k) *(                                                  &
             TWO*(Gamxxz(i,j,k) * gxyz(i,j,k) + Gamyxz(i,j,k) * gyyz(i,j,k) + Gamzxz(i,j,k) * gyzz(i,j,k)  + &
                  Gamxyz(i,j,k) * gxxz(i,j,k) + Gamyyz(i,j,k) * gxyz(i,j,k) + Gamzyz(i,j,k) * gxzz(i,j,k)) + &
                  Gamxyz(i,j,k) * gxzx(i,j,k) + Gamyyz(i,j,k) * gxzy(i,j,k) + Gamzyz(i,j,k) * gxzz(i,j,k)  + &
                  Gamxxz(i,j,k) * gyzx(i,j,k) + Gamyxz(i,j,k) * gyzy(i,j,k) + Gamzxz(i,j,k) * gyzz(i,j,k) )+ &
         gupxz(i,j,k) *(                                                  &
             TWO*(Gamxxz(i,j,k) * gxzz(i,j,k) + Gamyxz(i,j,k) * gyzz(i,j,k) + Gamzxz(i,j,k) * gzzz(i,j,k)  + &
                  Gamxzz(i,j,k) * gxxz(i,j,k) + Gamyzz(i,j,k) * gxyz(i,j,k) + Gamzzz(i,j,k) * gxzz(i,j,k)) + &
                  Gamxzz(i,j,k) * gxzx(i,j,k) + Gamyzz(i,j,k) * gxzy(i,j,k) + Gamzzz(i,j,k) * gxzz(i,j,k)  + &
                  Gamxxz(i,j,k) * gzzx(i,j,k) + Gamyxz(i,j,k) * gzzy(i,j,k) + Gamzxz(i,j,k) * gzzz(i,j,k) )+ &
         gupyy(i,j,k) *(                                                  &
             TWO*(Gamxyz(i,j,k) * gxyz(i,j,k) + Gamyyz(i,j,k) * gyyz(i,j,k) + Gamzyz(i,j,k) * gyzz(i,j,k)) + &
                  Gamxyz(i,j,k) * gyzx(i,j,k) + Gamyyz(i,j,k) * gyzy(i,j,k) + Gamzyz(i,j,k) * gyzz(i,j,k) )+ &
         gupyz(i,j,k) *(                                                  &
             TWO*(Gamxyz(i,j,k) * gxzz(i,j,k) + Gamyyz(i,j,k) * gyzz(i,j,k) + Gamzyz(i,j,k) * gzzz(i,j,k)  + &
                  Gamxzz(i,j,k) * gxyz(i,j,k) + Gamyzz(i,j,k) * gyyz(i,j,k) + Gamzzz(i,j,k) * gyzz(i,j,k)) + &
                  Gamxzz(i,j,k) * gyzx(i,j,k) + Gamyzz(i,j,k) * gyzy(i,j,k) + Gamzzz(i,j,k) * gyzz(i,j,k)  + &
                  Gamxyz(i,j,k) * gzzx(i,j,k) + Gamyyz(i,j,k) * gzzy(i,j,k) + Gamzyz(i,j,k) * gzzz(i,j,k) )+ &
         gupzz(i,j,k) *(                                                  &
             TWO*(Gamxzz(i,j,k) * gxzz(i,j,k) + Gamyzz(i,j,k) * gyzz(i,j,k) + Gamzzz(i,j,k) * gzzz(i,j,k)) + &
                  Gamxzz(i,j,k) * gzzx(i,j,k) + Gamyzz(i,j,k) * gzzy(i,j,k) + Gamzzz(i,j,k) * gzzz(i,j,k) )

        Rxy(i,j,k) = HALF*(     - Rxy(i,j,k)                                   + &
                   gxx(i,j,k) * Gamxy(i,j,k) +    gxy(i,j,k) * Gamyy(i,j,k) + gxz(i,j,k) * Gamzy(i,j,k)  + &
                   gxy(i,j,k) * Gamxx(i,j,k) +    gyy(i,j,k) * Gamyx(i,j,k) + gyz(i,j,k) * Gamzx(i,j,k)  + &
                 Gamxa(i,j,k) * gxyx(i,j,k) +  Gamya(i,j,k) * gyyx(i,j,k) +  Gamza(i,j,k) * gyzx(i,j,k)  + &
                 Gamxa(i,j,k) * gxxy(i,j,k) +  Gamya(i,j,k) * gxyy(i,j,k) +  Gamza(i,j,k) * gxzy(i,j,k) )+ &
         gupxx(i,j,k) *(                                                  &
                  Gamxxx(i,j,k) * gxxy(i,j,k) + Gamyxx(i,j,k) * gxyy(i,j,k) + Gamzxx(i,j,k) * gxzy(i,j,k)  + &
                  Gamxxy(i,j,k) * gxxx(i,j,k) + Gamyxy(i,j,k) * gxyx(i,j,k) + Gamzxy(i,j,k) * gxzx(i,j,k)  + &
                  Gamxxx(i,j,k) * gxyx(i,j,k) + Gamyxx(i,j,k) * gxyy(i,j,k) + Gamzxx(i,j,k) * gxyz(i,j,k) )+ &
         gupxy(i,j,k) *(                                                  &
                  Gamxxx(i,j,k) * gxyy(i,j,k) + Gamyxx(i,j,k) * gyyy(i,j,k) + Gamzxx(i,j,k) * gyzy(i,j,k)  + &
                  Gamxxy(i,j,k) * gxyx(i,j,k) + Gamyxy(i,j,k) * gyyx(i,j,k) + Gamzxy(i,j,k) * gyzx(i,j,k)  + &
                  Gamxxy(i,j,k) * gxyx(i,j,k) + Gamyxy(i,j,k) * gxyy(i,j,k) + Gamzxy(i,j,k) * gxyz(i,j,k)  + &
                  Gamxxy(i,j,k) * gxxy(i,j,k) + Gamyxy(i,j,k) * gxyy(i,j,k) + Gamzxy(i,j,k) * gxzy(i,j,k)  + &
                  Gamxyy(i,j,k) * gxxx(i,j,k) + Gamyyy(i,j,k) * gxyx(i,j,k) + Gamzyy(i,j,k) * gxzx(i,j,k)  + &
                  Gamxxx(i,j,k) * gyyx(i,j,k) + Gamyxx(i,j,k) * gyyy(i,j,k) + Gamzxx(i,j,k) * gyyz(i,j,k) )+ &
         gupxz(i,j,k) *(                                                  &
                  Gamxxx(i,j,k) * gxzy(i,j,k) + Gamyxx(i,j,k) * gyzy(i,j,k) + Gamzxx(i,j,k) * gzzy(i,j,k)  + &
                  Gamxxy(i,j,k) * gxzx(i,j,k) + Gamyxy(i,j,k) * gyzx(i,j,k) + Gamzxy(i,j,k) * gzzx(i,j,k)  + &
                  Gamxxz(i,j,k) * gxyx(i,j,k) + Gamyxz(i,j,k) * gxyy(i,j,k) + Gamzxz(i,j,k) * gxyz(i,j,k)  + &
                  Gamxxz(i,j,k) * gxxy(i,j,k) + Gamyxz(i,j,k) * gxyy(i,j,k) + Gamzxz(i,j,k) * gxzy(i,j,k)  + &
                  Gamxyz(i,j,k) * gxxx(i,j,k) + Gamyyz(i,j,k) * gxyx(i,j,k) + Gamzyz(i,j,k) * gxzx(i,j,k)  + &
                  Gamxxx(i,j,k) * gyzx(i,j,k) + Gamyxx(i,j,k) * gyzy(i,j,k) + Gamzxx(i,j,k) * gyzz(i,j,k) )+ &
         gupyy(i,j,k) *(                                                  &
                  Gamxxy(i,j,k) * gxyy(i,j,k) + Gamyxy(i,j,k) * gyyy(i,j,k) + Gamzxy(i,j,k) * gyzy(i,j,k)  + &
                  Gamxyy(i,j,k) * gxyx(i,j,k) + Gamyyy(i,j,k) * gyyx(i,j,k) + Gamzyy(i,j,k) * gyzx(i,j,k)  + &
                  Gamxxy(i,j,k) * gyyx(i,j,k) + Gamyxy(i,j,k) * gyyy(i,j,k) + Gamzxy(i,j,k) * gyyz(i,j,k) )+ &
         gupyz(i,j,k) *(                                                  &
                  Gamxxy(i,j,k) * gxzy(i,j,k) + Gamyxy(i,j,k) * gyzy(i,j,k) + Gamzxy(i,j,k) * gzzy(i,j,k)  + &
                  Gamxyy(i,j,k) * gxzx(i,j,k) + Gamyyy(i,j,k) * gyzx(i,j,k) + Gamzyy(i,j,k) * gzzx(i,j,k)  + &
                  Gamxxz(i,j,k) * gyyx(i,j,k) + Gamyxz(i,j,k) * gyyy(i,j,k) + Gamzxz(i,j,k) * gyyz(i,j,k)  + &
                  Gamxxz(i,j,k) * gxyy(i,j,k) + Gamyxz(i,j,k) * gyyy(i,j,k) + Gamzxz(i,j,k) * gyzy(i,j,k)  + &
                  Gamxyz(i,j,k) * gxyx(i,j,k) + Gamyyz(i,j,k) * gyyx(i,j,k) + Gamzyz(i,j,k) * gyzx(i,j,k)  + &
                  Gamxxy(i,j,k) * gyzx(i,j,k) + Gamyxy(i,j,k) * gyzy(i,j,k) + Gamzxy(i,j,k) * gyzz(i,j,k) )+ &
         gupzz(i,j,k) *(                                                  &
                  Gamxxz(i,j,k) * gxzy(i,j,k) + Gamyxz(i,j,k) * gyzy(i,j,k) + Gamzxz(i,j,k) * gzzy(i,j,k)  + &
                  Gamxyz(i,j,k) * gxzx(i,j,k) + Gamyyz(i,j,k) * gyzx(i,j,k) + Gamzyz(i,j,k) * gzzx(i,j,k)  + &
                  Gamxxz(i,j,k) * gyzx(i,j,k) + Gamyxz(i,j,k) * gyzy(i,j,k) + Gamzxz(i,j,k) * gyzz(i,j,k) )

        Rxz(i,j,k) = HALF*(     - Rxz(i,j,k)                                   + &
                   gxx(i,j,k) * Gamxz(i,j,k) +  gxy(i,j,k) * Gamyz(i,j,k) + gxz(i,j,k) * Gamzz(i,j,k)    + &
                   gxz(i,j,k) * Gamxx(i,j,k) +  gyz(i,j,k) * Gamyx(i,j,k) + gzz(i,j,k) * Gamzx(i,j,k)    + &
                 Gamxa(i,j,k) * gxzx(i,j,k) +  Gamya(i,j,k) * gyzx(i,j,k) +  Gamza(i,j,k) * gzzx(i,j,k)  + &
                 Gamxa(i,j,k) * gxxz(i,j,k) +  Gamya(i,j,k) * gxyz(i,j,k) +  Gamza(i,j,k) * gxzz(i,j,k) )+ &
         gupxx(i,j,k) *(                                                  &
                  Gamxxx(i,j,k) * gxxz(i,j,k) + Gamyxx(i,j,k) * gxyz(i,j,k) + Gamzxx(i,j,k) * gxzz(i,j,k)  + &
                  Gamxxz(i,j,k) * gxxx(i,j,k) + Gamyxz(i,j,k) * gxyx(i,j,k) + Gamzxz(i,j,k) * gxzx(i,j,k)  + &
                  Gamxxx(i,j,k) * gxzx(i,j,k) + Gamyxx(i,j,k) * gxzy(i,j,k) + Gamzxx(i,j,k) * gxzz(i,j,k) )+ &
         gupxy(i,j,k) *(                                                  &
                  Gamxxx(i,j,k) * gxyz(i,j,k) + Gamyxx(i,j,k) * gyyz(i,j,k) + Gamzxx(i,j,k) * gyzz(i,j,k)  + &
                  Gamxxz(i,j,k) * gxyx(i,j,k) + Gamyxz(i,j,k) * gyyx(i,j,k) + Gamzxz(i,j,k) * gyzx(i,j,k)  + &
                  Gamxxy(i,j,k) * gxzx(i,j,k) + Gamyxy(i,j,k) * gxzy(i,j,k) + Gamzxy(i,j,k) * gxzz(i,j,k)  + &
                  Gamxxy(i,j,k) * gxxz(i,j,k) + Gamyxy(i,j,k) * gxyz(i,j,k) + Gamzxy(i,j,k) * gxzz(i,j,k)  + &
                  Gamxyz(i,j,k) * gxxx(i,j,k) + Gamyyz(i,j,k) * gxyx(i,j,k) + Gamzyz(i,j,k) * gxzx(i,j,k)  + &
                  Gamxxx(i,j,k) * gyzx(i,j,k) + Gamyxx(i,j,k) * gyzy(i,j,k) + Gamzxx(i,j,k) * gyzz(i,j,k) )+ &
         gupxz(i,j,k) *(                                                  &
                  Gamxxx(i,j,k) * gxzz(i,j,k) + Gamyxx(i,j,k) * gyzz(i,j,k) + Gamzxx(i,j,k) * gzzz(i,j,k)  + &
                  Gamxxz(i,j,k) * gxzx(i,j,k) + Gamyxz(i,j,k) * gyzx(i,j,k) + Gamzxz(i,j,k) * gzzx(i,j,k)  + &
                  Gamxxz(i,j,k) * gxzx(i,j,k) + Gamyxz(i,j,k) * gxzy(i,j,k) + Gamzxz(i,j,k) * gxzz(i,j,k)  + &
                  Gamxxz(i,j,k) * gxxz(i,j,k) + Gamyxz(i,j,k) * gxyz(i,j,k) + Gamzxz(i,j,k) * gxzz(i,j,k)  + &
                  Gamxzz(i,j,k) * gxxx(i,j,k) + Gamyzz(i,j,k) * gxyx(i,j,k) + Gamzzz(i,j,k) * gxzx(i,j,k)  + &
                  Gamxxx(i,j,k) * gzzx(i,j,k) + Gamyxx(i,j,k) * gzzy(i,j,k) + Gamzxx(i,j,k) * gzzz(i,j,k) )+ &
         gupyy(i,j,k) *(                                                  &
                  Gamxxy(i,j,k) * gxyz(i,j,k) + Gamyxy(i,j,k) * gyyz(i,j,k) + Gamzxy(i,j,k) * gyzz(i,j,k)  + &
                  Gamxyz(i,j,k) * gxyx(i,j,k) + Gamyyz(i,j,k) * gyyx(i,j,k) + Gamzyz(i,j,k) * gyzx(i,j,k)  + &
                  Gamxxy(i,j,k) * gyzx(i,j,k) + Gamyxy(i,j,k) * gyzy(i,j,k) + Gamzxy(i,j,k) * gyzz(i,j,k) )+ &
         gupyz(i,j,k) *(                                                  &
                  Gamxxy(i,j,k) * gxzz(i,j,k) + Gamyxy(i,j,k) * gyzz(i,j,k) + Gamzxy(i,j,k) * gzzz(i,j,k)  + &
                  Gamxyz(i,j,k) * gxzx(i,j,k) + Gamyyz(i,j,k) * gyzx(i,j,k) + Gamzyz(i,j,k) * gzzx(i,j,k)  + &
                  Gamxxz(i,j,k) * gyzx(i,j,k) + Gamyxz(i,j,k) * gyzy(i,j,k) + Gamzxz(i,j,k) * gyzz(i,j,k)  + &
                  Gamxxz(i,j,k) * gxyz(i,j,k) + Gamyxz(i,j,k) * gyyz(i,j,k) + Gamzxz(i,j,k) * gyzz(i,j,k)  + &
                  Gamxzz(i,j,k) * gxyx(i,j,k) + Gamyzz(i,j,k) * gyyx(i,j,k) + Gamzzz(i,j,k) * gyzx(i,j,k)  + &
                  Gamxxy(i,j,k) * gzzx(i,j,k) + Gamyxy(i,j,k) * gzzy(i,j,k) + Gamzxy(i,j,k) * gzzz(i,j,k) )+ &
         gupzz(i,j,k) *(                                                  &
                  Gamxxz(i,j,k) * gxzz(i,j,k) + Gamyxz(i,j,k) * gyzz(i,j,k) + Gamzxz(i,j,k) * gzzz(i,j,k)  + &
                  Gamxzz(i,j,k) * gxzx(i,j,k) + Gamyzz(i,j,k) * gyzx(i,j,k) + Gamzzz(i,j,k) * gzzx(i,j,k)  + &
                  Gamxxz(i,j,k) * gzzx(i,j,k) + Gamyxz(i,j,k) * gzzy(i,j,k) + Gamzxz(i,j,k) * gzzz(i,j,k) )

        Ryz(i,j,k) = HALF*(     - Ryz(i,j,k)                                   + &
                   gxy(i,j,k) * Gamxz(i,j,k) + gyy(i,j,k) * Gamyz(i,j,k) + gyz(i,j,k) * Gamzz(i,j,k)     + &
                   gxz(i,j,k) * Gamxy(i,j,k) + gyz(i,j,k) * Gamyy(i,j,k) + gzz(i,j,k) * Gamzy(i,j,k)     + &
                 Gamxa(i,j,k) * gxzy(i,j,k) +  Gamya(i,j,k) * gyzy(i,j,k) +  Gamza(i,j,k) * gzzy(i,j,k)  + &
                 Gamxa(i,j,k) * gxyz(i,j,k) +  Gamya(i,j,k) * gyyz(i,j,k) +  Gamza(i,j,k) * gyzz(i,j,k) )+ &
         gupxx(i,j,k) *(                                                  &
                  Gamxxy(i,j,k) * gxxz(i,j,k) + Gamyxy(i,j,k) * gxyz(i,j,k) + Gamzxy(i,j,k) * gxzz(i,j,k)  + &
                  Gamxxz(i,j,k) * gxxy(i,j,k) + Gamyxz(i,j,k) * gxyy(i,j,k) + Gamzxz(i,j,k) * gxzy(i,j,k)  + &
                  Gamxxy(i,j,k) * gxzx(i,j,k) + Gamyxy(i,j,k) * gxzy(i,j,k) + Gamzxy(i,j,k) * gxzz(i,j,k) )+ &
         gupxy(i,j,k) *(                                                  &
                  Gamxxy(i,j,k) * gxyz(i,j,k) + Gamyxy(i,j,k) * gyyz(i,j,k) + Gamzxy(i,j,k) * gyzz(i,j,k)  + &
                  Gamxxz(i,j,k) * gxyy(i,j,k) + Gamyxz(i,j,k) * gyyy(i,j,k) + Gamzxz(i,j,k) * gyzy(i,j,k)  + &
                  Gamxyy(i,j,k) * gxzx(i,j,k) + Gamyyy(i,j,k) * gxzy(i,j,k) + Gamzyy(i,j,k) * gxzz(i,j,k)  + &
                  Gamxyy(i,j,k) * gxxz(i,j,k) + Gamyyy(i,j,k) * gxyz(i,j,k) + Gamzyy(i,j,k) * gxzz(i,j,k)  + &
                  Gamxyz(i,j,k) * gxxy(i,j,k) + Gamyyz(i,j,k) * gxyy(i,j,k) + Gamzyz(i,j,k) * gxzy(i,j,k)  + &
                  Gamxxy(i,j,k) * gyzx(i,j,k) + Gamyxy(i,j,k) * gyzy(i,j,k) + Gamzxy(i,j,k) * gyzz(i,j,k) )+ &
         gupxz(i,j,k) *(                                                  &
                  Gamxxy(i,j,k) * gxzz(i,j,k) + Gamyxy(i,j,k) * gyzz(i,j,k) + Gamzxy(i,j,k) * gzzz(i,j,k)  + &
                  Gamxxz(i,j,k) * gxzy(i,j,k) + Gamyxz(i,j,k) * gyzy(i,j,k) + Gamzxz(i,j,k) * gzzy(i,j,k)  + &
                  Gamxyz(i,j,k) * gxzx(i,j,k) + Gamyyz(i,j,k) * gxzy(i,j,k) + Gamzyz(i,j,k) * gxzz(i,j,k)  + &
                  Gamxyz(i,j,k) * gxxz(i,j,k) + Gamyyz(i,j,k) * gxyz(i,j,k) + Gamzyz(i,j,k) * gxzz(i,j,k)  + &
                  Gamxzz(i,j,k) * gxxy(i,j,k) + Gamyzz(i,j,k) * gxyy(i,j,k) + Gamzzz(i,j,k) * gxzy(i,j,k)  + &
                  Gamxxy(i,j,k) * gzzx(i,j,k) + Gamyxy(i,j,k) * gzzy(i,j,k) + Gamzxy(i,j,k) * gzzz(i,j,k) )+ &
         gupyy(i,j,k) *(                                                  &
                  Gamxyy(i,j,k) * gxyz(i,j,k) + Gamyyy(i,j,k) * gyyz(i,j,k) + Gamzyy(i,j,k) * gyzz(i,j,k)  + &
                  Gamxyz(i,j,k) * gxyy(i,j,k) + Gamyyz(i,j,k) * gyyy(i,j,k) + Gamzyz(i,j,k) * gyzy(i,j,k)  + &
                  Gamxyy(i,j,k) * gyzx(i,j,k) + Gamyyy(i,j,k) * gyzy(i,j,k) + Gamzyy(i,j,k) * gyzz(i,j,k) )+ &
         gupyz(i,j,k) *(                                                  &
                  Gamxyy(i,j,k) * gxzz(i,j,k) + Gamyyy(i,j,k) * gyzz(i,j,k) + Gamzyy(i,j,k) * gzzz(i,j,k)  + &
                  Gamxyz(i,j,k) * gxzy(i,j,k) + Gamyyz(i,j,k) * gyzy(i,j,k) + Gamzyz(i,j,k) * gzzy(i,j,k)  + &
                  Gamxyz(i,j,k) * gyzx(i,j,k) + Gamyyz(i,j,k) * gyzy(i,j,k) + Gamzyz(i,j,k) * gyzz(i,j,k)  + &
                  Gamxyz(i,j,k) * gxyz(i,j,k) + Gamyyz(i,j,k) * gyyz(i,j,k) + Gamzyz(i,j,k) * gyzz(i,j,k)  + &
                  Gamxzz(i,j,k) * gxyy(i,j,k) + Gamyzz(i,j,k) * gyyy(i,j,k) + Gamzzz(i,j,k) * gyzy(i,j,k)  + &
                  Gamxyy(i,j,k) * gzzx(i,j,k) + Gamyyy(i,j,k) * gzzy(i,j,k) + Gamzyy(i,j,k) * gzzz(i,j,k) )+ &
         gupzz(i,j,k) *(                                                  &
                  Gamxyz(i,j,k) * gxzz(i,j,k) + Gamyyz(i,j,k) * gyzz(i,j,k) + Gamzyz(i,j,k) * gzzz(i,j,k)  + &
                  Gamxzz(i,j,k) * gxzy(i,j,k) + Gamyzz(i,j,k) * gyzy(i,j,k) + Gamzzz(i,j,k) * gzzy(i,j,k)  + &
                  Gamxyz(i,j,k) * gzzx(i,j,k) + Gamyyz(i,j,k) * gzzy(i,j,k) + Gamzyz(i,j,k) * gzzz(i,j,k) )
      end do

    end do
  end do
!$OMP END DO NOWAIT
  call sft_trace("Rxx", ts0_thr, sft_get_ts(), 0, tid)
  ts0_thr = sft_get_ts()
!covariant second derivative of chi respect to tilted metric
! Batch chi + Lap sym_bd together
    call symmetry_bd_inner(OPT_SYMORD,ex,chi,fh_d,SSS)
    call symmetry_bd_inner(OPT_SYMORD,ex,Lap,fh_d2,SSS)
!$OMP BARRIER
    call fdderivs_inner(ex,fh_d,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,Symmetry)
!$OMP BARRIER
!$OMP DO COLLAPSE(1)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        fxx(i,j,k) = fxx(i,j,k) - Gamxxx(i,j,k) * chix(i,j,k) - Gamyxx(i,j,k) * chiy(i,j,k) - Gamzxx(i,j,k) * chiz(i,j,k)
        fxy(i,j,k) = fxy(i,j,k) - Gamxxy(i,j,k) * chix(i,j,k) - Gamyxy(i,j,k) * chiy(i,j,k) - Gamzxy(i,j,k) * chiz(i,j,k)
        fxz(i,j,k) = fxz(i,j,k) - Gamxxz(i,j,k) * chix(i,j,k) - Gamyxz(i,j,k) * chiy(i,j,k) - Gamzxz(i,j,k) * chiz(i,j,k)
        fyy(i,j,k) = fyy(i,j,k) - Gamxyy(i,j,k) * chix(i,j,k) - Gamyyy(i,j,k) * chiy(i,j,k) - Gamzyy(i,j,k) * chiz(i,j,k)
        fyz(i,j,k) = fyz(i,j,k) - Gamxyz(i,j,k) * chix(i,j,k) - Gamyyz(i,j,k) * chiy(i,j,k) - Gamzyz(i,j,k) * chiz(i,j,k)
        fzz(i,j,k) = fzz(i,j,k) - Gamxzz(i,j,k) * chix(i,j,k) - Gamyzz(i,j,k) * chiy(i,j,k) - Gamzzz(i,j,k) * chiz(i,j,k)
! Store D^l D_l chi - 3/(2*chi) D^l chi D_l chi in f

        f(i,j,k) =        gupxx(i,j,k) * ( fxx(i,j,k) - F3o2/chin1(i,j,k) * chix(i,j,k) * chix(i,j,k) ) + &
                           gupyy(i,j,k) * ( fyy(i,j,k) - F3o2/chin1(i,j,k) * chiy(i,j,k) * chiy(i,j,k) ) + &
                           gupzz(i,j,k) * ( fzz(i,j,k) - F3o2/chin1(i,j,k) * chiz(i,j,k) * chiz(i,j,k) ) + &
                     TWO * gupxy(i,j,k) * ( fxy(i,j,k) - F3o2/chin1(i,j,k) * chix(i,j,k) * chiy(i,j,k) ) + &
                     TWO * gupxz(i,j,k) * ( fxz(i,j,k) - F3o2/chin1(i,j,k) * chix(i,j,k) * chiz(i,j,k) ) + &
                     TWO * gupyz(i,j,k) * ( fyz(i,j,k) - F3o2/chin1(i,j,k) * chiy(i,j,k) * chiz(i,j,k) )
! Add chi part to Ricci tensor:

        Rxx(i,j,k) = Rxx(i,j,k) + (fxx(i,j,k) - chix(i,j,k)*chix(i,j,k)/chin1(i,j,k)/TWO + gxx(i,j,k) * f(i,j,k))/chin1(i,j,k)/TWO
        Ryy(i,j,k) = Ryy(i,j,k) + (fyy(i,j,k) - chiy(i,j,k)*chiy(i,j,k)/chin1(i,j,k)/TWO + gyy(i,j,k) * f(i,j,k))/chin1(i,j,k)/TWO
        Rzz(i,j,k) = Rzz(i,j,k) + (fzz(i,j,k) - chiz(i,j,k)*chiz(i,j,k)/chin1(i,j,k)/TWO + gzz(i,j,k) * f(i,j,k))/chin1(i,j,k)/TWO
        Rxy(i,j,k) = Rxy(i,j,k) + (fxy(i,j,k) - chix(i,j,k)*chiy(i,j,k)/chin1(i,j,k)/TWO + gxy(i,j,k) * f(i,j,k))/chin1(i,j,k)/TWO
        Rxz(i,j,k) = Rxz(i,j,k) + (fxz(i,j,k) - chix(i,j,k)*chiz(i,j,k)/chin1(i,j,k)/TWO + gxz(i,j,k) * f(i,j,k))/chin1(i,j,k)/TWO
        Ryz(i,j,k) = Ryz(i,j,k) + (fyz(i,j,k) - chiy(i,j,k)*chiz(i,j,k)/chin1(i,j,k)/TWO + gyz(i,j,k) * f(i,j,k))/chin1(i,j,k)/TWO
      end do
    end do
  end do
!$OMP END DO NOWAIT
  call sft_trace("fxx", ts0_thr, sft_get_ts(), 0, tid)
  ts0_thr = sft_get_ts()

! covariant second derivatives of the lapse — Lap fdderivs overlaps with chi DO above
    call fdderivs_inner(ex,fh_d2,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,symmetry)
!$OMP BARRIER
!$OMP DO COLLAPSE(1)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        gxxx(i,j,k) = (gupxx(i,j,k) * chix(i,j,k) + gupxy(i,j,k) * chiy(i,j,k) + gupxz(i,j,k) * chiz(i,j,k))/chin1(i,j,k)
        gxxy(i,j,k) = (gupxy(i,j,k) * chix(i,j,k) + gupyy(i,j,k) * chiy(i,j,k) + gupyz(i,j,k) * chiz(i,j,k))/chin1(i,j,k)
        gxxz(i,j,k) = (gupxz(i,j,k) * chix(i,j,k) + gupyz(i,j,k) * chiy(i,j,k) + gupzz(i,j,k) * chiz(i,j,k))/chin1(i,j,k)
! now get physical second kind of connection
        Gamxxx(i,j,k) = Gamxxx(i,j,k) - ( (chix(i,j,k) + chix(i,j,k))/chin1(i,j,k) - gxx(i,j,k) * gxxx(i,j,k) )*HALF
        Gamyxx(i,j,k) = Gamyxx(i,j,k) - (                                             - gxx(i,j,k) * gxxy(i,j,k) )*HALF
        Gamzxx(i,j,k) = Gamzxx(i,j,k) - (                                             - gxx(i,j,k) * gxxz(i,j,k) )*HALF
        Gamxyy(i,j,k) = Gamxyy(i,j,k) - (                                             - gyy(i,j,k) * gxxx(i,j,k) )*HALF
        Gamyyy(i,j,k) = Gamyyy(i,j,k) - ( (chiy(i,j,k) + chiy(i,j,k))/chin1(i,j,k) - gyy(i,j,k) * gxxy(i,j,k) )*HALF
        Gamzyy(i,j,k) = Gamzyy(i,j,k) - (                                             - gyy(i,j,k) * gxxz(i,j,k) )*HALF
        Gamxzz(i,j,k) = Gamxzz(i,j,k) - (                                             - gzz(i,j,k) * gxxx(i,j,k) )*HALF
        Gamyzz(i,j,k) = Gamyzz(i,j,k) - (                                             - gzz(i,j,k) * gxxy(i,j,k) )*HALF
        Gamzzz(i,j,k) = Gamzzz(i,j,k) - ( (chiz(i,j,k) + chiz(i,j,k))/chin1(i,j,k) - gzz(i,j,k) * gxxz(i,j,k) )*HALF
        Gamxxy(i,j,k) = Gamxxy(i,j,k) - (  chiy(i,j,k)        /chin1(i,j,k)         - gxy(i,j,k) * gxxx(i,j,k) )*HALF
        Gamyxy(i,j,k) = Gamyxy(i,j,k) - (                chix(i,j,k) /chin1(i,j,k)  - gxy(i,j,k) * gxxy(i,j,k) )*HALF
        Gamzxy(i,j,k) = Gamzxy(i,j,k) - (                                             - gxy(i,j,k) * gxxz(i,j,k) )*HALF
        Gamxxz(i,j,k) = Gamxxz(i,j,k) - (  chiz(i,j,k)        /chin1(i,j,k)         - gxz(i,j,k) * gxxx(i,j,k) )*HALF
        Gamyxz(i,j,k) = Gamyxz(i,j,k) - (                                             - gxz(i,j,k) * gxxy(i,j,k) )*HALF
        Gamzxz(i,j,k) = Gamzxz(i,j,k) - (                chix(i,j,k) /chin1(i,j,k)  - gxz(i,j,k) * gxxz(i,j,k) )*HALF
        Gamxyz(i,j,k) = Gamxyz(i,j,k) - (                                             - gyz(i,j,k) * gxxx(i,j,k) )*HALF
        Gamyyz(i,j,k) = Gamyyz(i,j,k) - (  chiz(i,j,k)        /chin1(i,j,k)         - gyz(i,j,k) * gxxy(i,j,k) )*HALF
        Gamzyz(i,j,k) = Gamzyz(i,j,k) - (                chiy(i,j,k) /chin1(i,j,k)  - gyz(i,j,k) * gxxz(i,j,k) )*HALF

        fxx(i,j,k) = fxx(i,j,k) - Gamxxx(i,j,k)*Lapx(i,j,k) - Gamyxx(i,j,k)*Lapy(i,j,k) - Gamzxx(i,j,k)*Lapz(i,j,k)
        fyy(i,j,k) = fyy(i,j,k) - Gamxyy(i,j,k)*Lapx(i,j,k) - Gamyyy(i,j,k)*Lapy(i,j,k) - Gamzyy(i,j,k)*Lapz(i,j,k)
        fzz(i,j,k) = fzz(i,j,k) - Gamxzz(i,j,k)*Lapx(i,j,k) - Gamyzz(i,j,k)*Lapy(i,j,k) - Gamzzz(i,j,k)*Lapz(i,j,k)
        fxy(i,j,k) = fxy(i,j,k) - Gamxxy(i,j,k)*Lapx(i,j,k) - Gamyxy(i,j,k)*Lapy(i,j,k) - Gamzxy(i,j,k)*Lapz(i,j,k)
        fxz(i,j,k) = fxz(i,j,k) - Gamxxz(i,j,k)*Lapx(i,j,k) - Gamyxz(i,j,k)*Lapy(i,j,k) - Gamzxz(i,j,k)*Lapz(i,j,k)
        fyz(i,j,k) = fyz(i,j,k) - Gamxyz(i,j,k)*Lapx(i,j,k) - Gamyyz(i,j,k)*Lapy(i,j,k) - Gamzyz(i,j,k)*Lapz(i,j,k)

! store D^i D_i Lap in trK_rhs upto chi
        trK_rhs(i,j,k) =    gupxx(i,j,k) * fxx(i,j,k) + gupyy(i,j,k) * fyy(i,j,k) + gupzz(i,j,k) * fzz(i,j,k) + &
               TWO* ( gupxy(i,j,k) * fxy(i,j,k) + gupxz(i,j,k) * fxz(i,j,k) + gupyz(i,j,k) * fyz(i,j,k) )
#if 1
!! follow bam code
        S(i,j,k) =  chin1(i,j,k) * ( gupxx(i,j,k) * Sxx(i,j,k) + gupyy(i,j,k) * Syy(i,j,k) + gupzz(i,j,k) * Szz(i,j,k) + &
           TWO * ( gupxy(i,j,k) * Sxy(i,j,k) + gupxz(i,j,k) * Sxz(i,j,k) + gupyz(i,j,k) * Syz(i,j,k) ) )
        f(i,j,k) = F2o3 * trK(i,j,k) * trK(i,j,k) -(&
             gupxx(i,j,k) * ( &
             gupxx(i,j,k) * Axx(i,j,k) * Axx(i,j,k) + gupyy(i,j,k) * Axy(i,j,k) * Axy(i,j,k) + gupzz(i,j,k) * Axz(i,j,k) * Axz(i,j,k) + &
             TWO * (gupxy(i,j,k) * Axx(i,j,k) * Axy(i,j,k) + gupxz(i,j,k) * Axx(i,j,k) * Axz(i,j,k) + gupyz(i,j,k) * Axy(i,j,k) * Axz(i,j,k)) ) + &
             gupyy(i,j,k) * ( &
             gupxx(i,j,k) * Axy(i,j,k) * Axy(i,j,k) + gupyy(i,j,k) * Ayy(i,j,k) * Ayy(i,j,k) + gupzz(i,j,k) * Ayz(i,j,k) * Ayz(i,j,k) + &
             TWO * (gupxy(i,j,k) * Axy(i,j,k) * Ayy(i,j,k) + gupxz(i,j,k) * Axy(i,j,k) * Ayz(i,j,k) + gupyz(i,j,k) * Ayy(i,j,k) * Ayz(i,j,k)) ) + &
             gupzz(i,j,k) * ( &
             gupxx(i,j,k) * Axz(i,j,k) * Axz(i,j,k) + gupyy(i,j,k) * Ayz(i,j,k) * Ayz(i,j,k) + gupzz(i,j,k) * Azz(i,j,k) * Azz(i,j,k) + &
             TWO * (gupxy(i,j,k) * Axz(i,j,k) * Ayz(i,j,k) + gupxz(i,j,k) * Axz(i,j,k) * Azz(i,j,k) + gupyz(i,j,k) * Ayz(i,j,k) * Azz(i,j,k)) ) + &
             TWO * ( &
             gupxy(i,j,k) * ( &
             gupxx(i,j,k) * Axx(i,j,k) * Axy(i,j,k) + gupyy(i,j,k) * Axy(i,j,k) * Ayy(i,j,k) + gupzz(i,j,k) * Axz(i,j,k) * Ayz(i,j,k) + &
             gupxy(i,j,k) * (Axx(i,j,k) * Ayy(i,j,k) + Axy(i,j,k) * Axy(i,j,k)) + &
             gupxz(i,j,k) * (Axx(i,j,k) * Ayz(i,j,k) + Axz(i,j,k) * Axy(i,j,k)) + &
             gupyz(i,j,k) * (Axy(i,j,k) * Ayz(i,j,k) + Axz(i,j,k) * Ayy(i,j,k)) ) + &
             gupxz(i,j,k) * ( &
             gupxx(i,j,k) * Axx(i,j,k) * Axz(i,j,k) + gupyy(i,j,k) * Axy(i,j,k) * Ayz(i,j,k) + gupzz(i,j,k) * Axz(i,j,k) * Azz(i,j,k) + &
             gupxy(i,j,k) * (Axx(i,j,k) * Ayz(i,j,k) + Axy(i,j,k) * Axz(i,j,k)) + &
             gupxz(i,j,k) * (Axx(i,j,k) * Azz(i,j,k) + Axz(i,j,k) * Axz(i,j,k)) + &
             gupyz(i,j,k) * (Axy(i,j,k) * Azz(i,j,k) + Axz(i,j,k) * Ayz(i,j,k)) ) + &
             gupyz(i,j,k) * ( &
             gupxx(i,j,k) * Axy(i,j,k) * Axz(i,j,k) + gupyy(i,j,k) * Ayy(i,j,k) * Ayz(i,j,k) + gupzz(i,j,k) * Ayz(i,j,k) * Azz(i,j,k) + &
             gupxy(i,j,k) * (Axy(i,j,k) * Ayz(i,j,k) + Ayy(i,j,k) * Axz(i,j,k)) + &
             gupxz(i,j,k) * (Axy(i,j,k) * Azz(i,j,k) + Ayz(i,j,k) * Axz(i,j,k)) + &
             gupyz(i,j,k) * (Ayy(i,j,k) * Azz(i,j,k) + Ayz(i,j,k) * Ayz(i,j,k)) ) )) -1.6d1*PI*rho(i,j,k) + EIGHT * PI * S(i,j,k)
        f(i,j,k) = - F1o3 *(  gupxx(i,j,k) * fxx(i,j,k) + gupyy(i,j,k) * fyy(i,j,k) + gupzz(i,j,k) * fzz(i,j,k) + &
               TWO* ( gupxy(i,j,k) * fxy(i,j,k) + gupxz(i,j,k) * fxz(i,j,k) + gupyz(i,j,k) * fyz(i,j,k) ) + alpn1(i,j,k)/chin1(i,j,k)*f(i,j,k))

        fxx(i,j,k) = alpn1(i,j,k) * (Rxx(i,j,k) - EIGHT * PI * Sxx(i,j,k)) - fxx(i,j,k)
        fxy(i,j,k) = alpn1(i,j,k) * (Rxy(i,j,k) - EIGHT * PI * Sxy(i,j,k)) - fxy(i,j,k)
        fxz(i,j,k) = alpn1(i,j,k) * (Rxz(i,j,k) - EIGHT * PI * Sxz(i,j,k)) - fxz(i,j,k)
        fyy(i,j,k) = alpn1(i,j,k) * (Ryy(i,j,k) - EIGHT * PI * Syy(i,j,k)) - fyy(i,j,k)
        fyz(i,j,k) = alpn1(i,j,k) * (Ryz(i,j,k) - EIGHT * PI * Syz(i,j,k)) - fyz(i,j,k)
        fzz(i,j,k) = alpn1(i,j,k) * (Rzz(i,j,k) - EIGHT * PI * Szz(i,j,k)) - fzz(i,j,k)
#else
! Add lapse and S_ij parts to Ricci tensor:

        fxx(i,j,k) = alpn1(i,j,k) * (Rxx(i,j,k) - EIGHT * PI * Sxx(i,j,k)) - fxx(i,j,k)
        fxy(i,j,k) = alpn1(i,j,k) * (Rxy(i,j,k) - EIGHT * PI * Sxy(i,j,k)) - fxy(i,j,k)
        fxz(i,j,k) = alpn1(i,j,k) * (Rxz(i,j,k) - EIGHT * PI * Sxz(i,j,k)) - fxz(i,j,k)
        fyy(i,j,k) = alpn1(i,j,k) * (Ryy(i,j,k) - EIGHT * PI * Syy(i,j,k)) - fyy(i,j,k)
        fyz(i,j,k) = alpn1(i,j,k) * (Ryz(i,j,k) - EIGHT * PI * Syz(i,j,k)) - fyz(i,j,k)
        fzz(i,j,k) = alpn1(i,j,k) * (Rzz(i,j,k) - EIGHT * PI * Szz(i,j,k)) - fzz(i,j,k)

! Compute trace-free part (note: chi^-1 and chi cancel!):

        f(i,j,k) = F1o3 *(  gupxx(i,j,k) * fxx(i,j,k) + gupyy(i,j,k) * fyy(i,j,k) + gupzz(i,j,k) * fzz(i,j,k) + &
               TWO* ( gupxy(i,j,k) * fxy(i,j,k) + gupxz(i,j,k) * fxz(i,j,k) + gupyz(i,j,k) * fyz(i,j,k) ) )
#endif
      end do
    end do
  end do
!$OMP END DO NOWAIT
  call sft_trace("gxxx", ts0_thr, sft_get_ts(), 0, tid)
  ts0_thr = sft_get_ts()

!$OMP BARRIER
!$OMP DO COLLAPSE(1)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
! Axx_rhs = trace-free projection — reads old fij and f, writes Aij_rhs
        Axx_rhs(i,j,k) = fxx(i,j,k) - gxx(i,j,k) * f(i,j,k)
        Ayy_rhs(i,j,k) = fyy(i,j,k) - gyy(i,j,k) * f(i,j,k)
        Azz_rhs(i,j,k) = fzz(i,j,k) - gzz(i,j,k) * f(i,j,k)
        Axy_rhs(i,j,k) = fxy(i,j,k) - gxy(i,j,k) * f(i,j,k)
        Axz_rhs(i,j,k) = fxz(i,j,k) - gxz(i,j,k) * f(i,j,k)
        Ayz_rhs(i,j,k) = fyz(i,j,k) - gyz(i,j,k) * f(i,j,k)
! Now: store A_il A^l_j into fij — overwrites fij
        fxx(i,j,k) =       gupxx(i,j,k) * Axx(i,j,k) * Axx(i,j,k) + gupyy(i,j,k) * Axy(i,j,k) * Axy(i,j,k) + gupzz(i,j,k) * Axz(i,j,k) * Axz(i,j,k) + &
             TWO * (gupxy(i,j,k) * Axx(i,j,k) * Axy(i,j,k) + gupxz(i,j,k) * Axx(i,j,k) * Axz(i,j,k) + gupyz(i,j,k) * Axy(i,j,k) * Axz(i,j,k))
        fyy(i,j,k) =       gupxx(i,j,k) * Axy(i,j,k) * Axy(i,j,k) + gupyy(i,j,k) * Ayy(i,j,k) * Ayy(i,j,k) + gupzz(i,j,k) * Ayz(i,j,k) * Ayz(i,j,k) + &
             TWO * (gupxy(i,j,k) * Axy(i,j,k) * Ayy(i,j,k) + gupxz(i,j,k) * Axy(i,j,k) * Ayz(i,j,k) + gupyz(i,j,k) * Ayy(i,j,k) * Ayz(i,j,k))
        fzz(i,j,k) =       gupxx(i,j,k) * Axz(i,j,k) * Axz(i,j,k) + gupyy(i,j,k) * Ayz(i,j,k) * Ayz(i,j,k) + gupzz(i,j,k) * Azz(i,j,k) * Azz(i,j,k) + &
             TWO * (gupxy(i,j,k) * Axz(i,j,k) * Ayz(i,j,k) + gupxz(i,j,k) * Axz(i,j,k) * Azz(i,j,k) + gupyz(i,j,k) * Ayz(i,j,k) * Azz(i,j,k))
        fxy(i,j,k) =       gupxx(i,j,k) * Axx(i,j,k) * Axy(i,j,k) + gupyy(i,j,k) * Axy(i,j,k) * Ayy(i,j,k) + gupzz(i,j,k) * Axz(i,j,k) * Ayz(i,j,k) + &
                    gupxy(i,j,k) *(Axx(i,j,k) * Ayy(i,j,k) + Axy(i,j,k) * Axy(i,j,k))                            + &
                    gupxz(i,j,k) *(Axx(i,j,k) * Ayz(i,j,k) + Axz(i,j,k) * Axy(i,j,k))                            + &
                    gupyz(i,j,k) *(Axy(i,j,k) * Ayz(i,j,k) + Axz(i,j,k) * Ayy(i,j,k))
        fxz(i,j,k) =       gupxx(i,j,k) * Axx(i,j,k) * Axz(i,j,k) + gupyy(i,j,k) * Axy(i,j,k) * Ayz(i,j,k) + gupzz(i,j,k) * Axz(i,j,k) * Azz(i,j,k) + &
                    gupxy(i,j,k) *(Axx(i,j,k) * Ayz(i,j,k) + Axy(i,j,k) * Axz(i,j,k))                            + &
                    gupxz(i,j,k) *(Axx(i,j,k) * Azz(i,j,k) + Axz(i,j,k) * Axz(i,j,k))                            + &
                    gupyz(i,j,k) *(Axy(i,j,k) * Azz(i,j,k) + Axz(i,j,k) * Ayz(i,j,k))
        fyz(i,j,k) =       gupxx(i,j,k) * Axy(i,j,k) * Axz(i,j,k) + gupyy(i,j,k) * Ayy(i,j,k) * Ayz(i,j,k) + gupzz(i,j,k) * Ayz(i,j,k) * Azz(i,j,k) + &
                    gupxy(i,j,k) *(Axy(i,j,k) * Ayz(i,j,k) + Ayy(i,j,k) * Axz(i,j,k))                            + &
                    gupxz(i,j,k) *(Axy(i,j,k) * Azz(i,j,k) + Ayz(i,j,k) * Axz(i,j,k))                            + &
                    gupyz(i,j,k) *(Ayy(i,j,k) * Azz(i,j,k) + Ayz(i,j,k) * Ayz(i,j,k))

        f(i,j,k) = chin1(i,j,k)
! store D^i D_i Lap in trK_rhs
        trK_rhs(i,j,k) = f(i,j,k)*trK_rhs(i,j,k)
      end do
    end do
  end do
!$OMP END DO NOWAIT
  call sft_trace("Axx_rhs_fxx", ts0_thr, sft_get_ts(), 0, tid)
  ts0_thr = sft_get_ts()

!$OMP BARRIER
!$OMP DO COLLAPSE(1)
  do k = 1, ex(3)
    do j = 1, ex(2)

      do i = 1, ex(1)
        Axx_rhs(i,j,k) =           f(i,j,k) * Axx_rhs(i,j,k)+ alpn1(i,j,k) * (trK(i,j,k) * Axx(i,j,k) - TWO * fxx(i,j,k))  + &
                 TWO * (  Axx(i,j,k) * betaxx(i,j,k) +   Axy(i,j,k) * betayx(i,j,k) +   Axz(i,j,k) * betazx(i,j,k) )- &
                   F2o3 * Axx(i,j,k) * div_beta(i,j,k)

        Ayy_rhs(i,j,k) =           f(i,j,k) * Ayy_rhs(i,j,k)+ alpn1(i,j,k) * (trK(i,j,k) * Ayy(i,j,k) - TWO * fyy(i,j,k))  + &
                 TWO * (  Axy(i,j,k) * betaxy(i,j,k) +   Ayy(i,j,k) * betayy(i,j,k) +   Ayz(i,j,k) * betazy(i,j,k) )- &
                   F2o3 * Ayy(i,j,k) * div_beta(i,j,k)

        Azz_rhs(i,j,k) =           f(i,j,k) * Azz_rhs(i,j,k)+ alpn1(i,j,k) * (trK(i,j,k) * Azz(i,j,k) - TWO * fzz(i,j,k))  + &
                 TWO * (  Axz(i,j,k) * betaxz(i,j,k) +   Ayz(i,j,k) * betayz(i,j,k) +   Azz(i,j,k) * betazz(i,j,k) )- &
                   F2o3 * Azz(i,j,k) * div_beta(i,j,k)

        Axy_rhs(i,j,k) =           f(i,j,k) * Axy_rhs(i,j,k)+ alpn1(i,j,k) *( trK(i,j,k) * Axy(i,j,k)  - TWO * fxy(i,j,k) )+ &
                          Axx(i,j,k) * betaxy(i,j,k)                  +   Axz(i,j,k) * betazy(i,j,k)  + &
                                               Ayy(i,j,k) * betayx(i,j,k) +   Ayz(i,j,k) * betazx(i,j,k)  + &
                   F1o3 * Axy(i,j,k) * div_beta(i,j,k)                -   Axy(i,j,k) * betazz(i,j,k)

        Ayz_rhs(i,j,k) =           f(i,j,k) * Ayz_rhs(i,j,k)+ alpn1(i,j,k) *( trK(i,j,k) * Ayz(i,j,k)  - TWO * fyz(i,j,k) )+ &
                          Axy(i,j,k) * betaxz(i,j,k) +   Ayy(i,j,k) * betayz(i,j,k)                   + &
                          Axz(i,j,k) * betaxy(i,j,k)                  +   Azz(i,j,k) * betazy(i,j,k)  + &
                   F1o3 * Ayz(i,j,k) * div_beta(i,j,k)                -   Ayz(i,j,k) * betaxx(i,j,k)
 
        Axz_rhs(i,j,k) =           f(i,j,k) * Axz_rhs(i,j,k)+ alpn1(i,j,k) *( trK(i,j,k) * Axz(i,j,k)  - TWO * fxz(i,j,k) )+ &
                          Axx(i,j,k) * betaxz(i,j,k) +   Axy(i,j,k) * betayz(i,j,k)                   + &
                                               Ayz(i,j,k) * betayx(i,j,k) +   Azz(i,j,k) * betazx(i,j,k)  + &
                   F1o3 * Axz(i,j,k) * div_beta(i,j,k)                -   Axz(i,j,k) * betayy(i,j,k)      !rhs for Aij
! Compute trace of S_ij (merged)
        S(i,j,k) =  f(i,j,k) * ( gupxx(i,j,k) * Sxx(i,j,k) + gupyy(i,j,k) * Syy(i,j,k) + gupzz(i,j,k) * Szz(i,j,k) + &
           TWO * ( gupxy(i,j,k) * Sxy(i,j,k) + gupxz(i,j,k) * Sxz(i,j,k) + gupyz(i,j,k) * Syz(i,j,k) ) )
      end do

    end do
  end do
!$OMP END DO NOWAIT
  call sft_trace("Axx_rhs_S", ts0_thr, sft_get_ts(), 0, tid)
  ts0_thr = sft_get_ts()

!$OMP BARRIER
!$OMP DO COLLAPSE(1)
  do k = 1, ex(3)
    do j = 1, ex(2)

      do i = 1, ex(1)
        trK_rhs(i,j,k) = - trK_rhs(i,j,k) + alpn1(i,j,k) *( F1o3 * trK(i,j,k) * trK(i,j,k)         + &
                      gupxx(i,j,k) * fxx(i,j,k) + gupyy(i,j,k) * fyy(i,j,k) + gupzz(i,j,k) * fzz(i,j,k)   + &
              TWO * ( gupxy(i,j,k) * fxy(i,j,k) + gupxz(i,j,k) * fxz(i,j,k) + gupyz(i,j,k) * fyz(i,j,k) ) + &
             FOUR * PI * ( rho(i,j,k) + S(i,j,k) ))                                !rhs for trK
        Lap_rhs(i,j,k) = -TWO*alpn1(i,j,k)*trK(i,j,k)                             !rhs for Lap (merged)
      end do

    end do
  end do
!$OMP END DO NOWAIT
  call sft_trace("trK_Lap_rhs", ts0_thr, sft_get_ts(), 0, tid)
  ts0_thr = sft_get_ts()

!!!! gauge variable part
#if (GAUGE == 0)
!$OMP BARRIER
!$OMP DO COLLAPSE(1)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
       betax_rhs(i,j,k) = FF*dtSfx(i,j,k)
        betay_rhs(i,j,k) = FF*dtSfy(i,j,k)
        betaz_rhs(i,j,k) = FF*dtSfz(i,j,k)

        dtSfx_rhs(i,j,k) = Gamx_rhs(i,j,k) - eta*dtSfx(i,j,k)
        dtSfy_rhs(i,j,k) = Gamy_rhs(i,j,k) - eta*dtSfy(i,j,k)
        dtSfz_rhs(i,j,k) = Gamz_rhs(i,j,k) - eta*dtSfz(i,j,k)
      end do
    end do
  end do
!$OMP END DO NOWAIT
  call sft_trace("betax_rhs", ts0_thr, sft_get_ts(), 0, tid)
  ts0_thr = sft_get_ts()
#elif (GAUGE == 1)
!$OMP BARRIER
!$OMP DO COLLAPSE(1)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        betax_rhs(i,j,k) = Gamx(i,j,k) - eta*betax(i,j,k)
        betay_rhs(i,j,k) = Gamy(i,j,k) - eta*betay(i,j,k)
        betaz_rhs(i,j,k) = Gamz(i,j,k) - eta*betaz(i,j,k)

        dtSfx_rhs(i,j,k) = ZEO
        dtSfy_rhs(i,j,k) = ZEO
        dtSfz_rhs(i,j,k) = ZEO
      end do
    end do
  end do
!$OMP END DO NOWAIT
  call sft_trace("betax_rhs", ts0_thr, sft_get_ts(), 0, tid)
  ts0_thr = sft_get_ts()
#elif (GAUGE == 2)
!$OMP BARRIER
!$OMP DO COLLAPSE(1)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        betax_rhs(i,j,k) = FF*dtSfx(i,j,k)
        betay_rhs(i,j,k) = FF*dtSfy(i,j,k)
        betaz_rhs(i,j,k) = FF*dtSfz(i,j,k)
      end do
    end do
  end do
!$OMP END DO NOWAIT
  call sft_trace("betax_rhs", ts0_thr, sft_get_ts(), 0, tid)
  ts0_thr = sft_get_ts()

!$OMP BARRIER
    call symmetry_bd_inner(OPT_SYMORD,ex,chi,fh_d,SSS)
!$OMP BARRIER
    call fderivs_inner(ex,fh_d,dtSfx_rhs,dtSfy_rhs,dtSfz_rhs,X,Y,Z,Symmetry)
!$OMP BARRIER
!$OMP DO COLLAPSE(1)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        reta(i,j,k) = gupxx(i,j,k) * dtSfx_rhs(i,j,k) * dtSfx_rhs(i,j,k) + gupyy(i,j,k) * dtSfy_rhs(i,j,k) * dtSfy_rhs(i,j,k) + gupzz(i,j,k) * dtSfz_rhs(i,j,k) * dtSfz_rhs(i,j,k) + &
             TWO * (gupxy(i,j,k) * dtSfx_rhs(i,j,k) * dtSfy_rhs(i,j,k) + gupxz(i,j,k) * dtSfx_rhs(i,j,k) * dtSfz_rhs(i,j,k) + gupyz(i,j,k) * dtSfy_rhs(i,j,k) * dtSfz_rhs(i,j,k))
        reta(i,j,k) = 1.31d0/2*dsqrt(reta(i,j,k)/chin1(i,j,k))/(1-dsqrt(chin1(i,j,k)))**2
        dtSfx_rhs(i,j,k) = Gamx_rhs(i,j,k) - reta(i,j,k)*dtSfx(i,j,k)
        dtSfy_rhs(i,j,k) = Gamy_rhs(i,j,k) - reta(i,j,k)*dtSfy(i,j,k)
        dtSfz_rhs(i,j,k) = Gamz_rhs(i,j,k) - reta(i,j,k)*dtSfz(i,j,k)
      end do
    end do
  end do
!$OMP END DO NOWAIT
  call sft_trace("reta", ts0_thr, sft_get_ts(), 0, tid)
  ts0_thr = sft_get_ts()
#elif (GAUGE == 3)
!$OMP BARRIER
!$OMP DO COLLAPSE(1)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        betax_rhs(i,j,k) = FF*dtSfx(i,j,k)
        betay_rhs(i,j,k) = FF*dtSfy(i,j,k)
        betaz_rhs(i,j,k) = FF*dtSfz(i,j,k)
      end do
    end do
  end do
!$OMP END DO NOWAIT
  call sft_trace("betax_rhs", ts0_thr, sft_get_ts(), 0, tid)
  ts0_thr = sft_get_ts()

!$OMP BARRIER
    call symmetry_bd_inner(OPT_SYMORD,ex,chi,fh_d,SSS)
!$OMP BARRIER
    call fderivs_inner(ex,fh_d,dtSfx_rhs,dtSfy_rhs,dtSfz_rhs,X,Y,Z,Symmetry)
!$OMP BARRIER
!$OMP DO COLLAPSE(1)
  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        reta(i,j,k) = gupxx(i,j,k) * dtSfx_rhs(i,j,k) * dtSfx_rhs(i,j,k) + gupyy(i,j,k) * dtSfy_rhs(i,j,k) * dtSfy_rhs(i,j,k) + gupzz(i,j,k) * dtSfz_rhs(i,j,k) * dtSfz_rhs(i,j,k) + &
             TWO * (gupxy(i,j,k) * dtSfx_rhs(i,j,k) * dtSfy_rhs(i,j,k) + gupxz(i,j,k) * dtSfx_rhs(i,j,k) * dtSfz_rhs(i,j,k) + gupyz(i,j,k) * dtSfy_rhs(i,j,k) * dtSfz_rhs(i,j,k))
        reta(i,j,k) = 1.31d0/2*dsqrt(reta(i,j,k)/chin1(i,j,k))/(1-chin1(i,j,k))**2
        dtSfx_rhs(i,j,k) = Gamx_rhs(i,j,k) - reta(i,j,k)*dtSfx(i,j,k)
        dtSfy_rhs(i,j,k) = Gamy_rhs(i,j,k) - reta(i,j,k)*dtSfy(i,j,k)
        dtSfz_rhs(i,j,k) = Gamz_rhs(i,j,k) - reta(i,j,k)*dtSfz(i,j,k)
      end do
    end do
  end do
!$OMP END DO NOWAIT
  call sft_trace("reta", ts0_thr, sft_get_ts(), 0, tid)
  ts0_thr = sft_get_ts()
#elif (GAUGE == 4)
!$OMP BARRIER
    call symmetry_bd_inner(OPT_SYMORD,ex,chi,fh_d,SSS)
!$OMP BARRIER
    call fderivs_inner(ex,fh_d,dtSfx_rhs,dtSfy_rhs,dtSfz_rhs,X,Y,Z,Symmetry)
!$OMP SINGLE
  reta = gupxx * dtSfx_rhs * dtSfx_rhs + gupyy * dtSfy_rhs * dtSfy_rhs + gupzz * dtSfz_rhs * dtSfz_rhs + &
       TWO * (gupxy * dtSfx_rhs * dtSfy_rhs + gupxz * dtSfx_rhs * dtSfz_rhs + gupyz * dtSfy_rhs * dtSfz_rhs)
  reta = 1.31d0/2*dsqrt(reta/chin1)/(1-dsqrt(chin1))**2
  betax_rhs = FF*Gamx - reta*betax
  betay_rhs = FF*Gamy - reta*betay
  betaz_rhs = FF*Gamz - reta*betaz

  dtSfx_rhs = ZEO
  dtSfy_rhs = ZEO
  dtSfz_rhs = ZEO
!$OMP END SINGLE
#elif (GAUGE == 5)
!$OMP BARRIER
    call symmetry_bd_inner(OPT_SYMORD,ex,chi,fh_d,SSS)
!$OMP BARRIER
    call fderivs_inner(ex,fh_d,dtSfx_rhs,dtSfy_rhs,dtSfz_rhs,X,Y,Z,Symmetry)
!$OMP SINGLE
  reta = gupxx * dtSfx_rhs * dtSfx_rhs + gupyy * dtSfy_rhs * dtSfy_rhs + gupzz * dtSfz_rhs * dtSfz_rhs + &
       TWO * (gupxy * dtSfx_rhs * dtSfy_rhs + gupxz * dtSfx_rhs * dtSfz_rhs + gupyz * dtSfy_rhs * dtSfz_rhs)
  reta = 1.31d0/2*dsqrt(reta/chin1)/(1-chin1)**2
  betax_rhs = FF*Gamx - reta*betax
  betay_rhs = FF*Gamy - reta*betay
  betaz_rhs = FF*Gamz - reta*betaz

  dtSfx_rhs = ZEO
  dtSfy_rhs = ZEO
  dtSfz_rhs = ZEO
!$OMP END SINGLE
#elif (GAUGE == 6)
!$OMP SINGLE
  if(BHN==2)then
   M = Mass(1)+Mass(2)
   A = 2.d0/M
   w1 = 1.2d1
   w2 = w1
   C1 = 1.d0/Mass(1) - A
   C2 = 1.d0/Mass(2) - A

   do k=1,ex(3)
   do j=1,ex(2)
   do i=1,ex(1)
     r1 = ((Porg(1)-X(i))**2+(Porg(2)-Y(j))**2+(Porg(3)-Z(k))**2)/ &
          ((Porg(1)-Porg(4))**2+(Porg(2)-Porg(5))**2+(Porg(3)-Porg(6))**2)
     r2 = ((Porg(4)-X(i))**2+(Porg(5)-Y(j))**2+(Porg(6)-Z(k))**2)/ &
          ((Porg(1)-Porg(4))**2+(Porg(2)-Porg(5))**2+(Porg(3)-Porg(6))**2)
     reta(i,j,k) = A + C1/(ONE+w1*r1) + C2/(ONE+w2*r2)
    enddo
    enddo
    enddo
  else
    write(*,*) "not support BH_num in Jason's form 1",BHN
  endif
  betax_rhs = FF*dtSfx
  betay_rhs = FF*dtSfy
  betaz_rhs = FF*dtSfz

  dtSfx_rhs = Gamx_rhs - reta*dtSfx
  dtSfy_rhs = Gamy_rhs - reta*dtSfy
  dtSfz_rhs = Gamz_rhs - reta*dtSfz
!$OMP END SINGLE
#elif (GAUGE == 7)
!$OMP SINGLE
  if(BHN==2)then
   M = Mass(1)+Mass(2)
   A = 2.d0/M
   w1 = 1.2d1
   w2 = w1
   C1 = 1.d0/Mass(1) - A
   C2 = 1.d0/Mass(2) - A

   do k=1,ex(3)
   do j=1,ex(2)
   do i=1,ex(1)
     r1 = ((Porg(1)-X(i))**2+(Porg(2)-Y(j))**2+(Porg(3)-Z(k))**2)/ &
          ((Porg(1)-Porg(4))**2+(Porg(2)-Porg(5))**2+(Porg(3)-Porg(6))**2)
     r2 = ((Porg(4)-X(i))**2+(Porg(5)-Y(j))**2+(Porg(6)-Z(k))**2)/ &
          ((Porg(1)-Porg(4))**2+(Porg(2)-Porg(5))**2+(Porg(3)-Porg(6))**2)
     reta(i,j,k) = A + C1*dexp(-w1*r1) + C2*dexp(-w2*r2)
    enddo
    enddo
    enddo
  else
    write(*,*) "not support BH_num in Jason's form 2",BHN
  endif
  betax_rhs = FF*dtSfx
  betay_rhs = FF*dtSfy
  betaz_rhs = FF*dtSfz

  dtSfx_rhs = Gamx_rhs - reta*dtSfx
  dtSfy_rhs = Gamy_rhs - reta*dtSfy
  dtSfz_rhs = Gamz_rhs - reta*dtSfz
!$OMP END SINGLE
#endif  

!!!!!!!!!advection term part (variable-level parallelism)

!$OMP END PARALLEL

!$OMP PARALLEL SECTIONS DEFAULT(SHARED) PRIVATE(ts0_thr, tid)
!$OMP SECTION
    tid = omp_get_thread_num()
    ts0_thr = sft_get_ts()
    call lopsided(ex,X,Y,Z,gxx,gxx_rhs,betax,betay,betaz,Symmetry,SSS)
    if(eps>0) call kodis(ex,X,Y,Z,dxx,gxx_rhs,SSS,Symmetry,eps)
    call sft_trace("advec_gxx", ts0_thr, sft_get_ts(), 0, tid)
!$OMP SECTION
    tid = omp_get_thread_num()
    ts0_thr = sft_get_ts()
    call lopsided(ex,X,Y,Z,gxy,gxy_rhs,betax,betay,betaz,Symmetry,AAS)
    if(eps>0) call kodis(ex,X,Y,Z,gxy,gxy_rhs,AAS,Symmetry,eps)
    call sft_trace("advec_gxy", ts0_thr, sft_get_ts(), 0, tid)
!$OMP SECTION
    tid = omp_get_thread_num()
    ts0_thr = sft_get_ts()
    call lopsided(ex,X,Y,Z,gxz,gxz_rhs,betax,betay,betaz,Symmetry,ASA)
    if(eps>0) call kodis(ex,X,Y,Z,gxz,gxz_rhs,ASA,Symmetry,eps)
    call sft_trace("advec_gxz", ts0_thr, sft_get_ts(), 0, tid)
!$OMP SECTION
    tid = omp_get_thread_num()
    ts0_thr = sft_get_ts()
    call lopsided(ex,X,Y,Z,gyy,gyy_rhs,betax,betay,betaz,Symmetry,SSS)
    if(eps>0) call kodis(ex,X,Y,Z,dyy,gyy_rhs,SSS,Symmetry,eps)
    call sft_trace("advec_gyy", ts0_thr, sft_get_ts(), 0, tid)
!$OMP SECTION
    tid = omp_get_thread_num()
    ts0_thr = sft_get_ts()
    call lopsided(ex,X,Y,Z,gyz,gyz_rhs,betax,betay,betaz,Symmetry,SAA)
    if(eps>0) call kodis(ex,X,Y,Z,gyz,gyz_rhs,SAA,Symmetry,eps)
    call sft_trace("advec_gyz", ts0_thr, sft_get_ts(), 0, tid)
!$OMP SECTION
    tid = omp_get_thread_num()
    ts0_thr = sft_get_ts()
    call lopsided(ex,X,Y,Z,gzz,gzz_rhs,betax,betay,betaz,Symmetry,SSS)
    if(eps>0) call kodis(ex,X,Y,Z,dzz,gzz_rhs,SSS,Symmetry,eps)
    call sft_trace("advec_gzz", ts0_thr, sft_get_ts(), 0, tid)
!$OMP SECTION
    tid = omp_get_thread_num()
    ts0_thr = sft_get_ts()
    call lopsided(ex,X,Y,Z,Axx,Axx_rhs,betax,betay,betaz,Symmetry,SSS)
    if(eps>0) call kodis(ex,X,Y,Z,Axx,Axx_rhs,SSS,Symmetry,eps)
    call sft_trace("advec_Axx", ts0_thr, sft_get_ts(), 0, tid)
!$OMP SECTION
    tid = omp_get_thread_num()
    ts0_thr = sft_get_ts()
    call lopsided(ex,X,Y,Z,Axy,Axy_rhs,betax,betay,betaz,Symmetry,AAS)
    if(eps>0) call kodis(ex,X,Y,Z,Axy,Axy_rhs,AAS,Symmetry,eps)
    call sft_trace("advec_Axy", ts0_thr, sft_get_ts(), 0, tid)
!$OMP SECTION
    tid = omp_get_thread_num()
    ts0_thr = sft_get_ts()
    call lopsided(ex,X,Y,Z,Axz,Axz_rhs,betax,betay,betaz,Symmetry,ASA)
    if(eps>0) call kodis(ex,X,Y,Z,Axz,Axz_rhs,ASA,Symmetry,eps)
    call sft_trace("advec_Axz", ts0_thr, sft_get_ts(), 0, tid)
!$OMP SECTION
    tid = omp_get_thread_num()
    ts0_thr = sft_get_ts()
    call lopsided(ex,X,Y,Z,Ayy,Ayy_rhs,betax,betay,betaz,Symmetry,SSS)
    if(eps>0) call kodis(ex,X,Y,Z,Ayy,Ayy_rhs,SSS,Symmetry,eps)
    call sft_trace("advec_Ayy", ts0_thr, sft_get_ts(), 0, tid)
!$OMP SECTION
    tid = omp_get_thread_num()
    ts0_thr = sft_get_ts()
    call lopsided(ex,X,Y,Z,Ayz,Ayz_rhs,betax,betay,betaz,Symmetry,SAA)
    if(eps>0) call kodis(ex,X,Y,Z,Ayz,Ayz_rhs,SAA,Symmetry,eps)
    call sft_trace("advec_Ayz", ts0_thr, sft_get_ts(), 0, tid)
!$OMP SECTION
    tid = omp_get_thread_num()
    ts0_thr = sft_get_ts()
    call lopsided(ex,X,Y,Z,Azz,Azz_rhs,betax,betay,betaz,Symmetry,SSS)
    if(eps>0) call kodis(ex,X,Y,Z,Azz,Azz_rhs,SSS,Symmetry,eps)
    call sft_trace("advec_Azz", ts0_thr, sft_get_ts(), 0, tid)
!$OMP SECTION
    tid = omp_get_thread_num()
    ts0_thr = sft_get_ts()
    call lopsided(ex,X,Y,Z,chi,chi_rhs,betax,betay,betaz,Symmetry,SSS)
    if(eps>0) call kodis(ex,X,Y,Z,chi,chi_rhs,SSS,Symmetry,eps)
    call sft_trace("advec_chi", ts0_thr, sft_get_ts(), 0, tid)
!$OMP SECTION
    tid = omp_get_thread_num()
    ts0_thr = sft_get_ts()
    call lopsided(ex,X,Y,Z,trK,trK_rhs,betax,betay,betaz,Symmetry,SSS)
    if(eps>0) call kodis(ex,X,Y,Z,trK,trK_rhs,SSS,Symmetry,eps)
    call sft_trace("advec_trK", ts0_thr, sft_get_ts(), 0, tid)
!$OMP SECTION
    tid = omp_get_thread_num()
    ts0_thr = sft_get_ts()
    call lopsided(ex,X,Y,Z,Gamx,Gamx_rhs,betax,betay,betaz,Symmetry,ASS)
    if(eps>0) call kodis(ex,X,Y,Z,Gamx,Gamx_rhs,ASS,Symmetry,eps)
    call sft_trace("advec_Gamx", ts0_thr, sft_get_ts(), 0, tid)
!$OMP SECTION
    tid = omp_get_thread_num()
    ts0_thr = sft_get_ts()
    call lopsided(ex,X,Y,Z,Gamy,Gamy_rhs,betax,betay,betaz,Symmetry,SAS)
    if(eps>0) call kodis(ex,X,Y,Z,Gamy,Gamy_rhs,SAS,Symmetry,eps)
    call sft_trace("advec_Gamy", ts0_thr, sft_get_ts(), 0, tid)
!$OMP SECTION
    tid = omp_get_thread_num()
    ts0_thr = sft_get_ts()
    call lopsided(ex,X,Y,Z,Gamz,Gamz_rhs,betax,betay,betaz,Symmetry,SSA)
    if(eps>0) call kodis(ex,X,Y,Z,Gamz,Gamz_rhs,SSA,Symmetry,eps)
    call sft_trace("advec_Gamz", ts0_thr, sft_get_ts(), 0, tid)
!$OMP SECTION
    tid = omp_get_thread_num()
    ts0_thr = sft_get_ts()
    call lopsided(ex,X,Y,Z,Lap,Lap_rhs,betax,betay,betaz,Symmetry,SSS)
    if(eps>0) call kodis(ex,X,Y,Z,Lap,Lap_rhs,SSS,Symmetry,eps)
    call sft_trace("advec_Lap", ts0_thr, sft_get_ts(), 0, tid)
#if (GAUGE == 0 || GAUGE == 1 || GAUGE == 2 || GAUGE == 3 || GAUGE == 4 || GAUGE == 5 || GAUGE == 6 || GAUGE == 7)
!$OMP SECTION
    tid = omp_get_thread_num()
    ts0_thr = sft_get_ts()
    call lopsided(ex,X,Y,Z,betax,betax_rhs,betax,betay,betaz,Symmetry,ASS)
    if(eps>0) call kodis(ex,X,Y,Z,betax,betax_rhs,ASS,Symmetry,eps)
    call sft_trace("advec_betax", ts0_thr, sft_get_ts(), 0, tid)
!$OMP SECTION
    tid = omp_get_thread_num()
    ts0_thr = sft_get_ts()
    call lopsided(ex,X,Y,Z,betay,betay_rhs,betax,betay,betaz,Symmetry,SAS)
    if(eps>0) call kodis(ex,X,Y,Z,betay,betay_rhs,SAS,Symmetry,eps)
    call sft_trace("advec_betay", ts0_thr, sft_get_ts(), 0, tid)
!$OMP SECTION
    tid = omp_get_thread_num()
    ts0_thr = sft_get_ts()
    call lopsided(ex,X,Y,Z,betaz,betaz_rhs,betax,betay,betaz,Symmetry,SSA)
    if(eps>0) call kodis(ex,X,Y,Z,betaz,betaz_rhs,SSA,Symmetry,eps)
    call sft_trace("advec_betaz", ts0_thr, sft_get_ts(), 0, tid)
#endif
#if (GAUGE == 0 || GAUGE == 2 || GAUGE == 3 || GAUGE == 6 || GAUGE == 7)
!$OMP SECTION
    tid = omp_get_thread_num()
    ts0_thr = sft_get_ts()
    call lopsided(ex,X,Y,Z,dtSfx,dtSfx_rhs,betax,betay,betaz,Symmetry,ASS)
    if(eps>0) call kodis(ex,X,Y,Z,dtSfx,dtSfx_rhs,ASS,Symmetry,eps)
    call sft_trace("advec_dtSfx", ts0_thr, sft_get_ts(), 0, tid)
!$OMP SECTION
    tid = omp_get_thread_num()
    ts0_thr = sft_get_ts()
    call lopsided(ex,X,Y,Z,dtSfy,dtSfy_rhs,betax,betay,betaz,Symmetry,SAS)
    if(eps>0) call kodis(ex,X,Y,Z,dtSfy,dtSfy_rhs,SAS,Symmetry,eps)
    call sft_trace("advec_dtSfy", ts0_thr, sft_get_ts(), 0, tid)
!$OMP SECTION
    tid = omp_get_thread_num()
    ts0_thr = sft_get_ts()
    call lopsided(ex,X,Y,Z,dtSfz,dtSfz_rhs,betax,betay,betaz,Symmetry,SSA)
    if(eps>0) call kodis(ex,X,Y,Z,dtSfz,dtSfz_rhs,SSA,Symmetry,eps)
    call sft_trace("advec_dtSfz", ts0_thr, sft_get_ts(), 0, tid)
#endif
!$OMP END PARALLEL SECTIONS

  if(co == 0)then
! ham_Res = trR + 2/3 * K^2 - A_ij * A^ij - 16 * PI * rho
! here trR is respect to physical metric
!$OMP PARALLEL PRIVATE(i,j,k,ts0_thr,tid) DEFAULT(SHARED)
  tid = omp_get_thread_num()
  ts0_thr = sft_get_ts()
!$OMP DO COLLAPSE(1)
  do k = 1, ex(3)
    do j = 1, ex(2)

      do i = 1, ex(1)
        ham_Res(i,j,k) = gupxx(i,j,k) * Rxx(i,j,k) + gupyy(i,j,k) * Ryy(i,j,k) + gupzz(i,j,k) * Rzz(i,j,k) + &
              TWO * ( gupxy(i,j,k) * Rxy(i,j,k) + gupxz(i,j,k) * Rxz(i,j,k) + gupyz(i,j,k) * Ryz(i,j,k) )

        ham_Res(i,j,k) = chin1(i,j,k) * ham_Res(i,j,k) + F2o3 * trK(i,j,k) * trK(i,j,k) -( &
             gupxx(i,j,k) * ( &
             gupxx(i,j,k) * Axx(i,j,k) * Axx(i,j,k) + gupyy(i,j,k) * Axy(i,j,k) * Axy(i,j,k) + gupzz(i,j,k) * Axz(i,j,k) * Axz(i,j,k) + &
             TWO * (gupxy(i,j,k) * Axx(i,j,k) * Axy(i,j,k) + gupxz(i,j,k) * Axx(i,j,k) * Axz(i,j,k) + gupyz(i,j,k) * Axy(i,j,k) * Axz(i,j,k)) ) + &
             gupyy(i,j,k) * ( &
             gupxx(i,j,k) * Axy(i,j,k) * Axy(i,j,k) + gupyy(i,j,k) * Ayy(i,j,k) * Ayy(i,j,k) + gupzz(i,j,k) * Ayz(i,j,k) * Ayz(i,j,k) + &
             TWO * (gupxy(i,j,k) * Axy(i,j,k) * Ayy(i,j,k) + gupxz(i,j,k) * Axy(i,j,k) * Ayz(i,j,k) + gupyz(i,j,k) * Ayy(i,j,k) * Ayz(i,j,k)) ) + &
             gupzz(i,j,k) * ( &
             gupxx(i,j,k) * Axz(i,j,k) * Axz(i,j,k) + gupyy(i,j,k) * Ayz(i,j,k) * Ayz(i,j,k) + gupzz(i,j,k) * Azz(i,j,k) * Azz(i,j,k) + &
             TWO * (gupxy(i,j,k) * Axz(i,j,k) * Ayz(i,j,k) + gupxz(i,j,k) * Axz(i,j,k) * Azz(i,j,k) + gupyz(i,j,k) * Ayz(i,j,k) * Azz(i,j,k)) ) + &
             TWO * ( &
             gupxy(i,j,k) * ( &
             gupxx(i,j,k) * Axx(i,j,k) * Axy(i,j,k) + gupyy(i,j,k) * Axy(i,j,k) * Ayy(i,j,k) + gupzz(i,j,k) * Axz(i,j,k) * Ayz(i,j,k) + &
             gupxy(i,j,k) * (Axx(i,j,k) * Ayy(i,j,k) + Axy(i,j,k) * Axy(i,j,k)) + &
             gupxz(i,j,k) * (Axx(i,j,k) * Ayz(i,j,k) + Axz(i,j,k) * Axy(i,j,k)) + &
             gupyz(i,j,k) * (Axy(i,j,k) * Ayz(i,j,k) + Axz(i,j,k) * Ayy(i,j,k)) ) + &
             gupxz(i,j,k) * ( &
             gupxx(i,j,k) * Axx(i,j,k) * Axz(i,j,k) + gupyy(i,j,k) * Axy(i,j,k) * Ayz(i,j,k) + gupzz(i,j,k) * Axz(i,j,k) * Azz(i,j,k) + &
             gupxy(i,j,k) * (Axx(i,j,k) * Ayz(i,j,k) + Axy(i,j,k) * Axz(i,j,k)) + &
             gupxz(i,j,k) * (Axx(i,j,k) * Azz(i,j,k) + Axz(i,j,k) * Axz(i,j,k)) + &
             gupyz(i,j,k) * (Axy(i,j,k) * Azz(i,j,k) + Axz(i,j,k) * Ayz(i,j,k)) ) + &
             gupyz(i,j,k) * ( &
             gupxx(i,j,k) * Axy(i,j,k) * Axz(i,j,k) + gupyy(i,j,k) * Ayy(i,j,k) * Ayz(i,j,k) + gupzz(i,j,k) * Ayz(i,j,k) * Azz(i,j,k) + &
             gupxy(i,j,k) * (Axy(i,j,k) * Ayz(i,j,k) + Ayy(i,j,k) * Axz(i,j,k)) + &
             gupxz(i,j,k) * (Axy(i,j,k) * Azz(i,j,k) + Ayz(i,j,k) * Axz(i,j,k)) + &
             gupyz(i,j,k) * (Ayy(i,j,k) * Azz(i,j,k) + Ayz(i,j,k) * Ayz(i,j,k)) ) )) - F16 * PI * rho(i,j,k)
      end do

    end do
  end do
!$OMP END DO NOWAIT
  call sft_trace("ham_Res", ts0_thr, sft_get_ts(), 0, tid)
  ts0_thr = sft_get_ts()

! mov_Res_j = gupkj*(-1/chi d_k chi*A_ij + D_k A_ij) - 2/3 d_j trK - 8 PI s_j where D respect to physical metric
! store D_i A_jk - 1/chi d_i chi*A_jk in gjk_i
! Batched Aij constraint derivatives — 6 sym_bd + 6 fderivs
    call symmetry_bd_inner(OPT_SYMORD,ex,Axx,fh_d,SSS)
    call symmetry_bd_inner(OPT_SYMORD,ex,Axy,fh_d2,AAS)
    call symmetry_bd_inner(OPT_SYMORD,ex,Axz,fh_d3,ASA)
    call symmetry_bd_inner(OPT_SYMORD,ex,Ayy,fh_d4,SSS)
    call symmetry_bd_inner(OPT_SYMORD,ex,Ayz,fh_d5,SAA)
    call symmetry_bd_inner(OPT_SYMORD,ex,Azz,fh_d6,SSS)
!$OMP BARRIER
    call fderivs_inner(ex,fh_d,gxxx,gxxy,gxxz,X,Y,Z,Symmetry)
    call fderivs_inner(ex,fh_d2,gxyx,gxyy,gxyz,X,Y,Z,Symmetry)
    call fderivs_inner(ex,fh_d3,gxzx,gxzy,gxzz,X,Y,Z,Symmetry)
    call fderivs_inner(ex,fh_d4,gyyx,gyyy,gyyz,X,Y,Z,Symmetry)
    call fderivs_inner(ex,fh_d5,gyzx,gyzy,gyzz,X,Y,Z,Symmetry)
    call fderivs_inner(ex,fh_d6,gzzx,gzzy,gzzz,X,Y,Z,Symmetry)
!$OMP BARRIER
!$OMP DO COLLAPSE(1)
  do k = 1, ex(3)
    do j = 1, ex(2)

      do i = 1, ex(1)
        gxxx(i,j,k) = gxxx(i,j,k) - (  Gamxxx(i,j,k) * Axx(i,j,k) + Gamyxx(i,j,k) * Axy(i,j,k) + Gamzxx(i,j,k) * Axz(i,j,k) &
                     + Gamxxx(i,j,k) * Axx(i,j,k) + Gamyxx(i,j,k) * Axy(i,j,k) + Gamzxx(i,j,k) * Axz(i,j,k)) - chix(i,j,k)*Axx(i,j,k)/chin1(i,j,k)
        gxyx(i,j,k) = gxyx(i,j,k) - (  Gamxxy(i,j,k) * Axx(i,j,k) + Gamyxy(i,j,k) * Axy(i,j,k) + Gamzxy(i,j,k) * Axz(i,j,k) &
                     + Gamxxx(i,j,k) * Axy(i,j,k) + Gamyxx(i,j,k) * Ayy(i,j,k) + Gamzxx(i,j,k) * Ayz(i,j,k)) - chix(i,j,k)*Axy(i,j,k)/chin1(i,j,k)
        gxzx(i,j,k) = gxzx(i,j,k) - (  Gamxxz(i,j,k) * Axx(i,j,k) + Gamyxz(i,j,k) * Axy(i,j,k) + Gamzxz(i,j,k) * Axz(i,j,k) &
                     + Gamxxx(i,j,k) * Axz(i,j,k) + Gamyxx(i,j,k) * Ayz(i,j,k) + Gamzxx(i,j,k) * Azz(i,j,k)) - chix(i,j,k)*Axz(i,j,k)/chin1(i,j,k)
        gyyx(i,j,k) = gyyx(i,j,k) - (  Gamxxy(i,j,k) * Axy(i,j,k) + Gamyxy(i,j,k) * Ayy(i,j,k) + Gamzxy(i,j,k) * Ayz(i,j,k) &
                     + Gamxxy(i,j,k) * Axy(i,j,k) + Gamyxy(i,j,k) * Ayy(i,j,k) + Gamzxy(i,j,k) * Ayz(i,j,k)) - chix(i,j,k)*Ayy(i,j,k)/chin1(i,j,k)
        gyzx(i,j,k) = gyzx(i,j,k) - (  Gamxxz(i,j,k) * Axy(i,j,k) + Gamyxz(i,j,k) * Ayy(i,j,k) + Gamzxz(i,j,k) * Ayz(i,j,k) &
                     + Gamxxy(i,j,k) * Axz(i,j,k) + Gamyxy(i,j,k) * Ayz(i,j,k) + Gamzxy(i,j,k) * Azz(i,j,k)) - chix(i,j,k)*Ayz(i,j,k)/chin1(i,j,k)
        gzzx(i,j,k) = gzzx(i,j,k) - (  Gamxxz(i,j,k) * Axz(i,j,k) + Gamyxz(i,j,k) * Ayz(i,j,k) + Gamzxz(i,j,k) * Azz(i,j,k) &
                     + Gamxxz(i,j,k) * Axz(i,j,k) + Gamyxz(i,j,k) * Ayz(i,j,k) + Gamzxz(i,j,k) * Azz(i,j,k)) - chix(i,j,k)*Azz(i,j,k)/chin1(i,j,k)
        gxxy(i,j,k) = gxxy(i,j,k) - (  Gamxxy(i,j,k) * Axx(i,j,k) + Gamyxy(i,j,k) * Axy(i,j,k) + Gamzxy(i,j,k) * Axz(i,j,k) &
                     + Gamxxy(i,j,k) * Axx(i,j,k) + Gamyxy(i,j,k) * Axy(i,j,k) + Gamzxy(i,j,k) * Axz(i,j,k)) - chiy(i,j,k)*Axx(i,j,k)/chin1(i,j,k)
        gxyy(i,j,k) = gxyy(i,j,k) - (  Gamxyy(i,j,k) * Axx(i,j,k) + Gamyyy(i,j,k) * Axy(i,j,k) + Gamzyy(i,j,k) * Axz(i,j,k) &
                     + Gamxxy(i,j,k) * Axy(i,j,k) + Gamyxy(i,j,k) * Ayy(i,j,k) + Gamzxy(i,j,k) * Ayz(i,j,k)) - chiy(i,j,k)*Axy(i,j,k)/chin1(i,j,k)
        gxzy(i,j,k) = gxzy(i,j,k) - (  Gamxyz(i,j,k) * Axx(i,j,k) + Gamyyz(i,j,k) * Axy(i,j,k) + Gamzyz(i,j,k) * Axz(i,j,k) &
                     + Gamxxy(i,j,k) * Axz(i,j,k) + Gamyxy(i,j,k) * Ayz(i,j,k) + Gamzxy(i,j,k) * Azz(i,j,k)) - chiy(i,j,k)*Axz(i,j,k)/chin1(i,j,k)
        gyyy(i,j,k) = gyyy(i,j,k) - (  Gamxyy(i,j,k) * Axy(i,j,k) + Gamyyy(i,j,k) * Ayy(i,j,k) + Gamzyy(i,j,k) * Ayz(i,j,k) &
                     + Gamxyy(i,j,k) * Axy(i,j,k) + Gamyyy(i,j,k) * Ayy(i,j,k) + Gamzyy(i,j,k) * Ayz(i,j,k)) - chiy(i,j,k)*Ayy(i,j,k)/chin1(i,j,k)
        gyzy(i,j,k) = gyzy(i,j,k) - (  Gamxyz(i,j,k) * Axy(i,j,k) + Gamyyz(i,j,k) * Ayy(i,j,k) + Gamzyz(i,j,k) * Ayz(i,j,k) &
                     + Gamxyy(i,j,k) * Axz(i,j,k) + Gamyyy(i,j,k) * Ayz(i,j,k) + Gamzyy(i,j,k) * Azz(i,j,k)) - chiy(i,j,k)*Ayz(i,j,k)/chin1(i,j,k)
        gzzy(i,j,k) = gzzy(i,j,k) - (  Gamxyz(i,j,k) * Axz(i,j,k) + Gamyyz(i,j,k) * Ayz(i,j,k) + Gamzyz(i,j,k) * Azz(i,j,k) &
                     + Gamxyz(i,j,k) * Axz(i,j,k) + Gamyyz(i,j,k) * Ayz(i,j,k) + Gamzyz(i,j,k) * Azz(i,j,k)) - chiy(i,j,k)*Azz(i,j,k)/chin1(i,j,k)
        gxxz(i,j,k) = gxxz(i,j,k) - (  Gamxxz(i,j,k) * Axx(i,j,k) + Gamyxz(i,j,k) * Axy(i,j,k) + Gamzxz(i,j,k) * Axz(i,j,k) &
                     + Gamxxz(i,j,k) * Axx(i,j,k) + Gamyxz(i,j,k) * Axy(i,j,k) + Gamzxz(i,j,k) * Axz(i,j,k)) - chiz(i,j,k)*Axx(i,j,k)/chin1(i,j,k)
        gxyz(i,j,k) = gxyz(i,j,k) - (  Gamxyz(i,j,k) * Axx(i,j,k) + Gamyyz(i,j,k) * Axy(i,j,k) + Gamzyz(i,j,k) * Axz(i,j,k) &
                     + Gamxxz(i,j,k) * Axy(i,j,k) + Gamyxz(i,j,k) * Ayy(i,j,k) + Gamzxz(i,j,k) * Ayz(i,j,k)) - chiz(i,j,k)*Axy(i,j,k)/chin1(i,j,k)
        gxzz(i,j,k) = gxzz(i,j,k) - (  Gamxzz(i,j,k) * Axx(i,j,k) + Gamyzz(i,j,k) * Axy(i,j,k) + Gamzzz(i,j,k) * Axz(i,j,k) &
                     + Gamxxz(i,j,k) * Axz(i,j,k) + Gamyxz(i,j,k) * Ayz(i,j,k) + Gamzxz(i,j,k) * Azz(i,j,k)) - chiz(i,j,k)*Axz(i,j,k)/chin1(i,j,k)
        gyyz(i,j,k) = gyyz(i,j,k) - (  Gamxyz(i,j,k) * Axy(i,j,k) + Gamyyz(i,j,k) * Ayy(i,j,k) + Gamzyz(i,j,k) * Ayz(i,j,k) &
                     + Gamxyz(i,j,k) * Axy(i,j,k) + Gamyyz(i,j,k) * Ayy(i,j,k) + Gamzyz(i,j,k) * Ayz(i,j,k)) - chiz(i,j,k)*Ayy(i,j,k)/chin1(i,j,k)
        gyzz(i,j,k) = gyzz(i,j,k) - (  Gamxzz(i,j,k) * Axy(i,j,k) + Gamyzz(i,j,k) * Ayy(i,j,k) + Gamzzz(i,j,k) * Ayz(i,j,k) &
                     + Gamxyz(i,j,k) * Axz(i,j,k) + Gamyyz(i,j,k) * Ayz(i,j,k) + Gamzyz(i,j,k) * Azz(i,j,k)) - chiz(i,j,k)*Ayz(i,j,k)/chin1(i,j,k)
        gzzz(i,j,k) = gzzz(i,j,k) - (  Gamxzz(i,j,k) * Axz(i,j,k) + Gamyzz(i,j,k) * Ayz(i,j,k) + Gamzzz(i,j,k) * Azz(i,j,k) &
                     + Gamxzz(i,j,k) * Axz(i,j,k) + Gamyzz(i,j,k) * Ayz(i,j,k) + Gamzzz(i,j,k) * Azz(i,j,k)) - chiz(i,j,k)*Azz(i,j,k)/chin1(i,j,k)
        movx_Res(i,j,k) = gupxx(i,j,k)*gxxx(i,j,k) + gupyy(i,j,k)*gxyy(i,j,k) + gupzz(i,j,k)*gxzz(i,j,k) &
                          +gupxy(i,j,k)*gxyx(i,j,k) + gupxz(i,j,k)*gxzx(i,j,k) + gupyz(i,j,k)*gxzy(i,j,k) &
                          +gupxy(i,j,k)*gxxy(i,j,k) + gupxz(i,j,k)*gxxz(i,j,k) + gupyz(i,j,k)*gxyz(i,j,k)
        movy_Res(i,j,k) = gupxx(i,j,k)*gxyx(i,j,k) + gupyy(i,j,k)*gyyy(i,j,k) + gupzz(i,j,k)*gyzz(i,j,k) &
                          +gupxy(i,j,k)*gyyx(i,j,k) + gupxz(i,j,k)*gyzx(i,j,k) + gupyz(i,j,k)*gyzy(i,j,k) &
                          +gupxy(i,j,k)*gxyy(i,j,k) + gupxz(i,j,k)*gxyz(i,j,k) + gupyz(i,j,k)*gyyz(i,j,k)
        movz_Res(i,j,k) = gupxx(i,j,k)*gxzx(i,j,k) + gupyy(i,j,k)*gyzy(i,j,k) + gupzz(i,j,k)*gzzz(i,j,k) &
                          +gupxy(i,j,k)*gyzx(i,j,k) + gupxz(i,j,k)*gzzx(i,j,k) + gupyz(i,j,k)*gzzy(i,j,k) &
                          +gupxy(i,j,k)*gxzy(i,j,k) + gupxz(i,j,k)*gxzz(i,j,k) + gupyz(i,j,k)*gyzz(i,j,k)
        movx_Res(i,j,k) = movx_Res(i,j,k) - F2o3*Kx(i,j,k) - F8*PI*sx(i,j,k)
        movy_Res(i,j,k) = movy_Res(i,j,k) - F2o3*Ky(i,j,k) - F8*PI*sy(i,j,k)
        movz_Res(i,j,k) = movz_Res(i,j,k) - F2o3*Kz(i,j,k) - F8*PI*sz(i,j,k)
      end do

    end do
  end do
!$OMP END DO NOWAIT
  call sft_trace("mov_Res", ts0_thr, sft_get_ts(), 0, tid)
  ts0_thr = sft_get_ts()
!$OMP END PARALLEL
  endif

#if (ABV == 1)
  call ricci_gamma(ex, X, Y, Z,                                      &
               chi,                                                  &
               dxx    ,   gxy    ,   gxz    ,   dyy    ,   gyz    ,   dzz,&
               Gamx   ,  Gamy    ,  Gamz    , &
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,&
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,&
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,&
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,&
               Symmetry)
  call constraint_bssn(ex, X, Y, Z,&
               chi,trK, &
               dxx,gxy,gxz,dyy,gyz,dzz, &
               Axx,Axy,Axz,Ayy,Ayz,Azz, &
               Gamx,Gamy,Gamz,&
               Lap,betax,betay,betaz,rho,Sx,Sy,Sz,&
               Gamxxx, Gamxxy, Gamxxz,Gamxyy, Gamxyz, Gamxzz, &
               Gamyxx, Gamyxy, Gamyxz,Gamyyy, Gamyyz, Gamyzz, &
               Gamzxx, Gamzxy, Gamzxz,Gamzyy, Gamzyz, Gamzzz, &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz, &
               ham_Res,movx_Res,movy_Res,movz_Res,Gmx_Res,Gmy_Res,Gmz_Res, &
               Symmetry)
#endif 
#if 0
#define i 2
if(Lev == 1)then
write(*,*) X(i),Y(i),Z(i)
write(*,*) Axx(i,i,i),Axy(i,i,i),Axz(i,i,i),Ayy(i,i,i),Ayz(i,i,i),Azz(i,i,i)
write(*,*) 1+Lap(i,i,i),dtSfx(i,i,i),dtSfy(i,i,i),dtSfz(i,i,i)
write(*,*) betax(i,i,i),betay(i,i,i),betaz(i,i,i)
write(*,*) 1+chi(i,i,i),Gamx(i,i,i),Gamy(i,i,i),Gamz(i,i,i)
write(*,*) gxx(i,i,i),gxy(i,i,i),gxz(i,i,i),gyy(i,i,i),gyz(i,i,i),gzz(i,i,i)
write(*,*) trK(i,i,i)
write(*,*) "====="
write(*,*) Axx_rhs(i,i,i),Axy_rhs(i,i,i),Axz_rhs(i,i,i),Ayy_rhs(i,i,i),Ayz_rhs(i,i,i),Azz_rhs(i,i,i)
write(*,*) Lap_rhs(i,i,i),dtSfx_rhs(i,i,i),dtSfy_rhs(i,i,i),dtSfz_rhs(i,i,i)
write(*,*) betax_rhs(i,i,i),betay_rhs(i,i,i),betaz_rhs(i,i,i)
write(*,*) chi_rhs(i,i,i),Gamx_rhs(i,i,i),Gamy_rhs(i,i,i),Gamz_rhs(i,i,i)
write(*,*) gxx_rhs(i,i,i),gxy_rhs(i,i,i),gxz_rhs(i,i,i),gyy_rhs(i,i,i),gyz_rhs(i,i,i),gzz_rhs(i,i,i)
write(*,*) trK_rhs(i,i,i)
endif
#undef i
!!stop
#endif

  gont = 0

  return

  end function compute_rhs_bssn_opt
#endif


  ! ----------------------------------------------------------------------
  ! Legacy implementation: ghost_width-aware (works for 4th, 6th, 8th order).
  ! Used by the wrapper above when ghost_width != 3 (e.g. Z4C constraint
  ! computation with 6th-order finite differences).
  ! ----------------------------------------------------------------------
  function compute_rhs_bssn_legacy(ex, T,X, Y, Z,                                     &
               chi    ,   trK    ,                                             &
               dxx    ,   gxy    ,   gxz    ,   dyy    ,   gyz    ,   dzz,     &
               Axx    ,   Axy    ,   Axz    ,   Ayy    ,   Ayz    ,   Azz,     &
               Gamx   ,  Gamy    ,  Gamz    ,                                  &
               Lap    ,  betax   ,  betay   ,  betaz   ,                       &
               dtSfx  ,  dtSfy   ,  dtSfz   ,                                  &
               chi_rhs,   trK_rhs,                                             &
               gxx_rhs,   gxy_rhs,   gxz_rhs,   gyy_rhs,   gyz_rhs,   gzz_rhs, &
               Axx_rhs,   Axy_rhs,   Axz_rhs,   Ayy_rhs,   Ayz_rhs,   Azz_rhs, &
               Gamx_rhs,  Gamy_rhs,  Gamz_rhs,                                 &
               Lap_rhs,  betax_rhs,  betay_rhs,  betaz_rhs,                    &
               dtSfx_rhs,  dtSfy_rhs,  dtSfz_rhs,                              &
               rho,Sx,Sy,Sz,Sxx,Sxy,Sxz,Syy,Syz,Szz,                           &
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,                      &
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,                      &
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,                      &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,                                        &
               ham_Res, movx_Res, movy_Res, movz_Res,                          &
                        Gmx_Res, Gmy_Res, Gmz_Res,                             &
               Symmetry,Lev,eps,co)  result(gont)
! calculate constraint violation when co=0               
  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3), Symmetry,Lev,co
  real*8, intent(in ):: T
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: chi,dxx,dyy,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: trK
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: gxy,gxz,gyz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamx,Gamy,Gamz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Lap, betax, betay, betaz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: dtSfx,  dtSfy,  dtSfz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: chi_rhs,trK_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gxx_rhs,gxy_rhs,gxz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gyy_rhs,gyz_rhs,gzz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Axx_rhs,Axy_rhs,Axz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Ayy_rhs,Ayz_rhs,Azz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamx_rhs,Gamy_rhs,Gamz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Lap_rhs, betax_rhs, betay_rhs, betaz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: dtSfx_rhs,dtSfy_rhs,dtSfz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: rho,Sx,Sy,Sz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Sxx,Sxy,Sxz,Syy,Syz,Szz
! when out, physical second kind of connection  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamzyy, Gamzyz, Gamzzz
! when out, physical Ricci tensor  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Rxx,Rxy,Rxz,Ryy,Ryz,Rzz
  real*8,intent(in) :: eps
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: ham_Res, movx_Res, movy_Res, movz_Res
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Gmx_Res, Gmy_Res, Gmz_Res
!  gont = 0: success; gont = 1: something wrong
  integer::gont

!~~~~~~> Other variables:

  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gyy,gzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: chix,chiy,chiz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxx,gxyx,gxzx,gyyx,gyzx,gzzx
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxy,gxyy,gxzy,gyyy,gyzy,gzzy
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxz,gxyz,gxzz,gyyz,gyzz,gzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Lapx,Lapy,Lapz
  real*8, dimension(ex(1),ex(2),ex(3)) :: betaxx,betaxy,betaxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: betayx,betayy,betayz
  real*8, dimension(ex(1),ex(2),ex(3)) :: betazx,betazy,betazz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamxx,Gamxy,Gamxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamyx,Gamyy,Gamyz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamzx,Gamzy,Gamzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Kx,Ky,Kz,div_beta,S
  real*8, dimension(ex(1),ex(2),ex(3)) :: f,fxx,fxy,fxz,fyy,fyz,fzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamxa,Gamya,Gamza,alpn1,chin1
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupyy,gupyz,gupzz

  real*8,dimension(3) ::SSS,AAS,ASA,SAA,ASS,SAS,SSA
  real*8            :: dX, dY, dZ, PI
  real*8, parameter :: ZEO = 0.d0,ONE = 1.D0, TWO = 2.D0, FOUR = 4.D0
  real*8, parameter :: EIGHT = 8.D0, HALF = 0.5D0, THR = 3.d0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  double precision,parameter::FF = 0.75d0,eta=2.d0
  real*8, parameter :: F1o3 = 1.D0/3.D0, F2o3 = 2.D0/3.D0,F3o2=1.5d0, F1o6 = 1.D0/6.D0
  real*8, parameter :: F16=1.6d1,F8=8.d0

#if (GAUGE == 2 || GAUGE == 3 || GAUGE == 4 || GAUGE == 5)
  real*8, dimension(ex(1),ex(2),ex(3)) :: reta
#endif

#if (GAUGE == 6 || GAUGE == 7)
  integer :: BHN,i,j,k
  real*8, dimension(9) :: Porg
  real*8, dimension(3) :: Mass
  real*8 :: r1,r2,M,A,w1,w2,C1,C2
  real*8, dimension(ex(1),ex(2),ex(3)) :: reta

  call getpbh(BHN,Porg,Mass)
#endif

!!! sanity check
  dX = sum(chi)+sum(trK)+sum(dxx)+sum(gxy)+sum(gxz)+sum(dyy)+sum(gyz)+sum(dzz) &
      +sum(Axx)+sum(Axy)+sum(Axz)+sum(Ayy)+sum(Ayz)+sum(Azz)                   &
      +sum(Gamx)+sum(Gamy)+sum(Gamz)                                           &
      +sum(Lap)+sum(betax)+sum(betay)+sum(betaz)
  if(dX.ne.dX) then
     if(sum(chi).ne.sum(chi))write(*,*)"bssn.f90: find NaN in chi"
     if(sum(trK).ne.sum(trK))write(*,*)"bssn.f90: find NaN in trk"
     if(sum(dxx).ne.sum(dxx))write(*,*)"bssn.f90: find NaN in dxx"
     if(sum(gxy).ne.sum(gxy))write(*,*)"bssn.f90: find NaN in gxy"
     if(sum(gxz).ne.sum(gxz))write(*,*)"bssn.f90: find NaN in gxz"
     if(sum(dyy).ne.sum(dyy))write(*,*)"bssn.f90: find NaN in dyy"
     if(sum(gyz).ne.sum(gyz))write(*,*)"bssn.f90: find NaN in gyz"
     if(sum(dzz).ne.sum(dzz))write(*,*)"bssn.f90: find NaN in dzz"
     if(sum(Axx).ne.sum(Axx))write(*,*)"bssn.f90: find NaN in Axx"
     if(sum(Axy).ne.sum(Axy))write(*,*)"bssn.f90: find NaN in Axy"
     if(sum(Axz).ne.sum(Axz))write(*,*)"bssn.f90: find NaN in Axz"
     if(sum(Ayy).ne.sum(Ayy))write(*,*)"bssn.f90: find NaN in Ayy"
     if(sum(Ayz).ne.sum(Ayz))write(*,*)"bssn.f90: find NaN in Ayz"
     if(sum(Azz).ne.sum(Azz))write(*,*)"bssn.f90: find NaN in Azz"
     if(sum(Gamx).ne.sum(Gamx))write(*,*)"bssn.f90: find NaN in Gamx"
     if(sum(Gamy).ne.sum(Gamy))write(*,*)"bssn.f90: find NaN in Gamy"
     if(sum(Gamz).ne.sum(Gamz))write(*,*)"bssn.f90: find NaN in Gamz"
     if(sum(Lap).ne.sum(Lap))write(*,*)"bssn.f90: find NaN in Lap"
     if(sum(betax).ne.sum(betax))write(*,*)"bssn.f90: find NaN in betax"
     if(sum(betay).ne.sum(betay))write(*,*)"bssn.f90: find NaN in betay"
     if(sum(betaz).ne.sum(betaz))write(*,*)"bssn.f90: find NaN in betaz"
     gont = 1
     return
  endif

  PI = dacos(-ONE)

  dX = X(2) - X(1)
  dY = Y(2) - Y(1)
  dZ = Z(2) - Z(1)

  alpn1 = Lap + ONE
  chin1 = chi + ONE
  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE

  call fderivs(ex,betax,betaxx,betaxy,betaxz,X,Y,Z,ANTI, SYM, SYM,Symmetry,Lev)
  call fderivs(ex,betay,betayx,betayy,betayz,X,Y,Z, SYM,ANTI, SYM,Symmetry,Lev)
  call fderivs(ex,betaz,betazx,betazy,betazz,X,Y,Z, SYM, SYM,ANTI,Symmetry,Lev)
  
  div_beta = betaxx + betayy + betazz
 
  call fderivs(ex,chi,chix,chiy,chiz,X,Y,Z,SYM,SYM,SYM,symmetry,Lev)

  chi_rhs = F2o3 *chin1*( alpn1 * trK - div_beta ) !rhs for chi

  call fderivs(ex,dxx,gxxx,gxxy,gxxz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)
  call fderivs(ex,gxy,gxyx,gxyy,gxyz,X,Y,Z,ANTI,ANTI,SYM ,Symmetry,Lev)
  call fderivs(ex,gxz,gxzx,gxzy,gxzz,X,Y,Z,ANTI,SYM ,ANTI,Symmetry,Lev)
  call fderivs(ex,dyy,gyyx,gyyy,gyyz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)
  call fderivs(ex,gyz,gyzx,gyzy,gyzz,X,Y,Z,SYM ,ANTI,ANTI,Symmetry,Lev)
  call fderivs(ex,dzz,gzzx,gzzy,gzzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)

  gxx_rhs = - TWO * alpn1 * Axx    -  F2o3 * gxx * div_beta          + &
              TWO *(  gxx * betaxx +   gxy * betayx +   gxz * betazx)

  gyy_rhs = - TWO * alpn1 * Ayy    -  F2o3 * gyy * div_beta          + &
              TWO *(  gxy * betaxy +   gyy * betayy +   gyz * betazy)

  gzz_rhs = - TWO * alpn1 * Azz    -  F2o3 * gzz * div_beta          + &
              TWO *(  gxz * betaxz +   gyz * betayz +   gzz * betazz)

  gxy_rhs = - TWO * alpn1 * Axy    +  F1o3 * gxy    * div_beta       + &
                      gxx * betaxy                  +   gxz * betazy + &
                                       gyy * betayx +   gyz * betazx   &
                                                    -   gxy * betazz

  gyz_rhs = - TWO * alpn1 * Ayz    +  F1o3 * gyz    * div_beta       + &
                      gxy * betaxz +   gyy * betayz                  + &
                      gxz * betaxy                  +   gzz * betazy   &
                                                    -   gyz * betaxx
 
  gxz_rhs = - TWO * alpn1 * Axz    +  F1o3 * gxz    * div_beta       + &
                      gxx * betaxz +   gxy * betayz                  + &
                                       gyz * betayx +   gzz * betazx   &
                                                    -   gxz * betayy     !rhs for gij

! invert tilted metric
  gupzz =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz - &
           gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz
  gupxx =   ( gyy * gzz - gyz * gyz ) / gupzz
  gupxy = - ( gxy * gzz - gyz * gxz ) / gupzz
  gupxz =   ( gxy * gyz - gyy * gxz ) / gupzz
  gupyy =   ( gxx * gzz - gxz * gxz ) / gupzz
  gupyz = - ( gxx * gyz - gxy * gxz ) / gupzz
  gupzz =   ( gxx * gyy - gxy * gxy ) / gupzz

  if(co == 0)then
! Gam^i_Res = Gam^i + gup^ij_,j
  Gmx_Res = Gamx - (gupxx*(gupxx*gxxx+gupxy*gxyx+gupxz*gxzx)&
                   +gupxy*(gupxx*gxyx+gupxy*gyyx+gupxz*gyzx)&
                   +gupxz*(gupxx*gxzx+gupxy*gyzx+gupxz*gzzx)&
                   +gupxx*(gupxy*gxxy+gupyy*gxyy+gupyz*gxzy)&
                   +gupxy*(gupxy*gxyy+gupyy*gyyy+gupyz*gyzy)&
                   +gupxz*(gupxy*gxzy+gupyy*gyzy+gupyz*gzzy)&
                   +gupxx*(gupxz*gxxz+gupyz*gxyz+gupzz*gxzz)&
                   +gupxy*(gupxz*gxyz+gupyz*gyyz+gupzz*gyzz)&
                   +gupxz*(gupxz*gxzz+gupyz*gyzz+gupzz*gzzz))
  Gmy_Res = Gamy - (gupxx*(gupxy*gxxx+gupyy*gxyx+gupyz*gxzx)&
                   +gupxy*(gupxy*gxyx+gupyy*gyyx+gupyz*gyzx)&
                   +gupxz*(gupxy*gxzx+gupyy*gyzx+gupyz*gzzx)&
                   +gupxy*(gupxy*gxxy+gupyy*gxyy+gupyz*gxzy)&
                   +gupyy*(gupxy*gxyy+gupyy*gyyy+gupyz*gyzy)&
                   +gupyz*(gupxy*gxzy+gupyy*gyzy+gupyz*gzzy)&
                   +gupxy*(gupxz*gxxz+gupyz*gxyz+gupzz*gxzz)&
                   +gupyy*(gupxz*gxyz+gupyz*gyyz+gupzz*gyzz)&
                   +gupyz*(gupxz*gxzz+gupyz*gyzz+gupzz*gzzz))
  Gmz_Res = Gamz - (gupxx*(gupxz*gxxx+gupyz*gxyx+gupzz*gxzx)&
                   +gupxy*(gupxz*gxyx+gupyz*gyyx+gupzz*gyzx)&
                   +gupxz*(gupxz*gxzx+gupyz*gyzx+gupzz*gzzx)&
                   +gupxy*(gupxz*gxxy+gupyz*gxyy+gupzz*gxzy)&
                   +gupyy*(gupxz*gxyy+gupyz*gyyy+gupzz*gyzy)&
                   +gupyz*(gupxz*gxzy+gupyz*gyzy+gupzz*gzzy)&
                   +gupxz*(gupxz*gxxz+gupyz*gxyz+gupzz*gxzz)&
                   +gupyz*(gupxz*gxyz+gupyz*gyyz+gupzz*gyzz)&
                   +gupzz*(gupxz*gxzz+gupyz*gyzz+gupzz*gzzz))
  endif

! second kind of connection
  Gamxxx =HALF*( gupxx*gxxx + gupxy*(TWO*gxyx - gxxy ) + gupxz*(TWO*gxzx - gxxz ))
  Gamyxx =HALF*( gupxy*gxxx + gupyy*(TWO*gxyx - gxxy ) + gupyz*(TWO*gxzx - gxxz ))
  Gamzxx =HALF*( gupxz*gxxx + gupyz*(TWO*gxyx - gxxy ) + gupzz*(TWO*gxzx - gxxz ))
 
  Gamxyy =HALF*( gupxx*(TWO*gxyy - gyyx ) + gupxy*gyyy + gupxz*(TWO*gyzy - gyyz ))
  Gamyyy =HALF*( gupxy*(TWO*gxyy - gyyx ) + gupyy*gyyy + gupyz*(TWO*gyzy - gyyz ))
  Gamzyy =HALF*( gupxz*(TWO*gxyy - gyyx ) + gupyz*gyyy + gupzz*(TWO*gyzy - gyyz ))

  Gamxzz =HALF*( gupxx*(TWO*gxzz - gzzx ) + gupxy*(TWO*gyzz - gzzy ) + gupxz*gzzz)
  Gamyzz =HALF*( gupxy*(TWO*gxzz - gzzx ) + gupyy*(TWO*gyzz - gzzy ) + gupyz*gzzz)
  Gamzzz =HALF*( gupxz*(TWO*gxzz - gzzx ) + gupyz*(TWO*gyzz - gzzy ) + gupzz*gzzz)

  Gamxxy =HALF*( gupxx*gxxy + gupxy*gyyx + gupxz*( gxzy + gyzx - gxyz ) )
  Gamyxy =HALF*( gupxy*gxxy + gupyy*gyyx + gupyz*( gxzy + gyzx - gxyz ) )
  Gamzxy =HALF*( gupxz*gxxy + gupyz*gyyx + gupzz*( gxzy + gyzx - gxyz ) )

  Gamxxz =HALF*( gupxx*gxxz + gupxy*( gxyz + gyzx - gxzy ) + gupxz*gzzx )
  Gamyxz =HALF*( gupxy*gxxz + gupyy*( gxyz + gyzx - gxzy ) + gupyz*gzzx )
  Gamzxz =HALF*( gupxz*gxxz + gupyz*( gxyz + gyzx - gxzy ) + gupzz*gzzx )

  Gamxyz =HALF*( gupxx*( gxyz + gxzy - gyzx ) + gupxy*gyyz + gupxz*gzzy )
  Gamyyz =HALF*( gupxy*( gxyz + gxzy - gyzx ) + gupyy*gyyz + gupyz*gzzy )
  Gamzyz =HALF*( gupxz*( gxyz + gxzy - gyzx ) + gupyz*gyyz + gupzz*gzzy )
! Raise indices of \tilde A_{ij} and store in R_ij

  Rxx =    gupxx * gupxx * Axx + gupxy * gupxy * Ayy + gupxz * gupxz * Azz + &
      TWO*(gupxx * gupxy * Axy + gupxx * gupxz * Axz + gupxy * gupxz * Ayz)

  Ryy =    gupxy * gupxy * Axx + gupyy * gupyy * Ayy + gupyz * gupyz * Azz + &
      TWO*(gupxy * gupyy * Axy + gupxy * gupyz * Axz + gupyy * gupyz * Ayz)

  Rzz =    gupxz * gupxz * Axx + gupyz * gupyz * Ayy + gupzz * gupzz * Azz + &
      TWO*(gupxz * gupyz * Axy + gupxz * gupzz * Axz + gupyz * gupzz * Ayz)

  Rxy =    gupxx * gupxy * Axx + gupxy * gupyy * Ayy + gupxz * gupyz * Azz + &
          (gupxx * gupyy       + gupxy * gupxy)* Axy                       + &
          (gupxx * gupyz       + gupxz * gupxy)* Axz                       + &
          (gupxy * gupyz       + gupxz * gupyy)* Ayz

  Rxz =    gupxx * gupxz * Axx + gupxy * gupyz * Ayy + gupxz * gupzz * Azz + &
          (gupxx * gupyz       + gupxy * gupxz)* Axy                       + &
          (gupxx * gupzz       + gupxz * gupxz)* Axz                       + &
          (gupxy * gupzz       + gupxz * gupyz)* Ayz

  Ryz =    gupxy * gupxz * Axx + gupyy * gupyz * Ayy + gupyz * gupzz * Azz + &
          (gupxy * gupyz       + gupyy * gupxz)* Axy                       + &
          (gupxy * gupzz       + gupyz * gupxz)* Axz                       + &
          (gupyy * gupzz       + gupyz * gupyz)* Ayz

! Right hand side for Gam^i without shift terms...
  call fderivs(ex,Lap,Lapx,Lapy,Lapz,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)
  call fderivs(ex,trK,Kx,Ky,Kz,X,Y,Z,SYM,SYM,SYM,symmetry,Lev)

   Gamx_rhs = - TWO * (   Lapx * Rxx +   Lapy * Rxy +   Lapz * Rxz ) + &
        TWO * alpn1 * (                                                &
        -F3o2/chin1 * (   chix * Rxx +   chiy * Rxy +   chiz * Rxz ) - &
              gupxx * (   F2o3 * Kx  +  EIGHT * PI * Sx            ) - &
              gupxy * (   F2o3 * Ky  +  EIGHT * PI * Sy            ) - &
              gupxz * (   F2o3 * Kz  +  EIGHT * PI * Sz            ) + &
                        Gamxxx * Rxx + Gamxyy * Ryy + Gamxzz * Rzz   + &
                TWO * ( Gamxxy * Rxy + Gamxxz * Rxz + Gamxyz * Ryz ) )

   Gamy_rhs = - TWO * (   Lapx * Rxy +   Lapy * Ryy +   Lapz * Ryz ) + &
        TWO * alpn1 * (                                                &
        -F3o2/chin1 * (   chix * Rxy +  chiy * Ryy +    chiz * Ryz ) - &
              gupxy * (   F2o3 * Kx  +  EIGHT * PI * Sx            ) - &
              gupyy * (   F2o3 * Ky  +  EIGHT * PI * Sy            ) - &
              gupyz * (   F2o3 * Kz  +  EIGHT * PI * Sz            ) + &
                        Gamyxx * Rxx + Gamyyy * Ryy + Gamyzz * Rzz   + &
                TWO * ( Gamyxy * Rxy + Gamyxz * Rxz + Gamyyz * Ryz ) )

   Gamz_rhs = - TWO * (   Lapx * Rxz +   Lapy * Ryz +   Lapz * Rzz ) + &
        TWO * alpn1 * (                                                &
        -F3o2/chin1 * (   chix * Rxz +  chiy * Ryz +    chiz * Rzz ) - &
              gupxz * (   F2o3 * Kx  +  EIGHT * PI * Sx            ) - &
              gupyz * (   F2o3 * Ky  +  EIGHT * PI * Sy            ) - &
              gupzz * (   F2o3 * Kz  +  EIGHT * PI * Sz            ) + &
                        Gamzxx * Rxx + Gamzyy * Ryy + Gamzzz * Rzz   + &
                TWO * ( Gamzxy * Rxy + Gamzxz * Rxz + Gamzyz * Ryz ) )

  call fdderivs(ex,betax,gxxx,gxyx,gxzx,gyyx,gyzx,gzzx,&
                X,Y,Z,ANTI,SYM, SYM ,Symmetry,Lev)
  call fdderivs(ex,betay,gxxy,gxyy,gxzy,gyyy,gyzy,gzzy,&
                X,Y,Z,SYM ,ANTI,SYM ,Symmetry,Lev)
  call fdderivs(ex,betaz,gxxz,gxyz,gxzz,gyyz,gyzz,gzzz,&
                X,Y,Z,SYM ,SYM, ANTI,Symmetry,Lev)

  fxx = gxxx + gxyy + gxzz
  fxy = gxyx + gyyy + gyzz
  fxz = gxzx + gyzy + gzzz

  Gamxa =       gupxx * Gamxxx + gupyy * Gamxyy + gupzz * Gamxzz + &
          TWO*( gupxy * Gamxxy + gupxz * Gamxxz + gupyz * Gamxyz )
  Gamya =       gupxx * Gamyxx + gupyy * Gamyyy + gupzz * Gamyzz + &
          TWO*( gupxy * Gamyxy + gupxz * Gamyxz + gupyz * Gamyyz )
  Gamza =       gupxx * Gamzxx + gupyy * Gamzyy + gupzz * Gamzzz + &
          TWO*( gupxy * Gamzxy + gupxz * Gamzxz + gupyz * Gamzyz )

  call fderivs(ex,Gamx,Gamxx,Gamxy,Gamxz,X,Y,Z,ANTI,SYM ,SYM ,Symmetry,Lev)
  call fderivs(ex,Gamy,Gamyx,Gamyy,Gamyz,X,Y,Z,SYM ,ANTI,SYM ,Symmetry,Lev)
  call fderivs(ex,Gamz,Gamzx,Gamzy,Gamzz,X,Y,Z,SYM ,SYM ,ANTI,Symmetry,Lev)

  Gamx_rhs =               Gamx_rhs +  F2o3 *  Gamxa * div_beta        - &
                     Gamxa * betaxx - Gamya * betaxy - Gamza * betaxz  + &
             F1o3 * (gupxx * fxx    + gupxy * fxy    + gupxz * fxz    ) + &
                     gupxx * gxxx   + gupyy * gyyx   + gupzz * gzzx    + &
              TWO * (gupxy * gxyx   + gupxz * gxzx   + gupyz * gyzx  )

  Gamy_rhs =               Gamy_rhs +  F2o3 *  Gamya * div_beta        - &
                     Gamxa * betayx - Gamya * betayy - Gamza * betayz  + &
             F1o3 * (gupxy * fxx    + gupyy * fxy    + gupyz * fxz    ) + &
                     gupxx * gxxy   + gupyy * gyyy   + gupzz * gzzy    + &
              TWO * (gupxy * gxyy   + gupxz * gxzy   + gupyz * gyzy  )

  Gamz_rhs =               Gamz_rhs +  F2o3 *  Gamza * div_beta        - &
                     Gamxa * betazx - Gamya * betazy - Gamza * betazz  + &
             F1o3 * (gupxz * fxx    + gupyz * fxy    + gupzz * fxz    ) + &
                     gupxx * gxxz   + gupyy * gyyz   + gupzz * gzzz    + &
              TWO * (gupxy * gxyz   + gupxz * gxzz   + gupyz * gyzz  )    !rhs for Gam^i

!first kind of connection stored in gij,k
  gxxx = gxx * Gamxxx + gxy * Gamyxx + gxz * Gamzxx
  gxyx = gxx * Gamxxy + gxy * Gamyxy + gxz * Gamzxy
  gxzx = gxx * Gamxxz + gxy * Gamyxz + gxz * Gamzxz
  gyyx = gxx * Gamxyy + gxy * Gamyyy + gxz * Gamzyy
  gyzx = gxx * Gamxyz + gxy * Gamyyz + gxz * Gamzyz
  gzzx = gxx * Gamxzz + gxy * Gamyzz + gxz * Gamzzz

  gxxy = gxy * Gamxxx + gyy * Gamyxx + gyz * Gamzxx
  gxyy = gxy * Gamxxy + gyy * Gamyxy + gyz * Gamzxy
  gxzy = gxy * Gamxxz + gyy * Gamyxz + gyz * Gamzxz
  gyyy = gxy * Gamxyy + gyy * Gamyyy + gyz * Gamzyy
  gyzy = gxy * Gamxyz + gyy * Gamyyz + gyz * Gamzyz
  gzzy = gxy * Gamxzz + gyy * Gamyzz + gyz * Gamzzz

  gxxz = gxz * Gamxxx + gyz * Gamyxx + gzz * Gamzxx
  gxyz = gxz * Gamxxy + gyz * Gamyxy + gzz * Gamzxy
  gxzz = gxz * Gamxxz + gyz * Gamyxz + gzz * Gamzxz
  gyyz = gxz * Gamxyy + gyz * Gamyyy + gzz * Gamzyy
  gyzz = gxz * Gamxyz + gyz * Gamyyz + gzz * Gamzyz
  gzzz = gxz * Gamxzz + gyz * Gamyzz + gzz * Gamzzz

!compute Ricci tensor for tilted metric
   call fdderivs(ex,dxx,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM ,SYM ,SYM ,symmetry,Lev)
   Rxx =   gupxx * fxx + gupyy * fyy + gupzz * fzz + &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) * TWO

   call fdderivs(ex,dyy,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM ,SYM ,SYM ,symmetry,Lev)
   Ryy =   gupxx * fxx + gupyy * fyy + gupzz * fzz + &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) * TWO

   call fdderivs(ex,dzz,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM ,SYM ,SYM ,symmetry,Lev)
   Rzz =   gupxx * fxx + gupyy * fyy + gupzz * fzz + &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) * TWO

   call fdderivs(ex,gxy,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,ANTI, ANTI,SYM ,symmetry,Lev)
   Rxy =   gupxx * fxx + gupyy * fyy + gupzz * fzz + &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) * TWO

   call fdderivs(ex,gxz,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,ANTI ,SYM ,ANTI,symmetry,Lev)
   Rxz =   gupxx * fxx + gupyy * fyy + gupzz * fzz + &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) * TWO

   call fdderivs(ex,gyz,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM ,ANTI ,ANTI,symmetry,Lev)
   Ryz =   gupxx * fxx + gupyy * fyy + gupzz * fzz + &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) * TWO

  Rxx =     - HALF * Rxx                                   + &
               gxx * Gamxx+ gxy * Gamyx   +    gxz * Gamzx + &
             Gamxa * gxxx +  Gamya * gxyx +  Gamza * gxzx  + &
   gupxx *(                                                  &
       TWO*(Gamxxx * gxxx + Gamyxx * gxyx + Gamzxx * gxzx) + &
            Gamxxx * gxxx + Gamyxx * gxxy + Gamzxx * gxxz )+ &
   gupxy *(                                                  &
       TWO*(Gamxxx * gxyx + Gamyxx * gyyx + Gamzxx * gyzx  + &
            Gamxxy * gxxx + Gamyxy * gxyx + Gamzxy * gxzx) + &
            Gamxxy * gxxx + Gamyxy * gxxy + Gamzxy * gxxz  + &
            Gamxxx * gxyx + Gamyxx * gxyy + Gamzxx * gxyz )+ &
   gupxz *(                                                  &
       TWO*(Gamxxx * gxzx + Gamyxx * gyzx + Gamzxx * gzzx  + &
            Gamxxz * gxxx + Gamyxz * gxyx + Gamzxz * gxzx) + &
            Gamxxz * gxxx + Gamyxz * gxxy + Gamzxz * gxxz  + &
            Gamxxx * gxzx + Gamyxx * gxzy + Gamzxx * gxzz )+ &
   gupyy *(                                                  &
       TWO*(Gamxxy * gxyx + Gamyxy * gyyx + Gamzxy * gyzx) + &
            Gamxxy * gxyx + Gamyxy * gxyy + Gamzxy * gxyz )+ &
   gupyz *(                                                  &
       TWO*(Gamxxy * gxzx + Gamyxy * gyzx + Gamzxy * gzzx  + &
            Gamxxz * gxyx + Gamyxz * gyyx + Gamzxz * gyzx) + &
            Gamxxz * gxyx + Gamyxz * gxyy + Gamzxz * gxyz  + &
            Gamxxy * gxzx + Gamyxy * gxzy + Gamzxy * gxzz )+ &
   gupzz *(                                                  &
       TWO*(Gamxxz * gxzx + Gamyxz * gyzx + Gamzxz * gzzx) + &
            Gamxxz * gxzx + Gamyxz * gxzy + Gamzxz * gxzz )

  Ryy =     - HALF * Ryy                                   + &
               gxy * Gamxy+  gyy * Gamyy  +  gyz * Gamzy   + &
             Gamxa * gxyy +  Gamya * gyyy +  Gamza * gyzy  + &
   gupxx *(                                                  &
       TWO*(Gamxxy * gxxy + Gamyxy * gxyy + Gamzxy * gxzy) + &
            Gamxxy * gxyx + Gamyxy * gxyy + Gamzxy * gxyz )+ &
   gupxy *(                                                  &
       TWO*(Gamxxy * gxyy + Gamyxy * gyyy + Gamzxy * gyzy  + &
            Gamxyy * gxxy + Gamyyy * gxyy + Gamzyy * gxzy) + &
            Gamxyy * gxyx + Gamyyy * gxyy + Gamzyy * gxyz  + &
            Gamxxy * gyyx + Gamyxy * gyyy + Gamzxy * gyyz )+ &
   gupxz *(                                                  &
       TWO*(Gamxxy * gxzy + Gamyxy * gyzy + Gamzxy * gzzy  + &
            Gamxyz * gxxy + Gamyyz * gxyy + Gamzyz * gxzy) + &
            Gamxyz * gxyx + Gamyyz * gxyy + Gamzyz * gxyz  + &
            Gamxxy * gyzx + Gamyxy * gyzy + Gamzxy * gyzz )+ &
   gupyy *(                                                  &
       TWO*(Gamxyy * gxyy + Gamyyy * gyyy + Gamzyy * gyzy) + &
            Gamxyy * gyyx + Gamyyy * gyyy + Gamzyy * gyyz )+ &
   gupyz *(                                                  &
       TWO*(Gamxyy * gxzy + Gamyyy * gyzy + Gamzyy * gzzy  + &
            Gamxyz * gxyy + Gamyyz * gyyy + Gamzyz * gyzy) + &
            Gamxyz * gyyx + Gamyyz * gyyy + Gamzyz * gyyz  + &
            Gamxyy * gyzx + Gamyyy * gyzy + Gamzyy * gyzz )+ &
   gupzz *(                                                  &
       TWO*(Gamxyz * gxzy + Gamyyz * gyzy + Gamzyz * gzzy) + &
            Gamxyz * gyzx + Gamyyz * gyzy + Gamzyz * gyzz )

  Rzz =     - HALF * Rzz                                   + &
               gxz * Gamxz+ gyz * Gamyz  +    gzz * Gamzz  + &
             Gamxa * gxzz +  Gamya * gyzz +  Gamza * gzzz  + &
   gupxx *(                                                  &
       TWO*(Gamxxz * gxxz + Gamyxz * gxyz + Gamzxz * gxzz) + &
            Gamxxz * gxzx + Gamyxz * gxzy + Gamzxz * gxzz )+ &
   gupxy *(                                                  &
       TWO*(Gamxxz * gxyz + Gamyxz * gyyz + Gamzxz * gyzz  + &
            Gamxyz * gxxz + Gamyyz * gxyz + Gamzyz * gxzz) + &
            Gamxyz * gxzx + Gamyyz * gxzy + Gamzyz * gxzz  + &
            Gamxxz * gyzx + Gamyxz * gyzy + Gamzxz * gyzz )+ &
   gupxz *(                                                  &
       TWO*(Gamxxz * gxzz + Gamyxz * gyzz + Gamzxz * gzzz  + &
            Gamxzz * gxxz + Gamyzz * gxyz + Gamzzz * gxzz) + &
            Gamxzz * gxzx + Gamyzz * gxzy + Gamzzz * gxzz  + &
            Gamxxz * gzzx + Gamyxz * gzzy + Gamzxz * gzzz )+ &
   gupyy *(                                                  &
       TWO*(Gamxyz * gxyz + Gamyyz * gyyz + Gamzyz * gyzz) + &
            Gamxyz * gyzx + Gamyyz * gyzy + Gamzyz * gyzz )+ &
   gupyz *(                                                  &
       TWO*(Gamxyz * gxzz + Gamyyz * gyzz + Gamzyz * gzzz  + &
            Gamxzz * gxyz + Gamyzz * gyyz + Gamzzz * gyzz) + &
            Gamxzz * gyzx + Gamyzz * gyzy + Gamzzz * gyzz  + &
            Gamxyz * gzzx + Gamyyz * gzzy + Gamzyz * gzzz )+ &
   gupzz *(                                                  &
       TWO*(Gamxzz * gxzz + Gamyzz * gyzz + Gamzzz * gzzz) + &
            Gamxzz * gzzx + Gamyzz * gzzy + Gamzzz * gzzz )

  Rxy = HALF*(     - Rxy                                   + &
               gxx * Gamxy +    gxy * Gamyy + gxz * Gamzy  + &
               gxy * Gamxx +    gyy * Gamyx + gyz * Gamzx  + &
             Gamxa * gxyx +  Gamya * gyyx +  Gamza * gyzx  + &
             Gamxa * gxxy +  Gamya * gxyy +  Gamza * gxzy )+ &
   gupxx *(                                                  &
            Gamxxx * gxxy + Gamyxx * gxyy + Gamzxx * gxzy  + &
            Gamxxy * gxxx + Gamyxy * gxyx + Gamzxy * gxzx  + &
            Gamxxx * gxyx + Gamyxx * gxyy + Gamzxx * gxyz )+ &
   gupxy *(                                                  &
            Gamxxx * gxyy + Gamyxx * gyyy + Gamzxx * gyzy  + &
            Gamxxy * gxyx + Gamyxy * gyyx + Gamzxy * gyzx  + &
            Gamxxy * gxyx + Gamyxy * gxyy + Gamzxy * gxyz  + &
            Gamxxy * gxxy + Gamyxy * gxyy + Gamzxy * gxzy  + &
            Gamxyy * gxxx + Gamyyy * gxyx + Gamzyy * gxzx  + &
            Gamxxx * gyyx + Gamyxx * gyyy + Gamzxx * gyyz )+ &
   gupxz *(                                                  &
            Gamxxx * gxzy + Gamyxx * gyzy + Gamzxx * gzzy  + &
            Gamxxy * gxzx + Gamyxy * gyzx + Gamzxy * gzzx  + &
            Gamxxz * gxyx + Gamyxz * gxyy + Gamzxz * gxyz  + &
            Gamxxz * gxxy + Gamyxz * gxyy + Gamzxz * gxzy  + &
            Gamxyz * gxxx + Gamyyz * gxyx + Gamzyz * gxzx  + &
            Gamxxx * gyzx + Gamyxx * gyzy + Gamzxx * gyzz )+ &
   gupyy *(                                                  &
            Gamxxy * gxyy + Gamyxy * gyyy + Gamzxy * gyzy  + &
            Gamxyy * gxyx + Gamyyy * gyyx + Gamzyy * gyzx  + &
            Gamxxy * gyyx + Gamyxy * gyyy + Gamzxy * gyyz )+ &
   gupyz *(                                                  &
            Gamxxy * gxzy + Gamyxy * gyzy + Gamzxy * gzzy  + &
            Gamxyy * gxzx + Gamyyy * gyzx + Gamzyy * gzzx  + &
            Gamxxz * gyyx + Gamyxz * gyyy + Gamzxz * gyyz  + &
            Gamxxz * gxyy + Gamyxz * gyyy + Gamzxz * gyzy  + &
            Gamxyz * gxyx + Gamyyz * gyyx + Gamzyz * gyzx  + &
            Gamxxy * gyzx + Gamyxy * gyzy + Gamzxy * gyzz )+ &
   gupzz *(                                                  &
            Gamxxz * gxzy + Gamyxz * gyzy + Gamzxz * gzzy  + &
            Gamxyz * gxzx + Gamyyz * gyzx + Gamzyz * gzzx  + &
            Gamxxz * gyzx + Gamyxz * gyzy + Gamzxz * gyzz )

  Rxz = HALF*(     - Rxz                                   + &
               gxx * Gamxz +  gxy * Gamyz + gxz * Gamzz    + &
               gxz * Gamxx +  gyz * Gamyx + gzz * Gamzx    + &
             Gamxa * gxzx +  Gamya * gyzx +  Gamza * gzzx  + &
             Gamxa * gxxz +  Gamya * gxyz +  Gamza * gxzz )+ &
   gupxx *(                                                  &
            Gamxxx * gxxz + Gamyxx * gxyz + Gamzxx * gxzz  + &
            Gamxxz * gxxx + Gamyxz * gxyx + Gamzxz * gxzx  + &
            Gamxxx * gxzx + Gamyxx * gxzy + Gamzxx * gxzz )+ &
   gupxy *(                                                  &
            Gamxxx * gxyz + Gamyxx * gyyz + Gamzxx * gyzz  + &
            Gamxxz * gxyx + Gamyxz * gyyx + Gamzxz * gyzx  + &
            Gamxxy * gxzx + Gamyxy * gxzy + Gamzxy * gxzz  + &
            Gamxxy * gxxz + Gamyxy * gxyz + Gamzxy * gxzz  + &
            Gamxyz * gxxx + Gamyyz * gxyx + Gamzyz * gxzx  + &
            Gamxxx * gyzx + Gamyxx * gyzy + Gamzxx * gyzz )+ &
   gupxz *(                                                  &
            Gamxxx * gxzz + Gamyxx * gyzz + Gamzxx * gzzz  + &
            Gamxxz * gxzx + Gamyxz * gyzx + Gamzxz * gzzx  + &
            Gamxxz * gxzx + Gamyxz * gxzy + Gamzxz * gxzz  + &
            Gamxxz * gxxz + Gamyxz * gxyz + Gamzxz * gxzz  + &
            Gamxzz * gxxx + Gamyzz * gxyx + Gamzzz * gxzx  + &
            Gamxxx * gzzx + Gamyxx * gzzy + Gamzxx * gzzz )+ &
   gupyy *(                                                  &
            Gamxxy * gxyz + Gamyxy * gyyz + Gamzxy * gyzz  + &
            Gamxyz * gxyx + Gamyyz * gyyx + Gamzyz * gyzx  + &
            Gamxxy * gyzx + Gamyxy * gyzy + Gamzxy * gyzz )+ &
   gupyz *(                                                  &
            Gamxxy * gxzz + Gamyxy * gyzz + Gamzxy * gzzz  + &
            Gamxyz * gxzx + Gamyyz * gyzx + Gamzyz * gzzx  + &
            Gamxxz * gyzx + Gamyxz * gyzy + Gamzxz * gyzz  + &
            Gamxxz * gxyz + Gamyxz * gyyz + Gamzxz * gyzz  + &
            Gamxzz * gxyx + Gamyzz * gyyx + Gamzzz * gyzx  + &
            Gamxxy * gzzx + Gamyxy * gzzy + Gamzxy * gzzz )+ &
   gupzz *(                                                  &
            Gamxxz * gxzz + Gamyxz * gyzz + Gamzxz * gzzz  + &
            Gamxzz * gxzx + Gamyzz * gyzx + Gamzzz * gzzx  + &
            Gamxxz * gzzx + Gamyxz * gzzy + Gamzxz * gzzz )

  Ryz = HALF*(     - Ryz                                   + &
               gxy * Gamxz + gyy * Gamyz + gyz * Gamzz     + &
               gxz * Gamxy + gyz * Gamyy + gzz * Gamzy     + &
             Gamxa * gxzy +  Gamya * gyzy +  Gamza * gzzy  + &
             Gamxa * gxyz +  Gamya * gyyz +  Gamza * gyzz )+ &
   gupxx *(                                                  &
            Gamxxy * gxxz + Gamyxy * gxyz + Gamzxy * gxzz  + &
            Gamxxz * gxxy + Gamyxz * gxyy + Gamzxz * gxzy  + &
            Gamxxy * gxzx + Gamyxy * gxzy + Gamzxy * gxzz )+ &
   gupxy *(                                                  &
            Gamxxy * gxyz + Gamyxy * gyyz + Gamzxy * gyzz  + &
            Gamxxz * gxyy + Gamyxz * gyyy + Gamzxz * gyzy  + &
            Gamxyy * gxzx + Gamyyy * gxzy + Gamzyy * gxzz  + &
            Gamxyy * gxxz + Gamyyy * gxyz + Gamzyy * gxzz  + &
            Gamxyz * gxxy + Gamyyz * gxyy + Gamzyz * gxzy  + &
            Gamxxy * gyzx + Gamyxy * gyzy + Gamzxy * gyzz )+ &
   gupxz *(                                                  &
            Gamxxy * gxzz + Gamyxy * gyzz + Gamzxy * gzzz  + &
            Gamxxz * gxzy + Gamyxz * gyzy + Gamzxz * gzzy  + &
            Gamxyz * gxzx + Gamyyz * gxzy + Gamzyz * gxzz  + &
            Gamxyz * gxxz + Gamyyz * gxyz + Gamzyz * gxzz  + &
            Gamxzz * gxxy + Gamyzz * gxyy + Gamzzz * gxzy  + &
            Gamxxy * gzzx + Gamyxy * gzzy + Gamzxy * gzzz )+ &
   gupyy *(                                                  &
            Gamxyy * gxyz + Gamyyy * gyyz + Gamzyy * gyzz  + &
            Gamxyz * gxyy + Gamyyz * gyyy + Gamzyz * gyzy  + &
            Gamxyy * gyzx + Gamyyy * gyzy + Gamzyy * gyzz )+ &
   gupyz *(                                                  &
            Gamxyy * gxzz + Gamyyy * gyzz + Gamzyy * gzzz  + &
            Gamxyz * gxzy + Gamyyz * gyzy + Gamzyz * gzzy  + &
            Gamxyz * gyzx + Gamyyz * gyzy + Gamzyz * gyzz  + &
            Gamxyz * gxyz + Gamyyz * gyyz + Gamzyz * gyzz  + &
            Gamxzz * gxyy + Gamyzz * gyyy + Gamzzz * gyzy  + &
            Gamxyy * gzzx + Gamyyy * gzzy + Gamzyy * gzzz )+ &
   gupzz *(                                                  &
            Gamxyz * gxzz + Gamyyz * gyzz + Gamzyz * gzzz  + &
            Gamxzz * gxzy + Gamyzz * gyzy + Gamzzz * gzzy  + &
            Gamxyz * gzzx + Gamyyz * gzzy + Gamzyz * gzzz )
!covariant second derivative of chi respect to tilted metric
  call fdderivs(ex,chi,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)

  fxx = fxx - Gamxxx * chix - Gamyxx * chiy - Gamzxx * chiz
  fxy = fxy - Gamxxy * chix - Gamyxy * chiy - Gamzxy * chiz
  fxz = fxz - Gamxxz * chix - Gamyxz * chiy - Gamzxz * chiz
  fyy = fyy - Gamxyy * chix - Gamyyy * chiy - Gamzyy * chiz
  fyz = fyz - Gamxyz * chix - Gamyyz * chiy - Gamzyz * chiz
  fzz = fzz - Gamxzz * chix - Gamyzz * chiy - Gamzzz * chiz
! Store D^l D_l chi - 3/(2*chi) D^l chi D_l chi in f

  f =        gupxx * ( fxx - F3o2/chin1 * chix * chix ) + &
             gupyy * ( fyy - F3o2/chin1 * chiy * chiy ) + &
             gupzz * ( fzz - F3o2/chin1 * chiz * chiz ) + &
       TWO * gupxy * ( fxy - F3o2/chin1 * chix * chiy ) + &
       TWO * gupxz * ( fxz - F3o2/chin1 * chix * chiz ) + &
       TWO * gupyz * ( fyz - F3o2/chin1 * chiy * chiz ) 
! Add chi part to Ricci tensor:

  Rxx = Rxx + (fxx - chix*chix/chin1/TWO + gxx * f)/chin1/TWO
  Ryy = Ryy + (fyy - chiy*chiy/chin1/TWO + gyy * f)/chin1/TWO
  Rzz = Rzz + (fzz - chiz*chiz/chin1/TWO + gzz * f)/chin1/TWO
  Rxy = Rxy + (fxy - chix*chiy/chin1/TWO + gxy * f)/chin1/TWO
  Rxz = Rxz + (fxz - chix*chiz/chin1/TWO + gxz * f)/chin1/TWO
  Ryz = Ryz + (fyz - chiy*chiz/chin1/TWO + gyz * f)/chin1/TWO

! covariant second derivatives of the lapse respect to physical metric
  call fdderivs(ex,Lap,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z, &
                SYM,SYM,SYM,symmetry,Lev)

  gxxx = (gupxx * chix + gupxy * chiy + gupxz * chiz)/chin1
  gxxy = (gupxy * chix + gupyy * chiy + gupyz * chiz)/chin1
  gxxz = (gupxz * chix + gupyz * chiy + gupzz * chiz)/chin1
! now get physical second kind of connection
  Gamxxx = Gamxxx - ( (chix + chix)/chin1 - gxx * gxxx )*HALF
  Gamyxx = Gamyxx - (                     - gxx * gxxy )*HALF
  Gamzxx = Gamzxx - (                     - gxx * gxxz )*HALF
  Gamxyy = Gamxyy - (                     - gyy * gxxx )*HALF
  Gamyyy = Gamyyy - ( (chiy + chiy)/chin1 - gyy * gxxy )*HALF
  Gamzyy = Gamzyy - (                     - gyy * gxxz )*HALF
  Gamxzz = Gamxzz - (                     - gzz * gxxx )*HALF
  Gamyzz = Gamyzz - (                     - gzz * gxxy )*HALF
  Gamzzz = Gamzzz - ( (chiz + chiz)/chin1 - gzz * gxxz )*HALF
  Gamxxy = Gamxxy - (  chiy        /chin1 - gxy * gxxx )*HALF
  Gamyxy = Gamyxy - (         chix /chin1 - gxy * gxxy )*HALF
  Gamzxy = Gamzxy - (                     - gxy * gxxz )*HALF
  Gamxxz = Gamxxz - (  chiz        /chin1 - gxz * gxxx )*HALF
  Gamyxz = Gamyxz - (                     - gxz * gxxy )*HALF
  Gamzxz = Gamzxz - (         chix /chin1 - gxz * gxxz )*HALF
  Gamxyz = Gamxyz - (                     - gyz * gxxx )*HALF
  Gamyyz = Gamyyz - (  chiz        /chin1 - gyz * gxxy )*HALF
  Gamzyz = Gamzyz - (         chiy /chin1 - gyz * gxxz )*HALF

  fxx = fxx - Gamxxx*Lapx - Gamyxx*Lapy - Gamzxx*Lapz
  fyy = fyy - Gamxyy*Lapx - Gamyyy*Lapy - Gamzyy*Lapz
  fzz = fzz - Gamxzz*Lapx - Gamyzz*Lapy - Gamzzz*Lapz
  fxy = fxy - Gamxxy*Lapx - Gamyxy*Lapy - Gamzxy*Lapz
  fxz = fxz - Gamxxz*Lapx - Gamyxz*Lapy - Gamzxz*Lapz
  fyz = fyz - Gamxyz*Lapx - Gamyyz*Lapy - Gamzyz*Lapz

! store D^i D_i Lap in trK_rhs upto chi
  trK_rhs =    gupxx * fxx + gupyy * fyy + gupzz * fzz + &
        TWO* ( gupxy * fxy + gupxz * fxz + gupyz * fyz )
#if 1        
!! follow bam code
  S =  chin1 * ( gupxx * Sxx + gupyy * Syy + gupzz * Szz + &
     TWO * ( gupxy * Sxy + gupxz * Sxz + gupyz * Syz ) )
  f = F2o3 * trK * trK -(&
       gupxx * ( &
       gupxx * Axx * Axx + gupyy * Axy * Axy + gupzz * Axz * Axz + &
       TWO * (gupxy * Axx * Axy + gupxz * Axx * Axz + gupyz * Axy * Axz) ) + &
       gupyy * ( &
       gupxx * Axy * Axy + gupyy * Ayy * Ayy + gupzz * Ayz * Ayz + &
       TWO * (gupxy * Axy * Ayy + gupxz * Axy * Ayz + gupyz * Ayy * Ayz) ) + &
       gupzz * ( &
       gupxx * Axz * Axz + gupyy * Ayz * Ayz + gupzz * Azz * Azz + &
       TWO * (gupxy * Axz * Ayz + gupxz * Axz * Azz + gupyz * Ayz * Azz) ) + &
       TWO * ( &
       gupxy * ( &
       gupxx * Axx * Axy + gupyy * Axy * Ayy + gupzz * Axz * Ayz + &
       gupxy * (Axx * Ayy + Axy * Axy) + &
       gupxz * (Axx * Ayz + Axz * Axy) + &
       gupyz * (Axy * Ayz + Axz * Ayy) ) + &
       gupxz * ( &
       gupxx * Axx * Axz + gupyy * Axy * Ayz + gupzz * Axz * Azz + &
       gupxy * (Axx * Ayz + Axy * Axz) + &
       gupxz * (Axx * Azz + Axz * Axz) + &
       gupyz * (Axy * Azz + Axz * Ayz) ) + &
       gupyz * ( &
       gupxx * Axy * Axz + gupyy * Ayy * Ayz + gupzz * Ayz * Azz + &
       gupxy * (Axy * Ayz + Ayy * Axz) + &
       gupxz * (Axy * Azz + Ayz * Axz) + &
       gupyz * (Ayy * Azz + Ayz * Ayz) ) )) -1.6d1*PI*rho + EIGHT * PI * S
  f = - F1o3 *(  gupxx * fxx + gupyy * fyy + gupzz * fzz + &
        TWO* ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) + alpn1/chin1*f)
  
  fxx = alpn1 * (Rxx - EIGHT * PI * Sxx) - fxx
  fxy = alpn1 * (Rxy - EIGHT * PI * Sxy) - fxy
  fxz = alpn1 * (Rxz - EIGHT * PI * Sxz) - fxz
  fyy = alpn1 * (Ryy - EIGHT * PI * Syy) - fyy
  fyz = alpn1 * (Ryz - EIGHT * PI * Syz) - fyz
  fzz = alpn1 * (Rzz - EIGHT * PI * Szz) - fzz
#else        
! Add lapse and S_ij parts to Ricci tensor:

  fxx = alpn1 * (Rxx - EIGHT * PI * Sxx) - fxx
  fxy = alpn1 * (Rxy - EIGHT * PI * Sxy) - fxy
  fxz = alpn1 * (Rxz - EIGHT * PI * Sxz) - fxz
  fyy = alpn1 * (Ryy - EIGHT * PI * Syy) - fyy
  fyz = alpn1 * (Ryz - EIGHT * PI * Syz) - fyz
  fzz = alpn1 * (Rzz - EIGHT * PI * Szz) - fzz

! Compute trace-free part (note: chi^-1 and chi cancel!):

  f = F1o3 *(  gupxx * fxx + gupyy * fyy + gupzz * fzz + &
        TWO* ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) )
#endif

  Axx_rhs = fxx - gxx * f
  Ayy_rhs = fyy - gyy * f
  Azz_rhs = fzz - gzz * f
  Axy_rhs = fxy - gxy * f
  Axz_rhs = fxz - gxz * f
  Ayz_rhs = fyz - gyz * f

! Now: store A_il A^l_j into fij:

  fxx =       gupxx * Axx * Axx + gupyy * Axy * Axy + gupzz * Axz * Axz + &
       TWO * (gupxy * Axx * Axy + gupxz * Axx * Axz + gupyz * Axy * Axz)
  fyy =       gupxx * Axy * Axy + gupyy * Ayy * Ayy + gupzz * Ayz * Ayz + &
       TWO * (gupxy * Axy * Ayy + gupxz * Axy * Ayz + gupyz * Ayy * Ayz)
  fzz =       gupxx * Axz * Axz + gupyy * Ayz * Ayz + gupzz * Azz * Azz + &
       TWO * (gupxy * Axz * Ayz + gupxz * Axz * Azz + gupyz * Ayz * Azz)
  fxy =       gupxx * Axx * Axy + gupyy * Axy * Ayy + gupzz * Axz * Ayz + &
              gupxy *(Axx * Ayy + Axy * Axy)                            + &
              gupxz *(Axx * Ayz + Axz * Axy)                            + &
              gupyz *(Axy * Ayz + Axz * Ayy)
  fxz =       gupxx * Axx * Axz + gupyy * Axy * Ayz + gupzz * Axz * Azz + &
              gupxy *(Axx * Ayz + Axy * Axz)                            + &
              gupxz *(Axx * Azz + Axz * Axz)                            + &
              gupyz *(Axy * Azz + Axz * Ayz)
  fyz =       gupxx * Axy * Axz + gupyy * Ayy * Ayz + gupzz * Ayz * Azz + &
              gupxy *(Axy * Ayz + Ayy * Axz)                            + &
              gupxz *(Axy * Azz + Ayz * Axz)                            + &
              gupyz *(Ayy * Azz + Ayz * Ayz)

  f = chin1
! store D^i D_i Lap in trK_rhs
  trK_rhs = f*trK_rhs
          
  Axx_rhs =           f * Axx_rhs+ alpn1 * (trK * Axx - TWO * fxx)  + &
           TWO * (  Axx * betaxx +   Axy * betayx +   Axz * betazx )- &
             F2o3 * Axx * div_beta

  Ayy_rhs =           f * Ayy_rhs+ alpn1 * (trK * Ayy - TWO * fyy)  + &
           TWO * (  Axy * betaxy +   Ayy * betayy +   Ayz * betazy )- &
             F2o3 * Ayy * div_beta

  Azz_rhs =           f * Azz_rhs+ alpn1 * (trK * Azz - TWO * fzz)  + &
           TWO * (  Axz * betaxz +   Ayz * betayz +   Azz * betazz )- &
             F2o3 * Azz * div_beta

  Axy_rhs =           f * Axy_rhs+ alpn1 *( trK * Axy  - TWO * fxy )+ &
                    Axx * betaxy                  +   Axz * betazy  + &
                                     Ayy * betayx +   Ayz * betazx  + &
             F1o3 * Axy * div_beta                -   Axy * betazz

  Ayz_rhs =           f * Ayz_rhs+ alpn1 *( trK * Ayz  - TWO * fyz )+ &
                    Axy * betaxz +   Ayy * betayz                   + &
                    Axz * betaxy                  +   Azz * betazy  + &
             F1o3 * Ayz * div_beta                -   Ayz * betaxx
 
  Axz_rhs =           f * Axz_rhs+ alpn1 *( trK * Axz  - TWO * fxz )+ &
                    Axx * betaxz +   Axy * betayz                   + &
                                     Ayz * betayx +   Azz * betazx  + &
             F1o3 * Axz * div_beta                -   Axz * betayy      !rhs for Aij

! Compute trace of S_ij

  S =  f * ( gupxx * Sxx + gupyy * Syy + gupzz * Szz + &
     TWO * ( gupxy * Sxy + gupxz * Sxz + gupyz * Syz ) )

  trK_rhs = - trK_rhs + alpn1 *( F1o3 * trK * trK         + &
                gupxx * fxx + gupyy * fyy + gupzz * fzz   + &
        TWO * ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) + &
       FOUR * PI * ( rho + S ))                                !rhs for trK
  
!!!! gauge variable part

  Lap_rhs = -TWO*alpn1*trK
#if (GAUGE == 0)
  betax_rhs = FF*dtSfx
  betay_rhs = FF*dtSfy
  betaz_rhs = FF*dtSfz

  dtSfx_rhs = Gamx_rhs - eta*dtSfx
  dtSfy_rhs = Gamy_rhs - eta*dtSfy
  dtSfz_rhs = Gamz_rhs - eta*dtSfz
#elif (GAUGE == 1)
  betax_rhs = Gamx - eta*betax
  betay_rhs = Gamy - eta*betay
  betaz_rhs = Gamz - eta*betaz

  dtSfx_rhs = ZEO
  dtSfy_rhs = ZEO
  dtSfz_rhs = ZEO
#elif (GAUGE == 2)
  betax_rhs = FF*dtSfx
  betay_rhs = FF*dtSfy
  betaz_rhs = FF*dtSfz

  call fderivs(ex,chi,dtSfx_rhs,dtSfy_rhs,dtSfz_rhs,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)
  reta = gupxx * dtSfx_rhs * dtSfx_rhs + gupyy * dtSfy_rhs * dtSfy_rhs + gupzz * dtSfz_rhs * dtSfz_rhs + &
       TWO * (gupxy * dtSfx_rhs * dtSfy_rhs + gupxz * dtSfx_rhs * dtSfz_rhs + gupyz * dtSfy_rhs * dtSfz_rhs)
  reta = 1.31d0/2*dsqrt(reta/chin1)/(1-dsqrt(chin1))**2
  dtSfx_rhs = Gamx_rhs - reta*dtSfx
  dtSfy_rhs = Gamy_rhs - reta*dtSfy
  dtSfz_rhs = Gamz_rhs - reta*dtSfz
#elif (GAUGE == 3)
  betax_rhs = FF*dtSfx
  betay_rhs = FF*dtSfy
  betaz_rhs = FF*dtSfz

  call fderivs(ex,chi,dtSfx_rhs,dtSfy_rhs,dtSfz_rhs,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)
  reta = gupxx * dtSfx_rhs * dtSfx_rhs + gupyy * dtSfy_rhs * dtSfy_rhs + gupzz * dtSfz_rhs * dtSfz_rhs + &
       TWO * (gupxy * dtSfx_rhs * dtSfy_rhs + gupxz * dtSfx_rhs * dtSfz_rhs + gupyz * dtSfy_rhs * dtSfz_rhs)
  reta = 1.31d0/2*dsqrt(reta/chin1)/(1-chin1)**2
  dtSfx_rhs = Gamx_rhs - reta*dtSfx
  dtSfy_rhs = Gamy_rhs - reta*dtSfy
  dtSfz_rhs = Gamz_rhs - reta*dtSfz
#elif (GAUGE == 4)
  call fderivs(ex,chi,dtSfx_rhs,dtSfy_rhs,dtSfz_rhs,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)
  reta = gupxx * dtSfx_rhs * dtSfx_rhs + gupyy * dtSfy_rhs * dtSfy_rhs + gupzz * dtSfz_rhs * dtSfz_rhs + &
       TWO * (gupxy * dtSfx_rhs * dtSfy_rhs + gupxz * dtSfx_rhs * dtSfz_rhs + gupyz * dtSfy_rhs * dtSfz_rhs)
  reta = 1.31d0/2*dsqrt(reta/chin1)/(1-dsqrt(chin1))**2
  betax_rhs = FF*Gamx - reta*betax
  betay_rhs = FF*Gamy - reta*betay
  betaz_rhs = FF*Gamz - reta*betaz

  dtSfx_rhs = ZEO
  dtSfy_rhs = ZEO
  dtSfz_rhs = ZEO
#elif (GAUGE == 5)
  call fderivs(ex,chi,dtSfx_rhs,dtSfy_rhs,dtSfz_rhs,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)
  reta = gupxx * dtSfx_rhs * dtSfx_rhs + gupyy * dtSfy_rhs * dtSfy_rhs + gupzz * dtSfz_rhs * dtSfz_rhs + &
       TWO * (gupxy * dtSfx_rhs * dtSfy_rhs + gupxz * dtSfx_rhs * dtSfz_rhs + gupyz * dtSfy_rhs * dtSfz_rhs)
  reta = 1.31d0/2*dsqrt(reta/chin1)/(1-chin1)**2
  betax_rhs = FF*Gamx - reta*betax
  betay_rhs = FF*Gamy - reta*betay
  betaz_rhs = FF*Gamz - reta*betaz

  dtSfx_rhs = ZEO
  dtSfy_rhs = ZEO
  dtSfz_rhs = ZEO
#elif (GAUGE == 6)
  if(BHN==2)then
   M = Mass(1)+Mass(2)
   A = 2.d0/M
   w1 = 1.2d1
   w2 = w1
   C1 = 1.d0/Mass(1) - A
   C2 = 1.d0/Mass(2) - A

   do k=1,ex(3)
   do j=1,ex(2)
   do i=1,ex(1)
     r1 = ((Porg(1)-X(i))**2+(Porg(2)-Y(j))**2+(Porg(3)-Z(k))**2)/ &
          ((Porg(1)-Porg(4))**2+(Porg(2)-Porg(5))**2+(Porg(3)-Porg(6))**2)
     r2 = ((Porg(4)-X(i))**2+(Porg(5)-Y(j))**2+(Porg(6)-Z(k))**2)/ &
          ((Porg(1)-Porg(4))**2+(Porg(2)-Porg(5))**2+(Porg(3)-Porg(6))**2)
     reta(i,j,k) = A + C1/(ONE+w1*r1) + C2/(ONE+w2*r2)
    enddo
    enddo
    enddo
  else
    write(*,*) "not support BH_num in Jason's form 1",BHN
  endif
  betax_rhs = FF*dtSfx
  betay_rhs = FF*dtSfy
  betaz_rhs = FF*dtSfz

  dtSfx_rhs = Gamx_rhs - reta*dtSfx
  dtSfy_rhs = Gamy_rhs - reta*dtSfy
  dtSfz_rhs = Gamz_rhs - reta*dtSfz
#elif (GAUGE == 7)
  if(BHN==2)then
   M = Mass(1)+Mass(2)
   A = 2.d0/M
   w1 = 1.2d1
   w2 = w1
   C1 = 1.d0/Mass(1) - A
   C2 = 1.d0/Mass(2) - A

   do k=1,ex(3)
   do j=1,ex(2)
   do i=1,ex(1)
     r1 = ((Porg(1)-X(i))**2+(Porg(2)-Y(j))**2+(Porg(3)-Z(k))**2)/ &
          ((Porg(1)-Porg(4))**2+(Porg(2)-Porg(5))**2+(Porg(3)-Porg(6))**2)
     r2 = ((Porg(4)-X(i))**2+(Porg(5)-Y(j))**2+(Porg(6)-Z(k))**2)/ &
          ((Porg(1)-Porg(4))**2+(Porg(2)-Porg(5))**2+(Porg(3)-Porg(6))**2)
     reta(i,j,k) = A + C1*dexp(-w1*r1) + C2*dexp(-w2*r2)
    enddo
    enddo
    enddo
  else
    write(*,*) "not support BH_num in Jason's form 2",BHN
  endif
  betax_rhs = FF*dtSfx
  betay_rhs = FF*dtSfy
  betaz_rhs = FF*dtSfz

  dtSfx_rhs = Gamx_rhs - reta*dtSfx
  dtSfy_rhs = Gamy_rhs - reta*dtSfy
  dtSfz_rhs = Gamz_rhs - reta*dtSfz
#endif  

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

  call lopsided(ex,X,Y,Z,gxx,gxx_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,gxy,gxy_rhs,betax,betay,betaz,Symmetry,AAS)
  call lopsided(ex,X,Y,Z,gxz,gxz_rhs,betax,betay,betaz,Symmetry,ASA)
  call lopsided(ex,X,Y,Z,gyy,gyy_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,gyz,gyz_rhs,betax,betay,betaz,Symmetry,SAA)
  call lopsided(ex,X,Y,Z,gzz,gzz_rhs,betax,betay,betaz,Symmetry,SSS)

  call lopsided(ex,X,Y,Z,Axx,Axx_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,Axy,Axy_rhs,betax,betay,betaz,Symmetry,AAS)
  call lopsided(ex,X,Y,Z,Axz,Axz_rhs,betax,betay,betaz,Symmetry,ASA)
  call lopsided(ex,X,Y,Z,Ayy,Ayy_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,Ayz,Ayz_rhs,betax,betay,betaz,Symmetry,SAA)
  call lopsided(ex,X,Y,Z,Azz,Azz_rhs,betax,betay,betaz,Symmetry,SSS)

  call lopsided(ex,X,Y,Z,chi,chi_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,trK,trK_rhs,betax,betay,betaz,Symmetry,SSS)

  call lopsided(ex,X,Y,Z,Gamx,Gamx_rhs,betax,betay,betaz,Symmetry,ASS)
  call lopsided(ex,X,Y,Z,Gamy,Gamy_rhs,betax,betay,betaz,Symmetry,SAS)
  call lopsided(ex,X,Y,Z,Gamz,Gamz_rhs,betax,betay,betaz,Symmetry,SSA)
!!
  call lopsided(ex,X,Y,Z,Lap,Lap_rhs,betax,betay,betaz,Symmetry,SSS)

#if (GAUGE == 0 || GAUGE == 1 || GAUGE == 2 || GAUGE == 3 || GAUGE == 4 || GAUGE == 5 || GAUGE == 6 || GAUGE == 7)
  call lopsided(ex,X,Y,Z,betax,betax_rhs,betax,betay,betaz,Symmetry,ASS)
  call lopsided(ex,X,Y,Z,betay,betay_rhs,betax,betay,betaz,Symmetry,SAS)
  call lopsided(ex,X,Y,Z,betaz,betaz_rhs,betax,betay,betaz,Symmetry,SSA)
#endif

#if (GAUGE == 0 || GAUGE == 2 || GAUGE == 3 || GAUGE == 6 || GAUGE == 7)
  call lopsided(ex,X,Y,Z,dtSfx,dtSfx_rhs,betax,betay,betaz,Symmetry,ASS)
  call lopsided(ex,X,Y,Z,dtSfy,dtSfy_rhs,betax,betay,betaz,Symmetry,SAS)
  call lopsided(ex,X,Y,Z,dtSfz,dtSfz_rhs,betax,betay,betaz,Symmetry,SSA)
#endif

  if(eps>0)then 
! usual Kreiss-Oliger dissipation      
  call kodis(ex,X,Y,Z,chi,chi_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,trK,trK_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,dxx,gxx_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,gxy,gxy_rhs,AAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,gxz,gxz_rhs,ASA,Symmetry,eps)
  call kodis(ex,X,Y,Z,dyy,gyy_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,gyz,gyz_rhs,SAA,Symmetry,eps)
  call kodis(ex,X,Y,Z,dzz,gzz_rhs,SSS,Symmetry,eps)
#if 0
#define i 42
#define j 40
#define k 40
if(Lev == 1)then
write(*,*) X(i),Y(j),Z(k)
write(*,*) "before",Axx_rhs(i,j,k)
endif
#undef i
#undef j
#undef k
!!stop
#endif
  call kodis(ex,X,Y,Z,Axx,Axx_rhs,SSS,Symmetry,eps)
#if 0
#define i 42
#define j 40
#define k 40
if(Lev == 1)then
write(*,*) X(i),Y(j),Z(k)
write(*,*) "after",Axx_rhs(i,j,k)
endif
#undef i
#undef j
#undef k
!!stop
#endif
  call kodis(ex,X,Y,Z,Axy,Axy_rhs,AAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Axz,Axz_rhs,ASA,Symmetry,eps)
  call kodis(ex,X,Y,Z,Ayy,Ayy_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Ayz,Ayz_rhs,SAA,Symmetry,eps)
  call kodis(ex,X,Y,Z,Azz,Azz_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Gamx,Gamx_rhs,ASS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Gamy,Gamy_rhs,SAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Gamz,Gamz_rhs,SSA,Symmetry,eps)

#if 1 
!! bam does not apply dissipation on gauge variables
  call kodis(ex,X,Y,Z,Lap,Lap_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,betax,betax_rhs,ASS,Symmetry,eps)
  call kodis(ex,X,Y,Z,betay,betay_rhs,SAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,betaz,betaz_rhs,SSA,Symmetry,eps)
#if (GAUGE == 0 || GAUGE == 2 || GAUGE == 3 || GAUGE == 6 || GAUGE == 7)
  call kodis(ex,X,Y,Z,dtSfx,dtSfx_rhs,ASS,Symmetry,eps)
  call kodis(ex,X,Y,Z,dtSfy,dtSfy_rhs,SAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,dtSfz,dtSfz_rhs,SSA,Symmetry,eps)
#endif
#endif

  endif

  if(co == 0)then
! ham_Res = trR + 2/3 * K^2 - A_ij * A^ij - 16 * PI * rho
! here trR is respect to physical metric
  ham_Res =   gupxx * Rxx + gupyy * Ryy + gupzz * Rzz + &
        TWO* ( gupxy * Rxy + gupxz * Rxz + gupyz * Ryz )

  ham_Res = chin1*ham_Res + F2o3 * trK * trK -(&
       gupxx * ( &
       gupxx * Axx * Axx + gupyy * Axy * Axy + gupzz * Axz * Axz + &
       TWO * (gupxy * Axx * Axy + gupxz * Axx * Axz + gupyz * Axy * Axz) ) + &
       gupyy * ( &
       gupxx * Axy * Axy + gupyy * Ayy * Ayy + gupzz * Ayz * Ayz + &
       TWO * (gupxy * Axy * Ayy + gupxz * Axy * Ayz + gupyz * Ayy * Ayz) ) + &
       gupzz * ( &
       gupxx * Axz * Axz + gupyy * Ayz * Ayz + gupzz * Azz * Azz + &
       TWO * (gupxy * Axz * Ayz + gupxz * Axz * Azz + gupyz * Ayz * Azz) ) + &
       TWO * ( &
       gupxy * ( &
       gupxx * Axx * Axy + gupyy * Axy * Ayy + gupzz * Axz * Ayz + &
       gupxy * (Axx * Ayy + Axy * Axy) + &
       gupxz * (Axx * Ayz + Axz * Axy) + &
       gupyz * (Axy * Ayz + Axz * Ayy) ) + &
       gupxz * ( &
       gupxx * Axx * Axz + gupyy * Axy * Ayz + gupzz * Axz * Azz + &
       gupxy * (Axx * Ayz + Axy * Axz) + &
       gupxz * (Axx * Azz + Axz * Axz) + &
       gupyz * (Axy * Azz + Axz * Ayz) ) + &
       gupyz * ( &
       gupxx * Axy * Axz + gupyy * Ayy * Ayz + gupzz * Ayz * Azz + &
       gupxy * (Axy * Ayz + Ayy * Axz) + &
       gupxz * (Axy * Azz + Ayz * Axz) + &
       gupyz * (Ayy * Azz + Ayz * Ayz) ) ))- F16 * PI * rho

! mov_Res_j = gupkj*(-1/chi d_k chi*A_ij + D_k A_ij) - 2/3 d_j trK - 8 PI s_j where D respect to physical metric
! store D_i A_jk - 1/chi d_i chi*A_jk in gjk_i
  call fderivs(ex,Axx,gxxx,gxxy,gxxz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)
  call fderivs(ex,Axy,gxyx,gxyy,gxyz,X,Y,Z,ANTI,ANTI,SYM ,Symmetry,0)
  call fderivs(ex,Axz,gxzx,gxzy,gxzz,X,Y,Z,ANTI,SYM ,ANTI,Symmetry,0)
  call fderivs(ex,Ayy,gyyx,gyyy,gyyz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)
  call fderivs(ex,Ayz,gyzx,gyzy,gyzz,X,Y,Z,SYM ,ANTI,ANTI,Symmetry,0)
  call fderivs(ex,Azz,gzzx,gzzy,gzzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)

  gxxx = gxxx - (  Gamxxx * Axx + Gamyxx * Axy + Gamzxx * Axz &
                 + Gamxxx * Axx + Gamyxx * Axy + Gamzxx * Axz) - chix*Axx/chin1
  gxyx = gxyx - (  Gamxxy * Axx + Gamyxy * Axy + Gamzxy * Axz &
                 + Gamxxx * Axy + Gamyxx * Ayy + Gamzxx * Ayz) - chix*Axy/chin1
  gxzx = gxzx - (  Gamxxz * Axx + Gamyxz * Axy + Gamzxz * Axz &
                 + Gamxxx * Axz + Gamyxx * Ayz + Gamzxx * Azz) - chix*Axz/chin1
  gyyx = gyyx - (  Gamxxy * Axy + Gamyxy * Ayy + Gamzxy * Ayz &
                 + Gamxxy * Axy + Gamyxy * Ayy + Gamzxy * Ayz) - chix*Ayy/chin1
  gyzx = gyzx - (  Gamxxz * Axy + Gamyxz * Ayy + Gamzxz * Ayz &
                 + Gamxxy * Axz + Gamyxy * Ayz + Gamzxy * Azz) - chix*Ayz/chin1
  gzzx = gzzx - (  Gamxxz * Axz + Gamyxz * Ayz + Gamzxz * Azz &
                 + Gamxxz * Axz + Gamyxz * Ayz + Gamzxz * Azz) - chix*Azz/chin1
  gxxy = gxxy - (  Gamxxy * Axx + Gamyxy * Axy + Gamzxy * Axz &
                 + Gamxxy * Axx + Gamyxy * Axy + Gamzxy * Axz) - chiy*Axx/chin1
  gxyy = gxyy - (  Gamxyy * Axx + Gamyyy * Axy + Gamzyy * Axz &
                 + Gamxxy * Axy + Gamyxy * Ayy + Gamzxy * Ayz) - chiy*Axy/chin1
  gxzy = gxzy - (  Gamxyz * Axx + Gamyyz * Axy + Gamzyz * Axz &
                 + Gamxxy * Axz + Gamyxy * Ayz + Gamzxy * Azz) - chiy*Axz/chin1
  gyyy = gyyy - (  Gamxyy * Axy + Gamyyy * Ayy + Gamzyy * Ayz &
                 + Gamxyy * Axy + Gamyyy * Ayy + Gamzyy * Ayz) - chiy*Ayy/chin1
  gyzy = gyzy - (  Gamxyz * Axy + Gamyyz * Ayy + Gamzyz * Ayz &
                 + Gamxyy * Axz + Gamyyy * Ayz + Gamzyy * Azz) - chiy*Ayz/chin1
  gzzy = gzzy - (  Gamxyz * Axz + Gamyyz * Ayz + Gamzyz * Azz &
                 + Gamxyz * Axz + Gamyyz * Ayz + Gamzyz * Azz) - chiy*Azz/chin1
  gxxz = gxxz - (  Gamxxz * Axx + Gamyxz * Axy + Gamzxz * Axz &
                 + Gamxxz * Axx + Gamyxz * Axy + Gamzxz * Axz) - chiz*Axx/chin1
  gxyz = gxyz - (  Gamxyz * Axx + Gamyyz * Axy + Gamzyz * Axz &
                 + Gamxxz * Axy + Gamyxz * Ayy + Gamzxz * Ayz) - chiz*Axy/chin1
  gxzz = gxzz - (  Gamxzz * Axx + Gamyzz * Axy + Gamzzz * Axz &
                 + Gamxxz * Axz + Gamyxz * Ayz + Gamzxz * Azz) - chiz*Axz/chin1
  gyyz = gyyz - (  Gamxyz * Axy + Gamyyz * Ayy + Gamzyz * Ayz &
                 + Gamxyz * Axy + Gamyyz * Ayy + Gamzyz * Ayz) - chiz*Ayy/chin1
  gyzz = gyzz - (  Gamxzz * Axy + Gamyzz * Ayy + Gamzzz * Ayz &
                 + Gamxyz * Axz + Gamyyz * Ayz + Gamzyz * Azz) - chiz*Ayz/chin1
  gzzz = gzzz - (  Gamxzz * Axz + Gamyzz * Ayz + Gamzzz * Azz &
                 + Gamxzz * Axz + Gamyzz * Ayz + Gamzzz * Azz) - chiz*Azz/chin1
movx_Res = gupxx*gxxx + gupyy*gxyy + gupzz*gxzz &
          +gupxy*gxyx + gupxz*gxzx + gupyz*gxzy &
          +gupxy*gxxy + gupxz*gxxz + gupyz*gxyz
movy_Res = gupxx*gxyx + gupyy*gyyy + gupzz*gyzz &
          +gupxy*gyyx + gupxz*gyzx + gupyz*gyzy &
          +gupxy*gxyy + gupxz*gxyz + gupyz*gyyz
movz_Res = gupxx*gxzx + gupyy*gyzy + gupzz*gzzz &
          +gupxy*gyzx + gupxz*gzzx + gupyz*gzzy &
          +gupxy*gxzy + gupxz*gxzz + gupyz*gyzz

movx_Res = movx_Res - F2o3*Kx - F8*PI*sx
movy_Res = movy_Res - F2o3*Ky - F8*PI*sy
movz_Res = movz_Res - F2o3*Kz - F8*PI*sz
  endif

#if (ABV == 1)
  call ricci_gamma(ex, X, Y, Z,                                      &
               chi,                                                  &
               dxx    ,   gxy    ,   gxz    ,   dyy    ,   gyz    ,   dzz,&
               Gamx   ,  Gamy    ,  Gamz    , &
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,&
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,&
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,&
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,&
               Symmetry)
  call constraint_bssn(ex, X, Y, Z,&
               chi,trK, &
               dxx,gxy,gxz,dyy,gyz,dzz, &
               Axx,Axy,Axz,Ayy,Ayz,Azz, &
               Gamx,Gamy,Gamz,&
               Lap,betax,betay,betaz,rho,Sx,Sy,Sz,&
               Gamxxx, Gamxxy, Gamxxz,Gamxyy, Gamxyz, Gamxzz, &
               Gamyxx, Gamyxy, Gamyxz,Gamyyy, Gamyyz, Gamyzz, &
               Gamzxx, Gamzxy, Gamzxz,Gamzyy, Gamzyz, Gamzzz, &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz, &
               ham_Res,movx_Res,movy_Res,movz_Res,Gmx_Res,Gmy_Res,Gmz_Res, &
               Symmetry)
#endif 
#if 0
#define i 2
if(Lev == 1)then
write(*,*) X(i),Y(i),Z(i)
write(*,*) Axx(i,i,i),Axy(i,i,i),Axz(i,i,i),Ayy(i,i,i),Ayz(i,i,i),Azz(i,i,i)
write(*,*) 1+Lap(i,i,i),dtSfx(i,i,i),dtSfy(i,i,i),dtSfz(i,i,i)
write(*,*) betax(i,i,i),betay(i,i,i),betaz(i,i,i)
write(*,*) 1+chi(i,i,i),Gamx(i,i,i),Gamy(i,i,i),Gamz(i,i,i)
write(*,*) gxx(i,i,i),gxy(i,i,i),gxz(i,i,i),gyy(i,i,i),gyz(i,i,i),gzz(i,i,i)
write(*,*) trK(i,i,i)
write(*,*) "====="
write(*,*) Axx_rhs(i,i,i),Axy_rhs(i,i,i),Axz_rhs(i,i,i),Ayy_rhs(i,i,i),Ayz_rhs(i,i,i),Azz_rhs(i,i,i)
write(*,*) Lap_rhs(i,i,i),dtSfx_rhs(i,i,i),dtSfy_rhs(i,i,i),dtSfz_rhs(i,i,i)
write(*,*) betax_rhs(i,i,i),betay_rhs(i,i,i),betaz_rhs(i,i,i)
write(*,*) chi_rhs(i,i,i),Gamx_rhs(i,i,i),Gamy_rhs(i,i,i),Gamz_rhs(i,i,i)
write(*,*) gxx_rhs(i,i,i),gxy_rhs(i,i,i),gxz_rhs(i,i,i),gyy_rhs(i,i,i),gyz_rhs(i,i,i),gzz_rhs(i,i,i)
write(*,*) trK_rhs(i,i,i)
endif
#undef i
!!stop
#endif

  gont = 0

  return

  end function compute_rhs_bssn_legacy



end function compute_rhs_bssn