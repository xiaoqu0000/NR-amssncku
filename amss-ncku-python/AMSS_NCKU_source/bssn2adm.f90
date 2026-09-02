
!-------------------------------------------------------------------------------!
! convert bssn variables to ADM variables                                       !
!-------------------------------------------------------------------------------!
  subroutine bssn2adm(ex,chi,trK, &
               gxx,gxy,gxz,gyy,gyz,gzz, &
               Axx,Axy,Axz,Ayy,Ayz,Azz, &
               adm_gxx,adm_gxy,adm_gxz,adm_gyy,adm_gyz,adm_gzz, &
               Kxx,Kxy,Kxz,Kyy,Kyz,Kzz)

  implicit none
!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3)
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::chi,trK
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::gxx,gxy,gxz,gyy,gyz,gzz
  double precision,intent(in),dimension(ex(1),ex(2),ex(3))::Axx,Axy,Axz,Ayy,Ayz,Azz

  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: adm_gxx,adm_gxy,adm_gxz,adm_gyy,adm_gyz,adm_gzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Kxx,Kxy,Kxz,Kyy,Kyz,Kzz

  real*8, parameter :: F1o3=1.d0/3.d0

  adm_gxx = gxx/chi
  adm_gxy = gxy/chi
  adm_gxz = gxz/chi
  adm_gyy = gyy/chi
  adm_gyz = gyz/chi
  adm_gzz = gzz/chi

  Kxx = Axx/chi+F1o3*trK*adm_gxx
  Kxy = Axy/chi+F1o3*trK*adm_gxy
  Kxz = Axz/chi+F1o3*trK*adm_gxz
  Kyy = Ayy/chi+F1o3*trK*adm_gyy
  Kyz = Ayz/chi+F1o3*trK*adm_gyz
  Kzz = Azz/chi+F1o3*trK*adm_gzz

  return

  end subroutine bssn2adm
