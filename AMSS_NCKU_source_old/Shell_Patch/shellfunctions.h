
#ifndef SHELLFUNCTIONS_H
#define SHELLFUNCTIONS_H

#ifdef fortran1
#define f_get_initial_nbhs_sh get_initial_nbhs_sh
#define f_xp_getxyz xp_getxyz
#define f_xm_getxyz xm_getxyz
#define f_yp_getxyz yp_getxyz
#define f_ym_getxyz ym_getxyz
#define f_zp_getxyz zp_getxyz
#define f_zm_getxyz zm_getxyz
#define f_xpm_getjacobian xpm_getjacobian
#define f_ypm_getjacobian ypm_getjacobian
#define f_zpm_getjacobian zpm_getjacobian
#define f_shellcordpar shellcordpar
#endif
#ifdef fortran2
#define f_get_initial_nbhs_sh GET_INITIAL_NBHS_SH
#define f_xp_getxyz XP_GETXYZ
#define f_xm_getxyz XM_GETXYZ
#define f_yp_getxyz YP_GETXYZ
#define f_ym_getxyz YM_GETXYZ
#define f_zp_getxyz ZP_GETXYZ
#define f_zm_getxyz ZM_GETXYZ
#define f_xpm_getjacobian XPM_GETJACOBIAN
#define f_ypm_getjacobian YPM_GETJACOBIAN
#define f_zpm_getjacobian ZPM_GETJACOBIAN
#define f_shellcordpar SHELLCORDPAR
#endif
#ifdef fortran3
#define f_get_initial_nbhs_sh get_initial_nbhs_sh_
#define f_xp_getxyz xp_getxyz_
#define f_xm_getxyz xm_getxyz_
#define f_yp_getxyz yp_getxyz_
#define f_ym_getxyz ym_getxyz_
#define f_zp_getxyz zp_getxyz_
#define f_zm_getxyz zm_getxyz_
#define f_xpm_getjacobian xpm_getjacobian_
#define f_ypm_getjacobian ypm_getjacobian_
#define f_zpm_getjacobian zpm_getjacobian_
#define f_shellcordpar shellcordpar_
#endif

extern "C"
{
	void f_get_initial_nbhs_sh(int *, double *, double *, double *,
							   double *, double *,
							   double *, double *, double *, double *, double *, double *,
							   double *, double *, double *, double *, double *, double *,
							   double *, double *, double *,
							   double *, double *, double *, double *,
							   double *, double *, double *,
							   double *, double *, double *, double *, int &);
}

extern "C"
{
	void f_xp_getxyz(int *, double *, double *, double *, double *, double *, double *);
}
extern "C"
{
	void f_xm_getxyz(int *, double *, double *, double *, double *, double *, double *);
}
extern "C"
{
	void f_yp_getxyz(int *, double *, double *, double *, double *, double *, double *);
}
extern "C"
{
	void f_ym_getxyz(int *, double *, double *, double *, double *, double *, double *);
}
extern "C"
{
	void f_zp_getxyz(int *, double *, double *, double *, double *, double *, double *);
}
extern "C"
{
	void f_zm_getxyz(int *, double *, double *, double *, double *, double *, double *);
}

extern "C"
{
	void f_xpm_getjacobian(int *, double *, double *, double *, double *, double *, double *,
						   double *, double *, double *, double *, double *, double *, double *, double *, double *,
						   double *, double *, double *, double *, double *, double *,
						   double *, double *, double *, double *, double *, double *,
						   double *, double *, double *, double *, double *, double *);
}
extern "C"
{
	void f_ypm_getjacobian(int *, double *, double *, double *, double *, double *, double *,
						   double *, double *, double *, double *, double *, double *, double *, double *, double *,
						   double *, double *, double *, double *, double *, double *,
						   double *, double *, double *, double *, double *, double *,
						   double *, double *, double *, double *, double *, double *);
}
extern "C"
{
	void f_zpm_getjacobian(int *, double *, double *, double *, double *, double *, double *,
						   double *, double *, double *, double *, double *, double *, double *, double *, double *,
						   double *, double *, double *, double *, double *, double *,
						   double *, double *, double *, double *, double *, double *,
						   double *, double *, double *, double *, double *, double *);
}

extern "C"
{
	void f_shellcordpar(double &, double &, double &, double &);
}

#endif /* SHELLFUNCTIONS_H */
