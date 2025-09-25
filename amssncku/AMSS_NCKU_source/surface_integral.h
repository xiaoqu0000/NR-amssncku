//$Id: surface_integral.h,v 1.9 2013/08/20 11:49:05 zjcao Exp $
#ifndef SURFACE_INTEGRAL_H
#define SURFACE_INTEGRAL_H

#ifdef newc
#include <iostream>
#include <iomanip>
#include <fstream>
#include <strstream>
#include <cmath>
using namespace std;
#else
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <strstream>
#include <math.h>
#endif

#include "cgh.h"
#include "ShellPatch.h"
#include "NullShellPatch.h"
#include "NullShellPatch2.h"
#include "var.h"
#include "monitor.h"

class surface_integral
{

private:
	int Symmetry, factor;
	int N_theta, N_phi; // Number of points in Theta & Phi directions
	double dphi, dcostheta;
	double *arcostheta, *wtcostheta;
	int n_tot; // size of arrays

	double *nx_g, *ny_g, *nz_g; // global list of unit normals
	int myrank, cpusize;

public:
	surface_integral(int iSymmetry);
	~surface_integral();

	void surf_Wave(double rex, int lev, cgh *GH, var *Rpsi4, var *Ipsi4,
				   int spinw, int maxl, int NN, double *RP, double *IP,
				   monitor *Monitor); // NN is the length of RP and IP
									  // this routine can only deal with the symmetry of Psi4
	void surf_Wave(double rex, int lev, ShellPatch *GH, var *Rpsi4, var *Ipsi4,
				   int spinw, int maxl, int NN, double *RP, double *IP,
				   monitor *Monitor);
	void surf_Wave(double rex, int lev, NullShellPatch *GH, var *Rpsi4, var *Ipsi4,
				   int spinw, int maxl, int NN, double *RP, double *IP,
				   monitor *Monitor);
	void surf_Wave(double rex, int lev, NullShellPatch2 *GH, var *Rpsi4, var *Ipsi4,
				   int spinw, int maxl, int NN, double *RP, double *IP,
				   monitor *Monitor);
	void surf_Wave(double rex, int lev, ShellPatch *GH,
				   var *Ex, var *Ey, var *Ez, var *Bx, var *By, var *Bz,
				   var *chi, var *gxx, var *gxy, var *gxz, var *gyy, var *gyz, var *gzz,
				   int spinw, int maxl, int NN, double *RP, double *IP,
				   monitor *Monitor); // NN is the length of RP and IP
	void surf_Wave(double rex, int lev, cgh *GH,
				   var *Ex, var *Ey, var *Ez, var *Bx, var *By, var *Bz,
				   var *chi, var *gxx, var *gxy, var *gxz, var *gyy, var *gyz, var *gzz,
				   int spinw, int maxl, int NN, double *RP, double *IP,
				   monitor *Monitor,
				   void (*funcs)(double &, double &, double &,
								 double &, double &, double &, double &, double &, double &, double &,
								 double &, double &, double &, double &, double &, double &,
								 double &, double &)); // NN is the length of RP and IP
	void surf_Wave(double rex, int lev, ShellPatch *GH,
				   var *Ex, var *Ey, var *Ez, var *Bx, var *By, var *Bz,
				   var *chi, var *gxx, var *gxy, var *gxz, var *gyy, var *gyz, var *gzz,
				   int spinw, int maxl, int NN, double *RP, double *IP,
				   monitor *Monitor,
				   void (*funcs)(double &, double &, double &,
								 double &, double &, double &, double &, double &, double &, double &,
								 double &, double &, double &, double &, double &, double &,
								 double &, double &)); // NN is the length of RP and IP
	void surf_MassPAng(double rex, int lev, cgh *GH, var *chi, var *trK,
					   var *gxx, var *gxy, var *gxz, var *gyy, var *gyz, var *gzz,
					   var *Axx, var *Axy, var *Axz, var *Ayy, var *Ayz, var *Azz,
					   var *Gmx, var *Gmy, var *Gmz,
					   var *Sfx_rhs, var *Sfy_rhs, var *Sfz_rhs,
					   double *Rout, monitor *Monitor);
	void surf_MassPAng(double rex, int lev, ShellPatch *GH, var *chi, var *trK,
					   var *gxx, var *gxy, var *gxz, var *gyy, var *gyz, var *gzz,
					   var *Axx, var *Axy, var *Axz, var *Ayy, var *Ayz, var *Azz,
					   var *Gmx, var *Gmy, var *Gmz,
					   var *Sfx_rhs, var *Sfy_rhs, var *Sfz_rhs,
					   double *Rout, monitor *Monitor);
	void surf_Wave(double rex, cgh *GH, ShellPatch *SH,
				   var *chi, var *trK,
				   var *gxx, var *gxy, var *gxz, var *gyy, var *gyz, var *gzz,
				   var *Axx, var *Axy, var *Axz, var *Ayy, var *Ayz, var *Azz,
				   var *chix, var *chiy, var *chiz,
				   var *trKx, var *trKy, var *trKz,
				   var *Axxx, var *Axxy, var *Axxz,
				   var *Axyx, var *Axyy, var *Axyz,
				   var *Axzx, var *Axzy, var *Axzz,
				   var *Ayyx, var *Ayyy, var *Ayyz,
				   var *Ayzx, var *Ayzy, var *Ayzz,
				   var *Azzx, var *Azzy, var *Azzz,
				   var *Gamxxx, var *Gamxxy, var *Gamxxz, var *Gamxyy, var *Gamxyz, var *Gamxzz,
				   var *Gamyxx, var *Gamyxy, var *Gamyxz, var *Gamyyy, var *Gamyyz, var *Gamyzz,
				   var *Gamzxx, var *Gamzxy, var *Gamzxz, var *Gamzyy, var *Gamzyz, var *Gamzzz,
				   var *Rxx, var *Rxy, var *Rxz, var *Ryy, var *Ryz, var *Rzz,
				   int spinw, int maxl, int NN, double *RP, double *IP,
				   monitor *Monitor);
	bool SR_Interp_Points(MyList<var> *VarList, cgh *GH, ShellPatch *SH,
						  int NN, double **XX, double *Shellf);

	void surf_MassPAng(double rex, int lev, cgh *GH, var *chi, var *trK,
					   var *gxx, var *gxy, var *gxz, var *gyy, var *gyz, var *gzz,
					   var *Axx, var *Axy, var *Axz, var *Ayy, var *Ayz, var *Azz,
					   var *Gmx, var *Gmy, var *Gmz,
					   var *Sfx_rhs, var *Sfy_rhs, var *Sfz_rhs, // temparay memory for mass^i
					   double *Rout, monitor *Monitor, MPI_Comm Comm_here);
	void surf_Wave(double rex, int lev, cgh *GH, var *Rpsi4, var *Ipsi4,
				   int spinw, int maxl, int NN, double *RP, double *IP,
				   monitor *Monitor, MPI_Comm Comm_here);
};
#endif /* SURFACE_INTEGRAL_H */
