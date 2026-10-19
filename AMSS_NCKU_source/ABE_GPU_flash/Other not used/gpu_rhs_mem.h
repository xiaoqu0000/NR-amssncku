#ifndef GPU_MEM_H_
#define GPU_MEM_H_
#include "macrodef.fh"
struct Meta
{
	//---------------in/out-------------------
	// int * ex;
	// int* Symmetry,Lev,co; //not array	//in
	// double *  T;				//not array	//in
	double *X, *Y, *Z;									 // in
	double *chi, *dxx, *dyy, *dzz;						 // inout
	double *trK;										 // in
	double *gxy, *gxz, *gyz;							 // in
	double *Axx, *Axy, *Axz, *Ayy, *Ayz, *Azz;			 // in
	double *Gamx, *Gamy, *Gamz;							 // in
	double *Lap, *betax, *betay, *betaz;				 // inout
	double *dtSfx, *dtSfy, *dtSfz;						 // in
	double *chi_rhs, *trK_rhs;							 // out
	double *gxx_rhs, *gxy_rhs, *gxz_rhs;				 // out
	double *gyy_rhs, *gyz_rhs, *gzz_rhs;				 // out
	double *Axx_rhs, *Axy_rhs, *Axz_rhs;				 // out
	double *Ayy_rhs, *Ayz_rhs, *Azz_rhs;				 // out
	double *Gamx_rhs, *Gamy_rhs, *Gamz_rhs;				 // out
	double *Lap_rhs, *betax_rhs, *betay_rhs, *betaz_rhs; // out
	double *dtSfx_rhs, *dtSfy_rhs, *dtSfz_rhs;			 // out
	double *rho, *Sx, *Sy, *Sz;							 // in
	double *Sxx, *Sxy, *Sxz, *Syy, *Syz, *Szz;			 // in

	// when out, physical second kind of connection  //out
	double *Gamxxx, *Gamxxy, *Gamxxz;
	double *Gamxyy, *Gamxyz, *Gamxzz;
	double *Gamyxx, *Gamyxy, *Gamyxz;
	double *Gamyyy, *Gamyyz, *Gamyzz;
	double *Gamzxx, *Gamzxy, *Gamzxz;
	double *Gamzyy, *Gamzyz, *Gamzzz;

	// when out, physical Ricci tensor
	double *Rxx, *Rxy, *Rxz, *Ryy, *Ryz, *Rzz; // out
	// double * eps;						//in
	double *ham_Res, *movx_Res, *movy_Res, *movz_Res; // inout
	double *Gmx_Res, *Gmy_Res, *Gmz_Res;			  // inout

	//---------------local-------------------

	double *gxx, *gyy, *gzz, *chix, *chiy, *chiz, *gxxx, *gxyx, *gxzx, *gyyx, *gyzx, *gzzx, *gxxy, *gxyy, *gxzy, *gyyy, *gyzy, *gzzy, *gxxz, *gxyz, *gxzz, *gyyz, *gyzz, *gzzz, *Lapx, *Lapy, *Lapz, *betaxx, *betaxy, *betaxz, *betayx, *betayy, *betayz, *betazx, *betazy, *betazz, *Gamxx, *Gamxy, *Gamxz, *Gamyx, *Gamyy, *Gamyz, *Gamzx, *Gamzy, *Gamzz, *Kx, *Ky, *Kz, *div_beta, *S, *f, *fxx, *fxy, *fxz, *fyy, *fyz, *fzz, *Gamxa, *Gamya, *Gamza, *alpn1, *chin1, *gupxx, *gupxy, *gupxz, *gupyy, *gupyz, *gupzz;

	//---------------subroutine----------------
	double *fh;
	double *fh2;

	/*double *SSS;
	double *AAS;
	double *ASA;
	double *SAA;
	double *ASS;
	double *SAS;
	double *SSA;*/
//---------------GAUGE--------------
#if (GAUGE == 2 || GAUGE == 3 || GAUGE == 4 || GAUGE == 5 || GAUGE == 6 || GAUGE == 7)
	double *reta;
#endif
};

// //------init constant memory---------

// // 1-----for compute_rhs_bssn---------
// extern __constant__ Meta metac;
// extern __constant__ int ex_c[3];
// extern __constant__ double T_c;
// extern __constant__ int Symmetry_c;
// extern __constant__ int Lev_c;
// extern __constant__ int co_c;
// extern __constant__ double eps_c;
// // local
// extern __constant__ double dX; // dX,dY,dZ
// extern __constant__ double dY;
// extern __constant__ double dZ;
// extern __constant__ double ZEO;
// extern __constant__ double ONE;
// extern __constant__ double TWO;
// extern __constant__ double FOUR;
// extern __constant__ double EIGHT;
// extern __constant__ double HALF;
// extern __constant__ double THR;
// extern __constant__ double SYM;
// extern __constant__ double ANTI ;
// extern __constant__ double FF ;
// extern __constant__ double eta;
// extern __constant__ double F1o3;
// extern __constant__ double F2o3;
// extern __constant__ double F3o2;
// extern __constant__ double F1o6;
// extern __constant__ double F8;
// extern __constant__ double F16 ;
// extern __constant__ double PI;
// /*__constant__ double SSS[3] = {1,1,1};
// __constant__ double AAS[3] = {-1,-1,1};
// __constant__ double ASA[3] = {-1,1,-1};
// __constant__ double SAA[3] = {1,-1,-1};
// __constant__ double ASS[3] = {-1,1,1};
// __constant__ double SAS[3] = {1,-1,1};
// __constant__ double SSA[3] = {1,1,-1};*/

// // 2--------for fderivs------------
// extern __constant__ int ijk_min[3];
// extern __constant__ int ijk_min2[3];
// extern __constant__ int ijk_min3[3];
// extern __constant__ int ijk_max[3];
// extern __constant__ double d12dxyz[3];
// extern __constant__ double d2dxyz[3];

// // 3--------for fdderivs------------
// extern __constant__ double Sdxdx;
// extern __constant__ double Sdydy;
// extern __constant__ double Sdzdz;
// extern __constant__ double Fdxdx;
// extern __constant__ double Fdydy;
// extern __constant__ double Fdzdz;
// extern __constant__ double Sdxdy;
// extern __constant__ double Sdxdz;
// extern __constant__ double Sdydz;
// extern __constant__ double Fdxdy;
// extern __constant__ double Fdxdz;
// extern __constant__ double Fdydz;

// // my own
// extern __constant__ int STEP_SIZE;
// /*__constant__ int MATRIX_SIZE;
// __constant__ int MATRIX_SIZE_FH;
// __constant__ int SQUARE_SIZE;
// __constant__ int SQUARE_SIZE_FH;
// __constant__ int LINE_SIZE_FH;*/

// extern __constant__ int _1D_SIZE[4]; // start from 0 !!
// extern __constant__ int _2D_SIZE[4]; ////start from 0 !!
// extern __constant__ int _3D_SIZE[4]; ////start from 0 !!

// #if (GAUGE == 6 || GAUGE == 7)
// __constant__ int BHN;
// __constant__ double Porg[9];
// __constant__ double Mass[3];
// __constant__ double /*r1,r2*/, M, A, /*w1,w2 (== 12)*/, C1, C2;
// #endif

/**/
#endif
