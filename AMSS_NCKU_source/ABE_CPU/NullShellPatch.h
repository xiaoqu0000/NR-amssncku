
#ifndef NULLSHELLPATCH_H
#define NULLSHELLPATCH_H

#ifdef newc
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <complex>
using namespace std;
#else
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#endif

#include <mpi.h>
#include "MyList.h"
#include "Block.h"
#include "Parallel.h"
#include "ShellPatch.h"
#include "var.h"
#include "macrodef.h" //need dim here; Vertex or Cell; ghost_width

#if (dim != 3)
#error NullShellPatch only supports 3 dimensional stuff yet
#endif

class xp_npatch : public ss_patch
{
public:
   xp_npatch(int ingfsi, int fngfsi, int *shapei, double *bboxi, int myranki) : ss_patch(ingfsi, fngfsi, shapei, bboxi, myranki) { sst = 2; };
};

class xm_npatch : public ss_patch
{
public:
   xm_npatch(int ingfsi, int fngfsi, int *shapei, double *bboxi, int myranki) : ss_patch(ingfsi, fngfsi, shapei, bboxi, myranki) { sst = 3; };
};
class yp_npatch : public ss_patch
{
public:
   yp_npatch(int ingfsi, int fngfsi, int *shapei, double *bboxi, int myranki) : ss_patch(ingfsi, fngfsi, shapei, bboxi, myranki) { sst = 4; };
};

class ym_npatch : public ss_patch
{
public:
   ym_npatch(int ingfsi, int fngfsi, int *shapei, double *bboxi, int myranki) : ss_patch(ingfsi, fngfsi, shapei, bboxi, myranki) { sst = 5; };
};
class zp_npatch : public ss_patch
{
public:
   zp_npatch(int ingfsi, int fngfsi, int *shapei, double *bboxi, int myranki) : ss_patch(ingfsi, fngfsi, shapei, bboxi, myranki) { sst = 0; };
};

class zm_npatch : public ss_patch
{
public:
   zm_npatch(int ingfsi, int fngfsi, int *shapei, double *bboxi, int myranki) : ss_patch(ingfsi, fngfsi, shapei, bboxi, myranki) { sst = 1; };
};

class NullShellPatch
{

public:
   struct pointstru
   {
      double gpox[dim]; // global cordinate
      double lpox[dim]; // local cordinate
      Block *Bg;
      int ssst; //-1: cardisian, others as sst of ss_patch source sst
      int tsst; //-1: cardisian, others as sst of ss_patch target sst
      double *coef;
      int *sind;
      int dumyd;            // the dimension which has common lines, only useful in interdata_packer
      complex<double> swtf; // exp(i gamma) of Eq.(26) of CQG 24 S327
   };

   var *FXZEO;
   var *gx, *gy, *gz;
   // we always assume the number of VarList = 2* the number of Varwt
   // so VarList must apear with pairs, either components of complex number or a fake pair
   var *beta, *W;
   var *Rnu, *Inu, *Rk, *Ik, *RB, *IB;
   var *RQ, *IQ, *RU, *IU, *RTheta, *ITheta;
   var *KK, *HKK, *KKx, *HKKx;
   var *RJo, *IJo, *omegao;
   var *RJ0, *IJ0, *omega0;
   var *RJ, *IJ, *omega;
   var *RJ1, *IJ1, *omega1;
   var *RJ_rhs, *IJ_rhs, *omega_rhs;

   var *quR1, *quR2, *quI1, *quI2;
   var *qlR1, *qlR2, *qlI1, *qlI2;
   var *gR, *gI;
   var *dquR1, *dquR2, *dquI1, *dquI2;
   var *bdquR1, *bdquR2, *bdquI1, *bdquI2;
   var *dgR, *dgI;
   var *bdgR, *bdgI;

   var *RNews, *INews;

   MyList<var> *StateList, *SynchList_pre, *SynchList_cor, *RHSList;
   MyList<var> *OldStateList, *DumpList, *CheckList;

   MyList<var> *betaList, *QUList, *WTheList, *TheList, *JrhsList, *J1List;
   int betawt[1], QUwt[2], WThewt[2];

   int myrank;
   int shape[dim]; // for (rho, sigma, X), for rho and sigma means number of points for every pi/2
   double Rmin, xmin, xmax;
   int Symmetry;
   int ingfs, fngfs;

   MyList<ss_patch> *PatL;

   MyList<pointstru> **ss_src, **ss_dst;
   MyList<pointstru> **cs_src, **cs_dst;

public:
   NullShellPatch(int *shapei, double Rmini, double xmini, double xmaxi, int Symmetry, int myranki);

   ~NullShellPatch();

   void destroypsuList(MyList<pointstru> *ct);
   void fill_symmetric_boundarybuffer(MyList<var> *VarList, int *Varwt);
   MyList<Block> *compose_sh(int cpusize);
   int getdumydimension(int acsst, int posst);
   void Setup_dyad();
   void Setup_Initial_Data(bool checkrun, double PhysTime);
   void eth_derivs(var *Rv, var *Iv, var *ethRv, var *ethIv, int s, int e);
   void eth_dderivs(var *Rv, var *Iv, var *ethRv, var *ethIv, int s, int e1, int e2);
   void getlocalpox_ss(int isst, double ix, double iy, double iz, int &sst, double &lx, double &ly, double &lz);
   void getlocalpox_fake(double x, double y, double z, int &sst, double &lx, double &ly, double &lz);
   void getlocalpox(double x, double y, double z, int &sst, double &lx, double &ly, double &lz);
   void getlocalpoxsst_ss(int isst, double ix, double iy, double iz, int lsst, double &lx, double &ly, double &lz);
   void getlocalpoxsst(double x, double y, double z, int sst, double &lx, double &ly, double &lz);
   void getglobalpox(double &x, double &y, double &z, int sst, double lx, double ly, double lz);
   complex<double> get_swtf(double *pox, int tsst, int ssst);
   void prolongpointstru(MyList<pointstru> *&psul, MyList<ss_patch> *sPpi, double DH[dim],
                         MyList<Patch> *Ppi, double CDH[dim], MyList<pointstru> *pss);
   bool prolongpointstru(MyList<pointstru> *&psul, bool ssyn, int tsst, MyList<ss_patch> *sPp, double DH[dim],
                         MyList<Patch> *Pp, double CDH[dim], double x, double y, double z, int Symmetry, int rank_in);
   bool prolongpointstru_ss(MyList<pointstru> *&psul, int tsst, MyList<ss_patch> *sPp, double DH[dim],
                            MyList<Patch> *Pp, double CDH[dim], double x, double y, double z, int Symmetry, int rank_in);
   void setupintintstuff(int cpusize, MyList<Patch> *CPatL, int Symmetry);
   void checkPatch();
   void checkBlock(int sst);
   double getdX(int dir);
   void shellname(char *sn, int i);
   void Dump_xyz(char *tag, double time, double dT);
   void Dump_Data(MyList<var> *DumpListi, char *tag, double time, double dT);
   void intertransfer(MyList<pointstru> **src, MyList<pointstru> **dst,
                      MyList<var> *VarList1 /* source */, MyList<var> *VarList2 /*target */,
                      int Symmetry, int *Varwt);
   int interdata_packer(double *data, MyList<pointstru> *src, MyList<pointstru> *dst, int rank_in, int dir,
                        MyList<var> *VarLists /* source */, MyList<var> *VarListd /* target */, int Symmetry, int *Varwt);
   void Synch(MyList<var> *VarList, int Symmetry, int *Varwt);
   void CS_Inter(MyList<var> *VarList, int Symmetry, int *Varwt);
   void check_pointstrul(MyList<pointstru> *pp, bool first_only);
   void check_pointstrul2(MyList<pointstru> *pp, int first_last_only);
   void matchcheck(MyList<Patch> *CPatL);
   void Interp_Points(MyList<var> *VarList,
                      int NN, double **XX, /*input global Cartesian coordinate*/
                      double *Shellf, int Symmetry);
   void Interp_Points_2D(MyList<var> *VarList,
                         int NN, double **XX, /*input global Cartesian coordinate*/
                         double *Shellf, int Symmetry);
   void Step(double dT, double PhysTime, monitor *ErrorMonitor);
   void Null_Boundary(double PhysTime);
   void HyperSlice(double dT, double PhysTime, monitor *ErrorMonitor, int RK_count);
   double News_Error_Check(double PhysTime, double dT, bool dp);
   double Error_Check(double PhysTime, double dT, bool dp);
   double EqTheta_Check(double PhysTime, double dT, bool dp);
   void Compute_News(double PhysTime, double dT, bool dp);
   void Check_News(double PhysTime, double dT, bool dp);
};

#endif /* NULLSHELLPATCH_H */
