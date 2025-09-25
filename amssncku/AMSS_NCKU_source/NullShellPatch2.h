
#ifndef NULLSHELLPATCH2_H
#define NULLSHELLPATCH2_H

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
#error NullShellPatch2 only supports 3 dimensional stuff yet
#endif

// x   x   x   x   x   o   *
//             *   o   x   x   x   x   x
// each side contribute an overlap points
// so we need half of that
#define overghost ((ghost_width + 1) / 2 + ghost_width)

class NullShellPatch2
{

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

public:
   struct pointstru
   {
      double gpox[dim]; // global cordinate
      double lpox[dim]; // local cordinate
      Block *Bg;
      int ssst; //-1: cardisian, others as sst of ss_patch source sst
      int tsst; //-1: cardisian, others as sst of ss_patch target sst
      double *coef;
      int *sind; // index position, considered dummy dimension already
      int dumyd; // the dimension which has common lines, only useful in interdata_packer
      double Jacob[2][2];
      int indz; // index position of r direction
   };

   var *gx, *gy, *gz;
   // surface variable
   var *g00, *g01, *p02, *p03, *g02, *g03;
   var *Theta22, *Theta23, *Theta33;

   // evolution variables
   var *g22o, *g23o, *g33o;
   var *g220, *g230, *g330;
   var *g22, *g23, *g33;
   var *g221, *g231, *g331;
   var *g22_rhs, *g23_rhs, *g33_rhs;

   var *RNews, *INews;
   var *omega, *dtomega;

   MyList<var> *StateList, *SynchList_pre, *SynchList_cor, *RHSList;
   MyList<var> *OldStateList, *DumpList, *CheckList;
   MyList<var> *NewsList;

   MyList<var> *g01List, *pg0AList, *g00List, *ThetaList;

   double **g01wt, **pg0Awt, **g00wt, **Thetawt;

   int myrank;
   int shape[dim]; // for (rho, sigma, X), for rho and sigma means number of points for every pi/2
   double Rmin, xmin, xmax;
   int Symmetry;
   int ingfs, fngfs;

   MyList<ss_patch> *PatL;

   MyList<pointstru> **ss_src, **ss_dst;
   MyList<pointstru> **cs_src, **cs_dst;

public:
   NullShellPatch2(int *shapei, double Rmini, double xmini, double xmaxi, int Symmetry, int myranki);

   ~NullShellPatch2();

   double getdX(int dir);
   void shellname(char *sn, int i);
   void destroypsuList(MyList<pointstru> *ct);
   MyList<Block> *compose_sh(int cpusize);
   void Dump_xyz(char *tag, double time, double dT);
   void Dump_Data(MyList<var> *DumpListi, char *tag, double time, double dT);
   void setupintintstuff(int cpusize, MyList<Patch> *CPatL, int Symmetry);
   void getlocalpox_ss(int isst, double ix, double iy, double iz, int &sst, double &lx, double &ly, double &lz);
   void getlocalpox_fake(double x, double y, double z, int &sst, double &lx, double &ly, double &lz);
   void getlocalpox(double x, double y, double z, int &sst, double &lx, double &ly, double &lz);
   void getlocalpoxsst_ss(int isst, double ix, double iy, double iz, int lsst, double &lx, double &ly, double &lz);
   void getlocalpoxsst(double x, double y, double z, int sst, double &lx, double &ly, double &lz);
   void getglobalpox(double &x, double &y, double &z, int sst, double lx, double ly, double lz);
   int getdumydimension(int acsst, int posst);
   void get_Jacob(double *pox, int tsst, int ssst, double J[2][2]);
   void prolongpointstru(MyList<pointstru> *&psul, MyList<ss_patch> *sPpi, double DH[dim],
                         MyList<Patch> *Ppi, double CDH[dim], MyList<pointstru> *pss);
   bool prolongpointstru(MyList<pointstru> *&psul, bool ssyn, int tsst, MyList<ss_patch> *sPp, double DH[dim],
                         MyList<Patch> *Pp, double CDH[dim], double x, double y, double z, int Symmetry, int rank_in, const int iz);
   bool prolongpointstru_ss(MyList<pointstru> *&psul, int tsst, MyList<ss_patch> *sPp, double DH[dim],
                            MyList<Patch> *Pp, double CDH[dim], double x, double y, double z, int Symmetry, int rank_in, const int iz);
   void Setup_Initial_Data(bool checkrun, double PhysTime);
   void Step(double dT, double PhysTime, monitor *ErrorMonitor);
   void HyperSlice(double dT, double PhysTime, monitor *ErrorMonitor, int RK_count);
   void Synch(MyList<var> *VarList, int Symmetry, double **Varwt, const short int svt);
   void fill_symmetric_boundarybuffer(MyList<var> *VarList, double **Varwt);
   void intertransfer(MyList<pointstru> **src, MyList<pointstru> **dst,
                      MyList<var> *VarList1 /* source */, MyList<var> *VarList2 /*target */,
                      int Symmetry, double **Varwt, const short int svt);
   int interdata_packer(double *data, MyList<pointstru> *src, MyList<pointstru> *dst, int rank_in, int dir,
                        MyList<var> *VarLists /* source */, MyList<var> *VarListd /* target */, int Symmetry, double **Varwt,
                        const short int svt);
   int interdata_packer_pre(double *data, MyList<pointstru> *src, MyList<pointstru> *dst, int rank_in, int dir,
                            MyList<var> *VarLists /* source */, MyList<var> *VarListd /* target */, int Symmetry, double **Varwt,
                            const short int svt);
   int interdata_packer_pot(double *data, MyList<pointstru> *src, MyList<pointstru> *dst, int rank_in, int dir,
                            MyList<var> *VarLists /* source */, MyList<var> *VarListd /* target */, int Symmetry, double **Varwt,
                            const short int svt);
   void check_pointstrul(MyList<pointstru> *pp, bool first_only);
   void checkBlock(int sst);
   void Null_Boundary(double PhysTime);
   void Compute_News(double PhysTime);
   void Interp_Points_2D(MyList<var> *VarList,
                         int NN, double **XX, /*input fake global Cartesian coordinate*/
                         double *Shellf, int Symmetry);
   double Error_Check(double PhysTime);
};

#endif /* NULLSHELLPATCH2_H */
