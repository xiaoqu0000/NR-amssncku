
#ifdef newc
#include <sstream>
#include <cstdio>
#include <map>
using namespace std;
#else
#include <stdio.h>
#include <map.h>
#endif

#include <time.h>

#include "macrodef.h"
#include "misc.h"
#include "Ansorg.h"
#include "fmisc.h"
#include "Parallel.h"
#include "Z4c_class.h"
#include "bssn_rhs.h"
#include "initial_puncture.h"
#include "enforce_algebra.h"
#include "rungekutta4_rout.h"
#include "sommerfeld_rout.h"
#include "getnp4.h"
#include "shellfunctions.h"
#include "cpbc.h"
#include "kodiss.h"
#include "parameters.h"

#ifdef With_AHF
#include "derivatives.h"
#include "myglobal.h"
#endif

//================================================================================================

// 定义 Z4c_class

// 它是继承了父类 bssn_class 的某些成员和方法，并对另一些成员和方法进行修改
// 修改的成员和方法在下面（以及头文件 Z4c_class.h ）中定义
// 其余的继承父类 bssn_class（头文件 bssn_class.h 中声明）

Z4c_class::Z4c_class(double Couranti, double StartTimei, double TotalTimei, 
                     double DumpTimei, double d2DumpTimei, 
                     double CheckTimei, double AnasTimei,
                     int Symmetryi, int checkruni, char *checkfilenamei, 
                     double numepssi, double numepsbi, double numepshi,
                     int a_levi, int maxli, int decni, double maxrexi, double drexi) 
                     : bssn_class(Couranti, StartTimei, TotalTimei, 
                                  DumpTimei, d2DumpTimei, CheckTimei, AnasTimei,
                                  Symmetryi, checkruni, checkfilenamei, numepssi, numepsbi, numepshi,
                                  a_levi, maxli, decni, maxrexi, drexi)
{
}

//================================================================================================



//================================================================================================

// 该成员函数用于对 Class 初始化

//================================================================================================

void Z4c_class::Initialize()
{
  TZo = new var("TZo", ngfs++, 1, 1, 1);
  TZ0 = new var("TZ0", ngfs++, 1, 1, 1);
  TZ = new var("TZ", ngfs++, 1, 1, 1);
  TZ1 = new var("TZ1", ngfs++, 1, 1, 1);
  TZ_rhs = new var("TZ_rhs", ngfs++, 1, 1, 1);

  if (myrank == 0)
    cout << "you have setted " << ngfs << " grid functions." << endl;

  OldStateList->insert(TZo);
  StateList->insert(TZ0);
  RHSList->insert(TZ_rhs);
  SynchList_pre->insert(TZ);
  SynchList_cor->insert(TZ1);
  // DumpList->insert(TZ0);
  ConstraintList->insert(TZ0);

  CheckPoint->addvariablelist(StateList);
  CheckPoint->addvariablelist(OldStateList);

  char pname[50];
  {
    map<string, string>::iterator iter = parameters::str_par.find("inputpar");
    if (iter != parameters::str_par.end())
    {
      strcpy(pname, (iter->second).c_str());
    }
    else
    {
      cout << "Error inputpar" << endl;
      exit(0);
    }
  }
  GH = new cgh(0, ngfs, Symmetry, pname, checkrun, ErrorMonitor);
  if (checkrun)
    CheckPoint->readcheck_cgh(PhysTime, GH, myrank, nprocs, Symmetry);
  else
    GH->compose_cgh(nprocs);

#ifdef WithShell
  SH = new ShellPatch(0, ngfs, pname, Symmetry, myrank, ErrorMonitor);
  if (!checkrun)
    SH->matchcheck(GH->PatL[0]);
  SH->compose_sh(nprocs);
  SH->setupcordtrans();
  SH->Dump_xyz(0, 0, 1);
  SH->setupintintstuff(nprocs, GH->PatL[0], Symmetry);

  if (checkrun)
    CheckPoint->readcheck_sh(SH, myrank);
#endif

  double h = GH->PatL[0]->data->blb->data->getdX(0);
  for (int i = 1; i < dim; i++)
    h = Mymin(h, GH->PatL[0]->data->blb->data->getdX(i));
  dT = Courant * h;

  if (checkrun)
  {
    CheckPoint->read_Black_Hole_position(BH_num_input, BH_num, Porg0, Pmom, Spin, Mass, Porgbr, Porg, Porg1, Porg_rhs);
  }
  else
  {
    PhysTime = StartTime;
    Setup_Black_Hole_position();
  }
}

//================================================================================================




//================================================================================================

// 该成员函数为析构函数，用于删除变量

//================================================================================================

Z4c_class::~Z4c_class()
{
  delete TZo;
  delete TZ0;
  delete TZ;
  delete TZ1;
  delete TZ_rhs;
}

//================================================================================================




//================================================================================================

// 该成员函数设定了时间演化过程中的单步时间演化

//================================================================================================

#define MRBD 0 // 0: fix BD for meshrefinement level; 1: sommerfeld_bam for them; 2: sommerfeld_yo for them

#ifndef CPBC
// for sommerfeld boundary

void Z4c_class::Step(int lev, int YN)
{
  double dT_lev = dT * pow(0.5, Mymax(lev, trfls));
#ifdef With_AHF
  AH_Step_Find(lev, dT_lev);
#endif
  bool BB = fgt(PhysTime, StartTime, dT_lev / 2);
  double ndeps = numepss;
  if (lev < GH->movls)
    ndeps = numepsb;
  double TRK4 = PhysTime;
  int iter_count = 0; // count RK4 substeps
  int pre = 0, cor = 1;
  int ERROR = 0;

  MyList<ss_patch> *sPp;
  // Predictor
  MyList<Patch> *Pp = GH->PatL[lev];
  while (Pp)
  {
    MyList<Block> *BP = Pp->data->blb;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
#if (AGM == 0)
        f_enforce_ga(cg->shape,
                     cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                     cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                     cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                     cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn]);
#endif

        if (f_compute_rhs_Z4c(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                              cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                              cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                              cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                              cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                              cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                              cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                              cg->fgfs[Lap0->sgfn], 
                              cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                              cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn],
                              cg->fgfs[TZ0->sgfn],
                              cg->fgfs[phi_rhs->sgfn], cg->fgfs[trK_rhs->sgfn],
                              cg->fgfs[gxx_rhs->sgfn], cg->fgfs[gxy_rhs->sgfn], cg->fgfs[gxz_rhs->sgfn],
                              cg->fgfs[gyy_rhs->sgfn], cg->fgfs[gyz_rhs->sgfn], cg->fgfs[gzz_rhs->sgfn],
                              cg->fgfs[Axx_rhs->sgfn], cg->fgfs[Axy_rhs->sgfn], cg->fgfs[Axz_rhs->sgfn],
                              cg->fgfs[Ayy_rhs->sgfn], cg->fgfs[Ayz_rhs->sgfn], cg->fgfs[Azz_rhs->sgfn],
                              cg->fgfs[Gmx_rhs->sgfn], cg->fgfs[Gmy_rhs->sgfn], cg->fgfs[Gmz_rhs->sgfn],
                              cg->fgfs[Lap_rhs->sgfn], 
                              cg->fgfs[Sfx_rhs->sgfn], cg->fgfs[Sfy_rhs->sgfn], cg->fgfs[Sfz_rhs->sgfn],
                              cg->fgfs[dtSfx_rhs->sgfn], cg->fgfs[dtSfy_rhs->sgfn], cg->fgfs[dtSfz_rhs->sgfn],
                              cg->fgfs[TZ_rhs->sgfn],
                              cg->fgfs[rho->sgfn], cg->fgfs[Sx->sgfn], cg->fgfs[Sy->sgfn], cg->fgfs[Sz->sgfn],
                              cg->fgfs[Sxx->sgfn], cg->fgfs[Sxy->sgfn], cg->fgfs[Sxz->sgfn], 
                              cg->fgfs[Syy->sgfn], cg->fgfs[Syz->sgfn], cg->fgfs[Szz->sgfn],
                              cg->fgfs[Gamxxx->sgfn], cg->fgfs[Gamxxy->sgfn], cg->fgfs[Gamxxz->sgfn],
                              cg->fgfs[Gamxyy->sgfn], cg->fgfs[Gamxyz->sgfn], cg->fgfs[Gamxzz->sgfn],
                              cg->fgfs[Gamyxx->sgfn], cg->fgfs[Gamyxy->sgfn], cg->fgfs[Gamyxz->sgfn],
                              cg->fgfs[Gamyyy->sgfn], cg->fgfs[Gamyyz->sgfn], cg->fgfs[Gamyzz->sgfn],
                              cg->fgfs[Gamzxx->sgfn], cg->fgfs[Gamzxy->sgfn], cg->fgfs[Gamzxz->sgfn],
                              cg->fgfs[Gamzyy->sgfn], cg->fgfs[Gamzyz->sgfn], cg->fgfs[Gamzzz->sgfn],
                              cg->fgfs[Rxx->sgfn], cg->fgfs[Rxy->sgfn], cg->fgfs[Rxz->sgfn], 
                              cg->fgfs[Ryy->sgfn], cg->fgfs[Ryz->sgfn], cg->fgfs[Rzz->sgfn],
                              cg->fgfs[Cons_Ham->sgfn],
                              cg->fgfs[Cons_Px->sgfn], cg->fgfs[Cons_Py->sgfn], cg->fgfs[Cons_Pz->sgfn],
                              cg->fgfs[Cons_Gx->sgfn], cg->fgfs[Cons_Gy->sgfn], cg->fgfs[Cons_Gz->sgfn],
                              Symmetry, lev, ndeps, pre))
        {
          cout << "find NaN in domain: (" 
               << cg->bbox[0] << ":" << cg->bbox[3] << "," 
               << cg->bbox[1] << ":" << cg->bbox[4] << ","
               << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
          ERROR = 1;
        }

        // rk4 substep and boundary
        {
          MyList<var> *varl0 = StateList, *varl = SynchList_pre, *varlrhs = RHSList; // we do not check the correspondence here
          while (varl0)
          {
#if (MRBD == 0)

#ifndef WithShell
            if (lev == 0) // sommerfeld indeed
              f_sommerfeld_routbam(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                   Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], 
                                   Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                                   cg->fgfs[varlrhs->data->sgfn],
                                   cg->fgfs[varl0->data->sgfn], 
                                   varl0->data->propspeed, varl0->data->SoA,
                                   Symmetry);
#endif

#elif (MRBD == 1)
            // sommerfeld indeed
            f_sommerfeld_routbam(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                 Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], 
                                 Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                                 cg->fgfs[varlrhs->data->sgfn],
                                 cg->fgfs[varl0->data->sgfn], 
                                 varl0->data->propspeed, varl0->data->SoA,
                                 Symmetry);
#endif
            f_rungekutta4_rout(cg->shape, dT_lev, 
                               cg->fgfs[varl0->data->sgfn], 
                               cg->fgfs[varl->data->sgfn], 
                               cg->fgfs[varlrhs->data->sgfn],
                               iter_count);

#if (MRBD == 0)

#ifndef WithShell
            if (lev > 0) // fix BD point
#endif
              f_sommerfeld_rout(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], 
                                Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                                dT_lev, 
                                cg->fgfs[phi0->sgfn],
                                cg->fgfs[Lap0->sgfn], 
                                cg->fgfs[varl0->data->sgfn], cg->fgfs[varl->data->sgfn], 
                                varl0->data->SoA,
                                Symmetry, cor);
#elif (MRBD == 2)
            f_sommerfeld_rout(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                              Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], 
                              Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                              dT_lev, 
                              cg->fgfs[phi0->sgfn],
                              cg->fgfs[Lap0->sgfn], 
                              cg->fgfs[varl0->data->sgfn], cg->fgfs[varl->data->sgfn], 
                              varl0->data->SoA,
                              Symmetry, pre);
#endif

            varl0 = varl0->next;
            varl = varl->next;
            varlrhs = varlrhs->next;
          }
        }
        f_lowerboundset(cg->shape, cg->fgfs[phi->sgfn], chitiny);
      }
      if (BP == Pp->data->ble)
        break;
      BP = BP->next;
    }
    Pp = Pp->next;
  }
  // check error information
  {
    int erh = ERROR;
    MPI_Allreduce(&erh, &ERROR, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  }
  if (ERROR)
  {
    Parallel::Dump_Data(GH->PatL[lev], StateList, 0, PhysTime, dT_lev);
    if (myrank == 0)
    {
      if (ErrorMonitor->outfile)
        ErrorMonitor->outfile << "find NaN in state variables at t = " << PhysTime 
                              << ", lev = " << lev << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

#ifdef WithShell
  // evolve Shell Patches
  if (lev == 0)
  {
    sPp = SH->PatL;
    while (sPp)
    {
      MyList<Block> *BP = sPp->data->blb;
      int fngfs = sPp->data->fngfs;
      while (BP)
      {
        Block *cg = BP->data;
        if (myrank == cg->rank)
        {
#if (AGM == 0)
          f_enforce_ga(cg->shape,
                       cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                       cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                       cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                       cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn]);
#endif

          if (f_compute_rhs_Z4c_ss(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                                   cg->fgfs[fngfs + ShellPatch::gx], 
                                   cg->fgfs[fngfs + ShellPatch::gy], 
                                   cg->fgfs[fngfs + ShellPatch::gz],
                                   cg->fgfs[fngfs + ShellPatch::drhodx], 
                                   cg->fgfs[fngfs + ShellPatch::drhody], 
                                   cg->fgfs[fngfs + ShellPatch::drhodz],
                                   cg->fgfs[fngfs + ShellPatch::dsigmadx], 
                                   cg->fgfs[fngfs + ShellPatch::dsigmady], 
                                   cg->fgfs[fngfs + ShellPatch::dsigmadz],
                                   cg->fgfs[fngfs + ShellPatch::dRdx], 
                                   cg->fgfs[fngfs + ShellPatch::dRdy], 
                                   cg->fgfs[fngfs + ShellPatch::dRdz],
                                   cg->fgfs[fngfs + ShellPatch::drhodxx], 
                                   cg->fgfs[fngfs + ShellPatch::drhodxy], 
                                   cg->fgfs[fngfs + ShellPatch::drhodxz],
                                   cg->fgfs[fngfs + ShellPatch::drhodyy], 
                                   cg->fgfs[fngfs + ShellPatch::drhodyz], 
                                   cg->fgfs[fngfs + ShellPatch::drhodzz],
                                   cg->fgfs[fngfs + ShellPatch::dsigmadxx], 
                                   cg->fgfs[fngfs + ShellPatch::dsigmadxy], 
                                   cg->fgfs[fngfs + ShellPatch::dsigmadxz],
                                   cg->fgfs[fngfs + ShellPatch::dsigmadyy], 
                                   cg->fgfs[fngfs + ShellPatch::dsigmadyz], 
                                   cg->fgfs[fngfs + ShellPatch::dsigmadzz],
                                   cg->fgfs[fngfs + ShellPatch::dRdxx], 
                                   cg->fgfs[fngfs + ShellPatch::dRdxy], 
                                   cg->fgfs[fngfs + ShellPatch::dRdxz],
                                   cg->fgfs[fngfs + ShellPatch::dRdyy], 
                                   cg->fgfs[fngfs + ShellPatch::dRdyz], 
                                   cg->fgfs[fngfs + ShellPatch::dRdzz],
                                   cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                                   cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                                   cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                                   cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                                   cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                                   cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                                   cg->fgfs[Lap0->sgfn], 
                                   cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                                   cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn],
                                   cg->fgfs[TZ0->sgfn],
                                   cg->fgfs[phi_rhs->sgfn], cg->fgfs[trK_rhs->sgfn],
                                   cg->fgfs[gxx_rhs->sgfn], cg->fgfs[gxy_rhs->sgfn], cg->fgfs[gxz_rhs->sgfn],
                                   cg->fgfs[gyy_rhs->sgfn], cg->fgfs[gyz_rhs->sgfn], cg->fgfs[gzz_rhs->sgfn],
                                   cg->fgfs[Axx_rhs->sgfn], cg->fgfs[Axy_rhs->sgfn], cg->fgfs[Axz_rhs->sgfn],
                                   cg->fgfs[Ayy_rhs->sgfn], cg->fgfs[Ayz_rhs->sgfn], cg->fgfs[Azz_rhs->sgfn],
                                   cg->fgfs[Gmx_rhs->sgfn], cg->fgfs[Gmy_rhs->sgfn], cg->fgfs[Gmz_rhs->sgfn],
                                   cg->fgfs[Lap_rhs->sgfn], 
                                   cg->fgfs[Sfx_rhs->sgfn], cg->fgfs[Sfy_rhs->sgfn], cg->fgfs[Sfz_rhs->sgfn],
                                   cg->fgfs[dtSfx_rhs->sgfn], cg->fgfs[dtSfy_rhs->sgfn], cg->fgfs[dtSfz_rhs->sgfn],
                                   cg->fgfs[TZ_rhs->sgfn],
                                   cg->fgfs[rho->sgfn], cg->fgfs[Sx->sgfn], cg->fgfs[Sy->sgfn], cg->fgfs[Sz->sgfn],
                                   cg->fgfs[Sxx->sgfn], cg->fgfs[Sxy->sgfn], cg->fgfs[Sxz->sgfn], 
                                   cg->fgfs[Syy->sgfn], cg->fgfs[Syz->sgfn], cg->fgfs[Szz->sgfn],
                                   cg->fgfs[Gamxxx->sgfn], cg->fgfs[Gamxxy->sgfn], cg->fgfs[Gamxxz->sgfn],
                                   cg->fgfs[Gamxyy->sgfn], cg->fgfs[Gamxyz->sgfn], cg->fgfs[Gamxzz->sgfn],
                                   cg->fgfs[Gamyxx->sgfn], cg->fgfs[Gamyxy->sgfn], cg->fgfs[Gamyxz->sgfn],
                                   cg->fgfs[Gamyyy->sgfn], cg->fgfs[Gamyyz->sgfn], cg->fgfs[Gamyzz->sgfn],
                                   cg->fgfs[Gamzxx->sgfn], cg->fgfs[Gamzxy->sgfn], cg->fgfs[Gamzxz->sgfn],
                                   cg->fgfs[Gamzyy->sgfn], cg->fgfs[Gamzyz->sgfn], cg->fgfs[Gamzzz->sgfn],
                                   cg->fgfs[Rxx->sgfn], cg->fgfs[Rxy->sgfn], cg->fgfs[Rxz->sgfn], 
                                   cg->fgfs[Ryy->sgfn], cg->fgfs[Ryz->sgfn], cg->fgfs[Rzz->sgfn],
                                   cg->fgfs[Cons_Ham->sgfn],
                                   cg->fgfs[Cons_Px->sgfn], cg->fgfs[Cons_Py->sgfn], cg->fgfs[Cons_Pz->sgfn],
                                   cg->fgfs[Cons_Gx->sgfn], cg->fgfs[Cons_Gy->sgfn], cg->fgfs[Cons_Gz->sgfn],
                                   Symmetry, lev, numepsh, sPp->data->sst, pre))
          {
            cout << "find NaN in Shell domain: sst = " << sPp->data->sst << ", (" 
                 << cg->bbox[0] << ":" << cg->bbox[3] << ","
                 << cg->bbox[1] << ":" << cg->bbox[4] << "," 
                 << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
            ERROR = 1;
          }

          // rk4 substep and boundary
          {
            MyList<var> *varl0 = StateList, *varl = SynchList_pre, *varlrhs = RHSList; // we do not check the correspondence here
            while (varl0)
            {
              // sommerfeld indeed for outter boudary while fix BD for inner boundary
              f_sommerfeld_routbam_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                      sPp->data->bbox[0], sPp->data->bbox[1], sPp->data->bbox[2], 
                                      sPp->data->bbox[3], sPp->data->bbox[4], sPp->data->bbox[5],
                                      cg->fgfs[varlrhs->data->sgfn],
                                      cg->fgfs[varl0->data->sgfn], 
                                      varl0->data->propspeed, varl0->data->SoA,
                                      Symmetry);

              f_rungekutta4_rout(cg->shape, dT_lev, 
                                 cg->fgfs[varl0->data->sgfn], 
                                 cg->fgfs[varl->data->sgfn], 
                                 cg->fgfs[varlrhs->data->sgfn],
                                 iter_count);

              varl0 = varl0->next;
              varl = varl->next;
              varlrhs = varlrhs->next;
            }
          }
          f_lowerboundset(cg->shape, cg->fgfs[phi->sgfn], chitiny);
        }
        if (BP == sPp->data->ble)
          break;
        BP = BP->next;
      }
      sPp = sPp->next;
    }
  }
  // check error information
  {
    int erh = ERROR;
    MPI_Allreduce(&erh, &ERROR, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  }
  if (ERROR)
  {
    SH->Dump_Data(StateList, 0, PhysTime, dT_lev);
    if (myrank == 0)
    {
      if (ErrorMonitor->outfile)
        ErrorMonitor->outfile << "find NaN in state variables on Shell Patches at t = " << PhysTime << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }
#endif

  Parallel::Sync(GH->PatL[lev], SynchList_pre, Symmetry);

#ifdef WithShell
  if (lev == 0)
  {
    clock_t prev_clock, curr_clock;
    if (myrank == 0)
      curr_clock = clock();
    SH->Synch(SynchList_pre, Symmetry);
    if (myrank == 0)
    {
      prev_clock = curr_clock;
      curr_clock = clock();
      cout << " Shell stuff synchronization used " 
           << (double)(curr_clock - prev_clock) / ((double)CLOCKS_PER_SEC) 
           << " seconds! " << endl;
    }
  }
#endif

  // for black hole position
  if (BH_num > 0 && lev == GH->levels - 1)
  {
    compute_Porg_rhs(Porg0, Porg_rhs, Sfx0, Sfy0, Sfz0, lev);
    for (int ithBH = 0; ithBH < BH_num; ithBH++)
    {
      f_rungekutta4_scalar(dT_lev, Porg0[ithBH][0], Porg[ithBH][0], Porg_rhs[ithBH][0], iter_count);
      f_rungekutta4_scalar(dT_lev, Porg0[ithBH][1], Porg[ithBH][1], Porg_rhs[ithBH][1], iter_count);
      f_rungekutta4_scalar(dT_lev, Porg0[ithBH][2], Porg[ithBH][2], Porg_rhs[ithBH][2], iter_count);
      if (Symmetry > 0)
        Porg[ithBH][2] = fabs(Porg[ithBH][2]);
      if (Symmetry == 2)
      {
        Porg[ithBH][0] = fabs(Porg[ithBH][0]);
        Porg[ithBH][1] = fabs(Porg[ithBH][1]);
      }
      if (!finite(Porg[ithBH][0]) || !finite(Porg[ithBH][1]) || !finite(Porg[ithBH][2]))
      {
        if (ErrorMonitor->outfile)
          ErrorMonitor->outfile << "predictor step finds NaN for BH's position from ("
                                << Porg0[ithBH][0] << "," << Porg0[ithBH][1] << "," << Porg0[ithBH][2] 
                                << ")" << endl;

        MyList<var> *DG_List = new MyList<var>(Sfx0);
        DG_List->insert(Sfy0);
        DG_List->insert(Sfz0);
        Parallel::Dump_Data(GH->PatL[lev], DG_List, 0, PhysTime, dT_lev);
        DG_List->clearList();
      }
    }
  }
  // data analysis part
  // Warning NOTE: the variables1 are used as temp storege room
  if (lev == a_lev)
  {
    AnalysisStuff(lev, dT_lev);
  }
  // corrector
  for (iter_count = 1; iter_count < 4; iter_count++)
  {
    // for RK4: t0, t0+dt/2, t0+dt/2, t0+dt;
    if (iter_count == 1 || iter_count == 3)
      TRK4 += dT_lev / 2;
    Pp = GH->PatL[lev];
    while (Pp)
    {
      MyList<Block> *BP = Pp->data->blb;
      while (BP)
      {
        Block *cg = BP->data;
        if (myrank == cg->rank)
        {
#if (AGM == 0)
          f_enforce_ga(cg->shape,
                       cg->fgfs[gxx->sgfn], cg->fgfs[gxy->sgfn], cg->fgfs[gxz->sgfn], 
                       cg->fgfs[gyy->sgfn], cg->fgfs[gyz->sgfn], cg->fgfs[gzz->sgfn],
                       cg->fgfs[Axx->sgfn], cg->fgfs[Axy->sgfn], cg->fgfs[Axz->sgfn], 
                       cg->fgfs[Ayy->sgfn], cg->fgfs[Ayz->sgfn], cg->fgfs[Azz->sgfn]);
#elif (AGM == 1)
          if (iter_count == 3)
            f_enforce_ga(cg->shape,
                         cg->fgfs[gxx->sgfn], cg->fgfs[gxy->sgfn], cg->fgfs[gxz->sgfn], 
                         cg->fgfs[gyy->sgfn], cg->fgfs[gyz->sgfn], cg->fgfs[gzz->sgfn],
                         cg->fgfs[Axx->sgfn], cg->fgfs[Axy->sgfn], cg->fgfs[Axz->sgfn], 
                         cg->fgfs[Ayy->sgfn], cg->fgfs[Ayz->sgfn], cg->fgfs[Azz->sgfn]);
#endif

          if (f_compute_rhs_Z4c(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                                cg->fgfs[phi->sgfn], cg->fgfs[trK->sgfn],
                                cg->fgfs[gxx->sgfn], cg->fgfs[gxy->sgfn], cg->fgfs[gxz->sgfn], 
                                cg->fgfs[gyy->sgfn], cg->fgfs[gyz->sgfn], cg->fgfs[gzz->sgfn],
                                cg->fgfs[Axx->sgfn], cg->fgfs[Axy->sgfn], cg->fgfs[Axz->sgfn], 
                                cg->fgfs[Ayy->sgfn], cg->fgfs[Ayz->sgfn], cg->fgfs[Azz->sgfn],
                                cg->fgfs[Gmx->sgfn], cg->fgfs[Gmy->sgfn], cg->fgfs[Gmz->sgfn],
                                cg->fgfs[Lap->sgfn], 
                                cg->fgfs[Sfx->sgfn], cg->fgfs[Sfy->sgfn], cg->fgfs[Sfz->sgfn],
                                cg->fgfs[dtSfx->sgfn], cg->fgfs[dtSfy->sgfn], cg->fgfs[dtSfz->sgfn],
                                cg->fgfs[TZ->sgfn],
                                cg->fgfs[phi1->sgfn], cg->fgfs[trK1->sgfn],
                                cg->fgfs[gxx1->sgfn], cg->fgfs[gxy1->sgfn], cg->fgfs[gxz1->sgfn],
                                cg->fgfs[gyy1->sgfn], cg->fgfs[gyz1->sgfn], cg->fgfs[gzz1->sgfn],
                                cg->fgfs[Axx1->sgfn], cg->fgfs[Axy1->sgfn], cg->fgfs[Axz1->sgfn],
                                cg->fgfs[Ayy1->sgfn], cg->fgfs[Ayz1->sgfn], cg->fgfs[Azz1->sgfn],
                                cg->fgfs[Gmx1->sgfn], cg->fgfs[Gmy1->sgfn], cg->fgfs[Gmz1->sgfn],
                                cg->fgfs[Lap1->sgfn], 
                                cg->fgfs[Sfx1->sgfn], cg->fgfs[Sfy1->sgfn], cg->fgfs[Sfz1->sgfn],
                                cg->fgfs[dtSfx1->sgfn], cg->fgfs[dtSfy1->sgfn], cg->fgfs[dtSfz1->sgfn],
                                cg->fgfs[TZ1->sgfn],
                                cg->fgfs[rho->sgfn], cg->fgfs[Sx->sgfn], cg->fgfs[Sy->sgfn], cg->fgfs[Sz->sgfn],
                                cg->fgfs[Sxx->sgfn], cg->fgfs[Sxy->sgfn], cg->fgfs[Sxz->sgfn], 
                                cg->fgfs[Syy->sgfn], cg->fgfs[Syz->sgfn], cg->fgfs[Szz->sgfn],
                                cg->fgfs[Gamxxx->sgfn], cg->fgfs[Gamxxy->sgfn], cg->fgfs[Gamxxz->sgfn],
                                cg->fgfs[Gamxyy->sgfn], cg->fgfs[Gamxyz->sgfn], cg->fgfs[Gamxzz->sgfn],
                                cg->fgfs[Gamyxx->sgfn], cg->fgfs[Gamyxy->sgfn], cg->fgfs[Gamyxz->sgfn],
                                cg->fgfs[Gamyyy->sgfn], cg->fgfs[Gamyyz->sgfn], cg->fgfs[Gamyzz->sgfn],
                                cg->fgfs[Gamzxx->sgfn], cg->fgfs[Gamzxy->sgfn], cg->fgfs[Gamzxz->sgfn],
                                cg->fgfs[Gamzyy->sgfn], cg->fgfs[Gamzyz->sgfn], cg->fgfs[Gamzzz->sgfn],
                                cg->fgfs[Rxx->sgfn], cg->fgfs[Rxy->sgfn], cg->fgfs[Rxz->sgfn], 
                                cg->fgfs[Ryy->sgfn], cg->fgfs[Ryz->sgfn], cg->fgfs[Rzz->sgfn],
                                cg->fgfs[Cons_Ham->sgfn],
                                cg->fgfs[Cons_Px->sgfn], cg->fgfs[Cons_Py->sgfn], cg->fgfs[Cons_Pz->sgfn],
                                cg->fgfs[Cons_Gx->sgfn], cg->fgfs[Cons_Gy->sgfn], cg->fgfs[Cons_Gz->sgfn],
                                Symmetry, lev, ndeps, cor))
          {
            cout << "find NaN in domain: (" 
                 << cg->bbox[0] << ":" << cg->bbox[3] << "," 
                 << cg->bbox[1] << ":" << cg->bbox[4] << ","
                 << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
            ERROR = 1;
          }
          // rk4 substep and boundary
          {
            MyList<var> *varl0 = StateList, *varl = SynchList_pre, *varl1 = SynchList_cor, *varlrhs = RHSList; 
            // we do not check the correspondence here
            
            while (varl0)
            {
#if (MRBD == 0)

#ifndef WithShell
              if (lev == 0) // sommerfeld indeed
                f_sommerfeld_routbam(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                     Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], 
                                     Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                                     cg->fgfs[varl1->data->sgfn],
                                     cg->fgfs[varl->data->sgfn], 
                                     varl0->data->propspeed, varl0->data->SoA,
                                     Symmetry);
#endif

#elif (MRBD == 1)
              // sommerfeld indeed
              f_sommerfeld_routbam(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                   Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], 
                                   Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                                   cg->fgfs[varl1->data->sgfn],
                                   cg->fgfs[varl->data->sgfn], 
                                   varl0->data->propspeed, varl0->data->SoA,
                                   Symmetry);
#endif

              f_rungekutta4_rout(cg->shape, dT_lev, 
                                 cg->fgfs[varl0->data->sgfn], 
                                 cg->fgfs[varl1->data->sgfn], 
                                 cg->fgfs[varlrhs->data->sgfn],
                                 iter_count);

#if (MRBD == 0)

#ifndef WithShell
              if (lev > 0) // fix BD point
#endif
                f_sommerfeld_rout(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                  Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], 
                                  Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                                  dT_lev, 
                                  cg->fgfs[phi0->sgfn],
                                  cg->fgfs[Lap0->sgfn], 
                                  cg->fgfs[varl0->data->sgfn], cg->fgfs[varl1->data->sgfn], 
                                  varl0->data->SoA,
                                  Symmetry, cor);
#elif (MRBD == 2)
              f_sommerfeld_rout(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], 
                                Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                                dT_lev, 
                                cg->fgfs[phi0->sgfn],
                                cg->fgfs[Lap0->sgfn], 
                                cg->fgfs[varl->data->sgfn], cg->fgfs[varl1->data->sgfn], 
                                varl0->data->SoA,
                                Symmetry, cor);
#endif

              varl0 = varl0->next;
              varl = varl->next;
              varl1 = varl1->next;
              varlrhs = varlrhs->next;
            }
          }
          f_lowerboundset(cg->shape, cg->fgfs[phi1->sgfn], chitiny);
        }
        if (BP == Pp->data->ble)
          break;
        BP = BP->next;
      }
      Pp = Pp->next;
    }

    // check error information
    {
      int erh = ERROR;
      MPI_Allreduce(&erh, &ERROR, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    }
    if (ERROR)
    {
      Parallel::Dump_Data(GH->PatL[lev], SynchList_pre, 0, PhysTime, dT_lev);
      if (myrank == 0)
      {
        if (ErrorMonitor->outfile)
          ErrorMonitor->outfile << "find NaN in RK4 substep#" << iter_count 
                                << " variables at t = " << PhysTime 
                                << ", lev = " << lev << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }

#ifdef WithShell
    // evolve Shell Patches
    if (lev == 0)
    {
      sPp = SH->PatL;
      while (sPp)
      {
        MyList<Block> *BP = sPp->data->blb;
        int fngfs = sPp->data->fngfs;
        while (BP)
        {
          Block *cg = BP->data;
          if (myrank == cg->rank)
          {
#if (AGM == 0)
            f_enforce_ga(cg->shape,
                         cg->fgfs[gxx->sgfn], cg->fgfs[gxy->sgfn], cg->fgfs[gxz->sgfn], 
                         cg->fgfs[gyy->sgfn], cg->fgfs[gyz->sgfn], cg->fgfs[gzz->sgfn],
                         cg->fgfs[Axx->sgfn], cg->fgfs[Axy->sgfn], cg->fgfs[Axz->sgfn], 
                         cg->fgfs[Ayy->sgfn], cg->fgfs[Ayz->sgfn], cg->fgfs[Azz->sgfn]);
#elif (AGM == 1)
            if (iter_count == 3)
              f_enforce_ga(cg->shape,
                           cg->fgfs[gxx->sgfn], cg->fgfs[gxy->sgfn], cg->fgfs[gxz->sgfn], 
                           cg->fgfs[gyy->sgfn], cg->fgfs[gyz->sgfn], cg->fgfs[gzz->sgfn],
                           cg->fgfs[Axx->sgfn], cg->fgfs[Axy->sgfn], cg->fgfs[Axz->sgfn], 
                           cg->fgfs[Ayy->sgfn], cg->fgfs[Ayz->sgfn], cg->fgfs[Azz->sgfn]);
#endif

            if (f_compute_rhs_Z4c_ss(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                                     cg->fgfs[fngfs + ShellPatch::gx], 
                                     cg->fgfs[fngfs + ShellPatch::gy], 
                                     cg->fgfs[fngfs + ShellPatch::gz],
                                     cg->fgfs[fngfs + ShellPatch::drhodx], 
                                     cg->fgfs[fngfs + ShellPatch::drhody], 
                                     cg->fgfs[fngfs + ShellPatch::drhodz],
                                     cg->fgfs[fngfs + ShellPatch::dsigmadx], 
                                     cg->fgfs[fngfs + ShellPatch::dsigmady], 
                                     cg->fgfs[fngfs + ShellPatch::dsigmadz],
                                     cg->fgfs[fngfs + ShellPatch::dRdx], 
                                     cg->fgfs[fngfs + ShellPatch::dRdy], 
                                     cg->fgfs[fngfs + ShellPatch::dRdz],
                                     cg->fgfs[fngfs + ShellPatch::drhodxx], 
                                     cg->fgfs[fngfs + ShellPatch::drhodxy], 
                                     cg->fgfs[fngfs + ShellPatch::drhodxz],
                                     cg->fgfs[fngfs + ShellPatch::drhodyy], 
                                     cg->fgfs[fngfs + ShellPatch::drhodyz], 
                                     cg->fgfs[fngfs + ShellPatch::drhodzz],
                                     cg->fgfs[fngfs + ShellPatch::dsigmadxx], 
                                     cg->fgfs[fngfs + ShellPatch::dsigmadxy], 
                                     cg->fgfs[fngfs + ShellPatch::dsigmadxz],
                                     cg->fgfs[fngfs + ShellPatch::dsigmadyy], 
                                     cg->fgfs[fngfs + ShellPatch::dsigmadyz], 
                                     cg->fgfs[fngfs + ShellPatch::dsigmadzz],
                                     cg->fgfs[fngfs + ShellPatch::dRdxx], 
                                     cg->fgfs[fngfs + ShellPatch::dRdxy], 
                                     cg->fgfs[fngfs + ShellPatch::dRdxz],
                                     cg->fgfs[fngfs + ShellPatch::dRdyy], 
                                     cg->fgfs[fngfs + ShellPatch::dRdyz], 
                                     cg->fgfs[fngfs + ShellPatch::dRdzz],
                                     cg->fgfs[phi->sgfn], cg->fgfs[trK->sgfn],
                                     cg->fgfs[gxx->sgfn], cg->fgfs[gxy->sgfn], cg->fgfs[gxz->sgfn], 
                                     cg->fgfs[gyy->sgfn], cg->fgfs[gyz->sgfn], cg->fgfs[gzz->sgfn],
                                     cg->fgfs[Axx->sgfn], cg->fgfs[Axy->sgfn], cg->fgfs[Axz->sgfn], 
                                     cg->fgfs[Ayy->sgfn], cg->fgfs[Ayz->sgfn], cg->fgfs[Azz->sgfn],
                                     cg->fgfs[Gmx->sgfn], cg->fgfs[Gmy->sgfn], cg->fgfs[Gmz->sgfn],
                                     cg->fgfs[Lap->sgfn], 
                                     cg->fgfs[Sfx->sgfn], cg->fgfs[Sfy->sgfn], cg->fgfs[Sfz->sgfn],
                                     cg->fgfs[dtSfx->sgfn], cg->fgfs[dtSfy->sgfn], cg->fgfs[dtSfz->sgfn],
                                     cg->fgfs[TZ->sgfn],
                                     cg->fgfs[phi1->sgfn], cg->fgfs[trK1->sgfn],
                                     cg->fgfs[gxx1->sgfn], cg->fgfs[gxy1->sgfn], cg->fgfs[gxz1->sgfn],
                                     cg->fgfs[gyy1->sgfn], cg->fgfs[gyz1->sgfn], cg->fgfs[gzz1->sgfn],
                                     cg->fgfs[Axx1->sgfn], cg->fgfs[Axy1->sgfn], cg->fgfs[Axz1->sgfn],
                                     cg->fgfs[Ayy1->sgfn], cg->fgfs[Ayz1->sgfn], cg->fgfs[Azz1->sgfn],
                                     cg->fgfs[Gmx1->sgfn], cg->fgfs[Gmy1->sgfn], cg->fgfs[Gmz1->sgfn],
                                     cg->fgfs[Lap1->sgfn], 
                                     cg->fgfs[Sfx1->sgfn], cg->fgfs[Sfy1->sgfn], cg->fgfs[Sfz1->sgfn],
                                     cg->fgfs[dtSfx1->sgfn], cg->fgfs[dtSfy1->sgfn], cg->fgfs[dtSfz1->sgfn],
                                     cg->fgfs[TZ1->sgfn],
                                     cg->fgfs[rho->sgfn], cg->fgfs[Sx->sgfn], cg->fgfs[Sy->sgfn], cg->fgfs[Sz->sgfn],
                                     cg->fgfs[Sxx->sgfn], cg->fgfs[Sxy->sgfn], cg->fgfs[Sxz->sgfn], 
                                     cg->fgfs[Syy->sgfn], cg->fgfs[Syz->sgfn], cg->fgfs[Szz->sgfn],
                                     cg->fgfs[Gamxxx->sgfn], cg->fgfs[Gamxxy->sgfn], cg->fgfs[Gamxxz->sgfn],
                                     cg->fgfs[Gamxyy->sgfn], cg->fgfs[Gamxyz->sgfn], cg->fgfs[Gamxzz->sgfn],
                                     cg->fgfs[Gamyxx->sgfn], cg->fgfs[Gamyxy->sgfn], cg->fgfs[Gamyxz->sgfn],
                                     cg->fgfs[Gamyyy->sgfn], cg->fgfs[Gamyyz->sgfn], cg->fgfs[Gamyzz->sgfn],
                                     cg->fgfs[Gamzxx->sgfn], cg->fgfs[Gamzxy->sgfn], cg->fgfs[Gamzxz->sgfn],
                                     cg->fgfs[Gamzyy->sgfn], cg->fgfs[Gamzyz->sgfn], cg->fgfs[Gamzzz->sgfn],
                                     cg->fgfs[Rxx->sgfn], cg->fgfs[Rxy->sgfn], cg->fgfs[Rxz->sgfn], 
                                     cg->fgfs[Ryy->sgfn], cg->fgfs[Ryz->sgfn], cg->fgfs[Rzz->sgfn],
                                     cg->fgfs[Cons_Ham->sgfn],
                                     cg->fgfs[Cons_Px->sgfn], cg->fgfs[Cons_Py->sgfn], cg->fgfs[Cons_Pz->sgfn],
                                     cg->fgfs[Cons_Gx->sgfn], cg->fgfs[Cons_Gy->sgfn], cg->fgfs[Cons_Gz->sgfn],
                                     Symmetry, lev, numepsh, sPp->data->sst, cor))
            {
              cout << "find NaN in Shell domain: sst = " << sPp->data->sst << ", (" 
                   << cg->bbox[0] << ":" << cg->bbox[3] << ","
                   << cg->bbox[1] << ":" << cg->bbox[4] << "," 
                   << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
              ERROR = 1;
            }
            // rk4 substep and boundary
            {
              MyList<var> *varl0 = StateList, *varl = SynchList_pre, *varl1 = SynchList_cor, *varlrhs = RHSList; 
              // we do not check the correspondence here
              
              while (varl0)
              {
                // sommerfeld indeed for outter boudary while fix BD for inner boundary
                f_sommerfeld_routbam_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                        sPp->data->bbox[0], sPp->data->bbox[1], sPp->data->bbox[2], 
                                        sPp->data->bbox[3], sPp->data->bbox[4], sPp->data->bbox[5],
                                        cg->fgfs[varl1->data->sgfn],
                                        cg->fgfs[varl->data->sgfn], 
                                        varl0->data->propspeed, varl0->data->SoA,
                                        Symmetry);

                f_rungekutta4_rout(cg->shape, dT_lev, 
                                   cg->fgfs[varl0->data->sgfn], 
                                   cg->fgfs[varl1->data->sgfn], 
                                   cg->fgfs[varlrhs->data->sgfn],
                                   iter_count);

                varl0 = varl0->next;
                varl = varl->next;
                varl1 = varl1->next;
                varlrhs = varlrhs->next;
              }
            }
            f_lowerboundset(cg->shape, cg->fgfs[phi1->sgfn], chitiny);
          }
          if (BP == sPp->data->ble)
            break;
          BP = BP->next;
        }
        sPp = sPp->next;
      }
    }
    // check error information
    {
      int erh = ERROR;
      MPI_Allreduce(&erh, &ERROR, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    }
    if (ERROR)
    {
      SH->Dump_Data(SynchList_pre, 0, PhysTime, dT_lev);
      if (myrank == 0)
      {
        if (ErrorMonitor->outfile)
          ErrorMonitor->outfile << "find NaN on Shell Patches in RK4 substep#" << iter_count 
                                << " variables at t = " << PhysTime << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }
#endif

    Parallel::Sync(GH->PatL[lev], SynchList_cor, Symmetry);

#ifdef WithShell
    if (lev == 0)
    {
      clock_t prev_clock, curr_clock;
      if (myrank == 0)
        curr_clock = clock();
      SH->Synch(SynchList_cor, Symmetry);
      if (myrank == 0)
      {
        prev_clock = curr_clock;
        curr_clock = clock();
        cout << " Shell stuff synchronization used " 
             << (double)(curr_clock - prev_clock) / ((double)CLOCKS_PER_SEC) 
             << " seconds! " << endl;
      }
    }
#endif
    // for black hole position
    if (BH_num > 0 && lev == GH->levels - 1)
    {
      compute_Porg_rhs(Porg, Porg1, Sfx, Sfy, Sfz, lev);
      for (int ithBH = 0; ithBH < BH_num; ithBH++)
      {
        f_rungekutta4_scalar(dT_lev, Porg0[ithBH][0], Porg1[ithBH][0], Porg_rhs[ithBH][0], iter_count);
        f_rungekutta4_scalar(dT_lev, Porg0[ithBH][1], Porg1[ithBH][1], Porg_rhs[ithBH][1], iter_count);
        f_rungekutta4_scalar(dT_lev, Porg0[ithBH][2], Porg1[ithBH][2], Porg_rhs[ithBH][2], iter_count);
        if (Symmetry > 0)
          Porg1[ithBH][2] = fabs(Porg1[ithBH][2]);
        if (Symmetry == 2)
        {
          Porg1[ithBH][0] = fabs(Porg1[ithBH][0]);
          Porg1[ithBH][1] = fabs(Porg1[ithBH][1]);
        }
        if (!finite(Porg1[ithBH][0]) || !finite(Porg1[ithBH][1]) || !finite(Porg1[ithBH][2]))
        {
          if (ErrorMonitor->outfile)
            ErrorMonitor->outfile << iter_count << " corrector step finds NaN for BH's position from ("
                                  << Porg[ithBH][0] << "," << Porg[ithBH][1] << "," << Porg[ithBH][2] 
                                  << ")" << endl;

          MyList<var> *DG_List = new MyList<var>(Sfx0);
          DG_List->insert(Sfx0);
          DG_List->insert(Sfy0);
          DG_List->insert(Sfz0);
          Parallel::Dump_Data(GH->PatL[lev], DG_List, 0, PhysTime, dT_lev);
          DG_List->clearList();
        }
      }
    }
    // swap time level
    if (iter_count < 3)
    {
      Pp = GH->PatL[lev];
      while (Pp)
      {
        MyList<Block> *BP = Pp->data->blb;
        while (BP)
        {
          Block *cg = BP->data;
          cg->swapList(SynchList_pre, SynchList_cor, myrank);
          if (BP == Pp->data->ble)
            break;
          BP = BP->next;
        }
        Pp = Pp->next;
      }
#ifdef WithShell
      if (lev == 0)
      {
        sPp = SH->PatL;
        while (sPp)
        {
          MyList<Block> *BP = sPp->data->blb;
          while (BP)
          {
            Block *cg = BP->data;
            cg->swapList(SynchList_pre, SynchList_cor, myrank);
            if (BP == sPp->data->ble)
              break;
            BP = BP->next;
          }
          sPp = sPp->next;
        }
      }
#endif
      // for black hole position
      if (BH_num > 0 && lev == GH->levels - 1)
      {
        for (int ithBH = 0; ithBH < BH_num; ithBH++)
        {
          Porg[ithBH][0] = Porg1[ithBH][0];
          Porg[ithBH][1] = Porg1[ithBH][1];
          Porg[ithBH][2] = Porg1[ithBH][2];
        }
      }
    }
  }

#if (RPS == 0)
  // mesh refinement boundary part
  RestrictProlong(lev, YN, BB);

#ifdef WithShell
  if (lev == 0)
  {
    clock_t prev_clock, curr_clock;
    if (myrank == 0)
      curr_clock = clock();
    SH->CS_Inter(SynchList_cor, Symmetry);
    if (myrank == 0)
    {
      prev_clock = curr_clock;
      curr_clock = clock();
      cout << " CS_Inter used " << (double)(curr_clock - prev_clock) / ((double)CLOCKS_PER_SEC) 
           << " seconds! " << endl;
    }
  }
#endif

#endif
  // note the data structure before update
  // SynchList_cor 1   -----------
  //
  // StateList     0   -----------
  //
  // OldStateList  old -----------
  // update
  Pp = GH->PatL[lev];
  while (Pp)
  {
    MyList<Block> *BP = Pp->data->blb;
    while (BP)
    {
      Block *cg = BP->data;
      cg->swapList(StateList, SynchList_cor, myrank);
      cg->swapList(OldStateList, SynchList_cor, myrank);
      if (BP == Pp->data->ble)
        break;
      BP = BP->next;
    }
    Pp = Pp->next;
  }
#ifdef WithShell
  if (lev == 0)
  {
    sPp = SH->PatL;
    while (sPp)
    {
      MyList<Block> *BP = sPp->data->blb;
      while (BP)
      {
        Block *cg = BP->data;
        cg->swapList(StateList, SynchList_cor, myrank);
        cg->swapList(OldStateList, SynchList_cor, myrank);
        if (BP == sPp->data->ble)
          break;
        BP = BP->next;
      }
      sPp = sPp->next;
    }
  }
#endif
  // for black hole position
  if (BH_num > 0 && lev == GH->levels - 1)
  {
    for (int ithBH = 0; ithBH < BH_num; ithBH++)
    {
      Porg0[ithBH][0] = Porg1[ithBH][0];
      Porg0[ithBH][1] = Porg1[ithBH][1];
      Porg0[ithBH][2] = Porg1[ithBH][2];
    }
  }
}
#else
// for constraint preserving boundary (CPBC)
#ifndef WithShell
#error "CPBC only supports Shell"
#endif

// 0: extroplate rhs, 1: extroplate variable
// 2: extroplate variable but before RHS calculation
#define EXTO 1

// #define SMOOTHSHELL

// change chi based on chitiny or not: 0: yes; 1: no
#define chinot 0
void Z4c_class::Step(int lev, int YN)
{
  //    Check_extrop();
  double dT_lev = dT * pow(0.5, Mymax(lev, trfls));
#ifdef With_AHF
  AH_Step_Find(lev, dT_lev);
#endif
  bool BB = fgt(PhysTime, StartTime, dT_lev / 2);
  double fbeps = -0.1;
  double ndeps = numepss;
  if (lev < GH->movls)
    ndeps = numepsb;
  double TRK4 = PhysTime;
  int iter_count = 0; // count RK4 substeps
  int pre = 0, cor = 1;
  int ERROR = 0;

  MyList<ss_patch> *sPp;
  // Predictor
  MyList<Patch> *Pp = GH->PatL[lev];
  while (Pp)
  {
    MyList<Block> *BP = Pp->data->blb;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
#if (AGM == 0)
        f_enforce_ga(cg->shape,
                     cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                     cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                     cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                     cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn]);
#endif

#if (chinot == 0)
        if (f_compute_rhs_Z4c(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                              cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                              cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                              cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                              cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                              cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                              cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                              cg->fgfs[Lap0->sgfn], 
                              cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                              cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn],
                              cg->fgfs[TZ0->sgfn],
                              cg->fgfs[phi_rhs->sgfn], cg->fgfs[trK_rhs->sgfn],
                              cg->fgfs[gxx_rhs->sgfn], cg->fgfs[gxy_rhs->sgfn], cg->fgfs[gxz_rhs->sgfn],
                              cg->fgfs[gyy_rhs->sgfn], cg->fgfs[gyz_rhs->sgfn], cg->fgfs[gzz_rhs->sgfn],
                              cg->fgfs[Axx_rhs->sgfn], cg->fgfs[Axy_rhs->sgfn], cg->fgfs[Axz_rhs->sgfn],
                              cg->fgfs[Ayy_rhs->sgfn], cg->fgfs[Ayz_rhs->sgfn], cg->fgfs[Azz_rhs->sgfn],
                              cg->fgfs[Gmx_rhs->sgfn], cg->fgfs[Gmy_rhs->sgfn], cg->fgfs[Gmz_rhs->sgfn],
                              cg->fgfs[Lap_rhs->sgfn], 
                              cg->fgfs[Sfx_rhs->sgfn], cg->fgfs[Sfy_rhs->sgfn], cg->fgfs[Sfz_rhs->sgfn],
                              cg->fgfs[dtSfx_rhs->sgfn], cg->fgfs[dtSfy_rhs->sgfn], cg->fgfs[dtSfz_rhs->sgfn],
                              cg->fgfs[TZ_rhs->sgfn],
                              cg->fgfs[rho->sgfn], cg->fgfs[Sx->sgfn], cg->fgfs[Sy->sgfn], cg->fgfs[Sz->sgfn],
                              cg->fgfs[Sxx->sgfn], cg->fgfs[Sxy->sgfn], cg->fgfs[Sxz->sgfn], 
                              cg->fgfs[Syy->sgfn], cg->fgfs[Syz->sgfn], cg->fgfs[Szz->sgfn],
                              cg->fgfs[Gamxxx->sgfn], cg->fgfs[Gamxxy->sgfn], cg->fgfs[Gamxxz->sgfn],
                              cg->fgfs[Gamxyy->sgfn], cg->fgfs[Gamxyz->sgfn], cg->fgfs[Gamxzz->sgfn],
                              cg->fgfs[Gamyxx->sgfn], cg->fgfs[Gamyxy->sgfn], cg->fgfs[Gamyxz->sgfn],
                              cg->fgfs[Gamyyy->sgfn], cg->fgfs[Gamyyz->sgfn], cg->fgfs[Gamyzz->sgfn],
                              cg->fgfs[Gamzxx->sgfn], cg->fgfs[Gamzxy->sgfn], cg->fgfs[Gamzxz->sgfn],
                              cg->fgfs[Gamzyy->sgfn], cg->fgfs[Gamzyz->sgfn], cg->fgfs[Gamzzz->sgfn],
                              cg->fgfs[Rxx->sgfn], cg->fgfs[Rxy->sgfn], cg->fgfs[Rxz->sgfn], 
                              cg->fgfs[Ryy->sgfn], cg->fgfs[Ryz->sgfn], cg->fgfs[Rzz->sgfn],
                              cg->fgfs[Cons_Ham->sgfn],
                              cg->fgfs[Cons_Px->sgfn], cg->fgfs[Cons_Py->sgfn], cg->fgfs[Cons_Pz->sgfn],
                              cg->fgfs[Cons_Gx->sgfn], cg->fgfs[Cons_Gy->sgfn], cg->fgfs[Cons_Gz->sgfn],
                              Symmetry, lev, ndeps, pre))
        {
          cout << "find NaN in domain: (" 
               << cg->bbox[0] << ":" << cg->bbox[3] << "," 
               << cg->bbox[1] << ":" << cg->bbox[4] << ","
               << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
          ERROR = 1;
        }
#else
        if (f_compute_rhs_Z4cnot(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                                 cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                                 cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                                 cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                                 cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                                 cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                                 cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                                 cg->fgfs[Lap0->sgfn], 
                                 cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                                 cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn],
                                 cg->fgfs[TZ0->sgfn],
                                 cg->fgfs[phi_rhs->sgfn], cg->fgfs[trK_rhs->sgfn],
                                 cg->fgfs[gxx_rhs->sgfn], cg->fgfs[gxy_rhs->sgfn], cg->fgfs[gxz_rhs->sgfn],
                                 cg->fgfs[gyy_rhs->sgfn], cg->fgfs[gyz_rhs->sgfn], cg->fgfs[gzz_rhs->sgfn],
                                 cg->fgfs[Axx_rhs->sgfn], cg->fgfs[Axy_rhs->sgfn], cg->fgfs[Axz_rhs->sgfn],
                                 cg->fgfs[Ayy_rhs->sgfn], cg->fgfs[Ayz_rhs->sgfn], cg->fgfs[Azz_rhs->sgfn],
                                 cg->fgfs[Gmx_rhs->sgfn], cg->fgfs[Gmy_rhs->sgfn], cg->fgfs[Gmz_rhs->sgfn],
                                 cg->fgfs[Lap_rhs->sgfn], 
                                 cg->fgfs[Sfx_rhs->sgfn], cg->fgfs[Sfy_rhs->sgfn], cg->fgfs[Sfz_rhs->sgfn],
                                 cg->fgfs[dtSfx_rhs->sgfn], cg->fgfs[dtSfy_rhs->sgfn], cg->fgfs[dtSfz_rhs->sgfn],
                                 cg->fgfs[TZ_rhs->sgfn],
                                 cg->fgfs[rho->sgfn], cg->fgfs[Sx->sgfn], cg->fgfs[Sy->sgfn], cg->fgfs[Sz->sgfn],
                                 cg->fgfs[Sxx->sgfn], cg->fgfs[Sxy->sgfn], cg->fgfs[Sxz->sgfn], 
                                 cg->fgfs[Syy->sgfn], cg->fgfs[Syz->sgfn], cg->fgfs[Szz->sgfn],
                                 cg->fgfs[Gamxxx->sgfn], cg->fgfs[Gamxxy->sgfn], cg->fgfs[Gamxxz->sgfn],
                                 cg->fgfs[Gamxyy->sgfn], cg->fgfs[Gamxyz->sgfn], cg->fgfs[Gamxzz->sgfn],
                                 cg->fgfs[Gamyxx->sgfn], cg->fgfs[Gamyxy->sgfn], cg->fgfs[Gamyxz->sgfn],
                                 cg->fgfs[Gamyyy->sgfn], cg->fgfs[Gamyyz->sgfn], cg->fgfs[Gamyzz->sgfn],
                                 cg->fgfs[Gamzxx->sgfn], cg->fgfs[Gamzxy->sgfn], cg->fgfs[Gamzxz->sgfn],
                                 cg->fgfs[Gamzyy->sgfn], cg->fgfs[Gamzyz->sgfn], cg->fgfs[Gamzzz->sgfn],
                                 cg->fgfs[Rxx->sgfn], cg->fgfs[Rxy->sgfn], cg->fgfs[Rxz->sgfn], 
                                 cg->fgfs[Ryy->sgfn], cg->fgfs[Ryz->sgfn], cg->fgfs[Rzz->sgfn],
                                 cg->fgfs[Cons_Ham->sgfn],
                                 cg->fgfs[Cons_Px->sgfn], cg->fgfs[Cons_Py->sgfn], cg->fgfs[Cons_Pz->sgfn],
                                 cg->fgfs[Cons_Gx->sgfn], cg->fgfs[Cons_Gy->sgfn], cg->fgfs[Cons_Gz->sgfn],
                                 Symmetry, lev, ndeps, pre, chitiny))
        {
          cout << "find NaN in domain: (" 
               << cg->bbox[0] << ":" << cg->bbox[3] << "," 
               << cg->bbox[1] << ":" << cg->bbox[4] << ","
               << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
          ERROR = 1;
        }
#endif
        // rk4 substep and boundary
        {
          MyList<var> *varl0 = StateList, *varl = SynchList_pre, *varlrhs = RHSList; 
          // we do not check the correspondence here
          
          while (varl0)
          {
#if (MRBD == 1)
            // sommerfeld indeed
            f_sommerfeld_routbam(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                 Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], 
                                 Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                                 cg->fgfs[varlrhs->data->sgfn],
                                 cg->fgfs[varl0->data->sgfn], 
                                 varl0->data->propspeed, varl0->data->SoA,
                                 Symmetry);
#endif
            f_rungekutta4_rout(cg->shape, dT_lev, 
                               cg->fgfs[varl0->data->sgfn], 
                               cg->fgfs[varl->data->sgfn], 
                               cg->fgfs[varlrhs->data->sgfn],
                               iter_count);
#if (MRBD == 0)
            f_sommerfeld_rout(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                              Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], 
                              Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                              dT_lev, 
                              cg->fgfs[phi0->sgfn],
                              cg->fgfs[Lap0->sgfn], 
                              cg->fgfs[varl0->data->sgfn], cg->fgfs[varl->data->sgfn], 
                              varl0->data->SoA,
                              Symmetry, cor);
#elif (MRBD == 2)
            f_sommerfeld_rout(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                              Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], 
                              Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                              dT_lev, 
                              cg->fgfs[phi0->sgfn],
                              cg->fgfs[Lap0->sgfn], 
                              cg->fgfs[varl0->data->sgfn], cg->fgfs[varl->data->sgfn], 
                              varl0->data->SoA,
                              Symmetry, pre);
#endif

            varl0 = varl0->next;
            varl = varl->next;
            varlrhs = varlrhs->next;
          }
        }
#if (chinot == 0)
        f_lowerboundset(cg->shape, cg->fgfs[phi->sgfn], chitiny);
#endif
      }
      if (BP == Pp->data->ble)
        break;
      BP = BP->next;
    }
    Pp = Pp->next;
  }
#if 0 
// check rhs
  {
         Parallel::Dump_Data(GH->PatL[lev],RHSList,0,PhysTime,dT_lev);
         if(myrank == 0)
	 {
            cout<<"check irhs for box"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	 }
  }
#endif
  // check error information
  {
    int erh = ERROR;
    MPI_Allreduce(&erh, &ERROR, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  }
  if (ERROR)
  {
    Parallel::Dump_Data(GH->PatL[lev], StateList, 0, PhysTime, dT_lev);
    if (myrank == 0)
    {
      if (ErrorMonitor->outfile)
        ErrorMonitor->outfile << "find NaN in state variables at t = " << PhysTime 
                              << ", lev = " << lev << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  // evolve Shell Patches
  if (lev == 0)
  {
    sPp = SH->PatL;
    while (sPp)
    {
      MyList<Block> *BP = sPp->data->blb;
      int fngfs = sPp->data->fngfs;
      while (BP)
      {
        Block *cg = BP->data;
        if (myrank == cg->rank)
        {
#if (EXTO == 2)
          // extroplate variable itself
          f_david_milton_extroplate_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                       cg->fgfs[TZ0->sgfn], cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                                       cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn],
                                       cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                                       cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn],
                                       cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                                       cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                                       cg->fgfs[Lap0->sgfn], 
                                       cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                                       cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn], 
                                       sPp->data->bbox[2], sPp->data->bbox[5]);
#endif

#if (AGM == 0)
          f_enforce_ga(cg->shape,
                       cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                       cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                       cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                       cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn]);
#endif

          if (f_compute_rhs_Z4c_ss(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                                   cg->fgfs[fngfs + ShellPatch::gx], 
                                   cg->fgfs[fngfs + ShellPatch::gy], 
                                   cg->fgfs[fngfs + ShellPatch::gz],
                                   cg->fgfs[fngfs + ShellPatch::drhodx], 
                                   cg->fgfs[fngfs + ShellPatch::drhody], 
                                   cg->fgfs[fngfs + ShellPatch::drhodz],
                                   cg->fgfs[fngfs + ShellPatch::dsigmadx], 
                                   cg->fgfs[fngfs + ShellPatch::dsigmady], 
                                   cg->fgfs[fngfs + ShellPatch::dsigmadz],
                                   cg->fgfs[fngfs + ShellPatch::dRdx], 
                                   cg->fgfs[fngfs + ShellPatch::dRdy], 
                                   cg->fgfs[fngfs + ShellPatch::dRdz],
                                   cg->fgfs[fngfs + ShellPatch::drhodxx], 
                                   cg->fgfs[fngfs + ShellPatch::drhodxy], 
                                   cg->fgfs[fngfs + ShellPatch::drhodxz],
                                   cg->fgfs[fngfs + ShellPatch::drhodyy], 
                                   cg->fgfs[fngfs + ShellPatch::drhodyz], 
                                   cg->fgfs[fngfs + ShellPatch::drhodzz],
                                   cg->fgfs[fngfs + ShellPatch::dsigmadxx], 
                                   cg->fgfs[fngfs + ShellPatch::dsigmadxy], 
                                   cg->fgfs[fngfs + ShellPatch::dsigmadxz],
                                   cg->fgfs[fngfs + ShellPatch::dsigmadyy], 
                                   cg->fgfs[fngfs + ShellPatch::dsigmadyz], 
                                   cg->fgfs[fngfs + ShellPatch::dsigmadzz],
                                   cg->fgfs[fngfs + ShellPatch::dRdxx], 
                                   cg->fgfs[fngfs + ShellPatch::dRdxy], 
                                   cg->fgfs[fngfs + ShellPatch::dRdxz],
                                   cg->fgfs[fngfs + ShellPatch::dRdyy], 
                                   cg->fgfs[fngfs + ShellPatch::dRdyz], 
                                   cg->fgfs[fngfs + ShellPatch::dRdzz],
                                   cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                                   cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                                   cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                                   cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                                   cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                                   cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                                   cg->fgfs[Lap0->sgfn], 
                                   cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                                   cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn],
                                   cg->fgfs[TZ0->sgfn],
                                   cg->fgfs[phi_rhs->sgfn], cg->fgfs[trK_rhs->sgfn],
                                   cg->fgfs[gxx_rhs->sgfn], cg->fgfs[gxy_rhs->sgfn], cg->fgfs[gxz_rhs->sgfn],
                                   cg->fgfs[gyy_rhs->sgfn], cg->fgfs[gyz_rhs->sgfn], cg->fgfs[gzz_rhs->sgfn],
                                   cg->fgfs[Axx_rhs->sgfn], cg->fgfs[Axy_rhs->sgfn], cg->fgfs[Axz_rhs->sgfn],
                                   cg->fgfs[Ayy_rhs->sgfn], cg->fgfs[Ayz_rhs->sgfn], cg->fgfs[Azz_rhs->sgfn],
                                   cg->fgfs[Gmx_rhs->sgfn], cg->fgfs[Gmy_rhs->sgfn], cg->fgfs[Gmz_rhs->sgfn],
                                   cg->fgfs[Lap_rhs->sgfn], 
                                   cg->fgfs[Sfx_rhs->sgfn], cg->fgfs[Sfy_rhs->sgfn], cg->fgfs[Sfz_rhs->sgfn],
                                   cg->fgfs[dtSfx_rhs->sgfn], cg->fgfs[dtSfy_rhs->sgfn], cg->fgfs[dtSfz_rhs->sgfn],
                                   cg->fgfs[TZ_rhs->sgfn],
                                   cg->fgfs[rho->sgfn], cg->fgfs[Sx->sgfn], cg->fgfs[Sy->sgfn], cg->fgfs[Sz->sgfn],
                                   cg->fgfs[Sxx->sgfn], cg->fgfs[Sxy->sgfn], cg->fgfs[Sxz->sgfn], 
                                   cg->fgfs[Syy->sgfn], cg->fgfs[Syz->sgfn], cg->fgfs[Szz->sgfn],
                                   cg->fgfs[Gamxxx->sgfn], cg->fgfs[Gamxxy->sgfn], cg->fgfs[Gamxxz->sgfn],
                                   cg->fgfs[Gamxyy->sgfn], cg->fgfs[Gamxyz->sgfn], cg->fgfs[Gamxzz->sgfn],
                                   cg->fgfs[Gamyxx->sgfn], cg->fgfs[Gamyxy->sgfn], cg->fgfs[Gamyxz->sgfn],
                                   cg->fgfs[Gamyyy->sgfn], cg->fgfs[Gamyyz->sgfn], cg->fgfs[Gamyzz->sgfn],
                                   cg->fgfs[Gamzxx->sgfn], cg->fgfs[Gamzxy->sgfn], cg->fgfs[Gamzxz->sgfn],
                                   cg->fgfs[Gamzyy->sgfn], cg->fgfs[Gamzyz->sgfn], cg->fgfs[Gamzzz->sgfn],
                                   cg->fgfs[Rxx->sgfn], cg->fgfs[Rxy->sgfn], cg->fgfs[Rxz->sgfn], 
                                   cg->fgfs[Ryy->sgfn], cg->fgfs[Ryz->sgfn], cg->fgfs[Rzz->sgfn],
                                   cg->fgfs[Cons_Ham->sgfn],
                                   cg->fgfs[Cons_Px->sgfn], cg->fgfs[Cons_Py->sgfn], cg->fgfs[Cons_Pz->sgfn],
                                   cg->fgfs[Cons_Gx->sgfn], cg->fgfs[Cons_Gy->sgfn], cg->fgfs[Cons_Gz->sgfn],
                                   Symmetry, lev, fbeps, sPp->data->sst, pre))
          {
            cout << "find NaN in Shell domain: sst = " << sPp->data->sst << ", (" 
                 << cg->bbox[0] << ":" << cg->bbox[3] << ","
                 << cg->bbox[1] << ":" << cg->bbox[4] << "," 
                 << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
            ERROR = 1;
          }

          // rk4 substep and boundary
          {
            // CPBC indeed for outter boudary while fix BD for inner boundary
            f_david_milton_cpbc_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                   cg->fgfs[fngfs + ShellPatch::gx], 
                                   cg->fgfs[fngfs + ShellPatch::gy], 
                                   cg->fgfs[fngfs + ShellPatch::gz],
                                   cg->fgfs[fngfs + ShellPatch::drhodx], 
                                   cg->fgfs[fngfs + ShellPatch::drhody], 
                                   cg->fgfs[fngfs + ShellPatch::drhodz],
                                   cg->fgfs[fngfs + ShellPatch::dsigmadx], 
                                   cg->fgfs[fngfs + ShellPatch::dsigmady], 
                                   cg->fgfs[fngfs + ShellPatch::dsigmadz],
                                   cg->fgfs[fngfs + ShellPatch::dRdx], 
                                   cg->fgfs[fngfs + ShellPatch::dRdy], 
                                   cg->fgfs[fngfs + ShellPatch::dRdz],
                                   cg->fgfs[fngfs + ShellPatch::drhodxx], 
                                   cg->fgfs[fngfs + ShellPatch::drhodxy], 
                                   cg->fgfs[fngfs + ShellPatch::drhodxz],
                                   cg->fgfs[fngfs + ShellPatch::drhodyy], 
                                   cg->fgfs[fngfs + ShellPatch::drhodyz], 
                                   cg->fgfs[fngfs + ShellPatch::drhodzz],
                                   cg->fgfs[fngfs + ShellPatch::dsigmadxx], 
                                   cg->fgfs[fngfs + ShellPatch::dsigmadxy], 
                                   cg->fgfs[fngfs + ShellPatch::dsigmadxz],
                                   cg->fgfs[fngfs + ShellPatch::dsigmadyy], 
                                   cg->fgfs[fngfs + ShellPatch::dsigmadyz], 
                                   cg->fgfs[fngfs + ShellPatch::dsigmadzz],
                                   cg->fgfs[fngfs + ShellPatch::dRdxx], 
                                   cg->fgfs[fngfs + ShellPatch::dRdxy], 
                                   cg->fgfs[fngfs + ShellPatch::dRdxz],
                                   cg->fgfs[fngfs + ShellPatch::dRdyy], 
                                   cg->fgfs[fngfs + ShellPatch::dRdyz], 
                                   cg->fgfs[fngfs + ShellPatch::dRdzz],
                                   sPp->data->bbox[0], sPp->data->bbox[1], sPp->data->bbox[2], 
                                   sPp->data->bbox[3], sPp->data->bbox[4], sPp->data->bbox[5],
                                   cg->fgfs[TZ0->sgfn], cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                                   cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                                   cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                                   cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                                   cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                                   cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                                   cg->fgfs[Lap0->sgfn], 
                                   cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                                   cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn],
                                   cg->fgfs[TZ_rhs->sgfn], cg->fgfs[phi_rhs->sgfn], cg->fgfs[trK_rhs->sgfn],
                                   cg->fgfs[gxx_rhs->sgfn], cg->fgfs[gxy_rhs->sgfn], cg->fgfs[gxz_rhs->sgfn],
                                   cg->fgfs[gyy_rhs->sgfn], cg->fgfs[gyz_rhs->sgfn], cg->fgfs[gzz_rhs->sgfn],
                                   cg->fgfs[Axx_rhs->sgfn], cg->fgfs[Axy_rhs->sgfn], cg->fgfs[Axz_rhs->sgfn],
                                   cg->fgfs[Ayy_rhs->sgfn], cg->fgfs[Ayz_rhs->sgfn], cg->fgfs[Azz_rhs->sgfn],
                                   cg->fgfs[Gmx_rhs->sgfn], cg->fgfs[Gmy_rhs->sgfn], cg->fgfs[Gmz_rhs->sgfn],
                                   cg->fgfs[Lap_rhs->sgfn], 
                                   cg->fgfs[Sfx_rhs->sgfn], cg->fgfs[Sfy_rhs->sgfn], cg->fgfs[Sfz_rhs->sgfn],
                                   cg->fgfs[dtSfx_rhs->sgfn], cg->fgfs[dtSfy_rhs->sgfn], cg->fgfs[dtSfz_rhs->sgfn],
                                   cg->fgfs[Gamxxx->sgfn], cg->fgfs[Gamxxy->sgfn], cg->fgfs[Gamxxz->sgfn],
                                   cg->fgfs[Gamxyy->sgfn], cg->fgfs[Gamxyz->sgfn], cg->fgfs[Gamxzz->sgfn],
                                   cg->fgfs[Gamyxx->sgfn], cg->fgfs[Gamyxy->sgfn], cg->fgfs[Gamyxz->sgfn],
                                   cg->fgfs[Gamyyy->sgfn], cg->fgfs[Gamyyz->sgfn], cg->fgfs[Gamyzz->sgfn],
                                   cg->fgfs[Gamzxx->sgfn], cg->fgfs[Gamzxy->sgfn], cg->fgfs[Gamzxz->sgfn],
                                   cg->fgfs[Gamzyy->sgfn], cg->fgfs[Gamzyz->sgfn], cg->fgfs[Gamzzz->sgfn],
                                   cg->fgfs[Rxx->sgfn], cg->fgfs[Rxy->sgfn], cg->fgfs[Rxz->sgfn], 
                                   cg->fgfs[Ryy->sgfn], cg->fgfs[Ryz->sgfn], cg->fgfs[Rzz->sgfn],
                                   cg->fgfs[Cons_Gx->sgfn], cg->fgfs[Cons_Gy->sgfn], cg->fgfs[Cons_Gz->sgfn],
#if (EXTO == 0)
                                   Symmetry, fbeps, sPp->data->sst);
            // extroplate rhs
            f_david_milton_extroplate_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                         cg->fgfs[TZ_rhs->sgfn], cg->fgfs[phi_rhs->sgfn], cg->fgfs[trK_rhs->sgfn],
                                         cg->fgfs[gxx_rhs->sgfn], cg->fgfs[gxy_rhs->sgfn], cg->fgfs[gxz_rhs->sgfn],
                                         cg->fgfs[gyy_rhs->sgfn], cg->fgfs[gyz_rhs->sgfn], cg->fgfs[gzz_rhs->sgfn],
                                         cg->fgfs[Axx_rhs->sgfn], cg->fgfs[Axy_rhs->sgfn], cg->fgfs[Axz_rhs->sgfn],
                                         cg->fgfs[Ayy_rhs->sgfn], cg->fgfs[Ayz_rhs->sgfn], cg->fgfs[Azz_rhs->sgfn],
                                         cg->fgfs[Gmx_rhs->sgfn], cg->fgfs[Gmy_rhs->sgfn], cg->fgfs[Gmz_rhs->sgfn],
                                         cg->fgfs[Lap_rhs->sgfn], 
                                         cg->fgfs[Sfx_rhs->sgfn], cg->fgfs[Sfy_rhs->sgfn], cg->fgfs[Sfz_rhs->sgfn],
                                         cg->fgfs[dtSfx_rhs->sgfn], cg->fgfs[dtSfy_rhs->sgfn], cg->fgfs[dtSfz_rhs->sgfn], 
                                         sPp->data->bbox[2], sPp->data->bbox[5]);

            MyList<var> *varl0 = StateList, *varl = SynchList_pre, *varlrhs = RHSList; 
            // we do not check the correspondence here
            
            while (varl0)
            {
              f_kodis_sh(cg->shape, cg->X[0], cg->X[1], cg->X[2], cg->fgfs[varl0->data->sgfn], cg->fgfs[varlrhs->data->sgfn],
                         varl0->data->SoA, Symmetry, numepsh, sPp->data->sst);
#elif (EXTO == 1 || EXTO == 2)
                                   Symmetry, numepsh, sPp->data->sst);

            MyList<var> *varl0 = StateList, *varl = SynchList_pre, *varlrhs = RHSList; 
            // we do not check the correspondence here
            
            while (varl0)
            {
#endif
              f_rungekutta4_rout(cg->shape, dT_lev, 
                                 cg->fgfs[varl0->data->sgfn], 
                                 cg->fgfs[varl->data->sgfn], 
                                 cg->fgfs[varlrhs->data->sgfn],
                                 iter_count);

              varl0 = varl0->next;
              varl = varl->next;
              varlrhs = varlrhs->next;
            }
#if (EXTO == 1)
            // extroplate variable itself
            f_david_milton_extroplate_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                         cg->fgfs[TZ->sgfn], cg->fgfs[phi->sgfn], cg->fgfs[trK->sgfn],
                                         cg->fgfs[gxx->sgfn], cg->fgfs[gxy->sgfn], cg->fgfs[gxz->sgfn],
                                         cg->fgfs[gyy->sgfn], cg->fgfs[gyz->sgfn], cg->fgfs[gzz->sgfn],
                                         cg->fgfs[Axx->sgfn], cg->fgfs[Axy->sgfn], cg->fgfs[Axz->sgfn],
                                         cg->fgfs[Ayy->sgfn], cg->fgfs[Ayz->sgfn], cg->fgfs[Azz->sgfn],
                                         cg->fgfs[Gmx->sgfn], cg->fgfs[Gmy->sgfn], cg->fgfs[Gmz->sgfn],
                                         cg->fgfs[Lap->sgfn], 
                                         cg->fgfs[Sfx->sgfn], cg->fgfs[Sfy->sgfn], cg->fgfs[Sfz->sgfn],
                                         cg->fgfs[dtSfx->sgfn], cg->fgfs[dtSfy->sgfn], cg->fgfs[dtSfz->sgfn], 
                                         sPp->data->bbox[2], sPp->data->bbox[5]);
#endif
          }
#if (chinot == 0)
          f_lowerboundset(cg->shape, cg->fgfs[phi->sgfn], chitiny);
#endif
        }
        if (BP == sPp->data->ble)
          break;
        BP = BP->next;
      }
      sPp = sPp->next;
    }
#if 0 
// check StateList
  {
         SH->Dump_Data(StateList,0,PhysTime,dT_lev);
         if(myrank == 0)
	 {
            cout<<"check StateList"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	 }
  }
#endif
#if 0
// check rhs  
  {
         SH->Dump_Data(RHSList,0,PhysTime,dT_lev);
         if(myrank == 0)
	 {
            cout<<"check rhs"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	 }
  }
#endif
#if 0 
// check SynchList_pre
  {
         SH->Dump_Data(SynchList_pre,0,PhysTime,dT_lev);
         if(myrank == 0)
	 {
            cout<<"check SynchList_pre"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	 }
  }
#endif
  }
  // check error information
  {
    int erh = ERROR;
    MPI_Allreduce(&erh, &ERROR, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  }
  if (ERROR)
  {
    SH->Dump_Data(StateList, 0, PhysTime, dT_lev);
    if (myrank == 0)
    {
      if (ErrorMonitor->outfile)
        ErrorMonitor->outfile << "find NaN in state variables on Shell Patches at t = " << PhysTime << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  Parallel::Sync(GH->PatL[lev], SynchList_pre, Symmetry);

  if (lev == 0)
  {
    clock_t prev_clock, curr_clock;
    if (myrank == 0)
      curr_clock = clock();
    SH->Synch(SynchList_pre, Symmetry);
    if (myrank == 0)
    {
      prev_clock = curr_clock;
      curr_clock = clock();
      cout << " Shell stuff synchronization used " 
           << (double)(curr_clock - prev_clock) / ((double)CLOCKS_PER_SEC) 
           << " seconds! " << endl;
    }

#ifdef SMOOTHSHELL
    // smooth Shell Patches
    if (lev == 0)
    {
      sPp = SH->PatL;
      while (sPp)
      {
        MyList<Block> *BP = sPp->data->blb;
        while (BP)
        {
          Block *cg = BP->data;
          if (myrank == cg->rank)
          {
            MyList<var> *varl = SynchList_pre;
            while (varl)
            {
              f_kodis_shcr(cg->shape, cg->X[0], cg->X[1], cg->X[2], 
                           cg->fgfs[varl->data->sgfn], cg->fgfs[varl->data->sgfn],
                           varl->data->SoA, Symmetry, numepsh, sPp->data->sst);
              varl = varl->next;
            }
          }
          if (BP == sPp->data->ble)
            break;
          BP = BP->next;
        }
        sPp = sPp->next;
      }
      SH->Synch(SynchList_pre, Symmetry);
    }
// end smooth
#endif

#if 0
// check SynchList_pre after Synch
  {
         SH->Dump_Data(SynchList_pre,0,PhysTime,dT_lev);
         if(myrank == 0)
	 {
            cout<<"check SynchList_pre"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	 }
  }
#endif
  }

  // for black hole position
  if (BH_num > 0 && lev == GH->levels - 1)
  {
    compute_Porg_rhs(Porg0, Porg_rhs, Sfx0, Sfy0, Sfz0, lev);
    for (int ithBH = 0; ithBH < BH_num; ithBH++)
    {
      f_rungekutta4_scalar(dT_lev, Porg0[ithBH][0], Porg[ithBH][0], Porg_rhs[ithBH][0], iter_count);
      f_rungekutta4_scalar(dT_lev, Porg0[ithBH][1], Porg[ithBH][1], Porg_rhs[ithBH][1], iter_count);
      f_rungekutta4_scalar(dT_lev, Porg0[ithBH][2], Porg[ithBH][2], Porg_rhs[ithBH][2], iter_count);
      if (Symmetry > 0)
        Porg[ithBH][2] = fabs(Porg[ithBH][2]);
      if (Symmetry == 2)
      {
        Porg[ithBH][0] = fabs(Porg[ithBH][0]);
        Porg[ithBH][1] = fabs(Porg[ithBH][1]);
      }
      if (!finite(Porg[ithBH][0]) || !finite(Porg[ithBH][1]) || !finite(Porg[ithBH][2]))
      {
        if (ErrorMonitor->outfile)
          ErrorMonitor->outfile << "predictor step finds NaN for BH's position from ("
                                << Porg0[ithBH][0] << "," << Porg0[ithBH][1] << "," << Porg0[ithBH][2] 
                                << ")" << endl;

        MyList<var> *DG_List = new MyList<var>(Sfx0);
        DG_List->insert(Sfx0);
        DG_List->insert(Sfy0);
        DG_List->insert(Sfz0);
        Parallel::Dump_Data(GH->PatL[lev], DG_List, 0, PhysTime, dT_lev);
        DG_List->clearList();
      }
    }
  }
  // data analysis part
  // Warning NOTE: the variables1 are used as temp storege room
  if (lev == a_lev)
  {
    AnalysisStuff(lev, dT_lev);
  }
  // corrector
  for (iter_count = 1; iter_count < 4; iter_count++)
  {
    // for RK4: t0, t0+dt/2, t0+dt/2, t0+dt;
    if (iter_count == 1 || iter_count == 3)
      TRK4 += dT_lev / 2;
    Pp = GH->PatL[lev];
    while (Pp)
    {
      MyList<Block> *BP = Pp->data->blb;
      while (BP)
      {
        Block *cg = BP->data;
        if (myrank == cg->rank)
        {
#if (AGM == 0)
          f_enforce_ga(cg->shape,
                       cg->fgfs[gxx->sgfn], cg->fgfs[gxy->sgfn], cg->fgfs[gxz->sgfn], 
                       cg->fgfs[gyy->sgfn], cg->fgfs[gyz->sgfn], cg->fgfs[gzz->sgfn],
                       cg->fgfs[Axx->sgfn], cg->fgfs[Axy->sgfn], cg->fgfs[Axz->sgfn], 
                       cg->fgfs[Ayy->sgfn], cg->fgfs[Ayz->sgfn], cg->fgfs[Azz->sgfn]);
#elif (AGM == 1)
          if (iter_count == 3)
            f_enforce_ga(cg->shape,
                         cg->fgfs[gxx->sgfn], cg->fgfs[gxy->sgfn], cg->fgfs[gxz->sgfn], 
                         cg->fgfs[gyy->sgfn], cg->fgfs[gyz->sgfn], cg->fgfs[gzz->sgfn],
                         cg->fgfs[Axx->sgfn], cg->fgfs[Axy->sgfn], cg->fgfs[Axz->sgfn], 
                         cg->fgfs[Ayy->sgfn], cg->fgfs[Ayz->sgfn], cg->fgfs[Azz->sgfn]);
#endif

#if (chinot == 0)
          if (f_compute_rhs_Z4c(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                                cg->fgfs[phi->sgfn], cg->fgfs[trK->sgfn],
                                cg->fgfs[gxx->sgfn], cg->fgfs[gxy->sgfn], cg->fgfs[gxz->sgfn], 
                                cg->fgfs[gyy->sgfn], cg->fgfs[gyz->sgfn], cg->fgfs[gzz->sgfn],
                                cg->fgfs[Axx->sgfn], cg->fgfs[Axy->sgfn], cg->fgfs[Axz->sgfn], 
                                cg->fgfs[Ayy->sgfn], cg->fgfs[Ayz->sgfn], cg->fgfs[Azz->sgfn],
                                cg->fgfs[Gmx->sgfn], cg->fgfs[Gmy->sgfn], cg->fgfs[Gmz->sgfn],
                                cg->fgfs[Lap->sgfn], 
                                cg->fgfs[Sfx->sgfn], cg->fgfs[Sfy->sgfn], cg->fgfs[Sfz->sgfn],
                                cg->fgfs[dtSfx->sgfn], cg->fgfs[dtSfy->sgfn], cg->fgfs[dtSfz->sgfn],
                                cg->fgfs[TZ->sgfn],
                                cg->fgfs[phi1->sgfn], cg->fgfs[trK1->sgfn],
                                cg->fgfs[gxx1->sgfn], cg->fgfs[gxy1->sgfn], cg->fgfs[gxz1->sgfn],
                                cg->fgfs[gyy1->sgfn], cg->fgfs[gyz1->sgfn], cg->fgfs[gzz1->sgfn],
                                cg->fgfs[Axx1->sgfn], cg->fgfs[Axy1->sgfn], cg->fgfs[Axz1->sgfn],
                                cg->fgfs[Ayy1->sgfn], cg->fgfs[Ayz1->sgfn], cg->fgfs[Azz1->sgfn],
                                cg->fgfs[Gmx1->sgfn], cg->fgfs[Gmy1->sgfn], cg->fgfs[Gmz1->sgfn],
                                cg->fgfs[Lap1->sgfn], 
                                cg->fgfs[Sfx1->sgfn], cg->fgfs[Sfy1->sgfn], cg->fgfs[Sfz1->sgfn],
                                cg->fgfs[dtSfx1->sgfn], cg->fgfs[dtSfy1->sgfn], cg->fgfs[dtSfz1->sgfn],
                                cg->fgfs[TZ1->sgfn],
                                cg->fgfs[rho->sgfn], 
                                cg->fgfs[Sx->sgfn], cg->fgfs[Sy->sgfn], cg->fgfs[Sz->sgfn],
                                cg->fgfs[Sxx->sgfn], cg->fgfs[Sxy->sgfn], cg->fgfs[Sxz->sgfn], 
                                cg->fgfs[Syy->sgfn], cg->fgfs[Syz->sgfn], cg->fgfs[Szz->sgfn],
                                cg->fgfs[Gamxxx->sgfn], cg->fgfs[Gamxxy->sgfn], cg->fgfs[Gamxxz->sgfn],
                                cg->fgfs[Gamxyy->sgfn], cg->fgfs[Gamxyz->sgfn], cg->fgfs[Gamxzz->sgfn],
                                cg->fgfs[Gamyxx->sgfn], cg->fgfs[Gamyxy->sgfn], cg->fgfs[Gamyxz->sgfn],
                                cg->fgfs[Gamyyy->sgfn], cg->fgfs[Gamyyz->sgfn], cg->fgfs[Gamyzz->sgfn],
                                cg->fgfs[Gamzxx->sgfn], cg->fgfs[Gamzxy->sgfn], cg->fgfs[Gamzxz->sgfn],
                                cg->fgfs[Gamzyy->sgfn], cg->fgfs[Gamzyz->sgfn], cg->fgfs[Gamzzz->sgfn],
                                cg->fgfs[Rxx->sgfn], cg->fgfs[Rxy->sgfn], cg->fgfs[Rxz->sgfn], 
                                cg->fgfs[Ryy->sgfn], cg->fgfs[Ryz->sgfn], cg->fgfs[Rzz->sgfn],
                                cg->fgfs[Cons_Ham->sgfn],
                                cg->fgfs[Cons_Px->sgfn], cg->fgfs[Cons_Py->sgfn], cg->fgfs[Cons_Pz->sgfn],
                                cg->fgfs[Cons_Gx->sgfn], cg->fgfs[Cons_Gy->sgfn], cg->fgfs[Cons_Gz->sgfn],
                                Symmetry, lev, ndeps, cor))
          {
            cout << "find NaN in domain: (" 
                 << cg->bbox[0] << ":" << cg->bbox[3] << "," 
                 << cg->bbox[1] << ":" << cg->bbox[4] << ","
                 << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
            ERROR = 1;
          }
#else
          if (f_compute_rhs_Z4cnot(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                                   cg->fgfs[phi->sgfn], cg->fgfs[trK->sgfn],
                                   cg->fgfs[gxx->sgfn], cg->fgfs[gxy->sgfn], cg->fgfs[gxz->sgfn], 
                                   cg->fgfs[gyy->sgfn], cg->fgfs[gyz->sgfn], cg->fgfs[gzz->sgfn],
                                   cg->fgfs[Axx->sgfn], cg->fgfs[Axy->sgfn], cg->fgfs[Axz->sgfn], 
                                   cg->fgfs[Ayy->sgfn], cg->fgfs[Ayz->sgfn], cg->fgfs[Azz->sgfn],
                                   cg->fgfs[Gmx->sgfn], cg->fgfs[Gmy->sgfn], cg->fgfs[Gmz->sgfn],
                                   cg->fgfs[Lap->sgfn], 
                                   cg->fgfs[Sfx->sgfn], cg->fgfs[Sfy->sgfn], cg->fgfs[Sfz->sgfn],
                                   cg->fgfs[dtSfx->sgfn], cg->fgfs[dtSfy->sgfn], cg->fgfs[dtSfz->sgfn],
                                   cg->fgfs[TZ->sgfn],
                                   cg->fgfs[phi1->sgfn], cg->fgfs[trK1->sgfn],
                                   cg->fgfs[gxx1->sgfn], cg->fgfs[gxy1->sgfn], cg->fgfs[gxz1->sgfn],
                                   cg->fgfs[gyy1->sgfn], cg->fgfs[gyz1->sgfn], cg->fgfs[gzz1->sgfn],
                                   cg->fgfs[Axx1->sgfn], cg->fgfs[Axy1->sgfn], cg->fgfs[Axz1->sgfn],
                                   cg->fgfs[Ayy1->sgfn], cg->fgfs[Ayz1->sgfn], cg->fgfs[Azz1->sgfn],
                                   cg->fgfs[Gmx1->sgfn], cg->fgfs[Gmy1->sgfn], cg->fgfs[Gmz1->sgfn],
                                   cg->fgfs[Lap1->sgfn], 
                                   cg->fgfs[Sfx1->sgfn], cg->fgfs[Sfy1->sgfn], cg->fgfs[Sfz1->sgfn],
                                   cg->fgfs[dtSfx1->sgfn], cg->fgfs[dtSfy1->sgfn], cg->fgfs[dtSfz1->sgfn],
                                   cg->fgfs[TZ1->sgfn],
                                   cg->fgfs[rho->sgfn], 
                                   cg->fgfs[Sx->sgfn], cg->fgfs[Sy->sgfn], cg->fgfs[Sz->sgfn],
                                   cg->fgfs[Sxx->sgfn], cg->fgfs[Sxy->sgfn], cg->fgfs[Sxz->sgfn], 
                                   cg->fgfs[Syy->sgfn], cg->fgfs[Syz->sgfn], cg->fgfs[Szz->sgfn],
                                   cg->fgfs[Gamxxx->sgfn], cg->fgfs[Gamxxy->sgfn], cg->fgfs[Gamxxz->sgfn],
                                   cg->fgfs[Gamxyy->sgfn], cg->fgfs[Gamxyz->sgfn], cg->fgfs[Gamxzz->sgfn],
                                   cg->fgfs[Gamyxx->sgfn], cg->fgfs[Gamyxy->sgfn], cg->fgfs[Gamyxz->sgfn],
                                   cg->fgfs[Gamyyy->sgfn], cg->fgfs[Gamyyz->sgfn], cg->fgfs[Gamyzz->sgfn],
                                   cg->fgfs[Gamzxx->sgfn], cg->fgfs[Gamzxy->sgfn], cg->fgfs[Gamzxz->sgfn],
                                   cg->fgfs[Gamzyy->sgfn], cg->fgfs[Gamzyz->sgfn], cg->fgfs[Gamzzz->sgfn],
                                   cg->fgfs[Rxx->sgfn], cg->fgfs[Rxy->sgfn], cg->fgfs[Rxz->sgfn], 
                                   cg->fgfs[Ryy->sgfn], cg->fgfs[Ryz->sgfn], cg->fgfs[Rzz->sgfn],
                                   cg->fgfs[Cons_Ham->sgfn],
                                   cg->fgfs[Cons_Px->sgfn], cg->fgfs[Cons_Py->sgfn], cg->fgfs[Cons_Pz->sgfn],
                                   cg->fgfs[Cons_Gx->sgfn], cg->fgfs[Cons_Gy->sgfn], cg->fgfs[Cons_Gz->sgfn],
                                   Symmetry, lev, ndeps, cor, chitiny))
          {
            cout << "find NaN in domain: (" 
                 << cg->bbox[0] << ":" << cg->bbox[3] << "," 
                 << cg->bbox[1] << ":" << cg->bbox[4] << ","
                 << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
            ERROR = 1;
          }
#endif
          // rk4 substep and boundary
          {
            MyList<var> *varl0 = StateList, *varl = SynchList_pre, *varl1 = SynchList_cor, *varlrhs = RHSList; 
            // we do not check the correspondence here
            
            while (varl0)
            {
#if (MRBD == 1)
              // sommerfeld indeed
              f_sommerfeld_routbam(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                   Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], 
                                   Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                                   cg->fgfs[varl1->data->sgfn],
                                   cg->fgfs[varl->data->sgfn], 
                                   varl0->data->propspeed, varl0->data->SoA,
                                   Symmetry);
#endif
              f_rungekutta4_rout(cg->shape, dT_lev, 
                                 cg->fgfs[varl0->data->sgfn], 
                                 cg->fgfs[varl1->data->sgfn], 
                                 cg->fgfs[varlrhs->data->sgfn],
                                 iter_count);
#if (MRBD == 0)
              f_sommerfeld_rout(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], 
                                Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                                dT_lev, 
                                cg->fgfs[phi0->sgfn],
                                cg->fgfs[Lap0->sgfn], 
                                cg->fgfs[varl0->data->sgfn], cg->fgfs[varl1->data->sgfn], 
                                varl0->data->SoA,
                                Symmetry, cor);
#elif (MRBD == 2)
              f_sommerfeld_rout(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], 
                                Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                                dT_lev, 
                                cg->fgfs[phi0->sgfn],
                                cg->fgfs[Lap0->sgfn], 
                                cg->fgfs[varl->data->sgfn], cg->fgfs[varl1->data->sgfn], 
                                varl0->data->SoA,
                                Symmetry, cor);
#endif

              varl0 = varl0->next;
              varl = varl->next;
              varl1 = varl1->next;
              varlrhs = varlrhs->next;
            }
          }
#if (chinot == 0)
          f_lowerboundset(cg->shape, cg->fgfs[phi1->sgfn], chitiny);
#endif
        }
        if (BP == Pp->data->ble)
          break;
        BP = BP->next;
      }
      Pp = Pp->next;
    }

    // check error information
    {
      int erh = ERROR;
      MPI_Allreduce(&erh, &ERROR, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    }
    if (ERROR)
    {
      Parallel::Dump_Data(GH->PatL[lev], SynchList_pre, 0, PhysTime, dT_lev);
      if (myrank == 0)
      {
        if (ErrorMonitor->outfile)
          ErrorMonitor->outfile << "find NaN in RK4 substep#" << iter_count 
                                << " variables at t = " << PhysTime 
                                << ", lev = " << lev << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }

    // evolve Shell Patches
    if (lev == 0)
    {
      sPp = SH->PatL;
      while (sPp)
      {
        MyList<Block> *BP = sPp->data->blb;
        int fngfs = sPp->data->fngfs;
        while (BP)
        {
          Block *cg = BP->data;
          if (myrank == cg->rank)
          {
#if (EXTO == 2)
            // extroplate variable itself
            f_david_milton_extroplate_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                         cg->fgfs[TZ->sgfn], cg->fgfs[phi->sgfn], cg->fgfs[trK->sgfn],
                                         cg->fgfs[gxx->sgfn], cg->fgfs[gxy->sgfn], cg->fgfs[gxz->sgfn],
                                         cg->fgfs[gyy->sgfn], cg->fgfs[gyz->sgfn], cg->fgfs[gzz->sgfn],
                                         cg->fgfs[Axx->sgfn], cg->fgfs[Axy->sgfn], cg->fgfs[Axz->sgfn],
                                         cg->fgfs[Ayy->sgfn], cg->fgfs[Ayz->sgfn], cg->fgfs[Azz->sgfn],
                                         cg->fgfs[Gmx->sgfn], cg->fgfs[Gmy->sgfn], cg->fgfs[Gmz->sgfn],
                                         cg->fgfs[Lap->sgfn], 
                                         cg->fgfs[Sfx->sgfn], cg->fgfs[Sfy->sgfn], cg->fgfs[Sfz->sgfn],
                                         cg->fgfs[dtSfx->sgfn], cg->fgfs[dtSfy->sgfn], cg->fgfs[dtSfz->sgfn], 
                                         sPp->data->bbox[2], sPp->data->bbox[5]);
#endif

#if (AGM == 0)
            f_enforce_ga(cg->shape,
                         cg->fgfs[gxx->sgfn], cg->fgfs[gxy->sgfn], cg->fgfs[gxz->sgfn], 
                         cg->fgfs[gyy->sgfn], cg->fgfs[gyz->sgfn], cg->fgfs[gzz->sgfn],
                         cg->fgfs[Axx->sgfn], cg->fgfs[Axy->sgfn], cg->fgfs[Axz->sgfn], 
                         cg->fgfs[Ayy->sgfn], cg->fgfs[Ayz->sgfn], cg->fgfs[Azz->sgfn]);
#elif (AGM == 1)
            if (iter_count == 3)
              f_enforce_ga(cg->shape,
                           cg->fgfs[gxx->sgfn], cg->fgfs[gxy->sgfn], cg->fgfs[gxz->sgfn], 
                           cg->fgfs[gyy->sgfn], cg->fgfs[gyz->sgfn], cg->fgfs[gzz->sgfn],
                           cg->fgfs[Axx->sgfn], cg->fgfs[Axy->sgfn], cg->fgfs[Axz->sgfn], 
                           cg->fgfs[Ayy->sgfn], cg->fgfs[Ayz->sgfn], cg->fgfs[Azz->sgfn]);
#endif

            if (f_compute_rhs_Z4c_ss(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                                     cg->fgfs[fngfs + ShellPatch::gx], 
                                     cg->fgfs[fngfs + ShellPatch::gy], 
                                     cg->fgfs[fngfs + ShellPatch::gz],
                                     cg->fgfs[fngfs + ShellPatch::drhodx], 
                                     cg->fgfs[fngfs + ShellPatch::drhody], 
                                     cg->fgfs[fngfs + ShellPatch::drhodz],
                                     cg->fgfs[fngfs + ShellPatch::dsigmadx], 
                                     cg->fgfs[fngfs + ShellPatch::dsigmady], 
                                     cg->fgfs[fngfs + ShellPatch::dsigmadz],
                                     cg->fgfs[fngfs + ShellPatch::dRdx], 
                                     cg->fgfs[fngfs + ShellPatch::dRdy], 
                                     cg->fgfs[fngfs + ShellPatch::dRdz],
                                     cg->fgfs[fngfs + ShellPatch::drhodxx], 
                                     cg->fgfs[fngfs + ShellPatch::drhodxy], 
                                     cg->fgfs[fngfs + ShellPatch::drhodxz],
                                     cg->fgfs[fngfs + ShellPatch::drhodyy], 
                                     cg->fgfs[fngfs + ShellPatch::drhodyz], 
                                     cg->fgfs[fngfs + ShellPatch::drhodzz],
                                     cg->fgfs[fngfs + ShellPatch::dsigmadxx], 
                                     cg->fgfs[fngfs + ShellPatch::dsigmadxy], 
                                     cg->fgfs[fngfs + ShellPatch::dsigmadxz],
                                     cg->fgfs[fngfs + ShellPatch::dsigmadyy], 
                                     cg->fgfs[fngfs + ShellPatch::dsigmadyz], 
                                     cg->fgfs[fngfs + ShellPatch::dsigmadzz],
                                     cg->fgfs[fngfs + ShellPatch::dRdxx], 
                                     cg->fgfs[fngfs + ShellPatch::dRdxy], 
                                     cg->fgfs[fngfs + ShellPatch::dRdxz],
                                     cg->fgfs[fngfs + ShellPatch::dRdyy], 
                                     cg->fgfs[fngfs + ShellPatch::dRdyz], 
                                     cg->fgfs[fngfs + ShellPatch::dRdzz],
                                     cg->fgfs[phi->sgfn], cg->fgfs[trK->sgfn],
                                     cg->fgfs[gxx->sgfn], cg->fgfs[gxy->sgfn], cg->fgfs[gxz->sgfn], 
                                     cg->fgfs[gyy->sgfn], cg->fgfs[gyz->sgfn], cg->fgfs[gzz->sgfn],
                                     cg->fgfs[Axx->sgfn], cg->fgfs[Axy->sgfn], cg->fgfs[Axz->sgfn], 
                                     cg->fgfs[Ayy->sgfn], cg->fgfs[Ayz->sgfn], cg->fgfs[Azz->sgfn],
                                     cg->fgfs[Gmx->sgfn], cg->fgfs[Gmy->sgfn], cg->fgfs[Gmz->sgfn],
                                     cg->fgfs[Lap->sgfn], 
                                     cg->fgfs[Sfx->sgfn], cg->fgfs[Sfy->sgfn], cg->fgfs[Sfz->sgfn],
                                     cg->fgfs[dtSfx->sgfn], cg->fgfs[dtSfy->sgfn], cg->fgfs[dtSfz->sgfn],
                                     cg->fgfs[TZ->sgfn],
                                     cg->fgfs[phi1->sgfn], cg->fgfs[trK1->sgfn],
                                     cg->fgfs[gxx1->sgfn], cg->fgfs[gxy1->sgfn], cg->fgfs[gxz1->sgfn],
                                     cg->fgfs[gyy1->sgfn], cg->fgfs[gyz1->sgfn], cg->fgfs[gzz1->sgfn],
                                     cg->fgfs[Axx1->sgfn], cg->fgfs[Axy1->sgfn], cg->fgfs[Axz1->sgfn],
                                     cg->fgfs[Ayy1->sgfn], cg->fgfs[Ayz1->sgfn], cg->fgfs[Azz1->sgfn],
                                     cg->fgfs[Gmx1->sgfn], cg->fgfs[Gmy1->sgfn], cg->fgfs[Gmz1->sgfn],
                                     cg->fgfs[Lap1->sgfn], 
                                     cg->fgfs[Sfx1->sgfn], cg->fgfs[Sfy1->sgfn], cg->fgfs[Sfz1->sgfn],
                                     cg->fgfs[dtSfx1->sgfn], cg->fgfs[dtSfy1->sgfn], cg->fgfs[dtSfz1->sgfn],
                                     cg->fgfs[TZ1->sgfn],
                                     cg->fgfs[rho->sgfn], 
                                     cg->fgfs[Sx->sgfn], cg->fgfs[Sy->sgfn], cg->fgfs[Sz->sgfn],
                                     cg->fgfs[Sxx->sgfn], cg->fgfs[Sxy->sgfn], cg->fgfs[Sxz->sgfn], 
                                     cg->fgfs[Syy->sgfn], cg->fgfs[Syz->sgfn], cg->fgfs[Szz->sgfn],
                                     cg->fgfs[Gamxxx->sgfn], cg->fgfs[Gamxxy->sgfn], cg->fgfs[Gamxxz->sgfn],
                                     cg->fgfs[Gamxyy->sgfn], cg->fgfs[Gamxyz->sgfn], cg->fgfs[Gamxzz->sgfn],
                                     cg->fgfs[Gamyxx->sgfn], cg->fgfs[Gamyxy->sgfn], cg->fgfs[Gamyxz->sgfn],
                                     cg->fgfs[Gamyyy->sgfn], cg->fgfs[Gamyyz->sgfn], cg->fgfs[Gamyzz->sgfn],
                                     cg->fgfs[Gamzxx->sgfn], cg->fgfs[Gamzxy->sgfn], cg->fgfs[Gamzxz->sgfn],
                                     cg->fgfs[Gamzyy->sgfn], cg->fgfs[Gamzyz->sgfn], cg->fgfs[Gamzzz->sgfn],
                                     cg->fgfs[Rxx->sgfn], cg->fgfs[Rxy->sgfn], cg->fgfs[Rxz->sgfn], 
                                     cg->fgfs[Ryy->sgfn], cg->fgfs[Ryz->sgfn], cg->fgfs[Rzz->sgfn],
                                     cg->fgfs[Cons_Ham->sgfn],
                                     cg->fgfs[Cons_Px->sgfn], cg->fgfs[Cons_Py->sgfn], cg->fgfs[Cons_Pz->sgfn],
                                     cg->fgfs[Cons_Gx->sgfn], cg->fgfs[Cons_Gy->sgfn], cg->fgfs[Cons_Gz->sgfn],
                                     Symmetry, lev, fbeps, sPp->data->sst, cor))
            {
              cout << "find NaN in Shell domain: sst = " << sPp->data->sst << ", (" 
                   << cg->bbox[0] << ":" << cg->bbox[3] << ","
                   << cg->bbox[1] << ":" << cg->bbox[4] << "," 
                   << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
              ERROR = 1;
            }
            // rk4 substep and boundary
            {
              // CPBC indeed for outter boudary while fix BD for inner boundary
              f_david_milton_cpbc_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                     cg->fgfs[fngfs + ShellPatch::gx], 
                                     cg->fgfs[fngfs + ShellPatch::gy], 
                                     cg->fgfs[fngfs + ShellPatch::gz],
                                     cg->fgfs[fngfs + ShellPatch::drhodx], 
                                     cg->fgfs[fngfs + ShellPatch::drhody], 
                                     cg->fgfs[fngfs + ShellPatch::drhodz],
                                     cg->fgfs[fngfs + ShellPatch::dsigmadx], 
                                     cg->fgfs[fngfs + ShellPatch::dsigmady], 
                                     cg->fgfs[fngfs + ShellPatch::dsigmadz],
                                     cg->fgfs[fngfs + ShellPatch::dRdx], 
                                     cg->fgfs[fngfs + ShellPatch::dRdy], 
                                     cg->fgfs[fngfs + ShellPatch::dRdz],
                                     cg->fgfs[fngfs + ShellPatch::drhodxx], 
                                     cg->fgfs[fngfs + ShellPatch::drhodxy], 
                                     cg->fgfs[fngfs + ShellPatch::drhodxz],
                                     cg->fgfs[fngfs + ShellPatch::drhodyy], 
                                     cg->fgfs[fngfs + ShellPatch::drhodyz], 
                                     cg->fgfs[fngfs + ShellPatch::drhodzz],
                                     cg->fgfs[fngfs + ShellPatch::dsigmadxx], 
                                     cg->fgfs[fngfs + ShellPatch::dsigmadxy], 
                                     cg->fgfs[fngfs + ShellPatch::dsigmadxz],
                                     cg->fgfs[fngfs + ShellPatch::dsigmadyy], 
                                     cg->fgfs[fngfs + ShellPatch::dsigmadyz], 
                                     cg->fgfs[fngfs + ShellPatch::dsigmadzz],
                                     cg->fgfs[fngfs + ShellPatch::dRdxx], 
                                     cg->fgfs[fngfs + ShellPatch::dRdxy], 
                                     cg->fgfs[fngfs + ShellPatch::dRdxz],
                                     cg->fgfs[fngfs + ShellPatch::dRdyy], 
                                     cg->fgfs[fngfs + ShellPatch::dRdyz], 
                                     cg->fgfs[fngfs + ShellPatch::dRdzz],
                                     sPp->data->bbox[0], sPp->data->bbox[1], sPp->data->bbox[2], 
                                     sPp->data->bbox[3], sPp->data->bbox[4], sPp->data->bbox[5],
                                     cg->fgfs[TZ->sgfn], cg->fgfs[phi->sgfn], cg->fgfs[trK->sgfn],
                                     cg->fgfs[gxx->sgfn], cg->fgfs[gxy->sgfn], cg->fgfs[gxz->sgfn], 
                                     cg->fgfs[gyy->sgfn], cg->fgfs[gyz->sgfn], cg->fgfs[gzz->sgfn],
                                     cg->fgfs[Axx->sgfn], cg->fgfs[Axy->sgfn], cg->fgfs[Axz->sgfn], 
                                     cg->fgfs[Ayy->sgfn], cg->fgfs[Ayz->sgfn], cg->fgfs[Azz->sgfn],
                                     cg->fgfs[Gmx->sgfn], cg->fgfs[Gmy->sgfn], cg->fgfs[Gmz->sgfn],
                                     cg->fgfs[Lap->sgfn], 
                                     cg->fgfs[Sfx->sgfn], cg->fgfs[Sfy->sgfn], cg->fgfs[Sfz->sgfn],
                                     cg->fgfs[dtSfx->sgfn], cg->fgfs[dtSfy->sgfn], cg->fgfs[dtSfz->sgfn],
                                     cg->fgfs[TZ1->sgfn], cg->fgfs[phi1->sgfn], cg->fgfs[trK1->sgfn],
                                     cg->fgfs[gxx1->sgfn], cg->fgfs[gxy1->sgfn], cg->fgfs[gxz1->sgfn],
                                     cg->fgfs[gyy1->sgfn], cg->fgfs[gyz1->sgfn], cg->fgfs[gzz1->sgfn],
                                     cg->fgfs[Axx1->sgfn], cg->fgfs[Axy1->sgfn], cg->fgfs[Axz1->sgfn],
                                     cg->fgfs[Ayy1->sgfn], cg->fgfs[Ayz1->sgfn], cg->fgfs[Azz1->sgfn],
                                     cg->fgfs[Gmx1->sgfn], cg->fgfs[Gmy1->sgfn], cg->fgfs[Gmz1->sgfn],
                                     cg->fgfs[Lap1->sgfn], 
                                     cg->fgfs[Sfx1->sgfn], cg->fgfs[Sfy1->sgfn], cg->fgfs[Sfz1->sgfn],
                                     cg->fgfs[dtSfx1->sgfn], cg->fgfs[dtSfy1->sgfn], cg->fgfs[dtSfz1->sgfn],
                                     cg->fgfs[Gamxxx->sgfn], cg->fgfs[Gamxxy->sgfn], cg->fgfs[Gamxxz->sgfn],
                                     cg->fgfs[Gamxyy->sgfn], cg->fgfs[Gamxyz->sgfn], cg->fgfs[Gamxzz->sgfn],
                                     cg->fgfs[Gamyxx->sgfn], cg->fgfs[Gamyxy->sgfn], cg->fgfs[Gamyxz->sgfn],
                                     cg->fgfs[Gamyyy->sgfn], cg->fgfs[Gamyyz->sgfn], cg->fgfs[Gamyzz->sgfn],
                                     cg->fgfs[Gamzxx->sgfn], cg->fgfs[Gamzxy->sgfn], cg->fgfs[Gamzxz->sgfn],
                                     cg->fgfs[Gamzyy->sgfn], cg->fgfs[Gamzyz->sgfn], cg->fgfs[Gamzzz->sgfn],
                                     cg->fgfs[Rxx->sgfn], cg->fgfs[Rxy->sgfn], cg->fgfs[Rxz->sgfn], 
                                     cg->fgfs[Ryy->sgfn], cg->fgfs[Ryz->sgfn], cg->fgfs[Rzz->sgfn],
                                     cg->fgfs[Cons_Gx->sgfn], cg->fgfs[Cons_Gy->sgfn], cg->fgfs[Cons_Gz->sgfn],
#if (EXTO == 0)
                                     Symmetry, fbeps, sPp->data->sst);
              // extroplate rhs
              f_david_milton_extroplate_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                           cg->fgfs[TZ1->sgfn], cg->fgfs[phi1->sgfn], cg->fgfs[trK1->sgfn],
                                           cg->fgfs[gxx1->sgfn], cg->fgfs[gxy1->sgfn], cg->fgfs[gxz1->sgfn],
                                           cg->fgfs[gyy1->sgfn], cg->fgfs[gyz1->sgfn], cg->fgfs[gzz1->sgfn],
                                           cg->fgfs[Axx1->sgfn], cg->fgfs[Axy1->sgfn], cg->fgfs[Axz1->sgfn],
                                           cg->fgfs[Ayy1->sgfn], cg->fgfs[Ayz1->sgfn], cg->fgfs[Azz1->sgfn],
                                           cg->fgfs[Gmx1->sgfn], cg->fgfs[Gmy1->sgfn], cg->fgfs[Gmz1->sgfn],
                                           cg->fgfs[Lap1->sgfn], 
                                           cg->fgfs[Sfx1->sgfn], cg->fgfs[Sfy1->sgfn], cg->fgfs[Sfz1->sgfn],
                                           cg->fgfs[dtSfx1->sgfn], cg->fgfs[dtSfy1->sgfn], cg->fgfs[dtSfz1->sgfn], 
                                           sPp->data->bbox[2], sPp->data->bbox[5]);

              MyList<var> *varl0 = StateList, *varl = SynchList_pre, *varl1 = SynchList_cor, *varlrhs = RHSList; 
              // we do not check the correspondence here
              
              while (varl0)
              {
                f_kodis_sh(cg->shape, cg->X[0], cg->X[1], cg->X[2], 
                           cg->fgfs[varl->data->sgfn], cg->fgfs[varl1->data->sgfn],
                           varl->data->SoA, Symmetry, numepsh, sPp->data->sst);
#elif (EXTO == 1 || EXTO == 2)
                           Symmetry, numepsh, sPp->data->sst);

              MyList<var> *varl0 = StateList, *varl = SynchList_pre, *varl1 = SynchList_cor, *varlrhs = RHSList; // we do not check the correspondence here
              while (varl0)
              {
#endif
                f_rungekutta4_rout(cg->shape, dT_lev, cg->fgfs[varl0->data->sgfn], cg->fgfs[varl1->data->sgfn], cg->fgfs[varlrhs->data->sgfn],
                                   iter_count);

                varl0 = varl0->next;
                varl = varl->next;
                varl1 = varl1->next;
                varlrhs = varlrhs->next;
              }
#if (EXTO == 1)
              // extroplate variable itself
              f_david_milton_extroplate_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                           cg->fgfs[TZ1->sgfn], cg->fgfs[phi1->sgfn], cg->fgfs[trK1->sgfn],
                                           cg->fgfs[gxx1->sgfn], cg->fgfs[gxy1->sgfn], cg->fgfs[gxz1->sgfn],
                                           cg->fgfs[gyy1->sgfn], cg->fgfs[gyz1->sgfn], cg->fgfs[gzz1->sgfn],
                                           cg->fgfs[Axx1->sgfn], cg->fgfs[Axy1->sgfn], cg->fgfs[Axz1->sgfn],
                                           cg->fgfs[Ayy1->sgfn], cg->fgfs[Ayz1->sgfn], cg->fgfs[Azz1->sgfn],
                                           cg->fgfs[Gmx1->sgfn], cg->fgfs[Gmy1->sgfn], cg->fgfs[Gmz1->sgfn],
                                           cg->fgfs[Lap1->sgfn], 
                                           cg->fgfs[Sfx1->sgfn], cg->fgfs[Sfy1->sgfn], cg->fgfs[Sfz1->sgfn],
                                           cg->fgfs[dtSfx1->sgfn], cg->fgfs[dtSfy1->sgfn], cg->fgfs[dtSfz1->sgfn], 
                                           sPp->data->bbox[2], sPp->data->bbox[5]);
#endif
            }
#if (chinot == 0)
            f_lowerboundset(cg->shape, cg->fgfs[phi1->sgfn], chitiny);
#endif
          }
          if (BP == sPp->data->ble)
            break;
          BP = BP->next;
        }
        sPp = sPp->next;
      }
    }
    // check error information
    {
      int erh = ERROR;
      MPI_Allreduce(&erh, &ERROR, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    }
    if (ERROR)
    {
      SH->Dump_Data(SynchList_pre, 0, PhysTime, dT_lev);
      if (myrank == 0)
      {
        if (ErrorMonitor->outfile)
          ErrorMonitor->outfile << "find NaN on Shell Patches in RK4 substep#" << iter_count 
                                << " variables at t = " << PhysTime << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }

    Parallel::Sync(GH->PatL[lev], SynchList_cor, Symmetry);

    if (lev == 0)
    {
      clock_t prev_clock, curr_clock;
      if (myrank == 0)
        curr_clock = clock();
      SH->Synch(SynchList_cor, Symmetry);
      if (myrank == 0)
      {
        prev_clock = curr_clock;
        curr_clock = clock();
        cout << " Shell stuff synchronization used " 
             << (double)(curr_clock - prev_clock) / ((double)CLOCKS_PER_SEC) 
             << " seconds! " << endl;
      }
    }

#ifdef SMOOTHSHELL
    // smooth Shell Patches
    if (lev == 0)
    {
      sPp = SH->PatL;
      while (sPp)
      {
        MyList<Block> *BP = sPp->data->blb;
        while (BP)
        {
          Block *cg = BP->data;
          if (myrank == cg->rank)
          {
            MyList<var> *varl = SynchList_cor;
            while (varl)
            {
              f_kodis_shcr(cg->shape, cg->X[0], cg->X[1], cg->X[2], 
                           cg->fgfs[varl->data->sgfn], cg->fgfs[varl->data->sgfn],
                           varl->data->SoA, Symmetry, numepsh, sPp->data->sst);
              varl = varl->next;
            }
          }
          if (BP == sPp->data->ble)
            break;
          BP = BP->next;
        }
        sPp = sPp->next;
      }
      SH->Synch(SynchList_cor, Symmetry);
    }
// end smooth
#endif

    // for black hole position
    if (BH_num > 0 && lev == GH->levels - 1)
    {
      compute_Porg_rhs(Porg, Porg1, Sfx, Sfy, Sfz, lev);
      for (int ithBH = 0; ithBH < BH_num; ithBH++)
      {
        f_rungekutta4_scalar(dT_lev, Porg0[ithBH][0], Porg1[ithBH][0], Porg_rhs[ithBH][0], iter_count);
        f_rungekutta4_scalar(dT_lev, Porg0[ithBH][1], Porg1[ithBH][1], Porg_rhs[ithBH][1], iter_count);
        f_rungekutta4_scalar(dT_lev, Porg0[ithBH][2], Porg1[ithBH][2], Porg_rhs[ithBH][2], iter_count);
        if (Symmetry > 0)
          Porg1[ithBH][2] = fabs(Porg1[ithBH][2]);
        if (Symmetry == 2)
        {
          Porg1[ithBH][0] = fabs(Porg1[ithBH][0]);
          Porg1[ithBH][1] = fabs(Porg1[ithBH][1]);
        }
        if (!finite(Porg1[ithBH][0]) || !finite(Porg1[ithBH][1]) || !finite(Porg1[ithBH][2]))
        {
          if (ErrorMonitor->outfile)
            ErrorMonitor->outfile << iter_count << " corrector step finds NaN for BH's position from ("
                                  << Porg[ithBH][0] << "," << Porg[ithBH][1] << "," << Porg[ithBH][2] 
                                  << ")" << endl;

          MyList<var> *DG_List = new MyList<var>(Sfx0);
          DG_List->insert(Sfx0);
          DG_List->insert(Sfy0);
          DG_List->insert(Sfz0);
          Parallel::Dump_Data(GH->PatL[lev], DG_List, 0, PhysTime, dT_lev);
          DG_List->clearList();
        }
      }
    }
    // swap time level
    if (iter_count < 3)
    {
      Pp = GH->PatL[lev];
      while (Pp)
      {
        MyList<Block> *BP = Pp->data->blb;
        while (BP)
        {
          Block *cg = BP->data;
          cg->swapList(SynchList_pre, SynchList_cor, myrank);
          if (BP == Pp->data->ble)
            break;
          BP = BP->next;
        }
        Pp = Pp->next;
      }

      if (lev == 0)
      {
        sPp = SH->PatL;
        while (sPp)
        {
          MyList<Block> *BP = sPp->data->blb;
          while (BP)
          {
            Block *cg = BP->data;
            cg->swapList(SynchList_pre, SynchList_cor, myrank);
            if (BP == sPp->data->ble)
              break;
            BP = BP->next;
          }
          sPp = sPp->next;
        }
      }
      // for black hole position
      if (BH_num > 0 && lev == GH->levels - 1)
      {
        for (int ithBH = 0; ithBH < BH_num; ithBH++)
        {
          Porg[ithBH][0] = Porg1[ithBH][0];
          Porg[ithBH][1] = Porg1[ithBH][1];
          Porg[ithBH][2] = Porg1[ithBH][2];
        }
      }
    }
  }

#if (RPS == 0)
  // mesh refinement boundary part
  RestrictProlong(lev, YN, BB);

  if (lev == 0)
  {
    clock_t prev_clock, curr_clock;
    if (myrank == 0)
      curr_clock = clock();
    SH->CS_Inter(SynchList_cor, Symmetry);
    if (myrank == 0)
    {
      prev_clock = curr_clock;
      curr_clock = clock();
      cout << " CS_Inter used " << (double)(curr_clock - prev_clock) / ((double)CLOCKS_PER_SEC) 
           << " seconds! " << endl;
    }
  }
#endif

  // note the data structure before update
  // SynchList_cor 1   -----------
  //
  // StateList     0   -----------
  //
  // OldStateList  old -----------
  // update
  Pp = GH->PatL[lev];
  while (Pp)
  {
    MyList<Block> *BP = Pp->data->blb;
    while (BP)
    {
      Block *cg = BP->data;
      cg->swapList(StateList, SynchList_cor, myrank);
      cg->swapList(OldStateList, SynchList_cor, myrank);
      if (BP == Pp->data->ble)
        break;
      BP = BP->next;
    }
    Pp = Pp->next;
  }

  if (lev == 0)
  {
    sPp = SH->PatL;
    while (sPp)
    {
      MyList<Block> *BP = sPp->data->blb;
      while (BP)
      {
        Block *cg = BP->data;
        cg->swapList(StateList, SynchList_cor, myrank);
        cg->swapList(OldStateList, SynchList_cor, myrank);
        if (BP == sPp->data->ble)
          break;
        BP = BP->next;
      }
      sPp = sPp->next;
    }
#if 0
// check StateList   
  {
         SH->Dump_Data(StateList,0,PhysTime,dT_lev);
         if(myrank == 0)
	 {
            cout<<"check StateList"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	 }
  }
#endif
  }

  // for black hole position
  if (BH_num > 0 && lev == GH->levels - 1)
  {
    for (int ithBH = 0; ithBH < BH_num; ithBH++)
    {
      Porg0[ithBH][0] = Porg1[ithBH][0];
      Porg0[ithBH][1] = Porg1[ithBH][1];
      Porg0[ithBH][2] = Porg1[ithBH][2];
    }
  }
#if 0
     if(lev>6)
	{     
	  char str[50];	
	  MyList<var> * DG_List=new MyList<var>(Cons_Ham);
       	  DG_List->insert(Cons_Px); DG_List->insert(Cons_Py); DG_List->insert(Cons_Px);
       	  DG_List->insert(Cons_Gx); DG_List->insert(Cons_Gy); DG_List->insert(Cons_Gx);
	  printf(str,"lao%d",lev);
	  Parallel::Dump_Data(GH->PatL[6],DG_List,str,PhysTime,dT_lev);
	  DG_List->clearList();
	}
#endif
}
#endif
#undef MRBD

//================================================================================================



//================================================================================================

// 该成员函数用于检查插值结果？

//================================================================================================

void Z4c_class::Check_extrop()
{
  MyList<ss_patch> *sPp;

  sPp = SH->PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        f_david_milton_extroplate_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                     cg->fgfs[TZ0->sgfn], cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                                     cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn],
                                     cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                                     cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn],
                                     cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                                     cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                                     cg->fgfs[Lap0->sgfn], 
                                     cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                                     cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn], 
                                     sPp->data->bbox[2], sPp->data->bbox[5]);
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }

  SH->Dump_Data(StateList, "extrop", 0, 1);
  if (myrank == 0)
    MPI_Abort(MPI_COMM_WORLD, 1);
}

//================================================================================================



//================================================================================================

// 该成员函数用于计算和输出约束违反

//================================================================================================

void Z4c_class::Constraint_Out()
{
  // 这里要跟父类中用同样的变量
  // 否则不会传入正确的时间
  LastConsOut += dT * pow(0.5, Mymax(0, trfls));
  
  if (LastConsOut >= AnasTime)
  // Constraint violation
  {
    // recompute least the constraint data lost for moved new grid
    for (int lev = 0; lev < GH->levels; lev++)
    {
      // make sure the data consistent for higher levels
      if (lev > 0)
      {
        double TRK4 = PhysTime;
        double ndeps = numepsb;
        int pre = 0;
        MyList<Patch> *Pp = GH->PatL[lev];
        while (Pp)
        {
          MyList<Block> *BP = Pp->data->blb;
          while (BP)
          {
            Block *cg = BP->data;
            if (myrank == cg->rank)
            {
              f_compute_rhs_Z4c(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                                cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                                cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                                cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                                cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                                cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                                cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                                cg->fgfs[Lap0->sgfn], 
                                cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                                cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn],
                                cg->fgfs[TZ0->sgfn],
                                cg->fgfs[phi_rhs->sgfn], cg->fgfs[trK_rhs->sgfn],
                                cg->fgfs[gxx_rhs->sgfn], cg->fgfs[gxy_rhs->sgfn], cg->fgfs[gxz_rhs->sgfn],
                                cg->fgfs[gyy_rhs->sgfn], cg->fgfs[gyz_rhs->sgfn], cg->fgfs[gzz_rhs->sgfn],
                                cg->fgfs[Axx_rhs->sgfn], cg->fgfs[Axy_rhs->sgfn], cg->fgfs[Axz_rhs->sgfn],
                                cg->fgfs[Ayy_rhs->sgfn], cg->fgfs[Ayz_rhs->sgfn], cg->fgfs[Azz_rhs->sgfn],
                                cg->fgfs[Gmx_rhs->sgfn], cg->fgfs[Gmy_rhs->sgfn], cg->fgfs[Gmz_rhs->sgfn],
                                cg->fgfs[Lap_rhs->sgfn], 
                                cg->fgfs[Sfx_rhs->sgfn], cg->fgfs[Sfy_rhs->sgfn], cg->fgfs[Sfz_rhs->sgfn],
                                cg->fgfs[dtSfx_rhs->sgfn], cg->fgfs[dtSfy_rhs->sgfn], cg->fgfs[dtSfz_rhs->sgfn],
                                cg->fgfs[TZ_rhs->sgfn],
                                cg->fgfs[rho->sgfn], cg->fgfs[Sx->sgfn], cg->fgfs[Sy->sgfn], cg->fgfs[Sz->sgfn],
                                cg->fgfs[Sxx->sgfn], cg->fgfs[Sxy->sgfn], cg->fgfs[Sxz->sgfn], 
                                cg->fgfs[Syy->sgfn], cg->fgfs[Syz->sgfn], cg->fgfs[Szz->sgfn],
                                cg->fgfs[Gamxxx->sgfn], cg->fgfs[Gamxxy->sgfn], cg->fgfs[Gamxxz->sgfn],
                                cg->fgfs[Gamxyy->sgfn], cg->fgfs[Gamxyz->sgfn], cg->fgfs[Gamxzz->sgfn],
                                cg->fgfs[Gamyxx->sgfn], cg->fgfs[Gamyxy->sgfn], cg->fgfs[Gamyxz->sgfn],
                                cg->fgfs[Gamyyy->sgfn], cg->fgfs[Gamyyz->sgfn], cg->fgfs[Gamyzz->sgfn],
                                cg->fgfs[Gamzxx->sgfn], cg->fgfs[Gamzxy->sgfn], cg->fgfs[Gamzxz->sgfn],
                                cg->fgfs[Gamzyy->sgfn], cg->fgfs[Gamzyz->sgfn], cg->fgfs[Gamzzz->sgfn],
                                cg->fgfs[Rxx->sgfn], cg->fgfs[Rxy->sgfn], cg->fgfs[Rxz->sgfn], 
                                cg->fgfs[Ryy->sgfn], cg->fgfs[Ryz->sgfn], cg->fgfs[Rzz->sgfn],
                                cg->fgfs[Cons_Ham->sgfn],
                                cg->fgfs[Cons_Px->sgfn], cg->fgfs[Cons_Py->sgfn], cg->fgfs[Cons_Pz->sgfn],
                                cg->fgfs[Cons_Gx->sgfn], cg->fgfs[Cons_Gy->sgfn], cg->fgfs[Cons_Gz->sgfn],
                                Symmetry, lev, ndeps, pre);
            }
            if (BP == Pp->data->ble)
              break;
            BP = BP->next;
          }
          Pp = Pp->next;
        }
      }
      Parallel::Sync(GH->PatL[lev], ConstraintList, Symmetry);
    }
#ifdef WithShell
    SH->Synch(ConstraintList, Symmetry);
#endif

    double ConV[8];

#ifdef WithShell
    ConV[0] = SH->L2Norm(Cons_Ham);
    ConV[1] = SH->L2Norm(Cons_Px);
    ConV[2] = SH->L2Norm(Cons_Py);
    ConV[3] = SH->L2Norm(Cons_Pz);
    ConV[4] = SH->L2Norm(Cons_Gx);
    ConV[5] = SH->L2Norm(Cons_Gy);
    ConV[6] = SH->L2Norm(Cons_Gz);
    ConV[7] = SH->L2Norm(TZ0);
    ConVMonitor->writefile(PhysTime, 8, ConV);
#endif
    for (int levi = 0; levi < GH->levels; levi++)
    {
      ConV[0] = Parallel::L2Norm(GH->PatL[levi]->data, Cons_Ham);
      ConV[1] = Parallel::L2Norm(GH->PatL[levi]->data, Cons_Px);
      ConV[2] = Parallel::L2Norm(GH->PatL[levi]->data, Cons_Py);
      ConV[3] = Parallel::L2Norm(GH->PatL[levi]->data, Cons_Pz);
      ConV[4] = Parallel::L2Norm(GH->PatL[levi]->data, Cons_Gx);
      ConV[5] = Parallel::L2Norm(GH->PatL[levi]->data, Cons_Gy);
      ConV[6] = Parallel::L2Norm(GH->PatL[levi]->data, Cons_Gz);
      ConV[7] = Parallel::L2Norm(GH->PatL[levi]->data, TZ0);
      ConVMonitor->writefile(PhysTime, 8, ConV);
      /*
        if(fabs(ConV[0])<0.00001)
        {
          MyList<var> * DG_List=new MyList<var>(Cons_Ham);
                DG_List->insert(Cons_Px); DG_List->insert(Cons_Py); DG_List->insert(Cons_Px);
                DG_List->insert(Cons_Gx); DG_List->insert(Cons_Gy); DG_List->insert(Cons_Gx);
          Parallel::Dump_Data(GH->PatL[levi],DG_List,"jiu",0,1);
          DG_List->clearList();
          if(myrank==0) MPI_Abort(MPI_COMM_WORLD,1);
        }
      */
    }
    
    LastConsOut = 0;
  }
}

//================================================================================================



//================================================================================================

// 该成员函数用于对约束数据进行插值

//================================================================================================

void Z4c_class::Interp_Constraint()
{
  // we do not support a_lev != 0 yet.
  if (a_lev > 0)
    return;

  for (int lev = 0; lev < GH->levels; lev++)
  {
    // make sure the data consistent for higher levels
    if (lev > 0)
    {
      double TRK4 = PhysTime;
      double ndeps = numepsb;
      int pre = 0;
      MyList<Patch> *Pp = GH->PatL[lev];
      while (Pp)
      {
        MyList<Block> *BP = Pp->data->blb;
        while (BP)
        {
          Block *cg = BP->data;
          if (myrank == cg->rank)
          {
            f_compute_rhs_Z4c(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                              cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                              cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                              cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                              cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                              cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                              cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                              cg->fgfs[Lap0->sgfn], 
                              cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                              cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn],
                              cg->fgfs[TZ0->sgfn],
                              cg->fgfs[phi_rhs->sgfn], cg->fgfs[trK_rhs->sgfn],
                              cg->fgfs[gxx_rhs->sgfn], cg->fgfs[gxy_rhs->sgfn], cg->fgfs[gxz_rhs->sgfn],
                              cg->fgfs[gyy_rhs->sgfn], cg->fgfs[gyz_rhs->sgfn], cg->fgfs[gzz_rhs->sgfn],
                              cg->fgfs[Axx_rhs->sgfn], cg->fgfs[Axy_rhs->sgfn], cg->fgfs[Axz_rhs->sgfn],
                              cg->fgfs[Ayy_rhs->sgfn], cg->fgfs[Ayz_rhs->sgfn], cg->fgfs[Azz_rhs->sgfn],
                              cg->fgfs[Gmx_rhs->sgfn], cg->fgfs[Gmy_rhs->sgfn], cg->fgfs[Gmz_rhs->sgfn],
                              cg->fgfs[Lap_rhs->sgfn], 
                              cg->fgfs[Sfx_rhs->sgfn], cg->fgfs[Sfy_rhs->sgfn], cg->fgfs[Sfz_rhs->sgfn],
                              cg->fgfs[dtSfx_rhs->sgfn], cg->fgfs[dtSfy_rhs->sgfn], cg->fgfs[dtSfz_rhs->sgfn],
                              cg->fgfs[TZ_rhs->sgfn],
                              cg->fgfs[rho->sgfn], cg->fgfs[Sx->sgfn], cg->fgfs[Sy->sgfn], cg->fgfs[Sz->sgfn],
                              cg->fgfs[Sxx->sgfn], cg->fgfs[Sxy->sgfn], cg->fgfs[Sxz->sgfn], 
                              cg->fgfs[Syy->sgfn], cg->fgfs[Syz->sgfn], cg->fgfs[Szz->sgfn],
                              cg->fgfs[Gamxxx->sgfn], cg->fgfs[Gamxxy->sgfn], cg->fgfs[Gamxxz->sgfn],
                              cg->fgfs[Gamxyy->sgfn], cg->fgfs[Gamxyz->sgfn], cg->fgfs[Gamxzz->sgfn],
                              cg->fgfs[Gamyxx->sgfn], cg->fgfs[Gamyxy->sgfn], cg->fgfs[Gamyxz->sgfn],
                              cg->fgfs[Gamyyy->sgfn], cg->fgfs[Gamyyz->sgfn], cg->fgfs[Gamyzz->sgfn],
                              cg->fgfs[Gamzxx->sgfn], cg->fgfs[Gamzxy->sgfn], cg->fgfs[Gamzxz->sgfn],
                              cg->fgfs[Gamzyy->sgfn], cg->fgfs[Gamzyz->sgfn], cg->fgfs[Gamzzz->sgfn],
                              cg->fgfs[Rxx->sgfn], cg->fgfs[Rxy->sgfn], cg->fgfs[Rxz->sgfn], 
                              cg->fgfs[Ryy->sgfn], cg->fgfs[Ryz->sgfn], cg->fgfs[Rzz->sgfn],
                              cg->fgfs[Cons_Ham->sgfn],
                              cg->fgfs[Cons_Px->sgfn], cg->fgfs[Cons_Py->sgfn], cg->fgfs[Cons_Pz->sgfn],
                              cg->fgfs[Cons_Gx->sgfn], cg->fgfs[Cons_Gy->sgfn], cg->fgfs[Cons_Gz->sgfn],
                              Symmetry, lev, ndeps, pre);
          }
          if (BP == Pp->data->ble)
            break;
          BP = BP->next;
        }
        Pp = Pp->next;
      }
    }
    Parallel::Sync(GH->PatL[lev], ConstraintList, Symmetry);
  }
#ifdef WithShell
  SH->Synch(ConstraintList, Symmetry);
#endif
  //    interpolate
  double *x1, *y1, *z1;
  const int n = 1000;
  double lmax, lmin, dd;
  lmin = 0;
#ifdef WithShell
  lmax = SH->Rrange[1];
#else
  lmax = GH->bbox[0][0][4];
#endif
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
  dd = (lmax - lmin) / (n - 1);
#else
#ifdef Cell
  dd = (lmax - lmin) / n;
#else
#error Not define Vertex nor Cell
#endif
#endif
  x1 = new double[n];
  y1 = new double[n];
  z1 = new double[n];
  for (int i = 0; i < n; i++)
  {
    x1[i] = 0;
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
    y1[i] = lmin + i * dd;
#else
#ifdef Cell
    y1[i] = lmin + (i + 0.5) * dd;
#else
#error Not define Vertex nor Cell
#endif
#endif
    z1[i] = 0;
  }

  int InList = 0;

  MyList<var> *varl = ConstraintList;
  while (varl)
  {
    InList++;
    varl = varl->next;
  }
  double *shellf;
  shellf = new double[n * InList];
  for (int i = 0; i < n; i++)
  {
    double XX[3];
    XX[0] = x1[i];
    XX[1] = y1[i];
    XX[2] = z1[i];
    bool fg = GH->Interp_One_Point(ConstraintList, XX, shellf + i * InList, Symmetry);
#ifdef WithShell
    if (!fg)
      fg = SH->Interp_One_Point(ConstraintList, XX, shellf + i * InList, Symmetry);
#endif
    if (!fg && myrank == 0)
    {
      cout << "bssn_class::Interp_Constraint meets wrong" << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  ofstream outfile;
  char filename[50];
  sprintf(filename, "%s/interp_constraint_%05d.dat", ErrorMonitor->out_dir.c_str(), int(PhysTime / dT + 0.5)); 
  // 0.5 for round off
  
  outfile.open(filename);
  outfile << "#  corrdinate, H_Res, Px_Res, Py_Res, Pz_Res, Gx_Res, Gy_Res, Gz_Res, ...." << endl;
  for (int i = 0; i < n; i++)
  {
    outfile << setw(10) << setprecision(10) << y1[i];
    for (int j = 0; j < InList; j++)
      outfile << " " << setw(16) << setprecision(15) << shellf[InList * i + j];
    outfile << endl;
  }

  delete[] shellf;
}

//================================================================================================



//================================================================================================

// 该成员函数用于计算约束违反

//================================================================================================

void Z4c_class::Compute_Constraint()
{
  double TRK4 = PhysTime;
  double ndeps = numepsb;
  int pre = 0;
  int lev;

  for (lev = 0; lev < GH->levels; lev++)
  {
    {
      MyList<Patch> *Pp = GH->PatL[lev];
      while (Pp)
      {
        MyList<Block> *BP = Pp->data->blb;
        while (BP)
        {
          Block *cg = BP->data;
          if (myrank == cg->rank)
          {
            f_compute_rhs_bssn(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                               cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                               cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                               cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                               cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                               cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                               cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                               cg->fgfs[Lap0->sgfn], 
                               cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                               cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn],
                               cg->fgfs[phi_rhs->sgfn], cg->fgfs[trK_rhs->sgfn],
                               cg->fgfs[gxx_rhs->sgfn], cg->fgfs[gxy_rhs->sgfn], cg->fgfs[gxz_rhs->sgfn],
                               cg->fgfs[gyy_rhs->sgfn], cg->fgfs[gyz_rhs->sgfn], cg->fgfs[gzz_rhs->sgfn],
                               cg->fgfs[Axx_rhs->sgfn], cg->fgfs[Axy_rhs->sgfn], cg->fgfs[Axz_rhs->sgfn],
                               cg->fgfs[Ayy_rhs->sgfn], cg->fgfs[Ayz_rhs->sgfn], cg->fgfs[Azz_rhs->sgfn],
                               cg->fgfs[Gmx_rhs->sgfn], cg->fgfs[Gmy_rhs->sgfn], cg->fgfs[Gmz_rhs->sgfn],
                               cg->fgfs[Lap_rhs->sgfn], 
                               cg->fgfs[Sfx_rhs->sgfn], cg->fgfs[Sfy_rhs->sgfn], cg->fgfs[Sfz_rhs->sgfn],
                               cg->fgfs[dtSfx_rhs->sgfn], cg->fgfs[dtSfy_rhs->sgfn], cg->fgfs[dtSfz_rhs->sgfn],
                               cg->fgfs[rho->sgfn], cg->fgfs[Sx->sgfn], cg->fgfs[Sy->sgfn], cg->fgfs[Sz->sgfn],
                               cg->fgfs[Sxx->sgfn], cg->fgfs[Sxy->sgfn], cg->fgfs[Sxz->sgfn], 
                               cg->fgfs[Syy->sgfn], cg->fgfs[Syz->sgfn], cg->fgfs[Szz->sgfn],
                               cg->fgfs[Gamxxx->sgfn], cg->fgfs[Gamxxy->sgfn], cg->fgfs[Gamxxz->sgfn],
                               cg->fgfs[Gamxyy->sgfn], cg->fgfs[Gamxyz->sgfn], cg->fgfs[Gamxzz->sgfn],
                               cg->fgfs[Gamyxx->sgfn], cg->fgfs[Gamyxy->sgfn], cg->fgfs[Gamyxz->sgfn],
                               cg->fgfs[Gamyyy->sgfn], cg->fgfs[Gamyyz->sgfn], cg->fgfs[Gamyzz->sgfn],
                               cg->fgfs[Gamzxx->sgfn], cg->fgfs[Gamzxy->sgfn], cg->fgfs[Gamzxz->sgfn],
                               cg->fgfs[Gamzyy->sgfn], cg->fgfs[Gamzyz->sgfn], cg->fgfs[Gamzzz->sgfn],
                               cg->fgfs[Rxx->sgfn], cg->fgfs[Rxy->sgfn], cg->fgfs[Rxz->sgfn], 
                               cg->fgfs[Ryy->sgfn], cg->fgfs[Ryz->sgfn], cg->fgfs[Rzz->sgfn],
                               cg->fgfs[Cons_Ham->sgfn],
                               cg->fgfs[Cons_Px->sgfn], cg->fgfs[Cons_Py->sgfn], cg->fgfs[Cons_Pz->sgfn],
                               cg->fgfs[Cons_Gx->sgfn], cg->fgfs[Cons_Gy->sgfn], cg->fgfs[Cons_Gz->sgfn],
                               Symmetry, lev, ndeps, pre);
          }
          if (BP == Pp->data->ble)
            break;
          BP = BP->next;
        }
        Pp = Pp->next;
      }
    }
    Parallel::Sync(GH->PatL[lev], ConstraintList, Symmetry);
  }
#ifdef WithShell
  lev = 0;
  {
    MyList<ss_patch> *sPp;
    sPp = SH->PatL;
    while (sPp)
    {
      MyList<Block> *BP = sPp->data->blb;
      int fngfs = sPp->data->fngfs;
      while (BP)
      {
        Block *cg = BP->data;
        if (myrank == cg->rank)
        {
          f_compute_rhs_Z4c_ss(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                               cg->fgfs[fngfs + ShellPatch::gx], 
                               cg->fgfs[fngfs + ShellPatch::gy], 
                               cg->fgfs[fngfs + ShellPatch::gz],
                               cg->fgfs[fngfs + ShellPatch::drhodx], 
                               cg->fgfs[fngfs + ShellPatch::drhody], 
                               cg->fgfs[fngfs + ShellPatch::drhodz],
                               cg->fgfs[fngfs + ShellPatch::dsigmadx], 
                               cg->fgfs[fngfs + ShellPatch::dsigmady], 
                               cg->fgfs[fngfs + ShellPatch::dsigmadz],
                               cg->fgfs[fngfs + ShellPatch::dRdx], 
                               cg->fgfs[fngfs + ShellPatch::dRdy], 
                               cg->fgfs[fngfs + ShellPatch::dRdz],
                               cg->fgfs[fngfs + ShellPatch::drhodxx], 
                               cg->fgfs[fngfs + ShellPatch::drhodxy], 
                               cg->fgfs[fngfs + ShellPatch::drhodxz],
                               cg->fgfs[fngfs + ShellPatch::drhodyy], 
                               cg->fgfs[fngfs + ShellPatch::drhodyz], 
                               cg->fgfs[fngfs + ShellPatch::drhodzz],
                               cg->fgfs[fngfs + ShellPatch::dsigmadxx], 
                               cg->fgfs[fngfs + ShellPatch::dsigmadxy], 
                               cg->fgfs[fngfs + ShellPatch::dsigmadxz],
                               cg->fgfs[fngfs + ShellPatch::dsigmadyy], 
                               cg->fgfs[fngfs + ShellPatch::dsigmadyz], 
                               cg->fgfs[fngfs + ShellPatch::dsigmadzz],
                               cg->fgfs[fngfs + ShellPatch::dRdxx], 
                               cg->fgfs[fngfs + ShellPatch::dRdxy], 
                               cg->fgfs[fngfs + ShellPatch::dRdxz],
                               cg->fgfs[fngfs + ShellPatch::dRdyy], 
                               cg->fgfs[fngfs + ShellPatch::dRdyz], 
                               cg->fgfs[fngfs + ShellPatch::dRdzz],
                               cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                               cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                               cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                               cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                               cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                               cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                               cg->fgfs[Lap0->sgfn], 
                               cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                               cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn],
                               cg->fgfs[TZ0->sgfn],
                               cg->fgfs[phi_rhs->sgfn], cg->fgfs[trK_rhs->sgfn],
                               cg->fgfs[gxx_rhs->sgfn], cg->fgfs[gxy_rhs->sgfn], cg->fgfs[gxz_rhs->sgfn],
                               cg->fgfs[gyy_rhs->sgfn], cg->fgfs[gyz_rhs->sgfn], cg->fgfs[gzz_rhs->sgfn],
                               cg->fgfs[Axx_rhs->sgfn], cg->fgfs[Axy_rhs->sgfn], cg->fgfs[Axz_rhs->sgfn],
                               cg->fgfs[Ayy_rhs->sgfn], cg->fgfs[Ayz_rhs->sgfn], cg->fgfs[Azz_rhs->sgfn],
                               cg->fgfs[Gmx_rhs->sgfn], cg->fgfs[Gmy_rhs->sgfn], cg->fgfs[Gmz_rhs->sgfn],
                               cg->fgfs[Lap_rhs->sgfn], 
                               cg->fgfs[Sfx_rhs->sgfn], cg->fgfs[Sfy_rhs->sgfn], cg->fgfs[Sfz_rhs->sgfn],
                               cg->fgfs[dtSfx_rhs->sgfn], cg->fgfs[dtSfy_rhs->sgfn], cg->fgfs[dtSfz_rhs->sgfn],
                               cg->fgfs[TZ_rhs->sgfn],
                               cg->fgfs[rho->sgfn], cg->fgfs[Sx->sgfn], cg->fgfs[Sy->sgfn], cg->fgfs[Sz->sgfn],
                               cg->fgfs[Sxx->sgfn], cg->fgfs[Sxy->sgfn], cg->fgfs[Sxz->sgfn], 
                               cg->fgfs[Syy->sgfn], cg->fgfs[Syz->sgfn], cg->fgfs[Szz->sgfn],
                               cg->fgfs[Gamxxx->sgfn], cg->fgfs[Gamxxy->sgfn], cg->fgfs[Gamxxz->sgfn],
                               cg->fgfs[Gamxyy->sgfn], cg->fgfs[Gamxyz->sgfn], cg->fgfs[Gamxzz->sgfn],
                               cg->fgfs[Gamyxx->sgfn], cg->fgfs[Gamyxy->sgfn], cg->fgfs[Gamyxz->sgfn],
                               cg->fgfs[Gamyyy->sgfn], cg->fgfs[Gamyyz->sgfn], cg->fgfs[Gamyzz->sgfn],
                               cg->fgfs[Gamzxx->sgfn], cg->fgfs[Gamzxy->sgfn], cg->fgfs[Gamzxz->sgfn],
                               cg->fgfs[Gamzyy->sgfn], cg->fgfs[Gamzyz->sgfn], cg->fgfs[Gamzzz->sgfn],
                               cg->fgfs[Rxx->sgfn], cg->fgfs[Rxy->sgfn], cg->fgfs[Rxz->sgfn], 
                               cg->fgfs[Ryy->sgfn], cg->fgfs[Ryz->sgfn], cg->fgfs[Rzz->sgfn],
                               cg->fgfs[Cons_Ham->sgfn],
                               cg->fgfs[Cons_Px->sgfn], cg->fgfs[Cons_Py->sgfn], cg->fgfs[Cons_Pz->sgfn],
                               cg->fgfs[Cons_Gx->sgfn], cg->fgfs[Cons_Gy->sgfn], cg->fgfs[Cons_Gz->sgfn],
                               Symmetry, lev, numepsh, sPp->data->sst, pre);
        }
        if (BP == sPp->data->ble)
          break;
        BP = BP->next;
      }
      sPp = sPp->next;
    }
  }
  SH->Synch(ConstraintList, Symmetry);
#endif
}

//================================================================================================

