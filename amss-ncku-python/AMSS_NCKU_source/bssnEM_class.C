
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
#include "bssnEM_class.h"
#include "bssn_rhs.h"
#include "empart.h"
#include "initial_puncture.h"
#include "initial_maxwell.h"
#include "enforce_algebra.h"
#include "rungekutta4_rout.h"
#include "sommerfeld_rout.h"
#include "getnp4.h"
#include "getnpem2.h"
#include "shellfunctions.h"
#include "parameters.h"

#ifdef With_AHF
#include "derivatives.h"
#include "myglobal.h"
#endif

//================================================================================================

// 定义 bssnEM_class

// 它是继承了父类 bssn_class 的某些成员和方法，并对另一些成员和方法进行修改
// 修改的成员和方法在下面（以及头文件 bssnEM_class.h ）中定义
// 其余的继承父类 bssn_class（头文件 bssn_class.h 中声明）

//================================================================================================

bssnEM_class::bssnEM_class(double Couranti, double StartTimei, double TotalTimei, 
                           double DumpTimei, double d2DumpTimei, double CheckTimei, double AnasTimei,
                           int Symmetryi, int checkruni, char *checkfilenamei, 
                           double numepssi, double numepsbi, double numepshi,
                           int a_levi, int maxli, int decni, double maxrexi, double drexi) 
                           : bssn_class(Couranti, StartTimei, TotalTimei, 
                                        DumpTimei, d2DumpTimei, CheckTimei, AnasTimei,
                                        Symmetryi, checkruni, checkfilenamei, numepssi, numepsbi, numepshi,
                                        a_levi, maxli, decni, maxrexi, drexi)
{
  // setup Monitors
  {
    char str[50];
    stringstream a_stream;
    a_stream.setf(ios::left);
    a_stream.str("");
    a_stream << setw(15) << "# time";
    for (int pl = 1; pl < maxl + 1; pl++)
      for (int pm = -pl; pm < pl + 1; pm++)
      {
        sprintf(str, "R%02dm%03d", pl, pm);
        a_stream << setw(16) << str;
        sprintf(str, "I%02dm%03d", pl, pm);
        a_stream << setw(16) << str;
      }
    Phi2Monitor = new monitor("bssn_phi2.dat", myrank, a_stream.str()); // myrank has been setup in bssn_class.C

    a_stream.clear();
    a_stream.str("");
    a_stream << setw(15) << "# time";
    for (int pl = 0; pl < maxl + 1; pl++)
      for (int pm = -pl; pm < pl + 1; pm++)
      {
        sprintf(str, "R%02dm%03d", pl, pm);
        a_stream << setw(16) << str;
        sprintf(str, "I%02dm%03d", pl, pm);
        a_stream << setw(16) << str;
      }
    Phi1Monitor = new monitor("bssn_phi1.dat", myrank, a_stream.str()); // myrank has been setup in bssn_class.C
  }
}

//================================================================================================



//================================================================================================

// 该成员函数用于对 Class 初始化

//================================================================================================

void bssnEM_class::Initialize()
{
  Exo = new var("Exo", ngfs++, -1, 1, 1);
  Eyo = new var("Eyo", ngfs++, 1, -1, 1);
  Ezo = new var("Ezo", ngfs++, 1, 1, -1);
  // note B is an axi vector
  Bxo = new var("Bxo", ngfs++, 1, -1, -1);
  Byo = new var("Byo", ngfs++, -1, 1, -1);
  Bzo = new var("Bzo", ngfs++, -1, -1, 1);
  Kpsio = new var("Kpsio", ngfs++, 1, 1, 1);
  Kphio = new var("Kphio", ngfs++, 1, 1, 1);

  Ex0 = new var("Ex0", ngfs++, -1, 1, 1);
  Ey0 = new var("Ey0", ngfs++, 1, -1, 1);
  Ez0 = new var("Ez0", ngfs++, 1, 1, -1);
  Bx0 = new var("Bx0", ngfs++, 1, -1, -1);
  By0 = new var("By0", ngfs++, -1, 1, -1);
  Bz0 = new var("Bz0", ngfs++, -1, -1, 1);
  Kpsi0 = new var("Kpsi0", ngfs++, 1, 1, 1);
  Kphi0 = new var("Kphi0", ngfs++, 1, 1, 1);

  Ex = new var("Ex", ngfs++, -1, 1, 1);
  Ey = new var("Ey", ngfs++, 1, -1, 1);
  Ez = new var("Ez", ngfs++, 1, 1, -1);
  Bx = new var("Bx", ngfs++, 1, -1, -1);
  By = new var("By", ngfs++, -1, 1, -1);
  Bz = new var("Bz", ngfs++, -1, -1, 1);
  Kpsi = new var("Kpsi", ngfs++, 1, 1, 1);
  Kphi = new var("Kphi", ngfs++, 1, 1, 1);

  Ex1 = new var("Ex1", ngfs++, -1, 1, 1);
  Ey1 = new var("Ey1", ngfs++, 1, -1, 1);
  Ez1 = new var("Ez1", ngfs++, 1, 1, -1);
  Bx1 = new var("Bx1", ngfs++, 1, -1, -1);
  By1 = new var("By1", ngfs++, -1, 1, -1);
  Bz1 = new var("Bz1", ngfs++, -1, -1, 1);
  Kpsi1 = new var("Kpsi1", ngfs++, 1, 1, 1);
  Kphi1 = new var("Kphi1", ngfs++, 1, 1, 1);

  Ex_rhs = new var("Ex_rhs", ngfs++, -1, 1, 1);
  Ey_rhs = new var("Ey_rhs", ngfs++, 1, -1, 1);
  Ez_rhs = new var("Ez_rhs", ngfs++, 1, 1, -1);
  Bx_rhs = new var("Bx_rhs", ngfs++, 1, -1, -1);
  By_rhs = new var("By_rhs", ngfs++, -1, 1, -1);
  Bz_rhs = new var("Bz_rhs", ngfs++, -1, -1, 1);
  Kpsi_rhs = new var("Kpsi_rhs", ngfs++, 1, 1, 1);
  Kphi_rhs = new var("Kphi_rhs", ngfs++, 1, 1, 1);

  qchar = new var("qchar", ngfs++, 1, 1, 1);
  Jx = new var("Jx", ngfs++, -1, 1, 1);
  Jy = new var("Jy", ngfs++, 1, -1, 1);
  Jz = new var("Jz", ngfs++, 1, 1, -1);

  Rphi2 = new var("Rphi2", ngfs++, 1, 1, 1);    // Etheta - Bphi in fact, so no symmetry at all
  Iphi2 = new var("Iphi2", ngfs++, -1, -1, -1); // Ephi - Btheta in fact, so no symmetry at all

  Rphi1 = new var("Rphi1", ngfs++, 1, 1, 1); // Er in fact
  Iphi1 = new var("Iphi1", ngfs++, 1, 1, 1); // Br in fact

  if (myrank == 0)
    cout << "you have setted " << ngfs << " grid functions." << endl;

  OldStateList->insert(Kpsio);
  OldStateList->insert(Kphio);
  OldStateList->insert(Exo);
  OldStateList->insert(Eyo);
  OldStateList->insert(Ezo);
  OldStateList->insert(Bxo);
  OldStateList->insert(Byo);
  OldStateList->insert(Bzo);

  StateList->insert(Kpsi0);
  StateList->insert(Kphi0);
  StateList->insert(Ex0);
  StateList->insert(Ey0);
  StateList->insert(Ez0);
  StateList->insert(Bx0);
  StateList->insert(By0);
  StateList->insert(Bz0);

  RHSList->insert(Kpsi_rhs);
  RHSList->insert(Kphi_rhs);
  RHSList->insert(Ex_rhs);
  RHSList->insert(Ey_rhs);
  RHSList->insert(Ez_rhs);
  RHSList->insert(Bx_rhs);
  RHSList->insert(By_rhs);
  RHSList->insert(Bz_rhs);

  SynchList_pre->insert(Kpsi);
  SynchList_pre->insert(Kphi);
  SynchList_pre->insert(Ex);
  SynchList_pre->insert(Ey);
  SynchList_pre->insert(Ez);
  SynchList_pre->insert(Bx);
  SynchList_pre->insert(By);
  SynchList_pre->insert(Bz);

  SynchList_cor->insert(Kpsi1);
  SynchList_cor->insert(Kphi1);
  SynchList_cor->insert(Ex1);
  SynchList_cor->insert(Ey1);
  SynchList_cor->insert(Ez1);
  SynchList_cor->insert(Bx1);
  SynchList_cor->insert(By1);
  SynchList_cor->insert(Bz1);

  DumpList->insert(Rphi2);
  DumpList->insert(Iphi2);
  DumpList->insert(Rphi1);
  DumpList->insert(Iphi1);
  DumpList->insert(Ex0);
  DumpList->insert(Bx0);

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

bssnEM_class::~bssnEM_class()
{
  delete Kpsio;
  delete Kphio;
  delete Exo;
  delete Eyo;
  delete Ezo;
  delete Bxo;
  delete Byo;
  delete Bzo;

  delete Kpsi0;
  delete Kphi0;
  delete Ex0;
  delete Ey0;
  delete Ez0;
  delete Bx0;
  delete By0;
  delete Bz0;

  delete Kpsi;
  delete Kphi;
  delete Ex;
  delete Ey;
  delete Ez;
  delete Bx;
  delete By;
  delete Bz;

  delete Kpsi1;
  delete Kphi1;
  delete Ex1;
  delete Ey1;
  delete Ez1;
  delete Bx1;
  delete By1;
  delete Bz1;

  delete Kpsi_rhs;
  delete Kphi_rhs;
  delete Ex_rhs;
  delete Ey_rhs;
  delete Ez_rhs;
  delete Bx_rhs;
  delete By_rhs;
  delete Bz_rhs;

  delete qchar;
  delete Jx;
  delete Jy;
  delete Jz;

  delete Rphi2;
  delete Iphi2;

  delete Rphi1;
  delete Iphi1;

  delete Phi2Monitor;

  delete Phi1Monitor;
}

//================================================================================================



//================================================================================================

// 该成员函数读入 Ansorg 方法求解的 TwoPuncture 初值

//================================================================================================

// Read initial data solved by Ansorg, PRD 70, 064011 (2004)

void bssnEM_class::Read_Ansorg()
{
  if (checkrun)
  {
    CheckPoint->readcheck_cgh(PhysTime, GH, myrank, nprocs, Symmetry);
#ifdef WithShell
    CheckPoint->readcheck_sh(SH, myrank);
#endif
  }
  else
  {
    if (myrank == 0)
      cout << "Read initial data from Ansorg's solver,"
           << " please be sure the input parameters for black holes are puncture parameters!!" << endl;
    char filename[50];
    {
      map<string, string>::iterator iter = parameters::str_par.find("inputpar");
      if (iter != parameters::str_par.end())
      {
        strcpy(filename, (iter->second).c_str());
      }
      else
      {
        cout << "Error inputpar" << endl;
        exit(0);
      }
    }
    int BH_NM;
    double *Porg_here, *Qchar;
    // read parameter from file
    {
      const int LEN = 256;
      char pline[LEN];
      string str, sgrp, skey, sval;
      int sind;
      ifstream inf(filename, ifstream::in);
      if (!inf.good() && myrank == 0)
      {
        if (ErrorMonitor->outfile)
          ErrorMonitor->outfile << "Can not open parameter file " << filename 
                                << " for inputing information of black holes" << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      for (int i = 1; inf.good(); i++)
      {
        inf.getline(pline, LEN);
        str = pline;

        int status = misc::parse_parts(str, sgrp, skey, sval, sind);
        if (status == -1)
        {
          if (ErrorMonitor->outfile)
            ErrorMonitor->outfile << "error reading parameter file " << filename 
                                  << " in line " << i << endl;
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else if (status == 0)
          continue;

        if (sgrp == "BSSN" && skey == "BH_num")
        {
          BH_NM = atoi(sval.c_str());
          break;
        }
      }
      inf.close();
    }

    Porg_here = new double[3 * BH_NM];
    Pmom = new double[3 * BH_NM];
    Spin = new double[3 * BH_NM];
    Mass = new double[BH_NM];
    Qchar = new double[BH_NM];
    // read parameter from file
    {
      const int LEN = 256;
      char pline[LEN];
      string str, sgrp, skey, sval;
      int sind;
      ifstream inf(filename, ifstream::in);
      if (!inf.good() && myrank == 0)
      {
        if (ErrorMonitor->outfile)
          ErrorMonitor->outfile << "Can not open parameter file " << filename
                                << " for inputing information of black holes" << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      for (int i = 1; inf.good(); i++)
      {
        inf.getline(pline, LEN);
        str = pline;

        int status = misc::parse_parts(str, sgrp, skey, sval, sind);
        if (status == -1)
        {
          if (ErrorMonitor->outfile)
            ErrorMonitor->outfile << "error reading parameter file " << filename 
                                  << " in line " << i << endl;
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else if (status == 0)
          continue;

        if (sgrp == "BSSN" && sind < BH_NM)
        {
          if (skey == "Mass")
            Mass[sind] = atof(sval.c_str());
          else if (skey == "Qchar")
          {
            Qchar[sind] = atof(sval.c_str());
            if (myrank == 0)
              cout << "black hole #" << sind << " has elctric charge " << Qchar[sind] << endl;
          }
          else if (skey == "Porgx")
            Porg_here[sind * 3] = atof(sval.c_str());
          else if (skey == "Porgy")
            Porg_here[sind * 3 + 1] = atof(sval.c_str());
          else if (skey == "Porgz")
            Porg_here[sind * 3 + 2] = atof(sval.c_str());
          else if (skey == "Spinx")
            Spin[sind * 3] = atof(sval.c_str());
          else if (skey == "Spiny")
            Spin[sind * 3 + 1] = atof(sval.c_str());
          else if (skey == "Spinz")
            Spin[sind * 3 + 2] = atof(sval.c_str());
          else if (skey == "Pmomx")
            Pmom[sind * 3] = atof(sval.c_str());
          else if (skey == "Pmomy")
            Pmom[sind * 3 + 1] = atof(sval.c_str());
          else if (skey == "Pmomz")
            Pmom[sind * 3 + 2] = atof(sval.c_str());
        }
      }
      inf.close();
    }
    int order = 6;
    Ansorg read_ansorg("Ansorg.psid", order);
    // set initial data
    for (int lev = 0; lev < GH->levels; lev++)
    {
      MyList<Patch> *Pp = GH->PatL[lev];
      while (Pp)
      {
        MyList<Block> *BL = Pp->data->blb;
        while (BL)
        {
          Block *cg = BL->data;
          if (myrank == cg->rank)
          {
            for (int k = 0; k < cg->shape[2]; k++)
              for (int j = 0; j < cg->shape[1]; j++)
                for (int i = 0; i < cg->shape[0]; i++)
                  cg->fgfs[phi0->sgfn][i + j * cg->shape[0] + k * cg->shape[0] * cg->shape[1]] =
                      read_ansorg.ps_u_at_xyz(cg->X[0][i], cg->X[1][j], cg->X[2][k]);

            f_get_ansorg_nbhs_em(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                 cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                                 cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                                 cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                                 cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                                 cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                                 cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                                 cg->fgfs[Lap0->sgfn], 
                                 cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                                 cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn],
                                 cg->fgfs[Ex0->sgfn], cg->fgfs[Ey0->sgfn], cg->fgfs[Ez0->sgfn],
                                 cg->fgfs[Bx0->sgfn], cg->fgfs[By0->sgfn], cg->fgfs[Bz0->sgfn],
                                 cg->fgfs[Kpsi0->sgfn], cg->fgfs[Kphi0->sgfn],
                                 Mass, Qchar, Porg_here, Pmom, Spin, BH_NM);
          }
          if (BL == Pp->data->ble)
            break;
          BL = BL->next;
        }
        Pp = Pp->next;
      }
    }
#ifdef WithShell
    // ShellPatch part
    MyList<ss_patch> *Pp = SH->PatL;
    while (Pp)
    {
      MyList<Block> *BL = Pp->data->blb;
      while (BL)
      {
        Block *cg = BL->data;
        if (myrank == cg->rank)
        {
          for (int k = 0; k < cg->shape[2]; k++)
            for (int j = 0; j < cg->shape[1]; j++)
              for (int i = 0; i < cg->shape[0]; i++)
                cg->fgfs[phi0->sgfn][i + j * cg->shape[0] + k * cg->shape[0] * cg->shape[1]] =
                    read_ansorg.ps_u_at_xyz(cg->fgfs[Pp->data->fngfs + ShellPatch::gx][i + j * cg->shape[0] + k * cg->shape[0] * cg->shape[1]],
                                            cg->fgfs[Pp->data->fngfs + ShellPatch::gy][i + j * cg->shape[0] + k * cg->shape[0] * cg->shape[1]],
                                            cg->fgfs[Pp->data->fngfs + ShellPatch::gz][i + j * cg->shape[0] + k * cg->shape[0] * cg->shape[1]]);

          f_get_ansorg_nbhs_ss_em(cg->shape, 
                                  cg->fgfs[Pp->data->fngfs + ShellPatch::gx], 
                                  cg->fgfs[Pp->data->fngfs + ShellPatch::gy],
                                  cg->fgfs[Pp->data->fngfs + ShellPatch::gz],
                                  cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                                  cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                                  cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                                  cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                                  cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                                  cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                                  cg->fgfs[Lap0->sgfn], 
                                  cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                                  cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn],
                                  cg->fgfs[Ex0->sgfn], cg->fgfs[Ey0->sgfn], cg->fgfs[Ez0->sgfn],
                                  cg->fgfs[Bx0->sgfn], cg->fgfs[By0->sgfn], cg->fgfs[Bz0->sgfn],
                                  cg->fgfs[Kpsi0->sgfn], cg->fgfs[Kphi0->sgfn],
                                  Mass, Qchar, Porg_here, Pmom, Spin, BH_NM);
        }
        if (BL == Pp->data->ble)
          break;
        BL = BL->next;
      }
      Pp = Pp->next;
    }
#endif

    delete[] Porg_here;
// dump read_in initial data
//   for(int lev=0;lev<GH->levels;lev++) Parallel::Dump_Data(GH->PatL[lev],StateList,0,PhysTime,dT);
// check initial constraint
#if 0 
    for(int lev=0;lev<GH->levels;lev++) Step(lev,0);
    if(myrank == 0) MPI_Abort(MPI_COMM_WORLD,1);
#endif
  }
}

//================================================================================================



//================================================================================================

// 该成员函数用于用解析函数设定带电磁场的数值相对论初值
// 看下面的描述，仅对对头碰撞 head on 的情况

//================================================================================================

// Set up initial data given by PRD 80, 104022 (2009)
void bssnEM_class::Setup_Initial_Data()
{
  if (checkrun)
  {
    CheckPoint->readcheck_cgh(PhysTime, GH, myrank, nprocs, Symmetry);
#ifdef WithShell
    CheckPoint->readcheck_sh(SH, myrank);
#endif
  }
  else
  {
    if (myrank == 0)
      cout << "Setup initial data for head on identical charge-mass ratio black holes." << endl;
    char filename[50];
    {
      map<string, string>::iterator iter = parameters::str_par.find("inputpar");
      if (iter != parameters::str_par.end())
      {
        strcpy(filename, (iter->second).c_str());
      }
      else
      {
        cout << "Error inputpar" << endl;
        exit(0);
      }
    }
    int BH_NM;
    double *Porg_here, *Qchar_here, *Pmom_here, *Spin_here, *Mass_here;
    // read parameter from file
    {
      const int LEN = 256;
      char pline[LEN];
      string str, sgrp, skey, sval;
      int sind;
      ifstream inf(filename, ifstream::in);
      if (!inf.good() && myrank == 0)
      {
        if (ErrorMonitor->outfile)
          ErrorMonitor->outfile << "Can not open parameter file " << filename 
                                << " for inputing information of black holes" << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      for (int i = 1; inf.good(); i++)
      {
        inf.getline(pline, LEN);
        str = pline;

        int status = misc::parse_parts(str, sgrp, skey, sval, sind);
        if (status == -1)
        {
          if (ErrorMonitor->outfile)
            ErrorMonitor->outfile << "error reading parameter file " << filename 
                                  << " in line " << i << endl;
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else if (status == 0)
          continue;

        if (sgrp == "BSSN" && skey == "BH_num")
        {
          BH_NM = atoi(sval.c_str());
          break;
        }
      }
      inf.close();
    }

    Porg_here = new double[3 * BH_NM];
    Pmom_here = new double[3 * BH_NM];
    Spin_here = new double[3 * BH_NM];
    Mass_here = new double[BH_NM];
    Qchar_here = new double[BH_NM];
    // read parameter from file
    {
      const int LEN = 256;
      char pline[LEN];
      string str, sgrp, skey, sval;
      int sind;
      ifstream inf(filename, ifstream::in);
      if (!inf.good() && myrank == 0)
      {
        if (ErrorMonitor->outfile)
          ErrorMonitor->outfile << "Can not open parameter file " << filename
                                << " for inputing information of black holes" << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      for (int i = 1; inf.good(); i++)
      {
        inf.getline(pline, LEN);
        str = pline;

        int status = misc::parse_parts(str, sgrp, skey, sval, sind);
        if (status == -1)
        {
          if (ErrorMonitor->outfile)
            ErrorMonitor->outfile << "error reading parameter file " << filename 
                                  << " in line " << i << endl;
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else if (status == 0)
          continue;

        if (sgrp == "BSSN" && sind < BH_NM)
        {
          if (skey == "Mass")
            Mass_here[sind] = atof(sval.c_str());
          else if (skey == "Qchar")
          {
            Qchar_here[sind] = atof(sval.c_str());
            if (myrank == 0)
              cout << "black hole #" << sind << " has elctric charge " << Qchar_here[sind] << endl;
          }
          else if (skey == "Porgx")
            Porg_here[sind * 3] = atof(sval.c_str());
          else if (skey == "Porgy")
            Porg_here[sind * 3 + 1] = atof(sval.c_str());
          else if (skey == "Porgz")
            Porg_here[sind * 3 + 2] = atof(sval.c_str());
          else if (skey == "Spinx")
            Spin_here[sind * 3] = atof(sval.c_str());
          else if (skey == "Spiny")
            Spin_here[sind * 3 + 1] = atof(sval.c_str());
          else if (skey == "Spinz")
            Spin_here[sind * 3 + 2] = atof(sval.c_str());
          else if (skey == "Pmomx")
            Pmom_here[sind * 3] = atof(sval.c_str());
          else if (skey == "Pmomy")
            Pmom_here[sind * 3 + 1] = atof(sval.c_str());
          else if (skey == "Pmomz")
            Pmom_here[sind * 3 + 2] = atof(sval.c_str());
        }
      }
      inf.close();
    }
    // set initial data
    for (int lev = 0; lev < GH->levels; lev++)
    {
      MyList<Patch> *Pp = GH->PatL[lev];
      while (Pp)
      {
        MyList<Block> *BL = Pp->data->blb;
        while (BL)
        {
          Block *cg = BL->data;
          if (myrank == cg->rank)
          {
            f_get_initial_nbhsem(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                 cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                                 cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                                 cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                                 cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                                 cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                                 cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                                 cg->fgfs[Lap0->sgfn], 
                                 cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                                 cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn],
                                 cg->fgfs[Ex0->sgfn], cg->fgfs[Ey0->sgfn], cg->fgfs[Ez0->sgfn],
                                 cg->fgfs[Bx0->sgfn], cg->fgfs[By0->sgfn], cg->fgfs[Bz0->sgfn],
                                 cg->fgfs[Kpsi0->sgfn], cg->fgfs[Kphi0->sgfn],
                                 Mass_here, Qchar_here, Porg_here, Pmom_here, Spin_here, BH_NM);
          }
          if (BL == Pp->data->ble)
            break;
          BL = BL->next;
        }
        Pp = Pp->next;
      }
    }
#ifdef WithShell
    // ShellPatch part
    MyList<ss_patch> *Pp = SH->PatL;
    while (Pp)
    {
      MyList<Block> *BL = Pp->data->blb;
      while (BL)
      {
        Block *cg = BL->data;
        if (myrank == cg->rank)
        {
          f_get_initial_nbhsem_ss(cg->shape, 
                                  cg->fgfs[Pp->data->fngfs + ShellPatch::gx], 
                                  cg->fgfs[Pp->data->fngfs + ShellPatch::gy],
                                  cg->fgfs[Pp->data->fngfs + ShellPatch::gz],
                                  cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                                  cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                                  cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                                  cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                                  cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                                  cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                                  cg->fgfs[Lap0->sgfn], 
                                  cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                                  cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn],
                                  cg->fgfs[Ex0->sgfn], cg->fgfs[Ey0->sgfn], cg->fgfs[Ez0->sgfn],
                                  cg->fgfs[Bx0->sgfn], cg->fgfs[By0->sgfn], cg->fgfs[Bz0->sgfn],
                                  cg->fgfs[Kpsi0->sgfn], cg->fgfs[Kphi0->sgfn],
                                  Mass_here, Qchar_here, Porg_here, Pmom_here, Spin_here, BH_NM);
        }
        if (BL == Pp->data->ble)
          break;
        BL = BL->next;
      }
      Pp = Pp->next;
    }
#endif

    delete[] Porg_here;
    delete[] Mass_here;
    delete[] Qchar_here;
    delete[] Pmom_here;
    delete[] Spin_here;
    // dump read_in initial data
    //   for(int lev=0;lev<GH->levels;lev++) Parallel::Dump_Data(GH->PatL[lev],StateList,0,PhysTime,dT);
  }
}

//================================================================================================



//================================================================================================

// 该成员函数设定了时间演化过程中的单步时间演化

//================================================================================================

void bssnEM_class::Step(int lev, int YN)
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

        if (
            f_compute_rhs_empart(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                 cg->fgfs[phi0->sgfn],
                                 cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                                 cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                                 cg->fgfs[Lap0->sgfn], 
                                 cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn], 
                                 cg->fgfs[trK0->sgfn],
                                 cg->fgfs[Ex0->sgfn], cg->fgfs[Ey0->sgfn], cg->fgfs[Ez0->sgfn], 
                                 cg->fgfs[Bx0->sgfn], cg->fgfs[By0->sgfn], cg->fgfs[Bz0->sgfn],
                                 cg->fgfs[Kpsi0->sgfn], cg->fgfs[Kphi0->sgfn], 
                                 cg->fgfs[Jx->sgfn], cg->fgfs[Jy->sgfn], cg->fgfs[Jz->sgfn], 
                                 cg->fgfs[qchar->sgfn],
                                 cg->fgfs[Ex_rhs->sgfn], cg->fgfs[Ey_rhs->sgfn], cg->fgfs[Ez_rhs->sgfn],
                                 cg->fgfs[Bx_rhs->sgfn], cg->fgfs[By_rhs->sgfn], cg->fgfs[Bz_rhs->sgfn],
                                 cg->fgfs[Kpsi_rhs->sgfn], cg->fgfs[Kphi_rhs->sgfn],
                                 cg->fgfs[rho->sgfn], 
                                 cg->fgfs[Sx->sgfn], cg->fgfs[Sy->sgfn], cg->fgfs[Sz->sgfn],
                                 cg->fgfs[Sxx->sgfn], cg->fgfs[Sxy->sgfn], cg->fgfs[Sxz->sgfn], 
                                 cg->fgfs[Syy->sgfn], cg->fgfs[Syz->sgfn], cg->fgfs[Szz->sgfn],
                                 Symmetry, lev, ndeps) ||
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
          MyList<var> *varl0 = StateList, *varl = SynchList_pre, *varlrhs = RHSList; 
          // we do not check the correspondence here
          
          while (varl0)
          {
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
            f_rungekutta4_rout(cg->shape, dT_lev, 
                               cg->fgfs[varl0->data->sgfn], 
                               cg->fgfs[varl->data->sgfn], 
                               cg->fgfs[varlrhs->data->sgfn],
                               iter_count);
#ifndef WithShell
            if (lev > 0) // fix BD point
#endif
              f_sommerfeld_rout(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], 
                                Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                                dT_lev, cg->fgfs[phi0->sgfn],
                                cg->fgfs[Lap0->sgfn], 
                                cg->fgfs[varl0->data->sgfn], cg->fgfs[varl->data->sgfn], 
                                varl0->data->SoA,
                                Symmetry, cor);

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

// check initial constraint
#if 0
    Parallel::Dump_Data(GH->PatL[lev],DumpList,0,PhysTime,dT_lev);
#endif

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

          if (
              f_compute_rhs_empart_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
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
                                      cg->fgfs[phi0->sgfn],
                                      cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                                      cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                                      cg->fgfs[Lap0->sgfn], 
                                      cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn], 
                                      cg->fgfs[trK0->sgfn],
                                      cg->fgfs[Ex0->sgfn], cg->fgfs[Ey0->sgfn], cg->fgfs[Ez0->sgfn], 
                                      cg->fgfs[Bx0->sgfn], cg->fgfs[By0->sgfn], cg->fgfs[Bz0->sgfn],
                                      cg->fgfs[Kpsi0->sgfn], cg->fgfs[Kphi0->sgfn], 
                                      cg->fgfs[Jx->sgfn], cg->fgfs[Jy->sgfn], cg->fgfs[Jz->sgfn], 
                                      cg->fgfs[qchar->sgfn],
                                      cg->fgfs[Ex_rhs->sgfn], cg->fgfs[Ey_rhs->sgfn], cg->fgfs[Ez_rhs->sgfn],
                                      cg->fgfs[Bx_rhs->sgfn], cg->fgfs[By_rhs->sgfn], cg->fgfs[Bz_rhs->sgfn],
                                      cg->fgfs[Kpsi_rhs->sgfn], cg->fgfs[Kphi_rhs->sgfn],
                                      cg->fgfs[rho->sgfn], 
                                      cg->fgfs[Sx->sgfn], cg->fgfs[Sy->sgfn], cg->fgfs[Sz->sgfn],
                                      cg->fgfs[Sxx->sgfn], cg->fgfs[Sxy->sgfn], cg->fgfs[Sxz->sgfn], 
                                      cg->fgfs[Syy->sgfn], cg->fgfs[Syz->sgfn], cg->fgfs[Szz->sgfn],
                                      Symmetry, lev, numepsh, sPp->data->sst) ||
              f_compute_rhs_bssn_ss(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
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
                                      cg->fgfs[varl0->data->sgfn], varl0->data->propspeed, varl0->data->SoA,
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
#if 1
            // falloff boundary condition
            {
              int n = 2;
              f_falloff_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                           sPp->data->bbox[0], sPp->data->bbox[1], sPp->data->bbox[2], 
                           sPp->data->bbox[3], sPp->data->bbox[4], sPp->data->bbox[5],
                           cg->fgfs[Ex->sgfn], n, Ex->SoA, Symmetry);
              f_falloff_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                           sPp->data->bbox[0], sPp->data->bbox[1], sPp->data->bbox[2], 
                           sPp->data->bbox[3], sPp->data->bbox[4], sPp->data->bbox[5],
                           cg->fgfs[Ey->sgfn], n, Ey->SoA, Symmetry);
              f_falloff_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                           sPp->data->bbox[0], sPp->data->bbox[1], sPp->data->bbox[2], 
                           sPp->data->bbox[3], sPp->data->bbox[4], sPp->data->bbox[5],
                           cg->fgfs[Ez->sgfn], n, Ez->SoA, Symmetry);
              f_falloff_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                           sPp->data->bbox[0], sPp->data->bbox[1], sPp->data->bbox[2], 
                           sPp->data->bbox[3], sPp->data->bbox[4], sPp->data->bbox[5],
                           cg->fgfs[Bx->sgfn], n, Bx->SoA, Symmetry);
              f_falloff_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                           sPp->data->bbox[0], sPp->data->bbox[1], sPp->data->bbox[2], 
                           sPp->data->bbox[3], sPp->data->bbox[4], sPp->data->bbox[5],
                           cg->fgfs[By->sgfn], n, By->SoA, Symmetry);
              f_falloff_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                           sPp->data->bbox[0], sPp->data->bbox[1], sPp->data->bbox[2], 
                           sPp->data->bbox[3], sPp->data->bbox[4], sPp->data->bbox[5],
                           cg->fgfs[Bz->sgfn], n, Bz->SoA, Symmetry);
              n = 3;
              f_falloff_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                           sPp->data->bbox[0], sPp->data->bbox[1], sPp->data->bbox[2], 
                           sPp->data->bbox[3], sPp->data->bbox[4], sPp->data->bbox[5],
                           cg->fgfs[Kpsi->sgfn], n, Kpsi->SoA, Symmetry);
              f_falloff_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                           sPp->data->bbox[0], sPp->data->bbox[1], sPp->data->bbox[2], 
                           sPp->data->bbox[3], sPp->data->bbox[4], sPp->data->bbox[5],
                           cg->fgfs[Kphi->sgfn], n, Kphi->SoA, Symmetry);
            }
#endif
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
    AnalysisStuff_EM(lev, dT_lev);
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

          if (
              f_compute_rhs_empart(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                   cg->fgfs[phi->sgfn],
                                   cg->fgfs[gxx->sgfn], cg->fgfs[gxy->sgfn], cg->fgfs[gxz->sgfn], 
                                   cg->fgfs[gyy->sgfn], cg->fgfs[gyz->sgfn], cg->fgfs[gzz->sgfn],
                                   cg->fgfs[Lap->sgfn], 
                                   cg->fgfs[Sfx->sgfn], cg->fgfs[Sfy->sgfn], cg->fgfs[Sfz->sgfn], 
                                   cg->fgfs[trK->sgfn],
                                   cg->fgfs[Ex->sgfn], cg->fgfs[Ey->sgfn], cg->fgfs[Ez->sgfn], 
                                   cg->fgfs[Bx->sgfn], cg->fgfs[By->sgfn], cg->fgfs[Bz->sgfn],
                                   cg->fgfs[Kpsi->sgfn], cg->fgfs[Kphi->sgfn], 
                                   cg->fgfs[Jx->sgfn], cg->fgfs[Jy->sgfn], cg->fgfs[Jz->sgfn], 
                                   cg->fgfs[qchar->sgfn],
                                   cg->fgfs[Ex1->sgfn], cg->fgfs[Ey1->sgfn], cg->fgfs[Ez1->sgfn], 
                                   cg->fgfs[Bx1->sgfn], cg->fgfs[By1->sgfn], cg->fgfs[Bz1->sgfn],
                                   cg->fgfs[Kpsi1->sgfn], cg->fgfs[Kphi1->sgfn],
                                   cg->fgfs[rho->sgfn], 
                                   cg->fgfs[Sx->sgfn], cg->fgfs[Sy->sgfn], cg->fgfs[Sz->sgfn],
                                   cg->fgfs[Sxx->sgfn], cg->fgfs[Sxy->sgfn], cg->fgfs[Sxz->sgfn], 
                                   cg->fgfs[Syy->sgfn], cg->fgfs[Syz->sgfn], cg->fgfs[Szz->sgfn],
                                   Symmetry, lev, ndeps) ||
              f_compute_rhs_bssn(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                                 cg->fgfs[phi->sgfn], cg->fgfs[trK->sgfn],
                                 cg->fgfs[gxx->sgfn], cg->fgfs[gxy->sgfn], cg->fgfs[gxz->sgfn], 
                                 cg->fgfs[gyy->sgfn], cg->fgfs[gyz->sgfn], cg->fgfs[gzz->sgfn],
                                 cg->fgfs[Axx->sgfn], cg->fgfs[Axy->sgfn], cg->fgfs[Axz->sgfn], 
                                 cg->fgfs[Ayy->sgfn], cg->fgfs[Ayz->sgfn], cg->fgfs[Azz->sgfn],
                                 cg->fgfs[Gmx->sgfn], cg->fgfs[Gmy->sgfn], cg->fgfs[Gmz->sgfn],
                                 cg->fgfs[Lap->sgfn], 
                                 cg->fgfs[Sfx->sgfn], cg->fgfs[Sfy->sgfn], cg->fgfs[Sfz->sgfn],
                                 cg->fgfs[dtSfx->sgfn], cg->fgfs[dtSfy->sgfn], cg->fgfs[dtSfz->sgfn],
                                 cg->fgfs[phi1->sgfn], cg->fgfs[trK1->sgfn],
                                 cg->fgfs[gxx1->sgfn], cg->fgfs[gxy1->sgfn], cg->fgfs[gxz1->sgfn],
                                 cg->fgfs[gyy1->sgfn], cg->fgfs[gyz1->sgfn], cg->fgfs[gzz1->sgfn],
                                 cg->fgfs[Axx1->sgfn], cg->fgfs[Axy1->sgfn], cg->fgfs[Axz1->sgfn],
                                 cg->fgfs[Ayy1->sgfn], cg->fgfs[Ayz1->sgfn], cg->fgfs[Azz1->sgfn],
                                 cg->fgfs[Gmx1->sgfn], cg->fgfs[Gmy1->sgfn], cg->fgfs[Gmz1->sgfn],
                                 cg->fgfs[Lap1->sgfn], cg->fgfs[Sfx1->sgfn], cg->fgfs[Sfy1->sgfn], cg->fgfs[Sfz1->sgfn],
                                 cg->fgfs[dtSfx1->sgfn], cg->fgfs[dtSfy1->sgfn], cg->fgfs[dtSfz1->sgfn],
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
              f_rungekutta4_rout(cg->shape, dT_lev, 
                                 cg->fgfs[varl0->data->sgfn], 
                                 cg->fgfs[varl1->data->sgfn], 
                                 cg->fgfs[varlrhs->data->sgfn],
                                 iter_count);

#ifndef WithShell
              if (lev > 0) // fix BD point
#endif
                f_sommerfeld_rout(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                  Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], 
                                  Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                                  dT_lev, cg->fgfs[phi0->sgfn],
                                  cg->fgfs[Lap0->sgfn], 
                                  cg->fgfs[varl0->data->sgfn], cg->fgfs[varl1->data->sgfn], 
                                  varl0->data->SoA,
                                  Symmetry, cor);

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

            if (
                f_compute_rhs_empart_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
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
                                        cg->fgfs[phi->sgfn],
                                        cg->fgfs[gxx->sgfn], cg->fgfs[gxy->sgfn], cg->fgfs[gxz->sgfn], 
                                        cg->fgfs[gyy->sgfn], cg->fgfs[gyz->sgfn], cg->fgfs[gzz->sgfn],
                                        cg->fgfs[Lap->sgfn], cg->fgfs[Sfx->sgfn], cg->fgfs[Sfy->sgfn], 
                                        cg->fgfs[Sfz->sgfn], cg->fgfs[trK->sgfn],
                                        cg->fgfs[Ex->sgfn], cg->fgfs[Ey->sgfn], cg->fgfs[Ez->sgfn], 
                                        cg->fgfs[Bx->sgfn], cg->fgfs[By->sgfn], cg->fgfs[Bz->sgfn],
                                        cg->fgfs[Kpsi->sgfn], cg->fgfs[Kphi->sgfn], 
                                        cg->fgfs[Jx->sgfn], cg->fgfs[Jy->sgfn], cg->fgfs[Jz->sgfn], 
                                        cg->fgfs[qchar->sgfn],
                                        cg->fgfs[Ex1->sgfn], cg->fgfs[Ey1->sgfn], cg->fgfs[Ez1->sgfn], 
                                        cg->fgfs[Bx1->sgfn], cg->fgfs[By1->sgfn], cg->fgfs[Bz1->sgfn],
                                        cg->fgfs[Kpsi1->sgfn], cg->fgfs[Kphi1->sgfn],
                                        cg->fgfs[rho->sgfn], 
                                        cg->fgfs[Sx->sgfn], cg->fgfs[Sy->sgfn], cg->fgfs[Sz->sgfn],
                                        cg->fgfs[Sxx->sgfn], cg->fgfs[Sxy->sgfn], cg->fgfs[Sxz->sgfn], 
                                        cg->fgfs[Syy->sgfn], cg->fgfs[Syz->sgfn], cg->fgfs[Szz->sgfn],
                                        Symmetry, lev, numepsh, sPp->data->sst) ||
                f_compute_rhs_bssn_ss(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
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
                                      cg->fgfs[phi1->sgfn], cg->fgfs[trK1->sgfn],
                                      cg->fgfs[gxx1->sgfn], cg->fgfs[gxy1->sgfn], cg->fgfs[gxz1->sgfn],
                                      cg->fgfs[gyy1->sgfn], cg->fgfs[gyz1->sgfn], cg->fgfs[gzz1->sgfn],
                                      cg->fgfs[Axx1->sgfn], cg->fgfs[Axy1->sgfn], cg->fgfs[Axz1->sgfn],
                                      cg->fgfs[Ayy1->sgfn], cg->fgfs[Ayz1->sgfn], cg->fgfs[Azz1->sgfn],
                                      cg->fgfs[Gmx1->sgfn], cg->fgfs[Gmy1->sgfn], cg->fgfs[Gmz1->sgfn],
                                      cg->fgfs[Lap1->sgfn], 
                                      cg->fgfs[Sfx1->sgfn], cg->fgfs[Sfy1->sgfn], cg->fgfs[Sfz1->sgfn],
                                      cg->fgfs[dtSfx1->sgfn], cg->fgfs[dtSfy1->sgfn], cg->fgfs[dtSfz1->sgfn],
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
#if 1
            // falloff boundary condition
            {
              int n = 2;
              f_falloff_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                           sPp->data->bbox[0], sPp->data->bbox[1], sPp->data->bbox[2], 
                           sPp->data->bbox[3], sPp->data->bbox[4], sPp->data->bbox[5],
                           cg->fgfs[Ex1->sgfn], n, Ex1->SoA, Symmetry);
              f_falloff_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                           sPp->data->bbox[0], sPp->data->bbox[1], sPp->data->bbox[2], 
                           sPp->data->bbox[3], sPp->data->bbox[4], sPp->data->bbox[5],
                           cg->fgfs[Ey1->sgfn], n, Ey1->SoA, Symmetry);
              f_falloff_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                           sPp->data->bbox[0], sPp->data->bbox[1], sPp->data->bbox[2], 
                           sPp->data->bbox[3], sPp->data->bbox[4], sPp->data->bbox[5],
                           cg->fgfs[Ez1->sgfn], n, Ez1->SoA, Symmetry);
              f_falloff_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                           sPp->data->bbox[0], sPp->data->bbox[1], sPp->data->bbox[2], 
                           sPp->data->bbox[3], sPp->data->bbox[4], sPp->data->bbox[5],
                           cg->fgfs[Bx1->sgfn], n, Bx1->SoA, Symmetry);
              f_falloff_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                           sPp->data->bbox[0], sPp->data->bbox[1], sPp->data->bbox[2], 
                           sPp->data->bbox[3], sPp->data->bbox[4], sPp->data->bbox[5],
                           cg->fgfs[By1->sgfn], n, By1->SoA, Symmetry);
              f_falloff_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                           sPp->data->bbox[0], sPp->data->bbox[1], sPp->data->bbox[2], 
                           sPp->data->bbox[3], sPp->data->bbox[4], sPp->data->bbox[5],
                           cg->fgfs[Bz1->sgfn], n, Bz1->SoA, Symmetry);
              n = 3;
              f_falloff_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                           sPp->data->bbox[0], sPp->data->bbox[1], sPp->data->bbox[2], 
                           sPp->data->bbox[3], sPp->data->bbox[4], sPp->data->bbox[5],
                           cg->fgfs[Kpsi1->sgfn], n, Kpsi1->SoA, Symmetry);
              f_falloff_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                           sPp->data->bbox[0], sPp->data->bbox[1], sPp->data->bbox[2], 
                           sPp->data->bbox[3], sPp->data->bbox[4], sPp->data->bbox[5],
                           cg->fgfs[Kphi1->sgfn], n, Kphi1->SoA, Symmetry);
            }
#endif
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
      cout << " CS_Inter used " 
           << (double)(curr_clock - prev_clock) / ((double)CLOCKS_PER_SEC) 
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

//================================================================================================



//================================================================================================

// 该成员函数用于计算电磁波 Phi2

//================================================================================================

void bssnEM_class::Compute_Phi2(int lev)
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
        f_getnpem2(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                   cg->fgfs[phi0->sgfn],
                   cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                   cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                   cg->fgfs[Ex0->sgfn], cg->fgfs[Ey0->sgfn], cg->fgfs[Ez0->sgfn], 
                   cg->fgfs[Bx0->sgfn], cg->fgfs[By0->sgfn], cg->fgfs[Bz0->sgfn],
                   cg->fgfs[Rphi2->sgfn], cg->fgfs[Iphi2->sgfn],
                   Symmetry);
        f_getnpem1(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                   cg->fgfs[phi0->sgfn],
                   cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                   cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                   cg->fgfs[Ex0->sgfn], cg->fgfs[Ey0->sgfn], cg->fgfs[Ez0->sgfn], 
                   cg->fgfs[Bx0->sgfn], cg->fgfs[By0->sgfn], cg->fgfs[Bz0->sgfn],
                   cg->fgfs[Rphi1->sgfn], cg->fgfs[Iphi1->sgfn],
                   Symmetry);
      }
      if (BP == Pp->data->ble)
        break;
      BP = BP->next;
    }
    Pp = Pp->next;
  }

#ifdef WithShell
  // ShellPatch part
  if (lev == 0)
  {
    MyList<ss_patch> *Pp = SH->PatL;
    while (Pp)
    {
      MyList<Block> *BL = Pp->data->blb;
      int fngfs = Pp->data->fngfs;
      while (BL)
      {
        Block *cg = BL->data;
        if (myrank == cg->rank)
        {
          f_getnpem2_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
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
                        cg->fgfs[phi0->sgfn],
                        cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                        cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                        cg->fgfs[Ex0->sgfn], cg->fgfs[Ey0->sgfn], cg->fgfs[Ez0->sgfn], 
                        cg->fgfs[Bx0->sgfn], cg->fgfs[By0->sgfn], cg->fgfs[Bz0->sgfn],
                        cg->fgfs[Rphi2->sgfn], cg->fgfs[Iphi2->sgfn],
                        Symmetry, Pp->data->sst);
          f_getnpem1_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
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
                        cg->fgfs[phi0->sgfn],
                        cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                        cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                        cg->fgfs[Ex0->sgfn], cg->fgfs[Ey0->sgfn], cg->fgfs[Ez0->sgfn], 
                        cg->fgfs[Bx0->sgfn], cg->fgfs[By0->sgfn], cg->fgfs[Bz0->sgfn],
                        cg->fgfs[Rphi1->sgfn], cg->fgfs[Iphi1->sgfn],
                        Symmetry, Pp->data->sst);
        }
        if (BL == Pp->data->ble)
          break;
        BL = BL->next;
      }
      Pp = Pp->next;
    }
  }
#endif

  MyList<var> *DG_List = new MyList<var>(Rphi2);
  DG_List->insert(Iphi2);
  DG_List->insert(Rphi1);
  DG_List->insert(Iphi1);
  Parallel::Sync(GH->PatL[lev], DG_List, Symmetry);

#ifdef WithShell
  if (lev == 0)
  {
    SH->Synch(DG_List, Symmetry);
  }
#endif
  DG_List->clearList();
}

//================================================================================================



//================================================================================================

// 该成员函数用于分析电磁场的数据

//================================================================================================

void bssnEM_class::AnalysisStuff_EM(int lev, double dT_lev)
{
  LastAnas += dT_lev;
  
  if (LastAnas >= AnasTime)
  {
    Compute_Phi2(lev);
    double *RP, *IP;
    int NN = 0;
    // for phi2
    for (int pl = 1; pl < maxl + 1; pl++)
      for (int pm = -pl; pm < pl + 1; pm++)
        NN++;
    RP = new double[NN];
    IP = new double[NN];
    double Rex = maxrex;
    for (int i = 0; i < decn; i++)
    {
#ifdef WithShell
      if (lev > 0 || Rex < GH->bbox[0][0][3])
      {
        //            Waveshell->surf_Wave(Rex,lev,GH, Rphi2, Iphi2,1,maxl,NN,RP,IP,ErrorMonitor);
        Waveshell->surf_Wave(Rex, lev, GH, 
                             Ex0, Ey0, Ez0, Bx0, By0, Bz0, phi0, 
                             gxx0, gxy0, gxz0, gyy0, gyz0, gzz0, 
                             1, maxl, NN, RP, IP, ErrorMonitor,
                             f_getnpem2_point);
      }
      else
      {
        //	    Waveshell->surf_Wave(Rex,lev,SH, Rphi2, Iphi2,1,maxl,NN,RP,IP,ErrorMonitor);
        //	    Waveshell->surf_Wave(Rex,lev,SH, Ex0,Ey0,Ez0,Bx0,By0,Bz0,phi0,gxx0,gxy0,gxz0,gyy0,gyz0,gzz0,1,maxl,NN,RP,IP,ErrorMonitor);
        Waveshell->surf_Wave(Rex, lev, SH, 
                             Ex0, Ey0, Ez0, Bx0, By0, Bz0, phi0, 
                             gxx0, gxy0, gxz0, gyy0, gyz0, gzz0, 
                             1, maxl, NN, RP, IP, ErrorMonitor,
                             f_getnpem2_point);
      }
#else
      //            Waveshell->surf_Wave(Rex,lev,GH, Rphi2, Iphi2,1,maxl,NN,RP,IP,ErrorMonitor);
      Waveshell->surf_Wave(Rex, lev, GH, 
                           Ex0, Ey0, Ez0, Bx0, By0, Bz0, phi0, 
                           gxx0, gxy0, gxz0, gyy0, gyz0, gzz0, 
                           1, maxl, NN, RP, IP, ErrorMonitor,
                           f_getnpem2_point);
#endif
      Phi2Monitor->writefile(PhysTime, NN, RP, IP);
      Rex = Rex - drex;
    }
    delete[] RP;
    delete[] IP;

    // for phi1
    NN = 0;
    for (int pl = 0; pl < maxl + 1; pl++)
      for (int pm = -pl; pm < pl + 1; pm++)
        NN++;
    RP = new double[NN];
    IP = new double[NN];
    Rex = maxrex;
    for (int i = 0; i < decn; i++)
    {
#ifdef WithShell
      if (lev > 0 || Rex < GH->bbox[0][0][3])
      {
        //            Waveshell->surf_Wave(Rex,lev,GH, Rphi1, Iphi1,0,maxl,NN,RP,IP,ErrorMonitor);
        Waveshell->surf_Wave(Rex, lev, GH, 
                             Ex0, Ey0, Ez0, Bx0, By0, Bz0, phi0, 
                             gxx0, gxy0, gxz0, gyy0, gyz0, gzz0, 
                             0, maxl, NN, RP, IP, ErrorMonitor,
                             f_getnpem1_point);
      }
      else
      {
        //	    Waveshell->surf_Wave(Rex,lev,SH, Rphi1, Iphi1,0,maxl,NN,RP,IP,ErrorMonitor);
        Waveshell->surf_Wave(Rex, lev, SH, 
                             Ex0, Ey0, Ez0, Bx0, By0, Bz0, phi0, 
                             gxx0, gxy0, gxz0, gyy0, gyz0, gzz0, 
                             0, maxl, NN, RP, IP, ErrorMonitor,
                             f_getnpem1_point);
      }
#else
      //            Waveshell->surf_Wave(Rex,lev,GH, Rphi1, Iphi1,0,maxl,NN,RP,IP,ErrorMonitor);
      Waveshell->surf_Wave(Rex, lev, GH, 
                           Ex0, Ey0, Ez0, Bx0, By0, Bz0, phi0, 
                           gxx0, gxy0, gxz0, gyy0, gyz0, gzz0, 
                           0, maxl, NN, RP, IP, ErrorMonitor,
                           f_getnpem1_point);
#endif
      Phi1Monitor->writefile(PhysTime, NN, RP, IP);
      Rex = Rex - drex;
    }
    delete[] RP;
    delete[] IP;
  }

  AnalysisStuff(lev, dT_lev); // LastAnas need and only need control here
  
  // 这是共享变量？每次分析完要归零？
  LastAnas = 0;
}

//================================================================================================



//================================================================================================

// 该成员函数用于对约束数据进行插值

//================================================================================================

void bssnEM_class::Interp_Constraint()
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
            f_compute_rhs_empart(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                 cg->fgfs[phi0->sgfn],
                                 cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                                 cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                                 cg->fgfs[Lap0->sgfn], 
                                 cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn], 
                                 cg->fgfs[trK0->sgfn],
                                 cg->fgfs[Ex0->sgfn], cg->fgfs[Ey0->sgfn], cg->fgfs[Ez0->sgfn], 
                                 cg->fgfs[Bx0->sgfn], cg->fgfs[By0->sgfn], cg->fgfs[Bz0->sgfn],
                                 cg->fgfs[Kpsi0->sgfn], cg->fgfs[Kphi0->sgfn], 
                                 cg->fgfs[Jx->sgfn], cg->fgfs[Jy->sgfn], cg->fgfs[Jz->sgfn], 
                                 cg->fgfs[qchar->sgfn],
                                 cg->fgfs[Ex_rhs->sgfn], cg->fgfs[Ey_rhs->sgfn], cg->fgfs[Ez_rhs->sgfn],
                                 cg->fgfs[Bx_rhs->sgfn], cg->fgfs[By_rhs->sgfn], cg->fgfs[Bz_rhs->sgfn],
                                 cg->fgfs[Kpsi_rhs->sgfn], cg->fgfs[Kphi_rhs->sgfn],
                                 cg->fgfs[rho->sgfn], 
                                 cg->fgfs[Sx->sgfn], cg->fgfs[Sy->sgfn], cg->fgfs[Sz->sgfn],
                                 cg->fgfs[Sxx->sgfn], cg->fgfs[Sxy->sgfn], cg->fgfs[Sxz->sgfn], 
                                 cg->fgfs[Syy->sgfn], cg->fgfs[Syz->sgfn], cg->fgfs[Szz->sgfn],
                                 Symmetry, lev, ndeps) ||
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

