
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
#include "bssnEScalar_class.h"
#include "bssn_rhs.h"
#include "initial_puncture.h"
#include "enforce_algebra.h"
#include "rungekutta4_rout.h"
#include "sommerfeld_rout.h"
#include "getnp4.h"
#include "shellfunctions.h"
#include "parameters.h"

#ifdef With_AHF
#include "derivatives.h"
#include "myglobal.h"
#endif

//================================================================================================

// 定义 bssnEScalar_class

// 它是继承了父类 bssn_class 的某些成员和方法，并对另一些成员和方法进行修改
// 修改的成员和方法在下面（以及头文件 bssnEScalar_class.h ）中定义
// 其余的继承父类 bssn_class（头文件 bssn_class.h 中声明）

//================================================================================================

bssnEScalar_class::bssnEScalar_class(double Couranti, double StartTimei, double TotalTimei, 
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
  // setup Monitors
  {
    char str[50];
    stringstream a_stream;
    a_stream.setf(ios::left);
    a_stream.str("");
    a_stream << setw(15) << "# time x y z maxs";
    MaxScalar_Monitor = new monitor("bssn_maxs.dat", myrank, a_stream.str()); 
    // myrank has been setup in bssn_class.C
  }
}

//================================================================================================



//================================================================================================

// 该成员函数用于对 Class 初始化

//================================================================================================

void bssnEScalar_class::Initialize()
{
  Sphio = new var("Sphio", ngfs++, 1, 1, 1);
  Spio = new var("Spio", ngfs++, 1, 1, 1);
  Sphi0 = new var("Sphi0", ngfs++, 1, 1, 1);
  Spi0 = new var("Spi0", ngfs++, 1, 1, 1);
  Sphi = new var("Sphi", ngfs++, 1, 1, 1);
  Spi = new var("Spi", ngfs++, 1, 1, 1);
  Sphi1 = new var("Sphi1", ngfs++, 1, 1, 1);
  Spi1 = new var("Spi1", ngfs++, 1, 1, 1);
  Sphi_rhs = new var("Sphi_rhs", ngfs++, 1, 1, 1);
  Spi_rhs = new var("Spi_rhs", ngfs++, 1, 1, 1);

  // constraint violation monitor variables
  Cons_fR = new var("Cons_fR", ngfs++, 1, 1, 1);

  if (myrank == 0)
    cout << "you have setted " << ngfs << " grid functions." << endl;

  OldStateList->insert(Sphio);
  OldStateList->insert(Spio);
  StateList->insert(Sphi0);
  StateList->insert(Spi0);
  RHSList->insert(Sphi_rhs);
  RHSList->insert(Spi_rhs);
  SynchList_pre->insert(Sphi);
  SynchList_pre->insert(Spi);
  SynchList_cor->insert(Sphi1);
  SynchList_cor->insert(Spi1);

  ConstraintList->insert(Cons_Gz);

  DumpList->insert(Sphi0);
  DumpList->insert(Spi0);
  DumpList->insert(Cons_fR);

  CheckPoint->addvariablelist(StateList);
  CheckPoint->addvariablelist(OldStateList);
  
    
  int myrank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    
  // read parameter from file
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

bssnEScalar_class::~bssnEScalar_class()
{
  delete Sphio;
  delete Spio;
  delete Sphi0;
  delete Spi0;
  delete Sphi;
  delete Spi;
  delete Sphi1;
  delete Spi1;
  delete Sphi_rhs;
  delete Spi_rhs;

  delete Cons_fR;

  delete MaxScalar_Monitor;
}

//================================================================================================



//================================================================================================

// 该成员函数读入 Ansorg 方法求解的 TwoPuncture 初值

//================================================================================================

// Read initial data solved by Ansorg, PRD 70, 064011 (2004)

void bssnEScalar_class::Read_Ansorg()
{
  if (!checkrun)
  {
    if (myrank == 0)
      cout << "Read initial data from Ansorg's solver,"
           << " please be sure the input parameters for black holes are puncture parameters!!" 
           << endl;
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
    double *Porg_here;
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
            ErrorMonitor->outfile << "error reading parameter file " 
                                  << filename << " in line " << i << endl;
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

            f_get_ansorg_nbhs_escalar(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                      cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                                      cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                                      cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                                      cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                                      cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                                      cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                                      cg->fgfs[Lap0->sgfn], 
                                      cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                                      cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn],
                                      cg->fgfs[Sphi0->sgfn], cg->fgfs[Spi0->sgfn],
                                      Mass, Porg_here, Pmom, Spin, BH_NM);
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

          f_get_ansorg_nbhs_ss_escalar(cg->shape, 
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
                                       cg->fgfs[Sphi0->sgfn], cg->fgfs[Spi0->sgfn],
                                       Mass, Porg_here, Pmom, Spin, BH_NM);
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
  }
}

//================================================================================================



//================================================================================================

// 该成员函数读入 Pablo Galaviz 的 Olliptic 程序求出的初值
// 但 Olliptic 程序已无人维护，不知如何运行

//================================================================================================

// Read initial data solved by Pablo's Olliptic Phys.Rev.D 82 024005 (2010)

void bssnEScalar_class::Read_Pablo()
{
  if (!checkrun)
  {
    if (myrank == 0)
      cout << "Read initial data from Pablo's solver,"
           << " please be sure the input parameters for black holes are puncture parameters!!" 
           << endl;
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
    double *Porg_here;
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
    bool flag = false;
    int DIM = dim;
    // set initial data
    for (int lev = 0; lev < GH->levels; lev++)
    {
      MyList<Patch> *Pp = GH->PatL[lev];
      int grd = 0;
      while (Pp)
      {
        double *databuffer = (double *)malloc(sizeof(double) 
                             * Pp->data->shape[0] * Pp->data->shape[1] * Pp->data->shape[2]);
        if (!databuffer)
        {
          cout << "bssnEScalar_class::Read_Pablo: on node# " << myrank 
               << ", out of memory when reading Pablo's data in" << endl;
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
        char filename[100];
        sprintf(filename, "Lev%02d-%02d.mgid_m", lev, grd);
        if (read_Pablo_file((int *)Pp->data->shape, databuffer, filename))
        {
          MyList<Block> *BL = Pp->data->blb;
          while (BL)
          {
            Block *cg = BL->data;
            if (myrank == cg->rank)
            {
              f_copy(DIM, cg->bbox, cg->bbox + DIM, cg->shape, cg->fgfs[phi0->sgfn],
                     Pp->data->bbox, Pp->data->bbox + DIM, Pp->data->shape, databuffer,
                     cg->bbox, cg->bbox + DIM);

              f_get_ansorg_nbhs_escalar(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                        cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                                        cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                                        cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                                        cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                                        cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                                        cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                                        cg->fgfs[Lap0->sgfn], 
                                        cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                                        cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn],
                                        cg->fgfs[Sphi0->sgfn], cg->fgfs[Spi0->sgfn],
                                        Mass, Porg_here, Pmom, Spin, BH_NM);
            }
            if (BL == Pp->data->ble)
              break;
            BL = BL->next;
          }
        }
        else
        {
          sprintf(filename, "Lev%02d-%02d.mgid", lev, grd);
          if (myrank == 0)
            write_Pablo_file((int *)Pp->data->shape, 
                             Pp->data->bbox[0], Pp->data->bbox[3], 
                             Pp->data->bbox[1], Pp->data->bbox[4],
                             Pp->data->bbox[2], Pp->data->bbox[5], 
                             filename);
          flag = true;
        }
        free(databuffer);
        Pp = Pp->next;
        grd++;
      }
    }

#ifdef WithShell
    // ShellPatch part
    MyList<ss_patch> *Pp = SH->PatL;
    while (Pp)
    {
      double *databuffer = (double *)malloc(sizeof(double) * Pp->data->shape[0] * Pp->data->shape[1] * Pp->data->shape[2]);
      if (!databuffer)
      {
        cout << "bssnEScalar_class::Read_Pablo: on node# " << myrank << ", out of memory when reading Pablo's data in" << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
      char filename[100], shn[10];
      SH->shellname(shn, Pp->data->sst);
      sprintf(filename, "LevSH-%s.mgid_m", shn);
      if (read_Pablo_file((int *)Pp->data->shape, databuffer, filename))
      {
        MyList<Block> *BL = Pp->data->blb;
        while (BL)
        {
          Block *cg = BL->data;
          if (myrank == cg->rank)
          {
            f_copy(DIM, cg->bbox, cg->bbox + DIM, cg->shape, cg->fgfs[phi0->sgfn],
                   Pp->data->bbox, Pp->data->bbox + DIM, Pp->data->shape, databuffer,
                   cg->bbox, cg->bbox + DIM);

            f_get_ansorg_nbhs_ss_escalar(cg->shape, 
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
                                         cg->fgfs[Sphi0->sgfn], cg->fgfs[Spi0->sgfn],
                                         Mass, Porg_here, Pmom, Spin, BH_NM);
          }
          if (BL == Pp->data->ble)
            break;
          BL = BL->next;
        }
      }
      else
      {
        sprintf(filename, "LevSH-%s.mgid", shn);
        if (myrank == 0)
          SH->write_Pablo_file_ss((int *)Pp->data->shape, 
                                  Pp->data->bbox[0], Pp->data->bbox[3], 
                                  Pp->data->bbox[1], Pp->data->bbox[4],
                                  Pp->data->bbox[2], Pp->data->bbox[5], 
                                  filename, Pp->data->sst);
        flag = true;
      }
      free(databuffer);
      Pp = Pp->next;
    }
#endif

    delete[] Porg_here;
    if (flag && myrank == 0)
      MPI_Abort(MPI_COMM_WORLD, 1);
    // dump read_in initial data
    for (int lev = 0; lev < GH->levels; lev++)
      Parallel::Dump_Data(GH->PatL[lev], StateList, 0, PhysTime, dT);
    SH->Dump_Data(StateList, 0, PhysTime, dT);
  }
}

//================================================================================================



//================================================================================================

// 该成员函数设定了时间演化过程中的单步时间演化

//================================================================================================

void bssnEScalar_class::Step(int lev, int YN)
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

        if (f_compute_rhs_bssn_escalar(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                                       cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                                       cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                                       cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                                       cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                                       cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                                       cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                                       cg->fgfs[Lap0->sgfn], 
                                       cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                                       cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn],
                                       cg->fgfs[Sphi0->sgfn], cg->fgfs[Spi0->sgfn],
                                       cg->fgfs[phi_rhs->sgfn], cg->fgfs[trK_rhs->sgfn],
                                       cg->fgfs[gxx_rhs->sgfn], cg->fgfs[gxy_rhs->sgfn], cg->fgfs[gxz_rhs->sgfn],
                                       cg->fgfs[gyy_rhs->sgfn], cg->fgfs[gyz_rhs->sgfn], cg->fgfs[gzz_rhs->sgfn],
                                       cg->fgfs[Axx_rhs->sgfn], cg->fgfs[Axy_rhs->sgfn], cg->fgfs[Axz_rhs->sgfn],
                                       cg->fgfs[Ayy_rhs->sgfn], cg->fgfs[Ayz_rhs->sgfn], cg->fgfs[Azz_rhs->sgfn],
                                       cg->fgfs[Gmx_rhs->sgfn], cg->fgfs[Gmy_rhs->sgfn], cg->fgfs[Gmz_rhs->sgfn],
                                       cg->fgfs[Lap_rhs->sgfn], 
                                       cg->fgfs[Sfx_rhs->sgfn], cg->fgfs[Sfy_rhs->sgfn], cg->fgfs[Sfz_rhs->sgfn],
                                       cg->fgfs[dtSfx_rhs->sgfn], cg->fgfs[dtSfy_rhs->sgfn], cg->fgfs[dtSfz_rhs->sgfn],
                                       cg->fgfs[Sphi_rhs->sgfn], cg->fgfs[Spi_rhs->sgfn],
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

          if (f_compute_rhs_bssn_escalar_ss(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
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
                                            cg->fgfs[Sphi0->sgfn], cg->fgfs[Spi0->sgfn],
                                            cg->fgfs[phi_rhs->sgfn], cg->fgfs[trK_rhs->sgfn],
                                            cg->fgfs[gxx_rhs->sgfn], cg->fgfs[gxy_rhs->sgfn], cg->fgfs[gxz_rhs->sgfn],
                                            cg->fgfs[gyy_rhs->sgfn], cg->fgfs[gyz_rhs->sgfn], cg->fgfs[gzz_rhs->sgfn],
                                            cg->fgfs[Axx_rhs->sgfn], cg->fgfs[Axy_rhs->sgfn], cg->fgfs[Axz_rhs->sgfn],
                                            cg->fgfs[Ayy_rhs->sgfn], cg->fgfs[Ayz_rhs->sgfn], cg->fgfs[Azz_rhs->sgfn],
                                            cg->fgfs[Gmx_rhs->sgfn], cg->fgfs[Gmy_rhs->sgfn], cg->fgfs[Gmz_rhs->sgfn],
                                            cg->fgfs[Lap_rhs->sgfn], 
                                            cg->fgfs[Sfx_rhs->sgfn], cg->fgfs[Sfy_rhs->sgfn], cg->fgfs[Sfz_rhs->sgfn],
                                            cg->fgfs[dtSfx_rhs->sgfn], cg->fgfs[dtSfy_rhs->sgfn], cg->fgfs[dtSfz_rhs->sgfn],
                                            cg->fgfs[Sphi_rhs->sgfn], cg->fgfs[Spi_rhs->sgfn],
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
    AnalysisStuff_EScalar(lev, dT_lev);
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

          if (f_compute_rhs_bssn_escalar(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                                         cg->fgfs[phi->sgfn], cg->fgfs[trK->sgfn],
                                         cg->fgfs[gxx->sgfn], cg->fgfs[gxy->sgfn], cg->fgfs[gxz->sgfn], 
                                         cg->fgfs[gyy->sgfn], cg->fgfs[gyz->sgfn], cg->fgfs[gzz->sgfn],
                                         cg->fgfs[Axx->sgfn], cg->fgfs[Axy->sgfn], cg->fgfs[Axz->sgfn], 
                                         cg->fgfs[Ayy->sgfn], cg->fgfs[Ayz->sgfn], cg->fgfs[Azz->sgfn],
                                         cg->fgfs[Gmx->sgfn], cg->fgfs[Gmy->sgfn], cg->fgfs[Gmz->sgfn],
                                         cg->fgfs[Lap->sgfn], 
                                         cg->fgfs[Sfx->sgfn], cg->fgfs[Sfy->sgfn], cg->fgfs[Sfz->sgfn],
                                         cg->fgfs[dtSfx->sgfn], cg->fgfs[dtSfy->sgfn], cg->fgfs[dtSfz->sgfn],
                                         cg->fgfs[Sphi->sgfn], cg->fgfs[Spi->sgfn],
                                         cg->fgfs[phi1->sgfn], cg->fgfs[trK1->sgfn],
                                         cg->fgfs[gxx1->sgfn], cg->fgfs[gxy1->sgfn], cg->fgfs[gxz1->sgfn],
                                         cg->fgfs[gyy1->sgfn], cg->fgfs[gyz1->sgfn], cg->fgfs[gzz1->sgfn],
                                         cg->fgfs[Axx1->sgfn], cg->fgfs[Axy1->sgfn], cg->fgfs[Axz1->sgfn],
                                         cg->fgfs[Ayy1->sgfn], cg->fgfs[Ayz1->sgfn], cg->fgfs[Azz1->sgfn],
                                         cg->fgfs[Gmx1->sgfn], cg->fgfs[Gmy1->sgfn], cg->fgfs[Gmz1->sgfn],
                                         cg->fgfs[Lap1->sgfn], 
                                         cg->fgfs[Sfx1->sgfn], cg->fgfs[Sfy1->sgfn], cg->fgfs[Sfz1->sgfn],
                                         cg->fgfs[dtSfx1->sgfn], cg->fgfs[dtSfy1->sgfn], cg->fgfs[dtSfz1->sgfn],
                                         cg->fgfs[Sphi1->sgfn], cg->fgfs[Spi1->sgfn],
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

            if (f_compute_rhs_bssn_escalar_ss(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
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
                                              cg->fgfs[Sphi->sgfn], cg->fgfs[Spi->sgfn],
                                              cg->fgfs[phi1->sgfn], cg->fgfs[trK1->sgfn],
                                              cg->fgfs[gxx1->sgfn], cg->fgfs[gxy1->sgfn], cg->fgfs[gxz1->sgfn],
                                              cg->fgfs[gyy1->sgfn], cg->fgfs[gyz1->sgfn], cg->fgfs[gzz1->sgfn],
                                              cg->fgfs[Axx1->sgfn], cg->fgfs[Axy1->sgfn], cg->fgfs[Axz1->sgfn],
                                              cg->fgfs[Ayy1->sgfn], cg->fgfs[Ayz1->sgfn], cg->fgfs[Azz1->sgfn],
                                              cg->fgfs[Gmx1->sgfn], cg->fgfs[Gmy1->sgfn], cg->fgfs[Gmz1->sgfn],
                                              cg->fgfs[Lap1->sgfn], 
                                              cg->fgfs[Sfx1->sgfn], cg->fgfs[Sfy1->sgfn], cg->fgfs[Sfz1->sgfn],
                                              cg->fgfs[dtSfx1->sgfn], cg->fgfs[dtSfy1->sgfn], cg->fgfs[dtSfz1->sgfn],
                                              cg->fgfs[Sphi1->sgfn], cg->fgfs[Spi1->sgfn],
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
              MyList<var> *varl0 = StateList, *varl = SynchList_pre, *varl1 = SynchList_cor, *varlrhs = RHSList; // we do not check the correspondence here
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

//================================================================================================



//================================================================================================

// 该成员函数用于计算引力波 Psi4

//================================================================================================

void bssnEScalar_class::Compute_Psi4(int lev)
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
#if (Psi4type == 0)
        // the input arguments Gamma^i_jk and R_ij do not need synch, because we do not need to derivate them
        f_getnp4scalar(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                       cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn], cg->fgfs[Sphi0->sgfn],
                       cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                       cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                       cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                       cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                       cg->fgfs[Gamxxx->sgfn], cg->fgfs[Gamxxy->sgfn], cg->fgfs[Gamxxz->sgfn],
                       cg->fgfs[Gamxyy->sgfn], cg->fgfs[Gamxyz->sgfn], cg->fgfs[Gamxzz->sgfn],
                       cg->fgfs[Gamyxx->sgfn], cg->fgfs[Gamyxy->sgfn], cg->fgfs[Gamyxz->sgfn],
                       cg->fgfs[Gamyyy->sgfn], cg->fgfs[Gamyyz->sgfn], cg->fgfs[Gamyzz->sgfn],
                       cg->fgfs[Gamzxx->sgfn], cg->fgfs[Gamzxy->sgfn], cg->fgfs[Gamzxz->sgfn],
                       cg->fgfs[Gamzyy->sgfn], cg->fgfs[Gamzyz->sgfn], cg->fgfs[Gamzzz->sgfn],
                       cg->fgfs[Rxx->sgfn], cg->fgfs[Rxy->sgfn], cg->fgfs[Rxz->sgfn], 
                       cg->fgfs[Ryy->sgfn], cg->fgfs[Ryz->sgfn], cg->fgfs[Rzz->sgfn],
                       cg->fgfs[Rpsi4->sgfn], cg->fgfs[Ipsi4->sgfn],
                       Symmetry);
#elif (Psi4type == 1)
        f_getnp4oldscalar(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                          cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn], cg->fgfs[Sphi0->sgfn],
                          cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                          cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                          cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                          cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                          cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                          cg->fgfs[Lap0->sgfn], 
                          cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                          cg->fgfs[Rpsi4->sgfn], cg->fgfs[Ipsi4->sgfn],
                          Symmetry);
#else
#error "not recognized Psi4type"
#endif
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
#if (Psi4type == 0)
          f_getnp4scalar_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
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
                            cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn], cg->fgfs[Sphi0->sgfn],
                            cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                            cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                            cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                            cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                            cg->fgfs[Gamxxx->sgfn], cg->fgfs[Gamxxy->sgfn], cg->fgfs[Gamxxz->sgfn],
                            cg->fgfs[Gamxyy->sgfn], cg->fgfs[Gamxyz->sgfn], cg->fgfs[Gamxzz->sgfn],
                            cg->fgfs[Gamyxx->sgfn], cg->fgfs[Gamyxy->sgfn], cg->fgfs[Gamyxz->sgfn],
                            cg->fgfs[Gamyyy->sgfn], cg->fgfs[Gamyyz->sgfn], cg->fgfs[Gamyzz->sgfn],
                            cg->fgfs[Gamzxx->sgfn], cg->fgfs[Gamzxy->sgfn], cg->fgfs[Gamzxz->sgfn],
                            cg->fgfs[Gamzyy->sgfn], cg->fgfs[Gamzyz->sgfn], cg->fgfs[Gamzzz->sgfn],
                            cg->fgfs[Rxx->sgfn], cg->fgfs[Rxy->sgfn], cg->fgfs[Rxz->sgfn], 
                            cg->fgfs[Ryy->sgfn], cg->fgfs[Ryz->sgfn], cg->fgfs[Rzz->sgfn],
                            cg->fgfs[Rpsi4->sgfn], cg->fgfs[Ipsi4->sgfn],
                            Symmetry, Pp->data->sst);
#elif (Psi4type == 1)
          f_getnp4oldscalar_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
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
                               cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn], cg->fgfs[Sphi0->sgfn],
                               cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                               cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                               cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                               cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                               cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                               cg->fgfs[Lap0->sgfn], 
                               cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                               cg->fgfs[Rpsi4->sgfn], cg->fgfs[Ipsi4->sgfn],
                               Symmetry, Pp->data->sst);
#else
#error "not recognized Psi4type"
#endif
        }
        if (BL == Pp->data->ble)
          break;
        BL = BL->next;
      }
      Pp = Pp->next;
    }
  }
#endif

  MyList<var> *DG_List = new MyList<var>(Rpsi4);
  DG_List->insert(Ipsi4);
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

// 该成员函数用于分析和检查标量场数据？

//================================================================================================

void bssnEScalar_class::AnalysisStuff_EScalar(int lev, double dT_lev)
{
  LastAnas += dT_lev;
  
  if (lev > 0)
  {
    cout << "AnalysisStuff_EScala only supports level 0, but lev = " << lev << endl;

    AnalysisStuff(lev, dT_lev);

    return;
  }

  if (LastAnas >= AnasTime)
  {
    MyList<var> *DG_List = new MyList<var>(Sphi0);
    double XX[3], maxs[1];

    double XXh[3], maxsh[1];
    for (int levh = GH->levels - 1; levh >= 0; levh--)
    {
      MyList<Patch> *Pp = GH->PatL[levh];

      maxsh[0] = -1; // for sure be rewriten
      while (Pp)
      {
        double XXhh[3], maxshh[1];
        Pp->data->Find_Maximum(DG_List, XXhh, maxshh);
        if (maxsh[0] < maxshh[0])
        {
          for (int i = 0; i < 3; i++)
            XXh[i] = XXhh[i];
          maxsh[0] = maxshh[0];
        }
        Pp = Pp->next;
      }

      if (levh == GH->levels - 1)
      {
        for (int i = 0; i < 3; i++)
          XX[i] = XXh[i];
        maxs[0] = maxsh[0];
      }
      else if (maxs[0] < maxsh[0])
      {
        bool fg = true;
        Pp = GH->PatL[levh + 1];

        while (Pp && fg)
        {
          if (Pp->data->Find_Point(XXh))
            fg = false; // we only take finner level
          Pp = Pp->next;
        }
        if (fg)
        {
          for (int i = 0; i < 3; i++)
            XX[i] = XXh[i];
          maxs[0] = maxsh[0];
        }
      }
    }

#ifdef WithShell
    SH->Find_Maximum(DG_List, XXh, maxsh);

    if (maxs[0] < maxsh[0])
    {
      bool fg = true;
      MyList<Patch> *Pp = GH->PatL[0];

      while (Pp && fg)
      {
        if (Pp->data->Find_Point(XXh))
          fg = false;
        Pp = Pp->next;
      }
      if (fg)
      {
        for (int i = 0; i < 3; i++)
          XX[i] = XXh[i];
        maxs[0] = maxsh[0];
      }
    }
#endif

    double RD[4];
    for (int i = 0; i < 3; i++)
      RD[i] = XX[i];
    RD[3] = maxs[0];
    MaxScalar_Monitor->writefile(PhysTime, 4, RD);

    DG_List->clearList();
  }

  AnalysisStuff(lev, dT_lev); // LastAnas need and only need control here
  
  LastAnas = 0;
}

//================================================================================================



//================================================================================================

// 该成员函数用于对约束数据进行插值

//================================================================================================

void bssnEScalar_class::Interp_Constraint()
{
  // we do not support a_lev != 0 yet.
  if (a_lev > 0)
    return;

  for (int lev = 0; lev < GH->levels; lev++)
  {
    // make sure the data consistent for higher levels
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
            if (lev > 0)
              f_compute_rhs_bssn_escalar(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                                         cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                                         cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                                         cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                                         cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                                         cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                                         cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                                         cg->fgfs[Lap0->sgfn], 
                                         cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                                         cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn],
                                         cg->fgfs[Sphi0->sgfn], cg->fgfs[Spi0->sgfn],
                                         cg->fgfs[phi_rhs->sgfn], cg->fgfs[trK_rhs->sgfn],
                                         cg->fgfs[gxx_rhs->sgfn], cg->fgfs[gxy_rhs->sgfn], cg->fgfs[gxz_rhs->sgfn],
                                         cg->fgfs[gyy_rhs->sgfn], cg->fgfs[gyz_rhs->sgfn], cg->fgfs[gzz_rhs->sgfn],
                                         cg->fgfs[Axx_rhs->sgfn], cg->fgfs[Axy_rhs->sgfn], cg->fgfs[Axz_rhs->sgfn],
                                         cg->fgfs[Ayy_rhs->sgfn], cg->fgfs[Ayz_rhs->sgfn], cg->fgfs[Azz_rhs->sgfn],
                                         cg->fgfs[Gmx_rhs->sgfn], cg->fgfs[Gmy_rhs->sgfn], cg->fgfs[Gmz_rhs->sgfn],
                                         cg->fgfs[Lap_rhs->sgfn], 
                                         cg->fgfs[Sfx_rhs->sgfn], cg->fgfs[Sfy_rhs->sgfn], cg->fgfs[Sfz_rhs->sgfn],
                                         cg->fgfs[dtSfx_rhs->sgfn], cg->fgfs[dtSfy_rhs->sgfn], cg->fgfs[dtSfz_rhs->sgfn],
                                         cg->fgfs[Sphi_rhs->sgfn], cg->fgfs[Spi_rhs->sgfn],
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
            f_compute_constraint_fr(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                    cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn], 
                                    cg->fgfs[rho->sgfn], cg->fgfs[Sphi0->sgfn],
                                    cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                                    cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                                    cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                                    cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                                    cg->fgfs[Rxx->sgfn], cg->fgfs[Rxy->sgfn], cg->fgfs[Rxz->sgfn], 
                                    cg->fgfs[Ryy->sgfn], cg->fgfs[Ryz->sgfn], cg->fgfs[Rzz->sgfn],
                                    cg->fgfs[Sxx->sgfn], cg->fgfs[Sxy->sgfn], cg->fgfs[Sxz->sgfn], 
                                    cg->fgfs[Syy->sgfn], cg->fgfs[Syz->sgfn], cg->fgfs[Szz->sgfn],
                                    cg->fgfs[Cons_fR->sgfn]);
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
  // ShellPatch part
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
          f_compute_constraint_fr(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                  cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn], 
                                  cg->fgfs[rho->sgfn], cg->fgfs[Sphi0->sgfn],
                                  cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                                  cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                                  cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                                  cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                                  cg->fgfs[Rxx->sgfn], cg->fgfs[Rxy->sgfn], cg->fgfs[Rxz->sgfn], 
                                  cg->fgfs[Ryy->sgfn], cg->fgfs[Ryz->sgfn], cg->fgfs[Rzz->sgfn],
                                  cg->fgfs[Sxx->sgfn], cg->fgfs[Sxy->sgfn], cg->fgfs[Sxz->sgfn], 
                                  cg->fgfs[Syy->sgfn], cg->fgfs[Syz->sgfn], cg->fgfs[Szz->sgfn],
                                  cg->fgfs[Cons_fR->sgfn]);
        }
        if (BL == Pp->data->ble)
          break;
        BL = BL->next;
      }
      Pp = Pp->next;
    }
  }

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
  outfile << "#  corrdinate, H_Res, Px_Res, Py_Res, Pz_Res, Gx_Res, Gy_Res, Gz_Res, fR_Res, ...." << endl;
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

// 该成员函数用于计算和输出约束违反

//================================================================================================

void bssnEScalar_class::Constraint_Out()
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
              if (lev > 0)
                f_compute_rhs_bssn_escalar(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                                           cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                                           cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                                           cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                                           cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                                           cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                                           cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                                           cg->fgfs[Lap0->sgfn], 
                                           cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                                           cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn],
                                           cg->fgfs[Sphi0->sgfn], cg->fgfs[Spi0->sgfn],
                                           cg->fgfs[phi_rhs->sgfn], cg->fgfs[trK_rhs->sgfn],
                                           cg->fgfs[gxx_rhs->sgfn], cg->fgfs[gxy_rhs->sgfn], cg->fgfs[gxz_rhs->sgfn],
                                           cg->fgfs[gyy_rhs->sgfn], cg->fgfs[gyz_rhs->sgfn], cg->fgfs[gzz_rhs->sgfn],
                                           cg->fgfs[Axx_rhs->sgfn], cg->fgfs[Axy_rhs->sgfn], cg->fgfs[Axz_rhs->sgfn],
                                           cg->fgfs[Ayy_rhs->sgfn], cg->fgfs[Ayz_rhs->sgfn], cg->fgfs[Azz_rhs->sgfn],
                                           cg->fgfs[Gmx_rhs->sgfn], cg->fgfs[Gmy_rhs->sgfn], cg->fgfs[Gmz_rhs->sgfn],
                                           cg->fgfs[Lap_rhs->sgfn], 
                                           cg->fgfs[Sfx_rhs->sgfn], cg->fgfs[Sfy_rhs->sgfn], cg->fgfs[Sfz_rhs->sgfn],
                                           cg->fgfs[dtSfx_rhs->sgfn], cg->fgfs[dtSfy_rhs->sgfn], cg->fgfs[dtSfz_rhs->sgfn],
                                           cg->fgfs[Sphi_rhs->sgfn], cg->fgfs[Spi_rhs->sgfn],
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
              f_compute_constraint_fr(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                      cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn], 
                                      cg->fgfs[rho->sgfn], cg->fgfs[Sphi0->sgfn],
                                      cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                                      cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                                      cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                                      cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                                      cg->fgfs[Rxx->sgfn], cg->fgfs[Rxy->sgfn], cg->fgfs[Rxz->sgfn], 
                                      cg->fgfs[Ryy->sgfn], cg->fgfs[Ryz->sgfn], cg->fgfs[Rzz->sgfn],
                                      cg->fgfs[Sxx->sgfn], cg->fgfs[Sxy->sgfn], cg->fgfs[Sxz->sgfn], 
                                      cg->fgfs[Syy->sgfn], cg->fgfs[Syz->sgfn], cg->fgfs[Szz->sgfn],
                                      cg->fgfs[Cons_fR->sgfn]);
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
    // ShellPatch part
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
            f_compute_constraint_fr(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                    cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn], 
                                    cg->fgfs[rho->sgfn], cg->fgfs[Sphi0->sgfn],
                                    cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                                    cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                                    cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                                    cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                                    cg->fgfs[Rxx->sgfn], cg->fgfs[Rxy->sgfn], cg->fgfs[Rxz->sgfn], 
                                    cg->fgfs[Ryy->sgfn], cg->fgfs[Ryz->sgfn], cg->fgfs[Rzz->sgfn],
                                    cg->fgfs[Sxx->sgfn], cg->fgfs[Sxy->sgfn], cg->fgfs[Sxz->sgfn], 
                                    cg->fgfs[Syy->sgfn], cg->fgfs[Syz->sgfn], cg->fgfs[Szz->sgfn],
                                    cg->fgfs[Cons_fR->sgfn]);
          }
          if (BL == Pp->data->ble)
            break;
          BL = BL->next;
        }
        Pp = Pp->next;
      }
    }

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
    ConV[7] = SH->L2Norm(Cons_fR);
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
      ConV[7] = Parallel::L2Norm(GH->PatL[levi]->data, Cons_fR);
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

// 读入标量张量理论的参数
// 小曲修改
// 一次读入多个值
// 原来函数是一个一个读入

//================================================================================================

extern "C"
{

#ifdef fortran1
  void set_escalar_parameter
#endif
#ifdef fortran2
      void SET_ESCALAR_PARAMETER
#endif
#ifdef fortran3
      void set_escalar_parameter_
#endif
  
  (double &a2, double &phi0, double &r0, double &sigma0, double &l2)
  {

    static bool file_status = true;  
    // 使用静态布尔类型控制，避免重复阅读参数文件？
    // 这种类型的变量似乎是共享的，一旦阅读了一此后，其它进程会记住它的状态
    // 参数文件读取完之后，file_status 会自动转换为 false

    static double aa2;
    static double ll2;
    static double pphi0;
    static double rr0;
    static double ssigma0;
    
    int myrank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    
    // read parameter from file
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
    
  if (file_status)
  {
    const int LEN = 256;
    char pline[LEN];
    string str, sgrp, skey, sval;
    int sind;
    ifstream inf(pname, ifstream::in);
    if (!inf.good() && myrank == 0)
    {
      cout << "Can not open parameter file " << pname << " for inputing information of EScalar" << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    for (int i = 1; inf.good(); i++)
    {
      inf.getline(pline, LEN);
      str = pline;

      int status = misc::parse_parts(str, sgrp, skey, sval, sind);
      if (status == -1)
      {
        cout << "error reading parameter file " << pname << " in line " << i << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
      else if (status == 0)
        continue;

      if (sgrp == "FR" && skey == "a2")
        aa2 = atof(sval.c_str());
      else if (sgrp == "FR" && skey == "l2")
        ll2 = atof(sval.c_str());
      else if (sgrp == "FR" && skey == "phi0")
        pphi0 = atof(sval.c_str());
      else if (sgrp == "FR" && skey == "r0")
        rr0 = atof(sval.c_str());
      else if (sgrp == "FR" && skey == "sigma0")
        ssigma0 = atof(sval.c_str());
    }
    
    inf.close(); // if not closed, it will fail when you try to open it next time.
    
    // 参数文件读取完之后，file_status 转换为 false
    file_status = false;
    
    int myrank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank == 0)
    {
      cout <<                                         endl;
      cout << " you have set a2     = " << aa2     << endl;
      cout << " you have set l2     = " << ll2     << endl;
      cout << " you have set phi0   = " << pphi0   << endl;
      cout << " you have set r0     = " << rr0     << endl;
      cout << " you have set sigma0 = " << ssigma0 << endl;
      cout <<                                         endl;
    }
  }

  a2     = aa2;
  phi0   = pphi0;
  r0     = rr0;
  sigma0 = ssigma0;
  l2     = ll2;
  }
}


// 原函数：一个一个完成读入，太麻烦了

extern "C"
{

#ifdef fortran1
  void seta2
#endif
#ifdef fortran2
      void SETA2
#endif
#ifdef fortran3
      void
      seta2_
#endif
      (double &a2)
  {
    static bool fga2 = true;
    static double aa2;

    if (fga2)
    {
      char s[1000], *t;
      FILE *fp;

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
      fp = fopen(pname, "r");
      if (!fp)
      {
        cout << "could not open " << pname << " for reading a2" << endl;
      }
      while (fgets(s, 1000, fp))
      {
        t = strstr(s, "FR::a2 ");
        if (t == s)
        {
          sscanf(s + 8, "%lf", &aa2);
          break;
        }
      }

      fclose(fp); // if not closed, it will fail when you try to open it next time.
      fga2 = false;

      int myrank = 0;
      MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
      if (myrank == 0)
      {
        printf("you have set a2 = %0.4lg\n", aa2);
      }
    }

    a2 = aa2;
  }
}

extern "C"
{

#ifdef fortran1
  void setphi0
#endif
#ifdef fortran2
      void SETPHI0
#endif
#ifdef fortran3
      void
      setphi0_
#endif
      (double &phi0)
  {
    static bool fgphi0 = true;
    static double pphi0;

    if (fgphi0)
    {
      char s[1000], *t;
      FILE *fp;

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
      fp = fopen(pname, "r");
      if (!fp)
      {
        cout << "could not open " << pname << " for reading phi0" << endl;
      }
      while (fgets(s, 1000, fp))
      {
        t = strstr(s, "FR::phi0 ");
        if (t == s)
        {
          sscanf(s + 10, "%lf", &pphi0);
          break;
        }
      }

      fclose(fp); // if not closed, it will fail when you try to open it next time.
      fgphi0 = false;

      int myrank = 0;
      MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
      if (myrank == 0)
      {
        printf("you have set phi0 = %0.4lg\n", pphi0);
      }
    }

    phi0 = pphi0;
  }
}

//================================================================================================

