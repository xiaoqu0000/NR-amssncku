
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
#include "fmisc.h"
#include "Parallel.h"
#include "scalar_class.h"
#include "scalar_rhs.h"
#include "initial_scalar.h"
#include "rungekutta4_rout.h"
#include "sommerfeld_rout.h"
#include "shellfunctions.h"
#include "parameters.h"

scalar_class::scalar_class(double Couranti, double StartTimei, double TotalTimei, double DumpTimei, double CheckTimei, double AnasTimei,
                           int Symmetryi, int checkruni, char *checkfilenamei, double numepssi, double numepsbi,
                           int a_levi) : Courant(Couranti), StartTime(StartTimei), TotalTime(TotalTimei), DumpTime(DumpTimei), CheckTime(CheckTimei), AnasTime(AnasTimei),
                                         Symmetry(Symmetryi), checkrun(checkruni), numepss(numepssi), numepsb(numepsbi),
                                         a_lev(a_levi)
{
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if (checkrun)
  {
  }
  else
  {
    PhysTime = StartTime;
  }
  // setup Monitors
  {
    stringstream a_stream;
    a_stream.setf(ios::left);
    a_stream << "# Error log information";
    ErrorMonitor = new monitor("Error.log", myrank, a_stream.str());
  }

  trfls = 0;
  // read parameter from file
  {
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
    const int LEN = 256;
    char pline[LEN];
    string str, sgrp, skey, sval;
    int sind;
    ifstream inf(filename, ifstream::in);
    if (!inf.good() && myrank == 0)
    {
      if (ErrorMonitor->outfile)
        ErrorMonitor->outfile << "Can not open parameter file " << filename << " for inputing information of black holes" << endl;
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
          ErrorMonitor->outfile << "error reading parameter file " << filename << " in line " << i << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
      else if (status == 0)
        continue;

      if (sgrp == "SCALAR" && skey == "time refinement start from level")
        trfls = atoi(sval.c_str());
    }
    inf.close();
  }
  // echo read-in information
  if (myrank == 0)
  {
    cout << "time refinement start from level #" << trfls << endl;
  }

  strcpy(checkfilename, checkfilenamei);

  int ngfs = 0;
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

  if (myrank == 0)
    cout << "you have setted " << ngfs << " grid functions." << endl;

  OldStateList = new MyList<var>(Sphio);
  OldStateList->insert(Spio);

  StateList = new MyList<var>(Sphi0);
  StateList->insert(Spi0);

  RHSList = new MyList<var>(Sphi_rhs);
  RHSList->insert(Spi_rhs);

  SynchList_pre = new MyList<var>(Sphi);
  SynchList_pre->insert(Spi);

  SynchList_cor = new MyList<var>(Sphi1);
  SynchList_cor->insert(Spi1);

  DumpList = new MyList<var>(Sphi0);
  DumpList->insert(Spi0);

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
  GH->compose_cgh(nprocs);
#ifdef WithShell
  SH = new ShellPatch(0, ngfs, pname, Symmetry, myrank, ErrorMonitor);
  SH->matchcheck(GH->PatL[0]);
  SH->compose_sh(nprocs);
  //   SH->compose_shr(nprocs);  //sh is faster than shr
  SH->setupcordtrans();
  SH->Dump_xyz(0, 0, 1);
  SH->setupintintstuff(nprocs, GH->PatL[0], Symmetry);
#else
  SH = 0;
#endif

  double h = GH->PatL[0]->data->blb->data->getdX(0);
  for (int i = 1; i < dim; i++)
    h = Mymin(h, GH->PatL[0]->data->blb->data->getdX(i));
  dT = Courant * h;
}
scalar_class::~scalar_class()
{
  StateList->clearList();
  RHSList->clearList();
  OldStateList->clearList();
  SynchList_pre->clearList();
  SynchList_cor->clearList();
  DumpList->clearList();

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

  delete GH;
#ifdef WithShell
  delete SH;
#endif

  delete ErrorMonitor;
}
void scalar_class::Setup_Initial_Data()
{
  if (checkrun)
  {
  }
  else
  {
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
    double R0, WD, A;
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
          ErrorMonitor->outfile << "Can not open parameter file " << filename << " for inputing information of black holes" << endl;
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
            ErrorMonitor->outfile << "error reading parameter file " << filename << " in line " << i << endl;
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else if (status == 0)
          continue;

        if (sgrp == "SCALAR")
        {
          if (skey == "center of Gauss")
            R0 = atof(sval.c_str());
          else if (skey == "width of Gauss")
            WD = atof(sval.c_str());
          else if (skey == "amplitude of Gauss")
            A = atof(sval.c_str());
        }
      }
      inf.close();
    }
    // echo read-in information
    if (myrank == 0)
    {
      cout << "Setup initial scalar with Gauss profile " << A << "*exp[-(r-" << R0 << ")^2/2/" << WD << "^2]" << endl;
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
            f_get_initial_scalar(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                 cg->fgfs[Sphi0->sgfn], cg->fgfs[Spi0->sgfn], R0, WD, A);
          }
          if (BL == Pp->data->ble)
            break;
          BL = BL->next;
        }
        Pp = Pp->next;
      }
    }

    for (int lev = 0; lev < GH->levels; lev++)
      Parallel::Dump_Data(GH->PatL[lev], DumpList, 0, PhysTime, dT);
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
          f_get_initial_scalar_sh(cg->shape, cg->fgfs[Pp->data->fngfs + ShellPatch::gx], cg->fgfs[Pp->data->fngfs + ShellPatch::gy],
                                  cg->fgfs[Pp->data->fngfs + ShellPatch::gz],
                                  cg->fgfs[Sphi0->sgfn], cg->fgfs[Spi0->sgfn], R0, WD, A);
        }
        if (BL == Pp->data->ble)
          break;
        BL = BL->next;
      }
      Pp = Pp->next;
    }
// dump read_in initial data
//   SH->Synch(GH->PatL[0],StateList,Symmetry);
//   for(int lev=0;lev<GH->levels;lev++) Parallel::Dump_Data(GH->PatL[lev],StateList,0,PhysTime,dT);
//   SH->Dump_Data(StateList,0,PhysTime,dT);
//   exit(0);
#endif
  }
}
void scalar_class::Evolve(int Steps)
{
  clock_t prev_clock, curr_clock;
  double LastDump = 0.0, LastCheck = 0.0;
  LastAnas = 0;

  double dT_mon = dT * pow(0.5, Mymax(0, trfls));

  for (int ncount = 1; ncount < Steps + 1; ncount++)
  {
    if (myrank == 0)
      curr_clock = clock();
    RecursiveStep(0);

    LastDump += dT_mon;
    LastCheck += dT_mon;

    if (LastDump >= DumpTime)
    {
      for (int lev = 0; lev < GH->levels; lev++)
        Parallel::Dump_Data(GH->PatL[lev], DumpList, 0, PhysTime, dT_mon);
#ifdef WithShell
      SH->Dump_Data(DumpList, 0, PhysTime, dT_mon);
#endif
      LastDump = 0;
    }
    if (myrank == 0)
    {
      prev_clock = curr_clock;
      curr_clock = clock();
      cout << " Timestep # " << ncount << ": integrating to time: " << PhysTime
           << " Computer used " << (double)(curr_clock - prev_clock) / ((double)CLOCKS_PER_SEC) << " seconds! " << endl;
    }
    if (PhysTime >= TotalTime)
      break;
  }
}
void scalar_class::RecursiveStep(int lev)
{
  int NoIterations = 1, YN;
  if (lev <= trfls)
    NoIterations = 1;
  else
    NoIterations = 2;

  for (int i = 0; i < NoIterations; i++)
  {
    //     if(myrank==0) cout<<"level now = "<<lev<<" NoIteration = "<<i<<endl;
    YN = (i == NoIterations - 1) ? 1 : 0; // 1: same time level for coarse level and fine level
    Step(lev, YN);
    if (lev < GH->levels - 1)
    {
      int lf = lev + 1;
      RecursiveStep(lf);
    }
    else
      PhysTime += dT * pow(0.5, lev);

    Parallel::Dump_Data(GH->PatL[lev], DumpList, 0, PhysTime, dT * pow(0.5, lev));
  }
}
#if 1
void scalar_class::Step(int lev, int YN)
{
  double dT_lev = dT * pow(0.5, Mymax(lev, trfls));
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
        if (f_compute_rhs_scalar(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                                 cg->fgfs[Sphi0->sgfn], cg->fgfs[Spi0->sgfn],
                                 cg->fgfs[Sphi_rhs->sgfn], cg->fgfs[Spi_rhs->sgfn],
                                 Symmetry, lev, ndeps))
        {
          cout << "find NaN in domain: (" << cg->bbox[0] << ":" << cg->bbox[3] << "," << cg->bbox[1] << ":" << cg->bbox[4] << ","
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
                                   Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                                   cg->fgfs[varlrhs->data->sgfn],
                                   cg->fgfs[varl0->data->sgfn], varl0->data->propspeed, varl0->data->SoA,
                                   Symmetry);

#endif
            f_rungekutta4_rout(cg->shape, dT_lev, cg->fgfs[varl0->data->sgfn], cg->fgfs[varl->data->sgfn], cg->fgfs[varlrhs->data->sgfn],
                               iter_count);
#ifndef WithShell
            if (lev > 0) // fix BD point
#endif
              f_sommerfeld_rout(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                                dT_lev, cg->fgfs[Sphi0->sgfn],
                                cg->fgfs[Spi0->sgfn], cg->fgfs[varl0->data->sgfn], cg->fgfs[varl->data->sgfn], varl0->data->SoA,
                                Symmetry, cor);

            varl0 = varl0->next;
            varl = varl->next;
            varlrhs = varlrhs->next;
          }
        }
      }
      if (BP == Pp->data->ble)
        break;
      BP = BP->next;
    }
    Pp = Pp->next;
  }

  //  if(lev==1) Parallel::Dump_Data(GH->PatL[lev],RHSList,0,PhysTime,dT_lev);
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
        ErrorMonitor->outfile << "find NaN in state variables at t = " << PhysTime << ", lev = " << lev << endl;
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
          if (f_compute_rhs_scalar_ss(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                                      cg->fgfs[fngfs + ShellPatch::gx], cg->fgfs[fngfs + ShellPatch::gy], cg->fgfs[fngfs + ShellPatch::gz],
                                      cg->fgfs[fngfs + ShellPatch::drhodx], cg->fgfs[fngfs + ShellPatch::drhody], cg->fgfs[fngfs + ShellPatch::drhodz],
                                      cg->fgfs[fngfs + ShellPatch::dsigmadx], cg->fgfs[fngfs + ShellPatch::dsigmady], cg->fgfs[fngfs + ShellPatch::dsigmadz],
                                      cg->fgfs[fngfs + ShellPatch::dRdx], cg->fgfs[fngfs + ShellPatch::dRdy], cg->fgfs[fngfs + ShellPatch::dRdz],
                                      cg->fgfs[fngfs + ShellPatch::drhodxx], cg->fgfs[fngfs + ShellPatch::drhodxy], cg->fgfs[fngfs + ShellPatch::drhodxz],
                                      cg->fgfs[fngfs + ShellPatch::drhodyy], cg->fgfs[fngfs + ShellPatch::drhodyz], cg->fgfs[fngfs + ShellPatch::drhodzz],
                                      cg->fgfs[fngfs + ShellPatch::dsigmadxx], cg->fgfs[fngfs + ShellPatch::dsigmadxy], cg->fgfs[fngfs + ShellPatch::dsigmadxz],
                                      cg->fgfs[fngfs + ShellPatch::dsigmadyy], cg->fgfs[fngfs + ShellPatch::dsigmadyz], cg->fgfs[fngfs + ShellPatch::dsigmadzz],
                                      cg->fgfs[fngfs + ShellPatch::dRdxx], cg->fgfs[fngfs + ShellPatch::dRdxy], cg->fgfs[fngfs + ShellPatch::dRdxz],
                                      cg->fgfs[fngfs + ShellPatch::dRdyy], cg->fgfs[fngfs + ShellPatch::dRdyz], cg->fgfs[fngfs + ShellPatch::dRdzz],
                                      cg->fgfs[Sphi0->sgfn], cg->fgfs[Spi0->sgfn],
                                      cg->fgfs[Sphi_rhs->sgfn], cg->fgfs[Spi_rhs->sgfn],
                                      Symmetry, lev, ndeps, sPp->data->sst))
          {
            cout << "find NaN in Shell domain: sst = " << sPp->data->sst << ", (" << cg->bbox[0] << ":" << cg->bbox[3] << ","
                 << cg->bbox[1] << ":" << cg->bbox[4] << "," << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
            ERROR = 1;
          }

          // rk4 substep and boundary
          {
            MyList<var> *varl0 = StateList, *varl = SynchList_pre, *varlrhs = RHSList; // we do not check the correspondence here
            while (varl0)
            {
              // sommerfeld indeed for outter boudary while fix BD for inner boundary
              f_sommerfeld_routbam_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                      sPp->data->bbox[0], sPp->data->bbox[1], sPp->data->bbox[2], sPp->data->bbox[3], sPp->data->bbox[4], sPp->data->bbox[5],
                                      cg->fgfs[varlrhs->data->sgfn],
                                      cg->fgfs[varl0->data->sgfn], varl0->data->propspeed, varl0->data->SoA,
                                      Symmetry);

              f_rungekutta4_rout(cg->shape, dT_lev, cg->fgfs[varl0->data->sgfn], cg->fgfs[varl->data->sgfn], cg->fgfs[varlrhs->data->sgfn],
                                 iter_count);

              varl0 = varl0->next;
              varl = varl->next;
              varlrhs = varlrhs->next;
            }
          }
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
      cout << " Shell stuff synchronization used " << (double)(curr_clock - prev_clock) / ((double)CLOCKS_PER_SEC) << " seconds! " << endl;
    }
  }
#endif
  // data analysis part
  // Warning NOTE: the variables1 are used as temp storege room
  if (lev == a_lev)
  {
    if (LastAnas >= AnasTime)
    {

      LastAnas = 0;
    }
    LastAnas += dT_lev;
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
          if (f_compute_rhs_scalar(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                                   cg->fgfs[Sphi->sgfn], cg->fgfs[Spi->sgfn],
                                   cg->fgfs[Sphi1->sgfn], cg->fgfs[Spi1->sgfn],
                                   Symmetry, lev, ndeps))
          {
            cout << "find NaN in domain: (" << cg->bbox[0] << ":" << cg->bbox[3] << "," << cg->bbox[1] << ":" << cg->bbox[4] << ","
                 << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
            ERROR = 1;
          }
          // rk4 substep and boundary
          {
            MyList<var> *varl0 = StateList, *varl = SynchList_pre, *varl1 = SynchList_cor, *varlrhs = RHSList; // we do not check the correspondence here
            while (varl0)
            {
#ifndef WithShell
              if (lev == 0) // sommerfeld indeed
                f_sommerfeld_routbam(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                     Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                                     cg->fgfs[varl1->data->sgfn],
                                     cg->fgfs[varl->data->sgfn], varl0->data->propspeed, varl0->data->SoA,
                                     Symmetry);
#endif
              f_rungekutta4_rout(cg->shape, dT_lev, cg->fgfs[varl0->data->sgfn], cg->fgfs[varl1->data->sgfn], cg->fgfs[varlrhs->data->sgfn],
                                 iter_count);

#ifndef WithShell
              if (lev > 0) // fix BD point
#endif
                f_sommerfeld_rout(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                  Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                                  dT_lev, cg->fgfs[Sphi0->sgfn],
                                  cg->fgfs[Spi0->sgfn], cg->fgfs[varl0->data->sgfn], cg->fgfs[varl1->data->sgfn], varl0->data->SoA,
                                  Symmetry, cor);

              varl0 = varl0->next;
              varl = varl->next;
              varl1 = varl1->next;
              varlrhs = varlrhs->next;
            }
          }
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
          ErrorMonitor->outfile << "find NaN in RK4 substep#" << iter_count << " variables at t = " << PhysTime << ", lev = " << lev << endl;
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
            if (f_compute_rhs_scalar_ss(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                                        cg->fgfs[fngfs + ShellPatch::gx], cg->fgfs[fngfs + ShellPatch::gy], cg->fgfs[fngfs + ShellPatch::gz],
                                        cg->fgfs[fngfs + ShellPatch::drhodx], cg->fgfs[fngfs + ShellPatch::drhody], cg->fgfs[fngfs + ShellPatch::drhodz],
                                        cg->fgfs[fngfs + ShellPatch::dsigmadx], cg->fgfs[fngfs + ShellPatch::dsigmady], cg->fgfs[fngfs + ShellPatch::dsigmadz],
                                        cg->fgfs[fngfs + ShellPatch::dRdx], cg->fgfs[fngfs + ShellPatch::dRdy], cg->fgfs[fngfs + ShellPatch::dRdz],
                                        cg->fgfs[fngfs + ShellPatch::drhodxx], cg->fgfs[fngfs + ShellPatch::drhodxy], cg->fgfs[fngfs + ShellPatch::drhodxz],
                                        cg->fgfs[fngfs + ShellPatch::drhodyy], cg->fgfs[fngfs + ShellPatch::drhodyz], cg->fgfs[fngfs + ShellPatch::drhodzz],
                                        cg->fgfs[fngfs + ShellPatch::dsigmadxx], cg->fgfs[fngfs + ShellPatch::dsigmadxy], cg->fgfs[fngfs + ShellPatch::dsigmadxz],
                                        cg->fgfs[fngfs + ShellPatch::dsigmadyy], cg->fgfs[fngfs + ShellPatch::dsigmadyz], cg->fgfs[fngfs + ShellPatch::dsigmadzz],
                                        cg->fgfs[fngfs + ShellPatch::dRdxx], cg->fgfs[fngfs + ShellPatch::dRdxy], cg->fgfs[fngfs + ShellPatch::dRdxz],
                                        cg->fgfs[fngfs + ShellPatch::dRdyy], cg->fgfs[fngfs + ShellPatch::dRdyz], cg->fgfs[fngfs + ShellPatch::dRdzz],
                                        cg->fgfs[Sphi->sgfn], cg->fgfs[Spi->sgfn],
                                        cg->fgfs[Sphi1->sgfn], cg->fgfs[Spi1->sgfn],
                                        Symmetry, lev, ndeps, sPp->data->sst))
            {
              cout << "find NaN in Shell domain: sst = " << sPp->data->sst << ", (" << cg->bbox[0] << ":" << cg->bbox[3] << ","
                   << cg->bbox[1] << ":" << cg->bbox[4] << "," << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
              ERROR = 1;
            }
            // rk4 substep and boundary
            {
              MyList<var> *varl0 = StateList, *varl = SynchList_pre, *varl1 = SynchList_cor, *varlrhs = RHSList; // we do not check the correspondence here
              while (varl0)
              {
                f_sommerfeld_routbam_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                        sPp->data->bbox[0], sPp->data->bbox[1], sPp->data->bbox[2], sPp->data->bbox[3], sPp->data->bbox[4], sPp->data->bbox[5],
                                        cg->fgfs[varl1->data->sgfn],
                                        cg->fgfs[varl->data->sgfn], varl0->data->propspeed, varl0->data->SoA,
                                        Symmetry);

                f_rungekutta4_rout(cg->shape, dT_lev, cg->fgfs[varl0->data->sgfn], cg->fgfs[varl1->data->sgfn], cg->fgfs[varlrhs->data->sgfn],
                                   iter_count);

                varl0 = varl0->next;
                varl = varl->next;
                varl1 = varl1->next;
                varlrhs = varlrhs->next;
              }
            }
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
          ErrorMonitor->outfile << "find NaN on Shell Patches in RK4 substep#" << iter_count << " variables at t = " << PhysTime << endl;
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
        cout << " Shell stuff synchronization used " << (double)(curr_clock - prev_clock) / ((double)CLOCKS_PER_SEC) << " seconds! " << endl;
      }
    }
#endif
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
    }
  }
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
      cout << " CS_Inter used " << (double)(curr_clock - prev_clock) / ((double)CLOCKS_PER_SEC) << " seconds! " << endl;
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
}
#else
// for check, using Euler method
void scalar_class::Step(int lev, int YN)
{
  double dT_lev = dT * pow(0.5, Mymax(lev, trfls));
  bool BB = fgt(PhysTime, StartTime, dT_lev / 2);
  double ndeps = numepss;
  if (lev < GH->movls)
    ndeps = numepsb;
  double TRK4 = PhysTime;
  int iter_count = 0; // count RK4 substeps
  int pre = 0, cor = 1;
  int ERROR = 0;

  MyList<ss_patch> *sPp;

  MyList<Patch> *Pp = GH->PatL[lev];
  while (Pp)
  {
    MyList<Block> *BP = Pp->data->blb;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        if (f_compute_rhs_scalar(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                                 cg->fgfs[Sphi0->sgfn], cg->fgfs[Spi0->sgfn],
                                 cg->fgfs[Sphi_rhs->sgfn], cg->fgfs[Spi_rhs->sgfn],
                                 Symmetry, lev, ndeps))
        {
          cout << "find NaN in domain: (" << cg->bbox[0] << ":" << cg->bbox[3] << "," << cg->bbox[1] << ":" << cg->bbox[4] << ","
               << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
          ERROR = 1;
        }

        // rk4 substep and boundary
        {
          MyList<var> *varl0 = StateList, *varl1 = SynchList_cor, *varlrhs = RHSList; // we do not check the correspondence here
          while (varl0)
          {
#ifndef WithShell
            if (lev == 0) // sommerfeld indeed
              f_sommerfeld_routbam(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                   Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                                   cg->fgfs[varl1->data->sgfn],
                                   cg->fgfs[varl0->data->sgfn], varl0->data->propspeed, varl0->data->SoA,
                                   Symmetry);
#endif
            f_euler_rout(cg->shape, dT_lev, cg->fgfs[varl0->data->sgfn], cg->fgfs[varl1->data->sgfn], cg->fgfs[varlrhs->data->sgfn]);

#ifndef WithShell
            if (lev > 0) // fix BD point
#endif
              f_sommerfeld_rout(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                                dT_lev, cg->fgfs[Sphi0->sgfn],
                                cg->fgfs[Spi0->sgfn], cg->fgfs[varl0->data->sgfn], cg->fgfs[varl1->data->sgfn], varl0->data->SoA,
                                Symmetry, cor);

            varl0 = varl0->next;
            varl1 = varl1->next;
            varlrhs = varlrhs->next;
          }
        }
      }
      if (BP == Pp->data->ble)
        break;
      BP = BP->next;
    }
    Pp = Pp->next;
  }

  //  if(lev==1) Parallel::Dump_Data(GH->PatL[lev],RHSList,0,PhysTime,dT_lev);
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
        ErrorMonitor->outfile << "find NaN in state variables at t = " << PhysTime << ", lev = " << lev << endl;
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
          if (f_compute_rhs_scalar_ss(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                                      cg->fgfs[fngfs + ShellPatch::gx], cg->fgfs[fngfs + ShellPatch::gy], cg->fgfs[fngfs + ShellPatch::gz],
                                      cg->fgfs[fngfs + ShellPatch::drhodx], cg->fgfs[fngfs + ShellPatch::drhody], cg->fgfs[fngfs + ShellPatch::drhodz],
                                      cg->fgfs[fngfs + ShellPatch::dsigmadx], cg->fgfs[fngfs + ShellPatch::dsigmady], cg->fgfs[fngfs + ShellPatch::dsigmadz],
                                      cg->fgfs[fngfs + ShellPatch::dRdx], cg->fgfs[fngfs + ShellPatch::dRdy], cg->fgfs[fngfs + ShellPatch::dRdz],
                                      cg->fgfs[fngfs + ShellPatch::drhodxx], cg->fgfs[fngfs + ShellPatch::drhodxy], cg->fgfs[fngfs + ShellPatch::drhodxz],
                                      cg->fgfs[fngfs + ShellPatch::drhodyy], cg->fgfs[fngfs + ShellPatch::drhodyz], cg->fgfs[fngfs + ShellPatch::drhodzz],
                                      cg->fgfs[fngfs + ShellPatch::dsigmadxx], cg->fgfs[fngfs + ShellPatch::dsigmadxy], cg->fgfs[fngfs + ShellPatch::dsigmadxz],
                                      cg->fgfs[fngfs + ShellPatch::dsigmadyy], cg->fgfs[fngfs + ShellPatch::dsigmadyz], cg->fgfs[fngfs + ShellPatch::dsigmadzz],
                                      cg->fgfs[fngfs + ShellPatch::dRdxx], cg->fgfs[fngfs + ShellPatch::dRdxy], cg->fgfs[fngfs + ShellPatch::dRdxz],
                                      cg->fgfs[fngfs + ShellPatch::dRdyy], cg->fgfs[fngfs + ShellPatch::dRdyz], cg->fgfs[fngfs + ShellPatch::dRdzz],
                                      cg->fgfs[Sphi0->sgfn], cg->fgfs[Spi0->sgfn],
                                      cg->fgfs[Sphi_rhs->sgfn], cg->fgfs[Spi_rhs->sgfn],
                                      Symmetry, lev, ndeps, sPp->data->sst))
          {
            cout << "find NaN in Shell domain: sst = " << sPp->data->sst << ", (" << cg->bbox[0] << ":" << cg->bbox[3] << ","
                 << cg->bbox[1] << ":" << cg->bbox[4] << "," << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
            ERROR = 1;
          }

          // euler step and boundary
          {
            MyList<var> *varl0 = StateList, *varl1 = SynchList_cor, *varlrhs = RHSList; // we do not check the correspondence here
            while (varl0)
            {
              f_sommerfeld_routbam_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                      sPp->data->bbox[0], sPp->data->bbox[1], sPp->data->bbox[2], sPp->data->bbox[3], sPp->data->bbox[4], sPp->data->bbox[5],
                                      cg->fgfs[varl1->data->sgfn],
                                      cg->fgfs[varl0->data->sgfn], varl0->data->propspeed, varl0->data->SoA,
                                      Symmetry);

              f_euler_rout(cg->shape, dT_lev, cg->fgfs[varl0->data->sgfn], cg->fgfs[varl1->data->sgfn], cg->fgfs[varlrhs->data->sgfn]);

              varl0 = varl0->next;
              varl1 = varl1->next;
              varlrhs = varlrhs->next;
            }
          }
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
      cout << " Shell stuff synchronization used " << (double)(curr_clock - prev_clock) / ((double)CLOCKS_PER_SEC) << " seconds! " << endl;
    }
  }
#endif
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
      cout << " CS_Inter used " << (double)(curr_clock - prev_clock) / ((double)CLOCKS_PER_SEC) << " seconds! " << endl;
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
}
#endif
void scalar_class::RestrictProlong(int lev, int YN, bool BB)
{
  if (lev > 0)
  {
    MyList<Patch> *Pp, *Ppc;
    if (lev > trfls && YN == 0) // time refinement levels and for intermediat time level
    {
      Pp = GH->PatL[lev - 1];
      while (Pp)
      {
        if (BB)
          Parallel::prepare_inter_time_level(Pp->data, StateList, OldStateList, SynchList_cor,
                                             SynchList_pre, 0); // use SynchList_pre as temporal storage space
        else
          Parallel::prepare_inter_time_level(Pp->data, StateList, OldStateList,
                                             SynchList_pre, 0); // use SynchList_pre as temporal storage space
        Pp = Pp->next;
      }

      Parallel::Restrict(GH->PatL[lev - 1], GH->PatL[lev], SynchList_cor, SynchList_pre, Symmetry);

      Parallel::Sync(GH->PatL[lev - 1], SynchList_pre, Symmetry);

      Ppc = GH->PatL[lev - 1];
      while (Ppc)
      {
        Pp = GH->PatL[lev];
        while (Pp)
        {
          Parallel::OutBdLow2Hi(Ppc->data, Pp->data, SynchList_pre, SynchList_cor, Symmetry);
          Pp = Pp->next;
        }
        Ppc = Ppc->next;
      }
    }
    else // no time refinement levels and for all same time levels
    {
      Parallel::Restrict(GH->PatL[lev - 1], GH->PatL[lev], SynchList_cor, StateList, Symmetry);

      Parallel::Sync(GH->PatL[lev - 1], StateList, Symmetry);

      Ppc = GH->PatL[lev - 1];
      while (Ppc)
      {
        Pp = GH->PatL[lev];
        while (Pp)
        {
          Parallel::OutBdLow2Hi(Ppc->data, Pp->data, StateList, SynchList_cor, Symmetry);
          Pp = Pp->next;
        }
        Ppc = Ppc->next;
      }
    }

    Parallel::Sync(GH->PatL[lev], SynchList_cor, Symmetry);
  }
}
void scalar_class::ProlongRestrict(int lev, int YN, bool BB)
{
  if (lev > 0)
  {
    MyList<Patch> *Pp, *Ppc;
    if (lev > trfls && YN == 0) // time refinement levels and for intermediat time level
    {
      Pp = GH->PatL[lev - 1];
      while (Pp)
      {
        if (BB)
          Parallel::prepare_inter_time_level(Pp->data, StateList, OldStateList, SynchList_cor,
                                             SynchList_pre, 0); // use SynchList_pre as temporal storage space
        else
          Parallel::prepare_inter_time_level(Pp->data, StateList, OldStateList,
                                             SynchList_pre, 0); // use SynchList_pre as temporal storage space
        Pp = Pp->next;
      }

      Ppc = GH->PatL[lev - 1];
      while (Ppc)
      {
        Pp = GH->PatL[lev];
        while (Pp)
        {
          Parallel::OutBdLow2Hi(Ppc->data, Pp->data, SynchList_pre, SynchList_cor, Symmetry);
          Pp = Pp->next;
        }
        Ppc = Ppc->next;
      }
    }
    else // no time refinement levels and for all same time levels
    {
      Ppc = GH->PatL[lev - 1];
      while (Ppc)
      {
        Pp = GH->PatL[lev];
        while (Pp)
        {
          Parallel::OutBdLow2Hi(Ppc->data, Pp->data, StateList, SynchList_cor, Symmetry);
          Pp = Pp->next;
        }
        Ppc = Ppc->next;
      }

      Parallel::Restrict(GH->PatL[lev - 1], GH->PatL[lev], SynchList_cor, StateList, Symmetry);

      Parallel::Sync(GH->PatL[lev - 1], StateList, Symmetry);
    }

    Parallel::Sync(GH->PatL[lev], SynchList_cor, Symmetry);
  }
}
