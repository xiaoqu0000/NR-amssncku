// $Id: testNull.C,v 1.8 2013/03/06 04:16:04 zjcao Exp $
#ifdef newc
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
using namespace std;
#else
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#endif

#include <mpi.h>

#include "misc.h"
#include "macrodef.h"
#include "NullShellPatch.h"
#include "monitor.h"
#include "surface_integral.h"

#define PI M_PI
//=======================================
int main(int argc, char *argv[])
{
  int myrank = 0, nprocs = 1;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  int checkrun;
  char checkfilename[50];
  int Steps;
  double StartTime, TotalTime;
  double AnasTime, DumpTime, CheckTime;
  double Courant;
  double numepss, numepsb;
  int Symmetry;
  int a_lev, maxl, decn;
  double maxrex, drex;

  int shapei[dim];
  double Rmin, xmin, xmax;

  // double RJerror[2];
  double RJerror;
  // read parameter from file
  {
    char filename[100] = "input.par";
    const int LEN = 256;
    char pline[LEN];
    string str, sgrp, skey, sval;
    int sind;
    ifstream inf(filename, ifstream::in);
    if (!inf.good())
    {
      cout << "Can not open parameter file " << filename
           << " for inputing information of Shell patches" << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    for (int i = 1; inf.good(); i++)
    {
      inf.getline(pline, LEN);
      str = pline;

      int status = misc::parse_parts(str, sgrp, skey, sval, sind);
      if (status == -1)
      {
        cout << "error reading parameter file " << filename << " in line " << i << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
      else if (status == 0)
        continue;

      if (sgrp == "BSSN")
      {
        if (skey == "Shell shape")
          shapei[sind] = atof(sval.c_str());
        else if (skey == "Rmin")
          Rmin = atof(sval.c_str());
        else if (skey == "xmin")
          xmin = atof(sval.c_str());
        else if (skey == "xmax")
          xmax = atof(sval.c_str());
      }
      if (sgrp == "ABE")
      {
        if (skey == "Symmetry")
          Symmetry = atoi(sval.c_str());
        else if (skey == "Courant")
          Courant = atof(sval.c_str());
        else if (skey == "DumpTime")
          DumpTime = atof(sval.c_str());
        else if (skey == "TotalTime")
          TotalTime = atof(sval.c_str());
        else if (skey == "AnalysisTime")
          AnasTime = atof(sval.c_str());
        else if (skey == "Max mode l")
          maxl = atoi(sval.c_str());
      }
    }
    inf.close();
  }

  monitor *ECmonitor, *NewsMonitor;
  // setup Monitors
  {
    stringstream a_stream;
    a_stream.setf(ios::left);
    a_stream << "# time L2norm_of_error";
    ECmonitor = new monitor("error.dat", myrank, a_stream.str());

    a_stream.clear();
    a_stream.str("");
    a_stream << setw(15) << "# time";
    char str[50];
    for (int pl = 2; pl < maxl + 1; pl++)
      for (int pm = -pl; pm < pl + 1; pm++)
      {
        sprintf(str, "R%02dm%03d", pl, pm);
        a_stream << setw(16) << str;
        sprintf(str, "I%02dm%03d", pl, pm);
        a_stream << setw(16) << str;
      }
    NewsMonitor = new monitor("null_news.dat", myrank, a_stream.str());
  }
  //===========================the computation body====================================================
  NullShellPatch *ADM;
  surface_integral *Waveshell;
  // setup sphere integration engine
  Waveshell = new surface_integral(Symmetry);

  ADM = new NullShellPatch(shapei, Rmin, xmin, xmax, Symmetry, myrank);
  ADM->compose_sh(nprocs);
  ADM->Setup_dyad();
  ADM->Dump_xyz(0, 0, 1);
  ADM->setupintintstuff(nprocs, 0, Symmetry);

  double PhysTime = 0, dT = Courant * PI / 4 / shapei[0];
  double LastDump = 0, LastAnas = 0;

  ADM->Setup_Initial_Data(false, PhysTime);
  while (PhysTime < TotalTime)
  {
    if (LastAnas >= AnasTime)
    {
      double *RP, *IP;
      int NN = 0;
      for (int pl = 2; pl < maxl + 1; pl++)
        for (int pm = -pl; pm < pl + 1; pm++)
          NN++;
      RP = new double[NN];
      IP = new double[NN];
//      ADM->Check_News(PhysTime,dT,false);
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
      Waveshell->surf_Wave(ADM->xmax, 0, ADM, ADM->RNews, ADM->INews, 2, maxl, NN, RP, IP, 0);
#else
#ifdef Cell
      Waveshell->surf_Wave(ADM->xmax - (ADM->getdX(2)) / 2.0, 0, ADM, ADM->RNews, ADM->INews, 2, maxl, NN, RP, IP, 0);
#else
#error Not define Vertex nor Cell
#endif
#endif
      NewsMonitor->writefile(PhysTime, NN, RP, IP);
      delete[] RP;
      delete[] IP;

      RJerror = ADM->Error_Check(PhysTime, dT, (LastDump >= DumpTime));
      // RJerror[1]=ADM->News_Error_Check(PhysTime,dT,(LastDump >= DumpTime));
      //     RJerror[0]=ADM->EqTheta_Check(PhysTime,dT,(LastDump >= DumpTime));

      ECmonitor->writefile(PhysTime, 1, &RJerror);

      LastAnas = 0;
    }

    if (LastDump >= DumpTime)
    {
      ADM->Dump_Data(ADM->StateList, 0, PhysTime, dT);
      LastDump = 0;
    }

    ADM->Step(dT, PhysTime, 0);
    PhysTime += dT;
    LastDump += dT;
    LastAnas += dT;
    if (myrank == 0)
      cout << "Time = " << PhysTime << endl;
    //     ADM->Dump_Data(ADM->StateList,0,PhysTime,dT);
  }

  ADM->Dump_Data(ADM->StateList, 0, PhysTime, dT);
  delete ADM;
  delete ECmonitor;
  delete NewsMonitor;
  delete Waveshell;
  //=======================caculation done=============================================================
  if (myrank == 0)
    cout << "===============================================================" << endl;
  if (myrank == 0)
    cout << "Simulation is successfully done!!" << endl;
  MPI_Finalize();

  exit(0);
}
