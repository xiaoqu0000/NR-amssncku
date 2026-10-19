
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
#include "microdef.h"
#include "scalar_class.h"

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
  double AnasTime, DumpTime, d2DumpTime, CheckTime;
  double Courant;
  double numepss, numepsb, numepsh;
  int Symmetry;
  int a_lev, maxl, decn;
  double maxrex, drex;
  // read parameter from file
  {
    const int LEN = 256;
    char pline[LEN];
    string str, sgrp, skey, sval;
    int sind;
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
    ifstream inf(pname, ifstream::in);
    if (!inf.good() && myrank == 0)
    {
      cout << "Can not open parameter file " << pname << endl;
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

      if (sgrp == "ABE")
      {
        if (skey == "checkrun")
          checkrun = atoi(sval.c_str());
        else if (skey == "checkfile")
          strcpy(checkfilename, sval.c_str());
        else if (skey == "Steps")
          Steps = atoi(sval.c_str());
        else if (skey == "StartTime")
          StartTime = atof(sval.c_str());
        else if (skey == "TotalTime")
          TotalTime = atof(sval.c_str());
        else if (skey == "DumpTime")
          DumpTime = atof(sval.c_str());
        else if (skey == "d2DumpTime")
          d2DumpTime = atof(sval.c_str());
        else if (skey == "CheckTime")
          CheckTime = atof(sval.c_str());
        else if (skey == "AnalysisTime")
          AnasTime = atof(sval.c_str());
        else if (skey == "Courant")
          Courant = atof(sval.c_str());
        else if (skey == "Symmetry")
          Symmetry = atoi(sval.c_str());
        else if (skey == "small dissipation")
          numepss = atof(sval.c_str());
        else if (skey == "big dissipation")
          numepsb = atof(sval.c_str());
        else if (skey == "shell dissipation")
          numepsh = atof(sval.c_str());
        else if (skey == "Analysis Level")
          a_lev = atoi(sval.c_str());
        else if (skey == "Max mode l")
          maxl = atoi(sval.c_str());
        else if (skey == "detector number")
          decn = atoi(sval.c_str());
        else if (skey == "farest detector position")
          maxrex = atof(sval.c_str());
        else if (skey == "detector distance")
          drex = atof(sval.c_str());
      }
    }
    inf.close();
  }
  // echo parameters
  if (myrank == 0)
  {
    cout << "///////////////////////////////////////////////////////////////" << endl;
#ifdef Cell
    cout << "Cell center numerical grid structure" << endl;
#endif
#ifdef Vertex
    cout << "Vertex center numerical grid structure" << endl;
#endif
    if (checkrun)
      cout << "                             checked run" << endl;
    else
      cout << "                                 new run" << endl;
    cout << "  simulation with cpu numbers = " << nprocs << endl;
    cout << "              simulation time = (" << StartTime << ", " << TotalTime << ")" << endl;
    cout << "simulation steps for this run = " << Steps << endl;
    cout << "               Courant number = " << Courant << endl;
    cout << "                   ghost zone = " << ghost_width << endl;
    cout << "                  buffer zone = " << buffer_width << endl;
    switch (Symmetry)
    {
    case 0:
      cout << "           Symmetry setting: No_Symmetry" << endl;
      break;
    case 1:
      cout << "           Symmetry setting: Equatorial" << endl;
      break;
    case 2:
      cout << "           Symmetry setting: Octant" << endl;
      break;
    default:
      cout << "OOOOps, not supported Symmetry setting!" << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    cout << "Courant = " << Courant << endl;
    cout << "artificial dissipation for shell patches = " << numepsh << endl;
    cout << "artificial dissipation for fixed levels = " << numepsb << endl;
    cout << "artificial dissipation for moving levels = " << numepss << endl;
    cout << "Dumpt Time = " << DumpTime << endl;
    cout << "Check Time = " << CheckTime << endl;
    cout << "Analysis Time = " << AnasTime << endl;
    cout << "Analysis level = " << a_lev << endl;
    cout << "checkfile = " << checkfilename << endl;
    switch (ghost_width)
    {
    case 2:
      cout << "second order finite difference is used" << endl;
      break;
    case 3:
      cout << "fourth order finite difference is used" << endl;
      break;
    case 4:
      cout << "sixth order finite difference is used" << endl;
      break;
    case 5:
      cout << "eighth order finite difference is used" << endl;
      break;
    default:
      cout << "Why are you using ghost width = " << ghost_width << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    cout << "///////////////////////////////////////////////////////////////" << endl;
  }
  //===========================the computation body====================================================
  scalar_class *ADM;

  ADM = new scalar_class(Courant, StartTime, TotalTime, DumpTime, CheckTime, AnasTime,
                         Symmetry, checkrun, checkfilename, numepss, numepsb,
                         a_lev);

  ADM->Setup_Initial_Data();

  ADM->Evolve(Steps);

  delete ADM;
  //=======================caculation done=============================================================
  if (myrank == 0)
    cout << "===============================================================" << endl;
  if (myrank == 0)
    cout << "Simulation is successfully done!!" << endl;
  MPI_Finalize();

  exit(0);
}
