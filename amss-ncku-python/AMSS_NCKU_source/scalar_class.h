
#ifndef SCALAR_CLASS_H
#define SCALAR_CLASS_H

#ifdef newc
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
using namespace std;
#else
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#endif

#include <mpi.h>

#include "cgh.h"
#include "ShellPatch.h"
#include "misc.h"
#include "var.h"
#include "MyList.h"
#include "monitor.h"

class scalar_class
{
protected:
     int myrank;
     cgh *GH;
     ShellPatch *SH;
     double PhysTime;

     int checkrun;
     char checkfilename[50];
     int Steps;
     double StartTime, TotalTime;
     double AnasTime, DumpTime, CheckTime;
     double LastAnas;
     double Courant;
     double numepss, numepsb;
     int Symmetry;
     int trfls, a_lev;

     double dT;

     var *Sphio, *Spio;
     var *Sphi0, *Spi0;
     var *Sphi, *Spi;
     var *Sphi1, *Spi1;
     var *Sphi_rhs, *Spi_rhs;

     MyList<var> *StateList, *SynchList_pre, *SynchList_cor, *RHSList;
     MyList<var> *OldStateList, *DumpList, *CheckList;

     monitor *ErrorMonitor;

public:
     scalar_class(double Couranti, double StartTimei, double TotalTimei, double DumpTimei, double CheckTimei, double AnasTimei,
                  int Symmetryi, int checkruni, char *checkfilenamei, double numepssi, double numepsbi,
                  int a_levi);
     ~scalar_class();
     void Setup_Initial_Data();
     void Evolve(int Steps);
     void RecursiveStep(int lev);
     void Step(int lev, int YN);
     void RestrictProlong(int lev, int YN, bool BB);
     void ProlongRestrict(int lev, int YN, bool BB);
};
#endif /* SCALAR_CLASS_H */
