
#ifndef BSSNEM_CLASS_H
#define BSSNEM_CLASS_H

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
#include "surface_integral.h"

#include "macrodef.h"

#ifdef USE_GPU
#include "bssn_gpu_class.h"
#else
#include "bssn_class.h"
#endif

class bssnEM_class : public bssn_class
{
public:
     bssnEM_class(double Couranti, double StartTimei, double TotalTimei, double DumpTimei, double d2DumpTimei, double CheckTimei, double AnasTimei,
                  int Symmetryi, int checkruni, char *checkfilenamei, double numepssi, double numepsbi, double numepshi,
                  int a_levi, int maxli, int decni, double maxrexi, double drexi);
     ~bssnEM_class();

     void Initialize();
     void Read_Ansorg();
     void Setup_Initial_Data();
     void Step(int lev, int YN);
     void Compute_Phi2(int lev);
     void AnalysisStuff_EM(int lev, double dT_lev);
     void Interp_Constraint();

protected:
     var *Exo, *Eyo, *Ezo, *Bxo, *Byo, *Bzo, *Kpsio, *Kphio;
     var *Ex0, *Ey0, *Ez0, *Bx0, *By0, *Bz0, *Kpsi0, *Kphi0;
     var *Ex, *Ey, *Ez, *Bx, *By, *Bz, *Kpsi, *Kphi;
     var *Ex1, *Ey1, *Ez1, *Bx1, *By1, *Bz1, *Kpsi1, *Kphi1;
     var *Ex_rhs, *Ey_rhs, *Ez_rhs, *Bx_rhs, *By_rhs, *Bz_rhs, *Kpsi_rhs, *Kphi_rhs;
     var *Jx, *Jy, *Jz, *qchar;
     var *Rphi2, *Iphi2;
     var *Rphi1, *Iphi1;

     monitor *Phi2Monitor;
     monitor *Phi1Monitor;
};
#endif /* BSSNEM_CLASS_H */
