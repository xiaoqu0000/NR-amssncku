
#ifndef BSSNESCALAR_CLASS_H
#define BSSNESCALAR_CLASS_H

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

class bssnEScalar_class : public bssn_class
{
public:
     bssnEScalar_class(double Couranti, double StartTimei, double TotalTimei, double DumpTimei, double d2DumpTimei, double CheckTimei, double AnasTimei,
                       int Symmetryi, int checkruni, char *checkfilenamei, double numepssi, double numepsbi, double numepshi,
                       int a_levi, int maxli, int decni, double maxrexi, double drexi);
     ~bssnEScalar_class();

     void Initialize();
     void Read_Ansorg();
     void Read_Pablo();
     void Compute_Psi4(int lev);
     void Step(int lev, int YN);
     void AnalysisStuff_EScalar(int lev, double dT_lev);
     void Interp_Constraint();
     void Constraint_Out(); 

protected:
     var *Sphio, *Spio;
     var *Sphi0, *Spi0;
     var *Sphi,  *Spi;
     var *Sphi1, *Spi1;
     var *Sphi_rhs, *Spi_rhs;

     var *Cons_fR;

     monitor *MaxScalar_Monitor;
};

#endif /* BSSNESCALAR_CLASS_H */

