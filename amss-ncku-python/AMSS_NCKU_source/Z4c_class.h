
#ifndef Z4c_CLASS_H
#define Z4c_CLASS_H

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

class Z4c_class : public bssn_class
{
public:
     Z4c_class(double Couranti, double StartTimei, double TotalTimei, double DumpTimei, double d2DumpTimei, double CheckTimei, double AnasTimei,
               int Symmetryi, int checkruni, char *checkfilenamei, double numepssi, double numepsbi, double numepshi,
               int a_levi, int maxli, int decni, double maxrexi, double drexi);
     ~Z4c_class();

     void Initialize();
     void Check_extrop();
     // Since we have set zero to variables at very begining
     // we can neglect TZ for initial data setting
     void Step(int lev, int YN);
     void Interp_Constraint();
     void Constraint_Out();
     void Compute_Constraint();

protected:
     var *TZo;
     var *TZ0;
     var *TZ;
     var *TZ1;
     var *TZ_rhs;
};
#endif /* Z4c_CLASS_H */
