
#ifndef BSSNEM_GPU_CLASS_H
#define BSSNEM_GPU_CLASS_H

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

#include "macrodef.h"
#include "cgh.h"
#include "ShellPatch.h"
#include "misc.h"
#include "var.h"
#include "MyList.h"
#include "monitor.h"
#include "surface_integral.h"
#include "checkpoint.h"

// GPU infrastructure
#include "gpu_manager.h"
#include "helper.h"

#include "bssn_gpu_class.h"

class bssnEM_class : public bssn_class
{
public:
    bssnEM_class(double Couranti, double StartTimei, double TotalTimei,
                 double DumpTimei, double d2DumpTimei, double CheckTimei, double AnasTimei,
                 int Symmetryi, int checkruni, char *checkfilenamei,
                 double numepssi, double numepsbi, double numepshi,
                 int a_levi, int maxli, int decni, double maxrexi, double drexi);
    ~bssnEM_class();

    void Initialize();
    void Read_Ansorg();
    void Setup_Initial_Data();
    void Step_GPU(int lev, int YN);
    void Compute_Phi2_GPU(int lev);
    void AnalysisStuff_EM(int lev, double dT_lev);
    void Interp_Constraint();

protected:
    // EM evolved fields — 5 time levels each
    var *Exo, *Eyo, *Ezo, *Bxo, *Byo, *Bzo, *Kpsio, *Kphio;
    var *Ex0, *Ey0, *Ez0, *Bx0, *By0, *Bz0, *Kpsi0, *Kphi0;
    var *Ex,  *Ey,  *Ez,  *Bx,  *By,  *Bz,  *Kpsi,  *Kphi;
    var *Ex1, *Ey1, *Ez1, *Bx1, *By1, *Bz1, *Kpsi1, *Kphi1;
    var *Ex_rhs, *Ey_rhs, *Ez_rhs, *Bx_rhs, *By_rhs, *Bz_rhs, *Kpsi_rhs, *Kphi_rhs;
    // EM source / diagnostic
    var *Jx, *Jy, *Jz, *qchar;
    var *Rphi2, *Iphi2;
    var *Rphi1, *Iphi1;

    monitor *Phi2Monitor;
    monitor *Phi1Monitor;
};

#endif /* BSSNEM_GPU_CLASS_H */
