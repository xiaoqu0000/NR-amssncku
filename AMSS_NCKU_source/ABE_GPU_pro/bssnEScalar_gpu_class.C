
#ifdef newc
#include <sstream>
#include <cstdio>
using namespace std;
#else
#include <stdio.h>
#endif

#include "macrodef.h"
#include "misc.h"
#include <time.h>
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

#ifdef With_AHF
#include "derivatives.h"
#include "myglobal.h"
#endif

#include "perf.h"
#include "derivatives.h"
#include "ricci_gamma.h"

// include GPU files
#include "gpu_manager.h"
#include "helper.h"
#include "macrodef.h"

#include "scalar_rhs.h"
#include "bssnEScalar_rhs_gpu.cuh"

// NOTE: This file provides the GPU implementation of bssnEScalar_class.
// It defines the constructor, destructor, Initialize(), and Step_GPU()
// for the BSSN-EScalar equation class when USE_GPU is enabled.
// The structure mirrors bssn_gpu_class.C (for pure BSSN) and
// bssn_step_gpu.C (for Step_GPU), adapted for the scalar field.

#if (PSTR != 0)
#warning "PSTR != 0 alternate Step_GPU implementations were removed during fused-kernel optimization; only PSTR=0 is provided."
#endif

//================================================================================================

// 构造函数，定义 bssnEScalar_gpu_class
// Constructor: mirrors bssn_gpu_class.C with additional EScalar variables

//================================================================================================

bssnEScalar_class::bssnEScalar_class(double Couranti, double StartTimei, double TotalTimei,
                                     double DumpTimei, double d2DumpTimei, double CheckTimei, double AnasTimei,
                                     int Symmetryi, int checkruni, char *checkfilenamei,
                                     double numepssi, double numepsbi, double numepshi,
                                     int a_levi, int maxli, int decni, double maxrexi, double drexi)
                                     : bssn_class(Couranti, StartTimei, TotalTimei,
                                                  DumpTimei, d2DumpTimei, CheckTimei, AnasTimei,
                                                  Symmetryi, checkruni, checkfilenamei, numepssi, numepsbi, numepshi,
                                                  a_levi, maxli, decni, maxrexi, drexi)
{
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  // setup the EScalar-specific monitor only.
  // ErrorMonitor / Psi4Monitor / BHMonitor / MAPMonitor / ConVMonitor are already
  // created by the base bssn_class GPU constructor. Re-creating them here makes
  // monitor() version-suffix the filenames (bssn_BH.dat -> 1_bssn_BH.dat, ...), so
  // the standard diagnostics never reach the expected bssn_*.dat outputs.
  {
    stringstream a_stream;
    a_stream.setf(ios::left);
    a_stream << setw(15) << "# time x y z maxs";
    MaxScalar_Monitor = new monitor("bssn_maxs.dat", myrank, a_stream.str());
  }
  // setup sphere integration engine
  Waveshell = new surface_integral(Symmetry);

  trfls = 0;
  chitiny = 0;
  // read parameter from file
  {
    char filename[50];
    strcpy(filename, "input.par");
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

      if (sgrp == "BSSN" && skey == "chitiny")
        chitiny = atof(sval.c_str());
      else if (sgrp == "BSSN" && skey == "time refinement start from level")
        trfls = atoi(sval.c_str());
#ifdef With_AHF
      else if (sgrp == "AHF" && skey == "AHfindevery")
        AHfindevery = atoi(sval.c_str());
      else if (sgrp == "AHF" && skey == "AHdumptime")
        AHdumptime = atof(sval.c_str());
#endif
    }
    inf.close();
  }
  if (myrank == 0)
  {
    cout << "chitiny = " << chitiny << endl;
    cout << "time refinement start from level #" << trfls << endl;
#ifdef With_AHF
    cout << "parameters for AHF:" << endl;
    cout << "AHfindevery = " << AHfindevery << endl;
    cout << "AHdumptime = " << AHdumptime << endl;
#endif
  }

  chitiny = chitiny - 1; // because we have subtracted one from chi

  strcpy(checkfilename, checkfilenamei);

  ngfs = 0;
  // --- BSSN OldState variables (same as bssn_gpu_class.C) ---
  phio = new var("phio", ngfs++, 1, 1, 1);  trKo = new var("trKo", ngfs++, 1, 1, 1);
  gxxo = new var("gxxo", ngfs++, 1, 1, 1);  gxyo = new var("gxyo", ngfs++, -1, -1, 1); gxzo = new var("gxzo", ngfs++, -1, 1, -1);
  gyyo = new var("gyyo", ngfs++, 1, 1, 1);  gyzo = new var("gyzo", ngfs++, 1, -1, -1); gzzo = new var("gzzo", ngfs++, 1, 1, 1);
  Axxo = new var("Axxo", ngfs++, 1, 1, 1);  Axyo = new var("Axyo", ngfs++, -1, -1, 1); Axzo = new var("Axzo", ngfs++, -1, 1, -1);
  Ayyo = new var("Ayyo", ngfs++, 1, 1, 1);  Ayzo = new var("Ayzo", ngfs++, 1, -1, -1); Azzo = new var("Azzo", ngfs++, 1, 1, 1);
  Gmxo = new var("Gmxo", ngfs++, -1, 1, 1); Gmyo = new var("Gmyo", ngfs++, 1, -1, 1);  Gmzo = new var("Gmzo", ngfs++, 1, 1, -1);
  Lapo = new var("Lapo", ngfs++, 1, 1, 1);
  Sfxo   = new var("Sfxo", ngfs++, -1, 1, 1);   Sfyo   = new var("Sfyo", ngfs++, 1, -1, 1);   Sfzo = new var("Sfzo", ngfs++, 1, 1, -1);
  dtSfxo = new var("dtSfxo", ngfs++, -1, 1, 1); dtSfyo = new var("dtSfyo", ngfs++, 1, -1, 1); dtSfzo = new var("dtSfzo", ngfs++, 1, 1, -1);
  // --- EScalar OldState variables ---
  Sphio = new var("Sphio", ngfs++, 1, 1, 1); 
  Spio = new var("Spio", ngfs++, 1, 1, 1);

  // --- BSSN State variables ---
  phi0 = new var("phi0", ngfs++, 1, 1, 1);  trK0 = new var("trK0", ngfs++, 1, 1, 1);
  gxx0 = new var("gxx0", ngfs++, 1, 1, 1);  gxy0 = new var("gxy0", ngfs++, -1, -1, 1); gxz0 = new var("gxz0", ngfs++, -1, 1, -1);
  gyy0 = new var("gyy0", ngfs++, 1, 1, 1);  gyz0 = new var("gyz0", ngfs++, 1, -1, -1); gzz0 = new var("gzz0", ngfs++, 1, 1, 1);
  Axx0 = new var("Axx0", ngfs++, 1, 1, 1);  Axy0 = new var("Axy0", ngfs++, -1, -1, 1); Axz0 = new var("Axz0", ngfs++, -1, 1, -1);
  Ayy0 = new var("Ayy0", ngfs++, 1, 1, 1);  Ayz0 = new var("Ayz0", ngfs++, 1, -1, -1); Azz0 = new var("Azz0", ngfs++, 1, 1, 1);
  Gmx0 = new var("Gmx0", ngfs++, -1, 1, 1); Gmy0 = new var("Gmy0", ngfs++, 1, -1, 1);  Gmz0 = new var("Gmz0", ngfs++, 1, 1, -1);
  Lap0 = new var("Lap0", ngfs++, 1, 1, 1);
  Sfx0   = new var("Sfx0", ngfs++, -1, 1, 1);   Sfy0   = new var("Sfy0", ngfs++, 1, -1, 1);   Sfz0   = new var("Sfz0", ngfs++, 1, 1, -1);
  dtSfx0 = new var("dtSfx0", ngfs++, -1, 1, 1); dtSfy0 = new var("dtSfy0", ngfs++, 1, -1, 1); dtSfz0 = new var("dtSfz0", ngfs++, 1, 1, -1);
  // --- EScalar State variables ---
  Sphi0 = new var("Sphi0", ngfs++, 1, 1, 1);
  Spi0 = new var("Spi0", ngfs++, 1, 1, 1);

  // --- BSSN SynchList_pre variables ---
  phi = new var("phi", ngfs++, 1, 1, 1);  trK = new var("trK", ngfs++, 1, 1, 1);
  gxx = new var("gxx", ngfs++, 1, 1, 1);  gxy = new var("gxy", ngfs++, -1, -1, 1); gxz = new var("gxz", ngfs++, -1, 1, -1);
  gyy = new var("gyy", ngfs++, 1, 1, 1);  gyz = new var("gyz", ngfs++, 1, -1, -1); gzz = new var("gzz", ngfs++, 1, 1, 1);
  Axx = new var("Axx", ngfs++, 1, 1, 1);  Axy = new var("Axy", ngfs++, -1, -1, 1); Axz = new var("Axz", ngfs++, -1, 1, -1);
  Ayy = new var("Ayy", ngfs++, 1, 1, 1);  Ayz = new var("Ayz", ngfs++, 1, -1, -1); Azz = new var("Azz", ngfs++, 1, 1, 1);
  Gmx = new var("Gmx", ngfs++, -1, 1, 1); Gmy = new var("Gmy", ngfs++, 1, -1, 1);  Gmz = new var("Gmz", ngfs++, 1, 1, -1);
  Lap = new var("Lap", ngfs++, 1, 1, 1);
  Sfx   = new var("Sfx", ngfs++, -1, 1, 1);   Sfy   = new var("Sfy", ngfs++, 1, -1, 1);   Sfz   = new var("Sfz", ngfs++, 1, 1, -1);
  dtSfx = new var("dtSfx", ngfs++, -1, 1, 1); dtSfy = new var("dtSfy", ngfs++, 1, -1, 1); dtSfz = new var("dtSfz", ngfs++, 1, 1, -1);
  // --- EScalar SynchList_pre variables ---
  Sphi = new var("Sphi", ngfs++, 1, 1, 1);
  Spi  = new var("Spi", ngfs++, 1, 1, 1);

  // --- BSSN SynchList_cor variables ---
  phi1 = new var("phi1", ngfs++, 1, 1, 1);  trK1 = new var("trK1", ngfs++, 1, 1, 1);
  gxx1 = new var("gxx1", ngfs++, 1, 1, 1);  gxy1 = new var("gxy1", ngfs++, -1, -1, 1); gxz1 = new var("gxz1", ngfs++, -1, 1, -1);
  gyy1 = new var("gyy1", ngfs++, 1, 1, 1);  gyz1 = new var("gyz1", ngfs++, 1, -1, -1); gzz1 = new var("gzz1", ngfs++, 1, 1, 1);
  Axx1 = new var("Axx1", ngfs++, 1, 1, 1);  Axy1 = new var("Axy1", ngfs++, -1, -1, 1); Axz1 = new var("Axz1", ngfs++, -1, 1, -1);
  Ayy1 = new var("Ayy1", ngfs++, 1, 1, 1);  Ayz1 = new var("Ayz1", ngfs++, 1, -1, -1); Azz1 = new var("Azz1", ngfs++, 1, 1, 1);
  Gmx1 = new var("Gmx1", ngfs++, -1, 1, 1); Gmy1 = new var("Gmy1", ngfs++, 1, -1, 1);  Gmz1 = new var("Gmz1", ngfs++, 1, 1, -1);
  Lap1 = new var("Lap1", ngfs++, 1, 1, 1);
  Sfx1   = new var("Sfx1", ngfs++, -1, 1, 1);   Sfy1   = new var("Sfy1", ngfs++, 1, -1, 1);   Sfz1   = new var("Sfz1", ngfs++, 1, 1, -1);
  dtSfx1 = new var("dtSfx1", ngfs++, -1, 1, 1); dtSfy1 = new var("dtSfy1", ngfs++, 1, -1, 1); dtSfz1 = new var("dtSfz1", ngfs++, 1, 1, -1);
  // --- EScalar SynchList_cor variables ---
  Sphi1 = new var("Sphi1", ngfs++, 1, 1, 1);
  Spi1  = new var("Spi1", ngfs++, 1, 1, 1);

  // --- BSSN RHS variables ---
  phi_rhs = new var("phi_rhs", ngfs++, 1, 1, 1); trK_rhs = new var("trK_rhs", ngfs++, 1, 1, 1);
  gxx_rhs = new var("gxx_rhs", ngfs++, 1, 1, 1); gxy_rhs = new var("gxy_rhs", ngfs++, -1, -1, 1); gxz_rhs = new var("gxz_rhs", ngfs++, -1, 1, -1);
  gyy_rhs = new var("gyy_rhs", ngfs++, 1, 1, 1); gyz_rhs = new var("gyz_rhs", ngfs++, 1, -1, -1); gzz_rhs = new var("gzz_rhs", ngfs++, 1, 1, 1);
  Axx_rhs = new var("Axx_rhs", ngfs++, 1, 1, 1);  Axy_rhs = new var("Axy_rhs", ngfs++, -1, -1, 1); Axz_rhs = new var("Axz_rhs", ngfs++, -1, 1, -1);
  Ayy_rhs = new var("Ayy_rhs", ngfs++, 1, 1, 1);  Ayz_rhs = new var("Ayz_rhs", ngfs++, 1, -1, -1); Azz_rhs = new var("Azz_rhs", ngfs++, 1, 1, 1);
  Gmx_rhs = new var("Gmx_rhs", ngfs++, -1, 1, 1); Gmy_rhs = new var("Gmy_rhs", ngfs++, 1, -1, 1);  Gmz_rhs = new var("Gmz_rhs", ngfs++, 1, 1, -1);
  Lap_rhs = new var("Lap_rhs", ngfs++, 1, 1, 1);
  Sfx_rhs   = new var("Sfx_rhs", ngfs++, -1, 1, 1);   Sfy_rhs   = new var("Sfy_rhs", ngfs++, 1, -1, 1);   Sfz_rhs   = new var("Sfz_rhs", ngfs++, 1, 1, -1);
  dtSfx_rhs = new var("dtSfx_rhs", ngfs++, -1, 1, 1); dtSfy_rhs = new var("dtSfy_rhs", ngfs++, 1, -1, 1); dtSfz_rhs = new var("dtSfz_rhs", ngfs++, 1, 1, -1);
  // --- EScalar RHS variables ---
  Sphi_rhs = new var("Sphi_rhs", ngfs++, 1, 1, 1);
  Spi_rhs  = new var("Spi_rhs", ngfs++, 1, 1, 1);

  // --- Matter source variables ---
  rho = new var("rho", ngfs++, 1, 1, 1);
  Sx  = new var("Sx", ngfs++, -1, 1, 1); Sy  = new var("Sy", ngfs++, 1, -1, 1);   Sz  = new var("Sz", ngfs++, 1, 1, -1);
  Sxx = new var("Sxx", ngfs++, 1, 1, 1); Sxy = new var("Sxy", ngfs++, -1, -1, 1); Sxz = new var("Sxz", ngfs++, -1, 1, -1);
  Syy = new var("Syy", ngfs++, 1, 1, 1); Syz = new var("Syz", ngfs++, 1, -1, -1); Szz = new var("Szz", ngfs++, 1, 1, 1);

  Gamxxx = new var("Gamxxx", ngfs++, -1, 1, 1); Gamxxy = new var("Gamxxy", ngfs++, 1, -1, 1);   Gamxxz = new var("Gamxxz", ngfs++, 1, 1, -1);
  Gamxyy = new var("Gamxyy", ngfs++, -1, 1, 1); Gamxyz = new var("Gamxyz", ngfs++, -1, -1, -1);
  Gamxzz = new var("Gamxzz", ngfs++, -1, 1, 1);
  Gamyxx = new var("Gamyxx", ngfs++, 1, -1, 1); Gamyxy = new var("Gamyxy", ngfs++, -1, 1, 1);   Gamyxz = new var("Gamyxz", ngfs++, -1, -1, -1);
  Gamyyy = new var("Gamyyy", ngfs++, 1, -1, 1); Gamyyz = new var("Gamyyz", ngfs++, 1, 1, -1);
  Gamyzz = new var("Gamyzz", ngfs++, 1, -1, 1);
  Gamzxx = new var("Gamzxx", ngfs++, 1, 1, -1); Gamzxy = new var("Gamzxy", ngfs++, -1, -1, -1); Gamzxz = new var("Gamzxz", ngfs++, -1, 1, 1);
  Gamzyy = new var("Gamzyy", ngfs++, 1, 1, -1); Gamzyz = new var("Gamzyz", ngfs++, 1, -1, 1); 
  Gamzzz = new var("Gamzzz", ngfs++, 1, 1, -1);

  Rxx = new var("Rxx", ngfs++, 1, 1, 1); Rxy = new var("Rxy", ngfs++, -1, -1, 1); Rxz = new var("Rxz", ngfs++, -1, 1, -1); 
  Ryy = new var("Ryy", ngfs++, 1, 1, 1); Ryz = new var("Ryz", ngfs++, 1, -1, -1); Rzz = new var("Rzz", ngfs++, 1, 1, 1);

  Rpsi4   = new var("Rpsi4", ngfs++, 1, 1, 1);
  Ipsi4   = new var("Ipsi4", ngfs++, -1, -1, -1);
  t1Rpsi4 = new var("t1Rpsi4", ngfs++, 1, 1, 1);  t1Ipsi4 = new var("t1Ipsi4", ngfs++, -1, -1, -1);
  t2Rpsi4 = new var("t2Rpsi4", ngfs++, 1, 1, 1);  t2Ipsi4 = new var("t2Ipsi4", ngfs++, -1, -1, -1);

  Cons_Ham = new var("Cons_Ham", ngfs++, 1, 1, 1);
  Cons_Px  = new var("Cons_Px", ngfs++, -1, 1, 1); Cons_Py = new var("Cons_Py", ngfs++, 1, -1, 1); Cons_Pz = new var("Cons_Pz", ngfs++, 1, 1, -1);
  Cons_Gx  = new var("Cons_Gx", ngfs++, -1, 1, 1); Cons_Gy = new var("Cons_Gy", ngfs++, 1, -1, 1); Cons_Gz = new var("Cons_Gz", ngfs++, 1, 1, -1);
  // --- EScalar constraint variable ---
  Cons_fR = new var("Cons_fR", ngfs++, 1, 1, 1);

#ifdef Point_Psi4
  phix = new var("phix", ngfs++, -1, 1, 1);   phiy = new var("phiy", ngfs++, 1, -1, 1);   phiz = new var("phiz", ngfs++, 1, 1, -1);
  trKx = new var("trKx", ngfs++, -1, 1, 1);   trKy = new var("trKy", ngfs++, 1, -1, 1);   trKz = new var("trKz", ngfs++, 1, 1, -1);
  Axxx = new var("Axxx", ngfs++, -1, 1, 1);   Axxy = new var("Axxy", ngfs++, 1, -1, 1);   Axxz = new var("Axxz", ngfs++, 1, 1, -1);
  Axyx = new var("Axyx", ngfs++, 1, -1, 1);   Axyy = new var("Axyy", ngfs++, -1, 1, 1);   Axyz = new var("Axyz", ngfs++, -1, -1, -1);
  Axzx = new var("Axzx", ngfs++, 1, 1, -1);   Axzy = new var("Axzy", ngfs++, -1, -1, -1); Axzz = new var("Axzz", ngfs++, -1, 1, 1);
  Ayyx = new var("Ayyx", ngfs++, -1, 1, 1);   Ayyy = new var("Ayyy", ngfs++, 1, -1, 1);   Ayyz = new var("Ayyz", ngfs++, 1, 1, -1);
  Ayzx = new var("Ayzx", ngfs++, -1, -1, -1); Ayzy = new var("Ayzy", ngfs++, 1, 1, -1);   Ayzz = new var("Ayzz", ngfs++, 1, -1, 1);
  Azzx = new var("Azzx", ngfs++, -1, 1, 1);   Azzy = new var("Azzy", ngfs++, 1, -1, 1);   Azzz = new var("Azzz", ngfs++, 1, 1, -1);
#endif

  // specific properspeed for 1+log slice
  {
    const double vl = sqrt(2);
    trKo->setpropspeed(vl); trK0->setpropspeed(vl); trK->setpropspeed(vl); trK1->setpropspeed(vl); trK_rhs->setpropspeed(vl);
    phio->setpropspeed(vl); phi0->setpropspeed(vl); phi->setpropspeed(vl); phi1->setpropspeed(vl); phi_rhs->setpropspeed(vl);
    Lapo->setpropspeed(vl); Lap0->setpropspeed(vl); Lap->setpropspeed(vl); Lap1->setpropspeed(vl); Lap_rhs->setpropspeed(vl);
  }

  // --- Build lists (BSSN + EScalar variables) ---
  OldStateList = new MyList<var>(phio);
  OldStateList->insert(trKo);
  OldStateList->insert(gxxo);   OldStateList->insert(gxyo);   OldStateList->insert(gxzo);
  OldStateList->insert(gyyo);   OldStateList->insert(gyzo);   OldStateList->insert(gzzo);
  OldStateList->insert(Axxo);   OldStateList->insert(Axyo);   OldStateList->insert(Axzo);
  OldStateList->insert(Ayyo);   OldStateList->insert(Ayzo);   OldStateList->insert(Azzo);
  OldStateList->insert(Gmxo);   OldStateList->insert(Gmyo);   OldStateList->insert(Gmzo);
  OldStateList->insert(Lapo);
  OldStateList->insert(Sfxo);   OldStateList->insert(Sfyo);   OldStateList->insert(Sfzo);
  OldStateList->insert(dtSfxo); OldStateList->insert(dtSfyo); OldStateList->insert(dtSfzo);
  OldStateList->insert(Sphio);  OldStateList->insert(Spio);

  StateList = new MyList<var>(phi0);
  StateList->insert(trK0);
  StateList->insert(gxx0);   StateList->insert(gxy0);   StateList->insert(gxz0);
  StateList->insert(gyy0);   StateList->insert(gyz0);   StateList->insert(gzz0);
  StateList->insert(Axx0);   StateList->insert(Axy0);   StateList->insert(Axz0);
  StateList->insert(Ayy0);   StateList->insert(Ayz0);   StateList->insert(Azz0);
  StateList->insert(Gmx0);   StateList->insert(Gmy0);   StateList->insert(Gmz0);
  StateList->insert(Lap0);
  StateList->insert(Sfx0);   StateList->insert(Sfy0);   StateList->insert(Sfz0);
  StateList->insert(dtSfx0); StateList->insert(dtSfy0); StateList->insert(dtSfz0);
  StateList->insert(Sphi0);  StateList->insert(Spi0);

  RHSList = new MyList<var>(phi_rhs);
  RHSList->insert(trK_rhs);
  RHSList->insert(gxx_rhs);   RHSList->insert(gxy_rhs);   RHSList->insert(gxz_rhs);
  RHSList->insert(gyy_rhs);   RHSList->insert(gyz_rhs);   RHSList->insert(gzz_rhs);
  RHSList->insert(Axx_rhs);   RHSList->insert(Axy_rhs);   RHSList->insert(Axz_rhs);
  RHSList->insert(Ayy_rhs);   RHSList->insert(Ayz_rhs);   RHSList->insert(Azz_rhs);
  RHSList->insert(Gmx_rhs);   RHSList->insert(Gmy_rhs);   RHSList->insert(Gmz_rhs);
  RHSList->insert(Lap_rhs);
  RHSList->insert(Sfx_rhs);   RHSList->insert(Sfy_rhs);   RHSList->insert(Sfz_rhs);
  RHSList->insert(dtSfx_rhs); RHSList->insert(dtSfy_rhs); RHSList->insert(dtSfz_rhs);
  RHSList->insert(Sphi_rhs);  RHSList->insert(Spi_rhs);

  SynchList_pre = new MyList<var>(phi);
  SynchList_pre->insert(trK);
  SynchList_pre->insert(gxx);   SynchList_pre->insert(gxy); SynchList_pre->insert(gxz);
  SynchList_pre->insert(gyy);   SynchList_pre->insert(gyz); SynchList_pre->insert(gzz);
  SynchList_pre->insert(Axx);   SynchList_pre->insert(Axy); SynchList_pre->insert(Axz);
  SynchList_pre->insert(Ayy);   SynchList_pre->insert(Ayz); SynchList_pre->insert(Azz);
  SynchList_pre->insert(Gmx);   SynchList_pre->insert(Gmy); SynchList_pre->insert(Gmz);
  SynchList_pre->insert(Lap);
  SynchList_pre->insert(Sfx);   SynchList_pre->insert(Sfy);   SynchList_pre->insert(Sfz);
  SynchList_pre->insert(dtSfx); SynchList_pre->insert(dtSfy); SynchList_pre->insert(dtSfz);
  SynchList_pre->insert(Sphi);  SynchList_pre->insert(Spi);

  SynchList_cor = new MyList<var>(phi1);
  SynchList_cor->insert(trK1);
  SynchList_cor->insert(gxx1);   SynchList_cor->insert(gxy1);   SynchList_cor->insert(gxz1);
  SynchList_cor->insert(gyy1);   SynchList_cor->insert(gyz1);   SynchList_cor->insert(gzz1);
  SynchList_cor->insert(Axx1);   SynchList_cor->insert(Axy1);   SynchList_cor->insert(Axz1);
  SynchList_cor->insert(Ayy1);   SynchList_cor->insert(Ayz1);   SynchList_cor->insert(Azz1);
  SynchList_cor->insert(Gmx1);   SynchList_cor->insert(Gmy1);   SynchList_cor->insert(Gmz1);
  SynchList_cor->insert(Lap1);
  SynchList_cor->insert(Sfx1);   SynchList_cor->insert(Sfy1);   SynchList_cor->insert(Sfz1);
  SynchList_cor->insert(dtSfx1); SynchList_cor->insert(dtSfy1); SynchList_cor->insert(dtSfz1);
  SynchList_cor->insert(Sphi1);  SynchList_cor->insert(Spi1);

  DumpList = new MyList<var>(phi0);
  DumpList->insert(trK0);
  DumpList->insert(gxx0);     DumpList->insert(gxy0);   DumpList->insert(gxz0);
  DumpList->insert(gyy0);     DumpList->insert(gyz0);   DumpList->insert(gzz0);
  DumpList->insert(Lap0);
  DumpList->insert(Rpsi4);    DumpList->insert(Ipsi4);
  DumpList->insert(Cons_Ham);
  DumpList->insert(Cons_Px);  DumpList->insert(Cons_Py); DumpList->insert(Cons_Pz);
  DumpList->insert(Cons_Gx);  DumpList->insert(Cons_Gy); DumpList->insert(Cons_Gz);
  DumpList->insert(Sphi0);    DumpList->insert(Spi0);
  DumpList->insert(Cons_fR);

  ConstraintList = new MyList<var>(Cons_Ham);
  ConstraintList->insert(Cons_Px); ConstraintList->insert(Cons_Py); ConstraintList->insert(Cons_Pz);
  ConstraintList->insert(Cons_Gx); ConstraintList->insert(Cons_Gy); ConstraintList->insert(Cons_Gz);
  ConstraintList->insert(Cons_fR);

  MiscList = new MyList<var>(rho);
  MiscList->insert(Sx);     MiscList->insert(Sy);     MiscList->insert(Sz);
  MiscList->insert(Sxx);    MiscList->insert(Sxy);    MiscList->insert(Sxz);
  MiscList->insert(Syy);    MiscList->insert(Syz);    MiscList->insert(Szz);
  MiscList->insert(Gamxxx); MiscList->insert(Gamxxy); MiscList->insert(Gamxxz);
  MiscList->insert(Gamxyy); MiscList->insert(Gamxyz); MiscList->insert(Gamxzz);
  MiscList->insert(Gamyxx); MiscList->insert(Gamyxy); MiscList->insert(Gamyxz);
  MiscList->insert(Gamyyy); MiscList->insert(Gamyyz); MiscList->insert(Gamyzz);
  MiscList->insert(Gamzxx); MiscList->insert(Gamzxy); MiscList->insert(Gamzxz);
  MiscList->insert(Gamzyy); MiscList->insert(Gamzyz); MiscList->insert(Gamzzz);
  MiscList->insert(Rxx);    MiscList->insert(Rxy);    MiscList->insert(Rxz);
  MiscList->insert(Ryy);    MiscList->insert(Ryz);    MiscList->insert(Rzz);

  DGList = new MyList<var>(Rpsi4);
  DGList->insert(Ipsi4);
}

//================================================================================================



//================================================================================================

// 析构函数，用于删除变量
// Destructor

//================================================================================================

bssnEScalar_class::~bssnEScalar_class()
{
#ifdef With_AHF
  AHList->clearList();
  AHDList->clearList();
  GaugeList->clearList();
  if (lastahdumpid) delete[] lastahdumpid;
  if (findeveryl) delete[] findeveryl;
  if (xc) { delete[] xc; delete[] yc; delete[] zc; delete[] xr; delete[] yr; delete[] zr; delete[] trigger; delete[] dumpid; delete[] dTT; }
  AHFinderDirect::AHFinderDirect_cleanup();
#endif

  StateList->clearList();
  RHSList->clearList();
  OldStateList->clearList();
  SynchList_pre->clearList();
  SynchList_cor->clearList();
  DumpList->clearList();
  ConstraintList->clearList();
  MiscList->clearList();
  DGList->clearList();

  // BSSN OldState
  delete phio; delete trKo;
  delete gxxo; delete gxyo; delete gxzo; delete gyyo; delete gyzo; delete gzzo;
  delete Axxo; delete Axyo; delete Axzo; delete Ayyo; delete Ayzo; delete Azzo;
  delete Gmxo; delete Gmyo; delete Gmzo;
  delete Lapo;
  delete Sfxo; delete Sfyo; delete Sfzo;
  delete dtSfxo; delete dtSfyo; delete dtSfzo;
  // EScalar OldState
  delete Sphio; delete Spio;

  // BSSN State
  delete phi0; delete trK0;
  delete gxx0; delete gxy0; delete gxz0; delete gyy0; delete gyz0; delete gzz0;
  delete Axx0; delete Axy0; delete Axz0; delete Ayy0; delete Ayz0; delete Azz0;
  delete Gmx0; delete Gmy0; delete Gmz0;
  delete Lap0;
  delete Sfx0; delete Sfy0; delete Sfz0;
  delete dtSfx0; delete dtSfy0; delete dtSfz0;
  // EScalar State
  delete Sphi0; delete Spi0;

  // BSSN SynchList_pre
  delete phi; delete trK;
  delete gxx; delete gxy; delete gxz; delete gyy; delete gyz; delete gzz;
  delete Axx; delete Axy; delete Axz; delete Ayy; delete Ayz; delete Azz;
  delete Gmx; delete Gmy; delete Gmz;
  delete Lap;
  delete Sfx; delete Sfy; delete Sfz;
  delete dtSfx; delete dtSfy; delete dtSfz;
  // EScalar SynchList_pre
  delete Sphi; delete Spi;

  // BSSN SynchList_cor
  delete phi1; delete trK1;
  delete gxx1; delete gxy1; delete gxz1; delete gyy1; delete gyz1; delete gzz1;
  delete Axx1; delete Axy1; delete Axz1; delete Ayy1; delete Ayz1; delete Azz1;
  delete Gmx1; delete Gmy1; delete Gmz1;
  delete Lap1;
  delete Sfx1; delete Sfy1; delete Sfz1;
  delete dtSfx1; delete dtSfy1; delete dtSfz1;
  // EScalar SynchList_cor
  delete Sphi1; delete Spi1;

  // BSSN RHS
  delete phi_rhs; delete trK_rhs;
  delete gxx_rhs; delete gxy_rhs; delete gxz_rhs; delete gyy_rhs; delete gyz_rhs; delete gzz_rhs;
  delete Axx_rhs; delete Axy_rhs; delete Axz_rhs; delete Ayy_rhs; delete Ayz_rhs; delete Azz_rhs;
  delete Gmx_rhs; delete Gmy_rhs; delete Gmz_rhs;
  delete Lap_rhs;
  delete Sfx_rhs; delete Sfy_rhs; delete Sfz_rhs;
  delete dtSfx_rhs; delete dtSfy_rhs; delete dtSfz_rhs;
  // EScalar RHS
  delete Sphi_rhs; delete Spi_rhs;

  // Matter
  delete rho; delete Sx; delete Sy; delete Sz;
  delete Sxx; delete Sxy; delete Sxz; delete Syy; delete Syz; delete Szz;

  // Connection
  delete Gamxxx; delete Gamxxy; delete Gamxxz; delete Gamxyy; delete Gamxyz; delete Gamxzz;
  delete Gamyxx; delete Gamyxy; delete Gamyxz; delete Gamyyy; delete Gamyyz; delete Gamyzz;
  delete Gamzxx; delete Gamzxy; delete Gamzxz; delete Gamzyy; delete Gamzyz; delete Gamzzz;

  // Ricci
  delete Rxx; delete Rxy; delete Rxz; delete Ryy; delete Ryz; delete Rzz;

  // Psi4
  delete Rpsi4; delete Ipsi4;
  delete t1Rpsi4; delete t1Ipsi4;
  delete t2Rpsi4; delete t2Ipsi4;

  // Constraints
  delete Cons_Ham; delete Cons_Px; delete Cons_Py; delete Cons_Pz;
  delete Cons_Gx; delete Cons_Gy; delete Cons_Gz;
  delete Cons_fR;

  delete ErrorMonitor; delete Psi4Monitor; delete BHMonitor; delete MAPMonitor; delete ConVMonitor;
  delete MaxScalar_Monitor;
  delete Waveshell; delete CheckPoint;
  delete GH;
#ifdef WithShell
  delete SH;
#endif
}

//================================================================================================



//================================================================================================

// 该成员函数用于对 Class 初始化
// Initialize

//================================================================================================

void bssnEScalar_class::Initialize()
{
  if (myrank == 0)
    cout << "you have setted " << ngfs << " grid functions." << endl;

  CheckPoint->addvariablelist(StateList);
  CheckPoint->addvariablelist(OldStateList);

  GH = new cgh(0, ngfs, Symmetry, "input.par", checkrun, ErrorMonitor);
  if (checkrun)
    CheckPoint->readcheck_cgh(PhysTime, GH, myrank, nprocs, Symmetry);
  else
    GH->compose_cgh(nprocs);

#ifdef WithShell
  SH = new ShellPatch(0, ngfs, "input.par", Symmetry, myrank, ErrorMonitor);
  SH->matchcheck(GH->PatL[0]);
  SH->compose_sh(nprocs);
  SH->setupcordtrans();
  SH->Dump_xyz(0, 0, 1);
  SH->setupintintstuff(nprocs, GH->PatL[0], Symmetry);
#ifdef USE_GPU
SH->build_gpu_block_table(GH->PatL[0]);
SH->build_gpu_pointstru_cache(SH->gpu_ss_cache, SH->ss_src, SH->ss_dst, nprocs, GH->PatL[0]);
SH->build_gpu_pointstru_cache(SH->gpu_csatc_cache, SH->csatc_src, SH->csatc_dst, nprocs, GH->PatL[0]);
SH->build_gpu_pointstru_cache(SH->gpu_csats_cache, SH->csats_src, SH->csats_dst, nprocs, GH->PatL[0]);
#endif

  if (checkrun)
    CheckPoint->readcheck_sh(SH, myrank);
#else
  SH = 0;
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



// NOTE: bssnEScalar_class::Step_GPU was moved to bssnEScalar_step_gpu.C
// to mirror the bssn_step_gpu.C / bssn_gpu_class.C split for the base class.


//================================================================================================
//
// Constraint_Out — GPU implementation
//
// Mirrors bssn_gpu_class.C::Constraint_Out (base BSSN). The CPU counterpart
// lives in bssnEScalar_class.C and is selected by `#ifndef USE_GPU`; this
// GPU body is always compiled into bssnEScalar_gpu_class.o and only ABEGPU
// links it (see makefile C++FILES_GPU).
//
// Strategy:
//   * For lev>0 re-run the EScalar RHS kernel on device so Cons_*, rho/S*
//     get refreshed on moved-grid ghost zones (mirrors what the CPU path did
//     with f_compute_rhs_bssn_escalar).
//   * For lev==0 trust the Cons_* the last Step_GPU RHS launch wrote.
//   * On every level compute the EScalar-only Cons_fR via
//     gpu_compute_constraint_fr_launch (4D Ricci -> f'(R) -> finalize).
//   * Cross-rank L2 reductions use Parallel::L2Norm_GPU.
//
// No StateList GPU->CPU sync, no CPU f_compute_rhs_bssn_escalar call, no
// CPU f_compute_constraint_fr. The previous fallback was the largest single
// host/device round-trip in the EScalar pipeline and accounted for the ~8s
// "Before Step" gap that fragmented GPU utilisation.
//
//================================================================================================



//================================================================================================

// 该成员函数用于计算和输出约束违反

//================================================================================================

void bssnEScalar_class::Constraint_Out()
{
    LastConsOut += dT * pow(0.5, Mymax(0, trfls));

    if (LastConsOut < AnasTime)
        return;

    if (use_gpu != 1) {
        // Defensive fallback: a GPU-built binary asked to run a CPU-only case
        // (use_gpu==0). The CPU bssnEScalar_class::Constraint_Out is excluded
        // by `#ifndef USE_GPU` in bssnEScalar_class.C, so we cannot delegate
        // — log and bail. ASC26 always runs use_gpu==1; if you hit this,
        // either rebuild ABE or extend this branch.
        if (myrank == 0 && ErrorMonitor && ErrorMonitor->outfile) {
            ErrorMonitor->outfile << "bssnEScalar_class::Constraint_Out (GPU TU) "
                                  << "called with use_gpu=0; CPU path not "
                                  << "linked in this binary." << endl;
        }
        return;
    }

    double TRK4 = PhysTime;
    double ndeps = numepsb;
    int pre = 0;

    for (int lev = 0; lev < GH->levels; lev++) {
        MyList<Patch> *Pp = GH->PatL[lev];
        while (Pp) {
            MyList<Block> *BP = Pp->data->blb;
            while (BP) {
                Block *cg = BP->data;
                if (myrank == cg->rank) {
                    if (lev > 0) {
                        // Refresh BSSN constraints + matter source on moved grid.
                        gpu_apply_symmetry_bd_fused_launch(
                            cg->stream, cg->shape[0], cg->shape[1], cg->shape[2],
                            cg->shape_pad[0], cg->shape_pad[1], cg->shape_pad[2], Symmetry,
                            cg->d_fgfs[phi0->sgfn], cg->d_fgfs[trK0->sgfn],
                            cg->d_fgfs[gxx0->sgfn], cg->d_fgfs[gxy0->sgfn], cg->d_fgfs[gxz0->sgfn],
                            cg->d_fgfs[gyy0->sgfn], cg->d_fgfs[gyz0->sgfn], cg->d_fgfs[gzz0->sgfn],
                            cg->d_fgfs[Axx0->sgfn], cg->d_fgfs[Axy0->sgfn], cg->d_fgfs[Axz0->sgfn],
                            cg->d_fgfs[Ayy0->sgfn], cg->d_fgfs[Ayz0->sgfn], cg->d_fgfs[Azz0->sgfn],
                            cg->d_fgfs[Gmx0->sgfn], cg->d_fgfs[Gmy0->sgfn], cg->d_fgfs[Gmz0->sgfn],
                            cg->d_fgfs[Lap0->sgfn],
                            cg->d_fgfs[Sfx0->sgfn], cg->d_fgfs[Sfy0->sgfn], cg->d_fgfs[Sfz0->sgfn],
                            cg->d_fgfs[dtSfx0->sgfn], cg->d_fgfs[dtSfy0->sgfn], cg->d_fgfs[dtSfz0->sgfn]
                        );
                        gpu_apply_scalar_symmetry_launch(
                            cg->stream, cg->shape[0], cg->shape[1], cg->shape[2],
                            cg->shape_pad[0], cg->shape_pad[1], cg->shape_pad[2],
                            cg->d_fgfs[Sphi0->sgfn], cg->d_fgfs[Spi0->sgfn]
                        );
                        gpu_compute_rhs_bssn_escalar_launch(
                            cg->stream,
                            cg->shape, cg->shape_pad, TRK4, cg->d_X[0], cg->d_X[1], cg->d_X[2],
                            cg->d_fgfs[phi0->sgfn], cg->d_fgfs[trK0->sgfn],
                            cg->d_fgfs[gxx0->sgfn], cg->d_fgfs[gxy0->sgfn], cg->d_fgfs[gxz0->sgfn],
                            cg->d_fgfs[gyy0->sgfn], cg->d_fgfs[gyz0->sgfn], cg->d_fgfs[gzz0->sgfn],
                            cg->d_fgfs[Axx0->sgfn], cg->d_fgfs[Axy0->sgfn], cg->d_fgfs[Axz0->sgfn],
                            cg->d_fgfs[Ayy0->sgfn], cg->d_fgfs[Ayz0->sgfn], cg->d_fgfs[Azz0->sgfn],
                            cg->d_fgfs[Gmx0->sgfn], cg->d_fgfs[Gmy0->sgfn], cg->d_fgfs[Gmz0->sgfn],
                            cg->d_fgfs[Lap0->sgfn],
                            cg->d_fgfs[Sfx0->sgfn], cg->d_fgfs[Sfy0->sgfn], cg->d_fgfs[Sfz0->sgfn],
                            cg->d_fgfs[dtSfx0->sgfn], cg->d_fgfs[dtSfy0->sgfn], cg->d_fgfs[dtSfz0->sgfn],
                            cg->d_fgfs[Sphi0->sgfn], cg->d_fgfs[Spi0->sgfn],
                            cg->d_fgfs[phi_rhs->sgfn], cg->d_fgfs[trK_rhs->sgfn],
                            cg->d_fgfs[gxx_rhs->sgfn], cg->d_fgfs[gxy_rhs->sgfn], cg->d_fgfs[gxz_rhs->sgfn],
                            cg->d_fgfs[gyy_rhs->sgfn], cg->d_fgfs[gyz_rhs->sgfn], cg->d_fgfs[gzz_rhs->sgfn],
                            cg->d_fgfs[Axx_rhs->sgfn], cg->d_fgfs[Axy_rhs->sgfn], cg->d_fgfs[Axz_rhs->sgfn],
                            cg->d_fgfs[Ayy_rhs->sgfn], cg->d_fgfs[Ayz_rhs->sgfn], cg->d_fgfs[Azz_rhs->sgfn],
                            cg->d_fgfs[Gmx_rhs->sgfn], cg->d_fgfs[Gmy_rhs->sgfn], cg->d_fgfs[Gmz_rhs->sgfn],
                            cg->d_fgfs[Lap_rhs->sgfn],
                            cg->d_fgfs[Sfx_rhs->sgfn], cg->d_fgfs[Sfy_rhs->sgfn], cg->d_fgfs[Sfz_rhs->sgfn],
                            cg->d_fgfs[dtSfx_rhs->sgfn], cg->d_fgfs[dtSfy_rhs->sgfn], cg->d_fgfs[dtSfz_rhs->sgfn],
                            cg->d_fgfs[Sphi_rhs->sgfn], cg->d_fgfs[Spi_rhs->sgfn],
                            cg->d_fgfs[rho->sgfn], cg->d_fgfs[Sx->sgfn], cg->d_fgfs[Sy->sgfn], cg->d_fgfs[Sz->sgfn],
                            cg->d_fgfs[Sxx->sgfn], cg->d_fgfs[Sxy->sgfn], cg->d_fgfs[Sxz->sgfn],
                            cg->d_fgfs[Syy->sgfn], cg->d_fgfs[Syz->sgfn], cg->d_fgfs[Szz->sgfn],
                            cg->d_fgfs[Gamxxx->sgfn], cg->d_fgfs[Gamxxy->sgfn], cg->d_fgfs[Gamxxz->sgfn],
                            cg->d_fgfs[Gamxyy->sgfn], cg->d_fgfs[Gamxyz->sgfn], cg->d_fgfs[Gamxzz->sgfn],
                            cg->d_fgfs[Gamyxx->sgfn], cg->d_fgfs[Gamyxy->sgfn], cg->d_fgfs[Gamyxz->sgfn],
                            cg->d_fgfs[Gamyyy->sgfn], cg->d_fgfs[Gamyyz->sgfn], cg->d_fgfs[Gamyzz->sgfn],
                            cg->d_fgfs[Gamzxx->sgfn], cg->d_fgfs[Gamzxy->sgfn], cg->d_fgfs[Gamzxz->sgfn],
                            cg->d_fgfs[Gamzyy->sgfn], cg->d_fgfs[Gamzyz->sgfn], cg->d_fgfs[Gamzzz->sgfn],
                            cg->d_fgfs[Rxx->sgfn], cg->d_fgfs[Rxy->sgfn], cg->d_fgfs[Rxz->sgfn],
                            cg->d_fgfs[Ryy->sgfn], cg->d_fgfs[Ryz->sgfn], cg->d_fgfs[Rzz->sgfn],
                            cg->d_fgfs[Cons_Ham->sgfn],
                            cg->d_fgfs[Cons_Px->sgfn], cg->d_fgfs[Cons_Py->sgfn], cg->d_fgfs[Cons_Pz->sgfn],
                            cg->d_fgfs[Cons_Gx->sgfn], cg->d_fgfs[Cons_Gy->sgfn], cg->d_fgfs[Cons_Gz->sgfn],
                            Symmetry, lev, ndeps, pre
                        );
                    }

                    // Cons_fR is EScalar-only — Step_GPU doesn't write it.
                    // Recompute on every level. Pooled scratch alloc is reused
                    // across patches (free-list shuffling, no actual cudaFree).
                    const size_t pad_n = (size_t)cg->shape_pad[0]
                                       * (size_t)cg->shape_pad[1]
                                       * (size_t)cg->shape_pad[2];
                    double *d_scratch_RR = GPUManager::getInstance()
                        .allocate_device_memory<double>(pad_n, cg->stream);
                    gpu_compute_constraint_fr_launch(
                        cg->stream, cg->shape, cg->shape_pad,
                        cg->d_fgfs[phi0->sgfn], cg->d_fgfs[trK0->sgfn],
                        cg->d_fgfs[rho->sgfn], cg->d_fgfs[Sphi0->sgfn],
                        cg->d_fgfs[gxx0->sgfn], cg->d_fgfs[gxy0->sgfn], cg->d_fgfs[gxz0->sgfn],
                        cg->d_fgfs[gyy0->sgfn], cg->d_fgfs[gyz0->sgfn], cg->d_fgfs[gzz0->sgfn],
                        cg->d_fgfs[Axx0->sgfn], cg->d_fgfs[Axy0->sgfn], cg->d_fgfs[Axz0->sgfn],
                        cg->d_fgfs[Ayy0->sgfn], cg->d_fgfs[Ayz0->sgfn], cg->d_fgfs[Azz0->sgfn],
                        cg->d_fgfs[Rxx->sgfn], cg->d_fgfs[Rxy->sgfn], cg->d_fgfs[Rxz->sgfn],
                        cg->d_fgfs[Ryy->sgfn], cg->d_fgfs[Ryz->sgfn], cg->d_fgfs[Rzz->sgfn],
                        cg->d_fgfs[Sxx->sgfn], cg->d_fgfs[Sxy->sgfn], cg->d_fgfs[Sxz->sgfn],
                        cg->d_fgfs[Syy->sgfn], cg->d_fgfs[Syz->sgfn], cg->d_fgfs[Szz->sgfn],
                        cg->d_fgfs[Cons_fR->sgfn], d_scratch_RR
                    );
                    GPUManager::getInstance().free_device_memory(d_scratch_RR, pad_n, cg->stream);
                }
                if (BP == Pp->data->ble) break;
                BP = BP->next;
            }
            Pp = Pp->next;
        }
        Parallel::Sync_GPU(GH->PatL[lev], ConstraintList, Symmetry);
    }

    double ConV[8];
    for (int levi = 0; levi < GH->levels; levi++) {
        ConV[0] = Parallel::L2Norm_GPU(GH->PatL[levi]->data, Cons_Ham);
        ConV[1] = Parallel::L2Norm_GPU(GH->PatL[levi]->data, Cons_Px);
        ConV[2] = Parallel::L2Norm_GPU(GH->PatL[levi]->data, Cons_Py);
        ConV[3] = Parallel::L2Norm_GPU(GH->PatL[levi]->data, Cons_Pz);
        ConV[4] = Parallel::L2Norm_GPU(GH->PatL[levi]->data, Cons_Gx);
        ConV[5] = Parallel::L2Norm_GPU(GH->PatL[levi]->data, Cons_Gy);
        ConV[6] = Parallel::L2Norm_GPU(GH->PatL[levi]->data, Cons_Gz);
        ConV[7] = Parallel::L2Norm_GPU(GH->PatL[levi]->data, Cons_fR);
        ConVMonitor->writefile(PhysTime, 8, ConV);
    }

    LastConsOut = 0;
}

//================================================================================================
