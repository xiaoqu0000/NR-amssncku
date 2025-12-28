
#ifdef newc
#include <typeinfo>
#include <sstream>
#include <cstdio>
#include <map>
#include <string>
#include <cstring> 
#include <iostream>
using namespace std;
#else
#include <stdio.h>
#include <map.h>
#include <string.h>
#endif

#include <time.h>

#include "macrodef.h"
#include "misc.h"
#include "Ansorg.h"
#include "fmisc.h"
#include "Parallel.h"
#include "bssn_class.h"
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

#include "perf.h"

#include "derivatives.h"
#include "ricci_gamma.h"

//================================================================================================

// 定义 bssn_class

//================================================================================================

bssn_class::bssn_class(double Couranti, double StartTimei, double TotalTimei, 
                       double DumpTimei, double d2DumpTimei, double CheckTimei, double AnasTimei,
                       int Symmetryi, int checkruni, char *checkfilenamei, 
                       double numepssi, double numepsbi, double numepshi,
                       int a_levi, int maxli, int decni, double maxrexi, double drexi) 
                       : Courant(Couranti), StartTime(StartTimei), TotalTime(TotalTimei), 
                         DumpTime(DumpTimei), d2DumpTime(d2DumpTimei), CheckTime(CheckTimei), AnasTime(AnasTimei),
                         Symmetry(Symmetryi), checkrun(checkruni), numepss(numepssi), numepsb(numepsbi), numepsh(numepshi),
#ifdef With_AHF
                         xc(0), yc(0), zc(0), xr(0), yr(0), zr(0), trigger(0), dTT(0), dumpid(0),
#endif
                         a_lev(a_levi), maxl(maxli), decn(decni), maxrex(maxrexi), drex(drexi),
                         CheckPoint(0)
                         // CheckPoint(0)
{
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  // setup Monitors
  {
    stringstream a_stream;
    a_stream.setf(ios::left);
    a_stream << "# Error log information";
    ErrorMonitor = new monitor("Error.log", myrank, a_stream.str());
    ErrorMonitor->print_message("Warning: we always assume intput parameter in cell center style.");

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
    Psi4Monitor = new monitor("bssn_psi4.dat", myrank, a_stream.str());

    a_stream.clear();
    a_stream.str("");
    a_stream << setw(15) << "# time";
    BHMonitor = new monitor("bssn_BH.dat", myrank, a_stream.str());

    a_stream.clear();
    a_stream.str("");
    a_stream << setw(15) << "# time ADMmass ADMPx ADMPy ADMPz ADMSx ADMSy ADMSz";
    MAPMonitor = new monitor("bssn_ADMQs.dat", myrank, a_stream.str());

    a_stream.clear();
    a_stream.str("");
    a_stream << setw(15) << "# time Ham Px Py Pz Gx Gy Gz";
    ConVMonitor = new monitor("bssn_constraint.dat", myrank, a_stream.str());
  }
  // setup sphere integration engine
  Waveshell = new surface_integral(Symmetry);

  trfls = 0;
  chitiny = 0;
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
    // echo information of lower bound of chi
    cout << " chitiny = " << chitiny << endl;
    cout << " time refinement start from level #" << trfls << endl;
#ifdef With_AHF
    cout << " parameters for AHF:" << endl;
    cout << " AHfindevery = " << AHfindevery << endl;
    cout << " AHdumptime = " << AHdumptime << endl;
#endif
  }

  chitiny = chitiny - 1; // because we have subtracted one from chi

  strcpy(checkfilename, checkfilenamei);

  ngfs = 0;
  phio = new var("phio", ngfs++, 1, 1, 1);
  trKo = new var("trKo", ngfs++, 1, 1, 1);
  gxxo = new var("gxxo", ngfs++, 1, 1, 1);
  gxyo = new var("gxyo", ngfs++, -1, -1, 1);
  gxzo = new var("gxzo", ngfs++, -1, 1, -1);
  gyyo = new var("gyyo", ngfs++, 1, 1, 1);
  gyzo = new var("gyzo", ngfs++, 1, -1, -1);
  gzzo = new var("gzzo", ngfs++, 1, 1, 1);
  Axxo = new var("Axxo", ngfs++, 1, 1, 1);
  Axyo = new var("Axyo", ngfs++, -1, -1, 1);
  Axzo = new var("Axzo", ngfs++, -1, 1, -1);
  Ayyo = new var("Ayyo", ngfs++, 1, 1, 1);
  Ayzo = new var("Ayzo", ngfs++, 1, -1, -1);
  Azzo = new var("Azzo", ngfs++, 1, 1, 1);
  Gmxo = new var("Gmxo", ngfs++, -1, 1, 1);
  Gmyo = new var("Gmyo", ngfs++, 1, -1, 1);
  Gmzo = new var("Gmzo", ngfs++, 1, 1, -1);
  Lapo = new var("Lapo", ngfs++, 1, 1, 1);
  Sfxo = new var("Sfxo", ngfs++, -1, 1, 1);
  Sfyo = new var("Sfyo", ngfs++, 1, -1, 1);
  Sfzo = new var("Sfzo", ngfs++, 1, 1, -1);
  dtSfxo = new var("dtSfxo", ngfs++, -1, 1, 1);
  dtSfyo = new var("dtSfyo", ngfs++, 1, -1, 1);
  dtSfzo = new var("dtSfzo", ngfs++, 1, 1, -1);

  phi0 = new var("phi0", ngfs++, 1, 1, 1);
  trK0 = new var("trK0", ngfs++, 1, 1, 1);
  gxx0 = new var("gxx0", ngfs++, 1, 1, 1);
  gxy0 = new var("gxy0", ngfs++, -1, -1, 1);
  gxz0 = new var("gxz0", ngfs++, -1, 1, -1);
  gyy0 = new var("gyy0", ngfs++, 1, 1, 1);
  gyz0 = new var("gyz0", ngfs++, 1, -1, -1);
  gzz0 = new var("gzz0", ngfs++, 1, 1, 1);
  Axx0 = new var("Axx0", ngfs++, 1, 1, 1);
  Axy0 = new var("Axy0", ngfs++, -1, -1, 1);
  Axz0 = new var("Axz0", ngfs++, -1, 1, -1);
  Ayy0 = new var("Ayy0", ngfs++, 1, 1, 1);
  Ayz0 = new var("Ayz0", ngfs++, 1, -1, -1);
  Azz0 = new var("Azz0", ngfs++, 1, 1, 1);
  Gmx0 = new var("Gmx0", ngfs++, -1, 1, 1);
  Gmy0 = new var("Gmy0", ngfs++, 1, -1, 1);
  Gmz0 = new var("Gmz0", ngfs++, 1, 1, -1);
  Lap0 = new var("Lap0", ngfs++, 1, 1, 1);
  Sfx0 = new var("Sfx0", ngfs++, -1, 1, 1);
  Sfy0 = new var("Sfy0", ngfs++, 1, -1, 1);
  Sfz0 = new var("Sfz0", ngfs++, 1, 1, -1);
  dtSfx0 = new var("dtSfx0", ngfs++, -1, 1, 1);
  dtSfy0 = new var("dtSfy0", ngfs++, 1, -1, 1);
  dtSfz0 = new var("dtSfz0", ngfs++, 1, 1, -1);

  phi = new var("phi", ngfs++, 1, 1, 1);
  trK = new var("trK", ngfs++, 1, 1, 1);
  gxx = new var("gxx", ngfs++, 1, 1, 1);
  gxy = new var("gxy", ngfs++, -1, -1, 1);
  gxz = new var("gxz", ngfs++, -1, 1, -1);
  gyy = new var("gyy", ngfs++, 1, 1, 1);
  gyz = new var("gyz", ngfs++, 1, -1, -1);
  gzz = new var("gzz", ngfs++, 1, 1, 1);
  Axx = new var("Axx", ngfs++, 1, 1, 1);
  Axy = new var("Axy", ngfs++, -1, -1, 1);
  Axz = new var("Axz", ngfs++, -1, 1, -1);
  Ayy = new var("Ayy", ngfs++, 1, 1, 1);
  Ayz = new var("Ayz", ngfs++, 1, -1, -1);
  Azz = new var("Azz", ngfs++, 1, 1, 1);
  Gmx = new var("Gmx", ngfs++, -1, 1, 1);
  Gmy = new var("Gmy", ngfs++, 1, -1, 1);
  Gmz = new var("Gmz", ngfs++, 1, 1, -1);
  Lap = new var("Lap", ngfs++, 1, 1, 1);
  Sfx = new var("Sfx", ngfs++, -1, 1, 1);
  Sfy = new var("Sfy", ngfs++, 1, -1, 1);
  Sfz = new var("Sfz", ngfs++, 1, 1, -1);
  dtSfx = new var("dtSfx", ngfs++, -1, 1, 1);
  dtSfy = new var("dtSfy", ngfs++, 1, -1, 1);
  dtSfz = new var("dtSfz", ngfs++, 1, 1, -1);

  phi1 = new var("phi1", ngfs++, 1, 1, 1);
  trK1 = new var("trK1", ngfs++, 1, 1, 1);
  gxx1 = new var("gxx1", ngfs++, 1, 1, 1);
  gxy1 = new var("gxy1", ngfs++, -1, -1, 1);
  gxz1 = new var("gxz1", ngfs++, -1, 1, -1);
  gyy1 = new var("gyy1", ngfs++, 1, 1, 1);
  gyz1 = new var("gyz1", ngfs++, 1, -1, -1);
  gzz1 = new var("gzz1", ngfs++, 1, 1, 1);
  Axx1 = new var("Axx1", ngfs++, 1, 1, 1);
  Axy1 = new var("Axy1", ngfs++, -1, -1, 1);
  Axz1 = new var("Axz1", ngfs++, -1, 1, -1);
  Ayy1 = new var("Ayy1", ngfs++, 1, 1, 1);
  Ayz1 = new var("Ayz1", ngfs++, 1, -1, -1);
  Azz1 = new var("Azz1", ngfs++, 1, 1, 1);
  Gmx1 = new var("Gmx1", ngfs++, -1, 1, 1);
  Gmy1 = new var("Gmy1", ngfs++, 1, -1, 1);
  Gmz1 = new var("Gmz1", ngfs++, 1, 1, -1);
  Lap1 = new var("Lap1", ngfs++, 1, 1, 1);
  Sfx1 = new var("Sfx1", ngfs++, -1, 1, 1);
  Sfy1 = new var("Sfy1", ngfs++, 1, -1, 1);
  Sfz1 = new var("Sfz1", ngfs++, 1, 1, -1);
  dtSfx1 = new var("dtSfx1", ngfs++, -1, 1, 1);
  dtSfy1 = new var("dtSfy1", ngfs++, 1, -1, 1);
  dtSfz1 = new var("dtSfz1", ngfs++, 1, 1, -1);

  phi_rhs = new var("phi_rhs", ngfs++, 1, 1, 1);
  trK_rhs = new var("trK_rhs", ngfs++, 1, 1, 1);
  gxx_rhs = new var("gxx_rhs", ngfs++, 1, 1, 1);
  gxy_rhs = new var("gxy_rhs", ngfs++, -1, -1, 1);
  gxz_rhs = new var("gxz_rhs", ngfs++, -1, 1, -1);
  gyy_rhs = new var("gyy_rhs", ngfs++, 1, 1, 1);
  gyz_rhs = new var("gyz_rhs", ngfs++, 1, -1, -1);
  gzz_rhs = new var("gzz_rhs", ngfs++, 1, 1, 1);
  Axx_rhs = new var("Axx_rhs", ngfs++, 1, 1, 1);
  Axy_rhs = new var("Axy_rhs", ngfs++, -1, -1, 1);
  Axz_rhs = new var("Axz_rhs", ngfs++, -1, 1, -1);
  Ayy_rhs = new var("Ayy_rhs", ngfs++, 1, 1, 1);
  Ayz_rhs = new var("Ayz_rhs", ngfs++, 1, -1, -1);
  Azz_rhs = new var("Azz_rhs", ngfs++, 1, 1, 1);
  Gmx_rhs = new var("Gmx_rhs", ngfs++, -1, 1, 1);
  Gmy_rhs = new var("Gmy_rhs", ngfs++, 1, -1, 1);
  Gmz_rhs = new var("Gmz_rhs", ngfs++, 1, 1, -1);
  Lap_rhs = new var("Lap_rhs", ngfs++, 1, 1, 1);
  Sfx_rhs = new var("Sfx_rhs", ngfs++, -1, 1, 1);
  Sfy_rhs = new var("Sfy_rhs", ngfs++, 1, -1, 1);
  Sfz_rhs = new var("Sfz_rhs", ngfs++, 1, 1, -1);
  dtSfx_rhs = new var("dtSfx_rhs", ngfs++, -1, 1, 1);
  dtSfy_rhs = new var("dtSfy_rhs", ngfs++, 1, -1, 1);
  dtSfz_rhs = new var("dtSfz_rhs", ngfs++, 1, 1, -1);

  rho = new var("rho", ngfs++, 1, 1, 1);
  Sx = new var("Sx", ngfs++, -1, 1, 1);
  Sy = new var("Sy", ngfs++, 1, -1, 1);
  Sz = new var("Sz", ngfs++, 1, 1, -1);
  Sxx = new var("Sxx", ngfs++, 1, 1, 1);
  Sxy = new var("Sxy", ngfs++, -1, -1, 1);
  Sxz = new var("Sxz", ngfs++, -1, 1, -1);
  Syy = new var("Syy", ngfs++, 1, 1, 1);
  Syz = new var("Syz", ngfs++, 1, -1, -1);
  Szz = new var("Szz", ngfs++, 1, 1, 1);

  Gamxxx = new var("Gamxxx", ngfs++, -1, 1, 1);
  Gamxxy = new var("Gamxxy", ngfs++, 1, -1, 1);
  Gamxxz = new var("Gamxxz", ngfs++, 1, 1, -1);
  Gamxyy = new var("Gamxyy", ngfs++, -1, 1, 1);
  Gamxyz = new var("Gamxyz", ngfs++, -1, -1, -1);
  Gamxzz = new var("Gamxzz", ngfs++, -1, 1, 1);
  Gamyxx = new var("Gamyxx", ngfs++, 1, -1, 1);
  Gamyxy = new var("Gamyxy", ngfs++, -1, 1, 1);
  Gamyxz = new var("Gamyxz", ngfs++, -1, -1, -1);
  Gamyyy = new var("Gamyyy", ngfs++, 1, -1, 1);
  Gamyyz = new var("Gamyyz", ngfs++, 1, 1, -1);
  Gamyzz = new var("Gamyzz", ngfs++, 1, -1, 1);
  Gamzxx = new var("Gamzxx", ngfs++, 1, 1, -1);
  Gamzxy = new var("Gamzxy", ngfs++, -1, -1, -1);
  Gamzxz = new var("Gamzxz", ngfs++, -1, 1, 1);
  Gamzyy = new var("Gamzyy", ngfs++, 1, 1, -1);
  Gamzyz = new var("Gamzyz", ngfs++, 1, -1, 1);
  Gamzzz = new var("Gamzzz", ngfs++, 1, 1, -1);

  Rxx = new var("Rxx", ngfs++, 1, 1, 1);
  Rxy = new var("Rxy", ngfs++, -1, -1, 1);
  Rxz = new var("Rxz", ngfs++, -1, 1, -1);
  Ryy = new var("Ryy", ngfs++, 1, 1, 1);
  Ryz = new var("Ryz", ngfs++, 1, -1, -1);
  Rzz = new var("Rzz", ngfs++, 1, 1, 1);

  // refer to PRD, 77, 024027 (2008)
  Rpsi4 = new var("Rpsi4", ngfs++, 1, 1, 1);
  Ipsi4 = new var("Ipsi4", ngfs++, -1, -1, -1);
  t1Rpsi4 = new var("t1Rpsi4", ngfs++, 1, 1, 1);
  t1Ipsi4 = new var("t1Ipsi4", ngfs++, -1, -1, -1);
  t2Rpsi4 = new var("t2Rpsi4", ngfs++, 1, 1, 1);
  t2Ipsi4 = new var("t2Ipsi4", ngfs++, -1, -1, -1);

  // constraint violation monitor variables
  Cons_Ham = new var("Cons_Ham", ngfs++, 1, 1, 1);
  Cons_Px = new var("Cons_Px", ngfs++, -1, 1, 1);
  Cons_Py = new var("Cons_Py", ngfs++, 1, -1, 1);
  Cons_Pz = new var("Cons_Pz", ngfs++, 1, 1, -1);
  Cons_Gx = new var("Cons_Gx", ngfs++, -1, 1, 1);
  Cons_Gy = new var("Cons_Gy", ngfs++, 1, -1, 1);
  Cons_Gz = new var("Cons_Gz", ngfs++, 1, 1, -1);

#ifdef Point_Psi4
  phix = new var("phix", ngfs++, -1, 1, 1);
  phiy = new var("phiy", ngfs++, 1, -1, 1);
  phiz = new var("phiz", ngfs++, 1, 1, -1);
  trKx = new var("trKx", ngfs++, -1, 1, 1);
  trKy = new var("trKy", ngfs++, 1, -1, 1);
  trKz = new var("trKz", ngfs++, 1, 1, -1);
  Axxx = new var("Axxx", ngfs++, -1, 1, 1);
  Axxy = new var("Axxy", ngfs++, 1, -1, 1);
  Axxz = new var("Axxz", ngfs++, 1, 1, -1);
  Axyx = new var("Axyx", ngfs++, 1, -1, 1);
  Axyy = new var("Axyy", ngfs++, -1, 1, 1);
  Axyz = new var("Axyz", ngfs++, -1, -1, -1);
  Axzx = new var("Axzx", ngfs++, 1, 1, -1);
  Axzy = new var("Axzy", ngfs++, -1, -1, -1);
  Axzz = new var("Axzz", ngfs++, -1, 1, 1);
  Ayyx = new var("Ayyx", ngfs++, -1, 1, 1);
  Ayyy = new var("Ayyy", ngfs++, 1, -1, 1);
  Ayyz = new var("Ayyz", ngfs++, 1, 1, -1);
  Ayzx = new var("Ayzx", ngfs++, -1, -1, -1);
  Ayzy = new var("Ayzy", ngfs++, 1, 1, -1);
  Ayzz = new var("Ayzz", ngfs++, 1, -1, 1);
  Azzx = new var("Azzx", ngfs++, -1, 1, 1);
  Azzy = new var("Azzy", ngfs++, 1, -1, 1);
  Azzz = new var("Azzz", ngfs++, 1, 1, -1);
#endif

  // specific properspeed for 1+log slice
  {
    const double vl = sqrt(2);
    trKo->setpropspeed(vl);
    trK0->setpropspeed(vl);
    trK->setpropspeed(vl);
    trK1->setpropspeed(vl);
    trK_rhs->setpropspeed(vl);

    phio->setpropspeed(vl);
    phi0->setpropspeed(vl);
    phi->setpropspeed(vl);
    phi1->setpropspeed(vl);
    phi_rhs->setpropspeed(vl);

    Lapo->setpropspeed(vl);
    Lap0->setpropspeed(vl);
    Lap->setpropspeed(vl);
    Lap1->setpropspeed(vl);
    Lap_rhs->setpropspeed(vl);
  }

  OldStateList = new MyList<var>(phio);
  OldStateList->insert(trKo);
  OldStateList->insert(gxxo);
  OldStateList->insert(gxyo);
  OldStateList->insert(gxzo);
  OldStateList->insert(gyyo);
  OldStateList->insert(gyzo);
  OldStateList->insert(gzzo);
  OldStateList->insert(Axxo);
  OldStateList->insert(Axyo);
  OldStateList->insert(Axzo);
  OldStateList->insert(Ayyo);
  OldStateList->insert(Ayzo);
  OldStateList->insert(Azzo);
  OldStateList->insert(Gmxo);
  OldStateList->insert(Gmyo);
  OldStateList->insert(Gmzo);
  OldStateList->insert(Lapo);
  OldStateList->insert(Sfxo);
  OldStateList->insert(Sfyo);
  OldStateList->insert(Sfzo);
  OldStateList->insert(dtSfxo);
  OldStateList->insert(dtSfyo);
  OldStateList->insert(dtSfzo);

  StateList = new MyList<var>(phi0);
  StateList->insert(trK0);
  StateList->insert(gxx0);
  StateList->insert(gxy0);
  StateList->insert(gxz0);
  StateList->insert(gyy0);
  StateList->insert(gyz0);
  StateList->insert(gzz0);
  StateList->insert(Axx0);
  StateList->insert(Axy0);
  StateList->insert(Axz0);
  StateList->insert(Ayy0);
  StateList->insert(Ayz0);
  StateList->insert(Azz0);
  StateList->insert(Gmx0);
  StateList->insert(Gmy0);
  StateList->insert(Gmz0);
  StateList->insert(Lap0);
  StateList->insert(Sfx0);
  StateList->insert(Sfy0);
  StateList->insert(Sfz0);
  StateList->insert(dtSfx0);
  StateList->insert(dtSfy0);
  StateList->insert(dtSfz0);

  RHSList = new MyList<var>(phi_rhs);
  RHSList->insert(trK_rhs);
  RHSList->insert(gxx_rhs);
  RHSList->insert(gxy_rhs);
  RHSList->insert(gxz_rhs);
  RHSList->insert(gyy_rhs);
  RHSList->insert(gyz_rhs);
  RHSList->insert(gzz_rhs);
  RHSList->insert(Axx_rhs);
  RHSList->insert(Axy_rhs);
  RHSList->insert(Axz_rhs);
  RHSList->insert(Ayy_rhs);
  RHSList->insert(Ayz_rhs);
  RHSList->insert(Azz_rhs);
  RHSList->insert(Gmx_rhs);
  RHSList->insert(Gmy_rhs);
  RHSList->insert(Gmz_rhs);
  RHSList->insert(Lap_rhs);
  RHSList->insert(Sfx_rhs);
  RHSList->insert(Sfy_rhs);
  RHSList->insert(Sfz_rhs);
  RHSList->insert(dtSfx_rhs);
  RHSList->insert(dtSfy_rhs);
  RHSList->insert(dtSfz_rhs);

  SynchList_pre = new MyList<var>(phi);
  SynchList_pre->insert(trK);
  SynchList_pre->insert(gxx);
  SynchList_pre->insert(gxy);
  SynchList_pre->insert(gxz);
  SynchList_pre->insert(gyy);
  SynchList_pre->insert(gyz);
  SynchList_pre->insert(gzz);
  SynchList_pre->insert(Axx);
  SynchList_pre->insert(Axy);
  SynchList_pre->insert(Axz);
  SynchList_pre->insert(Ayy);
  SynchList_pre->insert(Ayz);
  SynchList_pre->insert(Azz);
  SynchList_pre->insert(Gmx);
  SynchList_pre->insert(Gmy);
  SynchList_pre->insert(Gmz);
  SynchList_pre->insert(Lap);
  SynchList_pre->insert(Sfx);
  SynchList_pre->insert(Sfy);
  SynchList_pre->insert(Sfz);
  SynchList_pre->insert(dtSfx);
  SynchList_pre->insert(dtSfy);
  SynchList_pre->insert(dtSfz);

  SynchList_cor = new MyList<var>(phi1);
  SynchList_cor->insert(trK1);
  SynchList_cor->insert(gxx1);
  SynchList_cor->insert(gxy1);
  SynchList_cor->insert(gxz1);
  SynchList_cor->insert(gyy1);
  SynchList_cor->insert(gyz1);
  SynchList_cor->insert(gzz1);
  SynchList_cor->insert(Axx1);
  SynchList_cor->insert(Axy1);
  SynchList_cor->insert(Axz1);
  SynchList_cor->insert(Ayy1);
  SynchList_cor->insert(Ayz1);
  SynchList_cor->insert(Azz1);
  SynchList_cor->insert(Gmx1);
  SynchList_cor->insert(Gmy1);
  SynchList_cor->insert(Gmz1);
  SynchList_cor->insert(Lap1);
  SynchList_cor->insert(Sfx1);
  SynchList_cor->insert(Sfy1);
  SynchList_cor->insert(Sfz1);
  SynchList_cor->insert(dtSfx1);
  SynchList_cor->insert(dtSfy1);
  SynchList_cor->insert(dtSfz1);

  DumpList = new MyList<var>(phi0);
  DumpList->insert(trK0);
  DumpList->insert(gxx0);
  DumpList->insert(gxy0);
  DumpList->insert(gxz0);
  DumpList->insert(gyy0);
  DumpList->insert(gyz0);
  DumpList->insert(gzz0);
  // DumpList->insert(Axx0);
  // DumpList->insert(Axy0);
  // DumpList->insert(Axz0);
  // DumpList->insert(Ayy0);
  // DumpList->insert(Ayz0);
  // DumpList->insert(Azz0);
  // DumpList->insert(Gmx0);
  // DumpList->insert(Gmy0);
  // DumpList->insert(Gmz0);
  DumpList->insert(Lap0);
  // DumpList->insert(Sfx0);
  // DumpList->insert(Sfy0);
  // DumpList->insert(Sfz0);
  // DumpList->insert(dtSfx0);
  // DumpList->insert(dtSfy0);
  // DumpList->insert(dtSfz0);
  // DumpList->insert(Rpsi4);
  // DumpList->insert(Ipsi4);
  DumpList->insert(Cons_Ham);
  DumpList->insert(Cons_Px);
  DumpList->insert(Cons_Py);
  DumpList->insert(Cons_Pz);
  // DumpList->insert(Cons_Gx);
  // DumpList->insert(Cons_Gy);
  // DumpList->insert(Cons_Gz);

  ConstraintList = new MyList<var>(Cons_Ham);
  ConstraintList->insert(Cons_Px);
  ConstraintList->insert(Cons_Py);
  ConstraintList->insert(Cons_Pz);
  ConstraintList->insert(Cons_Gx);
  ConstraintList->insert(Cons_Gy);
  ConstraintList->insert(Cons_Gz);
#ifdef With_AHF
  // setup kinds of var list
  // List for AparentHorizonFinderDirect
  // special attension is payed to symmetry type
  // gij  gij,x  gij,y  gij,z
  AHList = new MyList<var>(gxx0);
  AHList->insert(Gamxxx);
  AHList->insert(Gamyxx);
  AHList->insert(Gamzxx);
  AHList->insert(gxy0);
  AHList->insert(Gamxxy);
  AHList->insert(Gamyxy);
  AHList->insert(Gamzxy);
  AHList->insert(gxz0);
  AHList->insert(Gamxxz);
  AHList->insert(Gamyxz);
  AHList->insert(Gamzxz);
  AHList->insert(gyy0);
  AHList->insert(Gamxyy);
  AHList->insert(Gamyyy);
  AHList->insert(Gamzyy);
  AHList->insert(gyz0);
  AHList->insert(Gamxyz);
  AHList->insert(Gamyyz);
  AHList->insert(Gamzyz);
  AHList->insert(gzz0);
  AHList->insert(Gamxzz);
  AHList->insert(Gamyzz);
  AHList->insert(Gamzzz);
  // phi  phi,x  phi,y  phi,z
  AHList->insert(phi0);
  AHList->insert(dtSfx_rhs);
  AHList->insert(dtSfy_rhs);
  AHList->insert(dtSfz_rhs);
  // Aij
  AHList->insert(Axx0);
  AHList->insert(Axy0);
  AHList->insert(Axz0);
  AHList->insert(Ayy0);
  AHList->insert(Ayz0);
  AHList->insert(Azz0);
  // trK
  AHList->insert(trK0);
  // gij,x  gij,y  gij,z
  AHDList = new MyList<var>(Gamxxx);
  AHDList->insert(Gamyxx);
  AHDList->insert(Gamzxx);
  AHDList->insert(Gamxxy);
  AHDList->insert(Gamyxy);
  AHDList->insert(Gamzxy);
  AHDList->insert(Gamxxz);
  AHDList->insert(Gamyxz);
  AHDList->insert(Gamzxz);
  AHDList->insert(Gamxyy);
  AHDList->insert(Gamyyy);
  AHDList->insert(Gamzyy);
  AHDList->insert(Gamxyz);
  AHDList->insert(Gamyyz);
  AHDList->insert(Gamzyz);
  AHDList->insert(Gamxzz);
  AHDList->insert(Gamyzz);
  AHDList->insert(Gamzzz);
  // phi,x  phi,y  phi,z
  AHDList->insert(dtSfx_rhs);
  AHDList->insert(dtSfy_rhs);
  AHDList->insert(dtSfz_rhs);

  GaugeList = new MyList<var>(Lap0);
  GaugeList->insert(Sfx0);
  GaugeList->insert(Sfy0);
  GaugeList->insert(Sfz0);
#endif
  
  // checkpoint class 在ubuntu24下遇到内存错误，以为是类型问题，结果发现不是
  // checkpoint class 中增加输出后，莫名其妙好了
  // 已查到原因，文件名宽度超出范围
  
  // checkpoint class 中的第一个变量类型为bool，而这里的变量类型是int，人为进行一次转换
  // bool checkrun00 = checkrun;
  // checkpoint class 中的第二个变量类型为const char，而这里的变量类型是char，人为进行一次转换
  // const char* checkfilename00 = checkfilename;

  CheckPoint = new checkpoint(checkrun, checkfilename, myrank);
  
  if (myrank==0) {
    cout << " BSSN class successfully created " << endl;
  }
}

//================================================================================================



//================================================================================================

// 该成员函数用于对 Class 初始化

//================================================================================================

void bssn_class::Initialize()
{
  if (myrank == 0)
    cout << " you have setted " << ngfs << " grid functions." << endl;

  CheckPoint->addvariablelist(StateList);
  CheckPoint->addvariablelist(OldStateList);

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
  SH->matchcheck(GH->PatL[0]);
  SH->compose_sh(nprocs);
  //   SH->compose_shr(nprocs);  //sh is faster than shr
  SH->setupcordtrans();
  SH->Dump_xyz(0, 0, 1);
  SH->setupintintstuff(nprocs, GH->PatL[0], Symmetry);

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
    setpbh(BH_num, Porg0, Mass, BH_num_input);
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

bssn_class::~bssn_class()
{
#ifdef With_AHF
  AHList->clearList();
  AHDList->clearList();
  GaugeList->clearList();
  if (lastahdumpid)
    delete[] lastahdumpid;
  if (findeveryl)
    delete[] findeveryl;

  if (xc)
  {
    delete[] xc;
    delete[] yc;
    delete[] zc;
    delete[] xr;
    delete[] yr;
    delete[] zr;
    delete[] trigger;
    delete[] dumpid;
    delete[] dTT;
  }

  AHFinderDirect::AHFinderDirect_cleanup();
#endif

  StateList->clearList();
  RHSList->clearList();
  OldStateList->clearList();
  SynchList_pre->clearList();
  SynchList_cor->clearList();
  DumpList->clearList();
  ConstraintList->clearList();

  delete phio;
  delete trKo;
  delete gxxo;
  delete gxyo;
  delete gxzo;
  delete gyyo;
  delete gyzo;
  delete gzzo;
  delete Axxo;
  delete Axyo;
  delete Axzo;
  delete Ayyo;
  delete Ayzo;
  delete Azzo;
  delete Gmxo;
  delete Gmyo;
  delete Gmzo;
  delete Lapo;
  delete Sfxo;
  delete Sfyo;
  delete Sfzo;
  delete dtSfxo;
  delete dtSfyo;
  delete dtSfzo;

  delete phi0;
  delete trK0;
  delete gxx0;
  delete gxy0;
  delete gxz0;
  delete gyy0;
  delete gyz0;
  delete gzz0;
  delete Axx0;
  delete Axy0;
  delete Axz0;
  delete Ayy0;
  delete Ayz0;
  delete Azz0;
  delete Gmx0;
  delete Gmy0;
  delete Gmz0;
  delete Lap0;
  delete Sfx0;
  delete Sfy0;
  delete Sfz0;
  delete dtSfx0;
  delete dtSfy0;
  delete dtSfz0;

  delete phi;
  delete trK;
  delete gxx;
  delete gxy;
  delete gxz;
  delete gyy;
  delete gyz;
  delete gzz;
  delete Axx;
  delete Axy;
  delete Axz;
  delete Ayy;
  delete Ayz;
  delete Azz;
  delete Gmx;
  delete Gmy;
  delete Gmz;
  delete Lap;
  delete Sfx;
  delete Sfy;
  delete Sfz;
  delete dtSfx;
  delete dtSfy;
  delete dtSfz;

  delete phi1;
  delete trK1;
  delete gxx1;
  delete gxy1;
  delete gxz1;
  delete gyy1;
  delete gyz1;
  delete gzz1;
  delete Axx1;
  delete Axy1;
  delete Axz1;
  delete Ayy1;
  delete Ayz1;
  delete Azz1;
  delete Gmx1;
  delete Gmy1;
  delete Gmz1;
  delete Lap1;
  delete Sfx1;
  delete Sfy1;
  delete Sfz1;
  delete dtSfx1;
  delete dtSfy1;
  delete dtSfz1;

  delete phi_rhs;
  delete trK_rhs;
  delete gxx_rhs;
  delete gxy_rhs;
  delete gxz_rhs;
  delete gyy_rhs;
  delete gyz_rhs;
  delete gzz_rhs;
  delete Axx_rhs;
  delete Axy_rhs;
  delete Axz_rhs;
  delete Ayy_rhs;
  delete Ayz_rhs;
  delete Azz_rhs;
  delete Gmx_rhs;
  delete Gmy_rhs;
  delete Gmz_rhs;
  delete Lap_rhs;
  delete Sfx_rhs;
  delete Sfy_rhs;
  delete Sfz_rhs;
  delete dtSfx_rhs;
  delete dtSfy_rhs;
  delete dtSfz_rhs;

  delete rho;
  delete Sx;
  delete Sy;
  delete Sz;
  delete Sxx;
  delete Sxy;
  delete Sxz;
  delete Syy;
  delete Syz;
  delete Szz;

  delete Gamxxx;
  delete Gamxxy;
  delete Gamxxz;
  delete Gamxyy;
  delete Gamxyz;
  delete Gamxzz;
  delete Gamyxx;
  delete Gamyxy;
  delete Gamyxz;
  delete Gamyyy;
  delete Gamyyz;
  delete Gamyzz;
  delete Gamzxx;
  delete Gamzxy;
  delete Gamzxz;
  delete Gamzyy;
  delete Gamzyz;
  delete Gamzzz;

  delete Rxx;
  delete Rxy;
  delete Rxz;
  delete Ryy;
  delete Ryz;
  delete Rzz;

  delete Rpsi4;
  delete Ipsi4;
  delete t1Rpsi4;
  delete t1Ipsi4;
  delete t2Rpsi4;
  delete t2Ipsi4;

  delete Cons_Ham;
  delete Cons_Px;
  delete Cons_Py;
  delete Cons_Pz;
  delete Cons_Gx;
  delete Cons_Gy;
  delete Cons_Gz;

#ifdef Point_Psi4
  delete phix;
  delete phiy;
  delete phiz;
  delete trKx;
  delete trKy;
  delete trKz;
  delete Axxx;
  delete Axxy;
  delete Axxz;
  delete Axyx;
  delete Axyy;
  delete Axyz;
  delete Axzx;
  delete Axzy;
  delete Axzz;
  delete Ayyx;
  delete Ayyy;
  delete Ayyz;
  delete Ayzx;
  delete Ayzy;
  delete Ayzz;
  delete Azzx;
  delete Azzy;
  delete Azzz;
#endif

  delete GH;
#ifdef WithShell
  delete SH;
#endif

  for (int i = 0; i < BH_num; i++)
  {
    delete[] Porg0[i];
    delete[] Porgbr[i];
    delete[] Porg[i];
    delete[] Porg1[i];
    delete[] Porg_rhs[i];
  }

  delete[] Porg0;
  delete[] Porgbr;
  delete[] Porg;
  delete[] Porg1;
  delete[] Porg_rhs;

  delete[] Mass;
  delete[] Spin;
  delete[] Pmom;

  delete ErrorMonitor;
  delete Psi4Monitor;
  delete BHMonitor;
  delete MAPMonitor;
  delete ConVMonitor;
  delete Waveshell;

  delete CheckPoint;
}

//================================================================================================



//================================================================================================

// 该成员函数用 Lousto 的解析方法求解初值

//================================================================================================

void bssn_class::Setup_Initial_Data_Lousto()
{
  if (!checkrun)
  {
    if (myrank == 0)
    {
      cout << endl;
      cout << " Setup initial data with Lousto's analytical formula. " << endl;
      cout << endl;
    }
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
    double *Porg_here, *Pmom_here, *Spin_here, *Mass_here;
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
            ErrorMonitor->outfile << "error reading parameter file " << filename << " in line " << i << endl;
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
    Pmom_here = new double[3 * BH_NM];
    Spin_here = new double[3 * BH_NM];
    Mass_here = new double[BH_NM];
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
            ErrorMonitor->outfile << "error reading parameter file " << filename << " in line " << i << endl;
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else if (status == 0)
          continue;

        if (sgrp == "BSSN" && sind < BH_NM)
        {
          if (skey == "Mass")
            Mass_here[sind] = atof(sval.c_str());
          else if (skey == "Porgx")
            Porg_here[sind * 3] = atof(sval.c_str());
          else if (skey == "Porgy")
            Porg_here[sind * 3 + 1] = atof(sval.c_str());
          else if (skey == "Porgz")
            Porg_here[sind * 3 + 2] = atof(sval.c_str());
          else if (skey == "Spinx")
            Spin_here[sind * 3] = atof(sval.c_str());
          else if (skey == "Spiny")
            Spin_here[sind * 3 + 1] = atof(sval.c_str());
          else if (skey == "Spinz")
            Spin_here[sind * 3 + 2] = atof(sval.c_str());
          else if (skey == "Pmomx")
            Pmom_here[sind * 3] = atof(sval.c_str());
          else if (skey == "Pmomy")
            Pmom_here[sind * 3 + 1] = atof(sval.c_str());
          else if (skey == "Pmomz")
            Pmom_here[sind * 3 + 2] = atof(sval.c_str());
        }
      }
      inf.close();
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
            // 用 Lousto 的解析公式求解初值
            f_get_lousto_nbhs(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                              cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                              cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                              cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                              cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                              cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                              cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                              cg->fgfs[Lap0->sgfn], 
                              cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                              cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn], 
                              Mass_here, Porg_here, Pmom_here, Spin_here, BH_NM);
          }
          if (BL == Pp->data->ble)
            break;
          BL = BL->next;
        }
        Pp = Pp->next;
      }
    }
    // dump read_in initial data
    for (int lev = 0; lev < GH->levels; lev++)
      Parallel::Dump_Data(GH->PatL[lev], StateList, 0, PhysTime, dT);
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
          f_get_initial_nbhs_sh(cg->shape, 
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
                                Mass_here, Porg_here, Pmom_here, Spin_here, BH_NM);
        }
        if (BL == Pp->data->ble)
          break;
        BL = BL->next;
      }
      Pp = Pp->next;
    }
    // dump read_in initial data
    SH->Dump_Data(StateList, 0, PhysTime, dT);
#endif

    delete[] Porg_here;
    delete[] Mass_here;
    delete[] Pmom_here;
    delete[] Spin_here;
    //   SH->Synch(GH->PatL[0],StateList,Symmetry);
    //   exit(0);
  }
}

//================================================================================================



//================================================================================================

// 该成员函数用曹老师的解析公式求解初值

//================================================================================================

void bssn_class::Setup_Initial_Data_Cao()
{
  if (!checkrun)
  {
    if (myrank == 0)
    {
      cout << endl;
      cout << " Setup initial data with Cao's analytical formula. " << endl;
      cout << endl;
    }
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
    double *Porg_here, *Pmom_here, *Spin_here, *Mass_here;
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
            ErrorMonitor->outfile << "error reading parameter file " << filename << " in line " << i << endl;
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
    Pmom_here = new double[3 * BH_NM];
    Spin_here = new double[3 * BH_NM];
    Mass_here = new double[BH_NM];
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
            ErrorMonitor->outfile << "error reading parameter file " << filename << " in line " << i << endl;
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else if (status == 0)
          continue;

        if (sgrp == "BSSN" && sind < BH_NM)
        {
          if (skey == "Mass")
            Mass_here[sind] = atof(sval.c_str());
          else if (skey == "Porgx")
            Porg_here[sind * 3] = atof(sval.c_str());
          else if (skey == "Porgy")
            Porg_here[sind * 3 + 1] = atof(sval.c_str());
          else if (skey == "Porgz")
            Porg_here[sind * 3 + 2] = atof(sval.c_str());
          else if (skey == "Spinx")
            Spin_here[sind * 3] = atof(sval.c_str());
          else if (skey == "Spiny")
            Spin_here[sind * 3 + 1] = atof(sval.c_str());
          else if (skey == "Spinz")
            Spin_here[sind * 3 + 2] = atof(sval.c_str());
          else if (skey == "Pmomx")
            Pmom_here[sind * 3] = atof(sval.c_str());
          else if (skey == "Pmomy")
            Pmom_here[sind * 3 + 1] = atof(sval.c_str());
          else if (skey == "Pmomz")
            Pmom_here[sind * 3 + 2] = atof(sval.c_str());
        }
      }
      inf.close();
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
            // 用曹老师的解析公式求解初值
            f_get_initial_nbhs(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                               cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                               cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                               cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                               cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                               cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                               cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                               cg->fgfs[Lap0->sgfn], 
                               cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                               cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn], 
                               Mass_here, Porg_here, Pmom_here, Spin_here, BH_NM);
          }
          if (BL == Pp->data->ble)
            break;
          BL = BL->next;
        }
        Pp = Pp->next;
      }
    }
    // dump read_in initial data
    for (int lev = 0; lev < GH->levels; lev++)
      Parallel::Dump_Data(GH->PatL[lev], StateList, 0, PhysTime, dT);
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
          f_get_initial_nbhs_sh(cg->shape, 
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
                                Mass_here, Porg_here, Pmom_here, Spin_here, BH_NM);
        }
        if (BL == Pp->data->ble)
          break;
        BL = BL->next;
      }
      Pp = Pp->next;
    }
    // dump read_in initial data
    SH->Dump_Data(StateList, 0, PhysTime, dT);
#endif

    delete[] Porg_here;
    delete[] Mass_here;
    delete[] Pmom_here;
    delete[] Spin_here;
    //   SH->Synch(GH->PatL[0],StateList,Symmetry);
    //   exit(0);
  }
}

//================================================================================================



//================================================================================================

// 该成员函数用解析方法求解的 KerrSchild 初值

//================================================================================================

void bssn_class::Setup_KerrSchild()
{
  if (!checkrun)
  {
    if (myrank == 0)
    {
      cout << endl;
      cout << " Setup initial data with Kerr-Schild formula. " << endl;
      cout << endl;
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
            f_get_initial_kerrschild(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                     cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                                     cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                                     cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                                     cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                                     cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                                     cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                                     cg->fgfs[Lap0->sgfn], 
                                     cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                                     cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn]);
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
      int lev = 0, fngfs = Pp->data->fngfs;

      MyList<Block> *BL = Pp->data->blb;
      while (BL)
      {
        Block *cg = BL->data;
        if (myrank == cg->rank)
        {
          f_get_initial_kerrschild_ss(cg->shape, 
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
                                      cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn]);
          /*
               f_fderivs_shc(cg->shape,cg->fgfs[phi0->sgfn],
                             cg->fgfs[Sfx_rhs->sgfn],
                             cg->fgfs[Sfy_rhs->sgfn],
                             cg->fgfs[Sfz_rhs->sgfn],
                             cg->X[0],cg->X[1],cg->X[2],
                             phi0->SoA[0],phi0->SoA[1],phi0->SoA[2],
                             Symmetry,lev,Pp->data->sst,
                             cg->fgfs[fngfs+ShellPatch::drhodx],
                             cg->fgfs[fngfs+ShellPatch::drhody],
                             cg->fgfs[fngfs+ShellPatch::drhodz],
                             cg->fgfs[fngfs+ShellPatch::dsigmadx],
                             cg->fgfs[fngfs+ShellPatch::dsigmady],
                             cg->fgfs[fngfs+ShellPatch::dsigmadz],
                             cg->fgfs[fngfs+ShellPatch::dRdx],
                             cg->fgfs[fngfs+ShellPatch::dRdy],
                             cg->fgfs[fngfs+ShellPatch::dRdz]);
               f_fdderivs_shc(cg->shape,cg->fgfs[phi0->sgfn],
                              cg->fgfs[Axx_rhs->sgfn],cg->fgfs[Axy_rhs->sgfn],cg->fgfs[Axz_rhs->sgfn],
                              cg->fgfs[Ayy_rhs->sgfn],cg->fgfs[Ayz_rhs->sgfn],cg->fgfs[Azz_rhs->sgfn],
                              cg->X[0],cg->X[1],cg->X[2],
                              phi0->SoA[0],phi0->SoA[1],phi0->SoA[2],
                              Symmetry,lev,Pp->data->sst,
                              cg->fgfs[fngfs+ShellPatch::drhodx],
                              cg->fgfs[fngfs+ShellPatch::drhody],
                              cg->fgfs[fngfs+ShellPatch::drhodz],
                              cg->fgfs[fngfs+ShellPatch::dsigmadx],
                              cg->fgfs[fngfs+ShellPatch::dsigmady],
                              cg->fgfs[fngfs+ShellPatch::dsigmadz],
                              cg->fgfs[fngfs+ShellPatch::dRdx],
                              cg->fgfs[fngfs+ShellPatch::dRdy],
                              cg->fgfs[fngfs+ShellPatch::dRdz],
                              cg->fgfs[fngfs+ShellPatch::drhodxx],
                              cg->fgfs[fngfs+ShellPatch::drhodxy],
                              cg->fgfs[fngfs+ShellPatch::drhodxz],
                              cg->fgfs[fngfs+ShellPatch::drhodyy],
                              cg->fgfs[fngfs+ShellPatch::drhodyz],
                              cg->fgfs[fngfs+ShellPatch::drhodzz],
                              cg->fgfs[fngfs+ShellPatch::dsigmadxx],
                              cg->fgfs[fngfs+ShellPatch::dsigmadxy],
                              cg->fgfs[fngfs+ShellPatch::dsigmadxz],
                              cg->fgfs[fngfs+ShellPatch::dsigmadyy],
                              cg->fgfs[fngfs+ShellPatch::dsigmadyz],
                              cg->fgfs[fngfs+ShellPatch::dsigmadzz],
                              cg->fgfs[fngfs+ShellPatch::dRdxx], 
                              cg->fgfs[fngfs+ShellPatch::dRdxy],
                              cg->fgfs[fngfs+ShellPatch::dRdxz],
                              cg->fgfs[fngfs+ShellPatch::dRdyy],
                              cg->fgfs[fngfs+ShellPatch::dRdyz],
                              cg->fgfs[fngfs+ShellPatch::dRdzz]);
          */
        }
        if (BL == Pp->data->ble)
          break;
        BL = BL->next;
      }
      Pp = Pp->next;
    }
#endif

    // dump read_in initial data
    //   SH->Synch(GH->PatL[0],StateList,Symmetry);
    //   for(int lev=0;lev<GH->levels;lev++) Parallel::Dump_Data(GH->PatL[lev],StateList,0,PhysTime,dT);
    //   SH->Dump_Data(StateList,0,PhysTime,dT);
    //   exit(0);

    /*
       {
             MyList<var> * DG_List=new MyList<var>(Sfx_rhs);
             DG_List->insert(Sfy_rhs); 
             DG_List->insert(Sfz_rhs);
             DG_List->insert(Axx_rhs); 
             DG_List->insert(Axy_rhs); 
             DG_List->insert(Axz_rhs);
             DG_List->insert(Ayy_rhs); 
             DG_List->insert(Ayz_rhs); 
             DG_List->insert(Azz_rhs);
             SH->Synch(DG_List,Symmetry);
             SH->Dump_Data(DG_List,0,PhysTime,dT);
             DG_List->clearList();
             exit(0);
       }
    */
  }
}

//================================================================================================



//================================================================================================

// 该成员函数读入 Pablo Galaviz 的 Olliptic 程序求出的初值

//================================================================================================

// Read initial data solved by Pablo's Olliptic Phys.Rev.D 82 024005 (2010)

//|----------------------------------------------------------------------------
//  read ASCII file with the style of Pablo
//|----------------------------------------------------------------------------
bool bssn_class::read_Pablo_file(int *ext, double *datain, char *filename)
{
  if (myrank == 0)
  {
    cout << endl;
    cout << " Setup initial data with Pablo_file. " << endl;
    cout << endl;
  }

  int nx = ext[0], ny = ext[1], nz = ext[2];
  int i, j, k;
  double x, y, z;
  //|--->open in put file
  ifstream infile;
  infile.open(filename);
  if (!infile)
  {
    cout << "bssn_class: read_Pablo_file can't open " << filename << " for input." << endl;
    return false;
  }
  for (k = 0; k < nz; k++)
    for (j = 0; j < ny; j++)
      for (i = 0; i < nx; i++)
      {
        infile >> x >> y >> z >> datain[i + j * nx + k * nx * ny];
      }

  infile.close();

  return true;
}

//================================================================================================



//================================================================================================

// 该成员函数写入 Pablo Galaviz 的 Olliptic 程序求出的初值

//================================================================================================

//|----------------------------------------------------------------------------
//  write ASCII file with the style of Pablo
//|----------------------------------------------------------------------------
void bssn_class::write_Pablo_file(int *ext, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax,
                                  char *filename)
{
  int nx = ext[0], ny = ext[1], nz = ext[2];
  int i, j, k;
  double *X, *Y, *Z;
  X = new double[nx];
  Y = new double[ny];
  Z = new double[nz];
  double dX, dY, dZ;
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
  dX = (xmax - xmin) / (nx - 1);
  for (i = 0; i < nx; i++)
    X[i] = xmin + i * dX;
  dY = (ymax - ymin) / (ny - 1);
  for (j = 0; j < ny; j++)
    Y[j] = ymin + j * dY;
  dZ = (zmax - zmin) / (nz - 1);
  for (k = 0; k < nz; k++)
    Z[k] = zmin + k * dZ;
#else
#ifdef Cell
  dX = (xmax - xmin) / nx;
  for (i = 0; i < nx; i++)
    X[i] = xmin + (i + 0.5) * dX;
  dY = (ymax - ymin) / ny;
  for (j = 0; j < ny; j++)
    Y[j] = ymin + (j + 0.5) * dY;
  dZ = (zmax - zmin) / nz;
  for (k = 0; k < nz; k++)
    Z[k] = zmin + (k + 0.5) * dZ;
#else
#error Not define Vertex nor Cell
#endif
#endif
  //|--->open out put file
  ofstream outfile;
  outfile.open(filename);
  if (!outfile)
  {
    cout << "bssn_class: write_Pablo_file can't open " << filename << " for output." << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  outfile.setf(ios::scientific, ios::floatfield);
  outfile.precision(16);
  for (k = 0; k < nz; k++)
    for (j = 0; j < ny; j++)
      for (i = 0; i < nx; i++)
      {
        outfile << X[i] << " " << Y[j] << " " << Z[k] << " "
                << 0 << endl;
      }
  outfile.close();

  delete[] X;
  delete[] Y;
  delete[] Z;
}

//================================================================================================



//================================================================================================

// 该成员函数读入 Ansorg 方法求解的 TwoPuncture 初值

//================================================================================================

// Read initial data solved by Ansorg, PRD 70, 064011 (2004)

void bssn_class::Read_Ansorg()
{
  if (!checkrun)
  {
    if (myrank == 0)
    {
      cout << endl;
      cout << " Read initial data from Ansorg's solver,"
           << " please be sure the input parameters for black holes are puncture parameters!! " << endl;
      cout << endl;
    }
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
    double *Porg_here, *Pmom_here, *Spin_here, *Mass_here;
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
            ErrorMonitor->outfile << "error reading parameter file " << filename << " in line " << i << endl;
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
    Pmom_here = new double[3 * BH_NM];
    Spin_here = new double[3 * BH_NM];
    Mass_here = new double[BH_NM];
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
            ErrorMonitor->outfile << "error reading parameter file " << filename << " in line " << i << endl;
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else if (status == 0)
          continue;

        if (sgrp == "BSSN" && sind < BH_NM)
        {
          if (skey == "Mass")
            Mass_here[sind] = atof(sval.c_str());
          else if (skey == "Porgx")
            Porg_here[sind * 3] = atof(sval.c_str());
          else if (skey == "Porgy")
            Porg_here[sind * 3 + 1] = atof(sval.c_str());
          else if (skey == "Porgz")
            Porg_here[sind * 3 + 2] = atof(sval.c_str());
          else if (skey == "Spinx")
            Spin_here[sind * 3] = atof(sval.c_str());
          else if (skey == "Spiny")
            Spin_here[sind * 3 + 1] = atof(sval.c_str());
          else if (skey == "Spinz")
            Spin_here[sind * 3 + 2] = atof(sval.c_str());
          else if (skey == "Pmomx")
            Pmom_here[sind * 3] = atof(sval.c_str());
          else if (skey == "Pmomy")
            Pmom_here[sind * 3 + 1] = atof(sval.c_str());
          else if (skey == "Pmomz")
            Pmom_here[sind * 3 + 2] = atof(sval.c_str());
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

            f_get_ansorg_nbhs(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                              cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                              cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                              cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                              cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                              cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                              cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                              cg->fgfs[Lap0->sgfn], 
                              cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                              cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn],
                              Mass_here, Porg_here, Pmom_here, Spin_here, BH_NM);
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

          f_get_ansorg_nbhs_ss(cg->shape, 
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
                               Mass_here, Porg_here, Pmom_here, Spin_here, BH_NM);
#if 0
// for check fderivs_sh
            f_fderivs_sh(cg->shape,cg->fgfs[Ayz0->sgfn],
			 cg->fgfs[Sfx0->sgfn],cg->fgfs[Sfy0->sgfn],cg->fgfs[Sfz0->sgfn],
                         cg->X[0],cg->X[1],cg->X[2],
                         Ayz0->SoA[0],Ayz0->SoA[1],Ayz0->SoA[2],
                         Symmetry,Pp->data->sst,Pp->data->sst);
#endif
#if 0
// for check fderivs_shc
            int fngfs = Pp->data->fngfs;
            f_fderivs_shc(cg->shape,cg->fgfs[Ayz0->sgfn],
			 cg->fgfs[Sfx0->sgfn],cg->fgfs[Sfy0->sgfn],cg->fgfs[Sfz0->sgfn],
                         cg->X[0],cg->X[1],cg->X[2],
                         Ayz0->SoA[0],Ayz0->SoA[1],Ayz0->SoA[2],
                         Symmetry,Pp->data->sst,Pp->data->sst,
                         cg->fgfs[fngfs+ShellPatch::drhodx],
                         cg->fgfs[fngfs+ShellPatch::drhody],
                         cg->fgfs[fngfs+ShellPatch::drhodz],
                         cg->fgfs[fngfs+ShellPatch::dsigmadx],
                         cg->fgfs[fngfs+ShellPatch::dsigmady],
                         cg->fgfs[fngfs+ShellPatch::dsigmadz],
                         cg->fgfs[fngfs+ShellPatch::dRdx],
                         cg->fgfs[fngfs+ShellPatch::dRdy],
                         cg->fgfs[fngfs+ShellPatch::dRdz]);
#endif
        }
        if (BL == Pp->data->ble)
          break;
        BL = BL->next;
      }
      Pp = Pp->next;
    }
#endif

    delete[] Porg_here;
    delete[] Mass_here;
    delete[] Pmom_here;
    delete[] Spin_here;

    Compute_Constraint();
    // dump read_in initial data
    for (int lev = 0; lev < GH->levels; lev++)
      Parallel::Dump_Data(GH->PatL[lev], DumpList, 0, PhysTime, dT);
#ifdef WithShell
    SH->Dump_Data(DumpList, 0, PhysTime, dT);
#endif
    //   if(myrank==0) MPI_Abort(MPI_COMM_WORLD,1);
  }
}

//================================================================================================



//================================================================================================

// 该成员函数设定了整个过程中的时间演化

//================================================================================================

void bssn_class::Evolve(int Steps)
{
  clock_t prev_clock, curr_clock;
  double LastDump = 0.0, LastCheck = 0.0, Last2dDump = 0.0;
  LastAnas = 0;
#if 0
//initial checkpoint for special uasge
     {
       CheckPoint->write_Black_Hole_position(BH_num_input,BH_num,Porg0,Porgbr,Mass);
       CheckPoint->writecheck_cgh(PhysTime,GH);
#ifdef WithShell   
       CheckPoint->writecheck_sh(PhysTime,SH);
#endif
       CheckPoint->write_bssn(LastDump,Last2dDump,LastAnas);
       misc::tillherecheck("complete initialization preparation"); // we need synchronization here
       if(myrank==0) MPI_Abort(MPI_COMM_WORLD,1);
     }
#endif
  // for step 0 constraint interpolation
  Interp_Constraint(true);

#ifdef With_AHF
  // setup apparent horizon finder direct of thornburg
  {
    HN_num = BH_num;
    for (int ia = 0; ia < BH_num; ia++)
      for (int ib = ia + 1; ib < BH_num; ib++)
        HN_num++;

    AHFinderDirect::AHFinderDirect_setup(AHList, GaugeList,
                                         this,
                                         Symmetry, HN_num, &PhysTime);

    lastahdumpid = new int[HN_num];
    findeveryl = new int[HN_num];
    xc = new double[HN_num];
    yc = new double[HN_num];
    zc = new double[HN_num];
    xr = new double[HN_num];
    yr = new double[HN_num];
    zr = new double[HN_num];
    dTT = new double[HN_num];
    trigger = new bool[HN_num];
    dumpid = new int[HN_num];

    for (int ihn = 0; ihn < HN_num; ihn++)
    {
      lastahdumpid[ihn] = 0;
      findeveryl[ihn] = AHfindevery;
    }
  }
#endif

  if (checkrun)
    CheckPoint->read_bssn(LastDump, Last2dDump, LastAnas);

  double dT_mon = dT * pow(0.5, Mymax(0, trfls));
  
  /*
  #ifdef With_AHF
  //initial apparent horizon finding
      {
         double gam;
         double massmin=Mass[0];
         for(int ihn=1;ihn<BH_num;ihn++) massmin=Mymin(massmin,Mass[ihn]);

         for(int ihn=0;ihn<BH_num;ihn++)
         {
           xc[ihn] = Porg0[ihn][0];
           yc[ihn] = Porg0[ihn][1];
           zc[ihn] = Porg0[ihn][2];
           gam = fabs(Pmom[ihn*3])/(Mass[ihn]);
           gam = sqrt(1-gam*gam);
           xr[ihn] = Mass[ihn]*gam;
           gam = fabs(Pmom[ihn*3+1])/(Mass[ihn]);
           gam = sqrt(1-gam*gam);
           yr[ihn] = Mass[ihn]*gam;
           gam = fabs(Pmom[ihn*3+2])/(Mass[ihn]);
           gam = sqrt(1-gam*gam);
           zr[ihn] = Mass[ihn]*gam;
           findeveryl[ihn] = findeveryl[ihn]*(Mymax(int(Mass[ihn]/massmin),1));
         }
         int ihn = BH_num;
         for(int ia=0;ia<BH_num;ia++)
           for(int ib=ia+1;ib<BH_num;ib++) {
             xc[ihn] = (Porg0[ia][0] + Porg0[ib][0])/2;
             yc[ihn] = (Porg0[ia][1] + Porg0[ib][1])/2;
             zc[ihn] = (Porg0[ia][2] + Porg0[ib][2])/2;
             xr[ihn] = yr[ihn] = zr[ihn] = (Mass[ia])+(Mass[ib]);
             findeveryl[ihn] = findeveryl[ihn]*(int(xr[ihn]/Mass[0]));
             ihn++;
           }

         AHFinderDirect::AHFinderDirect_enforcefind(HN_num,xc,yc,zc,xr,yr,zr);  
         // note rhs has been used as temp storage space

         delete[] xc; delete[] yc; delete[] zc; delete[] xr; delete[] yr; delete[] zr;
      }
  #endif
  */
  
  perf bssn_perf;
  size_t current_min, current_avg, current_max, peak_min, peak_avg, peak_max;

  for (int lev = 0; lev < GH->levels; lev++)
    GH->Lt[lev] = PhysTime;

  GH->settrfls(trfls);

  for (int ncount = 1; ncount < Steps + 1; ncount++)
  {
    // special for large mass ratio consideration
    //     if(fabs(Porg0[0][0]-Porg0[1][0])+fabs(Porg0[0][1]-Porg0[1][1])+fabs(Porg0[0][2]-Porg0[1][2])<1e-6) 
    //     { GH->levels=GH->movls; }

    if (myrank == 0)
      curr_clock = clock();
#if (PSTR == 0)
    RecursiveStep(0);
#elif (PSTR == 1 || PSTR == 2 || PSTR == 3)
    // data analysis part
    // Warning NOTE: the variables1 are used as temp storege room
    AnalysisStuff(a_lev, dT_mon);
    ParallelStep();
#endif

    // misc::tillherecheck("before Constraint_Out");

    Constraint_Out(); // this will affect the Dump_List

    LastDump += dT_mon;
    Last2dDump += dT_mon;
    LastCheck += dT_mon;

    // 每当超出 DumpTime 时间，输出相应的 2 进制数据
    if (LastDump >= DumpTime)
    {
      //       misc::tillherecheck("before Dump_Data");

      for (int lev = 0; lev < GH->levels; lev++)
        Parallel::Dump_Data(GH->PatL[lev], DumpList, 0, PhysTime, dT_mon);
#ifdef WithShell
      SH->Dump_Data(DumpList, 0, PhysTime, dT_mon);
#endif

      LastDump = 0;

      if (myrank == 0)
      {
        cout << " Dump done. " << endl;
      }
    }

    // 每当超出 d2DumpTime 时间，输出相应的 2 维数据
    if (Last2dDump >= d2DumpTime)
    {
      //       misc::tillherecheck("before 2dDump_Data");

      for (int lev = 0; lev < GH->levels; lev++)
        Parallel::d2Dump_Data(GH->PatL[lev], DumpList, 0, PhysTime, dT_mon);

      Last2dDump = 0;

      if (myrank == 0)
      {
        cout << " 2d Dump done. " << endl;
      }
    }

    if (myrank == 0)
    {
      prev_clock = curr_clock;
      curr_clock = clock();
      cout << endl;
      cout << " Timestep # " << ncount << ": integrating to time: " << PhysTime << "   "
           << " Computer used " << (double)(curr_clock - prev_clock) / ((double)CLOCKS_PER_SEC) 
           << " seconds! " << endl;
      // cout << endl;
    }

    if (PhysTime >= TotalTime)
      break;

#if (REGLEV == 1)
    GH->Regrid(Symmetry, BH_num, Porgbr, Porg0,
               SynchList_cor, OldStateList, StateList, SynchList_pre,
               fgt(PhysTime - dT_mon, StartTime, dT_mon / 2), ErrorMonitor);
#endif

#if (REGLEV == 0 && (PSTR == 1 || PSTR == 2))
//     GH->Regrid_fake(Symmetry,BH_num,Porgbr,Porg0,
//		SynchList_cor,OldStateList,StateList,SynchList_pre,
//		fgt(PhysTime-dT_mon,StartTime,dT_mon/2),ErrorMonitor);
#endif

    // 获取计算中用到的内存信息，主进程打印到屏幕上
    bssn_perf.MemoryUsage(&current_min, &current_avg, &current_max,
                          &peak_min, &peak_avg, &peak_max, nprocs);
    if (myrank == 0)
    {
      printf(" Memory usage: current %0.4lg/%0.4lg/%0.4lgMB, "
             "peak %0.4lg/%0.4lg/%0.4lgMB\n",
             (double)current_min / (1024.0 * 1024.0),
             (double)current_avg / (1024.0 * 1024.0),
             (double)current_max / (1024.0 * 1024.0),
             (double)peak_min / (1024.0 * 1024.0),
             (double)peak_avg / (1024.0 * 1024.0),
             (double)peak_max / (1024.0 * 1024.0));
      cout << endl;
    }
    
    // 输出每一步的 puncture 位置
    if (myrank == 0)
    {
      for (int i_count=0; i_count<BH_num; i_count++)
      {
        cout << " puncture position: no." 
             << setw(2) << setfill(' ')      << i_count                 
             << " = ("  << Porg0[i_count][0] << " " 
                        << Porg0[i_count][1] << " "
                        << Porg0[i_count][2] << ")" 
             << endl;
      }
      cout << endl;
      cout << " If you think the physical evolution time is enough for this simulation, please input 'stop' in the terminal to stop the MPI processes in the next evolution step ! " << endl;
      // cout << endl;
    }
    
    ////////////////////////////////////////////////////////////
    // 如果检测到有输入 "abort"，则退出 MPI 程序
    ////////////////////////////////////////////////////////////
    
    bool shouldAbort = false;
    
    // 只有进程 0 检查 stdin
    if (myrank == 0) {
        if (check_Stdin_Abort()) {
            shouldAbort = true;
        }
    }

    // 广播 abort 信号给所有进程
    int abortFlag = shouldAbort ? 1 : 0;
    MPI_Bcast(&abortFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (abortFlag && myrank == 0) {
        cout << endl;
        cout << " Aborting the MPI Process due to 'abort' command !! " << endl;
        cout << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);  // 退出 MPI 进程
        // MPI_Finalize();
    }
        
    ////////////////////////////////////////////////////////////

    // 每当超过 CheckTime，检查程序运行情况，并输出该时刻相应的数据
    if (LastCheck >= CheckTime)
    {
      LastCheck = 0;

      CheckPoint->write_Black_Hole_position(BH_num_input, BH_num, Porg0, Porgbr, Mass);
      CheckPoint->writecheck_cgh(PhysTime, GH);
#ifdef WithShell
      CheckPoint->writecheck_sh(PhysTime, SH);
#endif
      CheckPoint->write_bssn(LastDump, Last2dDump, LastAnas);
    }
  }
  /*
  #ifdef With_AHF
  // final apparent horizon finding
      {
         double gam;
         for(int ihn=0;ihn<BH_num;ihn++)
         {
           xc[ihn] = Porg0[ihn][0];
           yc[ihn] = Porg0[ihn][1];
           zc[ihn] = Porg0[ihn][2];
           gam = fabs(Pmom[ihn*3]  )/(Mass[ihn]);
           gam = sqrt(1-gam*gam);
           xr[ihn] = Mass[ihn]*gam;
           gam = fabs(Pmom[ihn*3+1])/(Mass[ihn]);
           gam = sqrt(1-gam*gam);
           yr[ihn] = Mass[ihn]*gam;
           gam = fabs(Pmom[ihn*3+2])/(Mass[ihn]);
           gam = sqrt(1-gam*gam);
           zr[ihn] = Mass[ihn]*gam;
         }
         int ihn = BH_num;
         for(int ia=0;ia<BH_num;ia++)
           for(int ib=ia+1;ib<BH_num;ib++)
           {
             xc[ihn] = (Porg0[ia][0] + Porg0[ib][0])/2;
             yc[ihn] = (Porg0[ia][1] + Porg0[ib][1])/2;
             zc[ihn] = (Porg0[ia][2] + Porg0[ib][2])/2;
             xr[ihn] = yr[ihn] = zr[ihn] = (Mass[ia])+(Mass[ib]);

             ihn++;
           }

         AHFinderDirect::AHFinderDirect_enforcefind(HN_num,xc,yc,zc,xr,yr,zr);  
         // note rhs has been used as temp storage space

         delete[] xc; delete[] yc; delete[] zc; delete[] xr; delete[] yr; delete[] zr;
      }
  #endif
  */
}

//================================================================================================




//================================================================================================

// 该成员函数设定了时间演化过程中的不同层网格间递归时间演化
// 仅仅针对 PSTR == 0 的情况

//================================================================================================

#if (PSTR == 0)
void bssn_class::RecursiveStep(int lev)
{
  double dT_lev = dT * pow(0.5, Mymax(lev, trfls));

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

#if (AGM == 2)
    if (GH->levels == 1)
    {
      Enforce_algcon(lev, 0);
    }
#endif

    GH->Lt[lev] += dT_lev;

    if (lev < GH->levels - 1)
    {
      int lf = lev + 1;
      RecursiveStep(lf);
    }
    else
      PhysTime += dT * pow(0.5, lev);

#if (AGM == 2)
    if (lev > 0)
    {
      Enforce_algcon(lev, 0);
      if (YN == 1)
        Enforce_algcon(lev - 1, 0);
    }
#endif

#if (RPS == 1)
    // mesh refinement boundary part
    //
    // till here the PhysTime has updated dT_lev
    //  if(myrank==0) cout<<"level now = "<<lev<<", "<<fgt(PhysTime-dT_lev,StartTime,dT_lev/2)<<endl;
    RestrictProlong(lev, YN, fgt(PhysTime - dT_lev, StartTime, dT_lev / 2), StateList, OldStateList, SynchList_cor);
    // RestrictProlong(lev,YN,false,StateList,OldStateList,SynchList_cor);

#ifdef WithShell
    if (lev == 0)
    {
      clock_t prev_clock, curr_clock;
      if (myrank == 0)
        curr_clock = clock();
      SH->CS_Inter(StateList, Symmetry);
      if (myrank == 0)
      {
        prev_clock = curr_clock;
        curr_clock = clock();
        cout << " CS_Inter used " << (double)(curr_clock - prev_clock) / ((double)CLOCKS_PER_SEC) << " seconds! " << endl;
      }
    }
#endif

#endif

    //     Parallel::Dump_Data(GH->PatL[lev],StateList,0,PhysTime,dT_lev);
  }

#if 0
    if(lev>0) Parallel::Restrict_after(GH->PatL[lev-1],GH->PatL[lev],StateList,StateList,Symmetry);
#endif

#if (REGLEV == 0)
  GH->Regrid_Onelevel(lev, Symmetry, BH_num, Porgbr, Porg0,
                      SynchList_cor, OldStateList, StateList, SynchList_pre,
                      fgt(PhysTime - dT_lev, StartTime, dT_lev / 2), ErrorMonitor);
#endif
}

//================================================================================================



//================================================================================================

// 该成员函数设定了时间演化过程中的不同层网格间递归时间演化
// 针对 PSTR == 1 和 PSTR == 2 的情况

//================================================================================================

#elif (PSTR == 1 || PSTR == 2)
void bssn_class::RecursiveStep(int lev)
{
  double dT_lev = dT * pow(0.5, Mymax(lev, trfls));

  int NoIterations = 1, YN;
  if (lev <= trfls)
    NoIterations = 1;
  else
    NoIterations = 2;

  for (int i = 0; i < NoIterations; i++)
  {
    //     if(myrank==0) cout<<"level now = "<<lev<<" NoIteration = "<<i<<endl;
    YN = (i == NoIterations - 1) ? 1 : 0; // 1: same time level for coarse level and fine level

    if (GH->mylev == lev)
    {
      Step(lev, YN);

#if (AGM == 2)
      if (GH->levels == 1)
      {
        Enforce_algcon(lev, 0);
      }
      else if (lev > 0)
      {
        Enforce_algcon(lev, 0);
        if (YN == 1)
          Enforce_algcon(lev - 1, 0);
      }
#endif
    }

    if (GH->mylev == lev || GH->mylev == lev - 1)
    {
#if (RPS == 1)
      // mesh refinement boundary part
      //
      // till here the PhysTime has updated dT_lev
      RestrictProlong(lev, YN, fgt(PhysTime - dT_lev, StartTime, dT_lev / 2), StateList, OldStateList, SynchList_cor);
#endif
    }

    GH->Lt[lev] += dT_lev;

    PhysTime += dT * pow(0.5, lev);
  }

  if (lev > 0)
  {
    int lc = lev - 1;
    RecursiveStep(lc);
  }
}

//================================================================================================



//================================================================================================

// ParallelStep 是对不同网格层进行时间演化（包含并行计算）

//================================================================================================

#if 1
void bssn_class::ParallelStep()
{
  int YN, NoIterations = 1;
  if (GH->mylev <= trfls)
    NoIterations = 1;
  else
    NoIterations = int(pow(2.0, GH->mylev - trfls));

  double dT_lev = dT * pow(0.5, Mymax(GH->mylev, trfls));

  for (int i = 0; i < NoIterations; i++)
  {
    YN = (i % 2) ? 1 : 0; // 1: same time level for coarse level and fine level

    Step(GH->mylev, YN);
#if (RPS == 1)
    // mesh refinement boundary part
    //
    // till here the PhysTime has updated dT_lev
    if (GH->mylev < GH->levels - 1)
    {
      RestrictProlong(GH->mylev + 1, 0, fgt(PhysTime - dT_lev / 2, StartTime, dT_lev / 4), StateList, OldStateList, SynchList_cor);
      RestrictProlong(GH->mylev + 1, 1, fgt(PhysTime - dT_lev / 2, StartTime, dT_lev / 4), StateList, OldStateList, SynchList_cor);
    }
    if (GH->mylev > 0)
      RestrictProlong(GH->mylev, YN, fgt(PhysTime - dT_lev, StartTime, dT_lev / 2), StateList, OldStateList, SynchList_cor);
#endif

#ifdef WithShell
    SHStep();
#if (RPS == 1)
    {
      clock_t prev_clock, curr_clock;
      if (myrank == 0)
        curr_clock = clock();
      SH->CS_Inter(StateList, Symmetry);
      if (myrank == 0)
      {
        prev_clock = curr_clock;
        curr_clock = clock();
        cout << " CS_Inter used " << (double)(curr_clock - prev_clock) / ((double)CLOCKS_PER_SEC) << " seconds! " << endl;
      }
    }
#endif
#endif
    PhysTime += dT_lev;
  }

  //  stringstream a_stream;
  //  a_stream.setf(ios::left);

  double *tporg, *tporgo;
  tporg = new double[3 * BH_num];
  tporgo = new double[3 * BH_num];
  {
    int lev = GH->mylev;
    MPI_Status status;
    // receive
    if (lev < GH->levels - 1)
    {
      if (myrank == GH->start_rank[lev])
      {
        MPI_Recv(tporgo, 3 * BH_num, MPI_DOUBLE, GH->start_rank[lev + 1], 1, MPI_COMM_WORLD, &status);
        //	 cout<<tporgo[0]<<","<<tporgo[1]<<","<<tporgo[2]<<endl;
      }
      else
        for (int i = 0; i < 3 * BH_num; i++)
          tporgo[i] = 0;
      MPI_Allreduce(tporgo, tporg, 3 * BH_num, MPI_DOUBLE, MPI_SUM, GH->Commlev[lev]);

      for (int i = 0; i < BH_num; i++)
        for (int j = 0; j < 3; j++)
          Porg0[i][j] = tporg[3 * i + j];

      //      if(myrank==GH->start_rank[lev]) cout<<Porg0[0][0]<<","<<Porg0[0][1]<<","<<Porg0[0][2]<<endl;
    }
    // send
    if (lev > 0 && myrank == GH->start_rank[lev])
    {
      for (int i = 0; i < BH_num; i++)
        for (int j = 0; j < 3; j++)
          tporg[3 * i + j] = Porg0[i][j];

      MPI_Send(tporg, 3 * BH_num, MPI_DOUBLE, GH->start_rank[lev - 1], 1, MPI_COMM_WORLD);
    }

    //   a_stream.clear();
    //   a_stream.str("");
    //   a_stream<<lev<<": after PBH Sync";
    //   misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],a_stream.str());
  }
  delete[] tporg;
  delete[] tporgo;
#if (REGLEV == 0)
  GH->Regrid_Onelevel(GH->mylev, Symmetry, BH_num, Porgbr, Porg0,
                      SynchList_cor, OldStateList, StateList, SynchList_pre,
                      fgt(PhysTime - dT_lev, StartTime, dT_lev / 2), ErrorMonitor);
#endif
}

//================================================================================================



//================================================================================================

// ParallelStep 是对不同网格层进行时间演化（包含并行计算）
// 这是另一份代码

//================================================================================================

#else
void bssn_class::ParallelStep()
{
  //  stringstream a_stream;
  //  a_stream.setf(ios::left);

  double *tporg, *tporgo;
  tporg = new double[3 * BH_num];
  tporgo = new double[3 * BH_num];

  int lev = GH->mylev;
  double dT_lev = dT * pow(0.5, Mymax(lev, trfls));
  double dT_levp1 = dT * pow(0.5, Mymax(lev + 1, trfls));
  double dT_levm1 = dT * pow(0.5, Mymax(lev - 1, trfls));

  int NoIterations = 1, YN;
  if (lev <= trfls)
    NoIterations = 1;
  else
    NoIterations = int(pow(2.0, lev - trfls));

  for (int i = 0; i < NoIterations; i++)
  {
    //     if(myrank==GH->start_rank[lev]) cout<<"level now = "<<lev<<" NoIteration = "<<i<<endl;

    if (NoIterations == 1)
      YN = 1;
    else if (i % 2 == 0)
      YN = 0;
    else
      YN = 1; // 1: same time level for coarse level and fine level

    //     a_stream<<lev<<": before calling Step";
    //     misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],a_stream.str());
    Step(lev, YN);

    //     a_stream.clear();
    //     a_stream.str("");
    //     a_stream<<lev<<": after calling Step";
    //     misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],a_stream.str());

#if (AGM == 2)
    if (GH->levels == 1)
    {
      Enforce_algcon(lev, 0);
    }
#endif

    GH->Lt[lev] += dT_lev;

    PhysTime += dT_lev;

#if (AGM == 2)
    if (lev > 0)
    {
      Enforce_algcon(lev, 0);
      if (YN == 1)
        Enforce_algcon(lev - 1, 0);
    }
#endif

#if (RPS == 1)
    // mesh refinement boundary part
    //
    // till here the PhysTime has updated dT_lev
    //   a_stream.clear();
    //   a_stream.str("");
    //   a_stream<<lev<<": before RestrictProlong";
    //   misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],a_stream.str());
    if (lev < GH->levels - 1)
    {
      if (lev + 1 <= trfls)
      {
        //	RestrictProlong_aux(lev,1,fgt(PhysTime-dT_lev,StartTime,dT_levp1/2),StateList,OldStateList,SynchList_cor);
        RestrictProlong(lev + 1, 1, fgt(PhysTime - dT_lev, StartTime, dT_levp1 / 2), StateList, OldStateList, SynchList_cor);
      }
      else
      {
        //        if(myrank==GH->start_rank[lev]) cout<<GH->mylev<<", "<<YN<<endl;
        //        misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],"between RestrictProlong");

        //	RestrictProlong_aux(lev,0,fgt(PhysTime-dT_lev,StartTime,dT_levp1/2),StateList,OldStateList,SynchList_cor);
        //	RestrictProlong_aux(lev,1,fgt(PhysTime-dT_levp1,StartTime,dT_levp1/2),StateList,OldStateList,SynchList_cor);
        RestrictProlong(lev + 1, 0, fgt(PhysTime - dT_lev, StartTime, dT_levp1 / 2), StateList, OldStateList, SynchList_cor);
        RestrictProlong(lev + 1, 1, fgt(PhysTime - dT_levp1, StartTime, dT_levp1 / 2), StateList, OldStateList, SynchList_cor);
      }
    }

    //   if(myrank==GH->start_rank[lev]) cout<<GH->mylev<<", "<<YN<<endl;
    //   a_stream.clear();
    //   a_stream.str("");
    //   a_stream<<lev<<": between RestrictProlong";
    //   misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],a_stream.str());

    RestrictProlong(lev, YN, fgt(PhysTime - dT_lev, StartTime, dT_lev / 2), StateList, OldStateList, SynchList_cor);
    // RestrictProlong(lev,YN,false,StateList,OldStateList,SynchList_cor);

//   if(myrank==GH->start_rank[lev]) cout<<GH->mylev<<endl;
//   a_stream.clear();
//   a_stream.str("");
//   a_stream<<lev<<": after RestrictProlong";
//   misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],a_stream.str());
#endif

    //     Parallel::Dump_Data(GH->PatL[lev],StateList,0,PhysTime,dT_lev);

    {
      MPI_Status status;
      // receive
      if (lev < GH->levels - 1)
      {
        if (myrank == GH->start_rank[lev])
        {
          MPI_Recv(tporgo, 3 * BH_num, MPI_DOUBLE, GH->start_rank[lev + 1], 1, MPI_COMM_WORLD, &status);
          //	 cout<<tporgo[0]<<","<<tporgo[1]<<","<<tporgo[2]<<endl;
        }
        else
          for (int i = 0; i < 3 * BH_num; i++)
            tporgo[i] = 0;
        MPI_Allreduce(tporgo, tporg, 3 * BH_num, MPI_DOUBLE, MPI_SUM, GH->Commlev[lev]);

        for (int i = 0; i < BH_num; i++)
          for (int j = 0; j < 3; j++)
            Porg0[i][j] = tporg[3 * i + j];

        //      if(myrank==GH->start_rank[lev]) cout<<Porg0[0][0]<<","<<Porg0[0][1]<<","<<Porg0[0][2]<<endl;
      }
      // send
      if (lev > 0 && YN == 1 && myrank == GH->start_rank[lev])
      {
        for (int i = 0; i < BH_num; i++)
          for (int j = 0; j < 3; j++)
            tporg[3 * i + j] = Porg0[i][j];

        MPI_Send(tporg, 3 * BH_num, MPI_DOUBLE, GH->start_rank[lev - 1], 1, MPI_COMM_WORLD);
      }

      //   a_stream.clear();
      //   a_stream.str("");
      //   a_stream<<lev<<": after PBH Sync";
      //   misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],a_stream.str());
    }
#if (REGLEV == 0)
    // for higher level
    if (lev < GH->levels - 1)
    {
      if (lev + 1 >= GH->movls)
      {
        //	       GH->Regrid_Onelevel_aux(lev,Symmetry,BH_num,Porgbr,Porg0,
        GH->Regrid_Onelevel(lev + 1, Symmetry, BH_num, Porgbr, Porg0,
                            SynchList_cor, OldStateList, StateList, SynchList_pre,
                            fgt(PhysTime - dT_levp1, StartTime, dT_levp1 / 2), ErrorMonitor);

        //               a_stream.clear();
        //               a_stream.str("");
        //               a_stream<<lev<<": after calling GH->Regrid_Onelevel_aux for higher level";
        //               misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],a_stream.str());
      }
    }

    // for this level
    if (YN == 1)
    {
      GH->Regrid_Onelevel(lev, Symmetry, BH_num, Porgbr, Porg0,
                          SynchList_cor, OldStateList, StateList, SynchList_pre,
                          fgt(PhysTime - dT_lev, StartTime, dT_lev / 2), ErrorMonitor);

      //               a_stream.clear();
      //               a_stream.str("");
      //               a_stream<<lev<<": after calling GH->Regrid_Onelevel";
      //               misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],a_stream.str());
    }

    // for lower level
    if (lev - 1 >= GH->movls)
    {
      if (lev - 1 <= trfls)
      {
        if (YN == 1)
        {
          //	   GH->Regrid_Onelevel_aux(lev-2,Symmetry,BH_num,Porgbr,Porg0,
          GH->Regrid_Onelevel(lev - 1, Symmetry, BH_num, Porgbr, Porg0,
                              SynchList_cor, OldStateList, StateList, SynchList_pre,
                              fgt(PhysTime - dT_lev, StartTime, dT_levm1 / 2), ErrorMonitor);

          //               a_stream.clear();
          //               a_stream.str("");
          //               a_stream<<lev<<": after calling GH->Regrid_Onelevel_aux for lower level";
          //               misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],a_stream.str());
        }
      }
      else
      {
        if (i % 4 == 3)
        {
          //	   GH->Regrid_Onelevel_aux(lev-2,Symmetry,BH_num,Porgbr,Porg0,
          GH->Regrid_Onelevel(lev - 1, Symmetry, BH_num, Porgbr, Porg0,
                              SynchList_cor, OldStateList, StateList, SynchList_pre,
                              fgt(PhysTime - dT_lev, StartTime, dT_levm1 / 2), ErrorMonitor);

          //               a_stream.clear();
          //               a_stream.str("");
          //               a_stream<<lev<<": after calling GH->Regrid_Onelevel_aux for lower level";
          //               misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],a_stream.str());
        }
      }
    }
#endif
  }

#ifdef WithShell
  SHStep();
  //               a_stream.clear();
  //               a_stream.str("");
  //               a_stream<<lev<<": after calling SHStep";
  //               misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],a_stream.str());

#if (RPS == 1)
  {
    clock_t prev_clock, curr_clock;
    if (myrank == 0)
      curr_clock = clock();
    SH->CS_Inter(StateList, Symmetry);
    if (myrank == 0)
    {
      prev_clock = curr_clock;
      curr_clock = clock();
      cout << " CS_Inter used " << (double)(curr_clock - prev_clock) / ((double)CLOCKS_PER_SEC) 
           << " seconds! " << endl;
    }
    //               a_stream.clear();
    //               a_stream.str("");
    //               a_stream<<lev<<": after calling shell cartisian sync";
    //               misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],a_stream.str());
  }
#endif

#endif

#if 0
    if(lev>0) Parallel::Restrict_after(GH->PatL[lev-1],GH->PatL[lev],StateList,StateList,Symmetry);
#endif

  delete[] tporg;
  delete[] tporgo;
}
#endif

//================================================================================================



//================================================================================================

// ParallelStep 是对不同网格层进行时间演化（包含并行计算）
// 这是又一份代码，针对 PSTR == 3 情形

//================================================================================================

#elif (PSTR == 3)
#warning "remember do not use Shell"
void bssn_class::ParallelStep()
{
  //  stringstream a_stream;
  //  a_stream.setf(ios::left);

  double *tporg, *tporgo;
  tporg = new double[3 * BH_num];
  tporgo = new double[3 * BH_num];

  int lev = GH->mylev;
  double dT_lev = dT * pow(0.5, Mymax(GH->levels - 1, trfls));
  if (lev == 1)
  {
    lev = GH->levels - 1;
    for (int i = 0; i < misc::MYpow2(lev); i++)
    {
      Step(lev, i % 2);
      PhysTime += dT_lev;
      //        if(myrank==nprocs-1) cout<<"OOO level now = "<<lev<<", "<<i%2<<", "<<fgt(PhysTime-dT_lev,StartTime,dT_lev/2)<<endl;
      RestrictProlong(lev, i % 2, fgt(PhysTime - dT_lev, StartTime, dT_lev / 2), StateList, OldStateList, SynchList_cor);
    }
  }
  else
  {
    lev = GH->levels - 2;
    for (int i = 1; i < misc::MYpow2(lev + 1); i++)
    {
      RecursiveStep(lev, i);
      PhysTime += dT_lev;
      if (i % 2 == 0)
      {
        //          if(myrank==0) cout<<"level now = "<<lev+1<<", "<<fgt(PhysTime-dT_lev,StartTime,dT_lev/2)<<endl;
        RestrictProlong(lev + 1, 1, fgt(PhysTime - dT_lev, StartTime, dT_lev / 2), StateList, OldStateList, SynchList_cor);
      }
    }
    PhysTime += dT_lev;
    //     if(myrank==0) cout<<"level now = "<<lev+1<<", "<<fgt(PhysTime-dT_lev,StartTime,dT_lev/2)<<endl;
    RestrictProlong(lev + 1, 1, fgt(PhysTime - dT_lev, StartTime, dT_lev / 2), StateList, OldStateList, SynchList_cor);
  }

  lev = GH->mylev;
  if (lev == -1)
    lev = 0;
  else
    lev = GH->levels - 1;

  {
    MPI_Status status;
    // receive
    if (lev == 0)
    {
      if (myrank == GH->start_rank[lev])
      {
        MPI_Recv(tporgo, 3 * BH_num, MPI_DOUBLE, GH->start_rank[GH->levels - 1], 1, MPI_COMM_WORLD, &status);
        //	 cout<<tporgo[0]<<","<<tporgo[1]<<","<<tporgo[2]<<endl;
      }
      else
        for (int i = 0; i < 3 * BH_num; i++)
          tporgo[i] = 0;
      MPI_Allreduce(tporgo, tporg, 3 * BH_num, MPI_DOUBLE, MPI_SUM, GH->Commlev[lev]);

      for (int i = 0; i < BH_num; i++)
        for (int j = 0; j < 3; j++)
          Porg0[i][j] = tporg[3 * i + j];

      //      if(myrank==GH->start_rank[lev]) cout<<Porg0[0][0]<<","<<Porg0[0][1]<<","<<Porg0[0][2]<<endl;
    }
    // send
    else if (myrank == GH->start_rank[lev])
    {
      for (int i = 0; i < BH_num; i++)
        for (int j = 0; j < 3; j++)
          tporg[3 * i + j] = Porg0[i][j];

      MPI_Send(tporg, 3 * BH_num, MPI_DOUBLE, GH->start_rank[0], 1, MPI_COMM_WORLD);
    }
  }

  delete[] tporg;
  delete[] tporgo;
}

//================================================================================================




//================================================================================================

// 该成员函数设定了时间演化过程中的不同层网格间递归时间演化

//================================================================================================

void bssn_class::RecursiveStep(int lev, int num) // in all 2^(lev+1)-1 steps
{
  if (trfls > 0)
    cout << "error: bssn_class::RecursiveStep does not support trfls > 0 yet" << endl;

  if (num / 2 * 2 == num)
    RecursiveStep(lev - 1, num / 2);
  else
  {
    Step(lev, 0);
    double dT_lev = dT * pow(0.5, Mymax(lev + 1, trfls));
    if (myrank == 0)
      cout << "level now = " << lev + 1 << ", " << (num - 1) % 2 << ", " 
           << fgt(PhysTime - dT_lev, StartTime, dT_lev / 2) << endl;
    RestrictProlong(lev + 1, (num - 1) % 2, fgt(PhysTime - dT_lev, StartTime, dT_lev / 2), StateList, OldStateList, SynchList_cor);
  }
}
#endif

//================================================================================================




//================================================================================================

// 该成员函数设定了时间演化过程中的每层网格的单步时间演化
// 针对 PSTR == 0 的情形

//================================================================================================

#if (PSTR == 0)
#if 1
void bssn_class::Step(int lev, int YN)
{
  setpbh(BH_num, Porg0, Mass, BH_num_input);

  double dT_lev = dT * pow(0.5, Mymax(lev, trfls));

// new code 2013-2-15, zjcao
#if (MAPBH == 1)
  // for black hole position
  if (BH_num > 0 && lev == GH->levels - 1)
  {
    compute_Porg_rhs(Porg0, Porg_rhs, Sfx0, Sfy0, Sfz0, lev);
    for (int ithBH = 0; ithBH < BH_num; ithBH++)
    {
      for (int ith = 0; ith < 3; ith++)
        Porg1[ithBH][ith] = Porg0[ithBH][ith] + Porg_rhs[ithBH][ith] * dT_lev;
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
          ErrorMonitor->outfile << "predictor step finds NaN for BH's position from ("
                                << Porg0[ithBH][0] << "," << Porg0[ithBH][1] << "," << Porg0[ithBH][2] << ")" << endl;

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
    AnalysisStuff(lev, dT_lev);
  }
#endif

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

        if (f_compute_rhs_bssn(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                               cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                               cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                               cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                               cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                               cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                               cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                               cg->fgfs[Lap0->sgfn], 
                               cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                               cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn],
                               cg->fgfs[phi_rhs->sgfn], cg->fgfs[trK_rhs->sgfn],
                               cg->fgfs[gxx_rhs->sgfn], cg->fgfs[gxy_rhs->sgfn], cg->fgfs[gxz_rhs->sgfn],
                               cg->fgfs[gyy_rhs->sgfn], cg->fgfs[gyz_rhs->sgfn], cg->fgfs[gzz_rhs->sgfn],
                               cg->fgfs[Axx_rhs->sgfn], cg->fgfs[Axy_rhs->sgfn], cg->fgfs[Axz_rhs->sgfn],
                               cg->fgfs[Ayy_rhs->sgfn], cg->fgfs[Ayz_rhs->sgfn], cg->fgfs[Azz_rhs->sgfn],
                               cg->fgfs[Gmx_rhs->sgfn], cg->fgfs[Gmy_rhs->sgfn], cg->fgfs[Gmz_rhs->sgfn],
                               cg->fgfs[Lap_rhs->sgfn], 
                               cg->fgfs[Sfx_rhs->sgfn], cg->fgfs[Sfy_rhs->sgfn], cg->fgfs[Sfz_rhs->sgfn],
                               cg->fgfs[dtSfx_rhs->sgfn], cg->fgfs[dtSfy_rhs->sgfn], cg->fgfs[dtSfz_rhs->sgfn],
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
#if (SommerType == 0)
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
                                dT_lev, 
                                cg->fgfs[phi0->sgfn],
                                cg->fgfs[Lap0->sgfn], 
                                cg->fgfs[varl0->data->sgfn], cg->fgfs[varl->data->sgfn], 
                                varl0->data->SoA,
                                Symmetry, cor);

#if (SommerType == 1)
#warning "shell part still bam type"
            if (lev == 0) // Shibata type sommerfeld
              f_sommerfeld_rout(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], 
                                Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                                dT_lev, 
                                cg->fgfs[phi0->sgfn],
                                cg->fgfs[Lap0->sgfn], 
                                cg->fgfs[varl0->data->sgfn], cg->fgfs[varl->data->sgfn], 
                                varl0->data->SoA,
                                Symmetry, pre);
#endif

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
#if (AGM == 0)
          f_enforce_ga(cg->shape,
                       cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                       cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                       cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                       cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn]);
#endif

          if (f_compute_rhs_bssn_ss(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
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
                                    cg->fgfs[phi_rhs->sgfn], cg->fgfs[trK_rhs->sgfn],
                                    cg->fgfs[gxx_rhs->sgfn], cg->fgfs[gxy_rhs->sgfn], cg->fgfs[gxz_rhs->sgfn],
                                    cg->fgfs[gyy_rhs->sgfn], cg->fgfs[gyz_rhs->sgfn], cg->fgfs[gzz_rhs->sgfn],
                                    cg->fgfs[Axx_rhs->sgfn], cg->fgfs[Axy_rhs->sgfn], cg->fgfs[Axz_rhs->sgfn],
                                    cg->fgfs[Ayy_rhs->sgfn], cg->fgfs[Ayz_rhs->sgfn], cg->fgfs[Azz_rhs->sgfn],
                                    cg->fgfs[Gmx_rhs->sgfn], cg->fgfs[Gmy_rhs->sgfn], cg->fgfs[Gmz_rhs->sgfn],
                                    cg->fgfs[Lap_rhs->sgfn], 
                                    cg->fgfs[Sfx_rhs->sgfn], cg->fgfs[Sfy_rhs->sgfn], cg->fgfs[Sfz_rhs->sgfn],
                                    cg->fgfs[dtSfx_rhs->sgfn], cg->fgfs[dtSfy_rhs->sgfn], cg->fgfs[dtSfz_rhs->sgfn],
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
            MyList<var> *varl0 = StateList, *varl = SynchList_pre, *varlrhs = RHSList; 
            // we do not check the correspondence here
            
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
#if 0 
// check rhs    
  {
         SH->Dump_Data(RHSList,0,PhysTime,dT_lev);
         if(myrank == 0)
	 {
            cout<<"check rhs"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	 }
  }
#endif
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

#if (MAPBH == 0)
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
                                << Porg0[ithBH][0] << "," << Porg0[ithBH][1] << "," << Porg0[ithBH][2] << ")" << endl;

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
    AnalysisStuff(lev, dT_lev);
  }
#endif

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

          if (f_compute_rhs_bssn(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                                 cg->fgfs[phi->sgfn], cg->fgfs[trK->sgfn],
                                 cg->fgfs[gxx->sgfn], cg->fgfs[gxy->sgfn], cg->fgfs[gxz->sgfn], 
                                 cg->fgfs[gyy->sgfn], cg->fgfs[gyz->sgfn], cg->fgfs[gzz->sgfn],
                                 cg->fgfs[Axx->sgfn], cg->fgfs[Axy->sgfn], cg->fgfs[Axz->sgfn], 
                                 cg->fgfs[Ayy->sgfn], cg->fgfs[Ayz->sgfn], cg->fgfs[Azz->sgfn],
                                 cg->fgfs[Gmx->sgfn], cg->fgfs[Gmy->sgfn], cg->fgfs[Gmz->sgfn],
                                 cg->fgfs[Lap->sgfn], 
                                 cg->fgfs[Sfx->sgfn], cg->fgfs[Sfy->sgfn], cg->fgfs[Sfz->sgfn],
                                 cg->fgfs[dtSfx->sgfn], cg->fgfs[dtSfy->sgfn], cg->fgfs[dtSfz->sgfn],
                                 cg->fgfs[phi1->sgfn], cg->fgfs[trK1->sgfn],
                                 cg->fgfs[gxx1->sgfn], cg->fgfs[gxy1->sgfn], cg->fgfs[gxz1->sgfn],
                                 cg->fgfs[gyy1->sgfn], cg->fgfs[gyz1->sgfn], cg->fgfs[gzz1->sgfn],
                                 cg->fgfs[Axx1->sgfn], cg->fgfs[Axy1->sgfn], cg->fgfs[Axz1->sgfn],
                                 cg->fgfs[Ayy1->sgfn], cg->fgfs[Ayz1->sgfn], cg->fgfs[Azz1->sgfn],
                                 cg->fgfs[Gmx1->sgfn], cg->fgfs[Gmy1->sgfn], cg->fgfs[Gmz1->sgfn],
                                 cg->fgfs[Lap1->sgfn], 
                                 cg->fgfs[Sfx1->sgfn], cg->fgfs[Sfy1->sgfn], cg->fgfs[Sfz1->sgfn],
                                 cg->fgfs[dtSfx1->sgfn], cg->fgfs[dtSfy1->sgfn], cg->fgfs[dtSfz1->sgfn],
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
            MyList<var> *varl0 = StateList, *varl = SynchList_pre, *varl1 = SynchList_cor, *varlrhs = RHSList; // we do not check the correspondence here
            while (varl0)
            {
#if (SommerType == 0)
#ifndef WithShell
              if (lev == 0) // sommerfeld indeed
                f_sommerfeld_routbam(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                     Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], 
                                     Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                                     cg->fgfs[varl1->data->sgfn],
                                     cg->fgfs[varl->data->sgfn], varl0->data->propspeed, varl0->data->SoA,
                                     Symmetry);
#endif
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
                                  dT_lev, 
                                  cg->fgfs[phi0->sgfn],
                                  cg->fgfs[Lap0->sgfn], 
                                  cg->fgfs[varl0->data->sgfn], cg->fgfs[varl1->data->sgfn], 
                                  varl0->data->SoA,
                                  Symmetry, cor);

#if (SommerType == 1)
              if (lev == 1) // shibata type sommerfeld
                f_sommerfeld_rout(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                  Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], 
                                  Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                                  dT_lev, 
                                  cg->fgfs[phi0->sgfn],
                                  cg->fgfs[Lap0->sgfn], 
                                  cg->fgfs[varl->data->sgfn], cg->fgfs[varl1->data->sgfn], 
                                  varl0->data->SoA,
                                  Symmetry, cor);
#endif

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

            if (f_compute_rhs_bssn_ss(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
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
                                      cg->fgfs[phi1->sgfn], cg->fgfs[trK1->sgfn],
                                      cg->fgfs[gxx1->sgfn], cg->fgfs[gxy1->sgfn], cg->fgfs[gxz1->sgfn],
                                      cg->fgfs[gyy1->sgfn], cg->fgfs[gyz1->sgfn], cg->fgfs[gzz1->sgfn],
                                      cg->fgfs[Axx1->sgfn], cg->fgfs[Axy1->sgfn], cg->fgfs[Axz1->sgfn],
                                      cg->fgfs[Ayy1->sgfn], cg->fgfs[Ayz1->sgfn], cg->fgfs[Azz1->sgfn],
                                      cg->fgfs[Gmx1->sgfn], cg->fgfs[Gmy1->sgfn], cg->fgfs[Gmz1->sgfn],
                                      cg->fgfs[Lap1->sgfn], 
                                      cg->fgfs[Sfx1->sgfn], cg->fgfs[Sfy1->sgfn], cg->fgfs[Sfz1->sgfn],
                                      cg->fgfs[dtSfx1->sgfn], cg->fgfs[dtSfy1->sgfn], cg->fgfs[dtSfz1->sgfn],
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
              MyList<var> *varl0 = StateList, *varl = SynchList_pre, *varl1 = SynchList_cor, *varlrhs = RHSList; 
              // we do not check the correspondence here
              
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
          ErrorMonitor->outfile << "find NaN on Shell Patches in RK4 substep#" 
                                << iter_count << " variables at t = " 
                                << PhysTime << endl;
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

#if (MAPBH == 0)
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

#if (MAPBH == 0)
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
#endif
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
#if 0
// check StateList   
  {
         SH->Dump_Data(StateList,0,PhysTime,dT_lev);
         if(myrank == 0)
	 {
            cout<<"check StateList"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	 }
  }
#endif
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

// 该成员函数设定了时间演化过程中的每层网格的单步时间演化（另一份）

//================================================================================================

// ICN for bam comparison

#else
void bssn_class::Step(int lev, int YN)
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

        if (f_compute_rhs_bssn(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                               cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                               cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                               cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                               cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                               cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                               cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                               cg->fgfs[Lap0->sgfn], 
                               cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                               cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn],
                               cg->fgfs[phi_rhs->sgfn], cg->fgfs[trK_rhs->sgfn],
                               cg->fgfs[gxx_rhs->sgfn], cg->fgfs[gxy_rhs->sgfn], cg->fgfs[gxz_rhs->sgfn],
                               cg->fgfs[gyy_rhs->sgfn], cg->fgfs[gyz_rhs->sgfn], cg->fgfs[gzz_rhs->sgfn],
                               cg->fgfs[Axx_rhs->sgfn], cg->fgfs[Axy_rhs->sgfn], cg->fgfs[Axz_rhs->sgfn],
                               cg->fgfs[Ayy_rhs->sgfn], cg->fgfs[Ayz_rhs->sgfn], cg->fgfs[Azz_rhs->sgfn],
                               cg->fgfs[Gmx_rhs->sgfn], cg->fgfs[Gmy_rhs->sgfn], cg->fgfs[Gmz_rhs->sgfn],
                               cg->fgfs[Lap_rhs->sgfn], 
                               cg->fgfs[Sfx_rhs->sgfn], cg->fgfs[Sfy_rhs->sgfn], cg->fgfs[Sfz_rhs->sgfn],
                               cg->fgfs[dtSfx_rhs->sgfn], cg->fgfs[dtSfy_rhs->sgfn], cg->fgfs[dtSfz_rhs->sgfn],
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
          MyList<var> *varl0 = StateList, *varl = SynchList_pre, *varlrhs = RHSList; 
          // we do not check the correspondence here
          while (varl0)
          {
#ifndef WithShell
            if (lev == 0) // sommerfeld indeed
              f_sommerfeld_routbam(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                   Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], 
                                   Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                                   cg->fgfs[varlrhs->data->sgfn],
                                   cg->fgfs[varl0->data->sgfn], varl0->data->propspeed, varl0->data->SoA,
                                   Symmetry);

#endif
            f_icn_rout(cg->shape, dT_lev, 
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
                                dT_lev, 
                                cg->fgfs[phi0->sgfn],
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

          if (f_compute_rhs_bssn_ss(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
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
                                    cg->fgfs[phi_rhs->sgfn], cg->fgfs[trK_rhs->sgfn],
                                    cg->fgfs[gxx_rhs->sgfn], cg->fgfs[gxy_rhs->sgfn], cg->fgfs[gxz_rhs->sgfn],
                                    cg->fgfs[gyy_rhs->sgfn], cg->fgfs[gyz_rhs->sgfn], cg->fgfs[gzz_rhs->sgfn],
                                    cg->fgfs[Axx_rhs->sgfn], cg->fgfs[Axy_rhs->sgfn], cg->fgfs[Axz_rhs->sgfn],
                                    cg->fgfs[Ayy_rhs->sgfn], cg->fgfs[Ayz_rhs->sgfn], cg->fgfs[Azz_rhs->sgfn],
                                    cg->fgfs[Gmx_rhs->sgfn], cg->fgfs[Gmy_rhs->sgfn], cg->fgfs[Gmz_rhs->sgfn],
                                    cg->fgfs[Lap_rhs->sgfn], 
                                    cg->fgfs[Sfx_rhs->sgfn], cg->fgfs[Sfy_rhs->sgfn], cg->fgfs[Sfz_rhs->sgfn],
                                    cg->fgfs[dtSfx_rhs->sgfn], cg->fgfs[dtSfy_rhs->sgfn], cg->fgfs[dtSfz_rhs->sgfn],
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

              f_icn_rout(cg->shape, dT_lev, 
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
#if 0 
// check rhs    
  {
         SH->Dump_Data(RHSList,0,PhysTime,dT_lev);
         if(myrank == 0)
	 {
            cout<<"check rhs"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	 }
  }
#endif
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
        ErrorMonitor->outfile << "find NaN in state variables on Shell Patches at t = " 
                              << PhysTime << endl;
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
      f_icn_scalar(dT_lev, Porg0[ithBH][0], Porg[ithBH][0], Porg_rhs[ithBH][0], iter_count);
      f_icn_scalar(dT_lev, Porg0[ithBH][1], Porg[ithBH][1], Porg_rhs[ithBH][1], iter_count);
      f_icn_scalar(dT_lev, Porg0[ithBH][2], Porg[ithBH][2], Porg_rhs[ithBH][2], iter_count);
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
    AnalysisStuff(lev, dT_lev);
  }
  // corrector
  for (iter_count = 1; iter_count < 3; iter_count++)
  {
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

          if (f_compute_rhs_bssn(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                                 cg->fgfs[phi->sgfn], cg->fgfs[trK->sgfn],
                                 cg->fgfs[gxx->sgfn], cg->fgfs[gxy->sgfn], cg->fgfs[gxz->sgfn], 
                                 cg->fgfs[gyy->sgfn], cg->fgfs[gyz->sgfn], cg->fgfs[gzz->sgfn],
                                 cg->fgfs[Axx->sgfn], cg->fgfs[Axy->sgfn], cg->fgfs[Axz->sgfn], 
                                 cg->fgfs[Ayy->sgfn], cg->fgfs[Ayz->sgfn], cg->fgfs[Azz->sgfn],
                                 cg->fgfs[Gmx->sgfn], cg->fgfs[Gmy->sgfn], cg->fgfs[Gmz->sgfn],
                                 cg->fgfs[Lap->sgfn], 
                                 cg->fgfs[Sfx->sgfn], cg->fgfs[Sfy->sgfn], cg->fgfs[Sfz->sgfn],
                                 cg->fgfs[dtSfx->sgfn], cg->fgfs[dtSfy->sgfn], cg->fgfs[dtSfz->sgfn],
                                 cg->fgfs[phi1->sgfn], cg->fgfs[trK1->sgfn],
                                 cg->fgfs[gxx1->sgfn], cg->fgfs[gxy1->sgfn], cg->fgfs[gxz1->sgfn],
                                 cg->fgfs[gyy1->sgfn], cg->fgfs[gyz1->sgfn], cg->fgfs[gzz1->sgfn],
                                 cg->fgfs[Axx1->sgfn], cg->fgfs[Axy1->sgfn], cg->fgfs[Axz1->sgfn],
                                 cg->fgfs[Ayy1->sgfn], cg->fgfs[Ayz1->sgfn], cg->fgfs[Azz1->sgfn],
                                 cg->fgfs[Gmx1->sgfn], cg->fgfs[Gmy1->sgfn], cg->fgfs[Gmz1->sgfn],
                                 cg->fgfs[Lap1->sgfn], 
                                 cg->fgfs[Sfx1->sgfn], cg->fgfs[Sfy1->sgfn], cg->fgfs[Sfz1->sgfn],
                                 cg->fgfs[dtSfx1->sgfn], cg->fgfs[dtSfy1->sgfn], cg->fgfs[dtSfz1->sgfn],
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
              f_icn_rout(cg->shape, dT_lev, 
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
                                  dT_lev, 
                                  cg->fgfs[phi0->sgfn],
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

            if (f_compute_rhs_bssn_ss(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
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
                                      cg->fgfs[phi1->sgfn], cg->fgfs[trK1->sgfn],
                                      cg->fgfs[gxx1->sgfn], cg->fgfs[gxy1->sgfn], cg->fgfs[gxz1->sgfn],
                                      cg->fgfs[gyy1->sgfn], cg->fgfs[gyz1->sgfn], cg->fgfs[gzz1->sgfn],
                                      cg->fgfs[Axx1->sgfn], cg->fgfs[Axy1->sgfn], cg->fgfs[Axz1->sgfn],
                                      cg->fgfs[Ayy1->sgfn], cg->fgfs[Ayz1->sgfn], cg->fgfs[Azz1->sgfn],
                                      cg->fgfs[Gmx1->sgfn], cg->fgfs[Gmy1->sgfn], cg->fgfs[Gmz1->sgfn],
                                      cg->fgfs[Lap1->sgfn], 
                                      cg->fgfs[Sfx1->sgfn], cg->fgfs[Sfy1->sgfn], cg->fgfs[Sfz1->sgfn],
                                      cg->fgfs[dtSfx1->sgfn], cg->fgfs[dtSfy1->sgfn], cg->fgfs[dtSfz1->sgfn],
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
              MyList<var> *varl0 = StateList, *varl = SynchList_pre, *varl1 = SynchList_cor, *varlrhs = RHSList; 
              // we do not check the correspondence here
              
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
        f_icn_scalar(dT_lev, Porg0[ithBH][0], Porg1[ithBH][0], Porg_rhs[ithBH][0], iter_count);
        f_icn_scalar(dT_lev, Porg0[ithBH][1], Porg1[ithBH][1], Porg_rhs[ithBH][1], iter_count);
        f_icn_scalar(dT_lev, Porg0[ithBH][2], Porg1[ithBH][2], Porg_rhs[ithBH][2], iter_count);
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
#if 0
// check StateList   
  {
         SH->Dump_Data(StateList,0,PhysTime,dT_lev);
         if(myrank == 0)
	 {
            cout<<"check StateList"<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	 }
  }
#endif
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
#endif

//================================================================================================



//================================================================================================

// 该成员函数设定了时间演化过程中的每层网格的单步时间演化
// 针对 PSTR == 0 以外的情形

//================================================================================================

#elif (PSTR == 1 || PSTR == 2 || PSTR == 3)
void bssn_class::Step(int lev, int YN)
{
  //   misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],"start Step");

  setpbh(BH_num, Porg0, Mass, BH_num_input);

  double dT_lev = dT * pow(0.5, Mymax(lev, trfls));

// new code 2013-2-15, zjcao
#if (MAPBH == 1)
  // for black hole position
  if (BH_num > 0 && lev == GH->levels - 1)
  {
    compute_Porg_rhs(Porg0, Porg_rhs, Sfx0, Sfy0, Sfz0, lev);
    for (int ithBH = 0; ithBH < BH_num; ithBH++)
    {
      for (int ith = 0; ith < 3; ith++)
        Porg1[ithBH][ith] = Porg0[ithBH][ith] + Porg_rhs[ithBH][ith] * dT_lev;
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
#endif

  //   misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],"before Predictor");

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

        if (f_compute_rhs_bssn(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                               cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                               cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                               cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                               cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                               cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                               cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                               cg->fgfs[Lap0->sgfn], 
                               cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                               cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn],
                               cg->fgfs[phi_rhs->sgfn], cg->fgfs[trK_rhs->sgfn],
                               cg->fgfs[gxx_rhs->sgfn], cg->fgfs[gxy_rhs->sgfn], cg->fgfs[gxz_rhs->sgfn],
                               cg->fgfs[gyy_rhs->sgfn], cg->fgfs[gyz_rhs->sgfn], cg->fgfs[gzz_rhs->sgfn],
                               cg->fgfs[Axx_rhs->sgfn], cg->fgfs[Axy_rhs->sgfn], cg->fgfs[Axz_rhs->sgfn],
                               cg->fgfs[Ayy_rhs->sgfn], cg->fgfs[Ayz_rhs->sgfn], cg->fgfs[Azz_rhs->sgfn],
                               cg->fgfs[Gmx_rhs->sgfn], cg->fgfs[Gmy_rhs->sgfn], cg->fgfs[Gmz_rhs->sgfn],
                               cg->fgfs[Lap_rhs->sgfn], 
                               cg->fgfs[Sfx_rhs->sgfn], cg->fgfs[Sfy_rhs->sgfn], cg->fgfs[Sfz_rhs->sgfn],
                               cg->fgfs[dtSfx_rhs->sgfn], cg->fgfs[dtSfy_rhs->sgfn], cg->fgfs[dtSfz_rhs->sgfn],
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
#if (SommerType == 0)
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
                                dT_lev, 
                                cg->fgfs[phi0->sgfn],
                                cg->fgfs[Lap0->sgfn], 
                                cg->fgfs[varl0->data->sgfn], cg->fgfs[varl->data->sgfn], 
                                varl0->data->SoA,
                                Symmetry, cor);

#if (SommerType == 1)
#warning "shell part still bam type"
            if (lev == 0) // Shibata type sommerfeld
              f_sommerfeld_rout(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], 
                                Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                                dT_lev, 
                                cg->fgfs[phi0->sgfn],
                                cg->fgfs[Lap0->sgfn], 
                                cg->fgfs[varl0->data->sgfn], cg->fgfs[varl->data->sgfn], 
                                varl0->data->SoA,
                                Symmetry, pre);
#endif

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

  //   misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],"after Predictor rhs calculation");

  // check error information
  {
    int erh = ERROR;
    MPI_Allreduce(&erh, &ERROR, 1, MPI_INT, MPI_SUM, GH->Commlev[lev]);
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

  //   misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],"before Predictor sync");

  Parallel::Sync(GH->PatL[lev], SynchList_pre, Symmetry);

#if (MAPBH == 0)
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
                                << Porg0[ithBH][0] << "," << Porg0[ithBH][1] << "," << Porg0[ithBH][2] << ")" << endl;

        MyList<var> *DG_List = new MyList<var>(Sfx0);
        DG_List->insert(Sfx0);
        DG_List->insert(Sfy0);
        DG_List->insert(Sfz0);
        Parallel::Dump_Data(GH->PatL[lev], DG_List, 0, PhysTime, dT_lev);
        DG_List->clearList();
      }
    }
  }
#endif

  //   misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],"before Corrector");

  // corrector
  for (iter_count = 1; iter_count < 4; iter_count++)
  {
    //   misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],"head of Corrector");

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

          if (f_compute_rhs_bssn(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                                 cg->fgfs[phi->sgfn], cg->fgfs[trK->sgfn],
                                 cg->fgfs[gxx->sgfn], cg->fgfs[gxy->sgfn], cg->fgfs[gxz->sgfn], 
                                 cg->fgfs[gyy->sgfn], cg->fgfs[gyz->sgfn], cg->fgfs[gzz->sgfn],
                                 cg->fgfs[Axx->sgfn], cg->fgfs[Axy->sgfn], cg->fgfs[Axz->sgfn], 
                                 cg->fgfs[Ayy->sgfn], cg->fgfs[Ayz->sgfn], cg->fgfs[Azz->sgfn],
                                 cg->fgfs[Gmx->sgfn], cg->fgfs[Gmy->sgfn], cg->fgfs[Gmz->sgfn],
                                 cg->fgfs[Lap->sgfn], 
                                 cg->fgfs[Sfx->sgfn], cg->fgfs[Sfy->sgfn], cg->fgfs[Sfz->sgfn],
                                 cg->fgfs[dtSfx->sgfn], cg->fgfs[dtSfy->sgfn], cg->fgfs[dtSfz->sgfn],
                                 cg->fgfs[phi1->sgfn], cg->fgfs[trK1->sgfn],
                                 cg->fgfs[gxx1->sgfn], cg->fgfs[gxy1->sgfn], cg->fgfs[gxz1->sgfn],
                                 cg->fgfs[gyy1->sgfn], cg->fgfs[gyz1->sgfn], cg->fgfs[gzz1->sgfn],
                                 cg->fgfs[Axx1->sgfn], cg->fgfs[Axy1->sgfn], cg->fgfs[Axz1->sgfn],
                                 cg->fgfs[Ayy1->sgfn], cg->fgfs[Ayz1->sgfn], cg->fgfs[Azz1->sgfn],
                                 cg->fgfs[Gmx1->sgfn], cg->fgfs[Gmy1->sgfn], cg->fgfs[Gmz1->sgfn],
                                 cg->fgfs[Lap1->sgfn], 
                                 cg->fgfs[Sfx1->sgfn], cg->fgfs[Sfy1->sgfn], cg->fgfs[Sfz1->sgfn],
                                 cg->fgfs[dtSfx1->sgfn], cg->fgfs[dtSfy1->sgfn], cg->fgfs[dtSfz1->sgfn],
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
#if (SommerType == 0)
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
                                  dT_lev, 
                                  cg->fgfs[phi0->sgfn],
                                  cg->fgfs[Lap0->sgfn], 
                                  cg->fgfs[varl0->data->sgfn], cg->fgfs[varl1->data->sgfn], 
                                  varl0->data->SoA,
                                  Symmetry, cor);

#if (SommerType == 1)
              if (lev == 1) // shibata type sommerfeld
                f_sommerfeld_rout(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                  Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], 
                                  Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                                  dT_lev, 
                                  cg->fgfs[phi0->sgfn],
                                  cg->fgfs[Lap0->sgfn], 
                                  cg->fgfs[varl->data->sgfn], cg->fgfs[varl1->data->sgfn], 
                                  varl0->data->SoA,
                                  Symmetry, cor);
#endif

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

    //   misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],"before Corrector error check");

    // check error information
    {
      int erh = ERROR;
      MPI_Allreduce(&erh, &ERROR, 1, MPI_INT, MPI_SUM, GH->Commlev[lev]);
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

    //    misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],"before Corrector sync");

    Parallel::Sync(GH->PatL[lev], SynchList_cor, Symmetry);

    //    misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],"after Corrector sync");

#if (MAPBH == 0)
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
//    misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],"after Corrector of black hole position");
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
      //    misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],"after pre cor swap");

#if (MAPBH == 0)
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
#endif
    }
    //    misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],"tail of corrector");
  }
#if (RPS == 0)
  // mesh refinement boundary part
  //   misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],"before RestrictProlong");
  RestrictProlong(lev, YN, BB);
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
  // for black hole position
  if (BH_num > 0 && lev == GH->levels - 1)
  {
    for (int ithBH = 0; ithBH < BH_num; ithBH++)
    {
      Porg0[ithBH][0] = Porg1[ithBH][0];
      Porg0[ithBH][1] = Porg1[ithBH][1];
      Porg0[ithBH][2] = Porg1[ithBH][2];
      //  if(myrank==GH->start_rank[lev]) 
      //     cout<<Porg0[ithBH][0]<<","<<Porg0[ithBH][1]<<","<<Porg0[ithBH][2]<<endl;
    }
  }

  //     if(myrank==GH->start_rank[lev]) cout<<GH->mylev<<endl;
  //     misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],"complet GH Step");
}

//================================================================================================



//================================================================================================

// 该成员函数设定了时间演化过程中球壳网格部分的单步时间演化

//================================================================================================

#ifdef WithShell
void bssn_class::SHStep()
{
  int lev = 0;
  // #if (PSTR == 1 || PSTR == 2)
  //    misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],"start Step");
  // #endif

  setpbh(BH_num, Porg0, Mass, BH_num_input);

  double dT_lev = dT * pow(0.5, Mymax(lev, trfls));

  // #if (PSTR == 1 || PSTR == 2)
  //    misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],"before Predictor");
  // #endif

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

        if (f_compute_rhs_bssn_ss(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
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
                                  cg->fgfs[phi_rhs->sgfn], cg->fgfs[trK_rhs->sgfn],
                                  cg->fgfs[gxx_rhs->sgfn], cg->fgfs[gxy_rhs->sgfn], cg->fgfs[gxz_rhs->sgfn],
                                  cg->fgfs[gyy_rhs->sgfn], cg->fgfs[gyz_rhs->sgfn], cg->fgfs[gzz_rhs->sgfn],
                                  cg->fgfs[Axx_rhs->sgfn], cg->fgfs[Axy_rhs->sgfn], cg->fgfs[Axz_rhs->sgfn],
                                  cg->fgfs[Ayy_rhs->sgfn], cg->fgfs[Ayz_rhs->sgfn], cg->fgfs[Azz_rhs->sgfn],
                                  cg->fgfs[Gmx_rhs->sgfn], cg->fgfs[Gmy_rhs->sgfn], cg->fgfs[Gmz_rhs->sgfn],
                                  cg->fgfs[Lap_rhs->sgfn], 
                                  cg->fgfs[Sfx_rhs->sgfn], cg->fgfs[Sfy_rhs->sgfn], cg->fgfs[Sfz_rhs->sgfn],
                                  cg->fgfs[dtSfx_rhs->sgfn], cg->fgfs[dtSfy_rhs->sgfn], cg->fgfs[dtSfz_rhs->sgfn],
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
          MyList<var> *varl0 = StateList, *varl = SynchList_pre, *varlrhs = RHSList; 
          // we do not check the correspondence here
          
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

#if (PSTR == 1 || PSTR == 2)
//   misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],"before Predictor's error check");
#endif
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

  // corrector
  for (iter_count = 1; iter_count < 4; iter_count++)
  {
    // for RK4: t0, t0+dt/2, t0+dt/2, t0+dt;
    if (iter_count == 1 || iter_count == 3)
      TRK4 += dT_lev / 2;

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

            if (f_compute_rhs_bssn_ss(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
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
                                      cg->fgfs[phi1->sgfn], cg->fgfs[trK1->sgfn],
                                      cg->fgfs[gxx1->sgfn], cg->fgfs[gxy1->sgfn], cg->fgfs[gxz1->sgfn],
                                      cg->fgfs[gyy1->sgfn], cg->fgfs[gyz1->sgfn], cg->fgfs[gzz1->sgfn],
                                      cg->fgfs[Axx1->sgfn], cg->fgfs[Axy1->sgfn], cg->fgfs[Axz1->sgfn],
                                      cg->fgfs[Ayy1->sgfn], cg->fgfs[Ayz1->sgfn], cg->fgfs[Azz1->sgfn],
                                      cg->fgfs[Gmx1->sgfn], cg->fgfs[Gmy1->sgfn], cg->fgfs[Gmz1->sgfn],
                                      cg->fgfs[Lap1->sgfn], 
                                      cg->fgfs[Sfx1->sgfn], cg->fgfs[Sfy1->sgfn], cg->fgfs[Sfz1->sgfn],
                                      cg->fgfs[dtSfx1->sgfn], cg->fgfs[dtSfy1->sgfn], cg->fgfs[dtSfz1->sgfn],
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
              MyList<var> *varl0 = StateList, *varl = SynchList_pre, *varl1 = SynchList_cor, *varlrhs = RHSList; 
              // we do not check the correspondence here
              
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
#if (RPS == 0)
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
  // note the data structure before update
  // SynchList_cor 1   -----------
  //
  // StateList     0   -----------
  //
  // OldStateList  old -----------
  // update
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
#endif

//================================================================================================



//================================================================================================

// 0: do not use mixing two levels data for OutBD; 1: do use

#define MIXOUTB 0
void bssn_class::RestrictProlong(int lev, int YN, bool BB,
                                 MyList<var> *SL, MyList<var> *OL, MyList<var> *corL)
// we assume
// StateList      1   -----------
//
// OldStateList   0   -----------
//
// SynchList_cor  old -----------
{
#if (PSTR == 1 || PSTR == 2)
//  stringstream a_stream;
//  a_stream.setf(ios::left);
#endif

  if (lev > 0)
  {
    MyList<Patch> *Pp, *Ppc;
    if (lev > trfls && YN == 0) // time refinement levels and for intermediat time level
    {
      Pp = GH->PatL[lev - 1];
      while (Pp)
      {
        if (BB)
          Parallel::prepare_inter_time_level(Pp->data, SL, OL, corL,
                                             SynchList_pre, 0); // use SynchList_pre as temporal storage space
        else
          Parallel::prepare_inter_time_level(Pp->data, SL, OL,
                                             SynchList_pre, 0); // use SynchList_pre as temporal storage space

#if (PSTR == 1 || PSTR == 2)
//	 Pp->data->checkPatch(0,GH->start_rank[GH->mylev]);
#endif
        Pp = Pp->next;
      }

#if (PSTR == 1 || PSTR == 2)
//       Pp=GH->PatL[lev];
//       while(Pp)
//       {
//	 Pp->data->checkPatch(0,GH->start_rank[GH->mylev]);
//	 Pp=Pp->next;
//       }

//       a_stream.clear();
//       a_stream.str("");
//       a_stream<<GH->mylev<<": 0 before Restrict";
//       misc::tillherecheck(GH->Commlev[GH->mylev],GH->start_rank[GH->mylev],a_stream.str());
#endif

#if (RPB == 0)
      Parallel::Restrict(GH->PatL[lev - 1], GH->PatL[lev], SL, SynchList_pre, Symmetry);
#elif (RPB == 1)
      //       Parallel::Restrict_bam(GH->PatL[lev-1],GH->PatL[lev],SL,SynchList_pre,Symmetry);
      Parallel::Restrict_bam(GH->PatL[lev - 1], GH->PatL[lev], SL, SynchList_pre, GH->rsul[lev], Symmetry);
#endif

#if (PSTR == 1 || PSTR == 2)
//       a_stream.clear();
//       a_stream.str("");
//       a_stream<<GH->mylev<<": 0 after Restrict";
//       misc::tillherecheck(GH->Commlev[GH->mylev],GH->start_rank[GH->mylev],a_stream.str());
#endif

      Parallel::Sync(GH->PatL[lev - 1], SynchList_pre, Symmetry);

#if (PSTR == 1 || PSTR == 2)
//       a_stream.clear();
//       a_stream.str("");
//       a_stream<<GH->mylev<<": 0 after Sync";
//       misc::tillherecheck(GH->Commlev[GH->mylev],GH->start_rank[GH->mylev],a_stream.str());
#endif

#if (RPB == 0)
      Ppc = GH->PatL[lev - 1];
      while (Ppc)
      {
        Pp = GH->PatL[lev];
        while (Pp)
        {
#if (MIXOUTB == 0)
          Parallel::OutBdLow2Hi(Ppc->data, Pp->data, SynchList_pre, SL, Symmetry);
#elif (MIXOUTB == 1)
          Parallel::OutBdLow2Himix(Ppc->data, Pp->data, SynchList_pre, SL, Symmetry);
#endif
          Pp = Pp->next;
        }
        Ppc = Ppc->next;
      }
#elif (RPB == 1)
      //       Parallel::OutBdLow2Hi_bam(GH->PatL[lev-1],GH->PatL[lev],SynchList_pre,SL,Symmetry);
      Parallel::OutBdLow2Hi_bam(GH->PatL[lev - 1], GH->PatL[lev], SynchList_pre, SL, GH->bdsul[lev], Symmetry);
#endif

#if (PSTR == 1 || PSTR == 2)
//       a_stream.clear();
//       a_stream.str("");
//       a_stream<<GH->mylev<<": 0 after OutBdLow2Hi";
//       misc::tillherecheck(GH->Commlev[GH->mylev],GH->start_rank[GH->mylev],a_stream.str());
#endif
    }
    else // no time refinement levels and for all same time levels
    {

#if (PSTR == 1 || PSTR == 2)
//       a_stream.clear();
//       a_stream.str("");
//       a_stream<<GH->mylev<<": 1 before Restrict";
//       misc::tillherecheck(GH->Commlev[GH->mylev],GH->start_rank[GH->mylev],a_stream.str());
#endif

#if (RPB == 0)
      Parallel::Restrict(GH->PatL[lev - 1], GH->PatL[lev], SL, SL, Symmetry);
#elif (RPB == 1)
      //       Parallel::Restrict_bam(GH->PatL[lev-1],GH->PatL[lev],SL,SL,Symmetry);
      Parallel::Restrict_bam(GH->PatL[lev - 1], GH->PatL[lev], SL, SL, GH->rsul[lev], Symmetry);
#endif

#if (PSTR == 1 || PSTR == 2)
//       a_stream.clear();
//       a_stream.str("");
//       a_stream<<GH->mylev<<": 1 before Sync";
//       misc::tillherecheck(GH->Commlev[GH->mylev],GH->start_rank[GH->mylev],a_stream.str());
#endif

      Parallel::Sync(GH->PatL[lev - 1], SL, Symmetry);

#if (PSTR == 1 || PSTR == 2)
//       a_stream.clear();
//       a_stream.str("");
//       a_stream<<GH->mylev<<": 1 after Sync";
//       misc::tillherecheck(GH->Commlev[GH->mylev],GH->start_rank[GH->mylev],a_stream.str());
#endif

#if (RPB == 0)
      Ppc = GH->PatL[lev - 1];
      while (Ppc)
      {
        Pp = GH->PatL[lev];
        while (Pp)
        {
#if (MIXOUTB == 0)
          Parallel::OutBdLow2Hi(Ppc->data, Pp->data, SL, SL, Symmetry);
#elif (MIXOUTB == 1)
          Parallel::OutBdLow2Himix(Ppc->data, Pp->data, SL, SL, Symmetry);
#endif
          Pp = Pp->next;
        }
        Ppc = Ppc->next;
      }
#elif (RPB == 1)
      //       Parallel::OutBdLow2Hi_bam(GH->PatL[lev-1],GH->PatL[lev],SL,SL,Symmetry);
      Parallel::OutBdLow2Hi_bam(GH->PatL[lev - 1], GH->PatL[lev], SL, SL, GH->bdsul[lev], Symmetry);
#endif

#if (PSTR == 1 || PSTR == 2)
//       a_stream.clear();
//       a_stream.str("");
//       a_stream<<GH->mylev<<": 1 after OutBdLow2Hi";
//       misc::tillherecheck(GH->Commlev[GH->mylev],GH->start_rank[GH->mylev],a_stream.str());
#endif
    }

    Parallel::Sync(GH->PatL[lev], SL, Symmetry);

#if (PSTR == 1 || PSTR == 2)
//    a_stream.clear();
//    a_stream.str("");
//    a_stream<<GH->mylev<<": after Sync";
//    misc::tillherecheck(GH->Commlev[GH->mylev],GH->start_rank[GH->mylev],a_stream.str());
#endif
  }
}

//================================================================================================



//================================================================================================

// auxiliary operation, input lev means original lev-1

void bssn_class::RestrictProlong_aux(int lev, int YN, bool BB,
                                     MyList<var> *SL, MyList<var> *OL, MyList<var> *corL)
// we assume
// StateList      1   -----------
//
// OldStateList   0   -----------
//
// SynchList_cor  old -----------
{
  //  misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],"starting RestrictProlong_aux");

  if (lev >= GH->levels - 1)
    return;
  lev = lev + 1;

  if (lev > 0)
  {
    MyList<Patch> *Pp, *Ppc;
    if (lev > trfls && YN == 0) // time refinement levels and for intermediat time level
    {
      Pp = GH->PatL[lev - 1];
      while (Pp)
      {
        if (BB)
          Parallel::prepare_inter_time_level(Pp->data, SL, OL, corL,
                                             SynchList_pre, 0); // use SynchList_pre as temporal storage space
        else
          Parallel::prepare_inter_time_level(Pp->data, SL, OL,
                                             SynchList_pre, 0); // use SynchList_pre as temporal storage space
        Pp = Pp->next;
      }

#if (RPB == 0)
      Parallel::Restrict(GH->PatL[lev - 1], GH->PatL[lev], SL, SynchList_pre, Symmetry);
#elif (RPB == 1)
      //       Parallel::Restrict_bam(GH->PatL[lev-1],GH->PatL[lev],SL,SynchList_pre,Symmetry);
      Parallel::Restrict_bam(GH->PatL[lev - 1], GH->PatL[lev], SL, SynchList_pre, GH->rsul[lev], Symmetry);
#endif

      Parallel::Sync(GH->PatL[lev - 1], SynchList_pre, Symmetry);

#if (RPB == 0)
      Ppc = GH->PatL[lev - 1];
      while (Ppc)
      {
        Pp = GH->PatL[lev];
        while (Pp)
        {
#if (MIXOUTB == 0)
          Parallel::OutBdLow2Hi(Ppc->data, Pp->data, SynchList_pre, SL, Symmetry);
#elif (MIXOUTB == 1)
          Parallel::OutBdLow2Himix(Ppc->data, Pp->data, SynchList_pre, SL, Symmetry);
#endif
          Pp = Pp->next;
        }
        Ppc = Ppc->next;
      }
#elif (RPB == 1)
      //       Parallel::OutBdLow2Hi_bam(GH->PatL[lev-1],GH->PatL[lev],SynchList_pre,SL,Symmetry);
      Parallel::OutBdLow2Hi_bam(GH->PatL[lev - 1], GH->PatL[lev], SynchList_pre, SL, GH->bdsul[lev], Symmetry);
#endif
    }
    else // no time refinement levels and for all same time levels
    {
#if (RPB == 0)
      Parallel::Restrict(GH->PatL[lev - 1], GH->PatL[lev], SL, SL, Symmetry);
#elif (RPB == 1)
      //       Parallel::Restrict_bam(GH->PatL[lev-1],GH->PatL[lev],SL,SL,Symmetry);
      Parallel::Restrict_bam(GH->PatL[lev - 1], GH->PatL[lev], SL, SL, GH->rsul[lev], Symmetry);
#endif

      Parallel::Sync(GH->PatL[lev - 1], SL, Symmetry);

#if (RPB == 0)
      Ppc = GH->PatL[lev - 1];
      while (Ppc)
      {
        Pp = GH->PatL[lev];
        while (Pp)
        {
#if (MIXOUTB == 0)
          Parallel::OutBdLow2Hi(Ppc->data, Pp->data, SL, SL, Symmetry);
#elif (MIXOUTB == 1)
          Parallel::OutBdLow2Himix(Ppc->data, Pp->data, SL, SL, Symmetry);
#endif
          Pp = Pp->next;
        }
        Ppc = Ppc->next;
      }
#elif (RPB == 1)
      //       Parallel::OutBdLow2Hi_bam(GH->PatL[lev-1],GH->PatL[lev],SL,SL,Symmetry);
      Parallel::OutBdLow2Hi_bam(GH->PatL[lev - 1], GH->PatL[lev], SL, SL, GH->bdsul[lev], Symmetry);
#endif
    }

    Parallel::Sync(GH->PatL[lev], SL, Symmetry);
  }
}

//================================================================================================



//================================================================================================

void bssn_class::RestrictProlong(int lev, int YN, bool BB)
{
  double dT_lev = dT * pow(0.5, Mymax(lev, trfls));
  // we assume  for fine
  // SynchList_cor 1   -----------
  //
  // StateList     0   -----------
  //
  // OldStateList  old -----------
  //            for coarse
  // StateList      1   -----------
  //
  // OldStateList   0   -----------
  //
  // SynchList_cor  old -----------
  if (lev > 0)
  {
    MyList<Patch> *Pp, *Ppc;
    if (lev > trfls && YN == 0) // time refinement levels and for intermediat time level
    {
      if (myrank == 0)
        cout << "/=: " << GH->Lt[lev - 1] << "," << GH->Lt[lev] + dT_lev << endl;
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

#if (RPB == 0)
      Parallel::Restrict(GH->PatL[lev - 1], GH->PatL[lev], SynchList_cor, SynchList_pre, Symmetry);
#elif (RPB == 1)
      //       Parallel::Restrict_bam(GH->PatL[lev-1],GH->PatL[lev],SynchList_cor,SynchList_pre,Symmetry);
      Parallel::Restrict_bam(GH->PatL[lev - 1], GH->PatL[lev], SynchList_cor, SynchList_pre, GH->rsul[lev], Symmetry);
#endif

      Parallel::Sync(GH->PatL[lev - 1], SynchList_pre, Symmetry);

#if (RPB == 0)
      Ppc = GH->PatL[lev - 1];
      while (Ppc)
      {
        Pp = GH->PatL[lev];
        while (Pp)
        {
#if (MIXOUTB == 0)
          Parallel::OutBdLow2Hi(Ppc->data, Pp->data, SynchList_pre, SynchList_cor, Symmetry);
#elif (MIXOUTB == 1)
          Parallel::OutBdLow2Himix(Ppc->data, Pp->data, SynchList_pre, SynchList_cor, Symmetry);
#endif
          Pp = Pp->next;
        }
        Ppc = Ppc->next;
      }
#elif (RPB == 1)
      //       Parallel::OutBdLow2Hi_bam(GH->PatL[lev-1],GH->PatL[lev],SynchList_pre,SynchList_cor,Symmetry);
      Parallel::OutBdLow2Hi_bam(GH->PatL[lev - 1], GH->PatL[lev], SynchList_pre, SynchList_cor, GH->bdsul[lev], Symmetry);
#endif
    }
    else // no time refinement levels and for all same time levels
    {
      if (myrank == 0)
        cout << "===: " << GH->Lt[lev - 1] << "," << GH->Lt[lev] + dT_lev << endl;
#if (RPB == 0)
      Parallel::Restrict(GH->PatL[lev - 1], GH->PatL[lev], SynchList_cor, StateList, Symmetry);
#elif (RPB == 1)
      //       Parallel::Restrict_bam(GH->PatL[lev-1],GH->PatL[lev],SynchList_cor,StateList,Symmetry);
      Parallel::Restrict_bam(GH->PatL[lev - 1], GH->PatL[lev], SynchList_cor, StateList, GH->rsul[lev], Symmetry);
#endif

      Parallel::Sync(GH->PatL[lev - 1], StateList, Symmetry);

#if (RPB == 0)
      Ppc = GH->PatL[lev - 1];
      while (Ppc)
      {
        Pp = GH->PatL[lev];
        while (Pp)
        {
#if (MIXOUTB == 0)
          Parallel::OutBdLow2Hi(Ppc->data, Pp->data, StateList, SynchList_cor, Symmetry);
#elif (MIXOUTB == 1)
          Parallel::OutBdLow2Himix(Ppc->data, Pp->data, StateList, SynchList_cor, Symmetry);
#endif
          Pp = Pp->next;
        }
        Ppc = Ppc->next;
      }
#elif (RPB == 1)
      //       Parallel::OutBdLow2Hi_bam(GH->PatL[lev-1],GH->PatL[lev],StateList,SynchList_cor,Symmetry);
      Parallel::OutBdLow2Hi_bam(GH->PatL[lev - 1], GH->PatL[lev], StateList, SynchList_cor, GH->bdsul[lev], Symmetry);
#endif
    }

    Parallel::Sync(GH->PatL[lev], SynchList_cor, Symmetry);
  }
}

//================================================================================================



//================================================================================================

void bssn_class::ProlongRestrict(int lev, int YN, bool BB)
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

#if (RPB == 0)
      Ppc = GH->PatL[lev - 1];
      while (Ppc)
      {
        Pp = GH->PatL[lev];
        while (Pp)
        {
#if (MIXOUTB == 0)
          Parallel::OutBdLow2Hi(Ppc->data, Pp->data, SynchList_pre, SynchList_cor, Symmetry);
#elif (MIXOUTB == 1)
          Parallel::OutBdLow2Himix(Ppc->data, Pp->data, SynchList_pre, SynchList_cor, Symmetry);
#endif
          Pp = Pp->next;
        }
        Ppc = Ppc->next;
      }
#elif (RPB == 1)
      //       Parallel::OutBdLow2Hi_bam(GH->PatL[lev-1],GH->PatL[lev],SynchList_pre,SynchList_cor,Symmetry);
      Parallel::OutBdLow2Hi_bam(GH->PatL[lev - 1], GH->PatL[lev], SynchList_pre, SynchList_cor, GH->bdsul[lev], Symmetry);
#endif
    }
    else // no time refinement levels and for all same time levels
    {
#if (RPB == 0)
      Ppc = GH->PatL[lev - 1];
      while (Ppc)
      {
        Pp = GH->PatL[lev];
        while (Pp)
        {
#if (MIXOUTB == 0)
          Parallel::OutBdLow2Hi(Ppc->data, Pp->data, StateList, SynchList_cor, Symmetry);
#elif (MIXOUTB == 1)
          Parallel::OutBdLow2Himix(Ppc->data, Pp->data, StateList, SynchList_cor, Symmetry);
#endif
          Pp = Pp->next;
        }
        Ppc = Ppc->next;
      }
#elif (RPB == 1)
      //       Parallel::OutBdLow2Hi_bam(GH->PatL[lev-1],GH->PatL[lev],StateList,SynchList_cor,Symmetry);
      Parallel::OutBdLow2Hi_bam(GH->PatL[lev - 1], GH->PatL[lev], StateList, SynchList_cor, GH->bdsul[lev], Symmetry);
#endif

#if 0
#if (RPB == 0)
       Parallel::Restrict(GH->PatL[lev-1],GH->PatL[lev],SynchList_cor,StateList,Symmetry);
#elif (RPB == 1)
//       Parallel::Restrict_bam(GH->PatL[lev-1],GH->PatL[lev],SynchList_cor,StateList,Symmetry);
       Parallel::Restrict_bam(GH->PatL[lev-1],GH->PatL[lev],SynchList_cor,StateList,GH->rsul[lev],Symmetry);
#endif
#else
      Parallel::Restrict_after(GH->PatL[lev - 1], GH->PatL[lev], SynchList_cor, StateList, Symmetry);
#endif
      Parallel::Sync(GH->PatL[lev - 1], StateList, Symmetry);
    }

    Parallel::Sync(GH->PatL[lev], SynchList_cor, Symmetry);
  }
}
#undef MIXOUTB

//================================================================================================



//================================================================================================

// 该成员函数用于计算引力波 Psi4

//================================================================================================

void bssn_class::Compute_Psi4(int lev)
{
  MyList<var> *DG_List = new MyList<var>(Rpsi4);
  DG_List->insert(Ipsi4);

#if 0 // test showes this operation does not help    
for(int ilev = GH->levels-1;ilev>=lev;ilev--)
{
  MyList<Patch> *Pp=GH->PatL[ilev];
#else
  MyList<Patch> *Pp = GH->PatL[lev];
#endif
  while (Pp)
  {
    MyList<Block> *BP = Pp->data->blb;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
#if (Psi4type == 0)
        if (0) // if Gamma^i_jk and R_ij can be reused from the rhs calculation
          f_ricci_gamma(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                        cg->fgfs[phi0->sgfn],
                        cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                        cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                        cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                        cg->fgfs[Gamxxx->sgfn], cg->fgfs[Gamxxy->sgfn], cg->fgfs[Gamxxz->sgfn],
                        cg->fgfs[Gamxyy->sgfn], cg->fgfs[Gamxyz->sgfn], cg->fgfs[Gamxzz->sgfn],
                        cg->fgfs[Gamyxx->sgfn], cg->fgfs[Gamyxy->sgfn], cg->fgfs[Gamyxz->sgfn],
                        cg->fgfs[Gamyyy->sgfn], cg->fgfs[Gamyyz->sgfn], cg->fgfs[Gamyzz->sgfn],
                        cg->fgfs[Gamzxx->sgfn], cg->fgfs[Gamzxy->sgfn], cg->fgfs[Gamzxz->sgfn],
                        cg->fgfs[Gamzyy->sgfn], cg->fgfs[Gamzyz->sgfn], cg->fgfs[Gamzzz->sgfn],
                        cg->fgfs[Rxx->sgfn], cg->fgfs[Rxy->sgfn], cg->fgfs[Rxz->sgfn], 
                        cg->fgfs[Ryy->sgfn], cg->fgfs[Ryz->sgfn], cg->fgfs[Rzz->sgfn],
                        Symmetry);
        // the input arguments Gamma^i_jk and R_ij do not need synch, because we do not need to derivate them
        f_getnp4(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                 cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
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
        f_getnp4old(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                    cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
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

#if 0  
    Parallel::Sync(GH->PatL[ilev],DG_List,Symmetry);
}
// because of double level data change, you can not do this in above loop
// prolong restrict Psi4
for(int ilev=GH->levels-1;ilev>lev;ilev--)
    RestrictProlong(ilev,1,false,DG_List,DG_List,DG_List);
#else
  Parallel::Sync(GH->PatL[lev], DG_List, Symmetry);
#endif

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
          if (0) // if Gamma^i_jk and R_ij can be reused from the rhs calculation
            f_ricci_gamma_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                             cg->fgfs[fngfs + 
                             ShellPatch::gx], cg->fgfs[fngfs + ShellPatch::gy], 
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
                             cg->fgfs[phi0->sgfn],
                             cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                             cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                             cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                             cg->fgfs[Gamxxx->sgfn], cg->fgfs[Gamxxy->sgfn], cg->fgfs[Gamxxz->sgfn],
                             cg->fgfs[Gamxyy->sgfn], cg->fgfs[Gamxyz->sgfn], cg->fgfs[Gamxzz->sgfn],
                             cg->fgfs[Gamyxx->sgfn], cg->fgfs[Gamyxy->sgfn], cg->fgfs[Gamyxz->sgfn],
                             cg->fgfs[Gamyyy->sgfn], cg->fgfs[Gamyyz->sgfn], cg->fgfs[Gamyzz->sgfn],
                             cg->fgfs[Gamzxx->sgfn], cg->fgfs[Gamzxy->sgfn], cg->fgfs[Gamzxz->sgfn],
                             cg->fgfs[Gamzyy->sgfn], cg->fgfs[Gamzyz->sgfn], cg->fgfs[Gamzzz->sgfn],
                             cg->fgfs[Rxx->sgfn], cg->fgfs[Rxy->sgfn], cg->fgfs[Rxz->sgfn], 
                             cg->fgfs[Ryy->sgfn], cg->fgfs[Ryz->sgfn], cg->fgfs[Rzz->sgfn],
                             Symmetry, lev, Pp->data->sst);

          f_getnp4_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
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
          f_getnp4old_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
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

    SH->Synch(DG_List, Symmetry);
#if 0     
// interpolate Psi4
     SH->CS_Inter(DG_List,Symmetry);
#endif
  }
#endif

  DG_List->clearList();

  //    misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],"end of Compute_Psi4");
}

//================================================================================================



//================================================================================================

// 该成员函数用于设定初始时刻黑洞的 Puncture 位置

//================================================================================================

void bssn_class::Setup_Black_Hole_position()
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
          ErrorMonitor->outfile << "error reading parameter file " << filename << " in line " << i << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
      else if (status == 0)
        continue;

      if (sgrp == "BSSN" && skey == "BH_num")
      {
        BH_num_input = BH_num = atoi(sval.c_str());
        break;
      }
    }
    inf.close();
  }
  // set up the data for black holes
  // these arrays will be deleted when bssn_class is deleted
  Pmom = new double[3 * BH_num];
  Spin = new double[3 * BH_num];
  Mass = new double[BH_num];
  Porg0 = new double *[BH_num];
  Porgbr = new double *[BH_num];
  Porg = new double *[BH_num];
  Porg1 = new double *[BH_num];
  Porg_rhs = new double *[BH_num];
  for (int i = 0; i < BH_num; i++)
  {
    Porg0[i] = new double[3];
    Porgbr[i] = new double[3];
    Porg[i] = new double[3];
    Porg1[i] = new double[3];
    Porg_rhs[i] = new double[3];
  }
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
          ErrorMonitor->outfile << "error reading parameter file " << filename << " in line " << i << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
      else if (status == 0)
        continue;

      if (sgrp == "BSSN" && sind < BH_num)
      {
        if (skey == "Mass")
          Mass[sind] = atof(sval.c_str());
        else if (skey == "Porgx")
          Porg0[sind][0] = atof(sval.c_str());
        else if (skey == "Porgy")
          Porg0[sind][1] = atof(sval.c_str());
        else if (skey == "Porgz")
          Porg0[sind][2] = atof(sval.c_str());
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
  // echo information of Black holes
  if (myrank == 0)
  {
    cout <<                                                              endl;
    cout << " initial information of " << BH_num << " Black Hole(s) " << endl;
    cout << setw(12) << "Mass"
         << setw(12) << "x"
         << setw(12) << "y"
         << setw(12) << "z"
         << setw(16) << "Px"
         << setw(16) << "Py"
         << setw(12) << "Pz"
         << setw(12) << "Sx"
         << setw(12) << "Sy"
         << setw(12) << "Sz" << endl;
    for (int i = 0; i < BH_num; i++)
    {
      cout << setw(12) << Mass[i]
           << setw(12) << Porg0[i][0]
           << setw(12) << Porg0[i][1]
           << setw(12) << Porg0[i][2]
           << setw(16) << Pmom[i * 3]
           << setw(16) << Pmom[i * 3 + 1]
           << setw(12) << Pmom[i * 3 + 2]
           << setw(12) << Spin[i * 3]
           << setw(12) << Spin[i * 3 + 1]
           << setw(12) << Spin[i * 3 + 2] << endl;
    }
  }

  int maxl = 1;
  int levels;
  int *grids;
  double bbox[6];
  // read parameter from file
  {
    const int LEN = 256;
    char pline[LEN];
    string str, sgrp, skey, sval;
    int sind1, sind2, sind3;
    ifstream inf(filename, ifstream::in);
    if (!inf.good() && myrank == 0)
    {
      cout << "bssn_class::Setup_Black_Hole_position: Can not open parameter file " << filename 
           << " for inputing information of black holes" << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    for (int i = 1; inf.good(); i++)
    {
      inf.getline(pline, LEN);
      str = pline;

      int status = misc::parse_parts(str, sgrp, skey, sval, sind1);
      if (status == -1)
      {
        cout << "error reading parameter file " << filename << " in line " << i << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
      else if (status == 0)
        continue;

      if (sgrp == "cgh" && skey == "levels")
      {
        levels = atoi(sval.c_str());
        break;
      }
    }
    inf.close();
  }
  grids = new int[levels];
  // read parameter from file
  {
    const int LEN = 256;
    char pline[LEN];
    string str, sgrp, skey, sval;
    int sind1, sind2, sind3;
    ifstream inf(filename, ifstream::in);
    if (!inf.good() && myrank == 0)
    {
      cout << "bssn_class::Setup_Black_Hole_position: Can not open parameter file " << filename 
           << " for inputing information of black holes" << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    for (int i = 1; inf.good(); i++)
    {
      inf.getline(pline, LEN);
      str = pline;

      int status = misc::parse_parts(str, sgrp, skey, sval, sind1, sind2, sind3);
      if (status == -1)
      {
        cout << "error reading parameter file " << filename << " in line " << i << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
      else if (status == 0)
        continue;

      if (sgrp == "cgh" && skey == "grids" && sind1 < levels)
        grids[sind1] = atoi(sval.c_str());
      if (sgrp == "cgh" && skey == "bbox" && sind1 == 0 && sind2 == 0)
        bbox[sind3] = atof(sval.c_str());
    }
    inf.close();
  }
  for (int i = 0; i < levels; i++)
    if (maxl < grids[i])
      maxl = grids[i];

  delete[] grids;

  if (BH_num > maxl)
  {
    int BH_numc = BH_num;
    for (int i = 0; i < BH_num; i++)
      if (Porg0[i][0] < bbox[0] || Porg0[i][0] > bbox[3] ||
          Porg0[i][1] < bbox[1] || Porg0[i][1] > bbox[4] ||
          Porg0[i][2] < bbox[2] || Porg0[i][2] > bbox[5])
      {
        delete[] Porg0[i];
        Porg0[i] = 0;
        BH_numc--;
      }

    if (BH_num > BH_numc)
    {
      maxl = BH_numc;
      int bhi;
      double *tmp;

      tmp = Pmom;
      Pmom = new double[3 * maxl];
      bhi = 0;
      for (int i = 0; i < BH_num; i++)
        if (Porg0[i])
        {
          for (int j = 0; j < 3; j++)
            Pmom[3 * bhi + j] = tmp[3 * i + j];
          bhi++;
        }
      delete[] tmp;

      tmp = Spin;
      Spin = new double[3 * maxl];
      bhi = 0;
      for (int i = 0; i < BH_num; i++)
        if (Porg0[i])
        {
          for (int j = 0; j < 3; j++)
            Spin[3 * bhi + j] = tmp[3 * i + j];
          bhi++;
        }
      delete[] tmp;

      tmp = Mass;
      Mass = new double[3 * maxl];
      bhi = 0;
      for (int i = 0; i < BH_num; i++)
        if (Porg0[i])
        {
          Mass[bhi] = tmp[i];
          bhi++;
        }
      delete[] tmp;

      double **ttmp;
      ttmp = Porg0;
      Porg0 = new double *[maxl];
      bhi = 0;
      for (int i = 0; i < BH_num; i++)
        if (ttmp[i])
        {
          Porg0[bhi] = ttmp[i];
          bhi++;
        }
      delete[] ttmp;

      for (int i = 0; i < BH_num; i++)
      {
        delete[] Porgbr[i];
        delete[] Porg[i];
        delete[] Porg1[i];
        delete[] Porg_rhs[i];
      }
      delete[] Porgbr;
      delete[] Porg;
      delete[] Porg1;
      delete[] Porg_rhs;

      BH_num = maxl;

      Porgbr = new double *[BH_num];
      Porg = new double *[BH_num];
      Porg1 = new double *[BH_num];
      Porg_rhs = new double *[BH_num];

      for (int i = 0; i < BH_num; i++)
      {
        Porgbr[i] = new double[3];
        Porg[i] = new double[3];
        Porg1[i] = new double[3];
        Porg_rhs[i] = new double[3];
      }
    }
  }

  for (int i = 0; i < BH_num; i++)
  {
    for (int j = 0; j < dim; j++)
      Porgbr[i][j] = Porg0[i][j];
  }

  setpbh(BH_num, Porg0, Mass, BH_num_input);
}

//================================================================================================



//================================================================================================

// 该成员函数用于计算黑洞的位置

//================================================================================================

#if 0
// old code 

void bssn_class::compute_Porg_rhs(double **BH_PS,double **BH_RHS,var *forx,var *fory,var *forz,int lev)
{
  const int InList = 3;

  MyList<var> * DG_List=new MyList<var>(forx);
  DG_List->insert(fory); DG_List->insert(forz);

  int n;
  double *x1,*y1,*z1;
  double *shellf;
  shellf=new double[3*BH_num];
  double *pox[3];
  for(int i=0;i<3;i++) pox[i] = new double[BH_num];
  for( n = 0; n < BH_num; n++)
  {
    pox[0][n] = BH_PS[n][0];
    pox[1][n] = BH_PS[n][1];
    pox[2][n] = BH_PS[n][2];
  }

  if(!Parallel::PatList_Interp_Points(GH->PatL[lev],DG_List,BH_num,pox,shellf,Symmetry))
  {
	  ErrorMonitor->outfile<<"fail to find black holes at t = "<<PhysTime<<endl;
          for( n = 0; n < BH_num; n++) 
              ErrorMonitor->outfile<<"(x,y,z) = ("<<pox[0][n]<<","<<pox[1][n]<<","<<pox[2][n]<<")"<<endl;
  }
  
  for( n = 0; n < BH_num; n++)
  {
    BH_RHS[n][0]=-shellf[3*n  ];
    BH_RHS[n][1]=-shellf[3*n+1];
    BH_RHS[n][2]=-shellf[3*n+2];
  }
   
  DG_List->clearList();
  delete[] shellf;
  for(int i=0;i<3;i++) delete[] pox[i];
}

#else

// new code considering diferent levels for different black hole

void bssn_class::compute_Porg_rhs(double **BH_PS, double **BH_RHS, var *forx, var *fory, var *forz, int ilev)
{
  const int InList = 3;

  MyList<var> *DG_List = new MyList<var>(forx);
  DG_List->insert(fory);
  DG_List->insert(forz);

  double *x1, *y1, *z1;
  double *shellf;
  shellf = new double[3];
  double *pox[3];
  for (int i = 0; i < 3; i++)
    pox[i] = new double[1];

  for (int n = 0; n < BH_num; n++)
  {
    pox[0][0] = BH_PS[n][0];
    pox[1][0] = BH_PS[n][1];
    pox[2][0] = BH_PS[n][2];

    int lev = ilev;

#if (PSTR == 0)
    while (!Parallel::PatList_Interp_Points(GH->PatL[lev], DG_List, 1, pox, shellf, Symmetry))
#elif (PSTR == 1 || PSTR == 2 || PSTR == 3)
    while (!Parallel::PatList_Interp_Points(GH->PatL[lev], DG_List, 1, pox, shellf, Symmetry, GH->Commlev[lev]))
#endif
    {
      lev--;
      if (lev < 0)
      {
        ErrorMonitor->outfile << "fail to find black holes at t = " << PhysTime << endl;
        for (n = 0; n < BH_num; n++)
          ErrorMonitor->outfile << "(x,y,z) = (" 
                                << pox[0][n] << "," << pox[1][n] << "," << pox[2][n] 
                                << ")" << endl;
        break;
      }
    }

    if (lev >= 0)
    {
      BH_RHS[n][0] = -shellf[0];
      BH_RHS[n][1] = -shellf[1];
      BH_RHS[n][2] = -shellf[2];
    }
  }

  DG_List->clearList();
  delete[] shellf;
  for (int i = 0; i < 3; i++)
    delete[] pox[i];
}
#endif

//================================================================================================



//================================================================================================

// 该成员函数用于计算引力波相关数据

//================================================================================================

void bssn_class::AnalysisStuff(int lev, double dT_lev)
{
  LastAnas += dT_lev;

  if (LastAnas >= AnasTime)
  {
#ifdef Point_Psi4
#error "not support parallel levels yet"
    // Gam_ijk and R_ij have been calculated in Interp_Constraint()
    double SYM = 1, ANT = -1;
    for (int levh = lev; levh < GH->levels; levh++)
    {
      MyList<Patch> *Pp = GH->PatL[levh];
      while (Pp)
      {
        MyList<Block> *BP = Pp->data->blb;
        while (BP)
        {
          Block *cg = BP->data;
          if (myrank == cg->rank)
          {
            f_fderivs(cg->shape, cg->fgfs[phi0->sgfn], 
                      cg->fgfs[phix->sgfn], cg->fgfs[phiy->sgfn], cg->fgfs[phiz->sgfn],
                      cg->X[0], cg->X[1], cg->X[2], 
                      SYM, SYM, SYM, Symmetry, levh);
            f_fderivs(cg->shape, cg->fgfs[trK0->sgfn], 
                      cg->fgfs[trKx->sgfn], cg->fgfs[trKy->sgfn], cg->fgfs[trKz->sgfn],
                      cg->X[0], cg->X[1], cg->X[2], 
                      SYM, SYM, SYM, Symmetry, levh);
            f_fderivs(cg->shape, cg->fgfs[Axx0->sgfn], 
                      cg->fgfs[Axxx->sgfn], cg->fgfs[Axxy->sgfn], cg->fgfs[Axxz->sgfn],
                      cg->X[0], cg->X[1], cg->X[2], 
                      SYM, SYM, SYM, Symmetry, levh);
            f_fderivs(cg->shape, cg->fgfs[Axy0->sgfn], 
                      cg->fgfs[Axyx->sgfn], cg->fgfs[Axyy->sgfn], cg->fgfs[Axyz->sgfn],
                      cg->X[0], cg->X[1], cg->X[2], 
                      ANT, ANT, SYM, Symmetry, levh);
            f_fderivs(cg->shape, cg->fgfs[Axz0->sgfn], 
                      cg->fgfs[Axzx->sgfn], cg->fgfs[Axzy->sgfn], cg->fgfs[Axzz->sgfn],
                      cg->X[0], cg->X[1], cg->X[2], 
                      ANT, SYM, ANT, Symmetry, levh);
            f_fderivs(cg->shape, cg->fgfs[Ayy0->sgfn], 
                      cg->fgfs[Ayyx->sgfn], cg->fgfs[Ayyy->sgfn], cg->fgfs[Ayyz->sgfn],
                      cg->X[0], cg->X[1], cg->X[2], 
                      SYM, SYM, SYM, Symmetry, levh);
            f_fderivs(cg->shape, cg->fgfs[Ayz0->sgfn], 
                      cg->fgfs[Ayzx->sgfn], cg->fgfs[Ayzy->sgfn], cg->fgfs[Ayzz->sgfn],
                      cg->X[0], cg->X[1], cg->X[2], 
                      SYM, ANT, ANT, Symmetry, levh);
            f_fderivs(cg->shape, cg->fgfs[Azz0->sgfn], 
                      cg->fgfs[Azzx->sgfn], cg->fgfs[Azzy->sgfn], cg->fgfs[Azzz->sgfn],
                      cg->X[0], cg->X[1], cg->X[2], 
                      SYM, SYM, SYM, Symmetry, levh);
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
              f_fderivs_shc(cg->shape, cg->fgfs[phi0->sgfn], 
                            cg->fgfs[phix->sgfn], cg->fgfs[phiy->sgfn], cg->fgfs[phiz->sgfn],
                            cg->X[0], cg->X[1], cg->X[2],
                            phi0->SoA[0], phi0->SoA[1], phi0->SoA[2], 
                            Symmetry, levh, Pp->data->sst,
                            cg->fgfs[fngfs + ShellPatch::drhodx], 
                            cg->fgfs[fngfs + ShellPatch::drhody], 
                            cg->fgfs[fngfs + ShellPatch::drhodz],
                            cg->fgfs[fngfs + ShellPatch::dsigmadx], 
                            cg->fgfs[fngfs + ShellPatch::dsigmady], 
                            cg->fgfs[fngfs + ShellPatch::dsigmadz],
                            cg->fgfs[fngfs + ShellPatch::dRdx], 
                            cg->fgfs[fngfs + ShellPatch::dRdy], 
                            cg->fgfs[fngfs + ShellPatch::dRdz]);
              f_fderivs_shc(cg->shape, cg->fgfs[trK0->sgfn], 
                            cg->fgfs[trKx->sgfn], cg->fgfs[trKy->sgfn], cg->fgfs[trKz->sgfn],
                            cg->X[0], cg->X[1], cg->X[2],
                            trK0->SoA[0], trK0->SoA[1], trK0->SoA[2], 
                            Symmetry, levh, Pp->data->sst,
                            cg->fgfs[fngfs + ShellPatch::drhodx], 
                            cg->fgfs[fngfs + ShellPatch::drhody], 
                            cg->fgfs[fngfs + ShellPatch::drhodz],
                            cg->fgfs[fngfs + ShellPatch::dsigmadx], 
                            cg->fgfs[fngfs + ShellPatch::dsigmady], 
                            cg->fgfs[fngfs + ShellPatch::dsigmadz],
                            cg->fgfs[fngfs + ShellPatch::dRdx], 
                            cg->fgfs[fngfs + ShellPatch::dRdy], 
                            cg->fgfs[fngfs + ShellPatch::dRdz]);
              f_fderivs_shc(cg->shape, cg->fgfs[Axx0->sgfn], 
                            cg->fgfs[Axxx->sgfn], cg->fgfs[Axxy->sgfn], cg->fgfs[Axxz->sgfn],
                            cg->X[0], cg->X[1], cg->X[2],
                            Axx0->SoA[0], Axx0->SoA[1], Axx0->SoA[2], 
                            Symmetry, levh, Pp->data->sst,
                            cg->fgfs[fngfs + ShellPatch::drhodx], 
                            cg->fgfs[fngfs + ShellPatch::drhody], 
                            cg->fgfs[fngfs + ShellPatch::drhodz],
                            cg->fgfs[fngfs + ShellPatch::dsigmadx], 
                            cg->fgfs[fngfs + ShellPatch::dsigmady], 
                            cg->fgfs[fngfs + ShellPatch::dsigmadz],
                            cg->fgfs[fngfs + ShellPatch::dRdx], 
                            cg->fgfs[fngfs + ShellPatch::dRdy], 
                            cg->fgfs[fngfs + ShellPatch::dRdz]);
              f_fderivs_shc(cg->shape, cg->fgfs[Axy0->sgfn], 
                            cg->fgfs[Axyx->sgfn], cg->fgfs[Axyy->sgfn], cg->fgfs[Axyz->sgfn],
                            cg->X[0], cg->X[1], cg->X[2],
                            Axy0->SoA[0], Axy0->SoA[1], Axy0->SoA[2], 
                            Symmetry, levh, Pp->data->sst,
                            cg->fgfs[fngfs + ShellPatch::drhodx], 
                            cg->fgfs[fngfs + ShellPatch::drhody], 
                            cg->fgfs[fngfs + ShellPatch::drhodz],
                            cg->fgfs[fngfs + ShellPatch::dsigmadx], 
                            cg->fgfs[fngfs + ShellPatch::dsigmady], 
                            cg->fgfs[fngfs + ShellPatch::dsigmadz],
                            cg->fgfs[fngfs + ShellPatch::dRdx], 
                            cg->fgfs[fngfs + ShellPatch::dRdy], 
                            cg->fgfs[fngfs + ShellPatch::dRdz]);
              f_fderivs_shc(cg->shape, cg->fgfs[Axz0->sgfn], 
                            cg->fgfs[Axzx->sgfn], cg->fgfs[Axzy->sgfn], cg->fgfs[Axzz->sgfn],
                            cg->X[0], cg->X[1], cg->X[2],
                            Axz0->SoA[0], Axz0->SoA[1], Axz0->SoA[2], 
                            Symmetry, levh, Pp->data->sst,
                            cg->fgfs[fngfs + ShellPatch::drhodx], 
                            cg->fgfs[fngfs + ShellPatch::drhody], 
                            cg->fgfs[fngfs + ShellPatch::drhodz],
                            cg->fgfs[fngfs + ShellPatch::dsigmadx], 
                            cg->fgfs[fngfs + ShellPatch::dsigmady], 
                            cg->fgfs[fngfs + ShellPatch::dsigmadz],
                            cg->fgfs[fngfs + ShellPatch::dRdx], 
                            cg->fgfs[fngfs + ShellPatch::dRdy], 
                            cg->fgfs[fngfs + ShellPatch::dRdz]);
              f_fderivs_shc(cg->shape, cg->fgfs[Ayy0->sgfn], 
                            cg->fgfs[Ayyx->sgfn], cg->fgfs[Ayyy->sgfn], cg->fgfs[Ayyz->sgfn],
                            cg->X[0], cg->X[1], cg->X[2],
                            Ayy0->SoA[0], Ayy0->SoA[1], Ayy0->SoA[2], 
                            Symmetry, levh, Pp->data->sst,
                            cg->fgfs[fngfs + ShellPatch::drhodx], 
                            cg->fgfs[fngfs + ShellPatch::drhody], 
                            cg->fgfs[fngfs + ShellPatch::drhodz],
                            cg->fgfs[fngfs + ShellPatch::dsigmadx], 
                            cg->fgfs[fngfs + ShellPatch::dsigmady], 
                            cg->fgfs[fngfs + ShellPatch::dsigmadz],
                            cg->fgfs[fngfs + ShellPatch::dRdx], 
                            cg->fgfs[fngfs + ShellPatch::dRdy], 
                            cg->fgfs[fngfs + ShellPatch::dRdz]);
              f_fderivs_shc(cg->shape, cg->fgfs[Ayz0->sgfn], 
                            cg->fgfs[Ayzx->sgfn], cg->fgfs[Ayzy->sgfn], cg->fgfs[Ayzz->sgfn],
                            cg->X[0], cg->X[1], cg->X[2],
                            Ayz0->SoA[0], Ayz0->SoA[1], Ayz0->SoA[2], 
                            Symmetry, levh, Pp->data->sst,
                            cg->fgfs[fngfs + ShellPatch::drhodx], 
                            cg->fgfs[fngfs + ShellPatch::drhody], 
                            cg->fgfs[fngfs + ShellPatch::drhodz],
                            cg->fgfs[fngfs + ShellPatch::dsigmadx], 
                            cg->fgfs[fngfs + ShellPatch::dsigmady], 
                            cg->fgfs[fngfs + ShellPatch::dsigmadz],
                            cg->fgfs[fngfs + ShellPatch::dRdx], 
                            cg->fgfs[fngfs + ShellPatch::dRdy], 
                            cg->fgfs[fngfs + ShellPatch::dRdz]);
              f_fderivs_shc(cg->shape, cg->fgfs[Azz0->sgfn], 
                            cg->fgfs[Azzx->sgfn], cg->fgfs[Azzy->sgfn], cg->fgfs[Azzz->sgfn],
                            cg->X[0], cg->X[1], cg->X[2],
                            Azz0->SoA[0], Azz0->SoA[1], Azz0->SoA[2], 
                            Symmetry, levh, Pp->data->sst,
                            cg->fgfs[fngfs + ShellPatch::drhodx], 
                            cg->fgfs[fngfs + ShellPatch::drhody], 
                            cg->fgfs[fngfs + ShellPatch::drhodz],
                            cg->fgfs[fngfs + ShellPatch::dsigmadx], 
                            cg->fgfs[fngfs + ShellPatch::dsigmady], 
                            cg->fgfs[fngfs + ShellPatch::dsigmadz],
                            cg->fgfs[fngfs + ShellPatch::dRdx], 
                            cg->fgfs[fngfs + ShellPatch::dRdy], 
                            cg->fgfs[fngfs + ShellPatch::dRdz]);
            }
            if (BL == Pp->data->ble)
              break;
            BL = BL->next;
          }
          Pp = Pp->next;
        }
      }
#endif
    }
#else
    Compute_Psi4(lev);
#endif
    double *RP, *IP, *RoutMAP;
    int NN = 0;
    for (int pl = 2; pl < maxl + 1; pl++)
      for (int pm = -pl; pm < pl + 1; pm++)
        NN++;
    RP = new double[NN];
    IP = new double[NN];
    RoutMAP = new double[7];
    double Rex = maxrex;
    for (int i = 0; i < decn; i++)
    {
#ifdef Point_Psi4
      Waveshell->surf_Wave(Rex, GH, SH,
                           phi, trK,
                           gxx0, gxy0, gxz0, gyy0, gyz0, gzz0,
                           Axx0, Axy0, Axz0, Ayy0, Ayz0, Azz0,
                           phix, phiy, phiz,
                           trKx, trKy, trKz,
                           Axxx, Axxy, Axxz,
                           Axyx, Axyy, Axyz,
                           Axzx, Axzy, Axzz,
                           Ayyx, Ayyy, Ayyz,
                           Ayzx, Ayzy, Ayzz,
                           Azzx, Azzy, Azzz,
                           Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
                           Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
                           Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
                           Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
                           2, maxl, NN, RP, IP, ErrorMonitor);
#ifdef WithShell
      if (lev > 0 || Rex < GH->bbox[0][0][3])
      {
        Waveshell->surf_MassPAng(Rex, lev, GH, phi0, trK0,
                                 gxx0, gxy0, gxz0, gyy0, gyz0, gzz0,
                                 Axx0, Axy0, Axz0, Ayy0, Ayz0, Azz0,
                                 Gmx0, Gmy0, Gmz0, Sfx1, Sfy1, Sfz1, // here we can not touch rhs variables, but 1 variables
                                 RoutMAP, ErrorMonitor);
      }
      else
      {
        Waveshell->surf_MassPAng(Rex, lev, SH, phi0, trK0,
                                 gxx0, gxy0, gxz0, gyy0, gyz0, gzz0,
                                 Axx0, Axy0, Axz0, Ayy0, Ayz0, Azz0,
                                 Gmx0, Gmy0, Gmz0, Sfx1, Sfy1, Sfz1, // here we can not touch rhs variables, but 1 variables
                                 RoutMAP, ErrorMonitor);
      }
#else
      Waveshell->surf_MassPAng(Rex, lev, GH, phi0, trK0,
                               gxx0, gxy0, gxz0, gyy0, gyz0, gzz0,
                               Axx0, Axy0, Axz0, Ayy0, Ayz0, Azz0,
                               Gmx0, Gmy0, Gmz0, Sfx1, Sfy1, Sfz1, // here we can not touch rhs variables, but 1 variables
                               RoutMAP, ErrorMonitor);
#endif
#else
//        misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],"before surface integral");
#ifdef WithShell
      if (lev > 0 || Rex < GH->bbox[0][0][3])
      {
        Waveshell->surf_Wave(Rex, lev, GH, Rpsi4, Ipsi4, 2, maxl, NN, RP, IP, ErrorMonitor);
        Waveshell->surf_MassPAng(Rex, lev, GH, phi0, trK0,
                                 gxx0, gxy0, gxz0, gyy0, gyz0, gzz0,
                                 Axx0, Axy0, Axz0, Ayy0, Ayz0, Azz0,
                                 Gmx0, Gmy0, Gmz0, Sfx1, Sfy1, Sfz1, // here we can not touch rhs variables, but 1 variables
                                 RoutMAP, ErrorMonitor);
      }
      else
      {
        Waveshell->surf_Wave(Rex, lev, SH, Rpsi4, Ipsi4, 2, maxl, NN, RP, IP, ErrorMonitor);
        Waveshell->surf_MassPAng(Rex, lev, SH, phi0, trK0,
                                 gxx0, gxy0, gxz0, gyy0, gyz0, gzz0,
                                 Axx0, Axy0, Axz0, Ayy0, Ayz0, Azz0,
                                 Gmx0, Gmy0, Gmz0, Sfx1, Sfy1, Sfz1, // here we can not touch rhs variables, but 1 variables
                                 RoutMAP, ErrorMonitor);
      }
#else
#if (PSTR == 0)
      Waveshell->surf_Wave(Rex, lev, GH, Rpsi4, Ipsi4, 2, maxl, NN, RP, IP, ErrorMonitor);
      Waveshell->surf_MassPAng(Rex, lev, GH, phi0, trK0,
                               gxx0, gxy0, gxz0, gyy0, gyz0, gzz0,
                               Axx0, Axy0, Axz0, Ayy0, Ayz0, Azz0,
                               Gmx0, Gmy0, Gmz0, Sfx1, Sfy1, Sfz1, // here we can not touch rhs variables, but 1 variables
                               RoutMAP, ErrorMonitor);
#elif (PSTR == 1 || PSTR == 2)
      Waveshell->surf_Wave(Rex, lev, GH, Rpsi4, Ipsi4, 2, maxl, NN, RP, IP, ErrorMonitor, GH->Commlev[lev]);
      //        misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],"after surf_Wave");
      Waveshell->surf_MassPAng(Rex, lev, GH, phi0, trK0,
                               gxx0, gxy0, gxz0, gyy0, gyz0, gzz0,
                               Axx0, Axy0, Axz0, Ayy0, Ayz0, Azz0,
                               Gmx0, Gmy0, Gmz0, Sfx1, Sfy1, Sfz1, // here we can not touch rhs variables, but 1 variables
                               RoutMAP, ErrorMonitor, GH->Commlev[lev]);
#endif
#endif
//        misc::tillherecheck(GH->Commlev[lev],GH->start_rank[lev],"end surface integral");
#endif
      if (i == 0)
      {
        ADMMass = RoutMAP[0];
      }
#if (PSTR == 1 || PSTR == 2)
      if (GH->start_rank[a_lev] > 0)
      {
        MPI_Status status;
        // receive
        if (myrank == 0)
        {
          MPI_Recv(RP, NN, MPI_DOUBLE, GH->start_rank[a_lev], 1, MPI_COMM_WORLD, &status);
          MPI_Recv(IP, NN, MPI_DOUBLE, GH->start_rank[a_lev], 2, MPI_COMM_WORLD, &status);
          MPI_Recv(RoutMAP, 7, MPI_DOUBLE, GH->start_rank[a_lev], 3, MPI_COMM_WORLD, &status);
        }
        // send
        if (myrank == GH->start_rank[a_lev])
        {
          MPI_Send(RP, NN, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
          MPI_Send(IP, NN, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
          MPI_Send(RoutMAP, 7, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
        }
      }
#endif
      Psi4Monitor->writefile(PhysTime, NN, RP, IP);
      MAPMonitor->writefile(PhysTime, 7, RoutMAP);
      Rex = Rex - drex;
    }
    delete[] RP;
    delete[] IP;
    delete[] RoutMAP;

    // black hole's position
    {
      double *pox;
      pox = new double[dim * BH_num];
      for (int bhi = 0; bhi < BH_num; bhi++)
        for (int i = 0; i < dim; i++)
          pox[dim * bhi + i] = Porg0[bhi][i];
      BHMonitor->writefile(PhysTime, dim * BH_num, pox);
      delete[] pox;
    }

    LastAnas = 0;
  }
}

//================================================================================================



//================================================================================================

// 该成员函数用于计算和输出约束违反

//================================================================================================

void bssn_class::Constraint_Out()
{
  LastConsOut += dT * pow(0.5, Mymax(0, trfls));

  if (LastConsOut >= AnasTime)
  // Constraint violation
  {
    // recompute least the constraint data lost for moved new grid
    for (int lev = 0; lev < GH->levels; lev++)
    {
      // make sure the data consistent for higher levels
      if (lev > 0) // if the constrait quantities can be reused from the step rhs calculation
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
              f_compute_rhs_bssn(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                                 cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                                 cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                                 cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                                 cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                                 cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                                 cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                                 cg->fgfs[Lap0->sgfn], 
                                 cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                                 cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn],
                                 cg->fgfs[phi_rhs->sgfn], cg->fgfs[trK_rhs->sgfn],
                                 cg->fgfs[gxx_rhs->sgfn], cg->fgfs[gxy_rhs->sgfn], cg->fgfs[gxz_rhs->sgfn],
                                 cg->fgfs[gyy_rhs->sgfn], cg->fgfs[gyz_rhs->sgfn], cg->fgfs[gzz_rhs->sgfn],
                                 cg->fgfs[Axx_rhs->sgfn], cg->fgfs[Axy_rhs->sgfn], cg->fgfs[Axz_rhs->sgfn],
                                 cg->fgfs[Ayy_rhs->sgfn], cg->fgfs[Ayz_rhs->sgfn], cg->fgfs[Azz_rhs->sgfn],
                                 cg->fgfs[Gmx_rhs->sgfn], cg->fgfs[Gmy_rhs->sgfn], cg->fgfs[Gmz_rhs->sgfn],
                                 cg->fgfs[Lap_rhs->sgfn], 
                                 cg->fgfs[Sfx_rhs->sgfn], cg->fgfs[Sfy_rhs->sgfn], cg->fgfs[Sfz_rhs->sgfn],
                                 cg->fgfs[dtSfx_rhs->sgfn], cg->fgfs[dtSfy_rhs->sgfn], cg->fgfs[dtSfz_rhs->sgfn],
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
    if (0) // if the constrait quantities can be reused from the step rhs calculation
    {
      MyList<ss_patch> *sPp;
      sPp = SH->PatL;
      while (sPp)
      {
        double TRK4 = PhysTime;
        int pre = 0;
        int lev = 0;
        MyList<Block> *BP = sPp->data->blb;
        int fngfs = sPp->data->fngfs;
        while (BP)
        {
          Block *cg = BP->data;
          if (myrank == cg->rank)
          {
            f_compute_rhs_bssn_ss(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
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
                                  cg->fgfs[phi_rhs->sgfn], cg->fgfs[trK_rhs->sgfn],
                                  cg->fgfs[gxx_rhs->sgfn], cg->fgfs[gxy_rhs->sgfn], cg->fgfs[gxz_rhs->sgfn],
                                  cg->fgfs[gyy_rhs->sgfn], cg->fgfs[gyz_rhs->sgfn], cg->fgfs[gzz_rhs->sgfn],
                                  cg->fgfs[Axx_rhs->sgfn], cg->fgfs[Axy_rhs->sgfn], cg->fgfs[Axz_rhs->sgfn],
                                  cg->fgfs[Ayy_rhs->sgfn], cg->fgfs[Ayz_rhs->sgfn], cg->fgfs[Azz_rhs->sgfn],
                                  cg->fgfs[Gmx_rhs->sgfn], cg->fgfs[Gmy_rhs->sgfn], cg->fgfs[Gmz_rhs->sgfn],
                                  cg->fgfs[Lap_rhs->sgfn], 
                                  cg->fgfs[Sfx_rhs->sgfn], cg->fgfs[Sfy_rhs->sgfn], cg->fgfs[Sfz_rhs->sgfn],
                                  cg->fgfs[dtSfx_rhs->sgfn], cg->fgfs[dtSfy_rhs->sgfn], cg->fgfs[dtSfz_rhs->sgfn],
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
                                  Symmetry, lev, numepsh, sPp->data->sst, pre);
          }
          if (BP == sPp->data->ble)
            break;
          BP = BP->next;
        }
        sPp = sPp->next;
      }
    }
    SH->Synch(ConstraintList, Symmetry);
#endif

    double ConV[7];
#if (PSTR == 1 || PSTR == 2)
    double ConV_h[7];
#endif

#ifdef WithShell
    ConV[0] = SH->L2Norm(Cons_Ham);
    ConV[1] = SH->L2Norm(Cons_Px);
    ConV[2] = SH->L2Norm(Cons_Py);
    ConV[3] = SH->L2Norm(Cons_Pz);
    ConV[4] = SH->L2Norm(Cons_Gx);
    ConV[5] = SH->L2Norm(Cons_Gy);
    ConV[6] = SH->L2Norm(Cons_Gz);
    ConVMonitor->writefile(PhysTime, 7, ConV);
#endif
    for (int levi = 0; levi < GH->levels; levi++)
    {
#if (PSTR == 0)
      ConV[0] = Parallel::L2Norm(GH->PatL[levi]->data, Cons_Ham);
      ConV[1] = Parallel::L2Norm(GH->PatL[levi]->data, Cons_Px);
      ConV[2] = Parallel::L2Norm(GH->PatL[levi]->data, Cons_Py);
      ConV[3] = Parallel::L2Norm(GH->PatL[levi]->data, Cons_Pz);
      ConV[4] = Parallel::L2Norm(GH->PatL[levi]->data, Cons_Gx);
      ConV[5] = Parallel::L2Norm(GH->PatL[levi]->data, Cons_Gy);
      ConV[6] = Parallel::L2Norm(GH->PatL[levi]->data, Cons_Gz);
#elif (PSTR == 1 || PSTR == 2)
      ConV[0] = Parallel::L2Norm(GH->PatL[levi]->data, Cons_Ham, GH->Commlev[levi]);
      ConV[1] = Parallel::L2Norm(GH->PatL[levi]->data, Cons_Px, GH->Commlev[levi]);
      ConV[2] = Parallel::L2Norm(GH->PatL[levi]->data, Cons_Py, GH->Commlev[levi]);
      ConV[3] = Parallel::L2Norm(GH->PatL[levi]->data, Cons_Pz, GH->Commlev[levi]);
      ConV[4] = Parallel::L2Norm(GH->PatL[levi]->data, Cons_Gx, GH->Commlev[levi]);
      ConV[5] = Parallel::L2Norm(GH->PatL[levi]->data, Cons_Gy, GH->Commlev[levi]);
      ConV[6] = Parallel::L2Norm(GH->PatL[levi]->data, Cons_Gz, GH->Commlev[levi]);
      //        misc::tillherecheck("before collect data to cpu0");
      // MPI_ALLREDUCE( sendbuf, recvbuf, count, datatype, op, comm), sendbu and recvbuf must be different
      if (levi > 0)
      {
        if (GH->mylev == levi && myrank == GH->start_rank[levi])
          for (int i = 0; i < 7; i++)
            ConV_h[i] = ConV[i];
        else
          for (int i = 0; i < 7; i++)
            ConV_h[i] = 0;
        MPI_Allreduce(ConV_h, ConV, 7, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      }
#endif
      ConVMonitor->writefile(PhysTime, 7, ConV);
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

    Interp_Constraint(false);

    LastConsOut = 0;
  }
}

//================================================================================================



//================================================================================================

// 该成员函数用于计算表观视界需要的导数

//================================================================================================

#ifdef With_AHF
void bssn_class::AH_Prepare_derivatives()
{
  double SYM = 1.0, ANT = -1.0;
  int ZEO = 0;

  for (int lev = 0; lev < GH->levels; lev++)
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
          f_fderivs(cg->shape, cg->fgfs[phi0->sgfn], 
                    cg->fgfs[dtSfx_rhs->sgfn], cg->fgfs[dtSfy_rhs->sgfn], cg->fgfs[dtSfz_rhs->sgfn],
                    cg->X[0], cg->X[1], cg->X[2], 
                    SYM, SYM, SYM, Symmetry, ZEO);
          f_fderivs(cg->shape, cg->fgfs[gxx0->sgfn], 
                    cg->fgfs[Gamxxx->sgfn], cg->fgfs[Gamyxx->sgfn], cg->fgfs[Gamzxx->sgfn],
                    cg->X[0], cg->X[1], cg->X[2], 
                    SYM, SYM, SYM, Symmetry, ZEO);
          f_fderivs(cg->shape, cg->fgfs[gxy0->sgfn], 
                    cg->fgfs[Gamxxy->sgfn], cg->fgfs[Gamyxy->sgfn], cg->fgfs[Gamzxy->sgfn],
                    cg->X[0], cg->X[1], cg->X[2], 
                    ANT, ANT, SYM, Symmetry, ZEO);
          f_fderivs(cg->shape, cg->fgfs[gxz0->sgfn], 
                    cg->fgfs[Gamxxz->sgfn], cg->fgfs[Gamyxz->sgfn], cg->fgfs[Gamzxz->sgfn],
                    cg->X[0], cg->X[1], cg->X[2], 
                    ANT, SYM, ANT, Symmetry, ZEO);
          f_fderivs(cg->shape, cg->fgfs[gyy0->sgfn], 
                    cg->fgfs[Gamxyy->sgfn], cg->fgfs[Gamyyy->sgfn], cg->fgfs[Gamzyy->sgfn],
                    cg->X[0], cg->X[1], cg->X[2], 
                    SYM, SYM, SYM, Symmetry, ZEO);
          f_fderivs(cg->shape, cg->fgfs[gyz0->sgfn], 
                    cg->fgfs[Gamxyz->sgfn], cg->fgfs[Gamyyz->sgfn], cg->fgfs[Gamzyz->sgfn],
                    cg->X[0], cg->X[1], cg->X[2], 
                    SYM, ANT, ANT, Symmetry, ZEO);
          f_fderivs(cg->shape, cg->fgfs[gzz0->sgfn], 
                    cg->fgfs[Gamxzz->sgfn], cg->fgfs[Gamyzz->sgfn], cg->fgfs[Gamzzz->sgfn],
                    cg->X[0], cg->X[1], cg->X[2], 
                    SYM, SYM, SYM, Symmetry, ZEO);
        }
        if (BP == Pp->data->ble)
          break;
        BP = BP->next;
      }
      Pp = Pp->next;
    }
    Parallel::Sync(GH->PatL[lev], AHDList, Symmetry);
  }
}

//================================================================================================



//================================================================================================

// 该成员函数用于对表观视界数据进行插值

//================================================================================================

bool bssn_class::AH_Interp_Points(MyList<var> *VarList,
                                  int NN, double **XX,
                                  double *Shellf, int Symmetryi)
{
  MyList<var> *varl;
  int num_var = 0;
  varl = VarList;
  while (varl)
  {
    num_var++;
    varl = varl->next;
  }

  double pox[3];
  for (int i = 0; i < NN; i++)
  {
    for (int j = 0; j < 3; j++)
      pox[j] = XX[j][i];
    int lev = GH->levels - 1;
    bool notfound = true;

    while (notfound)
    {
      if (lev < 0)
      {
#ifdef WithShell
        if (SH->Interp_One_Point(VarList, pox, Shellf + i * num_var, Symmetryi))
        {
          return true;
        }
        if (myrank == 0)
        {
          cout << " bssn_class::AH_Interp_Points: point (" 
               << pox[0] << "," << pox[1] << "," << pox[2] 
               << ") is out of cgh and shell domain!" << endl;
          if (ErrorMonitor->outfile)
            ErrorMonitor->outfile << " bssn_class::AH_Interp_Points: point (" 
                                  << pox[0] << "," << pox[1] << "," << pox[2] 
                                  << ") is out of cgh and shell domain!" << endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
#else
        if (myrank == 0)
        {
          cout << " bssn_class::AH_Interp_Points: point (" 
               << pox[0] << "," << pox[1] << "," << pox[2] 
               << ") is out of cgh domain!" << endl;
          if (ErrorMonitor->outfile)
            ErrorMonitor->outfile << " bssn_class::AH_Interp_Points: point (" 
                                  << pox[0] << "," << pox[1] << "," << pox[2] 
                                  << ") is out of cgh domain!" << endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
#endif
        return false;
      }
      MyList<Patch> *Pp = GH->PatL[lev];
      while (Pp)
      {
        if (Pp->data->Interp_ONE_Point(VarList, pox, Shellf + i * num_var, Symmetryi))
        {
          notfound = false;
          break;
        }
        Pp = Pp->next;
      }
      lev--;
    }
  }
  return true;
}

//================================================================================================



//================================================================================================

// 该成员函数用于计算表观视界

//================================================================================================

void bssn_class::AH_Step_Find(int lev, double dT_lev)
{
  if ((lev == GH->levels - 1))
  {
    int ncount = int(PhysTime / dT_lev);
    bool tf = false;
    for (int ihn = 0; ihn < HN_num; ihn++)
    {
      if (ncount % findeveryl[ihn] == 0)
      {
        tf = true;
        break;
      }
    }
    if (tf)
    {
      clock_t prev_clock, curr_clock;
      if (myrank == 0)
        prev_clock = clock();
      const int cdumpid = int(PhysTime / AHdumptime) + 1;
      for (int ihn = 0; ihn < HN_num; ihn++)
        dumpid[ihn] = cdumpid;

      double gam;
      for (int ihn = 0; ihn < BH_num; ihn++)
      {
        xc[ihn] = Porg0[ihn][0];
        yc[ihn] = Porg0[ihn][1];
        zc[ihn] = Porg0[ihn][2];
        gam = fabs(Pmom[ihn * 3]) / (Mass[ihn]);
        gam = sqrt(1 - gam * gam);
        xr[ihn] = Mass[ihn] * gam;
        gam = fabs(Pmom[ihn * 3 + 1]) / (Mass[ihn]);
        gam = sqrt(1 - gam * gam);
        yr[ihn] = Mass[ihn] * gam;
        gam = fabs(Pmom[ihn * 3 + 2]) / (Mass[ihn]);
        gam = sqrt(1 - gam * gam);
        zr[ihn] = Mass[ihn] * gam;
        dTT[ihn] = -1;

        if (ncount % findeveryl[ihn] == 0)
        {
          trigger[ihn] = true;
          dTT[ihn] = findeveryl[ihn] * dT_lev;
        }
        else
          trigger[ihn] = false;
        if (trigger[ihn] && (dumpid[ihn] > lastahdumpid[ihn]))
          lastahdumpid[ihn] = dumpid[ihn];
        else
          dumpid[ihn] = 0;
      }
      int ihn = BH_num;
      for (int ia = 0; ia < BH_num; ia++)
        for (int ib = ia + 1; ib < BH_num; ib++)
        {
          xc[ihn] = (Porg0[ia][0] + Porg0[ib][0]) / 2;
          yc[ihn] = (Porg0[ia][1] + Porg0[ib][1]) / 2;
          zc[ihn] = (Porg0[ia][2] + Porg0[ib][2]) / 2;

          xr[ihn] = yr[ihn] = zr[ihn] = Mass[ia] + Mass[ib];

          dTT[ihn] = -1;

          if (fabs(Porg0[ia][0] - Porg0[ib][0]) < 2 * xr[ihn] &&
              fabs(Porg0[ia][1] - Porg0[ib][1]) < 2 * xr[ihn] &&
              fabs(Porg0[ia][2] - Porg0[ib][2]) < 2 * xr[ihn] &&
              (ncount % findeveryl[ihn] == 0))
          {
            trigger[ihn] = true;
            dTT[ihn] = findeveryl[ihn] * dT_lev;
          }
          else
            trigger[ihn] = false;

          if (trigger[ihn] && (dumpid[ihn] > lastahdumpid[ihn]))
            lastahdumpid[ihn] = dumpid[ihn];
          else
            dumpid[ihn] = 0;

          ihn++;
        }
#if (ABEtype == 1)
      if (PhysTime > 10)
      {
        ihn--;
        trigger[ihn] = true;
        xr[ihn] = yr[ihn] = zr[ihn] = 50;
        //	if(myrank==0) for(ihn=0;ihn<HN_num;ihn++) cout<<"trigger#"<<ihn<<": "<<trigger[ihn]<<endl;
      }
#endif
      AHFinderDirect::AHFinderDirect_find_horizons(HN_num, dumpid,
                                                   xc, yc, zc, xr, yr, zr, trigger, dTT); 
      // note rhs and Gamijk have been used as temp storage space

      if (myrank == 0)
      {
        curr_clock = clock();
        cout << " Finding horizon used " 
             << (double)(curr_clock - prev_clock) / ((double)CLOCKS_PER_SEC) 
             << " seconds! " << endl;
      }
    }
  }
}
#endif

//================================================================================================



//================================================================================================

// 该成员函数用于对约束数据进行插值

//================================================================================================

void bssn_class::Interp_Constraint(bool infg)
{
  if (infg)
  {
    // we do not support a_lev != 0 yet.
    if (a_lev > 0)
      return;

    // recompute least the constraint data lost for moved new grid
    for (int lev = 0; lev < GH->levels; lev++)
    {
      // make sure the data consistent for higher levels
      if (lev > 0) // if the constrait quantities can be reused from the step rhs calculation
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
              f_compute_rhs_bssn(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                                 cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                                 cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                                 cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                                 cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                                 cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                                 cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                                 cg->fgfs[Lap0->sgfn], 
                                 cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                                 cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn],
                                 cg->fgfs[phi_rhs->sgfn], cg->fgfs[trK_rhs->sgfn],
                                 cg->fgfs[gxx_rhs->sgfn], cg->fgfs[gxy_rhs->sgfn], cg->fgfs[gxz_rhs->sgfn],
                                 cg->fgfs[gyy_rhs->sgfn], cg->fgfs[gyz_rhs->sgfn], cg->fgfs[gzz_rhs->sgfn],
                                 cg->fgfs[Axx_rhs->sgfn], cg->fgfs[Axy_rhs->sgfn], cg->fgfs[Axz_rhs->sgfn],
                                 cg->fgfs[Ayy_rhs->sgfn], cg->fgfs[Ayz_rhs->sgfn], cg->fgfs[Azz_rhs->sgfn],
                                 cg->fgfs[Gmx_rhs->sgfn], cg->fgfs[Gmy_rhs->sgfn], cg->fgfs[Gmz_rhs->sgfn],
                                 cg->fgfs[Lap_rhs->sgfn], 
                                 cg->fgfs[Sfx_rhs->sgfn], cg->fgfs[Sfy_rhs->sgfn], cg->fgfs[Sfz_rhs->sgfn],
                                 cg->fgfs[dtSfx_rhs->sgfn], cg->fgfs[dtSfy_rhs->sgfn], cg->fgfs[dtSfz_rhs->sgfn],
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
    if (0) // if the constrait quantities can be reused from the step rhs calculation
    {
      MyList<ss_patch> *sPp;
      sPp = SH->PatL;
      while (sPp)
      {
        double TRK4 = PhysTime;
        int pre = 0;
        int lev = 0;
        MyList<Block> *BP = sPp->data->blb;
        int fngfs = sPp->data->fngfs;
        while (BP)
        {
          Block *cg = BP->data;
          if (myrank == cg->rank)
          {
            f_compute_rhs_bssn_ss(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
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
                                  cg->fgfs[phi_rhs->sgfn], cg->fgfs[trK_rhs->sgfn],
                                  cg->fgfs[gxx_rhs->sgfn], cg->fgfs[gxy_rhs->sgfn], cg->fgfs[gxz_rhs->sgfn],
                                  cg->fgfs[gyy_rhs->sgfn], cg->fgfs[gyz_rhs->sgfn], cg->fgfs[gzz_rhs->sgfn],
                                  cg->fgfs[Axx_rhs->sgfn], cg->fgfs[Axy_rhs->sgfn], cg->fgfs[Axz_rhs->sgfn],
                                  cg->fgfs[Ayy_rhs->sgfn], cg->fgfs[Ayz_rhs->sgfn], cg->fgfs[Azz_rhs->sgfn],
                                  cg->fgfs[Gmx_rhs->sgfn], cg->fgfs[Gmy_rhs->sgfn], cg->fgfs[Gmz_rhs->sgfn],
                                  cg->fgfs[Lap_rhs->sgfn], 
                                  cg->fgfs[Sfx_rhs->sgfn], cg->fgfs[Sfy_rhs->sgfn], cg->fgfs[Sfz_rhs->sgfn],
                                  cg->fgfs[dtSfx_rhs->sgfn], cg->fgfs[dtSfy_rhs->sgfn], cg->fgfs[dtSfz_rhs->sgfn],
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
                                  Symmetry, lev, numepsh, sPp->data->sst, pre);
          }
          if (BP == sPp->data->ble)
            break;
          BP = BP->next;
        }
        sPp = sPp->next;
      }
    }
    SH->Synch(ConstraintList, Symmetry);
#endif
  }
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

  if (myrank == 0)
  {
    ofstream outfile;
    char filename[50];
    sprintf(filename, "%s/interp_constraint_%05d.dat", ErrorMonitor->out_dir.c_str(), int(PhysTime / dT + 0.5)); 
    // 0.5 for round off
    
    outfile.open(filename);
    outfile << "#  corrdinate, H_Res, Px_Res, Py_Res, Pz_Res, Gx_Res, Gy_Res, Gz_Res, ...." << endl;
    for (int i = 0; i < n; i++)
    {
      outfile << setw(10) << setprecision(10) << y1[i];
      for (int j = 0; j < InList; j++)
        outfile << " " << setw(16) << setprecision(15) << shellf[InList * i + j];
      outfile << endl;
    }
    outfile.close();
  }

  delete[] shellf;
}

//================================================================================================



//================================================================================================

// 该成员函数用于计算约束违反

//================================================================================================

void bssn_class::Compute_Constraint()
{
  double TRK4 = PhysTime;
  double ndeps = numepsb;
  int pre = 0;
  int lev;

  for (lev = 0; lev < GH->levels; lev++)
  {
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
            f_compute_rhs_bssn(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
                               cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                               cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                               cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                               cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                               cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                               cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                               cg->fgfs[Lap0->sgfn], 
                               cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                               cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn],
                               cg->fgfs[phi_rhs->sgfn], cg->fgfs[trK_rhs->sgfn],
                               cg->fgfs[gxx_rhs->sgfn], cg->fgfs[gxy_rhs->sgfn], cg->fgfs[gxz_rhs->sgfn],
                               cg->fgfs[gyy_rhs->sgfn], cg->fgfs[gyz_rhs->sgfn], cg->fgfs[gzz_rhs->sgfn],
                               cg->fgfs[Axx_rhs->sgfn], cg->fgfs[Axy_rhs->sgfn], cg->fgfs[Axz_rhs->sgfn],
                               cg->fgfs[Ayy_rhs->sgfn], cg->fgfs[Ayz_rhs->sgfn], cg->fgfs[Azz_rhs->sgfn],
                               cg->fgfs[Gmx_rhs->sgfn], cg->fgfs[Gmy_rhs->sgfn], cg->fgfs[Gmz_rhs->sgfn],
                               cg->fgfs[Lap_rhs->sgfn], 
                               cg->fgfs[Sfx_rhs->sgfn], cg->fgfs[Sfy_rhs->sgfn], cg->fgfs[Sfz_rhs->sgfn],
                               cg->fgfs[dtSfx_rhs->sgfn], cg->fgfs[dtSfy_rhs->sgfn], cg->fgfs[dtSfz_rhs->sgfn],
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
  // prolong restrict constraint quantities
  for (lev = GH->levels - 1; lev > 0; lev--)
    RestrictProlong(lev, 1, false, ConstraintList, ConstraintList, ConstraintList);

#ifdef WithShell
  lev = 0;
  {
    MyList<ss_patch> *sPp;
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
          f_compute_rhs_bssn_ss(cg->shape, TRK4, cg->X[0], cg->X[1], cg->X[2],
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
                                cg->fgfs[phi_rhs->sgfn], cg->fgfs[trK_rhs->sgfn],
                                cg->fgfs[gxx_rhs->sgfn], cg->fgfs[gxy_rhs->sgfn], cg->fgfs[gxz_rhs->sgfn],
                                cg->fgfs[gyy_rhs->sgfn], cg->fgfs[gyz_rhs->sgfn], cg->fgfs[gzz_rhs->sgfn],
                                cg->fgfs[Axx_rhs->sgfn], cg->fgfs[Axy_rhs->sgfn], cg->fgfs[Axz_rhs->sgfn],
                                cg->fgfs[Ayy_rhs->sgfn], cg->fgfs[Ayz_rhs->sgfn], cg->fgfs[Azz_rhs->sgfn],
                                cg->fgfs[Gmx_rhs->sgfn], cg->fgfs[Gmy_rhs->sgfn], cg->fgfs[Gmz_rhs->sgfn],
                                cg->fgfs[Lap_rhs->sgfn], 
                                cg->fgfs[Sfx_rhs->sgfn], cg->fgfs[Sfy_rhs->sgfn], cg->fgfs[Sfz_rhs->sgfn],
                                cg->fgfs[dtSfx_rhs->sgfn], cg->fgfs[dtSfy_rhs->sgfn], cg->fgfs[dtSfz_rhs->sgfn],
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
                                Symmetry, lev, numepsh, sPp->data->sst, pre);
        }
        if (BP == sPp->data->ble)
          break;
        BP = BP->next;
      }
      sPp = sPp->next;
    }
  }
  SH->Synch(ConstraintList, Symmetry);
  // interpolate constraint quantities
  SH->CS_Inter(ConstraintList, Symmetry);
#endif
}

//================================================================================================



//================================================================================================

void bssn_class::testRestrict()
{
  MyList<var> *DG_List = new MyList<var>(phi0);
  int lev = 0;
  double ZEO = 0, ONE = 1;
  MyList<Patch> *Pp = GH->PatL[lev];
  while (Pp)
  {
    MyList<Block> *BP = Pp->data->blb;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        f_set_value(cg->shape, cg->fgfs[phi0->sgfn], ZEO);
      }
      if (BP == Pp->data->ble)
        break;
      BP = BP->next;
    }
    Pp = Pp->next;
  }

  lev = 1;
  Pp = GH->PatL[lev];
  while (Pp)
  {
    MyList<Block> *BP = Pp->data->blb;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        f_set_value(cg->shape, cg->fgfs[phi0->sgfn], ONE);
      }
      if (BP == Pp->data->ble)
        break;
      BP = BP->next;
    }
    Pp = Pp->next;
  }

  Parallel::Restrict(GH->PatL[lev - 1], GH->PatL[lev], DG_List, DG_List, Symmetry);
  Parallel::Sync(GH->PatL[lev - 1], DG_List, Symmetry);

  Parallel::Dump_Data(GH->PatL[lev - 1], DG_List, 0, PhysTime, dT);
  Parallel::Dump_Data(GH->PatL[lev], DG_List, 0, PhysTime, dT);

  DG_List->clearList();
  exit(0);
}

//================================================================================================



//================================================================================================

void bssn_class::testOutBd()
{
  MyList<var> *DG_List = new MyList<var>(phi0);
  int lev = 1;
  double ZEO = 0, ONE = 1;
  MyList<Patch> *Pp = GH->PatL[lev];
  while (Pp)
  {
    MyList<Block> *BP = Pp->data->blb;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        f_set_value(cg->shape, cg->fgfs[phi0->sgfn], ZEO);
      }
      if (BP == Pp->data->ble)
        break;
      BP = BP->next;
    }
    Pp = Pp->next;
  }

  lev = 0;
  Pp = GH->PatL[lev];
  while (Pp)
  {
    MyList<Block> *BP = Pp->data->blb;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        f_set_value(cg->shape, cg->fgfs[phi0->sgfn], ONE);
      }
      if (BP == Pp->data->ble)
        break;
      BP = BP->next;
    }
    Pp = Pp->next;
  }

  lev = 1;
  MyList<Patch> *Ppc = GH->PatL[lev - 1];
  while (Ppc)
  {
    Pp = GH->PatL[lev];
    while (Pp)
    {
      Parallel::OutBdLow2Hi(Ppc->data, Pp->data, DG_List, DG_List, Symmetry);
      Pp = Pp->next;
    }
    Ppc = Ppc->next;
  }

  Parallel::Sync(GH->PatL[lev], DG_List, Symmetry);

  Parallel::Dump_Data(GH->PatL[lev], DG_List, 0, PhysTime, dT);
  Parallel::Dump_Data(GH->PatL[lev - 1], DG_List, 0, PhysTime, dT);

  DG_List->clearList();
  exit(0);
}

//================================================================================================



//================================================================================================

// 该成员函数用于检验无迹条件

//================================================================================================

void bssn_class::Enforce_algcon(int lev, int fg)
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
        if (fg == 0)
          f_enforce_ga(cg->shape,
                       cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                       cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                       cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                       cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn]);
        else
          f_enforce_ga(cg->shape,
                       cg->fgfs[gxx->sgfn], cg->fgfs[gxy->sgfn], cg->fgfs[gxz->sgfn], 
                       cg->fgfs[gyy->sgfn], cg->fgfs[gyz->sgfn], cg->fgfs[gzz->sgfn],
                       cg->fgfs[Axx->sgfn], cg->fgfs[Axy->sgfn], cg->fgfs[Axz->sgfn], 
                       cg->fgfs[Ayy->sgfn], cg->fgfs[Ayz->sgfn], cg->fgfs[Azz->sgfn]);
      }
      if (BP == Pp->data->ble)
        break;
      BP = BP->next;
    }
    Pp = Pp->next;
  }

#ifdef WithShell
  if (lev == 0)
  {
    MyList<ss_patch> *sPp = SH->PatL;
    while (sPp)
    {
      MyList<Block> *BP = sPp->data->blb;
      int fngfs = sPp->data->fngfs;
      while (BP)
      {
        Block *cg = BP->data;
        if (myrank == cg->rank)
        {
          if (fg == 0)
            f_enforce_ga(cg->shape,
                         cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], 
                         cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                         cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], 
                         cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn]);
          else
            f_enforce_ga(cg->shape,
                         cg->fgfs[gxx->sgfn], cg->fgfs[gxy->sgfn], cg->fgfs[gxz->sgfn], 
                         cg->fgfs[gyy->sgfn], cg->fgfs[gyz->sgfn], cg->fgfs[gzz->sgfn],
                         cg->fgfs[Axx->sgfn], cg->fgfs[Axy->sgfn], cg->fgfs[Axz->sgfn], 
                         cg->fgfs[Ayy->sgfn], cg->fgfs[Ayz->sgfn], cg->fgfs[Azz->sgfn]);
        }
        if (BP == sPp->data->ble)
          break;
        BP = BP->next;
      }
      sPp = sPp->next;
    }
  }
#endif
}

//================================================================================================



//================================================================================================

// 该成员函数用来监视是否在屏幕输入 abort

//================================================================================================

bool bssn_class::check_Stdin_Abort() 
{

    fd_set readfds;

    struct timeval timeout;

    FD_ZERO(&readfds);
    FD_SET(STDIN_FILENO, &readfds);

    // 设置超时为 0 → 非阻塞检查
    timeout.tv_sec = 0;
    timeout.tv_usec = 0;

    int activity = select(STDIN_FILENO + 1, &readfds, nullptr, nullptr, &timeout);

    if (activity > 0 && FD_ISSET(STDIN_FILENO, &readfds)) {
        string input_abort;
        if (cin >> input_abort) {
            if (input_abort == "stop") {
                return true;
            }
        }
    }

    return false;
}

//================================================================================================

