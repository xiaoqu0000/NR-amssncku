
#ifdef newc
#include <sstream>
#include <cstdio>
#include <cstring>
#include <map>
using namespace std;
#else
#include <stdio.h>
#include <string.h>
#include <map.h>
#endif

#include <time.h>

#include "macrodef.h"
#include "misc.h"
#include "Ansorg.h"
#include "fmisc.h"
#include "Parallel.h"
#include "bssnEM_gpu_class.h"
#include "bssn_rhs.h"
#include "empart.h"
#include "initial_puncture.h"
#include "initial_maxwell.h"
#include "enforce_algebra.h"
#include "rungekutta4_rout.h"
#include "sommerfeld_rout.h"
#include "getnp4.h"
#include "getnpem2.h"
#include "shellfunctions.h"
#include "parameters.h"

#include "gpu_manager.h"
#include "helper.h"

#ifdef With_AHF
#include "derivatives.h"
#include "myglobal.h"
#endif

//================================================================================================
// 构造函数，定义 bssnEM_class
// Constructor
//================================================================================================

bssnEM_class::bssnEM_class(double Couranti, double StartTimei, double TotalTimei,
                           double DumpTimei, double d2DumpTimei, double CheckTimei, double AnasTimei,
                           int Symmetryi, int checkruni, char *checkfilenamei,
                           double numepssi, double numepsbi, double numepshi,
                           int a_levi, int maxli, int decni, double maxrexi, double drexi)
                           : bssn_class(Couranti, StartTimei, TotalTimei,
                                        DumpTimei, d2DumpTimei, CheckTimei, AnasTimei,
                                        Symmetryi, checkruni, checkfilenamei, numepssi, numepsbi, numepshi,
                                        a_levi, maxli, decni, maxrexi, drexi)
{
    {
        char str[50];
        stringstream a_stream;
        a_stream.setf(ios::left);
        a_stream.str("");
        a_stream << setw(15) << "# time";
        for (int pl = 1; pl < maxl + 1; pl++)
            for (int pm = -pl; pm < pl + 1; pm++)
            {
                sprintf(str, "R%02dm%03d", pl, pm);
                a_stream << setw(16) << str;
                sprintf(str, "I%02dm%03d", pl, pm);
                a_stream << setw(16) << str;
            }
        Phi2Monitor = new monitor("bssn_phi2.dat", myrank, a_stream.str());

        a_stream.clear();
        a_stream.str("");
        a_stream << setw(15) << "# time";
        for (int pl = 0; pl < maxl + 1; pl++)
            for (int pm = -pl; pm < pl + 1; pm++)
            {
                sprintf(str, "R%02dm%03d", pl, pm);
                a_stream << setw(16) << str;
                sprintf(str, "I%02dm%03d", pl, pm);
                a_stream << setw(16) << str;
            }
        Phi1Monitor = new monitor("bssn_phi1.dat", myrank, a_stream.str());
    }
}

//================================================================================================



//================================================================================================
// 析构函数，用于删除变量
// Destructor
//================================================================================================

bssnEM_class::~bssnEM_class()
{
    delete Kpsio; delete Kphio; 
    delete Exo; delete Eyo; delete Ezo;
    delete Bxo; delete Byo; delete Bzo;

    delete Kpsi0; delete Kphi0; 
    delete Ex0; delete Ey0; delete Ez0;
    delete Bx0; delete By0; delete Bz0;

    delete Kpsi; delete Kphi;  
    delete Ex; delete Ey; delete Ez;
    delete Bx; delete By; delete Bz;

    delete Kpsi1; delete Kphi1; 
    delete Ex1; delete Ey1; delete Ez1;
    delete Bx1; delete By1; delete Bz1;

    delete Kpsi_rhs; delete Kphi_rhs;
    delete Ex_rhs; delete Ey_rhs; delete Ez_rhs;
    delete Bx_rhs; delete By_rhs; delete Bz_rhs;

    delete qchar; delete Jx; delete Jy; delete Jz;

    delete Rphi2; delete Iphi2;
    delete Rphi1; delete Iphi1;
    
    delete Phi2Monitor;
    delete Phi1Monitor;
}

//================================================================================================



//================================================================================================
// 该成员函数用于对 Class 初始化
// Initialize
//================================================================================================

void bssnEM_class::Initialize()
{
    Exo   = new var("Exo", ngfs++, -1,  1,  1);  Eyo = new var("Eyo", ngfs++,  1, -1,  1); Ezo = new var("Ezo", ngfs++,  1,  1, -1);
    Bxo   = new var("Bxo", ngfs++,  1, -1, -1);  Byo = new var("Byo", ngfs++, -1,  1, -1); Bzo = new var("Bzo", ngfs++, -1, -1,  1);
    Kpsio = new var("Kpsio", ngfs++,  1,  1,  1);
    Kphio = new var("Kphio", ngfs++,  1,  1,  1);

    Ex0   = new var("Ex0", ngfs++, -1,  1,  1);  Ey0 = new var("Ey0", ngfs++,  1, -1,  1); Ez0 = new var("Ez0", ngfs++,  1,  1, -1);
    Bx0   = new var("Bx0", ngfs++,  1, -1, -1);  By0 = new var("By0", ngfs++, -1,  1, -1); Bz0 = new var("Bz0", ngfs++, -1, -1,  1);
    Kpsi0 = new var("Kpsi0", ngfs++,  1,  1,  1);
    Kphi0 = new var("Kphi0", ngfs++,  1,  1,  1);

    Ex = new var("Ex", ngfs++, -1,  1,  1); Ey = new var("Ey", ngfs++,  1, -1,  1); Ez = new var("Ez", ngfs++,  1,  1, -1);
    Bx = new var("Bx", ngfs++,  1, -1, -1); By = new var("By", ngfs++, -1,  1, -1); Bz = new var("Bz", ngfs++, -1, -1,  1);
    Kpsi  = new var("Kpsi", ngfs++,  1,  1,  1);
    Kphi  = new var("Kphi", ngfs++,  1,  1,  1);

    Ex1 = new var("Ex1", ngfs++, -1,  1,  1); Ey1 = new var("Ey1", ngfs++,  1, -1,  1); Ez1 = new var("Ez1", ngfs++,  1,  1, -1);
    Bx1 = new var("Bx1", ngfs++,  1, -1, -1); By1 = new var("By1", ngfs++, -1,  1, -1); Bz1 = new var("Bz1", ngfs++, -1, -1,  1);
    Kpsi1 = new var("Kpsi1", ngfs++,  1,  1,  1);
    Kphi1 = new var("Kphi1", ngfs++,  1,  1,  1);

    Ex_rhs = new var("Ex_rhs", ngfs++, -1,  1,  1); Ey_rhs = new var("Ey_rhs", ngfs++,  1, -1,  1); Ez_rhs = new var("Ez_rhs", ngfs++,  1,  1, -1);
    Bx_rhs = new var("Bx_rhs", ngfs++,  1, -1, -1); By_rhs = new var("By_rhs", ngfs++, -1,  1, -1); Bz_rhs = new var("Bz_rhs", ngfs++, -1, -1,  1);
    Kpsi_rhs = new var("Kpsi_rhs", ngfs++,  1,  1,  1);
    Kphi_rhs = new var("Kphi_rhs", ngfs++,  1,  1,  1);

    qchar = new var("qchar", ngfs++,  1,  1,  1);
    Jx    = new var("Jx",    ngfs++, -1,  1,  1); Jy = new var("Jy", ngfs++,  1, -1,  1); Jz = new var("Jz", ngfs++,  1,  1, -1);

    Rphi2 = new var("Rphi2", ngfs++,  1,  1,  1); Iphi2 = new var("Iphi2", ngfs++, -1, -1, -1);
    Rphi1 = new var("Rphi1", ngfs++,  1,  1,  1); Iphi1 = new var("Iphi1", ngfs++,  1,  1,  1);

    if (myrank == 0)
        cout << "you have setted " << ngfs << " grid functions." << endl;

    OldStateList->insert(Kpsio); OldStateList->insert(Kphio);
    OldStateList->insert(Exo); OldStateList->insert(Eyo); OldStateList->insert(Ezo);
    OldStateList->insert(Bxo); OldStateList->insert(Byo); OldStateList->insert(Bzo);

    StateList->insert(Kpsi0); StateList->insert(Kphi0);
    StateList->insert(Ex0); StateList->insert(Ey0); StateList->insert(Ez0);
    StateList->insert(Bx0); StateList->insert(By0); StateList->insert(Bz0);

    RHSList->insert(Kpsi_rhs); RHSList->insert(Kphi_rhs);
    RHSList->insert(Ex_rhs); RHSList->insert(Ey_rhs); RHSList->insert(Ez_rhs);
    RHSList->insert(Bx_rhs); RHSList->insert(By_rhs); RHSList->insert(Bz_rhs);

    SynchList_pre->insert(Kpsi); SynchList_pre->insert(Kphi);
    SynchList_pre->insert(Ex); SynchList_pre->insert(Ey); SynchList_pre->insert(Ez);
    SynchList_pre->insert(Bx); SynchList_pre->insert(By); SynchList_pre->insert(Bz);

    SynchList_cor->insert(Kpsi1); SynchList_cor->insert(Kphi1);
    SynchList_cor->insert(Ex1); SynchList_cor->insert(Ey1); SynchList_cor->insert(Ez1);
    SynchList_cor->insert(Bx1); SynchList_cor->insert(By1); SynchList_cor->insert(Bz1);

    DumpList->insert(Rphi2); DumpList->insert(Iphi2);
    DumpList->insert(Rphi1); DumpList->insert(Iphi1);
    DumpList->insert(Ex0); DumpList->insert(Bx0);

    MiscList->insert(Jx); MiscList->insert(Jy); MiscList->insert(Jz); MiscList->insert(qchar);

    CheckPoint->addvariablelist(StateList);
    CheckPoint->addvariablelist(OldStateList);

    char pname[50];
    {
        map<string, string>::iterator iter = parameters::str_par.find("inputpar");
        if (iter != parameters::str_par.end())
            strcpy(pname, (iter->second).c_str());
        else { cout << "Error inputpar" << endl; exit(0); }
    }
    GH = new cgh(0, ngfs, Symmetry, pname, checkrun, ErrorMonitor);
    if (checkrun) CheckPoint->readcheck_cgh(PhysTime, GH, myrank, nprocs, Symmetry);
    else GH->compose_cgh(nprocs);

#ifdef WithShell
    SH = new ShellPatch(0, ngfs, pname, Symmetry, myrank, ErrorMonitor);
    SH->matchcheck(GH->PatL[0]); SH->compose_sh(nprocs);
    SH->setupcordtrans(); SH->Dump_xyz(0, 0, 1);
    SH->setupintintstuff(nprocs, GH->PatL[0], Symmetry);
#ifdef USE_GPU
SH->build_gpu_block_table(GH->PatL[0]);
SH->build_gpu_pointstru_cache(SH->gpu_ss_cache, SH->ss_src, SH->ss_dst, nprocs, GH->PatL[0]);
SH->build_gpu_pointstru_cache(SH->gpu_csatc_cache, SH->csatc_src, SH->csatc_dst, nprocs, GH->PatL[0]);
SH->build_gpu_pointstru_cache(SH->gpu_csats_cache, SH->csats_src, SH->csats_dst, nprocs, GH->PatL[0]);
#endif
    if (checkrun) CheckPoint->readcheck_sh(SH, myrank);
#endif

    double h = GH->PatL[0]->data->blb->data->getdX(0);
    for (int i = 1; i < dim; i++) h = Mymin(h, GH->PatL[0]->data->blb->data->getdX(i));
    dT = Courant * h;

    if (checkrun)
        CheckPoint->read_Black_Hole_position(BH_num_input, BH_num, Porg0, Pmom, Spin, Mass, Porgbr, Porg, Porg1, Porg_rhs);
    else { PhysTime = StartTime; Setup_Black_Hole_position(); }

    // Zero-initialize EM source terms (external current/charge = 0 for vacuum EM)
    for (int lev = 0; lev < GH->levels; lev++) {
        MyList<Patch> *Pp = GH->PatL[lev];
        while (Pp) {
            MyList<Block> *BL = Pp->data->blb;
            while (BL) {
                Block *cg = BL->data;
                if (myrank == cg->rank) {
                    int n = cg->shape[0] * cg->shape[1] * cg->shape[2];
                    memset(cg->fgfs[Jx->sgfn], 0, n * sizeof(double));
                    memset(cg->fgfs[Jy->sgfn], 0, n * sizeof(double));
                    memset(cg->fgfs[Jz->sgfn], 0, n * sizeof(double));
                    memset(cg->fgfs[qchar->sgfn], 0, n * sizeof(double));
                }
                if (BL == Pp->data->ble) break;
                BL = BL->next;
            }
            Pp = Pp->next;
        }
    }
}

//================================================================================================




//================================================================================================
// 该成员函数读入 Ansorg 方法求解的 TwoPuncture 初值
// Read_Ansorg
//================================================================================================

void bssnEM_class::Read_Ansorg()
{
    if (checkrun) {
        CheckPoint->readcheck_cgh(PhysTime, GH, myrank, nprocs, Symmetry);
#ifdef WithShell
        CheckPoint->readcheck_sh(SH, myrank);
#endif
    } else {
        if (myrank == 0)
            cout << "Read initial data from Ansorg's solver,"
                 << " please be sure the input parameters for black holes are puncture parameters!!" << endl;
        char filename[50];
        { map<string, string>::iterator iter = parameters::str_par.find("inputpar");
          if (iter != parameters::str_par.end()) strcpy(filename, (iter->second).c_str());
          else { cout << "Error inputpar" << endl; exit(0); } }
        int BH_NM; double *Porg_here, *Qchar;
        { const int LEN = 256; char pline[LEN]; string str, sgrp, skey, sval; int sind;
          ifstream inf(filename, ifstream::in);
          if (!inf.good() && myrank == 0) { if (ErrorMonitor->outfile) ErrorMonitor->outfile << "Can not open parameter file " << filename << endl; MPI_Abort(MPI_COMM_WORLD, 1); }
          for (int i = 1; inf.good(); i++) { inf.getline(pline, LEN); str = pline;
            int status = misc::parse_parts(str, sgrp, skey, sval, sind);
            if (status == -1) { if (ErrorMonitor->outfile) ErrorMonitor->outfile << "error reading parameter file " << filename << " in line " << i << endl; MPI_Abort(MPI_COMM_WORLD, 1); }
            else if (status == 0) continue;
            if (sgrp == "BSSN" && skey == "BH_num") { BH_NM = atoi(sval.c_str()); break; } } inf.close(); }
        Porg_here = new double[3*BH_NM]; Pmom = new double[3*BH_NM]; Spin = new double[3*BH_NM]; Mass = new double[BH_NM]; Qchar = new double[BH_NM];
        { const int LEN = 256; char pline[LEN]; string str, sgrp, skey, sval; int sind;
          ifstream inf(filename, ifstream::in);
          if (!inf.good() && myrank == 0) { if (ErrorMonitor->outfile) ErrorMonitor->outfile << "Can not open parameter file " << filename << endl; MPI_Abort(MPI_COMM_WORLD, 1); }
          for (int i = 1; inf.good(); i++) { inf.getline(pline, LEN); str = pline;
            int status = misc::parse_parts(str, sgrp, skey, sval, sind);
            if (status == -1) { if (ErrorMonitor->outfile) ErrorMonitor->outfile << "error reading parameter file " << filename << " in line " << i << endl; MPI_Abort(MPI_COMM_WORLD, 1); }
            else if (status == 0) continue;
            if (sgrp == "BSSN" && sind < BH_NM) {
              if (skey == "Mass") Mass[sind] = atof(sval.c_str());
              else if (skey == "Qchar") { Qchar[sind] = atof(sval.c_str()); if (myrank == 0) cout << "black hole #" << sind << " has electric charge " << Qchar[sind] << endl; }
              else if (skey == "Porgx") Porg_here[sind*3] = atof(sval.c_str());
              else if (skey == "Porgy") Porg_here[sind*3+1] = atof(sval.c_str());
              else if (skey == "Porgz") Porg_here[sind*3+2] = atof(sval.c_str());
              else if (skey == "Spinx") Spin[sind*3] = atof(sval.c_str());
              else if (skey == "Spiny") Spin[sind*3+1] = atof(sval.c_str());
              else if (skey == "Spinz") Spin[sind*3+2] = atof(sval.c_str());
              else if (skey == "Pmomx") Pmom[sind*3] = atof(sval.c_str());
              else if (skey == "Pmomy") Pmom[sind*3+1] = atof(sval.c_str());
              else if (skey == "Pmomz") Pmom[sind*3+2] = atof(sval.c_str());
            } } inf.close(); }
        int order = 6; Ansorg read_ansorg("Ansorg.psid", order);
        for (int lev = 0; lev < GH->levels; lev++) {
          MyList<Patch> *Pp = GH->PatL[lev];
          while (Pp) { MyList<Block> *BL = Pp->data->blb;
            while (BL) { Block *cg = BL->data; if (myrank == cg->rank) {
              for (int k = 0; k < cg->shape[2]; k++) for (int j = 0; j < cg->shape[1]; j++) for (int i = 0; i < cg->shape[0]; i++)
                cg->fgfs[phi0->sgfn][i+j*cg->shape[0]+k*cg->shape[0]*cg->shape[1]] = read_ansorg.ps_u_at_xyz(cg->X[0][i], cg->X[1][j], cg->X[2][k]);
              f_get_ansorg_nbhs_em(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                cg->fgfs[Lap0->sgfn], cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn],
                cg->fgfs[Ex0->sgfn], cg->fgfs[Ey0->sgfn], cg->fgfs[Ez0->sgfn], cg->fgfs[Bx0->sgfn], cg->fgfs[By0->sgfn], cg->fgfs[Bz0->sgfn],
                cg->fgfs[Kpsi0->sgfn], cg->fgfs[Kphi0->sgfn], Mass, Qchar, Porg_here, Pmom, Spin, BH_NM);
            } if (BL == Pp->data->ble) break; BL = BL->next; } Pp = Pp->next; } }
#ifdef WithShell
        /* shell patch - reserved */
#endif
        delete[] Porg_here;
    }
}

//================================================================================================



//================================================================================================
// 该成员函数用于用解析函数设定带电磁场的数值相对论初值
// Setup_Initial_Data
//================================================================================================

void bssnEM_class::Setup_Initial_Data()
{
    if (checkrun) {
        CheckPoint->readcheck_cgh(PhysTime, GH, myrank, nprocs, Symmetry);
#ifdef WithShell
        CheckPoint->readcheck_sh(SH, myrank);
#endif
    } else {
        if (myrank == 0) cout << "Setup initial data for head on identical charge-mass ratio black holes." << endl;
        char filename[50];
        { map<string, string>::iterator iter = parameters::str_par.find("inputpar");
          if (iter != parameters::str_par.end()) strcpy(filename, (iter->second).c_str());
          else { cout << "Error inputpar" << endl; exit(0); } }
        int BH_NM; double *Porg_here, *Qchar_here, *Pmom_here, *Spin_here, *Mass_here;
        { const int LEN = 256; char pline[LEN]; string str, sgrp, skey, sval; int sind;
          ifstream inf(filename, ifstream::in);
          if (!inf.good() && myrank == 0) { if (ErrorMonitor->outfile) ErrorMonitor->outfile << "Can not open parameter file " << filename << endl; MPI_Abort(MPI_COMM_WORLD, 1); }
          for (int i = 1; inf.good(); i++) { inf.getline(pline, LEN); str = pline;
            int status = misc::parse_parts(str, sgrp, skey, sval, sind);
            if (status == -1) { if (ErrorMonitor->outfile) ErrorMonitor->outfile << "error reading parameter file " << filename << " in line " << i << endl; MPI_Abort(MPI_COMM_WORLD, 1); }
            else if (status == 0) continue;
            if (sgrp == "BSSN" && skey == "BH_num") { BH_NM = atoi(sval.c_str()); break; } } inf.close(); }
        Porg_here = new double[3*BH_NM]; Pmom_here = new double[3*BH_NM]; Spin_here = new double[3*BH_NM]; Mass_here = new double[BH_NM]; Qchar_here = new double[BH_NM];
        { const int LEN = 256; char pline[LEN]; string str, sgrp, skey, sval; int sind;
          ifstream inf(filename, ifstream::in);
          if (!inf.good() && myrank == 0) { if (ErrorMonitor->outfile) ErrorMonitor->outfile << "Can not open parameter file " << filename << endl; MPI_Abort(MPI_COMM_WORLD, 1); }
          for (int i = 1; inf.good(); i++) { inf.getline(pline, LEN); str = pline;
            int status = misc::parse_parts(str, sgrp, skey, sval, sind);
            if (status == -1) { if (ErrorMonitor->outfile) ErrorMonitor->outfile << "error reading parameter file " << filename << " in line " << i << endl; MPI_Abort(MPI_COMM_WORLD, 1); }
            else if (status == 0) continue;
            if (sgrp == "BSSN" && sind < BH_NM) {
              if (skey == "Mass") Mass_here[sind] = atof(sval.c_str());
              else if (skey == "Qchar") { Qchar_here[sind] = atof(sval.c_str()); if (myrank == 0) cout << "black hole #" << sind << " has electric charge " << Qchar_here[sind] << endl; }
              else if (skey == "Porgx") Porg_here[sind*3] = atof(sval.c_str());
              else if (skey == "Porgy") Porg_here[sind*3+1] = atof(sval.c_str());
              else if (skey == "Porgz") Porg_here[sind*3+2] = atof(sval.c_str());
              else if (skey == "Spinx") Spin_here[sind*3] = atof(sval.c_str());
              else if (skey == "Spiny") Spin_here[sind*3+1] = atof(sval.c_str());
              else if (skey == "Spinz") Spin_here[sind*3+2] = atof(sval.c_str());
              else if (skey == "Pmomx") Pmom_here[sind*3] = atof(sval.c_str());
              else if (skey == "Pmomy") Pmom_here[sind*3+1] = atof(sval.c_str());
              else if (skey == "Pmomz") Pmom_here[sind*3+2] = atof(sval.c_str());
            } } inf.close(); }
        for (int lev = 0; lev < GH->levels; lev++) {
          MyList<Patch> *Pp = GH->PatL[lev];
          while (Pp) { MyList<Block> *BL = Pp->data->blb;
            while (BL) { Block *cg = BL->data; if (myrank == cg->rank) {
              f_get_initial_nbhsem(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                cg->fgfs[phi0->sgfn], cg->fgfs[trK0->sgfn],
                cg->fgfs[gxx0->sgfn], cg->fgfs[gxy0->sgfn], cg->fgfs[gxz0->sgfn], cg->fgfs[gyy0->sgfn], cg->fgfs[gyz0->sgfn], cg->fgfs[gzz0->sgfn],
                cg->fgfs[Axx0->sgfn], cg->fgfs[Axy0->sgfn], cg->fgfs[Axz0->sgfn], cg->fgfs[Ayy0->sgfn], cg->fgfs[Ayz0->sgfn], cg->fgfs[Azz0->sgfn],
                cg->fgfs[Gmx0->sgfn], cg->fgfs[Gmy0->sgfn], cg->fgfs[Gmz0->sgfn],
                cg->fgfs[Lap0->sgfn], cg->fgfs[Sfx0->sgfn], cg->fgfs[Sfy0->sgfn], cg->fgfs[Sfz0->sgfn],
                cg->fgfs[dtSfx0->sgfn], cg->fgfs[dtSfy0->sgfn], cg->fgfs[dtSfz0->sgfn],
                cg->fgfs[Ex0->sgfn], cg->fgfs[Ey0->sgfn], cg->fgfs[Ez0->sgfn], cg->fgfs[Bx0->sgfn], cg->fgfs[By0->sgfn], cg->fgfs[Bz0->sgfn],
                cg->fgfs[Kpsi0->sgfn], cg->fgfs[Kphi0->sgfn],
                Mass_here, Qchar_here, Porg_here, Pmom_here, Spin_here, BH_NM);
            } if (BL == Pp->data->ble) break; BL = BL->next; } Pp = Pp->next; } }
#ifdef WithShell
        /* shell patch - reserved */
#endif
        delete[] Porg_here; delete[] Mass_here; delete[] Qchar_here; delete[] Pmom_here; delete[] Spin_here;
    }
}

//================================================================================================



//================================================================================================
// 该成员函数设定了时间演化过程中的单步时间演化
// Step_GPU: Full standalone RK4 evolution with BSSN + EM RHS at each substep
//================================================================================================

// Structure mirrors bssn_step_gpu.C exactly (MAPBH=0/1, AGM=0/1, SommerType=0/1, RPS=0).
// At each substep (predictor + 3 correctors), we add EM symmetry BD and EM RHS
// AFTER enforce_ga and BEFORE BSSN symmetry BD, so the stress-energy tensor
// (rho, S_i, S_ij) is freshly computed before BSSN RHS reads it.
// The fused RK4 + Sommerfeld loop already processes ALL variables in the lists
// (including the 8 EM fields), so EM variables are integrated alongside BSSN ones.
//
// Variable mapping per substep:
//   Predictor (iter_count=0):  EM state=StateList (Ex0..Kphi0), EM RHS output=RHSList (Ex_rhs..)
//                              BSSN state=StateList (phi0..), BSSN RHS output=RHSList (phi_rhs..)
//   Corrector (iter_count=1..3): EM state=SynchList_pre (Ex..Kphi), EM RHS output=SynchList_cor (Ex1..)
//                                BSSN state=SynchList_pre (phi..), BSSN RHS output=SynchList_cor (phi1..)
//   Stress-energy (rho, S_i, S_ij) is in MiscList, written by EM RHS, read by BSSN RHS.

void bssnEM_class::Step_GPU(int lev, int YN)
{
    setpbh(BH_num, Porg0, Mass, BH_num_input);

    double dT_lev = dT * pow(0.5, Mymax(lev, trfls));

    // BH position RK4 predictor/corrector is unconditional in EM class
    // (matches CPU bssnEM_class.C which uses RK4 regardless of MAPBH)

#ifdef With_AHF
    AH_Step_Find(lev, dT_lev);
#endif

    bool BB = fgt(PhysTime, StartTime, dT_lev / 2);
    double ndeps = numepss;
    if (lev < GH->movls) ndeps = numepsb;
    double TRK4 = PhysTime;
    int iter_count = 0;
    int pre = 0, cor = 1;

    MyList<ss_patch> *sPp;
    // Predictor
    MyList<Patch> *Pp = GH->PatL[lev];
    while (Pp) {
        MyList<Block> *BP = Pp->data->blb;
        while (BP) {
            Block *cg = BP->data;
            if (myrank == cg->rank) {

#if (AGM == 0)
                gpu_enforce_ga_launch(
                    cg->stream, cg->shape, cg->shape_pad,
                    cg->d_fgfs[gxx0->sgfn], cg->d_fgfs[gxy0->sgfn], cg->d_fgfs[gxz0->sgfn], cg->d_fgfs[gyy0->sgfn], cg->d_fgfs[gyz0->sgfn], cg->d_fgfs[gzz0->sgfn],
                    cg->d_fgfs[Axx0->sgfn], cg->d_fgfs[Axy0->sgfn], cg->d_fgfs[Axz0->sgfn], cg->d_fgfs[Ayy0->sgfn], cg->d_fgfs[Ayz0->sgfn], cg->d_fgfs[Azz0->sgfn]
                );
#endif

                // Pre-fill EM and BSSN field ghost zones with symmetry-mirrored
                // values before EM RHS so the SMEM-based empart_rhs_kernel can
                // read all 19 input field ghost zones (Lap, chi, beta_i, dxx..dzz,
                // Kpsi/Kphi, Ex..Bz) by direct global load. Ordered as:
                //   sym_bd_empart  (8 EM fields)
                //   sym_bd_fused   (12 BSSN fields the EM RHS reads)
                //   empart_rhs     (consumes both)
                //   bssn_rhs       (also consumes BSSN ghost — still valid)
                // The original sequence had sym_bd_fused AFTER empart_rhs, which
                // was safe only because the old per-point kernel mirrored on
                // the fly inside d_fderivs_point_gw.
                gpu_apply_symmetry_bd_empart_launch(
                    cg->stream, cg->shape[0], cg->shape[1], cg->shape[2],
                    cg->shape_pad[0], cg->shape_pad[1], cg->shape_pad[2],
                    cg->d_fgfs[Ex0->sgfn], cg->d_fgfs[Ey0->sgfn], cg->d_fgfs[Ez0->sgfn],
                    cg->d_fgfs[Bx0->sgfn], cg->d_fgfs[By0->sgfn], cg->d_fgfs[Bz0->sgfn],
                    cg->d_fgfs[Kpsi0->sgfn], cg->d_fgfs[Kphi0->sgfn]
                );

                gpu_apply_symmetry_bd_fused_launch(
                    cg->stream, cg->shape[0], cg->shape[1], cg->shape[2], cg->shape_pad[0], cg->shape_pad[1], cg->shape_pad[2], Symmetry,
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

                // EM RHS (predictor): EM state from StateList, RHS output to RHSList
                // Also computes stress-energy tensor (rho, S_i, S_ij) into MiscList
                gpu_compute_rhs_empart_launch(
                    cg->stream,
                    cg->shape, cg->shape_pad, cg->d_X[0], cg->d_X[1], cg->d_X[2],
                    cg->d_fgfs[phi0->sgfn],
                    cg->d_fgfs[Jx->sgfn], cg->d_fgfs[Jy->sgfn], cg->d_fgfs[Jz->sgfn], cg->d_fgfs[qchar->sgfn],
                    cg->d_fgfs[gxx0->sgfn], cg->d_fgfs[gxy0->sgfn], cg->d_fgfs[gxz0->sgfn],
                    cg->d_fgfs[gyy0->sgfn], cg->d_fgfs[gyz0->sgfn], cg->d_fgfs[gzz0->sgfn],
                    cg->d_fgfs[Lap0->sgfn],
                    cg->d_fgfs[Sfx0->sgfn], cg->d_fgfs[Sfy0->sgfn], cg->d_fgfs[Sfz0->sgfn],
                    cg->d_fgfs[trK0->sgfn],
                    cg->d_fgfs[Ex0->sgfn], cg->d_fgfs[Ey0->sgfn], cg->d_fgfs[Ez0->sgfn],
                    cg->d_fgfs[Bx0->sgfn], cg->d_fgfs[By0->sgfn], cg->d_fgfs[Bz0->sgfn],
                    cg->d_fgfs[Kpsi0->sgfn], cg->d_fgfs[Kphi0->sgfn],
                    cg->d_fgfs[Ex_rhs->sgfn], cg->d_fgfs[Ey_rhs->sgfn], cg->d_fgfs[Ez_rhs->sgfn],
                    cg->d_fgfs[Bx_rhs->sgfn], cg->d_fgfs[By_rhs->sgfn], cg->d_fgfs[Bz_rhs->sgfn],
                    cg->d_fgfs[Kpsi_rhs->sgfn], cg->d_fgfs[Kphi_rhs->sgfn],
                    cg->d_fgfs[rho->sgfn],
                    cg->d_fgfs[Sx->sgfn], cg->d_fgfs[Sy->sgfn], cg->d_fgfs[Sz->sgfn],
                    cg->d_fgfs[Sxx->sgfn], cg->d_fgfs[Sxy->sgfn], cg->d_fgfs[Sxz->sgfn],
                    cg->d_fgfs[Syy->sgfn], cg->d_fgfs[Syz->sgfn], cg->d_fgfs[Szz->sgfn],
                    Symmetry, lev, ndeps
                );

                // BSSN RHS (predictor): state from StateList, RHS output to RHSList
                // Reads stress-energy from MiscList (freshly computed by EM RHS above)
                gpu_compute_rhs_bssn_launch(
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
                    cg->d_fgfs[phi_rhs->sgfn], cg->d_fgfs[trK_rhs->sgfn],
                    cg->d_fgfs[gxx_rhs->sgfn], cg->d_fgfs[gxy_rhs->sgfn], cg->d_fgfs[gxz_rhs->sgfn],
                    cg->d_fgfs[gyy_rhs->sgfn], cg->d_fgfs[gyz_rhs->sgfn], cg->d_fgfs[gzz_rhs->sgfn],
                    cg->d_fgfs[Axx_rhs->sgfn], cg->d_fgfs[Axy_rhs->sgfn], cg->d_fgfs[Axz_rhs->sgfn],
                    cg->d_fgfs[Ayy_rhs->sgfn], cg->d_fgfs[Ayz_rhs->sgfn], cg->d_fgfs[Azz_rhs->sgfn],
                    cg->d_fgfs[Gmx_rhs->sgfn], cg->d_fgfs[Gmy_rhs->sgfn], cg->d_fgfs[Gmz_rhs->sgfn],
                    cg->d_fgfs[Lap_rhs->sgfn],
                    cg->d_fgfs[Sfx_rhs->sgfn], cg->d_fgfs[Sfy_rhs->sgfn], cg->d_fgfs[Sfz_rhs->sgfn],
                    cg->d_fgfs[dtSfx_rhs->sgfn], cg->d_fgfs[dtSfy_rhs->sgfn], cg->d_fgfs[dtSfz_rhs->sgfn],
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
#ifdef With_AHF
                AH_Mark_RHS_Scratch(*cg);
#endif


                // Fused RK4 + Sommerfeld for ALL variables (BSSN + EM) in the lists
                MyList<var> *varl0 = StateList, *varl = SynchList_pre, *varlrhs = RHSList;
                while (varl0) {
                    FusedRK4Vars rk4_vars; rk4_vars.num_vars = 0;
                    FusedSommerfeldVars somm_vars; somm_vars.num_vars = 0;
#if (SommerType == 1)
                    double *shib_f_out[MAX_FUSED_VARS];
                    int shib_count = 0;
#endif

                    while (varl0 && rk4_vars.num_vars < MAX_FUSED_VARS) {
                        rk4_vars.f0[rk4_vars.num_vars] = cg->d_fgfs[varl0->data->sgfn];
                        rk4_vars.f1[rk4_vars.num_vars] = cg->d_fgfs[varl->data->sgfn];
                        rk4_vars.f_rhs[rk4_vars.num_vars] = cg->d_fgfs[varlrhs->data->sgfn];

                        somm_vars.f0[somm_vars.num_vars] = cg->d_fgfs[varl0->data->sgfn];
                        somm_vars.f_out[somm_vars.num_vars] = lev == 0 ? cg->d_fgfs[varlrhs->data->sgfn] : cg->d_fgfs[varl->data->sgfn];
                        somm_vars.propspeed[somm_vars.num_vars] = varl0->data->propspeed;
                        somm_vars.SoA0[somm_vars.num_vars] = varl0->data->SoA[0];
                        somm_vars.SoA1[somm_vars.num_vars] = varl0->data->SoA[1];
                        somm_vars.SoA2[somm_vars.num_vars] = varl0->data->SoA[2];
#if (SommerType == 1)
                        shib_f_out[shib_count++] = cg->d_fgfs[varl->data->sgfn];
#endif

                        rk4_vars.num_vars++;
                        somm_vars.num_vars++;

                        varl0 = varl0->next; varl = varl->next; varlrhs = varlrhs->next;
                    }

#if (SommerType == 0)
#ifndef WithShell
                    if (lev == 0) {
                        gpu_fused_sommerfeld_routbam_launch(
                            cg->stream, cg->shape, cg->shape_pad, cg->d_X[0], cg->d_X[1], cg->d_X[2],
                            Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                            somm_vars, Symmetry
                        );
                    }
#endif
#endif

                    gpu_fused_rungekutta4_rout_launch(cg->stream, cg->shape, cg->shape_pad, dT_lev, rk4_vars, iter_count);

#ifndef WithShell
                    if (lev > 0) {
#endif
                        gpu_fused_sommerfeld_rout_launch(
                            cg->stream, cg->shape, cg->shape_pad, cg->d_X[0], cg->d_X[1], cg->d_X[2],
                            Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                            dT_lev, cg->d_fgfs[phi0->sgfn], cg->d_fgfs[Lap0->sgfn],
                            somm_vars, Symmetry, cor
                        );
#ifndef WithShell
                    }
#endif

#if (SommerType == 1)
                    if (lev == 0) {
                        FusedSommerfeldVars somm_vars_shib;
                        somm_vars_shib.num_vars = somm_vars.num_vars;
                        for (int iv = 0; iv < somm_vars.num_vars; iv++) {
                            somm_vars_shib.f0[iv] = somm_vars.f0[iv];
                            somm_vars_shib.f_out[iv] = shib_f_out[iv];
                            somm_vars_shib.propspeed[iv] = somm_vars.propspeed[iv];
                            somm_vars_shib.SoA0[iv] = somm_vars.SoA0[iv];
                            somm_vars_shib.SoA1[iv] = somm_vars.SoA1[iv];
                            somm_vars_shib.SoA2[iv] = somm_vars.SoA2[iv];
                        }
                        gpu_fused_sommerfeld_rout_launch(
                            cg->stream, cg->shape, cg->shape_pad, cg->d_X[0], cg->d_X[1], cg->d_X[2],
                            Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                            dT_lev, cg->d_fgfs[phi0->sgfn], cg->d_fgfs[Lap0->sgfn],
                            somm_vars_shib, Symmetry, pre
                        );
                    }
#endif
                }

                gpu_lowerboundset_launch(cg->stream, cg->shape, cg->shape_pad, cg->d_fgfs[phi->sgfn], chitiny);

            }
            if (BP == Pp->data->ble) break;
            BP = BP->next;
        }
        Pp = Pp->next;
    }

#ifdef WithShell
    if (lev == 0)
    {
        int ERROR = 0;
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
                    gpu_enforce_ga_launch(
                        cg->stream, cg->shape, cg->shape_pad,
                        cg->d_fgfs[gxx0->sgfn], cg->d_fgfs[gxy0->sgfn], cg->d_fgfs[gxz0->sgfn], cg->d_fgfs[gyy0->sgfn], cg->d_fgfs[gyz0->sgfn], cg->d_fgfs[gzz0->sgfn],
                        cg->d_fgfs[Axx0->sgfn], cg->d_fgfs[Axy0->sgfn], cg->d_fgfs[Axz0->sgfn], cg->d_fgfs[Ayy0->sgfn], cg->d_fgfs[Ayz0->sgfn], cg->d_fgfs[Azz0->sgfn]
                    );
#endif

                    // [REMOVED] EM symmetry boundary - not in CPU version

                    // EM RHS (Shell predictor): shell coordinate transform + stress-energy
                    gpu_compute_rhs_empart_ss_launch(
                        cg->stream,
                        cg->shape, cg->shape_pad,
                        cg->d_fgfs[fngfs+ShellPatch::gx], cg->d_fgfs[fngfs+ShellPatch::gy], cg->d_fgfs[fngfs+ShellPatch::gz],
                        cg->d_fgfs[fngfs+ShellPatch::drhodx], cg->d_fgfs[fngfs+ShellPatch::drhody], cg->d_fgfs[fngfs+ShellPatch::drhodz],
                        cg->d_fgfs[fngfs+ShellPatch::dsigmadx], cg->d_fgfs[fngfs+ShellPatch::dsigmady], cg->d_fgfs[fngfs+ShellPatch::dsigmadz],
                        cg->d_fgfs[fngfs+ShellPatch::dRdx], cg->d_fgfs[fngfs+ShellPatch::dRdy], cg->d_fgfs[fngfs+ShellPatch::dRdz],
                        cg->d_fgfs[phi0->sgfn],
                        cg->d_fgfs[Jx->sgfn], cg->d_fgfs[Jy->sgfn], cg->d_fgfs[Jz->sgfn], cg->d_fgfs[qchar->sgfn],
                        cg->d_fgfs[gxx0->sgfn], cg->d_fgfs[gxy0->sgfn], cg->d_fgfs[gxz0->sgfn],
                        cg->d_fgfs[gyy0->sgfn], cg->d_fgfs[gyz0->sgfn], cg->d_fgfs[gzz0->sgfn],
                        cg->d_fgfs[Lap0->sgfn],
                        cg->d_fgfs[Sfx0->sgfn], cg->d_fgfs[Sfy0->sgfn], cg->d_fgfs[Sfz0->sgfn],
                        cg->d_fgfs[trK0->sgfn],
                        cg->d_fgfs[Ex0->sgfn], cg->d_fgfs[Ey0->sgfn], cg->d_fgfs[Ez0->sgfn],
                        cg->d_fgfs[Bx0->sgfn], cg->d_fgfs[By0->sgfn], cg->d_fgfs[Bz0->sgfn],
                        cg->d_fgfs[Kpsi0->sgfn], cg->d_fgfs[Kphi0->sgfn],
                        cg->d_fgfs[Ex_rhs->sgfn], cg->d_fgfs[Ey_rhs->sgfn], cg->d_fgfs[Ez_rhs->sgfn],
                        cg->d_fgfs[Bx_rhs->sgfn], cg->d_fgfs[By_rhs->sgfn], cg->d_fgfs[Bz_rhs->sgfn],
                        cg->d_fgfs[Kpsi_rhs->sgfn], cg->d_fgfs[Kphi_rhs->sgfn],
                        cg->d_fgfs[rho->sgfn],
                        cg->d_fgfs[Sx->sgfn], cg->d_fgfs[Sy->sgfn], cg->d_fgfs[Sz->sgfn],
                        cg->d_fgfs[Sxx->sgfn], cg->d_fgfs[Sxy->sgfn], cg->d_fgfs[Sxz->sgfn],
                        cg->d_fgfs[Syy->sgfn], cg->d_fgfs[Syz->sgfn], cg->d_fgfs[Szz->sgfn],
                        Symmetry, lev, numepsh, sPp->data->sst
                    );

                    gpu_compute_rhs_bssn_ss_launch(
                        cg->stream, cg->shape, cg->shape_pad, TRK4,
                        cg->d_X[0], cg->d_X[1], cg->d_X[2],
                        cg->d_fgfs[fngfs+ShellPatch::gx], cg->d_fgfs[fngfs+ShellPatch::gy], cg->d_fgfs[fngfs+ShellPatch::gz],
                        cg->d_fgfs[fngfs+ShellPatch::drhodx], cg->d_fgfs[fngfs+ShellPatch::drhody], cg->d_fgfs[fngfs+ShellPatch::drhodz],
                        cg->d_fgfs[fngfs+ShellPatch::dsigmadx], cg->d_fgfs[fngfs+ShellPatch::dsigmady], cg->d_fgfs[fngfs+ShellPatch::dsigmadz],
                        cg->d_fgfs[fngfs+ShellPatch::dRdx], cg->d_fgfs[fngfs+ShellPatch::dRdy], cg->d_fgfs[fngfs+ShellPatch::dRdz],
                        cg->d_fgfs[fngfs+ShellPatch::drhodxx], cg->d_fgfs[fngfs+ShellPatch::drhodxy], cg->d_fgfs[fngfs+ShellPatch::drhodxz],
                        cg->d_fgfs[fngfs+ShellPatch::drhodyy], cg->d_fgfs[fngfs+ShellPatch::drhodyz], cg->d_fgfs[fngfs+ShellPatch::drhodzz],
                        cg->d_fgfs[fngfs+ShellPatch::dsigmadxx], cg->d_fgfs[fngfs+ShellPatch::dsigmadxy], cg->d_fgfs[fngfs+ShellPatch::dsigmadxz],
                        cg->d_fgfs[fngfs+ShellPatch::dsigmadyy], cg->d_fgfs[fngfs+ShellPatch::dsigmadyz], cg->d_fgfs[fngfs+ShellPatch::dsigmadzz],
                        cg->d_fgfs[fngfs+ShellPatch::dRdxx], cg->d_fgfs[fngfs+ShellPatch::dRdxy], cg->d_fgfs[fngfs+ShellPatch::dRdxz],
                        cg->d_fgfs[fngfs+ShellPatch::dRdyy], cg->d_fgfs[fngfs+ShellPatch::dRdyz], cg->d_fgfs[fngfs+ShellPatch::dRdzz],
                        cg->d_fgfs[phi0->sgfn], cg->d_fgfs[trK0->sgfn],
                        cg->d_fgfs[gxx0->sgfn], cg->d_fgfs[gxy0->sgfn], cg->d_fgfs[gxz0->sgfn],
                        cg->d_fgfs[gyy0->sgfn], cg->d_fgfs[gyz0->sgfn], cg->d_fgfs[gzz0->sgfn],
                        cg->d_fgfs[Axx0->sgfn], cg->d_fgfs[Axy0->sgfn], cg->d_fgfs[Axz0->sgfn],
                        cg->d_fgfs[Ayy0->sgfn], cg->d_fgfs[Ayz0->sgfn], cg->d_fgfs[Azz0->sgfn],
                        cg->d_fgfs[Gmx0->sgfn], cg->d_fgfs[Gmy0->sgfn], cg->d_fgfs[Gmz0->sgfn],
                        cg->d_fgfs[Lap0->sgfn],
                        cg->d_fgfs[Sfx0->sgfn], cg->d_fgfs[Sfy0->sgfn], cg->d_fgfs[Sfz0->sgfn],
                        cg->d_fgfs[dtSfx0->sgfn], cg->d_fgfs[dtSfy0->sgfn], cg->d_fgfs[dtSfz0->sgfn],
                        cg->d_fgfs[phi_rhs->sgfn], cg->d_fgfs[trK_rhs->sgfn],
                        cg->d_fgfs[gxx_rhs->sgfn], cg->d_fgfs[gxy_rhs->sgfn], cg->d_fgfs[gxz_rhs->sgfn],
                        cg->d_fgfs[gyy_rhs->sgfn], cg->d_fgfs[gyz_rhs->sgfn], cg->d_fgfs[gzz_rhs->sgfn],
                        cg->d_fgfs[Axx_rhs->sgfn], cg->d_fgfs[Axy_rhs->sgfn], cg->d_fgfs[Axz_rhs->sgfn],
                        cg->d_fgfs[Ayy_rhs->sgfn], cg->d_fgfs[Ayz_rhs->sgfn], cg->d_fgfs[Azz_rhs->sgfn],
                        cg->d_fgfs[Gmx_rhs->sgfn], cg->d_fgfs[Gmy_rhs->sgfn], cg->d_fgfs[Gmz_rhs->sgfn],
                        cg->d_fgfs[Lap_rhs->sgfn],
                        cg->d_fgfs[Sfx_rhs->sgfn], cg->d_fgfs[Sfy_rhs->sgfn], cg->d_fgfs[Sfz_rhs->sgfn],
                        cg->d_fgfs[dtSfx_rhs->sgfn], cg->d_fgfs[dtSfy_rhs->sgfn], cg->d_fgfs[dtSfz_rhs->sgfn],
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
                        cg->d_fgfs[Cons_Ham->sgfn], cg->d_fgfs[Cons_Px->sgfn], cg->d_fgfs[Cons_Py->sgfn], cg->d_fgfs[Cons_Pz->sgfn],
                        cg->d_fgfs[Cons_Gx->sgfn], cg->d_fgfs[Cons_Gy->sgfn], cg->d_fgfs[Cons_Gz->sgfn],
                        Symmetry, lev, numepsh, sPp->data->sst, pre
                    );
#ifdef With_AHF
                    AH_Mark_RHS_Scratch(*cg);
#endif
                    cudaError_t err = cudaGetLastError();
                    if (err != cudaSuccess) {
                        cout << "GPU SS RHS error: " << cudaGetErrorString(err) << " sst=" << sPp->data->sst << endl;
                        ERROR = 1;
                    }

                    // Fused GPU sommerfeld + rungekutta4 for Shell patches (predictor)
                    {
                        MyList<var> *varl0 = StateList, *varl = SynchList_pre, *varlrhs = RHSList;
                        while (varl0) {
                            FusedRK4Vars rk4_vars; rk4_vars.num_vars = 0;
                            FusedSommerfeldVars somm_vars; somm_vars.num_vars = 0;

                            while (varl0 && rk4_vars.num_vars < MAX_FUSED_VARS) {
                                rk4_vars.f0[rk4_vars.num_vars] = cg->d_fgfs[varl0->data->sgfn];
                                rk4_vars.f1[rk4_vars.num_vars] = cg->d_fgfs[varl->data->sgfn];
                                rk4_vars.f_rhs[rk4_vars.num_vars] = cg->d_fgfs[varlrhs->data->sgfn];

                                somm_vars.f0[somm_vars.num_vars] = cg->d_fgfs[varl0->data->sgfn];
                                somm_vars.f_out[somm_vars.num_vars] = cg->d_fgfs[varlrhs->data->sgfn];
                                somm_vars.propspeed[somm_vars.num_vars] = varl0->data->propspeed;
                                somm_vars.SoA0[somm_vars.num_vars] = varl0->data->SoA[0];
                                somm_vars.SoA1[somm_vars.num_vars] = varl0->data->SoA[1];
                                somm_vars.SoA2[somm_vars.num_vars] = varl0->data->SoA[2];

                                rk4_vars.num_vars++;
                                somm_vars.num_vars++;

                                varl0 = varl0->next; varl = varl->next; varlrhs = varlrhs->next;
                            }

#if (SommerType == 0)
                            gpu_fused_sommerfeld_routbam_ss_launch(
                                cg->stream, cg->shape, cg->shape_pad, cg->d_X[0], cg->d_X[1], cg->d_X[2],
                                sPp->data->bbox[0], sPp->data->bbox[1], sPp->data->bbox[2], sPp->data->bbox[3], sPp->data->bbox[4], sPp->data->bbox[5],
                                somm_vars, Symmetry
                            );
#endif

                            gpu_fused_rungekutta4_rout_launch(cg->stream, cg->shape, cg->shape_pad, dT_lev, rk4_vars, iter_count);
                        }
                    }
                    // EM field falloff boundary on Shell patch (predictor)
                    // f ~ r^(-n): n=2 for E,B; n=3 for Kpsi,Kphi
                    gpu_falloff_ss_launch(cg->stream,
                        cg->shape[0], cg->shape[1], cg->shape[2],
                        cg->shape_pad[0], cg->shape_pad[1],
                        cg->d_fgfs[Ex->sgfn], cg->d_fgfs[fngfs+ShellPatch::gz], 2, sPp->data->bbox[5]);
                    gpu_falloff_ss_launch(cg->stream,
                        cg->shape[0], cg->shape[1], cg->shape[2],
                        cg->shape_pad[0], cg->shape_pad[1],
                        cg->d_fgfs[Ey->sgfn], cg->d_fgfs[fngfs+ShellPatch::gz], 2, sPp->data->bbox[5]);
                    gpu_falloff_ss_launch(cg->stream,
                        cg->shape[0], cg->shape[1], cg->shape[2],
                        cg->shape_pad[0], cg->shape_pad[1],
                        cg->d_fgfs[Ez->sgfn], cg->d_fgfs[fngfs+ShellPatch::gz], 2, sPp->data->bbox[5]);
                    gpu_falloff_ss_launch(cg->stream,
                        cg->shape[0], cg->shape[1], cg->shape[2],
                        cg->shape_pad[0], cg->shape_pad[1],
                        cg->d_fgfs[Bx->sgfn], cg->d_fgfs[fngfs+ShellPatch::gz], 2, sPp->data->bbox[5]);
                    gpu_falloff_ss_launch(cg->stream,
                        cg->shape[0], cg->shape[1], cg->shape[2],
                        cg->shape_pad[0], cg->shape_pad[1],
                        cg->d_fgfs[By->sgfn], cg->d_fgfs[fngfs+ShellPatch::gz], 2, sPp->data->bbox[5]);
                    gpu_falloff_ss_launch(cg->stream,
                        cg->shape[0], cg->shape[1], cg->shape[2],
                        cg->shape_pad[0], cg->shape_pad[1],
                        cg->d_fgfs[Bz->sgfn], cg->d_fgfs[fngfs+ShellPatch::gz], 2, sPp->data->bbox[5]);
                    gpu_falloff_ss_launch(cg->stream,
                        cg->shape[0], cg->shape[1], cg->shape[2],
                        cg->shape_pad[0], cg->shape_pad[1],
                        cg->d_fgfs[Kpsi->sgfn], cg->d_fgfs[fngfs+ShellPatch::gz], 3, sPp->data->bbox[5]);
                    gpu_falloff_ss_launch(cg->stream,
                        cg->shape[0], cg->shape[1], cg->shape[2],
                        cg->shape_pad[0], cg->shape_pad[1],
                        cg->d_fgfs[Kphi->sgfn], cg->d_fgfs[fngfs+ShellPatch::gz], 3, sPp->data->bbox[5]);
                    gpu_lowerboundset_launch(cg->stream, cg->shape, cg->shape_pad, cg->d_fgfs[phi->sgfn], chitiny);
                }
                if (BP == sPp->data->ble) break;
                BP = BP->next;
            }
            sPp = sPp->next;
        }

        int erh = ERROR;
        MPI_Allreduce(&erh, &ERROR, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
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
    }
#endif

    Parallel::Sync_GPU(GH->PatL[lev], SynchList_pre, Symmetry);

#ifdef WithShell
    if (lev == 0)
    {
        double prev_clock, curr_clock;
        GPUManager::getInstance().synchronize_all();
        if (myrank == 0) prev_clock = MPI_Wtime();
        SH->Synch_GPU(SynchList_pre, Symmetry);
        if (myrank == 0)
        {
            curr_clock = MPI_Wtime();
            cout << "Shell stuff synchronization used " << (curr_clock - prev_clock) << " seconds! (RK substep " << iter_count << ")" << endl;
        }
    }
#endif

    // for black hole position — RK4 predictor (unconditional for EM)
    if (BH_num > 0 && lev == GH->levels - 1)
    {
        compute_Porg_rhs(Porg0, Porg_rhs, Sfx0, Sfy0, Sfz0, lev);
        for (int ithBH = 0; ithBH < BH_num; ithBH++)
        {
            f_rungekutta4_scalar(dT_lev, Porg0[ithBH][0], Porg[ithBH][0], Porg_rhs[ithBH][0], iter_count);
            f_rungekutta4_scalar(dT_lev, Porg0[ithBH][1], Porg[ithBH][1], Porg_rhs[ithBH][1], iter_count);
            f_rungekutta4_scalar(dT_lev, Porg0[ithBH][2], Porg[ithBH][2], Porg_rhs[ithBH][2], iter_count);
            if (Symmetry > 0) Porg[ithBH][2] = fabs(Porg[ithBH][2]);
            if (Symmetry == 2) { Porg[ithBH][0] = fabs(Porg[ithBH][0]); Porg[ithBH][1] = fabs(Porg[ithBH][1]); }
            if (!finite(Porg[ithBH][0]) || !finite(Porg[ithBH][1]) || !finite(Porg[ithBH][2]))
            {
                if (ErrorMonitor->outfile)
                    ErrorMonitor->outfile << "predictor step finds NaN for BH's position from ("
                                          << Porg0[ithBH][0] << "," << Porg0[ithBH][1] << "," << Porg0[ithBH][2] << ")" << endl;

                MyList<var> *DG_List = new MyList<var>(Sfx0);
                DG_List->insert(Sfx0); DG_List->insert(Sfy0); DG_List->insert(Sfz0);
                Parallel::Dump_Data(GH->PatL[lev], DG_List, 0, PhysTime, dT_lev);
                DG_List->clearList();
            }
        }
    }
    if (lev == a_lev) { AnalysisStuff_EM(lev, dT_lev); }

    // corrector
    for (iter_count = 1; iter_count < 4; iter_count++) {
        if (iter_count == 1 || iter_count == 3) TRK4 += dT_lev / 2;
        Pp = GH->PatL[lev];
        while (Pp) {
            MyList<Block> *BP = Pp->data->blb;
            while (BP) {
                Block *cg = BP->data;
                if (myrank == cg->rank) {
#if (AGM == 0)
                    gpu_enforce_ga_launch(
                        cg->stream, cg->shape, cg->shape_pad,
                        cg->d_fgfs[gxx->sgfn], cg->d_fgfs[gxy->sgfn], cg->d_fgfs[gxz->sgfn], cg->d_fgfs[gyy->sgfn], cg->d_fgfs[gyz->sgfn], cg->d_fgfs[gzz->sgfn],
                        cg->d_fgfs[Axx->sgfn], cg->d_fgfs[Axy->sgfn], cg->d_fgfs[Axz->sgfn], cg->d_fgfs[Ayy->sgfn], cg->d_fgfs[Ayz->sgfn], cg->d_fgfs[Azz->sgfn]
                    );
#elif (AGM == 1)
                    if (iter_count == 3)
                        gpu_enforce_ga_launch(
                            cg->stream, cg->shape, cg->shape_pad,
                            cg->d_fgfs[gxx->sgfn], cg->d_fgfs[gxy->sgfn], cg->d_fgfs[gxz->sgfn], cg->d_fgfs[gyy->sgfn], cg->d_fgfs[gyz->sgfn], cg->d_fgfs[gzz->sgfn],
                            cg->d_fgfs[Axx->sgfn], cg->d_fgfs[Axy->sgfn], cg->d_fgfs[Axz->sgfn], cg->d_fgfs[Ayy->sgfn], cg->d_fgfs[Ayz->sgfn], cg->d_fgfs[Azz->sgfn]
                        );
#endif

                    // Pre-mirror EM + BSSN ghost zones before EM RHS so the
                    // SMEM-based empart_rhs_kernel can global-load both groups
                    // of ghost cells. See predictor for full rationale.
                    gpu_apply_symmetry_bd_empart_launch(
                        cg->stream, cg->shape[0], cg->shape[1], cg->shape[2],
                        cg->shape_pad[0], cg->shape_pad[1], cg->shape_pad[2],
                        cg->d_fgfs[Ex->sgfn], cg->d_fgfs[Ey->sgfn], cg->d_fgfs[Ez->sgfn],
                        cg->d_fgfs[Bx->sgfn], cg->d_fgfs[By->sgfn], cg->d_fgfs[Bz->sgfn],
                        cg->d_fgfs[Kpsi->sgfn], cg->d_fgfs[Kphi->sgfn]
                    );

                    gpu_apply_symmetry_bd_fused_launch(
                        cg->stream, cg->shape[0], cg->shape[1], cg->shape[2], cg->shape_pad[0], cg->shape_pad[1], cg->shape_pad[2], Symmetry,
                        cg->d_fgfs[phi->sgfn], cg->d_fgfs[trK->sgfn],
                        cg->d_fgfs[gxx->sgfn], cg->d_fgfs[gxy->sgfn], cg->d_fgfs[gxz->sgfn],
                        cg->d_fgfs[gyy->sgfn], cg->d_fgfs[gyz->sgfn], cg->d_fgfs[gzz->sgfn],
                        cg->d_fgfs[Axx->sgfn], cg->d_fgfs[Axy->sgfn], cg->d_fgfs[Axz->sgfn],
                        cg->d_fgfs[Ayy->sgfn], cg->d_fgfs[Ayz->sgfn], cg->d_fgfs[Azz->sgfn],
                        cg->d_fgfs[Gmx->sgfn], cg->d_fgfs[Gmy->sgfn], cg->d_fgfs[Gmz->sgfn],
                        cg->d_fgfs[Lap->sgfn],
                        cg->d_fgfs[Sfx->sgfn], cg->d_fgfs[Sfy->sgfn], cg->d_fgfs[Sfz->sgfn],
                        cg->d_fgfs[dtSfx->sgfn], cg->d_fgfs[dtSfy->sgfn], cg->d_fgfs[dtSfz->sgfn]
                    );

                    // EM RHS (corrector): EM state from SynchList_pre, RHS output to SynchList_cor
                    gpu_compute_rhs_empart_launch(
                        cg->stream,
                        cg->shape, cg->shape_pad, cg->d_X[0], cg->d_X[1], cg->d_X[2],
                        cg->d_fgfs[phi->sgfn],
                        cg->d_fgfs[Jx->sgfn], cg->d_fgfs[Jy->sgfn], cg->d_fgfs[Jz->sgfn], cg->d_fgfs[qchar->sgfn],
                        cg->d_fgfs[gxx->sgfn], cg->d_fgfs[gxy->sgfn], cg->d_fgfs[gxz->sgfn],
                        cg->d_fgfs[gyy->sgfn], cg->d_fgfs[gyz->sgfn], cg->d_fgfs[gzz->sgfn],
                        cg->d_fgfs[Lap->sgfn],
                        cg->d_fgfs[Sfx->sgfn], cg->d_fgfs[Sfy->sgfn], cg->d_fgfs[Sfz->sgfn],
                        cg->d_fgfs[trK->sgfn],
                        cg->d_fgfs[Ex->sgfn], cg->d_fgfs[Ey->sgfn], cg->d_fgfs[Ez->sgfn],
                        cg->d_fgfs[Bx->sgfn], cg->d_fgfs[By->sgfn], cg->d_fgfs[Bz->sgfn],
                        cg->d_fgfs[Kpsi->sgfn], cg->d_fgfs[Kphi->sgfn],
                        cg->d_fgfs[Ex1->sgfn], cg->d_fgfs[Ey1->sgfn], cg->d_fgfs[Ez1->sgfn],
                        cg->d_fgfs[Bx1->sgfn], cg->d_fgfs[By1->sgfn], cg->d_fgfs[Bz1->sgfn],
                        cg->d_fgfs[Kpsi1->sgfn], cg->d_fgfs[Kphi1->sgfn],
                        cg->d_fgfs[rho->sgfn],
                        cg->d_fgfs[Sx->sgfn], cg->d_fgfs[Sy->sgfn], cg->d_fgfs[Sz->sgfn],
                        cg->d_fgfs[Sxx->sgfn], cg->d_fgfs[Sxy->sgfn], cg->d_fgfs[Sxz->sgfn],
                        cg->d_fgfs[Syy->sgfn], cg->d_fgfs[Syz->sgfn], cg->d_fgfs[Szz->sgfn],
                        Symmetry, lev, ndeps
                    );

                    // BSSN RHS (corrector): state from SynchList_pre, RHS output to SynchList_cor
                    // (BSSN sym_bd already applied above, before EM RHS.)
                    gpu_compute_rhs_bssn_launch(
                        cg->stream,
                        cg->shape, cg->shape_pad, TRK4, cg->d_X[0], cg->d_X[1], cg->d_X[2],
                        cg->d_fgfs[phi->sgfn], cg->d_fgfs[trK->sgfn],
                        cg->d_fgfs[gxx->sgfn], cg->d_fgfs[gxy->sgfn], cg->d_fgfs[gxz->sgfn],
                        cg->d_fgfs[gyy->sgfn], cg->d_fgfs[gyz->sgfn], cg->d_fgfs[gzz->sgfn],
                        cg->d_fgfs[Axx->sgfn], cg->d_fgfs[Axy->sgfn], cg->d_fgfs[Axz->sgfn],
                        cg->d_fgfs[Ayy->sgfn], cg->d_fgfs[Ayz->sgfn], cg->d_fgfs[Azz->sgfn],
                        cg->d_fgfs[Gmx->sgfn], cg->d_fgfs[Gmy->sgfn], cg->d_fgfs[Gmz->sgfn],
                        cg->d_fgfs[Lap->sgfn],
                        cg->d_fgfs[Sfx->sgfn], cg->d_fgfs[Sfy->sgfn], cg->d_fgfs[Sfz->sgfn],
                        cg->d_fgfs[dtSfx->sgfn], cg->d_fgfs[dtSfy->sgfn], cg->d_fgfs[dtSfz->sgfn],
                        cg->d_fgfs[phi1->sgfn], cg->d_fgfs[trK1->sgfn],
                        cg->d_fgfs[gxx1->sgfn], cg->d_fgfs[gxy1->sgfn], cg->d_fgfs[gxz1->sgfn],
                        cg->d_fgfs[gyy1->sgfn], cg->d_fgfs[gyz1->sgfn], cg->d_fgfs[gzz1->sgfn],
                        cg->d_fgfs[Axx1->sgfn], cg->d_fgfs[Axy1->sgfn], cg->d_fgfs[Axz1->sgfn],
                        cg->d_fgfs[Ayy1->sgfn], cg->d_fgfs[Ayz1->sgfn], cg->d_fgfs[Azz1->sgfn],
                        cg->d_fgfs[Gmx1->sgfn], cg->d_fgfs[Gmy1->sgfn], cg->d_fgfs[Gmz1->sgfn],
                        cg->d_fgfs[Lap1->sgfn],
                        cg->d_fgfs[Sfx1->sgfn], cg->d_fgfs[Sfy1->sgfn], cg->d_fgfs[Sfz1->sgfn],
                        cg->d_fgfs[dtSfx1->sgfn], cg->d_fgfs[dtSfy1->sgfn], cg->d_fgfs[dtSfz1->sgfn],
                        cg->d_fgfs[rho->sgfn],
                        cg->d_fgfs[Sx->sgfn], cg->d_fgfs[Sy->sgfn], cg->d_fgfs[Sz->sgfn],
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
                        Symmetry, lev, ndeps, cor
                    );
#ifdef With_AHF
                    AH_Mark_RHS_Scratch(*cg);
#endif

                    // Fused RK4 + Sommerfeld for ALL variables (BSSN + EM)
                    MyList<var> *varl0 = StateList, *varl = SynchList_pre, *varl1 = SynchList_cor, *varlrhs = RHSList;
                    while (varl0) {
                        FusedRK4Vars rk4_vars;
                        rk4_vars.num_vars = 0;
                        FusedSommerfeldVars somm_vars;
                        somm_vars.num_vars = 0;
#if (SommerType == 1)
                        double *shib_f0[MAX_FUSED_VARS];
                        int shib_count = 0;
#endif

                        while (varl0 && rk4_vars.num_vars < MAX_FUSED_VARS) {
                            rk4_vars.f0[rk4_vars.num_vars] = cg->d_fgfs[varl0->data->sgfn];
                            rk4_vars.f1[rk4_vars.num_vars] = cg->d_fgfs[varl1->data->sgfn];
                            rk4_vars.f_rhs[rk4_vars.num_vars] = cg->d_fgfs[varlrhs->data->sgfn];

                            if (lev == 0) {
                                somm_vars.f0[somm_vars.num_vars] = cg->d_fgfs[varl->data->sgfn];
                                somm_vars.f_out[somm_vars.num_vars] = cg->d_fgfs[varl1->data->sgfn];
                            } else {
                                somm_vars.f0[somm_vars.num_vars] = cg->d_fgfs[varl0->data->sgfn];
                                somm_vars.f_out[somm_vars.num_vars] = cg->d_fgfs[varl1->data->sgfn];
                            }

                            somm_vars.propspeed[somm_vars.num_vars] = varl0->data->propspeed;
                            somm_vars.SoA0[somm_vars.num_vars] = varl0->data->SoA[0];
                            somm_vars.SoA1[somm_vars.num_vars] = varl0->data->SoA[1];
                            somm_vars.SoA2[somm_vars.num_vars] = varl0->data->SoA[2];
#if (SommerType == 1)
                            shib_f0[shib_count++] = cg->d_fgfs[varl->data->sgfn];
#endif

                            rk4_vars.num_vars++;
                            somm_vars.num_vars++;

                            varl0 = varl0->next;
                            varl = varl->next;
                            varl1 = varl1->next;
                            varlrhs = varlrhs->next;
                        }

#if (SommerType == 0)
#ifndef WithShell
                        if (lev == 0) {
                            gpu_fused_sommerfeld_routbam_launch(
                                cg->stream, cg->shape, cg->shape_pad, cg->d_X[0], cg->d_X[1], cg->d_X[2],
                                Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                                somm_vars, Symmetry
                            );
                        }
#endif
#endif

                        gpu_fused_rungekutta4_rout_launch(
                            cg->stream, cg->shape, cg->shape_pad, dT_lev,
                            rk4_vars, iter_count
                        );

#ifndef WithShell
                        if (lev > 0) {
#endif
                            gpu_fused_sommerfeld_rout_launch(
                                cg->stream, cg->shape, cg->shape_pad, cg->d_X[0], cg->d_X[1], cg->d_X[2],
                                Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                                dT_lev, cg->d_fgfs[phi0->sgfn], cg->d_fgfs[Lap0->sgfn],
                                somm_vars, Symmetry, cor
                            );
#ifndef WithShell
                        }
#endif

#if (SommerType == 1)
                        if (lev == 1) {
                            FusedSommerfeldVars somm_vars_shib;
                            somm_vars_shib.num_vars = somm_vars.num_vars;
                            for (int iv = 0; iv < somm_vars.num_vars; iv++) {
                                somm_vars_shib.f0[iv] = shib_f0[iv];
                                somm_vars_shib.f_out[iv] = somm_vars.f_out[iv];
                                somm_vars_shib.propspeed[iv] = somm_vars.propspeed[iv];
                                somm_vars_shib.SoA0[iv] = somm_vars.SoA0[iv];
                                somm_vars_shib.SoA1[iv] = somm_vars.SoA1[iv];
                                somm_vars_shib.SoA2[iv] = somm_vars.SoA2[iv];
                            }
                            gpu_fused_sommerfeld_rout_launch(
                                cg->stream, cg->shape, cg->shape_pad, cg->d_X[0], cg->d_X[1], cg->d_X[2],
                                Pp->data->bbox[0], Pp->data->bbox[1], Pp->data->bbox[2], Pp->data->bbox[3], Pp->data->bbox[4], Pp->data->bbox[5],
                                dT_lev, cg->d_fgfs[phi0->sgfn], cg->d_fgfs[Lap0->sgfn],
                                somm_vars_shib, Symmetry, cor
                            );
                        }
#endif
                    }

                    gpu_lowerboundset_launch(cg->stream, cg->shape, cg->shape_pad, cg->d_fgfs[phi1->sgfn], chitiny);
                }
                if (BP == Pp->data->ble) break;
                BP = BP->next;
            }
            Pp = Pp->next;
        }

#ifdef WithShell
        if (lev == 0)
        {
            int ERROR = 0;
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
                        gpu_enforce_ga_launch(
                            cg->stream, cg->shape, cg->shape_pad,
                            cg->d_fgfs[gxx->sgfn], cg->d_fgfs[gxy->sgfn], cg->d_fgfs[gxz->sgfn], cg->d_fgfs[gyy->sgfn], cg->d_fgfs[gyz->sgfn], cg->d_fgfs[gzz->sgfn],
                            cg->d_fgfs[Axx->sgfn], cg->d_fgfs[Axy->sgfn], cg->d_fgfs[Axz->sgfn], cg->d_fgfs[Ayy->sgfn], cg->d_fgfs[Ayz->sgfn], cg->d_fgfs[Azz->sgfn]
                        );
#elif (AGM == 1)
                        if (iter_count == 3)
                            gpu_enforce_ga_launch(
                                cg->stream, cg->shape, cg->shape_pad,
                                cg->d_fgfs[gxx->sgfn], cg->d_fgfs[gxy->sgfn], cg->d_fgfs[gxz->sgfn], cg->d_fgfs[gyy->sgfn], cg->d_fgfs[gyz->sgfn], cg->d_fgfs[gzz->sgfn],
                                cg->d_fgfs[Axx->sgfn], cg->d_fgfs[Axy->sgfn], cg->d_fgfs[Axz->sgfn], cg->d_fgfs[Ayy->sgfn], cg->d_fgfs[Ayz->sgfn], cg->d_fgfs[Azz->sgfn]
                            );
#endif

                        // [REMOVED] EM symmetry boundary - not in CPU version

                        // EM RHS (Shell corrector): shell coordinate transform + stress-energy
                        gpu_compute_rhs_empart_ss_launch(
                            cg->stream,
                            cg->shape, cg->shape_pad,
                            cg->d_fgfs[fngfs+ShellPatch::gx], cg->d_fgfs[fngfs+ShellPatch::gy], cg->d_fgfs[fngfs+ShellPatch::gz],
                            cg->d_fgfs[fngfs+ShellPatch::drhodx], cg->d_fgfs[fngfs+ShellPatch::drhody], cg->d_fgfs[fngfs+ShellPatch::drhodz],
                            cg->d_fgfs[fngfs+ShellPatch::dsigmadx], cg->d_fgfs[fngfs+ShellPatch::dsigmady], cg->d_fgfs[fngfs+ShellPatch::dsigmadz],
                            cg->d_fgfs[fngfs+ShellPatch::dRdx], cg->d_fgfs[fngfs+ShellPatch::dRdy], cg->d_fgfs[fngfs+ShellPatch::dRdz],
                            cg->d_fgfs[phi->sgfn],
                            cg->d_fgfs[Jx->sgfn], cg->d_fgfs[Jy->sgfn], cg->d_fgfs[Jz->sgfn], cg->d_fgfs[qchar->sgfn],
                            cg->d_fgfs[gxx->sgfn], cg->d_fgfs[gxy->sgfn], cg->d_fgfs[gxz->sgfn],
                            cg->d_fgfs[gyy->sgfn], cg->d_fgfs[gyz->sgfn], cg->d_fgfs[gzz->sgfn],
                            cg->d_fgfs[Lap->sgfn],
                            cg->d_fgfs[Sfx->sgfn], cg->d_fgfs[Sfy->sgfn], cg->d_fgfs[Sfz->sgfn],
                            cg->d_fgfs[trK->sgfn],
                            cg->d_fgfs[Ex->sgfn], cg->d_fgfs[Ey->sgfn], cg->d_fgfs[Ez->sgfn],
                            cg->d_fgfs[Bx->sgfn], cg->d_fgfs[By->sgfn], cg->d_fgfs[Bz->sgfn],
                            cg->d_fgfs[Kpsi->sgfn], cg->d_fgfs[Kphi->sgfn],
                            cg->d_fgfs[Ex1->sgfn], cg->d_fgfs[Ey1->sgfn], cg->d_fgfs[Ez1->sgfn],
                            cg->d_fgfs[Bx1->sgfn], cg->d_fgfs[By1->sgfn], cg->d_fgfs[Bz1->sgfn],
                            cg->d_fgfs[Kpsi1->sgfn], cg->d_fgfs[Kphi1->sgfn],
                            cg->d_fgfs[rho->sgfn],
                            cg->d_fgfs[Sx->sgfn], cg->d_fgfs[Sy->sgfn], cg->d_fgfs[Sz->sgfn],
                            cg->d_fgfs[Sxx->sgfn], cg->d_fgfs[Sxy->sgfn], cg->d_fgfs[Sxz->sgfn],
                            cg->d_fgfs[Syy->sgfn], cg->d_fgfs[Syz->sgfn], cg->d_fgfs[Szz->sgfn],
                            Symmetry, lev, numepsh, sPp->data->sst
                        );

                        gpu_compute_rhs_bssn_ss_launch(
                            cg->stream, cg->shape, cg->shape_pad, TRK4,
                            cg->d_X[0], cg->d_X[1], cg->d_X[2],
                            cg->d_fgfs[fngfs+ShellPatch::gx], cg->d_fgfs[fngfs+ShellPatch::gy], cg->d_fgfs[fngfs+ShellPatch::gz],
                            cg->d_fgfs[fngfs+ShellPatch::drhodx], cg->d_fgfs[fngfs+ShellPatch::drhody], cg->d_fgfs[fngfs+ShellPatch::drhodz],
                            cg->d_fgfs[fngfs+ShellPatch::dsigmadx], cg->d_fgfs[fngfs+ShellPatch::dsigmady], cg->d_fgfs[fngfs+ShellPatch::dsigmadz],
                            cg->d_fgfs[fngfs+ShellPatch::dRdx], cg->d_fgfs[fngfs+ShellPatch::dRdy], cg->d_fgfs[fngfs+ShellPatch::dRdz],
                            cg->d_fgfs[fngfs+ShellPatch::drhodxx], cg->d_fgfs[fngfs+ShellPatch::drhodxy], cg->d_fgfs[fngfs+ShellPatch::drhodxz],
                            cg->d_fgfs[fngfs+ShellPatch::drhodyy], cg->d_fgfs[fngfs+ShellPatch::drhodyz], cg->d_fgfs[fngfs+ShellPatch::drhodzz],
                            cg->d_fgfs[fngfs+ShellPatch::dsigmadxx], cg->d_fgfs[fngfs+ShellPatch::dsigmadxy], cg->d_fgfs[fngfs+ShellPatch::dsigmadxz],
                            cg->d_fgfs[fngfs+ShellPatch::dsigmadyy], cg->d_fgfs[fngfs+ShellPatch::dsigmadyz], cg->d_fgfs[fngfs+ShellPatch::dsigmadzz],
                            cg->d_fgfs[fngfs+ShellPatch::dRdxx], cg->d_fgfs[fngfs+ShellPatch::dRdxy], cg->d_fgfs[fngfs+ShellPatch::dRdxz],
                            cg->d_fgfs[fngfs+ShellPatch::dRdyy], cg->d_fgfs[fngfs+ShellPatch::dRdyz], cg->d_fgfs[fngfs+ShellPatch::dRdzz],
                            cg->d_fgfs[phi->sgfn], cg->d_fgfs[trK->sgfn],
                            cg->d_fgfs[gxx->sgfn], cg->d_fgfs[gxy->sgfn], cg->d_fgfs[gxz->sgfn],
                            cg->d_fgfs[gyy->sgfn], cg->d_fgfs[gyz->sgfn], cg->d_fgfs[gzz->sgfn],
                            cg->d_fgfs[Axx->sgfn], cg->d_fgfs[Axy->sgfn], cg->d_fgfs[Axz->sgfn],
                            cg->d_fgfs[Ayy->sgfn], cg->d_fgfs[Ayz->sgfn], cg->d_fgfs[Azz->sgfn],
                            cg->d_fgfs[Gmx->sgfn], cg->d_fgfs[Gmy->sgfn], cg->d_fgfs[Gmz->sgfn],
                            cg->d_fgfs[Lap->sgfn],
                            cg->d_fgfs[Sfx->sgfn], cg->d_fgfs[Sfy->sgfn], cg->d_fgfs[Sfz->sgfn],
                            cg->d_fgfs[dtSfx->sgfn], cg->d_fgfs[dtSfy->sgfn], cg->d_fgfs[dtSfz->sgfn],
                            cg->d_fgfs[phi1->sgfn], cg->d_fgfs[trK1->sgfn],
                            cg->d_fgfs[gxx1->sgfn], cg->d_fgfs[gxy1->sgfn], cg->d_fgfs[gxz1->sgfn],
                            cg->d_fgfs[gyy1->sgfn], cg->d_fgfs[gyz1->sgfn], cg->d_fgfs[gzz1->sgfn],
                            cg->d_fgfs[Axx1->sgfn], cg->d_fgfs[Axy1->sgfn], cg->d_fgfs[Axz1->sgfn],
                            cg->d_fgfs[Ayy1->sgfn], cg->d_fgfs[Ayz1->sgfn], cg->d_fgfs[Azz1->sgfn],
                            cg->d_fgfs[Gmx1->sgfn], cg->d_fgfs[Gmy1->sgfn], cg->d_fgfs[Gmz1->sgfn],
                            cg->d_fgfs[Lap1->sgfn],
                            cg->d_fgfs[Sfx1->sgfn], cg->d_fgfs[Sfy1->sgfn], cg->d_fgfs[Sfz1->sgfn],
                            cg->d_fgfs[dtSfx1->sgfn], cg->d_fgfs[dtSfy1->sgfn], cg->d_fgfs[dtSfz1->sgfn],
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
                            cg->d_fgfs[Cons_Ham->sgfn], cg->d_fgfs[Cons_Px->sgfn], cg->d_fgfs[Cons_Py->sgfn], cg->d_fgfs[Cons_Pz->sgfn],
                            cg->d_fgfs[Cons_Gx->sgfn], cg->d_fgfs[Cons_Gy->sgfn], cg->d_fgfs[Cons_Gz->sgfn],
                            Symmetry, lev, numepsh, sPp->data->sst, cor
                        );
#ifdef With_AHF
                        AH_Mark_RHS_Scratch(*cg);
#endif
                        cudaError_t err = cudaGetLastError();
                        if (err != cudaSuccess) {
                            cout << "GPU SS RHS error: " << cudaGetErrorString(err) << " sst=" << sPp->data->sst << endl;
                            ERROR = 1;
                        }

                        {
                            MyList<var> *varl0 = StateList, *varl = SynchList_pre, *varl1 = SynchList_cor, *varlrhs = RHSList;
                            while (varl0) {
                                FusedRK4Vars rk4_vars; rk4_vars.num_vars = 0;
                                FusedSommerfeldVars somm_vars; somm_vars.num_vars = 0;

                                while (varl0 && rk4_vars.num_vars < MAX_FUSED_VARS) {
                                    rk4_vars.f0[rk4_vars.num_vars] = cg->d_fgfs[varl0->data->sgfn];
                                    rk4_vars.f1[rk4_vars.num_vars] = cg->d_fgfs[varl1->data->sgfn];
                                    rk4_vars.f_rhs[rk4_vars.num_vars] = cg->d_fgfs[varlrhs->data->sgfn];

                                    somm_vars.f0[somm_vars.num_vars] = cg->d_fgfs[varl->data->sgfn];
                                    somm_vars.f_out[somm_vars.num_vars] = cg->d_fgfs[varl1->data->sgfn];
                                    somm_vars.propspeed[somm_vars.num_vars] = varl0->data->propspeed;
                                    somm_vars.SoA0[somm_vars.num_vars] = varl0->data->SoA[0];
                                    somm_vars.SoA1[somm_vars.num_vars] = varl0->data->SoA[1];
                                    somm_vars.SoA2[somm_vars.num_vars] = varl0->data->SoA[2];

                                    rk4_vars.num_vars++;
                                    somm_vars.num_vars++;

                                    varl0 = varl0->next; varl = varl->next; varl1 = varl1->next; varlrhs = varlrhs->next;
                                }

#if (SommerType == 0)
                                gpu_fused_sommerfeld_routbam_ss_launch(
                                    cg->stream, cg->shape, cg->shape_pad, cg->d_X[0], cg->d_X[1], cg->d_X[2],
                                    sPp->data->bbox[0], sPp->data->bbox[1], sPp->data->bbox[2], sPp->data->bbox[3], sPp->data->bbox[4], sPp->data->bbox[5],
                                    somm_vars, Symmetry
                                );
#endif

                                gpu_fused_rungekutta4_rout_launch(cg->stream, cg->shape, cg->shape_pad, dT_lev, rk4_vars, iter_count);
                            }
                        }
                        // EM field falloff boundary on Shell patch (corrector)
                        gpu_falloff_ss_launch(cg->stream,
                            cg->shape[0], cg->shape[1], cg->shape[2],
                            cg->shape_pad[0], cg->shape_pad[1],
                            cg->d_fgfs[Ex1->sgfn], cg->d_fgfs[fngfs+ShellPatch::gz], 2, sPp->data->bbox[5]);
                        gpu_falloff_ss_launch(cg->stream,
                            cg->shape[0], cg->shape[1], cg->shape[2],
                            cg->shape_pad[0], cg->shape_pad[1],
                            cg->d_fgfs[Ey1->sgfn], cg->d_fgfs[fngfs+ShellPatch::gz], 2, sPp->data->bbox[5]);
                        gpu_falloff_ss_launch(cg->stream,
                            cg->shape[0], cg->shape[1], cg->shape[2],
                            cg->shape_pad[0], cg->shape_pad[1],
                            cg->d_fgfs[Ez1->sgfn], cg->d_fgfs[fngfs+ShellPatch::gz], 2, sPp->data->bbox[5]);
                        gpu_falloff_ss_launch(cg->stream,
                            cg->shape[0], cg->shape[1], cg->shape[2],
                            cg->shape_pad[0], cg->shape_pad[1],
                            cg->d_fgfs[Bx1->sgfn], cg->d_fgfs[fngfs+ShellPatch::gz], 2, sPp->data->bbox[5]);
                        gpu_falloff_ss_launch(cg->stream,
                            cg->shape[0], cg->shape[1], cg->shape[2],
                            cg->shape_pad[0], cg->shape_pad[1],
                            cg->d_fgfs[By1->sgfn], cg->d_fgfs[fngfs+ShellPatch::gz], 2, sPp->data->bbox[5]);
                        gpu_falloff_ss_launch(cg->stream,
                            cg->shape[0], cg->shape[1], cg->shape[2],
                            cg->shape_pad[0], cg->shape_pad[1],
                            cg->d_fgfs[Bz1->sgfn], cg->d_fgfs[fngfs+ShellPatch::gz], 2, sPp->data->bbox[5]);
                        gpu_falloff_ss_launch(cg->stream,
                            cg->shape[0], cg->shape[1], cg->shape[2],
                            cg->shape_pad[0], cg->shape_pad[1],
                            cg->d_fgfs[Kpsi1->sgfn], cg->d_fgfs[fngfs+ShellPatch::gz], 3, sPp->data->bbox[5]);
                        gpu_falloff_ss_launch(cg->stream,
                            cg->shape[0], cg->shape[1], cg->shape[2],
                            cg->shape_pad[0], cg->shape_pad[1],
                            cg->d_fgfs[Kphi1->sgfn], cg->d_fgfs[fngfs+ShellPatch::gz], 3, sPp->data->bbox[5]);
                        gpu_lowerboundset_launch(cg->stream, cg->shape, cg->shape_pad, cg->d_fgfs[phi1->sgfn], chitiny);
                    }
                    if (BP == sPp->data->ble) break;
                    BP = BP->next;
                }
                sPp = sPp->next;
            }
            int erh = ERROR;
            MPI_Allreduce(&erh, &ERROR, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
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
        }
#endif

        Parallel::Sync_GPU(GH->PatL[lev], SynchList_cor, Symmetry);

#ifdef WithShell
        if (lev == 0)
        {
            double prev_clock, curr_clock;
            GPUManager::getInstance().synchronize_all();
            if (myrank == 0) prev_clock = MPI_Wtime();
            SH->Synch_GPU(SynchList_cor, Symmetry);
            if (myrank == 0)
            {
                curr_clock = MPI_Wtime();
                cout << "Shell stuff synchronization used " << (curr_clock - prev_clock) << " seconds! (RK substep " << iter_count << ")" << endl;
                if (iter_count == 3) cout << endl;
            }
        }
#endif

        // for black hole position — RK4 corrector (unconditional for EM)
        if (BH_num > 0 && lev == GH->levels - 1)
        {
            compute_Porg_rhs(Porg, Porg1, Sfx, Sfy, Sfz, lev);
            for (int ithBH = 0; ithBH < BH_num; ithBH++)
            {
                f_rungekutta4_scalar(dT_lev, Porg0[ithBH][0], Porg1[ithBH][0], Porg_rhs[ithBH][0], iter_count);
                f_rungekutta4_scalar(dT_lev, Porg0[ithBH][1], Porg1[ithBH][1], Porg_rhs[ithBH][1], iter_count);
                f_rungekutta4_scalar(dT_lev, Porg0[ithBH][2], Porg1[ithBH][2], Porg_rhs[ithBH][2], iter_count);
                if (Symmetry > 0) Porg1[ithBH][2] = fabs(Porg1[ithBH][2]);
                if (Symmetry == 2) { Porg1[ithBH][0] = fabs(Porg1[ithBH][0]); Porg1[ithBH][1] = fabs(Porg1[ithBH][1]); }
                if (!finite(Porg1[ithBH][0]) || !finite(Porg1[ithBH][1]) || !finite(Porg1[ithBH][2]))
                {
                    if (ErrorMonitor->outfile)
                        ErrorMonitor->outfile << iter_count << " corrector step finds NaN for BH's position from ("
                                              << Porg[ithBH][0] << "," << Porg[ithBH][1] << "," << Porg[ithBH][2] << ")" << endl;

                    MyList<var> *DG_List = new MyList<var>(Sfx0);
                    DG_List->insert(Sfx0); DG_List->insert(Sfy0); DG_List->insert(Sfz0);
                    Parallel::Dump_Data(GH->PatL[lev], DG_List, 0, PhysTime, dT_lev);
                    DG_List->clearList();
                }
            }
        }

        if (iter_count < 3) {
            Pp = GH->PatL[lev];
            while (Pp) {
                MyList<Block> *BP = Pp->data->blb;
                while (BP) {
                    Block *cg = BP->data;
                    cg->swapList(SynchList_pre, SynchList_cor, myrank);
                    if (BP == Pp->data->ble) break;
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
                        if (BP == sPp->data->ble) break;
                        BP = BP->next;
                    }
                    sPp = sPp->next;
                }
            }
#endif

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
    RestrictProlong(lev, YN, BB);

#ifdef WithShell
    if (lev == 0)
    {
        double prev_clock, curr_clock;
        if (myrank == 0) curr_clock = MPI_Wtime();
        SH->CS_Inter_GPU(SynchList_cor, Symmetry);
        if (myrank == 0)
        {
            prev_clock = curr_clock;
            curr_clock = MPI_Wtime();
            cout << "CS_Inter used " << (curr_clock - prev_clock) << " seconds!" << endl;
        }
    }
#endif
#endif

    Pp = GH->PatL[lev];
    while (Pp)
    {
        MyList<Block> *BP = Pp->data->blb;
        while (BP)
        {
            Block *cg = BP->data;
            cg->swapList(StateList, SynchList_cor, myrank);
            cg->swapList(OldStateList, SynchList_cor, myrank);
            if (BP == Pp->data->ble) break;
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
                if (BP == sPp->data->ble) break;
                BP = BP->next;
            }
            sPp = sPp->next;
        }
    }
#endif

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
// 该成员函数用于计算电磁波 Phi2
// Compute_Phi2_GPU — GPU NP scalar extraction
//================================================================================================

void bssnEM_class::Compute_Phi2_GPU(int lev)
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
                gpu_getnpem2_launch(
                    cg->stream, cg->shape, cg->shape_pad,
                    cg->d_X[0], cg->d_X[1], cg->d_X[2],
                    cg->d_fgfs[phi0->sgfn],
                    cg->d_fgfs[gxx0->sgfn], cg->d_fgfs[gxy0->sgfn], cg->d_fgfs[gxz0->sgfn],
                    cg->d_fgfs[gyy0->sgfn], cg->d_fgfs[gyz0->sgfn], cg->d_fgfs[gzz0->sgfn],
                    cg->d_fgfs[Ex0->sgfn], cg->d_fgfs[Ey0->sgfn], cg->d_fgfs[Ez0->sgfn],
                    cg->d_fgfs[Bx0->sgfn], cg->d_fgfs[By0->sgfn], cg->d_fgfs[Bz0->sgfn],
                    cg->d_fgfs[Rphi2->sgfn], cg->d_fgfs[Iphi2->sgfn],
                    Symmetry);

                gpu_getnpem1_launch(
                    cg->stream, cg->shape, cg->shape_pad,
                    cg->d_X[0], cg->d_X[1], cg->d_X[2],
                    cg->d_fgfs[phi0->sgfn],
                    cg->d_fgfs[gxx0->sgfn], cg->d_fgfs[gxy0->sgfn], cg->d_fgfs[gxz0->sgfn],
                    cg->d_fgfs[gyy0->sgfn], cg->d_fgfs[gyz0->sgfn], cg->d_fgfs[gzz0->sgfn],
                    cg->d_fgfs[Ex0->sgfn], cg->d_fgfs[Ey0->sgfn], cg->d_fgfs[Ez0->sgfn],
                    cg->d_fgfs[Bx0->sgfn], cg->d_fgfs[By0->sgfn], cg->d_fgfs[Bz0->sgfn],
                    cg->d_fgfs[Rphi1->sgfn], cg->d_fgfs[Iphi1->sgfn],
                    Symmetry);
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
            MyList<Block> *BL = sPp->data->blb;
            while (BL)
            {
                Block *cg = BL->data;
                if (myrank == cg->rank)
                {
                    gpu_getnpem2_ss_launch(
                        cg->stream, cg->shape, cg->shape_pad,
                        cg->d_X[0], cg->d_X[1], cg->d_X[2],
                        cg->d_fgfs[phi0->sgfn],
                        cg->d_fgfs[gxx0->sgfn], cg->d_fgfs[gxy0->sgfn], cg->d_fgfs[gxz0->sgfn],
                        cg->d_fgfs[gyy0->sgfn], cg->d_fgfs[gyz0->sgfn], cg->d_fgfs[gzz0->sgfn],
                        cg->d_fgfs[Ex0->sgfn], cg->d_fgfs[Ey0->sgfn], cg->d_fgfs[Ez0->sgfn],
                        cg->d_fgfs[Bx0->sgfn], cg->d_fgfs[By0->sgfn], cg->d_fgfs[Bz0->sgfn],
                        cg->d_fgfs[Rphi2->sgfn], cg->d_fgfs[Iphi2->sgfn],
                        Symmetry);

                    gpu_getnpem1_ss_launch(
                        cg->stream, cg->shape, cg->shape_pad,
                        cg->d_X[0], cg->d_X[1], cg->d_X[2],
                        cg->d_fgfs[phi0->sgfn],
                        cg->d_fgfs[gxx0->sgfn], cg->d_fgfs[gxy0->sgfn], cg->d_fgfs[gxz0->sgfn],
                        cg->d_fgfs[gyy0->sgfn], cg->d_fgfs[gyz0->sgfn], cg->d_fgfs[gzz0->sgfn],
                        cg->d_fgfs[Ex0->sgfn], cg->d_fgfs[Ey0->sgfn], cg->d_fgfs[Ez0->sgfn],
                        cg->d_fgfs[Bx0->sgfn], cg->d_fgfs[By0->sgfn], cg->d_fgfs[Bz0->sgfn],
                        cg->d_fgfs[Rphi1->sgfn], cg->d_fgfs[Iphi1->sgfn],
                        Symmetry);
                }
                if (BL == sPp->data->ble)
                    break;
                BL = BL->next;
            }
            sPp = sPp->next;
        }
    }
#endif

    MyList<var> *DG_List = new MyList<var>(Rphi2);
    DG_List->insert(Iphi2);
    DG_List->insert(Rphi1);
    DG_List->insert(Iphi1);
    Parallel::Sync_GPU(GH->PatL[lev], DG_List, Symmetry);
#ifdef WithShell
    if (lev == 0) SH->Synch_GPU(DG_List, Symmetry);
#endif
    DG_List->clearList();
}

//================================================================================================



//================================================================================================
// 该成员函数用于分析电磁场的数据
// AnalysisStuff_EM — GPU-aware analysis (delegates to Compute_Phi2_GPU)
//================================================================================================

void bssnEM_class::AnalysisStuff_EM(int lev, double dT_lev)
{
    LastAnas += dT_lev;

    if (LastAnas >= AnasTime)
    {
        Compute_Phi2_GPU(lev);
        double *RP, *IP;
        int NN = 0;
        for (int pl = 1; pl < maxl + 1; pl++)
            for (int pm = -pl; pm < pl + 1; pm++)
                NN++;
        RP = new double[NN];
        IP = new double[NN];
        double Rex = maxrex;
        for (int i = 0; i < decn; i++)
        {
#ifdef WithShell
            if (lev > 0 || Rex < GH->bbox[0][0][3])
            {
                Waveshell->gpu_surf_Wave(Rex, lev, GH,
                                     Ex0, Ey0, Ez0, Bx0, By0, Bz0, phi0,
                                     gxx0, gxy0, gxz0, gyy0, gyz0, gzz0,
                                     1, maxl, NN, RP, IP, ErrorMonitor,
                                     f_getnpem2_point);
            }
            else
            {
                Waveshell->gpu_surf_Wave(Rex, lev, SH,
                                     Ex0, Ey0, Ez0, Bx0, By0, Bz0, phi0,
                                     gxx0, gxy0, gxz0, gyy0, gyz0, gzz0,
                                     1, maxl, NN, RP, IP, ErrorMonitor,
                                     f_getnpem2_point);
            }
#else
            Waveshell->gpu_surf_Wave(Rex, lev, GH,
                                 Ex0, Ey0, Ez0, Bx0, By0, Bz0, phi0,
                                 gxx0, gxy0, gxz0, gyy0, gyz0, gzz0,
                                 1, maxl, NN, RP, IP, ErrorMonitor,
                                 f_getnpem2_point);
#endif
            Phi2Monitor->writefile(PhysTime, NN, RP, IP);
            Rex = Rex - drex;
        }
        delete[] RP;
        delete[] IP;

        NN = 0;
        for (int pl = 0; pl < maxl + 1; pl++)
            for (int pm = -pl; pm < pl + 1; pm++)
                NN++;
        RP = new double[NN];
        IP = new double[NN];
        Rex = maxrex;
        for (int i = 0; i < decn; i++)
        {
#ifdef WithShell
            if (lev > 0 || Rex < GH->bbox[0][0][3])
            {
                Waveshell->gpu_surf_Wave(Rex, lev, GH,
                                     Ex0, Ey0, Ez0, Bx0, By0, Bz0, phi0,
                                     gxx0, gxy0, gxz0, gyy0, gyz0, gzz0,
                                     0, maxl, NN, RP, IP, ErrorMonitor,
                                     f_getnpem1_point);
            }
            else
            {
                Waveshell->gpu_surf_Wave(Rex, lev, SH,
                                     Ex0, Ey0, Ez0, Bx0, By0, Bz0, phi0,
                                     gxx0, gxy0, gxz0, gyy0, gyz0, gzz0,
                                     0, maxl, NN, RP, IP, ErrorMonitor,
                                     f_getnpem1_point);
            }
#else
            Waveshell->gpu_surf_Wave(Rex, lev, GH,
                                 Ex0, Ey0, Ez0, Bx0, By0, Bz0, phi0,
                                 gxx0, gxy0, gxz0, gyy0, gyz0, gzz0,
                                 0, maxl, NN, RP, IP, ErrorMonitor,
                                 f_getnpem1_point);
#endif
            Phi1Monitor->writefile(PhysTime, NN, RP, IP);
            Rex = Rex - drex;
        }
        delete[] RP;
        delete[] IP;
    }

    AnalysisStuff(lev, dT_lev);

    LastAnas = 0;
}

//================================================================================================



//================================================================================================
// 该成员函数用于对约束数据进行插值
// Interp_Constraint — GPU RHS for EM + BSSN, CPU interpolation output
//================================================================================================

void bssnEM_class::Interp_Constraint()
{
    if (a_lev > 0) return;

    for (int lev = 0; lev < GH->levels; lev++)
    {
        if (lev > 0)
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
                        // SMEM-based EM RHS needs EM ghost zones pre-mirrored.
                        gpu_apply_symmetry_bd_empart_launch(
                            cg->stream, cg->shape[0], cg->shape[1], cg->shape[2],
                            cg->shape_pad[0], cg->shape_pad[1], cg->shape_pad[2],
                            cg->d_fgfs[Ex0->sgfn], cg->d_fgfs[Ey0->sgfn], cg->d_fgfs[Ez0->sgfn],
                            cg->d_fgfs[Bx0->sgfn], cg->d_fgfs[By0->sgfn], cg->d_fgfs[Bz0->sgfn],
                            cg->d_fgfs[Kpsi0->sgfn], cg->d_fgfs[Kphi0->sgfn]
                        );

                        gpu_compute_rhs_empart_launch(
                            cg->stream,
                            cg->shape, cg->shape_pad, cg->d_X[0], cg->d_X[1], cg->d_X[2],
                            cg->d_fgfs[phi0->sgfn],
                            cg->d_fgfs[Jx->sgfn], cg->d_fgfs[Jy->sgfn], cg->d_fgfs[Jz->sgfn], cg->d_fgfs[qchar->sgfn],
                            cg->d_fgfs[gxx0->sgfn], cg->d_fgfs[gxy0->sgfn], cg->d_fgfs[gxz0->sgfn],
                            cg->d_fgfs[gyy0->sgfn], cg->d_fgfs[gyz0->sgfn], cg->d_fgfs[gzz0->sgfn],
                            cg->d_fgfs[Lap0->sgfn],
                            cg->d_fgfs[Sfx0->sgfn], cg->d_fgfs[Sfy0->sgfn], cg->d_fgfs[Sfz0->sgfn],
                            cg->d_fgfs[trK0->sgfn],
                            cg->d_fgfs[Ex0->sgfn], cg->d_fgfs[Ey0->sgfn], cg->d_fgfs[Ez0->sgfn],
                            cg->d_fgfs[Bx0->sgfn], cg->d_fgfs[By0->sgfn], cg->d_fgfs[Bz0->sgfn],
                            cg->d_fgfs[Kpsi0->sgfn], cg->d_fgfs[Kphi0->sgfn],
                            cg->d_fgfs[Ex_rhs->sgfn], cg->d_fgfs[Ey_rhs->sgfn], cg->d_fgfs[Ez_rhs->sgfn],
                            cg->d_fgfs[Bx_rhs->sgfn], cg->d_fgfs[By_rhs->sgfn], cg->d_fgfs[Bz_rhs->sgfn],
                            cg->d_fgfs[Kpsi_rhs->sgfn], cg->d_fgfs[Kphi_rhs->sgfn],
                            cg->d_fgfs[rho->sgfn],
                            cg->d_fgfs[Sx->sgfn], cg->d_fgfs[Sy->sgfn], cg->d_fgfs[Sz->sgfn],
                            cg->d_fgfs[Sxx->sgfn], cg->d_fgfs[Sxy->sgfn], cg->d_fgfs[Sxz->sgfn],
                            cg->d_fgfs[Syy->sgfn], cg->d_fgfs[Syz->sgfn], cg->d_fgfs[Szz->sgfn],
                            Symmetry, lev, ndeps);

                        gpu_compute_rhs_bssn_launch(
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
                            cg->d_fgfs[phi_rhs->sgfn], cg->d_fgfs[trK_rhs->sgfn],
                            cg->d_fgfs[gxx_rhs->sgfn], cg->d_fgfs[gxy_rhs->sgfn], cg->d_fgfs[gxz_rhs->sgfn],
                            cg->d_fgfs[gyy_rhs->sgfn], cg->d_fgfs[gyz_rhs->sgfn], cg->d_fgfs[gzz_rhs->sgfn],
                            cg->d_fgfs[Axx_rhs->sgfn], cg->d_fgfs[Axy_rhs->sgfn], cg->d_fgfs[Axz_rhs->sgfn],
                            cg->d_fgfs[Ayy_rhs->sgfn], cg->d_fgfs[Ayz_rhs->sgfn], cg->d_fgfs[Azz_rhs->sgfn],
                            cg->d_fgfs[Gmx_rhs->sgfn], cg->d_fgfs[Gmy_rhs->sgfn], cg->d_fgfs[Gmz_rhs->sgfn],
                            cg->d_fgfs[Lap_rhs->sgfn],
                            cg->d_fgfs[Sfx_rhs->sgfn], cg->d_fgfs[Sfy_rhs->sgfn], cg->d_fgfs[Sfz_rhs->sgfn],
                            cg->d_fgfs[dtSfx_rhs->sgfn], cg->d_fgfs[dtSfy_rhs->sgfn], cg->d_fgfs[dtSfz_rhs->sgfn],
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
                            Symmetry, lev, ndeps, pre);
#ifdef With_AHF
                        AH_Mark_RHS_Scratch(*cg);
#endif
                    }
                    if (BP == Pp->data->ble) break;
                    BP = BP->next;
                }
                Pp = Pp->next;
            }
        }
        Parallel::Sync_GPU(GH->PatL[lev], ConstraintList, Symmetry);
    }
#ifdef WithShell
    SH->Synch_GPU(ConstraintList, Symmetry);
#endif

#warning "Interp_Constraint uses CPU path for interpolation output (GH->Interp_One_Point)"

    // Interpolation along x-axis — CPU path for output
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
    while (varl) { InList++; varl = varl->next; }
    double *shellf;
    shellf = new double[n * InList];
    for (int i = 0; i < n; i++)
    {
        double XX[3];
        XX[0] = x1[i]; XX[1] = y1[i]; XX[2] = z1[i];
        bool fg = GH->Interp_One_Point(ConstraintList, XX, shellf + i * InList, Symmetry);
#ifdef WithShell
        if (!fg) fg = SH->Interp_One_Point(ConstraintList, XX, shellf + i * InList, Symmetry);
#endif
        if (!fg && myrank == 0)
        {
            cout << "bssn_class::Interp_Constraint meets wrong" << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    ofstream outfile;
    char filename[50];
    sprintf(filename, "%s/interp_constraint_%05d.dat", ErrorMonitor->out_dir.c_str(), int(PhysTime / dT + 0.5));
    outfile.open(filename);
    outfile << "#  corrdinate, H_Res, Px_Res, Py_Res, Pz_Res, Gx_Res, Gy_Res, Gz_Res, ...." << endl;
    for (int i = 0; i < n; i++)
    {
        outfile << setw(10) << setprecision(10) << y1[i];
        for (int j = 0; j < InList; j++)
            outfile << " " << setw(16) << setprecision(15) << shellf[InList * i + j];
        outfile << endl;
    }

    delete[] shellf;
    delete[] x1; delete[] y1; delete[] z1;
}

//================================================================================================
