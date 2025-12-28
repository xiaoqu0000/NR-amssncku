
#ifdef newc
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <map>
using namespace std;
#else
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <map.h>
#endif

#include <mpi.h>

#include "macrodef.h"
#include "misc.h"
#include "cgh.h"
#include "Parallel.h"
#include "parameters.h"

//================================================================================================

// 定义 cgh class

//================================================================================================

cgh::cgh(int ingfsi, int fngfsi, int Symmetry, char *filename, int checkrun,
         monitor *ErrorMonitor) : ingfs(ingfsi), fngfs(fngfsi), trfls(0)
{
#if (PSTR == 1 || PSTR == 2 || PSTR == 3)
  Commlev = 0;
  start_rank = 0;
  end_rank = 0;
#endif

  if (!checkrun)
  {
    read_bbox(Symmetry, filename);
    sethandle(ErrorMonitor);
    for (int lev = 0; lev < levels; lev++)
      PatL[lev] = construct_patchlist(lev, Symmetry);
  }
}

//================================================================================================



//================================================================================================

// 该成员函数为析构函数，用于删除变量

//================================================================================================

cgh::~cgh()
{
  for (int lev = 0; lev < levels; lev++)
  {
    for (int grd = 0; grd < grids[lev]; grd++)
    {
      delete[] bbox[lev][grd];
      delete[] shape[lev][grd];
      delete[] handle[lev][grd];
    }
    delete[] bbox[lev];
    delete[] shape[lev];
    delete[] handle[lev];
    Parallel::KillBlocks(PatL[lev]);
    PatL[lev]->destroyList();
#if (RPB == 1)
    Parallel::destroypsuList_bam(bdsul[lev]);
    Parallel::destroypsuList_bam(rsul[lev]);
#endif
  }
  delete[] grids;
  delete[] Lt;
  delete[] bbox;
  delete[] shape;
  delete[] handle;
  delete[] PatL;
#if (RPB == 1)
  delete[] bdsul;
  delete[] rsul;
#endif

#if (PSTR == 1 || PSTR == 2 || PSTR == 3)
  for (int lev = 0; lev < levels; lev++)
  {
    MPI_Comm_free(&Commlev[lev]);
  }

  if (Commlev)
    delete[] Commlev;
  if (start_rank)
    delete[] start_rank;
  if (end_rank)
    delete[] end_rank;
#endif
  for (int lev = 0; lev < levels; lev++)
  {
    for (int ibh = 0; ibh < BH_num_in; ibh++)
      delete[] Porgls[lev][ibh];
    delete[] Porgls[lev];
  }
  delete[] Porgls;
}

//================================================================================================


//================================================================================================

// 该成员函数用于构建网格

//================================================================================================

#if (PSTR == 0)
void cgh::compose_cgh(int nprocs)
{
  for (int lev = 0; lev < levels; lev++)
  {
    checkPatchList(PatL[lev], false);
    Parallel::distribute(PatL[lev], nprocs, ingfs, fngfs, false);
#if (RPB == 1)
    // we need distributed box of PatL[lev] and PatL[lev-1]
    if (lev > 0)
    {
      Parallel::Constr_pointstr_OutBdLow2Hi(PatL[lev], PatL[lev - 1], bdsul[lev]);
      Parallel::Constr_pointstr_Restrict(PatL[lev], PatL[lev - 1], rsul[lev]);
    }
    else
    {
      bdsul[lev] = 0;
      rsul[lev] = 0;
    }
#endif
  }
}

//================================================================================================


//================================================================================================

// 该成员函数用于构建网格
// 针对 PSTR == 1 和 PSTR == 2 的情形

//================================================================================================

#elif (PSTR == 1 || PSTR == 2)
void cgh::compose_cgh(int nprocs)
{
  Commlev = new MPI_Comm[levels];
  construct_mylev(nprocs);
  for (int lev = 0; lev < levels; lev++)
  {
    MPI_Comm_split(MPI_COMM_WORLD, mylev, lev, &Commlev[lev]);
    checkPatchList(PatL[lev], false);
    Parallel::distribute(PatL[lev], end_rank[lev] - start_rank[lev] + 1, ingfs, fngfs, false, start_rank[lev], end_rank[lev]);
#if (RPB == 1)
#error "not support yet"
#endif
  }
  /* note different comm field has its own rank index
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    if(myrank==nprocs-1)
    {
        cout<<"myrank = "<<myrank<<", mylev = "<<mylev<<endl;
        MPI_Comm_rank(Commlev[levels-1],&myrank);
        cout<<myrank<<" :)"<<endl;
    }
  */
}

//================================================================================================

#if (PSTR == 1)
void cgh::construct_mylev(int nprocs)
{
  if (nprocs < levels)
  {
    cout << "Too few procs to use parallel level methods!" << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  start_rank = new int[levels];
  end_rank = new int[levels];

  int myrank;

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  int mp;
  mp = nprocs / levels;

  start_rank[0] = 0;
  end_rank[0] = mp - 1;
  for (int lev = 1; lev < levels - 1; lev++)
  {
    start_rank[lev] = end_rank[lev - 1] + 1;
    end_rank[lev] = end_rank[lev - 1] + mp;
  }
  start_rank[levels - 1] = end_rank[levels - 2] + 1;
  end_rank[levels - 1] = nprocs - 1;

  for (int lev = 0; lev < levels; lev++)
  {
    if (myrank >= start_rank[lev] && myrank <= end_rank[lev])
      mylev = lev;
  }
}
#elif (PSTR == 2)
void cgh::construct_mylev(int nprocs)
{
  if (nprocs < levels)
  {
    cout << "Too few procs to use parallel level methods!" << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  start_rank = new int[levels];
  end_rank = new int[levels];

  int myrank;

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  int mp;
  mp = nprocs / levels;

  start_rank[levels - 1] = 0;
  end_rank[levels - 1] = mp - 1;
  for (int lev = levels - 2; lev > 0; lev--)
  {
    start_rank[lev] = end_rank[lev - 1] + 1;
    end_rank[lev] = end_rank[lev - 1] + mp;
  }
  start_rank[0] = end_rank[1] + 1;
  end_rank[0] = nprocs - 1;

  for (int lev = levels - 1; lev >= 0; lev--)
  {
    if (myrank >= start_rank[lev] && myrank <= end_rank[lev])
      mylev = lev;
  }
}
#endif

#elif (PSTR == 3)
void cgh::construct_mylev(int nprocs)
{
  if (nprocs <= 1)
  {
    cout << " cgh::construct_mylev requires at least 2 procs" << endl;
    exit(0);
  }

  start_rank = new int[2];
  end_rank = new int[2];

  int myrank;

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  int mp;
  mp = nprocs / 2;

  // for other levels
  for (int lev = 0; lev < levels - 1; lev++)
  {
    start_rank[lev] = 0;
    end_rank[lev] = mp - 1;
  }
  // for finest level
  start_rank[levels - 1] = end_rank[0] + 1;
  end_rank[levels - 1] = nprocs - 1;

  if (myrank >= start_rank[0] && myrank <= end_rank[0])
    mylev = -1; // for other levels
  else
    mylev = 1; // for finest level
}


//-----------------------------------------------------------------------


void cgh::compose_cgh(int nprocs)
{
  Commlev = new MPI_Comm[levels];
  construct_mylev(nprocs);

  for (int lev = 0; lev < levels - 1; lev++)
  {
    MPI_Comm_split(MPI_COMM_WORLD, mylev, -1, &Commlev[lev]);
  }
  MPI_Comm_split(MPI_COMM_WORLD, mylev, 1, &Commlev[levels - 1]);

  for (int lev = 0; lev < levels; lev++)
  {
    checkPatchList(PatL[lev], false);
    Parallel::distribute(PatL[lev], end_rank[lev] - start_rank[lev] + 1, ingfs, fngfs, false, start_rank[lev], end_rank[lev]);
#if (RPB == 1)
#error "not support yet"
#endif
  }
}
#endif


void cgh::sethandle(monitor *ErrorMonitor)
{
  int BH_num;
  Porgls = new double **[levels];
  char filename[100];
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
    if (!inf.good() && ErrorMonitor && ErrorMonitor->outfile)
    {
      ErrorMonitor->outfile << "Can not open parameter file " << filename << " for inputing information of black holes" << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    for (int i = 1; inf.good(); i++)
    {
      inf.getline(pline, LEN);
      str = pline;

      int status = misc::parse_parts(str, sgrp, skey, sval, sind);
      if (status == -1)
      {
        if (ErrorMonitor && ErrorMonitor->outfile)
          ErrorMonitor->outfile << "error reading parameter file " << filename << " in line " << i << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
      else if (status == 0)
        continue;

      if (sgrp == "BSSN" && skey == "BH_num")
        BH_num = atoi(sval.c_str());
      else if (sgrp == "cgh" && skey == "moving levels start from")
      {
        movls = atoi(sval.c_str());
        movls = Mymin(movls, levels);
        movls = Mymax(0, movls);
      }
    }
    inf.close();
  }
  for (int lev = 0; lev < levels; lev++)
  {
    Porgls[lev] = new double *[BH_num];
    for (int i = 0; i < BH_num; i++)
      Porgls[lev][i] = new double[dim];
  }
  // read parameter from file
  {
    const int LEN = 256;
    char pline[LEN];
    string str, sgrp, skey, sval;
    int sind;
    ifstream inf(filename, ifstream::in);
    if (!inf.good() && ErrorMonitor && ErrorMonitor->outfile)
    {
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
        if (ErrorMonitor && ErrorMonitor->outfile)
          ErrorMonitor->outfile << "error reading parameter file " << filename << " in line " << i << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
      else if (status == 0)
        continue;

      if (sgrp == "BSSN" && sind < BH_num)
      {
        if (skey == "Porgx")
        {
          for (int lev = 0; lev < levels; lev++)
            Porgls[lev][sind][0] = atof(sval.c_str());
        }
        else if (skey == "Porgy")
        {
          for (int lev = 0; lev < levels; lev++)
            Porgls[lev][sind][1] = atof(sval.c_str());
        }
        else if (skey == "Porgz")
        {
          for (int lev = 0; lev < levels; lev++)
            Porgls[lev][sind][2] = atof(sval.c_str());
        }
      }
    }
    inf.close();
  }

  for (int lev = 0; lev < movls; lev++)
    for (int grd = 0; grd < grids[lev]; grd++)
      for (int i = 0; i < dim; i++)
        handle[lev][grd][i] = 0;

  if (movls < levels)
  {
    if (ErrorMonitor && ErrorMonitor->I_Print)
    {
      cout << endl;
      cout << " moving levels are lev #" << movls << "--" << levels - 1 << endl;
      cout << endl;
    }

    for (int lev = movls; lev < levels; lev++)
      for (int grd = 0; grd < grids[lev]; grd++)
      {
#if 0	
	 int bht=0;
	 for(int bhi=0;bhi<BH_num;bhi++)
	 {
	    bool flag=false;
    
	    for(int i=0;i<dim;i++)
               if(Porgls[0][bhi][i] < bbox[lev][grd][i] || Porgls[0][bhi][i] > bbox[lev][grd][i+dim]) {flag=true; break;}
	    if(flag) continue;
	    bht++;
	    if(bht==1)  for(int i=0;i<dim;i++) handle[lev][grd][i]=Porgls[0][bhi][i];
	    else if(ErrorMonitor && ErrorMonitor->outfile) 
	    {
               ErrorMonitor->outfile<<"cgh::sethandle: lev#"<<lev<<" grd#"<<grd<<" has too many black holes"<<endl;
               MPI_Abort(MPI_COMM_WORLD,1);
	    }
	 }
#else
        double xxc[dim], dis0, dis1;
        for (int i = 0; i < dim; i++)
          xxc[i] = (bbox[lev][grd][i] + bbox[lev][grd][i + dim]) / 2;
        int bht = 0;
        for (int bhi = 0; bhi < BH_num; bhi++)
        {
          if (bhi == 0)
          {
            dis0 = 0;
            for (int i = 0; i < dim; i++)
              dis0 += pow(Porgls[0][bhi][i] - xxc[i], 2);
            dis0 = sqrt(dis0);
          }
          else
          {
            dis1 = 0;
            for (int i = 0; i < dim; i++)
              dis1 += pow(Porgls[0][bhi][i] - xxc[i], 2);
            dis1 = sqrt(dis1);
            if (dis0 > dis1)
            {
              bht = bhi;
              dis0 = dis1;
            } // chose nearest one
          }
        }
        for (int i = 0; i < dim; i++)
          handle[lev][grd][i] = Porgls[0][bht][i];
#endif
      }
  }
  else if (ErrorMonitor && ErrorMonitor->I_Print)
  {
    if (levels > 1)
      cout << "fixed mesh refinement!" << endl;
    else
      cout << "unigrid simulation!" << endl;
  }

  BH_num_in = BH_num;
}
void cgh::checkPatchList(MyList<Patch> *PatL, bool buflog)
{
  while (PatL)
  {
    PatL->data->checkPatch(buflog);
    PatL = PatL->next;
  }
}


//================================================================================================

// 该成员函数用于将网格移动

//================================================================================================

void cgh::Regrid(int Symmetry, int BH_num, double **Porgbr, double **Porg0,
                 MyList<var> *OldList, MyList<var> *StateList,
                 MyList<var> *FutureList, MyList<var> *tmList, bool BB,
                 monitor *ErrorMonitor)
{
  // for moving part
  if (movls < levels)
  {
    bool tot_flag = false;
    bool *lev_flag;
    double **tmpPorg;
    tmpPorg = new double *[BH_num];
    for (int bhi = 0; bhi < BH_num; bhi++)
    {
      tmpPorg[bhi] = new double[dim];
      for (int i = 0; i < dim; i++)
        tmpPorg[bhi][i] = Porgbr[bhi][i];
    }
    lev_flag = new bool[levels - movls];
    for (int lev = movls; lev < levels; lev++)
    {
      lev_flag[lev - movls] = false;
      for (int grd = 0; grd < grids[lev]; grd++)
      {
        int flag;
        int do_every = 2;
        double dX = PatL[lev]->data->blb->data->getdX(0);
        double dY = PatL[lev]->data->blb->data->getdX(1);
        double dZ = PatL[lev]->data->blb->data->getdX(2);
        double rr;
        // make sure that the grid corresponds to the black hole
        int bhi = 0;
        for (bhi = 0; bhi < BH_num; bhi++)
        {
          // because finner level may also change Porgbr, so we need factor 2
          if (feq(Porgbr[bhi][0], handle[lev][grd][0], 2 * do_every * dX) &&
              feq(Porgbr[bhi][1], handle[lev][grd][1], 2 * do_every * dY) &&
              feq(Porgbr[bhi][2], handle[lev][grd][2], 2 * do_every * dZ))
            break;
        }
        if (bhi == BH_num)
        {
          // if the box has already touched the original point
          if (feq(0, bbox[lev][grd][0], dX / 2) &&
              feq(0, bbox[lev][grd][1], dY / 2) &&
              feq(0, bbox[lev][grd][2], dZ / 2))
            break;

          if (BH_num == 1)
          {
            bhi = 0;
            break;
          } // if only one black hole, it definitely match!

          if (ErrorMonitor->outfile)
          {
            ErrorMonitor->outfile << "cgh::Regrid: no black hole matches with grid lev#" << lev << " grd#" << grd
                                  << " with handle (" << handle[lev][grd][0] << "," << handle[lev][grd][1] << "," << handle[lev][grd][2] << ")" << endl;
            ErrorMonitor->outfile << "black holes' old positions:" << endl;
            for (bhi = 0; bhi < BH_num; bhi++)
              ErrorMonitor->outfile << "#" << bhi << ": (" << Porgbr[bhi][0] << "," << Porgbr[bhi][1] << "," << Porgbr[bhi][2] << ")" << endl;
            ErrorMonitor->outfile << "tolerance:" << endl;
            ErrorMonitor->outfile << "(" << 2 * do_every * dX << "," << 2 * do_every * dY << "," << 2 * do_every * dZ << ")" << endl;
            ErrorMonitor->outfile << "box lower boundary: (" << bbox[lev][grd][0] << "," << bbox[lev][grd][1] << "," << bbox[lev][grd][2] << ")" << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
          }

          delete[] lev_flag;
          for (bhi = 0; bhi < BH_num; bhi++)
            delete[] tmpPorg[bhi];
          delete[] tmpPorg;
          return;
        }
        // x direction
        rr = (Porg0[bhi][0] - handle[lev][grd][0]) / dX;
        if (rr > 0)
          flag = int(rr + 0.5) / do_every;
        else
          flag = int(rr - 0.5) / do_every;
        flag = flag * do_every;
        rr = bbox[lev][grd][0] + flag * dX;
        // pay attention to the symmetric case
        if (Symmetry == 2 && rr < 0)
          rr = -bbox[lev][grd][0];
        else
          rr = flag * dX;

        if (fabs(rr) > dX / 2)
        {
          lev_flag[lev - movls] = tot_flag = true;
          bbox[lev][grd][0] = bbox[lev][grd][0] + rr;
          bbox[lev][grd][3] = bbox[lev][grd][3] + rr;
          handle[lev][grd][0] += rr;
          tmpPorg[bhi][0] = Porg0[bhi][0];
        }

        // y direction
        rr = (Porg0[bhi][1] - handle[lev][grd][1]) / dY;
        if (rr > 0)
          flag = int(rr + 0.5) / do_every;
        else
          flag = int(rr - 0.5) / do_every;
        flag = flag * do_every;
        rr = bbox[lev][grd][1] + flag * dY;
        // pay attention to the symmetric case
        if (Symmetry == 2 && rr < 0)
          rr = -bbox[lev][grd][1];
        else
          rr = flag * dY;

        if (fabs(rr) > dY / 2)
        {
          lev_flag[lev - movls] = tot_flag = true;
          bbox[lev][grd][1] = bbox[lev][grd][1] + rr;
          bbox[lev][grd][4] = bbox[lev][grd][4] + rr;
          handle[lev][grd][1] += rr;
          tmpPorg[bhi][1] = Porg0[bhi][1];
        }

        // z direction
        rr = (Porg0[bhi][2] - handle[lev][grd][2]) / dZ;
        if (rr > 0)
          flag = int(rr + 0.5) / do_every;
        else
          flag = int(rr - 0.5) / do_every;
        flag = flag * do_every;
        rr = bbox[lev][grd][2] + flag * dZ;
        // pay attention to the symmetric case
        if (Symmetry > 0 && rr < 0)
          rr = -bbox[lev][grd][1];
        else
          rr = flag * dZ;

        if (fabs(rr) > dZ / 2)
        {
          lev_flag[lev - movls] = tot_flag = true;
          bbox[lev][grd][2] = bbox[lev][grd][2] + rr;
          bbox[lev][grd][5] = bbox[lev][grd][5] + rr;
          handle[lev][grd][2] += rr;
          tmpPorg[bhi][2] = Porg0[bhi][2];
        }
      }
      //   if(ErrorMonitor->outfile && lev_flag[lev-movls]) cout<<"lev#"<<lev<<"'s boxes moved"<<endl;
    }

    if (tot_flag)
    {
      int nprocs;
      MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
      recompose_cgh(nprocs, lev_flag, OldList, StateList, FutureList, tmList, Symmetry, BB);
      for (int bhi = 0; bhi < BH_num; bhi++)
      {
        for (int i = 0; i < dim; i++)
          Porgbr[bhi][i] = tmpPorg[bhi][i];
      }
    }

    delete[] lev_flag;
    for (int bhi = 0; bhi < BH_num; bhi++)
      delete[] tmpPorg[bhi];
    delete[] tmpPorg;
  }
}

//================================================================================================

#if (PSTR == 1 || PSTR == 2 || PSTR == 3)
#warning "Regrid is not implimented yet"
void cgh::Regrid_fake(int Symmetry, int BH_num, double **Porgbr, double **Porg0,
                      MyList<var> *OldList, MyList<var> *StateList,
                      MyList<var> *FutureList, MyList<var> *tmList, bool BB,
                      monitor *ErrorMonitor)
{
  // for moving part
  if (movls < levels)
  {
    bool tot_flag = false;
    bool *lev_flag;
    double **tmpPorg;
    tmpPorg = new double *[BH_num];
    for (int bhi = 0; bhi < BH_num; bhi++)
    {
      tmpPorg[bhi] = new double[dim];
      for (int i = 0; i < dim; i++)
        tmpPorg[bhi][i] = Porgbr[bhi][i];
    }
    lev_flag = new bool[levels - movls];
    for (int lev = movls; lev < levels; lev++)
    {
      lev_flag[lev - movls] = false;
      for (int grd = 0; grd < grids[lev]; grd++)
      {
        int flag;
        int do_every = 2;
        double dX = PatL[lev]->data->blb->data->getdX(0);
        double dY = PatL[lev]->data->blb->data->getdX(1);
        double dZ = PatL[lev]->data->blb->data->getdX(2);
        double rr;
        // make sure that the grid corresponds to the black hole
        int bhi = 0;
        for (bhi = 0; bhi < BH_num; bhi++)
        {
          // because finner level may also change Porgbr, so we need factor 2
          if (feq(Porgbr[bhi][0], handle[lev][grd][0], 2 * do_every * dX) &&
              feq(Porgbr[bhi][1], handle[lev][grd][1], 2 * do_every * dY) &&
              feq(Porgbr[bhi][2], handle[lev][grd][2], 2 * do_every * dZ))
            break;
        }
        if (bhi == BH_num)
        {
          // if the box has already touched the original point
          if (feq(0, bbox[lev][grd][0], dX / 2) &&
              feq(0, bbox[lev][grd][1], dY / 2) &&
              feq(0, bbox[lev][grd][2], dZ / 2))
            break;

          if (BH_num == 1)
          {
            bhi = 0;
            break;
          } // if only one black hole, it definitely match!

          if (ErrorMonitor->outfile)
          {
            ErrorMonitor->outfile << "cgh::Regrid: no black hole matches with grid lev#" << lev << " grd#" << grd
                                  << " with handle (" << handle[lev][grd][0] << "," << handle[lev][grd][1] << "," << handle[lev][grd][2] << ")" << endl;
            ErrorMonitor->outfile << "black holes' old positions:" << endl;
            for (bhi = 0; bhi < BH_num; bhi++)
              ErrorMonitor->outfile << "#" << bhi << ": (" << Porgbr[bhi][0] << "," << Porgbr[bhi][1] << "," << Porgbr[bhi][2] << ")" << endl;
            ErrorMonitor->outfile << "tolerance:" << endl;
            ErrorMonitor->outfile << "(" << 2 * do_every * dX << "," << 2 * do_every * dY << "," << 2 * do_every * dZ << ")" << endl;
            ErrorMonitor->outfile << "box lower boundary: (" << bbox[lev][grd][0] << "," << bbox[lev][grd][1] << "," << bbox[lev][grd][2] << ")" << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
          }

          delete[] lev_flag;
          for (bhi = 0; bhi < BH_num; bhi++)
            delete[] tmpPorg[bhi];
          delete[] tmpPorg;
          return;
        }
        // x direction
        rr = (Porg0[bhi][0] - handle[lev][grd][0]) / dX;
        if (rr > 0)
          flag = int(rr + 0.5) / do_every;
        else
          flag = int(rr - 0.5) / do_every;
        flag = flag * do_every;
        rr = bbox[lev][grd][0] + flag * dX;
        // pay attention to the symmetric case
        if (Symmetry == 2 && rr < 0)
          rr = -bbox[lev][grd][0];
        else
          rr = flag * dX;

        if (fabs(rr) > dX / 2)
        {
          lev_flag[lev - movls] = tot_flag = true;
          bbox[lev][grd][0] = bbox[lev][grd][0] + rr;
          bbox[lev][grd][3] = bbox[lev][grd][3] + rr;
          handle[lev][grd][0] += rr;
          tmpPorg[bhi][0] = Porg0[bhi][0];
        }

        // y direction
        rr = (Porg0[bhi][1] - handle[lev][grd][1]) / dY;
        if (rr > 0)
          flag = int(rr + 0.5) / do_every;
        else
          flag = int(rr - 0.5) / do_every;
        flag = flag * do_every;
        rr = bbox[lev][grd][1] + flag * dY;
        // pay attention to the symmetric case
        if (Symmetry == 2 && rr < 0)
          rr = -bbox[lev][grd][1];
        else
          rr = flag * dY;

        if (fabs(rr) > dY / 2)
        {
          lev_flag[lev - movls] = tot_flag = true;
          bbox[lev][grd][1] = bbox[lev][grd][1] + rr;
          bbox[lev][grd][4] = bbox[lev][grd][4] + rr;
          handle[lev][grd][1] += rr;
          tmpPorg[bhi][1] = Porg0[bhi][1];
        }

        // z direction
        rr = (Porg0[bhi][2] - handle[lev][grd][2]) / dZ;
        if (rr > 0)
          flag = int(rr + 0.5) / do_every;
        else
          flag = int(rr - 0.5) / do_every;
        flag = flag * do_every;
        rr = bbox[lev][grd][2] + flag * dZ;
        // pay attention to the symmetric case
        if (Symmetry > 0 && rr < 0)
          rr = -bbox[lev][grd][1];
        else
          rr = flag * dZ;

        if (fabs(rr) > dZ / 2)
        {
          lev_flag[lev - movls] = tot_flag = true;
          bbox[lev][grd][2] = bbox[lev][grd][2] + rr;
          bbox[lev][grd][5] = bbox[lev][grd][5] + rr;
          handle[lev][grd][2] += rr;
          tmpPorg[bhi][2] = Porg0[bhi][2];
        }
      }
      //   if(ErrorMonitor->outfile && lev_flag[lev-movls]) cout<<"lev#"<<lev<<"'s boxes moved"<<endl;
    }

    if (tot_flag)
    {
      int nprocs;
      MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
      recompose_cgh_fake(nprocs, lev_flag, OldList, StateList, FutureList, tmList, Symmetry, BB);
      for (int bhi = 0; bhi < BH_num; bhi++)
      {
        for (int i = 0; i < dim; i++)
          Porgbr[bhi][i] = tmpPorg[bhi][i];
      }
    }

    delete[] lev_flag;
    for (int bhi = 0; bhi < BH_num; bhi++)
      delete[] tmpPorg[bhi];
    delete[] tmpPorg;
  }
}
#endif


//================================================================================================

// 该成员函数用于重构网格

//================================================================================================

#if (PSTR == 0)
void cgh::recompose_cgh(int nprocs, bool *lev_flag,
                        MyList<var> *OldList, MyList<var> *StateList,
                        MyList<var> *FutureList, MyList<var> *tmList,
                        int Symmetry, bool BB)
{
  for (int lev = movls; lev < levels; lev++)
    if (lev_flag[lev - movls])
    {
      MyList<Patch> *tmPat = 0;
      tmPat = construct_patchlist(lev, Symmetry);
      // tmPat construction completes
      Parallel::distribute(tmPat, nprocs, ingfs, fngfs, false);
      //    checkPatchList(tmPat,true);
      bool CC = (lev > trfls);
      Parallel::fill_level_data(tmPat, PatL[lev], PatL[lev - 1], OldList, StateList, FutureList, tmList, Symmetry, BB, CC);

      Parallel::KillBlocks(PatL[lev]);
      PatL[lev]->destroyList();
      PatL[lev] = tmPat;
#if (RPB == 1)
      Parallel::destroypsuList_bam(bdsul[lev]);
      Parallel::destroypsuList_bam(rsul[lev]);
      Parallel::Constr_pointstr_OutBdLow2Hi(PatL[lev], PatL[lev - 1], bdsul[lev]);
      Parallel::Constr_pointstr_Restrict(PatL[lev], PatL[lev - 1], rsul[lev]);
#endif
    }
}
#elif (PSTR == 1 || PSTR == 2 || PSTR == 3)
#warning "recompose_cgh is not implimented yet"
void cgh::recompose_cgh(int nprocs, bool *lev_flag,
                        MyList<var> *OldList, MyList<var> *StateList,
                        MyList<var> *FutureList, MyList<var> *tmList,
                        int Symmetry, bool BB)
{
  for (int lev = movls; lev < levels; lev++)
    if (lev_flag[lev - movls])
    {
      MyList<Patch> *tmPat = 0;
      tmPat = construct_patchlist(lev, Symmetry);
      // tmPat construction completes
      Parallel::distribute(tmPat, end_rank[lev] - start_rank[lev] + 1, ingfs, fngfs, false, start_rank[lev], end_rank[lev]);
      //    checkPatchList(tmPat,true);
      bool CC = (lev > trfls);
      Parallel::fill_level_data(tmPat, PatL[lev], PatL[lev - 1], OldList, StateList, FutureList, tmList, Symmetry, BB, CC);

      Parallel::KillBlocks(PatL[lev]);
      PatL[lev]->destroyList();
      PatL[lev] = tmPat;
#if (RPB == 1)
#error "not support yet"
#endif
    }
}

//================================================================================================

void cgh::recompose_cgh_fake(int nprocs, bool *lev_flag,
                             MyList<var> *OldList, MyList<var> *StateList,
                             MyList<var> *FutureList, MyList<var> *tmList,
                             int Symmetry, bool BB)
{
  for (int lev = movls; lev < levels; lev++)
    if (lev_flag[lev - movls] && lev != mylev)
    {
      MyList<Patch> *tmPat = 0;
      tmPat = construct_patchlist(lev, Symmetry);
      // tmPat construction completes
      Parallel::distribute(tmPat, end_rank[lev] - start_rank[lev] + 1, ingfs, fngfs, false, start_rank[lev], end_rank[lev]);

      Parallel::KillBlocks(PatL[lev]);
      PatL[lev]->destroyList();
      PatL[lev] = tmPat;
    }
}
#endif

//================================================================================================

// 该成员函数从输入文件中读取网格信息

//================================================================================================

void cgh::read_bbox(int Symmetry, char *filename)
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  // read parameter from file
  {
    const int LEN = 256;
    char pline[LEN];
    string str, sgrp, skey, sval;
    int sind1, sind2, sind3;
    ifstream inf(filename, ifstream::in);
    if (!inf.good() && myrank == 0)
    {
      cout << "cgh::cgh: Can not open parameter file " << filename << " for inputing information of black holes" << endl;
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
  shape = new int **[levels];
  handle = new double **[levels];
  bbox = new double **[levels];
  PatL = new MyList<Patch> *[levels];
  Lt = new double[levels];
#if (RPB == 1)
  bdsul = new MyList<Parallel::pointstru_bam> *[levels];
  rsul = new MyList<Parallel::pointstru_bam> *[levels];
#endif
  // read parameter from file
  {
    const int LEN = 256;
    char pline[LEN];
    string str, sgrp, skey, sval;
    int sind1, sind2, sind3;
    ifstream inf(filename, ifstream::in);
    if (!inf.good() && myrank == 0)
    {
      cout << "cgh::cgh: Can not open parameter file " << filename << " for inputing information of black holes" << endl;
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
    }
    inf.close();
  }

  for (int sind1 = 0; sind1 < levels; sind1++)
  {
    shape[sind1] = new int *[grids[sind1]];
    handle[sind1] = new double *[grids[sind1]];
    bbox[sind1] = new double *[grids[sind1]];
    for (int sind2 = 0; sind2 < grids[sind1]; sind2++)
    {
      shape[sind1][sind2] = new int[dim];
      handle[sind1][sind2] = new double[dim];
      bbox[sind1][sind2] = new double[2 * dim];
    }
  }
  // read parameter from file
  {
    const int LEN = 256;
    char pline[LEN];
    string str, sgrp, skey, sval;
    int sind1, sind2, sind3;
    ifstream inf(filename, ifstream::in);
    if (!inf.good() && myrank == 0)
    {
      cout << "cgh::cgh: Can not open parameter file " << filename << " for inputing information of black holes" << endl;
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

      if (sgrp == "cgh" && sind1 < levels && sind2 < grids[sind1])
      {
        if (skey == "bbox")
          bbox[sind1][sind2][sind3] = atof(sval.c_str());
        else if (skey == "shape")
          shape[sind1][sind2][sind3] = atoi(sval.c_str());
      }
    }
    inf.close();
  }
// we always assume the input parameter is in cell center style
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
  for (int lev = 0; lev < levels; lev++)
    for (int grd = 0; grd < grids[lev]; grd++)
    {
      for (int i = 0; i < dim; i++)
      {

        shape[lev][grd][i] = shape[lev][grd][i] + 1;
      }
    }
#endif

  {

    // boxes align check
    double DH0[dim];
    for (int i = 0; i < dim; i++)
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
      DH0[i] = (bbox[0][0][i + dim] - bbox[0][0][i]) / (shape[0][0][i] - 1);
#else
#ifdef Cell
      DH0[i] = (bbox[0][0][i + dim] - bbox[0][0][i]) / shape[0][0][i];
#else
#error Not define Vertex nor Cell
#endif
#endif
    for (int lev = 0; lev < levels; lev++)
      for (int grd = 0; grd < grids[lev]; grd++)
        Parallel::aligncheck(bbox[0][0], bbox[lev][grd], lev, DH0, shape[lev][grd]);

#if 0 // we do not need it here, because we do it in construct_patchlist
// extend buffer points for shell overlap
#ifdef WithShell   
   for(int i=0;i<dim;i++)
   {
     shape[0][0][i] += buffer_width;
     bbox[0][0][i+dim] += DH0[i]*buffer_width;
   }
   if(Symmetry < 2 )
   {
     shape[0][0][0] += buffer_width;
     bbox[0][0][0] -= DH0[0]*buffer_width;
     shape[0][0][1] += buffer_width;
     bbox[0][0][1] -= DH0[1]*buffer_width;
     if(Symmetry < 1)
     {
       shape[0][0][2] += buffer_width;
       bbox[0][0][2] -= DH0[2]*buffer_width;
     }
   }
#endif
#endif
  }
  // print information of cgh
  if (myrank == 0)
  {
    cout << endl;
    cout << " cgh has levels: " << levels << endl;
    cout << endl;
    for (int lev = 0; lev < levels; lev++)
    {
      cout << " level #" << lev << " has boxes: " << grids[lev] << endl;
      for (int grd = 0; grd < grids[lev]; grd++)
      {
        cout << " #" << grd << " box is" << "  (" << bbox[lev][grd][0] << ":" << bbox[lev][grd][3]
             << "," << bbox[lev][grd][1] << ":" << bbox[lev][grd][4]
             << "," << bbox[lev][grd][2] << ":" << bbox[lev][grd][5]
             << ")." << endl;
      }
    }
  }
}

//================================================================================================


//================================================================================================

// 该成员函数用于生成所需要的网格信息

//================================================================================================

MyList<Patch> *cgh::construct_patchlist(int lev, int Symmetry)
{
  // Construct Patches
  MyList<Patch> *tmPat = 0;
  // construct box list
  MyList<Parallel::gridseg> *boxes = 0, *gs;

  /*
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (myrank == 0)
  {
      cout << " construct patchlist: " << " level = " << lev << ",  grids in this level = " << grids[lev] << endl;
  }
  */

  for (int grd = 0; grd < grids[lev]; grd++)
  {
    if (boxes)
    {
      gs->next = new MyList<Parallel::gridseg>;
      gs = gs->next;
      gs->data = new Parallel::gridseg;
    }
    else
    {
      boxes = gs = new MyList<Parallel::gridseg>;
      gs->data = new Parallel::gridseg;
    }
    for (int i = 0; i < dim; i++)
    {
      gs->data->llb[i] = bbox[lev][grd][i];
      gs->data->uub[i] = bbox[lev][grd][dim + i];
      gs->data->shape[i] = shape[lev][grd][i];
    }
    gs->data->Bg = 0;
    gs->next = 0;
  }

  // 对网格盒子进行合并（三个以上盒子会出bug）
  // Parallel::merge_gsl(boxes, ratio);
  if (grids[lev] < 3)
  {
    Parallel::merge_gsl(boxes, ratio);
  }

  // 当网格盒子相交时，对网格盒子进行重新切割
  // Parallel::cut_gsl(boxes);
  if (grids[lev] < 3)
  {
    Parallel::cut_gsl(boxes);
  }

  // 切割后增加新的ghost部分？
  // Parallel::add_ghost_touch(boxes);
  if (grids[lev] < 3)
  {
    Parallel::add_ghost_touch(boxes);
  }

  MyList<Patch> *gp;
  gs = boxes;
  while (gs)
  {
    double tbb[2 * dim];
    if (tmPat)
    {
      gp->next = new MyList<Patch>;
      gp = gp->next;
      for (int i = 0; i < dim; i++)
      {
        tbb[i] = gs->data->llb[i];
        tbb[dim + i] = gs->data->uub[i];
      }
#ifdef WithShell
      gp->data = new Patch(3, gs->data->shape, tbb, lev, true, Symmetry);
#else
      gp->data = new Patch(3, gs->data->shape, tbb, lev, (lev > 0), Symmetry);
#endif
    }
    else
    {
      tmPat = gp = new MyList<Patch>;
      for (int i = 0; i < dim; i++)
      {
        tbb[i] = gs->data->llb[i];
        tbb[dim + i] = gs->data->uub[i];
      }
#ifdef WithShell
      gp->data = new Patch(3, gs->data->shape, tbb, lev, true, Symmetry);
#else
      gp->data = new Patch(3, gs->data->shape, tbb, lev, (lev > 0), Symmetry);
#endif
    }
    gp->next = 0;

    gs = gs->next;
  }

  boxes->destroyList();

  return tmPat;
}

//================================================================================================


bool cgh::Interp_One_Point(MyList<var> *VarList,
                           double *XX, /*input global Cartesian coordinate*/
                           double *Shellf, int Symmetry)
{
  int lev = levels - 1;
  while (lev >= 0)
  {
    MyList<Patch> *Pp = PatL[lev];
    while (Pp)
    {
#if (PSTR == 0)
      if (Pp->data->Interp_ONE_Point(VarList, XX, Shellf, Symmetry))
        return true;
#elif (PSTR == 1 || PSTR == 2 || PSTR == 3)
      if (Pp->data->Interp_ONE_Point(VarList, XX, Shellf, Symmetry, Commlev[lev]))
        return true;
#endif
      Pp = Pp->next;
    }
    lev--;
  }
  return false;
}


void cgh::Regrid_Onelevel(int lev, int Symmetry, int BH_num, double **Porgbr, double **Porg0,
                          MyList<var> *OldList, MyList<var> *StateList,
                          MyList<var> *FutureList, MyList<var> *tmList, bool BB,
                          monitor *ErrorMonitor)
{
  if (lev < movls)
    return;

#if (0)
  // #if (PSTR == 1 || PSTR == 2)
  MyList<Patch> *Pp = PatL[lev];
  while (Pp)
  {
    Pp->data->checkPatch(0, start_rank[mylev]);
    Pp = Pp->next;
  }
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (myrank == start_rank[mylev])
  {
    cout << "out_rank = " << myrank << endl;
    for (int grd = 0; grd < grids[lev]; grd++)
    {
      cout << handle[lev][grd][0] << "," << handle[lev][grd][1] << "," << handle[lev][grd][2] << endl;
    }
    for (int bhi = 0; bhi < BH_num; bhi++)
    {
      cout << Porgls[lev][bhi][0] << "," << Porgls[lev][bhi][1] << "," << Porgls[lev][bhi][2] << endl;
      cout << Porg0[bhi][0] << "," << Porg0[bhi][1] << "," << Porg0[bhi][2] << endl;
    }
  }
#endif

  //   misc::tillherecheck(Commlev[lev],start_rank[lev],"start Regrid_Onelevel");
  // for moving part
  bool tot_flag = false;
  double **tmpPorg;
  tmpPorg = new double *[BH_num];
  for (int bhi = 0; bhi < BH_num; bhi++)
  {
    tmpPorg[bhi] = new double[dim];
    for (int i = 0; i < dim; i++)
      tmpPorg[bhi][i] = Porgls[lev][bhi][i];
  }

  for (int grd = 0; grd < grids[lev]; grd++)
  {
    int flag;
    int do_every = 2;
    double dX = PatL[lev]->data->blb->data->getdX(0);
    double dY = PatL[lev]->data->blb->data->getdX(1);
    double dZ = PatL[lev]->data->blb->data->getdX(2);
    double rr;
    // make sure that the grid corresponds to the black hole
    int bhi = 0;
    for (bhi = 0; bhi < BH_num; bhi++)
    {
      // because finner level may also change Porgbr, so we need factor 2
      // now I used Porgls
      if (feq(Porgls[lev][bhi][0], handle[lev][grd][0], 2 * do_every * dX) &&
          feq(Porgls[lev][bhi][1], handle[lev][grd][1], 2 * do_every * dY) &&
          feq(Porgls[lev][bhi][2], handle[lev][grd][2], 2 * do_every * dZ))
        break;
    }
    if (bhi == BH_num)
    {
      // if the box has already touched the original point
      if (feq(0, bbox[lev][grd][0], dX / 2) &&
          feq(0, bbox[lev][grd][1], dY / 2) &&
          feq(0, bbox[lev][grd][2], dZ / 2))
        break;

      if (BH_num == 1)
      {
        bhi = 0;
        break;
      } // if only one black hole, it definitely match!

      if (ErrorMonitor->outfile)
      {
        ErrorMonitor->outfile << "cgh::Regrid: no black hole matches with grid lev#" << lev << " grd#" << grd
                              << " with handle (" << handle[lev][grd][0] << "," << handle[lev][grd][1] << "," << handle[lev][grd][2] << ")" << endl;
        ErrorMonitor->outfile << "black holes' old positions:" << endl;
        for (bhi = 0; bhi < BH_num; bhi++)
          ErrorMonitor->outfile << "#" << bhi << ": (" << Porgls[lev][bhi][0] << "," << Porgls[lev][bhi][1] << ","
                                << Porgls[lev][bhi][2] << ")" << endl;
        ErrorMonitor->outfile << "tolerance:" << endl;
        ErrorMonitor->outfile << "(" << 2 * do_every * dX << "," << 2 * do_every * dY << "," << 2 * do_every * dZ << ")" << endl;
        ErrorMonitor->outfile << "box lower boundary: (" << bbox[lev][grd][0] << "," << bbox[lev][grd][1] << "," << bbox[lev][grd][2] << ")" << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      for (bhi = 0; bhi < BH_num; bhi++)
        delete[] tmpPorg[bhi];
      delete[] tmpPorg;
      return;
    }
    // x direction
    rr = (Porg0[bhi][0] - handle[lev][grd][0]) / dX;
    if (rr > 0)
      flag = int(rr + 0.5) / do_every;
    else
      flag = int(rr - 0.5) / do_every;
    flag = flag * do_every;
    rr = bbox[lev][grd][0] + flag * dX;
    // pay attention to the symmetric case
    if (Symmetry == 2 && rr < 0)
      rr = -bbox[lev][grd][0];
    else
      rr = flag * dX;

    if (fabs(rr) > dX / 2)
    {
      tot_flag = true;
      bbox[lev][grd][0] = bbox[lev][grd][0] + rr;
      bbox[lev][grd][3] = bbox[lev][grd][3] + rr;
      handle[lev][grd][0] += rr;
      tmpPorg[bhi][0] = Porg0[bhi][0];
    }

    // y direction
    rr = (Porg0[bhi][1] - handle[lev][grd][1]) / dY;
    if (rr > 0)
      flag = int(rr + 0.5) / do_every;
    else
      flag = int(rr - 0.5) / do_every;
    flag = flag * do_every;
    rr = bbox[lev][grd][1] + flag * dY;
    // pay attention to the symmetric case
    if (Symmetry == 2 && rr < 0)
      rr = -bbox[lev][grd][1];
    else
      rr = flag * dY;

    if (fabs(rr) > dY / 2)
    {
      tot_flag = true;
      bbox[lev][grd][1] = bbox[lev][grd][1] + rr;
      bbox[lev][grd][4] = bbox[lev][grd][4] + rr;
      handle[lev][grd][1] += rr;
      tmpPorg[bhi][1] = Porg0[bhi][1];
    }

    // z direction
    rr = (Porg0[bhi][2] - handle[lev][grd][2]) / dZ;
    if (rr > 0)
      flag = int(rr + 0.5) / do_every;
    else
      flag = int(rr - 0.5) / do_every;
    flag = flag * do_every;
    rr = bbox[lev][grd][2] + flag * dZ;
    // pay attention to the symmetric case
    if (Symmetry > 0 && rr < 0)
      rr = -bbox[lev][grd][1];
    else
      rr = flag * dZ;

    if (fabs(rr) > dZ / 2)
    {
      tot_flag = true;
      bbox[lev][grd][2] = bbox[lev][grd][2] + rr;
      bbox[lev][grd][5] = bbox[lev][grd][5] + rr;
      handle[lev][grd][2] += rr;
      tmpPorg[bhi][2] = Porg0[bhi][2];
    }
  }

  //   misc::tillherecheck(Commlev[lev],start_rank[lev],"after tot_flag check");

  if (tot_flag)
  {
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    //     misc::tillherecheck(Commlev[lev],start_rank[lev],"before recompose_cgh_Onelevel");

    recompose_cgh_Onelevel(nprocs, lev, OldList, StateList, FutureList, tmList, Symmetry, BB);

    //     misc::tillherecheck(Commlev[lev],start_rank[lev],"after recompose_cgh_Onelevel");

    for (int bhi = 0; bhi < BH_num; bhi++)
    {
      for (int i = 0; i < dim; i++)
        Porgls[lev][bhi][i] = tmpPorg[bhi][i];
    }

#if (PSTR == 1 || PSTR == 2)
//       MyList<Patch> *Pp=PatL[lev];
//       while(Pp)
//       {
//	 Pp->data->checkPatch(0,start_rank[mylev]);
//	 Pp=Pp->next;
//       }
#endif
  }

  for (int bhi = 0; bhi < BH_num; bhi++)
    delete[] tmpPorg[bhi];
  delete[] tmpPorg;
}


#if (PSTR == 0)
void cgh::recompose_cgh_Onelevel(int nprocs, int lev,
                                 MyList<var> *OldList, MyList<var> *StateList,
                                 MyList<var> *FutureList, MyList<var> *tmList,
                                 int Symmetry, bool BB)
{
  MyList<Patch> *tmPat = 0;
  tmPat = construct_patchlist(lev, Symmetry);
  // tmPat construction completes
  Parallel::distribute(tmPat, nprocs, ingfs, fngfs, false);
  //    checkPatchList(tmPat,true);
  bool CC = (lev > trfls);
  Parallel::fill_level_data(tmPat, PatL[lev], PatL[lev - 1], OldList, StateList, FutureList, tmList, Symmetry, BB, CC);

  Parallel::KillBlocks(PatL[lev]);
  PatL[lev]->destroyList();
  PatL[lev] = tmPat;
}
#elif (PSTR == 1 || PSTR == 2 || PSTR == 3)
#warning "recompose_cgh_Onelevel is not implimented yet"
void cgh::recompose_cgh_Onelevel(int nprocs, int lev,
                                 MyList<var> *OldList, MyList<var> *StateList,
                                 MyList<var> *FutureList, MyList<var> *tmList,
                                 int Symmetry, bool BB)
{
  MyList<Patch> *tmPat = 0;
  misc::tillherecheck(Commlev[lev], start_rank[lev], "before construct_patchlist");
  tmPat = construct_patchlist(lev, Symmetry);
  misc::tillherecheck(Commlev[lev], start_rank[lev], "after construct_patchlist");
  // tmPat construction completes
  Parallel::distribute(tmPat, end_rank[lev] - start_rank[lev] + 1, ingfs, fngfs, false, start_rank[lev], end_rank[lev]);
  misc::tillherecheck(Commlev[lev], start_rank[lev], "after distribute");
  //    checkPatchList(tmPat,true);
  bool CC = (lev > trfls);
  Parallel::fill_level_data(tmPat, PatL[lev], PatL[lev - 1], OldList, StateList, FutureList, tmList, Symmetry, BB, CC);
  misc::tillherecheck(Commlev[lev], start_rank[lev], "after fill_level_data");

  Parallel::KillBlocks(PatL[lev]);
  PatL[lev]->destroyList();
  PatL[lev] = tmPat;
}


// the input lev is lower level for regrid
void cgh::Regrid_Onelevel_aux(int lev, int Symmetry, int BH_num, double **Porgbr, double **Porg0,
                              MyList<var> *OldList, MyList<var> *StateList,
                              MyList<var> *FutureList, MyList<var> *tmList, bool BB,
                              monitor *ErrorMonitor)
{
  lev++;
  if (lev < movls)
    return;

  // for moving part
  bool tot_flag = false;
  double **tmpPorg;
  tmpPorg = new double *[BH_num];
  for (int bhi = 0; bhi < BH_num; bhi++)
  {
    tmpPorg[bhi] = new double[dim];
    for (int i = 0; i < dim; i++)
      tmpPorg[bhi][i] = Porgbr[bhi][i];
  }

  for (int grd = 0; grd < grids[lev]; grd++)
  {
    int flag;
    int do_every = 2;
    double dX = PatL[lev]->data->blb->data->getdX(0);
    double dY = PatL[lev]->data->blb->data->getdX(1);
    double dZ = PatL[lev]->data->blb->data->getdX(2);
    double rr;
    // make sure that the grid corresponds to the black hole
    int bhi = 0;
    for (bhi = 0; bhi < BH_num; bhi++)
    {
      // because finner level may also change Porgbr, so we need factor 2
      if (feq(Porgbr[bhi][0], handle[lev][grd][0], 2 * do_every * dX) &&
          feq(Porgbr[bhi][1], handle[lev][grd][1], 2 * do_every * dY) &&
          feq(Porgbr[bhi][2], handle[lev][grd][2], 2 * do_every * dZ))
        break;
    }
    if (bhi == BH_num)
    {
      // if the box has already touched the original point
      if (feq(0, bbox[lev][grd][0], dX / 2) &&
          feq(0, bbox[lev][grd][1], dY / 2) &&
          feq(0, bbox[lev][grd][2], dZ / 2))
        break;

      if (BH_num == 1)
      {
        bhi = 0;
        break;
      } // if only one black hole, it definitely match!

      if (ErrorMonitor->outfile)
      {
        ErrorMonitor->outfile << "cgh::Regrid: no black hole matches with grid lev#" << lev << " grd#" << grd
                              << " with handle (" << handle[lev][grd][0] << "," << handle[lev][grd][1] << "," << handle[lev][grd][2] << ")" << endl;
        ErrorMonitor->outfile << "black holes' old positions:" << endl;
        for (bhi = 0; bhi < BH_num; bhi++)
          ErrorMonitor->outfile << "#" << bhi << ": (" << Porgbr[bhi][0] << "," << Porgbr[bhi][1] << "," << Porgbr[bhi][2] << ")" << endl;
        ErrorMonitor->outfile << "tolerance:" << endl;
        ErrorMonitor->outfile << "(" << 2 * do_every * dX << "," << 2 * do_every * dY << "," << 2 * do_every * dZ << ")" << endl;
        ErrorMonitor->outfile << "box lower boundary: (" << bbox[lev][grd][0] << "," << bbox[lev][grd][1] << "," << bbox[lev][grd][2] << ")" << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      for (bhi = 0; bhi < BH_num; bhi++)
        delete[] tmpPorg[bhi];
      delete[] tmpPorg;
      return;
    }
    // x direction
    rr = (Porg0[bhi][0] - handle[lev][grd][0]) / dX;
    if (rr > 0)
      flag = int(rr + 0.5) / do_every;
    else
      flag = int(rr - 0.5) / do_every;
    flag = flag * do_every;
    rr = bbox[lev][grd][0] + flag * dX;
    // pay attention to the symmetric case
    if (Symmetry == 2 && rr < 0)
      rr = -bbox[lev][grd][0];
    else
      rr = flag * dX;

    if (fabs(rr) > dX / 2)
    {
      tot_flag = true;
      bbox[lev][grd][0] = bbox[lev][grd][0] + rr;
      bbox[lev][grd][3] = bbox[lev][grd][3] + rr;
      handle[lev][grd][0] += rr;
      tmpPorg[bhi][0] = Porg0[bhi][0];
    }

    // y direction
    rr = (Porg0[bhi][1] - handle[lev][grd][1]) / dY;
    if (rr > 0)
      flag = int(rr + 0.5) / do_every;
    else
      flag = int(rr - 0.5) / do_every;
    flag = flag * do_every;
    rr = bbox[lev][grd][1] + flag * dY;
    // pay attention to the symmetric case
    if (Symmetry == 2 && rr < 0)
      rr = -bbox[lev][grd][1];
    else
      rr = flag * dY;

    if (fabs(rr) > dY / 2)
    {
      tot_flag = true;
      bbox[lev][grd][1] = bbox[lev][grd][1] + rr;
      bbox[lev][grd][4] = bbox[lev][grd][4] + rr;
      handle[lev][grd][1] += rr;
      tmpPorg[bhi][1] = Porg0[bhi][1];
    }

    // z direction
    rr = (Porg0[bhi][2] - handle[lev][grd][2]) / dZ;
    if (rr > 0)
      flag = int(rr + 0.5) / do_every;
    else
      flag = int(rr - 0.5) / do_every;
    flag = flag * do_every;
    rr = bbox[lev][grd][2] + flag * dZ;
    // pay attention to the symmetric case
    if (Symmetry > 0 && rr < 0)
      rr = -bbox[lev][grd][1];
    else
      rr = flag * dZ;

    if (fabs(rr) > dZ / 2)
    {
      tot_flag = true;
      bbox[lev][grd][2] = bbox[lev][grd][2] + rr;
      bbox[lev][grd][5] = bbox[lev][grd][5] + rr;
      handle[lev][grd][2] += rr;
      tmpPorg[bhi][2] = Porg0[bhi][2];
    }
  }

  if (tot_flag)
  {
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    recompose_cgh_Onelevel(nprocs, lev, OldList, StateList, FutureList, tmList, Symmetry, BB);
  }

  for (int bhi = 0; bhi < BH_num; bhi++)
    delete[] tmpPorg[bhi];
  delete[] tmpPorg;
}
#endif


void cgh::settrfls(const int lev)
{
  trfls = lev;
}
