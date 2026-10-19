
#ifndef PARALLEL_BAM_H
#define PARALLEL_BAM_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <new>
using namespace std;

#include "var.h"
#include "MPatch.h"
#include "Block.h"
#include "MyList.h"
#include "macrodef.h"
namespace Parallel
{
  struct pointstru_bam
  {
    double pox[dim]; // cordinate
    Block *Bgs;      // interplate from
    Block *Bgd;      // interplate for
    double *coef;    // interpolation coefficients
    int sind[dim];   // interpolation starting array index
  };
  void destroypsuList_bam(MyList<pointstru_bam> *ct);
  void OutBdLow2Hi_bam(MyList<Patch> *PLc, MyList<Patch> *PLf,
                       MyList<var> *VarList1 /* source */, MyList<var> *VarList2 /* target */,
                       int Symmetry);
  void OutBdLow2Hi_bam(MyList<Patch> *PLc, MyList<Patch> *PLf,
                       MyList<var> *VarList1 /* source */, MyList<var> *VarList2 /* target */,
                       MyList<Parallel::pointstru_bam> *bdsul, int Symmetry);
  void Constr_pointstr_OutBdLow2Hi(MyList<Patch> *PLf, MyList<Patch> *PLc,
                                   MyList<Parallel::pointstru_bam> *&bdsul);
  void Restrict_bam(MyList<Patch> *PLc, MyList<Patch> *PLf,
                    MyList<var> *VarList1 /* source */, MyList<var> *VarList2 /* target */,
                    int Symmetry);
  void Restrict_bam(MyList<Patch> *PLc, MyList<Patch> *PLf,
                    MyList<var> *VarList1 /* source */, MyList<var> *VarList2 /* target */,
                    MyList<Parallel::pointstru_bam> *rsul, int Symmetry);
  void Constr_pointstr_Restrict(MyList<Patch> *PLf, MyList<Patch> *PLc,
                                MyList<Parallel::pointstru_bam> *&rsul);
  void intertransfer(MyList<Parallel::pointstru_bam> *&sul,
                     MyList<var> *VarList1 /* source */, MyList<var> *VarList2 /*target */,
                     int Symmetry);
  int interdata_packer(double *data, MyList<Parallel::pointstru_bam> *sul, int myrank, int node, int dir,
                       MyList<var> *VarLists /* source */, MyList<var> *VarListd /* target */, int Symmetry);
}
#endif /*PARALLEL_BAM_H */
