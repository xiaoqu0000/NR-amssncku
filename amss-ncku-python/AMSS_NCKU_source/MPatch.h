
#ifndef PATCH_H
#define PATCH_H

#include <mpi.h>
#include "MyList.h"
#include "Block.h"
#include "var.h"
#include "macrodef.h" //need dim here; Vertex or Cell; ghost_width

class Patch
{

public:
   int lev;
   int shape[dim];
   double bbox[2 * dim]; // this bbox includes buffer points
   MyList<Block> *blb, *ble;
   int lli[dim], uui[dim]; // denote the buffer points on each boundary

public:
   Patch() {};
   Patch(int DIM, int *shapei, double *bboxi, int levi, bool buflog, int Symmetry);

   ~Patch();

   void checkPatch(bool buflog);
   void checkPatch(bool buflog, const int out_rank);
   void checkBlock();
   void Interp_Points(MyList<var> *VarList,
                      int NN, double **XX,
                      double *Shellf, int Symmetry);
   bool Interp_ONE_Point(MyList<var> *VarList, double *XX,
                         double *Shellf, int Symmetry);
   double getdX(int dir);

   void Find_Maximum(MyList<var> *VarList, double *XX,
                     double *Shellf);

   bool Find_Point(double *XX);

   void Interp_Points(MyList<var> *VarList,
                      int NN, double **XX,
                      double *Shellf, int Symmetry, MPI_Comm Comm_here);
   bool Interp_ONE_Point(MyList<var> *VarList, double *XX,
                         double *Shellf, int Symmetry, MPI_Comm Comm_here);
   void Find_Maximum(MyList<var> *VarList, double *XX,
                     double *Shellf, MPI_Comm Comm_here);
};

#endif /* PATCH_H */
