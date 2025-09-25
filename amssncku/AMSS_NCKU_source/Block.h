
#ifndef BLOCK_H
#define BLOCK_H

#include <mpi.h>
#include "macrodef.h" //need dim here; Vertex or Cell
#include "var.h"
#include "MyList.h"
class Block
{

public:
   int shape[dim];
   double bbox[2 * dim];
   double *X[dim];
   int rank; // where the real data locate in
   int lev, cgpu;
   int ingfs, fngfs;
   int *(*igfs);
   double *(*fgfs);

public:
   Block() {};
   Block(int DIM, int *shapei, double *bboxi, int ranki, int ingfsi, int fngfs, int levi, const int cgpui = 0);

   ~Block();

   void checkBlock();

   double getdX(int dir);
   void swapList(MyList<var> *VarList1, MyList<var> *VarList2, int myrank);
};

#endif /* BLOCK_H */
