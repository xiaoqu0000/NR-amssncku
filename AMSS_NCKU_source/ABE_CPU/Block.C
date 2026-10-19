
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <new>
using namespace std;

#include "Block.h"
#include "misc.h"

Block::Block(int DIM, int *shapei, double *bboxi, int ranki, int ingfsi, int fngfsi, int levi, const int cgpui) : rank(ranki), ingfs(ingfsi), fngfs(fngfsi), lev(levi), cgpu(cgpui)
{
  for (int i = 0; i < dim; i++)
    X[i] = 0;

  if (DIM != dim)
  {
    cout << "dimension is not consistent in Block construction" << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  bool flag = false;
  for (int i = 0; i < dim; i++)
  {
    shape[i] = shapei[i];
    if (shape[i] <= 0)
      flag = true;
    bbox[i] = bboxi[i];
    bbox[dim + i] = bboxi[dim + i];
  }

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (flag)
  {
    cout << "myrank: " << myrank << ", on rank: " << rank << endl;
    cout << "error shape in Block construction: (" << shape[0] << "," << shape[1] << "," << shape[2] << ")" << endl;
    cout << "box boundary: (" << bbox[0] << ":" << bbox[3] << "," << bbox[1] << ":" << bbox[4] << "," << bbox[2] << ":" << bbox[5] << ")" << endl;
    cout << "belong to level " << lev << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

#ifndef FAKECHECK
  if (myrank == rank)
  {
    for (int i = 0; i < dim; i++)
    {
      X[i] = new double[shape[i]];
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
      double h = (bbox[dim + i] - bbox[i]) / (shape[i] - 1);
      for (int j = 0; j < shape[i]; j++)
        X[i][j] = bbox[i] + j * h;
#else
#ifdef Cell
      double h = (bbox[dim + i] - bbox[i]) / shape[i];
      for (int j = 0; j < shape[i]; j++)
        X[i][j] = bbox[i] + (j + 0.5) * h;
#else
#error Not define Vertex nor Cell
#endif
#endif
    }

    int nn = shape[0] * shape[1] * shape[2];
    fgfs = new double *[fngfs];
    for (int i = 0; i < fngfs; i++)
    {
      fgfs[i] = (double *)malloc(sizeof(double) * nn);
      if (!(fgfs[i]))
      {
        cout << "on node#" << rank << ", out of memory when constructing Block." << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
      memset(fgfs[i], 0, sizeof(double) * nn);
    }

    igfs = new int *[ingfs];
    for (int i = 0; i < ingfs; i++)
    {
      igfs[i] = (int *)malloc(sizeof(int) * nn);
      if (!(igfs[i]))
      {
        cout << "on node#" << rank << ", out of memory when constructing Block." << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
      memset(igfs[i], 0, sizeof(int) * nn);
    }
  }
#endif
}
Block::~Block()
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (myrank == rank)
  {
    for (int i = 0; i < dim; i++)
      delete[] X[i];
    for (int i = 0; i < ingfs; i++)
      free(igfs[i]);
    delete[] igfs;
    for (int i = 0; i < fngfs; i++)
      free(fgfs[i]);
    delete[] fgfs;
    X[0] = X[1] = X[2] = 0;
    igfs = 0;
    fgfs = 0;
  }
}
void Block::checkBlock()
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (myrank == 0)
  {
    cout << "belong to level " << lev << endl;
    cout << "shape: [";
    for (int i = 0; i < dim; i++)
    {
      cout << shape[i];
      if (i < dim - 1)
        cout << ",";
      else
        cout << "]";
    }
    cout << " resolution: [";
    for (int i = 0; i < dim; i++)
    {
      cout << getdX(i);
      if (i < dim - 1)
        cout << ",";
      else
        cout << "]" << endl;
    }
    cout << "locate on node " << rank << ", at (includes ghost zone):" << endl;
    cout << "(";
    for (int i = 0; i < dim; i++)
    {
      cout << bbox[i] << ":" << bbox[dim + i];
      if (i < dim - 1)
        cout << ",";
      else
        cout << ")" << endl;
    }
    cout << "has " << ingfs << " int type grids functions," << fngfs << " double type grids functions" << endl;
  }
}
double Block::getdX(int dir)
{
  if (dir < 0 || dir >= dim)
  {
    cout << "Block::getdX: error input dir = " << dir << ", this Block has direction (0," << dim - 1 << ")" << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  double h;
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
  if (shape[dir] == 1)
  {
    cout << "Block::getdX: for direction " << dir << ", this Block has only one point. Can not determine dX for vertex center grid." << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  h = (bbox[dim + dir] - bbox[dir]) / (shape[dir] - 1);
#else
#ifdef Cell
  h = (bbox[dim + dir] - bbox[dir]) / shape[dir];
#else
#error Not define Vertex nor Cell
#endif
#endif
  return h;
}
void Block::swapList(MyList<var> *VarList1, MyList<var> *VarList2, int myrank)
{
  if (rank == myrank)
  {
    MyList<var> *varl1 = VarList1, *varl2 = VarList2;
    while (varl1 && varl2)
    {
      misc::swap<double *>(fgfs[varl1->data->sgfn], fgfs[varl2->data->sgfn]);
      varl1 = varl1->next;
      varl2 = varl2->next;
    }
    if (varl1 || varl2)
    {
      cout << "error in Block::swaplist, var lists does not match." << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }
}
