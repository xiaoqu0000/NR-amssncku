
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <new>
using namespace std;

#include "NullShellPatch2.h"
#include "Parallel.h"
#include "fmisc.h"
#include "misc.h"
#include "shellfunctions.h"
#include "NullEvol.h"
#include "NullNews.h"
#include "initial_null2.h"
#include "rungekutta4_rout.h"
#include "kodiss.h"

#define PI M_PI

NullShellPatch2::NullShellPatch2(int *shapei, double Rmini, double xmini, double xmaxi, int Symmetryi, int myranki) : myrank(myranki), Rmin(Rmini), xmin(xmini), xmax(xmaxi), PatL(0), Symmetry(Symmetryi)
{
  for (int i = 0; i < dim; i++)
  {
    shape[i] = shapei[i];
// we always assume the input parameter is in cell center style
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
    shape[i] = shape[i] + 1;
#endif
  }

  if (myrank == 0)
  {
    cout << "null shell's range: r = [" << xmin * Rmin / (1 - xmin) << ":";
    if (xmax == 1)
      cout << "+Infty]" << endl;
    else
      cout << xmax * Rmin / (1 - xmax) << "]" << endl;
    cout << "                    x = [" << xmin << ":" << xmax << "]" << endl
         << "shape: " << shape[2] << endl
         << "resolution: [" << getdX(0) << "," << getdX(1) << "," << getdX(2) << "]" << endl;
  }
// in order to touch infinity, we always use vertex center in r direction
// for Cell center it is some fake as following
#ifdef Cell
#ifdef Vertex
#error Both Cell and Vertex are defined
#endif
  {
    double ht = (xmax - xmin) / shape[2];
    xmax = xmax + ht / 2;
    xmin = xmin - ht / 2;
    shape[2] = shape[2] + 1;
  }
#endif

  double bbox[2 * dim];
  int shape_here[dim];
  bbox[2] = xmin;
  bbox[5] = xmax;
  shape_here[2] = shape[2];

  switch (Symmetry)
  {
  case 0:
    for (int i = 0; i < 2; i++)
      shape_here[i] = shape[i] + 2 * overghost;
    bbox[0] = -PI / 4 - overghost * getdX(0);
    bbox[1] = -PI / 4 - overghost * getdX(1);
    bbox[3] = PI / 4 + overghost * getdX(0);
    bbox[4] = PI / 4 + overghost * getdX(1);
    PatL = new MyList<ss_patch>;
    PatL->data = new xp_npatch(ingfs, fngfs, shape_here, bbox, myrank);
    PatL->insert(new xm_npatch(ingfs, fngfs, shape_here, bbox, myrank));
    PatL->insert(new yp_npatch(ingfs, fngfs, shape_here, bbox, myrank));
    PatL->insert(new ym_npatch(ingfs, fngfs, shape_here, bbox, myrank));
    PatL->insert(new zp_npatch(ingfs, fngfs, shape_here, bbox, myrank));
    PatL->insert(new zm_npatch(ingfs, fngfs, shape_here, bbox, myrank));
    break;
  case 1:
    for (int i = 0; i < 2; i++)
      shape_here[i] = shape[i] + 2 * overghost;
    bbox[0] = -PI / 4 - overghost * getdX(0);
    bbox[1] = -PI / 4 - overghost * getdX(1);
    bbox[3] = PI / 4 + overghost * getdX(0);
    bbox[4] = PI / 4 + overghost * getdX(1);
    PatL = new MyList<ss_patch>;
    PatL->data = new zp_npatch(ingfs, fngfs, shape_here, bbox, myrank);
    shape_here[0] = shape[0] + 2 * overghost;
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
    shape_here[1] = (shape[1] + 1) / 2 + overghost;
#else
#ifdef Cell
    shape_here[1] = shape[1] / 2 + overghost;
#else
#error Not define Vertex nor Cell
#endif
#endif
    bbox[0] = -PI / 4 - overghost * getdX(0);
    shape_here[1] += ghost_width;
    bbox[1] = -ghost_width * getdX(1); // buffer points method to deal with boundary
    bbox[3] = PI / 4 + overghost * getdX(0);
    bbox[4] = PI / 4 + overghost * getdX(1);
    PatL->insert(new xp_npatch(ingfs, fngfs, shape_here, bbox, myrank));
    PatL->insert(new yp_npatch(ingfs, fngfs, shape_here, bbox, myrank));
    bbox[0] = -PI / 4 - overghost * getdX(0);
    bbox[1] = -PI / 4 - overghost * getdX(1);
    bbox[3] = PI / 4 + overghost * getdX(0);
    bbox[4] = ghost_width * getdX(1); // buffer points method to deal with boundary
    PatL->insert(new xm_npatch(ingfs, fngfs, shape_here, bbox, myrank));
    PatL->insert(new ym_npatch(ingfs, fngfs, shape_here, bbox, myrank));
    break;
  case 2:
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
    for (int i = 0; i < 2; i++)
      shape_here[i] = (shape[i] + 1) / 2 + overghost;
#else
#ifdef Cell
    for (int i = 0; i < 2; i++)
      shape_here[i] = shape[i] / 2 + overghost;
#else
#error Not define Vertex nor Cell
#endif
#endif
    shape_here[0] += ghost_width;
    shape_here[1] += ghost_width;
    bbox[0] = -ghost_width * getdX(0); // buffer points method to deal with boundary
    bbox[1] = -ghost_width * getdX(1); // buffer points method to deal with boundary
    bbox[3] = PI / 4 + overghost * getdX(0);
    bbox[4] = PI / 4 + overghost * getdX(1);
    PatL = new MyList<ss_patch>;
    PatL->data = new zp_npatch(ingfs, fngfs, shape_here, bbox, myrank);
    PatL->insert(new xp_npatch(ingfs, fngfs, shape_here, bbox, myrank));
    PatL->insert(new yp_npatch(ingfs, fngfs, shape_here, bbox, myrank));
    break;
  default:
    cout << "not recognized Symmetry type" << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  int ngfs = 0;
  gx = new var("gx", ngfs++, 1, 1, 1);
  gy = new var("gy", ngfs++, 1, 1, 1);
  gz = new var("gz", ngfs++, 1, 1, 1);

  g00 = new var("g00", ngfs++, 1, 1, 1);
  g01 = new var("g01", ngfs++, -1, 1, 1);
  p02 = new var("p02", ngfs++, 1, -1, 1);
  p03 = new var("p03", ngfs++, 1, 1, -1);
  g02 = new var("g02", ngfs++, 1, -1, 1);
  g03 = new var("g03", ngfs++, 1, 1, -1);
  Theta22 = new var("Theta22", ngfs++, 1, 1, 1);
  Theta23 = new var("Theta23", ngfs++, 1, -1, -1);
  Theta33 = new var("Theta33", ngfs++, 1, 1, 1);

  g22o = new var("g22o", ngfs++, 1, 1, 1);
  g23o = new var("g23o", ngfs++, 1, -1, -1);
  g33o = new var("g33o", ngfs++, 1, 1, 1);
  g220 = new var("g220", ngfs++, 1, 1, 1);
  g230 = new var("g230", ngfs++, 1, -1, -1);
  g330 = new var("g330", ngfs++, 1, 1, 1);
  g22 = new var("g22", ngfs++, 1, 1, 1);
  g23 = new var("g23", ngfs++, 1, -1, -1);
  g33 = new var("g33", ngfs++, 1, 1, 1);
  g221 = new var("g221", ngfs++, 1, 1, 1);
  g231 = new var("g231", ngfs++, 1, -1, -1);
  g331 = new var("g331", ngfs++, 1, 1, 1);
  g22_rhs = new var("g22_rhs", ngfs++, 1, 1, 1);
  g23_rhs = new var("g23_rhs", ngfs++, 1, -1, -1);
  g33_rhs = new var("g33_rhs", ngfs++, 1, 1, 1);

  RNews = new var("RNews", ngfs++, 1, 1, 1);
  INews = new var("INews", ngfs++, 1, 1, 1);
  omega = new var("omega", ngfs++, 1, 1, 1);
  dtomega = new var("dtomega", ngfs++, 1, 1, 1);

  DumpList = new MyList<var>(g220);
  DumpList->insert(g230);
  DumpList->insert(g330);

  OldStateList = new MyList<var>(g22o);
  OldStateList->insert(g23o);
  OldStateList->insert(g33o);
  StateList = new MyList<var>(g220);
  StateList->insert(g230);
  StateList->insert(g330);
  SynchList_pre = new MyList<var>(g22);
  SynchList_pre->insert(g23);
  SynchList_pre->insert(g33);
  RHSList = new MyList<var>(g22_rhs);
  RHSList->insert(g23_rhs);
  RHSList->insert(g33_rhs);
  SynchList_cor = new MyList<var>(g221);
  SynchList_cor->insert(g231);
  SynchList_cor->insert(g331);

  NewsList = new MyList<var>(RNews);
  NewsList->insert(INews);

  g01List = new MyList<var>(g01);
  g01wt = new double *[1];
  for (int ii = 0; ii < 1; ii++)
  {
    g01wt[ii] = new double[3];
    g01wt[ii][0] = g01wt[ii][1] = g01wt[ii][2] = 1;
  }

  pg0AList = new MyList<var>(p02);
  pg0AList->insert(p03);
  pg0AList->insert(g02);
  pg0AList->insert(g03);
  pg0Awt = new double *[4];
  for (int ii = 0; ii < 4; ii++)
  {
    pg0Awt[ii] = new double[3];
    pg0Awt[ii][0] = pg0Awt[ii][1] = pg0Awt[ii][2] = 1;
  }
  pg0Awt[0][0] = pg0Awt[1][1] = pg0Awt[2][0] = pg0Awt[3][1] = -1;

  g00List = new MyList<var>(g00);
  g00wt = new double *[1];
  for (int ii = 0; ii < 1; ii++)
  {
    g00wt[ii] = new double[3];
    g00wt[ii][0] = g00wt[ii][1] = g00wt[ii][2] = 1;
  }

  ThetaList = new MyList<var>(Theta22);
  ThetaList->insert(Theta23);
  ThetaList->insert(Theta33);
  Thetawt = new double *[3];
  for (int ii = 0; ii < 3; ii++)
  {
    Thetawt[ii] = new double[3];
    Thetawt[ii][0] = Thetawt[ii][1] = Thetawt[ii][2] = 1;
  }
  Thetawt[1][0] = Thetawt[1][1] = -1;

  ingfs = 0;
  fngfs = ngfs;
}
NullShellPatch2::~NullShellPatch2()
{
  int nprocs = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  for (int node = 0; node < nprocs; node++)
  {
    if (ss_src[node])
      destroypsuList(ss_src[node]);
    if (ss_dst[node])
      destroypsuList(ss_dst[node]);
    if (cs_src)
    {
      if (cs_src[node])
        destroypsuList(cs_src[node]);
      if (cs_dst[node])
        destroypsuList(cs_dst[node]);
    }
  }

  delete[] ss_src;
  delete[] ss_dst;
  if (cs_src)
  {
    delete[] cs_src;
    delete[] cs_dst;
  }

  while (PatL)
  {
    ss_patch *sPp = PatL->data;
    MyList<Block> *bg;
    while (sPp->blb)
    {
      if (sPp->blb == sPp->ble)
        break;
      bg = (sPp->blb->next) ? sPp->blb->next : 0;
      delete sPp->blb->data;
      delete sPp->blb;
      sPp->blb = bg;
    }
    if (sPp->ble)
    {
      delete sPp->ble->data;
      delete sPp->ble;
    }
    sPp->blb = sPp->ble = 0;
    PatL = PatL->next;
  }
  PatL->destroyList();

  StateList->clearList();
  SynchList_pre->clearList();
  SynchList_cor->clearList();
  RHSList->clearList();
  OldStateList->clearList();
  DumpList->clearList();
  CheckList->clearList();

  NewsList->clearList();

  g01List->clearList();
  g00List->clearList();
  pg0AList->clearList();
  ThetaList->clearList();

  delete gx;
  delete gy;
  delete gz;

  delete g00;
  delete g01;
  delete p02;
  delete p03;
  delete g02;
  delete g03;
  delete Theta22;
  delete Theta23;
  delete Theta33;

  delete g22o;
  delete g23o;
  delete g33o;
  delete g220;
  delete g230;
  delete g330;
  delete g22;
  delete g23;
  delete g33;
  delete g221;
  delete g231;
  delete g331;
  delete g22_rhs;
  delete g23_rhs;
  delete g33_rhs;

  delete RNews;
  delete INews;
  delete omega;
  delete dtomega;

  for (int ii = 0; ii < 1; ii++)
    delete[] g01wt[ii];
  delete[] g01wt;
  for (int ii = 0; ii < 4; ii++)
    delete[] pg0Awt[ii];
  delete[] pg0Awt;
  for (int ii = 0; ii < 1; ii++)
    delete[] g00wt[ii];
  delete[] g00wt;
  for (int ii = 0; ii < 3; ii++)
    delete[] Thetawt[ii];
  delete[] Thetawt;
}
double NullShellPatch2::getdX(int dir)
{
  if (dir < 0 || dir >= dim)
  {
    cout << "NullShellPatch::getdX: error input dir = " << dir << ", this Patch has direction (0," << dim - 1 << ")" << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  double h;
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
  if (shape[dir] == 1)
  {
    cout << "NullShellPatch::getdX: for direction " << dir << ", this Patch has only one point. Can not determine dX for vertex center grid." << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  if (dir < 2)
    h = PI / 2 / (shape[dir] - 1);
  else
    h = (xmax - xmin) / (shape[dir] - 1);
#else
#ifdef Cell
  if (dir < 2)
    h = PI / 2 / shape[dir];
  else
    h = (xmax - xmin) / shape[dir];
#else
#error Not define Vertex nor Cell
#endif
#endif
  return h;
}
void NullShellPatch2::destroypsuList(MyList<pointstru> *ct)
{
  MyList<pointstru> *n;
  while (ct)
  {
    n = ct->next;
    if (ct->data->coef)
    {
      delete[] ct->data->coef;
      delete[] ct->data->sind;
    }
    delete ct->data;
    delete ct;
    ct = n;
  }
}
void NullShellPatch2::shellname(char *sn, int i)
{
  switch (i)
  {
  case 0:
    sprintf(sn, "zp");
    return;
  case 1:
    sprintf(sn, "zm");
    return;
  case 2:
    sprintf(sn, "xp");
    return;
  case 3:
    sprintf(sn, "xm");
    return;
  case 4:
    sprintf(sn, "yp");
    return;
  case 5:
    sprintf(sn, "ym");
    return;
  }
}
MyList<Block> *NullShellPatch2::compose_sh(int cpusize)
{
  if (dim != 3)
  {
    cout << "distrivute: now we only support 3-dimension" << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  //  checkPatch();

  MyList<Block> *BlL = 0;

  int split_size, min_size, block_size = 0;

  int min_width = 2 * Mymax(ghost_width, buffer_width);
  int nxy[2], mmin_width[2], min_shape[2];

  MyList<ss_patch> *PLi = PatL;
  for (int i = 0; i < 2; i++)
    min_shape[i] = PLi->data->shape[i];
  PLi = PLi->next;
  while (PLi)
  {
    ss_patch *PP = PLi->data;
    for (int i = 0; i < 2; i++)
      min_shape[i] = Mymin(min_shape[i], PP->shape[i]);
    PLi = PLi->next;
  }

  for (int i = 0; i < 2; i++)
    mmin_width[i] = Mymin(min_width, min_shape[i]);

  min_size = mmin_width[0];
  for (int i = 1; i < 2; i++)
    min_size = min_size * mmin_width[i];

  PLi = PatL;
  while (PLi)
  {
    ss_patch *PP = PLi->data;
    //    PP->checkPatch(true);
    int bs = PP->shape[0];
    for (int i = 1; i < 2; i++)
      bs = bs * PP->shape[i];
    block_size = block_size + bs;
    PLi = PLi->next;
  }
  split_size = Mymax(min_size, block_size / cpusize);
  split_size = Mymax(1, split_size);

  int n_rank = 0;
  PLi = PatL;
  int reacpu = 0;
  while (PLi)
  {
    ss_patch *PP = PLi->data;

    reacpu += Parallel::partition2(nxy, split_size, mmin_width, cpusize, PP->shape); // r direction can not be splitted!! It's ode!

    Block *ng;
    int shape_here[3], ibbox_here[2 * 2];
    double bbox_here[2 * 3], dd;

    // ibbox : 0,...N-1
    for (int i = 0; i < nxy[0]; i++)
      for (int j = 0; j < nxy[1]; j++)
      {
        ibbox_here[0] = (PP->shape[0] * i) / nxy[0];
        ibbox_here[2] = (PP->shape[0] * (i + 1)) / nxy[0] - 1;
        ibbox_here[1] = (PP->shape[1] * j) / nxy[1];
        ibbox_here[3] = (PP->shape[1] * (j + 1)) / nxy[1] - 1;

        ibbox_here[0] = Mymax(0, ibbox_here[0] - ghost_width);
        ibbox_here[2] = Mymin(PP->shape[0] - 1, ibbox_here[2] + ghost_width);
        ibbox_here[1] = Mymax(0, ibbox_here[1] - ghost_width);
        ibbox_here[3] = Mymin(PP->shape[1] - 1, ibbox_here[3] + ghost_width);

        shape_here[0] = ibbox_here[2] - ibbox_here[0] + 1;
        shape_here[1] = ibbox_here[3] - ibbox_here[1] + 1;
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
        dd = (PP->bbox[3] - PP->bbox[0]) / (PP->shape[0] - 1);
        bbox_here[0] = PP->bbox[0] + ibbox_here[0] * dd;
        bbox_here[3] = PP->bbox[0] + ibbox_here[2] * dd;

        dd = (PP->bbox[4] - PP->bbox[1]) / (PP->shape[1] - 1);
        bbox_here[1] = PP->bbox[1] + ibbox_here[1] * dd;
        bbox_here[4] = PP->bbox[1] + ibbox_here[3] * dd;
#else
#ifdef Cell
        dd = (PP->bbox[3] - PP->bbox[0]) / PP->shape[0];
        bbox_here[0] = PP->bbox[0] + (ibbox_here[0]) * dd;
        bbox_here[3] = PP->bbox[0] + (ibbox_here[2] + 1) * dd;

        dd = (PP->bbox[4] - PP->bbox[1]) / PP->shape[1];
        bbox_here[1] = PP->bbox[1] + (ibbox_here[1]) * dd;
        bbox_here[4] = PP->bbox[1] + (ibbox_here[3] + 1) * dd;
#else
#error Not define Vertex nor Cell
#endif
#endif
        shape_here[2] = PP->shape[2];
        bbox_here[2] = PP->bbox[2];
        bbox_here[5] = PP->bbox[5];
        ng = new Block(dim, shape_here, bbox_here, n_rank++, ingfs, fngfs, 0); // delete through KillBlocks
        //	    ng->checkBlock();
        if (n_rank == cpusize)
          n_rank = 0;
        if (BlL)
          BlL->insert(ng);
        else
          BlL = new MyList<Block>(ng); // delete through KillBlocks

        // set PP->blb
        if (i == 0 && j == 0)
        {
          MyList<Block> *Bp = BlL;
          while (Bp->data != ng)
            Bp = Bp->next;
          PP->blb = Bp;
        }
      }
    // set PP->ble
    {
      MyList<Block> *Bp = BlL;
      while (Bp->data != ng)
        Bp = Bp->next;
      PP->ble = Bp;
    }
    PLi = PLi->next;
  }
  if (reacpu < cpusize * 2 / 3)
  {
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank == 0)
      cout << "NullShellPatch::distribute CAUSTION: uses essencially " << reacpu << " processors vs " << cpusize << " cpus run, your scientific computation scale is not as large as you estimate." << endl;
  }

  return BlL;
}
void NullShellPatch2::Dump_Data(MyList<var> *DumpListi, char *tag, double time, double dT)
{
  MyList<ss_patch> *PP = PatL;
  while (PP)
  {
    //   round at 4 and 5
    int ncount = int(time / dT + 0.5);

    MPI_Status sta;
    int DIM = 3;
    double llb[3], uub[3];
    double DX, DY, DZ;

    double *databuffer = 0;
    if (myrank == 0)
    {
      databuffer = (double *)malloc(sizeof(double) * PP->data->shape[0] * PP->data->shape[1] * PP->data->shape[2]);
      if (!databuffer)
      {
        cout << "NullShellPatch::Dump_Data: out of memory when dumping data." << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }

    MyList<var> *DumpList = DumpListi;
    while (DumpList)
    {
      var *VP = DumpList->data;

      MyList<Block> *Bp = PP->data->blb;
      while (Bp)
      {
        Block *BP = Bp->data;
        if (BP->rank == 0 && myrank == 0)
        {
          DX = BP->getdX(0);
          DY = BP->getdX(1);
          DZ = BP->getdX(2);
          llb[0] = (feq(BP->bbox[0], PP->data->bbox[0], DX / 2)) ? BP->bbox[0] : BP->bbox[0] + ghost_width * DX;
          llb[1] = (feq(BP->bbox[1], PP->data->bbox[1], DY / 2)) ? BP->bbox[1] : BP->bbox[1] + ghost_width * DY;
          llb[2] = (feq(BP->bbox[2], PP->data->bbox[2], DZ / 2)) ? BP->bbox[2] : BP->bbox[2] + ghost_width * DZ;
          uub[0] = (feq(BP->bbox[3], PP->data->bbox[3], DX / 2)) ? BP->bbox[3] : BP->bbox[3] - ghost_width * DX;
          uub[1] = (feq(BP->bbox[4], PP->data->bbox[4], DY / 2)) ? BP->bbox[4] : BP->bbox[4] - ghost_width * DY;
          uub[2] = (feq(BP->bbox[5], PP->data->bbox[5], DZ / 2)) ? BP->bbox[5] : BP->bbox[5] - ghost_width * DZ;
          f_copy(DIM, PP->data->bbox, PP->data->bbox + DIM, PP->data->shape, databuffer, BP->bbox, BP->bbox + DIM, BP->shape, BP->fgfs[VP->sgfn], llb, uub);
        }
        else
        {
          int nnn = (BP->shape[0]) * (BP->shape[1]) * (BP->shape[2]);
          if (myrank == 0)
          {
            double *bufferhere = (double *)malloc(sizeof(double) * nnn);
            if (!bufferhere)
            {
              cout << "on node#" << myrank << ", out of memory when dumping data." << endl;
              MPI_Abort(MPI_COMM_WORLD, 1);
            }
            MPI_Recv(bufferhere, nnn, MPI_DOUBLE, BP->rank, 0, MPI_COMM_WORLD, &sta);
            DX = BP->getdX(0);
            DY = BP->getdX(1);
            DZ = BP->getdX(2);
            llb[0] = (feq(BP->bbox[0], PP->data->bbox[0], DX / 2)) ? BP->bbox[0] : BP->bbox[0] + ghost_width * DX;
            llb[1] = (feq(BP->bbox[1], PP->data->bbox[1], DY / 2)) ? BP->bbox[1] : BP->bbox[1] + ghost_width * DY;
            llb[2] = (feq(BP->bbox[2], PP->data->bbox[2], DZ / 2)) ? BP->bbox[2] : BP->bbox[2] + ghost_width * DZ;
            uub[0] = (feq(BP->bbox[3], PP->data->bbox[3], DX / 2)) ? BP->bbox[3] : BP->bbox[3] - ghost_width * DX;
            uub[1] = (feq(BP->bbox[4], PP->data->bbox[4], DY / 2)) ? BP->bbox[4] : BP->bbox[4] - ghost_width * DY;
            uub[2] = (feq(BP->bbox[5], PP->data->bbox[5], DZ / 2)) ? BP->bbox[5] : BP->bbox[5] - ghost_width * DZ;
            f_copy(DIM, PP->data->bbox, PP->data->bbox + DIM, PP->data->shape, databuffer, BP->bbox, BP->bbox + DIM, BP->shape, bufferhere, llb, uub);
            free(bufferhere);
          }
          else if (myrank == BP->rank)
          {
            MPI_Send(BP->fgfs[VP->sgfn], nnn, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
          }
        }
        if (Bp == PP->data->ble)
          break;
        Bp = Bp->next;
      }
      if (myrank == 0)
      {
        char filename[100];
        char sn[3];
        shellname(sn, PP->data->sst);
        if (tag)
          sprintf(filename, "%s_LevSH-%s_%s_%05d.bin", tag, sn, VP->name, ncount);
        else
          sprintf(filename, "LevSH-%s_%s_%05d.bin", sn, VP->name, ncount);

        Parallel::writefile(time, PP->data->shape[0], PP->data->shape[1], PP->data->shape[2],
                            PP->data->bbox[0], PP->data->bbox[3], PP->data->bbox[1], PP->data->bbox[4],
                            PP->data->bbox[2], PP->data->bbox[5], filename, databuffer);
      }
      DumpList = DumpList->next;
    }

    if (myrank == 0)
      free(databuffer);

    PP = PP->next;
  }
}
// Now we dump the data including overlap points
void NullShellPatch2::Dump_xyz(char *tag, double time, double dT)
{
  MyList<var> *DumpListi = 0;
  DumpListi = new MyList<var>(gx);
  DumpListi->insert(gy);
  DumpListi->insert(gz);
  Dump_Data(DumpListi, tag, time, dT);
  DumpListi->clearList();
}
// setup interpatch interpolation stuffs
void NullShellPatch2::setupintintstuff(int cpusize, MyList<Patch> *CPatL, int Symmetry)
{
  const int hCS_width = 0; // do not input data from null shell to box
  const int hSC_width = 1; // do     input data from box to null shell
  int myrank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (myrank == 0)
    cout << "NullShellPatch2::setupintintstuff begines..." << endl;

  ss_src = new MyList<pointstru> *[cpusize];
  ss_dst = new MyList<pointstru> *[cpusize];

  if (!CPatL) // if characteristic evolve alone
  {
    cs_src = 0;
    cs_dst = 0;
  }
  else
  {
    cs_src = new MyList<pointstru> *[cpusize];
    cs_dst = new MyList<pointstru> *[cpusize];
  }

  MyList<pointstru> *ps, *ts;
  MyList<ss_patch> *sPp;
  MyList<Block> *Bgl;
  MyList<Patch> *Pp;
  Block *Bg;
  double CDH[dim], DH[dim], llb[dim], uub[dim];
  double x, y, z;

  for (int i = 0; i < dim; i++)
  {
    if (CPatL)
      CDH[i] = CPatL->data->getdX(i);
    DH[i] = getdX(i);
  }

  for (int i = 0; i < cpusize; i++)
  {
    ss_src[i] = 0;
    ss_dst[i] = 0;
    if (CPatL)
    {
      cs_src[i] = 0;
      cs_dst[i] = 0;
    }
  }

  sPp = PatL;
  while (sPp)
  {
    for (int iz = 0; iz < sPp->data->shape[2]; iz++)
      for (int is = 0; is < sPp->data->shape[1]; is++)
        for (int ir = 0; ir < sPp->data->shape[0]; ir++)
        {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
          x = sPp->data->bbox[0] + ir * DH[0];
          y = sPp->data->bbox[1] + is * DH[1];
          z = sPp->data->bbox[2] + iz * DH[2];
#else
#ifdef Cell
          x = sPp->data->bbox[0] + (ir + 0.5) * DH[0];
          y = sPp->data->bbox[1] + (is + 0.5) * DH[1];
          z = sPp->data->bbox[2] + (iz + 0.5) * DH[2];
#else
#error Not define Vertex nor Cell
#endif
#endif
          if (CPatL && z < sPp->data->bbox[2] + (hSC_width + 0.0001) * DH[2])
          {
            double gx, gy, gz;
            getglobalpox(gx, gy, gz, sPp->data->sst, x, y, z);
            bool flag = false;
            for (int i = 0; i < cpusize; i++)
            {
              flag = prolongpointstru(cs_src[i], false, sPp->data->sst, PatL, DH, CPatL, CDH, gx, gy, gz, Symmetry, i, iz);
              if (flag)
                break;
            }
            if (!flag)
            {
              CPatL->data->checkBlock();
              if (myrank == 0)
              {
                cout << "ShellPatch2::prolongpointstru fail to find cardisian source point for" << endl;
                cout << "sst = " << sPp->data->sst << " lx,ly,lz = " << x << "," << y << "," << z << endl;
                cout << "x,y,z = " << gx << "," << gy << "," << gz << endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
              }
            }
          }
          if (x < -PI / 4 - (overghost - ghost_width - 0.0001) * DH[0] || x > PI / 4 + (overghost - ghost_width - 0.0001) * DH[0] ||
              y < -PI / 4 - (overghost - ghost_width - 0.0001) * DH[1] || y > PI / 4 + (overghost - ghost_width - 0.0001) * DH[1])
          {
            double gx, gy, gz;
            if (z < 1 - 0.0001 * DH[2])
              getglobalpox(gx, gy, gz, sPp->data->sst, x, y, z);
            bool flag = true;
            if (flag)
            {
              flag = false;
              for (int i = 0; i < cpusize; i++)
              {
                if (z < 1 - 0.0001 * DH[2])
                  flag = prolongpointstru(ss_src[i], true, sPp->data->sst, PatL, DH, CPatL, CDH, gx, gy, gz, Symmetry, i, iz);
                else
                  flag = prolongpointstru_ss(ss_src[i], sPp->data->sst, PatL, DH, CPatL, CDH, x, y, z, Symmetry, i, iz);
                if (flag)
                  break;
              }
              if (!flag)
              {
                if (myrank == 0)
                {
                  // if you used Vertex grid please note x=1, try 0.999999 instead
                  cout << "NullShellPatch2::prolongpointstru fail to find shell source point for" << endl;
                  cout << "sst = " << sPp->data->sst << " lx,ly,lz = " << x << "," << y << "," << z << endl;
                  MPI_Abort(MPI_COMM_WORLD, 1);
                }
              }
            }
          }
        }
    sPp = sPp->next;
  }
  if (myrank == 0)
    cout << "NullShellPatch2::setupintintstuff ss_src completes" << endl;

  Pp = CPatL;
  while (Pp)
  {
    double llb[dim], uub[dim];
    if (Symmetry > 0)
      llb[2] = Pp->data->bbox[2] - 0.0001 * CDH[2];
    else
      llb[2] = Pp->data->bbox[2] + (hCS_width + 0.0001) * CDH[2];
    uub[2] = Pp->data->bbox[dim + 2] - (hCS_width + 0.0001) * CDH[2];
    for (int j = 0; j < 2; j++)
    {
      if (Symmetry > 1)
        llb[j] = Pp->data->bbox[j] - 0.0001 * CDH[j];
      else
        llb[j] = Pp->data->bbox[j] + (hCS_width + 0.0001) * CDH[j];
      uub[j] = Pp->data->bbox[dim + j] - (hCS_width + 0.0001) * CDH[j];
    }
    for (int iz = 0; iz < Pp->data->shape[2]; iz++)
      for (int iy = 0; iy < Pp->data->shape[1]; iy++)
        for (int ix = 0; ix < Pp->data->shape[0]; ix++)
        {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
          x = Pp->data->bbox[0] + ix * CDH[0];
          y = Pp->data->bbox[1] + iy * CDH[1];
          z = Pp->data->bbox[2] + iz * CDH[2];
#else
#ifdef Cell
          x = Pp->data->bbox[0] + (ix + 0.5) * CDH[0];
          y = Pp->data->bbox[1] + (iy + 0.5) * CDH[1];
          z = Pp->data->bbox[2] + (iz + 0.5) * CDH[2];
#else
#error Not define Vertex nor Cell
#endif
#endif
          if (x < llb[0] || x > uub[0] ||
              y < llb[1] || y > uub[1] ||
              z < llb[2] || z > uub[2])
          {
            int sst;
            double lx, ly, lz;
            bool flag = false;
            getlocalpox(x, y, z, sst, lx, ly, lz);
            for (int i = 0; i < cpusize; i++)
            {
              flag = prolongpointstru(cs_src[i], true, -1, PatL, DH, CPatL, CDH, x, y, z, Symmetry, i, -1);
              if (flag)
                break;
            }
            if (!flag)
            {
              if (myrank == 0)
              {
                cout << "ShellPatch2::prolongpointstru fail to find shell source point for" << endl;
                cout << "sst = -1, x,y,z = " << x << "," << y << "," << z << endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
              }
            }
          }
        }
    Pp = Pp->next;
  }
  if (myrank == 0)
    if (CPatL)
      cout << "NullShellPatch2::setupintintstuff cs_src completes" << endl;
    else
      cout << "NullShellPatch2::no cs_src exists" << endl;

  for (int i = 0; i < cpusize; i++)
  {
    ps = ss_src[i];
    while (ps)
    {
      ts = ps->next;
      prolongpointstru(ss_dst[i], PatL, DH, CPatL, CDH, ps); // ps may be insterted more here
      ps = ts;
    }

    if (CPatL)
    {
      ps = cs_src[i];
      while (ps)
      {
        ts = ps->next;
        prolongpointstru(cs_dst[i], PatL, DH, CPatL, CDH, ps); // ps may be insterted more here
        ps = ts;
      }
    }
  }
  if (myrank == 0)
    cout << "NullShellPatch2::setupintintstuff ss_dst and cs_dst complete" << endl;

  /*
    for(int i=0;i<cpusize;i++)
    {
       ps=ss_src[i];
       ts=ss_dst[i];
       while(ps)
       {
        if(myrank==0) cout<<"src:"<<endl;
        check_pointstrul(ps,1);
        if(myrank==0) cout<<"dst:"<<endl;
        check_pointstrul(ts,1);
        ps=ps->next;
        ts=ts->next;
       }
    }
    exit(0);
  */
}
// lz is x instead of r
void NullShellPatch2::getlocalpox(double x, double y, double z, int &sst, double &lx, double &ly, double &lz)
{
  double r;
  r = sqrt(x * x + y * y + z * z);
  lz = r / (r + Rmin);
  if (fabs(x) <= z && fabs(y) <= z)
  {
    sst = 0;
    lx = atan(x / z);
    ly = atan(y / z);
  }
  else if (fabs(x) <= -z && fabs(y) <= -z)
  {
    sst = 1;
    lx = atan(x / z);
    ly = atan(y / z);
  }
  else if (fabs(y) <= x && fabs(z) <= x)
  {
    sst = 2;
    lx = atan(y / x);
    ly = atan(z / x);
  }
  else if (fabs(y) <= -x && fabs(z) <= -x)
  {
    sst = 3;
    lx = atan(y / x);
    ly = atan(z / x);
  }
  else if (fabs(x) <= y && fabs(z) <= y)
  {
    sst = 4;
    lx = atan(x / y);
    ly = atan(z / y);
  }
  else if (fabs(x) <= -y && fabs(z) <= -y)
  {
    sst = 5;
    lx = atan(x / y);
    ly = atan(z / y);
  }
  else
  {
    cout << "NullShellPatch2::getlocalpox should not come here, something wrong" << endl;
  }
}
// lz is x instead of r
// using fake global coordinates to get local coordinate
void NullShellPatch2::getlocalpox_fake(double x, double y, double z, int &sst, double &lx, double &ly, double &lz)
{
  double r;
  r = sqrt(x * x + y * y + z * z);
  lz = r;
  if (fabs(x) <= z && fabs(y) <= z)
  {
    sst = 0;
    lx = atan(x / z);
    ly = atan(y / z);
  }
  else if (fabs(x) <= -z && fabs(y) <= -z)
  {
    sst = 1;
    lx = atan(x / z);
    ly = atan(y / z);
  }
  else if (fabs(y) <= x && fabs(z) <= x)
  {
    sst = 2;
    lx = atan(y / x);
    ly = atan(z / x);
  }
  else if (fabs(y) <= -x && fabs(z) <= -x)
  {
    sst = 3;
    lx = atan(y / x);
    ly = atan(z / x);
  }
  else if (fabs(x) <= y && fabs(z) <= y)
  {
    sst = 4;
    lx = atan(x / y);
    ly = atan(z / y);
  }
  else if (fabs(x) <= -y && fabs(z) <= -y)
  {
    sst = 5;
    lx = atan(x / y);
    ly = atan(z / y);
  }
  else
  {
    cout << "NullShellPatch2::getlocalpox should not come here, something wrong" << endl;
  }
}
// lz is x instead of r
// specially for usage from shell to shell
void NullShellPatch2::getlocalpox_ss(int isst, double ix, double iy, double iz, int &sst, double &lx, double &ly, double &lz)
{
  // fake global coordinate
  double r = 1, x, y, z;
  switch (isst)
  {
  case 0:
    x = tan(ix);
    y = tan(iy);
    z = r / sqrt(1 + x * x + y * y);
    x = z * x;
    y = z * y;
    break;
  case 1:
    x = tan(ix);
    y = tan(iy);
    z = -r / sqrt(1 + x * x + y * y);
    x = z * x;
    y = z * y;
    break;
  case 2:
    y = tan(ix);
    z = tan(iy);
    x = r / sqrt(1 + z * z + y * y);
    y = x * y;
    z = x * z;
    break;
  case 3:
    y = tan(ix);
    z = tan(iy);
    x = -r / sqrt(1 + z * z + y * y);
    y = x * y;
    z = x * z;
    break;
  case 4:
    x = tan(ix);
    z = tan(iy);
    y = r / sqrt(1 + x * x + z * z);
    x = y * x;
    z = y * z;
    break;
  case 5:
    x = tan(ix);
    z = tan(iy);
    y = -r / sqrt(1 + x * x + z * z);
    x = y * x;
    z = y * z;
    break;
  }

  // map with fake global coordinate
  if (fabs(x) <= z && fabs(y) <= z)
  {
    sst = 0;
    lx = atan(x / z);
    ly = atan(y / z);
  }
  else if (fabs(x) <= -z && fabs(y) <= -z)
  {
    sst = 1;
    lx = atan(x / z);
    ly = atan(y / z);
  }
  else if (fabs(y) <= x && fabs(z) <= x)
  {
    sst = 2;
    lx = atan(y / x);
    ly = atan(z / x);
  }
  else if (fabs(y) <= -x && fabs(z) <= -x)
  {
    sst = 3;
    lx = atan(y / x);
    ly = atan(z / x);
  }
  else if (fabs(x) <= y && fabs(z) <= y)
  {
    sst = 4;
    lx = atan(x / y);
    ly = atan(z / y);
  }
  else if (fabs(x) <= -y && fabs(z) <= -y)
  {
    sst = 5;
    lx = atan(x / y);
    ly = atan(z / y);
  }
  else
  {
    cout << "NullShellPatch2::getlocalpox should not come here, something wrong" << endl;
  }

  lz = iz;

  //     if(lx != lx) cout<<lx<<","<<ly<<","<<lz<<endl;
}
// lz is x instead of r
void NullShellPatch2::getlocalpoxsst(double x, double y, double z, int sst, double &lx, double &ly, double &lz)
{
  double r;
  r = sqrt(x * x + y * y + z * z);
  lz = r / (r + Rmin);
  switch (sst)
  {
  case -1:
    lx = x;
    ly = y;
    lz = z;
    break;
  case 0:
    lx = atan(x / z);
    ly = atan(y / z);
    break;
  case 1:
    lx = atan(x / z);
    ly = atan(y / z);
    break;
  case 2:
    lx = atan(y / x);
    ly = atan(z / x);
    break;
  case 3:
    lx = atan(y / x);
    ly = atan(z / x);
    break;
  case 4:
    lx = atan(x / y);
    ly = atan(z / y);
    break;
  case 5:
    lx = atan(x / y);
    ly = atan(z / y);
    break;
  default:
    cout << "NullShellPatch2::getlocalpoxsst should not come here, something wrong" << endl;
  }
}
// lz is x instead of r
// special for usage from shell to shell
void NullShellPatch2::getlocalpoxsst_ss(int isst, double ix, double iy, double iz, int lsst, double &lx, double &ly, double &lz)
{
  // fake global coordinate
  double r = 1, x, y, z;
  switch (isst)
  {
  case 0:
    x = tan(ix);
    y = tan(iy);
    z = r / sqrt(1 + x * x + y * y);
    x = z * x;
    y = z * y;
    break;
  case 1:
    x = tan(ix);
    y = tan(iy);
    z = -r / sqrt(1 + x * x + y * y);
    x = z * x;
    y = z * y;
    break;
  case 2:
    y = tan(ix);
    z = tan(iy);
    x = r / sqrt(1 + z * z + y * y);
    y = x * y;
    z = x * z;
    break;
  case 3:
    y = tan(ix);
    z = tan(iy);
    x = -r / sqrt(1 + z * z + y * y);
    y = x * y;
    z = x * z;
    break;
  case 4:
    x = tan(ix);
    z = tan(iy);
    y = r / sqrt(1 + x * x + z * z);
    x = y * x;
    z = y * z;
    break;
  case 5:
    x = tan(ix);
    z = tan(iy);
    y = -r / sqrt(1 + x * x + z * z);
    x = y * x;
    z = y * z;
    break;
  }

  // map with fake global coordinate
  switch (lsst)
  {
  case 0:
    lx = atan(x / z);
    ly = atan(y / z);
    break;
  case 1:
    lx = atan(x / z);
    ly = atan(y / z);
    break;
  case 2:
    lx = atan(y / x);
    ly = atan(z / x);
    break;
  case 3:
    lx = atan(y / x);
    ly = atan(z / x);
    break;
  case 4:
    lx = atan(x / y);
    ly = atan(z / y);
    break;
  case 5:
    lx = atan(x / y);
    ly = atan(z / y);
    break;
  default:
    cout << "NullShellPatch2::getlocalpoxsst_ss should not come here, something wrong" << endl;
  }

  lz = iz;
}
// lz is x instead of r
void NullShellPatch2::getglobalpox(double &x, double &y, double &z, int sst, double lx, double ly, double lz)
{
  double r = lz * Rmin / (1 - lz);
  switch (sst)
  {
  case 0:
    x = tan(lx);
    y = tan(ly);
    z = r / sqrt(1 + x * x + y * y);
    x = z * x;
    y = z * y;
    break;
  case 1:
    x = tan(lx);
    y = tan(ly);
    z = -r / sqrt(1 + x * x + y * y);
    x = z * x;
    y = z * y;
    break;
  case 2:
    y = tan(lx);
    z = tan(ly);
    x = r / sqrt(1 + z * z + y * y);
    y = x * y;
    z = x * z;
    break;
  case 3:
    y = tan(lx);
    z = tan(ly);
    x = -r / sqrt(1 + z * z + y * y);
    y = x * y;
    z = x * z;
    break;
  case 4:
    x = tan(lx);
    z = tan(ly);
    y = r / sqrt(1 + x * x + z * z);
    x = y * x;
    z = y * z;
    break;
  case 5:
    x = tan(lx);
    z = tan(ly);
    y = -r / sqrt(1 + x * x + z * z);
    x = y * x;
    z = y * z;
    break;
  }
}
void NullShellPatch2::checkBlock(int sst)
{
  if (myrank == 0)
  {
    cout << "checking shell patch sst = " << sst << endl;
    MyList<ss_patch> *Pp = PatL;
    while (Pp)
    {
      if (Pp->data->sst == sst)
      {
        MyList<Block> *BP = Pp->data->blb;
        while (BP)
        {
          BP->data->checkBlock();
          if (BP == Pp->data->ble)
            break;
          BP = BP->next;
        }
      }
      Pp = Pp->next;
    }
  }
}
void NullShellPatch2::check_pointstrul(MyList<pointstru> *pp, bool first_only)
{
  int myrank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (myrank == 0)
  {
    if (!pp)
      cout << "ShellPatch::check_pointstrul meets empty pointstru" << endl;
    else
      cout << "checking check_pointstrul..." << endl;
    while (pp)
    {
      if (pp->data->Bg)
        cout << "on node#" << pp->data->Bg->rank << endl;
      else
        cout << "virtual pointstru" << endl;
      cout << "source sst = " << pp->data->ssst << endl;
      cout << "target sst = " << pp->data->tsst << endl;
      cout << "dumy dimension = " << pp->data->dumyd << endl;
      cout << "global coordinates: (";
      for (int i = 0; i < dim; i++)
      {
        if (i < dim - 1)
          cout << pp->data->gpox[i] << ",";
        else
          cout << pp->data->gpox[i] << ")" << endl;
      }
      cout << "local coordinates: (";
      for (int i = 0; i < dim; i++)
      {
        if (i < dim - 1)
          cout << pp->data->lpox[i] << ",";
        else
          cout << pp->data->lpox[i] << ")" << endl;
      }
      if (first_only)
        return;
      pp = pp->next;
    }
  }
}
// for check
// used by _dst construction, so these x,y,z must coinside with grid point
// we have considered ghost points now
void NullShellPatch2::prolongpointstru(MyList<pointstru> *&psul, MyList<ss_patch> *sPpi, double DH[dim],
                                       MyList<Patch> *Ppi, double CDH[dim], MyList<pointstru> *pss)
{
  int n_dst = 0;
  MyList<ss_patch> *sPp = sPpi;
  MyList<Patch> *Pp = Ppi;
  MyList<Block> *Bgl;
  Block *Bg;
  double llb[dim], uub[dim];
  double lx, ly, lz, lsst;

  if (pss->data->tsst >= 0)
  {
    getlocalpoxsst(pss->data->gpox[0], pss->data->gpox[1], pss->data->gpox[2], pss->data->tsst,
                   lx, ly, lz);
    if (lx != lx)
      getlocalpoxsst_ss(pss->data->ssst, pss->data->lpox[0], pss->data->lpox[1], pss->data->lpox[2],
                        pss->data->tsst, lx, ly, lz);
    while (sPp)
    {
      if (sPp->data->sst == pss->data->tsst)
      {
        Bgl = sPp->data->blb;
        while (Bgl)
        {
          Bg = Bgl->data;
          {
            for (int j = 0; j < dim; j++)
            {
              llb[j] = Bg->bbox[j];
              uub[j] = Bg->bbox[j + dim];
            }

            if (lx > llb[0] - 0.1 * DH[0] && lx < uub[0] + 0.1 * DH[0] &&
                ly > llb[1] - 0.1 * DH[1] && ly < uub[1] + 0.1 * DH[1] &&
                lz > llb[2] - 0.1 * DH[2] && lz < uub[2] + 0.1 * DH[2])
            {
              MyList<pointstru> *ps = new MyList<pointstru>;
              ps->data = new pointstru;
              ps->next = 0;
              for (int i = 0; i < dim; i++)
                ps->data->gpox[i] = pss->data->gpox[i];
              ps->data->lpox[0] = lx;
              ps->data->lpox[1] = ly;
              ps->data->lpox[2] = lz;
              ps->data->ssst = pss->data->ssst;
              ps->data->tsst = sPp->data->sst;
              ps->data->dumyd = getdumydimension(ps->data->tsst, ps->data->ssst);
              ps->data->Bg = Bg;
              ps->data->coef = 0;
              ps->data->sind = 0;
              ps->data->indz = pss->data->indz;
              get_Jacob(ps->data->lpox, ps->data->tsst, ps->data->ssst, ps->data->Jacob);
              if (psul)
                psul->catList(ps);
              else
                psul = ps;
              n_dst++;
            }
          }
          if (Bgl == sPp->data->ble)
            break;
          Bgl = Bgl->next;
        }
      }
      sPp = sPp->next;
    }
  }
  else
  {
    if (pss->data->tsst != -1)
      cout << "somthing is wrong in NullShellPatch2::prolongpointstru" << endl;
    lx = pss->data->gpox[0];
    ly = pss->data->gpox[1];
    lz = pss->data->gpox[2];
    while (Pp)
    {
      Bgl = Pp->data->blb;
      while (Bgl)
      {
        Bg = Bgl->data;
        {
          for (int j = 0; j < dim; j++)
          {
            llb[j] = Bg->bbox[j];
            uub[j] = Bg->bbox[j + dim];
          }

          if (lx > llb[0] - 0.1 * CDH[0] && lx < uub[0] + 0.1 * CDH[0] &&
              ly > llb[1] - 0.1 * CDH[1] && ly < uub[1] + 0.1 * CDH[1] &&
              lz > llb[2] - 0.1 * CDH[2] && lz < uub[2] + 0.1 * CDH[2])
          {
            MyList<pointstru> *ps = new MyList<pointstru>;
            ps->data = new pointstru;
            ps->next = 0;
            for (int i = 0; i < dim; i++)
              ps->data->gpox[i] = pss->data->gpox[i];
            ps->data->lpox[0] = lx;
            ps->data->lpox[1] = ly;
            ps->data->lpox[2] = lz;
            ps->data->ssst = pss->data->ssst;
            ps->data->tsst = -1;
            ps->data->dumyd = getdumydimension(ps->data->tsst, ps->data->ssst);
            ps->data->Bg = Bg;
            ps->data->coef = 0;
            ps->data->sind = 0;
            ps->data->indz = pss->data->indz;
            for (int i = 0; i < 2; i++)
              for (int j = 0; j < 2; j++)
                ps->data->Jacob[i][j] = 0;
            if (psul)
              psul->catList(ps);
            else
              psul = ps;
            n_dst++;
          }
        }
        if (Bgl == Pp->data->ble)
          break;
        Bgl = Bgl->next;
      }
      Pp = Pp->next;
    }
  }
  // if n_dst > 0, that's because of ghost_points then prolong source list
  if (n_dst == 0)
  {
    int myrank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank == 0)
      cout << "NullShellPatch2::prolongpointstru fail to find target Block for pointstru:" << endl;
    check_pointstrul(pss, true);
    if (Pp == Ppi)
    {
      getlocalpoxsst(pss->data->gpox[0], pss->data->gpox[1], pss->data->gpox[2], pss->data->tsst,
                     lx, ly, lz);
      if (myrank == 0)
        cout << "sst = " << pss->data->tsst << ", lx,ly,lz = " << lx << "," << ly << "," << lz << endl;
      checkBlock(pss->data->tsst);
    }
    else
    {
      Pp = Ppi;
      while (Pp)
      {
        Pp->data->checkBlock();
        Pp = Pp->next;
      }
    }
    if (myrank == 0)
      MPI_Abort(MPI_COMM_WORLD, 1);
  }
  else
  {
    MyList<pointstru> *ts = 0;
    for (int i = 1; i < n_dst; i++)
    {
      MyList<pointstru> *ps = new MyList<pointstru>;
      ps->data = new pointstru;
      ps->next = (i == n_dst - 1) ? pss->next : 0;
      for (int i = 0; i < dim; i++)
      {
        ps->data->gpox[i] = pss->data->gpox[i];
        ps->data->lpox[i] = pss->data->lpox[i];
      }
      ps->data->ssst = pss->data->ssst;
      ps->data->tsst = pss->data->tsst;
      ps->data->dumyd = getdumydimension(ps->data->ssst, ps->data->tsst);
      ps->data->Bg = pss->data->Bg;
      ps->data->coef = 0;
      ps->data->sind = 0;
      ps->data->indz = pss->data->indz;
      for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
          ps->data->Jacob[i][j] = pss->data->Jacob[i][j];
      if (ts)
        ts->catList(ps);
      else
        ts = ps;
    }
    if (ts)
      pss->next = ts;
  }
}
// used by _src construction, so these x,y,z do not coinside with grid point
bool NullShellPatch2::prolongpointstru(MyList<pointstru> *&psul, bool ssyn, int tsst, MyList<ss_patch> *sPp, double DH[dim],
                                       MyList<Patch> *Pp, double CDH[dim], double x, double y, double z, int Symmetry, int rank_in,
                                       const int iz)
{
  MyList<Block> *Bgl;
  Block *Bg;
  double llb[dim], uub[dim];
  double lx, ly, lz;

  if (ssyn)
  {
    int sst;
    getlocalpox(x, y, z, sst, lx, ly, lz);
    while (sPp)
    {
      if (sPp->data->sst == sst)
      {
        Bgl = sPp->data->blb;
        while (Bgl)
        {
          Bg = Bgl->data;
          if (Bg->rank == rank_in)
          {
            for (int j = 0; j < 2; j++)
            {
              if (feq(Bg->bbox[j], -PI / 4 - overghost * DH[j], DH[j] / 2))
                llb[j] = -PI / 4;
              else if (feq(Bg->bbox[j], sPp->data->bbox[j], DH[j] / 2))
                llb[j] = Bg->bbox[j];
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
              else
                llb[j] = Bg->bbox[j] + (ghost_width - 1) * DH[j];
#else
#ifdef Cell
              else
                llb[j] = Bg->bbox[j] + ghost_width * DH[j];
#else
#error Not define Vertex nor Cell
#endif
#endif
              if (feq(Bg->bbox[dim + j], PI / 4 + overghost * DH[j], DH[j] / 2))
                uub[j] = PI / 4;
              else if (feq(Bg->bbox[dim + j], sPp->data->bbox[dim + j], DH[j] / 2))
                uub[j] = Bg->bbox[dim + j];
              else
                uub[j] = Bg->bbox[dim + j] - ghost_width * DH[j];
            }
            if (feq(Bg->bbox[2], sPp->data->bbox[2], DH[2] / 2))
              llb[2] = Bg->bbox[2];
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
            else
              llb[2] = Bg->bbox[2] + (ghost_width - 1) * DH[2];
#else
#ifdef Cell
            else
              llb[2] = Bg->bbox[2] + ghost_width * DH[2];
#else
#error Not define Vertex nor Cell
#endif
#endif
            if (feq(Bg->bbox[dim + 2], sPp->data->bbox[dim + 2], DH[2] / 2))
              uub[2] = Bg->bbox[dim + 2];
            else
              uub[2] = Bg->bbox[dim + 2] - ghost_width * DH[2];
            if (lx > llb[0] - 0.0001 * DH[0] && lx < uub[0] + 0.0001 * DH[0] &&
                ly > llb[1] - 0.0001 * DH[1] && ly < uub[1] + 0.0001 * DH[1] &&
                lz > llb[2] - 0.0001 * DH[2] && lz < uub[2] + 0.0001 * DH[2]) // even ghost_width-1 the region is like |----|----|
                                                                              //                                            ^
                                                                              // so for ^ point may miss for vertext center, so we use 0.0001
            {
              MyList<pointstru> *ps = new MyList<pointstru>;
              ps->data = new pointstru;
              ps->data->Bg = Bg;
              ps->data->gpox[0] = x;
              ps->data->gpox[1] = y;
              ps->data->gpox[2] = z;
              ps->data->lpox[0] = lx;
              ps->data->lpox[1] = ly;
              ps->data->lpox[2] = lz;
              ps->data->ssst = sPp->data->sst;
              ps->data->tsst = tsst;
              ps->data->dumyd = getdumydimension(ps->data->ssst, ps->data->tsst);
              ps->data->coef = 0;
              ps->data->sind = 0;
              ps->data->indz = iz;
              for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                  ps->data->Jacob[i][j] = 0;
              ps->next = 0;
              if (psul)
                psul->catList(ps);
              else
                psul = ps;
              return true;
            }
          }
          if (Bgl == sPp->data->ble)
            break;
          Bgl = Bgl->next;
        }
      }
      sPp = sPp->next;
    }
  }
  else
  {
    while (Pp)
    {
      Bgl = Pp->data->blb;
      while (Bgl)
      {
        Bg = Bgl->data;
        if (Bg->rank == rank_in)
        {
          for (int j = 0; j < dim; j++)
          {
            if (feq(Bg->bbox[j], Pp->data->bbox[j], CDH[j] / 2))
              llb[j] = Bg->bbox[j];
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
            else
              llb[j] = Bg->bbox[j] + (ghost_width - 1) * CDH[j];
#else
#ifdef Cell
            else
              llb[j] = Bg->bbox[j] + ghost_width * CDH[j];
#else
#error Not define Vertex nor Cell
#endif
#endif
            if (feq(Bg->bbox[dim + j], Pp->data->bbox[dim + j], CDH[j] / 2))
              uub[j] = Bg->bbox[dim + j];
            else
              uub[j] = Bg->bbox[dim + j] - ghost_width * CDH[j];
          }
          if (x > llb[0] - 0.0001 * CDH[0] && x < uub[0] + 0.0001 * CDH[0] &&
              y > llb[1] - 0.0001 * CDH[1] && y < uub[1] + 0.0001 * CDH[1] &&
              z > llb[2] - 0.0001 * CDH[2] && z < uub[2] + 0.0001 * CDH[2])
          {
            MyList<pointstru> *ps = new MyList<pointstru>;
            ps->data = new pointstru;
            ps->data->Bg = Bg;
            ps->data->gpox[0] = x;
            ps->data->gpox[1] = y;
            ps->data->gpox[2] = z;
            ps->data->lpox[0] = x;
            ps->data->lpox[1] = y;
            ps->data->lpox[2] = z;
            ps->data->ssst = -1;
            ps->data->tsst = tsst;
            ps->data->dumyd = getdumydimension(ps->data->ssst, ps->data->tsst);
            ps->data->coef = 0;
            ps->data->sind = 0;
            ps->data->indz = -1;
            for (int i = 0; i < 2; i++)
              for (int j = 0; j < 2; j++)
                ps->data->Jacob[i][j] = 0;
            ps->next = 0;
            if (psul)
              psul->catList(ps);
            else
              psul = ps;
            return true;
          }
        }
        if (Bgl == Pp->data->ble)
          break;
        Bgl = Bgl->next;
      }
      Pp = Pp->next;
    }
  }

  return false;
}
// used by _src construction, so these x,y,z do not coinside with grid point
// specially used from shell to shell
bool NullShellPatch2::prolongpointstru_ss(MyList<pointstru> *&psul, int tsst, MyList<ss_patch> *sPp, double DH[dim],
                                          MyList<Patch> *Pp, double CDH[dim], double x, double y, double z, int Symmetry, int rank_in, const int iz)
{
  MyList<Block> *Bgl;
  Block *Bg;
  double llb[dim], uub[dim];
  double lx, ly, lz;

  int sst;
  getlocalpox_ss(tsst, x, y, z, sst, lx, ly, lz);
  while (sPp)
  {
    if (sPp->data->sst == sst)
    {
      Bgl = sPp->data->blb;
      while (Bgl)
      {
        Bg = Bgl->data;
        if (Bg->rank == rank_in)
        {
          for (int j = 0; j < 2; j++)
          {
            if (feq(Bg->bbox[j], -PI / 4 - overghost * DH[j], DH[j] / 2))
              llb[j] = -PI / 4;
            else if (feq(Bg->bbox[j], sPp->data->bbox[j], DH[j] / 2))
              llb[j] = Bg->bbox[j];
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
            else
              llb[j] = Bg->bbox[j] + (ghost_width - 1) * DH[j];
#else
#ifdef Cell
            else
              llb[j] = Bg->bbox[j] + ghost_width * DH[j];
#else
#error Not define Vertex nor Cell
#endif
#endif
            if (feq(Bg->bbox[dim + j], PI / 4 + overghost * DH[j], DH[j] / 2))
              uub[j] = PI / 4;
            else if (feq(Bg->bbox[dim + j], sPp->data->bbox[dim + j], DH[j] / 2))
              uub[j] = Bg->bbox[dim + j];
            else
              uub[j] = Bg->bbox[dim + j] - ghost_width * DH[j];
          }
          if (feq(Bg->bbox[2], sPp->data->bbox[2], DH[2] / 2))
            llb[2] = Bg->bbox[2];
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
          else
            llb[2] = Bg->bbox[2] + (ghost_width - 1) * DH[2];
#else
#ifdef Cell
          else
            llb[2] = Bg->bbox[2] + ghost_width * DH[2];
#else
#error Not define Vertex nor Cell
#endif
#endif
          if (feq(Bg->bbox[dim + 2], sPp->data->bbox[dim + 2], DH[2] / 2))
            uub[2] = Bg->bbox[dim + 2];
          else
            uub[2] = Bg->bbox[dim + 2] - ghost_width * DH[2];
          if (lx > llb[0] - 0.0001 * DH[0] && lx < uub[0] + 0.0001 * DH[0] &&
              ly > llb[1] - 0.0001 * DH[1] && ly < uub[1] + 0.0001 * DH[1] &&
              lz > llb[2] - 0.0001 * DH[2] && lz < uub[2] + 0.0001 * DH[2]) // even ghost_width-1 the region is like |----|----|
                                                                            //                                            ^
                                                                            // so for ^ point may miss for vertext center, so we use 0.0001
          {
            MyList<pointstru> *ps = new MyList<pointstru>;
            ps->data = new pointstru;
            ps->data->Bg = Bg;
            ps->data->gpox[0] = 0; // global coordinate is not valid for r=infinity
            ps->data->gpox[1] = 0;
            ps->data->gpox[2] = 0;
            ps->data->lpox[0] = lx;
            ps->data->lpox[1] = ly;
            ps->data->lpox[2] = lz;
            ps->data->ssst = sPp->data->sst;
            ps->data->tsst = tsst;
            ps->data->dumyd = getdumydimension(ps->data->ssst, ps->data->tsst);
            ps->data->coef = 0;
            ps->data->sind = 0;
            ps->data->indz = iz;
            for (int i = 0; i < 2; i++)
              for (int j = 0; j < 2; j++)
                ps->data->Jacob[i][j] = 0;
            ps->next = 0;
            if (psul)
              psul->catList(ps);
            else
              psul = ps;
            return true;
          }
        }
        if (Bgl == sPp->data->ble)
          break;
        Bgl = Bgl->next;
      }
    }
    sPp = sPp->next;
  }

  return false;
}
// J[new][old] = d x_new/d x_old
void NullShellPatch2::get_Jacob(double *pox, int tsst, int ssst, double J[2][2])
{
  double rn = pox[0], sn = pox[1], ro, so;

  double cosro, sinro, cosso, sinso;
  if (tsst == 0 || tsst == 1) // z
  {
    if (ssst == 2 || ssst == 3) // x
    {
      ro = atan(tan(sn) / tan(rn));
      so = atan(1 / tan(rn));
      cosro = cos(ro);
      sinro = sin(ro);
      cosso = cos(so);
      sinso = sin(so);
      J[0][0] = 0;
      J[0][1] = -1;
      J[1][0] = cosso * cosso * sinro * sinro + cosro * cosro * sinso * sinso;
      J[1][1] = -cosro * sinro / J[1][0];
      J[1][0] = cosso * sinso / J[1][0];
    }
    else if (ssst == 4 || ssst == 5) // y
    {
      ro = atan(tan(rn) / tan(sn));
      so = atan(1 / tan(sn));
      cosro = cos(ro);
      sinro = sin(ro);
      cosso = cos(so);
      sinso = sin(so);
      J[0][0] = cosso * cosso * sinro * sinro + cosro * cosro * sinso * sinso;
      J[0][1] = -cosro * sinro / J[0][0];
      J[0][0] = cosso * sinso / J[0][0];
      J[1][0] = 0;
      J[1][1] = -1;
    }
    else
      cout << "Error in NullShellPatch2::get_Jacob 1" << endl;
  }
  else if (tsst == 2 || tsst == 3)
  {
    if (ssst == 0 || ssst == 1)
    {
      ro = atan(1 / tan(sn));
      so = atan(tan(rn) / tan(sn));
      cosro = cos(ro);
      sinro = sin(ro);
      cosso = cos(so);
      sinso = sin(so);
      J[0][0] = cosso * cosso * sinro * sinro + cosro * cosro * sinso * sinso;
      J[0][1] = cosro * sinro / J[0][0];
      J[0][0] = -cosso * sinso / J[0][0];
      J[1][0] = -1;
      J[1][1] = 0;
    }
    else if (ssst == 4 || ssst == 5)
    {
      ro = atan(1 / tan(rn));
      so = atan(tan(sn) / tan(rn));
      cosro = cos(ro);
      sinro = sin(ro);
      cosso = cos(so);
      sinso = sin(so);
      J[0][0] = -1;
      J[0][1] = 0;
      J[1][0] = cosso * cosso * sinro * sinro + cosro * cosro * sinso * sinso;
      J[1][1] = cosro * sinro / J[1][0];
      J[1][0] = -cosso * sinso / J[1][0];
    }
    else
      cout << "Error in NullShellPatch2::get_Jacob 2" << endl;
  }
  else if (tsst == 4 || tsst == 5)
  {
    if (ssst == 0 || ssst == 1)
    {
      ro = atan(tan(rn) / tan(sn));
      so = atan(1 / tan(sn));
      cosro = cos(ro);
      sinro = sin(ro);
      cosso = cos(so);
      sinso = sin(so);
      J[0][0] = cosso * cosso * sinro * sinro + cosro * cosro * sinso * sinso;
      J[0][1] = -cosro * sinro / J[0][0];
      J[0][0] = cosso * sinso / J[0][0];
      J[1][0] = 0;
      J[1][1] = -1;
    }
    else if (ssst == 2 || ssst == 3)
    {
      ro = atan(1 / tan(rn));
      so = atan(tan(sn) / tan(rn));
      cosro = cos(ro);
      sinro = sin(ro);
      cosso = cos(so);
      sinso = sin(so);
      J[0][0] = -1;
      J[0][1] = 0;
      J[1][0] = cosso * cosso * sinro * sinro + cosro * cosro * sinso * sinso;
      J[1][1] = cosro * sinro / J[1][0];
      J[1][0] = -cosso * sinso / J[1][0];
    }
    else
      cout << "Error in NullShellPatch2::get_Jacob 3" << endl;
  }
}
int NullShellPatch2::getdumydimension(int acsst, int posst) // -1 means no dumy dimension
{
  int dms;
  if (acsst == -1 || posst == -1)
    return -1;
  switch (acsst)
  {
  case 0:
  case 1:
    switch (posst)
    {
    case 0:
    case 1:
      cout << "error in NullShellPatch2::getdumydimension: acsst = " << acsst << ", posst = " << posst << endl;
      return -1;
    case 2:
    case 3:
      return 0;
    case 4:
    case 5:
      return 1;
    default:
      cout << "error in NullShellPatch2::getdumydimension: posst = " << posst << endl;
      return -1;
    }
  case 2:
  case 3:
    switch (posst)
    {
    case 0:
    case 1:
      return 1;
    case 2:
    case 3:
      cout << "error in NullShellPatch2::getdumydimension: acsst = " << acsst << ", posst = " << posst << endl;
      return -1;
    case 4:
    case 5:
      return 0;
    default:
      cout << "error in NullShellPatch2::getdumydimension: posst = " << posst << endl;
      return -1;
    }
  case 4:
  case 5:
    switch (posst)
    {
    case 0:
    case 1:
      return 1;
    case 2:
    case 3:
      return 0;
    case 4:
    case 5:
      cout << "error in NullShellPatch2::getdumydimension: acsst = " << acsst << ", posst = " << posst << endl;
      return -1;
    default:
      cout << "error in NullShellPatch2::getdumydimension: posst = " << posst << endl;
      return -1;
    }
  default:
    cout << "error in NullShellPatch2::getdumydimension: acsst = " << acsst << endl;
    return -1;
  }
}
void NullShellPatch2::Synch(MyList<var> *VarList, int Symmetry, double **Varwt, const short int svt)
{
  MyList<ss_patch> *Pp = PatL;
  while (Pp)
  {
    Pp->data->Sync(VarList, Symmetry);
    Pp = Pp->next;
  }

  // we need this before interpolation
  if (Symmetry > 0)
    fill_symmetric_boundarybuffer(VarList, Varwt);

  intertransfer(ss_src, ss_dst, VarList, VarList, Symmetry, Varwt, svt);

  // we need this here to correct conners
  if (Symmetry > 0)
    fill_symmetric_boundarybuffer(VarList, Varwt);
}
// Varwt: AoS of rho, sigma, x
void NullShellPatch2::fill_symmetric_boundarybuffer(MyList<var> *VarList, double **Varwt)
{
  MyList<var> *varl;
  int ind;
  double drho = getdX(0), dsigma = getdX(1);

  if (Symmetry == 0)
    return;
  else
  {
    MyList<ss_patch> *Pp = PatL;
    while (Pp)
    {
      MyList<Block> *BL = Pp->data->blb;
      while (BL)
      {
        Block *cg = BL->data;
        if (myrank == cg->rank)
        {
          varl = VarList;
          ind = 0;
          while (varl)
          {
            f_fill_symmetric_boundarybuffer2(cg->shape, cg->X[0], cg->X[1], cg->X[2], drho, dsigma,
                                             cg->fgfs[varl->data->sgfn],
                                             Symmetry, Pp->data->sst, Varwt[ind]); // defined in NullEvol2.f90
            varl = varl->next;
            ind++;
          }
        }
        if (BL == Pp->data->ble)
          break;
        BL = BL->next;
      }
      Pp = Pp->next;
    }
  }
}
void NullShellPatch2::intertransfer(MyList<pointstru> **src, MyList<pointstru> **dst,
                                    MyList<var> *VarList1 /* source */, MyList<var> *VarList2 /*target */,
                                    int Symmetry, double **Varwt, const short int svt)
{
  int myrank, cpusize;
  MPI_Comm_size(MPI_COMM_WORLD, &cpusize);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  int node;

  MPI_Request *reqs;
  MPI_Status *stats;
  reqs = new MPI_Request[2 * cpusize];
  stats = new MPI_Status[2 * cpusize];
  int req_no = 0;

  double **send_data, **rec_data;
  send_data = new double *[cpusize];
  rec_data = new double *[cpusize];
  int length;

  for (node = 0; node < cpusize; node++)
  {
    send_data[node] = rec_data[node] = 0;
    if (node == myrank)
    {
      if (length = interdata_packer(0, src[myrank], dst[myrank], node, PACK, VarList1, VarList2, Symmetry, Varwt, svt))
      {
        rec_data[node] = new double[length];
        if (!rec_data[node])
        {
          cout << "out of memory when new in short transfer, place 1" << endl;
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
        interdata_packer(rec_data[node], src[myrank], dst[myrank], node, PACK, VarList1, VarList2, Symmetry, Varwt, svt);
      }
    }
    else
    {
      // send from this cpu to cpu#node
      if (length = interdata_packer(0, src[myrank], dst[myrank], node, PACK, VarList1, VarList2, Symmetry, Varwt, svt))
      {
        send_data[node] = new double[length];
        if (!send_data[node])
        {
          cout << "out of memory when new in short transfer, place 2" << endl;
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
        interdata_packer(send_data[node], src[myrank], dst[myrank], node, PACK, VarList1, VarList2, Symmetry, Varwt, svt);
        MPI_Isend((void *)send_data[node], length, MPI_DOUBLE, node, 1, MPI_COMM_WORLD, reqs + req_no++);
      }
      // receive from cpu#node to this cpu
      if (length = interdata_packer(0, src[node], dst[node], node, UNPACK, VarList1, VarList2, Symmetry, Varwt, svt))
      {
        rec_data[node] = new double[length];
        if (!rec_data[node])
        {
          cout << "out of memory when new in short transfer, place 3" << endl;
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
        MPI_Irecv((void *)rec_data[node], length, MPI_DOUBLE, node, 1, MPI_COMM_WORLD, reqs + req_no++);
      }
    }
  }
  // wait for all requests to complete
  MPI_Waitall(req_no, reqs, stats);

  for (node = 0; node < cpusize; node++)
    if (rec_data[node])
      interdata_packer(rec_data[node], src[node], dst[node], node, UNPACK, VarList1, VarList2, Symmetry, Varwt, svt);

  for (node = 0; node < cpusize; node++)
  {
    if (send_data[node])
      delete[] send_data[node];
    if (rec_data[node])
      delete[] rec_data[node];
  }

  delete[] reqs;
  delete[] stats;
  delete[] send_data;
  delete[] rec_data;
}
//   PACK: prepare target data in 'data'
// UNPACK: copy target data from 'data' to corresponding numerical grids
int NullShellPatch2::interdata_packer(double *data, MyList<pointstru> *src, MyList<pointstru> *dst, int rank_in, int dir,
                                      MyList<var> *VarLists /* source */, MyList<var> *VarListd /* target */, int Symmetry, double **Varwt,
                                      const short int svt)
{
  int rev;
  rev = interdata_packer_pre(data, src, dst, rank_in, dir, VarLists, VarListd, Symmetry, Varwt, svt);
  if (dir == PACK)
    return rev;
  rev = interdata_packer_pot(data, src, dst, rank_in, dir, VarLists, VarListd, Symmetry, Varwt, svt);
  return rev;
}
int NullShellPatch2::interdata_packer_pre(double *data, MyList<pointstru> *src, MyList<pointstru> *dst, int rank_in, int dir,
                                          MyList<var> *VarLists /* source */, MyList<var> *VarListd /* target */, int Symmetry, double **Varwt,
                                          const short int svt)
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  int DIM = dim;
  int ordn = 2 * ghost_width;

  if (dir != PACK && dir != UNPACK)
  {
    cout << "error dir " << dir << " for data_packer " << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  int size_out = 0;

  if (!src || !dst)
    return size_out;

  MyList<var> *varls, *varld;

  varls = VarLists;
  varld = VarListd;
  while (varls && varld)
  {
    varls = varls->next;
    varld = varld->next;
  }

  if (varls || varld)
  {
    cout << "error in short data packer, var lists does not match." << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  while (src && dst)
  {
    if ((dir == PACK && dst->data->Bg->rank == rank_in && src->data->Bg->rank == myrank) ||
        (dir == UNPACK && src->data->Bg->rank == rank_in && dst->data->Bg->rank == myrank))
    {
      varls = VarLists;
      varld = VarListd;
      int vind = 1;
      bool flag = true;
      while (varls && varld)
      {
        if (data)
        {
          if (dir == PACK)
          {
            int DIMh = (src->data->dumyd == -1) ? dim : 1;
            if (src->data->coef == 0)
            {
              src->data->coef = new double[ordn * DIMh];
              src->data->sind = new int[dim];
              if (DIMh == 3)
              {
                for (int i = 0; i < DIMh; i++)
                {
                  double dd = src->data->Bg->getdX(i);
                  // 0.001 instead of 0.4 makes the point locate more center
                  src->data->sind[i] = int((src->data->lpox[i] - src->data->Bg->X[i][0]) / dd) - ordn / 2 + 1;
                  double h1, h2;
                  for (int j = 0; j < ordn; j++)
                  {
                    h1 = src->data->Bg->X[i][0] + (src->data->sind[i] + j) * dd;
                    src->data->coef[i * ordn + j] = 1;
                    for (int k = 0; k < j; k++)
                    {
                      h2 = src->data->Bg->X[i][0] + (src->data->sind[i] + k) * dd;
                      src->data->coef[i * ordn + j] *= (src->data->lpox[i] - h2) / (h1 - h2);
                    }
                    for (int k = j + 1; k < ordn; k++)
                    {
                      h2 = src->data->Bg->X[i][0] + (src->data->sind[i] + k) * dd;
                      src->data->coef[i * ordn + j] *= (src->data->lpox[i] - h2) / (h1 - h2);
                    }
                  }
                }
              }
              else
              {
                int actd = 1 - src->data->dumyd;
                double dd = src->data->Bg->getdX(actd);
                src->data->sind[0] = int((src->data->lpox[actd] - src->data->Bg->X[actd][0]) / dd) - ordn / 2 + 1;
                double h1, h2;
                for (int j = 0; j < ordn; j++)
                {
                  h1 = src->data->Bg->X[actd][0] + (src->data->sind[0] + j) * dd;
                  src->data->coef[j] = 1;
                  for (int k = 0; k < j; k++)
                  {
                    h2 = src->data->Bg->X[actd][0] + (src->data->sind[0] + k) * dd;
                    src->data->coef[j] *= (src->data->lpox[actd] - h2) / (h1 - h2);
                  }
                  for (int k = j + 1; k < ordn; k++)
                  {
                    h2 = src->data->Bg->X[actd][0] + (src->data->sind[0] + k) * dd;
                    src->data->coef[j] *= (src->data->lpox[actd] - h2) / (h1 - h2);
                  }
                }
                src->data->sind[2] = int((src->data->lpox[2] - src->data->Bg->X[2][0]) / src->data->Bg->getdX(2) + 0.001);
                if (!feq(src->data->Bg->X[2][src->data->sind[2]], src->data->lpox[2], src->data->Bg->getdX(2) / 2000))
                  cout << "error in NullShellPatch::interdata_packer point = " << src->data->lpox[2] << " != grid " << src->data->Bg->X[2][src->data->sind[2]] << endl;
                src->data->sind[1] = int((src->data->lpox[src->data->dumyd] - src->data->Bg->X[src->data->dumyd][0]) /
                                             src->data->Bg->getdX(src->data->dumyd) +
                                         0.001);
                if (!feq(src->data->Bg->X[src->data->dumyd][src->data->sind[1]], src->data->lpox[src->data->dumyd], src->data->Bg->getdX(src->data->dumyd) / 2000))
                  cout << "error in NullShellPatch::interdata_packer for dumy dimension point = "
                       << src->data->lpox[src->data->dumyd] << " != grid " << src->data->Bg->X[src->data->dumyd][src->data->sind[1]] << endl;
              }
            }
            // interpolate
            switch (DIMh)
            {
            case 3:
              f_global_interpind(src->data->Bg->shape, src->data->Bg->X[0], src->data->Bg->X[1], src->data->Bg->X[2],
                                 src->data->Bg->fgfs[varls->data->sgfn], data[size_out],
                                 src->data->lpox[0], src->data->lpox[1], src->data->lpox[2], ordn, varls->data->SoA, Symmetry,
                                 src->data->sind, src->data->coef, src->data->ssst);
              break;
            case 2:
              f_global_interpind2d(src->data->Bg->shape, src->data->Bg->X[0], src->data->Bg->X[1], src->data->Bg->X[2],
                                   src->data->Bg->fgfs[varls->data->sgfn], data[size_out],
                                   src->data->lpox[0], src->data->lpox[1], src->data->lpox[2], ordn, varls->data->SoA, Symmetry,
                                   src->data->sind, src->data->coef, src->data->ssst);
              break;
            case 1:
              f_global_interpind1d(src->data->Bg->shape, src->data->Bg->X[0], src->data->Bg->X[1], src->data->Bg->X[2],
                                   src->data->Bg->fgfs[varls->data->sgfn], data[size_out],
                                   src->data->lpox[0], src->data->lpox[1], src->data->lpox[2], ordn, varls->data->SoA, Symmetry,
                                   src->data->sind, src->data->coef, src->data->ssst, src->data->dumyd);
              break;
            default:
              cout << "NullShellPatch2::interdata_packer: not recognized DIM = " << DIMh << endl;
              MPI_Abort(MPI_COMM_WORLD, 1);
            }
          }
          if (dir == UNPACK) // from target data to corresponding grid
          {
            switch (svt)
            {
            case 1: // type(0,0)
              vind = 0;
              break;
            case 2: // type(0,1)
            {
              if (vind / 2 * 2 == vind)
              {
                double tmp[2];
                double Jon[2][2];
                Jon[0][0] = dst->data->Jacob[0][0];
                Jon[0][1] = dst->data->Jacob[0][1];
                Jon[1][0] = dst->data->Jacob[1][0];
                Jon[1][1] = dst->data->Jacob[1][1];

                tmp[0] = Jon[0][0] * Jon[1][1] - Jon[0][1] * Jon[1][0];
                tmp[1] = Jon[1][1] / tmp[0];
                Jon[0][1] = -Jon[0][1] / tmp[0];
                Jon[1][0] = -Jon[1][0] / tmp[0];
                Jon[1][1] = Jon[0][0] / tmp[0];
                Jon[0][0] = tmp[1];

                tmp[0] = data[size_out - 1];
                tmp[1] = data[size_out];
                data[size_out - 1] = Jon[0][0] * tmp[0] + Jon[1][0] * tmp[1];
                data[size_out] = Jon[0][1] * tmp[0] + Jon[1][1] * tmp[1];

                vind = 0;
              }
              break;
            }
            case 3: // symmetric type(0,2)
            {
              if (vind / 3 * 3 == vind)
              {
                double tmp[3];
                double Jon[2][2];
                Jon[0][0] = dst->data->Jacob[0][0];
                Jon[0][1] = dst->data->Jacob[0][1];
                Jon[1][0] = dst->data->Jacob[1][0];
                Jon[1][1] = dst->data->Jacob[1][1];
                tmp[0] = Jon[0][0] * Jon[1][1] - Jon[0][1] * Jon[1][0];
                tmp[1] = Jon[1][1] / tmp[0];
                Jon[0][1] = -Jon[0][1] / tmp[0];
                Jon[1][0] = -Jon[1][0] / tmp[0];
                Jon[1][1] = Jon[0][0] / tmp[0];
                Jon[0][0] = tmp[1];

                tmp[0] = data[size_out - 2];
                tmp[1] = data[size_out - 1];
                tmp[2] = data[size_out];
                data[size_out - 2] = Jon[0][0] * Jon[0][0] * tmp[0] + 2 * Jon[1][0] * Jon[0][0] * tmp[1] + Jon[1][0] * Jon[1][0] * tmp[2];
                data[size_out - 1] = Jon[0][0] * Jon[0][1] * tmp[0] + (Jon[1][0] * Jon[0][1] + Jon[0][0] * Jon[1][1]) * tmp[1] + Jon[1][0] * Jon[1][1] * tmp[2];
                data[size_out] = Jon[0][1] * Jon[0][1] * tmp[0] + 2 * Jon[1][1] * Jon[0][1] * tmp[1] + Jon[1][1] * Jon[1][1] * tmp[2];

                vind = 0;
              }
              break;
            }
            default:
            {
              cout << "NullShellPatch2::interdata_packer: not recognized svt = " << svt << endl;
              MPI_Abort(MPI_COMM_WORLD, 1);
            }
            }
          }
        }
        size_out += 1;
        vind += 1;
        varls = varls->next;
        varld = varld->next;
      }
    }
    dst = dst->next;
    src = src->next;
  }

  return size_out;
}
int NullShellPatch2::interdata_packer_pot(double *data, MyList<pointstru> *src, MyList<pointstru> *dst, int rank_in, int dir,
                                          MyList<var> *VarLists /* source */, MyList<var> *VarListd /* target */, int Symmetry, double **Varwt,
                                          const short int svt)
{
  if (dir != UNPACK)
    return 0;
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  int DIM = dim;
  int ordn = 2 * ghost_width;

  int size_out = 0;

  if (!src || !dst)
    return size_out;

  MyList<var> *varls, *varld;

  varls = VarLists;
  varld = VarListd;
  while (varls && varld)
  {
    varls = varls->next;
    varld = varld->next;
  }

  if (varls || varld)
  {
    cout << "error in short data packer, var lists does not match." << endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  while (src && dst)
  {
    if ((dir == UNPACK && src->data->Bg->rank == rank_in && dst->data->Bg->rank == myrank))
    {
      varls = VarLists;
      varld = VarListd;
      while (varls && varld)
      {
        if (data)
        {
          if (dir == UNPACK) // from target data to corresponding grid
          {
            f_pointcopy(DIM, dst->data->Bg->bbox, dst->data->Bg->bbox + dim, dst->data->Bg->shape, dst->data->Bg->fgfs[varld->data->sgfn],
                        dst->data->lpox[0], dst->data->lpox[1], dst->data->lpox[2], data[size_out]);
          }
        }
        size_out += 1;
        varls = varls->next;
        varld = varld->next;
      }
    }
    dst = dst->next;
    src = src->next;
  }

  return size_out;
}
void NullShellPatch2::Interp_Points_2D(MyList<var> *VarList,
                                       int NN, double **XX, /*input fake global Cartesian coordinate*/
                                       double *Shellf, int Symmetry)
{
  // NOTE: we do not Synchnize variables here, make sure of that before calling this routine
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  int ordn = 2 * ghost_width;
  MyList<var> *varl;
  int num_var = 0;
  varl = VarList;
  while (varl)
  {
    num_var++;
    varl = varl->next;
  }

  double *shellf;
  shellf = new double[NN * num_var];
  memset(shellf, 0, sizeof(double) * NN * num_var);

  // we use weight to monitor code, later some day we can move it for optimization
  int *weight;
  weight = new int[NN];
  memset(weight, 0, sizeof(int) * NN);

  double *DH, *llb, *uub;
  DH = new double[dim];

  for (int i = 0; i < dim; i++)
  {
    DH[i] = getdX(i);
  }
  llb = new double[dim];
  uub = new double[dim];

  for (int j = 0; j < NN; j++) // run along points
  {
    double pox[dim];
    int sst;
    getlocalpox_fake(XX[0][j], XX[1][j], XX[2][j], sst, pox[0], pox[1], pox[2]); // pox[2] is x indeed

    //    int indZ=int((pox[2]-xmin)/DH[2]);
    int indZ = shape[2]; // note we use index for Fortran
    MyList<ss_patch> *sPp = PatL;
    while (sPp->data->sst != sst)
      sPp = sPp->next;

    if (myrank == 0 && ((!sPp) || pox[2] < xmin - 0.0001 * DH[2] || pox[2] > xmax + 0.0001 * DH[2]))
    {
      cout << "NullShellPatch::Interp_Points: point gc = (";
      for (int k = 0; k < dim; k++)
      {
        cout << XX[k][j];
        if (k < dim - 1)
          cout << ",";
      }
      if (sPp)
      {
        cout << ") sst = " << sst << " lc = (";
        for (int k = 0; k < dim; k++)
        {
          cout << pox[k];
          if (k < dim - 1)
            cout << ",";
        }
      }
      cout << ") is out of the NullShellPatch." << endl;
      cout << "xmin = " << xmin << ", xmax = " << xmax << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (!sPp)
      return;

    MyList<Block> *Bp = sPp->data->blb;
    bool notfind = true;
    while (notfind && Bp) // run along Blocks
    {
      Block *BP = Bp->data;

      bool flag = true;
      for (int i = 0; i < dim; i++)
      {
// NOTE: our dividing structure is (exclude ghost)
// -1 0
//       1  2
// so (0,1) does not belong to any part for vertex structure
// here we put (0,0.5) to left part and (0.5,1) to right part
// BUT for cell structure the bbox is (-1.5,0.5) and (0.5,2.5), there is no missing region at all
//
// because of getlocalpox, pox will not goes into overghost region of ss_patch
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
        llb[i] = (feq(BP->bbox[i], sPp->data->bbox[i], DH[i] / 2)) ? BP->bbox[i] : BP->bbox[i] + (ghost_width - 0.5) * DH[i];
        uub[i] = (feq(BP->bbox[dim + i], sPp->data->bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] : BP->bbox[dim + i] - (ghost_width - 0.5) * DH[i];
#else
#ifdef Cell
        llb[i] = (feq(BP->bbox[i], sPp->data->bbox[i], DH[i] / 2)) ? BP->bbox[i] : BP->bbox[i] + ghost_width * DH[i];
        uub[i] = (feq(BP->bbox[dim + i], sPp->data->bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] : BP->bbox[dim + i] - ghost_width * DH[i];
#else
#error Not define Vertex nor Cell
#endif
#endif
        if (pox[i] - llb[i] < -DH[i] / 2 || pox[i] - uub[i] > DH[i] / 2)
        {
          flag = false;
          break;
        }
      }

      if (flag)
      {
        notfind = false;
        if (myrank == BP->rank)
        {
          //---> interpolation
          varl = VarList;
          int k = 0;
          while (varl) // run along variables
          {
            f_global_interp_ss_2d(BP->shape, BP->X[0], BP->X[1], indZ, BP->fgfs[varl->data->sgfn], shellf[j * num_var + k],
                                  pox[0], pox[1], ordn, varl->data->SoA, Symmetry, sst);
            varl = varl->next;
            k++;
          }
          weight[j] = 1;
        }
      }
      if (Bp == sPp->data->ble)
        break;
      Bp = Bp->next;
    }
  }

  MPI_Allreduce(shellf, Shellf, NN * num_var, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  int *Weight;
  Weight = new int[NN];
  MPI_Allreduce(weight, Weight, NN, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  for (int i = 0; i < NN; i++)
  {
    if (Weight[i] > 1)
    {
      if (myrank == 0)
        cout << "WARNING: NullShellPatch::Interp_Points meets multiple weight" << endl;
      for (int j = 0; j < num_var; j++)
        Shellf[j + i * num_var] = Shellf[j + i * num_var] / Weight[i];
    }
    else if (Weight[i] == 0 && myrank == 0)
    {
      cout << "ERROR: NullShellPatch::Interp_Points fails to find point (";
      for (int j = 0; j < dim; j++)
      {
        cout << XX[j][i];
        if (j < dim - 1)
          cout << ",";
        else
          cout << ")";
      }
      cout << " on NullShellPatch (" << xmin << ":" << xmax << ")" << endl;

      cout << "splited domains:" << endl;
      MyList<ss_patch> *sPp = PatL;
      while (sPp)
      {
        char sn[3];
        shellname(sn, sPp->data->sst);
        cout << "ss_patch " << sn << ":" << endl;
        MyList<Block> *Bp = sPp->data->blb;
        while (Bp)
        {
          Block *BP = Bp->data;

          for (int i = 0; i < dim; i++)
          {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
            llb[i] = (feq(BP->bbox[i], sPp->data->bbox[i], DH[i] / 2)) ? BP->bbox[i] : BP->bbox[i] + (ghost_width - 0.5) * DH[i];
            uub[i] = (feq(BP->bbox[dim + i], sPp->data->bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] : BP->bbox[dim + i] - (ghost_width - 0.5) * DH[i];
#else
#ifdef Cell
            llb[i] = (feq(BP->bbox[i], sPp->data->bbox[i], DH[i] / 2)) ? BP->bbox[i] : BP->bbox[i] + ghost_width * DH[i];
            uub[i] = (feq(BP->bbox[dim + i], sPp->data->bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] : BP->bbox[dim + i] - ghost_width * DH[i];
#else
#error Not define Vertex nor Cell
#endif
#endif
          }
          cout << "(";
          for (int j = 0; j < dim; j++)
          {
            cout << llb[j] << ":" << uub[j];
            if (j < dim - 1)
              cout << ",";
            else
              cout << ")" << endl;
          }
          if (Bp == sPp->data->ble)
            break;
          Bp = Bp->next;
        }
        sPp = sPp->next;
      }
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  delete[] shellf;
  delete[] weight;
  delete[] Weight;
  delete[] DH;
  delete[] llb;
  delete[] uub;
}
