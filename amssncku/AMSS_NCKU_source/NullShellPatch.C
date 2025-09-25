
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <new>
using namespace std;

#include "NullShellPatch.h"
#include "Parallel.h"
#include "fmisc.h"
#include "misc.h"
#include "shellfunctions.h"
#include "NullEvol.h"
#include "NullNews.h"
#include "initial_null.h"
#include "rungekutta4_rout.h"
#include "kodiss.h"

#define PI M_PI

// x   x   x   x   x   o   *
//             *   o   x   x   x   x   x
// each side contribute an overlap points
// so we need half of that
#define overghost ((ghost_width + 1) / 2 + ghost_width)

NullShellPatch::NullShellPatch(int *shapei, double Rmini, double xmini, double xmaxi, int Symmetryi, int myranki) : myrank(myranki), Rmin(Rmini), xmin(xmini), xmax(xmaxi), PatL(0), Symmetry(Symmetryi)
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
    cout << " null shell's range: r = [" << xmin * Rmin / (1 - xmin) << ":";
    if (xmax == 1)
      cout << " +Infty]" << endl;
    else
      cout << xmax * Rmin / (1 - xmax) << "]" << endl;
    cout << "                    x = [" << xmin << ":" << xmax << "]" << endl
         << " shape: " << shape[2] << endl
         << " resolution: [" << getdX(0) << "," << getdX(1) << "," << getdX(2) << "]" << endl;
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
  FXZEO = new var("FXZEO", ngfs++, 1, 1, 1);
  gx = new var("gx", ngfs++, 1, 1, 1);
  gy = new var("gy", ngfs++, 1, 1, 1);
  gz = new var("gz", ngfs++, 1, 1, 1);
  // every thing is taken as scalar
  beta = new var("beta", ngfs++, 1, 1, 1);
  W = new var("W", ngfs++, 1, 1, 1);
  KK = new var("KK", ngfs++, 1, 1, 1);
  HKK = new var("HKK", ngfs++, 1, 1, 1);
  KKx = new var("KKx", ngfs++, 1, 1, 1);
  HKKx = new var("HKKx", ngfs++, 1, 1, 1);
  Rnu = new var("Rnu", ngfs++, 1, 1, 1);
  Inu = new var("Inu", ngfs++, 1, 1, 1);
  Rk = new var("Rk", ngfs++, 1, 1, 1);
  Ik = new var("Ik", ngfs++, 1, 1, 1);
  RB = new var("RB", ngfs++, 1, 1, 1);
  IB = new var("IB", ngfs++, 1, 1, 1);
  RQ = new var("RQ", ngfs++, 1, 1, 1);
  IQ = new var("IQ", ngfs++, 1, 1, 1);
  RU = new var("RU", ngfs++, 1, 1, 1);
  IU = new var("IU", ngfs++, 1, 1, 1);
  RTheta = new var("RTheta", ngfs++, 1, 1, 1);
  ITheta = new var("ITheta", ngfs++, 1, 1, 1);
  RJo = new var("RJo", ngfs++, 1, 1, 1);
  IJo = new var("IJo", ngfs++, 1, 1, 1);
  omegao = new var("omegao", ngfs++, 1, 1, 1);
  RJ0 = new var("RJ0", ngfs++, 1, 1, 1);
  IJ0 = new var("IJ0", ngfs++, 1, 1, 1);
  omega0 = new var("omega0", ngfs++, 1, 1, 1);
  RJ = new var("RJ", ngfs++, 1, 1, 1);
  IJ = new var("IJ", ngfs++, 1, 1, 1);
  omega = new var("omega", ngfs++, 1, 1, 1);
  RJ1 = new var("RJ1", ngfs++, 1, 1, 1);
  IJ1 = new var("IJ1", ngfs++, 1, 1, 1);
  omega1 = new var("omega1", ngfs++, 1, 1, 1);
  RJ_rhs = new var("RJ_rhs", ngfs++, 1, 1, 1);
  IJ_rhs = new var("IJ_rhs", ngfs++, 1, 1, 1);
  omega_rhs = new var("omega_rhs", ngfs++, 1, 1, 1);

  quR1 = new var("quR1", ngfs++, 1, 1, 1);
  quI1 = new var("quI1", ngfs++, 1, 1, 1);
  quR2 = new var("quR2", ngfs++, 1, 1, 1);
  quI2 = new var("quI2", ngfs++, 1, 1, 1);
  qlR1 = new var("qlR1", ngfs++, 1, 1, 1);
  qlI1 = new var("qlI1", ngfs++, 1, 1, 1);
  qlR2 = new var("qlR2", ngfs++, 1, 1, 1);
  qlI2 = new var("qlI2", ngfs++, 1, 1, 1);
  gR = new var("gR", ngfs++, 1, 1, 1);
  gI = new var("gI", ngfs++, 1, 1, 1);

  dquR1 = new var("dquR1", ngfs++, 1, 1, 1);
  dquI1 = new var("dquI1", ngfs++, 1, 1, 1);
  dquR2 = new var("dquR2", ngfs++, 1, 1, 1);
  dquI2 = new var("dquI2", ngfs++, 1, 1, 1);
  bdquR1 = new var("bdquR1", ngfs++, 1, 1, 1);
  bdquI1 = new var("bdquI1", ngfs++, 1, 1, 1);
  bdquR2 = new var("bdquR2", ngfs++, 1, 1, 1);
  bdquI2 = new var("bdquI2", ngfs++, 1, 1, 1);
  dgR = new var("dgR", ngfs++, 1, 1, 1);
  dgI = new var("dgI", ngfs++, 1, 1, 1);
  bdgR = new var("bdgR", ngfs++, 1, 1, 1);
  bdgI = new var("bdgI", ngfs++, 1, 1, 1);

  RNews = new var("RNews", ngfs++, 1, 1, 1);
  INews = new var("INews", ngfs++, 1, 1, 1);

  DumpList = new MyList<var>(RJ0);
  DumpList->insert(IJ0);

  betaList = new MyList<var>(beta);
  betaList->insert(beta);
  betawt[0] = 0;
  QUList = new MyList<var>(RQ);
  QUList->insert(IQ);
  QUList->insert(RU);
  QUList->insert(IU);
  QUwt[0] = QUwt[1] = 1;
  WTheList = new MyList<var>(W);
  WTheList->insert(W);
  WTheList->insert(RTheta);
  WTheList->insert(ITheta);
  WThewt[0] = 0;
  WThewt[1] = 2;

  TheList = new MyList<var>(RTheta);
  TheList->insert(ITheta);

  OldStateList = new MyList<var>(RJo);
  OldStateList->insert(IJo);
  OldStateList->insert(omegao);
  StateList = new MyList<var>(RJ0);
  StateList->insert(IJ0);
  StateList->insert(omega0);
  SynchList_pre = new MyList<var>(RJ);
  SynchList_pre->insert(IJ);
  SynchList_pre->insert(omega);
  RHSList = new MyList<var>(RJ_rhs);
  RHSList->insert(IJ_rhs);
  RHSList->insert(omega_rhs);
  SynchList_cor = new MyList<var>(RJ1);
  SynchList_cor->insert(IJ1);
  SynchList_cor->insert(omega1);

  JrhsList = new MyList<var>(RJ_rhs);
  JrhsList->insert(IJ_rhs);
  J1List = new MyList<var>(RJ1);
  J1List->insert(IJ1);

  ingfs = 0;
  fngfs = ngfs;
}
NullShellPatch::~NullShellPatch()
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
  betaList->clearList();
  QUList->clearList();
  WTheList->clearList();
  TheList->clearList();
  JrhsList->clearList();
  J1List->clearList();

  delete FXZEO;
  delete gx;
  delete gy;
  delete gz;
  delete beta;
  delete W;
  delete Rnu;
  delete Inu;
  delete Rk;
  delete Ik;
  delete RB;
  delete IB;
  delete RQ;
  delete IQ;
  delete RU;
  delete IU;
  delete RTheta;
  delete ITheta;
  delete KK;
  delete HKK;
  delete KKx;
  delete HKKx;

  delete RJo;
  delete IJo;
  delete omegao;
  delete RJ0;
  delete IJ0;
  delete omega0;
  delete RJ;
  delete IJ;
  delete omega;
  delete RJ1;
  delete IJ1;
  delete omega1;
  delete RJ_rhs;
  delete IJ_rhs;
  delete omega_rhs;

  delete quR1;
  delete quR2;
  delete quI1;
  delete quI2;
  delete qlR1;
  delete qlR2;
  delete qlI1;
  delete qlI2;
  delete gR;
  delete gI;
  delete dquR1;
  delete dquR2;
  delete dquI1;
  delete dquI2;
  delete bdquR1;
  delete bdquR2;
  delete bdquI1;
  delete bdquI2;
  delete dgR;
  delete dgI;
  delete bdgR;
  delete bdgI;

  delete RNews;
  delete INews;
}
void NullShellPatch::destroypsuList(MyList<pointstru> *ct)
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
// the number of VarList = 2* the number of Varwt
void NullShellPatch::fill_symmetric_boundarybuffer(MyList<var> *VarList, int *Varwt)
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
            f_fill_symmetric_boundarybuffer(cg->shape, cg->X[0], cg->X[1], cg->X[2], drho, dsigma,
                                            cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                                            cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                                            cg->fgfs[varl->data->sgfn], cg->fgfs[varl->next->data->sgfn],
                                            Symmetry, Pp->data->sst, Varwt[ind]);
            varl = varl->next;
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
MyList<Block> *NullShellPatch::compose_sh(int cpusize)
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
int NullShellPatch::getdumydimension(int acsst, int posst) // -1 means no dumy dimension
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
      cout << "error in NullShellPatch::getdumydimension: acsst = " << acsst << ", posst = " << posst << endl;
      return -1;
    case 2:
    case 3:
      return 0;
    case 4:
    case 5:
      return 1;
    default:
      cout << "error in NullShellPatch::getdumydimension: posst = " << posst << endl;
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
      cout << "error in NullShellPatch::getdumydimension: acsst = " << acsst << ", posst = " << posst << endl;
      return -1;
    case 4:
    case 5:
      return 0;
    default:
      cout << "error in NullShellPatch::getdumydimension: posst = " << posst << endl;
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
      cout << "error in NullShellPatch::getdumydimension: acsst = " << acsst << ", posst = " << posst << endl;
      return -1;
    default:
      cout << "error in NullShellPatch::getdumydimension: posst = " << posst << endl;
      return -1;
    }
  default:
    cout << "error in NullShellPatch::getdumydimension: acsst = " << acsst << endl;
    return -1;
  }
}
void NullShellPatch::Setup_dyad()
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
        f_setup_dyad(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                     cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                     cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                     cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn],
                     cg->fgfs[dquR1->sgfn], cg->fgfs[dquR2->sgfn], cg->fgfs[dquI1->sgfn], cg->fgfs[dquI2->sgfn],
                     cg->fgfs[bdquR1->sgfn], cg->fgfs[bdquR2->sgfn], cg->fgfs[bdquI1->sgfn], cg->fgfs[bdquI2->sgfn],
                     cg->fgfs[dgR->sgfn], cg->fgfs[dgI->sgfn], cg->fgfs[bdgR->sgfn], cg->fgfs[bdgI->sgfn],
                     cg->fgfs[gx->sgfn], cg->fgfs[gy->sgfn], cg->fgfs[gz->sgfn],
                     Pp->data->sst, Rmin);
      }
      if (BL == Pp->data->ble)
        break;
      BL = BL->next;
    }
    Pp = Pp->next;
  }
}
void NullShellPatch::Setup_Initial_Data(bool checkrun, double PhysTime)
{
  if (checkrun)
  {
  }
  else
  {
    double one = 1.0;
    MyList<ss_patch> *Pp = PatL;
    while (Pp)
    {
      MyList<Block> *BL = Pp->data->blb;
      while (BL)
      {
        Block *cg = BL->data;
        if (myrank == cg->rank)
        {

          f_get_exact_null(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                           cg->fgfs[RJ0->sgfn], cg->fgfs[IJ0->sgfn], Pp->data->sst, Rmin, PhysTime,
                           cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                           cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                           cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn],
                           cg->fgfs[dquR1->sgfn], cg->fgfs[dquR2->sgfn], cg->fgfs[dquI1->sgfn], cg->fgfs[dquI2->sgfn],
                           cg->fgfs[bdquR1->sgfn], cg->fgfs[bdquR2->sgfn], cg->fgfs[bdquI1->sgfn], cg->fgfs[bdquI2->sgfn],
                           cg->fgfs[dgR->sgfn], cg->fgfs[dgI->sgfn], cg->fgfs[bdgR->sgfn], cg->fgfs[bdgI->sgfn]);
          //	   f_get_initial_null(cg->shape,cg->X[0],cg->X[1],cg->X[2],
          //                              cg->fgfs[RJ0->sgfn],cg->fgfs[IJ0->sgfn],Pp->data->sst,Rmin);
          //	   f_set_value(cg->shape,cg->fgfs[omega0->sgfn],one);
          f_get_exact_omega(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                            cg->fgfs[omega0->sgfn], Pp->data->sst, Rmin, PhysTime);
        }
        if (BL == Pp->data->ble)
          break;
        BL = BL->next;
      }
      Pp = Pp->next;
    }
    int Varwt[1];
    MyList<var> *DG_List;
#if 0     
     eth_derivs(RJ0,IJ0,RJ1,IJ1,0,1);     
     Varwt[0]=1; 
     DG_List=new MyList<var>(RJ1); DG_List->insert(IJ1);
     Synch(DG_List,Symmetry,Varwt);
     eth_derivs(RJ1,IJ1,RJ0,IJ0,1,1); 
     DG_List->clearList();  // after this DG_List = 0
#elif 0
    eth_dderivs(RJ1, IJ1, RJ0, IJ0, 0, 1, 1);
#endif
    DG_List = new MyList<var>(RJ0);
    DG_List->insert(IJ0);
    Varwt[0] = 2;
    Synch(DG_List, Symmetry, Varwt);

    Dump_Data(DG_List, 0, 0, 1);
    DG_List->clearList();
  }
}
void NullShellPatch::eth_derivs(var *Rv, var *Iv, var *ethRv, var *ethIv, int s, int e)
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
        f_eth_derivs(cg->shape, cg->X[0], cg->X[1], cg->fgfs[Rv->sgfn], cg->fgfs[Iv->sgfn],
                     cg->fgfs[ethRv->sgfn], cg->fgfs[ethIv->sgfn], s, e,
                     cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                     cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn]);
      }
      if (BL == Pp->data->ble)
        break;
      BL = BL->next;
    }
    Pp = Pp->next;
  }

  int Varwt[1];
  MyList<var> *DG_List;
  DG_List = new MyList<var>(ethRv);
  DG_List->insert(ethIv);
  Varwt[0] = s + e;
  Synch(DG_List, Symmetry, Varwt);
  DG_List->clearList();
}
void NullShellPatch::eth_dderivs(var *Rv, var *Iv, var *ethRv, var *ethIv, int s, int e1, int e2)
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
        f_eth_dderivs(cg->shape, cg->X[0], cg->X[1], cg->fgfs[Rv->sgfn], cg->fgfs[Iv->sgfn],
                      cg->fgfs[ethRv->sgfn], cg->fgfs[ethIv->sgfn], s, e1, e2,
                      cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                      cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn],
                      cg->fgfs[dquR1->sgfn], cg->fgfs[dquR2->sgfn], cg->fgfs[dquI1->sgfn], cg->fgfs[dquI2->sgfn],
                      cg->fgfs[bdquR1->sgfn], cg->fgfs[bdquR2->sgfn], cg->fgfs[bdquI1->sgfn], cg->fgfs[bdquI2->sgfn],
                      cg->fgfs[dgR->sgfn], cg->fgfs[dgI->sgfn], cg->fgfs[bdgR->sgfn], cg->fgfs[bdgI->sgfn]);
      }
      if (BL == Pp->data->ble)
        break;
      BL = BL->next;
    }
    Pp = Pp->next;
  }
  int Varwt[1];
  MyList<var> *DG_List;
  DG_List = new MyList<var>(ethRv);
  DG_List->insert(ethIv);
  Varwt[0] = s + e1 + e2;
  Synch(DG_List, Symmetry, Varwt);
  DG_List->clearList();
}
// lz is x instead of r
void NullShellPatch::getlocalpox(double x, double y, double z, int &sst, double &lx, double &ly, double &lz)
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
    cout << "NullShellPatch::getlocalpox should not come here, something wrong" << endl;
  }
}
// lz is x instead of r
// using fake global coordinates to get local coordinate
void NullShellPatch::getlocalpox_fake(double x, double y, double z, int &sst, double &lx, double &ly, double &lz)
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
    cout << "NullShellPatch::getlocalpox should not come here, something wrong" << endl;
  }
}
// lz is x instead of r
// specially for usage from shell to shell
void NullShellPatch::getlocalpox_ss(int isst, double ix, double iy, double iz, int &sst, double &lx, double &ly, double &lz)
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
    cout << "NullShellPatch::getlocalpox should not come here, something wrong" << endl;
  }

  lz = iz;

  //     if(lx != lx) cout<<lx<<","<<ly<<","<<lz<<endl;
}
// lz is x instead of r
void NullShellPatch::getlocalpoxsst(double x, double y, double z, int sst, double &lx, double &ly, double &lz)
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
    cout << "NullShellPatch::getlocalpoxsst should not come here, something wrong" << endl;
  }
}
// lz is x instead of r
// special for usage from shell to shell
void NullShellPatch::getlocalpoxsst_ss(int isst, double ix, double iy, double iz, int lsst, double &lx, double &ly, double &lz)
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
    cout << "NullShellPatch::getlocalpoxsst_ss should not come here, something wrong" << endl;
  }

  lz = iz;
}
// lz is x instead of r
void NullShellPatch::getglobalpox(double &x, double &y, double &z, int sst, double lx, double ly, double lz)
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
// we assume pox is the coordinate of target point
#if 1
complex<double> NullShellPatch::get_swtf(double *pox, int tsst, int ssst)
{
  double rn = pox[0], sn = pox[1], ro, so;
  double tcn, tsn, tco, tso;
  tcn = sqrt((1 - sin(rn) * sin(sn)) / 2);
  tsn = sqrt((1 + sin(rn) * sin(sn)) / 2);
  // upper a
  complex<double> qan[2];
  qan[0] = complex<double>(tsn, tcn);
  qan[1] = complex<double>(tsn, -tcn);
  qan[0] = 2.0 * tcn * tsn / cos(sn) * qan[0];
  qan[1] = 2.0 * tcn * tsn / cos(rn) * qan[1];
  if (tsst == 1 || tsst == 3 || tsst == 4)
  {
    qan[0] = conj(qan[0]);
    qan[1] = conj(qan[1]);
  }

  complex<double> qao[2];
  complex<double> gont;

  double J[2][2];
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
      tco = sqrt((1 - sin(ro) * sin(so)) / 2);
      tso = sqrt((1 + sin(ro) * sin(so)) / 2);
      // upper a
      qao[1] = complex<double>(tso, -tco);
      qao[1] = 2.0 * tco * tso / cos(ro) * qao[1];
      if (ssst == 1 || ssst == 3 || ssst == 4)
      {
        qao[1] = conj(qao[1]);
      }
      gont = -qan[0] / qao[1];
    }
    else if (ssst == 4 || ssst == 5) // y
    {
      ro = atan(tan(rn) / tan(sn));
      so = atan(1 / tan(sn));
      cosro = cos(ro);
      sinro = sin(ro);
      cosso = cos(so);
      sinso = sin(so);
      tco = sqrt((1 - sin(ro) * sin(so)) / 2);
      tso = sqrt((1 + sin(ro) * sin(so)) / 2);
      // upper a
      qao[1] = complex<double>(tso, -tco);
      qao[1] = 2.0 * tco * tso / cos(ro) * qao[1];
      if (ssst == 1 || ssst == 3 || ssst == 4)
      {
        qao[1] = conj(qao[1]);
      }
      gont = -qan[1] / qao[1];
    }
    else
      cout << "Error in NullShellPatch::get_swtf 1" << endl;
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
      tco = sqrt((1 - sin(ro) * sin(so)) / 2);
      tso = sqrt((1 + sin(ro) * sin(so)) / 2);
      // upper a
      qao[0] = complex<double>(tso, tco);
      qao[0] = 2.0 * tco * tso / cos(so) * qao[0];
      if (ssst == 1 || ssst == 3 || ssst == 4)
      {
        qao[0] = conj(qao[0]);
      }
      gont = -qan[1] / qao[0];
    }
    else if (ssst == 4 || ssst == 5)
    {
      ro = atan(1 / tan(rn));
      so = atan(tan(sn) / tan(rn));
      cosro = cos(ro);
      sinro = sin(ro);
      cosso = cos(so);
      sinso = sin(so);
      tco = sqrt((1 - sin(ro) * sin(so)) / 2);
      tso = sqrt((1 + sin(ro) * sin(so)) / 2);
      // upper a
      qao[0] = complex<double>(tso, tco);
      qao[0] = 2.0 * tco * tso / cos(so) * qao[0];
      if (ssst == 1 || ssst == 3 || ssst == 4)
      {
        qao[0] = conj(qao[0]);
      }
      gont = -qan[0] / qao[0];
    }
    else
      cout << "Error in NullShellPatch::get_swtf 2" << endl;
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
      tco = sqrt((1 - sin(ro) * sin(so)) / 2);
      tso = sqrt((1 + sin(ro) * sin(so)) / 2);
      // upper a
      qao[1] = complex<double>(tso, -tco);
      qao[1] = 2.0 * tco * tso / cos(ro) * qao[1];
      if (ssst == 1 || ssst == 3 || ssst == 4)
      {
        qao[1] = conj(qao[1]);
      }
      gont = -qan[1] / qao[1];
    }
    else if (ssst == 2 || ssst == 3)
    {
      ro = atan(1 / tan(rn));
      so = atan(tan(sn) / tan(rn));
      cosro = cos(ro);
      sinro = sin(ro);
      cosso = cos(so);
      sinso = sin(so);
      tco = sqrt((1 - sin(ro) * sin(so)) / 2);
      tso = sqrt((1 + sin(ro) * sin(so)) / 2);
      // upper a
      qao[0] = complex<double>(tso, tco);
      qao[0] = 2.0 * tco * tso / cos(so) * qao[0];
      if (ssst == 1 || ssst == 3 || ssst == 4)
      {
        qao[0] = conj(qao[0]);
      }
      gont = -qan[0] / qao[0];
    }
    else
      cout << "Error in NullShellPatch::get_swtf 3" << endl;
  }

  return gont;
}
#else
// #define DEBUG
complex<double> NullShellPatch::get_swtf(double *pox, int tsst, int ssst)
{
  double rn = pox[0], sn = pox[1], ro, so;
  double tcn, tsn, tco, tso;
  tcn = sqrt((1 - sin(rn) * sin(sn)) / 2);
  tsn = sqrt((1 + sin(rn) * sin(sn)) / 2);
#ifdef DEBUG
  // upper a
  complex<double> qan[2];
  qan[0] = complex<double>(tsn, tcn);
  qan[1] = complex<double>(tsn, -tcn);
  qan[0] = 2.0 * tcn * tsn / cos(sn) * qan[0];
  qan[1] = 2.0 * tcn * tsn / cos(rn) * qan[1];
  if (tsst == 1 || tsst == 3 || tsst == 4)
  {
    qan[0] = conj(qan[0]);
    qan[1] = conj(qan[1]);
  }
#endif
  // lower bar a
  complex<double> lan[2];
  lan[0] = complex<double>(tcn, -tsn);
  lan[1] = complex<double>(tcn, tsn);
  lan[0] = cos(sn) / 4.0 / tcn / tcn / tsn / tsn * lan[0];
  lan[1] = cos(rn) / 4.0 / tcn / tcn / tsn / tsn * lan[1];

  if (tsst == 1 || tsst == 3 || tsst == 4)
  {
    lan[0] = conj(lan[0]);
    lan[1] = conj(lan[1]);
  }

  complex<double> gont = complex<double>(2, 0);

  double J[2][2];
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
      cout << "Error in NullShellPatch::get_swtf 1" << endl;
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
      cout << "Error in NullShellPatch::get_swtf 2" << endl;
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
      cout << "Error in NullShellPatch::get_swtf 3" << endl;
  }
  tco = sqrt((1 - sin(ro) * sin(so)) / 2);
  tso = sqrt((1 + sin(ro) * sin(so)) / 2);

  complex<double> qao[2];
  // upper a
  qao[0] = complex<double>(tso, tco);
  qao[1] = complex<double>(tso, -tco);
  qao[0] = 2.0 * tco * tso / cos(so) * qao[0];
  qao[1] = 2.0 * tco * tso / cos(ro) * qao[1];
  if (ssst == 1 || ssst == 3 || ssst == 4)
  {
    qao[0] = conj(qao[0]);
    qao[1] = conj(qao[1]);
  }

  gont /= J[0][0] * lan[0] * qao[0] + J[0][1] * lan[0] * qao[1] + J[1][0] * lan[1] * qao[0] + J[1][1] * lan[1] * qao[1];

#ifdef DEBUG

  complex<double> lao[2];
  // lower bar a
  lao[0] = complex<double>(tco, -tso);
  lao[1] = complex<double>(tco, tso);
  lao[0] = cos(so) / 4.0 / tco / tco / tso / tso * lao[0];
  lao[1] = cos(ro) / 4.0 / tco / tco / tso / tso * lao[1];
  if (ssst == 1 || ssst == 3 || ssst == 4)
  {
    lao[0] = conj(lao[0]);
    lao[1] = conj(lao[1]);
  }

  static bool f1 = true, f2 = true, f3 = true, f4 = true;
  static bool f5 = true, f6 = true, f7 = true, f8 = true;
  static bool f9 = true, f10 = true, f11 = true, f12 = true;
  double hn11, hn12, hn22;
  double ho11, ho12, ho22;
  if (f1 && tsst == 0 && ssst == 2)
  {
    cout << "tsst = " << tsst << ", ssst = " << ssst << endl;
    cout << tan(rn) << "," << tan(sn) << "," << lan[0] * qan[0] + lan[1] * qan[1] << endl;
    cout << tan(ro) << "," << tan(so) << "," << lao[0] * qao[0] + lao[1] * qao[1] << endl;
    cout << "x+ -> z+; g -> x+; g -> z+" << endl;
    double the = atan(sqrt(tan(rn) * tan(rn) + tan(sn) * tan(sn))), phi = atan2(tan(sn), tan(rn));
    cout << (J[0][0] * qao[0] + J[0][1] * qao[1]) / qan[0] << ","
         << complex<double>(cos(phi), sin(phi) * cos(the)) / (cos(phi) * cos(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / (J[0][0] * qao[0] + J[0][1] * qao[1]) << ","
         << complex<double>(cos(phi), sin(phi) * cos(the)) / (cos(phi) * cos(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / qan[0]
         << endl;
    cout << (J[1][0] * qao[0] + J[1][1] * qao[1]) / qan[1] << ","
         << complex<double>(sin(phi), -cos(phi) * cos(the)) / (sin(phi) * sin(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / (J[1][0] * qao[0] + J[1][1] * qao[1]) << ","
         << complex<double>(sin(phi), -cos(phi) * cos(the)) / (sin(phi) * sin(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / qan[1]
         << endl;
    ho11 = pow(1 - sin(ro) * sin(ro) * sin(so) * sin(so), -2);
    ho12 = -0.25 * sin(2 * ro) * sin(2 * so) * ho11;
    ho22 = cos(ro) * cos(ro) * ho11;
    ho11 = cos(so) * cos(so) * ho11;
    hn11 = pow(1 - sin(rn) * sin(rn) * sin(sn) * sin(sn), -2);
    hn12 = -0.25 * sin(2 * rn) * sin(2 * sn) * hn11;
    hn22 = cos(rn) * cos(rn) * hn11;
    hn11 = cos(sn) * cos(sn) * hn11;
    cout << ho11 << "," << ho12 << "," << ho22 << endl;
    cout << hn11 * J[0][0] * J[0][0] + hn12 * J[0][0] * J[1][0] + hn12 * J[1][0] * J[0][0] + hn22 * J[1][0] * J[1][0] << ","
         << hn11 * J[0][0] * J[0][1] + hn12 * J[0][0] * J[1][1] + hn12 * J[1][0] * J[0][1] + hn22 * J[1][0] * J[1][1] << ","
         << hn11 * J[0][1] * J[0][1] + hn12 * J[0][1] * J[1][1] + hn12 * J[1][1] * J[0][1] + hn22 * J[1][1] * J[1][1] << endl;
    cout << "swtf = " << gont << endl;
    f1 = false;
  }
  else if (f2 && tsst == 0 && ssst == 3)
  {
    cout << "tsst = " << tsst << ", ssst = " << ssst << endl;
    cout << tan(rn) << "," << tan(sn) << "," << lan[0] * qan[0] + lan[1] * qan[1] << endl;
    cout << tan(ro) << "," << tan(so) << "," << lao[0] * qao[0] + lao[1] * qao[1] << endl;
    cout << "x- -> z+; g -> x-; g -> z+" << endl;
    double the = atan(sqrt(tan(rn) * tan(rn) + tan(sn) * tan(sn))), phi = atan2(tan(sn), tan(rn));
    cout << (J[0][0] * qao[0] + J[0][1] * qao[1]) / qan[0] << ","
         << complex<double>(cos(phi), sin(phi) * cos(the)) / (cos(phi) * cos(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / (J[0][0] * qao[0] + J[0][1] * qao[1]) << ","
         << complex<double>(cos(phi), sin(phi) * cos(the)) / (cos(phi) * cos(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / qan[0]
         << endl;
    cout << (J[1][0] * qao[0] + J[1][1] * qao[1]) / qan[1] << ","
         << complex<double>(sin(phi), -cos(phi) * cos(the)) / (sin(phi) * sin(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / (J[1][0] * qao[0] + J[1][1] * qao[1]) << ","
         << complex<double>(sin(phi), -cos(phi) * cos(the)) / (sin(phi) * sin(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / qan[1]
         << endl;
    ho11 = pow(1 - sin(ro) * sin(ro) * sin(so) * sin(so), -2);
    ho12 = -0.25 * sin(2 * ro) * sin(2 * so) * ho11;
    ho22 = cos(ro) * cos(ro) * ho11;
    ho11 = cos(so) * cos(so) * ho11;
    hn11 = pow(1 - sin(rn) * sin(rn) * sin(sn) * sin(sn), -2);
    hn12 = -0.25 * sin(2 * rn) * sin(2 * sn) * hn11;
    hn22 = cos(rn) * cos(rn) * hn11;
    hn11 = cos(sn) * cos(sn) * hn11;
    cout << ho11 << "," << ho12 << "," << ho22 << endl;
    cout << hn11 * J[0][0] * J[0][0] + hn12 * J[0][0] * J[1][0] + hn12 * J[1][0] * J[0][0] + hn22 * J[1][0] * J[1][0] << ","
         << hn11 * J[0][0] * J[0][1] + hn12 * J[0][0] * J[1][1] + hn12 * J[1][0] * J[0][1] + hn22 * J[1][0] * J[1][1] << ","
         << hn11 * J[0][1] * J[0][1] + hn12 * J[0][1] * J[1][1] + hn12 * J[1][1] * J[0][1] + hn22 * J[1][1] * J[1][1] << endl;
    cout << "swtf = " << gont << endl;
    f2 = false;
  }
  else if (f3 && tsst == 0 && ssst == 4)
  {
    cout << "tsst = " << tsst << ", ssst = " << ssst << endl;
    cout << tan(rn) << "," << tan(sn) << "," << lan[0] * qan[0] + lan[1] * qan[1] << endl;
    cout << tan(ro) << "," << tan(so) << "," << lao[0] * qao[0] + lao[1] * qao[1] << endl;
    cout << "y+ -> z+; g -> y+; g -> z+" << endl;
    double the = atan(sqrt(tan(rn) * tan(rn) + tan(sn) * tan(sn))), phi = atan2(tan(sn), tan(rn));
    cout << (J[0][0] * qao[0] + J[0][1] * qao[1]) / qan[0] << ","
         << complex<double>(cos(phi), sin(phi) * cos(the)) / (cos(phi) * cos(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / (J[0][0] * qao[0] + J[0][1] * qao[1]) << ","
         << complex<double>(cos(phi), sin(phi) * cos(the)) / (cos(phi) * cos(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / qan[0]
         << endl;
    cout << (J[1][0] * qao[0] + J[1][1] * qao[1]) / qan[1] << ","
         << complex<double>(sin(phi), -cos(phi) * cos(the)) / (sin(phi) * sin(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / (J[1][0] * qao[0] + J[1][1] * qao[1]) << ","
         << complex<double>(sin(phi), -cos(phi) * cos(the)) / (sin(phi) * sin(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / qan[1]
         << endl;
    ho12 = -0.25 * sin(2 * ro) * sin(2 * so) * ho11;
    ho22 = cos(ro) * cos(ro) * ho11;
    ho11 = cos(so) * cos(so) * ho11;
    hn11 = pow(1 - sin(rn) * sin(rn) * sin(sn) * sin(sn), -2);
    hn12 = -0.25 * sin(2 * rn) * sin(2 * sn) * hn11;
    hn22 = cos(rn) * cos(rn) * hn11;
    hn11 = cos(sn) * cos(sn) * hn11;
    cout << ho11 << "," << ho12 << "," << ho22 << endl;
    cout << hn11 * J[0][0] * J[0][0] + hn12 * J[0][0] * J[1][0] + hn12 * J[1][0] * J[0][0] + hn22 * J[1][0] * J[1][0] << ","
         << hn11 * J[0][0] * J[0][1] + hn12 * J[0][0] * J[1][1] + hn12 * J[1][0] * J[0][1] + hn22 * J[1][0] * J[1][1] << ","
         << hn11 * J[0][1] * J[0][1] + hn12 * J[0][1] * J[1][1] + hn12 * J[1][1] * J[0][1] + hn22 * J[1][1] * J[1][1] << endl;
    cout << "swtf = " << gont << endl;
    f3 = false;
  }
  else if (f4 && tsst == 0 && ssst == 5)
  {
    cout << "tsst = " << tsst << ", ssst = " << ssst << endl;
    cout << tan(rn) << "," << tan(sn) << "," << lan[0] * qan[0] + lan[1] * qan[1] << endl;
    cout << tan(ro) << "," << tan(so) << "," << lao[0] * qao[0] + lao[1] * qao[1] << endl;
    cout << "y- -> z+; g -> y-; g -> z+" << endl;
    double the = atan(sqrt(tan(rn) * tan(rn) + tan(sn) * tan(sn))), phi = atan2(tan(sn), tan(rn));
    cout << (J[0][0] * qao[0] + J[0][1] * qao[1]) / qan[0] << ","
         << complex<double>(cos(phi), sin(phi) * cos(the)) / (cos(phi) * cos(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / (J[0][0] * qao[0] + J[0][1] * qao[1]) << ","
         << complex<double>(cos(phi), sin(phi) * cos(the)) / (cos(phi) * cos(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / qan[0]
         << endl;
    cout << (J[1][0] * qao[0] + J[1][1] * qao[1]) / qan[1] << ","
         << complex<double>(sin(phi), -cos(phi) * cos(the)) / (sin(phi) * sin(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / (J[1][0] * qao[0] + J[1][1] * qao[1]) << ","
         << complex<double>(sin(phi), -cos(phi) * cos(the)) / (sin(phi) * sin(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / qan[1]
         << endl;
    ho11 = pow(1 - sin(ro) * sin(ro) * sin(so) * sin(so), -2);
    ho12 = -0.25 * sin(2 * ro) * sin(2 * so) * ho11;
    ho22 = cos(ro) * cos(ro) * ho11;
    ho11 = cos(so) * cos(so) * ho11;
    hn11 = pow(1 - sin(rn) * sin(rn) * sin(sn) * sin(sn), -2);
    hn12 = -0.25 * sin(2 * rn) * sin(2 * sn) * hn11;
    hn22 = cos(rn) * cos(rn) * hn11;
    hn11 = cos(sn) * cos(sn) * hn11;
    cout << ho11 << "," << ho12 << "," << ho22 << endl;
    cout << hn11 * J[0][0] * J[0][0] + hn12 * J[0][0] * J[1][0] + hn12 * J[1][0] * J[0][0] + hn22 * J[1][0] * J[1][0] << ","
         << hn11 * J[0][0] * J[0][1] + hn12 * J[0][0] * J[1][1] + hn12 * J[1][0] * J[0][1] + hn22 * J[1][0] * J[1][1] << ","
         << hn11 * J[0][1] * J[0][1] + hn12 * J[0][1] * J[1][1] + hn12 * J[1][1] * J[0][1] + hn22 * J[1][1] * J[1][1] << endl;
    cout << "swtf = " << gont << endl;
    f4 = false;
  }
  else if (f5 && tsst == 1 && ssst == 2)
  {
    cout << "tsst = " << tsst << ", ssst = " << ssst << endl;
    cout << tan(rn) << "," << tan(sn) << "," << lan[0] * qan[0] + lan[1] * qan[1] << endl;
    cout << tan(ro) << "," << tan(so) << "," << lao[0] * qao[0] + lao[1] * qao[1] << endl;
    cout << "x+ -> z-; g -> x+; g -> z-" << endl;
    double the = PI - atan(sqrt(tan(rn) * tan(rn) + tan(sn) * tan(sn))), phi = atan2(tan(sn), tan(rn));
    cout << (J[0][0] * qao[0] + J[0][1] * qao[1]) / qan[0] << ","
         << complex<double>(cos(phi), sin(phi) * cos(the)) / (cos(phi) * cos(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / (J[0][0] * qao[0] + J[0][1] * qao[1]) << ","
         << complex<double>(cos(phi), sin(phi) * cos(the)) / (cos(phi) * cos(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / qan[0]
         << endl;
    cout << (J[1][0] * qao[0] + J[1][1] * qao[1]) / qan[1] << ","
         << complex<double>(sin(phi), -cos(phi) * cos(the)) / (sin(phi) * sin(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / (J[1][0] * qao[0] + J[1][1] * qao[1]) << ","
         << complex<double>(sin(phi), -cos(phi) * cos(the)) / (sin(phi) * sin(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / qan[1]
         << endl;
    ho11 = pow(1 - sin(ro) * sin(ro) * sin(so) * sin(so), -2);
    ho12 = -0.25 * sin(2 * ro) * sin(2 * so) * ho11;
    ho22 = cos(ro) * cos(ro) * ho11;
    ho11 = cos(so) * cos(so) * ho11;
    hn11 = pow(1 - sin(rn) * sin(rn) * sin(sn) * sin(sn), -2);
    hn12 = -0.25 * sin(2 * rn) * sin(2 * sn) * hn11;
    hn22 = cos(rn) * cos(rn) * hn11;
    hn11 = cos(sn) * cos(sn) * hn11;
    cout << ho11 << "," << ho12 << "," << ho22 << endl;
    cout << hn11 * J[0][0] * J[0][0] + hn12 * J[0][0] * J[1][0] + hn12 * J[1][0] * J[0][0] + hn22 * J[1][0] * J[1][0] << ","
         << hn11 * J[0][0] * J[0][1] + hn12 * J[0][0] * J[1][1] + hn12 * J[1][0] * J[0][1] + hn22 * J[1][0] * J[1][1] << ","
         << hn11 * J[0][1] * J[0][1] + hn12 * J[0][1] * J[1][1] + hn12 * J[1][1] * J[0][1] + hn22 * J[1][1] * J[1][1] << endl;
    cout << "swtf = " << gont << endl;
    f5 = false;
  }
  else if (f6 && tsst == 1 && ssst == 3)
  {
    cout << "tsst = " << tsst << ", ssst = " << ssst << endl;
    cout << tan(rn) << "," << tan(sn) << "," << lan[0] * qan[0] + lan[1] * qan[1] << endl;
    cout << tan(ro) << "," << tan(so) << "," << lao[0] * qao[0] + lao[1] * qao[1] << endl;
    cout << "x- -> z-; g -> x-; g -> z-" << endl;
    double the = PI - atan(sqrt(tan(rn) * tan(rn) + tan(sn) * tan(sn))), phi = atan2(tan(sn), tan(rn));
    cout << (J[0][0] * qao[0] + J[0][1] * qao[1]) / qan[0] << ","
         << complex<double>(cos(phi), sin(phi) * cos(the)) / (cos(phi) * cos(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / (J[0][0] * qao[0] + J[0][1] * qao[1]) << ","
         << complex<double>(cos(phi), sin(phi) * cos(the)) / (cos(phi) * cos(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / qan[0]
         << endl;
    cout << (J[1][0] * qao[0] + J[1][1] * qao[1]) / qan[1] << ","
         << complex<double>(sin(phi), -cos(phi) * cos(the)) / (sin(phi) * sin(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / (J[1][0] * qao[0] + J[1][1] * qao[1]) << ","
         << complex<double>(sin(phi), -cos(phi) * cos(the)) / (sin(phi) * sin(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / qan[1]
         << endl;
    ho11 = pow(1 - sin(ro) * sin(ro) * sin(so) * sin(so), -2);
    ho12 = -0.25 * sin(2 * ro) * sin(2 * so) * ho11;
    ho22 = cos(ro) * cos(ro) * ho11;
    ho11 = cos(so) * cos(so) * ho11;
    hn11 = pow(1 - sin(rn) * sin(rn) * sin(sn) * sin(sn), -2);
    hn12 = -0.25 * sin(2 * rn) * sin(2 * sn) * hn11;
    hn22 = cos(rn) * cos(rn) * hn11;
    hn11 = cos(sn) * cos(sn) * hn11;
    cout << ho11 << "," << ho12 << "," << ho22 << endl;
    cout << hn11 * J[0][0] * J[0][0] + hn12 * J[0][0] * J[1][0] + hn12 * J[1][0] * J[0][0] + hn22 * J[1][0] * J[1][0] << ","
         << hn11 * J[0][0] * J[0][1] + hn12 * J[0][0] * J[1][1] + hn12 * J[1][0] * J[0][1] + hn22 * J[1][0] * J[1][1] << ","
         << hn11 * J[0][1] * J[0][1] + hn12 * J[0][1] * J[1][1] + hn12 * J[1][1] * J[0][1] + hn22 * J[1][1] * J[1][1] << endl;
    cout << "swtf = " << gont << endl;
    f6 = false;
  }
  else if (f7 && tsst == 1 && ssst == 4)
  {
    cout << "tsst = " << tsst << ", ssst = " << ssst << endl;
    cout << tan(rn) << "," << tan(sn) << "," << lan[0] * qan[0] + lan[1] * qan[1] << endl;
    cout << tan(ro) << "," << tan(so) << "," << lao[0] * qao[0] + lao[1] * qao[1] << endl;
    cout << "y+ -> z-; g -> y+; g -> z-" << endl;
    double the = PI - atan(sqrt(tan(rn) * tan(rn) + tan(sn) * tan(sn))), phi = atan2(tan(sn), tan(rn));
    cout << (J[0][0] * qao[0] + J[0][1] * qao[1]) / qan[0] << ","
         << complex<double>(cos(phi), sin(phi) * cos(the)) / (cos(phi) * cos(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / (J[0][0] * qao[0] + J[0][1] * qao[1]) << ","
         << complex<double>(cos(phi), sin(phi) * cos(the)) / (cos(phi) * cos(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / qan[0]
         << endl;
    cout << (J[1][0] * qao[0] + J[1][1] * qao[1]) / qan[1] << ","
         << complex<double>(sin(phi), -cos(phi) * cos(the)) / (sin(phi) * sin(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / (J[1][0] * qao[0] + J[1][1] * qao[1]) << ","
         << complex<double>(sin(phi), -cos(phi) * cos(the)) / (sin(phi) * sin(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / qan[1]
         << endl;
    ho11 = pow(1 - sin(ro) * sin(ro) * sin(so) * sin(so), -2);
    ho12 = -0.25 * sin(2 * ro) * sin(2 * so) * ho11;
    ho22 = cos(ro) * cos(ro) * ho11;
    ho11 = cos(so) * cos(so) * ho11;
    hn11 = pow(1 - sin(rn) * sin(rn) * sin(sn) * sin(sn), -2);
    hn12 = -0.25 * sin(2 * rn) * sin(2 * sn) * hn11;
    hn22 = cos(rn) * cos(rn) * hn11;
    hn11 = cos(sn) * cos(sn) * hn11;
    cout << ho11 << "," << ho12 << "," << ho22 << endl;
    cout << hn11 * J[0][0] * J[0][0] + hn12 * J[0][0] * J[1][0] + hn12 * J[1][0] * J[0][0] + hn22 * J[1][0] * J[1][0] << ","
         << hn11 * J[0][0] * J[0][1] + hn12 * J[0][0] * J[1][1] + hn12 * J[1][0] * J[0][1] + hn22 * J[1][0] * J[1][1] << ","
         << hn11 * J[0][1] * J[0][1] + hn12 * J[0][1] * J[1][1] + hn12 * J[1][1] * J[0][1] + hn22 * J[1][1] * J[1][1] << endl;
    cout << "swtf = " << gont << endl;
    f7 = false;
  }
  else if (f8 && tsst == 1 && ssst == 5)
  {
    cout << "tsst = " << tsst << ", ssst = " << ssst << endl;
    cout << tan(rn) << "," << tan(sn) << "," << lan[0] * qan[0] + lan[1] * qan[1] << endl;
    cout << tan(ro) << "," << tan(so) << "," << lao[0] * qao[0] + lao[1] * qao[1] << endl;
    cout << "y- -> z-; g -> y-; g -> z-" << endl;
    double the = PI - atan(sqrt(tan(rn) * tan(rn) + tan(sn) * tan(sn))), phi = atan2(tan(sn), tan(rn));
    cout << (J[0][0] * qao[0] + J[0][1] * qao[1]) / qan[0] << ","
         << complex<double>(cos(phi), sin(phi) * cos(the)) / (cos(phi) * cos(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / (J[0][0] * qao[0] + J[0][1] * qao[1]) << ","
         << complex<double>(cos(phi), sin(phi) * cos(the)) / (cos(phi) * cos(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / qan[0]
         << endl;
    cout << (J[1][0] * qao[0] + J[1][1] * qao[1]) / qan[1] << ","
         << complex<double>(sin(phi), -cos(phi) * cos(the)) / (sin(phi) * sin(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / (J[1][0] * qao[0] + J[1][1] * qao[1]) << ","
         << complex<double>(sin(phi), -cos(phi) * cos(the)) / (sin(phi) * sin(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / qan[1]
         << endl;
    ho11 = pow(1 - sin(ro) * sin(ro) * sin(so) * sin(so), -2);
    ho12 = -0.25 * sin(2 * ro) * sin(2 * so) * ho11;
    ho22 = cos(ro) * cos(ro) * ho11;
    ho11 = cos(so) * cos(so) * ho11;
    hn11 = pow(1 - sin(rn) * sin(rn) * sin(sn) * sin(sn), -2);
    hn12 = -0.25 * sin(2 * rn) * sin(2 * sn) * hn11;
    hn22 = cos(rn) * cos(rn) * hn11;
    hn11 = cos(sn) * cos(sn) * hn11;
    cout << ho11 << "," << ho12 << "," << ho22 << endl;
    cout << hn11 * J[0][0] * J[0][0] + hn12 * J[0][0] * J[1][0] + hn12 * J[1][0] * J[0][0] + hn22 * J[1][0] * J[1][0] << ","
         << hn11 * J[0][0] * J[0][1] + hn12 * J[0][0] * J[1][1] + hn12 * J[1][0] * J[0][1] + hn22 * J[1][0] * J[1][1] << ","
         << hn11 * J[0][1] * J[0][1] + hn12 * J[0][1] * J[1][1] + hn12 * J[1][1] * J[0][1] + hn22 * J[1][1] * J[1][1] << endl;
    cout << "swtf = " << gont << endl;
    f8 = false;
  }
  else if (f9 && tsst == 2 && ssst == 0)
  {
    cout << "tsst = " << tsst << ", ssst = " << ssst << endl;
    cout << tan(rn) << "," << tan(sn) << "," << lan[0] * qan[0] + lan[1] * qan[1] << endl;
    cout << tan(ro) << "," << tan(so) << "," << lao[0] * qao[0] + lao[1] * qao[1] << endl;
    cout << "z+ -> x+; g -> z+; g -> x+" << endl;
    double the = atan(sqrt(tan(rn) * tan(rn) + 1) / tan(sn)), phi = rn;
    if (the < 0)
      the = PI + the;
    cout << (J[0][0] * qao[0] + J[0][1] * qao[1]) / qan[0] << ","
         << complex<double>(0, -1) / sin(the) / (J[0][0] * qao[0] + J[0][1] * qao[1]) << ","
         << complex<double>(0, -1) / sin(the) / qan[0]
         << endl;
    cout << (J[1][0] * qao[0] + J[1][1] * qao[1]) / qan[1] << ","
         << complex<double>(-cos(phi), -sin(phi) * cos(the)) / (cos(phi) * cos(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / (J[1][0] * qao[0] + J[1][1] * qao[1]) << ","
         << complex<double>(-cos(phi), -sin(phi) * cos(the)) / (cos(phi) * cos(phi) * sin(the) * sin(the) + cos(the) * cos(the)) / qan[1]
         << endl;
    ho11 = pow(1 - sin(ro) * sin(ro) * sin(so) * sin(so), -2);
    ho12 = -0.25 * sin(2 * ro) * sin(2 * so) * ho11;
    ho22 = cos(ro) * cos(ro) * ho11;
    ho11 = cos(so) * cos(so) * ho11;
    hn11 = pow(1 - sin(rn) * sin(rn) * sin(sn) * sin(sn), -2);
    hn12 = -0.25 * sin(2 * rn) * sin(2 * sn) * hn11;
    hn22 = cos(rn) * cos(rn) * hn11;
    hn11 = cos(sn) * cos(sn) * hn11;
    cout << ho11 << "," << ho12 << "," << ho22 << endl;
    cout << hn11 * J[0][0] * J[0][0] + hn12 * J[0][0] * J[1][0] + hn12 * J[1][0] * J[0][0] + hn22 * J[1][0] * J[1][0] << ","
         << hn11 * J[0][0] * J[0][1] + hn12 * J[0][0] * J[1][1] + hn12 * J[1][0] * J[0][1] + hn22 * J[1][0] * J[1][1] << ","
         << hn11 * J[0][1] * J[0][1] + hn12 * J[0][1] * J[1][1] + hn12 * J[1][1] * J[0][1] + hn22 * J[1][1] * J[1][1] << endl;
    cout << "swtf = " << gont << endl;
    f9 = false;
  }
  else if (f10 && tsst == 2 && ssst == 1)
  {
    cout << "tsst = " << tsst << ", ssst = " << ssst << endl;
    cout << tan(rn) << "," << tan(sn) << "," << lan[0] * qan[0] + lan[1] * qan[1] << endl;
    cout << tan(ro) << "," << tan(so) << "," << lao[0] * qao[0] + lao[1] * qao[1] << endl;
    cout << (J[0][0] * qao[0] + J[0][1] * qao[1]) / qan[0] << endl;
    cout << (J[1][0] * qao[0] + J[1][1] * qao[1]) / qan[1] << endl;
    ho11 = pow(1 - sin(ro) * sin(ro) * sin(so) * sin(so), -2);
    ho12 = -0.25 * sin(2 * ro) * sin(2 * so) * ho11;
    ho22 = cos(ro) * cos(ro) * ho11;
    ho11 = cos(so) * cos(so) * ho11;
    hn11 = pow(1 - sin(rn) * sin(rn) * sin(sn) * sin(sn), -2);
    hn12 = -0.25 * sin(2 * rn) * sin(2 * sn) * hn11;
    hn22 = cos(rn) * cos(rn) * hn11;
    hn11 = cos(sn) * cos(sn) * hn11;
    cout << ho11 << "," << ho12 << "," << ho22 << endl;
    cout << hn11 * J[0][0] * J[0][0] + hn12 * J[0][0] * J[1][0] + hn12 * J[1][0] * J[0][0] + hn22 * J[1][0] * J[1][0] << ","
         << hn11 * J[0][0] * J[0][1] + hn12 * J[0][0] * J[1][1] + hn12 * J[1][0] * J[0][1] + hn22 * J[1][0] * J[1][1] << ","
         << hn11 * J[0][1] * J[0][1] + hn12 * J[0][1] * J[1][1] + hn12 * J[1][1] * J[0][1] + hn22 * J[1][1] * J[1][1] << endl;
    cout << "swtf = " << gont << endl;
    f10 = false;
  }
  else if (f11 && tsst == 2 && ssst == 4)
  {
    cout << "tsst = " << tsst << ", ssst = " << ssst << endl;
    cout << tan(rn) << "," << tan(sn) << "," << lan[0] * qan[0] + lan[1] * qan[1] << endl;
    cout << tan(ro) << "," << tan(so) << "," << lao[0] * qao[0] + lao[1] * qao[1] << endl;
    cout << (J[0][0] * qao[0] + J[0][1] * qao[1]) / qan[0] << endl;
    cout << (J[1][0] * qao[0] + J[1][1] * qao[1]) / qan[1] << endl;
    ho11 = pow(1 - sin(ro) * sin(ro) * sin(so) * sin(so), -2);
    ho12 = -0.25 * sin(2 * ro) * sin(2 * so) * ho11;
    ho22 = cos(ro) * cos(ro) * ho11;
    ho11 = cos(so) * cos(so) * ho11;
    hn11 = pow(1 - sin(rn) * sin(rn) * sin(sn) * sin(sn), -2);
    hn12 = -0.25 * sin(2 * rn) * sin(2 * sn) * hn11;
    hn22 = cos(rn) * cos(rn) * hn11;
    hn11 = cos(sn) * cos(sn) * hn11;
    cout << ho11 << "," << ho12 << "," << ho22 << endl;
    cout << hn11 * J[0][0] * J[0][0] + hn12 * J[0][0] * J[1][0] + hn12 * J[1][0] * J[0][0] + hn22 * J[1][0] * J[1][0] << ","
         << hn11 * J[0][0] * J[0][1] + hn12 * J[0][0] * J[1][1] + hn12 * J[1][0] * J[0][1] + hn22 * J[1][0] * J[1][1] << ","
         << hn11 * J[0][1] * J[0][1] + hn12 * J[0][1] * J[1][1] + hn12 * J[1][1] * J[0][1] + hn22 * J[1][1] * J[1][1] << endl;
    cout << "swtf = " << gont << endl;
    f11 = false;
  }
  else if (f12 && tsst == 2 && ssst == 5)
  {
    cout << "tsst = " << tsst << ", ssst = " << ssst << endl;
    cout << tan(rn) << "," << tan(sn) << "," << lan[0] * qan[0] + lan[1] * qan[1] << endl;
    cout << tan(ro) << "," << tan(so) << "," << lao[0] * qao[0] + lao[1] * qao[1] << endl;
    cout << (J[0][0] * qao[0] + J[0][1] * qao[1]) / qan[0] << endl;
    cout << (J[1][0] * qao[0] + J[1][1] * qao[1]) / qan[1] << endl;
    ho11 = pow(1 - sin(ro) * sin(ro) * sin(so) * sin(so), -2);
    ho12 = -0.25 * sin(2 * ro) * sin(2 * so) * ho11;
    ho22 = cos(ro) * cos(ro) * ho11;
    ho11 = cos(so) * cos(so) * ho11;
    hn11 = pow(1 - sin(rn) * sin(rn) * sin(sn) * sin(sn), -2);
    hn12 = -0.25 * sin(2 * rn) * sin(2 * sn) * hn11;
    hn22 = cos(rn) * cos(rn) * hn11;
    hn11 = cos(sn) * cos(sn) * hn11;
    cout << ho11 << "," << ho12 << "," << ho22 << endl;
    cout << hn11 * J[0][0] * J[0][0] + hn12 * J[0][0] * J[1][0] + hn12 * J[1][0] * J[0][0] + hn22 * J[1][0] * J[1][0] << ","
         << hn11 * J[0][0] * J[0][1] + hn12 * J[0][0] * J[1][1] + hn12 * J[1][0] * J[0][1] + hn22 * J[1][0] * J[1][1] << ","
         << hn11 * J[0][1] * J[0][1] + hn12 * J[0][1] * J[1][1] + hn12 * J[1][1] * J[0][1] + hn22 * J[1][1] * J[1][1] << endl;
    cout << "swtf = " << gont << endl;
    f12 = false;
  }

#endif

  return gont;
}
#endif
// for check
// used by _dst construction, so these x,y,z must coinside with grid point
// we have considered ghost points now
void NullShellPatch::prolongpointstru(MyList<pointstru> *&psul, MyList<ss_patch> *sPpi, double DH[dim],
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
              ps->data->swtf = get_swtf(ps->data->lpox, ps->data->tsst, ps->data->ssst);
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
      cout << "somthing is wrong in NullShellPatch::prolongpointstru" << endl;
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
            ps->data->swtf = 1;
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
      cout << "NullShellPatch::prolongpointstru fail to find target Block for pointstru:" << endl;
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
      ps->data->swtf = pss->data->swtf;
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
bool NullShellPatch::prolongpointstru(MyList<pointstru> *&psul, bool ssyn, int tsst, MyList<ss_patch> *sPp, double DH[dim],
                                      MyList<Patch> *Pp, double CDH[dim], double x, double y, double z, int Symmetry, int rank_in)
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
              ps->data->swtf = 1;
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
            ps->data->swtf = 1;
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
bool NullShellPatch::prolongpointstru_ss(MyList<pointstru> *&psul, int tsst, MyList<ss_patch> *sPp, double DH[dim],
                                         MyList<Patch> *Pp, double CDH[dim], double x, double y, double z, int Symmetry, int rank_in)
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
            ps->data->swtf = 1;
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
// setup interpatch interpolation stuffs
void NullShellPatch::setupintintstuff(int cpusize, MyList<Patch> *CPatL, int Symmetry)
{
  const int hCS_width = 0; // do not input data from null shell to box
  const int hSC_width = 1; // do     input data from box to null shell
  int myrank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (myrank == 0)
    cout << "NullShellPatch::setupintintstuff begines..." << endl;

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
              flag = prolongpointstru(cs_src[i], false, sPp->data->sst, PatL, DH, CPatL, CDH, gx, gy, gz, Symmetry, i);
              if (flag)
                break;
            }
            if (!flag)
            {
              CPatL->data->checkBlock();
              if (myrank == 0)
              {
                cout << "ShellPatch::prolongpointstru fail to find cardisian source point for" << endl;
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
                  flag = prolongpointstru(ss_src[i], true, sPp->data->sst, PatL, DH, CPatL, CDH, gx, gy, gz, Symmetry, i);
                else
                  flag = prolongpointstru_ss(ss_src[i], sPp->data->sst, PatL, DH, CPatL, CDH, x, y, z, Symmetry, i);
                if (flag)
                  break;
              }
              if (!flag)
              {
                if (myrank == 0)
                {
                  // if you used Vertex grid please note x=1, try 0.999999 instead
                  cout << "NullShellPatch::prolongpointstru fail to find shell source point for" << endl;
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
    cout << "NullShellPatch::setupintintstuff ss_src completes" << endl;

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
              flag = prolongpointstru(cs_src[i], true, -1, PatL, DH, CPatL, CDH, x, y, z, Symmetry, i);
              if (flag)
                break;
            }
            if (!flag)
            {
              if (myrank == 0)
              {
                cout << "ShellPatch::prolongpointstru fail to find shell source point for" << endl;
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
      cout << "NullShellPatch::setupintintstuff cs_src completes" << endl;
    else
      cout << "NullShellPatch::no cs_src exists" << endl;

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
    cout << "NullShellPatch::setupintintstuff ss_dst and cs_dst complete" << endl;

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
void NullShellPatch::checkPatch()
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (myrank == 0)
  {
    cout << " belong to NullShell Patchs " << endl;
    MyList<ss_patch> *Pp = PatL;
    while (Pp)
    {
      cout << " shape: [";
      for (int i = 0; i < dim; i++)
      {
        cout << Pp->data->shape[i];
        if (i < dim - 1)
          cout << ",";
        else
          cout << "]" << endl;
      }
      cout << " range:" << "(";
      for (int i = 0; i < dim; i++)
      {
        cout << Pp->data->bbox[i] << ":" << Pp->data->bbox[dim + i];
        if (i < dim - 1)
          cout << ",";
        else
          cout << ")" << endl;
      }
      Pp = Pp->next;
    }
  }
}
void NullShellPatch::checkBlock(int sst)
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
double NullShellPatch::getdX(int dir)
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
void NullShellPatch::shellname(char *sn, int i)
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
// Now we dump the data including overlap points
void NullShellPatch::Dump_xyz(char *tag, double time, double dT)
{
  MyList<var> *DumpListi = 0;
  DumpListi = new MyList<var>(gx);
  DumpListi->insert(gy);
  DumpListi->insert(gz);
  Dump_Data(DumpListi, tag, time, dT);
  DumpListi->clearList();
}
void NullShellPatch::Dump_Data(MyList<var> *DumpListi, char *tag, double time, double dT)
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
void NullShellPatch::intertransfer(MyList<pointstru> **src, MyList<pointstru> **dst,
                                   MyList<var> *VarList1 /* source */, MyList<var> *VarList2 /*target */,
                                   int Symmetry, int *Varwt)
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
      if (length = interdata_packer(0, src[myrank], dst[myrank], node, PACK, VarList1, VarList2, Symmetry, Varwt))
      {
        rec_data[node] = new double[length];
        if (!rec_data[node])
        {
          cout << "out of memory when new in short transfer, place 1" << endl;
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
        interdata_packer(rec_data[node], src[myrank], dst[myrank], node, PACK, VarList1, VarList2, Symmetry, Varwt);
      }
    }
    else
    {
      // send from this cpu to cpu#node
      if (length = interdata_packer(0, src[myrank], dst[myrank], node, PACK, VarList1, VarList2, Symmetry, Varwt))
      {
        send_data[node] = new double[length];
        if (!send_data[node])
        {
          cout << "out of memory when new in short transfer, place 2" << endl;
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
        interdata_packer(send_data[node], src[myrank], dst[myrank], node, PACK, VarList1, VarList2, Symmetry, Varwt);
        MPI_Isend((void *)send_data[node], length, MPI_DOUBLE, node, 1, MPI_COMM_WORLD, reqs + req_no++);
      }
      // receive from cpu#node to this cpu
      if (length = interdata_packer(0, src[node], dst[node], node, UNPACK, VarList1, VarList2, Symmetry, Varwt))
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
      interdata_packer(rec_data[node], src[node], dst[node], node, UNPACK, VarList1, VarList2, Symmetry, Varwt);

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
int NullShellPatch::interdata_packer(double *data, MyList<pointstru> *src, MyList<pointstru> *dst, int rank_in, int dir,
                                     MyList<var> *VarLists /* source */, MyList<var> *VarListd /* target */, int Symmetry, int *Varwt)
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
      int vind = 0;
      bool flag = true;
      while (varls && varld)
      {
        if (data)
        {
          if (dir == PACK)
          {
            /*
                     f_global_interp(src->data->Bg->shape,src->data->Bg->X[0],src->data->Bg->X[1],src->data->Bg->X[2],
                   src->data->Bg->fgfs[varls->data->sgfn],data[size_out],
                   src->data->lpox[0],src->data->lpox[1],src->data->lpox[2],ordn,varls->data->SoA,Symmetry);
            */
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
              cout << "NullShellPatch::interdata_packer: not recognized DIM = " << DIMh << endl;
              MPI_Abort(MPI_COMM_WORLD, 1);
            }
          }
          if (dir == UNPACK) // from target data to corresponding grid
          {
            if (Varwt[vind / 2] != 0) // we always assume 2 time number relation
            {
              if (flag)
              {
                complex<double> rtp = complex<double>(data[size_out], data[size_out + 1]);
                rtp = rtp * pow(dst->data->swtf, Varwt[vind / 2]); // note we only stored the factor in dst
                data[size_out] = rtp.real();
                data[size_out + 1] = rtp.imag();
              }
              flag = !flag; // on-off method
            }
            // if(dst->data->tsst==2 && fabs(dst->data->lpox[0]+0.02617993878)<0.00001 && fabs(dst->data->lpox[2]-0.510417)<0.00001)cout<<varld->data->name<<endl;
            f_pointcopy(DIM, dst->data->Bg->bbox, dst->data->Bg->bbox + dim, dst->data->Bg->shape, dst->data->Bg->fgfs[varld->data->sgfn],
                        dst->data->lpox[0], dst->data->lpox[1], dst->data->lpox[2], data[size_out]);
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
void NullShellPatch::Synch(MyList<var> *VarList, int Symmetry, int *Varwt)
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

  intertransfer(ss_src, ss_dst, VarList, VarList, Symmetry, Varwt);

  // we need this here to correct conners
  if (Symmetry > 0)
    fill_symmetric_boundarybuffer(VarList, Varwt);
}
void NullShellPatch::check_pointstrul(MyList<pointstru> *pp, bool first_only)
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
void NullShellPatch::check_pointstrul2(MyList<pointstru> *pp, int first_last_only)
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
      if (first_last_only == 2)
      {
        if (pp->next == 0)
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
        }
      }
      else
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
        if (first_last_only == 1)
          return;
      }
      pp = pp->next;
    }
  }
}
void NullShellPatch::matchcheck(MyList<Patch> *CPatL)
{
  double cbd = CPatL->data->bbox[dim];
  for (int i = 1; i < dim; i++)
    cbd = Mymin(cbd, CPatL->data->bbox[dim + i]);
  cbd = cbd - xmin * Rmin / (1 - xmin);
  double dr, dc;
  dc = CPatL->data->getdX(0);
  dr = getdX(2);
  for (int i = 1; i < dim; i++)
  {
    dc = Mymax(dc, CPatL->data->getdX(i));
    //    dr = Mymax(dr,getdX(i));
  }

  int ir, ic;
  ir = int(cbd / dr);
  ic = int(cbd / dc);
  if (Mymin(ir, ic) < 3 * ghost_width)
  {
    int myrank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank == 0)
    {
      cout << "NullShell Patches insert too shallow:" << endl;
      cout << "distantance between these two boundaries is " << cbd << ", spatial step is " << Mymax(dc, dr) << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }
}
void NullShellPatch::Interp_Points(MyList<var> *VarList,
                                   int NN, double **XX, /*input global Cartesian coordinate*/
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
    getlocalpox(XX[0][j], XX[1][j], XX[2][j], sst, pox[0], pox[1], pox[2]); // pox[2] is x indeed

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
            f_global_interp_ss(BP->shape, BP->X[0], BP->X[1], BP->X[2], BP->fgfs[varl->data->sgfn], shellf[j * num_var + k],
                               pox[0], pox[1], pox[2], ordn, varl->data->SoA, Symmetry, sst);
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
void NullShellPatch::Interp_Points_2D(MyList<var> *VarList,
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

    int indZ = int((pox[2] - xmin) / DH[2]);
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
void NullShellPatch::Step(double dT, double PhysTime, monitor *ErrorMonitor)
{
  int iter_count = 0; // count RK4 substeps
  int pre = 0, cor = 1;
  int ERROR = 0;
  double TT = PhysTime;
  double neps = 0.05;
  MyList<ss_patch> *sPp;

  // Predictor
  HyperSlice(dT, TT, ErrorMonitor, iter_count);
  {
    sPp = PatL;
    while (sPp)
    {
      MyList<Block> *BP = sPp->data->blb;
      while (BP)
      {
        Block *cg = BP->data;
        //          cg->swapList(TheList,JrhsList,myrank);
        if (myrank == cg->rank)
        {
          // rhs calculation
          f_array_copy(cg->shape, cg->fgfs[RJ_rhs->sgfn], cg->fgfs[RTheta->sgfn]);
          f_array_copy(cg->shape, cg->fgfs[IJ_rhs->sgfn], cg->fgfs[ITheta->sgfn]);
          f_kodis_shor(cg->shape, cg->X[0], cg->X[1], cg->X[2], cg->fgfs[RJ0->sgfn], cg->fgfs[RJ_rhs->sgfn],
                       RJ0->SoA, Symmetry, neps, sPp->data->sst);
          f_kodis_shor(cg->shape, cg->X[0], cg->X[1], cg->X[2], cg->fgfs[IJ0->sgfn], cg->fgfs[IJ_rhs->sgfn],
                       RJ0->SoA, Symmetry, neps, sPp->data->sst);
          f_omega_rhs(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                      cg->fgfs[omega0->sgfn], cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn], cg->fgfs[omega_rhs->sgfn],
                      cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                      cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn]);
          f_rungekutta4_rout(cg->shape, dT, cg->fgfs[RJ0->sgfn], cg->fgfs[RJ->sgfn], cg->fgfs[RJ_rhs->sgfn],
                             iter_count);
          f_rungekutta4_rout(cg->shape, dT, cg->fgfs[IJ0->sgfn], cg->fgfs[IJ->sgfn], cg->fgfs[IJ_rhs->sgfn],
                             iter_count);
          f_rungekutta4_rout(cg->shape, dT, cg->fgfs[omega0->sgfn], cg->fgfs[omega->sgfn], cg->fgfs[omega_rhs->sgfn],
                             iter_count);
        }
        if (BP == sPp->data->ble)
          break;
        BP = BP->next;
      }
      sPp = sPp->next;
    }
  }
  /*
      {
        char str[50];
        sprintf(str,"rk%d",iter_count);
        Dump_Data(SynchList_pre,str,PhysTime,dT);
        Dump_Data(RHSList,str,PhysTime,dT);
      }
  */
  // no nedd to synchronize J, because Theta has already been synchnized previously
  int Varwt[1];
  MyList<var> *DG_List;
  DG_List = new MyList<var>(omega);
  DG_List->insert(FXZEO);
  Varwt[0] = 0;
  Synch(DG_List, Symmetry, Varwt);
  DG_List->clearList();

  Compute_News(PhysTime, dT, false); // put here because after step J and omega are at t+dt, while other variables at t

  // corrector
  for (iter_count = 1; iter_count < 4; iter_count++)
  {

    // for RK4: t0, t0+dt/2, t0+dt/2, t0+dt;
    if (iter_count == 1 || iter_count == 3)
      TT += dT / 2;
    HyperSlice(dT, TT, ErrorMonitor, iter_count);
    {
      sPp = PatL;
      while (sPp)
      {
        MyList<Block> *BP = sPp->data->blb;
        while (BP)
        {
          Block *cg = BP->data;
          //          cg->swapList(TheList,J1List,myrank);
          if (myrank == cg->rank)
          {
            // rhs calculation
            f_array_copy(cg->shape, cg->fgfs[RJ1->sgfn], cg->fgfs[RTheta->sgfn]);
            f_array_copy(cg->shape, cg->fgfs[IJ1->sgfn], cg->fgfs[ITheta->sgfn]);
            f_kodis_shor(cg->shape, cg->X[0], cg->X[1], cg->X[2], cg->fgfs[RJ0->sgfn], cg->fgfs[RJ1->sgfn],
                         RJ0->SoA, Symmetry, neps, sPp->data->sst);
            f_kodis_shor(cg->shape, cg->X[0], cg->X[1], cg->X[2], cg->fgfs[IJ0->sgfn], cg->fgfs[IJ1->sgfn],
                         RJ0->SoA, Symmetry, neps, sPp->data->sst);
            f_omega_rhs(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                        cg->fgfs[omega->sgfn], cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn], cg->fgfs[omega1->sgfn],
                        cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                        cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn]);
            f_rungekutta4_rout(cg->shape, dT, cg->fgfs[RJ0->sgfn], cg->fgfs[RJ1->sgfn], cg->fgfs[RJ_rhs->sgfn],
                               iter_count);
            f_rungekutta4_rout(cg->shape, dT, cg->fgfs[IJ0->sgfn], cg->fgfs[IJ1->sgfn], cg->fgfs[IJ_rhs->sgfn],
                               iter_count);
            f_rungekutta4_rout(cg->shape, dT, cg->fgfs[omega0->sgfn], cg->fgfs[omega1->sgfn], cg->fgfs[omega_rhs->sgfn],
                               iter_count);
          }
          if (iter_count < 3)
            cg->swapList(SynchList_cor, SynchList_pre, myrank);
          else
          {
            cg->swapList(StateList, SynchList_cor, myrank);
            cg->swapList(OldStateList, SynchList_cor, myrank);
          }
          if (BP == sPp->data->ble)
            break;
          BP = BP->next;
        }
        sPp = sPp->next;
      }
    }

    int Varwt[1];
    MyList<var> *DG_List;
    DG_List = new MyList<var>(omega0);
    DG_List->insert(FXZEO);
    Varwt[0] = 0;
    Synch(DG_List, Symmetry, Varwt);
    DG_List->clearList();

    /*
      {
          char str[50];
          sprintf(str,"rk%d",iter_count);
          Dump_Data(SynchList_cor,str,PhysTime,dT);
      }
    */
  }
}
void NullShellPatch::Null_Boundary(double PhysTime)
{
  MyList<ss_patch> *sPp;

  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        f_get_null_boundary(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                            cg->fgfs[beta->sgfn], cg->fgfs[RQ->sgfn], cg->fgfs[IQ->sgfn],
                            cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn],
                            cg->fgfs[W->sgfn], cg->fgfs[RTheta->sgfn], cg->fgfs[ITheta->sgfn],
                            cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                            cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                            cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn],
                            cg->fgfs[dquR1->sgfn], cg->fgfs[dquR2->sgfn], cg->fgfs[dquI1->sgfn], cg->fgfs[dquI2->sgfn],
                            cg->fgfs[bdquR1->sgfn], cg->fgfs[bdquR2->sgfn], cg->fgfs[bdquI1->sgfn], cg->fgfs[bdquI2->sgfn],
                            cg->fgfs[dgR->sgfn], cg->fgfs[dgI->sgfn], cg->fgfs[bdgR->sgfn], cg->fgfs[bdgI->sgfn],
                            PhysTime, Rmin, sPp->data->sst);
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }

  int Varwt[3];
  MyList<var> *DG_List;
  DG_List = new MyList<var>(RU);
  DG_List->insert(IU);
  Varwt[0] = 1;
  DG_List->insert(RQ);
  DG_List->insert(IQ);
  Varwt[1] = 1;
  DG_List->insert(RTheta);
  DG_List->insert(ITheta);
  Varwt[2] = 2;

  Synch(DG_List, Symmetry, Varwt);
  //     Dump_Data(DG_List,0,0,1);
  DG_List->clearList();
}
#if 1
// real evolve
void NullShellPatch::HyperSlice(double dT, double PhysTime, monitor *ErrorMonitor, int RK_count)
{
  int ERROR = 0;
  Null_Boundary(PhysTime);

  int spin, e;

  MyList<ss_patch> *sPp;

  // evolve beta
  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        if (RK_count == 0)
        {
          f_calculate_K(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                        cg->fgfs[RJ0->sgfn], cg->fgfs[IJ0->sgfn], cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn], cg->fgfs[KKx->sgfn], cg->fgfs[HKKx->sgfn]);
          if (f_NullEvol_beta(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                              cg->fgfs[RJ0->sgfn], cg->fgfs[IJ0->sgfn], cg->fgfs[beta->sgfn], cg->fgfs[KKx->sgfn], cg->fgfs[HKKx->sgfn]))
          {
            cout << "find NaN in NullShell domain: sst = " << sPp->data->sst << ", (" << cg->bbox[0] << ":" << cg->bbox[3] << ","
                 << cg->bbox[1] << ":" << cg->bbox[4] << "," << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
            ERROR = 1;
          }
        }
        else
        {
          f_calculate_K(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                        cg->fgfs[RJ->sgfn], cg->fgfs[IJ->sgfn], cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn], cg->fgfs[KKx->sgfn], cg->fgfs[HKKx->sgfn]);
          if (f_NullEvol_beta(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                              cg->fgfs[RJ->sgfn], cg->fgfs[IJ->sgfn], cg->fgfs[beta->sgfn], cg->fgfs[KKx->sgfn], cg->fgfs[HKKx->sgfn]))
          {
            cout << "find NaN in NullShell domain: sst = " << sPp->data->sst << ", (" << cg->bbox[0] << ":" << cg->bbox[3] << ","
                 << cg->bbox[1] << ":" << cg->bbox[4] << "," << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
            ERROR = 1;
          }
        }
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }
  // check error information
  {
    int erh = ERROR;
    MPI_Allreduce(&erh, &ERROR, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  }
  if (ERROR)
  {
    Dump_Data(betaList, 0, PhysTime, dT);
    if (myrank == 0)
    {
      if (ErrorMonitor && ErrorMonitor->outfile)
        ErrorMonitor->outfile << "find NaN in beta on NullShell Patches at t = " << PhysTime << endl;
      else
        cout << "find NaN in beta on NullShell Patches at t = " << PhysTime << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  Synch(betaList, Symmetry, betawt);
  // get nu, k and B
  spin = 2;
  e = -1;
  if (RK_count == 0)
    eth_derivs(RJ0, IJ0, Rnu, Inu, spin, e);
  else
    eth_derivs(RJ, IJ, Rnu, Inu, spin, e);
  spin = 0;
  e = 1;
  eth_derivs(KK, FXZEO, Rk, Ik, spin, e);
  eth_derivs(beta, FXZEO, RB, IB, spin, e);

  // evolve Q and U
  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        if (RK_count == 0)
        {
          if (f_NullEvol_Q(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                           cg->fgfs[RJ0->sgfn], cg->fgfs[IJ0->sgfn],
                           cg->fgfs[Rk->sgfn], cg->fgfs[Ik->sgfn],
                           cg->fgfs[Rnu->sgfn], cg->fgfs[Inu->sgfn],
                           cg->fgfs[RB->sgfn], cg->fgfs[IB->sgfn],
                           cg->fgfs[RQ->sgfn], cg->fgfs[IQ->sgfn],
                           cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn], cg->fgfs[KKx->sgfn], cg->fgfs[HKKx->sgfn],
                           cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                           cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                           cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn]) ||
              // since we do not need derivetive of Q, we can deal with U together here
              // at this stage Q has been updated already
              f_NullEvol_U(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                           cg->fgfs[RJ0->sgfn], cg->fgfs[IJ0->sgfn],
                           cg->fgfs[RQ->sgfn], cg->fgfs[IQ->sgfn], cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn],
                           cg->fgfs[beta->sgfn], cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn], Rmin))
          {
            cout << "find NaN in NullShell domain: sst = " << sPp->data->sst << ", (" << cg->bbox[0] << ":" << cg->bbox[3] << ","
                 << cg->bbox[1] << ":" << cg->bbox[4] << "," << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
            ERROR = 1;
          }
        }
        else
        {
          if (f_NullEvol_Q(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                           cg->fgfs[RJ->sgfn], cg->fgfs[IJ->sgfn],
                           cg->fgfs[Rk->sgfn], cg->fgfs[Ik->sgfn],
                           cg->fgfs[Rnu->sgfn], cg->fgfs[Inu->sgfn],
                           cg->fgfs[RB->sgfn], cg->fgfs[IB->sgfn],
                           cg->fgfs[RQ->sgfn], cg->fgfs[IQ->sgfn],
                           cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn], cg->fgfs[KKx->sgfn], cg->fgfs[HKKx->sgfn],
                           cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                           cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                           cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn]) ||
              // since we do not need derivetive of Q, we can deal with U together here
              // at this stage Q has been updated already
              f_NullEvol_U(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                           cg->fgfs[RJ->sgfn], cg->fgfs[IJ->sgfn],
                           cg->fgfs[RQ->sgfn], cg->fgfs[IQ->sgfn], cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn],
                           cg->fgfs[beta->sgfn], cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn], Rmin))
          {
            cout << "find NaN in NullShell domain: sst = " << sPp->data->sst << ", (" << cg->bbox[0] << ":" << cg->bbox[3] << ","
                 << cg->bbox[1] << ":" << cg->bbox[4] << "," << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
            ERROR = 1;
          }
        }
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }
  // check error information
  {
    int erh = ERROR;
    MPI_Allreduce(&erh, &ERROR, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  }
  if (ERROR)
  {
    Dump_Data(QUList, 0, PhysTime, dT);
    if (myrank == 0)
    {
      if (ErrorMonitor && ErrorMonitor->outfile)
        ErrorMonitor->outfile << "find NaN in Q and/or U on NullShell Patches at t = " << PhysTime << endl;
      else
        cout << "find NaN in Q and/or U on NullShell Patches at t = " << PhysTime << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  Synch(QUList, Symmetry, QUwt);

  // evolve W and Theta
  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        if (RK_count == 0)
        {
          if (f_NullEvol_W(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                           cg->fgfs[RJ0->sgfn], cg->fgfs[IJ0->sgfn],
                           cg->fgfs[RB->sgfn], cg->fgfs[IB->sgfn],
                           cg->fgfs[Rnu->sgfn], cg->fgfs[Inu->sgfn],
                           cg->fgfs[Rk->sgfn], cg->fgfs[Ik->sgfn],
                           cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn],
                           cg->fgfs[RQ->sgfn], cg->fgfs[IQ->sgfn],
                           cg->fgfs[W->sgfn], cg->fgfs[beta->sgfn], cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn], Rmin,
                           cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                           cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                           cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn]) ||
              // since we do not need derivetive of W, we can deal with Theta together here
              // at this stage W has been updated already
              f_NullEvol_Theta(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                               cg->fgfs[RJ0->sgfn], cg->fgfs[IJ0->sgfn],
                               cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn],
                               cg->fgfs[beta->sgfn],
                               cg->fgfs[RB->sgfn], cg->fgfs[IB->sgfn],
                               cg->fgfs[Rnu->sgfn], cg->fgfs[Inu->sgfn],
                               cg->fgfs[Rk->sgfn], cg->fgfs[Ik->sgfn],
                               cg->fgfs[RTheta->sgfn], cg->fgfs[ITheta->sgfn],
                               cg->fgfs[W->sgfn],
                               cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn], cg->fgfs[KKx->sgfn], cg->fgfs[HKKx->sgfn],
                               Rmin,
                               cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                               cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                               cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn]))
          {
            cout << "find NaN in NullShell domain: sst = " << sPp->data->sst << ", (" << cg->bbox[0] << ":" << cg->bbox[3] << ","
                 << cg->bbox[1] << ":" << cg->bbox[4] << "," << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
            ERROR = 1;
          }
        }
        else
        {
          if (f_NullEvol_W(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                           cg->fgfs[RJ->sgfn], cg->fgfs[IJ->sgfn],
                           cg->fgfs[RB->sgfn], cg->fgfs[IB->sgfn],
                           cg->fgfs[Rnu->sgfn], cg->fgfs[Inu->sgfn],
                           cg->fgfs[Rk->sgfn], cg->fgfs[Ik->sgfn],
                           cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn],
                           cg->fgfs[RQ->sgfn], cg->fgfs[IQ->sgfn],
                           cg->fgfs[W->sgfn], cg->fgfs[beta->sgfn], cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn], Rmin,
                           cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                           cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                           cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn]) ||
              // since we do not need derivetive of W, we can deal with Theta together here
              // at this stage W has been updated already
              f_NullEvol_Theta(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                               cg->fgfs[RJ->sgfn], cg->fgfs[IJ->sgfn],
                               cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn],
                               cg->fgfs[beta->sgfn],
                               cg->fgfs[RB->sgfn], cg->fgfs[IB->sgfn],
                               cg->fgfs[Rnu->sgfn], cg->fgfs[Inu->sgfn],
                               cg->fgfs[Rk->sgfn], cg->fgfs[Ik->sgfn],
                               cg->fgfs[RTheta->sgfn], cg->fgfs[ITheta->sgfn],
                               cg->fgfs[W->sgfn],
                               cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn], cg->fgfs[KKx->sgfn], cg->fgfs[HKKx->sgfn],
                               Rmin,
                               cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                               cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                               cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn]))
          {
            cout << "find NaN in NullShell domain: sst = " << sPp->data->sst << ", (" << cg->bbox[0] << ":" << cg->bbox[3] << ","
                 << cg->bbox[1] << ":" << cg->bbox[4] << "," << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
            ERROR = 1;
          }
        }
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }
  // check error information
  {
    int erh = ERROR;
    MPI_Allreduce(&erh, &ERROR, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  }
  if (ERROR)
  {
    Dump_Data(WTheList, 0, PhysTime, dT);
    if (myrank == 0)
    {
      if (ErrorMonitor && ErrorMonitor->outfile)
        ErrorMonitor->outfile << "find NaN in W and/or Theta on NullShell Patches at t = " << PhysTime << endl;
      else
        cout << "find NaN in W and/or Theta on NullShell Patches at t = " << PhysTime << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  Synch(WTheList, Symmetry, WThewt);
}
#else
#if 0
//For check, give all surface varialbes
//check J evolve only
void NullShellPatch::HyperSlice(double dT,double PhysTime,monitor *ErrorMonitor,int RK_count)
{ 
    int ERROR=0;	

    int spin,e;

    MyList<ss_patch> *sPp;

    sPp=PatL;
    while(sPp)
    {
      MyList<Block> *BP=sPp->data->blb;
      int fngfs = sPp->data->fngfs;
      while(BP)
      {
        Block *cg=BP->data;
        if(myrank == cg->rank) 
        {
/*
          f_get_exact_null_theta(cg->shape,cg->X[0],cg->X[1],cg->X[2],
			            cg->fgfs[RTheta->sgfn],cg->fgfs[ITheta->sgfn],sPp->data->sst,Rmin,PhysTime,
		   	            cg->fgfs[quR1->sgfn],cg->fgfs[quR2->sgfn],cg->fgfs[quI1->sgfn],cg->fgfs[quI2->sgfn],
		   	            cg->fgfs[qlR1->sgfn],cg->fgfs[qlR2->sgfn],cg->fgfs[qlI1->sgfn],cg->fgfs[qlI2->sgfn],
			            cg->fgfs[gR->sgfn],cg->fgfs[gI->sgfn],
			            cg->fgfs[dquR1->sgfn],cg->fgfs[dquR2->sgfn],cg->fgfs[dquI1->sgfn],cg->fgfs[dquI2->sgfn],
			            cg->fgfs[bdquR1->sgfn],cg->fgfs[bdquR2->sgfn],cg->fgfs[bdquI1->sgfn],cg->fgfs[bdquI2->sgfn],
			            cg->fgfs[dgR->sgfn],cg->fgfs[dgI->sgfn],cg->fgfs[bdgR->sgfn],cg->fgfs[bdgI->sgfn]);
*/	  
	   f_get_null_boundary_c(cg->shape,cg->X[0],cg->X[1],cg->X[2],
			            cg->fgfs[beta->sgfn],cg->fgfs[RQ->sgfn],cg->fgfs[IQ->sgfn],
				    cg->fgfs[RU->sgfn],cg->fgfs[IU->sgfn],
			            cg->fgfs[W->sgfn],cg->fgfs[RTheta->sgfn],cg->fgfs[ITheta->sgfn],
		   	            cg->fgfs[quR1->sgfn],cg->fgfs[quR2->sgfn],cg->fgfs[quI1->sgfn],cg->fgfs[quI2->sgfn],
		   	            cg->fgfs[qlR1->sgfn],cg->fgfs[qlR2->sgfn],cg->fgfs[qlI1->sgfn],cg->fgfs[qlI2->sgfn],
			            cg->fgfs[gR->sgfn],cg->fgfs[gI->sgfn],
			            cg->fgfs[dquR1->sgfn],cg->fgfs[dquR2->sgfn],cg->fgfs[dquI1->sgfn],cg->fgfs[dquI2->sgfn],
			            cg->fgfs[bdquR1->sgfn],cg->fgfs[bdquR2->sgfn],cg->fgfs[bdquI1->sgfn],cg->fgfs[bdquI2->sgfn],
			            cg->fgfs[dgR->sgfn],cg->fgfs[dgI->sgfn],cg->fgfs[bdgR->sgfn],cg->fgfs[bdgI->sgfn],
				    PhysTime,Rmin,sPp->data->sst);

	}
        if(BP==sPp->data->ble) break;
        BP=BP->next;
      }
      sPp=sPp->next;
    }
}
#elif 0
// For check Theta calculation with given Theta_x
void NullShellPatch::HyperSlice(double dT, double PhysTime, monitor *ErrorMonitor, int RK_count)
{
  int ERROR = 0;

  int spin, e;

  MyList<ss_patch> *sPp;

  // calculate K
  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        f_get_null_boundary_c(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                              cg->fgfs[beta->sgfn], cg->fgfs[RQ->sgfn], cg->fgfs[IQ->sgfn],
                              cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn],
                              cg->fgfs[W->sgfn], cg->fgfs[RTheta->sgfn], cg->fgfs[ITheta->sgfn],
                              cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                              cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                              cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn],
                              cg->fgfs[dquR1->sgfn], cg->fgfs[dquR2->sgfn], cg->fgfs[dquI1->sgfn], cg->fgfs[dquI2->sgfn],
                              cg->fgfs[bdquR1->sgfn], cg->fgfs[bdquR2->sgfn], cg->fgfs[bdquI1->sgfn], cg->fgfs[bdquI2->sgfn],
                              cg->fgfs[dgR->sgfn], cg->fgfs[dgI->sgfn], cg->fgfs[bdgR->sgfn], cg->fgfs[bdgI->sgfn],
                              PhysTime, Rmin, sPp->data->sst);
        if (RK_count == 0)
        {
          f_calculate_K(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                        cg->fgfs[RJ0->sgfn], cg->fgfs[IJ0->sgfn], cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn], cg->fgfs[KKx->sgfn], cg->fgfs[HKKx->sgfn]);
        }
        else
        {
          f_calculate_K(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                        cg->fgfs[RJ->sgfn], cg->fgfs[IJ->sgfn], cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn], cg->fgfs[KKx->sgfn], cg->fgfs[HKKx->sgfn]);
        }
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }

  if (0)
  {
    int Varwt[3];
    MyList<var> *DG_List;
    DG_List = new MyList<var>(RU);
    DG_List->insert(IU);
    Varwt[0] = 1;
    DG_List->insert(RQ);
    DG_List->insert(IQ);
    Varwt[1] = 1;
    DG_List->insert(RTheta);
    DG_List->insert(ITheta);
    Varwt[2] = 2;

    Synch(DG_List, Symmetry, Varwt);
    DG_List->clearList();
  }

  // get nu, k and B
  spin = 2;
  e = -1;
  if (RK_count == 0)
    eth_derivs(RJ0, IJ0, Rnu, Inu, spin, e);
  else
    eth_derivs(RJ, IJ, Rnu, Inu, spin, e);
  spin = 0;
  e = 1;
  eth_derivs(KK, FXZEO, Rk, Ik, spin, e);
  eth_derivs(beta, FXZEO, RB, IB, spin, e);

  // evolve Theta
  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        if (RK_count == 0)
        {
          if (f_NullEvol_Theta_givenx(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                      cg->fgfs[RJ0->sgfn], cg->fgfs[IJ0->sgfn],
                                      cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn],
                                      cg->fgfs[beta->sgfn],
                                      cg->fgfs[RB->sgfn], cg->fgfs[IB->sgfn],
                                      cg->fgfs[Rnu->sgfn], cg->fgfs[Inu->sgfn],
                                      cg->fgfs[Rk->sgfn], cg->fgfs[Ik->sgfn],
                                      cg->fgfs[RTheta->sgfn], cg->fgfs[ITheta->sgfn],
                                      cg->fgfs[W->sgfn],
                                      cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn], cg->fgfs[KKx->sgfn], cg->fgfs[HKKx->sgfn],
                                      Rmin,
                                      cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                                      cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                                      cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn],
                                      cg->fgfs[dquR1->sgfn], cg->fgfs[dquR2->sgfn], cg->fgfs[dquI1->sgfn], cg->fgfs[dquI2->sgfn],
                                      cg->fgfs[bdquR1->sgfn], cg->fgfs[bdquR2->sgfn], cg->fgfs[bdquI1->sgfn], cg->fgfs[bdquI2->sgfn],
                                      cg->fgfs[dgR->sgfn], cg->fgfs[dgI->sgfn], cg->fgfs[bdgR->sgfn], cg->fgfs[bdgI->sgfn],
                                      PhysTime, sPp->data->sst))
          {
            cout << "find NaN in NullShell domain: sst = " << sPp->data->sst << ", (" << cg->bbox[0] << ":" << cg->bbox[3] << ","
                 << cg->bbox[1] << ":" << cg->bbox[4] << "," << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
            ERROR = 1;
          }
        }
        else
        {
          if (f_NullEvol_Theta_givenx(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                      cg->fgfs[RJ->sgfn], cg->fgfs[IJ->sgfn],
                                      cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn],
                                      cg->fgfs[beta->sgfn],
                                      cg->fgfs[RB->sgfn], cg->fgfs[IB->sgfn],
                                      cg->fgfs[Rnu->sgfn], cg->fgfs[Inu->sgfn],
                                      cg->fgfs[Rk->sgfn], cg->fgfs[Ik->sgfn],
                                      cg->fgfs[RTheta->sgfn], cg->fgfs[ITheta->sgfn],
                                      cg->fgfs[W->sgfn],
                                      cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn], cg->fgfs[KKx->sgfn], cg->fgfs[HKKx->sgfn],
                                      Rmin,
                                      cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                                      cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                                      cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn],
                                      cg->fgfs[dquR1->sgfn], cg->fgfs[dquR2->sgfn], cg->fgfs[dquI1->sgfn], cg->fgfs[dquI2->sgfn],
                                      cg->fgfs[bdquR1->sgfn], cg->fgfs[bdquR2->sgfn], cg->fgfs[bdquI1->sgfn], cg->fgfs[bdquI2->sgfn],
                                      cg->fgfs[dgR->sgfn], cg->fgfs[dgI->sgfn], cg->fgfs[bdgR->sgfn], cg->fgfs[bdgI->sgfn],
                                      PhysTime, sPp->data->sst))
          {
            cout << "find NaN in NullShell domain: sst = " << sPp->data->sst << ", (" << cg->bbox[0] << ":" << cg->bbox[3] << ","
                 << cg->bbox[1] << ":" << cg->bbox[4] << "," << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
            ERROR = 1;
          }
        }
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }
  // check error information
  {
    int erh = ERROR;
    MPI_Allreduce(&erh, &ERROR, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  }
  if (ERROR)
  {
    Dump_Data(WTheList, 0, PhysTime, dT);
    if (myrank == 0)
    {
      if (ErrorMonitor && ErrorMonitor->outfile)
        ErrorMonitor->outfile << "find NaN in W and/or Theta on NullShell Patches at t = " << PhysTime << endl;
      else
        cout << "find NaN in W and/or Theta on NullShell Patches at t = " << PhysTime << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  Synch(WTheList, Symmetry, WThewt);
}
#elif 0
// For check Theta calculation
void NullShellPatch::HyperSlice(double dT, double PhysTime, monitor *ErrorMonitor, int RK_count)
{
  int ERROR = 0;

  int spin, e;

  MyList<ss_patch> *sPp;

  // calculate K
  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        f_get_null_boundary_c(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                              cg->fgfs[beta->sgfn], cg->fgfs[RQ->sgfn], cg->fgfs[IQ->sgfn],
                              cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn],
                              cg->fgfs[W->sgfn], cg->fgfs[RTheta->sgfn], cg->fgfs[ITheta->sgfn],
                              cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                              cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                              cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn],
                              cg->fgfs[dquR1->sgfn], cg->fgfs[dquR2->sgfn], cg->fgfs[dquI1->sgfn], cg->fgfs[dquI2->sgfn],
                              cg->fgfs[bdquR1->sgfn], cg->fgfs[bdquR2->sgfn], cg->fgfs[bdquI1->sgfn], cg->fgfs[bdquI2->sgfn],
                              cg->fgfs[dgR->sgfn], cg->fgfs[dgI->sgfn], cg->fgfs[bdgR->sgfn], cg->fgfs[bdgI->sgfn],
                              PhysTime, Rmin, sPp->data->sst);
        if (RK_count == 0)
        {
          f_calculate_K(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                        cg->fgfs[RJ0->sgfn], cg->fgfs[IJ0->sgfn], cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn], cg->fgfs[KKx->sgfn], cg->fgfs[HKKx->sgfn]);
        }
        else
        {
          f_calculate_K(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                        cg->fgfs[RJ->sgfn], cg->fgfs[IJ->sgfn], cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn], cg->fgfs[KKx->sgfn], cg->fgfs[HKKx->sgfn]);
        }
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }

  if (0)
  {
    int Varwt[3];
    MyList<var> *DG_List;
    DG_List = new MyList<var>(RU);
    DG_List->insert(IU);
    Varwt[0] = 1;
    DG_List->insert(RQ);
    DG_List->insert(IQ);
    Varwt[1] = 1;
    DG_List->insert(RTheta);
    DG_List->insert(ITheta);
    Varwt[2] = 2;

    Synch(DG_List, Symmetry, Varwt);
    DG_List->clearList();
  }

  // get nu, k and B
  spin = 2;
  e = -1;
  if (RK_count == 0)
    eth_derivs(RJ0, IJ0, Rnu, Inu, spin, e);
  else
    eth_derivs(RJ, IJ, Rnu, Inu, spin, e);
  spin = 0;
  e = 1;
  eth_derivs(KK, FXZEO, Rk, Ik, spin, e);
  eth_derivs(beta, FXZEO, RB, IB, spin, e);

  // evolve Theta
  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        if (RK_count == 0)
        {
          if (f_NullEvol_Theta(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                               cg->fgfs[RJ0->sgfn], cg->fgfs[IJ0->sgfn],
                               cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn],
                               cg->fgfs[beta->sgfn],
                               cg->fgfs[RB->sgfn], cg->fgfs[IB->sgfn],
                               cg->fgfs[Rnu->sgfn], cg->fgfs[Inu->sgfn],
                               cg->fgfs[Rk->sgfn], cg->fgfs[Ik->sgfn],
                               cg->fgfs[RTheta->sgfn], cg->fgfs[ITheta->sgfn],
                               cg->fgfs[W->sgfn],
                               cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn], cg->fgfs[KKx->sgfn], cg->fgfs[HKKx->sgfn],
                               Rmin,
                               cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                               cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                               cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn]))
          {
            cout << "find NaN in NullShell domain: sst = " << sPp->data->sst << ", (" << cg->bbox[0] << ":" << cg->bbox[3] << ","
                 << cg->bbox[1] << ":" << cg->bbox[4] << "," << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
            ERROR = 1;
          }
        }
        else
        {
          if (f_NullEvol_Theta(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                               cg->fgfs[RJ->sgfn], cg->fgfs[IJ->sgfn],
                               cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn],
                               cg->fgfs[beta->sgfn],
                               cg->fgfs[RB->sgfn], cg->fgfs[IB->sgfn],
                               cg->fgfs[Rnu->sgfn], cg->fgfs[Inu->sgfn],
                               cg->fgfs[Rk->sgfn], cg->fgfs[Ik->sgfn],
                               cg->fgfs[RTheta->sgfn], cg->fgfs[ITheta->sgfn],
                               cg->fgfs[W->sgfn],
                               cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn], cg->fgfs[KKx->sgfn], cg->fgfs[HKKx->sgfn],
                               Rmin,
                               cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                               cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                               cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn]))
          {
            cout << "find NaN in NullShell domain: sst = " << sPp->data->sst << ", (" << cg->bbox[0] << ":" << cg->bbox[3] << ","
                 << cg->bbox[1] << ":" << cg->bbox[4] << "," << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
            ERROR = 1;
          }
        }
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }
  // check error information
  {
    int erh = ERROR;
    MPI_Allreduce(&erh, &ERROR, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  }
  if (ERROR)
  {
    Dump_Data(WTheList, 0, PhysTime, dT);
    if (myrank == 0)
    {
      if (ErrorMonitor && ErrorMonitor->outfile)
        ErrorMonitor->outfile << "find NaN in W and/or Theta on NullShell Patches at t = " << PhysTime << endl;
      else
        cout << "find NaN in W and/or Theta on NullShell Patches at t = " << PhysTime << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  Synch(WTheList, Symmetry, WThewt);
}
#elif 0
// For check W and Theta calculation
void NullShellPatch::HyperSlice(double dT, double PhysTime, monitor *ErrorMonitor, int RK_count)
{
  int ERROR = 0;

  int spin, e;

  MyList<ss_patch> *sPp;

  // calculate K
  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        f_get_null_boundary_c(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                              cg->fgfs[beta->sgfn], cg->fgfs[RQ->sgfn], cg->fgfs[IQ->sgfn],
                              cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn],
                              cg->fgfs[W->sgfn], cg->fgfs[RTheta->sgfn], cg->fgfs[ITheta->sgfn],
                              cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                              cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                              cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn],
                              cg->fgfs[dquR1->sgfn], cg->fgfs[dquR2->sgfn], cg->fgfs[dquI1->sgfn], cg->fgfs[dquI2->sgfn],
                              cg->fgfs[bdquR1->sgfn], cg->fgfs[bdquR2->sgfn], cg->fgfs[bdquI1->sgfn], cg->fgfs[bdquI2->sgfn],
                              cg->fgfs[dgR->sgfn], cg->fgfs[dgI->sgfn], cg->fgfs[bdgR->sgfn], cg->fgfs[bdgI->sgfn],
                              PhysTime, Rmin, sPp->data->sst);
        if (RK_count == 0)
        {
          f_calculate_K(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                        cg->fgfs[RJ0->sgfn], cg->fgfs[IJ0->sgfn], cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn], cg->fgfs[KKx->sgfn], cg->fgfs[HKKx->sgfn]);
        }
        else
        {
          f_calculate_K(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                        cg->fgfs[RJ->sgfn], cg->fgfs[IJ->sgfn], cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn], cg->fgfs[KKx->sgfn], cg->fgfs[HKKx->sgfn]);
        }
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }

  {
    int Varwt[3];
    MyList<var> *DG_List;
    DG_List = new MyList<var>(RU);
    DG_List->insert(IU);
    Varwt[0] = 1;
    DG_List->insert(RQ);
    DG_List->insert(IQ);
    Varwt[1] = 1;
    DG_List->insert(RTheta);
    DG_List->insert(ITheta);
    Varwt[2] = 2;

    Synch(DG_List, Symmetry, Varwt);
    DG_List->clearList();
  }

  // get nu, k and B
  spin = 2;
  e = -1;
  if (RK_count == 0)
    eth_derivs(RJ0, IJ0, Rnu, Inu, spin, e);
  else
    eth_derivs(RJ, IJ, Rnu, Inu, spin, e);
  spin = 0;
  e = 1;
  eth_derivs(KK, FXZEO, Rk, Ik, spin, e);
  eth_derivs(beta, FXZEO, RB, IB, spin, e);

  // evolve W and Theta
  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        if (RK_count == 0)
        {
          if (f_NullEvol_W(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                           cg->fgfs[RJ0->sgfn], cg->fgfs[IJ0->sgfn],
                           cg->fgfs[RB->sgfn], cg->fgfs[IB->sgfn],
                           cg->fgfs[Rnu->sgfn], cg->fgfs[Inu->sgfn],
                           cg->fgfs[Rk->sgfn], cg->fgfs[Ik->sgfn],
                           cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn],
                           cg->fgfs[RQ->sgfn], cg->fgfs[IQ->sgfn],
                           cg->fgfs[W->sgfn], cg->fgfs[beta->sgfn], cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn], Rmin,
                           cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                           cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                           cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn]) ||
              // since we do not need derivetive of W, we can deal with Theta together here
              // at this stage W has been updated already
              f_NullEvol_Theta(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                               cg->fgfs[RJ0->sgfn], cg->fgfs[IJ0->sgfn],
                               cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn],
                               cg->fgfs[beta->sgfn],
                               cg->fgfs[RB->sgfn], cg->fgfs[IB->sgfn],
                               cg->fgfs[Rnu->sgfn], cg->fgfs[Inu->sgfn],
                               cg->fgfs[Rk->sgfn], cg->fgfs[Ik->sgfn],
                               cg->fgfs[RTheta->sgfn], cg->fgfs[ITheta->sgfn],
                               cg->fgfs[W->sgfn],
                               cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn], cg->fgfs[KKx->sgfn], cg->fgfs[HKKx->sgfn],
                               Rmin,
                               cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                               cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                               cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn]))
          {
            cout << "find NaN in NullShell domain: sst = " << sPp->data->sst << ", (" << cg->bbox[0] << ":" << cg->bbox[3] << ","
                 << cg->bbox[1] << ":" << cg->bbox[4] << "," << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
            ERROR = 1;
          }
        }
        else
        {
          if (f_NullEvol_W(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                           cg->fgfs[RJ->sgfn], cg->fgfs[IJ->sgfn],
                           cg->fgfs[RB->sgfn], cg->fgfs[IB->sgfn],
                           cg->fgfs[Rnu->sgfn], cg->fgfs[Inu->sgfn],
                           cg->fgfs[Rk->sgfn], cg->fgfs[Ik->sgfn],
                           cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn],
                           cg->fgfs[RQ->sgfn], cg->fgfs[IQ->sgfn],
                           cg->fgfs[W->sgfn], cg->fgfs[beta->sgfn], cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn], Rmin,
                           cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                           cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                           cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn]) ||
              // since we do not need derivetive of W, we can deal with Theta together here
              // at this stage W has been updated already
              f_NullEvol_Theta(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                               cg->fgfs[RJ->sgfn], cg->fgfs[IJ->sgfn],
                               cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn],
                               cg->fgfs[beta->sgfn],
                               cg->fgfs[RB->sgfn], cg->fgfs[IB->sgfn],
                               cg->fgfs[Rnu->sgfn], cg->fgfs[Inu->sgfn],
                               cg->fgfs[Rk->sgfn], cg->fgfs[Ik->sgfn],
                               cg->fgfs[RTheta->sgfn], cg->fgfs[ITheta->sgfn],
                               cg->fgfs[W->sgfn],
                               cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn], cg->fgfs[KKx->sgfn], cg->fgfs[HKKx->sgfn],
                               Rmin,
                               cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                               cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                               cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn]))
          {
            cout << "find NaN in NullShell domain: sst = " << sPp->data->sst << ", (" << cg->bbox[0] << ":" << cg->bbox[3] << ","
                 << cg->bbox[1] << ":" << cg->bbox[4] << "," << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
            ERROR = 1;
          }
        }
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }
  // check error information
  {
    int erh = ERROR;
    MPI_Allreduce(&erh, &ERROR, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  }
  if (ERROR)
  {
    Dump_Data(WTheList, 0, PhysTime, dT);
    if (myrank == 0)
    {
      if (ErrorMonitor && ErrorMonitor->outfile)
        ErrorMonitor->outfile << "find NaN in W and/or Theta on NullShell Patches at t = " << PhysTime << endl;
      else
        cout << "find NaN in W and/or Theta on NullShell Patches at t = " << PhysTime << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  Synch(WTheList, Symmetry, WThewt);
}
#elif 1
// For check Q, U, W and Theta calculation
void NullShellPatch::HyperSlice(double dT, double PhysTime, monitor *ErrorMonitor, int RK_count)
{
  int ERROR = 0;

  int spin, e;

  MyList<ss_patch> *sPp;

  // calculate K
  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        f_get_null_boundary_c(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                              cg->fgfs[beta->sgfn], cg->fgfs[RQ->sgfn], cg->fgfs[IQ->sgfn],
                              cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn],
                              cg->fgfs[W->sgfn], cg->fgfs[RTheta->sgfn], cg->fgfs[ITheta->sgfn],
                              cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                              cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                              cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn],
                              cg->fgfs[dquR1->sgfn], cg->fgfs[dquR2->sgfn], cg->fgfs[dquI1->sgfn], cg->fgfs[dquI2->sgfn],
                              cg->fgfs[bdquR1->sgfn], cg->fgfs[bdquR2->sgfn], cg->fgfs[bdquI1->sgfn], cg->fgfs[bdquI2->sgfn],
                              cg->fgfs[dgR->sgfn], cg->fgfs[dgI->sgfn], cg->fgfs[bdgR->sgfn], cg->fgfs[bdgI->sgfn],
                              PhysTime, Rmin, sPp->data->sst);
        if (RK_count == 0)
        {
          f_calculate_K(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                        cg->fgfs[RJ0->sgfn], cg->fgfs[IJ0->sgfn], cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn], cg->fgfs[KKx->sgfn], cg->fgfs[HKKx->sgfn]);
        }
        else
        {
          f_calculate_K(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                        cg->fgfs[RJ->sgfn], cg->fgfs[IJ->sgfn], cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn], cg->fgfs[KKx->sgfn], cg->fgfs[HKKx->sgfn]);
        }
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }

  {
    int Varwt[3];
    MyList<var> *DG_List;
    DG_List = new MyList<var>(RU);
    DG_List->insert(IU);
    Varwt[0] = 1;
    DG_List->insert(RQ);
    DG_List->insert(IQ);
    Varwt[1] = 1;
    DG_List->insert(RTheta);
    DG_List->insert(ITheta);
    Varwt[2] = 2;

    Synch(DG_List, Symmetry, Varwt);
    DG_List->clearList();
  }

  // get nu, k and B
  spin = 2;
  e = -1;
  if (RK_count == 0)
    eth_derivs(RJ0, IJ0, Rnu, Inu, spin, e);
  else
    eth_derivs(RJ, IJ, Rnu, Inu, spin, e);
  spin = 0;
  e = 1;
  eth_derivs(KK, FXZEO, Rk, Ik, spin, e);
  eth_derivs(beta, FXZEO, RB, IB, spin, e);

  // evolve Q and U
  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        if (RK_count == 0)
        {
          if (f_NullEvol_Q(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                           cg->fgfs[RJ0->sgfn], cg->fgfs[IJ0->sgfn],
                           cg->fgfs[Rk->sgfn], cg->fgfs[Ik->sgfn],
                           cg->fgfs[Rnu->sgfn], cg->fgfs[Inu->sgfn],
                           cg->fgfs[RB->sgfn], cg->fgfs[IB->sgfn],
                           cg->fgfs[RQ->sgfn], cg->fgfs[IQ->sgfn],
                           cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn], cg->fgfs[KKx->sgfn], cg->fgfs[HKKx->sgfn],
                           cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                           cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                           cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn]) ||
              // since we do not need derivetive of Q, we can deal with U together here
              // at this stage Q has been updated already
              f_NullEvol_U(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                           cg->fgfs[RJ0->sgfn], cg->fgfs[IJ0->sgfn],
                           cg->fgfs[RQ->sgfn], cg->fgfs[IQ->sgfn], cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn],
                           cg->fgfs[beta->sgfn], cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn], Rmin))
          {
            cout << "find NaN in NullShell domain: sst = " << sPp->data->sst << ", (" << cg->bbox[0] << ":" << cg->bbox[3] << ","
                 << cg->bbox[1] << ":" << cg->bbox[4] << "," << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
            ERROR = 1;
          }
        }
        else
        {
          if (f_NullEvol_Q(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                           cg->fgfs[RJ->sgfn], cg->fgfs[IJ->sgfn],
                           cg->fgfs[Rk->sgfn], cg->fgfs[Ik->sgfn],
                           cg->fgfs[Rnu->sgfn], cg->fgfs[Inu->sgfn],
                           cg->fgfs[RB->sgfn], cg->fgfs[IB->sgfn],
                           cg->fgfs[RQ->sgfn], cg->fgfs[IQ->sgfn],
                           cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn], cg->fgfs[KKx->sgfn], cg->fgfs[HKKx->sgfn],
                           cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                           cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                           cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn]) ||
              // since we do not need derivetive of Q, we can deal with U together here
              // at this stage Q has been updated already
              f_NullEvol_U(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                           cg->fgfs[RJ->sgfn], cg->fgfs[IJ->sgfn],
                           cg->fgfs[RQ->sgfn], cg->fgfs[IQ->sgfn], cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn],
                           cg->fgfs[beta->sgfn], cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn], Rmin))
          {
            cout << "find NaN in NullShell domain: sst = " << sPp->data->sst << ", (" << cg->bbox[0] << ":" << cg->bbox[3] << ","
                 << cg->bbox[1] << ":" << cg->bbox[4] << "," << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
            ERROR = 1;
          }
        }
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }
  // check error information
  {
    int erh = ERROR;
    MPI_Allreduce(&erh, &ERROR, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  }
  if (ERROR)
  {
    Dump_Data(QUList, 0, PhysTime, dT);
    Dump_Data(SynchList_pre, 0, PhysTime, dT);
    if (myrank == 0)
    {
      if (ErrorMonitor && ErrorMonitor->outfile)
        ErrorMonitor->outfile << "find NaN in Q and/or U on NullShell Patches at t = " << PhysTime << endl;
      else
        cout << "find NaN in Q and/or U on NullShell Patches at t = " << PhysTime << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  Synch(QUList, Symmetry, QUwt);

  // evolve W and Theta
  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        if (RK_count == 0)
        {
          if (f_NullEvol_W(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                           cg->fgfs[RJ0->sgfn], cg->fgfs[IJ0->sgfn],
                           cg->fgfs[RB->sgfn], cg->fgfs[IB->sgfn],
                           cg->fgfs[Rnu->sgfn], cg->fgfs[Inu->sgfn],
                           cg->fgfs[Rk->sgfn], cg->fgfs[Ik->sgfn],
                           cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn],
                           cg->fgfs[RQ->sgfn], cg->fgfs[IQ->sgfn],
                           cg->fgfs[W->sgfn], cg->fgfs[beta->sgfn], cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn], Rmin,
                           cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                           cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                           cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn]) ||
              // since we do not need derivetive of W, we can deal with Theta together here
              // at this stage W has been updated already
              f_NullEvol_Theta(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                               cg->fgfs[RJ0->sgfn], cg->fgfs[IJ0->sgfn],
                               cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn],
                               cg->fgfs[beta->sgfn],
                               cg->fgfs[RB->sgfn], cg->fgfs[IB->sgfn],
                               cg->fgfs[Rnu->sgfn], cg->fgfs[Inu->sgfn],
                               cg->fgfs[Rk->sgfn], cg->fgfs[Ik->sgfn],
                               cg->fgfs[RTheta->sgfn], cg->fgfs[ITheta->sgfn],
                               cg->fgfs[W->sgfn],
                               cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn], cg->fgfs[KKx->sgfn], cg->fgfs[HKKx->sgfn],
                               Rmin,
                               cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                               cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                               cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn]))
          {
            cout << "find NaN in NullShell domain: sst = " << sPp->data->sst << ", (" << cg->bbox[0] << ":" << cg->bbox[3] << ","
                 << cg->bbox[1] << ":" << cg->bbox[4] << "," << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
            ERROR = 1;
          }
        }
        else
        {
          if (f_NullEvol_W(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                           cg->fgfs[RJ->sgfn], cg->fgfs[IJ->sgfn],
                           cg->fgfs[RB->sgfn], cg->fgfs[IB->sgfn],
                           cg->fgfs[Rnu->sgfn], cg->fgfs[Inu->sgfn],
                           cg->fgfs[Rk->sgfn], cg->fgfs[Ik->sgfn],
                           cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn],
                           cg->fgfs[RQ->sgfn], cg->fgfs[IQ->sgfn],
                           cg->fgfs[W->sgfn], cg->fgfs[beta->sgfn], cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn], Rmin,
                           cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                           cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                           cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn]) ||
              // since we do not need derivetive of W, we can deal with Theta together here
              // at this stage W has been updated already
              f_NullEvol_Theta(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                               cg->fgfs[RJ->sgfn], cg->fgfs[IJ->sgfn],
                               cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn],
                               cg->fgfs[beta->sgfn],
                               cg->fgfs[RB->sgfn], cg->fgfs[IB->sgfn],
                               cg->fgfs[Rnu->sgfn], cg->fgfs[Inu->sgfn],
                               cg->fgfs[Rk->sgfn], cg->fgfs[Ik->sgfn],
                               cg->fgfs[RTheta->sgfn], cg->fgfs[ITheta->sgfn],
                               cg->fgfs[W->sgfn],
                               cg->fgfs[KK->sgfn], cg->fgfs[HKK->sgfn], cg->fgfs[KKx->sgfn], cg->fgfs[HKKx->sgfn],
                               Rmin,
                               cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                               cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                               cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn]))
          {
            cout << "find NaN in NullShell domain: sst = " << sPp->data->sst << ", (" << cg->bbox[0] << ":" << cg->bbox[3] << ","
                 << cg->bbox[1] << ":" << cg->bbox[4] << "," << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
            ERROR = 1;
          }
        }
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }

  // check error information
  {
    int erh = ERROR;
    MPI_Allreduce(&erh, &ERROR, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  }
  if (ERROR)
  {
    Dump_Data(WTheList, 0, PhysTime, dT);
    Dump_Data(QUList, 0, PhysTime, dT);
    if (myrank == 0)
    {
      if (ErrorMonitor && ErrorMonitor->outfile)
        ErrorMonitor->outfile << "find NaN in W and/or Theta on NullShell Patches at t = " << PhysTime << endl;
      else
        cout << "find NaN in W and/or Theta on NullShell Patches at t = " << PhysTime << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  Synch(WTheList, Symmetry, WThewt);
}
#endif
#endif
#if 1
// need evolve step
// 0: real L2 norm; 1: root mean squar
#define L2m 0
double NullShellPatch::Error_Check(double PhysTime, double dT, bool dp)
{
  MyList<ss_patch> *sPp;

  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        f_get_exact_null(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                         cg->fgfs[RJ1->sgfn], cg->fgfs[IJ1->sgfn], sPp->data->sst, Rmin, PhysTime,
                         cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                         cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                         cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn],
                         cg->fgfs[dquR1->sgfn], cg->fgfs[dquR2->sgfn], cg->fgfs[dquI1->sgfn], cg->fgfs[dquI2->sgfn],
                         cg->fgfs[bdquR1->sgfn], cg->fgfs[bdquR2->sgfn], cg->fgfs[bdquI1->sgfn], cg->fgfs[bdquI2->sgfn],
                         cg->fgfs[dgR->sgfn], cg->fgfs[dgI->sgfn], cg->fgfs[bdgR->sgfn], cg->fgfs[bdgI->sgfn]);
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }

  if (0)
  {
    int Varwt[1];
    MyList<var> *DG_List;
    DG_List = new MyList<var>(RJ1);
    DG_List->insert(IJ1);
    Varwt[0] = 2;
    Synch(DG_List, Symmetry, Varwt);

    if (dp)
    {
      DG_List->insert(RJ0);
      DG_List->insert(IJ0);
      Dump_Data(DG_List, 0, PhysTime, dT);
    }
    DG_List->clearList();
  }

  double tvf, dtvf = 0;
  int tN, dtN = 0;
  int BDW = ghost_width, OBDW = overghost;
  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        f_array_subtract(cg->shape, cg->fgfs[RJ1->sgfn], cg->fgfs[RJ0->sgfn]);
#if (L2m == 0)
        f_l2normhelper_sh(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                          sPp->data->bbox[0], sPp->data->bbox[1], sPp->data->bbox[2], sPp->data->bbox[3], sPp->data->bbox[4], sPp->data->bbox[5],
                          cg->fgfs[RJ1->sgfn], tvf, BDW, OBDW, Symmetry);
#elif (L2m == 1)
        f_l2normhelper_sh_rms(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                              sPp->data->bbox[0], sPp->data->bbox[1], sPp->data->bbox[2], sPp->data->bbox[3], sPp->data->bbox[4], sPp->data->bbox[5],
                              cg->fgfs[RJ1->sgfn], tvf, BDW, OBDW, Symmetry, dtN);
        dtN += dtN;
#endif

        dtvf += tvf;
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }

  MPI_Allreduce(&dtvf, &tvf, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#if (L2m == 0)
  tvf = sqrt(tvf);
#elif (L2m == 1)
  MPI_Allreduce(&dtN, &tN, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  tvf = sqrt(tvf / tN);
#endif
#if 0 
    {
     MyList<var> * DG_List;    
     DG_List=new MyList<var>(RJ1); DG_List->insert(IJ1);
     
     Dump_Data(DG_List,0,0,1);
     DG_List->clearList(); 
     if(myrank==0) MPI_Abort(MPI_COMM_WORLD,1);
    }
#endif

  return tvf;
}
#else
// only check Theta calculation, do not need Evolve step
double NullShellPatch::Error_Check(double PhysTime, double dT, bool dp)
{
  MyList<ss_patch> *sPp;

  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        f_get_exact_null(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                         cg->fgfs[RJ0->sgfn], cg->fgfs[IJ0->sgfn], sPp->data->sst, Rmin, PhysTime,
                         cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                         cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                         cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn],
                         cg->fgfs[dquR1->sgfn], cg->fgfs[dquR2->sgfn], cg->fgfs[dquI1->sgfn], cg->fgfs[dquI2->sgfn],
                         cg->fgfs[bdquR1->sgfn], cg->fgfs[bdquR2->sgfn], cg->fgfs[bdquI1->sgfn], cg->fgfs[bdquI2->sgfn],
                         cg->fgfs[dgR->sgfn], cg->fgfs[dgI->sgfn], cg->fgfs[bdgR->sgfn], cg->fgfs[bdgI->sgfn]);
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }

  {
    int Varwt[1];
    MyList<var> *DG_List;
    DG_List = new MyList<var>(RJ0);
    DG_List->insert(IJ0);
    Varwt[0] = 2;
    Synch(DG_List, Symmetry, Varwt);
    DG_List->clearList();
  }

  HyperSlice(dT, PhysTime, 0, 0);

  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        f_get_null_boundary_c(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                              cg->fgfs[beta->sgfn], cg->fgfs[RQ->sgfn], cg->fgfs[IQ->sgfn],
                              cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn],
                              cg->fgfs[W->sgfn], cg->fgfs[RJ1->sgfn], cg->fgfs[IJ1->sgfn],
                              cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                              cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                              cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn],
                              cg->fgfs[dquR1->sgfn], cg->fgfs[dquR2->sgfn], cg->fgfs[dquI1->sgfn], cg->fgfs[dquI2->sgfn],
                              cg->fgfs[bdquR1->sgfn], cg->fgfs[bdquR2->sgfn], cg->fgfs[bdquI1->sgfn], cg->fgfs[bdquI2->sgfn],
                              cg->fgfs[dgR->sgfn], cg->fgfs[dgI->sgfn], cg->fgfs[bdgR->sgfn], cg->fgfs[bdgI->sgfn],
                              PhysTime, Rmin, sPp->data->sst);
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }

  {
    int Varwt[1];
    MyList<var> *DG_List;
    DG_List = new MyList<var>(RJ1);
    DG_List->insert(IJ1);
    Varwt[0] = 2;
    Synch(DG_List, Symmetry, Varwt);

    if (dp)
    {
      DG_List->insert(RTheta);
      DG_List->insert(ITheta);
      Dump_Data(DG_List, 0, PhysTime, dT);
    }
    DG_List->clearList();
  }

  double tvf, dtvf = 0;
  int BDW = ghost_width, OBDW = overghost;
  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        f_array_subtract(cg->shape, cg->fgfs[RJ1->sgfn], cg->fgfs[RTheta->sgfn]);

        f_l2normhelper_sh(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                          sPp->data->bbox[0], sPp->data->bbox[1], sPp->data->bbox[2], sPp->data->bbox[3], sPp->data->bbox[4], sPp->data->bbox[5],
                          cg->fgfs[RJ1->sgfn], tvf, BDW, OBDW, Symmetry);
        dtvf += tvf;
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }

  MPI_Allreduce(&dtvf, &tvf, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  tvf = sqrt(tvf);

  return tvf;
}
#endif
double NullShellPatch::EqTheta_Check(double PhysTime, double dT, bool dp)
{
  int ERROR = 0;

  MyList<ss_patch> *sPp;

  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        f_get_exact_null(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                         cg->fgfs[RJ0->sgfn], cg->fgfs[IJ0->sgfn], sPp->data->sst, Rmin, PhysTime,
                         cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                         cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                         cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn],
                         cg->fgfs[dquR1->sgfn], cg->fgfs[dquR2->sgfn], cg->fgfs[dquI1->sgfn], cg->fgfs[dquI2->sgfn],
                         cg->fgfs[bdquR1->sgfn], cg->fgfs[bdquR2->sgfn], cg->fgfs[bdquI1->sgfn], cg->fgfs[bdquI2->sgfn],
                         cg->fgfs[dgR->sgfn], cg->fgfs[dgI->sgfn], cg->fgfs[bdgR->sgfn], cg->fgfs[bdgI->sgfn]);
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }

  {
    int Varwt[1];
    MyList<var> *DG_List;
    DG_List = new MyList<var>(RJ0);
    DG_List->insert(IJ0);
    Varwt[0] = 2;
    Synch(DG_List, Symmetry, Varwt);
    DG_List->clearList();
  }

  HyperSlice(dT, PhysTime, 0, 0);

  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        f_get_null_boundary_c(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                              cg->fgfs[beta->sgfn], cg->fgfs[RQ->sgfn], cg->fgfs[IQ->sgfn],
                              cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn],
                              cg->fgfs[W->sgfn], cg->fgfs[RTheta->sgfn], cg->fgfs[ITheta->sgfn],
                              cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                              cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                              cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn],
                              cg->fgfs[dquR1->sgfn], cg->fgfs[dquR2->sgfn], cg->fgfs[dquI1->sgfn], cg->fgfs[dquI2->sgfn],
                              cg->fgfs[bdquR1->sgfn], cg->fgfs[bdquR2->sgfn], cg->fgfs[bdquI1->sgfn], cg->fgfs[bdquI2->sgfn],
                              cg->fgfs[dgR->sgfn], cg->fgfs[dgI->sgfn], cg->fgfs[bdgR->sgfn], cg->fgfs[bdgI->sgfn],
                              PhysTime, Rmin, sPp->data->sst);
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }

  {
    int Varwt[1];
    MyList<var> *DG_List;
    DG_List = new MyList<var>(RTheta);
    DG_List->insert(ITheta);
    Varwt[0] = 2;
    Synch(DG_List, Symmetry, Varwt);

    DG_List->clearList();
  }

  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        if (f_Eq_Theta(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                       cg->fgfs[RJ0->sgfn], cg->fgfs[IJ0->sgfn],
                       cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn],
                       cg->fgfs[beta->sgfn],
                       cg->fgfs[RB->sgfn], cg->fgfs[IB->sgfn],
                       cg->fgfs[Rnu->sgfn], cg->fgfs[Inu->sgfn],
                       cg->fgfs[Rk->sgfn], cg->fgfs[Ik->sgfn],
                       cg->fgfs[RTheta->sgfn], cg->fgfs[ITheta->sgfn],
                       cg->fgfs[W->sgfn],
                       Rmin,
                       cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                       cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                       cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn]))
        /*		if(f_Eq_Theta_2(cg->shape,cg->X[0],cg->X[1],cg->X[2],
                     cg->fgfs[RJ0->sgfn],cg->fgfs[IJ0->sgfn],
                     cg->fgfs[RU->sgfn],cg->fgfs[IU->sgfn],
                               cg->fgfs[beta->sgfn],
                   cg->fgfs[RB->sgfn],cg->fgfs[IB->sgfn],
                   cg->fgfs[Rnu->sgfn],cg->fgfs[Inu->sgfn],
                   cg->fgfs[Rk->sgfn],cg->fgfs[Ik->sgfn],
                   cg->fgfs[RTheta->sgfn],cg->fgfs[ITheta->sgfn],
                               cg->fgfs[W->sgfn],
                   Rmin,
                               cg->fgfs[qlR1->sgfn],cg->fgfs[qlR2->sgfn],cg->fgfs[qlI1->sgfn],cg->fgfs[qlI2->sgfn],
                               cg->fgfs[quR1->sgfn],cg->fgfs[quR2->sgfn],cg->fgfs[quI1->sgfn],cg->fgfs[quI2->sgfn],
                         cg->fgfs[gR->sgfn],cg->fgfs[gI->sgfn],PhysTime,sPp->data->sst))	*/
        {
          cout << "find NaN in NullShell domain: sst = " << sPp->data->sst << ", (" << cg->bbox[0] << ":" << cg->bbox[3] << ","
               << cg->bbox[1] << ":" << cg->bbox[4] << "," << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
          ERROR = 1;
        }
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }
  // check error information
  {
    int erh = ERROR;
    MPI_Allreduce(&erh, &ERROR, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  }
  if (ERROR)
  {
    Dump_Data(WTheList, 0, PhysTime, dT);
    if (myrank == 0)
    {
      cout << "find NaN in W and/or Theta on NullShell Patches at t = " << PhysTime << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  Synch(WTheList, Symmetry, WThewt);

  if (dp)
  {
    MyList<var> *DG_List;
    DG_List = new MyList<var>(RTheta);
    DG_List->insert(ITheta);
    Dump_Data(DG_List, 0, PhysTime, dT);
    DG_List->clearList();
  }

  double tvf, dtvf = 0;
  int BDW = ghost_width, OBDW = overghost;
  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        f_l2normhelper_sh(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                          sPp->data->bbox[0], sPp->data->bbox[1], sPp->data->bbox[2], sPp->data->bbox[3], sPp->data->bbox[4], sPp->data->bbox[5],
                          cg->fgfs[RTheta->sgfn], tvf, BDW, OBDW, Symmetry);
        dtvf += tvf;
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }

  MPI_Allreduce(&dtvf, &tvf, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  tvf = sqrt(tvf);

  return tvf;
}
void NullShellPatch::Compute_News(double PhysTime, double dT, bool dp)
{
  MyList<ss_patch> *sPp;

  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
// for check
#if 0
     		f_get_exact_omega(cg->shape,cg->X[0],cg->X[1],cg->X[2],
                                  cg->fgfs[omega0->sgfn],sPp->data->sst,Rmin,PhysTime);
#endif
#if 1
        f_drive_null_news(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                          cg->fgfs[RJ0->sgfn], cg->fgfs[IJ0->sgfn], cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn],
                          cg->fgfs[RTheta->sgfn], cg->fgfs[ITheta->sgfn],
                          cg->fgfs[omega0->sgfn], cg->fgfs[beta->sgfn],
                          cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                          cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                          cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn],
                          cg->fgfs[dquR1->sgfn], cg->fgfs[dquR2->sgfn], cg->fgfs[dquI1->sgfn], cg->fgfs[dquI2->sgfn],
                          cg->fgfs[bdquR1->sgfn], cg->fgfs[bdquR2->sgfn], cg->fgfs[bdquI1->sgfn], cg->fgfs[bdquI2->sgfn],
                          cg->fgfs[dgR->sgfn], cg->fgfs[dgI->sgfn], cg->fgfs[bdgR->sgfn], cg->fgfs[bdgI->sgfn],
                          cg->fgfs[RNews->sgfn], cg->fgfs[INews->sgfn], Rmin, sPp->data->sst);
#else
        f_drive_null_news_diff(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                               cg->fgfs[RJ0->sgfn], cg->fgfs[IJ0->sgfn], cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn],
                               cg->fgfs[RTheta->sgfn], cg->fgfs[ITheta->sgfn],
                               cg->fgfs[omega0->sgfn], cg->fgfs[beta->sgfn],
                               cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                               cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                               cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn],
                               cg->fgfs[dquR1->sgfn], cg->fgfs[dquR2->sgfn], cg->fgfs[dquI1->sgfn], cg->fgfs[dquI2->sgfn],
                               cg->fgfs[bdquR1->sgfn], cg->fgfs[bdquR2->sgfn], cg->fgfs[bdquI1->sgfn], cg->fgfs[bdquI2->sgfn],
                               cg->fgfs[dgR->sgfn], cg->fgfs[dgI->sgfn], cg->fgfs[bdgR->sgfn], cg->fgfs[bdgI->sgfn],
                               cg->fgfs[RNews->sgfn], cg->fgfs[INews->sgfn], Rmin, sPp->data->sst, PhysTime);
#endif
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }

  {
    int Varwt[1];
    MyList<var> *DG_List;
    DG_List = new MyList<var>(RNews);
    DG_List->insert(INews);
    Varwt[0] = 2;
    Synch(DG_List, Symmetry, Varwt);
    DG_List->clearList();
  }
}
#if 1
// evolve omega
void NullShellPatch::Check_News(double PhysTime, double dT, bool dp)
{
  MyList<ss_patch> *sPp;

  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        f_get_exact_null(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                         cg->fgfs[RJ0->sgfn], cg->fgfs[IJ0->sgfn], sPp->data->sst, Rmin, PhysTime,
                         cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                         cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                         cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn],
                         cg->fgfs[dquR1->sgfn], cg->fgfs[dquR2->sgfn], cg->fgfs[dquI1->sgfn], cg->fgfs[dquI2->sgfn],
                         cg->fgfs[bdquR1->sgfn], cg->fgfs[bdquR2->sgfn], cg->fgfs[bdquI1->sgfn], cg->fgfs[bdquI2->sgfn],
                         cg->fgfs[dgR->sgfn], cg->fgfs[dgI->sgfn], cg->fgfs[bdgR->sgfn], cg->fgfs[bdgI->sgfn]);

        f_get_null_boundary_c(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                              cg->fgfs[beta->sgfn], cg->fgfs[RQ->sgfn], cg->fgfs[IQ->sgfn],
                              cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn],
                              cg->fgfs[W->sgfn], cg->fgfs[RTheta->sgfn], cg->fgfs[ITheta->sgfn],
                              cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                              cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                              cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn],
                              cg->fgfs[dquR1->sgfn], cg->fgfs[dquR2->sgfn], cg->fgfs[dquI1->sgfn], cg->fgfs[dquI2->sgfn],
                              cg->fgfs[bdquR1->sgfn], cg->fgfs[bdquR2->sgfn], cg->fgfs[bdquI1->sgfn], cg->fgfs[bdquI2->sgfn],
                              cg->fgfs[dgR->sgfn], cg->fgfs[dgI->sgfn], cg->fgfs[bdgR->sgfn], cg->fgfs[bdgI->sgfn],
                              PhysTime, Rmin, sPp->data->sst);

        f_drive_null_news(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                          cg->fgfs[RJ0->sgfn], cg->fgfs[IJ0->sgfn], cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn],
                          cg->fgfs[RTheta->sgfn], cg->fgfs[ITheta->sgfn],
                          cg->fgfs[omega0->sgfn], cg->fgfs[beta->sgfn],
                          cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                          cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                          cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn],
                          cg->fgfs[dquR1->sgfn], cg->fgfs[dquR2->sgfn], cg->fgfs[dquI1->sgfn], cg->fgfs[dquI2->sgfn],
                          cg->fgfs[bdquR1->sgfn], cg->fgfs[bdquR2->sgfn], cg->fgfs[bdquI1->sgfn], cg->fgfs[bdquI2->sgfn],
                          cg->fgfs[dgR->sgfn], cg->fgfs[dgI->sgfn], cg->fgfs[bdgR->sgfn], cg->fgfs[bdgI->sgfn],
                          cg->fgfs[RNews->sgfn], cg->fgfs[INews->sgfn], Rmin, sPp->data->sst);
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }

  {
    int Varwt[1];
    MyList<var> *DG_List;
    DG_List = new MyList<var>(RNews);
    DG_List->insert(INews);
    Varwt[0] = 2;
    Synch(DG_List, Symmetry, Varwt);
    DG_List->clearList();
  }
  // evolve omega
  int iter_count = 0; // count RK4 substeps
  int pre = 0, cor = 1;
  int ERROR = 0;
  double TT = PhysTime;

  // Predictor
  {
    sPp = PatL;
    while (sPp)
    {
      MyList<Block> *BP = sPp->data->blb;
      while (BP)
      {
        Block *cg = BP->data;
        cg->swapList(TheList, JrhsList, myrank);
        if (myrank == cg->rank)
        {
#if 1
          f_get_exact_omegau(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                             cg->fgfs[omega_rhs->sgfn],
                             cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                             cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                             cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn],
                             cg->fgfs[dquR1->sgfn], cg->fgfs[dquR2->sgfn], cg->fgfs[dquI1->sgfn], cg->fgfs[dquI2->sgfn],
                             cg->fgfs[bdquR1->sgfn], cg->fgfs[bdquR2->sgfn], cg->fgfs[bdquI1->sgfn], cg->fgfs[bdquI2->sgfn],
                             cg->fgfs[dgR->sgfn], cg->fgfs[dgI->sgfn], cg->fgfs[bdgR->sgfn], cg->fgfs[bdgI->sgfn],
                             PhysTime, Rmin, sPp->data->sst);
#if 0	     
             f_euler_rout(cg->shape, dT,cg->fgfs[omega0->sgfn],cg->fgfs[omega_rhs->sgfn]);
	     PhysTime += dT;
             f_get_exact_omega(cg->shape,cg->X[0],cg->X[1],cg->X[2],
                                  cg->fgfs[omega->sgfn],sPp->data->sst,Rmin,PhysTime);
	     PhysTime -= dT;
	     if(sPp->data->sst==0 && cg->X[0][0] < -PI/4 && cg->X[1][0] < -PI/4)
	     {
		     int hi=cg->shape[0]/2-1,hj=cg->shape[1]/2-1,hk=cg->shape[2]-1;
		     int hg=hi+hj*cg->shape[0]+hk*cg->shape[0]*cg->shape[1];
		     cout<<cg->fgfs[omega->sgfn][hg]-1<<","<<cg->fgfs[omega0->sgfn][hg]-1<<endl;
	     }
#endif
#else
          f_omega_rhs(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                      cg->fgfs[omega0->sgfn], cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn], cg->fgfs[omega_rhs->sgfn],
                      cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                      cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn]);
#endif
          f_rungekutta4_rout(cg->shape, dT, cg->fgfs[omega0->sgfn], cg->fgfs[omega->sgfn], cg->fgfs[omega_rhs->sgfn],
                             iter_count);
        }
        if (BP == sPp->data->ble)
          break;
        BP = BP->next;
      }
      sPp = sPp->next;
    }
  }

  // corrector
  for (iter_count = 1; iter_count < 4; iter_count++)
  {
    // for RK4: t0, t0+dt/2, t0+dt/2, t0+dt;
    if (iter_count == 1 || iter_count == 3)
      TT += dT / 2;
    {
      sPp = PatL;
      while (sPp)
      {
        MyList<Block> *BP = sPp->data->blb;
        while (BP)
        {
          Block *cg = BP->data;
          cg->swapList(TheList, J1List, myrank);
          if (myrank == cg->rank)
          {
            f_get_exact_null(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                             cg->fgfs[RJ0->sgfn], cg->fgfs[IJ0->sgfn], sPp->data->sst, Rmin, TT,
                             cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                             cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                             cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn],
                             cg->fgfs[dquR1->sgfn], cg->fgfs[dquR2->sgfn], cg->fgfs[dquI1->sgfn], cg->fgfs[dquI2->sgfn],
                             cg->fgfs[bdquR1->sgfn], cg->fgfs[bdquR2->sgfn], cg->fgfs[bdquI1->sgfn], cg->fgfs[bdquI2->sgfn],
                             cg->fgfs[dgR->sgfn], cg->fgfs[dgI->sgfn], cg->fgfs[bdgR->sgfn], cg->fgfs[bdgI->sgfn]);

            f_get_null_boundary_c(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                  cg->fgfs[beta->sgfn], cg->fgfs[RQ->sgfn], cg->fgfs[IQ->sgfn],
                                  cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn],
                                  cg->fgfs[W->sgfn], cg->fgfs[RTheta->sgfn], cg->fgfs[ITheta->sgfn],
                                  cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                                  cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                                  cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn],
                                  cg->fgfs[dquR1->sgfn], cg->fgfs[dquR2->sgfn], cg->fgfs[dquI1->sgfn], cg->fgfs[dquI2->sgfn],
                                  cg->fgfs[bdquR1->sgfn], cg->fgfs[bdquR2->sgfn], cg->fgfs[bdquI1->sgfn], cg->fgfs[bdquI2->sgfn],
                                  cg->fgfs[dgR->sgfn], cg->fgfs[dgI->sgfn], cg->fgfs[bdgR->sgfn], cg->fgfs[bdgI->sgfn],
                                  TT, Rmin, sPp->data->sst);
#if 1
            f_get_exact_omegau(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                               cg->fgfs[omega1->sgfn],
                               cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                               cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                               cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn],
                               cg->fgfs[dquR1->sgfn], cg->fgfs[dquR2->sgfn], cg->fgfs[dquI1->sgfn], cg->fgfs[dquI2->sgfn],
                               cg->fgfs[bdquR1->sgfn], cg->fgfs[bdquR2->sgfn], cg->fgfs[bdquI1->sgfn], cg->fgfs[bdquI2->sgfn],
                               cg->fgfs[dgR->sgfn], cg->fgfs[dgI->sgfn], cg->fgfs[bdgR->sgfn], cg->fgfs[bdgI->sgfn],
                               PhysTime, Rmin, sPp->data->sst);
#else
            f_omega_rhs(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                        cg->fgfs[omega->sgfn], cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn], cg->fgfs[omega1->sgfn],
                        cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                        cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn]);
#endif
            f_rungekutta4_rout(cg->shape, dT, cg->fgfs[omega0->sgfn], cg->fgfs[omega1->sgfn], cg->fgfs[omega_rhs->sgfn],
                               iter_count);
          }
          if (iter_count < 3)
            cg->swapList(SynchList_cor, SynchList_pre, myrank);
          else
          {
            cg->swapList(StateList, SynchList_cor, myrank);
            cg->swapList(OldStateList, SynchList_cor, myrank);
          }
          if (BP == sPp->data->ble)
            break;
          BP = BP->next;
        }
        sPp = sPp->next;
      }
    }

    int Varwt[1];
    MyList<var> *DG_List;
    DG_List = new MyList<var>(omega0);
    DG_List->insert(FXZEO);
    Varwt[0] = 0;
    Synch(DG_List, Symmetry, Varwt);
    DG_List->clearList();
  }
#if 0
    {
      sPp=PatL;
      while(sPp)
      {
        MyList<Block> *BP=sPp->data->blb;
        while(BP)
        {
          Block *cg=BP->data;
          cg->swapList(TheList,J1List,myrank);
  	  if(myrank == cg->rank) 
	  {
	     PhysTime += dT;
             f_get_exact_omega(cg->shape,cg->X[0],cg->X[1],cg->X[2],
                                  cg->fgfs[omega->sgfn],sPp->data->sst,Rmin,PhysTime);
	     PhysTime -= dT;
	     if(sPp->data->sst==0 && cg->X[0][0] < -PI/4 && cg->X[1][0] < -PI/4)
	     {
		     int hi=cg->shape[0]/2-1,hj=cg->shape[1]/2-1,hk=cg->shape[2]-1;
		     int hg=hi+hj*cg->shape[0]+hk*cg->shape[0]*cg->shape[1];
		     cout<<cg->fgfs[omega->sgfn][hg]-1<<","<<cg->fgfs[omega0->sgfn][hg]-1<<endl;
	     }
	  }
          if(BP==sPp->data->ble) break;
          BP=BP->next;
        }
        sPp=sPp->next;
      }
    }
#endif

#if 0
// dump omega for check
{
     MyList<var> * DG_List; 
     DG_List=new MyList<var>(omega0);
     Dump_Data(DG_List,"evo",PhysTime,dT);

    sPp=PatL;
    while(sPp)
    {
      MyList<Block> *BP=sPp->data->blb;
      int fngfs = sPp->data->fngfs;
      while(BP)
      {
        Block *cg=BP->data;
        if(myrank == cg->rank) 
        {
     		f_get_exact_omega(cg->shape,cg->X[0],cg->X[1],cg->X[2],
                                  cg->fgfs[omega0->sgfn],sPp->data->sst,Rmin,TT);
	}
        if(BP==sPp->data->ble) break;
        BP=BP->next;
      }
      sPp=sPp->next;
    }

     Dump_Data(DG_List,"exa",PhysTime,dT);
     DG_List->clearList(); 

     if(TT>0.5 && myrank==0) MPI_Abort(MPI_COMM_WORLD,1);
}
#endif
}
#else
// given omega
void NullShellPatch::Check_News(double PhysTime, double dT, bool dp)
{
  MyList<ss_patch> *sPp;

  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        f_get_exact_null(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                         cg->fgfs[RJ0->sgfn], cg->fgfs[IJ0->sgfn], sPp->data->sst, Rmin, PhysTime,
                         cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                         cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                         cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn],
                         cg->fgfs[dquR1->sgfn], cg->fgfs[dquR2->sgfn], cg->fgfs[dquI1->sgfn], cg->fgfs[dquI2->sgfn],
                         cg->fgfs[bdquR1->sgfn], cg->fgfs[bdquR2->sgfn], cg->fgfs[bdquI1->sgfn], cg->fgfs[bdquI2->sgfn],
                         cg->fgfs[dgR->sgfn], cg->fgfs[dgI->sgfn], cg->fgfs[bdgR->sgfn], cg->fgfs[bdgI->sgfn]);

        f_get_null_boundary_c(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                              cg->fgfs[beta->sgfn], cg->fgfs[RQ->sgfn], cg->fgfs[IQ->sgfn],
                              cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn],
                              cg->fgfs[W->sgfn], cg->fgfs[RTheta->sgfn], cg->fgfs[ITheta->sgfn],
                              cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                              cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                              cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn],
                              cg->fgfs[dquR1->sgfn], cg->fgfs[dquR2->sgfn], cg->fgfs[dquI1->sgfn], cg->fgfs[dquI2->sgfn],
                              cg->fgfs[bdquR1->sgfn], cg->fgfs[bdquR2->sgfn], cg->fgfs[bdquI1->sgfn], cg->fgfs[bdquI2->sgfn],
                              cg->fgfs[dgR->sgfn], cg->fgfs[dgI->sgfn], cg->fgfs[bdgR->sgfn], cg->fgfs[bdgI->sgfn],
                              PhysTime, Rmin, sPp->data->sst);

        f_get_exact_omega(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                          cg->fgfs[omega0->sgfn], sPp->data->sst, Rmin, PhysTime);

        f_drive_null_news(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                          cg->fgfs[RJ0->sgfn], cg->fgfs[IJ0->sgfn], cg->fgfs[RU->sgfn], cg->fgfs[IU->sgfn],
                          cg->fgfs[RTheta->sgfn], cg->fgfs[ITheta->sgfn],
                          cg->fgfs[omega0->sgfn], cg->fgfs[beta->sgfn],
                          cg->fgfs[qlR1->sgfn], cg->fgfs[qlR2->sgfn], cg->fgfs[qlI1->sgfn], cg->fgfs[qlI2->sgfn],
                          cg->fgfs[quR1->sgfn], cg->fgfs[quR2->sgfn], cg->fgfs[quI1->sgfn], cg->fgfs[quI2->sgfn],
                          cg->fgfs[gR->sgfn], cg->fgfs[gI->sgfn],
                          cg->fgfs[dquR1->sgfn], cg->fgfs[dquR2->sgfn], cg->fgfs[dquI1->sgfn], cg->fgfs[dquI2->sgfn],
                          cg->fgfs[bdquR1->sgfn], cg->fgfs[bdquR2->sgfn], cg->fgfs[bdquI1->sgfn], cg->fgfs[bdquI2->sgfn],
                          cg->fgfs[dgR->sgfn], cg->fgfs[dgI->sgfn], cg->fgfs[bdgR->sgfn], cg->fgfs[bdgI->sgfn],
                          cg->fgfs[RNews->sgfn], cg->fgfs[INews->sgfn], Rmin, sPp->data->sst);
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }

  {
    int Varwt[1];
    MyList<var> *DG_List;
    DG_List = new MyList<var>(RNews);
    DG_List->insert(INews);
    Varwt[0] = 2;
    Synch(DG_List, Symmetry, Varwt);
    DG_List->clearList();
  }
}
#endif
double NullShellPatch::News_Error_Check(double PhysTime, double dT, bool dp)
{
  MyList<ss_patch> *sPp;

  double tvf, dtvf = 0;
  int BDW = ghost_width, OBDW = overghost;
  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        f_l2normhelper_sh(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                          sPp->data->bbox[0], sPp->data->bbox[1], sPp->data->bbox[2], sPp->data->bbox[3], sPp->data->bbox[4], sPp->data->bbox[5],
                          cg->fgfs[RNews->sgfn], tvf, BDW, OBDW, Symmetry);
        dtvf += tvf;
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }

  MPI_Allreduce(&dtvf, &tvf, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  tvf = sqrt(tvf);

  return tvf;
}
