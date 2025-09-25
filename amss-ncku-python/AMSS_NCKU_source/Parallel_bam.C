
#include "Parallel.h"
#include "fmisc.h"
#include "prolongrestrict.h"
#include "misc.h"

void Parallel::OutBdLow2Hi_bam(MyList<Patch> *PLc, MyList<Patch> *PLf,
                               MyList<var> *VarList1 /* source */, MyList<var> *VarList2 /* target */,
                               int Symmetry)
{
   MyList<Parallel::pointstru_bam> *bdsul;
   Constr_pointstr_OutBdLow2Hi(PLf, PLc, bdsul);

   intertransfer(bdsul, VarList1, VarList2, Symmetry);

   destroypsuList_bam(bdsul);
}
void Parallel::Restrict_bam(MyList<Patch> *PLc, MyList<Patch> *PLf,
                            MyList<var> *VarList1 /* source */, MyList<var> *VarList2 /* target */,
                            int Symmetry)
{
   MyList<Parallel::pointstru_bam> *rsul;
   Constr_pointstr_Restrict(PLf, PLc, rsul);

   intertransfer(rsul, VarList1, VarList2, Symmetry);

   destroypsuList_bam(rsul);
}
void Parallel::OutBdLow2Hi_bam(MyList<Patch> *PLc, MyList<Patch> *PLf,
                               MyList<var> *VarList1 /* source */, MyList<var> *VarList2 /* target */,
                               MyList<Parallel::pointstru_bam> *bdsul, int Symmetry)
{
   intertransfer(bdsul, VarList1, VarList2, Symmetry);
}
void Parallel::Restrict_bam(MyList<Patch> *PLc, MyList<Patch> *PLf,
                            MyList<var> *VarList1 /* source */, MyList<var> *VarList2 /* target */,
                            MyList<Parallel::pointstru_bam> *rsul, int Symmetry)
{
   intertransfer(rsul, VarList1, VarList2, Symmetry);
}
void Parallel::Constr_pointstr_OutBdLow2Hi(MyList<Patch> *PLf, MyList<Patch> *PLc,
                                           MyList<Parallel::pointstru_bam> *&bdsul)
{
   MyList<Patch> *PL;

   MyList<Parallel::pointstru_bam> *ps;
   bdsul = 0;

   // find out points
   PL = PLf;
   while (PL)
   {
      double dx, dy, dz;

      dx = PL->data->blb->data->getdX(0);
      dy = PL->data->blb->data->getdX(1);
      dz = PL->data->blb->data->getdX(2);

      double uub[3], llb[3];

      llb[0] = PL->data->bbox[0] + PL->data->lli[0] * dx;
      llb[1] = PL->data->bbox[1] + PL->data->lli[1] * dy;
      llb[2] = PL->data->bbox[2] + PL->data->lli[2] * dz;
      uub[0] = PL->data->bbox[3] - PL->data->uui[0] * dx;
      uub[1] = PL->data->bbox[4] - PL->data->uui[1] * dy;
      uub[2] = PL->data->bbox[5] - PL->data->uui[2] * dz;

      double x, y, z;

      for (int i = 0; i < PL->data->shape[0]; i++)
      {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
         x = PL->data->bbox[0] + i * dx;
#else
#ifdef Cell
         x = PL->data->bbox[0] + (0.5 + i) * dx;
#else
#error Not define Vertex nor Cell
#endif
#endif
         for (int j = 0; j < PL->data->shape[1]; j++)
         {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
            y = PL->data->bbox[1] + j * dy;
#else
#ifdef Cell
            y = PL->data->bbox[1] + (0.5 + j) * dy;
#else
#error Not define Vertex nor Cell
#endif
#endif
            for (int k = 0; k < PL->data->shape[2]; k++)
            {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
               z = PL->data->bbox[2] + k * dz;
#else
#ifdef Cell
               z = PL->data->bbox[2] + (0.5 + k) * dz;
#else
#error Not define Vertex nor Cell
#endif
#endif
               if (!(llb[0] - TINY < x && uub[0] + TINY > x &&
                     llb[1] - TINY < y && uub[1] + TINY > y &&
                     llb[2] - TINY < z && uub[2] + TINY > z)) // not in the inner part
               {
                  if (bdsul)
                  {
                     ps->next = new MyList<Parallel::pointstru_bam>;
                     ps = ps->next;
                     ps->data = new Parallel::pointstru_bam;
                  }
                  else
                  {
                     bdsul = ps = new MyList<Parallel::pointstru_bam>;
                     ps->data = new Parallel::pointstru_bam;
                  }

                  ps->data->pox[0] = x;
                  ps->data->pox[1] = y;
                  ps->data->pox[2] = z;
                  ps->data->Bgs = 0;
                  ps->data->Bgd = 0;
                  ps->data->coef = 0;

                  ps->next = 0;
               }
            }
         }
      }

      PL = PL->next;
   }

   // find out blocks
   ps = bdsul;
   while (ps)
   {
      double x, y, z;
      x = ps->data->pox[0];
      y = ps->data->pox[1];
      z = ps->data->pox[2];
      bool flag;
      // find target block
      flag = true;
      PL = PLf;
      while (flag && PL)
      {
         MyList<Block> *BP = PL->data->blb;
         while (flag && BP)
         {
            double llb[3], uub[3];

            for (int i = 0; i < dim; i++)
            {
               double DH = BP->data->getdX(i);
               uub[i] = (feq(BP->data->bbox[dim + i], PL->data->bbox[dim + i], DH / 2)) ? BP->data->bbox[dim + i] : BP->data->bbox[dim + i] - ghost_width * DH;
               llb[i] = (feq(BP->data->bbox[i], PL->data->bbox[i], DH / 2)) ? BP->data->bbox[i] : BP->data->bbox[i] + ghost_width * DH;
            }

            if (llb[0] - TINY < x && uub[0] + TINY > x &&
                llb[1] - TINY < y && uub[1] + TINY > y &&
                llb[2] - TINY < z && uub[2] + TINY > z)
            {
               ps->data->Bgd = BP->data;
               flag = false;
            }

            if (BP == PL->data->ble)
               break;
            BP = BP->next;
         }
         PL = PL->next;
      }
      if (flag)
      {
         cout << "error in Parallel::Constr_pointstr_OutBdLow2Hi 2" << endl;
         MPI_Abort(MPI_COMM_WORLD, 1);
      }
      // find source block
      flag = true;
      PL = PLc;
      while (flag && PL)
      {
         MyList<Block> *BP = PL->data->blb;
         while (flag && BP)
         {
            double llb[3], uub[3];

            for (int i = 0; i < dim; i++)
            {
               double DH = BP->data->getdX(i);
               uub[i] = (feq(BP->data->bbox[dim + i], PL->data->bbox[dim + i], DH / 2)) ? BP->data->bbox[dim + i] : BP->data->bbox[dim + i] - ghost_width * DH;
               llb[i] = (feq(BP->data->bbox[i], PL->data->bbox[i], DH / 2)) ? BP->data->bbox[i] : BP->data->bbox[i] + ghost_width * DH;
            }

            if (llb[0] - TINY < x && uub[0] + TINY > x &&
                llb[1] - TINY < y && uub[1] + TINY > y &&
                llb[2] - TINY < z && uub[2] + TINY > z)
            {
               ps->data->Bgs = BP->data;
               flag = false;
            }

            if (BP == PL->data->ble)
               break;
            BP = BP->next;
         }
         PL = PL->next;
      }
      if (flag)
      {
         cout << "error in Parallel::Constr_pointstr_OutBdLow2Hi 3" << endl;
         MPI_Abort(MPI_COMM_WORLD, 1);
      }

      ps = ps->next;
   }
}
void Parallel::Constr_pointstr_Restrict(MyList<Patch> *PLf, MyList<Patch> *PLc,
                                        MyList<Parallel::pointstru_bam> *&rsul)
{
   MyList<Parallel::gridseg> *gdlf = 0, *gs;
   MyList<Patch> *PL = PLf;
   while (PL)
   {
      if (gdlf)
      {
         gs->next = new MyList<Parallel::gridseg>;
         gs = gs->next;
         gs->data = new Parallel::gridseg;
      }
      else
      {
         gdlf = gs = new MyList<Parallel::gridseg>;
         gs->data = new Parallel::gridseg;
      }

      gs->next = 0;

      for (int i = 0; i < dim; i++)
      {
         double DH = PL->data->blb->data->getdX(i);

         gs->data->llb[i] = PL->data->bbox[i] + PL->data->lli[i] * DH;
         gs->data->uub[i] = PL->data->bbox[dim + i] - PL->data->uui[i] * DH;
      }

      PL = PL->next;
   }

   MyList<Parallel::pointstru_bam> *ps;
   rsul = 0;

   // find out points
   gs = gdlf;
   while (gs)
   {
      PL = PLc;
      bool flag = true;
      while (flag)
      {
         if (!PL)
         {
            int myrank;
            MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
            if (myrank == 0)
            {
               cout << "error in Parallel::Constr_pointstr_Restrict: fail to find grid segment [" << gs->data->llb[0] << ":" << gs->data->uub[0] << ","
                    << gs->data->llb[1] << ":" << gs->data->uub[1] << ","
                    << gs->data->llb[2] << ":" << gs->data->uub[2] << "]"
                    << endl;
               PL = PLc;
               while (PL)
               {
                  PL->data->checkPatch(0);
                  PL = PL->next;
               }
            }

            misc::tillherecheck("for wait.");
            MPI_Abort(MPI_COMM_WORLD, 1);
         }
         if (gs->data->llb[0] > PL->data->bbox[0] - TINY && gs->data->uub[0] < PL->data->bbox[3] + TINY &&
             gs->data->llb[1] > PL->data->bbox[1] - TINY && gs->data->uub[1] < PL->data->bbox[4] + TINY &&
             gs->data->llb[2] > PL->data->bbox[2] - TINY && gs->data->uub[2] < PL->data->bbox[5] + TINY)
            flag = false;

         if (flag)
            PL = PL->next;
      }

      double dx, dy, dz;

      dx = PL->data->blb->data->getdX(0);
      dy = PL->data->blb->data->getdX(1);
      dz = PL->data->blb->data->getdX(2);

      double x, y, z;

      for (int i = 0; i < PL->data->shape[0]; i++)
      {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
         x = PL->data->bbox[0] + i * dx;
#else
#ifdef Cell
         x = PL->data->bbox[0] + (0.5 + i) * dx;
#else
#error Not define Vertex nor Cell
#endif
#endif
         for (int j = 0; j < PL->data->shape[1]; j++)
         {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
            y = PL->data->bbox[1] + j * dy;
#else
#ifdef Cell
            y = PL->data->bbox[1] + (0.5 + j) * dy;
#else
#error Not define Vertex nor Cell
#endif
#endif
            for (int k = 0; k < PL->data->shape[2]; k++)
            {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
               z = PL->data->bbox[2] + k * dz;
#else
#ifdef Cell
               z = PL->data->bbox[2] + (0.5 + k) * dz;
#else
#error Not define Vertex nor Cell
#endif
#endif
               if (gs->data->llb[0] - TINY < x && gs->data->uub[0] + TINY > x &&
                   gs->data->llb[1] - TINY < y && gs->data->uub[1] + TINY > y &&
                   gs->data->llb[2] - TINY < z && gs->data->uub[2] + TINY > z) // in the inner part
               {
                  if (rsul)
                  {
                     ps->next = new MyList<Parallel::pointstru_bam>;
                     ps = ps->next;
                     ps->data = new Parallel::pointstru_bam;
                  }
                  else
                  {
                     rsul = ps = new MyList<Parallel::pointstru_bam>;
                     ps->data = new Parallel::pointstru_bam;
                  }

                  ps->data->pox[0] = x;
                  ps->data->pox[1] = y;
                  ps->data->pox[2] = z;
                  ps->data->Bgs = 0;
                  ps->data->Bgd = 0;
                  ps->data->coef = 0;

                  ps->next = 0;
               }
            }
         }
      }

      gs = gs->next;
   }

   gdlf->destroyList();

   // find out blocks
   ps = rsul;
   while (ps)
   {
      double x, y, z;
      x = ps->data->pox[0];
      y = ps->data->pox[1];
      z = ps->data->pox[2];
      bool flag;
      // find source block
      flag = true;
      PL = PLf;
      while (flag && PL)
      {
         MyList<Block> *BP = PL->data->blb;
         while (flag && BP)
         {
            double llb[3], uub[3];

            for (int i = 0; i < dim; i++)
            {
               double DH = BP->data->getdX(i);
               uub[i] = (feq(BP->data->bbox[dim + i], PL->data->bbox[dim + i], DH / 2)) ? BP->data->bbox[dim + i] : BP->data->bbox[dim + i] - ghost_width * DH;
               llb[i] = (feq(BP->data->bbox[i], PL->data->bbox[i], DH / 2)) ? BP->data->bbox[i] : BP->data->bbox[i] + ghost_width * DH;
            }

            if (llb[0] - TINY < x && uub[0] + TINY > x &&
                llb[1] - TINY < y && uub[1] + TINY > y &&
                llb[2] - TINY < z && uub[2] + TINY > z)
            {
               ps->data->Bgs = BP->data;
               flag = false;
            }

            if (BP == PL->data->ble)
               break;
            BP = BP->next;
         }
         PL = PL->next;
      }
      if (flag)
      {
         cout << "error in Parallel::Constr_pointstr_Restrict 2" << endl;
         MPI_Abort(MPI_COMM_WORLD, 1);
      }
      // find target block
      flag = true;
      PL = PLc;
      while (flag && PL)
      {
         MyList<Block> *BP = PL->data->blb;
         while (flag && BP)
         {
            double llb[3], uub[3];

            for (int i = 0; i < dim; i++)
            {
               double DH = BP->data->getdX(i);
               uub[i] = (feq(BP->data->bbox[dim + i], PL->data->bbox[dim + i], DH / 2)) ? BP->data->bbox[dim + i] : BP->data->bbox[dim + i] - ghost_width * DH;
               llb[i] = (feq(BP->data->bbox[i], PL->data->bbox[i], DH / 2)) ? BP->data->bbox[i] : BP->data->bbox[i] + ghost_width * DH;
            }

            if (llb[0] - TINY < x && uub[0] + TINY > x &&
                llb[1] - TINY < y && uub[1] + TINY > y &&
                llb[2] - TINY < z && uub[2] + TINY > z)
            {
               ps->data->Bgd = BP->data;
               flag = false;
            }

            if (BP == PL->data->ble)
               break;
            BP = BP->next;
         }
         PL = PL->next;
      }
      if (flag)
      {
         cout << "error in Parallel::Constr_pointstr_Restrict 3" << endl;
         MPI_Abort(MPI_COMM_WORLD, 1);
      }

      ps = ps->next;
   }
}

void Parallel::intertransfer(MyList<Parallel::pointstru_bam> *&sul,
                             MyList<var> *VarList1 /* source */, MyList<var> *VarList2 /*target */,
                             int Symmetry)
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
         // myrank: local; node : remote
         if (length = interdata_packer(0, sul, myrank, node, PACK, VarList1, VarList2, Symmetry))
         {
            rec_data[node] = new double[length];
            if (!rec_data[node])
            {
               cout << "Parallel::intertransfer: out of memory when new in short transfer, place 1" << endl;
               MPI_Abort(MPI_COMM_WORLD, 1);
            }
            interdata_packer(rec_data[node], sul, myrank, node, PACK, VarList1, VarList2, Symmetry);
         }
      }
      else
      {
         // send from this cpu to cpu#node
         if (length = interdata_packer(0, sul, myrank, node, PACK, VarList1, VarList2, Symmetry))
         {
            send_data[node] = new double[length];
            if (!send_data[node])
            {
               cout << "Parallel::intertransfer: out of memory when new in short transfer, place 2" << endl;
               MPI_Abort(MPI_COMM_WORLD, 1);
            }
            interdata_packer(send_data[node], sul, myrank, node, PACK, VarList1, VarList2, Symmetry);
            MPI_Isend((void *)send_data[node], length, MPI_DOUBLE, node, 1, MPI_COMM_WORLD, reqs + req_no++);
         }
         // receive from cpu#node to this cpu
         if (length = interdata_packer(0, sul, myrank, node, UNPACK, VarList1, VarList2, Symmetry))
         {
            rec_data[node] = new double[length];
            if (!rec_data[node])
            {
               cout << "Parallel::intertransfer: out of memory when new in short transfer, place 3" << endl;
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
         interdata_packer(rec_data[node], sul, myrank, node, UNPACK, VarList1, VarList2, Symmetry);

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
int Parallel::interdata_packer(double *data, MyList<Parallel::pointstru_bam> *sul, int myrank, int node, int dir,
                               MyList<var> *VarLists /* source */, MyList<var> *VarListd /* target */, int Symmetry)
{
   int DIM = dim;
   int ordn = 2 * ghost_width;

   if (dir != PACK && dir != UNPACK)
   {
      cout << "Parallel::interdata_packer: error dir " << dir << " for data_packer " << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
   }

   int size_out = 0;

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

   while (sul)
   {
      if ((dir == PACK && sul->data->Bgs->rank == myrank && sul->data->Bgd->rank == node) ||
          (dir == UNPACK && sul->data->Bgd->rank == myrank && sul->data->Bgs->rank == node))
      {
         varls = VarLists;
         varld = VarListd;
         while (varls && varld)
         {
            if (data)
            {
               if (dir == PACK)
               {
                  //         f_global_interp(sul->data->Bgs->shape,sul->data->Bgs->X[0],sul->data->Bgs->X[1],sul->data->Bgs->X[2],
                  //			 sul->data->Bgs->fgfs[varls->data->sgfn],data[size_out],
                  //			 sul->data->pox[0],sul->data->pox[1],sul->data->pox[2],ordn,varls->data->SoA,Symmetry);
                  if (sul->data->coef == 0)
                  {
                     sul->data->coef = new double[ordn * dim];
                     for (int i = 0; i < dim; i++)
                     {
                        double dd = sul->data->Bgs->getdX(i);
                        sul->data->sind[i] = int((sul->data->pox[i] - sul->data->Bgs->X[i][0]) / dd) - ordn / 2 + 1;
                        double h1, h2;
                        for (int j = 0; j < ordn; j++)
                        {
                           h1 = sul->data->Bgs->X[i][0] + (sul->data->sind[i] + j) * dd;
                           sul->data->coef[i * ordn + j] = 1;
                           for (int k = 0; k < j; k++)
                           {
                              h2 = sul->data->Bgs->X[i][0] + (sul->data->sind[i] + k) * dd;
                              sul->data->coef[i * ordn + j] *= (sul->data->pox[i] - h2) / (h1 - h2);
                           }
                           for (int k = j + 1; k < ordn; k++)
                           {
                              h2 = sul->data->Bgs->X[i][0] + (sul->data->sind[i] + k) * dd;
                              sul->data->coef[i * ordn + j] *= (sul->data->pox[i] - h2) / (h1 - h2);
                           }
                        }
                     }
                  }
                  int sst = -1;
                  f_global_interpind(sul->data->Bgs->shape, sul->data->Bgs->X[0], sul->data->Bgs->X[1], sul->data->Bgs->X[2],
                                     sul->data->Bgs->fgfs[varls->data->sgfn], data[size_out],
                                     sul->data->pox[0], sul->data->pox[1], sul->data->pox[2], ordn, varls->data->SoA, Symmetry,
                                     sul->data->sind, sul->data->coef, sst);
               }
               if (dir == UNPACK) // from target data to corresponding grid
                  f_pointcopy(DIM, sul->data->Bgd->bbox, sul->data->Bgd->bbox + dim, sul->data->Bgd->shape, sul->data->Bgd->fgfs[varld->data->sgfn],
                              sul->data->pox[0], sul->data->pox[1], sul->data->pox[2], data[size_out]);
            }
            size_out += 1;
            varls = varls->next;
            varld = varld->next;
         }
      }
      sul = sul->next;
   }

   return size_out;
}
void Parallel::destroypsuList_bam(MyList<pointstru_bam> *ct)
{
   MyList<pointstru_bam> *n;
   while (ct)
   {
      n = ct->next;
      if (ct->data->coef)
         delete[] ct->data->coef;
      delete ct->data;
      delete ct;
      ct = n;
   }
}
