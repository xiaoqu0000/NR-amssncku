// Z4c_gpu_class.C
//
// GPU-aware overrides for Z4c_class. Mirrors bssn_gpu_class.C structurally
// (the parent class is bssn_class compiled in its GPU variant). We only
// override the pieces that differ from BSSN:
//   * Step(int,int) dispatches to Step_GPU (defined in Z4c_step_gpu.C)
//   * Constraint_Out() outputs 8 ConV entries (BSSN's seven plus TZ0 L2-norm)
//
// Compute_Constraint / Interp_Constraint / Check_extrop are inherited from
// bssn_class. The ConstraintList includes TZ0, so it gets MPI-synced and
// (for Interp_Constraint) interpolated alongside the other constraints
// automatically; the BSSN-only GPU constraint kernel does not write to TZ0
// but TZ0 already holds the Theta value evolved by the Z4C RHS.

#ifdef newc
#include <sstream>
#include <cstdio>
#else
#include <stdio.h>
#endif

#include "macrodef.h"
#include "misc.h"
#include "Parallel.h"
#include "Z4c_class.h"
#include "bssn_rhs.h"
#include "Z4c_rhs.h"
#include "fmisc.h"
#include "mpi_error_batch.h"
#include "monitor.h"

// ---------------------------------------------------------------------------
// Step(int,int) — dispatch CPU/GPU path.
// Under USE_GPU build use_gpu==1 is set by Get_runtime_envirment(); we
// therefore always run Step_GPU. If use_gpu is somehow 0 we abort because
// the CPU body of Z4c_class::Step is unavailable in a GPU build (it lives
// under #ifndef USE_GPU in Z4c_class.C).
// ---------------------------------------------------------------------------
void Z4c_class::Step(int lev, int YN)
{
  if (use_gpu == 1) {
    Step_GPU(lev, YN);
    return;
  }
  if (myrank == 0) {
    cout << "Z4c_class::Step called with use_gpu=0 in a USE_GPU build -- "
         << "CPU Z4C step is not compiled in. Aborting." << endl;
  }
  MPI_Abort(MPI_COMM_WORLD, 1);
}

// ---------------------------------------------------------------------------
// Constraint_Out() — Z4C variant outputting 8 entries.
//
//   ConV[0..6] : Hamiltonian, momentum (Px,Py,Pz), Gamma (Gx,Gy,Gz)
//   ConV[7]    : TZ0 L2-norm (Z4C-specific Theta constraint monitor)
//
// Structure mirrors bssn_class::Constraint_Out (the GPU variant compiled in
// bssn_gpu_class.C). The bulk recomputation uses the Z4C RHS kernel so that
// TZ-dependent constraint contributions on refined levels are consistent.
// ---------------------------------------------------------------------------
void Z4c_class::Constraint_Out()
{
  LastConsOut += dT * pow(0.5, Mymax(0, trfls));

  if (LastConsOut < AnasTime)
    return;

  // Recompute least the constraint data lost for moved new grid (lev > 0).
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
            gpu_apply_symmetry_bd_fused_launch(
              cg->stream, cg->shape[0], cg->shape[1], cg->shape[2],
              cg->shape_pad[0], cg->shape_pad[1], cg->shape_pad[2], Symmetry,
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
            gpu_apply_symmetry_bd(
              cg->stream,
              cg->shape[0], cg->shape[1], cg->shape[2],
              cg->shape_pad[0], cg->shape_pad[1], cg->shape_pad[2],
              cg->d_fgfs[TZ0->sgfn],
              TZ0->SoA[0], TZ0->SoA[1], TZ0->SoA[2]
            );

            gpu_compute_rhs_z4c_launch(
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
              cg->d_fgfs[TZ0->sgfn],
              cg->d_fgfs[phi_rhs->sgfn], cg->d_fgfs[trK_rhs->sgfn],
              cg->d_fgfs[gxx_rhs->sgfn], cg->d_fgfs[gxy_rhs->sgfn], cg->d_fgfs[gxz_rhs->sgfn],
              cg->d_fgfs[gyy_rhs->sgfn], cg->d_fgfs[gyz_rhs->sgfn], cg->d_fgfs[gzz_rhs->sgfn],
              cg->d_fgfs[Axx_rhs->sgfn], cg->d_fgfs[Axy_rhs->sgfn], cg->d_fgfs[Axz_rhs->sgfn],
              cg->d_fgfs[Ayy_rhs->sgfn], cg->d_fgfs[Ayz_rhs->sgfn], cg->d_fgfs[Azz_rhs->sgfn],
              cg->d_fgfs[Gmx_rhs->sgfn], cg->d_fgfs[Gmy_rhs->sgfn], cg->d_fgfs[Gmz_rhs->sgfn],
              cg->d_fgfs[Lap_rhs->sgfn],
              cg->d_fgfs[Sfx_rhs->sgfn], cg->d_fgfs[Sfy_rhs->sgfn], cg->d_fgfs[Sfz_rhs->sgfn],
              cg->d_fgfs[dtSfx_rhs->sgfn], cg->d_fgfs[dtSfy_rhs->sgfn], cg->d_fgfs[dtSfz_rhs->sgfn],
              cg->d_fgfs[TZ_rhs->sgfn],
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
          }
          if (BP == Pp->data->ble)
            break;
          BP = BP->next;
        }
        Pp = Pp->next;
      }
    }
    Parallel::Sync_GPU(GH->PatL[lev], ConstraintList, Symmetry);
  }

  double ConV[8];
  for (int levi = 0; levi < GH->levels; levi++)
  {
    ConV[0] = Parallel::L2Norm_GPU(GH->PatL[levi]->data, Cons_Ham);
    ConV[1] = Parallel::L2Norm_GPU(GH->PatL[levi]->data, Cons_Px);
    ConV[2] = Parallel::L2Norm_GPU(GH->PatL[levi]->data, Cons_Py);
    ConV[3] = Parallel::L2Norm_GPU(GH->PatL[levi]->data, Cons_Pz);
    ConV[4] = Parallel::L2Norm_GPU(GH->PatL[levi]->data, Cons_Gx);
    ConV[5] = Parallel::L2Norm_GPU(GH->PatL[levi]->data, Cons_Gy);
    ConV[6] = Parallel::L2Norm_GPU(GH->PatL[levi]->data, Cons_Gz);
    ConV[7] = Parallel::L2Norm_GPU(GH->PatL[levi]->data, TZ0);
    ConVMonitor->writefile(PhysTime, 8, ConV);
  }

  Interp_Constraint(false);

  LastConsOut = 0;
}
