

#include "macrodef.h"
#ifdef With_AHF

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>

#include "cctk.h"

#include "config.h"
#include "stdc.h"
#include "util.h"
#include "array.h"
#include "cpm_map.h"
#include "linear_map.h"

#include "coords.h"
#include "tgrid.h"
#include "fd_grid.h"
#include "patch.h"
#include "patch_edge.h"
#include "patch_interp.h"
#include "ghost_zone.h"
#include "patch_system.h"

#include "Jacobian.h"

#include "gfns.h"
#include "gr.h"

#include "horizon_sequence.h"
#include "BH_diagnostics.h"
#include "myglobal.h"

namespace AHFinderDirect
{
	void recentering(patch_system &ps, double max_x, double max_y, double max_z,
					 double min_x, double min_y, double min_z,
					 double centroid_x, double centroid_y, double centroid_z);
	extern struct state state;

	void AHFinderDirect_find_horizons(int HN, int *dumpid,
									  double *xc, double *yc, double *zc, double *xr, double *yr, double *zr,
									  bool *trigger, double *dT)
	{
		const int my_proc = state.my_proc;
		horizon_sequence &hs = *state.my_hs;
		if (my_proc == 0 && hs.N_horizons() != HN)
		{
			cout << "input number " << HN << " != " << "number of wanted horizons " << hs.N_horizons() << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		state.ADM->AH_Prepare_derivatives();

		for (int hn = hs.init_hn(); hs.is_genuine(); hn = hs.next_hn())
		{
			int ihn = hs.get_hn();
			assert(ihn > 0 && ihn <= HN);
			ihn = ihn - 1;

			struct AH_data &AH_data = *state.AH_data_array[hn];

			AH_data.find_trigger = trigger[ihn];
			if (AH_data.find_trigger)
			{
				if (AH_data.found_flag)
					AH_data.initial_find_flag = false;
				else if (AH_data.recentering_flag == false)
				{
					patch_system &ps = *AH_data.ps_ptr;
					recentering(ps, xc[ihn] + xr[ihn] / 2, yc[ihn] + yr[ihn] / 2, zc[ihn] + zr[ihn] / 2,
								xc[ihn] - xr[ihn] / 2, yc[ihn] - yr[ihn] / 2, zc[ihn] - zr[ihn] / 2,
								xc[ihn], yc[ihn], zc[ihn]);
					setup_initial_guess(ps, xc[ihn], yc[ihn], zc[ihn], xr[ihn], yr[ihn], zr[ihn]);
					AH_data.initial_find_flag = true;
				}
				else
					AH_data.stop_finding == true;
			}

		} // end for hn

		Newton(state.N_procs, state.N_active_procs, my_proc,
			   *state.my_hs, state.AH_data_array,
			   state.isb, dumpid, dT);
	}

	void AHFinderDirect_enforcefind(int HN,
									double *xc, double *yc, double *zc, double *xr, double *yr, double *zr)
	{
		const int my_proc = state.my_proc;
		horizon_sequence &hs = *state.my_hs;
		if (my_proc == 0 && hs.N_horizons() != HN)
		{
			cout << "input number " << HN << " != " << "number of wanted horizons " << hs.N_horizons() << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		bool *trigger;
		int *dumpid;
		double *dTT;
		trigger = new bool[HN];
		dumpid = new int[HN];
		dTT = new double[HN];
		for (int ihn = 0; ihn < HN; ihn++)
		{
			trigger[ihn] = true;
			dumpid[ihn] = 1;
			dTT[ihn] = 1;
		}

		for (int hn = hs.init_hn(); hs.is_genuine(); hn = hs.next_hn())
		{
			int ihn = hs.get_hn();
			assert(ihn > 0 && ihn <= HN);

			struct AH_data &AH_data = *state.AH_data_array[hn];

			AH_data.find_trigger = true;
			AH_data.stop_finding = false;
			AH_data.found_flag = false;
			AH_data.recentering_flag = false;
			AH_data.initial_find_flag = true;

		} // end for hn

		AHFinderDirect_find_horizons(HN, dumpid, xc, yc, zc, xr, yr, zr, trigger, dTT);

		delete[] trigger;
		delete[] dumpid;
		delete[] dTT;
	}
} // namespace AHFinderDirect
#endif
