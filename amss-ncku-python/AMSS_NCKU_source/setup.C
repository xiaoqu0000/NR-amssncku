#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include <mpi.h>

#include "util_Table.h"
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
#include "driver.h"
using namespace std;

#include "myglobal.h"
#include "bssn_class.h"

namespace AHFinderDirect
{
	struct state state;

	using jtutil::error_exit;

	namespace
	{
		int allocate_horizons_to_processor(int N_procs, int my_proc,
										   int N_horizons, bool multiproc_flag,
										   horizon_sequence &my_hs)
		{
			const int N_active_procs = multiproc_flag ? Mymin(N_procs, N_horizons)
													  : 1;
			// Implementation note:
			// We allocate the horizons to active processors in round-robin order.
			//
			int proc = 0;
			for (int hn = 1; hn <= N_horizons; ++hn)
			{
				if (proc == my_proc)
					my_hs.append_hn(hn);
				if (++proc >= N_active_procs)
					proc = 0;
			}

			return N_active_procs;
		}
	}

	extern struct state state;

	void AHFinderDirect_setup(MyList<var> *AHList, MyList<var> *GaugeList, bssn_class *ADM,
							  int Symmetry, int HN, double *PhysTime)
	{
		enum patch_system::patch_system_type ps_type;

		switch (Symmetry)
		{
		case 2:
			ps_type = patch_system::patch_system__plus_xyz_octant_mirrored;
			break;
		case 1:
			ps_type = patch_system::patch_system__plus_z_hemisphere;
			break;
		case 0:
			ps_type = patch_system::patch_system__full_sphere;
			break;
		default:
			jtutil::error_exit(ERROR_EXIT, "** Symmetry=%d is not support by AHFD yet.", Symmetry);
		}

		int nprocs = 1, myrank = 0;
		MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
		MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

		state.PhysTime = PhysTime; // Synchonize the PhysTime
		state.Symmetry = Symmetry;
		state.AHList = AHList;
		state.GaugeList = GaugeList;
		state.ADM = ADM;
		state.N_procs = nprocs;
		state.my_proc = myrank;

		state.N_horizons = HN;

		//
		// (genuine) horizon sequence for this processor
		//
		state.my_hs = new horizon_sequence(state.N_horizons);
		horizon_sequence &hs = *state.my_hs;

		const bool multiproc_flag = true;
		state.N_active_procs = allocate_horizons_to_processor(state.N_procs, state.my_proc,
															  state.N_horizons, multiproc_flag,
															  hs);

		// ... horizon numbers run from 1 to N_horizons inclusive
		//     so the array size is N_horizons+1
		state.AH_data_array = new AH_data *[HN + 1];
		for (int hn = 0; hn <= HN; ++hn)
		{
			state.AH_data_array[hn] = NULL;
		}

		int NNP = 0, NNP_out;
		for (int hn = 1; hn <= hs.N_horizons(); ++hn)
		{
			const bool genuine_flag = hs.is_hn_genuine(hn);
			state.AH_data_array[hn] = new AH_data;
			struct AH_data &AH_data = *state.AH_data_array[hn];

			AH_data.recentering_flag = false;
			AH_data.stop_finding = false;

			// create the patch system
			AH_data.ps_ptr = new patch_system(0, 0, 0, // just dummy set, we will recenter it when setting initial guess
											  ps_type, 2, 1,
											  20, 1,
											  //			      (genuine_flag ? 53 : 0),
											  (genuine_flag ? gfns::nominal_max_gfn
															: gfns::skeletal_nominal_max_gfn),
											  -1, -1,
											  1, 1,
											  1, 1,
											  true, false);
			patch_system &ps = *AH_data.ps_ptr;

			if (genuine_flag)
				ps.set_gridfn_to_constant(1.0, gfns::gfn__one);

			AH_data.Jac_ptr = genuine_flag ? new Jacobian(ps) : NULL;

			AH_data.surface_expansion = 0;

			AH_data.initial_find_flag = genuine_flag;

			AH_data.found_flag = false;
			AH_data.BH_diagnostics_fileptr = NULL;

			NNP = Mymax(NNP, AH_data.ps_ptr->N_grid_points());
		} // end of for hn

		MPI_Allreduce(&NNP, &NNP_out, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

		state.Data = new double[NNP_out * 35];
		state.oX = new double[NNP_out];
		state.oY = new double[NNP_out];
		state.oZ = new double[NNP_out];
	}
	void AHFinderDirect_cleanup()
	{
		horizon_sequence &hs = *state.my_hs;
		for (int hn = 1; hn <= hs.N_horizons(); ++hn)
		{
			struct AH_data &AH_data = *state.AH_data_array[hn];
			if (AH_data.ps_ptr)
				delete AH_data.ps_ptr;
			if (AH_data.Jac_ptr)
				delete AH_data.Jac_ptr;
			delete state.AH_data_array[hn];
		} // end of for hn
		delete[] state.AH_data_array;
		delete state.my_hs;
		delete[] state.oX;
		delete[] state.oY;
		delete[] state.oZ;
		delete[] state.Data;
	}
} // namespace AHFinderDirect
