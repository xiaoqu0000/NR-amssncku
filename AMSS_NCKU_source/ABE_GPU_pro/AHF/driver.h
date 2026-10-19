#ifndef DRIVER_H
#define DRIVER_H
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>

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

namespace AHFinderDirect
{
	struct iteration_status_buffers
	{
		int *hn_buffer;
		int *iteration_buffer;
		enum expansion_status *expansion_status_buffer;
		fp *mean_horizon_radius_buffer;
		fp *Theta_infinity_norm_buffer;
		bool *found_horizon_buffer;

		jtutil::array2d<CCTK_REAL> *send_buffer_ptr;
		jtutil::array2d<CCTK_REAL> *receive_buffer_ptr;

		iteration_status_buffers()
			: hn_buffer(NULL), iteration_buffer(NULL),
			  expansion_status_buffer(NULL),
			  mean_horizon_radius_buffer(NULL),
			  Theta_infinity_norm_buffer(NULL),
			  found_horizon_buffer(NULL),
			  send_buffer_ptr(NULL), receive_buffer_ptr(NULL)
		{
		}
	};

	//
	// This struct holds interprocessor-communication buffers for broadcasting
	// the BH diagnostics and horizon shape from the processor which finds a
	// given horizon, to all processors.
	//
	struct horizon_buffers
	{
		int N_buffer;
		double *send_buffer;
		double *receive_buffer;

		horizon_buffers()
			: N_buffer(0),
			  send_buffer(NULL),
			  receive_buffer(NULL)
		{
		}
	};
	//
	struct AH_data
	{
		patch_system *ps_ptr;
		Jacobian *Jac_ptr;
		double surface_expansion;

		bool initial_find_flag;
		bool recentering_flag, stop_finding, find_trigger;

		bool found_flag; // did we find this horizon (successfully)

		struct BH_diagnostics BH_diagnostics;
		FILE *BH_diagnostics_fileptr;

		// interprocessor-communication buffers
		// for this horizon's BH diagnostics and (optionally) horizon shape
		struct horizon_buffers horizon_buffers;
	};

	// initial_guess.cc
	void setup_initial_guess(patch_system &ps,
							 fp x_center, fp y_center, fp z_center,
							 fp x_radius, fp y_radius, fp z_radius);

	// Newton.cc
	void Newton(int N_procs, int N_active_procs, int my_proc,
				horizon_sequence &hs, struct AH_data *const AH_data_array[],
				struct iteration_status_buffers &isb, int *dumpid, double *);

} // namespace AHFinderDirect
#endif /*     DRIVER_H    */
