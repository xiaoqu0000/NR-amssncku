//$Id: Newton.C,v 1.1 2012/04/03 10:49:44 zjcao Exp $

#include "macrodef.h"
#ifdef With_AHF

#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <float.h>
#include <math.h>
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
#include "myglobal.h"

namespace AHFinderDirect
{
	extern struct state state;
	using jtutil::error_exit;

	void recentering(patch_system &ps, double max_x, double max_y, double max_z,
					 double min_x, double min_y, double min_z,
					 double centroid_x, double centroid_y, double centroid_z)
	{
		fp ox = ps.origin_x();
		fp oy = ps.origin_y();
		fp oz = ps.origin_z();

		const fp CTR_TOLERENCE = .45;
		bool center = (abs(max_x + min_x - 2.0 * ox) < CTR_TOLERENCE * (max_x - min_x)) &&
					  (abs(max_y + min_y - 2.0 * oy) < CTR_TOLERENCE * (max_y - min_y)) &&
					  (abs(max_z + min_z - 2.0 * oz) < CTR_TOLERENCE * (max_z - min_z));

		if (!center)
		{

			for (int pn = 0; pn < ps.N_patches(); ++pn)
			{
				patch &p = ps.ith_patch(pn);

				for (int irho = p.min_irho(); irho <= p.max_irho(); ++irho)
					for (int isigma = p.min_isigma(); isigma <= p.max_isigma(); ++isigma)
					{

						p.ghosted_gridfn(gfns::gfn__h, irho, isigma) =
							sqrt(jtutil::pow2(p.gridfn(gfns::gfn__global_x, irho, isigma) - centroid_x) +
								 jtutil::pow2(p.gridfn(gfns::gfn__global_y, irho, isigma) - centroid_y) +
								 jtutil::pow2(p.gridfn(gfns::gfn__global_z, irho, isigma) - centroid_z));
					}
			}

			ps.recentering(centroid_x, centroid_y, centroid_z);
		}
	}

	namespace
	{
		bool broadcast_status(int N_procs, int N_active_procs,
							  int my_proc, bool my_active_flag,
							  int hn, int iteration,
							  enum expansion_status expansion_status,
							  fp mean_horizon_radius, fp infinity_norm,
							  bool found_this_horizon, bool I_need_more_iterations,
							  struct iteration_status_buffers &isb);

		void Newton_step(patch_system &ps,
						 fp mean_horizon_radius, fp max_allowable_Delta_h_over_h);

		void save_oldh(patch_system &ps);

		int interpolate_alsh(patch_system *ps_ptr)
		{
			int status = 1;

#define CAST_PTR_OR_NULL(type_, ptr_) \
	(ps_ptr == NULL) ? NULL : static_cast<type_>(ptr_)

			//
			// ***** interpolation points *****
			//
			const int N_interp_points = (ps_ptr == NULL) ? 0 : ps_ptr->N_grid_points();
			double *interp_coords[3] = {
				CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__global_x)),
				CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__global_y)),
				CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__global_z)),
			};

			double *const output_arrays[] = {
				CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__global_xx)), // Lapse-1
				CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__global_xy)), // Sfx
				CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__global_xz)), // Sfy
				CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__global_yy)), // Sfz
			};

			const int N_output_arrays_dim = sizeof(output_arrays) / sizeof(output_arrays[0]);
			const int N_output_arrays_use = N_output_arrays_dim;

			double *Data, *oX, *oY, *oZ;

			int s;
			int Npts = 0;
			for (int ncpu = 0; ncpu < state.N_procs; ncpu++)
			{

				if (state.my_proc == ncpu)
					Npts = N_interp_points;

				MPI_Bcast(&Npts, 1, MPI_INT, ncpu, MPI_COMM_WORLD);

				if (Npts != 0)
				{
					Data = new double[Npts * N_output_arrays_use];

					oX = new double[Npts];
					oY = new double[Npts];
					oZ = new double[Npts];
					if (state.my_proc == ncpu)
					{
						memcpy(oX, interp_coords[0], Npts * sizeof(double));
						memcpy(oY, interp_coords[1], Npts * sizeof(double));
						memcpy(oZ, interp_coords[2], Npts * sizeof(double));
					}
					MPI_Bcast(oX, Npts, MPI_DOUBLE, ncpu, MPI_COMM_WORLD);
					MPI_Bcast(oY, Npts, MPI_DOUBLE, ncpu, MPI_COMM_WORLD);
					MPI_Bcast(oZ, Npts, MPI_DOUBLE, ncpu, MPI_COMM_WORLD);

					// each cpu calls interpolator
					s = globalInterpGFLlash(
						oX, oY, oZ, Npts,
						Data); // 1 succuss; 0 fail

					if (state.my_proc == ncpu)
					{
						status = s;

						if (status == 1)
						{
							for (int ngf = 0; ngf < N_output_arrays_use; ngf++)
							{
								memcpy(output_arrays[ngf], Data + ngf * N_interp_points,
									   sizeof(double) * N_interp_points);
							}
						}
					}

					delete[] oX;
					delete[] oY;
					delete[] oZ;
					delete[] Data;
				}
			}

			return status;
		}

	}

	//******************************************************************************
	void Newton(int N_procs, int N_active_procs, int my_proc,
				horizon_sequence &hs, struct AH_data *const AH_data_array[],
				struct iteration_status_buffers &isb, int *dumpid, double *dT)
	{
		const bool my_active_flag = hs.has_genuine_horizons();
		const int N_horizons = hs.N_horizons();

		for (int hn = hs.init_hn();; hn = hs.next_hn()) // hn always =0 for cpu who has no patch_system
		{
			bool horizon_is_genuine = hs.is_genuine();
			const bool there_is_another_genuine_horizon = hs.is_next_genuine();

			struct AH_data *AH_data_ptr = horizon_is_genuine ? AH_data_array[hn] : NULL;

			horizon_is_genuine = horizon_is_genuine && AH_data_ptr->find_trigger && !AH_data_ptr->stop_finding;
			if (horizon_is_genuine)
				cout << "being finding horizon #" << hn << endl;
			patch_system *const ps_ptr = horizon_is_genuine ? AH_data_ptr->ps_ptr : NULL;
			Jacobian *const Jac_ptr = horizon_is_genuine ? AH_data_ptr->Jac_ptr : NULL;
			const double add_to_expansion = horizon_is_genuine ? -AH_data_ptr->surface_expansion : 0.0;
			const int max_iterations = horizon_is_genuine
										   ? (AH_data_ptr->initial_find_flag ? 80 : 20)
										   : INT_MAX;

			if (horizon_is_genuine)
				save_oldh(*ps_ptr);

			for (int iteration = 1;; ++iteration)
			{
				if (horizon_is_genuine && iteration == max_iterations)
					cout << "AHfinder: fail to find horizon #" << hn
						 << " with Newton iteration " << iteration << " steps!!!" << endl;
				jtutil::norm<fp> Theta_norms;

				const enum expansion_status raw_expansion_status = expansion(ps_ptr, add_to_expansion,
																			 (iteration == 1), true, &Theta_norms);

				const bool Theta_is_ok = (raw_expansion_status == expansion_success);
				const bool norms_are_ok = horizon_is_genuine && Theta_is_ok;

				//
				// have we found this horizon?
				// if so, compute and output BH diagnostics
				//
				const bool found_this_horizon = norms_are_ok && (Theta_norms.infinity_norm() <= 1e-11);

				if (horizon_is_genuine)
					AH_data_ptr->found_flag = found_this_horizon;

				if (horizon_is_genuine && found_this_horizon)
					cout << "found horizon #" << hn << " with " << iteration << " steps!!!" << endl;
				//
				// see if the expansion is too big
				// (if so, we'll give up on this horizon)
				//
				const bool expansion_is_too_large = norms_are_ok && (Theta_norms.infinity_norm() > 1e10);

				//
				// compute the mean horizon radius, and if it's too large,
				// then pretend expansion() returned a "surface too large" error status
				//
				jtutil::norm<fp> h_norms;
				if (horizon_is_genuine)
					then ps_ptr->ghosted_gridfn_norms(gfns::gfn__h, h_norms);
				const fp mean_horizon_radius = horizon_is_genuine ? h_norms.mean()
																  : 0.0;
				const bool horizon_is_too_large = (mean_horizon_radius > 1e10);

				const enum expansion_status effective_expansion_status = horizon_is_too_large ? expansion_failure__surface_too_large
																							  : raw_expansion_status;

				//
				// see if we need more iterations (either on this or another horizon)
				//

				// does *this* horizon need more iterations?
				// i.e. has this horizon's Newton iteration not yet converged?
				const bool this_horizon_needs_more_iterations = horizon_is_genuine && Theta_is_ok && !found_this_horizon && !expansion_is_too_large && !horizon_is_too_large && (iteration < max_iterations);

				// do I (this processor) need to do more iterations
				// on this or a following horizon?
				const bool I_need_more_iterations = this_horizon_needs_more_iterations || there_is_another_genuine_horizon;

				//
				// broadcast iteration status from each active processor
				// to all processors, and inclusive-or the "we need more iterations"
				// flags to see if *any* (active) processor needs more iterations
				//
				const bool any_proc_needs_more_iterations = broadcast_status(N_procs, N_active_procs,
																			 my_proc, my_active_flag,
																			 hn, iteration, effective_expansion_status,
																			 mean_horizon_radius,
																			 (norms_are_ok ? Theta_norms.infinity_norm() : 0.0),
																			 found_this_horizon, I_need_more_iterations,
																			 isb);
				// set found-this-horizon flags
				// for all active processors' non-dummy horizons
				for (int found_proc = 0; found_proc < N_active_procs; ++found_proc)
				{
					const int found_hn = isb.hn_buffer[found_proc];
					if (found_hn > 0)
						AH_data_array[found_hn]->found_flag = isb.found_horizon_buffer[found_proc];
				}

				//
				// prepare lapse and shift
				{
					int ff = 0, fft = 0;
					if (found_this_horizon && dumpid[hn - 1] > 0 && dT[hn - 1] > 0)
						fft = 1;
					MPI_Allreduce(&fft, &ff, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

					if (ff)
					{
						if ((interpolate_alsh(ps_ptr) == 0) && (state.my_proc == 0))
							cout << "interpolation of lapse and shift for AH failed." << endl;
					}
				}

				if (found_this_horizon)
				{
					struct BH_diagnostics &BH_diagnostics = AH_data_ptr->BH_diagnostics;
					// output data
					if (dumpid[hn - 1] > 0)
					{
						char filename[100];
						sprintf(filename, "ah%02d_%05d.dat", hn, dumpid[hn - 1]);
						if (dT[hn - 1] > 0)
						{
							// gridfunction xx,xy,xz,yy,yz,zz will be used as temp storage
							BH_diagnostics.compute_signature(*ps_ptr, dT[hn - 1]);
							ps_ptr->print_gridfn_with_xyz(gfns::gfn__global_zz, true, gfns::gfn__h, filename);
						}
						else
							ps_ptr->print_ghosted_gridfn_with_xyz(gfns::gfn__h, true, gfns::gfn__h, filename, false);
					}

					BH_diagnostics.compute(*ps_ptr); // gridfunction xx,xy,xz,yy,yz,zz changed

					if (AH_data_ptr->BH_diagnostics_fileptr == NULL)
						AH_data_ptr->BH_diagnostics_fileptr = BH_diagnostics.setup_output_file(N_horizons, hn);
					BH_diagnostics.output(AH_data_ptr->BH_diagnostics_fileptr, (*state.PhysTime));

					// recentering
					recentering(*ps_ptr, (AH_data_ptr->BH_diagnostics).max_x, (AH_data_ptr->BH_diagnostics).max_y, (AH_data_ptr->BH_diagnostics).max_z,
								(AH_data_ptr->BH_diagnostics).min_x, (AH_data_ptr->BH_diagnostics).min_y, (AH_data_ptr->BH_diagnostics).min_z,
								(AH_data_ptr->BH_diagnostics).centroid_x, (AH_data_ptr->BH_diagnostics).centroid_y, (AH_data_ptr->BH_diagnostics).centroid_z);
					AH_data_ptr->recentering_flag = true;
				}

				//
				// are all processors done with all their genuine horizons?
				// or if this is a single-processor run, are we done with this horizon?
				//
				if (!any_proc_needs_more_iterations)
					return; // *** NORMAL RETURN ***

				//
				// compute the Jacobian matrix
				// *** this is a synchronous operation across all processors ***
				//

				const enum expansion_status
					Jacobian_status = expansion_Jacobian(this_horizon_needs_more_iterations ? ps_ptr : NULL,
														 this_horizon_needs_more_iterations ? Jac_ptr : NULL,
														 add_to_expansion,
														 (iteration == 1),
														 false);
				const bool Jacobian_is_ok = (Jacobian_status == expansion_success);

				//
				// skip to the next horizon unless
				// this is a genuine Jacobian computation, and it went ok
				//
				if (!(this_horizon_needs_more_iterations && Jacobian_is_ok))
					break; // *** LOOP EXIT ***

				//
				// compute the Newton step
				//
				Jac_ptr->solve_linear_system(gfns::gfn__Theta, gfns::gfn__Delta_h, false);

				Newton_step(*ps_ptr, mean_horizon_radius, 0.1);

				// end of this Newton iteration
			}

			// end of this horizon
		}

		// we should never get to here
		assert(false);
	}

	//******************************************************************************
	//******************************************************************************
	//******************************************************************************
	namespace
	{
		bool broadcast_status(int N_procs, int N_active_procs,
							  int my_proc, bool my_active_flag,
							  int hn, int iteration,
							  enum expansion_status effective_expansion_status,
							  fp mean_horizon_radius, fp infinity_norm,
							  bool found_this_horizon, bool I_need_more_iterations,
							  struct iteration_status_buffers &isb)
		{
			assert(my_proc >= 0);
			assert(my_proc < N_procs);

			enum
			{
				buffer_var__hn = 0,	   // also encodes found_this_horizon flag
									   // in sign: +=true, -=false
				buffer_var__iteration, // also encodes I_need_more_iterations flag
									   // in sign: +=true, -=false
				buffer_var__expansion_status,
				buffer_var__mean_horizon_radius,
				buffer_var__Theta_infinity_norm,
				N_buffer_vars // no comma
			};

			//
			// allocate buffers if this is the first use
			//
			if (isb.hn_buffer == NULL)
				then
				{
					isb.hn_buffer = new int[N_active_procs];
					isb.iteration_buffer = new int[N_active_procs];
					isb.expansion_status_buffer = new enum expansion_status[N_active_procs];
					isb.mean_horizon_radius_buffer = new fp[N_active_procs];
					isb.Theta_infinity_norm_buffer = new fp[N_active_procs];
					isb.found_horizon_buffer = new bool[N_active_procs];

					isb.send_buffer_ptr = new jtutil::array2d<double>(0, N_active_procs - 1,
																	  0, N_buffer_vars - 1);
					isb.receive_buffer_ptr = new jtutil::array2d<double>(0, N_active_procs - 1,
																		 0, N_buffer_vars - 1);
				}
			jtutil::array2d<double> &send_buffer = *isb.send_buffer_ptr;
			jtutil::array2d<double> &receive_buffer = *isb.receive_buffer_ptr;

			//
			// pack this processor's values into the reduction buffer
			//
			jtutil::zero_C_array(send_buffer.N_array(), send_buffer.data_array());
			if (my_active_flag)
				then
				{
					assert(send_buffer.is_valid_i(my_proc));
					assert(hn >= 0);	   // encoding scheme assumes this
					assert(iteration > 0); // encoding scheme assumes this
					send_buffer(my_proc, buffer_var__hn) = found_this_horizon ? +hn : -hn;
					send_buffer(my_proc, buffer_var__iteration) = I_need_more_iterations ? +iteration : -iteration;
					send_buffer(my_proc, buffer_var__expansion_status) = int(effective_expansion_status);
					send_buffer(my_proc, buffer_var__mean_horizon_radius) = mean_horizon_radius;
					send_buffer(my_proc, buffer_var__Theta_infinity_norm) = infinity_norm;
				}

			const int reduction_status = MPI_Allreduce(static_cast<void *>(send_buffer.data_array()),
													   static_cast<void *>(receive_buffer.data_array()),
													   send_buffer.N_array(),
													   MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);

			// if (reduction_status < 0)
			if (reduction_status != MPI_SUCCESS)
				then CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
								"broadcast_status(): error status %d from reduction!",
								reduction_status); /*NOTREACHED*/

			//
			// unpack the reduction buffer back to the high-level result buffers and
			// compute the inclusive-or of the broadcast I_need_more_iterations flags
			//
			bool any_proc_needs_more_iterations = false;
			for (int proc = 0; proc < N_active_procs; ++proc)
			{
				const int hn_temp = static_cast<int>(
					receive_buffer(proc, buffer_var__hn));
				isb.hn_buffer[proc] = jtutil::abs(hn_temp);
				isb.found_horizon_buffer[proc] = (hn_temp > 0);

				const int iteration_temp = static_cast<int>(
					receive_buffer(proc, buffer_var__iteration));
				isb.iteration_buffer[proc] = jtutil::abs(iteration_temp);
				const bool proc_needs_more_iterations = (iteration_temp > 0);
				any_proc_needs_more_iterations |= proc_needs_more_iterations;

				isb.expansion_status_buffer[proc] = static_cast<enum expansion_status>(
					static_cast<int>(
						receive_buffer(proc, buffer_var__expansion_status)));

				isb.mean_horizon_radius_buffer[proc] = receive_buffer(proc, buffer_var__mean_horizon_radius);
				isb.Theta_infinity_norm_buffer[proc] = receive_buffer(proc, buffer_var__Theta_infinity_norm);
			}

			return any_proc_needs_more_iterations;
		}
	}
	//
	// This function takes the Newton step, scaling it down if it's too large.
	//
	// Arguments:
	// ps = The patch system containing the gridfns h and Delta_h.
	// mean_horizon_radius = ||h||_mean
	// max_allowable_Delta_h_over_h = The maximum allowable
	//				     ||Delta_h||_infinity / ||h||_mean
	//				  Any step over this is internally clamped
	//				  (scaled down) to this size.
	//
	namespace
	{
		void Newton_step(patch_system &ps,
						 fp mean_horizon_radius, fp max_allowable_Delta_h_over_h)
		{
			//
			// compute scale factor (1 for small steps, <1 for large steps)
			//

			const fp max_allowable_Delta_h = max_allowable_Delta_h_over_h * mean_horizon_radius;

			jtutil::norm<fp> Delta_h_norms;
			ps.gridfn_norms(gfns::gfn__Delta_h, Delta_h_norms);
			const fp max_Delta_h = Delta_h_norms.infinity_norm();

			const fp scale = (max_Delta_h <= max_allowable_Delta_h)
								 ? 1.0
								 : max_allowable_Delta_h / max_Delta_h;

			//
			// take the Newton step (scaled if necessary)
			//
			for (int pn = 0; pn < ps.N_patches(); ++pn)
			{
				patch &p = ps.ith_patch(pn);

				for (int irho = p.min_irho(); irho <= p.max_irho(); ++irho)
				{
					for (int isigma = p.min_isigma();
						 isigma <= p.max_isigma();
						 ++isigma)
					{
						p.ghosted_gridfn(gfns::gfn__h, irho, isigma) -= scale * p.gridfn(gfns::gfn__Delta_h, irho, isigma);
					}
				}
			}
		}
		void save_oldh(patch_system &ps)
		{
			for (int pn = 0; pn < ps.N_patches(); ++pn)
			{
				patch &p = ps.ith_patch(pn);

				for (int irho = p.min_irho(); irho <= p.max_irho(); ++irho)
				{
					for (int isigma = p.min_isigma();
						 isigma <= p.max_isigma();
						 ++isigma)
					{
						p.gridfn(gfns::gfn__oldh, irho, isigma) = p.ghosted_gridfn(gfns::gfn__h, irho, isigma);
					}
				}
			}
		}
	}

	//******************************************************************************

} // namespace AHFinderDirect
#endif
