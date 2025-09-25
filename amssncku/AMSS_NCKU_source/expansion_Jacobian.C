

#include "macrodef.h"
#ifdef With_AHF

#include <stdio.h>
#include <assert.h>
#include <math.h>

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

namespace AHFinderDirect
{
	using jtutil::error_exit;

	namespace
	{

		void expansion_Jacobian_partial_SD(patch_system &ps, Jacobian &Jac,
										   bool print_msg_flag);

		void add_ghost_zone_Jacobian(const patch_system &ps,
									 Jacobian &Jac,
									 fp mol,
									 const patch &xp, const ghost_zone &xmgz,
									 int x_II,
									 int xm_irho, int xm_isigma);

		enum expansion_status
		expansion_Jacobian_dr_FD(patch_system *ps_ptr, Jacobian *Jac_ptr, fp add_to_expansion,
								 bool initial_flag,
								 bool print_msg_flag);
	}

	//******************************************************************************

	//
	// If ps_ptr != NULL and Jac_ptr != NULL, this function computes the
	// Jacobian matrix J[Theta(h)] of the expansion Theta(h).  We assume
	// that Theta(h) has already been computed.
	//
	// If ps_ptr == NULL and Jac_ptr == NULL, this function does a dummy
	// computation, in which only any expansion() (and hence geometry
	// interpolator) calls are done, these with the number of interpolation
	// points set to 0 and all the output array pointers set to NULL.
	//
	// It's illegal for one but not both of ps_ptr and Jac_ptr to be NULL.
	//
	// Arguments:
	// ps_ptr --> The patch system, or == NULL to do (only) a dummy computation.
	// Jac_ptr --> The Jacobian, or == NULL to do (only) a dummy computation.
	// add_to_expansion = A real number to add to the expansion.
	//
	// Results:
	// This function returns a status code indicating whether the computation
	// succeeded or failed, and if the latter, what caused the failure.
	//
	enum expansion_status
	expansion_Jacobian(patch_system *ps_ptr, Jacobian *Jac_ptr,
					   fp add_to_expansion,
					   bool initial_flag,
					   bool print_msg_flag /* = false */)
	{
		const bool active_flag = (ps_ptr != NULL) && (Jac_ptr != NULL);
		enum expansion_status status;

		if (active_flag)
			then expansion_Jacobian_partial_SD(*ps_ptr, *Jac_ptr,
											   print_msg_flag);
		// this function looks at ps_ptr and Jac_ptr (non-NULL vs NULL)
		// to choose a normal vs dummy computation
		{
			status = expansion_Jacobian_dr_FD(ps_ptr, Jac_ptr, add_to_expansion,
											  initial_flag,
											  print_msg_flag);
			if (status != expansion_success)
				then return status; // *** ERROR RETURN ***
		}

		return expansion_success; // *** NORMAL RETURN ***
	}
	//
	// This function computes the partial derivative terms in the Jacobian
	// matrix of the expansion Theta(h), by symbolic differentiation from
	// the Jacobian coefficient (angular) gridfns.  The Jacobian is traversed
	// by rows, using equation (25) of my 1996 apparent horizon finding paper.
	//
	// Inputs (angular gridfns, on ghosted grid):
	//	h			# shape of trial surface
	//	Theta			# Theta(h) assumed to already be computed
	//	partial_Theta_wrt_partial_d_h	# Jacobian coefficients
	//	partial_Theta_wrt_partial_dd_h	# (also assumed to already be computed)
	//
	// Outputs:
	//	The Jacobian matrix is stored in the Jacobian object Jac.
	//
	namespace
	{
		void expansion_Jacobian_partial_SD(patch_system &ps, Jacobian &Jac,
										   bool print_msg_flag)
		{
			Jac.zero_matrix();
			ps.compute_synchronize_Jacobian();

			for (int xpn = 0; xpn < ps.N_patches(); ++xpn)
			{
				patch &xp = ps.ith_patch(xpn);

				for (int x_irho = xp.min_irho(); x_irho <= xp.max_irho(); ++x_irho)
				{
					for (int x_isigma = xp.min_isigma(); x_isigma <= xp.max_isigma(); ++x_isigma)
					{
						//
						// compute the main Jacobian terms for this grid point, i.e.
						//	partial Theta(this point x, Jacobian row II)
						//	---------------------------------------------
						//	partial h(other points y, Jacobian column JJ)
						//

						// Jacobian row index
						const int II = ps.gpn_of_patch_irho_isigma(xp, x_irho, x_isigma);

						// Jacobian coefficients for this point
						const fp Jacobian_coeff_rho = xp.gridfn(gfns::gfn__partial_Theta_wrt_partial_d_h_1,
																x_irho, x_isigma);
						const fp Jacobian_coeff_sigma = xp.gridfn(gfns::gfn__partial_Theta_wrt_partial_d_h_2,
																  x_irho, x_isigma);
						const fp Jacobian_coeff_rho_rho = xp.gridfn(gfns::gfn__partial_Theta_wrt_partial_dd_h_11,
																	x_irho, x_isigma);
						const fp Jacobian_coeff_rho_sigma = xp.gridfn(gfns::gfn__partial_Theta_wrt_partial_dd_h_12,
																	  x_irho, x_isigma);
						const fp Jacobian_coeff_sigma_sigma = xp.gridfn(gfns::gfn__partial_Theta_wrt_partial_dd_h_22,
																		x_irho, x_isigma);

						// partial_rho, partial_rho_rho
						{
							for (int m_irho = xp.molecule_min_m();
								 m_irho <= xp.molecule_max_m();
								 ++m_irho)
							{
								const int xm_irho = x_irho + m_irho;
								const fp Jac_rho = Jacobian_coeff_rho * xp.partial_rho_coeff(m_irho);
								const fp Jac_rho_rho = Jacobian_coeff_rho_rho * xp.partial_rho_rho_coeff(m_irho);
								const fp Jac_sum = Jac_rho + Jac_rho_rho;
								if (xp.is_in_nominal_grid(xm_irho, x_isigma))
									then
									{
										const int xm_JJ = Jac.II_of_patch_irho_isigma(xp, xm_irho, x_isigma);
										Jac.sum_into_element(II, xm_JJ, Jac_sum);
									}
								else
									add_ghost_zone_Jacobian(ps, Jac,
															Jac_sum,
															xp, xp.minmax_rho_ghost_zone(m_irho < 0),
															II, xm_irho, x_isigma);
							}
						}

						// partial_sigma, partial_sigma_sigma
						{
							for (int m_isigma = xp.molecule_min_m();
								 m_isigma <= xp.molecule_max_m();
								 ++m_isigma)
							{
								const int xm_isigma = x_isigma + m_isigma;
								const fp Jac_sigma = Jacobian_coeff_sigma * xp.partial_sigma_coeff(m_isigma);
								const fp Jac_sigma_sigma = Jacobian_coeff_sigma_sigma * xp.partial_sigma_sigma_coeff(m_isigma);
								const fp Jac_sum = Jac_sigma + Jac_sigma_sigma;
								if (xp.is_in_nominal_grid(x_irho, xm_isigma))
									then
									{
										const int xm_JJ = Jac.II_of_patch_irho_isigma(xp, x_irho, xm_isigma);
										Jac.sum_into_element(II, xm_JJ, Jac_sum);
									}
								else
									add_ghost_zone_Jacobian(ps, Jac,
															Jac_sum,
															xp, xp.minmax_sigma_ghost_zone(m_isigma < 0),
															II, x_irho, xm_isigma);
							}
						}

						// partial_rho_sigma
						{
							for (int m_irho = xp.molecule_min_m();
								 m_irho <= xp.molecule_max_m();
								 ++m_irho)
							{
								for (int m_isigma = xp.molecule_min_m();
									 m_isigma <= xp.molecule_max_m();
									 ++m_isigma)
								{
									const int xm_irho = x_irho + m_irho;
									const int xm_isigma = x_isigma + m_isigma;
									const fp Jac_rho_sigma = Jacobian_coeff_rho_sigma * xp.partial_rho_sigma_coeff(m_irho, m_isigma);
									if (xp.is_in_nominal_grid(xm_irho, xm_isigma))
										then
										{
											const int xm_JJ = Jac.II_of_patch_irho_isigma(xp, xm_irho, xm_isigma);
											Jac.sum_into_element(II, xm_JJ, Jac_rho_sigma);
										}
									else
									{
										const ghost_zone &xmgz = xp.corner_ghost_zone_containing_point(m_irho < 0, m_isigma < 0,
																									   xm_irho, xm_isigma);
										add_ghost_zone_Jacobian(ps, Jac,
																Jac_rho_sigma,
																xp, xmgz,
																II, xm_irho, xm_isigma);
									}
								}
							}
						}
					}
				}
			}
		}
	}

	//******************************************************************************

	//
	// This function adds the ghost-zone Jacobian dependency contributions
	// for a single ghost-zone point, to a Jacobian matrix.
	//
	// Arguments:
	// ps = The patch system.
	// Jac = (out) The Jacobian matrix.
	// mol = The molecule coefficient.
	// xp = The patch containing the center point of the molecule.
	// xmgz = If the x+m point is in a ghost zone, this must be that ghost zone.
	//	  If the x+m point is not in a ghost zone, this argument is ignored.
	// x_II = The Jacobian row of the x point.
	// xm_(irho,isigma) = The coordinates (in xp) of the x+m point of the molecule.

	namespace
	{
		void add_ghost_zone_Jacobian(const patch_system &ps,
									 Jacobian &Jac,
									 fp mol,
									 const patch &xp, const ghost_zone &xmgz,
									 int x_II,
									 int xm_irho, int xm_isigma)
		{
			const patch_edge &xme = xmgz.my_edge();
			const int xm_iperp = xme.iperp_of_irho_isigma(xm_irho, xm_isigma);
			const int xm_ipar = xme.ipar_of_irho_isigma(xm_irho, xm_isigma);

			// FIXME: this won't change from one call to another
			//        ==> it would be more efficient to reuse the same buffer
			//            across multiple calls on this function
			int global_min_ym, global_max_ym;
			ps.synchronize_Jacobian_global_minmax_ym(global_min_ym, global_max_ym);
			jtutil::array1d<fp> Jacobian_buffer(global_min_ym, global_max_ym);

			// on what other points y does this molecule point xm depend
			// via the patch_system::synchronize() operation?
			int y_iperp;
			int y_posn, min_ym, max_ym;
			const patch_edge &ye = ps.synchronize_Jacobian(xmgz,
														   xm_iperp, xm_ipar,
														   y_iperp,
														   y_posn, min_ym, max_ym,
														   Jacobian_buffer);
			patch &yp = ye.my_patch();

			// add the Jacobian contributions from the ym points
			for (int ym = min_ym; ym <= max_ym; ++ym)
			{
				const int y_ipar = y_posn + ym;
				const int y_irho = ye.irho_of_iperp_ipar(y_iperp, y_ipar);
				const int y_isigma = ye.isigma_of_iperp_ipar(y_iperp, y_ipar);
				const int y_JJ = Jac.II_of_patch_irho_isigma(yp, y_irho, y_isigma);
				Jac.sum_into_element(x_II, y_JJ, mol * Jacobian_buffer(ym));
			}
		}
	}

	//******************************************************************************

	//
	// If ps_ptr != NULL and Jac_ptr != NULL, this function sums the d/dr
	// terms into the Jacobian matrix of the expansion Theta(h), computing
	// those terms by finite differencing.
	//
	// If ps_ptr == NULL and Jac_ptr == NULL, this function does a dummy
	// computation, in which only any expansion() (and hence geometry
	// interpolator) calls are done, these with the number of interpolation
	// points set to 0 and all the output array pointers set to NULL.
	//
	// It's illegal for one but not both of ps_ptr and Jac_ptr to be NULL.
	//
	// The basic algorithm is that
	//	Jac += diag[ (Theta(h+epsilon) - Theta(h)) / epsilon ]
	//
	// Inputs (angular gridfns, on ghosted grid):
	//	h			# shape of trial surface
	//	Theta			# Theta(h) assumed to already be computed
	//
	// Outputs:
	//	Jac += d/dr terms
	//
	// Results:
	// This function returns a status code indicating whether the computation
	// succeeded or failed, and if the latter, what caused the failure.
	//
	namespace
	{
		enum expansion_status
		expansion_Jacobian_dr_FD(patch_system *ps_ptr, Jacobian *Jac_ptr, fp add_to_expansion,
								 bool initial_flag,
								 bool print_msg_flag)
		{
			const bool active_flag = (ps_ptr != NULL) && (Jac_ptr != NULL);

			const double epsilon = 1e-6;
			// compute Theta(h+epsilon)
			if (active_flag)
				then
				{
					ps_ptr->gridfn_copy(gfns::gfn__Theta, gfns::gfn__save_Theta);
					ps_ptr->add_to_ghosted_gridfn(epsilon, gfns::gfn__h);
				}
			const enum expansion_status status = expansion(ps_ptr, add_to_expansion,
														   initial_flag);
			if (status != expansion_success)
				then return status; // *** ERROR RETURN ***

			if (active_flag)
				then
				{
					for (int pn = 0; pn < ps_ptr->N_patches(); ++pn)
					{
						patch &p = ps_ptr->ith_patch(pn);
						for (int irho = p.min_irho(); irho <= p.max_irho(); ++irho)
						{
							for (int isigma = p.min_isigma();
								 isigma <= p.max_isigma();
								 ++isigma)
							{
								const int II = ps_ptr->gpn_of_patch_irho_isigma(p, irho, isigma);
								const fp old_Theta = p.gridfn(gfns::gfn__save_Theta,
															  irho, isigma);
								const fp new_Theta = p.gridfn(gfns::gfn__Theta,
															  irho, isigma);
								const fp d_dr_term = (new_Theta - old_Theta) / epsilon;
								Jac_ptr->sum_into_element(II, II, d_dr_term);
							}
						}
					}

					// restore h and Theta
					ps_ptr->add_to_ghosted_gridfn(-epsilon, gfns::gfn__h);
					ps_ptr->gridfn_copy(gfns::gfn__save_Theta, gfns::gfn__Theta);
				}

			return expansion_success; // *** NORMAL RETURN ***
		}
	}

	//******************************************************************************

} // namespace AHFinderDirect
#endif
