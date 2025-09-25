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

namespace AHFinderDirect
{
	int lagrange_interp(double coor_orin, double dx, double *gf,
						int PTS, double ipx, double *out, int *mposn, double *Jac,
						int ORD) // ORD-1 order lagrange interpolation
	{
		assert(PTS >= ORD);
		int mi, mf;

		double *L, *x;
		L = new double[PTS];
		x = new double[PTS];
		int i, j, k;

		//-- Determine molecular range
		//   for odd points, say 5, the molecular is
		//             |
		//   +-----+---x-+-----+-----+
		//
		mi = jtutil::round<double>::ceiling((ipx - coor_orin) / dx) - ORD / 2;
		mf = mi + ORD;
		if (mi < 0)
		{
			mi = 0;
			mf = ORD;
		}
		else if (mf > PTS)
		{
			mf = PTS;
			mi = PTS - ORD;
		}

		//-- Setup coordinate by input origin, dx
		for (j = mi; j < mf; j++)
			x[j] = coor_orin + j * dx;

		//-- Lagrange basis function
		*out = 0;
		for (i = mi; i < mf; i++)
		{
			L[i] = 1.0;
			for (k = mi; k < mf; k++)
				if (k != i)
				{
					L[i] *= (ipx - x[k]) / (x[i] - x[k]);
				}
			*out += *(gf + i) * L[i];
			*Jac = L[i];
			Jac++;
		}

		*mposn = mi;

		delete[] L;
		delete[] x;

		return 0; // Normal retrun
	}

	using jtutil::error_exit;

	patch_interp::patch_interp(const patch_edge &my_edge_in,
							   int min_iperp_in, int max_iperp_in,
							   const jtutil::array1d<int> &min_parindex_array_in,
							   const jtutil::array1d<int> &max_parindex_array_in,
							   const jtutil::array2d<fp> &interp_par_in,
							   bool ok_to_use_min_par_ghost_zone,
							   bool ok_to_use_max_par_ghost_zone,
							   int interp_handle_in, int interp_par_table_handle_in)
		: my_patch_(my_edge_in.my_patch()),
		  my_edge_(my_edge_in),
		  min_gfn_(my_patch().ghosted_min_gfn()),
		  max_gfn_(my_patch().ghosted_max_gfn()),
		  ok_to_use_min_par_ghost_zone_(ok_to_use_min_par_ghost_zone),
		  ok_to_use_max_par_ghost_zone_(ok_to_use_max_par_ghost_zone),
		  min_iperp_(min_iperp_in), max_iperp_(max_iperp_in),
		  min_ipar_(ok_to_use_min_par_ghost_zone
						? my_edge_in.min_ipar_with_corners()
						: my_edge_in.min_ipar_without_corners()),
		  max_ipar_(ok_to_use_max_par_ghost_zone
						? my_edge_in.max_ipar_with_corners()
						: my_edge_in.max_ipar_without_corners()),
		  min_parindex_array_(min_parindex_array_in),
		  max_parindex_array_(max_parindex_array_in),
		  interp_par_(interp_par_in),
		  interp_handle_(interp_handle_in),
		  interp_par_table_handle_(1),
		  gridfn_coord_origin_(my_edge().par_map().fp_of_int(min_ipar_)),
		  gridfn_coord_delta_(my_edge().par_map().delta_fp()),
		  gridfn_data_ptrs_(min_gfn_, max_gfn_),
		  interp_data_buffer_ptrs_(min_gfn_, max_gfn_) // no comma
	{
		int status;

		const CCTK_INT stride = my_edge().ghosted_par_stride();

		status = 0;
		if (status < 0)
			then error_exit(ERROR_EXIT,
							"***** patch_interp::patch_interp():\n"
							"        can't set gridfn stride in interpolator parmameter table!\n"
							"        error status=%d\n",
							status); /*NOTREACHED*/
	}

	patch_interp::~patch_interp()
	{
	}

	void patch_interp::interpolate(int ghosted_min_gfn_to_interp,
								   int ghosted_max_gfn_to_interp,
								   jtutil::array3d<fp> &data_buffer,
								   jtutil::array2d<CCTK_INT> &posn_buffer,
								   jtutil::array3d<fp> &Jacobian_buffer)
		const

	{
		int status;

		const int N_dims = 1;
		const int N_gridfns = jtutil::how_many_in_range(ghosted_min_gfn_to_interp,
														ghosted_max_gfn_to_interp);
		const CCTK_INT N_gridfn_data_points = jtutil::how_many_in_range(min_ipar(), max_ipar());

		//--  Jacobian
		const int Jacobian_interp_point_stride = Jacobian_buffer.subscript_stride_j();

		//
		// do the interpolations at each iperp
		//
		for (int iperp = min_iperp(); iperp <= max_iperp(); ++iperp)
		{
			//
			// interpolation-point coordinates
			//
			const int min_parindex = min_parindex_array_(iperp);
			const int max_parindex = max_parindex_array_(iperp);
			const CCTK_INT N_interp_points = jtutil::how_many_in_range(min_parindex, max_parindex);
			const fp *const interp_coords_ptr = &interp_par_(iperp, min_parindex);
			const void *const interp_coords[N_dims] = {static_cast<const void *>(interp_coords_ptr)};

			//
			// pointers to gridfn data to interpolate, and to result buffer
			//
			for (int ghosted_gfn = ghosted_min_gfn_to_interp;
				 ghosted_gfn <= ghosted_max_gfn_to_interp;
				 ++ghosted_gfn)
			{
				// set up data pointer to --> (iperp,min_ipar) gridfn
				const int start_irho = my_edge().irho_of_iperp_ipar(iperp, min_ipar());
				const int start_isigma = my_edge().isigma_of_iperp_ipar(iperp, min_ipar());
				gridfn_data_ptrs_(ghosted_gfn) = static_cast<const void *>(
					&my_patch()
						 .ghosted_gridfn(ghosted_gfn,
										 start_irho, start_isigma));
				interp_data_buffer_ptrs_(ghosted_gfn) = static_cast<void *>(
					&data_buffer(ghosted_gfn, iperp, min_parindex));
			}
			const void *const *const gridfn_data = &gridfn_data_ptrs_(ghosted_min_gfn_to_interp);
			void *const *const interp_buffer = &interp_data_buffer_ptrs_(ghosted_min_gfn_to_interp);

			//--  molecule position
			CCTK_POINTER molecule_posn_ptrs[N_dims] = {static_cast<CCTK_POINTER>(&posn_buffer(iperp, min_parindex))};
			//--  Jacobian
			CCTK_POINTER const Jacobian_ptrs[1] //[N_gridfns]
				= {static_cast<CCTK_POINTER>(
					&Jacobian_buffer(iperp, min_parindex, 0))};
			// Jacobian_buffer has continuous memory allocation.

			const CCTK_INT stride = my_edge().ghosted_par_stride();
			double y[N_gridfn_data_points];

			for (int i = 0; i < N_gridfn_data_points; i++)
			{
				y[i] = *((double *)(*gridfn_data) + stride * i);
			}

			const int ORD = 6;
			double Jac[ORD];
			int posn; // of molecular, starting from 0
			for (int i = 0; i < N_interp_points; i++)
			{
				status = lagrange_interp(gridfn_coord_origin_, gridfn_coord_delta_,
										 y, N_gridfn_data_points,
										 *((double *)interp_coords[0] + i), ((double *)(*interp_buffer) + i),
										 &posn, Jac, ORD);

				*((int *)molecule_posn_ptrs[0] + i) = posn + 2;

				memcpy((double *)(Jacobian_ptrs[0]) + Jacobian_buffer.min_k() +
						   Jacobian_interp_point_stride * i,
					   Jac, sizeof(Jac));
			}

			// convert the molecule positions from  parindex-min_ipar
			// to  parindex  values (again, cf comments on array subscripting
			// at the start of "patch_interp.hh")
			for (int parindex = min_parindex;
				 parindex <= max_parindex;
				 ++parindex)
			{
				posn_buffer(iperp, parindex) += min_ipar();
			}

			if (status < 0)
				then error_exit(ERROR_EXIT,
								"***** patch_interp::interpolate():\n"
								"        error return %d from interpolator at iperp=%d of [%d,%d]!\n"
								"        my_patch()=\"%s\" my_edge()=\"%s\"\n",
								status, iperp, min_iperp(), max_iperp(),
								my_patch().name(), my_edge().name()); /*NOTREACHED*/

		} // end for iperp
	}

	void patch_interp::verify_Jacobian_sparsity_pattern_ok()
		const
	{
		CCTK_INT MSS_is_fn_of_interp_coords = 0, MSS_is_fn_of_input_array_values = 0;
		CCTK_INT Jacobian_is_fn_of_input_array_values = 0;

		//
		// verify that we grok the Jacobian sparsity pattern
		//
		if (MSS_is_fn_of_interp_coords || MSS_is_fn_of_input_array_values || Jacobian_is_fn_of_input_array_values)
			then error_exit(ERROR_EXIT,
							"***** patch_interp::verify_Jacobian_sparsity_pattern_ok():\n"
							"        implementation restriction: we only grok Jacobians with\n"
							"        fixed-sized hypercube-shaped molecules, independent of\n"
							"        the interpolation coordinates and the floating-point values!\n"
							"        MSS_is_fn_of_interp_coords=(int)%d (we only grok 0)\n"
							"        MSS_is_fn_of_input_array_values=(int)%d (we only grok 0)\n"
							"        Jacobian_is_fn_of_input_array_values=(int)%d (we only grok 0)\n",
							MSS_is_fn_of_interp_coords,
							MSS_is_fn_of_input_array_values,
							Jacobian_is_fn_of_input_array_values);
	}

	//******************************************************************************

	//
	// This function queries the interpolator to get the [min,max] ipar m
	// coordinates of the interpolation molecules.
	//
	// (This API implicitly assumes that the Jacobian sparsity is one which
	// is "ok" as verified by  verify_Jacobian_sparsity_pattern_ok() .)
	//
	void patch_interp::molecule_minmax_ipar_m(int &min_ipar_m, int &max_ipar_m)
		const
	{
		min_ipar_m = -2;
		max_ipar_m = 3;
	}

	//******************************************************************************

	//
	// This function queries the interpolator at each iperp to find out the
	// molecule ipar positions (which we implicitly assume to be independent
	// of ghosted_gfn), and stores these in  posn_buffer(iperp, parindex) .
	//
	// (This API implicitly assumes that the Jacobian sparsity is one which
	// is "ok" as verified by  verify_Jacobian_sparsity_pattern_ok() .)
	//
	void patch_interp::molecule_posn(jtutil::array2d<CCTK_INT> &posn_buffer)
		const
	{
		const int N_dims = 1;
		int status;

		for (int iperp = min_iperp(); iperp <= max_iperp(); ++iperp)
		{
			const int min_parindex = min_parindex_array_(iperp);
			const int max_parindex = max_parindex_array_(iperp);

			// set up the molecule-position query in the parameter table
			CCTK_POINTER molecule_posn_ptrs[N_dims] = {static_cast<CCTK_POINTER>(&posn_buffer(iperp, min_parindex))};
			status = 0; // Util_TableSetPointerArray(interp_par_table_handle_, N_dims,
						//               molecule_posn_ptrs, "molecule_positions");

			if (status < 0)
				then error_exit(ERROR_EXIT,
								"***** patch_interp::molecule_posn():\n"
								"        can't set molecule position query\n"
								"        in interpolator parmameter table at iperp=%d of [%d,%d]!\n"
								"        error status=%d\n",
								iperp, min_iperp(), max_iperp(),
								status); /*NOTREACHED*/

			for (int parindex = min_parindex;
				 parindex <= max_parindex;
				 ++parindex)
			{
				posn_buffer(iperp, parindex) += min_ipar();
			}
		}
	}

	void patch_interp::Jacobian(jtutil::array3d<fp> &Jacobian_buffer)
		const
	{
		const int N_dims = 1;
		const int N_gridfns = 1;

		int status1, status2;

		//
		// set Jacobian stride info in parameter table
		//
		const int Jacobian_interp_point_stride = Jacobian_buffer.subscript_stride_j();

		status1 = 0;

		status2 = 0;

		if ((status1 < 0) || (status2 < 0))
			then error_exit(ERROR_EXIT,
							"***** patch_interp::Jacobian():\n"
							"        can't set Jacobian stride info in interpolator parmameter table!\n"
							"        error status1=%d status2=%d\n",
							status1, status2);

		//
		// query the Jacobians at each iperp
		//
		for (int iperp = min_iperp(); iperp <= max_iperp(); ++iperp)
		{
			const int min_parindex = min_parindex_array_(iperp);
			const int max_parindex = max_parindex_array_(iperp);

			//
			// set up the Jacobian query in the parameter table
			//
			CCTK_POINTER const Jacobian_ptrs[N_gridfns] = {static_cast<CCTK_POINTER>(
				&Jacobian_buffer(iperp, min_parindex, 0))};
		}
	}
} // namespace AHFinderDirect
