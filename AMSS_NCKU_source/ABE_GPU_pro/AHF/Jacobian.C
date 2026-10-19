#include <stdlib.h>
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
#include "ilucg.h"
// all the code in this file is inside this namespace
namespace AHFinderDirect
{
	// this represents a single element stored in the matrix for
	// sort_row_into_column_order()  and  sort_row_into_column_order__cmp()
	struct matrix_element
	{
		int JA;
		fp A;
	};

	Jacobian::Jacobian(patch_system &ps)
		: ps_(ps),
		  N_rows_(ps.N_grid_points()),
		  N_nonzeros_(0), current_N_rows_(0), N_nonzeros_allocated_(0),
		  IA_(new integer[N_rows_ + 1]), JA_(NULL), A_(NULL),
		  itemp_(NULL), rtemp_(NULL)
	{
		IO_ = 1;
		zero_matrix();
	}

	Jacobian::~Jacobian()
	{
		if (A_)
			delete[] A_;
		if (JA_)
			delete[] JA_;
		if (IA_)
			delete[] IA_;
		if (rtemp_)
			delete[] rtemp_;
		if (itemp_)
			delete[] itemp_;
	}

	double Jacobian::element(int II, int JJ)
		const
	{
		const int posn = find_element(II, JJ);
		return (posn >= 0) ? A_[posn] : 0.0;
	}

	void Jacobian::zero_matrix()
	{

		N_nonzeros_ = 0;
		current_N_rows_ = 0;
		IA_[0] = IO_;
	}

	void Jacobian::set_element(int II, int JJ, fp value)
	{
		const int posn = find_element(II, JJ);
		if (posn >= 0)
			then A_[posn] = value;
		else
			insert_element(II, JJ, value);
	}

	void Jacobian::sum_into_element(int II, int JJ, fp value)
	{
		const int posn = find_element(II, JJ);
		if (posn >= 0)
			then A_[posn] += value;
		else
			insert_element(II, JJ, value);
	}

	int Jacobian::find_element(int II, int JJ)
		const
	{
		if (II >= current_N_rows_)
			then return -1; // this row not defined yet

		const int start = IA_[II] - IO_;
		const int stop = IA_[II + 1] - IO_;
		for (int posn = start; posn < stop; ++posn)
		{
			if (JA_[posn] - IO_ == JJ)
				then return posn; // found
		}

		return -1; // not found
	}

	int Jacobian::insert_element(int II, int JJ, double value)
	{
		if (!((II == current_N_rows_ - 1) || (II == current_N_rows_)))
		{
			printf(
				"***** row_sparse_Jacobian::insert_element(II=%d, JJ=%d, value=%g):\n"
				"        attempt to insert element elsewhere than {last row, last row+1}!\n"
				"        N_rows_=%d   current_N_rows_=%d   IO_=%d\n"
				"        N_nonzeros_=%d   N_nonzeros_allocated_=%d\n",
				II, JJ, double(value),
				N_rows_, current_N_rows_, IO_,
				N_nonzeros_, N_nonzeros_allocated_);
			abort();
		}

		// start a new row if necessary
		if (II == current_N_rows_)
			then
			{
				assert(current_N_rows_ < N_rows_);
				IA_[current_N_rows_ + 1] = IA_[current_N_rows_];
				++current_N_rows_;
			}

		// insert into current row
		assert(II == current_N_rows_ - 1);
		if (IA_[II + 1] - IO_ >= N_nonzeros_allocated_)
			then grow_arrays();
		const int posn = IA_[II + 1] - IO_;
		assert(posn < N_nonzeros_allocated_);
		JA_[posn] = JJ + IO_;
		A_[posn] = value;
		++IA_[II + 1];
		++N_nonzeros_;

		return posn;
	}

	void Jacobian::grow_arrays()
	{
		N_nonzeros_allocated_ += base_growth_amount + (N_nonzeros_allocated_ >> 1);

		int *const new_JA = new int[N_nonzeros_allocated_];
		double *const new_A = new double[N_nonzeros_allocated_];
		for (int posn = 0; posn < N_nonzeros_; ++posn)
		{
			new_JA[posn] = JA_[posn];
			new_A[posn] = A_[posn];
		}
		delete[] A_;
		delete[] JA_;
		JA_ = new_JA;
		A_ = new_A;
	}

	int compare_matrix_elements(const void *x, const void *y)
	{
		const struct matrix_element *const px = static_cast<const struct matrix_element *>(x);
		const struct matrix_element *const py = static_cast<const struct matrix_element *>(y);

		return px->JA - py->JA;
	}

	void Jacobian::sort_each_row_into_column_order()
	{
		// buffer must be big enough to hold the largest row
		int max_N_in_row = 0;
		{
			for (int II = 0; II < N_rows_; ++II)
			{
				max_N_in_row = max(max_N_in_row, IA_[II + 1] - IA_[II]);
			}
		}

		// contiguous buffer for sorting
		struct matrix_element *const buffer = new struct matrix_element[max_N_in_row];

		{
			for (int II = 0; II < N_rows_; ++II)
			{
				const int N_in_row = IA_[II + 1] - IA_[II];

				// copy this row's JA_[] and A_[] values to the buffer
				const int start = IA_[II] - IO_;
				for (int p = 0; p < N_in_row; ++p)
				{
					const int posn = start + p;
					buffer[p].JA = JA_[posn];
					buffer[p].A = A_[posn];
				}

				// sort the buffer
				qsort(static_cast<void *>(buffer), N_in_row, sizeof(buffer[0]),
					  &compare_matrix_elements);

				// copy the buffer values back to this row's JA_[] and A_[]
				for (int p = 0; p < N_in_row; ++p)
				{
					const int posn = start + p;
					JA_[posn] = buffer[p].JA;
					A_[posn] = buffer[p].A;
				}
			}
		}

		delete[] buffer;
	}

	double Jacobian::solve_linear_system(int rhs_gfn, int x_gfn, bool print_msg_flag)
	{
		assert(IO_ == Fortran_index_origin);
		assert(current_N_rows_ == N_rows_);

		if (itemp_ == NULL)
			then
			{
				itemp_ = new int[3 * N_rows_ + 3 * N_nonzeros_ + 2];
				rtemp_ = new double[4 * N_rows_ + N_nonzeros_];
			}

		// initial guess = all zeros
		double *x = ps_.gridfn_data(x_gfn);
		for (int II = 0; II < N_rows_; ++II)
		{
			x[II] = 0.0;
		}

		const int N = N_rows_;
		const double *rhs = ps_.gridfn_data(rhs_gfn);
		const double eps = 1e-10;
		const int max_iterations = N_rows_;
		int istatus;

		// the actual linear solution
		f_ilucg(N,
				IA_, JA_, A_,
				rhs, x,
				itemp_, rtemp_,
				eps, max_iterations,
				istatus);

		if (istatus < 0)
		{
			printf(
				"***** row_sparse_Jacobian__ILUCG::solve_linear_system(rhs_gfn=%d, x_gfn=%d):\n"
				"        error return from [sd]ilucg() routine!\n"
				"        istatus=%d < 0 ==> bad matrix structure, eg. zero diagonal element!\n",
				rhs_gfn, x_gfn,
				int(istatus));
			abort();
		}

		return -1.0;
	}

} // namespace AHFinderDirect
