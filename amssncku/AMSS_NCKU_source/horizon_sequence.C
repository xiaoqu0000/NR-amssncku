#include <stdio.h>
#include <assert.h>

#include "stdc.h"
#include "util.h"

#include "horizon_sequence.h"

namespace AHFinderDirect
{

	horizon_sequence::horizon_sequence(int N_horizons_in)
		: N_horizons_(N_horizons_in),
		  my_N_horizons_(0), // sequence starts out empty
		  posn_(-1),
		  my_hn_(new int[N_horizons_in])
	{
	}

	horizon_sequence::~horizon_sequence()
	{
		delete[] my_hn_;
	}
	//
	// This function appends  hn  to the sequence.  It returns the new value
	// of my_N_horizons().
	//
	int horizon_sequence::append_hn(int hn)
	{
		assert(hn > 0);						  // can only append genuine horizons
		assert(my_N_horizons_ < N_horizons_); // make sure there's space for it
		my_hn_[my_N_horizons_++] = hn;
		posn_ = 0;
		return my_N_horizons_;
	}

	//******************************************************************************

	//
	// This function computes the internal position immediately following
	// a given internal position in the sequence.
	//
	// Arguments:
	// p = (in) The current internal position, with posn_ semantics
	//
	// Results:
	// This function returns the next internal position after p.
	//
	int horizon_sequence::next_posn(int pos)
		const
	{
		return (pos < 0)					? pos - 1
			   : (pos + 1 < my_N_horizons_) ? pos + 1
											: -1;
	}

	//******************************************************************************

	//
	// This function determines whether or not a given  hn  is genuine.
	//
	bool horizon_sequence::is_hn_genuine(int hn)
		const
	{
		for (int pos = 0; pos < my_N_horizons_; ++pos)
		{
			if (my_hn_[pos] == hn)
				then return true;
		}

		return false;
	}

	//******************************************************************************

} // namespace AHFinderDirect
