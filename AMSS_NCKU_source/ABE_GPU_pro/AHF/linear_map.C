#include <assert.h>
#include <stdio.h>

#include "stdc.h"
#include "util.h"
#include "linear_map.h"

namespace AHFinderDirect
{
	namespace jtutil
	{
		template <typename fp_t>
		linear_map<fp_t>::linear_map(int min_int_in, int max_int_in,
									 fp_t min_fp_in, fp_t delta_fp_in, fp_t max_fp_in)
			: delta_(delta_fp_in), inverse_delta_(1.0 / delta_fp_in),
			  min_int_(min_int_in), max_int_(max_int_in)
		{
			constructor_common(min_fp_in, max_fp_in);
		}

		template <typename fp_t>
		linear_map<fp_t>::linear_map(const linear_map<fp_t> &lm_in,
									 int min_int_in, int max_int_in) // subrange
			: delta_(lm_in.delta_fp()), inverse_delta_(lm_in.inverse_delta_fp()),
			  min_int_(min_int_in), max_int_(max_int_in)
		{
			if (!(is_in_range(min_int_in) && is_in_range(max_int_in)))
				then error_exit(ERROR_EXIT,
								"***** linear_map<fp_t>::linear_map:\n"
								"        min_int_in=%d and/or max_int_in=%d\n"
								"        aren't in integer range [%d,%d] of existing linear_map!\n",
								min_int_, max_int_,
								lm_in.min_int(), lm_in.max_int()); /*NOTREACHED*/

			constructor_common(lm_in.fp_of_int_unchecked(min_int_in),
							   lm_in.fp_of_int_unchecked(max_int_in));
		}

		//******************************************************************************

		//
		// This function does the common argument validation and setup for
		// all the constructors of class  linear_map<fp_t>:: .
		//
		template <typename fp_t>
		void linear_map<fp_t>::constructor_common(fp_t min_fp_in, fp_t max_fp_in)
		// assumes
		//	min_int_, max_int_, delta_, inverse_delta_
		// are already initialized
		// ==> ok to use min_int(), max_int(), delta_fp(), inverse_delta_fp()
		// ... other class members *not* yet initialized
		{
			origin_ = 0.0; // temp value
			origin_ = min_fp_in - fp_of_int_unchecked(min_int());

			// this should be guaranteed by the above calculation
			assert(fuzzy<fp_t>::EQ(fp_of_int_unchecked(min_int()), min_fp_in));

			// this is a test of the consistency of the input arguments
			if (fuzzy<fp_t>::NE(fp_of_int_unchecked(max_int()), max_fp_in))
				then error_exit(ERROR_EXIT,
								"***** linear_map<fp_t>::linear_map:\n"
								"        int range [%d,%d]\n"
								"        and fp range [%g(%g)%g]\n"
								"        are (fuzzily) inconsistent!\n",
								min_int(), max_int(),
								double(min_fp_in), double(delta_fp()), double(max_fp_in));
			/*NOTREACHED*/
		}

		//******************************************************************************

		//
		// This function converts  fp  --> int  coordinate, returning the result
		// as an fp (which need not be fuzzily integral).
		//
		template <typename fp_t>
		fp_t linear_map<fp_t>::fp_int_of_fp(fp_t x)
			const
		{
			if (!is_in_range(x))
				then error_exit(ERROR_EXIT,
								"***** linear_map<fp_t>::fp_int_of_fp:\n"
								"        fp value x=%g is (fuzzily) outside the grid!\n"
								"        {min(delta)max}_fp = %g(%g)%g\n",
								double(x),
								double(min_fp()), double(delta_fp()), double(max_fp()));
			/*NOTREACHED*/

			return inverse_delta_ * (x - origin_);
		}

		//******************************************************************************

		//
		// This function converts  fp  --> int  and checks that the result is
		// fuzzily integral.  (The  nia  argument specifies what to do if the
		// result *isn't* fuzzily integral.)
		//
		// FIXME:
		// Having to explicitly specify the namespace for jtutil::round<fp_t>::
		// is ++ugly. :(
		//
		template <typename fp_t>
		int linear_map<fp_t>::int_of_fp(fp_t x, noninteger_action nia /* = nia_error */)
			const
		{
			const fp_t fp_int = fp_int_of_fp(x);

			if (fuzzy<fp_t>::is_integer(fp_int))
				then
				{
					// x is (fuzzily) a grid point ==> return that
					return jtutil::round<fp_t>::to_integer(fp_int); // *** EARLY RETURN ***
				}

			// get to here ==> x isn't (fuzzily) a grid point
			static const char *const noninteger_msg =
				"%s linear_map<fp_t>::int_of_fp:\n"
				"        x=%g isn't (fuzzily) a grid point!\n"
				"        {min(delta)max}_fp() = %g(%g)%g\n";
			switch (nia)
			{
			case nia_error:
				error_exit(ERROR_EXIT,
						   noninteger_msg,
						   "*****",
						   double(x),
						   double(min_fp()), double(delta_fp()), double(max_fp()));
				/*NOTREACHED*/

			case nia_warning:
				printf(noninteger_msg,
					   "---",
					   double(x),
					   double(min_fp()), double(delta_fp()), double(max_fp()));
				// fall through

			case nia_round:
				return jtutil::round<fp_t>::to_integer(fp_int); // *** EARLY RETURN ***

			case nia_floor:
				return jtutil::round<fp_t>::floor(fp_int); // *** EARLY RETURN ***

			case nia_ceiling:
				return jtutil::round<fp_t>::ceiling(fp_int); // *** EARLY RETURN ***

			default:
				error_exit(PANIC_EXIT,
						   "***** linear_map<fp_t>::int_of_fp: illegal nia=(int)%d\n"
						   "                                   (this should never happen!)\n",
						   int(nia)); /*NOTREACHED*/
			}
			return 0; // dummy return to quiet gcc
					  // (which doesn't grok that error_exit() never returns)
		}

		//******************************************************************************

		//
		// This function converts "delta" spacings in the fp coordinate to
		// corresponding "delta" spacings in the int coordinate, and checks that
		// the result is fuzzily integral.  (The  nia  argument specifies what to
		// do if the result *isn't* fuzzily integral.)
		//
		// FIXME:
		// Having to explicitly specify the namespace for jtutil::round<fp_t>::
		// is ++ugly. :(
		//
		template <typename fp_t>
		int linear_map<fp_t>::delta_int_of_delta_fp(fp_t delta_x, noninteger_action nia /* = nia_error */)
			const
		{
			const fp_t fp_delta_int = inverse_delta_ * delta_x;

			if (fuzzy<fp_t>::is_integer(fp_delta_int))
				then
				{
					// delta_x is (fuzzily) an integer number of grid spacings
					// ==> return that
					return jtutil::round<fp_t>::to_integer(fp_delta_int);
					// *** EARLY RETURN ***
				}

			// get to here ==> delta_x isn't (fuzzily) an integer number of grid spacings
			static const char *const noninteger_msg =
				"%s linear_map<fp_t>::delta_int_of_delta_fp:\n"
				"        delta_x=%g isn't (fuzzily) an integer number of grid spacings!\n"
				"        {min(delta)max}_fp() = %g(%g)%g\n";
			switch (nia)
			{
			case nia_error:
				error_exit(ERROR_EXIT,
						   noninteger_msg,
						   "*****",
						   double(delta_x),
						   double(min_fp()), double(delta_fp()), double(max_fp()));
				/*NOTREACHED*/

			case nia_warning:
				printf(noninteger_msg,
					   "---",
					   double(delta_x),
					   double(min_fp()), double(delta_fp()), double(max_fp()));
				// fall through

			case nia_round:
				return jtutil::round<fp_t>::to_integer(fp_delta_int);
				// *** EARLY RETURN ***

			case nia_floor:
				return jtutil::round<fp_t>::floor(fp_delta_int); // *** EARLY RETURN ***

			case nia_ceiling:
				return jtutil::round<fp_t>::ceiling(fp_delta_int);
				// *** EARLY RETURN ***

			default:
				error_exit(PANIC_EXIT,
						   "***** linear_map<fp_t>::delta_int_of_delta_fp: illegal nia=(int)%d\n"
						   "                                               (this should never happen!)\n",
						   int(nia)); /*NOTREACHED*/
			}
			return 0; // dummy return to quiet gcc
					  // (which doesn't grok that error_exit() never returns)
		}

		//******************************************************************************
		//******************************************************************************
		//******************************************************************************

		//
		// ***** template instantiation *****
		//

		template class linear_map<float>;
		template class linear_map<double>;

		//******************************************************************************
		//******************************************************************************
		//******************************************************************************

	} // namespace jtutil
} // namespace AHFinderDirect
