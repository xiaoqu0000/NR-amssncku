#include <math.h>
#include <assert.h>
#include <stdlib.h>

#include "util.h"

namespace AHFinderDirect
{
	namespace jtutil
	{

		template <typename fp_t>
		norm<fp_t>::norm()
			: N_(0L),
			  sum_(0.0), sum2_(0.0),
			  max_abs_value_(0.0), min_abs_value_(0.0),
			  max_value_(0.0), min_value_(0.0)
		{
		}

		template <typename fp_t>
		void norm<fp_t>::reset()
		{
			N_ = 0L;
			sum_ = 0.0;
			sum2_ = 0.0;
			max_abs_value_ = 0.0;
			min_abs_value_ = 0.0;
			max_value_ = 0.0;
			min_value_ = 0.0;
		}

		template <typename fp_t>
		void norm<fp_t>::data(fp_t x)
		{
			sum_ += x;
			sum2_ += x * x;

			const fp_t abs_x = jtutil::abs<fp_t>(x);
			max_abs_value_ = jtutil::tmax(max_abs_value_, abs_x);
			min_abs_value_ = (N_ == 0) ? abs_x : jtutil::tmin(min_abs_value_, abs_x);

			min_value_ = (N_ == 0) ? x : jtutil::tmin(min_value_, x);
			max_value_ = (N_ == 0) ? x : jtutil::tmax(max_value_, x);

			++N_;
		}

		template <typename fp_t>
		fp_t norm<fp_t>::mean() const { return sum_ / fp_t(N_); }
		template <typename fp_t>
		fp_t norm<fp_t>::two_norm() const { return sqrt(sum2_); }
		template <typename fp_t>
		fp_t norm<fp_t>::rms_norm() const
		{
			assert(is_nonempty());
			return sqrt(sum2_ / fp_t(N_));
		}

		template class jtutil::norm<float>;
		template class jtutil::norm<double>;

		//******************************************************************************
		//******************************************************************************
		//******************************************************************************

	} // namespace jtutil
} // namespace AHFinderDirect
