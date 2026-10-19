#include <stdlib.h>
#include <stdio.h>

#include "stdc.h"
#include "util.h"

namespace AHFinderDirect
{
  namespace jtutil
  {
    template <typename fp_t>
    bool fuzzy<fp_t>::EQ(fp_t x, fp_t y)
    {
      fp_t max_abs = jtutil::tmax(jtutil::abs(x), jtutil::abs(y));
      fp_t epsilon = jtutil::tmax(tolerance_, tolerance_ * max_abs);

      return jtutil::abs(x - y) <= epsilon;
    }

    //******************************************************************************

    template <typename fp_t>
    bool fuzzy<fp_t>::is_integer(fp_t x)
    {
      int i = round<fp_t>::to_integer(x);
      return EQ(x, fp_t(i));
    }

    //******************************************************************************

    template <typename fp_t>
    int fuzzy<fp_t>::floor(fp_t x)
    {
      return fuzzy<fp_t>::is_integer(x)
                 ? round<fp_t>::to_integer(x)
                 : round<fp_t>::floor(x);
    }

    //******************************************************************************

    template <typename fp_t>
    int fuzzy<fp_t>::ceiling(fp_t x)
    {
      return fuzzy<fp_t>::is_integer(x)
                 ? round<fp_t>::to_integer(x)
                 : round<fp_t>::ceiling(x);
    }
    template <>
    float fuzzy<float>::tolerance_ = 1.0e-5; // about 100 * FLT_EPSILON

    template <>
    double fuzzy<double>::tolerance_ = 1.0e-12; // about 1e4 * DBL_EPSILON

    // template instantiations
    template class fuzzy<float>;
    template class fuzzy<double>;

    //******************************************************************************
    //******************************************************************************
    //******************************************************************************

  } // namespace jtutil
} // namespace AHFinderDirect
