#include <stdlib.h>

#include "stdc.h"
#include "util.h"

namespace AHFinderDirect
{
       namespace jtutil
       {
              template <typename fp_t>
              int round<fp_t>::to_integer(fp_t x)
              {
                     return (x >= 0.0)
                                ? int(x + 0.5)      // eg 3.6 --> int(4.1) = 4
                                : -int((-x) + 0.5); // eg -3.6 --> - int(4.1) = -4
              }

              template <typename fp_t>
              int round<fp_t>::floor(fp_t x)
              {
                     return (x >= 0.0)
                                ? int(x)
                                : -ceiling(-x);
              }

              template <typename fp_t>
              int round<fp_t>::ceiling(fp_t x)
              {
                     return (x >= 0.0)
                                ? int(x) + (x != fp_t(int(x)))
                                : -floor(-x);
              }

              template class round<float>;
              template class round<double>;

       } // namespace jtutil
} // namespace AHFinderDirect
