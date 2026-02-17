#include <math.h>
#include <stdlib.h>

#include "cctk.h"

#include "stdc.h"
#include "util.h"

namespace AHFinderDirect
{
	namespace jtutil
	{
		double signum(double x)
		{
			if (x == 0.0)
				then return 0.0;
			else
				return (x > 0.0) ? 1.0 : -1.0;
		}
		double hypot3(double x, double y, double z)
		{
			return sqrt(x * x + y * y + z * z);
		}
		double arctan_xy(double x, double y)
		{
			return ((x == 0.0) && (y == 0.0)) ? 0.0 : atan2(y, x);
		}
		double modulo_reduce(double x, double xmod, double xmin, double xmax)
		{
			double xx = x;

			while (fuzzy<double>::LT(xx, xmin))
			{
				xx += xmod;
			}

			while (fuzzy<double>::GT(xx, xmax))
			{
				xx -= xmod;
			}

			if (!(fuzzy<double>::GE(xx, xmin) && fuzzy<double>::LE(xx, xmax)))
				then error_exit(ERROR_EXIT,
								"***** modulo_reduce(): no modulo value is fuzzily within specified range!\n"
								"                       x = %g   xmod = %g\n"
								"                       [xmin,xmax] = [%g,%g]\n"
								"                       ==> xx = %g\n",
								x, xmod,
								xmin, xmax,
								xx); /*NOTREACHED*/

			return xx;
		}
		template <typename fp_t>
		void zero_C_array(int N, fp_t array[])
		{
			for (int i = 0; i < N; ++i)
			{
				array[i] = 0;
			}
		}

		template void zero_C_array<CCTK_REAL>(int, CCTK_REAL[]);

	} // namespace jtutil
} // namespace AHFinderDirect
