#ifndef AHFINDERDIRECT__STDC_H
#define AHFINDERDIRECT__STDC_H

#define then /* empty */

#ifdef M_PI
#define PI M_PI
#endif

#define iabs(x_) abs(x_)

namespace AHFinderDirect
{
	namespace jtutil
	{

		int error_exit(int msg_level, const char *format, ...);

#define ERROR_EXIT (-1)
#define PANIC_EXIT (-2)
	}
}

#endif /* AHFINDERDIRECT__STDC_H */
