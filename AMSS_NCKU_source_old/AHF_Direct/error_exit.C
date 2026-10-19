#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#include "cctk.h"

#include "config.h"
#include "stdc.h"

namespace AHFinderDirect
{
	namespace jtutil
	{
		int error_exit(int msg_level, const char *format, ...)
		{
			const int N_buffer = 2000;
			char buffer[N_buffer];

			va_list ap;
			va_start(ap, format);
			vsnprintf(buffer, N_buffer, format, ap);
			va_end(ap);

			const int len = strlen(buffer);
			if ((len > 0) && (buffer[len - 1] == '\n'))
				then buffer[len - 1] = '\0';

			CCTK_VWarn(msg_level, __LINE__, __FILE__, CCTK_THORNSTRING, "%s", buffer);

			// if we got here, evidently  msg_level  wasn't drastic enough
			abort(); /*NOTREACHED*/
		}

		//******************************************************************************

	} // namespace jtutil
} // namespace AHFinderDirect
