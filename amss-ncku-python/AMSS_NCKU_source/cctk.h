#ifndef _CCTK_H_
#define _CCTK_H_ 1

/* Grab the main configuration info. */
#include "cctk_Config.h"

#define CCTK_THORNSTRING "AHFinderDirect"

/* Include the constants */
#include "cctk_Constants.h"

/* get the definition of ptrdiff_t */
#include <stddef.h>
int CCTK_VInfo(const char *thorn, const char *format, ...);
int CCTK_VWarn(int level,
               int line,
               const char *file,
               const char *thorn,
               const char *format,
               ...);
#define CCTK_ERROR_INTERP_GHOST_SIZE_TOO_SMALL (-1001)
#ifdef __cplusplus
#define HAVE_INLINE
#else
#ifndef inline
#define HAVE_INLINE
#endif
#endif

#define CCTK_PRINTSEPARATOR \
  printf("--------------------------------------------------------------------------------\n");

#define _DECLARE_CCTK_ARGUMENTS _DECLARE_CCTK_CARGUMENTS
#define _DECLARE_CCTK_CARGUMENTS          \
  ptrdiff_t cctki_dummy_int;              \
  CCTK_REAL cctk_time = cctkGH->PhysTime; \
  int cctk_iteration = 1;                 \
  int cctk_dim = 3;

#define CCTK_EQUALS(a, b) (CCTK_Equals((a), (b)))

#define CCTK_PASS_CTOC cctkGH

#define CCTK_ORIGIN_SPACE(x) (cctk_origin_space[x] + cctk_delta_space[x] / cctk_levfac[x] * cctk_levoff[x] / cctk_levoffdenom[x])
#define CCTK_DELTA_SPACE(x) (cctk_delta_space[x] / cctk_levfac[x])
#define CCTK_DELTA_TIME (cctk_delta_time / cctk_timefac)
#define CCTK_LSSH(stag, dim) cctk_lssh[(stag) + CCTK_NSTAGGER * (dim)]
#define CCTK_LSSH_IDX(stag, dim) ((stag) + CCTK_NSTAGGER * (dim))

#define CCTK_WARN(a, b) CCTK_Warn(a, __LINE__, __FILE__, CCTK_THORNSTRING, b)

#define CCTK_MALLOC(s) CCTKi_Malloc(s, __LINE__, __FILE__)
#define CCTK_FREE(p) CCTKi_Free(p)

#define CCTK_INFO(a) CCTK_Info(CCTK_THORNSTRING, (a))
#define CCTK_PARAMWARN(a) CCTK_ParamWarn(CCTK_THORNSTRING, (a))

#endif
