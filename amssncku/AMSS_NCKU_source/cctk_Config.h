#ifndef _CCTK_CONFIG_H_
#define _CCTK_CONFIG_H_

#define STDC_HEADERS 1

#define CCTK_FCALL 

#define HAVE_GETHOSTBYNAME 1
#define HAVE_GETOPT_LONG_ONLY 1
#define HAVE_CRYPT 1
#define HAVE_FINITE 1
#define HAVE_ISNAN 1
#define HAVE_ISINF 1
#define HAVE_MKSTEMP 1
#define HAVE_VA_COPY 1

/* Do we have mode_t ? */
#define HAVE_MODE_T 1

#define HAVE_SOCKLEN_T 1
#ifdef HAVE_SOCKLEN_T
#  define CCTK_SOCKLEN_T socklen_t
#else
#  define CCTK_SOCKLEN_T int
#endif

#define HAVE_TIME_H 1
#define HAVE_SYS_IOCTL_H 1
#define HAVE_SYS_SOCKET_H 1
#define HAVE_SYS_TIME_H 1
#define HAVE_SYS_TYPES_H 1
#define HAVE_UNISTD_H 1
#define HAVE_STRING_H 1
#define HAVE_ASSERT_H 1
#define HAVE_TGMATH_H 1
#define HAVE_SYS_STAT_H 1
#define HAVE_GETOPT_H 1
#define HAVE_REGEX_H 1
#define HAVE_NETINET_IN_H 1
#define HAVE_NETDB_H 1
#define HAVE_ARPA_INET_H 1
#define HAVE_CRYPT_H 1
#define HAVE_DIRENT_H 1
#define HAVE_SIGNAL_H 1
#define HAVE_MALLOC_H 1
#define HAVE_MALLINFO 1
#define HAVE_MALLOPT 1
#define HAVE_M_MMAP_THRESHOLD_VALUE 1

#define TIME_WITH_SYS_TIME 1

#define HAVE_VECTOR 1
#define HAVE_VECTOR_H 1

#define GETTIMEOFDAY_NEEDS_TIMEZONE 1

#define CCTK_CACHELINE_BYTES 64
#define CCTK_CACHE_SIZE 1024*1024

#define NULL_DEVICE "/dev/null"

#define CCTK_BUILD_OS "linux-gnu"
#define CCTK_BUILD_CPU "x86_64"
#define CCTK_BUILD_VENDOR "unknown"

#define SIZEOF_SHORT_INT 2
#define SIZEOF_INT 4
#define SIZEOF_LONG_INT 8
#define SIZEOF_LONG_LONG 8
#define SIZEOF_LONG_DOUBLE 16
#define SIZEOF_DOUBLE 8
#define SIZEOF_FLOAT 4
#define SIZEOF_CHAR_P 8

#define CCTK_REAL_PRECISION_8 1

#define CCTK_INTEGER_PRECISION_4 1

#define HAVE_CCTK_INT8 1
#define HAVE_CCTK_INT4 1
#define HAVE_CCTK_INT2 1
#define HAVE_CCTK_INT1 1

#define HAVE_CCTK_REAL16 1
#define HAVE_CCTK_REAL8 1
#define HAVE_CCTK_REAL4 1

#define CCTK_INT8 long int
#define CCTK_INT4 int
#define CCTK_INT2 short int
#define CCTK_INT1 signed char

#define CCTK_REAL16 long double
#define CCTK_REAL8 double
#define CCTK_REAL4 float

#ifndef __cplusplus

#ifdef CCTK_C_RESTRICT
#define restrict CCTK_C_RESTRICT
#endif

/* Allow the use of CCTK_RESTRICT as a qualifier always. */
#ifdef CCTK_C_RESTRICT
#define CCTK_RESTRICT CCTK_C_RESTRICT
#else
#define CCTK_RESTRICT restrict
#endif

#ifdef HAVE_CCTK_C_BOOL
#define CCTK_HAVE_C_BOOL
#endif

#endif /* ! defined __cplusplus */
/****************************************************************************/

/****************************************************************************/
/* C++ specific stuff */
/****************************************************************************/
#ifdef __cplusplus

/* Some C++ compilers don't have bool ! */
#define HAVE_CCTK_CXX_BOOL 1

#ifndef HAVE_CCTK_CXX_BOOL
typedef enum {false, true} bool;
#else
/* deprecated in beta15 */
#define CCTK_HAVE_CXX_BOOL
#endif

/* Some C++ compilers recognise the restrict keyword */
#define CCTK_CXX_RESTRICT __restrict__

/* Since this is non-standard leave commented out for the moment */
#if 0
/* Define to empty if the keyword does not work. */
#ifdef CCTK_CXX_RESTRICT
#define restrict CCTK_CXX_RESTRICT
#endif
#endif

/* Allow the use of CCTK_RESTRICT as a qualifier always. */
#ifdef CCTK_CXX_RESTRICT
#define CCTK_RESTRICT CCTK_CXX_RESTRICT
#else
#define CCTK_RESTRICT restrict
#endif

#endif /* __cplusplus */
/****************************************************************************/

#ifdef FCODE

#define HAVE_CCTK_FORTRAN_REAL4 1
#define HAVE_CCTK_FORTRAN_REAL8 1
#define HAVE_CCTK_FORTRAN_REAL16 1

#define HAVE_CCTK_FORTRAN_COMPLEX8 1
#define HAVE_CCTK_FORTRAN_COMPLEX16 1
#define HAVE_CCTK_FORTRAN_COMPLEX32 1

#endif /* FCODE */

/* Now include the code to pick an appropriate precison for reals and ints */
#include "cctk_Types.h"

#endif /* _CCTK_CONFIG_H_ */
