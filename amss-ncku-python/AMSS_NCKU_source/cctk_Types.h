#ifndef _CCTK_TYPES_H_
#define _CCTK_TYPES_H_

#ifndef _CCTK_CONFIG_H_
#include "cctk_Config.h"
#endif

typedef void *CCTK_POINTER;
typedef const void *CCTK_POINTER_TO_CONST;
typedef void (*CCTK_FPOINTER)(void);
#define HAVE_CCTK_POINTER 1
#define HAVE_CCTK_POINTER_TO_CONST 1
#define HAVE_CCTK_FPOINTER 1

/* Character types */
typedef char CCTK_CHAR;
typedef const char * CCTK_STRING;
#define HAVE_CCTK_CHAR 1
#define HAVE_CCTK_STRING 1

/* Structures for complex types */

#ifdef HAVE_CCTK_REAL16
#define HAVE_CCTK_COMPLEX32 1
typedef struct CCTK_COMPLEX32
{
  CCTK_REAL16 Re;
  CCTK_REAL16 Im;
#ifdef __cplusplus
  CCTK_REAL16 real() const { return Re; }
  CCTK_REAL16 imag() const { return Im; }
#endif
} CCTK_COMPLEX32;
#endif

#ifdef HAVE_CCTK_REAL8
#define HAVE_CCTK_COMPLEX16 1
typedef struct CCTK_COMPLEX16
{
  CCTK_REAL8 Re;
  CCTK_REAL8 Im;
#ifdef __cplusplus
  CCTK_REAL8 real() const { return Re; }
  CCTK_REAL8 imag() const { return Im; }
#endif
} CCTK_COMPLEX16;
#endif

#ifdef HAVE_CCTK_REAL4
#define HAVE_CCTK_COMPLEX8 1
typedef struct CCTK_COMPLEX8
{
  CCTK_REAL4 Re;
  CCTK_REAL4 Im;
#ifdef __cplusplus
  CCTK_REAL4 real() const { return Re; }
  CCTK_REAL4 imag() const { return Im; }
#endif
} CCTK_COMPLEX8;
#endif

/* Small positive integer type */
typedef unsigned char CCTK_BYTE;
#define HAVE_CCTK_BYTE 1

/* Define stuff for fortran. */
#ifdef FCODE

#define CCTK_POINTER          integer*SIZEOF_CHAR_P
#define CCTK_POINTER_TO_CONST integer*SIZEOF_CHAR_P
/* TODO: add autoconf for determining the size of function pointers */
#define CCTK_FPOINTER         integer*SIZEOF_CHAR_P
#define HAVE_CCTK_POINTER 1
#define HAVE_CCTK_POINTER_TO_CONST 1
#define HAVE_CCTK_FPOINTER 1

/* Character types */
/* A single character does not exist in Fortran; in Fortran, all
   character types are strings.  Hence we do not define CCTK_CHAR.  */
/* #define CCTK_CHAR   CHARACTER */
/* #define HAVE_CCTK_CHAR 1 */
/* This is a C-string, i.e., only a pointer */
#define CCTK_STRING CCTK_POINTER_TO_CONST
#define HAVE_CCTK_STRING 1

#ifdef HAVE_CCTK_INT8
#define CCTK_INT8 INTEGER*8
#endif
#ifdef HAVE_CCTK_INT4
#define CCTK_INT4 INTEGER*4
#endif
#ifdef HAVE_CCTK_INT2
#define CCTK_INT2 INTEGER*2
#endif
#ifdef HAVE_CCTK_INT1
#define CCTK_INT1 INTEGER*1
#endif

#ifdef HAVE_CCTK_REAL16
#define CCTK_REAL16 REAL*16
#define HAVE_CCTK_COMPLEX32 1
#define CCTK_COMPLEX32  COMPLEX*32
#endif

#ifdef HAVE_CCTK_REAL8
#define CCTK_REAL8  REAL*8
#define HAVE_CCTK_COMPLEX16 1
#define CCTK_COMPLEX16  COMPLEX*16
#endif

#ifdef HAVE_CCTK_REAL4
#define CCTK_REAL4  REAL*4
#define HAVE_CCTK_COMPLEX8 1
#define CCTK_COMPLEX8   COMPLEX*8
#endif

/* Should be unsigned, but Fortran doesn't have that */
#define CCTK_BYTE INTEGER*1
#define HAVE_CCTK_BYTE 1

#endif /*FCODE */

/* Now pick the types based upon the precision variable. */

/* Floating point precision */
#ifdef CCTK_REAL_PRECISION_16
#define CCTK_REAL_PRECISION 16
#define CCTK_REAL CCTK_REAL16
#endif

#ifdef CCTK_REAL_PRECISION_8
#define CCTK_REAL_PRECISION 8
#define CCTK_REAL CCTK_REAL8
#endif

#ifdef CCTK_REAL_PRECISION_4
#define CCTK_REAL_PRECISION 4
#define CCTK_REAL CCTK_REAL4
#endif

/* Integer precision */

#ifdef CCTK_INTEGER_PRECISION_8
#define CCTK_INTEGER_PRECISION 8
#define CCTK_INT CCTK_INT8
#endif

#ifdef CCTK_INTEGER_PRECISION_4
#define CCTK_INTEGER_PRECISION 4
#define CCTK_INT CCTK_INT4
#endif

#ifdef CCTK_INTEGER_PRECISION_2
#define CCTK_INTEGER_PRECISION 2
#define CCTK_INT CCTK_INT2
#endif

#ifdef CCTK_INTEGER_PRECISION_1
#define CCTK_INTEGER_PRECISION 1
#define CCTK_INT CCTK_INT1
#endif

/* Complex precision */
#ifdef CCTK_REAL_PRECISION_16
#define CCTK_COMPLEX_PRECISION 32
#define CCTK_COMPLEX CCTK_COMPLEX32
#endif

#ifdef CCTK_REAL_PRECISION_8
#define CCTK_COMPLEX_PRECISION 16
#define CCTK_COMPLEX CCTK_COMPLEX16
#endif

#ifdef CCTK_REAL_PRECISION_4
#define CCTK_COMPLEX_PRECISION 8
#define CCTK_COMPLEX CCTK_COMPLEX8
#endif

#endif /*_CCTK_TYPES_H_ */

