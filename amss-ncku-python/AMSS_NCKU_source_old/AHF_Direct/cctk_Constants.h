#ifndef _CCTK_CONSTANTS_H_
#define _CCTK_CONSTANTS_H_

#define CCTK_VARIABLE_VOID             100
#define CCTK_VARIABLE_BYTE             101
#define CCTK_VARIABLE_INT              102
#define CCTK_VARIABLE_INT1             103
#define CCTK_VARIABLE_INT2             104
#define CCTK_VARIABLE_INT4             105
#define CCTK_VARIABLE_INT8             106
#define CCTK_VARIABLE_REAL             107
#define CCTK_VARIABLE_REAL4            108
#define CCTK_VARIABLE_REAL8            109
#define CCTK_VARIABLE_REAL16           110
#define CCTK_VARIABLE_COMPLEX          111
#define CCTK_VARIABLE_COMPLEX8         112
#define CCTK_VARIABLE_COMPLEX16        113
#define CCTK_VARIABLE_COMPLEX32        114
#define CCTK_VARIABLE_CHAR             115
#define CCTK_VARIABLE_STRING           116
#define CCTK_VARIABLE_POINTER          117
#define CCTK_VARIABLE_POINTER_TO_CONST 118
#define CCTK_VARIABLE_FPOINTER         119

/* DEPRECATED IN BETA 12 */
#define CCTK_VARIABLE_FN_POINTER CCTK_VARIABLE_FPOINTER

/* steerable status of parameters */
#define CCTK_STEERABLE_NEVER   200
#define CCTK_STEERABLE_ALWAYS  201
#define CCTK_STEERABLE_RECOVER 202

/* number of staggerings */
#define CCTK_NSTAGGER      3

/* group distributions */
#define CCTK_DISTRIB_CONSTANT 301
#define CCTK_DISTRIB_DEFAULT  302

/* group types */
#define CCTK_SCALAR 401
#define CCTK_GF     402
#define CCTK_ARRAY  403

/* group scopes */
#define CCTK_PRIVATE   501
#define CCTK_PROTECTED 502
#define CCTK_PUBLIC    503

/* constants for CCTK_TraverseString() */
#define CCTK_VAR          601
#define CCTK_GROUP        602
#define CCTK_GROUP_OR_VAR 603


#endif /* _CCTK_CONSTANTS_ */

