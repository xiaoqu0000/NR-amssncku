
#ifndef MICRODEF_H
#define MICRODEF_H

#include "microdef.fh"

// application parameters

/// ****
// sommerfeld boundary type
// 0: bam, 1: shibata
#define SommerType 0

/// ****
// for Using Gauss-Legendre quadrature in theta direction
#define GaussInt

/// ****
// 0: BSSN vacuum
// 1: coupled to scalar field
// 2: Z4c vacuum
// 3: coupled to Maxwell field
//
#define ABEtype 2

/// ****
// using Apparent Horizon Finder
//#define With_AHF

/// ****
// Psi4 calculation method
// 0: EB method
// 1: 4-D method
//
#define Psi4type 0

/// ****
// for Using point psi4 or not
//#define Point_Psi4

/// ****
// RestrictProlong in Step (0) or after Step (1)
#define RPS 1

/// ****
// Enforce algebra constraint
// for every RK4 sub step: 0
// only when iter_count == 3: 1
// after routine Step: 2
#define AGM 0

/// ****
// Restrict Prolong using BAM style 1 or old style 0
#define RPB 0

/// ****
// 1: move Analysis out ot 4 sub steps and treat PBH with Euler method
#define MAPBH 1

/// ****
// parallel structure, 0: level by level, 1: considering all levels, 2: as 1 but reverse the CPU order, 3: Frank's scheme
#define PSTR 0

/// ****
// regrid for every level or for all levels at a time
// 0: for every level; 1: for all
#define REGLEV 0

/// ****
// use gpu or not
//#define USE_GPU

/// ****
// use checkpoint for every process
//#define CHECKDETAIL

/// ****
// use FakeCheckPrepare to write CheckPoint
//#define FAKECHECK
////================================================================
//  some basic parameters for numerical calculation
#define dim 3

//#define Cell or Vertex in "microdef.fh"

// ******
// buffer point number for mesh refinement interface
#define buffer_width 6

// ******
// buffer point number shell-box interface, on shell
#define SC_width buffer_width
// buffer point number shell-box interface, on box
#define CS_width (2*buffer_width)

#if(buffer_width < ghost_width)
#error we always assume buffer_width>ghost_width
#endif

#define PACK 1
#define UNPACK 2

#define Mymax(a,b) (((a) > (b)) ? (a) : (b))
#define Mymin(a,b) (((a) < (b)) ? (a) : (b))

#define feq(a,b,d) (fabs(a-b)<d)
#define flt(a,b,d) ((a-b)<d)
#define fgt(a,b,d) ((a-b)>d)

#define TINY 1e-10

#endif   /* MICRODEF_H */
