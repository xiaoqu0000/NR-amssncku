/* $Id: comm.h,v 1.1.1.1 2017/09/25 02:35:55 zjcao Exp $ */

#define Cell

#include "TwoPunctures.h"

//|----------------------------------------------------------------------------
//  write ASCII file with the style of Pablo
//|----------------------------------------------------------------------------
void write_Pablo_file(const int *ext,const double xmin,const double xmax,const double ymin,const double ymax,const double zmin,const double zmax,
		     const char *filename, TwoPunctures *ADM,
		     const int n1,const int n2,const int n3);
