//$Id: lagpolint.h,v 1.1.1.1 2017/09/25 02:35:56 zjcao Exp $
#ifndef LAGPOLINT_H
#define LAGPOLINT_H

#ifdef fortran1
#define f_polint polint
#endif
#ifdef fortran2
#define f_polint POLINT
#endif
#ifdef fortran3
#define f_polint polint_
#endif

extern "C" { void f_polint(double *,double *,double &,double &,double &,int &); }

#endif    /* LAGPOLINT_H */
