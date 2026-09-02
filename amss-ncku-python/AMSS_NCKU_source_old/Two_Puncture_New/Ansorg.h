//$Id: Ansorg.h,v 1.1.1.1 2017/09/25 02:35:55 zjcao Exp $
#ifndef Ansorg_H
#define Ansorg_H

#ifdef newc
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
using namespace std;
#else
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#endif

#include <mpi.h>

#define PI M_PI

   class Ansorg
{
  protected:
      int n1,n2,n3,ntotal;
      int order;
      double *coordA, *coordB, *coordphi;
      int ps_rxx,ps_rxy,ps_ryx,ps_ryy;
      double ps_b,ps_dx;
      double PIh;
      double *pu_ps;
      int myrank;
  public:
       Ansorg(char*filename,int orderi);
       ~Ansorg();
       double ps_u_at_xyz(double x, double y, double z);
       void write_Pablo_file(const int *ext,const double xmin,const double xmax,const double ymin,const double ymax,const double zmin,const double zmax,
		     const char *filename);
       void set_ABp();
       void xyz_to_ABp(double x, double y, double z,
		       double *A, double *B, double *phi);
       double interpolate_tri_bar(double x, double y, double z,
			   int n1, int n2, int n3, 
			   double *x1, double *x2, double *x3, double *yp);
       int find_point_bisection(double x, int n, double *xp, int o);
       void barycentric_omega(int n, int s, double *x, double *omega);
       double barycentric(double x0, int n, int s, double *x, double *y,
		   double *omega);

  private:
       bool zheadon;
};
#endif    /* Ansorg_H */
