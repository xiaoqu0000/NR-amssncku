/* $Id: write_Pablo_file.C,v 1.1.1.1 2017/09/25 02:35:55 zjcao Exp $ */
#ifdef newc
#include <sstream>
#include <cstdio>
#include <map>
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
#include <stdio.h>
#include <map.h>
#endif

#include "comm.h"

#include "../For_Debug/tillherecheck.h"

//|----------------------------------------------------------------------------
//  write ASCII file with the style of Pablo
//|----------------------------------------------------------------------------
void TwoPunctures::write_Pablo_file(const int *ext,
		                    const double xmin,const double xmax,
				    const double ymin,const double ymax,
				    const double zmin,const double zmax,
	       			    const char *filename)
{
    int nx = ext[0],ny = ext[1],nz = ext[2];
    int i,j,k;
    double *X, *Y, *Z;
    X = new double[nx];
    Y = new double[ny];
    Z = new double[nz];
    double dX,dY,dZ;
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
    dX = (xmax - xmin)/(nx-1);
    for (i = 0; i < nx; i++) X[i] = xmin + i*dX;
    dY = (ymax - ymin)/(ny-1);
    for (j = 0; j < ny; j++) Y[j] = ymin + j*dY;
    dZ = (zmax - zmin)/(nz-1);
    for (k = 0; k < nz; k++) Z[k] = zmin + k*dZ;
#else
#ifdef Cell
    dX = (xmax - xmin)/nx;
    for (i = 0; i < nx; i++) X[i] = xmin + (i+0.5)*dX;
    dY = (ymax - ymin)/ny;
    for (j = 0; j < ny; j++) Y[j] = ymin + (j+0.5)*dY;
    dZ = (zmax - zmin)/nz;
    for (k = 0; k < nz; k++) Z[k] = zmin + (k+0.5)*dZ;
#else
#error Not define Vertex nor Cell
#endif  
#endif
//|--->open out put file  
    ofstream outfile; 
    outfile.open(filename);
    if(!outfile)
    {
      cout << "write_Pablo_file can't open " << filename << " for output." << endl;
      exit(0);
    }

    outfile<<"#using center of mass coordinate and the BBHs locate on x-axis"<<endl;
    outfile.setf(ios::scientific,ios::floatfield);
    outfile.precision( 16 );

    double D=2*par_b,x1,x2;
    x1= D*target_M_minus/(target_M_plus+target_M_minus);
    x2=-D*target_M_plus/(target_M_plus+target_M_minus);

    for(k=0;k<nz;k++)
      for(j=0;j<ny;j++)
	for(i=0;i<nx;i++)
	{
// using center of mass coordinate
	   outfile<<X[i]     <<" "<<Y[j]     <<" "<<Z[k]     <<" "
		  <<PunctIntPolAtArbitPosition(0,1,npoints_A, npoints_B, npoints_phi,v,X[i]-(x1+x2)/2,Y[j],Z[k])<<endl;
	}

    outfile.close();

    delete[] X; delete[] Y; delete[] Z;
}
void TwoPunctures::write_Pablo_file_F(const int *ext,
		                    const double xmin,const double xmax,
				    const double ymin,const double ymax,
				    const double zmin,const double zmax,
	       			    const char *filename)
{
    int nx = ext[0],ny = ext[1],nz = ext[2];
    int i,j,k;
    double *X, *Y, *Z;
    X = new double[nx];
    Y = new double[ny];
    Z = new double[nz];
    double dX,dY,dZ;
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
    dX = (xmax - xmin)/(nx-1);
    for (i = 0; i < nx; i++) X[i] = xmin + i*dX;
    dY = (ymax - ymin)/(ny-1);
    for (j = 0; j < ny; j++) Y[j] = ymin + j*dY;
    dZ = (zmax - zmin)/(nz-1);
    for (k = 0; k < nz; k++) Z[k] = zmin + k*dZ;
#else
#ifdef Cell
    dX = (xmax - xmin)/nx;
    for (i = 0; i < nx; i++) X[i] = xmin + (i+0.5)*dX;
    dY = (ymax - ymin)/ny;
    for (j = 0; j < ny; j++) Y[j] = ymin + (j+0.5)*dY;
    dZ = (zmax - zmin)/nz;
    for (k = 0; k < nz; k++) Z[k] = zmin + (k+0.5)*dZ;
#else
#error Not define Vertex nor Cell
#endif  
#endif
//|--->open out put file  
    ofstream outfile; 
    outfile.open(filename);
    if(!outfile)
    {
      cout << "write_Pablo_file can't open " << filename << " for output." << endl;
      exit(0);
    }

    outfile<<"#using center of mass coordinate and the BBHs locate on x-axis"<<endl;
    outfile.setf(ios::scientific,ios::floatfield);
    outfile.precision( 16 );

    double D=2*par_b,x1,x2;
    x1= D*target_M_minus/(target_M_plus+target_M_minus);
    x2=-D*target_M_plus/(target_M_plus+target_M_minus);

    for(k=0;k<nz;k++)
      for(j=0;j<ny;j++)
	for(i=0;i<nx;i++)
	{
// using center of mass coordinate
	   outfile<<X[i]     <<" "<<Y[j]     <<" "<<Z[k]     <<" "
		  <<PunctIntPolAtArbitPosition_F(0,1,npoints_A, npoints_B, npoints_phi,F,X[i]-(x1+x2)/2,Y[j],Z[k])<<endl;
	}

    outfile.close();

    delete[] X; delete[] Y; delete[] Z;
}
/* -------------------------------------------------------------------------*/
/* Calculates the value of residual at an arbitrary position (x,y,z)*/
// note here (x,y,z) means BBHs locate at (\pm b,0,0)
double TwoPunctures::PunctIntPolAtArbitPosition_F(int ivar, int nvar, int n1,
			    int n2, int n3, double *Ff, double x, double y,
			    double z)
{
  double xs, ys, zs, rs2, phi, X, R, A, B, aux1, aux2, result, Ui;

  xs = x / par_b;
  ys = y / par_b;
  zs = z / par_b;
  rs2 = ys * ys + zs * zs;
  phi = atan2 (z, y);
  if (phi < 0)  phi += 2 * Pi;

  aux1 = 0.5 * (xs * xs + rs2 - 1);
  aux2 = sqrt (aux1 * aux1 + rs2);
  X = asinh (sqrt (aux1 + aux2));
  R = asin (min(1.0, sqrt (-aux1 + aux2)));
  if (x < 0) R = Pi - R;

  A = 2 * tanh (0.5 * X) - 1;
  B = tan (0.5 * R - Piq);

  result = PunctEvalAtArbitPosition (F, ivar, A, B, phi, nvar, n1, n2, n3);

  Ui = (A - 1) * result;

  return Ui;
}
