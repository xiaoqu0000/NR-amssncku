
//-----------------------------------------------------------------------
// Read binary files and do fancy things with them...
//-----------------------------------------------------------------------
#ifdef newc
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>
using namespace std;
#else
#include <math.h>
#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fstream.h>
#endif

#include "microdef.fh"

int main(int argc, char *argv[])
{
  //
  // USE:  DataCT flag file1 [ file2 ]
  //
  // where: - flag can be XY,XZ,YZ
  //
  void set_fname(char *fname);

  if (argc < 3)
  {
    cout << "\aUsage: DataCT flag binaryfile1 [ binaryfile2 ] \n "
         << "   where: - flag can be XY,XZ,YZ"
         << endl;
    exit(1);
  }
  ifstream infile1;
  infile1.open(argv[2]);
  if (!infile1)
  {
    cerr << "\a Can't open " << argv[2] << " for input." << endl;
    exit(1);
  }

  /* read properties of the binary file */
  double time;
  int nx, ny, nz;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  infile1.seekg(0, ios::beg);
  infile1.read((char *)&time, sizeof(double));
  infile1.read((char *)&nx, sizeof(int));
  infile1.read((char *)&ny, sizeof(int));
  infile1.read((char *)&nz, sizeof(int));
  infile1.read((char *)&xmin, sizeof(double));
  infile1.read((char *)&xmax, sizeof(double));
  infile1.read((char *)&ymin, sizeof(double));
  infile1.read((char *)&ymax, sizeof(double));
  infile1.read((char *)&zmin, sizeof(double));
  infile1.read((char *)&zmax, sizeof(double));

  /* get rid of any 4 character suffix */
  set_fname(argv[2]);

  /* sanity check */
  if (nx != ny || nx != nz)
  {
    cout << "\n"
         << endl;
    cout << " nx, ny and nz do not agree! Using a symmetry?... ";
    cout << "\n"
         << endl;
  }

  cout << "\n Reading file : " << argv[2] << endl;
  cout << "\n  Time        : " << time << endl;
  cout << "  Dimensions  : " << setw(16) << nx << setw(16) << ny << setw(16) << nz << endl;
  cout << "   xmin, xmax : " << setw(16) << xmin << setw(16) << xmax << endl;
  cout << "   ymin, ymax : " << setw(16) << ymin << setw(16) << ymax << endl;
  cout << "   zmin, zmax : " << setw(16) << zmin << setw(16) << zmax << endl;
  cout << "\n";

  double *data;
  data = new double[nx * ny * nz];
  int i = 0, j = 0, k = 0;
  infile1.read((char *)data, nx * ny * nz * sizeof(double));
  infile1.close();
  //
  //
  // if second file given, open second file and subtract from first one!
  //
  //
  if (argc == 4)
  {
    infile1.open(argv[3]);
    if (!infile1)
    {
      cerr << "\a Can't open " << argv[3] << " for input." << endl;
      exit(1);
    }
    double *indata;
    indata = new double[nx * ny * nz];
    // read in header
    infile1.seekg(0, ios::beg);
    int nxin, nyin, nzin;
    infile1.read((char *)&time, sizeof(double));
    infile1.read((char *)&nxin, sizeof(int));
    infile1.read((char *)&nyin, sizeof(int));
    infile1.read((char *)&nzin, sizeof(int));
    infile1.read((char *)&xmin, sizeof(double));
    infile1.read((char *)&xmax, sizeof(double));
    infile1.read((char *)&ymin, sizeof(double));
    infile1.read((char *)&ymax, sizeof(double));
    infile1.read((char *)&zmin, sizeof(double));
    infile1.read((char *)&zmax, sizeof(double));
    if (nxin != nx || nyin != ny || nzin != nz)
    {
      cerr << "\a Number of indices do not agree! " << endl;
      exit(1);
    }
    cout << " Comparing with data at time " << time << "\n"
         << endl;
    infile1.read((char *)indata, nx * ny * nz * sizeof(double));
    infile1.close();
    for (i = 0; i < nx * ny * nz; i++)
      data[i] -= indata[i];
  }

  double *X, *Y, *Z;
  X = new double[nx];
  Y = new double[ny];
  Z = new double[nz];
  double dd;
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
  dd = (xmax - xmin) / (nx - 1);
  for (i = 0; i < nx; i++)
    X[i] = xmin + i * dd;
  dd = (ymax - ymin) / (ny - 1);
  for (j = 0; j < ny; j++)
    Y[j] = ymin + j * dd;
  dd = (zmax - zmin) / (nz - 1);
  for (k = 0; k < nz; k++)
    Z[k] = zmin + k * dd;
#else
#ifdef Cell
  dd = (xmax - xmin) / nx;
  for (i = 0; i < nx; i++)
    X[i] = xmin + (i + 0.5) * dd;
  dd = (ymax - ymin) / ny;
  for (j = 0; j < ny; j++)
    Y[j] = ymin + (j + 0.5) * dd;
  dd = (zmax - zmin) / nz;
  for (k = 0; k < nz; k++)
    Z[k] = zmin + (k + 0.5) * dd;
#else
#error Not define Vertex nor Cell
#endif
#endif

  int ext[3];
  ext[0] = nx;
  ext[1] = ny;
  ext[2] = nz;
  void writefile(int *ext, double *XX, double *YY, double *ZZ, double *datain,
                 char *filename, const char *flag);
  writefile(ext, X, Y, Z, data, argv[2], argv[1]);

  delete[] data;
  delete[] X;
  delete[] Y;
  delete[] Z;
}

/*-----------------------------------*/
/* get rid of any 4 character suffix */
/*-----------------------------------*/
void set_fname(char *fname)
{
  int len = strlen(fname) - 4;
  char *n_fname;
  n_fname = new char[len];

  for (int i = 0; i < len; ++i)
  {
    n_fname[i] = fname[i];
    //     cout << n_fname[i] << " " << i << endl;
  }
  n_fname[len] = '\0';

  //      cout << "n_fname: " << n_fname << " fname: " << fname << ", "
  // 	  << len << endl;

  strcpy(fname, n_fname); /* Send back the old pointer */
  delete n_fname;
}
//|----------------------------------------------------------------------------
//  writefile
//|----------------------------------------------------------------------------
void writefile(int *ext, double *XX, double *YY, double *ZZ, double *datain,
               char *filename, const char *flag)
{
  int nx = ext[0], ny = ext[1], nz = ext[2];
  int i, j, k;
  char filename_h[50];
  //|--->open out put file
  ofstream outfile;

  if (!strcmp(flag, "YZ"))
  {
    for (i = 0; i < nx; i++)
    {
      sprintf(filename_h, "%s_%d.dat", filename, i);
      outfile.open(filename_h);
      outfile << "# CT along X at " << i << endl;
      for (k = 0; k < nz; k++)
      {
        for (j = 0; j < ny; j++)
        {
          outfile << setw(10) << setprecision(10) << YY[j] << " "
                  << setw(10) << setprecision(10) << ZZ[k] << " "
                  << datain[i + j * nx + k * nx * ny] << " "
                  << endl;
        }
        outfile << "\n"; /* blanck line for gnuplot */
      }
      outfile.close();
    }
  }
  else if (!strcmp(flag, "XZ"))
  {
    for (j = 0; j < ny; j++)
    {
      sprintf(filename_h, "%s_%d.dat", filename, j);
      outfile.open(filename_h);
      outfile << "# CT along Y at " << j << endl;
      for (k = 0; k < nz; k++)
      {
        for (i = 0; i < nx; i++)
        {
          outfile << setw(10) << setprecision(10) << XX[i] << " "
                  << setw(10) << setprecision(10) << ZZ[k] << " "
                  << datain[i + j * nx + k * nx * ny] << " "
                  << endl;
        }
        outfile << "\n"; /* blanck line for gnuplot */
      }
      outfile.close();
    }
  }
  else if (!strcmp(flag, "XY"))
  {
    for (k = 0; k < nz; k++)
    {
      sprintf(filename_h, "%s_%d.dat", filename, k);
      outfile.open(filename_h);
      outfile << "# CT along Z at " << k << endl;
      for (j = 0; j < ny; j++)
      {
        for (i = 0; i < nx; i++)
        {
          outfile << setw(10) << setprecision(10) << XX[i] << " "
                  << setw(10) << setprecision(10) << YY[j] << " "
                  << datain[i + j * nx + k * nx * ny] << " "
                  << endl;
        }
        outfile << "\n"; /* blanck line for gnuplot */
      }
      outfile.close();
    }
  }
  else
  {
    cout << "In output_data: not recognized flag-->" << flag << endl;
    exit(0);
  }
}
