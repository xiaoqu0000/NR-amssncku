
#ifdef newc
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cmath>
using namespace std;
#else
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#endif
/* Linear equation solution by Gauss-Jordan elimination.
a[0..n-1][0..n-1] is the input matrix. b[0..n-1] is input
containing the right-hand side vectors. On output a is
replaced by its matrix inverse, and b is replaced by the
corresponding set of solution vectors */

int gaussj(double *a, double *b, int n)
{
  double swap;

  int *indxc, *indxr, *ipiv;
  indxc = new int[n];
  indxr = new int[n];
  ipiv = new int[n];

  int i, icol, irow, j, k, l, ll;
  double big, dum, pivinv, temp;

  for (j = 0; j < n; j++)
    ipiv[j] = 0;
  for (i = 0; i < n; i++)
  {
    big = 0.0;
    for (j = 0; j < n; j++)
      if (ipiv[j] != 1)
        for (k = 0; k < n; k++)
        {
          if (ipiv[k] == 0)
          {
            if (fabs(a[j * n + k]) >= big)
            {
              big = fabs(a[j * n + k]);
              irow = j;
              icol = k;
            }
          }
          else if (ipiv[k] > 1)
          {
            cout << "gaussj: Singular Matrix-1" << endl;
            for (int ii = 0; ii < n; ii++)
            {
              for (int jj = 0; jj < n; jj++)
                cout << a[ii * n + jj] << " ";
              cout << endl;
            }
            return 1; // error return
          }
        }

    ipiv[icol] = ipiv[icol] + 1;
    if (irow != icol)
    {
      for (l = 0; l < n; l++)
      {
        swap = a[irow * n + l];
        a[irow * n + l] = a[icol * n + l];
        a[icol * n + l] = swap;
      }

      swap = b[irow];
      b[irow] = b[icol];
      b[icol] = swap;
    }

    indxr[i] = irow;
    indxc[i] = icol;

    if (a[icol * n + icol] == 0.0)
    {
      cout << "gaussj: Singular Matrix-2" << endl;
      for (int ii = 0; ii < n; ii++)
      {
        for (int jj = 0; jj < n; jj++)
          cout << a[ii * n + jj] << " ";
        cout << endl;
      }
      return 1; // error return
    }

    pivinv = 1.0 / a[icol * n + icol];
    a[icol * n + icol] = 1.0;
    for (l = 0; l < n; l++)
      a[icol * n + l] *= pivinv;
    b[icol] *= pivinv;
    for (ll = 0; ll < n; ll++)
      if (ll != icol)
      {
        dum = a[ll * n + icol];
        a[ll * n + icol] = 0.0;
        for (l = 0; l < n; l++)
          a[ll * n + l] -= a[icol * n + l] * dum;
        b[ll] -= b[icol] * dum;
      }
  }

  for (l = n - 1; l >= 0; l--)
  {
    if (indxr[l] != indxc[l])
      for (k = 0; k < n; k++)
      {
        swap = a[k * n + indxr[l]];
        a[k * n + indxr[l]] = a[k * n + indxc[l]];
        a[k * n + indxc[l]] = swap;
      }
  }

  delete[] indxc;
  delete[] indxr;
  delete[] ipiv;

  return 0;
}
// for check usage
/*
int main()
{
  double *A,*b;
  A=new double[9];
  b=new double[3];

  A[0]=0.5; A[1]=1.0/3; A[2]=1;
  A[3]=1;   A[4]=5.0/3; A[5]=3;
  A[6]=2;   A[7]=4.0/3; A[8]=5;

  b[0]=1; b[1]=3; b[2]=2;

  cout<<"initial data:"<<endl;
  for(int i=0;i<3;i++) cout<<A[i*3]<<" "<<A[i*3+1]<<" "<<A[i*3+2]<<" "<<b[i]<<endl;

  gaussj(A, b, 3);

  cout<<"final data:"<<endl;
  for(int i=0;i<3;i++) cout<<A[i*3]<<" "<<A[i*3+1]<<" "<<A[i*3+2]<<" "<<b[i]<<endl;

  delete[] A; delete[] b;
}
*/
