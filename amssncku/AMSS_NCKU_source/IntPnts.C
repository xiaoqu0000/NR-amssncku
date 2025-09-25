//$Id: IntPnts.C,v 1.1 2012/04/03 10:49:42 zjcao Exp $

#include "macrodef.h"
#ifdef With_AHF

#include <math.h>
#include <stdio.h>

#include <iostream>
using namespace std;

#include "myglobal.h"

namespace AHFinderDirect
{
  extern struct state state;
  int globalInterpGFL(double *X, double *Y, double *Z, int Ns,
                      double *Data)
  {
    if (Ns == 0)
      return 0;
    int n;
    double *pox[3];
    for (int i = 0; i < 3; i++)
      pox[i] = new double[Ns];
    for (n = 0; n < Ns; n++)
    {
      pox[0][n] = X[n];
      pox[1][n] = Y[n];
      pox[2][n] = Z[n];
    }

    const int InList = 35;

    double *datap;
    datap = new double[Ns * InList];
    if (!(state.ADM->AH_Interp_Points(state.AHList, Ns, pox, datap, state.Symmetry)))
      return 0;
    // reform data
    for (int pnt = 0; pnt < Ns; pnt++)
      for (int ii = 0; ii < InList; ii++)
      {
        if (ii == 0 || ii == 12 || ii == 20)
          Data[pnt + ii * Ns] = datap[ii + pnt * InList] + 1;
        else if (ii == 24) // from chi-1 to psi
          Data[pnt + ii * Ns] = pow(datap[ii + pnt * InList] + 1, -0.25);
        else if (ii == 25 || ii == 26 || ii == 27) // from chi,i to psi,i
          Data[pnt + ii * Ns] = -pow(datap[24 + pnt * InList] + 1, -1.25) / 4 * datap[ii + pnt * InList];
        else
          Data[pnt + ii * Ns] = datap[ii + pnt * InList];
      }
    delete[] datap;

    delete[] pox[0];
    delete[] pox[1];
    delete[] pox[2];

    return 1;
  }
  // inerpolate lapse and shift
  int globalInterpGFLlash(double *X, double *Y, double *Z, int Ns,
                          double *Data)
  {
    if (Ns == 0)
      return 0;
    int n;
    double *pox[3];
    for (int i = 0; i < 3; i++)
      pox[i] = new double[Ns];
    for (n = 0; n < Ns; n++)
    {
      pox[0][n] = X[n];
      pox[1][n] = Y[n];
      pox[2][n] = Z[n];
    }

    double SYM = 1.0, ANT = -1.0;
    const int InList = 4;

    double *datap;
    datap = new double[Ns * InList];
    state.ADM->AH_Interp_Points(state.GaugeList, Ns, pox, datap, state.Symmetry);
    // reform data
    for (int pnt = 0; pnt < Ns; pnt++)
      for (int ii = 0; ii < InList; ii++)
        Data[pnt + ii * Ns] = datap[ii + pnt * InList];

    delete[] datap;
    delete[] pox[0];
    delete[] pox[1];
    delete[] pox[2];

    return 1;
  }

} // namespace AHFinderDirect
#endif
