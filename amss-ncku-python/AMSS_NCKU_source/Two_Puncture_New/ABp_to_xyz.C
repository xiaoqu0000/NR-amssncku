
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

#include "TwoPunctures.h"

#include "../For_Debug/tillherecheck.h"

//|----------------------------------------------------------------------------
//  write ASCII file with the style of Pablo
//|----------------------------------------------------------------------------
void TwoPunctures::ABp_to_xyz(const double A, const double B, const double phi,
                              double &x, double &y, double &z)
{
  double At = 0.5 * (A + 1);
  double X = 2 * atanh(At);
  double R = Pih + 2 * atan(B);

  complex<double> C = complex<double>(X, R);
  complex<double> c = cosh(C) * par_b; /* c=b*cosh(C)*/

  x = real(c);
  double r = imag(c);

  double sin_phi = sin(phi);
  double cos_phi = cos(phi);

  y = r * cos_phi;
  z = r * sin_phi;

  return;
}
