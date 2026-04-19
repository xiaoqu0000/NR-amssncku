//$Id: PN_orbit.C,v 1.1.1.1 2017/09/25 02:35:55 zjcao Exp $
#ifdef newc
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <complex>
using namespace std;
#else
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#endif

#include "misc.h"
#include "lagpolint.h"
#define PI M_PI

//----------------------------------------------------------------------------------------------------
// 该函数设置两体后牛顿方程的右边项
//----------------------------------------------------------------------------------------------------
// PRD 74, 104005 (2006)
// we do not consider spin here
// we have assumed M = 1, G = 1

void compute_rhs(const int N, double *dyn, double *dyn_rhs, double eta)
{
      double qx, qy, qz, px, py, pz;
      qx = dyn[0];
      qy = dyn[1];
      qz = dyn[2];
      px = dyn[3];
      py = dyn[4];
      pz = dyn[5];

      double rqx, rqy, rqz, rpx, rpy, rpz;
#if 1
      // diff(H,y) part
      double MapleGenVar1 = px + (3.0 * eta - 1.0) * (px * px + py * py + pz * pz) * px / 2.0 - (2.0 * (3.0 + eta) * px + 2.0 * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * qx) / sqrt(qx * qx + qy * qy + qz * qz) / 2.0 + 3.0 / 8.0 * (1.0 - 5.0 * eta + 5.0 * eta * eta) * pow(px * px + py * py + pz * pz, 2.0) * px + (4.0 * (5.0 - 20.0 * eta - 3.0 * eta * eta) * (px * px + py * py + pz * pz) * px - 4.0 * eta * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * (px * px + py * py + pz * pz) * qx - 4.0 * eta * eta * pow(qx * px + qy * py + qz * pz, 2.0) / (qx * qx + qy * qy + qz * qz) * px - 12.0 * eta * eta * pow(qx * px + qy * py + qz * pz, 3.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * qx) / sqrt(qx * qx + qy * qy + qz * qz) / 8.0;
      double MapleGenVar2 = MapleGenVar1 + (2.0 * (5.0 + 8.0 * eta) * px + 6.0 * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * qx) / (qx * qx + qy * qy + qz * qz) / 2.0 + (-5.0 + 35.0 * eta - 70.0 * eta * eta + 35.0 * eta * eta * eta) * pow(px * px + py * py + pz * pz, 3.0) * px / 16.0;
      double MapleGenVar3 = MapleGenVar2 + (6.0 * (-7.0 + 42.0 * eta - 53.0 * eta * eta - 5.0 * eta * eta * eta) * pow(px * px + py * py + pz * pz, 2.0) * px + 2.0 * (2.0 - 3.0 * eta) * eta * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * pow(px * px + py * py + pz * pz, 2.0) * qx + 4.0 * (2.0 - 3.0 * eta) * eta * eta * pow(qx * px + qy * py + qz * pz, 2.0) / (qx * qx + qy * qy + qz * qz) * (px * px + py * py + pz * pz) * px + 12.0 * (1.0 - eta) * eta * eta * pow(qx * px + qy * py + qz * pz, 3.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * (px * px + py * py + pz * pz) * qx + 6.0 * (1.0 - eta) * eta * eta * pow(qx * px + qy * py + qz * pz, 4.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * px - 30.0 * eta * eta * eta * pow(qx * px + qy * py + qz * pz, 5.0) / pow(qx * qx + qy * qy + qz * qz, 3.0) * qx) / sqrt(qx * qx + qy * qy + qz * qz) / 16.0;
      rqx = MapleGenVar3 + ((-27.0 + 136.0 * eta + 109.0 * eta * eta) * (px * px + py * py + pz * pz) * px / 4.0 + (17.0 + 30.0 * eta) * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * (px * px + py * py + pz * pz) * qx / 8.0 + (17.0 + 30.0 * eta) * eta * pow(qx * px + qy * py + qz * pz, 2.0) / (qx * qx + qy * qy + qz * qz) * px / 8.0 + (5.0 + 43.0 * eta) * eta * pow(qx * px + qy * py + qz * pz, 3.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * qx / 3.0) / (qx * qx + qy * qy + qz * qz) + (2.0 * (-25.0 / 8.0 + (0.3141592653589793E1 * 0.3141592653589793E1 / 64.0 - 335.0 / 48.0) * eta - 23.0 / 8.0 * eta * eta) * px + 2.0 * (-85.0 / 16.0 - 3.0 / 64.0 * 0.3141592653589793E1 * 0.3141592653589793E1 - 7.0 / 4.0 * eta) * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * qx) / sqrt(pow(qx * qx + qy * qy + qz * qz, 3.0));

      MapleGenVar1 = py + (3.0 * eta - 1.0) * (px * px + py * py + pz * pz) * py / 2.0 - (2.0 * (3.0 + eta) * py + 2.0 * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * qy) / sqrt(qx * qx + qy * qy + qz * qz) / 2.0 + 3.0 / 8.0 * (1.0 - 5.0 * eta + 5.0 * eta * eta) * pow(px * px + py * py + pz * pz, 2.0) * py + (4.0 * (5.0 - 20.0 * eta - 3.0 * eta * eta) * (px * px + py * py + pz * pz) * py - 4.0 * eta * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * (px * px + py * py + pz * pz) * qy - 4.0 * eta * eta * pow(qx * px + qy * py + qz * pz, 2.0) / (qx * qx + qy * qy + qz * qz) * py - 12.0 * eta * eta * pow(qx * px + qy * py + qz * pz, 3.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * qy) / sqrt(qx * qx + qy * qy + qz * qz) / 8.0;
      MapleGenVar2 = MapleGenVar1 + (2.0 * (5.0 + 8.0 * eta) * py + 6.0 * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * qy) / (qx * qx + qy * qy + qz * qz) / 2.0 + (-5.0 + 35.0 * eta - 70.0 * eta * eta + 35.0 * eta * eta * eta) * pow(px * px + py * py + pz * pz, 3.0) * py / 16.0;
      MapleGenVar3 = MapleGenVar2 + (6.0 * (-7.0 + 42.0 * eta - 53.0 * eta * eta - 5.0 * eta * eta * eta) * pow(px * px + py * py + pz * pz, 2.0) * py + 2.0 * (2.0 - 3.0 * eta) * eta * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * pow(px * px + py * py + pz * pz, 2.0) * qy + 4.0 * (2.0 - 3.0 * eta) * eta * eta * pow(qx * px + qy * py + qz * pz, 2.0) / (qx * qx + qy * qy + qz * qz) * (px * px + py * py + pz * pz) * py + 12.0 * (1.0 - eta) * eta * eta * pow(qx * px + qy * py + qz * pz, 3.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * (px * px + py * py + pz * pz) * qy + 6.0 * (1.0 - eta) * eta * eta * pow(qx * px + qy * py + qz * pz, 4.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * py - 30.0 * eta * eta * eta * pow(qx * px + qy * py + qz * pz, 5.0) / pow(qx * qx + qy * qy + qz * qz, 3.0) * qy) / sqrt(qx * qx + qy * qy + qz * qz) / 16.0;
      rqy = MapleGenVar3 + ((-27.0 + 136.0 * eta + 109.0 * eta * eta) * (px * px + py * py + pz * pz) * py / 4.0 + (17.0 + 30.0 * eta) * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * (px * px + py * py + pz * pz) * qy / 8.0 + (17.0 + 30.0 * eta) * eta * pow(qx * px + qy * py + qz * pz, 2.0) / (qx * qx + qy * qy + qz * qz) * py / 8.0 + (5.0 + 43.0 * eta) * eta * pow(qx * px + qy * py + qz * pz, 3.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * qy / 3.0) / (qx * qx + qy * qy + qz * qz) + (2.0 * (-25.0 / 8.0 + (0.3141592653589793E1 * 0.3141592653589793E1 / 64.0 - 335.0 / 48.0) * eta - 23.0 / 8.0 * eta * eta) * py + 2.0 * (-85.0 / 16.0 - 3.0 / 64.0 * 0.3141592653589793E1 * 0.3141592653589793E1 - 7.0 / 4.0 * eta) * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * qy) / sqrt(pow(qx * qx + qy * qy + qz * qz, 3.0));

      MapleGenVar1 = pz + (3.0 * eta - 1.0) * (px * px + py * py + pz * pz) * pz / 2.0 - (2.0 * (3.0 + eta) * pz + 2.0 * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * qz) / sqrt(qx * qx + qy * qy + qz * qz) / 2.0 + 3.0 / 8.0 * (1.0 - 5.0 * eta + 5.0 * eta * eta) * pow(px * px + py * py + pz * pz, 2.0) * pz + (4.0 * (5.0 - 20.0 * eta - 3.0 * eta * eta) * (px * px + py * py + pz * pz) * pz - 4.0 * eta * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * (px * px + py * py + pz * pz) * qz - 4.0 * eta * eta * pow(qx * px + qy * py + qz * pz, 2.0) / (qx * qx + qy * qy + qz * qz) * pz - 12.0 * eta * eta * pow(qx * px + qy * py + qz * pz, 3.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * qz) / sqrt(qx * qx + qy * qy + qz * qz) / 8.0;
      MapleGenVar2 = MapleGenVar1 + (2.0 * (5.0 + 8.0 * eta) * pz + 6.0 * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * qz) / (qx * qx + qy * qy + qz * qz) / 2.0 + (-5.0 + 35.0 * eta - 70.0 * eta * eta + 35.0 * eta * eta * eta) * pow(px * px + py * py + pz * pz, 3.0) * pz / 16.0;
      MapleGenVar3 = MapleGenVar2 + (6.0 * (-7.0 + 42.0 * eta - 53.0 * eta * eta - 5.0 * eta * eta * eta) * pow(px * px + py * py + pz * pz, 2.0) * pz + 2.0 * (2.0 - 3.0 * eta) * eta * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * pow(px * px + py * py + pz * pz, 2.0) * qz + 4.0 * (2.0 - 3.0 * eta) * eta * eta * pow(qx * px + qy * py + qz * pz, 2.0) / (qx * qx + qy * qy + qz * qz) * (px * px + py * py + pz * pz) * pz + 12.0 * (1.0 - eta) * eta * eta * pow(qx * px + qy * py + qz * pz, 3.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * (px * px + py * py + pz * pz) * qz + 6.0 * (1.0 - eta) * eta * eta * pow(qx * px + qy * py + qz * pz, 4.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * pz - 30.0 * eta * eta * eta * pow(qx * px + qy * py + qz * pz, 5.0) / pow(qx * qx + qy * qy + qz * qz, 3.0) * qz) / sqrt(qx * qx + qy * qy + qz * qz) / 16.0;
      rqz = MapleGenVar3 + ((-27.0 + 136.0 * eta + 109.0 * eta * eta) * (px * px + py * py + pz * pz) * pz / 4.0 + (17.0 + 30.0 * eta) * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * (px * px + py * py + pz * pz) * qz / 8.0 + (17.0 + 30.0 * eta) * eta * pow(qx * px + qy * py + qz * pz, 2.0) / (qx * qx + qy * qy + qz * qz) * pz / 8.0 + (5.0 + 43.0 * eta) * eta * pow(qx * px + qy * py + qz * pz, 3.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * qz / 3.0) / (qx * qx + qy * qy + qz * qz) + (2.0 * (-25.0 / 8.0 + (0.3141592653589793E1 * 0.3141592653589793E1 / 64.0 - 335.0 / 48.0) * eta - 23.0 / 8.0 * eta * eta) * pz + 2.0 * (-85.0 / 16.0 - 3.0 / 64.0 * 0.3141592653589793E1 * 0.3141592653589793E1 - 7.0 / 4.0 * eta) * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * qz) / sqrt(pow(qx * qx + qy * qy + qz * qz, 3.0));

      MapleGenVar2 = 1 / (sqrt(pow(qx * qx + qy * qy + qz * qz, 3.0))) * qx - (2.0 * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * px - 2.0 * eta * pow(qx * px + qy * py + qz * pz, 2.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * qx) / sqrt(qx * qx + qy * qy + qz * qz) / 2.0 + ((3.0 + eta) * (px * px + py * py + pz * pz) + eta * pow(qx * px + qy * py + qz * pz, 2.0) / (qx * qx + qy * qy + qz * qz)) / sqrt(pow(qx * qx + qy * qy + qz * qz, 3.0)) * qx / 2.0 - 1 / (pow(qx * qx + qy * qy + qz * qz, 2.0)) * qx;
      MapleGenVar3 = MapleGenVar2 + (-4.0 * eta * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * (px * px + py * py + pz * pz) * px + 4.0 * eta * eta * pow(qx * px + qy * py + qz * pz, 2.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * (px * px + py * py + pz * pz) * qx - 12.0 * eta * eta * pow(qx * px + qy * py + qz * pz, 3.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * px + 12.0 * eta * eta * pow(qx * px + qy * py + qz * pz, 4.0) / pow(qx * qx + qy * qy + qz * qz, 3.0) * qx) / sqrt(qx * qx + qy * qy + qz * qz) / 8.0;
      MapleGenVar1 = MapleGenVar3 - ((5.0 - 20.0 * eta - 3.0 * eta * eta) * pow(px * px + py * py + pz * pz, 2.0) - 2.0 * eta * eta * pow(qx * px + qy * py + qz * pz, 2.0) / (qx * qx + qy * qy + qz * qz) * (px * px + py * py + pz * pz) - 3.0 * eta * eta * pow(qx * px + qy * py + qz * pz, 4.0) / pow(qx * qx + qy * qy + qz * qz, 2.0)) / sqrt(pow(qx * qx + qy * qy + qz * qz, 3.0)) * qx / 8.0 + (6.0 * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * px - 6.0 * eta * pow(qx * px + qy * py + qz * pz, 2.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * qx) / (qx * qx + qy * qy + qz * qz) / 2.0 - ((5.0 + 8.0 * eta) * (px * px + py * py + pz * pz) + 3.0 * eta * pow(qx * px + qy * py + qz * pz, 2.0) / (qx * qx + qy * qy + qz * qz)) / pow(qx * qx + qy * qy + qz * qz, 2.0) * qx;
      MapleGenVar3 = MapleGenVar1 + 3.0 / 4.0 * (1.0 + 3.0 * eta) / sqrt(pow(qx * qx + qy * qy + qz * qz, 5.0)) * qx;
      double MapleGenVar4 = MapleGenVar3;
      double MapleGenVar6 = (2.0 * (2.0 - 3.0 * eta) * eta * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * pow(px * px + py * py + pz * pz, 2.0) * px - 2.0 * (2.0 - 3.0 * eta) * eta * eta * pow(qx * px + qy * py + qz * pz, 2.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * pow(px * px + py * py + pz * pz, 2.0) * qx + 12.0 * (1.0 - eta) * eta * eta * pow(qx * px + qy * py + qz * pz, 3.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * (px * px + py * py + pz * pz) * px - 12.0 * (1.0 - eta) * eta * eta * pow(qx * px + qy * py + qz * pz, 4.0) / pow(qx * qx + qy * qy + qz * qz, 3.0) * (px * px + py * py + pz * pz) * qx - 30.0 * eta * eta * eta * pow(qx * px + qy * py + qz * pz, 5.0) / pow(qx * qx + qy * qy + qz * qz, 3.0) * px + 30.0 * eta * eta * eta * pow(qx * px + qy * py + qz * pz, 6.0) / pow(qx * qx + qy * qy + qz * qz, 4.0) * qx) / sqrt(qx * qx + qy * qy + qz * qz) / 16.0;
      double MapleGenVar7 = -((-7.0 + 42.0 * eta - 53.0 * eta * eta - 5.0 * eta * eta * eta) * pow(px * px +
                                                                                                       py * py + pz * pz,
                                                                                                   3.0) +
                              (2.0 - 3.0 * eta) * eta * eta * pow(qx * px + qy * py + qz * pz, 2.0) / (qx * qx + qy * qy + qz * qz) * pow(px * px + py * py + pz * pz, 2.0) + 3.0 * (1.0 - eta) * eta * eta * pow(qx * px + qy * py + qz * pz, 4.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * (px * px + py * py + pz * pz) - 5.0 * eta * eta * eta * pow(qx * px + qy * py + qz * pz, 6.0) / pow(qx * qx + qy * qy + qz * qz, 3.0)) /
                            sqrt(pow(qx * qx + qy * qy + qz * qz, 3.0)) *
                            qx / 16.0;
      double MapleGenVar5 = MapleGenVar6 + MapleGenVar7;
      MapleGenVar2 = MapleGenVar4 + MapleGenVar5;
      MapleGenVar3 = MapleGenVar2 + ((17.0 + 30.0 * eta) * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * (px * px + py * py + pz * pz) * px / 8.0 - (17.0 + 30.0 * eta) * eta * pow(qx * px + qy * py + qz * pz, 2.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * (px * px + py * py + pz * pz) * qx / 8.0 + (5.0 + 43.0 * eta) * eta * pow(qx * px + qy * py + qz * pz, 3.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * px / 3.0 - (5.0 + 43.0 * eta) * eta * pow(qx * px + qy * py + qz * pz, 4.0) / pow(qx * qx + qy * qy + qz * qz, 3.0) * qx / 3.0) / (qx * qx + qy * qy + qz * qz) - 2.0 * ((-27.0 + 136.0 * eta + 109.0 * eta * eta) * pow(px * px + py * py + pz * pz, 2.0) / 16.0 + (17.0 + 30.0 * eta) * eta * pow(qx * px + qy * py + qz * pz, 2.0) / (qx * qx + qy * qy + qz * qz) * (px * px + py * py + pz * pz) / 16.0 + (5.0 + 43.0 * eta) * eta * pow(qx * px + qy * py + qz * pz, 4.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) / 12.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * qx;
      rpx = MapleGenVar3 + (2.0 * (-85.0 / 16.0 - 3.0 / 64.0 * 0.3141592653589793E1 * 0.3141592653589793E1 - 7.0 / 4.0 * eta) * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * px - 2.0 * (-85.0 / 16.0 - 3.0 / 64.0 * 0.3141592653589793E1 * 0.3141592653589793E1 - 7.0 / 4.0 * eta) * eta * pow(qx * px + qy * py + qz * pz, 2.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * qx) / sqrt(pow(qx * qx + qy * qy + qz * qz, 3.0)) - 3.0 * ((-25.0 / 8.0 + (0.3141592653589793E1 * 0.3141592653589793E1 / 64.0 - 335.0 / 48.0) * eta - 23.0 / 8.0 * eta * eta) * (px * px + py * py + pz * pz) + (-85.0 / 16.0 - 3.0 / 64.0 * 0.3141592653589793E1 * 0.3141592653589793E1 - 7.0 / 4.0 * eta) * eta * pow(qx * px + qy * py + qz * pz, 2.0) / (qx * qx + qy * qy + qz * qz)) / sqrt(pow(qx * qx + qy * qy + qz * qz, 5.0)) * qx - 4.0 * (1.0 / 8.0 + (109.0 / 12.0 - 21.0 / 32.0 * 0.3141592653589793E1 * 0.3141592653589793E1) * eta) / pow(qx * qx + qy * qy + qz * qz, 3.0) * qx;

      MapleGenVar2 = 1 / (sqrt(pow(qx * qx + qy * qy + qz * qz, 3.0))) * qy - (2.0 * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * py - 2.0 * eta * pow(qx * px + qy * py + qz * pz, 2.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * qy) / sqrt(qx * qx + qy * qy + qz * qz) / 2.0 + ((3.0 + eta) * (px * px + py * py + pz * pz) + eta * pow(qx * px + qy * py + qz * pz, 2.0) / (qx * qx + qy * qy + qz * qz)) / sqrt(pow(qx * qx + qy * qy + qz * qz, 3.0)) * qy / 2.0 - 1 / (pow(qx * qx + qy * qy + qz * qz, 2.0)) * qy;
      MapleGenVar3 = MapleGenVar2 + (-4.0 * eta * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * (px * px + py * py + pz * pz) * py + 4.0 * eta * eta * pow(qx * px + qy * py + qz * pz, 2.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * (px * px + py * py + pz * pz) * qy - 12.0 * eta * eta * pow(qx * px + qy * py + qz * pz, 3.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * py + 12.0 * eta * eta * pow(qx * px + qy * py + qz * pz, 4.0) / pow(qx * qx + qy * qy + qz * qz, 3.0) * qy) / sqrt(qx * qx + qy * qy + qz * qz) / 8.0;
      MapleGenVar1 = MapleGenVar3 - ((5.0 - 20.0 * eta - 3.0 * eta * eta) * pow(px * px + py * py + pz * pz, 2.0) - 2.0 * eta * eta * pow(qx * px + qy * py + qz * pz, 2.0) / (qx * qx + qy * qy + qz * qz) * (px * px + py * py + pz * pz) - 3.0 * eta * eta * pow(qx * px + qy * py + qz * pz, 4.0) / pow(qx * qx + qy * qy + qz * qz, 2.0)) / sqrt(pow(qx * qx + qy * qy + qz * qz, 3.0)) * qy / 8.0 + (6.0 * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * py - 6.0 * eta * pow(qx * px + qy * py + qz * pz, 2.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * qy) / (qx * qx + qy * qy + qz * qz) / 2.0 - ((5.0 + 8.0 * eta) * (px * px + py * py + pz * pz) + 3.0 * eta * pow(qx * px + qy * py + qz * pz, 2.0) / (qx * qx + qy * qy + qz * qz)) / pow(qx * qx + qy * qy + qz * qz, 2.0) * qy;
      MapleGenVar3 = MapleGenVar1 + 3.0 / 4.0 * (1.0 + 3.0 * eta) / sqrt(pow(qx * qx + qy * qy + qz * qz, 5.0)) * qy;
      MapleGenVar4 = MapleGenVar3;
      MapleGenVar6 = (2.0 * (2.0 - 3.0 * eta) * eta * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * pow(px * px + py * py + pz * pz, 2.0) * py - 2.0 * (2.0 - 3.0 * eta) * eta * eta * pow(qx * px + qy * py + qz * pz, 2.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * pow(px * px + py * py + pz * pz, 2.0) * qy + 12.0 * (1.0 - eta) * eta * eta * pow(qx * px + qy * py + qz * pz, 3.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * (px * px + py * py + pz * pz) * py - 12.0 * (1.0 - eta) * eta * eta * pow(qx * px + qy * py + qz * pz, 4.0) / pow(qx * qx + qy * qy + qz * qz, 3.0) * (px * px + py * py + pz * pz) * qy - 30.0 * eta * eta * eta * pow(qx * px + qy * py + qz * pz, 5.0) / pow(qx * qx + qy * qy + qz * qz, 3.0) * py + 30.0 * eta * eta * eta * pow(qx * px + qy * py + qz * pz, 6.0) / pow(qx * qx + qy * qy + qz * qz, 4.0) * qy) / sqrt(qx * qx + qy * qy + qz * qz) / 16.0;
      MapleGenVar7 = -((-7.0 + 42.0 * eta - 53.0 * eta * eta - 5.0 * eta * eta * eta) * pow(px * px +
                                                                                                py * py + pz * pz,
                                                                                            3.0) +
                       (2.0 - 3.0 * eta) * eta * eta * pow(qx * px + qy * py + qz * pz, 2.0) / (qx * qx + qy * qy + qz * qz) * pow(px * px + py * py + pz * pz, 2.0) + 3.0 * (1.0 - eta) * eta * eta * pow(qx * px + qy * py + qz * pz, 4.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * (px * px + py * py + pz * pz) - 5.0 * eta * eta * eta * pow(qx * px + qy * py + qz * pz, 6.0) / pow(qx * qx + qy * qy + qz * qz, 3.0)) /
                     sqrt(pow(qx * qx + qy * qy + qz * qz, 3.0)) *
                     qy / 16.0;
      MapleGenVar5 = MapleGenVar6 + MapleGenVar7;
      MapleGenVar2 = MapleGenVar4 + MapleGenVar5;
      MapleGenVar3 = MapleGenVar2 + ((17.0 + 30.0 * eta) * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * (px * px + py * py + pz * pz) * py / 8.0 - (17.0 + 30.0 * eta) * eta * pow(qx * px + qy * py + qz * pz, 2.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * (px * px + py * py + pz * pz) * qy / 8.0 + (5.0 + 43.0 * eta) * eta * pow(qx * px + qy * py + qz * pz, 3.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * py / 3.0 - (5.0 + 43.0 * eta) * eta * pow(qx * px + qy * py + qz * pz, 4.0) / pow(qx * qx + qy * qy + qz * qz, 3.0) * qy / 3.0) / (qx * qx + qy * qy + qz * qz) - 2.0 * ((-27.0 + 136.0 * eta + 109.0 * eta * eta) * pow(px * px + py * py + pz * pz, 2.0) / 16.0 + (17.0 + 30.0 * eta) * eta * pow(qx * px + qy * py + qz * pz, 2.0) / (qx * qx + qy * qy + qz * qz) * (px * px + py * py + pz * pz) / 16.0 + (5.0 + 43.0 * eta) * eta * pow(qx * px + qy * py + qz * pz, 4.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) / 12.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * qy;
      rpy = MapleGenVar3 + (2.0 * (-85.0 / 16.0 - 3.0 / 64.0 * 0.3141592653589793E1 * 0.3141592653589793E1 - 7.0 / 4.0 * eta) * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * py - 2.0 * (-85.0 / 16.0 - 3.0 / 64.0 * 0.3141592653589793E1 * 0.3141592653589793E1 - 7.0 / 4.0 * eta) * eta * pow(qx * px + qy * py + qz * pz, 2.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * qy) / sqrt(pow(qx * qx + qy * qy + qz * qz, 3.0)) - 3.0 * ((-25.0 / 8.0 + (0.3141592653589793E1 * 0.3141592653589793E1 / 64.0 - 335.0 / 48.0) * eta - 23.0 / 8.0 * eta * eta) * (px * px + py * py + pz * pz) + (-85.0 / 16.0 - 3.0 / 64.0 * 0.3141592653589793E1 * 0.3141592653589793E1 - 7.0 / 4.0 * eta) * eta * pow(qx * px + qy * py + qz * pz, 2.0) / (qx * qx + qy * qy + qz * qz)) / sqrt(pow(qx * qx + qy * qy + qz * qz, 5.0)) * qy - 4.0 * (1.0 / 8.0 + (109.0 / 12.0 - 21.0 / 32.0 * 0.3141592653589793E1 * 0.3141592653589793E1) * eta) / pow(qx * qx + qy * qy + qz * qz, 3.0) * qy;

      MapleGenVar2 = 1 / (sqrt(pow(qx * qx + qy * qy + qz * qz, 3.0))) * qz - (2.0 * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * pz - 2.0 * eta * pow(qx * px + qy * py + qz * pz, 2.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * qz) / sqrt(qx * qx + qy * qy + qz * qz) / 2.0 + ((3.0 + eta) * (px * px + py * py + pz * pz) + eta * pow(qx * px + qy * py + qz * pz, 2.0) / (qx * qx + qy * qy + qz * qz)) / sqrt(pow(qx * qx + qy * qy + qz * qz, 3.0)) * qz / 2.0 - 1 / (pow(qx * qx + qy * qy + qz * qz, 2.0)) * qz;
      MapleGenVar3 = MapleGenVar2 + (-4.0 * eta * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * (px * px + py * py + pz * pz) * pz + 4.0 * eta * eta * pow(qx * px + qy * py + qz * pz, 2.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * (px * px + py * py + pz * pz) * qz - 12.0 * eta * eta * pow(qx * px + qy * py + qz * pz, 3.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * pz + 12.0 * eta * eta * pow(qx * px + qy * py + qz * pz, 4.0) / pow(qx * qx + qy * qy + qz * qz, 3.0) * qz) / sqrt(qx * qx + qy * qy + qz * qz) / 8.0;
      MapleGenVar1 = MapleGenVar3 - ((5.0 - 20.0 * eta - 3.0 * eta * eta) * pow(px * px + py * py + pz * pz, 2.0) - 2.0 * eta * eta * pow(qx * px + qy * py + qz * pz, 2.0) / (qx * qx + qy * qy + qz * qz) * (px * px + py * py + pz * pz) - 3.0 * eta * eta * pow(qx * px + qy * py + qz * pz, 4.0) / pow(qx * qx + qy * qy + qz * qz, 2.0)) / sqrt(pow(qx * qx + qy * qy + qz * qz, 3.0)) * qz / 8.0 + (6.0 * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * pz - 6.0 * eta * pow(qx * px + qy * py + qz * pz, 2.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * qz) / (qx * qx + qy * qy + qz * qz) / 2.0 - ((5.0 + 8.0 * eta) * (px * px + py * py + pz * pz) + 3.0 * eta * pow(qx * px + qy * py + qz * pz, 2.0) / (qx * qx + qy * qy + qz * qz)) / pow(qx * qx + qy * qy + qz * qz, 2.0) * qz;
      MapleGenVar3 = MapleGenVar1 + 3.0 / 4.0 * (1.0 + 3.0 * eta) / sqrt(pow(qx * qx + qy * qy + qz * qz, 5.0)) * qz;
      MapleGenVar4 = MapleGenVar3;
      MapleGenVar6 = (2.0 * (2.0 - 3.0 * eta) * eta * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * pow(px * px + py * py + pz * pz, 2.0) * pz - 2.0 * (2.0 - 3.0 * eta) * eta * eta * pow(qx * px + qy * py + qz * pz, 2.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * pow(px * px + py * py + pz * pz, 2.0) * qz + 12.0 * (1.0 - eta) * eta * eta * pow(qx * px + qy * py + qz * pz, 3.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * (px * px + py * py + pz * pz) * pz - 12.0 * (1.0 - eta) * eta * eta * pow(qx * px + qy * py + qz * pz, 4.0) / pow(qx * qx + qy * qy + qz * qz, 3.0) * (px * px + py * py + pz * pz) * qz - 30.0 * eta * eta * eta * pow(qx * px + qy * py + qz * pz, 5.0) / pow(qx * qx + qy * qy + qz * qz, 3.0) * pz + 30.0 * eta * eta * eta * pow(qx * px + qy * py + qz * pz, 6.0) / pow(qx * qx + qy * qy + qz * qz, 4.0) * qz) / sqrt(qx * qx + qy * qy + qz * qz) / 16.0;
      MapleGenVar7 = -((-7.0 + 42.0 * eta - 53.0 * eta * eta - 5.0 * eta * eta * eta) * pow(px * px +
                                                                                                py * py + pz * pz,
                                                                                            3.0) +
                       (2.0 - 3.0 * eta) * eta * eta * pow(qx * px + qy * py + qz * pz, 2.0) / (qx * qx + qy * qy + qz * qz) * pow(px * px + py * py + pz * pz, 2.0) + 3.0 * (1.0 - eta) * eta * eta * pow(qx * px + qy * py + qz * pz, 4.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * (px * px + py * py + pz * pz) - 5.0 * eta * eta * eta * pow(qx * px + qy * py + qz * pz, 6.0) / pow(qx * qx + qy * qy + qz * qz, 3.0)) /
                     sqrt(pow(qx * qx + qy * qy + qz * qz, 3.0)) *
                     qz / 16.0;
      MapleGenVar5 = MapleGenVar6 + MapleGenVar7;
      MapleGenVar2 = MapleGenVar4 + MapleGenVar5;
      MapleGenVar3 = MapleGenVar2 + ((17.0 + 30.0 * eta) * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * (px * px + py * py + pz * pz) * pz / 8.0 - (17.0 + 30.0 * eta) * eta * pow(qx * px + qy * py + qz * pz, 2.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * (px * px + py * py + pz * pz) * qz / 8.0 + (5.0 + 43.0 * eta) * eta * pow(qx * px + qy * py + qz * pz, 3.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * pz / 3.0 - (5.0 + 43.0 * eta) * eta * pow(qx * px + qy * py + qz * pz, 4.0) / pow(qx * qx + qy * qy + qz * qz, 3.0) * qz / 3.0) / (qx * qx + qy * qy + qz * qz) - 2.0 * ((-27.0 + 136.0 * eta + 109.0 * eta * eta) * pow(px * px + py * py + pz * pz, 2.0) / 16.0 + (17.0 + 30.0 * eta) * eta * pow(qx * px + qy * py + qz * pz, 2.0) / (qx * qx + qy * qy + qz * qz) * (px * px + py * py + pz * pz) / 16.0 + (5.0 + 43.0 * eta) * eta * pow(qx * px + qy * py + qz * pz, 4.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) / 12.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * qz;
      rpz = MapleGenVar3 + (2.0 * (-85.0 / 16.0 - 3.0 / 64.0 * 0.3141592653589793E1 * 0.3141592653589793E1 - 7.0 / 4.0 * eta) * eta * (qx * px + qy * py + qz * pz) / (qx * qx + qy * qy + qz * qz) * pz - 2.0 * (-85.0 / 16.0 - 3.0 / 64.0 * 0.3141592653589793E1 * 0.3141592653589793E1 - 7.0 / 4.0 * eta) * eta * pow(qx * px + qy * py + qz * pz, 2.0) / pow(qx * qx + qy * qy + qz * qz, 2.0) * qz) / sqrt(pow(qx * qx + qy * qy + qz * qz, 3.0)) - 3.0 * ((-25.0 / 8.0 + (0.3141592653589793E1 * 0.3141592653589793E1 / 64.0 - 335.0 / 48.0) * eta - 23.0 / 8.0 * eta * eta) * (px * px + py * py + pz * pz) + (-85.0 / 16.0 - 3.0 / 64.0 * 0.3141592653589793E1 * 0.3141592653589793E1 - 7.0 / 4.0 * eta) * eta * pow(qx * px + qy * py + qz * pz, 2.0) / (qx * qx + qy * qy + qz * qz)) / sqrt(pow(qx * qx + qy * qy + qz * qz, 5.0)) * qz - 4.0 * (1.0 / 8.0 + (109.0 / 12.0 - 21.0 / 32.0 * 0.3141592653589793E1 * 0.3141592653589793E1) * eta) / pow(qx * qx + qy * qy + qz * qz, 3.0) * qz;
#else
      {
            double mu = eta;
            double ata = eta;

            double pp = px * px + py * py + pz * pz;
            double sq = sqrt(qx * qx + qy * qy + qz * qz);
            double nx = qx / sq;
            double ny = qy / sq;
            double nz = qz / sq;
            double np = px * nx + py * ny + pz * nz;

            double coe_p_newton = 1.0;
            double coe_n_newton = 0.0;

            double coe_p_1pn = 0.5 * (3.0 * ata - 1.0) * pp - (3.0 + ata) / sq;
            double coe_n_1pn = -1.0 * ata * np / sq;

            double coe_p_2pn = 0.375 * (1 - 5.0 * ata + 5.0 * ata * ata) * pp * pp + 0.5 / sq * (5.0 - 20.0 * ata - 3.0 * ata * ata) * pp - 0.5 * ata * ata * np * np / sq + (5.0 + 8.0 * ata) / pow(sq, 2);
            double coe_n_2pn = -0.5 / sq * ata * ata * np * pp - 1.5 / sq * ata * ata * pow(np, 3) + 3.0 * ata * np / (sq * sq);

            double coe_p_3pn = 0.0625 * (-5.0 + 35.0 * ata - 70.0 * ata * ata + 35.0 * pow(ata, 3)) * pow(pp, 3) + 0.375 / sq * (-7.0 + 42.0 * ata - 53.0 * pow(ata, 2) - 5.0 * pow(ata, 3)) * pow(pp, 2) + 0.25 / sq * (2.0 - 3.0 * ata) * ata * ata * np * np * pp + 0.375 / sq * (1.0 - ata) * ata * ata * pow(np, 4) + 0.25 / (sq * sq) * (-27.0 + 136.0 * ata + 109.0 * ata * ata) * pp + 0.125 / (sq * sq) * (17.0 + 30.0 * ata) * ata * np * np + 2.0 / pow(sq, 3) * (-3.125 + (0.015625 * PI * PI - 335.0 / 48.0) * ata - 2.875 * pow(ata, 2));

            double coe_n_3pn = 0.125 / sq * (2.0 - 3.0 * ata) * ata * ata * pp * pp * np + 0.75 / sq * (1.0 - ata) * ata * ata * pow(np, 3) * pp - 1.875 / sq * pow(ata, 3) * pow(np, 5) + 0.125 / (sq * sq) * (17.0 + 30.0 * ata) * ata * np * pp + 1.0 / (3.0 * sq * sq) * (5.0 + 43.0 * ata) * ata * pow(np, 3) + 2.0 * ata / pow(sq, 3) * (-5.3125 - 0.046875 * PI * PI - 1.75 * ata) * np;

            double coe_p = coe_p_newton + coe_p_1pn + coe_p_2pn + coe_p_3pn;
            double coe_n = coe_n_newton + coe_n_1pn + coe_n_2pn + coe_n_3pn;

            rqx = coe_p * px + coe_n * nx;
            rqy = coe_p * py + coe_n * ny;
            rqz = coe_p * pz + coe_n * nz;

            double q_crs_px = qy * pz - qz * py;
            double q_crs_py = qz * px - qx * pz;
            double q_crs_pz = qx * py - qy * px;

            double bx = q_crs_py * qz - q_crs_pz * qy;
            double by = q_crs_pz * qx - q_crs_px * qz;
            double bz = q_crs_px * qy - q_crs_py * qx;

            bx = bx / pow(sq, 3);
            by = by / pow(sq, 3);
            bz = bz / pow(sq, 3);

            coe_n_newton = 1.0 / pow(sq, 2);
            double coe_b_newton = 0.0;

            coe_n_1pn = 0.5 / (sq * sq) * ((3 + ata) * pp + ata * np * np) - 1.0 / pow(sq, 3);
            double coe_b_1pn = -1.0 / sq * ata * np;

            coe_n_2pn = -0.125 / (sq * sq) * ((5.0 - 20.0 * ata - 3.0 * ata * ata) * pp * pp - 2.0 * ata * ata * np * np * pp - 3.0 * ata * ata * pow(np, 4)) - 1.0 / pow(sq, 3) * ((5.0 + 8.0 * ata) * pp + 3.0 * ata * np * np) + 0.75 / pow(sq, 4) * (1.0 + 3.0 * ata);

            double coe_b_2pn = -0.5 / sq * (ata * ata * np * pp + 3.0 * ata * ata * pow(np, 3)) + 3.0 / (sq * sq) * ata * np;

            coe_n_3pn = 0.0625 * ((-7.0 + 42.0 * ata - 53.0 * pow(ata, 2) - 5.0 * pow(ata, 3)) * pow(pp, 3) + (2.0 - 3.0 * ata) * ata * ata * np * np * pp * pp + 3.0 * (1.0 - ata) * ata * ata * pow(np, 4) * pp - 5.0 * pow(ata, 3) * pow(np, 6)) * (-1.0 / pow(sq, 2)) + (0.0625 * (-27.0 + 136.0 * ata + 109.0 * ata * ata) * pp * pp + 0.0625 * (17.0 + 30.0 * ata) * ata * np * np * pp + 1.0 / 12.0 * (5.0 + 43.0 * ata) * ata * pow(np, 4)) * (-2.0 / pow(sq, 3)) + ((-3.125 + (0.015625 * PI * PI - 335.0 / 48.0) * ata - 2.875 * ata * ata) * pp + (-5.3125 - 0.046875 * PI * PI - 1.75 * ata) * ata * np * np) * (-3.0 / pow(sq, 4)) + (0.125 + (109.0 / 12.0 - 0.65625 * PI * PI) * ata) * (-4.0 / pow(sq, 5));

            double coe_b_3pn = 0.0625 / sq * ((2.0 - 3.0 * ata) * ata * ata * pp * pp * 2.0 * np + 12.0 * (1.0 - ata) * ata * ata * pow(np, 3) * pp - 30.0 * pow(ata, 3) * pow(np, 5)) + 1.0 / (sq * sq) * (0.125 * (17.0 + 30.0 * ata) * ata * np * pp + 1.0 / 3.0 * (5.0 + 43.0 * ata) * ata * pow(np, 3)) + 2.0 / pow(sq, 3) * (-5.3125 - 0.046875 * PI * PI - 1.75 * ata) * ata * np;

            coe_n = coe_n_newton + coe_n_1pn + coe_n_2pn + coe_n_3pn;
            double coe_b = coe_b_newton + coe_b_1pn + coe_b_2pn + coe_b_3pn;

            rpx = coe_n * nx + coe_b * bx;
            rpy = coe_n * ny + coe_b * by;
            rpz = coe_n * nz + coe_b * bz;
      }
#endif

      // add F
      double vx, vy, vz, ome, vome;
      // ome = |v x r|/r^2
      vx = rqy * qz - rqz * qy;
      vy = rqz * qx - rqx * qz;
      vz = rqx * qy - rqy * qx;

      ome = sqrt(vx * vx + vy * vy + vz * vz) / (qx * qx + qy * qy + qz * qz);
      // |L| = |q x P| = nu * |q x p|
      double L;
      vx = qy * pz - qz * py;
      vy = qz * px - qx * pz;
      vz = qx * py - qy * px;
      L = eta * sqrt(vx * vx + vy * vy + vz * vz);

      vome = pow(ome, 1.0 / 3);

      double dedt;

      // PRD 74, 104005 (2006), Eq.(3.15)
      double f2, f3, f4, f5, f6, fl6, f7;
      const double gammae = 0.5772156649;
      f2 = -1247 / 336.0 - 35 / 12.0 * eta;
      f3 = 4 * PI;
      f4 = -44711 / 9072.0 + 9271 / 504.0 * eta + 65 / 18.0 * eta * eta;
      f5 = -(8191 / 672.0 + 583 / 24.0 * eta) * PI;
      f6 = 6643739519 / 69854400.0 + 16 / 3.0 * PI * PI - 1712 / 105.0 * gammae + (-134543 / 7776.0 + 41 / 48.0 * PI * PI) * eta - 94403 / 3024.0 * eta * eta - 775 / 324.0 * eta * eta * eta;
      fl6 = -1712 / 105.0;
      f7 = (-16285 / 504.0 + 214745 / 1728.0 * eta + 193385 / 3024.0 * eta * eta) * PI;

      dedt = -32 / 5.0 * eta * eta * pow(vome, 10) * (1 + f2 * vome * vome + f3 * pow(vome, 3) + f4 * pow(vome, 4) + f5 * pow(vome, 5) + f6 * pow(vome, 6) + fl6 * pow(vome, 6) * log(4 * vome) + f7 * pow(vome, 7));

      // 根据论文 PRD 74, 104005 (2006) 中 (3.14)
      double ome_B;
      ome_B = dedt / ome / L;

      double rpx_rhs, rpy_rhs, rpz_rhs;
      rpx_rhs = -rpx + ome_B * px;
      rpy_rhs = -rpy + ome_B * py;
      rpz_rhs = -rpz + ome_B * pz;

      dyn_rhs[0] = rqx;
      dyn_rhs[1] = rqy;
      dyn_rhs[2] = rqz;
      dyn_rhs[3] = rpx_rhs;
      dyn_rhs[4] = rpy_rhs;
      dyn_rhs[5] = rpz_rhs;
}
//----------------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------------------
// 该函数设定 PN_Orbit 的输入
//----------------------------------------------------------------------------------------------------
void setid(double &x, double &y, double &z, double &px, double &py, double &pz, double nu)
{
      double r;
      // read parameter from file
      {
            char filename[50];
            strcpy(filename, "PN_Orbit.par");
            const int LEN = 256;
            char pline[LEN];
            string str, sgrp, skey, sval;
            int sind;
            ifstream inf(filename, ifstream::in);
            if (!inf.good())
            {
                  cout << "Can not open parameter file " << filename << " for inputing information of black holes" << endl;
                  exit(0);
            }

            for (int i = 1; inf.good(); i++)
            {
                  inf.getline(pline, LEN);
                  str = pline;

                  int status = misc::parse_parts(str, sgrp, skey, sval, sind);
                  if (status == -1)
                  {
                        cout << "error reading parameter file " << filename << " in line " << i << endl;
                        exit(0);
                  }
                  else if (status == 0)
                        continue;

                  if (sgrp == "PN" && skey == "initial coordinate separation")
                        r = atof(sval.c_str());
            }
            inf.close();
      }

      x = r;
      y = 0;
      z = 0;

      // PRD 77,024027 (2008) Eq.(65)
      double p = sqrt(1 / r);
      p = p + 2 * pow(p, 3) + 1 / 16.0 * (42 - 43 * nu) * pow(p, 5) + 1 / 128.0 * (480 + (163 * PI * PI - 4556) * nu + 104 * nu * nu) * pow(p, 7);

      px = 0;
      py = p;
      pz = 0;
}
//----------------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------------------
// 该函数设定切向动量与径向动量
//----------------------------------------------------------------------------------------------------
void trdecompose(double x, double y, double z,
                 double px, double py, double pz,
                 double &pr, double &pt)
{
      double r = sqrt(x * x + y * y + z * z);

      double vx, vy, vz;
      vx = y * pz - z * py;
      vy = z * px - x * pz;
      vz = x * py - y * px;
      pt = sqrt(vx * vx + vy * vy + vz * vz) / r;
      pr = (px * x + py * y + pz * z) / r;
}
//----------------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------------------
// push several elements backward in an array
void pushback(int N, int m, double *a)
{
      for (int i = 0; i + m < N; i++)
            a[i] = a[i + m];
}
//----------------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------------------
// 此为程序的主函数
//----------------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{

      double nu;
      double dT;
      double wr;
      const int ord = 6;

      // 读取参数文件
      // read parameter from file
      {
            char filename[50];
            strcpy(filename, "PN_Orbit.par");
            const int LEN = 256;
            char pline[LEN];
            string str, sgrp, skey, sval;
            int sind;
            ifstream inf(filename, ifstream::in);
            if (!inf.good())
            {
                  cout << "Can not open parameter file " << filename << " for inputing information of black holes" << endl;
                  exit(0);
            }

            for (int i = 1; inf.good(); i++)
            {
                  inf.getline(pline, LEN);
                  str = pline;

                  int status = misc::parse_parts(str, sgrp, skey, sval, sind);
                  if (status == -1)
                  {
                        cout << "error reading parameter file " << filename << " in line " << i << endl;
                        exit(0);
                  }
                  else if (status == 0)
                        continue;

                  if (sgrp == "PN" && skey == "symmetric mass ratio")
                        nu = atof(sval.c_str());
                  else if (sgrp == "PN" && skey == "dT")
                        dT = atof(sval.c_str());
                  else if (sgrp == "PN" && skey == "wanted r")
                        wr = atof(sval.c_str());
            }
            inf.close();
      }

      //   initial data
      // 根据参数文件设定 PN_Orbit 的初值

      cout << "M = 1; nu = " << nu << endl;

      int N = 6;
      double *dyn0, *dyn1, *dyn_rhs;
      double *dyn;
      dyn0 = new double[N];
      dyn1 = new double[N];
      dyn_rhs = new double[N];
      dyn = new double[N];

      // 利用 setid 函数设定 PN_Orbit 的输入
      setid(dyn0[0], dyn0[1], dyn0[2], dyn0[3], dyn0[4], dyn0[5], nu);

      // 屏幕输出初始的位置和动量
      cout << "initial data are: x y z px py pz" << endl;
      cout << "        x  = " << dyn0[0] << endl
           << "        y  = " << dyn0[1] << endl
           << "        z  = " << dyn0[2] << endl
           << "       px  = " << dyn0[3] << endl
           << "       py  = " << dyn0[4] << endl
           << "       pz  = " << dyn0[5] << endl
           << "       Px  = " << dyn0[3] * nu << endl
           << "       Py  = " << dyn0[4] * nu << endl
           << "       Pz  = " << dyn0[5] * nu << endl;

      ofstream outfile;
      outfile.open("PNorbit.dat");
      outfile.setf(ios::scientific, ios::floatfield);
      outfile.precision(16);
      double TotalTime = 1000000, dumptime = 5;
      int Steps = int(TotalTime / dT), dump_every = int(dumptime / dT);
      if (dump_every < 1)
            dump_every = 1;

      cout << "evolve with time step dT = " << dT << endl;
      outfile << "# time  x y z px py pz r" << endl;

      double RR[ord], PT[ord], PR[ord];

      // 下面使用龙格库塔法计算每步迭代过程中的位置和动量

      for (int ncount = 0; ncount < Steps; ncount++)
      {
            int iter_count = 0;

            // 用龙格库塔法计算下一步的轨道位置和动量

            // Predictor
            compute_rhs(N, dyn0, dyn_rhs, nu);
            misc::rungekutta4(N, dT, dyn0, dyn, dyn_rhs, iter_count);
            // Corrector
            for (iter_count = 1; iter_count < 4; iter_count++)
            {
                  compute_rhs(N, dyn, dyn1, nu);
                  misc::rungekutta4(N, dT, dyn0, dyn1, dyn_rhs, iter_count);
                  if (iter_count < 3)
                  {
                        double *tt = dyn;
                        dyn = dyn1;
                        dyn1 = tt;
                  }
                  else
                  {
                        double *tt = dyn0;
                        dyn0 = dyn1;
                        dyn1 = tt;
                  }
            }

            // 根据轨道位置和动量，计算出双星系统的切向动量和径向动量
            double pt, pr;
            trdecompose(dyn0[0], dyn0[1], dyn0[2], dyn0[3], dyn0[4], dyn0[5], pr, pt);

            if (ncount < ord)
            {
                  RR[ncount] = sqrt(dyn0[0] * dyn0[0] + dyn0[1] * dyn0[1] + dyn0[2] * dyn0[2]);
                  PT[ncount] = pt;
                  PR[ncount] = pr;
            }
            else
            {
                  pushback(ord, 1, RR);
                  pushback(ord, 1, PT);
                  pushback(ord, 1, PR);

                  RR[ord - 1] = sqrt(dyn0[0] * dyn0[0] + dyn0[1] * dyn0[1] + dyn0[2] * dyn0[2]);
                  PT[ord - 1] = pt;
                  PR[ord - 1] = pr;

                  // 如果半径 RR 达到设定值 wr，则输出最后的结果
                  if ((RR[ord / 2] - wr) * (RR[ord / 2 - 1] - wr) < 0)
                  {
                        double dy;
                        int nord = ord;
                        f_polint(RR, PT, wr, pt, dy, nord);
                        f_polint(RR, PR, wr, pr, dy, nord);

                        cout << "succussfully find at r = " << wr << endl;
                        cout << "                    pt = " << pt << endl;
                        cout << "                    pr = " << pr << endl;
                        cout << "                    Pt = nu*pt = " << nu * pt << endl;
                        cout << "                    Pr = nu*pr = " << nu * pr << endl;
                        exit(1);
                  }
            }

            // 每经过 dump_every 步，向输出文件中写入数据
            if (ncount / dump_every * dump_every == ncount)
            {

                  outfile << ncount * dT;
                  for (int i = 0; i < N; i++)
                        outfile << " " << dyn0[i];
                  outfile << " " << sqrt(dyn0[0] * dyn0[0] + dyn0[1] * dyn0[1] + dyn0[2] * dyn0[2]);
                  outfile << " " << nu * pr << " " << nu * pt;
                  outfile << endl;
            }

            // 检查是否算出了 NaN
            bool ffgg = finite(dyn0[0]);
            for (int i = 1; i < N; i++)
                  ffgg = (ffgg && (finite(dyn0[i])));
            if (!ffgg)
            {
                  cout << "find NaN at t = " << ncount * dT << endl;
                  cout << "        x  = " << dyn0[0] << endl
                       << "        y  = " << dyn0[1] << endl
                       << "        z  = " << dyn0[2] << endl
                       << "       px  = " << dyn0[3] << endl
                       << "       py  = " << dyn0[4] << endl
                       << "       pz  = " << dyn0[5] << endl;
                  break;
            }
            
            // 如果达到最大迭代步数，输出相关信息
            if (ncount == Steps - 1)
            {
                  ncount = ord;
                  cout << Steps << " steps have been done" << endl;
            }
      }
      outfile.close();

      delete[] dyn0, delete[] dyn1;
      delete[] dyn_rhs;
      delete[] dyn;
}
//----------------------------------------------------------------------------------------------------
// 程序主函数结束
//----------------------------------------------------------------------------------------------------
