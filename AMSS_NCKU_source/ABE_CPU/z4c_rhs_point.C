

// Z4c rhs without advection term
#ifdef newc
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
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

#include "macrodef.fh"

#define Power(x, y) (pow((double)(x), (double)(y)))
#define Sqrt(x) sqrt(x)
#define Log(x) log((double)(x))
#define pow2(x) ((x) * (x))
#define pow3(x) ((x) * (x) * (x))
#define pow4(x) ((x) * (x) * (x) * (x))
#define pow2inv(x) (1.0 / ((x) * (x)))

#define Cal(x, y, z) ((x) ? (y) : (z))

#define Tan(x) tan(x)
#define ArcTan(x) atan(x)
#define Sin(x) sin(x)
#define Cos(x) cos(x)
#define Csc(x) (1. / sin(x))
#define Abs(x) (fabs(x))
#define sqrt2 (sqrt(2))
#define Tanh(x) tanh(x)
#define Sech(x) (1 / cosh(x))

extern "C"
{

#ifdef fortran1
  void z4c_rhs_point
#endif
#ifdef fortran2
      void Z4C_RHS_POINT
#endif
#ifdef fortran3
      void
      z4c_rhs_point_
#endif
      (double &A11,
       double &A12,
       double &A13,
       double &A22,
       double &A23,
       double &A33,
       double &alpha,
       double &B1,
       double &B2,
       double &B3,
       double &beta1,
       double &beta2,
       double &beta3,
       double &chi,
       double &chiDivFloor,
       double &da1,
       double &dA111,
       double &dA112,
       double &dA113,
       double &dA122,
       double &dA123,
       double &dA133,
       double &da2,
       double &dA211,
       double &dA212,
       double &dA213,
       double &dA222,
       double &dA223,
       double &dA233,
       double &da3,
       double &dA311,
       double &dA312,
       double &dA313,
       double &dA322,
       double &dA323,
       double &dA333,
       double &db11,
       double &dB11,
       double &db12,
       double &dB12,
       double &db13,
       double &dB13,
       double &db21,
       double &dB21,
       double &db22,
       double &dB22,
       double &db23,
       double &dB23,
       double &db31,
       double &dB31,
       double &db32,
       double &dB32,
       double &db33,
       double &dB33,
       double &dchi1,
       double &dchi2,
       double &dchi3,
       double &dda11,
       double &dda12,
       double &dda13,
       double &dda22,
       double &dda23,
       double &dda33,
       double &ddb111,
       double &ddb112,
       double &ddb113,
       double &ddb121,
       double &ddb122,
       double &ddb123,
       double &ddb131,
       double &ddb132,
       double &ddb133,
       double &ddb221,
       double &ddb222,
       double &ddb223,
       double &ddb231,
       double &ddb232,
       double &ddb233,
       double &ddb331,
       double &ddb332,
       double &ddb333,
       double &ddchi11,
       double &ddchi12,
       double &ddchi13,
       double &ddchi22,
       double &ddchi23,
       double &ddchi33,
       double &deldelg1111,
       double &deldelg1112,
       double &deldelg1113,
       double &deldelg1122,
       double &deldelg1123,
       double &deldelg1133,
       double &deldelg1211,
       double &deldelg1212,
       double &deldelg1213,
       double &deldelg1222,
       double &deldelg1223,
       double &deldelg1233,
       double &deldelg1311,
       double &deldelg1312,
       double &deldelg1313,
       double &deldelg1322,
       double &deldelg1323,
       double &deldelg1333,
       double &deldelg2211,
       double &deldelg2212,
       double &deldelg2213,
       double &deldelg2222,
       double &deldelg2223,
       double &deldelg2233,
       double &deldelg2311,
       double &deldelg2312,
       double &deldelg2313,
       double &deldelg2322,
       double &deldelg2323,
       double &deldelg2333,
       double &deldelg3311,
       double &deldelg3312,
       double &deldelg3313,
       double &deldelg3322,
       double &deldelg3323,
       double &deldelg3333,
       double &delG11,
       double &delg111,
       double &delg112,
       double &delg113,
       double &delG12,
       double &delg122,
       double &delg123,
       double &delG13,
       double &delg133,
       double &delG21,
       double &delg211,
       double &delg212,
       double &delg213,
       double &delG22,
       double &delg222,
       double &delg223,
       double &delG23,
       double &delg233,
       double &delG31,
       double &delg311,
       double &delg312,
       double &delg313,
       double &delG32,
       double &delg322,
       double &delg323,
       double &delG33,
       double &delg333,
       double &dKhat1,
       double &dKhat2,
       double &dKhat3,
       double &dTheta1,
       double &dTheta2,
       double &dTheta3,
       double &G1,
       double &g11,
       double &g12,
       double &g13,
       double &G2,
       double &g22,
       double &g23,
       double &G3,
       double &g33,
       double &kappa1,
       double &kappa2,
       double &Khat,
       double &rA11,
       double &rA12,
       double &rA13,
       double &rA22,
       double &rA23,
       double &rA33,
       double &rchi,
       double &rG1,
       double &rg11,
       double &rg12,
       double &rg13,
       double &rG2,
       double &rg22,
       double &rg23,
       double &rG3,
       double &rg33,
       double &rKhat,
       double &rTheta,
       double &Theta)
  {

    double AA11;
    double AA12;
    double AA13;
    double AA22;
    double AA23;
    double AA33;
    double Ainv11;
    double Ainv12;
    double Ainv13;
    double Ainv22;
    double Ainv23;
    double Ainv33;
    double cAA;
    double cdda11;
    double cdda12;
    double cdda13;
    double cdda22;
    double cdda23;
    double cdda33;
    double cddf11;
    double cddf12;
    double cddf13;
    double cddf22;
    double cddf23;
    double cddf33;
    double chiguard;
    double chiguarded;
    double chipsipower;
    double ddf11;
    double ddf12;
    double ddf13;
    double ddf22;
    double ddf23;
    double ddf33;
    double detginv;
    double df1;
    double df2;
    double df3;
    double dGd11;
    double dGd12;
    double dGd13;
    double dGd21;
    double dGd22;
    double dGd23;
    double dGd31;
    double dGd32;
    double dGd33;
    double dginv111;
    double dginv112;
    double dginv113;
    double dginv122;
    double dginv123;
    double dginv133;
    double dginv211;
    double dginv212;
    double dginv213;
    double dginv222;
    double dginv223;
    double dginv233;
    double dginv311;
    double dginv312;
    double dginv313;
    double dginv322;
    double dginv323;
    double dginv333;
    double divAinv1;
    double divAinv2;
    double divAinv3;
    double divbeta;
    double dK1;
    double dK2;
    double dK3;
    double dphi1;
    double dphi2;
    double dphi3;
    double dZ11;
    double DZ11;
    double dZ12;
    double DZ12;
    double dZ13;
    double DZ13;
    double dZ21;
    double DZ21;
    double dZ22;
    double DZ22;
    double dZ23;
    double DZ23;
    double dZ31;
    double DZ31;
    double dZ32;
    double DZ32;
    double dZ33;
    double DZ33;
    double dZinv11;
    double DZinv11;
    double dZinv12;
    double DZinv12;
    double dZinv13;
    double DZinv13;
    double dZinv21;
    double DZinv21;
    double dZinv22;
    double DZinv22;
    double dZinv23;
    double DZinv23;
    double dZinv31;
    double DZinv31;
    double dZinv32;
    double DZinv32;
    double dZinv33;
    double DZinv33;
    double DZsym11;
    double DZsym12;
    double DZsym13;
    double DZsym21;
    double DZsym22;
    double DZsym23;
    double DZsym31;
    double DZsym32;
    double DZsym33;
    double f;
    double ff;
    double gamma111;
    double gamma112;
    double gamma113;
    double gamma122;
    double gamma123;
    double gamma133;
    double gamma211;
    double gamma212;
    double gamma213;
    double gamma222;
    double gamma223;
    double gamma233;
    double gamma311;
    double gamma312;
    double gamma313;
    double gamma322;
    double gamma323;
    double gamma333;
    double gammado111;
    double gammado112;
    double gammado113;
    double gammado122;
    double gammado123;
    double gammado133;
    double gammado211;
    double gammado212;
    double gammado213;
    double gammado222;
    double gammado223;
    double gammado233;
    double gammado311;
    double gammado312;
    double gammado313;
    double gammado322;
    double gammado323;
    double gammado333;
    double gammaF111;
    double gammaF112;
    double gammaF113;
    double gammaF121;
    double gammaF122;
    double gammaF123;
    double gammaF131;
    double gammaF132;
    double gammaF133;
    double gammaF211;
    double gammaF212;
    double gammaF213;
    double gammaF221;
    double gammaF222;
    double gammaF223;
    double gammaF231;
    double gammaF232;
    double gammaF233;
    double gammaF311;
    double gammaF312;
    double gammaF313;
    double gammaF321;
    double gammaF322;
    double gammaF323;
    double gammaF331;
    double gammaF332;
    double gammaF333;
    double Gd1;
    double Gd2;
    double Gd3;
    double Gfromg1;
    double Gfromg2;
    double Gfromg3;
    double ginv11;
    double ginv12;
    double ginv13;
    double ginv22;
    double ginv23;
    double ginv33;
    double Hhat;
    double K;
    double lieA11;
    double lieA12;
    double lieA13;
    double lieA22;
    double lieA23;
    double lieA33;
    double liechi;
    double lieg11;
    double lieg12;
    double lieg13;
    double lieg22;
    double lieg23;
    double lieg33;
    double oochipsipower;
    double ootddivbeta1;
    double ootddivbeta2;
    double ootddivbeta3;
    double pseudolieG1;
    double pseudolieG2;
    double pseudolieG3;
    double psim4;
    double R11;
    double R12;
    double R13;
    double R22;
    double R23;
    double R33;
    double Rhat;
    double Rphi11;
    double Rphi12;
    double Rphi13;
    double Rphi22;
    double Rphi23;
    double Rphi33;
    double totdivbeta;
    double trcdda;
    double trcddf;
    double trDZsym;
    double Z1;
    double Z2;
    double Z3;
    double Zinv1;
    double Zinv2;
    double Zinv3;

    chipsipower =
        -4.;

    K =
        Khat + 2. * Theta;

    dK1 =
        dKhat1 + 2. * dTheta1;

    dK2 =
        dKhat2 + 2. * dTheta2;

    dK3 =
        dKhat3 + 2. * dTheta3;

    detginv =
        1 / (2. * g12 * g13 * g23 - g33 * pow2(g12) + g22 * (g11 * g33 - pow2(g13)) -
             g11 * pow2(g23));

    ginv11 =
        detginv * (g22 * g33 - pow2(g23));

    ginv12 =
        detginv * (g13 * g23 - g12 * g33);

    ginv13 =
        detginv * (-(g13 * g22) + g12 * g23);

    ginv22 =
        detginv * (g11 * g33 - pow2(g13));

    ginv23 =
        detginv * (g12 * g13 - g11 * g23);

    ginv33 =
        detginv * (g11 * g22 - pow2(g12));

    dginv111 =
        -2. * (delg123 * ginv12 * ginv13 + ginv11 * (delg112 * ginv12 + delg113 * ginv13)) -
        delg111 * pow2(ginv11) - delg122 * pow2(ginv12) - delg133 * pow2(ginv13);

    dginv112 =
        -(ginv11 * (delg111 * ginv12 + delg112 * ginv22 + delg113 * ginv23)) -
        ginv12 * (delg113 * ginv13 + delg122 * ginv22 + delg123 * ginv23) -
        ginv13 * (delg123 * ginv22 + delg133 * ginv23) - delg112 * pow2(ginv12);

    dginv113 =
        -(ginv11 * (delg111 * ginv13 + delg112 * ginv23 + delg113 * ginv33)) -
        ginv12 * (delg112 * ginv13 + delg122 * ginv23 + delg123 * ginv33) -
        ginv13 * (delg123 * ginv23 + delg133 * ginv33) - delg113 * pow2(ginv13);

    dginv122 =
        -2. * (delg123 * ginv22 * ginv23 + ginv12 * (delg112 * ginv22 + delg113 * ginv23)) -
        delg111 * pow2(ginv12) - delg122 * pow2(ginv22) - delg133 * pow2(ginv23);

    dginv123 =
        -(ginv13 * (delg112 * ginv22 + delg113 * ginv23)) - delg133 * ginv23 * ginv33 -
        ginv12 * (delg111 * ginv13 + delg112 * ginv23 + delg113 * ginv33) -
        ginv22 * (delg122 * ginv23 + delg123 * ginv33) - delg123 * pow2(ginv23);

    dginv133 =
        -2. * (delg123 * ginv23 * ginv33 + ginv13 * (delg112 * ginv23 + delg113 * ginv33)) -
        delg111 * pow2(ginv13) - delg122 * pow2(ginv23) - delg133 * pow2(ginv33);

    dginv211 =
        -2. * (delg223 * ginv12 * ginv13 + ginv11 * (delg212 * ginv12 + delg213 * ginv13)) -
        delg211 * pow2(ginv11) - delg222 * pow2(ginv12) - delg233 * pow2(ginv13);

    dginv212 =
        -(ginv11 * (delg211 * ginv12 + delg212 * ginv22 + delg213 * ginv23)) -
        ginv12 * (delg213 * ginv13 + delg222 * ginv22 + delg223 * ginv23) -
        ginv13 * (delg223 * ginv22 + delg233 * ginv23) - delg212 * pow2(ginv12);

    dginv213 =
        -(ginv11 * (delg211 * ginv13 + delg212 * ginv23 + delg213 * ginv33)) -
        ginv12 * (delg212 * ginv13 + delg222 * ginv23 + delg223 * ginv33) -
        ginv13 * (delg223 * ginv23 + delg233 * ginv33) - delg213 * pow2(ginv13);

    dginv222 =
        -2. * (delg223 * ginv22 * ginv23 + ginv12 * (delg212 * ginv22 + delg213 * ginv23)) -
        delg211 * pow2(ginv12) - delg222 * pow2(ginv22) - delg233 * pow2(ginv23);

    dginv223 =
        -(ginv13 * (delg212 * ginv22 + delg213 * ginv23)) - delg233 * ginv23 * ginv33 -
        ginv12 * (delg211 * ginv13 + delg212 * ginv23 + delg213 * ginv33) -
        ginv22 * (delg222 * ginv23 + delg223 * ginv33) - delg223 * pow2(ginv23);

    dginv233 =
        -2. * (delg223 * ginv23 * ginv33 + ginv13 * (delg212 * ginv23 + delg213 * ginv33)) -
        delg211 * pow2(ginv13) - delg222 * pow2(ginv23) - delg233 * pow2(ginv33);

    dginv311 =
        -2. * (delg323 * ginv12 * ginv13 + ginv11 * (delg312 * ginv12 + delg313 * ginv13)) -
        delg311 * pow2(ginv11) - delg322 * pow2(ginv12) - delg333 * pow2(ginv13);

    dginv312 =
        -(ginv11 * (delg311 * ginv12 + delg312 * ginv22 + delg313 * ginv23)) -
        ginv12 * (delg313 * ginv13 + delg322 * ginv22 + delg323 * ginv23) -
        ginv13 * (delg323 * ginv22 + delg333 * ginv23) - delg312 * pow2(ginv12);

    dginv313 =
        -(ginv11 * (delg311 * ginv13 + delg312 * ginv23 + delg313 * ginv33)) -
        ginv12 * (delg312 * ginv13 + delg322 * ginv23 + delg323 * ginv33) -
        ginv13 * (delg323 * ginv23 + delg333 * ginv33) - delg313 * pow2(ginv13);

    dginv322 =
        -2. * (delg323 * ginv22 * ginv23 + ginv12 * (delg312 * ginv22 + delg313 * ginv23)) -
        delg311 * pow2(ginv12) - delg322 * pow2(ginv22) - delg333 * pow2(ginv23);

    dginv323 =
        -(ginv13 * (delg312 * ginv22 + delg313 * ginv23)) - delg333 * ginv23 * ginv33 -
        ginv12 * (delg311 * ginv13 + delg312 * ginv23 + delg313 * ginv33) -
        ginv22 * (delg322 * ginv23 + delg323 * ginv33) - delg323 * pow2(ginv23);

    dginv333 =
        -2. * (delg323 * ginv23 * ginv33 + ginv13 * (delg312 * ginv23 + delg313 * ginv33)) -
        delg311 * pow2(ginv13) - delg322 * pow2(ginv23) - delg333 * pow2(ginv33);

    gammado111 =
        0.5 * delg111;

    gammado112 =
        0.5 * delg211;

    gammado113 =
        0.5 * delg311;

    gammado122 =
        -0.5 * delg122 + delg212;

    gammado123 =
        0.5 * (-delg123 + delg213 + delg312);

    gammado133 =
        -0.5 * delg133 + delg313;

    gammado211 =
        delg112 - 0.5 * delg211;

    gammado212 =
        0.5 * delg122;

    gammado213 =
        0.5 * (delg123 - delg213 + delg312);

    gammado222 =
        0.5 * delg222;

    gammado223 =
        0.5 * delg322;

    gammado233 =
        -0.5 * delg233 + delg323;

    gammado311 =
        delg113 - 0.5 * delg311;

    gammado312 =
        0.5 * (delg123 + delg213 - delg312);

    gammado313 =
        0.5 * delg133;

    gammado322 =
        delg223 - 0.5 * delg322;

    gammado323 =
        0.5 * delg233;

    gammado333 =
        0.5 * delg333;

    gamma111 =
        gammado111 * ginv11 + gammado211 * ginv12 + gammado311 * ginv13;

    gamma112 =
        gammado112 * ginv11 + gammado212 * ginv12 + gammado312 * ginv13;

    gamma113 =
        gammado113 * ginv11 + gammado213 * ginv12 + gammado313 * ginv13;

    gamma122 =
        gammado122 * ginv11 + gammado222 * ginv12 + gammado322 * ginv13;

    gamma123 =
        gammado123 * ginv11 + gammado223 * ginv12 + gammado323 * ginv13;

    gamma133 =
        gammado133 * ginv11 + gammado233 * ginv12 + gammado333 * ginv13;

    gamma211 =
        gammado111 * ginv12 + gammado211 * ginv22 + gammado311 * ginv23;

    gamma212 =
        gammado112 * ginv12 + gammado212 * ginv22 + gammado312 * ginv23;

    gamma213 =
        gammado113 * ginv12 + gammado213 * ginv22 + gammado313 * ginv23;

    gamma222 =
        gammado122 * ginv12 + gammado222 * ginv22 + gammado322 * ginv23;

    gamma223 =
        gammado123 * ginv12 + gammado223 * ginv22 + gammado323 * ginv23;

    gamma233 =
        gammado133 * ginv12 + gammado233 * ginv22 + gammado333 * ginv23;

    gamma311 =
        gammado111 * ginv13 + gammado211 * ginv23 + gammado311 * ginv33;

    gamma312 =
        gammado112 * ginv13 + gammado212 * ginv23 + gammado312 * ginv33;

    gamma313 =
        gammado113 * ginv13 + gammado213 * ginv23 + gammado313 * ginv33;

    gamma322 =
        gammado122 * ginv13 + gammado222 * ginv23 + gammado322 * ginv33;

    gamma323 =
        gammado123 * ginv13 + gammado223 * ginv23 + gammado323 * ginv33;

    gamma333 =
        gammado133 * ginv13 + gammado233 * ginv23 + gammado333 * ginv33;

    Gfromg1 =
        gamma111 * ginv11 + gamma122 * ginv22 +
        2. * (gamma112 * ginv12 + gamma113 * ginv13 + gamma123 * ginv23) + gamma133 * ginv33;

    Gfromg2 =
        gamma211 * ginv11 + gamma222 * ginv22 +
        2. * (gamma212 * ginv12 + gamma213 * ginv13 + gamma223 * ginv23) + gamma233 * ginv33;

    Gfromg3 =
        gamma311 * ginv11 + gamma322 * ginv22 +
        2. * (gamma312 * ginv12 + gamma313 * ginv13 + gamma323 * ginv23) + gamma333 * ginv33;

    R11 =
        delG11 * g11 + delG12 * g12 + delG13 * g13 + gammado111 * Gfromg1 +
        gammado112 * Gfromg2 + gammado113 * Gfromg3 +
        (-0.5 * deldelg1111 + 3. * gamma111 * gammado111 +
         2. * (gamma211 * gammado112 + gamma311 * gammado113) +
         gamma211 * gammado211 + gamma311 * gammado311) *
            ginv11 +
        (-deldelg1211 + 3. * (gamma112 * gammado111 + gamma111 * gammado112) +
         2. * (gamma212 * gammado112 + gamma312 * gammado113 +
               gamma211 * gammado122 + gamma311 * gammado123) +
         gamma212 * gammado211 +
         gamma211 * gammado212 + gamma312 * gammado311 + gamma311 * gammado312) *
            ginv12 +
        (-deldelg1311 + 3. * (gamma113 * gammado111 + gamma111 * gammado113) +
         2. * (gamma213 * gammado112 + gamma313 * gammado113 +
               gamma211 * gammado123 + gamma311 * gammado133) +
         gamma213 * gammado211 +
         gamma211 * gammado213 + gamma313 * gammado311 + gamma311 * gammado313) *
            ginv13 +
        (-0.5 * deldelg2211 + 3. * gamma112 * gammado112 +
         2. * (gamma212 * gammado122 + gamma312 * gammado123) +
         gamma212 * gammado212 + gamma312 * gammado312) *
            ginv22 +
        (-deldelg2311 + 3. * (gamma113 * gammado112 + gamma112 * gammado113) +
         2. * (gamma213 * gammado122 + (gamma212 + gamma313) * gammado123 +
               gamma312 * gammado133) +
         gamma213 * gammado212 + gamma212 * gammado213 +
         gamma313 * gammado312 + gamma312 * gammado313) *
            ginv23 +
        (-0.5 * deldelg3311 + 3. * gamma113 * gammado113 +
         2. * (gamma213 * gammado123 + gamma313 * gammado133) + gamma213 * gammado213 +
         gamma313 * gammado313) *
            ginv33;

    R12 =
        0.5 * (delG21 * g11 + (delG11 + delG22) * g12 + delG23 * g13 + delG12 * g22 +
               delG13 * g23 + (gammado112 + gammado211) * Gfromg1 +
               (gammado122 + gammado212) * Gfromg2 + (gammado123 + gammado213) * Gfromg3) +
        (-0.5 * deldelg1112 + gamma112 * gammado111 +
         (gamma111 + gamma212) * gammado112 + gamma312 * gammado113 +
         gamma111 * gammado211 + 2. * gamma211 * gammado212 +
         gamma311 * (gammado213 + gammado312)) *
            ginv11 +
        (-deldelg1212 + gamma122 * gammado111 +
         (2. * gamma112 + gamma222) * gammado112 + gamma322 * gammado113 +
         (gamma111 + gamma212) * gammado122 + gamma112 * gammado211 +
         (gamma111 + 2. * gamma212) * gammado212 + 2. * gamma211 * gammado222 +
         gamma312 * (gammado123 + gammado213 + gammado312) +
         gamma311 * (gammado223 + gammado322)) *
            ginv12 +
        (-deldelg1312 + gamma123 * gammado111 + (gamma113 + gamma223) * gammado112 +
         (gamma112 + gamma323) * gammado113 + (gamma111 + gamma212) * gammado123 +
         gamma312 * gammado133 + gamma113 * gammado211 +
         (gamma111 + gamma313) * gammado213 +
         2. * (gamma213 * gammado212 + gamma211 * gammado223) +
         gamma313 * gammado312 + gamma311 * (gammado233 + gammado323)) *
            ginv13 +
        (-0.5 * deldelg2212 + gamma122 * gammado112 +
         (gamma112 + gamma222) * gammado122 + gamma322 * gammado123 +
         gamma112 * gammado212 + 2. * gamma212 * gammado222 +
         gamma312 * (gammado223 + gammado322)) *
            ginv22 +
        (-deldelg2312 + gamma123 * gammado112 + gamma122 * gammado113 +
         (gamma113 + gamma223) * gammado122 +
         (gamma112 + gamma222 + gamma323) * gammado123 + gamma322 * gammado133 +
         gamma113 * gammado212 + gamma112 * gammado213 +
         2. * (gamma213 * gammado222 + gamma212 * gammado223) +
         gamma313 * (gammado223 + gammado322) +
         gamma312 * (gammado233 + gammado323)) *
            ginv23 +
        (-0.5 * deldelg3312 + gamma123 * gammado113 +
         (gamma113 + gamma223) * gammado123 + gamma323 * gammado133 +
         gamma113 * gammado213 + 2. * gamma213 * gammado223 +
         gamma313 * (gammado233 + gammado323)) *
            ginv33;

    R13 =
        0.5 * (delG31 * g11 + delG32 * g12 + (delG11 + delG33) * g13 + delG12 * g23 +
               delG13 * g33 + (gammado113 + gammado311) * Gfromg1 +
               (gammado123 + gammado312) * Gfromg2 + (gammado133 + gammado313) * Gfromg3) +
        (-0.5 * deldelg1113 + gamma113 * gammado111 + gamma213 * gammado112 +
         (gamma111 + gamma313) * gammado113 + gamma111 * gammado311 +
         gamma211 * (gammado213 + gammado312) + 2. * gamma311 * gammado313) *
            ginv11 +
        (-deldelg1213 + gamma123 * gammado111 + (gamma113 + gamma223) * gammado112 +
         (gamma112 + gamma323) * gammado113 + gamma213 * gammado122 +
         (gamma111 + gamma313) * gammado123 + gamma112 * gammado311 +
         gamma111 * gammado312 + gamma212 * (gammado213 + gammado312) +
         gamma211 * (gammado223 + gammado322) +
         2. * (gamma312 * gammado313 + gamma311 * gammado323)) *
            ginv12 +
        (-deldelg1313 + gamma133 * gammado111 + gamma233 * gammado112 +
         (2. * gamma113 + gamma333) * gammado113 +
         (gamma111 + gamma313) * gammado133 + gamma113 * gammado311 +
         gamma213 * (gammado123 + gammado213 + gammado312) +
         (gamma111 + 2. * gamma313) * gammado313 +
         gamma211 * (gammado233 + gammado323) + 2. * gamma311 * gammado333) *
            ginv13 +
        (-0.5 * deldelg2213 + gamma123 * gammado112 + gamma223 * gammado122 +
         (gamma112 + gamma323) * gammado123 + gamma112 * gammado312 +
         gamma212 * (gammado223 + gammado322) + 2. * gamma312 * gammado323) *
            ginv22 +
        (-deldelg2313 + gamma133 * gammado112 + gamma123 * gammado113 +
         gamma233 * gammado122 + (gamma113 + gamma223 + gamma333) * gammado123 +
         (gamma112 + gamma323) * gammado133 + gamma113 * gammado312 +
         gamma112 * gammado313 + gamma213 * (gammado223 + gammado322) +
         gamma212 * (gammado233 + gammado323) +
         2. * (gamma313 * gammado323 + gamma312 * gammado333)) *
            ginv23 +
        (-0.5 * deldelg3313 + gamma133 * gammado113 + gamma233 * gammado123 +
         (gamma113 + gamma333) * gammado133 + gamma113 * gammado313 +
         gamma213 * (gammado233 + gammado323) + 2. * gamma313 * gammado333) *
            ginv33;

    R22 =
        delG21 * g12 + delG22 * g22 + delG23 * g23 + gammado212 * Gfromg1 +
        gammado222 * Gfromg2 + gammado223 * Gfromg3 +
        (-0.5 * deldelg1122 + gamma112 * (gammado112 + 2. * gammado211) +
         3. * gamma212 * gammado212 + gamma312 * (2. * gammado213 + gammado312)) *
            ginv11 +
        (-deldelg1222 + gamma122 * (gammado112 + 2. * gammado211) +
         gamma112 * (gammado122 + 2. * gammado212) +
         3. * (gamma222 * gammado212 + gamma212 * gammado222) +
         2. * (gamma322 * gammado213 + gamma312 * gammado223) +
         gamma322 * gammado312 + gamma312 * gammado322) *
            ginv12 +
        (-deldelg1322 + gamma123 * (gammado112 + 2. * gammado211) +
         gamma112 * (gammado123 + 2. * gammado213) +
         3. * (gamma223 * gammado212 + gamma212 * gammado223) +
         2. * (gamma323 * gammado213 + gamma312 * gammado233) +
         gamma323 * gammado312 + gamma312 * gammado323) *
            ginv13 +
        (-0.5 * deldelg2222 + gamma122 * (gammado122 + 2. * gammado212) +
         3. * gamma222 * gammado222 + gamma322 * (2. * gammado223 + gammado322)) *
            ginv22 +
        (-deldelg2322 + gamma123 * (gammado122 + 2. * gammado212) +
         gamma122 * (gammado123 + 2. * gammado213) +
         3. * (gamma223 * gammado222 + gamma222 * gammado223) +
         2. * (gamma323 * gammado223 + gamma322 * gammado233) +
         gamma323 * gammado322 + gamma322 * gammado323) *
            ginv23 +
        (-0.5 * deldelg3322 + gamma123 * (gammado123 + 2. * gammado213) +
         3. * gamma223 * gammado223 + gamma323 * (2. * gammado233 + gammado323)) *
            ginv33;

    R23 =
        0.5 * (delG31 * g12 + delG21 * g13 + delG32 * g22 + (delG22 + delG33) * g23 +
               delG23 * g33 + (gammado213 + gammado312) * Gfromg1 +
               (gammado223 + gammado322) * Gfromg2 + (gammado233 + gammado323) * Gfromg3) +
        (-0.5 * deldelg1123 + gamma113 * gammado211 + gamma213 * gammado212 +
         (gamma212 + gamma313) * gammado213 +
         gamma112 * (gammado113 + gammado311) + gamma212 * gammado312 +
         2. * gamma312 * gammado313) *
            ginv11 +
        (-deldelg1223 + gamma123 * gammado211 + (gamma113 + gamma223) * gammado212 +
         (gamma222 + gamma323) * gammado213 + gamma213 * gammado222 +
         (gamma212 + gamma313) * gammado223 +
         gamma122 * (gammado113 + gammado311) + gamma222 * gammado312 +
         gamma112 * (gammado123 + gammado312) + gamma212 * gammado322 +
         2. * (gamma322 * gammado313 + gamma312 * gammado323)) *
            ginv12 +
        (-deldelg1323 + gamma133 * gammado211 + gamma233 * gammado212 +
         (gamma113 + gamma223 + gamma333) * gammado213 + gamma213 * gammado223 +
         (gamma212 + gamma313) * gammado233 +
         gamma123 * (gammado113 + gammado311) + gamma223 * gammado312 +
         gamma112 * (gammado133 + gammado313) + gamma212 * gammado323 +
         2. * (gamma323 * gammado313 + gamma312 * gammado333)) *
            ginv13 +
        (-0.5 * deldelg2223 + gamma123 * gammado212 + gamma223 * gammado222 +
         (gamma222 + gamma323) * gammado223 +
         gamma122 * (gammado123 + gammado312) + gamma222 * gammado322 +
         2. * gamma322 * gammado323) *
            ginv22 +
        (-deldelg2323 + gamma133 * gammado212 + gamma233 * gammado222 +
         (2. * gamma223 + gamma333) * gammado223 +
         (gamma222 + gamma323) * gammado233 +
         gamma123 * (gammado123 + gammado213 + gammado312) +
         gamma122 * (gammado133 + gammado313) + gamma223 * gammado322 +
         (gamma222 + 2. * gamma323) * gammado323 + 2. * gamma322 * gammado333) *
            ginv23 +
        (-0.5 * deldelg3323 + gamma133 * gammado213 + gamma233 * gammado223 +
         (gamma223 + gamma333) * gammado233 +
         gamma123 * (gammado133 + gammado313) + gamma223 * gammado323 +
         2. * gamma323 * gammado333) *
            ginv33;

    R33 =
        delG31 * g13 + delG32 * g23 + delG33 * g33 + gammado313 * Gfromg1 +
        gammado323 * Gfromg2 + gammado333 * Gfromg3 +
        (-0.5 * deldelg1133 + gamma113 * (gammado113 + 2. * gammado311) +
         gamma213 * (gammado213 + 2. * gammado312) + 3. * gamma313 * gammado313) *
            ginv11 +
        (-deldelg1233 + gamma123 * (gammado113 + 2. * gammado311) +
         gamma113 * (gammado123 + 2. * gammado312) +
         gamma223 * (gammado213 + 2. * gammado312) +
         gamma213 * (gammado223 + 2. * gammado322) +
         3. * (gamma323 * gammado313 + gamma313 * gammado323)) *
            ginv12 +
        (-deldelg1333 + gamma133 * (gammado113 + 2. * gammado311) +
         gamma233 * (gammado213 + 2. * gammado312) +
         gamma113 * (gammado133 + 2. * gammado313) +
         gamma213 * (gammado233 + 2. * gammado323) +
         3. * (gamma333 * gammado313 + gamma313 * gammado333)) *
            ginv13 +
        (-0.5 * deldelg2233 + gamma123 * (gammado123 + 2. * gammado312) +
         gamma223 * (gammado223 + 2. * gammado322) + 3. * gamma323 * gammado323) *
            ginv22 +
        (-deldelg2333 + gamma133 * (gammado123 + 2. * gammado312) +
         gamma123 * (gammado133 + 2. * gammado313) +
         gamma233 * (gammado223 + 2. * gammado322) +
         gamma223 * (gammado233 + 2. * gammado323) +
         3. * (gamma333 * gammado323 + gamma323 * gammado333)) *
            ginv23 +
        (-0.5 * deldelg3333 + gamma133 * (gammado133 + 2. * gammado313) +
         gamma233 * (gammado233 + 2. * gammado323) + 3. * gamma333 * gammado333) *
            ginv33;

    chiguard =
        chiDivFloor;

    chiguarded =
        chi;

    if (chiguarded < chiguard)
      chiguarded = chiguard;

    ff =
        chiguarded;

    oochipsipower =
        1 / chipsipower;

    f =
        oochipsipower * log(ff);

    psim4 =
        exp(-4. * f);

    df1 =
        (dchi1 * oochipsipower) / chiguarded;

    df2 =
        (dchi2 * oochipsipower) / chiguarded;

    df3 =
        (dchi3 * oochipsipower) / chiguarded;

    ddf11 =
        (ddchi11 * oochipsipower) / chiguarded - chipsipower * pow2(df1);

    ddf12 =
        -(chipsipower * df1 * df2) + (ddchi12 * oochipsipower) / chiguarded;

    ddf13 =
        -(chipsipower * df1 * df3) + (ddchi13 * oochipsipower) / chiguarded;

    ddf22 =
        (ddchi22 * oochipsipower) / chiguarded - chipsipower * pow2(df2);

    ddf23 =
        -(chipsipower * df2 * df3) + (ddchi23 * oochipsipower) / chiguarded;

    ddf33 =
        (ddchi33 * oochipsipower) / chiguarded - chipsipower * pow2(df3);

    cddf11 =
        ddf11 - df1 * gamma111 - df2 * gamma211 - df3 * gamma311;

    cddf12 =
        ddf12 - df1 * gamma112 - df2 * gamma212 - df3 * gamma312;

    cddf13 =
        ddf13 - df1 * gamma113 - df2 * gamma213 - df3 * gamma313;

    cddf22 =
        ddf22 - df1 * gamma122 - df2 * gamma222 - df3 * gamma322;

    cddf23 =
        ddf23 - df1 * gamma123 - df2 * gamma223 - df3 * gamma323;

    cddf33 =
        ddf33 - df1 * gamma133 - df2 * gamma233 - df3 * gamma333;

    trcddf =
        cddf11 * ginv11 + cddf22 * ginv22 +
        2. * (cddf12 * ginv12 + cddf13 * ginv13 + cddf23 * ginv23) + cddf33 * ginv33;

    Rphi11 =
        -2. * (cddf11 + g11 * trcddf) + (4. - 4. * g11 * ginv11) * pow2(df1) -
        g11 * (8. * (df1 * (df2 * ginv12 + df3 * ginv13) + df2 * df3 * ginv23) +
               4. * (ginv22 * pow2(df2) + ginv33 * pow2(df3)));

    Rphi12 =
        df1 * df2 * (4. - 8. * g12 * ginv12) - 2. * (cddf12 + g12 * trcddf) -
        g12 * (8. * df3 * (df1 * ginv13 + df2 * ginv23) +
               4. * (ginv11 * pow2(df1) + ginv22 * pow2(df2) + ginv33 * pow2(df3)));

    Rphi13 =
        df1 * (4. * df3 - 8. * df2 * g13 * ginv12) - 2. * (cddf13 + g13 * trcddf) -
        g13 * (8. * df3 * (df1 * ginv13 + df2 * ginv23) +
               4. * (ginv11 * pow2(df1) + ginv22 * pow2(df2) + ginv33 * pow2(df3)));

    Rphi22 =
        -2. * (cddf22 + g22 * trcddf) + (4. - 4. * g22 * ginv22) * pow2(df2) -
        g22 * (8. * (df1 * (df2 * ginv12 + df3 * ginv13) + df2 * df3 * ginv23) +
               4. * (ginv11 * pow2(df1) + ginv33 * pow2(df3)));

    Rphi23 =
        df2 * (-8. * df1 * g23 * ginv12 + df3 * (4. - 8. * g23 * ginv23)) -
        2. * (cddf23 + g23 * trcddf) - g23 * (8. * df1 * df3 * ginv13 + 4. * (ginv11 * pow2(df1) + ginv22 * pow2(df2) + ginv33 * pow2(df3)));

    Rphi33 =
        -2. * (cddf33 + g33 * trcddf) - g33 * (8. * (df1 * (df2 * ginv12 + df3 * ginv13) + df2 * df3 * ginv23) + 4. * (ginv11 * pow2(df1) + ginv22 * pow2(df2))) +
        (4. - 4. * g33 * ginv33) * pow2(df3);

    cdda11 =
        dda11 - da2 * gamma211 - da3 * gamma311 +
        da1 * (-gamma111 + df1 * (-4. + 2. * g11 * ginv11)) +
        2. * g11 * ((da2 * df1 + da1 * df2) * ginv12 + (da3 * df1 + da1 * df3) * ginv13 + da2 * df2 * ginv22 + (da3 * df2 + da2 * df3) * ginv23 + da3 * df3 * ginv33);

    cdda12 =
        dda12 - da1 * gamma112 - da2 * gamma212 - da3 * gamma312 +
        2. * (-(da2 * df1) - da1 * df2 + g12 * (da1 * df1 * ginv11 + (da2 * df1 + da1 * df2) * ginv12 + (da3 * df1 + da1 * df3) * ginv13 + da2 * df2 * ginv22 + (da3 * df2 + da2 * df3) * ginv23 + da3 * df3 * ginv33));

    cdda13 =
        dda13 - da1 * gamma113 - da2 * gamma213 - da3 * gamma313 +
        2. * (-(da3 * df1) - da1 * df3 + g13 * (da1 * df1 * ginv11 + (da2 * df1 + da1 * df2) * ginv12 + (da3 * df1 + da1 * df3) * ginv13 + da2 * df2 * ginv22 + (da3 * df2 + da2 * df3) * ginv23 + da3 * df3 * ginv33));

    cdda22 =
        dda22 - da1 * gamma122 - da2 * (4. * df2 + gamma222) - da3 * gamma322 +
        2. * g22 * (da1 * df1 * ginv11 + (da2 * df1 + da1 * df2) * ginv12 + (da3 * df1 + da1 * df3) * ginv13 + da2 * df2 * ginv22 + (da3 * df2 + da2 * df3) * ginv23 + da3 * df3 * ginv33);

    cdda23 =
        dda23 - da1 * gamma123 - da2 * gamma223 - da3 * gamma323 +
        2. * (-(da3 * df2) - da2 * df3 + g23 * (da1 * df1 * ginv11 + (da2 * df1 + da1 * df2) * ginv12 + (da3 * df1 + da1 * df3) * ginv13 + da2 * df2 * ginv22 + (da3 * df2 + da2 * df3) * ginv23 + da3 * df3 * ginv33));

    cdda33 =
        dda33 - da1 * gamma133 - da2 * gamma233 - da3 * (4. * df3 + gamma333) +
        2. * g33 * (da1 * df1 * ginv11 + (da2 * df1 + da1 * df2) * ginv12 + (da3 * df1 + da1 * df3) * ginv13 + da2 * df2 * ginv22 + (da3 * df2 + da2 * df3) * ginv23 + da3 * df3 * ginv33);

    trcdda =
        (cdda11 * ginv11 + cdda22 * ginv22 +
         2. * (cdda12 * ginv12 + cdda13 * ginv13 + cdda23 * ginv23) + cdda33 * ginv33) *
        psim4;

    AA11 =
        2. * (A11 * (A12 * ginv12 + A13 * ginv13) + A12 * A13 * ginv23) + ginv11 * pow2(A11) +
        ginv22 * pow2(A12) + ginv33 * pow2(A13);

    AA12 =
        (A12 * A13 + A11 * A23) * ginv13 + A12 * (A11 * ginv11 + A22 * ginv22) +
        (A13 * A22 + A12 * A23) * ginv23 + A13 * A23 * ginv33 + ginv12 * (A11 * A22 + pow2(A12));

    AA13 =
        (A12 * A13 + A11 * A23) * ginv12 + A12 * A23 * ginv22 + (A13 * A23 + A12 * A33) * ginv23 +
        A13 * (A11 * ginv11 + A33 * ginv33) + ginv13 * (A11 * A33 + pow2(A13));

    AA22 =
        2. * (A12 * (A22 * ginv12 + A23 * ginv13) + A22 * A23 * ginv23) + ginv11 * pow2(A12) +
        ginv22 * pow2(A22) + ginv33 * pow2(A23);

    AA23 =
        A12 * A13 * ginv11 + (A13 * A22 + A12 * A23) * ginv12 + (A13 * A23 + A12 * A33) * ginv13 +
        A23 * (A22 * ginv22 + A33 * ginv33) + ginv23 * (A22 * A33 + pow2(A23));

    AA33 =
        2. * (A13 * (A23 * ginv12 + A33 * ginv13) + A23 * A33 * ginv23) + ginv11 * pow2(A13) +
        ginv22 * pow2(A23) + ginv33 * pow2(A33);

    cAA =
        AA11 * ginv11 + AA22 * ginv22 + 2. * (AA12 * ginv12 + AA13 * ginv13 + AA23 * ginv23) +
        AA33 * ginv33;

    Ainv11 =
        2. * (A23 * ginv12 * ginv13 + ginv11 * (A12 * ginv12 + A13 * ginv13)) +
        A11 * pow2(ginv11) + A22 * pow2(ginv12) + A33 * pow2(ginv13);

    Ainv12 =
        ginv11 * (A11 * ginv12 + A12 * ginv22 + A13 * ginv23) +
        ginv12 * (A13 * ginv13 + A22 * ginv22 + A23 * ginv23) +
        ginv13 * (A23 * ginv22 + A33 * ginv23) + A12 * pow2(ginv12);

    Ainv13 =
        ginv11 * (A11 * ginv13 + A12 * ginv23 + A13 * ginv33) +
        ginv12 * (A12 * ginv13 + A22 * ginv23 + A23 * ginv33) +
        ginv13 * (A23 * ginv23 + A33 * ginv33) + A13 * pow2(ginv13);

    Ainv22 =
        2. * (A23 * ginv22 * ginv23 + ginv12 * (A12 * ginv22 + A13 * ginv23)) +
        A11 * pow2(ginv12) + A22 * pow2(ginv22) + A33 * pow2(ginv23);

    Ainv23 =
        ginv13 * (A12 * ginv22 + A13 * ginv23) + A33 * ginv23 * ginv33 +
        ginv12 * (A11 * ginv13 + A12 * ginv23 + A13 * ginv33) +
        ginv22 * (A22 * ginv23 + A23 * ginv33) + A23 * pow2(ginv23);

    Ainv33 =
        2. * (A23 * ginv23 * ginv33 + ginv13 * (A12 * ginv23 + A13 * ginv33)) +
        A11 * pow2(ginv13) + A22 * pow2(ginv23) + A33 * pow2(ginv33);

    divAinv1 =
        (-1.5 * (Ainv11 * dchi1 + Ainv12 * dchi2 + Ainv13 * dchi3)) / chiguarded +
        Ainv11 * gamma111 + Ainv22 * gamma122 +
        2. * (Ainv12 * gamma112 + Ainv13 * gamma113 + Ainv23 * gamma123) +
        Ainv33 * gamma133 - (0.66666666666666666667 * dKhat1 + 0.33333333333333333333 * dTheta1) * ginv11 -
        (0.66666666666666666667 * dKhat2 + 0.33333333333333333333 * dTheta2) * ginv12 -
        (0.66666666666666666667 * dKhat3 + 0.33333333333333333333 * dTheta3) * ginv13;

    divAinv2 =
        (-1.5 * (Ainv12 * dchi1 + Ainv22 * dchi2 + Ainv23 * dchi3)) / chiguarded +
        Ainv11 * gamma211 + Ainv22 * gamma222 +
        2. * (Ainv12 * gamma212 + Ainv13 * gamma213 + Ainv23 * gamma223) +
        Ainv33 * gamma233 - (0.66666666666666666667 * dKhat1 + 0.33333333333333333333 * dTheta1) * ginv12 -
        (0.66666666666666666667 * dKhat2 + 0.33333333333333333333 * dTheta2) * ginv22 -
        (0.66666666666666666667 * dKhat3 + 0.33333333333333333333 * dTheta3) * ginv23;

    divAinv3 =
        (-1.5 * (Ainv13 * dchi1 + Ainv23 * dchi2 + Ainv33 * dchi3)) / chiguarded +
        Ainv11 * gamma311 + Ainv22 * gamma322 +
        2. * (Ainv12 * gamma312 + Ainv13 * gamma313 + Ainv23 * gamma323) +
        Ainv33 * gamma333 - (0.66666666666666666667 * dKhat1 + 0.33333333333333333333 * dTheta1) * ginv13 -
        (0.66666666666666666667 * dKhat2 + 0.33333333333333333333 * dTheta2) * ginv23 -
        (0.66666666666666666667 * dKhat3 + 0.33333333333333333333 * dTheta3) * ginv33;

    Rhat =
        psim4 * (ginv11 * (R11 + Rphi11) + ginv22 * (R22 + Rphi22) +
                 2. * (ginv12 * (R12 + Rphi12) + ginv13 * (R13 + Rphi13) +
                       ginv23 * (R23 + Rphi23)) +
                 ginv33 * (R33 + Rphi33));

    Hhat =
        -cAA + Rhat + 0.66666666666666666667 * pow2(K);

    divbeta =
        db11 + db22 + db33;

    totdivbeta =
        0.66666666666666666667 * divbeta;

    ootddivbeta1 =
        0.33333333333333333333 * (ddb111 + ddb122 + ddb133);

    ootddivbeta2 =
        0.33333333333333333333 * (ddb121 + ddb222 + ddb233);

    ootddivbeta3 =
        0.33333333333333333333 * (ddb131 + ddb232 + ddb333);

    lieg11 =
        2. * (db11 * g11 + db12 * g12 + db13 * g13) - g11 * totdivbeta;

    lieg12 =
        db21 * g11 + db23 * g13 + db12 * g22 + db13 * g23 + g12 * (db11 + db22 - totdivbeta);

    lieg13 =
        db31 * g11 + db32 * g12 + db12 * g23 + db13 * g33 + g13 * (db11 + db33 - totdivbeta);

    lieg22 =
        2. * (db21 * g12 + db22 * g22 + db23 * g23) - g22 * totdivbeta;

    lieg23 =
        db31 * g12 + db21 * g13 + db32 * g22 + db23 * g33 + g23 * (db22 + db33 - totdivbeta);

    lieg33 =
        2. * (db31 * g13 + db32 * g23 + db33 * g33) - g33 * totdivbeta;

    lieA11 =
        2. * (A11 * db11 + A12 * db12 + A13 * db13) - A11 * totdivbeta;

    lieA12 =
        A22 * db12 + A23 * db13 + A11 * db21 + A13 * db23 + A12 * (db11 + db22 - totdivbeta);

    lieA13 =
        A23 * db12 + A33 * db13 + A11 * db31 + A12 * db32 + A13 * (db11 + db33 - totdivbeta);

    lieA22 =
        2. * (A12 * db21 + A22 * db22 + A23 * db23) - A22 * totdivbeta;

    lieA23 =
        A13 * db21 + A33 * db23 + A12 * db31 + A22 * db32 + A23 * (db22 + db33 - totdivbeta);

    lieA33 =
        2. * (A13 * db31 + A23 * db32 + A33 * db33) - A33 * totdivbeta;

    liechi =
        0.16666666666666666667 * chiguarded * chipsipower * divbeta;

    pseudolieG1 =
        -(db11 * Gfromg1) - db21 * Gfromg2 - db31 * Gfromg3 + ddb221 * ginv22 +
        2. * ddb231 * ginv23 + ddb331 * ginv33 + ginv11 * (ddb111 + ootddivbeta1) +
        ginv12 * (2. * ddb121 + ootddivbeta2) + ginv13 * (2. * ddb131 + ootddivbeta3) +
        Gfromg1 * totdivbeta;

    pseudolieG2 =
        -(db12 * Gfromg1) - db22 * Gfromg2 - db32 * Gfromg3 + ddb112 * ginv11 +
        2. * ddb132 * ginv13 + ddb332 * ginv33 + ginv12 * (2. * ddb122 + ootddivbeta1) +
        ginv22 * (ddb222 + ootddivbeta2) + ginv23 * (2. * ddb232 + ootddivbeta3) +
        Gfromg2 * totdivbeta;

    pseudolieG3 =
        -(db13 * Gfromg1) - db23 * Gfromg2 - db33 * Gfromg3 + ddb113 * ginv11 +
        2. * ddb123 * ginv12 + ddb223 * ginv22 + ginv13 * (2. * ddb133 + ootddivbeta1) +
        ginv23 * (2. * ddb233 + ootddivbeta2) + ginv33 * (ddb333 + ootddivbeta3) +
        Gfromg3 * totdivbeta;

    rg11 =
        -2. * A11 * alpha + lieg11;

    rg12 =
        -2. * A12 * alpha + lieg12;

    rg13 =
        -2. * A13 * alpha + lieg13;

    rg22 =
        -2. * A22 * alpha + lieg22;

    rg23 =
        -2. * A23 * alpha + lieg23;

    rg33 =
        -2. * A33 * alpha + lieg33;

    rA11 =
        lieA11 + alpha * (-2. * AA11 + A11 * K + psim4 * R11 - 0.33333333333333333333 * g11 * Rhat) + psim4 * (-cdda11 + alpha * Rphi11) +
        0.33333333333333333333 * g11 * trcdda;

    rA12 =
        lieA12 + alpha * (-2. * AA12 + A12 * K + psim4 * R12 - 0.33333333333333333333 * g12 * Rhat) + psim4 * (-cdda12 + alpha * Rphi12) +
        0.33333333333333333333 * g12 * trcdda;

    rA13 =
        lieA13 + alpha * (-2. * AA13 + A13 * K + psim4 * R13 - 0.33333333333333333333 * g13 * Rhat) + psim4 * (-cdda13 + alpha * Rphi13) +
        0.33333333333333333333 * g13 * trcdda;

    rA22 =
        lieA22 + alpha * (-2. * AA22 + A22 * K + psim4 * R22 - 0.33333333333333333333 * g22 * Rhat) + psim4 * (-cdda22 + alpha * Rphi22) +
        0.33333333333333333333 * g22 * trcdda;

    rA23 =
        lieA23 + alpha * (-2. * AA23 + A23 * K + psim4 * R23 - 0.33333333333333333333 * g23 * Rhat) + psim4 * (-cdda23 + alpha * Rphi23) +
        0.33333333333333333333 * g23 * trcdda;

    rA33 =
        lieA33 + alpha * (-2. * AA33 + A33 * K + psim4 * R33 - 0.33333333333333333333 * g33 * Rhat) + psim4 * (-cdda33 + alpha * Rphi33) +
        0.33333333333333333333 * g33 * trcdda;

    rG1 =
        -2. * (Ainv11 * da1 + Ainv12 * da2 + Ainv13 * da3) +
        alpha * (2. * divAinv1 + 2. * (-G1 + Gfromg1) * kappa1) + pseudolieG1;

    rG2 =
        -2. * (Ainv12 * da1 + Ainv22 * da2 + Ainv23 * da3) +
        alpha * (2. * divAinv2 + 2. * (-G2 + Gfromg2) * kappa1) + pseudolieG2;

    rG3 =
        -2. * (Ainv13 * da1 + Ainv23 * da2 + Ainv33 * da3) +
        alpha * (2. * divAinv3 + 2. * (-G3 + Gfromg3) * kappa1) + pseudolieG3;

    rKhat =
        -trcdda + alpha * (cAA + kappa1 * (Theta - kappa2 * Theta) +
                           0.33333333333333333333 * pow2(K));

    rchi =
        -0.16666666666666666667 * alpha * chiguarded * chipsipower * K + liechi;

    rTheta =
        alpha * (0.5 * Hhat - kappa1 * (2. + kappa2) * Theta);

#if 0   
// this part is for CCZ4
dginv111
=
-2.*(delg123*ginv12*ginv13 + ginv11*(delg112*ginv12 + delg113*ginv13)) - 
  delg111*pow2(ginv11) - delg122*pow2(ginv12) - delg133*pow2(ginv13)
;

dginv112
=
-(ginv11*(delg111*ginv12 + delg112*ginv22 + delg113*ginv23)) - 
  ginv12*(delg113*ginv13 + delg122*ginv22 + delg123*ginv23) - 
  ginv13*(delg123*ginv22 + delg133*ginv23) - delg112*pow2(ginv12)
;

dginv113
=
-(ginv11*(delg111*ginv13 + delg112*ginv23 + delg113*ginv33)) - 
  ginv12*(delg112*ginv13 + delg122*ginv23 + delg123*ginv33) - 
  ginv13*(delg123*ginv23 + delg133*ginv33) - delg113*pow2(ginv13)
;

dginv122
=
-2.*(delg123*ginv22*ginv23 + ginv12*(delg112*ginv22 + delg113*ginv23)) - 
  delg111*pow2(ginv12) - delg122*pow2(ginv22) - delg133*pow2(ginv23)
;

dginv123
=
-(ginv13*(delg112*ginv22 + delg113*ginv23)) - delg133*ginv23*ginv33 - 
  ginv12*(delg111*ginv13 + delg112*ginv23 + delg113*ginv33) - 
  ginv22*(delg122*ginv23 + delg123*ginv33) - delg123*pow2(ginv23)
;

dginv133
=
-2.*(delg123*ginv23*ginv33 + ginv13*(delg112*ginv23 + delg113*ginv33)) - 
  delg111*pow2(ginv13) - delg122*pow2(ginv23) - delg133*pow2(ginv33)
;

dginv211
=
-2.*(delg223*ginv12*ginv13 + ginv11*(delg212*ginv12 + delg213*ginv13)) - 
  delg211*pow2(ginv11) - delg222*pow2(ginv12) - delg233*pow2(ginv13)
;

dginv212
=
-(ginv11*(delg211*ginv12 + delg212*ginv22 + delg213*ginv23)) - 
  ginv12*(delg213*ginv13 + delg222*ginv22 + delg223*ginv23) - 
  ginv13*(delg223*ginv22 + delg233*ginv23) - delg212*pow2(ginv12)
;

dginv213
=
-(ginv11*(delg211*ginv13 + delg212*ginv23 + delg213*ginv33)) - 
  ginv12*(delg212*ginv13 + delg222*ginv23 + delg223*ginv33) - 
  ginv13*(delg223*ginv23 + delg233*ginv33) - delg213*pow2(ginv13)
;

dginv222
=
-2.*(delg223*ginv22*ginv23 + ginv12*(delg212*ginv22 + delg213*ginv23)) - 
  delg211*pow2(ginv12) - delg222*pow2(ginv22) - delg233*pow2(ginv23)
;

dginv223
=
-(ginv13*(delg212*ginv22 + delg213*ginv23)) - delg233*ginv23*ginv33 - 
  ginv12*(delg211*ginv13 + delg212*ginv23 + delg213*ginv33) - 
  ginv22*(delg222*ginv23 + delg223*ginv33) - delg223*pow2(ginv23)
;

dginv233
=
-2.*(delg223*ginv23*ginv33 + ginv13*(delg212*ginv23 + delg213*ginv33)) - 
  delg211*pow2(ginv13) - delg222*pow2(ginv23) - delg233*pow2(ginv33)
;

dginv311
=
-2.*(delg323*ginv12*ginv13 + ginv11*(delg312*ginv12 + delg313*ginv13)) - 
  delg311*pow2(ginv11) - delg322*pow2(ginv12) - delg333*pow2(ginv13)
;

dginv312
=
-(ginv11*(delg311*ginv12 + delg312*ginv22 + delg313*ginv23)) - 
  ginv12*(delg313*ginv13 + delg322*ginv22 + delg323*ginv23) - 
  ginv13*(delg323*ginv22 + delg333*ginv23) - delg312*pow2(ginv12)
;

dginv313
=
-(ginv11*(delg311*ginv13 + delg312*ginv23 + delg313*ginv33)) - 
  ginv12*(delg312*ginv13 + delg322*ginv23 + delg323*ginv33) - 
  ginv13*(delg323*ginv23 + delg333*ginv33) - delg313*pow2(ginv13)
;

dginv322
=
-2.*(delg323*ginv22*ginv23 + ginv12*(delg312*ginv22 + delg313*ginv23)) - 
  delg311*pow2(ginv12) - delg322*pow2(ginv22) - delg333*pow2(ginv23)
;

dginv323
=
-(ginv13*(delg312*ginv22 + delg313*ginv23)) - delg333*ginv23*ginv33 - 
  ginv12*(delg311*ginv13 + delg312*ginv23 + delg313*ginv33) - 
  ginv22*(delg322*ginv23 + delg323*ginv33) - delg323*pow2(ginv23)
;

dginv333
=
-2.*(delg323*ginv23*ginv33 + ginv13*(delg312*ginv23 + delg313*ginv33)) - 
  delg311*pow2(ginv13) - delg322*pow2(ginv23) - delg333*pow2(ginv33)
;

dphi1
=
(-0.25*dchi1)/chiguarded
;

dphi2
=
(-0.25*dchi2)/chiguarded
;

dphi3
=
(-0.25*dchi3)/chiguarded
;

gammaF111
=
gamma111 + dphi1*(4. - 2.*g11*ginv11) - 2.*g11*(dphi2*ginv12 + dphi3*ginv13)
;

gammaF112
=
gamma112 + dphi2*(2. - 2.*g12*ginv12) - 2.*g12*(dphi1*ginv11 + dphi3*ginv13)
;

gammaF113
=
gamma113 - 2.*g13*(dphi1*ginv11 + dphi2*ginv12) + dphi3*(2. - 2.*g13*ginv13)
;

gammaF121
=
gamma112 + dphi2*(2. - 2.*g12*ginv12) - 2.*g12*(dphi1*ginv11 + dphi3*ginv13)
;

gammaF122
=
gamma122 - 2.*g22*(dphi1*ginv11 + dphi2*ginv12 + dphi3*ginv13)
;

gammaF123
=
gamma123 - 2.*g23*(dphi1*ginv11 + dphi2*ginv12 + dphi3*ginv13)
;

gammaF131
=
gamma113 - 2.*g13*(dphi1*ginv11 + dphi2*ginv12) + dphi3*(2. - 2.*g13*ginv13)
;

gammaF132
=
gamma123 - 2.*g23*(dphi1*ginv11 + dphi2*ginv12 + dphi3*ginv13)
;

gammaF133
=
gamma133 - 2.*g33*(dphi1*ginv11 + dphi2*ginv12 + dphi3*ginv13)
;

gammaF211
=
gamma211 - 2.*g11*(dphi1*ginv12 + dphi2*ginv22 + dphi3*ginv23)
;

gammaF212
=
gamma212 + dphi1*(2. - 2.*g12*ginv12) - 2.*g12*(dphi2*ginv22 + dphi3*ginv23)
;

gammaF213
=
gamma213 - 2.*g13*(dphi1*ginv12 + dphi2*ginv22 + dphi3*ginv23)
;

gammaF221
=
gamma212 + dphi1*(2. - 2.*g12*ginv12) - 2.*g12*(dphi2*ginv22 + dphi3*ginv23)
;

gammaF222
=
gamma222 + dphi2*(4. - 2.*g22*ginv22) - 2.*g22*(dphi1*ginv12 + dphi3*ginv23)
;

gammaF223
=
gamma223 - 2.*g23*(dphi1*ginv12 + dphi2*ginv22) + dphi3*(2. - 2.*g23*ginv23)
;

gammaF231
=
gamma213 - 2.*g13*(dphi1*ginv12 + dphi2*ginv22 + dphi3*ginv23)
;

gammaF232
=
gamma223 - 2.*g23*(dphi1*ginv12 + dphi2*ginv22) + dphi3*(2. - 2.*g23*ginv23)
;

gammaF233
=
gamma233 - 2.*g33*(dphi1*ginv12 + dphi2*ginv22 + dphi3*ginv23)
;

gammaF311
=
gamma311 - 2.*g11*(dphi1*ginv13 + dphi2*ginv23 + dphi3*ginv33)
;

gammaF312
=
gamma312 - 2.*g12*(dphi1*ginv13 + dphi2*ginv23 + dphi3*ginv33)
;

gammaF313
=
gamma313 + dphi1*(2. - 2.*g13*ginv13) - 2.*g13*(dphi2*ginv23 + dphi3*ginv33)
;

gammaF321
=
gamma312 - 2.*g12*(dphi1*ginv13 + dphi2*ginv23 + dphi3*ginv33)
;

gammaF322
=
gamma322 - 2.*g22*(dphi1*ginv13 + dphi2*ginv23 + dphi3*ginv33)
;

gammaF323
=
gamma323 + dphi2*(2. - 2.*g23*ginv23) - 2.*g23*(dphi1*ginv13 + dphi3*ginv33)
;

gammaF331
=
gamma313 + dphi1*(2. - 2.*g13*ginv13) - 2.*g13*(dphi2*ginv23 + dphi3*ginv33)
;

gammaF332
=
gamma323 + dphi2*(2. - 2.*g23*ginv23) - 2.*g23*(dphi1*ginv13 + dphi3*ginv33)
;

gammaF333
=
gamma333 - 2.*g33*(dphi1*ginv13 + dphi2*ginv23) + dphi3*(4. - 2.*g33*ginv33)
;

Gd1
=
ginv11*((2.*delg112 + delg211)*ginv12 + (2.*delg113 + delg311)*ginv13 + 
     delg212*ginv22 + (delg213 + delg312)*ginv23 + delg313*ginv33) + 
  ginv12*((2.*delg123 + delg213 + delg312)*ginv13 + delg222*ginv22 + 
     (delg223 + delg322)*ginv23 + delg323*ginv33) + 
  ginv13*(delg223*ginv22 + (delg233 + delg323)*ginv23 + delg333*ginv33) + 
  delg111*pow2(ginv11) + (delg122 + delg212)*pow2(ginv12) + 
  (delg133 + delg313)*pow2(ginv13)
;

Gd2
=
ginv11*(delg111*ginv12 + delg112*ginv22 + delg113*ginv23) + 
  ginv13*((delg123 + delg312)*ginv22 + (delg133 + delg313)*ginv23) + 
  delg333*ginv23*ginv33 + ginv12*
   ((delg113 + delg311)*ginv13 + (delg122 + 2.*delg212)*ginv22 + 
     (delg123 + 2.*delg213 + delg312)*ginv23 + delg313*ginv33) + 
  ginv22*((2.*delg223 + delg322)*ginv23 + delg323*ginv33) + 
  (delg112 + delg211)*pow2(ginv12) + delg222*pow2(ginv22) + 
  (delg233 + delg323)*pow2(ginv23)
;

Gd3
=
(delg233 + 2.*delg323)*ginv23*ginv33 + 
  ginv11*(delg111*ginv13 + delg112*ginv23 + delg113*ginv33) + 
  ginv12*((delg112 + delg211)*ginv13 + (delg122 + delg212)*ginv23 + 
     (delg123 + delg213)*ginv33) + 
  ginv22*(delg222*ginv23 + delg223*ginv33) + 
  ginv13*(delg212*ginv22 + (delg123 + delg213 + 2.*delg312)*ginv23 + 
     (delg133 + 2.*delg313)*ginv33) + (delg113 + delg311)*pow2(ginv13) + 
  (delg223 + delg322)*pow2(ginv23) + delg333*pow2(ginv33)
;

dGd11
=
(delg212*dginv111 + delg222*dginv112 + delg223*dginv113)*ginv22 + 
  ((delg213 + delg312)*dginv111 + (delg223 + delg322)*dginv112 + 
     (delg233 + delg323)*dginv113)*ginv23 + 
  (delg313*dginv111 + delg323*dginv112 + delg333*dginv113)*ginv33 + 
  ginv11*(delg211*dginv112 + delg311*dginv113 + 
     2.*(delg111*dginv111 + delg112*dginv112 + delg113*dginv113) + 
     delg212*dginv122 + (delg213 + delg312)*dginv123 + delg313*dginv133 + 
     (2.*deldelg1112 + deldelg1211)*ginv12 + 
     (2.*deldelg1113 + deldelg1311)*ginv13 + deldelg1212*ginv22 + 
     (deldelg1213 + deldelg1312)*ginv23 + deldelg1313*ginv33) + 
  ginv12*((2.*delg112 + delg211)*dginv111 + (delg213 + delg312)*dginv113 + 
     2.*((delg122 + delg212)*dginv112 + delg123*dginv113) + 
     delg222*dginv122 + (delg223 + delg322)*dginv123 + delg323*dginv133 + 
     (2.*deldelg1123 + deldelg1213 + deldelg1312)*ginv13 + 
     deldelg1222*ginv22 + (deldelg1223 + deldelg1322)*ginv23 + 
     deldelg1323*ginv33) + ginv13*
   ((2.*delg113 + delg311)*dginv111 + 
     (2.*delg123 + delg213 + delg312)*dginv112 + 
     2.*(delg133 + delg313)*dginv113 + delg223*dginv122 + 
     (delg233 + delg323)*dginv123 + delg333*dginv133 + deldelg1223*ginv22 + 
     (deldelg1233 + deldelg1323)*ginv23 + deldelg1333*ginv33) + 
  deldelg1111*pow2(ginv11) + (deldelg1122 + deldelg1212)*pow2(ginv12) + 
  (deldelg1133 + deldelg1313)*pow2(ginv13)
;

dGd12
=
ginv11*(delg111*dginv112 + delg112*dginv122 + delg113*dginv123 + 
     deldelg1111*ginv12 + deldelg1112*ginv22 + deldelg1113*ginv23) + 
  ginv13*((delg113 + delg311)*dginv112 + (delg123 + delg312)*dginv122 + 
     (delg133 + delg313)*dginv123 + (deldelg1123 + deldelg1312)*ginv22 + 
     (deldelg1133 + deldelg1313)*ginv23) + 
  (delg313*dginv112 + delg323*dginv122 + delg333*dginv123)*ginv33 + 
  ginv12*(delg111*dginv111 + (delg113 + delg311)*dginv113 + 
     delg122*dginv122 + (delg123 + delg312)*dginv123 + 
     2.*((delg112 + delg211)*dginv112 + delg212*dginv122 + 
        delg213*dginv123) + delg313*dginv133 + 
     (deldelg1113 + deldelg1311)*ginv13 + 
     (deldelg1122 + 2.*deldelg1212)*ginv22 + 
     (deldelg1123 + 2.*deldelg1213 + deldelg1312)*ginv23 + 
     deldelg1313*ginv33) + ginv22*
   (delg112*dginv111 + (delg122 + 2.*delg212)*dginv112 + 
     (delg123 + delg312)*dginv113 + delg322*dginv123 + 
     2.*(delg222*dginv122 + delg223*dginv123) + delg323*dginv133 + 
     (2.*deldelg1223 + deldelg1322)*ginv23 + deldelg1323*ginv33) + 
  ginv23*(delg113*dginv111 + (delg123 + 2.*delg213 + delg312)*dginv112 + 
     (delg133 + delg313)*dginv113 + (2.*delg223 + delg322)*dginv122 + 
     2.*(delg233 + delg323)*dginv123 + delg333*dginv133 + deldelg1333*ginv33\
) + (deldelg1112 + deldelg1211)*pow2(ginv12) + deldelg1222*pow2(ginv22) + 
  (deldelg1233 + deldelg1323)*pow2(ginv23)
;

dGd13
=
(delg113*dginv111 + (delg123 + delg213)*dginv112 + 
     (delg133 + 2.*delg313)*dginv113 + delg223*dginv122 + 
     (delg233 + 2.*delg323)*dginv123 + 2.*delg333*dginv133)*ginv33 + 
  ginv11*(delg111*dginv113 + delg112*dginv123 + delg113*dginv133 + 
     deldelg1111*ginv13 + deldelg1112*ginv23 + deldelg1113*ginv33) + 
  ginv12*((delg112 + delg211)*dginv113 + (delg122 + delg212)*dginv123 + 
     (delg123 + delg213)*dginv133 + (deldelg1112 + deldelg1211)*ginv13 + 
     (deldelg1122 + deldelg1212)*ginv23 + (deldelg1123 + deldelg1213)*ginv33\
) + ginv22*(delg212*dginv113 + delg222*dginv123 + delg223*dginv133 + 
     deldelg1222*ginv23 + deldelg1223*ginv33) + 
  ginv13*(delg111*dginv111 + (delg112 + delg211)*dginv112 + 
     delg212*dginv122 + (delg123 + delg213)*dginv123 + delg133*dginv133 + 
     2.*((delg113 + delg311)*dginv113 + delg312*dginv123 + 
        delg313*dginv133) + deldelg1212*ginv22 + 
     (deldelg1123 + deldelg1213 + 2.*deldelg1312)*ginv23 + 
     (deldelg1133 + 2.*deldelg1313)*ginv33) + 
  ginv23*(delg112*dginv111 + (delg122 + delg212)*dginv112 + 
     (delg123 + delg213 + 2.*delg312)*dginv113 + delg222*dginv122 + 
     delg233*dginv133 + 2.*((delg223 + delg322)*dginv123 + 
        delg323*dginv133) + (deldelg1233 + 2.*deldelg1323)*ginv33) + 
  (deldelg1113 + deldelg1311)*pow2(ginv13) + 
  (deldelg1223 + deldelg1322)*pow2(ginv23) + deldelg1333*pow2(ginv33)
;

dGd21
=
(delg212*dginv211 + delg222*dginv212 + delg223*dginv213)*ginv22 + 
  ((delg213 + delg312)*dginv211 + (delg223 + delg322)*dginv212 + 
     (delg233 + delg323)*dginv213)*ginv23 + 
  (delg313*dginv211 + delg323*dginv212 + delg333*dginv213)*ginv33 + 
  ginv11*(delg211*dginv212 + delg311*dginv213 + 
     2.*(delg111*dginv211 + delg112*dginv212 + delg113*dginv213) + 
     delg212*dginv222 + (delg213 + delg312)*dginv223 + delg313*dginv233 + 
     (2.*deldelg1212 + deldelg2211)*ginv12 + 
     (2.*deldelg1213 + deldelg2311)*ginv13 + deldelg2212*ginv22 + 
     (deldelg2213 + deldelg2312)*ginv23 + deldelg2313*ginv33) + 
  ginv12*((2.*delg112 + delg211)*dginv211 + (delg213 + delg312)*dginv213 + 
     2.*((delg122 + delg212)*dginv212 + delg123*dginv213) + 
     delg222*dginv222 + (delg223 + delg322)*dginv223 + delg323*dginv233 + 
     (2.*deldelg1223 + deldelg2213 + deldelg2312)*ginv13 + 
     deldelg2222*ginv22 + (deldelg2223 + deldelg2322)*ginv23 + 
     deldelg2323*ginv33) + ginv13*
   ((2.*delg113 + delg311)*dginv211 + 
     (2.*delg123 + delg213 + delg312)*dginv212 + 
     2.*(delg133 + delg313)*dginv213 + delg223*dginv222 + 
     (delg233 + delg323)*dginv223 + delg333*dginv233 + deldelg2223*ginv22 + 
     (deldelg2233 + deldelg2323)*ginv23 + deldelg2333*ginv33) + 
  deldelg1211*pow2(ginv11) + (deldelg1222 + deldelg2212)*pow2(ginv12) + 
  (deldelg1233 + deldelg2313)*pow2(ginv13)
;

dGd22
=
ginv11*(delg111*dginv212 + delg112*dginv222 + delg113*dginv223 + 
     deldelg1211*ginv12 + deldelg1212*ginv22 + deldelg1213*ginv23) + 
  ginv13*((delg113 + delg311)*dginv212 + (delg123 + delg312)*dginv222 + 
     (delg133 + delg313)*dginv223 + (deldelg1223 + deldelg2312)*ginv22 + 
     (deldelg1233 + deldelg2313)*ginv23) + 
  (delg313*dginv212 + delg323*dginv222 + delg333*dginv223)*ginv33 + 
  ginv12*(delg111*dginv211 + (delg113 + delg311)*dginv213 + 
     delg122*dginv222 + (delg123 + delg312)*dginv223 + 
     2.*((delg112 + delg211)*dginv212 + delg212*dginv222 + 
        delg213*dginv223) + delg313*dginv233 + 
     (deldelg1213 + deldelg2311)*ginv13 + 
     (deldelg1222 + 2.*deldelg2212)*ginv22 + 
     (deldelg1223 + 2.*deldelg2213 + deldelg2312)*ginv23 + 
     deldelg2313*ginv33) + ginv22*
   (delg112*dginv211 + (delg122 + 2.*delg212)*dginv212 + 
     (delg123 + delg312)*dginv213 + delg322*dginv223 + 
     2.*(delg222*dginv222 + delg223*dginv223) + delg323*dginv233 + 
     (2.*deldelg2223 + deldelg2322)*ginv23 + deldelg2323*ginv33) + 
  ginv23*(delg113*dginv211 + (delg123 + 2.*delg213 + delg312)*dginv212 + 
     (delg133 + delg313)*dginv213 + (2.*delg223 + delg322)*dginv222 + 
     2.*(delg233 + delg323)*dginv223 + delg333*dginv233 + deldelg2333*ginv33\
) + (deldelg1212 + deldelg2211)*pow2(ginv12) + deldelg2222*pow2(ginv22) + 
  (deldelg2233 + deldelg2323)*pow2(ginv23)
;

dGd23
=
(delg113*dginv211 + (delg123 + delg213)*dginv212 + 
     (delg133 + 2.*delg313)*dginv213 + delg223*dginv222 + 
     (delg233 + 2.*delg323)*dginv223 + 2.*delg333*dginv233)*ginv33 + 
  ginv11*(delg111*dginv213 + delg112*dginv223 + delg113*dginv233 + 
     deldelg1211*ginv13 + deldelg1212*ginv23 + deldelg1213*ginv33) + 
  ginv12*((delg112 + delg211)*dginv213 + (delg122 + delg212)*dginv223 + 
     (delg123 + delg213)*dginv233 + (deldelg1212 + deldelg2211)*ginv13 + 
     (deldelg1222 + deldelg2212)*ginv23 + (deldelg1223 + deldelg2213)*ginv33\
) + ginv22*(delg212*dginv213 + delg222*dginv223 + delg223*dginv233 + 
     deldelg2222*ginv23 + deldelg2223*ginv33) + 
  ginv13*(delg111*dginv211 + (delg112 + delg211)*dginv212 + 
     delg212*dginv222 + (delg123 + delg213)*dginv223 + delg133*dginv233 + 
     2.*((delg113 + delg311)*dginv213 + delg312*dginv223 + 
        delg313*dginv233) + deldelg2212*ginv22 + 
     (deldelg1223 + deldelg2213 + 2.*deldelg2312)*ginv23 + 
     (deldelg1233 + 2.*deldelg2313)*ginv33) + 
  ginv23*(delg112*dginv211 + (delg122 + delg212)*dginv212 + 
     (delg123 + delg213 + 2.*delg312)*dginv213 + delg222*dginv222 + 
     delg233*dginv233 + 2.*((delg223 + delg322)*dginv223 + 
        delg323*dginv233) + (deldelg2233 + 2.*deldelg2323)*ginv33) + 
  (deldelg1213 + deldelg2311)*pow2(ginv13) + 
  (deldelg2223 + deldelg2322)*pow2(ginv23) + deldelg2333*pow2(ginv33)
;

dGd31
=
(delg212*dginv311 + delg222*dginv312 + delg223*dginv313)*ginv22 + 
  ((delg213 + delg312)*dginv311 + (delg223 + delg322)*dginv312 + 
     (delg233 + delg323)*dginv313)*ginv23 + 
  (delg313*dginv311 + delg323*dginv312 + delg333*dginv313)*ginv33 + 
  ginv11*(delg211*dginv312 + delg311*dginv313 + 
     2.*(delg111*dginv311 + delg112*dginv312 + delg113*dginv313) + 
     delg212*dginv322 + (delg213 + delg312)*dginv323 + delg313*dginv333 + 
     (2.*deldelg1312 + deldelg2311)*ginv12 + 
     (2.*deldelg1313 + deldelg3311)*ginv13 + deldelg2312*ginv22 + 
     (deldelg2313 + deldelg3312)*ginv23 + deldelg3313*ginv33) + 
  ginv12*((2.*delg112 + delg211)*dginv311 + (delg213 + delg312)*dginv313 + 
     2.*((delg122 + delg212)*dginv312 + delg123*dginv313) + 
     delg222*dginv322 + (delg223 + delg322)*dginv323 + delg323*dginv333 + 
     (2.*deldelg1323 + deldelg2313 + deldelg3312)*ginv13 + 
     deldelg2322*ginv22 + (deldelg2323 + deldelg3322)*ginv23 + 
     deldelg3323*ginv33) + ginv13*
   ((2.*delg113 + delg311)*dginv311 + 
     (2.*delg123 + delg213 + delg312)*dginv312 + 
     2.*(delg133 + delg313)*dginv313 + delg223*dginv322 + 
     (delg233 + delg323)*dginv323 + delg333*dginv333 + deldelg2323*ginv22 + 
     (deldelg2333 + deldelg3323)*ginv23 + deldelg3333*ginv33) + 
  deldelg1311*pow2(ginv11) + (deldelg1322 + deldelg2312)*pow2(ginv12) + 
  (deldelg1333 + deldelg3313)*pow2(ginv13)
;

dGd32
=
ginv11*(delg111*dginv312 + delg112*dginv322 + delg113*dginv323 + 
     deldelg1311*ginv12 + deldelg1312*ginv22 + deldelg1313*ginv23) + 
  ginv13*((delg113 + delg311)*dginv312 + (delg123 + delg312)*dginv322 + 
     (delg133 + delg313)*dginv323 + (deldelg1323 + deldelg3312)*ginv22 + 
     (deldelg1333 + deldelg3313)*ginv23) + 
  (delg313*dginv312 + delg323*dginv322 + delg333*dginv323)*ginv33 + 
  ginv12*(delg111*dginv311 + (delg113 + delg311)*dginv313 + 
     delg122*dginv322 + (delg123 + delg312)*dginv323 + 
     2.*((delg112 + delg211)*dginv312 + delg212*dginv322 + 
        delg213*dginv323) + delg313*dginv333 + 
     (deldelg1313 + deldelg3311)*ginv13 + 
     (deldelg1322 + 2.*deldelg2312)*ginv22 + 
     (deldelg1323 + 2.*deldelg2313 + deldelg3312)*ginv23 + 
     deldelg3313*ginv33) + ginv22*
   (delg112*dginv311 + (delg122 + 2.*delg212)*dginv312 + 
     (delg123 + delg312)*dginv313 + delg322*dginv323 + 
     2.*(delg222*dginv322 + delg223*dginv323) + delg323*dginv333 + 
     (2.*deldelg2323 + deldelg3322)*ginv23 + deldelg3323*ginv33) + 
  ginv23*(delg113*dginv311 + (delg123 + 2.*delg213 + delg312)*dginv312 + 
     (delg133 + delg313)*dginv313 + (2.*delg223 + delg322)*dginv322 + 
     2.*(delg233 + delg323)*dginv323 + delg333*dginv333 + deldelg3333*ginv33\
) + (deldelg1312 + deldelg2311)*pow2(ginv12) + deldelg2322*pow2(ginv22) + 
  (deldelg2333 + deldelg3323)*pow2(ginv23)
;

dGd33
=
(delg113*dginv311 + (delg123 + delg213)*dginv312 + 
     (delg133 + 2.*delg313)*dginv313 + delg223*dginv322 + 
     (delg233 + 2.*delg323)*dginv323 + 2.*delg333*dginv333)*ginv33 + 
  ginv11*(delg111*dginv313 + delg112*dginv323 + delg113*dginv333 + 
     deldelg1311*ginv13 + deldelg1312*ginv23 + deldelg1313*ginv33) + 
  ginv12*((delg112 + delg211)*dginv313 + (delg122 + delg212)*dginv323 + 
     (delg123 + delg213)*dginv333 + (deldelg1312 + deldelg2311)*ginv13 + 
     (deldelg1322 + deldelg2312)*ginv23 + (deldelg1323 + deldelg2313)*ginv33\
) + ginv22*(delg212*dginv313 + delg222*dginv323 + delg223*dginv333 + 
     deldelg2322*ginv23 + deldelg2323*ginv33) + 
  ginv13*(delg111*dginv311 + (delg112 + delg211)*dginv312 + 
     delg212*dginv322 + (delg123 + delg213)*dginv323 + delg133*dginv333 + 
     2.*((delg113 + delg311)*dginv313 + delg312*dginv323 + 
        delg313*dginv333) + deldelg2312*ginv22 + 
     (deldelg1323 + deldelg2313 + 2.*deldelg3312)*ginv23 + 
     (deldelg1333 + 2.*deldelg3313)*ginv33) + 
  ginv23*(delg112*dginv311 + (delg122 + delg212)*dginv312 + 
     (delg123 + delg213 + 2.*delg312)*dginv313 + delg222*dginv322 + 
     delg233*dginv333 + 2.*((delg223 + delg322)*dginv323 + 
        delg323*dginv333) + (deldelg2333 + 2.*deldelg3323)*ginv33) + 
  (deldelg1313 + deldelg3311)*pow2(ginv13) + 
  (deldelg2323 + deldelg3322)*pow2(ginv23) + deldelg3333*pow2(ginv33)
;

Zinv1
=
0.5*(G1 - Gd1)
;

Zinv2
=
0.5*(G2 - Gd2)
;

Zinv3
=
0.5*(G3 - Gd3)
;

dZinv11
=
0.5*(delG11 - dGd11)
;

dZinv12
=
0.5*(delG12 - dGd12)
;

dZinv13
=
0.5*(delG13 - dGd13)
;

dZinv21
=
0.5*(delG21 - dGd21)
;

dZinv22
=
0.5*(delG22 - dGd22)
;

dZinv23
=
0.5*(delG23 - dGd23)
;

dZinv31
=
0.5*(delG31 - dGd31)
;

dZinv32
=
0.5*(delG32 - dGd32)
;

dZinv33
=
0.5*(delG33 - dGd33)
;

Z1
=
g11*Zinv1 + g12*Zinv2 + g13*Zinv3
;

Z2
=
g12*Zinv1 + g22*Zinv2 + g23*Zinv3
;

Z3
=
g13*Zinv1 + g23*Zinv2 + g33*Zinv3
;

dZ11
=
dZinv11*g11 + dZinv12*g12 + dZinv13*g13 + delg111*Zinv1 + delg112*Zinv2 + 
  delg113*Zinv3
;

dZ12
=
dZinv11*g12 + dZinv12*g22 + dZinv13*g23 + delg112*Zinv1 + delg122*Zinv2 + 
  delg123*Zinv3
;

dZ13
=
dZinv11*g13 + dZinv12*g23 + dZinv13*g33 + delg113*Zinv1 + delg123*Zinv2 + 
  delg133*Zinv3
;

dZ21
=
dZinv21*g11 + dZinv22*g12 + dZinv23*g13 + delg211*Zinv1 + delg212*Zinv2 + 
  delg213*Zinv3
;

dZ22
=
dZinv21*g12 + dZinv22*g22 + dZinv23*g23 + delg212*Zinv1 + delg222*Zinv2 + 
  delg223*Zinv3
;

dZ23
=
dZinv21*g13 + dZinv22*g23 + dZinv23*g33 + delg213*Zinv1 + delg223*Zinv2 + 
  delg233*Zinv3
;

dZ31
=
dZinv31*g11 + dZinv32*g12 + dZinv33*g13 + delg311*Zinv1 + delg312*Zinv2 + 
  delg313*Zinv3
;

dZ32
=
dZinv31*g12 + dZinv32*g22 + dZinv33*g23 + delg312*Zinv1 + delg322*Zinv2 + 
  delg323*Zinv3
;

dZ33
=
dZinv31*g13 + dZinv32*g23 + dZinv33*g33 + delg313*Zinv1 + delg323*Zinv2 + 
  delg333*Zinv3
;

DZinv11
=
dZinv11 + gammaF111*Zinv1 + gammaF112*Zinv2 + gammaF113*Zinv3
;

DZinv12
=
dZinv12 + gammaF211*Zinv1 + gammaF212*Zinv2 + gammaF213*Zinv3
;

DZinv13
=
dZinv13 + gammaF311*Zinv1 + gammaF312*Zinv2 + gammaF313*Zinv3
;

DZinv21
=
dZinv21 + gammaF121*Zinv1 + gammaF122*Zinv2 + gammaF123*Zinv3
;

DZinv22
=
dZinv22 + gammaF221*Zinv1 + gammaF222*Zinv2 + gammaF223*Zinv3
;

DZinv23
=
dZinv23 + gammaF321*Zinv1 + gammaF322*Zinv2 + gammaF323*Zinv3
;

DZinv31
=
dZinv31 + gammaF131*Zinv1 + gammaF132*Zinv2 + gammaF133*Zinv3
;

DZinv32
=
dZinv32 + gammaF231*Zinv1 + gammaF232*Zinv2 + gammaF233*Zinv3
;

DZinv33
=
dZinv33 + gammaF331*Zinv1 + gammaF332*Zinv2 + gammaF333*Zinv3
;

DZ11
=
dZ11 - gammaF111*Z1 - gammaF211*Z2 - gammaF311*Z3
;

DZ12
=
dZ12 - gammaF112*Z1 - gammaF212*Z2 - gammaF312*Z3
;

DZ13
=
dZ13 - gammaF113*Z1 - gammaF213*Z2 - gammaF313*Z3
;

DZ21
=
dZ21 - gammaF121*Z1 - gammaF221*Z2 - gammaF321*Z3
;

DZ22
=
dZ22 - gammaF122*Z1 - gammaF222*Z2 - gammaF322*Z3
;

DZ23
=
dZ23 - gammaF123*Z1 - gammaF223*Z2 - gammaF323*Z3
;

DZ31
=
dZ31 - gammaF131*Z1 - gammaF231*Z2 - gammaF331*Z3
;

DZ32
=
dZ32 - gammaF132*Z1 - gammaF232*Z2 - gammaF332*Z3
;

DZ33
=
dZ33 - gammaF133*Z1 - gammaF233*Z2 - gammaF333*Z3
;

DZsym11
=
2.*DZ11
;

DZsym12
=
DZ12 + DZ21
;

DZsym13
=
DZ13 + DZ31
;

DZsym21
=
DZ12 + DZ21
;

DZsym22
=
2.*DZ22
;

DZsym23
=
DZ23 + DZ32
;

DZsym31
=
DZ13 + DZ31
;

DZsym32
=
DZ23 + DZ32
;

DZsym33
=
2.*DZ33
;

trDZsym
=
(DZsym11*ginv11 + (DZsym12 + DZsym21)*ginv12 + (DZsym13 + DZsym31)*ginv13 + 
    DZsym22*ginv22 + (DZsym23 + DZsym32)*ginv23 + DZsym33*ginv33)*psim4
;

rA11
=
rA11 + alpha*(-2.*A11*Theta + chi*
      (DZsym11 - 0.33333333333333333333*g11*trDZsym))
;

rA12
=
rA12 + alpha*(-2.*A12*Theta + chi*
      (DZsym21 - 0.33333333333333333333*g12*trDZsym))
;

rA13
=
rA13 + alpha*(-2.*A13*Theta + chi*
      (DZsym31 - 0.33333333333333333333*g13*trDZsym))
;

rA22
=
rA22 + alpha*(-2.*A22*Theta + chi*
      (DZsym22 - 0.33333333333333333333*g22*trDZsym))
;

rA23
=
rA23 + alpha*(-2.*A23*Theta + chi*
      (DZsym32 - 0.33333333333333333333*g23*trDZsym))
;

rA33
=
rA33 + alpha*(-2.*A33*Theta + chi*
      (DZsym33 - 0.33333333333333333333*g33*trDZsym))
;

rTheta
=
alpha*(DZinv11 + DZinv22 + DZinv33) + rTheta - da1*Zinv1 - da2*Zinv2 - 
  da3*Zinv3
;

rG1
=
rG1 - ginv11*(1.3333333333333333333*alpha*K*Z1 + 
     2.*(da1*Theta + alpha*kappa1*Z1)) - 
  ginv12*(1.3333333333333333333*alpha*K*Z2 + 
     2.*(da2*Theta + alpha*kappa1*Z2)) - 
  ginv13*(1.3333333333333333333*alpha*K*Z3 + 
     2.*(da3*Theta + alpha*kappa1*Z3))
;

rG2
=
rG2 - ginv12*(1.3333333333333333333*alpha*K*Z1 + 
     2.*(da1*Theta + alpha*kappa1*Z1)) - 
  ginv22*(1.3333333333333333333*alpha*K*Z2 + 
     2.*(da2*Theta + alpha*kappa1*Z2)) - 
  ginv23*(1.3333333333333333333*alpha*K*Z3 + 
     2.*(da3*Theta + alpha*kappa1*Z3))
;

rG3
=
rG3 - ginv13*(1.3333333333333333333*alpha*K*Z1 + 
     2.*(da1*Theta + alpha*kappa1*Z1)) - 
  ginv23*(1.3333333333333333333*alpha*K*Z2 + 
     2.*(da2*Theta + alpha*kappa1*Z2)) - 
  ginv33*(1.3333333333333333333*alpha*K*Z3 + 
     2.*(da3*Theta + alpha*kappa1*Z3))
;
#endif

  } /* function */
}
