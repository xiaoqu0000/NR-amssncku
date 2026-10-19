

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


#define Power(x,y) (pow((double) (x), (double) (y)))
#define Sqrt(x)    sqrt(x)
#define Log(x)     log((double) (x))
#define pow2(x)    ((x)*(x))
#define pow3(x)    ((x)*(x)*(x))
#define pow4(x)    ((x)*(x)*(x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))

#define Cal(x,y,z) ((x)?(y):(z))

#define Tan(x)     tan(x)
#define ArcTan(x)  atan(x)
#define Sin(x)     sin(x)
#define Cos(x)     cos(x)
#define Csc(x)     (1./sin(x))
#define Abs(x)     (fabs(x))
#define sqrt2      (sqrt(2))
#define Tanh(x)    tanh(x)
#define Sech(x)    (1/cosh(x))


extern "C" {

#ifdef fortran1
void cpbc_point
#endif	
#ifdef fortran2
void CPBC_POINT
#endif
#ifdef fortran3
void cpbc_point_
#endif
(double & r,double & xp,double & yp,double & zp,
		 double & Theta,double & chi,double & Khat,
		 double & g11,double & g12,double & g13,double & g22,double & g23,double & g33,
		 double & A11,double & A12,double & A13,double & A22,double & A23,double & A33,
		 double & G1,double & G2,double & G3,
		 double & alpha,double & beta1,double & beta2,double & beta3,
		 double & da1,double & da2,double & da3,
		 double & dda11,double & dda12,double & dda13,double & dda22,double & dda23,double & dda33,
		 double & db11,double & db21,double & db31,
		 double & db12,double & db22,double & db32,
		 double & db13,double & db23,double & db33,
		 double & ddb111,double & ddb121,double & ddb131,double & ddb221,double & ddb231,double & ddb331,
		 double & ddb112,double & ddb122,double & ddb132,double & ddb222,double & ddb232,double & ddb332,
		 double & ddb113,double & ddb123,double & ddb133,double & ddb223,double & ddb233,double & ddb333,
		 double & dchi1,double & dchi2,double & dchi3,
		 double & ddchi11,double & ddchi12,double & ddchi13,double & ddchi22,double & ddchi23,double & ddchi33,
		 double & dg111,double & dg112,double & dg113,double & dg122,double & dg123,double & dg133,
		 double & dg211,double & dg212,double & dg213,double & dg222,double & dg223,double & dg233,
		 double & dg311,double & dg312,double & dg313,double & dg322,double & dg323,double & dg333,
		 double & ddg1111,double & ddg1211,double & ddg1311,double & ddg2211,double & ddg2311,double & ddg3311,
		 double & ddg1112,double & ddg1212,double & ddg1312,double & ddg2212,double & ddg2312,double & ddg3312,
		 double & ddg1113,double & ddg1213,double & ddg1313,double & ddg2213,double & ddg2313,double & ddg3313,
		 double & ddg1122,double & ddg1222,double & ddg1322,double & ddg2222,double & ddg2322,double & ddg3322,
		 double & ddg1123,double & ddg1223,double & ddg1323,double & ddg2223,double & ddg2323,double & ddg3323,
		 double & ddg1133,double & ddg1233,double & ddg1333,double & ddg2233,double & ddg2333,double & ddg3333,
		 double & dKhat1,double & dKhat2,double & dKhat3,
		 double & dA111,double & dA112,double & dA113,double & dA122,double & dA123,double & dA133,
		 double & dA211,double & dA212,double & dA213,double & dA222,double & dA223,double & dA233,
		 double & dA311,double & dA312,double & dA313,double & dA322,double & dA323,double & dA333,
		 double & dG11,double & dG21,double & dG31,
		 double & dG12,double & dG22,double & dG32,
		 double & dG13,double & dG23,double & dG33,
		 double & dTheta1,double & dTheta2,double & dTheta3,
		 double & rKhat,double & rTheta,
		 double & rA11,double & rA12,double & rA13,double & rA22,double & rA23,double & rA33,
		 double & rG1,double & rG2,double & rG3,
		 double &kappa1,double &kappa2,double &shiftdriver)
{

double AA11;
double AA12;
double AA13;
double AA21;
double AA22;
double AA23;
double AA31;
double AA32;
double AA33;
double ADMginv11;
double ADMginv12;
double ADMginv13;
double ADMginv22;
double ADMginv23;
double ADMginv33;
double Ainv11;
double Ainv12;
double Ainv13;
double Ainv22;
double Ainv23;
double Ainv33;
double betaA1;
double betaA2;
double betaA3;
double betas;
double cdA111;
double cdA112;
double cdA113;
double cdA122;
double cdA123;
double cdA133;
double cdA211;
double cdA212;
double cdA213;
double cdA222;
double cdA223;
double cdA233;
double cdA311;
double cdA312;
double cdA313;
double cdA322;
double cdA323;
double cdA333;
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
const double chipsipower = -4;
double Dalpha;
double DbetaA1;
double DbetaA2;
double DbetaA3;
double Dbetas;
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
double DGamA1;
double DGamA2;
double DGamA3;
double DGams;
double dGfromgdu11;
double dGfromgdu12;
double dGfromgdu13;
double dGfromgdu21;
double dGfromgdu22;
double dGfromgdu23;
double dGfromgdu31;
double dGfromgdu32;
double dGfromgdu33;
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
double divbeta;
double DK;
double dK1;
double dK2;
double dK3;
double DKhat;
double DTheta;
double f;
double ff;
double gADM11;
double gADM12;
double gADM13;
double gADM21;
double gADM22;
double gADM23;
double gADM31;
double gADM32;
double gADM33;
double GamA1;
double GamA2;
double GamA3;
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
double Gams;
double Gfromg1;
double Gfromg2;
double Gfromg3;
double ginv11;
double ginv12;
double ginv13;
double ginv22;
double ginv23;
double ginv33;
const bool givehPsi0 = false;
const double hPsi0para = 0;
const double hPsi0parb = 0;
const double hPsi0parc = 0;
double ImhPsi0;
double K;
double lieA11;
double lieA12;
double lieA13;
double lieA22;
double lieA23;
double lieA33;
double lieg11;
double lieg12;
double lieg13;
double lieg22;
double lieg23;
double lieg33;
double lienK;
double lienKhat;
double lienTheta;
double modshatARG;
double muL;
double muStilde;
double oochipsipower;
double oomodshat;
double psim4;
double qdd11;
double qdd12;
double qdd13;
double qdd22;
double qdd23;
double qdd33;
double qPhysuudd1111;
double qPhysuudd1112;
double qPhysuudd1113;
double qPhysuudd1122;
double qPhysuudd1123;
double qPhysuudd1133;
double qPhysuudd1211;
double qPhysuudd1212;
double qPhysuudd1213;
double qPhysuudd1222;
double qPhysuudd1223;
double qPhysuudd1233;
double qPhysuudd1311;
double qPhysuudd1312;
double qPhysuudd1313;
double qPhysuudd1322;
double qPhysuudd1323;
double qPhysuudd1333;
double qPhysuudd2211;
double qPhysuudd2212;
double qPhysuudd2213;
double qPhysuudd2222;
double qPhysuudd2223;
double qPhysuudd2233;
double qPhysuudd2311;
double qPhysuudd2312;
double qPhysuudd2313;
double qPhysuudd2322;
double qPhysuudd2323;
double qPhysuudd2333;
double qPhysuudd3311;
double qPhysuudd3312;
double qPhysuudd3313;
double qPhysuudd3322;
double qPhysuudd3323;
double qPhysuudd3333;
double qud11;
double qud12;
double qud13;
double qud21;
double qud22;
double qud23;
double qud31;
double qud32;
double qud33;
double quu11;
double quu12;
double quu13;
double quu22;
double quu23;
double quu33;
double R11;
double R12;
double R13;
double R22;
double R23;
double R33;
double rACABTF11;
double rACABTF12;
double rACABTF13;
double rACABTF22;
double rACABTF23;
double rACABTF33;
double rACqq;
double rACsA1;
double rACsA2;
double rACsA3;
double rACss;
double RehPsi0;
double Rf11;
double Rf12;
double Rf13;
double Rf22;
double Rf23;
double Rf33;
double rGamA1;
double rGamA2;
double rGamA3;
double rGams;
double Rhat;
double Rphi11;
double Rphi12;
double Rphi13;
double Rphi22;
double Rphi23;
double Rphi33;
double sdotv;
double sdotw;
double sdown1;
double sdown2;
double sdown3;
double shat1;
double shat2;
double shat3;
double sup1;
double sup2;
double sup3;
const double time = 0;
double totdivbeta;
double trcdda;
double trcddf;
double vbetaA;
double vbetas;
double vd1;
double vd2;
double vd3;
double vdotv;
double vdotw;
double vu1;
double vu2;
double vu3;
double wd1;
double wd2;
double wd3;
double wdotw;
double wu1;
double wu2;
double wu3;

shat1=xp/r;shat2=yp/r;shat3=zp/r; 

#if 0
// my code
detginv
=
1/(2.*g12*g13*g23 - g33*pow2(g12) + g22*(g11*g33 - pow2(g13)) - 
    g11*pow2(g23))
;

ginv11
=
detginv*(g22*g33 - pow2(g23))
;

ginv12
=
detginv*(g13*g23 - g12*g33)
;

ginv13
=
detginv*(-(g13*g22) + g12*g23)
;

ginv22
=
detginv*(g11*g33 - pow2(g13))
;

ginv23
=
detginv*(g12*g13 - g11*g23)
;

ginv33
=
detginv*(g11*g22 - pow2(g12))
;

ADMginv11
=
chi*ginv11
;

ADMginv12
=
chi*ginv12
;

ADMginv13
=
chi*ginv13
;

ADMginv22
=
chi*ginv22
;

ADMginv23
=
chi*ginv23
;

ADMginv33
=
chi*ginv33
;

modshatARG
=
2.*(ADMginv23*shat2*shat3 + shat1*(ADMginv12*shat2 + ADMginv13*shat3)) + 
  ADMginv11*pow2(shat1) + ADMginv22*pow2(shat2) + ADMginv33*pow2(shat3)
;


if (modshatARG<0.00001) {                           
      printf("modshat is wrong (%e)\n",modshatARG);
      modshatARG = 0.00001;
    }oomodshat
=
1/sqrt(modshatARG)
;

sdown1
=
oomodshat*shat1
;

sdown2
=
oomodshat*shat2
;

sdown3
=
oomodshat*shat3
;

sup1
=
ADMginv11*sdown1 + ADMginv12*sdown2 + ADMginv13*sdown3
;

sup2
=
ADMginv12*sdown1 + ADMginv22*sdown2 + ADMginv23*sdown3
;

sup3
=
ADMginv13*sdown1 + ADMginv23*sdown2 + ADMginv33*sdown3
;

qud11
=
1. - sdown1*sup1
;

qud12
=
-(sdown2*sup1)
;

qud13
=
-(sdown3*sup1)
;

qud21
=
-(sdown1*sup2)
;

qud22
=
1. - sdown2*sup2
;

qud23
=
-(sdown3*sup2)
;

qud31
=
-(sdown1*sup3)
;

qud32
=
-(sdown2*sup3)
;

qud33
=
1. - sdown3*sup3
;

qdd11
=
g11/chi - pow2(sdown1)
;

qdd12
=
g12/chi - sdown1*sdown2
;

qdd13
=
g13/chi - sdown1*sdown3
;

qdd22
=
g22/chi - pow2(sdown2)
;

qdd23
=
g23/chi - sdown2*sdown3
;

qdd33
=
g33/chi - pow2(sdown3)
;

quu11
=
ADMginv11 - pow2(sup1)
;

quu12
=
ADMginv12 - sup1*sup2
;

quu13
=
ADMginv13 - sup1*sup3
;

quu22
=
ADMginv22 - pow2(sup2)
;

quu23
=
ADMginv23 - sup2*sup3
;

quu33
=
ADMginv33 - pow2(sup3)
;

qPhysuudd1111
=
-0.5*qdd11*quu11 + pow2(qud11)
;

qPhysuudd1112
=
qud11*qud12 - 0.5*qdd12*quu11
;

qPhysuudd1113
=
qud11*qud13 - 0.5*qdd13*quu11
;

qPhysuudd1122
=
-0.5*qdd22*quu11 + pow2(qud12)
;

qPhysuudd1123
=
qud12*qud13 - 0.5*qdd23*quu11
;

qPhysuudd1133
=
-0.5*qdd33*quu11 + pow2(qud13)
;

qPhysuudd1211
=
qud11*qud21 - 0.5*qdd11*quu12
;

qPhysuudd1212
=
0.5*(qud12*qud21 + qud11*qud22 - qdd12*quu12)
;

qPhysuudd1213
=
0.5*(qud13*qud21 + qud11*qud23 - qdd13*quu12)
;

qPhysuudd1222
=
qud12*qud22 - 0.5*qdd22*quu12
;

qPhysuudd1223
=
0.5*(qud13*qud22 + qud12*qud23 - qdd23*quu12)
;

qPhysuudd1233
=
qud13*qud23 - 0.5*qdd33*quu12
;

qPhysuudd1311
=
qud11*qud31 - 0.5*qdd11*quu13
;

qPhysuudd1312
=
0.5*(qud12*qud31 + qud11*qud32 - qdd12*quu13)
;

qPhysuudd1313
=
0.5*(qud13*qud31 + qud11*qud33 - qdd13*quu13)
;

qPhysuudd1322
=
qud12*qud32 - 0.5*qdd22*quu13
;

qPhysuudd1323
=
0.5*(qud13*qud32 + qud12*qud33 - qdd23*quu13)
;

qPhysuudd1333
=
qud13*qud33 - 0.5*qdd33*quu13
;

qPhysuudd2211
=
-0.5*qdd11*quu22 + pow2(qud21)
;

qPhysuudd2212
=
qud21*qud22 - 0.5*qdd12*quu22
;

qPhysuudd2213
=
qud21*qud23 - 0.5*qdd13*quu22
;

qPhysuudd2222
=
-0.5*qdd22*quu22 + pow2(qud22)
;

qPhysuudd2223
=
qud22*qud23 - 0.5*qdd23*quu22
;

qPhysuudd2233
=
-0.5*qdd33*quu22 + pow2(qud23)
;

qPhysuudd2311
=
qud21*qud31 - 0.5*qdd11*quu23
;

qPhysuudd2312
=
0.5*(qud22*qud31 + qud21*qud32 - qdd12*quu23)
;

qPhysuudd2313
=
0.5*(qud23*qud31 + qud21*qud33 - qdd13*quu23)
;

qPhysuudd2322
=
qud22*qud32 - 0.5*qdd22*quu23
;

qPhysuudd2323
=
0.5*(qud23*qud32 + qud22*qud33 - qdd23*quu23)
;

qPhysuudd2333
=
qud23*qud33 - 0.5*qdd33*quu23
;

qPhysuudd3311
=
-0.5*qdd11*quu33 + pow2(qud31)
;

qPhysuudd3312
=
qud31*qud32 - 0.5*qdd12*quu33
;

qPhysuudd3313
=
qud31*qud33 - 0.5*qdd13*quu33
;

qPhysuudd3322
=
-0.5*qdd22*quu33 + pow2(qud32)
;

qPhysuudd3323
=
qud32*qud33 - 0.5*qdd23*quu33
;

qPhysuudd3333
=
-0.5*qdd33*quu33 + pow2(qud33)
;

muL
=
2./alpha
;

muStilde
=
1/chi
;

vbetas
=
2.*sqrt(0.33333333333333333333*muStilde)
;

vbetaA
=
sqrt(muStilde)
;

K
=
Khat + 2.*Theta
;

dK1
=
dKhat1 + 2.*dTheta1
;

dK2
=
dKhat2 + 2.*dTheta2
;

dK3
=
dKhat3 + 2.*dTheta3
;

dginv111
=
-2.*(dg123*ginv12*ginv13 + ginv11*(dg112*ginv12 + dg113*ginv13)) - 
  dg111*pow2(ginv11) - dg122*pow2(ginv12) - dg133*pow2(ginv13)
;

dginv112
=
-(ginv11*(dg111*ginv12 + dg112*ginv22 + dg113*ginv23)) - 
  ginv12*(dg113*ginv13 + dg122*ginv22 + dg123*ginv23) - 
  ginv13*(dg123*ginv22 + dg133*ginv23) - dg112*pow2(ginv12)
;

dginv113
=
-(ginv11*(dg111*ginv13 + dg112*ginv23 + dg113*ginv33)) - 
  ginv12*(dg112*ginv13 + dg122*ginv23 + dg123*ginv33) - 
  ginv13*(dg123*ginv23 + dg133*ginv33) - dg113*pow2(ginv13)
;

dginv122
=
-2.*(dg123*ginv22*ginv23 + ginv12*(dg112*ginv22 + dg113*ginv23)) - 
  dg111*pow2(ginv12) - dg122*pow2(ginv22) - dg133*pow2(ginv23)
;

dginv123
=
-(ginv13*(dg112*ginv22 + dg113*ginv23)) - dg133*ginv23*ginv33 - 
  ginv12*(dg111*ginv13 + dg112*ginv23 + dg113*ginv33) - 
  ginv22*(dg122*ginv23 + dg123*ginv33) - dg123*pow2(ginv23)
;

dginv133
=
-2.*(dg123*ginv23*ginv33 + ginv13*(dg112*ginv23 + dg113*ginv33)) - 
  dg111*pow2(ginv13) - dg122*pow2(ginv23) - dg133*pow2(ginv33)
;

dginv211
=
-2.*(dg223*ginv12*ginv13 + ginv11*(dg212*ginv12 + dg213*ginv13)) - 
  dg211*pow2(ginv11) - dg222*pow2(ginv12) - dg233*pow2(ginv13)
;

dginv212
=
-(ginv11*(dg211*ginv12 + dg212*ginv22 + dg213*ginv23)) - 
  ginv12*(dg213*ginv13 + dg222*ginv22 + dg223*ginv23) - 
  ginv13*(dg223*ginv22 + dg233*ginv23) - dg212*pow2(ginv12)
;

dginv213
=
-(ginv11*(dg211*ginv13 + dg212*ginv23 + dg213*ginv33)) - 
  ginv12*(dg212*ginv13 + dg222*ginv23 + dg223*ginv33) - 
  ginv13*(dg223*ginv23 + dg233*ginv33) - dg213*pow2(ginv13)
;

dginv222
=
-2.*(dg223*ginv22*ginv23 + ginv12*(dg212*ginv22 + dg213*ginv23)) - 
  dg211*pow2(ginv12) - dg222*pow2(ginv22) - dg233*pow2(ginv23)
;

dginv223
=
-(ginv13*(dg212*ginv22 + dg213*ginv23)) - dg233*ginv23*ginv33 - 
  ginv12*(dg211*ginv13 + dg212*ginv23 + dg213*ginv33) - 
  ginv22*(dg222*ginv23 + dg223*ginv33) - dg223*pow2(ginv23)
;

dginv233
=
-2.*(dg223*ginv23*ginv33 + ginv13*(dg212*ginv23 + dg213*ginv33)) - 
  dg211*pow2(ginv13) - dg222*pow2(ginv23) - dg233*pow2(ginv33)
;

dginv311
=
-2.*(dg323*ginv12*ginv13 + ginv11*(dg312*ginv12 + dg313*ginv13)) - 
  dg311*pow2(ginv11) - dg322*pow2(ginv12) - dg333*pow2(ginv13)
;

dginv312
=
-(ginv11*(dg311*ginv12 + dg312*ginv22 + dg313*ginv23)) - 
  ginv12*(dg313*ginv13 + dg322*ginv22 + dg323*ginv23) - 
  ginv13*(dg323*ginv22 + dg333*ginv23) - dg312*pow2(ginv12)
;

dginv313
=
-(ginv11*(dg311*ginv13 + dg312*ginv23 + dg313*ginv33)) - 
  ginv12*(dg312*ginv13 + dg322*ginv23 + dg323*ginv33) - 
  ginv13*(dg323*ginv23 + dg333*ginv33) - dg313*pow2(ginv13)
;

dginv322
=
-2.*(dg323*ginv22*ginv23 + ginv12*(dg312*ginv22 + dg313*ginv23)) - 
  dg311*pow2(ginv12) - dg322*pow2(ginv22) - dg333*pow2(ginv23)
;

dginv323
=
-(ginv13*(dg312*ginv22 + dg313*ginv23)) - dg333*ginv23*ginv33 - 
  ginv12*(dg311*ginv13 + dg312*ginv23 + dg313*ginv33) - 
  ginv22*(dg322*ginv23 + dg323*ginv33) - dg323*pow2(ginv23)
;

dginv333
=
-2.*(dg323*ginv23*ginv33 + ginv13*(dg312*ginv23 + dg313*ginv33)) - 
  dg311*pow2(ginv13) - dg322*pow2(ginv23) - dg333*pow2(ginv33)
;

gammado111
=
0.5*dg111
;

gammado112
=
0.5*dg211
;

gammado113
=
0.5*dg311
;

gammado122
=
-0.5*dg122 + dg212
;

gammado123
=
0.5*(-dg123 + dg213 + dg312)
;

gammado133
=
-0.5*dg133 + dg313
;

gammado211
=
dg112 - 0.5*dg211
;

gammado212
=
0.5*dg122
;

gammado213
=
0.5*(dg123 - dg213 + dg312)
;

gammado222
=
0.5*dg222
;

gammado223
=
0.5*dg322
;

gammado233
=
-0.5*dg233 + dg323
;

gammado311
=
dg113 - 0.5*dg311
;

gammado312
=
0.5*(dg123 + dg213 - dg312)
;

gammado313
=
0.5*dg133
;

gammado322
=
dg223 - 0.5*dg322
;

gammado323
=
0.5*dg233
;

gammado333
=
0.5*dg333
;

gamma111
=
gammado111*ginv11 + gammado211*ginv12 + gammado311*ginv13
;

gamma112
=
gammado112*ginv11 + gammado212*ginv12 + gammado312*ginv13
;

gamma113
=
gammado113*ginv11 + gammado213*ginv12 + gammado313*ginv13
;

gamma122
=
gammado122*ginv11 + gammado222*ginv12 + gammado322*ginv13
;

gamma123
=
gammado123*ginv11 + gammado223*ginv12 + gammado323*ginv13
;

gamma133
=
gammado133*ginv11 + gammado233*ginv12 + gammado333*ginv13
;

gamma211
=
gammado111*ginv12 + gammado211*ginv22 + gammado311*ginv23
;

gamma212
=
gammado112*ginv12 + gammado212*ginv22 + gammado312*ginv23
;

gamma213
=
gammado113*ginv12 + gammado213*ginv22 + gammado313*ginv23
;

gamma222
=
gammado122*ginv12 + gammado222*ginv22 + gammado322*ginv23
;

gamma223
=
gammado123*ginv12 + gammado223*ginv22 + gammado323*ginv23
;

gamma233
=
gammado133*ginv12 + gammado233*ginv22 + gammado333*ginv23
;

gamma311
=
gammado111*ginv13 + gammado211*ginv23 + gammado311*ginv33
;

gamma312
=
gammado112*ginv13 + gammado212*ginv23 + gammado312*ginv33
;

gamma313
=
gammado113*ginv13 + gammado213*ginv23 + gammado313*ginv33
;

gamma322
=
gammado122*ginv13 + gammado222*ginv23 + gammado322*ginv33
;

gamma323
=
gammado123*ginv13 + gammado223*ginv23 + gammado323*ginv33
;

gamma333
=
gammado133*ginv13 + gammado233*ginv23 + gammado333*ginv33
;

Gfromg1
=
gamma111*ginv11 + gamma122*ginv22 + 
  2.*(gamma112*ginv12 + gamma113*ginv13 + gamma123*ginv23) + gamma133*ginv33
;

Gfromg2
=
gamma211*ginv11 + gamma222*ginv22 + 
  2.*(gamma212*ginv12 + gamma213*ginv13 + gamma223*ginv23) + gamma233*ginv33
;

Gfromg3
=
gamma311*ginv11 + gamma322*ginv22 + 
  2.*(gamma312*ginv12 + gamma313*ginv13 + gamma323*ginv23) + gamma333*ginv33
;

dGfromgdu11
=
-((dg122*(4.*dg112 + dg211) + 2.*dg112*dg212 + dg111*dg222)*
     Power(ginv12,3)) - (dg133*(4.*dg113 + dg311) + 2.*dg113*dg313 + 
     dg111*dg333)*Power(ginv13,3) - 2.*Power(ginv11,3)*pow2(dg111) + 
  (ddg1111 - dg111*((8.*dg112 + 2.*dg211)*ginv12 + 
        (8.*dg113 + 2.*dg311)*ginv13) - 
     (dg113*(4.*dg112 + dg211) + dg112*dg311 + dg111*(dg213 + dg312))*
      ginv23 - ginv22*(dg112*dg211 + dg111*dg212 + 2.*pow2(dg112)) - 
     ginv33*(dg113*dg311 + dg111*dg313 + 2.*pow2(dg113)))*pow2(ginv11) + 
  (ddg1122 + ddg1212 - (dg123*(8.*dg112 + 2.*dg211) + 
        dg113*(4.*dg122 + 2.*dg212) + dg122*dg311 + 
        2.*(dg111*dg223 + dg112*(dg213 + dg312)) + dg111*dg322)*ginv13 - 
     (dg123*(4.*dg122 + 2.*dg212) + 
        2.*(dg113*dg222 + dg122*(dg213 + dg312) + dg112*(dg223 + dg322)))*
      ginv23 - ginv22*(3.*(dg122*dg212 + dg112*dg222) + 2.*pow2(dg122)) - 
     ginv33*(dg123*(dg213 + dg312) + dg122*dg313 + dg113*(dg223 + dg322) + 
        dg112*dg323 + 2.*pow2(dg123)))*pow2(ginv12) + 
  (ddg1133 + ddg1313 - (dg133*(4.*dg123 + 2.*(dg213 + dg312)) + 
        2.*(dg123*dg313 + dg113*(dg233 + dg323) + dg112*dg333))*ginv23 - 
     ginv22*(dg133*dg212 + dg113*dg223 + dg123*(dg213 + dg312) + 
        dg112*(dg233 + dg323) + 2.*pow2(dg123)) - 
     ginv33*(3.*(dg133*dg313 + dg113*dg333) + 2.*pow2(dg133)))*pow2(ginv13) \
+ ginv13*(ddg1333*ginv33 + ginv22*
      (ddg1223 - (dg133*dg222 + dg123*(4.*dg223 + dg322) + 
           dg122*(dg233 + dg323))*ginv23 - 
        (dg133*dg223 + dg123*(dg233 + 2.*dg323))*ginv33) + 
     ginv23*(ddg1233 + ddg1323 - 
        (dg133*(2.*dg233 + 3.*dg323) + 3.*dg123*dg333)*ginv33) - 
     (dg123*dg222 + dg122*dg223)*pow2(ginv22) - 
     (dg133*dg322 + 2.*(dg133*dg223 + dg123*(dg233 + dg323)) + 
        dg122*dg333)*pow2(ginv23) - 2.*dg133*dg333*pow2(ginv33)) + 
  ginv11*(ddg1313*ginv33 + ginv12*
      (2.*ddg1112 + ddg1211 - 
        (dg113*(12.*dg112 + 3.*dg211) + 3.*dg112*dg311 + 
           dg111*(8.*dg123 + 3.*(dg213 + dg312)))*ginv13 - 
        (dg122*(4.*dg112 + dg211) + 6.*dg112*dg212 + dg111*dg222)*ginv22 - 
        (dg123*dg211 + dg122*dg311 + 
           4.*(dg113*(dg122 + dg212) + dg112*(dg123 + dg213 + dg312)) + 
           dg111*(dg223 + dg322))*ginv23 - 
        (dg123*dg311 + dg113*(4.*dg123 + 2.*(dg213 + dg312)) + 
           2.*dg112*dg313 + dg111*dg323)*ginv33) + 
     ginv22*(ddg1212 - (dg113*dg222 + 2.*(dg123*dg212 + dg112*dg223) + 
           dg122*(dg213 + dg312) + dg112*dg322)*ginv23 - 
        (dg113*dg223 + dg123*(dg213 + dg312) + dg112*dg323)*ginv33) + 
     ginv13*(2.*ddg1113 + ddg1311 - 
        (dg123*(4.*dg112 + dg211) + dg111*dg223 + 
           2.*(dg113*dg212 + dg112*(dg213 + dg312)))*ginv22 - 
        (dg133*dg211 + dg123*dg311 + 
           4.*(dg113*(dg123 + dg213 + dg312) + dg112*(dg133 + dg313)) + 
           dg111*(dg233 + dg323))*ginv23 - 
        (dg133*(4.*dg113 + dg311) + 6.*dg113*dg313 + dg111*dg333)*ginv33) + 
     ginv23*(ddg1213 + ddg1312 - 
        (dg133*(dg213 + dg312) + 2.*dg123*dg313 + 
           dg113*(dg233 + 2.*dg323) + dg112*dg333)*ginv33) - 
     (3.*dg112*dg211 + dg111*(4.*dg122 + 3.*dg212) + 6.*pow2(dg112))*
      pow2(ginv12) - (3.*dg113*dg311 + dg111*(4.*dg133 + 3.*dg313) + 
        6.*pow2(dg113))*pow2(ginv13) - 
     (dg122*dg212 + dg112*dg222)*pow2(ginv22) - 
     (dg133*dg212 + dg123*(dg213 + dg312) + dg122*dg313 + 
        dg113*(dg223 + dg322) + dg112*(dg233 + dg323))*pow2(ginv23) - 
     (dg133*dg313 + dg113*dg333)*pow2(ginv33)) + 
  ginv12*(ddg1323*ginv33 + ginv22*
      (ddg1222 - (3.*(dg123*dg222 + dg122*dg223) + 2.*dg122*dg322)*ginv23 - 
        (dg123*(2.*dg223 + dg322) + dg122*dg323)*ginv33) + 
     ginv23*(ddg1223 + ddg1322 - 
        (dg133*(dg223 + dg322) + dg123*(dg233 + 4.*dg323) + dg122*dg333)*
         ginv33) + ginv13*(2.*ddg1123 + ddg1213 + ddg1312 - 
        (dg113*dg222 + 4.*(dg123*(dg122 + dg212) + dg112*dg223) + 
           dg122*(dg213 + dg312) + dg112*dg322)*ginv22 - 
        (dg133*(4.*dg123 + dg213 + dg312) + 4.*dg123*dg313 + 
           dg113*(dg233 + 4.*dg323) + dg112*dg333)*ginv33 - 
        ginv23*(2.*(dg133*dg212 + dg112*dg233 + dg122*dg313 + 
              dg113*dg322) + 4.*
            (dg122*dg133 + dg113*dg223 + dg123*(dg213 + dg312) + 
              dg112*dg323 + pow2(dg123)))) - 
     (dg133*(4.*dg112 + dg211) + dg113*(8.*dg123 + 2.*(dg213 + dg312)) + 
        2.*(dg123*dg311 + dg112*dg313) + dg111*(dg233 + 2.*dg323))*
      pow2(ginv13) - 2.*dg122*dg222*pow2(ginv22) - 
     (dg133*dg222 + 2.*dg123*(dg223 + dg322) + dg122*(dg233 + 2.*dg323))*
      pow2(ginv23) - (dg133*dg323 + dg123*dg333)*pow2(ginv33))
;

dGfromgdu12
=
-((dg133*dg322 + 2.*(dg133*dg223 + dg123*(dg233 + dg323)) + dg122*dg333)*
     Power(ginv23,3)) - 2.*(dg122*dg222*Power(ginv22,3) + 
     Power(ginv12,3)*(dg112*dg211 + dg111*(dg122 + dg212) + pow2(dg112)) + 
     (dg111*(dg112*ginv22 + dg113*ginv23) + ginv12*pow2(dg111))*pow2(ginv11)\
) + (ddg1112 + ddg1211 - (4.*(dg112*dg113 + dg111*dg123) + 
        2.*(dg113*dg211 + dg112*dg311 + dg111*(dg213 + dg312)))*ginv13 - 
     (dg122*(6.*dg112 + 2.*dg211) + 6.*dg112*dg212 + 2.*dg111*dg222)*
      ginv22 - (4.*(dg113*(dg122 + dg212) + dg112*(dg123 + dg213)) + 
        dg122*dg311 + 2.*(dg123*dg211 + dg111*dg223 + dg112*dg312) + 
        dg111*dg322)*ginv23 - 
     (dg123*dg311 + dg113*(2.*(dg123 + dg213) + dg312) + dg112*dg313 + 
        dg111*dg323)*ginv33)*pow2(ginv12) - 
  ((2.*(dg113*dg123 + dg112*dg133) + dg123*dg311 + dg113*dg312 + 
        dg112*dg313 + dg111*dg323)*ginv22 + 
     (dg133*(4.*dg113 + dg311) + 2.*dg113*dg313 + dg111*dg333)*ginv23)*
   pow2(ginv13) + (ddg1222 - (4.*(dg123*dg222 + dg122*dg223) + 
        2.*dg122*dg322)*ginv23 - 
     (dg123*(2.*dg223 + dg322) + dg122*dg323)*ginv33)*pow2(ginv22) + 
  (ddg1233 + ddg1323 - (dg133*(2.*dg233 + 3.*dg323) + 3.*dg123*dg333)*
      ginv33)*pow2(ginv23) + ginv11*
   (ginv23*(ddg1113 - 2.*dg113*(dg133 + dg313)*ginv33) + 
     ginv22*(ddg1112 - (dg112*(4.*dg123 + 2.*dg213) + 
           2.*(dg113*(dg122 + dg212) + dg112*dg312))*ginv23 - 
        (dg113*(2.*dg123 + dg312) + dg112*dg313)*ginv33) + 
     ginv12*(ddg1111 - dg111*(6.*dg113 + 2.*dg311)*ginv13 - 
        (dg113*(8.*dg112 + 2.*dg211) + dg112*dg311 + 
           dg111*(2.*(dg123 + dg213) + dg312))*ginv23 - 
        ginv22*(2.*(dg112*dg211 + dg111*(dg122 + dg212)) + 
           6.*pow2(dg112)) - ginv33*
         (dg113*dg311 + dg111*dg313 + 2.*pow2(dg113))) - 
     ginv13*((dg112*(4.*dg113 + dg311) + dg111*(2.*dg123 + dg312))*
         ginv22 + ginv23*(dg113*dg311 + dg111*(2.*dg133 + dg313) + 
           4.*pow2(dg113))) - dg111*(6.*dg112 + 2.*dg211)*pow2(ginv12) - 
     2.*dg112*(dg122 + dg212)*pow2(ginv22) - 
     (2.*(dg112*dg133 + dg113*(dg123 + dg213)) + dg113*dg312 + dg112*dg313)*
      pow2(ginv23)) + ginv13*(ginv22*
      (ddg1123 + ddg1312 - (dg133*(2.*dg123 + dg312) + 
           2.*(dg123*dg313 + dg113*dg323) + dg112*dg333)*ginv33 - 
        ginv23*(2.*(dg133*(dg122 + dg212) + dg123*dg213 + dg113*dg223 + 
              dg112*dg233) + dg122*dg313 + dg113*dg322 + 
           4.*(dg123*dg312 + dg112*dg323 + pow2(dg123)))) + 
     ginv23*(ddg1133 + ddg1313 - 
        ginv33*(3.*(dg133*dg313 + dg113*dg333) + 2.*pow2(dg133))) - 
     (2.*(dg123*(dg122 + dg212) + dg112*dg223) + dg122*dg312 + 
        dg112*dg322)*pow2(ginv22) - 
     (dg133*(4.*dg123 + 2.*(dg213 + dg312)) + 
        2.*(dg123*dg313 + dg113*(dg233 + dg323) + dg112*dg333))*pow2(ginv23)\
) + ginv23*(ddg1333*ginv33 - 2.*dg133*dg333*pow2(ginv33)) + 
  ginv12*(ddg1313*ginv33 + ginv13*
      (ddg1113 + ddg1311 - (2.*
            (dg123*dg211 + dg113*(dg122 + dg212) + dg111*dg223) + 
           dg122*dg311 + dg112*(8.*dg123 + 2.*dg213 + 4.*dg312) + 
           dg111*dg322)*ginv22 - 
        (dg133*(4.*dg112 + 2.*dg211) + 
           dg113*(8.*dg123 + 4.*(dg213 + dg312)) + 4.*dg112*dg313 + 
           2.*(dg123*dg311 + dg111*(dg233 + dg323)))*ginv23 - 
        (dg133*(2.*dg113 + dg311) + 4.*dg113*dg313 + dg111*dg333)*ginv33) + 
     ginv23*(ddg1123 + 2.*ddg1213 + ddg1312 - 
        (2.*(dg133*(dg123 + dg213) + dg113*dg233) + dg133*dg312 + 
           4.*(dg123*dg313 + dg113*dg323) + dg112*dg333)*ginv33) + 
     ginv22*(ddg1122 + 2.*ddg1212 - 
        (4.*(dg122*dg213 + dg113*dg222) + 
           6.*(dg123*(dg122 + dg212) + dg112*dg223) + 
           3.*(dg122*dg312 + dg112*dg322))*ginv23 - 
        ginv33*(dg122*dg313 + dg113*dg322 + 
           2.*(dg113*dg223 + dg123*(dg213 + dg312) + dg112*dg323 + 
              pow2(dg123)))) - 
     2.*(dg113*dg311 + dg111*(dg133 + dg313) + pow2(dg113))*pow2(ginv13) - 
     (4.*(dg122*dg212 + dg112*dg222) + 2.*pow2(dg122))*pow2(ginv22) - 
     (4.*(dg123*dg213 + dg113*dg223) + 
        2.*(dg133*(dg122 + dg212) + dg123*dg312 + dg122*dg313 + 
           dg113*dg322 + dg112*(dg233 + dg323) + pow2(dg123)))*pow2(ginv23) \
- (dg133*dg313 + dg113*dg333)*pow2(ginv33)) + 
  ginv22*(ddg1323*ginv33 + ginv23*
      (2.*ddg1223 + ddg1322 - (2.*(dg133*dg223 + dg123*dg233) + 
           dg133*dg322 + 6.*dg123*dg323 + dg122*dg333)*ginv33) - 
     (2.*(dg133*dg222 + dg122*dg233) + dg123*(6.*dg223 + 3.*dg322) + 
        3.*dg122*dg323)*pow2(ginv23) - 
     (dg133*dg323 + dg123*dg333)*pow2(ginv33))
;

dGfromgdu13
=
-((dg133*dg222 + 2.*dg123*(dg223 + dg322) + dg122*(dg233 + 2.*dg323))*
     Power(ginv23,3)) - 2.*(dg133*dg333*Power(ginv33,3) + 
     Power(ginv13,3)*(dg113*dg311 + dg111*(dg133 + dg313) + pow2(dg113)) + 
     (dg111*(dg112*ginv23 + dg113*ginv33) + ginv13*pow2(dg111))*pow2(ginv11)\
) - ((dg122*(4.*dg112 + dg211) + 2.*dg112*dg212 + dg111*dg222)*ginv23 + 
     (2.*(dg113*dg122 + dg112*dg123) + dg123*dg211 + dg113*dg212 + 
        dg112*dg213 + dg111*dg223)*ginv33 + 
     2.*ginv13*(dg112*dg211 + dg111*(dg122 + dg212) + pow2(dg112)))*
   pow2(ginv12) + (ddg1113 + ddg1311 - 
     (dg123*(2.*dg112 + dg211) + dg113*dg212 + dg111*dg223 + 
        dg112*(dg213 + 2.*dg312))*ginv22 - 
     (dg133*dg211 + 2.*(dg113*dg213 + dg123*dg311) + 
        4.*(dg113*(dg123 + dg312) + dg112*(dg133 + dg313)) + 
        dg111*(dg233 + 2.*dg323))*ginv23 - 
     (dg133*(6.*dg113 + 2.*dg311) + 6.*dg113*dg313 + 2.*dg111*dg333)*ginv33\
)*pow2(ginv13) - (2.*dg122*dg222*ginv23 + 
     (dg123*dg222 + dg122*dg223)*ginv33)*pow2(ginv22) + 
  (ddg1223 + ddg1322 - (3.*(dg133*dg223 + dg123*dg233) + 6.*dg123*dg323 + 
        2.*(dg133*dg322 + dg122*dg333))*ginv33)*pow2(ginv23) + 
  ddg1333*pow2(ginv33) + ginv11*
   (ddg1113*ginv33 - ginv22*(2.*dg112*(dg122 + dg212)*ginv23 + 
        (dg113*dg212 + dg112*(2.*dg123 + dg213))*ginv33) + 
     ginv23*(ddg1112 - (dg113*(4.*dg123 + 2.*dg213) + 
           2.*(dg113*dg312 + dg112*(dg133 + dg313)))*ginv33) - 
     ginv12*(dg111*(6.*dg112 + 2.*dg211)*ginv13 + 
        (dg113*(4.*dg112 + dg211) + dg111*(2.*dg123 + dg213))*ginv33 + 
        ginv23*(dg112*dg211 + dg111*(2.*dg122 + dg212) + 4.*pow2(dg112))) + 
     ginv13*(ddg1111 - (dg113*(8.*dg112 + dg211) + 2.*dg112*dg311 + 
           dg111*(dg213 + 2.*(dg123 + dg312)))*ginv23 - 
        ginv22*(dg112*dg211 + dg111*dg212 + 2.*pow2(dg112)) - 
        ginv33*(2.*(dg113*dg311 + dg111*(dg133 + dg313)) + 6.*pow2(dg113))) \
- dg111*(6.*dg113 + 2.*dg311)*pow2(ginv13) - 
     (dg113*dg212 + dg112*dg213 + 
        2.*(dg113*dg122 + dg112*(dg123 + dg312)))*pow2(ginv23) - 
     2.*dg113*(dg133 + dg313)*pow2(ginv33)) + 
  ginv12*((ddg1123 + ddg1213)*ginv33 + 
     ginv13*(ddg1112 + ddg1211 - 
        (dg122*(2.*dg112 + dg211) + 4.*dg112*dg212 + dg111*dg222)*ginv22 - 
        (dg123*(8.*dg112 + 2.*dg211) + 
           4.*(dg113*(dg122 + dg212) + dg112*(dg213 + dg312)) + 
           2.*(dg122*dg311 + dg111*(dg223 + dg322)))*ginv23 - 
        (dg133*(2.*dg112 + dg211) + 
           dg113*(8.*dg123 + 4.*dg213 + 2.*dg312) + 
           2.*(dg123*dg311 + dg112*dg313) + dg111*(dg233 + 2.*dg323))*
         ginv33) - ginv22*((dg122*dg213 + dg113*dg222 + 
           2.*(dg123*(dg122 + dg212) + dg112*dg223))*ginv33 + 
        ginv23*(3.*(dg122*dg212 + dg112*dg222) + 2.*pow2(dg122))) + 
     ginv23*(ddg1122 + ddg1212 - 
        ginv33*(dg133*(2.*dg122 + dg212) + 
           2.*(dg123*dg312 + dg122*dg313 + dg113*dg322) + 
           dg112*(dg233 + 2.*dg323) + 
           4.*(dg123*dg213 + dg113*dg223 + pow2(dg123)))) - 
     (4.*(dg112*dg113 + dg111*dg123) + 
        2.*(dg113*dg211 + dg112*dg311 + dg111*(dg213 + dg312)))*
      pow2(ginv13) - (dg123*(4.*dg122 + 2.*dg212) + 
        2.*(dg113*dg222 + dg122*(dg213 + dg312) + dg112*(dg223 + dg322)))*
      pow2(ginv23) - (dg133*(2.*dg123 + dg213) + 2.*dg123*dg313 + 
        dg113*(dg233 + 2.*dg323))*pow2(ginv33)) + 
  ginv22*(ddg1223*ginv33 + ginv23*
      (ddg1222 - (dg133*dg222 + dg123*(6.*dg223 + 2.*dg322) + 
           dg122*(dg233 + 2.*dg323))*ginv33) - 
     (3.*(dg123*dg222 + dg122*dg223) + 2.*dg122*dg322)*pow2(ginv23) - 
     (dg133*dg223 + dg123*(dg233 + 2.*dg323))*pow2(ginv33)) + 
  ginv23*((ddg1233 + 2.*ddg1323)*ginv33 - 
     (dg133*(2.*dg233 + 4.*dg323) + 4.*dg123*dg333)*pow2(ginv33)) + 
  ginv13*((ddg1133 + 2.*ddg1313)*ginv33 + 
     ginv23*(ddg1123 + ddg1213 + 2.*ddg1312 - 
        (dg133*(6.*dg123 + 3.*dg213 + 4.*dg312) + 6.*dg123*dg313 + 
           dg113*(3.*dg233 + 6.*dg323) + 4.*dg112*dg333)*ginv33) + 
     ginv22*(ddg1212 - (dg123*(2.*dg122 + 4.*dg212) + dg113*dg222 + 
           dg122*(dg213 + 2.*dg312) + dg112*(4.*dg223 + 2.*dg322))*ginv23 - 
        ginv33*(dg133*dg212 + dg112*(dg233 + 2.*dg323) + 
           2.*(dg113*dg223 + dg123*(dg213 + dg312) + pow2(dg123)))) - 
     (dg122*dg212 + dg112*dg222)*pow2(ginv22) - 
     (4.*(dg123*dg312 + dg112*dg323) + 
        2.*(dg133*(dg122 + dg212) + dg123*dg213 + dg112*dg233 + 
           dg122*dg313 + dg113*(dg223 + dg322) + pow2(dg123)))*pow2(ginv23) \
- (4.*(dg133*dg313 + dg113*dg333) + 2.*pow2(dg133))*pow2(ginv33))
;

dGfromgdu21
=
-((dg233*dg311 + 2.*(dg113*dg233 + dg213*(dg133 + dg313)) + dg211*dg333)*
     Power(ginv13,3)) - 2.*(dg111*dg211*Power(ginv11,3) + 
     Power(ginv12,3)*(dg122*dg212 + (dg112 + dg211)*dg222 + pow2(dg212))) + 
  (ddg1211 - (4.*(dg113*dg211 + dg111*dg213) + 2.*dg211*dg311)*ginv13 - 
     2.*(dg112 + dg211)*dg212*ginv22 - 
     (2.*(dg113*dg212 + (dg112 + dg211)*dg213) + dg212*dg311 + 
        dg211*dg312)*ginv23 - 
     (dg213*(2.*dg113 + dg311) + dg211*dg313)*ginv33 - 
     ginv12*(4.*(dg112*dg211 + dg111*dg212) + 2.*pow2(dg211)))*pow2(ginv11) \
+ (ddg1222 + ddg2212 - (4.*(dg212*(dg123 + dg213) + 
           (dg112 + dg211)*dg223) + dg222*dg311 + 
        2.*(dg122*dg213 + dg113*dg222 + dg212*dg312) + dg211*dg322)*ginv13 \
- (2.*dg122 + 6.*dg212)*dg222*ginv22 - 
     ((2.*dg122 + 4.*dg212)*dg223 + 
        dg222*(4.*dg213 + 2.*(dg123 + dg312)) + 2.*dg212*dg322)*ginv23 - 
     (dg223*(2.*(dg123 + dg213) + dg312) + dg222*dg313 + dg213*dg322 + 
        dg212*dg323)*ginv33)*pow2(ginv12) + 
  (ddg1233 + ddg2313 - (2.*((dg123 + dg213)*dg223 + dg212*dg233) + 
        dg223*dg312 + dg212*dg323)*ginv22 - 
     (dg233*(4.*dg213 + 2.*dg312) + 
        2.*(dg123*dg233 + dg223*(dg133 + dg313) + dg213*dg323 + 
           dg212*dg333))*ginv23 - 
     (dg233*(2.*dg133 + 3.*dg313) + 3.*dg213*dg333)*ginv33)*pow2(ginv13) + 
  ginv11*(ddg2313*ginv33 + ginv22*
      (ddg2212 - (dg222*(2.*dg213 + dg312) + dg212*(4.*dg223 + dg322))*
         ginv23 - (dg223*(2.*dg213 + dg312) + dg212*dg323)*ginv33) + 
     ginv23*(ddg2213 + ddg2312 - 
        (dg233*(2.*dg213 + dg312) + 2.*(dg223*dg313 + dg213*dg323) + 
           dg212*dg333)*ginv33) + 
     ginv13*(2.*ddg1213 + ddg2311 - 
        (2.*(dg112 + dg211)*dg223 + 
           dg212*(4.*dg213 + 2.*(dg123 + dg312)))*ginv22 - 
        (2.*(dg133*dg213 + dg113*dg233) + dg233*dg311 + 6.*dg213*dg313 + 
           dg211*dg333)*ginv33 - 
        ginv23*(2.*(dg133*dg212 + dg123*dg213 + dg113*dg223 + 
              (dg112 + dg211)*dg233) + dg223*dg311 + dg211*dg323 + 
           4.*(dg213*dg312 + dg212*dg313 + pow2(dg213)))) + 
     ginv12*(2.*ddg1212 + ddg2211 - 
        (6.*(dg113*dg212 + dg112*dg213) + 4.*dg111*dg223 + 
           3.*dg212*dg311 + dg211*(4.*dg123 + 6.*dg213 + 3.*dg312))*ginv13 \
- (2.*(dg123*dg212 + dg122*dg213 + dg113*dg222 + 
              (dg112 + dg211)*dg223) + dg222*dg311 + 
           dg212*(8.*dg213 + 4.*dg312) + dg211*dg322)*ginv23 - 
        ginv22*(2.*(dg122*dg212 + (dg112 + dg211)*dg222) + 
           6.*pow2(dg212)) - ginv33*
         (dg223*dg311 + dg211*dg323 + 
           2.*(dg113*dg223 + dg213*(dg123 + dg312) + dg212*dg313 + 
              pow2(dg213)))) - 
     (6.*dg112*dg212 + dg211*(2.*dg122 + 6.*dg212) + 2.*dg111*dg222)*
      pow2(ginv12) - (2.*(dg133*dg211 + dg111*dg233) + 
        dg213*(6.*dg113 + 3.*dg311) + 3.*dg211*dg313)*pow2(ginv13) - 
     2.*dg212*dg222*pow2(ginv22) - 
     (2.*(dg213*dg223 + dg212*dg233) + dg223*dg312 + dg222*dg313 + 
        dg213*dg322 + dg212*dg323)*pow2(ginv23) - 
     (dg233*dg313 + dg213*dg333)*pow2(ginv33)) + 
  ginv12*(ddg2323*ginv33 + ginv13*
      (2.*ddg1223 + ddg2213 + ddg2312 - 
        (2.*((dg123 + dg213)*dg222 + dg122*dg223) + dg222*dg312 + 
           dg212*(8.*dg223 + dg322))*ginv22 - 
        (dg223*(8.*dg213 + 4.*(dg123 + dg312)) + 
           2.*(dg122*dg233 + dg222*(dg133 + dg313) + dg213*dg322) + 
           4.*dg212*(dg233 + dg323))*ginv23 - 
        (2.*(dg133*dg223 + (dg123 + dg213)*dg233) + dg233*dg312 + 
           4.*(dg223*dg313 + dg213*dg323) + dg212*dg333)*ginv33) + 
     ginv23*(ddg2223 + ddg2322 - 
        (dg233*(2.*dg223 + dg322) + 4.*dg223*dg323 + dg222*dg333)*ginv33) + 
     ginv22*(ddg2222 - dg222*(6.*dg223 + 2.*dg322)*ginv23 - 
        ginv33*(dg223*dg322 + dg222*dg323 + 2.*pow2(dg223))) - 
     (4.*(dg123*dg213 + dg113*dg223) + 
        2.*((dg112 + dg211)*dg233 + dg223*dg311 + dg213*dg312 + 
           dg212*(dg133 + dg313) + dg211*dg323 + pow2(dg213)))*pow2(ginv13) \
- 2.*(pow2(dg222)*pow2(ginv22) + 
        (dg223*dg322 + dg222*(dg233 + dg323) + pow2(dg223))*pow2(ginv23)) - 
     (dg233*dg323 + dg223*dg333)*pow2(ginv33)) + 
  ginv13*(ddg2333*ginv33 + ginv22*
      (ddg2223 - 2.*dg223*(dg233 + dg323)*ginv33 - 
        ginv23*(dg223*dg322 + dg222*(2.*dg233 + dg323) + 4.*pow2(dg223))) + 
     ginv23*(ddg2233 + ddg2323 - 
        ginv33*(3.*(dg233*dg323 + dg223*dg333) + 2.*pow2(dg233))) - 
     (dg233*(4.*dg223 + dg322) + 2.*dg223*dg323 + dg222*dg333)*
      pow2(ginv23) - 2.*(dg222*dg223*pow2(ginv22) + dg233*dg333*pow2(ginv33))\
)
;

dGfromgdu22
=
-((2.*dg112*dg212 + dg211*(dg122 + 4.*dg212) + dg111*dg222)*
     Power(ginv12,3)) - (dg233*(4.*dg223 + dg322) + 2.*dg223*dg323 + 
     dg222*dg333)*Power(ginv23,3) - 2.*Power(ginv22,3)*pow2(dg222) - 
  (2.*dg111*dg211*ginv12 + (dg112*dg211 + dg111*dg212)*ginv22 + 
     (dg113*dg211 + dg111*dg213)*ginv23)*pow2(ginv11) + 
  (ddg1212 + ddg2211 - (2.*(dg123*dg211 + dg112*dg213 + dg111*dg223 + 
           dg212*(dg113 + dg311)) + dg211*(4.*dg213 + 2.*dg312))*ginv13 - 
     (2.*(dg123*dg212 + dg122*dg213 + dg113*dg222 + dg112*dg223) + 
        dg222*dg311 + dg212*(8.*dg213 + 2.*dg312) + 
        dg211*(4.*dg223 + dg322))*ginv23 - 
     ginv22*(4.*dg211*dg222 + 3.*(dg122*dg212 + dg112*dg222) + 
        6.*pow2(dg212)) - ginv33*
      (dg223*(dg113 + dg311) + dg213*(dg123 + dg312) + dg212*dg313 + 
        dg211*dg323 + 2.*pow2(dg213)))*pow2(ginv12) - 
  ((dg112*dg233 + dg223*(dg113 + dg311) + dg213*(dg123 + dg312) + 
        dg212*(dg133 + dg313) + dg211*dg323)*ginv22 + 
     (dg233*dg311 + 2.*(dg113*dg233 + dg213*(dg133 + dg313)) + 
        dg211*dg333)*ginv23)*pow2(ginv13) + 
  (ddg2222 - dg222*(8.*dg223 + 2.*dg322)*ginv23 - 
     ginv33*(dg223*dg322 + dg222*dg323 + 2.*pow2(dg223)))*pow2(ginv22) + 
  (ddg2233 + ddg2323 - ginv33*
      (3.*(dg233*dg323 + dg223*dg333) + 2.*pow2(dg233)))*pow2(ginv23) + 
  ginv13*(ginv22*(ddg1223 + ddg2312 - 
        (dg122*dg233 + dg222*(dg133 + dg313) + dg213*dg322 + 
           4.*(dg223*(dg123 + dg213 + dg312) + dg212*(dg233 + dg323)))*
         ginv23 - (dg233*(dg123 + dg312) + dg223*(dg133 + 2.*dg313) + 
           2.*dg213*dg323 + dg212*dg333)*ginv33) + 
     ginv23*(ddg1233 + ddg2313 - 
        (dg233*(2.*dg133 + 3.*dg313) + 3.*dg213*dg333)*ginv33) - 
     ((dg122 + 4.*dg212)*dg223 + dg222*(dg123 + dg312) + dg212*dg322)*
      pow2(ginv22) - (dg233*(4.*dg213 + 2.*dg312) + 
        2.*(dg123*dg233 + dg223*(dg133 + dg313) + dg213*dg323 + 
           dg212*dg333))*pow2(ginv23)) + 
  ginv11*(-(ginv13*((2.*(dg113*dg212 + dg112*dg213) + dg111*dg223 + 
             dg212*dg311 + dg211*(dg123 + dg312))*ginv22 + 
          (dg111*dg233 + dg213*(4.*dg113 + dg311) + dg211*(dg133 + dg313))*
           ginv23)) + ginv12*(ddg1211 - 
        (3.*(dg113*dg211 + dg111*dg213) + 2.*dg211*dg311)*ginv13 - 
        (6.*dg112*dg212 + dg211*(dg122 + 4.*dg212) + dg111*dg222)*ginv22 - 
        (4.*(dg113*dg212 + dg112*dg213) + dg111*dg223 + dg212*dg311 + 
           dg211*(dg123 + 4.*dg213 + dg312))*ginv23 - 
        (dg213*(2.*dg113 + dg311) + dg211*dg313)*ginv33) + 
     ginv22*(ddg1212 - (dg122*dg213 + dg113*dg222 + 2.*dg112*dg223 + 
           dg212*(4.*dg213 + 2.*(dg123 + dg312)))*ginv23 - 
        (dg113*dg223 + dg213*(dg123 + dg312) + dg212*dg313)*ginv33) + 
     ginv23*(ddg1213 - (dg113*dg233 + dg213*(dg133 + 2.*dg313))*ginv33) - 
     (3.*(dg112*dg211 + dg111*dg212) + 2.*pow2(dg211))*pow2(ginv12) - 
     (dg122*dg212 + dg112*dg222 + 2.*pow2(dg212))*pow2(ginv22) - 
     (dg113*dg223 + dg112*dg233 + dg213*(dg123 + dg312) + 
        dg212*(dg133 + dg313) + 2.*pow2(dg213))*pow2(ginv23)) + 
  ginv23*(ddg2333*ginv33 - 2.*dg233*dg333*pow2(ginv33)) + 
  ginv12*(ddg2313*ginv33 + ginv22*
      (ddg1222 + 2.*ddg2212 - 
        ((3.*dg122 + 12.*dg212)*dg223 + 
           dg222*(8.*dg213 + 3.*(dg123 + dg312)) + 3.*dg212*dg322)*ginv23 \
- (dg223*(4.*dg213 + 2.*(dg123 + dg312)) + dg222*dg313 + dg213*dg322 + 
           2.*dg212*dg323)*ginv33) + 
     ginv23*(ddg1223 + 2.*ddg2213 + ddg2312 - 
        (dg233*(dg123 + 4.*dg213 + dg312) + dg223*(dg133 + 4.*dg313) + 
           4.*dg213*dg323 + dg212*dg333)*ginv33) + 
     ginv13*(ddg1213 + ddg2311 - 
        (dg122*dg213 + dg222*(dg113 + dg311) + 
           4.*((dg112 + dg211)*dg223 + dg212*(dg123 + dg213 + dg312)) + 
           dg211*dg322)*ginv22 - 
        (dg233*(dg113 + dg311) + dg213*(dg133 + 4.*dg313) + dg211*dg333)*
         ginv33 - ginv23*(2.*(dg133*dg212 + dg112*dg233 + dg223*dg311 + 
              dg211*dg323) + 4.*
            (dg113*dg223 + dg211*dg233 + dg213*(dg123 + dg312) + 
              dg212*dg313 + pow2(dg213)))) - 
     (dg111*dg233 + 2.*dg213*(dg113 + dg311) + dg211*(dg133 + 2.*dg313))*
      pow2(ginv13) - (2.*dg122 + 8.*dg212)*dg222*pow2(ginv22) - 
     ((dg122 + 4.*dg212)*dg233 + dg223*(8.*dg213 + 2.*(dg123 + dg312)) + 
        dg222*(dg133 + 2.*dg313) + 2.*(dg213*dg322 + dg212*dg323))*
      pow2(ginv23) - (dg233*dg313 + dg213*dg333)*pow2(ginv33)) + 
  ginv22*(ddg2323*ginv33 + ginv23*
      (2.*ddg2223 + ddg2322 - (dg233*(4.*dg223 + dg322) + 6.*dg223*dg323 + 
           dg222*dg333)*ginv33) - 
     (3.*dg223*dg322 + dg222*(4.*dg233 + 3.*dg323) + 6.*pow2(dg223))*
      pow2(ginv23) - (dg233*dg323 + dg223*dg333)*pow2(ginv33))
;

dGfromgdu23
=
-((dg111*dg233 + 2.*dg213*(dg113 + dg311) + dg211*(dg133 + 2.*dg313))*
     Power(ginv13,3)) - (2.*dg111*dg211*ginv13 + 
     (dg112*dg211 + dg111*dg212)*ginv23 + 
     (dg113*dg211 + dg111*dg213)*ginv33)*pow2(ginv11) - 
  ((2.*dg112*dg212 + dg211*(dg122 + 4.*dg212) + dg111*dg222)*ginv13 + 
     (dg122*dg213 + dg212*(dg123 + 2.*dg213) + dg113*dg222 + 
        (dg112 + 2.*dg211)*dg223)*ginv33 + 
     2.*ginv23*(dg122*dg212 + (dg112 + dg211)*dg222 + pow2(dg212)))*
   pow2(ginv12) + (ddg1213 + ddg2311 - 
     ((dg112 + 2.*dg211)*dg223 + dg212*(dg123 + 2.*(dg213 + dg312)))*
      ginv22 - (3.*(dg133*dg213 + dg113*dg233) + 6.*dg213*dg313 + 
        2.*(dg233*dg311 + dg211*dg333))*ginv33 - 
     ginv23*(4.*(dg213*dg312 + dg212*dg313) + 
        2.*(dg133*dg212 + dg123*dg213 + (dg112 + dg211)*dg233 + 
           dg223*(dg113 + dg311) + dg211*dg323 + pow2(dg213))))*pow2(ginv13) \
- 2.*(dg233*dg333*Power(ginv33,3) + 
     Power(ginv23,3)*(dg223*dg322 + dg222*(dg233 + dg323) + pow2(dg223)) + 
     (dg222*dg223*ginv33 + ginv23*pow2(dg222))*pow2(ginv22)) + 
  (ddg2223 + ddg2322 - (dg233*(6.*dg223 + 2.*dg322) + 6.*dg223*dg323 + 
        2.*dg222*dg333)*ginv33)*pow2(ginv23) + ddg2333*pow2(ginv33) + 
  ginv11*(ddg1213*ginv33 + ginv13*
      (ddg1211 - 2.*(dg112 + dg211)*dg212*ginv22 - 
        (4.*(dg113*dg212 + dg112*dg213) + dg111*dg223 + 2.*dg212*dg311 + 
           dg211*(dg123 + 2.*(dg213 + dg312)))*ginv23 - 
        (dg111*dg233 + dg213*(6.*dg113 + 2.*dg311) + 
           dg211*(dg133 + 2.*dg313))*ginv33) - 
     ginv12*((4.*dg112*dg212 + dg211*(dg122 + 2.*dg212) + dg111*dg222)*
         ginv23 + (dg211*(dg123 + 2.*dg213) + 
           2.*(dg113*dg212 + dg112*dg213) + dg111*dg223)*ginv33 + 
        ginv13*(3.*(dg112*dg211 + dg111*dg212) + 2.*pow2(dg211))) - 
     ginv22*((dg212*(dg123 + 2.*dg213) + dg112*dg223)*ginv33 + 
        ginv23*(dg122*dg212 + dg112*dg222 + 2.*pow2(dg212))) + 
     ginv23*(ddg1212 - ginv33*
         (dg112*dg233 + dg212*(dg133 + 2.*dg313) + 
           2.*(dg113*dg223 + dg213*(dg123 + dg312) + pow2(dg213)))) - 
     (3.*(dg113*dg211 + dg111*dg213) + 2.*dg211*dg311)*pow2(ginv13) - 
     (dg122*dg213 + dg113*dg222 + dg112*dg223 + 
        dg212*(dg123 + 2.*(dg213 + dg312)))*pow2(ginv23) - 
     (dg113*dg233 + dg213*(dg133 + 2.*dg313))*pow2(ginv33)) + 
  ginv22*(ddg2223*ginv33 + ginv23*
      (ddg2222 - ginv33*(2.*(dg223*dg322 + dg222*(dg233 + dg323)) + 
           6.*pow2(dg223))) - dg222*(6.*dg223 + 2.*dg322)*pow2(ginv23) - 
     2.*dg223*(dg233 + dg323)*pow2(ginv33)) + 
  ginv12*((ddg1223 + ddg2213)*ginv33 - 
     ginv22*((2.*dg122 + 6.*dg212)*dg222*ginv23 + 
        ((dg123 + 2.*dg213)*dg222 + (dg122 + 4.*dg212)*dg223)*ginv33) + 
     ginv23*(ddg1222 + ddg2212 - 
        ((dg122 + 2.*dg212)*dg233 + 
           dg223*(4.*dg123 + 8.*dg213 + 2.*dg312) + 
           dg222*(dg133 + 2.*dg313) + 2.*(dg213*dg322 + dg212*dg323))*
         ginv33) + ginv13*(ddg1212 + ddg2211 - 
        (4.*(dg112 + dg211)*dg223 + 
           dg212*(8.*dg213 + 4.*(dg123 + dg312)) + 
           2.*(dg122*dg213 + dg222*(dg113 + dg311) + dg211*dg322))*ginv23 \
- ginv22*(dg122*dg212 + (dg112 + 2.*dg211)*dg222 + 4.*pow2(dg212)) - 
        ginv33*((dg112 + 2.*dg211)*dg233 + dg212*(dg133 + 2.*dg313) + 
           2.*(dg223*dg311 + dg213*dg312 + dg211*dg323) + 
           4.*(dg123*dg213 + dg113*dg223 + pow2(dg213)))) - 
     (2.*(dg123*dg211 + dg112*dg213 + dg111*dg223 + 
           dg212*(dg113 + dg311)) + dg211*(4.*dg213 + 2.*dg312))*
      pow2(ginv13) - ((2.*dg122 + 4.*dg212)*dg223 + 
        dg222*(4.*dg213 + 2.*(dg123 + dg312)) + 2.*dg212*dg322)*
      pow2(ginv23) - ((dg123 + 2.*dg213)*dg233 + 
        dg223*(dg133 + 2.*dg313) + 2.*dg213*dg323)*pow2(ginv33)) + 
  ginv13*((ddg1233 + 2.*ddg2313)*ginv33 + 
     ginv22*(ddg2212 - ((dg122 + 8.*dg212)*dg223 + 
           dg222*(dg123 + 2.*(dg213 + dg312)) + 2.*dg212*dg322)*ginv23 - 
        (dg223*(4.*dg213 + 2.*(dg123 + dg312)) + 2.*dg212*(dg233 + dg323))*
         ginv33) + ginv23*(ddg1223 + ddg2213 + 2.*ddg2312 - 
        (3.*(dg133*dg223 + dg123*dg233) + dg233*(6.*dg213 + 4.*dg312) + 
           6.*(dg223*dg313 + dg213*dg323) + 4.*dg212*dg333)*ginv33) - 
     2.*dg212*dg222*pow2(ginv22) - 
     ((dg122 + 4.*dg212)*dg233 + dg223*(2.*dg123 + 4.*(dg213 + dg312)) + 
        dg222*(dg133 + 2.*dg313) + 2.*dg213*dg322 + 4.*dg212*dg323)*
      pow2(ginv23) - (dg233*(2.*dg133 + 4.*dg313) + 4.*dg213*dg333)*
      pow2(ginv33)) + ginv23*((ddg2233 + 2.*ddg2323)*ginv33 - 
     (4.*(dg233*dg323 + dg223*dg333) + 2.*pow2(dg233))*pow2(ginv33))
;

dGfromgdu31
=
-((dg222*dg311 + dg211*dg322 + 2.*((dg122 + dg212)*dg312 + dg112*dg322))*
     Power(ginv12,3)) - 2.*(dg111*dg311*Power(ginv11,3) + 
     Power(ginv13,3)*(dg133*dg313 + (dg113 + dg311)*dg333 + pow2(dg313))) + 
  (ddg1311 - ((4.*dg112 + 2.*dg211)*dg311 + 4.*dg111*dg312)*ginv12 - 
     (dg212*dg311 + (2.*dg112 + dg211)*dg312)*ginv22 - 
     (dg311*(dg213 + 2.*dg312) + dg211*dg313 + 
        2.*(dg113*dg312 + dg112*dg313))*ginv23 - 
     2.*(dg113 + dg311)*dg313*ginv33 - 
     ginv13*(4.*(dg113*dg311 + dg111*dg313) + 2.*pow2(dg311)))*pow2(ginv11) \
+ (ddg1322 + ddg2312 - (2.*dg122*dg322 + 3.*(dg222*dg312 + dg212*dg322))*
      ginv22 - ((2.*dg213 + 4.*dg312)*dg322 + 
        2.*(dg223*dg312 + dg222*dg313 + dg123*dg322 + 
           (dg122 + dg212)*dg323))*ginv23 - 
     (dg313*(dg223 + 2.*dg322) + (dg213 + 2.*(dg123 + dg312))*dg323)*
      ginv33 - ginv13*(4.*(dg123*dg312 + dg112*dg323) + 
        2.*(dg213*dg312 + (dg122 + dg212)*dg313 + dg113*dg322 + 
           dg311*(dg223 + dg322) + dg211*dg323 + pow2(dg312))))*pow2(ginv12) \
+ (ddg1333 + ddg3313 - (dg233*dg312 + dg223*dg313 + 
        (dg213 + 2.*(dg123 + dg312))*dg323 + dg212*dg333)*ginv22 - 
     (2.*(dg233*dg313 + dg133*dg323 + (dg123 + dg213)*dg333) + 
        4.*(dg313*dg323 + dg312*dg333))*ginv23 - 
     (2.*dg133 + 6.*dg313)*dg333*ginv33)*pow2(ginv13) + 
  ginv11*(ddg3313*ginv33 + ginv22*
      (ddg2312 - (dg222*dg313 + dg213*dg322 + 
           2.*(dg312*(dg223 + dg322) + dg212*dg323))*ginv23 - 
        (dg223*dg313 + (dg213 + 2.*dg312)*dg323)*ginv33) + 
     ginv23*(ddg2313 + ddg3312 - 
        (dg313*(dg233 + 4.*dg323) + (dg213 + 2.*dg312)*dg333)*ginv33) + 
     ginv12*(2.*ddg1312 + ddg2311 - 
        (dg311*(4.*dg123 + 3.*dg213 + 6.*dg312) + 3.*dg211*dg313 + 
           6.*(dg113*dg312 + dg112*dg313) + 4.*dg111*dg323)*ginv13 - 
        (dg222*dg311 + (2.*dg122 + 6.*dg212)*dg312 + 
           (2.*dg112 + dg211)*dg322)*ginv22 - 
        (4.*dg312*dg313 + 2.*((dg123 + dg213)*dg313 + 
              (dg113 + dg311)*dg323))*ginv33 - 
        ginv23*((2.*dg123 + 4.*dg213)*dg312 + dg311*(dg223 + 2.*dg322) + 
           dg211*dg323 + 2.*(dg122*dg313 + dg113*dg322 + dg112*dg323) + 
           4.*(dg212*dg313 + pow2(dg312)))) + 
     ginv13*(2.*ddg1313 + ddg3311 - 
        ((4.*dg213 + 8.*dg312)*dg313 + dg311*(dg233 + 2.*dg323) + 
           dg211*dg333 + 2.*(dg133*dg312 + dg123*dg313 + dg113*dg323 + 
              dg112*dg333))*ginv23 - 
        ginv22*(dg223*dg311 + dg211*dg323 + 
           2.*((dg123 + dg213)*dg312 + dg212*dg313 + dg112*dg323 + 
              pow2(dg312))) - 
        ginv33*(2.*(dg133*dg313 + (dg113 + dg311)*dg333) + 6.*pow2(dg313))) \
- ((2.*dg122 + 3.*dg212)*dg311 + (6.*dg112 + 3.*dg211)*dg312 + 
        2.*dg111*dg322)*pow2(ginv12) - 
     (6.*dg113*dg313 + dg311*(2.*dg133 + 6.*dg313) + 2.*dg111*dg333)*
      pow2(ginv13) - (dg222*dg312 + dg212*dg322)*pow2(ginv22) - 
     (dg313*(dg223 + 2.*dg322) + dg213*dg323 + dg312*(dg233 + 2.*dg323) + 
        dg212*dg333)*pow2(ginv23) - 2.*dg313*dg333*pow2(ginv33)) + 
  ginv12*(ddg3323*ginv33 + ginv13*
      (2.*ddg1323 + ddg2313 + ddg3312 - 
        (dg222*dg313 + (2.*dg123 + dg213)*dg322 + 
           dg312*(4.*dg223 + 2.*dg322) + (2.*dg122 + 4.*dg212)*dg323)*
         ginv22 - ((4.*dg213 + 8.*dg312)*dg323 + 
           4.*(dg313*(dg223 + dg322) + dg123*dg323) + 
           2.*(dg233*dg312 + dg133*dg322 + (dg122 + dg212)*dg333))*ginv23 \
- (dg313*(dg233 + 8.*dg323) + (dg213 + 2.*dg312)*dg333 + 
           2.*(dg133*dg323 + dg123*dg333))*ginv33) + 
     ginv22*(ddg2322 - 2.*(dg223 + dg322)*dg323*ginv33 - 
        ginv23*(3.*(dg223*dg322 + dg222*dg323) + 2.*pow2(dg322))) + 
     ginv23*(ddg2323 + ddg3322 - 
        ginv33*(dg233*dg323 + (dg223 + 2.*dg322)*dg333 + 4.*pow2(dg323))) - 
     (dg311*(dg233 + 4.*dg323) + 
        4.*((dg123 + dg312)*dg313 + dg113*dg323) + dg211*dg333 + 
        2.*(dg133*dg312 + dg213*dg313 + dg112*dg333))*pow2(ginv13) - 
     (2.*dg223*dg323 + dg322*(dg233 + 4.*dg323) + dg222*dg333)*
      pow2(ginv23) - 2.*(dg222*dg322*pow2(ginv22) + 
        dg323*dg333*pow2(ginv33))) + 
  ginv13*(ddg3333*ginv33 + ginv23*
      (ddg2333 + ddg3323 - (2.*dg233 + 6.*dg323)*dg333*ginv33) + 
     ginv22*(ddg2323 - (4.*dg223*dg323 + dg322*(dg233 + 2.*dg323) + 
           dg222*dg333)*ginv23 - 
        ginv33*(dg233*dg323 + dg223*dg333 + 2.*pow2(dg323))) - 
     (dg223*dg322 + dg222*dg323)*pow2(ginv22) - 
     2.*((dg233*dg323 + (dg223 + dg322)*dg333 + pow2(dg323))*pow2(ginv23) + 
        pow2(dg333)*pow2(ginv33)))
;

dGfromgdu32
=
-(((dg122 + 2.*dg212)*dg311 + 2.*(dg112 + dg211)*dg312 + dg111*dg322)*
     Power(ginv12,3)) - 2.*(dg222*dg322*Power(ginv22,3) + 
     Power(ginv23,3)*(dg233*dg323 + (dg223 + dg322)*dg333 + pow2(dg323))) - 
  (2.*dg111*dg311*ginv12 + (dg112*dg311 + dg111*dg312)*ginv22 + 
     (dg113*dg311 + dg111*dg313)*ginv23)*pow2(ginv11) + 
  (ddg1312 + ddg2311 - (4.*dg311*dg312 + 
        2.*((dg123 + dg213)*dg311 + dg113*dg312 + 
           (dg112 + dg211)*dg313 + dg111*dg323))*ginv13 - 
     ((3.*dg122 + 6.*dg212)*dg312 + 3.*dg112*dg322 + 
        2.*(dg222*dg311 + dg211*dg322))*ginv22 - 
     ((dg123 + 2.*(dg213 + dg312))*dg313 + (dg113 + 2.*dg311)*dg323)*
      ginv33 - ginv23*(4.*(dg213*dg312 + dg212*dg313) + 
        2.*(dg123*dg312 + dg122*dg313 + dg113*dg322 + 
           dg311*(dg223 + dg322) + (dg112 + dg211)*dg323 + pow2(dg312))))*
   pow2(ginv12) - ((dg123*dg313 + dg312*(dg133 + 2.*dg313) + 
        (dg113 + 2.*dg311)*dg323 + dg112*dg333)*ginv22 + 
     2.*ginv23*(dg133*dg313 + (dg113 + dg311)*dg333 + pow2(dg313)))*
   pow2(ginv13) + (ddg2322 - 2.*(dg223 + dg322)*dg323*ginv33 - 
     ginv23*(4.*(dg223*dg322 + dg222*dg323) + 2.*pow2(dg322)))*pow2(ginv22) \
+ (ddg2333 + ddg3323 - (2.*dg233 + 6.*dg323)*dg333*ginv33)*pow2(ginv23) + 
  ginv11*(-(ginv13*((dg311*(dg123 + 2.*dg312) + 
             2.*(dg113*dg312 + dg112*dg313) + dg111*dg323)*ginv22 + 
          (4.*dg113*dg313 + dg311*(dg133 + 2.*dg313) + dg111*dg333)*ginv23)\
) + ginv12*(ddg1311 - ((dg122 + 2.*dg212)*dg311 + 
           (6.*dg112 + 2.*dg211)*dg312 + dg111*dg322)*ginv22 - 
        (dg311*(dg123 + 2.*(dg213 + dg312)) + 2.*dg211*dg313 + 
           4.*(dg113*dg312 + dg112*dg313) + dg111*dg323)*ginv23 - 
        2.*(dg113 + dg311)*dg313*ginv33 - 
        ginv13*(3.*(dg113*dg311 + dg111*dg313) + 2.*pow2(dg311))) + 
     ginv22*(ddg1312 - ((dg123 + 2.*dg312)*dg313 + dg113*dg323)*ginv33 - 
        ginv23*(dg122*dg313 + dg113*dg322 + 
           2.*((dg123 + dg213)*dg312 + dg212*dg313 + dg112*dg323 + 
              pow2(dg312)))) + 
     ginv23*(ddg1313 - ginv33*
         (dg133*dg313 + dg113*dg333 + 2.*pow2(dg313))) - 
     ((3.*dg112 + 2.*dg211)*dg311 + 3.*dg111*dg312)*pow2(ginv12) - 
     ((dg122 + 2.*dg212)*dg312 + dg112*dg322)*pow2(ginv22) - 
     (dg133*dg312 + (dg123 + 2.*(dg213 + dg312))*dg313 + dg113*dg323 + 
        dg112*dg333)*pow2(ginv23)) + 
  ginv13*(ginv23*(ddg1333 + ddg3313 - (2.*dg133 + 6.*dg313)*dg333*ginv33) + 
     ginv22*(ddg1323 + ddg3312 - 
        (dg133*dg322 + (4.*dg123 + 2.*dg213 + 8.*dg312)*dg323 + 
           dg122*dg333 + 2.*(dg233*dg312 + dg313*(dg223 + dg322) + 
              dg212*dg333))*ginv23 - 
        ((dg133 + 4.*dg313)*dg323 + (dg123 + 2.*dg312)*dg333)*ginv33) - 
     (dg123*dg322 + dg122*dg323 + 
        2.*(dg312*(dg223 + dg322) + dg212*dg323))*pow2(ginv22) - 
     (2.*(dg233*dg313 + dg133*dg323 + (dg123 + dg213)*dg333) + 
        4.*(dg313*dg323 + dg312*dg333))*pow2(ginv23)) + 
  ginv12*(ddg3313*ginv33 + ginv22*
      (ddg1322 + 2.*ddg2312 - 
        (4.*(dg222*dg313 + dg213*dg322) + 
           3.*(dg123*dg322 + dg122*dg323) + 
           6.*(dg312*(dg223 + dg322) + dg212*dg323))*ginv23 - 
        ((2.*dg213 + 4.*dg312)*dg323 + 
           2.*(dg313*(dg223 + dg322) + dg123*dg323))*ginv33) + 
     ginv23*(ddg1323 + 2.*ddg2313 + ddg3312 - 
        (dg133*dg323 + dg313*(2.*dg233 + 8.*dg323) + 
           (dg123 + 2.*(dg213 + dg312))*dg333)*ginv33) + 
     ginv13*(ddg1313 + ddg3311 - 
        (8.*dg312*dg313 + 4.*
            ((dg123 + dg213)*dg313 + (dg113 + dg311)*dg323) + 
           2.*(dg233*dg311 + dg133*dg312 + (dg112 + dg211)*dg333))*ginv23 \
- ginv22*(dg122*dg313 + dg113*dg322 + 
           2.*(dg213*dg312 + dg212*dg313 + dg311*(dg223 + dg322) + 
              dg211*dg323) + 4.*(dg123*dg312 + dg112*dg323 + pow2(dg312))) \
- ginv33*(dg133*dg313 + (dg113 + 2.*dg311)*dg333 + 4.*pow2(dg313))) - 
     (2.*dg113*dg313 + dg311*(dg133 + 4.*dg313) + dg111*dg333)*
      pow2(ginv13) - (2.*dg122*dg322 + 4.*(dg222*dg312 + dg212*dg322))*
      pow2(ginv22) - (dg133*dg322 + 
        4.*(dg313*(dg223 + dg322) + (dg213 + dg312)*dg323) + 
        dg122*dg333 + 2.*(dg233*dg312 + dg123*dg323 + dg212*dg333))*
      pow2(ginv23) - 2.*dg313*dg333*pow2(ginv33)) + 
  ginv22*(ddg3323*ginv33 + ginv23*
      (2.*ddg2323 + ddg3322 - 
        ginv33*(2.*(dg233*dg323 + (dg223 + dg322)*dg333) + 6.*pow2(dg323))) \
- (6.*dg223*dg323 + dg322*(2.*dg233 + 6.*dg323) + 2.*dg222*dg333)*
      pow2(ginv23) - 2.*dg323*dg333*pow2(ginv33)) + 
  ginv23*(ddg3333*ginv33 - 2.*pow2(dg333)*pow2(ginv33))
;

dGfromgdu33
=
-((2.*dg113*dg313 + dg311*(dg133 + 4.*dg313) + dg111*dg333)*
     Power(ginv13,3)) - (2.*dg223*dg323 + dg322*(dg233 + 4.*dg323) + 
     dg222*dg333)*Power(ginv23,3) - 2.*Power(ginv33,3)*pow2(dg333) - 
  (2.*dg111*dg311*ginv13 + (dg112*dg311 + dg111*dg312)*ginv23 + 
     (dg113*dg311 + dg111*dg313)*ginv33)*pow2(ginv11) - 
  (((dg122 + 2.*dg212)*dg311 + 2.*(dg112 + dg211)*dg312 + dg111*dg322)*
      ginv13 + (dg222*dg311 + dg211*dg322 + 
        2.*((dg122 + dg212)*dg312 + dg112*dg322))*ginv23 + 
     (dg223*dg311 + (dg123 + dg213)*dg312 + (dg122 + dg212)*dg313 + 
        dg113*dg322 + (dg112 + dg211)*dg323)*ginv33)*pow2(ginv12) + 
  (ddg1313 + ddg3311 - ((2.*dg213 + 8.*dg312)*dg313 + 
        dg311*(dg233 + 4.*dg323) + dg211*dg333 + 
        2.*(dg133*dg312 + dg123*dg313 + dg113*dg323 + dg112*dg333))*ginv23 \
- ginv22*(dg223*dg311 + (dg123 + dg213)*dg312 + dg212*dg313 + 
        (dg112 + dg211)*dg323 + 2.*pow2(dg312)) - 
     ginv33*(4.*dg311*dg333 + 3.*(dg133*dg313 + dg113*dg333) + 
        6.*pow2(dg313)))*pow2(ginv13) - 
  (2.*dg222*dg322*ginv23 + (dg223*dg322 + dg222*dg323)*ginv33)*
   pow2(ginv22) + (ddg2323 + ddg3322 - 
     ginv33*(4.*dg322*dg333 + 3.*(dg233*dg323 + dg223*dg333) + 
        6.*pow2(dg323)))*pow2(ginv23) + ddg3333*pow2(ginv33) + 
  ginv13*((ddg1333 + 2.*ddg3313)*ginv33 + 
     ginv22*(ddg2312 - (dg222*dg313 + (dg123 + dg213)*dg322 + 
           dg122*dg323 + 4.*(dg312*(dg223 + dg322) + dg212*dg323))*ginv23 \
- (dg312*(dg233 + 4.*dg323) + 2.*(dg223*dg313 + (dg123 + dg213)*dg323) + 
           dg212*dg333)*ginv33) + 
     ginv23*(ddg1323 + ddg2313 + 2.*ddg3312 - 
        (12.*dg313*dg323 + (3.*dg213 + 8.*dg312)*dg333 + 
           3.*(dg233*dg313 + dg133*dg323 + dg123*dg333))*ginv33) - 
     (dg222*dg312 + dg212*dg322)*pow2(ginv22) - 
     ((dg133 + 4.*dg313)*dg322 + (2.*dg213 + 8.*dg312)*dg323 + 
        dg122*dg333 + 2.*(dg233*dg312 + dg223*dg313 + dg123*dg323 + 
           dg212*dg333))*pow2(ginv23) - 
     (2.*dg133 + 8.*dg313)*dg333*pow2(ginv33)) + 
  ginv23*((ddg2333 + 2.*ddg3323)*ginv33 - 
     (2.*dg233 + 8.*dg323)*dg333*pow2(ginv33)) + 
  ginv12*((ddg1323 + ddg2313)*ginv33 - 
     ginv22*((2.*dg122*dg322 + 3.*(dg222*dg312 + dg212*dg322))*ginv23 + 
        (dg222*dg313 + (dg123 + dg213)*dg322 + dg122*dg323 + 
           2.*(dg223*dg312 + dg212*dg323))*ginv33) + 
     ginv23*(ddg1322 + ddg2312 - 
        (dg233*dg312 + dg133*dg322 + 
           4.*(dg313*(dg223 + dg322) + (dg123 + dg213 + dg312)*dg323) + 
           (dg122 + dg212)*dg333)*ginv33) + 
     ginv13*(ddg1312 + ddg2311 - 
        (dg222*dg311 + (dg122 + 4.*dg212)*dg312 + (dg112 + dg211)*dg322)*
         ginv22 - (dg133*dg312 + dg311*(dg233 + 4.*dg323) + 
           4.*((dg123 + dg213 + dg312)*dg313 + dg113*dg323) + 
           (dg112 + dg211)*dg333)*ginv33 - 
        ginv23*(2.*(dg223*dg311 + dg122*dg313 + dg113*dg322 + 
              dg211*dg323) + 4.*
            ((dg123 + dg213)*dg312 + dg212*dg313 + dg311*dg322 + 
              dg112*dg323 + pow2(dg312)))) - 
     (4.*dg311*dg312 + 2.*((dg123 + dg213)*dg311 + dg113*dg312 + 
           (dg112 + dg211)*dg313 + dg111*dg323))*pow2(ginv13) - 
     ((2.*dg213 + 4.*dg312)*dg322 + 
        2.*(dg223*dg312 + dg222*dg313 + dg123*dg322 + 
           (dg122 + dg212)*dg323))*pow2(ginv23) - 
     (dg133*dg323 + dg313*(dg233 + 4.*dg323) + (dg123 + dg213)*dg333)*
      pow2(ginv33)) + ginv11*(ddg1313*ginv33 - 
     ginv12*(((3.*dg112 + 2.*dg211)*dg311 + 3.*dg111*dg312)*ginv13 + 
        ((dg122 + dg212)*dg311 + (4.*dg112 + dg211)*dg312 + dg111*dg322)*
         ginv23 + ((dg123 + dg213)*dg311 + dg211*dg313 + 
           2.*(dg113*dg312 + dg112*dg313) + dg111*dg323)*ginv33) - 
     ginv22*(((dg122 + 2.*dg212)*dg312 + dg112*dg322)*ginv23 + 
        ((dg123 + dg213)*dg312 + dg212*dg313 + dg112*dg323)*ginv33) + 
     ginv13*(ddg1311 - (dg212*dg311 + (2.*dg112 + dg211)*dg312)*ginv22 - 
        ((dg123 + dg213)*dg311 + 4.*(dg113 + dg311)*dg312 + 
           (4.*dg112 + dg211)*dg313 + dg111*dg323)*ginv23 - 
        (6.*dg113*dg313 + dg311*(dg133 + 4.*dg313) + dg111*dg333)*ginv33) + 
     ginv23*(ddg1312 - (dg312*(dg133 + 4.*dg313) + 
           2.*((dg123 + dg213)*dg313 + dg113*dg323) + dg112*dg333)*ginv33) \
- (3.*(dg113*dg311 + dg111*dg313) + 2.*pow2(dg311))*pow2(ginv13) - 
     ((dg123 + dg213)*dg312 + (dg122 + dg212)*dg313 + dg113*dg322 + 
        dg112*dg323 + 2.*pow2(dg312))*pow2(ginv23) - 
     (dg133*dg313 + dg113*dg333 + 2.*pow2(dg313))*pow2(ginv33)) + 
  ginv22*(ddg2323*ginv33 + ginv23*
      (ddg2322 - (6.*dg223*dg323 + dg322*(dg233 + 4.*dg323) + dg222*dg333)*
         ginv33) - (3.*(dg223*dg322 + dg222*dg323) + 2.*pow2(dg322))*
      pow2(ginv23) - (dg233*dg323 + dg223*dg333 + 2.*pow2(dg323))*
      pow2(ginv33))
;

R11
=
dG11*g11 + dG12*g12 + dG13*g13 + gammado111*Gfromg1 + gammado112*Gfromg2 + 
  gammado113*Gfromg3 + (-0.5*ddg1111 + 3.*gamma111*gammado111 + 
     2.*(gamma211*gammado112 + gamma311*gammado113) + 
     gamma211*gammado211 + gamma311*gammado311)*ginv11 + 
  (-ddg1211 + 3.*(gamma112*gammado111 + gamma111*gammado112) + 
     2.*(gamma212*gammado112 + gamma312*gammado113 + 
        gamma211*gammado122 + gamma311*gammado123) + gamma212*gammado211 + 
     gamma211*gammado212 + gamma312*gammado311 + gamma311*gammado312)*ginv12 \
+ (-ddg1311 + 3.*(gamma113*gammado111 + gamma111*gammado113) + 
     2.*(gamma213*gammado112 + gamma313*gammado113 + 
        gamma211*gammado123 + gamma311*gammado133) + gamma213*gammado211 + 
     gamma211*gammado213 + gamma313*gammado311 + gamma311*gammado313)*ginv13 \
+ (-0.5*ddg2211 + 3.*gamma112*gammado112 + 
     2.*(gamma212*gammado122 + gamma312*gammado123) + 
     gamma212*gammado212 + gamma312*gammado312)*ginv22 + 
  (-ddg2311 + 3.*(gamma113*gammado112 + gamma112*gammado113) + 
     2.*(gamma213*gammado122 + (gamma212 + gamma313)*gammado123 + 
        gamma312*gammado133) + gamma213*gammado212 + gamma212*gammado213 + 
     gamma313*gammado312 + gamma312*gammado313)*ginv23 + 
  (-0.5*ddg3311 + 3.*gamma113*gammado113 + 
     2.*(gamma213*gammado123 + gamma313*gammado133) + gamma213*gammado213 + 
     gamma313*gammado313)*ginv33
;

R12
=
0.5*(dG21*g11 + (dG11 + dG22)*g12 + dG23*g13 + dG12*g22 + dG13*g23 + 
     (gammado112 + gammado211)*Gfromg1 + 
     (gammado122 + gammado212)*Gfromg2 + (gammado123 + gammado213)*Gfromg3) \
+ (-0.5*ddg1112 + gamma112*gammado111 + (gamma111 + gamma212)*gammado112 + 
     gamma312*gammado113 + gamma111*gammado211 + 2.*gamma211*gammado212 + 
     gamma311*(gammado213 + gammado312))*ginv11 + 
  (-ddg1212 + gamma122*gammado111 + (2.*gamma112 + gamma222)*gammado112 + 
     gamma322*gammado113 + (gamma111 + gamma212)*gammado122 + 
     gamma112*gammado211 + (gamma111 + 2.*gamma212)*gammado212 + 
     2.*gamma211*gammado222 + 
     gamma312*(gammado123 + gammado213 + gammado312) + 
     gamma311*(gammado223 + gammado322))*ginv12 + 
  (-ddg1312 + gamma123*gammado111 + (gamma113 + gamma223)*gammado112 + 
     (gamma112 + gamma323)*gammado113 + (gamma111 + gamma212)*gammado123 + 
     gamma312*gammado133 + gamma113*gammado211 + 
     (gamma111 + gamma313)*gammado213 + 
     2.*(gamma213*gammado212 + gamma211*gammado223) + 
     gamma313*gammado312 + gamma311*(gammado233 + gammado323))*ginv13 + 
  (-0.5*ddg2212 + gamma122*gammado112 + (gamma112 + gamma222)*gammado122 + 
     gamma322*gammado123 + gamma112*gammado212 + 2.*gamma212*gammado222 + 
     gamma312*(gammado223 + gammado322))*ginv22 + 
  (-ddg2312 + gamma123*gammado112 + gamma122*gammado113 + 
     (gamma113 + gamma223)*gammado122 + 
     (gamma112 + gamma222 + gamma323)*gammado123 + gamma322*gammado133 + 
     gamma113*gammado212 + gamma112*gammado213 + 
     2.*(gamma213*gammado222 + gamma212*gammado223) + 
     gamma313*(gammado223 + gammado322) + 
     gamma312*(gammado233 + gammado323))*ginv23 + 
  (-0.5*ddg3312 + gamma123*gammado113 + (gamma113 + gamma223)*gammado123 + 
     gamma323*gammado133 + gamma113*gammado213 + 2.*gamma213*gammado223 + 
     gamma313*(gammado233 + gammado323))*ginv33
;

R13
=
0.5*(dG31*g11 + dG32*g12 + (dG11 + dG33)*g13 + dG12*g23 + dG13*g33 + 
     (gammado113 + gammado311)*Gfromg1 + 
     (gammado123 + gammado312)*Gfromg2 + (gammado133 + gammado313)*Gfromg3) \
+ (-0.5*ddg1113 + gamma113*gammado111 + gamma213*gammado112 + 
     (gamma111 + gamma313)*gammado113 + gamma111*gammado311 + 
     gamma211*(gammado213 + gammado312) + 2.*gamma311*gammado313)*ginv11 + 
  (-ddg1213 + gamma123*gammado111 + (gamma113 + gamma223)*gammado112 + 
     (gamma112 + gamma323)*gammado113 + gamma213*gammado122 + 
     (gamma111 + gamma313)*gammado123 + gamma112*gammado311 + 
     gamma111*gammado312 + gamma212*(gammado213 + gammado312) + 
     gamma211*(gammado223 + gammado322) + 
     2.*(gamma312*gammado313 + gamma311*gammado323))*ginv12 + 
  (-ddg1313 + gamma133*gammado111 + gamma233*gammado112 + 
     (2.*gamma113 + gamma333)*gammado113 + 
     (gamma111 + gamma313)*gammado133 + gamma113*gammado311 + 
     gamma213*(gammado123 + gammado213 + gammado312) + 
     (gamma111 + 2.*gamma313)*gammado313 + 
     gamma211*(gammado233 + gammado323) + 2.*gamma311*gammado333)*ginv13 + 
  (-0.5*ddg2213 + gamma123*gammado112 + gamma223*gammado122 + 
     (gamma112 + gamma323)*gammado123 + gamma112*gammado312 + 
     gamma212*(gammado223 + gammado322) + 2.*gamma312*gammado323)*ginv22 + 
  (-ddg2313 + gamma133*gammado112 + gamma123*gammado113 + 
     gamma233*gammado122 + (gamma113 + gamma223 + gamma333)*gammado123 + 
     (gamma112 + gamma323)*gammado133 + gamma113*gammado312 + 
     gamma112*gammado313 + gamma213*(gammado223 + gammado322) + 
     gamma212*(gammado233 + gammado323) + 
     2.*(gamma313*gammado323 + gamma312*gammado333))*ginv23 + 
  (-0.5*ddg3313 + gamma133*gammado113 + gamma233*gammado123 + 
     (gamma113 + gamma333)*gammado133 + gamma113*gammado313 + 
     gamma213*(gammado233 + gammado323) + 2.*gamma313*gammado333)*ginv33
;

R22
=
dG21*g12 + dG22*g22 + dG23*g23 + gammado212*Gfromg1 + gammado222*Gfromg2 + 
  gammado223*Gfromg3 + (-0.5*ddg1122 + 
     gamma112*(gammado112 + 2.*gammado211) + 3.*gamma212*gammado212 + 
     gamma312*(2.*gammado213 + gammado312))*ginv11 + 
  (-ddg1222 + gamma122*(gammado112 + 2.*gammado211) + 
     gamma112*(gammado122 + 2.*gammado212) + 
     3.*(gamma222*gammado212 + gamma212*gammado222) + 
     2.*(gamma322*gammado213 + gamma312*gammado223) + 
     gamma322*gammado312 + gamma312*gammado322)*ginv12 + 
  (-ddg1322 + gamma123*(gammado112 + 2.*gammado211) + 
     gamma112*(gammado123 + 2.*gammado213) + 
     3.*(gamma223*gammado212 + gamma212*gammado223) + 
     2.*(gamma323*gammado213 + gamma312*gammado233) + 
     gamma323*gammado312 + gamma312*gammado323)*ginv13 + 
  (-0.5*ddg2222 + gamma122*(gammado122 + 2.*gammado212) + 
     3.*gamma222*gammado222 + gamma322*(2.*gammado223 + gammado322))*ginv22 \
+ (-ddg2322 + gamma123*(gammado122 + 2.*gammado212) + 
     gamma122*(gammado123 + 2.*gammado213) + 
     3.*(gamma223*gammado222 + gamma222*gammado223) + 
     2.*(gamma323*gammado223 + gamma322*gammado233) + 
     gamma323*gammado322 + gamma322*gammado323)*ginv23 + 
  (-0.5*ddg3322 + gamma123*(gammado123 + 2.*gammado213) + 
     3.*gamma223*gammado223 + gamma323*(2.*gammado233 + gammado323))*ginv33
;

R23
=
0.5*(dG31*g12 + dG21*g13 + dG32*g22 + (dG22 + dG33)*g23 + dG23*g33 + 
     (gammado213 + gammado312)*Gfromg1 + 
     (gammado223 + gammado322)*Gfromg2 + (gammado233 + gammado323)*Gfromg3) \
+ (-0.5*ddg1123 + gamma113*gammado211 + gamma213*gammado212 + 
     (gamma212 + gamma313)*gammado213 + 
     gamma112*(gammado113 + gammado311) + gamma212*gammado312 + 
     2.*gamma312*gammado313)*ginv11 + 
  (-ddg1223 + gamma123*gammado211 + (gamma113 + gamma223)*gammado212 + 
     (gamma222 + gamma323)*gammado213 + gamma213*gammado222 + 
     (gamma212 + gamma313)*gammado223 + 
     gamma122*(gammado113 + gammado311) + gamma222*gammado312 + 
     gamma112*(gammado123 + gammado312) + gamma212*gammado322 + 
     2.*(gamma322*gammado313 + gamma312*gammado323))*ginv12 + 
  (-ddg1323 + gamma133*gammado211 + gamma233*gammado212 + 
     (gamma113 + gamma223 + gamma333)*gammado213 + gamma213*gammado223 + 
     (gamma212 + gamma313)*gammado233 + 
     gamma123*(gammado113 + gammado311) + gamma223*gammado312 + 
     gamma112*(gammado133 + gammado313) + gamma212*gammado323 + 
     2.*(gamma323*gammado313 + gamma312*gammado333))*ginv13 + 
  (-0.5*ddg2223 + gamma123*gammado212 + gamma223*gammado222 + 
     (gamma222 + gamma323)*gammado223 + 
     gamma122*(gammado123 + gammado312) + gamma222*gammado322 + 
     2.*gamma322*gammado323)*ginv22 + 
  (-ddg2323 + gamma133*gammado212 + gamma233*gammado222 + 
     (2.*gamma223 + gamma333)*gammado223 + 
     (gamma222 + gamma323)*gammado233 + 
     gamma123*(gammado123 + gammado213 + gammado312) + 
     gamma122*(gammado133 + gammado313) + gamma223*gammado322 + 
     (gamma222 + 2.*gamma323)*gammado323 + 2.*gamma322*gammado333)*ginv23 + 
  (-0.5*ddg3323 + gamma133*gammado213 + gamma233*gammado223 + 
     (gamma223 + gamma333)*gammado233 + 
     gamma123*(gammado133 + gammado313) + gamma223*gammado323 + 
     2.*gamma323*gammado333)*ginv33
;

R33
=
dG31*g13 + dG32*g23 + dG33*g33 + gammado313*Gfromg1 + gammado323*Gfromg2 + 
  gammado333*Gfromg3 + (-0.5*ddg1133 + 
     gamma113*(gammado113 + 2.*gammado311) + 
     gamma213*(gammado213 + 2.*gammado312) + 3.*gamma313*gammado313)*ginv11 \
+ (-ddg1233 + gamma123*(gammado113 + 2.*gammado311) + 
     gamma113*(gammado123 + 2.*gammado312) + 
     gamma223*(gammado213 + 2.*gammado312) + 
     gamma213*(gammado223 + 2.*gammado322) + 
     3.*(gamma323*gammado313 + gamma313*gammado323))*ginv12 + 
  (-ddg1333 + gamma133*(gammado113 + 2.*gammado311) + 
     gamma233*(gammado213 + 2.*gammado312) + 
     gamma113*(gammado133 + 2.*gammado313) + 
     gamma213*(gammado233 + 2.*gammado323) + 
     3.*(gamma333*gammado313 + gamma313*gammado333))*ginv13 + 
  (-0.5*ddg2233 + gamma123*(gammado123 + 2.*gammado312) + 
     gamma223*(gammado223 + 2.*gammado322) + 3.*gamma323*gammado323)*ginv22 \
+ (-ddg2333 + gamma133*(gammado123 + 2.*gammado312) + 
     gamma123*(gammado133 + 2.*gammado313) + 
     gamma233*(gammado223 + 2.*gammado322) + 
     gamma223*(gammado233 + 2.*gammado323) + 
     3.*(gamma333*gammado323 + gamma323*gammado333))*ginv23 + 
  (-0.5*ddg3333 + gamma133*(gammado133 + 2.*gammado313) + 
     gamma233*(gammado233 + 2.*gammado323) + 3.*gamma333*gammado333)*ginv33
;

ff
=
chi
;

oochipsipower
=
1/chipsipower
;

f
=
oochipsipower*log(ff)
;

psim4
=
exp(-4.*f)
;

df1
=
(dchi1*oochipsipower)/chi
;

df2
=
(dchi2*oochipsipower)/chi
;

df3
=
(dchi3*oochipsipower)/chi
;

ddf11
=
(ddchi11*oochipsipower)/chi - chipsipower*pow2(df1)
;

ddf12
=
-(chipsipower*df1*df2) + (ddchi12*oochipsipower)/chi
;

ddf13
=
-(chipsipower*df1*df3) + (ddchi13*oochipsipower)/chi
;

ddf22
=
(ddchi22*oochipsipower)/chi - chipsipower*pow2(df2)
;

ddf23
=
-(chipsipower*df2*df3) + (ddchi23*oochipsipower)/chi
;

ddf33
=
(ddchi33*oochipsipower)/chi - chipsipower*pow2(df3)
;

cddf11
=
ddf11 - df1*gamma111 - df2*gamma211 - df3*gamma311
;

cddf12
=
ddf12 - df1*gamma112 - df2*gamma212 - df3*gamma312
;

cddf13
=
ddf13 - df1*gamma113 - df2*gamma213 - df3*gamma313
;

cddf22
=
ddf22 - df1*gamma122 - df2*gamma222 - df3*gamma322
;

cddf23
=
ddf23 - df1*gamma123 - df2*gamma223 - df3*gamma323
;

cddf33
=
ddf33 - df1*gamma133 - df2*gamma233 - df3*gamma333
;

trcddf
=
cddf11*ginv11 + cddf22*ginv22 + 
  2.*(cddf12*ginv12 + cddf13*ginv13 + cddf23*ginv23) + cddf33*ginv33
;

Rphi11
=
-2.*(cddf11 + g11*trcddf) + (4. - 4.*g11*ginv11)*pow2(df1) - 
  g11*(8.*(df1*(df2*ginv12 + df3*ginv13) + df2*df3*ginv23) + 
     4.*(ginv22*pow2(df2) + ginv33*pow2(df3)))
;

Rphi12
=
df1*df2*(4. - 8.*g12*ginv12) - 2.*(cddf12 + g12*trcddf) - 
  g12*(8.*df3*(df1*ginv13 + df2*ginv23) + 
     4.*(ginv11*pow2(df1) + ginv22*pow2(df2) + ginv33*pow2(df3)))
;

Rphi13
=
df1*(4.*df3 - 8.*df2*g13*ginv12) - 2.*(cddf13 + g13*trcddf) - 
  g13*(8.*df3*(df1*ginv13 + df2*ginv23) + 
     4.*(ginv11*pow2(df1) + ginv22*pow2(df2) + ginv33*pow2(df3)))
;

Rphi22
=
-2.*(cddf22 + g22*trcddf) + (4. - 4.*g22*ginv22)*pow2(df2) - 
  g22*(8.*(df1*(df2*ginv12 + df3*ginv13) + df2*df3*ginv23) + 
     4.*(ginv11*pow2(df1) + ginv33*pow2(df3)))
;

Rphi23
=
df2*(-8.*df1*g23*ginv12 + df3*(4. - 8.*g23*ginv23)) - 
  2.*(cddf23 + g23*trcddf) - g23*
   (8.*df1*df3*ginv13 + 4.*(ginv11*pow2(df1) + ginv22*pow2(df2) + 
        ginv33*pow2(df3)))
;

Rphi33
=
-2.*(cddf33 + g33*trcddf) - g33*
   (8.*(df1*(df2*ginv12 + df3*ginv13) + df2*df3*ginv23) + 
     4.*(ginv11*pow2(df1) + ginv22*pow2(df2))) + 
  (4. - 4.*g33*ginv33)*pow2(df3)
;

Rf11
=
R11 + Rphi11
;

Rf12
=
R12 + Rphi12
;

Rf13
=
R13 + Rphi13
;

Rf22
=
R22 + Rphi22
;

Rf23
=
R23 + Rphi23
;

Rf33
=
R33 + Rphi33
;

Rhat
=
psim4*(ginv11*Rf11 + ginv22*Rf22 + 
    2.*(ginv12*Rf12 + ginv13*Rf13 + ginv23*Rf23) + ginv33*Rf33)
;

cdda11
=
dda11 - da2*gamma211 - da3*gamma311 + 
  da1*(-gamma111 + df1*(-4. + 2.*g11*ginv11)) + 
  2.*g11*((da2*df1 + da1*df2)*ginv12 + (da3*df1 + da1*df3)*ginv13 + 
     da2*df2*ginv22 + (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33)
;

cdda12
=
dda12 - da1*gamma112 - da2*gamma212 - da3*gamma312 + 
  2.*(-(da2*df1) - da1*df2 + g12*
      (da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
        (da3*df1 + da1*df3)*ginv13 + da2*df2*ginv22 + 
        (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33))
;

cdda13
=
dda13 - da1*gamma113 - da2*gamma213 - da3*gamma313 + 
  2.*(-(da3*df1) - da1*df3 + g13*
      (da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
        (da3*df1 + da1*df3)*ginv13 + da2*df2*ginv22 + 
        (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33))
;

cdda22
=
dda22 - da1*gamma122 - da2*(4.*df2 + gamma222) - da3*gamma322 + 
  2.*g22*(da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
     (da3*df1 + da1*df3)*ginv13 + da2*df2*ginv22 + 
     (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33)
;

cdda23
=
dda23 - da1*gamma123 - da2*gamma223 - da3*gamma323 + 
  2.*(-(da3*df2) - da2*df3 + g23*
      (da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
        (da3*df1 + da1*df3)*ginv13 + da2*df2*ginv22 + 
        (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33))
;

cdda33
=
dda33 - da1*gamma133 - da2*gamma233 - da3*(4.*df3 + gamma333) + 
  2.*g33*(da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
     (da3*df1 + da1*df3)*ginv13 + da2*df2*ginv22 + 
     (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33)
;

trcdda
=
(cdda11*ginv11 + cdda22*ginv22 + 
    2.*(cdda12*ginv12 + cdda13*ginv13 + cdda23*ginv23) + cdda33*ginv33)*psim4
;

AA11
=
2.*(A11*(A12*ginv12 + A13*ginv13) + A12*A13*ginv23) + ginv11*pow2(A11) + 
  ginv22*pow2(A12) + ginv33*pow2(A13)
;

AA12
=
(A12*A13 + A11*A23)*ginv13 + A12*(A11*ginv11 + A22*ginv22) + 
  (A13*A22 + A12*A23)*ginv23 + A13*A23*ginv33 + ginv12*(A11*A22 + pow2(A12))
;

AA13
=
(A12*A13 + A11*A23)*ginv12 + A12*A23*ginv22 + (A13*A23 + A12*A33)*ginv23 + 
  A13*(A11*ginv11 + A33*ginv33) + ginv13*(A11*A33 + pow2(A13))
;

AA21
=
(A12*A13 + A11*A23)*ginv13 + A12*(A11*ginv11 + A22*ginv22) + 
  (A13*A22 + A12*A23)*ginv23 + A13*A23*ginv33 + ginv12*(A11*A22 + pow2(A12))
;

AA22
=
2.*(A12*(A22*ginv12 + A23*ginv13) + A22*A23*ginv23) + ginv11*pow2(A12) + 
  ginv22*pow2(A22) + ginv33*pow2(A23)
;

AA23
=
A12*A13*ginv11 + (A13*A22 + A12*A23)*ginv12 + (A13*A23 + A12*A33)*ginv13 + 
  A23*(A22*ginv22 + A33*ginv33) + ginv23*(A22*A33 + pow2(A23))
;

AA31
=
(A12*A13 + A11*A23)*ginv12 + A12*A23*ginv22 + (A13*A23 + A12*A33)*ginv23 + 
  A13*(A11*ginv11 + A33*ginv33) + ginv13*(A11*A33 + pow2(A13))
;

AA32
=
A12*A13*ginv11 + (A13*A22 + A12*A23)*ginv12 + (A13*A23 + A12*A33)*ginv13 + 
  A23*(A22*ginv22 + A33*ginv33) + ginv23*(A22*A33 + pow2(A23))
;

AA33
=
2.*(A13*(A23*ginv12 + A33*ginv13) + A23*A33*ginv23) + ginv11*pow2(A13) + 
  ginv22*pow2(A23) + ginv33*pow2(A33)
;

Ainv11
=
2.*(A23*ginv12*ginv13 + ginv11*(A12*ginv12 + A13*ginv13)) + 
  A11*pow2(ginv11) + A22*pow2(ginv12) + A33*pow2(ginv13)
;

Ainv12
=
ginv11*(A11*ginv12 + A12*ginv22 + A13*ginv23) + 
  ginv12*(A13*ginv13 + A22*ginv22 + A23*ginv23) + 
  ginv13*(A23*ginv22 + A33*ginv23) + A12*pow2(ginv12)
;

Ainv13
=
ginv11*(A11*ginv13 + A12*ginv23 + A13*ginv33) + 
  ginv12*(A12*ginv13 + A22*ginv23 + A23*ginv33) + 
  ginv13*(A23*ginv23 + A33*ginv33) + A13*pow2(ginv13)
;

Ainv22
=
2.*(A23*ginv22*ginv23 + ginv12*(A12*ginv22 + A13*ginv23)) + 
  A11*pow2(ginv12) + A22*pow2(ginv22) + A33*pow2(ginv23)
;

Ainv23
=
ginv13*(A12*ginv22 + A13*ginv23) + A33*ginv23*ginv33 + 
  ginv12*(A11*ginv13 + A12*ginv23 + A13*ginv33) + 
  ginv22*(A22*ginv23 + A23*ginv33) + A23*pow2(ginv23)
;

Ainv33
=
2.*(A23*ginv23*ginv33 + ginv13*(A12*ginv23 + A13*ginv33)) + 
  A11*pow2(ginv13) + A22*pow2(ginv23) + A33*pow2(ginv33)
;

cdA111
=
dA111 - 2.*(A11*gamma111 + A12*gamma211 + A13*gamma311)
;

cdA112
=
dA112 - A11*gamma112 - A22*gamma211 - A12*(gamma111 + gamma212) - 
  A23*gamma311 - A13*gamma312
;

cdA113
=
dA113 - A11*gamma113 - A23*gamma211 - A12*gamma213 - A33*gamma311 - 
  A13*(gamma111 + gamma313)
;

cdA122
=
dA122 - 2.*(A12*gamma112 + A22*gamma212 + A23*gamma312)
;

cdA123
=
dA123 - A13*gamma112 - A12*gamma113 - A22*gamma213 - A33*gamma312 - 
  A23*(gamma212 + gamma313)
;

cdA133
=
dA133 - 2.*(A13*gamma113 + A23*gamma213 + A33*gamma313)
;

cdA211
=
dA211 - 2.*(A11*gamma112 + A12*gamma212 + A13*gamma312)
;

cdA212
=
dA212 - A11*gamma122 - A22*gamma212 - A12*(gamma112 + gamma222) - 
  A23*gamma312 - A13*gamma322
;

cdA213
=
dA213 - A11*gamma123 - A23*gamma212 - A12*gamma223 - A33*gamma312 - 
  A13*(gamma112 + gamma323)
;

cdA222
=
dA222 - 2.*(A12*gamma122 + A22*gamma222 + A23*gamma322)
;

cdA223
=
dA223 - A13*gamma122 - A12*gamma123 - A22*gamma223 - A33*gamma322 - 
  A23*(gamma222 + gamma323)
;

cdA233
=
dA233 - 2.*(A13*gamma123 + A23*gamma223 + A33*gamma323)
;

cdA311
=
dA311 - 2.*(A11*gamma113 + A12*gamma213 + A13*gamma313)
;

cdA312
=
dA312 - A11*gamma123 - A22*gamma213 - A12*(gamma113 + gamma223) - 
  A23*gamma313 - A13*gamma323
;

cdA313
=
dA313 - A11*gamma133 - A23*gamma213 - A12*gamma233 - A33*gamma313 - 
  A13*(gamma113 + gamma333)
;

cdA322
=
dA322 - 2.*(A12*gamma123 + A22*gamma223 + A23*gamma323)
;

cdA323
=
dA323 - A13*gamma123 - A12*gamma133 - A22*gamma233 - A33*gamma323 - 
  A23*(gamma223 + gamma333)
;

cdA333
=
dA333 - 2.*(A13*gamma133 + A23*gamma233 + A33*gamma333)
;

divbeta
=
db11 + db22 + db33
;

totdivbeta
=
0.66666666666666666667*divbeta
;

lieg11
=
beta1*dg111 + beta2*dg211 + beta3*dg311 + 
  2.*(db11*g11 + db12*g12 + db13*g13) - g11*totdivbeta
;

lieg12
=
beta1*dg112 + beta2*dg212 + beta3*dg312 + db21*g11 + db23*g13 + db12*g22 + 
  db13*g23 + g12*(db11 + db22 - totdivbeta)
;

lieg13
=
beta1*dg113 + beta2*dg213 + beta3*dg313 + db31*g11 + db32*g12 + db12*g23 + 
  db13*g33 + g13*(db11 + db33 - totdivbeta)
;

lieg22
=
beta1*dg122 + beta2*dg222 + beta3*dg322 + 
  2.*(db21*g12 + db22*g22 + db23*g23) - g22*totdivbeta
;

lieg23
=
beta1*dg123 + beta2*dg223 + beta3*dg323 + db31*g12 + db21*g13 + db32*g22 + 
  db23*g33 + g23*(db22 + db33 - totdivbeta)
;

lieg33
=
beta1*dg133 + beta2*dg233 + beta3*dg333 + 
  2.*(db31*g13 + db32*g23 + db33*g33) - g33*totdivbeta
;

lieA11
=
beta1*dA111 + beta2*dA211 + beta3*dA311 + 
  2.*(A11*db11 + A12*db12 + A13*db13) - A11*totdivbeta
;

lieA12
=
beta1*dA112 + beta2*dA212 + beta3*dA312 + A22*db12 + A23*db13 + A11*db21 + 
  A13*db23 + A12*(db11 + db22 - totdivbeta)
;

lieA13
=
beta1*dA113 + beta2*dA213 + beta3*dA313 + A23*db12 + A33*db13 + A11*db31 + 
  A12*db32 + A13*(db11 + db33 - totdivbeta)
;

lieA22
=
beta1*dA122 + beta2*dA222 + beta3*dA322 + 
  2.*(A12*db21 + A22*db22 + A23*db23) - A22*totdivbeta
;

lieA23
=
beta1*dA123 + beta2*dA223 + beta3*dA323 + A13*db21 + A33*db23 + A12*db31 + 
  A22*db32 + A23*(db22 + db33 - totdivbeta)
;

lieA33
=
beta1*dA133 + beta2*dA233 + beta3*dA333 + 
  2.*(A13*db31 + A23*db32 + A33*db33) - A33*totdivbeta
;

betas
=
beta1*sdown1 + beta2*sdown2 + beta3*sdown3
;

Dbetas
=
(db11*sdown1 + db12*sdown2 + db13*sdown3)*sup1 + 
  (db21*sdown1 + db22*sdown2 + db23*sdown3)*sup2 + 
  (db31*sdown1 + db32*sdown2 + db33*sdown3)*sup3
;

Dalpha
=
da1*sup1 + da2*sup2 + da3*sup3
;

DKhat
=
dKhat1*sup1 + dKhat2*sup2 + dKhat3*sup3
;

DK
=
dK1*sup1 + dK2*sup2 + dK3*sup3
;

DTheta
=
dTheta1*sup1 + dTheta2*sup2 + dTheta3*sup3
;

Gams
=
G1*sdown1 + G2*sdown2 + G3*sdown3
;

DGams
=
(dG11*sdown1 + dG12*sdown2 + dG13*sdown3)*sup1 + 
  (dG21*sdown1 + dG22*sdown2 + dG23*sdown3)*sup2 + 
  (dG31*sdown1 + dG32*sdown2 + dG33*sdown3)*sup3
;

GamA1
=
G1*qud11 + G2*qud12 + G3*qud13
;

GamA2
=
G1*qud21 + G2*qud22 + G3*qud23
;

GamA3
=
G1*qud31 + G2*qud32 + G3*qud33
;

DGamA1
=
(dG11*qud11 + dG12*qud12 + dG13*qud13)*sup1 + 
  (dG21*qud11 + dG22*qud12 + dG23*qud13)*sup2 + 
  (dG31*qud11 + dG32*qud12 + dG33*qud13)*sup3
;

DGamA2
=
(dG11*qud21 + dG12*qud22 + dG13*qud23)*sup1 + 
  (dG21*qud21 + dG22*qud22 + dG23*qud23)*sup2 + 
  (dG31*qud21 + dG32*qud22 + dG33*qud23)*sup3
;

DGamA3
=
(dG11*qud31 + dG12*qud32 + dG13*qud33)*sup1 + 
  (dG21*qud31 + dG22*qud32 + dG23*qud33)*sup2 + 
  (dG31*qud31 + dG32*qud32 + dG33*qud33)*sup3
;

betaA1
=
beta1*qud11 + beta2*qud12 + beta3*qud13
;

betaA2
=
beta1*qud21 + beta2*qud22 + beta3*qud23
;

betaA3
=
beta1*qud31 + beta2*qud32 + beta3*qud33
;

DbetaA1
=
(db11*qud11 + db12*qud12 + db13*qud13)*sup1 + 
  (db21*qud11 + db22*qud12 + db23*qud13)*sup2 + 
  (db31*qud11 + db32*qud12 + db33*qud13)*sup3
;

DbetaA2
=
(db11*qud21 + db12*qud22 + db13*qud23)*sup1 + 
  (db21*qud21 + db22*qud22 + db23*qud23)*sup2 + 
  (db31*qud21 + db32*qud22 + db33*qud23)*sup3
;

DbetaA3
=
(db11*qud31 + db12*qud32 + db13*qud33)*sup1 + 
  (db21*qud31 + db22*qud32 + db23*qud33)*sup2 + 
  (db31*qud31 + db32*qud32 + db33*qud33)*sup3
;

lienKhat
=
-((DKhat + Khat/r)*sqrt(muL))
;

lienTheta
=
-DTheta - (kappa1*(2. + kappa2) + 1/r)*Theta
;

lienK
=
lienKhat + 2.*lienTheta
;

rKhat
=
beta1*dKhat1 + beta2*dKhat2 + beta3*dKhat3 + alpha*lienKhat
;

#if 0
// David's new version
rGams
=
(beta1*dG11 + beta2*dG21 + beta3*dG31 + 
     (ddb111*quu11 + ddb221*quu22 + 
        2.*(ddb121*quu12 + ddb131*quu13 + ddb231*quu23) + ddb331*quu33)/chi\
)*sdown1 + (beta1*dG12 + beta2*dG22 + beta3*dG32 + 
     (ddb112*quu11 + ddb222*quu22 + 
        2.*(ddb122*quu12 + ddb132*quu13 + ddb232*quu23) + ddb332*quu33)/chi\
)*sdown2 + (beta1*dG13 + beta2*dG23 + beta3*dG33 + 
     (ddb113*quu11 + ddb223*quu22 + 
        2.*(ddb123*quu12 + ddb133*quu13 + ddb233*quu23) + ddb333*quu33)/chi\
)*sdown3 - ((ddb111*qud11 + ddb112*qud12 + ddb113*qud13 + ddb121*qud21 + 
        ddb122*qud22 + ddb123*qud23 + ddb131*qud31 + ddb132*qud32 + 
        ddb133*qud33)*sup1 + (ddb121*qud11 + ddb122*qud12 + 
        ddb123*qud13 + ddb221*qud21 + ddb222*qud22 + ddb223*qud23 + 
        ddb231*qud31 + ddb232*qud32 + ddb233*qud33)*sup2 + 
     (ddb131*qud11 + ddb132*qud12 + ddb133*qud13 + ddb231*qud21 + 
        ddb232*qud22 + ddb233*qud23 + ddb331*qud31 + ddb332*qud32 + 
        ddb333*qud33)*sup3)/chi - (dG11 + dG22 + dG33)*vbetas + 
  2.*((0.33333333333333333333*alpha*
        (dTheta1*sup1 + dTheta2*sup2 + dTheta3*sup3))/(chi + chi*vbetas) + 
     ((db11 + db22 + db33)*shiftdriver)/(vbetaA*sqrt(3.))) + 
  (1.3333333333333333333*alpha*(dKhat1*sup1 + dKhat2*sup2 + dKhat3*sup3)*
     sqrt(muL))/(chi*(vbetas + sqrt(muL)))
;
#else
//David's old version
rGams
=
shiftdriver*((beta1*db11 + beta2*db21)*(db12*sdown2 + db13*sdown3) + 
     2.*beta1*((beta2*ddb121 + beta3*ddb131)*sdown1 + 
        (beta2*ddb122 + beta3*ddb132)*sdown2 + 
        (beta2*ddb123 + beta3*ddb133)*sdown3) + 
     sdown1*(db21*(beta1*db12 + beta2*(db11 + db22) + beta3*db32) + 
        db31*(beta1*db13 + beta2*db23 + beta3*(db11 + db33)) + 
        beta2*(2.*beta3*ddb231 + dG21) + beta3*dG31 + ddb111*pow2(beta1) + 
        ddb221*pow2(beta2) + ddb331*pow2(beta3) + beta1*(dG11 + pow2(db11))\
) + sdown2*(db12*(beta1*db22 + beta3*db31) + 
        db32*(beta1*db13 + beta2*db23 + beta3*(db22 + db33)) + beta1*dG12 + 
        beta3*dG32 + ddb112*pow2(beta1) + ddb222*pow2(beta2) + 
        ddb332*pow2(beta3) + beta2*(2.*beta3*ddb232 + dG22 + pow2(db22)))) - 
  ((beta1*db11 + beta2*db21 + beta3*db31)*sdown1 + 
     (beta2*db22 + beta3*db32)*sdown2 + beta2*db23*sdown3 + 
     beta1*(db12*sdown2 + db13*sdown3))*pow2(shiftdriver) + 
  sdown3*(shiftdriver*((beta1*db12 + beta2*db22)*db23 + beta1*dG13 + 
        beta2*dG23 + ddb113*pow2(beta1) + ddb223*pow2(beta2) + 
        ddb333*pow2(beta3) + beta3*
         (db13*db31 + db23*db32 + 2.*beta2*ddb233 + dG33 + pow2(db33))) + 
     db33*((beta1*db13 + beta2*db23)*shiftdriver - beta3*pow2(shiftdriver)))
;
#endif

rTheta
=
beta1*dTheta1 + beta2*dTheta2 + beta3*dTheta3 + alpha*lienTheta
;

rACss
=
2.*((A23*alpha*K + lieA23)*sup2*sup3 + 
     sup1*((A12*alpha*K + lieA12)*sup2 + A13*alpha*K*sup3) + 
     psim4*((-cdda23 + alpha*Rf23)*sup2*sup3 + 
        sup1*((-cdda12 + alpha*Rf12)*sup2 - cdda13*sup3))) + 
  0.66666666666666666667*(g13*sup1 + g23*sup2)*sup3*trcdda + 
  sup1*(2.*(-(AA31*alpha) + lieA13)*sup3 + 
     0.66666666666666666667*g12*sup2*trcdda) + 
  (lieA11 + psim4*(-cdda11 + alpha*Rf11) + 
     0.33333333333333333333*g11*(-(alpha*Rhat) + trcdda))*pow2(sup1) + 
  (lieA22 - cdda22*psim4 + alpha*
      (A22*K + psim4*Rf22 - 0.33333333333333333333*g22*Rhat) + 
     0.33333333333333333333*g22*trcdda)*pow2(sup2) + 
  (lieA33 - cdda33*psim4 + alpha*
      (A33*K + psim4*Rf33 - 0.33333333333333333333*g33*Rhat) + 
     0.33333333333333333333*g33*trcdda)*pow2(sup3) + 
  alpha*(ginv11*((-2.*cdA111*chi + 3.*A11*dchi1)*sup1 + 
        (-2.*cdA112*chi + 3.*A12*dchi1)*sup2 + 
        (-2.*cdA113*chi + 3.*A13*dchi1)*sup3) + 
     ginv22*((-2.*cdA212*chi + 3.*A12*dchi2)*sup1 + 
        (-2.*cdA222*chi + 3.*A22*dchi2)*sup2 + 
        (-2.*cdA223*chi + 3.*A23*dchi2)*sup3) + 
     ginv33*((-2.*cdA313*chi + 3.*A13*dchi3)*sup1 + 
        (-2.*cdA323*chi + 3.*A23*dchi3)*sup2 + 
        (-2.*cdA333*chi + 3.*A33*dchi3)*sup3) + 
     chi*(-2.*DTheta + 1.3333333333333333333*
         (dK1*sup1 + dK2*sup2 + dK3*sup3)) + 
     ginv12*((-2.*cdA212*chi + 3.*A12*dchi2)*sup2 + 
        (-2.*cdA213*chi + 3.*A13*dchi2)*sup3 - 
        2.*chi*((cdA112 + cdA211)*sup1 + cdA122*sup2 + cdA123*sup3) + 
        3.*((A12*dchi1 + A11*dchi2)*sup1 + dchi1*(A22*sup2 + A23*sup3))) + 
     ginv13*((-2.*cdA312*chi + 3.*A12*dchi3)*sup2 + 
        (-2.*cdA313*chi + 3.*A13*dchi3)*sup3 - 
        2.*chi*((cdA113 + cdA311)*sup1 + cdA123*sup2 + cdA133*sup3) + 
        3.*((A13*dchi1 + A11*dchi3)*sup1 + dchi1*(A23*sup2 + A33*sup3))) + 
     ginv23*((-2.*cdA322*chi + 3.*A22*dchi3)*sup2 + 
        (-2.*cdA323*chi + 3.*A23*dchi3)*sup3 - 
        2.*chi*((cdA213 + cdA312)*sup1 + cdA223*sup2 + cdA233*sup3) + 
        3.*((A13*dchi2 + A12*dchi3)*sup1 + dchi2*(A23*sup2 + A33*sup3))) + 
     (0.33333333333333333333*((dG11 - dGfromgdu11)*qud11 + 
           (dG12 - dGfromgdu12)*qud12 + (dG13 - dGfromgdu13)*qud13 + 
           (dG21 - dGfromgdu21)*qud21 + (dG22 - dGfromgdu22)*qud22 + 
           (dG23 - dGfromgdu23)*qud23 + (dG31 - dGfromgdu31)*qud31 + 
           (dG32 - dGfromgdu32)*qud32 + (dG33 - dGfromgdu33)*qud33) + 
        kappa1*((G1 - Gfromg1)*sdown1 + (G2 - Gfromg2)*sdown2 + 
           (G3 - Gfromg3)*sdown3) + 
        0.66666666666666666667*
         ((dGfromgdu21*sdown1 + dGfromgdu22*sdown2)*sup2 + 
           sdown3*((-dG13 + dGfromgdu13)*sup1 - dG23*sup2 - dG33*sup3) + 
           sdown1*((-dG11 + dGfromgdu11)*sup1 - dG21*sup2 - dG31*sup3 + 
              dGfromgdu31*sup3) + 
           sdown2*((-dG12 + dGfromgdu12)*sup1 - dG22*sup2 - dG32*sup3 + 
              dGfromgdu32*sup3)))*pow2(chi) + 
     0.66666666666666666667*sup2*
      (-(Rhat*(g12*sup1 + g23*sup3)) + dGfromgdu23*sdown3*pow2(chi)) + 
     sup3*((2.*psim4*Rf13 - 0.66666666666666666667*g13*Rhat)*sup1 + 
        0.66666666666666666667*dGfromgdu33*sdown3*pow2(chi)) + 
     (-2.*AA11 + A11*K)*pow2(sup1) - 
     2.*((AA23 + AA32)*sup2*sup3 + sup1*((AA12 + AA21)*sup2 + AA13*sup3) + 
        AA22*pow2(sup2) + AA33*pow2(sup3)))
;

rACqq
=
chi*(-((4.*(A12*Ainv12 + A13*Ainv13 + A23*Ainv23) + 
          2.*(A11*Ainv11 + A22*Ainv22 + A33*Ainv33))*alpha) + 
     Ainv11*lieg11 + Ainv22*lieg22 + 
     2.*(Ainv12*lieg12 + Ainv13*lieg13 + Ainv23*lieg23) + Ainv33*lieg33) - 
  rACss
;

rGamA1
=
-(((dG21*qud11 + dG22*qud12 + dG23*qud13)*sup2 + 
       (dG31*qud11 + dG32*qud12 + dG33*qud13)*sup3)*vbetaA) + 
  qud11*(beta2*dG21 + beta3*dG31 + 
     (1.3333333333333333333*ddb111*quu11 + 
        2.3333333333333333333*(ddb121*quu12 + ddb131*quu13) + 
        ddb221*quu22 + ddb331*quu33 + 
        (shiftdriver*(db11*sup1 + db31*sup3))/vbetaA)/chi + 
     dG11*(beta1 - sup1*vbetaA)) + 
  qud12*(beta2*dG22 + beta3*dG32 + 
     (1.3333333333333333333*ddb112*quu11 + 
        2.3333333333333333333*(ddb122*quu12 + ddb132*quu13) + 
        ddb222*quu22 + 2.*ddb232*quu23 + ddb332*quu33 + 
        (shiftdriver*(db12*sup1 + db22*sup2 + db32*sup3))/vbetaA)/chi + 
     dG12*(beta1 - sup1*vbetaA)) + 
  qud13*(beta2*dG23 + beta3*dG33 + 
     (1.3333333333333333333*ddb113*quu11 + 
        2.3333333333333333333*(ddb123*quu12 + ddb133*quu13) + 
        ddb223*quu22 + 2.*ddb233*quu23 + ddb333*quu33 + 
        (shiftdriver*(db13*sup1 + db23*sup2 + db33*sup3))/vbetaA)/chi + 
     dG13*(beta1 - sup1*vbetaA)) + 
  (0.33333333333333333333*((ddb121*qud21 + ddb122*qud22 + ddb123*qud23 + 
           ddb131*qud31 + ddb132*qud32 + ddb133*qud33)*quu11 + 
        (ddb221*qud21 + ddb223*qud23 + ddb231*qud31 + ddb232*qud32 + 
           ddb233*qud33)*quu12 + 
        (ddb231*qud21 + ddb232*qud22 + ddb233*qud23 + ddb331*qud31 + 
           ddb332*qud32)*quu13) - 
     alpha*((1.3333333333333333333*dKhat1 + 
           0.66666666666666666667*dTheta1)*quu11 + 
        1.3333333333333333333*(dKhat2*quu12 + dKhat3*quu13)) + 
     1.3333333333333333333*((ddb132*quu13*sdown2 + ddb113*quu11*sdown3)*
         sup1 + (quu13*(ddb231*sdown1 + ddb232*sdown2) + 
           quu12*(ddb222*sdown2 + ddb223*sdown3))*sup2 + 
        (quu12*(ddb232*sdown2 + ddb233*sdown3) + 
           quu13*(ddb331*sdown1 + ddb332*sdown2 + ddb333*sdown3))*sup3 + 
        sdown1*((ddb121*quu12 + ddb131*quu13)*sup1 + ddb221*quu12*sup2 + 
           ddb131*quu11*sup3) + 
        sdown2*((ddb112*quu11 + ddb122*quu12)*sup1 + 
           quu11*(ddb122*sup2 + ddb132*sup3)) + 
        sdown3*((ddb123*quu12 + ddb133*quu13)*sup1 + 
           quu11*(ddb123*sup2 + ddb133*sup3))) + 
     qud11*(2.*ddb231*quu23 + (db21*shiftdriver*sup2)/vbetaA) - 
     (((db11*quu11 + db21*quu12)*sdown1 + 
          (db12*quu11 + db22*quu12 + db32*quu13)*sdown2 + 
          (db13*quu11 + db23*quu12 + db33*quu13)*sdown3)*shiftdriver)/
      vbetaA + ((dG22*quu12 + dG32*quu13)*sdown2 + 
        (dG13*quu11 + dG23*quu12)*sdown3)*vbetaA + 
     quu11*(1.3333333333333333333*sdown1*(ddb111*sup1 + ddb121*sup2) + 
        (dG11*sdown1 + dG12*sdown2)*vbetaA) + 
     quu12*(-0.66666666666666666667*alpha*dTheta2 + 
        0.33333333333333333333*ddb222*qud22 + 
        sdown1*(1.3333333333333333333*ddb231*sup3 + dG21*vbetaA)) + 
     quu13*(-0.66666666666666666667*alpha*dTheta3 + 
        0.33333333333333333333*ddb333*qud33 - 
        (db31*sdown1*shiftdriver)/vbetaA + dG31*sdown1*vbetaA + 
        sdown3*(1.3333333333333333333*ddb233*sup2 + dG33*vbetaA)))/chi
;

rGamA2
=
-(((dG21*qud21 + dG22*qud22 + dG23*qud23)*sup2 + 
       (dG31*qud21 + dG32*qud22 + dG33*qud23)*sup3)*vbetaA) + 
  qud21*(beta2*dG21 + beta3*dG31 + 
     (ddb111*quu11 + 2.*ddb131*quu13 + 
        1.3333333333333333333*ddb221*quu22 + 
        2.3333333333333333333*(ddb121*quu12 + ddb231*quu23) + 
        ddb331*quu33 + (shiftdriver*(db11*sup1 + db31*sup3))/vbetaA)/chi + 
     dG11*(beta1 - sup1*vbetaA)) + 
  qud22*(beta2*dG22 + beta3*dG32 + 
     (ddb112*quu11 + 2.*ddb132*quu13 + 
        1.3333333333333333333*ddb222*quu22 + 
        2.3333333333333333333*(ddb122*quu12 + ddb232*quu23) + 
        ddb332*quu33 + (shiftdriver*(db12*sup1 + db22*sup2 + db32*sup3))/
         vbetaA)/chi + dG12*(beta1 - sup1*vbetaA)) + 
  qud23*(beta2*dG23 + beta3*dG33 + 
     (ddb113*quu11 + 2.*ddb133*quu13 + 
        1.3333333333333333333*ddb223*quu22 + 
        2.3333333333333333333*(ddb123*quu12 + ddb233*quu23) + 
        ddb333*quu33 + (shiftdriver*(db13*sup1 + db23*sup2 + db33*sup3))/
         vbetaA)/chi + dG13*(beta1 - sup1*vbetaA)) + 
  (0.33333333333333333333*((ddb111*qud11 + ddb112*qud12 + ddb113*qud13 + 
           ddb131*qud31 + ddb132*qud32 + ddb133*qud33)*quu12 + 
        (ddb121*qud11 + ddb123*qud13 + ddb231*qud31 + ddb232*qud32 + 
           ddb233*qud33)*quu22 + 
        (ddb131*qud11 + ddb132*qud12 + ddb133*qud13 + ddb331*qud31 + 
           ddb332*qud32)*quu23) - 
     alpha*((1.3333333333333333333*dKhat1 + 
           0.66666666666666666667*dTheta1)*quu12 + 
        1.3333333333333333333*(dKhat2*quu22 + dKhat3*quu23)) + 
     1.3333333333333333333*((ddb132*quu23*sdown2 + ddb113*quu12*sdown3)*
         sup1 + (quu23*(ddb231*sdown1 + ddb232*sdown2) + 
           quu22*(ddb222*sdown2 + ddb223*sdown3))*sup2 + 
        (quu22*(ddb232*sdown2 + ddb233*sdown3) + 
           quu23*(ddb331*sdown1 + ddb332*sdown2 + ddb333*sdown3))*sup3 + 
        sdown1*((ddb121*quu22 + ddb131*quu23)*sup1 + ddb221*quu22*sup2 + 
           ddb131*quu12*sup3) + 
        sdown2*((ddb112*quu12 + ddb122*quu22)*sup1 + 
           quu12*(ddb122*sup2 + ddb132*sup3)) + 
        sdown3*((ddb123*quu22 + ddb133*quu23)*sup1 + 
           quu12*(ddb123*sup2 + ddb133*sup3))) - 
     (((db11*quu12 + db21*quu22)*sdown1 + 
          (db12*quu12 + db22*quu22 + db32*quu23)*sdown2 + 
          (db13*quu12 + db23*quu22 + db33*quu23)*sdown3)*shiftdriver)/
      vbetaA + (db21*qud21*shiftdriver*sup2)/vbetaA + 
     ((dG22*quu22 + dG32*quu23)*sdown2 + (dG13*quu12 + dG23*quu22)*sdown3)*
      vbetaA + quu12*(1.3333333333333333333*sdown1*
         (ddb111*sup1 + ddb121*sup2) + (dG11*sdown1 + dG12*sdown2)*vbetaA) \
+ quu22*(-0.66666666666666666667*alpha*dTheta2 + 
        0.33333333333333333333*ddb122*qud12 + 
        sdown1*(1.3333333333333333333*ddb231*sup3 + dG21*vbetaA)) + 
     quu23*(-0.66666666666666666667*alpha*dTheta3 + 
        0.33333333333333333333*ddb333*qud33 - 
        (db31*sdown1*shiftdriver)/vbetaA + dG31*sdown1*vbetaA + 
        sdown3*(1.3333333333333333333*ddb233*sup2 + dG33*vbetaA)))/chi
;

rGamA3
=
-(((dG21*qud31 + dG22*qud32 + dG23*qud33)*sup2 + 
       (dG31*qud31 + dG32*qud32 + dG33*qud33)*sup3)*vbetaA) + 
  qud31*(beta2*dG21 + beta3*dG31 + 
     (ddb111*quu11 + 2.*ddb121*quu12 + ddb221*quu22 + 
        2.3333333333333333333*(ddb131*quu13 + ddb231*quu23) + 
        1.3333333333333333333*ddb331*quu33 + 
        (shiftdriver*(db11*sup1 + db31*sup3))/vbetaA)/chi + 
     dG11*(beta1 - sup1*vbetaA)) + 
  qud32*(beta2*dG22 + beta3*dG32 + 
     (ddb112*quu11 + 2.*ddb122*quu12 + ddb222*quu22 + 
        2.3333333333333333333*(ddb132*quu13 + ddb232*quu23) + 
        1.3333333333333333333*ddb332*quu33 + 
        (shiftdriver*(db12*sup1 + db22*sup2 + db32*sup3))/vbetaA)/chi + 
     dG12*(beta1 - sup1*vbetaA)) + 
  qud33*(beta2*dG23 + beta3*dG33 + 
     (ddb113*quu11 + 2.*ddb123*quu12 + ddb223*quu22 + 
        2.3333333333333333333*(ddb133*quu13 + ddb233*quu23) + 
        1.3333333333333333333*ddb333*quu33 + 
        (shiftdriver*(db13*sup1 + db23*sup2 + db33*sup3))/vbetaA)/chi + 
     dG13*(beta1 - sup1*vbetaA)) + 
  (0.33333333333333333333*((ddb111*qud11 + ddb112*qud12 + ddb113*qud13 + 
           ddb121*qud21 + ddb122*qud22 + ddb123*qud23)*quu13 + 
        (ddb121*qud11 + ddb123*qud13 + ddb221*qud21 + ddb222*qud22 + 
           ddb223*qud23)*quu23 + 
        (ddb131*qud11 + ddb132*qud12 + ddb133*qud13 + ddb231*qud21 + 
           ddb232*qud22)*quu33) - 
     alpha*((1.3333333333333333333*dKhat1 + 
           0.66666666666666666667*dTheta1)*quu13 + 
        1.3333333333333333333*(dKhat2*quu23 + dKhat3*quu33)) + 
     1.3333333333333333333*((ddb132*quu33*sdown2 + ddb113*quu13*sdown3)*
         sup1 + (quu33*(ddb231*sdown1 + ddb232*sdown2) + 
           quu23*(ddb222*sdown2 + ddb223*sdown3))*sup2 + 
        (quu23*(ddb232*sdown2 + ddb233*sdown3) + 
           quu33*(ddb331*sdown1 + ddb332*sdown2 + ddb333*sdown3))*sup3 + 
        sdown1*((ddb121*quu23 + ddb131*quu33)*sup1 + ddb221*quu23*sup2 + 
           ddb131*quu13*sup3) + 
        sdown2*((ddb112*quu13 + ddb122*quu23)*sup1 + 
           quu13*(ddb122*sup2 + ddb132*sup3)) + 
        sdown3*((ddb123*quu23 + ddb133*quu33)*sup1 + 
           quu13*(ddb123*sup2 + ddb133*sup3))) - 
     (((db11*quu13 + db21*quu23)*sdown1 + 
          (db12*quu13 + db22*quu23 + db32*quu33)*sdown2 + 
          (db13*quu13 + db23*quu23 + db33*quu33)*sdown3)*shiftdriver)/
      vbetaA + (db21*qud31*shiftdriver*sup2)/vbetaA + 
     ((dG22*quu23 + dG32*quu33)*sdown2 + (dG13*quu13 + dG23*quu23)*sdown3)*
      vbetaA + quu13*(1.3333333333333333333*sdown1*
         (ddb111*sup1 + ddb121*sup2) + (dG11*sdown1 + dG12*sdown2)*vbetaA) \
+ quu33*(-0.66666666666666666667*alpha*dTheta3 + 
        ddb233*(0.33333333333333333333*qud23 + 
           1.3333333333333333333*sdown3*sup2) - 
        (db31*sdown1*shiftdriver)/vbetaA + 
        (dG31*sdown1 + dG33*sdown3)*vbetaA) + 
     quu23*(-0.66666666666666666667*alpha*dTheta2 + 
        0.33333333333333333333*ddb122*qud12 + 
        sdown1*(1.3333333333333333333*ddb231*sup3 + dG21*vbetaA)))/chi
;

rACsA1
=
(qud11*(lieA11 + alpha*chi*Rf11) + 
     qud21*(lieA12 + alpha*(-2.*AA12 + chi*Rf12)) + 
     qud31*(lieA13 + alpha*(-2.*AA13 + chi*Rf13)))*sup1 + 
  qud11*((-(cdda11*chi) + A11*alpha*K)*sup1 + lieA12*sup2 + 
     (A13*alpha*K + lieA13)*sup3 + 
     alpha*((-2.*AA21 + A12*K)*sup2 - 2.*AA31*sup3)) + 
  qud21*((-(cdda12*chi) + A12*alpha*K)*sup1 + lieA22*sup2 + 
     (A23*alpha*K + lieA23)*sup3 + 
     alpha*((-2.*AA22 + A22*K)*sup2 - 2.*AA32*sup3)) + 
  qud31*((-(cdda13*chi) + A13*alpha*K)*sup1 + lieA23*sup2 + 
     (A33*alpha*K + lieA33)*sup3 + 
     alpha*((-2.*AA23 + A23*K)*sup2 - 2.*AA33*sup3)) + 
  alpha*(ginv11*((-(cdA111*chi) + 1.5*A11*dchi1)*qud11 + 
        (-(cdA112*chi) + 1.5*A12*dchi1)*qud21 + 
        (-(cdA113*chi) + 1.5*A13*dchi1)*qud31) + 
     ginv22*((-(cdA212*chi) + 1.5*A12*dchi2)*qud11 + 
        (-(cdA222*chi) + 1.5*A22*dchi2)*qud21 + 
        (-(cdA223*chi) + 1.5*A23*dchi2)*qud31) + 
     ginv33*((-(cdA313*chi) + 1.5*A13*dchi3)*qud11 + 
        (-(cdA323*chi) + 1.5*A23*dchi3)*qud21 + 
        (-(cdA333*chi) + 1.5*A33*dchi3)*qud31) + 
     chi*((0.66666666666666666667*dK1 - dTheta1)*qud11 + 
        (0.66666666666666666667*dK2 - dTheta2)*qud21 + 
        (0.66666666666666666667*dK3 - dTheta3)*qud31) + 
     ginv12*((-(cdA212*chi) + 1.5*A12*dchi2)*qud21 + 
        (-(cdA213*chi) + 1.5*A13*dchi2)*qud31 - 
        chi*((cdA112 + cdA211)*qud11 + cdA122*qud21 + cdA123*qud31) + 
        1.5*((A12*dchi1 + A11*dchi2)*qud11 + dchi1*(A22*qud21 + A23*qud31))\
) + ginv13*((-(cdA312*chi) + 1.5*A12*dchi3)*qud21 + 
        (-(cdA313*chi) + 1.5*A13*dchi3)*qud31 - 
        chi*((cdA113 + cdA311)*qud11 + cdA123*qud21 + cdA133*qud31) + 
        1.5*((A13*dchi1 + A11*dchi3)*qud11 + dchi1*(A23*qud21 + A33*qud31))\
) + ginv23*((-(cdA322*chi) + 1.5*A22*dchi3)*qud21 + 
        (-(cdA323*chi) + 1.5*A23*dchi3)*qud31 - 
        chi*((cdA213 + cdA312)*qud11 + cdA223*qud21 + cdA233*qud31) + 
        1.5*((A13*dchi2 + A12*dchi3)*qud11 + dchi2*(A23*qud21 + A33*qud31))\
) + 0.5*(kappa1*((G1 - Gfromg1)*qdd11 + (G2 - Gfromg2)*qdd12 + 
           (G3 - Gfromg3)*qdd13) - dG13*qdd13*sup1 - dG21*qdd11*sup2 + 
        (dGfromgdu22*qdd12 - dG23*qdd13)*sup2 + 
        (dGfromgdu31*qdd11 + dGfromgdu32*qdd12 - dG33*qdd13)*sup3 + 
        qdd11*((-dG11 + dGfromgdu11)*sup1 + dGfromgdu21*sup2 - 
           dG31*sup3) + qdd12*
         ((-dG12 + dGfromgdu12)*sup1 - dG22*sup2 - dG32*sup3))*pow2(chi) + 
     sup1*(-2.*AA11*qud11 + 0.5*dGfromgdu13*qdd13*pow2(chi))) + 
  sup2*(chi*(-(cdda12*qud11) - cdda22*qud21 - cdda23*qud31 + 
        alpha*qud21*Rf22) + alpha*
      (chi*(qud11*Rf12 + qud31*Rf23) + 0.5*dGfromgdu23*qdd13*pow2(chi))) + 
  sup3*(chi*(-(cdda13*qud11) - cdda23*qud21 - cdda33*qud31 + 
        alpha*qud21*Rf23) + alpha*
      (chi*(qud11*Rf13 + qud31*Rf33) + 0.5*dGfromgdu33*qdd13*pow2(chi)))
;

rACsA2
=
(qud12*(lieA11 + alpha*chi*Rf11) + 
     qud22*(lieA12 + alpha*(-2.*AA12 + chi*Rf12)) + 
     qud32*(lieA13 + alpha*(-2.*AA13 + chi*Rf13)))*sup1 + 
  qud12*((-(cdda11*chi) + A11*alpha*K)*sup1 + lieA12*sup2 + 
     (A13*alpha*K + lieA13)*sup3 + 
     alpha*((-2.*AA21 + A12*K)*sup2 - 2.*AA31*sup3)) + 
  qud22*((-(cdda12*chi) + A12*alpha*K)*sup1 + lieA22*sup2 + 
     (A23*alpha*K + lieA23)*sup3 + 
     alpha*((-2.*AA22 + A22*K)*sup2 - 2.*AA32*sup3)) + 
  qud32*((-(cdda13*chi) + A13*alpha*K)*sup1 + lieA23*sup2 + 
     (A33*alpha*K + lieA33)*sup3 + 
     alpha*((-2.*AA23 + A23*K)*sup2 - 2.*AA33*sup3)) + 
  alpha*(ginv11*((-(cdA111*chi) + 1.5*A11*dchi1)*qud12 + 
        (-(cdA112*chi) + 1.5*A12*dchi1)*qud22 + 
        (-(cdA113*chi) + 1.5*A13*dchi1)*qud32) + 
     ginv22*((-(cdA212*chi) + 1.5*A12*dchi2)*qud12 + 
        (-(cdA222*chi) + 1.5*A22*dchi2)*qud22 + 
        (-(cdA223*chi) + 1.5*A23*dchi2)*qud32) + 
     ginv33*((-(cdA313*chi) + 1.5*A13*dchi3)*qud12 + 
        (-(cdA323*chi) + 1.5*A23*dchi3)*qud22 + 
        (-(cdA333*chi) + 1.5*A33*dchi3)*qud32) + 
     chi*((0.66666666666666666667*dK1 - dTheta1)*qud12 + 
        (0.66666666666666666667*dK2 - dTheta2)*qud22 + 
        (0.66666666666666666667*dK3 - dTheta3)*qud32) + 
     ginv12*((-(cdA212*chi) + 1.5*A12*dchi2)*qud22 + 
        (-(cdA213*chi) + 1.5*A13*dchi2)*qud32 - 
        chi*((cdA112 + cdA211)*qud12 + cdA122*qud22 + cdA123*qud32) + 
        1.5*((A12*dchi1 + A11*dchi2)*qud12 + dchi1*(A22*qud22 + A23*qud32))\
) + ginv13*((-(cdA312*chi) + 1.5*A12*dchi3)*qud22 + 
        (-(cdA313*chi) + 1.5*A13*dchi3)*qud32 - 
        chi*((cdA113 + cdA311)*qud12 + cdA123*qud22 + cdA133*qud32) + 
        1.5*((A13*dchi1 + A11*dchi3)*qud12 + dchi1*(A23*qud22 + A33*qud32))\
) + ginv23*((-(cdA322*chi) + 1.5*A22*dchi3)*qud22 + 
        (-(cdA323*chi) + 1.5*A23*dchi3)*qud32 - 
        chi*((cdA213 + cdA312)*qud12 + cdA223*qud22 + cdA233*qud32) + 
        1.5*((A13*dchi2 + A12*dchi3)*qud12 + dchi2*(A23*qud22 + A33*qud32))\
) + 0.5*(kappa1*((G1 - Gfromg1)*qdd12 + (G2 - Gfromg2)*qdd22 + 
           (G3 - Gfromg3)*qdd23) - dG13*qdd23*sup1 - dG21*qdd12*sup2 + 
        (dGfromgdu22*qdd22 - dG23*qdd23)*sup2 + 
        (dGfromgdu31*qdd12 + dGfromgdu32*qdd22 - dG33*qdd23)*sup3 + 
        qdd12*((-dG11 + dGfromgdu11)*sup1 + dGfromgdu21*sup2 - 
           dG31*sup3) + qdd22*
         ((-dG12 + dGfromgdu12)*sup1 - dG22*sup2 - dG32*sup3))*pow2(chi) + 
     sup1*(-2.*AA11*qud12 + 0.5*dGfromgdu13*qdd23*pow2(chi))) + 
  sup2*(chi*(-(cdda12*qud12) - cdda22*qud22 - cdda23*qud32 + 
        alpha*qud22*Rf22) + alpha*
      (chi*(qud12*Rf12 + qud32*Rf23) + 0.5*dGfromgdu23*qdd23*pow2(chi))) + 
  sup3*(chi*(-(cdda13*qud12) - cdda23*qud22 - cdda33*qud32 + 
        alpha*qud22*Rf23) + alpha*
      (chi*(qud12*Rf13 + qud32*Rf33) + 0.5*dGfromgdu33*qdd23*pow2(chi)))
;

rACsA3
=
(qud13*(lieA11 + alpha*chi*Rf11) + 
     qud23*(lieA12 + alpha*(-2.*AA12 + chi*Rf12)) + 
     qud33*(lieA13 + alpha*(-2.*AA13 + chi*Rf13)))*sup1 + 
  qud13*((-(cdda11*chi) + A11*alpha*K)*sup1 + lieA12*sup2 + 
     (A13*alpha*K + lieA13)*sup3 + 
     alpha*((-2.*AA21 + A12*K)*sup2 - 2.*AA31*sup3)) + 
  qud23*((-(cdda12*chi) + A12*alpha*K)*sup1 + lieA22*sup2 + 
     (A23*alpha*K + lieA23)*sup3 + 
     alpha*((-2.*AA22 + A22*K)*sup2 - 2.*AA32*sup3)) + 
  qud33*((-(cdda13*chi) + A13*alpha*K)*sup1 + lieA23*sup2 + 
     (A33*alpha*K + lieA33)*sup3 + 
     alpha*((-2.*AA23 + A23*K)*sup2 - 2.*AA33*sup3)) + 
  alpha*(ginv11*((-(cdA111*chi) + 1.5*A11*dchi1)*qud13 + 
        (-(cdA112*chi) + 1.5*A12*dchi1)*qud23 + 
        (-(cdA113*chi) + 1.5*A13*dchi1)*qud33) + 
     ginv22*((-(cdA212*chi) + 1.5*A12*dchi2)*qud13 + 
        (-(cdA222*chi) + 1.5*A22*dchi2)*qud23 + 
        (-(cdA223*chi) + 1.5*A23*dchi2)*qud33) + 
     ginv33*((-(cdA313*chi) + 1.5*A13*dchi3)*qud13 + 
        (-(cdA323*chi) + 1.5*A23*dchi3)*qud23 + 
        (-(cdA333*chi) + 1.5*A33*dchi3)*qud33) + 
     chi*((0.66666666666666666667*dK1 - dTheta1)*qud13 + 
        (0.66666666666666666667*dK2 - dTheta2)*qud23 + 
        (0.66666666666666666667*dK3 - dTheta3)*qud33) + 
     ginv12*((-(cdA212*chi) + 1.5*A12*dchi2)*qud23 + 
        (-(cdA213*chi) + 1.5*A13*dchi2)*qud33 - 
        chi*((cdA112 + cdA211)*qud13 + cdA122*qud23 + cdA123*qud33) + 
        1.5*((A12*dchi1 + A11*dchi2)*qud13 + dchi1*(A22*qud23 + A23*qud33))\
) + ginv13*((-(cdA312*chi) + 1.5*A12*dchi3)*qud23 + 
        (-(cdA313*chi) + 1.5*A13*dchi3)*qud33 - 
        chi*((cdA113 + cdA311)*qud13 + cdA123*qud23 + cdA133*qud33) + 
        1.5*((A13*dchi1 + A11*dchi3)*qud13 + dchi1*(A23*qud23 + A33*qud33))\
) + ginv23*((-(cdA322*chi) + 1.5*A22*dchi3)*qud23 + 
        (-(cdA323*chi) + 1.5*A23*dchi3)*qud33 - 
        chi*((cdA213 + cdA312)*qud13 + cdA223*qud23 + cdA233*qud33) + 
        1.5*((A13*dchi2 + A12*dchi3)*qud13 + dchi2*(A23*qud23 + A33*qud33))\
) + 0.5*(kappa1*((G1 - Gfromg1)*qdd13 + (G2 - Gfromg2)*qdd23 + 
           (G3 - Gfromg3)*qdd33) - dG13*qdd33*sup1 - dG21*qdd13*sup2 + 
        (dGfromgdu22*qdd23 - dG23*qdd33)*sup2 + 
        (dGfromgdu31*qdd13 + dGfromgdu32*qdd23 - dG33*qdd33)*sup3 + 
        qdd13*((-dG11 + dGfromgdu11)*sup1 + dGfromgdu21*sup2 - 
           dG31*sup3) + qdd23*
         ((-dG12 + dGfromgdu12)*sup1 - dG22*sup2 - dG32*sup3))*pow2(chi) + 
     sup1*(-2.*AA11*qud13 + 0.5*dGfromgdu13*qdd33*pow2(chi))) + 
  sup2*(chi*(-(cdda12*qud13) - cdda22*qud23 - cdda23*qud33 + 
        alpha*qud23*Rf22) + alpha*
      (chi*(qud13*Rf12 + qud33*Rf23) + 0.5*dGfromgdu23*qdd33*pow2(chi))) + 
  sup3*(chi*(-(cdda13*qud13) - cdda23*qud23 - cdda33*qud33 + 
        alpha*qud23*Rf23) + alpha*
      (chi*(qud13*Rf13 + qud33*Rf33) + 0.5*dGfromgdu33*qdd33*pow2(chi)))
;

rACABTF11
=
-(qPhysuudd1211*(2.*cdda12*chi + alpha*(AA21 + cdA112*sup1))) + 
  qPhysuudd3311*(-(cdda33*chi) + lieA33 + 
     alpha*(0.66666666666666666667*A33*K + cdA313*sup1 + cdA323*sup2)) + 
  qPhysuudd1111*(-(cdda11*chi) + lieA11 + 
     alpha*(-AA11 + 0.66666666666666666667*A11*K + cdA112*sup2 + 
        cdA113*sup3)) + qPhysuudd1211*
   (2.*lieA12 + alpha*(-AA12 + 1.3333333333333333333*A12*K + cdA211*sup1 + 
        cdA122*sup2 + cdA123*sup3)) + 
  qPhysuudd1311*(2.*(-(cdda13*chi) + lieA13) + 
     alpha*(-AA31 + 1.3333333333333333333*A13*K + cdA311*sup1 + 
        cdA123*sup2 + cdA133*sup3)) + 
  qPhysuudd2211*(-(cdda22*chi) + lieA22 + 
     alpha*(0.66666666666666666667*A22*K + cdA212*sup1 + cdA223*sup3)) + 
  qPhysuudd2311*(2.*(-(cdda23*chi) + lieA23) + 
     alpha*(-AA32 + 1.3333333333333333333*A23*K + cdA213*sup1 + 
        cdA322*sup2 + cdA233*sup3)) - 
  alpha*(AA13*qPhysuudd1311 + AA22*qPhysuudd2211 + AA23*qPhysuudd2311 + 
     AA33*qPhysuudd3311 + qPhysuudd1111*(cdA211*sup2 + cdA311*sup3)) + 
  alpha*(-((2.*cdA213*qPhysuudd1311 + 
          (0.5*(A12*dchi1*qPhysuudd1111 + A23*dchi3*qPhysuudd3311))/chi)*
        sup2) - qPhysuudd3311*((cdA133 + (0.5*A13*dchi3)/chi)*sup1 + 
        cdA233*sup2) - 2.*cdA312*qPhysuudd1211*sup3 + 
     qPhysuudd1211*((-cdA212 + (0.5*A12*dchi2)/chi)*sup2 + cdA213*sup3) + 
     qPhysuudd1311*((-cdA113 + (0.5*A13*dchi1)/chi)*sup1 + cdA312*sup2 - 
        cdA313*sup3) - qPhysuudd2211*
      ((cdA122 + (0.5*A12*dchi2)/chi)*sup1 + cdA322*sup3) + 
     qPhysuudd2311*((cdA312 + (A23*dchi1)/chi)*sup1 + 
        (0.5*A23*dchi2*sup2)/chi - cdA323*sup3) - 
     qPhysuudd2311*((2.*cdA123 + (0.5*A13*dchi2)/chi)*sup1 + cdA223*sup2 + 
        (0.5*A33*dchi2*sup3)/chi) + 
     ((-0.5*A22*dchi1*qPhysuudd1211 + A13*dchi2*qPhysuudd1311)*sup2 + 
        (A12*dchi3*qPhysuudd1211 - 
           0.5*dchi1*(A13*qPhysuudd1111 + A23*qPhysuudd1211))*sup3 + 
        0.5*(((A12*dchi1 - A11*dchi2)*qPhysuudd1211 - 
              dchi3*(A11*qPhysuudd1311 + A12*qPhysuudd2311) + 
              dchi1*(A22*qPhysuudd2211 + A33*qPhysuudd3311))*sup1 + 
           (-((A23*dchi1 + A12*dchi3)*qPhysuudd1311) - 
              A22*dchi3*qPhysuudd2311 + 
              dchi2*(A11*qPhysuudd1111 + A33*qPhysuudd3311))*sup2 + 
           (-(A33*dchi1*qPhysuudd1311) + 
              A13*(-(dchi2*qPhysuudd1211) + dchi3*qPhysuudd1311) + 
              dchi3*(A11*qPhysuudd1111 + A22*qPhysuudd2211) + 
              A23*(-(dchi2*qPhysuudd2211) + dchi3*qPhysuudd2311))*sup3))/chi)
;

rACABTF12
=
-(qPhysuudd1212*(2.*cdda12*chi + alpha*(AA21 + cdA112*sup1))) + 
  qPhysuudd3312*(-(cdda33*chi) + lieA33 + 
     alpha*(0.66666666666666666667*A33*K + cdA313*sup1 + cdA323*sup2)) + 
  qPhysuudd1112*(-(cdda11*chi) + lieA11 + 
     alpha*(-AA11 + 0.66666666666666666667*A11*K + cdA112*sup2 + 
        cdA113*sup3)) + qPhysuudd1212*
   (2.*lieA12 + alpha*(-AA12 + 1.3333333333333333333*A12*K + cdA211*sup1 + 
        cdA122*sup2 + cdA123*sup3)) + 
  qPhysuudd1312*(2.*(-(cdda13*chi) + lieA13) + 
     alpha*(-AA31 + 1.3333333333333333333*A13*K + cdA311*sup1 + 
        cdA123*sup2 + cdA133*sup3)) + 
  qPhysuudd2212*(-(cdda22*chi) + lieA22 + 
     alpha*(0.66666666666666666667*A22*K + cdA212*sup1 + cdA223*sup3)) + 
  qPhysuudd2312*(2.*(-(cdda23*chi) + lieA23) + 
     alpha*(-AA32 + 1.3333333333333333333*A23*K + cdA213*sup1 + 
        cdA322*sup2 + cdA233*sup3)) - 
  alpha*(AA13*qPhysuudd1312 + AA22*qPhysuudd2212 + AA23*qPhysuudd2312 + 
     AA33*qPhysuudd3312 + qPhysuudd1112*(cdA211*sup2 + cdA311*sup3)) + 
  alpha*(-((2.*cdA213*qPhysuudd1312 + 
          (0.5*(A12*dchi1*qPhysuudd1112 + A23*dchi3*qPhysuudd3312))/chi)*
        sup2) - qPhysuudd3312*((cdA133 + (0.5*A13*dchi3)/chi)*sup1 + 
        cdA233*sup2) - 2.*cdA312*qPhysuudd1212*sup3 + 
     qPhysuudd1212*((-cdA212 + (0.5*A12*dchi2)/chi)*sup2 + cdA213*sup3) + 
     qPhysuudd1312*((-cdA113 + (0.5*A13*dchi1)/chi)*sup1 + cdA312*sup2 - 
        cdA313*sup3) - qPhysuudd2212*
      ((cdA122 + (0.5*A12*dchi2)/chi)*sup1 + cdA322*sup3) + 
     qPhysuudd2312*((cdA312 + (A23*dchi1)/chi)*sup1 + 
        (0.5*A23*dchi2*sup2)/chi - cdA323*sup3) - 
     qPhysuudd2312*((2.*cdA123 + (0.5*A13*dchi2)/chi)*sup1 + cdA223*sup2 + 
        (0.5*A33*dchi2*sup3)/chi) + 
     ((-0.5*A22*dchi1*qPhysuudd1212 + A13*dchi2*qPhysuudd1312)*sup2 + 
        (A12*dchi3*qPhysuudd1212 - 
           0.5*dchi1*(A13*qPhysuudd1112 + A23*qPhysuudd1212))*sup3 + 
        0.5*(((A12*dchi1 - A11*dchi2)*qPhysuudd1212 - 
              dchi3*(A11*qPhysuudd1312 + A12*qPhysuudd2312) + 
              dchi1*(A22*qPhysuudd2212 + A33*qPhysuudd3312))*sup1 + 
           (-((A23*dchi1 + A12*dchi3)*qPhysuudd1312) - 
              A22*dchi3*qPhysuudd2312 + 
              dchi2*(A11*qPhysuudd1112 + A33*qPhysuudd3312))*sup2 + 
           (-(A33*dchi1*qPhysuudd1312) + 
              A13*(-(dchi2*qPhysuudd1212) + dchi3*qPhysuudd1312) + 
              dchi3*(A11*qPhysuudd1112 + A22*qPhysuudd2212) + 
              A23*(-(dchi2*qPhysuudd2212) + dchi3*qPhysuudd2312))*sup3))/chi)
;

rACABTF13
=
-(qPhysuudd1213*(2.*cdda12*chi + alpha*(AA21 + cdA112*sup1))) + 
  qPhysuudd3313*(-(cdda33*chi) + lieA33 + 
     alpha*(0.66666666666666666667*A33*K + cdA313*sup1 + cdA323*sup2)) + 
  qPhysuudd1113*(-(cdda11*chi) + lieA11 + 
     alpha*(-AA11 + 0.66666666666666666667*A11*K + cdA112*sup2 + 
        cdA113*sup3)) + qPhysuudd1213*
   (2.*lieA12 + alpha*(-AA12 + 1.3333333333333333333*A12*K + cdA211*sup1 + 
        cdA122*sup2 + cdA123*sup3)) + 
  qPhysuudd1313*(2.*(-(cdda13*chi) + lieA13) + 
     alpha*(-AA31 + 1.3333333333333333333*A13*K + cdA311*sup1 + 
        cdA123*sup2 + cdA133*sup3)) + 
  qPhysuudd2213*(-(cdda22*chi) + lieA22 + 
     alpha*(0.66666666666666666667*A22*K + cdA212*sup1 + cdA223*sup3)) + 
  qPhysuudd2313*(2.*(-(cdda23*chi) + lieA23) + 
     alpha*(-AA32 + 1.3333333333333333333*A23*K + cdA213*sup1 + 
        cdA322*sup2 + cdA233*sup3)) - 
  alpha*(AA13*qPhysuudd1313 + AA22*qPhysuudd2213 + AA23*qPhysuudd2313 + 
     AA33*qPhysuudd3313 + qPhysuudd1113*(cdA211*sup2 + cdA311*sup3)) + 
  alpha*(-((2.*cdA213*qPhysuudd1313 + 
          (0.5*(A12*dchi1*qPhysuudd1113 + A23*dchi3*qPhysuudd3313))/chi)*
        sup2) - qPhysuudd3313*((cdA133 + (0.5*A13*dchi3)/chi)*sup1 + 
        cdA233*sup2) - 2.*cdA312*qPhysuudd1213*sup3 + 
     qPhysuudd1213*((-cdA212 + (0.5*A12*dchi2)/chi)*sup2 + cdA213*sup3) + 
     qPhysuudd1313*((-cdA113 + (0.5*A13*dchi1)/chi)*sup1 + cdA312*sup2 - 
        cdA313*sup3) - qPhysuudd2213*
      ((cdA122 + (0.5*A12*dchi2)/chi)*sup1 + cdA322*sup3) + 
     qPhysuudd2313*((cdA312 + (A23*dchi1)/chi)*sup1 + 
        (0.5*A23*dchi2*sup2)/chi - cdA323*sup3) - 
     qPhysuudd2313*((2.*cdA123 + (0.5*A13*dchi2)/chi)*sup1 + cdA223*sup2 + 
        (0.5*A33*dchi2*sup3)/chi) + 
     ((-0.5*A22*dchi1*qPhysuudd1213 + A13*dchi2*qPhysuudd1313)*sup2 + 
        (A12*dchi3*qPhysuudd1213 - 
           0.5*dchi1*(A13*qPhysuudd1113 + A23*qPhysuudd1213))*sup3 + 
        0.5*(((A12*dchi1 - A11*dchi2)*qPhysuudd1213 - 
              dchi3*(A11*qPhysuudd1313 + A12*qPhysuudd2313) + 
              dchi1*(A22*qPhysuudd2213 + A33*qPhysuudd3313))*sup1 + 
           (-((A23*dchi1 + A12*dchi3)*qPhysuudd1313) - 
              A22*dchi3*qPhysuudd2313 + 
              dchi2*(A11*qPhysuudd1113 + A33*qPhysuudd3313))*sup2 + 
           (-(A33*dchi1*qPhysuudd1313) + 
              A13*(-(dchi2*qPhysuudd1213) + dchi3*qPhysuudd1313) + 
              dchi3*(A11*qPhysuudd1113 + A22*qPhysuudd2213) + 
              A23*(-(dchi2*qPhysuudd2213) + dchi3*qPhysuudd2313))*sup3))/chi)
;

rACABTF22
=
-(qPhysuudd1222*(2.*cdda12*chi + alpha*(AA21 + cdA112*sup1))) + 
  qPhysuudd3322*(-(cdda33*chi) + lieA33 + 
     alpha*(0.66666666666666666667*A33*K + cdA313*sup1 + cdA323*sup2)) + 
  qPhysuudd1122*(-(cdda11*chi) + lieA11 + 
     alpha*(-AA11 + 0.66666666666666666667*A11*K + cdA112*sup2 + 
        cdA113*sup3)) + qPhysuudd1222*
   (2.*lieA12 + alpha*(-AA12 + 1.3333333333333333333*A12*K + cdA211*sup1 + 
        cdA122*sup2 + cdA123*sup3)) + 
  qPhysuudd1322*(2.*(-(cdda13*chi) + lieA13) + 
     alpha*(-AA31 + 1.3333333333333333333*A13*K + cdA311*sup1 + 
        cdA123*sup2 + cdA133*sup3)) + 
  qPhysuudd2222*(-(cdda22*chi) + lieA22 + 
     alpha*(0.66666666666666666667*A22*K + cdA212*sup1 + cdA223*sup3)) + 
  qPhysuudd2322*(2.*(-(cdda23*chi) + lieA23) + 
     alpha*(-AA32 + 1.3333333333333333333*A23*K + cdA213*sup1 + 
        cdA322*sup2 + cdA233*sup3)) - 
  alpha*(AA13*qPhysuudd1322 + AA22*qPhysuudd2222 + AA23*qPhysuudd2322 + 
     AA33*qPhysuudd3322 + qPhysuudd1122*(cdA211*sup2 + cdA311*sup3)) + 
  alpha*(-((2.*cdA213*qPhysuudd1322 + 
          (0.5*(A12*dchi1*qPhysuudd1122 + A23*dchi3*qPhysuudd3322))/chi)*
        sup2) - qPhysuudd3322*((cdA133 + (0.5*A13*dchi3)/chi)*sup1 + 
        cdA233*sup2) - 2.*cdA312*qPhysuudd1222*sup3 + 
     qPhysuudd1222*((-cdA212 + (0.5*A12*dchi2)/chi)*sup2 + cdA213*sup3) + 
     qPhysuudd1322*((-cdA113 + (0.5*A13*dchi1)/chi)*sup1 + cdA312*sup2 - 
        cdA313*sup3) - qPhysuudd2222*
      ((cdA122 + (0.5*A12*dchi2)/chi)*sup1 + cdA322*sup3) + 
     qPhysuudd2322*((cdA312 + (A23*dchi1)/chi)*sup1 + 
        (0.5*A23*dchi2*sup2)/chi - cdA323*sup3) - 
     qPhysuudd2322*((2.*cdA123 + (0.5*A13*dchi2)/chi)*sup1 + cdA223*sup2 + 
        (0.5*A33*dchi2*sup3)/chi) + 
     ((-0.5*A22*dchi1*qPhysuudd1222 + A13*dchi2*qPhysuudd1322)*sup2 + 
        (A12*dchi3*qPhysuudd1222 - 
           0.5*dchi1*(A13*qPhysuudd1122 + A23*qPhysuudd1222))*sup3 + 
        0.5*(((A12*dchi1 - A11*dchi2)*qPhysuudd1222 - 
              dchi3*(A11*qPhysuudd1322 + A12*qPhysuudd2322) + 
              dchi1*(A22*qPhysuudd2222 + A33*qPhysuudd3322))*sup1 + 
           (-((A23*dchi1 + A12*dchi3)*qPhysuudd1322) - 
              A22*dchi3*qPhysuudd2322 + 
              dchi2*(A11*qPhysuudd1122 + A33*qPhysuudd3322))*sup2 + 
           (-(A33*dchi1*qPhysuudd1322) + 
              A13*(-(dchi2*qPhysuudd1222) + dchi3*qPhysuudd1322) + 
              dchi3*(A11*qPhysuudd1122 + A22*qPhysuudd2222) + 
              A23*(-(dchi2*qPhysuudd2222) + dchi3*qPhysuudd2322))*sup3))/chi)
;

rACABTF23
=
-(qPhysuudd1223*(2.*cdda12*chi + alpha*(AA21 + cdA112*sup1))) + 
  qPhysuudd3323*(-(cdda33*chi) + lieA33 + 
     alpha*(0.66666666666666666667*A33*K + cdA313*sup1 + cdA323*sup2)) + 
  qPhysuudd1123*(-(cdda11*chi) + lieA11 + 
     alpha*(-AA11 + 0.66666666666666666667*A11*K + cdA112*sup2 + 
        cdA113*sup3)) + qPhysuudd1223*
   (2.*lieA12 + alpha*(-AA12 + 1.3333333333333333333*A12*K + cdA211*sup1 + 
        cdA122*sup2 + cdA123*sup3)) + 
  qPhysuudd1323*(2.*(-(cdda13*chi) + lieA13) + 
     alpha*(-AA31 + 1.3333333333333333333*A13*K + cdA311*sup1 + 
        cdA123*sup2 + cdA133*sup3)) + 
  qPhysuudd2223*(-(cdda22*chi) + lieA22 + 
     alpha*(0.66666666666666666667*A22*K + cdA212*sup1 + cdA223*sup3)) + 
  qPhysuudd2323*(2.*(-(cdda23*chi) + lieA23) + 
     alpha*(-AA32 + 1.3333333333333333333*A23*K + cdA213*sup1 + 
        cdA322*sup2 + cdA233*sup3)) - 
  alpha*(AA13*qPhysuudd1323 + AA22*qPhysuudd2223 + AA23*qPhysuudd2323 + 
     AA33*qPhysuudd3323 + qPhysuudd1123*(cdA211*sup2 + cdA311*sup3)) + 
  alpha*(-((2.*cdA213*qPhysuudd1323 + 
          (0.5*(A12*dchi1*qPhysuudd1123 + A23*dchi3*qPhysuudd3323))/chi)*
        sup2) - qPhysuudd3323*((cdA133 + (0.5*A13*dchi3)/chi)*sup1 + 
        cdA233*sup2) - 2.*cdA312*qPhysuudd1223*sup3 + 
     qPhysuudd1223*((-cdA212 + (0.5*A12*dchi2)/chi)*sup2 + cdA213*sup3) + 
     qPhysuudd1323*((-cdA113 + (0.5*A13*dchi1)/chi)*sup1 + cdA312*sup2 - 
        cdA313*sup3) - qPhysuudd2223*
      ((cdA122 + (0.5*A12*dchi2)/chi)*sup1 + cdA322*sup3) + 
     qPhysuudd2323*((cdA312 + (A23*dchi1)/chi)*sup1 + 
        (0.5*A23*dchi2*sup2)/chi - cdA323*sup3) - 
     qPhysuudd2323*((2.*cdA123 + (0.5*A13*dchi2)/chi)*sup1 + cdA223*sup2 + 
        (0.5*A33*dchi2*sup3)/chi) + 
     ((-0.5*A22*dchi1*qPhysuudd1223 + A13*dchi2*qPhysuudd1323)*sup2 + 
        (A12*dchi3*qPhysuudd1223 - 
           0.5*dchi1*(A13*qPhysuudd1123 + A23*qPhysuudd1223))*sup3 + 
        0.5*(((A12*dchi1 - A11*dchi2)*qPhysuudd1223 - 
              dchi3*(A11*qPhysuudd1323 + A12*qPhysuudd2323) + 
              dchi1*(A22*qPhysuudd2223 + A33*qPhysuudd3323))*sup1 + 
           (-((A23*dchi1 + A12*dchi3)*qPhysuudd1323) - 
              A22*dchi3*qPhysuudd2323 + 
              dchi2*(A11*qPhysuudd1123 + A33*qPhysuudd3323))*sup2 + 
           (-(A33*dchi1*qPhysuudd1323) + 
              A13*(-(dchi2*qPhysuudd1223) + dchi3*qPhysuudd1323) + 
              dchi3*(A11*qPhysuudd1123 + A22*qPhysuudd2223) + 
              A23*(-(dchi2*qPhysuudd2223) + dchi3*qPhysuudd2323))*sup3))/chi)
;

rACABTF33
=
-(qPhysuudd1233*(2.*cdda12*chi + alpha*(AA21 + cdA112*sup1))) + 
  qPhysuudd3333*(-(cdda33*chi) + lieA33 + 
     alpha*(0.66666666666666666667*A33*K + cdA313*sup1 + cdA323*sup2)) + 
  qPhysuudd1133*(-(cdda11*chi) + lieA11 + 
     alpha*(-AA11 + 0.66666666666666666667*A11*K + cdA112*sup2 + 
        cdA113*sup3)) + qPhysuudd1233*
   (2.*lieA12 + alpha*(-AA12 + 1.3333333333333333333*A12*K + cdA211*sup1 + 
        cdA122*sup2 + cdA123*sup3)) + 
  qPhysuudd1333*(2.*(-(cdda13*chi) + lieA13) + 
     alpha*(-AA31 + 1.3333333333333333333*A13*K + cdA311*sup1 + 
        cdA123*sup2 + cdA133*sup3)) + 
  qPhysuudd2233*(-(cdda22*chi) + lieA22 + 
     alpha*(0.66666666666666666667*A22*K + cdA212*sup1 + cdA223*sup3)) + 
  qPhysuudd2333*(2.*(-(cdda23*chi) + lieA23) + 
     alpha*(-AA32 + 1.3333333333333333333*A23*K + cdA213*sup1 + 
        cdA322*sup2 + cdA233*sup3)) - 
  alpha*(AA13*qPhysuudd1333 + AA22*qPhysuudd2233 + AA23*qPhysuudd2333 + 
     AA33*qPhysuudd3333 + qPhysuudd1133*(cdA211*sup2 + cdA311*sup3)) + 
  alpha*(-((2.*cdA213*qPhysuudd1333 + 
          (0.5*(A12*dchi1*qPhysuudd1133 + A23*dchi3*qPhysuudd3333))/chi)*
        sup2) - qPhysuudd3333*((cdA133 + (0.5*A13*dchi3)/chi)*sup1 + 
        cdA233*sup2) - 2.*cdA312*qPhysuudd1233*sup3 + 
     qPhysuudd1233*((-cdA212 + (0.5*A12*dchi2)/chi)*sup2 + cdA213*sup3) + 
     qPhysuudd1333*((-cdA113 + (0.5*A13*dchi1)/chi)*sup1 + cdA312*sup2 - 
        cdA313*sup3) - qPhysuudd2233*
      ((cdA122 + (0.5*A12*dchi2)/chi)*sup1 + cdA322*sup3) + 
     qPhysuudd2333*((cdA312 + (A23*dchi1)/chi)*sup1 + 
        (0.5*A23*dchi2*sup2)/chi - cdA323*sup3) - 
     qPhysuudd2333*((2.*cdA123 + (0.5*A13*dchi2)/chi)*sup1 + cdA223*sup2 + 
        (0.5*A33*dchi2*sup3)/chi) + 
     ((-0.5*A22*dchi1*qPhysuudd1233 + A13*dchi2*qPhysuudd1333)*sup2 + 
        (A12*dchi3*qPhysuudd1233 - 
           0.5*dchi1*(A13*qPhysuudd1133 + A23*qPhysuudd1233))*sup3 + 
        0.5*(((A12*dchi1 - A11*dchi2)*qPhysuudd1233 - 
              dchi3*(A11*qPhysuudd1333 + A12*qPhysuudd2333) + 
              dchi1*(A22*qPhysuudd2233 + A33*qPhysuudd3333))*sup1 + 
           (-((A23*dchi1 + A12*dchi3)*qPhysuudd1333) - 
              A22*dchi3*qPhysuudd2333 + 
              dchi2*(A11*qPhysuudd1133 + A33*qPhysuudd3333))*sup2 + 
           (-(A33*dchi1*qPhysuudd1333) + 
              A13*(-(dchi2*qPhysuudd1233) + dchi3*qPhysuudd1333) + 
              dchi3*(A11*qPhysuudd1133 + A22*qPhysuudd2233) + 
              A23*(-(dchi2*qPhysuudd2233) + dchi3*qPhysuudd2333))*sup3))/chi)
;


if (givehPsi0) { 

gADM11
=
g11/chi
;

gADM12
=
g12/chi
;

gADM13
=
g13/chi
;

gADM21
=
g12/chi
;

gADM22
=
g22/chi
;

gADM23
=
g23/chi
;

gADM31
=
g13/chi
;

gADM32
=
g23/chi
;

gADM33
=
g33/chi
;

vu1
=
-yp
;

vu2
=
xp
;

vu3
=
0
;

wu1
=
((-(ADMginv13*sup2) + ADMginv12*sup3)*vu1 + 
    (ADMginv13*sup1 - ADMginv11*sup3)*vu2 + 
    (-(ADMginv12*sup1) + ADMginv11*sup2)*vu3)/Power(chi,1.5)
;

wu2
=
((-(ADMginv23*sup2) + ADMginv22*sup3)*vu1 + 
    (ADMginv23*sup1 - ADMginv12*sup3)*vu2 + 
    (-(ADMginv22*sup1) + ADMginv12*sup2)*vu3)/Power(chi,1.5)
;

wu3
=
((-(ADMginv33*sup2) + ADMginv23*sup3)*vu1 + 
    (ADMginv33*sup1 - ADMginv13*sup3)*vu2 + 
    (-(ADMginv23*sup1) + ADMginv13*sup2)*vu3)/Power(chi,1.5)
;

sdotv
=
(gADM11*sup1 + gADM21*sup2 + gADM31*sup3)*vu1 + 
  (gADM12*sup1 + gADM22*sup2 + gADM32*sup3)*vu2 + 
  (gADM13*sup1 + gADM23*sup2 + gADM33*sup3)*vu3
;

vu1
=
-(sdotv*sup1) + vu1
;

vu2
=
-(sdotv*sup2) + vu2
;

vu3
=
-(sdotv*sup3) + vu3
;

vdotv
=
(gADM31*vu1 + (gADM23 + gADM32)*vu2)*vu3 + 
  vu1*((gADM12 + gADM21)*vu2 + gADM13*vu3) + gADM11*pow2(vu1) + 
  gADM22*pow2(vu2) + gADM33*pow2(vu3)
;

vu1
=
vu1/Sqrt(vdotv)
;

vu2
=
vu2/Sqrt(vdotv)
;

vu3
=
vu3/Sqrt(vdotv)
;

sdotw
=
(gADM11*sup1 + gADM21*sup2 + gADM31*sup3)*wu1 + 
  (gADM12*sup1 + gADM22*sup2 + gADM32*sup3)*wu2 + 
  (gADM13*sup1 + gADM23*sup2 + gADM33*sup3)*wu3
;

vdotw
=
(gADM11*vu1 + gADM21*vu2 + gADM31*vu3)*wu1 + 
  (gADM12*vu1 + gADM22*vu2 + gADM32*vu3)*wu2 + 
  (gADM13*vu1 + gADM23*vu2 + gADM33*vu3)*wu3
;

wu1
=
-(sdotw*sup1) - vdotw*vu1 + wu1
;

wu2
=
-(sdotw*sup2) - vdotw*vu2 + wu2
;

wu3
=
-(sdotw*sup3) - vdotw*vu3 + wu3
;

wdotw
=
(gADM31*wu1 + (gADM23 + gADM32)*wu2)*wu3 + 
  wu1*((gADM12 + gADM21)*wu2 + gADM13*wu3) + gADM11*pow2(wu1) + 
  gADM22*pow2(wu2) + gADM33*pow2(wu3)
;

wu1
=
wu1/Sqrt(wdotw)
;

wu2
=
wu2/Sqrt(wdotw)
;

wu3
=
wu3/Sqrt(wdotw)
;

vd1
=
gADM11*vu1 + gADM12*vu2 + gADM13*vu3
;

vd2
=
gADM21*vu1 + gADM22*vu2 + gADM23*vu3
;

vd3
=
gADM31*vu1 + gADM32*vu2 + gADM33*vu3
;

wd1
=
gADM11*wu1 + gADM12*wu2 + gADM13*wu3
;

wd2
=
gADM21*wu1 + gADM22*wu2 + gADM23*wu3
;

wd3
=
gADM31*wu1 + gADM32*wu2 + gADM33*wu3
;

RehPsi0
=
Power(2.7182818284590452354,pow2(hPsi0parb)*
    (2.*hPsi0parc*time - pow2(hPsi0parc) - pow2(time)))*hPsi0para
;

ImhPsi0
=
0
;

rACABTF11
=
rACABTF11 + alpha*chi*(2.*ImhPsi0*vd1*wd1 + RehPsi0*(pow2(vd1) - pow2(wd1)))
;

rACABTF12
=
rACABTF12 + alpha*chi*(vd2*(RehPsi0*vd1 + ImhPsi0*wd1) + 
     (ImhPsi0*vd1 - RehPsi0*wd1)*wd2)
;

rACABTF13
=
rACABTF13 + alpha*chi*(vd3*(RehPsi0*vd1 + ImhPsi0*wd1) + 
     (ImhPsi0*vd1 - RehPsi0*wd1)*wd3)
;

rACABTF22
=
rACABTF22 + alpha*chi*(2.*ImhPsi0*vd2*wd2 + RehPsi0*(pow2(vd2) - pow2(wd2)))
;

rACABTF23
=
rACABTF23 + alpha*chi*(vd3*(RehPsi0*vd2 + ImhPsi0*wd2) + 
     (ImhPsi0*vd2 - RehPsi0*wd2)*wd3)
;

rACABTF33
=
rACABTF33 + alpha*chi*(2.*ImhPsi0*vd3*wd3 + RehPsi0*(pow2(vd3) - pow2(wd3)))
;


 }  

rA11
=
rACABTF11 + 0.5*qdd11*rACqq + 2.*
   (qud11*rACsA1 + qud21*rACsA2 + qud31*rACsA3)*sdown1 + rACss*pow2(sdown1)
;

rA12
=
rACABTF12 + 0.5*qdd12*rACqq + (qud11*rACsA1 + qud21*rACsA2 + qud31*rACsA3)*
   sdown2 + sdown1*(qud12*rACsA1 + qud22*rACsA2 + qud32*rACsA3 + 
     rACss*sdown2)
;

rA13
=
rACABTF13 + 0.5*qdd13*rACqq + (qud11*rACsA1 + qud21*rACsA2 + qud31*rACsA3)*
   sdown3 + sdown1*(qud13*rACsA1 + qud23*rACsA2 + qud33*rACsA3 + 
     rACss*sdown3)
;

rA22
=
rACABTF22 + 0.5*qdd22*rACqq + 2.*
   (qud12*rACsA1 + qud22*rACsA2 + qud32*rACsA3)*sdown2 + rACss*pow2(sdown2)
;

rA23
=
rACABTF23 + 0.5*qdd23*rACqq + (qud12*rACsA1 + qud22*rACsA2 + qud32*rACsA3)*
   sdown3 + sdown2*(qud13*rACsA1 + qud23*rACsA2 + qud33*rACsA3 + 
     rACss*sdown3)
;

rA33
=
rACABTF33 + 0.5*qdd33*rACqq + 2.*
   (qud13*rACsA1 + qud23*rACsA2 + qud33*rACsA3)*sdown3 + rACss*pow2(sdown3)
;

rG1
=
qud11*rGamA1 + qud12*rGamA2 + qud13*rGamA3 + rGams*sup1
;

rG2
=
qud21*rGamA1 + qud22*rGamA2 + qud23*rGamA3 + rGams*sup2
;

rG3
=
qud31*rGamA1 + qud32*rGamA2 + qud33*rGamA3 + rGams*sup3
;
#else
// code adapted from David 2012-8-18

detginv
=
1/(2.*g12*g13*g23 + g11*g22*g33 - 
    g33*pow2(g12) - g22*pow2(g13) - 
    g11*pow2(g23))
;

ginv11
=
detginv*(g22*g33 - pow2(g23))
;

ginv12
=
detginv*(g13*g23 - g12*g33)
;

ginv13
=
detginv*(-(g13*g22) + g12*g23)
;

ginv22
=
detginv*(g11*g33 - pow2(g13))
;

ginv23
=
detginv*(g12*g13 - g11*g23)
;

ginv33
=
detginv*(g11*g22 - pow2(g12))
;

ADMginv11
=
ginv11*chi
;

ADMginv12
=
ginv12*chi
;

ADMginv13
=
ginv13*chi
;

ADMginv22
=
ginv22*chi
;

ADMginv23
=
ginv23*chi
;

ADMginv33
=
ginv33*chi
;

modshatARG
=
2.*(ADMginv23*shat2*shat3 + shat1*(ADMginv12*shat2 + ADMginv13*shat3)) + 
  ADMginv11*pow2(shat1) + ADMginv22*pow2(shat2) + ADMginv33*pow2(shat3)
;


if (modshatARG<0.00001) {                           
      printf("modshat is wrong (%e)\n",modshatARG);
      modshatARG = 0.00001;
    }oomodshat
=
1/sqrt(modshatARG)
;

sdown1
=
oomodshat*shat1
;

sdown2
=
oomodshat*shat2
;

sdown3
=
oomodshat*shat3
;

sup1
=
ADMginv11*sdown1 + ADMginv12*sdown2 + ADMginv13*sdown3
;

sup2
=
ADMginv12*sdown1 + ADMginv22*sdown2 + ADMginv23*sdown3
;

sup3
=
ADMginv13*sdown1 + ADMginv23*sdown2 + ADMginv33*sdown3
;

qud11
=
1. - sdown1*sup1
;

qud12
=
-(sdown2*sup1)
;

qud13
=
-(sdown3*sup1)
;

qud21
=
-(sdown1*sup2)
;

qud22
=
1. - sdown2*sup2
;

qud23
=
-(sdown3*sup2)
;

qud31
=
-(sdown1*sup3)
;

qud32
=
-(sdown2*sup3)
;

qud33
=
1. - sdown3*sup3
;

qdd11
=
g11/chi - pow2(sdown1)
;

qdd12
=
-(sdown1*sdown2) + g12/chi
;

qdd13
=
-(sdown1*sdown3) + g13/chi
;

qdd22
=
g22/chi - pow2(sdown2)
;

qdd23
=
-(sdown2*sdown3) + g23/chi
;

qdd33
=
g33/chi - pow2(sdown3)
;

quu11
=
ADMginv11 - pow2(sup1)
;

quu12
=
ADMginv12 - sup1*sup2
;

quu13
=
ADMginv13 - sup1*sup3
;

quu22
=
ADMginv22 - pow2(sup2)
;

quu23
=
ADMginv23 - sup2*sup3
;

quu33
=
ADMginv33 - pow2(sup3)
;

qPhysuudd1111
=
-0.5*qdd11*quu11 + pow2(qud11)
;

qPhysuudd1112
=
qud11*qud12 - 0.5*qdd12*quu11
;

qPhysuudd1113
=
qud11*qud13 - 0.5*qdd13*quu11
;

qPhysuudd1122
=
-0.5*qdd22*quu11 + pow2(qud12)
;

qPhysuudd1123
=
qud12*qud13 - 0.5*qdd23*quu11
;

qPhysuudd1133
=
-0.5*qdd33*quu11 + pow2(qud13)
;

qPhysuudd1211
=
qud11*qud21 - 0.5*qdd11*quu12
;

qPhysuudd1212
=
0.5*(qud12*qud21 + qud11*qud22 - qdd12*quu12)
;

qPhysuudd1213
=
0.5*(qud13*qud21 + qud11*qud23 - qdd13*quu12)
;

qPhysuudd1222
=
qud12*qud22 - 0.5*qdd22*quu12
;

qPhysuudd1223
=
0.5*(qud13*qud22 + qud12*qud23 - qdd23*quu12)
;

qPhysuudd1233
=
qud13*qud23 - 0.5*qdd33*quu12
;

qPhysuudd1311
=
qud11*qud31 - 0.5*qdd11*quu13
;

qPhysuudd1312
=
0.5*(qud12*qud31 + qud11*qud32 - qdd12*quu13)
;

qPhysuudd1313
=
0.5*(qud13*qud31 + qud11*qud33 - qdd13*quu13)
;

qPhysuudd1322
=
qud12*qud32 - 0.5*qdd22*quu13
;

qPhysuudd1323
=
0.5*(qud13*qud32 + qud12*qud33 - qdd23*quu13)
;

qPhysuudd1333
=
qud13*qud33 - 0.5*qdd33*quu13
;

qPhysuudd2211
=
-0.5*qdd11*quu22 + pow2(qud21)
;

qPhysuudd2212
=
qud21*qud22 - 0.5*qdd12*quu22
;

qPhysuudd2213
=
qud21*qud23 - 0.5*qdd13*quu22
;

qPhysuudd2222
=
-0.5*qdd22*quu22 + pow2(qud22)
;

qPhysuudd2223
=
qud22*qud23 - 0.5*qdd23*quu22
;

qPhysuudd2233
=
-0.5*qdd33*quu22 + pow2(qud23)
;

qPhysuudd2311
=
qud21*qud31 - 0.5*qdd11*quu23
;

qPhysuudd2312
=
0.5*(qud22*qud31 + qud21*qud32 - qdd12*quu23)
;

qPhysuudd2313
=
0.5*(qud23*qud31 + qud21*qud33 - qdd13*quu23)
;

qPhysuudd2322
=
qud22*qud32 - 0.5*qdd22*quu23
;

qPhysuudd2323
=
0.5*(qud23*qud32 + qud22*qud33 - qdd23*quu23)
;

qPhysuudd2333
=
qud23*qud33 - 0.5*qdd33*quu23
;

qPhysuudd3311
=
-0.5*qdd11*quu33 + pow2(qud31)
;

qPhysuudd3312
=
qud31*qud32 - 0.5*qdd12*quu33
;

qPhysuudd3313
=
qud31*qud33 - 0.5*qdd13*quu33
;

qPhysuudd3322
=
-0.5*qdd22*quu33 + pow2(qud32)
;

qPhysuudd3323
=
qud32*qud33 - 0.5*qdd23*quu33
;

qPhysuudd3333
=
-0.5*qdd33*quu33 + pow2(qud33)
;

muL
=
2./alpha
;

muStilde
=
1/chi
;

vbetas
=
2.*sqrt(0.33333333333333333333*muStilde)
;

vbetaA
=
sqrt(muStilde)
;

K
=
Khat + 2.*Theta
;

dK1
=
dKhat1 + 2.*dTheta1
;

dK2
=
dKhat2 + 2.*dTheta2
;

dK3
=
dKhat3 + 2.*dTheta3
;

dginv111
=
-2.*(dg123*ginv12*ginv13 + ginv11*(dg112*ginv12 + dg113*ginv13)) - 
  dg111*pow2(ginv11) - dg122*pow2(ginv12) - dg133*pow2(ginv13)
;

dginv112
=
-(ginv11*(dg111*ginv12 + dg112*ginv22 + dg113*ginv23)) - 
  ginv12*(dg113*ginv13 + dg122*ginv22 + dg123*ginv23) - 
  ginv13*(dg123*ginv22 + dg133*ginv23) - dg112*pow2(ginv12)
;

dginv113
=
-(ginv11*(dg111*ginv13 + dg112*ginv23 + dg113*ginv33)) - 
  ginv12*(dg112*ginv13 + dg122*ginv23 + dg123*ginv33) - 
  ginv13*(dg123*ginv23 + dg133*ginv33) - dg113*pow2(ginv13)
;

dginv122
=
-2.*(dg123*ginv22*ginv23 + ginv12*(dg112*ginv22 + dg113*ginv23)) - 
  dg111*pow2(ginv12) - dg122*pow2(ginv22) - dg133*pow2(ginv23)
;

dginv123
=
-(ginv13*(dg112*ginv22 + dg113*ginv23)) - dg133*ginv23*ginv33 - 
  ginv12*(dg111*ginv13 + dg112*ginv23 + dg113*ginv33) - 
  ginv22*(dg122*ginv23 + dg123*ginv33) - dg123*pow2(ginv23)
;

dginv133
=
-2.*(dg123*ginv23*ginv33 + ginv13*(dg112*ginv23 + dg113*ginv33)) - 
  dg111*pow2(ginv13) - dg122*pow2(ginv23) - dg133*pow2(ginv33)
;

dginv211
=
-2.*(dg223*ginv12*ginv13 + ginv11*(dg212*ginv12 + dg213*ginv13)) - 
  dg211*pow2(ginv11) - dg222*pow2(ginv12) - dg233*pow2(ginv13)
;

dginv212
=
-(ginv11*(dg211*ginv12 + dg212*ginv22 + dg213*ginv23)) - 
  ginv12*(dg213*ginv13 + dg222*ginv22 + dg223*ginv23) - 
  ginv13*(dg223*ginv22 + dg233*ginv23) - dg212*pow2(ginv12)
;

dginv213
=
-(ginv11*(dg211*ginv13 + dg212*ginv23 + dg213*ginv33)) - 
  ginv12*(dg212*ginv13 + dg222*ginv23 + dg223*ginv33) - 
  ginv13*(dg223*ginv23 + dg233*ginv33) - dg213*pow2(ginv13)
;

dginv222
=
-2.*(dg223*ginv22*ginv23 + ginv12*(dg212*ginv22 + dg213*ginv23)) - 
  dg211*pow2(ginv12) - dg222*pow2(ginv22) - dg233*pow2(ginv23)
;

dginv223
=
-(ginv13*(dg212*ginv22 + dg213*ginv23)) - dg233*ginv23*ginv33 - 
  ginv12*(dg211*ginv13 + dg212*ginv23 + dg213*ginv33) - 
  ginv22*(dg222*ginv23 + dg223*ginv33) - dg223*pow2(ginv23)
;

dginv233
=
-2.*(dg223*ginv23*ginv33 + ginv13*(dg212*ginv23 + dg213*ginv33)) - 
  dg211*pow2(ginv13) - dg222*pow2(ginv23) - dg233*pow2(ginv33)
;

dginv311
=
-2.*(dg323*ginv12*ginv13 + ginv11*(dg312*ginv12 + dg313*ginv13)) - 
  dg311*pow2(ginv11) - dg322*pow2(ginv12) - dg333*pow2(ginv13)
;

dginv312
=
-(ginv11*(dg311*ginv12 + dg312*ginv22 + dg313*ginv23)) - 
  ginv12*(dg313*ginv13 + dg322*ginv22 + dg323*ginv23) - 
  ginv13*(dg323*ginv22 + dg333*ginv23) - dg312*pow2(ginv12)
;

dginv313
=
-(ginv11*(dg311*ginv13 + dg312*ginv23 + dg313*ginv33)) - 
  ginv12*(dg312*ginv13 + dg322*ginv23 + dg323*ginv33) - 
  ginv13*(dg323*ginv23 + dg333*ginv33) - dg313*pow2(ginv13)
;

dginv322
=
-2.*(dg323*ginv22*ginv23 + ginv12*(dg312*ginv22 + dg313*ginv23)) - 
  dg311*pow2(ginv12) - dg322*pow2(ginv22) - dg333*pow2(ginv23)
;

dginv323
=
-(ginv13*(dg312*ginv22 + dg313*ginv23)) - dg333*ginv23*ginv33 - 
  ginv12*(dg311*ginv13 + dg312*ginv23 + dg313*ginv33) - 
  ginv22*(dg322*ginv23 + dg323*ginv33) - dg323*pow2(ginv23)
;

dginv333
=
-2.*(dg323*ginv23*ginv33 + ginv13*(dg312*ginv23 + dg313*ginv33)) - 
  dg311*pow2(ginv13) - dg322*pow2(ginv23) - dg333*pow2(ginv33)
;

gammado111
=
0.5*dg111
;

gammado112
=
0.5*dg211
;

gammado113
=
0.5*dg311
;

gammado122
=
-0.5*dg122 + dg212
;

gammado123
=
0.5*(-dg123 + dg213 + dg312)
;

gammado133
=
-0.5*dg133 + dg313
;

gammado211
=
dg112 - 0.5*dg211
;

gammado212
=
0.5*dg122
;

gammado213
=
0.5*(dg123 - dg213 + dg312)
;

gammado222
=
0.5*dg222
;

gammado223
=
0.5*dg322
;

gammado233
=
-0.5*dg233 + dg323
;

gammado311
=
dg113 - 0.5*dg311
;

gammado312
=
0.5*(dg123 + dg213 - dg312)
;

gammado313
=
0.5*dg133
;

gammado322
=
dg223 - 0.5*dg322
;

gammado323
=
0.5*dg233
;

gammado333
=
0.5*dg333
;

gamma111
=
gammado111*ginv11 + gammado211*ginv12 + gammado311*ginv13
;

gamma112
=
gammado112*ginv11 + gammado212*ginv12 + gammado312*ginv13
;

gamma113
=
gammado113*ginv11 + gammado213*ginv12 + gammado313*ginv13
;

gamma122
=
gammado122*ginv11 + gammado222*ginv12 + gammado322*ginv13
;

gamma123
=
gammado123*ginv11 + gammado223*ginv12 + gammado323*ginv13
;

gamma133
=
gammado133*ginv11 + gammado233*ginv12 + gammado333*ginv13
;

gamma211
=
gammado111*ginv12 + gammado211*ginv22 + gammado311*ginv23
;

gamma212
=
gammado112*ginv12 + gammado212*ginv22 + gammado312*ginv23
;

gamma213
=
gammado113*ginv12 + gammado213*ginv22 + gammado313*ginv23
;

gamma222
=
gammado122*ginv12 + gammado222*ginv22 + gammado322*ginv23
;

gamma223
=
gammado123*ginv12 + gammado223*ginv22 + gammado323*ginv23
;

gamma233
=
gammado133*ginv12 + gammado233*ginv22 + gammado333*ginv23
;

gamma311
=
gammado111*ginv13 + gammado211*ginv23 + gammado311*ginv33
;

gamma312
=
gammado112*ginv13 + gammado212*ginv23 + gammado312*ginv33
;

gamma313
=
gammado113*ginv13 + gammado213*ginv23 + gammado313*ginv33
;

gamma322
=
gammado122*ginv13 + gammado222*ginv23 + gammado322*ginv33
;

gamma323
=
gammado123*ginv13 + gammado223*ginv23 + gammado323*ginv33
;

gamma333
=
gammado133*ginv13 + gammado233*ginv23 + gammado333*ginv33
;

Gfromg1
=
gamma111*ginv11 + gamma122*ginv22 + 
  2.*(gamma112*ginv12 + gamma113*ginv13 + gamma123*ginv23) + gamma133*ginv33
;

Gfromg2
=
gamma211*ginv11 + gamma222*ginv22 + 
  2.*(gamma212*ginv12 + gamma213*ginv13 + gamma223*ginv23) + gamma233*ginv33
;

Gfromg3
=
gamma311*ginv11 + gamma322*ginv22 + 
  2.*(gamma312*ginv12 + gamma313*ginv13 + gamma323*ginv23) + gamma333*ginv33
;

dGfromgdu11
=
(ddg1111 - dg111*((8.*dg112 + 2.*dg211)*ginv12 + 
        (8.*dg113 + 2.*dg311)*ginv13) - 
     (dg113*(4.*dg112 + dg211) + dg112*dg311 + dg111*(dg213 + dg312))*
      ginv23 - ginv22*(dg112*dg211 + dg111*dg212 + 2.*pow2(dg112)) - 
     ginv33*(dg113*dg311 + dg111*dg313 + 2.*pow2(dg113)))*pow2(ginv11) + 
  (ddg1122 + ddg1212 - (dg123*(8.*dg112 + 2.*dg211) + 
        dg113*(4.*dg122 + 2.*dg212) + dg122*dg311 + 
        2.*(dg111*dg223 + dg112*(dg213 + dg312)) + dg111*dg322)*ginv13 - 
     (dg123*(4.*dg122 + 2.*dg212) + 
        2.*(dg113*dg222 + dg122*(dg213 + dg312) + dg112*(dg223 + dg322)))*
      ginv23 - ginv22*(3.*(dg122*dg212 + dg112*dg222) + 2.*pow2(dg122)) - 
     ginv33*(dg123*(dg213 + dg312) + dg122*dg313 + dg113*(dg223 + dg322) + 
        dg112*dg323 + 2.*pow2(dg123)))*pow2(ginv12) + 
  (ddg1133 + ddg1313 - (dg133*(4.*dg123 + 2.*(dg213 + dg312)) + 
        2.*(dg123*dg313 + dg113*(dg233 + dg323) + dg112*dg333))*ginv23 - 
     ginv22*(dg133*dg212 + dg113*dg223 + dg123*(dg213 + dg312) + 
        dg112*(dg233 + dg323) + 2.*pow2(dg123)) - 
     ginv33*(3.*(dg133*dg313 + dg113*dg333) + 2.*pow2(dg133)))*pow2(ginv13) \
+ ginv13*(ddg1333*ginv33 + ginv22*
      (ddg1223 - (dg133*dg222 + dg123*(4.*dg223 + dg322) + 
           dg122*(dg233 + dg323))*ginv23 - 
        (dg133*dg223 + dg123*(dg233 + 2.*dg323))*ginv33) + 
     ginv23*(ddg1233 + ddg1323 - 
        (dg133*(2.*dg233 + 3.*dg323) + 3.*dg123*dg333)*ginv33) - 
     (dg123*dg222 + dg122*dg223)*pow2(ginv22) - 
     (dg133*dg322 + 2.*(dg133*dg223 + dg123*(dg233 + dg323)) + 
        dg122*dg333)*pow2(ginv23) - 2.*dg133*dg333*pow2(ginv33)) + 
  ginv11*(ddg1313*ginv33 + ginv12*
      (2.*ddg1112 + ddg1211 - (dg113*(12.*dg112 + 3.*dg211) + 
           3.*dg112*dg311 + dg111*(8.*dg123 + 3.*(dg213 + dg312)))*ginv13 \
- (dg122*(4.*dg112 + dg211) + 6.*dg112*dg212 + dg111*dg222)*ginv22 - 
        (dg123*dg211 + dg122*dg311 + 
           4.*(dg113*(dg122 + dg212) + dg112*(dg123 + dg213 + dg312)) + 
           dg111*(dg223 + dg322))*ginv23 - 
        (dg123*dg311 + dg113*(4.*dg123 + 2.*(dg213 + dg312)) + 
           2.*dg112*dg313 + dg111*dg323)*ginv33) + 
     ginv22*(ddg1212 - (dg113*dg222 + 2.*(dg123*dg212 + dg112*dg223) + 
           dg122*(dg213 + dg312) + dg112*dg322)*ginv23 - 
        (dg113*dg223 + dg123*(dg213 + dg312) + dg112*dg323)*ginv33) + 
     ginv13*(2.*ddg1113 + ddg1311 - 
        (dg123*(4.*dg112 + dg211) + dg111*dg223 + 
           2.*(dg113*dg212 + dg112*(dg213 + dg312)))*ginv22 - 
        (dg133*dg211 + dg123*dg311 + 
           4.*(dg113*(dg123 + dg213 + dg312) + dg112*(dg133 + dg313)) + 
           dg111*(dg233 + dg323))*ginv23 - 
        (dg133*(4.*dg113 + dg311) + 6.*dg113*dg313 + dg111*dg333)*ginv33) + 
     ginv23*(ddg1213 + ddg1312 - 
        (dg133*(dg213 + dg312) + 2.*dg123*dg313 + 
           dg113*(dg233 + 2.*dg323) + dg112*dg333)*ginv33) - 
     (3.*dg112*dg211 + dg111*(4.*dg122 + 3.*dg212) + 6.*pow2(dg112))*
      pow2(ginv12) - (3.*dg113*dg311 + dg111*(4.*dg133 + 3.*dg313) + 
        6.*pow2(dg113))*pow2(ginv13) - 
     (dg122*dg212 + dg112*dg222)*pow2(ginv22) - 
     (dg133*dg212 + dg123*(dg213 + dg312) + dg122*dg313 + 
        dg113*(dg223 + dg322) + dg112*(dg233 + dg323))*pow2(ginv23) - 
     (dg133*dg313 + dg113*dg333)*pow2(ginv33)) + 
  ginv12*(ddg1323*ginv33 + ginv22*
      (ddg1222 - (3.*(dg123*dg222 + dg122*dg223) + 2.*dg122*dg322)*
         ginv23 - (dg123*(2.*dg223 + dg322) + dg122*dg323)*ginv33) + 
     ginv23*(ddg1223 + ddg1322 - 
        (dg133*(dg223 + dg322) + dg123*(dg233 + 4.*dg323) + dg122*dg333)*
         ginv33) + ginv13*(2.*ddg1123 + ddg1213 + ddg1312 - 
        (dg113*dg222 + 4.*(dg123*(dg122 + dg212) + dg112*dg223) + 
           dg122*(dg213 + dg312) + dg112*dg322)*ginv22 - 
        (dg133*(4.*dg123 + dg213 + dg312) + 4.*dg123*dg313 + 
           dg113*(dg233 + 4.*dg323) + dg112*dg333)*ginv33 - 
        ginv23*(2.*(dg133*dg212 + dg112*dg233 + dg122*dg313 + 
              dg113*dg322) + 4.*
            (dg122*dg133 + dg113*dg223 + dg123*(dg213 + dg312) + 
              dg112*dg323 + pow2(dg123)))) - 
     (dg133*(4.*dg112 + dg211) + dg113*(8.*dg123 + 2.*(dg213 + dg312)) + 
        2.*(dg123*dg311 + dg112*dg313) + dg111*(dg233 + 2.*dg323))*
      pow2(ginv13) - 2.*dg122*dg222*pow2(ginv22) - 
     (dg133*dg222 + 2.*dg123*(dg223 + dg322) + dg122*(dg233 + 2.*dg323))*
      pow2(ginv23) - (dg133*dg323 + dg123*dg333)*pow2(ginv33)) - 
  2.*pow2(dg111)*pow3(ginv11) - 
  (dg122*(4.*dg112 + dg211) + 2.*dg112*dg212 + dg111*dg222)*pow3(ginv12) - 
  (dg133*(4.*dg113 + dg311) + 2.*dg113*dg313 + dg111*dg333)*pow3(ginv13)
;

dGfromgdu12
=
(ddg1112 + ddg1211 - (4.*(dg112*dg113 + dg111*dg123) + 
        2.*(dg113*dg211 + dg112*dg311 + dg111*(dg213 + dg312)))*ginv13 - 
     (dg122*(6.*dg112 + 2.*dg211) + 6.*dg112*dg212 + 2.*dg111*dg222)*
      ginv22 - (4.*(dg113*(dg122 + dg212) + dg112*(dg123 + dg213)) + 
        dg122*dg311 + 2.*(dg123*dg211 + dg111*dg223 + dg112*dg312) + 
        dg111*dg322)*ginv23 - (dg123*dg311 + 
        dg113*(2.*(dg123 + dg213) + dg312) + dg112*dg313 + dg111*dg323)*
      ginv33)*pow2(ginv12) - ((2.*(dg113*dg123 + dg112*dg133) + 
        dg123*dg311 + dg113*dg312 + dg112*dg313 + dg111*dg323)*ginv22 + 
     (dg133*(4.*dg113 + dg311) + 2.*dg113*dg313 + dg111*dg333)*ginv23)*
   pow2(ginv13) + (ddg1222 - (4.*(dg123*dg222 + dg122*dg223) + 
        2.*dg122*dg322)*ginv23 - 
     (dg123*(2.*dg223 + dg322) + dg122*dg323)*ginv33)*pow2(ginv22) + 
  (ddg1233 + ddg1323 - (dg133*(2.*dg233 + 3.*dg323) + 3.*dg123*dg333)*
      ginv33)*pow2(ginv23) + ginv11*
   (ginv23*(ddg1113 - 2.*dg113*(dg133 + dg313)*ginv33) + 
     ginv22*(ddg1112 - (dg112*(4.*dg123 + 2.*dg213) + 
           2.*(dg113*(dg122 + dg212) + dg112*dg312))*ginv23 - 
        (dg113*(2.*dg123 + dg312) + dg112*dg313)*ginv33) + 
     ginv12*(ddg1111 - dg111*(6.*dg113 + 2.*dg311)*ginv13 - 
        (dg113*(8.*dg112 + 2.*dg211) + dg112*dg311 + 
           dg111*(2.*(dg123 + dg213) + dg312))*ginv23 - 
        ginv22*(2.*(dg112*dg211 + dg111*(dg122 + dg212)) + 
           6.*pow2(dg112)) - ginv33*
         (dg113*dg311 + dg111*dg313 + 2.*pow2(dg113))) - 
     ginv13*((dg112*(4.*dg113 + dg311) + dg111*(2.*dg123 + dg312))*
         ginv22 + ginv23*(dg113*dg311 + dg111*(2.*dg133 + dg313) + 
           4.*pow2(dg113))) - dg111*(6.*dg112 + 2.*dg211)*pow2(ginv12) - 
     2.*dg112*(dg122 + dg212)*pow2(ginv22) - 
     (2.*(dg112*dg133 + dg113*(dg123 + dg213)) + dg113*dg312 + dg112*dg313)*
      pow2(ginv23)) + ginv13*(ginv22*
      (ddg1123 + ddg1312 - (dg133*(2.*dg123 + dg312) + 
           2.*(dg123*dg313 + dg113*dg323) + dg112*dg333)*ginv33 - 
        ginv23*(2.*(dg133*(dg122 + dg212) + dg123*dg213 + dg113*dg223 + 
              dg112*dg233) + dg122*dg313 + dg113*dg322 + 
           4.*(dg123*dg312 + dg112*dg323 + pow2(dg123)))) + 
     ginv23*(ddg1133 + ddg1313 - 
        ginv33*(3.*(dg133*dg313 + dg113*dg333) + 2.*pow2(dg133))) - 
     (2.*(dg123*(dg122 + dg212) + dg112*dg223) + dg122*dg312 + 
        dg112*dg322)*pow2(ginv22) - 
     (dg133*(4.*dg123 + 2.*(dg213 + dg312)) + 
        2.*(dg123*dg313 + dg113*(dg233 + dg323) + dg112*dg333))*pow2(ginv23)\
) + ginv23*(ddg1333*ginv33 - 2.*dg133*dg333*pow2(ginv33)) + 
  ginv12*(ddg1313*ginv33 + ginv13*
      (ddg1113 + ddg1311 - (2.*
            (dg123*dg211 + dg113*(dg122 + dg212) + dg111*dg223) + 
           dg122*dg311 + dg112*(8.*dg123 + 2.*dg213 + 4.*dg312) + 
           dg111*dg322)*ginv22 - 
        (dg133*(4.*dg112 + 2.*dg211) + 
           dg113*(8.*dg123 + 4.*(dg213 + dg312)) + 4.*dg112*dg313 + 
           2.*(dg123*dg311 + dg111*(dg233 + dg323)))*ginv23 - 
        (dg133*(2.*dg113 + dg311) + 4.*dg113*dg313 + dg111*dg333)*ginv33) + 
     ginv23*(ddg1123 + 2.*ddg1213 + ddg1312 - 
        (2.*(dg133*(dg123 + dg213) + dg113*dg233) + dg133*dg312 + 
           4.*(dg123*dg313 + dg113*dg323) + dg112*dg333)*ginv33) + 
     ginv22*(ddg1122 + 2.*ddg1212 - 
        (4.*(dg122*dg213 + dg113*dg222) + 
           6.*(dg123*(dg122 + dg212) + dg112*dg223) + 
           3.*(dg122*dg312 + dg112*dg322))*ginv23 - 
        ginv33*(dg122*dg313 + dg113*dg322 + 
           2.*(dg113*dg223 + dg123*(dg213 + dg312) + dg112*dg323 + 
              pow2(dg123)))) - 
     2.*(dg113*dg311 + dg111*(dg133 + dg313) + pow2(dg113))*pow2(ginv13) - 
     (4.*(dg122*dg212 + dg112*dg222) + 2.*pow2(dg122))*pow2(ginv22) - 
     (4.*(dg123*dg213 + dg113*dg223) + 
        2.*(dg133*(dg122 + dg212) + dg123*dg312 + dg122*dg313 + 
           dg113*dg322 + dg112*(dg233 + dg323) + pow2(dg123)))*pow2(ginv23) \
- (dg133*dg313 + dg113*dg333)*pow2(ginv33)) + 
  ginv22*(ddg1323*ginv33 + ginv23*
      (2.*ddg1223 + ddg1322 - (2.*(dg133*dg223 + dg123*dg233) + 
           dg133*dg322 + 6.*dg123*dg323 + dg122*dg333)*ginv33) - 
     (2.*(dg133*dg222 + dg122*dg233) + dg123*(6.*dg223 + 3.*dg322) + 
        3.*dg122*dg323)*pow2(ginv23) - 
     (dg133*dg323 + dg123*dg333)*pow2(ginv33)) - 
  2.*((dg111*(dg112*ginv22 + dg113*ginv23) + ginv12*pow2(dg111))*
      pow2(ginv11) + (dg112*dg211 + dg111*(dg122 + dg212) + pow2(dg112))*
      pow3(ginv12) + dg122*dg222*pow3(ginv22)) - 
  (dg133*dg322 + 2.*(dg133*dg223 + dg123*(dg233 + dg323)) + dg122*dg333)*
   pow3(ginv23)
;

dGfromgdu13
=
-(((dg122*(4.*dg112 + dg211) + 2.*dg112*dg212 + dg111*dg222)*ginv23 + 
       (2.*(dg113*dg122 + dg112*dg123) + dg123*dg211 + dg113*dg212 + 
          dg112*dg213 + dg111*dg223)*ginv33 + 
       2.*ginv13*(dg112*dg211 + dg111*(dg122 + dg212) + pow2(dg112)))*
     pow2(ginv12)) + (ddg1113 + ddg1311 - 
     (dg123*(2.*dg112 + dg211) + dg113*dg212 + dg111*dg223 + 
        dg112*(dg213 + 2.*dg312))*ginv22 - 
     (dg133*dg211 + 2.*(dg113*dg213 + dg123*dg311) + 
        4.*(dg113*(dg123 + dg312) + dg112*(dg133 + dg313)) + 
        dg111*(dg233 + 2.*dg323))*ginv23 - 
     (dg133*(6.*dg113 + 2.*dg311) + 6.*dg113*dg313 + 2.*dg111*dg333)*ginv33\
)*pow2(ginv13) - (2.*dg122*dg222*ginv23 + 
     (dg123*dg222 + dg122*dg223)*ginv33)*pow2(ginv22) + 
  (ddg1223 + ddg1322 - (3.*(dg133*dg223 + dg123*dg233) + 6.*dg123*dg323 + 
        2.*(dg133*dg322 + dg122*dg333))*ginv33)*pow2(ginv23) + 
  ddg1333*pow2(ginv33) + ginv11*
   (ddg1113*ginv33 - ginv22*(2.*dg112*(dg122 + dg212)*ginv23 + 
        (dg113*dg212 + dg112*(2.*dg123 + dg213))*ginv33) + 
     ginv23*(ddg1112 - (dg113*(4.*dg123 + 2.*dg213) + 
           2.*(dg113*dg312 + dg112*(dg133 + dg313)))*ginv33) - 
     ginv12*(dg111*(6.*dg112 + 2.*dg211)*ginv13 + 
        (dg113*(4.*dg112 + dg211) + dg111*(2.*dg123 + dg213))*ginv33 + 
        ginv23*(dg112*dg211 + dg111*(2.*dg122 + dg212) + 4.*pow2(dg112))) + 
     ginv13*(ddg1111 - (dg113*(8.*dg112 + dg211) + 2.*dg112*dg311 + 
           dg111*(dg213 + 2.*(dg123 + dg312)))*ginv23 - 
        ginv22*(dg112*dg211 + dg111*dg212 + 2.*pow2(dg112)) - 
        ginv33*(2.*(dg113*dg311 + dg111*(dg133 + dg313)) + 6.*pow2(dg113))) \
- dg111*(6.*dg113 + 2.*dg311)*pow2(ginv13) - 
     (dg113*dg212 + dg112*dg213 + 
        2.*(dg113*dg122 + dg112*(dg123 + dg312)))*pow2(ginv23) - 
     2.*dg113*(dg133 + dg313)*pow2(ginv33)) + 
  ginv12*((ddg1123 + ddg1213)*ginv33 + 
     ginv13*(ddg1112 + ddg1211 - 
        (dg122*(2.*dg112 + dg211) + 4.*dg112*dg212 + dg111*dg222)*ginv22 - 
        (dg123*(8.*dg112 + 2.*dg211) + 
           4.*(dg113*(dg122 + dg212) + dg112*(dg213 + dg312)) + 
           2.*(dg122*dg311 + dg111*(dg223 + dg322)))*ginv23 - 
        (dg133*(2.*dg112 + dg211) + 
           dg113*(8.*dg123 + 4.*dg213 + 2.*dg312) + 
           2.*(dg123*dg311 + dg112*dg313) + dg111*(dg233 + 2.*dg323))*
         ginv33) - ginv22*((dg122*dg213 + dg113*dg222 + 
           2.*(dg123*(dg122 + dg212) + dg112*dg223))*ginv33 + 
        ginv23*(3.*(dg122*dg212 + dg112*dg222) + 2.*pow2(dg122))) + 
     ginv23*(ddg1122 + ddg1212 - 
        ginv33*(dg133*(2.*dg122 + dg212) + 
           2.*(dg123*dg312 + dg122*dg313 + dg113*dg322) + 
           dg112*(dg233 + 2.*dg323) + 
           4.*(dg123*dg213 + dg113*dg223 + pow2(dg123)))) - 
     (4.*(dg112*dg113 + dg111*dg123) + 
        2.*(dg113*dg211 + dg112*dg311 + dg111*(dg213 + dg312)))*
      pow2(ginv13) - (dg123*(4.*dg122 + 2.*dg212) + 
        2.*(dg113*dg222 + dg122*(dg213 + dg312) + dg112*(dg223 + dg322)))*
      pow2(ginv23) - (dg133*(2.*dg123 + dg213) + 2.*dg123*dg313 + 
        dg113*(dg233 + 2.*dg323))*pow2(ginv33)) + 
  ginv22*(ddg1223*ginv33 + ginv23*
      (ddg1222 - (dg133*dg222 + dg123*(6.*dg223 + 2.*dg322) + 
           dg122*(dg233 + 2.*dg323))*ginv33) - 
     (3.*(dg123*dg222 + dg122*dg223) + 2.*dg122*dg322)*pow2(ginv23) - 
     (dg133*dg223 + dg123*(dg233 + 2.*dg323))*pow2(ginv33)) + 
  ginv23*((ddg1233 + 2.*ddg1323)*ginv33 - 
     (dg133*(2.*dg233 + 4.*dg323) + 4.*dg123*dg333)*pow2(ginv33)) + 
  ginv13*((ddg1133 + 2.*ddg1313)*ginv33 + 
     ginv23*(ddg1123 + ddg1213 + 2.*ddg1312 - 
        (dg133*(6.*dg123 + 3.*dg213 + 4.*dg312) + 6.*dg123*dg313 + 
           dg113*(3.*dg233 + 6.*dg323) + 4.*dg112*dg333)*ginv33) + 
     ginv22*(ddg1212 - (dg123*(2.*dg122 + 4.*dg212) + dg113*dg222 + 
           dg122*(dg213 + 2.*dg312) + dg112*(4.*dg223 + 2.*dg322))*ginv23 \
- ginv33*(dg133*dg212 + dg112*(dg233 + 2.*dg323) + 
           2.*(dg113*dg223 + dg123*(dg213 + dg312) + pow2(dg123)))) - 
     (dg122*dg212 + dg112*dg222)*pow2(ginv22) - 
     (4.*(dg123*dg312 + dg112*dg323) + 
        2.*(dg133*(dg122 + dg212) + dg123*dg213 + dg112*dg233 + 
           dg122*dg313 + dg113*(dg223 + dg322) + pow2(dg123)))*pow2(ginv23) \
- (4.*(dg133*dg313 + dg113*dg333) + 2.*pow2(dg133))*pow2(ginv33)) - 
  (dg133*dg222 + 2.*dg123*(dg223 + dg322) + dg122*(dg233 + 2.*dg323))*
   pow3(ginv23) - 2.*((dg111*(dg112*ginv23 + dg113*ginv33) + 
        ginv13*pow2(dg111))*pow2(ginv11) + 
     (dg113*dg311 + dg111*(dg133 + dg313) + pow2(dg113))*pow3(ginv13) + 
     dg133*dg333*pow3(ginv33))
;

dGfromgdu21
=
(ddg1211 - (4.*(dg113*dg211 + dg111*dg213) + 2.*dg211*dg311)*ginv13 - 
     2.*(dg112 + dg211)*dg212*ginv22 - 
     (2.*(dg113*dg212 + (dg112 + dg211)*dg213) + dg212*dg311 + 
        dg211*dg312)*ginv23 - (dg213*(2.*dg113 + dg311) + dg211*dg313)*
      ginv33 - ginv12*(4.*(dg112*dg211 + dg111*dg212) + 2.*pow2(dg211)))*
   pow2(ginv11) + (ddg1222 + ddg2212 - 
     (4.*(dg212*(dg123 + dg213) + (dg112 + dg211)*dg223) + dg222*dg311 + 
        2.*(dg122*dg213 + dg113*dg222 + dg212*dg312) + dg211*dg322)*ginv13 \
- (2.*dg122 + 6.*dg212)*dg222*ginv22 - 
     ((2.*dg122 + 4.*dg212)*dg223 + 
        dg222*(4.*dg213 + 2.*(dg123 + dg312)) + 2.*dg212*dg322)*ginv23 - 
     (dg223*(2.*(dg123 + dg213) + dg312) + dg222*dg313 + dg213*dg322 + 
        dg212*dg323)*ginv33)*pow2(ginv12) + 
  (ddg1233 + ddg2313 - (2.*((dg123 + dg213)*dg223 + dg212*dg233) + 
        dg223*dg312 + dg212*dg323)*ginv22 - 
     (dg233*(4.*dg213 + 2.*dg312) + 
        2.*(dg123*dg233 + dg223*(dg133 + dg313) + dg213*dg323 + 
           dg212*dg333))*ginv23 - 
     (dg233*(2.*dg133 + 3.*dg313) + 3.*dg213*dg333)*ginv33)*pow2(ginv13) + 
  ginv11*(ddg2313*ginv33 + ginv22*
      (ddg2212 - (dg222*(2.*dg213 + dg312) + dg212*(4.*dg223 + dg322))*
         ginv23 - (dg223*(2.*dg213 + dg312) + dg212*dg323)*ginv33) + 
     ginv23*(ddg2213 + ddg2312 - 
        (dg233*(2.*dg213 + dg312) + 2.*(dg223*dg313 + dg213*dg323) + 
           dg212*dg333)*ginv33) + 
     ginv13*(2.*ddg1213 + ddg2311 - 
        (2.*(dg112 + dg211)*dg223 + 
           dg212*(4.*dg213 + 2.*(dg123 + dg312)))*ginv22 - 
        (2.*(dg133*dg213 + dg113*dg233) + dg233*dg311 + 6.*dg213*dg313 + 
           dg211*dg333)*ginv33 - 
        ginv23*(2.*(dg133*dg212 + dg123*dg213 + dg113*dg223 + 
              (dg112 + dg211)*dg233) + dg223*dg311 + dg211*dg323 + 
           4.*(dg213*dg312 + dg212*dg313 + pow2(dg213)))) + 
     ginv12*(2.*ddg1212 + ddg2211 - 
        (6.*(dg113*dg212 + dg112*dg213) + 4.*dg111*dg223 + 
           3.*dg212*dg311 + dg211*(4.*dg123 + 6.*dg213 + 3.*dg312))*ginv13 \
- (2.*(dg123*dg212 + dg122*dg213 + dg113*dg222 + 
              (dg112 + dg211)*dg223) + dg222*dg311 + 
           dg212*(8.*dg213 + 4.*dg312) + dg211*dg322)*ginv23 - 
        ginv22*(2.*(dg122*dg212 + (dg112 + dg211)*dg222) + 
           6.*pow2(dg212)) - ginv33*
         (dg223*dg311 + dg211*dg323 + 
           2.*(dg113*dg223 + dg213*(dg123 + dg312) + dg212*dg313 + 
              pow2(dg213)))) - 
     (6.*dg112*dg212 + dg211*(2.*dg122 + 6.*dg212) + 2.*dg111*dg222)*
      pow2(ginv12) - (2.*(dg133*dg211 + dg111*dg233) + 
        dg213*(6.*dg113 + 3.*dg311) + 3.*dg211*dg313)*pow2(ginv13) - 
     2.*dg212*dg222*pow2(ginv22) - 
     (2.*(dg213*dg223 + dg212*dg233) + dg223*dg312 + dg222*dg313 + 
        dg213*dg322 + dg212*dg323)*pow2(ginv23) - 
     (dg233*dg313 + dg213*dg333)*pow2(ginv33)) + 
  ginv12*(ddg2323*ginv33 + ginv13*
      (2.*ddg1223 + ddg2213 + ddg2312 - 
        (2.*((dg123 + dg213)*dg222 + dg122*dg223) + dg222*dg312 + 
           dg212*(8.*dg223 + dg322))*ginv22 - 
        (dg223*(8.*dg213 + 4.*(dg123 + dg312)) + 
           2.*(dg122*dg233 + dg222*(dg133 + dg313) + dg213*dg322) + 
           4.*dg212*(dg233 + dg323))*ginv23 - 
        (2.*(dg133*dg223 + (dg123 + dg213)*dg233) + dg233*dg312 + 
           4.*(dg223*dg313 + dg213*dg323) + dg212*dg333)*ginv33) + 
     ginv23*(ddg2223 + ddg2322 - 
        (dg233*(2.*dg223 + dg322) + 4.*dg223*dg323 + dg222*dg333)*ginv33) + 
     ginv22*(ddg2222 - dg222*(6.*dg223 + 2.*dg322)*ginv23 - 
        ginv33*(dg223*dg322 + dg222*dg323 + 2.*pow2(dg223))) - 
     (4.*(dg123*dg213 + dg113*dg223) + 
        2.*((dg112 + dg211)*dg233 + dg223*dg311 + dg213*dg312 + 
           dg212*(dg133 + dg313) + dg211*dg323 + pow2(dg213)))*pow2(ginv13) \
- 2.*(pow2(dg222)*pow2(ginv22) + 
        (dg223*dg322 + dg222*(dg233 + dg323) + pow2(dg223))*pow2(ginv23)) - 
     (dg233*dg323 + dg223*dg333)*pow2(ginv33)) + 
  ginv13*(ddg2333*ginv33 + ginv22*
      (ddg2223 - 2.*dg223*(dg233 + dg323)*ginv33 - 
        ginv23*(dg223*dg322 + dg222*(2.*dg233 + dg323) + 4.*pow2(dg223))) + 
     ginv23*(ddg2233 + ddg2323 - 
        ginv33*(3.*(dg233*dg323 + dg223*dg333) + 2.*pow2(dg233))) - 
     (dg233*(4.*dg223 + dg322) + 2.*dg223*dg323 + dg222*dg333)*
      pow2(ginv23) - 2.*(dg222*dg223*pow2(ginv22) + 
        dg233*dg333*pow2(ginv33))) - 
  2.*(dg111*dg211*pow3(ginv11) + 
     (dg122*dg212 + (dg112 + dg211)*dg222 + pow2(dg212))*pow3(ginv12)) - 
  (dg233*dg311 + 2.*(dg113*dg233 + dg213*(dg133 + dg313)) + dg211*dg333)*
   pow3(ginv13)
;

dGfromgdu22
=
-((2.*dg111*dg211*ginv12 + (dg112*dg211 + dg111*dg212)*ginv22 + 
       (dg113*dg211 + dg111*dg213)*ginv23)*pow2(ginv11)) + 
  (ddg1212 + ddg2211 - (2.*(dg123*dg211 + dg112*dg213 + dg111*dg223 + 
           dg212*(dg113 + dg311)) + dg211*(4.*dg213 + 2.*dg312))*ginv13 - 
     (2.*(dg123*dg212 + dg122*dg213 + dg113*dg222 + dg112*dg223) + 
        dg222*dg311 + dg212*(8.*dg213 + 2.*dg312) + 
        dg211*(4.*dg223 + dg322))*ginv23 - 
     ginv22*(4.*dg211*dg222 + 3.*(dg122*dg212 + dg112*dg222) + 
        6.*pow2(dg212)) - ginv33*
      (dg223*(dg113 + dg311) + dg213*(dg123 + dg312) + dg212*dg313 + 
        dg211*dg323 + 2.*pow2(dg213)))*pow2(ginv12) - 
  ((dg112*dg233 + dg223*(dg113 + dg311) + dg213*(dg123 + dg312) + 
        dg212*(dg133 + dg313) + dg211*dg323)*ginv22 + 
     (dg233*dg311 + 2.*(dg113*dg233 + dg213*(dg133 + dg313)) + 
        dg211*dg333)*ginv23)*pow2(ginv13) + 
  (ddg2222 - dg222*(8.*dg223 + 2.*dg322)*ginv23 - 
     ginv33*(dg223*dg322 + dg222*dg323 + 2.*pow2(dg223)))*pow2(ginv22) + 
  (ddg2233 + ddg2323 - ginv33*(3.*(dg233*dg323 + dg223*dg333) + 
        2.*pow2(dg233)))*pow2(ginv23) + 
  ginv13*(ginv22*(ddg1223 + ddg2312 - 
        (dg122*dg233 + dg222*(dg133 + dg313) + dg213*dg322 + 
           4.*(dg223*(dg123 + dg213 + dg312) + dg212*(dg233 + dg323)))*
         ginv23 - (dg233*(dg123 + dg312) + dg223*(dg133 + 2.*dg313) + 
           2.*dg213*dg323 + dg212*dg333)*ginv33) + 
     ginv23*(ddg1233 + ddg2313 - 
        (dg233*(2.*dg133 + 3.*dg313) + 3.*dg213*dg333)*ginv33) - 
     ((dg122 + 4.*dg212)*dg223 + dg222*(dg123 + dg312) + dg212*dg322)*
      pow2(ginv22) - (dg233*(4.*dg213 + 2.*dg312) + 
        2.*(dg123*dg233 + dg223*(dg133 + dg313) + dg213*dg323 + 
           dg212*dg333))*pow2(ginv23)) + 
  ginv11*(-(ginv13*((2.*(dg113*dg212 + dg112*dg213) + dg111*dg223 + 
             dg212*dg311 + dg211*(dg123 + dg312))*ginv22 + 
          (dg111*dg233 + dg213*(4.*dg113 + dg311) + dg211*(dg133 + dg313))*
           ginv23)) + ginv12*(ddg1211 - 
        (3.*(dg113*dg211 + dg111*dg213) + 2.*dg211*dg311)*ginv13 - 
        (6.*dg112*dg212 + dg211*(dg122 + 4.*dg212) + dg111*dg222)*ginv22 - 
        (4.*(dg113*dg212 + dg112*dg213) + dg111*dg223 + dg212*dg311 + 
           dg211*(dg123 + 4.*dg213 + dg312))*ginv23 - 
        (dg213*(2.*dg113 + dg311) + dg211*dg313)*ginv33) + 
     ginv22*(ddg1212 - (dg122*dg213 + dg113*dg222 + 2.*dg112*dg223 + 
           dg212*(4.*dg213 + 2.*(dg123 + dg312)))*ginv23 - 
        (dg113*dg223 + dg213*(dg123 + dg312) + dg212*dg313)*ginv33) + 
     ginv23*(ddg1213 - (dg113*dg233 + dg213*(dg133 + 2.*dg313))*ginv33) - 
     (3.*(dg112*dg211 + dg111*dg212) + 2.*pow2(dg211))*pow2(ginv12) - 
     (dg122*dg212 + dg112*dg222 + 2.*pow2(dg212))*pow2(ginv22) - 
     (dg113*dg223 + dg112*dg233 + dg213*(dg123 + dg312) + 
        dg212*(dg133 + dg313) + 2.*pow2(dg213))*pow2(ginv23)) + 
  ginv23*(ddg2333*ginv33 - 2.*dg233*dg333*pow2(ginv33)) + 
  ginv12*(ddg2313*ginv33 + ginv22*
      (ddg1222 + 2.*ddg2212 - ((3.*dg122 + 12.*dg212)*dg223 + 
           dg222*(8.*dg213 + 3.*(dg123 + dg312)) + 3.*dg212*dg322)*ginv23 \
- (dg223*(4.*dg213 + 2.*(dg123 + dg312)) + dg222*dg313 + dg213*dg322 + 
           2.*dg212*dg323)*ginv33) + 
     ginv23*(ddg1223 + 2.*ddg2213 + ddg2312 - 
        (dg233*(dg123 + 4.*dg213 + dg312) + dg223*(dg133 + 4.*dg313) + 
           4.*dg213*dg323 + dg212*dg333)*ginv33) + 
     ginv13*(ddg1213 + ddg2311 - 
        (dg122*dg213 + dg222*(dg113 + dg311) + 
           4.*((dg112 + dg211)*dg223 + dg212*(dg123 + dg213 + dg312)) + 
           dg211*dg322)*ginv22 - 
        (dg233*(dg113 + dg311) + dg213*(dg133 + 4.*dg313) + dg211*dg333)*
         ginv33 - ginv23*(2.*(dg133*dg212 + dg112*dg233 + dg223*dg311 + 
              dg211*dg323) + 4.*
            (dg113*dg223 + dg211*dg233 + dg213*(dg123 + dg312) + 
              dg212*dg313 + pow2(dg213)))) - 
     (dg111*dg233 + 2.*dg213*(dg113 + dg311) + dg211*(dg133 + 2.*dg313))*
      pow2(ginv13) - (2.*dg122 + 8.*dg212)*dg222*pow2(ginv22) - 
     ((dg122 + 4.*dg212)*dg233 + dg223*(8.*dg213 + 2.*(dg123 + dg312)) + 
        dg222*(dg133 + 2.*dg313) + 2.*(dg213*dg322 + dg212*dg323))*
      pow2(ginv23) - (dg233*dg313 + dg213*dg333)*pow2(ginv33)) + 
  ginv22*(ddg2323*ginv33 + ginv23*
      (2.*ddg2223 + ddg2322 - (dg233*(4.*dg223 + dg322) + 
           6.*dg223*dg323 + dg222*dg333)*ginv33) - 
     (3.*dg223*dg322 + dg222*(4.*dg233 + 3.*dg323) + 6.*pow2(dg223))*
      pow2(ginv23) - (dg233*dg323 + dg223*dg333)*pow2(ginv33)) - 
  (2.*dg112*dg212 + dg211*(dg122 + 4.*dg212) + dg111*dg222)*pow3(ginv12) - 
  2.*pow2(dg222)*pow3(ginv22) - 
  (dg233*(4.*dg223 + dg322) + 2.*dg223*dg323 + dg222*dg333)*pow3(ginv23)
;

dGfromgdu23
=
-((2.*dg111*dg211*ginv13 + (dg112*dg211 + dg111*dg212)*ginv23 + 
       (dg113*dg211 + dg111*dg213)*ginv33)*pow2(ginv11)) - 
  ((2.*dg112*dg212 + dg211*(dg122 + 4.*dg212) + dg111*dg222)*ginv13 + 
     (dg122*dg213 + dg212*(dg123 + 2.*dg213) + dg113*dg222 + 
        (dg112 + 2.*dg211)*dg223)*ginv33 + 
     2.*ginv23*(dg122*dg212 + (dg112 + dg211)*dg222 + pow2(dg212)))*
   pow2(ginv12) + (ddg1213 + ddg2311 - 
     ((dg112 + 2.*dg211)*dg223 + dg212*(dg123 + 2.*(dg213 + dg312)))*
      ginv22 - (3.*(dg133*dg213 + dg113*dg233) + 6.*dg213*dg313 + 
        2.*(dg233*dg311 + dg211*dg333))*ginv33 - 
     ginv23*(4.*(dg213*dg312 + dg212*dg313) + 
        2.*(dg133*dg212 + dg123*dg213 + (dg112 + dg211)*dg233 + 
           dg223*(dg113 + dg311) + dg211*dg323 + pow2(dg213))))*pow2(ginv13) \
+ (ddg2223 + ddg2322 - (dg233*(6.*dg223 + 2.*dg322) + 6.*dg223*dg323 + 
        2.*dg222*dg333)*ginv33)*pow2(ginv23) + ddg2333*pow2(ginv33) + 
  ginv11*(ddg1213*ginv33 + ginv13*
      (ddg1211 - 2.*(dg112 + dg211)*dg212*ginv22 - 
        (4.*(dg113*dg212 + dg112*dg213) + dg111*dg223 + 2.*dg212*dg311 + 
           dg211*(dg123 + 2.*(dg213 + dg312)))*ginv23 - 
        (dg111*dg233 + dg213*(6.*dg113 + 2.*dg311) + 
           dg211*(dg133 + 2.*dg313))*ginv33) - 
     ginv12*((4.*dg112*dg212 + dg211*(dg122 + 2.*dg212) + dg111*dg222)*
         ginv23 + (dg211*(dg123 + 2.*dg213) + 
           2.*(dg113*dg212 + dg112*dg213) + dg111*dg223)*ginv33 + 
        ginv13*(3.*(dg112*dg211 + dg111*dg212) + 2.*pow2(dg211))) - 
     ginv22*((dg212*(dg123 + 2.*dg213) + dg112*dg223)*ginv33 + 
        ginv23*(dg122*dg212 + dg112*dg222 + 2.*pow2(dg212))) + 
     ginv23*(ddg1212 - ginv33*(dg112*dg233 + dg212*(dg133 + 2.*dg313) + 
           2.*(dg113*dg223 + dg213*(dg123 + dg312) + pow2(dg213)))) - 
     (3.*(dg113*dg211 + dg111*dg213) + 2.*dg211*dg311)*pow2(ginv13) - 
     (dg122*dg213 + dg113*dg222 + dg112*dg223 + 
        dg212*(dg123 + 2.*(dg213 + dg312)))*pow2(ginv23) - 
     (dg113*dg233 + dg213*(dg133 + 2.*dg313))*pow2(ginv33)) + 
  ginv22*(ddg2223*ginv33 + ginv23*
      (ddg2222 - ginv33*(2.*(dg223*dg322 + dg222*(dg233 + dg323)) + 
           6.*pow2(dg223))) - dg222*(6.*dg223 + 2.*dg322)*pow2(ginv23) - 
     2.*dg223*(dg233 + dg323)*pow2(ginv33)) + 
  ginv12*((ddg1223 + ddg2213)*ginv33 - 
     ginv22*((2.*dg122 + 6.*dg212)*dg222*ginv23 + 
        ((dg123 + 2.*dg213)*dg222 + (dg122 + 4.*dg212)*dg223)*ginv33) + 
     ginv23*(ddg1222 + ddg2212 - 
        ((dg122 + 2.*dg212)*dg233 + 
           dg223*(4.*dg123 + 8.*dg213 + 2.*dg312) + 
           dg222*(dg133 + 2.*dg313) + 2.*(dg213*dg322 + dg212*dg323))*
         ginv33) + ginv13*(ddg1212 + ddg2211 - 
        (4.*(dg112 + dg211)*dg223 + 
           dg212*(8.*dg213 + 4.*(dg123 + dg312)) + 
           2.*(dg122*dg213 + dg222*(dg113 + dg311) + dg211*dg322))*ginv23 \
- ginv22*(dg122*dg212 + (dg112 + 2.*dg211)*dg222 + 4.*pow2(dg212)) - 
        ginv33*((dg112 + 2.*dg211)*dg233 + dg212*(dg133 + 2.*dg313) + 
           2.*(dg223*dg311 + dg213*dg312 + dg211*dg323) + 
           4.*(dg123*dg213 + dg113*dg223 + pow2(dg213)))) - 
     (2.*(dg123*dg211 + dg112*dg213 + dg111*dg223 + 
           dg212*(dg113 + dg311)) + dg211*(4.*dg213 + 2.*dg312))*
      pow2(ginv13) - ((2.*dg122 + 4.*dg212)*dg223 + 
        dg222*(4.*dg213 + 2.*(dg123 + dg312)) + 2.*dg212*dg322)*
      pow2(ginv23) - ((dg123 + 2.*dg213)*dg233 + 
        dg223*(dg133 + 2.*dg313) + 2.*dg213*dg323)*pow2(ginv33)) + 
  ginv13*((ddg1233 + 2.*ddg2313)*ginv33 + 
     ginv22*(ddg2212 - ((dg122 + 8.*dg212)*dg223 + 
           dg222*(dg123 + 2.*(dg213 + dg312)) + 2.*dg212*dg322)*ginv23 - 
        (dg223*(4.*dg213 + 2.*(dg123 + dg312)) + 2.*dg212*(dg233 + dg323))*
         ginv33) + ginv23*(ddg1223 + ddg2213 + 2.*ddg2312 - 
        (3.*(dg133*dg223 + dg123*dg233) + dg233*(6.*dg213 + 4.*dg312) + 
           6.*(dg223*dg313 + dg213*dg323) + 4.*dg212*dg333)*ginv33) - 
     2.*dg212*dg222*pow2(ginv22) - 
     ((dg122 + 4.*dg212)*dg233 + dg223*(2.*dg123 + 4.*(dg213 + dg312)) + 
        dg222*(dg133 + 2.*dg313) + 2.*dg213*dg322 + 4.*dg212*dg323)*
      pow2(ginv23) - (dg233*(2.*dg133 + 4.*dg313) + 4.*dg213*dg333)*
      pow2(ginv33)) + ginv23*((ddg2233 + 2.*ddg2323)*ginv33 - 
     (4.*(dg233*dg323 + dg223*dg333) + 2.*pow2(dg233))*pow2(ginv33)) - 
  (dg111*dg233 + 2.*dg213*(dg113 + dg311) + dg211*(dg133 + 2.*dg313))*
   pow3(ginv13) - 2.*((dg222*dg223*ginv33 + ginv23*pow2(dg222))*
      pow2(ginv22) + (dg223*dg322 + dg222*(dg233 + dg323) + pow2(dg223))*
      pow3(ginv23) + dg233*dg333*pow3(ginv33))
;

dGfromgdu31
=
(ddg1311 - ((4.*dg112 + 2.*dg211)*dg311 + 4.*dg111*dg312)*ginv12 - 
     (dg212*dg311 + (2.*dg112 + dg211)*dg312)*ginv22 - 
     (dg311*(dg213 + 2.*dg312) + dg211*dg313 + 
        2.*(dg113*dg312 + dg112*dg313))*ginv23 - 
     2.*(dg113 + dg311)*dg313*ginv33 - 
     ginv13*(4.*(dg113*dg311 + dg111*dg313) + 2.*pow2(dg311)))*pow2(ginv11) \
+ (ddg1322 + ddg2312 - (2.*dg122*dg322 + 3.*(dg222*dg312 + dg212*dg322))*
      ginv22 - ((2.*dg213 + 4.*dg312)*dg322 + 
        2.*(dg223*dg312 + dg222*dg313 + dg123*dg322 + 
           (dg122 + dg212)*dg323))*ginv23 - 
     (dg313*(dg223 + 2.*dg322) + (dg213 + 2.*(dg123 + dg312))*dg323)*
      ginv33 - ginv13*(4.*(dg123*dg312 + dg112*dg323) + 
        2.*(dg213*dg312 + (dg122 + dg212)*dg313 + dg113*dg322 + 
           dg311*(dg223 + dg322) + dg211*dg323 + pow2(dg312))))*pow2(ginv12) \
+ (ddg1333 + ddg3313 - (dg233*dg312 + dg223*dg313 + 
        (dg213 + 2.*(dg123 + dg312))*dg323 + dg212*dg333)*ginv22 - 
     (2.*(dg233*dg313 + dg133*dg323 + (dg123 + dg213)*dg333) + 
        4.*(dg313*dg323 + dg312*dg333))*ginv23 - 
     (2.*dg133 + 6.*dg313)*dg333*ginv33)*pow2(ginv13) + 
  ginv11*(ddg3313*ginv33 + ginv22*
      (ddg2312 - (dg222*dg313 + dg213*dg322 + 
           2.*(dg312*(dg223 + dg322) + dg212*dg323))*ginv23 - 
        (dg223*dg313 + (dg213 + 2.*dg312)*dg323)*ginv33) + 
     ginv23*(ddg2313 + ddg3312 - 
        (dg313*(dg233 + 4.*dg323) + (dg213 + 2.*dg312)*dg333)*ginv33) + 
     ginv12*(2.*ddg1312 + ddg2311 - 
        (dg311*(4.*dg123 + 3.*dg213 + 6.*dg312) + 3.*dg211*dg313 + 
           6.*(dg113*dg312 + dg112*dg313) + 4.*dg111*dg323)*ginv13 - 
        (dg222*dg311 + (2.*dg122 + 6.*dg212)*dg312 + 
           (2.*dg112 + dg211)*dg322)*ginv22 - 
        (4.*dg312*dg313 + 2.*((dg123 + dg213)*dg313 + 
              (dg113 + dg311)*dg323))*ginv33 - 
        ginv23*((2.*dg123 + 4.*dg213)*dg312 + dg311*(dg223 + 2.*dg322) + 
           dg211*dg323 + 2.*(dg122*dg313 + dg113*dg322 + dg112*dg323) + 
           4.*(dg212*dg313 + pow2(dg312)))) + 
     ginv13*(2.*ddg1313 + ddg3311 - 
        ((4.*dg213 + 8.*dg312)*dg313 + dg311*(dg233 + 2.*dg323) + 
           dg211*dg333 + 2.*(dg133*dg312 + dg123*dg313 + dg113*dg323 + 
              dg112*dg333))*ginv23 - 
        ginv22*(dg223*dg311 + dg211*dg323 + 
           2.*((dg123 + dg213)*dg312 + dg212*dg313 + dg112*dg323 + 
              pow2(dg312))) - ginv33*
         (2.*(dg133*dg313 + (dg113 + dg311)*dg333) + 6.*pow2(dg313))) - 
     ((2.*dg122 + 3.*dg212)*dg311 + (6.*dg112 + 3.*dg211)*dg312 + 
        2.*dg111*dg322)*pow2(ginv12) - 
     (6.*dg113*dg313 + dg311*(2.*dg133 + 6.*dg313) + 2.*dg111*dg333)*
      pow2(ginv13) - (dg222*dg312 + dg212*dg322)*pow2(ginv22) - 
     (dg313*(dg223 + 2.*dg322) + dg213*dg323 + dg312*(dg233 + 2.*dg323) + 
        dg212*dg333)*pow2(ginv23) - 2.*dg313*dg333*pow2(ginv33)) + 
  ginv12*(ddg3323*ginv33 + ginv13*
      (2.*ddg1323 + ddg2313 + ddg3312 - 
        (dg222*dg313 + (2.*dg123 + dg213)*dg322 + 
           dg312*(4.*dg223 + 2.*dg322) + (2.*dg122 + 4.*dg212)*dg323)*
         ginv22 - ((4.*dg213 + 8.*dg312)*dg323 + 
           4.*(dg313*(dg223 + dg322) + dg123*dg323) + 
           2.*(dg233*dg312 + dg133*dg322 + (dg122 + dg212)*dg333))*ginv23 \
- (dg313*(dg233 + 8.*dg323) + (dg213 + 2.*dg312)*dg333 + 
           2.*(dg133*dg323 + dg123*dg333))*ginv33) + 
     ginv22*(ddg2322 - 2.*(dg223 + dg322)*dg323*ginv33 - 
        ginv23*(3.*(dg223*dg322 + dg222*dg323) + 2.*pow2(dg322))) + 
     ginv23*(ddg2323 + ddg3322 - 
        ginv33*(dg233*dg323 + (dg223 + 2.*dg322)*dg333 + 4.*pow2(dg323))) - 
     (dg311*(dg233 + 4.*dg323) + 
        4.*((dg123 + dg312)*dg313 + dg113*dg323) + dg211*dg333 + 
        2.*(dg133*dg312 + dg213*dg313 + dg112*dg333))*pow2(ginv13) - 
     (2.*dg223*dg323 + dg322*(dg233 + 4.*dg323) + dg222*dg333)*
      pow2(ginv23) - 2.*(dg222*dg322*pow2(ginv22) + 
        dg323*dg333*pow2(ginv33))) + 
  ginv13*(ddg3333*ginv33 + ginv23*
      (ddg2333 + ddg3323 - (2.*dg233 + 6.*dg323)*dg333*ginv33) + 
     ginv22*(ddg2323 - (4.*dg223*dg323 + dg322*(dg233 + 2.*dg323) + 
           dg222*dg333)*ginv23 - 
        ginv33*(dg233*dg323 + dg223*dg333 + 2.*pow2(dg323))) - 
     (dg223*dg322 + dg222*dg323)*pow2(ginv22) - 
     2.*((dg233*dg323 + (dg223 + dg322)*dg333 + pow2(dg323))*pow2(ginv23) + 
        pow2(dg333)*pow2(ginv33))) - 
  (dg222*dg311 + dg211*dg322 + 2.*((dg122 + dg212)*dg312 + dg112*dg322))*
   pow3(ginv12) - 2.*(dg111*dg311*pow3(ginv11) + 
     (dg133*dg313 + (dg113 + dg311)*dg333 + pow2(dg313))*pow3(ginv13))
;

dGfromgdu32
=
-((2.*dg111*dg311*ginv12 + (dg112*dg311 + dg111*dg312)*ginv22 + 
       (dg113*dg311 + dg111*dg313)*ginv23)*pow2(ginv11)) + 
  (ddg1312 + ddg2311 - (4.*dg311*dg312 + 
        2.*((dg123 + dg213)*dg311 + dg113*dg312 + 
           (dg112 + dg211)*dg313 + dg111*dg323))*ginv13 - 
     ((3.*dg122 + 6.*dg212)*dg312 + 3.*dg112*dg322 + 
        2.*(dg222*dg311 + dg211*dg322))*ginv22 - 
     ((dg123 + 2.*(dg213 + dg312))*dg313 + (dg113 + 2.*dg311)*dg323)*
      ginv33 - ginv23*(4.*(dg213*dg312 + dg212*dg313) + 
        2.*(dg123*dg312 + dg122*dg313 + dg113*dg322 + 
           dg311*(dg223 + dg322) + (dg112 + dg211)*dg323 + pow2(dg312))))*
   pow2(ginv12) - ((dg123*dg313 + dg312*(dg133 + 2.*dg313) + 
        (dg113 + 2.*dg311)*dg323 + dg112*dg333)*ginv22 + 
     2.*ginv23*(dg133*dg313 + (dg113 + dg311)*dg333 + pow2(dg313)))*
   pow2(ginv13) + (ddg2322 - 2.*(dg223 + dg322)*dg323*ginv33 - 
     ginv23*(4.*(dg223*dg322 + dg222*dg323) + 2.*pow2(dg322)))*pow2(ginv22) \
+ (ddg2333 + ddg3323 - (2.*dg233 + 6.*dg323)*dg333*ginv33)*pow2(ginv23) + 
  ginv11*(-(ginv13*((dg311*(dg123 + 2.*dg312) + 
             2.*(dg113*dg312 + dg112*dg313) + dg111*dg323)*ginv22 + 
          (4.*dg113*dg313 + dg311*(dg133 + 2.*dg313) + dg111*dg333)*ginv23)\
) + ginv12*(ddg1311 - ((dg122 + 2.*dg212)*dg311 + 
           (6.*dg112 + 2.*dg211)*dg312 + dg111*dg322)*ginv22 - 
        (dg311*(dg123 + 2.*(dg213 + dg312)) + 2.*dg211*dg313 + 
           4.*(dg113*dg312 + dg112*dg313) + dg111*dg323)*ginv23 - 
        2.*(dg113 + dg311)*dg313*ginv33 - 
        ginv13*(3.*(dg113*dg311 + dg111*dg313) + 2.*pow2(dg311))) + 
     ginv22*(ddg1312 - ((dg123 + 2.*dg312)*dg313 + dg113*dg323)*ginv33 - 
        ginv23*(dg122*dg313 + dg113*dg322 + 
           2.*((dg123 + dg213)*dg312 + dg212*dg313 + dg112*dg323 + 
              pow2(dg312)))) + 
     ginv23*(ddg1313 - ginv33*
         (dg133*dg313 + dg113*dg333 + 2.*pow2(dg313))) - 
     ((3.*dg112 + 2.*dg211)*dg311 + 3.*dg111*dg312)*pow2(ginv12) - 
     ((dg122 + 2.*dg212)*dg312 + dg112*dg322)*pow2(ginv22) - 
     (dg133*dg312 + (dg123 + 2.*(dg213 + dg312))*dg313 + dg113*dg323 + 
        dg112*dg333)*pow2(ginv23)) + 
  ginv13*(ginv23*(ddg1333 + ddg3313 - (2.*dg133 + 6.*dg313)*dg333*ginv33) + 
     ginv22*(ddg1323 + ddg3312 - 
        (dg133*dg322 + (4.*dg123 + 2.*dg213 + 8.*dg312)*dg323 + 
           dg122*dg333 + 2.*(dg233*dg312 + dg313*(dg223 + dg322) + 
              dg212*dg333))*ginv23 - 
        ((dg133 + 4.*dg313)*dg323 + (dg123 + 2.*dg312)*dg333)*ginv33) - 
     (dg123*dg322 + dg122*dg323 + 
        2.*(dg312*(dg223 + dg322) + dg212*dg323))*pow2(ginv22) - 
     (2.*(dg233*dg313 + dg133*dg323 + (dg123 + dg213)*dg333) + 
        4.*(dg313*dg323 + dg312*dg333))*pow2(ginv23)) + 
  ginv12*(ddg3313*ginv33 + ginv22*
      (ddg1322 + 2.*ddg2312 - (4.*(dg222*dg313 + dg213*dg322) + 
           3.*(dg123*dg322 + dg122*dg323) + 
           6.*(dg312*(dg223 + dg322) + dg212*dg323))*ginv23 - 
        ((2.*dg213 + 4.*dg312)*dg323 + 
           2.*(dg313*(dg223 + dg322) + dg123*dg323))*ginv33) + 
     ginv23*(ddg1323 + 2.*ddg2313 + ddg3312 - 
        (dg133*dg323 + dg313*(2.*dg233 + 8.*dg323) + 
           (dg123 + 2.*(dg213 + dg312))*dg333)*ginv33) + 
     ginv13*(ddg1313 + ddg3311 - 
        (8.*dg312*dg313 + 4.*
            ((dg123 + dg213)*dg313 + (dg113 + dg311)*dg323) + 
           2.*(dg233*dg311 + dg133*dg312 + (dg112 + dg211)*dg333))*ginv23 \
- ginv22*(dg122*dg313 + dg113*dg322 + 
           2.*(dg213*dg312 + dg212*dg313 + dg311*(dg223 + dg322) + 
              dg211*dg323) + 4.*(dg123*dg312 + dg112*dg323 + pow2(dg312))) \
- ginv33*(dg133*dg313 + (dg113 + 2.*dg311)*dg333 + 4.*pow2(dg313))) - 
     (2.*dg113*dg313 + dg311*(dg133 + 4.*dg313) + dg111*dg333)*
      pow2(ginv13) - (2.*dg122*dg322 + 4.*(dg222*dg312 + dg212*dg322))*
      pow2(ginv22) - (dg133*dg322 + 
        4.*(dg313*(dg223 + dg322) + (dg213 + dg312)*dg323) + 
        dg122*dg333 + 2.*(dg233*dg312 + dg123*dg323 + dg212*dg333))*
      pow2(ginv23) - 2.*dg313*dg333*pow2(ginv33)) + 
  ginv22*(ddg3323*ginv33 + ginv23*
      (2.*ddg2323 + ddg3322 - ginv33*
         (2.*(dg233*dg323 + (dg223 + dg322)*dg333) + 6.*pow2(dg323))) - 
     (6.*dg223*dg323 + dg322*(2.*dg233 + 6.*dg323) + 2.*dg222*dg333)*
      pow2(ginv23) - 2.*dg323*dg333*pow2(ginv33)) + 
  ginv23*(ddg3333*ginv33 - 2.*pow2(dg333)*pow2(ginv33)) - 
  ((dg122 + 2.*dg212)*dg311 + 2.*(dg112 + dg211)*dg312 + dg111*dg322)*
   pow3(ginv12) - 2.*(dg222*dg322*pow3(ginv22) + 
     (dg233*dg323 + (dg223 + dg322)*dg333 + pow2(dg323))*pow3(ginv23))
;

dGfromgdu33
=
-((2.*dg111*dg311*ginv13 + (dg112*dg311 + dg111*dg312)*ginv23 + 
       (dg113*dg311 + dg111*dg313)*ginv33)*pow2(ginv11)) - 
  (((dg122 + 2.*dg212)*dg311 + 2.*(dg112 + dg211)*dg312 + dg111*dg322)*
      ginv13 + (dg222*dg311 + dg211*dg322 + 
        2.*((dg122 + dg212)*dg312 + dg112*dg322))*ginv23 + 
     (dg223*dg311 + (dg123 + dg213)*dg312 + (dg122 + dg212)*dg313 + 
        dg113*dg322 + (dg112 + dg211)*dg323)*ginv33)*pow2(ginv12) + 
  (ddg1313 + ddg3311 - ((2.*dg213 + 8.*dg312)*dg313 + 
        dg311*(dg233 + 4.*dg323) + dg211*dg333 + 
        2.*(dg133*dg312 + dg123*dg313 + dg113*dg323 + dg112*dg333))*ginv23 \
- ginv22*(dg223*dg311 + (dg123 + dg213)*dg312 + dg212*dg313 + 
        (dg112 + dg211)*dg323 + 2.*pow2(dg312)) - 
     ginv33*(4.*dg311*dg333 + 3.*(dg133*dg313 + dg113*dg333) + 
        6.*pow2(dg313)))*pow2(ginv13) - 
  (2.*dg222*dg322*ginv23 + (dg223*dg322 + dg222*dg323)*ginv33)*
   pow2(ginv22) + (ddg2323 + ddg3322 - 
     ginv33*(4.*dg322*dg333 + 3.*(dg233*dg323 + dg223*dg333) + 
        6.*pow2(dg323)))*pow2(ginv23) + ddg3333*pow2(ginv33) + 
  ginv13*((ddg1333 + 2.*ddg3313)*ginv33 + 
     ginv22*(ddg2312 - (dg222*dg313 + (dg123 + dg213)*dg322 + 
           dg122*dg323 + 4.*(dg312*(dg223 + dg322) + dg212*dg323))*ginv23 \
- (dg312*(dg233 + 4.*dg323) + 2.*(dg223*dg313 + (dg123 + dg213)*dg323) + 
           dg212*dg333)*ginv33) + 
     ginv23*(ddg1323 + ddg2313 + 2.*ddg3312 - 
        (12.*dg313*dg323 + (3.*dg213 + 8.*dg312)*dg333 + 
           3.*(dg233*dg313 + dg133*dg323 + dg123*dg333))*ginv33) - 
     (dg222*dg312 + dg212*dg322)*pow2(ginv22) - 
     ((dg133 + 4.*dg313)*dg322 + (2.*dg213 + 8.*dg312)*dg323 + 
        dg122*dg333 + 2.*(dg233*dg312 + dg223*dg313 + dg123*dg323 + 
           dg212*dg333))*pow2(ginv23) - 
     (2.*dg133 + 8.*dg313)*dg333*pow2(ginv33)) + 
  ginv23*((ddg2333 + 2.*ddg3323)*ginv33 - 
     (2.*dg233 + 8.*dg323)*dg333*pow2(ginv33)) + 
  ginv12*((ddg1323 + ddg2313)*ginv33 - 
     ginv22*((2.*dg122*dg322 + 3.*(dg222*dg312 + dg212*dg322))*ginv23 + 
        (dg222*dg313 + (dg123 + dg213)*dg322 + dg122*dg323 + 
           2.*(dg223*dg312 + dg212*dg323))*ginv33) + 
     ginv23*(ddg1322 + ddg2312 - 
        (dg233*dg312 + dg133*dg322 + 
           4.*(dg313*(dg223 + dg322) + (dg123 + dg213 + dg312)*dg323) + 
           (dg122 + dg212)*dg333)*ginv33) + 
     ginv13*(ddg1312 + ddg2311 - 
        (dg222*dg311 + (dg122 + 4.*dg212)*dg312 + (dg112 + dg211)*dg322)*
         ginv22 - (dg133*dg312 + dg311*(dg233 + 4.*dg323) + 
           4.*((dg123 + dg213 + dg312)*dg313 + dg113*dg323) + 
           (dg112 + dg211)*dg333)*ginv33 - 
        ginv23*(2.*(dg223*dg311 + dg122*dg313 + dg113*dg322 + 
              dg211*dg323) + 4.*
            ((dg123 + dg213)*dg312 + dg212*dg313 + dg311*dg322 + 
              dg112*dg323 + pow2(dg312)))) - 
     (4.*dg311*dg312 + 2.*((dg123 + dg213)*dg311 + dg113*dg312 + 
           (dg112 + dg211)*dg313 + dg111*dg323))*pow2(ginv13) - 
     ((2.*dg213 + 4.*dg312)*dg322 + 
        2.*(dg223*dg312 + dg222*dg313 + dg123*dg322 + 
           (dg122 + dg212)*dg323))*pow2(ginv23) - 
     (dg133*dg323 + dg313*(dg233 + 4.*dg323) + (dg123 + dg213)*dg333)*
      pow2(ginv33)) + ginv11*(ddg1313*ginv33 - 
     ginv12*(((3.*dg112 + 2.*dg211)*dg311 + 3.*dg111*dg312)*ginv13 + 
        ((dg122 + dg212)*dg311 + (4.*dg112 + dg211)*dg312 + dg111*dg322)*
         ginv23 + ((dg123 + dg213)*dg311 + dg211*dg313 + 
           2.*(dg113*dg312 + dg112*dg313) + dg111*dg323)*ginv33) - 
     ginv22*(((dg122 + 2.*dg212)*dg312 + dg112*dg322)*ginv23 + 
        ((dg123 + dg213)*dg312 + dg212*dg313 + dg112*dg323)*ginv33) + 
     ginv13*(ddg1311 - (dg212*dg311 + (2.*dg112 + dg211)*dg312)*ginv22 - 
        ((dg123 + dg213)*dg311 + 4.*(dg113 + dg311)*dg312 + 
           (4.*dg112 + dg211)*dg313 + dg111*dg323)*ginv23 - 
        (6.*dg113*dg313 + dg311*(dg133 + 4.*dg313) + dg111*dg333)*ginv33) + 
     ginv23*(ddg1312 - (dg312*(dg133 + 4.*dg313) + 
           2.*((dg123 + dg213)*dg313 + dg113*dg323) + dg112*dg333)*ginv33) \
- (3.*(dg113*dg311 + dg111*dg313) + 2.*pow2(dg311))*pow2(ginv13) - 
     ((dg123 + dg213)*dg312 + (dg122 + dg212)*dg313 + dg113*dg322 + 
        dg112*dg323 + 2.*pow2(dg312))*pow2(ginv23) - 
     (dg133*dg313 + dg113*dg333 + 2.*pow2(dg313))*pow2(ginv33)) + 
  ginv22*(ddg2323*ginv33 + ginv23*
      (ddg2322 - (6.*dg223*dg323 + dg322*(dg233 + 4.*dg323) + dg222*dg333)*
         ginv33) - (3.*(dg223*dg322 + dg222*dg323) + 2.*pow2(dg322))*
      pow2(ginv23) - (dg233*dg323 + dg223*dg333 + 2.*pow2(dg323))*
      pow2(ginv33)) - (2.*dg113*dg313 + dg311*(dg133 + 4.*dg313) + 
     dg111*dg333)*pow3(ginv13) - 
  (2.*dg223*dg323 + dg322*(dg233 + 4.*dg323) + dg222*dg333)*pow3(ginv23) - 
  2.*pow2(dg333)*pow3(ginv33)
;

R11
=
gammado111*Gfromg1 + gammado112*Gfromg2 + gammado113*Gfromg3 + 
  (-0.5*ddg1111 + 3.*gamma111*gammado111 + 
     2.*(gamma211*gammado112 + gamma311*gammado113) + 
     gamma211*gammado211 + gamma311*gammado311)*ginv11 + 
  (-ddg1211 + 3.*(gamma112*gammado111 + gamma111*gammado112) + 
     2.*(gamma212*gammado112 + gamma312*gammado113 + 
        gamma211*gammado122 + gamma311*gammado123) + gamma212*gammado211 + 
     gamma211*gammado212 + gamma312*gammado311 + gamma311*gammado312)*ginv12 \
+ (-ddg1311 + 3.*(gamma113*gammado111 + gamma111*gammado113) + 
     2.*(gamma213*gammado112 + gamma313*gammado113 + 
        gamma211*gammado123 + gamma311*gammado133) + gamma213*gammado211 + 
     gamma211*gammado213 + gamma313*gammado311 + gamma311*gammado313)*ginv13 \
+ (-0.5*ddg2211 + 3.*gamma112*gammado112 + 
     2.*(gamma212*gammado122 + gamma312*gammado123) + 
     gamma212*gammado212 + gamma312*gammado312)*ginv22 + 
  (-ddg2311 + 3.*(gamma113*gammado112 + gamma112*gammado113) + 
     2.*(gamma213*gammado122 + (gamma212 + gamma313)*gammado123 + 
        gamma312*gammado133) + gamma213*gammado212 + gamma212*gammado213 + 
     gamma313*gammado312 + gamma312*gammado313)*ginv23 + 
  (-0.5*ddg3311 + 3.*gamma113*gammado113 + 
     2.*(gamma213*gammado123 + gamma313*gammado133) + 
     gamma213*gammado213 + gamma313*gammado313)*ginv33 + dG11*g11 + 
  dG12*g12 + dG13*g13
;

R12
=
(-0.5*ddg1112 + gamma112*gammado111 + (gamma111 + gamma212)*gammado112 + 
     gamma312*gammado113 + gamma111*gammado211 + 2.*gamma211*gammado212 + 
     gamma311*(gammado213 + gammado312))*ginv11 + 
  (-ddg1212 + gamma122*gammado111 + (2.*gamma112 + gamma222)*gammado112 + 
     gamma322*gammado113 + (gamma111 + gamma212)*gammado122 + 
     gamma112*gammado211 + (gamma111 + 2.*gamma212)*gammado212 + 
     2.*gamma211*gammado222 + gamma312*
      (gammado123 + gammado213 + gammado312) + 
     gamma311*(gammado223 + gammado322))*ginv12 + 
  (-ddg1312 + gamma123*gammado111 + (gamma113 + gamma223)*gammado112 + 
     (gamma112 + gamma323)*gammado113 + (gamma111 + gamma212)*gammado123 + 
     gamma312*gammado133 + gamma113*gammado211 + 
     (gamma111 + gamma313)*gammado213 + 
     2.*(gamma213*gammado212 + gamma211*gammado223) + 
     gamma313*gammado312 + gamma311*(gammado233 + gammado323))*ginv13 + 
  (-0.5*ddg2212 + gamma122*gammado112 + (gamma112 + gamma222)*gammado122 + 
     gamma322*gammado123 + gamma112*gammado212 + 2.*gamma212*gammado222 + 
     gamma312*(gammado223 + gammado322))*ginv22 + 
  (-ddg2312 + gamma123*gammado112 + gamma122*gammado113 + 
     (gamma113 + gamma223)*gammado122 + 
     (gamma112 + gamma222 + gamma323)*gammado123 + gamma322*gammado133 + 
     gamma113*gammado212 + gamma112*gammado213 + 
     2.*(gamma213*gammado222 + gamma212*gammado223) + 
     gamma313*(gammado223 + gammado322) + 
     gamma312*(gammado233 + gammado323))*ginv23 + 
  (-0.5*ddg3312 + gamma123*gammado113 + (gamma113 + gamma223)*gammado123 + 
     gamma323*gammado133 + gamma113*gammado213 + 2.*gamma213*gammado223 + 
     gamma313*(gammado233 + gammado323))*ginv33 + 
  0.5*((gammado112 + gammado211)*Gfromg1 + 
     (gammado122 + gammado212)*Gfromg2 + (gammado123 + gammado213)*Gfromg3 + 
     dG21*g11 + (dG11 + dG22)*g12 + dG23*g13 + 
     dG12*g22 + dG13*g23)
;

R13
=
(-0.5*ddg1113 + gamma113*gammado111 + gamma213*gammado112 + 
     (gamma111 + gamma313)*gammado113 + gamma111*gammado311 + 
     gamma211*(gammado213 + gammado312) + 2.*gamma311*gammado313)*ginv11 + 
  (-ddg1213 + gamma123*gammado111 + (gamma113 + gamma223)*gammado112 + 
     (gamma112 + gamma323)*gammado113 + gamma213*gammado122 + 
     (gamma111 + gamma313)*gammado123 + gamma112*gammado311 + 
     gamma111*gammado312 + gamma212*(gammado213 + gammado312) + 
     gamma211*(gammado223 + gammado322) + 
     2.*(gamma312*gammado313 + gamma311*gammado323))*ginv12 + 
  (-ddg1313 + gamma133*gammado111 + gamma233*gammado112 + 
     (2.*gamma113 + gamma333)*gammado113 + 
     (gamma111 + gamma313)*gammado133 + gamma113*gammado311 + 
     gamma213*(gammado123 + gammado213 + gammado312) + 
     (gamma111 + 2.*gamma313)*gammado313 + 
     gamma211*(gammado233 + gammado323) + 2.*gamma311*gammado333)*ginv13 + 
  (-0.5*ddg2213 + gamma123*gammado112 + gamma223*gammado122 + 
     (gamma112 + gamma323)*gammado123 + gamma112*gammado312 + 
     gamma212*(gammado223 + gammado322) + 2.*gamma312*gammado323)*ginv22 + 
  (-ddg2313 + gamma133*gammado112 + gamma123*gammado113 + 
     gamma233*gammado122 + (gamma113 + gamma223 + gamma333)*gammado123 + 
     (gamma112 + gamma323)*gammado133 + gamma113*gammado312 + 
     gamma112*gammado313 + gamma213*(gammado223 + gammado322) + 
     gamma212*(gammado233 + gammado323) + 
     2.*(gamma313*gammado323 + gamma312*gammado333))*ginv23 + 
  (-0.5*ddg3313 + gamma133*gammado113 + gamma233*gammado123 + 
     (gamma113 + gamma333)*gammado133 + gamma113*gammado313 + 
     gamma213*(gammado233 + gammado323) + 2.*gamma313*gammado333)*ginv33 + 
  0.5*((gammado113 + gammado311)*Gfromg1 + 
     (gammado123 + gammado312)*Gfromg2 + (gammado133 + gammado313)*Gfromg3 + 
     dG31*g11 + dG32*g12 + (dG11 + dG33)*g13 + 
     dG12*g23 + dG13*g33)
;

R22
=
gammado212*Gfromg1 + gammado222*Gfromg2 + gammado223*Gfromg3 + 
  (-0.5*ddg1122 + gamma112*(gammado112 + 2.*gammado211) + 
     3.*gamma212*gammado212 + gamma312*(2.*gammado213 + gammado312))*ginv11 \
+ (-ddg1222 + gamma122*(gammado112 + 2.*gammado211) + 
     gamma112*(gammado122 + 2.*gammado212) + 
     3.*(gamma222*gammado212 + gamma212*gammado222) + 
     2.*(gamma322*gammado213 + gamma312*gammado223) + 
     gamma322*gammado312 + gamma312*gammado322)*ginv12 + 
  (-ddg1322 + gamma123*(gammado112 + 2.*gammado211) + 
     gamma112*(gammado123 + 2.*gammado213) + 
     3.*(gamma223*gammado212 + gamma212*gammado223) + 
     2.*(gamma323*gammado213 + gamma312*gammado233) + 
     gamma323*gammado312 + gamma312*gammado323)*ginv13 + 
  (-0.5*ddg2222 + gamma122*(gammado122 + 2.*gammado212) + 
     3.*gamma222*gammado222 + gamma322*(2.*gammado223 + gammado322))*ginv22 \
+ (-ddg2322 + gamma123*(gammado122 + 2.*gammado212) + 
     gamma122*(gammado123 + 2.*gammado213) + 
     3.*(gamma223*gammado222 + gamma222*gammado223) + 
     2.*(gamma323*gammado223 + gamma322*gammado233) + 
     gamma323*gammado322 + gamma322*gammado323)*ginv23 + 
  (-0.5*ddg3322 + gamma123*(gammado123 + 2.*gammado213) + 
     3.*gamma223*gammado223 + gamma323*(2.*gammado233 + gammado323))*ginv33 \
+ dG21*g12 + dG22*g22 + dG23*g23
;

R23
=
(-0.5*ddg1123 + gamma113*gammado211 + gamma213*gammado212 + 
     (gamma212 + gamma313)*gammado213 + 
     gamma112*(gammado113 + gammado311) + gamma212*gammado312 + 
     2.*gamma312*gammado313)*ginv11 + 
  (-ddg1223 + gamma123*gammado211 + (gamma113 + gamma223)*gammado212 + 
     (gamma222 + gamma323)*gammado213 + gamma213*gammado222 + 
     (gamma212 + gamma313)*gammado223 + 
     gamma122*(gammado113 + gammado311) + gamma222*gammado312 + 
     gamma112*(gammado123 + gammado312) + gamma212*gammado322 + 
     2.*(gamma322*gammado313 + gamma312*gammado323))*ginv12 + 
  (-ddg1323 + gamma133*gammado211 + gamma233*gammado212 + 
     (gamma113 + gamma223 + gamma333)*gammado213 + gamma213*gammado223 + 
     (gamma212 + gamma313)*gammado233 + 
     gamma123*(gammado113 + gammado311) + gamma223*gammado312 + 
     gamma112*(gammado133 + gammado313) + gamma212*gammado323 + 
     2.*(gamma323*gammado313 + gamma312*gammado333))*ginv13 + 
  (-0.5*ddg2223 + gamma123*gammado212 + gamma223*gammado222 + 
     (gamma222 + gamma323)*gammado223 + 
     gamma122*(gammado123 + gammado312) + gamma222*gammado322 + 
     2.*gamma322*gammado323)*ginv22 + 
  (-ddg2323 + gamma133*gammado212 + gamma233*gammado222 + 
     (2.*gamma223 + gamma333)*gammado223 + 
     (gamma222 + gamma323)*gammado233 + 
     gamma123*(gammado123 + gammado213 + gammado312) + 
     gamma122*(gammado133 + gammado313) + gamma223*gammado322 + 
     (gamma222 + 2.*gamma323)*gammado323 + 2.*gamma322*gammado333)*ginv23 + 
  (-0.5*ddg3323 + gamma133*gammado213 + gamma233*gammado223 + 
     (gamma223 + gamma333)*gammado233 + 
     gamma123*(gammado133 + gammado313) + gamma223*gammado323 + 
     2.*gamma323*gammado333)*ginv33 + 
  0.5*((gammado213 + gammado312)*Gfromg1 + 
     (gammado223 + gammado322)*Gfromg2 + (gammado233 + gammado323)*Gfromg3 + 
     dG31*g12 + dG21*g13 + dG32*g22 + 
     (dG22 + dG33)*g23 + dG23*g33)
;

R33
=
gammado313*Gfromg1 + gammado323*Gfromg2 + gammado333*Gfromg3 + 
  (-0.5*ddg1133 + gamma113*(gammado113 + 2.*gammado311) + 
     gamma213*(gammado213 + 2.*gammado312) + 3.*gamma313*gammado313)*ginv11 \
+ (-ddg1233 + gamma123*(gammado113 + 2.*gammado311) + 
     gamma113*(gammado123 + 2.*gammado312) + 
     gamma223*(gammado213 + 2.*gammado312) + 
     gamma213*(gammado223 + 2.*gammado322) + 
     3.*(gamma323*gammado313 + gamma313*gammado323))*ginv12 + 
  (-ddg1333 + gamma133*(gammado113 + 2.*gammado311) + 
     gamma233*(gammado213 + 2.*gammado312) + 
     gamma113*(gammado133 + 2.*gammado313) + 
     gamma213*(gammado233 + 2.*gammado323) + 
     3.*(gamma333*gammado313 + gamma313*gammado333))*ginv13 + 
  (-0.5*ddg2233 + gamma123*(gammado123 + 2.*gammado312) + 
     gamma223*(gammado223 + 2.*gammado322) + 3.*gamma323*gammado323)*ginv22 \
+ (-ddg2333 + gamma133*(gammado123 + 2.*gammado312) + 
     gamma123*(gammado133 + 2.*gammado313) + 
     gamma233*(gammado223 + 2.*gammado322) + 
     gamma223*(gammado233 + 2.*gammado323) + 
     3.*(gamma333*gammado323 + gamma323*gammado333))*ginv23 + 
  (-0.5*ddg3333 + gamma133*(gammado133 + 2.*gammado313) + 
     gamma233*(gammado233 + 2.*gammado323) + 3.*gamma333*gammado333)*ginv33 \
+ dG31*g13 + dG32*g23 + dG33*g33
;

ff
=
chi
;

oochipsipower
=
1/chipsipower
;

f
=
oochipsipower*log(ff)
;

psim4
=
exp(-4.*f)
;

df1
=
(dchi1*oochipsipower)/chi
;

df2
=
(dchi2*oochipsipower)/chi
;

df3
=
(dchi3*oochipsipower)/chi
;

ddf11
=
(ddchi11*oochipsipower)/chi - chipsipower*pow2(df1)
;

ddf12
=
-(chipsipower*df1*df2) + (ddchi12*oochipsipower)/chi
;

ddf13
=
-(chipsipower*df1*df3) + (ddchi13*oochipsipower)/chi
;

ddf22
=
(ddchi22*oochipsipower)/chi - chipsipower*pow2(df2)
;

ddf23
=
-(chipsipower*df2*df3) + (ddchi23*oochipsipower)/chi
;

ddf33
=
(ddchi33*oochipsipower)/chi - chipsipower*pow2(df3)
;

cddf11
=
ddf11 - df1*gamma111 - df2*gamma211 - df3*gamma311
;

cddf12
=
ddf12 - df1*gamma112 - df2*gamma212 - df3*gamma312
;

cddf13
=
ddf13 - df1*gamma113 - df2*gamma213 - df3*gamma313
;

cddf22
=
ddf22 - df1*gamma122 - df2*gamma222 - df3*gamma322
;

cddf23
=
ddf23 - df1*gamma123 - df2*gamma223 - df3*gamma323
;

cddf33
=
ddf33 - df1*gamma133 - df2*gamma233 - df3*gamma333
;

trcddf
=
cddf11*ginv11 + cddf22*ginv22 + 
  2.*(cddf12*ginv12 + cddf13*ginv13 + cddf23*ginv23) + cddf33*ginv33
;

Rphi11
=
-2.*(cddf11 + trcddf*g11) + (4. - 4.*ginv11*g11)*pow2(df1) - 
  g11*(8.*(df1*(df2*ginv12 + df3*ginv13) + df2*df3*ginv23) + 
     4.*(ginv22*pow2(df2) + ginv33*pow2(df3)))
;

Rphi12
=
df1*df2*(4. - 8.*ginv12*g12) - 2.*(cddf12 + trcddf*g12) - 
  g12*(8.*df3*(df1*ginv13 + df2*ginv23) + 
     4.*(ginv11*pow2(df1) + ginv22*pow2(df2) + ginv33*pow2(df3)))
;

Rphi13
=
df1*(4.*df3 - 8.*df2*ginv12*g13) - 2.*(cddf13 + trcddf*g13) - 
  g13*(8.*df3*(df1*ginv13 + df2*ginv23) + 
     4.*(ginv11*pow2(df1) + ginv22*pow2(df2) + ginv33*pow2(df3)))
;

Rphi22
=
-2.*(cddf22 + trcddf*g22) + (4. - 4.*ginv22*g22)*pow2(df2) - 
  g22*(8.*(df1*(df2*ginv12 + df3*ginv13) + df2*df3*ginv23) + 
     4.*(ginv11*pow2(df1) + ginv33*pow2(df3)))
;

Rphi23
=
df2*(4.*df3 - 8.*df1*ginv12*g23) - 2.*(cddf23 + trcddf*g23) - 
  g23*(8.*df3*(df1*ginv13 + df2*ginv23) + 
     4.*(ginv11*pow2(df1) + ginv22*pow2(df2) + ginv33*pow2(df3)))
;

Rphi33
=
-2.*(cddf33 + trcddf*g33) - 
  g33*(8.*(df1*(df2*ginv12 + df3*ginv13) + df2*df3*ginv23) + 
     4.*(ginv11*pow2(df1) + ginv22*pow2(df2))) + 
  (4. - 4.*ginv33*g33)*pow2(df3)
;

Rf11
=
R11 + Rphi11
;

Rf12
=
R12 + Rphi12
;

Rf13
=
R13 + Rphi13
;

Rf22
=
R22 + Rphi22
;

Rf23
=
R23 + Rphi23
;

Rf33
=
R33 + Rphi33
;

Rhat
=
psim4*(ginv11*Rf11 + ginv22*Rf22 + 
    2.*(ginv12*Rf12 + ginv13*Rf13 + ginv23*Rf23) + ginv33*Rf33)
;

cdda11
=
dda11 - da2*gamma211 - da3*gamma311 + 
  2.*((da2*df1 + da1*df2)*ginv12 + (da3*df1 + da1*df3)*ginv13 + 
     da2*df2*ginv22 + (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33)*g11 \
+ da1*(-4.*df1 - gamma111 + 2.*df1*ginv11*g11)
;

cdda12
=
dda12 - 2.*(da2*df1 + da1*df2) - da1*gamma112 - da2*gamma212 - 
  da3*gamma312 + 2.*(da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
     (da3*df1 + da1*df3)*ginv13 + da2*df2*ginv22 + 
     (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33)*g12
;

cdda13
=
dda13 - 2.*(da3*df1 + da1*df3) - da1*gamma113 - da2*gamma213 - 
  da3*gamma313 + 2.*(da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
     (da3*df1 + da1*df3)*ginv13 + da2*df2*ginv22 + 
     (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33)*g13
;

cdda22
=
dda22 - da1*gamma122 - da3*gamma322 + 
  2.*(da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
     (da3*df1 + da1*df3)*ginv13 + (da3*df2 + da2*df3)*ginv23 + 
     da3*df3*ginv33)*g22 + 
  da2*(-4.*df2 - gamma222 + 2.*df2*ginv22*g22)
;

cdda23
=
dda23 - 2.*(da3*df2 + da2*df3) - da1*gamma123 - da2*gamma223 - 
  da3*gamma323 + 2.*(da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
     (da3*df1 + da1*df3)*ginv13 + da2*df2*ginv22 + 
     (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33)*g23
;

cdda33
=
dda33 - da1*gamma133 - da2*gamma233 + 
  2.*(da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
     (da3*df1 + da1*df3)*ginv13 + da2*df2*ginv22 + 
     (da3*df2 + da2*df3)*ginv23)*g33 + 
  da3*(-4.*df3 - gamma333 + 2.*df3*ginv33*g33)
;

dda12
=
dda12 - 2.*(da2*df1 + da1*df2) - da1*gamma112 - da2*gamma212 - 
  da3*gamma312 + 2.*(da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
     (da3*df1 + da1*df3)*ginv13 + da2*df2*ginv22 + 
     (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33)*g12
;

dda13
=
dda13 - 2.*(da3*df1 + da1*df3) - da1*gamma113 - da2*gamma213 - 
  da3*gamma313 + 2.*(da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
     (da3*df1 + da1*df3)*ginv13 + da2*df2*ginv22 + 
     (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33)*g13
;

dda23
=
dda23 - 2.*(da3*df2 + da2*df3) - da1*gamma123 - da2*gamma223 - 
  da3*gamma323 + 2.*(da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
     (da3*df1 + da1*df3)*ginv13 + da2*df2*ginv22 + 
     (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33)*g23
;

trcdda
=
(cdda11*ginv11 + (cdda12 + dda12)*ginv12 + (cdda13 + dda13)*ginv13 + 
    cdda22*ginv22 + (cdda23 + dda23)*ginv23 + cdda33*ginv33)*psim4
;

AA11
=
2.*(ginv23*A12*A13 + 
     A11*(ginv12*A12 + ginv13*A13)) + ginv11*pow2(A11) + 
  ginv22*pow2(A12) + ginv33*pow2(A13)
;

AA12
=
A12*(ginv11*A11 + ginv22*A22) + ginv33*A13*A23 + 
  ginv13*(A12*A13 + A11*A23) + 
  ginv23*(A13*A22 + A12*A23) + 
  ginv12*(A11*A22 + pow2(A12))
;

AA13
=
ginv22*A12*A23 + ginv12*(A12*A13 + A11*A23) + 
  A13*(ginv11*A11 + ginv33*A33) + 
  ginv23*(A13*A23 + A12*A33) + 
  ginv13*(A11*A33 + pow2(A13))
;

AA21
=
A12*(ginv11*A11 + ginv22*A22) + ginv33*A13*A23 + 
  ginv13*(A12*A13 + A11*A23) + 
  ginv23*(A13*A22 + A12*A23) + 
  ginv12*(A11*A22 + pow2(A12))
;

AA22
=
2.*(ginv23*A22*A23 + 
     A12*(ginv12*A22 + ginv13*A23)) + ginv11*pow2(A12) + 
  ginv22*pow2(A22) + ginv33*pow2(A23)
;

AA23
=
ginv11*A12*A13 + ginv12*(A13*A22 + A12*A23) + 
  A23*(ginv22*A22 + ginv33*A33) + 
  ginv13*(A13*A23 + A12*A33) + 
  ginv23*(A22*A33 + pow2(A23))
;

AA31
=
ginv22*A12*A23 + ginv12*(A12*A13 + A11*A23) + 
  A13*(ginv11*A11 + ginv33*A33) + 
  ginv23*(A13*A23 + A12*A33) + 
  ginv13*(A11*A33 + pow2(A13))
;

AA32
=
ginv11*A12*A13 + ginv12*(A13*A22 + A12*A23) + 
  A23*(ginv22*A22 + ginv33*A33) + 
  ginv13*(A13*A23 + A12*A33) + 
  ginv23*(A22*A33 + pow2(A23))
;

AA33
=
2.*(ginv23*A23*A33 + 
     A13*(ginv12*A23 + ginv13*A33)) + ginv11*pow2(A13) + 
  ginv22*pow2(A23) + ginv33*pow2(A33)
;

Ainv11
=
2.*(ginv11*(ginv12*A12 + ginv13*A13) + ginv12*ginv13*A23) + 
  A11*pow2(ginv11) + A22*pow2(ginv12) + A33*pow2(ginv13)
;

Ainv12
=
ginv11*(ginv12*A11 + ginv22*A12 + ginv23*A13) + 
  ginv12*(ginv13*A13 + ginv22*A22 + ginv23*A23) + 
  ginv13*(ginv22*A23 + ginv23*A33) + A12*pow2(ginv12)
;

Ainv13
=
ginv11*(ginv13*A11 + ginv23*A12 + ginv33*A13) + 
  ginv12*(ginv13*A12 + ginv23*A22 + ginv33*A23) + 
  ginv13*(ginv23*A23 + ginv33*A33) + A13*pow2(ginv13)
;

Ainv22
=
2.*(ginv12*(ginv22*A12 + ginv23*A13) + ginv22*ginv23*A23) + 
  A11*pow2(ginv12) + A22*pow2(ginv22) + A33*pow2(ginv23)
;

Ainv23
=
ginv13*(ginv22*A12 + ginv23*A13) + 
  ginv12*(ginv13*A11 + ginv23*A12 + ginv33*A13) + 
  ginv22*(ginv23*A22 + ginv33*A23) + ginv23*ginv33*A33 + 
  A23*pow2(ginv23)
;

Ainv33
=
2.*(ginv13*(ginv23*A12 + ginv33*A13) + ginv23*ginv33*A23) + 
  A11*pow2(ginv13) + A22*pow2(ginv23) + A33*pow2(ginv33)
;

cdA111
=
dA111 - 2.*(gamma111*A11 + gamma211*A12 + gamma311*A13)
;

cdA112
=
dA112 - gamma112*A11 - (gamma111 + gamma212)*A12 - 
  gamma312*A13 - gamma211*A22 - gamma311*A23
;

cdA113
=
dA113 - gamma113*A11 - gamma213*A12 - 
  (gamma111 + gamma313)*A13 - gamma211*A23 - gamma311*A33
;

cdA122
=
dA122 - 2.*(gamma112*A12 + gamma212*A22 + gamma312*A23)
;

cdA123
=
dA123 - gamma113*A12 - gamma112*A13 - gamma213*A22 - 
  (gamma212 + gamma313)*A23 - gamma312*A33
;

cdA133
=
dA133 - 2.*(gamma113*A13 + gamma213*A23 + gamma313*A33)
;

cdA211
=
dA211 - 2.*(gamma112*A11 + gamma212*A12 + gamma312*A13)
;

cdA212
=
dA212 - gamma122*A11 - (gamma112 + gamma222)*A12 - 
  gamma322*A13 - gamma212*A22 - gamma312*A23
;

cdA213
=
dA213 - gamma123*A11 - gamma223*A12 - 
  (gamma112 + gamma323)*A13 - gamma212*A23 - gamma312*A33
;

cdA222
=
dA222 - 2.*(gamma122*A12 + gamma222*A22 + gamma322*A23)
;

cdA223
=
dA223 - gamma123*A12 - gamma122*A13 - gamma223*A22 - 
  (gamma222 + gamma323)*A23 - gamma322*A33
;

cdA233
=
dA233 - 2.*(gamma123*A13 + gamma223*A23 + gamma323*A33)
;

cdA311
=
dA311 - 2.*(gamma113*A11 + gamma213*A12 + gamma313*A13)
;

cdA312
=
dA312 - gamma123*A11 - (gamma113 + gamma223)*A12 - 
  gamma323*A13 - gamma213*A22 - gamma313*A23
;

cdA313
=
dA313 - gamma133*A11 - gamma233*A12 - 
  (gamma113 + gamma333)*A13 - gamma213*A23 - gamma313*A33
;

cdA322
=
dA322 - 2.*(gamma123*A12 + gamma223*A22 + gamma323*A23)
;

cdA323
=
dA323 - gamma133*A12 - gamma123*A13 - gamma233*A22 - 
  (gamma223 + gamma333)*A23 - gamma323*A33
;

cdA333
=
dA333 - 2.*(gamma133*A13 + gamma233*A23 + gamma333*A33)
;

divbeta
=
db11 + db22 + db33
;

totdivbeta
=
0.66666666666666666667*divbeta
;

lieg11
=
dg111*beta1 + dg211*beta2 + dg311*beta3 + 
  (2.*db11 - totdivbeta)*g11 + 2.*(db12*g12 + db13*g13)
;

lieg12
=
dg112*beta1 + dg212*beta2 + dg312*beta3 + db21*g11 + 
  (db11 + db22 - totdivbeta)*g12 + db23*g13 + db12*g22 + 
  db13*g23
;

lieg13
=
dg113*beta1 + dg213*beta2 + dg313*beta3 + db31*g11 + 
  db32*g12 + (db11 + db33 - totdivbeta)*g13 + db12*g23 + 
  db13*g33
;

lieg22
=
dg122*beta1 + dg222*beta2 + dg322*beta3 - 
  totdivbeta*g22 + 2.*(db21*g12 + db22*g22 + db23*g23)
;

lieg23
=
dg123*beta1 + dg223*beta2 + dg323*beta3 + db31*g12 + 
  db21*g13 + db32*g22 + (db22 + db33 - totdivbeta)*g23 + 
  db23*g33
;

lieg33
=
dg133*beta1 + dg233*beta2 + dg333*beta3 - 
  totdivbeta*g33 + 2.*(db31*g13 + db32*g23 + db33*g33)
;

lieA11
=
(2.*db11 - totdivbeta)*A11 + 2.*(db12*A12 + db13*A13) + 
  dA111*beta1 + dA211*beta2 + dA311*beta3
;

lieA12
=
db21*A11 + (db11 + db22 - totdivbeta)*A12 + db23*A13 + 
  db12*A22 + db13*A23 + dA112*beta1 + dA212*beta2 + 
  dA312*beta3
;

lieA13
=
db31*A11 + db32*A12 + (db11 + db33 - totdivbeta)*A13 + 
  db12*A23 + db13*A33 + dA113*beta1 + dA213*beta2 + 
  dA313*beta3
;

lieA22
=
-(totdivbeta*A22) + 2.*(db21*A12 + db22*A22 + 
     db23*A23) + dA122*beta1 + dA222*beta2 + dA322*beta3
;

lieA23
=
db31*A12 + db21*A13 + db32*A22 + 
  (db22 + db33 - totdivbeta)*A23 + db23*A33 + dA123*beta1 + 
  dA223*beta2 + dA323*beta3
;

lieA33
=
-(totdivbeta*A33) + 2.*(db31*A13 + db32*A23 + 
     db33*A33) + dA133*beta1 + dA233*beta2 + dA333*beta3
;

betas
=
sdown1*beta1 + sdown2*beta2 + sdown3*beta3
;

Dbetas
=
(db11*sdown1 + db12*sdown2 + db13*sdown3)*sup1 + 
  (db21*sdown1 + db22*sdown2 + db23*sdown3)*sup2 + 
  (db31*sdown1 + db32*sdown2 + db33*sdown3)*sup3
;

Dalpha
=
da1*sup1 + da2*sup2 + da3*sup3
;

DKhat
=
dKhat1*sup1 + dKhat2*sup2 + dKhat3*sup3
;

DK
=
dK1*sup1 + dK2*sup2 + dK3*sup3
;

DTheta
=
dTheta1*sup1 + dTheta2*sup2 + dTheta3*sup3
;

Gams
=
sdown1*G1 + sdown2*G2 + sdown3*G3
;

DGams
=
(dG11*sdown1 + dG12*sdown2 + dG13*sdown3)*sup1 + 
  (dG21*sdown1 + dG22*sdown2 + dG23*sdown3)*sup2 + 
  (dG31*sdown1 + dG32*sdown2 + dG33*sdown3)*sup3
;

GamA1
=
qud11*G1 + qud12*G2 + qud13*G3
;

GamA2
=
qud21*G1 + qud22*G2 + qud23*G3
;

GamA3
=
qud31*G1 + qud32*G2 + qud33*G3
;

DGamA1
=
(dG11*qud11 + dG12*qud12 + dG13*qud13)*sup1 + 
  (dG21*qud11 + dG22*qud12 + dG23*qud13)*sup2 + 
  (dG31*qud11 + dG32*qud12 + dG33*qud13)*sup3
;

DGamA2
=
(dG11*qud21 + dG12*qud22 + dG13*qud23)*sup1 + 
  (dG21*qud21 + dG22*qud22 + dG23*qud23)*sup2 + 
  (dG31*qud21 + dG32*qud22 + dG33*qud23)*sup3
;

DGamA3
=
(dG11*qud31 + dG12*qud32 + dG13*qud33)*sup1 + 
  (dG21*qud31 + dG22*qud32 + dG23*qud33)*sup2 + 
  (dG31*qud31 + dG32*qud32 + dG33*qud33)*sup3
;

betaA1
=
qud11*beta1 + qud12*beta2 + qud13*beta3
;

betaA2
=
qud21*beta1 + qud22*beta2 + qud23*beta3
;

betaA3
=
qud31*beta1 + qud32*beta2 + qud33*beta3
;

DbetaA1
=
(db11*qud11 + db12*qud12 + db13*qud13)*sup1 + 
  (db21*qud11 + db22*qud12 + db23*qud13)*sup2 + 
  (db31*qud11 + db32*qud12 + db33*qud13)*sup3
;

DbetaA2
=
(db11*qud21 + db12*qud22 + db13*qud23)*sup1 + 
  (db21*qud21 + db22*qud22 + db23*qud23)*sup2 + 
  (db31*qud21 + db32*qud22 + db33*qud23)*sup3
;

DbetaA3
=
(db11*qud31 + db12*qud32 + db13*qud33)*sup1 + 
  (db21*qud31 + db22*qud32 + db23*qud33)*sup2 + 
  (db31*qud31 + db32*qud32 + db33*qud33)*sup3
;

lienKhat
=
-((DKhat + Khat/r)*sqrt(muL))
;

lienTheta
=
-DTheta - (kappa1*(2. + kappa2) + 1/r)*Theta
;

lienK
=
lienKhat + 2.*lienTheta
;

rKhat
=
lienKhat*alpha + dKhat1*beta1 + dKhat2*beta2 + 
  dKhat3*beta3
;

rGams
=
-(((db11*sdown1 + db12*sdown2)*beta1 + 
       (db21*sdown1 + db22*sdown2 + db23*sdown3)*beta2 + 
       db31*sdown1*beta3)*pow2(shiftdriver)) + 
  beta3*(2.*ddb231*sdown1*shiftdriver*beta2 + 
     sdown2*(2.*ddb132*shiftdriver*beta1 - db32*pow2(shiftdriver)) + 
     sdown3*(shiftdriver*(dG33 + 2.*ddb133*beta1) - 
        db33*pow2(shiftdriver))) + 
  sdown3*(db13*(db21*shiftdriver*beta2 - 
        beta1*pow2(shiftdriver)) + 
     shiftdriver*((db12*db23 + db13*(db11 + db33) + dG13)*beta1 + 
        (db23*db33 + dG23)*beta2 + 
        beta3*(db13*db31 + db23*db32 + pow2(db33)) + 
        ddb113*pow2(beta1))) + 
  shiftdriver*((dG22*sdown2 + db22*(db21*sdown1 + db23*sdown3) + 
        2.*ddb123*sdown3*beta1)*beta2 + 
     2.*((ddb232*sdown2 + ddb233*sdown3)*beta2*beta3 + 
        beta1*((ddb121*sdown1 + ddb122*sdown2)*beta2 + 
           ddb131*sdown1*beta3)) + 
     sdown2*((db13*db32 + dG12)*beta1 + 
        (db32*(db22 + db33) + dG32)*beta3 + 
        db12*((db11 + db22)*beta1 + db31*beta3) + 
        beta2*(db12*db21 + db23*db32 + pow2(db22))) + 
     (ddb111*sdown1 + ddb112*sdown2)*pow2(beta1) + 
     (ddb221*sdown1 + ddb222*sdown2 + ddb223*sdown3)*pow2(beta2) + 
     (ddb332*sdown2 + ddb333*sdown3)*pow2(beta3) + 
     sdown1*((db11*db21 + db23*db31 + dG21)*beta2 + 
        (db21*db32 + db31*(db11 + db33) + dG31)*beta3 + 
        beta1*(db12*db21 + db13*db31 + dG11 + pow2(db11)) + 
        ddb331*pow2(beta3)))
;

rTheta
=
lienTheta*alpha + dTheta1*beta1 + dTheta2*beta2 + 
  dTheta3*beta3
;

rACss
=
sup1*(2.*lieA13*sup3 + 1.3333333333333333333*dK1*alpha*chi + 
     sup2*(-(cdda12*psim4) + 2.*(lieA12 - AA12*alpha) + 
        0.66666666666666666667*trcdda*g12)) + 
  sup3*(2.*((psim4*Rf13*sup1 - AA23*sup2)*alpha + 
        sup2*(lieA23 + (-AA32 + psim4*Rf23)*alpha)) + 
     1.3333333333333333333*dK3*alpha*chi + 
     sup1*(-(dda13*psim4) - 2.*AA13*alpha + 
        0.66666666666666666667*trcdda*g13)) + 
  (lieA11 - cdda11*psim4 - 2.*AA11*alpha + 
     0.33333333333333333333*trcdda*g11)*pow2(sup1) + 
  (lieA22 - 2.*AA22*alpha + 0.33333333333333333333*trcdda*g22)*
   pow2(sup2) - psim4*((cdda23 + dda23)*sup2*sup3 + 
     sup1*(dda12*sup2 + cdda13*sup3) + cdda22*pow2(sup2)) + 
  (lieA33 - cdda33*psim4 - 2.*AA33*alpha + 
     0.33333333333333333333*trcdda*g33)*pow2(sup3) - 
  alpha*(sup2*(2.*AA21*sup1 + 
        0.66666666666666666667*Rhat*sup3*g23) + 
     0.33333333333333333333*(Rhat*g11*pow2(sup1) + 
        dGfromgdu11*qud11*pow2(chi))) + 
  alpha*(1.3333333333333333333*dK2*sup2*chi + 
     2.*(sup1*(-(AA31*sup3) + K*sup2*A12) + K*sup2*sup3*A23 - 
        DTheta*chi) + ginv11*
      (3.*dchi1*(sup1*A11 + sup2*A12 + sup3*A13) - 
        2.*(cdA111*sup1 + cdA112*sup2 + cdA113*sup3)*chi) + 
     ginv12*(3.*(sup1*(dchi2*A11 + dchi1*A12) + 
           dchi2*(sup2*A12 + sup3*A13) + 
           dchi1*(sup2*A22 + sup3*A23)) - 
        2.*((cdA112 + cdA211)*sup1 + (cdA122 + cdA212)*sup2 + 
           (cdA123 + cdA213)*sup3)*chi) + 
     ginv22*(3.*dchi2*(sup1*A12 + sup2*A22 + sup3*A23) - 
        2.*(cdA212*sup1 + cdA222*sup2 + cdA223*sup3)*chi) + 
     ginv13*(3.*(dchi3*(sup1*A11 + sup2*A12) + 
           (dchi1*sup1 + dchi3*sup3)*A13 + 
           dchi1*(sup2*A23 + sup3*A33)) - 
        2.*((cdA113 + cdA311)*sup1 + (cdA123 + cdA312)*sup2 + 
           (cdA133 + cdA313)*sup3)*chi) + 
     ginv23*(3.*(sup1*(dchi3*A12 + dchi2*A13) + 
           sup2*(dchi3*A22 + dchi2*A23) + 
           sup3*(dchi3*A23 + dchi2*A33)) - 
        2.*((cdA213 + cdA312)*sup1 + (cdA223 + cdA322)*sup2 + 
           (cdA233 + cdA323)*sup3)*chi) + 
     ginv33*(3.*dchi3*(sup1*A13 + sup2*A23 + sup3*A33) - 
        2.*(cdA313*sup1 + cdA323*sup2 + cdA333*sup3)*chi) - 
     0.66666666666666666667*Rhat*sup1*sup3*g13 + 
     psim4*(2.*Rf12*sup1*sup2 + Rf11*pow2(sup1) + Rf22*pow2(sup2) + 
        Rf33*pow2(sup3)) + K*(2.*sup1*sup3*A13 + 
        A11*pow2(sup1) + A22*pow2(sup2) + A33*pow2(sup3)) + 
     (0.33333333333333333333*dG11*qud11 - 
        sdown3*(Gfromg3*kappa1 + 0.66666666666666666667*dG33*sup3) + 
        sdown1*(0.66666666666666666667*dGfromgdu11*sup1 + 
           kappa1*G1) + kappa1*
         (-(Gfromg1*sdown1) - Gfromg2*sdown2 + sdown2*G2 + 
           sdown3*G3))*pow2(chi) + 
     0.33333333333333333333*(-(Rhat*
           (g22*pow2(sup2) + g33*pow2(sup3))) + 
        ((dG12 - dGfromgdu12)*qud12 + (dG13 - dGfromgdu13)*qud13 + 
           (dG21 - dGfromgdu21)*qud21 + (dG22 - dGfromgdu22)*qud22 + 
           (dG23 - dGfromgdu23)*qud23 + (dG31 - dGfromgdu31)*qud31 + 
           (dG32 - dGfromgdu32)*qud32 + (dG33 - dGfromgdu33)*qud33)*
         pow2(chi))) + 0.66666666666666666667*
   (sup2*(sup3*trcdda*g23 + 
        (-(dG21*sdown1) - dG22*sdown2 + dGfromgdu22*sdown2 + 
           dGfromgdu23*sdown3)*alpha*pow2(chi)) + 
     alpha*((-(sdown3*(dG13*sup1 + dG23*sup2)) + 
           (-(dG31*sdown1) - dG32*sdown2 + dGfromgdu32*sdown2 + 
              dGfromgdu33*sdown3)*sup3 + 
           sdown1*(dGfromgdu21*sup2 + dGfromgdu31*sup3))*pow2(chi) + 
        sup1*(-(Rhat*sup2*g12) + 
           (-(dG11*sdown1) - dG12*sdown2 + dGfromgdu12*sdown2 + 
              dGfromgdu13*sdown3)*pow2(chi))))
;

rACqq
=
-rACss + (Ainv22*lieg22 + 2.*(Ainv12*lieg12 + Ainv13*lieg13 + 
        Ainv23*lieg23) - (2.*Ainv22*A22 + 
        4.*(Ainv12*A12 + Ainv13*A13 + Ainv23*A23))*
      alpha + Ainv11*(lieg11 - 2.*A11*alpha) + 
     Ainv33*(lieg33 - 2.*A33*alpha))*chi
;

rGamA1
=
-(((dG11*qud11 + dG12*qud12 + dG13*qud13)*sup1 + 
       (dG22*qud12 + dG23*qud13)*sup2 + (dG32*qud12 + dG33*qud13)*sup3 + 
       qud11*(dG21*sup2 + dG31*sup3))*vbetaA) + 
  (dG11*qud11 + dG12*qud12 + dG13*qud13)*beta1 + 
  (dG21*qud11 + dG22*qud12 + dG23*qud13)*beta2 + 
  (dG31*qud11 + dG32*qud12 + dG33*qud13)*beta3 - 
  ((((db11*quu11 + db21*quu12 + db31*quu13)*sdown1 + 
          (db12*quu11 + db22*quu12 + db32*quu13)*sdown2 + 
          (db13*quu11 + db23*quu12)*sdown3)*shiftdriver)/vbetaA + 
     (0.66666666666666666667*dTheta1*quu11 + 
        (1.3333333333333333333*dKhat2 + 0.66666666666666666667*dTheta2)*
         quu12)*alpha + quu13*
      ((db33*sdown3*shiftdriver)/vbetaA + 
        1.3333333333333333333*dKhat3*alpha))/chi + 
  (2.3333333333333333333*((ddb121*qud11 + ddb122*qud12 + ddb123*qud13)*
         quu12 + (ddb131*qud11 + ddb132*qud12 + ddb133*qud13)*quu13) + 
     0.33333333333333333333*((ddb122*qud22 + ddb123*qud23 + 
           ddb131*qud31 + ddb132*qud32)*quu11 + 
        (ddb221*qud21 + ddb222*qud22 + ddb223*qud23 + ddb231*qud31 + 
           ddb232*qud32 + ddb233*qud33)*quu12 + 
        (ddb231*qud21 + ddb232*qud22 + ddb233*qud23 + ddb331*qud31 + 
           ddb332*qud32 + ddb333*qud33)*quu13) + 
     (ddb221*qud11 + ddb222*qud12 + ddb223*qud13)*quu22 + 
     2.*(ddb231*qud11 + ddb232*qud12 + ddb233*qud13)*quu23 + 
     (ddb331*qud11 + ddb332*qud12 + ddb333*qud13)*quu33 + 
     1.3333333333333333333*((ddb111*qud11 + ddb112*qud12)*quu11 + 
        (ddb132*quu13*sdown2 + ddb113*quu11*sdown3)*sup1 + 
        (ddb232*quu13*sdown2 + ddb123*quu11*sdown3)*sup2 + 
        (ddb332*quu13*sdown2 + ddb133*quu11*sdown3)*sup3 + 
        sdown2*((ddb112*quu11 + ddb122*quu12)*sup1 + 
           (ddb122*quu11 + ddb222*quu12)*sup2 + 
           (ddb132*quu11 + ddb232*quu12)*sup3) + 
        sdown1*((ddb121*quu12 + ddb131*quu13)*sup1 + 
           (ddb221*quu12 + ddb231*quu13)*sup2 + 
           (ddb231*quu12 + ddb331*quu13)*sup3) + 
        sdown3*((ddb123*quu12 + ddb133*quu13)*sup1 + 
           (ddb223*quu12 + ddb233*quu13)*sup2 + 
           (ddb233*quu12 + ddb333*quu13)*sup3)) + 
     (shiftdriver*((db11*qud11 + db12*qud12 + db13*qud13)*sup1 + 
          (db21*qud11 + db22*qud12 + db23*qud13)*sup2 + 
          (db31*qud11 + db32*qud12 + db33*qud13)*sup3))/vbetaA + 
     ((dG21*quu12 + dG31*quu13)*sdown1 + 
        (dG12*quu11 + dG22*quu12 + dG32*quu13)*sdown2 + 
        (dG13*quu11 + dG23*quu12 + dG33*quu13)*sdown3)*vbetaA - 
     0.66666666666666666667*dTheta3*quu13*alpha + 
     quu11*(0.33333333333333333333*(ddb121*qud21 + ddb133*qud33) + 
        dG11*sdown1*vbetaA + 1.3333333333333333333*
         (ddb113*qud13 + sdown1*(ddb111*sup1 + ddb121*sup2 + ddb131*sup3) - 
           dKhat1*alpha)))/chi
;

rGamA2
=
-(((dG11*qud21 + dG12*qud22 + dG13*qud23)*sup1 + 
       (dG22*qud22 + dG23*qud23)*sup2 + (dG32*qud22 + dG33*qud23)*sup3 + 
       qud21*(dG21*sup2 + dG31*sup3))*vbetaA) + 
  (dG11*qud21 + dG12*qud22 + dG13*qud23)*beta1 + 
  (dG21*qud21 + dG22*qud22 + dG23*qud23)*beta2 + 
  (dG31*qud21 + dG32*qud22 + dG33*qud23)*beta3 - 
  ((((db11*quu12 + db21*quu22 + db31*quu23)*sdown1 + 
          (db12*quu12 + db22*quu22 + db32*quu23)*sdown2 + 
          (db13*quu12 + db23*quu22)*sdown3)*shiftdriver)/vbetaA + 
     (0.66666666666666666667*dTheta1*quu12 + 
        (1.3333333333333333333*dKhat2 + 0.66666666666666666667*dTheta2)*
         quu22)*alpha + quu23*
      ((db33*sdown3*shiftdriver)/vbetaA + 
        1.3333333333333333333*dKhat3*alpha))/chi + 
  ((ddb111*qud21 + ddb112*qud22)*quu11 + 
     2.*(ddb131*qud21 + ddb132*qud22 + ddb133*qud23)*quu13 + 
     (1.3333333333333333333*ddb223*qud23 + 
        0.33333333333333333333*(ddb121*qud11 + ddb231*qud31))*quu22 + 
     2.3333333333333333333*((ddb121*qud21 + ddb122*qud22)*quu12 + 
        (ddb231*qud21 + ddb232*qud22)*quu23) + 
     0.33333333333333333333*((ddb112*qud12 + ddb113*qud13 + 
           ddb132*qud32 + ddb133*qud33)*quu12 + 
        (ddb122*qud12 + ddb123*qud13 + ddb232*qud32 + ddb233*qud33)*
         quu22 + (ddb132*qud12 + ddb133*qud13 + ddb332*qud32 + 
           ddb333*qud33)*quu23) + 
     (ddb331*qud21 + ddb332*qud22 + ddb333*qud23)*quu33 + 
     1.3333333333333333333*((ddb221*qud21 + ddb222*qud22)*quu22 + 
        (ddb132*quu23*sdown2 + ddb113*quu12*sdown3)*sup1 + 
        (ddb232*quu23*sdown2 + ddb123*quu12*sdown3)*sup2 + 
        (ddb332*quu23*sdown2 + ddb133*quu12*sdown3)*sup3 + 
        sdown2*((ddb112*quu12 + ddb122*quu22)*sup1 + 
           (ddb122*quu12 + ddb222*quu22)*sup2 + 
           (ddb132*quu12 + ddb232*quu22)*sup3) + 
        sdown1*((ddb121*quu22 + ddb131*quu23)*sup1 + 
           (ddb221*quu22 + ddb231*quu23)*sup2 + 
           (ddb231*quu22 + ddb331*quu23)*sup3) + 
        sdown3*((ddb123*quu22 + ddb133*quu23)*sup1 + 
           (ddb223*quu22 + ddb233*quu23)*sup2 + 
           (ddb233*quu22 + ddb333*quu23)*sup3)) + 
     qud23*(ddb113*quu11 + (db33*shiftdriver*sup3)/vbetaA) + 
     (shiftdriver*((db11*qud21 + db12*qud22 + db13*qud23)*sup1 + 
          (db21*qud21 + db22*qud22 + db23*qud23)*sup2 + 
          (db31*qud21 + db32*qud22)*sup3))/vbetaA + 
     ((dG21*quu22 + dG31*quu23)*sdown1 + 
        (dG12*quu12 + dG22*quu22 + dG32*quu23)*sdown2 + 
        (dG13*quu12 + dG23*quu22 + dG33*quu23)*sdown3)*vbetaA + 
     quu23*(2.3333333333333333333*ddb233*qud23 + 
        0.33333333333333333333*(ddb131*qud11 + ddb331*qud31) - 
        0.66666666666666666667*dTheta3*alpha) + 
     quu12*(2.3333333333333333333*ddb123*qud23 + 
        0.33333333333333333333*(ddb111*qud11 + ddb131*qud31) + 
        dG11*sdown1*vbetaA + 1.3333333333333333333*
         (sdown1*(ddb111*sup1 + ddb121*sup2 + ddb131*sup3) - 
           dKhat1*alpha)))/chi
;

rGamA3
=
-(((dG11*qud31 + dG12*qud32 + dG13*qud33)*sup1 + 
       (dG22*qud32 + dG23*qud33)*sup2 + (dG32*qud32 + dG33*qud33)*sup3 + 
       qud31*(dG21*sup2 + dG31*sup3))*vbetaA) + 
  (dG11*qud31 + dG12*qud32 + dG13*qud33)*beta1 + 
  (dG21*qud31 + dG22*qud32 + dG23*qud33)*beta2 + 
  (dG31*qud31 + dG32*qud32 + dG33*qud33)*beta3 - 
  ((((db11*quu13 + db21*quu23 + db31*quu33)*sdown1 + 
          (db12*quu13 + db22*quu23 + db32*quu33)*sdown2 + 
          (db13*quu13 + db23*quu23)*sdown3)*shiftdriver)/vbetaA + 
     (0.66666666666666666667*dTheta1*quu13 + 
        (1.3333333333333333333*dKhat2 + 0.66666666666666666667*dTheta2)*
         quu23)*alpha + quu33*
      ((db33*sdown3*shiftdriver)/vbetaA + 
        1.3333333333333333333*dKhat3*alpha))/chi + 
  ((ddb111*qud31 + ddb112*qud32)*quu11 + 
     2.*(ddb122*qud32 + ddb123*qud33)*quu12 + 
     (ddb222*qud32 + ddb223*qud33)*quu22 + 
     qud31*(2.*ddb121*quu12 + ddb221*quu22) + 
     (0.33333333333333333333*(ddb121*qud11 + ddb223*qud23) + 
        2.3333333333333333333*ddb231*qud31)*quu23 + 
     2.3333333333333333333*((ddb132*qud32 + ddb133*qud33)*quu13 + 
        (ddb232*qud32 + ddb233*qud33)*quu23) + 
     0.33333333333333333333*((ddb112*qud12 + ddb113*qud13 + 
           ddb121*qud21 + ddb122*qud22)*quu13 + 
        (ddb122*qud12 + ddb123*qud13 + ddb221*qud21 + ddb222*qud22)*
         quu23 + (ddb132*qud12 + ddb133*qud13 + ddb231*qud21 + 
           ddb232*qud22)*quu33) + 
     1.3333333333333333333*((ddb332*qud32 + ddb333*qud33)*quu33 + 
        (ddb132*quu33*sdown2 + ddb113*quu13*sdown3)*sup1 + 
        (ddb232*quu33*sdown2 + ddb123*quu13*sdown3)*sup2 + 
        (ddb332*quu33*sdown2 + ddb133*quu13*sdown3)*sup3 + 
        sdown2*((ddb112*quu13 + ddb122*quu23)*sup1 + 
           (ddb122*quu13 + ddb222*quu23)*sup2 + 
           (ddb132*quu13 + ddb232*quu23)*sup3) + 
        sdown1*((ddb121*quu23 + ddb131*quu33)*sup1 + 
           (ddb221*quu23 + ddb231*quu33)*sup2 + 
           (ddb231*quu23 + ddb331*quu33)*sup3) + 
        sdown3*((ddb123*quu23 + ddb133*quu33)*sup1 + 
           (ddb223*quu23 + ddb233*quu33)*sup2 + 
           (ddb233*quu23 + ddb333*quu33)*sup3)) + 
     qud33*(ddb113*quu11 + (db33*shiftdriver*sup3)/vbetaA) + 
     (shiftdriver*((db11*qud31 + db12*qud32 + db13*qud33)*sup1 + 
          (db21*qud31 + db22*qud32 + db23*qud33)*sup2 + 
          (db31*qud31 + db32*qud32)*sup3))/vbetaA + 
     ((dG21*quu23 + dG31*quu33)*sdown1 + 
        (dG12*quu13 + dG22*quu23 + dG32*quu33)*sdown2 + 
        (dG13*quu13 + dG23*quu23 + dG33*quu33)*sdown3)*vbetaA + 
     quu33*(0.33333333333333333333*(ddb131*qud11 + ddb233*qud23) + 
        1.3333333333333333333*ddb331*qud31 - 
        0.66666666666666666667*dTheta3*alpha) + 
     quu13*(0.33333333333333333333*(ddb111*qud11 + ddb123*qud23) + 
        ddb131*(2.3333333333333333333*qud31 + 
           1.3333333333333333333*sdown1*sup3) + dG11*sdown1*vbetaA + 
        1.3333333333333333333*(sdown1*(ddb111*sup1 + ddb121*sup2) - 
           dKhat1*alpha)))/chi
;

rACsA1
=
-2.*((AA12*qud21 + AA13*qud31)*sup1 + (AA22*qud21 + AA23*qud31)*sup2 + 
     (AA32*qud21 + AA33*qud31)*sup3)*alpha - 
  ((cdda12*qud21 + cdda13*qud31)*sup1 + 
     (cdda22*qud21 + cdda23*qud31)*sup2 + dda23*qud21*sup3)*chi + 
  (-(cdda33*qud31*sup3) + 0.66666666666666666667*dK1*qud11*alpha)*
   chi + sup1*(qud11*(lieA11 - 2.*AA11*alpha - cdda11*chi) + 
     qud21*(lieA12 + Rf12*alpha*chi) + 
     qud31*(lieA13 + Rf13*alpha*chi)) + 
  sup2*(lieA23*qud31 + qud11*(lieA12 - 2.*AA21*alpha - 
        dda12*chi + Rf12*alpha*chi) + 
     qud21*(lieA22 + Rf22*alpha*chi)) + 
  sup3*(qud11*(lieA13 - 2.*AA31*alpha - dda13*chi) + 
     qud21*(lieA23 + Rf23*alpha*chi) + 
     qud31*(lieA33 + alpha*(K*A33 + Rf33*chi)) - 
     0.5*dG33*qdd13*alpha*pow2(chi)) + 
  alpha*(K*(sup1*(qud11*A11 + qud21*A12 + qud31*A13) + 
        qud11*(sup2*A12 + sup3*A13) + qud21*sup3*A23 + 
        sup2*(qud21*A22 + qud31*A23)) + 
     (-(dTheta1*qud11) - dTheta2*qud21 - dTheta3*qud31 + 
        0.66666666666666666667*(dK2*qud21 + dK3*qud31) + qud31*Rf23*sup2 + 
        qud11*(Rf11*sup1 + Rf13*sup3))*chi + 
     ginv11*(1.5*dchi1*(qud11*A11 + qud21*A12 + qud31*A13) - 
        (cdA111*qud11 + cdA112*qud21 + cdA113*qud31)*chi) + 
     ginv12*(1.5*(qud11*(dchi2*A11 + dchi1*A12) + 
           dchi2*(qud21*A12 + qud31*A13) + 
           dchi1*(qud21*A22 + qud31*A23)) - 
        ((cdA112 + cdA211)*qud11 + (cdA122 + cdA212)*qud21 + 
           (cdA123 + cdA213)*qud31)*chi) + 
     ginv22*(1.5*dchi2*(qud11*A12 + qud21*A22 + qud31*A23) - 
        (cdA212*qud11 + cdA222*qud21 + cdA223*qud31)*chi) + 
     ginv13*(1.5*(dchi3*(qud11*A11 + qud21*A12) + 
           (dchi1*qud11 + dchi3*qud31)*A13 + 
           dchi1*(qud21*A23 + qud31*A33)) - 
        ((cdA113 + cdA311)*qud11 + (cdA123 + cdA312)*qud21 + 
           (cdA133 + cdA313)*qud31)*chi) + 
     ginv23*(1.5*(qud11*(dchi3*A12 + dchi2*A13) + 
           qud21*(dchi3*A22 + dchi2*A23) + 
           qud31*(dchi3*A23 + dchi2*A33)) - 
        ((cdA213 + cdA312)*qud11 + (cdA223 + cdA322)*qud21 + 
           (cdA233 + cdA323)*qud31)*chi) + 
     ginv33*(1.5*dchi3*(qud11*A13 + qud21*A23 + qud31*A33) - 
        (cdA313*qud11 + cdA323*qud21 + cdA333*qud31)*chi) + 
     0.5*((-(dG11*qdd11) - dG12*qdd12 + dGfromgdu12*qdd12 + 
           dGfromgdu13*qdd13)*sup1 + 
        (-(dG21*qdd11) - dG22*qdd12 + dGfromgdu22*qdd12 + 
           dGfromgdu23*qdd13)*sup2 - 
        qdd13*(Gfromg3*kappa1 + dG13*sup1 + dG23*sup2) + 
        (-(dG31*qdd11) - dG32*qdd12 + dGfromgdu32*qdd12 + 
           dGfromgdu33*qdd13)*sup3 + 
        qdd11*(dGfromgdu11*sup1 + dGfromgdu21*sup2 + dGfromgdu31*sup3) + 
        kappa1*(-(Gfromg1*qdd11) - Gfromg2*qdd12 + qdd11*G1 + 
           qdd12*G2 + qdd13*G3))*pow2(chi))
;

rACsA2
=
-2.*((AA12*qud22 + AA13*qud32)*sup1 + (AA22*qud22 + AA23*qud32)*sup2 + 
     (AA32*qud22 + AA33*qud32)*sup3)*alpha - 
  ((cdda12*qud22 + cdda13*qud32)*sup1 + 
     (cdda22*qud22 + cdda23*qud32)*sup2 + dda23*qud22*sup3)*chi + 
  (-(cdda33*qud32*sup3) + 0.66666666666666666667*dK1*qud12*alpha)*
   chi + sup1*(qud12*(lieA11 - 2.*AA11*alpha - cdda11*chi) + 
     qud22*(lieA12 + Rf12*alpha*chi) + 
     qud32*(lieA13 + Rf13*alpha*chi)) + 
  sup2*(lieA23*qud32 + qud12*(lieA12 - 2.*AA21*alpha - 
        dda12*chi + Rf12*alpha*chi) + 
     qud22*(lieA22 + Rf22*alpha*chi)) + 
  sup3*(qud12*(lieA13 - 2.*AA31*alpha - dda13*chi) + 
     qud22*(lieA23 + Rf23*alpha*chi) + 
     qud32*(lieA33 + alpha*(K*A33 + Rf33*chi)) - 
     0.5*dG33*qdd23*alpha*pow2(chi)) + 
  alpha*(K*(sup1*(qud12*A11 + qud22*A12 + qud32*A13) + 
        qud12*(sup2*A12 + sup3*A13) + qud22*sup3*A23 + 
        sup2*(qud22*A22 + qud32*A23)) + 
     (-(dTheta1*qud12) - dTheta2*qud22 - dTheta3*qud32 + 
        0.66666666666666666667*(dK2*qud22 + dK3*qud32) + qud32*Rf23*sup2 + 
        qud12*(Rf11*sup1 + Rf13*sup3))*chi + 
     ginv11*(1.5*dchi1*(qud12*A11 + qud22*A12 + qud32*A13) - 
        (cdA111*qud12 + cdA112*qud22 + cdA113*qud32)*chi) + 
     ginv12*(1.5*(qud12*(dchi2*A11 + dchi1*A12) + 
           dchi2*(qud22*A12 + qud32*A13) + 
           dchi1*(qud22*A22 + qud32*A23)) - 
        ((cdA112 + cdA211)*qud12 + (cdA122 + cdA212)*qud22 + 
           (cdA123 + cdA213)*qud32)*chi) + 
     ginv22*(1.5*dchi2*(qud12*A12 + qud22*A22 + qud32*A23) - 
        (cdA212*qud12 + cdA222*qud22 + cdA223*qud32)*chi) + 
     ginv13*(1.5*(dchi3*(qud12*A11 + qud22*A12) + 
           (dchi1*qud12 + dchi3*qud32)*A13 + 
           dchi1*(qud22*A23 + qud32*A33)) - 
        ((cdA113 + cdA311)*qud12 + (cdA123 + cdA312)*qud22 + 
           (cdA133 + cdA313)*qud32)*chi) + 
     ginv23*(1.5*(qud12*(dchi3*A12 + dchi2*A13) + 
           qud22*(dchi3*A22 + dchi2*A23) + 
           qud32*(dchi3*A23 + dchi2*A33)) - 
        ((cdA213 + cdA312)*qud12 + (cdA223 + cdA322)*qud22 + 
           (cdA233 + cdA323)*qud32)*chi) + 
     ginv33*(1.5*dchi3*(qud12*A13 + qud22*A23 + qud32*A33) - 
        (cdA313*qud12 + cdA323*qud22 + cdA333*qud32)*chi) + 
     0.5*((-(dG11*qdd12) - dG12*qdd22 + dGfromgdu12*qdd22 + 
           dGfromgdu13*qdd23)*sup1 + 
        (-(dG21*qdd12) - dG22*qdd22 + dGfromgdu22*qdd22 + 
           dGfromgdu23*qdd23)*sup2 - 
        qdd23*(Gfromg3*kappa1 + dG13*sup1 + dG23*sup2) + 
        (-(dG31*qdd12) - dG32*qdd22 + dGfromgdu32*qdd22 + 
           dGfromgdu33*qdd23)*sup3 + 
        qdd12*(dGfromgdu11*sup1 + dGfromgdu21*sup2 + dGfromgdu31*sup3) + 
        kappa1*(-(Gfromg1*qdd12) - Gfromg2*qdd22 + qdd12*G1 + 
           qdd22*G2 + qdd23*G3))*pow2(chi))
;

rACsA3
=
-2.*((AA12*qud23 + AA13*qud33)*sup1 + (AA22*qud23 + AA23*qud33)*sup2 + 
     (AA32*qud23 + AA33*qud33)*sup3)*alpha - 
  ((cdda12*qud23 + cdda13*qud33)*sup1 + 
     (cdda22*qud23 + cdda23*qud33)*sup2 + dda23*qud23*sup3)*chi + 
  (-(cdda33*qud33*sup3) + 0.66666666666666666667*dK1*qud13*alpha)*
   chi + sup1*(qud13*(lieA11 - 2.*AA11*alpha - cdda11*chi) + 
     qud23*(lieA12 + Rf12*alpha*chi) + 
     qud33*(lieA13 + Rf13*alpha*chi)) + 
  sup2*(lieA23*qud33 + qud13*(lieA12 - 2.*AA21*alpha - 
        dda12*chi + Rf12*alpha*chi) + 
     qud23*(lieA22 + Rf22*alpha*chi)) + 
  sup3*(qud13*(lieA13 - 2.*AA31*alpha - dda13*chi) + 
     qud23*(lieA23 + Rf23*alpha*chi) + 
     qud33*(lieA33 + alpha*(K*A33 + Rf33*chi)) - 
     0.5*dG33*qdd33*alpha*pow2(chi)) + 
  alpha*(K*(sup1*(qud13*A11 + qud23*A12 + qud33*A13) + 
        qud13*(sup2*A12 + sup3*A13) + qud23*sup3*A23 + 
        sup2*(qud23*A22 + qud33*A23)) + 
     (-(dTheta1*qud13) - dTheta2*qud23 - dTheta3*qud33 + 
        0.66666666666666666667*(dK2*qud23 + dK3*qud33) + qud33*Rf23*sup2 + 
        qud13*(Rf11*sup1 + Rf13*sup3))*chi + 
     ginv11*(1.5*dchi1*(qud13*A11 + qud23*A12 + qud33*A13) - 
        (cdA111*qud13 + cdA112*qud23 + cdA113*qud33)*chi) + 
     ginv12*(1.5*(qud13*(dchi2*A11 + dchi1*A12) + 
           dchi2*(qud23*A12 + qud33*A13) + 
           dchi1*(qud23*A22 + qud33*A23)) - 
        ((cdA112 + cdA211)*qud13 + (cdA122 + cdA212)*qud23 + 
           (cdA123 + cdA213)*qud33)*chi) + 
     ginv22*(1.5*dchi2*(qud13*A12 + qud23*A22 + qud33*A23) - 
        (cdA212*qud13 + cdA222*qud23 + cdA223*qud33)*chi) + 
     ginv13*(1.5*(dchi3*(qud13*A11 + qud23*A12) + 
           (dchi1*qud13 + dchi3*qud33)*A13 + 
           dchi1*(qud23*A23 + qud33*A33)) - 
        ((cdA113 + cdA311)*qud13 + (cdA123 + cdA312)*qud23 + 
           (cdA133 + cdA313)*qud33)*chi) + 
     ginv23*(1.5*(qud13*(dchi3*A12 + dchi2*A13) + 
           qud23*(dchi3*A22 + dchi2*A23) + 
           qud33*(dchi3*A23 + dchi2*A33)) - 
        ((cdA213 + cdA312)*qud13 + (cdA223 + cdA322)*qud23 + 
           (cdA233 + cdA323)*qud33)*chi) + 
     ginv33*(1.5*dchi3*(qud13*A13 + qud23*A23 + qud33*A33) - 
        (cdA313*qud13 + cdA323*qud23 + cdA333*qud33)*chi) + 
     0.5*((-(dG11*qdd13) - dG12*qdd23 + dGfromgdu12*qdd23 + 
           dGfromgdu13*qdd33)*sup1 + 
        (-(dG21*qdd13) - dG22*qdd23 + dGfromgdu22*qdd23 + 
           dGfromgdu23*qdd33)*sup2 - 
        qdd33*(Gfromg3*kappa1 + dG13*sup1 + dG23*sup2) + 
        (-(dG31*qdd13) - dG32*qdd23 + dGfromgdu32*qdd23 + 
           dGfromgdu33*qdd33)*sup3 + 
        qdd13*(dGfromgdu11*sup1 + dGfromgdu21*sup2 + dGfromgdu31*sup3) + 
        kappa1*(-(Gfromg1*qdd13) - Gfromg2*qdd23 + qdd13*G1 + 
           qdd23*G2 + qdd33*G3))*pow2(chi))
;

rACABTF11
=
2.*(lieA12*qPhysuudd1211 + lieA13*qPhysuudd1311 + 
     qPhysuudd2311*(lieA23 - cdA123*sup1*alpha)) + 
  qPhysuudd1111*(lieA11 + alpha*
      (-AA11 - cdA211*sup2 + sup2*(cdA112 - (0.5*dchi1*A12)/chi))) \
+ alpha*(qPhysuudd1111*(cdA113*sup3 + 
        0.66666666666666666667*K*A11) + 
     1.3333333333333333333*K*(qPhysuudd1211*A12 + 
        qPhysuudd1311*A13 + qPhysuudd2311*A23) + 
     qPhysuudd3311*(-(cdA233*sup2) + 0.66666666666666666667*K*A33) + 
     sup3*(cdA123*qPhysuudd1211 + 
        qPhysuudd1111*(-cdA311 + (0.5*dchi3*A11)/chi)) + 
     qPhysuudd2211*A22*(0.66666666666666666667*K + 
        (0.5*dchi3*sup3)/chi) + 
     sup2*(-2.*cdA213*qPhysuudd1311 - cdA223*qPhysuudd2311 + 
        cdA322*qPhysuudd2311 + cdA323*qPhysuudd3311 + 
        (0.5*dchi2*qPhysuudd1111*A11)/chi) + 
     (dchi3*(-0.5*qPhysuudd1311*sup2 + qPhysuudd1211*sup3)*A12 + 
        (-0.5*dchi3*qPhysuudd3311*sup1 + dchi2*qPhysuudd1311*sup2)*
         A13 + sup1*(-0.5*dchi2*qPhysuudd1211*A11 + 
           dchi1*qPhysuudd2311*A23) + 
        0.5*((-(dchi3*qPhysuudd2311*sup1) + dchi2*qPhysuudd1211*sup2)*
            A12 + dchi3*qPhysuudd1311*sup3*A13 - 
           (dchi1*qPhysuudd1211 + dchi3*qPhysuudd2311)*sup2*A22 + 
           sup1*((dchi1*qPhysuudd1211 - dchi2*qPhysuudd2211)*A12 + 
              (dchi1*qPhysuudd1311 - dchi2*qPhysuudd2311)*A13 + 
              dchi1*qPhysuudd2211*A22) - 
           (dchi3*qPhysuudd3311*sup2 + dchi1*qPhysuudd1211*sup3)*
            A23 + ((-(dchi1*qPhysuudd1311) + dchi2*qPhysuudd2311)*
               sup2 + (-(dchi2*qPhysuudd2211) + dchi3*qPhysuudd2311)*sup3\
)*A23 + qPhysuudd3311*(dchi1*sup1 + dchi2*sup2)*A33 - 
           sup3*((dchi1*qPhysuudd1111 + dchi2*qPhysuudd1211)*A13 + 
              (dchi1*qPhysuudd1311 + dchi2*qPhysuudd2311)*A33)))/
      chi) - cdda11*qPhysuudd1111*chi + 
  qPhysuudd1211*(((-cdA112 + cdA211)*sup1 + (cdA122 - cdA212)*sup2 + 
        (cdA213 - 2.*cdA312)*sup3)*alpha - cdda12*chi) + 
  qPhysuudd1311*(((-cdA113 + cdA311)*sup1 + (cdA123 + cdA312)*sup2 + 
        (cdA133 - cdA313)*sup3)*alpha - cdda13*chi) + 
  qPhysuudd2211*(lieA22 + (-AA22 - cdA122*sup1 + cdA212*sup1 + 
        cdA223*sup3 - cdA322*sup3)*alpha - cdda22*chi) + 
  qPhysuudd2311*(((cdA213 + cdA312)*sup1 + (cdA233 - cdA323)*sup3)*
      alpha - cdda23*chi) + 
  qPhysuudd3311*(lieA33 + (-AA33 - cdA133*sup1 + cdA313*sup1)*alpha - 
     cdda33*chi) - qPhysuudd1211*
   ((AA12 + AA21)*alpha + dda12*chi) - 
  qPhysuudd1311*(alpha*(AA13 + AA31 + 
        (0.5*dchi3*sup1*A11)/chi) + dda13*chi) - 
  qPhysuudd2311*((AA23 + AA32)*alpha + dda23*chi)
;

rACABTF12
=
2.*(lieA12*qPhysuudd1212 + lieA13*qPhysuudd1312 + 
     qPhysuudd2312*(lieA23 - cdA123*sup1*alpha)) + 
  qPhysuudd1112*(lieA11 + alpha*
      (-AA11 - cdA211*sup2 + sup2*(cdA112 - (0.5*dchi1*A12)/chi))) \
+ alpha*(qPhysuudd1112*(cdA113*sup3 + 
        0.66666666666666666667*K*A11) + 
     1.3333333333333333333*K*(qPhysuudd1212*A12 + 
        qPhysuudd1312*A13 + qPhysuudd2312*A23) + 
     qPhysuudd3312*(-(cdA233*sup2) + 0.66666666666666666667*K*A33) + 
     sup3*(cdA123*qPhysuudd1212 + 
        qPhysuudd1112*(-cdA311 + (0.5*dchi3*A11)/chi)) + 
     qPhysuudd2212*A22*(0.66666666666666666667*K + 
        (0.5*dchi3*sup3)/chi) + 
     sup2*(-2.*cdA213*qPhysuudd1312 - cdA223*qPhysuudd2312 + 
        cdA322*qPhysuudd2312 + cdA323*qPhysuudd3312 + 
        (0.5*dchi2*qPhysuudd1112*A11)/chi) + 
     (dchi3*(-0.5*qPhysuudd1312*sup2 + qPhysuudd1212*sup3)*A12 + 
        (-0.5*dchi3*qPhysuudd3312*sup1 + dchi2*qPhysuudd1312*sup2)*
         A13 + sup1*(-0.5*dchi2*qPhysuudd1212*A11 + 
           dchi1*qPhysuudd2312*A23) + 
        0.5*((-(dchi3*qPhysuudd2312*sup1) + dchi2*qPhysuudd1212*sup2)*
            A12 + dchi3*qPhysuudd1312*sup3*A13 - 
           (dchi1*qPhysuudd1212 + dchi3*qPhysuudd2312)*sup2*A22 + 
           sup1*((dchi1*qPhysuudd1212 - dchi2*qPhysuudd2212)*A12 + 
              (dchi1*qPhysuudd1312 - dchi2*qPhysuudd2312)*A13 + 
              dchi1*qPhysuudd2212*A22) - 
           (dchi3*qPhysuudd3312*sup2 + dchi1*qPhysuudd1212*sup3)*
            A23 + ((-(dchi1*qPhysuudd1312) + dchi2*qPhysuudd2312)*
               sup2 + (-(dchi2*qPhysuudd2212) + dchi3*qPhysuudd2312)*sup3\
)*A23 + qPhysuudd3312*(dchi1*sup1 + dchi2*sup2)*A33 - 
           sup3*((dchi1*qPhysuudd1112 + dchi2*qPhysuudd1212)*A13 + 
              (dchi1*qPhysuudd1312 + dchi2*qPhysuudd2312)*A33)))/
      chi) - cdda11*qPhysuudd1112*chi + 
  qPhysuudd1212*(((-cdA112 + cdA211)*sup1 + (cdA122 - cdA212)*sup2 + 
        (cdA213 - 2.*cdA312)*sup3)*alpha - cdda12*chi) + 
  qPhysuudd1312*(((-cdA113 + cdA311)*sup1 + (cdA123 + cdA312)*sup2 + 
        (cdA133 - cdA313)*sup3)*alpha - cdda13*chi) + 
  qPhysuudd2212*(lieA22 + (-AA22 - cdA122*sup1 + cdA212*sup1 + 
        cdA223*sup3 - cdA322*sup3)*alpha - cdda22*chi) + 
  qPhysuudd2312*(((cdA213 + cdA312)*sup1 + (cdA233 - cdA323)*sup3)*
      alpha - cdda23*chi) + 
  qPhysuudd3312*(lieA33 + (-AA33 - cdA133*sup1 + cdA313*sup1)*alpha - 
     cdda33*chi) - qPhysuudd1212*
   ((AA12 + AA21)*alpha + dda12*chi) - 
  qPhysuudd1312*(alpha*(AA13 + AA31 + 
        (0.5*dchi3*sup1*A11)/chi) + dda13*chi) - 
  qPhysuudd2312*((AA23 + AA32)*alpha + dda23*chi)
;

rACABTF13
=
2.*(lieA12*qPhysuudd1213 + lieA13*qPhysuudd1313 + 
     qPhysuudd2313*(lieA23 - cdA123*sup1*alpha)) + 
  qPhysuudd1113*(lieA11 + alpha*
      (-AA11 - cdA211*sup2 + sup2*(cdA112 - (0.5*dchi1*A12)/chi))) \
+ alpha*(qPhysuudd1113*(cdA113*sup3 + 
        0.66666666666666666667*K*A11) + 
     1.3333333333333333333*K*(qPhysuudd1213*A12 + 
        qPhysuudd1313*A13 + qPhysuudd2313*A23) + 
     qPhysuudd3313*(-(cdA233*sup2) + 0.66666666666666666667*K*A33) + 
     sup3*(cdA123*qPhysuudd1213 + 
        qPhysuudd1113*(-cdA311 + (0.5*dchi3*A11)/chi)) + 
     qPhysuudd2213*A22*(0.66666666666666666667*K + 
        (0.5*dchi3*sup3)/chi) + 
     sup2*(-2.*cdA213*qPhysuudd1313 - cdA223*qPhysuudd2313 + 
        cdA322*qPhysuudd2313 + cdA323*qPhysuudd3313 + 
        (0.5*dchi2*qPhysuudd1113*A11)/chi) + 
     (dchi3*(-0.5*qPhysuudd1313*sup2 + qPhysuudd1213*sup3)*A12 + 
        (-0.5*dchi3*qPhysuudd3313*sup1 + dchi2*qPhysuudd1313*sup2)*
         A13 + sup1*(-0.5*dchi2*qPhysuudd1213*A11 + 
           dchi1*qPhysuudd2313*A23) + 
        0.5*((-(dchi3*qPhysuudd2313*sup1) + dchi2*qPhysuudd1213*sup2)*
            A12 + dchi3*qPhysuudd1313*sup3*A13 - 
           (dchi1*qPhysuudd1213 + dchi3*qPhysuudd2313)*sup2*A22 + 
           sup1*((dchi1*qPhysuudd1213 - dchi2*qPhysuudd2213)*A12 + 
              (dchi1*qPhysuudd1313 - dchi2*qPhysuudd2313)*A13 + 
              dchi1*qPhysuudd2213*A22) - 
           (dchi3*qPhysuudd3313*sup2 + dchi1*qPhysuudd1213*sup3)*
            A23 + ((-(dchi1*qPhysuudd1313) + dchi2*qPhysuudd2313)*
               sup2 + (-(dchi2*qPhysuudd2213) + dchi3*qPhysuudd2313)*sup3\
)*A23 + qPhysuudd3313*(dchi1*sup1 + dchi2*sup2)*A33 - 
           sup3*((dchi1*qPhysuudd1113 + dchi2*qPhysuudd1213)*A13 + 
              (dchi1*qPhysuudd1313 + dchi2*qPhysuudd2313)*A33)))/
      chi) - cdda11*qPhysuudd1113*chi + 
  qPhysuudd1213*(((-cdA112 + cdA211)*sup1 + (cdA122 - cdA212)*sup2 + 
        (cdA213 - 2.*cdA312)*sup3)*alpha - cdda12*chi) + 
  qPhysuudd1313*(((-cdA113 + cdA311)*sup1 + (cdA123 + cdA312)*sup2 + 
        (cdA133 - cdA313)*sup3)*alpha - cdda13*chi) + 
  qPhysuudd2213*(lieA22 + (-AA22 - cdA122*sup1 + cdA212*sup1 + 
        cdA223*sup3 - cdA322*sup3)*alpha - cdda22*chi) + 
  qPhysuudd2313*(((cdA213 + cdA312)*sup1 + (cdA233 - cdA323)*sup3)*
      alpha - cdda23*chi) + 
  qPhysuudd3313*(lieA33 + (-AA33 - cdA133*sup1 + cdA313*sup1)*alpha - 
     cdda33*chi) - qPhysuudd1213*
   ((AA12 + AA21)*alpha + dda12*chi) - 
  qPhysuudd1313*(alpha*(AA13 + AA31 + 
        (0.5*dchi3*sup1*A11)/chi) + dda13*chi) - 
  qPhysuudd2313*((AA23 + AA32)*alpha + dda23*chi)
;

rACABTF22
=
2.*(lieA12*qPhysuudd1222 + lieA13*qPhysuudd1322 + 
     qPhysuudd2322*(lieA23 - cdA123*sup1*alpha)) + 
  qPhysuudd1122*(lieA11 + alpha*
      (-AA11 - cdA211*sup2 + sup2*(cdA112 - (0.5*dchi1*A12)/chi))) \
+ alpha*(qPhysuudd1122*(cdA113*sup3 + 
        0.66666666666666666667*K*A11) + 
     1.3333333333333333333*K*(qPhysuudd1222*A12 + 
        qPhysuudd1322*A13 + qPhysuudd2322*A23) + 
     qPhysuudd3322*(-(cdA233*sup2) + 0.66666666666666666667*K*A33) + 
     sup3*(cdA123*qPhysuudd1222 + 
        qPhysuudd1122*(-cdA311 + (0.5*dchi3*A11)/chi)) + 
     qPhysuudd2222*A22*(0.66666666666666666667*K + 
        (0.5*dchi3*sup3)/chi) + 
     sup2*(-2.*cdA213*qPhysuudd1322 - cdA223*qPhysuudd2322 + 
        cdA322*qPhysuudd2322 + cdA323*qPhysuudd3322 + 
        (0.5*dchi2*qPhysuudd1122*A11)/chi) + 
     (dchi3*(-0.5*qPhysuudd1322*sup2 + qPhysuudd1222*sup3)*A12 + 
        (-0.5*dchi3*qPhysuudd3322*sup1 + dchi2*qPhysuudd1322*sup2)*
         A13 + sup1*(-0.5*dchi2*qPhysuudd1222*A11 + 
           dchi1*qPhysuudd2322*A23) + 
        0.5*((-(dchi3*qPhysuudd2322*sup1) + dchi2*qPhysuudd1222*sup2)*
            A12 + dchi3*qPhysuudd1322*sup3*A13 - 
           (dchi1*qPhysuudd1222 + dchi3*qPhysuudd2322)*sup2*A22 + 
           sup1*((dchi1*qPhysuudd1222 - dchi2*qPhysuudd2222)*A12 + 
              (dchi1*qPhysuudd1322 - dchi2*qPhysuudd2322)*A13 + 
              dchi1*qPhysuudd2222*A22) - 
           (dchi3*qPhysuudd3322*sup2 + dchi1*qPhysuudd1222*sup3)*
            A23 + ((-(dchi1*qPhysuudd1322) + dchi2*qPhysuudd2322)*
               sup2 + (-(dchi2*qPhysuudd2222) + dchi3*qPhysuudd2322)*sup3\
)*A23 + qPhysuudd3322*(dchi1*sup1 + dchi2*sup2)*A33 - 
           sup3*((dchi1*qPhysuudd1122 + dchi2*qPhysuudd1222)*A13 + 
              (dchi1*qPhysuudd1322 + dchi2*qPhysuudd2322)*A33)))/
      chi) - cdda11*qPhysuudd1122*chi + 
  qPhysuudd1222*(((-cdA112 + cdA211)*sup1 + (cdA122 - cdA212)*sup2 + 
        (cdA213 - 2.*cdA312)*sup3)*alpha - cdda12*chi) + 
  qPhysuudd1322*(((-cdA113 + cdA311)*sup1 + (cdA123 + cdA312)*sup2 + 
        (cdA133 - cdA313)*sup3)*alpha - cdda13*chi) + 
  qPhysuudd2222*(lieA22 + (-AA22 - cdA122*sup1 + cdA212*sup1 + 
        cdA223*sup3 - cdA322*sup3)*alpha - cdda22*chi) + 
  qPhysuudd2322*(((cdA213 + cdA312)*sup1 + (cdA233 - cdA323)*sup3)*
      alpha - cdda23*chi) + 
  qPhysuudd3322*(lieA33 + (-AA33 - cdA133*sup1 + cdA313*sup1)*alpha - 
     cdda33*chi) - qPhysuudd1222*
   ((AA12 + AA21)*alpha + dda12*chi) - 
  qPhysuudd1322*(alpha*(AA13 + AA31 + 
        (0.5*dchi3*sup1*A11)/chi) + dda13*chi) - 
  qPhysuudd2322*((AA23 + AA32)*alpha + dda23*chi)
;

rACABTF23
=
2.*(lieA12*qPhysuudd1223 + lieA13*qPhysuudd1323 + 
     qPhysuudd2323*(lieA23 - cdA123*sup1*alpha)) + 
  qPhysuudd1123*(lieA11 + alpha*
      (-AA11 - cdA211*sup2 + sup2*(cdA112 - (0.5*dchi1*A12)/chi))) \
+ alpha*(qPhysuudd1123*(cdA113*sup3 + 
        0.66666666666666666667*K*A11) + 
     1.3333333333333333333*K*(qPhysuudd1223*A12 + 
        qPhysuudd1323*A13 + qPhysuudd2323*A23) + 
     qPhysuudd3323*(-(cdA233*sup2) + 0.66666666666666666667*K*A33) + 
     sup3*(cdA123*qPhysuudd1223 + 
        qPhysuudd1123*(-cdA311 + (0.5*dchi3*A11)/chi)) + 
     qPhysuudd2223*A22*(0.66666666666666666667*K + 
        (0.5*dchi3*sup3)/chi) + 
     sup2*(-2.*cdA213*qPhysuudd1323 - cdA223*qPhysuudd2323 + 
        cdA322*qPhysuudd2323 + cdA323*qPhysuudd3323 + 
        (0.5*dchi2*qPhysuudd1123*A11)/chi) + 
     (dchi3*(-0.5*qPhysuudd1323*sup2 + qPhysuudd1223*sup3)*A12 + 
        (-0.5*dchi3*qPhysuudd3323*sup1 + dchi2*qPhysuudd1323*sup2)*
         A13 + sup1*(-0.5*dchi2*qPhysuudd1223*A11 + 
           dchi1*qPhysuudd2323*A23) + 
        0.5*((-(dchi3*qPhysuudd2323*sup1) + dchi2*qPhysuudd1223*sup2)*
            A12 + dchi3*qPhysuudd1323*sup3*A13 - 
           (dchi1*qPhysuudd1223 + dchi3*qPhysuudd2323)*sup2*A22 + 
           sup1*((dchi1*qPhysuudd1223 - dchi2*qPhysuudd2223)*A12 + 
              (dchi1*qPhysuudd1323 - dchi2*qPhysuudd2323)*A13 + 
              dchi1*qPhysuudd2223*A22) - 
           (dchi3*qPhysuudd3323*sup2 + dchi1*qPhysuudd1223*sup3)*
            A23 + ((-(dchi1*qPhysuudd1323) + dchi2*qPhysuudd2323)*
               sup2 + (-(dchi2*qPhysuudd2223) + dchi3*qPhysuudd2323)*sup3\
)*A23 + qPhysuudd3323*(dchi1*sup1 + dchi2*sup2)*A33 - 
           sup3*((dchi1*qPhysuudd1123 + dchi2*qPhysuudd1223)*A13 + 
              (dchi1*qPhysuudd1323 + dchi2*qPhysuudd2323)*A33)))/
      chi) - cdda11*qPhysuudd1123*chi + 
  qPhysuudd1223*(((-cdA112 + cdA211)*sup1 + (cdA122 - cdA212)*sup2 + 
        (cdA213 - 2.*cdA312)*sup3)*alpha - cdda12*chi) + 
  qPhysuudd1323*(((-cdA113 + cdA311)*sup1 + (cdA123 + cdA312)*sup2 + 
        (cdA133 - cdA313)*sup3)*alpha - cdda13*chi) + 
  qPhysuudd2223*(lieA22 + (-AA22 - cdA122*sup1 + cdA212*sup1 + 
        cdA223*sup3 - cdA322*sup3)*alpha - cdda22*chi) + 
  qPhysuudd2323*(((cdA213 + cdA312)*sup1 + (cdA233 - cdA323)*sup3)*
      alpha - cdda23*chi) + 
  qPhysuudd3323*(lieA33 + (-AA33 - cdA133*sup1 + cdA313*sup1)*alpha - 
     cdda33*chi) - qPhysuudd1223*
   ((AA12 + AA21)*alpha + dda12*chi) - 
  qPhysuudd1323*(alpha*(AA13 + AA31 + 
        (0.5*dchi3*sup1*A11)/chi) + dda13*chi) - 
  qPhysuudd2323*((AA23 + AA32)*alpha + dda23*chi)
;

rACABTF33
=
2.*(lieA12*qPhysuudd1233 + lieA13*qPhysuudd1333 + 
     qPhysuudd2333*(lieA23 - cdA123*sup1*alpha)) + 
  qPhysuudd1133*(lieA11 + alpha*
      (-AA11 - cdA211*sup2 + sup2*(cdA112 - (0.5*dchi1*A12)/chi))) \
+ alpha*(qPhysuudd1133*(cdA113*sup3 + 
        0.66666666666666666667*K*A11) + 
     1.3333333333333333333*K*(qPhysuudd1233*A12 + 
        qPhysuudd1333*A13 + qPhysuudd2333*A23) + 
     qPhysuudd3333*(-(cdA233*sup2) + 0.66666666666666666667*K*A33) + 
     sup3*(cdA123*qPhysuudd1233 + 
        qPhysuudd1133*(-cdA311 + (0.5*dchi3*A11)/chi)) + 
     qPhysuudd2233*A22*(0.66666666666666666667*K + 
        (0.5*dchi3*sup3)/chi) + 
     sup2*(-2.*cdA213*qPhysuudd1333 - cdA223*qPhysuudd2333 + 
        cdA322*qPhysuudd2333 + cdA323*qPhysuudd3333 + 
        (0.5*dchi2*qPhysuudd1133*A11)/chi) + 
     (dchi3*(-0.5*qPhysuudd1333*sup2 + qPhysuudd1233*sup3)*A12 + 
        (-0.5*dchi3*qPhysuudd3333*sup1 + dchi2*qPhysuudd1333*sup2)*
         A13 + sup1*(-0.5*dchi2*qPhysuudd1233*A11 + 
           dchi1*qPhysuudd2333*A23) + 
        0.5*((-(dchi3*qPhysuudd2333*sup1) + dchi2*qPhysuudd1233*sup2)*
            A12 + dchi3*qPhysuudd1333*sup3*A13 - 
           (dchi1*qPhysuudd1233 + dchi3*qPhysuudd2333)*sup2*A22 + 
           sup1*((dchi1*qPhysuudd1233 - dchi2*qPhysuudd2233)*A12 + 
              (dchi1*qPhysuudd1333 - dchi2*qPhysuudd2333)*A13 + 
              dchi1*qPhysuudd2233*A22) - 
           (dchi3*qPhysuudd3333*sup2 + dchi1*qPhysuudd1233*sup3)*
            A23 + ((-(dchi1*qPhysuudd1333) + dchi2*qPhysuudd2333)*
               sup2 + (-(dchi2*qPhysuudd2233) + dchi3*qPhysuudd2333)*sup3\
)*A23 + qPhysuudd3333*(dchi1*sup1 + dchi2*sup2)*A33 - 
           sup3*((dchi1*qPhysuudd1133 + dchi2*qPhysuudd1233)*A13 + 
              (dchi1*qPhysuudd1333 + dchi2*qPhysuudd2333)*A33)))/
      chi) - cdda11*qPhysuudd1133*chi + 
  qPhysuudd1233*(((-cdA112 + cdA211)*sup1 + (cdA122 - cdA212)*sup2 + 
        (cdA213 - 2.*cdA312)*sup3)*alpha - cdda12*chi) + 
  qPhysuudd1333*(((-cdA113 + cdA311)*sup1 + (cdA123 + cdA312)*sup2 + 
        (cdA133 - cdA313)*sup3)*alpha - cdda13*chi) + 
  qPhysuudd2233*(lieA22 + (-AA22 - cdA122*sup1 + cdA212*sup1 + 
        cdA223*sup3 - cdA322*sup3)*alpha - cdda22*chi) + 
  qPhysuudd2333*(((cdA213 + cdA312)*sup1 + (cdA233 - cdA323)*sup3)*
      alpha - cdda23*chi) + 
  qPhysuudd3333*(lieA33 + (-AA33 - cdA133*sup1 + cdA313*sup1)*alpha - 
     cdda33*chi) - qPhysuudd1233*
   ((AA12 + AA21)*alpha + dda12*chi) - 
  qPhysuudd1333*(alpha*(AA13 + AA31 + 
        (0.5*dchi3*sup1*A11)/chi) + dda13*chi) - 
  qPhysuudd2333*((AA23 + AA32)*alpha + dda23*chi)
;


if (givehPsi0) { 

gADM11
=
g11/chi
;

gADM12
=
g12/chi
;

gADM13
=
g13/chi
;

gADM21
=
g12/chi
;

gADM22
=
g22/chi
;

gADM23
=
g23/chi
;

gADM31
=
g13/chi
;

gADM32
=
g23/chi
;

gADM33
=
g33/chi
;

vu1
=
-yp
;

vu2
=
xp
;

vu3
=
0
;

wu1
=
((-(ADMginv13*sup2) + ADMginv12*sup3)*vu1 + 
    (ADMginv13*sup1 - ADMginv11*sup3)*vu2 + 
    (-(ADMginv12*sup1) + ADMginv11*sup2)*vu3)/Power(chi,1.5)
;

wu2
=
((-(ADMginv23*sup2) + ADMginv22*sup3)*vu1 + 
    (ADMginv23*sup1 - ADMginv12*sup3)*vu2 + 
    (-(ADMginv22*sup1) + ADMginv12*sup2)*vu3)/Power(chi,1.5)
;

wu3
=
((-(ADMginv33*sup2) + ADMginv23*sup3)*vu1 + 
    (ADMginv33*sup1 - ADMginv13*sup3)*vu2 + 
    (-(ADMginv23*sup1) + ADMginv13*sup2)*vu3)/Power(chi,1.5)
;

sdotv
=
(gADM11*sup1 + gADM21*sup2 + gADM31*sup3)*vu1 + 
  (gADM12*sup1 + gADM22*sup2 + gADM32*sup3)*vu2 + 
  (gADM13*sup1 + gADM23*sup2 + gADM33*sup3)*vu3
;

vu1
=
-(sdotv*sup1) + vu1
;

vu2
=
-(sdotv*sup2) + vu2
;

vu3
=
-(sdotv*sup3) + vu3
;

vdotv
=
(gADM31*vu1 + (gADM23 + gADM32)*vu2)*vu3 + 
  vu1*((gADM12 + gADM21)*vu2 + gADM13*vu3) + gADM11*pow2(vu1) + 
  gADM22*pow2(vu2) + gADM33*pow2(vu3)
;

vu1
=
vu1/Sqrt(vdotv)
;

vu2
=
vu2/Sqrt(vdotv)
;

vu3
=
vu3/Sqrt(vdotv)
;

sdotw
=
(gADM11*sup1 + gADM21*sup2 + gADM31*sup3)*wu1 + 
  (gADM12*sup1 + gADM22*sup2 + gADM32*sup3)*wu2 + 
  (gADM13*sup1 + gADM23*sup2 + gADM33*sup3)*wu3
;

vdotw
=
(gADM11*vu1 + gADM21*vu2 + gADM31*vu3)*wu1 + 
  (gADM12*vu1 + gADM22*vu2 + gADM32*vu3)*wu2 + 
  (gADM13*vu1 + gADM23*vu2 + gADM33*vu3)*wu3
;

wu1
=
-(sdotw*sup1) - vdotw*vu1 + wu1
;

wu2
=
-(sdotw*sup2) - vdotw*vu2 + wu2
;

wu3
=
-(sdotw*sup3) - vdotw*vu3 + wu3
;

wdotw
=
(gADM31*wu1 + (gADM23 + gADM32)*wu2)*wu3 + 
  wu1*((gADM12 + gADM21)*wu2 + gADM13*wu3) + gADM11*pow2(wu1) + 
  gADM22*pow2(wu2) + gADM33*pow2(wu3)
;

wu1
=
wu1/Sqrt(wdotw)
;

wu2
=
wu2/Sqrt(wdotw)
;

wu3
=
wu3/Sqrt(wdotw)
;

vd1
=
gADM11*vu1 + gADM12*vu2 + gADM13*vu3
;

vd2
=
gADM21*vu1 + gADM22*vu2 + gADM23*vu3
;

vd3
=
gADM31*vu1 + gADM32*vu2 + gADM33*vu3
;

wd1
=
gADM11*wu1 + gADM12*wu2 + gADM13*wu3
;

wd2
=
gADM21*wu1 + gADM22*wu2 + gADM23*wu3
;

wd3
=
gADM31*wu1 + gADM32*wu2 + gADM33*wu3
;

RehPsi0
=
Power(2.7182818284590452354,pow2(hPsi0parb)*
    (2.*hPsi0parc*time - pow2(hPsi0parc) - pow2(time)))*hPsi0para
;

ImhPsi0
=
0
;

rACABTF11
=
rACABTF11 + alpha*chi*
   (2.*ImhPsi0*vd1*wd1 + RehPsi0*(pow2(vd1) - pow2(wd1)))
;

rACABTF12
=
rACABTF12 + (vd2*(RehPsi0*vd1 + ImhPsi0*wd1) + 
     (ImhPsi0*vd1 - RehPsi0*wd1)*wd2)*alpha*chi
;

rACABTF13
=
rACABTF13 + (vd3*(RehPsi0*vd1 + ImhPsi0*wd1) + 
     (ImhPsi0*vd1 - RehPsi0*wd1)*wd3)*alpha*chi
;

rACABTF22
=
rACABTF22 + alpha*chi*
   (2.*ImhPsi0*vd2*wd2 + RehPsi0*(pow2(vd2) - pow2(wd2)))
;

rACABTF23
=
rACABTF23 + (vd3*(RehPsi0*vd2 + ImhPsi0*wd2) + 
     (ImhPsi0*vd2 - RehPsi0*wd2)*wd3)*alpha*chi
;

rACABTF33
=
rACABTF33 + alpha*chi*
   (2.*ImhPsi0*vd3*wd3 + RehPsi0*(pow2(vd3) - pow2(wd3)))
;


 }  

rA11
=
rACABTF11 + 0.5*qdd11*rACqq + 2.*
   (qud11*rACsA1 + qud21*rACsA2 + qud31*rACsA3)*sdown1 + rACss*pow2(sdown1)
;

rA12
=
rACABTF12 + 0.5*qdd12*rACqq + (qud11*rACsA1 + qud21*rACsA2 + qud31*rACsA3)*
   sdown2 + sdown1*(qud12*rACsA1 + qud22*rACsA2 + qud32*rACsA3 + 
     rACss*sdown2)
;

rA13
=
rACABTF13 + 0.5*qdd13*rACqq + (qud11*rACsA1 + qud21*rACsA2 + qud31*rACsA3)*
   sdown3 + sdown1*(qud13*rACsA1 + qud23*rACsA2 + qud33*rACsA3 + 
     rACss*sdown3)
;

rA22
=
rACABTF22 + 0.5*qdd22*rACqq + 2.*
   (qud12*rACsA1 + qud22*rACsA2 + qud32*rACsA3)*sdown2 + rACss*pow2(sdown2)
;

rA23
=
rACABTF23 + 0.5*qdd23*rACqq + (qud12*rACsA1 + qud22*rACsA2 + qud32*rACsA3)*
   sdown3 + sdown2*(qud13*rACsA1 + qud23*rACsA2 + qud33*rACsA3 + 
     rACss*sdown3)
;

rA33
=
rACABTF33 + 0.5*qdd33*rACqq + 2.*
   (qud13*rACsA1 + qud23*rACsA2 + qud33*rACsA3)*sdown3 + rACss*pow2(sdown3)
;

rG1
=
qud11*rGamA1 + qud12*rGamA2 + qud13*rGamA3 + rGams*sup1
;

rG2
=
qud21*rGamA1 + qud22*rGamA2 + qud23*rGamA3 + rGams*sup2
;

rG3
=
qud31*rGamA1 + qud32*rGamA2 + qud33*rGamA3 + rGams*sup3
;

#if 0 
rG1 -= kappa1*(G1-Gfromg1);
rG2 -= kappa1*(G2-Gfromg2);
rG3 -= kappa1*(G3-Gfromg3);

rA11 -= kappa1*A11/r;
rA12 -= kappa1*A12/r;
rA13 -= kappa1*A13/r;
rA22 -= kappa1*A22/r;
rA23 -= kappa1*A23/r;
rA33 -= kappa1*A33/r;
#endif

#endif
}  /* function */
// f and tof are uper index
#ifdef fortran1
void decompose2p1_1
#endif	
#ifdef fortran2
void DECOMPOSE2P1_1
#endif
#ifdef fortran3
void decompose2p1_1_
#endif
(double & r,double & xp,double & yp,double & zp,double & chi,
		 double & g11,double & g12,double & g13,double & g22,double & g23,double & g33,
		 double & f1,double & f2,double & f3,double & tofs,double & tof1,double & tof2,double & tof3)
{
double ADMginv11;
double ADMginv12;
double ADMginv13;
double ADMginv22;
double ADMginv23;
double ADMginv33;
double detginv;
double ginv11;
double ginv12;
double ginv13;
double ginv22;
double ginv23;
double ginv33;
double modshatARG;
double oomodshat;
double qud11;
double qud12;
double qud13;
double qud21;
double qud22;
double qud23;
double qud31;
double qud32;
double qud33;
double sdown1;
double sdown2;
double sdown3;
double shat1;
double shat2;
double shat3;
double sup1;
double sup2;
double sup3;

shat1=xp/r;shat2=yp/r;shat3=zp/r; 

detginv
=
1/(2.*g12*g13*g23 - g33*pow2(g12) + g22*(g11*g33 - pow2(g13)) - 
    g11*pow2(g23))
;

ginv11
=
detginv*(g22*g33 - pow2(g23))
;

ginv12
=
detginv*(g13*g23 - g12*g33)
;

ginv13
=
detginv*(-(g13*g22) + g12*g23)
;

ginv22
=
detginv*(g11*g33 - pow2(g13))
;

ginv23
=
detginv*(g12*g13 - g11*g23)
;

ginv33
=
detginv*(g11*g22 - pow2(g12))
;

ADMginv11
=
chi*ginv11
;

ADMginv12
=
chi*ginv12
;

ADMginv13
=
chi*ginv13
;

ADMginv22
=
chi*ginv22
;

ADMginv23
=
chi*ginv23
;

ADMginv33
=
chi*ginv33
;

modshatARG
=
2.*(ADMginv23*shat2*shat3 + shat1*(ADMginv12*shat2 + ADMginv13*shat3)) + 
  ADMginv11*pow2(shat1) + ADMginv22*pow2(shat2) + ADMginv33*pow2(shat3)
;


if (modshatARG<0.00001) {                           
      printf("modshat is wrong (%e)\n",modshatARG);
      modshatARG = 0.00001;
    }oomodshat
=
1/sqrt(modshatARG)
;

sdown1
=
oomodshat*shat1
;

sdown2
=
oomodshat*shat2
;

sdown3
=
oomodshat*shat3
;

sup1
=
ADMginv11*sdown1 + ADMginv12*sdown2 + ADMginv13*sdown3
;

sup2
=
ADMginv12*sdown1 + ADMginv22*sdown2 + ADMginv23*sdown3
;

sup3
=
ADMginv13*sdown1 + ADMginv23*sdown2 + ADMginv33*sdown3
;

qud11
=
1. - sdown1*sup1
;

qud12
=
-(sdown2*sup1)
;

qud13
=
-(sdown3*sup1)
;

qud21
=
-(sdown1*sup2)
;

qud22
=
1. - sdown2*sup2
;

qud23
=
-(sdown3*sup2)
;

qud31
=
-(sdown1*sup3)
;

qud32
=
-(sdown2*sup3)
;

qud33
=
1. - sdown3*sup3
;

tofs
=
f1*sdown1 + f2*sdown2 + f3*sdown3
;

tof1
=
f1*qud11 + f2*qud12 + f3*qud13
;

tof2
=
f1*qud21 + f2*qud22 + f3*qud23
;

tof3
=
f1*qud31 + f2*qud32 + f3*qud33
;
} /* function */
// f and tof are lower index
#ifdef fortran1
void decompose2p1_2
#endif	
#ifdef fortran2
void DECOMPOSE2P1_2
#endif
#ifdef fortran3
void decompose2p1_2_
#endif
(double & r,double & xp,double & yp,double & zp,double & chi,
		 double & g11,double & g12,double & g13,double & g22,double & g23,double & g33,
		 double & f11,double & f12,double & f13,double & f22,double & f23,double & f33,
		 double & tofqq,double & tofss,double & tofs1,double & tofs2,double & tofs3,
		 double & tof11,double & tof12,double & tof13,double & tof22,double & tof23,double & tof33)
{
double ADMginv11;
double ADMginv12;
double ADMginv13;
double ADMginv22;
double ADMginv23;
double ADMginv33;
double detginv;
double ginv11;
double ginv12;
double ginv13;
double ginv22;
double ginv23;
double ginv33;
double modshatARG;
double oomodshat;
double qdd11;
double qdd12;
double qdd13;
double qdd22;
double qdd23;
double qdd33;
double qPhysuudd1111;
double qPhysuudd1112;
double qPhysuudd1113;
double qPhysuudd1122;
double qPhysuudd1123;
double qPhysuudd1133;
double qPhysuudd1211;
double qPhysuudd1212;
double qPhysuudd1213;
double qPhysuudd1222;
double qPhysuudd1223;
double qPhysuudd1233;
double qPhysuudd1311;
double qPhysuudd1312;
double qPhysuudd1313;
double qPhysuudd1322;
double qPhysuudd1323;
double qPhysuudd1333;
double qPhysuudd2211;
double qPhysuudd2212;
double qPhysuudd2213;
double qPhysuudd2222;
double qPhysuudd2223;
double qPhysuudd2233;
double qPhysuudd2311;
double qPhysuudd2312;
double qPhysuudd2313;
double qPhysuudd2322;
double qPhysuudd2323;
double qPhysuudd2333;
double qPhysuudd3311;
double qPhysuudd3312;
double qPhysuudd3313;
double qPhysuudd3322;
double qPhysuudd3323;
double qPhysuudd3333;
double qud11;
double qud12;
double qud13;
double qud21;
double qud22;
double qud23;
double qud31;
double qud32;
double qud33;
double quu11;
double quu12;
double quu13;
double quu22;
double quu23;
double quu33;
double sdown1;
double sdown2;
double sdown3;
double shat1;
double shat2;
double shat3;
double sup1;
double sup2;
double sup3;

shat1=xp/r;shat2=yp/r;shat3=zp/r; 

detginv
=
1/(2.*g12*g13*g23 - g33*pow2(g12) + g22*(g11*g33 - pow2(g13)) - 
    g11*pow2(g23))
;

ginv11
=
detginv*(g22*g33 - pow2(g23))
;

ginv12
=
detginv*(g13*g23 - g12*g33)
;

ginv13
=
detginv*(-(g13*g22) + g12*g23)
;

ginv22
=
detginv*(g11*g33 - pow2(g13))
;

ginv23
=
detginv*(g12*g13 - g11*g23)
;

ginv33
=
detginv*(g11*g22 - pow2(g12))
;

ADMginv11
=
chi*ginv11
;

ADMginv12
=
chi*ginv12
;

ADMginv13
=
chi*ginv13
;

ADMginv22
=
chi*ginv22
;

ADMginv23
=
chi*ginv23
;

ADMginv33
=
chi*ginv33
;

modshatARG
=
2.*(ADMginv23*shat2*shat3 + shat1*(ADMginv12*shat2 + ADMginv13*shat3)) + 
  ADMginv11*pow2(shat1) + ADMginv22*pow2(shat2) + ADMginv33*pow2(shat3)
;


if (modshatARG<0.00001) {                           
      printf("modshat is wrong (%e)\n",modshatARG);
      modshatARG = 0.00001;
    }oomodshat
=
1/sqrt(modshatARG)
;

sdown1
=
oomodshat*shat1
;

sdown2
=
oomodshat*shat2
;

sdown3
=
oomodshat*shat3
;

sup1
=
ADMginv11*sdown1 + ADMginv12*sdown2 + ADMginv13*sdown3
;

sup2
=
ADMginv12*sdown1 + ADMginv22*sdown2 + ADMginv23*sdown3
;

sup3
=
ADMginv13*sdown1 + ADMginv23*sdown2 + ADMginv33*sdown3
;

qud11
=
1. - sdown1*sup1
;

qud12
=
-(sdown2*sup1)
;

qud13
=
-(sdown3*sup1)
;

qud21
=
-(sdown1*sup2)
;

qud22
=
1. - sdown2*sup2
;

qud23
=
-(sdown3*sup2)
;

qud31
=
-(sdown1*sup3)
;

qud32
=
-(sdown2*sup3)
;

qud33
=
1. - sdown3*sup3
;

qdd11
=
g11/chi - pow2(sdown1)
;

qdd12
=
g12/chi - sdown1*sdown2
;

qdd13
=
g13/chi - sdown1*sdown3
;

qdd22
=
g22/chi - pow2(sdown2)
;

qdd23
=
g23/chi - sdown2*sdown3
;

qdd33
=
g33/chi - pow2(sdown3)
;

quu11
=
ADMginv11 - pow2(sup1)
;

quu12
=
ADMginv12 - sup1*sup2
;

quu13
=
ADMginv13 - sup1*sup3
;

quu22
=
ADMginv22 - pow2(sup2)
;

quu23
=
ADMginv23 - sup2*sup3
;

quu33
=
ADMginv33 - pow2(sup3)
;

qPhysuudd1111
=
-0.5*qdd11*quu11 + pow2(qud11)
;

qPhysuudd1112
=
qud11*qud12 - 0.5*qdd12*quu11
;

qPhysuudd1113
=
qud11*qud13 - 0.5*qdd13*quu11
;

qPhysuudd1122
=
-0.5*qdd22*quu11 + pow2(qud12)
;

qPhysuudd1123
=
qud12*qud13 - 0.5*qdd23*quu11
;

qPhysuudd1133
=
-0.5*qdd33*quu11 + pow2(qud13)
;

qPhysuudd1211
=
qud11*qud21 - 0.5*qdd11*quu12
;

qPhysuudd1212
=
0.5*(qud12*qud21 + qud11*qud22 - qdd12*quu12)
;

qPhysuudd1213
=
0.5*(qud13*qud21 + qud11*qud23 - qdd13*quu12)
;

qPhysuudd1222
=
qud12*qud22 - 0.5*qdd22*quu12
;

qPhysuudd1223
=
0.5*(qud13*qud22 + qud12*qud23 - qdd23*quu12)
;

qPhysuudd1233
=
qud13*qud23 - 0.5*qdd33*quu12
;

qPhysuudd1311
=
qud11*qud31 - 0.5*qdd11*quu13
;

qPhysuudd1312
=
0.5*(qud12*qud31 + qud11*qud32 - qdd12*quu13)
;

qPhysuudd1313
=
0.5*(qud13*qud31 + qud11*qud33 - qdd13*quu13)
;

qPhysuudd1322
=
qud12*qud32 - 0.5*qdd22*quu13
;

qPhysuudd1323
=
0.5*(qud13*qud32 + qud12*qud33 - qdd23*quu13)
;

qPhysuudd1333
=
qud13*qud33 - 0.5*qdd33*quu13
;

qPhysuudd2211
=
-0.5*qdd11*quu22 + pow2(qud21)
;

qPhysuudd2212
=
qud21*qud22 - 0.5*qdd12*quu22
;

qPhysuudd2213
=
qud21*qud23 - 0.5*qdd13*quu22
;

qPhysuudd2222
=
-0.5*qdd22*quu22 + pow2(qud22)
;

qPhysuudd2223
=
qud22*qud23 - 0.5*qdd23*quu22
;

qPhysuudd2233
=
-0.5*qdd33*quu22 + pow2(qud23)
;

qPhysuudd2311
=
qud21*qud31 - 0.5*qdd11*quu23
;

qPhysuudd2312
=
0.5*(qud22*qud31 + qud21*qud32 - qdd12*quu23)
;

qPhysuudd2313
=
0.5*(qud23*qud31 + qud21*qud33 - qdd13*quu23)
;

qPhysuudd2322
=
qud22*qud32 - 0.5*qdd22*quu23
;

qPhysuudd2323
=
0.5*(qud23*qud32 + qud22*qud33 - qdd23*quu23)
;

qPhysuudd2333
=
qud23*qud33 - 0.5*qdd33*quu23
;

qPhysuudd3311
=
-0.5*qdd11*quu33 + pow2(qud31)
;

qPhysuudd3312
=
qud31*qud32 - 0.5*qdd12*quu33
;

qPhysuudd3313
=
qud31*qud33 - 0.5*qdd13*quu33
;

qPhysuudd3322
=
-0.5*qdd22*quu33 + pow2(qud32)
;

qPhysuudd3323
=
qud32*qud33 - 0.5*qdd23*quu33
;

qPhysuudd3333
=
-0.5*qdd33*quu33 + pow2(qud33)
;

tofss
=
2.*(f23*sup2*sup3 + sup1*(f12*sup2 + f13*sup3)) + f11*pow2(sup1) + 
  f22*pow2(sup2) + f33*pow2(sup3)
;

tofqq
=
f12*quu12 + f13*quu13 + f23*quu23 + 0.5*(f11*quu11 + f22*quu22 + f33*quu33)
;

tofs1
=
(f11*qud11 + f12*qud21 + f13*qud31)*sup1 + 
  (f12*qud11 + f22*qud21 + f23*qud31)*sup2 + 
  (f13*qud11 + f23*qud21 + f33*qud31)*sup3
;

tofs2
=
(f11*qud12 + f12*qud22 + f13*qud32)*sup1 + 
  (f12*qud12 + f22*qud22 + f23*qud32)*sup2 + 
  (f13*qud12 + f23*qud22 + f33*qud32)*sup3
;

tofs3
=
(f11*qud13 + f12*qud23 + f13*qud33)*sup1 + 
  (f12*qud13 + f22*qud23 + f23*qud33)*sup2 + 
  (f13*qud13 + f23*qud23 + f33*qud33)*sup3
;

tof11
=
f11*qPhysuudd1111 + f22*qPhysuudd2211 + 
  2.*(f12*qPhysuudd1211 + f13*qPhysuudd1311 + f23*qPhysuudd2311) + 
  f33*qPhysuudd3311
;

tof12
=
f11*qPhysuudd1112 + f22*qPhysuudd2212 + 
  2.*(f12*qPhysuudd1212 + f13*qPhysuudd1312 + f23*qPhysuudd2312) + 
  f33*qPhysuudd3312
;

tof13
=
f11*qPhysuudd1113 + f22*qPhysuudd2213 + 
  2.*(f12*qPhysuudd1213 + f13*qPhysuudd1313 + f23*qPhysuudd2313) + 
  f33*qPhysuudd3313
;

tof22
=
f11*qPhysuudd1122 + f22*qPhysuudd2222 + 
  2.*(f12*qPhysuudd1222 + f13*qPhysuudd1322 + f23*qPhysuudd2322) + 
  f33*qPhysuudd3322
;

tof23
=
f11*qPhysuudd1123 + f22*qPhysuudd2223 + 
  2.*(f12*qPhysuudd1223 + f13*qPhysuudd1323 + f23*qPhysuudd2323) + 
  f33*qPhysuudd3323
;

tof33
=
f11*qPhysuudd1133 + f22*qPhysuudd2233 + 
  2.*(f12*qPhysuudd1233 + f13*qPhysuudd1333 + f23*qPhysuudd2333) + 
  f33*qPhysuudd3333
;
} /*function */
// f and tof are uper index
#ifdef fortran1
void compose2p1_1
#endif	
#ifdef fortran2
void COMPOSE2P1_1
#endif
#ifdef fortran3
void compose2p1_1_
#endif
(double & r,double & xp,double & yp,double & zp,double & chi,
		 double & g11,double & g12,double & g13,double & g22,double & g23,double & g33,
		 double & f1,double & f2,double & f3,double & tofs,double & tof1,double & tof2,double & tof3)
{
double ADMginv11;
double ADMginv12;
double ADMginv13;
double ADMginv22;
double ADMginv23;
double ADMginv33;
double detginv;
double ginv11;
double ginv12;
double ginv13;
double ginv22;
double ginv23;
double ginv33;
double modshatARG;
double oomodshat;
double qud11;
double qud12;
double qud13;
double qud21;
double qud22;
double qud23;
double qud31;
double qud32;
double qud33;
double sdown1;
double sdown2;
double sdown3;
double shat1;
double shat2;
double shat3;
double sup1;
double sup2;
double sup3;

shat1=xp/r;shat2=yp/r;shat3=zp/r; 

detginv
=
1/(2.*g12*g13*g23 - g33*pow2(g12) + g22*(g11*g33 - pow2(g13)) - 
    g11*pow2(g23))
;

ginv11
=
detginv*(g22*g33 - pow2(g23))
;

ginv12
=
detginv*(g13*g23 - g12*g33)
;

ginv13
=
detginv*(-(g13*g22) + g12*g23)
;

ginv22
=
detginv*(g11*g33 - pow2(g13))
;

ginv23
=
detginv*(g12*g13 - g11*g23)
;

ginv33
=
detginv*(g11*g22 - pow2(g12))
;

ADMginv11
=
chi*ginv11
;

ADMginv12
=
chi*ginv12
;

ADMginv13
=
chi*ginv13
;

ADMginv22
=
chi*ginv22
;

ADMginv23
=
chi*ginv23
;

ADMginv33
=
chi*ginv33
;

modshatARG
=
2.*(ADMginv23*shat2*shat3 + shat1*(ADMginv12*shat2 + ADMginv13*shat3)) + 
  ADMginv11*pow2(shat1) + ADMginv22*pow2(shat2) + ADMginv33*pow2(shat3)
;


if (modshatARG<0.00001) {                           
      printf("modshat is wrong (%e)\n",modshatARG);
      modshatARG = 0.00001;
    }oomodshat
=
1/sqrt(modshatARG)
;

sdown1
=
oomodshat*shat1
;

sdown2
=
oomodshat*shat2
;

sdown3
=
oomodshat*shat3
;

sup1
=
ADMginv11*sdown1 + ADMginv12*sdown2 + ADMginv13*sdown3
;

sup2
=
ADMginv12*sdown1 + ADMginv22*sdown2 + ADMginv23*sdown3
;

sup3
=
ADMginv13*sdown1 + ADMginv23*sdown2 + ADMginv33*sdown3
;

qud11
=
1. - sdown1*sup1
;

qud12
=
-(sdown2*sup1)
;

qud13
=
-(sdown3*sup1)
;

qud21
=
-(sdown1*sup2)
;

qud22
=
1. - sdown2*sup2
;

qud23
=
-(sdown3*sup2)
;

qud31
=
-(sdown1*sup3)
;

qud32
=
-(sdown2*sup3)
;

qud33
=
1. - sdown3*sup3
;

f1
=
qud11*tof1 + qud12*tof2 + qud13*tof3 + sup1*tofs
;

f2
=
qud21*tof1 + qud22*tof2 + qud23*tof3 + sup2*tofs
;

f3
=
qud31*tof1 + qud32*tof2 + qud33*tof3 + sup3*tofs
;
}  /* function */
// f and tof are lower index
#ifdef fortran1
void compose2p1_2
#endif	
#ifdef fortran2
void COMPOSE2P1_2
#endif
#ifdef fortran3
void compose2p1_2_
#endif
(double & r,double & xp,double & yp,double & zp,double & chi,
		 double & g11,double & g12,double & g13,double & g22,double & g23,double & g33,
		 double & f11,double & f12,double & f13,double & f22,double & f23,double & f33,
		 double & tofqq,double & tofss,double & tofs1,double & tofs2,double & tofs3,
		 double & tof11,double & tof12,double & tof13,double & tof22,double & tof23,double & tof33)
{
double ADMginv11;
double ADMginv12;
double ADMginv13;
double ADMginv22;
double ADMginv23;
double ADMginv33;
double detginv;
double ginv11;
double ginv12;
double ginv13;
double ginv22;
double ginv23;
double ginv33;
double modshatARG;
double oomodshat;
double qdd11;
double qdd12;
double qdd13;
double qdd22;
double qdd23;
double qdd33;
double qPhysuudd1111;
double qPhysuudd1112;
double qPhysuudd1113;
double qPhysuudd1122;
double qPhysuudd1123;
double qPhysuudd1133;
double qPhysuudd1211;
double qPhysuudd1212;
double qPhysuudd1213;
double qPhysuudd1222;
double qPhysuudd1223;
double qPhysuudd1233;
double qPhysuudd1311;
double qPhysuudd1312;
double qPhysuudd1313;
double qPhysuudd1322;
double qPhysuudd1323;
double qPhysuudd1333;
double qPhysuudd2211;
double qPhysuudd2212;
double qPhysuudd2213;
double qPhysuudd2222;
double qPhysuudd2223;
double qPhysuudd2233;
double qPhysuudd2311;
double qPhysuudd2312;
double qPhysuudd2313;
double qPhysuudd2322;
double qPhysuudd2323;
double qPhysuudd2333;
double qPhysuudd3311;
double qPhysuudd3312;
double qPhysuudd3313;
double qPhysuudd3322;
double qPhysuudd3323;
double qPhysuudd3333;
double qud11;
double qud12;
double qud13;
double qud21;
double qud22;
double qud23;
double qud31;
double qud32;
double qud33;
double quu11;
double quu12;
double quu13;
double quu22;
double quu23;
double quu33;
double sdown1;
double sdown2;
double sdown3;
double shat1;
double shat2;
double shat3;
double sup1;
double sup2;
double sup3;



shat1
=
0
;

shat2
=
0
;

shat3
=
0
;


shat1=xp/r;shat2=yp/r;shat3=zp/r; 

detginv
=
1/(2.*g12*g13*g23 - g33*pow2(g12) + g22*(g11*g33 - pow2(g13)) - 
    g11*pow2(g23))
;

ginv11
=
detginv*(g22*g33 - pow2(g23))
;

ginv12
=
detginv*(g13*g23 - g12*g33)
;

ginv13
=
detginv*(-(g13*g22) + g12*g23)
;

ginv22
=
detginv*(g11*g33 - pow2(g13))
;

ginv23
=
detginv*(g12*g13 - g11*g23)
;

ginv33
=
detginv*(g11*g22 - pow2(g12))
;

ADMginv11
=
chi*ginv11
;

ADMginv12
=
chi*ginv12
;

ADMginv13
=
chi*ginv13
;

ADMginv22
=
chi*ginv22
;

ADMginv23
=
chi*ginv23
;

ADMginv33
=
chi*ginv33
;

modshatARG
=
2.*(ADMginv23*shat2*shat3 + shat1*(ADMginv12*shat2 + ADMginv13*shat3)) + 
  ADMginv11*pow2(shat1) + ADMginv22*pow2(shat2) + ADMginv33*pow2(shat3)
;


if (modshatARG<0.00001) {                           
      printf("modshat is wrong (%e)\n",modshatARG);
      modshatARG = 0.00001;
    }oomodshat
=
1/sqrt(modshatARG)
;

sdown1
=
oomodshat*shat1
;

sdown2
=
oomodshat*shat2
;

sdown3
=
oomodshat*shat3
;

sup1
=
ADMginv11*sdown1 + ADMginv12*sdown2 + ADMginv13*sdown3
;

sup2
=
ADMginv12*sdown1 + ADMginv22*sdown2 + ADMginv23*sdown3
;

sup3
=
ADMginv13*sdown1 + ADMginv23*sdown2 + ADMginv33*sdown3
;

qud11
=
1. - sdown1*sup1
;

qud12
=
-(sdown2*sup1)
;

qud13
=
-(sdown3*sup1)
;

qud21
=
-(sdown1*sup2)
;

qud22
=
1. - sdown2*sup2
;

qud23
=
-(sdown3*sup2)
;

qud31
=
-(sdown1*sup3)
;

qud32
=
-(sdown2*sup3)
;

qud33
=
1. - sdown3*sup3
;

qdd11
=
g11/chi - pow2(sdown1)
;

qdd12
=
g12/chi - sdown1*sdown2
;

qdd13
=
g13/chi - sdown1*sdown3
;

qdd22
=
g22/chi - pow2(sdown2)
;

qdd23
=
g23/chi - sdown2*sdown3
;

qdd33
=
g33/chi - pow2(sdown3)
;

quu11
=
ADMginv11 - pow2(sup1)
;

quu12
=
ADMginv12 - sup1*sup2
;

quu13
=
ADMginv13 - sup1*sup3
;

quu22
=
ADMginv22 - pow2(sup2)
;

quu23
=
ADMginv23 - sup2*sup3
;

quu33
=
ADMginv33 - pow2(sup3)
;

qPhysuudd1111
=
-0.5*qdd11*quu11 + pow2(qud11)
;

qPhysuudd1112
=
qud11*qud12 - 0.5*qdd12*quu11
;

qPhysuudd1113
=
qud11*qud13 - 0.5*qdd13*quu11
;

qPhysuudd1122
=
-0.5*qdd22*quu11 + pow2(qud12)
;

qPhysuudd1123
=
qud12*qud13 - 0.5*qdd23*quu11
;

qPhysuudd1133
=
-0.5*qdd33*quu11 + pow2(qud13)
;

qPhysuudd1211
=
qud11*qud21 - 0.5*qdd11*quu12
;

qPhysuudd1212
=
0.5*(qud12*qud21 + qud11*qud22 - qdd12*quu12)
;

qPhysuudd1213
=
0.5*(qud13*qud21 + qud11*qud23 - qdd13*quu12)
;

qPhysuudd1222
=
qud12*qud22 - 0.5*qdd22*quu12
;

qPhysuudd1223
=
0.5*(qud13*qud22 + qud12*qud23 - qdd23*quu12)
;

qPhysuudd1233
=
qud13*qud23 - 0.5*qdd33*quu12
;

qPhysuudd1311
=
qud11*qud31 - 0.5*qdd11*quu13
;

qPhysuudd1312
=
0.5*(qud12*qud31 + qud11*qud32 - qdd12*quu13)
;

qPhysuudd1313
=
0.5*(qud13*qud31 + qud11*qud33 - qdd13*quu13)
;

qPhysuudd1322
=
qud12*qud32 - 0.5*qdd22*quu13
;

qPhysuudd1323
=
0.5*(qud13*qud32 + qud12*qud33 - qdd23*quu13)
;

qPhysuudd1333
=
qud13*qud33 - 0.5*qdd33*quu13
;

qPhysuudd2211
=
-0.5*qdd11*quu22 + pow2(qud21)
;

qPhysuudd2212
=
qud21*qud22 - 0.5*qdd12*quu22
;

qPhysuudd2213
=
qud21*qud23 - 0.5*qdd13*quu22
;

qPhysuudd2222
=
-0.5*qdd22*quu22 + pow2(qud22)
;

qPhysuudd2223
=
qud22*qud23 - 0.5*qdd23*quu22
;

qPhysuudd2233
=
-0.5*qdd33*quu22 + pow2(qud23)
;

qPhysuudd2311
=
qud21*qud31 - 0.5*qdd11*quu23
;

qPhysuudd2312
=
0.5*(qud22*qud31 + qud21*qud32 - qdd12*quu23)
;

qPhysuudd2313
=
0.5*(qud23*qud31 + qud21*qud33 - qdd13*quu23)
;

qPhysuudd2322
=
qud22*qud32 - 0.5*qdd22*quu23
;

qPhysuudd2323
=
0.5*(qud23*qud32 + qud22*qud33 - qdd23*quu23)
;

qPhysuudd2333
=
qud23*qud33 - 0.5*qdd33*quu23
;

qPhysuudd3311
=
-0.5*qdd11*quu33 + pow2(qud31)
;

qPhysuudd3312
=
qud31*qud32 - 0.5*qdd12*quu33
;

qPhysuudd3313
=
qud31*qud33 - 0.5*qdd13*quu33
;

qPhysuudd3322
=
-0.5*qdd22*quu33 + pow2(qud32)
;

qPhysuudd3323
=
qud32*qud33 - 0.5*qdd23*quu33
;

qPhysuudd3333
=
-0.5*qdd33*quu33 + pow2(qud33)
;

// my equations
#if 0
f11
=
qPhysuudd1111*tof11 + qPhysuudd2211*tof22 + 
  2.*(qPhysuudd1211*tof12 + qPhysuudd1311*tof13 + qPhysuudd2311*tof23) + 
  qPhysuudd3311*tof33 + qdd11*tofqq + 
  1.*sdown1*(qud11*tofs1 + qud21*tofs2 + qud31*tofs3) + tofss*pow2(sdown1)
;

f12
=
qPhysuudd1112*tof11 + qPhysuudd2212*tof22 + 
  2.*(qPhysuudd1212*tof12 + qPhysuudd1312*tof13 + qPhysuudd2312*tof23) + 
  qPhysuudd3312*tof33 + qdd12*tofqq + 
  0.5*((qud12*sdown1 + qud11*sdown2)*tofs1 + 
     (qud22*sdown1 + qud21*sdown2)*tofs2 + 
     (qud32*sdown1 + qud31*sdown2)*tofs3) + sdown1*sdown2*tofss
;

f13
=
qPhysuudd1113*tof11 + qPhysuudd2213*tof22 + 
  2.*(qPhysuudd1213*tof12 + qPhysuudd1313*tof13 + qPhysuudd2313*tof23) + 
  qPhysuudd3313*tof33 + qdd13*tofqq + 
  0.5*((qud13*sdown1 + qud11*sdown3)*tofs1 + 
     (qud23*sdown1 + qud21*sdown3)*tofs2 + 
     (qud33*sdown1 + qud31*sdown3)*tofs3) + sdown1*sdown3*tofss
;

f22
=
qPhysuudd1122*tof11 + qPhysuudd2222*tof22 + 
  2.*(qPhysuudd1222*tof12 + qPhysuudd1322*tof13 + qPhysuudd2322*tof23) + 
  qPhysuudd3322*tof33 + qdd22*tofqq + 
  1.*sdown2*(qud12*tofs1 + qud22*tofs2 + qud32*tofs3) + tofss*pow2(sdown2)
;

f23
=
qPhysuudd1123*tof11 + qPhysuudd2223*tof22 + 
  2.*(qPhysuudd1223*tof12 + qPhysuudd1323*tof13 + qPhysuudd2323*tof23) + 
  qPhysuudd3323*tof33 + qdd23*tofqq + 
  0.5*((qud13*sdown2 + qud12*sdown3)*tofs1 + 
     (qud23*sdown2 + qud22*sdown3)*tofs2 + 
     (qud33*sdown2 + qud32*sdown3)*tofs3) + sdown2*sdown3*tofss
;

f33
=
qPhysuudd1133*tof11 + qPhysuudd2233*tof22 + 
  2.*(qPhysuudd1233*tof12 + qPhysuudd1333*tof13 + qPhysuudd2333*tof23) + 
  qPhysuudd3333*tof33 + qdd33*tofqq + 
  1.*sdown3*(qud13*tofs1 + qud23*tofs2 + qud33*tofs3) + tofss*pow2(sdown3)
;
// David's equations
#else
f11
=
tof11 + 0.5*qdd11*tofqq + 2.*sdown1*
   (qud11*tofs1 + qud21*tofs2 + qud31*tofs3) + tofss*pow2(sdown1)
;

f12
=
tof12 + 0.5*qdd12*tofqq + (qud12*sdown1 + qud11*sdown2)*tofs1 + 
  (qud22*sdown1 + qud21*sdown2)*tofs2 + 
  (qud32*sdown1 + qud31*sdown2)*tofs3 + sdown1*sdown2*tofss
;

f13
=
tof13 + 0.5*qdd13*tofqq + (qud13*sdown1 + qud11*sdown3)*tofs1 + 
  (qud23*sdown1 + qud21*sdown3)*tofs2 + 
  (qud33*sdown1 + qud31*sdown3)*tofs3 + sdown1*sdown3*tofss
;

f22
=
tof22 + 0.5*qdd22*tofqq + 2.*sdown2*
   (qud12*tofs1 + qud22*tofs2 + qud32*tofs3) + tofss*pow2(sdown2)
;

f23
=
tof23 + 0.5*qdd23*tofqq + (qud13*sdown2 + qud12*sdown3)*tofs1 + 
  (qud23*sdown2 + qud22*sdown3)*tofs2 + 
  (qud33*sdown2 + qud32*sdown3)*tofs3 + sdown2*sdown3*tofss
;

f33
=
tof33 + 0.5*qdd33*tofqq + 2.*sdown3*
   (qud13*tofs1 + qud23*tofs2 + qud33*tofs3) + tofss*pow2(sdown3)
;
#endif

}   /* function */
#ifdef fortran1
void racqq_point
#endif	
#ifdef fortran2
void RACQQ_POINT
#endif
#ifdef fortran3
void racqq_point_
#endif
(double &A11,
double &A12,
double &A13,
double &A22,
double &A23,
double &A33,
double &alpha,
double &beta1,
double &beta2,
double &beta3,
double &chi,
double &db11,
double &db12,
double &db13,
double &db21,
double &db22,
double &db23,
double &db31,
double &db32,
double &db33,
double &dg111,
double &dg112,
double &dg113,
double &dg122,
double &dg123,
double &dg133,
double &dg211,
double &dg212,
double &dg213,
double &dg222,
double &dg223,
double &dg233,
double &dg311,
double &dg312,
double &dg313,
double &dg322,
double &dg323,
double &dg333,
double &g11,
double &g12,
double &g13,
double &g22,
double &g23,
double &g33,
double &rACqq,
double &rACss)
{

double Ainv11;
double Ainv12;
double Ainv13;
double Ainv22;
double Ainv23;
double Ainv33;
double detginv;
double divbeta;
double ginv11;
double ginv12;
double ginv13;
double ginv22;
double ginv23;
double ginv33;
double lieg11;
double lieg12;
double lieg13;
double lieg22;
double lieg23;
double lieg33;
double totdivbeta;



detginv
=
1/(2.*g12*g13*g23 - g33*pow2(g12) + g22*(g11*g33 - pow2(g13)) - 
    g11*pow2(g23))
;

ginv11
=
detginv*(g22*g33 - pow2(g23))
;

ginv12
=
detginv*(g13*g23 - g12*g33)
;

ginv13
=
detginv*(-(g13*g22) + g12*g23)
;

ginv22
=
detginv*(g11*g33 - pow2(g13))
;

ginv23
=
detginv*(g12*g13 - g11*g23)
;

ginv33
=
detginv*(g11*g22 - pow2(g12))
;

divbeta
=
db11 + db22 + db33
;

totdivbeta
=
0.66666666666666666667*divbeta
;

Ainv11
=
2.*(A23*ginv12*ginv13 + ginv11*(A12*ginv12 + A13*ginv13)) + 
  A11*pow2(ginv11) + A22*pow2(ginv12) + A33*pow2(ginv13)
;

Ainv12
=
ginv11*(A11*ginv12 + A12*ginv22 + A13*ginv23) + 
  ginv12*(A13*ginv13 + A22*ginv22 + A23*ginv23) + 
  ginv13*(A23*ginv22 + A33*ginv23) + A12*pow2(ginv12)
;

Ainv13
=
ginv11*(A11*ginv13 + A12*ginv23 + A13*ginv33) + 
  ginv12*(A12*ginv13 + A22*ginv23 + A23*ginv33) + 
  ginv13*(A23*ginv23 + A33*ginv33) + A13*pow2(ginv13)
;

Ainv22
=
2.*(A23*ginv22*ginv23 + ginv12*(A12*ginv22 + A13*ginv23)) + 
  A11*pow2(ginv12) + A22*pow2(ginv22) + A33*pow2(ginv23)
;

Ainv23
=
ginv13*(A12*ginv22 + A13*ginv23) + A33*ginv23*ginv33 + 
  ginv12*(A11*ginv13 + A12*ginv23 + A13*ginv33) + 
  ginv22*(A22*ginv23 + A23*ginv33) + A23*pow2(ginv23)
;

Ainv33
=
2.*(A23*ginv23*ginv33 + ginv13*(A12*ginv23 + A13*ginv33)) + 
  A11*pow2(ginv13) + A22*pow2(ginv23) + A33*pow2(ginv33)
;

lieg11
=
beta1*dg111 + beta2*dg211 + beta3*dg311 + 
  2.*(db11*g11 + db12*g12 + db13*g13) - g11*totdivbeta
;

lieg12
=
beta1*dg112 + beta2*dg212 + beta3*dg312 + db21*g11 + db23*g13 + db12*g22 + 
  db13*g23 + g12*(db11 + db22 - totdivbeta)
;

lieg13
=
beta1*dg113 + beta2*dg213 + beta3*dg313 + db31*g11 + db32*g12 + db12*g23 + 
  db13*g33 + g13*(db11 + db33 - totdivbeta)
;

lieg22
=
beta1*dg122 + beta2*dg222 + beta3*dg322 + 
  2.*(db21*g12 + db22*g22 + db23*g23) - g22*totdivbeta
;

lieg23
=
beta1*dg123 + beta2*dg223 + beta3*dg323 + db31*g12 + db21*g13 + db32*g22 + 
  db23*g33 + g23*(db22 + db33 - totdivbeta)
;

lieg33
=
beta1*dg133 + beta2*dg233 + beta3*dg333 + 
  2.*(db31*g13 + db32*g23 + db33*g33) - g33*totdivbeta
;

rACqq
=
chi*(-((4.*(A12*Ainv12 + A13*Ainv13 + A23*Ainv23) + 
          2.*(A11*Ainv11 + A22*Ainv22 + A33*Ainv33))*alpha) + 
     Ainv11*lieg11 + Ainv22*lieg22 + 
     2.*(Ainv12*lieg12 + Ainv13*lieg13 + Ainv23*lieg23) + Ainv33*lieg33) - 
  rACss
;

}  /* function */
#ifdef fortran1
void rkhat_point
#endif	
#ifdef fortran2
void RKHAT_POINT
#endif
#ifdef fortran3
void rkhat_point_
#endif
(double &alpha,
double &beta1,
double &beta2,
double &beta3,
double &chi,
double &dKhat1,
double &dKhat2,
double &dKhat3,
double &dTheta1,
double &dTheta2,
double &dTheta3,
double &g11,
double &g12,
double &g13,
double &g22,
double &g23,
double &g33,
double &kappa1,
double &kappa2,
double &Khat,
double &r,
double &rKhat,
double &Theta,
double &xp,
double &yp,
double &zp)
{

double ADMginv11;
double ADMginv12;
double ADMginv13;
double ADMginv22;
double ADMginv23;
double ADMginv33;
double detginv;
double DKhat;
double DTheta;
double ginv11;
double ginv12;
double ginv13;
double ginv22;
double ginv23;
double ginv33;
double lienK;
double lienKhat;
double lienTheta;
double modshatARG;
double muL;
double oomodshat;
double sdown1;
double sdown2;
double sdown3;
double shat1;
double shat2;
double shat3;
double sup1;
double sup2;
double sup3;

shat1=xp/r;shat2=yp/r;shat3=zp/r; 

detginv
=
1/(2.*g12*g13*g23 - g33*pow2(g12) + g22*(g11*g33 - pow2(g13)) - 
    g11*pow2(g23))
;

ginv11
=
detginv*(g22*g33 - pow2(g23))
;

ginv12
=
detginv*(g13*g23 - g12*g33)
;

ginv13
=
detginv*(-(g13*g22) + g12*g23)
;

ginv22
=
detginv*(g11*g33 - pow2(g13))
;

ginv23
=
detginv*(g12*g13 - g11*g23)
;

ginv33
=
detginv*(g11*g22 - pow2(g12))
;

ADMginv11
=
chi*ginv11
;

ADMginv12
=
chi*ginv12
;

ADMginv13
=
chi*ginv13
;

ADMginv22
=
chi*ginv22
;

ADMginv23
=
chi*ginv23
;

ADMginv33
=
chi*ginv33
;

modshatARG
=
2.*(ADMginv23*shat2*shat3 + shat1*(ADMginv12*shat2 + ADMginv13*shat3)) + 
  ADMginv11*pow2(shat1) + ADMginv22*pow2(shat2) + ADMginv33*pow2(shat3)
;


if (modshatARG<0.00001) {                           
      printf("modshat is wrong (%e)\n",modshatARG);
      modshatARG = 0.00001;
    }oomodshat
=
1/sqrt(modshatARG)
;

sdown1
=
oomodshat*shat1
;

sdown2
=
oomodshat*shat2
;

sdown3
=
oomodshat*shat3
;

sup1
=
ADMginv11*sdown1 + ADMginv12*sdown2 + ADMginv13*sdown3
;

sup2
=
ADMginv12*sdown1 + ADMginv22*sdown2 + ADMginv23*sdown3
;

sup3
=
ADMginv13*sdown1 + ADMginv23*sdown2 + ADMginv33*sdown3
;

muL
=
2./alpha
;

DKhat
=
dKhat1*sup1 + dKhat2*sup2 + dKhat3*sup3
;

DTheta
=
dTheta1*sup1 + dTheta2*sup2 + dTheta3*sup3
;

lienKhat
=
-((DKhat + Khat/r)*sqrt(muL))
;

lienTheta
=
-DTheta - (kappa1*(2. + kappa2) + 1/r)*Theta
;

lienK
=
lienKhat + 2.*lienTheta
;

rKhat
=
beta1*dKhat1 + beta2*dKhat2 + beta3*dKhat3 + alpha*lienKhat
;

}  /* function */
#ifdef fortran1
void rtheta_point
#endif	
#ifdef fortran2
void RTHETA_POINT
#endif
#ifdef fortran3
void rtheta_point_
#endif
(double &alpha,
double &beta1,
double &beta2,
double &beta3,
double &chi,
double &dTheta1,
double &dTheta2,
double &dTheta3,
double &g11,
double &g12,
double &g13,
double &g22,
double &g23,
double &g33,
double &kappa1,
double &kappa2,
double &r,
double &rTheta,
double &Theta,
double &xp,
double &yp,
double &zp)
{

double ADMginv11;
double ADMginv12;
double ADMginv13;
double ADMginv22;
double ADMginv23;
double ADMginv33;
double detginv;
double DTheta;
double ginv11;
double ginv12;
double ginv13;
double ginv22;
double ginv23;
double ginv33;
double lienTheta;
double modshatARG;
double oomodshat;
double sdown1;
double sdown2;
double sdown3;
double shat1;
double shat2;
double shat3;
double sup1;
double sup2;
double sup3;



shat1
=
0
;

shat2
=
0
;

shat3
=
0
;


shat1=xp/r;shat2=yp/r;shat3=zp/r; 

detginv
=
1/(2.*g12*g13*g23 - g33*pow2(g12) + g22*(g11*g33 - pow2(g13)) - 
    g11*pow2(g23))
;

ginv11
=
detginv*(g22*g33 - pow2(g23))
;

ginv12
=
detginv*(g13*g23 - g12*g33)
;

ginv13
=
detginv*(-(g13*g22) + g12*g23)
;

ginv22
=
detginv*(g11*g33 - pow2(g13))
;

ginv23
=
detginv*(g12*g13 - g11*g23)
;

ginv33
=
detginv*(g11*g22 - pow2(g12))
;

ADMginv11
=
chi*ginv11
;

ADMginv12
=
chi*ginv12
;

ADMginv13
=
chi*ginv13
;

ADMginv22
=
chi*ginv22
;

ADMginv23
=
chi*ginv23
;

ADMginv33
=
chi*ginv33
;

modshatARG
=
2.*(ADMginv23*shat2*shat3 + shat1*(ADMginv12*shat2 + ADMginv13*shat3)) + 
  ADMginv11*pow2(shat1) + ADMginv22*pow2(shat2) + ADMginv33*pow2(shat3)
;


if (modshatARG<0.00001) {                           
      printf("modshat is wrong (%e)\n",modshatARG);
      modshatARG = 0.00001;
    }oomodshat
=
1/sqrt(modshatARG)
;

sdown1
=
oomodshat*shat1
;

sdown2
=
oomodshat*shat2
;

sdown3
=
oomodshat*shat3
;

sup1
=
ADMginv11*sdown1 + ADMginv12*sdown2 + ADMginv13*sdown3
;

sup2
=
ADMginv12*sdown1 + ADMginv22*sdown2 + ADMginv23*sdown3
;

sup3
=
ADMginv13*sdown1 + ADMginv23*sdown2 + ADMginv33*sdown3
;

DTheta
=
dTheta1*sup1 + dTheta2*sup2 + dTheta3*sup3
;

lienTheta
=
-DTheta - (kappa1*(2. + kappa2) + 1/r)*Theta
;

rTheta
=
beta1*dTheta1 + beta2*dTheta2 + beta3*dTheta3 + alpha*lienTheta
;

}  /* function */

#ifdef fortran1
void rgam_point
#endif	
#ifdef fortran2
void RGAM_POINT
#endif
#ifdef fortran3
void rgam_point_
#endif
(double &alpha,
double &beta1,
double &beta2,
double &beta3,
double &chi,
double &db11,
double &db12,
double &db13,
double &db21,
double &db22,
double &db23,
double &db31,
double &db32,
double &db33,
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
double &dG11,
double &dG12,
double &dG13,
double &dG21,
double &dG22,
double &dG23,
double &dG31,
double &dG32,
double &dG33,
double &dKhat1,
double &dKhat2,
double &dKhat3,
double &dTheta1,
double &dTheta2,
double &dTheta3,
double &g11,
double &g12,
double &g13,
double &g22,
double &g23,
double &g33,
double &r,
double &rGamA1,
double &rGamA2,
double &rGamA3,
double &rGams,
double &shiftdriver,
double &xp,
double &yp,
double &zp)
{

double ADMginv11;
double ADMginv12;
double ADMginv13;
double ADMginv22;
double ADMginv23;
double ADMginv33;
double detginv;
double ginv11;
double ginv12;
double ginv13;
double ginv22;
double ginv23;
double ginv33;
double modshatARG;
double muL;
double muStilde;
double oomodshat;
double qdd11;
double qdd12;
double qdd13;
double qdd22;
double qdd23;
double qdd33;
double qud11;
double qud12;
double qud13;
double qud21;
double qud22;
double qud23;
double qud31;
double qud32;
double qud33;
double quu11;
double quu12;
double quu13;
double quu22;
double quu23;
double quu33;
double sdown1;
double sdown2;
double sdown3;
double shat1;
double shat2;
double shat3;
double sup1;
double sup2;
double sup3;
double vbetaA;
double vbetas;

shat1=xp/r;shat2=yp/r;shat3=zp/r; 

detginv
=
1/(2.*g12*g13*g23 - g33*pow2(g12) + g22*(g11*g33 - pow2(g13)) - 
    g11*pow2(g23))
;

ginv11
=
detginv*(g22*g33 - pow2(g23))
;

ginv12
=
detginv*(g13*g23 - g12*g33)
;

ginv13
=
detginv*(-(g13*g22) + g12*g23)
;

ginv22
=
detginv*(g11*g33 - pow2(g13))
;

ginv23
=
detginv*(g12*g13 - g11*g23)
;

ginv33
=
detginv*(g11*g22 - pow2(g12))
;

ADMginv11
=
chi*ginv11
;

ADMginv12
=
chi*ginv12
;

ADMginv13
=
chi*ginv13
;

ADMginv22
=
chi*ginv22
;

ADMginv23
=
chi*ginv23
;

ADMginv33
=
chi*ginv33
;

modshatARG
=
2.*(ADMginv23*shat2*shat3 + shat1*(ADMginv12*shat2 + ADMginv13*shat3)) + 
  ADMginv11*pow2(shat1) + ADMginv22*pow2(shat2) + ADMginv33*pow2(shat3)
;


if (modshatARG<0.00001) {                           
      printf("modshat is wrong (%e)\n",modshatARG);
      modshatARG = 0.00001;
    }oomodshat
=
1/sqrt(modshatARG)
;

sdown1
=
oomodshat*shat1
;

sdown2
=
oomodshat*shat2
;

sdown3
=
oomodshat*shat3
;

sup1
=
ADMginv11*sdown1 + ADMginv12*sdown2 + ADMginv13*sdown3
;

sup2
=
ADMginv12*sdown1 + ADMginv22*sdown2 + ADMginv23*sdown3
;

sup3
=
ADMginv13*sdown1 + ADMginv23*sdown2 + ADMginv33*sdown3
;

qud11
=
1. - sdown1*sup1
;

qud12
=
-(sdown2*sup1)
;

qud13
=
-(sdown3*sup1)
;

qud21
=
-(sdown1*sup2)
;

qud22
=
1. - sdown2*sup2
;

qud23
=
-(sdown3*sup2)
;

qud31
=
-(sdown1*sup3)
;

qud32
=
-(sdown2*sup3)
;

qud33
=
1. - sdown3*sup3
;

qdd11
=
g11/chi - pow2(sdown1)
;

qdd12
=
g12/chi - sdown1*sdown2
;

qdd13
=
g13/chi - sdown1*sdown3
;

qdd22
=
g22/chi - pow2(sdown2)
;

qdd23
=
g23/chi - sdown2*sdown3
;

qdd33
=
g33/chi - pow2(sdown3)
;

quu11
=
ADMginv11 - pow2(sup1)
;

quu12
=
ADMginv12 - sup1*sup2
;

quu13
=
ADMginv13 - sup1*sup3
;

quu22
=
ADMginv22 - pow2(sup2)
;

quu23
=
ADMginv23 - sup2*sup3
;

quu33
=
ADMginv33 - pow2(sup3)
;

muL
=
2./alpha
;

muStilde
=
1/chi
;

vbetas
=
2.*sqrt(0.33333333333333333333*muStilde)
;

vbetaA
=
sqrt(muStilde)
;

rGams
=
(beta1*dG11 + beta2*dG21 + beta3*dG31 + 
     (ddb111*quu11 + ddb221*quu22 + 
        2.*(ddb121*quu12 + ddb131*quu13 + ddb231*quu23) + ddb331*quu33)/chi\
)*sdown1 + (beta1*dG12 + beta2*dG22 + beta3*dG32 + 
     (ddb112*quu11 + ddb222*quu22 + 
        2.*(ddb122*quu12 + ddb132*quu13 + ddb232*quu23) + ddb332*quu33)/chi\
)*sdown2 + (beta1*dG13 + beta2*dG23 + beta3*dG33 + 
     (ddb113*quu11 + ddb223*quu22 + 
        2.*(ddb123*quu12 + ddb133*quu13 + ddb233*quu23) + ddb333*quu33)/chi\
)*sdown3 - ((ddb111*qud11 + ddb112*qud12 + ddb113*qud13 + ddb121*qud21 + 
        ddb122*qud22 + ddb123*qud23 + ddb131*qud31 + ddb132*qud32 + 
        ddb133*qud33)*sup1 + (ddb121*qud11 + ddb122*qud12 + 
        ddb123*qud13 + ddb221*qud21 + ddb222*qud22 + ddb223*qud23 + 
        ddb231*qud31 + ddb232*qud32 + ddb233*qud33)*sup2 + 
     (ddb131*qud11 + ddb132*qud12 + ddb133*qud13 + ddb231*qud21 + 
        ddb232*qud22 + ddb233*qud23 + ddb331*qud31 + ddb332*qud32 + 
        ddb333*qud33)*sup3)/chi - (dG11 + dG22 + dG33)*vbetas + 
  2.*((0.33333333333333333333*alpha*
        (dTheta1*sup1 + dTheta2*sup2 + dTheta3*sup3))/(chi + chi*vbetas) + 
     ((db11 + db22 + db33)*shiftdriver)/(vbetaA*sqrt(3.))) + 
  (1.3333333333333333333*alpha*(dKhat1*sup1 + dKhat2*sup2 + dKhat3*sup3)*
     sqrt(muL))/(chi*(vbetas + sqrt(muL)))
;

rGamA1
=
-(((dG21*qud11 + dG22*qud12 + dG23*qud13)*sup2 + 
       (dG31*qud11 + dG32*qud12 + dG33*qud13)*sup3)*vbetaA) + 
  qud11*(beta2*dG21 + beta3*dG31 + 
     (1.3333333333333333333*ddb111*quu11 + 
        2.3333333333333333333*(ddb121*quu12 + ddb131*quu13) + 
        ddb221*quu22 + ddb331*quu33 + 
        (shiftdriver*(db11*sup1 + db31*sup3))/vbetaA)/chi + 
     dG11*(beta1 - sup1*vbetaA)) + 
  qud12*(beta2*dG22 + beta3*dG32 + 
     (1.3333333333333333333*ddb112*quu11 + 
        2.3333333333333333333*(ddb122*quu12 + ddb132*quu13) + 
        ddb222*quu22 + 2.*ddb232*quu23 + ddb332*quu33 + 
        (shiftdriver*(db12*sup1 + db22*sup2 + db32*sup3))/vbetaA)/chi + 
     dG12*(beta1 - sup1*vbetaA)) + 
  qud13*(beta2*dG23 + beta3*dG33 + 
     (1.3333333333333333333*ddb113*quu11 + 
        2.3333333333333333333*(ddb123*quu12 + ddb133*quu13) + 
        ddb223*quu22 + 2.*ddb233*quu23 + ddb333*quu33 + 
        (shiftdriver*(db13*sup1 + db23*sup2 + db33*sup3))/vbetaA)/chi + 
     dG13*(beta1 - sup1*vbetaA)) + 
  (0.33333333333333333333*((ddb121*qud21 + ddb122*qud22 + ddb123*qud23 + 
           ddb131*qud31 + ddb132*qud32 + ddb133*qud33)*quu11 + 
        (ddb221*qud21 + ddb223*qud23 + ddb231*qud31 + ddb232*qud32 + 
           ddb233*qud33)*quu12 + 
        (ddb231*qud21 + ddb232*qud22 + ddb233*qud23 + ddb331*qud31 + 
           ddb332*qud32)*quu13) - 
     alpha*((1.3333333333333333333*dKhat1 + 
           0.66666666666666666667*dTheta1)*quu11 + 
        1.3333333333333333333*(dKhat2*quu12 + dKhat3*quu13)) + 
     1.3333333333333333333*((ddb132*quu13*sdown2 + ddb113*quu11*sdown3)*
         sup1 + (quu13*(ddb231*sdown1 + ddb232*sdown2) + 
           quu12*(ddb222*sdown2 + ddb223*sdown3))*sup2 + 
        (quu12*(ddb232*sdown2 + ddb233*sdown3) + 
           quu13*(ddb331*sdown1 + ddb332*sdown2 + ddb333*sdown3))*sup3 + 
        sdown1*((ddb121*quu12 + ddb131*quu13)*sup1 + ddb221*quu12*sup2 + 
           ddb131*quu11*sup3) + 
        sdown2*((ddb112*quu11 + ddb122*quu12)*sup1 + 
           quu11*(ddb122*sup2 + ddb132*sup3)) + 
        sdown3*((ddb123*quu12 + ddb133*quu13)*sup1 + 
           quu11*(ddb123*sup2 + ddb133*sup3))) + 
     qud11*(2.*ddb231*quu23 + (db21*shiftdriver*sup2)/vbetaA) - 
     (((db11*quu11 + db21*quu12)*sdown1 + 
          (db12*quu11 + db22*quu12 + db32*quu13)*sdown2 + 
          (db13*quu11 + db23*quu12 + db33*quu13)*sdown3)*shiftdriver)/
      vbetaA + ((dG22*quu12 + dG32*quu13)*sdown2 + 
        (dG13*quu11 + dG23*quu12)*sdown3)*vbetaA + 
     quu11*(1.3333333333333333333*sdown1*(ddb111*sup1 + ddb121*sup2) + 
        (dG11*sdown1 + dG12*sdown2)*vbetaA) + 
     quu12*(-0.66666666666666666667*alpha*dTheta2 + 
        0.33333333333333333333*ddb222*qud22 + 
        sdown1*(1.3333333333333333333*ddb231*sup3 + dG21*vbetaA)) + 
     quu13*(-0.66666666666666666667*alpha*dTheta3 + 
        0.33333333333333333333*ddb333*qud33 - 
        (db31*sdown1*shiftdriver)/vbetaA + dG31*sdown1*vbetaA + 
        sdown3*(1.3333333333333333333*ddb233*sup2 + dG33*vbetaA)))/chi
;

rGamA2
=
-(((dG21*qud21 + dG22*qud22 + dG23*qud23)*sup2 + 
       (dG31*qud21 + dG32*qud22 + dG33*qud23)*sup3)*vbetaA) + 
  qud21*(beta2*dG21 + beta3*dG31 + 
     (ddb111*quu11 + 2.*ddb131*quu13 + 
        1.3333333333333333333*ddb221*quu22 + 
        2.3333333333333333333*(ddb121*quu12 + ddb231*quu23) + 
        ddb331*quu33 + (shiftdriver*(db11*sup1 + db31*sup3))/vbetaA)/chi + 
     dG11*(beta1 - sup1*vbetaA)) + 
  qud22*(beta2*dG22 + beta3*dG32 + 
     (ddb112*quu11 + 2.*ddb132*quu13 + 
        1.3333333333333333333*ddb222*quu22 + 
        2.3333333333333333333*(ddb122*quu12 + ddb232*quu23) + 
        ddb332*quu33 + (shiftdriver*(db12*sup1 + db22*sup2 + db32*sup3))/
         vbetaA)/chi + dG12*(beta1 - sup1*vbetaA)) + 
  qud23*(beta2*dG23 + beta3*dG33 + 
     (ddb113*quu11 + 2.*ddb133*quu13 + 
        1.3333333333333333333*ddb223*quu22 + 
        2.3333333333333333333*(ddb123*quu12 + ddb233*quu23) + 
        ddb333*quu33 + (shiftdriver*(db13*sup1 + db23*sup2 + db33*sup3))/
         vbetaA)/chi + dG13*(beta1 - sup1*vbetaA)) + 
  (0.33333333333333333333*((ddb111*qud11 + ddb112*qud12 + ddb113*qud13 + 
           ddb131*qud31 + ddb132*qud32 + ddb133*qud33)*quu12 + 
        (ddb121*qud11 + ddb123*qud13 + ddb231*qud31 + ddb232*qud32 + 
           ddb233*qud33)*quu22 + 
        (ddb131*qud11 + ddb132*qud12 + ddb133*qud13 + ddb331*qud31 + 
           ddb332*qud32)*quu23) - 
     alpha*((1.3333333333333333333*dKhat1 + 
           0.66666666666666666667*dTheta1)*quu12 + 
        1.3333333333333333333*(dKhat2*quu22 + dKhat3*quu23)) + 
     1.3333333333333333333*((ddb132*quu23*sdown2 + ddb113*quu12*sdown3)*
         sup1 + (quu23*(ddb231*sdown1 + ddb232*sdown2) + 
           quu22*(ddb222*sdown2 + ddb223*sdown3))*sup2 + 
        (quu22*(ddb232*sdown2 + ddb233*sdown3) + 
           quu23*(ddb331*sdown1 + ddb332*sdown2 + ddb333*sdown3))*sup3 + 
        sdown1*((ddb121*quu22 + ddb131*quu23)*sup1 + ddb221*quu22*sup2 + 
           ddb131*quu12*sup3) + 
        sdown2*((ddb112*quu12 + ddb122*quu22)*sup1 + 
           quu12*(ddb122*sup2 + ddb132*sup3)) + 
        sdown3*((ddb123*quu22 + ddb133*quu23)*sup1 + 
           quu12*(ddb123*sup2 + ddb133*sup3))) - 
     (((db11*quu12 + db21*quu22)*sdown1 + 
          (db12*quu12 + db22*quu22 + db32*quu23)*sdown2 + 
          (db13*quu12 + db23*quu22 + db33*quu23)*sdown3)*shiftdriver)/
      vbetaA + (db21*qud21*shiftdriver*sup2)/vbetaA + 
     ((dG22*quu22 + dG32*quu23)*sdown2 + (dG13*quu12 + dG23*quu22)*sdown3)*
      vbetaA + quu12*(1.3333333333333333333*sdown1*
         (ddb111*sup1 + ddb121*sup2) + (dG11*sdown1 + dG12*sdown2)*vbetaA) \
+ quu22*(-0.66666666666666666667*alpha*dTheta2 + 
        0.33333333333333333333*ddb122*qud12 + 
        sdown1*(1.3333333333333333333*ddb231*sup3 + dG21*vbetaA)) + 
     quu23*(-0.66666666666666666667*alpha*dTheta3 + 
        0.33333333333333333333*ddb333*qud33 - 
        (db31*sdown1*shiftdriver)/vbetaA + dG31*sdown1*vbetaA + 
        sdown3*(1.3333333333333333333*ddb233*sup2 + dG33*vbetaA)))/chi
;

rGamA3
=
-(((dG21*qud31 + dG22*qud32 + dG23*qud33)*sup2 + 
       (dG31*qud31 + dG32*qud32 + dG33*qud33)*sup3)*vbetaA) + 
  qud31*(beta2*dG21 + beta3*dG31 + 
     (ddb111*quu11 + 2.*ddb121*quu12 + ddb221*quu22 + 
        2.3333333333333333333*(ddb131*quu13 + ddb231*quu23) + 
        1.3333333333333333333*ddb331*quu33 + 
        (shiftdriver*(db11*sup1 + db31*sup3))/vbetaA)/chi + 
     dG11*(beta1 - sup1*vbetaA)) + 
  qud32*(beta2*dG22 + beta3*dG32 + 
     (ddb112*quu11 + 2.*ddb122*quu12 + ddb222*quu22 + 
        2.3333333333333333333*(ddb132*quu13 + ddb232*quu23) + 
        1.3333333333333333333*ddb332*quu33 + 
        (shiftdriver*(db12*sup1 + db22*sup2 + db32*sup3))/vbetaA)/chi + 
     dG12*(beta1 - sup1*vbetaA)) + 
  qud33*(beta2*dG23 + beta3*dG33 + 
     (ddb113*quu11 + 2.*ddb123*quu12 + ddb223*quu22 + 
        2.3333333333333333333*(ddb133*quu13 + ddb233*quu23) + 
        1.3333333333333333333*ddb333*quu33 + 
        (shiftdriver*(db13*sup1 + db23*sup2 + db33*sup3))/vbetaA)/chi + 
     dG13*(beta1 - sup1*vbetaA)) + 
  (0.33333333333333333333*((ddb111*qud11 + ddb112*qud12 + ddb113*qud13 + 
           ddb121*qud21 + ddb122*qud22 + ddb123*qud23)*quu13 + 
        (ddb121*qud11 + ddb123*qud13 + ddb221*qud21 + ddb222*qud22 + 
           ddb223*qud23)*quu23 + 
        (ddb131*qud11 + ddb132*qud12 + ddb133*qud13 + ddb231*qud21 + 
           ddb232*qud22)*quu33) - 
     alpha*((1.3333333333333333333*dKhat1 + 
           0.66666666666666666667*dTheta1)*quu13 + 
        1.3333333333333333333*(dKhat2*quu23 + dKhat3*quu33)) + 
     1.3333333333333333333*((ddb132*quu33*sdown2 + ddb113*quu13*sdown3)*
         sup1 + (quu33*(ddb231*sdown1 + ddb232*sdown2) + 
           quu23*(ddb222*sdown2 + ddb223*sdown3))*sup2 + 
        (quu23*(ddb232*sdown2 + ddb233*sdown3) + 
           quu33*(ddb331*sdown1 + ddb332*sdown2 + ddb333*sdown3))*sup3 + 
        sdown1*((ddb121*quu23 + ddb131*quu33)*sup1 + ddb221*quu23*sup2 + 
           ddb131*quu13*sup3) + 
        sdown2*((ddb112*quu13 + ddb122*quu23)*sup1 + 
           quu13*(ddb122*sup2 + ddb132*sup3)) + 
        sdown3*((ddb123*quu23 + ddb133*quu33)*sup1 + 
           quu13*(ddb123*sup2 + ddb133*sup3))) - 
     (((db11*quu13 + db21*quu23)*sdown1 + 
          (db12*quu13 + db22*quu23 + db32*quu33)*sdown2 + 
          (db13*quu13 + db23*quu23 + db33*quu33)*sdown3)*shiftdriver)/
      vbetaA + (db21*qud31*shiftdriver*sup2)/vbetaA + 
     ((dG22*quu23 + dG32*quu33)*sdown2 + (dG13*quu13 + dG23*quu23)*sdown3)*
      vbetaA + quu13*(1.3333333333333333333*sdown1*
         (ddb111*sup1 + ddb121*sup2) + (dG11*sdown1 + dG12*sdown2)*vbetaA) \
+ quu33*(-0.66666666666666666667*alpha*dTheta3 + 
        ddb233*(0.33333333333333333333*qud23 + 
           1.3333333333333333333*sdown3*sup2) - 
        (db31*sdown1*shiftdriver)/vbetaA + 
        (dG31*sdown1 + dG33*sdown3)*vbetaA) + 
     quu23*(-0.66666666666666666667*alpha*dTheta2 + 
        0.33333333333333333333*ddb122*qud12 + 
        sdown1*(1.3333333333333333333*ddb231*sup3 + dG21*vbetaA)))/chi
;

}  /* function */
#ifdef fortran1
void ra_point
#endif	
#ifdef fortran2
void RA_POINT
#endif
#ifdef fortran3
void ra_point_
#endif
(double &A11,
double &A12,
double &A13,
double &A22,
double &A23,
double &A33,
double &alpha,
double &beta1,
double &beta2,
double &beta3,
double &chi,
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
double &db12,
double &db13,
double &db21,
double &db22,
double &db23,
double &db31,
double &db32,
double &db33,
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
double &ddg1111,
double &ddg1112,
double &ddg1113,
double &ddg1122,
double &ddg1123,
double &ddg1133,
double &ddg1211,
double &ddg1212,
double &ddg1213,
double &ddg1222,
double &ddg1223,
double &ddg1233,
double &ddg1311,
double &ddg1312,
double &ddg1313,
double &ddg1322,
double &ddg1323,
double &ddg1333,
double &ddg2211,
double &ddg2212,
double &ddg2213,
double &ddg2222,
double &ddg2223,
double &ddg2233,
double &ddg2311,
double &ddg2312,
double &ddg2313,
double &ddg2322,
double &ddg2323,
double &ddg2333,
double &ddg3311,
double &ddg3312,
double &ddg3313,
double &ddg3322,
double &ddg3323,
double &ddg3333,
double &dG11,
double &dg111,
double &dg112,
double &dg113,
double &dG12,
double &dg122,
double &dg123,
double &dG13,
double &dg133,
double &dG21,
double &dg211,
double &dg212,
double &dg213,
double &dG22,
double &dg222,
double &dg223,
double &dG23,
double &dg233,
double &dG31,
double &dg311,
double &dg312,
double &dg313,
double &dG32,
double &dg322,
double &dg323,
double &dG33,
double &dg333,
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
double &Khat,
double &r,
double &rACABTF11,
double &rACABTF12,
double &rACABTF13,
double &rACABTF22,
double &rACABTF23,
double &rACABTF33,
double &rACsA1,
double &rACsA2,
double &rACsA3,
double &rACss,
double &Theta,
double &xp,
double &yp,
double &zp)
{

double AA11;
double AA12;
double AA13;
double AA21;
double AA22;
double AA23;
double AA31;
double AA32;
double AA33;
double ADMginv11;
double ADMginv12;
double ADMginv13;
double ADMginv22;
double ADMginv23;
double ADMginv33;
double cdA111;
double cdA112;
double cdA113;
double cdA122;
double cdA123;
double cdA133;
double cdA211;
double cdA212;
double cdA213;
double cdA222;
double cdA223;
double cdA233;
double cdA311;
double cdA312;
double cdA313;
double cdA322;
double cdA323;
double cdA333;
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
double dGfromgdu11;
double dGfromgdu12;
double dGfromgdu13;
double dGfromgdu21;
double dGfromgdu22;
double dGfromgdu23;
double dGfromgdu31;
double dGfromgdu32;
double dGfromgdu33;
double divbeta;
double dK1;
double dK2;
double dK3;
double DTheta;
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
double Gfromg1;
double Gfromg2;
double Gfromg3;
double ginv11;
double ginv12;
double ginv13;
double ginv22;
double ginv23;
double ginv33;
double K;
double lieA11;
double lieA12;
double lieA13;
double lieA22;
double lieA23;
double lieA33;
double modshatARG;
double oochipsipower;
double oomodshat;
double psim4;
double qdd11;
double qdd12;
double qdd13;
double qdd22;
double qdd23;
double qdd33;
double qPhysuudd1111;
double qPhysuudd1112;
double qPhysuudd1113;
double qPhysuudd1122;
double qPhysuudd1123;
double qPhysuudd1133;
double qPhysuudd1211;
double qPhysuudd1212;
double qPhysuudd1213;
double qPhysuudd1222;
double qPhysuudd1223;
double qPhysuudd1233;
double qPhysuudd1311;
double qPhysuudd1312;
double qPhysuudd1313;
double qPhysuudd1322;
double qPhysuudd1323;
double qPhysuudd1333;
double qPhysuudd2211;
double qPhysuudd2212;
double qPhysuudd2213;
double qPhysuudd2222;
double qPhysuudd2223;
double qPhysuudd2233;
double qPhysuudd2311;
double qPhysuudd2312;
double qPhysuudd2313;
double qPhysuudd2322;
double qPhysuudd2323;
double qPhysuudd2333;
double qPhysuudd3311;
double qPhysuudd3312;
double qPhysuudd3313;
double qPhysuudd3322;
double qPhysuudd3323;
double qPhysuudd3333;
double qud11;
double qud12;
double qud13;
double qud21;
double qud22;
double qud23;
double qud31;
double qud32;
double qud33;
double quu11;
double quu12;
double quu13;
double quu22;
double quu23;
double quu33;
double R11;
double R12;
double R13;
double R22;
double R23;
double R33;
double Rf11;
double Rf12;
double Rf13;
double Rf22;
double Rf23;
double Rf33;
double Rhat;
double Rphi11;
double Rphi12;
double Rphi13;
double Rphi22;
double Rphi23;
double Rphi33;
double sdown1;
double sdown2;
double sdown3;
double shat1;
double shat2;
double shat3;
double sup1;
double sup2;
double sup3;
double totdivbeta;
double trcdda;
double trcddf;



chipsipower
=
-4.
;

shat1=xp/r;shat2=yp/r;shat3=zp/r; 

detginv
=
1/(2.*g12*g13*g23 - g33*pow2(g12) + g22*(g11*g33 - pow2(g13)) - 
    g11*pow2(g23))
;

ginv11
=
detginv*(g22*g33 - pow2(g23))
;

ginv12
=
detginv*(g13*g23 - g12*g33)
;

ginv13
=
detginv*(-(g13*g22) + g12*g23)
;

ginv22
=
detginv*(g11*g33 - pow2(g13))
;

ginv23
=
detginv*(g12*g13 - g11*g23)
;

ginv33
=
detginv*(g11*g22 - pow2(g12))
;

ADMginv11
=
chi*ginv11
;

ADMginv12
=
chi*ginv12
;

ADMginv13
=
chi*ginv13
;

ADMginv22
=
chi*ginv22
;

ADMginv23
=
chi*ginv23
;

ADMginv33
=
chi*ginv33
;

modshatARG
=
2.*(ADMginv23*shat2*shat3 + shat1*(ADMginv12*shat2 + ADMginv13*shat3)) + 
  ADMginv11*pow2(shat1) + ADMginv22*pow2(shat2) + ADMginv33*pow2(shat3)
;


if (modshatARG<0.00001) {                           
      printf("modshat is wrong (%e)\n",modshatARG);
      modshatARG = 0.00001;
    }oomodshat
=
1/sqrt(modshatARG)
;

sdown1
=
oomodshat*shat1
;

sdown2
=
oomodshat*shat2
;

sdown3
=
oomodshat*shat3
;

sup1
=
ADMginv11*sdown1 + ADMginv12*sdown2 + ADMginv13*sdown3
;

sup2
=
ADMginv12*sdown1 + ADMginv22*sdown2 + ADMginv23*sdown3
;

sup3
=
ADMginv13*sdown1 + ADMginv23*sdown2 + ADMginv33*sdown3
;

qud11
=
1. - sdown1*sup1
;

qud12
=
-(sdown2*sup1)
;

qud13
=
-(sdown3*sup1)
;

qud21
=
-(sdown1*sup2)
;

qud22
=
1. - sdown2*sup2
;

qud23
=
-(sdown3*sup2)
;

qud31
=
-(sdown1*sup3)
;

qud32
=
-(sdown2*sup3)
;

qud33
=
1. - sdown3*sup3
;

qdd11
=
g11/chi - pow2(sdown1)
;

qdd12
=
g12/chi - sdown1*sdown2
;

qdd13
=
g13/chi - sdown1*sdown3
;

qdd22
=
g22/chi - pow2(sdown2)
;

qdd23
=
g23/chi - sdown2*sdown3
;

qdd33
=
g33/chi - pow2(sdown3)
;

quu11
=
ADMginv11 - pow2(sup1)
;

quu12
=
ADMginv12 - sup1*sup2
;

quu13
=
ADMginv13 - sup1*sup3
;

quu22
=
ADMginv22 - pow2(sup2)
;

quu23
=
ADMginv23 - sup2*sup3
;

quu33
=
ADMginv33 - pow2(sup3)
;

qPhysuudd1111
=
-0.5*qdd11*quu11 + pow2(qud11)
;

qPhysuudd1112
=
qud11*qud12 - 0.5*qdd12*quu11
;

qPhysuudd1113
=
qud11*qud13 - 0.5*qdd13*quu11
;

qPhysuudd1122
=
-0.5*qdd22*quu11 + pow2(qud12)
;

qPhysuudd1123
=
qud12*qud13 - 0.5*qdd23*quu11
;

qPhysuudd1133
=
-0.5*qdd33*quu11 + pow2(qud13)
;

qPhysuudd1211
=
qud11*qud21 - 0.5*qdd11*quu12
;

qPhysuudd1212
=
0.5*(qud12*qud21 + qud11*qud22 - qdd12*quu12)
;

qPhysuudd1213
=
0.5*(qud13*qud21 + qud11*qud23 - qdd13*quu12)
;

qPhysuudd1222
=
qud12*qud22 - 0.5*qdd22*quu12
;

qPhysuudd1223
=
0.5*(qud13*qud22 + qud12*qud23 - qdd23*quu12)
;

qPhysuudd1233
=
qud13*qud23 - 0.5*qdd33*quu12
;

qPhysuudd1311
=
qud11*qud31 - 0.5*qdd11*quu13
;

qPhysuudd1312
=
0.5*(qud12*qud31 + qud11*qud32 - qdd12*quu13)
;

qPhysuudd1313
=
0.5*(qud13*qud31 + qud11*qud33 - qdd13*quu13)
;

qPhysuudd1322
=
qud12*qud32 - 0.5*qdd22*quu13
;

qPhysuudd1323
=
0.5*(qud13*qud32 + qud12*qud33 - qdd23*quu13)
;

qPhysuudd1333
=
qud13*qud33 - 0.5*qdd33*quu13
;

qPhysuudd2211
=
-0.5*qdd11*quu22 + pow2(qud21)
;

qPhysuudd2212
=
qud21*qud22 - 0.5*qdd12*quu22
;

qPhysuudd2213
=
qud21*qud23 - 0.5*qdd13*quu22
;

qPhysuudd2222
=
-0.5*qdd22*quu22 + pow2(qud22)
;

qPhysuudd2223
=
qud22*qud23 - 0.5*qdd23*quu22
;

qPhysuudd2233
=
-0.5*qdd33*quu22 + pow2(qud23)
;

qPhysuudd2311
=
qud21*qud31 - 0.5*qdd11*quu23
;

qPhysuudd2312
=
0.5*(qud22*qud31 + qud21*qud32 - qdd12*quu23)
;

qPhysuudd2313
=
0.5*(qud23*qud31 + qud21*qud33 - qdd13*quu23)
;

qPhysuudd2322
=
qud22*qud32 - 0.5*qdd22*quu23
;

qPhysuudd2323
=
0.5*(qud23*qud32 + qud22*qud33 - qdd23*quu23)
;

qPhysuudd2333
=
qud23*qud33 - 0.5*qdd33*quu23
;

qPhysuudd3311
=
-0.5*qdd11*quu33 + pow2(qud31)
;

qPhysuudd3312
=
qud31*qud32 - 0.5*qdd12*quu33
;

qPhysuudd3313
=
qud31*qud33 - 0.5*qdd13*quu33
;

qPhysuudd3322
=
-0.5*qdd22*quu33 + pow2(qud32)
;

qPhysuudd3323
=
qud32*qud33 - 0.5*qdd23*quu33
;

qPhysuudd3333
=
-0.5*qdd33*quu33 + pow2(qud33)
;

K
=
Khat + 2.*Theta
;

dK1
=
dKhat1 + 2.*dTheta1
;

dK2
=
dKhat2 + 2.*dTheta2
;

dK3
=
dKhat3 + 2.*dTheta3
;

gammado111
=
0.5*dg111
;

gammado112
=
0.5*dg211
;

gammado113
=
0.5*dg311
;

gammado122
=
-0.5*dg122 + dg212
;

gammado123
=
0.5*(-dg123 + dg213 + dg312)
;

gammado133
=
-0.5*dg133 + dg313
;

gammado211
=
dg112 - 0.5*dg211
;

gammado212
=
0.5*dg122
;

gammado213
=
0.5*(dg123 - dg213 + dg312)
;

gammado222
=
0.5*dg222
;

gammado223
=
0.5*dg322
;

gammado233
=
-0.5*dg233 + dg323
;

gammado311
=
dg113 - 0.5*dg311
;

gammado312
=
0.5*(dg123 + dg213 - dg312)
;

gammado313
=
0.5*dg133
;

gammado322
=
dg223 - 0.5*dg322
;

gammado323
=
0.5*dg233
;

gammado333
=
0.5*dg333
;

gamma111
=
gammado111*ginv11 + gammado211*ginv12 + gammado311*ginv13
;

gamma112
=
gammado112*ginv11 + gammado212*ginv12 + gammado312*ginv13
;

gamma113
=
gammado113*ginv11 + gammado213*ginv12 + gammado313*ginv13
;

gamma122
=
gammado122*ginv11 + gammado222*ginv12 + gammado322*ginv13
;

gamma123
=
gammado123*ginv11 + gammado223*ginv12 + gammado323*ginv13
;

gamma133
=
gammado133*ginv11 + gammado233*ginv12 + gammado333*ginv13
;

gamma211
=
gammado111*ginv12 + gammado211*ginv22 + gammado311*ginv23
;

gamma212
=
gammado112*ginv12 + gammado212*ginv22 + gammado312*ginv23
;

gamma213
=
gammado113*ginv12 + gammado213*ginv22 + gammado313*ginv23
;

gamma222
=
gammado122*ginv12 + gammado222*ginv22 + gammado322*ginv23
;

gamma223
=
gammado123*ginv12 + gammado223*ginv22 + gammado323*ginv23
;

gamma233
=
gammado133*ginv12 + gammado233*ginv22 + gammado333*ginv23
;

gamma311
=
gammado111*ginv13 + gammado211*ginv23 + gammado311*ginv33
;

gamma312
=
gammado112*ginv13 + gammado212*ginv23 + gammado312*ginv33
;

gamma313
=
gammado113*ginv13 + gammado213*ginv23 + gammado313*ginv33
;

gamma322
=
gammado122*ginv13 + gammado222*ginv23 + gammado322*ginv33
;

gamma323
=
gammado123*ginv13 + gammado223*ginv23 + gammado323*ginv33
;

gamma333
=
gammado133*ginv13 + gammado233*ginv23 + gammado333*ginv33
;

Gfromg1
=
gamma111*ginv11 + gamma122*ginv22 + 
  2.*(gamma112*ginv12 + gamma113*ginv13 + gamma123*ginv23) + gamma133*ginv33
;

Gfromg2
=
gamma211*ginv11 + gamma222*ginv22 + 
  2.*(gamma212*ginv12 + gamma213*ginv13 + gamma223*ginv23) + gamma233*ginv33
;

Gfromg3
=
gamma311*ginv11 + gamma322*ginv22 + 
  2.*(gamma312*ginv12 + gamma313*ginv13 + gamma323*ginv23) + gamma333*ginv33
;

dGfromgdu11
=
-((dg122*(4.*dg112 + dg211) + 2.*dg112*dg212 + dg111*dg222)*
     Power(ginv12,3)) - (dg133*(4.*dg113 + dg311) + 2.*dg113*dg313 + 
     dg111*dg333)*Power(ginv13,3) - 2.*Power(ginv11,3)*pow2(dg111) + 
  (ddg1111 - dg111*((8.*dg112 + 2.*dg211)*ginv12 + 
        (8.*dg113 + 2.*dg311)*ginv13) - 
     (dg113*(4.*dg112 + dg211) + dg112*dg311 + dg111*(dg213 + dg312))*
      ginv23 - ginv22*(dg112*dg211 + dg111*dg212 + 2.*pow2(dg112)) - 
     ginv33*(dg113*dg311 + dg111*dg313 + 2.*pow2(dg113)))*pow2(ginv11) + 
  (ddg1122 + ddg1212 - (dg123*(8.*dg112 + 2.*dg211) + 
        dg113*(4.*dg122 + 2.*dg212) + dg122*dg311 + 
        2.*(dg111*dg223 + dg112*(dg213 + dg312)) + dg111*dg322)*ginv13 - 
     (dg123*(4.*dg122 + 2.*dg212) + 
        2.*(dg113*dg222 + dg122*(dg213 + dg312) + dg112*(dg223 + dg322)))*
      ginv23 - ginv22*(3.*(dg122*dg212 + dg112*dg222) + 2.*pow2(dg122)) - 
     ginv33*(dg123*(dg213 + dg312) + dg122*dg313 + dg113*(dg223 + dg322) + 
        dg112*dg323 + 2.*pow2(dg123)))*pow2(ginv12) + 
  (ddg1133 + ddg1313 - (dg133*(4.*dg123 + 2.*(dg213 + dg312)) + 
        2.*(dg123*dg313 + dg113*(dg233 + dg323) + dg112*dg333))*ginv23 - 
     ginv22*(dg133*dg212 + dg113*dg223 + dg123*(dg213 + dg312) + 
        dg112*(dg233 + dg323) + 2.*pow2(dg123)) - 
     ginv33*(3.*(dg133*dg313 + dg113*dg333) + 2.*pow2(dg133)))*pow2(ginv13) \
+ ginv13*(ddg1333*ginv33 + ginv22*
      (ddg1223 - (dg133*dg222 + dg123*(4.*dg223 + dg322) + 
           dg122*(dg233 + dg323))*ginv23 - 
        (dg133*dg223 + dg123*(dg233 + 2.*dg323))*ginv33) + 
     ginv23*(ddg1233 + ddg1323 - 
        (dg133*(2.*dg233 + 3.*dg323) + 3.*dg123*dg333)*ginv33) - 
     (dg123*dg222 + dg122*dg223)*pow2(ginv22) - 
     (dg133*dg322 + 2.*(dg133*dg223 + dg123*(dg233 + dg323)) + 
        dg122*dg333)*pow2(ginv23) - 2.*dg133*dg333*pow2(ginv33)) + 
  ginv11*(ddg1313*ginv33 + ginv12*
      (2.*ddg1112 + ddg1211 - 
        (dg113*(12.*dg112 + 3.*dg211) + 3.*dg112*dg311 + 
           dg111*(8.*dg123 + 3.*(dg213 + dg312)))*ginv13 - 
        (dg122*(4.*dg112 + dg211) + 6.*dg112*dg212 + dg111*dg222)*ginv22 - 
        (dg123*dg211 + dg122*dg311 + 
           4.*(dg113*(dg122 + dg212) + dg112*(dg123 + dg213 + dg312)) + 
           dg111*(dg223 + dg322))*ginv23 - 
        (dg123*dg311 + dg113*(4.*dg123 + 2.*(dg213 + dg312)) + 
           2.*dg112*dg313 + dg111*dg323)*ginv33) + 
     ginv22*(ddg1212 - (dg113*dg222 + 2.*(dg123*dg212 + dg112*dg223) + 
           dg122*(dg213 + dg312) + dg112*dg322)*ginv23 - 
        (dg113*dg223 + dg123*(dg213 + dg312) + dg112*dg323)*ginv33) + 
     ginv13*(2.*ddg1113 + ddg1311 - 
        (dg123*(4.*dg112 + dg211) + dg111*dg223 + 
           2.*(dg113*dg212 + dg112*(dg213 + dg312)))*ginv22 - 
        (dg133*dg211 + dg123*dg311 + 
           4.*(dg113*(dg123 + dg213 + dg312) + dg112*(dg133 + dg313)) + 
           dg111*(dg233 + dg323))*ginv23 - 
        (dg133*(4.*dg113 + dg311) + 6.*dg113*dg313 + dg111*dg333)*ginv33) + 
     ginv23*(ddg1213 + ddg1312 - 
        (dg133*(dg213 + dg312) + 2.*dg123*dg313 + 
           dg113*(dg233 + 2.*dg323) + dg112*dg333)*ginv33) - 
     (3.*dg112*dg211 + dg111*(4.*dg122 + 3.*dg212) + 6.*pow2(dg112))*
      pow2(ginv12) - (3.*dg113*dg311 + dg111*(4.*dg133 + 3.*dg313) + 
        6.*pow2(dg113))*pow2(ginv13) - 
     (dg122*dg212 + dg112*dg222)*pow2(ginv22) - 
     (dg133*dg212 + dg123*(dg213 + dg312) + dg122*dg313 + 
        dg113*(dg223 + dg322) + dg112*(dg233 + dg323))*pow2(ginv23) - 
     (dg133*dg313 + dg113*dg333)*pow2(ginv33)) + 
  ginv12*(ddg1323*ginv33 + ginv22*
      (ddg1222 - (3.*(dg123*dg222 + dg122*dg223) + 2.*dg122*dg322)*ginv23 - 
        (dg123*(2.*dg223 + dg322) + dg122*dg323)*ginv33) + 
     ginv23*(ddg1223 + ddg1322 - 
        (dg133*(dg223 + dg322) + dg123*(dg233 + 4.*dg323) + dg122*dg333)*
         ginv33) + ginv13*(2.*ddg1123 + ddg1213 + ddg1312 - 
        (dg113*dg222 + 4.*(dg123*(dg122 + dg212) + dg112*dg223) + 
           dg122*(dg213 + dg312) + dg112*dg322)*ginv22 - 
        (dg133*(4.*dg123 + dg213 + dg312) + 4.*dg123*dg313 + 
           dg113*(dg233 + 4.*dg323) + dg112*dg333)*ginv33 - 
        ginv23*(2.*(dg133*dg212 + dg112*dg233 + dg122*dg313 + 
              dg113*dg322) + 4.*
            (dg122*dg133 + dg113*dg223 + dg123*(dg213 + dg312) + 
              dg112*dg323 + pow2(dg123)))) - 
     (dg133*(4.*dg112 + dg211) + dg113*(8.*dg123 + 2.*(dg213 + dg312)) + 
        2.*(dg123*dg311 + dg112*dg313) + dg111*(dg233 + 2.*dg323))*
      pow2(ginv13) - 2.*dg122*dg222*pow2(ginv22) - 
     (dg133*dg222 + 2.*dg123*(dg223 + dg322) + dg122*(dg233 + 2.*dg323))*
      pow2(ginv23) - (dg133*dg323 + dg123*dg333)*pow2(ginv33))
;

dGfromgdu12
=
-((dg133*dg322 + 2.*(dg133*dg223 + dg123*(dg233 + dg323)) + dg122*dg333)*
     Power(ginv23,3)) - 2.*(dg122*dg222*Power(ginv22,3) + 
     Power(ginv12,3)*(dg112*dg211 + dg111*(dg122 + dg212) + pow2(dg112)) + 
     (dg111*(dg112*ginv22 + dg113*ginv23) + ginv12*pow2(dg111))*pow2(ginv11)\
) + (ddg1112 + ddg1211 - (4.*(dg112*dg113 + dg111*dg123) + 
        2.*(dg113*dg211 + dg112*dg311 + dg111*(dg213 + dg312)))*ginv13 - 
     (dg122*(6.*dg112 + 2.*dg211) + 6.*dg112*dg212 + 2.*dg111*dg222)*
      ginv22 - (4.*(dg113*(dg122 + dg212) + dg112*(dg123 + dg213)) + 
        dg122*dg311 + 2.*(dg123*dg211 + dg111*dg223 + dg112*dg312) + 
        dg111*dg322)*ginv23 - 
     (dg123*dg311 + dg113*(2.*(dg123 + dg213) + dg312) + dg112*dg313 + 
        dg111*dg323)*ginv33)*pow2(ginv12) - 
  ((2.*(dg113*dg123 + dg112*dg133) + dg123*dg311 + dg113*dg312 + 
        dg112*dg313 + dg111*dg323)*ginv22 + 
     (dg133*(4.*dg113 + dg311) + 2.*dg113*dg313 + dg111*dg333)*ginv23)*
   pow2(ginv13) + (ddg1222 - (4.*(dg123*dg222 + dg122*dg223) + 
        2.*dg122*dg322)*ginv23 - 
     (dg123*(2.*dg223 + dg322) + dg122*dg323)*ginv33)*pow2(ginv22) + 
  (ddg1233 + ddg1323 - (dg133*(2.*dg233 + 3.*dg323) + 3.*dg123*dg333)*
      ginv33)*pow2(ginv23) + ginv11*
   (ginv23*(ddg1113 - 2.*dg113*(dg133 + dg313)*ginv33) + 
     ginv22*(ddg1112 - (dg112*(4.*dg123 + 2.*dg213) + 
           2.*(dg113*(dg122 + dg212) + dg112*dg312))*ginv23 - 
        (dg113*(2.*dg123 + dg312) + dg112*dg313)*ginv33) + 
     ginv12*(ddg1111 - dg111*(6.*dg113 + 2.*dg311)*ginv13 - 
        (dg113*(8.*dg112 + 2.*dg211) + dg112*dg311 + 
           dg111*(2.*(dg123 + dg213) + dg312))*ginv23 - 
        ginv22*(2.*(dg112*dg211 + dg111*(dg122 + dg212)) + 
           6.*pow2(dg112)) - ginv33*
         (dg113*dg311 + dg111*dg313 + 2.*pow2(dg113))) - 
     ginv13*((dg112*(4.*dg113 + dg311) + dg111*(2.*dg123 + dg312))*
         ginv22 + ginv23*(dg113*dg311 + dg111*(2.*dg133 + dg313) + 
           4.*pow2(dg113))) - dg111*(6.*dg112 + 2.*dg211)*pow2(ginv12) - 
     2.*dg112*(dg122 + dg212)*pow2(ginv22) - 
     (2.*(dg112*dg133 + dg113*(dg123 + dg213)) + dg113*dg312 + dg112*dg313)*
      pow2(ginv23)) + ginv13*(ginv22*
      (ddg1123 + ddg1312 - (dg133*(2.*dg123 + dg312) + 
           2.*(dg123*dg313 + dg113*dg323) + dg112*dg333)*ginv33 - 
        ginv23*(2.*(dg133*(dg122 + dg212) + dg123*dg213 + dg113*dg223 + 
              dg112*dg233) + dg122*dg313 + dg113*dg322 + 
           4.*(dg123*dg312 + dg112*dg323 + pow2(dg123)))) + 
     ginv23*(ddg1133 + ddg1313 - 
        ginv33*(3.*(dg133*dg313 + dg113*dg333) + 2.*pow2(dg133))) - 
     (2.*(dg123*(dg122 + dg212) + dg112*dg223) + dg122*dg312 + 
        dg112*dg322)*pow2(ginv22) - 
     (dg133*(4.*dg123 + 2.*(dg213 + dg312)) + 
        2.*(dg123*dg313 + dg113*(dg233 + dg323) + dg112*dg333))*pow2(ginv23)\
) + ginv23*(ddg1333*ginv33 - 2.*dg133*dg333*pow2(ginv33)) + 
  ginv12*(ddg1313*ginv33 + ginv13*
      (ddg1113 + ddg1311 - (2.*
            (dg123*dg211 + dg113*(dg122 + dg212) + dg111*dg223) + 
           dg122*dg311 + dg112*(8.*dg123 + 2.*dg213 + 4.*dg312) + 
           dg111*dg322)*ginv22 - 
        (dg133*(4.*dg112 + 2.*dg211) + 
           dg113*(8.*dg123 + 4.*(dg213 + dg312)) + 4.*dg112*dg313 + 
           2.*(dg123*dg311 + dg111*(dg233 + dg323)))*ginv23 - 
        (dg133*(2.*dg113 + dg311) + 4.*dg113*dg313 + dg111*dg333)*ginv33) + 
     ginv23*(ddg1123 + 2.*ddg1213 + ddg1312 - 
        (2.*(dg133*(dg123 + dg213) + dg113*dg233) + dg133*dg312 + 
           4.*(dg123*dg313 + dg113*dg323) + dg112*dg333)*ginv33) + 
     ginv22*(ddg1122 + 2.*ddg1212 - 
        (4.*(dg122*dg213 + dg113*dg222) + 
           6.*(dg123*(dg122 + dg212) + dg112*dg223) + 
           3.*(dg122*dg312 + dg112*dg322))*ginv23 - 
        ginv33*(dg122*dg313 + dg113*dg322 + 
           2.*(dg113*dg223 + dg123*(dg213 + dg312) + dg112*dg323 + 
              pow2(dg123)))) - 
     2.*(dg113*dg311 + dg111*(dg133 + dg313) + pow2(dg113))*pow2(ginv13) - 
     (4.*(dg122*dg212 + dg112*dg222) + 2.*pow2(dg122))*pow2(ginv22) - 
     (4.*(dg123*dg213 + dg113*dg223) + 
        2.*(dg133*(dg122 + dg212) + dg123*dg312 + dg122*dg313 + 
           dg113*dg322 + dg112*(dg233 + dg323) + pow2(dg123)))*pow2(ginv23) \
- (dg133*dg313 + dg113*dg333)*pow2(ginv33)) + 
  ginv22*(ddg1323*ginv33 + ginv23*
      (2.*ddg1223 + ddg1322 - (2.*(dg133*dg223 + dg123*dg233) + 
           dg133*dg322 + 6.*dg123*dg323 + dg122*dg333)*ginv33) - 
     (2.*(dg133*dg222 + dg122*dg233) + dg123*(6.*dg223 + 3.*dg322) + 
        3.*dg122*dg323)*pow2(ginv23) - 
     (dg133*dg323 + dg123*dg333)*pow2(ginv33))
;

dGfromgdu13
=
-((dg133*dg222 + 2.*dg123*(dg223 + dg322) + dg122*(dg233 + 2.*dg323))*
     Power(ginv23,3)) - 2.*(dg133*dg333*Power(ginv33,3) + 
     Power(ginv13,3)*(dg113*dg311 + dg111*(dg133 + dg313) + pow2(dg113)) + 
     (dg111*(dg112*ginv23 + dg113*ginv33) + ginv13*pow2(dg111))*pow2(ginv11)\
) - ((dg122*(4.*dg112 + dg211) + 2.*dg112*dg212 + dg111*dg222)*ginv23 + 
     (2.*(dg113*dg122 + dg112*dg123) + dg123*dg211 + dg113*dg212 + 
        dg112*dg213 + dg111*dg223)*ginv33 + 
     2.*ginv13*(dg112*dg211 + dg111*(dg122 + dg212) + pow2(dg112)))*
   pow2(ginv12) + (ddg1113 + ddg1311 - 
     (dg123*(2.*dg112 + dg211) + dg113*dg212 + dg111*dg223 + 
        dg112*(dg213 + 2.*dg312))*ginv22 - 
     (dg133*dg211 + 2.*(dg113*dg213 + dg123*dg311) + 
        4.*(dg113*(dg123 + dg312) + dg112*(dg133 + dg313)) + 
        dg111*(dg233 + 2.*dg323))*ginv23 - 
     (dg133*(6.*dg113 + 2.*dg311) + 6.*dg113*dg313 + 2.*dg111*dg333)*ginv33\
)*pow2(ginv13) - (2.*dg122*dg222*ginv23 + 
     (dg123*dg222 + dg122*dg223)*ginv33)*pow2(ginv22) + 
  (ddg1223 + ddg1322 - (3.*(dg133*dg223 + dg123*dg233) + 6.*dg123*dg323 + 
        2.*(dg133*dg322 + dg122*dg333))*ginv33)*pow2(ginv23) + 
  ddg1333*pow2(ginv33) + ginv11*
   (ddg1113*ginv33 - ginv22*(2.*dg112*(dg122 + dg212)*ginv23 + 
        (dg113*dg212 + dg112*(2.*dg123 + dg213))*ginv33) + 
     ginv23*(ddg1112 - (dg113*(4.*dg123 + 2.*dg213) + 
           2.*(dg113*dg312 + dg112*(dg133 + dg313)))*ginv33) - 
     ginv12*(dg111*(6.*dg112 + 2.*dg211)*ginv13 + 
        (dg113*(4.*dg112 + dg211) + dg111*(2.*dg123 + dg213))*ginv33 + 
        ginv23*(dg112*dg211 + dg111*(2.*dg122 + dg212) + 4.*pow2(dg112))) + 
     ginv13*(ddg1111 - (dg113*(8.*dg112 + dg211) + 2.*dg112*dg311 + 
           dg111*(dg213 + 2.*(dg123 + dg312)))*ginv23 - 
        ginv22*(dg112*dg211 + dg111*dg212 + 2.*pow2(dg112)) - 
        ginv33*(2.*(dg113*dg311 + dg111*(dg133 + dg313)) + 6.*pow2(dg113))) \
- dg111*(6.*dg113 + 2.*dg311)*pow2(ginv13) - 
     (dg113*dg212 + dg112*dg213 + 
        2.*(dg113*dg122 + dg112*(dg123 + dg312)))*pow2(ginv23) - 
     2.*dg113*(dg133 + dg313)*pow2(ginv33)) + 
  ginv12*((ddg1123 + ddg1213)*ginv33 + 
     ginv13*(ddg1112 + ddg1211 - 
        (dg122*(2.*dg112 + dg211) + 4.*dg112*dg212 + dg111*dg222)*ginv22 - 
        (dg123*(8.*dg112 + 2.*dg211) + 
           4.*(dg113*(dg122 + dg212) + dg112*(dg213 + dg312)) + 
           2.*(dg122*dg311 + dg111*(dg223 + dg322)))*ginv23 - 
        (dg133*(2.*dg112 + dg211) + 
           dg113*(8.*dg123 + 4.*dg213 + 2.*dg312) + 
           2.*(dg123*dg311 + dg112*dg313) + dg111*(dg233 + 2.*dg323))*
         ginv33) - ginv22*((dg122*dg213 + dg113*dg222 + 
           2.*(dg123*(dg122 + dg212) + dg112*dg223))*ginv33 + 
        ginv23*(3.*(dg122*dg212 + dg112*dg222) + 2.*pow2(dg122))) + 
     ginv23*(ddg1122 + ddg1212 - 
        ginv33*(dg133*(2.*dg122 + dg212) + 
           2.*(dg123*dg312 + dg122*dg313 + dg113*dg322) + 
           dg112*(dg233 + 2.*dg323) + 
           4.*(dg123*dg213 + dg113*dg223 + pow2(dg123)))) - 
     (4.*(dg112*dg113 + dg111*dg123) + 
        2.*(dg113*dg211 + dg112*dg311 + dg111*(dg213 + dg312)))*
      pow2(ginv13) - (dg123*(4.*dg122 + 2.*dg212) + 
        2.*(dg113*dg222 + dg122*(dg213 + dg312) + dg112*(dg223 + dg322)))*
      pow2(ginv23) - (dg133*(2.*dg123 + dg213) + 2.*dg123*dg313 + 
        dg113*(dg233 + 2.*dg323))*pow2(ginv33)) + 
  ginv22*(ddg1223*ginv33 + ginv23*
      (ddg1222 - (dg133*dg222 + dg123*(6.*dg223 + 2.*dg322) + 
           dg122*(dg233 + 2.*dg323))*ginv33) - 
     (3.*(dg123*dg222 + dg122*dg223) + 2.*dg122*dg322)*pow2(ginv23) - 
     (dg133*dg223 + dg123*(dg233 + 2.*dg323))*pow2(ginv33)) + 
  ginv23*((ddg1233 + 2.*ddg1323)*ginv33 - 
     (dg133*(2.*dg233 + 4.*dg323) + 4.*dg123*dg333)*pow2(ginv33)) + 
  ginv13*((ddg1133 + 2.*ddg1313)*ginv33 + 
     ginv23*(ddg1123 + ddg1213 + 2.*ddg1312 - 
        (dg133*(6.*dg123 + 3.*dg213 + 4.*dg312) + 6.*dg123*dg313 + 
           dg113*(3.*dg233 + 6.*dg323) + 4.*dg112*dg333)*ginv33) + 
     ginv22*(ddg1212 - (dg123*(2.*dg122 + 4.*dg212) + dg113*dg222 + 
           dg122*(dg213 + 2.*dg312) + dg112*(4.*dg223 + 2.*dg322))*ginv23 - 
        ginv33*(dg133*dg212 + dg112*(dg233 + 2.*dg323) + 
           2.*(dg113*dg223 + dg123*(dg213 + dg312) + pow2(dg123)))) - 
     (dg122*dg212 + dg112*dg222)*pow2(ginv22) - 
     (4.*(dg123*dg312 + dg112*dg323) + 
        2.*(dg133*(dg122 + dg212) + dg123*dg213 + dg112*dg233 + 
           dg122*dg313 + dg113*(dg223 + dg322) + pow2(dg123)))*pow2(ginv23) \
- (4.*(dg133*dg313 + dg113*dg333) + 2.*pow2(dg133))*pow2(ginv33))
;

dGfromgdu21
=
-((dg233*dg311 + 2.*(dg113*dg233 + dg213*(dg133 + dg313)) + dg211*dg333)*
     Power(ginv13,3)) - 2.*(dg111*dg211*Power(ginv11,3) + 
     Power(ginv12,3)*(dg122*dg212 + (dg112 + dg211)*dg222 + pow2(dg212))) + 
  (ddg1211 - (4.*(dg113*dg211 + dg111*dg213) + 2.*dg211*dg311)*ginv13 - 
     2.*(dg112 + dg211)*dg212*ginv22 - 
     (2.*(dg113*dg212 + (dg112 + dg211)*dg213) + dg212*dg311 + 
        dg211*dg312)*ginv23 - 
     (dg213*(2.*dg113 + dg311) + dg211*dg313)*ginv33 - 
     ginv12*(4.*(dg112*dg211 + dg111*dg212) + 2.*pow2(dg211)))*pow2(ginv11) \
+ (ddg1222 + ddg2212 - (4.*(dg212*(dg123 + dg213) + 
           (dg112 + dg211)*dg223) + dg222*dg311 + 
        2.*(dg122*dg213 + dg113*dg222 + dg212*dg312) + dg211*dg322)*ginv13 \
- (2.*dg122 + 6.*dg212)*dg222*ginv22 - 
     ((2.*dg122 + 4.*dg212)*dg223 + 
        dg222*(4.*dg213 + 2.*(dg123 + dg312)) + 2.*dg212*dg322)*ginv23 - 
     (dg223*(2.*(dg123 + dg213) + dg312) + dg222*dg313 + dg213*dg322 + 
        dg212*dg323)*ginv33)*pow2(ginv12) + 
  (ddg1233 + ddg2313 - (2.*((dg123 + dg213)*dg223 + dg212*dg233) + 
        dg223*dg312 + dg212*dg323)*ginv22 - 
     (dg233*(4.*dg213 + 2.*dg312) + 
        2.*(dg123*dg233 + dg223*(dg133 + dg313) + dg213*dg323 + 
           dg212*dg333))*ginv23 - 
     (dg233*(2.*dg133 + 3.*dg313) + 3.*dg213*dg333)*ginv33)*pow2(ginv13) + 
  ginv11*(ddg2313*ginv33 + ginv22*
      (ddg2212 - (dg222*(2.*dg213 + dg312) + dg212*(4.*dg223 + dg322))*
         ginv23 - (dg223*(2.*dg213 + dg312) + dg212*dg323)*ginv33) + 
     ginv23*(ddg2213 + ddg2312 - 
        (dg233*(2.*dg213 + dg312) + 2.*(dg223*dg313 + dg213*dg323) + 
           dg212*dg333)*ginv33) + 
     ginv13*(2.*ddg1213 + ddg2311 - 
        (2.*(dg112 + dg211)*dg223 + 
           dg212*(4.*dg213 + 2.*(dg123 + dg312)))*ginv22 - 
        (2.*(dg133*dg213 + dg113*dg233) + dg233*dg311 + 6.*dg213*dg313 + 
           dg211*dg333)*ginv33 - 
        ginv23*(2.*(dg133*dg212 + dg123*dg213 + dg113*dg223 + 
              (dg112 + dg211)*dg233) + dg223*dg311 + dg211*dg323 + 
           4.*(dg213*dg312 + dg212*dg313 + pow2(dg213)))) + 
     ginv12*(2.*ddg1212 + ddg2211 - 
        (6.*(dg113*dg212 + dg112*dg213) + 4.*dg111*dg223 + 
           3.*dg212*dg311 + dg211*(4.*dg123 + 6.*dg213 + 3.*dg312))*ginv13 \
- (2.*(dg123*dg212 + dg122*dg213 + dg113*dg222 + 
              (dg112 + dg211)*dg223) + dg222*dg311 + 
           dg212*(8.*dg213 + 4.*dg312) + dg211*dg322)*ginv23 - 
        ginv22*(2.*(dg122*dg212 + (dg112 + dg211)*dg222) + 
           6.*pow2(dg212)) - ginv33*
         (dg223*dg311 + dg211*dg323 + 
           2.*(dg113*dg223 + dg213*(dg123 + dg312) + dg212*dg313 + 
              pow2(dg213)))) - 
     (6.*dg112*dg212 + dg211*(2.*dg122 + 6.*dg212) + 2.*dg111*dg222)*
      pow2(ginv12) - (2.*(dg133*dg211 + dg111*dg233) + 
        dg213*(6.*dg113 + 3.*dg311) + 3.*dg211*dg313)*pow2(ginv13) - 
     2.*dg212*dg222*pow2(ginv22) - 
     (2.*(dg213*dg223 + dg212*dg233) + dg223*dg312 + dg222*dg313 + 
        dg213*dg322 + dg212*dg323)*pow2(ginv23) - 
     (dg233*dg313 + dg213*dg333)*pow2(ginv33)) + 
  ginv12*(ddg2323*ginv33 + ginv13*
      (2.*ddg1223 + ddg2213 + ddg2312 - 
        (2.*((dg123 + dg213)*dg222 + dg122*dg223) + dg222*dg312 + 
           dg212*(8.*dg223 + dg322))*ginv22 - 
        (dg223*(8.*dg213 + 4.*(dg123 + dg312)) + 
           2.*(dg122*dg233 + dg222*(dg133 + dg313) + dg213*dg322) + 
           4.*dg212*(dg233 + dg323))*ginv23 - 
        (2.*(dg133*dg223 + (dg123 + dg213)*dg233) + dg233*dg312 + 
           4.*(dg223*dg313 + dg213*dg323) + dg212*dg333)*ginv33) + 
     ginv23*(ddg2223 + ddg2322 - 
        (dg233*(2.*dg223 + dg322) + 4.*dg223*dg323 + dg222*dg333)*ginv33) + 
     ginv22*(ddg2222 - dg222*(6.*dg223 + 2.*dg322)*ginv23 - 
        ginv33*(dg223*dg322 + dg222*dg323 + 2.*pow2(dg223))) - 
     (4.*(dg123*dg213 + dg113*dg223) + 
        2.*((dg112 + dg211)*dg233 + dg223*dg311 + dg213*dg312 + 
           dg212*(dg133 + dg313) + dg211*dg323 + pow2(dg213)))*pow2(ginv13) \
- 2.*(pow2(dg222)*pow2(ginv22) + 
        (dg223*dg322 + dg222*(dg233 + dg323) + pow2(dg223))*pow2(ginv23)) - 
     (dg233*dg323 + dg223*dg333)*pow2(ginv33)) + 
  ginv13*(ddg2333*ginv33 + ginv22*
      (ddg2223 - 2.*dg223*(dg233 + dg323)*ginv33 - 
        ginv23*(dg223*dg322 + dg222*(2.*dg233 + dg323) + 4.*pow2(dg223))) + 
     ginv23*(ddg2233 + ddg2323 - 
        ginv33*(3.*(dg233*dg323 + dg223*dg333) + 2.*pow2(dg233))) - 
     (dg233*(4.*dg223 + dg322) + 2.*dg223*dg323 + dg222*dg333)*
      pow2(ginv23) - 2.*(dg222*dg223*pow2(ginv22) + dg233*dg333*pow2(ginv33))\
)
;

dGfromgdu22
=
-((2.*dg112*dg212 + dg211*(dg122 + 4.*dg212) + dg111*dg222)*
     Power(ginv12,3)) - (dg233*(4.*dg223 + dg322) + 2.*dg223*dg323 + 
     dg222*dg333)*Power(ginv23,3) - 2.*Power(ginv22,3)*pow2(dg222) - 
  (2.*dg111*dg211*ginv12 + (dg112*dg211 + dg111*dg212)*ginv22 + 
     (dg113*dg211 + dg111*dg213)*ginv23)*pow2(ginv11) + 
  (ddg1212 + ddg2211 - (2.*(dg123*dg211 + dg112*dg213 + dg111*dg223 + 
           dg212*(dg113 + dg311)) + dg211*(4.*dg213 + 2.*dg312))*ginv13 - 
     (2.*(dg123*dg212 + dg122*dg213 + dg113*dg222 + dg112*dg223) + 
        dg222*dg311 + dg212*(8.*dg213 + 2.*dg312) + 
        dg211*(4.*dg223 + dg322))*ginv23 - 
     ginv22*(4.*dg211*dg222 + 3.*(dg122*dg212 + dg112*dg222) + 
        6.*pow2(dg212)) - ginv33*
      (dg223*(dg113 + dg311) + dg213*(dg123 + dg312) + dg212*dg313 + 
        dg211*dg323 + 2.*pow2(dg213)))*pow2(ginv12) - 
  ((dg112*dg233 + dg223*(dg113 + dg311) + dg213*(dg123 + dg312) + 
        dg212*(dg133 + dg313) + dg211*dg323)*ginv22 + 
     (dg233*dg311 + 2.*(dg113*dg233 + dg213*(dg133 + dg313)) + 
        dg211*dg333)*ginv23)*pow2(ginv13) + 
  (ddg2222 - dg222*(8.*dg223 + 2.*dg322)*ginv23 - 
     ginv33*(dg223*dg322 + dg222*dg323 + 2.*pow2(dg223)))*pow2(ginv22) + 
  (ddg2233 + ddg2323 - ginv33*
      (3.*(dg233*dg323 + dg223*dg333) + 2.*pow2(dg233)))*pow2(ginv23) + 
  ginv13*(ginv22*(ddg1223 + ddg2312 - 
        (dg122*dg233 + dg222*(dg133 + dg313) + dg213*dg322 + 
           4.*(dg223*(dg123 + dg213 + dg312) + dg212*(dg233 + dg323)))*
         ginv23 - (dg233*(dg123 + dg312) + dg223*(dg133 + 2.*dg313) + 
           2.*dg213*dg323 + dg212*dg333)*ginv33) + 
     ginv23*(ddg1233 + ddg2313 - 
        (dg233*(2.*dg133 + 3.*dg313) + 3.*dg213*dg333)*ginv33) - 
     ((dg122 + 4.*dg212)*dg223 + dg222*(dg123 + dg312) + dg212*dg322)*
      pow2(ginv22) - (dg233*(4.*dg213 + 2.*dg312) + 
        2.*(dg123*dg233 + dg223*(dg133 + dg313) + dg213*dg323 + 
           dg212*dg333))*pow2(ginv23)) + 
  ginv11*(-(ginv13*((2.*(dg113*dg212 + dg112*dg213) + dg111*dg223 + 
             dg212*dg311 + dg211*(dg123 + dg312))*ginv22 + 
          (dg111*dg233 + dg213*(4.*dg113 + dg311) + dg211*(dg133 + dg313))*
           ginv23)) + ginv12*(ddg1211 - 
        (3.*(dg113*dg211 + dg111*dg213) + 2.*dg211*dg311)*ginv13 - 
        (6.*dg112*dg212 + dg211*(dg122 + 4.*dg212) + dg111*dg222)*ginv22 - 
        (4.*(dg113*dg212 + dg112*dg213) + dg111*dg223 + dg212*dg311 + 
           dg211*(dg123 + 4.*dg213 + dg312))*ginv23 - 
        (dg213*(2.*dg113 + dg311) + dg211*dg313)*ginv33) + 
     ginv22*(ddg1212 - (dg122*dg213 + dg113*dg222 + 2.*dg112*dg223 + 
           dg212*(4.*dg213 + 2.*(dg123 + dg312)))*ginv23 - 
        (dg113*dg223 + dg213*(dg123 + dg312) + dg212*dg313)*ginv33) + 
     ginv23*(ddg1213 - (dg113*dg233 + dg213*(dg133 + 2.*dg313))*ginv33) - 
     (3.*(dg112*dg211 + dg111*dg212) + 2.*pow2(dg211))*pow2(ginv12) - 
     (dg122*dg212 + dg112*dg222 + 2.*pow2(dg212))*pow2(ginv22) - 
     (dg113*dg223 + dg112*dg233 + dg213*(dg123 + dg312) + 
        dg212*(dg133 + dg313) + 2.*pow2(dg213))*pow2(ginv23)) + 
  ginv23*(ddg2333*ginv33 - 2.*dg233*dg333*pow2(ginv33)) + 
  ginv12*(ddg2313*ginv33 + ginv22*
      (ddg1222 + 2.*ddg2212 - 
        ((3.*dg122 + 12.*dg212)*dg223 + 
           dg222*(8.*dg213 + 3.*(dg123 + dg312)) + 3.*dg212*dg322)*ginv23 \
- (dg223*(4.*dg213 + 2.*(dg123 + dg312)) + dg222*dg313 + dg213*dg322 + 
           2.*dg212*dg323)*ginv33) + 
     ginv23*(ddg1223 + 2.*ddg2213 + ddg2312 - 
        (dg233*(dg123 + 4.*dg213 + dg312) + dg223*(dg133 + 4.*dg313) + 
           4.*dg213*dg323 + dg212*dg333)*ginv33) + 
     ginv13*(ddg1213 + ddg2311 - 
        (dg122*dg213 + dg222*(dg113 + dg311) + 
           4.*((dg112 + dg211)*dg223 + dg212*(dg123 + dg213 + dg312)) + 
           dg211*dg322)*ginv22 - 
        (dg233*(dg113 + dg311) + dg213*(dg133 + 4.*dg313) + dg211*dg333)*
         ginv33 - ginv23*(2.*(dg133*dg212 + dg112*dg233 + dg223*dg311 + 
              dg211*dg323) + 4.*
            (dg113*dg223 + dg211*dg233 + dg213*(dg123 + dg312) + 
              dg212*dg313 + pow2(dg213)))) - 
     (dg111*dg233 + 2.*dg213*(dg113 + dg311) + dg211*(dg133 + 2.*dg313))*
      pow2(ginv13) - (2.*dg122 + 8.*dg212)*dg222*pow2(ginv22) - 
     ((dg122 + 4.*dg212)*dg233 + dg223*(8.*dg213 + 2.*(dg123 + dg312)) + 
        dg222*(dg133 + 2.*dg313) + 2.*(dg213*dg322 + dg212*dg323))*
      pow2(ginv23) - (dg233*dg313 + dg213*dg333)*pow2(ginv33)) + 
  ginv22*(ddg2323*ginv33 + ginv23*
      (2.*ddg2223 + ddg2322 - (dg233*(4.*dg223 + dg322) + 6.*dg223*dg323 + 
           dg222*dg333)*ginv33) - 
     (3.*dg223*dg322 + dg222*(4.*dg233 + 3.*dg323) + 6.*pow2(dg223))*
      pow2(ginv23) - (dg233*dg323 + dg223*dg333)*pow2(ginv33))
;

dGfromgdu23
=
-((dg111*dg233 + 2.*dg213*(dg113 + dg311) + dg211*(dg133 + 2.*dg313))*
     Power(ginv13,3)) - (2.*dg111*dg211*ginv13 + 
     (dg112*dg211 + dg111*dg212)*ginv23 + 
     (dg113*dg211 + dg111*dg213)*ginv33)*pow2(ginv11) - 
  ((2.*dg112*dg212 + dg211*(dg122 + 4.*dg212) + dg111*dg222)*ginv13 + 
     (dg122*dg213 + dg212*(dg123 + 2.*dg213) + dg113*dg222 + 
        (dg112 + 2.*dg211)*dg223)*ginv33 + 
     2.*ginv23*(dg122*dg212 + (dg112 + dg211)*dg222 + pow2(dg212)))*
   pow2(ginv12) + (ddg1213 + ddg2311 - 
     ((dg112 + 2.*dg211)*dg223 + dg212*(dg123 + 2.*(dg213 + dg312)))*
      ginv22 - (3.*(dg133*dg213 + dg113*dg233) + 6.*dg213*dg313 + 
        2.*(dg233*dg311 + dg211*dg333))*ginv33 - 
     ginv23*(4.*(dg213*dg312 + dg212*dg313) + 
        2.*(dg133*dg212 + dg123*dg213 + (dg112 + dg211)*dg233 + 
           dg223*(dg113 + dg311) + dg211*dg323 + pow2(dg213))))*pow2(ginv13) \
- 2.*(dg233*dg333*Power(ginv33,3) + 
     Power(ginv23,3)*(dg223*dg322 + dg222*(dg233 + dg323) + pow2(dg223)) + 
     (dg222*dg223*ginv33 + ginv23*pow2(dg222))*pow2(ginv22)) + 
  (ddg2223 + ddg2322 - (dg233*(6.*dg223 + 2.*dg322) + 6.*dg223*dg323 + 
        2.*dg222*dg333)*ginv33)*pow2(ginv23) + ddg2333*pow2(ginv33) + 
  ginv11*(ddg1213*ginv33 + ginv13*
      (ddg1211 - 2.*(dg112 + dg211)*dg212*ginv22 - 
        (4.*(dg113*dg212 + dg112*dg213) + dg111*dg223 + 2.*dg212*dg311 + 
           dg211*(dg123 + 2.*(dg213 + dg312)))*ginv23 - 
        (dg111*dg233 + dg213*(6.*dg113 + 2.*dg311) + 
           dg211*(dg133 + 2.*dg313))*ginv33) - 
     ginv12*((4.*dg112*dg212 + dg211*(dg122 + 2.*dg212) + dg111*dg222)*
         ginv23 + (dg211*(dg123 + 2.*dg213) + 
           2.*(dg113*dg212 + dg112*dg213) + dg111*dg223)*ginv33 + 
        ginv13*(3.*(dg112*dg211 + dg111*dg212) + 2.*pow2(dg211))) - 
     ginv22*((dg212*(dg123 + 2.*dg213) + dg112*dg223)*ginv33 + 
        ginv23*(dg122*dg212 + dg112*dg222 + 2.*pow2(dg212))) + 
     ginv23*(ddg1212 - ginv33*
         (dg112*dg233 + dg212*(dg133 + 2.*dg313) + 
           2.*(dg113*dg223 + dg213*(dg123 + dg312) + pow2(dg213)))) - 
     (3.*(dg113*dg211 + dg111*dg213) + 2.*dg211*dg311)*pow2(ginv13) - 
     (dg122*dg213 + dg113*dg222 + dg112*dg223 + 
        dg212*(dg123 + 2.*(dg213 + dg312)))*pow2(ginv23) - 
     (dg113*dg233 + dg213*(dg133 + 2.*dg313))*pow2(ginv33)) + 
  ginv22*(ddg2223*ginv33 + ginv23*
      (ddg2222 - ginv33*(2.*(dg223*dg322 + dg222*(dg233 + dg323)) + 
           6.*pow2(dg223))) - dg222*(6.*dg223 + 2.*dg322)*pow2(ginv23) - 
     2.*dg223*(dg233 + dg323)*pow2(ginv33)) + 
  ginv12*((ddg1223 + ddg2213)*ginv33 - 
     ginv22*((2.*dg122 + 6.*dg212)*dg222*ginv23 + 
        ((dg123 + 2.*dg213)*dg222 + (dg122 + 4.*dg212)*dg223)*ginv33) + 
     ginv23*(ddg1222 + ddg2212 - 
        ((dg122 + 2.*dg212)*dg233 + 
           dg223*(4.*dg123 + 8.*dg213 + 2.*dg312) + 
           dg222*(dg133 + 2.*dg313) + 2.*(dg213*dg322 + dg212*dg323))*
         ginv33) + ginv13*(ddg1212 + ddg2211 - 
        (4.*(dg112 + dg211)*dg223 + 
           dg212*(8.*dg213 + 4.*(dg123 + dg312)) + 
           2.*(dg122*dg213 + dg222*(dg113 + dg311) + dg211*dg322))*ginv23 \
- ginv22*(dg122*dg212 + (dg112 + 2.*dg211)*dg222 + 4.*pow2(dg212)) - 
        ginv33*((dg112 + 2.*dg211)*dg233 + dg212*(dg133 + 2.*dg313) + 
           2.*(dg223*dg311 + dg213*dg312 + dg211*dg323) + 
           4.*(dg123*dg213 + dg113*dg223 + pow2(dg213)))) - 
     (2.*(dg123*dg211 + dg112*dg213 + dg111*dg223 + 
           dg212*(dg113 + dg311)) + dg211*(4.*dg213 + 2.*dg312))*
      pow2(ginv13) - ((2.*dg122 + 4.*dg212)*dg223 + 
        dg222*(4.*dg213 + 2.*(dg123 + dg312)) + 2.*dg212*dg322)*
      pow2(ginv23) - ((dg123 + 2.*dg213)*dg233 + 
        dg223*(dg133 + 2.*dg313) + 2.*dg213*dg323)*pow2(ginv33)) + 
  ginv13*((ddg1233 + 2.*ddg2313)*ginv33 + 
     ginv22*(ddg2212 - ((dg122 + 8.*dg212)*dg223 + 
           dg222*(dg123 + 2.*(dg213 + dg312)) + 2.*dg212*dg322)*ginv23 - 
        (dg223*(4.*dg213 + 2.*(dg123 + dg312)) + 2.*dg212*(dg233 + dg323))*
         ginv33) + ginv23*(ddg1223 + ddg2213 + 2.*ddg2312 - 
        (3.*(dg133*dg223 + dg123*dg233) + dg233*(6.*dg213 + 4.*dg312) + 
           6.*(dg223*dg313 + dg213*dg323) + 4.*dg212*dg333)*ginv33) - 
     2.*dg212*dg222*pow2(ginv22) - 
     ((dg122 + 4.*dg212)*dg233 + dg223*(2.*dg123 + 4.*(dg213 + dg312)) + 
        dg222*(dg133 + 2.*dg313) + 2.*dg213*dg322 + 4.*dg212*dg323)*
      pow2(ginv23) - (dg233*(2.*dg133 + 4.*dg313) + 4.*dg213*dg333)*
      pow2(ginv33)) + ginv23*((ddg2233 + 2.*ddg2323)*ginv33 - 
     (4.*(dg233*dg323 + dg223*dg333) + 2.*pow2(dg233))*pow2(ginv33))
;

dGfromgdu31
=
-((dg222*dg311 + dg211*dg322 + 2.*((dg122 + dg212)*dg312 + dg112*dg322))*
     Power(ginv12,3)) - 2.*(dg111*dg311*Power(ginv11,3) + 
     Power(ginv13,3)*(dg133*dg313 + (dg113 + dg311)*dg333 + pow2(dg313))) + 
  (ddg1311 - ((4.*dg112 + 2.*dg211)*dg311 + 4.*dg111*dg312)*ginv12 - 
     (dg212*dg311 + (2.*dg112 + dg211)*dg312)*ginv22 - 
     (dg311*(dg213 + 2.*dg312) + dg211*dg313 + 
        2.*(dg113*dg312 + dg112*dg313))*ginv23 - 
     2.*(dg113 + dg311)*dg313*ginv33 - 
     ginv13*(4.*(dg113*dg311 + dg111*dg313) + 2.*pow2(dg311)))*pow2(ginv11) \
+ (ddg1322 + ddg2312 - (2.*dg122*dg322 + 3.*(dg222*dg312 + dg212*dg322))*
      ginv22 - ((2.*dg213 + 4.*dg312)*dg322 + 
        2.*(dg223*dg312 + dg222*dg313 + dg123*dg322 + 
           (dg122 + dg212)*dg323))*ginv23 - 
     (dg313*(dg223 + 2.*dg322) + (dg213 + 2.*(dg123 + dg312))*dg323)*
      ginv33 - ginv13*(4.*(dg123*dg312 + dg112*dg323) + 
        2.*(dg213*dg312 + (dg122 + dg212)*dg313 + dg113*dg322 + 
           dg311*(dg223 + dg322) + dg211*dg323 + pow2(dg312))))*pow2(ginv12) \
+ (ddg1333 + ddg3313 - (dg233*dg312 + dg223*dg313 + 
        (dg213 + 2.*(dg123 + dg312))*dg323 + dg212*dg333)*ginv22 - 
     (2.*(dg233*dg313 + dg133*dg323 + (dg123 + dg213)*dg333) + 
        4.*(dg313*dg323 + dg312*dg333))*ginv23 - 
     (2.*dg133 + 6.*dg313)*dg333*ginv33)*pow2(ginv13) + 
  ginv11*(ddg3313*ginv33 + ginv22*
      (ddg2312 - (dg222*dg313 + dg213*dg322 + 
           2.*(dg312*(dg223 + dg322) + dg212*dg323))*ginv23 - 
        (dg223*dg313 + (dg213 + 2.*dg312)*dg323)*ginv33) + 
     ginv23*(ddg2313 + ddg3312 - 
        (dg313*(dg233 + 4.*dg323) + (dg213 + 2.*dg312)*dg333)*ginv33) + 
     ginv12*(2.*ddg1312 + ddg2311 - 
        (dg311*(4.*dg123 + 3.*dg213 + 6.*dg312) + 3.*dg211*dg313 + 
           6.*(dg113*dg312 + dg112*dg313) + 4.*dg111*dg323)*ginv13 - 
        (dg222*dg311 + (2.*dg122 + 6.*dg212)*dg312 + 
           (2.*dg112 + dg211)*dg322)*ginv22 - 
        (4.*dg312*dg313 + 2.*((dg123 + dg213)*dg313 + 
              (dg113 + dg311)*dg323))*ginv33 - 
        ginv23*((2.*dg123 + 4.*dg213)*dg312 + dg311*(dg223 + 2.*dg322) + 
           dg211*dg323 + 2.*(dg122*dg313 + dg113*dg322 + dg112*dg323) + 
           4.*(dg212*dg313 + pow2(dg312)))) + 
     ginv13*(2.*ddg1313 + ddg3311 - 
        ((4.*dg213 + 8.*dg312)*dg313 + dg311*(dg233 + 2.*dg323) + 
           dg211*dg333 + 2.*(dg133*dg312 + dg123*dg313 + dg113*dg323 + 
              dg112*dg333))*ginv23 - 
        ginv22*(dg223*dg311 + dg211*dg323 + 
           2.*((dg123 + dg213)*dg312 + dg212*dg313 + dg112*dg323 + 
              pow2(dg312))) - 
        ginv33*(2.*(dg133*dg313 + (dg113 + dg311)*dg333) + 6.*pow2(dg313))) \
- ((2.*dg122 + 3.*dg212)*dg311 + (6.*dg112 + 3.*dg211)*dg312 + 
        2.*dg111*dg322)*pow2(ginv12) - 
     (6.*dg113*dg313 + dg311*(2.*dg133 + 6.*dg313) + 2.*dg111*dg333)*
      pow2(ginv13) - (dg222*dg312 + dg212*dg322)*pow2(ginv22) - 
     (dg313*(dg223 + 2.*dg322) + dg213*dg323 + dg312*(dg233 + 2.*dg323) + 
        dg212*dg333)*pow2(ginv23) - 2.*dg313*dg333*pow2(ginv33)) + 
  ginv12*(ddg3323*ginv33 + ginv13*
      (2.*ddg1323 + ddg2313 + ddg3312 - 
        (dg222*dg313 + (2.*dg123 + dg213)*dg322 + 
           dg312*(4.*dg223 + 2.*dg322) + (2.*dg122 + 4.*dg212)*dg323)*
         ginv22 - ((4.*dg213 + 8.*dg312)*dg323 + 
           4.*(dg313*(dg223 + dg322) + dg123*dg323) + 
           2.*(dg233*dg312 + dg133*dg322 + (dg122 + dg212)*dg333))*ginv23 \
- (dg313*(dg233 + 8.*dg323) + (dg213 + 2.*dg312)*dg333 + 
           2.*(dg133*dg323 + dg123*dg333))*ginv33) + 
     ginv22*(ddg2322 - 2.*(dg223 + dg322)*dg323*ginv33 - 
        ginv23*(3.*(dg223*dg322 + dg222*dg323) + 2.*pow2(dg322))) + 
     ginv23*(ddg2323 + ddg3322 - 
        ginv33*(dg233*dg323 + (dg223 + 2.*dg322)*dg333 + 4.*pow2(dg323))) - 
     (dg311*(dg233 + 4.*dg323) + 
        4.*((dg123 + dg312)*dg313 + dg113*dg323) + dg211*dg333 + 
        2.*(dg133*dg312 + dg213*dg313 + dg112*dg333))*pow2(ginv13) - 
     (2.*dg223*dg323 + dg322*(dg233 + 4.*dg323) + dg222*dg333)*
      pow2(ginv23) - 2.*(dg222*dg322*pow2(ginv22) + 
        dg323*dg333*pow2(ginv33))) + 
  ginv13*(ddg3333*ginv33 + ginv23*
      (ddg2333 + ddg3323 - (2.*dg233 + 6.*dg323)*dg333*ginv33) + 
     ginv22*(ddg2323 - (4.*dg223*dg323 + dg322*(dg233 + 2.*dg323) + 
           dg222*dg333)*ginv23 - 
        ginv33*(dg233*dg323 + dg223*dg333 + 2.*pow2(dg323))) - 
     (dg223*dg322 + dg222*dg323)*pow2(ginv22) - 
     2.*((dg233*dg323 + (dg223 + dg322)*dg333 + pow2(dg323))*pow2(ginv23) + 
        pow2(dg333)*pow2(ginv33)))
;

dGfromgdu32
=
-(((dg122 + 2.*dg212)*dg311 + 2.*(dg112 + dg211)*dg312 + dg111*dg322)*
     Power(ginv12,3)) - 2.*(dg222*dg322*Power(ginv22,3) + 
     Power(ginv23,3)*(dg233*dg323 + (dg223 + dg322)*dg333 + pow2(dg323))) - 
  (2.*dg111*dg311*ginv12 + (dg112*dg311 + dg111*dg312)*ginv22 + 
     (dg113*dg311 + dg111*dg313)*ginv23)*pow2(ginv11) + 
  (ddg1312 + ddg2311 - (4.*dg311*dg312 + 
        2.*((dg123 + dg213)*dg311 + dg113*dg312 + 
           (dg112 + dg211)*dg313 + dg111*dg323))*ginv13 - 
     ((3.*dg122 + 6.*dg212)*dg312 + 3.*dg112*dg322 + 
        2.*(dg222*dg311 + dg211*dg322))*ginv22 - 
     ((dg123 + 2.*(dg213 + dg312))*dg313 + (dg113 + 2.*dg311)*dg323)*
      ginv33 - ginv23*(4.*(dg213*dg312 + dg212*dg313) + 
        2.*(dg123*dg312 + dg122*dg313 + dg113*dg322 + 
           dg311*(dg223 + dg322) + (dg112 + dg211)*dg323 + pow2(dg312))))*
   pow2(ginv12) - ((dg123*dg313 + dg312*(dg133 + 2.*dg313) + 
        (dg113 + 2.*dg311)*dg323 + dg112*dg333)*ginv22 + 
     2.*ginv23*(dg133*dg313 + (dg113 + dg311)*dg333 + pow2(dg313)))*
   pow2(ginv13) + (ddg2322 - 2.*(dg223 + dg322)*dg323*ginv33 - 
     ginv23*(4.*(dg223*dg322 + dg222*dg323) + 2.*pow2(dg322)))*pow2(ginv22) \
+ (ddg2333 + ddg3323 - (2.*dg233 + 6.*dg323)*dg333*ginv33)*pow2(ginv23) + 
  ginv11*(-(ginv13*((dg311*(dg123 + 2.*dg312) + 
             2.*(dg113*dg312 + dg112*dg313) + dg111*dg323)*ginv22 + 
          (4.*dg113*dg313 + dg311*(dg133 + 2.*dg313) + dg111*dg333)*ginv23)\
) + ginv12*(ddg1311 - ((dg122 + 2.*dg212)*dg311 + 
           (6.*dg112 + 2.*dg211)*dg312 + dg111*dg322)*ginv22 - 
        (dg311*(dg123 + 2.*(dg213 + dg312)) + 2.*dg211*dg313 + 
           4.*(dg113*dg312 + dg112*dg313) + dg111*dg323)*ginv23 - 
        2.*(dg113 + dg311)*dg313*ginv33 - 
        ginv13*(3.*(dg113*dg311 + dg111*dg313) + 2.*pow2(dg311))) + 
     ginv22*(ddg1312 - ((dg123 + 2.*dg312)*dg313 + dg113*dg323)*ginv33 - 
        ginv23*(dg122*dg313 + dg113*dg322 + 
           2.*((dg123 + dg213)*dg312 + dg212*dg313 + dg112*dg323 + 
              pow2(dg312)))) + 
     ginv23*(ddg1313 - ginv33*
         (dg133*dg313 + dg113*dg333 + 2.*pow2(dg313))) - 
     ((3.*dg112 + 2.*dg211)*dg311 + 3.*dg111*dg312)*pow2(ginv12) - 
     ((dg122 + 2.*dg212)*dg312 + dg112*dg322)*pow2(ginv22) - 
     (dg133*dg312 + (dg123 + 2.*(dg213 + dg312))*dg313 + dg113*dg323 + 
        dg112*dg333)*pow2(ginv23)) + 
  ginv13*(ginv23*(ddg1333 + ddg3313 - (2.*dg133 + 6.*dg313)*dg333*ginv33) + 
     ginv22*(ddg1323 + ddg3312 - 
        (dg133*dg322 + (4.*dg123 + 2.*dg213 + 8.*dg312)*dg323 + 
           dg122*dg333 + 2.*(dg233*dg312 + dg313*(dg223 + dg322) + 
              dg212*dg333))*ginv23 - 
        ((dg133 + 4.*dg313)*dg323 + (dg123 + 2.*dg312)*dg333)*ginv33) - 
     (dg123*dg322 + dg122*dg323 + 
        2.*(dg312*(dg223 + dg322) + dg212*dg323))*pow2(ginv22) - 
     (2.*(dg233*dg313 + dg133*dg323 + (dg123 + dg213)*dg333) + 
        4.*(dg313*dg323 + dg312*dg333))*pow2(ginv23)) + 
  ginv12*(ddg3313*ginv33 + ginv22*
      (ddg1322 + 2.*ddg2312 - 
        (4.*(dg222*dg313 + dg213*dg322) + 
           3.*(dg123*dg322 + dg122*dg323) + 
           6.*(dg312*(dg223 + dg322) + dg212*dg323))*ginv23 - 
        ((2.*dg213 + 4.*dg312)*dg323 + 
           2.*(dg313*(dg223 + dg322) + dg123*dg323))*ginv33) + 
     ginv23*(ddg1323 + 2.*ddg2313 + ddg3312 - 
        (dg133*dg323 + dg313*(2.*dg233 + 8.*dg323) + 
           (dg123 + 2.*(dg213 + dg312))*dg333)*ginv33) + 
     ginv13*(ddg1313 + ddg3311 - 
        (8.*dg312*dg313 + 4.*
            ((dg123 + dg213)*dg313 + (dg113 + dg311)*dg323) + 
           2.*(dg233*dg311 + dg133*dg312 + (dg112 + dg211)*dg333))*ginv23 \
- ginv22*(dg122*dg313 + dg113*dg322 + 
           2.*(dg213*dg312 + dg212*dg313 + dg311*(dg223 + dg322) + 
              dg211*dg323) + 4.*(dg123*dg312 + dg112*dg323 + pow2(dg312))) \
- ginv33*(dg133*dg313 + (dg113 + 2.*dg311)*dg333 + 4.*pow2(dg313))) - 
     (2.*dg113*dg313 + dg311*(dg133 + 4.*dg313) + dg111*dg333)*
      pow2(ginv13) - (2.*dg122*dg322 + 4.*(dg222*dg312 + dg212*dg322))*
      pow2(ginv22) - (dg133*dg322 + 
        4.*(dg313*(dg223 + dg322) + (dg213 + dg312)*dg323) + 
        dg122*dg333 + 2.*(dg233*dg312 + dg123*dg323 + dg212*dg333))*
      pow2(ginv23) - 2.*dg313*dg333*pow2(ginv33)) + 
  ginv22*(ddg3323*ginv33 + ginv23*
      (2.*ddg2323 + ddg3322 - 
        ginv33*(2.*(dg233*dg323 + (dg223 + dg322)*dg333) + 6.*pow2(dg323))) \
- (6.*dg223*dg323 + dg322*(2.*dg233 + 6.*dg323) + 2.*dg222*dg333)*
      pow2(ginv23) - 2.*dg323*dg333*pow2(ginv33)) + 
  ginv23*(ddg3333*ginv33 - 2.*pow2(dg333)*pow2(ginv33))
;

dGfromgdu33
=
-((2.*dg113*dg313 + dg311*(dg133 + 4.*dg313) + dg111*dg333)*
     Power(ginv13,3)) - (2.*dg223*dg323 + dg322*(dg233 + 4.*dg323) + 
     dg222*dg333)*Power(ginv23,3) - 2.*Power(ginv33,3)*pow2(dg333) - 
  (2.*dg111*dg311*ginv13 + (dg112*dg311 + dg111*dg312)*ginv23 + 
     (dg113*dg311 + dg111*dg313)*ginv33)*pow2(ginv11) - 
  (((dg122 + 2.*dg212)*dg311 + 2.*(dg112 + dg211)*dg312 + dg111*dg322)*
      ginv13 + (dg222*dg311 + dg211*dg322 + 
        2.*((dg122 + dg212)*dg312 + dg112*dg322))*ginv23 + 
     (dg223*dg311 + (dg123 + dg213)*dg312 + (dg122 + dg212)*dg313 + 
        dg113*dg322 + (dg112 + dg211)*dg323)*ginv33)*pow2(ginv12) + 
  (ddg1313 + ddg3311 - ((2.*dg213 + 8.*dg312)*dg313 + 
        dg311*(dg233 + 4.*dg323) + dg211*dg333 + 
        2.*(dg133*dg312 + dg123*dg313 + dg113*dg323 + dg112*dg333))*ginv23 \
- ginv22*(dg223*dg311 + (dg123 + dg213)*dg312 + dg212*dg313 + 
        (dg112 + dg211)*dg323 + 2.*pow2(dg312)) - 
     ginv33*(4.*dg311*dg333 + 3.*(dg133*dg313 + dg113*dg333) + 
        6.*pow2(dg313)))*pow2(ginv13) - 
  (2.*dg222*dg322*ginv23 + (dg223*dg322 + dg222*dg323)*ginv33)*
   pow2(ginv22) + (ddg2323 + ddg3322 - 
     ginv33*(4.*dg322*dg333 + 3.*(dg233*dg323 + dg223*dg333) + 
        6.*pow2(dg323)))*pow2(ginv23) + ddg3333*pow2(ginv33) + 
  ginv13*((ddg1333 + 2.*ddg3313)*ginv33 + 
     ginv22*(ddg2312 - (dg222*dg313 + (dg123 + dg213)*dg322 + 
           dg122*dg323 + 4.*(dg312*(dg223 + dg322) + dg212*dg323))*ginv23 \
- (dg312*(dg233 + 4.*dg323) + 2.*(dg223*dg313 + (dg123 + dg213)*dg323) + 
           dg212*dg333)*ginv33) + 
     ginv23*(ddg1323 + ddg2313 + 2.*ddg3312 - 
        (12.*dg313*dg323 + (3.*dg213 + 8.*dg312)*dg333 + 
           3.*(dg233*dg313 + dg133*dg323 + dg123*dg333))*ginv33) - 
     (dg222*dg312 + dg212*dg322)*pow2(ginv22) - 
     ((dg133 + 4.*dg313)*dg322 + (2.*dg213 + 8.*dg312)*dg323 + 
        dg122*dg333 + 2.*(dg233*dg312 + dg223*dg313 + dg123*dg323 + 
           dg212*dg333))*pow2(ginv23) - 
     (2.*dg133 + 8.*dg313)*dg333*pow2(ginv33)) + 
  ginv23*((ddg2333 + 2.*ddg3323)*ginv33 - 
     (2.*dg233 + 8.*dg323)*dg333*pow2(ginv33)) + 
  ginv12*((ddg1323 + ddg2313)*ginv33 - 
     ginv22*((2.*dg122*dg322 + 3.*(dg222*dg312 + dg212*dg322))*ginv23 + 
        (dg222*dg313 + (dg123 + dg213)*dg322 + dg122*dg323 + 
           2.*(dg223*dg312 + dg212*dg323))*ginv33) + 
     ginv23*(ddg1322 + ddg2312 - 
        (dg233*dg312 + dg133*dg322 + 
           4.*(dg313*(dg223 + dg322) + (dg123 + dg213 + dg312)*dg323) + 
           (dg122 + dg212)*dg333)*ginv33) + 
     ginv13*(ddg1312 + ddg2311 - 
        (dg222*dg311 + (dg122 + 4.*dg212)*dg312 + (dg112 + dg211)*dg322)*
         ginv22 - (dg133*dg312 + dg311*(dg233 + 4.*dg323) + 
           4.*((dg123 + dg213 + dg312)*dg313 + dg113*dg323) + 
           (dg112 + dg211)*dg333)*ginv33 - 
        ginv23*(2.*(dg223*dg311 + dg122*dg313 + dg113*dg322 + 
              dg211*dg323) + 4.*
            ((dg123 + dg213)*dg312 + dg212*dg313 + dg311*dg322 + 
              dg112*dg323 + pow2(dg312)))) - 
     (4.*dg311*dg312 + 2.*((dg123 + dg213)*dg311 + dg113*dg312 + 
           (dg112 + dg211)*dg313 + dg111*dg323))*pow2(ginv13) - 
     ((2.*dg213 + 4.*dg312)*dg322 + 
        2.*(dg223*dg312 + dg222*dg313 + dg123*dg322 + 
           (dg122 + dg212)*dg323))*pow2(ginv23) - 
     (dg133*dg323 + dg313*(dg233 + 4.*dg323) + (dg123 + dg213)*dg333)*
      pow2(ginv33)) + ginv11*(ddg1313*ginv33 - 
     ginv12*(((3.*dg112 + 2.*dg211)*dg311 + 3.*dg111*dg312)*ginv13 + 
        ((dg122 + dg212)*dg311 + (4.*dg112 + dg211)*dg312 + dg111*dg322)*
         ginv23 + ((dg123 + dg213)*dg311 + dg211*dg313 + 
           2.*(dg113*dg312 + dg112*dg313) + dg111*dg323)*ginv33) - 
     ginv22*(((dg122 + 2.*dg212)*dg312 + dg112*dg322)*ginv23 + 
        ((dg123 + dg213)*dg312 + dg212*dg313 + dg112*dg323)*ginv33) + 
     ginv13*(ddg1311 - (dg212*dg311 + (2.*dg112 + dg211)*dg312)*ginv22 - 
        ((dg123 + dg213)*dg311 + 4.*(dg113 + dg311)*dg312 + 
           (4.*dg112 + dg211)*dg313 + dg111*dg323)*ginv23 - 
        (6.*dg113*dg313 + dg311*(dg133 + 4.*dg313) + dg111*dg333)*ginv33) + 
     ginv23*(ddg1312 - (dg312*(dg133 + 4.*dg313) + 
           2.*((dg123 + dg213)*dg313 + dg113*dg323) + dg112*dg333)*ginv33) \
- (3.*(dg113*dg311 + dg111*dg313) + 2.*pow2(dg311))*pow2(ginv13) - 
     ((dg123 + dg213)*dg312 + (dg122 + dg212)*dg313 + dg113*dg322 + 
        dg112*dg323 + 2.*pow2(dg312))*pow2(ginv23) - 
     (dg133*dg313 + dg113*dg333 + 2.*pow2(dg313))*pow2(ginv33)) + 
  ginv22*(ddg2323*ginv33 + ginv23*
      (ddg2322 - (6.*dg223*dg323 + dg322*(dg233 + 4.*dg323) + dg222*dg333)*
         ginv33) - (3.*(dg223*dg322 + dg222*dg323) + 2.*pow2(dg322))*
      pow2(ginv23) - (dg233*dg323 + dg223*dg333 + 2.*pow2(dg323))*
      pow2(ginv33))
;

R11
=
dG11*g11 + dG12*g12 + dG13*g13 + gammado111*Gfromg1 + gammado112*Gfromg2 + 
  gammado113*Gfromg3 + (-0.5*ddg1111 + 3.*gamma111*gammado111 + 
     2.*(gamma211*gammado112 + gamma311*gammado113) + 
     gamma211*gammado211 + gamma311*gammado311)*ginv11 + 
  (-ddg1211 + 3.*(gamma112*gammado111 + gamma111*gammado112) + 
     2.*(gamma212*gammado112 + gamma312*gammado113 + 
        gamma211*gammado122 + gamma311*gammado123) + gamma212*gammado211 + 
     gamma211*gammado212 + gamma312*gammado311 + gamma311*gammado312)*ginv12 \
+ (-ddg1311 + 3.*(gamma113*gammado111 + gamma111*gammado113) + 
     2.*(gamma213*gammado112 + gamma313*gammado113 + 
        gamma211*gammado123 + gamma311*gammado133) + gamma213*gammado211 + 
     gamma211*gammado213 + gamma313*gammado311 + gamma311*gammado313)*ginv13 \
+ (-0.5*ddg2211 + 3.*gamma112*gammado112 + 
     2.*(gamma212*gammado122 + gamma312*gammado123) + 
     gamma212*gammado212 + gamma312*gammado312)*ginv22 + 
  (-ddg2311 + 3.*(gamma113*gammado112 + gamma112*gammado113) + 
     2.*(gamma213*gammado122 + (gamma212 + gamma313)*gammado123 + 
        gamma312*gammado133) + gamma213*gammado212 + gamma212*gammado213 + 
     gamma313*gammado312 + gamma312*gammado313)*ginv23 + 
  (-0.5*ddg3311 + 3.*gamma113*gammado113 + 
     2.*(gamma213*gammado123 + gamma313*gammado133) + gamma213*gammado213 + 
     gamma313*gammado313)*ginv33
;

R12
=
0.5*(dG21*g11 + (dG11 + dG22)*g12 + dG23*g13 + dG12*g22 + dG13*g23 + 
     (gammado112 + gammado211)*Gfromg1 + 
     (gammado122 + gammado212)*Gfromg2 + (gammado123 + gammado213)*Gfromg3) \
+ (-0.5*ddg1112 + gamma112*gammado111 + (gamma111 + gamma212)*gammado112 + 
     gamma312*gammado113 + gamma111*gammado211 + 2.*gamma211*gammado212 + 
     gamma311*(gammado213 + gammado312))*ginv11 + 
  (-ddg1212 + gamma122*gammado111 + (2.*gamma112 + gamma222)*gammado112 + 
     gamma322*gammado113 + (gamma111 + gamma212)*gammado122 + 
     gamma112*gammado211 + (gamma111 + 2.*gamma212)*gammado212 + 
     2.*gamma211*gammado222 + 
     gamma312*(gammado123 + gammado213 + gammado312) + 
     gamma311*(gammado223 + gammado322))*ginv12 + 
  (-ddg1312 + gamma123*gammado111 + (gamma113 + gamma223)*gammado112 + 
     (gamma112 + gamma323)*gammado113 + (gamma111 + gamma212)*gammado123 + 
     gamma312*gammado133 + gamma113*gammado211 + 
     (gamma111 + gamma313)*gammado213 + 
     2.*(gamma213*gammado212 + gamma211*gammado223) + 
     gamma313*gammado312 + gamma311*(gammado233 + gammado323))*ginv13 + 
  (-0.5*ddg2212 + gamma122*gammado112 + (gamma112 + gamma222)*gammado122 + 
     gamma322*gammado123 + gamma112*gammado212 + 2.*gamma212*gammado222 + 
     gamma312*(gammado223 + gammado322))*ginv22 + 
  (-ddg2312 + gamma123*gammado112 + gamma122*gammado113 + 
     (gamma113 + gamma223)*gammado122 + 
     (gamma112 + gamma222 + gamma323)*gammado123 + gamma322*gammado133 + 
     gamma113*gammado212 + gamma112*gammado213 + 
     2.*(gamma213*gammado222 + gamma212*gammado223) + 
     gamma313*(gammado223 + gammado322) + 
     gamma312*(gammado233 + gammado323))*ginv23 + 
  (-0.5*ddg3312 + gamma123*gammado113 + (gamma113 + gamma223)*gammado123 + 
     gamma323*gammado133 + gamma113*gammado213 + 2.*gamma213*gammado223 + 
     gamma313*(gammado233 + gammado323))*ginv33
;

R13
=
0.5*(dG31*g11 + dG32*g12 + (dG11 + dG33)*g13 + dG12*g23 + dG13*g33 + 
     (gammado113 + gammado311)*Gfromg1 + 
     (gammado123 + gammado312)*Gfromg2 + (gammado133 + gammado313)*Gfromg3) \
+ (-0.5*ddg1113 + gamma113*gammado111 + gamma213*gammado112 + 
     (gamma111 + gamma313)*gammado113 + gamma111*gammado311 + 
     gamma211*(gammado213 + gammado312) + 2.*gamma311*gammado313)*ginv11 + 
  (-ddg1213 + gamma123*gammado111 + (gamma113 + gamma223)*gammado112 + 
     (gamma112 + gamma323)*gammado113 + gamma213*gammado122 + 
     (gamma111 + gamma313)*gammado123 + gamma112*gammado311 + 
     gamma111*gammado312 + gamma212*(gammado213 + gammado312) + 
     gamma211*(gammado223 + gammado322) + 
     2.*(gamma312*gammado313 + gamma311*gammado323))*ginv12 + 
  (-ddg1313 + gamma133*gammado111 + gamma233*gammado112 + 
     (2.*gamma113 + gamma333)*gammado113 + 
     (gamma111 + gamma313)*gammado133 + gamma113*gammado311 + 
     gamma213*(gammado123 + gammado213 + gammado312) + 
     (gamma111 + 2.*gamma313)*gammado313 + 
     gamma211*(gammado233 + gammado323) + 2.*gamma311*gammado333)*ginv13 + 
  (-0.5*ddg2213 + gamma123*gammado112 + gamma223*gammado122 + 
     (gamma112 + gamma323)*gammado123 + gamma112*gammado312 + 
     gamma212*(gammado223 + gammado322) + 2.*gamma312*gammado323)*ginv22 + 
  (-ddg2313 + gamma133*gammado112 + gamma123*gammado113 + 
     gamma233*gammado122 + (gamma113 + gamma223 + gamma333)*gammado123 + 
     (gamma112 + gamma323)*gammado133 + gamma113*gammado312 + 
     gamma112*gammado313 + gamma213*(gammado223 + gammado322) + 
     gamma212*(gammado233 + gammado323) + 
     2.*(gamma313*gammado323 + gamma312*gammado333))*ginv23 + 
  (-0.5*ddg3313 + gamma133*gammado113 + gamma233*gammado123 + 
     (gamma113 + gamma333)*gammado133 + gamma113*gammado313 + 
     gamma213*(gammado233 + gammado323) + 2.*gamma313*gammado333)*ginv33
;

R22
=
dG21*g12 + dG22*g22 + dG23*g23 + gammado212*Gfromg1 + gammado222*Gfromg2 + 
  gammado223*Gfromg3 + (-0.5*ddg1122 + 
     gamma112*(gammado112 + 2.*gammado211) + 3.*gamma212*gammado212 + 
     gamma312*(2.*gammado213 + gammado312))*ginv11 + 
  (-ddg1222 + gamma122*(gammado112 + 2.*gammado211) + 
     gamma112*(gammado122 + 2.*gammado212) + 
     3.*(gamma222*gammado212 + gamma212*gammado222) + 
     2.*(gamma322*gammado213 + gamma312*gammado223) + 
     gamma322*gammado312 + gamma312*gammado322)*ginv12 + 
  (-ddg1322 + gamma123*(gammado112 + 2.*gammado211) + 
     gamma112*(gammado123 + 2.*gammado213) + 
     3.*(gamma223*gammado212 + gamma212*gammado223) + 
     2.*(gamma323*gammado213 + gamma312*gammado233) + 
     gamma323*gammado312 + gamma312*gammado323)*ginv13 + 
  (-0.5*ddg2222 + gamma122*(gammado122 + 2.*gammado212) + 
     3.*gamma222*gammado222 + gamma322*(2.*gammado223 + gammado322))*ginv22 \
+ (-ddg2322 + gamma123*(gammado122 + 2.*gammado212) + 
     gamma122*(gammado123 + 2.*gammado213) + 
     3.*(gamma223*gammado222 + gamma222*gammado223) + 
     2.*(gamma323*gammado223 + gamma322*gammado233) + 
     gamma323*gammado322 + gamma322*gammado323)*ginv23 + 
  (-0.5*ddg3322 + gamma123*(gammado123 + 2.*gammado213) + 
     3.*gamma223*gammado223 + gamma323*(2.*gammado233 + gammado323))*ginv33
;

R23
=
0.5*(dG31*g12 + dG21*g13 + dG32*g22 + (dG22 + dG33)*g23 + dG23*g33 + 
     (gammado213 + gammado312)*Gfromg1 + 
     (gammado223 + gammado322)*Gfromg2 + (gammado233 + gammado323)*Gfromg3) \
+ (-0.5*ddg1123 + gamma113*gammado211 + gamma213*gammado212 + 
     (gamma212 + gamma313)*gammado213 + 
     gamma112*(gammado113 + gammado311) + gamma212*gammado312 + 
     2.*gamma312*gammado313)*ginv11 + 
  (-ddg1223 + gamma123*gammado211 + (gamma113 + gamma223)*gammado212 + 
     (gamma222 + gamma323)*gammado213 + gamma213*gammado222 + 
     (gamma212 + gamma313)*gammado223 + 
     gamma122*(gammado113 + gammado311) + gamma222*gammado312 + 
     gamma112*(gammado123 + gammado312) + gamma212*gammado322 + 
     2.*(gamma322*gammado313 + gamma312*gammado323))*ginv12 + 
  (-ddg1323 + gamma133*gammado211 + gamma233*gammado212 + 
     (gamma113 + gamma223 + gamma333)*gammado213 + gamma213*gammado223 + 
     (gamma212 + gamma313)*gammado233 + 
     gamma123*(gammado113 + gammado311) + gamma223*gammado312 + 
     gamma112*(gammado133 + gammado313) + gamma212*gammado323 + 
     2.*(gamma323*gammado313 + gamma312*gammado333))*ginv13 + 
  (-0.5*ddg2223 + gamma123*gammado212 + gamma223*gammado222 + 
     (gamma222 + gamma323)*gammado223 + 
     gamma122*(gammado123 + gammado312) + gamma222*gammado322 + 
     2.*gamma322*gammado323)*ginv22 + 
  (-ddg2323 + gamma133*gammado212 + gamma233*gammado222 + 
     (2.*gamma223 + gamma333)*gammado223 + 
     (gamma222 + gamma323)*gammado233 + 
     gamma123*(gammado123 + gammado213 + gammado312) + 
     gamma122*(gammado133 + gammado313) + gamma223*gammado322 + 
     (gamma222 + 2.*gamma323)*gammado323 + 2.*gamma322*gammado333)*ginv23 + 
  (-0.5*ddg3323 + gamma133*gammado213 + gamma233*gammado223 + 
     (gamma223 + gamma333)*gammado233 + 
     gamma123*(gammado133 + gammado313) + gamma223*gammado323 + 
     2.*gamma323*gammado333)*ginv33
;

R33
=
dG31*g13 + dG32*g23 + dG33*g33 + gammado313*Gfromg1 + gammado323*Gfromg2 + 
  gammado333*Gfromg3 + (-0.5*ddg1133 + 
     gamma113*(gammado113 + 2.*gammado311) + 
     gamma213*(gammado213 + 2.*gammado312) + 3.*gamma313*gammado313)*ginv11 \
+ (-ddg1233 + gamma123*(gammado113 + 2.*gammado311) + 
     gamma113*(gammado123 + 2.*gammado312) + 
     gamma223*(gammado213 + 2.*gammado312) + 
     gamma213*(gammado223 + 2.*gammado322) + 
     3.*(gamma323*gammado313 + gamma313*gammado323))*ginv12 + 
  (-ddg1333 + gamma133*(gammado113 + 2.*gammado311) + 
     gamma233*(gammado213 + 2.*gammado312) + 
     gamma113*(gammado133 + 2.*gammado313) + 
     gamma213*(gammado233 + 2.*gammado323) + 
     3.*(gamma333*gammado313 + gamma313*gammado333))*ginv13 + 
  (-0.5*ddg2233 + gamma123*(gammado123 + 2.*gammado312) + 
     gamma223*(gammado223 + 2.*gammado322) + 3.*gamma323*gammado323)*ginv22 \
+ (-ddg2333 + gamma133*(gammado123 + 2.*gammado312) + 
     gamma123*(gammado133 + 2.*gammado313) + 
     gamma233*(gammado223 + 2.*gammado322) + 
     gamma223*(gammado233 + 2.*gammado323) + 
     3.*(gamma333*gammado323 + gamma323*gammado333))*ginv23 + 
  (-0.5*ddg3333 + gamma133*(gammado133 + 2.*gammado313) + 
     gamma233*(gammado233 + 2.*gammado323) + 3.*gamma333*gammado333)*ginv33
;

ff
=
chi
;

oochipsipower
=
1/chipsipower
;

f
=
oochipsipower*log(ff)
;

psim4
=
exp(-4.*f)
;

df1
=
(dchi1*oochipsipower)/chi
;

df2
=
(dchi2*oochipsipower)/chi
;

df3
=
(dchi3*oochipsipower)/chi
;

ddf11
=
(ddchi11*oochipsipower)/chi - chipsipower*pow2(df1)
;

ddf12
=
-(chipsipower*df1*df2) + (ddchi12*oochipsipower)/chi
;

ddf13
=
-(chipsipower*df1*df3) + (ddchi13*oochipsipower)/chi
;

ddf22
=
(ddchi22*oochipsipower)/chi - chipsipower*pow2(df2)
;

ddf23
=
-(chipsipower*df2*df3) + (ddchi23*oochipsipower)/chi
;

ddf33
=
(ddchi33*oochipsipower)/chi - chipsipower*pow2(df3)
;

cddf11
=
ddf11 - df1*gamma111 - df2*gamma211 - df3*gamma311
;

cddf12
=
ddf12 - df1*gamma112 - df2*gamma212 - df3*gamma312
;

cddf13
=
ddf13 - df1*gamma113 - df2*gamma213 - df3*gamma313
;

cddf22
=
ddf22 - df1*gamma122 - df2*gamma222 - df3*gamma322
;

cddf23
=
ddf23 - df1*gamma123 - df2*gamma223 - df3*gamma323
;

cddf33
=
ddf33 - df1*gamma133 - df2*gamma233 - df3*gamma333
;

trcddf
=
cddf11*ginv11 + cddf22*ginv22 + 
  2.*(cddf12*ginv12 + cddf13*ginv13 + cddf23*ginv23) + cddf33*ginv33
;

Rphi11
=
-2.*(cddf11 + g11*trcddf) + (4. - 4.*g11*ginv11)*pow2(df1) - 
  g11*(8.*(df1*(df2*ginv12 + df3*ginv13) + df2*df3*ginv23) + 
     4.*(ginv22*pow2(df2) + ginv33*pow2(df3)))
;

Rphi12
=
df1*df2*(4. - 8.*g12*ginv12) - 2.*(cddf12 + g12*trcddf) - 
  g12*(8.*df3*(df1*ginv13 + df2*ginv23) + 
     4.*(ginv11*pow2(df1) + ginv22*pow2(df2) + ginv33*pow2(df3)))
;

Rphi13
=
df1*(4.*df3 - 8.*df2*g13*ginv12) - 2.*(cddf13 + g13*trcddf) - 
  g13*(8.*df3*(df1*ginv13 + df2*ginv23) + 
     4.*(ginv11*pow2(df1) + ginv22*pow2(df2) + ginv33*pow2(df3)))
;

Rphi22
=
-2.*(cddf22 + g22*trcddf) + (4. - 4.*g22*ginv22)*pow2(df2) - 
  g22*(8.*(df1*(df2*ginv12 + df3*ginv13) + df2*df3*ginv23) + 
     4.*(ginv11*pow2(df1) + ginv33*pow2(df3)))
;

Rphi23
=
df2*(-8.*df1*g23*ginv12 + df3*(4. - 8.*g23*ginv23)) - 
  2.*(cddf23 + g23*trcddf) - g23*
   (8.*df1*df3*ginv13 + 4.*(ginv11*pow2(df1) + ginv22*pow2(df2) + 
        ginv33*pow2(df3)))
;

Rphi33
=
-2.*(cddf33 + g33*trcddf) - g33*
   (8.*(df1*(df2*ginv12 + df3*ginv13) + df2*df3*ginv23) + 
     4.*(ginv11*pow2(df1) + ginv22*pow2(df2))) + 
  (4. - 4.*g33*ginv33)*pow2(df3)
;

Rf11
=
R11 + Rphi11
;

Rf12
=
R12 + Rphi12
;

Rf13
=
R13 + Rphi13
;

Rf22
=
R22 + Rphi22
;

Rf23
=
R23 + Rphi23
;

Rf33
=
R33 + Rphi33
;

Rhat
=
psim4*(ginv11*Rf11 + ginv22*Rf22 + 
    2.*(ginv12*Rf12 + ginv13*Rf13 + ginv23*Rf23) + ginv33*Rf33)
;

cdda11
=
dda11 - da2*gamma211 - da3*gamma311 + 
  da1*(-gamma111 + df1*(-4. + 2.*g11*ginv11)) + 
  2.*g11*((da2*df1 + da1*df2)*ginv12 + (da3*df1 + da1*df3)*ginv13 + 
     da2*df2*ginv22 + (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33)
;

cdda12
=
dda12 - da1*gamma112 - da2*gamma212 - da3*gamma312 + 
  2.*(-(da2*df1) - da1*df2 + g12*
      (da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
        (da3*df1 + da1*df3)*ginv13 + da2*df2*ginv22 + 
        (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33))
;

cdda13
=
dda13 - da1*gamma113 - da2*gamma213 - da3*gamma313 + 
  2.*(-(da3*df1) - da1*df3 + g13*
      (da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
        (da3*df1 + da1*df3)*ginv13 + da2*df2*ginv22 + 
        (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33))
;

cdda22
=
dda22 - da1*gamma122 - da2*(4.*df2 + gamma222) - da3*gamma322 + 
  2.*g22*(da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
     (da3*df1 + da1*df3)*ginv13 + da2*df2*ginv22 + 
     (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33)
;

cdda23
=
dda23 - da1*gamma123 - da2*gamma223 - da3*gamma323 + 
  2.*(-(da3*df2) - da2*df3 + g23*
      (da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
        (da3*df1 + da1*df3)*ginv13 + da2*df2*ginv22 + 
        (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33))
;

cdda33
=
dda33 - da1*gamma133 - da2*gamma233 - da3*(4.*df3 + gamma333) + 
  2.*g33*(da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
     (da3*df1 + da1*df3)*ginv13 + da2*df2*ginv22 + 
     (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33)
;

trcdda
=
(cdda11*ginv11 + cdda22*ginv22 + 
    2.*(cdda12*ginv12 + cdda13*ginv13 + cdda23*ginv23) + cdda33*ginv33)*psim4
;

AA11
=
2.*(A11*(A12*ginv12 + A13*ginv13) + A12*A13*ginv23) + ginv11*pow2(A11) + 
  ginv22*pow2(A12) + ginv33*pow2(A13)
;

AA12
=
(A12*A13 + A11*A23)*ginv13 + A12*(A11*ginv11 + A22*ginv22) + 
  (A13*A22 + A12*A23)*ginv23 + A13*A23*ginv33 + ginv12*(A11*A22 + pow2(A12))
;

AA13
=
(A12*A13 + A11*A23)*ginv12 + A12*A23*ginv22 + (A13*A23 + A12*A33)*ginv23 + 
  A13*(A11*ginv11 + A33*ginv33) + ginv13*(A11*A33 + pow2(A13))
;

AA21
=
(A12*A13 + A11*A23)*ginv13 + A12*(A11*ginv11 + A22*ginv22) + 
  (A13*A22 + A12*A23)*ginv23 + A13*A23*ginv33 + ginv12*(A11*A22 + pow2(A12))
;

AA22
=
2.*(A12*(A22*ginv12 + A23*ginv13) + A22*A23*ginv23) + ginv11*pow2(A12) + 
  ginv22*pow2(A22) + ginv33*pow2(A23)
;

AA23
=
A12*A13*ginv11 + (A13*A22 + A12*A23)*ginv12 + (A13*A23 + A12*A33)*ginv13 + 
  A23*(A22*ginv22 + A33*ginv33) + ginv23*(A22*A33 + pow2(A23))
;

AA31
=
(A12*A13 + A11*A23)*ginv12 + A12*A23*ginv22 + (A13*A23 + A12*A33)*ginv23 + 
  A13*(A11*ginv11 + A33*ginv33) + ginv13*(A11*A33 + pow2(A13))
;

AA32
=
A12*A13*ginv11 + (A13*A22 + A12*A23)*ginv12 + (A13*A23 + A12*A33)*ginv13 + 
  A23*(A22*ginv22 + A33*ginv33) + ginv23*(A22*A33 + pow2(A23))
;

AA33
=
2.*(A13*(A23*ginv12 + A33*ginv13) + A23*A33*ginv23) + ginv11*pow2(A13) + 
  ginv22*pow2(A23) + ginv33*pow2(A33)
;

cdA111
=
dA111 - 2.*(A11*gamma111 + A12*gamma211 + A13*gamma311)
;

cdA112
=
dA112 - A11*gamma112 - A22*gamma211 - A12*(gamma111 + gamma212) - 
  A23*gamma311 - A13*gamma312
;

cdA113
=
dA113 - A11*gamma113 - A23*gamma211 - A12*gamma213 - A33*gamma311 - 
  A13*(gamma111 + gamma313)
;

cdA122
=
dA122 - 2.*(A12*gamma112 + A22*gamma212 + A23*gamma312)
;

cdA123
=
dA123 - A13*gamma112 - A12*gamma113 - A22*gamma213 - A33*gamma312 - 
  A23*(gamma212 + gamma313)
;

cdA133
=
dA133 - 2.*(A13*gamma113 + A23*gamma213 + A33*gamma313)
;

cdA211
=
dA211 - 2.*(A11*gamma112 + A12*gamma212 + A13*gamma312)
;

cdA212
=
dA212 - A11*gamma122 - A22*gamma212 - A12*(gamma112 + gamma222) - 
  A23*gamma312 - A13*gamma322
;

cdA213
=
dA213 - A11*gamma123 - A23*gamma212 - A12*gamma223 - A33*gamma312 - 
  A13*(gamma112 + gamma323)
;

cdA222
=
dA222 - 2.*(A12*gamma122 + A22*gamma222 + A23*gamma322)
;

cdA223
=
dA223 - A13*gamma122 - A12*gamma123 - A22*gamma223 - A33*gamma322 - 
  A23*(gamma222 + gamma323)
;

cdA233
=
dA233 - 2.*(A13*gamma123 + A23*gamma223 + A33*gamma323)
;

cdA311
=
dA311 - 2.*(A11*gamma113 + A12*gamma213 + A13*gamma313)
;

cdA312
=
dA312 - A11*gamma123 - A22*gamma213 - A12*(gamma113 + gamma223) - 
  A23*gamma313 - A13*gamma323
;

cdA313
=
dA313 - A11*gamma133 - A23*gamma213 - A12*gamma233 - A33*gamma313 - 
  A13*(gamma113 + gamma333)
;

cdA322
=
dA322 - 2.*(A12*gamma123 + A22*gamma223 + A23*gamma323)
;

cdA323
=
dA323 - A13*gamma123 - A12*gamma133 - A22*gamma233 - A33*gamma323 - 
  A23*(gamma223 + gamma333)
;

cdA333
=
dA333 - 2.*(A13*gamma133 + A23*gamma233 + A33*gamma333)
;

divbeta
=
db11 + db22 + db33
;

totdivbeta
=
0.66666666666666666667*divbeta
;

lieA11
=
beta1*dA111 + beta2*dA211 + beta3*dA311 + 
  2.*(A11*db11 + A12*db12 + A13*db13) - A11*totdivbeta
;

lieA12
=
beta1*dA112 + beta2*dA212 + beta3*dA312 + A22*db12 + A23*db13 + A11*db21 + 
  A13*db23 + A12*(db11 + db22 - totdivbeta)
;

lieA13
=
beta1*dA113 + beta2*dA213 + beta3*dA313 + A23*db12 + A33*db13 + A11*db31 + 
  A12*db32 + A13*(db11 + db33 - totdivbeta)
;

lieA22
=
beta1*dA122 + beta2*dA222 + beta3*dA322 + 
  2.*(A12*db21 + A22*db22 + A23*db23) - A22*totdivbeta
;

lieA23
=
beta1*dA123 + beta2*dA223 + beta3*dA323 + A13*db21 + A33*db23 + A12*db31 + 
  A22*db32 + A23*(db22 + db33 - totdivbeta)
;

lieA33
=
beta1*dA133 + beta2*dA233 + beta3*dA333 + 
  2.*(A13*db31 + A23*db32 + A33*db33) - A33*totdivbeta
;

DTheta
=
dTheta1*sup1 + dTheta2*sup2 + dTheta3*sup3
;

rACss
=
2.*((A23*alpha*K + lieA23)*sup2*sup3 + 
     sup1*((A12*alpha*K + lieA12)*sup2 + A13*alpha*K*sup3) + 
     psim4*((-cdda23 + alpha*Rf23)*sup2*sup3 + 
        sup1*((-cdda12 + alpha*Rf12)*sup2 - cdda13*sup3))) + 
  0.66666666666666666667*(g13*sup1 + g23*sup2)*sup3*trcdda + 
  sup1*(2.*(-(AA31*alpha) + lieA13)*sup3 + 
     0.66666666666666666667*g12*sup2*trcdda) + 
  (lieA11 + psim4*(-cdda11 + alpha*Rf11) + 
     0.33333333333333333333*g11*(-(alpha*Rhat) + trcdda))*pow2(sup1) + 
  (lieA22 - cdda22*psim4 + alpha*
      (A22*K + psim4*Rf22 - 0.33333333333333333333*g22*Rhat) + 
     0.33333333333333333333*g22*trcdda)*pow2(sup2) + 
  (lieA33 - cdda33*psim4 + alpha*
      (A33*K + psim4*Rf33 - 0.33333333333333333333*g33*Rhat) + 
     0.33333333333333333333*g33*trcdda)*pow2(sup3) + 
  alpha*(ginv11*((-2.*cdA111*chi + 3.*A11*dchi1)*sup1 + 
        (-2.*cdA112*chi + 3.*A12*dchi1)*sup2 + 
        (-2.*cdA113*chi + 3.*A13*dchi1)*sup3) + 
     ginv22*((-2.*cdA212*chi + 3.*A12*dchi2)*sup1 + 
        (-2.*cdA222*chi + 3.*A22*dchi2)*sup2 + 
        (-2.*cdA223*chi + 3.*A23*dchi2)*sup3) + 
     ginv33*((-2.*cdA313*chi + 3.*A13*dchi3)*sup1 + 
        (-2.*cdA323*chi + 3.*A23*dchi3)*sup2 + 
        (-2.*cdA333*chi + 3.*A33*dchi3)*sup3) + 
     chi*(-2.*DTheta + 1.3333333333333333333*
         (dK1*sup1 + dK2*sup2 + dK3*sup3)) + 
     ginv12*((-2.*cdA212*chi + 3.*A12*dchi2)*sup2 + 
        (-2.*cdA213*chi + 3.*A13*dchi2)*sup3 - 
        2.*chi*((cdA112 + cdA211)*sup1 + cdA122*sup2 + cdA123*sup3) + 
        3.*((A12*dchi1 + A11*dchi2)*sup1 + dchi1*(A22*sup2 + A23*sup3))) + 
     ginv13*((-2.*cdA312*chi + 3.*A12*dchi3)*sup2 + 
        (-2.*cdA313*chi + 3.*A13*dchi3)*sup3 - 
        2.*chi*((cdA113 + cdA311)*sup1 + cdA123*sup2 + cdA133*sup3) + 
        3.*((A13*dchi1 + A11*dchi3)*sup1 + dchi1*(A23*sup2 + A33*sup3))) + 
     ginv23*((-2.*cdA322*chi + 3.*A22*dchi3)*sup2 + 
        (-2.*cdA323*chi + 3.*A23*dchi3)*sup3 - 
        2.*chi*((cdA213 + cdA312)*sup1 + cdA223*sup2 + cdA233*sup3) + 
        3.*((A13*dchi2 + A12*dchi3)*sup1 + dchi2*(A23*sup2 + A33*sup3))) + 
     (0.33333333333333333333*((dG11 - dGfromgdu11)*qud11 + 
           (dG12 - dGfromgdu12)*qud12 + (dG13 - dGfromgdu13)*qud13 + 
           (dG21 - dGfromgdu21)*qud21 + (dG22 - dGfromgdu22)*qud22 + 
           (dG23 - dGfromgdu23)*qud23 + (dG31 - dGfromgdu31)*qud31 + 
           (dG32 - dGfromgdu32)*qud32 + (dG33 - dGfromgdu33)*qud33) + 
        kappa1*((G1 - Gfromg1)*sdown1 + (G2 - Gfromg2)*sdown2 + 
           (G3 - Gfromg3)*sdown3) + 
        0.66666666666666666667*
         ((dGfromgdu21*sdown1 + dGfromgdu22*sdown2)*sup2 + 
           sdown3*((-dG13 + dGfromgdu13)*sup1 - dG23*sup2 - dG33*sup3) + 
           sdown1*((-dG11 + dGfromgdu11)*sup1 - dG21*sup2 - dG31*sup3 + 
              dGfromgdu31*sup3) + 
           sdown2*((-dG12 + dGfromgdu12)*sup1 - dG22*sup2 - dG32*sup3 + 
              dGfromgdu32*sup3)))*pow2(chi) + 
     0.66666666666666666667*sup2*
      (-(Rhat*(g12*sup1 + g23*sup3)) + dGfromgdu23*sdown3*pow2(chi)) + 
     sup3*((2.*psim4*Rf13 - 0.66666666666666666667*g13*Rhat)*sup1 + 
        0.66666666666666666667*dGfromgdu33*sdown3*pow2(chi)) + 
     (-2.*AA11 + A11*K)*pow2(sup1) - 
     2.*((AA23 + AA32)*sup2*sup3 + sup1*((AA12 + AA21)*sup2 + AA13*sup3) + 
        AA22*pow2(sup2) + AA33*pow2(sup3)))
;

rACsA1
=
(qud11*(lieA11 + alpha*chi*Rf11) + 
     qud21*(lieA12 + alpha*(-2.*AA12 + chi*Rf12)) + 
     qud31*(lieA13 + alpha*(-2.*AA13 + chi*Rf13)))*sup1 + 
  qud11*((-(cdda11*chi) + A11*alpha*K)*sup1 + lieA12*sup2 + 
     (A13*alpha*K + lieA13)*sup3 + 
     alpha*((-2.*AA21 + A12*K)*sup2 - 2.*AA31*sup3)) + 
  qud21*((-(cdda12*chi) + A12*alpha*K)*sup1 + lieA22*sup2 + 
     (A23*alpha*K + lieA23)*sup3 + 
     alpha*((-2.*AA22 + A22*K)*sup2 - 2.*AA32*sup3)) + 
  qud31*((-(cdda13*chi) + A13*alpha*K)*sup1 + lieA23*sup2 + 
     (A33*alpha*K + lieA33)*sup3 + 
     alpha*((-2.*AA23 + A23*K)*sup2 - 2.*AA33*sup3)) + 
  alpha*(ginv11*((-(cdA111*chi) + 1.5*A11*dchi1)*qud11 + 
        (-(cdA112*chi) + 1.5*A12*dchi1)*qud21 + 
        (-(cdA113*chi) + 1.5*A13*dchi1)*qud31) + 
     ginv22*((-(cdA212*chi) + 1.5*A12*dchi2)*qud11 + 
        (-(cdA222*chi) + 1.5*A22*dchi2)*qud21 + 
        (-(cdA223*chi) + 1.5*A23*dchi2)*qud31) + 
     ginv33*((-(cdA313*chi) + 1.5*A13*dchi3)*qud11 + 
        (-(cdA323*chi) + 1.5*A23*dchi3)*qud21 + 
        (-(cdA333*chi) + 1.5*A33*dchi3)*qud31) + 
     chi*((0.66666666666666666667*dK1 - dTheta1)*qud11 + 
        (0.66666666666666666667*dK2 - dTheta2)*qud21 + 
        (0.66666666666666666667*dK3 - dTheta3)*qud31) + 
     ginv12*((-(cdA212*chi) + 1.5*A12*dchi2)*qud21 + 
        (-(cdA213*chi) + 1.5*A13*dchi2)*qud31 - 
        chi*((cdA112 + cdA211)*qud11 + cdA122*qud21 + cdA123*qud31) + 
        1.5*((A12*dchi1 + A11*dchi2)*qud11 + dchi1*(A22*qud21 + A23*qud31))\
) + ginv13*((-(cdA312*chi) + 1.5*A12*dchi3)*qud21 + 
        (-(cdA313*chi) + 1.5*A13*dchi3)*qud31 - 
        chi*((cdA113 + cdA311)*qud11 + cdA123*qud21 + cdA133*qud31) + 
        1.5*((A13*dchi1 + A11*dchi3)*qud11 + dchi1*(A23*qud21 + A33*qud31))\
) + ginv23*((-(cdA322*chi) + 1.5*A22*dchi3)*qud21 + 
        (-(cdA323*chi) + 1.5*A23*dchi3)*qud31 - 
        chi*((cdA213 + cdA312)*qud11 + cdA223*qud21 + cdA233*qud31) + 
        1.5*((A13*dchi2 + A12*dchi3)*qud11 + dchi2*(A23*qud21 + A33*qud31))\
) + 0.5*(kappa1*((G1 - Gfromg1)*qdd11 + (G2 - Gfromg2)*qdd12 + 
           (G3 - Gfromg3)*qdd13) - dG13*qdd13*sup1 - dG21*qdd11*sup2 + 
        (dGfromgdu22*qdd12 - dG23*qdd13)*sup2 + 
        (dGfromgdu31*qdd11 + dGfromgdu32*qdd12 - dG33*qdd13)*sup3 + 
        qdd11*((-dG11 + dGfromgdu11)*sup1 + dGfromgdu21*sup2 - 
           dG31*sup3) + qdd12*
         ((-dG12 + dGfromgdu12)*sup1 - dG22*sup2 - dG32*sup3))*pow2(chi) + 
     sup1*(-2.*AA11*qud11 + 0.5*dGfromgdu13*qdd13*pow2(chi))) + 
  sup2*(chi*(-(cdda12*qud11) - cdda22*qud21 - cdda23*qud31 + 
        alpha*qud21*Rf22) + alpha*
      (chi*(qud11*Rf12 + qud31*Rf23) + 0.5*dGfromgdu23*qdd13*pow2(chi))) + 
  sup3*(chi*(-(cdda13*qud11) - cdda23*qud21 - cdda33*qud31 + 
        alpha*qud21*Rf23) + alpha*
      (chi*(qud11*Rf13 + qud31*Rf33) + 0.5*dGfromgdu33*qdd13*pow2(chi)))
;

rACsA2
=
(qud12*(lieA11 + alpha*chi*Rf11) + 
     qud22*(lieA12 + alpha*(-2.*AA12 + chi*Rf12)) + 
     qud32*(lieA13 + alpha*(-2.*AA13 + chi*Rf13)))*sup1 + 
  qud12*((-(cdda11*chi) + A11*alpha*K)*sup1 + lieA12*sup2 + 
     (A13*alpha*K + lieA13)*sup3 + 
     alpha*((-2.*AA21 + A12*K)*sup2 - 2.*AA31*sup3)) + 
  qud22*((-(cdda12*chi) + A12*alpha*K)*sup1 + lieA22*sup2 + 
     (A23*alpha*K + lieA23)*sup3 + 
     alpha*((-2.*AA22 + A22*K)*sup2 - 2.*AA32*sup3)) + 
  qud32*((-(cdda13*chi) + A13*alpha*K)*sup1 + lieA23*sup2 + 
     (A33*alpha*K + lieA33)*sup3 + 
     alpha*((-2.*AA23 + A23*K)*sup2 - 2.*AA33*sup3)) + 
  alpha*(ginv11*((-(cdA111*chi) + 1.5*A11*dchi1)*qud12 + 
        (-(cdA112*chi) + 1.5*A12*dchi1)*qud22 + 
        (-(cdA113*chi) + 1.5*A13*dchi1)*qud32) + 
     ginv22*((-(cdA212*chi) + 1.5*A12*dchi2)*qud12 + 
        (-(cdA222*chi) + 1.5*A22*dchi2)*qud22 + 
        (-(cdA223*chi) + 1.5*A23*dchi2)*qud32) + 
     ginv33*((-(cdA313*chi) + 1.5*A13*dchi3)*qud12 + 
        (-(cdA323*chi) + 1.5*A23*dchi3)*qud22 + 
        (-(cdA333*chi) + 1.5*A33*dchi3)*qud32) + 
     chi*((0.66666666666666666667*dK1 - dTheta1)*qud12 + 
        (0.66666666666666666667*dK2 - dTheta2)*qud22 + 
        (0.66666666666666666667*dK3 - dTheta3)*qud32) + 
     ginv12*((-(cdA212*chi) + 1.5*A12*dchi2)*qud22 + 
        (-(cdA213*chi) + 1.5*A13*dchi2)*qud32 - 
        chi*((cdA112 + cdA211)*qud12 + cdA122*qud22 + cdA123*qud32) + 
        1.5*((A12*dchi1 + A11*dchi2)*qud12 + dchi1*(A22*qud22 + A23*qud32))\
) + ginv13*((-(cdA312*chi) + 1.5*A12*dchi3)*qud22 + 
        (-(cdA313*chi) + 1.5*A13*dchi3)*qud32 - 
        chi*((cdA113 + cdA311)*qud12 + cdA123*qud22 + cdA133*qud32) + 
        1.5*((A13*dchi1 + A11*dchi3)*qud12 + dchi1*(A23*qud22 + A33*qud32))\
) + ginv23*((-(cdA322*chi) + 1.5*A22*dchi3)*qud22 + 
        (-(cdA323*chi) + 1.5*A23*dchi3)*qud32 - 
        chi*((cdA213 + cdA312)*qud12 + cdA223*qud22 + cdA233*qud32) + 
        1.5*((A13*dchi2 + A12*dchi3)*qud12 + dchi2*(A23*qud22 + A33*qud32))\
) + 0.5*(kappa1*((G1 - Gfromg1)*qdd12 + (G2 - Gfromg2)*qdd22 + 
           (G3 - Gfromg3)*qdd23) - dG13*qdd23*sup1 - dG21*qdd12*sup2 + 
        (dGfromgdu22*qdd22 - dG23*qdd23)*sup2 + 
        (dGfromgdu31*qdd12 + dGfromgdu32*qdd22 - dG33*qdd23)*sup3 + 
        qdd12*((-dG11 + dGfromgdu11)*sup1 + dGfromgdu21*sup2 - 
           dG31*sup3) + qdd22*
         ((-dG12 + dGfromgdu12)*sup1 - dG22*sup2 - dG32*sup3))*pow2(chi) + 
     sup1*(-2.*AA11*qud12 + 0.5*dGfromgdu13*qdd23*pow2(chi))) + 
  sup2*(chi*(-(cdda12*qud12) - cdda22*qud22 - cdda23*qud32 + 
        alpha*qud22*Rf22) + alpha*
      (chi*(qud12*Rf12 + qud32*Rf23) + 0.5*dGfromgdu23*qdd23*pow2(chi))) + 
  sup3*(chi*(-(cdda13*qud12) - cdda23*qud22 - cdda33*qud32 + 
        alpha*qud22*Rf23) + alpha*
      (chi*(qud12*Rf13 + qud32*Rf33) + 0.5*dGfromgdu33*qdd23*pow2(chi)))
;

rACsA3
=
(qud13*(lieA11 + alpha*chi*Rf11) + 
     qud23*(lieA12 + alpha*(-2.*AA12 + chi*Rf12)) + 
     qud33*(lieA13 + alpha*(-2.*AA13 + chi*Rf13)))*sup1 + 
  qud13*((-(cdda11*chi) + A11*alpha*K)*sup1 + lieA12*sup2 + 
     (A13*alpha*K + lieA13)*sup3 + 
     alpha*((-2.*AA21 + A12*K)*sup2 - 2.*AA31*sup3)) + 
  qud23*((-(cdda12*chi) + A12*alpha*K)*sup1 + lieA22*sup2 + 
     (A23*alpha*K + lieA23)*sup3 + 
     alpha*((-2.*AA22 + A22*K)*sup2 - 2.*AA32*sup3)) + 
  qud33*((-(cdda13*chi) + A13*alpha*K)*sup1 + lieA23*sup2 + 
     (A33*alpha*K + lieA33)*sup3 + 
     alpha*((-2.*AA23 + A23*K)*sup2 - 2.*AA33*sup3)) + 
  alpha*(ginv11*((-(cdA111*chi) + 1.5*A11*dchi1)*qud13 + 
        (-(cdA112*chi) + 1.5*A12*dchi1)*qud23 + 
        (-(cdA113*chi) + 1.5*A13*dchi1)*qud33) + 
     ginv22*((-(cdA212*chi) + 1.5*A12*dchi2)*qud13 + 
        (-(cdA222*chi) + 1.5*A22*dchi2)*qud23 + 
        (-(cdA223*chi) + 1.5*A23*dchi2)*qud33) + 
     ginv33*((-(cdA313*chi) + 1.5*A13*dchi3)*qud13 + 
        (-(cdA323*chi) + 1.5*A23*dchi3)*qud23 + 
        (-(cdA333*chi) + 1.5*A33*dchi3)*qud33) + 
     chi*((0.66666666666666666667*dK1 - dTheta1)*qud13 + 
        (0.66666666666666666667*dK2 - dTheta2)*qud23 + 
        (0.66666666666666666667*dK3 - dTheta3)*qud33) + 
     ginv12*((-(cdA212*chi) + 1.5*A12*dchi2)*qud23 + 
        (-(cdA213*chi) + 1.5*A13*dchi2)*qud33 - 
        chi*((cdA112 + cdA211)*qud13 + cdA122*qud23 + cdA123*qud33) + 
        1.5*((A12*dchi1 + A11*dchi2)*qud13 + dchi1*(A22*qud23 + A23*qud33))\
) + ginv13*((-(cdA312*chi) + 1.5*A12*dchi3)*qud23 + 
        (-(cdA313*chi) + 1.5*A13*dchi3)*qud33 - 
        chi*((cdA113 + cdA311)*qud13 + cdA123*qud23 + cdA133*qud33) + 
        1.5*((A13*dchi1 + A11*dchi3)*qud13 + dchi1*(A23*qud23 + A33*qud33))\
) + ginv23*((-(cdA322*chi) + 1.5*A22*dchi3)*qud23 + 
        (-(cdA323*chi) + 1.5*A23*dchi3)*qud33 - 
        chi*((cdA213 + cdA312)*qud13 + cdA223*qud23 + cdA233*qud33) + 
        1.5*((A13*dchi2 + A12*dchi3)*qud13 + dchi2*(A23*qud23 + A33*qud33))\
) + 0.5*(kappa1*((G1 - Gfromg1)*qdd13 + (G2 - Gfromg2)*qdd23 + 
           (G3 - Gfromg3)*qdd33) - dG13*qdd33*sup1 - dG21*qdd13*sup2 + 
        (dGfromgdu22*qdd23 - dG23*qdd33)*sup2 + 
        (dGfromgdu31*qdd13 + dGfromgdu32*qdd23 - dG33*qdd33)*sup3 + 
        qdd13*((-dG11 + dGfromgdu11)*sup1 + dGfromgdu21*sup2 - 
           dG31*sup3) + qdd23*
         ((-dG12 + dGfromgdu12)*sup1 - dG22*sup2 - dG32*sup3))*pow2(chi) + 
     sup1*(-2.*AA11*qud13 + 0.5*dGfromgdu13*qdd33*pow2(chi))) + 
  sup2*(chi*(-(cdda12*qud13) - cdda22*qud23 - cdda23*qud33 + 
        alpha*qud23*Rf22) + alpha*
      (chi*(qud13*Rf12 + qud33*Rf23) + 0.5*dGfromgdu23*qdd33*pow2(chi))) + 
  sup3*(chi*(-(cdda13*qud13) - cdda23*qud23 - cdda33*qud33 + 
        alpha*qud23*Rf23) + alpha*
      (chi*(qud13*Rf13 + qud33*Rf33) + 0.5*dGfromgdu33*qdd33*pow2(chi)))
;

rACABTF11
=
-(qPhysuudd1211*(2.*cdda12*chi + alpha*(AA21 + cdA112*sup1))) + 
  qPhysuudd3311*(-(cdda33*chi) + lieA33 + 
     alpha*(0.66666666666666666667*A33*K + cdA313*sup1 + cdA323*sup2)) + 
  qPhysuudd1111*(-(cdda11*chi) + lieA11 + 
     alpha*(-AA11 + 0.66666666666666666667*A11*K + cdA112*sup2 + 
        cdA113*sup3)) + qPhysuudd1211*
   (2.*lieA12 + alpha*(-AA12 + 1.3333333333333333333*A12*K + cdA211*sup1 + 
        cdA122*sup2 + cdA123*sup3)) + 
  qPhysuudd1311*(2.*(-(cdda13*chi) + lieA13) + 
     alpha*(-AA31 + 1.3333333333333333333*A13*K + cdA311*sup1 + 
        cdA123*sup2 + cdA133*sup3)) + 
  qPhysuudd2211*(-(cdda22*chi) + lieA22 + 
     alpha*(0.66666666666666666667*A22*K + cdA212*sup1 + cdA223*sup3)) + 
  qPhysuudd2311*(2.*(-(cdda23*chi) + lieA23) + 
     alpha*(-AA32 + 1.3333333333333333333*A23*K + cdA213*sup1 + 
        cdA322*sup2 + cdA233*sup3)) - 
  alpha*(AA13*qPhysuudd1311 + AA22*qPhysuudd2211 + AA23*qPhysuudd2311 + 
     AA33*qPhysuudd3311 + qPhysuudd1111*(cdA211*sup2 + cdA311*sup3)) + 
  alpha*(-((2.*cdA213*qPhysuudd1311 + 
          (0.5*(A12*dchi1*qPhysuudd1111 + A23*dchi3*qPhysuudd3311))/chi)*
        sup2) - qPhysuudd3311*((cdA133 + (0.5*A13*dchi3)/chi)*sup1 + 
        cdA233*sup2) - 2.*cdA312*qPhysuudd1211*sup3 + 
     qPhysuudd1211*((-cdA212 + (0.5*A12*dchi2)/chi)*sup2 + cdA213*sup3) + 
     qPhysuudd1311*((-cdA113 + (0.5*A13*dchi1)/chi)*sup1 + cdA312*sup2 - 
        cdA313*sup3) - qPhysuudd2211*
      ((cdA122 + (0.5*A12*dchi2)/chi)*sup1 + cdA322*sup3) + 
     qPhysuudd2311*((cdA312 + (A23*dchi1)/chi)*sup1 + 
        (0.5*A23*dchi2*sup2)/chi - cdA323*sup3) - 
     qPhysuudd2311*((2.*cdA123 + (0.5*A13*dchi2)/chi)*sup1 + cdA223*sup2 + 
        (0.5*A33*dchi2*sup3)/chi) + 
     ((-0.5*A22*dchi1*qPhysuudd1211 + A13*dchi2*qPhysuudd1311)*sup2 + 
        (A12*dchi3*qPhysuudd1211 - 
           0.5*dchi1*(A13*qPhysuudd1111 + A23*qPhysuudd1211))*sup3 + 
        0.5*(((A12*dchi1 - A11*dchi2)*qPhysuudd1211 - 
              dchi3*(A11*qPhysuudd1311 + A12*qPhysuudd2311) + 
              dchi1*(A22*qPhysuudd2211 + A33*qPhysuudd3311))*sup1 + 
           (-((A23*dchi1 + A12*dchi3)*qPhysuudd1311) - 
              A22*dchi3*qPhysuudd2311 + 
              dchi2*(A11*qPhysuudd1111 + A33*qPhysuudd3311))*sup2 + 
           (-(A33*dchi1*qPhysuudd1311) + 
              A13*(-(dchi2*qPhysuudd1211) + dchi3*qPhysuudd1311) + 
              dchi3*(A11*qPhysuudd1111 + A22*qPhysuudd2211) + 
              A23*(-(dchi2*qPhysuudd2211) + dchi3*qPhysuudd2311))*sup3))/chi)
;

rACABTF12
=
-(qPhysuudd1212*(2.*cdda12*chi + alpha*(AA21 + cdA112*sup1))) + 
  qPhysuudd3312*(-(cdda33*chi) + lieA33 + 
     alpha*(0.66666666666666666667*A33*K + cdA313*sup1 + cdA323*sup2)) + 
  qPhysuudd1112*(-(cdda11*chi) + lieA11 + 
     alpha*(-AA11 + 0.66666666666666666667*A11*K + cdA112*sup2 + 
        cdA113*sup3)) + qPhysuudd1212*
   (2.*lieA12 + alpha*(-AA12 + 1.3333333333333333333*A12*K + cdA211*sup1 + 
        cdA122*sup2 + cdA123*sup3)) + 
  qPhysuudd1312*(2.*(-(cdda13*chi) + lieA13) + 
     alpha*(-AA31 + 1.3333333333333333333*A13*K + cdA311*sup1 + 
        cdA123*sup2 + cdA133*sup3)) + 
  qPhysuudd2212*(-(cdda22*chi) + lieA22 + 
     alpha*(0.66666666666666666667*A22*K + cdA212*sup1 + cdA223*sup3)) + 
  qPhysuudd2312*(2.*(-(cdda23*chi) + lieA23) + 
     alpha*(-AA32 + 1.3333333333333333333*A23*K + cdA213*sup1 + 
        cdA322*sup2 + cdA233*sup3)) - 
  alpha*(AA13*qPhysuudd1312 + AA22*qPhysuudd2212 + AA23*qPhysuudd2312 + 
     AA33*qPhysuudd3312 + qPhysuudd1112*(cdA211*sup2 + cdA311*sup3)) + 
  alpha*(-((2.*cdA213*qPhysuudd1312 + 
          (0.5*(A12*dchi1*qPhysuudd1112 + A23*dchi3*qPhysuudd3312))/chi)*
        sup2) - qPhysuudd3312*((cdA133 + (0.5*A13*dchi3)/chi)*sup1 + 
        cdA233*sup2) - 2.*cdA312*qPhysuudd1212*sup3 + 
     qPhysuudd1212*((-cdA212 + (0.5*A12*dchi2)/chi)*sup2 + cdA213*sup3) + 
     qPhysuudd1312*((-cdA113 + (0.5*A13*dchi1)/chi)*sup1 + cdA312*sup2 - 
        cdA313*sup3) - qPhysuudd2212*
      ((cdA122 + (0.5*A12*dchi2)/chi)*sup1 + cdA322*sup3) + 
     qPhysuudd2312*((cdA312 + (A23*dchi1)/chi)*sup1 + 
        (0.5*A23*dchi2*sup2)/chi - cdA323*sup3) - 
     qPhysuudd2312*((2.*cdA123 + (0.5*A13*dchi2)/chi)*sup1 + cdA223*sup2 + 
        (0.5*A33*dchi2*sup3)/chi) + 
     ((-0.5*A22*dchi1*qPhysuudd1212 + A13*dchi2*qPhysuudd1312)*sup2 + 
        (A12*dchi3*qPhysuudd1212 - 
           0.5*dchi1*(A13*qPhysuudd1112 + A23*qPhysuudd1212))*sup3 + 
        0.5*(((A12*dchi1 - A11*dchi2)*qPhysuudd1212 - 
              dchi3*(A11*qPhysuudd1312 + A12*qPhysuudd2312) + 
              dchi1*(A22*qPhysuudd2212 + A33*qPhysuudd3312))*sup1 + 
           (-((A23*dchi1 + A12*dchi3)*qPhysuudd1312) - 
              A22*dchi3*qPhysuudd2312 + 
              dchi2*(A11*qPhysuudd1112 + A33*qPhysuudd3312))*sup2 + 
           (-(A33*dchi1*qPhysuudd1312) + 
              A13*(-(dchi2*qPhysuudd1212) + dchi3*qPhysuudd1312) + 
              dchi3*(A11*qPhysuudd1112 + A22*qPhysuudd2212) + 
              A23*(-(dchi2*qPhysuudd2212) + dchi3*qPhysuudd2312))*sup3))/chi)
;

rACABTF13
=
-(qPhysuudd1213*(2.*cdda12*chi + alpha*(AA21 + cdA112*sup1))) + 
  qPhysuudd3313*(-(cdda33*chi) + lieA33 + 
     alpha*(0.66666666666666666667*A33*K + cdA313*sup1 + cdA323*sup2)) + 
  qPhysuudd1113*(-(cdda11*chi) + lieA11 + 
     alpha*(-AA11 + 0.66666666666666666667*A11*K + cdA112*sup2 + 
        cdA113*sup3)) + qPhysuudd1213*
   (2.*lieA12 + alpha*(-AA12 + 1.3333333333333333333*A12*K + cdA211*sup1 + 
        cdA122*sup2 + cdA123*sup3)) + 
  qPhysuudd1313*(2.*(-(cdda13*chi) + lieA13) + 
     alpha*(-AA31 + 1.3333333333333333333*A13*K + cdA311*sup1 + 
        cdA123*sup2 + cdA133*sup3)) + 
  qPhysuudd2213*(-(cdda22*chi) + lieA22 + 
     alpha*(0.66666666666666666667*A22*K + cdA212*sup1 + cdA223*sup3)) + 
  qPhysuudd2313*(2.*(-(cdda23*chi) + lieA23) + 
     alpha*(-AA32 + 1.3333333333333333333*A23*K + cdA213*sup1 + 
        cdA322*sup2 + cdA233*sup3)) - 
  alpha*(AA13*qPhysuudd1313 + AA22*qPhysuudd2213 + AA23*qPhysuudd2313 + 
     AA33*qPhysuudd3313 + qPhysuudd1113*(cdA211*sup2 + cdA311*sup3)) + 
  alpha*(-((2.*cdA213*qPhysuudd1313 + 
          (0.5*(A12*dchi1*qPhysuudd1113 + A23*dchi3*qPhysuudd3313))/chi)*
        sup2) - qPhysuudd3313*((cdA133 + (0.5*A13*dchi3)/chi)*sup1 + 
        cdA233*sup2) - 2.*cdA312*qPhysuudd1213*sup3 + 
     qPhysuudd1213*((-cdA212 + (0.5*A12*dchi2)/chi)*sup2 + cdA213*sup3) + 
     qPhysuudd1313*((-cdA113 + (0.5*A13*dchi1)/chi)*sup1 + cdA312*sup2 - 
        cdA313*sup3) - qPhysuudd2213*
      ((cdA122 + (0.5*A12*dchi2)/chi)*sup1 + cdA322*sup3) + 
     qPhysuudd2313*((cdA312 + (A23*dchi1)/chi)*sup1 + 
        (0.5*A23*dchi2*sup2)/chi - cdA323*sup3) - 
     qPhysuudd2313*((2.*cdA123 + (0.5*A13*dchi2)/chi)*sup1 + cdA223*sup2 + 
        (0.5*A33*dchi2*sup3)/chi) + 
     ((-0.5*A22*dchi1*qPhysuudd1213 + A13*dchi2*qPhysuudd1313)*sup2 + 
        (A12*dchi3*qPhysuudd1213 - 
           0.5*dchi1*(A13*qPhysuudd1113 + A23*qPhysuudd1213))*sup3 + 
        0.5*(((A12*dchi1 - A11*dchi2)*qPhysuudd1213 - 
              dchi3*(A11*qPhysuudd1313 + A12*qPhysuudd2313) + 
              dchi1*(A22*qPhysuudd2213 + A33*qPhysuudd3313))*sup1 + 
           (-((A23*dchi1 + A12*dchi3)*qPhysuudd1313) - 
              A22*dchi3*qPhysuudd2313 + 
              dchi2*(A11*qPhysuudd1113 + A33*qPhysuudd3313))*sup2 + 
           (-(A33*dchi1*qPhysuudd1313) + 
              A13*(-(dchi2*qPhysuudd1213) + dchi3*qPhysuudd1313) + 
              dchi3*(A11*qPhysuudd1113 + A22*qPhysuudd2213) + 
              A23*(-(dchi2*qPhysuudd2213) + dchi3*qPhysuudd2313))*sup3))/chi)
;

rACABTF22
=
-(qPhysuudd1222*(2.*cdda12*chi + alpha*(AA21 + cdA112*sup1))) + 
  qPhysuudd3322*(-(cdda33*chi) + lieA33 + 
     alpha*(0.66666666666666666667*A33*K + cdA313*sup1 + cdA323*sup2)) + 
  qPhysuudd1122*(-(cdda11*chi) + lieA11 + 
     alpha*(-AA11 + 0.66666666666666666667*A11*K + cdA112*sup2 + 
        cdA113*sup3)) + qPhysuudd1222*
   (2.*lieA12 + alpha*(-AA12 + 1.3333333333333333333*A12*K + cdA211*sup1 + 
        cdA122*sup2 + cdA123*sup3)) + 
  qPhysuudd1322*(2.*(-(cdda13*chi) + lieA13) + 
     alpha*(-AA31 + 1.3333333333333333333*A13*K + cdA311*sup1 + 
        cdA123*sup2 + cdA133*sup3)) + 
  qPhysuudd2222*(-(cdda22*chi) + lieA22 + 
     alpha*(0.66666666666666666667*A22*K + cdA212*sup1 + cdA223*sup3)) + 
  qPhysuudd2322*(2.*(-(cdda23*chi) + lieA23) + 
     alpha*(-AA32 + 1.3333333333333333333*A23*K + cdA213*sup1 + 
        cdA322*sup2 + cdA233*sup3)) - 
  alpha*(AA13*qPhysuudd1322 + AA22*qPhysuudd2222 + AA23*qPhysuudd2322 + 
     AA33*qPhysuudd3322 + qPhysuudd1122*(cdA211*sup2 + cdA311*sup3)) + 
  alpha*(-((2.*cdA213*qPhysuudd1322 + 
          (0.5*(A12*dchi1*qPhysuudd1122 + A23*dchi3*qPhysuudd3322))/chi)*
        sup2) - qPhysuudd3322*((cdA133 + (0.5*A13*dchi3)/chi)*sup1 + 
        cdA233*sup2) - 2.*cdA312*qPhysuudd1222*sup3 + 
     qPhysuudd1222*((-cdA212 + (0.5*A12*dchi2)/chi)*sup2 + cdA213*sup3) + 
     qPhysuudd1322*((-cdA113 + (0.5*A13*dchi1)/chi)*sup1 + cdA312*sup2 - 
        cdA313*sup3) - qPhysuudd2222*
      ((cdA122 + (0.5*A12*dchi2)/chi)*sup1 + cdA322*sup3) + 
     qPhysuudd2322*((cdA312 + (A23*dchi1)/chi)*sup1 + 
        (0.5*A23*dchi2*sup2)/chi - cdA323*sup3) - 
     qPhysuudd2322*((2.*cdA123 + (0.5*A13*dchi2)/chi)*sup1 + cdA223*sup2 + 
        (0.5*A33*dchi2*sup3)/chi) + 
     ((-0.5*A22*dchi1*qPhysuudd1222 + A13*dchi2*qPhysuudd1322)*sup2 + 
        (A12*dchi3*qPhysuudd1222 - 
           0.5*dchi1*(A13*qPhysuudd1122 + A23*qPhysuudd1222))*sup3 + 
        0.5*(((A12*dchi1 - A11*dchi2)*qPhysuudd1222 - 
              dchi3*(A11*qPhysuudd1322 + A12*qPhysuudd2322) + 
              dchi1*(A22*qPhysuudd2222 + A33*qPhysuudd3322))*sup1 + 
           (-((A23*dchi1 + A12*dchi3)*qPhysuudd1322) - 
              A22*dchi3*qPhysuudd2322 + 
              dchi2*(A11*qPhysuudd1122 + A33*qPhysuudd3322))*sup2 + 
           (-(A33*dchi1*qPhysuudd1322) + 
              A13*(-(dchi2*qPhysuudd1222) + dchi3*qPhysuudd1322) + 
              dchi3*(A11*qPhysuudd1122 + A22*qPhysuudd2222) + 
              A23*(-(dchi2*qPhysuudd2222) + dchi3*qPhysuudd2322))*sup3))/chi)
;

rACABTF23
=
-(qPhysuudd1223*(2.*cdda12*chi + alpha*(AA21 + cdA112*sup1))) + 
  qPhysuudd3323*(-(cdda33*chi) + lieA33 + 
     alpha*(0.66666666666666666667*A33*K + cdA313*sup1 + cdA323*sup2)) + 
  qPhysuudd1123*(-(cdda11*chi) + lieA11 + 
     alpha*(-AA11 + 0.66666666666666666667*A11*K + cdA112*sup2 + 
        cdA113*sup3)) + qPhysuudd1223*
   (2.*lieA12 + alpha*(-AA12 + 1.3333333333333333333*A12*K + cdA211*sup1 + 
        cdA122*sup2 + cdA123*sup3)) + 
  qPhysuudd1323*(2.*(-(cdda13*chi) + lieA13) + 
     alpha*(-AA31 + 1.3333333333333333333*A13*K + cdA311*sup1 + 
        cdA123*sup2 + cdA133*sup3)) + 
  qPhysuudd2223*(-(cdda22*chi) + lieA22 + 
     alpha*(0.66666666666666666667*A22*K + cdA212*sup1 + cdA223*sup3)) + 
  qPhysuudd2323*(2.*(-(cdda23*chi) + lieA23) + 
     alpha*(-AA32 + 1.3333333333333333333*A23*K + cdA213*sup1 + 
        cdA322*sup2 + cdA233*sup3)) - 
  alpha*(AA13*qPhysuudd1323 + AA22*qPhysuudd2223 + AA23*qPhysuudd2323 + 
     AA33*qPhysuudd3323 + qPhysuudd1123*(cdA211*sup2 + cdA311*sup3)) + 
  alpha*(-((2.*cdA213*qPhysuudd1323 + 
          (0.5*(A12*dchi1*qPhysuudd1123 + A23*dchi3*qPhysuudd3323))/chi)*
        sup2) - qPhysuudd3323*((cdA133 + (0.5*A13*dchi3)/chi)*sup1 + 
        cdA233*sup2) - 2.*cdA312*qPhysuudd1223*sup3 + 
     qPhysuudd1223*((-cdA212 + (0.5*A12*dchi2)/chi)*sup2 + cdA213*sup3) + 
     qPhysuudd1323*((-cdA113 + (0.5*A13*dchi1)/chi)*sup1 + cdA312*sup2 - 
        cdA313*sup3) - qPhysuudd2223*
      ((cdA122 + (0.5*A12*dchi2)/chi)*sup1 + cdA322*sup3) + 
     qPhysuudd2323*((cdA312 + (A23*dchi1)/chi)*sup1 + 
        (0.5*A23*dchi2*sup2)/chi - cdA323*sup3) - 
     qPhysuudd2323*((2.*cdA123 + (0.5*A13*dchi2)/chi)*sup1 + cdA223*sup2 + 
        (0.5*A33*dchi2*sup3)/chi) + 
     ((-0.5*A22*dchi1*qPhysuudd1223 + A13*dchi2*qPhysuudd1323)*sup2 + 
        (A12*dchi3*qPhysuudd1223 - 
           0.5*dchi1*(A13*qPhysuudd1123 + A23*qPhysuudd1223))*sup3 + 
        0.5*(((A12*dchi1 - A11*dchi2)*qPhysuudd1223 - 
              dchi3*(A11*qPhysuudd1323 + A12*qPhysuudd2323) + 
              dchi1*(A22*qPhysuudd2223 + A33*qPhysuudd3323))*sup1 + 
           (-((A23*dchi1 + A12*dchi3)*qPhysuudd1323) - 
              A22*dchi3*qPhysuudd2323 + 
              dchi2*(A11*qPhysuudd1123 + A33*qPhysuudd3323))*sup2 + 
           (-(A33*dchi1*qPhysuudd1323) + 
              A13*(-(dchi2*qPhysuudd1223) + dchi3*qPhysuudd1323) + 
              dchi3*(A11*qPhysuudd1123 + A22*qPhysuudd2223) + 
              A23*(-(dchi2*qPhysuudd2223) + dchi3*qPhysuudd2323))*sup3))/chi)
;

rACABTF33
=
-(qPhysuudd1233*(2.*cdda12*chi + alpha*(AA21 + cdA112*sup1))) + 
  qPhysuudd3333*(-(cdda33*chi) + lieA33 + 
     alpha*(0.66666666666666666667*A33*K + cdA313*sup1 + cdA323*sup2)) + 
  qPhysuudd1133*(-(cdda11*chi) + lieA11 + 
     alpha*(-AA11 + 0.66666666666666666667*A11*K + cdA112*sup2 + 
        cdA113*sup3)) + qPhysuudd1233*
   (2.*lieA12 + alpha*(-AA12 + 1.3333333333333333333*A12*K + cdA211*sup1 + 
        cdA122*sup2 + cdA123*sup3)) + 
  qPhysuudd1333*(2.*(-(cdda13*chi) + lieA13) + 
     alpha*(-AA31 + 1.3333333333333333333*A13*K + cdA311*sup1 + 
        cdA123*sup2 + cdA133*sup3)) + 
  qPhysuudd2233*(-(cdda22*chi) + lieA22 + 
     alpha*(0.66666666666666666667*A22*K + cdA212*sup1 + cdA223*sup3)) + 
  qPhysuudd2333*(2.*(-(cdda23*chi) + lieA23) + 
     alpha*(-AA32 + 1.3333333333333333333*A23*K + cdA213*sup1 + 
        cdA322*sup2 + cdA233*sup3)) - 
  alpha*(AA13*qPhysuudd1333 + AA22*qPhysuudd2233 + AA23*qPhysuudd2333 + 
     AA33*qPhysuudd3333 + qPhysuudd1133*(cdA211*sup2 + cdA311*sup3)) + 
  alpha*(-((2.*cdA213*qPhysuudd1333 + 
          (0.5*(A12*dchi1*qPhysuudd1133 + A23*dchi3*qPhysuudd3333))/chi)*
        sup2) - qPhysuudd3333*((cdA133 + (0.5*A13*dchi3)/chi)*sup1 + 
        cdA233*sup2) - 2.*cdA312*qPhysuudd1233*sup3 + 
     qPhysuudd1233*((-cdA212 + (0.5*A12*dchi2)/chi)*sup2 + cdA213*sup3) + 
     qPhysuudd1333*((-cdA113 + (0.5*A13*dchi1)/chi)*sup1 + cdA312*sup2 - 
        cdA313*sup3) - qPhysuudd2233*
      ((cdA122 + (0.5*A12*dchi2)/chi)*sup1 + cdA322*sup3) + 
     qPhysuudd2333*((cdA312 + (A23*dchi1)/chi)*sup1 + 
        (0.5*A23*dchi2*sup2)/chi - cdA323*sup3) - 
     qPhysuudd2333*((2.*cdA123 + (0.5*A13*dchi2)/chi)*sup1 + cdA223*sup2 + 
        (0.5*A33*dchi2*sup3)/chi) + 
     ((-0.5*A22*dchi1*qPhysuudd1233 + A13*dchi2*qPhysuudd1333)*sup2 + 
        (A12*dchi3*qPhysuudd1233 - 
           0.5*dchi1*(A13*qPhysuudd1133 + A23*qPhysuudd1233))*sup3 + 
        0.5*(((A12*dchi1 - A11*dchi2)*qPhysuudd1233 - 
              dchi3*(A11*qPhysuudd1333 + A12*qPhysuudd2333) + 
              dchi1*(A22*qPhysuudd2233 + A33*qPhysuudd3333))*sup1 + 
           (-((A23*dchi1 + A12*dchi3)*qPhysuudd1333) - 
              A22*dchi3*qPhysuudd2333 + 
              dchi2*(A11*qPhysuudd1133 + A33*qPhysuudd3333))*sup2 + 
           (-(A33*dchi1*qPhysuudd1333) + 
              A13*(-(dchi2*qPhysuudd1233) + dchi3*qPhysuudd1333) + 
              dchi3*(A11*qPhysuudd1133 + A22*qPhysuudd2233) + 
              A23*(-(dchi2*qPhysuudd2233) + dchi3*qPhysuudd2333))*sup3))/chi)
;

}  /* function */

}
