//$Id: misc.C,v 1.1.1.1 2017/09/25 02:35:55 zjcao Exp $
#ifdef newc
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <strstream>
#include <cmath>
using namespace std;
#else
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <string.h>
#include <math.h>
#endif
#include <mpi.h>

#include "misc.h"

#define PI M_PI

// pick out value from input string
int misc::parse_parts(string str, string& sgrp, string& skey, string& sval, int& ind)
{
      int pos1, pos2;
      string s0;
      
      ind = 0;
      
      // remove comments
      str = str.substr(0, str.find("#") );
      if ( rTrim(str).empty() ) return 0;   // continue;
      
      // parse {group, key, val}
      pos1 = str.find("::");  pos2 = str.find("=");
      if (pos1 == string::npos || pos2 == string::npos) return -1;
      
      s0 = str.substr(     0, pos1       );  sgrp = lTrim( s0 );
      s0 = str.substr(pos1+2, pos2-pos1-2);  skey = rTrim( s0 );
      s0 = str.substr(pos2+1)             ;  sval = Trim( s0 );  

      pos1 = sval.find("\"");  pos2 = sval.rfind("\"");
      if ( pos1 != string::npos ) {
         sval = sval.substr(1, pos2-1); 
      }

      pos1 = skey.find("[");  pos2 = skey.find("]");
      if ( pos1 != string::npos ) {
         s0   = skey.substr(0,pos1);
         ind = atoi( skey.substr(pos1+1 ,pos2-pos1-1).c_str() );
         skey  = s0; 
      }
      
      return 1;
}
int misc::parse_parts(string str, string& sgrp, string& skey, string& sval, int& ind1, int& ind2)
{
      int pos1, pos2;
      string s0,s1;
      
      ind1 = ind2 = 0;
      
      // remove comments
      str = str.substr(0, str.find("#") );
      if ( rTrim(str).empty() ) return 0;   // continue;
      
      // parse {group, key, val}
      pos1 = str.find("::");  pos2 = str.find("=");
      if (pos1 == string::npos || pos2 == string::npos) return -1;
      
      s0 = str.substr(     0, pos1       );  sgrp = lTrim( s0 );
      s0 = str.substr(pos1+2, pos2-pos1-2);  skey = rTrim( s0 );
      s0 = str.substr(pos2+1)             ;  sval = Trim( s0 );  

      pos1 = sval.find("\"");  pos2 = sval.rfind("\"");
      if ( pos1 != string::npos ) {
         sval = sval.substr(1, pos2-1); 
      }
      
      pos1 = skey.find("[");  pos2 = skey.find("]");
      if ( pos1 != string::npos ) {
         s0   = skey.substr(0,pos1);
	 s1   = skey.substr(pos2+1);
         ind1 = atoi( skey.substr(pos1+1 ,pos2-pos1-1).c_str() );
         skey  = s0; 
      }

      pos1 = s1.find("[");  pos2 = s1.find("]");
      if ( pos1 != string::npos ) {
         s0   = s1.substr(pos2+1);
         ind2 = atoi( s1.substr(pos1+1 ,pos2-pos1-1).c_str() );
      }
      
      return 1;
}
int misc::parse_parts(string str, string& sgrp, string& skey, string& sval, int& ind1, int& ind2, int& ind3)
{
      int pos1, pos2;
      string s0,s1;
      
      ind1 = ind2 = ind3 = 0;
      
      // remove comments
      str = str.substr(0, str.find("#") );
      if ( rTrim(str).empty() ) return 0;   // continue;
      
      // parse {group, key, val}
      pos1 = str.find("::");  pos2 = str.find("=");
      if (pos1 == string::npos || pos2 == string::npos) return -1;
      
      s0 = str.substr(     0, pos1       );  sgrp = lTrim( s0 );
      s0 = str.substr(pos1+2, pos2-pos1-2);  skey = rTrim( s0 );
      s0 = str.substr(pos2+1)             ;  sval = Trim( s0 );  

      pos1 = sval.find("\"");  pos2 = sval.rfind("\"");
      if ( pos1 != string::npos ) {
         sval = sval.substr(1, pos2-1); 
      }
      
      pos1 = skey.find("[");  pos2 = skey.find("]");
      if ( pos1 != string::npos ) {
         s0   = skey.substr(0,pos1);
	 s1   = skey.substr(pos2+1);
         ind1 = atoi( skey.substr(pos1+1 ,pos2-pos1-1).c_str() );
         skey  = s0; 
      }

      pos1 = s1.find("[");  pos2 = s1.find("]");
      if ( pos1 != string::npos ) {
         s0   = s1.substr(pos2+1);
         ind2 = atoi( s1.substr(pos1+1 ,pos2-pos1-1).c_str() );
      }

      pos1 = s0.find("[");  pos2 = s0.find("]");
      if ( pos1 != string::npos ) {
         ind3 = atoi( s0.substr(pos1+1 ,pos2-pos1-1).c_str() );
      }
      
      return 1;
}
void misc::rungekutta4(const int N, const double dT,double *f0,
		        double *f1,double *f_rhs,const int RK4)
{
  const double F1o6=1.0/6,HLF=0.5,TWO=2;
  switch(RK4)
  {
    case 0:
     for(int i=0;i<N;i++)f1[i] = f0[i]+HLF*dT*f_rhs[i];
     break;
    case 1:
     for(int i=0;i<N;i++){
       f_rhs[i] = f_rhs[i]+TWO*f1[i];
       f1[i] = f0[i]+HLF*dT*f1[i];
       }
     break;
    case 2:
     for(int i=0;i<N;i++){
       f_rhs[i] = f_rhs[i] + TWO*f1[i];
       f1[i] = f0[i]+dT*f1[i];
       }
     break;
    case 3:
     for(int i=0;i<N;i++)f1[i] = f0[i]+F1o6*dT*(f1[i]+f_rhs[i]);
     break;
    default:
     cout<<__FILE__<<"("<<__LINE__<<"), "<<__func__<<": something is wrong in RK4 counting!!"<<endl;
  }
}
