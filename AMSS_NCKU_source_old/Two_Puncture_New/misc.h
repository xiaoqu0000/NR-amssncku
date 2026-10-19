//$Id: misc.h,v 1.1.1.1 2017/09/25 02:35:55 zjcao Exp $
#ifndef MISC_H
#define MISC_H

#ifdef newc
#include <algorithm>   
#include <functional> 
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <complex>
using namespace std;
#else
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#endif

#include <mpi.h>

namespace misc
{
inline string&  lTrim(string   &ss)  {   
    string::iterator  p=find_if(ss.begin(),ss.end(),not1(ptr_fun<int,int>(isspace)));   
    ss.erase(ss.begin(),p);   
    return  ss;   
}   
inline string&  rTrim(string   &ss)  {    
    string::reverse_iterator  p=find_if(ss.rbegin(),ss.rend(),not1(ptr_fun<int,int>(isspace)));   
    ss.erase(p.base(),ss.end());   
    return   ss;   
}   
inline string& Trim(string   &st)  {   
    lTrim(rTrim(st));   
    return   st;   
}
int parse_parts(string str, string& sgrp, string& skey, string& sval, int& ind);
int parse_parts(string str, string& sgrp, string& skey, string& sval, int& ind1, int& ind2);
int parse_parts(string str, string& sgrp, string& skey, string& sval, int& ind1, int& ind2, int& ind3);
void rungekutta4(const int N, const double dT,double *f0,
		        double *f1,double *f_rhs,const int RK4);
}
#endif   /* MISC_H */
