
#ifndef PARAMETERS_H
#define PARAMETERS_H

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
#include <map>
using namespace std;
#else
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <map.h>
#endif

#include <mpi.h>

namespace parameters
{
    extern map<string,int> int_par;
    extern map<string,double> dou_par;
    extern map<string,string> str_par;
}
#endif   /* PARAMETERS_H */
