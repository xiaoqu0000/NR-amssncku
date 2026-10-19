
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include <mpi.h>

#include "myglobal.h"

int CCTK_VInfo(const char *thorn, const char *format, ...)
{
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
   if (myrank !=0) return 0;
   
   va_list ap;
   va_start (ap, format);
   fprintf (stdout, "INFO (%s): ", thorn);
   vfprintf (stdout, format, ap);
   fprintf (stdout, "\n");
   va_end (ap);
   return 0;
}
int CCTK_VWarn (int level,
                int line,
                const char *file,
                const char *thorn,
                const char *format,
                ...) 
{  
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
   if (myrank !=0) return 0;
   
   va_list ap;
   va_start (ap, format);
   fprintf (stdout, "WARN (%s): ", thorn);
   vfprintf (stdout, format, ap);
   fprintf (stdout, "\n");
   va_end (ap);
   return 0;
}
