
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <new>
using namespace std;

#include "NullShellPatch2.h"
#include "Parallel.h"
#include "fmisc.h"
#include "misc.h"
#include "shellfunctions.h"
#include "NullEvol.h"
#include "NullNews.h"
#include "initial_null2.h"
#include "rungekutta4_rout.h"
#include "kodiss.h"

#define PI M_PI

#if 0
// for RT
void NullShellPatch2::Setup_Initial_Data(bool checkrun,double PhysTime)
{ 
  if(checkrun)
  {
  }
  else
  {
     MyList<ss_patch> *Pp=PatL;
     while(Pp)
     {
      MyList<Block> *BL=Pp->data->blb;
      while(BL)
      {
       Block *cg=BL->data;
       if(myrank == cg->rank) 
       {
          f_get_initial_null2(cg->shape,cg->X[0],cg->X[1],cg->X[2],
			            cg->fgfs[g220->sgfn],cg->fgfs[g230->sgfn],cg->fgfs[g330->sgfn],
				    Pp->data->sst,Rmin);
// for Theta_AB	  
	  f_get_gauge_g00_K(cg->shape,cg->X[0],cg->X[1],cg->X[2],
			  	   cg->fgfs[g220->sgfn],cg->fgfs[g230->sgfn],cg->fgfs[g330->sgfn],
                                   cg->fgfs[Theta22->sgfn],cg->fgfs[Theta23->sgfn],cg->fgfs[Theta33->sgfn],
			  	   cg->fgfs[g00->sgfn],Rmin);
       }
       if(BL == Pp->data->ble) break;
       BL=BL->next;
      }
      Pp=Pp->next;
     }
//Synchronize K
    Synch(g00List,Symmetry,g00wt,1);
    Pp=PatL;
    int IONE=1;
    while(Pp)
    {
      MyList<Block> *BP=Pp->data->blb;
      int fngfs = Pp->data->fngfs;
      while(BP)
      {
        Block *cg=BP->data;
        if(myrank == cg->rank) 
        {
	   f_get_gauge_g00(cg->shape,cg->X[0],cg->X[1],cg->X[2],
			  	   cg->fgfs[g220->sgfn],cg->fgfs[g230->sgfn],cg->fgfs[g330->sgfn],
                                   cg->fgfs[Theta22->sgfn],cg->fgfs[Theta23->sgfn],cg->fgfs[Theta33->sgfn],
			  	   cg->fgfs[g00->sgfn],Rmin,IONE);
	}
        if(BP==Pp->data->ble) break;
        BP=BP->next;
      }
      Pp=Pp->next;
    }
    Synch(ThetaList,Symmetry,Thetawt,3);
  }
}
#else
void NullShellPatch2::Setup_Initial_Data(bool checkrun, double PhysTime)
{
  if (checkrun)
  {
  }
  else
  {
    MyList<ss_patch> *Pp = PatL;
    while (Pp)
    {
      MyList<Block> *BL = Pp->data->blb;
      while (BL)
      {
        Block *cg = BL->data;
        if (myrank == cg->rank)
        {
          f_get_initial_null3(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                              cg->fgfs[g220->sgfn], cg->fgfs[g230->sgfn], cg->fgfs[g330->sgfn],
                              Pp->data->sst, Rmin);
        }
        if (BL == Pp->data->ble)
          break;
        BL = BL->next;
      }
      Pp = Pp->next;
    }
  }
}
#endif
void NullShellPatch2::Step(double dT, double PhysTime, monitor *ErrorMonitor)
{
  int iter_count = 0; // count RK4 substeps
  int pre = 0, cor = 1;
  int ERROR = 0;
  double TT = PhysTime;
  double neps = -0.05;
  MyList<ss_patch> *sPp;

  // Predictor
  HyperSlice(dT, TT, ErrorMonitor, iter_count);
  {
    sPp = PatL;
    while (sPp)
    {
      MyList<Block> *BP = sPp->data->blb;
      while (BP)
      {
        Block *cg = BP->data;
        if (myrank == cg->rank)
        {
          // rhs calculation
          f_array_copy(cg->shape, cg->fgfs[g22_rhs->sgfn], cg->fgfs[Theta22->sgfn]);
          f_array_copy(cg->shape, cg->fgfs[g23_rhs->sgfn], cg->fgfs[Theta23->sgfn]);
          f_array_copy(cg->shape, cg->fgfs[g33_rhs->sgfn], cg->fgfs[Theta33->sgfn]);
          f_kodis_sh(cg->shape, cg->X[0], cg->X[1], cg->X[2], cg->fgfs[g220->sgfn], cg->fgfs[g22_rhs->sgfn],
                     Thetawt[0], Symmetry, neps, sPp->data->sst);
          f_kodis_sh(cg->shape, cg->X[0], cg->X[1], cg->X[2], cg->fgfs[g230->sgfn], cg->fgfs[g23_rhs->sgfn],
                     Thetawt[1], Symmetry, neps, sPp->data->sst);
          f_kodis_sh(cg->shape, cg->X[0], cg->X[1], cg->X[2], cg->fgfs[g330->sgfn], cg->fgfs[g33_rhs->sgfn],
                     Thetawt[2], Symmetry, neps, sPp->data->sst);
          f_rungekutta4_rout(cg->shape, dT, cg->fgfs[g220->sgfn], cg->fgfs[g22->sgfn], cg->fgfs[g22_rhs->sgfn],
                             iter_count);
          f_rungekutta4_rout(cg->shape, dT, cg->fgfs[g230->sgfn], cg->fgfs[g23->sgfn], cg->fgfs[g23_rhs->sgfn],
                             iter_count);
          f_rungekutta4_rout(cg->shape, dT, cg->fgfs[g330->sgfn], cg->fgfs[g33->sgfn], cg->fgfs[g33_rhs->sgfn],
                             iter_count);
        }
        if (BP == sPp->data->ble)
          break;
        BP = BP->next;
      }
      sPp = sPp->next;
    }
  }
  Synch(SynchList_pre, Symmetry, Thetawt, 3);
  //    Synch(SynchList_pre,Symmetry,g00wt,1);

  // corrector
  for (iter_count = 1; iter_count < 4; iter_count++)
  {
    // for RK4: t0, t0+dt/2, t0+dt/2, t0+dt;
    if (iter_count == 1 || iter_count == 3)
      TT += dT / 2;
    HyperSlice(dT, TT, ErrorMonitor, iter_count);
    {
      sPp = PatL;
      while (sPp)
      {
        MyList<Block> *BP = sPp->data->blb;
        while (BP)
        {
          Block *cg = BP->data;
          if (myrank == cg->rank)
          {
            // rhs calculation
            f_array_copy(cg->shape, cg->fgfs[g221->sgfn], cg->fgfs[Theta22->sgfn]);
            f_array_copy(cg->shape, cg->fgfs[g231->sgfn], cg->fgfs[Theta23->sgfn]);
            f_array_copy(cg->shape, cg->fgfs[g331->sgfn], cg->fgfs[Theta33->sgfn]);
            f_kodis_sh(cg->shape, cg->X[0], cg->X[1], cg->X[2], cg->fgfs[g22->sgfn], cg->fgfs[g221->sgfn],
                       Thetawt[0], Symmetry, neps, sPp->data->sst);
            f_kodis_sh(cg->shape, cg->X[0], cg->X[1], cg->X[2], cg->fgfs[g23->sgfn], cg->fgfs[g231->sgfn],
                       Thetawt[1], Symmetry, neps, sPp->data->sst);
            f_kodis_sh(cg->shape, cg->X[0], cg->X[1], cg->X[2], cg->fgfs[g33->sgfn], cg->fgfs[g331->sgfn],
                       Thetawt[2], Symmetry, neps, sPp->data->sst);
            f_rungekutta4_rout(cg->shape, dT, cg->fgfs[g220->sgfn], cg->fgfs[g221->sgfn], cg->fgfs[g22_rhs->sgfn],
                               iter_count);
            f_rungekutta4_rout(cg->shape, dT, cg->fgfs[g230->sgfn], cg->fgfs[g231->sgfn], cg->fgfs[g23_rhs->sgfn],
                               iter_count);
            f_rungekutta4_rout(cg->shape, dT, cg->fgfs[g330->sgfn], cg->fgfs[g331->sgfn], cg->fgfs[g33_rhs->sgfn],
                               iter_count);
          }
          if (iter_count < 3)
            cg->swapList(SynchList_cor, SynchList_pre, myrank);
          else
          {
            cg->swapList(StateList, SynchList_cor, myrank);
            cg->swapList(OldStateList, SynchList_cor, myrank);
          }
          if (BP == sPp->data->ble)
            break;
          BP = BP->next;
        }
        sPp = sPp->next;
      }
    }
    if (iter_count < 3)
      Synch(SynchList_pre, Symmetry, Thetawt, 3);
    else
      Synch(StateList, Symmetry, Thetawt, 3);
    //  if( iter_count < 3 ) Synch(SynchList_pre,Symmetry,g00wt,1);
    //  else                 Synch(StateList,Symmetry,g00wt,1);
  }
}
// really ODEs, so we do not need Synch in this routine at all
#if 0
void NullShellPatch2::HyperSlice(double dT,double PhysTime,monitor *ErrorMonitor,int RK_count)
{ 	
    int ERROR=0;
    Null_Boundary(PhysTime);
#if 1
    MyList<ss_patch> *sPp;

// evolve g01 
    sPp=PatL;
    while(sPp)
    {
      MyList<Block> *BP=sPp->data->blb;
      int fngfs = sPp->data->fngfs;
      while(BP)
      {
        Block *cg=BP->data;
        if(myrank == cg->rank) 
        {
	   if(RK_count==0)
	   {
		if(f_NullEvol_g01(cg->shape,cg->X[0],cg->X[1],cg->X[2],
			  	  cg->fgfs[g220->sgfn],cg->fgfs[g230->sgfn],cg->fgfs[g330->sgfn],
				  cg->fgfs[g01->sgfn],Rmin))		    
	       	{
			cout<<"find NaN of g01 in NullShell domain: sst = "<<sPp->data->sst<<", ("<<cg->bbox[0]<<":"<<cg->bbox[3]<<","
			<<cg->bbox[1]<<":"<<cg->bbox[4]<<","<<cg->bbox[2]<<":"<<cg->bbox[5]<<")"<<endl;
			ERROR = 1;
	        }
	   }
	   else
	   {
		if(f_NullEvol_g01(cg->shape,cg->X[0],cg->X[1],cg->X[2],
			  	  cg->fgfs[g22->sgfn],cg->fgfs[g23->sgfn],cg->fgfs[g33->sgfn],
				  cg->fgfs[g01->sgfn],Rmin))		    
	       	{
			cout<<"find NaN of g01 in NullShell domain: sst = "<<sPp->data->sst<<", ("<<cg->bbox[0]<<":"<<cg->bbox[3]<<","
			<<cg->bbox[1]<<":"<<cg->bbox[4]<<","<<cg->bbox[2]<<":"<<cg->bbox[5]<<")"<<endl;
			ERROR = 1;
	        }
	   }
	}
        if(BP==sPp->data->ble) break;
        BP=BP->next;
      }
      sPp=sPp->next;
    }
    Synch(g01List,Symmetry,g01wt,1);
    if(RK_count==3) Dump_Data(g01List,0,PhysTime,dT);
//check error information
  {int erh=ERROR;MPI_Allreduce(&erh,&ERROR,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD); }
  if(ERROR)
  {
	 if(RK_count==0)  Dump_Data(StateList,0,PhysTime,dT);
	 else             Dump_Data(SynchList_pre,0,PhysTime,dT);
         if(myrank == 0)
	 {
            if(ErrorMonitor && ErrorMonitor->outfile) 
	       ErrorMonitor->outfile<<"find NaN in beta on NullShell Patches at t = "<<PhysTime<<endl;
	    else
	                        cout<<"find NaN in beta on NullShell Patches at t = "<<PhysTime<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	 }
  }

// evolve p02, p03, g02 and g03
    sPp=PatL;
    while(sPp)
    {
      MyList<Block> *BP=sPp->data->blb;
      int fngfs = sPp->data->fngfs;
      while(BP)
      {
        Block *cg=BP->data;
        if(myrank == cg->rank) 
        {
	   if(RK_count==0)
	   {
		if(f_NullEvol_pg0A(cg->shape,cg->X[0],cg->X[1],cg->X[2],
			  	   cg->fgfs[g220->sgfn],cg->fgfs[g230->sgfn],cg->fgfs[g330->sgfn],
			  	   cg->fgfs[g01->sgfn],
			  	   cg->fgfs[p02->sgfn],cg->fgfs[p03->sgfn],
			  	   cg->fgfs[g02->sgfn],cg->fgfs[g03->sgfn],Rmin))	    
	       	{
			cout<<"find NaN of pg0A in NullShell domain: sst = "<<sPp->data->sst<<", ("<<cg->bbox[0]<<":"<<cg->bbox[3]<<","
			<<cg->bbox[1]<<":"<<cg->bbox[4]<<","<<cg->bbox[2]<<":"<<cg->bbox[5]<<")"<<endl;
			ERROR = 1;
	        }
	   }
	   else
	   {
		if(f_NullEvol_pg0A(cg->shape,cg->X[0],cg->X[1],cg->X[2],
			  	   cg->fgfs[g22->sgfn],cg->fgfs[g23->sgfn],cg->fgfs[g33->sgfn],
			  	   cg->fgfs[g01->sgfn],
			  	   cg->fgfs[p02->sgfn],cg->fgfs[p03->sgfn],
			  	   cg->fgfs[g02->sgfn],cg->fgfs[g03->sgfn],Rmin))		    
	       	{
			cout<<"find NaN of pg0A in NullShell domain: sst = "<<sPp->data->sst<<", ("<<cg->bbox[0]<<":"<<cg->bbox[3]<<","
			<<cg->bbox[1]<<":"<<cg->bbox[4]<<","<<cg->bbox[2]<<":"<<cg->bbox[5]<<")"<<endl;
			ERROR = 1;
	        }
	   }
	}
        if(BP==sPp->data->ble) break;
        BP=BP->next;
      }
      sPp=sPp->next;
    }
    Synch(pg0AList,Symmetry,pg0Awt,2);
    if(RK_count==3) Dump_Data(pg0AList,0,PhysTime,dT);
//check error information
  {int erh=ERROR;MPI_Allreduce(&erh,&ERROR,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD); }
  if(ERROR)
  {
         Dump_Data(g01List,0,PhysTime,dT);
         if(myrank == 0)
	 {
            if(ErrorMonitor && ErrorMonitor->outfile) 
	       ErrorMonitor->outfile<<"find NaN in beta on NullShell Patches at t = "<<PhysTime<<endl;
	    else
	                        cout<<"find NaN in beta on NullShell Patches at t = "<<PhysTime<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	 }
  }
#if 0
// for gauge variable g00
  {
    sPp=PatL;
    while(sPp)
    {
      MyList<Block> *BP=sPp->data->blb;
      int fngfs = sPp->data->fngfs;
      while(BP)
      {
        Block *cg=BP->data;
        if(myrank == cg->rank) 
        {
	   if(RK_count==0)
	   {
		f_get_gauge_g00_real(cg->shape,cg->X[0],cg->X[1],cg->X[2],
			  	   cg->fgfs[g220->sgfn],cg->fgfs[g230->sgfn],cg->fgfs[g330->sgfn],
                                   cg->fgfs[Theta22->sgfn],cg->fgfs[Theta23->sgfn],cg->fgfs[Theta33->sgfn],
			  	   cg->fgfs[g00->sgfn],Rmin);
	   }
	   else
	   {
	       f_get_gauge_g00_real(cg->shape,cg->X[0],cg->X[1],cg->X[2],
			  	   cg->fgfs[g22->sgfn],cg->fgfs[g23->sgfn],cg->fgfs[g33->sgfn],
                                   cg->fgfs[Theta22->sgfn],cg->fgfs[Theta23->sgfn],cg->fgfs[Theta33->sgfn],
			  	   cg->fgfs[g00->sgfn],Rmin);
	   }
	}
        if(BP==sPp->data->ble) break;
        BP=BP->next;
      }
      sPp=sPp->next;
    }
    Synch(g00List,Symmetry,g00wt,1);
    if(RK_count==3) Dump_Data(g00List,0,PhysTime,dT); 
  }
// evolve ThetaAB
    sPp=PatL;
    while(sPp)
    {
      MyList<Block> *BP=sPp->data->blb;
      int fngfs = sPp->data->fngfs;
      while(BP)
      {
        Block *cg=BP->data;
        if(myrank == cg->rank) 
        {
	   if(RK_count==0)
	   {
		if(f_NullEvol_Theta2(cg->shape,cg->X[0],cg->X[1],cg->X[2],
			  	   cg->fgfs[g220->sgfn],cg->fgfs[g230->sgfn],cg->fgfs[g330->sgfn],
			  	   cg->fgfs[g00->sgfn],cg->fgfs[g01->sgfn],
			  	   cg->fgfs[g02->sgfn],cg->fgfs[g03->sgfn],
			  	   cg->fgfs[p02->sgfn],cg->fgfs[p03->sgfn],
			  	   cg->fgfs[Theta22->sgfn],cg->fgfs[Theta23->sgfn],cg->fgfs[Theta33->sgfn],Rmin))		    
	       	{
			cout<<"find NaN of ThetaAB in NullShell domain: sst = "<<sPp->data->sst<<", ("<<cg->bbox[0]<<":"<<cg->bbox[3]<<","
			<<cg->bbox[1]<<":"<<cg->bbox[4]<<","<<cg->bbox[2]<<":"<<cg->bbox[5]<<")"<<endl;
			ERROR = 1;
	        }
	   }
	   else
	   {
		if(f_NullEvol_Theta2(cg->shape,cg->X[0],cg->X[1],cg->X[2],
			  	   cg->fgfs[g22->sgfn],cg->fgfs[g23->sgfn],cg->fgfs[g33->sgfn],
			  	   cg->fgfs[g00->sgfn],cg->fgfs[g01->sgfn],
			  	   cg->fgfs[g02->sgfn],cg->fgfs[g03->sgfn],
			  	   cg->fgfs[p02->sgfn],cg->fgfs[p03->sgfn],
			  	   cg->fgfs[Theta22->sgfn],cg->fgfs[Theta23->sgfn],cg->fgfs[Theta33->sgfn],Rmin))		    
	       	{
			cout<<"find NaN of ThetaAB in NullShell domain: sst = "<<sPp->data->sst<<", ("<<cg->bbox[0]<<":"<<cg->bbox[3]<<","
			<<cg->bbox[1]<<":"<<cg->bbox[4]<<","<<cg->bbox[2]<<":"<<cg->bbox[5]<<")"<<endl;
			ERROR = 1;
	        }
	   }
	}
        if(BP==sPp->data->ble) break;
        BP=BP->next;
      }
      sPp=sPp->next;
    }
    Synch(ThetaList,Symmetry,Thetawt,3);
    if(RK_count==3) Dump_Data(ThetaList,0,PhysTime,dT);
//check error information
  {int erh=ERROR;MPI_Allreduce(&erh,&ERROR,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD); }
  if(ERROR)
  {
         Dump_Data(pg0AList,0,PhysTime,dT);
         if(myrank == 0)
	 {
            if(ErrorMonitor && ErrorMonitor->outfile) 
	       ErrorMonitor->outfile<<"find NaN in beta on NullShell Patches at t = "<<PhysTime<<endl;
	    else
	                        cout<<"find NaN in beta on NullShell Patches at t = "<<PhysTime<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	 }
  }
#elif 1
// evolve ThetaAB and g00
    sPp=PatL;
    while(sPp)
    {
      MyList<Block> *BP=sPp->data->blb;
      int fngfs = sPp->data->fngfs;
      while(BP)
      {
        Block *cg=BP->data;
        if(myrank == cg->rank) 
        {
	   if(RK_count==0)
	   {
		if(f_NullEvol_Thetag00(cg->shape,cg->X[0],cg->X[1],cg->X[2],
			  	   cg->fgfs[g220->sgfn],cg->fgfs[g230->sgfn],cg->fgfs[g330->sgfn],
			  	   cg->fgfs[g00->sgfn],cg->fgfs[g01->sgfn],
			  	   cg->fgfs[g02->sgfn],cg->fgfs[g03->sgfn],
			  	   cg->fgfs[p02->sgfn],cg->fgfs[p03->sgfn],
			  	   cg->fgfs[Theta22->sgfn],cg->fgfs[Theta23->sgfn],cg->fgfs[Theta33->sgfn],Rmin))		    
	       	{
			cout<<"find NaN of ThetaAB in NullShell domain: sst = "<<sPp->data->sst<<", ("<<cg->bbox[0]<<":"<<cg->bbox[3]<<","
			<<cg->bbox[1]<<":"<<cg->bbox[4]<<","<<cg->bbox[2]<<":"<<cg->bbox[5]<<")"<<endl;
			ERROR = 1;
	        }
	   }
	   else
	   {
		if(f_NullEvol_Thetag00(cg->shape,cg->X[0],cg->X[1],cg->X[2],
			  	   cg->fgfs[g22->sgfn],cg->fgfs[g23->sgfn],cg->fgfs[g33->sgfn],
			  	   cg->fgfs[g00->sgfn],cg->fgfs[g01->sgfn],
			  	   cg->fgfs[g02->sgfn],cg->fgfs[g03->sgfn],
			  	   cg->fgfs[p02->sgfn],cg->fgfs[p03->sgfn],
			  	   cg->fgfs[Theta22->sgfn],cg->fgfs[Theta23->sgfn],cg->fgfs[Theta33->sgfn],Rmin))		    
	       	{
			cout<<"find NaN of ThetaAB in NullShell domain: sst = "<<sPp->data->sst<<", ("<<cg->bbox[0]<<":"<<cg->bbox[3]<<","
			<<cg->bbox[1]<<":"<<cg->bbox[4]<<","<<cg->bbox[2]<<":"<<cg->bbox[5]<<")"<<endl;
			ERROR = 1;
	        }
	   }
	}
        if(BP==sPp->data->ble) break;
        BP=BP->next;
      }
      sPp=sPp->next;
    }
    Synch(ThetaList,Symmetry,Thetawt,3);
    if(RK_count==3) Dump_Data(ThetaList,0,PhysTime,dT);
    Synch(g00List,Symmetry,g00wt,1);
    if(RK_count==3) Dump_Data(g00List,0,PhysTime,dT); 
//check error information
  {int erh=ERROR;MPI_Allreduce(&erh,&ERROR,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD); }
  if(ERROR)
  {
         Dump_Data(pg0AList,0,PhysTime,dT);
         if(myrank == 0)
	 {
            if(ErrorMonitor && ErrorMonitor->outfile) 
	       ErrorMonitor->outfile<<"find NaN in beta on NullShell Patches at t = "<<PhysTime<<endl;
	    else
	                        cout<<"find NaN in beta on NullShell Patches at t = "<<PhysTime<<endl;
	    MPI_Abort(MPI_COMM_WORLD,1);
	 }
  }
#endif
#endif
}
#else
void NullShellPatch2::HyperSlice(double dT, double PhysTime, monitor *ErrorMonitor, int RK_count)
{
  int ERROR = 0;
  Null_Boundary(PhysTime);

  MyList<ss_patch> *sPp;

  // evolve g01
  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        if (RK_count == 0)
        {
          if (f_NullEvol_g01(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                             cg->fgfs[g220->sgfn], cg->fgfs[g230->sgfn], cg->fgfs[g330->sgfn],
                             cg->fgfs[g01->sgfn], Rmin))
          {
            cout << "find NaN of g01 in NullShell domain: sst = " << sPp->data->sst << ", (" << cg->bbox[0] << ":" << cg->bbox[3] << ","
                 << cg->bbox[1] << ":" << cg->bbox[4] << "," << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
            ERROR = 1;
          }
        }
        else
        {
          if (f_NullEvol_g01(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                             cg->fgfs[g22->sgfn], cg->fgfs[g23->sgfn], cg->fgfs[g33->sgfn],
                             cg->fgfs[g01->sgfn], Rmin))
          {
            cout << "find NaN of g01 in NullShell domain: sst = " << sPp->data->sst << ", (" << cg->bbox[0] << ":" << cg->bbox[3] << ","
                 << cg->bbox[1] << ":" << cg->bbox[4] << "," << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
            ERROR = 1;
          }
        }
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }
  Synch(g01List, Symmetry, g01wt, 1);
  //    if(RK_count==3) Dump_Data(g01List,0,PhysTime,dT);
  // check error information
  {
    int erh = ERROR;
    MPI_Allreduce(&erh, &ERROR, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  }
  if (ERROR)
  {
    if (RK_count == 0)
      Dump_Data(StateList, 0, PhysTime, dT);
    else
      Dump_Data(SynchList_pre, 0, PhysTime, dT);
    if (myrank == 0)
    {
      if (ErrorMonitor && ErrorMonitor->outfile)
        ErrorMonitor->outfile << "find NaN in beta on NullShell Patches at t = " << PhysTime << endl;
      else
        cout << "find NaN in beta on NullShell Patches at t = " << PhysTime << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  // evolve p02, p03, g02 and g03
  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        if (RK_count == 0)
        {
          if (f_NullEvol_pg0A(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                              cg->fgfs[g220->sgfn], cg->fgfs[g230->sgfn], cg->fgfs[g330->sgfn],
                              cg->fgfs[g01->sgfn],
                              cg->fgfs[p02->sgfn], cg->fgfs[p03->sgfn],
                              cg->fgfs[g02->sgfn], cg->fgfs[g03->sgfn], Rmin))
          {
            cout << "find NaN of pg0A in NullShell domain: sst = " << sPp->data->sst << ", (" << cg->bbox[0] << ":" << cg->bbox[3] << ","
                 << cg->bbox[1] << ":" << cg->bbox[4] << "," << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
            ERROR = 1;
          }
        }
        else
        {
          if (f_NullEvol_pg0A(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                              cg->fgfs[g22->sgfn], cg->fgfs[g23->sgfn], cg->fgfs[g33->sgfn],
                              cg->fgfs[g01->sgfn],
                              cg->fgfs[p02->sgfn], cg->fgfs[p03->sgfn],
                              cg->fgfs[g02->sgfn], cg->fgfs[g03->sgfn], Rmin))
          {
            cout << "find NaN of pg0A in NullShell domain: sst = " << sPp->data->sst << ", (" << cg->bbox[0] << ":" << cg->bbox[3] << ","
                 << cg->bbox[1] << ":" << cg->bbox[4] << "," << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
            ERROR = 1;
          }
        }
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }
  Synch(pg0AList, Symmetry, pg0Awt, 2);
  //    if(RK_count==3) Dump_Data(pg0AList,0,PhysTime,dT);
  // check error information
  {
    int erh = ERROR;
    MPI_Allreduce(&erh, &ERROR, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  }
  if (ERROR)
  {
    Dump_Data(g01List, 0, PhysTime, dT);
    if (myrank == 0)
    {
      if (ErrorMonitor && ErrorMonitor->outfile)
        ErrorMonitor->outfile << "find NaN in beta on NullShell Patches at t = " << PhysTime << endl;
      else
        cout << "find NaN in beta on NullShell Patches at t = " << PhysTime << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }
  // for gauge variable g00
  {
    sPp = PatL;
    while (sPp)
    {
      MyList<Block> *BP = sPp->data->blb;
      int fngfs = sPp->data->fngfs;
      while (BP)
      {
        Block *cg = BP->data;
        if (myrank == cg->rank)
        {
          f_get_g00_with_t(PhysTime, cg->shape, cg->X[0], cg->X[1], cg->X[2],
                           cg->fgfs[g00->sgfn], Rmin, sPp->data->sst);
        }
        if (BP == sPp->data->ble)
          break;
        BP = BP->next;
      }
      sPp = sPp->next;
    }
    //    if(RK_count==3) Dump_Data(g00List,0,PhysTime,dT);
  }
  // evolve ThetaAB
  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        if (RK_count == 0)
        {
          if (f_NullEvol_Theta2(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                cg->fgfs[g220->sgfn], cg->fgfs[g230->sgfn], cg->fgfs[g330->sgfn],
                                cg->fgfs[g00->sgfn], cg->fgfs[g01->sgfn],
                                cg->fgfs[g02->sgfn], cg->fgfs[g03->sgfn],
                                cg->fgfs[p02->sgfn], cg->fgfs[p03->sgfn],
                                cg->fgfs[Theta22->sgfn], cg->fgfs[Theta23->sgfn], cg->fgfs[Theta33->sgfn], Rmin))
          {
            cout << "find NaN of ThetaAB in NullShell domain: sst = " << sPp->data->sst << ", (" << cg->bbox[0] << ":" << cg->bbox[3] << ","
                 << cg->bbox[1] << ":" << cg->bbox[4] << "," << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
            ERROR = 1;
          }
        }
        else
        {
          if (f_NullEvol_Theta2(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                cg->fgfs[g22->sgfn], cg->fgfs[g23->sgfn], cg->fgfs[g33->sgfn],
                                cg->fgfs[g00->sgfn], cg->fgfs[g01->sgfn],
                                cg->fgfs[g02->sgfn], cg->fgfs[g03->sgfn],
                                cg->fgfs[p02->sgfn], cg->fgfs[p03->sgfn],
                                cg->fgfs[Theta22->sgfn], cg->fgfs[Theta23->sgfn], cg->fgfs[Theta33->sgfn], Rmin))
          {
            cout << "find NaN of ThetaAB in NullShell domain: sst = " << sPp->data->sst << ", (" << cg->bbox[0] << ":" << cg->bbox[3] << ","
                 << cg->bbox[1] << ":" << cg->bbox[4] << "," << cg->bbox[2] << ":" << cg->bbox[5] << ")" << endl;
            ERROR = 1;
          }
        }
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }
  Synch(ThetaList, Symmetry, Thetawt, 3);
  //    if(RK_count==3) Dump_Data(ThetaList,0,PhysTime,dT);

  // check error information
  {
    int erh = ERROR;
    MPI_Allreduce(&erh, &ERROR, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  }
  if (ERROR)
  {
    Dump_Data(pg0AList, 0, PhysTime, dT);
    if (myrank == 0)
    {
      if (ErrorMonitor && ErrorMonitor->outfile)
        ErrorMonitor->outfile << "find NaN in beta on NullShell Patches at t = " << PhysTime << endl;
      else
        cout << "find NaN in beta on NullShell Patches at t = " << PhysTime << endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }
}
#endif
#if 0
void NullShellPatch2::Null_Boundary(double PhysTime)
{
    MyList<ss_patch> *sPp;

    sPp=PatL;
    while(sPp)
    {
      MyList<Block> *BP=sPp->data->blb;
      int fngfs = sPp->data->fngfs;
      while(BP)
      {
        Block *cg=BP->data;
        if(myrank == cg->rank) 
        {
		f_get_null_boundary2(cg->shape,cg->X[0],cg->X[1],cg->X[2],
			            cg->fgfs[g22->sgfn],cg->fgfs[g23->sgfn],cg->fgfs[g33->sgfn],
				    cg->fgfs[g01->sgfn],
				    cg->fgfs[p02->sgfn],cg->fgfs[p03->sgfn],
				    cg->fgfs[g02->sgfn],cg->fgfs[g03->sgfn],
		   	            cg->fgfs[Theta22->sgfn],cg->fgfs[Theta23->sgfn],cg->fgfs[Theta33->sgfn],Rmin);
// for Theta_AB	  
	  f_get_gauge_g00_K(cg->shape,cg->X[0],cg->X[1],cg->X[2],
			  	   cg->fgfs[g220->sgfn],cg->fgfs[g230->sgfn],cg->fgfs[g330->sgfn],
                                   cg->fgfs[Theta22->sgfn],cg->fgfs[Theta23->sgfn],cg->fgfs[Theta33->sgfn],
			  	   cg->fgfs[g00->sgfn],Rmin);
	}
        if(BP==sPp->data->ble) break;
        BP=BP->next;
      }
      sPp=sPp->next;
    }
// boundary for Theta_AB   
//Synchronize K
    Synch(g00List,Symmetry,g00wt,1);
    sPp=PatL;
    int IZEO=1;
    while(sPp)
    {
      MyList<Block> *BP=sPp->data->blb;
      int fngfs = sPp->data->fngfs;
      while(BP)
      {
        Block *cg=BP->data;
        if(myrank == cg->rank) 
        {
	   f_get_gauge_g00(cg->shape,cg->X[0],cg->X[1],cg->X[2],
			  	   cg->fgfs[g220->sgfn],cg->fgfs[g230->sgfn],cg->fgfs[g330->sgfn],
                                   cg->fgfs[Theta22->sgfn],cg->fgfs[Theta23->sgfn],cg->fgfs[Theta33->sgfn],
			  	   cg->fgfs[g00->sgfn],Rmin,IZEO);
	}
        if(BP==sPp->data->ble) break;
        BP=BP->next;
      }
      sPp=sPp->next;
    }
    Synch(ThetaList,Symmetry,Thetawt,3); 
    //Synch(ThetaList,Symmetry,g00wt,1); 
// boundary condition is independent of angular direction, do not need synch
//    Synch(pg0AList,Symmetry,pg0Awt,2,-1);
//    Synch(g00List,Symmetry,g00wt,1,-1); 
//    Synch(ThetaList,Symmetry,Thetawt,3,-1); 
}
#else
void NullShellPatch2::Null_Boundary(double PhysTime)
{
  MyList<ss_patch> *sPp;

  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        f_get_null_boundary3(PhysTime, cg->shape, cg->X[0], cg->X[1], cg->X[2],
                             cg->fgfs[g22->sgfn], cg->fgfs[g23->sgfn], cg->fgfs[g33->sgfn],
                             cg->fgfs[g01->sgfn],
                             cg->fgfs[p02->sgfn], cg->fgfs[p03->sgfn],
                             cg->fgfs[g02->sgfn], cg->fgfs[g03->sgfn],
                             cg->fgfs[Theta22->sgfn], cg->fgfs[Theta23->sgfn], cg->fgfs[Theta33->sgfn], Rmin, sPp->data->sst);
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }
  /*
  // check Synch
      Synch(g01List,Symmetry,g01wt,1);
      Dump_Data(g01List,0,PhysTime,1);
      Synch(pg0AList,Symmetry,pg0Awt,2);
      Dump_Data(pg0AList,0,PhysTime,1);
      Synch(StateList,Symmetry,Thetawt,3);
      Dump_Data(StateList,0,PhysTime,1);
      Synch(ThetaList,Symmetry,Thetawt,3);
      Dump_Data(ThetaList,0,PhysTime,1);
      if(myrank==0) MPI_Abort(MPI_COMM_WORLD,1);
  */
}
// 0: real L2 norm; 1: root mean squar
#define L2m 0
double NullShellPatch2::Error_Check(double PhysTime)
{
  MyList<ss_patch> *sPp;

  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        f_get_null_boundary3(PhysTime, cg->shape, cg->X[0], cg->X[1], cg->X[2],
                             cg->fgfs[g221->sgfn], cg->fgfs[g231->sgfn], cg->fgfs[g331->sgfn],
                             cg->fgfs[g01->sgfn],
                             cg->fgfs[p02->sgfn], cg->fgfs[p03->sgfn],
                             cg->fgfs[g22_rhs->sgfn], cg->fgfs[g03->sgfn],
                             cg->fgfs[Theta22->sgfn], cg->fgfs[Theta23->sgfn], cg->fgfs[Theta33->sgfn], Rmin, sPp->data->sst);
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }

  double tvf, dtvf = 0;
  int tN, dtN = 0;
  int BDW = ghost_width, OBDW = overghost;
  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        f_array_subtract(cg->shape, cg->fgfs[g22_rhs->sgfn], cg->fgfs[g02->sgfn]);
#if (L2m == 0)
        f_l2normhelper_sh(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                          sPp->data->bbox[0], sPp->data->bbox[1], sPp->data->bbox[2], sPp->data->bbox[3], sPp->data->bbox[4], sPp->data->bbox[5],
                          cg->fgfs[g22_rhs->sgfn], tvf, BDW, OBDW, Symmetry);
#elif (L2m == 1)
        f_l2normhelper_sh_rms(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                              sPp->data->bbox[0], sPp->data->bbox[1], sPp->data->bbox[2], sPp->data->bbox[3], sPp->data->bbox[4], sPp->data->bbox[5],
                              cg->fgfs[g22_rhs->sgfn], tvf, BDW, OBDW, Symmetry, dtN);
        dtN += dtN;
#endif

        dtvf += tvf;
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }

  //    Dump_Data(RHSList,0,PhysTime,1);
  //    Dump_Data(ThetaList,0,PhysTime,1);
  //    if(myrank==0) MPI_Abort(MPI_COMM_WORLD,1);

  MPI_Allreduce(&dtvf, &tvf, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#if (L2m == 0)
  tvf = sqrt(tvf);
#elif (L2m == 1)
  MPI_Allreduce(&dtN, &tN, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  tvf = sqrt(tvf / tN);
#endif

  return tvf;
}
#undef L2m
#endif

void NullShellPatch2::Compute_News(double PhysTime)
{
  MyList<ss_patch> *sPp;

// get omega and dtomega
// for RT
#if 0 
    sPp=PatL;
    while(sPp)
    {
      MyList<Block> *BP=sPp->data->blb;
      int fngfs = sPp->data->fngfs;
      while(BP)
      {
        Block *cg=BP->data;
        if(myrank == cg->rank) 
        {
	   f_get_omega_and_dtomega_pre(cg->shape,cg->X[0],cg->X[1],cg->X[2],
				   cg->fgfs[g220->sgfn],cg->fgfs[g230->sgfn],cg->fgfs[g330->sgfn],
			  	   cg->fgfs[omega->sgfn],cg->fgfs[dtomega->sgfn],Rmin);
	} 
        if(BP==sPp->data->ble) break;
        BP=BP->next;
      }
      sPp=sPp->next;
    }
// Synch    
    {
    MyList<var> * DG_List;
    DG_List=new MyList<var>(omega);
    Synch(DG_List,Symmetry,g00wt,1);
    DG_List->clearList(); 
    DG_List=new MyList<var>(dtomega);
    Synch(DG_List,Symmetry,g00wt,1);
    DG_List->clearList(); 
    }
// get dtomega    
    sPp=PatL;
    while(sPp)
    {
      MyList<Block> *BP=sPp->data->blb;
      int fngfs = sPp->data->fngfs;
      while(BP)
      {
        Block *cg=BP->data;
        if(myrank == cg->rank) 
        {
	   f_get_dtomega(cg->shape,cg->X[0],cg->X[1],cg->X[2],
				   cg->fgfs[g220->sgfn],cg->fgfs[g230->sgfn],cg->fgfs[g330->sgfn],
			  	   cg->fgfs[omega->sgfn],cg->fgfs[dtomega->sgfn],Rmin);
	} 
        if(BP==sPp->data->ble) break;
        BP=BP->next;
      }
      sPp=sPp->next;
    }
// Synch    
    {
    MyList<var> * DG_List;
    DG_List=new MyList<var>(dtomega);
    Synch(DG_List,Symmetry,g00wt,1);
    DG_List->clearList(); 
    }
#else
  // for linear wave
  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        f_get_omega_and_dtomega_LN(PhysTime, cg->shape, cg->X[0], cg->X[1], cg->X[2],
                                   cg->fgfs[omega->sgfn], cg->fgfs[dtomega->sgfn], Rmin, sPp->data->sst);
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }
  // Synch
  {
    MyList<var> *DG_List;
    DG_List = new MyList<var>(omega);
    Synch(DG_List, Symmetry, g00wt, 1);
    DG_List->clearList();
    DG_List = new MyList<var>(dtomega);
    Synch(DG_List, Symmetry, g00wt, 1);
    DG_List->clearList();
  }
#endif
  // calculate News
  sPp = PatL;
  while (sPp)
  {
    MyList<Block> *BP = sPp->data->blb;
    int fngfs = sPp->data->fngfs;
    while (BP)
    {
      Block *cg = BP->data;
      if (myrank == cg->rank)
      {
        f_get_null_news2(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                         cg->fgfs[omega->sgfn], cg->fgfs[dtomega->sgfn],
                         cg->fgfs[g00->sgfn], cg->fgfs[g01->sgfn],
                         cg->fgfs[g02->sgfn], cg->fgfs[g03->sgfn],
                         cg->fgfs[g220->sgfn], cg->fgfs[g230->sgfn], cg->fgfs[g330->sgfn],
                         cg->fgfs[Theta22->sgfn], cg->fgfs[Theta23->sgfn], cg->fgfs[Theta33->sgfn],
                         cg->fgfs[RNews->sgfn], cg->fgfs[INews->sgfn], Rmin, sPp->data->sst);
      }
      if (BP == sPp->data->ble)
        break;
      BP = BP->next;
    }
    sPp = sPp->next;
  }
}
