#ifndef MYGLOBAL_H
#define MYGLOBAL_H

#include "var.h"
#include "MyList.h"

#ifdef USE_GPU
#include "bssn_gpu_class.h"
#else
#include "bssn_class.h"
#endif

#include "driver.h"

namespace AHFinderDirect
{

	int globalInterpGFL(double *X, double *Y, double *Z, int Ns,
						double *Data);

	int globalInterpGFLlash(double *X, double *Y, double *Z, int Ns,
							double *Data);

	void AHFinderDirect_setup(MyList<var> *AHList, MyList<var> *GaugeList, bssn_class *ADM,
							  int Symmetry, int HN, double *PhysTime);

	void AHFinderDirect_cleanup();

	void AHFinderDirect_find_horizons(int HN, int *dumpid,
									  double *xc, double *yc, double *zc, double *xr, double *yr, double *zr,
									  bool *trigger, double *);

	void AHFinderDirect_enforcefind(int HN,
									double *xc, double *yc, double *zc, double *xr, double *yr, double *zr);
	//
	struct state
	{
		int N_procs; // total number of processors
		int my_proc; // processor number of this processor
					 // (0 to N_procs-1)

		int Symmetry;
		double *PhysTime;

		MyList<var> *AHList;
		MyList<var> *GaugeList;

		bssn_class *ADM;

		int N_horizons; // total number of genuine horizons
						// being searched for
		int N_active_procs; // total number of active processors
							// (the active processors are processor
							//  numbers 0 to N_active_procs-1)

		struct iteration_status_buffers isb;

		horizon_sequence *my_hs;

		struct AH_data **AH_data_array;

		double *Data, *oX, *oY, *oZ;
	};
}
#endif /* MYGLOBAL_H */
