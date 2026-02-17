
#ifndef FMISC_H
#define FMISC_H

#ifdef fortran1
#define f_interp_2 interp_2
#define f_pointcopy pointcopy
#define f_copy copy
#define f_global_interp global_interp
#define f_global_interp_ss global_interp_ss
#define f_global_interp_ss_2d global_interp_ss_2d
#define f_global_interpind global_interpind
#define f_global_interpind2d global_interpind2d
#define f_global_interpind1d global_interpind1d
#define f_l2normhelper l2normhelper
#define f_l2normhelper_sh l2normhelper_sh
#define f_l2normhelper_sh_rms l2normhelper_sh_rms
#define f_average average
#define f_average3 average3
#define f_average2 average2
#define f_average2p average2p
#define f_average2m average2m
#define f_lowerboundset lowerboundset
#define f_set_value set_value
#define f_add_value add_value
#define f_array_add array_add
#define f_array_copy array_copy
#define f_array_subtract array_subtract
#define f_fft four1
#define f_find_maximum find_maximum
#define f_polint polint
#define f_d2dump d2dump
#endif
#ifdef fortran2
#define f_interp_2 INTERP_2
#define f_pointcopy POINTCOPY
#define f_copy COPY
#define f_global_interp GLOBAL_INTERP
#define f_global_interp_ss GLOBAL_INTERP_SS
#define f_global_interp_ss_2d GLOBAL_INTERP_SS_2D
#define f_global_interpind GLOBAL_INTERPIND
#define f_global_interpind2d GLOBAL_INTERPIND2D
#define f_global_interpind1d GLOBAL_INTERPIND1D
#define f_l2normhelper L2NORMHELPER
#define f_l2normhelper_sh L2NORMHELPER_SH
#define f_l2normhelper_sh_rms L2NORMHELPER_SH_RMS
#define f_average AVERAGE
#define f_average3 AVERAGE3
#define f_average2 AVERAGE2
#define f_average2p AVERAGE2P
#define f_average2m AVERAGE2M
#define f_lowerboundset LOWERBOUNDSET
#define f_set_value SET_VALU
#define f_add_value ADD_VALUE
#define f_array_add ARRAY_ADD
#define f_array_copy ARRAY_COPY
#define f_array_subtract ARRAY_SUBTRACT
#define f_fft FOUR1
#define f_find_maximum FIND_MAXIMUM
#define f_polint POLINT
#define f_d2dump D2DUMP
#endif
#ifdef fortran3
#define f_interp_2 interp_2_
#define f_pointcopy pointcopy_
#define f_copy copy_
#define f_global_interp global_interp_
#define f_global_interp_ss global_interp_ss_
#define f_global_interp_ss_2d global_interp_ss_2d_
#define f_global_interpind global_interpind_
#define f_global_interpind2d global_interpind2d_
#define f_global_interpind1d global_interpind1d_
#define f_l2normhelper l2normhelper_
#define f_l2normhelper_sh l2normhelper_sh_
#define f_l2normhelper_sh_rms l2normhelper_sh_rms_
#define f_average average_
#define f_average3 average3_
#define f_average2 average2_
#define f_average2p average2p_
#define f_average2m average2m_
#define f_lowerboundset lowerboundset_
#define f_set_value set_value_
#define f_add_value add_value_
#define f_array_add array_add_
#define f_array_copy array_copy_
#define f_array_subtract array_subtract_
#define f_fft four1_
#define f_find_maximum find_maximum_
#define f_polint polint_
#define f_d2dump d2dump_
#endif

extern "C"
{
	void f_pointcopy(int &,
					 double *, double *, int *, double *,
					 double &, double &, double &, double &);
}

extern "C"
{
	void f_copy(int &,
				double *, double *, int *, double *,
				double *, double *, int *, double *,
				double *, double *);
}

extern "C"
{
	void f_global_interp(int *, double *, double *, double *,
						 double *, double &,
						 double &, double &, double &,
						 int &, double *, int &);
}

extern "C"
{
	void f_global_interp_ss(int *, double *, double *, double *,
							double *, double &,
							double &, double &, double &,
							int &, double *, int &, int &);
}

extern "C"
{
	void f_global_interp_ss_2d(int *, double *, double *, int &,
							   double *, double &,
							   double &, double &,
							   int &, double *, int &, int &);
}

extern "C"
{
	void f_global_interpind(int *, double *, double *, double *,
							double *, double &,
							double &, double &, double &,
							int &, double *, int &,
							int *, double *, int &);
}

extern "C"
{
	void f_global_interpind2d(int *, double *, double *, double *,
							  double *, double &,
							  double &, double &, double &,
							  int &, double *, int &,
							  int *, double *, int &);
}

extern "C"
{
	void f_global_interpind1d(int *, double *, double *, double *,
							  double *, double &,
							  double &, double &, double &,
							  int &, double *, int &,
							  int *, double *, int &, int &);
}

extern "C"
{
	void f_l2normhelper(int *, double *, double *, double *,
						double &, double &, double &,
						double &, double &, double &,
						double *, double &, int &);
}

extern "C"
{
	void f_l2normhelper_sh(int *, double *, double *, double *,
						   double &, double &, double &,
						   double &, double &, double &,
						   double *, double &, int &, int &, int &);
}

extern "C"
{
	void f_l2normhelper_sh_rms(int *, double *, double *, double *,
							   double &, double &, double &,
							   double &, double &, double &,
							   double *, double &, int &, int &, int &, int &);
}

extern "C"
{
	void f_average(int *, double *, double *, double *);
}

extern "C"
{
	void f_average3(int *, double *, double *, double *);
}

extern "C"
{
	void f_average2(int *, double *, double *, double *, double *);
}

extern "C"
{
	void f_average2p(int *, double *, double *, double *, double *);
}

extern "C"
{
	void f_average2m(int *, double *, double *, double *, double *);
}

extern "C"
{
	void f_lowerboundset(int *, double *, double &);
}

#if 0
extern "C" { void f_interp_2( int *, double *,
	                    double &, int &,int &,int &,int &,int &,int &,
	                    double &,double &,double &,
			    int &, double *,
			    int &
			    );}
#endif

extern "C"
{
	void f_set_value(int *, double *, double &);
}
extern "C"
{
	void f_add_value(int *, double *, double &);
}
extern "C"
{
	void f_array_add(int *, double *, double *);
}
extern "C"
{
	void f_array_copy(int *, double *, double *);
}
extern "C"
{
	void f_array_subtract(int *, double *, double *);
}

extern "C"
{
	void f_fft(double *, int &, int &);
}

extern "C"
{
	void f_find_maximum(int *,
						double *, double *, double *, double *,
						double &, double *, int *, int *);
}

extern "C"
{
	void f_polint(double *, double *, double &, double &, double &, int &);
}

extern "C"
{
	void f_d2dump(int &, double *, double *, int *, double *, double *, int &, double *);
}
#endif /* FMISC_H */
