
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

#ifdef USE_GPU
#include <cuda_runtime.h>
#include "MPatch.h"
#include "MyList.h"
#include "var.h"

void gpu_lowerboundset_launch(
    cudaStream_t &stream,
    int ex[3], int pex[3],
    double* d_chi0, double TINNY
);

void gpu_pack_launch(
    cudaStream_t stream, const double* d_src_3d, double* d_dst_1d,
    int src_pnx, int src_pny, int dst_nx, int dst_ny, int dst_nz,
    int off_x, int off_y, int off_z
);

void gpu_unpack_launch(
    cudaStream_t stream, const double* d_src_1d, double* d_dst_3d,
    int dst_pnx, int dst_pny, int src_nx, int src_ny, int src_nz,
    int off_x, int off_y, int off_z
);

void gpu_set_value_launch(cudaStream_t stream, const int ex[3], const int pex[3], double* d_data, double value);
void gpu_add_value_launch(cudaStream_t stream, const int ex[3], const int pex[3], double* d_data, double value);
void gpu_array_copy_launch(cudaStream_t stream, const int ex[3], const int pex[3], double* d_dst, const double* d_src);
void gpu_array_add_launch(cudaStream_t stream, const int ex[3], const int pex[3], double* d_dst, const double* d_src);
void gpu_array_subtract_launch(cudaStream_t stream, const int ex[3], const int pex[3], double* d_dst, const double* d_src);

void gpu_scan_nan_launch(
    cudaStream_t stream, const int ex[3], const int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_f, int* d_count, double* d_first_xyzf
);

void gpu_pointcopy_cell_launch(
    cudaStream_t stream,
    const int ext_out[3], const int pext_out[3],
    double* d_data_out,
    const double llbout[3], const double uubout[3],
    double xx, double yy, double zz, double dv
);

void gpu_copy_cell_launch(
    cudaStream_t stream,
    const int ext_out[3], const int pext_out[3], double* d_data_out,
    const double llbout[3], const double uubout[3],
    const int ext_in[3], const int pext_in[3], const double* d_data_in,
    const double llbin[3], const double uubin[3],
    const double lcopy[3], const double ucopy[3]
);

void gpu_average_launch(cudaStream_t stream, const int ex[3], const int pex[3], const double* d_f1, const double* d_f2, double* d_fout);
void gpu_average3_launch(cudaStream_t stream, const int ex[3], const int pex[3], const double* d_f1, const double* d_f2, double* d_fout);
void gpu_average2_launch(cudaStream_t stream, const int ex[3], const int pex[3], const double* d_f1, const double* d_f2, const double* d_f3, double* d_fout);
void gpu_average2p_launch(cudaStream_t stream, const int ex[3], const int pex[3], const double* d_f1, const double* d_f2, const double* d_f3, double* d_fout);
void gpu_average2m_launch(cudaStream_t stream, const int ex[3], const int pex[3], const double* d_f1, const double* d_f2, const double* d_f3, double* d_fout);

void gpu_interp_2_launch(
    cudaStream_t stream,
    const int ex[3], const int pex[3],
    const double* d_f, double* d_f_int,
    int il, int jl, int kl,
    double Dx, double Dy, double Dz,
    int ordn, const double SoA[3]
);

void gpu_polint_launch(cudaStream_t stream, const double* d_xa, const double* d_ya, double x, double* d_y, double* d_dy, int ordn);

void gpu_polin2_launch(
    cudaStream_t stream,
    const double* d_x1a, const double* d_x2a, const double* d_ya,
    double x1, double x2, double* d_y, double* d_dy, int ordn
);

void gpu_polin3_launch(
    cudaStream_t stream,
    const double* d_x1a, const double* d_x2a, const double* d_x3a, const double* d_ya,
    double x1, double x2, double x3, double* d_y, double* d_dy, int ordn
);

void gpu_d2dump_launch(
    cudaStream_t stream, const int ex[3], const int pex[3],
    const double* d_data_in, double* d_data_out,
    const double llb[3], const double uub[3],
    int gord, const double SoA[3]
);

void gpu_global_interp_launch(
    cudaStream_t stream,
    int NN, int DIM,
    double* d_XX_0, double* d_XX_1, double* d_XX_2,
    int shape_0, int shape_1, int shape_2,
    int pshape_0, int pshape_1, int pshape_2,
    double* d_X_0, double* d_X_1, double* d_X_2,
    double* d_field,
    double llb_0, double llb_1, double llb_2,
    double uub_0, double uub_1, double uub_2,
    double DH_0, double DH_1, double DH_2,
    int ordn, double SoA_0, double SoA_1, double SoA_2,
    int Symmetry, int var_idx, int num_var,
    double* d_shellf, int* d_weight
);

void gpu_global_interp_scalar_launch(
    cudaStream_t stream,
    const int ex[3], const int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_f, double* d_f_int,
    double x1, double y1, double z1,
    int ordn, const double SoA[3], int symmetry
);

void gpu_global_interp_ss_launch(
    cudaStream_t stream,
    const int ex[3], const int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_f, double* d_f_int,
    double x1, double y1, double z1,
    int ordn, const double SoA[3], int sst
);

void gpu_global_interp_ss_2d_launch(
    cudaStream_t stream,
    const int ex[3], const int pex[3],
    const double* d_X, const double* d_Y,
    const double* d_f, double* d_f_int,
    double x1, double y1, int indZ,
    int ordn, const double SoA[3], int sst
);

void gpu_l2normhelper_launch(
    cudaStream_t stream,
    const int* ex, const int* pex,
    const double* X, const double* Y, const double* Z,
    double xmin, double ymin, double zmin,
    double xmax, double ymax, double zmax,
    const double* d_f, double* d_sum, int gw
);

void gpu_l2normhelper_sh_launch(
    cudaStream_t stream,
    const int* ex, const int* pex,
    const double* X, const double* Y, const double* Z,
    double xmin, double ymin, double zmin,
    double xmax, double ymax, double zmax,
    const double* d_f, double* d_sum, int gw, int ogw, int Symmetry
);

void gpu_l2normhelper_sh_rms_launch(
    cudaStream_t stream,
    const int* ex, const int* pex,
    const double* X, const double* Y, const double* Z,
    double xmin, double ymin, double zmin,
    double xmax, double ymax, double zmax,
    const double* d_f, double* d_sum, int gw, int ogw, int Symmetry, int* d_count
);

void gpu_global_interpind_launch(
    cudaStream_t stream,
    const int ex[3], const int pex[3],
    const double* d_f, double* d_f_int,
    int ordn, const double SoA[3],
    const int* d_inds, const double* d_coef, int sst
);

void gpu_global_interpind2d_launch(
    cudaStream_t stream,
    const int ex[3], const int pex[3],
    const double* d_f, double* d_f_int,
    int ordn, const double SoA[3],
    const int* d_inds, const double* d_coef, int sst
);

void gpu_global_interpind1d_launch(
    cudaStream_t stream,
    const int ex[3], const int pex[3],
    const double* d_f, double* d_f_int,
    int ordn, const double SoA[3],
    const int* d_inds, const double* d_coef, int sst, int dumyd
);

void gpu_find_maximum_launch(
    cudaStream_t stream,
    const int ex[3], const int pex[3],
    const double* d_X, const double* d_Y, const double* d_Z,
    const double* d_fun,
    double* d_val, double* d_pos,
    const int llb[3], const int uub[3],
    int* d_tmp_idx
);

void gpu_ffact_launch(cudaStream_t stream, int n, double* d_out);
void gpu_wigner_d_launch(cudaStream_t stream, int l, int m, int s, double costheta, double* d_out);
void gpu_yslm_launch(
    cudaStream_t stream, int s, int l, int m, double theta, double phi,
    double* d_real_out, double* d_imag_out
);

void gpu_global_interp_amr_launch(
    cudaStream_t stream,
    int active_count, int DIM,
    int* d_active_indices,
    double* d_XX_0, double* d_XX_1, double* d_XX_2,
    int shape_0, int shape_1, int shape_2,
    int pshape_0, int pshape_1, int pshape_2,
    double* d_X_0, double* d_X_1, double* d_X_2,
    double* d_field,
    int ordn, double SoA_0, double SoA_1, double SoA_2,
    int Symmetry, int var_idx, int num_var,
    double* d_shellf
);

void gpu_apply_symmetry_bd(
    cudaStream_t stream,
    int extc0, int extc1, int extc2,
    int pextc0, int pextc1, int pextc2,
    double* __restrict__ func,
    double SoA0, double SoA1, double SoA2
);

void gpu_apply_symmetry_bd_fused_launch(
    cudaStream_t stream,
    int extc0, int extc1, int extc2,
    int pextc0, int pextc1, int pextc2,
    int sym_type,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Gamx, double* Gamy, double* Gamz, double* Lap,
    double* betax, double* betay, double* betaz,
    double* dtSfx, double* dtSfy, double* dtSfz
);

void gpu_apply_scalar_symmetry_launch(
    cudaStream_t stream,
    int extc0, int extc1, int extc2,
    int pextc0, int pextc1, int pextc2,
    double* Sphi, double* Spi
);

#include "gpu_fusion.h"

void gpu_apply_symmetry_bd_vars_launch(
    cudaStream_t stream,
    int extc0, int extc1, int extc2,
    int pextc0, int pextc1, int pextc2,
    int sym_type,
    FusedVars vars
);

void gpu_fused_pack_launch(
    cudaStream_t stream,
    FusedVars vars, double* d_dst_1d,
    int src_pnx, int src_pny,
    int dst_nx, int dst_ny, int dst_nz,
    int off_x, int off_y, int off_z
);

void gpu_fused_unpack_launch(
    cudaStream_t stream, 
    const double* d_src_1d, FusedVars vars,
    int dst_pnx, int dst_pny,
    int src_nx, int src_ny, int src_nz,
    int off_x, int off_y, int off_z
);

void gpu_fused_time_interp_3_launch(cudaStream_t stream, const int ext[3], const int pext[3], FusedTimeVars3 vars, int tindex);

void gpu_fused_time_interp_4_launch(cudaStream_t stream, const int ext[3], const int pext[3], FusedTimeVars4 vars, int tindex);

void gpu_fused_global_interp_launch(
    cudaStream_t stream,
    int NN, int DIM,
    double* d_XX_0, double* d_XX_1, double* d_XX_2,
    int shape_0, int shape_1, int shape_2,
    int pshape_0, int pshape_1, int pshape_2,
    double* d_X_0, double* d_X_1, double* d_X_2,
    FusedVars vars,
    double llb_0, double llb_1, double llb_2,
    double uub_0, double uub_1, double uub_2,
    double DH_0, double DH_1, double DH_2,
    int ordn, int Symmetry,
    int base_var_idx, int total_num_var,
    double* d_shellf, int* d_weight
);

void gpu_fused_global_interp_amr_launch(
    cudaStream_t stream,
    int active_count, int DIM,
    int* d_active_indices,
    double* d_XX_0, double* d_XX_1, double* d_XX_2,
    int shape_0, int shape_1, int shape_2,
    int pshape_0, int pshape_1, int pshape_2,
    double* d_X_0, double* d_X_1, double* d_X_2,
    FusedVars vars,
    int ordn, int Symmetry,
    int base_var_idx, int total_num_var,
    double* d_shellf
);

void gpu_batch_interp_ss_launch(
    cudaStream_t stream,
    int NN, int DIM, int sst,
    double* d_XX_0, double* d_XX_1, double* d_XX_2,
    int shape_0, int shape_1, int shape_2,
    int pshape_0, int pshape_1, int pshape_2,
    double* d_X_0, double* d_X_1, double* d_X_2,
    double* d_field,
    double llb_0, double llb_1, double llb_2,
    double uub_0, double uub_1, double uub_2,
    double DH_0, double DH_1, double DH_2,
    int ordn, double SoA_0, double SoA_1, double SoA_2,
    int Symmetry, int var_idx, int num_var,
    double* d_shellf, int* d_weight,
    int* d_point_sst
);

#endif /* USE_GPU */

#endif /* FMISC_H */
