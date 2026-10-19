#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <sys/time.h>

#ifdef RESULT_CHECK
#include <fstream>
#endif

#include "macrodef.h"
#include "fmisc.h"
#include "bssn_gpu_class.h"
#include "Z4c_gpu_class.h"
#include "bssn_rhs.h"
#include "bssn_rhs_shell_gpu.cuh"
#include "ahf_derivatives_gpu.h"
#include "ahf_gpu.h"
#ifdef With_AHF
#include "myglobal.h"
#endif
#include "enforce_algebra.h"
#include "rungekutta4_rout.h"
#include "sommerfeld_rout.h"

#include "gpu_runtime_coordinator.h"
#include "data_motion_helper.h"

namespace {

void abort_nan_black_hole_predictor(bssn_class& solver, int lev, int black_hole_index, double dT_lev, double **target_positions) {
	if (solver.ErrorMonitor->outfile)
		solver.ErrorMonitor->outfile << "predictor step finds NaN for BH's position from ("
			<< solver.Porg0[black_hole_index][0] << "," << solver.Porg0[black_hole_index][1] << "," << solver.Porg0[black_hole_index][2] << ")" << endl;
	cout << "predictor step finds NaN for BH's position from ("
		 << solver.Porg0[black_hole_index][0] << "," << solver.Porg0[black_hole_index][1] << "," << solver.Porg0[black_hole_index][2] << ")" << endl;
	cout << "to ("
		 << target_positions[black_hole_index][0] << "," << target_positions[black_hole_index][1] << "," << target_positions[black_hole_index][2] << ")" << endl;
	cout << "with velocity ("
		 << solver.Porg_rhs[black_hole_index][0] << "," << solver.Porg_rhs[black_hole_index][1] << "," << solver.Porg_rhs[black_hole_index][2] << ")" << endl;
	MPI_Abort(MPI_COMM_WORLD, 1);
	MyList<var> *dg_list = new MyList<var>(solver.Sfx0);
	dg_list->insert(solver.Sfx0);
	dg_list->insert(solver.Sfy0);
	dg_list->insert(solver.Sfz0);
	Parallel::Dump_Data(solver.GH->PatL[lev], dg_list, 0, solver.PhysTime, dT_lev);
	dg_list->clearList();
}

void apply_black_hole_position_symmetry(bssn_class& solver, double **positions, int black_hole_index) {
	if (solver.Symmetry > 0)
		positions[black_hole_index][2] = fabs(positions[black_hole_index][2]);
	if (solver.Symmetry == 2) {
		positions[black_hole_index][0] = fabs(positions[black_hole_index][0]);
		positions[black_hole_index][1] = fabs(positions[black_hole_index][1]);
	}
}

bool black_hole_position_is_finite(double **positions, int black_hole_index) {
	return finite(positions[black_hole_index][0]) &&
	       finite(positions[black_hole_index][1]) &&
	       finite(positions[black_hole_index][2]);
}

void advance_black_hole_euler_predictor(bssn_class& solver, int lev, double dT_lev) {
	if (solver.BH_num <= 0 || lev != solver.GH->levels - 1)
		return;

	solver.compute_Porg_rhs(solver.Porg0, solver.Porg_rhs, solver.Sfx0, solver.Sfy0, solver.Sfz0, lev);
	for (int black_hole_index = 0; black_hole_index < solver.BH_num; ++black_hole_index) {
		for (int axis = 0; axis < 3; ++axis)
			solver.Porg1[black_hole_index][axis] = solver.Porg0[black_hole_index][axis] + solver.Porg_rhs[black_hole_index][axis] * dT_lev;
		apply_black_hole_position_symmetry(solver, solver.Porg1, black_hole_index);
		if (!black_hole_position_is_finite(solver.Porg1, black_hole_index))
			abort_nan_black_hole_predictor(solver, lev, black_hole_index, dT_lev, solver.Porg1);
	}
}

void advance_black_hole_rk4_predictor(bssn_class& solver, int lev, double dT_lev) {
	if (solver.BH_num <= 0 || lev != solver.GH->levels - 1)
		return;

	int iter_count = 0;
	solver.compute_Porg_rhs(solver.Porg0, solver.Porg_rhs, solver.Sfx0, solver.Sfy0, solver.Sfz0, lev);
	for (int black_hole_index = 0; black_hole_index < solver.BH_num; ++black_hole_index) {
		for (int axis = 0; axis < 3; ++axis)
			f_rungekutta4_scalar(dT_lev, solver.Porg0[black_hole_index][axis], solver.Porg[black_hole_index][axis], solver.Porg_rhs[black_hole_index][axis], iter_count);
		apply_black_hole_position_symmetry(solver, solver.Porg, black_hole_index);
		if (!black_hole_position_is_finite(solver.Porg, black_hole_index))
			abort_nan_black_hole_predictor(solver, lev, black_hole_index, dT_lev, solver.Porg);
	}
}

void abort_nan_black_hole_corrector(bssn_class& solver, int lev, int black_hole_index, double dT_lev, int iter_count) {
	if (solver.ErrorMonitor->outfile)
		solver.ErrorMonitor->outfile << iter_count << " corrector step finds NaN for BH's position from ("
			<< solver.Porg[black_hole_index][0] << "," << solver.Porg[black_hole_index][1] << "," << solver.Porg[black_hole_index][2] << ")" << endl;
	cout << iter_count << " corrector step finds NaN for BH's position from ("
		 << solver.Porg[black_hole_index][0] << "," << solver.Porg[black_hole_index][1] << "," << solver.Porg[black_hole_index][2] << ")" << endl;
	MPI_Abort(MPI_COMM_WORLD, 1);
	MyList<var> *dg_list = new MyList<var>(solver.Sfx0);
	dg_list->insert(solver.Sfx0);
	dg_list->insert(solver.Sfy0);
	dg_list->insert(solver.Sfz0);
	Parallel::Dump_Data(solver.GH->PatL[lev], dg_list, 0, solver.PhysTime, dT_lev);
	dg_list->clearList();
}

void advance_black_hole_rk4_corrector(bssn_class& solver, int lev, double dT_lev, int iter_count) {
	if (solver.BH_num <= 0 || lev != solver.GH->levels - 1)
		return;

	solver.compute_Porg_rhs(solver.Porg, solver.Porg1, solver.Sfx, solver.Sfy, solver.Sfz, lev);
	for (int black_hole_index = 0; black_hole_index < solver.BH_num; ++black_hole_index) {
		for (int axis = 0; axis < 3; ++axis)
			f_rungekutta4_scalar(dT_lev, solver.Porg0[black_hole_index][axis], solver.Porg1[black_hole_index][axis], solver.Porg_rhs[black_hole_index][axis], iter_count);
		apply_black_hole_position_symmetry(solver, solver.Porg1, black_hole_index);
		if (!black_hole_position_is_finite(solver.Porg1, black_hole_index))
			abort_nan_black_hole_corrector(solver, lev, black_hole_index, dT_lev, iter_count);
	}
}

void prepare_black_hole_corrector_for_next_stage(bssn_class& solver, int lev) {
	if (solver.BH_num <= 0 || lev != solver.GH->levels - 1)
		return;

	for (int black_hole_index = 0; black_hole_index < solver.BH_num; ++black_hole_index) {
		solver.Porg[black_hole_index][0] = solver.Porg1[black_hole_index][0];
		solver.Porg[black_hole_index][1] = solver.Porg1[black_hole_index][1];
		solver.Porg[black_hole_index][2] = solver.Porg1[black_hole_index][2];
	}
}

void trace_block_metric_det_summary(
	bssn_class& solver,
	const char* phase,
	int lev,
	int yn,
	Block* block,
	var* gxx_var,
	var* gxy_var,
	var* gxz_var,
	var* gyy_var,
	var* gyz_var,
	var* gzz_var
);

void apply_stage_updates(
	bssn_class& solver,
	Patch* patch,
	Block* block,
	double dT_lev,
	int lev,
	int iter_count,
	MyList<var> *base_state_list,
	MyList<var> *bam_source_list,
	MyList<var> *bam_target_list,
	MyList<var> *rk_target_list,
	MyList<var> *rhs_list,
	int cor
) {
#ifdef WithShell
	const bool use_bam_boundary = false;
	const bool use_fix_boundary = true;
#else
	const bool use_bam_boundary = (lev == 0);
	const bool use_fix_boundary = (lev > 0);
#endif
	std::vector<double*> base_fields;
	std::vector<double*> rk_fields;
	std::vector<double*> rhs_fields;
	base_fields.reserve(32);
	rk_fields.reserve(32);
	rhs_fields.reserve(32);

	MyList<var> *base_node = base_state_list;
	MyList<var> *bam_source_node = bam_source_list;
	MyList<var> *bam_target_node = bam_target_list;
	MyList<var> *rk_target_node = rk_target_list;
	MyList<var> *rhs_node = rhs_list;
	while (base_node) {
		if (use_bam_boundary) {
			gpu_sommerfeld_routbam_launch(
				block->stream,
				block->shape,
				block->d_X[0], block->d_X[1], block->d_X[2],
				patch->bbox[0], patch->bbox[1], patch->bbox[2], patch->bbox[3], patch->bbox[4], patch->bbox[5],
				block->d_fgfs[bam_target_node->data->sgfn],
				block->d_fgfs[bam_source_node->data->sgfn], base_node->data->propspeed, base_node->data->SoA,
				solver.Symmetry
			);
		}
		base_fields.push_back(block->d_fgfs[base_node->data->sgfn]);
		rk_fields.push_back(block->d_fgfs[rk_target_node->data->sgfn]);
		rhs_fields.push_back(block->d_fgfs[rhs_node->data->sgfn]);

		base_node = base_node->next;
		bam_source_node = bam_source_node->next;
		bam_target_node = bam_target_node->next;
		rk_target_node = rk_target_node->next;
		rhs_node = rhs_node->next;
	}

	if (!base_fields.empty()) {
		gpu_rungekutta4_rout_batch_launch(
			block->stream,
			block->shape,
			dT_lev,
			base_fields,
			rk_fields,
			rhs_fields,
			iter_count
		);
	}

	if (use_fix_boundary) {
		base_node = base_state_list;
		rk_target_node = rk_target_list;
		while (base_node) {
			gpu_sommerfeld_rout_launch(
				block->stream,
				block->shape,
				block->d_X[0], block->d_X[1], block->d_X[2],
				patch->bbox[0], patch->bbox[1], patch->bbox[2], patch->bbox[3], patch->bbox[4], patch->bbox[5],
				dT_lev, block->d_fgfs[solver.phi0->sgfn],
				block->d_fgfs[solver.Lap0->sgfn], block->d_fgfs[base_node->data->sgfn], block->d_fgfs[rk_target_node->data->sgfn], base_node->data->SoA,
				solver.Symmetry, cor
			);

			base_node = base_node->next;
			rk_target_node = rk_target_node->next;
		}
	}
}

void launch_bssn_rhs_for_block(
	bssn_class& solver,
	Block* block,
	double trk4,
	int lev,
	double ndeps,
	int co,
	var* chi_var,
	var* trK_var,
	var* gxx_var,
	var* gxy_var,
	var* gxz_var,
	var* gyy_var,
	var* gyz_var,
	var* gzz_var,
	var* Axx_var,
	var* Axy_var,
	var* Axz_var,
	var* Ayy_var,
	var* Ayz_var,
	var* Azz_var,
	var* Gamx_var,
	var* Gamy_var,
	var* Gamz_var,
	var* Lap_var,
	var* betax_var,
	var* betay_var,
	var* betaz_var,
	var* dtSfx_var,
	var* dtSfy_var,
	var* dtSfz_var,
	var* chi_rhs_var,
	var* trK_rhs_var,
	var* gxx_rhs_var,
	var* gxy_rhs_var,
	var* gxz_rhs_var,
	var* gyy_rhs_var,
	var* gyz_rhs_var,
	var* gzz_rhs_var,
	var* Axx_rhs_var,
	var* Axy_rhs_var,
	var* Axz_rhs_var,
	var* Ayy_rhs_var,
	var* Ayz_rhs_var,
	var* Azz_rhs_var,
	var* Gamx_rhs_var,
	var* Gamy_rhs_var,
	var* Gamz_rhs_var,
	var* Lap_rhs_var,
	var* betax_rhs_var,
	var* betay_rhs_var,
	var* betaz_rhs_var,
	var* dtSfx_rhs_var,
	var* dtSfy_rhs_var,
	var* dtSfz_rhs_var
) {
		gpu_compute_rhs_bssn_launch(
			block->stream,
			block->shape, trk4, block->d_X[0], block->d_X[1], block->d_X[2],
			block->d_fgfs[chi_var->sgfn], block->d_fgfs[trK_var->sgfn],
			block->d_fgfs[gxx_var->sgfn], block->d_fgfs[gxy_var->sgfn], block->d_fgfs[gxz_var->sgfn],
			block->d_fgfs[gyy_var->sgfn], block->d_fgfs[gyz_var->sgfn], block->d_fgfs[gzz_var->sgfn],
			block->d_fgfs[Axx_var->sgfn], block->d_fgfs[Axy_var->sgfn], block->d_fgfs[Axz_var->sgfn],
			block->d_fgfs[Ayy_var->sgfn], block->d_fgfs[Ayz_var->sgfn], block->d_fgfs[Azz_var->sgfn],
			block->d_fgfs[Gamx_var->sgfn], block->d_fgfs[Gamy_var->sgfn], block->d_fgfs[Gamz_var->sgfn],
			block->d_fgfs[Lap_var->sgfn],
			block->d_fgfs[betax_var->sgfn], block->d_fgfs[betay_var->sgfn], block->d_fgfs[betaz_var->sgfn],
			block->d_fgfs[dtSfx_var->sgfn], block->d_fgfs[dtSfy_var->sgfn], block->d_fgfs[dtSfz_var->sgfn],
			block->d_fgfs[chi_rhs_var->sgfn], block->d_fgfs[trK_rhs_var->sgfn],
			block->d_fgfs[gxx_rhs_var->sgfn], block->d_fgfs[gxy_rhs_var->sgfn], block->d_fgfs[gxz_rhs_var->sgfn],
			block->d_fgfs[gyy_rhs_var->sgfn], block->d_fgfs[gyz_rhs_var->sgfn], block->d_fgfs[gzz_rhs_var->sgfn],
			block->d_fgfs[Axx_rhs_var->sgfn], block->d_fgfs[Axy_rhs_var->sgfn], block->d_fgfs[Axz_rhs_var->sgfn],
			block->d_fgfs[Ayy_rhs_var->sgfn], block->d_fgfs[Ayz_rhs_var->sgfn], block->d_fgfs[Azz_rhs_var->sgfn],
			block->d_fgfs[Gamx_rhs_var->sgfn], block->d_fgfs[Gamy_rhs_var->sgfn], block->d_fgfs[Gamz_rhs_var->sgfn],
			block->d_fgfs[Lap_rhs_var->sgfn],
			block->d_fgfs[betax_rhs_var->sgfn], block->d_fgfs[betay_rhs_var->sgfn], block->d_fgfs[betaz_rhs_var->sgfn],
			block->d_fgfs[dtSfx_rhs_var->sgfn], block->d_fgfs[dtSfy_rhs_var->sgfn], block->d_fgfs[dtSfz_rhs_var->sgfn],
			block->d_fgfs[solver.rho->sgfn], block->d_fgfs[solver.Sx->sgfn], block->d_fgfs[solver.Sy->sgfn], block->d_fgfs[solver.Sz->sgfn],
			block->d_fgfs[solver.Sxx->sgfn], block->d_fgfs[solver.Sxy->sgfn], block->d_fgfs[solver.Sxz->sgfn],
			block->d_fgfs[solver.Syy->sgfn], block->d_fgfs[solver.Syz->sgfn], block->d_fgfs[solver.Szz->sgfn],
			block->d_fgfs[solver.Gamxxx->sgfn], block->d_fgfs[solver.Gamxxy->sgfn], block->d_fgfs[solver.Gamxxz->sgfn],
			block->d_fgfs[solver.Gamxyy->sgfn], block->d_fgfs[solver.Gamxyz->sgfn], block->d_fgfs[solver.Gamxzz->sgfn],
			block->d_fgfs[solver.Gamyxx->sgfn], block->d_fgfs[solver.Gamyxy->sgfn], block->d_fgfs[solver.Gamyxz->sgfn],
			block->d_fgfs[solver.Gamyyy->sgfn], block->d_fgfs[solver.Gamyyz->sgfn], block->d_fgfs[solver.Gamyzz->sgfn],
			block->d_fgfs[solver.Gamzxx->sgfn], block->d_fgfs[solver.Gamzxy->sgfn], block->d_fgfs[solver.Gamzxz->sgfn],
			block->d_fgfs[solver.Gamzyy->sgfn], block->d_fgfs[solver.Gamzyz->sgfn], block->d_fgfs[solver.Gamzzz->sgfn],
			block->d_fgfs[solver.Rxx->sgfn], block->d_fgfs[solver.Rxy->sgfn], block->d_fgfs[solver.Rxz->sgfn],
			block->d_fgfs[solver.Ryy->sgfn], block->d_fgfs[solver.Ryz->sgfn], block->d_fgfs[solver.Rzz->sgfn],
			block->d_fgfs[solver.Cons_Ham->sgfn],
			block->d_fgfs[solver.Cons_Px->sgfn], block->d_fgfs[solver.Cons_Py->sgfn], block->d_fgfs[solver.Cons_Pz->sgfn],
			block->d_fgfs[solver.Cons_Gx->sgfn], block->d_fgfs[solver.Cons_Gy->sgfn], block->d_fgfs[solver.Cons_Gz->sgfn],
			solver.Symmetry, lev, ndeps, co
		);
}

void run_predictor_block_step(bssn_class& solver, Patch* patch, Block* block, double trk4, double dT_lev, int lev, double ndeps, int pre, int cor, int iter_count, int yn) {
	trace_block_metric_det_summary(
		solver, "predictor_before_enforce_g0", lev, yn, block,
		solver.gxx0, solver.gxy0, solver.gxz0, solver.gyy0, solver.gyz0, solver.gzz0
	);
	gpu_enforce_ga_launch(
		block->stream, block->shape,
		block->d_fgfs[solver.gxx0->sgfn], block->d_fgfs[solver.gxy0->sgfn], block->d_fgfs[solver.gxz0->sgfn], block->d_fgfs[solver.gyy0->sgfn], block->d_fgfs[solver.gyz0->sgfn], block->d_fgfs[solver.gzz0->sgfn],
		block->d_fgfs[solver.Axx0->sgfn], block->d_fgfs[solver.Axy0->sgfn], block->d_fgfs[solver.Axz0->sgfn], block->d_fgfs[solver.Ayy0->sgfn], block->d_fgfs[solver.Ayz0->sgfn], block->d_fgfs[solver.Azz0->sgfn]
	);
	trace_block_metric_det_summary(
		solver, "predictor_after_enforce_g0", lev, yn, block,
		solver.gxx0, solver.gxy0, solver.gxz0, solver.gyy0, solver.gyz0, solver.gzz0
	);

	solver.Launch_RHS_GPU(block, trk4, ndeps, lev, pre, true);

	apply_stage_updates(solver, patch, block, dT_lev, lev, iter_count, solver.StateList, solver.StateList, solver.RHSList, solver.SynchList_pre, solver.RHSList, cor);
	gpu_lowerboundset_launch(block->stream, block->shape, block->d_fgfs[solver.phi->sgfn], solver.chitiny);
	trace_block_metric_det_summary(
		solver, "predictor_after_stage_pre", lev, yn, block,
		solver.gxx, solver.gxy, solver.gxz, solver.gyy, solver.gyz, solver.gzz
	);
}

void run_corrector_block_step(bssn_class& solver, Patch* patch, Block* block, double trk4, double dT_lev, double ndeps, int lev, int cor, int iter_count, int yn) {
	char phase[96];
	snprintf(phase, sizeof(phase), "corrector_before_enforce_iter%d", iter_count);
	trace_block_metric_det_summary(
		solver, phase, lev, yn, block,
		solver.gxx, solver.gxy, solver.gxz, solver.gyy, solver.gyz, solver.gzz
	);
	gpu_enforce_ga_launch(
		block->stream, block->shape,
		block->d_fgfs[solver.gxx->sgfn], block->d_fgfs[solver.gxy->sgfn], block->d_fgfs[solver.gxz->sgfn], block->d_fgfs[solver.gyy->sgfn], block->d_fgfs[solver.gyz->sgfn], block->d_fgfs[solver.gzz->sgfn],
		block->d_fgfs[solver.Axx->sgfn], block->d_fgfs[solver.Axy->sgfn], block->d_fgfs[solver.Axz->sgfn], block->d_fgfs[solver.Ayy->sgfn], block->d_fgfs[solver.Ayz->sgfn], block->d_fgfs[solver.Azz->sgfn]
	);
	snprintf(phase, sizeof(phase), "corrector_after_enforce_iter%d", iter_count);
	trace_block_metric_det_summary(
		solver, phase, lev, yn, block,
		solver.gxx, solver.gxy, solver.gxz, solver.gyy, solver.gyz, solver.gzz
	);

	solver.Launch_RHS_GPU(block, trk4, ndeps, lev, cor, false);

	apply_stage_updates(solver, patch, block, dT_lev, lev, iter_count, solver.StateList, solver.SynchList_pre, solver.SynchList_cor, solver.SynchList_cor, solver.RHSList, cor);
	gpu_lowerboundset_launch(block->stream, block->shape, block->d_fgfs[solver.phi1->sgfn], solver.chitiny);
	snprintf(phase, sizeof(phase), "corrector_after_stage_cor_iter%d", iter_count);
	trace_block_metric_det_summary(
		solver, phase, lev, yn, block,
		solver.gxx1, solver.gxy1, solver.gxz1, solver.gyy1, solver.gyz1, solver.gzz1
		);
}

double* shell_geom_field(Block* block, int fngfs, int field_index) {
	return block->d_fgfs[fngfs + field_index];
}

void launch_bssn_shell_rhs_for_block(
	bssn_class& solver,
	ss_patch* shell_patch,
	Block* block,
	double trk4,
	int lev,
	double eps,
	int co,
	var* chi_var,
	var* trK_var,
	var* gxx_var,
	var* gxy_var,
	var* gxz_var,
	var* gyy_var,
	var* gyz_var,
	var* gzz_var,
	var* Axx_var,
	var* Axy_var,
	var* Axz_var,
	var* Ayy_var,
	var* Ayz_var,
	var* Azz_var,
	var* Gamx_var,
	var* Gamy_var,
	var* Gamz_var,
	var* Lap_var,
	var* betax_var,
	var* betay_var,
	var* betaz_var,
	var* dtSfx_var,
	var* dtSfy_var,
	var* dtSfz_var,
	var* chi_rhs_var,
	var* trK_rhs_var,
	var* gxx_rhs_var,
	var* gxy_rhs_var,
	var* gxz_rhs_var,
	var* gyy_rhs_var,
	var* gyz_rhs_var,
	var* gzz_rhs_var,
	var* Axx_rhs_var,
	var* Axy_rhs_var,
	var* Axz_rhs_var,
	var* Ayy_rhs_var,
	var* Ayz_rhs_var,
	var* Azz_rhs_var,
	var* Gamx_rhs_var,
	var* Gamy_rhs_var,
	var* Gamz_rhs_var,
	var* Lap_rhs_var,
	var* betax_rhs_var,
	var* betay_rhs_var,
	var* betaz_rhs_var,
	var* dtSfx_rhs_var,
	var* dtSfy_rhs_var,
	var* dtSfz_rhs_var,
	var* TZ_var = 0,
	var* TZ_out_var = 0
) {
	const int fngfs = shell_patch->fngfs;
	gpu_compute_rhs_bssn_ss_launch(
		block->stream,
		block->shape, trk4,
		block->d_X[0], block->d_X[1], block->d_X[2],
		shell_geom_field(block, fngfs, ShellPatch::gx),
		shell_geom_field(block, fngfs, ShellPatch::gy),
		shell_geom_field(block, fngfs, ShellPatch::gz),
		shell_geom_field(block, fngfs, ShellPatch::drhodx),
		shell_geom_field(block, fngfs, ShellPatch::drhody),
		shell_geom_field(block, fngfs, ShellPatch::drhodz),
		shell_geom_field(block, fngfs, ShellPatch::dsigmadx),
		shell_geom_field(block, fngfs, ShellPatch::dsigmady),
		shell_geom_field(block, fngfs, ShellPatch::dsigmadz),
		shell_geom_field(block, fngfs, ShellPatch::dRdx),
		shell_geom_field(block, fngfs, ShellPatch::dRdy),
		shell_geom_field(block, fngfs, ShellPatch::dRdz),
		shell_geom_field(block, fngfs, ShellPatch::drhodxx),
		shell_geom_field(block, fngfs, ShellPatch::drhodxy),
		shell_geom_field(block, fngfs, ShellPatch::drhodxz),
		shell_geom_field(block, fngfs, ShellPatch::drhodyy),
		shell_geom_field(block, fngfs, ShellPatch::drhodyz),
		shell_geom_field(block, fngfs, ShellPatch::drhodzz),
		shell_geom_field(block, fngfs, ShellPatch::dsigmadxx),
		shell_geom_field(block, fngfs, ShellPatch::dsigmadxy),
		shell_geom_field(block, fngfs, ShellPatch::dsigmadxz),
		shell_geom_field(block, fngfs, ShellPatch::dsigmadyy),
		shell_geom_field(block, fngfs, ShellPatch::dsigmadyz),
		shell_geom_field(block, fngfs, ShellPatch::dsigmadzz),
		shell_geom_field(block, fngfs, ShellPatch::dRdxx),
		shell_geom_field(block, fngfs, ShellPatch::dRdxy),
		shell_geom_field(block, fngfs, ShellPatch::dRdxz),
		shell_geom_field(block, fngfs, ShellPatch::dRdyy),
		shell_geom_field(block, fngfs, ShellPatch::dRdyz),
		shell_geom_field(block, fngfs, ShellPatch::dRdzz),
		block->d_fgfs[chi_var->sgfn], block->d_fgfs[trK_var->sgfn],
		block->d_fgfs[gxx_var->sgfn], block->d_fgfs[gxy_var->sgfn], block->d_fgfs[gxz_var->sgfn],
		block->d_fgfs[gyy_var->sgfn], block->d_fgfs[gyz_var->sgfn], block->d_fgfs[gzz_var->sgfn],
		block->d_fgfs[Axx_var->sgfn], block->d_fgfs[Axy_var->sgfn], block->d_fgfs[Axz_var->sgfn],
		block->d_fgfs[Ayy_var->sgfn], block->d_fgfs[Ayz_var->sgfn], block->d_fgfs[Azz_var->sgfn],
		block->d_fgfs[Gamx_var->sgfn], block->d_fgfs[Gamy_var->sgfn], block->d_fgfs[Gamz_var->sgfn],
		block->d_fgfs[Lap_var->sgfn],
		block->d_fgfs[betax_var->sgfn], block->d_fgfs[betay_var->sgfn], block->d_fgfs[betaz_var->sgfn],
		block->d_fgfs[dtSfx_var->sgfn], block->d_fgfs[dtSfy_var->sgfn], block->d_fgfs[dtSfz_var->sgfn],
		block->d_fgfs[chi_rhs_var->sgfn], block->d_fgfs[trK_rhs_var->sgfn],
		block->d_fgfs[gxx_rhs_var->sgfn], block->d_fgfs[gxy_rhs_var->sgfn], block->d_fgfs[gxz_rhs_var->sgfn],
		block->d_fgfs[gyy_rhs_var->sgfn], block->d_fgfs[gyz_rhs_var->sgfn], block->d_fgfs[gzz_rhs_var->sgfn],
		block->d_fgfs[Axx_rhs_var->sgfn], block->d_fgfs[Axy_rhs_var->sgfn], block->d_fgfs[Axz_rhs_var->sgfn],
		block->d_fgfs[Ayy_rhs_var->sgfn], block->d_fgfs[Ayz_rhs_var->sgfn], block->d_fgfs[Azz_rhs_var->sgfn],
		block->d_fgfs[Gamx_rhs_var->sgfn], block->d_fgfs[Gamy_rhs_var->sgfn], block->d_fgfs[Gamz_rhs_var->sgfn],
		block->d_fgfs[Lap_rhs_var->sgfn],
		block->d_fgfs[betax_rhs_var->sgfn], block->d_fgfs[betay_rhs_var->sgfn], block->d_fgfs[betaz_rhs_var->sgfn],
		block->d_fgfs[dtSfx_rhs_var->sgfn], block->d_fgfs[dtSfy_rhs_var->sgfn], block->d_fgfs[dtSfz_rhs_var->sgfn],
		block->d_fgfs[solver.rho->sgfn], block->d_fgfs[solver.Sx->sgfn], block->d_fgfs[solver.Sy->sgfn], block->d_fgfs[solver.Sz->sgfn],
		block->d_fgfs[solver.Sxx->sgfn], block->d_fgfs[solver.Sxy->sgfn], block->d_fgfs[solver.Sxz->sgfn],
		block->d_fgfs[solver.Syy->sgfn], block->d_fgfs[solver.Syz->sgfn], block->d_fgfs[solver.Szz->sgfn],
		block->d_fgfs[solver.Gamxxx->sgfn], block->d_fgfs[solver.Gamxxy->sgfn], block->d_fgfs[solver.Gamxxz->sgfn],
		block->d_fgfs[solver.Gamxyy->sgfn], block->d_fgfs[solver.Gamxyz->sgfn], block->d_fgfs[solver.Gamxzz->sgfn],
		block->d_fgfs[solver.Gamyxx->sgfn], block->d_fgfs[solver.Gamyxy->sgfn], block->d_fgfs[solver.Gamyxz->sgfn],
		block->d_fgfs[solver.Gamyyy->sgfn], block->d_fgfs[solver.Gamyyz->sgfn], block->d_fgfs[solver.Gamyzz->sgfn],
		block->d_fgfs[solver.Gamzxx->sgfn], block->d_fgfs[solver.Gamzxy->sgfn], block->d_fgfs[solver.Gamzxz->sgfn],
		block->d_fgfs[solver.Gamzyy->sgfn], block->d_fgfs[solver.Gamzyz->sgfn], block->d_fgfs[solver.Gamzzz->sgfn],
		block->d_fgfs[solver.Rxx->sgfn], block->d_fgfs[solver.Rxy->sgfn], block->d_fgfs[solver.Rxz->sgfn],
		block->d_fgfs[solver.Ryy->sgfn], block->d_fgfs[solver.Ryz->sgfn], block->d_fgfs[solver.Rzz->sgfn],
		block->d_fgfs[solver.Cons_Ham->sgfn],
		block->d_fgfs[solver.Cons_Px->sgfn], block->d_fgfs[solver.Cons_Py->sgfn], block->d_fgfs[solver.Cons_Pz->sgfn],
		block->d_fgfs[solver.Cons_Gx->sgfn], block->d_fgfs[solver.Cons_Gy->sgfn], block->d_fgfs[solver.Cons_Gz->sgfn],
		solver.Symmetry, lev, eps, shell_patch->sst, co,
		TZ_var ? block->d_fgfs[TZ_var->sgfn] : 0,
		TZ_out_var ? block->d_fgfs[TZ_out_var->sgfn] : 0
	);
}

void launch_z4c_rhs_for_block(
	bssn_class& solver,
	Block* block,
	double trk4,
	int lev,
	double ndeps,
	int co,
	var* chi_var,
	var* trK_var,
	var* gxx_var,
	var* gxy_var,
	var* gxz_var,
	var* gyy_var,
	var* gyz_var,
	var* gzz_var,
	var* Axx_var,
	var* Axy_var,
	var* Axz_var,
	var* Ayy_var,
	var* Ayz_var,
	var* Azz_var,
	var* Gamx_var,
	var* Gamy_var,
	var* Gamz_var,
	var* Lap_var,
	var* betax_var,
	var* betay_var,
	var* betaz_var,
	var* dtSfx_var,
	var* dtSfy_var,
	var* dtSfz_var,
	var* TZ_var,
	var* chi_rhs_var,
	var* trK_rhs_var,
	var* gxx_rhs_var,
	var* gxy_rhs_var,
	var* gxz_rhs_var,
	var* gyy_rhs_var,
	var* gyz_rhs_var,
	var* gzz_rhs_var,
	var* Axx_rhs_var,
	var* Axy_rhs_var,
	var* Axz_rhs_var,
	var* Ayy_rhs_var,
	var* Ayz_rhs_var,
	var* Azz_rhs_var,
	var* Gamx_rhs_var,
	var* Gamy_rhs_var,
	var* Gamz_rhs_var,
	var* Lap_rhs_var,
	var* betax_rhs_var,
	var* betay_rhs_var,
	var* betaz_rhs_var,
	var* dtSfx_rhs_var,
	var* dtSfy_rhs_var,
	var* dtSfz_rhs_var,
	var* TZ_rhs_var
) {
	gpu_compute_rhs_z4c_launch(
		block->stream,
		block->shape, trk4, block->d_X[0], block->d_X[1], block->d_X[2],
		block->d_fgfs[chi_var->sgfn], block->d_fgfs[trK_var->sgfn],
		block->d_fgfs[gxx_var->sgfn], block->d_fgfs[gxy_var->sgfn], block->d_fgfs[gxz_var->sgfn],
		block->d_fgfs[gyy_var->sgfn], block->d_fgfs[gyz_var->sgfn], block->d_fgfs[gzz_var->sgfn],
		block->d_fgfs[Axx_var->sgfn], block->d_fgfs[Axy_var->sgfn], block->d_fgfs[Axz_var->sgfn],
		block->d_fgfs[Ayy_var->sgfn], block->d_fgfs[Ayz_var->sgfn], block->d_fgfs[Azz_var->sgfn],
		block->d_fgfs[Gamx_var->sgfn], block->d_fgfs[Gamy_var->sgfn], block->d_fgfs[Gamz_var->sgfn],
		block->d_fgfs[Lap_var->sgfn],
		block->d_fgfs[betax_var->sgfn], block->d_fgfs[betay_var->sgfn], block->d_fgfs[betaz_var->sgfn],
		block->d_fgfs[dtSfx_var->sgfn], block->d_fgfs[dtSfy_var->sgfn], block->d_fgfs[dtSfz_var->sgfn],
		block->d_fgfs[TZ_var->sgfn],
		block->d_fgfs[chi_rhs_var->sgfn], block->d_fgfs[trK_rhs_var->sgfn],
		block->d_fgfs[gxx_rhs_var->sgfn], block->d_fgfs[gxy_rhs_var->sgfn], block->d_fgfs[gxz_rhs_var->sgfn],
		block->d_fgfs[gyy_rhs_var->sgfn], block->d_fgfs[gyz_rhs_var->sgfn], block->d_fgfs[gzz_rhs_var->sgfn],
		block->d_fgfs[Axx_rhs_var->sgfn], block->d_fgfs[Axy_rhs_var->sgfn], block->d_fgfs[Axz_rhs_var->sgfn],
		block->d_fgfs[Ayy_rhs_var->sgfn], block->d_fgfs[Ayz_rhs_var->sgfn], block->d_fgfs[Azz_rhs_var->sgfn],
		block->d_fgfs[Gamx_rhs_var->sgfn], block->d_fgfs[Gamy_rhs_var->sgfn], block->d_fgfs[Gamz_rhs_var->sgfn],
		block->d_fgfs[Lap_rhs_var->sgfn],
		block->d_fgfs[betax_rhs_var->sgfn], block->d_fgfs[betay_rhs_var->sgfn], block->d_fgfs[betaz_rhs_var->sgfn],
		block->d_fgfs[dtSfx_rhs_var->sgfn], block->d_fgfs[dtSfy_rhs_var->sgfn], block->d_fgfs[dtSfz_rhs_var->sgfn],
		block->d_fgfs[TZ_rhs_var->sgfn],
		block->d_fgfs[solver.rho->sgfn], block->d_fgfs[solver.Sx->sgfn], block->d_fgfs[solver.Sy->sgfn], block->d_fgfs[solver.Sz->sgfn],
		block->d_fgfs[solver.Sxx->sgfn], block->d_fgfs[solver.Sxy->sgfn], block->d_fgfs[solver.Sxz->sgfn],
		block->d_fgfs[solver.Syy->sgfn], block->d_fgfs[solver.Syz->sgfn], block->d_fgfs[solver.Szz->sgfn],
		block->d_fgfs[solver.Gamxxx->sgfn], block->d_fgfs[solver.Gamxxy->sgfn], block->d_fgfs[solver.Gamxxz->sgfn],
		block->d_fgfs[solver.Gamxyy->sgfn], block->d_fgfs[solver.Gamxyz->sgfn], block->d_fgfs[solver.Gamxzz->sgfn],
		block->d_fgfs[solver.Gamyxx->sgfn], block->d_fgfs[solver.Gamyxy->sgfn], block->d_fgfs[solver.Gamyxz->sgfn],
		block->d_fgfs[solver.Gamyyy->sgfn], block->d_fgfs[solver.Gamyyz->sgfn], block->d_fgfs[solver.Gamyzz->sgfn],
		block->d_fgfs[solver.Gamzxx->sgfn], block->d_fgfs[solver.Gamzxy->sgfn], block->d_fgfs[solver.Gamzxz->sgfn],
		block->d_fgfs[solver.Gamzyy->sgfn], block->d_fgfs[solver.Gamzyz->sgfn], block->d_fgfs[solver.Gamzzz->sgfn],
		block->d_fgfs[solver.Rxx->sgfn], block->d_fgfs[solver.Rxy->sgfn], block->d_fgfs[solver.Rxz->sgfn],
		block->d_fgfs[solver.Ryy->sgfn], block->d_fgfs[solver.Ryz->sgfn], block->d_fgfs[solver.Rzz->sgfn],
		block->d_fgfs[solver.Cons_Ham->sgfn],
		block->d_fgfs[solver.Cons_Px->sgfn], block->d_fgfs[solver.Cons_Py->sgfn], block->d_fgfs[solver.Cons_Pz->sgfn],
		block->d_fgfs[solver.Cons_Gx->sgfn], block->d_fgfs[solver.Cons_Gy->sgfn], block->d_fgfs[solver.Cons_Gz->sgfn],
		solver.Symmetry, lev, ndeps, co
	);
}

void run_z4c_predictor_block_step(
	bssn_class& solver,
	Patch* patch,
	Block* block,
	double trk4,
	double dT_lev,
	int lev,
	double ndeps,
	int pre,
	int cor,
	int iter_count,
	int yn,
	var* TZ0_var,
	var* TZ_var,
	var* TZ_rhs_var
) {
	trace_block_metric_det_summary(
		solver, "z4c_predictor_before_enforce_g0", lev, yn, block,
		solver.gxx0, solver.gxy0, solver.gxz0, solver.gyy0, solver.gyz0, solver.gzz0
	);
	gpu_enforce_ga_launch(
		block->stream, block->shape,
		block->d_fgfs[solver.gxx0->sgfn], block->d_fgfs[solver.gxy0->sgfn], block->d_fgfs[solver.gxz0->sgfn], block->d_fgfs[solver.gyy0->sgfn], block->d_fgfs[solver.gyz0->sgfn], block->d_fgfs[solver.gzz0->sgfn],
		block->d_fgfs[solver.Axx0->sgfn], block->d_fgfs[solver.Axy0->sgfn], block->d_fgfs[solver.Axz0->sgfn], block->d_fgfs[solver.Ayy0->sgfn], block->d_fgfs[solver.Ayz0->sgfn], block->d_fgfs[solver.Azz0->sgfn]
	);
	trace_block_metric_det_summary(
		solver, "z4c_predictor_after_enforce_g0", lev, yn, block,
		solver.gxx0, solver.gxy0, solver.gxz0, solver.gyy0, solver.gyz0, solver.gzz0
	);

	launch_z4c_rhs_for_block(
		solver, block, trk4, lev, ndeps, pre,
		solver.phi0, solver.trK0,
		solver.gxx0, solver.gxy0, solver.gxz0, solver.gyy0, solver.gyz0, solver.gzz0,
		solver.Axx0, solver.Axy0, solver.Axz0, solver.Ayy0, solver.Ayz0, solver.Azz0,
		solver.Gmx0, solver.Gmy0, solver.Gmz0,
		solver.Lap0,
		solver.Sfx0, solver.Sfy0, solver.Sfz0,
		solver.dtSfx0, solver.dtSfy0, solver.dtSfz0,
		TZ0_var,
		solver.phi_rhs, solver.trK_rhs,
		solver.gxx_rhs, solver.gxy_rhs, solver.gxz_rhs, solver.gyy_rhs, solver.gyz_rhs, solver.gzz_rhs,
		solver.Axx_rhs, solver.Axy_rhs, solver.Axz_rhs, solver.Ayy_rhs, solver.Ayz_rhs, solver.Azz_rhs,
		solver.Gmx_rhs, solver.Gmy_rhs, solver.Gmz_rhs,
		solver.Lap_rhs,
		solver.Sfx_rhs, solver.Sfy_rhs, solver.Sfz_rhs,
		solver.dtSfx_rhs, solver.dtSfy_rhs, solver.dtSfz_rhs,
		TZ_rhs_var
	);

	apply_stage_updates(solver, patch, block, dT_lev, lev, iter_count, solver.StateList, solver.StateList, solver.RHSList, solver.SynchList_pre, solver.RHSList, cor);
	gpu_lowerboundset_launch(block->stream, block->shape, block->d_fgfs[solver.phi->sgfn], solver.chitiny);
	trace_block_metric_det_summary(
		solver, "z4c_predictor_after_stage_pre", lev, yn, block,
		solver.gxx, solver.gxy, solver.gxz, solver.gyy, solver.gyz, solver.gzz
	);
}

void run_z4c_corrector_block_step(
	bssn_class& solver,
	Patch* patch,
	Block* block,
	double trk4,
	double dT_lev,
	double ndeps,
	int lev,
	int cor,
	int iter_count,
	int yn,
	var* TZ_var,
	var* TZ1_var,
	var* TZ_rhs_var
) {
	char phase[96];
	snprintf(phase, sizeof(phase), "z4c_corrector_before_enforce_iter%d", iter_count);
	trace_block_metric_det_summary(
		solver, phase, lev, yn, block,
		solver.gxx, solver.gxy, solver.gxz, solver.gyy, solver.gyz, solver.gzz
	);
	gpu_enforce_ga_launch(
		block->stream, block->shape,
		block->d_fgfs[solver.gxx->sgfn], block->d_fgfs[solver.gxy->sgfn], block->d_fgfs[solver.gxz->sgfn], block->d_fgfs[solver.gyy->sgfn], block->d_fgfs[solver.gyz->sgfn], block->d_fgfs[solver.gzz->sgfn],
		block->d_fgfs[solver.Axx->sgfn], block->d_fgfs[solver.Axy->sgfn], block->d_fgfs[solver.Axz->sgfn], block->d_fgfs[solver.Ayy->sgfn], block->d_fgfs[solver.Ayz->sgfn], block->d_fgfs[solver.Azz->sgfn]
	);
	snprintf(phase, sizeof(phase), "z4c_corrector_after_enforce_iter%d", iter_count);
	trace_block_metric_det_summary(
		solver, phase, lev, yn, block,
		solver.gxx, solver.gxy, solver.gxz, solver.gyy, solver.gyz, solver.gzz
	);

	launch_z4c_rhs_for_block(
		solver, block, trk4, lev, ndeps, cor,
		solver.phi, solver.trK,
		solver.gxx, solver.gxy, solver.gxz, solver.gyy, solver.gyz, solver.gzz,
		solver.Axx, solver.Axy, solver.Axz, solver.Ayy, solver.Ayz, solver.Azz,
		solver.Gmx, solver.Gmy, solver.Gmz,
		solver.Lap,
		solver.Sfx, solver.Sfy, solver.Sfz,
		solver.dtSfx, solver.dtSfy, solver.dtSfz,
		TZ_var,
		solver.phi1, solver.trK1,
		solver.gxx1, solver.gxy1, solver.gxz1, solver.gyy1, solver.gyz1, solver.gzz1,
		solver.Axx1, solver.Axy1, solver.Axz1, solver.Ayy1, solver.Ayz1, solver.Azz1,
		solver.Gmx1, solver.Gmy1, solver.Gmz1,
		solver.Lap1,
		solver.Sfx1, solver.Sfy1, solver.Sfz1,
		solver.dtSfx1, solver.dtSfy1, solver.dtSfz1,
		TZ1_var
	);

	apply_stage_updates(solver, patch, block, dT_lev, lev, iter_count, solver.StateList, solver.SynchList_pre, solver.SynchList_cor, solver.SynchList_cor, solver.RHSList, cor);
	gpu_lowerboundset_launch(block->stream, block->shape, block->d_fgfs[solver.phi1->sgfn], solver.chitiny);
	snprintf(phase, sizeof(phase), "z4c_corrector_after_stage_cor_iter%d", iter_count);
	trace_block_metric_det_summary(
		solver, phase, lev, yn, block,
		solver.gxx1, solver.gxy1, solver.gxz1, solver.gyy1, solver.gyz1, solver.gzz1
	);
	(void)TZ_rhs_var;
}

void run_z4c_predictor_pass(
	bssn_class& solver,
	int lev,
	double trk4,
	double dT_lev,
	double ndeps,
	int pre,
	int cor,
	int iter_count,
	int yn,
	var* TZ0_var,
	var* TZ_var,
	var* TZ_rhs_var
) {
	MyList<Patch> *patch_node = solver.GH->PatL[lev];
	while (patch_node) {
		MyList<Block> *block_node = patch_node->data->blb;
		while (block_node) {
			Block *block = block_node->data;
			if (solver.myrank == block->rank)
				run_z4c_predictor_block_step(solver, patch_node->data, block, trk4, dT_lev, lev, ndeps, pre, cor, iter_count, yn, TZ0_var, TZ_var, TZ_rhs_var);
			if (block_node == patch_node->data->ble)
				break;
			block_node = block_node->next;
		}
		patch_node = patch_node->next;
	}
	(void)TZ_var;
}

void run_z4c_corrector_pass(
	bssn_class& solver,
	int lev,
	double trk4,
	double dT_lev,
	double ndeps,
	int cor,
	int iter_count,
	int yn,
	var* TZ_var,
	var* TZ1_var,
	var* TZ_rhs_var
) {
	MyList<Patch> *patch_node = solver.GH->PatL[lev];
	while (patch_node) {
		MyList<Block> *block_node = patch_node->data->blb;
		while (block_node) {
			Block *block = block_node->data;
			if (solver.myrank == block->rank)
				run_z4c_corrector_block_step(solver, patch_node->data, block, trk4, dT_lev, ndeps, lev, cor, iter_count, yn, TZ_var, TZ1_var, TZ_rhs_var);
			if (block_node == patch_node->data->ble)
				break;
			block_node = block_node->next;
		}
		patch_node = patch_node->next;
	}
}

void run_predictor_pass(bssn_class& solver, int lev, double trk4, double dT_lev, double ndeps, int pre, int cor, int iter_count, int yn) {
	MyList<Patch> *patch_node = solver.GH->PatL[lev];
	while (patch_node) {
		MyList<Block> *block_node = patch_node->data->blb;
		while (block_node) {
			Block *block = block_node->data;
			if (solver.myrank == block->rank)
				run_predictor_block_step(solver, patch_node->data, block, trk4, dT_lev, lev, ndeps, pre, cor, iter_count, yn);
			if (block_node == patch_node->data->ble)
				break;
			block_node = block_node->next;
		}
		patch_node = patch_node->next;
	}
}

void run_corrector_pass(bssn_class& solver, int lev, double trk4, double dT_lev, double ndeps, int cor, int iter_count, int yn) {
	MyList<Patch> *patch_node = solver.GH->PatL[lev];
	while (patch_node) {
		MyList<Block> *block_node = patch_node->data->blb;
		while (block_node) {
			Block *block = block_node->data;
			if (solver.myrank == block->rank)
				run_corrector_block_step(solver, patch_node->data, block, trk4, dT_lev, ndeps, lev, cor, iter_count, yn);
			if (block_node == patch_node->data->ble)
				break;
			block_node = block_node->next;
		}
		patch_node = patch_node->next;
	}
}

void apply_shell_stage_updates(
	bssn_class& solver,
	ss_patch* shell_patch,
	Block* block,
	double dT_lev,
	int iter_count,
	MyList<var> *base_state_list,
	MyList<var> *bam_source_list,
	MyList<var> *bam_target_list,
	MyList<var> *rk_target_list,
	MyList<var> *rhs_list
) {
	MyList<var> *base_node = base_state_list;
	MyList<var> *bam_source_node = bam_source_list;
	MyList<var> *bam_target_node = bam_target_list;
	MyList<var> *rk_target_node = rk_target_list;
	MyList<var> *rhs_node = rhs_list;
	while (base_node) {
		gpu_sommerfeld_routbam_ss_launch(
			block->stream,
			block->shape,
			block->d_X[0], block->d_X[1], block->d_X[2],
			shell_patch->bbox[0], shell_patch->bbox[1], shell_patch->bbox[2],
			shell_patch->bbox[3], shell_patch->bbox[4], shell_patch->bbox[5],
			block->d_fgfs[bam_target_node->data->sgfn],
			block->d_fgfs[bam_source_node->data->sgfn],
			base_node->data->propspeed, base_node->data->SoA,
			solver.Symmetry
		);
		gpu_rungekutta4_rout_launch(
			block->stream,
			block->shape, dT_lev,
			block->d_fgfs[base_node->data->sgfn],
			block->d_fgfs[rk_target_node->data->sgfn],
			block->d_fgfs[rhs_node->data->sgfn],
			iter_count
		);
		base_node = base_node->next;
		bam_source_node = bam_source_node->next;
		bam_target_node = bam_target_node->next;
		rk_target_node = rk_target_node->next;
		rhs_node = rhs_node->next;
	}
}

void launch_bssn_shell_rhs_for_stage(
	bssn_class& solver,
	ss_patch* shell_patch,
	Block* block,
	double trk4,
	int lev,
	double eps,
	int co,
	bool predictor,
	var* TZ_var = 0,
	var* TZ_out_var = 0
) {
	launch_bssn_shell_rhs_for_block(
		solver, shell_patch, block, trk4, lev, eps, co,
		predictor ? solver.phi0 : solver.phi,
		predictor ? solver.trK0 : solver.trK,
		predictor ? solver.gxx0 : solver.gxx,
		predictor ? solver.gxy0 : solver.gxy,
		predictor ? solver.gxz0 : solver.gxz,
		predictor ? solver.gyy0 : solver.gyy,
		predictor ? solver.gyz0 : solver.gyz,
		predictor ? solver.gzz0 : solver.gzz,
		predictor ? solver.Axx0 : solver.Axx,
		predictor ? solver.Axy0 : solver.Axy,
		predictor ? solver.Axz0 : solver.Axz,
		predictor ? solver.Ayy0 : solver.Ayy,
		predictor ? solver.Ayz0 : solver.Ayz,
		predictor ? solver.Azz0 : solver.Azz,
		predictor ? solver.Gmx0 : solver.Gmx,
		predictor ? solver.Gmy0 : solver.Gmy,
		predictor ? solver.Gmz0 : solver.Gmz,
		predictor ? solver.Lap0 : solver.Lap,
		predictor ? solver.Sfx0 : solver.Sfx,
		predictor ? solver.Sfy0 : solver.Sfy,
		predictor ? solver.Sfz0 : solver.Sfz,
		predictor ? solver.dtSfx0 : solver.dtSfx,
		predictor ? solver.dtSfy0 : solver.dtSfy,
		predictor ? solver.dtSfz0 : solver.dtSfz,
		predictor ? solver.phi_rhs : solver.phi1,
		predictor ? solver.trK_rhs : solver.trK1,
		predictor ? solver.gxx_rhs : solver.gxx1,
		predictor ? solver.gxy_rhs : solver.gxy1,
		predictor ? solver.gxz_rhs : solver.gxz1,
		predictor ? solver.gyy_rhs : solver.gyy1,
		predictor ? solver.gyz_rhs : solver.gyz1,
		predictor ? solver.gzz_rhs : solver.gzz1,
		predictor ? solver.Axx_rhs : solver.Axx1,
		predictor ? solver.Axy_rhs : solver.Axy1,
		predictor ? solver.Axz_rhs : solver.Axz1,
		predictor ? solver.Ayy_rhs : solver.Ayy1,
		predictor ? solver.Ayz_rhs : solver.Ayz1,
		predictor ? solver.Azz_rhs : solver.Azz1,
		predictor ? solver.Gmx_rhs : solver.Gmx1,
		predictor ? solver.Gmy_rhs : solver.Gmy1,
		predictor ? solver.Gmz_rhs : solver.Gmz1,
		predictor ? solver.Lap_rhs : solver.Lap1,
		predictor ? solver.Sfx_rhs : solver.Sfx1,
		predictor ? solver.Sfy_rhs : solver.Sfy1,
		predictor ? solver.Sfz_rhs : solver.Sfz1,
		predictor ? solver.dtSfx_rhs : solver.dtSfx1,
		predictor ? solver.dtSfy_rhs : solver.dtSfy1,
		predictor ? solver.dtSfz_rhs : solver.dtSfz1,
		TZ_var,
		TZ_out_var
	);
}

void run_shell_predictor_pass(bssn_class& solver, int lev, double trk4, double dT_lev, int pre, int iter_count, int yn) {
	(void)yn;
#ifdef WithShell
	if (lev != 0 || !solver.SH)
		return;

	MyList<ss_patch> *shell_node = solver.SH->PatL;
	while (shell_node) {
		MyList<Block> *block_node = shell_node->data->blb;
		while (block_node) {
			Block *block = block_node->data;
			if (solver.myrank == block->rank) {
				gpu_enforce_ga_launch(
					block->stream, block->shape,
					block->d_fgfs[solver.gxx0->sgfn], block->d_fgfs[solver.gxy0->sgfn], block->d_fgfs[solver.gxz0->sgfn],
					block->d_fgfs[solver.gyy0->sgfn], block->d_fgfs[solver.gyz0->sgfn], block->d_fgfs[solver.gzz0->sgfn],
					block->d_fgfs[solver.Axx0->sgfn], block->d_fgfs[solver.Axy0->sgfn], block->d_fgfs[solver.Axz0->sgfn],
					block->d_fgfs[solver.Ayy0->sgfn], block->d_fgfs[solver.Ayz0->sgfn], block->d_fgfs[solver.Azz0->sgfn]
				);
				launch_bssn_shell_rhs_for_stage(solver, shell_node->data, block, trk4, lev, solver.numepsh, pre, true);
				apply_shell_stage_updates(
					solver, shell_node->data, block, dT_lev, iter_count,
					solver.StateList, solver.StateList, solver.RHSList, solver.SynchList_pre, solver.RHSList
				);
				gpu_lowerboundset_launch(block->stream, block->shape, block->d_fgfs[solver.phi->sgfn], solver.chitiny);
			}
			if (block_node == shell_node->data->ble)
				break;
			block_node = block_node->next;
		}
		shell_node = shell_node->next;
	}
#else
	(void)solver;
	(void)lev;
	(void)trk4;
	(void)dT_lev;
	(void)pre;
	(void)iter_count;
#endif
}

void run_shell_corrector_pass(bssn_class& solver, int lev, double trk4, double dT_lev, int cor, int iter_count, int yn) {
	(void)yn;
#ifdef WithShell
	if (lev != 0 || !solver.SH)
		return;

	MyList<ss_patch> *shell_node = solver.SH->PatL;
	while (shell_node) {
		MyList<Block> *block_node = shell_node->data->blb;
		while (block_node) {
			Block *block = block_node->data;
			if (solver.myrank == block->rank) {
				gpu_enforce_ga_launch(
					block->stream, block->shape,
					block->d_fgfs[solver.gxx->sgfn], block->d_fgfs[solver.gxy->sgfn], block->d_fgfs[solver.gxz->sgfn],
					block->d_fgfs[solver.gyy->sgfn], block->d_fgfs[solver.gyz->sgfn], block->d_fgfs[solver.gzz->sgfn],
					block->d_fgfs[solver.Axx->sgfn], block->d_fgfs[solver.Axy->sgfn], block->d_fgfs[solver.Axz->sgfn],
					block->d_fgfs[solver.Ayy->sgfn], block->d_fgfs[solver.Ayz->sgfn], block->d_fgfs[solver.Azz->sgfn]
				);
				launch_bssn_shell_rhs_for_stage(solver, shell_node->data, block, trk4, lev, solver.numepsh, cor, false);
				apply_shell_stage_updates(
					solver, shell_node->data, block, dT_lev, iter_count,
					solver.StateList, solver.SynchList_pre, solver.SynchList_cor, solver.SynchList_cor, solver.RHSList
				);
				gpu_lowerboundset_launch(block->stream, block->shape, block->d_fgfs[solver.phi1->sgfn], solver.chitiny);
			}
			if (block_node == shell_node->data->ble)
				break;
			block_node = block_node->next;
		}
		shell_node = shell_node->next;
	}
#else
	(void)solver;
	(void)lev;
	(void)trk4;
	(void)dT_lev;
	(void)cor;
	(void)iter_count;
#endif
}

void run_z4c_shell_predictor_pass(
	bssn_class& solver,
	int lev,
	double trk4,
	double dT_lev,
	int pre,
	int iter_count,
	int yn,
	var* TZ0_var,
	var* TZ_var,
	var* TZ_rhs_var
) {
	(void)yn;
	(void)TZ_var;
#ifdef WithShell
	if (lev != 0 || !solver.SH)
		return;

	MyList<ss_patch> *shell_node = solver.SH->PatL;
	while (shell_node) {
		MyList<Block> *block_node = shell_node->data->blb;
		while (block_node) {
			Block *block = block_node->data;
			if (solver.myrank == block->rank) {
				gpu_enforce_ga_launch(
					block->stream, block->shape,
					block->d_fgfs[solver.gxx0->sgfn], block->d_fgfs[solver.gxy0->sgfn], block->d_fgfs[solver.gxz0->sgfn],
					block->d_fgfs[solver.gyy0->sgfn], block->d_fgfs[solver.gyz0->sgfn], block->d_fgfs[solver.gzz0->sgfn],
					block->d_fgfs[solver.Axx0->sgfn], block->d_fgfs[solver.Axy0->sgfn], block->d_fgfs[solver.Axz0->sgfn],
					block->d_fgfs[solver.Ayy0->sgfn], block->d_fgfs[solver.Ayz0->sgfn], block->d_fgfs[solver.Azz0->sgfn]
				);
				launch_bssn_shell_rhs_for_stage(
					solver, shell_node->data, block, trk4, lev, solver.numepsh, pre, true,
					TZ0_var, TZ_rhs_var
				);
				apply_shell_stage_updates(
					solver, shell_node->data, block, dT_lev, iter_count,
					solver.StateList, solver.StateList, solver.RHSList, solver.SynchList_pre, solver.RHSList
				);
				gpu_lowerboundset_launch(block->stream, block->shape, block->d_fgfs[solver.phi->sgfn], solver.chitiny);
			}
			if (block_node == shell_node->data->ble)
				break;
			block_node = block_node->next;
		}
		shell_node = shell_node->next;
	}
#else
	(void)solver;
	(void)lev;
	(void)trk4;
	(void)dT_lev;
	(void)pre;
	(void)iter_count;
#endif
}

void run_z4c_shell_corrector_pass(
	bssn_class& solver,
	int lev,
	double trk4,
	double dT_lev,
	int cor,
	int iter_count,
	int yn,
	var* TZ_var,
	var* TZ1_var,
	var* TZ_rhs_var
) {
	(void)yn;
	(void)TZ_rhs_var;
#ifdef WithShell
	if (lev != 0 || !solver.SH)
		return;

	MyList<ss_patch> *shell_node = solver.SH->PatL;
	while (shell_node) {
		MyList<Block> *block_node = shell_node->data->blb;
		while (block_node) {
			Block *block = block_node->data;
			if (solver.myrank == block->rank) {
				gpu_enforce_ga_launch(
					block->stream, block->shape,
					block->d_fgfs[solver.gxx->sgfn], block->d_fgfs[solver.gxy->sgfn], block->d_fgfs[solver.gxz->sgfn],
					block->d_fgfs[solver.gyy->sgfn], block->d_fgfs[solver.gyz->sgfn], block->d_fgfs[solver.gzz->sgfn],
					block->d_fgfs[solver.Axx->sgfn], block->d_fgfs[solver.Axy->sgfn], block->d_fgfs[solver.Axz->sgfn],
					block->d_fgfs[solver.Ayy->sgfn], block->d_fgfs[solver.Ayz->sgfn], block->d_fgfs[solver.Azz->sgfn]
				);
				launch_bssn_shell_rhs_for_stage(
					solver, shell_node->data, block, trk4, lev, solver.numepsh, cor, false,
					TZ_var, TZ1_var
				);
				apply_shell_stage_updates(
					solver, shell_node->data, block, dT_lev, iter_count,
					solver.StateList, solver.SynchList_pre, solver.SynchList_cor, solver.SynchList_cor, solver.RHSList
				);
				gpu_lowerboundset_launch(block->stream, block->shape, block->d_fgfs[solver.phi1->sgfn], solver.chitiny);
			}
			if (block_node == shell_node->data->ble)
				break;
			block_node = block_node->next;
		}
		shell_node = shell_node->next;
	}
#else
	(void)solver;
	(void)lev;
	(void)trk4;
	(void)dT_lev;
	(void)cor;
	(void)iter_count;
#endif
}

void reduce_stage_error(int& error_flag) {
	int reduced_error = error_flag;
	MPI_Allreduce(&reduced_error, &error_flag, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
}

void abort_stage_error(bssn_class& solver, int lev, double dT_lev, MyList<var> *dump_list, const char* message_prefix, int iter_count) {
	Parallel::Dump_Data(solver.GH->PatL[lev], dump_list, 0, solver.PhysTime, dT_lev);
	if (solver.myrank == 0) {
		if (solver.ErrorMonitor->outfile) {
			if (iter_count >= 0)
				solver.ErrorMonitor->outfile << message_prefix << iter_count << " variables at t = " << solver.PhysTime << ", lev = " << lev << endl;
			else
				solver.ErrorMonitor->outfile << message_prefix << solver.PhysTime << ", lev = " << lev << endl;
		}
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
}

void swap_patch_block_lists(bssn_class& solver, int lev, MyList<var> *left_list, MyList<var> *right_list) {
	MyList<Patch> *patch_node = solver.GH->PatL[lev];
	while (patch_node) {
		MyList<Block> *block_node = patch_node->data->blb;
		while (block_node) {
			block_node->data->swapList(left_list, right_list, solver.myrank);
			if (block_node == patch_node->data->ble)
				break;
			block_node = block_node->next;
		}
		patch_node = patch_node->next;
	}
}

void swap_shell_block_lists(bssn_class& solver, MyList<var> *left_list, MyList<var> *right_list) {
#ifdef WithShell
	if (!solver.SH)
		return;

	MyList<ss_patch> *shell_node = solver.SH->PatL;
	while (shell_node) {
		MyList<Block> *block_node = shell_node->data->blb;
		while (block_node) {
			block_node->data->swapList(left_list, right_list, solver.myrank);
			if (block_node == shell_node->data->ble)
				break;
			block_node = block_node->next;
		}
		shell_node = shell_node->next;
	}
#else
	(void)solver;
	(void)left_list;
	(void)right_list;
#endif
}

void mark_gpu_modified_list(bssn_class& solver, int lev, MyList<var> *variable_list) {
	if (!variable_list)
		return;

	MyList<Patch> *patch_node = solver.GH->PatL[lev];
	while (patch_node) {
		MyList<Block> *block_node = patch_node->data->blb;
		while (block_node) {
			Block *block = block_node->data;
			if (solver.myrank == block->rank) {
				MyList<var> *var_node = variable_list;
				while (var_node) {
					const int sgfn = var_node->data->sgfn;
					if (sgfn >= 0 && sgfn < block->fngfs)
						block->mark_gpu_modified(sgfn);
					var_node = var_node->next;
				}
			}
			if (block_node == patch_node->data->ble)
				break;
			block_node = block_node->next;
		}
		patch_node = patch_node->next;
	}
}

void mark_shell_gpu_modified_list(bssn_class& solver, MyList<var> *variable_list) {
#ifdef WithShell
	if (!variable_list || !solver.SH)
		return;

	MyList<ss_patch> *shell_node = solver.SH->PatL;
	while (shell_node) {
		MyList<Block> *block_node = shell_node->data->blb;
		while (block_node) {
			Block *block = block_node->data;
			if (solver.myrank == block->rank) {
				MyList<var> *var_node = variable_list;
				while (var_node) {
					const int sgfn = var_node->data->sgfn;
					if (sgfn >= 0 && sgfn < block->fngfs)
						block->mark_gpu_modified(sgfn);
					var_node = var_node->next;
				}
			}
			if (block_node == shell_node->data->ble)
				break;
			block_node = block_node->next;
		}
		shell_node = shell_node->next;
	}
#else
	(void)solver;
	(void)variable_list;
#endif
}

void mark_bssn_step_device_outputs(bssn_class& solver, int lev) {
	mark_gpu_modified_list(solver, lev, solver.StateList);
	mark_gpu_modified_list(solver, lev, solver.OldStateList);
	mark_gpu_modified_list(solver, lev, solver.RHSList);
	mark_gpu_modified_list(solver, lev, solver.SynchList_pre);
	mark_gpu_modified_list(solver, lev, solver.SynchList_cor);
	mark_gpu_modified_list(solver, lev, solver.MiscList);
	mark_gpu_modified_list(solver, lev, solver.ConstraintList);
}

void mark_bssn_shell_step_device_outputs(bssn_class& solver) {
	mark_shell_gpu_modified_list(solver, solver.StateList);
	mark_shell_gpu_modified_list(solver, solver.OldStateList);
	mark_shell_gpu_modified_list(solver, solver.RHSList);
	mark_shell_gpu_modified_list(solver, solver.SynchList_pre);
	mark_shell_gpu_modified_list(solver, solver.SynchList_cor);
	mark_shell_gpu_modified_list(solver, solver.MiscList);
	mark_shell_gpu_modified_list(solver, solver.ConstraintList);
}

void mark_z4c_step_device_outputs(Z4c_gpu_class& solver, int lev) {
	mark_gpu_modified_list(solver, lev, solver.StateList);
	mark_gpu_modified_list(solver, lev, solver.OldStateList);
	mark_gpu_modified_list(solver, lev, solver.RHSList);
	mark_gpu_modified_list(solver, lev, solver.SynchList_pre);
	mark_gpu_modified_list(solver, lev, solver.SynchList_cor);
	mark_gpu_modified_list(solver, lev, solver.ConstraintList);
}

void commit_black_hole_predictor(bssn_class& solver, int lev) {
	if (solver.BH_num <= 0 || lev != solver.GH->levels - 1)
		return;

	for (int black_hole_index = 0; black_hole_index < solver.BH_num; ++black_hole_index) {
		solver.Porg0[black_hole_index][0] = solver.Porg1[black_hole_index][0];
		solver.Porg0[black_hole_index][1] = solver.Porg1[black_hole_index][1];
		solver.Porg0[black_hole_index][2] = solver.Porg1[black_hole_index][2];
	}
}

int step_trace_env_int(const char* name, int fallback) {
	const char* value = getenv(name);
	return value ? atoi(value) : fallback;
}

bool step_nonfinite_trace_in_scope(int lev, int yn) {
	if (!getenv("AMSS_TRAJECTORY_PREFIX_TRACE_STEP_NONFINITE"))
		return false;
	const int focus_lev = step_trace_env_int("AMSS_TRAJECTORY_PREFIX_TRACE_STEP_NONFINITE_LEVEL", -1);
	if (focus_lev >= 0 && focus_lev != lev)
		return false;
	const int focus_yn = step_trace_env_int("AMSS_TRAJECTORY_PREFIX_TRACE_STEP_NONFINITE_YN", -1);
	if (focus_yn >= 0 && focus_yn != yn)
		return false;
	return true;
}

bool step_det_trace_in_scope(int lev, int yn) {
	if (!getenv("AMSS_TRAJECTORY_PREFIX_TRACE_STEP_DET"))
		return false;
	const int focus_lev = step_trace_env_int("AMSS_TRAJECTORY_PREFIX_TRACE_STEP_DET_LEVEL", -1);
	if (focus_lev >= 0 && focus_lev != lev)
		return false;
	const int focus_yn = step_trace_env_int("AMSS_TRAJECTORY_PREFIX_TRACE_STEP_DET_YN", -1);
	if (focus_yn >= 0 && focus_yn != yn)
		return false;
	return true;
}

double step_trace_env_double(const char* name, double fallback) {
	const char* value = getenv(name);
	return value ? atof(value) : fallback;
}

double metric_det_from_components(double gxx, double gxy, double gxz, double gyy, double gyz, double gzz) {
	return gxx * (gyy * gzz - gyz * gyz)
	     - gxy * (gxy * gzz - gyz * gxz)
	     + gxz * (gxy * gyz - gyy * gxz);
}

void trace_block_metric_det_summary(
	bssn_class& solver,
	const char* phase,
	int lev,
	int yn,
	Block* block,
	var* gxx_var,
	var* gxy_var,
	var* gxz_var,
	var* gyy_var,
	var* gyz_var,
	var* gzz_var
) {
	if (!step_det_trace_in_scope(lev, yn) || !block)
		return;

	var* vars[6] = {gxx_var, gxy_var, gxz_var, gyy_var, gyz_var, gzz_var};
	for (int component = 0; component < 6; ++component) {
		if (!vars[component])
			return;
		const int sgfn = vars[component]->sgfn;
		if (sgfn < 0 || sgfn >= block->fngfs)
			return;
	}

	const int sgxx = gxx_var->sgfn;
	const int sgxy = gxy_var->sgfn;
	const int sgxz = gxz_var->sgfn;
	const int sgyy = gyy_var->sgfn;
	const int sgyz = gyz_var->sgfn;
	const int sgzz = gzz_var->sgfn;
	const int nx = block->shape[0];
	const int ny = block->shape[1];
	const int nz = block->shape[2];
	const int count = nx * ny * nz;
	const double det_floor = step_trace_env_double("AMSS_TRAJECTORY_PREFIX_TRACE_STEP_DET_FLOOR", 0.0);
	double* metric_host[6];
	const int sgfns[6] = {sgxx, sgxy, sgxz, sgyy, sgyz, sgzz};
	for (int component = 0; component < 6; ++component) {
		metric_host[component] = new double[count];
		GpuRuntimeCoordinator::instance().copy_to_host(
			metric_host[component], block->d_fgfs[sgfns[component]], count, block->stream
		);
	}
	GpuRuntimeCoordinator::wait_for_stream(block->stream);

	int nonfinite_metric_count = 0;
	int nonfinite_det_count = 0;
	int nonpositive_det_count = 0;
	int first_i = -1;
	int first_j = -1;
	int first_k = -1;
	double first_x = 0.0;
	double first_y = 0.0;
	double first_z = 0.0;
	double first_det = 0.0;
	double first_metric[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	double min_det = 0.0;
	double max_det = 0.0;
	bool have_det = false;

	for (int idx = 0; idx < count; ++idx) {
		const double gxx = metric_host[0][idx] + 1.0;
		const double gxy = metric_host[1][idx];
		const double gxz = metric_host[2][idx];
		const double gyy = metric_host[3][idx] + 1.0;
		const double gyz = metric_host[4][idx];
		const double gzz = metric_host[5][idx] + 1.0;
		const bool metric_finite = finite(gxx) && finite(gxy) && finite(gxz) && finite(gyy) && finite(gyz) && finite(gzz);
		if (!metric_finite)
			++nonfinite_metric_count;
		const double det = metric_det_from_components(gxx, gxy, gxz, gyy, gyz, gzz);
		if (!finite(det))
			++nonfinite_det_count;
		else {
			if (!have_det) {
				min_det = det;
				max_det = det;
				have_det = true;
			} else {
				if (det < min_det)
					min_det = det;
				if (det > max_det)
					max_det = det;
			}
			if (det <= det_floor)
				++nonpositive_det_count;
		}

		const bool bad_metric_det = !metric_finite || !finite(det) || det <= det_floor;
		if (!bad_metric_det || first_i >= 0)
			continue;

		const int i = idx % nx;
		const int j = (idx / nx) % ny;
		const int k = idx / (nx * ny);
		first_i = i;
		first_j = j;
		first_k = k;
		first_x = block->X[0][i];
		first_y = block->X[1][j];
		first_z = block->X[2][k];
		first_det = det;
		first_metric[0] = gxx;
		first_metric[1] = gxy;
		first_metric[2] = gxz;
		first_metric[3] = gyy;
		first_metric[4] = gyz;
		first_metric[5] = gzz;
	}

	if (!getenv("AMSS_TRAJECTORY_PREFIX_TRACE_STEP_DET_PRINT_CLEAN") &&
	    nonfinite_metric_count == 0 && nonfinite_det_count == 0 && nonpositive_det_count == 0) {
		for (int component = 0; component < 6; ++component)
			delete[] metric_host[component];
		return;
	}

	cout << " step_metric_det"
	     << " phase=" << phase
	     << " lev=" << lev
	     << " yn=" << yn
	     << " phys_time=" << solver.PhysTime;
	if (lev >= 0 && lev < solver.GH->levels)
		cout << " level_time=" << solver.GH->Lt[lev];
	cout << " shape=(" << nx << "," << ny << "," << nz << ")"
	     << " det_floor=" << det_floor
	     << " nonfinite_metric_count=" << nonfinite_metric_count
	     << " nonfinite_det_count=" << nonfinite_det_count
	     << " nonpositive_det_count=" << nonpositive_det_count
	     << " min_det=" << min_det
	     << " max_det=" << max_det
	     << " first_ijk=(" << first_i << "," << first_j << "," << first_k << ")"
	     << " first_xyz=(" << first_x << "," << first_y << "," << first_z << ")"
	     << " first_det=" << first_det
	     << " first_metric=(" << first_metric[0] << "," << first_metric[1] << "," << first_metric[2] << ","
	     << first_metric[3] << "," << first_metric[4] << "," << first_metric[5] << ")"
	     << endl;

	for (int component = 0; component < 6; ++component)
		delete[] metric_host[component];
}

void trace_step_nonfinite_summary(bssn_class& solver, const char* phase, int lev, int yn, MyList<var>* variable_list) {
	if (!step_nonfinite_trace_in_scope(lev, yn) || !solver.GH || !variable_list)
		return;

	for (MyList<Patch>* patch_node = solver.GH->PatL[lev]; patch_node; patch_node = patch_node->next) {
		for (MyList<Block>* block_node = patch_node->data->blb; block_node; block_node = block_node->next) {
			Block* block = block_node->data;
			if (solver.myrank == block->rank) {
				for (MyList<var>* var_node = variable_list; var_node; var_node = var_node->next) {
					const int sgfn = var_node->data->sgfn;
					if (sgfn < 0 || sgfn >= block->fngfs)
						continue;
					block->mark_gpu_modified(sgfn);
					block->require_on_cpu(sgfn);
				}
			}
			if (block_node == patch_node->data->ble)
				break;
		}
	}
	GpuRuntimeCoordinator::instance().wait_for_device();

	int total_nonfinite = 0;
	int variables_with_nonfinite = 0;
	int first_block = -1;
	int first_sgfn = -1;
	const char* first_var = "none";
	int first_i = -1;
	int first_j = -1;
	int first_k = -1;
	double first_x = 0.0;
	double first_y = 0.0;
	double first_z = 0.0;
	double first_value = 0.0;

	int block_id = 0;
	for (MyList<Patch>* patch_node = solver.GH->PatL[lev]; patch_node; patch_node = patch_node->next) {
		for (MyList<Block>* block_node = patch_node->data->blb; block_node; block_node = block_node->next, ++block_id) {
			Block* block = block_node->data;
			if (solver.myrank == block->rank) {
				const int nx = block->shape[0];
				const int ny = block->shape[1];
				const int count = block->shape[0] * block->shape[1] * block->shape[2];
				for (MyList<var>* var_node = variable_list; var_node; var_node = var_node->next) {
					const int sgfn = var_node->data->sgfn;
					if (sgfn < 0 || sgfn >= block->fngfs)
						continue;
					int var_nonfinite = 0;
					for (int idx = 0; idx < count; ++idx) {
						const double value = block->fgfs[sgfn][idx];
						if (finite(value))
							continue;
						++var_nonfinite;
						if (first_block < 0) {
							const int i = idx % nx;
							const int j = (idx / nx) % ny;
							const int k = idx / (nx * ny);
							first_block = block_id;
							first_sgfn = sgfn;
							first_var = var_node->data->name;
							first_i = i;
							first_j = j;
							first_k = k;
							first_x = block->X[0][i];
							first_y = block->X[1][j];
							first_z = block->X[2][k];
							first_value = value;
						}
					}
					if (var_nonfinite > 0) {
						++variables_with_nonfinite;
						total_nonfinite += var_nonfinite;
					}
				}
			}
			if (block_node == patch_node->data->ble)
				break;
		}
	}

	cout << " step_nonfinite"
	     << " phase=" << phase
	     << " lev=" << lev
	     << " yn=" << yn
	     << " phys_time=" << solver.PhysTime;
	if (lev >= 0 && lev < solver.GH->levels)
		cout << " level_time=" << solver.GH->Lt[lev];
	cout << " total_nonfinite=" << total_nonfinite
	     << " variables_with_nonfinite=" << variables_with_nonfinite
	     << " first_block=" << first_block
	     << " first_var=" << first_var
	     << " first_sgfn=" << first_sgfn
	     << " first_ijk=(" << first_i << "," << first_j << "," << first_k << ")"
	     << " first_xyz=(" << first_x << "," << first_y << "," << first_z << ")"
	     << " first_value=" << first_value
	     << endl;
	}

	}

void bssn_class::Launch_RHS_GPU(Block *block, double trk4, double ndeps, int lev, int co, bool predictor) {
	launch_bssn_rhs_for_block(
		*this, block, trk4, lev, ndeps, co,
		predictor ? phi0 : phi,
		predictor ? trK0 : trK,
		predictor ? gxx0 : gxx,
		predictor ? gxy0 : gxy,
		predictor ? gxz0 : gxz,
		predictor ? gyy0 : gyy,
		predictor ? gyz0 : gyz,
		predictor ? gzz0 : gzz,
		predictor ? Axx0 : Axx,
		predictor ? Axy0 : Axy,
		predictor ? Axz0 : Axz,
		predictor ? Ayy0 : Ayy,
		predictor ? Ayz0 : Ayz,
		predictor ? Azz0 : Azz,
		predictor ? Gmx0 : Gmx,
		predictor ? Gmy0 : Gmy,
		predictor ? Gmz0 : Gmz,
		predictor ? Lap0 : Lap,
		predictor ? Sfx0 : Sfx,
		predictor ? Sfy0 : Sfy,
		predictor ? Sfz0 : Sfz,
		predictor ? dtSfx0 : dtSfx,
		predictor ? dtSfy0 : dtSfy,
		predictor ? dtSfz0 : dtSfz,
		predictor ? phi_rhs : phi1,
		predictor ? trK_rhs : trK1,
		predictor ? gxx_rhs : gxx1,
		predictor ? gxy_rhs : gxy1,
		predictor ? gxz_rhs : gxz1,
		predictor ? gyy_rhs : gyy1,
		predictor ? gyz_rhs : gyz1,
		predictor ? gzz_rhs : gzz1,
		predictor ? Axx_rhs : Axx1,
		predictor ? Axy_rhs : Axy1,
		predictor ? Axz_rhs : Axz1,
		predictor ? Ayy_rhs : Ayy1,
		predictor ? Ayz_rhs : Ayz1,
		predictor ? Azz_rhs : Azz1,
		predictor ? Gmx_rhs : Gmx1,
		predictor ? Gmy_rhs : Gmy1,
		predictor ? Gmz_rhs : Gmz1,
		predictor ? Lap_rhs : Lap1,
		predictor ? Sfx_rhs : Sfx1,
		predictor ? Sfy_rhs : Sfy1,
		predictor ? Sfz_rhs : Sfz1,
		predictor ? dtSfx_rhs : dtSfx1,
		predictor ? dtSfy_rhs : dtSfy1,
		predictor ? dtSfz_rhs : dtSfz1
	);
}

#ifdef With_AHF
void bssn_class::AH_Prepare_derivatives_GPU()
{
	for (int lev = 0; lev < GH->levels; ++lev) {
		MyList<Patch> *Pp = GH->PatL[lev];
		while (Pp) {
			MyList<Block> *BP = Pp->data->blb;
			while (BP) {
				Block *cg = BP->data;
				if (myrank == cg->rank) {
					gpu_ah_prepare_derivatives_launch(
						cg->stream,
						cg->shape, cg->d_X[0], cg->d_X[1], cg->d_X[2],
						cg->d_fgfs[phi0->sgfn],
						cg->d_fgfs[dtSfx_rhs->sgfn], cg->d_fgfs[dtSfy_rhs->sgfn], cg->d_fgfs[dtSfz_rhs->sgfn],
						cg->d_fgfs[gxx0->sgfn],
						cg->d_fgfs[Gamxxx->sgfn], cg->d_fgfs[Gamyxx->sgfn], cg->d_fgfs[Gamzxx->sgfn],
						cg->d_fgfs[gxy0->sgfn],
						cg->d_fgfs[Gamxxy->sgfn], cg->d_fgfs[Gamyxy->sgfn], cg->d_fgfs[Gamzxy->sgfn],
						cg->d_fgfs[gxz0->sgfn],
						cg->d_fgfs[Gamxxz->sgfn], cg->d_fgfs[Gamyxz->sgfn], cg->d_fgfs[Gamzxz->sgfn],
						cg->d_fgfs[gyy0->sgfn],
						cg->d_fgfs[Gamxyy->sgfn], cg->d_fgfs[Gamyyy->sgfn], cg->d_fgfs[Gamzyy->sgfn],
						cg->d_fgfs[gyz0->sgfn],
						cg->d_fgfs[Gamxyz->sgfn], cg->d_fgfs[Gamyyz->sgfn], cg->d_fgfs[Gamzyz->sgfn],
						cg->d_fgfs[gzz0->sgfn],
						cg->d_fgfs[Gamxzz->sgfn], cg->d_fgfs[Gamyzz->sgfn], cg->d_fgfs[Gamzzz->sgfn],
						Symmetry);
				}
				if (BP == Pp->data->ble)
					break;
				BP = BP->next;
			}
			Pp = Pp->next;
		}
		GpuRuntimeCoordinator::instance().wait_for_device();
		Parallel::Sync_GPU(GH->PatL[lev], AHDList, Symmetry);
	}
}

bool bssn_class::AH_Interp_Points_GPU(MyList<var> *VarList,
	                                  int NN, double **XX,
	                                  double *Shellf, int Symmetryi)
{
	int num_var = 0;
	for (MyList<var> *varl = VarList; varl; varl = varl->next)
		++num_var;

	if (NN <= 0 || num_var <= 0)
		return true;

#ifdef WithShell
	if (SH) {
		// Horizon sample points live around the punctures, deep inside the
		// Cartesian levels, so the Cartesian-only interpolation below is
		// sufficient. There is no ShellPatch fallback yet: a point escaping
		// the Cartesian domain is reported by the fails-to-find check below.
		static bool shell_interp_warned = false;
		if (!shell_interp_warned) {
			shell_interp_warned = true;
			if (myrank == 0)
				cout << "AH_Interp_Points_GPU: ShellPatch fallback is not implemented; "
				     << "horizon points outside the Cartesian grid will abort." << endl;
			if (ErrorMonitor && ErrorMonitor->outfile)
				ErrorMonitor->outfile << "AH_Interp_Points_GPU: ShellPatch fallback is not implemented; "
				                      << "horizon points outside the Cartesian grid will abort." << endl;
		}
	}
#endif

	cudaStream_t sample_stream = porg_sample_stream;
	ensure_porg_sample_capacity(NN, num_var);

	int *active_points = h_porg_active_bh;
	int *next_active_points = h_porg_next_active_bh;
	int *point_weight = h_porg_point_weight;
	for (int n = 0; n < NN; ++n)
		active_points[n] = n;

	int active_count = NN;
	for (int lev = GH->levels - 1; lev >= 0 && active_count > 0; --lev) {
		for (int active_index = 0; active_index < active_count; ++active_index) {
			const int n = active_points[active_index];
			h_porg_points[active_index] = XX[0][n];
			h_porg_points[active_count + active_index] = XX[1][n];
			h_porg_points[2 * active_count + active_index] = XX[2][n];
		}

		GpuRuntimeCoordinator::instance().copy_to_device(
			h_porg_points, d_porg_points, 3 * active_count, sample_stream);

		double *d_pox[3] = {
			d_porg_points,
			d_porg_points + active_count,
			d_porg_points + 2 * active_count
		};

		Parallel::PatList_Interp_Points_GPU(
			sample_stream,
			GH->PatL[lev],
			VarList,
			active_count,
			d_pox,
			d_porg_shellf,
			h_porg_shellf,
			point_weight,
			Symmetryi
		);

		int next_active_count = 0;
		for (int active_index = 0; active_index < active_count; ++active_index) {
			const int n = active_points[active_index];
			if (point_weight[active_index] > 0) {
				for (int variable_index = 0; variable_index < num_var; ++variable_index)
					Shellf[n * num_var + variable_index] = h_porg_shellf[active_index * num_var + variable_index];
			}
			else {
				next_active_points[next_active_count] = n;
				++next_active_count;
			}
		}

		int *tmp_active_points = active_points;
		active_points = next_active_points;
		next_active_points = tmp_active_points;
		active_count = next_active_count;
	}

	if (active_count > 0) {
		for (int active_index = 0; active_index < active_count; ++active_index) {
			const int n = active_points[active_index];
			if (ErrorMonitor && ErrorMonitor->outfile)
				ErrorMonitor->outfile << "AH_Interp_Points_GPU fails to find point ("
				                      << XX[0][n] << "," << XX[1][n] << "," << XX[2][n]
				                      << ") at t = " << PhysTime << endl;
			if (myrank == 0)
				cout << "AH_Interp_Points_GPU fails to find point ("
				     << XX[0][n] << "," << XX[1][n] << "," << XX[2][n] << ")" << endl;
		}
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	GpuRuntimeCoordinator::wait_for_stream(sample_stream);
	return true;
}

void bssn_class::AH_Step_Find_GPU(int lev, double dT_lev)
{
	if (!(lev == GH->levels - 1))
		return;

	int ncount = int(PhysTime / dT_lev);
	bool tf = false;
	for (int ihn = 0; ihn < HN_num; ihn++) {
		if (ncount % findeveryl[ihn] == 0) {
			tf = true;
			break;
		}
	}
	if (!tf)
		return;

	double prev_clock, curr_clock;
	if (myrank == 0)
		prev_clock = MPI_Wtime();

	const int cdumpid = int(PhysTime / AHdumptime) + 1;
	for (int ihn = 0; ihn < HN_num; ihn++)
		dumpid[ihn] = cdumpid;

	double gam;
	for (int ihn = 0; ihn < BH_num; ihn++) {
		xc[ihn] = Porg0[ihn][0];
		yc[ihn] = Porg0[ihn][1];
		zc[ihn] = Porg0[ihn][2];
		gam = fabs(Pmom[ihn * 3]) / Mass[ihn];
		gam = sqrt(1 - gam * gam);
		xr[ihn] = Mass[ihn] * gam;
		gam = fabs(Pmom[ihn * 3 + 1]) / Mass[ihn];
		gam = sqrt(1 - gam * gam);
		yr[ihn] = Mass[ihn] * gam;
		gam = fabs(Pmom[ihn * 3 + 2]) / Mass[ihn];
		gam = sqrt(1 - gam * gam);
		zr[ihn] = Mass[ihn] * gam;
		dTT[ihn] = -1;

		if (ncount % findeveryl[ihn] == 0) {
			trigger[ihn] = true;
			dTT[ihn] = findeveryl[ihn] * dT_lev;
		}
		else {
			trigger[ihn] = false;
		}

		if (trigger[ihn] && dumpid[ihn] > lastahdumpid[ihn])
			lastahdumpid[ihn] = dumpid[ihn];
		else
			dumpid[ihn] = 0;
	}

	int ihn = BH_num;
	for (int ia = 0; ia < BH_num; ia++) {
		for (int ib = ia + 1; ib < BH_num; ib++) {
			xc[ihn] = (Porg0[ia][0] + Porg0[ib][0]) / 2;
			yc[ihn] = (Porg0[ia][1] + Porg0[ib][1]) / 2;
			zc[ihn] = (Porg0[ia][2] + Porg0[ib][2]) / 2;
			xr[ihn] = yr[ihn] = zr[ihn] = Mass[ia] + Mass[ib];
			dTT[ihn] = -1;

			if (fabs(Porg0[ia][0] - Porg0[ib][0]) < 2 * xr[ihn] &&
			    fabs(Porg0[ia][1] - Porg0[ib][1]) < 2 * xr[ihn] &&
			    fabs(Porg0[ia][2] - Porg0[ib][2]) < 2 * xr[ihn] &&
			    ncount % findeveryl[ihn] == 0) {
				trigger[ihn] = true;
				dTT[ihn] = findeveryl[ihn] * dT_lev;
			}
			else {
				trigger[ihn] = false;
			}

			if (trigger[ihn] && dumpid[ihn] > lastahdumpid[ihn])
				lastahdumpid[ihn] = dumpid[ihn];
			else
				dumpid[ihn] = 0;

			ihn++;
		}
	}
#if (ABEtype == 1)
	if (PhysTime > 10) {
		ihn--;
		trigger[ihn] = true;
		xr[ihn] = yr[ihn] = zr[ihn] = 50;
	}
#endif

	(void)&bssn_class::AH_Interp_Points_GPU;
	(void)&gpu_ahf_radial_residual_newton_launch;

	AHFinderDirect::AHFinderDirect_find_horizons(HN_num, dumpid,
		xc, yc, zc, xr, yr, zr, trigger, dTT);

	if (myrank == 0)
	{
		curr_clock = MPI_Wtime();
		cout << " Finding horizon used " << (curr_clock - prev_clock)
			<< " seconds!" << endl;
		cout << endl;
	}
}

void bssn_class::AH_Step_Find(int lev, double dT_lev)
{
	(void)lev;
	(void)dT_lev;
	if (myrank == 0) {
		cout << "GPU ABE build received CPU AH_Step_Find; use AH_Step_Find_GPU instead." << endl;
		if (ErrorMonitor && ErrorMonitor->outfile)
			ErrorMonitor->outfile << "GPU ABE build received CPU AH_Step_Find; use AH_Step_Find_GPU instead." << endl;
	}
	MPI_Abort(MPI_COMM_WORLD, 1);
}
#endif

bool bssn_class::Use_RK4_Puncture_Position_GPU() const {
#if defined(MAPBH) && (MAPBH == 0)
	return true;
#else
	return false;
#endif
}

bool Z4c_gpu_class::Use_RK4_Puncture_Position_GPU() const {
	return true;
}

void bssn_class::Step_GPU(int lev, int YN) {
	setpbh(BH_num, Porg0, Mass, BH_num_input);
	const bool use_rk4_puncture_position = Use_RK4_Puncture_Position_GPU();

	double dT_lev = dT * pow(0.5, Mymax(lev, trfls));
#ifdef With_AHF
	AH_Step_Find_GPU(lev, dT_lev);
#endif
	if (!use_rk4_puncture_position) {
		advance_black_hole_euler_predictor(*this, lev, dT_lev);
	}

	if (!use_rk4_puncture_position && lev == a_lev) {
		AnalysisStuff(lev, dT_lev);
	}

	bool BB = fgt(PhysTime, StartTime, dT_lev / 2);
	double ndeps = numepss;
	if (lev < GH->movls)
		ndeps = numepsb;
	double TRK4 = PhysTime;
	int iter_count = 0;
	int pre = 0, cor = 1;
	int ERROR = 0;

	trace_step_nonfinite_summary(*this, "entry_state", lev, YN, StateList);
	run_predictor_pass(*this, lev, TRK4, dT_lev, ndeps, pre, cor, iter_count, YN);
#ifdef WithShell
	run_shell_predictor_pass(*this, lev, TRK4, dT_lev, pre, iter_count, YN);
#endif
	GpuRuntimeCoordinator::instance().wait_for_device();
	trace_step_nonfinite_summary(*this, "after_predictor_rhs", lev, YN, RHSList);
	trace_step_nonfinite_summary(*this, "after_predictor_pre", lev, YN, SynchList_pre);
	reduce_stage_error(ERROR);
	if (ERROR)
		abort_stage_error(*this, lev, dT_lev, StateList, "find NaN in state variables at t = ", -1);
	Parallel::Sync_GPU(GH->PatL[lev], SynchList_pre, Symmetry);
#ifdef WithShell
	if (lev == 0)
	{
		double prev_clock, curr_clock;
		if (myrank == 0)
			prev_clock = MPI_Wtime();
		SH->Synch_GPU(SynchList_pre, Symmetry);
		if (myrank == 0)
		{
			curr_clock = MPI_Wtime();
			cout << "Shell stuff synchronization used " << (curr_clock - prev_clock) << " seconds! (RK substep " << iter_count << ")" << endl;
		}
	}
#endif
	if (use_rk4_puncture_position) {
		advance_black_hole_rk4_predictor(*this, lev, dT_lev);
		if (lev == a_lev)
			AnalysisStuff(lev, dT_lev);
	}
	trace_step_nonfinite_summary(*this, "after_predictor_sync_pre", lev, YN, SynchList_pre);
	for (iter_count = 1; iter_count < 4; iter_count++) {
		if (iter_count == 1 || iter_count == 3)
			TRK4 += dT_lev / 2;
		char phase[96];
		run_corrector_pass(*this, lev, TRK4, dT_lev, ndeps, cor, iter_count, YN);
#ifdef WithShell
		run_shell_corrector_pass(*this, lev, TRK4, dT_lev, cor, iter_count, YN);
#endif
		GpuRuntimeCoordinator::instance().wait_for_device();
		snprintf(phase, sizeof(phase), "after_corrector_rhs_iter%d", iter_count);
		trace_step_nonfinite_summary(*this, phase, lev, YN, RHSList);
		snprintf(phase, sizeof(phase), "after_corrector_cor_iter%d", iter_count);
		trace_step_nonfinite_summary(*this, phase, lev, YN, SynchList_cor);
		reduce_stage_error(ERROR);
		if (ERROR)
			abort_stage_error(*this, lev, dT_lev, SynchList_pre, "find NaN in RK4 substep#", iter_count);
		Parallel::Sync_GPU(GH->PatL[lev], SynchList_cor, Symmetry);
#ifdef WithShell
		if (lev == 0)
		{
			double prev_clock, curr_clock;
			if (myrank == 0)
				prev_clock = MPI_Wtime();
			SH->Synch_GPU(SynchList_cor, Symmetry);
			if (myrank == 0)
			{
				curr_clock = MPI_Wtime();
				cout << "Shell stuff synchronization used " << (curr_clock - prev_clock) << " seconds! (RK substep " << iter_count << ")" << endl;
				if (iter_count == 3) cout << endl;
			}
		}
#endif
		if (use_rk4_puncture_position)
			advance_black_hole_rk4_corrector(*this, lev, dT_lev, iter_count);
		snprintf(phase, sizeof(phase), "after_corrector_sync_cor_iter%d", iter_count);
		trace_step_nonfinite_summary(*this, phase, lev, YN, SynchList_cor);
		if (iter_count < 3) {
			swap_patch_block_lists(*this, lev, SynchList_pre, SynchList_cor);
#ifdef WithShell
			if (lev == 0)
				swap_shell_block_lists(*this, SynchList_pre, SynchList_cor);
#endif
			if (use_rk4_puncture_position)
				prepare_black_hole_corrector_for_next_stage(*this, lev);
			snprintf(phase, sizeof(phase), "after_stage_swap_pre_iter%d", iter_count);
			trace_step_nonfinite_summary(*this, phase, lev, YN, SynchList_pre);
		}
	}
	trace_step_nonfinite_summary(*this, "before_final_swap_state", lev, YN, StateList);
	trace_step_nonfinite_summary(*this, "before_final_swap_cor", lev, YN, SynchList_cor);
	swap_patch_block_lists(*this, lev, StateList, SynchList_cor);
#ifdef WithShell
	if (lev == 0)
		swap_shell_block_lists(*this, StateList, SynchList_cor);
#endif
	trace_step_nonfinite_summary(*this, "after_state_swap_state", lev, YN, StateList);
	swap_patch_block_lists(*this, lev, OldStateList, SynchList_cor);
#ifdef WithShell
	if (lev == 0)
		swap_shell_block_lists(*this, OldStateList, SynchList_cor);
#endif
	trace_step_nonfinite_summary(*this, "after_old_swap_state", lev, YN, StateList);
	mark_bssn_step_device_outputs(*this, lev);
#ifdef WithShell
	if (lev == 0)
		mark_bssn_shell_step_device_outputs(*this);
#endif
	commit_black_hole_predictor(*this, lev);
}

void Z4c_gpu_class::Step_GPU(int lev, int YN) {
	setpbh(BH_num, Porg0, Mass, BH_num_input);
	const bool use_rk4_puncture_position = Use_RK4_Puncture_Position_GPU();

	double dT_lev = dT * pow(0.5, Mymax(lev, trfls));
#ifdef With_AHF
	AH_Step_Find_GPU(lev, dT_lev);
#endif

	if (!use_rk4_puncture_position) {
		advance_black_hole_euler_predictor(*this, lev, dT_lev);
	}

	if (!use_rk4_puncture_position && lev == a_lev) {
		AnalysisStuff(lev, dT_lev);
	}

	double ndeps = numepss;
	if (lev < GH->movls)
		ndeps = numepsb;
	double TRK4 = PhysTime;
	int iter_count = 0;
	int pre = 0, cor = 1;
	int ERROR = 0;

	trace_step_nonfinite_summary(*this, "z4c_entry_state", lev, YN, StateList);
	run_z4c_predictor_pass(*this, lev, TRK4, dT_lev, ndeps, pre, cor, iter_count, YN, TZ0, TZ, TZ_rhs);
#ifdef WithShell
	run_z4c_shell_predictor_pass(*this, lev, TRK4, dT_lev, pre, iter_count, YN, TZ0, TZ, TZ_rhs);
#endif
	GpuRuntimeCoordinator::instance().wait_for_device();
	trace_step_nonfinite_summary(*this, "z4c_after_predictor_rhs", lev, YN, RHSList);
	trace_step_nonfinite_summary(*this, "z4c_after_predictor_pre", lev, YN, SynchList_pre);
	reduce_stage_error(ERROR);
	if (ERROR)
		abort_stage_error(*this, lev, dT_lev, StateList, "find NaN in Z4C state variables at t = ", -1);
	Parallel::Sync_GPU(GH->PatL[lev], SynchList_pre, Symmetry);
#ifdef WithShell
	if (lev == 0)
	{
		double prev_clock, curr_clock;
		if (myrank == 0)
			prev_clock = MPI_Wtime();
		SH->Synch_GPU(SynchList_pre, Symmetry);
		if (myrank == 0)
		{
			curr_clock = MPI_Wtime();
			cout << "Shell stuff synchronization used " << (curr_clock - prev_clock) << " seconds! (RK substep " << iter_count << ")" << endl;
		}
	}
#endif
	if (use_rk4_puncture_position) {
		advance_black_hole_rk4_predictor(*this, lev, dT_lev);
		if (lev == a_lev)
			AnalysisStuff(lev, dT_lev);
	}
	trace_step_nonfinite_summary(*this, "z4c_after_predictor_sync_pre", lev, YN, SynchList_pre);
	for (iter_count = 1; iter_count < 4; iter_count++) {
		if (iter_count == 1 || iter_count == 3)
			TRK4 += dT_lev / 2;
		char phase[96];
		run_z4c_corrector_pass(*this, lev, TRK4, dT_lev, ndeps, cor, iter_count, YN, TZ, TZ1, TZ_rhs);
#ifdef WithShell
		run_z4c_shell_corrector_pass(*this, lev, TRK4, dT_lev, cor, iter_count, YN, TZ, TZ1, TZ_rhs);
#endif
		GpuRuntimeCoordinator::instance().wait_for_device();
		snprintf(phase, sizeof(phase), "z4c_after_corrector_rhs_iter%d", iter_count);
		trace_step_nonfinite_summary(*this, phase, lev, YN, RHSList);
		snprintf(phase, sizeof(phase), "z4c_after_corrector_cor_iter%d", iter_count);
		trace_step_nonfinite_summary(*this, phase, lev, YN, SynchList_cor);
		reduce_stage_error(ERROR);
		if (ERROR)
			abort_stage_error(*this, lev, dT_lev, SynchList_pre, "find NaN in Z4C RK4 substep#", iter_count);
		Parallel::Sync_GPU(GH->PatL[lev], SynchList_cor, Symmetry);
#ifdef WithShell
		if (lev == 0)
		{
			double prev_clock, curr_clock;
			if (myrank == 0)
				prev_clock = MPI_Wtime();
			SH->Synch_GPU(SynchList_cor, Symmetry);
			if (myrank == 0)
			{
				curr_clock = MPI_Wtime();
				cout << "Shell stuff synchronization used " << (curr_clock - prev_clock) << " seconds! (RK substep " << iter_count << ")" << endl;
				if (iter_count == 3) cout << endl;
			}
		}
#endif
		if (use_rk4_puncture_position)
			advance_black_hole_rk4_corrector(*this, lev, dT_lev, iter_count);
		snprintf(phase, sizeof(phase), "z4c_after_corrector_sync_cor_iter%d", iter_count);
		trace_step_nonfinite_summary(*this, phase, lev, YN, SynchList_cor);
		if (iter_count < 3) {
			swap_patch_block_lists(*this, lev, SynchList_pre, SynchList_cor);
#ifdef WithShell
			if (lev == 0)
				swap_shell_block_lists(*this, SynchList_pre, SynchList_cor);
#endif
			if (use_rk4_puncture_position)
				prepare_black_hole_corrector_for_next_stage(*this, lev);
			snprintf(phase, sizeof(phase), "z4c_after_stage_swap_pre_iter%d", iter_count);
			trace_step_nonfinite_summary(*this, phase, lev, YN, SynchList_pre);
		}
	}
	trace_step_nonfinite_summary(*this, "z4c_before_final_swap_state", lev, YN, StateList);
	trace_step_nonfinite_summary(*this, "z4c_before_final_swap_cor", lev, YN, SynchList_cor);
	swap_patch_block_lists(*this, lev, StateList, SynchList_cor);
#ifdef WithShell
	if (lev == 0)
		swap_shell_block_lists(*this, StateList, SynchList_cor);
#endif
	trace_step_nonfinite_summary(*this, "z4c_after_state_swap_state", lev, YN, StateList);
	swap_patch_block_lists(*this, lev, OldStateList, SynchList_cor);
#ifdef WithShell
	if (lev == 0)
		swap_shell_block_lists(*this, OldStateList, SynchList_cor);
#endif
	trace_step_nonfinite_summary(*this, "z4c_after_old_swap_state", lev, YN, StateList);
	mark_z4c_step_device_outputs(*this, lev);
#ifdef WithShell
	if (lev == 0)
		mark_bssn_shell_step_device_outputs(*this);
#endif
	commit_black_hole_predictor(*this, lev);
}
