

include makefile.inc

.SUFFIXES: .o .f90 .C .for .cu

.f90.o:
	$(f90) $(f90appflags) -c $< -o $@

.C.o:
	${CXX} $(CXXAPPFLAGS) -c $< $(filein) -o $@

.cpp.o:
	${CXX} $(CXXAPPFLAGS) -c $< $(filein) -o $@

.for.o:
	$(f77) -c $< -o $@

.cu.o:
	$(Cu) $(CUDA_APP_FLAGS) -c $< -o $@ $(CUDA_LIB_PATH)

# Input files
C++FILES_GPU = ABE.o Ansorg.o Block.o misc.o monitor.o Parallel.o MPatch.o var.o\
           cgh.o surface_integral.o ShellPatch.o\
	   bssnEScalar_class.o perf.o Z4c_class.o NullShellPatch.o\
	   bssnEM_gpu_class.o cpbc_util.o z4c_rhs_point.o checkpoint.o\
           Parallel_bam.o scalar_class.o transpbh.o NullShellPatch2.o\
	   NullShellPatch2_Evo.o \
	   bssn_gpu_class.o bssn_step_gpu.o bssn_macro.o writefile_f.o\
	   bssnEScalar_gpu_class.o bssnEScalar_step_gpu.o\
	   Z4c_step_gpu.o Z4c_gpu_class.o\
	   Parallel_GPU.o helper.o ShellPatch_GPU.o

F90FILES = enforce_algebra.o fmisc.o initial_puncture.o prolongrestrict.o\
	   prolongrestrict_cell.o prolongrestrict_vertex.o\
	   rungekutta4_rout.o bssn_rhs.o diff_new.o kodiss.o kodiss_sh.o\
	   lopsidediff.o sommerfeld_rout.o getnp4.o diff_new_sh.o\
	   shellfunctions.o bssn_rhs_ss.o Set_Rho_ADM.o\
           getnp4EScalar.o bssnEScalar_rhs.o bssn_constraint.o ricci_gamma.o\
           fadmquantites_bssn.o Z4c_rhs.o Z4c_rhs_ss.o point_diff_new_sh.o\
	   cpbc.o getnp4old.o NullEvol.o initial_null.o initial_maxwell.o\
	   getnpem2.o empart.o NullNews.o fourdcurvature.o\
	   bssn2adm.o adm_constraint.o adm_ricci_gamma.o\
	   scalar_rhs.o initial_scalar.o NullEvol2.o initial_null2.o\
	   NullNews2.o tool_f.o

F77FILES = zbesh.o

AHFDOBJS = AHF/expansion.o AHF/expansion_Jacobian.o AHF/patch.o AHF/coords.o AHF/patch_info.o AHF/patch_interp.o AHF/patch_system.o \
AHF/tgrid.o AHF/fd_grid.o AHF/ghost_zone.o AHF/array.o AHF/round.o AHF/norm.o AHF/fuzzy.o AHF/error_exit.o AHF/miscfp.o \
AHF/linear_map.o AHF/cpm_map.o AHF/BH_diagnostics.o AHF/setup.o AHF/horizon_sequence.o AHF/find_horizons.o \
AHF/initial_guess.o AHF/Newton.o AHF/Jacobian.o AHF/ilucg.o AHF/IntPnts0.o AHF/IntPnts.o AHF_GPU/ahf_gpu.o

TwoPunctureFILES = ./TwoPuncture_OMP/TwoPunctureABE.o ./TwoPuncture_OMP/TwoPunctures.o

PNOrbit_FILES = ./PN_Orbit/PN_orbit.o ./PN_Orbit/misc.o ./PN_Orbit/lagpolint.o

CUDAFILES = fmisc_gpu.o bssn_rhs_gpu.o bssn_rhs_ss_gpu.o bssnEScalar_rhs_gpu.o Z4c_rhs_gpu.o Z4c_rhs_ss_gpu.o bssn_constraint_gpu.o prolongrestrict_cell_gpu.o prolongrestrict_vertex_gpu.o gpu_manager.o enforce_algebra_gpu.o sommerfeld_rout_gpu.o rungekutta4_rout_gpu.o fadmquantites_bssn_gpu.o MPatch_gpu.o ShellPatch_gpu.o surface_integral_gpu.o getnp4_gpu.o diff_new_sh_gpu.o cpbc_gpu.o ricci_gamma_gpu.o getnp4EScalar_gpu.o getnpem2_gpu.o empart_gpu.o adm_ricci_gamma_gpu.o adm_constraint_gpu.o fourdcurvature_gpu.o Set_Rho_ADM_gpu.o

CUDAFILES += AHF_GPU/ahf_derivatives_gpu.o AHF_GPU/ahf_sampling_gpu.o AHF_GPU/ahf_surface_gpu.o AHF_GPU/ahf_ilucg_gpu.o
AHFDOBJS += AHF_GPU/ahf_ordered_plan.o AHF_GPU/ahf_diagnostics_plan.o
CUDAFILES += AHF_GPU/ahf_ordered_csr.o AHF_GPU/ahf_newton_gpu.o AHF_GPU/ahf_diagnostics_gpu.o
AHF_GPU_LIBS = -lcusparse -lcusolver
# Preserve the evolution shell compilation recipe.
ShellPatch_gpu.o:
	$(Cu) $(CUDA_APP_FLAGS) -rdc=true -c $(@:.o=.cu) -o $@ $(CUDA_LIB_PATH)
$(AHFDOBJS): CXXAPPFLAGS += -I$(CUDA_DIR)/include
AHF_SURFACE_CONSUMERS = bssn_gpu_class.o AHF/setup.o AHF/find_horizons.o AHF_GPU/ahf_gpu.o AHF_GPU/ahf_surface_gpu.o AHF_GPU/ahf_sampling_gpu.o AHF_GPU/ahf_derivatives_gpu.o
$(AHF_SURFACE_CONSUMERS): AHF_GPU/ahf_gpu.h AHF_GPU/ahf_native_types.h AHF_GPU/ahf_ordered_csr.h AHF_GPU/ahf_newton_gpu.h AHF_GPU/ahf_diagnostics_gpu.h
AHF_GPU/ahf_ilucg_gpu.o AHF_GPU/ahf_surface_gpu.o AHF_GPU/ahf_gpu.o: AHF_GPU/ahf_ilucg_gpu.h
AHF_STRICT_CUDA = AHF_GPU/ahf_sampling_gpu.o AHF_GPU/ahf_surface_gpu.o AHF_GPU/ahf_ilucg_gpu.o AHF_GPU/ahf_ordered_csr.o AHF_GPU/ahf_newton_gpu.o AHF_GPU/ahf_diagnostics_gpu.o
AHF_STRICT_HOST = AHF_GPU/ahf_gpu.o AHF_GPU/ahf_ordered_plan.o AHF_GPU/ahf_diagnostics_plan.o
# AHF 严格浮点标志默认关闭（提速）；需要与旧构建逐位对齐时：make AHF_STRICT_FP=1
ifeq ($(AHF_STRICT_FP),1)
$(AHF_STRICT_CUDA): private CUDA_APP_FLAGS += --fmad=false --prec-div=true --prec-sqrt=true -Xcompiler=-ffp-contract=off,-fno-fast-math
$(AHF_STRICT_HOST): private CXXAPPFLAGS += -ffp-contract=off -fno-fast-math
endif
$(AHF_STRICT_CUDA) $(AHF_STRICT_HOST): $(wildcard AHF_GPU/*.h AHF_GPU/*.cuh)
AHF_GPU/ahf_gpu.o: AHF/patch.h AHF/fd_grid.h AHF/ghost_zone.h
AHF/expansion.o AHF/BH_diagnostics.o AHF/setup.o AHF/find_horizons.o \
AHF/initial_guess.o AHF/Newton.o AHF/IntPnts0.o AHF/IntPnts.o \
AHF_GPU/ahf_gpu.o AHF_GPU/ahf_sampling_gpu.o surface_integral_gpu.o: ShellPatch.h

# file dependences
$(C++FILES_GPU) $(F90FILES) $(AHFDOBJS) $(CUDAFILES): macrodef.fh

$(C++FILES_GPU): Block.h enforce_algebra.h fmisc.h initial_puncture.h macrodef.h\
	     misc.h monitor.h MyList.h Parallel.h MPatch.h prolongrestrict.h\
	     rungekutta4_rout.h var.h bssn_rhs.h sommerfeld_rout.h\
	     cgh.h surface_integral.h ShellPatch.h shellfunctions.h perf.h\
             fadmquantites_bssn.h cpbc.h getnp4.h initial_null.h NullEvol.h\
	     NullShellPatch.h initial_maxwell.h bssnEM_class.h getnpem2.h\
	     empart.h NullNews.h kodiss.h Parallel_bam.h ricci_gamma.h\
             initial_null2.h NullShellPatch2.h \
             bssn_gpu_class.h bssnEScalar_class.h bssn_macro.h helper.h \
			 gpu_fusion.h cpbc.h cpbc_gpu.cuh bssnEScalar_rhs_gpu.cuh
             
$(AHFDOBJS): AHF/cctk.h AHF/cctk_Config.h AHF/cctk_Types.h AHF/cctk_Constants.h AHF/myglobal.h

$(C++FILES_GPU) $(AHFDOBJS) $(CUDAFILES): macrodef.h gpu_info.h

TwoPunctureFILES: ./TwoPuncture_OMP/TwoPunctures.h

$(CUDAFILES): kodiss.h derivatives.h gpu_manager.h fmisc.h rungekutta4_rout.h sommerfeld_rout.h derivatives.h \
			  fmisc_gpu.cuh diff_new_gpu.cuh diff_new_sh_gpu.cuh kodiss_sh_gpu.cuh misc_gpu.cuh lopsidediff_gpu.cuh kodiss_gpu.cuh rhs_smem_common.cuh bssnEScalar_rhs_gpu.cuh \
			  bssn_constraint.h bssn_constraint_gpu.cuh adm_ricci_gamma_gpu.cuh bssn2adm_gpu.cuh gpu_common.cuh \
			  gpu_fusion.h cpbc.h cpbc_gpu.cuh ricci_gamma.h Z4c_rhs.h Z4c_rhs_gpu.cuh Z4c_rhs_ss_gpu.cuh point_diff_new_sh_gpu.cuh getnpem2.h getnp4.h

misc.o : zbesh.o

gpu_link.o: $(CUDAFILES)
	$(Cu) -arch=sm_$(CUDA_COMPUTE_CAPABILITY) -dlink $(CUDAFILES) -o $@

# projects
ABEGPU: $(C++FILES_GPU) $(F90FILES) $(F77FILES) $(AHFDOBJS) $(CUDAFILES) gpu_link.o
	$(CLINKER) $(CXXAPPFLAGS) -o $@ $(C++FILES_GPU) $(F90FILES) $(F77FILES) $(AHFDOBJS) $(CUDAFILES) gpu_link.o $(LDLIBS) $(AHF_GPU_LIBS)

TwoPunctureABE: $(TwoPunctureFILES)
	$(CLINKER) $(CXXAPPFLAGS) -o $@ $(TwoPunctureFILES) $(LDLIBS)

PNOrbit: $(PNOrbit_FILES)
	$(CLINKER) $(CXXAPPFLAGS) -o $@ $(PNOrbit_FILES) $(LDLIBS)

C++FILES_GPU_NOTEST = $(filter-out ABE.o,$(C++FILES_GPU))
test_em_kernels: test_em_kernels.o $(C++FILES_GPU_NOTEST) $(F90FILES) $(F77FILES) $(AHFDOBJS) $(CUDAFILES) gpu_link.o
	$(CLINKER) $(CXXAPPFLAGS) -o $@ test_em_kernels.o $(C++FILES_GPU_NOTEST) $(F90FILES) $(F77FILES) $(AHFDOBJS) $(CUDAFILES) gpu_link.o $(LDLIBS) $(AHF_GPU_LIBS)

test_full_step: test_full_step.o $(C++FILES_GPU_NOTEST) $(F90FILES) $(F77FILES) $(AHFDOBJS) $(CUDAFILES) gpu_link.o
	$(CLINKER) $(CXXAPPFLAGS) -o $@ test_full_step.o $(C++FILES_GPU_NOTEST) $(F90FILES) $(F77FILES) $(AHFDOBJS) $(CUDAFILES) gpu_link.o $(LDLIBS) $(AHF_GPU_LIBS)

clean:
	rm *.o gpu_link.o ABEGPU TwoPunctureABE PNOrbit test_em_kernels test_full_step make.log -f
	rm ./AHF/*.o ./AHF_GPU/*.o -f
	rm ./PN_Orbit/*.o -f
	rm ./TwoPuncture_OMP/*.o -f
