

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
	   bssnEScalar_gpu_class.o perf.o Z4c_gpu_class.o NullShellPatch.o\
	   bssnEM_gpu_class.o cpbc_util.o z4c_rhs_point.o checkpoint.o\
           Parallel_bam.o transpbh.o NullShellPatch2.o\
	   NullShellPatch2_Evo.o \
	   bssn_gpu_class.o bssn_step_gpu.o bssn_macro.o writefile_f.o\
	   Parallel_GPU.o data_motion_helper.o

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
AHF/initial_guess.o AHF/Newton.o AHF/Jacobian.o AHF/ilucg.o AHF/IntPnts0.o AHF/IntPnts.o

TwoPunctureFILES = ./TwoPuncture_OMP/TwoPunctureABE.o ./TwoPuncture_OMP/TwoPunctures.o 

PNOrbit_FILES = ./PN_Orbit/PN_orbit.o ./PN_Orbit/misc.o ./PN_Orbit/lagpolint.o

CUDAFILES = fmisc_gpu.o bssn_rhs_gpu.o bssn_rhs_shell_gpu.o prolongrestrict_cell_gpu.o prolongrestrict_vertex_gpu.o gpu_runtime_coordinator.o enforce_algebra_gpu.o sommerfeld_rout_gpu.o rungekutta4_rout_gpu.o fadmquantites_bssn_gpu.o MPatch_gpu.o surface_integral_gpu.o getnp4_gpu.o getnpem2_gpu.o initial_data_gpu.o shell_patch_gpu.o ahf_derivatives_gpu.o ahf_gpu.o

# file dependences
$(C++FILES_GPU) $(F90FILES) $(AHFDOBJS) $(CUDAFILES): macrodef.fh

$(C++FILES_GPU): Block.h enforce_algebra.h fmisc.h initial_puncture.h macrodef.h\
	     misc.h monitor.h MyList.h Parallel.h MPatch.h prolongrestrict.h\
	     rungekutta4_rout.h var.h bssn_rhs.h sommerfeld_rout.h\
	     cgh.h surface_integral.h ShellPatch.h shellfunctions.h perf.h\
             fadmquantites_bssn.h cpbc.h getnp4.h initial_null.h NullEvol.h\
	     NullShellPatch.h initial_maxwell.h bssnEM_gpu_class.h getnpem2.h\
	     empart.h NullNews.h kodiss.h Parallel_bam.h ricci_gamma.h\
             initial_null2.h NullShellPatch2.h \
	             bssn_gpu_class.h bssn_macro.h data_motion_helper.h initial_data_gpu.h
             
$(AHFDOBJS): AHF/cctk.h AHF/cctk_Config.h AHF/cctk_Types.h AHF/cctk_Constants.h AHF/myglobal.h

$(C++FILES_GPU) $(AHFDOBJS) $(CUDAFILES): macrodef.h

TwoPunctureFILES: ./TwoPuncture_OMP/TwoPunctures.h

$(CUDAFILES): kodiss.h derivatives.h gpu_runtime_coordinator.h fmisc.h rungekutta4_rout.h sommerfeld_rout.h derivatives.h \
				  fmisc_gpu.cuh diff_new_gpu.cuh misc_gpu.cuh lopsidediff_gpu.cuh kodiss_gpu.cuh bssn_rhs_gpu.cuh bssn_rhs_shell_gpu.cuh getnpem2_gpu.cuh initial_data_gpu.h shell_patch_gpu.h AHF_GPU/ahf_derivatives_gpu.h

ahf_gpu.o: AHF_GPU/ahf_gpu.cu AHF_GPU/ahf_gpu.h
	$(Cu) $(CUDA_APP_FLAGS) -c "$<" -o $@ $(CUDA_LIB_PATH)

ahf_derivatives_gpu.o: AHF_GPU/ahf_derivatives_gpu.cu AHF_GPU/ahf_derivatives_gpu.h
	$(Cu) $(CUDA_APP_FLAGS) -c "$<" -o $@ $(CUDA_LIB_PATH)

misc.o : zbesh.o

bssn_rhs_shell_gpu.o: bssn_rhs_shell_gpu.cu bssn_rhs_shell_gpu.cuh bssn_rhs.h
	$(Cu) $(CUDA_APP_FLAGS) --fmad=false -c $< -o $@ $(CUDA_LIB_PATH)

gpu_link.o: $(CUDAFILES)
	$(Cu) -arch=sm_$(CUDA_COMPUTE_CAPABILITY) -dlink $(CUDAFILES) -o $@

# projects
ABEGPU: $(C++FILES_GPU) $(F90FILES) $(F77FILES) $(AHFDOBJS) $(CUDAFILES) gpu_link.o
	$(CLINKER) $(CXXAPPFLAGS) -o $@ $(C++FILES_GPU) $(F90FILES) $(F77FILES) $(AHFDOBJS) $(CUDAFILES) gpu_link.o $(LDLIBS)

TwoPunctureABE: $(TwoPunctureFILES)
	$(CLINKER) $(CXXAPPFLAGS) -o $@ $(TwoPunctureFILES) $(LDLIBS)

PNOrbit: $(PNOrbit_FILES)
	$(CLINKER) $(CXXAPPFLAGS) -o $@ $(PNOrbit_FILES) $(LDLIBS)

clean:
	rm *.o gpu_link.o ABEGPU TwoPunctureABE PNOrbit make.log -f
	rm ./AHF/*.o -f
	rm ./PN_Orbit/*.o -f
	rm ./TwoPuncture_OMP/*.o -f
