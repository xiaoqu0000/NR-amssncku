

include makefile.inc

.SUFFIXES: .o .f90 .C .for .cu

.f90.o:
	$(f90) $(f90appflags) -c $< -o $@

.C.o:
	${CXX} $(CXXAPPFLAGS) -c $< $(filein) -o $@

.for.o:
	$(f77) -c $< -o $@

.cu.o:
	$(Cu) $(CUDA_APP_FLAGS) -c $< -o $@ $(CUDA_LIB_PATH)

# Input files
C++FILES = ABE.o Ansorg.o Block.o misc.o monitor.o Parallel.o MPatch.o var.o\
           cgh.o bssn_class.o surface_integral.o ShellPatch.o\
	   bssnEScalar_class.o perf.o Z4c_class.o NullShellPatch.o\
	   bssnEM_class.o cpbc_util.o z4c_rhs_point.o checkpoint.o\
           Parallel_bam.o scalar_class.o transpbh.o NullShellPatch2.o\
	   NullShellPatch2_Evo.o writefile_f.o
	   
C++FILES_GPU = ABE.o Ansorg.o Block.o misc.o monitor.o Parallel.o MPatch.o var.o\
           cgh.o surface_integral.o ShellPatch.o\
	   bssnEScalar_class.o perf.o Z4c_class.o NullShellPatch.o\
	   bssnEM_class.o cpbc_util.o z4c_rhs_point.o checkpoint.o\
           Parallel_bam.o scalar_class.o transpbh.o NullShellPatch2.o\
	   NullShellPatch2_Evo.o \
	   bssn_gpu_class.o bssn_step_gpu.o bssn_macro.o writefile_f.o

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

AHFDOBJS = expansion.o expansion_Jacobian.o patch.o coords.o patch_info.o patch_interp.o patch_system.o \
tgrid.o fd_grid.o ghost_zone.o array.o round.o norm.o fuzzy.o error_exit.o miscfp.o \
linear_map.o cpm_map.o BH_diagnostics.o setup.o horizon_sequence.o find_horizons.o \
initial_guess.o Newton.o Jacobian.o ilucg.o IntPnts0.o IntPnts.o

TwoPunctureFILES = TwoPunctureABE.o TwoPunctures.o 

CUDAFILES = bssn_gpu.o bssn_gpu_rhs_ss.o

# file dependences
$(C++FILES) $(C++FILESGPU) $(F90FILES) $(AHFDOBJS) $(CUDAFILES): macrodef.fh

$(C++FILES): Block.h enforce_algebra.h fmisc.h initial_puncture.h macrodef.h\
	     misc.h monitor.h MyList.h Parallel.h MPatch.h prolongrestrict.h\
	     rungekutta4_rout.h var.h bssn_class.h bssn_rhs.h sommerfeld_rout.h\
	     cgh.h surface_integral.h ShellPatch.h shellfunctions.h perf.h\
             fadmquantites_bssn.h cpbc.h getnp4.h initial_null.h NullEvol.h\
	     NullShellPatch.h initial_maxwell.h bssnEM_class.h getnpem2.h\
	     empart.h NullNews.h kodiss.h Parallel_bam.h ricci_gamma.h\
             initial_null2.h NullShellPatch2.h 
             
$(C++FILES_GPU): Block.h enforce_algebra.h fmisc.h initial_puncture.h macrodef.h\
	     misc.h monitor.h MyList.h Parallel.h MPatch.h prolongrestrict.h\
	     rungekutta4_rout.h var.h bssn_rhs.h sommerfeld_rout.h\
	     cgh.h surface_integral.h ShellPatch.h shellfunctions.h perf.h\
             fadmquantites_bssn.h cpbc.h getnp4.h initial_null.h NullEvol.h\
	     NullShellPatch.h initial_maxwell.h bssnEM_class.h getnpem2.h\
	     empart.h NullNews.h kodiss.h Parallel_bam.h ricci_gamma.h\
             initial_null2.h NullShellPatch2.h \
             bssn_gpu_class.h bssn_macro.h
             
$(AHFDOBJS): cctk.h cctk_Config.h cctk_Types.h cctk_Constants.h myglobal.h

$(C++FILES) $(C++FILES_GPU) $(AHFDOBJS) $(CUDAFILES): macrodef.h

TwoPunctureFILES: TwoPunctures.h

$(CUDAFILES): bssn_gpu.h gpu_mem.h gpu_rhsSS_mem.h

misc.o : zbesh.o

# projects
ABE: $(C++FILES) $(F90FILES) $(F77FILES) $(AHFDOBJS) 
	$(CLINKER) $(CXXAPPFLAGS) -o $@ $(C++FILES) $(F90FILES) $(F77FILES) $(AHFDOBJS) $(LDLIBS)
	
ABEGPU: $(C++FILES_GPU) $(F90FILES) $(F77FILES) $(AHFDOBJS) $(CUDAFILES)
	$(CLINKER) $(CXXAPPFLAGS) -o $@ $(C++FILES_GPU) $(F90FILES) $(F77FILES) $(AHFDOBJS) $(CUDAFILES) $(LDLIBS)

TwoPunctureABE: $(TwoPunctureFILES)
	$(CLINKER) $(CXXAPPFLAGS) -o $@ $(TwoPunctureFILES) $(LDLIBS)

clean:
	rm *.o ABE ABEGPU TwoPunctureABE make.log -f
