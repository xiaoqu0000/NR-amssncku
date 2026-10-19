

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
	   NullShellPatch2_Evo.o writefile_f.o\
	   worker_pool.o interp_shm.o

F90FILES = sft_trace_mod.o enforce_algebra.o fmisc.o fmisc_interp.o initial_puncture.o prolongrestrict.o\
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

# file dependences
$(C++FILES) $(F90FILES) $(AHFDOBJS): macrodef.fh

$(C++FILES): Block.h enforce_algebra.h fmisc.h initial_puncture.h macrodef.h\
	     misc.h monitor.h MyList.h Parallel.h MPatch.h prolongrestrict.h\
	     rungekutta4_rout.h var.h bssn_class.h bssn_rhs.h sommerfeld_rout.h\
	     cgh.h surface_integral.h ShellPatch.h shellfunctions.h perf.h\
             fadmquantites_bssn.h cpbc.h getnp4.h initial_null.h NullEvol.h\
	     NullShellPatch.h initial_maxwell.h bssnEM_class.h getnpem2.h\
	     empart.h NullNews.h kodiss.h Parallel_bam.h ricci_gamma.h\
             initial_null2.h NullShellPatch2.h 
             
$(AHFDOBJS): AHF/cctk.h AHF/cctk_Config.h AHF/cctk_Types.h AHF/cctk_Constants.h AHF/myglobal.h

$(C++FILES) $(AHFDOBJS): macrodef.h

TwoPunctureFILES: ./TwoPuncture_OMP/TwoPunctures.h

# fmisc_interp.f90 含 interp 小循环代码，用 precise 浮点模型编译（不用 fast=2）
fmisc_interp.o: fmisc_interp.f90
	$(f90) $(f90flags_interp) -c $< -o $@

# Fortran module dependencies (module must be compiled before its users)
sft_trace_mod.o: sft_trace_mod.f90

bssn_rhs.o: sft_trace_mod.o
rungekutta4_rout.o: sft_trace_mod.o
diff_new.o: sft_trace_mod.o
enforce_algebra.o: sft_trace_mod.o
lopsidediff.o: sft_trace_mod.o
fmisc.o: sft_trace_mod.o
kodiss.o: sft_trace_mod.o
sommerfeld_rout.o: sft_trace_mod.o

misc.o : zbesh.o

# projects
ABE: $(C++FILES) $(F90FILES) $(F77FILES) $(AHFDOBJS) 
	$(CLINKER) $(CXXAPPFLAGS) -o $@ $(C++FILES) $(F90FILES) $(F77FILES) $(AHFDOBJS) $(LDLIBS)

TwoPunctureABE: $(TwoPunctureFILES)
	$(CLINKER) $(CXXAPPFLAGS) -o $@ $(TwoPunctureFILES) $(LDLIBS)
	
PNOrbit: $(PNOrbit_FILES)
	$(CLINKER) $(CXXAPPFLAGS) -o $@ $(PNOrbit_FILES) $(LDLIBS)

clean:
	rm *.o ABE TwoPunctureABE PNOrbit make.log -f
	rm ./AHF/*.o -f
	rm ./PN_Orbit/*.o -f
	rm ./TwoPuncture_OMP/*.o -f
