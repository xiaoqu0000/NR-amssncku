

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
C++FILES = ABE.o \
           ./Initial_Data_Solver/Ansorg.o \
           ./cgh/Block.o \
           ./misc/misc.o \
           ./Monitor/monitor.o \
           ./Parallel/Parallel.o \
           ./Patch/MPatch.o \
           ./Variable/var.o \
           ./cgh/cgh.o \
           ./BSSN/bssn_class.o \
           ./Surface_Integral/surface_integral.o \
           ./Shell_Patch/ShellPatch.o \
	   ./Scalar/bssnEScalar_class.o \
	   ./System_Performance/perf.o \
	   ./Z4C/Z4c_class.o \
	   ./Null_Evolve/NullShellPatch.o \
	   ./BSSN/bssnEM_class.o \
	   ./Z4C/cpbc_util.o \
	   ./Z4C/z4c_rhs_point.o \
	   ./Check_Point/checkpoint.o \
           ./Parallel/Parallel_bam.o \
           ./Scalar/scalar_class.o \
           ./BSSN/transpbh.o \
           ./Null_Evolve/NullShellPatch2.o \
	   ./Null_Evolve/NullShellPatch2_Evo.o \
	   ./Read_and_Write/writefile_f.o
	   
C++FILES_GPU = ABE.o \
               ./Initial_Data_Solver/Ansorg.o \
               ./cgh/Block.o \
               ./misc/misc.o \
               ./Monitor/monitor.o \
               ./Parallel/Parallel.o \
               ./Patch/MPatch.o \
               ./Variable/var.o \
               ./cgh/cgh.o \
               ./Surface_Integral/surface_integral.o \
               ./Shell_Patch/ShellPatch.o \
	       ./Scalar/bssnEScalar_class.o \
	       ./System_Performance/perf.o \
	       ./Z4C/Z4c_class.o \
	       ./Null_Evolve/NullShellPatch.o \
	       ./BSSN/bssnEM_class.o \
	       ./Z4C/cpbc_util.o \
	       ./Z4C/z4c_rhs_point.o \
	       ./Check_Point/checkpoint.o \
               ./Parallel/Parallel_bam.o \
               ./Scalar/scalar_class.o \
               ./BSSN/transpbh.o \
               ./Null_Evolve/NullShellPatch2.o \
	       ./Null_Evolve/NullShellPatch2_Evo.o \
	       ./BSSN_GPU/bssn_gpu_class.o \
	       ./BSSN_GPU/bssn_step_gpu.o \
	       ./BSSN_GPU/bssn_macro.o \
	       ./Read_and_Write/writefile_f.o

F90FILES = ./BSSN/enforce_algebra.o \
           ./misc/fmisc.o \
           ./Initial_Data_Solver/initial_puncture.o \
           ./BSSN/prolongrestrict.o \
	   ./BSSN/prolongrestrict_cell.o \
	   ./BSSN/prolongrestrict_vertex.o \
	   ./Runge_Kutta/rungekutta4_rout.o \
	   ./BSSN/bssn_rhs.o \
	   ./Derivative/diff_new.o \
	   ./KO_dissipation/kodiss.o \
	   ./KO_dissipation/kodiss_sh.o \
	   ./BSSN/lopsidediff.o \
	   ./BSSN/sommerfeld_rout.o \
	   ./Psi4/getnp4.o \
	   ./Derivative/diff_new_sh.o \
	   ./Shell_Patch/shellfunctions.o \
	   ./BSSN/bssn_rhs_ss.o \
	   ./Scalar/Set_Rho_ADM.o \
           ./Psi4/getnp4EScalar.o \
           ./Scalar/bssnEScalar_rhs.o \
           ./BSSN/bssn_constraint.o \
           ./Psi4/ricci_gamma.o \
           ./BSSN/fadmquantites_bssn.o \
           ./Z4C/Z4c_rhs.o \
           ./Z4C/Z4c_rhs_ss.o \
           ./Derivative/point_diff_new_sh.o \
	   ./Z4C/cpbc.o \
	   ./Psi4/getnp4old.o \
	   ./Null_Evolve/NullEvol.o \
	   ./Initial_Data_Solver/initial_null.o \
	   ./Initial_Data_Solver/initial_maxwell.o\
	   ./Psi4/getnpem2.o \
	   ./BSSN/empart.o \
	   ./Null_Evolve/NullNews.o \
	   ./BSSN/fourdcurvature.o \
	   ./BSSN/bssn2adm.o \
	   ./BSSN/adm_constraint.o \
	   ./Psi4/adm_ricci_gamma.o\
	   ./Scalar/scalar_rhs.o \
	   ./Initial_Data_Solver/initial_scalar.o \
	   ./Null_Evolve/NullEvol2.o \
	   ./Initial_Data_Solver/initial_null2.o \
	   ./Null_Evolve/NullNews2.o \
	   ./Read_and_Write/tool_f.o

F77FILES = ./Special_Function/zbesh.o 

CUDAFILES = ./BSSN_GPU/bssn_gpu.o \
            ./BSSN_GPU/bssn_gpu_rhs_ss.o

AHFDirect = ./AHF_Direct/expansion.o  ./AHF_Direct/expansion_Jacobian.o  \
            ./AHF_Direct/patch.o      ./AHF_Direct/coords.o  \
            ./AHF_Direct/patch_info.o ./AHF_Direct/patch_interp.o  ./AHF_Direct/patch_system.o \
            ./AHF_Direct/tgrid.o      ./AHF_Direct/fd_grid.o \
            ./AHF_Direct/ghost_zone.o ./AHF_Direct/array.o  \
            ./AHF_Direct/round.o      ./AHF_Direct/norm.o   \
            ./AHF_Direct/fuzzy.o      ./AHF_Direct/error_exit.o \
            ./AHF_Direct/miscfp.o     ./AHF_Direct/linear_map.o \
            ./AHF_Direct/cpm_map.o    ./AHF_Direct/BH_diagnostics.o \
            ./AHF_Direct/setup.o      ./AHF_Direct/horizon_sequence.o ./AHF_Direct/find_horizons.o \
            ./AHF_Direct/initial_guess.o \
            ./AHF_Direct/Newton.o     ./AHF_Direct/Jacobian.o \
            ./AHF_Direct/ilucg.o      ./AHF_Direct/IntPnts0.o ./AHF_Direct/IntPnts.o

TwoPunctureFILES = ./Two_Puncture/TwoPunctureABE.o ./Two_Puncture/TwoPunctures.o 

# file dependences
$(C++FILES) $(C++FILESGPU) $(F90FILES) $(AHFDIRECT_OBJS) $(CUDAFILES): macrodef.fh

$(C++FILES): ./cgh/Block.h \
             ./BSSN/enforce_algebra.h \
             ./misc/fmisc.h \
             ./Initial_Data_Solver/initial_puncture.h \
             macrodef.h \
	     ./misc/misc.h \
	     ./Monitor/monitor.h \
	     ./Variable/MyList.h \
	     ./Variable/parameters.h \
	     ./Parallel/Parallel.h \
	     ./Patch/MPatch.h \
	     ./BSSN/prolongrestrict.h \
	     ./Runge_Kutta/rungekutta4_rout.h \
	     ./Variable/var.h \
	     ./BSSN/bssn_class.h \
	     ./BSSN/bssn_rhs.h \
	     ./BSSN/sommerfeld_rout.h \
	     ./cgh/cgh.h \
	     ./Surface_Integral/surface_integral.h \
	     ./Shell_Patch/ShellPatch.h \
	     ./Shell_Patch/shellfunctions.h \
	     ./System_Performance/perf.h \
             ./BSSN/fadmquantites_bssn.h \
             ./Z4C/cpbc.h \
             ./Psi4/getnp4.h \
             ./Initial_Data_Solver/initial_null.h \
             ./Null_Evolve/NullEvol.h \
	     ./Null_Evolve/NullShellPatch.h \
	     ./Initial_Data_Solver/initial_maxwell.h \
	     ./BSSN/bssnEM_class.h \
	     ./Psi4/getnpem2.h \
	     ./BSSN/empart.h \
	     ./Null_Evolve/NullNews.h \
	     ./KO_dissipation/kodiss.h \
	     ./Parallel/Parallel_bam.h \
	     ./Psi4/ricci_gamma.h\
             ./Initial_Data_Solver/initial_null2.h \
             ./Null_Evolve/NullShellPatch2.h 
             
$(C++FILES_GPU): ./cgh/Block.h \
                 ./BSSN/enforce_algebra.h \
                 ./misc/fmisc.h \
                 ./Initial_Data_Solver/initial_puncture.h \
                 macrodef.h \
	         ./misc/misc.h \
	         ./Monitor/monitor.h \
	         ./Variable/MyList.h \
	         ./Variable/parameters.h \
	         ./Parallel/Parallel.h \
	         ./Patch/MPatch.h \
	         ./BSSN/prolongrestrict.h \
	         ./Runge_Kutta/rungekutta4_rout.h \
	         ./Variable/var.h \
	         ./BSSN/bssn_rhs.h \
	         ./BSSN/sommerfeld_rout.h\
	         ./cgh/cgh.h \
	         ./Surface_Integral/surface_integral.h \
	         ./Shell_Patch/ShellPatch.h \
	         ./Shell_Patch/shellfunctions.h \
	         ./System_Performance/perf.h \
                 ./BSSN/fadmquantites_bssn.h \
                 ./Z4C/cpbc.h \
                 ./Psi4/getnp4.h \
                 ./Initial_Data_Solver/initial_null.h \
                 ./Null_Evolve/NullEvol.h\
	         ./Null_Evolve/NullShellPatch.h \
	         ./Initial_Data_Solver/initial_maxwell.h \
	         ./BSSN/bssnEM_class.h \
	         ./Psi4/getnpem2.h \
	         ./BSSN/empart.h \
	         ./Null_Evolve/NullNews.h \
	         ./KO_dissipation/kodiss.h \
	         ./Parallel/Parallel_bam.h \
	         ./Psi4/ricci_gamma.h \
                 ./Initial_Data_Solver/initial_null2.h \
                 ./Null_Evolve/NullShellPatch2.h \
                 ./BSSN_GPU/bssn_gpu_class.h \
                 ./BSSN_GPU/bssn_macro.h
              
$(BSSN_GPU_OBJS): ./BSSN_GPU/bssn_gpu_class.h
             
$(Z4C_OBJS): ./Z4C/Z4c_class.h ./Z4C/cpbc.h 
             
$(ID_OBJS): ./Initial_Data_Solver/Ansorg.h       ./Initial_Data_Solver/initial_puncture.h ./Initial_Data_Solver/initial_maxwell.h\
            ./Initial_Data_Solver/initial_null.h ./Initial_Data_Solver/initial_null2.h
            
$(NULL_OBJS): ./Null_Evolve/NullEvol.h       ./Null_Evolve/NullNews.h \
	      ./Null_Evolve/NullShellPatch.h ./Null_Evolve/NullShellPatch2.h
             
$(AHFDirect): ./AHF_Direct/array.h             ./AHF_Direct/cctk.h       ./AHF_Direct/config.h      ./AHF_Direct/cctk_Types.h      ./AHF_Direct/cctk_Constants.h\
              ./AHF_Direct/cctk_Config.h       ./AHF_Direct/coords.h     ./AHF_Direct/cpm_map.h     ./AHF_Direct/driver.h\
              ./AHF_Direct/fd_grid.h           ./AHF_Direct/gfns.h       ./AHF_Direct/ghost_zone.h  ./AHF_Direct/gr.h\
              ./AHF_Direct/horizon_sequence.h  ./AHF_Direct/ilucg.h      ./AHF_Direct/linear_map.h  ./AHF_Direct/BH_diagnostics.h  ./AHF_Direct/myglobal.h\
              ./AHF_Direct/patch.h             ./AHF_Direct/patch_edge.h ./AHF_Direct/patch_info.h  ./AHF_Direct/patch_interp.h    ./AHF_Direct/patch_system.h\
              ./AHF_Direct/stdc.h              ./AHF_Direct/tgrid.h      ./AHF_Direct/util.h        ./AHF_Direct/util_String.h     ./AHF_Direct/util_Table.h

$(C++FILES) $(CUDAFILES) $(AHFDirect): macrodef.h

TwoPunctureFILES: ./Two_Puncture/TwoPunctures.h

$(CUDAFILES): ./BSSN_GPU/bssn_gpu.h  ./BSSN_GPU/gpu_mem.h  ./BSSN_GPU/gpu_rhsSS_mem.h

./misc/misc.o: ./Special_Function/zbesh.o

# projects
ABE: ABE.o $(C++FILES) $(F90FILES) $(F77FILES) $(AHFDirect) 
	$(CLINKER) $(CXXAPPFLAGS) -o $@  $(C++FILES) $(F90FILES) $(F77FILES) $(AHFDirect) $(LDLIBS)
	
ABEGPU: $(C++FILES_GPU) $(F90FILES) $(F77FILES) $(CUDAFILES) $(AHFDirect) 
	$(CLINKER) $(CXXAPPFLAGS) -o $@ $(C++FILES_GPU) $(F90FILES) $(F77FILES) $(CUDAFILES) $(AHFDirect) $(LDLIBS)

TwoPunctureABE: $(TwoPunctureFILES)
	$(CLINKER) $(CXXAPPFLAGS) -o $@ $(TwoPunctureFILES) $(LDLIBS)

clean:
	rm *.o ABE ABEGPU TwoPunctureABE make.log -f
