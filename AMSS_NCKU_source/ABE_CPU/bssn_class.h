
#ifndef BSSN_CLASS_H
#define BSSN_CLASS_H

#ifdef newc
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <vector>
using namespace std;
#else
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#endif

#include <mpi.h>

#include "macrodef.h"
#include "cgh.h"
#include "ShellPatch.h"
#include "misc.h"
#include "var.h"
#include "MyList.h"
#include "monitor.h"
#include "surface_integral.h"
#include "checkpoint.h"
// Propagate the BSSN timing switch to the sft_timer layer.
// Must appear before worker_pool.h is included.
#ifndef BSSN_TIMING_ENABLED
#  ifndef SFT_TIMER_DISABLED
#    define SFT_TIMER_DISABLED
#  endif
#endif
#include "worker_pool.h"

#ifdef BSSN_TIMING_ENABLED

// BSSN_TIMED: wrap a block, accumulate MPI_Wtime into timing field,
// and emit a Chrome Trace event with pid=MPI-rank, tid=0.
// Call from the main (serial) thread.
#define BSSN_TIMED(name, field, ...)                                        \
  do {                                                                       \
    auto _bssn_ts0 = sft_timer::get_trace_timestamp();                       \
    double _bssn_t0 = MPI_Wtime();                                           \
    __VA_ARGS__;                                                             \
    current_iter_timing.field += MPI_Wtime() - _bssn_t0;                    \
    sft_timer::add_kernel_trace(name, _bssn_ts0,                             \
        sft_timer::get_trace_timestamp(), myrank, 0);                        \
  } while (0)

// BSSN_TRACE: pure trace only, no IterationTimingInfo update.
// For new instrumentation points that have no corresponding field.
#define BSSN_TRACE(name, ...)                                                \
  do {                                                                       \
    auto _bssn_ts0 = sft_timer::get_trace_timestamp();                       \
    __VA_ARGS__;                                                             \
    sft_timer::add_kernel_trace(name, _bssn_ts0,                             \
        sft_timer::get_trace_timestamp(), myrank, 0);                        \
  } while (0)

// BSSN_SCOPED_TRACE: RAII guard — drops a trace span for the enclosing scope.
// Usage: place at the TOP of the block you want to trace.
//   { BSSN_SCOPED_TRACE("rhs_predictor");
//     ... 40-line rhs call ...
//   }
// Works with multi-line expressions where BSSN_TIMED is too unwieldy.
struct BssnScopedTrace {
  const char* name_;
  int         pid_;
  uint64_t    ts0_;
  BssnScopedTrace(const char* name, int pid)
    : name_(name), pid_(pid), ts0_(sft_timer::get_trace_timestamp()) {}
  ~BssnScopedTrace() {
    sft_timer::add_kernel_trace(name_, ts0_, sft_timer::get_trace_timestamp(), pid_, 0);
  }
};
// Concatenation helper for unique variable names per line
#define _BSSN_CONCAT2(a,b) a##b
#define _BSSN_CONCAT(a,b)  _BSSN_CONCAT2(a,b)
#define BSSN_SCOPED_TRACE(name) \
  BssnScopedTrace _BSSN_CONCAT(_bssn_scope_, __LINE__)((name), myrank)

// BSSN_TRACE_DUR: emit a synthetic trace span ending NOW from a pre-measured
// duration (in seconds, as returned by MPI_Wtime differences).
// Used for sub-timings copied from sub-objects (e.g. Waveshell->time_xxx).
#define BSSN_TRACE_DUR(name, dur_sec)                                        \
  do {                                                                       \
    auto _dur_ts1 = sft_timer::get_trace_timestamp();                        \
    uint64_t _dur_cyc = static_cast<uint64_t>((dur_sec) *                    \
        (sft_timer::get_trace_timestamp(), 1.0) /* force init */);           \
    (void)_dur_cyc; /* see note below */                                     \
    sft_timer::add_kernel_trace(name, _dur_ts1, _dur_ts1, myrank, 0);       \
  } while (0)

// BSSN_TIMED_OMP: to be used INSIDE an already-open #pragma omp parallel region.
// Each OMP thread records its own span independently (tid = omp_get_thread_num()).
// pid should be the MPI rank so traces from different ranks stay on separate rows.
// Example:
//   #pragma omp parallel for
//   for (int i = 0; i < n; i++) {
//     BSSN_TIMED_OMP("bssn_rhs", myrank, omp_get_thread_num(), kernel(i));
//   }
#include <omp.h>
#define BSSN_TIMED_OMP(name, pid, tid, ...)                                  \
  do {                                                                       \
    auto _omp_ts0 = sft_timer::get_trace_timestamp();                        \
    __VA_ARGS__;                                                             \
    sft_timer::add_kernel_trace(name, _omp_ts0,                              \
        sft_timer::get_trace_timestamp(), (pid), (tid));                     \
  } while (0)

// Direct field helpers — use these for raw start/stop patterns that don't
// fit the BSSN_TIMED wrapper (e.g. multiple code paths writing one field).
// BSSN_FIELD_ACCUM(field, expr) : current_iter_timing.field += expr
// BSSN_FIELD_SET(field, expr)   : current_iter_timing.field  = expr
#define BSSN_FIELD_ACCUM(field, ...) (current_iter_timing.field += (__VA_ARGS__))
#define BSSN_FIELD_SET(field, ...)   (current_iter_timing.field  = (__VA_ARGS__))

#else  // !BSSN_TIMING_ENABLED — zero-overhead stubs ─────────────────────────
// Wrapped code still executes; all measurement is stripped at compile time.
#define BSSN_TIMED(name, field, ...)        do { __VA_ARGS__; } while (0)
#define BSSN_TRACE(name, ...)               do { __VA_ARGS__; } while (0)
#define BSSN_SCOPED_TRACE(name)             (void)0
#define BSSN_TRACE_DUR(name, dur_sec)       do { (void)(dur_sec); } while (0)
#define BSSN_TIMED_OMP(name, pid, tid, ...) do { __VA_ARGS__; } while (0)
#define BSSN_FIELD_ACCUM(field, ...)        do {} while (0)
#define BSSN_FIELD_SET(field, ...)          do {} while (0)
#endif  // BSSN_TIMING_ENABLED ───────────────────────────────────────────────

extern void setpbh(int iBHN, double **iPBH, double *iMass, int rBHN);

#ifdef BSSN_TIMING_ENABLED
struct IterationTimingInfo
{
       double iteration_wall;
       double evolve_kernel;
       double analysis;
       double constraint;
       double dump_binary;
       double dump_2d;
       double checkpoint;
       double regrid;
       double memory_probe;
       double abort_check;
       double kernel_step;
       double kernel_enforce;
       double kernel_restrict;
       double kernel_restrict_prepare;    // prepare_inter_time_level (time-refine branch only)
       double kernel_restrict_restrict;   // Restrict / Restrict_bam
       double kernel_restrict_sync_coarse;// Sync on coarse level (lev-1)
       double kernel_restrict_outbd;      // OutBdLow2Hi / OutBdLow2Hi_bam
       double kernel_restrict_sync_fine;  // final Sync on fine level (lev)
       double kernel_regrid_local;
       double kernel_step_rhs;
       double kernel_step_rhs_predictor; // Added for detailed timing
       double kernel_step_rhs_corrector; // Added for detailed timing
       double kernel_step_rungekutta;
       double kernel_step_boundary;
       double kernel_step_sync;
       double kernel_step_shell_rhs;
       double kernel_step_shell_rungekutta;
       double kernel_step_shell_boundary;
       double kernel_step_shell_sync;
       double kernel_step_bh;
       double kernel_step_analysis;
       double kernel_step_analysis_waveshell; // Added for detailed timing
       double kernel_step_analysis_fderivs; // Added for detailed timing
       double kernel_step_analysis_psi4; // Added for detailed timing
       double kernel_step_analysis_waveshell_surf_wave;
       double kernel_step_analysis_waveshell_surf_masspang;
       double kernel_step_analysis_waveshell_interp;
       double kernel_step_analysis_waveshell_int;
       double kernel_step_analysis_waveshell_comm;
       double kernel_step_analysis_waveshell_prepare;
       // Fine-grained timing for "other" breakdown
       double kernel_step_enforce_ga;
       double kernel_step_lowerbound;
       double kernel_step_error_check;
       double kernel_step_loop_overhead;
};
#else
struct IterationTimingInfo {};
#endif  // BSSN_TIMING_ENABLED

class bssn_class
{
public:
       int ngfs;
       int nprocs, myrank;
       cgh *GH;
       ShellPatch *SH;
       double PhysTime;

       int checkrun;
       char checkfilename[50];
       int Steps;
       double StartTime, TotalTime;
       double AnasTime, DumpTime, d2DumpTime, CheckTime;
       double LastAnas, LastConsOut;
       double Courant;
       double numepss, numepsb, numepsh;
       int Symmetry;
       int maxl, decn;
       double maxrex, drex;
       int trfls, a_lev;

       double dT;
       double chitiny;

       double **Porg0, **Porgbr, **Porg, **Porg1, **Porg_rhs;
       int BH_num, BH_num_input;
       double *Mass, *Pmom, *Spin;
       double ADMMass;

       var *phio, *trKo;
       var *gxxo, *gxyo, *gxzo, *gyyo, *gyzo, *gzzo;
       var *Axxo, *Axyo, *Axzo, *Ayyo, *Ayzo, *Azzo;
       var *Gmxo, *Gmyo, *Gmzo;
       var *Lapo, *Sfxo, *Sfyo, *Sfzo;
       var *dtSfxo, *dtSfyo, *dtSfzo;

       var *phi0, *trK0;
       var *gxx0, *gxy0, *gxz0, *gyy0, *gyz0, *gzz0;
       var *Axx0, *Axy0, *Axz0, *Ayy0, *Ayz0, *Azz0;
       var *Gmx0, *Gmy0, *Gmz0;
       var *Lap0, *Sfx0, *Sfy0, *Sfz0;
       var *dtSfx0, *dtSfy0, *dtSfz0;

       var *phi, *trK;
       var *gxx, *gxy, *gxz, *gyy, *gyz, *gzz;
       var *Axx, *Axy, *Axz, *Ayy, *Ayz, *Azz;
       var *Gmx, *Gmy, *Gmz;
       var *Lap, *Sfx, *Sfy, *Sfz;
       var *dtSfx, *dtSfy, *dtSfz;

       var *phi1, *trK1;
       var *gxx1, *gxy1, *gxz1, *gyy1, *gyz1, *gzz1;
       var *Axx1, *Axy1, *Axz1, *Ayy1, *Ayz1, *Azz1;
       var *Gmx1, *Gmy1, *Gmz1;
       var *Lap1, *Sfx1, *Sfy1, *Sfz1;
       var *dtSfx1, *dtSfy1, *dtSfz1;

       var *phi_rhs, *trK_rhs;
       var *gxx_rhs, *gxy_rhs, *gxz_rhs, *gyy_rhs, *gyz_rhs, *gzz_rhs;
       var *Axx_rhs, *Axy_rhs, *Axz_rhs, *Ayy_rhs, *Ayz_rhs, *Azz_rhs;
       var *Gmx_rhs, *Gmy_rhs, *Gmz_rhs;
       var *Lap_rhs, *Sfx_rhs, *Sfy_rhs, *Sfz_rhs;
       var *dtSfx_rhs, *dtSfy_rhs, *dtSfz_rhs;

       var *rho, *Sx, *Sy, *Sz, *Sxx, *Sxy, *Sxz, *Syy, *Syz, *Szz;

       var *Gamxxx, *Gamxxy, *Gamxxz, *Gamxyy, *Gamxyz, *Gamxzz;
       var *Gamyxx, *Gamyxy, *Gamyxz, *Gamyyy, *Gamyyz, *Gamyzz;
       var *Gamzxx, *Gamzxy, *Gamzxz, *Gamzyy, *Gamzyz, *Gamzzz;

       var *Rxx, *Rxy, *Rxz, *Ryy, *Ryz, *Rzz;

       var *Rpsi4, *Ipsi4;
       var *t1Rpsi4, *t1Ipsi4, *t2Rpsi4, *t2Ipsi4;

       var *Cons_Ham, *Cons_Px, *Cons_Py, *Cons_Pz, *Cons_Gx, *Cons_Gy, *Cons_Gz;

#ifdef Point_Psi4
       var *phix, *phiy, *phiz;
       var *trKx, *trKy, *trKz;
       var *Axxx, *Axxy, *Axxz;
       var *Axyx, *Axyy, *Axyz;
       var *Axzx, *Axzy, *Axzz;
       var *Ayyx, *Ayyy, *Ayyz;
       var *Ayzx, *Ayzy, *Ayzz;
       var *Azzx, *Azzy, *Azzz;
#endif
       // FIXME: uc = StateList, up = OldStateList, upp = SynchList_cor; so never touch these three data
       MyList<var> *StateList, *SynchList_pre, *SynchList_cor, *RHSList;
       MyList<var> *OldStateList, *DumpList;
       MyList<var> *ConstraintList;

       monitor *ErrorMonitor, *Psi4Monitor, *BHMonitor, *MAPMonitor;
       monitor *ConVMonitor;
       surface_integral *Waveshell;
       checkpoint *CheckPoint;

       IterationTimingInfo current_iter_timing;
       IterationTimingInfo total_iter_timing;
       int completed_iterations;
       bool total_timing_printed;
public:
       bssn_class(double Couranti, double StartTimei, double TotalTimei, double DumpTimei, double d2DumpTimei, double CheckTimei, double AnasTimei,
                  int Symmetryi, int checkruni, char *checkfilenamei, double numepssi, double numepsbi, double numepshi,
                  int a_levi, int maxli, int decni, double maxrexi, double drexi);
       ~bssn_class();

       void Evolve(int Steps);
       void RecursiveStep(int lev);
#if (PSTR == 3)
       void RecursiveStep(int lev, int num);
#endif
#if (PSTR == 1 || PSTR == 2 || PSTR == 3)
       void ParallelStep();
       void SHStep();
#endif
       void RestrictProlong(int lev, int YN, bool BB, MyList<var> *SL, MyList<var> *OL, MyList<var> *corL);
       void RestrictProlong_aux(int lev, int YN, bool BB, MyList<var> *SL, MyList<var> *OL, MyList<var> *corL);
       void RestrictProlong(int lev, int YN, bool BB);
       void ProlongRestrict(int lev, int YN, bool BB);
       void Setup_Black_Hole_position();
       void compute_Porg_rhs(double **BH_PS, double **BH_RHS, var *forx, var *fory, var *forz, int lev);
       bool read_Pablo_file(int *ext, double *datain, char *filename);
       void write_Pablo_file(int *ext, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax,
                             char *filename);
       void AnalysisStuff(int lev, double dT_lev);
       void Setup_KerrSchild();
       void Enforce_algcon(int lev, int fg);

       void testRestrict();
       void testOutBd();
       
       bool check_Stdin_Stop(); 

       virtual void Setup_Initial_Data_Cao();
       virtual void Setup_Initial_Data_Lousto();
       virtual void Initialize();
       virtual void Read_Ansorg();
       virtual void Read_Pablo() {};
       virtual void Compute_Psi4(int lev);
       virtual void Step(int lev, int YN);
       virtual void Interp_Constraint(bool infg);
       virtual void Constraint_Out();
       virtual void Compute_Constraint();

       void reset_iteration_timing();
       void finalize_iteration_timing(int stepIndex, double physicalTime, bool isTerminal, const char *reason);
       void print_total_timing_summary(const char *reason);

#ifdef With_AHF
protected:
       MyList<var> *AHList, *AHDList, *GaugeList;
       int AHfindevery;
       double AHdumptime;
       int *lastahdumpid, HN_num; // number of possible horizons
       int *findeveryl;
       double *xc, *yc, *zc, *xr, *yr, *zr;
       bool *trigger;
       double *dTT;
       int *dumpid;

public:
       void AH_Prepare_derivatives();
       bool AH_Interp_Points(MyList<var> *VarList,
                             int NN, double **XX,
                             double *Shellf, int Symmetryi);
       void AH_Step_Find(int lev, double dT_lev);
#endif
};
#endif /* BSSN_CLASS_H */
