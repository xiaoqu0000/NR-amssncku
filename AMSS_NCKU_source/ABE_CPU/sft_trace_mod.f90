!> @file sft_trace_mod.f90
!> @brief Fortran ISO_C_BINDING interface to sft_timer Chrome-Trace API.
!>
!> Allows Fortran OpenMP kernels to emit per-thread trace events visible
!> alongside C++ events in Perfetto / chrome://tracing.
!>
!> Usage pattern inside an OMP region:
!>
!>   use sft_trace_mod
!>   use omp_lib
!>   integer(8) :: ts0
!>   integer    :: tid
!>
!>   !$OMP PARALLEL PRIVATE(ts0, tid)
!>     tid = omp_get_thread_num()
!>     ts0 = sft_get_ts()
!>     !$OMP DO PRIVATE(i,j,k) COLLAPSE(1)
!>       ...kernel...
!>     !$OMP END DO
!>     call sft_add_trace("bssn_rhs_opt", ts0, sft_get_ts(), mpi_rank, tid)
!>   !$OMP END PARALLEL
!>
!> pid (mpi_rank) should be passed in from the calling routine so that
!> different MPI processes appear as separate rows in Perfetto.
!>
!> When BSSN_TIMING_ENABLED is NOT defined (default), all routines become
!> compile-time no-ops: sft_get_ts() returns the constant 0 and sft_trace()
!> is an empty subroutine that the compiler will inline and eliminate.

module sft_trace_mod
  use iso_c_binding
  implicit none

#ifdef BSSN_TIMING_ENABLED
  ! ── enabled: real ISO_C_BINDING interface to C++ sft_timer ───────────────
  interface

    ! Get the current RDTSC timestamp.
    ! Call before and after a computation to bracket it.
    function sft_get_ts() result(ts) bind(C, name="sft_c_get_ts")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t) :: ts
    end function sft_get_ts

    ! Add one trace event for the span [ts0, ts1].
    ! name     : label string (Fortran passes fixed-length; trailing spaces are trimmed)
    ! name_len : len_trim(name) - pass explicitly to avoid zero-termination issues
    ! ts0, ts1 : timestamps from sft_get_ts()
    ! pid      : MPI rank (becomes the "process" row in Perfetto)
    ! tid      : OMP thread id (becomes the "thread" row within the process)
    subroutine sft_add_trace(name, name_len, ts0, ts1, pid, tid) &
        bind(C, name="sft_c_add_trace")
      use iso_c_binding, only: c_char, c_int, c_int64_t
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int),   value, intent(in) :: name_len
      integer(c_int64_t), value, intent(in) :: ts0, ts1
      integer(c_int),   value, intent(in) :: pid, tid
    end subroutine sft_add_trace

  end interface
#endif  /* BSSN_TIMING_ENABLED */

contains

#ifdef BSSN_TIMING_ENABLED
  ! ── enabled: wrapper that computes len_trim and appends null char ─────────
  !> Convenience wrapper: automatically computes len_trim and adds null char.
  !> Usage: call sft_trace("bssn_rhs_opt", ts0, sft_get_ts(), mpi_rank, tid)
  subroutine sft_trace(name, ts0, ts1, pid, tid)
    character(len=*), intent(in)    :: name
    integer(c_int64_t), intent(in)  :: ts0, ts1
    integer(c_int), intent(in)      :: pid, tid

    ! Convert Fortran string to C-compatible call
    call sft_add_trace(name // c_null_char, int(len_trim(name), c_int), &
                       ts0, ts1, pid, tid)
  end subroutine sft_trace

#else
  ! ── disabled: no-op stubs ── BSSN_TIMING_ENABLED not set ─────────────────
  ! sft_get_ts() returns the compile-time constant 0.
  ! The compiler will propagate ts0 = 0 and eliminate the sft_trace() call.
  pure function sft_get_ts() result(ts)
    integer(c_int64_t) :: ts
    ts = 0_c_int64_t
  end function sft_get_ts

  subroutine sft_trace(name, ts0, ts1, pid, tid)
    character(len=*), intent(in)    :: name
    integer(c_int64_t), intent(in)  :: ts0, ts1
    integer(c_int), intent(in)      :: pid, tid
    ! Intentionally empty — inlined and eliminated by the compiler.
  end subroutine sft_trace
#endif  /* BSSN_TIMING_ENABLED */

end module sft_trace_mod
