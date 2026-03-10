module io_logger

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Centralized I/O tracing module for recording subroutine inputs/outputs.
  ! All logging calls in the codebase are gated by #ifdef IO_TRACE so there
  ! is zero overhead when IO_TRACE is not defined at compile time.
  !
  ! Directory layout produced at runtime:
  !   io_trace/
  !     call_000001__SoilAlbedo/
  !       inputs.txt
  !       outputs.txt
  !     call_000002__MLCanopyFluxes/
  !       inputs.txt
  !       outputs.txt
  !     ...
  !
  ! A global call counter ensures that nested / sequential calls are
  ! numbered in chronological order.
  !-----------------------------------------------------------------------

  use shr_kind_mod,  only : r8 => shr_kind_r8
  use iso_c_binding, only : c_int, c_char, c_null_char
  implicit none
  private

  ! ---- module state ----
  integer, save          :: global_call_counter = 0
  character(len=256), save :: trace_base_dir    = 'io_trace'
  logical, save          :: trace_initialized   = .false.

  ! ---- per-subroutine call limiting ----
  integer, parameter     :: max_traces_per_sub = 2
  integer, parameter     :: max_tracked_subs   = 128
  integer, save          :: num_tracked_subs   = 0
  character(len=64), save :: tracked_sub_names(128)
  integer, save          :: tracked_sub_counts(128) = 0

  ! ---- POSIX mkdir via C interop (avoids fork-per-call overhead) ----
  interface
    integer(c_int) function c_mkdir(path, mode) bind(C, name="mkdir")
      import :: c_int, c_char
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int), value :: mode
    end function c_mkdir
  end interface

  ! ---- public API ----
  public :: io_trace_init
  public :: io_trace_begin
  public :: io_trace_end
  public :: io_trace_open_stage
  public :: io_trace_close_stage

  ! scalars
  public :: log_r8, log_int, log_char

  ! assumed-shape arrays (for pointer members, assumed-shape args, etc.)
  public :: log_r8_arr1d, log_r8_arr2d, log_r8_arr3d
  public :: log_int_arr1d

contains

  ! =====================================================================
  ! Internal helpers
  ! =====================================================================

  integer function find_sub_index(sub_name)
    character(len=*), intent(in) :: sub_name
    integer :: i
    find_sub_index = 0
    do i = 1, num_tracked_subs
      if (trim(tracked_sub_names(i)) == trim(sub_name)) then
        find_sub_index = i
        return
      end if
    end do
  end function

  subroutine make_dir(path)
    character(len=*), intent(in) :: path
    integer(c_int) :: ierr
    ierr = c_mkdir(trim(path) // c_null_char, int(o'755', c_int))
    ! ierr /= 0 with EEXIST is benign
  end subroutine

  ! =====================================================================
  ! Trace lifecycle
  ! =====================================================================

  subroutine io_trace_init()
    if (trace_initialized) return
    call make_dir(trace_base_dir)
    trace_initialized = .true.
  end subroutine

  subroutine io_trace_begin(sub_name, call_id)
    character(len=*), intent(in)  :: sub_name
    integer,          intent(out) :: call_id
    character(len=512) :: dir_path
    integer :: idx

    if (.not. trace_initialized) call io_trace_init()

    ! --- per-subroutine call counting ---
    idx = find_sub_index(sub_name)
    if (idx == 0) then
      ! first time seeing this subroutine
      num_tracked_subs = num_tracked_subs + 1
      idx = num_tracked_subs
      tracked_sub_names(idx) = sub_name
      tracked_sub_counts(idx) = 0
    end if
    tracked_sub_counts(idx) = tracked_sub_counts(idx) + 1

    if (tracked_sub_counts(idx) > max_traces_per_sub) then
      call_id = -1   ! sentinel: skip all logging for this call
      return
    end if

    global_call_counter = global_call_counter + 1
    call_id = global_call_counter

    write(dir_path, '(a,"/call_",i6.6,"__",a)') &
      trim(trace_base_dir), call_id, trim(sub_name)
    call make_dir(dir_path)
  end subroutine

  subroutine io_trace_end(call_id)
    integer, intent(in) :: call_id
    ! placeholder for future per-call cleanup
  end subroutine

  subroutine io_trace_open_stage(call_id, sub_name, stage, funit)
    integer,          intent(in)  :: call_id
    character(len=*), intent(in)  :: sub_name, stage
    integer,          intent(out) :: funit
    character(len=512) :: fpath

    if (call_id < 0) then
      funit = -1
      return
    end if

    write(fpath, '(a,"/call_",i6.6,"__",a,"/",a,".txt")') &
      trim(trace_base_dir), call_id, trim(sub_name), trim(stage)
    open(newunit=funit, file=trim(fpath), status='replace', action='write')
  end subroutine

  subroutine io_trace_close_stage(funit)
    integer, intent(in) :: funit
    if (funit < 0) return
    close(funit)
  end subroutine

  ! =====================================================================
  ! Scalar loggers
  ! =====================================================================

  subroutine log_r8(funit, name, val)
    integer,          intent(in) :: funit
    character(len=*), intent(in) :: name
    real(r8),         intent(in) :: val
    if (funit < 0) return
    write(funit, '(a," = ",es23.15e3)') trim(name), val
  end subroutine

  subroutine log_int(funit, name, val)
    integer,          intent(in) :: funit
    character(len=*), intent(in) :: name
    integer,          intent(in) :: val
    if (funit < 0) return
    write(funit, '(a," = ",i0)') trim(name), val
  end subroutine

  subroutine log_char(funit, name, val)
    integer,          intent(in) :: funit
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: val
    if (funit < 0) return
    write(funit, '(a," = ",a)') trim(name), trim(val)
  end subroutine

  ! =====================================================================
  ! Array loggers  (assumed-shape – works with pointers & allocatables)
  ! =====================================================================

  subroutine log_r8_arr1d(funit, name, arr)
    integer,          intent(in) :: funit
    character(len=*), intent(in) :: name
    real(r8),         intent(in) :: arr(:)
    integer :: i, n
    if (funit < 0) return
    n = size(arr)
    write(funit, '(a," [r8_1d, n=",i0,"]")') trim(name), n
    do i = 1, n
      write(funit, '(2x,es23.15e3)') arr(i)
    end do
    write(funit, '(a)') '---'
  end subroutine

  subroutine log_r8_arr2d(funit, name, arr)
    integer,          intent(in) :: funit
    character(len=*), intent(in) :: name
    real(r8),         intent(in) :: arr(:,:)
    integer :: i, j, m, n
    if (funit < 0) return
    m = size(arr, 1)
    n = size(arr, 2)
    write(funit, '(a," [r8_2d, shape=",i0,"x",i0,"]")') trim(name), m, n
    do j = 1, n
      do i = 1, m
        write(funit, '(2x,es23.15e3)', advance='no') arr(i,j)
      end do
      write(funit, *)
    end do
    write(funit, '(a)') '---'
  end subroutine

  subroutine log_r8_arr3d(funit, name, arr)
    integer,          intent(in) :: funit
    character(len=*), intent(in) :: name
    real(r8),         intent(in) :: arr(:,:,:)
    integer :: i, j, l, m, n, k
    if (funit < 0) return
    m = size(arr, 1)
    n = size(arr, 2)
    k = size(arr, 3)
    write(funit, '(a," [r8_3d, shape=",i0,"x",i0,"x",i0,"]")') trim(name), m, n, k
    do l = 1, k
      write(funit, '(a,i0)') '  k=', l
      do j = 1, n
        do i = 1, m
          write(funit, '(2x,es23.15e3)', advance='no') arr(i,j,l)
        end do
        write(funit, *)
      end do
    end do
    write(funit, '(a)') '---'
  end subroutine

  subroutine log_int_arr1d(funit, name, arr)
    integer,          intent(in) :: funit
    character(len=*), intent(in) :: name
    integer,          intent(in) :: arr(:)
    integer :: i, n
    if (funit < 0) return
    n = size(arr)
    write(funit, '(a," [int_1d, n=",i0,"]")') trim(name), n
    do i = 1, n
      write(funit, '(2x,i0)') arr(i)
    end do
    write(funit, '(a)') '---'
  end subroutine

end module io_logger
