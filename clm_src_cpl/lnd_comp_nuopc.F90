module lnd_comp_nuopc

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Interface of the active land model component of CESM (CLM, Community Land Model)
  ! with the main CESM driver.
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use decompMod, only : bounds_type
#ifdef IO_TRACE
  use io_logger
#endif
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: InitializeRealize  ! CLM initialization
  public :: ModelAdvance       ! CLM run phase
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine InitializeRealize (bounds)
    !
    ! !DESCRIPTION:
    ! Initialize land surface model
    !
    ! !USES:
    use clm_initializeMod, only : initialize1, initialize2
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds
    !---------------------------------------------------------------------

    call initialize1 ()
    call initialize2 (bounds)

  end subroutine InitializeRealize

  !-----------------------------------------------------------------------
  subroutine ModelAdvance (bounds, time_indx, fin1, fin2)
    !
    ! !DESCRIPTION:
    ! Run CLM model
    !
    ! !USES:
    use clm_driver, only : clm_drv
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: time_indx            ! Time index from reference date (0Z January 1 of current year, when calday = 1.000
    character(len=256) :: fin1, fin2            ! File name
    !
    ! !LOCAL VARIABLES:
#ifdef IO_TRACE
    integer :: io_call_id, io_unit
#endif
    !---------------------------------------------------------------------

#ifdef IO_TRACE
    call io_trace_begin("ModelAdvance", io_call_id)
    call io_trace_open_stage(io_call_id, "ModelAdvance", "inputs", io_unit)
    call log_int(io_unit, "bounds%begp", bounds%begp)
    call log_int(io_unit, "bounds%endp", bounds%endp)
    call log_int(io_unit, "time_indx", time_indx)
    call log_char(io_unit, "fin1", fin1)
    call log_char(io_unit, "fin2", fin2)
    call io_trace_close_stage(io_unit)
#endif

    call clm_drv (bounds, time_indx, fin1, fin2)

#ifdef IO_TRACE
    call io_trace_open_stage(io_call_id, "ModelAdvance", "outputs", io_unit)
    call log_char(io_unit, "status", "completed")
    call io_trace_close_stage(io_unit)
    call io_trace_end(io_call_id)
#endif

  end subroutine ModelAdvance

end module lnd_comp_nuopc
