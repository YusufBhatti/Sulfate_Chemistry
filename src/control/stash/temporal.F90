! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    Routine: TEMPORAL -------------------------------------------------
!
!    Purpose: Control routine to handle temporal processing options
!             within STASH.  Its input and output arguments look like
!             1D arrays (ie. all the data should be in contiguous areas
!             of memory).  Lower level service routines are called to
!             perform the individual processing options.
!
!    Programming standard: UM Doc Paper 3
!
!    External documentation:
!      Unified Model Doc Paper C4 - Storage handling and diagnostic
!                                   system (STASH)
!
!    Interface and arguments: ------------------------------------------
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: STASH

SUBROUTINE temporal(variable,RESULT,fsize,extra_size,             &
  control,control_size,                                            &
  timestep,icode,cmessage,start,amdi)
!

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE stparam_mod, ONLY: st_proc_no_code, st_replace_code, st_accum_code,&
    st_start_time_code, st_time_mean_code, st_period_code,             &
    st_infinite_time, st_freq_code, st_gridpoint_code, st_max_code,    &
    st_min_code, st_time_series_code, st_append_traj_code,             &
    st_time_series_mean, block_size
USE sterr_mod, ONLY: st_not_supported, unknown_processing
USE ereport_mod, ONLY: ereport
USE umPrintMgr, ONLY: umPrint, umMessage
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE
!
INTEGER :: fsize                 ! IN  size of arrays
REAL :: variable(fsize)          ! IN  data array
REAL :: RESULT(fsize)            ! OUT output array
INTEGER :: extra_size            ! IN  size of extra data
INTEGER :: control_size          ! IN  size of control
INTEGER :: control(control_size) ! IN  control
INTEGER :: timestep              ! IN  present value of timestep
INTEGER :: icode                 ! OUT error code
CHARACTER(LEN=errormessagelength) :: cmessage     ! OUT error message
REAL :: amdi                     ! IN  missing data indicator
LOGICAL :: start                 ! OUT true if start timestep

! ----------------------------------------------------------------------
!
! Local variables
!
INTEGER :: masking        ! control for masking/treating missing data
INTEGER :: proc_code      ! value of processing code
INTEGER :: mask_code      ! value of masking code
REAL :: divisor           ! divisor for the time mean (1/period)
INTEGER :: mod_period     ! timesteps since start modulo period.
INTEGER :: start_time     ! value of start time
INTEGER :: i              ! loop counter
LOGICAL :: last_ts        ! true if end timestep
INTEGER :: proc_size      ! size of data to be processed

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='TEMPORAL'
! ---------------------------------------------------------------------
!  1. Set processing option code and select appropriate service routine
!
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
proc_size=fsize-extra_size
proc_code=control(st_proc_no_code)
!
!  Replace (null processing)
!
IF (proc_code == st_replace_code) THEN
  DO i=1,fsize
    RESULT(i)=variable(i)
  END DO
  start=(control(st_start_time_code) == timestep)
  !
  !  Mean/accumulation
  !
ELSE IF (proc_code == st_accum_code .OR.                           &
        proc_code == st_time_mean_code) THEN
  start_time=control(st_start_time_code)
  IF (control(st_period_code) == st_infinite_time) THEN
    start=(timestep == start_time)
    last_ts=.FALSE.
  ELSE
    mod_period=MOD(timestep-start_time,control(st_period_code))
    start=(mod_period == 0)
    last_ts=(mod_period == (control(st_period_code)-              &
                        control(st_freq_code)))
  END IF

  mask_code=control(st_gridpoint_code)
  mask_code=MOD(mask_code,block_size)
  masking=mask_code

  IF (start) THEN      ! first timestep.
    DO i=1,fsize
      RESULT(i)=variable(i)
    END DO
  ELSE
    ! DEPENDS ON: staccum
    CALL staccum(variable,RESULT,proc_size,masking)
    DO i=proc_size+1,fsize
      RESULT(i)=variable(i) ! copy over the extra data (if any)
    END DO
  END IF
  !  Normalise at end of mean period
  IF (last_ts .AND. proc_code == st_time_mean_code) THEN
    divisor=(REAL(control(st_freq_code))/                         &
             REAL(control(st_period_code)))
    ! If field is masked test for MDI, otherwise don't
    IF (masking > 1) THEN
      DO i=1,proc_size
        IF (RESULT(i) /= amdi) THEN
          RESULT(i)=RESULT(i)*divisor
        END IF
      END DO
    ELSE
      DO i=1,proc_size
        RESULT(i)=RESULT(i)*divisor
      END DO
    END IF
  END IF
  !
  !  Maximum
  !
ELSE IF (proc_code == st_max_code) THEN
  start_time=control(st_start_time_code)
  mod_period=MOD(timestep-start_time,control(st_period_code))
  start=(mod_period == 0)
  IF (start) THEN
    DO i=1,fsize
      RESULT(i)=variable(i)
    END DO
  ELSE
    mask_code=control(st_gridpoint_code)
    mask_code=MOD(mask_code,block_size)
    masking=mask_code
    ! DEPENDS ON: stmax
    CALL stmax(variable,RESULT,proc_size,masking)
    DO i=proc_size+1,fsize
      RESULT(i)=variable(i) ! copy over the extra data (if any)
    END DO
  END IF
  !
  !  Minimum
  !
ELSE IF (proc_code == st_min_code) THEN
  start_time=control(st_start_time_code)
  mod_period=MOD(timestep-start_time,control(st_period_code))
  start=(mod_period == 0)
  IF (start) THEN
    DO i=1,fsize
      RESULT(i)=variable(i)
    END DO
  ELSE
    mask_code=control(st_gridpoint_code)
    mask_code=MOD(mask_code,block_size)
    masking=mask_code
    ! DEPENDS ON: stmin
    CALL stmin(variable,RESULT,proc_size,masking)
    DO i=proc_size+1,fsize
      RESULT(i)=variable(i) ! copy over the extra data (if any)
    END DO
  END IF
  !
  !  Timeseries (append)
  !
ELSE IF (proc_code == st_time_series_code) THEN
  DO i=1,fsize
    ! Note that on start timestep this will include the extra data
    RESULT(i)=variable(i)
  END DO
  start_time=control(st_start_time_code)
  mod_period=MOD(timestep-start_time,control(st_period_code))
  start=(mod_period == 0)
  last_ts=(mod_period == (control(st_period_code)-                &
                      control(st_freq_code)))
  !
  !  Append trajectories
  !
ELSE IF (proc_code == st_append_traj_code) THEN
  start_time=control(st_start_time_code)
  mod_period=MOD(timestep-start_time,control(st_period_code))
  start=(mod_period == 0)
  last_ts=(mod_period == (control(st_period_code)-                &
                      control(st_freq_code)))
  icode=st_not_supported
  cmessage='TEMPORAL: does not support append trajectories'
  GO TO 9999
  !
  !  Timeseries (append) - option 8 daily mean
  !
ELSE IF (proc_code == st_time_series_mean) THEN

  DO i=1,fsize
    ! Note that on start timestep this will include the extra data
    RESULT(i)=variable(i)
  END DO
  start_time=control(st_start_time_code)
  mod_period=MOD(timestep-start_time,control(st_period_code))
  start=(mod_period == 0)
  last_ts=(mod_period == (control(st_period_code)-                &
                      control(st_freq_code)))
  !
  !  Error condition
  !
ELSE
  icode=unknown_processing
  WRITE(cmessage,'(A,1x,i5)')'TEMPORAL: unknown processing code', &
                              proc_code
  GO TO 9999
END IF
!
9999  CONTINUE   ! jump for errors
!
IF (icode >  0) THEN
  WRITE(umMessage,'(A,A,I5)')                                     &
  'TEMPORAL: Error in temporal processing:',', code ',icode
  CALL umPrint(umMessage,src='temporal')
  CALL ereport('TEMPORAL',icode,cmessage)
END IF
!
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE temporal
