! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    Routine: INCRTIME -------------------------------------------------
!
!    Purpose: Increments the model time by one atmosphere timestep.
!             Also updates timestamps in dump LOOKUP headers of
!             PROGNOSTIC fields (diagnostic LOOKUP headers
!             are updated exclusively by STASH).
!
!    Programming standard: UM Doc Paper 3
!
!    External documentation: On-line UM document C0 - The top-level
!                            control system
!
!    -------------------------------------------------------------------
!    Interface and arguments: ------------------------------------------
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Top Level

SUBROUTINE incrtime(                                              &
      internal_model,icode,cmessage )
!

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE IAU_mod, ONLY: &
    l_iau

USE umPrintMgr
USE UM_ParParams
USE Control_Max_Sizes
USE lookup_addresses
USE dump_headers_mod, ONLY: a_fixhd, a_lookup
USE dynamics_testing_mod, ONLY:  L_Backwards
USE nlstgen_mod, ONLY: steps_per_periodim, secs_per_periodim
USE submodel_mod, ONLY: atmos_sm
USE nlstcall_mod, ONLY: model_analysis_hrs, &
                         model_analysis_mins, &
                         l_fastrun, lcal360

USE history, ONLY: h_stepim

USE nlsizes_namelist_mod, ONLY:                                        &
    a_len1_coldepc, a_len1_flddepc, a_len1_levdepc, a_len1_rowdepc,    &
    a_len2_coldepc, a_len2_flddepc, a_len2_levdepc, a_len2_lookup,     &
    a_len2_rowdepc, a_len_cfi1, a_len_cfi2, a_len_cfi3, a_len_extcnst, &
    a_len_inthd, a_len_realhd, a_prog_lookup, len1_lookup,             &
    len_dumphist, len_fixhd, mpp_len1_lookup

USE model_time_mod, ONLY:                                                 &
    basis_time_days, basis_time_secs, data_minus_basis_hrs, forecast_hrs, &
    i_day, i_day_number, i_hour, i_minute, i_month, i_second, i_year,     &
    previous_time, stepim
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!
!   Arguments
!
INTEGER ::    internal_model    ! IN : internal model identifier
INTEGER ::    icode             ! OUT: Error return code
CHARACTER(LEN=errormessagelength) :: cmessage  ! OUT: Error return message
!
! ----------------------------------------------------------------------
!  Local variables
!
INTEGER ::                                                        &
    elapsed_days,                                                &
                             ! Elapsed days  since basis time
    elapsed_secs,                                                &
                             ! Elapsed secs  since basis time
    elapsed_days_prev,                                           &
                             ! Elapsed days, end of previous step
    elapsed_secs_prev,                                           &
                             ! Elapsed secs, end of previous step
    i                       ! Loop index
INTEGER ::                                                        &
                             ! Local scalars of internal model
 step                                                            &
                             !  arrays.
, steps_per_period                                                &
, secs_per_period

! model_analysis_hrs replaced by model_analysis_mins in cntlall.h -
! Requires ELAPSED_HRS changed to REAL
REAL :: elapsed_hrs             ! Elapsed hours since basis time

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INCRTIME'

!
! ----------------------------------------------------------------------
!  1. General timestep, increment STEP by one and update atmos
!     elapsed seconds (integer) relative to basis time
!
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
step            =stepim(internal_model)
steps_per_period=steps_per_periodim(internal_model)
secs_per_period =secs_per_periodim(internal_model)

! DEPENDS ON: stp2time
CALL stp2time(step,steps_per_period,secs_per_period,            &
             elapsed_days_prev,elapsed_secs_prev)
step = step+1
! DEPENDS ON: stp2time
CALL stp2time(step,steps_per_period,secs_per_period,            &
             elapsed_days,elapsed_secs)

stepim(internal_model) = step
h_stepim(internal_model) = step

!
!  1.1 If integrating backwards, negate elapsed times.
!
IF (L_Backwards .AND. internal_model  ==  atmos_sm) THEN
  elapsed_days_prev = -elapsed_days_prev
  elapsed_secs_prev = -elapsed_secs_prev
  elapsed_days      = -elapsed_days
  elapsed_secs      = -elapsed_secs
END IF

!
!  1.2 Set FORECAST_HRS - the number of hours relative to the current
!                         data time.
!
model_analysis_hrs =  REAL(model_analysis_mins)/60.0
elapsed_hrs = REAL(elapsed_secs)/3600 + elapsed_days*24

IF ((l_fastrun .OR. l_iau) .AND.             &
    elapsed_hrs >= model_analysis_hrs) THEN
  ! At or beyond data time reset step:
  forecast_hrs = elapsed_hrs - model_analysis_hrs
ELSE
  ! Data time of initial dump still valid:
  forecast_hrs = elapsed_hrs - data_minus_basis_hrs
END IF

IF ( PrintStatus  >  PrStatus_Oper ) THEN
  CALL umPrint( '',src='incrtime' )
  WRITE(umMessage,'(A,I10)')   'incrtime: ELAPSED SECS       ', &
      elapsed_secs
  CALL umPrint(umMessage,src='incrtime')
  WRITE(umMessage,'(A,F10.5)') 'incrtime: ELAPSED HRS        ', &
      elapsed_hrs
  CALL umPrint(umMessage,src='incrtime')
  WRITE(umMessage,'(A,F10.5)') 'incrtime: FORECAST HRS       ', &
      forecast_hrs
  CALL umPrint(umMessage,src='incrtime')
  WRITE(umMessage,'(A,F10.5)') 'incrtime: model_analysis_hrs ', &
      model_analysis_hrs
  CALL umPrint(umMessage,src='incrtime')
  CALL umPrint('',src='incrtime')
END IF

! ----------------------------------------------------------------------
!  2. Convert elapsed seconds since basis time to calendar time/date
!
! DEPENDS ON: sec2time
CALL sec2time(elapsed_days_prev,elapsed_secs_prev,                &
             basis_time_days,basis_time_secs,                    &
             previous_time(1),previous_time(2),previous_time(3), &
             previous_time(4),previous_time(5),previous_time(6), &
             previous_time(7),lcal360)
! DEPENDS ON: sec2time
CALL sec2time(elapsed_days,elapsed_secs,                          &
             basis_time_days,basis_time_secs,                    &
             i_year,i_month,i_day,i_hour,i_minute,i_second,      &
             i_day_number,lcal360)
! ----------------------------------------------------------------------
!  3. Copy date/time information into the dump header and update
!     VALIDITY TIME and FORECAST PERIOD in prognostic field LOOKUP
!     headers.
!
a_fixhd(28) = i_year
a_fixhd(29) = i_month
a_fixhd(30) = i_day
a_fixhd(31) = i_hour
a_fixhd(32) = i_minute
a_fixhd(33) = i_second
a_fixhd(34) = i_day_number
DO i=1,a_prog_lookup
  a_lookup(lbyr  ,i)=i_year
  a_lookup(lbmon ,i)=i_month
  a_lookup(lbdat ,i)=i_day
  a_lookup(lbhr  ,i)=i_hour
  a_lookup(lbmin ,i)=i_minute
  a_lookup(lbsec ,i)=i_second
  a_lookup(lbft,i)=forecast_hrs
END DO

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)

RETURN
! ----------------------------------------------------------------------
END SUBROUTINE incrtime
