! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
MODULE check_nlstcgen_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CHECK_NLSTCGEN_MOD'

CONTAINS

SUBROUTINE check_nlstcgen()

! Description:
!   Subroutine to apply logic controls and set control variables based on the
!   options selected in the nlstcgen namelist.
!
!   Converts the dumping frequency to timesteps if it was specified in
!   hours or days.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

! Code Description:
!   Language: FORTRAN 90


! Dr Hook Modules
USE yomhook,           ONLY: lhook, dr_hook
USE parkind1,          ONLY: jprb, jpim
USE nlstcall_subs_mod, ONLY: totime
USE submodel_mod,      ONLY: atmos_im
USE nlstgen_mod,       ONLY: dumpfreqim, i_dump_output, dumptimesim,          &
                             dump_frequency_units, dumptimes_len1, ppxm,      &
                             greg_dump_freq, secs_per_periodim,               &
                             steps_per_periodim, dump_packim, mean_reftimeim, &
                             meanfreqim, ppselectim
USE nlstcall_mod,      ONLY: lclimrealyr
USE conversions_mod,   ONLY: isec_per_day
USE chk_opts_mod,      ONLY: chk_var, def_src

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CHECK_NLSTCGEN'
INTEGER                       :: outtime = 0
INTEGER                       :: i
INTEGER                       :: days_per_period
INTEGER                       :: steps_per_day

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
def_src = RoutineName

! Where possible, check that sensible namelist values have been chosen
CALL chk_var( i_dump_output,        'i_dump_output',        '[1,2,3]' )
IF (i_dump_output > 1) THEN
  CALL chk_var( dump_frequency_units, 'dump_frequency_units', '[1,2,3]' )
  CALL chk_var( dump_packim(1),       'dump_packim(1)',       '[1,2,3]' )
  ! For the Climate Meaning system to work, the user must set year, month, day 
  ! of the reference time, but hour, min & sec must be zero.
  CALL chk_var( mean_reftimeim(4,1),  'mean_reftimeim(hr)',   '[0]' )
  CALL chk_var( mean_reftimeim(5,1),  'mean_reftimeim(min)',  '[0]' )
  CALL chk_var( mean_reftimeim(6,1),  'mean_reftimeim(sec)',  '[0]' )
  CALL chk_var( meanfreqim(1,1),  'meanfreqim(1)',  '[>=0]' )
  CALL chk_var( meanfreqim(2,1),  'meanfreqim(2)',  '[>=0]' )
  CALL chk_var( meanfreqim(3,1),  'meanfreqim(3)',  '[>=0]' )
  CALL chk_var( meanfreqim(4,1),  'meanfreqim(4)',  '[>=0]' )
  CALL chk_var( ppselectim(1,1),  'ppselectim(1)',  '[0:1]' )
  CALL chk_var( ppselectim(2,1),  'ppselectim(2)',  '[0:1]' )
  CALL chk_var( ppselectim(3,1),  'ppselectim(3)',  '[0:1]' )
  CALL chk_var( ppselectim(4,1),  'ppselectim(4)',  '[0:1]' )
  ! Only check ppxm if any of ppselectim == 1
  IF (SUM(ppselectim(:, 1)) > 0) THEN
    CALL chk_var( ppxm,           'ppxm',           '[0,1,2,4,5]' )
  END IF
END IF
CALL chk_var( secs_per_periodim(1), 'secs_per_periodim',  '[86400:1036800]' )
CALL chk_var( steps_per_periodim(1),'steps_per_periodim',   '[1:86400]' )


! No need to loop over second index of arrays as only dealing with atmos model
IF (i_dump_output == 2) THEN ! Regular freq. dumps
  CALL totime(dumpfreqim(atmos_im), dump_frequency_units, outtime)
  dumpfreqim(atmos_im) = outtime
ELSE IF (i_dump_output == 3) THEN ! Irregular freq. dumps
  DO i = 1, dumptimes_len1
    CALL totime(dumptimesim(i, atmos_im), dump_frequency_units, outtime)
    dumptimesim(i, atmos_im) = outtime
  END DO
END IF

! Set up variables for dumping using the gregorian calendar
! Historically Gregorian calendar runs have required a dumping period of
! daily or more frequently. This is set by dumpfreqim, which controls both
! the dumping and accumulation of partial sums for climate meaning.
! In the case where lclimrealyr is set to true, the variable greg_dump_freq
! controls the production of restart dumps, however dumpfreqim controls the
! climate meaning system. To ensure that the climate meaning system is
! unaffected by the flexible dumping system, once greg_dump_freq has been
! set, dumpfreqim is overwritten to have a value corresponding to one day.
days_per_period = secs_per_periodim(1) / isec_per_day
steps_per_day = steps_per_periodim(1) / days_per_period

IF (lclimrealyr) THEN
  IF (dumpfreqim(atmos_im) <= steps_per_day) THEN
    greg_dump_freq = dumpfreqim(atmos_im)
  ELSE
    greg_dump_freq = dumpfreqim(atmos_im)
    dumpfreqim(atmos_im) = steps_per_day
  END IF
END IF

def_src = ''
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_nlstcgen

END MODULE check_nlstcgen_mod

