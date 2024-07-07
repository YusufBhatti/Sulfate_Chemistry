! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE time_utils_mod

USE missing_data_mod, ONLY: imdi

! DrHook Modules
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

! Description:
!   Time format utilities. Convert a fieldsfile-header style timedate to an
!   integer number of seconds since an epoch. Wraps a UM subroutine. Note that
!   the calendar must be set using set_calendar before you can convert.
!   This  supports both Gregorian and 360-day calendars, with the setting being
!   read from the fixed header of the input file (fixed_header(8)).
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.

INTEGER :: calendar = imdi

! These values correspond to the header values for fixed_header(8) "Indicator
! for calendar" according to UMDP F3.
INTEGER, PARAMETER :: cal_Gregorian  = 1
INTEGER, PARAMETER :: cal_360day     = 2

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'TIME_UTILS_MOD'

CONTAINS
  
!-------------------------------------------------------------------------------
INTEGER FUNCTION timedate_to_seconds(input_timedate)
! Convert a timedate array into integer seconds

USE conversions_mod, ONLY: isec_per_day
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE

INTEGER, INTENT(IN) :: input_timedate(6)

INTEGER, PARAMETER :: basis_time_days = 0
INTEGER, PARAMETER :: basis_time_secs = 0

LOGICAL :: lcal360

INTEGER :: elapsed_time_secs, elapsed_time_days

INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER :: RoutineName='TIMEDATE_TO_SECONDS'

! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

SELECT CASE (calendar)
  CASE(cal_Gregorian)
    lcal360 = .FALSE.
  CASE(cal_360day)
    lcal360 = .TRUE.
  CASE DEFAULT
    cmessage = 'Calendar not defined'
    icode = 1
    CALL ereport(RoutineName, icode, cmessage)
END SELECT

  
! DEPENDS ON: time2sec    
CALL time2sec(input_timedate(1), input_timedate(2), input_timedate(3),   &
              input_timedate(4), input_timedate(5), input_timedate(6),   &
              basis_time_days, basis_time_secs, elapsed_time_days,       &
              elapsed_time_secs, lcal360) 

timedate_to_seconds = elapsed_time_days * isec_per_day + elapsed_time_secs

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION timedate_to_seconds 

!-------------------------------------------------------------------------------

SUBROUTINE set_calendar(value)
! Set the calendar type

USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE

INTEGER, INTENT(IN) :: value

INTEGER :: icode
CHARACTER(LEN=*), PARAMETER :: routinename='SET_CALENDAR'
CHARACTER(LEN=errormessagelength) :: cmessage

! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

SELECT CASE (value)
  CASE(cal_Gregorian, cal_360day)
    CONTINUE
  CASE DEFAULT
    WRITE(cmessage,'(A,I0)') 'Invalid calendar defined: ',value
    icode = 1
    CALL ereport(routinename, icode, cmessage)  
END SELECT

IF (calendar /= imdi .AND. calendar /= value) THEN
  ! We're changing the calendar. This is probably unwise so abort.
  WRITE(cmessage,'(A,I0,A,I0)') 'Attempted to change calendar from ', calendar,&
                                ' to ', value
  icode = 1
  CALL ereport(routinename, icode, cmessage)
END IF
calendar = value

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_calendar
!-------------------------------------------------------------------------------

END MODULE time_utils_mod

