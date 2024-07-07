! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE totimp_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='TOTIMP_MOD'

CONTAINS

! Convert multiple of specified time period to no. of timesteps
! Function Interface:

INTEGER FUNCTION totimp(period,time_unit,mdl)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE nlstgen_mod, ONLY: steps_per_periodim, secs_per_periodim, dumpfreqim
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE cstash_mod, ONLY:  &
    stsh_timesteps, stsh_hours, stsh_days, stsh_dumps, stsh_minutes, &
    stsh_seconds

IMPLICIT NONE

! Description: Get no. of timesteps from time period information taken
!  from STASH profiles to convert to no. of timesteps required to
!  generate STASH list contents for controlling diagnostic output times.
!
! Method: Simple conversion of specified time periods. Illegal
!  combinations are returned as -999 to be trapped by calling routine.
!
! Note: Negative periods are allowed since the namelist allows time profiles
! with a start and end date indicating when the diagnostic is to be
! written out: Model runs stopped by the operator and restarted may result
! in a start date for the diagnostic which is earlier than the model basis
! time.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Stash
!
!  Code description:
!    FORTRAN 90
!    Written to UMDP 003 programming standards.
!

! Function arguments:
!   Scalar arguments with intent(in):
INTEGER, INTENT(IN)        :: period     ! multiples of time period
INTEGER, INTENT(IN)        :: time_unit  ! descriptor for time period
INTEGER, INTENT(IN)        :: mdl        ! internal model

! Local scalars:
INTEGER, PARAMETER::      totimp_error=-999 ! error marker
REAL :: fac

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='TOTIMP'

!- End of Header ----------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

fac= REAL(steps_per_periodim(mdl))/REAL(secs_per_periodim(mdl))

IF (time_unit == stsh_timesteps) THEN           ! timesteps
  totimp=period
ELSE IF (time_unit == stsh_hours) THEN          ! hours
  ! Convert (negative or positive) hours into timesteps using nearest integer.
  totimp=NINT(period*fac*3600.0)

ELSE IF (time_unit == stsh_days) THEN           ! days
  ! Convert (negative or positive) days into timesteps using nearest integer.
  totimp=NINT(period*fac*3600.0*24.0)

ELSE IF (time_unit == stsh_dumps) THEN          ! dump periods
  IF (DUMPFREQim(mdl) == 0) THEN
    WRITE(umMessage,*)'TOTIMP:IRREGULAR DUMPS FOR DUMP FREQUENCY'
    CALL umPrint(umMessage,src='totimp')
    totimp= totimp_error
  ELSE
    totimp=period*dumpfreqim(mdl)
  END IF

ELSE IF (time_unit == stsh_minutes) THEN        ! minutes
  ! Convert (negative or positive) minutes to timesteps using nearest integer.
  totimp=NINT(period*fac*60.0)

ELSE IF (time_unit == stsh_seconds) THEN        ! seconds
  ! Convert (negative or positive) seconds to timesteps using nearest integer.
  totimp=NINT(period*fac)

ELSE                                            ! illegal unit
  totimp= totimp_error
  WRITE(umMessage,*)'TOTIMP: UNEXPECTED TIME UNIT=',time_unit
  CALL umPrint(umMessage,src='totimp')
END IF
! Note: the special case of period  ==  -1 (indefinite) is trapped by
! the calling routine, otherwise it would be necessary to include lines
! ELSE IF(PERIOD  ==  -1) THEN TOTIMP= -1 here.


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END FUNCTION totimp

END MODULE totimp_mod
