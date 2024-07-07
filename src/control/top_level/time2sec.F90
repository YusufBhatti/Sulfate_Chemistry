! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    Routine: TIME2SEC -------------------------------------------------
!
!    Purpose: Converts from calendar date/time to an integer number
!             of elapsed seconds since the model basis time, using the
!             absolute calendar zero point as a reference.  30-day
!             month or standard calendar may be used.
!             NB: BASIS_TIME_SECS is the number of seconds from the
!             calendar zero point to the basis time for the model, and
!             is calculated in INITTIME.
!
!    Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!
!    Logical components covered: S62
!
!    Project task: S62
!
!    External documentation: On-line UM document C0 - The top-level
!                            control system
!
!    -------------------------------------------------------------------
!    Interface and arguments: ------------------------------------------
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Top Level

SUBROUTINE time2sec (i_year,i_month,i_day,i_hour,i_minute,i_second&
,                    basis_time_days,basis_time_secs              &
,                    elapsed_days,elapsed_secs,lcal360)
!

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE cdaydata_mod, ONLY: days_per_4c, days_per_c, days_per_4y,     &
                        days_per_y, days_in_month, days_to_month
USE conversions_mod, ONLY: isec_per_day

IMPLICIT NONE

LOGICAL :: lcal360
!
INTEGER ::                                                        &
     i_year,                                                      &
                             ! IN  - model time (years)
     i_month,                                                     &
                             ! IN  - model time (months)
     i_day,                                                       &
                             ! IN  - model time (days)
     i_hour,                                                      &
                             ! IN  - model time (hours)
     i_minute,                                                    &
                             ! IN  - model time (minutes)
     i_second,                                                    &
                             ! IN  - model time (seconds)
     basis_time_days,                                             &
                             ! IN  - whole days to basis time
     basis_time_secs,                                             &
                             ! IN  - secs in day at basis time
     elapsed_days,                                                &
                             ! OUT - elapsed days since basis time
     elapsed_secs            ! OUT - elapsed secs in part of day
!                                  !       relative to basis time
!
! ----------------------------------------------------------------------
!  Local variables
!
INTEGER ::                                                        &
       year                                                       &
                   ! years
,      day         ! number of days since calendar zero

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='TIME2SEC'
!
! ----------------------------------------------------------------------
!  1. Add up days from time zero to specified time
!
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
IF (lcal360) THEN
  !
  !  1.1 30-day month (360 day year) calendar
  !
  day = 360*i_year+30*(i_month-1)+i_day-1
  !
ELSE
  !
  !  1.2 Gregorian calendar
  !
  !      If leap year and after 28 February, adjust day number by one
  IF (MOD(i_year,4) == 0   .AND.                                    &
     (MOD(i_year,400) == 0 .OR. MOD(i_year,100) /= 0) .AND.         &
      i_month >  2) THEN
    day = i_day+1
  ELSE
    day = i_day
  END IF
  !      Add on days in the preceding months in the year
  day = day + days_to_month(i_month) - 1
  year = i_year - 1
  !      Add on days up to the specified year
  day =  day+(year/400)*days_per_4c
  year = year-(year/400)*400
  day =  day+(year/100)*days_per_c
  year = year-(year/100)*100
  day =  day+(year/4)*days_per_4y
  year = year-(year/4)*4
  day =  day+year*days_per_y
  !
END IF       ! LCAL360
! ----------------------------------------------------------------------
!  2. Convert days, hours and minutes to days/secs since calendar zero,
!      and subtract basis time in days/secs to get elapsed time since
!      basis, converted to whole days and +ve no of secs in partial day
!
elapsed_days=day-basis_time_days
elapsed_secs=3600*i_hour+60*i_minute+i_second-basis_time_secs
IF (elapsed_secs >= isec_per_day) THEN
  elapsed_days=elapsed_days+elapsed_secs/isec_per_day
  elapsed_secs=MOD(elapsed_secs,isec_per_day)
ELSE IF (elapsed_secs <  0) THEN
  elapsed_days=elapsed_days+(elapsed_secs+1-isec_per_day)/isec_per_day
  elapsed_secs=MOD(elapsed_secs,isec_per_day)+isec_per_day
END IF
!
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
! ----------------------------------------------------------------------
END SUBROUTINE time2sec
