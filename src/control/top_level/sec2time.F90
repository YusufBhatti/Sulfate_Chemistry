! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    Routine: SEC2TIM0 -------------------------------------------------
!
!    Purpose: Converts from an integer number of elapsed seconds since
!             the model basis time to a calendar date/time, using the
!             absolute calendar zero point as a reference.  30-day
!             month or standard calendar may be used.
!             NB: BASIS_TIME_SECS is the number of seconds from the
!             calendar zero point to the basis time for the run, and
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

SUBROUTINE sec2time(elapsed_days,elapsed_secs                     &
,                   basis_time_days,basis_time_secs               &
,                   i_year,i_month,i_day,i_hour,i_minute,i_second &
,                   i_day_number,lcal360)
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
     elapsed_days,                                                &
                             ! IN  - elapsed days since basis time
     elapsed_secs,                                                &
                             ! IN  - elapsed secs in part of day
     basis_time_days,                                             &
                             ! IN  - whole days to basis time
     basis_time_secs,                                             &
                             ! IN  - secs in day at basis time
!                                  !       relative to calendar zero
           i_second,                                                    &
                                   ! OUT - model time (seconds)
           i_minute,                                                    &
                                   ! OUT - model time (minutes)
           i_hour,                                                      &
                                   ! OUT - model time (hours)
           i_day,                                                       &
                                   ! OUT - model time (days)
           i_month,                                                     &
                                   ! OUT - model time (months)
           i_year,                                                      &
                                   ! OUT - model time (years)
           i_day_number            ! OUT - model time (day number)

!
! ----------------------------------------------------------------------
!  Local variables
!
INTEGER ::                                                        &
       n_second      ! number of seconds since calendar zero

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SEC2TIME'
!
! ----------------------------------------------------------------------
!  1. Add elapsed time to basis time in days/seconds to get elapsed
!      since calendar zero, and convert to hours, minutes, seconds and
!      total days since calendar zero
!
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
n_second   = basis_time_secs+elapsed_secs
i_day      = basis_time_days+elapsed_days+n_second/isec_per_day
IF (n_second  >=  0) THEN
  n_second = MOD(n_second,isec_per_day)
ELSE
  n_second = MOD(n_second,isec_per_day)
  ! Correction for fact that rounding was upwards.
  IF (n_second  /=  0) THEN
    i_day  = i_day - 1
    n_second = n_second + isec_per_day
  END IF
END IF
i_hour   = MOD(n_second/3600,24)
i_minute = MOD(n_second/60  ,60)
i_second = MOD(n_second,60)
! ----------------------------------------------------------------------
!  2. Convert day number to date
!
IF (lcal360) THEN
  !
  !  2.1 30-day month (360 day year) calendar
  !
  i_year  = i_day/360
  i_month = MOD(i_day/30,12)+1
  i_day   = MOD(i_day,30)+1
  i_day_number = i_day+30*(i_month-1)

ELSE
  !
  !  2.2 Gregorian calendar
  !
  i_year = (i_day/days_per_4c)*400
  i_day = i_day-(i_day/days_per_4c)*days_per_4c
  !      Catch special case 31 Dec in leap years
  IF (i_day == 4*days_per_c) THEN
    i_year = i_year+400
    i_day = days_per_y+1
  ELSE
    i_year = i_year+(i_day/days_per_c)*100
    i_day = i_day-(i_day/days_per_c)*days_per_c
    i_year = i_year+(i_day/days_per_4y)*4
    i_day = i_day-(i_day/days_per_4y)*days_per_4y
    IF (i_day == 4*days_per_y) THEN
      i_year = i_year+4
      i_day = days_per_y+1
    ELSE
      i_year = i_year+(i_day/days_per_y) + 1
      i_day = i_day-(i_day/days_per_y)*days_per_y + 1
    END IF
  END IF
  i_day_number = i_day
  !      Find month/day from day no in year
  i_month = 1
  DO WHILE ((i_month  <=  12) .AND.                                 &
    (i_day  >   days_in_month(i_month)))
    i_day = i_day-days_in_month(i_month)
    i_month = i_month+1
  END DO
  !      Adjust if leap year and after February
  IF (i_month >  2 .AND. MOD(i_year,4) == 0 .AND.                   &
      (MOD(i_year,400) == 0 .OR. MOD(i_year,100) /= 0)) THEN
    i_day = i_day-1
    IF (i_day == 0) THEN
      i_month = i_month-1
      i_day = days_in_month(i_month)
      IF (i_month == 2) i_day=29
    END IF
  END IF
END IF  !  LCAL360
! ----------------------------------------------------------------------
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE sec2time
