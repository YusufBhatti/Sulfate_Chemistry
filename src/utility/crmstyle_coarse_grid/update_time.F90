! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Update date time by a given increment

MODULE update_time_mod

IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UPDATE_TIME_MOD'

CONTAINS


SUBROUTINE update_time(in_date, inc_date, out_date)


USE missing_data_mod, ONLY: rmdi

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Description:
! Update date and time by a given increment
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utilty - crmstyle_coarse_grid
!
! Code Description:
!   Language:           Fortran 90
!   This code is written to UM programming standards version 8.3.


INTEGER, INTENT(IN)  ::   &
  in_date(6)              & ! input  date  (year,mon, day, hour, min, sec)
 ,inc_date(6)               ! increment to date

INTEGER, INTENT(OUT)  ::  &
  out_date(6)               ! out date (year,mon, day, hour, min, sec)

! local variables

INTEGER, PARAMETER ::                                                         &
   days_per_mon(12) = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

INTEGER             ::     &
  inct_yr                  & ! increment to year
 ,inct_mon                 & ! increment to month
 ,inct_day                 & ! increment to day
 ,inct_hr                  & ! increment to hour
 ,inct_min                   ! increment to minute


! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UPDATE_TIME'

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------
! update date and time

! seconds

out_date(6) = in_date(6) + inc_date(6)

IF (out_date(6) >= 60) THEN
  out_date(6) = out_date(6) - 60     ! assume input is always less than 60
  inct_min = inc_date(5) + 1
ELSE
  inct_min = inc_date(5)
END IF

! minutes

out_date(5) = in_date(5) + inct_min

IF (out_date(5) >= 60) THEN
  ! assume input is always less than 59 to start
  out_date(5) = out_date(5) - 60
  inct_hr  = inc_date(4) + 1
ELSE
  inct_hr = inc_date(4)
END IF

! hours

out_date(4) = in_date(4) + inct_hr

IF (out_date(4) >= 24) THEN
  ! assume input is always less than 24 to start
  out_date(4)  = out_date(4) - 24
  inct_day = inc_date(3) + 1
ELSE
  inct_day = inc_date(3)
END IF

! days

out_date(3) = in_date(3) + inct_day

! Depends on calendar being used - setup for real calendar but no leap year

IF (out_date(3) > days_per_mon(in_date(2)) ) THEN
  out_date(3)  = out_date(3) - days_per_mon(in_date(2))
  inct_mon = inc_date(2) + 1
ELSE
  inct_mon = inc_date(2)
END IF

! months

out_date(2) = in_date(2) + inct_mon

IF (out_date(2) > 12) THEN
  ! assume input is always less than 13 to start
  out_date(2) = out_date(2) - 12
  inct_yr  = inc_date(1) + 1
ELSE
  inct_yr  = inc_date(1)
END IF

! year first

out_date(1) = in_date(1) + inct_yr

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------

RETURN
END SUBROUTINE update_time
END MODULE update_time_mod
