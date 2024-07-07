! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Returns the day of the week (Sunday=1, ..., Saturday=7) given the date.
!
! Method:
!   Based on algorithm from Calendar FAQ, v. 2.6
!   (modified 24 June 2003) Part 1/3
!   http://www.faqs.org/faqs/calendars/faq/part1/index.html
!   Assumes Gregorian calendar
!   a = (14 - month) / 12
!   y = year - a
!   m = month + 12*a - 2
!   For Julian calendar:
!      d = (5 + day + y + y/4 + (31*m)/12) mod 7
!   For Gregorian calendar:
!      d = (day + y + y/4 - y/100 + y/400 + (31*m)/12) mod 7
!   We add 1 to the result so that Sunday=1, ..., Saturday=7. This way
!   we can use it to index arrays with lower and upper bounds equal to
!   1 and 7, respectively.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
! Code description:
!   Language: Fortran 95.
!   This code is written to UMDP3 programming standards.
!
! --------------------------------------------------------------------------
!
MODULE day_of_week_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DAY_OF_WEEK_MOD'

CONTAINS

INTEGER FUNCTION day_of_week (i_day, i_month, i_year)

USE nlstcall_mod,          ONLY: lcal360
USE ereport_mod,           ONLY: ereport

USE parkind1,              ONLY: jpim,  jprb  ! DrHook
USE yomhook,               ONLY: lhook, dr_hook

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: i_day, i_month, i_year

! Local Variables
INTEGER :: a, y , m
INTEGER :: ierror = 0  ! error flag (0 = OK)
INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1

REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DAY_OF_WEEK'

CHARACTER (LEN=errormessagelength) :: cmessage ! error message for ereprot


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Report error if the model uses 360-day calendar
IF (lcal360) THEN
  ierror   = 1
  cmessage = 'This routine only works for a Gregorian calendar'
  CALL ereport ('DAY_OF_WEEK', ierror, cmessage)
END IF

! See header for explanation of the calculations
a = (14 - i_month) / 12
y = i_year - a
m = i_month + 12*a - 2
day_of_week = MOD ((i_day + y + y/4 - y/100 + y/400 + (31*m)/12), 7) + 1

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END FUNCTION day_of_week

END MODULE day_of_week_mod
