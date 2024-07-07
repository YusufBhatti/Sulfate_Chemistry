! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

SUBROUTINE add_period_to_date(  year,  month,  day,  hour,  minute,  secs, &
                              i_year,i_month,i_day,i_hour,i_minute,i_secs, &
                              lcal360 )

USE cdaydata_mod, ONLY: days_per_4c, days_per_c, days_per_4y, days_per_y,  &
                        days_in_month, days_to_month
USE umPrintMgr, ONLY: umPrint, umMessage
USE errormessagelength_mod, ONLY: errormessagelength
USE ereport_mod, ONLY: ereport

IMPLICIT NONE


! Routine: add_period_to_date
!
! Purpose: Adds a period to a date, returning the new date.
!
! Supports Gregorian and 360-day calendar.
!
! Dates should be valid dates in the Common Era.
! Periods can be expressed in any positive number of each unit.
! For example, the i_secs increment is not limited to 0-59.
! No negative values are allowed in the period field
!
! Return icode: 0 - success. 1 - date invalid. 2 - increment invalid
!
! Gregorian calendar issues:
! The following examples illustrate how variable month lengths are handled
! When adding 1 month to 15th March, the result is 15th April.
! When adding 1 month to 31st March, the result is 1st May because there
! is no 31st April.
! When adding 1 month to 31st January, the result is 2nd or 3rd March
! depending on whether it is a leap year or not.
! When adding a year to 29th February, the result is 1st March.
!
! Programming standard: UM Doc Paper 3, version 8.2 (25/3/2009)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
! External documentation: On-line UM document C0 - The top-level
!                         control system
!

! Input date and output result
INTEGER, INTENT(INOUT) :: year
INTEGER, INTENT(INOUT) :: month
INTEGER, INTENT(INOUT) :: day
INTEGER, INTENT(INOUT) :: hour
INTEGER, INTENT(INOUT) :: minute
INTEGER, INTENT(INOUT) :: secs

! Time period to add
INTEGER, INTENT(IN)    :: i_year
INTEGER, INTENT(IN)    :: i_month
INTEGER, INTENT(IN)    :: i_day
INTEGER, INTENT(IN)    :: i_hour
INTEGER, INTENT(IN)    :: i_minute
INTEGER, INTENT(IN)    :: i_secs

! Is this a 360-day calendar or a Gregorian calendar
LOGICAL, INTENT(IN)    :: lcal360

! Local variables

INTEGER :: days_in_month_this_year(12)
INTEGER :: icode
LOGICAL :: valid_date

CHARACTER(LEN=18)  :: routinename = 'ADD_PERIOD_TO_DATE'
CHARACTER(LEN=errormessagelength) :: cmessage

! cdaydata.h needed for integer, parameter :: days_in_month(12)

icode = 0

! Check input date is valid
valid_date = .TRUE.
IF (year < 1 .OR. month < 1 .OR. month > 12 .OR. hour < 0 .OR. hour > 23 &
   .OR. minute < 0 .OR. minute > 59 .OR. secs < 0 .OR. secs > 59     &
   .OR. day < 0) THEN
  ! Illegal input date
  valid_date = .FALSE.
END IF
IF (lcal360) THEN
  IF (day > 30) THEN
    valid_date = .FALSE.
  END IF
ELSE
  days_in_month_this_year = days_in_month
  IF (MOD(year,4) == 0   .AND.                                    &
        (MOD(year,400) == 0 .OR. MOD(year,100) /= 0)) THEN
    days_in_month_this_year(2) = 29
  END IF

  IF (day > days_in_month_this_year(month)) THEN
    valid_date = .FALSE.
  END IF
END IF

IF (.NOT. valid_date) THEN
  WRITE(cmessage,'(A,I5,5I3)')                                           &
    'Input date is invalid: Date ',                                      &
    year, month, day, hour, minute, secs
  icode = 1
  GO TO 9999
END IF

! Check increment date contains no negative values
IF (i_year < 0 .OR. i_month < 0 .OR. i_day < 0 .OR. i_hour < 0 .OR. &
    i_minute < 0 .OR. i_secs < 0) THEN
  WRITE(cmessage,'(A,I5,5I3)')                                           &
    'Negative increments not allowed: Increment ',                       &
     i_year, i_month, i_day, i_hour, i_minute, i_secs
  icode = 2
  GO TO 9999
END IF

! Add each element of the increment to the corresponding date element
secs = secs + i_secs
minute = minute + i_minute
hour = hour + i_hour
day = day + i_day
month = month + i_month
year = year + i_year

! Round down each period in turn
IF (secs >= 60) THEN
  minute = minute + (secs/60)
  secs = MOD(secs,60)
END IF

IF (minute >= 60) THEN
  hour = hour + (minute/60)
  minute = MOD(minute,60)
END IF

IF (hour >= 24) THEN
  day = day + (hour/24)
  hour = MOD(hour,24)
END IF

! Days can be dealt with now in a 360-day calendar,
! but Gregorian calendar needs months to be considered before days.
IF (lcal360) THEN
  ! 360-day calendar
  IF (day > 30) THEN
    ! Example: If day=60 and month=1 we would want day=30 and month=2
    ! not day=0 and month=3
    month = month + ((day-1)/30)
    day = MOD(day,30)
    IF (day == 0) THEN
      day = 30
    END IF
  END IF
END IF

IF (month > 12) THEN
  ! Example:
  ! If month=24 and year=2001 we would want month=12 and year=2002
  ! not month=0 and year=2003
  year = year + ((month-1)/12)
  month = MOD(month,12)
  IF (month == 0) THEN
    month = 12
  END IF
END IF


IF (.NOT. lcal360) THEN
  ! Gregorian calendar
  days_in_month_this_year = days_in_month
  IF (MOD(year,4) == 0   .AND.                                    &
         (MOD(year,400) == 0 .OR. MOD(year,100) /= 0)) THEN
    days_in_month_this_year(2) = 29
  END IF

  ! Increment days and decrement months till less than one month
  ! worth of days, taking into account length of each month
  DO WHILE (day > days_in_month_this_year(month))
    day = day - days_in_month_this_year(month)
    month = month + 1

    ! End of year, need to recalculate whether next year is leap year
    IF (month == 13) THEN
      year = year + 1
      month = 1
      IF (MOD(year,4) == 0   .AND.                                    &
             (MOD(year,400) == 0 .OR. MOD(year,100) /= 0)) THEN
        days_in_month_this_year(2) = 29
      ELSE
        days_in_month_this_year(2) = days_in_month(2)
      END IF
    END IF
  END DO
END IF

9999 CONTINUE
IF (icode  /=  0) THEN
  CALL ereport(routinename, icode, cmessage)
END IF

RETURN
END SUBROUTINE add_period_to_date
