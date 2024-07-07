! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Top level control

!  Description:

! Fortran module to provide subroutines required for the population of the
! namelist nlstcall

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

! Code Description:
!   Language: FORTRAN 90

MODULE nlstcall_subs_mod

IMPLICIT NONE

INTEGER, PARAMETER :: hours = 1
INTEGER, PARAMETER :: days = 2
INTEGER, PARAMETER :: timesteps = 3
INTEGER, PARAMETER :: real_months = 4

CONTAINS

SUBROUTINE totime(intime,UNIT,outtime)

  ! Subroutine is a direct replacement for the UMUI totime function
  ! The subroutine takes a time (intime) in units (unit) of hours
  ! 'H', days 'D', or timesteps 'T', and returns that time in units of
  ! timesteps (outtime).

USE nlstgen_mod, ONLY: secs_per_periodim, &
                        steps_per_periodim
USE conversions_mod, ONLY: isec_per_day

IMPLICIT NONE

INTEGER :: intime,outtime
INTEGER :: days_per_period
INTEGER :: UNIT

days_per_period = secs_per_periodim(1) / isec_per_day

IF (UNIT == hours) THEN
  outtime = (intime*steps_per_periodim(1)) / (days_per_period * 24)
ELSE IF (UNIT == days) THEN
  outtime = (intime*steps_per_periodim(1)) / days_per_period
ELSE IF (UNIT == timesteps) THEN
  outtime = intime
END IF

END SUBROUTINE totime

END MODULE nlstcall_subs_mod
