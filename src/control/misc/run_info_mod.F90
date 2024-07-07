! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Misc.
!
! Variables for full UM run information

MODULE run_info

USE get_wallclock_time_mod, ONLY: get_wallclock_time
 
IMPLICIT NONE

PRIVATE

REAL, PROTECTED, PUBLIC :: start_time  !! time at start of UM_SHELL

REAL, PUBLIC :: time_start_atmstep !! time at start of current call to ATMSTEP

PUBLIC :: set_start_time

CONTAINS

SUBROUTINE set_start_time()

USE umPrintMgr, ONLY: umPrint
USE IOS_common, ONLY: IOS_start_time

IMPLICIT NONE


LOGICAL :: after_first_call = .FALSE.

IF (after_first_call) THEN

  CALL umPrint('Start time already set. Doing nothing.',src='set_start_time')

ELSE

! Not done via parallel region as start_time/ios_start_time are shared
  start_time =  get_wallclock_time()
  IOS_start_time = start_time
  after_first_call = .TRUE.

END IF

END SUBROUTINE set_start_time

END MODULE run_info
