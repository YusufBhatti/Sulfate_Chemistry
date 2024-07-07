! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large scale precip (CASIM)
MODULE mphys_die

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='MPHYS_DIE'

! Integer constants for various CASIM errors
! These are used within the CASIM repository

INTEGER, PARAMETER :: incorrect_opt = 1
INTEGER, PARAMETER :: bad_values    = 2
INTEGER, PARAMETER :: warn          = -1

! Length of a standard CASIM message
! This is deliberately shorter than errormessagelength
! to allow for additional CASIM information to be added.
INTEGER, PARAMETER :: std_msg_len = errormessagelength - 150 

! Standard CASIM message, used exclusively within CASIM
! The same variable is available in mphys_die on the 
! CASIM repository for CASIM runs not made with the UM.
! It is preferable to use std_msg within the CASIM code
! but revert to using cmessage or ummessage within 
! UM routines.
CHARACTER(len=std_msg_len) :: std_msg = ''

CONTAINS

SUBROUTINE throw_mphys_error(itype, casim_routine, info)

! If modifying the subroutine or argument list, ensure that the 
! CASIM version of mphys_die (on the CASIM repository) 
! is also modified to give the same result

USE ereport_mod, ONLY: ereport
USE umprintmgr,  ONLY: newline

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER, INTENT(IN) :: itype 
! Type of error: 
! 1:  Incorrect specification of options
! 2:  Bad values found
! 3+:  Unknown error
! <0: Warning, code prints out a warning to ereport and continues

CHARACTER(LEN=*), INTENT(IN) :: casim_routine
! CASIM Routine causing the error or warning

CHARACTER(LEN=std_msg_len), INTENT(IN) :: info 
! Additional error or warning information

! Local variables
INTEGER, PARAMETER :: um_error_flag = 100
INTEGER, PARAMETER :: um_warning   = -100
! Positive value for call to ereport, negative for warning

CHARACTER(LEN=errormessagelength) :: message_str = ''
! Error or Warning message string

CHARACTER(LEN=*), PARAMETER   :: RoutineName='THROW_MPHYS_ERROR'

INTEGER :: errcode ! Error code

!----------------------------------------------------------------------------

IF ( itype > 0 ) THEN

!  --------------------------------------------
!                    ERROR
!  Code must output an error message and abort
!  --------------------------------------------

  IF ( itype == incorrect_opt ) THEN
    message_str = 'ERROR in CASIM microphysics:' //newline//&
      'Incorrect specification of options.'      //newline//&
      'Additional information:'                  //newline//&
      TRIM(info)
 
  ELSE IF ( itype == bad_values ) THEN

    message_str = 'ERROR in CASIM microphysics:' //newline//&
      'Bad values found.'                        //newline//&
      'Additional information:'                  //newline//&
      TRIM(info)

  ELSE

!   Unknown error. Just report a general error message

    message_str = 'ERROR in CASIM microphysics:' //newline//&
      'Additional information:'                  //newline//&
      TRIM(info)

  END IF ! itype

  errcode = um_error_flag

  CALL ereport(casim_routine, errcode, message_str)

ELSE IF ( itype < 0 ) THEN

!  --------------------------------------------
!                  WARNING
!  Output a warning via ereport and continue
!  --------------------------------------------

  message_str = 'WARNING from CASIM microphysics:' //newline//&
                 TRIM(info)

  errcode = um_warning

  CALL ereport(casim_routine, errcode, message_str)

END IF

END SUBROUTINE throw_mphys_error


SUBROUTINE mphys_message(casim_routine, casim_message)

USE umprintmgr,  ONLY: umPrint, newline, umMessage

IMPLICIT NONE

CHARACTER(*), INTENT(IN) :: casim_routine
! CASIM Routine causing the error

CHARACTER(*), INTENT(IN) :: casim_message 
! Additional error information

!----------------------------------------------------------------------------

WRITE(umMessage, '(A)') 'CASIM Message | '//TRIM(casim_routine)// &
                         TRIM(casim_message)

CALL umPrint( umMessage, src='mphys_message' )

END SUBROUTINE mphys_message

END MODULE mphys_die
