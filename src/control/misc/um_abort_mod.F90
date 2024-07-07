! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Misc.

MODULE um_abort_mod

USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: ERROR_UNIT

IMPLICIT NONE 

PRIVATE

PUBLIC :: set_abort_handler, um_abort

! Interface for an abort handler
ABSTRACT INTERFACE
SUBROUTINE handler_sub(status)
IMPLICIT NONE
INTEGER, INTENT(IN) :: status
END SUBROUTINE handler_sub
END INTERFACE

! Pointer to the current handler which will be called when calling um_abort
PROCEDURE(handler_sub), POINTER :: abort_handler => NULL()

CONTAINS

! Called to update the target for the handler, to provide a custom abort
! function to be called by um_abort
!-------------------------------------------------------------------------------
SUBROUTINE set_abort_handler(handler)
IMPLICIT NONE 
PROCEDURE(handler_sub) :: handler
abort_handler => handler
END SUBROUTINE set_abort_handler

! Called from the code to abort the program; will call whichever handler is
! currently set as the abort_handler
!-------------------------------------------------------------------------------
SUBROUTINE um_abort(status)
IMPLICIT NONE 
INTEGER, INTENT(IN) :: status
IF (ASSOCIATED(abort_handler)) THEN
  CALL abort_handler(status)
ELSE
  ! If no handler was been registered, STOP the program; note this should never
  ! be called by the UM or any other application which utilises multiple
  ! processes (its inclusion is purely for libraries using this module)
  WRITE(ERROR_UNIT, "(A)") &
    "um_abort_mod: NO ABORT HANDLER SET, CALLING 'STOP'"
  STOP 
END IF
END SUBROUTINE um_abort
!-------------------------------------------------------------------------------

END MODULE um_abort_mod


