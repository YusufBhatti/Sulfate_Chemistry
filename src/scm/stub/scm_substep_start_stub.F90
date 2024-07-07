! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Signal start of sub-step to SCM output diagnostic system

SUBROUTINE scm_substep_start(substep_number)

IMPLICIT NONE

! Description:
!   In order to be able to output diagnostics that are defined
!   within sub-stepped parts of the model, the diagnostic system
!   needs to know which parts of the model are sub-stepped and
!   which sub-step the model is currently on, if any. This routine
!   is called to let the system know that the specified sub-step
!   has just started. The system will assume that it is on the
!   specified sub-step until this routine is called again to
!   indicate the start of the next sub-step, or SCM_SUBSTEPPING_END
!   is called to indicate that the sub-stepping has now ended. Any
!   call to SCMoutput from a sub-stepped part of the code that is
!   not delimited by calls to these two routines (the first inside
!   the sub-stepping loop and the second outside once it has
!   finished) will result in an error message.

! Method:

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model


! Code Description:
!   Language: Fortran90

INTEGER :: substep_number

RETURN

END SUBROUTINE scm_substep_start

