! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routine: EXITCHEK -------------------------------------------------
!
! Purpose: Checks for end-of-run condition and returns a logical to
!          the top level.  There is now only one reason for stopping,
!          namely:
!          (i)   Model has completed the required integration;
!
! Programming standard: UM Doc Paper 3
!
! Logical components covered: C0
!
! Project task: C0
!
! External documentation: On-line UM document C0 - The top-level
!                         control system
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

SUBROUTINE exitchek                                               &
        ( internal_model, lexitNOW )

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE Control_Max_Sizes
USE submodel_mod, ONLY: atmos_im
USE model_time_mod, ONLY: &
    stepim, target_end_stepim

IMPLICIT NONE

INTEGER :: internal_model ! In  - id of current internal model
LOGICAL ::   lexitNOW  ! Out - True/False flag for stopping
!
! ----------------------------------------------------------------------
!  Local variables
!
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EXITCHEK'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

!
! ----------------------------------------------------------------------
!  1. Check for completed run
!
lexitNOW = ( stepim(atmos_im) >= target_end_stepim(atmos_im) )
! ----------------------------------------------------------------------
!
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE exitchek
