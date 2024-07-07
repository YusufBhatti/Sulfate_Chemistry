#if defined(MCT)
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************

SUBROUTINE oasis3_advance_date()


! Description: This routine advances the prism date and calculates
!              the appropriate time for use in  put/get operations.
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Coupling
!=============================================================


USE um_types

USE OASIS3_atmos_init_mod, ONLY: prism_nsec, prism_timestep
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='OASIS3_ADVANCE_DATE'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! Advance the date/time as used by prism
prism_nsec = prism_nsec + NINT(prism_timestep)

!=======================================================================

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE oasis3_advance_date
#endif
