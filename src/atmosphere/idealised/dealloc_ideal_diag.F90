! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! 
! Deallocate diagnostic arrays used by idealised UM

MODULE dealloc_ideal_diag_mod

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Deallocate diagnostic arrays required by the idealised UM
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: Fortran 95.
!   This code is written to UMDP standards.
!------------------------------------------------------------------------------

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DEALLOC_IDEAL_DIAG_MOD'

CONTAINS

! Subroutine Interface:
SUBROUTINE dealloc_ideal_diag( )

USE idealised_diag_mod, ONLY:                                             &
  dt_inc_ideal_um, dq_inc_ideal_um, du_inc_ideal_um, dv_inc_ideal_um,     &
  dtheta_inc_ideal_um, dcolqdt_ideal_um, de_cvt_ideal_um, de_u2_ideal_um, &
  de_v2_ideal_um

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

!------------------------------------------------------------------------------
! Local variables

! Variables required for Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DEALLOC_IDEAL_DIAG'

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------
! Deallocate in reverse order to allocation

DEALLOCATE ( dtheta_inc_ideal_um )
DEALLOCATE ( dv_inc_ideal_um )
DEALLOCATE ( du_inc_ideal_um )
DEALLOCATE ( dq_inc_ideal_um )
DEALLOCATE ( dt_inc_ideal_um )
DEALLOCATE ( de_v2_ideal_um )
DEALLOCATE ( de_u2_ideal_um )
DEALLOCATE ( de_cvt_ideal_um )
DEALLOCATE ( dcolqdt_ideal_um )

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE dealloc_ideal_diag

END MODULE dealloc_ideal_diag_mod
