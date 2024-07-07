! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
! ALLOCATE full grid arrays
MODULE reset_field_flags_mod

IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RESET_FIELD_FLAGS_MOD'

CONTAINS

SUBROUTINE reset_field_flags( )

USE field_flags_mod, ONLY:                                                  &
   l_u, l_v, l_w, l_theta, l_ptheta, l_q, l_qcl, l_qcf, l_qrain, l_qgraup,  &
   l_dt1, l_dt2, l_dt4, l_dt9, l_dt12, l_dt30, l_dq4, l_dq9, l_dq12, l_dq30,&
   l_dqcl4, l_dqcl9, l_dqcl12, l_dqcl30,                                    &
   l_dqcf4, l_dqcf3, l_dqcf12, l_dqcf30, l_got_fields,                      &
   l_rain, l_snow, l_precip, l_sh, l_lh, l_zh, l_tstar, l_pstar

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!   ALLOCATE and initialise  arrays required by convection for a full model
!   timestep on the first convection substep of a model timestep.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utility - crmstyle_coarse_grid
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3
! ------------------------------------------------------------------------------
! Subroutine arguments  - None

!-------------------------------------------------------------------------------
! Local variables

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RESET_FIELD_FLAGS'

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------
! reset all field flags to false for new time

l_u  = .FALSE.
l_v  = .FALSE.
l_w  = .FALSE.

l_theta  = .FALSE.
l_q       = .FALSE.
l_qcl     = .FALSE.
l_qcf     = .FALSE.
l_qrain   = .FALSE.
l_qgraup  = .FALSE.
l_ptheta  = .FALSE.

l_dt1  = .FALSE.
l_dt2  = .FALSE.
l_dt4  = .FALSE.
l_dt9  = .FALSE.
l_dt12  = .FALSE.
l_dt30  = .FALSE.

l_dq4  = .FALSE.
l_dq9  = .FALSE.
l_dq12  = .FALSE.
l_dq30  = .FALSE.

l_dqcl4  = .FALSE.
l_dqcl9  = .FALSE.
l_dqcl12  = .FALSE.
l_dqcl30  = .FALSE.

l_dqcf4  = .FALSE.
l_dqcf3  = .FALSE.
l_dqcf12  = .FALSE.
l_dqcf30  = .FALSE.

l_rain = .FALSE.
l_snow = .FALSE.
l_precip = .FALSE.
l_sh = .FALSE.
l_lh = .FALSE.
l_zh = .FALSE.
l_tstar = .FALSE.
l_pstar = .FALSE.

l_got_fields = .FALSE.
!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------
RETURN
END SUBROUTINE reset_field_flags
END MODULE reset_field_flags_mod
