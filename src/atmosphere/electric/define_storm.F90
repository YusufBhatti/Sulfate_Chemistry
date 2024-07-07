! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Purpose: Within the electric scheme, this routine defines whether
! a storm exists within each model column based on a given definition
! of 'a storm'. It then sets a 2D logical (storm_field) with the
! .TRUE. condition being used to determine a storm point.

! The current definition for a storm is based entirely
! on graupel water path (gwp) but the code is flexible to add
! other definitions in future.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Electric

MODULE define_storm_mod

USE atm_fields_bounds_mod,  ONLY: tdims
USE electric_constants_mod, ONLY: i, j, gwp_thresh

! Dr Hook modules
USE yomhook,                ONLY: lhook, dr_hook
USE parkind1,               ONLY: jprb, jpim

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DEFINE_STORM_MOD'

CONTAINS

SUBROUTINE define_storm( gwp, storm_field, nspts )

IMPLICIT NONE

REAL, INTENT(IN)       :: gwp(          tdims%i_start : tdims%i_end,           &
                                        tdims%j_start : tdims%j_end )

LOGICAL, INTENT(INOUT) :: storm_field ( tdims%i_start : tdims%i_end,           &
                                        tdims%j_start : tdims%j_end )

INTEGER, INTENT(OUT)   :: nspts ! number of storm points

!==============================================================
! Local variables
!==============================================================

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DEFINE_STORM'

!==================================================================
! Start the subroutine
!==================================================================

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! This is a very simple definition so far ...

nspts = 0

DO j = tdims%j_start, tdims%j_end

  DO i = tdims%i_start, tdims%i_end

    IF ( gwp(i,j) > gwp_thresh ) THEN

      storm_field(i,j) = .TRUE.

      nspts = nspts + 1

    END IF

  END DO ! i

END DO  ! j

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE define_storm

END MODULE define_storm_mod
