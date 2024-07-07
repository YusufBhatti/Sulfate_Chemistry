! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Set all tile fractions to bare soil

MODULE Rcf_init_tile_frac_mod

!  Subroutine Rcf_init_tile_frac
!
! Description:
!   Set all tile fractions to bare soil
!
! Method:
!   Over-write the input tile fractions with 0 for all tiles except
!   bare soil, which is set to 1
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_INIT_TILE_FRAC_MOD'

CONTAINS

SUBROUTINE Rcf_init_tile_frac(current_field_out)

USE jules_surface_types_mod, ONLY: soil

USE Rcf_Field_Type_Mod, ONLY: field_type

USE yomhook,   ONLY: lhook, dr_hook

USE parkind1,  ONLY: jprb, jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ), INTENT(INOUT)   :: current_field_out

! Local Variables/Parameters
INTEGER               :: i                   ! Looper
INTEGER               :: k                   ! Looper 2

CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_INIT_TILE_FRAC'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!--------------------------------------------------------------------
! Set everything to be tile 8 (bare soil)
!--------------------------------------------------------------------
DO k = 1, current_field_out % levels

  IF ( k == soil ) THEN

    DO i = 1, current_field_out % level_size
      current_field_out % Data( i, k ) = 1.0
    END DO

  ELSE !all other surface types

    DO i = 1, current_field_out % level_size
      current_field_out % Data( i, k ) = 0.0
    END DO

  END IF !surface type

END DO !surface types
      
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_init_tile_frac

END MODULE Rcf_init_tile_frac_mod
