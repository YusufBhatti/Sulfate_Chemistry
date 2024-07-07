! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE Rcf_Calc_2D_CCA_Mod

!  Subroutine Rcf_Calc_2D_CCA - calculates 2D cca from 3D cca
!
! Description: Calculates a 2D convective cloud amount from the 2D
!              convective cloud amount array.
!
! Method: If an anvil is detected (ie there is more convective cloud
!         at the top than the base), the base value is divided by
!         the tower_factor. Otherwise, the value is the one at
!         the cloud base. Note that the tower_factor used to genenerate
!         the 3D field is not available, so the most common one
!         in practice is used.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!   Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Reconfiguration

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_CALC_2D_CCA_MOD'

CONTAINS

SUBROUTINE Rcf_Calc_2D_CCA( cca_3d, ccb, cct, cca_2d )

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ), INTENT(IN)    :: cca_3d
TYPE( field_type ), INTENT(IN)    :: ccb
TYPE( field_type ), INTENT(IN)    :: cct
TYPE( field_type ), INTENT(INOUT) :: cca_2d

! Local variables
INTEGER                           :: i          ! Looper

! The most commonly used value for tower_factor is set as a parameter
! as it is not supplied to the reconfiguration
REAL, PARAMETER                   :: tower_factor = 0.25

  ! min_val is the amount of difference in convective cloud which
  ! constitutes an anvil.
  ! A two stage solution is used as some compilers will not promote parameters
  ! when '-r8' is used. This ensures that TINY() is passed a variable of the
  ! default real kind.
REAL            :: arbitraryReal = 0.1
REAL, PARAMETER :: min_val = TINY(arbitraryReal)

! max_2d_cca is the maximum amount of 2d cca permitted
REAL, PARAMETER                   :: max_2d_cca = 0.5

CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_CALC_2D_CCA'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! Loop over all points
DO i=1, cca_3d % level_size

  ! Do we have any cca to convert?
  IF ( ccb % Data_Int(i,1) == 0 .AND. cct % Data_Int(i,1) == 0) THEN

    cca_2d % DATA(i,1) = 0.0

    ! Do we have an anvil?
  ELSE IF ((cca_3d % DATA( i, cct % Data_Int(i,1) - 1) -        &
       cca_3d % DATA( i, ccb % Data_Int(i,1))) >  min_val) THEN
    cca_2d % DATA(i,1) = cca_3d % DATA(i, ccb % Data_Int(i,1))  &
                       / tower_factor

  ELSE
    cca_2d % DATA(i,1) = cca_3d % DATA(i, ccb % Data_Int(i,1))
  END IF

  ! Ensure that the maximum 2d cca isn't breached.
  IF (cca_2d % DATA(i,1) > max_2d_cca) THEN
    cca_2d % DATA(i,1) = max_2d_cca
  END IF

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Calc_2D_CCA
END MODULE Rcf_Calc_2D_CCA_Mod
