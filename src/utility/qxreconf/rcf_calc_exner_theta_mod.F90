! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT******************************

!  Calculates exner at theta levels from values on rho levels.

MODULE Rcf_Calc_exner_theta_mod

!  Subroutine Rcf_Calc_exner_theta
!
! Description: Calculates exner at theta levels.
!
! Method: Identical to that previously used for pressure instead of
!         exner, see "A semi-Implicit scheme for the Unified Model".
!          F.R. Division working paper No 154.
!          M. J. P. Cullen and T. Davies.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_CALC_EXNER_THETA_MOD'

CONTAINS

SUBROUTINE Rcf_calc_exner_theta( &
      p_field,               &! Intent(IN) Number of points in field
      model_levels,          &! Intent(IN) Number of model levels
      r_theta_levels,        &! Intent(IN) Height at theta levels
      r_rho_levels,          &! Intent(IN) Height at rho levels
      exner_rho,             &! Intent(IN) Exner at rho levels
      exner_theta_levels     &! Intent(OUT) Exner at theta levels
                               )
USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.
INTEGER, INTENT(IN)     :: p_field       ! Number of points in field
INTEGER, INTENT(IN)     :: model_levels  ! number of model levels.

REAL, INTENT(IN)        :: exner_rho(p_field,model_levels+1)
REAL, INTENT(IN)        :: r_theta_levels(p_field,0:model_levels)
REAL, INTENT(IN)        :: r_rho_levels(p_field, model_levels)

! Arguments with Intent OUT. ie: Output variables.
REAL, INTENT(OUT)       :: exner_theta_levels(p_field, model_levels)

! Local Variables.
INTEGER                 :: i, k      ! Loop indices

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RCF_CALC_EXNER_THETA'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! --------------------------------------------------------------------
! Section 1.   Calculate pressure at desired theta levels.
! --------------------------------------------------------------------
DO k = 1, model_levels - 1
  DO i = 1, p_field
    exner_theta_levels(i,k) =                                         &
        (exner_rho(i,k) *(r_rho_levels(i,k+1) - r_theta_levels(i,k))  &
       + exner_rho(i,k+1) * (r_theta_levels(i,k) - r_rho_levels(i,k)))&
                         / (r_rho_levels(i,k+1) - r_rho_levels(i,k) )
  END DO
END DO


k = model_levels

!  extra pressure level above top theta level is same height above
!  as is the pressure below - hence weights are 0.5 and there is
!  no need to store the r_rho_level for the extra pressure

DO i = 1, p_field
  exner_theta_levels(i,k) = 0.5 * ( exner_rho(i,k) + exner_rho(i,k+1) )
END DO


! End of routine
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Calc_exner_theta

END MODULE Rcf_Calc_exner_theta_mod
