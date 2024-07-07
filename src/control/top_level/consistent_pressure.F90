! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Given the exner field at rho levels derive all other pressure fields
!
! Called by atm_step.
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Top Level

SUBROUTINE Consistent_Pressure (                                  &
           exner_rho_levels,                                      &
           offx,offy,halo_i,halo_J,                               &
           row_length,rows,model_levels,                          &
           r_theta_levels, r_rho_levels, rho,                     &
           p, pstar, p_theta_levels,exner_theta_levels)

USE planet_constants_mod, ONLY: recip_kappa, p_zero
USE mpp_conf_mod,         ONLY: swap_field_is_scalar

USE atm_fields_bounds_mod

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE Field_Types
IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels
INTEGER, INTENT(IN) :: offx
INTEGER, INTENT(IN) :: offy
INTEGER, INTENT(IN) :: halo_i
INTEGER, INTENT(IN) :: halo_j

REAL, INTENT(IN) :: exner_rho_levels                              &
                           (pdims_s%i_start:pdims_s%i_end,        &
                            pdims_s%j_start:pdims_s%j_end,        &
                            pdims_s%k_start:pdims_s%k_end + 1)

REAL, INTENT(IN) :: rho    (pdims_s%i_start:pdims_s%i_end,        &
                            pdims_s%j_start:pdims_s%j_end,        &
                            pdims_s%k_start:pdims_s%k_end)

REAL, INTENT(IN) :: r_rho_levels                                  &
                           (pdims_l%i_start:pdims_l%i_end,        &
                            pdims_l%j_start:pdims_l%j_end,        &
                            pdims_l%k_start:pdims_l%k_end)

REAL, INTENT(IN) :: r_theta_levels                                &
                           (tdims_l%i_start:tdims_l%i_end,        &
                            tdims_l%j_start:tdims_l%j_end,        &
                                          0:tdims_l%k_end)

REAL, INTENT(OUT) :: p     (pdims_s%i_start:pdims_s%i_end,        &
                            pdims_s%j_start:pdims_s%j_end,        &
                            pdims_s%k_start:pdims_s%k_end + 1)

REAL, INTENT(OUT) :: pstar(row_length,rows)

REAL, INTENT(OUT) :: p_theta_levels                               &
                           (tdims_s%i_start:tdims_s%i_end,        &
                            tdims_s%j_start:tdims_s%j_end,        &
                            tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(OUT) :: exner_theta_levels                           &
                           (tdims_s%i_start:tdims_s%i_end,        &
                            tdims_s%j_start:tdims_s%j_end,        &
                            tdims_s%k_start:tdims_s%k_end)

INTEGER :: i,j,k
LOGICAL, PARAMETER :: l_include_halos = .FALSE.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CONSISTENT_PRESSURE'


! Calculate pressure from Exner.
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
DO k = pdims_s%k_start, pdims_s%k_end + 1
  DO j = pdims_s%j_start, pdims_s%j_end
    !CDIR NODEP
    DO i = pdims_s%i_start, pdims_s%i_end
      p(i,j,k)= (exner_rho_levels(i,j,k)**recip_kappa) * p_zero
    END DO
  END DO
END DO

! Halos updated
! DEPENDS ON: swap_bounds
CALL Swap_Bounds(p,                                               &
                 row_length, rows, model_levels+1,                &
                 offx, offy, fld_type_p, swap_field_is_scalar)

! DEPENDS ON: calc_p_star
CALL Calc_P_star(                                                 &
                   r_theta_levels, r_rho_levels,                  &
                   p, rho,                                        &
                   row_length, rows, model_levels,                &
                   offx, offy, halo_i, halo_j,                    &
                   pstar )

! DEPENDS ON: calc_exner_at_theta
CALL Calc_Exner_at_theta(                                         &
                   r_theta_levels, r_rho_levels,                  &
                   exner_rho_levels,                              &
                   row_length, rows, model_levels,                &
                   offx, offy, halo_i, halo_j,                    &
                   exner_theta_levels, l_include_halos)

! Calculate pressure from Exner at theta levels.
! DEPENDS ON: calc_p_from_exner
CALL Calc_P_from_Exner(                                           &
                   p_theta_levels,                                &
                   row_length, rows,                              &
                   tdims_s%k_len,                                 &
                   offx, offy,                                    &
                   exner_theta_levels, l_include_halos)

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE Consistent_pressure
