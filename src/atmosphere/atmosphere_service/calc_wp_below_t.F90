! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE calc_wp_below_t_mod

! Purpose: calculates any water path in units of kg m-2

! This is similar in nature to the calculations performed
! for Section 30 (Climate diagnostics), the code of which can
! be found in vert_eng_massq (within energy_correction)
! However, the code in calc_wp is much cheaper than vert_eng_massq
! as it only performs one integral for one liquid species only.
! Therefore calc_wp is recommended for simple calculations involving
! one or two moist variables.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Atmosphere Service

USE atm_fields_bounds_mod, ONLY: tdims

! Dr Hook modules
USE yomhook,               ONLY: lhook, dr_hook
USE parkind1,              ONLY: jprb, jpim

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CALC_WP_BELOW_T_MOD'

CONTAINS

SUBROUTINE calc_wp_below_t( wc, rhodz, t, t_thresh, wp_tot )

IMPLICIT NONE

!-------------------------------
! Subroutine arguments
!-------------------------------

REAL, INTENT(IN)::    wc(tdims%i_start : tdims%i_end,                         &
                         tdims%j_start : tdims%j_end,                         &
                                     1 : tdims%k_end)
! Any water mixing ratio (q, qcl, qcf, qrain, qgraup)- units kg/kg

REAL, INTENT(IN):: rhodz(tdims%i_start : tdims%i_end,                         &
                         tdims%j_start : tdims%j_end,                         &
                                     1 : tdims%k_end)
! Air density * depth of model level (kg m-2)

REAL, INTENT(IN):: t(tdims%i_start : tdims%i_end,                             &
                     tdims%j_start : tdims%j_end,                             &
                                 1 : tdims%k_end)
! Temperature (K)

REAL :: t_thresh
! Temperature threshold (K)

REAL, INTENT(OUT) :: wp_tot(tdims%i_start : tdims%i_end,                      &
                            tdims%j_start : tdims%j_end)

!-------------------
! Local variables
!-------------------

INTEGER :: i, j, k ! Loop indices

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_WP_BELOW_T'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

wp_tot(:,:) = 0.0

DO k = 1, tdims%k_end

  DO j = tdims%j_start, tdims%j_end

    DO i = tdims%i_start, tdims%i_end

      IF (wc(i,j,k) > 0.0 .AND. t(i,j,k) <= t_thresh) THEN
         ! only include if water content is positive and temperature
         ! below threshold

        wp_tot(i,j) = wp_tot(i,j) + (rhodz(i,j,k) * wc(i,j,k) )

      END IF

    END DO ! i

  END DO ! j

END DO ! k

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE calc_wp_below_t

END MODULE calc_wp_below_t_mod

