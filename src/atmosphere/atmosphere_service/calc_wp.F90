! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE calc_wp_mod

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

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CALC_WP_MOD'

CONTAINS

SUBROUTINE calc_wp( wc, rhodz, wp_tot )

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
! Air density * depth of model level

REAL, INTENT(OUT) :: wp_tot(tdims%i_start : tdims%i_end,                      &
                            tdims%j_start : tdims%j_end)
! Output water path

!-------------------
! Local variables
!-------------------

INTEGER :: i, j, k ! Loop indices

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_WP'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

wp_tot(:,:) = 0.0

DO k = 1, tdims%k_end

  DO j = tdims%j_start, tdims%j_end

    DO i = tdims%i_start, tdims%i_end

      IF (wc(i,j,k) > 0.0) THEN ! only include if water content is positive

        wp_tot(i,j) = wp_tot(i,j) + (rhodz(i,j,k) * wc(i,j,k) )

      END IF

    END DO ! i

  END DO ! j

END DO ! k

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE calc_wp

END MODULE calc_wp_mod

