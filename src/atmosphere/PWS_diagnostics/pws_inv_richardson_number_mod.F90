! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Purpose: Calculates Inverse Richardson number, stash section 20
!
!   Programming standard; Unified Model Documentation Paper No. 3
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS_diagnostics

MODULE pws_inv_richardson_number_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                      &
                         ModuleName='PWS_INV_RICHARDSON_NUMBER_MOD'

CONTAINS

SUBROUTINE pws_inv_richardson_number(u, v, theta, p_theta_levels)


USE yomhook,                    ONLY: lhook, dr_hook
USE parkind1,                   ONLY: jprb, jpim

USE pws_diags_mod,              ONLY: inv_richardson_press_levels,           &
                                      inv_richardson_press,                  &
                                      pws_inv_richardson

USE atm_fields_bounds_mod

USE model_domain_mod,           ONLY: model_type, mt_lam, l_regular
USE ereport_mod,                ONLY: ereport
USE missing_data_mod,           ONLY: rmdi, imdi
USE nlsizes_namelist_mod,       ONLY: row_length, rows, model_levels
USE pws_brunt_vaisala_freq_mod, ONLY: pws_brunt_vaisala_freq
USE pws_vert_wind_shear_sq_mod, ONLY: pws_vert_wind_shear_sq

IMPLICIT NONE

REAL, INTENT(IN)  :: u(udims_s%i_start:udims_s%i_end,                        &
                       udims_s%j_start:udims_s%j_end,                        &
                       udims_s%k_start:udims_s%k_end)

REAL, INTENT(IN)  :: v(vdims_s%i_start:vdims_s%i_end,                        &
                       vdims_s%j_start:vdims_s%j_end,                        &
                       vdims_s%k_start:vdims_s%k_end)

REAL, INTENT(IN)  :: theta(tdims_s%i_start:tdims_s%i_end,                    &
                           tdims_s%j_start:tdims_s%j_end,                    &
                           tdims_s%k_start:tdims_s%k_end)


REAL, INTENT(IN)  :: p_theta_levels(tdims_s%i_start:tdims_s%i_end,           &
                                    tdims_s%j_start:tdims_s%j_end,           &
                                    tdims_s%k_start:tdims_s%k_end)


! Local variables

! Wind shear squared on pressure levels
REAL :: wind_shear_squared_plevs(row_length, rows, inv_richardson_press_levels)

! Length scale
REAL :: dz(row_length, rows, inv_richardson_press_levels)

! Brunt-Vaisala Frequency on pressure levels
REAL :: bv_freq(row_length, rows, inv_richardson_press_levels)

! Masks - false if the point is NOT valid
LOGICAL :: mask_windshear(row_length, rows, inv_richardson_press_levels)
LOGICAL :: mask_bvfreq(row_length, rows, inv_richardson_press_levels)


INTEGER :: i,j,l ! Loop counters
INTEGER :: ErrorStatus

REAL, PARAMETER :: small_number = TINY(1.0)

CHARACTER (LEN=*), PARAMETER :: RoutineName='PWS_INV_RICHARDSON_NUMBER'

! dr_hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (.NOT. l_regular) THEN
  ErrorStatus=100
  CALL EReport( RoutineName, ErrorStatus,                                    &
                "Cannot calculate Inverse Richardson number - " //           &
                "var res and/or rotated grids are not supported" )
END IF


! Calculate vertical wind shear squared
CALL pws_vert_wind_shear_sq(inv_richardson_press_levels,                     &
                            inv_richardson_press, u, v, p_theta_levels,      &
                            wind_shear_squared_plevs, mask_windshear)

! Calculate Brunt-Vaisala Frequency and lengthscale
CALL pws_brunt_vaisala_freq(inv_richardson_press_levels,                     &
                            inv_richardson_press, theta, p_theta_levels,     &
                            bv_freq, dz, mask_bvfreq)


! Finally calculate the Inverse Richardson Number
pws_inv_richardson = rmdi

DO l = 1, inv_richardson_press_levels
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      IF (mask_windshear(i,j,l) .AND. mask_bvfreq(i,j,l) .AND.               &
          bv_freq(i,j,l) > small_number) THEN
        pws_inv_richardson(i,j,l) = wind_shear_squared_plevs(i,j,l) /        &
                                    bv_freq(i,j,l)
      END IF
    END DO
  END DO
END DO


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE pws_inv_richardson_number


END MODULE pws_inv_richardson_number_mod
