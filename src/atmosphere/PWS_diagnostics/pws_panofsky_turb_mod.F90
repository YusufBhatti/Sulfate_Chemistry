! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Purpose: Calculates Colson-Panofsky turbulence, stash section 20
!
!   Programming standard; Unified Model Documentation Paper No. 3
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS_diagnostics

MODULE pws_colson_panofsky_turb_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='PWS_COLSON_PANOFSKY_TURB_MOD'

CONTAINS

SUBROUTINE pws_colson_panofsky_turb(u, v, theta, p_theta_levels)


USE yomhook,                    ONLY: lhook, dr_hook
USE parkind1,                   ONLY: jprb, jpim

USE pws_diags_mod,              ONLY: panofsky_turb_press_levels,            &
                                      panofsky_turb_press,                   &
                                      pws_panofsky_turb

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
REAL :: wind_shear_squared_plevs(row_length, rows, panofsky_turb_press_levels)

! Length scale
REAL :: dz(row_length, rows, panofsky_turb_press_levels)

! Brunt-Vaisala Frequency on pressure levels
REAL :: bv_freq(row_length,                                                  &
                rows,                                                        &
                panofsky_turb_press_levels)

! Masks - false if the point is NOT valid
LOGICAL :: mask_windshear(row_length, rows, panofsky_turb_press_levels)
LOGICAL :: mask_bvfreq(row_length, rows, panofsky_turb_press_levels)


! Interpolation:
! Model level corresponding to pressure level we're interested in
INTEGER :: level
! Fraction of the way we are through the layer
REAL :: alpha

REAL, PARAMETER :: tuning_parameter = 2.0

INTEGER :: i,j,l ! Loop counters
INTEGER :: ErrorStatus

CHARACTER (LEN=*), PARAMETER :: RoutineName='PWS_COLSON_PANOFSKY_TURB'

! dr_hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (.NOT. l_regular) THEN
  ErrorStatus=100
  CALL EReport( RoutineName, ErrorStatus,                                    &
                "Cannot calculate Colson-Panofsky turbulence - " //          &
                "var res and/or rotated grids are not supported" )
END IF


! Calculate vertical wind shear squared
CALL pws_vert_wind_shear_sq(panofsky_turb_press_levels, panofsky_turb_press, &
                        u, v, p_theta_levels, wind_shear_squared_plevs,      &
                        mask_windshear)

! Calculate Brunt-Vaisala Frequency and lengthscale
CALL pws_brunt_vaisala_freq(panofsky_turb_press_levels, panofsky_turb_press, &
                        theta, p_theta_levels, bv_freq, dz, mask_bvfreq)


!   Calculate Colson-Panofsky Index using
!      CP_Index = (length_scale^2) * (S_sq) * (1 - (N_sq/S_sq)/ Ri_crit)
!   
!   Ref: Colson and Panofsky (1965), 'An index of clear-air turbulence.' 
!        Quart. J.Roy.Meteor.Soc., 91, pp 507-513.
!

pws_panofsky_turb = rmdi

DO l = 1, panofsky_turb_press_levels
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      IF (mask_bvfreq(i,j,l) .AND. mask_windshear(i,j,l)) THEN
        pws_panofsky_turb(i,j,l) = dz(i,j,l)**2 *                            &
               wind_shear_squared_plevs(i,j,l) *                             &
                 (1.0 - ( tuning_parameter * (bv_freq(i,j,l)/                &
                  wind_shear_squared_plevs(i,j,l))))
      END IF
    END DO
  END DO
END DO



IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE pws_colson_panofsky_turb


END MODULE pws_colson_panofsky_turb_mod
