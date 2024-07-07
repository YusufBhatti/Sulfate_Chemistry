! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Set parameters used in parametrization of fractional standard
! deviation (FSD) of subgrid water content.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!------------------------------------------------------------------------------
MODULE set_fsd_parameters_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SET_FSD_PARAMETERS_MOD'

CONTAINS

SUBROUTINE set_fsd_parameters(                                                 &
! Model dimensions.
  row_length, rows, n_cca_levels, offx, offy,                                  &
! Properties of clouds
  cca_dp, cca_md, cca_sh, xx_cos_theta_latitude,                               &
! Model switches
  l_rad_step_diag, l_rad_step_prog, at_extremity,                              &
! Error information
  Error_code  )

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE level_heights_mod, ONLY: r_theta_levels
USE mphys_inputs_mod, ONLY: l_fsd_generator
USE atm_fields_bounds_mod, ONLY: tdims
USE field_types, ONLY: fld_type_p
USE UM_ParParams, ONLY: PNorth, PSouth
USE mpp_conf_mod, ONLY: swap_field_is_scalar

USE model_domain_mod, ONLY: model_type, mt_global
USE rad_input_mod, ONLY: l_fsd_eff_res, i_fsd
USE cderived_mod, ONLY: delta_lambda, delta_phi
USE conversions_mod, ONLY: pi_over_180
USE fsd_parameters_mod, ONLY: f_arr, f_arr_c, f_cons,                          &
  ip_fsd_constant, ip_fsd_param, ip_fsd_regime, ip_fsd_regime_no_sh,           &
  ip_fsd_boutle, ip_fsd_regime_smooth, ip_fsd_regime_smooth_no_sh,             &
  fsd_eff_lam, fsd_eff_phi

IMPLICIT NONE


INTEGER, INTENT(INOUT) :: Error_code
!   Error flag

! Dimensions of arrays
INTEGER, INTENT(IN) ::                                                         &
  row_length,                                                                  &
!   Number of points on a row
  rows,                                                                        &
!   Number of rows
  n_cca_levels,                                                                &
!   Number of convective cloud amount levels
  offx,                                                                        &
!   Size of small halo in i
  offy
!   Size of small halo in j

! Convective cloud cover
REAL, INTENT(IN) ::                                                            &
  cca_dp(row_length, rows, n_cca_levels),                                      &
!   Convective cloud amount from deep scheme
  cca_md(row_length, rows, n_cca_levels),                                      &
!   Convective cloud amount from mid-level scheme
  cca_sh(row_length, rows, n_cca_levels),                                      &
!   Convective cloud amount from shallow sheme
  xx_cos_theta_latitude(1-offx:row_length+offx,1-offy:rows+offy)
!   Finite volume cosine of latitude.

! Model Switches
LOGICAL, INTENT(IN) ::                                                         &
  l_rad_step_diag,                                                             &
!   True if fast radiation timestep
  l_rad_step_prog,                                                             &
!   True if slow radiation timestep
  at_extremity(4)
!   Indicates if this processor is at north, south
!   east or west of the processor grid


! Local variables.
REAL, ALLOCATABLE :: conv_ind(:,:,:)
!   Convection indicator derived from convective cloud
!   fraction. Required for regime dependent FSD parametrization.

REAL, ALLOCATABLE :: conv_ind_temp(:,:,:)
!   Convection indicator derived from convective cloud
!   fraction. Temporary value required for smoothing

REAL ::                                                                        &
  conv_thick_part,                                                             &
!   Part of FSD parametrization that may vary with layer
!   thickness and regime.
  x_in_km
!   Grid-box size in km

INTEGER :: i,j,k
!   Loop counters

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SET_FSD_PARAMETERS'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (l_fsd_eff_res) THEN
  ! USE effective resolution in FSD parametrization
  ! N96 - 1.875 x 1.25 degrees
  fsd_eff_lam = 1.875*pi_over_180
  fsd_eff_phi = 1.25*pi_over_180
ELSE
  ! USE real grid-length in FSD parametrization
  fsd_eff_lam = delta_lambda
  fsd_eff_phi = delta_phi
END IF

IF ((l_rad_step_prog) .OR. (l_rad_step_diag)                                   &
    .OR. (l_fsd_generator)) THEN

  IF (i_fsd /= ip_fsd_constant) THEN
    ALLOCATE(f_arr(3, tdims%i_start:tdims%i_end,                               &
                      tdims%j_start:tdims%j_end, 1:tdims%k_end))
  END IF

  ! Calculate convection indicator for those parametrizations that require one.
  IF ((i_fsd == ip_fsd_regime)                                                 &
    .OR. (i_fsd == ip_fsd_regime_no_sh)                                        &
    .OR. (i_fsd == ip_fsd_regime_smooth)                                       &
    .OR. (i_fsd == ip_fsd_regime_smooth_no_sh)) THEN

    IF ((i_fsd == ip_fsd_regime_smooth)                                        &
      .OR. (i_fsd == ip_fsd_regime_smooth_no_sh)) THEN
      ALLOCATE(conv_ind_temp(tdims%i_start-1:tdims%i_end+1,                    &
                      tdims%j_start-1:tdims%j_end+1, 1:tdims%k_end))
      conv_ind_temp = 0.0
      ALLOCATE(conv_ind(tdims%i_start-2:tdims%i_end+2,                         &
                        tdims%j_start-2:tdims%j_end+2,                         &
                        1:tdims%k_end))
      conv_ind = 0.0

    ELSE
      ALLOCATE(conv_ind(tdims%i_start-1:tdims%i_end+1,                         &
                        tdims%j_start-1:tdims%j_end+1,                         &
                        1:tdims%k_end))
      conv_ind = 0.0

    END IF

    IF (i_fsd == ip_fsd_regime) THEN
      DO k=1, n_cca_levels
        DO j=1, rows
          DO i=1, row_length
            IF (cca_dp(i,j,k) > 0.0) conv_ind(i,j,k) = 1.0
            IF (cca_md(i,j,k) > 0.0) conv_ind(i,j,k) = 1.0
            IF (cca_sh(i,j,k) > 0.0) conv_ind(i,j,k) = 1.0
          END DO
        END DO
      END DO

    ELSE IF (i_fsd == ip_fsd_regime_smooth) THEN
      DO k=1, n_cca_levels
        DO j=1, rows
          DO i=1, row_length
            IF (cca_dp(i,j,k) > 0.0) conv_ind_temp(i,j,k) = 1.0
            IF (cca_md(i,j,k) > 0.0) conv_ind_temp(i,j,k) = 1.0
            IF (cca_sh(i,j,k) > 0.0) conv_ind_temp(i,j,k) = 1.0
          END DO
        END DO
      END DO

      CALL swap_bounds(conv_ind_temp, row_length, rows, n_cca_levels,          &
                      offx, offy, fld_type_p, swap_field_is_scalar)

      DO k=1, n_cca_levels
        DO j=0, rows+1
          DO i=0, row_length+1
            IF (conv_ind_temp(i,j,k) > 0.0) THEN
              conv_ind(i-1:i+1,j-1:j+1,k) = 1.0
            END IF
          END DO
        END DO
      END DO

    ELSE IF (i_fsd == ip_fsd_regime_no_sh) THEN
      DO k=1, n_cca_levels
        DO j=1, rows
          DO i=1, row_length
            IF (cca_dp(i,j,k) > 0.0) conv_ind(i,j,k) = 1.0
            IF (cca_md(i,j,k) > 0.0) conv_ind(i,j,k) = 1.0
          END DO
        END DO
      END DO

    ELSE IF (i_fsd == ip_fsd_regime_smooth_no_sh) THEN
      DO k=1, n_cca_levels
        DO j=1, rows
          DO i=1, row_length
            IF (cca_dp(i,j,k) > 0.0) conv_ind_temp(i,j,k) = 1.0
            IF (cca_md(i,j,k) > 0.0) conv_ind_temp(i,j,k) = 1.0
          END DO
        END DO
      END DO

      CALL swap_bounds(conv_ind_temp, row_length, rows, n_cca_levels,          &
                      offx, offy, fld_type_p, swap_field_is_scalar)

      DO k=1, n_cca_levels
        DO j=0, rows+1
          DO i=0, row_length+1
            IF (conv_ind_temp(i,j,k) > 0.0) THEN
              conv_ind(i-1:i+1,j-1:j+1,k) = 1.0
            END IF
          END DO
        END DO
      END DO

    END IF

  END IF ! convection indicator required

  IF (i_fsd == ip_fsd_param) THEN
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          conv_thick_part = 1.17 * ((r_theta_levels(i, j, k)                   &
            -r_theta_levels(i, j, k-1)) * 0.001) ** 0.11
          f_arr(1,i,j,k) = 0.12*conv_thick_part
          f_arr(2,i,j,k) = 0.23*conv_thick_part
          f_arr(3,i,j,k) = 0.05*conv_thick_part
        END DO
      END DO
    END DO
    f_cons(1) = 0.016
    f_cons(2) = 2.76
    f_cons(3) = -0.09

  ELSE IF ((i_fsd == ip_fsd_regime)                                            &
    .OR. (i_fsd == ip_fsd_regime_no_sh)) THEN
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          x_in_km = 0.001*SQRT ( r_theta_levels(i, j, k) * fsd_eff_lam         &
                                * r_theta_levels(i, j, k) * fsd_eff_phi        &
                                * xx_cos_theta_latitude(i,j)  )
          IF (conv_ind(i,j,k) > 0) THEN
            conv_thick_part = 2.81 * (x_in_km**(-0.12))                        &
              *((r_theta_levels(i, j, k)-r_theta_levels(i, j, k-1))            &
              * 0.001) ** 0.07
          ELSE
            conv_thick_part = 1.14 * (x_in_km**0.002)                          &
              *((r_theta_levels(i, j, k)-r_theta_levels(i, j, k-1))            &
              * 0.001) ** 0.12
          END IF
          f_arr(1,i,j,k) = 0.12*conv_thick_part
          f_arr(2,i,j,k) = 0.23*conv_thick_part
          f_arr(3,i,j,k) = 0.05*conv_thick_part
        END DO
      END DO
    END DO
    f_cons(1) = 0.016
    f_cons(2) = 2.76
    f_cons(3) = -0.09

  ELSE IF ((i_fsd == ip_fsd_regime_smooth)                                     &
    .OR. (i_fsd == ip_fsd_regime_smooth_no_sh)) THEN
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          x_in_km = 0.001*SQRT ( r_theta_levels(i, j, k) * fsd_eff_lam         &
                                * r_theta_levels(i, j, k) * fsd_eff_phi        &
                                * xx_cos_theta_latitude(i,j)  )
          IF (conv_ind(i,j,k) > 0) THEN
            conv_thick_part = 2.73 * (x_in_km**(-0.12))                        &
              *((r_theta_levels(i, j, k)-r_theta_levels(i, j, k-1))            &
              * 0.001) ** 0.08
          ELSE
            conv_thick_part = 1.12 * (x_in_km**0.002)                          &
              *((r_theta_levels(i, j, k)-r_theta_levels(i, j, k-1))            &
              * 0.001) ** 0.13
          END IF
          f_arr(1,i,j,k) = 0.12*conv_thick_part
          f_arr(2,i,j,k) = 0.23*conv_thick_part
          f_arr(3,i,j,k) = 0.05*conv_thick_part
        END DO
      END DO
    END DO
    f_cons(1) = 0.016
    f_cons(2) = 2.76
    f_cons(3) = -0.09

  ELSE IF (i_fsd == ip_fsd_boutle) THEN
    f_arr(1,:,:,:) = 0.11
    f_arr(2,:,:,:) = 0.45
    f_arr(3,:,:,:) = 0.25
    f_cons(1) = 0.06
    f_cons(2) = 1.50
    f_cons(3) = -0.17

  END IF

  IF (ALLOCATED(conv_ind)) DEALLOCATE(conv_ind)
  IF (ALLOCATED(conv_ind_temp)) DEALLOCATE(conv_ind_temp)

END IF ! Rad timestep or using param in microphysics.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_fsd_parameters

END MODULE set_fsd_parameters_mod
