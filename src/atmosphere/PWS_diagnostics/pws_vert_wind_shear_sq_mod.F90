! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate vertical wind shear squared on pressure levels
!
! Description:
!  Routine to calculate vertical wind shear squared (S_sq) on standard
! (pressure) levels using model level u and v fields.
!
! Method:
!  1) Calculate du/dz and dv/dz across model rho levels at u and v gridpoints.
!  2) Interpolate du/dz and dv/dz back to pressure gridpoint.
!  3) Calculate wind shear squared at this point
!  4) Determine height of midpoint between rho levels at pressure gridpoints
!     (this is the actual position where wind shear sq has been calculated)
!  5) For each standard (pressure) level, for each gridpoint, determine the
!     upper and lower wind shear squared values and their heights.
!  6) Calculate the wind shear squared on the standard level by interpolating
!     between these upper and lower values.
!
!   Programming standard; Unified Model Documentation Paper No. 3
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS_diagnostics

MODULE pws_vert_wind_shear_sq_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PWS_VERT_WIND_SHEAR_SQ_MOD'

CONTAINS

SUBROUTINE pws_vert_wind_shear_sq(num_press_lvls, press_lvls, u, v,      &
                              p_theta_levels, wind_shear_squared_plevs, mask)


USE atm_fields_bounds_mod
USE level_heights_mod,     ONLY: r_theta_levels
USE missing_data_mod,      ONLY: rmdi, imdi
USE nlsizes_namelist_mod,  ONLY: row_length, rows, model_levels
USE u_to_p_mod,            ONLY: u_to_p
USE v_to_p_mod,            ONLY: v_to_p
USE field_types,           ONLY: fld_type_p, fld_type_u, fld_type_v
USE halo_exchange,         ONLY: swap_bounds
USE mpp_conf_mod,          ONLY: swap_field_is_scalar
USE um_parvars,            ONLY: at_extremity, offx, offy
USE yomhook,               ONLY: lhook, dr_hook
USE parkind1,              ONLY: jprb, jpim

IMPLICIT NONE

! Inputs

! Pressure levels to derive BV-FREQ on
INTEGER, INTENT(IN) :: num_press_lvls
REAL, INTENT(IN) :: press_lvls(num_press_lvls)

REAL, INTENT(IN)  :: u(udims_s%i_start:udims_s%i_end,                    &
                       udims_s%j_start:udims_s%j_end,                    &
                       udims_s%k_start:udims_s%k_end)

REAL, INTENT(IN)  :: v(vdims_s%i_start:vdims_s%i_end,                    &
                       vdims_s%j_start:vdims_s%j_end,                    &
                       vdims_s%k_start:vdims_s%k_end)



REAL, INTENT(IN)  :: p_theta_levels(tdims_s%i_start:tdims_s%i_end,       &
                                    tdims_s%j_start:tdims_s%j_end,       &
                                    tdims_s%k_start:tdims_s%k_end)

! Outputs

! Wind shear squared on pressure levels
REAL, INTENT(OUT) :: wind_shear_squared_plevs(row_length, rows, num_press_lvls)

! Mask containing invalid values if FALSE
LOGICAL, INTENT(OUT) :: mask(row_length, rows, num_press_lvls)


! Local variables

! Rho level du/dz
REAL :: dudz(udims_s%i_start:udims_s%i_end,                              &
                       udims_s%j_start:udims_s%j_end,                    &
                       udims_s%k_start:udims_s%k_end)

! Rho level dv/dz
REAL :: dvdz(vdims_s%i_start:vdims_s%i_end,                              &
                       vdims_s%j_start:vdims_s%j_end,                    &
                       vdims_s%k_start:vdims_s%k_end)

! Rho level du/dz on pressure points
REAL :: dudz_Ppoints(row_length, rows, model_levels)

! Rho level dv/dz on pressure points
REAL :: dvdz_Ppoints(row_length, rows, model_levels)

! Wind shear squared on model levels
REAL :: wind_shear_squared(row_length, rows, model_levels)



! Interpolation:
! Model level corresponding to pressure level we're interested in
INTEGER :: level
! Fraction of the way we are through the layer
REAL :: alpha


INTEGER :: i,j,k,l ! Loop counters

CHARACTER (LEN=*), PARAMETER :: RoutineName='PWS_VERT_WIND_SHEAR_SQ'

! dr_hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

mask = .TRUE.


dudz = rmdi
dudz_Ppoints = rmdi

DO k = 1, model_levels-1
  DO j = udims%j_start,udims%j_end
    DO i = udims%i_start,udims%i_end
      dudz(i,j,k) = (u(i,j,k+1) - u(i,j,k)) /                            &
             (r_theta_levels(i,j,k+1) - r_theta_levels(i,j,k))
    END DO
  END DO
END DO


CALL swap_bounds(dudz, udims%i_len, udims%j_len, model_levels, offx,     &
                 offy, fld_type_u,swap_field_is_scalar)

CALL u_to_p(dudz(:,:,:), udims_s%i_start, udims_s%i_end,                 &
                         udims_s%j_start, udims_s%j_end,                 &
                         tdims%i_start,   tdims%i_end,                   &
                         tdims%j_start,   tdims%j_end,                   &
                         model_levels, at_extremity, dudz_Ppoints(:,:,:))


dvdz = rmdi
dvdz_Ppoints = rmdi

DO k = 1, model_levels-1
  DO j = vdims%j_start,vdims%j_end
    DO i = vdims%i_start,vdims%i_end
      dvdz(i,j,k) = (v(i,j,k+1) - v(i,j,k)) / &
             (r_theta_levels(i,j,k+1) - r_theta_levels(i,j,k))
    END DO
  END DO
END DO


CALL swap_bounds(dvdz, vdims%i_len, vdims%j_len, model_levels, offx,     &
                 offy, fld_type_v, swap_field_is_scalar)

CALL v_to_p(dvdz(:,:,:), vdims_s%i_start, vdims_s%i_end,                 &
                         vdims_s%j_start, vdims_s%j_end,                 &
                         tdims%i_start,   tdims%i_end,                   &
                         tdims%j_start,   tdims%j_end,                   &
                         model_levels, at_extremity, dvdz_Ppoints(:,:,:))


! Calculate wind shear squared on model levels
wind_shear_squared = dudz_Ppoints**2 + dvdz_Ppoints **2


! Interpolate wind shear squared onto requested pressure levels
DO l = 1, num_press_lvls

  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      level = imdi
vert: DO k = tdims%k_start, tdims%k_end - 1
        IF ((press_lvls(l) >  p_theta_levels(i,j,k)  .AND.               &
             press_lvls(l) <= p_theta_levels(i,j,k+1) ) .OR.             &
            (press_lvls(l) <  p_theta_levels(i,j,k)  .AND.               &
             press_lvls(l) >= p_theta_levels(i,j,k+1) ) ) THEN
              level = k
              alpha = (press_lvls(l) - p_theta_levels(i,j,level))        &
                  /  (p_theta_levels(i,j,level+1) - p_theta_levels(i,j,level))
              EXIT vert
        END IF
      END DO vert
      IF (level == imdi) THEN
        mask(i,j,l) = .FALSE.
      ELSE
        wind_shear_squared_plevs(i,j,l) =                                & 
                wind_shear_squared(i,j,level) + alpha *                  &
           (wind_shear_squared(i,j,level+1)  - wind_shear_squared(i,j,level))

        ! Avoid a potential division by zero later on
        IF (wind_shear_squared_plevs(i,j,l) == 0.0) THEN
          mask(i,j,l) = .FALSE.
        END IF

      END IF
      
    END DO
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE pws_vert_wind_shear_sq

END MODULE pws_vert_wind_shear_sq_mod
