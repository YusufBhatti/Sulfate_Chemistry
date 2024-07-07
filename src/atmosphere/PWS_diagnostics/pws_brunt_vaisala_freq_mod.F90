! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate Brunt-Vaisala frequency squared on pressure levels
!
! Method:
!   1) Calculate d(potential temperature)/dz (here dtheta/dz) for each model
!      slab from input theta and height model level fields
!   2) Determine pressure (position) of the mid point between model (theta)
!      levels
!   3) Then for each standard (pressure) level and each gridpoint:
!           i) determine the potential temperature theta by finding the model
!              levels above and below the standard (pressure) level and
!              interpolating the values of potential temperature on these
!              levels.
!          ii) determine which model slabs are above and below it by comparing
!              mid point pressure values. Hence determine the dtheta/dz values
!              of the mid points above and below the standard pressure level
!              and interpolate these to find dtheta/dz on the standard
!              (pressure) level.
!         iii) calculate Brunt Vaisala frequency as
!              N_sq = (g/theta)*(dtheta/dz)
!   4) Determine the length scale as the vertical distance between the 
!      model level above and the model level below each standard level 
!      gridpoint.
!
!   Programming standard; Unified Model Documentation Paper No. 3
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS_diagnostics

MODULE pws_brunt_vaisala_freq_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PWS_BRUNT_VAISALA_FREQ_MOD'

CONTAINS

SUBROUTINE pws_brunt_vaisala_freq(num_press_lvls, press_lvls, theta,         &
                              p_theta_levels, bv_freq_p_levs, dz, mask)

USE atm_fields_bounds_mod
USE level_heights_mod,     ONLY: r_theta_levels
USE missing_data_mod,      ONLY: rmdi, imdi
USE planet_constants_mod,  ONLY: g
USE nlsizes_namelist_mod,  ONLY: row_length, rows, model_levels
USE yomhook,               ONLY: lhook, dr_hook
USE parkind1,              ONLY: jprb, jpim

IMPLICIT NONE

! Inputs

! Pressure levels to derive BV-FREQ on
INTEGER, INTENT(IN) :: num_press_lvls
REAL, INTENT(IN) :: press_lvls(num_press_lvls)

REAL, INTENT(IN)  :: theta(tdims_s%i_start:tdims_s%i_end,                    &
                           tdims_s%j_start:tdims_s%j_end,                    &
                           tdims_s%k_start:tdims_s%k_end)


REAL, INTENT(IN)  :: p_theta_levels(tdims_s%i_start:tdims_s%i_end,           &
                                    tdims_s%j_start:tdims_s%j_end,           &
                                    tdims_s%k_start:tdims_s%k_end)

! Outputs

! Length scale
REAL, INTENT(OUT) :: dz(row_length, rows, num_press_lvls)

! Brunt-Vaisala Frequency on pressure levels
REAL, INTENT(OUT) :: bv_freq_p_levs(row_length, rows, num_press_lvls)

! Mask containing invalid values if FALSE
LOGICAL, INTENT(OUT) :: mask(row_length, rows, num_press_lvls)


! Local variables

! Rho level d(theta)/dZ
REAL :: dthetadZ(row_length, rows, model_levels)

! Pressure level d(theta)/dZ
REAL :: dthetadZ_plevs(row_length, rows, num_press_lvls)

! Theta on pressure levels
REAL :: theta_plevs(row_length, rows, num_press_lvls)


! Interpolation:
! Model level corresponding to pressure level we're interested in
INTEGER :: level
! Fraction of the way we are through the layer
REAL :: alpha


INTEGER :: i,j,k,l ! Loop counters

CHARACTER (LEN=*), PARAMETER :: RoutineName='PWS_BRUNT_VAISALA_FREQ'

! dr_hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


theta_plevs = rmdi
dthetadz = rmdi
dthetadZ_plevs = rmdi
dz = rmdi
bv_freq_p_levs = rmdi
mask = .TRUE.

DO k = 1, model_levels - 1
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dthetadz(i,j,k) = (theta(i,j,k+1) - theta(i,j,k)) /                    &
             (r_theta_levels(i,j,k+1) - r_theta_levels(i,j,k))
    END DO
  END DO
END DO

! Interpolate theta and d(theta)/dz onto requested pressure levels, then
! calculate B-V frequency.
DO l = 1, num_press_lvls

  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
       level = imdi
vert:  DO k = 1, model_levels - 1
        IF ((press_lvls(l) >  p_theta_levels(i,j,k)     .AND.                &
             press_lvls(l) <= p_theta_levels(i,j,k+1) ) .OR.                 &
            (press_lvls(l) <  p_theta_levels(i,j,k)     .AND.                &
             press_lvls(l) >= p_theta_levels(i,j,k+1) ) ) THEN
              level = k
              dz(i,j,l) = r_theta_levels(i,j,k+1) - r_theta_levels(i,j,k)
              alpha = (press_lvls(l) - p_theta_levels(i,j,k)) /              &
                      (p_theta_levels(i,j,k+1) - p_theta_levels(i,j,k))
              EXIT vert
        END IF
      END DO vert

      IF (level == imdi) THEN
        mask(i,j,l) = .FALSE.
      ELSE 
        dthetadZ_plevs(i,j,l) =                                              & 
                dthetadZ(i,j,level) + alpha *                                &
                (dthetadZ(i,j,level+1)  - dthetadZ(i,j,level))

        theta_plevs(i,j,l) =                                                 & 
                theta(i,j,level) + alpha *                                   &
                (theta(i,j,level+1)  - theta(i,j,level))

        bv_freq_p_levs(i,j,l) = (g/theta_plevs(i,j,l)) *                     &
                                                      (dthetadZ_plevs(i,j,l))
      END IF
    END DO
  END DO
END DO


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE pws_brunt_vaisala_freq

END MODULE pws_brunt_vaisala_freq_mod
