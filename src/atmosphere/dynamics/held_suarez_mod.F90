! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE held_suarez_mod
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='HELD_SUAREZ_MOD'

CONTAINS

SUBROUTINE eg_held_suarez(                                        &

! in data fields.
  u, v, theta,  exner_theta_levels, exner,                        &
! in/out
  theta_star,  r_u, r_v,                                          &
  L_HeldSuarez, L_HeldSuarez1_drag,                               &
! error information
  error_code  )

USE planet_constants_mod,      ONLY: recip_kappa
USE conversions_mod,           ONLY: rsec_per_day
USE water_constants_mod
USE parkind1,                  ONLY: jpim, jprb       !DrHook
USE yomhook,                   ONLY: lhook, dr_hook   !DrHook
USE level_heights_mod,         ONLY: r_theta_levels
USE timestep_mod,              ONLY: timestep
USE atm_fields_bounds_mod
USE horiz_grid_mod

IMPLICIT NONE
!
! Description: Initialises star state for ENDGame when there is no
!              physics. Optionally applies Held-Suarez test of
!              dynamical core (BAMS 75, 1825-1830).
!
!
! Method:
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics
!
! Code description:
!   Language: Fortran 90.
!   This code is written to programming standards in UMDP3.

! Subroutine arguments


! Data arrays
REAL, INTENT (INOUT) ::                                               &
  theta(             tdims_s%i_start:tdims_s%i_end,                   &
                     tdims_s%j_start:tdims_s%j_end,                   &
                     tdims_s%k_start:tdims_s%k_end )

REAL, INTENT (INOUT) ::                                               &
  theta_star(        tdims%i_start:tdims%i_end,                       &
                     tdims%j_start:tdims%j_end,                       &
                     tdims%k_start:tdims%k_end),                      &
  r_u(               udims_s%i_start:udims_s%i_end,                   &
                     udims_s%j_start:udims_s%j_end,                   &
                     udims_s%k_start:udims_s%k_end),                  &
  r_v(               vdims_s%i_start:vdims_s%i_end,                   &
                     vdims_s%j_start:vdims_s%j_end,                   &
                     vdims_s%k_start:vdims_s%k_end),                  &
  u(                 udims_s%i_start:udims_s%i_end,                   &
                     udims_s%j_start:udims_s%j_end,                   &
                     udims_s%k_start:udims_s%k_end),                  &
  v(                 vdims_s%i_start:vdims_s%i_end,                   &
                     vdims_s%j_start:vdims_s%j_end ,                  &
                     vdims_s%k_start:vdims_s%k_end),                  &
  exner_theta_levels(tdims_s%i_start:tdims_s%i_end,                   &
                     tdims_s%j_start:tdims_s%j_end,                   &
                     tdims_s%k_start:tdims_s%k_end),                  &
  exner(             pdims_s%i_start:pdims_s%i_end,                   &
                     pdims_s%j_start:pdims_s%j_end,                   &
                     pdims_s%k_start:pdims_s%k_end+1)

INTEGER :: error_code

LOGICAL :: L_HeldSuarez, L_HeldSuarez1_drag

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_HELD_SUAREZ'

! Variables for Held Suarez test case:
INTEGER :: i, j,k
REAL    :: temp1, temp2, SuHe_sigma_cutoff, friction_level,           &
           base_frictional_timescale, SuHe_pole_equ_deltaT,           &
           SuHe_static_stab, SuHe_level_weight,                       &
           SuHe_newtonian_timescale_ka,                               &
           SuHe_newtonian_timescale_ks, newtonian_timescale, sigma

REAL    ::                                                            &
           theta_eq(tdims_s%i_start:tdims_s%i_end,                    &
                    tdims_s%j_start:tdims_s%j_end,                    &
                    tdims_s%k_start:tdims_s%k_end)


! 1.0 Start of subroutine code: perform the calculation.
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (L_HeldSuarez) THEN
  ! ----------------------------------------------------------------------
  ! Section 1. Initialise parameters
  ! ----------------------------------------------------------------------
  !
  !     k_f = 1 1/day in H&S
  !
  !     here: 1/(60*60*24) 1/s = 1./86400. 1/s = 1./rsec_per_day
  base_frictional_timescale = 1.0/rsec_per_day

  !     sigma_b
  SuHe_sigma_cutoff = 0.7

  SuHe_pole_equ_deltaT = 60.0
  SuHe_static_stab     = 10.0

  ! 1/40 day^-1 = 1/(40*86400) 1/s = 1/(N*rsec_per_day)
  SuHe_newtonian_timescale_ka = 1.0/(40.0*rsec_per_day)
  SuHe_newtonian_timescale_ks = 1.0/(4.0*rsec_per_day)

  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

        ! calculate equilibrium potential temperature (theta_eq)
        temp1 = 200.0 / exner_theta_levels(i,j,k)
        temp2 = 315.0                                             &
                - SuHe_pole_equ_deltaT * Snxi2_p(j) * Snxi2_p(j)  &
                - SuHe_static_stab * recip_kappa                  &
                  * LOG(exner_theta_levels(i,j,k))                &
                  * Csxi2_p(j) * Csxi2_p(j)

        theta_eq(i,j,k) = MAX(temp1, temp2)

      END DO
    END DO
  END DO

  L_HeldSuarez_drag_if:  IF (L_HeldSuarez1_drag) THEN
    !   add to increment field
    DO k = udims%k_start, udims%k_end
      DO j = udims%j_start, udims%j_end
        DO i = udims%i_start, udims%i_end

          sigma= 0.5*( (exner(i,j,k) /                            &
                      exner_theta_levels(i,j,0))**recip_kappa +   &
                     (exner(i+1,j,k) /                            &
                      exner_theta_levels(i+1,j,0))**recip_kappa )

          temp1 = ( sigma - SuHe_sigma_cutoff)  /                 &
                ( 1.0   - SuHe_sigma_cutoff)

          SuHe_level_weight = MAX(0.0, temp1)
          friction_level = base_frictional_timescale              &
                         * SuHe_level_weight

          r_u(i,j,k) = r_u(i,j,k) - timestep * friction_level     &
                                  * u(i,j,k)
        END DO
      END DO
    END DO

    !   add on to increment field
    DO k = vdims%k_start, vdims%k_end
      DO j = vdims%j_start, vdims%j_end
        DO i = vdims%i_start, vdims%i_end

          sigma= 0.5*( (exner(i,j,k) /                            &
                      exner_theta_levels(i,j,0))**recip_kappa +   &
                     (exner(i,j+1,k) /                            &
                      exner_theta_levels(i,j+1,0))**recip_kappa )

          temp1 = ( sigma - SuHe_sigma_cutoff)  /                 &
                ( 1.0   - SuHe_sigma_cutoff)

          SuHe_level_weight = MAX(0.0, temp1)
          friction_level = base_frictional_timescale              &
                         * SuHe_level_weight

          r_v(i,j,k) = r_v(i,j,k) - timestep * friction_level     &
                                   * v(i,j,k)
        END DO
      END DO
    END DO

  END IF L_HeldSuarez_drag_if

  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

        sigma=(exner_theta_levels(i,j,k) /                        &
               exner_theta_levels(i,j,0))**recip_kappa

        temp1 = ( sigma - SuHe_sigma_cutoff)  /                   &
                ( 1.0   - SuHe_sigma_cutoff)

        SuHe_level_weight = MAX(0.0, temp1)
        ! calculate relaxation term.
        newtonian_timescale = SuHe_newtonian_timescale_ka         &
                                + ( SuHe_newtonian_timescale_ks - &
                                    SuHe_newtonian_timescale_ka ) &
                                * Csxi2_p(j) ** 4                 &
                                * SuHe_level_weight

        theta_star(i,j,k) = theta_star(i,j,k)                     &
                            - timestep * newtonian_timescale *    &
                            (theta(i,j,k) - theta_eq(i,j,k))
      END DO
    END DO
  END DO

END IF ! L_HeldSuarez

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_held_suarez


END MODULE

