! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE held_suarez_forcing_mod

  USE planet_constants_mod
  USE conversions_mod,           ONLY: rsec_per_day
  USE atm_fields_bounds_mod
  USE horiz_grid_mod
  USE timestep_mod,              ONLY: timestep
  USE model_domain_mod

  IMPLICIT NONE
  
  ! Description: 
  !   Module containing the subroutines/functions required for
  !   an idealised Held-Suarez test. Including the derivation
  !   of the required equilibrium potential temperature,
  !   relaxation timescale and friction.
  !
  ! Method:
  !   The values are derived as described in Held & Suarez (1995)
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
  ! This file belongs in section: IDEALISED
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !   This code is written to UMDP3 programming standards.
  
  CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='HELD_SUAREZ_FORCING_MOD'

CONTAINS
  
  REAL FUNCTION held_suarez_theta(i, j, k, exner_theta_levels)

    USE parkind1,                  ONLY: jpim, jprb       !DrHook
    USE yomhook,                   ONLY: lhook, dr_hook   !DrHook

    IMPLICIT NONE

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='HELD_SUAREZ_THETA'

    ! INPUT
    INTEGER :: i,j,k
    REAL    ::   exner_theta_levels(tdims_s%i_start:tdims_s%i_end,       &
         tdims_s%j_start:tdims_s%j_end,                                  &
         tdims_s%k_start:tdims_s%k_end) 
    ! LOCAL VARIABLES
    REAL :: T_min, T_surf  ! Minimum/Stratospheric & surface temperature
    REAL :: DT_eq_pole     ! Equator-Pole Temp diff
    REAL :: Static_Stab    ! Static Stability temperature
    REAL :: theta1, theta2 ! Potential temperatures
    REAL :: theta_out
    ! Function to calculate the Held Suarez prescribed
    ! potential temperature (Held & Suarez, 1995)

    ! 1.0 Start of function code: perform the calculation.
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ! Held-Suarez Temperature profile
    Static_Stab=10.0
    DT_eq_pole=60.0
    T_min=200.0
    T_surf=315.0
    theta1 = T_min / exner_theta_levels(i,j,k)
    theta2 = T_surf                                                      &
         - DT_eq_pole * Snxi2_p(j) * Snxi2_p(j)                          &
         - Static_Stab * recip_kappa                                     &
         * LOG(exner_theta_levels(i,j,k))                                &
         * Csxi2_p(j) * Csxi2_p(j)
    theta_out = MAX(theta1, theta2)
    held_suarez_theta=theta_out

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END FUNCTION held_suarez_theta

  REAL FUNCTION held_suarez_recip_newt(i, j, k, exner_theta_levels)
    USE parkind1,                  ONLY: jpim, jprb       !DrHook
    USE yomhook,                   ONLY: lhook, dr_hook   !DrHook

    IMPLICIT NONE

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='HELD_SUAREZ_RECIP_NEWT'

    ! INPUT
    INTEGER :: i,j,k
    REAL    ::   exner_theta_levels(tdims_s%i_start:tdims_s%i_end,      &
         tdims_s%j_start:tdims_s%j_end,                                 &
         tdims_s%k_start:tdims_s%k_end) 

    ! LOCAL VARIABLES
    ! The level weighting for Held-Suarez Relaxation scheme
    REAL :: sigma          ! The relative pressure of the level
    REAL :: SuHe_level_weight, temp_weight
    REAL :: SuHe_sigma_cutoff
    REAL :: SuHe_newtonian_timescale_ka, SuHe_newtonian_timescale_ks
    REAL :: recip_tscale_out

    ! Function to calculate the Held Suarez prescribed
    ! newtonian relaxation timescale (reciprocal) 
    ! (Held & Suarez, 1995)

    ! 1.0 Start of function code: perform the calculation.
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ! Held-Suarez (Held & Suarez, 1995)
    !     sigma_b
    SuHe_sigma_cutoff = 0.7

    ! 1/40 day**-1 = 1/(40*86400) 1/s = 1/(N*rsec_per_day)
    SuHe_newtonian_timescale_ka = 1.0/(40.0*rsec_per_day)
    SuHe_newtonian_timescale_ks = 1.0/(4.0*rsec_per_day)
    ! Derive sigma level
    sigma=(exner_theta_levels(i,j,k) /                                   &
         exner_theta_levels(i,j,0))**recip_kappa
    ! Create the level weighting
    temp_weight = ( sigma - SuHe_sigma_cutoff)  /                        &
         ( 1.0   - SuHe_sigma_cutoff)            
    SuHe_level_weight = MAX(0.0, temp_weight)
    ! Calculate the timescale
    recip_tscale_out = SuHe_newtonian_timescale_ka                       &
         + ( SuHe_newtonian_timescale_ks -                               &
         SuHe_newtonian_timescale_ka )                                   &
         * Csxi2_p(j) ** 4                                               &
         * SuHe_level_weight
    held_suarez_recip_newt=recip_tscale_out

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END FUNCTION held_suarez_recip_newt

  SUBROUTINE held_suarez_friction(                                       &
                                ! In data fields
       u,v, exner_theta_levels, exner,                                   &
                                ! In/out data fields
       r_u, r_v,                                                         &
       fric_num) 

    USE bl_option_mod, ONLY: fric_HeldSuarez_1, fric_HeldSuarez_2
    USE parkind1,                  ONLY: jpim, jprb       !DrHook
    USE yomhook,                   ONLY: lhook, dr_hook   !DrHook

    IMPLICIT NONE

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='HELD_SUAREZ_FRICTION'

    ! Data arrays
    REAL, INTENT (INOUT) ::                                              &
         r_u(               udims_s%i_start:udims_s%i_end,               &
         udims_s%j_start:udims_s%j_end,                                  &
         udims_s%k_start:udims_s%k_end),                                 &
         r_v(               vdims_s%i_start:vdims_s%i_end,               &
         vdims_s%j_start:vdims_s%j_end,                                  &
         vdims_s%k_start:vdims_s%k_end),                                 &
         u(                 udims_s%i_start:udims_s%i_end,               &
         udims_s%j_start:udims_s%j_end,                                  &
         udims_s%k_start:udims_s%k_end),                                 &
         v(                 vdims_s%i_start:vdims_s%i_end,               &
         vdims_s%j_start:vdims_s%j_end,                                  &
         vdims_s%k_start:vdims_s%k_end),                                 &
         exner_theta_levels(tdims_s%i_start:tdims_s%i_end,               &
         tdims_s%j_start:tdims_s%j_end,                                  &
         tdims_s%k_start:tdims_s%k_end),                                 &
         exner(             pdims_s%i_start:pdims_s%i_end,               &
         pdims_s%j_start:pdims_s%j_end,                                  &
         pdims_s%k_start:pdims_s%k_end+1)  
    INTEGER :: fric_num

    ! Local Variables for the friction application
    INTEGER :: i, j,k

    REAL :: SuHe_sigma_cutoff
    REAL :: base_frictional_timescale
    REAL :: sigma
    REAL :: temp_weight
    REAL :: SuHe_level_weight
    REAL :: friction_level

    ! 1.0 Start of subroutine code: perform the calculation.
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ! Set parameters
    !     k_f = 1 1/day in H&S
    !
    !     here: 1/(60*60*24) 1/s = 1./86400. 1/s = 1./rsec_per_day
    base_frictional_timescale = 1.0/rsec_per_day
    !     sigma_b
    SuHe_sigma_cutoff = 0.7
    !   add to increment u field
    DO k = udims%k_start, udims%k_end 
       DO j = udims%j_start, udims%j_end 
          DO i = udims%i_start, udims%i_end 
             sigma= 0.5*( (exner(i,j,k) /                                &
                  exner_theta_levels(i,j,0))**recip_kappa +              &
                  (exner(i+1,j,k) /                                      &
                  exner_theta_levels(i+1,j,0))**recip_kappa )

             temp_weight = ( sigma - SuHe_sigma_cutoff)  /               &
                  ( 1.0   - SuHe_sigma_cutoff)

             SuHe_level_weight = MAX(0.0, temp_weight)
             friction_level = base_frictional_timescale                  &
                  * SuHe_level_weight
             IF (fric_num == fric_HeldSuarez_1) THEN    
                r_u(i,j,k) = r_u(i,j,k) - timestep * friction_level      &
                     * u(i,j,k)
             ELSE IF (fric_num == fric_HeldSuarez_2) THEN
                r_u(i,j,k) = r_u(i,j,k) - timestep * friction_level      &
                     * (r_u(i,j,k)+u(i,j,k))
             END IF
          END DO
       END DO
    END DO


    !   add on to increment v field
    DO k = vdims%k_start, vdims%k_end 
       DO j = vdims%j_start, vdims%j_end 
          DO i = vdims%i_start, vdims%i_end             
             sigma= 0.5*( (exner(i,j,k) /                                &
                  exner_theta_levels(i,j,0))**recip_kappa +              &
                  (exner(i,j+1,k) /                                      &
                  exner_theta_levels(i,j+1,0))**recip_kappa )

             temp_weight = ( sigma - SuHe_sigma_cutoff)  /               &
                  ( 1.0   - SuHe_sigma_cutoff)

             SuHe_level_weight = MAX(0.0, temp_weight)
             friction_level = base_frictional_timescale                  &  
                  * SuHe_level_weight

             IF (fric_num == fric_HeldSuarez_1) THEN    
                r_v(i,j,k) = r_v(i,j,k) - timestep * friction_level      &
                     * v(i,j,k)
             ELSE IF (fric_num == fric_HeldSuarez_2) THEN
                r_v(i,j,k) = r_v(i,j,k) - timestep * friction_level      &
                     * (r_v(i,j,k)+v(i,j,k))
             END IF
          END DO
       END DO
    END DO

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END SUBROUTINE held_suarez_friction

END MODULE held_suarez_forcing_mod
