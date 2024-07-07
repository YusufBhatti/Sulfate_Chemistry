! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE eg_idl_theta_forcing_mod

  ! Modules for supported setups
  USE held_suarez_forcing_mod
  USE planet_forcing_mod

  USE umPrintMgr
  USE ereport_mod,               ONLY: ereport

  IMPLICIT NONE

  ! Description: 
  !  Module containing subroutines to assign, and employ
  !  the selected temperature forcing (actually force
  !  potential temperature)
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
  ! This file belongs in section: IDEALISED
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !   This code is written to UMDP3 programming standards.

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_IDL_THETA_FORCING_MOD'

CONTAINS

  SUBROUTINE eg_idl_theta_forcing(                                           &
                                ! in data fields.
       theta, exner_theta_levels, exner,                                     &
                                ! in/out
       theta_star,                                                           &
                                ! Profile switches, and anciliary info
       tforce_number, trelax_number,                                         &
       t_surface,                                                            &
                                ! error information
       error_code  )

    ! Description: This forces the potential temperature towards
    !              a selected (tforce_number) temperature profile
    !              using a specified newtonian timescale scheme 
    !              (trelax_number)

    USE parkind1,                  ONLY: jpim, jprb       !DrHook
    USE yomhook,                   ONLY: lhook, dr_hook   !DrHook
    USE timestep_mod,              ONLY: timestep
    USE atm_fields_bounds_mod
    USE horiz_grid_mod
    USE tforce_mod,                ONLY: tf_none
    USE trelax_mod,                ONLY: tr_none
    USE model_domain_mod

    IMPLICIT NONE

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_IDL_THETA_FORCING'

    ! Subroutine arguments

    ! Data arrays
    REAL, INTENT (INOUT) ::                                                  &
         theta(             tdims_s%i_start:tdims_s%i_end,                   &
         tdims_s%j_start:tdims_s%j_end,                                      &
         tdims_s%k_start:tdims_s%k_end )

    REAL, INTENT (INOUT) ::                                                  &
         theta_star(tdims%i_start:tdims%i_end,                               &
                    tdims%j_start:tdims%j_end,                               &
                    tdims%k_start:tdims%k_end),                              &
         exner_theta_levels(tdims_s%i_start:tdims_s%i_end,                   &
         tdims_s%j_start:tdims_s%j_end,                                      &
         tdims_s%k_start:tdims_s%k_end),                                     &
         exner(             pdims_s%i_start:pdims_s%i_end,                   &
         pdims_s%j_start:pdims_s%j_end,                                      &
         pdims_s%k_start:pdims_s%k_end+1)

    INTEGER :: tforce_number, trelax_number
    REAL    :: t_surface

    INTEGER :: error_code

    ! Local Variables for the temperature profile settings
    ! and the relaxation
    INTEGER :: i, j,k

    REAL    ::                                                               &
         theta_eq(tdims_s%i_start:tdims_s%i_end,                             &
         tdims_s%j_start:tdims_s%j_end,                                      &
         tdims_s%k_start:tdims_s%k_end),                                     &
         newtonian_timescale(tdims_s%i_start:tdims_s%i_end,                  &
         tdims_s%j_start:tdims_s%j_end,                                      &
         tdims_s%k_start:tdims_s%k_end)

    ! 1.0 Start of subroutine code: perform the calculation.
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
    ! ----------------------------------------------------------------------
    ! Section 1. Create the Equlibrium Theta and Newtonian Timescale
    ! ----------------------------------------------------------------------
    ! Initialise the forcing variables
    newtonian_timescale=0.0
    theta_eq=theta
    ! If we are not forcing do nothing
    IF (tforce_number == tf_none .OR. trelax_number == tr_none) THEN
       ! No forcing 
       theta_eq=theta
       newtonian_timescale=0.0
    ELSE
       ! Set the equilibrium theta and newtonian timescale
       DO k = tdims%k_start, tdims%k_end 
          DO j = tdims%j_start, tdims%j_end 
             DO i = tdims%i_start, tdims%i_end 
                theta_eq(i,j,k)=equilibrium_theta(tforce_number, i, j, k,    &
                     exner_theta_levels, t_surface)
                newtonian_timescale(i,j,k)=recip_relaxation_tscale(          &
                     trelax_number, i, j, k, exner_theta_levels, theta)
                ! ----------------------------------------------------------
                ! Section 1.2 Perform Relaxation
                ! ----------------------------------------------------------
                ! Check we do not over-correct
                IF (timestep * newtonian_timescale(i,j,k) > 1.0)             &
                     newtonian_timescale(i,j,k)=1.0/timestep
                ! Apply parameterised heating/cooling
                theta_star(i,j,k) = theta_star(i,j,k)                        &
                     - timestep * newtonian_timescale(i,j,k) *               &
                     (theta(i,j,k) - theta_eq(i,j,k))
             END DO
          END DO
       END DO
    END IF

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END SUBROUTINE eg_idl_theta_forcing

  ! ----------------------------------------------------------------------
  ! FUNCTIONS: 1. Calculate the equilibrium temperature
  ! ----------------------------------------------------------------------

  REAL FUNCTION equilibrium_theta(tforce_number, i, j, k,                    &
       exner_theta_levels, t_surface)
    ! This function returns the equilibrium
    ! potential temperature from the
    ! profile=tforce_number, at position
    ! i,j,k
    USE tforce_mod

    USE parkind1,                  ONLY: jpim, jprb       !DrHook
    USE yomhook,                   ONLY: lhook, dr_hook   !DrHook

    IMPLICIT NONE

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='EQUILIBRIUM_THETA'

    ! INPUT
    INTEGER :: tforce_number,err_code 
    INTEGER :: i,j,k
    REAL    :: t_surface
    REAL    ::   exner_theta_levels(tdims_s%i_start:tdims_s%i_end,           &
         tdims_s%j_start:tdims_s%j_end,                                      &
         tdims_s%k_start:tdims_s%k_end) 
    ! OUTPUT
    REAL :: theta_out

    ! 1.0 Start of function code: perform the calculation.
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ! Set the required potential temperature
    IF (tforce_number == tf_isothermal) THEN
       ! Isothermal Sphere
       theta_out=t_surface/exner_theta_levels(i,j,k)
    ELSE IF (tforce_number == tf_HeldSuarez) THEN
       theta_out=held_suarez_theta(i, j, k,                                  &
            exner_theta_levels)
    ELSE IF (l_planet_forcing(tforce_number)) THEN
       theta_out=planet_forcing_theta(tforce_number, i, j, k,                &
            exner_theta_levels)
    ELSE
       WRITE(umMessage,'(2A)') '-T-P profile selection not supported:',      &
            tforce_number
       CALL umPrint(umMessage,src='eg_idl_theta_forcing_mod.F90')
       err_code = 1
       CALL ereport("eg_idl_theta_forcing_mod", err_code,                    &
            "T-P forcing profile not supported" )
    END IF
    equilibrium_theta=theta_out

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END FUNCTION equilibrium_theta

  ! ----------------------------------------------------------------------
  ! FUNCTIONS: 2. Derive relaxation timescale
  ! ----------------------------------------------------------------------

  REAL FUNCTION recip_relaxation_tscale(trelax_number, i, j, k,              &
       exner_theta_levels, theta)
    ! This function returns the relaxation timescale
    ! for the trelax_number, at position
    ! i,j,k
    USE trelax_mod

    USE parkind1,                  ONLY: jpim, jprb       !DrHook
    USE yomhook,                   ONLY: lhook, dr_hook   !DrHook

    IMPLICIT NONE

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='RECIP_RELAXATION_TSCALE'

    ! INPUT
    INTEGER :: trelax_number,err_code 
    INTEGER :: i,j,k
    REAL    :: exner_theta_levels(tdims_s%i_start:tdims_s%i_end,             &
         tdims_s%j_start:tdims_s%j_end,                                      &
         tdims_s%k_start:tdims_s%k_end) 
    REAL    :: theta(tdims_s%i_start:tdims_s%i_end,                          &
         tdims_s%j_start:tdims_s%j_end,                                      &
         tdims_s%k_start:tdims_s%k_end )
    ! OUTPUT
    REAL :: recip_tscale_out

    ! 1.0 Start of function code: perform the calculation.
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    IF (trelax_number == tr_HeldSuarez) THEN
       recip_tscale_out=held_suarez_recip_newt(i, j, k,                      &
            exner_theta_levels)
    ELSE IF (l_planet_relax(trelax_number)) THEN
       recip_tscale_out=planet_recip_newt(trelax_number, i, j, k,            &
            exner_theta_levels)
    ELSE
       WRITE(umMessage,'(2A)') '-Radiative timescale not supported:',        &
            trelax_number
       CALL umPrint(umMessage,src='eg_idl_theta_forcing_mod.F90')
       err_code = 1
       CALL ereport("eg_idl_theta_forcing_mod", err_code,                    &
            "Radiative Timescale not supported" )
    END IF
    recip_relaxation_tscale=recip_tscale_out

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END FUNCTION recip_relaxation_tscale

END MODULE eg_idl_theta_forcing_mod
