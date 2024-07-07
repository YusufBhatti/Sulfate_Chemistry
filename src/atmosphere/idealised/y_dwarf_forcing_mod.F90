! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE y_dwarf_forcing_mod

  USE planet_constants_mod
  USE atm_fields_bounds_mod
  USE model_domain_mod
  USE horiz_grid_mod
  USE gravity_mod, ONLY: g_theta

  IMPLICIT NONE

  ! Description: 
  !   Module containing the subroutines/functions required for
  !   an idealised Y-Dwarf test. Including the derivation
  !   of the required equilibrium potential temperature,
  !   relaxation timescale and friction.
  !
  ! Method:
  !   T-P profile derived using ATMO (University of Exeter)
  !   Published in Tremblin et al (2015). Relaxation 
  !   timescale is the generic one from Showman & Guillot (2002)
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
  ! This file belongs in section: IDEALISED
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !   This code is written to UMDP3 programming standards.

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='Y_DWARF_FORCING_MOD'

CONTAINS

  REAL FUNCTION y_dwarf_theta(i, j, k, exner_theta_levels)
    USE parkind1,                  ONLY: jpim, jprb       !DrHook
    USE yomhook,                   ONLY: lhook, dr_hook   !DrHook

    IMPLICIT NONE
    
    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='Y_DWARF_THETA'
    
    ! INPUT
    INTEGER :: i,j,k
    REAL    ::   exner_theta_levels(tdims_s%i_start:tdims_s%i_end,       &
         tdims_s%j_start:tdims_s%j_end,                  &
         tdims_s%k_start:tdims_s%k_end) 
    
    REAL :: t_bd
    REAL :: theta_out

    ! Function to calculate the Y Dwarf potential temperature
    ! from Pascal Tremblin

    ! 1.0 Start of function code: perform the calculation.
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ! Calculate the temperature:
    t_bd=y_dwarf_temperature(exner_theta_levels(i,j,k))
    ! Convert to potential temperature
    theta_out=t_bd / exner_theta_levels(i,j,k)
    y_dwarf_theta=theta_out

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END FUNCTION y_dwarf_theta

  REAL FUNCTION y_dwarf_recip_newt(i, j, k, exner_theta_levels)
    USE parkind1,                  ONLY: jpim, jprb       !DrHook
    USE yomhook,                   ONLY: lhook, dr_hook   !DrHook

    IMPLICIT NONE
    
    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='Y_DWARF_RECIP_NEWT'
    
    ! INPUT
    INTEGER :: i,j,k
    REAL    ::   exner_theta_levels(tdims_s%i_start:tdims_s%i_end,       &
         tdims_s%j_start:tdims_s%j_end,                  &
         tdims_s%k_start:tdims_s%k_end) 

    ! LOCAL VARIABLES
    ! This should be moved to a module
    REAL, PARAMETER :: stefan_boltzman=5.67E-8 ! SI: W/m**2/K**4
    REAL :: pressure, temperature
    REAL :: recip_tscale_out

    ! Function to calculate the estimated relaxation timescale
    ! for Newtonian cooling from Showman and Guillot (2002)

    ! 1.0 Start of function code: perform the calculation.
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ! Calculate Pressure
    pressure=(exner_theta_levels(i,j,k)**recip_kappa)*p_zero
    ! Also need the temperature at this point
    temperature=y_dwarf_temperature(exner_theta_levels(i,j,k))
    ! Now calcualte the Showman and Guillot (2002) estimated
    ! relaxation timescale
    ! Gravity is taken from gravity_mod, and Stefan Boltzman is 
    ! defined locally (the latter should be changed)
    recip_tscale_out=(pressure/g_theta(i,j,k))*                          &
         (cp/(4.0*stefan_boltzman*temperature**3.0))
    ! Set output
    y_dwarf_recip_newt=recip_tscale_out

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END FUNCTION y_dwarf_recip_newt

  ! ----------------------------------------------------------------------
  ! FUNCTIONS: Derive Y Dwarf Temperature
  ! ----------------------------------------------------------------------

  REAL FUNCTION y_dwarf_temperature(exner_theta_levels)
    ! Returns Pascal's polynomial fit T-P from Atmo
    USE parkind1,                  ONLY: jpim, jprb       !DrHook
    USE yomhook,                   ONLY: lhook, dr_hook   !DrHook

    IMPLICIT NONE
    
    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='Y_DWARF_TEMPERATURE'
    
    ! INPUT
    REAL, INTENT (IN) :: exner_theta_levels

    ! LOCAL VARIABLES
    REAL :: log_sigma, pressure
    REAL :: P_high, P_low
    ! OUTPUT
    REAL :: t_bd

    ! 1.0 Start of function code: perform the calculation.
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ! Fit valid from 1 to 1e7 Pa
    P_high=1.0e7
    P_low=1.0

    ! First construct the pressure variable
    pressure=(exner_theta_levels**recip_kappa)*p_zero
    IF (pressure < P_low) THEN
       log_sigma=LOG10(P_low/1.0e5)
    ELSE IF (pressure > P_high) THEN
       log_sigma=LOG10(P_high/1.0e5)
    ELSE
       log_sigma=LOG10(pressure/1.0e5)
    END IF
    ! Calculate Temperature
    t_bd=(527.802318                              &
         +357.598964*log_sigma                    &
         +143.297285*(log_sigma**2.0)             &
         +30.6223854*(log_sigma**3.0)             &
         -3.45983112*(log_sigma**4.0)             &
         -9.66674914*(log_sigma**5.0)             &
         -3.37366445*(log_sigma**6.0)             &
         +0.814368708*(log_sigma**7.0)            &
         +0.700909211*(log_sigma**8.0)            &
         +0.140662132*(log_sigma**9.0)            &
         +0.00930073262*(log_sigma**10.0))     
    ! Set output
    y_dwarf_temperature=t_bd

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END FUNCTION y_dwarf_temperature

END MODULE y_dwarf_forcing_mod
