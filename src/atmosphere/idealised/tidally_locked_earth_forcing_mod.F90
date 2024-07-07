! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE tidally_locked_earth_forcing_mod

  USE planet_constants_mod
  USE atm_fields_bounds_mod
  USE horiz_grid_mod
  USE model_domain_mod
  
  IMPLICIT NONE

  ! Description: 
  !   Module containing the subroutines/functions required for
  !   an idealised Held-Suarez Tidally-Locked Earthtest. 
  !   Including the derivation
  !   of the required equilibrium potential temperature,
  !   relaxation timescale and friction.
  !
  ! Method:
  !   The values are derived as described in Heng et al (2011)
  !   based on the test of Merlis & Schneider (2011)
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
  ! This file belongs in section: IDEALISED
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !   This code is written to UMDP3 programming standards.

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
    ModuleName='TIDALLY_LOCKED_EARTH_FORCING_MOD'

CONTAINS
  
  REAL FUNCTION tidally_locked_earth_theta(i, j, k, exner_theta_levels)
    USE parkind1,                  ONLY: jpim, jprb       !DrHook
    USE yomhook,                   ONLY: lhook, dr_hook   !DrHook

    IMPLICIT NONE

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='TIDALLY_LOCKED_EARTH_THETA'

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
    ! Function to calculate the Tidally-Locked Earth prescribed
    ! potential temperature (Merlis & Schneider, 2010)

    ! 1.0 Start of function code: perform the calculation.
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ! Set local variables
    Static_Stab=10.0
    DT_eq_pole=60.0
    T_min=200.0
    T_surf=315.0
    ! Used cos(x-180)=-cos(x)
    ! Missing Final pressure term as using 
    ! potential temperature
    theta1 = T_min / exner_theta_levels(i,j,k)
    theta2 = T_surf - DT_eq_pole                            &
         * Csxi1_p(i)                                       &
         * Csxi2_p(j)                                       &
         - Static_Stab                                      &
         * LOG(exner_theta_levels(i,j,k))                   &
         * recip_kappa                                      &
         * Csxi2_p(j)                                       &
         * Csxi2_p(j)          
    theta_out = MAX(theta1, theta2)
    tidally_locked_earth_theta=theta_out

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END FUNCTION tidally_locked_earth_theta
    
END MODULE tidally_locked_earth_forcing_mod
