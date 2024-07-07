! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
MODULE earth_like_forcing_mod
  
  USE planet_constants_mod
  USE conversions_mod,           ONLY: pi
  USE level_heights_mod,         ONLY: r_theta_levels
  USE atm_fields_bounds_mod
  USE horiz_grid_mod
  USE model_domain_mod  

  IMPLICIT NONE
    
  ! Description: 
  !   Module containing the subroutines/functions required for
  !   an idealised Earth-Like test. Including the derivation
  !   of the required equilibrium potential temperature,
  !   relaxation timescale and friction.
  !
  ! Method:
  !   The values are derived as described in Menou & Rauscher (2009)
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
  ! This file belongs in section: IDEALISED
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !   This code is written to UMDP3 programming standards.  

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EARTH_LIKE_FORCING_MOD'

CONTAINS
  
  REAL FUNCTION earth_like_theta(i, j, k, exner_theta_levels)
    USE parkind1,                  ONLY: jpim, jprb       !DrHook
    USE yomhook,                   ONLY: lhook, dr_hook   !DrHook

    IMPLICIT NONE
    
    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='EARTH_LIKE_THETA'
    
    ! INPUT
    INTEGER :: i,j,k
    REAL    ::   exner_theta_levels(tdims_s%i_start:tdims_s%i_end,       &
         tdims_s%j_start:tdims_s%j_end,                  &
         tdims_s%k_start:tdims_s%k_end) 
    ! LOCAL VARIABLES
    REAL :: T_min, T_surf  ! Minimum/Stratospheric & surface temperature
    REAL :: DT_eq_pole     ! Equator-Pole Temp diff
    REAL :: Static_Stab    ! Static Stability temperature
    REAL :: z_strat        ! Height (m) of Stratosphere
    REAL :: gamma_trop     ! Lapse Rate
    INTEGER :: strat_loc   ! The location of the Stratosphere (index)
    REAL :: sigma_strat    ! Sigma value of the Stratosphere
    REAL :: sigma          ! The relative pressure of the level
    REAL :: T_vert         ! Stratospheric Temperature profile
    REAL :: beta_trop      ! Horizontal temperature gradient in troposphere
    REAL :: theta_out
    ! Function to calculate Earth-Like 
    ! potential temperature (Menou & Rauscher, 2009)
    
    ! 1.0 Start of function code: perform the calculation.
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
    
    ! Set local variables
    z_strat=1.2e4
    gamma_trop=6.5e-3
    Static_Stab=2.0
    T_min=212.0
    T_surf=288.0
    DT_eq_pole=60.0
    ! Find location and sigma at the stratosphere
    strat_loc=MINLOC(ABS(r_theta_levels(i,j,:)-Planet_radius-z_strat),1)
    ! Exner=(P/P_0)**kappa
    ! Sigma=(Exner/Exner_surf)**(1.0/kappa)
    sigma_strat=(exner_theta_levels(i,j,strat_loc)/&
         exner_theta_levels(i,j,0))**recip_kappa
    ! sigma_strat should be ~0.22 for Earth
    ! We require Sigma coordinates for this forcing
    ! NOTE: Exner=(P/P0)**1/kappa /= Sigma = P/P_surf!
    sigma=(exner_theta_levels(i,j,k)/&
         exner_theta_levels(i,j,0))**recip_kappa
    IF (r_theta_levels(i,j,k)-Planet_radius <= z_strat) THEN
       ! In the Stratosphere
       beta_trop=SIN(pi*(sigma-sigma_strat)/(2.0*(1.0-sigma_strat)))
       T_vert=T_surf-gamma_trop*                                             &
            (z_strat+(r_theta_levels(i,j,k)-Planet_radius-z_strat)*0.5) +    &
            ((gamma_trop*(r_theta_levels(i,j,k)-                             &
            Planet_radius-z_strat)*0.5)**2.0+Static_Stab**2.0)**0.5
    ELSE
       ! Above Stratosphere
       beta_trop=0.0
       T_vert=T_surf-gamma_trop*z_strat+Static_Stab
    END IF
    ! Recall Using Potential TEMPERATURE
    ! Therefore, must convert the temperature to potential 
    ! temperature using the exner function (exner*theta=Temp)
    ! Used cos(x-180)=-cos(x)
    theta_out=(T_vert + beta_trop * DT_eq_pole                    &
         * ((1.0/3.0)-Snxi2_p(j)*Snxi2_p(j)))                  &
          /exner_theta_levels(i,j,k)
    earth_like_theta=theta_out

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END FUNCTION earth_like_theta
  
  REAL FUNCTION earth_like_recip_newt()
    USE parkind1,                  ONLY: jpim, jprb       !DrHook
    USE yomhook,                   ONLY: lhook, dr_hook   !DrHook

    IMPLICIT NONE
    
    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='EARTH_LIKE_RECIP_NEWT'
    
    ! LOCAL VARIABLE
    REAL :: recip_tscale_out
    
    ! Function to calculate the Earth-Like  (Menou & Rauscher, 2009)
    ! newtonian relaxation timescale (reciprocal)
    
    ! 1.0 Start of function code: perform the calculation.
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
    
    ! Set as a constant
    ! Heng et al (2011) use constant tau_rad
    ! of 15 days, therefore  1.296e6 seconds
    recip_tscale_out=1.0/1.296e6
    earth_like_recip_newt=recip_tscale_out

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END FUNCTION earth_like_recip_newt
  
END MODULE earth_like_forcing_mod
