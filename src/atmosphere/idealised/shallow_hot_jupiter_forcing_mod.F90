! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Temperature and wind forcing fuctions for shallow hot jupiter
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3

MODULE shallow_hot_jupiter_forcing_mod
  
  USE planet_constants_mod
  USE conversions_mod,           ONLY: pi
  USE level_heights_mod,         ONLY: r_theta_levels
  USE atm_fields_bounds_mod
  USE horiz_grid_mod
  USE model_domain_mod
  
  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
    ModuleName='SHALLOW_HOT_JUPITER_FORCING_MOD'

CONTAINS
  
  REAL FUNCTION shallow_hot_jupiter_theta(i, j, k, exner_theta_levels)
    USE parkind1,                  ONLY: jpim, jprb       !DrHook
    USE yomhook,                   ONLY: lhook, dr_hook   !DrHook

    IMPLICIT NONE
!
! Description:
!   Function to calculate Shallow Hot Jupiter potential temperature
!   (Menou & Rauscher, 2009)

! Subroutine arguments    
    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='SHALLOW_HOT_JUPITER_THETA'
    
    ! INPUT
    INTEGER :: i,j,k
    REAL    ::   exner_theta_levels(tdims_s%i_start:tdims_s%i_end,       &
                                    tdims_s%j_start:tdims_s%j_end,       &
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
    ! 
    
    ! 1.0 Start of function code: perform the calculation.
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
    ! Set local variables
    z_strat=2.0e6
    gamma_trop=2.0e-4
    Static_Stab=10.0
    T_min=1210.0
    T_surf=1600.0
    DT_eq_pole=300.0
    ! Find location and sigma at the stratosphere
    strat_loc=MINLOC(ABS(r_theta_levels(i,j,:)-Planet_radius-z_strat),1)
    ! Exner=(P/P_0)**kappa
    ! Sigma=(Exner/Exner_surf)**(1.0/kappa)
    sigma_strat=(exner_theta_levels(i,j,strat_loc)/&
         exner_theta_levels(i,j,0))**recip_kappa
    ! sigma_strat should be ~0.12 for Shallow Hot Jupiter
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
    theta_out=(T_vert - beta_trop * DT_eq_pole                    &
         * Csxi1_p(i)                                             &
         * Csxi2_p(j))                                            &
         /exner_theta_levels(i,j,k)
    shallow_hot_jupiter_theta=theta_out
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END FUNCTION shallow_hot_jupiter_theta
  
  REAL FUNCTION shallow_hot_jupiter_recip_newt()
    USE parkind1,                  ONLY: jpim, jprb       !DrHook
    USE yomhook,                   ONLY: lhook, dr_hook   !DrHook
    
    IMPLICIT NONE
!
! Description:
!   Function to calculate the Shallow Hot Jupiter (Menou & Rauscher, 2009)
! newtonian relaxation timescale (reciprocal)

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='SHALLOW_HOT_JUPITER_RECIP_NEWT'
    
    ! LOCAL VARIABLE
    REAL :: recip_tscale_out
    
    ! 1.0 Start of function code: perform the calculation.
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
    
    ! Set as a constant
    ! Heng et al (2011) use pi/Rotation period~1.731 days, 
    ! therefore 1.49558e5 seconds    
    recip_tscale_out=1.0/1.49558e5

    shallow_hot_jupiter_recip_newt=recip_tscale_out
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END FUNCTION shallow_hot_jupiter_recip_newt
  
END MODULE shallow_hot_jupiter_forcing_mod
