! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE gj1214b_forcing_mod

  USE planet_constants_mod, ONLY: p_zero, recip_kappa, cp, r
  USE atm_fields_bounds_mod, ONLY : tdims_s
  USE tforce_mod, ONLY: tf_GJ1214b, tf_GJ1214b_dT800
  USE trelax_mod, ONLY: tr_GJ1214b
  USE horiz_grid_mod, ONLY: Csxi1_p, Csxi2_p
  USE ereport_mod,               ONLY: ereport
  USE parkind1,                  ONLY: jpim, jprb       !DrHook
  USE yomhook,                   ONLY: lhook, dr_hook   !DrHook

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: gj1214b_theta, gj1214b_recip_newt, gj1214b_mean_temp

  ! Description: 
  !   Module containing the subroutines/functions required for
  !   GJ 1214b, a warm Neptune. Including the derivation
  !   of the required equilibrium potential temperature,
  !   and relaxation timescale.
  !
  ! Method:
  !   The values are derived as in Zhang and Showman 2017
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
  ! This file belongs in section: IDEALISED
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !   This code is written to UMDP3 programming standards.

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='GJ1214B_FORCING_MOD'

CONTAINS

  ! ----------------------------------------------------------------------
  ! FUNCTION: Get GJ1214b equilibrium potential temperature profile
  ! ----------------------------------------------------------------------  
  REAL FUNCTION gj1214b_theta(i, j, k, exner_theta_levels, tforce_number)

    IMPLICIT NONE    

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='GJ1214B_THETA'
    
    ! INPUT
    INTEGER :: i,j,k,tforce_number
    REAL    ::   exner_theta_levels(tdims_s%i_start:tdims_s%i_end,           &
         tdims_s%j_start:tdims_s%j_end,                                      &
         tdims_s%k_start:tdims_s%k_end) 

    ! LOCAL VARIABLES
    REAL :: Tmean  ! The mean temperature
    REAL :: dTeq   ! The difference in the equilibrium temperature
    REAL :: Tnight ! Nightside temperature
    REAL :: theta_out
    
    ! Function to calculate the GJ1214b equilibrium potential temperature
    ! profile following the method of Zhang and Showman 2017

    ! 1.0 Start of function code: perform the calculation.
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ! Calculate mean temperature T0
    ! (See Fig 2, Zhang and Showman 2017)
    Tmean = gj1214b_mean_temp(exner_theta_levels(i,j,k))
    
    ! Difference in equilibrium temperature - pressure dependent
    dTeq = gj1214b_delta_temp_eq(exner_theta_levels(i,j,k),tforce_number)
    
    ! Nightside temperature profile
    Tnight = Tmean - (dTeq/2.)
    
    ! Equilibrium (potential) temperature
    ! Equal to Tnight for nightside but depends
    ! on lat and lon on dayside
    ! Dayside is 90<longitude<270
    IF (Csxi1_p(i) <= 0.0) THEN
      theta_out = Tnight -                                                   &
          dTeq * Csxi1_p(i) * Csxi2_p(j)                                     
      theta_out = theta_out/exner_theta_levels(i,j,k)
    ELSE
      theta_out = Tnight/exner_theta_levels(i,j,k)
    END IF
    gj1214b_theta = theta_out
    
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END FUNCTION gj1214b_theta

  ! ----------------------------------------------------------------------
  ! FUNCTION: Get GJ1214b radiative timescale
  ! ----------------------------------------------------------------------  

  REAL FUNCTION gj1214b_recip_newt(i, j, k, exner_theta_levels)

    IMPLICIT NONE
    
    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='GJ1214B_RECIP_NEWT'
    
    ! INPUT
    INTEGER :: i,j,k
    REAL    ::   exner_theta_levels(tdims_s%i_start:tdims_s%i_end,           &
         tdims_s%j_start:tdims_s%j_end,                                      &
         tdims_s%k_start:tdims_s%k_end) 

    ! LOCAL VARIABLES
    REAL :: recip_tscale_out

    ! Function to calculate the GJ1214b prescribed
    ! newtonian relaxation timescale, following the method of 
    ! Zhang and Showman 2017

    ! 1.0 Start of function code: perform the calculation.
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    recip_tscale_out=recip_newt_tscale_gj1214b(exner_theta_levels(i,j,k))

    gj1214b_recip_newt=recip_tscale_out

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END FUNCTION gj1214b_recip_newt

  ! ----------------------------------------------------------------------
  ! FUNCTION: Derive GJ1214b Mean Temperature
  ! ----------------------------------------------------------------------  
  REAL FUNCTION gj1214b_mean_temp(exner_theta_levels)
    ! Function which returns the mean temperature for 
    ! GJ1214b, from Zhang and Showman 2017
    
    IMPLICIT NONE

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='GJ1214B_MEAN_TEMP'
    
    ! INPUT 
    REAL, INTENT(IN) :: exner_theta_levels
    
    ! LOCAL VARIABLES
    REAL :: pressure, ptrans
    ! OUTPUT
    REAL :: Tmean

    ! 1.0 Start of function code: perform the calculation.
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle) 

    ! First construct the pressure variable
    pressure=(exner_theta_levels**recip_kappa)*p_zero
    
    ! Transition pressure between low and high pressure polynomial fit
    ptrans = 1.0e6
    
    ! Use polynomial fit to Zhang and Showman 2017 mean temperature
    ! Valid between 10 Pa and 8e7 Pa.
    IF (pressure>=ptrans) THEN
      Tmean = 7.44675684e+02 +                                               &
          7.28917654e-06*(pressure) -                                        &
          1.06804684e-13*(pressure**2) +                                    &
          1.93590145e-21*(pressure**3) -                                     &
          2.22910784e-29*(pressure**4) +                                    &
          1.02507321e-37*(pressure**5)
    ELSE
      Tmean = 5.29329330e+02 +                                               &
          1.33581984e-03*(pressure) -                                        &
          4.19764028e-09*(pressure**2) +                                    &
          6.70345780e-15*(pressure**3) -                                     &
          4.75819938e-21*(pressure**4) +                                    &
          1.15162580e-27*(pressure**5)  
    END IF
    gj1214b_mean_temp=Tmean

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END FUNCTION gj1214b_mean_temp
  
  ! ----------------------------------------------------------------------
  ! FUNCTIONS: Compute the difference in equilibrium temperature
  ! ----------------------------------------------------------------------  
  REAL FUNCTION gj1214b_delta_temp_eq(exner_theta_levels,tforce_number)
    ! Function which returns the difference in equilibrium temperature for 
    ! GJ1214b, from Zhang and Showman 2017
    
    IMPLICIT NONE

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='GJ1214B_DELTA_TEMP_EQ'
    
    ! INPUT 
    INTEGER, INTENT(IN) :: tforce_number
    REAL, INTENT(IN) :: exner_theta_levels
    
    ! LOCAL VARIABLES
    REAL :: pressure, Pmax, Pmin
    ! OUTPUT
    REAL :: dTeq

    ! 1.0 Start of function code: perform the calculation.
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ! First construct the pressure variable
    pressure=(exner_theta_levels**recip_kappa)*p_zero
    
    ! Pressure limits of fit
    Pmin = 10.
    Pmax = 200e5
    
    dTeq = 0.
    ! dTeq is 600 K at top of atmosphere and decays linearly in logP
    ! Constant value of 600 K below Pmin and zero for above Pmax
    IF (pressure<=Pmin) THEN
      IF (tforce_number==tf_GJ1214b) THEN
        dTeq = 600.
      ELSE IF (tforce_number==tf_GJ1214b_dt800) THEN
        dTeq = 800.
      END IF
    ELSE IF (pressure<Pmax) THEN
      IF (tforce_number==tf_GJ1214b) THEN
        dTeq = -95.22252718*(LOG10(pressure)) +                              &
            695.22252718
      ELSE IF (tforce_number==tf_GJ1214b_dT800) THEN
        dTeq = -126.98412698*(LOG10(pressure)) +                             &
            926.98412698
      END IF
    END IF
    gj1214b_delta_temp_eq=dTeq

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END FUNCTION gj1214b_delta_temp_eq

  ! ----------------------------------------------------------------------
  ! FUNCTIONS: Compute the radiative timescale for GJ1214b
  ! ----------------------------------------------------------------------  
  REAL FUNCTION recip_newt_tscale_gj1214b(exner_theta_levels)
    ! Function which returns the radiative timescale for 
    ! GJ1214b, from Zhang and Showman 2017
    
    IMPLICIT NONE

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='RECIP_NEWT_TSCALE_GJ1214B'

    ! INPUT 
    REAL :: exner_theta_levels
    
    ! LOCAL VARIABLES
    REAL :: pressure, Pmin, Pmax
    REAL :: tscale_h2, tscale
    REAL, PARAMETER :: r_h2 = 41158.4 ! J/K/kg
    REAL, PARAMETER :: cp_h2 = 1.4569e4 ! J/kg/K
    ! OUTPUT
    REAL :: recip_newt_tscale

    ! 1.0 Start of function code: perform the calculation.
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ! First construct the pressure variable
    pressure=(exner_theta_levels**recip_kappa)*p_zero

    ! Set pressure limits
    Pmin = 1.0e2
    Pmax = 1.0e6
    
    ! Molecular hydrogen radiative timescale - pressure dependent
    IF (pressure<Pmin) THEN
      tscale_h2 = 1.0e4
    ELSE IF (pressure>Pmax) THEN
      tscale_h2 = 1.0e7
    ELSE 
      tscale_h2 = 10**(5./2.)                                                &
        *pressure**(3./4.)
    END IF

    ! Radiative timescale
    ! Scaled from H2 timescale
    tscale = tscale_h2*(cp/cp_h2)                                            &
        *(r_h2/r)
    
    ! Reciprocal radiative timescale
    recip_newt_tscale_gj1214b = 1.0/tscale
    
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END FUNCTION recip_newt_tscale_gj1214b

END MODULE gj1214b_forcing_mod
