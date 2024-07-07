! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE hd209458b_forcing_mod

  USE planet_constants_mod
  USE atm_fields_bounds_mod
  USE tforce_mod, ONLY: tf_HD209458b_Heng,                                   &
       tf_HD209458b_Heng_smooth,                                             &
       tf_HD209458b_iro
  USE trelax_mod, ONLY: tr_HD209458b_Iro
  USE model_domain_mod
  USE horiz_grid_mod
  USE umPrintMgr
  USE ereport_mod,               ONLY: ereport

  IMPLICIT NONE

  ! Description: 
  !   Module containing the subroutines/functions required for
  !   an idealised HD 209458b test. Including the derivation
  !   of the required equilibrium potential temperature,
  !   relaxation timescale and friction.
  !
  ! Method:
  !   The values are derived as described in Heng et al (2011)
  !   and adapted in Mayne et al (2014)
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
  ! This file belongs in section: IDEALISED
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !   This code is written to UMDP3 programming standards.

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='HD209458B_FORCING_MOD'

CONTAINS

  REAL FUNCTION hd209458b_theta(i, j, k, exner_theta_levels, tforce_number)
    USE parkind1,                  ONLY: jpim, jprb       !DrHook
    USE yomhook,                   ONLY: lhook, dr_hook   !DrHook

    IMPLICIT NONE    

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='HD209458B_THETA'
    
    ! INPUT
    INTEGER :: i,j,k
    REAL    ::   exner_theta_levels(tdims_s%i_start:tdims_s%i_end,           &
         tdims_s%j_start:tdims_s%j_end,                                      &
         tdims_s%k_start:tdims_s%k_end) 
    INTEGER :: tforce_number

    ! LOCAL VARIABLES
    REAL :: T_day, T_night ! The day/night side temperatures of Hot Jupiters
    REAL :: theta_out
    ! Function to calculate the HD209458b 
    ! potential temperature (Heng et al, 2011)
    ! Supports either the original or a smoothed profile


    ! 1.0 Start of function code: perform the calculation.
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ! Set the day and night side T-P profiles
    IF (tforce_number == tf_HD209458b_Heng .OR.                              &
         tforce_number == tf_HD209458b_Heng_smooth) THEN
       T_night=hd209458b_night_temp(exner_theta_levels(i,j,k),               &
            tforce_number)   
       T_day=hd209458b_day_temp(exner_theta_levels(i,j,k),                   &
            tforce_number)     
    ELSE IF (tforce_number == tf_HD209458b_iro) THEN
       T_night=iro_HD209458b_temp(exner_theta_levels(i,j,k),                 &
            tforce_number)
       T_day=T_night
    END IF

    ! Now apply the temperature, dependent on the location
    ! Note there is an error in Heng et al (2011) Eqn 26
    ! 90<phi<270 should be 90<theta<270
    ! Day side is 90<theta(longitude)<270
    ! Therefore, cos_theta_longitude <=0.0 
    ! Again using potential temperature so convert by /exner
    IF (Csxi1_p(i) <= 0.0) THEN
       theta_out=((T_night**4.0 -                                            &
            (T_day**4.0 - T_night**4.0)                                      &
            * Csxi1_p(i) * Csxi2_p(j))**0.25)                                &
            / exner_theta_levels(i,j,k) 
    ELSE
       theta_out=T_night / exner_theta_levels(i,j,k)
    END IF
    hd209458b_theta=theta_out

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END FUNCTION hd209458b_theta

  REAL FUNCTION hd209458b_recip_newt(i, j, k, exner_theta_levels,            &
       trelax_number)
    USE parkind1,                  ONLY: jpim, jprb       !DrHook
    USE yomhook,                   ONLY: lhook, dr_hook   !DrHook

    IMPLICIT NONE
    
    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='HD209458B_RECIP_NEWT'
    
    ! INPUT
    INTEGER :: i,j,k
    REAL    ::   exner_theta_levels(tdims_s%i_start:tdims_s%i_end,           &
         tdims_s%j_start:tdims_s%j_end,                                      &
         tdims_s%k_start:tdims_s%k_end) 
    INTEGER :: trelax_number, err_code

    ! LOCAL VARIABLES
    REAL :: recip_tscale_out

    ! Function to calculate the HD209458b prescribed
    ! newtonian relaxation timescale (reciprocal) 

    ! 1.0 Start of function code: perform the calculation.
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    IF (trelax_number == tr_HD209458b_Iro) THEN
       recip_tscale_out=recip_newt_tscale_iro(exner_theta_levels(i,j,k))
    ELSE
       WRITE(umMessage,'(A,I3)')                                             &
            '-HD 209458b relaxation timescale not supported:',               &
            trelax_number
       CALL umPrint(umMessage,src='hd209458b_forcing_mod.F90')
       err_code = 1
       CALL ereport("hd209458b_forcing_mod", err_code,                       &
            "Relaxation timescale not supported" )
    END IF
    hd209458b_recip_newt=recip_tscale_out

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END FUNCTION hd209458b_recip_newt

  ! ----------------------------------------------------------------------
  ! FUNCTIONS: Derive HD209458b Night-Side Temperature
  ! ----------------------------------------------------------------------

  REAL FUNCTION hd209458b_night_temp(exner_theta_levels, tforce_number)
    ! Function which returns the night side equilibrium temperature
    ! for HD209458b, smoothed or original version
    USE parkind1,                  ONLY: jpim, jprb       !DrHook
    USE yomhook,                   ONLY: lhook, dr_hook   !DrHook
    
    IMPLICIT NONE

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='HD209458B_NIGHT_TEMP'
    
    ! INPUT
    REAL, INTENT (IN) :: exner_theta_levels
    INTEGER :: tforce_number, err_code

    ! LOCAL VARIABLES
    REAL :: log_sigma, pressure, log_sigma_inactive
    REAL :: P_low, P_high
    REAL :: T_night_active, T_night_inactive
    REAL :: Temp_low, alpha
    ! OUTPUT
    REAL :: T_night

    ! 1.0 Start of function code: perform the calculation.
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ! First construct the pressure variable
    pressure=(exner_theta_levels**recip_kappa)*p_zero
    ! Set some constants for the smoothing
    IF (tforce_number ==  tf_HD209458b_Heng_smooth) THEN
       Temp_low=250.0
       alpha=0.10
    END IF

    ! Now derive night side temperature
    IF (tforce_number == tf_HD209458b_Heng .OR.                              &
         tforce_number == tf_HD209458b_Heng_smooth) THEN
       ! From Heng et al (2011)
       ! These are valid from 1microbar (0.1Pa) to 3488 bar (3488e5 Pa)
       ! Switch to inactive region at 10bar (1e6Pa)
       ! Set the pressure limits
       P_high=1.0e6
       P_low=0.1
       IF (tforce_number == tf_HD209458b_Heng_smooth)                        &
            P_low=100.0
       ! Calculate the log_sigma required for T_night_active
       IF (pressure >= P_high) THEN
          log_sigma=LOG10(P_high/1.0e5)
       ELSE IF (pressure <= P_low) THEN
          log_sigma=LOG10(P_low/1.0e5)
       ELSE
          log_sigma=LOG10(pressure/1.0e5)
       END IF
       T_night_active=(1388.2145+267.66586*log_sigma                         &
            -215.53357*(log_sigma**2.0)                                      &
            +61.814807*(log_sigma**3.0)                                      &
            +135.68661*(log_sigma**4.0)                                      &
            +2.0149044*(log_sigma**5.0)                                      &
            -40.907246*(log_sigma**6.0)                                      &
            -19.015628*(log_sigma**7.0)                                      &
            -3.8771634*(log_sigma**8.0)                                      &
            -0.38413901*(log_sigma**9.0)                                     &
            -0.015089084*(log_sigma**10.0))     
       ! Now calculate the Inactive Region
       IF (pressure > P_high) THEN
          ! Use a different log_sigma to avoid confusion
          log_sigma_inactive=LOG10(pressure/1.0e5)
          T_night_inactive=                                                  & 
               (5529.7168-6869.6504*log_sigma_inactive                       &
               +4142.7231*(log_sigma_inactive**2.0)                          &
               -936.23053*(log_sigma_inactive**3.0)                          &
               +87.120975*(log_sigma_inactive**4.0))  
       ELSE 
          T_night_inactive=0.0
       END IF
       ! Now construct the output temperature
       IF (pressure > P_high) THEN
          IF (tforce_number == tf_HD209458b_Heng) THEN
             T_night=T_night_inactive
          ELSE IF (tforce_number == tf_HD209458b_Heng_smooth) THEN
             T_night=T_night_active+100.0*                                   &
                  (1.0-EXP(-1.0*(LOG10(pressure)-LOG10(P_high))))
          END IF
       ELSE IF (pressure < P_low) THEN
          IF (tforce_number == tf_HD209458b_Heng) THEN
             T_night=T_night_active
          ELSE IF (tforce_number == tf_HD209458b_Heng_smooth) THEN
             T_night=                                                        &
                  MAX(T_night_active*                                        &
                  EXP(alpha*(LOG10(pressure)-LOG10(P_low))), Temp_low)
          END IF
       ELSE
          T_night=T_night_active     
       END IF
    ELSE
       WRITE(umMessage,'(A,I3)')                                             &
            '-HD 209458b nightside T-P not supported:',                      &
            tforce_number
       CALL umPrint(umMessage,src='hd209458b_forcing_mod.F90')
       err_code = 1
       CALL ereport("hd209458b_forcing_mod", err_code,                       &
            "HD 209458b nightside T-P profile not supported" )
    END IF
    hd209458B_night_temp=T_night

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END FUNCTION hd209458b_night_temp

! ----------------------------------------------------------------------
! FUNCTIONS: Derive HD209458b Day-Side Temperature
! ----------------------------------------------------------------------

  REAL FUNCTION hd209458b_day_temp(exner_theta_levels, tforce_number)
    ! Function which returns the day side equilibrium temperature
    ! for HD209458b
    USE parkind1,                  ONLY: jpim, jprb       !DrHook
    USE yomhook,                   ONLY: lhook, dr_hook   !DrHook
    
    IMPLICIT NONE

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='HD209458B_DAY_TEMP'
    
    ! INPUT
    REAL, INTENT (IN) :: exner_theta_levels
    INTEGER :: tforce_number,err_code 

    ! LOCAL VARIABLES
    REAL :: log_sigma, pressure, log_sigma_inactive
    REAL :: T_day_active, T_day_inactive
    REAL :: P_low, P_high
    REAL :: Temp_low, alpha
    ! OUTPUT
    REAL :: T_day

    ! 1.0 Start of function code: perform the calculation.
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ! First construct the pressure variable
    pressure=(exner_theta_levels**recip_kappa)*p_zero
    ! Set some constants for the smoothing
    IF (tforce_number ==  tf_HD209458b_Heng_smooth) THEN
       Temp_low=1000.0
       alpha=0.015
    END IF

    ! Now derive day side temperature
    IF (tforce_number == tf_HD209458b_Heng .OR.                              &
         tforce_number == tf_HD209458b_Heng_smooth) THEN
       ! From Heng et al (2011)
       ! These are valid from 1microbar (0.1Pa) to 3488 bar (3488e5 Pa)
       ! Switch to inactive region at 10bar (1e6Pa)
       ! Set the pressure limits
       P_high=1.0e6
       P_low=0.1
       IF (tforce_number == tf_HD209458b_Heng_smooth)                        &
            P_low=100.0
       ! Calculate the log_sigma required for T_day_active
       IF (pressure >= P_high) THEN
          log_sigma=LOG10(P_high/1.0e5)
       ELSE IF (pressure <= P_low) THEN
          log_sigma=LOG10(P_low/1.0e5)
       ELSE
          log_sigma=LOG10(pressure/1.0e5)
       END IF
       T_day_active=(2149.9581+4.1395571*(log_sigma)                         &
            -186.24851*(log_sigma**2.0)                                      &
            +135.52524*(log_sigma**3.0)                                      &
            +106.20433*(log_sigma**4.0)                                      &
            -35.851966*(log_sigma**5.0)                                      &
            -50.022826*(log_sigma**6.0)                                      &
            -18.462489*(log_sigma**7.0)                                      &
            -3.3319965*(log_sigma**8.0)                                      &
            -0.30295925*(log_sigma**9.0)                                     &
            -0.011122316*(log_sigma**10.0))  
       ! Now calculate the Inactive Region
       IF (pressure > P_high) THEN
          ! Use a different log_sigma to avoid confusion
          log_sigma_inactive=LOG10(pressure/1.0e5)
          T_day_inactive=                                                    &
               (5529.7168-6869.6504*log_sigma_inactive                       &
               +4142.7231*(log_sigma_inactive**2.0)                          &
               -936.23053*(log_sigma_inactive**3.0)                          &
               +87.120975*(log_sigma_inactive**4.0))  
       ELSE 
          T_day_inactive=0.0
       END IF
       ! Now construct the output temperature
       IF (pressure > P_high) THEN
          IF (tforce_number == tf_HD209458b_Heng) THEN
             T_day=T_day_inactive
          ELSE IF (tforce_number == tf_HD209458b_Heng_smooth) THEN
             T_day=T_day_active-120.0*                                       &
                  (1.0-EXP(-1.0*(LOG10(pressure)-LOG10(P_high))))
          END IF
       ELSE IF (pressure < P_low) THEN
          IF (tforce_number == tf_HD209458b_Heng) THEN
             T_day=T_day_active
          ELSE IF (tforce_number == tf_HD209458b_Heng_smooth) THEN
             T_day=                                                          &
                  MAX(T_day_active*                                          &
                  EXP(alpha*(LOG10(pressure)-LOG10(P_low))), Temp_low)
          END IF
       ELSE
          T_day=T_day_active     
       END IF
    ELSE
       WRITE(umMessage,'(A,I3)')                                             &
            '-HD 209458b dayside T-P not supported:',                        &
            tforce_number
       CALL umPrint(umMessage,src='hd209458b_forcing_mod.F90')
       err_code =1
       CALL ereport("hd209458b_forcing_mod", err_code,                       &
            "HD 209458b dayside T-P profile not supported" )
    END IF
    hd209458b_day_temp=T_day

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END FUNCTION hd209458b_day_temp

! ----------------------------------------------------------------------
! FUNCTIONS: IRO Temp HD209458b
! ----------------------------------------------------------------------

  REAL FUNCTION iro_HD209458b_temp(exner_theta_levels, tforce_number)
    ! Function which returns the Temperature for the radiative 
    ! equilibrium solution of Iro et al (2005) for HD209458b
    USE parkind1,                  ONLY: jpim, jprb       !DrHook
    USE yomhook,                   ONLY: lhook, dr_hook   !DrHook

    IMPLICIT NONE
    
    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='IRO_HD209458B_TEMP'
    
    ! INPUT
    REAL, INTENT (IN) :: exner_theta_levels
    INTEGER :: tforce_number,err_code  

    ! LOCAL VARIABLES
    REAL :: log_sigma, pressure
    REAL :: P_low, P_high
    ! OUTPUT
    REAL :: T_iro

    ! 1.0 Start of function code: perform the calculation.
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
    ! First construct the pressure variable
    pressure=(exner_theta_levels**recip_kappa)*p_zero
    IF (tforce_number == tf_HD209458b_iro) THEN
       ! Now set validity limits
       ! These are valid from 1microbar (0.1Pa) to 3488 bar (3488e5 Pa)
       P_high=3488.0e5
       P_low=0.1
       IF (pressure >= P_high) THEN
          log_sigma=LOG10(P_high/1.0e5)
       ELSE IF (pressure <= P_low) THEN
          log_sigma=LOG10(P_low/1.0e5)
       ELSE
          log_sigma=LOG10(pressure/1.0e5)
       END IF
       ! Now calculate the Iro Temperature
       T_iro=(1696.6986+132.23180*(log_sigma)                                &
            -174.30459*(log_sigma**2.0)                                      &
            +12.579612*(log_sigma**3.0)                                      &
            +59.513639*(log_sigma**4.0)                                      &
            +9.6706522*(log_sigma**5.0)                                      &
            -4.1136048*(log_sigma**6.0)                                      &
            -1.0632301*(log_sigma**7.0)                                      &
            +0.064400203*(log_sigma**8.0)                                    &
            +0.035974396*(log_sigma**9.0)                                    &
            +0.0025740066*(log_sigma**10.0))  
    ELSE
       WRITE(umMessage,'(A,I3)')                                             &
            '-HD 209458b IRO T-P not supported:',                            &
            tforce_number
       CALL umPrint(umMessage,src='hd209458b_forcing_mod.F90')
       err_code = 1
       CALL ereport("hd209458b_forcing_mod", err_code,                       &
            "HD 209458b IRO T-P profile not supported" )
    END IF
    iro_HD209458b_temp=T_iro

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END FUNCTION iro_HD209458b_temp

! ----------------------------------------------------------------------
! FUNCTIONS: Derive HD209458b  Relaxation Timescale Iro et al (2005)
! ----------------------------------------------------------------------

  REAL FUNCTION recip_newt_tscale_iro(exner_theta_levels)
    ! Function which returns the reciprocal of the relaxation timescale
    ! or radiative timescale at the pressure (exner_theta_levels) for
    ! the profile of Iro et al (2005) for HD209458b
    USE parkind1,                  ONLY: jpim, jprb       !DrHook
    USE yomhook,                   ONLY: lhook, dr_hook   !DrHook

    IMPLICIT NONE

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='RECIP_NEWT_TSCALE_IRO'

    ! INPUT
    REAL, INTENT (INOUT) :: exner_theta_levels

    ! LOCAL VARIABLES
    REAL :: log_sigma, pressure
    REAL :: P_low, P_high
    ! OUTPUT
    REAL :: recip_newt_tscale

    ! 1.0 Start of function code: perform the calculation.
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ! First construct the pressure variable
    pressure=(exner_theta_levels**recip_kappa)*p_zero

    ! P>10bar, inactive core
    ! Heng et al (2011) gives log(tau_rad)
    ! We need 1/T_rad, s**-1 (so conversion included)
    ! Heng et al (2011) States that the fit for relaxation time
    ! is valid for 10 microbar (1Pa) <P<8.5bar (8.5e5Pa)
    ! I am not sure what this means you should do between 
    ! 8.5 and 10 bar
    P_high=1.0e6 ! 10 bar
    P_low=1.0 ! 1e-4 bar
    IF (pressure < P_high) THEN
       ! If the pressure gets to low cap the relaxation time
       ! at the value at 10microbar
       IF (pressure < P_low) pressure=P_low
       ! Construct the pressure variable
       ! From Heng et al (2011) =log(P/1bar) (1bar=1.0e5 pa)
       log_sigma=LOG10(pressure/1.0e5)
       recip_newt_tscale=1.0/(10.0**(5.4659686+1.4940124*log_sigma           &
            +0.66079196*(log_sigma**2.0)+0.16475329*(log_sigma**3.0)         &
            +0.014241552*(log_sigma**4.0)))
    ELSE
       recip_newt_tscale=0.0
    END IF
    recip_newt_tscale_iro=recip_newt_tscale

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END FUNCTION recip_newt_tscale_iro

END MODULE hd209458b_forcing_mod
