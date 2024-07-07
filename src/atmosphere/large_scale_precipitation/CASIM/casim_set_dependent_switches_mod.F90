! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Cloud Aerosol Interacting Microphysics (CASIM)
! Set variables in module casim_switches dependent on values in
! run_precip namelist

MODULE casim_set_dependent_switches_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='CASIM_SET_DEPENDENT_SWITCHES_MOD'

CONTAINS

SUBROUTINE casim_set_dependent_switches

USE mphys_inputs_mod, ONLY: l_mcr_qrain, l_mcr_qcf2, l_mcr_qgraup,             &
                            l_psd, l_psd_global, casim_moments_choice,         &
                            casim_aerosol_process_level,                       &
                            casim_aerosol_couple_choice,                       &
                            l_separate_process_rain
                          
USE casim_switches, ONLY: l_mp_cloudnumber, l_mp_rainnumber, l_mp_rain3mom,    &
                          l_mp_icenumber, l_mp_snownumber, l_mp_snow3mom,      &
                          l_mp_graupnumber, l_mp_graup3mom,                    &
                          l_mp_activesolliquid, l_mp_activesolrain,            &
                          l_mp_activeinsolice, l_mp_activesolice,              &
                          l_mp_activeinsolliquid, l_mp_activesolnumber,        &
                          l_mp_activeinsolnumber, l_casim_warm_only,           &
                          n_casim_progs, casim_moments_option,                 &
                          l_fix_aerosol, l_tracer_aerosol, l_ukca_aerosol,     &
                          l_ukca_feeding_in, l_ukca_feeding_out,               &
                          cloud_mom, rain_mom, ice_mom, snow_mom, graup_mom,   &
                          no_moments, single_moment, double_moment,            &
                          triple_moment, no_processing, passive_processing,    &
                          full_processing, passive_ice_only,                   &
                          passive_liquid_only

USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,    ONLY: lhook, dr_hook
USE parkind1,   ONLY: jprb, jpim
USE umprintmgr, ONLY: newline

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!   This routine overrides default values held in casim_switches
!   with values depending on the run_precip namelist input. Called from readlsta
!   reconfiguration or scm_shell after run_precip namelist read in.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Precipitation (CASIM)
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.1.
!
! ------------------------------------------------------------------------------

! Local parameters for the different settings are contained below

! Parameters for casim_moments_choice
INTEGER, PARAMETER :: all_single_moment = 0
INTEGER, PARAMETER :: all_double_moment = 1
INTEGER, PARAMETER :: warm_cloud2_rain3 = 2
INTEGER, PARAMETER :: warm_cloud1_rain2 = 3
INTEGER, PARAMETER :: warm_cloud1_rain1 = 4
INTEGER, PARAMETER :: warm_cloud1_rain3 = 5
INTEGER, PARAMETER :: warm_cloud2_rain2 = 6
INTEGER, PARAMETER :: double_no_graupel = 7
INTEGER, PARAMETER :: all_moments_on    = 8

! Parameters for casim_aerosol_couple_choice
INTEGER, PARAMETER :: fixed_aerosol      = 0
INTEGER, PARAMETER :: tracer_aerosol     = 1
INTEGER, PARAMETER :: ukca_aerosol_in    = 2
INTEGER, PARAMETER :: ukca_aerosol_inout = 3

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CASIM_SET_DEPENDENT_SWITCHES'

INTEGER :: errorstatus            ! Return code : 0 Normal Exit : >0 Error
CHARACTER(LEN=errormessagelength) :: cmessage
                                  ! Error message if Errorstatus /=0

LOGICAL ::  l_process             ! True if aerosol processing
LOGICAL ::  l_passivenumbers      ! True if aerosol processing and passive
                                  ! numbers in liquid
LOGICAL ::  l_passivenumbers_ice  ! True if aerosol processing and passive
                                  ! numbers in ice

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!-------------------------------------------------------------------------------
! 1. Initialization of variables
!-------------------------------------------------------------------------------
! Initialise error status to zero (normal/no error)
! Initialise number of extra casim prognostics to 0
errorstatus   = 0
n_casim_progs = 0

! Set all moments to zero initially and then modify depending on
! casim_moments_choice
cloud_mom = 0
rain_mom  = 0
ice_mom   = 0
snow_mom  = 0
graup_mom = 0

!-------------------------------------------------------------------------------
! 2.1 Select moments for microphysical species dependent on casim_moments_choice

! Note: casim_moments_option is moments for cloud, rain, ice, snow, graupel
! It is required as an integer in this form for the CASIM code itself
!-------------------------------------------------------------------------------

SELECT CASE (casim_moments_choice)

CASE (all_single_moment)
  ! Equivalent to current 3D setup in convective-scale models
  casim_moments_option = 11111
  cloud_mom            = 1
  rain_mom             = 1
  ice_mom              = 1
  snow_mom             = 1
  graup_mom            = 1

CASE (all_double_moment)
  ! All species active and double-moment
  casim_moments_option = 22222
  cloud_mom            = 2
  rain_mom             = 2
  ice_mom              = 2
  snow_mom             = 2
  graup_mom            = 2

CASE (warm_cloud2_rain3)
  ! Warm rain only
  ! Double-moment cloud
  ! Triple-moment rain
  casim_moments_option = 23000
  cloud_mom            = 2
  rain_mom             = 3
  l_casim_warm_only    = .TRUE.

CASE (warm_cloud1_rain2)
  ! Warm rain only
  ! Single-moment cloud
  ! Double-moment rain
  casim_moments_option = 12000
  cloud_mom            = 1
  rain_mom             = 2
  l_casim_warm_only    = .TRUE.

CASE (warm_cloud1_rain1)
  ! Warm rain only
  ! Single-moment cloud
  ! Single-moment rain
  casim_moments_option = 11000
  cloud_mom            = 1
  rain_mom             = 1
  l_casim_warm_only    = .TRUE.

CASE (warm_cloud1_rain3)
  ! Warm rain only
  ! Single-moment cloud
  ! Triple-moment rain
  casim_moments_option = 13000
  cloud_mom            = 1
  rain_mom             = 3
  l_casim_warm_only    = .TRUE.

CASE (warm_cloud2_rain2)
  ! Warm rain only
  ! Double-moment cloud
  ! Double-moment rain
  casim_moments_option = 22000
  cloud_mom            = 2
  rain_mom             = 2
  l_casim_warm_only    = .TRUE.

CASE (double_no_graupel)
  ! No graupel
  ! Double-moment cloud
  ! Double-moment rain
  ! Double-moment ice
  ! Double-moment snow
  casim_moments_option = 22220
  cloud_mom            = 2
  rain_mom             = 2
  ice_mom              = 2
  snow_mom             = 2

CASE (all_moments_on)
  ! All microphysics prognostics on
  ! Useful for prognostic variable testing
  casim_moments_option = 23233
  cloud_mom            = 2
  rain_mom             = 3
  ice_mom              = 2
  snow_mom             = 3
  graup_mom            = 3

CASE DEFAULT
  ! If this occurs then the user has clearly set up the namelist
  ! incorrectly. Throw an error message and stop the UM immediately
  errorstatus = 1

  WRITE(cmessage, '(A,I0)')                                                    &
    'Variable casim_moments_choice in the run_precip namelist has '// newline//&
    'been set incorrectly. Please check namelist input and correct'// newline//&
    'It should take an integer value >= 0 and agree with one of   '// newline//&
    'the available choices: see documentation for details.'        // newline//&
    'casim_moments_choice =', casim_moments_choice

    CALL ereport('casim_set_dependent_switches', errorstatus, cmessage)

END SELECT
!-------------------------------------------------------------------------------
! 2.2 Select level of aerosol coupling based on casim_aerosol_couple_choice
!-------------------------------------------------------------------------------

SELECT CASE (casim_aerosol_couple_choice)

! N.B. Most defensive programming and checking for casim aerosol options is 
! performed in init_casim_run due to ensuring other information required
! has been loaded in by other namelists (which may be after this point)

CASE (fixed_aerosol)
  ! Use a fixed value read from mphys_inputs_mod
  l_fix_aerosol = .TRUE. ! Other related switches already are false

CASE (tracer_aerosol)
  ! Use tracer aerosol data
  l_tracer_aerosol = .TRUE. ! Other related switches already are false

CASE (ukca_aerosol_in)
  ! Use UKCA aerosol, inwards only
  l_ukca_aerosol    = .TRUE. 
  l_ukca_feeding_in = .TRUE. ! UKCA set to feed in only. Other switches
                             ! are already false

CASE (ukca_aerosol_inout)
  ! Use UKCA aerosol in and out
  l_ukca_aerosol     = .TRUE. 
  l_ukca_feeding_in  = .TRUE. 
  l_ukca_feeding_out = .TRUE. ! UKCA set to feed in and out. Other switches
                              ! are already false

CASE DEFAULT 
  ! If this occurs then the user has clearly set up the namelist
  ! incorrectly. Throw an error message and stop the UM immediately
  errorstatus = 1
  WRITE(cmessage, '(A,I0)')                                                    &
    'Variable casim_aerosol_couple_choice in the run_precip '      // newline//&
    'namelist has been set incorrectly. Please check namelist'     // newline//&
    'input and correct. It should take an integer value >= 0 and'  // newline//&
    'agree with one of the available choices: see documentation '  // newline//&
    'for details. casim_aerosol_couple_choice =', casim_aerosol_couple_choice

    CALL ereport('casim_set_dependent_switches', errorstatus, cmessage)

END SELECT

IF ( l_ukca_aerosol ) THEN

    ! UKCA has not yet been coupled to CASIM. The rose settings should
    ! have prevented this from happening, but if not the model needs
    ! to throw an error immediately

    errorstatus = 1
    cmessage    = 'CASIM has been set to be coupled to UKCA    ' //newline//&
                  'aerosol. This option is not yet available,  ' //newline//&
                  'so for your model to run, please choose an  ' //newline//&
                  'alternative aerosol input option in the     ' //newline//&
                  'CASIM settings of the run_precip namelist.'
  
    CALL ereport('casim_set_dependent_switches', errorstatus, cmessage)

END IF ! l_ukca_aerosol

!-------------------------------------------------------------------------------
! 2.3 Select moments for microphysical species dependent on casim_aerosol_choice
!-------------------------------------------------------------------------------

SELECT CASE (casim_aerosol_process_level)

CASE (no_processing)
  l_process            = .FALSE.
  l_passivenumbers     = .FALSE.
  l_passivenumbers_ice = .FALSE.

CASE (passive_processing)
  l_process            = .TRUE.
  l_passivenumbers     = .TRUE.
  l_passivenumbers_ice = .TRUE.

CASE (full_processing)
  l_process            = .TRUE.
  l_passivenumbers     = .FALSE.
  l_passivenumbers_ice = .FALSE.

CASE (passive_ice_only)
  l_process            = .TRUE.
  l_passivenumbers     = .FALSE.
  l_passivenumbers_ice = .TRUE.

CASE (passive_liquid_only)
  l_process            = .TRUE.
  l_passivenumbers     = .TRUE.
  l_passivenumbers_ice = .FALSE.

CASE DEFAULT
  ! If this occurs then the user has clearly set up the namelist
  ! incorrectly. Throw an error message and stop the UM immediately
  errorstatus = 1
  WRITE(cmessage, '(A,I0)' )                                                   &
    'Variable casim_aerosol_process_level in the run_precip       '// newline//&
    'namelist has been set incorrectly. Please check namelist     '// newline//&
    'input and correct It should take an integer value >= 0 and   '// newline//&
    'agree with one of the available choices: see documentation   '// newline//&
    'for details. casim_aerosol_process_level =', casim_aerosol_process_level
  

  CALL ereport('casim_set_dependent_switches', errorstatus, cmessage)

END SELECT

! Set up prognostic switches based on level of aerosol procesing selected
! and the equivalent on the CASIM repository

IF (l_process) THEN

  IF (l_casim_warm_only) THEN

    l_mp_activesolliquid   = .TRUE.
    l_mp_activeinsolliquid = .FALSE.
    l_mp_activesolice      = .FALSE.
    l_mp_activeinsolice    = .FALSE.

    n_casim_progs = n_casim_progs + 1

  ELSE ! l_casim_warm_only

    l_mp_activesolliquid   = .TRUE.
    l_mp_activeinsolliquid = .TRUE.
    l_mp_activesolice      = .TRUE.
    l_mp_activeinsolice    = .TRUE.

    n_casim_progs = n_casim_progs + 4
 
  END IF ! l_casim_warm_only

  IF ( l_passivenumbers ) THEN
    l_mp_activesolnumber = .TRUE.
    n_casim_progs        = n_casim_progs + 1
  ELSE
    l_mp_activesolnumber = .FALSE.
  END IF ! l_passivenumbers
    
  IF ( l_passivenumbers_ice ) THEN
    l_mp_activeinsolnumber = .TRUE.
    n_casim_progs          = n_casim_progs + 1
  ELSE
    l_mp_activeinsolnumber = .FALSE.
  END IF ! l_passivenumbers_ice

  IF (l_separate_process_rain) THEN
    ! Rain is in its own category for processing
    l_mp_activesolrain = .TRUE.
    n_casim_progs      = n_casim_progs + 1
  ELSE
    ! No extra prognostic required for rain
    l_mp_activesolrain = .FALSE.
  END IF ! l_separate_process_rain

ELSE ! l_process

  l_mp_activesolliquid    = .FALSE.
  l_mp_activesolrain      = .FALSE.
  l_mp_activeinsolice     = .FALSE.
  l_mp_activesolice       = .FALSE.
  l_mp_activeinsolliquid  = .FALSE.
  l_mp_activesolnumber    = .FALSE.
  l_mp_activeinsolnumber  = .FALSE.

  ! ( Number of CASIM prognostics does not need increasing )

END IF ! l_process

!-------------------------------------------------------------------------------
! 3. User input error checking and correction

! While the existing (3D) microphysics and CASIM exist side by side in the UM
! code, options set for the existing microphysics will not work for CASIM.

! The GUI should check inputs and produce a warning, but adding some checking
! code here is a useful 'belt and braces' approach to ensure there are no
! unwanted effects of the code, especially if the user just modifies the 
! namelist input, without using the GUI.
!-------------------------------------------------------------------------------

IF ( rain_mom == no_moments .AND. l_mcr_qrain ) THEN

! Correct conflict between no-rain and prognostic rain

  errorstatus = -100
  cmessage    =                                                     newline//&
  'Prognostic rain has been switched on (variable l_mcr_qrain   '// newline//&
  'set to .TRUE. in the run_precip namelist), but this conflicts'// newline//&
  'with your choice of a no-rain simulation in CASIM. This has  '// newline//&
  'been corrected. To prevent this warning appearing, set       '// newline//&
  'l_mcr_qrain to .FALSE.'

  CALL ereport('casim_set_dependent_switches', errorstatus, cmessage)

  l_mcr_qrain = .FALSE. ! Correct the error

ELSE IF ( rain_mom > no_moments                         .AND.                  &
          casim_moments_choice > all_single_moment      .AND.                  &
          .NOT. l_mcr_qrain )                           THEN

! Correct conflict between rain-included and no prognostic rain
! Do not correct if casim_moments_choice set to 0.

  errorstatus = -100
  cmessage    =                                                     newline//&
  'Prognostic rain has been switched off (variable l_mcr_qrain  '// newline//&
  'set to .FALSE. in the run_precip namelist), but this         '// newline//&
  'conflicts with your choice of a simulation in CASIM including'// newline//&
  'rain. This has been corrected. To prevent this warning       '// newline//&
  'appearing, set l_mcr_qrain to .TRUE.'

  CALL ereport('casim_set_dependent_switches', errorstatus, cmessage)


  l_mcr_qrain = .TRUE. ! Correct the error

END IF ! rain_mom and l_mcr_qrain

IF ( graup_mom == no_moments .AND. l_mcr_qgraup ) THEN

! Correct conflict between no graupel in CASIM and prognostic graupel

  errorstatus = -100
  cmessage    =                                                     newline//&
  'Prognostic graupel has been switched on (l_mcr_qgraup set to '// newline//&
  '.TRUE. in the run_precip namelist), but this conflicts with  '// newline//&
  'your choice of a no-graupel simulation in CASIM. This has    '// newline//&
  'been corrected. To prevent this warning appearing, set       '// newline//&
  'l_mcr_qgraup to .FALSE.'

  CALL ereport('casim_set_dependent_switches', errorstatus, cmessage)

  l_mcr_qgraup = .FALSE. ! Correct the error

ELSE IF ( graup_mom > no_moments                         .AND.                 &
          casim_moments_choice > all_single_moment       .AND.                 &
          .NOT. l_mcr_qgraup )                           THEN

! Correct conflict between graupel-included and no prognostic graupel
! Do not correct if casim_moments_choice set to 0.

  errorstatus = -100
  cmessage    =                                                     newline//&
  'Prognostic graupel has been switched off (l_mcr_qgraup set to'// newline//&
  '.FALSE. in the run_precip namelist), but this conflicts with '// newline//&
  'your choice of a graupel-included simulation in CASIM. This  '// newline//&
  'has been corrected. To prevent this warning appearing, set   '// newline//&
  'l_mcr_qgraup to .TRUE.'

  CALL ereport('casim_set_dependent_switches', errorstatus, cmessage)

  l_mcr_qgraup = .TRUE. ! Correct the error

END IF ! graup_mom and l_mcr_qgraup

! N.B. CASIM uses the second ice cloud prognostic, which is in mphys_inputs_mod
! but not explicitly set in the run_precip namelist as it has not been used in
! any Wilson and Ballard microphysics jobs in the operational model. Rather than
! include it back in the namelist and have an extra option for users to worry
! about getting correct, it is far simpler to just set it on here, as it will 
! always be .TRUE. for CASIM simulations where ice_mom > 0.
IF ( ice_mom > no_moments ) l_mcr_qcf2 = .TRUE.

! The second ice cloud prognostic which we have just turned on does not work 
! with the generic ice particle size distribution (l_psd or l_psd_global set 
! to .TRUE.). Check whether the user has either of these turned on and if so, 
! throw a warning message and correct.

IF ( l_psd .OR. l_psd_global ) THEN

  errorstatus = -100
  cmessage    =                                                     newline//&
  'The generic ice particle size distribution has been turned on'// newline//&
  'in your run (l_psd and possibly l_psd_global set to .TRUE.   '// newline//&
  'in the run_precip namelist), but this does not work with     '// newline//&
  'CASIM. This has been corrected. To remove this warning please'// newline//&
  'set both of these logicals to .FALSE.'

  CALL ereport('casim_set_dependent_switches', errorstatus, cmessage)

  ! Now correct the error
  l_psd        = .FALSE.
  l_psd_global = .FALSE.

END IF ! l_psd or l_psd_global

!------------------------------------------------------------------------------
! 4.1 Determine (liquid) cloud moment and initialise switches
!------------------------------------------------------------------------------

IF (cloud_mom == double_moment) THEN
  l_mp_cloudnumber = .TRUE.
  n_casim_progs    = n_casim_progs + 1
END IF

!-------------------------------------------------------------------------------
! 4.2 Determine ice cloud moment and initialise switches
!------------------------------------------------------------------------------

IF (ice_mom == double_moment) THEN
  l_mp_icenumber = .TRUE.
  n_casim_progs  = n_casim_progs + 1
END IF

!------------------------------------------------------------------------------
! 4.3 Determine rain moment and initialise switches
!------------------------------------------------------------------------------
SELECT CASE (rain_mom)

CASE (single_moment) ! Single-moment rain

  l_mp_rainnumber = .FALSE.
  l_mp_rain3mom   = .FALSE.

CASE (double_moment) ! Double-moment rain

  l_mp_rainnumber = .TRUE.
  l_mp_rain3mom   = .FALSE.
  n_casim_progs   = n_casim_progs + 1

CASE (triple_moment) ! Triple-moment rain

  l_mp_rainnumber = .TRUE.
  l_mp_rain3mom   = .TRUE.
  n_casim_progs   = n_casim_progs + 2

END SELECT

!------------------------------------------------------------------------------
! 4.4 Determine snow moment and initialise switches
!------------------------------------------------------------------------------
SELECT CASE (snow_mom)

CASE (no_moments) ! No snow

  l_mp_snownumber = .FALSE.
  l_mp_snow3mom   = .FALSE.

CASE (single_moment) ! Single-moment snow

  l_mp_snownumber = .FALSE.
  l_mp_snow3mom   = .FALSE.

CASE (double_moment) ! Double-moment snow

  l_mp_snownumber = .TRUE.
  l_mp_snow3mom   = .FALSE.
  n_casim_progs   = n_casim_progs + 1

CASE (triple_moment) ! Triple-moment snow

  l_mp_snownumber = .TRUE.
  l_mp_snow3mom   = .TRUE.
  n_casim_progs   = n_casim_progs + 2

END SELECT

!------------------------------------------------------------------------------
! 4.5 Determine graupel moment and initialise switches
!------------------------------------------------------------------------------
SELECT CASE (graup_mom)

CASE (no_moments) ! No graupel

  l_mp_graupnumber = .FALSE.
  l_mp_graup3mom   = .FALSE.

CASE (single_moment) ! Single-moment graupel

  l_mp_graupnumber = .FALSE.
  l_mp_graup3mom   = .FALSE.

CASE (double_moment) ! Double-moment graupel

  l_mp_graupnumber = .TRUE.
  l_mp_graup3mom   = .FALSE.
  n_casim_progs    = n_casim_progs + 1

CASE (triple_moment) ! Triple-moment graupel

  l_mp_graupnumber = .TRUE.
  l_mp_graup3mom   = .TRUE.
  n_casim_progs    = n_casim_progs + 2

END SELECT

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE casim_set_dependent_switches

!==============================================================================
! Break between subroutines
!==============================================================================

SUBROUTINE casim_print_dependent_switches

USE casim_switches, ONLY: l_mp_cloudnumber, l_mp_rainnumber, l_mp_rain3mom,    &
                          l_mp_icenumber, l_mp_snownumber, l_mp_snow3mom,      &
                          l_mp_graupnumber, l_mp_graup3mom,                    &
                          l_mp_activesolliquid, l_mp_activesolrain,            &
                          l_mp_activeinsolice, l_mp_activesolice,              &
                          l_mp_activeinsolliquid, l_mp_activesolnumber,        &
                          l_mp_activeinsolnumber, n_casim_progs,               &
                          casim_moments_option, l_fix_aerosol,                 &
                          l_tracer_aerosol, l_ukca_aerosol, l_ukca_feeding_in, &
                          l_ukca_feeding_out

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE umprintmgr, ONLY: umprint, ummessage, PrNorm

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CASIM_PRINT_DEPENDENT_SWITCHES'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umprint('Model run includes CASIM Microphysics',                          &
              src='casim_set_dependent_switches_mod')

CALL umprint('Switches set for CASIM based on inputs from namelist run_precip',&
              src='casim_set_dependent_switches_mod')

WRITE(ummessage, '(A, L1)')'l_mp_cloudnumber = ', l_mp_cloudnumber
CALL umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

WRITE(ummessage, '(A, L1)')'l_mp_rainnumber = ', l_mp_rainnumber
CALL umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

WRITE(ummessage, '(A, L1)')'l_mp_rain3mom = ', l_mp_rain3mom
CALL umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

WRITE(ummessage, '(A, L1)')'l_mp_icenumber = ', l_mp_icenumber
CALL umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

WRITE(ummessage, '(A, L1)')'l_mp_snownumber = ', l_mp_snownumber
CALL umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

WRITE(ummessage, '(A, L1)')'l_mp_snow3mom = ', l_mp_snow3mom
CALL umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

WRITE(ummessage, '(A, L1)')'l_mp_graupnumber = ', l_mp_graupnumber
CALL umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

WRITE(ummessage, '(A, L1)')'l_mp_graup3mom = ', l_mp_graup3mom
CALL umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

CALL umprint('- - - - - - - - - - - - - - - - - - - - - - - - - - ',           &
             src='casim_set_dependent_switches_mod')

WRITE(ummessage, '(A, L1)')'l_mp_activesolliquid = ', l_mp_activesolliquid
CALL umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

WRITE(ummessage, '(A, L1)')'l_mp_activesolrain = ', l_mp_activesolrain
CALL umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

WRITE(ummessage, '(A, L1)')'l_mp_activeinsolice = ', l_mp_activeinsolice
CALL umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

WRITE(ummessage, '(A, L1)')'l_mp_activesolice = ', l_mp_activesolice
CALL umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

WRITE(ummessage, '(A, L1)')'l_mp_activeinsolliquid = ', l_mp_activeinsolliquid
CALL umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

WRITE(ummessage, '(A, L1)')'l_mp_activesolnumber = ', l_mp_activesolnumber
CALL umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

WRITE(ummessage, '(A, L1)')'l_mp_activeinsolnumber = ', l_mp_activeinsolnumber
CALL umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

WRITE(ummessage,'(A,I5)')'casim_moments_option = ', casim_moments_option
CALL umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

CALL umprint('- - - - - - - - - - - - - - - - - - - - - - - - - - ',           &
             src='casim_set_dependent_switches_mod')

WRITE(ummessage,'(A,I2)')'Number of additional CASIM prognostics = ',          &
                          n_casim_progs
CALL umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

CALL umprint('- - - - - - - - - - - - - - - - - - - - - - - - - - ',           &
             src='casim_set_dependent_switches_mod')

WRITE(ummessage, '(A, L1)')'l_fix_aerosol = ', l_fix_aerosol
CALL umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

WRITE(ummessage, '(A, L1)')'l_tracer_aerosol = ', l_tracer_aerosol
CALL umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

WRITE(ummessage, '(A, L1)')'l_ukca_aerosol = ', l_ukca_aerosol
CALL umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

WRITE(ummessage, '(A, L1)')'l_ukca_feeding_in = ', l_ukca_feeding_in
CALL umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

WRITE(ummessage, '(A, L1)')'l_ukca_feeding_out = ', l_ukca_feeding_out
CALL umprint(ummessage, level=PrNorm, src='casim_set_dependent_switches_mod')

CALL umprint('- - - - - - end of CASIM switches - - - - - -',                  &
    src='casim_set_dependent_switches_mod')

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE casim_print_dependent_switches

END MODULE casim_set_dependent_switches_mod
