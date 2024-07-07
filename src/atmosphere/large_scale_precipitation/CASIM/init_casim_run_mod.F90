! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!----------------------------------------------------------------------
! Purpose: Sets up all physics settings for a CASIM run
!----------------------------------------------------------------------
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Precipitation (CASIM)
!----------------------------------------------------------------------

MODULE init_casim_run_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='INIT_CASIM_RUN_MOD'

CONTAINS

SUBROUTINE init_casim_run()

! CASIM modules (CASIM repository)
USE casim_parent_mod,             ONLY: casim_parent, parent_um
USE generic_diagnostic_variables, ONLY: casdiags
USE mphys_switches,               ONLY: set_mphys_switches, max_step_length,   &
                                        l_warm, process_level, l_separate_rain,&
                                        nq_l,   nq_r,  nq_i,  nq_s,  nq_g,     &
                                        i_ActiveSolLiquid   => i_am4,          &
                                        i_ActiveSolRain     => i_am5,          &
                                        i_ActiveInsolIce    => i_am7,          &
                                        i_ActiveSolIce      => i_am8,          &
                                        i_ActiveInsolLiquid => i_am9,          &
                                        i_ActiveSolNumber   => i_an11,         &
                                        i_ActiveInSolNumber => i_an12

USE initialize,                   ONLY: mphys_init

USE passive_fields,               ONLY: rhcrit_1d

! Microphysics modules (UM repository)
USE casim_switches,               ONLY: its, ite, jts, jte, kts, kte,          &
                                        ils, ile, jls, jle, kls, kle,          &
                                        irs, ire, jrs, jre, krs, kre,          &
                                        casim_moments_option, n_casim_tracers, &
                                        l_casim_warm_only, l_micro_in_rim,     &
                                        l_mp_CloudNumber, l_mp_RainNumber,     &
                                        l_mp_Rain3mom,    l_mp_IceNumber,      &
                                        l_mp_SnowNumber,  l_mp_Snow3mom,       &
                                        l_mp_GraupNumber, l_mp_Graup3mom,      &
                                        l_mp_ActiveSolLiquid,                  &
                                        l_mp_ActiveSolRain,                    &
                                        l_mp_ActiveInsolIce, l_mp_ActiveSolIce,&
                                        l_mp_ActiveInsolLiquid,                &
                                        l_mp_ActiveSolNumber,                  &
                                        l_mp_ActiveInSolNumber

USE mphys_inputs_mod,             ONLY: i_mcr_iter,  i_mcr_iter_niters,        &
                                        i_mcr_iter_tstep, timestep_mp_in,      &
                                        niters_mp, casim_aerosol_couple_choice,&
                                        casim_aerosol_process_level,           &
                                        casim_aerosol_option,                  &
                                        l_separate_process_rain, l_mcr_qcf2,   &
                                        l_mcr_qgraup

USE cloud_inputs_mod,             ONLY: i_rhcpt, i_cld_vn, rhcrit
USE pc2_constants_mod,            ONLY: i_cld_off, rhcpt_off

USE submodel_mod,                 ONLY: atmos_im
USE model_time_mod,               ONLY: secs_per_stepim

! Dr Hook Modules
USE yomhook,                      ONLY: lhook, dr_hook
USE parkind1,                     ONLY: jprb, jpim

! Error reporting modules
USE errormessagelength_mod,       ONLY: errormessagelength
USE ereport_mod,                  ONLY: ereport

! MPP and bounds modules
USE umPrintMgr,                   ONLY: umprint, umMessage, newline,           &
                                        prstatus_oper
USE um_parvars,                   ONLY: at_extremity
USE um_ParParams,                 ONLY: PNorth, PEast, PSouth, PWest
USE lbc_mod,                      ONLY: rimwidtha
USE rimtypes,                     ONLY: rima_type_norm
USE atm_fields_bounds_mod,        ONLY: tdims
USE atmos_max_sizes,              ONLY: model_levels_max
USE nlsizes_namelist_mod,         ONLY: tr_vars
USE model_domain_mod,             ONLY: model_type, mt_single_column

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER       :: RoutineName='INIT_CASIM_RUN'
CHARACTER(LEN=errormessagelength) :: cmessage           ! Error message

REAL :: timestep      ! Model timestep [s]
REAL :: test_mphys_ts ! Trial microphysics timestep [s]

INTEGER  :: errcode ! Error code
INTEGER, PARAMETER :: warn      = -1 ! Produce UM Warning
INTEGER, PARAMETER :: run_abort = 1  ! Abort UM

INTEGER, PARAMETER :: tracer_aerosol = 1
INTEGER, PARAMETER :: pe0            = 0

INTEGER :: nlevels ! Number of CASIM model levels
INTEGER :: k       ! Loop counter in the vertical (z) dimension

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!---------------------------------------------------------------------
! Start of subroutine
!---------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

errcode = 0

!---------------------------------------------------------------------
! Check that the SCM is not in use as CASIM does not yet work with the
! SCM at this UM version
!---------------------------------------------------------------------
IF ( model_type == mt_single_column ) THEN

  errcode = run_abort
  WRITE(cmessage,FMT='(A)')                                                 &
      'Model configuration set for single column model and   '//newline//   &
      'CASIM microphysics. This option is not yet available. '//newline//   &
      'Switch off CASIM microphysics (l_casim = False in the '//newline//   &
      'run_precip namelist) or set model type to be a        '//newline//   &
      'non-single column configuration.'

  CALL ereport(RoutineName, errcode, cmessage)

END IF ! model_type

!---------------------------------------------------------------------
! Check that the options set by run_cloud namelist
! do not conflict with CASIM. If they do conflict, there is no
! point in doing anything here other than just calling ereport
!---------------------------------------------------------------------
IF ( i_cld_vn /= i_cld_off) THEN

  errcode = run_abort
  WRITE(cmessage,FMT='(A,I0,A,I0,A)')                                       & 
        'Variable i_cld_vn = ', i_cld_vn,   newline //                      &
        'This will not work with CASIM microphysics. '//newline //          &
        'Please set to ', i_cld_off, ' in the run_cloud namelist.'

  CALL ereport(RoutineName, errcode, cmessage)

END IF ! i_cld_vn

IF ( i_rhcpt /= rhcpt_off) THEN

  errcode = run_abort
  WRITE(cmessage,FMT='(A,I0,A,I0,A)')                                       & 
        'Variable i_rhcpt = ', i_rhcpt,   newline //                        &
        'This will not work with CASIM microphysics. '//newline //          &
        'Please set to ', rhcpt_off, ' in the run_cloud namelist.'

  CALL ereport(RoutineName, errcode, cmessage)

END IF ! i_rhcpt

!---------------------------------------------------------------------
! Check that options set by tracers and UKCA do not conflict with
! CASIM. If they do conflict, just call ereport
!---------------------------------------------------------------------

IF ( casim_aerosol_couple_choice == tracer_aerosol ) THEN

  IF ( tr_vars < n_casim_tracers ) THEN

    ! Abort the UM as there are insufficient tracers
    errcode = run_abort

    WRITE(cmessage,FMT='(A,I0,A,I0,A)')                                     & 
    'Error in CASIM configuration: Number of tracers set = ', tr_vars,      &
     newline //                                                             &
    'CASIM is expecting the number of tracers to be ', n_casim_tracers,     &
     newline //                                                             &
    'Please increase the number of tracers to meet this requirement.'

    CALL ereport(RoutineName, errcode, cmessage)

  ELSE IF ( tr_vars > n_casim_tracers ) THEN
    ! Produce a warning to say that the excess tracers will be ignored
    errcode = warn

    WRITE(cmessage,FMT='(A,I0,A,I0,A)')                                     & 
    'CASIM configuration: Number of tracers set = ', tr_vars,  newline//    &
    'CASIM will only use the first ', n_casim_tracers, ' tracers.'          &
     // newline //                                                          &
    'Excess tracers will be ignored by CASIM.'

    CALL ereport(RoutineName, errcode, cmessage)

  END IF ! tr_vars

ELSE IF ( casim_aerosol_couple_choice /= tracer_aerosol                     &
          .AND. tr_vars > 0 ) THEN

  ! Produce a warning to say that tracers will be ignored competeley
  errcode = warn

  WRITE(cmessage,FMT='(A,I0,A,I0,A)')                                       & 
  'Number of tracers set in the CASIM configuration = ', tr_vars,           &
   newline //                                                               &
  'Variable casim_aerosol_couple_choice = ', casim_aerosol_couple_choice,   &
   newline //                                                               &
  'This configuration means that the tracer variables will be ' //newline// &
  'ignored completely by CASIM.'

  CALL ereport(RoutineName, errcode, cmessage)

END IF ! casim_aerosol_couple_choice

!---------------------------------------------------------------------
! Basic setup for CASIM
!---------------------------------------------------------------------

! Work out model timestep
timestep = secs_per_stepim(atmos_im) ! timestep in seconds

! Tell CASIM that its parent model is the UM. This allows for any UM-specific
! operations to take place within CASIM. 
casim_parent = parent_um

! Tell CASIM that it needs to output rain and snowfall rates
! Required on every single UM timestep for JULES
casdiags % l_surface_rain  = .TRUE.
casdiags % l_surface_snow  = .TRUE.
casdiags % l_surface_graup = .TRUE.

!---------------------------------------------------------------------
! Set up microphysics timestep
!---------------------------------------------------------------------
IF (i_mcr_iter == i_mcr_iter_tstep .AND. timestep_mp_in > 0 ) THEN
  ! User has input a microphysics timestep, so we need to use
  ! this to work out the max step length in seconds.

  ! As timestep_mp_in is an INTEGER, it first needs to be made
  ! into REAL before passing to set_mphys_switches.

  test_mphys_ts = REAL(timestep_mp_in)

  ! If that timestep is longer than the model timestep, this needs
  ! to be reset to the model timestep

  IF ( test_mphys_ts > timestep ) THEN
    ! Set max step length to be model timestep
    max_step_length = timestep

    ! Now generate a warning for the user to let them know it has
    ! been reset to the model timestep

    errcode = warn
    WRITE(cmessage,FMT='(A,I0,A,F14.2,A,F14.2,A)')                      &
             'User has requested a ', timestep_mp_in,                   &
             ' second microphysics iterative '                          &
             // newline // ' timestep ' // newline //                   &
             'This is longer than the model timestep (',                &
             timestep, ' seconds) ' // newline //                       &
           '  Running with a single microphysics iteration such that'   &
          // newline // ' max_step_length = ', max_step_length, ' seconds'

    CALL ereport(RoutineName, errcode, cmessage)

  ELSE
    ! Timestep OK
    max_step_length = REAL(test_mphys_ts) 

    WRITE(umMessage,FMT='(A,F14.2,A)')                                  &
           'Microphysics is set up with max_step_length = ',            &
            max_step_length, ' seconds'

    CALL umPrint(umMessage,pe=pe0,level=prstatus_oper)

  END IF ! test_mphys_ts > timestep

ELSE IF (i_mcr_iter == i_mcr_iter_niters .AND. niters_mp > 0 ) THEN

  ! Calculate the max step length as a function of the number of
  ! iterations

  max_step_length = timestep / REAL(niters_mp)

  WRITE(umMessage,FMT='(A,I0,A,F14.2,A)')                               &
           'Microphysics is using ', niters_mp, ' iterations:'          &
         // newline // ' max_step_length = ', max_step_length, ' seconds'
  CALL umPrint(umMessage,pe=pe0,level=prstatus_oper)

ELSE ! test on i_mcr_iter

  ! User has not input a microphysics timestep or it is a wrong
  ! input, so default to fixed value (120 s)

  max_step_length = 120.0

  WRITE(umMessage,FMT='(A,F14.2,A)')                                    &
           'Microphysics is using a default value:'                     &
         // newline // ' max_step_length = ', max_step_length, ' seconds'
  CALL umPrint(umMessage,pe=pe0,level=prstatus_oper)

END IF  !i_mcr_iter

!---------------------------------------------------------------------
! Set up microphysics dimensions
!---------------------------------------------------------------------
! Set up dimensions of the parent model, (the UM) based on tdims
its = tdims%i_start
ite = tdims%i_end
jts = tdims%j_start
jte = tdims%j_end
kts = 1
kte = tdims%k_end

! We don't want to do any microphysics in the lateral boundary
! rim for a nested model...
ils = its
ile = ite
jls = jts
jle = jte
kls = kts
kle = kte - 2

irs = its
ire = ite
jrs = jts
jre = jte
krs = kts
kre = kte - 2

IF (at_extremity(PNorth)) jre = jte - rimwidtha(rima_type_norm)
IF (at_extremity(PEast))  ire = ite - rimwidtha(rima_type_norm)
IF (at_extremity(PSouth)) jrs = jts + rimwidtha(rima_type_norm)
IF (at_extremity(PWest))  irs = its + rimwidtha(rima_type_norm)

IF (l_micro_in_rim) THEN
  jle = jte
  ile = ite
  jls = jts
  ils = its
ELSE
  jle = jre
  ile = ire
  jls = jrs
  ils = irs
END IF ! l_micro_in_rim

!---------------------------------------------------------------------
! Set up microphysics switches
!---------------------------------------------------------------------

IF (l_casim_warm_only) l_warm = .TRUE.

! Set aerosol processing level directly on the CASIM side using the value
! obtained from the run_precip namelist.
process_level = casim_aerosol_process_level

! Set up separate rain processing category for active aerosol on the
! CASIM side.
l_separate_rain = l_separate_process_rain

! Call set mphys_switches with the options passed in from the run_precip
! namelist directly to CASIM.
CALL set_mphys_switches(casim_moments_option, casim_aerosol_option)

!---------------------------------------------------------------------
! Initialise and allocate the space required (CASIM repository)
!---------------------------------------------------------------------
CALL mphys_init( its, ite, jts, jte, kts, kte, ils, ile, jls, jle, kls, kle,   &
                 l_tendency=.FALSE. )

!---------------------------------------------------------------------
! Setup RHCrit
!---------------------------------------------------------------------
! Variable rhcrit_1d will have been allocated with dimensions (kts:kte)
! in the mphys_init routine called above.
nlevels = kte - kts + 1

IF ( nlevels > model_levels_max ) THEN

  ! Something has gone very wrong with the rhcrit input provided.
  ! Abort the model immediately as this will just produce nonsense
  ! in CASIM. 

  errcode = run_abort

  WRITE(cmessage,FMT='(A,I0,A,I0,A)')                                         &
        'While creating RHCrit array for CASIM,'                 //newline//  &
        'Number of CASIM model levels was found to be ', nlevels , newline//  &
        'While model_levels_max is ', model_levels_max,            newline//  &
        'Number of CASIM model levels should be < model_levels_max'

  CALL ereport(RoutineName, errcode, cmessage)

END IF

! Loop over k to generate RHCrit on the CASIM side.
! This must be done as a loop as UM RHCrit will be an array
! of size (model_levels_max), which will usually be bigger
! than nlevels. 
DO k = kts, kte

  IF ( rhcrit(k) > 0.0 ) THEN

    ! Rhcrit input sensible, can pass this value from 
    ! cloud_inputs_mod directly to CASIM 
    rhcrit_1d(k) = rhcrit(k)

  ELSE  ! rhcrit <= 0

    ! Rhcrit value incorrect (and possibly missing data)
    ! Throw an error

    errcode = run_abort
    WRITE(cmessage,FMT='(A,F14.2,A)')                                         &
          'Rhcrit input from run_cloud namelist =', rhcrit(k), newline//      &
          'This is lower than expected and may be missing data' //newline//   &
          'Please correct in the run_cloud namelist.'

    CALL ereport(RoutineName, errcode, cmessage)

  END IF ! rhcrit > 0

END DO ! k

!---------------------------------------------------------------------
! Perform a cross-check of CASIM and UM switches
!---------------------------------------------------------------------

!---------------------------------------------------------------------
! Cross-checking of options relating to moments
!---------------------------------------------------------------------

! Check for liquid or rain being not set. If they are, raise an error

IF (nq_l <= 0 .OR. nq_r <= 0) THEN
  errcode = run_abort

  WRITE(cmessage, FMT = '(A, I0, A, I0, A)')                                &
        'Liquid cloud or rain configuration in UM does not match'//newline//&
        'expected configuration in CASIM. nq_l = ', nq_l, ' and '//newline//&
        'nq_r = ', nq_r, '. Both variables should be above zero '//newline//&
        'Please correct in your configuration set up.'

  CALL ereport(RoutineName, errcode, cmessage)
 
END IF ! nq_l / nq_r <= 0

IF (nq_l > 1 .AND. .NOT. l_mp_CloudNumber ) THEN
  errcode = run_abort
  WRITE(cmessage, FMT = '(A, I0, A)')                                       &
        'Liquid cloud configuration in UM does not match        '//newline//&
        'expected configuration in CASIM. nq_l = ', nq_l, ' and '//newline//&
        'l_mp_CloudNumber is .FALSE. Please correct in    '      //newline//&
        'your configuration set up.'

  CALL ereport(RoutineName, errcode, cmessage)

END IF

IF (nq_r > 1 .AND. .NOT. l_mp_RainNumber ) THEN
  errcode = run_abort
  WRITE(cmessage, FMT = '(A, I0, A)')                                       &
        'Rain configuration in UM does not match expected       '//newline//&
        'configuration in CASIM. nq_r = ', nq_r, ' and '         //newline//&
        'l_mp_RainNumber is .FALSE. Please correct in     '      //newline//&
        'your configuration set up.'

  CALL ereport(RoutineName, errcode, cmessage)

END IF

IF (nq_r > 2 .AND. .NOT. l_mp_Rain3mom ) THEN
  errcode = run_abort
  WRITE(cmessage, FMT = '(A, I0, A)')                                      &
        'Rain configuration in UM does not match expected      '//newline//&
        'configuration in CASIM. nq_r = ', nq_r, ' and '        //newline//&
        'l_mp_Rain3mom is .FALSE. Please correct in your  '     //newline//&
        'configuration set up.'

  CALL ereport(RoutineName, errcode, cmessage)

END IF

IF ( .NOT. l_casim_warm_only ) THEN

  IF ( nq_s <= 0 .AND. nq_i <= 0 ) THEN
    errcode = run_abort

    WRITE(cmessage, FMT = '(A, I0, A, I0, A)')                             &
        'Ice cloud or snow configuration in UM does not match'  //newline//&
        'expected configuration in CASIM. nq_s = ', nq_s, ' and'//newline//&
        'nq_i = ', nq_i, '. Both variables should be above zero'//newline//&
        'Please correct in your configuration set up'

    CALL ereport(RoutineName, errcode, cmessage)

  END IF !nq_s / nq_i < 0

  IF ( nq_i > 0 .AND. .NOT. l_mcr_qcf2 ) THEN
    errcode = run_abort
    WRITE(cmessage, FMT = '(A)')                                           &
          'Logical l_mcr_qcf2 is not set up correctly for      '//newline//&
          'CASIM runs with cloud ice. Please correct in your   '//newline//&
          'configuration set up'

    CALL ereport(RoutineName, errcode, cmessage)

  END IF ! nq_i > 0

  IF ( nq_g > 0 .AND. .NOT. l_mcr_qgraup ) THEN

    errcode = run_abort
    WRITE(cmessage, FMT = '(A)')                                            &
          'Logical l_mcr_qgraup is not set up correctly for     '//newline//&
          'CASIM runs with graupel. Please correct in your      '//newline//&
          'configuration set up'

    CALL ereport(RoutineName, errcode, cmessage)

  END IF ! nq_g > 0

  IF (nq_i > 1 .AND. .NOT. l_mp_IceNumber ) THEN
    errcode = run_abort

    WRITE(cmessage, FMT = '(A, I0, A)')                                      &
          'Ice cloud configuration in UM does not match expected '//newline//&
          'configuration in CASIM. nq_i = ', nq_i, ' and '        //newline//&
          'l_mp_IceNumber is .FALSE. Please correct in your      '//newline//&
          'configuration set up.'

      CALL ereport(RoutineName, errcode, cmessage)

  END IF ! nq_i > 1

  IF (nq_s > 1 .AND. .NOT. l_mp_SnowNumber ) THEN
    errcode = run_abort

    WRITE(cmessage, FMT = '(A, I0, A)')                                      &
          'Snow configuration in UM does not match expected '     //newline//&
          'expected configuration in CASIM. nq_s = ', nq_s, ' and'//newline//&
          'l_mp_SnowNumber is .FALSE. Please correct in your '    //newline//&
          'configuration set up.'

    CALL ereport(RoutineName, errcode, cmessage)

  END IF ! nq_s > 1

  IF (nq_s > 2 .AND. .NOT. l_mp_Snow3mom ) THEN
    errcode = run_abort

    WRITE(cmessage, FMT = '(A, I0, A)')                                      &
          'Snow configuration in UM does not match expected '     //newline//&
          'expected configuration in CASIM. nq_s = ', nq_s, ' and'//newline//&
          'l_mp_Snow3mom is .FALSE. Please correct in your '      //newline//&
          'configuration set up.'


    CALL ereport(RoutineName, errcode, cmessage)

  END IF ! nq_s > 2

  IF (nq_g > 1 .AND. .NOT. l_mp_GraupNumber ) THEN
    errcode = run_abort

    WRITE(cmessage, FMT = '(A, I0, A)')                                      &
          'Graupel configuration in UM does not match expected '  //newline//&
          'expected configuration in CASIM. nq_g = ', nq_g, ' and'//newline//&
          'l_mp_GraupNumber is .FALSE. Please correct in your '   //newline//&
          'configuration set up.'

    CALL ereport(RoutineName, errcode, cmessage)

  END IF ! nq_g > 1

  IF (nq_g > 2 .AND. .NOT. l_mp_Graup3mom ) THEN
    errcode = run_abort

    WRITE(cmessage, FMT = '(A, I0, A)')                                      &
          'Graupel configuration in UM does not match expected '  //newline//&
          'expected configuration in CASIM. nq_g = ', nq_g, ' and'//newline//&
          'l_mp_Graup3mom is .FALSE. Please correct in your '     //newline//&
          'configuration set up.'

    CALL ereport(RoutineName, errcode, cmessage)

  END IF ! nq_g > 2

END IF ! .NOT. l_casim_warm_only

!---------------------------------------------------------------------
! Cross-checking of options relating to aerosols
!---------------------------------------------------------------------

IF ( i_ActiveSolLiquid > 0 .AND. .NOT. l_mp_ActiveSolLiquid ) THEN
  errcode = run_abort

  WRITE(cmessage, FMT = '(A, I0, A)')                                        &
        'Error in CASIM configuration set up: '                   //newline//&
        'i_ActiveSolLiquid = ', i_ActiveSolLiquid, ' and logical '//newline//&
        'l_mp_ActiveSolLiquid is False. Please correct in your  ' //newline//&
        'configuration set up.'

  CALL ereport(RoutineName, errcode, cmessage)

END IF

IF ( i_ActiveSolRain > 0 .AND. .NOT. l_mp_ActiveSolRain ) THEN
  errcode = run_abort

  WRITE(cmessage, FMT = '(A, I0, A)')                                        &
        'Error in CASIM configuration set up: '                   //newline//&
        'i_ActiveSolRain = ', i_ActiveSolRain, ' and logical '    //newline//&
        'l_mp_ActiveSolRain is False. Please correct in your  '   //newline//& 
        'configuration set up.'

  CALL ereport(RoutineName, errcode, cmessage)

END IF

IF ( i_ActiveInsolIce > 0 .AND. .NOT. l_mp_ActiveInsolIce ) THEN
  errcode = run_abort

  WRITE(cmessage, FMT = '(A, I0, A)')                                        &
        'Error in CASIM configuration set up: '                   //newline//&
        'i_ActiveInsolIce = ', i_ActiveInsolIce, ' and logical '  //newline//&
        'l_mp_ActiveInsolIce is False. Please correct in your  '  //newline//&
        'configuration set up.'

  CALL ereport(RoutineName, errcode, cmessage)

END IF
                                        
IF ( i_ActiveSolIce > 0 .AND. .NOT. l_mp_ActiveSolIce ) THEN
  errcode = run_abort

  WRITE(cmessage, FMT = '(A, I0, A)')                                        &
        'Error in CASIM configuration set up: '                   //newline//&
        'i_ActiveSolIce = ', i_ActiveSolIce, 'and logical '       //newline//&
        'l_mp_ActiveSolIce is False. Please correct in your  '    //newline//&
        'configuration set up.'

  CALL ereport(RoutineName, errcode, cmessage)

END IF

IF ( i_ActiveInsolLiquid > 0 .AND. .NOT. l_mp_ActiveInsolLiquid ) THEN
  errcode = run_abort

  WRITE(cmessage, FMT = '(A, I0, A)')                                         &
     'Error in CASIM configuration set up: '                       //newline//&
     'i_ActiveInsolLiquid = ', i_ActiveInsolLiquid, ' and logical '//newline//&
     'l_mp_ActiveInsolLiquid is False. Please correct in your  '   //newline//&
     'configuration set up.'

  CALL ereport(RoutineName, errcode, cmessage)

END IF

IF ( i_ActiveSolNumber > 0 .AND. .NOT. l_mp_ActiveSolNumber ) THEN
  errcode = run_abort

  WRITE(cmessage, FMT = '(A, I0, A)')                                        &
        'Error in CASIM configuration set up: '                   //newline//&
        'i_ActiveSolNumber = ', i_ActiveSolNumber, ' and logical' //newline//&
        'l_mp_ActiveSolNumber is False. Please correct in your  ' //newline//&
        'configuration set up.'

  CALL ereport(RoutineName, errcode, cmessage)

END IF                               
                                        
IF ( i_ActiveInSolNumber > 0 .AND. .NOT. l_mp_ActiveInSolNumber ) THEN
  errcode = run_abort

  WRITE(cmessage, FMT = '(A, I0, A)')                                         &
      'Error in CASIM configuration set up: '                      //newline//&
      'i_ActiveInSolNumber = ', i_ActiveInSolNumber, ' and logical'//newline//&
      'l_mp_ActiveInSolNumber is False. Please correct in your   ' //newline//&
      'configuration set up.'

  CALL ereport(RoutineName, errcode, cmessage)

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE init_casim_run
END MODULE init_casim_run_mod
