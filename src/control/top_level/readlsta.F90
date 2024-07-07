! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description: Read run-time control namelist information configurable
!   in the gui, as required to set up parametrization constants and
!   logical switches needed by physics and dynamics schemes for the
!   Atmosphere model.
!
! Method:  Sequential read of namelists. 
!          Optional print out of namelists
!          Check input namelists
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
! System component: Control Atmos
!

! Subroutine Interface:
SUBROUTINE Readlsta()

! Modules for AC scheme:
USE acdiag_namel_mod,                 ONLY: read_nml_acdiag, print_nlist_acdiag
USE acp_namel_mod,                    ONLY: read_nml_acp, check_nml_acp, &
                                             print_nlist_acp

USE atmos_max_sizes,                  ONLY: model_levels_max

! module for Boundary Layer options (run_BL namelist)
USE bl_option_mod,                    ONLY: &
    check_run_bl,                           &
    print_nlist_run_bl,                     &
    read_nml_run_bl

USE calc_pmsl_inputs_mod,             ONLY: &
    check_nml_run_calc_pmsl,                &
    print_nlist_run_calc_pmsl,              &
    read_nml_run_calc_pmsl

! module for carbon_options namelist
USE carbon_options_mod,               ONLY: &
    check_carbon_options,                   &
    print_nlist_carbon_options,             &
    read_nml_carbon_options

USE casim_set_dependent_switches_mod, ONLY: &
    casim_set_dependent_switches,           &
    casim_print_dependent_switches

USE casim_switches,                   ONLY: n_casim_progs

USE cderived_mod,                     ONLY: h_split_defaults

USE clmchfcg_scenario_mod,            ONLY: &
    read_nml_clmchfcg,                      &
    print_nlist_clmchfcg,                   &
    clmchfcg_rates

! module for RUN_CLOUD namelist
USE cloud_inputs_mod,                 ONLY: &
    i_cld_vn,                               &
    print_nlist_RUN_Cloud,                  &
    check_run_cloud,                        &
    read_nml_run_cloud

USE Control_Max_Sizes

USE coradoca,                         ONLY: &
    coradoca_defaults,                      &
    print_nlist_radfcdia,                   &
    read_nml_radfcdia

! Module for COSP
USE cosp_input_mod,                   ONLY: &
    print_nlist_run_cosp,                   &
    read_nml_run_cosp

! Module for coupling control namelist:
USE coupling_control_mod,             ONLY: &
    check_nml_coupling_control,             &
    read_nml_coupling_control,              &
    print_nlist_coupling_control,           &
    l_oasis

USE cv_param_mod                       ! Access to all variables

USE cv_run_mod                         ! Access to all variables + subroutines

! module for mineral dust scheme options
USE dust_parameters_mod,              ONLY: &
    dust_parameters_check,                  &
    dust_parameters_load,                   &
    dust_size_dist_initialise,              &
    print_nlist_run_dust,                   &
    read_nml_run_dust

USE dynamics_input_mod,               ONLY: &
    check_run_dyn,                          &
    run_dyn,                                &
    read_nml_run_dyn,                       &
    print_nlist_run_dyn,                    &
    numcycles

USE dynamics_testing_mod

! module for RUN_ELECTRIC namelist
USE electric_inputs_mod,              ONLY: &
    check_run_electric,                     &
    print_nlist_run_electric,               &
    read_nml_run_electric

! module for RUN_Eng_Corr namelist
USE eng_corr_inputs_mod,              ONLY: &
    print_nlist_run_eng_corr,               &
    read_nml_run_eng_corr,                  &
    check_run_eng_corr

USE ereport_mod,                      ONLY: ereport

USE errormessagelength_mod,           ONLY: errormessagelength

USE file_manager,                     ONLY: get_file_unit_by_id

USE free_tracers_inputs_mod,          ONLY: read_nml_run_free_tracers,  &
                                            print_nlist_run_free_tracers

! Module for GWD scheme:
USE g_wave_input_mod

! Module for general physics options
USE gen_phys_inputs_mod,              ONLY: &
    read_nml_gen_phys_inputs,               &
    gen_phys_inputs,                        &
    print_nlist_gen_phys_inputs

! module for GLOMAP_CLIM options
USE glomap_clim_option_mod,           ONLY: &
    print_nlist_run_glomap_clim,            &
    read_nml_run_glomap_clim,               &
    check_run_glomap_clim

! Module for IAU scheme:
USE IAU_mod,                          ONLY: &
    l_iau,                                  &
    read_nml_iau_nl

USE idealise_run_mod,                 ONLY: &
    read_nml_idealised,                     &
    check_nlist_idealised,                  &
    print_nlist_idealised

USE jules_surface_mod,                ONLY: &
    l_flake_model,                          &
    Limit_ObukhovL,                         &
    cor_mo_iter

USE jules_sea_seaice_mod,             ONLY: &
    l_ssice_albedo,                         &
    alpham,                                 &
    ssalpham,                               &
    alphac,                                 &
    ssalphac,                               &
    dtice,                                  &
    ssdtice

USE lam_config_inputs_mod,            ONLY: &
    check_nml_lam_config,                   &
    print_nlist_lam_config,                 &
    read_nml_lam_config

! Module for lbc options
USE lbc_read_data_mod,                ONLY: &
    lbc_options,                            &
    read_nml_lbc_options,                   &
    print_nlist_lbc_options

USE lw_rad_input_mod

USE model_domain_mod,                 ONLY: &
    check_nml_model_domain,                 &
    read_nml_model_domain,                  &
    print_nlist_model_domain

! module for RUN_PRECIP namelist
USE mphys_inputs_mod,                 ONLY: &
    check_run_precip,                       &
    print_nlist_run_precip,                 &
    read_nml_run_precip,                    &
    l_casim

! module for RUN_Murk namelist
USE murk_inputs_mod,                  ONLY: &
    print_nlist_run_murk,                   &
    read_nml_run_murk

USE nlsizes_namelist_mod,             ONLY: land_field

USE nlstcall_mod,                     ONLY: l_fastrun

! Module for Nudging scheme:
USE nudging_input_mod,                ONLY: &
    print_nlist_run_nudging,                &
    read_nml_run_nudging

USE ozone_inputs_mod,                 ONLY: &
    print_nlist_run_ozone,                  &
    check_run_ozone,                        &
    read_nml_run_ozone

USE parkind1,                         ONLY: &
    jprb,                                   &
    jpim

USE pc2_constants_mod,                ONLY: i_cld_pc2

USE planet_constants_mod,             ONLY: &
    planet_constants,                       &
    set_planet_constants,                   &
    read_nml_planet_constants,              &
    print_nlist_planet_constants

USE problem_mod,                      ONLY: standard

USE rad_input_mod,                    ONLY: &
    check_run_radiation,                    &
    print_nlist_run_radiation,              &
    read_nml_run_radiation,                 &
    l_radiation

! module for RUN_RIVERS namelist
USE river_inputs_mod,                 ONLY: &
    print_nlist_run_rivers,                 &
    read_nml_run_rivers

! module for aerosol emissions options
USE run_aerosol_mod,                  ONLY: &
    print_nlist_run_aerosol,                &
    read_nml_run_aerosol

USE science_fixes_mod

USE set_rad_steps_mod,                ONLY: set_a_radstep

USE sl_input_mod,                     ONLY: &
    print_nlist_run_sl,                     &
    check_run_sl,                           &
    read_nml_run_sl

! Module for stochastic physics UM section 35 : SKEB2
USE stochastic_physics_run_mod,       ONLY: &
    print_nlist_run_stochastic,             &
    check_run_stochastic,                   &
    read_nml_run_stochastic

! Module for reading JULES namelists:
USE surf_couple_read_namelists_mod,   ONLY: surf_couple_read_namelists

USE sw_rad_input_mod

! module for FV-TRACK (Cyclone tracking)
USE track_mod,                        ONLY: &
    print_nlist_run_track,                  &
    read_nml_run_track

! module for Segments options (tuning_segments namelist)
USE tuning_segments_mod,              ONLY: &
    check_tuning_segments,                  &
    print_nlist_tuning_segments,            &
    read_nml_tuning_segments

USE turb_diff_mod,                    ONLY: &
    run_diffusion,                          &
    check_run_diffusion,                    &
    print_nlist_run_diffusion,              &
    read_nml_run_diffusion

USE ukca_init_mod,                    ONLY: ukca_init

! module for UKCA options
USE ukca_option_mod,                  ONLY: &
    print_nlist_RUN_ukca,                   &
    read_nml_run_ukca,                      &
    l_ukca

USE UM_ParCore,                       ONLY: mype

USE umPrintMgr,                       ONLY: &
    PrintStatus,                            &
    PrStatus_Normal,                        &
    umPrint,                                &
    umMessage

USE yomhook,                          ONLY: &
    lhook,                                  &
    dr_hook

! module for EasyAerosol options
USE easyaerosol_option_mod,           ONLY: &
    print_nlist_easyaerosol,                &
    read_nml_easyaerosol

USE cv_set_dependent_switches_mod, ONLY: cv_set_dependent_switches

USE qsat_mod,                         ONLY: init_qsat_switches

IMPLICIT NONE


! Local parameters:
CHARACTER(LEN=*) :: RoutineName
PARAMETER (   RoutineName='READLSTA')

! Local scalars:
INTEGER ::                                                                    &
 ErrorStatus      ! Return code : 0 Normal Exit : >0 Error
CHARACTER(LEN=errormessagelength) ::                                          &
 CMessage         ! Error message if Errorstatus >0

INTEGER :: level ,i, j, k, jj
INTEGER :: atmoscntl_unit  ! unit no. for ATMOSCNTL file
INTEGER :: shared_unit     ! unit no. for SHARED file

LOGICAL :: l_print_namelist = .FALSE.  ! print out namelist entries.

CHARACTER(LEN=50000) :: lineBuffer

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! -------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

ErrorStatus = 0

! Retrieve the unit number of the main namelist files
atmoscntl_unit = get_file_unit_by_id("atmoscntl", handler="fortran")
shared_unit    = get_file_unit_by_id("shared",    handler="fortran")

! height values to split model levels into l/m/h cloud 
! not an input at present.
CALL h_split_defaults()

IF (PrintStatus >= PrStatus_Normal .AND. mype == 0) THEN
  WRITE(lineBuffer,'(A)') '******************** ' // RoutineName //           &
                          ': Atmosphere run-time constants *******************'
  CALL umPrint(lineBuffer,src=RoutineName)
  l_print_namelist = .TRUE.
END IF

! -----------------------------------------
! read in scientific fixes namelist
! controls what fixes are applied
! each of these fixes are anticipated
! to become the default code in the future
! -----------------------------------------

CALL read_nml_temp_fixes(shared_unit)
IF (l_print_namelist) CALL print_nlist_temp_fixes()
CALL warn_temp_fixes()

! Carbon options
CALL read_nml_carbon_options(shared_unit)
IF (l_print_namelist) CALL print_nlist_carbon_options()
CALL check_carbon_options()

! Coupling options
CALL read_nml_coupling_control(shared_unit)
IF (l_print_namelist) CALL print_nlist_coupling_control()
CALL check_nml_coupling_control()

! Model domain
CALL read_nml_model_domain(shared_unit)
IF (l_print_namelist) CALL print_nlist_model_domain()
CALL check_nml_model_domain()

! Planet constants
CALL read_nml_planet_constants(shared_unit)
CALL set_planet_constants()
! Print needs to be after set otherwise g etc look unset for Earth simulations.
IF (l_print_namelist) CALL print_nlist_planet_constants()

! Mineral dust modelling
! Initialise sizes in case not set in namelist
CALL dust_size_dist_initialise
CALL read_nml_run_dust(shared_unit)
IF (l_print_namelist) CALL print_nlist_RUN_Dust()
CALL dust_parameters_load
CALL dust_parameters_check

! FV-track for cyclone tracking diagnostics
CALL read_nml_run_track(atmoscntl_unit)
IF (l_print_namelist) CALL print_nlist_run_track()

! GLOMAP_CLIM Sub-model
CALL read_nml_run_glomap_clim(shared_unit)
IF (l_print_namelist) CALL print_nlist_RUN_Glomap_Clim()

! UKCA Sub-model
CALL read_nml_run_ukca(shared_unit)
IF (l_print_namelist) CALL print_nlist_RUN_UKCA() 

! Gravity wave drag physics
CALL read_nml_run_gwd(shared_unit)
IF (l_print_namelist) CALL print_nlist_RUN_GWD()

! Murk aerosol physics
! This namelist must be read before the run_precip namelist
! and check_run_precip is called as precip error checking 
! depends on run_murk switches (see check_run_precip for details) 
CALL read_nml_run_murk(shared_unit)
IF (l_print_namelist) CALL print_nlist_RUN_Murk()

! Convection physics
CALL read_nml_run_convection(shared_unit)
IF (l_print_namelist) CALL print_nlist_RUN_Convection() 
CALL check_run_convection()

! Boundary layer physics
CALL read_nml_run_bl(shared_unit)
IF (l_print_namelist) CALL print_nlist_RUN_BL()
CALL check_run_bl()

! River routing
CALL read_nml_run_rivers(shared_unit)
IF (l_print_namelist) CALL print_nlist_RUN_RIVERS()

! Large scale precipitation physics
CALL read_nml_run_precip(shared_unit)
IF (l_print_namelist) CALL print_nlist_RUN_Precip()
CALL check_run_precip()

! Radiation physics
CALL read_nml_run_radiation(shared_unit)
IF (l_print_namelist) CALL print_nlist_RUN_Radiation()
CALL check_run_radiation()

! Large scale cloud physics
CALL read_nml_run_cloud(shared_unit)
IF (l_print_namelist) CALL print_nlist_RUN_Cloud()
CALL check_run_cloud()

! Aerosol Modelling
CALL read_nml_run_aerosol(shared_unit)
IF (l_print_namelist) CALL print_nlist_RUN_Aerosol()
! run_aerosol_check called in readsize because it need model_levels

! LAM configuration
CALL read_nml_lam_config(shared_unit)
IF (l_print_namelist) CALL print_nlist_lam_config()
CALL check_nml_lam_config()

! Ozone
CALL read_nml_run_ozone(shared_unit)
IF (l_print_namelist) CALL print_nlist_run_ozone()
CALL check_run_ozone()

! Free tracers
CALL read_nml_run_free_tracers(shared_unit)
IF (l_print_namelist) CALL print_nlist_run_free_tracers()

! Energy correction physics
CALL read_nml_run_eng_corr(shared_unit)
IF (l_print_namelist) CALL print_nlist_RUN_Eng_Corr()
CALL check_run_eng_corr()

! Calc pmsl
CALL read_nml_run_calc_pmsl(atmoscntl_unit)
IF (l_print_namelist) CALL print_nlist_run_calc_pmsl()
CALL check_nml_run_calc_pmsl()

! General physics
CALL read_nml_gen_phys_inputs(shared_unit)
IF (l_print_namelist) CALL print_nlist_gen_phys_inputs()

!Set new qsat switches before any calls to qsat
CALL init_qsat_switches()

! LBC options
CALL read_nml_lbc_options(atmoscntl_unit)
IF (l_print_namelist) CALL print_nlist_lbc_options()

! UM nudging
CALL read_nml_run_nudging(atmoscntl_unit)
IF (l_print_namelist) CALL print_nlist_RUN_Nudging()

! Generalised integration and GCR dynamics
CALL read_nml_run_dyn(shared_unit)
IF (l_print_namelist) CALL print_nlist_RUN_Dyn()
CALL check_run_dyn()

! Generalised integration and GCR dynamics
CALL read_nml_run_dyntest(shared_unit)
IF (l_print_namelist) CALL print_nlist_RUN_Dyntest()
CALL check_run_dyntest()

! Idealised Model if required
! - problem_number is set in run_dyntest
IF (problem_number /= standard) THEN
  CALL read_nml_idealised()
  IF (l_print_namelist) CALL print_nlist_idealised()
  CALL check_nlist_idealised()
END IF

! Semi-Lagrangian advection dynamics
CALL read_nml_run_sl(atmoscntl_unit)
IF (l_print_namelist) CALL print_nlist_RUN_SL()
CALL check_run_sl()

! Diffusion, divergence damping and filtering dynamics
CALL read_nml_run_diffusion(atmoscntl_unit)
IF (l_print_namelist) CALL print_nlist_RUN_Diffusion()
CALL check_run_diffusion()

! Call to COSP
CALL read_nml_run_cosp(atmoscntl_unit)
IF (l_print_namelist) CALL print_nlist_RUN_COSP()

! Diagnostic double call to radiation
CALL read_nml_radfcdia(atmoscntl_unit)
IF (l_print_namelist) CALL print_nlist_RADFCDIA()

! EasyAerosol
CALL read_nml_easyaerosol(shared_unit)
IF (l_print_namelist) CALL print_nlist_easyaerosol()

!-------------------------------------------------------
! Set other dependent convective switches valid for whole run
CALL cv_set_dependent_switches

! Set model timesteps per radiation timestep (valid for whole run), and
! the radiation timestep length
IF (l_radiation)  CALL set_a_radstep()

!-------------------------------------------------------
! Set up switches based on CASIM, valid for whole run.
! Only call this routine if CASIM is switched on.
IF ( l_casim ) THEN
  CALL casim_set_dependent_switches
  IF (l_print_namelist) CALL casim_print_dependent_switches
ELSE
! Explicitly set the number of extra CASIM prognostics
! to zero to avoid affecting moisture array size and
! other variables
  n_casim_progs = 0
END IF

!-------------------------------------------------------
! Check that if L_GLOMAP_CLIM_RADAER is selected certain
! other switches are not set so that they conflict.
!-------------------------------------------------------
CALL check_run_glomap_clim()


IF (i_cld_vn == i_cld_pc2) THEN
  ! So that PC2 does NOT affect Section 5 diagnostics
  l_dcpl_cld4pc2=.TRUE.
END IF

! Options for the shortwave radiation
CALL sw_input

! Options for the longwave radiation
CALL lw_input

CALL read_nml_clmchfcg(atmoscntl_unit)
CALL clmchfcg_rates()
IF (l_print_namelist) CALL print_nlist_clmchfcg()

CALL coradoca_defaults

! ACP, ACDIAG
CALL read_nml_acp(atmoscntl_unit)
IF (l_print_namelist) CALL print_nlist_acp()
CALL check_nml_acp()
CALL read_nml_acdiag(atmoscntl_unit)
IF (l_print_namelist) CALL print_nlist_acdiag()

! Read the JULES namelists
CALL surf_couple_read_namelists("UM",shared_unit,atmoscntl_unit)

! Stochastic physics
CALL read_nml_run_stochastic(shared_unit)
IF (l_print_namelist) CALL print_nlist_RUN_Stochastic()
CALL check_run_stochastic()

! Electric physics
CALL read_nml_run_electric(shared_unit)
IF (l_print_namelist) CALL print_nlist_run_electric()
CALL check_run_electric()

!-----------------------------------------------------
! Below is where we provide a location for inter namelist checks
!-----------------------------------------------------

! If UKCA is active, check UKCA logicals are consistent and
! set internal UKCA values from UKCA namelists
! Needs to be after CLASSIC aerosol namelist is read as
! well as after UKCA and after gen_phys_inputs
IF (l_ukca) THEN
  CALL ukca_init()
END IF

! In the coupled model case, set alpham/c to ssalpham/c and dtice to ssdtice
! if l_ssice_albedo == T. Note, ssalpham/c, ssdtice accessed via rad_input_mod
! Do this after JULES namelist reads, as l_ssice_albedo is in jules_sea_seaice
IF (l_oasis .AND. l_ssice_albedo) THEN
  alpham = ssalpham
  alphac = ssalphac
  dtice  = ssdtice
END IF

!--------------------------------------------------------

! Read IAU namelist:
CALL read_nml_iau_nl(atmoscntl_unit)

! Turn off IAU if in FASTRUN mode:
IF (l_fastrun) THEN
  l_iau = .FALSE.
END IF

! Segments
CALL read_nml_tuning_segments(atmoscntl_unit)
IF (l_print_namelist) CALL print_nlist_tuning_segments()
CALL check_tuning_segments()

! Check error condition
IF (ErrorStatus >  0) THEN

  CALL Ereport(RoutineName,ErrorStatus,Cmessage)
END IF

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE Readlsta
