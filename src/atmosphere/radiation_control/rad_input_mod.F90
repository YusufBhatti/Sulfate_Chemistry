! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Input control for radiation.

! Description:
!   Module containing input switches/settings as used by the radiation code.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control

! Method:
!   Switches are initialised to false and read in from the
!   namelist. The module may then be used directly where the switches
!   are needed within the radiation code.

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3

MODULE rad_input_mod

USE missing_data_mod, ONLY: imdi, rmdi
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! ----------------
! Control options
! ----------------

LOGICAL :: l_radiation = .FALSE.    !  F: Turns off radiation code

LOGICAL :: l_use_dust = .FALSE.     !  Use mineral dust in rad calculations
LOGICAL :: l_use_biogenic = .FALSE. !  Use biogenic aerosol in radiation code

! Use SO4 aerosol from sulphur cycle for direct/indirect effect
! in radiation, the latter for both SW and LW.
LOGICAL :: l_use_sulpc_direct = .FALSE.
LOGICAL :: l_use_sulpc_indirect_sw = .FALSE.
LOGICAL :: l_use_sulpc_indirect_lw = .FALSE.

! Indirect radiative effect of sea-salt
LOGICAL :: l_use_seasalt_indirect = .FALSE.

! Direct radiative effect of sea-salt
LOGICAL :: l_use_seasalt_direct = .FALSE.

LOGICAL :: l_use_soot_direct = .FALSE.   ! direct radiative effects of soot
LOGICAL :: l_use_soot_indirect = .FALSE. ! indirect effects of soot

! Use biomass aerosol for direct/indirect effect in radiation.
LOGICAL :: l_use_bmass_direct = .FALSE.
LOGICAL :: l_use_bmass_indirect = .FALSE.

! Use fossil-fuel organic carbon aerosol for direct/indirect
! effect in radiation
LOGICAL :: l_use_ocff_direct = .FALSE.
LOGICAL :: l_use_ocff_indirect = .FALSE.

! Use ammonium nitrate aerosol for direct/indirect effect in radiation
LOGICAL :: l_use_nitrate_direct = .FALSE.
LOGICAL :: l_use_nitrate_indirect = .FALSE.

! Use the same number of cloud droplets for 1st and 2nd indirect effects
LOGICAL :: l_consistent_cdnc = .FALSE.

! Use aerosol climatologies in radiation instead of prognostic variables
! Set on a species by species basis
LOGICAL :: l_use_arclbiom = .FALSE. ! biomass burning aerosol
LOGICAL :: l_use_arclblck = .FALSE. ! black carbon
LOGICAL :: l_use_arclsslt = .FALSE. ! sea salt
LOGICAL :: l_use_arclsulp = .FALSE. ! sulpahtes
LOGICAL :: l_use_arcldust = .FALSE. ! mineral dust
LOGICAL :: l_use_arclocff = .FALSE. ! organic carbon (fossil fuel)
LOGICAL :: l_use_arcldlta = .FALSE. ! delta aerosol

! Use droplet number from n_drop_pot array:
LOGICAL :: l_use_ndrop = .FALSE.

! Use Liu 2008 spectral broadening in calculation of effective radius
LOGICAL :: l_use_liu_spec = .FALSE.
! Parameters for tuning the beta function in the Liu (2008) scheme
REAL :: aparam = rmdi
REAL :: bparam = rmdi

LOGICAL :: L_rad_deg = .FALSE.  ! controls the use of spatial degradation
!                                 of radiation calc.

! Logicals for different radiation packages
LOGICAL :: l_forcing     = .FALSE. ! Calculate radiative forcings
LOGICAL :: l_radiance    = .FALSE. ! Calculate radiances
LOGICAL :: l_timestep    = .FALSE. ! Use new timestepping scheme
LOGICAL :: l_rad_perturb = .FALSE. ! Use the perturbation version of
                                   ! the radiative time-stepping
! Control integer for different radiation packages used by check_run_radiation()
! routine to set l_forcing, l_radiance, l_timestep and l_rad_perturb
INTEGER :: i_rad_extra_call = imdi

! Use a solar zenith angle correction based on optical depth
LOGICAL :: l_rad_szacor = .FALSE.
!                       Needed in glue_rad-rad_ctl3c.
!                       Switch for the solar zenith angle correction to surface
!                       fluxes using the change in optical depth.

! Scale the condensed water content to simulate
! inhomogeneous clouds
LOGICAL :: l_inhom_cloud = .FALSE.

! Orography correction to SW radiation
LOGICAL :: l_use_orog_corr  = .FALSE.  !  Find gradients from mean orog
LOGICAL :: l_use_grad_corr  = .FALSE.  !  Use ancillary X & Y gradients
! Correction for skyview factor in LW and direct SW
LOGICAL :: l_use_skyview    = .FALSE.
LOGICAL :: l_orog_unfilt    = .FALSE.  !  Use unfiltered ancillary orog
! Control integer used by check_radiaiton to set the orography correction
! and skyview logicals
INTEGER :: i_rad_topography = imdi

! ----------------------------

INTEGER :: h_swbands    = imdi   ! Number of shortwave radiation bands
INTEGER :: h_lwbands    = imdi   ! Number of longwave radiation bands

! Number of advection steps per prognostic/diagnostic SW and LW step.
!'Prognostic' and 'Diagnostic' refer to the frequency of the calls
! to radiation code in the Unified Model.
! In the case of time stepping prognostic and diagnostic refer to the
! slow and fast radiative timestep respectively. In the case of radiative
! forcing they refer to the prognostic and diagnostic calls to radiation.

INTEGER :: a_sw_radstep_diag = imdi
! Number of advection steps per 'fast' SW step (3C)
INTEGER :: a_lw_radstep_diag = imdi
! Number of advection steps per 'fast' LW step (3C)
INTEGER :: a_sw_radstep_prog = imdi
! Number of advection steps per 'slow' LW step (3C)
INTEGER :: a_lw_radstep_prog = imdi
! Number of advection steps per 'slow' LW step (3C)

INTEGER :: i_ozone_int = imdi ! Option for interpolation of ozone


! The following three switches REMOVED from run_radiation NL (ROSE project)
! They are set in check_run_radiation, dependent on settings of cusack_aero and
! cusack_aero_hgt, which have been added to the NL

! True if climatological aerosol is included.
LOGICAL :: L_climat_aerosol = .FALSE.

! True to use real boundary layer heights to specify the boundary
! layer aerosol.
LOGICAL :: L_clim_aero_hgt = .FALSE.

! Flag to use HadGEM1 setting for climatological aerosols
LOGICAL :: L_HadGEM1_Clim_Aero = .FALSE.

! These two switches ADDED to NL (ROSE project)
INTEGER :: cusack_aero     = imdi
INTEGER :: cusack_aero_hgt = imdi

LOGICAL :: lrad_ccrad = .FALSE.
!             Allows access to ccrad code and the logicals
!             lrad_ovrlap and lrad_ccw_scav

! Convert zonal mean ozone to field
LOGICAL :: lexpand_ozone

! convert zonal mean tpps ozone to field
LOGICAL :: lexpand_tpps_ozone

! Tropopause-based Ozone Scheme
LOGICAL :: l_use_tpps_ozone = .FALSE.  !  Use TPPS ozone scheme

!-----------------------------------------------------------
! run_radiation namelists
! ----------------------------------------------------------

LOGICAL :: l_sec_var  = .FALSE.     ! true if using time varying astronomy

LOGICAL :: l_rad_ovrlap   = .FALSE.
!                           Requires l_ccrad=.TRUE.
!                           Allows Convective and LS Cloud to overlap
!                           for radiative impacts.
!                           (Experimental, defaulted to FALSE,
!                           requires a hand-edit to change)
!             (THIS IS EXPERIMENTAL AND USED FOR DEVELOPMENT ONLY).
!             Current convective/large-scale cloud fractions in the
!             radiation scheme are mutally exclusive. This assumes CCA
!             and the large-scale to overlap with CCA taking dominance.
!             I.E. Large-scale cloud fraction must exceed the convective
!             cloud fraction before having any presence.

LOGICAL :: l_rad_ccw_scav = .FALSE.
!                           Requires l_ccrad=.TRUE. .AND. l_rad_ovrlap=.TRUE.
!                           Allows Convective Cloud Water (CCW) to
!                           compensate for LS Cloud water in overlapping
!                           LS/CCA fractions.
!                           (Experimental, defaulted to FALSE, requires a
!                           hand-edit to change)
!             (THIS IS EXPERIMENTAL AND USED FOR DEVELOPMENT ONLY)
!             Allowing the CCA to negate large-scale fractions of lower
!             values means that the large-scale cloud water in the
!             overlapping fraction is lost. This switch will scavenge
!             the large-scale cloud water from the overlapping fraction
!             and combine it with the convective cloud water to
!             conpensate.

LOGICAL :: l_rad_use_clim_volc=.FALSE.
!                           If .TRUE. use climatological volcanic
!                           eruption code in climatological aerosol
!                           code


LOGICAL :: l_rad_snow_emis = .FALSE.
!          Switch to adjust the emissivity in radiation for snow cover.
!          This should eventually be moved to the surface scheme, but
!          JULES cannot currently cope with distinct emissivities
!          for snow-covered surfaces, so the switch currently acts
!          only in radiation and logically belongs here for the present.
LOGICAL :: l_t_land_nosnow = .FALSE.
!          Switch for emissivity of snow used in averaging the
!          surface temperature. Setting this switch to .TRUE. is
!          deprecated and it is included only for historical reasons.
LOGICAL :: l_quad_t_coast = .FALSE.
!          Switch for quadratic averaging of the surface temperature at
!          coastal points. .FALSE. is deprecated.
LOGICAL :: l_t_rad_solid = .FALSE.
!          Switch to use common soid temperature at coastal points with
!          sea-ice. .TRUE. is deprecated.

LOGICAL :: l_t_bdy_surf = .FALSE.
!          Take the temperature of the air just above the surface as
!          the temperature at the surface.

! ------------------------------------------

! number of components of clouds
INTEGER,PARAMETER:: npd_cloud_component=4

INTEGER :: aero_bl_levels = imdi
!                          Common number of layers taken to be
!                          occupied by the boundary-layer
!                          aerosol if the boundary layer
!                          depth is not used to determine the
!                          number separately at each grid-point
!                          In previous versions of the code,
!                          this was taken to be BL_LEVELS

INTEGER :: clim_rad_volc_eruption_year = imdi  ! Climatological volcano
!                                                eruption year

INTEGER :: clim_rad_volc_eruption_month = imdi    ! Climatological volcano
!                                                eruption month

INTEGER :: rad_mcica_sampling = imdi   ! Version of McICA used (was 1)
!             Needed in open_cloud_gen. Selects the version of McICA
!             used to sample the generated cloud:
!                               0 = full sampling
!                               1 = single sampling
!                               2 = optimal sampling

! --------------------------------

REAL    :: rad_mcica_sigma = rmdi  ! Normalised cloud condensate standard
!                                    deviation for the cloud generator.
!                                    Needed in open_cloud_gen.

REAL    :: two_d_fsd_factor = rmdi ! Ratio between fsd in 2D (required by UM) 
                                   ! and 1D (as parametrized).

LOGICAL :: l_fsd_eff_res = .FALSE. ! Use an effective resolution rather than
                                   ! the actual grid-length in the FSD
                                   ! parametrization

INTEGER :: i_cloud_representation = imdi
INTEGER :: i_inhom = imdi
INTEGER :: i_overlap = imdi
INTEGER :: i_fsd = imdi
INTEGER :: i_cloud_representation_2 = imdi
INTEGER :: i_inhom_2 = imdi
INTEGER :: i_overlap_2 = imdi
INTEGER :: i_fsd_2 = imdi

! Mass Mixing Ratios (MMR) of minor Gases
REAL :: co2_mmr    = rmdi ! CO2 concentration (if constant)
REAL :: n2ommr     = rmdi ! N2O mmr
REAL :: ch4mmr     = rmdi ! CH4 mmr
REAL :: c11mmr     = rmdi ! CFC11 mmr
REAL :: c12mmr     = rmdi ! CFC12 mmr
REAL :: o2mmr      = rmdi ! O2 mmr
REAL :: so2mmr     = rmdi ! SO2 mmr
REAL :: c113mmr    = rmdi ! CFC113 mmr
REAL :: c114mmr    = rmdi ! CFC114 mmr
REAL :: hcfc22mmr  = rmdi ! HCFC22 mmr
REAL :: hfc125mmr  = rmdi ! HFC125 mmr
REAL :: hfc134ammr = rmdi ! HFC134A mmr

! Scaling factors to simulate inhomogeneous cloud.
REAL    :: inhom_cloud_sw(npd_cloud_component) = rmdi
REAL    :: inhom_cloud_lw(npd_cloud_component) = rmdi


! Decorrelation pressure scale for large scale cloud
REAL    :: dp_corr_strat = rmdi

! Decorrelation pressure scale for convective cloud
REAL    :: dp_corr_conv  = rmdi


REAL    :: clim_rad_volc_eruption_weight = rmdi
! Eruption weighting factor for idealised volcanic aerosol.
! 1.0 is an average 20th century tropical explosive eruption.

REAL    :: aeroscl_csk_clim(5) = (/ rmdi, rmdi, rmdi, rmdi, rmdi /)
! Scalings for aerosols in Cusack's climatology

! Number of radiation prognostic/diagnostic timesteps per day
INTEGER :: i_sw_radstep_perday_prog = imdi
INTEGER :: i_lw_radstep_perday_prog = imdi
INTEGER :: i_sw_radstep_perday_diag = imdi
INTEGER :: i_lw_radstep_perday_diag = imdi

! Ozone tracer as input to radiation scheme
LOGICAL :: l_use_cariolle   = .FALSE.
LOGICAL :: l_use_ozoneinrad = .FALSE.

! Use abundances from Burrows & Sharp, ApJ, 1999 (for hot Jupiters)
LOGICAL :: l_BS1999_abundances = .FALSE.

! Calculate layer masses using the hydrostatic approximation
LOGICAL :: l_hydrostatic_mass = .FALSE.

! Calculate layer heat capacities including the moisture
LOGICAL :: l_moist_heat_capacity = .FALSE.

! Create an extra top layer for radiation
LOGICAL :: l_extra_top = .FALSE.

NAMELIST/RUN_Radiation/ &
       cusack_aero, cusack_aero_hgt, aeroscl_csk_clim, co2_mmr,        &
       l_sec_var, inhom_cloud_sw, inhom_cloud_lw, dp_corr_strat,       &
       rad_mcica_sampling, rad_mcica_sigma, two_d_fsd_factor,          &
       l_fsd_eff_res, dp_corr_conv, aero_bl_levels,                    &
       l_rad_use_clim_volc, clim_rad_volc_eruption_year,               &
       clim_rad_volc_eruption_month, clim_rad_volc_eruption_weight,    &
       n2ommr, ch4mmr, so2mmr, c11mmr, c12mmr,                         &
       o2mmr, c113mmr, c114mmr, hcfc22mmr, hfc125mmr, hfc134ammr,      &
       i_cloud_representation, i_cloud_representation_2,               &
       i_inhom, i_inhom_2, i_overlap, i_overlap_2, i_fsd, i_fsd_2,     &
       l_rad_snow_emis, l_t_land_nosnow, l_quad_t_coast, l_t_rad_solid,&
       l_t_bdy_surf,                                                   &
       h_swbands, h_lwbands, i_rad_extra_call, i_rad_topography,       &
       l_radiation, l_rad_deg,                                         &
       l_rad_szacor, i_sw_radstep_perday_prog,                         &
       i_lw_radstep_perday_prog,  i_sw_radstep_perday_diag,            &
       i_lw_radstep_perday_diag, l_use_cariolle, l_use_ozoneinrad,     &
       i_ozone_int, l_consistent_cdnc, l_use_sulpc_direct,             &
       l_use_sulpc_indirect_sw, l_use_sulpc_indirect_lw,               &
       l_use_seasalt_direct, l_use_seasalt_indirect, l_use_biogenic,   &
       l_use_nitrate_direct, l_use_nitrate_indirect, l_use_soot_direct,&
       l_use_soot_indirect, l_use_bmass_direct, l_use_bmass_indirect,  &
       l_use_ocff_direct, l_use_ocff_indirect, l_use_dust,             &
       l_use_arclbiom, l_use_arclblck, l_use_arcldlta, l_use_arcldust, &
       l_use_arclocff, l_use_arclsslt, l_use_arclsulp,                 &
       l_BS1999_abundances, l_hydrostatic_mass, l_moist_heat_capacity, &
       l_extra_top, l_use_liu_spec, aparam, bparam

! Logical variables to control whether different cca type progonostics
! are required by the cloud generator
LOGICAL :: l_cca_dp_prog = .FALSE.
LOGICAL :: l_cca_md_prog = .FALSE.
LOGICAL :: l_cca_sh_prog = .FALSE.

! The number of calls to SW/LW radiation
INTEGER :: n_swcall = imdi
INTEGER :: n_lwcall = imdi

! Parameters for the values that i_rad_extra_call can have:
INTEGER, PARAMETER :: ip_single_call           = 0
INTEGER, PARAMETER :: ip_diagnostic_call       = 1
INTEGER, PARAMETER :: ip_increment_call        = 2
INTEGER, PARAMETER :: ip_radiance_call         = 3

! Timestep number of the first step with a radiation call; assumed to be 1
! in the full UM, but may be set otherwise in the SCM:
INTEGER :: it_rad1 = 1

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RAD_INPUT_MOD'

CONTAINS

SUBROUTINE check_run_radiation()

! Description:
!   Subroutine to apply logic controls and set control variables based on the
!   options selected in the run_radiation namelist.

USE ereport_mod,  ONLY: ereport
USE fsd_parameters_mod, ONLY: ip_fsd_regime, ip_fsd_regime_no_sh,      &
                        ip_fsd_regime_smooth,                          &
                        ip_fsd_regime_smooth_no_sh
USE max_calls, ONLY: npd_swcall, npd_lwcall

IMPLICIT NONE

INTEGER                       :: icode         ! used for ereport
CHARACTER (LEN=errormessagelength)   :: cmessage      ! used for ereport
CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'CHECK_RUN_RADIATION'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set logicals for different radiation packages based on a control
! integer choice in namelist
SELECT CASE (i_rad_extra_call)
  CASE (ip_single_call)
    l_forcing     = .FALSE.
    l_timestep    = .FALSE.
    l_rad_perturb = .FALSE.
    l_radiance    = .FALSE.
    n_lwcall      = 1
    n_swcall      = 1
  CASE (ip_diagnostic_call)
    l_forcing     = .TRUE.
    l_timestep    = .FALSE.
    l_rad_perturb = .FALSE.
    l_radiance    = .FALSE.
    n_lwcall      = npd_lwcall
    n_swcall      = npd_swcall
  CASE (ip_increment_call)
    l_forcing     = .FALSE.
    l_timestep    = .TRUE.
    l_rad_perturb = .TRUE.
    l_radiance    = .FALSE.
    n_lwcall      = 2
    n_swcall      = 2
  CASE (ip_radiance_call)
    l_forcing     = .FALSE.
    l_timestep    = .FALSE.
    l_rad_perturb = .FALSE.
    l_radiance    = .TRUE.
    n_lwcall      = npd_lwcall
    n_swcall      = npd_swcall
  CASE DEFAULT
    WRITE (cmessage,'(A,I1,A)')                                         &
            'i_rad_extra_call value invalid, default to ',              &
            ip_single_call, ': single call to radiation'
    icode = -100
    CALL ereport(RoutineName, icode, cmessage)
END SELECT

! Set logicals for different radiation topography based on a control
! integer choice in namelist
IF (i_rad_topography == 0) THEN
  l_use_orog_corr   = .FALSE.
  l_use_grad_corr   = .FALSE.
  l_use_skyview     = .FALSE.
  l_orog_unfilt     = .FALSE.
ELSE IF (i_rad_topography == 1) THEN
  l_use_orog_corr   = .TRUE.
  l_use_grad_corr   = .FALSE.
  l_use_skyview     = .FALSE.
  l_orog_unfilt     = .FALSE.
ELSE IF (i_rad_topography == 2) THEN
  l_use_orog_corr   = .FALSE.
  l_use_grad_corr   = .TRUE.
  l_use_skyview     = .FALSE.
  l_orog_unfilt     = .FALSE.
ELSE IF (i_rad_topography == 3) THEN
  l_use_orog_corr   = .TRUE.
  l_use_grad_corr   = .FALSE.
  l_use_skyview     = .TRUE.
  l_orog_unfilt     = .FALSE.
ELSE IF (i_rad_topography == 4) THEN
  l_use_orog_corr   = .FALSE.
  l_use_grad_corr   = .TRUE.
  l_use_skyview     = .TRUE.
  l_orog_unfilt     = .TRUE.
ELSE
  WRITE (cmessage,'(A58)') 'i_rad_topography value invalid, default to ' &
                            // '0: flat surface'
  icode = -100
  CALL ereport(RoutineName, icode, cmessage)
END IF

! Warn if duplicate effects selected

IF (l_use_arclbiom .AND. l_use_bmass_direct) THEN
  WRITE (cmessage,'(A)') 'arclbiom and bmass_direct should not both be true'
  icode = -110
  CALL ereport(RoutineName, icode, cmessage)
END IF

IF (l_use_arclblck .AND. l_use_soot_direct) THEN
  WRITE (cmessage,'(A)') 'arclblck and soot_direct should not both be true'
  icode = -110
  CALL ereport(RoutineName, icode, cmessage)
END IF

IF (l_use_arclocff .AND. l_use_ocff_direct) THEN
  WRITE (cmessage,'(A)') 'arclocff and ocff_direct should not both be true'
  icode = -110
  CALL ereport(RoutineName, icode, cmessage)
END IF

IF (l_use_arclsslt .AND. l_use_seasalt_direct) THEN
  WRITE (cmessage,'(A)') 'arclsslt and seasalt_direct should not both be true'
  icode = -110
  CALL ereport(RoutineName, icode, cmessage)
END IF

IF (l_use_arclsulp .AND. l_use_sulpc_direct) THEN
  WRITE (cmessage,'(A)') 'arclsulp and sulpc_direct should not both be true'
  icode = -110
  CALL ereport(RoutineName, icode, cmessage)
END IF


! Set logicals for whether cca from different convection types are
! required by the cloud generator
IF ((i_fsd == ip_fsd_regime) .OR. (i_fsd == ip_fsd_regime_smooth)) THEN
  l_cca_sh_prog = .TRUE.
  l_cca_dp_prog = .TRUE.
  l_cca_md_prog = .TRUE.
ELSE IF ((i_fsd == ip_fsd_regime_no_sh)  .OR.                            &
         (i_fsd == ip_fsd_regime_smooth_no_sh)) THEN
  l_cca_dp_prog = .TRUE.
  l_cca_md_prog = .TRUE.
END IF

! Should also check validity of aerosol radiative effect choices,
! but this needs to be done after run_aerosol has been read.

! Set radiation aerosol switches
IF (cusack_aero==2 .OR.  cusack_aero==3    ) l_climat_aerosol    = .TRUE.
IF (cusack_aero==3 .AND. cusack_aero_hgt==1) l_clim_aero_hgt     = .TRUE.
IF (cusack_aero==2)                          l_HadGEM1_clim_aero = .TRUE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_run_radiation


SUBROUTINE print_nlist_run_radiation()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_RUN_RADIATION'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist run_radiation', &
    src='rad_input_mod')

WRITE(lineBuffer,*)' cusack_aero = ',cusack_aero
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' cusack_aero_hgt = ',cusack_aero_hgt
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' co2_mmr = ',co2_mmr
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' l_sec_var = ',l_sec_var
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' inhom_cloud_sw = ',inhom_cloud_sw
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' inhom_cloud_lw = ',inhom_cloud_lw
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' dp_corr_strat = ',dp_corr_strat
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' rad_mcica_sampling = ',rad_mcica_sampling
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' rad_mcica_sigma = ',rad_mcica_sigma
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' two_d_fsd_factor = ',two_d_fsd_factor
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' l_fsd_eff_res = ',l_fsd_eff_res
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' dp_corr_conv = ',dp_corr_conv
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' aero_bl_levels = ',aero_bl_levels
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' l_rad_use_clim_volc = ',l_rad_use_clim_volc
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' clim_rad_volc_eruption_year = ', &
    clim_rad_volc_eruption_year
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' clim_rad_volc_eruption_month = ', &
    clim_rad_volc_eruption_month
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' clim_rad_volc_eruption_weight = ', &
    clim_rad_volc_eruption_weight
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' n2ommr = ',n2ommr
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' ch4mmr = ',ch4mmr
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' c11mmr = ',c11mmr
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' c12mmr = ',c12mmr
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' o2mmr = ',o2mmr
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' so2mmr = ',so2mmr
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' c113mmr = ',c113mmr
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' c114mmr = ',c114mmr
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' hcfc22mmr = ',hcfc22mmr
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' hfc125mmr = ',hfc125mmr
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' hfc134ammr = ',hfc134ammr
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' l_use_sulpc_direct = ',l_use_sulpc_direct
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' l_use_sulpc_indirect_sw = ',l_use_sulpc_indirect_sw
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' l_use_sulpc_indirect_lw = ',l_use_sulpc_indirect_lw
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' l_use_seasalt_direct = ',l_use_seasalt_direct
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' l_use_seasalt_indirect = ',l_use_seasalt_indirect
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' l_use_biogenic = ',l_use_biogenic
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*) ' l_use_ocff_direct = ',l_use_ocff_direct
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*) ' l_use_ocff_indirect = ',l_use_ocff_indirect
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*) ' l_use_nitrate_direct = ',l_use_nitrate_direct
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*) ' l_use_nitrate_indirect = ',l_use_nitrate_indirect
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*) ' l_consistent_cdnc = ',l_consistent_cdnc
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*) ' l_use_soot_direct = ',l_use_soot_direct
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*) ' l_use_soot_indirect = ',l_use_soot_indirect
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*) ' l_use_bmass_direct = ',l_use_bmass_direct
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*) ' l_use_bmass_indirect = ',l_use_bmass_indirect
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*) ' l_use_dust = ',l_use_dust
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*) ' l_use_arclbiom = ',l_use_arclbiom
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*) ' l_use_arclblck = ',l_use_arclblck
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*) ' l_use_arcldlta = ',l_use_arcldlta
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*) ' l_use_arcldust = ',l_use_arcldust
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*) ' l_use_arclocff = ',l_use_arclocff
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*) ' l_use_arclsslt = ',l_use_arclsslt
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*) ' l_use_arclsulp = ',l_use_arclsulp
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*) ' l_use_cariolle = ',l_use_cariolle
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*) ' l_bs1999_abundances = ',l_bs1999_abundances
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,'(A,L7)') ' l_hydrostatic_mass = ',l_hydrostatic_mass
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,'(A,L7)') ' l_moist_heat_capacity = ',l_moist_heat_capacity
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,'(A,L7)') ' l_extra_top = ',l_extra_top
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,'(A,L7)') ' l_use_liu_spec = ',l_use_liu_spec
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' aparam = ',aparam
CALL umPrint(lineBuffer,src='rad_input_mod')
WRITE(lineBuffer,*)' bparam = ',bparam
CALL umPrint(lineBuffer,src='rad_input_mod')

CALL umPrint('- - - - - - end of namelist - - - - - -', &
    src='rad_input_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_run_radiation

#if !defined(LFRIC)
SUBROUTINE read_nml_run_radiation(unit_in)

USE um_parcore, ONLY: mype

USE check_iostat_mod, ONLY: check_iostat

USE setup_namelist, ONLY: setup_nml_type

IMPLICIT NONE

INTEGER,INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_RUN_RADIATION'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int = 23
INTEGER, PARAMETER :: n_real = 19 + 5 + 2 * npd_cloud_component
INTEGER, PARAMETER :: n_log = 41

TYPE my_namelist
  SEQUENCE
  INTEGER :: cusack_aero
  INTEGER :: cusack_aero_hgt
  INTEGER :: rad_mcica_sampling
  INTEGER :: aero_bl_levels
  INTEGER :: clim_rad_volc_eruption_year
  INTEGER :: clim_rad_volc_eruption_month
  INTEGER :: i_cloud_representation
  INTEGER :: i_cloud_representation_2
  INTEGER :: i_inhom
  INTEGER :: i_inhom_2
  INTEGER :: i_overlap
  INTEGER :: i_overlap_2
  INTEGER :: i_fsd
  INTEGER :: i_fsd_2
  INTEGER :: h_swbands
  INTEGER :: h_lwbands
  INTEGER :: i_rad_extra_call
  INTEGER :: i_rad_topography
  INTEGER :: i_sw_radstep_perday_prog
  INTEGER :: i_lw_radstep_perday_prog
  INTEGER :: i_sw_radstep_perday_diag
  INTEGER :: i_lw_radstep_perday_diag
  INTEGER :: i_ozone_int
  REAL :: aeroscl_csk_clim(5)
  REAL :: co2_mmr
  REAL :: inhom_cloud_sw(npd_cloud_component)
  REAL :: inhom_cloud_lw(npd_cloud_component)
  REAL :: dp_corr_strat
  REAL :: rad_mcica_sigma
  REAL :: two_d_fsd_factor
  REAL :: dp_corr_conv
  REAL :: clim_rad_volc_eruption_weight
  REAL :: n2ommr
  REAL :: ch4mmr
  REAL :: c11mmr
  REAL :: c12mmr
  REAL :: o2mmr
  REAL :: so2mmr
  REAL :: c113mmr
  REAL :: c114mmr
  REAL :: hcfc22mmr
  REAL :: hfc125mmr
  REAL :: hfc134ammr
  REAL :: aparam
  REAL :: bparam
  LOGICAL :: l_sec_var
  LOGICAL :: l_rad_use_clim_volc
  LOGICAL :: l_rad_snow_emis
  LOGICAL :: l_t_land_nosnow
  LOGICAL :: l_quad_t_coast
  LOGICAL :: l_t_rad_solid
  LOGICAL :: l_t_bdy_surf
  LOGICAL :: l_radiation
  LOGICAL :: l_rad_deg
  LOGICAL :: l_rad_szacor
  LOGICAL :: l_use_cariolle
  LOGICAL :: l_use_ozoneinrad
  LOGICAL :: l_consistent_cdnc
  LOGICAL :: l_use_sulpc_direct
  LOGICAL :: l_use_sulpc_indirect_sw
  LOGICAL :: l_use_sulpc_indirect_lw
  LOGICAL :: l_use_seasalt_direct
  LOGICAL :: l_use_seasalt_indirect
  LOGICAL :: l_use_biogenic
  LOGICAL :: l_use_nitrate_direct
  LOGICAL :: l_use_nitrate_indirect
  LOGICAL :: l_use_soot_direct
  LOGICAL :: l_use_soot_indirect
  LOGICAL :: l_use_bmass_direct
  LOGICAL :: l_use_bmass_indirect
  LOGICAL :: l_use_ocff_direct
  LOGICAL :: l_use_ocff_indirect
  LOGICAL :: l_use_dust
  LOGICAL :: l_use_arclbiom
  LOGICAL :: l_use_arclblck
  LOGICAL :: l_use_arcldlta
  LOGICAL :: l_use_arcldust
  LOGICAL :: l_use_arclocff
  LOGICAL :: l_use_arclsslt
  LOGICAL :: l_use_arclsulp
  LOGICAL :: l_bs1999_abundances
  LOGICAL :: l_fsd_eff_res
  LOGICAL :: l_hydrostatic_mass
  LOGICAL :: l_moist_heat_capacity
  LOGICAL :: l_extra_top
  LOGICAL :: l_use_liu_spec
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,     &
                    n_real_in=n_real, n_log_in=n_log)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=RUN_Radiation, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist RUN_Radiation", iomessage)

  my_nml % cusack_aero        = cusack_aero
  my_nml % cusack_aero_hgt    = cusack_aero_hgt
  my_nml % rad_mcica_sampling = rad_mcica_sampling
  my_nml % aero_bl_levels     = aero_bl_levels
  my_nml % clim_rad_volc_eruption_year  = clim_rad_volc_eruption_year
  my_nml % clim_rad_volc_eruption_month = clim_rad_volc_eruption_month
  my_nml % i_cloud_representation   = i_cloud_representation
  my_nml % i_cloud_representation_2 = i_cloud_representation_2
  my_nml % i_inhom            = i_inhom
  my_nml % i_inhom_2          = i_inhom_2
  my_nml % i_overlap          = i_overlap
  my_nml % i_overlap_2        = i_overlap_2
  my_nml % i_fsd              = i_fsd
  my_nml % i_fsd_2            = i_fsd_2
  my_nml % h_swbands          = h_swbands
  my_nml % h_lwbands          = h_lwbands
  my_nml % i_rad_extra_call   = i_rad_extra_call
  my_nml % i_rad_topography   = i_rad_topography
  my_nml % i_sw_radstep_perday_prog = i_sw_radstep_perday_prog
  my_nml % i_lw_radstep_perday_prog = i_lw_radstep_perday_prog
  my_nml % i_sw_radstep_perday_diag = i_sw_radstep_perday_diag
  my_nml % i_lw_radstep_perday_diag = i_lw_radstep_perday_diag
  my_nml % i_ozone_int        = i_ozone_int
  ! end of integers
  my_nml % aeroscl_csk_clim = aeroscl_csk_clim
  my_nml % co2_mmr          = co2_mmr
  my_nml % inhom_cloud_sw   = inhom_cloud_sw
  my_nml % inhom_cloud_lw   = inhom_cloud_lw
  my_nml % dp_corr_strat    = dp_corr_strat
  my_nml % rad_mcica_sigma  = rad_mcica_sigma
  my_nml % two_d_fsd_factor = two_d_fsd_factor
  my_nml % dp_corr_conv     = dp_corr_conv
  my_nml % clim_rad_volc_eruption_weight = clim_rad_volc_eruption_weight
  my_nml % n2ommr           = n2ommr
  my_nml % ch4mmr           = ch4mmr
  my_nml % c11mmr           = c11mmr
  my_nml % c12mmr           = c12mmr
  my_nml % o2mmr            = o2mmr
  my_nml % so2mmr           = so2mmr
  my_nml % c113mmr          = c113mmr
  my_nml % c114mmr          = c114mmr
  my_nml % hcfc22mmr        = hcfc22mmr
  my_nml % hfc125mmr        = hfc125mmr
  my_nml % hfc134ammr       = hfc134ammr
  my_nml % aparam           = aparam
  my_nml % bparam           = bparam
  ! end of reals
  my_nml % l_sec_var           = l_sec_var
  my_nml % l_rad_use_clim_volc = l_rad_use_clim_volc
  my_nml % l_rad_snow_emis     = l_rad_snow_emis
  my_nml % l_t_land_nosnow     = l_t_land_nosnow
  my_nml % l_quad_t_coast      = l_quad_t_coast
  my_nml % l_t_rad_solid       = l_t_rad_solid
  my_nml % l_t_bdy_surf        = l_t_bdy_surf
  my_nml % l_radiation         = l_radiation
  my_nml % l_rad_deg           = l_rad_deg
  my_nml % l_rad_szacor        = l_rad_szacor
  my_nml % l_use_cariolle      = l_use_cariolle
  my_nml % l_use_ozoneinrad    = l_use_ozoneinrad
  my_nml % l_consistent_cdnc   = l_consistent_cdnc
  my_nml % l_use_sulpc_direct  = l_use_sulpc_direct
  my_nml % l_use_sulpc_indirect_sw = l_use_sulpc_indirect_sw
  my_nml % l_use_sulpc_indirect_lw = l_use_sulpc_indirect_lw
  my_nml % l_use_seasalt_direct    = l_use_seasalt_direct
  my_nml % l_use_seasalt_indirect  = l_use_seasalt_indirect
  my_nml % l_use_biogenic          = l_use_biogenic
  my_nml % l_use_nitrate_direct    = l_use_nitrate_direct
  my_nml % l_use_nitrate_indirect  = l_use_nitrate_indirect
  my_nml % l_use_soot_direct       = l_use_soot_direct
  my_nml % l_use_soot_indirect     = l_use_soot_indirect
  my_nml % l_use_bmass_direct      = l_use_bmass_direct
  my_nml % l_use_bmass_indirect    = l_use_bmass_indirect
  my_nml % l_use_ocff_direct       = l_use_ocff_direct
  my_nml % l_use_ocff_indirect     = l_use_ocff_indirect
  my_nml % l_use_dust              = l_use_dust
  my_nml % l_use_arclbiom          = l_use_arclbiom
  my_nml % l_use_arclblck          = l_use_arclblck
  my_nml % l_use_arcldlta          = l_use_arcldlta
  my_nml % l_use_arcldust          = l_use_arcldust
  my_nml % l_use_arclocff          = l_use_arclocff
  my_nml % l_use_arclsslt          = l_use_arclsslt
  my_nml % l_use_arclsulp          = l_use_arclsulp
  my_nml % l_bs1999_abundances     = l_bs1999_abundances
  my_nml % l_fsd_eff_res           = l_fsd_eff_res
  my_nml % l_hydrostatic_mass      = l_hydrostatic_mass
  my_nml % l_moist_heat_capacity   = l_moist_heat_capacity
  my_nml % l_extra_top             = l_extra_top
  my_nml % l_use_liu_spec          = l_use_liu_spec
END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  cusack_aero        = my_nml % cusack_aero
  cusack_aero_hgt    = my_nml % cusack_aero_hgt
  rad_mcica_sampling = my_nml % rad_mcica_sampling
  aero_bl_levels     = my_nml % aero_bl_levels
  clim_rad_volc_eruption_year  = my_nml % clim_rad_volc_eruption_year
  clim_rad_volc_eruption_month = my_nml % clim_rad_volc_eruption_month
  i_cloud_representation   = my_nml % i_cloud_representation
  i_cloud_representation_2 = my_nml % i_cloud_representation_2
  i_inhom            = my_nml % i_inhom
  i_inhom_2          = my_nml % i_inhom_2
  i_overlap          = my_nml % i_overlap
  i_overlap_2        = my_nml % i_overlap_2
  i_fsd              = my_nml % i_fsd
  i_fsd_2            = my_nml % i_fsd_2
  h_swbands          = my_nml % h_swbands
  h_lwbands          = my_nml % h_lwbands
  i_rad_extra_call   = my_nml % i_rad_extra_call
  i_rad_topography   = my_nml % i_rad_topography
  i_sw_radstep_perday_prog = my_nml % i_sw_radstep_perday_prog
  i_lw_radstep_perday_prog = my_nml % i_lw_radstep_perday_prog
  i_sw_radstep_perday_diag = my_nml % i_sw_radstep_perday_diag
  i_lw_radstep_perday_diag = my_nml % i_lw_radstep_perday_diag
  i_ozone_int        = my_nml % i_ozone_int
  ! end of integers
  aeroscl_csk_clim = my_nml % aeroscl_csk_clim
  co2_mmr          = my_nml % co2_mmr
  inhom_cloud_sw   = my_nml % inhom_cloud_sw
  inhom_cloud_lw   = my_nml % inhom_cloud_lw
  dp_corr_strat    = my_nml % dp_corr_strat
  rad_mcica_sigma  = my_nml % rad_mcica_sigma
  two_d_fsd_factor = my_nml % two_d_fsd_factor
  dp_corr_conv     = my_nml % dp_corr_conv
  clim_rad_volc_eruption_weight = my_nml % clim_rad_volc_eruption_weight
  n2ommr           = my_nml % n2ommr
  ch4mmr           = my_nml % ch4mmr
  c11mmr           = my_nml % c11mmr
  c12mmr           = my_nml % c12mmr
  o2mmr            = my_nml % o2mmr
  so2mmr           = my_nml % so2mmr
  c113mmr          = my_nml % c113mmr
  c114mmr          = my_nml % c114mmr
  hcfc22mmr        = my_nml % hcfc22mmr
  hfc125mmr        = my_nml % hfc125mmr
  hfc134ammr       = my_nml % hfc134ammr
  aparam           = my_nml % aparam
  bparam           = my_nml % bparam
  ! end of reals
  l_sec_var           = my_nml % l_sec_var
  l_rad_use_clim_volc = my_nml % l_rad_use_clim_volc
  l_rad_snow_emis     = my_nml % l_rad_snow_emis
  l_t_land_nosnow     = my_nml % l_t_land_nosnow
  l_quad_t_coast      = my_nml % l_quad_t_coast
  l_t_rad_solid       = my_nml % l_t_rad_solid
  l_t_bdy_surf        = my_nml % l_t_bdy_surf
  l_radiation         = my_nml % l_radiation
  l_rad_deg           = my_nml % l_rad_deg
  l_rad_szacor        = my_nml % l_rad_szacor
  l_use_cariolle      = my_nml % l_use_cariolle
  l_use_ozoneinrad    = my_nml % l_use_ozoneinrad
  l_consistent_cdnc   = my_nml % l_consistent_cdnc
  l_use_sulpc_direct  = my_nml % l_use_sulpc_direct
  l_use_sulpc_indirect_sw = my_nml % l_use_sulpc_indirect_sw
  l_use_sulpc_indirect_lw = my_nml % l_use_sulpc_indirect_lw
  l_use_seasalt_direct    = my_nml % l_use_seasalt_direct
  l_use_seasalt_indirect  = my_nml % l_use_seasalt_indirect
  l_use_biogenic          = my_nml % l_use_biogenic
  l_use_nitrate_direct    = my_nml % l_use_nitrate_direct
  l_use_nitrate_indirect  = my_nml % l_use_nitrate_indirect
  l_use_soot_direct       = my_nml % l_use_soot_direct
  l_use_soot_indirect     = my_nml % l_use_soot_indirect
  l_use_bmass_direct      = my_nml % l_use_bmass_direct
  l_use_bmass_indirect    = my_nml % l_use_bmass_indirect
  l_use_ocff_direct       = my_nml % l_use_ocff_direct
  l_use_ocff_indirect     = my_nml % l_use_ocff_indirect
  l_use_dust              = my_nml % l_use_dust
  l_use_arclbiom          = my_nml % l_use_arclbiom
  l_use_arclblck          = my_nml % l_use_arclblck
  l_use_arcldlta          = my_nml % l_use_arcldlta
  l_use_arcldust          = my_nml % l_use_arcldust
  l_use_arclocff          = my_nml % l_use_arclocff
  l_use_arclsslt          = my_nml % l_use_arclsslt
  l_use_arclsulp          = my_nml % l_use_arclsulp
  l_bs1999_abundances     = my_nml % l_bs1999_abundances
  l_fsd_eff_res           = my_nml % l_fsd_eff_res
  l_hydrostatic_mass      = my_nml % l_hydrostatic_mass
  l_moist_heat_capacity   = my_nml % l_moist_heat_capacity
  l_extra_top             = my_nml % l_extra_top
  l_use_liu_spec          = my_nml % l_use_liu_spec
END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_run_radiation
#endif

END MODULE rad_input_mod

