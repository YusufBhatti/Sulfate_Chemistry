! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine microphys_ctl: CASIM VERSION
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Precipitation (CASIM)

MODULE casim_ctl_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CASIM_CTL_MOD'

CONTAINS

SUBROUTINE casim_ctl(                                                   &
! Primary fields passed in
     t_n, q_n, qcl_n, qcf_n, qcf2_n, qrain_n, qgraup_n, w,              &
     p_layer_centres, rho_r2, exner_theta_levels,                       &
     flash_pot, tracers,                                                &
! diagnostic info
     stashwork4, stashwork21,                                           &
! Increment fields passed in/out
     t_inc, q_inc, qcl_inc, qcf_inc, qcf2_inc, qrain_inc, qgraup_inc,   &
! Fields required elsewhere
     ls_rain, ls_snow, ls_graup, micro_tends, n_drop_pot )


USE nlsizes_namelist_mod,  ONLY: bl_levels, tr_vars
USE um_parvars,            ONLY: at_extremity
USE stash_array_mod,       ONLY: sf
USE model_domain_mod,      ONLY: model_type, mt_single_column

! General atmosphere modules
USE planet_constants_mod,  ONLY: lcrcp, lsrcp
USE mphys_inputs_mod,      ONLY: casim_aerosol_process_level,           &
                                 casim_aerosol_option, l_mcr_qgraup
USE mphys_air_density_mod, ONLY: mphys_air_density
USE electric_inputs_mod,   ONLY: l_use_electric
USE electric_main_mod,     ONLY: electric_main
USE gen_phys_inputs_mod,   ONLY: l_mr_physics
USE timestep_mod,          ONLY: timestep, recip_timestep
USE rad_input_mod,         ONLY: l_consistent_cdnc

! Grid bounds module
USE atm_fields_bounds_mod, ONLY: tdims, tdims_s, wdims_s, pdims_s
USE micro_main,            ONLY: shipway_microphysics

USE casim_switches,        ONLY:  l_mp_CloudNumber, l_mp_RainNumber,        &
                                  l_mp_Rain3mom,    l_mp_IceNumber,         &
                                  l_mp_SnowNumber,  l_mp_Snow3mom,          &
                                  l_mp_GraupNumber, l_mp_Graup3mom,         &
                                  l_mp_ActiveSolLiquid, l_mp_ActiveSolRain, &
                                  l_mp_ActiveInsolIce, l_mp_ActiveSolIce,   &
                                  l_mp_ActiveInsolLiquid,                   &
                                  l_mp_ActiveSolNumber,                     &
                                  l_mp_ActiveInSolNumber,                   &
                                  its, ite, jts, jte, kts, kte,             &
                                  ils, ile, jls, jle, kls, kle,             &
                                  irs, ire, jrs, jre, krs, kre,             &
                                  l_fix_aerosol, l_casim_warm_only,         &
                                  rain_mom, snow_mom, graup_mom,            &
                                  single_moment, double_moment,             &
                                  triple_moment, no_processing

USE casim_prognostics,     ONLY:  CloudNumber, RainNumber, Rain3mom,        &
                                  IceNumber, SnowNumber,                    &
                                  Snow3mom, GraupNumber, Graup3mom,         &
                                  CloudNumber_inc, RainNumber_inc,          &
                                  Rain3mom_inc, IceNumber_inc,              &
                                  SnowNumber_inc, Snow3mom_inc,             &
                                  GraupNumber_inc, Graup3mom_inc,           &
                                  ActiveSolLiquid, ActiveSolRain,           &
                                  ActiveInsolIce, ActiveSolIce,             &
                                  ActiveInsolLiquid, ActiveSolNumber,       &
                                  ActiveInsolNumber,                        &
                                  dActiveSolLiquid, dActiveSolRain,         &
                                  dActiveInsolIce, dActiveSolIce,           &
                                  dActiveInsolLiquid, dActiveSolNumber,     &
                                  dActiveInsolNumber

USE diagnostics_casim_mod,  ONLY: diagnostics_casim

USE generic_diagnostic_variables, ONLY: allocate_diagnostic_space,          &
                                        deallocate_diagnostic_space,        &
                                        casdiags

! Dr Hook modules
USE yomhook,               ONLY: lhook, dr_hook
USE parkind1,              ONLY: jprb, jpim

IMPLICIT NONE

!------------------------------------------------------------------------------
! Subroutine arguments
!------------------------------------------------------------------------------
! Input of temperature [K]
REAL, INTENT(IN) ::  t_n( tdims%i_start : tdims%i_end,                         &
                          tdims%j_start : tdims%j_end,                         &
                                      1 : tdims%k_end )
! Input of vapour mixing ratio [kg kg-1]
REAL, INTENT(IN) ::  q_n( tdims%i_start : tdims%i_end,                         &
                          tdims%j_start : tdims%j_end,                         &
                                      1 : tdims%k_end )

! Input of cloud liquid mixing ratio [kg kg-1]
REAL, INTENT(IN) :: qcl_n( tdims%i_start : tdims%i_end,                        &
                           tdims%j_start : tdims%j_end,                        &
                                       1 : tdims%k_end )

! Input of cloud ice aggregate mixing ratio [kg kg-1]
REAL, INTENT(IN) :: qcf_n( tdims%i_start : tdims%i_end,                        &
                           tdims%j_start : tdims%j_end,                        &
                                       1 : tdims%k_end )

! Input of cloud ice crystal mixing ratio [kg kg-1]
REAL, INTENT(IN) :: qcf2_n( tdims%i_start : tdims%i_end,                       &
                            tdims%j_start : tdims%j_end,                       &
                                        1 : tdims%k_end )

! Input of rain mixing ratio [kg kg-1]
REAL, INTENT(IN) :: qrain_n( tdims%i_start : tdims%i_end,                      &
                             tdims%j_start : tdims%j_end,                      &
                                         1 : tdims%k_end )

! Input of graupel mixing ratio [kg kg-1]
REAL, INTENT(IN) :: qgraup_n( tdims%i_start : tdims%i_end,                     &
                              tdims%j_start : tdims%j_end,                     &
                                          1 : tdims%k_end )

! Vertical component of velocity [m s-1]
REAL, INTENT(IN) :: w( wdims_s%i_start : wdims_s%i_end,                        &
                       wdims_s%j_start : wdims_s%j_end,                        &
                       wdims_s%k_start : wdims_s%k_end ) 

! Pressure at layer centres [Pa]
REAL, INTENT(IN) ::  p_layer_centres( tdims%i_start : tdims%i_end,             &
                                      tdims%j_start : tdims%j_end,             &
                                      tdims%k_start : tdims%k_end )

! Exner function on theta levels
REAL, INTENT(IN) :: exner_theta_levels( tdims_s%i_start : tdims_s%i_end,       &
                                        tdims_s%j_start : tdims_s%j_end,       &
                                                      1 : tdims_s%k_end )
! Air density * (Earth radius ^ 2) [kg m]                               
REAL, INTENT(IN) :: rho_r2(  pdims_s%i_start : pdims_s%i_end,                  &
                             pdims_s%j_start : pdims_s%j_end,                  &
                             pdims_s%k_start : pdims_s%k_end )

! Lightning Flash Potential [] - for use with lightning scheme
REAL, INTENT(INOUT) :: flash_pot(  tdims%i_start : tdims%i_end,                &
                                   tdims%j_start : tdims%j_end,                &
                                   tdims%k_start : tdims%k_end )

! Free tracers []
REAL, INTENT(INOUT) ::  tracers( tdims_s%i_start : tdims_s%i_end,              &
                                 tdims_s%j_start : tdims_s%j_end,              &
                                 tdims_s%k_start : tdims_s%k_end, tr_vars )

! Diagnostics STASH sections
REAL, INTENT(INOUT) ::  stashwork4(*)     ! STASH workspace (Microphysics)
REAL, INTENT(INOUT) ::  stashwork21(*)    ! STASH workspace (Electrification)

! Temperature increment [K]
REAL, INTENT(INOUT) ::  t_inc( tdims%i_start : tdims%i_end,                    &
                               tdims%j_start : tdims%j_end,                    &
                                           1 : tdims%k_end )

! Vapour mixing ratio increment [kg kg-1]
REAL, INTENT(INOUT) ::  q_inc( tdims%i_start : tdims%i_end,                    &
                               tdims%j_start : tdims%j_end,                    &
                                           1 : tdims%k_end )

! Cloud Liquid water mixing ratio increment [kg kg-1]
REAL, INTENT(INOUT) ::  qcl_inc( tdims%i_start : tdims%i_end,                  &
                                 tdims%j_start : tdims%j_end,                  &
                                             1 : tdims%k_end )

! Cloud Ice aggregate mixing ratio increment [kg kg-1]
REAL, INTENT(INOUT) ::  qcf_inc( tdims%i_start : tdims%i_end,                  &
                                 tdims%j_start : tdims%j_end,                  &
                                             1 : tdims%k_end )

! Cloud Ice crystal mixing ratio increment [kg kg-1]
REAL, INTENT(INOUT) ::  qcf2_inc( tdims%i_start : tdims%i_end,                 &
                                  tdims%j_start : tdims%j_end,                 &
                                              1 : tdims%k_end )

! Rain mixing ratio increment [kg kg-1]
REAL, INTENT(INOUT) ::  qrain_inc( tdims%i_start : tdims%i_end,                &
                                   tdims%j_start : tdims%j_end,                &
                                               1 : tdims%k_end )
! Graupel mixing ratio increment [kg kg-1]
REAL, INTENT(INOUT) ::  qgraup_inc( tdims%i_start : tdims%i_end,               &
                                    tdims%j_start : tdims%j_end,               &
                                                1 : tdims%k_end )

! Rainfall rate to pass out to other sections (e.g. Land surface) [kg m-2 s-1]
REAL, INTENT(OUT) :: ls_rain( tdims%i_start : tdims%i_end,                     &
                              tdims%j_start : tdims%j_end )

! Snowfall rate to pass out to other sections (e.g. Land surface) [kg m-2 s-1]
REAL, INTENT(OUT) :: ls_snow( tdims%i_start : tdims%i_end,                     &
                              tdims%j_start : tdims%j_end )


! This should be made intent out on upgrade to vn10.7


REAL, INTENT(OUT) :: ls_graup( tdims%i_start : tdims%i_end,                    &
                               tdims%j_start : tdims%j_end )

! Microphysical tendencies for boundary layer within BL levels
! [TL, K/s; QW, kg/kg/s]
! The number 2 refers to the two categories
! 1 is TL and 2 is QW
REAL, INTENT(OUT) ::  micro_tends( tdims%i_start : tdims%i_end,                &
                                   tdims%j_start : tdims%j_end,                &
                                   2, bl_levels )

! Potential cloud droplet number concentration, for use in Radiation
! [m-3]               
REAL, INTENT(OUT) ::  n_drop_pot( tdims%i_start : tdims%i_end,                 &
                                  tdims%j_start : tdims%j_end,                 &
                                              1 : tdims%k_end )

!------------------------------------------------------------------------------
! Local Variables
!------------------------------------------------------------------------------

INTEGER, PARAMETER :: stash_sec = 4 ! STASH section for microphysics
INTEGER, PARAMETER :: izero     = 0 ! Zero STASH flag item

! Fixed accumulation mode quantities
REAL, PARAMETER :: fixed_accum_sol_mass   = 1.5e-9
REAL, PARAMETER :: fixed_accum_sol_number = 100.0e6

! Dr Hook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER   :: RoutineName='CASIM_CTL'

! loop counters
INTEGER :: i, j, k

! Fixed value of TKE [m2 s-2]
REAL, PARAMETER  :: fixed_tke_value = 0.1

! Local working temperature [K]
REAL  ::  t_work( tdims%i_start : tdims%i_end,                                 &
                  tdims%j_start : tdims%j_end,                                 &
                              1 : tdims%k_end )

! Local working cloud ice aggregate content [kg kg-1]
REAL, ALLOCATABLE :: qcf_work(:,:,:)  

! Local working second ice content [kg kg-1]
REAL, ALLOCATABLE :: qcf2_work(:,:,:)

! Local working graupel content [kg kg-1]
REAL, ALLOCATABLE :: qgraup_work(:,:,:)

! Dry density [kg m-3]
REAL :: rhodz_dry( tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end )

! Moist density [kg m-3]
REAL :: rhodz_moist( tdims%i_start : tdims%i_end,                             &
                     tdims%j_start : tdims%j_end,                             &
                                 1 : tdims%k_end )

! Layer thickness [m]
REAL :: deltaz( tdims%i_start : tdims%i_end,                                  &
                tdims%j_start : tdims%j_end,                                  &
                            1 : tdims%k_end )

! -----------------------------------------------------------------------------
! Local CASIM microphysics variables
! These are all of the column-first form (z, x, y) hence
! are defined on these dimensions appropriately using tdims
! They are then passed in/out of the CASIM microphysics scheme
! and copied back to the (x, y, z) form that is native to the UM. 
!------------------------------------------------------------------------------
! Height on theta levels [m]
REAL :: height(               1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Height on rho levels [m]
REAL :: height_rho(           0 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Theta on CASIM dimensions [k]
REAL :: th_casim(             1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Pressure on CASIM dimensions [Pa]
REAL :: p_casim(              1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Layer thickness on CASIM dimensions [m]
REAL :: dz_casim(             1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Vapour mixing ratio on CASIM dimensions [kg kg-1]
REAL :: qv_casim(             1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Cloud liquid content mixing ratio on CASIM dimensions [kg kg-1]
REAL :: qc_casim(             1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Cloud number concentration on CASIM dimensions [m-3]
REAL :: nc_casim(             1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Rain mixing ratio on CASIM dimensions [kg kg-1]
REAL :: qr_casim(             1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Rain number on CASIM dimensions [m-3]
REAL :: nr_casim(             1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Rain third moment on CASIM dimensions 
REAL :: m3r_casim(            1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Air density on CASIM dimensions [kg m-3]
REAL :: rho_casim(            1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Exner function on CASIM dimensions
REAL :: pii_casim(            1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Vertical velocity component on CASIM dimensions [m s-1]
REAL :: w_casim(              1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Turbulent kinetic energy dissipation rate on CASIM dimensions [m2 s-3]
REAL :: tke_casim(            1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Ice mixing ratio on CASIM dimensions [kg kg-1]
REAL :: qi_casim(             1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Snow mixing ratio on CASIM dimensions [kg kg-1]
REAL :: qs_casim(             1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Graupel mixing ratio on CASIM dimensions [kg kg-1]
REAL :: qg_casim(             1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Ice number on CASIM dimensions [m-3]
REAL :: ni_casim(             1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Snow number on CASIM dimensions [m-3]
REAL :: ns_casim(             1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Graupel number on CASIM dimensions [m-3]
REAL :: ng_casim(             1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Snow third moment on CASIM dimensions
REAL :: m3s_casim(            1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Graupel third moment on CASIM dimensions
REAL :: m3g_casim(            1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Change in Theta on CASIM dimensions [k]
REAL :: dth_casim(            1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Change in Vapour mixing ratio on CASIM dimensions [kg kg-1]
REAL :: dqv_casim(            1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Change in Cloud liquid content mixing ratio on CASIM dimensions [kg kg-1]
REAL :: dqc_casim(            1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Change in Cloud number concentration on CASIM dimensions [m-3]
REAL :: dnc_casim(            1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Change in Rain mixing ratio on CASIM dimensions [kg kg-1]
REAL :: dqr_casim(            1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Change in Rain number on CASIM dimensions [m-3]
REAL :: dnr_casim(            1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Change in Rain third moment on CASIM dimensions 
REAL :: dm3r_casim(           1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Change in Ice mixing ratio on CASIM dimensions [kg kg-1]
REAL :: dqi_casim(            1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Change in Snow mixing ratio on CASIM dimensions [kg kg-1]
REAL :: dqs_casim(            1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Change in Graupel mixing ratio on CASIM dimensions [kg kg-1]
REAL :: dqg_casim(            1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Change in Ice number on CASIM dimensions [m-3]
REAL :: dni_casim(            1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Change in Snow number on CASIM dimensions [m-3]
REAL :: dns_casim(            1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Change in Graupel number on CASIM dimensions [m-3]
REAL :: dng_casim(            1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Change in Snow third moment on CASIM dimensions
REAL :: dm3s_casim(           1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Change in Graupel third moment on CASIM dimensions
REAL :: dm3g_casim(           1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Activated soluble liquid on CASIM dimensions
REAL :: ActSolLiq_casim(      1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Change in Activated soluble liquid on CASIM dimensions
REAL :: dActSolLiq_casim(     1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Activated soluble rain on CASIM dimensions
REAL :: ActSolRain_casim(     1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Change in Activated soluble rain on CASIM dimensions
REAL :: dActSolRain_casim(    1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Activated insoluble ice on CASIM dimensions
REAL :: ActInsolIce_casim(    1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Change in Activated insoluble ice on CASIM dimensions
REAL :: dActInsolIce_casim(   1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Activated soluble ice on CASIM dimensions
REAL :: ActSolIce_casim(      1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Change in Activated soluble ice on CASIM dimensions
REAL :: dActSolIce_casim(     1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Activated insoluble liquid on CASIM dimensions
REAL :: ActInsolLiq_casim(    1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Change in Activated insoluble liquid on CASIM dimensions
REAL :: dActInsolLiq_casim(   1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Activated soluble number on CASIM dimensions
REAL :: ActSolNumber_casim(   1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Activated insoluble number on CASIM dimensions
REAL :: ActInsolNumber_casim( 1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Change in Activated soluble number on CASIM dimensions
REAL :: dActSolNumber_casim(  1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Change in Activated insoluble number on CASIM dimensions
REAL :: dActInsolNumber_casim(1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Accumulation mode aerosol mass
REAL :: AccumSolMass         (1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Accumulation model aerosol number
REAL :: AccumSolNumber       (1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Aitken mode aerosol mass
REAL  :: AitkenSolMass       (1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Aitken model aerosol number
REAL :: AitkenSolNumber      (1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )


! Coarse mode aerosol mass
REAL  :: CoarseSolMass       (1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Coarse mode aerosol number
REAL :: CoarseSolNumber      (1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Coarse mode Dust mass
REAL :: CoarseDustMass       (1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Coarse mode Dust number
REAL :: CoarseDustNumber     (1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Accumulation mode dust mass
REAL  :: AccumDustMass       (1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Accumulation mode dust number
REAL :: AccumDustNumber      (1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Change in Accumulation mode aerosol mass
REAL :: dAccumSolMass        (1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Change in Accumulation model aerosol number
REAL :: dAccumSolNumber      (1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Change in Aitken mode aerosol mass
REAL  :: dAitkenSolMass      (1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Change in Aitken model aerosol number
REAL :: dAitkenSolNumber     (1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )


! Change in Coarse mode aerosol mass
REAL  :: dCoarseSolMass      (1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Change in Coarse mode aerosol number
REAL :: dCoarseSolNumber     (1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Change in Coarse mode Dust mass
REAL :: dCoarseDustMass      (1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Change in Coarse mode Dust number
REAL :: dCoarseDustNumber    (1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Change in Accumulation mode dust mass
REAL  :: dAccumDustMass      (1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

! Change in Accumulation mode dust number
REAL :: dAccumDustNumber     (1 : tdims%k_end,                                &
                  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )

!==============================================================================
! End of Declarations and start of the subroutine doing something
!==============================================================================

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Ensure all diagnostics from casdiags are set to their appropriate STASH
! flag before calling the routine to allocate the diagnostic space for
! non-SCM runs only where diagnostic output is required on this timestep.

IF ( model_type /= mt_single_column .AND. sf(izero, stash_sec) ) THEN

  ! Set switch for radar reflectivity diagnostics if any are in use.
  ! N.B. casdiags % l_radar will be set to .FALSE. in the CASIM routine
  ! deallocate_diagnostic_space, which is called at the end of casim_ctl.  
  IF (sf(110, stash_sec) .OR. sf(111, stash_sec) .OR. sf(112, stash_sec) .OR. &
      sf(113, stash_sec) .OR. sf(114, stash_sec) .OR. sf(115, stash_sec) .OR. &
      sf(116, stash_sec) .OR. sf(117, stash_sec) .OR. sf(118, stash_sec) ) THEN

    casdiags % l_radar = .TRUE.

  END IF ! radar reflectivity stash codes 110-118.

  ! Set logical variables for process rate diagnostics
  ! All process rates diagnostic flags are set to .FALSE. in the
  ! deallocate_diagnostic_space routine.
  ! All the following diagnostics have units of [kg/kg/s].
  casdiags % l_phomc = sf(240, stash_sec) ! Homogeneous nucleation rate
  casdiags % l_pinuc = sf(241, stash_sec) ! Heterogeneous nucleation rate
  casdiags % l_pidep = sf(243, stash_sec) ! Deposition rate ice crystals
  casdiags % l_psdep = sf(245, stash_sec) ! Deposition rate snow aggregates
  casdiags % l_piacw = sf(247, stash_sec) ! Riming rate ice crystals
  casdiags % l_psacw = sf(248, stash_sec) ! Riming rate snow aggregates
  casdiags % l_psacr = sf(250, stash_sec) ! Snow aggregate-rain capture rate
  casdiags % l_pisub = sf(251, stash_sec) ! Evaporation (sublimation) of 
                                          ! melting ice crystals
  casdiags % l_pssub = sf(252, stash_sec) ! Evaporation (sublimation) of 
                                          ! melting snow aggregates.
  casdiags % l_pimlt = sf(253, stash_sec) ! Melting rate of ice crystals
  casdiags % l_psmlt = sf(254, stash_sec) ! Melting rate of snow aggregates
  casdiags % l_psaut = sf(255, stash_sec) ! Snow autoconversion rate
  casdiags % l_psaci = sf(256, stash_sec) ! Snow capture of ice crystals
  casdiags % l_praut = sf(257, stash_sec) ! Rain autoconversion rate
  casdiags % l_pracw = sf(258, stash_sec) ! Rain accretion rate
  casdiags % l_prevp = sf(259, stash_sec) ! Rain evaporation rate
  casdiags % l_pgacw = sf(261, stash_sec) ! Graupel accretion rate
  casdiags % l_pgacs = sf(262, stash_sec) ! Graupel-snow capture rate
  casdiags % l_pgmlt = sf(263, stash_sec) ! Graupel melting rate
  casdiags % l_pgsub = sf(264, stash_sec) ! Graupel evap. (sublimation) rate
  casdiags % l_psedi = sf(265, stash_sec) ! Ice crystal sedimentation rate
  casdiags % l_pseds = sf(266, stash_sec) ! Snow aggregate sedimentation rate
  casdiags % l_psedr = sf(267, stash_sec) ! Rain sedimentation rate
  casdiags % l_psedg = sf(268, stash_sec) ! Graupel sedimentation rate
  casdiags % l_psedl = sf(269, stash_sec) ! Liquid sedimentation rate
  casdiags % l_phomr = sf(271, stash_sec) ! Homogeneous freezing of rain
  casdiags % l_pcond = sf(325, stash_sec) ! Condensation/Evaporation rate
  casdiags % l_psedi = sf(336, stash_sec) ! Ice cloud sublimation rate

  ! All the following diagnostics have units of [Number/kg/s].
  casdiags % l_nhomc = sf(350, stash_sec) ! Homo. Freezing of Cloud Number Tend.
  casdiags % l_nhomr = sf(351, stash_sec) ! Homo. Freezing of Rain Number Tend.
  casdiags % l_nihal = sf(352, stash_sec) ! Ice number tendency: Hallett-Mossop
  casdiags % l_ninuc = sf(353, stash_sec) ! Ice number tendency: Ice Nucleation
  casdiags % l_nsedi = sf(354, stash_sec) ! Ice number tendency: sedimentation
  casdiags % l_nseds = sf(355, stash_sec) ! Snow number tendency: sedimentation
  casdiags % l_nsedg = sf(356, stash_sec) ! Graupel number tend.: sedimentation

END IF ! model type /= mt_single_column.

CALL allocate_diagnostic_space(its, ite, jts, jte, kts, kte)

!----------------------------------------------------------------------------
! Allocate additional ice working variables if required by lightning scheme
! later
!----------------------------------------------------------------------------

IF (l_use_electric .AND. l_mcr_qgraup .AND. .NOT. l_casim_warm_only ) THEN

  ALLOCATE ( qcf_work(    tdims%i_start : tdims%i_end,                  &
                          tdims%j_start : tdims%j_end,                  &
                                      1 : tdims%k_end ) )

  ALLOCATE ( qcf2_work(   tdims%i_start : tdims%i_end,                  &
                          tdims%j_start : tdims%j_end,                  &
                                      1 : tdims%k_end ) )

  ALLOCATE ( qgraup_work(  tdims%i_start : tdims%i_end,                 &
                           tdims%j_start : tdims%j_end,                 &
                                       1 : tdims%k_end ) )

END IF

!-----------------------------------------------------------------------------
! Configure cloud number if required by the radiation scheme
!----------------------------------------------------------------------------

IF ( l_consistent_cdnc ) THEN
  IF ( l_mp_CloudNumber ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP SHARED( tdims, n_drop_pot, CloudNumber)                                 &
!$OMP PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          n_drop_pot(i,j,k) = CloudNumber(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF
END IF ! l_consistent_cdnc

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP SHARED( tdims, height, dAitkenSolMass, dAitkenSolNumber, dAccumSolMass, &
!$OMP        dAccumSolNumber, dCoarseSolMass, dCoarseSolNumber,               &
!$OMP        dActSolLiq_casim, dCoarseDustMass, dCoarseDustNumber,            &
!$OMP        dActInsolIce_casim, dActSolIce_casim, dActInsolLiq_casim,        &
!$OMP        dActSolNumber_casim, dActInsolNumber_casim, dActSolRain_casim,   &
!$OMP        ActSolLiq_casim, ActSolRain_casim, ActInsolIce_casim,            &
!$OMP        ActSolIce_casim, ActInsolLiq_casim, ActSolNumber_casim,          &
!$OMP        ActInsolNumber_casim  )                                          &
!$OMP PRIVATE( i, j, k )
DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end
    DO k = 1, tdims%k_end
      height(k,i,j)                = 0.0
      dAitkenSolMass(k,i,j)        = 0.0
      dAitkenSolNumber(k,i,j)      = 0.0
      dAccumSolMass(k,i,j)         = 0.0
      dAccumSolNumber(k,i,j)       = 0.0
      dCoarseSolMass(k,i,j)        = 0.0
      dCoarseSolNumber(k,i,j)      = 0.0
      dActSolLiq_casim(k,i,j)      = 0.0
      dCoarseDustMass(k,i,j)       = 0.0
      dCoarseDustNumber(k,i,j)     = 0.0
      dActInsolIce_casim(k,i,j)    = 0.0
      dActSolIce_casim(k,i,j)      = 0.0
      dActInsolLiq_casim(k,i,j)    = 0.0
      dActSolNumber_casim(k,i,j)   = 0.0
      dActInsolNumber_casim(k,i,j) = 0.0
      dActSolRain_casim(k,i,j)     = 0.0
      ActSolLiq_casim(k,i,j)       = 0.0
      ActSolRain_casim(k,i,j)      = 0.0
      ActInsolIce_casim(k,i,j)     = 0.0
      ActSolIce_casim(k,i,j)       = 0.0
      ActInsolLiq_casim(k,i,j)     = 0.0
      ActSolNumber_casim(k,i,j)    = 0.0
      ActInsolNumber_casim(k,i,j)  = 0.0
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

height_rho(:,:,:) = 0.0

!-------------------------------------------------------------------------
! Initialise large scale rain, snow and graupel rates to zero.
!-------------------------------------------------------------------------
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, ls_rain, ls_snow, ls_graup )                        &
!$OMP PRIVATE( i, j )
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      ls_rain(i,j)  = 0.0
      ls_snow(i,j)  = 0.0
      ls_graup(i,j) = 0.0
    END DO
  END DO
!$OMP END PARALLEL DO

!-------------------------------------------------------------------------
! Work out the air density for use in CASIM microphysics
!-------------------------------------------------------------------------

CALL mphys_air_density( rho_r2, q_n, qcl_n, qcf_n, qcf2_n, qrain_n, qgraup_n, &
                        rhodz_dry, rhodz_moist, deltaz )

IF (l_mr_physics) THEN

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, dz_casim, deltaz, rho_casim, rhodz_dry, rhodz_moist)&
!$OMP PRIVATE( i, j, k )
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dz_casim(k,i,j)  = deltaz(i,j,k)
        rho_casim(k,i,j) = rhodz_dry(i,j,k) / dz_casim(k,i,j)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

ELSE ! l_mr_physics

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, dz_casim, deltaz, rho_casim, rhodz_dry, rhodz_moist)&
!$OMP PRIVATE( i, j, k )
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dz_casim(k,i,j)  = deltaz(i,j,k)
        rho_casim(k,i,j) = rhodz_moist(i,j,k) / dz_casim(k,i,j)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

END IF ! l_mr_physics

!------------------------------------------------------------------------
! Include fixed aerosol number, if required
!------------------------------------------------------------------------
IF ( l_fix_aerosol .AND. casim_aerosol_option > 0 ) THEN
! no transport of aerosol,  set CASIM aerosol accumulation mode to a
! fixed value all others set to 0.0

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( tdims, AitkenSolMass, AitkenSolNumber, CoarseSolMass,    &
!$OMP         CoarseSolNumber, AccumDustMass, AccumDustNumber,         &
!$OMP         CoarseDustMass,  CoarseDustNumber,  AccumSolMass,        &
!$OMP         AccumSolNumber )                                         &
!$OMP PRIVATE( i, j, k )
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

        AitkenSolMass    (k,i,j)  = 0.0
        AitkenSolNumber  (k,i,j)  = 0.0
        CoarseSolMass    (k,i,j)  = 0.0
        CoarseSolNumber  (k,i,j)  = 0.0
        AccumDustMass    (k,i,j)  = 0.0
        AccumDustNumber  (k,i,j)  = 0.0
        CoarseDustMass   (k,i,j)  = 0.0
        CoarseDustNumber (k,i,j)  = 0.0
        AccumSolMass     (k,i,j)  = fixed_accum_sol_mass
        AccumSolNumber   (k,i,j)  = fixed_accum_sol_number

      END DO ! i
    END DO   ! j
  END DO     ! k
!$OMP END PARALLEL DO
ELSE IF (casim_aerosol_option == 0 ) THEN
! Set all aerosol species to zero (failsafe option)

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( tdims, AitkenSolMass, AitkenSolNumber, CoarseSolMass,    &
!$OMP         CoarseSolNumber, AccumDustMass, AccumDustNumber,         &
!$OMP         CoarseDustMass,  CoarseDustNumber,  AccumSolMass,        &
!$OMP         AccumSolNumber )                                         &
!$OMP PRIVATE( i, j, k )
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

        AitkenSolMass    (k,i,j)  = 0.0
        AitkenSolNumber  (k,i,j)  = 0.0
        CoarseSolMass    (k,i,j)  = 0.0
        CoarseSolNumber  (k,i,j)  = 0.0
        AccumDustMass    (k,i,j)  = 0.0
        AccumDustNumber  (k,i,j)  = 0.0
        CoarseDustMass   (k,i,j)  = 0.0
        CoarseDustNumber (k,i,j)  = 0.0
        AccumSolMass     (k,i,j)  = 0.0
        AccumSolNumber   (k,i,j)  = 0.0

      END DO ! i
    END DO   ! j
  END DO     ! k
END IF ! l_fix_aerosol

!------------------------------------------------------------------------
! Convert UM (i,j,k) dimensions to CASIM (k,i,j) dimensions for various
! input fields
!------------------------------------------------------------------------
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, th_casim, t_n, exner_theta_levels, dth_casim,       &
!$OMP         t_inc, pii_casim, p_casim, p_layer_centres, w_casim, w,    &
!$OMP         tke_casim, qv_casim, q_n, dqv_casim, q_inc, qc_casim,      &
!$OMP         qcl_n, dqc_casim, qcl_inc, qr_casim, qrain_n, dqr_casim,   &
!$OMP         qrain_inc )                                                &
!$OMP PRIVATE( i, j, k )
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      th_casim(k,i,j)  = t_n(i,j,k)/exner_theta_levels(i,j,k)
      dth_casim(k,i,j) = t_inc(i,j,k)/exner_theta_levels(i,j,k)
      pii_casim(k,i,j) = exner_theta_levels(i,j,k)
      p_casim(k,i,j)   = p_layer_centres(i,j,k)
      w_casim(k,i,j)   = w(i,j,k)
      tke_casim(k,i,j) = fixed_tke_value
      qv_casim(k,i,j)  = q_n(i,j,k)
      dqv_casim(k,i,j) = q_inc(i,j,k)
      qc_casim(k,i,j)  = qcl_n(i,j,k)
      dqc_casim(k,i,j) = qcl_inc(i,j,k)
      qr_casim(k,i,j)  = qrain_n(i,j,k)
      dqr_casim(k,i,j) = qrain_inc(i,j,k)
    END DO ! i
  END DO   ! j
END DO     ! k
!$OMP END PARALLEL DO

!-----------------------------------------------------------------------
! Convert additional cloud number, if present
!-----------------------------------------------------------------------
IF (l_mp_CloudNumber) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, nc_casim, CloudNumber, dnc_casim, CloudNumber_inc ) &
!$OMP PRIVATE( i, j, k )
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        nc_casim(k,i,j)  = CloudNumber(i,j,k)
        dnc_casim(k,i,j) = CloudNumber_inc(i,j,k)
      END DO ! i
    END DO   ! j
  END DO     ! k
!$OMP END PARALLEL DO
END IF ! l_mp_CloudNumber

!-----------------------------------------------------------------------
! Convert additional rain moments, if present
!-----------------------------------------------------------------------
IF ( rain_mom == triple_moment ) THEN

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, nr_casim,  RainNumber, dnr_casim,  RainNumber_inc,  &
!$OMP         m3r_casim, Rain3mom, dm3r_casim, Rain3mom_inc )            &
!$OMP PRIVATE( i, j, k )
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        nr_casim(k,i,j)   = RainNumber(i,j,k)
        dnr_casim(k,i,j)  = RainNumber_inc(i,j,k)
        m3r_casim(k,i,j)  = Rain3mom(i,j,k)
        dm3r_casim(k,i,j) = Rain3mom_inc(i,j,k)
      END DO ! i
    END DO   ! j
  END DO     ! k
!$OMP END PARALLEL DO

ELSE IF ( rain_mom == double_moment ) THEN

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, nr_casim,  RainNumber, dnr_casim,  RainNumber_inc ) &
!$OMP PRIVATE( i, j, k )
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        nr_casim(k,i,j)  = RainNumber(i,j,k)
        dnr_casim(k,i,j) = RainNumber_inc(i,j,k)
      END DO ! i
    END DO   ! j
  END DO     ! k
!$OMP END PARALLEL DO

END IF ! rain_mom (double or triple)

!-----------------------------------------------------------------------
! Convert ice phase variables to CASIM dimensions
!-----------------------------------------------------------------------
IF (.NOT. l_casim_warm_only) THEN

! First do the ice and snow which are always present
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, qi_casim, qs_casim, qcf_n, qcf2_n, dqi_casim,       &
!$OMP         qcf2_inc, dqs_casim, qcf_inc )                             &
!$OMP PRIVATE( i, j, k )
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        qi_casim(k,i,j)  = qcf2_n(i,j,k) 
        qs_casim(k,i,j)  = qcf_n(i,j,k)
        dqi_casim(k,i,j) = qcf2_inc(i,j,k)
        dqs_casim(k,i,j) = qcf_inc(i,j,k)
      END DO ! i
    END DO   ! j
  END DO     ! k
!$OMP END PARALLEL DO

!-----------------------------------------------------------------------
! Convert additional ice and snow moments, if present
!-----------------------------------------------------------------------
  IF (l_mp_IceNumber ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, ni_casim, IceNumber, dni_casim, IceNumber_inc )     &
!$OMP PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          ni_casim(k,i,j)  = IceNumber(i,j,k)
          dni_casim(k,i,j) = IceNumber_inc(i,j,k)
        END DO ! i
      END DO   ! j
    END DO     ! k
!$OMP END PARALLEL DO
  END IF

  IF ( snow_mom == triple_moment ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, ns_casim, SnowNumber, m3s_casim, Snow3mom,          &
!$OMP         dns_casim, SnowNumber_inc, dm3s_casim, Snow3mom_inc )      &
!$OMP PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          ns_casim(k,i,j)   = SnowNumber(i,j,k)
          m3s_casim(k,i,j)  = Snow3mom(i,j,k)
          dns_casim(k,i,j)  = SnowNumber_inc(i,j,k)
          dm3s_casim(k,i,j) = Snow3mom_inc(i,j,k)
        END DO ! i
      END DO   ! j
    END DO     ! k
!$OMP END PARALLEL DO

  ELSE IF ( snow_mom == double_moment ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, ns_casim, SnowNumber, dns_casim, SnowNumber_inc  )  &
!$OMP PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          ns_casim(k,i,j)   = SnowNumber(i,j,k)
          dns_casim(k,i,j)  = SnowNumber_inc(i,j,k)
        END DO ! i
      END DO   ! j
    END DO     ! k
!$OMP END PARALLEL DO

  END IF ! snow_mom (double or triple)

  IF ( graup_mom == triple_moment ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, qg_casim, qgraup_n, ng_casim, GraupNumber,          &
!$OMP         dqg_casim, qgraup_inc, dng_casim, GraupNumber_inc,         &
!$OMP         m3g_casim, Graup3mom, dm3g_casim, Graup3mom_inc )          &
!$OMP PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qg_casim(k,i,j)   = qgraup_n(i,j,k)
          ng_casim(k,i,j)   = GraupNumber(i,j,k)
          m3g_casim(k,i,j)  = Graup3mom(i,j,k)
          dqg_casim(k,i,j)  = qgraup_inc(i,j,k)
          dng_casim(k,i,j)  = GraupNumber_inc(i,j,k)
          dm3g_casim(k,i,j) = Graup3mom_inc(i,j,k)
        END DO ! i
      END DO   ! j
    END DO     ! k
!$OMP END PARALLEL DO
  ELSE IF ( graup_mom == double_moment ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, qg_casim, qgraup_n, ng_casim, GraupNumber,          &
!$OMP         dqg_casim, qgraup_inc, dng_casim, GraupNumber_inc )        &
!$OMP PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qg_casim(k,i,j)   = qgraup_n(i,j,k)
          ng_casim(k,i,j)   = GraupNumber(i,j,k)
          dqg_casim(k,i,j)  = qgraup_inc(i,j,k)
          dng_casim(k,i,j)  = GraupNumber_inc(i,j,k)
        END DO ! i
      END DO   ! j
    END DO     ! k
!$OMP END PARALLEL DO
  ELSE IF ( graup_mom == single_moment ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, qg_casim, qgraup_n, dqg_casim, qgraup_inc )         &
!$OMP PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qg_casim(k,i,j)   = qgraup_n(i,j,k)
          dqg_casim(k,i,j)  = qgraup_inc(i,j,k)
        END DO ! i
      END DO   ! j
    END DO     ! k
!$OMP END PARALLEL DO
  END IF ! graup_mom (double or triple)

END IF ! .NOT. l_casim_warm_only

!------------------------------------------------------------------------------
! For aerosol processing runs, produce aerosol processing prognostics on
! CASIM dimensions.
!------------------------------------------------------------------------------

IF ( casim_aerosol_process_level > no_processing ) THEN

  IF ( l_mp_ActiveSolLiquid ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, ActSolLiq_casim, ActiveSolLiquid )                  &
!$OMP PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          ActSolLiq_casim(k,i,j) = ActiveSolLiquid(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF ! l_mp_ActiveSolLiquid

  IF ( l_mp_ActiveSolRain ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, ActSolRain_casim, ActiveSolRain )                   &
!$OMP PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          ActSolRain_casim(k,i,j) = ActiveSolRain(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  IF ( l_mp_ActiveInsolIce ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, ActInsolIce_casim, ActiveInsolIce )                 &
!$OMP PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          ActInsolIce_casim(k,i,j) = ActiveInsolIce(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  IF ( l_mp_ActiveSolIce ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, ActSolIce_casim, ActiveSolIce )                     &
!$OMP PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
            ActSolIce_casim(k,i,j) = ActiveSolIce(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  IF ( l_mp_ActiveInsolLiquid ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, ActInsolLiq_casim, ActiveInsolLiquid )              &
!$OMP PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
            ActInsolLiq_casim(k,i,j) = ActiveInsolLiquid(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  IF ( l_mp_ActiveSolNumber ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, ActSolNumber_casim, ActiveSolNumber )               &
!$OMP PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          ActSolNumber_casim(k,i,j) = ActiveSolNumber(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  IF ( l_mp_ActiveInSolNumber ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, ActInsolNumber_casim, ActiveInsolNumber )           &
!$OMP PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
            ActInsolNumber_casim(k,i,j)  = ActiveInsolNumber(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
   END IF

 END IF  ! casim_aerosol_process_level > no_processing

! ------------------------------------------------------------------------------
! Call the Shipway Microphysics code on CASIM side
! ------------------------------------------------------------------------------

  CALL shipway_microphysics( its, ite, jts, jte, kts, kte,  timestep,          &
                             qv_casim, qc_casim, qr_casim, nc_casim, nr_casim, &
                             m3r_casim, qi_casim, qs_casim, qg_casim, ni_casim,&
                             ns_casim, ng_casim, m3s_casim, m3g_casim,         &
                             th_casim,                                         &
                             AitkenSolMass, AitkenSolNumber, AccumSolMass,     &
                             AccumSolNumber, CoarseSolMass, CoarseSolNumber,   &
                             ActSolLiq_casim, ActSolRain_casim, CoarseDustMass,&
                             CoarseDustNumber, ActInsolIce_casim,              &
                             ActSolIce_casim, ActInsolLiq_casim, AccumDustMass,&
                             AccumDustNumber, ActSolNumber_casim,              &
                             ActInsolNumber_casim, pii_casim, p_casim,         &
                             rho_casim, w_casim, tke_casim, height_rho,        &
                             height, dz_casim,                                 &
!                   input variables above  || in/out variables below 
                             dqv_casim, dqc_casim,  dqr_casim, dnc_casim,      &
                             dnr_casim, dm3r_casim, dqi_casim, dqs_casim,      &
                             dqg_casim, dni_casim, dns_casim,  dng_casim,      &
                             dm3s_casim, dm3g_casim, dth_casim,                &
                             dAitkenSolMass, dAitkenSolNumber, dAccumSolMass,  &
                             dAccumSolNumber, dCoarseSolMass, dCoarseSolNumber,&
                             dActSolLiq_casim,   dActSolRain_casim,            &
                             dCoarseDustMass,    dCoarseDustNumber,            &
                             dActInsolIce_casim, dActSolIce_casim,             &
                             dActInsolLiq_casim, dAccumDustMass,               &
                             dAccumDustNumber,   dActSolNumber_casim,          &
                             dActInsolNumber_casim,                            &
                             ils, ile,  jls, jle, kls, kle, l_tendency=.FALSE. )

!------------------------------------------------------------------------------
! Calculation of increment fields from CASIM dimensions back to UM dimensions
!------------------------------------------------------------------------------

IF ( rain_mom == triple_moment ) THEN

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                          &
!$OMP SHARED( tdims, RainNumber_inc, dnr_casim, Rain3mom_inc, dm3r_casim, &
!$OMP         t_inc, dth_casim, exner_theta_levels, q_inc, dqv_casim,     &
!$OMP         qcl_inc, dqc_casim, qrain_inc, dqr_casim, t_work, t_n  )    &
!$OMP PRIVATE( i, j, k )
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        RainNumber_inc(i,j,k)  = RainNumber_inc(i,j,k) + dnr_casim(k,i,j)
        Rain3mom_inc(i,j,k)    = Rain3mom_inc(i,j,k)   + dm3r_casim(k,i,j)
        t_inc(i,j,k)           = t_inc(i,j,k)          + dth_casim(k,i,j) *   &
                                 exner_theta_levels(i,j,k)

        t_work(i,j,k)    = t_n(i,j,k)       + t_inc(i,j,k)
        q_inc(i,j,k)     = q_inc(i,j,k)     + dqv_casim(k,i,j)
        qcl_inc(i,j,k)   = qcl_inc(i,j,k)   + dqc_casim(k,i,j)
        qrain_inc(i,j,k) = qrain_inc(i,j,k) + dqr_casim(k,i,j)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

ELSE IF ( rain_mom == double_moment ) THEN

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                          &
!$OMP SHARED( tdims, RainNumber_inc, dnr_casim, t_inc, dth_casim,         &
!$OMP         exner_theta_levels, q_inc, dqv_casim, qcl_inc, dqc_casim,   &
!$OMP         qrain_inc, dqr_casim, t_work, t_n )                         &
!$OMP PRIVATE( i, j, k )
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        RainNumber_inc(i,j,k)  = RainNumber_inc(i,j,k) + dnr_casim(k,i,j)
        t_inc(i,j,k)           = t_inc(i,j,k)          + dth_casim(k,i,j) *   &
                                 exner_theta_levels(i,j,k)

        t_work(i,j,k)    = t_n(i,j,k)       + t_inc(i,j,k)
        q_inc(i,j,k)     = q_inc(i,j,k)     + dqv_casim(k,i,j)
        qcl_inc(i,j,k)   = qcl_inc(i,j,k)   + dqc_casim(k,i,j)
        qrain_inc(i,j,k) = qrain_inc(i,j,k) + dqr_casim(k,i,j)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

ELSE IF ( rain_mom == single_moment ) THEN

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, t_inc, dth_casim, exner_theta_levels, q_inc,        &
!$OMP         dqv_casim, qcl_inc, dqc_casim, qrain_inc, dqr_casim,       &
!$OMP         t_work, t_n )                                              &
!$OMP PRIVATE( i, j, k )
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        t_inc(i,j,k)     = t_inc(i,j,k)     + dth_casim(k,i,j) *              &
                           exner_theta_levels(i,j,k)

        t_work(i,j,k)    = t_n(i,j,k)       + t_inc(i,j,k)
        q_inc(i,j,k)     = q_inc(i,j,k)     + dqv_casim(k,i,j)
        qcl_inc(i,j,k)   = qcl_inc(i,j,k)   + dqc_casim(k,i,j)
        qrain_inc(i,j,k) = qrain_inc(i,j,k) + dqr_casim(k,i,j)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

END IF ! rain_mom

IF ( l_mp_CloudNumber ) THEN

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, CloudNumber_inc, dnc_casim )                        &
!$OMP PRIVATE( i, j, k )
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        CloudNumber_inc(i,j,k) = CloudNumber_inc(i,j,k) + dnc_casim(k,i,j)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

END IF ! l_mp_CloudNumber

!------------------------------------------------------------------------------
! Calculation of ice microphysics increment fields from CASIM dimensions back 
! to UM dimensions
!------------------------------------------------------------------------------

IF ( .NOT. l_casim_warm_only ) THEN

  IF ( snow_mom == triple_moment ) THEN

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, qcf_inc, dqs_casim, qcf2_inc, dqi_casim,            &
!$OMP         SnowNumber_inc, dns_casim, Snow3mom_inc, dm3s_casim )      &
!$OMP PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qcf_inc(i,j,k)        = qcf_inc(i,j,k)        + dqs_casim(k,i,j)
          qcf2_inc(i,j,k)       = qcf2_inc(i,j,k)       + dqi_casim(k,i,j)
          SnowNumber_inc(i,j,k) = SnowNumber_inc(i,j,k) + dns_casim(k,i,j)
          Snow3mom_inc(i,j,k)   = Snow3mom_inc(i,j,k)   + dm3s_casim(k,i,j)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

  ELSE IF ( snow_mom == double_moment ) THEN

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, qcf_inc, dqs_casim, qcf2_inc, dqi_casim,            &
!$OMP         SnowNumber_inc, dns_casim )                                &
!$OMP PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qcf_inc(i,j,k)        = qcf_inc(i,j,k)        + dqs_casim(k,i,j)
          qcf2_inc(i,j,k)       = qcf2_inc(i,j,k)       + dqi_casim(k,i,j)
          SnowNumber_inc(i,j,k) = SnowNumber_inc(i,j,k) + dns_casim(k,i,j)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

  ELSE IF ( snow_mom == single_moment ) THEN

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, qcf_inc, dqs_casim, qcf2_inc, dqi_casim )           &
!$OMP PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qcf_inc(i,j,k)  = qcf_inc(i,j,k)   + dqs_casim(k,i,j)
          qcf2_inc(i,j,k) = qcf2_inc(i,j,k)  + dqi_casim(k,i,j)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

  END IF ! snow_mom

  IF ( l_mp_IceNumber ) THEN

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, IceNumber_inc, dni_casim )                          &
!$OMP PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          IceNumber_inc(i,j,k) = IceNumber_inc(i,j,k) + dni_casim(k,i,j)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

  END IF ! l_mp_icenumber

  IF ( graup_mom == triple_moment ) THEN

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, qgraup_inc, dqg_casim, GraupNumber_inc, dng_casim,  &
!$OMP         Graup3mom_inc, dm3g_casim )                                &
!$OMP PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qgraup_inc(i,j,k)      = qgraup_inc(i,j,k)      + dqg_casim(k,i,j)
          GraupNumber_inc(i,j,k) = GraupNumber_inc(i,j,k) + dng_casim(k,i,j)
          Graup3mom_inc(i,j,k)   = Graup3mom_inc(i,j,k)   + dm3g_casim(k,i,j)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

  ELSE IF ( graup_mom == double_moment ) THEN

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, qgraup_inc, dqg_casim, GraupNumber_inc, dng_casim ) &
!$OMP PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qgraup_inc(i,j,k)      = qgraup_inc(i,j,k)      + dqg_casim(k,i,j)
          GraupNumber_inc(i,j,k) = GraupNumber_inc(i,j,k) + dng_casim(k,i,j)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

  ELSE IF ( graup_mom == single_moment ) THEN

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, qgraup_inc, dqg_casim )                             &
!$OMP PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qgraup_inc(i,j,k) = qgraup_inc(i,j,k) + dqg_casim(k,i,j)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

  END IF ! graup_mom

END IF ! l_casim_warm_only

!------------------------------------------------------------------------------
! Pass back increments in Aerosol processing switches from CASIM dimensions
! to UM dimensions
!------------------------------------------------------------------------------

IF ( casim_aerosol_process_level > no_processing ) THEN

  IF ( l_mp_ActiveSolLiquid ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, dActiveSolLiquid, dActSolLiq_casim )                &
!$OMP PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          dActiveSolLiquid(i,j,k) = dActiveSolLiquid(i,j,k) +                  &
                                    dActSolLiq_casim(k,i,j)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF ! l_mp_ActiveSolLiquid

  IF ( l_mp_ActiveInsolIce ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, dActiveInsolIce, dActInsolIce_casim )               &
!$OMP PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          dActiveInsolIce(i,j,k) = dActiveInsolIce(i,j,k) +                    &
                                   dActInsolIce_casim(k,i,j)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF ! l_mp_ActiveSolLiquid

  IF ( l_mp_ActiveSolIce ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, dActiveSolIce, dActSolIce_casim )                   &
!$OMP PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          dActiveSolIce(i,j,k) = dActiveSolIce(i,j,k) +                        &
                                 dActSolIce_casim(k,i,j)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF ! l_mp_ActiveSolIce

  IF ( l_mp_ActiveInsolLiquid ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, dActiveInsolLiquid, dActInsolLiq_casim )            &
!$OMP PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          dActiveInsolLiquid(i,j,k) = dActiveInsolLiquid(i,j,k) +              &
                                      dActInsolLiq_casim(k,i,j)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF ! l_mp_ActiveInsolLiquid

  IF ( l_mp_ActiveSolRain ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, dActiveSolRain, dActSolRain_casim )                 &
!$OMP PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          dActiveSolRain(i,j,k) = dActiveSolRain(i,j,k) +                      &
                                  dActSolRain_casim(k,i,j)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF ! l_mp_ActiveSolRain

  IF ( l_mp_ActiveInsolNumber ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, dActiveInsolNumber, dActInsolNumber_casim )         &
!$OMP PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          dActiveInsolNumber(i,j,k) = dActiveInsolNumber(i,j,k) +              &
                                      dActInsolNumber_casim(k,i,j)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF ! l_mp_ActiveInsolNumber

  IF ( l_mp_ActiveSolNumber ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, dActiveSolNumber, dActSolNumber_casim )             &
!$OMP PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end        
          dActiveSolNumber(i,j,k) = dActiveSolNumber(i,j,k) +                  &
                                    dActSolNumber_casim(k,i,j)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF ! l_mp_ActiveSolNumber

END IF ! casim_aerosol_process_level > no_processing

!------------------------------------------------------------------------------
! Write precipitation diagnostics into quantities to send to other UM sections
!------------------------------------------------------------------------------
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP SHARED( tdims, ls_rain, ls_snow, ls_graup, casdiags )              &
!$OMP PRIVATE( i, j )
DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end
    ls_rain(i,j)  = casdiags % SurfaceRainR(i,j)
    ls_snow(i,j)  = casdiags % SurfaceSnowR(i,j)
    ls_graup(i,j) = casdiags % SurfaceGraupR(i,j)
  END DO
END DO
!$OMP END PARALLEL DO

!------------------------------------------------------------------------------
! Write micro_tends out for the boundary layer section 
! These are microphys increments to TL (micro_tends (:,:,1,:)) and 
! QW (micro_tends (:,:,2,:)) .
!------------------------------------------------------------------------------
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                           &
!$OMP SHARED( bl_levels, tdims, micro_tends, t_inc, lcrcp, lsrcp, qcl_inc, &
!$OMP         qrain_inc, qcf_inc, qcf2_inc, qgraup_inc, q_inc,             &
!$OMP         recip_timestep )                                             &
!$OMP PRIVATE( i, j, k )
DO k = 1, bl_levels
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end

       ! Recall that:
       ! lcrcp = lc/cp
       ! lsrcp = (lc+lf)/cp
       ! Micro tends (:,:,1,:) is t_inc - ( lcrcp * liquid increments ) 
       !                                - ( lsrcp * ice increments    )
       micro_tends(i,j,1,k) = t_inc(i,j,k)                                     &
                            - ( lcrcp * ( qcl_inc(i,j,k) + qrain_inc(i,j,k)) ) &
                            - ( lsrcp * ( qcf_inc(i,j,k) + qcf2_inc(i,j,k) +   &
                                          qgraup_inc(i,j,k) ) )


       ! Micro tends (:,:,2,:) is all moist mass variable increments
       micro_tends(i,j,2,k) = q_inc(i,j,k)     + qcl_inc(i,j,k)    +           &
                              qrain_inc(i,j,k) + qcf_inc(i,j,k)    +           &
                              qcf2_inc(i,j,k)  + qgraup_inc(i,j,k)

       ! Each tendency needs to be divided by timestep length in order to 
       ! convert to a rate. Multiplication by recip_timestep is quicker.
       micro_tends(i,j,1,k) = micro_tends(i,j,1,k) * recip_timestep
       micro_tends(i,j,2,k) = micro_tends(i,j,2,k) * recip_timestep

    END DO
  END DO
END DO
!$OMP END PARALLEL DO

!------------------------------------------------------------------------------
! Produce Section 4 diagnostics appropriately, if requested
!------------------------------------------------------------------------------
IF (model_type /= mt_single_column .AND. sf(izero, stash_sec) ) THEN

  CALL diagnostics_casim( stashwork4, t_n, q_n, qcl_n, qcf2_n, t_inc, q_inc,  &
                          qcl_inc, qcf_inc, qrain_inc, qgraup_inc, qcf2_inc,  &
                          p_layer_centres, deltaz )

END IF ! model_type and sf

!------------------------------------------------------------------------------
! Call Lightning scheme, if required and prognostic graupel is in use
! with ice cloud in use as well
!------------------------------------------------------------------------------
IF (l_use_electric .AND. l_mcr_qgraup .AND. .NOT. l_casim_warm_only) THEN

  ! Need to calculate working qcf, qcf2 and qgraup based on updates
  ! from CASIM. Working temperature already calculated above.

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                            &
!$OMP SHARED( tdims, qcf_work, qcf_n, qcf_inc, qcf2_work, qcf2_n, qcf2_inc, &
!$OMP         qgraup_work, qgraup_n, qgraup_inc )                           &
!$OMP PRIVATE( i, j, k )
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        qcf_work(i,j,k)    = qcf_n(i,j,k)    +  qcf_inc(i,j,k)
        qcf2_work(i,j,k)   = qcf2_n(i,j,k)   +  qcf2_inc(i,j,k)
        qgraup_work(i,j,k) = qgraup_n(i,j,k) +  qgraup_inc(i,j,k)
      END DO ! i
    END DO ! j
  END DO ! k
!$OMP END PARALLEL DO    

! -------------------------------------------------------------------
! Call electric main in section atmosphere/electric
! ------------------------------------------------------------------

  CALL electric_main( qcf_work, qcf2_work, qgraup_work, rhodz_dry,        &
                      rhodz_moist, t_work, w, at_extremity, stashwork21,  &
                      flash_pot(:,:,1 : tdims%k_end ) )

  ! Deallocate electric-specific working arrays.
  DEALLOCATE ( qgraup_work )
  DEALLOCATE ( qcf2_work   )
  DEALLOCATE ( qcf_work    )

END IF ! l_use_electric ( lightning scheme )

!-------------------------------------------------------------------
! Final deallocations of any diagnostic space used.
!-------------------------------------------------------------------
CALL deallocate_diagnostic_space()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE casim_ctl

END MODULE casim_ctl_mod
