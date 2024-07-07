! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Interface to Atmos Physics parametrizations before S-L advection.
!
! Subroutine Interface:
SUBROUTINE Atmos_Physics1(                                              &

! Parallel variables
        global_row_length, global_rows, n_proc, n_procx, n_procy        &
      , g_rows, g_row_length                                            &

! model dimensions
      , row_length, rows, n_rows, land_points                           &
      , bl_levels, dst_levels, dsm_levels, Ozone_levels, cloud_levels   &
      , land_ice_points, soil_points, n_cca_levels, ntiles              &
      , salt_dim1, salt_dim2, salt_dim3, tr_levels, tr_ukca             &
      , cdnc_dim1, cdnc_dim2, cdnc_dim3                                 &
      , co2_dim_len, co2_dim_row, co2_dim_lev                           &
      , n_arcl_species, n_arcl_compnts, i_arcl_compnts                  &

! model switches
      , l_lbc_old                                                       &
      , l_ukca_chem, l_ukca_set_trace_gases                             &
      , l_ukca_strat, l_ukca_strattrop                                  &
      , l_ukca_prescribech4, l_use_arcl                                 &

! model Parameters
      , rhcrit                                                          &
      , min_trop_level, max_trop_level                                  &

! Position of greenhouse gases in tracer_ukca array
      , ngrgas, grgas_addr                                              &
! parameter for stochastic physics random parameters
      , m_ci                                                            &

! in coordinate information
      , delta_lambda, delta_phi, lat_rot_NP, long_rot_NP                &

! in time stepping information.
      , val_year, val_day_number, val_hour, val_minute                  &
      , val_second, previous_time                                       &

! diagnostic info
      , STASHwork1,STASHwork2,STASHwork4,STASHwork6,STASHwork14         &
      , STASHwork21                                                     &
!
! Additional variables for SCM diagnostics
      , nSCMDpkgs, L_SCMDiags                                           &
!
! in data fields.
      , theta, q, qcl, qcf, qcf2, qrain, qgraup, rho, u, v, w           &
      , p, p_star, exner_rho_levels, exner_theta_levels, land_sea_mask  &
      , p_theta_levels, fland_ctile, frac_control, ukca_cdnc            &

! ancillary fields and fields needed to be kept from timestep to
! timestep

      , land_index,rgrain,soot,canht,ntml,cumulus                       &
      , ice_fract,ice_fract_cat,ice_thick_cat                           &
      , cca_dp, cca_md, cca_sh                                          &
      , cca, ccb, cct, cclwp, ccw, lcbase, totalppn                     &
      , t_surf, tstar_land_ctile, tstar_sea_ctile, tstar_sice_ctile     &
      , sice_alb_ctile,land_alb_ctile,snow_depth,snow_depth_sea_cat     &
      , ozone, SW_incs, LW_incs, dirpar_inc                             &
      , O3_trop_level, O3_trop_height, T_trop_level, T_trop_height      &
      , zh, sd_orog_land , orog_grad_xx_land, orog_grad_xy_land         &
      , orog_grad_yy_land, area_cloud_fraction, cf, cfl, cff            &
      , aerosol_em, arcl                                                &
      , albsoil, albobs_sw, albobs_vis, albobs_nir, lai, snow_tile      &
      , tile_frac, tstar_tile, z0_tile                                  &
      , dOLR_rts, LW_down, SW_tile_rts, es_space_interp, rad_mask       &
      , cos_zenith_angle                                                &
      , easyaerosol_sw, easyaerosol_lw, easyaerosol_cdnc                &

! Variables for COSP
      , cosp_crain_3d,cosp_csnow_3d                                     &

! IN/OUT JULES prognostics
      , snowdepth_p,lake_h_ice_p, z0m_sea, chloro_sea                   &
! In variable storing BL w-variance diagnostic
      , bl_w_var                                                        &
! in/out
      , theta_star, q_star, qcl_star, qcf_star, qcf2_star, qrain_star   &
      , qgraup_star, cf_star, cfl_star, cff_star, u_inc, v_inc          &
      , energy_correction, sum_eng_fluxes, sum_moist_flux, aerosol      &
      , flash_pot                                                       &
      , dust_div1,dust_div2,dust_div3,dust_div4,dust_div5,dust_div6     &
      , so2, so4_aitken, so4_accu, so4_diss, nh3                        &
      , soot_new, soot_agd, soot_cld, bmass_new, bmass_agd, bmass_cld   &
      , ocff_new, ocff_agd, ocff_cld, nitr_acc, nitr_diss               &
      , co2, tracer, tracer_ukca, biogenic, asteps_since_triffid        &
      , ukca_radaer                                                     &
! out fields
      , ls_rain, ls_rainfrac, ls_snow, ls_graup, micro_tends            &
      , unscaled_dry_rho, photosynth_act_rad, rad_hr, dOLR, SW_tile     &

! error information
      , Error_code  )

USE swap_bounds_mv_mod,   ONLY: swap_bounds_mv
USE dynamics_input_mod,   ONLY: l_check_moist_inc
USE dynamics_testing_mod, ONLY: l_idealised_data

USE planet_constants_mod, ONLY: kappa, p_zero, cp

USE water_constants_mod, ONLY: lc, lf, tfs

USE atm_fields_bounds_mod, ONLY: pdims, pdims_s, tdims, tdims_s,       &
                                 tdims_l, udims, udims_s, vdims,       &
                                 vdims_s, wdims_s

USE cv_run_mod, ONLY: l_pc2_diag_sh

USE cloud_inputs_mod,  ONLY: l_micro_eros,  i_rhcpt,                    &
                             falliceshear_method, i_cld_vn
USE pc2_constants_mod, ONLY: real_shear, rhcpt_horiz_var, i_cld_pc2

USE ukca_trace_gas_mixratio, ONLY: ukca_set_trace_gas_mixratio

USE timestep_mod, ONLY: timestep_number, timestep, recip_timestep

USE rad_input_mod

USE nlstcall_mod, ONLY: ltimer

! Get flags for whether this is a prog / diag radiation timestep
USE set_rad_steps_mod, ONLY: l_rad_step_prog, l_rad_step_diag

USE clmchfcg_scenario_mod, ONLY: clim_fcg_nyears, clim_fcg_years,       &
                                 clim_fcg_levls, clim_fcg_rates,        &
     s_co2, s_n2o, s_ch4, s_cfc11, s_cfc12, s_cfc113, s_cfc114,         &
     s_hcfc22, s_hfc125, s_hfc134a, lenscen
! histories/scenarios of climate change forcings

USE arcl_mod, ONLY: npd_arcl_species, npd_arcl_compnts

USE NI_gwd_ctl_mod, ONLY: NI_gwd_ctl

USE physics_tendencies_mod,  ONLY: init_slowphys_tendencies,            &
                                   dtheta_ph1,                          &
                                   l_retain_slow_tendencies,            &
                                   l_retain_ph1_tendencies

USE trignometric_mod, ONLY:                                             &
    cos_theta_latitude, sin_theta_latitude,                             &
    cos_theta_longitude, sin_theta_longitude,                           &
    FV_cos_theta_latitude, true_latitude, true_longitude

USE global_2d_sums_mod, ONLY: global_2d_sums

USE g_wave_input_mod, ONLY: l_gwd, l_use_ussp

! JULES
USE jules_surface_mod, ONLY: l_flake_model
USE jules_sea_seaice_mod, ONLY: nice, nice_use, l_ctile
USE jules_radiation_mod, ONLY: l_snow_albedo
USE jules_vegetation_mod, ONLY: l_triffid
USE prognostics, ONLY:                                                  &
  snowdepth_surft
USE lake_mod, ONLY:                                                     &
  lake_h_ice_gb

USE ancil_info, ONLY:                                                   &
    ssi_pts, sea_pts, sice_pts                                          &
  , sice_pts_ncat, ssi_index, sea_index, sice_index                     &
  , sice_index_ncat, fssi_ij, sea_frac, sice_frac, sice_frac_ncat

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ukca_radaer_struct_mod
USE swapable_field_mod,    ONLY: swapable_field_pointer_type

USE level_heights_mod, ONLY:                                            &
    eta_theta_levels, eta_rho_levels,                                   &
    r_theta_levels, r_rho_levels

USE jules_surface_types_mod
USE stash_array_mod, ONLY: sf_calc, sf
USE ereport_mod, ONLY: ereport
USE UM_ParVars
USE cosp_types_mod, ONLY: cosp_config, cosp_subgrid, cosp_sghydro
USE cosp_init_mod
USE cosp_main_mod
USE ukca_feedback_mod, ONLY: p_o3, p_ch4, p_n2o, p_f11, p_f12,          &
                             p_f113, p_f22, p_h2os

USE uv_p_pnts_mod,   ONLY: uv_p_pnts, uv_p_pnts_halo
USE set_seasalt_mod, ONLY: set_seasalt_4A


USE add_eng_corr_mod, ONLY: add_eng_corr
USE dust_parameters_mod, ONLY: l_dust

USE carbon_options_mod, ONLY: l_co2_interactive
USE gen_phys_inputs_mod, ONLY: l_use_methox, l_mr_physics
USE ukca_option_mod, ONLY:                                              &
     l_ukca, l_ukca_aie1, l_ukca_radaer, l_ukca_radaer_sustrat,         &
     i_ukca_scenario
USE glomap_clim_option_mod, ONLY: &
    l_glomap_clim_radaer,         &
    l_glomap_clim_radaer_sustrat, &
    l_glomap_clim_aie1
USE run_aerosol_mod, ONLY: l_sulpc_so2, l_soot, l_ocff, l_biomass,      &
                           l_sulpc_nh3, l_nitrate
USE ukca_scenario_ctl_mod, ONLY: i_ukca_scenario_um
USE mphys_inputs_mod, ONLY: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain,      &
                            l_rain, l_use_seasalt_autoconv, l_casim

USE easyaerosol_option_mod, ONLY:                                       &
     l_easyaerosol_sw, l_easyaerosol_lw, l_easyaerosol_cdnc,            &
     l_easyaerosol_autoconv
USE def_easyaerosol, ONLY: t_easyaerosol_rad, t_easyaerosol_cdnc

USE casim_prognostics, ONLY: casim_alloc_increments, casim_prognostics_update

USE eng_corr_inputs_mod, ONLY: l_emcorr, lflux_reset
USE murk_inputs_mod, ONLY: l_murk, l_murk_source, l_murk_bdry,          &
                           l_murk_rad, murk_source_scale

USE cosp_input_mod, ONLY: l_cosp

USE diagnostics_pc2checks_mod, ONLY: diagnostics_pc2checks
USE u_to_p_mod, ONLY: u_to_p
USE v_to_p_mod, ONLY: v_to_p
USE trsrce_mod, ONLY: trsrce
USE check_dmoist_inc_mod, ONLY: check_dmoist_inc
USE set_fsd_parameters_mod, ONLY: set_fsd_parameters
USE umPrintMgr, ONLY: umPrint, umMessage

USE cosp_variable_mod, ONLY: cosp_gbx
USE null_cosp_gridbox_mod, ONLY: allocate_null_gbx, free_null_gbx


USE s_scmop_mod,   ONLY: default_streams                                &
  , t_avg, d_wet, d_all, scmdiag_pc2
USE scmoutput_mod, ONLY: scmoutput

USE errormessagelength_mod, ONLY: errormessagelength

USE model_domain_mod, ONLY: model_type, mt_global, mt_lam, mt_single_column

USE nlsizes_namelist_mod, ONLY: model_levels, tr_vars, super_array_size
USE swap_bounds_2d_mv_mod, ONLY: swap_bounds_2d_mv
USE mpp_conf_mod,  ONLY: swap_field_is_scalar

! CASIM
USE casim_ctl_mod, ONLY: casim_ctl

USE flux_diag_mod, ONLY: flux_diag
USE close_cloud_gen_mod, ONLY: close_cloud_gen
USE gas_calc_mod, ONLY: gas_calc
USE open_cloud_gen_mod, ONLY: open_cloud_gen
USE rad_ctl_mod, ONLY: rad_ctl
USE tropin_mod, ONLY: tropin
USE pc2_checks_mod, ONLY: pc2_checks

USE atm_step_local, ONLY: L_print_L2norms, L_tracer
USE atmos_phys1_norm_mod
USE atmos_start_norm_mod
USE sl_tracer_norm_mod
USE turb_diff_mod, ONLY: norm_lev_start, norm_lev_end
USE lam_config_inputs_mod, ONLY: n_rims_to_do

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE
!
! Description: This version interfaces to physics schemes in sequence:
!    energy correction              (optional)
!    microphysics (cloud and large scale precipitation schemes)
!    radiation
!    gravity wave drag
!
!          CALLed before Semi-Lagrangian in atmosphere timestep.
! Method:
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
!   Language: FORTRAN 90 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component: Control Atmos
!
! Declarations:
!
! Subroutine arguments ===============================================

! arguments with intent in. ie: input variables.

! Parallel setup variables
INTEGER ::                                                              &
  global_row_length                                                     &
                     ! number of points on a row
, global_rows                                                           &
                     ! NUMBER OF global rows
, n_proc                                                                &
             ! Total number of processors
, n_procx                                                               &
             ! Number of processors in longitude
, n_procy                                                               &
             ! Number of processors in latitude
, g_rows (0:n_proc-1)                                                   &
, g_row_length (0:n_proc-1)


! Model dimensions
INTEGER ::                                                              &
  row_length                                                            &
, rows                                                                  &
, n_rows                                                                &
, salt_dim1                                                             &
                !
, salt_dim2                                                             &
                ! Dimensions of sea-salt aerosol arrays
, salt_dim3                                                             &
                !
, land_points                                                           &
                ! IN No.of land points being processed, can be 0.
, bl_levels                                                             &
, dst_levels                                                            &
                ! number of deep soil temperature levels
, dsm_levels                                                            &
                ! number of deep soil moisture levels
, Ozone_levels                                                          &
, cloud_levels                                                          &
, land_ice_points                                                       &
                ! number of land ice points
, soil_points                                                           &
                ! number of soil points
, n_cca_levels                                                          &
                ! Number of levels for conv cloud
, tr_levels                                                             &
                ! number of tracer levels
, tr_ukca                                                               &
                ! number of UKCA tracers
                ! amount: 1 for 2D, nlevs for 3D.
, ntiles                                                                &
, co2_dim_len                                                           &
                !\ For dimension 3-D CO2 field to be passed
, co2_dim_row                                                           &
                !/ to rad_ctl
, co2_dim_lev   !/

LOGICAL ::                                                              &
  l_lbc_old                                                             &
                !  false for new lbc treatment
, l_ukca_chem                                                           &
          ! T for ukca chemistry
, l_ukca_set_trace_gases                                                &
          ! Switch for prescribing atmospheric gases in UKCA
, l_ukca_strat                                                          &
          ! True for stratospheric chemistry scheme
, l_ukca_strattrop                                                      &
          ! True for stratospheric+tropospheric chemistry scheme
, l_ukca_prescribech4
          ! Switch for prescribing surface ch4 in UKCA

REAL ::                                                                 &
              ! Intent(IN)
  m_ci        ! variable to modify ice fall speed for LSPCON3C
              !  for stochastic physics random parameters

INTEGER ::                                                              &
  min_trop_level                                                        &
                  ! Lowest permitted level for the tropopause
!                       ! used for radiative purposes.
      , max_trop_level  ! Highest permitted level for the tropopause
!                       ! used for radiative purposes.


REAL ::                                                                 &
  rhcrit(model_levels) ! IN Critical relative humidity.
                       ! the values need to be tuned
                       ! for the given set of levels.

! Diagnostics info
REAL ::                                                                 &
 stashwork1(*)                                                          &
                   ! STASH workspace for section 1 (SW rad)
,stashwork2(*)                                                          &
                   ! STASH workspace for section 2 (LW rad)
,stashwork4(*)                                                          &
                   ! STASH workspace for section 4 (LS precip)
,stashwork6(*)                                                          &
                   ! STASH workspace for section 6 (gw drag)
,stashwork14(*)                                                         &
                   ! STASH workspace for section 14 (Energy cor)
,stashwork21(*)
                   ! STASH workspace for section 21 (electric)

! Data arrays

! Small halo
REAL, INTENT (INOUT) :: u  (udims_s%i_start:udims_s%i_end,              &
                            udims_s%j_start:udims_s%j_end,              &
                            udims_s%k_start:udims_s%k_end)
REAL, INTENT (INOUT) :: v  (vdims_s%i_start:vdims_s%i_end,              &
                            vdims_s%j_start:vdims_s%j_end,              &
                            vdims_s%k_start:vdims_s%k_end)
REAL, INTENT (INOUT) :: rho(pdims_s%i_start:pdims_s%i_end,              &
                            pdims_s%j_start:pdims_s%j_end,              &
                            pdims_s%k_start:pdims_s%k_end)
REAL, INTENT (INOUT) :: p  (pdims_s%i_start:pdims_s%i_end,              &
                            pdims_s%j_start:pdims_s%j_end,              &
                            pdims_s%k_start:pdims_s%k_end)
REAL, INTENT (INOUT) :: p_theta_levels                                  &
                           (tdims_s%i_start:tdims_s%i_end,              &
                            tdims_s%j_start:tdims_s%j_end,              &
                            tdims_s%k_start:tdims_s%k_end)
REAL, INTENT (INOUT) :: theta                                           &
                           (tdims_s%i_start:tdims_s%i_end,              &
                            tdims_s%j_start:tdims_s%j_end,              &
                            tdims_s%k_start:tdims_s%k_end)
REAL, INTENT (INOUT) :: exner_rho_levels                                &
                           (pdims_s%i_start:pdims_s%i_end,              &
                            pdims_s%j_start:pdims_s%j_end,              &
                            pdims_s%k_start:pdims_s%k_end + 1)
REAL, INTENT (INOUT) :: exner_theta_levels                              &
                           (tdims_s%i_start:tdims_s%i_end,              &
                            tdims_s%j_start:tdims_s%j_end,              &
                            tdims_s%k_start:tdims_s%k_end)

!No halo
REAL, INTENT (INOUT) :: p_star (pdims%i_start:pdims%i_end,              &
                                pdims%j_start:pdims%j_end)

!Large halo
REAL, INTENT (INOUT) :: q  (tdims_l%i_start:tdims_l%i_end,              &
                            tdims_l%j_start:tdims_l%j_end,              &
                            tdims_l%k_start:tdims_l%k_end)
REAL, INTENT (INOUT) :: qcl(tdims_l%i_start:tdims_l%i_end,              &
                            tdims_l%j_start:tdims_l%j_end,              &
                            tdims_l%k_start:tdims_l%k_end)
REAL, INTENT (INOUT) :: qcf(tdims_l%i_start:tdims_l%i_end,              &
                            tdims_l%j_start:tdims_l%j_end,              &
                            tdims_l%k_start:tdims_l%k_end)
REAL, INTENT (INOUT) :: qcf2                                            &
                           (tdims_l%i_start:tdims_l%i_end,              &
                            tdims_l%j_start:tdims_l%j_end,              &
                            tdims_l%k_start:tdims_l%k_end)
REAL, INTENT (INOUT) :: qrain                                           &
                           (tdims_l%i_start:tdims_l%i_end,              &
                            tdims_l%j_start:tdims_l%j_end,              &
                            tdims_l%k_start:tdims_l%k_end)
REAL, INTENT (INOUT) :: qgraup                                          &
                           (tdims_l%i_start:tdims_l%i_end,              &
                            tdims_l%j_start:tdims_l%j_end,              &
                            tdims_l%k_start:tdims_l%k_end)
REAL, INTENT (INOUT) :: cf (tdims_l%i_start:tdims_l%i_end,              &
                            tdims_l%j_start:tdims_l%j_end,              &
                            tdims_l%k_start:tdims_l%k_end)
REAL, INTENT (INOUT) :: cfl(tdims_l%i_start:tdims_l%i_end,              &
                            tdims_l%j_start:tdims_l%j_end,              &
                            tdims_l%k_start:tdims_l%k_end)
REAL, INTENT (INOUT) :: cff(tdims_l%i_start:tdims_l%i_end,              &
                            tdims_l%j_start:tdims_l%j_end,              &
                            tdims_l%k_start:tdims_l%k_end)

REAL ::                                                                 &
  energy_correction

! ancillary arrays and fields required to be saved from timestep to
! timestep.

REAL ::                                                                 &
  T_surf(row_length, rows)
! The following are only used if coastal tiling is switched on:
REAL ::                                                                 &
 fland_ctile(land_points)                                               &
                             ! IN Land fraction on land points.
,tstar_land_ctile(row_length,rows)                                      &
!                                  ! IN Land mean sfc temperature (K)
      ,tstar_sea_ctile(row_length,rows)                                 &
!                                  ! IN Open sea sfc temperature (K).
      ,tstar_sice_ctile(row_length,rows,nice_use)                       &
!                                  ! IN Sea-ice sfc temperature (K).
      ,land_alb_ctile(row_length,rows)                                  &
!                                  ! INOUT Mean land albedo.
      ,sice_alb_ctile(row_length,rows)
!                                  ! INOUT Sea-ice albedo.

LOGICAL ::                                                              &
  land_sea_mask(row_length, rows)                                       &
, rad_mask(row_length, rows)
!  A mask which ensures a chequerboard pattern of radiation calculations
!  over the whole domain (not just one PE)

INTEGER ::                                                              &
  land_index (land_points)                                              &
                                ! set from land_sea_mask
, ntml(row_length, rows)

LOGICAL ::                                                              &
  cumulus (row_length, rows) ! bl convection flag

REAL ::                                                                 &
  snow_depth (row_length, rows)                                         &
                                ! snow/qrclim.snow.(month)
, snow_depth_sea_cat (row_length, rows, nice_use)                       &
                                ! snow depth on sea ice
, sd_orog_land (land_points)                                            &
                             ! orog/qrparm.orog.stdev
, orog_grad_xx_land(land_points)                                        &
                                 ! orog/qrparm.orog.sigmaxx
, orog_grad_xy_land(land_points)                                        &
                                 ! orog/qrparm.orog.sigmaxy
, orog_grad_yy_land(land_points)                                        &
                                 ! orog/qrparm.orog.sigmayy
, zh(row_length, rows)                                                  &
                       ! boundary layer height
, aerosol_em(tdims%i_start:tdims%i_end,                                 &
             tdims%j_start:tdims%j_end,                                 &
             tdims%k_start:tdims%k_end)

REAL ::                                                                 &
  snowdepth_p(land_points,ntiles)                                       &
                            ! Snow depth on ground on tiles (m)
, lake_h_ice_p(land_points)
                            ! FLake lake ice thickness (m)

REAL ::                                                                 &
  ozone(row_length, rows, ozone_levels)                                 &
, O3_trop_level(row_length,rows)                                        &
, O3_trop_height(row_length,rows)                                       &
, cos_zenith_angle(row_length, rows)                                    &
, T_trop_level(row_length,rows)                                         &
, T_trop_height(row_length,rows)                                        &
, SW_incs(row_length, rows, 0:model_levels+1)                           &
, LW_incs(row_length, rows, 0:model_levels)                             &
, dirpar_inc(row_length, rows)                                          &
, soot(row_length, rows)                                                &
                                 ! Snow soot
, ice_fract (row_length, rows)                                          &
                                 ! ice/qrclim.ice.(month)
, ice_fract_cat (row_length, rows, nice_use)                            &
                                 ! ice fraction per cat
, ice_thick_cat (row_length, rows, nice_use)
                                 ! effective ice thickness per cat

REAL ::                                                                 &
  albsoil(land_points)                                                  &
, albobs_sw(land_points)                                                &
, albobs_vis(land_points)                                               &
, albobs_nir(land_points)                                               &
, lai(land_points, npft)                                                &
, canht(land_points, npft)                                              &
, rgrain(land_points, ntiles)                                           &
, snow_tile(land_points, ntiles)                                        &
, tile_frac(land_points, ntype)                                         &
, tstar_tile(land_points, ntiles)                                       &
, z0_tile(land_points, ntiles)                                          &
, z0m_sea(row_length, rows)                                             &
, chloro_sea(row_length, rows)                                          &
, dOLR_rts(row_length, rows)                                            &
                                   ! TOA - surface upward LW
, LW_down(row_length, rows)                                             &
                                   ! Surface downward LW
, SW_tile_rts(land_points, ntiles)                                      &
                                   ! Surface net SW on land tiles
, es_space_interp(4, row_length, rows)
!              Coeffs for spatial interpolation of radiation quantities
!

! Aerosol climatology for NWP

      ! Number of requested species within the climatology
INTEGER :: n_arcl_species

! Corresponding number of requested components
INTEGER :: n_arcl_compnts

! Model switch for each species
LOGICAL :: l_use_arcl(npd_arcl_species)

! Mass-mixing ratios
REAL :: arcl (tdims%i_start:tdims%i_end,                                &
              tdims%j_start:tdims%j_end,                                &
                          1:tdims%k_end,                                &
                          1:n_arcl_compnts)

! Array index of each component
INTEGER :: i_arcl_compnts(npd_arcl_compnts)

! EasyAerosol distributions
TYPE (t_easyaerosol_rad), INTENT(IN) :: easyaerosol_sw
TYPE (t_easyaerosol_rad), INTENT(IN) :: easyaerosol_lw
TYPE (t_easyaerosol_cdnc), INTENT(IN) :: easyaerosol_cdnc

! Convection Scheme

REAL    :: ccw(row_length, rows, model_levels)
INTEGER :: lcbase(row_length, rows)

REAL ::                                                                 &
  cca (row_length, rows, n_cca_levels)                                  &
, cca_dp (tdims%i_start:tdims%i_end,                                    &
          tdims%j_start:tdims%j_end, n_cca_levels)                      &
, cca_md (tdims%i_start:tdims%i_end,                                    &
          tdims%j_start:tdims%j_end, n_cca_levels)                      &
, cca_sh (tdims%i_start:tdims%i_end,                                    &
          tdims%j_start:tdims%j_end, n_cca_levels)                      &
, cclwp(row_length, rows) ! condensed water path (KG/M**2)

INTEGER ::                                                              &
  ccb (row_length, rows)                                                &
, cct (row_length, rows)

! Co-ordinate arrays

REAL ::                                                                 &
  delta_lambda                                                          &
, delta_phi

! time information for current timestep
INTEGER ::                                                              &
  val_year                                                              &
, val_day_number                                                        &
, val_hour                                                              &
, val_minute                                                            &
, val_second                                                            &
,         previous_time(7)                                              &
,     asteps_since_triffid   ! INOUT  Number of atmospheric
!                                  !        timesteps since last call
!                                  !        to TRIFFID.
! UKCA_RADAER structure
TYPE (ukca_radaer_struct) :: ukca_radaer

! Diagnostic variables
REAL ::                                                                 &
  lat_rot_NP                                                            &
, long_rot_NP

REAL ::                                                                 &
  area_cloud_fraction(row_length, rows, model_levels)

! Variables for COSP
! 3D convective rainfall flux COSP
REAL,INTENT(IN) :: cosp_crain_3d(row_length,rows,model_levels)
! 3D convective snowfall flux COSP
REAL,INTENT(IN) :: cosp_csnow_3d(row_length,rows,model_levels)

! arguments with intent in/out. ie: input variables changed on output.

REAL, INTENT (INOUT) ::                                                 &
  sum_eng_fluxes(row_length,rows)                                       &
                                   ! sum atmosphere fluxes
, sum_moist_flux(row_length,rows)  ! sum moist fluxes

REAL, INTENT (INOUT) :: theta_star                                      &
                           (tdims%i_start:tdims%i_end,                  &
                            tdims%j_start:tdims%j_end,                  &
                            tdims%k_start:tdims%k_end)
REAL, INTENT (INOUT) :: q_star                                          &
                           (tdims%i_start:tdims%i_end,                  &
                            tdims%j_start:tdims%j_end,                  &
                            tdims%k_start:tdims%k_end)
REAL, INTENT (INOUT) :: qcl_star                                        &
                           (tdims%i_start:tdims%i_end,                  &
                            tdims%j_start:tdims%j_end,                  &
                            tdims%k_start:tdims%k_end)
REAL, INTENT (INOUT) :: qcf_star                                        &
                           (tdims%i_start:tdims%i_end,                  &
                            tdims%j_start:tdims%j_end,                  &
                            tdims%k_start:tdims%k_end)
REAL, INTENT (INOUT) :: qcf2_star                                       &
                           (tdims%i_start:tdims%i_end,                  &
                            tdims%j_start:tdims%j_end,                  &
                            tdims%k_start:tdims%k_end)
REAL, INTENT (INOUT) :: qrain_star                                      &
                           (tdims%i_start:tdims%i_end,                  &
                            tdims%j_start:tdims%j_end,                  &
                            tdims%k_start:tdims%k_end)
REAL, INTENT (INOUT) :: qgraup_star                                     &
                           (tdims%i_start:tdims%i_end,                  &
                            tdims%j_start:tdims%j_end,                  &
                            tdims%k_start:tdims%k_end)
REAL, INTENT (INOUT) :: cf_star                                         &
                           (tdims%i_start:tdims%i_end,                  &
                            tdims%j_start:tdims%j_end,                  &
                            tdims%k_start:tdims%k_end)
REAL, INTENT (INOUT) :: cfl_star                                        &
                           (tdims%i_start:tdims%i_end,                  &
                            tdims%j_start:tdims%j_end,                  &
                            tdims%k_start:tdims%k_end)
REAL, INTENT (INOUT) :: cff_star                                        &
                           (tdims%i_start:tdims%i_end,                  &
                            tdims%j_start:tdims%j_end,                  &
                            tdims%k_start:tdims%k_end)
REAL, INTENT (INOUT) :: u_inc                                           &
                           (udims_s%i_start:udims_s%i_end,              &
                            udims_s%j_start:udims_s%j_end,              &
                            udims_s%k_start:udims_s%k_end)
REAL, INTENT (INOUT) :: v_inc                                           &
                           (vdims_s%i_start:vdims_s%i_end,              &
                            vdims_s%j_start:vdims_s%j_end,              &
                            vdims_s%k_start:vdims_s%k_end)

REAL, INTENT(INOUT) ::  aerosol  (tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(INOUT) ::  flash_pot(tdims%i_start:tdims%i_end,            &
                                  tdims%j_start:tdims%j_end,            &
                                  tdims%k_start:tdims%k_end)

REAL, INTENT(INOUT) ::  dust_div1(tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(INOUT) ::  dust_div2(tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(INOUT) ::  dust_div3(tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(INOUT) ::  dust_div4(tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(INOUT) ::  dust_div5(tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(INOUT) ::  dust_div6(tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(INOUT) ::  so2      (tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  so4_aitken                                      &
                                 (tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  so4_accu (tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  so4_diss (tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  nh3      (tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  soot_new (tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  soot_agd (tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  soot_cld (tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  bmass_new(tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  bmass_agd(tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  bmass_cld(tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  ocff_new (tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  ocff_agd (tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  ocff_cld (tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  nitr_acc (tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  nitr_diss(tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  co2      (tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::  tracer   (tdims_s%i_start:tdims_s%i_end,  &
                                  tdims_s%j_start:tdims_s%j_end,  &
                                  tdims_s%k_start:tdims_s%k_end,  &
                                  tr_vars)

REAL, INTENT(INOUT) ::  tracer_ukca                               &
                                 (tdims_s%i_start:tdims_s%i_end,  &
                                  tdims_s%j_start:tdims_s%j_end,  &
                                  tdims_s%k_start:tdims_s%k_end,  &
                                  tr_ukca)

REAL, INTENT(IN) :: biogenic (tdims%i_start:tdims%i_end,                &
                              tdims%j_start:tdims%j_end,                &
                                          1:tdims%k_end)

REAL, INTENT(IN) :: bl_w_var(tdims%i_start : tdims%i_end,               &
                             tdims%j_start : tdims%j_end,               &
                                         1 : tdims%k_end)

REAL, INTENT(IN) :: totalppn (tdims%i_start:tdims%i_end,                &
                              tdims%j_start:tdims%j_end)

! arguments with intent out. ie: output variables.

! arrays holding information to be passed between physics
! routines.

REAL, INTENT(OUT) ::                                                    &
  ls_rain(row_length, rows)                                             &
, ls_snow(row_length, rows)                                             &
, ls_graup(row_length, rows)                                            &
, micro_tends(row_length, rows, 2, bl_levels)                           &
!                          ! Tendencies from microphys within BL levels
!                          ! (TL, K/s; QW, kg/kg/s)
      , ls_rainfrac(land_points)
                        ! Rain fraction array on land points to be
                        ! passed to atmos_physics2

! Radiation fields 1. SW & common with LW.
REAL ::                                                                 &
  photosynth_act_rad(row_length, rows)                                  &
                                       ! Net downward
!                                 shortwave radiation in band 1 (w/m2).
      , rad_hr(row_length, rows, 2, bl_levels)                          &
                                               !
!                               ! BL radiative (LW,SW) heating rates
      , dOLR(row_length, rows)                                          &
                                    ! TOA - surface upward LW
      , SW_tile(land_points, ntiles)! Surface net SW on land tiles

! Position of greenhouse gases in tracer_ukca array
INTEGER, INTENT(IN) :: ngrgas
INTEGER, INTENT(IN) :: grgas_addr(ngrgas)

! UKCA cloud drop number concentration dimensions
INTEGER, INTENT(IN) :: cdnc_dim1
INTEGER, INTENT(IN) :: cdnc_dim2
INTEGER, INTENT(IN) :: cdnc_dim3

!  Additional variables for SCM diagnostics (dummy in full UM)
INTEGER ::                                                              &
  nSCMDpkgs             ! No of SCM diagnostics packages

LOGICAL ::                                                              &
  L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages

CHARACTER(LEN=errormessagelength) ::                                    &
   cmessage
INTEGER ::                                                              &
  Error_code

REAL, INTENT(IN):: w(wdims_s%i_start:wdims_s%i_end,                     &
                     wdims_s%j_start:wdims_s%j_end,                     &
                     wdims_s%k_start:wdims_s%k_end)


! === End of arguments ==============================================

! Local parameters:
CHARACTER(LEN=*) :: RoutineName
PARAMETER (   RoutineName='ATMOS_PHYSICS1')
REAL ::                                                                 &
  amp                                                                   &
                  ! amplitude of diurnal variation in tracers
, tau_decay                                                             &
                  ! time constant for decay of tracer to clim
, clim_murk_land                                                        &
                  ! climatological murk value over land points
, clim_murk_sea   ! climatological murk value over sea points
PARAMETER (                                                             &
            amp=0.7                                                     &
,           tau_decay=1.728e5                                           &
,           clim_murk_land=25.0                                         &
,           clim_murk_sea=12.5                                          &
           )

! Local scalars:
REAL :: z0_sea    ! Roughness length over sea (m)
REAL :: p1
REAL :: windspeed_1
REAL :: windspeed_10m(salt_dim1, salt_dim2)

! loop counters
INTEGER ::                                                              &
  i, j, k                                                               &
, l, n

! local variables
INTEGER ::                                                              &
  rhc_row_length                                                        &
                  ! Row length for RHcrit array
, rhc_rows                                                              &
                  ! Row number for RHcrit array
, lspice_dim1,lspice_dim2,lspice_dim3                                   &
!                       ! Required for array dimensions MCR_CTL2
      , sulp_dim1                                                       &
                                ! dimensions for sulphate arrays in
      , dust_dim1                                                       &
                                ! dimensions for mineral dust arrays i
      , dust_dim2                                                       &
                                ! in rad_ctl
      , biogenic_dim1                                                   &
                                ! dimensions of biogenic arrays in
      , biogenic_dim2                                                   &
                                !
      , sulp_dim2                                                       &
                                !   RAD_CTL
      , soot_dim1, soot_dim2                                            &
                                ! dimensions of soot arrays in RAD_CTL
      , bmass_dim1, bmass_dim2                                          &
                                ! dimensions of biomass arrays in radn
      , ocff_dim1, ocff_dim2                                            &
                                ! dimensions of OCFF arrays in radiation
      , nitrate_dim1, nitrate_dim2                                      &
                                ! dimensions of nitrate arrays in radn
      , arcl_dim1, arcl_dim2                                            &
                     ! dimensions of aerosol clim for NWP arrays in radn
      , ukca_dim1, ukca_dim2                                            &
                     ! dimensions of UKCA_RADAER arrays in RAD_CTL
      , i_start                                                         &
                        ! Row start point for polar row tidying
      , istat           ! Status (error code) indicator

! Local data arrays

REAL ::                                                                 &
  u_1(salt_dim1, salt_dim2)                                             &
, v_1(salt_dim1, salt_dim2)                                             &
, co2_3D(co2_dim_len, co2_dim_row, co2_dim_lev)                         &
, u_1_mean                                                              &
, v_1_mean                                                              &
, frac_control(land_points,ntype)

! local arrays and scalar for a 10m windspeed on the p grid on sea points:
! Needed for radiation, sea salt and DMS:
REAL ::                                                                 &
  height_theta(pdims%i_len, pdims%j_len, pdims%k_len),                  &
!          theta level centre height above surface
        height_rho_1,                                                          &
!          theta level 1 rho centre height above surface
        u_l1p(pdims%i_len, pdims%j_len),                                &
!          u wind on level 1, on the p grid
        v_l1p(pdims%i_len, pdims%j_len),                                &
!          v wind on level 1, on the p grid
        ws_l1p,                                                                &
!          wind speed on level 1, on the p grid
        p_l1_sea,                                                              &
!          10m, in logspace with roughness, on the p grid, for sea points
        ws_10m_sea(pdims%i_len, pdims%j_len)
!          wind speed at 10m, on the p grid, calculated for sea points

REAL :: u_1_arr(1), v_1_arr(1)

REAL :: T_n   (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
                           1:tdims%k_end)
REAL :: T_n2  (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
                           1:tdims%k_end)
REAL :: q_n   (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
                           1:tdims%k_end)
REAL :: qcl_n (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
                           1:tdims%k_end)
REAL :: qcf_n (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
                           1:tdims%k_end)
REAL :: cf_n  (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
                           1:tdims%k_end)
REAL :: cfl_n (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
                           1:tdims%k_end)
REAL :: cff_n (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
                           1:tdims%k_end)
REAL :: u_on_p(pdims%i_start:pdims%i_end,                               &
               pdims%j_start:pdims%j_end,                               &
               pdims%k_start:pdims%k_end)
REAL :: v_on_p(pdims%i_start:pdims%i_end,                               &
               pdims%j_start:pdims%j_end,                               &
               pdims%k_start:pdims%k_end)
REAL :: theta_inc (tdims%i_start:tdims%i_end,                           &
                   tdims%j_start:tdims%j_end,                           &
                               1:tdims%k_end)
REAL :: q_inc     (tdims%i_start:tdims%i_end,                           &
                   tdims%j_start:tdims%j_end,                           &
                               1:tdims%k_end)
REAL :: qcl_inc   (tdims%i_start:tdims%i_end,                           &
                   tdims%j_start:tdims%j_end,                           &
                               1:tdims%k_end)
REAL :: qcf_inc   (tdims%i_start:tdims%i_end,                           &
                   tdims%j_start:tdims%j_end,                           &
                               1:tdims%k_end)
REAL :: cf_inc    (tdims%i_start:tdims%i_end,                           &
                   tdims%j_start:tdims%j_end,                           &
                               1:tdims%k_end)
REAL :: cfl_inc   (tdims%i_start:tdims%i_end,                           &
                   tdims%j_start:tdims%j_end,                           &
                               1:tdims%k_end)
REAL :: cff_inc   (tdims%i_start:tdims%i_end,                           &
                   tdims%j_start:tdims%j_end,                           &
                               1:tdims%k_end)

! Dummy array to pass to pc2_turbulence_ctl.
REAL :: dummy_zeros(tdims%i_start:tdims%i_end,                          &
                    tdims%j_start:tdims%j_end,                          &
                    tdims%k_start:tdims%k_end)


! Potential droplet number (includes values where cloud not present)
REAL :: n_drop_pot(tdims%i_start : tdims%i_end,                         &
                   tdims%j_start : tdims%j_end,                         &
                               1 : tdims%k_end )

! Scaled aerosol emissions for murk
REAL :: aerosol_em_scaled(tdims%i_start:tdims%i_end,                    &
                          tdims%j_start:tdims%j_end)

! COSP variables
TYPE(cosp_config) :: cosp_cfg
TYPE(cosp_subgrid) :: cosp_sgx  ! Subgrid inputs
TYPE(cosp_sghydro) :: cosp_sgh  ! Subgrid hydrometeors
REAL :: cosp_frac_agg(tdims%i_start : tdims%i_end,                      &
                      tdims%j_start : tdims%j_end,                      &
                      1 : tdims%k_end )
INTEGER :: cosp_npoints
LOGICAL :: L_cosp_call, L_cosp_run

! Local logical which is set by tests dependent on model_type
LOGICAL :: l_RequestedDiags

! Local Arrays to store additional microphysics fields if in use
REAL, ALLOCATABLE ::                                 &
  qcf2_n(:,:,:),   qrain_n(:,:,:),   qgraup_n (:,:,:)                   &
, qcf2_inc(:,:,:), qrain_inc(:,:,:), qgraup_inc (:,:,:)                 &
, T_inc_diag(:,:,:),  q_inc_diag(:,:,:), qcl_inc_diag(:,:,:)            &
, qcf_inc_diag (:,:,:), cf_inc_diag(:,:,:), cfl_inc_diag(:,:,:)         &
, cff_inc_diag(:,:,:)

REAL ::                                                                 &
  p_layer_boundaries(row_length, rows, 0:model_levels)                  &
        ! pressure at layer boundaries. Same as p except at
        ! bottom level = pstar, and at top = 0.
, p_layer_centres(row_length, rows, 0:model_levels)               
        ! pressure at layer centres. Same as p_theta_levels
        !except bottom level = pstar, and at top = 0.


REAL :: moist                                                     &
                           (tdims_l%i_start:tdims_l%i_end,        &
                            tdims_l%j_start:tdims_l%j_end,        &
                            tdims_l%k_start:tdims_l%k_end)

! holds total moisture for conversion from wet to dry density
REAL :: unscaled_dry_rho(pdims_s%i_start:pdims_s%i_end,                 &
                         pdims_s%j_start:pdims_s%j_end,                 &
                         pdims_s%k_start:pdims_s%k_end)

! unscaled dry density
REAL :: weight1
REAL :: weight2
REAL :: weight3
REAL :: temp

REAL, TARGET, ALLOCATABLE ::                                            &
  ext_p_layer_centres(:,:,:),                                           &
  ext_tl(:,:,:),                                                        &
  ext_ql(:,:,:),                                                        &
  ext_qcf(:,:,:),                                                       &
  ext_ice_frac(:,:),                                                    &
  ext_land_frac(:,:)

INTEGER :: i_field   ! field counter for multivar swapbounds
TYPE(swapable_field_pointer_type) :: fields_to_swap(4)
                     ! mv swapbounds

REAL ::                                                                 &
 fland(land_points)                                                     &
                             ! Land fraction on land points.
,tstar_sea(row_length,rows)                                             &
                             ! Open sea sfc temperature (K).
,tstar_sice_cat(row_length,rows,nice_use)                               &
                             ! Sea-ice sfc temperature (K).
,land_alb(row_length,rows)                                              &
                             ! Mean land albedo.
,sice_alb(row_length,rows)   ! Sea-ice albedo.

LOGICAL ::                                                              &
 land0p5(row_length, rows)

! Diagnostics controlled by Diagnostic switches

! fields output by ls_cld not using stash flag

! Energy correction work variables

REAL ::                                                           &
  tot_precip_scaled_1(row_length, rows)                           &
, tot_precip_scaled_2(row_length, rows)                           &
, lclf

! Checking moisture increments
REAL, ALLOCATABLE ::                                                    &
  ep1(:,:)             & ! method 1
 ,ep2(:,:)             & ! method 2
 ,ep3(:,:)             & ! method 3
 ,ep4(:,:)             & ! copy of stash diagnostics i.e. precip
 ,q_methane(:,:,:)     & ! Methane increment
 ,dqt(:,:,:)           & ! total increment to moist in units (kg/kg/s)
 ,rho_con(:,:,:)         ! Density *r*r (kg/m) either wet density (specfic)
                         ! or dry density (mixing ratio)

! Temporary logicals
LOGICAL ::                                                              &
  L_zero_boundaries                                                     &
, L_zero_halos                                                          &
, L_use_dirpar

INTEGER ::                                                              &
  level

REAL, ALLOCATABLE :: grgas_field(:,:,:,:)

!     Global topmost cloudy level
INTEGER :: global_cloud_top = 0

!     Level of tropopause
INTEGER :: trindx(tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end)

!     SCM Dummy variables to keep call to tropin consistent.
REAL :: scm_dummy_1d(1,1), scm_dummy_2d(1,1,0:model_levels)

REAL ::                                                                 &
  sea_salt_film(salt_dim1, salt_dim2, salt_dim3),                       &
!          Film-mode sea-salt aerosol number concentration
        sea_salt_jet(salt_dim1, salt_dim2, salt_dim3),                  &
!          Jet-mode sea-salt aerosol number concentration
        height(salt_dim1, salt_dim2, salt_dim3)
!          Layer-centre height above surface

REAL :: ukca_cdnc(cdnc_dim1, cdnc_dim2, cdnc_dim3)
!          Cloud droplet number from UKCA-MODE

REAL :: land_fract(tdims%i_start:tdims%i_end,                           &
                   tdims%j_start:tdims%j_end)


LOGICAL :: l_pc2_prod_qcl_mp = .FALSE. ! Apply dqcl from mp calc

CHARACTER(LEN=60) :: scheme_name  ! description for moisture checking

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


! needed for endgame for fv_cos_theta_latitude vs cos_theta_latitude
REAL, POINTER :: xx_cos_theta_latitude (:,:)

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

xx_cos_theta_latitude => cos_theta_latitude

! ----------------------------------------------------------------------
! Section 0. Norms printing before physics1
! ----------------------------------------------------------------------
IF ( L_print_L2norms ) THEN
  WRITE(umMessage,'(A)') ' ***  L2 norms before atmos_physics1 ***'
  CALL umPrint(umMessage,src='ATMOS_PHYSICS1')
  CALL atmos_start_norm(                                                &
                        norm_lev_start, norm_lev_end, n_rims_to_do,     &
                        exner_rho_levels, rho, u, v, w,                 &
                        theta, q, qcl, qcf, qcf2, qrain, qgraup,        &
                        cf, cfl, cff, .TRUE., .FALSE., .FALSE. )

  IF ( L_tracer ) THEN
    WRITE(umMessage, '(A)')                                             &
                   ' *** L2 norms of tracers before atmos_physics1 ***'
    CALL umPrint(umMessage,src='ATMOS_PHYSICS1')
    CALL sl_tracer1_norm(                                               &
                         super_array_size,                              &
                         norm_lev_start, norm_lev_end, n_rims_to_do,    &
                         co2, aerosol, dust_div1, dust_div2,            &
                         dust_div3, dust_div4, dust_div5, dust_div6,    &
                         soot_new, soot_agd, soot_cld,                  &
                         bmass_new, bmass_agd, bmass_cld,               &
                         ocff_new, ocff_agd, ocff_cld,                  &
                         so2, so4_aitken, so4_accu, so4_diss,           &
                         nh3, nitr_acc, nitr_diss,                      &
                         tracer_ukca, tr_ukca, ozone,                   &
                         .TRUE., .FALSE., .FALSE. )
  END IF  !  L_tracer
END IF !  L_print_L2norms

! ----------------------------------------------------------------------
! Section INI. Initialisation of variables.
! ----------------------------------------------------------------------

z0_sea = 2.5e-04    ! Roughness length over sea (m)
p1     = LOG(10.0/z0_sea)


IF (L_mcr_qcf2) THEN  ! Second cloud ice variable in use
  ALLOCATE ( qcf2_inc(row_length, rows, model_levels) )
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(i,j,k)                          &
!$OMP SHARED(model_levels,rows,row_length,qcf2_inc)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        qcf2_inc(i,j,k) = 0.0
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  ALLOCATE ( qcf2_inc(1,1,1) )
END IF

IF (L_mcr_qrain) THEN  ! Prognostic rain in use
  ALLOCATE ( qrain_inc(row_length, rows, model_levels) )
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(i,j,k)                          &
!$OMP SHARED(model_levels,rows,row_length,qrain_inc)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        qrain_inc(i,j,k) = 0.0
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  ALLOCATE ( qrain_inc(1,1,1) )
END IF

IF (L_mcr_qgraup) THEN  ! Prognostic graupel in use
  ALLOCATE ( qgraup_inc(row_length, rows, model_levels) )
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(i,j,k)                          &
!$OMP SHARED(model_levels,rows,row_length,qgraup_inc)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        qgraup_inc(i,j,k) = 0.0
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  ALLOCATE ( qgraup_inc(1,1,1) )
END IF

IF ( l_casim ) CALL casim_alloc_increments

! Allocate arrays required for moisture conservation checking
IF (l_check_moist_inc .AND. model_type == mt_global .AND.               &
     (l_rain .OR. l_use_methox)) THEN

  ! ep4 holds diagnostic eg total large-scale precipitation or
  ! 0 for methane
  ALLOCATE (ep1(row_length, rows) )
  ALLOCATE (ep2(row_length, rows) )
  ALLOCATE (ep3(row_length, rows) )
  ALLOCATE (ep4(row_length, rows) )
  ALLOCATE (rho_con(pdims_s%i_start:pdims_s%i_end,                      &
                      pdims_s%j_start:pdims_s%j_end,                    &
                      pdims_s%k_start:pdims_s%k_end) )
  ALLOCATE (dqt(row_length, rows, model_levels) )
END IF

IF (l_triffid) THEN
  ! Increment counter for number of atmosphere timesteps since last
  ! call to TRIFFID vegetation model
  asteps_since_triffid = asteps_since_triffid + 1
END IF

! MORUSES prognostics mapped in surf_couple_initialise as some are required
! by JULES before now.

! map JULES prognostics to module prognostics before use
snowdepth_surft = snowdepth_p
IF (l_flake_model) THEN
  lake_h_ice_gb = lake_h_ice_p
END IF

! Allocate arrays of slow physics tendencies if SPT is TRUE
IF (l_retain_slow_tendencies) CALL init_slowphys_tendencies()

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, l, weight2, weight1,   &
!$OMP weight3, temp)

! set p at layer boundaries.
! NB: some arrays have haloes but are unset, if never used they will
!     be removed.

!$OMP DO SCHEDULE(STATIC)
DO j = 1, rows
  DO i = 1, row_length
    p_layer_boundaries(i,j,0) = p_star(i,j)
    p_layer_centres(i,j,0) = p_star(i,j)
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels - 1
  DO j = 1, rows
    DO i = 1, row_length
      p_layer_boundaries(i,j,k) = p(i,j,k+1)
      p_layer_centres(i,j,k) = p_theta_levels(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

k = model_levels

!$OMP DO SCHEDULE(STATIC)
DO j = 1, rows
  DO i = 1, row_length
    p_layer_boundaries(i,j,model_levels) = 0.0
    p_layer_centres(i,j,model_levels) =                                 &
                           p_theta_levels(i,j,model_levels)
  END DO
END DO
!$OMP END DO NOWAIT

! Height of theta levels above the surface (without orography)
!$OMP DO SCHEDULE(STATIC)
DO k = pdims%k_start, pdims%k_end
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      height_theta(i,j,k) = r_theta_levels(i,j,k) - r_theta_levels(i,j,0)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

! For idealised forcing, use the increments generated in eg_idl_forcing:
! - convert the potential temperature increment into a temperature increment;
! - wind increments are in u_inc & v_inc (initialised here if not idealised).
IF (l_idealised_data) THEN

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        theta_inc(i,j,k) = theta_star(i,j,k) * exner_theta_levels(i,j,k)
        q_inc(i,j,k)     = q_star(i,j,k)
        qcl_inc(i,j,k)   = 0.0
        qcf_inc(i,j,k)   = 0.0
        cf_inc(i,j,k)    = 0.0
        cfl_inc(i,j,k)   = 0.0
        cff_inc(i,j,k)   = 0.0
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

ELSE

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        theta_inc(i,j,k) = 0.0
        q_inc(i,j,k)     = 0.0
        qcl_inc(i,j,k)   = 0.0
        qcf_inc(i,j,k)   = 0.0
        cf_inc(i,j,k)    = 0.0
        cfl_inc(i,j,k)   = 0.0
        cff_inc(i,j,k)   = 0.0
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO k = udims_s%k_start, udims_s%k_end
    DO j = udims_s%j_start, udims_s%j_end
      DO i = udims_s%i_start, udims_s%i_end
        u_inc(i,j,k) = 0.0
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO k = vdims_s%k_start, vdims_s%k_end
    DO j = vdims_s%j_start, vdims_s%j_end
      DO i = vdims_s%i_start, vdims_s%i_end
        v_inc(i,j,k) = 0.0
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

END IF

!$OMP DO SCHEDULE(STATIC)
DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end
    land_fract(i,j) = 0.0
  END DO
END DO
!$OMP END DO

! ----------------------------------------------------------------------
! Set Coastal tiling dependent prognostics:
! ----------------------------------------------------------------------
IF (l_ctile) THEN

!$OMP DO SCHEDULE(STATIC)
  DO j = 1, rows
    DO i = 1, row_length
      tstar_sea(i,j)=tstar_sea_ctile(i,j)
      land_alb(i,j)=land_alb_ctile(i,j)
      sice_alb(i,j)=sice_alb_ctile(i,j)
      land0p5(i,j)=.FALSE.
    END DO
  END DO
!$OMP END DO
!wait land0p5

  DO k = 1, nice_use
!$OMP DO SCHEDULE(STATIC)
    DO j = 1, rows
      DO i = 1, row_length
        tstar_sice_cat(i,j,k)=tstar_sice_ctile(i,j,k)
      END DO
    END DO
!$OMP END DO NOWAIT
  END DO

!$OMP DO SCHEDULE(STATIC)
  DO l=1,land_points
    j=(land_index(l)-1)/row_length + 1
    i = land_index(l) - (j-1)*row_length
    fland(l)=fland_ctile(l)
    IF (fland(l) >= 0.5)land0p5(i,j)=.TRUE.
    land_fract(i,j) = fland(l)
  END DO
!$OMP END DO NOWAIT

ELSE
  ! No coastal tiling so land fraction is just 1 for all land points.
!$OMP DO SCHEDULE(STATIC)
  DO l=1,land_points
    j=(land_index(l)-1)/row_length + 1
    i = land_index(l) - (j-1)*row_length
    fland(l) = 1.0
    land_fract(i,j) = 1.0
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  !(Note, nice_use must = 1 if not using coastal tiling)
  DO j = 1, rows
    DO i = 1, row_length
      IF (.NOT. land_sea_mask(i,j)) THEN
        IF (ice_fract(i,j) <= 0.0) THEN
          tstar_sea(i,j)=t_surf(i,j)
          tstar_sice_cat(i,j,1)=t_surf(i,j)
        ELSE
          tstar_sea(i,j)=tfs
          tstar_sice_cat(i,j,1)=(t_surf(i,j)                            &
            -(1.0-ice_fract(i,j))*tstar_sea(i,j))                       &
                /ice_fract(i,j)
        END IF
      ELSE
        tstar_sea(i,j)=t_surf(i,j)
        tstar_sice_cat(i,j,1)=t_surf(i,j)
      END IF
      land0p5(i,j)=land_sea_mask(i,j)

      land_alb(i,j)=rmdi
      sice_alb(i,j)=rmdi
    END DO
  END DO
!$OMP END DO  NOWAIT
END IF

!$OMP END PARALLEL

!-----------------------------------------------------------------------
! Initialise sea and sea-ice variables for JULES that require a land
! mask to be available  (also used by MOSES in the radiation scheme)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Set up index for sea and sea-ice
!-----------------------------------------------------------------------
ssi_pts = 0
ssi_index(:)=0
DO j=1,rows
  DO i=1,row_length
    ssi_pts=ssi_pts + 1
    IF ( land_fract(i,j) < 1.0 ) THEN
      ssi_index(ssi_pts) = (j-1)*row_length + i
    END IF
    fssi_ij(i,j)=1.0 - land_fract(i,j)
  END DO
END DO

!-----------------------------------------------------------------------
! Allocate space for and initialise sea and sea-ice indices.
!-----------------------------------------------------------------------

sea_pts = 0
sice_pts = 0
sea_index(:)=0
sice_index(:)=0
sice_frac(:)=0.0
sea_frac(:)=0.0
DO l=1,ssi_pts
  j=(ssi_index(l)-1)/row_length + 1
  i = ssi_index(l) - (j-1)*row_length
  IF (ssi_index(l) > 0) THEN
    IF (ice_fract(i,j) > 0.0) THEN
      sice_pts=sice_pts+1
      sice_index(sice_pts)=l
      sice_frac(l)=ice_fract(i,j)
    END IF
    IF (ice_fract(i,j) < 1.0) THEN
      sea_pts=sea_pts+1
      sea_index(sea_pts)=l
      sea_frac(l)=1.0 - sice_frac(l)
    END IF
  END IF
END DO

! NOTE these settings use NICE_USE not NICE so are suitable for use here and
! in the explicit part of the surface exchange code.  They are then updated
! using NICE before the call to the implicit part of the surface exchange
! code
sice_pts_ncat(:)=0
sice_index_ncat(:,:)=0
sice_frac_ncat(:,:)=0.0
DO n=1,nice_use
  DO l=1,ssi_pts
    j=(ssi_index(l)-1)/row_length + 1
    i = ssi_index(l) - (j-1)*row_length
    IF (ssi_index(l) > 0) THEN
      IF (ice_fract_cat(i,j,n) > 0.0) THEN
        sice_pts_ncat(n)=sice_pts_ncat(n)+1
        sice_index_ncat(sice_pts_ncat(n),n)=l
        sice_frac_ncat(l,n)=ice_fract_cat(i,j,n)
      END IF
    END IF
  END DO
END DO

! Interpolate the u and v winds to produce a 10m wind speed (for the radiation
! scheme, sea salt and DMS emission):
!
! calculate level 1 winds on p points using subroutine
! that deals with grids and halos
CALL uv_p_pnts_halo(u(:,:,udims%k_start), v(:,:,vdims%k_start),         &
                      cos_theta_longitude, sin_theta_longitude,         &
                      global_row_length, gc_proc_row_group,             &
                      u_l1p, v_l1p)

! Calculate level 1 windspeed (NB on rho levels)
! Adjust level 1 windspeed to be 10m windspeed

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, ws_l1p, p_l1_sea,         &
!$OMP height_rho_1)                                                     &
!$OMP SHARED(pdims,u_l1p,v_l1p,z0m_sea,ws_10m_sea,r_rho_levels,         &
!$OMP     r_theta_levels,halo_j,halo_i)                                 &
!$OMP SCHEDULE(STATIC)   
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end

    height_rho_1 = r_rho_levels(i-halo_i, j-halo_j, 1) -                 &
                        r_theta_levels(i-halo_i, j-halo_j, 0)

    ws_l1p = SQRT(u_l1p(i,j)**2 + v_l1p(i,j)**2)

    p_l1_sea = LOG(10.0 / z0m_sea(i,j))

    ws_10m_sea(i,j) = ws_l1p * p_l1_sea /                                &
                      (LOG(height_rho_1 / z0m_sea(i,j)) )
  END DO
END DO
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! Section ENG.1  Add energy correction increments to temperature
! ----------------------------------------------------------------------
IF (Ltimer) CALL timer ('AP1 Energy Correct.',5)

IF (L_emcorr) THEN

  IF (Lflux_reset) THEN
    ! reinitialise net flux field at beginning of energy correction period
    DO j = 1, rows
      DO i = 1, row_length
        sum_eng_fluxes(i,j)=0.0
        sum_moist_flux(i,j)=0.0
      END DO
    END DO
  END IF

  ! Add energy correction increments every timestep.
  ! This is a temperature increment

  CALL add_eng_corr (energy_correction,Theta_inc,                       &
                         STASHwork14)

END IF    ! (L_emcorr)

IF (Ltimer) CALL timer ('AP1 Energy Correct.',6)

! ----------------------------------------------------------------------
! Section Set-up time-level n
! ----------------------------------------------------------------------

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, j, k )                          &
!$OMP SHARED(tdims,T_n,theta,exner_theta_levels,q_n,qcl_n,qcf_n,  &
!$OMP    cf_n,cfl_n,cff_n,q,qcl,qcf,cf,cfl,cff,T_n2)

!$OMP DO SCHEDULE(STATIC)   
DO k =             1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      T_n(i,j,k) = theta(i,j,k) * exner_theta_levels(i,j,k)
      T_n2(i,j,k) = theta(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)   
DO k =             1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      q_n(i,j,k) = q(i,j,k)
      qcl_n(i,j,k) = qcl(i,j,k)
      qcf_n(i,j,k) = qcf(i,j,k)
      cf_n(i,j,k)  = cf(i,j,k)
      cfl_n(i,j,k) = cfl(i,j,k)
      cff_n(i,j,k) = cff(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

IF (falliceshear_method == real_shear) THEN

!$OMP PARALLEL DEFAULT(NONE)                                     &
!$OMP SHARED(u,udims_s,pdims,model_levels,at_extremity,u_on_p,   &
!$OMP    vdims_s,v,v_on_p)

  CALL u_to_p (u, udims_s%i_start,udims_s%i_end,                  &
               udims_s%j_start,udims_s%j_end,                     &
               pdims%i_start,pdims%i_end,                         &
               pdims%j_start,pdims%j_end,                         &
               model_levels, at_extremity, u_on_p)

  CALL v_to_p (v, vdims_s%i_start,vdims_s%i_end,                        &
               vdims_s%j_start,vdims_s%j_end,                           &
               pdims%i_start,pdims%i_end,                               &
               pdims%j_start,pdims%j_end,                               &
               model_levels, at_extremity, v_on_p)

!$OMP END PARALLEL

END IF ! falliceshear_method

IF (l_co2_interactive) THEN
  DO k = 1, co2_dim_lev
    DO j = 1, co2_dim_row
      DO i = 1, co2_dim_len
        co2_3D(i,j,k) = co2(i,j,k)
      END DO
    END DO
  END DO
ELSE
  co2_3D(1,1,1) = 0.0
END IF

IF (L_mcr_qcf2) THEN  ! Second cloud ice variable in use
  ALLOCATE ( qcf2_n(row_length, rows, model_levels) )
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(i,j,k)                          &
!$OMP SHARED(model_levels,rows,row_length,qcf2_n,qcf2)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        qcf2_n(i,j,k) = qcf2(i, j, k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  ALLOCATE ( qcf2_n(1,1,1) )
END IF

IF (L_mcr_qrain) THEN  ! Prognostic rain in use
  ALLOCATE ( qrain_n(row_length, rows, model_levels) )
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(i,j,k)                          &
!$OMP SHARED(model_levels,rows,row_length,qrain_n,qrain)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        qrain_n(i,j,k) = qrain(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  ALLOCATE ( qrain_n(1,1,1) )
END IF

IF (L_mcr_qgraup) THEN  ! Prognostic graupel in use
  ALLOCATE ( qgraup_n(row_length, rows, model_levels) )
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(i,j,k)                          &
!$OMP SHARED(model_levels,rows,row_length,qgraup_n,qgraup)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        qgraup_n(i,j,k) = qgraup(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  ALLOCATE ( qgraup_n(1,1,1) )
END IF


IF (model_type == mt_lam .AND. L_lbc_old) THEN

  ! In the LAM set qcl_n, qcf_n and cloud fractions to zero on model
  ! boundaries. This avoids failures due to inconsistences in these fields.
  ! As the increments on the boundary are purely from the boundary
  ! conditions we are free to do anything sensible to these values.
  ! Note that this only applies to the variables in the physics.

  L_zero_boundaries=.TRUE.
  L_zero_halos=.FALSE.

  ! DEPENDS ON: zero_lateral_boundaries
  CALL zero_lateral_boundaries(                                         &
   row_length,rows,0,0,model_levels,fld_type_p,qcl_n,                   &
   1, at_extremity,                                                     &
   L_zero_boundaries,L_zero_halos)

  ! DEPENDS ON: zero_lateral_boundaries
  CALL zero_lateral_boundaries(                                         &
   row_length,rows,0,0,model_levels,fld_type_p,qcf_n,                   &
   1, at_extremity,                                                     &
   L_zero_boundaries,L_zero_halos)

  IF (L_mcr_qcf2)                                                       &
                  ! prognostic second cloud ice in use
  ! DEPENDS ON: zero_lateral_boundaries
           CALL zero_lateral_boundaries(                                &
            row_length,rows,0,0,model_levels,fld_type_p,qcf2_n,         &
            1, at_extremity,                                            &
            L_zero_boundaries,L_zero_halos)

  IF (L_mcr_qrain)                                                      &
                   ! prognostic rain in use
  ! DEPENDS ON: zero_lateral_boundaries
           CALL zero_lateral_boundaries(                                &
            row_length,rows,0,0,model_levels,fld_type_p,qrain_n,        &
            1, at_extremity,                                            &
            L_zero_boundaries,L_zero_halos)

  IF (L_mcr_qgraup)                                                     &
                    ! prognostic graupel in use
  ! DEPENDS ON: zero_lateral_boundaries
           CALL zero_lateral_boundaries(                                &
            row_length,rows,0,0,model_levels,fld_type_p,qgraup_n,       &
            1, at_extremity,                                            &
            L_zero_boundaries,L_zero_halos)

  ! DEPENDS ON: zero_lateral_boundaries
  CALL zero_lateral_boundaries(                                         &
   row_length,rows,0,0,model_levels,fld_type_p,                         &
   area_cloud_fraction,                                                 &
   1, at_extremity,                                                     &
   L_zero_boundaries,L_zero_halos)

  ! DEPENDS ON: zero_lateral_boundaries
  CALL zero_lateral_boundaries(                                         &
   row_length,rows,0,0,model_levels,fld_type_p,                         &
   cf_n,                                                                &
   1, at_extremity,                                                     &
   L_zero_boundaries,L_zero_halos)

  ! DEPENDS ON: zero_lateral_boundaries
  CALL zero_lateral_boundaries(                                         &
   row_length,rows,0,0,model_levels,fld_type_p,                         &
   cfl_n,                                                               &
   1, at_extremity,                                                     &
   L_zero_boundaries,L_zero_halos)

  ! DEPENDS ON: zero_lateral_boundaries
  CALL zero_lateral_boundaries(                                         &
   row_length,rows,0,0,model_levels,fld_type_p,                         &
   cff_n,                                                               &
   1, at_extremity,                                                     &
   L_zero_boundaries,L_zero_halos)


END IF !  model_type == mt_lam .and. L_lbc_old


! ----------------------------------------------------------------------
! Section Communications.
! ----------------------------------------------------------------------

SELECT CASE (model_type)

  !       If climatological aerosols are present, or if any diagnostics
  !       for fluxes at the tropopause are requested, or if any
  !       stratospheric aerosol optical depth diagnostics are requested,
  !       or if aerosol properties are to switch to sulphuric acid in
  !       the stratosphere, then it is necessary to find the tropopause.

  !       This sets a logical based different tests as SCM has no
  !       Stash

CASE (mt_single_column)
  l_RequestedDiags =                                                    &
    ( l_climat_aerosol .OR.                                             &
      ( l_ukca_radaer .AND. l_ukca_radaer_sustrat ) .OR.                &
      ( l_glomap_clim_radaer .AND. l_glomap_clim_radaer_sustrat ) )

CASE DEFAULT
  l_RequestedDiags =                                                    &
    ( l_climat_aerosol .OR.                                             &
      ( l_ukca_radaer .AND. l_ukca_radaer_sustrat ) .OR.                &
      ( l_glomap_clim_radaer .AND. l_glomap_clim_radaer_sustrat ) .OR.  &
      sf_calc(237,1) .OR. sf_calc(238,1) .OR.                           &
      sf_calc(237,2) .OR. sf_calc(238,2) .OR.                           &
      sf_calc(251,2) .OR. sf_calc(252,2) .OR.                           &
      sf_calc(253,2) .OR. sf_calc(254,2) .OR.                           &
      sf_calc(255,2) .OR. sf_calc(256,2) )

END SELECT

!     Find the tropopause level
!     -------------------------
IF ( (l_rad_step_diag .OR. l_rad_step_prog) .AND.                       &
     l_RequestedDiags ) THEN

  !       The variables min_trop_level and max_trop_level
  !       index rho-levels in the control routines,
  !       which start from a zeroth level at the surface.
  !       Layer boundaries in physical routines are indexed with
  !       the convention that the lowest is indexed by 1.
  !       Since the first rho-level is omitted from the physics
  !       grid, the two methods of indexing refer to the same
  !       horizontal level from the second rho-level upwards,
  !       so there is actually no need to adjust these variables
  !       for the change of indexing convention.

  !       Set SCM dummy values to zero
  scm_dummy_1d(:,:)   = 0.0
  scm_dummy_2d(:,:,:) = 0.0

  CALL tropin (t_n, exner_rho_levels,                                   &
    exner_theta_levels(:,:,1:tdims%k_end),                              &
    row_length, rows, model_levels, offx, offy,                         &
    at_extremity,scm_dummy_1d,scm_dummy_2d,                             &
    min_trop_level, max_trop_level, trindx )

END IF

CALL set_fsd_parameters(                                                &

!       Model dimensions.
        row_length, rows, n_cca_levels, offx, offy,                     &

!       Properties of clouds
        cca_dp, cca_md, cca_sh, xx_cos_theta_latitude,                  &

!       Model switches
        l_rad_step_diag, l_rad_step_prog, at_extremity,                 &

!       Error information
        Error_code  )

!     Generate sub-grid cloud field
!     -----------------------------
CALL open_cloud_gen (                                                   &

!       Parallel variables
        global_row_length, global_rows,                                 &
        mype,n_proc, at_extremity,                                      &

!       Model dimensions.
        row_length, rows,                                               &
        row_length*rows, p_layer_centres, offx, offy,                   &

!       Properties of clouds
        area_cloud_fraction, dp_corr_strat, cct, global_cloud_top,      &
        cca, n_cca_levels, ccw, qcl_n, qcf_n, xx_cos_theta_latitude,    &

!       Model switches
        l_rad_step_diag, l_rad_step_prog,                               &

!       Time stepping information.
        val_year, val_day_number, val_hour, val_minute, val_second,     &

!       Error information
        Error_code  )

!     Diagnostic RHcrit
!     -----------------
IF (i_rhcpt == rhcpt_horiz_var) THEN

  !       Dimension diagnostic 3D RHcrit array
  rhc_row_length = tdims%i_len
  rhc_rows       = tdims%j_len

  ALLOCATE(ext_p_layer_centres(0:rhc_row_length+1,0:rhc_rows+1,         &
                                         0:tdims%k_end))
  ALLOCATE(ext_tl(0:rhc_row_length+1, 0:rhc_rows+1, 1: tdims%k_end))
  ALLOCATE(ext_ql(0:rhc_row_length+1, 0:rhc_rows+1, 1: tdims%k_end))
  ALLOCATE(ext_qcf(0:rhc_row_length+1,0:rhc_rows+1, 1: tdims%k_end))
  ALLOCATE(ext_ice_frac(0:rhc_row_length+1,0:rhc_rows+1))
  ALLOCATE(ext_land_frac(0:rhc_row_length+1,0:rhc_rows+1))

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, tdims%k_end
    DO j = 1, rhc_rows
      DO i = 1, rhc_row_length

        ext_p_layer_centres(i,j,k) = p_layer_centres(i,j,k)

        IF (i_cld_vn == i_cld_pc2) THEN
          ext_tl(i,j,k) = t_n(i,j,k)                                    &
                          - (lc * qcl_n(i,j,k) ) / cp
          ext_ql(i,j,k) = q_n(i,j,k) + qcl_n(i,j,k)
        ELSE  ! i_cld_pc2
          ext_tl(i,j,k) = t_n(i,j,k)
          ext_ql(i,j,k) = q_n(i,j,k)
        END IF  ! i_cld_pc2
        ext_qcf(i,j,k) = qcf_n(i,j,k)

      END DO ! i
    END DO ! j
  END DO ! k
!$OMP END DO

  !       Rhc_rows_do2:
!$OMP DO SCHEDULE(STATIC)
  DO j = 1, rhc_rows

    !         Rhc_rowlen_do2:
    DO i = 1, rhc_row_length

      ext_p_layer_centres(i,j,0) = p_layer_centres(i,j,0)
      ext_ice_frac(i,j) = ice_fract(i,j)
      ext_land_frac(i,j) = 0.0

    END DO ! Rhc_rowlen_do2

  END DO ! Rhc_rows_do2
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, land_points

    j = (land_index(k)-1)/( tdims%i_len) + 1
    i = land_index(k) - (j-1)*( tdims%i_len)
    ext_land_frac(i,j) = fland(k)

  END DO
!$OMP END DO

!$OMP END PARALLEL
  !       Synchronize haloes.

  i_field = 0

  i_field = i_field + 1
  fields_to_swap(i_field) % field    => ext_p_layer_centres(:,:,:)
  fields_to_swap(i_field) % field_type  =  fld_type_p
  fields_to_swap(i_field) % levels      =  tdims%k_end+1
  fields_to_swap(i_field) % rows        =  rhc_rows
  fields_to_swap(i_field) % vector      =  .FALSE.

  i_field = i_field + 1
  fields_to_swap(i_field) % field       => ext_tl(:,:,:)
  fields_to_swap(i_field) % field_type  =  fld_type_p
  fields_to_swap(i_field) % levels      =  tdims%k_end
  fields_to_swap(i_field) % rows        =  rhc_rows
  fields_to_swap(i_field) % vector      =  .FALSE.

  i_field = i_field + 1
  fields_to_swap(i_field) % field       => ext_ql(:,:,:)
  fields_to_swap(i_field) % field_type  =  fld_type_p
  fields_to_swap(i_field) % levels      =  tdims%k_end
  fields_to_swap(i_field) % rows        =  rhc_rows
  fields_to_swap(i_field) % vector      =  .FALSE.

  i_field = i_field + 1
  fields_to_swap(i_field) % field       => ext_qcf(:,:,:)
  fields_to_swap(i_field) % field_type  =  fld_type_p
  fields_to_swap(i_field) % levels      =  tdims%k_end
  fields_to_swap(i_field) % rows        =  rhc_rows
  fields_to_swap(i_field) % vector      =  .FALSE.

  CALL swap_bounds_mv( fields_to_swap, i_field,                         &
                   rhc_row_length, 1, 1)

  i_field = 0
  i_field = i_field + 1
  fields_to_swap(i_field) % field_2d    => ext_land_frac(:,:)
  fields_to_swap(i_field) % field_type  =  fld_type_p
  fields_to_swap(i_field) % levels      =  1
  fields_to_swap(i_field) % rows        =  rhc_rows
  fields_to_swap(i_field) % vector      =  .FALSE.

  i_field = i_field + 1
  fields_to_swap(i_field) % field_2d    => ext_ice_frac(:,:)
  fields_to_swap(i_field) % field_type  =  fld_type_p
  fields_to_swap(i_field) % levels      =  1
  fields_to_swap(i_field) % rows        =  rhc_rows
  fields_to_swap(i_field) % vector      =  .FALSE.

  CALL swap_bounds_2d_mv( fields_to_swap, i_field,                      &
                   rhc_row_length, 1, 1)

ELSE ! i_rhcpt = rhcpt_off or rhcpt_tke_based

  !       RHcrit will be a 1D parametrized array input from user interface
  rhc_row_length = 1
  rhc_rows = 1

  ALLOCATE(ext_p_layer_centres(1,1,1))
  ALLOCATE(ext_tl(1,1,1))
  ALLOCATE(ext_ql(1,1,1))
  ALLOCATE(ext_qcf(1,1,1))
  ALLOCATE(ext_ice_frac(1,1))
  ALLOCATE(ext_land_frac(1,1))

END IF ! i_rhcpt


!     Set sea salt arrays
!     -------------------
sea_salt: IF (((l_use_seasalt_direct .OR. l_use_seasalt_indirect) .AND. &
           (l_rad_step_diag .OR. l_rad_step_prog)) .OR.                 &
            l_use_seasalt_autoconv) THEN


  ! calculate level 1 winds on p points using subroutine
  CALL uv_p_pnts(u(0:row_length-1,1:rows,1),v(1:row_length,0:n_rows-1,1),&
                 cos_theta_longitude,sin_theta_longitude,             &
                 global_row_length,gc_proc_row_group,u_1,v_1)

  ! Calculate level 1 windspeed (NB on rho levels)
  ! Adjust level 1 windspeed to be 10m windspeed

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, windspeed_1)              &
!$OMP SHARED(rows,row_length,u_1,v_1,z0_sea,p1,windspeed_10m,           &
!$OMP     height_theta)                                                 &
!$OMP SCHEDULE(STATIC)   
  DO j = 1, rows
    DO i = 1, row_length
      windspeed_1=SQRT(u_1(i,j)**2+v_1(i,j)**2)
      windspeed_10m(i,j)=windspeed_1*p1/(LOG(height_theta(i,j,1)/z0_sea))
    END DO
  END DO
!$OMP END PARALLEL DO

  CALL set_seasalt_4A(windspeed_10m, height_theta, land_fract,        &
                  ice_fract,row_length, rows, model_levels,           &
                  salt_dim1, salt_dim2, salt_dim3,                    &
                  bl_levels, sea_salt_film, sea_salt_jet)

ELSE ! l_use_seasalt_direct etc

  sea_salt_film(1, 1, 1) = 0.0
  sea_salt_jet( 1, 1, 1) = 0.0

END IF sea_salt ! l_use_seasalt_direct etc


!  ^^ Keep this last in the communications section
!     as call to set_seasalt has a load imbalance.

!-----------------------------------------------------------------------
! COSP initialization before calling mycrophysics and radiation
!-----------------------------------------------------------------------
l_cosp_call = .FALSE.
IF (l_cosp) THEN
  CALL cosp_init(l_cosp,l_rad_step_prog,l_radiation,row_length,         &
     rows,model_levels,cosp_crain_3d,cosp_csnow_3d,p,q_n,               &
     T_n,t_surf,p_star,p_theta_levels,r_theta_levels,r_rho_levels,      &
     l_cosp_call,l_cosp_run,cosp_npoints,cosp_cfg,                      &
     cosp_gbx,cosp_sgx,cosp_sgh)
ELSE
  CALL allocate_null_gbx()
END IF

! ----------------------------------------------------------------------
! Section Microphysics. Call microphys_ctl routine
! ----------------------------------------------------------------------

IF (Ltimer) CALL timer ('AP1 Microphys (AP1M)',5)
!
IF (L_rain) THEN
  !
  ! Only allow dynamic allocation of full space for arrays for 3D
  ! precipitation diagnostics if they are being used. Otherwise save space
  ! and give them a minimum size of 1 by 1.
  SELECT CASE (model_type)

  CASE (mt_single_column)
    lspice_dim1 = row_length
    lspice_dim2 = rows
    lspice_dim3 = model_levels

  CASE DEFAULT
    IF ( sf(222,4) .OR. sf(223,4) .OR. sf(224,4) .OR. sf(225,4) .OR.    &
         sf(226,4) .OR. sf(227,4) .OR. sf(323,4) .OR. l_dust    .OR.    &
         l_sulpc_so2 .OR. l_soot  .OR. l_biomass .OR. l_ocff    .OR.    &
         L_cosp_call) THEN
      lspice_dim1 = row_length
      lspice_dim2 = rows
      lspice_dim3 = model_levels
    ELSE
      lspice_dim1 = 1
      lspice_dim2 = 1
      lspice_dim3 = 1
    END IF

  END SELECT ! model_type

  IF ( l_casim ) THEN

    CALL casim_ctl (                                                   &
    ! Primary fields passed in
       T_n, q_n, qcl_n, qcf_n, qcf2_n, qrain_n, qgraup_n, w,           &
       p_layer_centres, rho, exner_theta_levels,                       &
       flash_pot, tracer,                                              &
       stashwork4, stashwork21,                                        &
    ! Increment fields passed in/out
       Theta_inc, q_inc, qcl_inc, qcf_inc,                             &
       qcf2_inc, qrain_inc, qgraup_inc,                                &
    ! Fields required elsewhere
       ls_rain, ls_snow, ls_graup, micro_tends, n_drop_pot )

  ELSE ! l_casim

    ! DEPENDS ON: microphys_ctl
    CALL microphys_ctl (                                                &

    ! Parallel variables
          halo_i, halo_j, offx, offy, global_row_length,                &
          at_extremity,                                                 &

  ! model dimensions.
          row_length, rows, rhc_row_length, rhc_rows, n_rows,           &
          land_points, bl_levels,                                       &
          lspice_dim1,lspice_dim2,lspice_dim3,                          &
          salt_dim1, salt_dim2, salt_dim3,                              &
          cdnc_dim1, cdnc_dim2, cdnc_dim3,                              &
          n_arcl_species, n_arcl_compnts, i_arcl_compnts,               &

  ! Model switches
          L_sulpc_so2, L_sulpc_nh3, L_soot, L_biomass, L_ocff, L_nitrate,&
          L_cosp_call,                                                  &

  ! Model parameters
          rhcrit,                                                       &

  ! Primary fields passed in
          T_n, T_n2, q_n, qcl_n, qcf_n, qcf2_n, qrain_n, qgraup_n,      &
          cf_n, cfl_n, cff_n,                                           &
          u_on_p, v_on_p, w,                                            &
          snow_depth,                                                   &
          land_sea_mask, ice_fract,                                     &
          p_layer_centres, p_layer_boundaries,                          &
          rho, aerosol, flash_pot,                                      &
          ukca_cdnc, easyaerosol_cdnc,                                  &
          dust_div1, dust_div2, dust_div3,                              &
          dust_div4, dust_div5, dust_div6,                              &
          so2, nh3, so4_aitken, so4_accu, so4_diss,                     &
          soot_agd, soot_cld, bmass_agd, bmass_cld,                     &
          ocff_agd, ocff_cld, nitr_acc, nitr_diss, biogenic,            &
          sea_salt_film, sea_salt_jet, arcl,                            &

  ! Other fields passed in
          ntml, cumulus,                                                &
          fland, land_index,                                            &
  ! For seeder feeder scheme
          sd_orog_land,                                                 &
  ! Variables for stochastic physics random parameters
          m_ci,                                                         &
  ! diagnostic info
          stashwork4, stashwork21,                                      &
  !
  ! SCM diagnostics switches (dummy in full UM)
          nSCMDpkgs, L_SCMDiags,                                        &

  ! Increment fields passed in/out
          Theta_inc, q_inc, qcl_inc, qcf_inc,                           &
          qcf2_inc, qrain_inc, qgraup_inc,                              &
          cf_inc, cfl_inc, cff_inc,                                     &

  ! Fields required elsewhere
          ls_rain, ls_rainfrac, ls_snow, ls_graup, micro_tends,         &
          cosp_gbx, cosp_frac_agg, n_drop_pot,                          &
  ! Field for Rh crit parametrization
          ext_p_layer_centres,                                          &
          ext_tl,                                                       &
          ext_ql,                                                       &
          ext_qcf,                                                      &
          ext_ice_frac,                                                 &
          ext_land_frac,                                                &

  ! BL w-variance required for mixed phase cloud calculation
          bl_w_var,                                                     &
  ! error information
          Error_code  )

  END IF ! l_casim

ELSE ! l_rain
  ls_rain(:,:) = 0.0
  ls_snow(:,:) = 0.0
  ls_graup(:,:)= 0.0
  micro_tends(:,:,:,:) = 0.0
END IF !l_rain

! Deallocate arrays used for RH crit parametrization
DEALLOCATE(ext_land_frac)
DEALLOCATE(ext_ice_frac)
DEALLOCATE(ext_qcf)
DEALLOCATE(ext_ql)
DEALLOCATE(ext_tl)
DEALLOCATE(ext_p_layer_centres)

IF (Ltimer) CALL timer ('AP1 Microphys (AP1M)',6)

!-------------------------------------------------
! check moisture conservation for microphysics
!-------------------------------------------------
IF (l_check_moist_inc .AND. model_type == mt_global) THEN

  ! ep4 holds total large-scale precipitation
  DO j=1,rows
    DO i=1,row_length
      ep4(i,j) = ls_rain(i,j) + ls_snow(i,j)
    END DO
  END DO

  ! total moist increment

  DO k=1,Model_levels
    DO j=1,rows
      DO i=1,row_length

        dqt(i,j,k) = q_inc(i,j,k)+ qcl_inc(i,j,k)+qcf_inc(i,j,k)
        IF (l_mcr_qrain ) THEN
          dqt(i,j,k) = dqt(i,j,k)+qrain_inc(i,j,k)
        END IF
        IF (l_mcr_qgraup ) THEN
          dqt(i,j,k) = dqt(i,j,k)+qgraup_inc(i,j,k)
        END IF
        IF (l_mcr_qcf2 ) THEN
          dqt(i,j,k) = dqt(i,j,k)+qcf2_inc(i,j,k)
        END IF
        dqt(i,j,k) = dqt(i,j,k)*recip_timestep
      END DO
    END DO
  END DO


  ! Need wet or dry density depending on specific or mixing ratios
  ! In this subroutine rho holds a wet density, unscaled_dry_rho
  ! a dry density

  IF (l_mr_physics) THEN
    DO k=1,Model_levels
      DO j=1,rows
        DO i=1,row_length
          rho_con(i,j,k) = unscaled_dry_rho(i,j,k)*r_rho_levels(i,j,k)  &
                                              *r_rho_levels(i,j,k)
        END DO
      END DO
    END DO
  ELSE   ! want wet density
    DO k=1,Model_levels
      DO j=1,rows
        DO i=1,row_length
          rho_con(i,j,k) = rho(i,j,k)
        END DO
      END DO
    END DO
  END IF    ! test on l_mr_physics

  scheme_name = 'Micro-physics'
  CALL check_dmoist_inc(row_length, rows, model_levels,                 &
                      delta_lambda,delta_phi,timestep,                  &
                      p_layer_boundaries,rho_con,                       &
                      dqt,                                              &
                      ep4, scheme_name,                                 &
                      ep1, ep2, ep3)
END IF

! ----------------------------------------------------------------------
! Section Turbulence. Call pc2_turbulence_ctl routine
! ----------------------------------------------------------------------
! Earlier versions of PC2 had the erosion term included here,
! parallel with the rest of the slow physics. We have since decided to
! move the erosion to be parallel with the PC2 response to convection
! since this results in improved numerical balances in shallow
! convection.
! In the case of PC2 using diagnostic shallow cloud it is better
! to do the erosion term here so the call is now under a switch for this.

IF (i_cld_vn == i_cld_pc2 .AND. l_pc2_diag_sh .AND. .NOT. l_micro_eros ) THEN

  ! Because of the switches used for this pc2 call, a dummy array is
  ! needed in place of a variable that is not used here.
  dummy_zeros(:,:,:)=0.0
  !DEPENDS ON: pc2_turbulence_ctl
  CALL pc2_turbulence_ctl (                                             &

  ! Primary fields passed in, unchanged on exit
            T_n, q_n, qcl_n, cf_n, cfl_n, cff_n                         &
        ,   p_layer_centres(1,1,1),                                     &
            dummy_zeros,                                                &
  ! diagnostic info
            STASHwork4                                                  &
  !
  ! SCM diagnostics switches (dummy in full UM)
        ,   nSCMDpkgs, L_SCMDiags                                       &

  ! Increment fields passed in/out, updated on exit
        ,   Theta_inc, q_inc, qcl_inc, cf_inc, cfl_inc, l_pc2_prod_qcl_mp)

END IF  ! i_cld_pc2 and l_pc2_diag_sh
!
! ----------------------------------------------------------------------
! Section RAD   Radiation scheme.
!               This incorporates radiation code for non-radiation
!               timesteps which is in CLDCTL1.dk in the UM.
!-----------------------------------------------------------------------

IF (Ltimer) CALL timer ('AP1 Radiation (AP1R)',5)
!  Set dimensions of mineral dust arrays for passing to rad_ctl
IF (l_dust) THEN
  dust_dim1 = row_length*rows
  dust_dim2 = model_levels
ELSE
  dust_dim1 = 1
  dust_dim2 = 1
END IF
!  Set dimensions of _SULPHATE arrays for passing to RAD_CTL
IF (l_sulpc_so2) THEN
  sulp_dim1 = rows*row_length
  sulp_dim2 = model_levels
ELSE
  sulp_dim1 = 1
  sulp_dim2 = 1
END IF
!  Set dimensions of soot arrays for passing to RAD_CTL
IF (l_soot) THEN
  soot_dim1 = rows*row_length
  soot_dim2 = model_levels
ELSE
  soot_dim1 = 1
  soot_dim2 = 1
END IF
!  Set dimensions of array for aerosol climatology for NWP
IF (n_arcl_species > 0) THEN
  arcl_dim1 = rows*row_length
  arcl_dim2 = model_levels
ELSE
  arcl_dim1 = 1
  arcl_dim2 = 1
END IF
!  Set dimensions of biomass arrays for passing to RAD_CTL
IF (l_biomass) THEN
  bmass_dim1 = rows*row_length
  bmass_dim2 = model_levels
ELSE
  bmass_dim1 = 1
  bmass_dim2 = 1
END IF
!  Set dimensions of biogenic array for passing to RAD_CTL
IF (l_use_biogenic) THEN
  biogenic_dim1 = rows*row_length
  biogenic_dim2 = model_levels
ELSE
  biogenic_dim1 = 1
  biogenic_dim2 = 1
END IF
!  Set dimensions of OCFF array for passing to RAD_CTL
IF (l_ocff) THEN
  ocff_dim1 = rows*row_length
  ocff_dim2 = model_levels
ELSE
  ocff_dim1 = 1
  ocff_dim2 = 1
END IF
!  Set dimensions of nitrate aerosol array for passing to RAD_CTL
IF (L_nitrate) THEN
  nitrate_dim1 = rows*row_length
  nitrate_dim2 = model_levels
ELSE
  nitrate_dim1 = 1
  nitrate_dim2 = 1
END IF
!  Set dimensions of UKCA_RADAER arrays for passing to RAD_CTL
IF ( l_ukca_radaer .OR. l_glomap_clim_radaer ) THEN
  ukca_dim1 = rows*row_length
  ukca_dim2 = model_levels
ELSE
  ukca_dim1 = 1
  ukca_dim2 = 1
END IF

!     Code to calculate mixing ratio of well-mixed greenhouse gases
!       if scenarios for their time variation have been prescribed.

IF ( clim_fcg_nyears(s_co2) > 0 ) THEN  !  CO2 level calculated
  CALL gas_calc ( co2_mmr,                                              &
    clim_fcg_nyears(s_co2),  clim_fcg_years(1,s_co2),                   &
    clim_fcg_levls(1,s_co2), clim_fcg_rates(1,s_co2),                   &
    lenscen, error_code)
  IF ( error_code /= 0 ) THEN
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END IF

IF ( clim_fcg_nyears(s_n2o) > 0 ) THEN  !  same for N2O
  CALL gas_calc ( n2ommr,                                               &
     clim_fcg_nyears(s_n2o),  clim_fcg_years(1,s_n2o),                  &
     clim_fcg_levls(1,s_n2o), clim_fcg_rates(1,s_n2o),                  &
     lenscen, error_code)
  IF ( error_code /= 0 ) THEN
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END IF

IF ( clim_fcg_nyears(s_ch4) > 0 ) THEN  !  CH4
  CALL gas_calc ( ch4mmr,                                               &
     clim_fcg_nyears(s_ch4),  clim_fcg_years(1,s_ch4),                  &
     clim_fcg_levls(1,s_ch4), clim_fcg_rates(1,s_ch4),                  &
     lenscen, error_code)
  IF ( error_code /= 0 ) THEN
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END IF

IF ( clim_fcg_nyears(s_cfc11) > 0 ) THEN   !  "CFC11"
  CALL gas_calc ( c11mmr,                                               &
     clim_fcg_nyears(s_cfc11),  clim_fcg_years(1,s_cfc11),              &
     clim_fcg_levls(1,s_cfc11), clim_fcg_rates(1,s_cfc11),              &
     lenscen, error_code)
  IF ( error_code /= 0 ) THEN
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END IF

IF ( clim_fcg_nyears(s_cfc12) > 0 ) THEN   !  "CFC12"
  CALL gas_calc ( c12mmr,                                               &
     clim_fcg_nyears(s_cfc12),  clim_fcg_years(1,s_cfc12),              &
     clim_fcg_levls(1,s_cfc12), clim_fcg_rates(1,s_cfc12),              &
     lenscen, error_code)
  IF ( error_code /= 0 ) THEN
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END IF

IF ( clim_fcg_nyears(s_cfc113) > 0 ) THEN  !  "CFC113"
  CALL gas_calc ( c113mmr,                                              &
     clim_fcg_nyears(s_cfc113),  clim_fcg_years(1,s_cfc113),            &
     clim_fcg_levls(1,s_cfc113), clim_fcg_rates(1,s_cfc113),            &
     lenscen, error_code)
  IF ( error_code /= 0 ) THEN
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END IF

IF ( clim_fcg_nyears(s_cfc114) > 0 ) THEN  !  "CFC114"
  CALL gas_calc ( c114mmr,                                              &
     clim_fcg_nyears(s_cfc114),  clim_fcg_years(1,s_cfc114),            &
     clim_fcg_levls(1,s_cfc114), clim_fcg_rates(1,s_cfc114),            &
     lenscen, error_code)
  IF ( error_code /= 0 ) THEN
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END IF

IF ( clim_fcg_nyears(s_hcfc22) > 0 ) THEN  !  "HCFC22"
  CALL gas_calc ( hcfc22mmr,                                            &
     clim_fcg_nyears(s_hcfc22),  clim_fcg_years(1,s_hcfc22),            &
     clim_fcg_levls(1,s_hcfc22), clim_fcg_rates(1,s_hcfc22),            &
     lenscen, error_code)
  IF ( error_code /= 0 ) THEN
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END IF

IF ( clim_fcg_nyears(s_hfc125) > 0 ) THEN  !  "HFC125"
  CALL gas_calc ( hfc125mmr,                                            &
     clim_fcg_nyears(s_hfc125),  clim_fcg_years(1,s_hfc125),            &
     clim_fcg_levls(1,s_hfc125), clim_fcg_rates(1,s_hfc125),            &
     lenscen, error_code)
  IF ( error_code /= 0 ) THEN
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END IF

IF ( clim_fcg_nyears(s_hfc134a) > 0 ) THEN  ! "HFC134a"
  CALL gas_calc ( hfc134ammr,                                           &
     clim_fcg_nyears(s_hfc134a),  clim_fcg_years(1,s_hfc134a),          &
     clim_fcg_levls(1,s_hfc134a), clim_fcg_rates(1,s_hfc134a),          &
     lenscen, error_code)
  IF ( error_code /= 0 ) THEN
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END IF

! Copy greenhouse gases into greenhouse gas array
ALLOCATE(grgas_field(row_length, rows, model_levels, ngrgas))

!$OMP PARALLEL DEFAULT(NONE) SHARED(ngrgas,grgas_addr,grgas_field,   &
!$OMP tracer_ukca,row_length,rows,model_levels,l_ukca_strat,         &
!$OMP l_ukca_strattrop,ozone_levels,c11mmr,c12mmr,ozone,ch4mmr,      &
!$OMP n2ommr,c113mmr,hcfc22mmr,q_n)                                  &
!$OMP PRIVATE(i,j,k)

! check to see if UKCA is specifying trace gases, and fill fields
!$OMP DO SCHEDULE(STATIC)
DO i=1,ngrgas
  ! Make sure the field if used is all positive.
  IF (grgas_addr(i) > 0) THEN
  
    grgas_field(:,:,:,i) =                   &
      tracer_ukca(1:row_length,1:rows,1:model_levels,grgas_addr(i))

    ! For some configurations CFC-11 is held as a lumped species, so
    !  rescale field by the ratio of CFC-11/field at the surface
    IF (i == p_f11 .AND.                    &
      (L_ukca_strat .OR. L_ukca_strattrop)) THEN
      DO k=model_levels,1,-1
        grgas_field(:,:,k,p_f11) = grgas_field(:,:,k,p_f11) *       &
          c11mmr / grgas_field(:,:,1,p_f11)
      END DO
    END IF

    ! For some configurations CFC-12 is held as a lumped species, so
    !  rescale field by the ratio of CFC-12/field at the surface
    IF (i == p_f12 .AND.                    &
      (L_ukca_strat .OR. L_ukca_strattrop)) THEN
      DO k=model_levels,1,-1
        grgas_field(:,:,k,p_f12) = grgas_field(:,:,k,p_f12) *       &
          c12mmr / grgas_field(:,:,1,p_f12)
      END DO
    END IF
    
    grgas_field(:,:,:,i) = MAX(grgas_field(:,:,:,i),0.0)
  ELSE
    ! Set to RMDI if not used.
    grgas_field(:,:,:,i) = rmdi
  END IF
END DO    ! Loop over grgas species
!$OMP END DO

IF (grgas_addr(p_o3) <= 0) THEN                                        
!$OMP DO SCHEDULE(STATIC)
  DO k=1,ozone_levels
    DO j = 1, rows
      DO i = 1, row_length
        grgas_field(i,j,k,p_o3) = ozone(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (grgas_addr(p_ch4) <= 0) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k=1,model_levels
    DO j = 1, rows
      DO i = 1, row_length
        grgas_field(i,j,k,p_ch4) = ch4mmr          ! global constant
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (grgas_addr(p_n2o) <= 0) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k=1,model_levels
    DO j = 1, rows
      DO i = 1, row_length
        grgas_field(i,j,k,p_n2o) = n2ommr          ! global constant
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (grgas_addr(p_f11) <= 0) THEN                                     
!$OMP DO SCHEDULE(STATIC)
  DO k=1,model_levels
    DO j = 1, rows
      DO i = 1, row_length
        grgas_field(i,j,k,p_f11) = c11mmr          ! global constant
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (grgas_addr(p_f12) <= 0) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k=1,model_levels
    DO j = 1, rows
      DO i = 1, row_length
        grgas_field(i,j,k,p_f12) = c12mmr          ! global constant
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (grgas_addr(p_f113) <= 0) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k=1,model_levels
    DO j = 1, rows
      DO i = 1, row_length
        grgas_field(i,j,k,p_f113) = c113mmr        ! global constant
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (grgas_addr(p_f22) <= 0) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k=1,model_levels
    DO j = 1, rows
      DO i = 1, row_length
        grgas_field(i,j,k,p_f22) = hcfc22mmr       ! global constant
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

!$OMP DO SCHEDULE(STATIC)
DO k=1,model_levels
  DO j = 1, rows
    DO i = 1, row_length
      grgas_field(i,j,k,p_h2os) =                          &
        q_n(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

!  Write trace gas mixing ratios to module for use in UKCA. The value can be
!  either that which was passed to this routine or the value calculated from
!  the time interpolation, depending on whether L_CLMCHFCG is true or false.
IF (L_ukca .AND.  (((i_ukca_scenario == i_ukca_scenario_um) .OR.        &
    L_ukca_set_trace_gases) .OR. L_ukca_prescribech4)) THEN
  IF (L_ukca_prescribech4 .AND. (ch4mmr == rmdi)) THEN
    cmessage = 'Missing value for ch4_mix_ratio'
    Error_code = 1
    CALL umPrint(                                                       &
        'Set value in Gas MMRs subsection of run_radiation in gui ',    &
        src='atmos_physics1')
    CALL ereport('Atmos_Physics1',Error_code,cmessage)
  END IF
  CALL ukca_set_trace_gas_mixratio(                                     &
        ch4mmr, co2_mmr, n2ommr, o2mmr,                                 &
        c11mmr, c12mmr, c113mmr, c114mmr,                               &
        hcfc22mmr, hfc125mmr, hfc134ammr)
END IF

!     Set number of cloud droplets to be used for the 1st indirect effect
IF (l_consistent_cdnc) THEN
  !       Use the n_drop_pot array from microphysics in radiation
  l_use_ndrop = .TRUE.
ELSE IF (l_ukca_aie1 .OR. l_glomap_clim_aie1) THEN
  !       Use the UKCA cdnc array in radiation
  n_drop_pot = ukca_cdnc
  l_use_ndrop = .TRUE.
END IF

IF (l_radiation) THEN

  CALL rad_ctl (                                                        &

  ! Parallel variables
          offx, offy, at_extremity, n_proc, global_cloud_top            &

  ! model dimensions.
        , row_length, rows                                              &
        , bl_levels, Ozone_levels, cloud_levels, N_cca_levels           &
        , ntiles, land_points, dust_dim1, dust_dim2                     &
        , biogenic_dim1                                                 &
        , biogenic_dim2, sulp_dim1, sulp_dim2, soot_dim1, soot_dim2     &
        , bmass_dim1, bmass_dim2, salt_dim1, salt_dim2, salt_dim3       &
        , co2_dim_len, co2_dim_row, co2_dim_lev, arcl_dim1, arcl_dim2   &
        , n_arcl_species, n_arcl_compnts, i_arcl_compnts                &
        , ocff_dim1, ocff_dim2, nitrate_dim1, nitrate_dim2              &
        , ukca_dim1, ukca_dim2, ukca_radaer%n_mode, ukca_radaer%n_cpnt  &

  ! Model switches
        , L_emcorr                                                      &
        , L_snow_albedo, l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup          &
        , l_murk_rad                                                    &
        , L_co2_interactive                                             &
        , l_use_arcl                                                    &
        , l_dust, l_sulpc_so2, l_soot, l_biomass, l_ocff, l_nitrate     &
        , l_cosp_call                                                   &
        , l_easyaerosol_sw, l_easyaerosol_lw, l_easyaerosol_cdnc        &

  ! in coordinate information
        , rho, trindx                                                   &

  ! in time stepping information.
        , timestep_number, previous_time                                &

  ! diagnostic info
             ,                                                          &
         STASHwork1,                                                    &
         STASHwork2                                                     &
  !
  ! SCM diagnostics switches (dummy in full UM)
        , nSCMDpkgs, L_SCMDiags                                         &

  ! in data fields.
        , p_star                                                        &
        , p_layer_boundaries, p_layer_centres                           &
        , land_sea_mask, fland, land0p5, sd_orog_land                   &
        , t_surf, tstar_sea, tstar_sice_cat, area_cloud_fraction        &
        , dust_div1(tdims_s%i_start,tdims_s%j_start,1)                  &
        , dust_div2(tdims_s%i_start,tdims_s%j_start,1)                  &
        , dust_div3(tdims_s%i_start,tdims_s%j_start,1)                  &
        , dust_div4(tdims_s%i_start,tdims_s%j_start,1)                  &
        , dust_div5(tdims_s%i_start,tdims_s%j_start,1)                  &
        , dust_div6(tdims_s%i_start,tdims_s%j_start,1)                  &
        , biogenic, so4_aitken(tdims_s%i_start,tdims_s%j_start,1)       &
        , so4_accu(tdims_s%i_start,tdims_s%j_start,1)                   &
        , so4_diss(tdims_s%i_start,tdims_s%j_start,1)                   &
        , soot_new(tdims_s%i_start,tdims_s%j_start,1)                   &
        , soot_agd(tdims_s%i_start,tdims_s%j_start,1)                   &
        , bmass_new(tdims_s%i_start,tdims_s%j_start,1)                  &
        , bmass_agd(tdims_s%i_start,tdims_s%j_start,1)                  &
        , bmass_cld(tdims_s%i_start,tdims_s%j_start,1)                  &
        , ocff_new(tdims_s%i_start,tdims_s%j_start,1)                   &
        , ocff_agd(tdims_s%i_start,tdims_s%j_start,1)                   &
        , ocff_cld(tdims_s%i_start,tdims_s%j_start,1)                   &
        , nitr_acc(tdims_s%i_start,tdims_s%j_start,1)                   &
        , nitr_diss(tdims_s%i_start,tdims_s%j_start,1)                  &
        , aerosol(1:row_length, 1:rows, 1:model_levels), arcl, ukca_radaer&
        , sea_salt_film, sea_salt_jet, ws_10m_sea, co2_3D               &
        , frac_control, n_drop_pot                                      &
        , easyaerosol_sw, easyaerosol_lw, easyaerosol_cdnc              &
  ! chemical greenhouse gas fields
        , ngrgas, grgas_field                                           &

  ! ancillary fields and fields needed to be kept from timestep to
  ! timestep
        , snow_depth, snow_depth_sea_cat, ice_fract, ice_fract_cat      &
        , ice_thick_cat, chloro_sea, rgrain, soot, canht                &
        , cca, ccb, cct, cclwp, ccw, lcbase                             &
  ! chemical ozone; replaces climatological ozone
        , grgas_field(:,:,1:ozone_levels,p_o3)                          &
        , SW_incs, LW_incs                                              &
        , O3_trop_level, O3_trop_height                                 &
        , T_trop_level, T_trop_height, zh, land_index, albsoil, albobs_sw&
        , albobs_vis, albobs_nir, lai, snow_tile, tile_frac, tstar_tile &
        , z0_tile, dOLR_rts, LW_down, SW_tile_rts                       &
        , land_alb,sice_alb                                             &
        , es_space_interp, rad_mask                                     &

  ! in/out
  ! chemical water replaces hydrological water
        , T_n, grgas_field(:,:,1:model_levels,p_h2os)                   &
        , qcl_n, qcf_n, cf_n, cfl_n, cff_n                              &
        , qcf2_n, qrain_n, qgraup_n                                     &
        , Theta_inc, q_inc, qcl_inc, cf_inc, cfl_inc                    &
        , sum_eng_fluxes, cos_zenith_angle                              &

  ! out.
        , photosynth_act_rad, rad_hr, dOLR, SW_tile                     &
  ! COSP arguments
        , cosp_gbx, cosp_sgx, cosp_sgh                                  &
  ! error information
        , Error_code  )

  ! Check error condition
  IF (Error_code > 0) THEN

    CALL ereport(RoutineName, Error_code,                               &
      "Error on return from radiation code.")
  END IF

ELSE
  photosynth_act_rad = 0.0
  rad_hr = 0.0
  dolr = 0.0
  sw_tile = 0.0
END IF

IF (Ltimer) CALL timer ('AP1 Radiation (AP1R)',6)

CALL close_cloud_gen

DEALLOCATE(grgas_field)

!-----------------------------------------------------------------------
!     Call to COSP
!-----------------------------------------------------------------------

IF (L_cosp_call) THEN
  CALL cosp_main(at_extremity,L_cosp_run,row_length,rows,               &
          n_rows,model_levels,qrain,qgraup,cosp_frac_agg(1,1,1),        &
          cosp_cfg,cosp_gbx,cosp_sgx,cosp_sgh,                          &
          STASHwork2)
  CALL free_cosp_gridbox(cosp_gbx)
  CALL free_cosp_subgrid(cosp_sgx)
  CALL free_cosp_sghydro(cosp_sgh)
END IF

IF (.NOT. l_cosp) THEN
  CALL free_null_gbx()
END IF

! ----------------------------------------------------------------------
! Set Coastal tiling dependent prognostics for output to D1 array:
! ----------------------------------------------------------------------

IF (l_ctile) THEN
  DO j = 1, rows
    DO i = 1, row_length
      land_alb_ctile(i,j)=land_alb(i,j)
      sice_alb_ctile(i,j)=sice_alb(i,j)
    END DO
  END DO
END IF

IF (Error_code  ==  0) THEN

  ! Convert temperature held in t_n to potential temperature.
  ! Convert increment after call to gravity wave drag, to allow for
  ! inclusion of dissipative heating term.
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k)                        &
!$OMP SHARED(tdims,T_n, Theta)                                          &
!$OMP SCHEDULE(STATIC)   
  DO k =             1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        T_n(i,j,k) = Theta(i,j,k)
      END DO
    END DO
  END DO
!$END PARALLEL DO
END IF ! on error code equal to zero

! ----------------------------------------------------------------------
! Section CNV.2 Energy correction code
! ----------------------------------------------------------------------

IF (Ltimer) CALL timer ('AP1 Conv Eng Corr',5)
IF (L_emcorr .AND. Error_code  ==  0) THEN

  ! Add convective + large scale, rain and snow, at the surface to the
  ! diabatic heating for use in the energy correction
  ! procedure.
  ! Scale variables by conversion factor so that only one call is required

  lclf = lc + lf

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, j)                           &
!$OMP SHARED(tot_precip_scaled_1,tot_precip_scaled_2, ls_rain,       &
!$OMP ls_snow,lclf,rows,row_length)                            

!$OMP DO SCHEDULE(STATIC)   
  DO j = 1, rows
    DO i = 1, row_length
      tot_precip_scaled_1(i,j) = ls_rain(i,j) * lc +              &
                                 ls_snow(i,j) * lclf
    END DO
  END DO
!$OMP END DO NOWAIT

  ! moist fluxes
!$OMP DO SCHEDULE(STATIC)   
  DO j = 1, rows
    DO i = 1, row_length
      tot_precip_scaled_2(i,j) = -ls_rain(i,j) -ls_snow(i,j)
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

  CALL flux_diag(tot_precip_scaled_1, xx_cos_theta_latitude,      &
                 row_length, rows ,offx, offy, 1.0,               &
                 sum_eng_fluxes,timestep)

  CALL flux_diag(tot_precip_scaled_2, xx_cos_theta_latitude,      &
                 row_length, rows ,offx, offy, 1.0,               &
                 sum_moist_flux,timestep)

END IF   ! L_emcorr

IF (Ltimer) CALL timer ('AP1 Conv Eng Corr',6)



! ----------------------------------------------------------------------
! Section GWD.1 Call gravity wave drag
! l_gwd      = orographic GWD
! l_use_ussp = middle atmosphere non-orographic GWD scheme
! ----------------------------------------------------------------------

IF ( l_gwd .OR. l_use_ussp ) THEN

  IF (error_code  ==  0 ) THEN
    IF (Ltimer) CALL timer ('AP1 G-wave drag',5)

    ! reset p at layer boundaries for spectral GWD code
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, j, k)                        &
!$OMP SHARED(model_levels,rows,row_length, p_theta_levels,           &
!$OMP p_layer_boundaries)                            

!$OMP DO SCHEDULE(STATIC)   
    DO k = 1, model_levels - 1
      DO j = 1, rows
        DO i = 1, row_length
          p_layer_boundaries(i,j,k) = p_theta_levels(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)   
    DO j = 1, rows
      DO i = 1, row_length
        p_layer_boundaries(i,j,model_levels) = 0.0
      END DO
    END DO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

    CALL NI_gwd_ctl(                                                &
                     halo_i, halo_j, offx, offy                     &
  ,  global_row_length,n_proc, n_procy, gc_proc_row_group           &
  ,                  at_extremity, neighbour                        &

    ! model dimensions.
          ,                  row_length, rows, n_rows, land_points      &
          ,                  model_levels                               &
    ! trig arrays
          ,                  sin_theta_longitude, sin_theta_latitude    &

    ! in coordinate information
          ,                  delta_lambda,delta_phi,true_latitude       &
          ,                  exner_theta_levels                         &

    ! in time stepping information.
          ,                  timestep, timestep_number                  &

    ! diagnostic info
               ,                                                        &
           STASHwork6                                                   &
    ! SCM diagnostics (dummy in full UM)
          ,                  nSCMDpkgs, L_SCMDiags                      &
    ! in data fields.
          ,                  u, v, totalppn, land_sea_mask              &
          ,                  p_layer_boundaries                         &
          ,                  rho, t_n, sd_orog_land                     &
          ,                  orog_grad_xx_land, orog_grad_xy_land       &
          ,                  orog_grad_yy_land, land_index              &

    ! in/out
          ,                  u_inc, v_inc,Theta_inc                     &

    ! error information
          , Error_code  )

    IF (Ltimer) CALL timer ('AP1 G-wave drag',6)

  END IF

END IF



IF (Error_code  ==  0) THEN

! Now convert temperature increment to theta increment
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k)                        &
!$OMP SHARED(tdims,Theta_inc, exner_theta_levels)                       &
!$OMP SCHEDULE(STATIC)   
  DO k =             1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        Theta_inc(i,j,k) = Theta_inc(i,j,k) /                           &
                           exner_theta_levels(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

END IF! on error code equal to zero


!-----------------------------------------------------------------------
! Tracer Source and Boundary updating where applicable
!-----------------------------------------------------------------------
! Aerosol
IF ( L_murk_source ) THEN

  DO level = 1, model_levels
    aerosol_em_scaled(:,:) = murk_source_scale * aerosol_em(:,:,level) 
    CALL trsrce(                                                        &
       rows, row_length, offx, offy, 0, 0                               &
,      theta, q_n, qcl_n, qcf_n, exner_rho_levels, rho                  &
,      aerosol(:,:,level), aerosol_em_scaled                            &
,      level, timestep, val_hour, val_minute, amp                       &
)
  END DO
END IF

IF ( L_murk_bdry ) THEN
  ! DEPENDS ON: trbdry
  CALL trbdry(                                                          &
       row_length, rows, n_rows, model_levels                           &
,      offx, offy, at_extremity                                         &
,      p, u, v                                                          &
,      aerosol, timestep                                                &
)
END IF

IF ( L_murk ) THEN
  IF (.NOT. L_murk_bdry) THEN
    ! Set the external halos if not specified by the
    ! formula for UK mes.

    ! DEPENDS ON: swap_bounds
    CALL Swap_Bounds(                                                   &
         aerosol, row_length, rows,                                     &
         model_levels, offx, offy, fld_type_p, swap_field_is_scalar)

  END IF

END IF

!---------------------------------------------------------------------
! Section METHOX Call methane oxidation
!---------------------------------------------------------------------

IF (l_use_methox) THEN

  IF (Ltimer) CALL timer ('AP1 NI_methox',5)

  !-------------------------------------------------
  ! check moisture conservation for methane oxidation
  !-------------------------------------------------
  IF (l_check_moist_inc .AND. model_type == mt_global) THEN

    ALLOCATE (q_methane(row_length, rows, model_levels) )

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)         &
!$OMP SHARED(model_levels,rows, row_length,q_inc, q_methane)
    DO k=1,Model_levels
      DO j=1,rows
        DO i=1,row_length
          q_methane(i,j,k) = q_inc(i,j,k)   ! copy before
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  IF (error_code == 0) THEN
     
    !    call methane oxidation directly (no interface routine needed)
     
    ! DEPENDS ON: ni_methox
    CALL NI_methox( row_length, rows, eta_theta_levels, timestep,       &
                    STASHwork4, q_n, q_inc )

  END IF

  !-------------------------------------------------
  ! check moisture conservation for methane oxidation
  !-------------------------------------------------
  IF (l_check_moist_inc .AND. model_type == mt_global) THEN
    ! no diagnostic value of source/sink
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,k) SHARED(rows,row_length,       &
!$OMP model_levels,ep4,q_inc,q_methane,dqt,recip_timestep,l_mr_physics, &
!$OMP rho_con,unscaled_dry_rho,r_rho_levels,rho)
!$OMP DO SCHEDULE(STATIC)
    DO j=1,rows
      DO i=1,row_length
        ep4(i,j) = 0.0
      END DO
    END DO
!$OMP END DO
!$OMP DO SCHEDULE(STATIC)
    DO k=1,Model_levels
      DO j=1,rows
        DO i=1,row_length
          ! increment from methane oxidation to q
          dqt(i,j,k) = (q_inc(i,j,k) - q_methane(i,j,k))*recip_timestep
        END DO
      END DO
    END DO
!$OMP END DO

    IF (l_mr_physics) THEN
!$OMP DO SCHEDULE(STATIC)
      DO k=1,Model_levels
        DO j=1,rows
          DO i=1,row_length
            rho_con(i,j,k) = unscaled_dry_rho(i,j,k)*r_rho_levels(i,j,k)&
                                              *r_rho_levels(i,j,k)
          END DO
        END DO
      END DO
!$OMP END DO
    ELSE   ! want wet density
!$OMP DO SCHEDULE(STATIC)
      DO k=1,Model_levels
        DO j=1,rows
          DO i=1,row_length
            rho_con(i,j,k) = rho(i,j,k)
          END DO
        END DO
      END DO
!$OMP END DO
    END IF

!$OMP END PARALLEL

    scheme_name = 'Methane oxidation'
    CALL check_dmoist_inc(row_length, rows, model_levels,               &
                      delta_lambda,delta_phi,timestep,                  &
                      p_layer_boundaries,rho_con,                       &
                      dqt,                                              &
                      ep4, scheme_name,                                 &
                      ep1, ep2, ep3)

    DEALLOCATE (q_methane)

  END IF  ! moisture check

  IF (Ltimer) CALL timer ('AP1 NI_methox',6)
END IF

! ----------------------------------------------------------------------
! Copy increment arrays into star locations
! ----------------------------------------------------------------------

! In the LAM set physics increments on boundary to zero.
IF (model_type == mt_lam .AND. L_lbc_old) THEN

  L_zero_boundaries=.TRUE.
  L_zero_halos=.FALSE.

  ! DEPENDS ON: zero_lateral_boundaries
  CALL zero_lateral_boundaries(                                         &
   row_length,rows,0,0,model_levels,fld_type_p,theta_inc,               &
   1, at_extremity,                                                     &
   L_zero_boundaries,L_zero_halos)

  ! DEPENDS ON: zero_lateral_boundaries
  CALL zero_lateral_boundaries(                                         &
   row_length,rows,0,0,model_levels,fld_type_p,q_inc,                   &
   1, at_extremity,                                                     &
   L_zero_boundaries,L_zero_halos)

  ! DEPENDS ON: zero_lateral_boundaries
  CALL zero_lateral_boundaries(                                         &
   row_length,rows,0,0,model_levels,fld_type_p,qcl_inc,                 &
   1, at_extremity,                                                     &
   L_zero_boundaries,L_zero_halos)

  ! DEPENDS ON: zero_lateral_boundaries
  CALL zero_lateral_boundaries(                                         &
   row_length,rows,0,0,model_levels,fld_type_p,qcf_inc,                 &
   1, at_extremity,                                                     &
   L_zero_boundaries,L_zero_halos)

  IF (L_mcr_qcf2)                                                       &
                  ! prognostic second cloud ice in use
  ! DEPENDS ON: zero_lateral_boundaries
           CALL zero_lateral_boundaries(                                &
            row_length,rows,0,0,model_levels,fld_type_p,qcf2_inc,       &
            1, at_extremity,                                            &
            L_zero_boundaries,L_zero_halos)

  IF (L_mcr_qrain)                                                      &
                   ! prognostic rain in use
  ! DEPENDS ON: zero_lateral_boundaries
           CALL zero_lateral_boundaries(                                &
            row_length,rows,0,0,model_levels,fld_type_p,qrain_inc,      &
            1, at_extremity,                                            &
            L_zero_boundaries,L_zero_halos)

  IF (L_mcr_qgraup)                                                     &
                    ! prognostic graupel in use
  ! DEPENDS ON: zero_lateral_boundaries
           CALL zero_lateral_boundaries(                                &
            row_length,rows,0,0,model_levels,fld_type_p,qgraup_inc,     &
            1, at_extremity,                                            &
            L_zero_boundaries,L_zero_halos)

  ! DEPENDS ON: zero_lateral_boundaries
  CALL zero_lateral_boundaries(                                         &
   row_length,rows,0,0,model_levels,fld_type_p,cf_inc,                  &
   1, at_extremity,                                                     &
   L_zero_boundaries,L_zero_halos)

  ! DEPENDS ON: zero_lateral_boundaries
  CALL zero_lateral_boundaries(                                         &
   row_length,rows,0,0,model_levels,fld_type_p,cfl_inc,                 &
   1, at_extremity,                                                     &
   L_zero_boundaries,L_zero_halos)

  ! DEPENDS ON: zero_lateral_boundaries
  CALL zero_lateral_boundaries(                                         &
   row_length,rows,0,0,model_levels,fld_type_p,cff_inc,                 &
   1, at_extremity,                                                     &
   L_zero_boundaries,L_zero_halos)

  ! DEPENDS ON: zero_lateral_boundaries
  CALL zero_lateral_boundaries(                                         &
   row_length,n_ROWS,offx,offy,model_levels,fld_type_v,v_inc,           &
   1, at_extremity,                                                     &
   L_zero_boundaries,L_zero_halos)

  ! DEPENDS ON: zero_lateral_boundaries
  CALL zero_lateral_boundaries(                                         &
   row_length,rows,offx,offy,model_levels,fld_type_u,u_inc,             &
   1, at_extremity,                                                     &
   L_zero_boundaries,L_zero_halos)

END IF  !  model_type == mt_lam .and. L_lbc_old



IF (i_cld_vn == i_cld_pc2) THEN
  !
  ALLOCATE( T_inc_diag  (row_length,rows,model_levels) )
  ALLOCATE( q_inc_diag  (row_length,rows,model_levels) )
  ALLOCATE( qcl_inc_diag(row_length,rows,model_levels) )
  ALLOCATE( qcf_inc_diag(row_length,rows,model_levels) )
  ALLOCATE( cf_inc_diag (row_length,rows,model_levels) )
  ALLOCATE( cfl_inc_diag(row_length,rows,model_levels) )
  ALLOCATE( cff_inc_diag(row_length,rows,model_levels) )

  !
  ! Now call a consistency check for the moisture and cloud fields.
  ! The catch is that _inc variables hold the increments and not the
  ! full variables yet, so temporarily form them as _inc, check them,
  ! then rewrite them as increments.
  !
  ! 1. Create full variables (not increments).
  !    (theta_inc and T_n hold potential temperature)
  !
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k)                   &
!$OMP SHARED(model_levels,rows,row_length,theta_inc,t_n,           &
!$OMP    exner_theta_levels,q_inc,q_n,qcl_inc,qcl_n,qcf_inc,qcf_n, &
!$OMP    cf_inc,cf_n,cfl_inc,cfl_n,cff_inc,cff_n,t_inc_diag,       &
!$OMP    q_inc_diag,qcl_inc_diag,qcf_inc_diag,cf_inc_diag,         &
!$OMP    cfl_inc_diag,cff_inc_diag)                                &
!$OMP SCHEDULE(STATIC)   
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        theta_inc(i,j,k) = (t_n(i,j,k)   + theta_inc(i,j,k))            &
                           * exner_theta_levels(i,j,k)
        ! theta_inc now holds temperature, t_n still holds potential temp
        q_inc(i,j,k)     = q_n(i,j,k)   + q_inc(i,j,k)
        qcl_inc(i,j,k)   = qcl_n(i,j,k) + qcl_inc(i,j,k)
        qcf_inc(i,j,k)   = qcf_n(i,j,k) + qcf_inc(i,j,k)
        cf_inc(i,j,k)    = cf_n(i,j,k)  + cf_inc(i,j,k)
        cfl_inc(i,j,k)   = cfl_n(i,j,k) + cfl_inc(i,j,k)
        cff_inc(i,j,k)   = cff_n(i,j,k) + cff_inc(i,j,k)

        t_inc_diag(i,j,k)   = theta_inc(i,j,k)
        q_inc_diag(i,j,k)   = q_inc(i,j,k)
        qcl_inc_diag(i,j,k) = qcl_inc(i,j,k)
        qcf_inc_diag(i,j,k) = qcf_inc(i,j,k)
        cf_inc_diag(i,j,k)  = cf_inc(i,j,k)
        cfl_inc_diag(i,j,k) = cfl_inc(i,j,k)
        cff_inc_diag(i,j,k) = cff_inc(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  !
  ! 2. Call the checking routine (needs to use temperature, not
  !    potential temperature)
  !
  CALL pc2_checks(p_layer_centres(1,1,1)                                &
                 ,theta_inc, cf_inc, cfl_inc                            &
                 ,cff_inc, q_inc, qcl_inc, qcf_inc                      &
                 ,l_mr_physics)
  !
  ! Now form the increment diagnostics in t_inc_diag etc.
  !
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k)                   &
!$OMP SHARED(model_levels,rows,row_length,theta_inc, q_inc,        &
!$OMP    qcl_inc,qcf_inc, cf_inc,cfl_inc,cff_inc,t_inc_diag,       &
!$OMP    q_inc_diag,qcl_inc_diag,qcf_inc_diag,cf_inc_diag,         &
!$OMP    cfl_inc_diag,cff_inc_diag)                                &
!$OMP SCHEDULE(STATIC)   
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        t_inc_diag(i,j,k)   = theta_inc(i,j,k) - t_inc_diag(i,j,k)
        q_inc_diag(i,j,k)   = q_inc(i,j,k)     - q_inc_diag(i,j,k)
        qcl_inc_diag(i,j,k) = qcl_inc(i,j,k)   - qcl_inc_diag(i,j,k)
        qcf_inc_diag(i,j,k) = qcf_inc(i,j,k)   - qcf_inc_diag(i,j,k)
        cf_inc_diag(i,j,k)  = cf_inc(i,j,k)    - cf_inc_diag(i,j,k)
        cfl_inc_diag(i,j,k) = cfl_inc(i,j,k)   - cfl_inc_diag(i,j,k)
        cff_inc_diag(i,j,k) = cff_inc(i,j,k)   - cff_inc_diag(i,j,k)
      END DO ! i
    END DO ! j
  END DO ! k
!$OMP END PARALLEL DO
  !
  ! 3. Update fields to produce the net increments which are written to
  !    the _star variables.
  !

  SELECT CASE (model_type)
  CASE (mt_single_column)
    !-----------------------------------------------------------------------
    !     SCM PC2 Diagnostics Package
    !-----------------------------------------------------------------------
    IF (L_SCMDiags(SCMDiag_PC2)) THEN

      !           Stash 4,141
      CALL SCMoutput(T_inc_diag,'dt_pc2ck',                             &
           'Temperature inc PC2 checks in AP1','K',                     &
            t_avg,d_all,default_streams,'',RoutineName)

      !           Stash 4,142
      CALL SCMoutput(q_inc_diag,'dq_pc2ck',                             &
           'Specific humidity inc PC2 checks in AP1','kg/kg',           &
            t_avg,d_wet,default_streams,'',RoutineName)

      !           Stash 4,143
      CALL SCMoutput(qcl_inc_diag,                                      &
           'dqcl_pc2ck','QCL inc PC2 checks in AP1','kg/kg',            &
           t_avg,d_wet,default_streams,'',RoutineName)

      !           Stash 4,144
      CALL SCMoutput(qcf_inc_diag,                                      &
           'dqcf_pc2ck','QCF inc PC2 checks in AP1','kg/kg',            &
           t_avg,d_wet,default_streams,'',RoutineName)

      !           Stash 4,152
      CALL SCMoutput(cf_inc_diag,'dbcf_pc2ck',                          &
           'Bulk cloud fraction inc PC2 checks in AP1','Fraction',      &
            t_avg,d_wet,default_streams,'',RoutineName)

      !           Stash 4,153
      CALL SCMoutput(cfl_inc_diag,'dcfl_pc2ck',                         &
           'Liquid cloud fraction inc PC2 checks in AP1','Fraction',    &
           t_avg,d_wet,default_streams,'',RoutineName)

      !           Stash 4,154
      CALL SCMoutput(cff_inc_diag,'dcff_pc2ck',                         &
           'Ice cloud fraction inc PC2 checks in AP1','Fraction',       &
           t_avg,d_wet,default_streams,'',RoutineName)

    END IF ! L_SCMDiags(SCMDiag_PC2)

  CASE DEFAULT

    !
    ! 4. Call diagnostic writing subroutine.
    ! Not called for SCM.
    !
    ! Check that checking diagnostics requested this timestep

    IF (sf(0,4)) THEN
      CALL diagnostics_pc2checks(                                       &
                       row_length, rows                                 &
,                      mype, at_extremity                               &
,                      T_inc_diag,q_inc_diag,qcl_inc_diag               &
,                      qcf_inc_diag,cf_inc_diag,cfl_inc_diag            &
,                      cff_inc_diag                                     &
,                                                                       &
 STASHwork4                                                             &
       )

    END IF   ! on error_code and sf(0,4)

  END SELECT ! model_type

  !
  ! 5. Deallocate fields
  !
  DEALLOCATE( T_inc_diag )
  DEALLOCATE( q_inc_diag )
  DEALLOCATE( qcl_inc_diag )
  DEALLOCATE( qcf_inc_diag )
  DEALLOCATE(  cf_inc_diag )
  DEALLOCATE( cfl_inc_diag )
  DEALLOCATE( cff_inc_diag )
  !
  ! 6. Update fields to produce the net increments which are written to
  !    the _star variables.
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k)                   &
!$OMP SHARED(model_levels,rows,row_length,theta_inc,t_n,           &
!$OMP    exner_theta_levels,q_inc,q_n,qcl_inc,qcl_n,qcf_inc,qcf_n, &
!$OMP    cf_inc,cf_n,cfl_inc,cfl_n,cff_inc,cff_n,theta_star,       &
!$OMP    q_star,qcl_star,qcf_star,cf_star,cfl_star,cff_star)       &
!$OMP SCHEDULE(STATIC)   
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        theta_star(i,j,k) = theta_inc(i,j,k)                            &
                           /exner_theta_levels(i,j,k) - t_n(i,j,k)
        ! theta_star now holds the potential temperature increment
        q_star(i,j,k)     = q_inc(i,j,k)     - q_n(i,j,k)
        qcl_star(i,j,k)   = qcl_inc(i,j,k)   - qcl_n(i,j,k)
        qcf_star(i,j,k)   = qcf_inc(i,j,k)   - qcf_n(i,j,k)
        cf_star(i,j,k)    = cf_inc(i,j,k)    - cf_n(i,j,k)
        cfl_star(i,j,k)   = cfl_inc(i,j,k)   - cfl_n(i,j,k)
        cff_star(i,j,k)   = cff_inc(i,j,k)   - cff_n(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

ELSE  ! i_cld_pc2

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k)                   &
!$OMP SHARED(model_levels,rows,row_length,theta_inc,               &
!$OMP    q_inc,qcl_inc,qcf_inc,cf_inc,cfl_inc,cff_inc,theta_star,  &
!$OMP    q_star,qcl_star,qcf_star,cf_star,cfl_star,cff_star)       &
!$OMP SCHEDULE(STATIC)   
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        theta_star(i,j,k) = theta_inc(i,j,k)
        q_star  (i,j,k) = q_inc  (i,j,k)
        qcl_star(i,j,k) = qcl_inc(i,j,k)
        qcf_star(i,j,k) = qcf_inc(i,j,k)
        cf_star (i,j,k) = cf_inc (i,j,k)
        cfl_star(i,j,k) = cfl_inc(i,j,k)
        cff_star(i,j,k) = cff_inc(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

END IF  ! i_cld_pc2

! Get inc from slow physics for physics_tendencies_mod
IF (l_retain_ph1_tendencies) THEN
  DO k = 1,tdims%k_end
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
        dtheta_ph1(i,j,k) = theta_star(i,j,k)
      END DO
    END DO
  END DO
END IF

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(k)                            &
!$OMP SHARED(l_mcr_qcf2,qcf2_star,qcf2_inc,qrain_star,qrain_inc,   &
!$OMP    l_mcr_qrain,l_mcr_qgraup,qgraup_star,qgraup_inc,          &
!$OMP    model_levels,rows,row_length)

! prognostic second cloud ice in use
IF (L_mcr_qcf2) THEN
!$OMP DO SCHEDULE(STATIC)   
  DO k = 1, model_levels
    qcf2_star(1:row_length, 1:rows, k) =                      &
                                             qcf2_inc(:,:,k)
  END DO
!$OMP END DO NOWAIT
END IF

! prognostic rain in use
IF (L_mcr_qrain) THEN
!$OMP  DO SCHEDULE(STATIC)   
  DO k = 1, model_levels
    qrain_star(1:row_length, 1:rows, k) =                     &
                                             qrain_inc(:,:,k)
  END DO
!$OMP END DO NOWAIT
END IF

! prognostic graupel in use
IF (L_mcr_qgraup) THEN
!$OMP  DO SCHEDULE(STATIC)   
  DO k = 1, model_levels
    qgraup_star(1:row_length, 1:rows, k) =                    &
                                            qgraup_inc(:,:,k)
  END DO
!$OMP END DO NOWAIT
END IF

!$OMP END PARALLEL


IF (l_casim) CALL casim_prognostics_update

! Deallocate additional microphysics variables
DEALLOCATE ( qcf2_inc )
DEALLOCATE ( qrain_inc )
DEALLOCATE ( qgraup_inc )
DEALLOCATE ( qcf2_n )
DEALLOCATE ( qrain_n )
DEALLOCATE ( qgraup_n )

! Deallocate arrays required for moisture conservation checking
IF (l_check_moist_inc .AND. model_type == mt_global .AND.               &
     (l_rain .OR. l_use_methox)) THEN
  DEALLOCATE (ep1)
  DEALLOCATE (ep2)
  DEALLOCATE (ep3)
  DEALLOCATE (ep4)
  DEALLOCATE (rho_con)
  DEALLOCATE (dqt)
END IF

IF ( L_print_L2norms ) THEN
  WRITE(umMessage, '(A)') '** L2 norms of increments after atmos_physics1 **'
  CALL umPrint(umMessage,src='ATMOS_PHYSICS1')
  WRITE(umMessage,'(A,L1)') ' Mixing ratio physics, l_mr_physics = '    &
                                            , l_mr_physics
  CALL umPrint(umMessage,src='ATMOS_PHYSICS1')
  CALL atmos_phys1_norm(                                                &
                        norm_lev_start, norm_lev_end, n_rims_to_do,     &
                        theta_star, q_star, qcl_star, qcf_star,         &
                        qrain_star, qcf2_star, qgraup_star,             &
                        cf_star, cfl_star, cff_star, u_inc, v_inc,      &
                        ls_rain, ls_snow,  dolr,                        &
                        .TRUE., .FALSE., .FALSE. )
  IF ( L_tracer ) THEN
    WRITE(umMessage, '(A)')                                             &
                   ' *** L2 norms of tracers after atmos_physics1 ***'
    CALL umPrint(umMessage,src='ATMOS_PHYSICS1')
    CALL sl_tracer1_norm(                                               &
                         super_array_size,                              &
                         norm_lev_start, norm_lev_end, n_rims_to_do,    &
                         co2, aerosol, dust_div1, dust_div2,            &
                         dust_div3, dust_div4, dust_div5, dust_div6,    &
                         soot_new, soot_agd, soot_cld,                  &
                         bmass_new, bmass_agd, bmass_cld,               &
                         ocff_new, ocff_agd, ocff_cld,                  &
                         so2, so4_aitken, so4_accu, so4_diss,           &
                         nh3, nitr_acc, nitr_diss,                      &
                         tracer_ukca, tr_ukca, ozone,                   &
                         .TRUE., .FALSE., .FALSE. )
  END IF  !  L_tracer
END IF !  L_print_L2norms

! end of routine Atmos_physics1
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE Atmos_Physics1

