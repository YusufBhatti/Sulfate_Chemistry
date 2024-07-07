! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Interface to Atmos Physics parametrizations after S-L advection.
!
! Subroutine Interface:
SUBROUTINE Atmos_Physics2(                                              &
! Parallel variables
        global_row_length, global_rows, n_proc, n_procx, n_procy        &
      , g_rows, g_row_length, NumCycles, CycleNo                        &

! model dimensions.
      , row_length, rows, n_rows, land_points                           &
      , bl_levels, st_levels, sm_levels, cloud_levels                   &
      , land_ice_points, soil_points, n_cca_levels, ntiles, tr_levels   &
      , first_constant_r_rho_level, dim_cs1, dim_cs2                    &

! Model switches
      , L_dry, L_lbc_old                                                &
      , l_cal360                                                        &

! Model Parameters
      , rhcrit, co2_mmr, tr_vars, tr_ukca                               &

! in coordinate information
      , unscaled_dry_rho                                                &
      , delta_lambda, delta_phi                                         &
! This following line is not used by Endgame.
      , dlambda_p, dphi_p, wt_lambda_p, wt_lambda_u, wt_phi_p, wt_phi_v &
! The following line is used by Endgame.
      , lat_rot_NP, long_rot_NP , f3_at_u                               &

! in time stepping information.
      , val_year, val_day_number, val_hour, val_minute                  &
      , val_second,                                                     &

! River routing
       aocpl_row_length,aocpl_p_rows,xpa,xua,xva,ypa,yua,yva,           &
       g_p_field, g_r_field, a_steps_since_riv, river_row_length,       &
       river_rows, global_river_row_length, global_river_rows,          &
       river_vel, river_mcoef, i_river_vn,                              &
! Add inland basin outflow to arguments
       trivdir, trivseq, twatstor,inlandout_atm,                        &

!  Add lake evaporation:
        acc_lake_evap,                                                  &

! Grid-to-grid river routing
        r_area, slope, flowobs1, r_inext, r_jnext, r_land,              &
        substore, surfstore, flowin, bflowin,                           &
!
! diagnostic info
       STASHwork3,STASHwork5,STASHwork8,STASHwork9,STASHwork19,         &
       STASHwork26                                                      &

! SCM Diagnostics (dummy values in full UM)
      , nscmdpkgs, l_scmdiags, conv_mode, l_emcorr_opt                  &

! in data fields.
      , theta, q, qcl, qcf, qrain, qgraup, qcf2                         &
      , rho_wet_rsq, u, v, w, etadot, p, p_star                         &
      , exner_rho_levels, exner_theta_levels                            &
      , land_sea_mask, p_theta_levels                                   &
! Mixing ratio prognostics; passed in separately to the q fields so
! that any convection scheme which works with mixing ratios
! can just use them instead of having to convert the q-fields back to
! mixing ratios, which would seem perverse!
      , m_v, m_cl, m_cf, m_cf2, m_r, m_gr                               &

! ancillary fields and fields needed to be kept from timestep to
! timestep

      , land_index, land_ice_index, soil_index, canopy_gb, snow_depth   &
      , hcon, hcap, smvccl, smvcwt, smvcst, z0m_soil, sthf, sthu        &
      , sil_orog_land, ho2r2_orog, sd_orog, di, ice_fract               &
      , u_0, v_0, u_0_p, v_0_p                                          &
      , cca0_dp, cca0_md, cca0_sh                                       &
      , cca0, ccb0, cct0, cclwp0, ccw_out, lcbase_out, t_soil           &
      , ti, ti_gb, k_sice_ml                                            &
      , tstar, z0msea, ice_fract_ncat,di_ncat,satcon,sathh,clapp        &
      , soil_layer_moisture, t1_sd, q1_sd, zh, ddmfx                    &
      , area_cloud_fraction, bulk_cloud_fraction_halos                  &
      , cloud_fraction_liquid_halos, cloud_fraction_frozen_halos        &
      , ls_rain, ls_rainfrac, ls_snow, ls_graup, micro_tends, totalppn  &
      , photosynth_act_rad, rad_hr                                      &
      , soil_clay,soil_silt,soil_sand,dust_mrel1,dust_mrel2,dust_mrel3  &
      , dust_mrel4,dust_mrel5,dust_mrel6                                &
      , so2_hilem, so2_em, nh3_em, dms_em, soot_hilem, soot_em          &
      , ocff_hilem, ocff_em, co2_emits, co2flux                         &
      , deep_flag, past_precip, past_conv_ht                            &

! in/out
      , theta_star, q_star, qcl_star, qcf_star, qrain_star, qgraup_star &
      , qcf2_star, cf_star, cfl_star                                    &
      , cff_star, R_u, R_v, R_w, sum_eng_fluxes, sum_moist_flux         &

! In/Out tracer fields
      , aerosol, free_tracers, tracer_ukca                              &
      , dust_div1,dust_div2,dust_div3,dust_div4,dust_div5,dust_div6     &

      , so2, dms, so4_aitken, so4_accu, so4_diss, nh3, soot_new         &
      , soot_aged, soot_cld, bmass_new, bmass_aged, bmass_cld           &
      , ocff_new, ocff_aged, ocff_cld, nitr_acc, nitr_diss, co2         &

! IN/OUT RIVERS
      , tot_surf_runoff, tot_sub_runoff                                 &
!
! out fields
      , ntml, cumulus, nbdsc, ntdsc ,rhcpt                              &
      , rhc_row_length, rhc_rows                                        &
      , zlcl_mixed                                                      &
! Additional variables for MOSES II
      , frac, frac_disturb, canht_ft, lai_ft, canopy, catch, catch_snow &
      , snow_grnd, snow_tile, z0_tile, z0h_tile_bare                    &
      , tstar_tile, tsurf_elev_surft, infil_tile, rgrain                &
      , cs, gs, co2_dim_row, co2_dim_len                                &
      , asteps_since_triffid, a_step                                    &
      , g_leaf_acc, g_leaf_phen_acc, npp_ft_acc, resp_w_ft_acc          &
      , resp_s_acc, land_pts_trif, npft_trif, olr, lw_down, sw_tile     &
      , fland_ctile,tstar_land_ctile,tstar_sea_ctile                    &
      , tstar_sice_cat_ctile,tstar_sice_ctile                           &
      , albsoil,cos_zenith_angle                                        &

! INOUT variables for TKE based turbulence schemes
      , e_trb, tsq_trb, qsq_trb, cov_trb, zhpar_shcu                    &
! OUT variable to store BL w-variance diagnostic
      , bl_w_var                                                        &
! IN/OUT convection prognostics
      , conv_prog_1, conv_prog_2, conv_prog_3, conv_prog_precip         &
      , ux_ccp, uy_ccp, um_ccp, g_ccp, h_ccp, riso_ccp, rdir_ccp        &
! Additional variables required for large-scale hydrology:
      , fexp,gamtot,ti_mean,ti_sig,fsat,fwetl,zw                        &
      , sthzw,a_fsat,c_fsat,a_fwet,c_fwet                               &
! JULES 2 Prognostics
      , snowdepth_p, rho_snow_grnd_p                                    &
      , nsnow_p                                                         &
      , ds_p, sice_p, sliq_p, tsnowlayer_p, rho_snow_p, rgrainl_p       &
! FLake lake scheme prognostics
      , lake_depth_p, lake_fetch_p, lake_t_mean_p, lake_t_mxl_p         &
      , lake_t_ice_p, lake_h_mxl_p, lake_h_ice_p,  lake_shape_p         &
      , lake_g_dt_p                                                     &
! Cariolle ozone and associated ancillaries
      , ozone_tracer                                                    &
!
! Additional screen-level variables
      , TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans                          &

! Variables required for COSP (out)
      , cosp_crain_3d, cosp_csnow_3d                                    &

! Variables for SCM and idealised UM
      , flux_e,flux_h,ustar_in,L_spec_z0,z0m_scm,z0h_scm,               &
!
! error information
        Error_code  )


USE swap_bounds_mv_mod, ONLY: swap_bounds_mv
USE swap_bounds_2d_mv_mod, ONLY: swap_bounds_2d_mv
USE dynamics_input_mod, ONLY: l_check_moist_inc

USE atmos_physics2_alloc_mod

USE nlsizes_namelist_mod, ONLY: model_levels, super_array_size

USE gen_phys_inputs_mod, ONLY: l_mr_physics

USE bl_option_mod, ONLY: i_bl_vn, i_bl_vn_1a, i_bl_vn_0,                &
                         on, l_quick_ap2, off, l_conv_tke, max_tke,     &
                         l_calc_tau_at_p
USE bl_diags_mod, ONLY: BL_diag
USE sf_diags_mod, ONLY: sf_diag
USE bl_pert_theta_mod, ONLY: bl_pert_theta
USE eg_idl_friction_mod, ONLY: eg_idl_friction
USE atm_fields_bounds_mod, ONLY: pdims, pdims_s,                        &
                                 tdims, tdims_s, tdims_l,               &
                                 udims, udims_s, vdims, vdims_s,        &
                                 wdims, wdims_l, wdims_s

USE vertnamelist_mod, ONLY: z_top_of_model

USE planet_constants_mod, ONLY: kappa, p_zero, g

USE water_constants_mod, ONLY: lc, lf, tfs, rho_water, tm

USE timestep_mod, ONLY: timestep_number, timestep, recip_timestep

USE level_heights_mod, ONLY:                                            &
                r_theta_levels, r_rho_levels,                           &
                eta_theta_levels, eta_rho_levels

USE rad_input_mod, ONLY: a_sw_radstep_diag, a_lw_radstep_diag,          &
                         a_sw_radstep_prog, a_lw_radstep_prog,          &
                         l_use_cariolle, l_cca_dp_prog,                 &
                         l_cca_md_prog, l_cca_sh_prog

USE trignometric_mod, ONLY: sec_theta_latitude,                         &
                      sin_theta_longitude, cos_theta_longitude,         &
                      cos_theta_latitude,                               &
                      FV_cos_theta_latitude

USE g_wave_input_mod, ONLY:                                             &
    l_use_ussp

USE cv_run_mod,  ONLY:                                                  &
    i_convection_vn, i_convection_vn_5a, i_convection_vn_6a,            &
    cape_bottom, cape_top, l_ccrad,                                     &
    rad_cloud_decay_opt, l_mom, l_rediagnosis, l_pc2_diag_sh,           &
    l_param_conv, cldbase_opt_dp, cldbase_opt_md, i_cv_llcs,            &
    l_jules_flux, cnv_cold_pools, ccp_off

USE cv_param_mod,  ONLY:                                                &
    rad_decay_off, rad_decay_full_timestep, rad_decay_conv_substep

USE cv_stash_flg_mod,  ONLY:                                            &
  ! subroutine
   set_convection_output_flags,                                         &
  ! variables
    flg_up_flx, flg_dwn_flx,                                            &
    l_apply_diag

USE swapable_field_mod, ONLY:                                           &
    swapable_field_pointer_type

USE nlstcall_mod, ONLY: ltimer

USE jules_surface_mod, ONLY: l_aggregate, l_flake_model, l_urban2t,     &
                             IScrnTDiag, formdrag, explicit_stress,     &
                             ISrfExCnvGust, IP_SrfExWithCnv

USE jules_vegetation_mod, ONLY: l_triffid, l_landuse, l_trif_crop,      &
                                l_nitrogen

USE jules_sea_seaice_mod, ONLY: nice, nice_use, l_use_dtstar_sea,       &
                           l_sice_multilayers, l_ctile, buddy_sea

USE jules_snow_mod, ONLY: nsmax, rho_snow_const

USE p_s_parms, ONLY:   smvccl_levs => smvccl_soilt,                     &
                       smvcwt_levs => smvcwt_soilt,                     &
                       smvcst_levs => smvcst_soilt

USE surf_couple_extra_mod, ONLY: surf_couple_extra

! FLake model
USE lake_mod, ONLY:  h_snow_min_flk                                     &
                   , u_s_lake_gb                                        &
                   , surf_ht_flux_lake_ij                               &
                   , surf_ht_flux_lk_gb                                 &
                   , sw_down_gb                                         &
                   , lake_depth_gb                                      &
                   , lake_fetch_gb                                      &
                   , coriolis_param_gb                                  &
                   , lake_albedo_gb                                     &
                   , lake_t_snow_gb                                     &
                   , lake_t_ice_gb                                      &
                   , lake_t_mean_gb                                     &
                   , lake_t_mxl_gb                                      &
                   , lake_shape_factor_gb                               &
                   , lake_h_snow_gb                                     &
                   , lake_h_ice_gb                                      &
                   , lake_h_mxl_gb                                      &
                   , lake_t_sfc_gb                                      &
                   , ts1_lake_gb                                        &
                   , g_dt_gb                                            &
                   , h_snow_sw_att                                      &
                   , trap_frozen                                        &
                   , trap_unfrozen

USE ancil_info, ONLY:                                                   &
 ssi_pts                                                                &
,sea_pts                                                                &
! arrays which JULES needs to allocate
   ,ssi_index                                                           &
   ,sea_index                                                           &
   ,sice_pts_ncat                                                       &
   ,sice_index_ncat                                                     &
   ,fssi_ij                                                             &
   ,sea_frac                                                            &
   ,sice_frac_ncat

USE jules_surface_types_mod

USE missing_data_mod, ONLY: rmdi
! arrays which JULES needs to allocate
USE prognostics, ONLY:                                                  &
  snowdepth_surft                                                       &
 ,rho_snow_grnd_surft                                                   &
 ,nsnow_surft                                                           &
 ,sice_surft                                                            &
 ,sliq_surft                                                            &
 ,tsnow_surft                                                           &
 ,rho_snow_surft                                                        &
 ,rgrainl_surft                                                         &
 ,ds_surft                                                              &
 ,wood_prod_fast_gb                                                     &
 ,wood_prod_med_gb                                                      &
 ,wood_prod_slow_gb                                                     &
 ,frac_agr_prev_gb                                                      &
 ,frac_past_prev_gb                                                     &
 ,ns_pool_gb                                                            &
 ,n_inorg_soilt_lyrs                                                    &
 ,n_inorg_gb                                                            &
 ,triffid_co2_gb

USE trif_vars_mod, ONLY:                                                &
  frac_past_gb

! Import UM versions of prognostics from D1 array for passing to
! equivalent JULES variables
USE atm_fields_real_mod, ONLY:                                    &
  wood_prod_fast_d1                                               &
 ,wood_prod_med_d1                                                &
 ,wood_prod_slow_d1                                               &
 ,disturb_veg_prev                                                &
 ,pasture_frac_d1                                                 &
 ,pasture_frac_prev_d1                                            &
 ,agr_crop_frac_d1                                                &
 ,agr_crop_frac_prev_d1                                           &
 ,soil_nitro1                                                     &
 ,soil_nitro2                                                     &
 ,soil_nitro3                                                     &
 ,soil_nitro4                                                     &
 ,soil_inorgnit                                                   &
 ,nitrogen_deposition_d1                                          &
 ,hgt, hwr, wrr, disp, ztm, albwl, albrd, emisw, emisr            &
 ,triffid_co2_d1

USE trif_vars_mod, ONLY: deposition_n_gb

! Convective diagnostic output arrays
USE cv_diagnostic_array_mod, ONLY:                                      &
  cape_out, up_flux, dwn_flux                                           &
 ,theta_incr_diag_conv, conv_rain_3d, conv_snow_3d

USE c_kappai, ONLY: kappai,de

USE ni_conv_ctl_mod, ONLY: ni_conv_ctl
USE conv_diag_5a_mod, ONLY: conv_diag_5a
USE conv_diag_6a_mod, ONLY: conv_diag_6a
USE other_conv_ctl_mod, ONLY: other_conv_ctl
USE conv_cold_pools_mod, ONLY: conv_cold_pools

USE check_dmoist_inc_mod, ONLY: check_dmoist_inc

! Copy of arrays from Convection required by SKEB2
! Dev Note: request owners of APP code to store orig arrays in modules
USE stochastic_physics_run_mod,  ONLY:                                  &
    skeb2_up_flux, skeb2_dwn_flux, skeb2_cape,                          &
    i_pert_theta, pert_theta_mag, pert_theta_and_moist

! Copy of Tendencies for SPT and other schemes
USE physics_tendencies_mod,  ONLY:                                      &
    init_convection_tendencies, du_conv, dv_conv,                       &
    l_retain_conv_all_tendencies, l_retain_conv_mom_tendencies

USE dust_parameters_mod, ONLY: ndiv, ndivh, l_dust, l_dust_diag

USE carbon_options_mod, ONLY: l_co2_interactive, l_co2_emits

USE run_aerosol_mod, ONLY: l_sulpc_so2, l_sulpc_dms, l_sulpc_nh3,       &
                           l_soot, l_ocff, l_biomass, l_nitrate

USE ukca_option_mod, ONLY: l_ukca
USE cloud_inputs_mod, ONLY: i_cld_vn, l_pc2_reset
USE pc2_constants_mod, ONLY: i_cld_pc2
USE river_inputs_mod, ONLY: l_rivers
USE eng_corr_inputs_mod, ONLY: l_emcorr
USE murk_inputs_mod, ONLY: l_murk, l_murk_advect, l_murk_source
USE free_tracers_inputs_mod, ONLY: l_bl_tracer_mix
USE cosp_input_mod, ONLY: l_cosp
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE UM_ParParams, ONLY: pnorth, psouth
USE Field_Types
USE ereport_mod, ONLY: ereport

USE horiz_grid_mod, ONLY: intw_rho2w,intw_u2p,intw_v2p
USE metric_terms_mod

USE compute_chunk_size_mod, ONLY: compute_chunk_size

USE p_to_t_mod,      ONLY: p_to_t
USE p_to_u_mod,      ONLY: p_to_u
USE p_to_u_land_mod, ONLY: p_to_u_land
USE p_to_u_sea_mod,  ONLY: p_to_u_sea
USE p_to_v_mod,      ONLY: p_to_v
USE p_to_v_land_mod, ONLY: p_to_v_land
USE p_to_v_sea_mod,  ONLY: p_to_v_sea
USE u_to_p_mod,      ONLY: u_to_p
USE v_to_p_mod,      ONLY: v_to_p

USE stash_array_mod, ONLY: sf

USE s_scmop_mod,  ONLY: default_streams                                 &
  , t_inst, t_avg, t_mult, d_sl, d_bl, d_all, d_tile, d_land            &
  , scmdiag_gen, scmdiag_bl, scmdiag_surf, scmdiag_land, scmdiag_lsp    &
  , scmdiag_conv, scmdiag_incs, scmdiag_mlsnow
USE scmoutput_mod, ONLY: scmoutput
USE scm_diag_conv_mod, ONLY: scm_diag_conv

USE model_domain_mod, ONLY: model_type, mt_global, mt_lam, mt_single_column
USE urban_param, ONLY:                                                     &
    hgt_gb, hwr_gb, wrr_gb, disp_gb, ztm_gb, albwl_gb, albrd_gb,           &
    emisw_gb, emisr_gb
USE turb_diff_mod,        ONLY: l_leonard_term
USE leonard_incs_mod,     ONLY: u_inc_leonard, v_inc_leonard, w_inc_leonard, &
                                thetal_inc_leonard, qw_inc_leonard
USE leonard_term_ctl_mod, ONLY: leonard_term_ctl

USE flux_diag_mod, ONLY: flux_diag
USE cv_alloc_diag_array_mod, ONLY: cv_alloc_diag_array
USE cv_dealloc_diag_array_mod, ONLY: cv_dealloc_diag_array
USE diagnostics_conv_mod, ONLY: diagnostics_conv
USE bdy_expl3_mod, ONLY: bdy_expl3
USE ni_bl_ctl_mod, ONLY: ni_bl_ctl
USE ni_imp_ctl_mod, ONLY: ni_imp_ctl

USE umPrintMgr, ONLY: umPrint, umMessage, printstatus, PrStatus_Oper
USE errormessagelength_mod, ONLY: errormessagelength

USE atm_step_local, ONLY: L_print_L2norms, L_tracer
USE atmos_phys2_norm_mod
USE sl_tracer_norm_mod
USE turb_diff_mod, ONLY: norm_lev_start, norm_lev_end
USE lam_config_inputs_mod, ONLY: n_rims_to_do

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE
!
! Description: This version interfaces to physics schemes in sequence:
!    convection                     (optional)
!    boundary layer
!    convection                     (optional)
!    hydrology
!    river routing                  (optional)
!
!          CALLed after Semi-Lagrangian in atmosphere timestep.
! Method:
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.3 programming standards.
!
! System component: Control Atmos
!
! Declarations:
!
! arguments with intent in. ie: input variables.

! Parallel setup variables
INTEGER ::                                                              &
             ! Size of small halo in j.
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
, g_row_length (0:n_proc-1)                                             &
, NumCycles                                                             &
             ! Number of cycles (iterations) for iterative SISL
, CycleNo
             ! Sweep number

! Parameters

! Model dimensions
INTEGER ::                                                              &
  row_length                                                            &
, rows                                                                  &
, n_rows                                                                &
, land_points                                                           &
              ! IN No.of land points being processed, can be 0.
, bl_levels                                                             &
, st_levels                                                             &
              ! number of deep soil temperature levels
, sm_levels                                                             &
              ! number of deep soil moisture levels
, cloud_levels                                                          &
, land_ice_points                                                       &
                  ! number of land ice points
, soil_points                                                           &
                  ! number of soil points
, n_cca_levels                                                          &
                ! Number of levels for conv cloud
                ! amount: 1 for 2D, nlevs for 3D.
, ntiles                                                                &
                ! No. of land-surface tiles ( MOSES II )
, tr_levels                                                             &
                ! No. of free tracer levels
, first_constant_r_rho_level                                            &
                             ! 1st rho level on which r constant
, dim_cs1, dim_cs2  ! soil carbon dimensions

! Model switches

LOGICAL ::                                                              &
  L_lbc_old                                                             &
               !  false for new lbc treatment
, L_dry       ! true if model to be run with no moisture


LOGICAL ::                                                              &
  l_cal360                                                              &
                  ! true if using 360 day calender
, l_combcld_cca0
                    ! Combined cld diag uses Sec 0 convective cld
!

! model parameters
REAL ::                                                                 &
  rhcrit(model_levels)                                                  &
                      ! IN Critical relative humidity.
                            ! the values need to be tuned
                            ! for the given set of levels.
, co2_mmr
                  ! set equal to co2_start

INTEGER ::                                                              &
  tr_vars                                                               &
                    ! number of free tracer variables
, tr_ukca
                    ! number of ukca tracer variables

! RIVER routing

INTEGER ::                                                              &
 aocpl_p_rows, aocpl_row_length                                         &
, g_p_field                                                             &
                            ! IN size of global ATMOS field
, g_r_field                                                             &
                            ! IN Size of global river field
, river_row_length                                                      &
                            ! IN local river row length
, river_rows                                                            &
                            ! IN local river rows
, global_river_row_length                                               &
                            ! IN global river row length
, global_river_rows                                                     &
                                  ! IN GLOBAL river rows
, a_steps_since_riv         ! IN No. Physics timsteps since last
!                                 ! call to river routing

! Local parameters:
INTEGER ::                                                              &
       swap_levels                                                      &
                            ! no. levels for SWAPBOUNDS
, gather_pe_trip                                                        &
                            ! pe River routing to be run on
, info                                                                  &
                        ! Return code from MPP
, icode                 ! Return code : 0 Normal Exit : >0 Error

CHARACTER (LEN=errormessagelength) :: cmessage      ! used for ereport

PARAMETER(swap_levels=1)              ! by definition for A- T

! data to regrid runoff from ATMOS to river routing grid
REAL ::                                                                 &
 xpa(aocpl_row_length+1)                                                &
                         ! IN Atmosphere TP long coordinates
,xua(0:aocpl_row_length)                                                &
                         ! IN Atmosphere U long coordinates
,xva(aocpl_row_length+1)                                                &
                         ! IN Atmosphere V long coordinates
,ypa(aocpl_p_rows)                                                      &
                         ! IN Atmosphere TP lat coordinates
,yua(aocpl_p_rows)                                                      &
                         ! IN Atmosphere U lat coordinates
,yva(0:aocpl_p_rows)     ! IN Atmosphere V lat coordinates
! Data to run river routing
REAL ::                                                                 &
 trivdir(river_row_length, river_rows)                                  &
                                        ! IN River direction file
,trivseq(river_row_length, river_rows)                                  &
                                        ! IN River sequence file
,twatstor(river_row_length, river_rows)                                 &
                                        ! IN/OUT Water storage
!                                             ! file (Kg)
      ,river_vel                                                        &
                                              ! IN river velocity (m/s)
      ,river_mcoef
                                              ! IN meander coefficient

INTEGER ::                                                              &
 i_river_vn
                                        ! IN river model type

!
REAL ::                                                                 &
r_area(row_length,rows),                                                &
!               accumulated areas file
       r_inext(row_length,rows),                                        &
!               x-coordinate of downstream grid point
       r_jnext(row_length,rows),                                        &
!               y-coordinate of downstream grid point
        slope(row_length,rows),                                         &
!             slopes (not used yet)
        flowobs1(row_length,rows),                                      &
!             initialisation for flows
       r_land(row_length,rows),                                         &
!             land/river/sea
       substore(row_length,rows),                                       &
!          routing sub_surface store (mm)
       surfstore(row_length,rows),                                      &
!          routing surface store (mm)
       flowin(row_length,rows),                                         &
!          surface lateral inflow (mm)
       bflowin(row_length,rows)
!          sub-surface lateral inflow (mm)
! Diagnostics info
REAL ::                                                                 &
 stashwork3(*)                                                          &
                   ! STASH workspace for section 3 (b layer)
,stashwork5(*)                                                          &
                   ! STASH workspace for section 5 (convection)
,stashwork8(*)                                                          &
                   ! STASH workspace for section 8 (hydrology)
,STASHwork9(*)                                                          &
                   ! STASH workspace for section 9 (LS Cloud)
,STASHwork19(*)                                                         &
                   ! STASH workspace for section 19 (Veg)
,STASHwork26(*)    ! STASH workspace for sect. 26 (River routing)
!
! Data arrays
REAL :: u      (udims_s%i_start:udims_s%i_end,                          &
                udims_s%j_start:udims_s%j_end,                          &
                udims_s%k_start:udims_s%k_end)
REAL :: v      (vdims_s%i_start:vdims_s%i_end,                          &
                vdims_s%j_start:vdims_s%j_end,                          &
                vdims_s%k_start:vdims_s%k_end)
REAL :: w      (wdims_s%i_start:wdims_s%i_end,                          &
                wdims_s%j_start:wdims_s%j_end,                          &
                wdims_s%k_start:wdims_s%k_end)
REAL :: etadot (wdims_s%i_start:wdims_s%i_end,                          &
                wdims_s%j_start:wdims_s%j_end,                          &
                wdims_s%k_start:wdims_s%k_end)
REAL :: rho_wet_rsq(pdims_s%i_start:pdims_s%i_end,                      &
                    pdims_s%j_start:pdims_s%j_end,                      &
                    pdims_s%k_start:pdims_s%k_end)
REAL :: p      (pdims_s%i_start:pdims_s%i_end,                          &
                pdims_s%j_start:pdims_s%j_end,                          &
                pdims_s%k_start:pdims_s%k_end)
REAL :: p_theta_levels                                                  &
               (tdims_s%i_start:tdims_s%i_end,                          &
                tdims_s%j_start:tdims_s%j_end,                          &
                tdims_s%k_start:tdims_s%k_end)
REAL :: theta  (tdims_s%i_start:tdims_s%i_end,                          &
                tdims_s%j_start:tdims_s%j_end,                          &
                tdims_s%k_start:tdims_s%k_end)
REAL :: exner_rho_levels                                                &
               (pdims_s%i_start:pdims_s%i_end,                          &
                pdims_s%j_start:pdims_s%j_end,                          &
                pdims_s%k_start:pdims_s%k_end + 1)
REAL :: exner_theta_levels                                              &
               (tdims_s%i_start:tdims_s%i_end,                          &
                tdims_s%j_start:tdims_s%j_end,                          &
                tdims_s%k_start:tdims_s%k_end)

REAL :: p_star (pdims%i_start:pdims%i_end,                              &
                pdims%j_start:pdims%j_end)

REAL :: q      (tdims_l%i_start:tdims_l%i_end,                          &
                tdims_l%j_start:tdims_l%j_end,                          &
                tdims_l%k_start:tdims_l%k_end)
REAL :: qcl    (tdims_l%i_start:tdims_l%i_end,                          &
                tdims_l%j_start:tdims_l%j_end,                          &
                tdims_l%k_start:tdims_l%k_end)
REAL :: qcf    (tdims_l%i_start:tdims_l%i_end,                          &
                tdims_l%j_start:tdims_l%j_end,                          &
                tdims_l%k_start:tdims_l%k_end)

REAL, INTENT (IN) ::                                                    &
qrain  (tdims_l%i_start:tdims_l%i_end,     & ! prognostic rain if present
        tdims_l%j_start:tdims_l%j_end,                                  &
        tdims_l%k_start:tdims_l%k_end)                                  &
,qgraup (tdims_l%i_start:tdims_l%i_end,     & ! prognostic graupel
        tdims_l%j_start:tdims_l%j_end,     & ! if present
        tdims_l%k_start:tdims_l%k_end)                                  &
,qcf2   (tdims_l%i_start:tdims_l%i_end,     & ! 2nd ice type if present
        tdims_l%j_start:tdims_l%j_end,                                  &
        tdims_l%k_start:tdims_l%k_end)

! Start-of-timestep mixing ratios of water species:
                    ! Water vapour
REAL, INTENT(IN) :: m_v   ( tdims_s%i_start:tdims_s%i_end,              &
                            tdims_s%j_start:tdims_s%j_end,              &
                            tdims_s%k_start:tdims_s%k_end )
                    ! Cloud liquid water
REAL, INTENT(IN) :: m_cl  ( tdims_s%i_start:tdims_s%i_end,              &
                            tdims_s%j_start:tdims_s%j_end,              &
                            tdims_s%k_start:tdims_s%k_end )
                    ! Cloud ice and snow
REAL, INTENT(IN) :: m_cf  ( tdims_s%i_start:tdims_s%i_end,              &
                            tdims_s%j_start:tdims_s%j_end,              &
                            tdims_s%k_start:tdims_s%k_end )
                    ! 2nd cloud ice category (optional)
REAL, INTENT(IN) :: m_cf2 ( tdims_s%i_start:tdims_s%i_end,              &
                            tdims_s%j_start:tdims_s%j_end,              &
                            tdims_s%k_start:tdims_s%k_end )
                    ! Prognostic rain water (optional)
REAL, INTENT(IN) :: m_r   ( tdims_s%i_start:tdims_s%i_end,              &
                            tdims_s%j_start:tdims_s%j_end,              &
                            tdims_s%k_start:tdims_s%k_end )
                    ! Prognostic graupel (optional)
REAL, INTENT(IN) :: m_gr  ( tdims_s%i_start:tdims_s%i_end,              &
                            tdims_s%j_start:tdims_s%j_end,              &
                            tdims_s%k_start:tdims_s%k_end )

REAL :: unscaled_dry_rho                                                &
                    (pdims_s%i_start:pdims_s%i_end,                     &
                     pdims_s%j_start:pdims_s%j_end,                     &
                     pdims_s%k_start:pdims_s%k_end)
       ! unscaled dry density

! Convection prognostics
REAL, INTENT (INOUT) :: conv_prog_1(      tdims_s%i_start:tdims_s%i_end,&
                                          tdims_s%j_start:tdims_s%j_end,&
                                          tdims_s%k_start:tdims_s%k_end )
REAL, INTENT (INOUT) :: conv_prog_2(      tdims_s%i_start:tdims_s%i_end,&
                                          tdims_s%j_start:tdims_s%j_end,&
                                          tdims_s%k_start:tdims_s%k_end )
REAL, INTENT (INOUT) :: conv_prog_3(      tdims_s%i_start:tdims_s%i_end,&
                                          tdims_s%j_start:tdims_s%j_end,&
                                          tdims_s%k_start:tdims_s%k_end )
REAL, INTENT (INOUT) :: conv_prog_precip( tdims_s%i_start:tdims_s%i_end,&
                                          tdims_s%j_start:tdims_s%j_end,&
                                          tdims_s%k_start:tdims_s%k_end )
! convective cold-pool prognostics
REAL, INTENT (INOUT) ::     ux_ccp( pdims%i_start:pdims%i_end,          &
                                    pdims%j_start:pdims%j_end)          &
!            x-component of vector sum of c.c.p. front velocities (m/s)
!
                        ,   uy_ccp( pdims%i_start:pdims%i_end,          &
                                    pdims%j_start:pdims%j_end)          &
!            y-component of vector sum of c.c.p. front velocities (m/s)
!
                        ,   um_ccp( pdims%i_start:pdims%i_end,          &
                                    pdims%j_start:pdims%j_end)          &
!            scalar sum of c.c.p. front speeds (m/s)
!
                        ,    g_ccp( pdims%i_start:pdims%i_end,          &
                                    pdims%j_start:pdims%j_end)          &
!            gridbox c.c.p. reduced gravity (m/s^2)
!
                        ,    h_ccp( pdims%i_start:pdims%i_end,          &
                                    pdims%j_start:pdims%j_end)          &
!            gridbox c.c.p. height (m)
!
                        , riso_ccp( pdims%i_start:pdims%i_end,          &
                                    pdims%j_start:pdims%j_end)          &
!             remain counter (isotropic) (dimensionless)
!
                        , rdir_ccp( pdims%i_start:pdims%i_end,          &
                                    pdims%j_start:pdims%j_end)
!             remain counter (directed) (dimensionless)

! ancillary arrays and fields required to be saved from timestep to
! timestep.

REAL ::                                                                 &
  tstar(row_length, rows)

LOGICAL ::                                                              &
  land_sea_mask(row_length, rows)

INTEGER ::                                                              &
  land_index (land_points)                                              &
                                ! set from land_sea_mask
, land_ice_index (land_points)                                          &
                                ! Array of land ice points.
, soil_index(land_points)       ! Array of soil points.

REAL :: u_0  (udims%i_start:udims%i_end,                                &
              udims%j_start:udims%j_end)
                          ! set to zero
REAL :: v_0  (vdims%i_start:vdims%i_end,                                &
              vdims%j_start:vdims%j_end)
                          ! set to zero
REAL :: u_0_p(pdims%i_start:pdims%i_end,                                &
              pdims%j_start:pdims%j_end)
                            ! set to zero
REAL :: v_0_p(pdims%i_start:pdims%i_end,                                &
              pdims%j_start:pdims%j_end)
                          ! set to zero

REAL :: totalppn(tdims%i_start:tdims%i_end,                             &
                 tdims%j_start:tdims%j_end)

REAL ::                                                                 &
  hcon (land_points)                                                    &
                       ! soil/qrparm.soil.hcond
, hcap (land_points)                                                    &
                       ! soil/qrparm.soil.hcap
, smvccl (land_points)                                                  &
                       ! soil/qrparm.soil.crit
, smvcwt (land_points)                                                  &
                       ! soil/qrparm.soil.wilt
, smvcst (land_points)                                                  &
                       ! soil/qrparm.soil.satn
, z0m_soil(land_points)                                                 &
                       ! soil/qrparm.soil momentum z0
, sthf(land_points,sm_levels)                                           &
                                ! IN Frozen soil moisture content
!                of each layer as a fraction of saturation.
      , sthu(land_points,sm_levels)                                     &
                                      ! IN Unfrozen soil moisture
!                content of each layer as a fraction of saturation.
      , canopy_gb (land_points)                                         &
                                ! set to zero.
      , snow_depth (row_length, rows) ! snow/qrclim.snow.(month)

REAL ::                                                                 &
  ice_fract (row_length, rows)                                          &
                               ! ice/qrclim.ice.(month)
, di(row_length, rows)                                                  &
                        ! ice/qrclim.ice_thick.(month)
, ice_fract_ncat(row_length, rows, nice)                                &
, di_ncat(row_length, rows, nice)                                       &
, z0msea(row_length, rows)                                              &
                             ! Sea surface roughness
, z0m_scm(row_length, rows)                                             &
                             ! Fixed sea surface roughness
                             ! length(m) for MOMENTUM,
                             ! used in SCM
, z0h_scm(row_length, rows)  ! Fixed sea surface roughness
                             ! length(m) for HEAT,
                             ! used in SCM

REAL ::                                                                 &
  sil_orog_land (land_points)                                           &
                             ! orog/qrparm.orog.as
, ho2r2_orog (land_points)                                              &
                             ! orog/qrparm.orog.h2root2
, sd_orog (land_points)                                                 &
                             ! orog/qrparm.orog.stdev
, t_soil(land_points,sm_levels)                                         &
                             ! slt/qrclim.slt_pm(lev).(month)
, ti_gb(row_length, rows)                                               &
                             ! sea ice sfc layer temp (ice mean)
, ti(row_length, rows, nice)                                            &
                             ! sea ice sfc layer temp on categories
, k_sice_ml(row_length, rows, nice)                                     &
                             ! sea ice effective conductivity in
                             !  sfc layer on categories (W/m2/K)
                             !  (only set if l_sice_multilayers=T)
, t1_sd(row_length, rows)                                               &
                          ! set to zero initially
, q1_sd(row_length, rows)                                               &
                          ! set to zero initially
, zh(row_length, rows)                                                  &
                       ! boundary layer height
, zh_prev(row_length, rows)                                             &
                            ! boundary layer height from previous
!                                 ! timestep
      , ddmfx(row_length, rows)
!                                 ! Convective downdraught mass-flux
!                                 ! at cloud-base

REAL ::                                                                 &
  clapp(land_points)                                                    &
                       !  qrparm.soil.bwag ?
!                               Clapp-Hornberger exponent.
      , satcon(land_points)                                             &
                                 !  qrparm.soil.satcon
      , sathh(land_points)                                              &
                             !  soil water suction
      , soil_layer_moisture(land_points, sm_levels)
                             !  qrclim.smc_pm(lev).(month)

! CLoud fields

! local variables.
INTEGER ::                                                              &
  rhc_row_length                                                        &
                  ! Row length for RHcrit array
, rhc_rows        ! Row number for RHcrit array

REAL :: area_cloud_fraction (tdims%i_start:tdims%i_end,                 &
                             tdims%j_start:tdims%j_end,                 &
                                         1:tdims%k_end)
REAL :: bulk_cloud_fraction_halos                                       &
                            (tdims_l%i_start:tdims_l%i_end,             &
                             tdims_l%j_start:tdims_l%j_end,             &
                             tdims_l%k_start:tdims_l%k_end)
REAL :: cloud_fraction_liquid_halos                                     &
                            (tdims_l%i_start:tdims_l%i_end,             &
                             tdims_l%j_start:tdims_l%j_end,             &
                             tdims_l%k_start:tdims_l%k_end)
REAL :: cloud_fraction_frozen_halos                                     &
                            (tdims_l%i_start:tdims_l%i_end,             &
                             tdims_l%j_start:tdims_l%j_end,             &
                             tdims_l%k_start:tdims_l%k_end)
REAL :: rhcpt (rhc_row_length, rhc_rows,1:tdims%k_end)

! Rain fields
REAL ::                                                                 &
  ls_rain(row_length, rows)                                             &
, ls_snow(row_length, rows)                                             &
, ls_graup(row_length, rows)                                            &
, micro_tends(row_length, rows, 2, bl_levels)                           &
!                          ! Tendencies from microphys within BL levels
!                          ! (TL, K/s; QW, kg/kg/s)
      , conv_rain(row_length, rows)                                     &
      , conv_snow(row_length, rows)

REAL, INTENT(INOUT) ::                                                  &
  ls_rainfrac(land_points)
                     ! Rain fraction array on land points

! height of lcl in a well-mixed BL (types 3 or 4), 0 otherwise
REAL :: zlcl_mixed(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                     ! passed from atmos_physics1

! Radiation fields
REAL ::                                                                 &
  photosynth_act_rad(row_length, rows)                                  &
                                       ! Net downward
!                                 shortwave radiation in band 1 (w/m2).
      , rad_hr(row_length, rows, 2, bl_levels)
                                               !
!                               ! BL radiative (LW,SW) heating rates

! Fields for mineral dust source flux calculations
REAL, INTENT(IN) ::                                                     &
  soil_clay ( row_length, rows )                                        &
, soil_silt ( row_length, rows )                                        &
, soil_sand ( row_length, rows )                                        &
, dust_mrel1 ( row_length, rows )                                       &
, dust_mrel2 ( row_length, rows )                                       &
, dust_mrel3 ( row_length, rows )                                       &
, dust_mrel4 ( row_length, rows )                                       &
, dust_mrel5 ( row_length, rows )                                       &
, dust_mrel6 ( row_length, rows )

! Tracer emissions
REAL, INTENT(IN) ::                                                     &
  so2_hilem ( row_length, rows )                                        &
, so2_em    ( row_length, rows )                                        &
, nh3_em    ( row_length, rows )                                        &
, dms_em    ( row_length, rows )                                        &
, soot_hilem( row_length, rows )                                        &
, soot_em   ( row_length, rows )                                        &
, ocff_hilem( row_length, rows )                                        &
, ocff_em   ( row_length, rows )

! CO2 fields
REAL ::                                                                 &
  co2_emits ( row_length, rows )                                        &
, co2flux ( row_length, rows )

! Fields holding information on past history of convection. Note if
! prognostics not requested these fields are unset and there is no
! space allocated for them at the top level of the UM.
! Fields are reduced over time or reset if convection occurs.
REAL, INTENT(INOUT) ::                                                  &
 deep_flag(row_length,rows)   &  ! Value between 0.0 and 1.0
                                 ! 1 if deep convection last timestep
                                 ! Reduced according to a decay period
                                 ! if no deep convection
,past_precip(row_length,rows) &  ! Rate of convective precip (kg/m2/s)
                                 ! decayed over time.
,past_conv_ht(row_length,rows)   ! Height of convection (m)

! Convection
REAL ::                                                                 &
  cca (row_length, rows, n_cca_levels)                                  &
, cclwp(row_length, rows) ! condensed water path (KG/M**2)

INTEGER ::                                                              &
  ccb (row_length, rows)                                                &
, cct (row_length, rows)


REAL ::                                                                 &
     ! local vertical co-ordinate information
  delta_lambda                                                          &
, delta_phi

REAL ::                                                                 &
     !VarRes horizontal co-ordinate information
  dlambda_p(1-halo_i:row_length+halo_i)                                 &
, dphi_p(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)              &
, wt_lambda_p(1-halo_i:row_length+halo_i)                               &
, wt_lambda_u(1-halo_i:row_length+halo_i)                               &
, wt_phi_p(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)            &
, wt_phi_v(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)


REAL,INTENT(IN) ::                                                      &
 f3_at_u(1-offx:row_length+offx, 1-offy:rows+offy)

! time information for current timestep
INTEGER ::                                                              &
  val_year                                                              &
, val_day_number                                                        &
, val_hour                                                              &
, val_minute                                                            &
, val_second

! Diagnostic variables
REAL ::                                                                 &
  lat_rot_NP                                                            &
, long_rot_NP

! Variables to be used by COSP
!  Convective rainfall
REAL,INTENT(OUT) :: cosp_crain_3d(row_length,rows,model_levels)
!  Convective snowfall
REAL,INTENT(OUT) :: cosp_csnow_3d(row_length,rows,model_levels)

! arguments with intent in/out. ie: input variables changed on output.


REAL :: theta_star(tdims%i_start:tdims%i_end,                           &
                   tdims%j_start:tdims%j_end,                           &
                   tdims%k_start:tdims%k_end)
REAL :: q_star    (tdims%i_start:tdims%i_end,                           &
                   tdims%j_start:tdims%j_end,                           &
                   tdims%k_start:tdims%k_end)
REAL :: qcl_star  (tdims%i_start:tdims%i_end,                           &
                   tdims%j_start:tdims%j_end,                           &
                   tdims%k_start:tdims%k_end)
REAL :: qcf_star  (tdims%i_start:tdims%i_end,                           &
                   tdims%j_start:tdims%j_end,                           &
                   tdims%k_start:tdims%k_end)

! Extra prognostics -
! Passed down so convection can do conversions to/from mixing to specific
! accurately if any are present in the run.
! Also new convection schemes in other_conv_ctl may modify them.
REAL, INTENT (INOUT) ::                                                 &
qrain_star(tdims%i_start:tdims%i_end,     & ! prognostic rain
           tdims%j_start:tdims%j_end,     & ! if present
           tdims%k_start:tdims%k_end)                                   &
,qgraup_star(tdims%i_start:tdims%i_end,   & ! prognostic graupel
           tdims%j_start:tdims%j_end,     & ! if present
           tdims%k_start:tdims%k_end)                                   &
,qcf2_star (tdims%i_start:tdims%i_end,    & ! 2nd ice type if present
           tdims%j_start:tdims%j_end,                                   &
           tdims%k_start:tdims%k_end)

REAL :: cf_star   (tdims%i_start:tdims%i_end,                           &
                   tdims%j_start:tdims%j_end,                           &
                   tdims%k_start:tdims%k_end)
REAL :: cfl_star  (tdims%i_start:tdims%i_end,                           &
                   tdims%j_start:tdims%j_end,                           &
                   tdims%k_start:tdims%k_end)
REAL :: cff_star  (tdims%i_start:tdims%i_end,                           &
                   tdims%j_start:tdims%j_end,                           &
                   tdims%k_start:tdims%k_end)
REAL, TARGET :: R_u(udims_s%i_start:udims_s%i_end,                      &
                    udims_s%j_start:udims_s%j_end,                      &
                    udims_s%k_start:udims_s%k_end)
REAL, TARGET :: R_v(vdims_s%i_start:vdims_s%i_end,                      &
                    vdims_s%j_start:vdims_s%j_end,                      &
                    vdims_s%k_start:vdims_s%k_end)
REAL :: R_w       (wdims%i_start:wdims%i_end,                           &
                   wdims%j_start:wdims%j_end,                           &
                   wdims%k_start:wdims%k_end)

REAL :: sum_eng_fluxes(row_length, rows)
REAL :: sum_moist_flux(row_length, rows)

! Used in SCM for prescribed surface flux forcing
REAL ::                                                                 &
  flux_e(row_length, rows)                                              &
                             ! Surface latent heat flux (W/m^2)
, flux_h(row_length, rows)                                              &
                             ! Surface sensible heat flux (W/m^2)
, ustar_in(row_length, rows) ! Surf. friction velocity (m/s)

!    IN logicals for surface forcing
LOGICAL ::                                                              &
  L_spec_z0,  &  ! T if roughness lengths have been specified
  l_aero_classic ! T if any aerosol scheme is active

! Tracer variables
REAL, INTENT(INOUT) ::                                                  &
  aerosol     ( tdims_s%i_start:tdims_s%i_end,                          &
                tdims_s%j_start:tdims_s%j_end,                          &
                tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) ::                                                  &
  free_tracers( tdims_s%i_start:tdims_s%i_end,                          &
                tdims_s%j_start:tdims_s%j_end,                          &
                tdims_s%k_start:tdims_s%k_end, tr_vars)
REAL, INTENT(INOUT) ::                                                  &
  tracer_ukca( tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end,tr_ukca)

REAL, INTENT(INOUT) ::                                                  &
  dust_div1   ( tdims_s%i_start:tdims_s%i_end,                          &
                tdims_s%j_start:tdims_s%j_end,                          &
                tdims_s%k_start:tdims_s%k_end )
REAL, INTENT(INOUT) ::                                                  &
  dust_div2   ( tdims_s%i_start:tdims_s%i_end,                          &
                tdims_s%j_start:tdims_s%j_end,                          &
                tdims_s%k_start:tdims_s%k_end )
REAL, INTENT(INOUT) ::                                                  &
  dust_div3   ( tdims_s%i_start:tdims_s%i_end,                          &
                tdims_s%j_start:tdims_s%j_end,                          &
                tdims_s%k_start:tdims_s%k_end )
REAL, INTENT(INOUT) ::                                                  &
  dust_div4   ( tdims_s%i_start:tdims_s%i_end,                          &
                tdims_s%j_start:tdims_s%j_end,                          &
                tdims_s%k_start:tdims_s%k_end )
REAL, INTENT(INOUT) ::                                                  &
  dust_div5   ( tdims_s%i_start:tdims_s%i_end,                          &
                tdims_s%j_start:tdims_s%j_end,                          &
                tdims_s%k_start:tdims_s%k_end )
REAL, INTENT(INOUT) ::                                                  &
  dust_div6   ( tdims_s%i_start:tdims_s%i_end,                          &
                tdims_s%j_start:tdims_s%j_end,                          &
                tdims_s%k_start:tdims_s%k_end )


REAL, INTENT(INOUT) ::                                                  &
  so2         ( tdims_s%i_start:tdims_s%i_end,                          &
                tdims_s%j_start:tdims_s%j_end,                          &
                tdims_s%k_start:tdims_s%k_end )
REAL, INTENT(INOUT) ::                                                  &
  dms         ( tdims_s%i_start:tdims_s%i_end,                          &
                tdims_s%j_start:tdims_s%j_end,                          &
                tdims_s%k_start:tdims_s%k_end )
REAL, INTENT(INOUT) ::                                                  &
  so4_aitken  ( tdims_s%i_start:tdims_s%i_end,                          &
                tdims_s%j_start:tdims_s%j_end,                          &
                tdims_s%k_start:tdims_s%k_end )
REAL, INTENT(INOUT) ::                                                  &
  so4_accu    ( tdims_s%i_start:tdims_s%i_end,                          &
                tdims_s%j_start:tdims_s%j_end,                          &
                tdims_s%k_start:tdims_s%k_end )
REAL, INTENT(INOUT) ::                                                  &
  so4_diss    ( tdims_s%i_start:tdims_s%i_end,                          &
                tdims_s%j_start:tdims_s%j_end,                          &
                tdims_s%k_start:tdims_s%k_end )
REAL, INTENT(INOUT) ::                                                  &
  nh3         ( tdims_s%i_start:tdims_s%i_end,                          &
                tdims_s%j_start:tdims_s%j_end,                          &
                tdims_s%k_start:tdims_s%k_end )

REAL, INTENT(INOUT) ::                                                  &
  soot_new    ( tdims_s%i_start:tdims_s%i_end,                          &
                tdims_s%j_start:tdims_s%j_end,                          &
                tdims_s%k_start:tdims_s%k_end )
REAL, INTENT(INOUT) ::                                                  &
  soot_aged   ( tdims_s%i_start:tdims_s%i_end,                          &
                tdims_s%j_start:tdims_s%j_end,                          &
                tdims_s%k_start:tdims_s%k_end )
REAL, INTENT(INOUT) ::                                                  &
  soot_cld    ( tdims_s%i_start:tdims_s%i_end,                          &
                tdims_s%j_start:tdims_s%j_end,                          &
                tdims_s%k_start:tdims_s%k_end )
REAL, INTENT(INOUT) ::                                                  &
  bmass_new   ( tdims_s%i_start:tdims_s%i_end,                          &
                tdims_s%j_start:tdims_s%j_end,                          &
                tdims_s%k_start:tdims_s%k_end )
REAL, INTENT(INOUT) ::                                                  &
  bmass_aged  ( tdims_s%i_start:tdims_s%i_end,                          &
                tdims_s%j_start:tdims_s%j_end,                          &
                tdims_s%k_start:tdims_s%k_end )
REAL, INTENT(INOUT) ::                                                  &
  bmass_cld   ( tdims_s%i_start:tdims_s%i_end,                          &
                tdims_s%j_start:tdims_s%j_end,                          &
                tdims_s%k_start:tdims_s%k_end )
REAL, INTENT(INOUT) ::                                                  &
  ocff_new    ( tdims_s%i_start:tdims_s%i_end,                          &
                tdims_s%j_start:tdims_s%j_end,                          &
                tdims_s%k_start:tdims_s%k_end )
REAL, INTENT(INOUT) ::                                                  &
  ocff_aged   ( tdims_s%i_start:tdims_s%i_end,                          &
                tdims_s%j_start:tdims_s%j_end,                          &
                tdims_s%k_start:tdims_s%k_end )
REAL, INTENT(INOUT) ::                                                  &
  ocff_cld    ( tdims_s%i_start:tdims_s%i_end,                          &
                tdims_s%j_start:tdims_s%j_end,                          &
                tdims_s%k_start:tdims_s%k_end )
REAL, INTENT(INOUT) ::                                                  &
  nitr_acc    ( tdims_s%i_start:tdims_s%i_end,                          &
                tdims_s%j_start:tdims_s%j_end,                          &
                tdims_s%k_start:tdims_s%k_end )
REAL, INTENT(INOUT) ::                                                  &
  nitr_diss   ( tdims_s%i_start:tdims_s%i_end,                          &
                tdims_s%j_start:tdims_s%j_end,                          &
                tdims_s%k_start:tdims_s%k_end )
REAL, INTENT(INOUT) ::                                                  &
  co2         ( tdims_s%i_start:tdims_s%i_end,                          &
                tdims_s%j_start:tdims_s%j_end,                          &
                tdims_s%k_start:tdims_s%k_end )

! Add definitions for the cariolle scheme.

REAL, INTENT(INOUT) ::                                                  &
   ozone_tracer(tdims_s%i_start:tdims_s%i_end,                          &
                tdims_s%j_start:tdims_s%j_end,                          &
                tdims_s%k_start:tdims_s%k_end  )

! arguments with intent out. ie: output variables.
INTEGER ::                                                              &
  ntml (row_length, rows)                                               &
, nbdsc(row_length, rows)                                               &
, ntdsc(row_length, rows)

LOGICAL ::                                                              &
  cumulus (row_length, rows) ! bl convection flag

!     Variables for screen-level diagnostics
REAL, INTENT(INOUT) :: TScrnDcl_SSI(row_length,rows)
!                           !    Decoupled screen-level temperature
!                           !    over sea or sea-ice
REAL, INTENT(INOUT) :: TScrnDcl_TILE(land_points,ntiles)
!                           !    Decoupled screen-level temperature
!                           !    over land tiles
REAL, INTENT(INOUT) :: tStbTrans(row_length,rows)
!                           !    Time since the transition

INTEGER ::                                                              &
  Error_code

INTEGER ::                                                              &
 co2_dim_len                                                            &
                             ! IN Length of a CO2 field row.
,co2_dim_row                                                            &
                             ! IN Number of CO2 field rows.
,land_pts_trif                                                          &
                             ! IN For dimensioning land fields
,npft_trif                                                              &
                             ! IN For dimensioning PFT fields
!                                  !    available only with TRIFFID.
!                                  !    Set to NPFT when TRIFFID on,
!                                  !    set to 1 when TRIFFID off.
      ,asteps_since_triffid                                             &
                                   ! IN Number of atmospheric
!                                  !    timesteps since last call
!                                  !    to TRIFFID.
      ,a_step
                                   ! IN Atmospheric timestep number.

REAL ::                                                                 &
 frac(land_points,ntype)                                                &
                                ! IN Fractions of surface types.
,frac_disturb(land_points)                                              &
                              ! IN Fraction of gridbox in which
!                                   !    vegetation is disturbed.
      ,canht_ft(land_points,npft)                                       &
                                      ! IN Canopy height (m)
      ,lai_ft(land_points,npft)                                         &
                                      ! IN Leaf area index
      ,canopy(land_points,ntiles)                                       &
                                      ! IN Surface/canopy water for
!                                  !    snow-free land tiles (kg/m2)
      ,catch(land_points,ntiles)                                        &
                                      ! IN Surface/canopy water capacity
!                                  !    of snow-free land tiles (kg/m2).
      ,catch_snow(land_points,ntiles)                                   &
                                   ! IN Coniferous canopy snow capacity
!                                  !    (kg/m2).
      ,snow_tile(land_points,ntiles)                                    &
                                      ! IN Lying snow on tiles (kg/m2)
      ,z0_tile(land_points,ntiles)                                      &
                                      ! IN Tile roughness lengths (m).
      ,z0h_tile_bare(land_points,ntiles)                                &
                                      ! IN Tile thermal roughness
                                      ! lengths without snow (m).
                                      ! c.f. z0h_tile which is after
                                      ! adjustment for snow.
      ,tstar_tile(land_points,ntiles)                                   &
                                      ! IN Surface tile temperatures
      ,tsurf_elev_surft(land_points,ntiles)                             &
                                      ! IN Temperature of elevated
                                      ! subsurface tiles (K)
      ,infil_tile(land_points,ntiles)                                   &
!                          ! IN Maximum surface infiltration
      ,rgrain(land_points,ntiles)                                       &
                                 ! INOUT Snow grain size (microns).
      ,snow_grnd(land_points,ntiles)                                    &
                                   ! INOUT Snow below canopy (kg/m2).
      ,gs(land_points)                                                  &
                                      ! INOUT "Stomatal" conductance to
!                                  !        evaporation (m/s).
      ,cs(land_points,dim_cs1)                                          &
                                   ! IN Soil carbon (kg C/m2).
      ,g_leaf_acc(land_points,npft)                                     &
                                      ! INOUT Accumulated G_LEAF
!g_leaf_acc

      ,g_leaf_phen_acc(land_points,npft)                                &
!                                  ! INOUT Accumulated leaf turnover
!                                  !       rate including phenology.
      ,npp_ft_acc(land_pts_trif,npft_trif)                              &
!                                  ! INOUT Accumulated NPP_FT
      ,resp_w_ft_acc(land_pts_trif,npft_trif)                           &
!                                  ! INOUT Accum RESP_W_FT
      ,resp_s_acc(land_pts_trif,dim_cs1)                                &
                                        ! INOUT Accumulated RESP_S
      ,olr(row_length,rows)                                             &
                                ! IN    TOA - surface upward LW on
!                               !       last radiation timestep
!                               ! OUT   Corrected TOA outward LW
      ,lw_down(row_length,rows)                                         &
                                   ! IN Surface downward LW radiation
!                                  !    (W/m2).
      ,sw_tile(land_points,ntiles)    ! IN Surface net SW radiation on
!                                  !    land tiles (W/m2).
! The following are only used if coastal tiling is switched on:
REAL ::                                                                 &
 fland_ctile(land_points)                                               &
                             ! IN Land fraction on land tiles.
,tstar_land_ctile(row_length,rows)                                      &
!                                  ! INOUT Land mean sfc temperature (K)
      ,tstar_sea_ctile(row_length,rows)                                 &
!                                  ! IN    Open sea sfc temperature (K).
      ,tstar_sice_cat_ctile(row_length,rows,nice_use)                   &
!                                  ! INOUT Category seaice sfc temp (K).
      ,tstar_sice_ctile(row_length,rows)                                &
!                                  ! INOUT Sea-ice sfc temperature (K).
      ,albsoil(land_points)                                             &
                                   ! Soil albedo.
      ,cos_zenith_angle(row_length,rows)
!                                  ! Cosine of the zenith angle

! INOUT variables for TKE based turbulence schemes
REAL ::                                                                 &
  e_trb(tdims%i_start:tdims%i_end,                                      &
        tdims%j_start:tdims%j_end, bl_levels)                           &
!                   ! TKE defined on theta levels K-1
      , tsq_trb(tdims%i_start:tdims%i_end,                              &
                tdims%j_start:tdims%j_end, bl_levels)                   &
!                   ! Self covariance of liquid potential temperature
!                   ! (thetal'**2) defined on theta levels K-1
      , qsq_trb(tdims%i_start:tdims%i_end,                              &
                tdims%j_start:tdims%j_end, bl_levels)                   &
!                   ! Self covariance of total water
!                   ! (qw'**2) defined on theta levels K-1
      , cov_trb(tdims%i_start:tdims%i_end,                              &
                tdims%j_start:tdims%j_end, bl_levels)                   &
!                   ! Correlation between thetal and qw
!                   ! (thetal'qw') defined on theta levels K-1
      , zhpar_shcu(row_length, rows)
!                   ! Height of mixed layer used to evaluate
!                   ! the non-gradient buoyancy flux

! OUT variable to store BL w-variance diagnostic
REAL, INTENT(OUT) :: bl_w_var(tdims%i_start : tdims%i_end,              &
                              tdims%j_start : tdims%j_end,              &
                                          1 : tdims%k_end)

! Additional variables required for large-scale hydrology:
REAL ::                                                                 &
 fexp(land_points)                                                      &
                      ! IN Decay factor in Sat. Conductivity
!                           !    in water table layer.
      ,gamtot(land_points)                                              &
                            ! IN Integrated complete Gamma function.
      ,ti_mean(land_points)                                             &
                            ! IN Mean topographic index.
      ,ti_sig(land_points)  ! IN Standard dev. of topographic index.
!                           !    in water table layer.
REAL ::                                                                 &
 fsat(land_points)                                                      &
                      ! INOUT Surface saturation fraction.
,fwetl(land_points)                                                     &
                      ! INOUT Wetland fraction.
,zw(land_points)                                                        &
                      ! INOUT Water table depth (m).
,sthzw(land_points)                                                     &
                      ! INOUT soil moist fract. in deep-zw layer.
,a_fsat(land_points)                                                    &
                      ! IN Fitting parameter for Fsat in LSH model
,c_fsat(land_points)                                                    &
                      ! IN Fitting parameter for Fsat in LSH model
,a_fwet(land_points)                                                    &
                      ! IN Fitting parameter for Fwet in LSH model
,c_fwet(land_points)  ! IN Fitting parameter for Fwet in LSH model


! Additional variables for SCM diagnostics which are dummy in full UM
INTEGER :: conv_mode
INTEGER :: nscmdpkgs             ! No of SCM diagnostics packages
LOGICAL :: l_scmdiags(nscmdpkgs) ! Logicals for SCM diagnostics packages
LOGICAL :: l_emcorr_opt

! loop counters
INTEGER ::                                                              &
  i, j, k, l, m, n                                                      &
, j_begin, j_end

LOGICAL ::                                                              &
  L_poles   !  include poles in etadot calc (false for LAMs)

! Diagnostic switches
! a) hydrology
LOGICAL ::                                                              &
  stf_sub_surf_roff                                                     &
, smlt

! Local parameters:
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'ATMOS_PHYSICS2'

! local variables

! Local data arrays
REAL :: T_latest             (tdims%i_start:tdims%i_end,                &
                              tdims%j_start:tdims%j_end,                &
                                          1:tdims%k_end)
REAL :: q_latest             (tdims%i_start:tdims%i_end,                &
                              tdims%j_start:tdims%j_end,                &
                                          1:tdims%k_end)
REAL :: qcl_latest           (tdims%i_start:tdims%i_end,                &
                              tdims%j_start:tdims%j_end,                &
                                          1:tdims%k_end)
REAL :: qcf_latest           (tdims%i_start:tdims%i_end,                &
                              tdims%j_start:tdims%j_end,                &
                                          1:tdims%k_end)
REAL :: cf_latest            (tdims%i_start:tdims%i_end,                &
                              tdims%j_start:tdims%j_end,                &
                                          1:tdims%k_end)
REAL :: cfl_latest           (tdims%i_start:tdims%i_end,                &
                              tdims%j_start:tdims%j_end,                &
                                          1:tdims%k_end)
REAL :: cff_latest           (tdims%i_start:tdims%i_end,                &
                              tdims%j_start:tdims%j_end,                &
                                          1:tdims%k_end)
REAL :: bulk_cloud_fraction  (tdims%i_start:tdims%i_end,                &
                              tdims%j_start:tdims%j_end,                &
                                          1:tdims%k_end)
REAL :: cloud_fraction_liquid(tdims%i_start:tdims%i_end,                &
                              tdims%j_start:tdims%j_end,                &
                                          1:tdims%k_end)
REAL :: cloud_fraction_frozen(tdims%i_start:tdims%i_end,                &
                              tdims%j_start:tdims%j_end,                &
                                          1:tdims%k_end)

REAL :: l_s_poles(row_length, first_constant_r_rho_level-1)
REAL :: l_n_poles(row_length, first_constant_r_rho_level-1)

!-------------------------------------------------------------
REAL ::                                                                 &
  p_layer_boundaries(row_length, rows, 0:model_levels)                  &
        ! pressure at layer boundaries. Same as p except at
        ! bottom level = pstar, and at top = 0.
, p_layer_centres(row_length, rows, 0:model_levels)                     &
        ! pressure at layer centres. Same as p_theta_levels
        !except bottom level = pstar, and at top = 0.
, exner_layer_boundaries(row_length, rows, 0:model_levels)              &
, exner_layer_centres(row_length, rows, 0:model_levels)

! Local arrays used by both conv_diag and convection. Calculations moved
! from conv_diag & convection to save CPU.

REAL :: z_theta(row_length, rows, model_levels)
                  ! height of theta levels above surface (m)
REAL :: z_rho  (row_length, rows, model_levels)
                  ! height of rho levels above surface (m)
REAL :: rho_wet(row_length, rows, model_levels)
                  ! wet density on rho levels (kg/m3)
REAL :: rho_wet_tq(row_length, rows, model_levels-1)
                  ! wet density on theta levels (not top) (kg/m3)
REAL :: rho_dry(row_length, rows, model_levels)
                  ! dry density on rho levels (kg/m3)
REAL :: rho_dry_theta(row_length, rows, model_levels-1)
                  ! dry density on theta levels (kg/m3)
REAL ::                                                                 &
  etadot_copy(row_length,rows,0:model_levels)
! Local arrays holding information to be passed between physics
! routines.

REAL ::                                                                 &
  ls_rain_land(land_points)                                             &
, ls_snow_land(land_points)                                             &
, ls_graup_land(land_points)                                            &
, conv_rain_land(land_points)                                           &
, conv_snow_land(land_points)                                           &
, con_rainfrac(land_points)

! Local arrays for moisture conservation checking ep1-3 could be become
! diagnostics if required.
REAL, ALLOCATABLE ::                                                    &
  ep1(:,:)               & ! method 1
 ,ep2(:,:)               & ! method 2
 ,ep3(:,:)               & ! method 3
 ,dqt(:,:,:)             & ! total moisture increment
 ,rho_dry_r2(:,:,:)        ! rho dry * r*r

INTEGER ::                                                              &
  lcbase (row_length, rows)                                             &
                               ! lowest convective cloud base leve
, lctop (row_length, rows)     ! lowest convective cloud top level

!                                        ! shallow convection
REAL ::                                                                 &
  lcca(row_length, rows)                                                &
                          ! lowest convective cloud amount (%)
, ccw(row_length, rows, model_levels)
                            ! convective cloud liquid water
                            ! (G/KG) on model levels

!=================================================================
! DO NOT REMOVE OR WRITE TO if l_ccrad=.FALSE.
!
! Intermediary versions of ccw lcbase. These are used when
! l_ccrad=.TRUE.. They exist because no space is allocated in the
! D1 array in atm_step when l_ccrad=.FALSE., without them
! overwriting will occur. They should not be removed unless an
! alternative in atm_step can be found.

REAL    :: ccw_out    (row_length, rows, model_levels)
INTEGER :: lcbase_out (row_length, rows)

! DO NOT REMOVE OR WRITE TO if l_ccrad=.FALSE.
!=================================================================

REAL :: cca0    (row_length,rows,n_cca_levels)
REAL :: cca0_dp (row_length,rows,n_cca_levels)
REAL :: cca0_md (row_length,rows,n_cca_levels)
REAL :: cca0_sh (row_length,rows,n_cca_levels)
REAL :: ccw0    (row_length,rows,model_levels)
REAL :: cca0_2d (row_length,rows)
REAL :: cclwp0  (row_length,rows)

INTEGER :: ccb0    (row_length,rows)
INTEGER :: cct0    (row_length,rows)
INTEGER :: lcbase0 (row_length, rows)

INTEGER ::                                                              &
  ntpar(row_length, rows)                                               &
                               ! top of diagnostic parcel ascent
, nlcl(row_length, rows)       ! lifting condensation level

! Convective type array ::
INTEGER :: conv_type(row_length, rows)
                ! Integer index describing convective type:
                !    0=no convection
                !    1=non-precipitating shallow
                !    2=drizzling shallow
                !    3=warm congestus
                !    ...
                !    8=deep convection

LOGICAL ::                                                              &
l_shallow(row_length, rows)     & ! Logical indicator of shallow
                                  ! shallow
, l_congestus(row_length, rows)  & ! Logical indicator of congestus
, l_congestus2(row_length, rows) & ! Logical indicator of congestus
, l_mid_level(row_length, rows)  & ! Logical indicator of mid-level convection
, l_pc2_diag_sh_pts(row_length, rows)  ! PC2 is using diagnostic
                                      ! shallow convection

REAL ::                                                                 &
 CIN_undilute(row_length, rows)                                         &
                                 ! undilute CIN from parcel ascent
                                 ! (m2/s2)
,CAPE_undilute(row_length, rows) ! undilute CAPE from parcel
                                 ! ascent (m2/s2)

REAL ::                                                                 &
 delthvu(row_length, rows)     & ! buoyancy integral
, zhpar(row_length, rows)       & ! height of ntpar
, dzh(row_length, rows)         & ! inversion thickness
, qcl_inv_top(row_length,rows)  & ! Parcel water content at inversion top
, zlcl(row_length, rows)        & ! height of nlcl
, zlcl_uv(row_length,rows)      & ! height of nlcl for uv grid
, ql_ad(row_length,rows)        & ! adiabatic liquid water content
                                 ! at inversion (kg/kg)
, entrain_coef(row_length,rows) & ! entrainment factor
, qsat_lcl(row_length,rows)       ! qsat at cloud base (kg/kg)

REAL ::                                                                 &
  wstar(row_length, rows)                                               &
                               ! surface-based mixed layer
!                                    ! velocity scale
      , wthvs(row_length, rows)                                         &
                                     ! surface buoyancy flux
      , w_max(row_length, rows) ! max w in column

REAL, ALLOCATABLE :: theta_star_surf(:,:) ! surface flux potential temp scale
REAL, ALLOCATABLE :: qv_star_surf(:,:)    ! surface flux moisture scale

REAL ::                                                                 &
  sub_surf_roff(land_points)                                            &
                              ! sub-surface runoff
, snomlt_sub_htf(land_points)
                              ! subsurface snowmelt heat flux
                              ! after hydrol

REAL ::                                                                 &
  snow_melt(land_points)                                                &
                          ! snowmelt (kg/m2/s).
, surf_roff(land_points)                                                &
                          ! surface runoff (kg/m2/s).
, tot_tfall(land_points)  ! total throughfall (kg/m2/s).

! Fields passed from BDY_LAYR to IMP_SOLVER
REAL ::                                                                 &
   rhokh (row_length, rows, bl_levels)

! Additional fields
REAL ::                                                                 &
 tstar_land(row_length,rows)                                            &
                             ! Land mean sfc temperature (K)
,tstar_sea(row_length,rows)                                             &
                             ! Open sea sfc temperature (K).
,tstar_sice_cat(row_length,rows,nice_use)                               &
                             ! Sea-ice category sfc temperature (K).
,tstar_sice(row_length,rows)                                            &
                             ! Sea-ice sfc temperature (K).
,tstar_ssi(row_length,rows)  ! Sea mean sfc temperature (K).
! diagnostics
REAL ::                                                                 &
      ! output from bdy_layr.
  fqw(row_length, rows, bl_levels)                                      &
, ftl(row_length, rows, bl_levels)                                      &
, taux(udims%i_start:udims%i_end,udims%j_start:udims%j_end,             &
       0:bl_levels-1)                                                   &
, tauy(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,             &
       0:bl_levels-1)                                                   &
, fb_surf(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                          ! Surface flux buoyancy over density (m^2/s^3)

! data required for tracer mixing :
REAL ::                                                                 &
  drydep2(row_length,rows,ndiv)
                                !dry dep though grav. set.

! Fields passed out of boundary layer into hydrology
REAL ::                                                                 &
       !
  snowmelt(row_length, rows)                                            &
                                  !output from sf_evap.
, ext(land_points,sm_levels) ! Extraction of water from each
!                                    soil layer (kg/m2/s).

REAL ::                                                                 &
       ! ( needed for soil_hft )
  surf_ht_flux_land(row_length, rows)                                   &
                                       !
, surf_ht_flux_ld(land_points)                                          &
                                       !
, snomlt_surf_htf(row_length, rows)

REAL ::                                                                 &
  cca_2d(row_length, rows)                                              &
, lclf

REAL ::                                                           &
  tot_precip_scaled_1(row_length, rows)                           &
, tot_precip_scaled_2(row_length, rows)

! logicals
LOGICAL ::                                                              &
  L_zero_boundaries                                                     &
, L_zero_halos
!
LOGICAL ::                                                              &
  L_scrn                                                                &
                           ! Logical to control output
                           !    of screen level T,Q,QCL,QCF
, L_plsp                   ! Logical to control output
                           !    of Probability of LS Precip

LOGICAL ::                                                              &
  L_calc_dxek                                                           &
               ! Switch for calculation of condensate increment
, L_q_interact ! Switch allows overwriting of parcel variables
!                      when calculating condensate increments.
LOGICAL ::                                                              &
  L_cape_opt_345 !Logical to store whether cldbase_opt_dp has value
                 ! 3, 4, 5 now also 6.

! Additional variables for MOSES II
REAL ::                                                                 &
  ftl_tile(land_points,ntiles)                                          &
                               ! Surface FTL for land tiles
, radnet_sea(row_length,rows)                                           &
                               ! Surface net radiation on
                               !     open sea (W/m2)
, radnet_sice(row_length,rows,nice_use)                                 &
                               ! Surface net radiation on
                               !     sea-ice (W/m2)
, fqw_tile(land_points,ntiles)                                          &
                               ! Surface FQW for land tiles
, epot_tile(land_points,ntiles)                                         &
                               ! OUT Local EPOT for land tiles.
, fqw_ice(row_length,rows,nice_use)                                     &
                               ! Surface FQW for sea-ice
, ftl_ice(row_length,rows,nice_use)                                     &
                               ! Surface FTL for sea-ice
, taux_land(udims%i_start:udims%i_end,udims%j_start:udims%j_end)        &
                               ! Taux over land part of grid box.
, taux_ssi(udims%i_start:udims%i_end,udims%j_start:udims%j_end)         &
                               ! Taux over sea part of grid box.
, tauy_land(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)        &
                               ! Tauy over land part of grid box.
, tauy_ssi(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)

REAL ::                                                                 &
 ei_tile(land_points,ntiles)                                            &
                             ! EI for land tiles
,ecan_tile(land_points,ntiles)                                          &
                          ! ECAN for land tiles
,melt_tile(land_points,ntiles)                                          &
!                               ! Snowmelt on tiles (kg/m2/s).
,snow_soil_htf(land_points,ntiles)                                      &
                                ! Tiled snow->soil heat flux (W/m2)
,snow_smb_surft(land_points,ntiles)                                     &
                                ! Tiled snow mass change (kg/m2/s)
,surf_htf_tile(land_points,ntiles)
!                               ! Net downward surface heat flux
!                               ! on tiles (W/m2)

INTEGER ::                                                              &
 phenol_call                                                            &
              ! indicates whether phenology is to be called
,triffid_call                                                           &
              ! indicates whether TRIFFID is to be called
,nstep_trif   ! Number of atmospheric timesteps between calls to
!                   ! TRIFFID vegetation model.

! Additional variables for JULES

REAL ::                                                                 &
 dtstar_sea(row_length,rows)
REAL ::                                                                 &
 dtstar_sice(row_length,rows,nice_use)
REAL ::                                                                 &
 dtstar_tile(land_points,ntiles)
                             ! Change in TSTAR over timestep
!                                  ! for land tiles
REAL ::                                                                 &
 snowdepth_p(    land_points,ntiles)                                    &
                       ! Snow depth on ground on tiles (m)
,rho_snow_grnd_p(land_points,ntiles)                                    &
                       ! Snowpack bulk density (kg/m3)
,nsnow_p(        land_points,ntiles)                                    &
                       ! Number of snow layers on ground on tiles
                       ! NOTE that this is converted to an integer.
,ds_p(        land_points,ntiles,nsmax)                                 &
                       ! Snow layer thickness (m)
,sice_p(      land_points,ntiles,nsmax)                                 &
                       ! Snow layer ice mass on tiles (Kg/m2)
,sliq_p(      land_points,ntiles,nsmax)                                 &
                       ! Snow layer liquid mass on tiles (Kg/m2)
,tsnowlayer_p(land_points,ntiles,nsmax)                                 &
                       ! Snow layer temperature (K)
,rho_snow_p(  land_points,ntiles,nsmax)                                 &
                       ! Snow layer densities (kg/m3)
,rgrainl_p(   land_points,ntiles,nsmax)
                       ! Snow layer grain size on tiles (microns)
! prognostic FLake fields
REAL ::                                                                 &
  lake_depth_p( land_points)                                            &
, lake_fetch_p( land_points)                                            &
, lake_t_mean_p(land_points)                                            &
, lake_t_mxl_p( land_points)                                            &
, lake_t_ice_p( land_points)                                            &
, lake_h_mxl_p( land_points)                                            &
, lake_h_ice_p( land_points)                                            &
, lake_shape_p( land_points)                                            &
, lake_g_dt_p(  land_points)

! Sea ice fields
REAL ::                                                                 &
 ice_fract_cat_use(row_length, rows, nice_use)                          &
                         ! Sea ice fraction passed to ni_bl_ctl
,k_sice(row_length, rows, nice)
                         ! Sea ice effective conductivity

! River routing:
REAL :: tot_surf_runoff(land_points)! Accumulated runoff over
REAL :: tot_sub_runoff(land_points) ! river routing timestep
                                   ! (Kg/m2/s)
REAL :: inlandout_atm(land_points)         ! INLAND BASIN FLOW
!                                  land points only kg/m2/s

! Water conservation
! Remove global lake evaporation from soil moisture store
REAL :: acc_lake_evap(row_length,rows)
!                                       ! Accumulated lake evap over
!                                       ! river routing timestep
!                                       ! (Kg/m2)

REAL ::                                                                 &
 w_copy(row_length,rows,0:model_levels)  ! copy of w to pass to
                                         ! conv_diag, BL?

!  Workspace :-

!  Local scalars :-

! Variables for multivariate swapbounds
INTEGER :: i_field
TYPE(swapable_field_pointer_type) :: fields_to_swap(7)

REAL ::                                                                 &
  mag_vector_np (model_levels)                                          &
, dir_vector_np (model_levels)                                          &
, mag_vector_sp (model_levels)                                          &
, dir_vector_sp (model_levels)

REAL ::                                                                 &
 pptn_rate(row_length,rows)                                             &
                              ! Total precipitation
                              ! (convective+large scale) (kg/m2/s)
,accum_pptn(row_length,rows)                                            &
                              ! Accumulated total precip (kg/m2)
,tot_rain(row_length,rows)                                              &
                              ! Total rain (kg/m2/s)
,tot_snow(row_length,rows)    ! Total snow (kg/m2/s)

! Extra variable for JULES snow scheme
LOGICAL, PARAMETER :: stf_hf_snow_melt = .TRUE.

REAL, ALLOCATABLE ::                     &
  r_u_p(:,:,:)                           & ! r_u on P-grid
 ,r_v_p(:,:,:)                           & ! r_v on P-grid
 ,ustar_p(:,:,:)                         & ! u+r_u on P-grid
 ,vstar_p(:,:,:)                           ! v+r_v on P-grid

! Work arrays for convection wind variables
REAL, TARGET, ALLOCATABLE ::             & ! Convective momentum transport
  dubydt_conv_p(:,:,:)                   & ! tendencies on p-grid, calculated
 ,dvbydt_conv_p(:,:,:)                     ! in the convection scheme.
REAL, ALLOCATABLE ::                     &
  dubydt_conv_u(:,:,:)                   & ! du/dt on u grid
 ,dvbydt_conv_v(:,:,:)                   & ! dv/dt on v grid
 ,work_u_halo(:,:,:)                     & ! work array for u
 ,work_v_halo(:,:,:)                       ! work array for v

! "steering level" properties for convective downdraughts
REAL :: u_steer(  pdims%i_start:pdims%i_end,                            &
                  pdims%j_start:pdims%j_end)                            &
!                                   ! U wind at cloud base
      , v_steer(  pdims%i_start:pdims%i_end,                            &
                  pdims%j_start:pdims%j_end)                            &
!                                   ! V wind at cloud base
      , cca_steer(pdims%i_start:pdims%i_end,                            &
                  pdims%j_start:pdims%j_end)
!                                   ! cnv cloud amnt at cloud base [0-1]

! New varibles added for message passing
REAL, TARGET ::                                                         &
       flandfac(pdims_s%i_start:pdims_s%i_end,                          &
                pdims_s%j_start:pdims_s%j_end),                         &
       fseafac(pdims_s%i_start:pdims_s%i_end,                           &
              pdims_s%j_start:pdims_s%j_end),                           &
       rhokm_land(pdims_s%i_start:pdims_s%i_end,                        &
                  pdims_s%j_start:pdims_s%j_end),                       &
       rhokm_ssi(pdims_s%i_start:pdims_s%i_end,                         &
                 pdims_s%j_start:pdims_s%j_end),                        &
       cdr10m(pdims_s%i_start:pdims_s%i_end,                            &
              pdims_s%j_start:pdims_s%j_end),                           &
       tau_fd_x(pdims_s%i_start:pdims_s%i_end,                          &
                pdims_s%j_start:pdims_s%j_end,0:bl_levels-1),           &
       tau_fd_y(pdims_s%i_start:pdims_s%i_end,                          &
                pdims_s%j_start:pdims_s%j_end,0:bl_levels-1)

REAL ::                                                                 &
  flandfac_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),      &
  flandfac_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),      &
  fseafac_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),       &
  fseafac_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),       &
  taux_fd_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,        &
            0:bl_levels-1),                                             &
  tauy_fd_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,        &
            0:bl_levels-1),                                             &
  rhokm_u_land(udims%i_start:udims%i_end,udims%j_start:udims%j_end),    &
  rhokm_u_ssi(udims%i_start:udims%i_end,udims%j_start:udims%j_end),     &
  rhokm_v_land(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),    &
  rhokm_v_ssi(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)

! 1A BL Scheme uses rho * gamma for non-gradient stress
REAL, TARGET, ALLOCATABLE :: rhogamu(:,:,:),                            &
                             rhogamv(:,:,:)
REAL, ALLOCATABLE :: rhogamu_u(:,:,:),                                  &
                     rhogamv_v(:,:,:)
! Other BL schemes use dimensionless function for non-gradient stress
REAL, TARGET, ALLOCATABLE :: f_ngstress(:,:,:)
REAL, ALLOCATABLE :: f_ngstress_u(:,:,:),                               &
                     f_ngstress_v(:,:,:)

! Required for convective sub-stepping with rediagnosis

LOGICAL ::                                                              &
  cumulus_copy(row_length, rows)  & ! Copy of cumulus from conv_diag
 ,no_cumulus(row_length, rows)      ! Points changed from cumulus to
                                    ! non-cumulus by BL routine

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=60) :: scheme_name    ! description for moisture checking

! Error reporting

CHARACTER(LEN=256)       :: message
INTEGER                  :: errorstatus

INTEGER, PARAMETER :: i_river_vn_1A = 1
INTEGER, PARAMETER :: i_river_vn_2A = 2


! Local
LOGICAL :: l_emcorr_tmp      ! Tmp scm logical
LOGICAL :: l_CallConvection
LOGICAL :: l_jules_call ! flag for whether the BL scheme is being called
                        ! for jules only or for the main scheme

LOGICAL :: l_extra_call = .FALSE.  ! true this is an additional
                                    ! call to conv_diag within a timestep

LOGICAL :: l_calc_at_p   ! Flag input to bdy_expl3 to indicate whether
                         ! arguments are on p-grid or native u and v grids.

! needed for ENDGame etadot calc below
REAL :: u_at_w,v_at_w

! Variables used with FLake and ML-snow
REAL :: dwsw_sub_snow   ! (W/m2) remaining SW flux under the snowpack
REAL :: dhf_surf_minus_soil(land_points)
                        ! (W/m2) heat flux difference across the snowpack

! needed for endgame for fv_cos_theta_latitude vs cos_theta_latitude
REAL, POINTER :: xx_cos_theta_latitude (:,:)

INTEGER :: ompt_start   ! start of thread's work
INTEGER :: ompt_end     ! end of thread's work

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

xx_cos_theta_latitude => cos_theta_latitude

! ----------------------------------------------------------------------
! Section INI. Initialisation of variables.
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Communication Section
! ----------------------------------------------------------------------

! allocate arrays if this is the first time AP2 has been called
IF ( .NOT. ALLOCATED(uhalo) ) THEN
  CALL atmos_physics2_alloc(land_points,ntiles,ndiv,ndivh,npft,ntype,   &
                            dim_cs1,dim_cs2,sm_levels,bl_levels,nice_use)
END IF

! only call on 1st cycle or if not fast running
IF ( cycleno == 1 .OR. .NOT. l_quick_ap2 ) THEN

  !   Include surface currents in coupling.
  !   Was not included for old HadGEM1.
  !   SURFACE CURRENT INTERPOLATION START-
  !   interpolate U_0,V_0 onto U_0_P,V_0_P

!$OMP PARALLEL DEFAULT(NONE)                                           &
!$OMP SHARED( udims, vdims, uhalo, u_0, vhalo, v_0 )                   &
!$OMP PRIVATE( j )
!$OMP DO SCHEDULE(STATIC)
  DO j=udims%j_start, udims%j_end
    uhalo(udims%i_start:udims%i_end,j) = u_0(udims%i_start:udims%i_end,j)
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO j=vdims%j_start, vdims%j_end
    vhalo(vdims%i_start:vdims%i_end,j) = v_0(vdims%i_start:vdims%i_end,j)
  END DO
!$OMP END DO
!$OMP END PARALLEL

  i_field = 0

  i_field = i_field + 1
  fields_to_swap(i_field) % field_2d   => uhalo(:,:)
  fields_to_swap(i_field) % field_type = fld_type_u
  fields_to_swap(i_field) % levels     = 1
  fields_to_swap(i_field) % rows       = rows
  fields_to_swap(i_field) % vector     = .TRUE.

  i_field = i_field + 1
  fields_to_swap(i_field) % field_2d   => vhalo(:,:)
  fields_to_swap(i_field) % field_type = fld_type_v
  fields_to_swap(i_field) % levels     = 1
  fields_to_swap(i_field) % rows       = n_rows
  fields_to_swap(i_field) % vector     = .TRUE.

  !   Need to call swap bounds as halo points not set
  CALL swap_bounds_2d_mv(fields_to_swap, i_field, row_length,           &
                         offx, offy)

END IF ! ( cycleno == 1 .OR. .NOT. l_quick_ap2 )

CALL v_to_p(vhalo,                                                      &
                vdims_s%i_start,vdims_s%i_end,                          &
                vdims_s%j_start,vdims_s%j_end,                          &
                pdims%i_start,pdims%i_end,                              &
                pdims%j_start,pdims%j_end,                              &
                1, at_extremity,v_0_p)

CALL u_to_p(uhalo,                                                      &
                udims_s%i_start,udims_s%i_end,                          &
                udims_s%j_start,udims_s%j_end,                          &
                pdims%i_start,pdims%i_end,                              &
                pdims%j_start,pdims%j_end,                              &
                1, at_extremity,u_0_p)

! only call on 1st cycle or if not fast running
IF ( cycleno == 1 .OR. .NOT. l_quick_ap2 ) THEN
  !
  !   Calculate etadot, requires comms for ND only
  !
  L_poles = .FALSE.
  j_begin = 1
  j_end = rows
  IF (model_type == mt_global) THEN
    IF (at_extremity(PSouth)) j_begin = 2
    IF (at_extremity(PNorth)) j_end = rows-1
    L_poles = .TRUE.
  END IF ! model_type == mt_global

  SELECT CASE (model_type)

  CASE (mt_single_column)
    DO j=1, rows
      DO i=1, row_length
        DO k=0, model_levels
          etadot_copy(i,j,k) = w(i,j,k) / z_top_of_model
        END DO
      END DO
    END DO

  CASE DEFAULT

    !------------------------------------------------------------------------
    !       Initialise etadot: Will be computed by the solver at timestep>1
    !------------------------------------------------------------------------
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(k,j,i,u_at_w,v_at_w)                  &
!$OMP SHARED(model_levels,pdims,intw_rho2w,intw_u2p,u,intw_v2p,v,          &
!$OMP    etadot_copy,w, h3_p_eta,dxi1_xi3,h1_p_eta,dxi2_xi3,h2_p_eta,      &
!$OMP    deta_xi3_theta)

!$OMP DO SCHEDULE(STATIC)
    DO k = 1, model_levels-1
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end

          u_at_w = intw_rho2w(k,1)*( intw_u2p(i,1)*u(i-1,j,k+1) +     &
                                     intw_u2p(i,2)*u(i,j,k+1) ) +     &
                   intw_rho2w(k,2)*( intw_u2p(i,1)*u(i-1,j,k) +       &
                                       intw_u2p(i,2)*u(i,j,k) )

          v_at_w = intw_rho2w(k,1)*( intw_v2p(j,1)*v(i,j-1,k+1) +     &
                                     intw_v2p(j,2)*v(i,j,k+1) ) +     &
                   intw_rho2w(k,2)*( intw_v2p(j,1)*v(i,j-1,k) +       &
                                       intw_v2p(j,2)*v(i,j,k) )

          etadot_copy(i,j,k) = ( w(i,j,k)/h3_p_eta(i,j,k) -           &
                             u_at_w*dxi1_xi3(i,j,k)/                  &
                                           h1_p_eta(i,j,k) -          &
                             v_at_w*dxi2_xi3(i,j,k)/                  &
                                           h2_p_eta(i,j,k) ) /        &
                                             deta_xi3_theta(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        etadot_copy(i,j,0) = 0.0
        etadot_copy(i,j,model_levels) = 0.0
      END DO
    END DO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

  END SELECT ! model_type

!$OMP  PARALLEL DEFAULT(NONE)                                 &
!$OMP SHARED(u,udims_s,pdims,model_levels,at_extremity,u_p,   &
!$OMP    vdims_s,v,v_p)

  !   u and v winds on all levels copied to p-grid, again only requires
  !   comms for ND runs
  CALL u_to_p(u,                                                        &
              udims_s%i_start,udims_s%i_end,                            &
              udims_s%j_start,udims_s%j_end,                            &
              pdims%i_start,pdims%i_end,                                &
              pdims%j_start,pdims%j_end,                                &
              model_levels, at_extremity,u_p)

  CALL v_to_p(v,                                                        &
              vdims_s%i_start,vdims_s%i_end,                            &
              vdims_s%j_start,vdims_s%j_end,                            &
              pdims%i_start,pdims%i_end,                                &
              pdims%j_start,pdims%j_end,                                &
              model_levels, at_extremity,v_p)

!$OMP END PARALLEL

END IF ! ( cycleno == 1 .OR. .NOT. l_quick_ap2 )

! Winds required for convection only if convective momentum being used
IF ( l_param_conv .AND. l_mom ) THEN

  ! Convection either needs r_u or u_star on the p-grid, depending
  ! on whether doing rediagnosis.  Only allocate space for the one
  ! that's needed.
  IF (l_rediagnosis) THEN
    ALLOCATE( r_u_p ( pdims%i_start:pdims%i_end,                        &
                      pdims%j_start:pdims%j_end,                        &
                      pdims%k_start:pdims%k_end ) )
    ALLOCATE( r_v_p ( pdims%i_start:pdims%i_end,                        &
                      pdims%j_start:pdims%j_end,                        &
                      pdims%k_start:pdims%k_end ) )
    ALLOCATE( ustar_p(1,1,1) )
    ALLOCATE( vstar_p(1,1,1) )
  ELSE
    ALLOCATE( r_u_p(1,1,1) )
    ALLOCATE( r_v_p(1,1,1) )
    ALLOCATE( ustar_p ( pdims%i_start:pdims%i_end,                      &
                        pdims%j_start:pdims%j_end,                      &
                        pdims%k_start:pdims%k_end ) )
    ALLOCATE( vstar_p ( pdims%i_start:pdims%i_end,                      &
                        pdims%j_start:pdims%j_end,                      &
                        pdims%k_start:pdims%k_end ) )
  END IF

  SELECT CASE (model_type)

  CASE DEFAULT

    !     CALL swapbounds on r_u and r_v.
    i_field = 0
    i_field = i_field + 1
    fields_to_swap(i_field) % field       => r_u(:,:,:)
    fields_to_swap(i_field) % field_type  =  fld_type_u
    fields_to_swap(i_field) % levels      =  model_levels
    fields_to_swap(i_field) % rows        =  rows
    fields_to_swap(i_field) % vector      =  .TRUE.

    i_field = i_field + 1
    fields_to_swap(i_field) % field       => r_v(:,:,:)
    fields_to_swap(i_field) % field_type  =  fld_type_v
    fields_to_swap(i_field) % levels      =  model_levels
    fields_to_swap(i_field) % rows        =  n_rows
    fields_to_swap(i_field) % vector      =  .TRUE.

    CALL swap_bounds_mv(fields_to_swap, i_field, row_length,            &
                         offx, offy)

    IF (l_rediagnosis) THEN

      CALL u_to_p (r_u, udims_s%i_start,udims_s%i_end,                  &
                    udims_s%j_start,udims_s%j_end,                      &
                    pdims%i_start,pdims%i_end,                          &
                    pdims%j_start,pdims%j_end,                          &
                    model_levels, at_extremity, r_u_p)

      CALL v_to_p (r_v, vdims_s%i_start,vdims_s%i_end,                  &
                    vdims_s%j_start,vdims_s%j_end,                      &
                    pdims%i_start,pdims%i_end,                          &
                    pdims%j_start,pdims%j_end,                          &
                    model_levels, at_extremity, r_v_p)


    ELSE     ! non rediagnosis option

      ALLOCATE(work_u_halo(udims_s%i_start:udims_s%i_end,               &
                           udims_s%j_start:udims_s%j_end,               &
                           udims_s%k_start:udims_s%k_end) )
      ALLOCATE(work_v_halo(vdims_s%i_start:vdims_s%i_end,               &
                           vdims_s%j_start:vdims_s%j_end,               &
                           vdims_s%k_start:vdims_s%k_end) )

      ! Work out ustar=u+du on u grid and then interpolate to p_grid

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(k,j,i,ompt_start,ompt_end)    &
!$OMP SHARED(u,udims_s,pdims,model_levels,at_extremity,r_v,        &
!$OMP    ustar_p,vdims_s,v,vstar_p,r_u,work_u_halo,work_v_halo)

      ! The following loops and sub routines are all parallelised using
      ! ompt_start and ompt_end as calculated by compute_chunk_size
      ! since the variables udims, vdims, pdims all have k_start = 1
      ! and k_end = model_levels

      ompt_start = 1
      ompt_end   = model_levels

      ! only call compute_chunk_size if compiling with OMP
      ! The procedure call is protected by the optional compile
      ! sentinel
!$ CALL compute_chunk_size(1,model_levels,ompt_start,ompt_end)

      DO k = ompt_start,ompt_end
        DO j = udims_s%j_start, udims_s%j_end
          DO i = udims_s%i_start,udims_s%i_end
            work_u_halo(i,j,k) = u(i,j,k) + r_u(i,j,k)
          END DO ! i
        END DO ! j
      END DO ! k

      CALL u_to_p (work_u_halo,                                     &
                    udims_s%i_start,udims_s%i_end,                  &
                    udims_s%j_start,udims_s%j_end,                  &
                    pdims%i_start,pdims%i_end,                      &
                    pdims%j_start,pdims%j_end,                      &
                    model_levels, at_extremity, ustar_p)

      DO k = ompt_start,ompt_end
        DO j = vdims_s%j_start, vdims_s%j_end
          DO i = vdims_s%i_start,vdims_s%i_end
            work_v_halo(i,j,k) = v(i,j,k) + r_v(i,j,k)
          END DO ! i
        END DO ! j
      END DO ! k

      CALL v_to_p (work_v_halo,                                     &
                    vdims_s%i_start,vdims_s%i_end,                  &
                    vdims_s%j_start,vdims_s%j_end,                  &
                    pdims%i_start,pdims%i_end,                      &
                    pdims%j_start,pdims%j_end,                      &
                    model_levels, at_extremity, vstar_p)

!$OMP END PARALLEL

      DEALLOCATE(work_v_halo)
      DEALLOCATE(work_u_halo)

    END IF ! l_rediagnosis

    !     SCM version u and v on p grid  (No problems with v having less rows)

  CASE (mt_single_column)
    ! SCM version u and v on p grid
    ! (No problems with v having less rows)
    IF (l_rediagnosis) THEN

      DO k = 1, tdims%k_end
        DO j = 1, rows
          DO i = 1, row_length
            r_u_p(i,j,k) = r_u(i,j,k)
            r_v_p(i,j,k) = r_v(i,j,k)
          END DO ! i
        END DO ! j
      END DO ! k

    ELSE

      DO k = 1, tdims%k_end
        DO j = 1, rows
          DO i = 1, row_length
            ustar_p(i,j,k) = u(i,j,k) + r_u(i,j,k)
            vstar_p(i,j,k) = v(i,j,k) + r_v(i,j,k)
          END DO ! i
        END DO ! j
      END DO ! k

    END IF ! test on l_rediagnosis

  END SELECT ! model_type

  ! Allocate arrays for convective momentum tendencies on p-grid (with haloes
  ! for interpolation onto u/v grids later)
  ALLOCATE( dubydt_conv_p( pdims_s%i_start:pdims_s%i_end,           &
                           pdims_s%j_start:pdims_s%j_end,           &
                           pdims_s%k_start:pdims_s%k_end ) )
  ALLOCATE( dvbydt_conv_p( pdims_s%i_start:pdims_s%i_end,           &
                           pdims_s%j_start:pdims_s%j_end,           &
                           pdims_s%k_start:pdims_s%k_end ) )

END IF  ! l_param_conv .AND. l_mom

! Calculate extended halo land fraction
! Surface currents and level 1 winds needed for coastal tiling only

! only call on 1st cycle or if not fast running
IF ( cycleno == 1 .OR. .NOT. l_quick_ap2 ) THEN

  IF ( l_ctile ) THEN
    DO l=1,land_points
      fland(l)=fland_ctile(l)
    END DO
  ELSE
    DO l=1,land_points
      fland(l)=1.0
    END DO
  END IF ! L_CTILE

  !-----------------------------------------------------------------------
  ! Expand land fraction to global field:
  !-----------------------------------------------------------------------
!$OMP PARALLEL DEFAULT(NONE) PRIVATE( i, j, k )                         &
!$OMP SHARED( rows, row_length, bl_levels, flandg, u_px, u_p, v_px, v_p,&
!$OMP         u_0_px, u_0_p, v_0_px, v_0_p )

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, bl_levels
    DO j = 1, rows
      DO i = 1, row_length
        u_px(i,j,k) = u_p(i,j,k)
        v_px(i,j,k) = v_p(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO j=1,rows
    DO i=1,row_length
      flandg(i,j)=0.0
      u_0_px(i,j) = u_0_p(i,j)
      v_0_px(i,j) = v_0_p(i,j)
    END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL

  DO l=1,land_points
    j=(land_index(l)-1)/row_length + 1
    i = land_index(l) - (j-1)*row_length
    flandg(i,j)=fland(l)
  END DO

  IF (model_type /= mt_single_column) THEN
    i_field = 0
    i_field = i_field + 1
    fields_to_swap(i_field) % field_2d   => flandg(:,:)
    fields_to_swap(i_field) % field_type = fld_type_p
    fields_to_swap(i_field) % levels     = 1
    fields_to_swap(i_field) % rows       = rows
    fields_to_swap(i_field) % vector     = .FALSE.
  END IF ! model_type

  IF (model_type /= mt_single_column) THEN

    IF (l_ctile .AND. buddy_sea == on) THEN

      i_field = i_field + 1
      fields_to_swap(i_field) % field_2d   => u_0_px(:,:)
      fields_to_swap(i_field) % field_type = fld_type_p
      fields_to_swap(i_field) % levels     = 1
      fields_to_swap(i_field) % rows       = rows
      fields_to_swap(i_field) % vector     = .FALSE.

      i_field = i_field + 1
      fields_to_swap(i_field) % field_2d   => v_0_px(:,:)
      fields_to_swap(i_field) % field_type = fld_type_p
      fields_to_swap(i_field) % levels     = 1
      fields_to_swap(i_field) % rows       = rows
      fields_to_swap(i_field) % vector     = .FALSE.

    END IF  ! test on buddy_sea switch

    CALL swap_bounds_2d_mv(fields_to_swap, i_field, row_length,         &
                           offx , offy )

    IF (l_ctile .AND. buddy_sea == on) THEN

      i_field = 0
      i_field = i_field + 1
      fields_to_swap(i_field) % field      => u_px(:,:,:)
      fields_to_swap(i_field) % field_type = fld_type_p
      fields_to_swap(i_field) % levels     = bl_levels
      fields_to_swap(i_field) % rows       = rows
      fields_to_swap(i_field) % vector     = .FALSE.

      i_field = i_field + 1
      fields_to_swap(i_field) % field      => v_px(:,:,:)
      fields_to_swap(i_field) % field_type = fld_type_p
      fields_to_swap(i_field) % levels     = bl_levels
      fields_to_swap(i_field) % rows       = rows
      fields_to_swap(i_field) % vector     = .FALSE.


      CALL swap_bounds_mv(fields_to_swap, i_field, row_length,          &
                           offx , offy )


    END IF  ! test on buddy_sea switch

  END IF ! model_type

END IF ! ( cycleno == 1 .OR. .NOT. l_quick_ap2 )

! only call on 1st cycle
IF ( cycleno == 1 ) THEN

  ! Represent propagation of convective cold pools.
  IF ( cnv_cold_pools > ccp_off ) THEN

    DO j =   pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        ! steering wind fields to be added later
        u_steer(  i,j) = 0.0
        v_steer(  i,j) = 0.0
        ! only use cloud amount if the cloud level is non-zero
        IF (ccb0(i,j) > 0) THEN
          cca_steer(i,j) = cca0(i,j,ccb0(i,j))
        ELSE
          cca_steer(i,j) = 0.0
        END IF
      END DO
    END DO

    IF (ISrfExCnvGust == IP_SrfExWithCnv) THEN

      CALL conv_cold_pools( global_row_length, global_rows, ddmfx,      &
                            u_steer, v_steer, cca_steer,                &
                            ux_ccp, uy_ccp, um_ccp, g_ccp, h_ccp,       &
                            riso_ccp, rdir_ccp )

    ELSE
      cmessage =                                                        &
       'ISrfExCnvGust must be set for CCP, resetting Cnv_Cold_Pools = 0'
      cnv_cold_pools = ccp_off
      icode=-35
      CALL ereport(RoutineName, icode, cmessage)
    END IF

  END IF ! ( cnv_cold_pools > ccp_off )

END IF ! ( cycleno == 1 )


! ----------------------------------------------------------------------
! End of communication Section
! ----------------------------------------------------------------------

L_apply_diag = CycleNo == NumCycles

SELECT CASE (model_type)
CASE DEFAULT
  ! Hydrology
  !   Sub-surface runoff must be calculated if running
  !   with river routing or model will fail.
  stf_sub_surf_roff = sf(205,8) .OR. sf(235,8) .OR. l_rivers
  smlt = sf(202,8)

CASE (mt_single_column)
  stf_sub_surf_roff=.TRUE.
  smlt=.TRUE.

  IF (.NOT. l_emcorr_opt) THEN
    l_emcorr_tmp = l_emcorr
    l_emcorr = l_emcorr_opt
  END IF ! l_emcorr_opt

END SELECT ! model_type

!-----------------------------------------------------------------------
! map JULES prognostics to module prognostics
!-----------------------------------------------------------------------
snowdepth_surft = snowdepth_p
IF ( nsmax > 0 ) THEN
  rho_snow_grnd_surft = rho_snow_grnd_p
  nsnow_surft         = nsnow_p
  ds_surft            = ds_p
  sice_surft          = sice_p
  sliq_surft          = sliq_p
  tsnow_surft         = tsnowlayer_p
  rho_snow_surft      = rho_snow_p
  rgrainl_surft       = rgrainl_p
ELSE
  rho_snow_grnd_surft = rho_snow_const
  nsnow_surft         = 0
END IF

! Give the JULES prognostics the values of the UM prognostics.
! If it is the first timestep of the run JULES will be passed
! the values read in from the start dump. If not there will
! be no change because the UM prognostics are updated to the
! JULES prognostics after SURF_COUPLE_EXTRA

! JULES prognostic = UM prognostic

IF (l_triffid) THEN

  ! Wood product pool prognostics
  IF (l_landuse) THEN
    wood_prod_fast_gb = wood_prod_fast_d1
    wood_prod_med_gb  = wood_prod_med_d1
    wood_prod_slow_gb = wood_prod_slow_d1
  END IF

  ! Nitrogen prognostics
  IF (l_nitrogen) THEN
    ns_pool_gb(:,1,1) =  soil_nitro1(:)
    ns_pool_gb(:,1,2) =  soil_nitro2(:)
    ns_pool_gb(:,1,3) =  soil_nitro3(:)
    ns_pool_gb(:,1,4) =  soil_nitro4(:)
    deposition_n_gb(:) = nitrogen_deposition_d1(:)
    ! Variable n_inorg_soilt_lyrs was introduced with the layered soil carbon
    ! scheme. When L_layeredC =.FALSE., as is currently the case in all UM+JULES
    ! runs because layered soil carbon is available only in standalone JULES at
    ! present, the top layer of that variable is used. However the existing
    ! inorganic N prognostic, n_inorg_gb, is also still used as the sum
    ! over all levels. Both need to be initialised with the inorganic soil N
    ! prognostic (stash 446) read in from the start dumps (soil_inorgnit),
    ! otherwise the inorganic N pool gets reset to 0 at the start of each run.
    n_inorg_gb(:)             = soil_inorgnit(:)
    n_inorg_soilt_lyrs(:,1,1) = soil_inorgnit(:)
  END IF ! l_nitrogen

  ! Disturbed fraction or crop and pasture fraction prognostics:
  !   If L_TRIF_CROP is TRUE then the disturbed fraction variable is crop
  !     fraction only, with a separate prognostic for pasture.
  !   If L_TRIF_CROP is FALSE then there is only one agricultural fraction as
  !     in HadGEM2ES and that is only used if L_LANDUSE is true.
  IF (l_trif_crop) THEN
    ! The disturbed fraction variable is used for crop fraction
    frac_agr_prev_gb  = agr_crop_frac_prev_d1
    ! Pasture fractions
    frac_past_gb      = pasture_frac_d1
    frac_past_prev_gb = pasture_frac_prev_d1
  ELSE
    ! The disturbed fraction variable contains the disturbed fraction
    ! as always
    IF (l_landuse) THEN
      frac_agr_prev_gb  = disturb_veg_prev
    END IF
  END IF

  IF ( l_co2_interactive ) THEN

      !
      ! If TRIFFID has not yet been called during this run then it won't have
      ! calculated triffid_co2_gb, which is the variable passed to 
      ! BL_TRMIX_DD for adding to the atmosphere every timestep. Therefore we 
      ! need to give triffid_co2_gb the prognostic version read in from the 
      ! start dump (triffid_co2_d1), which was calculated by TRIFFID at the 
      ! very end of the previous CRUN.
      !
      ! If TRIFFID has been called this run, triffid_co2_d1 would have 
      ! been updated with the new value of triffid_co2_gb after 
      ! SURF_COUPLE_EXTRA, for storing in the start dump, so this line 
      ! effectively does nothing.    
      !

    triffid_co2_gb(:) =  triffid_co2_d1(:)

  END IF ! l_co2_interactive

END IF ! L_TRIFFID

IF ( l_nitrogen ) THEN
  ! Give JULES prognostics N the values of the UM prognostics.
  ns_pool_gb(:,1,1) =  soil_nitro1(:)
  ns_pool_gb(:,1,2) =  soil_nitro2(:)
  ns_pool_gb(:,1,3) =  soil_nitro3(:)
  ns_pool_gb(:,1,4) =  soil_nitro4(:)
  n_inorg_gb(:)     = soil_inorgnit(:)
  deposition_n_gb(:) = nitrogen_deposition_d1(:)
END IF ! l_nitrogen

IF ( l_flake_model                         &
     .AND. (.NOT. l_aggregate   )           &
     .AND. (land_points > 0    ) ) THEN
  !
  ! map FLake prognostics to module prognostics
  lake_depth_gb        = lake_depth_p
  lake_fetch_gb        = lake_fetch_p
  lake_t_mean_gb       = lake_t_mean_p
  lake_t_mxl_gb        = lake_t_mxl_p
  lake_t_ice_gb        = lake_t_ice_p
  lake_h_mxl_gb        = lake_h_mxl_p
  lake_h_ice_gb        = lake_h_ice_p
  lake_shape_factor_gb = lake_shape_p
  g_dt_gb              = lake_g_dt_p

  ! initialise FLake variables needed in surface exchange
  DO l = 1, land_points
    ! set surface T based on T* of the lake tile
    lake_t_sfc_gb(l) = tstar_tile(l, lake)

    ! set the FLake snow depth,
    ! BUT depending on the presence of ice
    lake_h_snow_gb(l) = 0.0
    IF (lake_h_ice_gb(l) > 0.0) THEN
      lake_h_snow_gb(l) = snowdepth_surft(l, lake)
    END IF

    ! set the snow surface T based on T* of the lake tile
    lake_t_snow_gb(l) = tstar_tile(l, lake)
  END DO

END IF ! l_flake_model

IF (l_param_conv) THEN     ! A convection scheme is being called

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(k,j,i)                  &
!$OMP SHARED(model_levels,rows,row_length,ccw_out,ccw0,      &
!$OMP    lcbase_out,lcbase0,l_ccrad,                         &
!$OMP    rad_cloud_decay_opt,ccb0,cct0,cca0_2d,cclwp0,       &
!$OMP    cca0,l_cca_dp_prog,cca0_dp,l_cca_md_prog,           &
!$OMP    cca0_md,l_cca_sh_prog,cca0_sh,ccb,cct,cca_2d,       &
!$OMP    cclwp,lcbase,lctop,lcca,ccw,cca,n_cca_levels)

  ! Copy ccw/lcbase values for Radiation if l_ccrad in use.
  ! See comments at declaration of ccw_out/lcbase_out
  IF (l_ccrad) THEN
!$OMP DO SCHEDULE(STATIC)
    DO k=1, model_levels
      DO j=1, rows
        DO i=1, row_length
          ccw0(i,j,k) = ccw_out(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
    DO j=1, rows
      DO i=1, row_length
        lcbase0(i,j)  = lcbase_out(i,j)
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF ! l_ccrad

  IF (rad_cloud_decay_opt == rad_decay_off) THEN
    ! Zero section 0 convective cloud diagnostics going to radiation
!$OMP DO SCHEDULE(STATIC)
    DO j = 1, rows
      DO i = 1, row_length
        ccb0(i,j)    = 0
        cct0(i,j)    = 0
        cca0_2d(i,j) = 0.0
        cclwp0(i,j)  = 0.0
        lcbase0(i,j) = 0
      END DO
    END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
    DO k=1, model_levels
      DO j=1, rows
        DO i=1, row_length
          ccw0(i,j,k) = 0.0
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
    DO k=1, n_cca_levels
      DO j=1, rows
        DO i=1, row_length
          cca0(i,j,k) = 0.0
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

    IF (l_cca_dp_prog) THEN
!$OMP DO SCHEDULE(STATIC)
      DO k=1, n_cca_levels
        DO j=1, rows
          DO i=1, row_length
            cca0_dp(i,j,k) = 0.0
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT
    END IF
    IF (l_cca_md_prog) THEN
!$OMP DO SCHEDULE(STATIC)
      DO k=1, n_cca_levels
        DO j=1, rows
          DO i=1, row_length
            cca0_md(i,j,k) = 0.0
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT
    END IF
    IF (l_cca_sh_prog) THEN
!$OMP DO SCHEDULE(STATIC)
      DO k=1, n_cca_levels
        DO j=1, rows
          DO i=1, row_length
            cca0_sh(i,j,k) = 0.0
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT
    END IF
  END IF ! Rad decay disabled

!$OMP DO SCHEDULE(STATIC)
  DO j = 1, rows
    DO i = 1, row_length
      ccb(i,j)    = 0
      cct(i,j)    = 0
      cca_2d(i,j) = 0.0
      cclwp(i,j)  = 0.0
      lcbase(i,j) = 0
      lctop(i,j)  = 0
      lcca(i,j)   = 0.0
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO k=1, model_levels
    DO j=1, rows
      DO i=1, row_length
        ccw(i,j,k) = 0.0
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO k=1, n_cca_levels
    DO j=1, rows
      DO i=1, row_length
        cca(i,j,k) = 0.0
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


ELSE   ! no call to a convection scheme

  ! There should be no testing against convection namelist switches as these
  ! should have no impact.
  ! All convective cloud fields should be zero if they are passed to other 
  ! parts of the code. These are the fields with zero in their name.
  ! Other convective cloud fields should not be in use.
  ! Note if l_param_conv = .false. it is assumed l_ccrad=.false. so
  ! ccw_out and lcbase_out should not be required or set.

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(k,j,i)                  &
!$OMP SHARED(model_levels,rows,row_length, n_cca_levels,     &
!$OMP    lcbase0,lcbase,conv_rain,conv_snow,                 &
!$OMP    ccb0,cct0,cca0_2d,cclwp0,cca0,ccw0,cca)

!$OMP DO SCHEDULE(STATIC)
  DO j = 1, rows
    DO i = 1, row_length
      conv_rain(i,j) = 0.0
      conv_snow(i,j) = 0.0
      ccb0(i,j)    = 0
      cct0(i,j)    = 0
      cca0_2d(i,j) = 0.0
      cclwp0(i,j)  = 0.0
      lcbase0(i,j) = 0
      lcbase(i,j)  = 0        ! needs to be set as used in ni_imp_ctl
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO k=1, model_levels
    DO j=1, rows
      DO i=1, row_length
        ccw0(i,j,k) = 0.0
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO k=1, n_cca_levels
    DO j=1, rows
      DO i=1, row_length
        cca0(i,j,k) = 0.0
        cca(i,j,k)  = 0.0      ! Needs to be set as used in  ni_imp_ctl
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

END IF    ! l_param_conv

! set temporary logicals to disabled un-called physics.
IF (i_bl_vn == i_bl_vn_0 .AND. l_param_conv) THEN
  WRITE(umMessage,*)' convection on and boundary layer off,'
  CALL umPrint(umMessage,src='atmos_physics2')
  WRITE(umMessage,*)' not a sensible choice'
  CALL umPrint(umMessage,src='atmos_physics2')
  WRITE(umMessage,*)' will try and run but results may be garbage'
  CALL umPrint(umMessage,src='atmos_physics2')
END IF


!$OMP PARALLEL DEFAULT(NONE)                                            &
!$OMP SHARED(p, p_layer_boundaries, p_layer_centres, p_theta_levels,    &
!$OMP exner_layer_boundaries, exner_layer_centres, exner_rho_levels,    &
!$OMP exner_theta_levels, model_levels, rows, row_length,               &
!$OMP z_theta, z_rho, r_theta_levels, r_rho_levels, rho_wet, rho_dry,   &
!$OMP rho_wet_rsq, unscaled_dry_rho, rho_wet_tq, halo_i, halo_j,        &
!$OMP rho_dry_theta, w_copy, w, zh_prev, zh, p_star, p_zero, kappa,     &
!$OMP bulk_cloud_fraction, bulk_cloud_fraction_halos,                   &
!$OMP cloud_fraction_liquid, cloud_fraction_liquid_halos,               &
!$OMP cloud_fraction_frozen, cloud_fraction_frozen_halos, tdims)        &
!$OMP PRIVATE(i,j,k)

! set p at layer boundaries.
! NB: some arrays have haloes but are unset, if never used they will
!     be removed.
!$OMP DO SCHEDULE(STATIC)
DO j = 1, rows
  DO i = 1, row_length
    zh_prev(i,j) = zh(i,j)  ! make a copy of zh, as passed in
    p_layer_boundaries(i,j,0) = p_star(i,j)
    p_layer_centres(i,j,0) = p_star(i,j)
    exner_layer_boundaries(i,j,0) = (p_layer_boundaries(i,j,0)/         &
                                     p_zero)**kappa
    exner_layer_centres(i,j,0) = exner_layer_boundaries(i,j,0)

  END DO
END DO
!$OMP END DO NOWAIT

! k=model_levels never appears on the LHS, and is not modified by this loop.
! There are instances of p( , ,k+1) and exner_rho_levels( , ,k+1) in the first
! loop, but those arrays are not modified by the second. Hence nowait valid.
!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels - 1
  DO j = 1, rows
    DO i = 1, row_length
      p_layer_boundaries(i,j,k) = p(i,j,k+1)
      p_layer_centres(i,j,k)    = p_theta_levels(i,j,k)
      exner_layer_boundaries(i,j,k) = exner_rho_levels(i,j,k+1)
      exner_layer_centres(i,j,k)    = exner_theta_levels(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

k = model_levels

!$OMP DO SCHEDULE(STATIC)
DO j = 1, rows
  DO i = 1, row_length
    p_layer_boundaries(i,j,model_levels) = 0.0
    p_layer_centres(i,j,model_levels)    = p_theta_levels(i,j,model_levels)
    exner_layer_boundaries(i,j,k) = 0.0
    exner_layer_centres(i,j,k)    = exner_theta_levels(i,j,k)
  END DO
END DO
!$OMP END DO NOWAIT
!
! Variables required by conv_diag and convection_control
! Heights of model levels, density; all without halos
!

!$OMP DO SCHEDULE(STATIC)
DO k =                 1, tdims%k_end
  DO j =   tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      z_theta(i,j,k) = r_theta_levels(i,j,k) - r_theta_levels(i,j,0)
      z_rho(i,j,k)   = r_rho_levels(i,j,k)   - r_theta_levels(i,j,0)
      rho_wet(i,j,k) = rho_wet_rsq(i,j,k) /( r_rho_levels(i,j,k)        &
                                        *r_rho_levels(i,j,k) )
      rho_dry(i,j,k) = unscaled_dry_rho(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO
!Implicit barrier

! density on theta levels (not top level).

! Following routine is thread-safe and can be called from a parallel region
CALL p_to_t(row_length,rows, halo_i, halo_j,0,0, model_levels-1,             &
            r_theta_levels(:,:,0:model_levels),   &
            r_rho_levels(:,:,1:model_levels),     &
            rho_wet, rho_wet_tq)

! dry density on theta levels (not top level).

! Following routine is thread-safe and can be called from a parallel region
CALL p_to_t(row_length,rows, halo_i,halo_j,0,0, model_levels-1,             &
            r_theta_levels(:,:,0:model_levels),   &
            r_rho_levels(:,:,1:model_levels),     &
            rho_dry, rho_dry_theta)

!$OMP DO SCHEDULE(STATIC)
DO k =                 1, tdims%k_end
  DO j =   tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      bulk_cloud_fraction(i,j,k)   = bulk_cloud_fraction_halos(i,j,k)
      cloud_fraction_liquid(i,j,k) = cloud_fraction_liquid_halos(i,j,k)
      cloud_fraction_frozen(i,j,k) = cloud_fraction_frozen_halos(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

! Copy of w for use in conv_diag and other routines.

!$OMP DO SCHEDULE(STATIC)
DO k=0,model_levels
  DO j=1,rows
    DO i=1,row_length
      w_copy(i,j,k) = w(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

IF (model_type == mt_lam .AND. L_lbc_old ) THEN

  ! zero qcl, qcf and bulk_cloud_fraction on LAM boundaries to avoid
  ! inconsistencies between cloud fields and cloud amounts.
  ! As results on boundaries do not affect model results this can be done.

  L_zero_boundaries=.TRUE.
  L_zero_halos=.FALSE.

  ! DEPENDS ON: zero_lateral_boundaries
  CALL zero_lateral_boundaries(                                         &
   row_length,rows,halo_i,halo_j,model_levels,fld_type_p,qcl,           &
   1, at_extremity,                                                     &
   L_zero_boundaries,L_zero_halos)

  ! DEPENDS ON: zero_lateral_boundaries
  CALL zero_lateral_boundaries(                                         &
   row_length,rows,halo_i,halo_j,model_levels,fld_type_p,qcf,           &
   1, at_extremity,                                                     &
   L_zero_boundaries,L_zero_halos)

  ! DEPENDS ON: zero_lateral_boundaries
  CALL zero_lateral_boundaries(                                         &
   row_length,rows,0,0,model_levels,fld_type_p,                         &
   bulk_cloud_fraction,                                                 &
   1, at_extremity,                                                     &
   L_zero_boundaries,L_zero_halos)

END IF !  model_type == mt_lam .AND. L_lbc_old

! ----------------------------------------------------------------------
! Set Coastal tiling dependent prognostics:

! ----------------------------------------------------------------------
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k)

IF ( l_ctile ) THEN

!$OMP DO SCHEDULE(STATIC)
  DO j = 1, rows
    DO i = 1, row_length
      tstar_land(i,j)=tstar_land_ctile(i,j)
      tstar_sea(i,j)=tstar_sea_ctile(i,j)
      tstar_sice(i,j)=tstar_sice_ctile(i,j)
      IF (ice_fract(i,j) <= 0.0) THEN
        tstar_ssi(i,j) = tstar_sea(i,j)
      ELSE
        tstar_ssi(i,j)=ice_fract(i,j)*tstar_sice(i,j)                   &
                +(1.0-ice_fract(i,j))*tstar_sea(i,j)
      END IF
    END DO
  END DO
!$OMP END DO NOWAIT

  DO k = 1, nice_use
!$OMP DO SCHEDULE(STATIC)
    DO j = 1, rows
      DO i = 1, row_length
        tstar_sice_cat(i,j,k)=tstar_sice_cat_ctile(i,j,k)
      END DO
    END DO
!$OMP END DO NOWAIT
  END DO

ELSE ! l_ctile

!$OMP DO SCHEDULE(STATIC)
  DO j = 1, rows
    DO i = 1, row_length
      tstar_land(i,j)=tstar(i,j)
      tstar_ssi(i,j)=tstar(i,j)

      IF (.NOT. land_sea_mask(i,j)) THEN
        IF (ice_fract(i,j) <= 0.0) THEN
          tstar_sea(i,j)=tstar(i,j)
          tstar_sice(i,j)=tstar(i,j)
        ELSE
          tstar_sea(i,j)=tfs
          tstar_sice(i,j)=(tstar(i,j)                                   &
            -(1.0-ice_fract(i,j))*tstar_sea(i,j))/ice_fract(i,j)
        END IF
      ELSE
        tstar_sea(i,j)=tstar(i,j)
        tstar_sice(i,j)=tstar(i,j)
      END IF

    END DO
  END DO
!$OMP END DO NOWAIT

  !   Note that nice_use=1 if not using coastal tiling
  DO k = 1, nice_use
!$OMP DO SCHEDULE(STATIC)
    DO j = 1, rows
      DO i = 1, row_length
        tstar_sice_cat(i,j,k)=tstar_sice(i,j)
      END DO
    END DO
!$OMP END DO NOWAIT
  END DO

END IF ! l_ctile

! Set up sea ice fraction and conductivity fields depending on science choices

IF (l_sice_multilayers) THEN  !Note nice_use=nice in this case
  DO k = 1, nice
!$OMP DO SCHEDULE(STATIC)
    DO j = 1, rows
      DO i = 1, row_length
        k_sice(i,j,k) = k_sice_ml(i,j,k)  ! Received from sea ice model
      END DO
    END DO
!$OMP END DO NOWAIT
  END DO
ELSE    ! Note nice may not equal nice_use
  DO k = 1, nice
!$OMP DO SCHEDULE(STATIC)
    DO j = 1, rows
      DO i = 1, row_length
        ! Set sea ice effective conductivity to constant value
        k_sice(i,j,k) = 2.0*kappai/de
      END DO
    END DO
!$OMP END DO NOWAIT
  END DO
END IF

IF (nice_use > 1) THEN   ! Use categories fully
  DO k = 1, nice_use
!$OMP DO SCHEDULE(STATIC)
    DO j = 1, rows
      DO i = 1, row_length
        ice_fract_cat_use(i,j,k) = ice_fract_ncat(i,j,k)
      END DO
    END DO
!$OMP END DO NOWAIT
  END DO
ELSE
!$OMP DO SCHEDULE(STATIC)
  DO j = 1, rows
    DO i = 1, row_length
      ice_fract_cat_use(i,j,1) = ice_fract(i,j)
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

!$OMP END PARALLEL

!----------------------------------
! compute surface fluxes here
!----------------------------------
SELECT CASE (model_type)

CASE DEFAULT
  L_plsp=sf(281,3) .OR. sf(282,3) .OR. sf(283,3)
  L_scrn=L_plsp

CASE (mt_single_column)
  L_plsp=.TRUE.
  L_scrn=L_plsp

END SELECT ! model_type

IF ( i_bl_vn /= i_bl_vn_0 ) THEN
  !
  !   only call on 1st cycle or if not fast running
  IF ( cycleno == 1 .OR. .NOT. l_quick_ap2 ) THEN
    ! set local switch for whether model is using CLASSIC aerosol scheme
    l_aero_classic = (l_sulpc_so2 .OR. l_soot .OR. l_biomass .OR.       &
                      l_ocff .OR. l_nitrate .OR. l_dust)

    ! if cycleno > 1 and not l_quick_ap2 (already tested above), then we
    ! need to reset sice_frac_ncat to ice_fract as was done in atmos_phys1
    IF ( cycleno > 1 ) THEN
      sice_pts_ncat(:)=0
      sice_index_ncat(:,:)=0
      sice_frac_ncat(:,:)=0.0
      DO n=1,nice_use
        DO l=1,ssi_pts
          j=(ssi_index(l)-1)/row_length + 1
          i = ssi_index(l) - (j-1)*row_length
          IF (ssi_index(l) > 0) THEN
            IF (ice_fract(i,j) > 0.0) THEN
              sice_pts_ncat(n)=sice_pts_ncat(n)+1
              sice_index_ncat(sice_pts_ncat(n),n)=l
              sice_frac_ncat(l,n)=ice_fract(i,j)
            END IF
          END IF
        END DO
      END DO
    END IF ! cycleno > 1

    !     allocate variables depending on which bl scheme is used
    IF ( i_bl_vn == i_bl_vn_1a ) THEN
      ALLOCATE(rhogamu(pdims_s%i_start:pdims_s%i_end,                   &
                       pdims_s%j_start:pdims_s%j_end ,0:bl_levels-1))
      ALLOCATE(rhogamv(pdims_s%i_start:pdims_s%i_end,                   &
                       pdims_s%j_start:pdims_s%j_end ,0:bl_levels-1))
      ALLOCATE(f_ngstress(1,1,1))
    ELSE
      ALLOCATE(rhogamu(1,1,1))
      ALLOCATE(rhogamv(1,1,1))
      ALLOCATE(f_ngstress(pdims_s%i_start:pdims_s%i_end,                &
                          pdims_s%j_start:pdims_s%j_end,1:bl_levels-1))
    END IF ! i_bl_vn

    !     Allocation of arrays for neutral wind diagnostics.
    !     At least initially, we expect equivalent neutral diagnostics
    !     to be required as a block.

    SELECT CASE (model_type)

    CASE (mt_single_column)
      sf_diag%suv10m_n = .FALSE.

    CASE DEFAULT
      sf_diag%suv10m_n = ( sf(368,3) .OR. sf(369,3) .OR.                &
                   sf(370,3) .OR. sf(371,3) .OR.                        &
                   sf(365,3) .OR. sf(366,3) .OR. sf(367,3) ) .AND.      &
                   (cycleno == numcycles .OR. l_quick_ap2)

    END SELECT
    !
    IF (sf_diag%suv10m_n) THEN
      !       Full sizes required for diagnostics arrays for neutral winds.
      !       (These are declared in atmos_physics2_alloc, but are not
      !       allocated there since the relevent STASH flags may not be
      !       available on the first timestep.)
      ALLOCATE(cdr10m_n(pdims_s%i_start:pdims_s%i_end,                  &
                        pdims_s%j_start:pdims_s%j_end))
      ALLOCATE(cdr10m_n_u(udims%i_start:udims%i_end,                    &
                          udims%j_start:udims%j_end))
      ALLOCATE(cdr10m_n_v(vdims%i_start:vdims%i_end,                    &
                          vdims%j_start:vdims%j_end))
      ALLOCATE(cd10m_n(pdims_s%i_start:pdims_s%i_end,                   &
                       pdims_s%j_start:pdims_s%j_end))
      ALLOCATE(cd10m_n_u(udims%i_start:udims%i_end,                     &
                         udims%j_start:udims%j_end))
      ALLOCATE(cd10m_n_v(vdims%i_start:vdims%i_end,                     &
                         vdims%j_start:vdims%j_end))
    ELSE
      !       Minimal sizes are alllocated to ensure valid addresses are
      !       present in subroutine calls.
      ALLOCATE(cdr10m_n(1,1))
      ALLOCATE(cdr10m_n_u(1,1))
      ALLOCATE(cdr10m_n_v(1,1))
      ALLOCATE(cd10m_n(1,1))
      ALLOCATE(cd10m_n_u(1,1))
      ALLOCATE(cd10m_n_v(1,1))
    END IF

    IF (l_jules_flux) THEN
      l_jules_call = .TRUE.

      CALL NI_bl_ctl (                                                  &
      ! IN parameters for SISL scheme
        CycleNo,l_jules_call,                                           &
      ! IN time stepping information
        val_year, val_day_number, val_hour, val_minute, val_second,     &
      ! IN model dimensions.
        land_points, ntiles, bl_levels,                                 &
      ! IN switches
        L_scrn, L_aero_classic,                                         &
      ! IN data fields.
        p, p_layer_centres, rho_wet_rsq,rho_wet,rho_dry, u_p, v_p,      &
        u_px, v_px, u_0_px, v_0_px,                                     &
        land_sea_mask,q,qcl,qcf,p_star,theta,exner_theta_levels,rad_hr, &
        micro_tends, soil_layer_moisture, rho_wet_tq, z_rho, z_theta,   &
      ! IN ancillary fields and fields needed to be kept from tstep to tstep
        hcon,smvccl_levs,smvcwt_levs,smvcst_levs,sthf,sthu,sil_orog_land,&
      !-------------------------------------------------------------------------
        ho2r2_orog, sd_orog, ice_fract_cat_use, k_sice(:,:,1:nice_use), &
        land_index, photosynth_act_rad,                                 &
        soil_clay,soil_sand,dust_mrel1,dust_mrel2,                      &
        dust_mrel3,dust_mrel4,dust_mrel5,dust_mrel6,                    &
      ! IN additional variables for JULES
        canopy, catch, catch_snow, snow_tile, z0_tile, z0h_tile_bare,   &
        z0m_soil, lw_down, sw_tile, tstar_tile, tsurf_elev_surft,       &
        co2(1:co2_dim_len,1:co2_dim_row,1),                             &
        asteps_since_triffid,                                           &
        cs,frac,canht_ft,lai_ft,fland,flandg,albsoil,cos_zenith_angle,  &
      ! IN everything not covered so far
        t_soil,ti_gb,                                                   &
        ti,tstar,zh_prev,ddmfx,bulk_cloud_fraction,zhpar,zlcl,          &
      ! IN SCM namelist data
        L_spec_z0, z0m_scm, z0h_scm, flux_e, flux_h, ustar_in,          &
      ! SCM diagnostics and STASH
        nSCMDpkgs, L_SCMDiags, BL_diag, sf_diag,                        &
      ! INOUT data
        gs,z0msea,w_copy,etadot_copy,tstar_sea,tstar_sice_cat,zh,dzh,   &
        cumulus,ntml,ntpar,l_shallow,error_code,                        &
      ! INOUT additional variables for JULES
        g_leaf_acc,npp_ft_acc,resp_w_ft_acc,resp_s_acc,                 &
      ! INOUT variables for TKE based turbulence schemes
        e_trb, tsq_trb, qsq_trb, cov_trb, zhpar_shcu,                   &
      ! INOUT variables from bdy_expl1 needed elsewhere
        bq_gb, bt_gb, dtrdz_charney_grid,rdz_charney_grid,              &
        dtrdz_u, dtrdz_v, rdz_u, rdz_v, k_blend_tq, k_blend_uv,         &
      ! INOUT variables from Jules needed elsewhere
        flandfac,fseafac,rhokm_land,rhokm_ssi,cdr10m,cdr10m_n,cd10m_n,  &
        fqw, ftl, rib_gb, vshr, z0m_eff_gb, z0h_eff_gb, r_b_dust,       &
        rho_aresist,aresist,resist_b, rhokm,rhokh,                      &
      ! INOUT diagnostics required for soil moisture nudging scheme :
        wt_ext,                                                         &
      ! INOUT variables required in IMP_SOLVER
        alpha1_sea, alpha1_sice, ashtf_sea, ashtf, uStarGBM,            &
      ! INOUT additional variables for JULES
        ftl_tile,radnet_sea,radnet_sice,rib_tile,rho_aresist_tile,      &
        aresist_tile,resist_b_tile,alpha1,ashtf_tile,fqw_tile,epot_tile,&
        fqw_ice,ftl_ice,fraca,resfs,resft,rhokh_tile,rhokh_sice,rhokh_sea,&
        z0hssi,z0h_tile,z0m_gb,z0mssi,z0m_tile,chr1p5m,chr1p5m_sice,smc,&
        gpp,npp,resp_p,g_leaf,gpp_ft,npp_ft,resp_p_ft,resp_s,resp_s_tot,&
        resp_w_ft,gc,canhc_tile,wt_ext_tile,flake,tile_index,tile_pts,  &
        tile_frac,fsmc,vshr_land,vshr_ssi,tstar_land,tstar_ssi,dtstar_tile,&
        dtstar_sea,dtstar_sice,hcons,emis_tile,emis_soil,t1_sd,q1_sd,fb_surf,&
      ! OUT variables for message passing
        tau_fd_x, tau_fd_y, rhogamu, rhogamv, f_ngstress,               &
      ! OUT diagnostics (done after implicit solver)
        zht, zhnl, shallowc,cu_over_orog,bl_type_1,bl_type_2,bl_type_3, &
        bl_type_4,bl_type_5,bl_type_6, bl_type_7, bl_w_var,             &
      ! OUT variables required for mineral dust scheme
        dust_flux,dust_emiss_frac, u_s_t_tile,u_s_t_dry_tile,           &
        u_s_std_tile, kent, we_lim, t_frac, zrzi,                       &
        kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhsc,               &
      ! OUT fields
        nbdsc,ntdsc,wstar,wthvs,uw0,vw0,taux_p,tauy_p,rhcpt, rib_ssi    &
        )
    END IF !l_jules_flux

  END IF !l_quick_ap2

END IF !bl_vn > 0
! ---------------------------------------------
! CALL CONV_DIAG to diagnose convection
! ---------------------------------------------

! latest values needed for substepping/fully sequential BL.

! Determine value of L_cape_opt_345
L_cape_opt_345 = .FALSE.   ! default whether a convection scheme called or not

IF (l_param_conv) THEN

  ! Use  convection switches to decide the value of  L_cape_opt_345
  IF (i_convection_vn == i_convection_vn_5a .OR.                          &
         i_convection_vn == i_convection_vn_6a ) THEN
    L_cape_opt_345 = ( (cldbase_opt_dp == 3) .OR. (cldbase_opt_md == 3) .OR. &
                       (cldbase_opt_dp == 4) .OR. (cldbase_opt_md == 4) .OR. &
                       (cldbase_opt_dp == 5) .OR. (cldbase_opt_md == 5) .OR. &
                       (cldbase_opt_dp == 6) .OR. (cldbase_opt_md == 6) )
  END IF
END IF

! Start of OpenMP parallel region
!$OMP PARALLEL DEFAULT(NONE)                                                 &
!$OMP SHARED(rows, row_length,                                               &
!$OMP ntml, ntpar, nlcl, cumulus, l_shallow, l_mid_level, delthvu,           &
!$OMP ql_ad, zhpar, dzh, qcl_inv_top, zlcl, zlcl_uv,                         &
!$OMP cape_bottom, L_cape_opt_345, cape_top, w_copy, conv_type,              &
!$OMP no_cumulus, w_max)                                                     &
!$OMP PRIVATE(i,j,k)

! Initialise conv_diag output arrays
!$OMP DO SCHEDULE(STATIC)
DO j = 1, rows
  DO i = 1, row_length
    ntml(i,j)       = 1
    ntpar(i,j)      = 1
    nlcl(i,j)       = 1
    cumulus(i,j)    = .FALSE.
    no_cumulus(i,j) = .FALSE.
    l_shallow(i,j)  = .FALSE.
    l_mid_level(i,j)= .FALSE.
    conv_type(i,j)  = 0
    delthvu(i,j)    = 0.0
    ql_ad(i,j)      = 0.0
    zhpar(i,j)      = 0.0
    dzh(i,j)        = 0.0
    qcl_inv_top(i,j)= 0.0
    zlcl(i,j)       = 0.0
    zlcl_uv(i,j)    = 0.0

    ! Initialise the w_max array.
    w_max(i,j)      = 0.0

  END DO
END DO
!$OMP END DO NOWAIT

IF (L_cape_opt_345) THEN
  !   Find w_max for each column. The w_max array is initialised just
  !   before the start of the OpenMP parallel region.
  DO k =  cape_bottom, cape_top
!$OMP DO SCHEDULE(STATIC)
    DO j = 1, rows
      DO i = 1, row_length
        w_max(i,j) = MAX(w_max(i,j), w_copy(i,j,k))
      END DO
    END DO
!$OMP END DO NOWAIT
  END DO
END IF  ! L_cape_opt_345

! End of OpenMP parallel region
!$OMP END PARALLEL
!

! only call on 1st cycle or if not fast running
IF ( cycleno == 1 .OR. .NOT. l_quick_ap2 ) THEN
  ! Uses beginning of timestep values for Theta, q, qcl, qcf ,u and v etc

  SELECT CASE ( i_convection_vn )

  CASE ( i_convection_vn_5a )

    CALL conv_diag_5a(                                                  &

    !     IN Parallel variables
            row_length, rows                                            &

    !     IN model dimensions.
          , bl_levels                                                   &
          , p, p_layer_centres(1,1,1),exner_rho_levels                  &
          , rho_wet, rho_wet_tq, z_theta, z_rho                         &

    !     IN Model switches
          , l_extra_call                                                &
          , no_cumulus                                                  &

    !     IN cloud data
          , qcf(1:row_length,1:rows,1:tdims%k_end)                      &
          , qcl(1:row_length,1:rows,1:tdims%k_end), bulk_cloud_fraction &

    !     IN everything not covered so far :

          , p_star, q(1:row_length,1:rows,1:tdims%k_end)                &
          , theta(tdims%i_start:tdims%i_end,                            &
              tdims%j_start:tdims%j_end,1:tdims%k_end)                  &
          , exner_theta_levels(tdims%i_start:tdims%i_end,               &
              tdims%j_start:tdims%j_end, 1:tdims%k_end)                 &
          , u_p, v_p, u_0_p, v_0_p                                      &
          , tstar_land, tstar_sea, tstar_sice, z0msea                   &
          , flux_e, flux_h, ustar_in, L_spec_z0, z0m_scm, z0h_scm       &
          , tstar, land_sea_mask, flandg, ice_fract, timestep           &
          , w_copy, w_max, deep_flag, past_precip, past_conv_ht         &
    !
    !     IN surface fluxes
          , fb_surf, ustarGBM                                           &
    !     SCM Diagnostics (dummy values in full UM)
          , nSCMDpkgs,L_SCMDiags                                        &

    !     OUT data required elsewhere in UM system :
          , zh,zhpar,dzh,qcl_inv_top,zlcl,zlcl_uv,delthvu,ql_ad         &
          , ntml,ntpar,nlcl                                             &
          , cumulus,l_shallow,l_congestus,l_congestus2                  &
          , conv_type                                                   &
          , CIN_undilute,CAPE_undilute, wstar, wthvs                    &
          , entrain_coef, qsat_lcl                                      &
          , Error_code                                                  &
           )

  ! Preferred method of convective diagonsis. Used by 6A scheme and called
  ! by any other convection schemes which do not require diagnosis info.
  ! The BL scheme requires information from a conv_diag call.
  ! Other convection schemes with number >9 are allowed
  ! Currently only option 11 is an allowed value but the test below allows
  ! other schemes to be tested in a branch.

  CASE ( i_convection_vn_6a, 10:20 )

    CALL conv_diag_6a(                                                  &

    !     IN Parallel variables
            row_length, rows                                            &

    !     IN model dimensions.
          , bl_levels                                                   &
          , p, p_layer_centres(1,1,1),exner_rho_levels                  &
          , rho_wet, rho_wet_tq, z_theta, z_rho                         &

    !     IN Model switches
          , l_extra_call                                                &
          , no_cumulus                                                  &

    !     IN cloud data
          , qcf(1:row_length,1:rows,1:tdims%k_end)                      &
          , qcl(1:row_length,1:rows,1:tdims%k_end), bulk_cloud_fraction &

    !     IN everything not covered so far :

          , p_star, q(1:row_length,1:rows,1:tdims%k_end)                &
          , theta(tdims%i_start:tdims%i_end,                            &
              tdims%j_start:tdims%j_end,1:tdims%k_end)                  &
          , exner_theta_levels(tdims%i_start:tdims%i_end,               &
              tdims%j_start:tdims%j_end, 1:tdims%k_end)                 &
          , u_p, v_p, u_0_p, v_0_p                                      &
          , tstar_land, tstar_sea, tstar_sice, z0msea                   &
          , flux_e, flux_h, ustar_in, L_spec_z0, z0m_scm, z0h_scm       &
          , tstar, land_sea_mask, flandg, ice_fract                     &
          , w_copy, w_max, deep_flag, past_precip, past_conv_ht         &
          , conv_prog_precip                                            &
          , g_ccp, h_ccp                                                &
    !
    !     IN surface fluxes
          , fb_surf, ustarGBM                                           &
    !     SCM Diagnostics (dummy values in full UM)
          , nSCMDpkgs,L_SCMDiags                                        &

    !     OUT data required elsewhere in UM system :
          , zh,zhpar,dzh,qcl_inv_top,zlcl,zlcl_uv,delthvu,ql_ad         &
          , ntml,ntpar,nlcl                                             &
          , cumulus,l_shallow,l_congestus,l_congestus2                  &
          , conv_type                                                   &
          , CIN_undilute,CAPE_undilute, wstar, wthvs                    &
          , entrain_coef, qsat_lcl                                      &
          , Error_code                                                  &
            )

  CASE DEFAULT ! i_convection_vn

    errorstatus = 10
    WRITE (message,'(A)') 'Convection scheme version value not recognised'
    WRITE (message,'(A,I6)') '   i_convection_vn = ',i_convection_vn
    CALL Ereport ( RoutineName, errorstatus, message)

  END SELECT ! i_convection_vn

  IF (l_quick_ap2) THEN
    !     save outputs for second EG cycle
    !     N.B. any new OUT data added to conv_diag will need to be saved here
    !     and restored below. The array size will also need to be increased in
    !     atmos_physics2_alloc.

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,i) SCHEDULE(STATIC)
    DO j = 1, rows
      DO i = 1, row_length
        conv_diag_reals(i,j,1)=zh(i,j)
        conv_diag_reals(i,j,2)=zhpar(i,j)
        conv_diag_reals(i,j,3)=zlcl(i,j)
        conv_diag_reals(i,j,4)=zlcl_uv(i,j)
        conv_diag_reals(i,j,5)=delthvu(i,j)
        conv_diag_reals(i,j,6)=ql_ad(i,j)
        conv_diag_reals(i,j,7)=cin_undilute(i,j)
        conv_diag_reals(i,j,8)=cape_undilute(i,j)
        conv_diag_reals(i,j,9)=wstar(i,j)
        conv_diag_reals(i,j,10)=wthvs(i,j)
        conv_diag_reals(i,j,11)=entrain_coef(i,j)
        conv_diag_reals(i,j,12)=qsat_lcl(i,j)
        conv_diag_reals(i,j,13)=dzh(i,j)
        conv_diag_reals(i,j,14)=qcl_inv_top(i,j)

        conv_diag_ints(i,j,1)=ntml(i,j)
        conv_diag_ints(i,j,2)=ntpar(i,j)
        conv_diag_ints(i,j,3)=nlcl(i,j)
        conv_diag_ints(i,j,4)=conv_type(i,j)

        conv_diag_logs(i,j,1)=cumulus(i,j)
        conv_diag_logs(i,j,2)=l_shallow(i,j)
        conv_diag_logs(i,j,3)=l_congestus(i,j)
        conv_diag_logs(i,j,4)=l_congestus2(i,j)
      END DO
    END DO
!$OMP END PARALLEL DO

  END IF ! l_quick_ap2

ELSE ! ( cycleno == 1 .OR. .NOT. l_quick_ap2 )

  ! restore outputs on second EG cycle

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,i) SCHEDULE(STATIC)
  DO j = 1, rows
    DO i = 1, row_length
      zh(i,j)=conv_diag_reals(i,j,1)
      zhpar(i,j)=conv_diag_reals(i,j,2)
      zlcl(i,j)=conv_diag_reals(i,j,3)
      zlcl_uv(i,j)=conv_diag_reals(i,j,4)
      delthvu(i,j)=conv_diag_reals(i,j,5)
      ql_ad(i,j)=conv_diag_reals(i,j,6)
      cin_undilute(i,j)=conv_diag_reals(i,j,7)
      cape_undilute(i,j)=conv_diag_reals(i,j,8)
      wstar(i,j)=conv_diag_reals(i,j,9)
      wthvs(i,j)=conv_diag_reals(i,j,10)
      entrain_coef(i,j)=conv_diag_reals(i,j,11)
      qsat_lcl(i,j)=conv_diag_reals(i,j,12)
      dzh(i,j)=conv_diag_reals(i,j,13)
      qcl_inv_top(i,j)=conv_diag_reals(i,j,14)

      ntml(i,j)=conv_diag_ints(i,j,1)
      ntpar(i,j)=conv_diag_ints(i,j,2)
      nlcl(i,j)=conv_diag_ints(i,j,3)
      conv_type(i,j)=conv_diag_ints(i,j,4)

      cumulus(i,j)=conv_diag_logs(i,j,1)
      l_shallow(i,j)=conv_diag_logs(i,j,2)
      l_congestus(i,j)=conv_diag_logs(i,j,3)
      l_congestus2(i,j)=conv_diag_logs(i,j,4)
    END DO
  END DO
!$OMP END PARALLEL DO

END IF ! ( cycleno == 1 .OR. .NOT. l_quick_ap2 )

! Make a copy of the cumulus array so that can tell which points are altered
! by the boundary layer scheme.

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( rows, row_length, cumulus_copy, cumulus, wstar, wthvs )  &
!$OMP PRIVATE( i, j )
DO j = 1, rows
  DO i = 1, row_length
    cumulus_copy(i,j) = cumulus(i,j)

    ! Initialise arrays which will be passed from BL to convection
    wstar(i,j) = 0.0
    wthvs(i,j) = 0.0
  END DO
END DO
!$OMP END PARALLEL DO

!----------------------------------------------------------------------
! Section BL CALL Explicit part of Boundary Layer scheme.
! ---------------------------------------------------------------------

IF (Ltimer) CALL timer ('AP2 Boundary Layer',5)
IF (Ltimer) CALL timer ('AP2 Explicit BL',5)

IF ( i_bl_vn /= i_bl_vn_0 ) THEN
  !
  !   only call on 1st cycle or if not fast running
  IF ( cycleno == 1 .OR. .NOT. l_quick_ap2 ) THEN

    l_jules_call = .FALSE.

    CALL NI_bl_ctl (                                                    &
    ! IN parameters for SISL scheme
      CycleNo,l_jules_call,                                             &
    ! IN time stepping information
      val_year, val_day_number, val_hour, val_minute, val_second,       &
    ! IN model dimensions.
      land_points, ntiles, bl_levels,                                   &
    ! IN switches
      L_scrn, L_aero_classic,                                           &
    ! IN data fields.
      p, p_layer_centres, rho_wet_rsq,rho_wet,rho_dry, u_p, v_p,        &
      u_px, v_px, u_0_px, v_0_px,                                       &
      land_sea_mask,q,qcl,qcf,p_star,theta,exner_theta_levels,rad_hr,   &
      micro_tends, soil_layer_moisture, rho_wet_tq, z_rho, z_theta,     &
    ! IN ancillary fields and fields needed to be kept from tstep to tstep
      hcon,smvccl_levs,smvcwt_levs,smvcst_levs,sthf,sthu,sil_orog_land, &
    !-------------------------------------------------------------------------
      ho2r2_orog, sd_orog, ice_fract_cat_use, k_sice(:,:,1:nice_use),   &
      land_index, photosynth_act_rad,                                   &
      soil_clay,soil_sand,dust_mrel1,dust_mrel2,                        &
      dust_mrel3,dust_mrel4,dust_mrel5,dust_mrel6,                      &
    ! IN additional variables for JULES
      canopy, catch, catch_snow, snow_tile, z0_tile, z0h_tile_bare,     &
      z0m_soil, lw_down, sw_tile, tstar_tile, tsurf_elev_surft,         &
      co2(1:co2_dim_len,1:co2_dim_row,1),                               &
      asteps_since_triffid,                                             &
      cs,frac,canht_ft,lai_ft,fland,flandg,albsoil,cos_zenith_angle,    &
    ! IN everything not covered so far
      t_soil,ti_gb,                                                     &
      ti,tstar,zh_prev,ddmfx,bulk_cloud_fraction,zhpar,zlcl,            &
    ! IN SCM namelist data
      L_spec_z0, z0m_scm, z0h_scm, flux_e, flux_h, ustar_in,            &
    ! SCM diagnostics and STASH
      nSCMDpkgs, L_SCMDiags, BL_diag, sf_diag,                          &
    ! INOUT data
      gs,z0msea,w_copy,etadot_copy,tstar_sea,tstar_sice_cat,zh,dzh,     &
      cumulus,ntml,ntpar,l_shallow,error_code,                          &
    ! INOUT additional variables for JULES
      g_leaf_acc,npp_ft_acc,resp_w_ft_acc,resp_s_acc,                   &
    ! INOUT variables for TKE based turbulence schemes
      e_trb, tsq_trb, qsq_trb, cov_trb, zhpar_shcu,                     &
    ! INOUT variables from bdy_expl1 needed elsewhere
      bq_gb, bt_gb, dtrdz_charney_grid,rdz_charney_grid,                &
      dtrdz_u, dtrdz_v, rdz_u, rdz_v, k_blend_tq, k_blend_uv,           &
    ! INOUT variables from Jules needed elsewhere
      flandfac,fseafac,rhokm_land,rhokm_ssi,cdr10m,cdr10m_n,cd10m_n,    &
      fqw, ftl, rib_gb, vshr, z0m_eff_gb, z0h_eff_gb, r_b_dust,         &
      rho_aresist,aresist,resist_b, rhokm,rhokh,                        &
    ! INOUT diagnostics required for soil moisture nudging scheme :
      wt_ext,                                                           &
    ! INOUT variables required in IMP_SOLVER
      alpha1_sea, alpha1_sice, ashtf_sea, ashtf, uStarGBM,              &
    ! INOUT additional variables for JULES
      ftl_tile,radnet_sea,radnet_sice,rib_tile,rho_aresist_tile,        &
      aresist_tile,resist_b_tile,alpha1,ashtf_tile,fqw_tile,epot_tile,  &
      fqw_ice,ftl_ice,fraca,resfs,resft,rhokh_tile,rhokh_sice,rhokh_sea,&
      z0hssi,z0h_tile,z0m_gb,z0mssi,z0m_tile,chr1p5m,chr1p5m_sice,smc,  &
      gpp,npp,resp_p,g_leaf,gpp_ft,npp_ft,resp_p_ft,resp_s,resp_s_tot,  &
      resp_w_ft,gc,canhc_tile,wt_ext_tile,flake,tile_index,tile_pts,    &
      tile_frac,fsmc,vshr_land,vshr_ssi,tstar_land,tstar_ssi,dtstar_tile,&
      dtstar_sea,dtstar_sice,hcons,emis_tile,emis_soil,t1_sd,q1_sd,fb_surf,&
    ! OUT variables for message passing
      tau_fd_x, tau_fd_y, rhogamu, rhogamv, f_ngstress,                 &
    ! OUT diagnostics (done after implicit solver)
      zht, zhnl, shallowc,cu_over_orog,bl_type_1,bl_type_2,bl_type_3,   &
      bl_type_4,bl_type_5,bl_type_6, bl_type_7, bl_w_var,               &
    ! OUT variables required for mineral dust scheme
      dust_flux,dust_emiss_frac, u_s_t_tile,u_s_t_dry_tile,             &
      u_s_std_tile, kent, we_lim, t_frac, zrzi,                         &
      kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhsc,                 &
    ! OUT fields
      nbdsc,ntdsc,wstar,wthvs,uw0,vw0,taux_p,tauy_p,rhcpt, rib_ssi      &
      )

    IF (l_quick_ap2) THEN
      ! save outputs for second EG cycle
      ! any new OUT data added to ni_bl_ctl which gets subsequently modified
      ! (i.e. is INOUT to conv_ctl or imp_ctl) will need to be saved here
      ! and restored below. The array size will also need to be increased in
      ! atmos_physics2_alloc

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,i,k)

!$OMP DO SCHEDULE(STATIC)
      DO j = 1, rows
        DO i = 1, row_length
          bl_ctl_2d(i,j,1)=zh(i,j)    !INOUT to bl_ctl, gets reset in atm_step
          bl_ctl_2d(i,j,2)=z0msea(i,j)!INOUT to bl_ctl, gets reset in atm_step
          bl_ctl_2d(i,j,3)=wstar(i,j) !INOUT to conv_ctl, needs reset
          bl_ctl_2d(i,j,4)=wthvs(i,j) !INOUT to conv_ctl, needs reset

          bl_ctl_int2d(i,j,1)=ntml(i,j)!INOUT to conv_ctl, needs reset

          bl_ctl_log2d(i,j,1)=cumulus(i,j)
          !INOUT to conv_ctl, needs reset
          bl_ctl_log2d(i,j,2)=l_shallow(i,j)
          !INOUT to conv_ctl, needs reset
        END DO
      END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
      DO k = 1, bl_levels
        DO j = 1, rows
          DO i = 1, row_length
            bl_ctl_3d(i,j,k,1)=fqw(i,j,k)  !INOUT to imp_ctl, needs reset
            bl_ctl_3d(i,j,k,2)=ftl(i,j,k)  !INOUT to imp_ctl, needs reset
            bl_ctl_3d(i,j,k,3)=rhokh(i,j,k)!INOUT to imp_ctl, needs reset
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT

      IF (l_use_dtstar_sea) THEN
!$OMP DO SCHEDULE(STATIC)
        DO j = 1, rows
          DO i = 1, row_length
            sea_save(i,j,1)=radnet_sea(i,j) !INOUT to imp_ctl, needs reset
            sea_save(i,j,2)=dtstar_sea(i,j) !INOUT to imp_ctl, needs reset
          END DO
        END DO
!$OMP END DO NOWAIT
      END IF

      DO k = 1, nice_use
!$OMP DO SCHEDULE(STATIC)
        DO j = 1, rows
          DO i = 1, row_length
            sice_save(i,j,k,1)=radnet_sice(i,j,k)!INOUT to imp_ctl, needs reset
            sice_save(i,j,k,2)=fqw_ice(i,j,k)    !INOUT to imp_ctl, needs reset
            sice_save(i,j,k,3)=ftl_ice(i,j,k)    !INOUT to imp_ctl, needs reset
            sice_save(i,j,k,4)=dtstar_sice(i,j,k)!INOUT to imp_ctl, needs reset
          END DO
        END DO
!$OMP END DO NOWAIT
      END DO

!$OMP END PARALLEL

!$OMP PARALLEL DO IF(ntiles > 1) DEFAULT(NONE) SCHEDULE(STATIC)        &
!$OMP SHARED ( ntiles, land_points, tile_save, ftl_tile, fqw_tile,     &
!$OMP          epot_tile, dtstar_tile )                                &
!$OMP PRIVATE( i, j )
      DO j = 1, ntiles
        DO i = 1, land_points
          tile_save(i,j,1)=ftl_tile(i,j)     !INOUT to imp_ctl, needs reset
          tile_save(i,j,2)=fqw_tile(i,j)     !INOUT to imp_ctl, needs reset
          tile_save(i,j,3)=epot_tile(i,j)    !INOUT to imp_ctl, needs reset
          tile_save(i,j,4)=dtstar_tile(i,j)  !INOUT to imp_ctl, needs reset
        END DO
      END DO
!$OMP END PARALLEL DO

      land_save(:,1)=gs !INOUT to bl_ctl, gets reset in atm_step

    END IF ! l_quick_ap2

  ELSE ! ( cycleno == 1 .OR. .NOT. l_quick_ap2 )

    !     restore outputs on second EG cycle
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,i,k)

!$OMP DO SCHEDULE(STATIC)
    DO j = 1, rows
      DO i = 1, row_length
        zh(i,j)     = bl_ctl_2d(i,j,1)
        z0msea(i,j) = bl_ctl_2d(i,j,2)
        wstar(i,j)  = bl_ctl_2d(i,j,3)
        wthvs(i,j)  = bl_ctl_2d(i,j,4)

        ntml(i,j) = bl_ctl_int2d(i,j,1)

        cumulus(i,j)   = bl_ctl_log2d(i,j,1)
        l_shallow(i,j) = bl_ctl_log2d(i,j,2)
      END DO
    END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
    DO k = 1, bl_levels
      DO j = 1, rows
        DO i = 1, row_length
          fqw(i,j,k)   = bl_ctl_3d(i,j,k,1)
          ftl(i,j,k)   = bl_ctl_3d(i,j,k,2)
          rhokh(i,j,k) = bl_ctl_3d(i,j,k,3)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

    IF (l_use_dtstar_sea) THEN
!$OMP DO SCHEDULE(STATIC)
      DO j = 1, rows
        DO i = 1, row_length
          radnet_sea(i,j) = sea_save(i,j,1)
          dtstar_sea(i,j) = sea_save(i,j,2)
        END DO
      END DO
!$OMP END DO NOWAIT
    END IF

    DO k = 1, nice_use
!$OMP DO SCHEDULE(STATIC)
      DO j = 1, rows
        DO i = 1, row_length
          radnet_sice(i,j,k) = sice_save(i,j,k,1)
          fqw_ice(i,j,k)     = sice_save(i,j,k,2)
          ftl_ice(i,j,k)     = sice_save(i,j,k,3)
          dtstar_sice(i,j,k) = sice_save(i,j,k,4)
        END DO
      END DO
!$OMP END DO NOWAIT
    END DO

!$OMP END PARALLEL

!$OMP PARALLEL DO IF(ntiles > 1) DEFAULT(NONE) SCHEDULE(STATIC)        &
!$OMP SHARED ( ntiles, land_points, tile_save, ftl_tile, fqw_tile,     &
!$OMP          epot_tile, dtstar_tile )                                &
!$OMP PRIVATE( i, j )
    DO j = 1, ntiles
      DO i = 1, land_points
        ftl_tile(i,j)    = tile_save(i,j,1)
        fqw_tile(i,j)    = tile_save(i,j,2)
        epot_tile(i,j)   = tile_save(i,j,3)
        dtstar_tile(i,j) = tile_save(i,j,4)
      END DO
    END DO
!$OMP END PARALLEL DO

    gs = land_save(:,1)

  END IF ! ( cycleno == 1 .OR. .NOT. l_quick_ap2 )

END IF ! i_bl_vn

IF (Ltimer) CALL timer ('AP2 Explicit BL',6)

IF (Ltimer) CALL timer ('AP2 Boundary Layer',6)


!------------------------------------------------------------------------
! Calculate Leonard term fluxes and increments
!------------------------------------------------------------------------
IF ( l_leonard_term .AND. i_bl_vn /= i_bl_vn_0 ) THEN
  ! Note: the Leonard terms code will not work if ni_bl_ctl isn't called,
  ! due to dependence on the dt/(dz rho r^2) arrays 
  ! dtrdz_charney_grid, dtrdz_u, dtrdz_v,
  ! which are calculated in the explicit BL scheme code.

  IF ( cycleno == 1 .OR. .NOT. l_quick_ap2 ) THEN

    CALL leonard_term_ctl( etadot, u, v, w, theta, q, qcl, qcf,         &
                           rho_wet_rsq, exner_theta_levels,             &
                           exner_rho_levels,                            &
                           dtrdz_charney_grid, dtrdz_u, dtrdz_v,        &
                           u_inc_leonard, v_inc_leonard, w_inc_leonard, &
                           thetal_inc_leonard, qw_inc_leonard,          &
                           stashwork3 )

  END IF

  ! Add Leonard term increments to the latest fields

  ! Increments to 3-D winds:
  DO k = 1, bl_levels
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        r_u(i,j,k)        = r_u(i,j,k)        + u_inc_leonard(i,j,k)
      END DO
    END DO
  END DO
  DO k = 1, bl_levels
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        r_v(i,j,k)        = r_v(i,j,k)        + v_inc_leonard(i,j,k)
      END DO
    END DO
  END DO
  DO k = 1, bl_levels
    DO j = wdims%j_start, wdims%j_end
      DO i = wdims%i_start, wdims%i_end
        r_w(i,j,k)        = r_w(i,j,k)        + w_inc_leonard(i,j,k)
      END DO
    END DO
  END DO

  ! Add increments of theta_l / qw onto theta and q.  The appropriate
  ! cloud response to the Leonard terms is done in ni_imp_ctl.
  DO k = 1, bl_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        theta_star(i,j,k) = theta_star(i,j,k) + thetal_inc_leonard(i,j,k)
      END DO
    END DO
  END DO
  DO k = 1, bl_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        q_star(i,j,k)     = q_star(i,j,k)     + qw_inc_leonard(i,j,k)
      END DO
    END DO
  END DO

END IF  ! ( l_leonard_term .AND. i_bl_vn /= i_bl_vn_0 )


!------------------------------------------------------------------------
! Which cumulus points has BL changed?
!  cumulus_copy - cumulus as diagnosed by conv_diag
!  cumulus      - cumulus array after BL
!  no_cumulus   - .true. for those points which have been changed from
!                  .true. to .false.
! Note the BL scheme has good reasons for changing the diagnosis of
! cumulus to .false. for a small number of locations. The information
! used to make these changes is not available to the convection scheme.
!------------------------------------------------------------------------
!$OMP PARALLEL DEFAULT(NONE)                                           &
!$OMP SHARED( rows, row_length, no_cumulus, cumulus, i_convection_vn,  &
!$OMP         cumulus_copy, l_pc2_diag_sh, l_pc2_diag_sh_pts,          &
!$OMP         l_shallow, bl_type_6 )                                   &
!$OMP PRIVATE( i, j )
!$OMP DO SCHEDULE(STATIC)
DO j= 1,rows
  DO i=1,row_length
    no_cumulus(i,j) = .NOT. cumulus(i,j) .AND. cumulus_copy(i,j)
  END DO
END DO
!$OMP END DO NOWAIT

! PC2 required mask
IF (i_convection_vn == i_convection_vn_5a .OR.                          &
    i_convection_vn == i_convection_vn_6a ) THEN
      ! Calculate for each point whether we want to carry
      ! convective cloud information for shallow conv in PC2
  IF ( l_pc2_diag_sh) THEN
    ! Carry convective cloud information
!$OMP DO SCHEDULE(STATIC)
    DO j = 1, rows
      DO i = 1, row_length
        l_pc2_diag_sh_pts(i,j) = l_shallow(i,j)
      END DO
    END DO
!$OMP END DO NOWAIT
  ELSE
    ! Do not carry convective cloud information
!$OMP DO SCHEDULE(STATIC)
    DO j = 1, rows
      DO i = 1, row_length
        l_pc2_diag_sh_pts(i,j) = .FALSE.
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF ! l_pc2_diag_sh

ELSE ! i_convection_vn
      ! Calculate for each point whether we want to carry
      ! convective cloud information for shallow conv in PC2
!$OMP DO SCHEDULE(STATIC)
  DO j = 1, rows
    DO i = 1, row_length
      IF (l_pc2_diag_sh .AND. l_shallow(i,j) .AND.                      &
           bl_type_6(i,j) == 1.0) THEN
        ! Carry convective cloud information
        l_pc2_diag_sh_pts(i,j) = .TRUE.
      ELSE
        ! Do not carry convective cloud information
        l_pc2_diag_sh_pts(i,j) = .FALSE.
      END IF
    END DO
  END DO
!$OMP END DO NOWAIT
END IF ! i_convection_vn
!$OMP END PARALLEL

! ----------------------------------------------------------------------
! Section CNV.1 CALL Convection scheme.
! ----------------------------------------------------------------------
l_CallConvection = .FALSE.
SELECT CASE (model_type)

CASE (mt_single_column)
  IF ( l_param_conv .AND. (conv_mode < 2) ) l_CallConvection = .TRUE.

CASE DEFAULT
  IF ( l_param_conv ) l_CallConvection = .TRUE.

END SELECT

IF (l_CallConvection) THEN
  IF ( Error_code == 0 ) THEN

    IF (Ltimer) CALL timer ('AP2 Convection',5)

    !     Set up logicals for condensate increment calculation.
    !     The PC2 Convection code is implemented in three phases:
    !      1. OFF-OFF No extra diagnostic space reserved, no calculations.
    !      2.  ON-OFF NEW diagnostics calculated, no existing fields touched.
    !      3.  ON-ON  Full prognostic interactions and accompanying diagnostics.
    !
    !  Other parts of the PC2 code use i_cld_PC2 to perform overwriting calcs,
    !  then L_PC2_RESET to write back the original fields for non-interacting
    !  mode if required. Hence use of different logicals here.
    !
    L_calc_dxek  = ( i_cld_vn == i_cld_pc2 )
    L_q_interact = ( i_cld_vn == i_cld_pc2  .AND.  ( .NOT. L_pc2_reset ) )

    !-----------------------------------------------------------------------
    ! Setup flags to control convection output diagnostics
    !-----------------------------------------------------------------------
    !       Apply diags at last cycle only
    CALL set_convection_output_flags()

    !     Convection diagnostics needed for SKEB2
    flg_up_flx       = .TRUE.
    flg_dwn_flx      = .TRUE.

    ! Allocate SPT arrays
    IF (l_retain_conv_all_tendencies) THEN
      CALL init_convection_tendencies
    END IF

    !---------------------------------------------------------------------------
    !     Allocate arrays required for convective diagnostics
    !---------------------------------------------------------------------------

    CALL cv_alloc_diag_array( row_length, rows, l_cosp,                 &
                              l_dust, l_sulpc_so2, l_sulpc_nh3, l_soot, &
                              l_biomass, l_ocff, l_nitrate )

    IF (i_convection_vn < 10) THEN

      !----------------------------------------------------------------------
      ! Met Office supported mass flux convection scheme
      ! Based on Gregory-Rowntree
      ! Current allowed variants 5A and 6A schemes
      !----------------------------------------------------------------------

      CALL NI_conv_ctl (                                                &
      !     Parallel variables
            at_extremity, n_proc                                        &
          , delta_lambda,delta_phi                                      &
      ! parameters for cycling physics-dynamics
          , NumCycles, CycleNo                                          &
      !     Model dimensions.
          , row_length, rows                                            &
          , rows*row_length                                             &
          , bl_levels, n_cca_levels                                     &
          , tr_vars, tr_ukca                                            &

      !     Model switches
          , L_calc_dxek, L_q_interact, l_spec_z0                        &

      !     IN coordinate information
          ,z_rho, z_theta                                               &
      !     SCM/Idealised UM
          , flux_e,flux_h,ustar_in,z0m_scm,z0h_scm                      &
          , tstar, zh, u_0_p, v_0_p, fb_surf, uStarGBM                  &
          , ls_rain, ls_snow                                            &
      !
      !     SCM diagnostics (dummy in full UM)
          , nSCMDpkgs,L_SCMDiags, conv_mode                             &

      !     IN data fields.
          , rho_wet_rsq, rho_wet, rho_wet_tq, rho_dry, rho_dry_theta    &
          , u_p,v_p, ustar_p, vstar_p,w_copy, p, p_star, exner_rho_levels&
          , land_sea_mask, flandg, ice_fract                            &
          , tstar_land, tstar_sea, tstar_sice, z0msea                   &
          , p_layer_boundaries, p_layer_centres                         &
          , exner_layer_boundaries, exner_layer_centres                 &
          , t1_sd, q1_sd, exner_theta_levels                            &
          , uw0, vw0, w_max, zlcl, zlcl_uv, zhpar, dzh, qcl_inv_top     &
          , entrain_coef, conv_type                                     &
          , cumulus, l_shallow, l_congestus, l_congestus2, l_mid_level  &
          , l_pc2_diag_sh_pts                                           &
          , no_cumulus                                                  &
          , ntml, ntpar                                                 &
          , wstar, wthvs, delthvu, ql_ad, qsat_lcl, ftl ,fqw            &
          , shallowc, cu_over_orog, cape_undilute, cin_undilute         &
          , deep_flag, past_precip, past_conv_ht                        &

      !     IN start of time step prognostic values
          , theta, q, qcl, qcf, qrain, qgraup, qcf2                     &
          , cloud_fraction_liquid_halos                                 &
          , cloud_fraction_frozen_halos, bulk_cloud_fraction_halos      &
      !     IN cold-pool fields
          , g_ccp, h_ccp                                                &
      !     IN/OUT primary fields with all increments so far added on
          , theta_star, q_star, qcl_star, qcf_star                      &
          , qcf2_star, qrain_star, qgraup_star                          &
          , cfl_star, cff_star, cf_star                                 &
      !     OUT convective momentum transport tendencies on p-grid
          , dubydt_conv_p, dvbydt_conv_p                                &
      !     IN/OUT values with all increments added
          , conv_prog_1, conv_prog_2, conv_prog_3, conv_prog_precip     &
      !     IN total wind increments before convection, on p-grid
          , r_u_p, r_v_p                                                &
      !     IN/OUT tracers
          , aerosol                                                     &
          , dust_div1,dust_div2,dust_div3,dust_div4,dust_div5,dust_div6 &
          , so2, so4_aitken, so4_accu, so4_diss                         &
          , dms, nh3, soot_new, soot_aged, soot_cld, bmass_new, bmass_aged&
          , bmass_cld, ocff_new, ocff_aged, ocff_cld, nitr_acc, nitr_diss&
          , co2,  free_tracers, tracer_ukca                             &
          , ozone_tracer                                                &

      ! out fields
          , cca0_dp, cca0_md, cca0_sh                                   &
          , cca0, ccw0, ccb0, cct0, cclwp0, lcbase0, cca0_2d, lctop, lcca&
          , cca , ccw , ccb , cct , cclwp , lcbase,  cca_2d             &
          , conv_rain, conv_snow, ddmfx                                 &

      !     ERROR information
          , Error_code  )

    ELSE
      !----------------------------------------------------------------------
      ! Top level control to call alternative convection schemes
      ! Must have i_convection_vn > 9
      ! None at present - expect Betts-Miller to be added here
      !----------------------------------------------------------------------
      CALL other_conv_ctl(                                                   &
                numcycles, cycleno,                                          &
                row_length, rows, n_cca_levels,                              &
                tr_levels, tr_vars, tr_ukca,                                 &
                land_points, nscmdpkgs,                                      &
                ntml,                                                        &
             ! Logical in
                l_calc_dxek, l_q_interact, l_scmdiags,                       &
                land_sea_mask,                                               &
             ! Real in...
             ! Turbulence fields
                zh, dzh, taux_p, tauy_p, ftl, fqw, rhokm, rhokh, bl_w_var,   &
             ! Density, pressure, exner
                rho_wet_rsq, rho_wet, rho_wet_tq, rho_dry, rho_dry_theta,    &
                p_star, p, exner_rho_levels, exner_theta_levels,             &
                p_layer_boundaries, p_layer_centres,                         &
                exner_layer_boundaries, exner_layer_centres,                 &
             ! in Prognostics at beginning of time step
                theta, q, qcl, qcf, qrain, qgraup, qcf2,                     &
                cloud_fraction_liquid_halos, cloud_fraction_frozen_halos,    &
                bulk_cloud_fraction_halos,                                   &
                u_p, v_p, w,                                                 &
                m_v, m_cl, m_cf, m_cf2, m_r, m_gr,                           &
             ! Integer inout (cloud level info)
                ccb0, cct0, lcbase0, ccb, cct, lcbase, lctop,                &
             ! In/Out: fields updated with all increments so far
             ! For safety, the k=0 level (which should not be altered by
             ! convection) is not passed in.
                theta_star(:,:,1:tdims%k_end), q_star(:,:,1:tdims%k_end),    &
                qcl_star(:,:,1:tdims%k_end), qcf_star(:,:,1:tdims%k_end),    &
                qcf2_star(:,:,1:tdims%k_end), qrain_star(:,:,1:tdims%k_end), &
                qgraup_star(:,:,1:tdims%k_end),                              &
                cfl_star(:,:,1:tdims%k_end), cff_star(:,:,1:tdims%k_end),    &
                cf_star(:,:,1:tdims%k_end),                                  &
             ! In: winds updated with increments so far, interpolated to p-grid
                ustar_p, vstar_p,                                            &
             ! Out: convective momentum transport tendencies on p-grid
                dubydt_conv_p, dvbydt_conv_p,                                &
             ! Cloud prognostics & diagostics
                cca0_dp, cca0_md, cca0_sh, cca0, ccw0, cclwp0, cca0_2d,      &
                lcca, cca,  ccw, cclwp, cca_2d,                              &
             ! Other output fields
                conv_rain, conv_snow,                                        &
             ! Tracers
                aerosol, dust_div1, dust_div2,                               &
                dust_div3, dust_div4, dust_div5, dust_div6,                  &
                so2, so4_aitken, so4_accu, so4_diss,                         &
                dms, nh3, soot_new, soot_aged, soot_cld,                     &
                bmass_new, bmass_aged, bmass_cld,                            &
                ocff_new, ocff_aged, ocff_cld, nitr_acc, nitr_diss,          &
                co2,  free_tracers, tracer_ukca, ozone_tracer,               &
             ! Real out - may be expected by surface scheme
                ddmfx                                                        &
                )

    END IF   ! on convection scheme control

    !     Arrays from Convection required by SKEB2
    !     Passed through module stochastic_physics_run_mod
    !     Deallocation occurs in stph_skeb2
    IF (.NOT. ALLOCATED(skeb2_up_flux)) THEN
      ALLOCATE (skeb2_up_flux(row_length, rows, model_levels))
    END IF
    IF (.NOT. ALLOCATED(skeb2_dwn_flux)) THEN
      ALLOCATE (skeb2_dwn_flux(row_length, rows, model_levels))
    END IF
    IF (.NOT. ALLOCATED(skeb2_cape)) THEN
      ALLOCATE (skeb2_cape(row_length, rows))
    END IF

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(k,j,i)                          &
!$OMP SHARED(model_levels,rows,row_length,skeb2_up_flux,up_flux,     &
!$OMP    skeb2_dwn_flux,dwn_flux,skeb2_cape,cape_out,l_cosp,         &
!$OMP    cosp_crain_3d, cosp_csnow_3d, conv_rain_3d, conv_snow_3d)

!$OMP DO SCHEDULE(STATIC)
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          skeb2_up_flux(i,j,k) = up_flux(i,j,k)
          skeb2_dwn_flux(i,j,k) = dwn_flux(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
    DO j = 1, rows
      DO i = 1, row_length
        skeb2_cape(i,j) = cape_out(i,j)
      END DO
    END DO
!$OMP END DO NOWAIT

    !     Copy variables required by COSP
    IF ( L_cosp ) THEN
!$OMP DO SCHEDULE(STATIC)
      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
            cosp_crain_3d(i,j,k) = conv_rain_3d(i,j,k)
            cosp_csnow_3d(i,j,k) = conv_snow_3d(i,j,k)
            IF (cosp_crain_3d(i,j,k) < 0.0) cosp_crain_3d(i,j,k) = 0.0
            IF (cosp_csnow_3d(i,j,k) < 0.0) cosp_csnow_3d(i,j,k) = 0.0
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT
    END IF ! L_cosp

!$OMP END PARALLEL

    IF (model_type == mt_single_column) THEN

      ! DEPENDS ON: sub_data
      CALL sub_data(                                                    &
           row_length, rows,                                            &
           bl_levels, st_levels, sm_levels,                             &
           ' After convect, before hydr        ', a_step,               &
           val_year, val_day_number,                                    &
           u, v, theta, q, qcl, qcf, bulk_cloud_fraction_halos,         &
           p, rho_wet_rsq, exner_rho_levels, exner_theta_levels,        &
           t_soil, smc,canopy, snow_depth, tstar, zh, z0msea,           &
           cca, ccb, cct, soil_layer_moisture )

    END IF ! model_type

    IF ( l_mom ) THEN

      ! Deallocate temporary copies of r_u and u_star on p-grid
      ! that were needed for the CMT calc in the convection scheme
      DEALLOCATE( vstar_p )
      DEALLOCATE( ustar_p )
      DEALLOCATE( r_v_p )
      DEALLOCATE( r_u_p )

    END IF ! l_mom

    IF (Ltimer) CALL timer ('AP2 Convection',6)

  END IF ! on error code equal to zero


  IF ( Error_code  ==  0 ) THEN  ! Check error code still zero

    IF (Ltimer) CALL timer ('AP2 Convection',5)

    IF (Ltimer) CALL timer ('AP2 Convection',6)

    ! ----------------------------------------------------------------------
    ! Section CNV.2 Energy correction code
    ! ----------------------------------------------------------------------

    IF ( CycleNo == NumCycles .AND. L_emcorr ) THEN

      IF (Ltimer) CALL timer ('AP2 Conv Eng Corr',5)

      !     Add convective rain and snow, at the surface to the
      !     diabatic heating for use in the energy correction
      !     procedure.
      !     Scale variables by conversion factor so that only one call
      !     is required

      lclf = lc + lf

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, j)                           &
!$OMP SHARED(tot_precip_scaled_1,tot_precip_scaled_2, conv_rain,     &
!$OMP conv_snow,lclf,rows,row_length)

!$OMP DO SCHEDULE(STATIC)
      DO j = 1, rows
        DO i = 1, row_length
          tot_precip_scaled_1(i,j) =  conv_rain(i,j) * lc +     &
                                      conv_snow(i,j) * lclf
        END DO
      END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
      DO j = 1, rows
        DO i = 1, row_length
          tot_precip_scaled_2(i,j) = -conv_rain(i,j)-conv_snow(i,j)
        END DO
      END DO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

      CALL flux_diag(tot_precip_scaled_1, xx_cos_theta_latitude,        &
                     row_length, rows ,offx,offy, 1.0,                  &
                     sum_eng_fluxes,    timestep)

      CALL flux_diag(tot_precip_scaled_2, xx_cos_theta_latitude,        &
                     row_length, rows ,offx,offy, 1.0,                  &
                     sum_moist_flux,    timestep)

      IF (Ltimer) CALL timer ('AP2 Conv Eng Corr',6)

    END IF   ! L_emcorr

  END IF ! Error_code

END IF ! l_param_conv convection option

! ----------------------------------------------------------------------
! Section BL CALL implicit solver
! ---------------------------------------------------------------------
IF (Ltimer) CALL timer ('AP2 Boundary Layer',5)

IF (Ltimer) CALL timer ('AP2 Implicit BL',5)

IF ( Error_code  ==  0 ) THEN
  !--------------------------------------------------------------------------
  ! Second communication section before ni_imp_ctl
  !--------------------------------------------------------------------------
  IF (model_type /= mt_single_column) THEN

    !     swap bounds for 3d fields first
    i_field = 0

    IF ( i_bl_vn /= i_bl_vn_0 ) THEN
      !     only call on 1st cycle or if not fast running
      IF ( cycleno == 1 .OR. .NOT. l_quick_ap2 ) THEN

        i_field = i_field + 1
        fields_to_swap(i_field) % field       => rhokm(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_p
        fields_to_swap(i_field) % levels      =  bl_levels
        fields_to_swap(i_field) % rows        =  rows
        fields_to_swap(i_field) % vector      =  .FALSE.

        IF (formdrag ==  explicit_stress) THEN

          i_field = i_field + 1
          fields_to_swap(i_field) % field       => tau_fd_x(:,:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  bl_levels
          fields_to_swap(i_field) % rows        =  rows
          fields_to_swap(i_field) % vector      =  .FALSE.

          i_field = i_field + 1
          fields_to_swap(i_field) % field       => tau_fd_y(:,:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  bl_levels
          fields_to_swap(i_field) % rows        =  rows
          fields_to_swap(i_field) % vector      =  .FALSE.

        END IF

        IF ( i_bl_vn == i_bl_vn_1a ) THEN
          i_field = i_field + 1
          fields_to_swap(i_field) % field       => rhogamu(:,:,1:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  bl_levels-1
          fields_to_swap(i_field) % rows        =  rows
          fields_to_swap(i_field) % vector      =  .FALSE.

          i_field = i_field + 1
          fields_to_swap(i_field) % field       => rhogamv(:,:,1:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  bl_levels-1
          fields_to_swap(i_field) % rows        =  rows
          fields_to_swap(i_field) % vector      =  .FALSE.
        ELSE
          i_field = i_field + 1
          fields_to_swap(i_field) % field       => f_ngstress(:,:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  bl_levels-1
          fields_to_swap(i_field) % rows        =  rows
          fields_to_swap(i_field) % vector      =  .FALSE.
        END IF ! i_bl_vn

      END IF ! ( cycleno == 1 .OR. .NOT. l_quick_ap2 )

    END IF ! i_bl_vn

    !     Convection winds fields requiring swap bound if CMT being used

    IF ( l_param_conv .AND. l_mom ) THEN

      i_field = i_field + 1
      fields_to_swap(i_field) % field       => dubydt_conv_p(:,:,:)
      fields_to_swap(i_field) % field_type  =  fld_type_p
      fields_to_swap(i_field) % levels      =  model_levels
      fields_to_swap(i_field) % rows        =  rows
      fields_to_swap(i_field) % vector      =  .FALSE.

      i_field = i_field + 1
      fields_to_swap(i_field) % field       => dvbydt_conv_p(:,:,:)
      fields_to_swap(i_field) % field_type  =  fld_type_p
      fields_to_swap(i_field) % levels      =  model_levels
      fields_to_swap(i_field) % rows        =  rows
      fields_to_swap(i_field) % vector      =  .FALSE.

    END IF ! l_param_conv .AND. l_mom

    CALL swap_bounds_mv(fields_to_swap, i_field, row_length,            &
                        offx , offy )
  END IF ! model_type

  IF ( i_bl_vn /= i_bl_vn_0 ) THEN
    !   only call on 1st cycle or if not fast running
    IF ( cycleno == 1 .OR. .NOT. l_quick_ap2 ) THEN


      SELECT CASE (model_type)
      CASE DEFAULT
        !         then 2d fields
        i_field = 0
        i_field = i_field + 1
        fields_to_swap(i_field) % field_2d   => rhokm_land(:,:)
        fields_to_swap(i_field) % field_type = fld_type_p
        fields_to_swap(i_field) % levels     = 1
        fields_to_swap(i_field) % rows       = rows
        fields_to_swap(i_field) % vector     = .FALSE.

        i_field = i_field + 1
        fields_to_swap(i_field) % field_2d   => rhokm_ssi(:,:)
        fields_to_swap(i_field) % field_type = fld_type_p
        fields_to_swap(i_field) % levels     = 1
        fields_to_swap(i_field) % rows       = rows
        fields_to_swap(i_field) % vector     = .FALSE.

        IF (l_ctile .AND. buddy_sea == on) THEN
          !           Interpolate wind speed factors to u and v columns

          i_field = i_field + 1
          fields_to_swap(i_field) % field_2d   => flandfac(:,:)
          fields_to_swap(i_field) % field_type = fld_type_p
          fields_to_swap(i_field) % levels     = 1
          fields_to_swap(i_field) % rows       = rows
          fields_to_swap(i_field) % vector     = .FALSE.

          i_field = i_field + 1
          fields_to_swap(i_field) % field_2d   => fseafac(:,:)
          fields_to_swap(i_field) % field_type = fld_type_p
          fields_to_swap(i_field) % levels     = 1
          fields_to_swap(i_field) % rows       = rows
          fields_to_swap(i_field) % vector     = .FALSE.

        END IF !c_tile and buddy_sea

        IF (sf_diag%su10 .OR. sf_diag%sv10) THEN
          i_field = i_field + 1
          fields_to_swap(i_field) % field_2d   => cdr10m(:,:)
          fields_to_swap(i_field) % field_type = fld_type_p
          fields_to_swap(i_field) % levels     = 1
          fields_to_swap(i_field) % rows       = rows
          fields_to_swap(i_field) % vector     = .FALSE.
        END IF ! (su10 .OR. sv10)

        IF (sf_diag%suv10m_n) THEN

          i_field = i_field + 1
          fields_to_swap(i_field) % field_2d   => cdr10m_n(:,:)
          fields_to_swap(i_field) % field_type = fld_type_p
          fields_to_swap(i_field) % levels     = 1
          fields_to_swap(i_field) % rows       = rows
          fields_to_swap(i_field) % vector     = .FALSE.

          i_field = i_field + 1
          fields_to_swap(i_field) % field_2d   => cd10m_n(:,:)
          fields_to_swap(i_field) % field_type = fld_type_p
          fields_to_swap(i_field) % levels     = 1
          fields_to_swap(i_field) % rows       = rows
          fields_to_swap(i_field) % vector     = .FALSE.

        END IF ! suv10m_n


        CALL swap_bounds_2d_mv(fields_to_swap, i_field, row_length,     &
                               offx , offy )

!$OMP PARALLEL DEFAULT(NONE)                                           &
!$OMP SHARED(rhokm, rhokm_u, rhokm_v, pdims_s, udims, vdims, bl_levels)
        CALL p_to_u(rhokm,                                              &
                     pdims_s%i_start,pdims_s%i_end,                     &
                     pdims_s%j_start,pdims_s%j_end,                     &
                     udims%i_start,udims%i_end,                         &
                     udims%j_start,udims%j_end,                         &
                     0,bl_levels-1, rhokm_u)

        CALL p_to_v(rhokm,                                              &
                     pdims_s%i_start,pdims_s%i_end,                     &
                     pdims_s%j_start,pdims_s%j_end,                     &
                     vdims%i_start,vdims%i_end,                         &
                     vdims%j_start,vdims%j_end,                         &
                     0,bl_levels-1, rhokm_v)
!$OMP END PARALLEL

        CALL p_to_u_land(rhokm_land, flandg,                            &
                     pdims_s%i_start,pdims_s%i_end,                     &
                     pdims_s%j_start,pdims_s%j_end,                     &
                     udims%i_start,udims%i_end,                         &
                     udims%j_start,udims%j_end,                         &
                     1,1, rhokm_u_land)

        CALL p_to_u_sea(rhokm_ssi, flandg,                              &
                     pdims_s%i_start,pdims_s%i_end,                     &
                     pdims_s%j_start,pdims_s%j_end,                     &
                     udims%i_start,udims%i_end,                         &
                     udims%j_start,udims%j_end,                         &
                     1,1, rhokm_u_ssi)

        CALL p_to_u(flandg,                                             &
                     pdims_s%i_start,pdims_s%i_end,                     &
                     pdims_s%j_start,pdims_s%j_end,                     &
                     udims%i_start,udims%i_end,                         &
                     udims%j_start,udims%j_end,                         &
                     1,1,flandg_u)

        CALL p_to_v_land(rhokm_land, flandg,                            &
                     pdims_s%i_start,pdims_s%i_end,                     &
                     pdims_s%j_start,pdims_s%j_end,                     &
                     vdims%i_start,vdims%i_end,                         &
                     vdims%j_start,vdims%j_end,                         &
                     1,1, rhokm_v_land)

        CALL p_to_v_sea(rhokm_ssi, flandg,                              &
                     pdims_s%i_start,pdims_s%i_end,                     &
                     pdims_s%j_start,pdims_s%j_end,                     &
                     vdims%i_start,vdims%i_end,                         &
                     vdims%j_start,vdims%j_end,                         &
                     1,1, rhokm_v_ssi)

        CALL p_to_v(flandg,                                             &
                     pdims_s%i_start,pdims_s%i_end,                     &
                     pdims_s%j_start,pdims_s%j_end,                     &
                     vdims%i_start,vdims%i_end,                         &
                     vdims%j_start,vdims%j_end,                         &
                     1,1,flandg_v)

        IF (l_ctile .AND. buddy_sea == on) THEN

          CALL p_to_u(flandfac,                                         &
                    pdims_s%i_start,pdims_s%i_end,                      &
                    pdims_s%j_start,pdims_s%j_end,                      &
                    udims%i_start,udims%i_end,                          &
                    udims%j_start,udims%j_end,                          &
                    1,1, flandfac_u)

          CALL p_to_u(fseafac,                                          &
                    pdims_s%i_start,pdims_s%i_end,                      &
                    pdims_s%j_start,pdims_s%j_end,                      &
                    udims%i_start,udims%i_end,                          &
                    udims%j_start,udims%j_end,                          &
                    1,1, fseafac_u)

          CALL p_to_v(flandfac,                                         &
                    pdims_s%i_start,pdims_s%i_end,                      &
                    pdims_s%j_start,pdims_s%j_end,                      &
                    vdims%i_start,vdims%i_end,                          &
                    vdims%j_start,vdims%j_end,                          &
                    1,1, flandfac_v)

          CALL p_to_v(fseafac,                                          &
                    pdims_s%i_start,pdims_s%i_end,                      &
                    pdims_s%j_start,pdims_s%j_end,                      &
                    vdims%i_start,vdims%i_end,                          &
                    vdims%j_start,vdims%j_end,                          &
                    1,1, fseafac_v)

        END IF !c_tile and buddy_sea

      CASE (mt_single_column)

        DO k = 0, bl_levels-1
          DO j = 1, rows
            DO i = 1, row_length
              rhokm_u(i,j,k) = rhokm(i,j,k)
              rhokm_v(i,j,k) = rhokm(i,j,k)
            END DO
          END DO
        END DO

        DO j = 1, rows
          DO i = 1, row_length
            rhokm_u_land(i,j) = rhokm_land(i,j)
            rhokm_v_land(i,j) = rhokm_land(i,j)
            rhokm_u_ssi(i,j)  = rhokm_ssi(i,j)
            rhokm_v_ssi(i,j)  = rhokm_ssi(i,j)
            flandg_u(i,j)     = flandg(i,j)
            flandg_v(i,j)     = flandg(i,j)
            flandfac_u(i,j)   = 1.0
            flandfac_v(i,j)   = 1.0
            fseafac_u(i,j)    = 1.0
            fseafac_v(i,j)    = 1.0
          END DO
        END DO

      END SELECT ! model_type

      SELECT CASE (model_type)

      CASE DEFAULT
        IF (sf_diag%su10) THEN
          CALL p_to_u(cdr10m,                                           &
                  pdims_s%i_start,pdims_s%i_end,                        &
                  pdims_s%j_start,pdims_s%j_end,                        &
                  udims%i_start,udims%i_end,                            &
                  udims%j_start,udims%j_end,                            &
                  1,1, cdr10m_u)

        END IF ! su10

        IF (sf_diag%sv10) THEN
          CALL p_to_v(cdr10m,                                           &
                  pdims_s%i_start,pdims_s%i_end,                        &
                  pdims_s%j_start,pdims_s%j_end,                        &
                  vdims%i_start,vdims%i_end,                            &
                  vdims%j_start,vdims%j_end,                            &
                  1,1, cdr10m_v)

        END IF ! sv10

        IF (sf_diag%suv10m_n) THEN

          CALL p_to_u(cdr10m_n,                                         &
                    pdims_s%i_start,pdims_s%i_end,                      &
                    pdims_s%j_start,pdims_s%j_end,                      &
                    udims%i_start,udims%i_end,                          &
                    udims%j_start,udims%j_end,                          &
                    1,1, cdr10m_n_u)

          CALL p_to_u(cd10m_n,                                          &
                    pdims_s%i_start,pdims_s%i_end,                      &
                    pdims_s%j_start,pdims_s%j_end,                      &
                    udims%i_start,udims%i_end,                          &
                    udims%j_start,udims%j_end,                          &
                    1,1, cd10m_n_u)

          CALL p_to_v(cdr10m_n,                                         &
                    pdims_s%i_start,pdims_s%i_end,                      &
                    pdims_s%j_start,pdims_s%j_end,                      &
                    vdims%i_start,vdims%i_end,                          &
                    vdims%j_start,vdims%j_end,                          &
                    1,1, cdr10m_n_v)

          CALL p_to_v(cd10m_n,                                          &
                    pdims_s%i_start,pdims_s%i_end,                      &
                    pdims_s%j_start,pdims_s%j_end,                      &
                    vdims%i_start,vdims%i_end,                          &
                    vdims%j_start,vdims%j_end,                          &
                    1,1, cd10m_n_v)

        END IF ! IF (suv10m_n)


      CASE (mt_single_column)

        DO j = 1, rows
          DO i = 1, row_length
            cdr10m_u(i,j) = cdr10m(i,j)
            cdr10m_v(i,j) = cdr10m(i,j)
          END DO
        END DO

        IF (sf_diag%suv10m_n) THEN
          DO j = 1, rows
            DO i = 1, row_length
              cdr10m_n_u(i,j) = cdr10m_n(i,j)
              cdr10m_n_v(i,j) = cdr10m_n(i,j)
              cd10m_n_u(i,j) = cd10m_n(i,j)
              cd10m_n_v(i,j) = cd10m_n(i,j)
            END DO
          END DO
        END IF

      END SELECT ! model_type

      IF ( i_bl_vn == i_bl_vn_1a ) THEN
        ALLOCATE(rhogamu_u(udims%i_start:udims%i_end,                   &
                           udims%j_start:udims%j_end,                   &
                           1:bl_levels-1))
        ALLOCATE(rhogamv_v(vdims%i_start:vdims%i_end,                   &
                           vdims%j_start:vdims%j_end,                   &
                           1:bl_levels-1))
        ALLOCATE(f_ngstress_u(1,1,1))
        ALLOCATE(f_ngstress_v(1,1,1))
      ELSE
        ALLOCATE(rhogamu_u(1,1,1))
        ALLOCATE(rhogamv_v(1,1,1))
        ALLOCATE(f_ngstress_u(udims%i_start:udims%i_end,                &
                              udims%j_start:udims%j_end,                &
                              1:bl_levels-1))
        ALLOCATE(f_ngstress_v(vdims%i_start:vdims%i_end,                &
                              vdims%j_start:vdims%j_end,                &
                              1:bl_levels-1))
      END IF ! i_bl_vn

      SELECT CASE (model_type)

      CASE DEFAULT
!$OMP PARALLEL DEFAULT(NONE)                                            &
!$OMP SHARED( formdrag, tau_fd_x, pdims_s, udims, bl_levels, taux_fd_u, &
!$OMP         tau_fd_y, vdims, tauy_fd_v, i_bl_vn, rhogamu, rhogamu_u,  &
!$OMP         rhogamv, rhogamv_v, f_ngstress, f_ngstress_u,             &
!$OMP         f_ngstress_v )
        IF ( formdrag == explicit_stress ) THEN
          CALL p_to_u (tau_fd_x,                                        &
                     pdims_s%i_start,pdims_s%i_end,                     &
                     pdims_s%j_start,pdims_s%j_end,                     &
                     udims%i_start,udims%i_end,                         &
                     udims%j_start,udims%j_end,                         &
                     0,bl_levels-1,taux_fd_u)
          CALL p_to_v (tau_fd_y,                                        &
                     pdims_s%i_start,pdims_s%i_end,                     &
                     pdims_s%j_start,pdims_s%j_end,                     &
                     vdims%i_start,vdims%i_end,                         &
                     vdims%j_start,vdims%j_end,                         &
                     0,bl_levels-1,tauy_fd_v)
        END IF ! formdrag

        IF ( i_bl_vn == i_bl_vn_1a ) THEN
          CALL p_to_u(rhogamu(:,:,1:bl_levels-1),                       &
                    pdims_s%i_start,pdims_s%i_end,                      &
                    pdims_s%j_start,pdims_s%j_end,                      &
                    udims%i_start,udims%i_end,                          &
                    udims%j_start,udims%j_end,                          &
                    1,bl_levels-1,rhogamu_u)

          CALL p_to_v(rhogamv(:,:,1:bl_levels-1),                       &
                    pdims_s%i_start,pdims_s%i_end,                      &
                    pdims_s%j_start,pdims_s%j_end,                      &
                    vdims%i_start,vdims%i_end,                          &
                    vdims%j_start,vdims%j_end,                          &
                    1,bl_levels-1,rhogamv_v)
        ELSE
          CALL p_to_u (f_ngstress,                                      &
                     pdims_s%i_start,pdims_s%i_end,                     &
                     pdims_s%j_start,pdims_s%j_end,                     &
                     udims%i_start,udims%i_end,                         &
                     udims%j_start,udims%j_end,                         &
                     1,bl_levels-1,f_ngstress_u)

          CALL p_to_v (f_ngstress,                                      &
                     pdims_s%i_start,pdims_s%i_end,                     &
                     pdims_s%j_start,pdims_s%j_end,                     &
                     vdims%i_start,vdims%i_end,                         &
                     vdims%j_start,vdims%j_end,                         &
                     1,bl_levels-1,f_ngstress_v)
        END IF ! i_bl_vn
!$OMP END PARALLEL


      CASE (mt_single_column)

        IF ( i_bl_vn == i_bl_vn_1a ) THEN
          DO k = 1, bl_levels-1
            DO j = 1, rows
              DO i = 1, row_length
                rhogamu_u(i,j,k) = rhogamu(i,j,k)
                rhogamv_v(i,j,k) = rhogamv(i,j,k)
              END DO
            END DO
          END DO
        ELSE
          DO k = 1, bl_levels-1
            DO j = 1, rows
              DO i = 1, row_length
                f_ngstress_v(i,j,k) = f_ngstress(i,j,k)
                f_ngstress_u(i,j,k) = f_ngstress(i,j,k)
              END DO
            END DO
          END DO
        END IF  ! i_bl_vn

        IF ( formdrag == explicit_stress ) THEN
          DO k = 0, bl_levels-1
            DO j = 1, rows
              DO i = 1, row_length
                taux_fd_u(i,j,k) = tau_fd_x(i,j,k)
                tauy_fd_v(i,j,k) = tau_fd_y(i,j,k)
              END DO
            END DO
          END DO
        END IF ! formdrag

      END SELECT ! model_type

      !       deallocate temp variables for message passing
      DEALLOCATE(rhogamu)
      DEALLOCATE(rhogamv)
      DEALLOCATE(f_ngstress)

    END IF ! ( cycleno == 1 .OR. .NOT. l_quick_ap2 )

  END IF ! i_bl_vn

  !-------------------------------------------------------------------------
  !   Rest of convection section regridding fields after swap_bounds and
  !   outputing convection diagnostics to stash
  !-------------------------------------------------------------------------

  IF (l_param_conv) THEN

    IF (Ltimer) CALL timer ('AP2 Convection',5)

    IF ( l_mom ) THEN
      IF ( .NOT. ( model_type==mt_single_column .AND. conv_mode==1 ) ) THEN
        ! If this is an SCM job running the convection scheme in
        ! diagnostics-only mode (conv_mode=1), don't add on the CMT increments.


        ALLOCATE(dubydt_conv_u(udims%i_start:udims%i_end,               &
                               udims%j_start:udims%j_end,               &
                               udims%k_start:udims%k_end) )
        ALLOCATE(dvbydt_conv_v(vdims%i_start:vdims%i_end,               &
                               vdims%j_start:vdims%j_end,               &
                               vdims%k_start:vdims%k_end) )

        ! CONV have set dubydt_conv_p, dvbydt_conv_p to be on tdims grid.
!$OMP PARALLEL DEFAULT(SHARED)                                          &
!$OMP SHARED( dubydt_conv_p, pdims_s, udims, dubydt_conv_u,             &
!$OMP         dvbydt_conv_p, vdims, dvbydt_conv_v,                      &
!$OMP         r_u, timestep, r_v, du_conv, dv_conv )                    &
!$OMP PRIVATE( i, j, k )
        ! interpolate to u grid
        CALL p_to_u (dubydt_conv_p,                                     &
                     pdims_s%i_start,pdims_s%i_end,                     &
                     pdims_s%j_start,pdims_s%j_end,                     &
                     udims%i_start,udims%i_end,                         &
                     udims%j_start,udims%j_end,                         &
                     1,udims%k_end, dubydt_conv_u)

        ! interpolate to v grid
        CALL p_to_v (dvbydt_conv_p,                                     &
                     pdims_s%i_start,pdims_s%i_end,                     &
                     pdims_s%j_start,pdims_s%j_end,                     &
                     vdims%i_start,vdims%i_end,                         &
                     vdims%j_start,vdims%j_end,                         &
                     1,vdims%k_end, dvbydt_conv_v)
!$OMP BARRIER

        ! u: add on to increment field

!$OMP DO SCHEDULE(STATIC)
        DO k =               1,udims%k_end
          DO j=  udims%j_start,udims%j_end
            DO i=udims%i_start,udims%i_end
              r_u(i,j,k) = r_u(i,j,k) + dubydt_conv_u(i,j,k)*timestep
            END DO
          END DO
        END DO
!$OMP END DO NOWAIT

        ! v: add on to increment field

!$OMP DO SCHEDULE(STATIC)
        DO k =               1,vdims%k_end
          DO j=  vdims%j_start,vdims%j_end
            DO i=vdims%i_start,vdims%i_end
              r_v(i,j,k) = r_v(i,j,k) + dvbydt_conv_v(i,j,k)*timestep
            END DO
          END DO
        END DO
!$OMP END DO NOWAIT

        ! Copy increments of convective momentum if requested
        IF (l_retain_conv_mom_tendencies) THEN
!$OMP DO SCHEDULE(STATIC)
          DO k = 1, udims%k_end
            DO j = udims%j_start,udims%j_end
              DO i = udims%i_start,udims%i_end
                du_conv(i,j,k) = dubydt_conv_u(i,j,k)*timestep
              END DO
            END DO
          END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
          DO k = 1, vdims%k_end
            DO j = vdims%j_start,vdims%j_end
              DO i = vdims%i_start,vdims%i_end
                dv_conv(i,j,k) = dvbydt_conv_v(i,j,k)*timestep
              END DO
            END DO
          END DO
!$OMP END DO NOWAIT
        END IF   !end loop for l_retain_conv_mom_tendencies

!$OMP END PARALLEL

      END IF  ! not SCM running conv in diags-only mode

    END IF ! l_mom

    ! provide some estimate of TKE in convective plumes, based on
    ! the mass flux and convective cloud area
    IF (l_conv_tke .AND. bl_diag%l_tke .AND. l_apply_diag) THEN
      DO k = 1, bl_levels-1
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            ! TKE diag is offset up by a level
            bl_diag%tke(i,j,k+1) = MIN(max_tke,MAX(bl_diag%tke(i,j,k+1),&
                 ( up_flux(i,j,k) / ( g*rho_wet_tq(i,j,k)*              &
                   MIN(0.5,MAX(0.05,cca_2d(i,j))) ) )**2))
                   ! 0.5 and 0.05 are used here as plausible max and min
                   ! values of CCA to prevent numerical problems
          END DO
        END DO
      END DO
    END IF

    ! Output section 5 (convection) stash diagnostics. Moved to atmos_physics2
    ! so that swap_bounds calls for wind increments can be combined with
    ! BL swap_bounds.

    IF ( l_apply_diag ) THEN
      IF (model_type == mt_single_column) THEN

        ! SCM diagnostics
        CALL scm_diag_conv( n_cca_levels, rows, row_length,                 &
                            ntml, ntpar,                                    &
                            ccb, cct, lcbase, lctop, conv_type,             &
                            nscmdpkgs, l_scmdiags,                          &
                            cumulus, l_shallow, l_congestus, l_mid_level,   &
                            wstar, wthvs, cclwp, cca_2d, lcca, deep_flag,   &
                            past_precip, past_conv_ht,                      &
                            dubydt_conv_u, dvbydt_conv_v,                   &
                            theta_star, exner_theta_levels,                 &
                            cca, ccw  )

      ELSE  ! model_type

        IF ( sf(0,5) ) THEN
          CALL diagnostics_conv(                                        &
                         row_length, rows, at_extremity                 &
,                          u, v, p, r_u, r_v                            &
,                          dubydt_conv_p, dvbydt_conv_p                 &
,                          dubydt_conv_u, dvbydt_conv_v                 &
,                          theta_star, q_star                           &
,                          exner_theta_levels                           &
,                          ls_rain, ls_snow                             &
,                          ccw, conv_rain, conv_snow                    &
,                          cca_2d, cca, ccb, cct                        &
,                          cu_over_orog, cape_undilute, cin_undilute    &
,                          lcbase, lctop, lcca, n_cca_levels            &
,                          l_dust                                       &
,                          conv_type                                    &
,                          timestep                                     &
,                                                                       &
                         STASHwork5 )
        END IF  ! test on sf(0,5)

      END IF  ! model_type
    END IF ! l_apply_diag

    ! Deallocate convective momentum transport tendencies
    IF ( l_mom ) THEN
      DEALLOCATE(dvbydt_conv_v)
      DEALLOCATE(dubydt_conv_u)
      DEALLOCATE(dvbydt_conv_p)
      DEALLOCATE(dubydt_conv_p)
    END IF

    !     Deallocate convective diagnostics arrays allocated earlier
    CALL cv_dealloc_diag_array( )

    ! If running the convection scheme in diagnostics-only mode,
    ! zero any convection outputs now that we've output the diagnostics
    IF ( model_type==mt_single_column .AND. conv_mode==1 ) THEN

      DO j = 1, rows
        DO i = 1, row_length
          conv_rain(i,j) = 0.0
          conv_snow(i,j) = 0.0
          ddmfx(i,j) = 0.0
          cclwp0(i,j) = 0.0
          cca0_2d(i,j) = 0.0
        END DO
      END DO

      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
            cca0(i,j,k) = 0.0
            ccw0(i,j,k) = 0.0
          END DO
        END DO
      END DO

    END IF

    IF (Ltimer) CALL timer ('AP2 Convection',6)

  END IF ! l_param_conv  test

  !------------------------------------------------------------------------
  !   End of second communication section
  !------------------------------------------------------------------------

  IF ( i_bl_vn /= i_bl_vn_0 ) THEN
    !   only call on 1st cycle or if not fast running
    IF ( cycleno == 1 .OR. .NOT. l_quick_ap2 ) THEN

      l_calc_at_p = .FALSE.  ! Flag for calling this routine on p-grid
      CALL bdy_expl3 (                                                 &
      ! IN grid related variables
           bl_levels, l_calc_at_p,                                     &
           udims, vdims, udims_s, vdims_s, pdims,                      &
           udims, vdims,                                               &
      ! IN SCM diags
           nSCMDpkgs,L_SCMDiags,                                       &
      ! IN variables used in flux calculations
           u, v, u_0, v_0, rhokm_u_land, rhokm_v_land, flandfac_u,     &
           flandfac_v, rhokm_u_ssi, rhokm_v_ssi, fseafac_u, fseafac_v, &
           flandg_u, flandg_v, zhnl, rdz_u, rdz_v, rhokm_u, rhokm_v,   &
           taux_fd_u, tauy_fd_v,                                       &
           rhogamu_u, rhogamv_v, f_ngstress_u, f_ngstress_v,           &
      ! OUT explicit momentum fluxes
           taux_land, tauy_land, taux_ssi, tauy_ssi, taux, tauy        &
                   )

      IF (l_quick_ap2) THEN
        ! save outputs for second EG cycle
        ! N.B. any new OUT data added to bdy_expl3 will need to be saved here
        ! and restored below. The array size will also need to be increased in
        ! atmos_physics2_alloc
!$OMP PARALLEL DEFAULT(NONE)                                           &
!$OMP SHARED( bl_levels, udims, vdims, bdy_expl3_u3d, taux,            &
!$OMP         bdy_expl3_v3d, tauy, bdy_expl3_u2d, taux_land,           &
!$OMP         taux_ssi, bdy_expl3_v2d, tauy_land, tauy_ssi )           &
!$OMP PRIVATE( i, j, k )
!$OMP DO SCHEDULE(STATIC)
        DO k = 0, bl_levels-1
          DO j = udims%j_start, udims%j_end
            DO i = udims%i_start, udims%i_end
              bdy_expl3_u3d(i,j,k)=taux(i,j,k)    !INOUT to imp_ctl, needs reset
            END DO
          END DO

          DO j = vdims%j_start, vdims%j_end
            DO i = vdims%i_start, vdims%i_end
              bdy_expl3_v3d(i,j,k)=tauy(i,j,k)    !INOUT to imp_ctl, needs reset
            END DO
          END DO
        END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
        DO j = udims%j_start, udims%j_end
          DO i = udims%i_start, udims%i_end
            bdy_expl3_u2d(i,j,1)=taux_land(i,j)   !INOUT to imp_ctl, needs reset
            bdy_expl3_u2d(i,j,2)=taux_ssi(i,j)    !INOUT to imp_ctl, needs reset
          END DO
        END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
        DO j = vdims%j_start, vdims%j_end
          DO i = vdims%i_start, vdims%i_end
            bdy_expl3_v2d(i,j,1)=tauy_land(i,j)   !INOUT to imp_ctl, needs reset
            bdy_expl3_v2d(i,j,2)=tauy_ssi(i,j)    !INOUT to imp_ctl, needs reset
          END DO
        END DO
!$OMP END DO
!$OMP END PARALLEL

      END IF

      DEALLOCATE(rhogamu_u)
      DEALLOCATE(rhogamv_v)
      DEALLOCATE(f_ngstress_u)
      DEALLOCATE(f_ngstress_v)

    ELSE    ! ( cycleno == 1 .OR. .NOT. l_quick_ap2 )
      ! restore outputs on second EG cycle
!$OMP PARALLEL DEFAULT(NONE)                                           &
!$OMP SHARED( bl_levels, udims, vdims, bdy_expl3_u3d, taux,            &
!$OMP         bdy_expl3_v3d, tauy, bdy_expl3_u2d, taux_land,           &
!$OMP         taux_ssi, bdy_expl3_v2d, tauy_land, tauy_ssi )           &
!$OMP PRIVATE( i, j, k )
!$OMP DO SCHEDULE(STATIC)
      DO k = 0, bl_levels-1
        DO j = udims%j_start, udims%j_end
          DO i = udims%i_start, udims%i_end
            taux(i,j,k)=bdy_expl3_u3d(i,j,k)
          END DO
        END DO

        DO j = vdims%j_start, vdims%j_end
          DO i = vdims%i_start, vdims%i_end
            tauy(i,j,k)=bdy_expl3_v3d(i,j,k)
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
      DO j = udims%j_start, udims%j_end
        DO i = udims%i_start, udims%i_end
          taux_land(i,j)=bdy_expl3_u2d(i,j,1)
          taux_ssi(i,j)=bdy_expl3_u2d(i,j,2)
        END DO
      END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
      DO j = vdims%j_start, vdims%j_end
        DO i = vdims%i_start, vdims%i_end
          tauy_land(i,j)=bdy_expl3_v2d(i,j,1)
          tauy_ssi(i,j)=bdy_expl3_v2d(i,j,2)
        END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL

    END IF ! ( cycleno == 1 .OR. .NOT. l_quick_ap2 )
    !-----------------------------------------------------------------------
    ! Reset index for sea-ice categories
    !-----------------------------------------------------------------------
    ! These indices are set up initially in atmos_physics1 using nice_use as
    ! there is the option of only partially using the categories in the
    ! radation and explicit surface scheme.
    ! Here they are reset using nice as all the categories are always used in
    ! the implicit surface scheme.

    ! only call on 1st cycle or if not fast running
    IF ( cycleno == 1 .OR. .NOT. l_quick_ap2 ) THEN
      sice_pts_ncat(:)=0
      sice_index_ncat(:,:)=0
      sice_frac_ncat(:,:)=0.0
      DO n=1,nice
        DO l=1,ssi_pts
          j=(ssi_index(l)-1)/row_length + 1
          i = ssi_index(l) - (j-1)*row_length
          IF (ssi_index(l) > 0) THEN
            IF (ice_fract_ncat(i,j,n) > 0.0) THEN
              sice_pts_ncat(n)=sice_pts_ncat(n)+1
              sice_index_ncat(sice_pts_ncat(n),n)=l
              sice_frac_ncat(l,n)=ice_fract_ncat(i,j,n)
            END IF
          END IF
        END DO
      END DO
    END IF !l_quick_ap2

    IF (model_type == mt_single_column) THEN
      ALLOCATE(theta_star_surf(1,1))
      ALLOCATE(qv_star_surf(1,1))
    ELSE
      IF ( i_pert_theta > pert_theta_mag .OR.                           &
           (sf(487,3) .AND. l_apply_diag) ) THEN
        ALLOCATE(theta_star_surf(tdims%i_start:tdims%i_end,             &
                                 tdims%j_start:tdims%j_end))
      ELSE
        ALLOCATE(theta_star_surf(1,1))
      END IF
      IF ( i_pert_theta == pert_theta_and_moist .OR.                    &
           (sf(488,3) .AND. l_apply_diag) ) THEN
        ALLOCATE(qv_star_surf(tdims%i_start:tdims%i_end,                &
                              tdims%j_start:tdims%j_end))
      ELSE
        ALLOCATE(qv_star_surf(1,1))
      END IF
    END IF

    ! Start of OpenMP parallel region
!$OMP  PARALLEL DEFAULT(NONE)                                           &
!$OMP  SHARED(T_latest, theta_star,                                     &
!$OMP  exner_theta_levels, q_latest, q_star, qcl_latest,                &
!$OMP  qcl_star, qcf_latest, qcf_star, cf_latest, cf_star, cfl_latest,  &
!$OMP  cfl_star, cff_latest, cff_star, tdims)                           &
!$OMP  PRIVATE(i,j,k)

!$OMP DO SCHEDULE(STATIC)
    DO k =                 1, tdims%k_end
      DO j =   tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          T_latest(i,j,k) = theta_star(i,j,k) * exner_theta_levels(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
    DO k =                 1, tdims%k_end
      DO j =   tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          q_latest(i,j,k)   = q_star(i,j,k)
          qcl_latest(i,j,k) = qcl_star(i,j,k)
          qcf_latest(i,j,k) = qcf_star(i,j,k)
          cf_latest(i,j,k)  = cf_star(i,j,k)
          cfl_latest(i,j,k) = cfl_star(i,j,k)
          cff_latest(i,j,k) = cff_star(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

    ! End of OpenMP parallel region
!$OMP END PARALLEL

    CALL NI_imp_ctl (                                                   &
    ! IN model dimensions.
        rhc_row_length, rhc_rows, land_points                           &
      , ntiles, bl_levels                                               &
      , cloud_levels, n_cca_levels                                      &
    ! IN Model switches
      , CycleNo                                                         &
    ! IN model Parameters
      , tr_vars                                                         &
    ! IN trig arrays
      , xx_cos_theta_latitude                                           &
    ! IN data fields.
      , p_layer_centres, p_layer_boundaries, rho_wet_rsq, rho_wet_tq    &
      , u, v, w                                                         &
      , land_sea_mask, q, qcl, qcf, p_star, theta, exner_theta_levels   &
    ! IN ancillary fields and fields needed to be kept from tstep to tstep
      , sil_orog_land, ho2r2_orog                                       &
      , ice_fract, di_ncat, ice_fract_ncat, k_sice, u_0, v_0, land_index&
      , cca, lcbase, ccb0, cct0                                         &
      , ls_rain, ls_snow, conv_rain, conv_snow, L_scrn, L_plsp          &
    ! IN variables required from BDY_LAYR
      , alpha1_sea, alpha1_sice, ashtf_sea, ashtf, bq_gb, bt_gb         &
      , dtrdz_charney_grid, rdz_charney_grid, dtrdz_u, dtrdz_v          &
      , rdz_u, rdz_v, cdr10m_u, cdr10m_v, cdr10m_n_u, cdr10m_n_v        &
      , cd10m_n_u, cd10m_n_v, z_theta                                   &
      , k_blend_tq, k_blend_uv, uStarGBM, rhokm, rhokm_u, rhokm_v       &
    ! IN diagnostics (started or from) BDY_LAYR
      , rib_gb,zlcl,zht,zhnl,dzh,qcl_inv_top,zh                         &
      , bl_type_1,bl_type_2,bl_type_3,bl_type_4,bl_type_5,bl_type_6     &
      , bl_type_7, z0m_gb, z0m_eff_gb, z0h_eff_gb                       &
      , ntml, cumulus, l_pc2_diag_sh_pts                                &
    ! IN data required for tracer mixing :
      , rho_aresist,aresist,r_b_dust                                    &
      , kent, we_lim, t_frac, zrzi                                      &
      , kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc                      &
      , zhsc,z_rho,dust_flux,dust_emiss_frac                            &
      , u_s_t_tile,u_s_t_dry_tile,u_s_std_tile                          &
    ! IN additional variables for JULES. Now includes lai_ft, canht_ft.
      , tile_pts,tile_index,tile_frac,canopy                            &
      , alpha1,fraca,rhokh_tile,smc,chr1p5m,resfs,z0hssi,z0mssi         &
      , canhc_tile,flake,wt_ext_tile,lw_down,lai_ft,canht_ft            &
      , sw_tile,ashtf_tile,gc,aresist_tile                              &
      , resft,rhokh_sice,rhokh_sea,z0h_tile,z0m_tile                    &
      , chr1p5m_sice                                                    &
      , fland, flandg, flandg_u,flandg_v,vshr_land,vshr_ssi             &
      , emis_tile, t_soil, snow_tile, rib_ssi                           &
    ! IN JULES variables for STASH
      , gs,gpp,npp,resp_p,gpp_ft,npp_ft,resp_p_ft,resp_s                &
      , resp_s_tot,cs                                                   &
      , rib_tile,fsmc,catch,g_leaf                                      &
      , co2_emits, co2flux                                              &
    ! IN additional variables for soil moisture nudging scheme
      , wt_ext,                                                         &
    ! INOUT diagnostic info
        STASHwork3, STASHwork9                                          &
    ! SCM Diagnostics (dummy in full UM) & bl diags
      , nSCMDpkgs, L_SCMDiags, bl_diag, sf_diag                         &
    ! INOUT (Note ti and ti_gb are IN only if l_sice_multilayers=T)
      , TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans                          &
      , cca0,ccw0,cca0_2d,fqw,ftl,taux,tauy, rhokh                      &
      , fqw_ice,ftl_ice,dtstar_tile,dtstar_sea,dtstar_sice,ti           &
      , area_cloud_fraction, bulk_cloud_fraction                        &
      , T_latest, q_latest, qcl_latest, qcf_latest                      &
      , cf_latest, cfl_latest, cff_latest                               &
      , R_u, R_v, R_w, cloud_fraction_liquid, cloud_fraction_frozen     &
      , sum_eng_fluxes,sum_moist_flux, rhcpt                            &
    ! INOUT tracer fields
      , aerosol, free_tracers,  resist_b,  resist_b_tile                &
      , dust_div1,dust_div2,dust_div3,dust_div4,dust_div5,dust_div6     &
      , drydep2, so2, dms, so4_aitken, so4_accu, so4_diss, nh3          &
      , soot_new, soot_aged, soot_cld, bmass_new, bmass_aged            &
      , bmass_cld, ocff_new, ocff_aged, ocff_cld, nitr_acc, nitr_diss   &
      , co2, ozone_tracer                                               &
    ! INOUT additional variables for JULES
      , tstar_tile,fqw_tile,epot_tile,ftl_tile                          &
      , radnet_sea,radnet_sice,olr,tstar_sice_cat,tstar_ssi             &
      , tstar_sea,taux_land,taux_ssi,tauy_land,tauy_ssi,Error_code      &
    ! OUT fields
      , surf_ht_flux_land, zlcl_mixed                                   &
      , theta_star_surf, qv_star_surf                                   &
    ! OUT additional variables for JULES
      , tstar, ti_gb, ext, snowmelt,tstar_land,tstar_sice, ei_tile      &
      , ecan_tile,melt_tile, surf_htf_tile                              &
      )
    !

    !     Deallocate diagnostic workspace.
    IF ( (CycleNo == NumCycles) .OR. (.NOT. l_quick_ap2) ) THEN
      DEALLOCATE(cdr10m_n)
      DEALLOCATE(cdr10m_n_u)
      DEALLOCATE(cdr10m_n_v)
      DEALLOCATE(cd10m_n)
      DEALLOCATE(cd10m_n_u)
      DEALLOCATE(cd10m_n_v)
    END IF
    !

    !-----------------------------------------------------------
    ! check moisture conservation for Boundary layer
    ! only call on final cycle
    !-----------------------------------------------------------
    IF (l_check_moist_inc .AND. cycleno == NumCycles                    &
                  .AND. model_type == mt_global) THEN
      ALLOCATE (ep1(row_length, rows) )
      ALLOCATE (ep2(row_length, rows) )
      ALLOCATE (ep3(row_length, rows) )
      ALLOCATE (dqt(row_length, rows,model_levels) )
      ! Sum up increments from Boundary layer scheme and convert to rate
      ! At present there are no increments to qrain, qcf2 or qgraup possible
      ! from the boundary layer scheme

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( model_levels, rows, row_length, dqt, bl_diag,            &
!$OMP         recip_timestep )                                         &
!$OMP PRIVATE( i, j, k )
      DO k=1,model_levels
        DO j=1,rows
          DO i=1,row_length
            dqt(i,j,k) = (bl_diag%q_incr(i,j,k)+bl_diag%qcl_incr(i,j,k)+&
                           bl_diag%qcf_incr(i,j,k)) *recip_timestep
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO
      scheme_name = 'Boundary layer'

      IF (l_mr_physics) THEN
        ! increments are mixing ratios need dry density
        ALLOCATE (rho_dry_r2(row_length, rows,model_levels) )

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( model_levels, rows, row_length, rho_dry_r2,              &
!$OMP         unscaled_dry_rho, r_rho_levels )                         &
!$OMP PRIVATE( i, j, k )
        DO k=1,model_levels
          DO j=1,rows
            DO i=1,row_length
              rho_dry_r2(i,j,k) = unscaled_dry_rho(i,j,k)*              &
                             r_rho_levels(i,j,k)*r_rho_levels(i,j,k)
            END DO
          END DO
        END DO
!$OMP END PARALLEL DO

        CALL check_dmoist_inc(row_length, rows, model_levels,           &
                      delta_lambda,delta_phi,timestep,                  &
                      p_layer_boundaries,rho_dry_r2,                    &
                      dqt,                                              &
                      fqw(1,1,1),  scheme_name,                         &
                      ep1, ep2, ep3)

        DEALLOCATE(rho_dry_r2)

      ELSE
        ! Use wet density
        CALL check_dmoist_inc(row_length, rows, model_levels,           &
                      delta_lambda,delta_phi,timestep,                  &
                      p_layer_boundaries,rho_wet_rsq,                   &
                      dqt,                                              &
                      fqw(1,1,1), scheme_name,                          &
                      ep1, ep2, ep3)

      END IF
      DEALLOCATE(dqt)
      DEALLOCATE(ep3)
      DEALLOCATE(ep2)
      DEALLOCATE(ep1)

    END IF
    ! Deallocation of some BL_diag moved from ni_imp_ctl as required at
    ! this level if checking moisture conservation

    DEALLOCATE (BL_diag%qcf_incr)
    DEALLOCATE (BL_diag%qcl_incr)
    DEALLOCATE (BL_diag%q_incr)

    ! ----------------------------------------------------------------------
    ! Section RESTORE: Copy _latest fields into _star locations
    ! Area cloud fraction has been done in NI_IMP_CTL
    ! ----------------------------------------------------------------------

      ! Start of OpenMP parallel region
!$OMP  PARALLEL DEFAULT(NONE)                                           &
!$OMP  SHARED(theta_star, T_latest,                                     &
!$OMP  exner_theta_levels, q_star, q_latest, qcl_star,                  &
!$OMP  qcl_latest, qcf_star, qcf_latest, cf_star, cf_latest,            &
!$OMP  cfl_star, cfl_latest, cff_star, cff_latest,                      &
!$OMP  tdims)                                                           &
!$OMP  PRIVATE(i,j,k)

!$OMP DO SCHEDULE(STATIC)
    DO k =                 1, tdims%k_end
      DO j =   tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          theta_star(i,j,k) = T_latest(i,j,k)/exner_theta_levels(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
    DO k =                 1, tdims%k_end
      DO j =   tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          q_star(i,j,k)   = q_latest(i,j,k)
          qcl_star(i,j,k) = qcl_latest(i,j,k)
          qcf_star(i,j,k) = qcf_latest(i,j,k)
          cf_star(i,j,k)  = cf_latest(i,j,k)
          cfl_star(i,j,k) = cfl_latest(i,j,k)
          cff_star(i,j,k) = cff_latest(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

    ! End of OpenMP parallel region
!$OMP END PARALLEL

  ELSE ! i_bl_vn

    DO k = 1, model_levels
      rhcpt(1,1,k) = rhcrit(k)
    END DO

    CALL eg_idl_friction(u, v, exner_theta_levels, exner_rho_levels,    &
                         theta_star, q_star, p_theta_levels,            &
                         r_u, r_v, error_code, l_mr_physics)

  END IF ! i_bl_vn

END IF ! error_code == 0

IF (Ltimer) CALL timer ('AP2 Implicit BL',6)
IF (Ltimer) CALL timer ('AP2 Boundary Layer',6)

IF ( error_code  ==  0 ) THEN

  !   Reset halo variables for writing back into D1 array on exit from
  !   this subroutine.
  IF (i_cld_vn /= i_cld_pc2 ) THEN

!$OMP  PARALLEL DEFAULT(NONE)                                           &
!$OMP  SHARED(                                                          &
!$OMP  bulk_cloud_fraction_halos, bulk_cloud_fraction,                  &
!$OMP  cloud_fraction_liquid_halos, cloud_fraction_liquid,              &
!$OMP  cloud_fraction_frozen, cloud_fraction_frozen_halos,              &
!$OMP  tdims)                                                           &
!$OMP  PRIVATE(i,j,k)

!$OMP DO SCHEDULE(STATIC)
    DO k =                 1, tdims%k_end
      DO j =   tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          bulk_cloud_fraction_halos(i,j,k)   = bulk_cloud_fraction(i,j,k)
          cloud_fraction_liquid_halos(i,j,k) = cloud_fraction_liquid(i,j,k)
          cloud_fraction_frozen_halos(i,j,k) = cloud_fraction_frozen(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

    ! End of OpenMP parallel region
!$OMP END PARALLEL

  END IF  ! .not. i_cld_pc2

END IF ! on error code equal to zero

IF (i_bl_vn /= i_bl_vn_0) THEN

   ! random perturbations of theta_star
  IF (i_pert_theta > off) THEN
    CALL bl_pert_theta(cycleno, theta_star, theta_star_surf,            &
                       q_star, qv_star_surf, ntml, cumulus,             &
                       z_theta, STASHwork3)
  END IF

  DEALLOCATE(theta_star_surf)
  DEALLOCATE(qv_star_surf)

END IF

! Only do land surface on last cycle
IF ( CycleNo == NumCycles ) THEN

  IF ( .NOT. l_triffid) THEN
    DO j = 1,soil_points
      i = soil_index(j)
      cs_ch4(i)=cs(i,1)
    END DO
  END IF

  SELECT CASE (model_type)

  CASE (mt_single_column)
    sf_diag%l_snice = .FALSE.

  CASE DEFAULT
    sf_diag%l_snice = ( sf(578,8) .OR. sf(579,8) .OR.  &
                 sf(580,8) .OR. sf(581,8) .OR.         &
                 sf(582,8) .OR. sf(583,8) ) .AND.      &
                 (cycleno == numcycles .OR. l_quick_ap2)

  END SELECT
  !
  IF (sf_diag%l_snice) THEN
    ALLOCATE(sf_diag%snice_smb_surft(land_points,ntiles))
    ALLOCATE(sf_diag%snice_m_surft(land_points,ntiles))
    ALLOCATE(sf_diag%snice_freez_surft(land_points,ntiles))
    ALLOCATE(sf_diag%snice_runoff_surft(land_points,ntiles))
    ALLOCATE(sf_diag%snice_sicerate_surft(land_points,ntiles))
    ALLOCATE(sf_diag%snice_sliqrate_surft(land_points,ntiles))
  ELSE
    !       Minimal sizes are alllocated to ensure valid addresses are
    !       present in subroutine calls.
    ALLOCATE(sf_diag%snice_smb_surft(1,1))
    ALLOCATE(sf_diag%snice_m_surft(1,1))
    ALLOCATE(sf_diag%snice_freez_surft(1,1))
    ALLOCATE(sf_diag%snice_runoff_surft(1,1))
    ALLOCATE(sf_diag%snice_sicerate_surft(1,1))
    ALLOCATE(sf_diag%snice_sliqrate_surft(1,1))
  END IF

  CALL surf_couple_extra(                                               &
  !IN
   land_points, row_length, rows, river_row_length, river_rows, land_index,&
   ls_rain, conv_rain, ls_snow, ls_graup, conv_snow, surf_ht_flux_land, &
   cca_2d, smlt, ntiles, tile_pts, tile_index, tile_frac,               &
   ei_tile, surf_htf_tile,                                              &
   tstar_tile, hcons,                                                   &
   land_ice_points, land_ice_index, soil_points, soil_index, ext,       &
   stf_sub_surf_roff,                                                   &
   ecan_tile, fexp, gamtot, ti_mean, ti_sig, cs_ch4, p_star,            &
   a_fsat, c_fsat, a_fwet, c_fwet,                                      &
   ntype, fqw_tile,                                                     &
   halo_i, halo_j, model_levels,                                        &
   delta_lambda, delta_phi, xx_cos_theta_latitude, i_river_vn,          &
   aocpl_row_length, aocpl_p_rows, xpa, xua, xva, ypa, yua, yva,        &
   g_p_field, g_r_field, n_proc, global_row_length, global_rows,        &
   global_river_row_length, global_river_rows, flandg, river_vel, river_mcoef,&
   trivdir, trivseq, r_area, slope, flowobs1, r_inext, r_jnext, r_land, &
   substore, surfstore, flowin, bflowin, smvcst, smvcwt,                &
   a_step, n_rows, offx, offy, n_procx, n_procy, g_rows, g_row_length,  &
   at_extremity, frac_disturb, satcon, soil_clay, resp_s, npp, z0m_soil,&
   !INOUT
   a_steps_since_riv, melt_tile, t_soil, tsurf_elev_surft,              &
   rgrain, snow_grnd, snow_tile,                                        &
   soil_layer_moisture, sthf, sthu, canopy, fsat, fwetl, zw, sthzw,     &
   snow_depth, snowmelt, ls_rainfrac,                                   &
   tot_surf_runoff, tot_sub_runoff, acc_lake_evap, twatstor,            &
   asteps_since_triffid, g_leaf_acc, g_leaf_phen_acc, npp_ft_acc, resp_s_acc,&
   resp_w_ft_acc, cs, frac, lai_ft, canht_ft, catch_snow, catch, infil_tile,&
   !OUT- mostly for SCM output calls below
   inlandout_atm, surf_ht_flux_ld, snow_melt,                           &
   snomlt_sub_htf, dhf_surf_minus_soil,                                 &
   canopy_gb, snomlt_surf_htf, smc, sub_surf_roff, surf_roff,           &
   tot_tfall, z0_tile, z0h_tile_bare, snow_soil_htf,                    &
   land_sea_mask                                                        &
   )

  IF ( l_nitrogen ) THEN
    ! JULES has updated the N-pools, now give those updated values to the
    ! UM equivalents.
    ! UM prognostic = JULES prognostic
    soil_nitro1(:)  = ns_pool_gb(:,1,1)
    soil_nitro2(:)  = ns_pool_gb(:,1,2)
    soil_nitro3(:)  = ns_pool_gb(:,1,3)
    soil_nitro4(:)  = ns_pool_gb(:,1,4)
    soil_inorgnit(:)  = n_inorg_gb(:)
    nitrogen_deposition_d1(:) = deposition_n_gb(:)
  END IF   !l_nitrogen

  ! JULES has updated the wood product pools and previous agricultural
  ! fraction. Now give those updated values to the UM equivalents. IF TRIFFID
  ! has not been called this timestep it will effectively do nothing, passing
  ! non-updated values from the JULES module to the UM.

  IF (l_triffid) THEN

    IF (l_landuse) THEN
      wood_prod_fast_d1 = wood_prod_fast_gb
      wood_prod_med_d1  = wood_prod_med_gb
      wood_prod_slow_d1 = wood_prod_slow_gb
    END IF

    IF ( l_nitrogen ) THEN
      soil_nitro1(:)            = ns_pool_gb(:,1,1)
      soil_nitro2(:)            = ns_pool_gb(:,1,2)
      soil_nitro3(:)            = ns_pool_gb(:,1,3)
      soil_nitro4(:)            = ns_pool_gb(:,1,4)
      soil_inorgnit(:)          = n_inorg_gb(:)
      nitrogen_deposition_d1(:) = deposition_n_gb(:)
    END IF

    IF (l_trif_crop) THEN
      ! ...we have both pasture and crop fractions; the crop fraction
      ! makes use of the disturbed fraction variables.
      agr_crop_frac_prev_d1 = frac_agr_prev_gb
      pasture_frac_d1       = frac_past_gb
      pasture_frac_prev_d1  = frac_past_prev_gb
    ELSE
      ! ...we have just the disturbed fraction variable
      IF (l_landuse) THEN
        disturb_veg_prev      = frac_agr_prev_gb
      END IF
    END IF ! l_trif_crop

    ! If TRIFFID has just been called it will have updated triffid_co2_gb,
    ! the CO2 flux from exudates + decay of the wood product pools and + 
    ! crop harvest, combined. Copy that to the equivalent UM variable for 
    ! writing to the start dumps so that it can be added to the atmosphere 
    ! by BL_TRMIX_DD next time that is called (which might be the start of 
    ! the next CRUN). 

    IF ( l_co2_interactive ) THEN
      triffid_co2_d1(:) =  triffid_co2_gb(:)
    END IF ! l_co2_interactive

  END IF ! l_triffid



  ! MORUSES prognostics are not currently updated, but some may be in the
  ! future to account for snow.
  IF ( l_urban2t ) THEN
    hgt   = hgt_gb
    hwr   = hwr_gb
    wrr   = wrr_gb
    disp  = disp_gb
    ztm   = ztm_gb
    albwl = albwl_gb
    albrd = albrd_gb
    emisw = emisw_gb
    emisr = emisr_gb
  END IF

  !   Map module prognostics back to JULES prognostics
  snowdepth_p = snowdepth_surft
  IF ( nsmax > 0 ) THEN
    rho_snow_grnd_p = rho_snow_grnd_surft
    nsnow_p         = nsnow_surft
    ds_p            = ds_surft
    sice_p          = sice_surft
    sliq_p          = sliq_surft
    tsnowlayer_p    = tsnow_surft
    rho_snow_p      = rho_snow_surft
    rgrainl_p       = rgrainl_surft
  END IF ! NSMAX

  DEALLOCATE(sf_diag%snice_sliqrate_surft)
  DEALLOCATE(sf_diag%snice_sicerate_surft)
  DEALLOCATE(sf_diag%snice_runoff_surft)
  DEALLOCATE(sf_diag%snice_freez_surft)
  DEALLOCATE(sf_diag%snice_m_surft)
  DEALLOCATE(sf_diag%snice_smb_surft)

  !-------------------------------------------------------------------------
  !End of JULES coupling point- FLake and land diagnostic calls
  !-------------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !     call to the FLake interface
  !-----------------------------------------------------------------------
  IF (       l_flake_model                &
       .AND. (.NOT. l_aggregate)           &
       .AND. (land_points > 0 ) ) THEN

    DO k=1,tile_pts(lake)
      l = tile_index(k,lake)
      j=(land_index(l)-1)/row_length + 1
      i = land_index(l) - (j-1)*row_length

      !       U* of lake obtained by matching surface stress
      u_s_lake_gb(l) =   u_s_std_tile(l, lake)                          &
                    * SQRT( rho_wet_tq(i,j,1) / rho_water )

      !       Downwelling SW on lake tile
      sw_down_gb(l) = sw_tile(l,lake) / (1.0 - lake_albedo_gb(l))

      !       Take the net SW flux out of the surface heat flux
      !       since this is done separately within FLake.
      surf_ht_flux_lk_gb(l) =   surf_ht_flux_lake_ij(i,j) - sw_tile(l,lake)

      IF ( (nsmax > 0) .AND. (lake_h_snow_gb(l) > EPSILON(1.0)) ) THEN

        ! For the new snow scheme, FLake is forced with zero snow
        ! and the forcing fluxes are those at the bottom of the snowpack.
        ! The order of the following calculations is important.

        ! attenuated DWSW below snow
        dwsw_sub_snow = sw_down_gb(l) * EXP( -lake_h_snow_gb(l)/h_snow_sw_att )

        ! following last calculation, set snow depth to zero
        lake_h_snow_gb(l) = 0.0

        ! heat flux into the lake becomes
        ! the sub-snow value minus remaining DWSW

        surf_ht_flux_lk_gb(l) =   surf_ht_flux_lake_ij(i,j)    &
                                - dhf_surf_minus_soil(l)       &
                                - dwsw_sub_snow

        ! now overwrite the DWSW with the remaining sub-snow amount
        sw_down_gb(l) = dwsw_sub_snow

      END IF

    END DO

    trap_frozen   = 0
    trap_unfrozen = 0

    ! DEPENDS ON: flake_interface
    CALL flake_interface( land_points                                   &
                         ,tile_pts(lake)                                &
                         ,tile_index(:,lake)                            &
                         ,u_s_lake_gb                                   &
                         ,surf_ht_flux_lk_gb                            &
                         ,sw_down_gb                                    &
                         ,lake_depth_gb                                 &
                         ,lake_fetch_gb                                 &
                         ,coriolis_param_gb                             &
                         ,timestep                                      &
                         ,lake_albedo_gb                                &
                         ,lake_t_snow_gb                                &
                         ,lake_t_ice_gb                                 &
                         ,lake_t_mean_gb                                &
                         ,lake_t_mxl_gb                                 &
                         ,lake_shape_factor_gb                          &
                         ,lake_h_snow_gb                                &
                         ,lake_h_ice_gb                                 &
                         ,lake_h_mxl_gb                                 &
                         ,lake_t_sfc_gb                                 &
                         ,ts1_lake_gb                                   &
                         ,g_dt_gb                                       &
                         ,trap_frozen                                   &
                         ,trap_unfrozen )

    IF ( printstatus > PrStatus_Oper ) THEN
      IF ( trap_frozen > 0 ) THEN
        WRITE(umMessage,*)                                              &
          'AP2-FLake: # zero-divide (  frozen) avoided =',              &
           trap_frozen
        CALL umPrint(umMessage,src='atmos_physics2')
      END IF
      IF ( trap_unfrozen > 0 ) THEN
        WRITE(umMessage,*)                                              &
          'AP2-FLake: # zero-divide (unfrozen) avoided =',              &
          trap_unfrozen
        CALL umPrint(umMessage,src='atmos_physics2')
      END IF
    END IF

    !       copy FROM FLake module variables TO Prognostics
    lake_depth_p  = lake_depth_gb
    lake_fetch_p  = lake_fetch_gb
    lake_t_mean_p = lake_t_mean_gb
    lake_t_mxl_p  = lake_t_mxl_gb
    lake_t_ice_p  = lake_t_ice_gb
    lake_h_mxl_p  = lake_h_mxl_gb
    lake_h_ice_p  = lake_h_ice_gb
    lake_shape_p  = lake_shape_factor_gb
    lake_g_dt_p   = g_dt_gb

  END IF ! Flake

  !-------------------------------------------------------------------------
  !   Set Coastal tiling dependent prognostics for output to D1 array:
  !-------------------------------------------------------------------------
  IF ( L_ctile ) THEN

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(k,j,i)                                &
!$OMP SHARED (rows,row_length,tstar_land_ctile,tstar_land,tstar_sice_ctile,&
!$OMP    tstar_sice,nice_use,tstar_sice_cat_ctile,tstar_sice_cat,tstar_sea,&
!$OMP    tstar_sea_ctile)

!$OMP DO SCHEDULE(STATIC)
    DO j = 1, rows
      DO i = 1, row_length
        tstar_land_ctile(i,j)=tstar_land(i,j)
        tstar_sea_ctile(i,j)=tstar_sea(i,j)
        tstar_sice_ctile(i,j)=tstar_sice(i,j)
      END DO
    END DO
!$OMP END DO NOWAIT

    DO k = 1, nice_use
!$OMP DO SCHEDULE(STATIC)
      DO j = 1, rows
        DO i = 1, row_length
          tstar_sice_cat_ctile(i,j,k)=tstar_sice_cat(i,j,k)
        END DO
      END DO
!$OMP END DO NOWAIT
    END DO

!$OMP END PARALLEL

  END IF ! L_ctile

END IF ! CycleNo == NumCycles

IF (ALLOCATED(sf_diag%t1p5m_surft)) THEN
  DEALLOCATE(sf_diag%t1p5m_surft)
END IF

IF (ALLOCATED(sf_diag%q1p5m_surft)) THEN
  DEALLOCATE(sf_diag%q1p5m_surft)
END IF


IF (l_param_conv) THEN
  ! This code is not required and should not be executed if no convection 
  ! scheme is called. The setting of l_ccrad should not be requried so by 
  ! default is .false. 
  ! Copy ccw/lcbase values for Radiation if l_ccrad in use.
  ! See comments at declaration of ccw_out/lcbase_out
  IF (l_ccrad) THEN

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(k,j,i)                  &
!$OMP SHARED(model_levels,rows,row_length,ccw_out,ccw0,      &
!$OMP    lcbase_out,lcbase0)

!$OMP DO SCHEDULE(STATIC)
    DO k=1, model_levels
      DO j=1, rows
        DO i=1, row_length
          ccw_out(i,j,k) = ccw0(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
    DO j=1, rows
      DO i=1, row_length
        lcbase_out(i,j)  = lcbase0(i,j)
      END DO
    END DO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

  END IF ! l_ccrad

END IF ! l_param_conv


!
! ----------------------------------------------------------------------+-------
! Whether model_type is scm or not, replicate total precipitation rate for dump
! to pass as start-of-timestep input to s6 (gw_ussp), but check if dump space
! has been requested. For code as introduced l_use_ussp is sufficient
! but more rigorous to define l_totppn_softs if requirement moves beyond
! the simple (l_totppn_softs == l_use_ussp).
! Note that conv_rain/snow default to 0 (above) if Convection switched off.
! ----------------------------------------------------------------------+-------
tppsofts_if1: IF (l_use_ussp) THEN  

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(j,i)                               &
!$OMP SHARED (tdims,totalppn,ls_rain,conv_rain,ls_snow,conv_snow)

!$OMP DO SCHEDULE(STATIC)
  Rows_tppsoftsdo1: DO j = tdims%j_start,tdims%j_end
    Row_length_tppsoftsdo1: DO i = tdims%i_start,tdims%i_end
      totalppn(i,j) = (ls_rain(i,j) + conv_rain(i,j) +                  &
                       ls_snow(i,j) + conv_snow(i,j)  )
    END DO  Row_length_tppsoftsdo1
  END DO  Rows_tppsoftsdo1
!$OMP END DO NOWAIT

!$OMP END PARALLEL

END IF  tppsofts_if1

IF (model_type == mt_single_column) THEN
  !-------------------------------------------------------------------------
  !   Output some SCM diagnostics
  !-------------------------------------------------------------------------

      ! The SCM output diagnostic system needs to know that
      ! sub-stepping has ended
  ! DEPENDS ON: scm_substepping_end
  CALL SCM_substepping_end

  !-----------------------------------------------------------------------
  !   SCM Boundary Layer Diagnostics Package
  !-----------------------------------------------------------------------
  IF (L_SCMDiags(SCMDiag_bl)) THEN

    CALL SCMoutput(rib_gb,'rib',                                        &
         'Richardson number (lowest layer)','-',                        &
         t_avg,d_sl,default_streams,'',RoutineName)

    CALL SCMoutput(taux,'taux',                                         &
         'Wind stress x','N/m2',                                        &
         t_avg,d_bl,default_streams,'',RoutineName)

    CALL SCMoutput(tauy,'tauy',                                         &
         'Wind stress y','M/m2',                                        &
         t_avg,d_bl,default_streams,'',RoutineName)

    CALL SCMoutput(fqw,'fqT',                                           &
         'Turbulent moisture flux','kg/m2/s',                           &
         t_avg,d_bl,default_streams,'',RoutineName)

    CALL SCMoutput(ftl,'ftl',                                           &
         'Turbulent sensible heat flux','W/m2',                         &
         t_avg,d_bl,default_streams,'',RoutineName)

    IF ( l_calc_tau_at_p ) THEN
      ! Explicit wind stresses on the theta-grid
      CALL SCMoutput(taux_p,'taux_p',                                   &
         'Explicit wind stress x at theta','N/m2',                      &
         t_avg,d_bl,default_streams,'',RoutineName)
      CALL SCMoutput(tauy_p,'tauy_p',                                   &
         'Explicit wind stress y at theta','N/m2',                      &
         t_avg,d_bl,default_streams,'',RoutineName)
    END IF

  END IF ! L_SCMDiags(SCMDiag_bl)
  !
  !---------------------------------------------------------------------
  !   SCM Land Points Diagnostics Package
  !---------------------------------------------------------------------
  IF (L_SCMDiags(SCMDiag_land)) THEN

    CALL SCMoutput(ftl_tile,'shf_tile',                                 &
         'Tile sensible heat flux','W/m2',                              &
         t_avg,d_tile,default_streams,'',RoutineName)

    CALL SCMoutput(surf_ht_flux_land,'surf_ht_flux_ld',                 &
         'Net downward heat flux into land fraction','W/m2',            &
         t_avg,d_sl,default_streams,'',RoutineName)

    CALL SCMoutput(snomlt_sub_htf,'snomlt_sub_htf',                     &
         'Snow melt heat flux into sub-surface','W/m2',                 &
         t_avg,d_land,default_streams,'',RoutineName)

    CALL SCMoutput(snow_melt,'snomlt_hyd',                              &
         'Snow melt (hydrology scheme)','kg/m2/day',                    &
         t_mult,d_land,default_streams,'sec_day',RoutineName)

    CALL SCMoutput(snowmelt,'snomlt_bl',                                &
         'Snow melt (boundary layer scheme)','kg/m2/day',               &
         t_mult,d_land,default_streams,'sec_day',RoutineName)

    CALL SCMoutput(tot_tfall,'thro_fall',                               &
         'Throughfall','kg/m2/dat',                                     &
         t_mult,d_land,default_streams,'sec_day',RoutineName)

    CALL SCMoutput(surf_roff,'surf_roff',                               &
         'Surface runoff','kg/m2/day',                                  &
         t_mult,d_land,default_streams,'sec_day',RoutineName)

    CALL SCMoutput(sub_surf_roff,'Sub_roff',                            &
         'Sub-surface runoff','kg/m2/day',                              &
         t_mult,d_land,default_streams,'sec_day',RoutineName)

    CALL SCMoutput(gpp,'gpp',                                           &
         'Gross primary productivity','kg C/m2/s',                      &
         t_mult,d_land,default_streams,'oneKsecday',RoutineName)

    CALL SCMoutput(npp,'npp',                                           &
         'Net primary productivity','kg C/m2/s',                        &
         t_mult,d_land,default_streams,'oneKsecday',RoutineName)

    CALL SCMoutput(resp_p,'resp_p',                                     &
         'Plant respiration','kg/m2/s',                                 &
         t_mult,d_land,default_streams,'oneKsecday',RoutineName)

    CALL SCMoutput(snomlt_surf_htf,'snomlt_surf_htf',                   &
         'Snowmelt heatflux (boundary layer scheme)','W/m2',            &
         t_avg,d_sl,default_streams,'',RoutineName)

  END IF ! L_SCMDiags(SCMDiag_land)
  !

  !---------------------------------------------------------------------
  !   SCM General OR Surface Diagnostics Package
  !---------------------------------------------------------------------
  IF (L_SCMDiags(SCMDiag_gen)                                           &
      .OR. L_SCMDiags(SCMDiag_surf)) THEN

    CALL SCMoutput(ftl(1,1,1),'sens_ht',                                &
         'Surface sensible heat flux','W/m2',                           &
         t_avg,d_sl,default_streams,'',RoutineName)

  END IF ! L_SCMDiags(SCMDiag_surf)
  !
  !---------------------------------------------------------------------
  !   SCM Convection OR Increments Diagnostics Package
  !---------------------------------------------------------------------
  IF (l_param_conv) THEN
    IF (L_SCMDiags(SCMDiag_conv)                                        &
         .OR. L_SCMDiags(SCMDiag_incs)) THEN

      CALL SCMoutput(theta_incr_diag_conv,'dth_conv',                   &
           'Convective increment - theta','K',                          &
           t_avg,d_all,default_streams,'',RoutineName)

    END IF ! L_SCMDiags(SCMDiag_conv)
  END IF
  !
  !---------------------------------------------------------------------
  !   SCM General OR Large Scale Precipitation Diagnostics Package
  !---------------------------------------------------------------------
  IF (L_SCMDiags(SCMDiag_gen) .OR. L_SCMDiags(SCMDiag_lsp)) THEN

    CALL SCMoutput(ls_rain,'ls_rain',                                   &
         'Large scale rainfall rate','kg/m2/day',                       &
         t_mult,d_sl,default_streams,'sec_day',RoutineName)

    CALL SCMoutput(ls_snow,'ls_snow',                                   &
         'Large scale snowfall rate','kg/m2/day',                       &
         t_mult,d_sl,default_streams,'sec_day',RoutineName)

    CALL SCMoutput(ls_rain,'ls_rain_inst',                              &
         'Large scale rainfall rate','kg/m2/s',                         &
         t_inst,d_sl,default_streams,'',RoutineName)

    CALL SCMoutput(ls_snow,'ls_snow_inst',                              &
         'Large scale snowfall rate','kg/m2/s',                         &
         t_inst,d_sl,default_streams,'',RoutineName)

  END IF ! L_SCMDiags(SCMDiag_gen) .OR. L_SCMDiags(SCMDiag_lsp)
  !
  !-----------------------------------------------------------------------
  !   SCM General OR Convection Diagnostics Package
  !---------------------------------------------------------------------
  IF (L_SCMDiags(SCMDiag_gen) .OR. L_SCMDiags(SCMDiag_conv)) THEN

    CALL SCMoutput(conv_rain,'conv_rain',                               &
         'Convective rainfall','kg/m2/day',                             &
         t_mult,d_sl,default_streams,'sec_day',RoutineName)

    CALL SCMoutput(conv_snow,'conv_snow',                               &
         'Convective snowfall','kg/m2/day',                             &
         t_mult,d_sl,default_streams,'sec_day',RoutineName)

    CALL SCMoutput(conv_rain,'conv_rain_inst',                          &
         'Convective rainfall','kg/m2/s',                               &
         t_inst,d_sl,default_streams,'',RoutineName)

    CALL SCMoutput(conv_snow,'conv_snow_inst',                          &
         'Convective snowfall','kg/m2/s',                               &
         t_inst,d_sl,default_streams,'',RoutineName)

  END IF ! L_SCMDiags(SCMDiag_gen) .OR. L_SCMDiags(SCMDiag_conv)
  !
  !---------------------------------------------------------------------
  !   SCM Large Scale Precipitation OR Convection OR General
  !   Diagnostics Package
  !---------------------------------------------------------------------
  IF (L_SCMDiags(SCMDiag_lsp) .OR. L_SCMDiags(SCMDiag_conv)             &
      .OR. L_SCMDiags(SCMDiag_gen)) THEN

    ! Calculate precipitation rates
    DO j=1,rows
      DO i=1,row_length
        tot_rain(i,j)   = ls_rain(i,j)   + conv_rain(i,j)
        tot_snow(i,j)   = ls_snow(i,j)   + conv_snow(i,j)
        pptn_rate(i,j)  = ls_rain(i,j)   + ls_snow(i,j)                 &
                        + conv_rain(i,j) + conv_snow(i,j)
        accum_pptn(i,j) = timestep * pptn_rate(i,j)
      END DO ! i
    END DO ! j

    CALL SCMoutput(pptn_rate,'tot_precip_avg',                          &
         'Precipitation rate','kg/m2/s',                                &
         t_avg,d_sl,default_streams,'',RoutineName)

    CALL SCMoutput(accum_pptn,'tot_precip_acc',                         &
         'Accumulated Precipitation','kg/m2/s',                         &
         t_avg,d_sl,default_streams,'',RoutineName)

    !     Stash 5,216
    CALL SCMoutput(pptn_rate,'tot_precip',                              &
         'Total precipitation rate','kg/m2/day',                        &
         t_mult,d_sl,default_streams,'sec_day',RoutineName)

    CALL SCMoutput(pptn_rate,'tot_precip_inst',                         &
         'Total precipitation rate','kg/m2/s',                          &
         t_inst,d_sl,default_streams,'',RoutineName)

    !     Stash 5,214
    CALL SCMoutput(tot_rain,'tot_rain',                                 &
         'Total rainfall rate','kg/m2/day',                             &
         t_mult,d_sl,default_streams,'sec_day',RoutineName)

    CALL SCMoutput(tot_rain,'tot_rain_inst',                            &
         'Total rainfall rate','kg/m2/s',                               &
         t_inst,d_sl,default_streams,'',RoutineName)

    !     Stash 5,215
    CALL SCMoutput(tot_snow,'tot_snow',                                 &
         'Total snowfall rate','kg/m2/day',                             &
         t_mult,d_sl,default_streams,'sec_day',RoutineName)

    CALL SCMoutput(tot_snow,'tot_snow_inst',                            &
         'Total snowfall rate','kg/m2/s',                               &
         t_inst,d_sl,default_streams,'',RoutineName)

  END IF ! L_SCMDiags(SCMDiag_lsp) .OR. L_SCMDiags(SCMDiag_conv)
         ! .OR. L_SCMDiags(SCMDiag_gen)

  IF (.NOT. l_emcorr_opt) THEN
    l_emcorr = l_emcorr_tmp
  END IF

END IF ! model_type

! Deallocated here as needed later in the timestep than all the rest
! of the convection diags, which are deallocated by the call to
! cv_dealloc_diag_array.
IF (l_CallConvection .AND. Error_code == 0 ) DEALLOCATE( theta_incr_diag_conv )

IF ( L_print_L2norms ) THEN
  WRITE(umMessage,'(A, I2, A)') ' **  cycleno = ', cycleno              &
                        , ' **  L2 norms after atmos_physics2 **'
  CALL umPrint(umMessage,src='ATMOS_PHYSICS2')
  CALL atmos_phys2_norm(                                                &
                        norm_lev_start, norm_lev_end, n_rims_to_do,     &
                        theta_star, q_star, qcl_star, qcf_star,         &
                        qrain_star, qcf2_star, qgraup_star,             &
                        cf_star, cfl_star, cff_star,                    &
                        r_u, r_v, r_w, .TRUE., .FALSE., .FALSE. )
  IF ( L_tracer ) THEN
    WRITE(umMessage,'(A, I2, A)') ' **  cycleno = ', cycleno            &
                   ,' **  L2 norms of tracers after atmos_physics2  **'
    CALL umPrint(umMessage,src='atm_step_4A')
    CALL sl_tracer1_norm(                                               &
                         super_array_size,                              &
                         norm_lev_start, norm_lev_end, n_rims_to_do,    &
                         co2, aerosol, dust_div1, dust_div2,            &
                         dust_div3, dust_div4, dust_div5, dust_div6,    &
                         soot_new, soot_aged, soot_cld,                 &
                         bmass_new, bmass_aged, bmass_cld,              &
                         ocff_new, ocff_aged, ocff_cld,                 &
                         so2, so4_aitken, so4_accu, so4_diss,           &
                         nh3, nitr_acc, nitr_diss,                      &
                         tracer_ukca, tr_ukca, ozone_tracer,            &
                         .TRUE., .FALSE., .FALSE. )
 
    IF ( tr_vars > 0 )  CALL free_tracer_norm ( free_tracers, tr_vars,  &
                                          norm_lev_start, norm_lev_end, &
                                          n_rims_to_do,                 &
                                          .TRUE., .FALSE., .FALSE. )
 
    IF ( l_Sulpc_so2 .AND. l_sulpc_dms ) CALL dms_norm( dms,            &
                                          norm_lev_start, norm_lev_end, &
                                          n_rims_to_do,                 &
                                          .TRUE., .FALSE., .FALSE. )
  END IF  !  L_tracer
END IF !  L_print_L2norms

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE Atmos_Physics2
