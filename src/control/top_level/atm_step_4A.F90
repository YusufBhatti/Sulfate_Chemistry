! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Perform a 1-timestep integration of the Atmosphere Model
!
! Subroutine Interface:

SUBROUTINE atm_step_4A(                                               &
! River routing
 land, cumulus, nbdsc, ntdsc, ntml,                                   &
 ccb,cct,                                                             &
 g_p_field,                                                           &
 g_r_field,                                                           &
 obs_flag,obs,obs_flag_len,obs_len,                                   &
 ngrgas,grgas_addr)

USE sl_tracer1_mod
USE tr_reset_mod, ONLY: tr_reset_4A
USE TR_Set_Phys_mod, ONLY: tr_set_phys_4A
USE atm_fields_bounds_mod, ONLY : array_dims,                           &
                                  o3dims2,                              &
                                  pdims, pdims_s,                       &
                                  tdims, tdims_l, tdims_s,              &
                                  udims_l, udims_s, vdims_l, vdims_s,   &
                                  wdims, wdims_l, wdims_s
USE dump_headers_mod, ONLY: rh_z_top_theta, rh_tot_m_init,            &
                            rh_tot_mass_init, rh_tot_energy_init,     &
                            rh_energy_corr, a_inthd, a_realhd,        &
                            ih_stochastic_flag, ih_stph_n1,           &
                            ih_stph_n2, ih_sp_seed,                   &
                            allocate_sp_coefficients,                 &
                            assign_a_flddepc_input_values
USE rad_input_mod
USE jules_sea_seaice_mod, ONLY: l_ctile
USE cv_cntl_mod, ONLY: lcv_3d_ccw
USE level_heights_mod
USE trignometric_mod
USE cderived_mod, ONLY: delta_lambda, delta_phi, base_lambda,  &
                        base_phi, lat_rot_NP, long_rot_NP
USE dyn_coriolis_mod
USE dyn_var_res_mod
USE diff_coeff_mod
USE turb_diff_mod
USE metric_terms_mod
USE eg_calc_p_star_mod
USE rad_mask_trop_mod
USE rot_coeff_mod
USE swapable_field_mod, ONLY: swapable_field_pointer_type
USE timestep_mod
USE atm_step_local
USE atm_step_timestep_mod, ONLY: atm_step_info
! Logicals and stash diagnostic types needed for stochastic physics
USE stochastic_physics_run_mod,  ONLY:                                  &
    l_skeb2, l_rp2, l_spt,                                              &
    i_rp_scheme, i_rp2b, lai_mult_rp,                                   &
    rhcrit_max, rhcrit_min, m_ci, m_ci_max, m_ci_min,                   &
    par_mezcla_max, par_mezcla_min, par_mezcla,                         &
    g0_rp, g0_rp_max, g0_rp_min, charnock_max, charnock_min,            &
    ricrit_rp, ricrit_rp_max, ricrit_rp_min,                            &
    lambda_min_rp, lambda_min_rp_max, lambda_min_rp_min,                &
    a_ent_1_rp, a_ent_1_rp_max, a_ent_1_rp_min,                         &
    g1_rp, g1_rp_max, g1_rp_min, i_pert_theta,                          &
    stph_seed_present

USE jules_snow_mod, ONLY: nsmax
USE arcl_mod, ONLY: npd_arcl_species, npd_arcl_compnts
USE iau_mod, ONLY: l_iau, iau_firstcallts, iau_lastcallts

USE planet_constants_mod, ONLY: g, planet_radius, recip_epsilon
USE planet_suite_mod, ONLY: nsteps_consv_print

USE g_wave_input_mod
USE bl_option_mod, ONLY: l_quick_ap2, off, i_bl_vn, i_bl_vn_1a, i_bl_vn_0
USE mym_option_mod, ONLY: l_adv_turb_field
USE jules_sea_seaice_mod, ONLY: nice, nice_use, charnock
USE submodel_mod, ONLY: atmos_sm, atmos_im
USE stash_array_mod, ONLY: stash_maxlen, sf
USE jules_surface_types_mod, ONLY: ntype, npft

USE water_constants_mod


USE ukca_cdnc_mod, ONLY: ukca_cdnc, cdnc_dim1, cdnc_dim2, cdnc_dim3
USE ukca_option_mod, ONLY: l_ukca, l_ukca_chem, l_ukca_mode,      &
                     l_ukca_set_trace_gases, l_ukca_strat,        &
                     l_ukca_strattrop, l_ukca_prescribech4,       &
                     l_ukca_trop, l_ukca_tropisop, l_ukca_raq,    &
                     l_ukca, l_ukca_aie1, l_ukca_aie2,            &
                     l_ukca_radaer, l_ukca_radaer_sustrat,        &
                     l_ukca_offline, l_conserve_ukca_with_tr,     &
                     l_ukca_raqaero, l_ukca_plume_scav,           &
                     l_ukca_h2o_feedback
USE glomap_clim_option_mod, ONLY: &
    l_glomap_mode_clim,           &
    l_glomap_clim_arg_act,        &
    l_glomap_clim_aie1,           &
    l_glomap_clim_aie2,           &
    l_glomap_clim_radaer,         &
    l_glomap_clim_radaer_sustrat, &
    i_glomap_clim_setup,          &
    i_gc_sussocbc_5mode
USE ukca_main1_mod, ONLY: ukca_main1
USE ukca_d1_defs,   ONLY: item1_plume_diags, nmax_plume_diags,    &
                          l_ukca_plume_diags

USE ukca_radaer_get_mod,        ONLY: ukca_radaer_get
USE ukca_radaer_init_mod,       ONLY: ukca_radaer_init
USE ukca_radaer_saved_mod,      ONLY: ukca_radaer
USE ukca_mode_setup,            ONLY: ukca_mode_sussbcoc_5mode
USE allocate_ukca_cdnc_mod,     ONLY: allocate_ukca_cdnc
USE glomap_clim_radaer_get_mod, ONLY: glomap_clim_radaer_get

! If using Priestley Conservation scheme
USE eg_correct_tracers_priestley_mod
USE eg_correct_moisture_priestley_mod
USE eg_correct_tracers_ukca_mod
USE eg_correct_thetav_priestley_mod

!      USE ancil_info, ONLY: nsmax

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

! deallocate FV-TRACK arrays
USE track_mod, ONLY: destroy_track_arr

USE eg_destroy_vert_damp_mod, ONLY: eg_destroy_vert_damp
USE eg_dry_static_adj_mod, ONLY: eg_dry_static_adj
USE um_types
USE eg_helmholtz_mod

USE eg_idl_forcing_mod
USE external_force_2_mod, ONLY: external_force_2
USE eg_r_mod
USE eg_r_s_mod
USE eg_f1sp_mod
USE eg_f1sp_inc_mod
USE eg_thetav_theta_mod
USE eg_sisl_resetcon_mod
USE eg_set_helm_lhs_mod
USE eg_sisl_init_mod
USE eg_q_to_mix_mod
USE eg_sl_helmholtz_mod
USE eg_sl_helmholtz_inc_mod
USE eg_sl_full_wind_mod
USE eg_sl_rho_mod
USE eg_rho_pseudo_lbflux_mod, ONLY: eg_rho_pseudo_lbflux
USE eg_sl_thermo_mod
USE eg_sl_moisture_mod
USE eg_sl_turb_mod, ONLY: eg_sl_turb
USE eg_sl_casim_mod, ONLY: eg_sl_casim
USE eg_moisture_pseudo_lbflux_mod, ONLY: eg_moisture_pseudo_lbflux
USE eg_idl_set_init_mod
USE eg_coriolis_star_mod
USE eg_dx_diags_mod
USE eg_set_helmholtz_mod
USE eg_mix_to_q_mod
USE eg_v_at_poles_mod
USE eg_NI_filter_Ctl_mod
USE eg_NI_filter_incs_Ctl_mod
USE eg_dxout_mod
USE horiz_grid_mod
USE ref_pro_mod

USE eg_total_mass_region_mod, ONLY: eg_total_mass_region

USE eg_balance_lbc_values_mod

USE eg_alpha_mod
USE eg_alpha_ramp_mod
USE set_metric_terms_4A_mod

USE Control_Max_Sizes
USE Atmos_Max_Sizes, ONLY: model_levels_max
USE UM_ParVars
USE UM_ParCore,      ONLY: mype, nproc
USE UM_ParParams,    ONLY: halo_type_extended, halo_type_no_halo
USE Field_Types,     ONLY: fld_type_p, fld_type_r, fld_type_u, fld_type_v, &
                           fld_type_w
USE ereport_mod, ONLY: ereport
USE umPrintMgr
USE set_star_zero_level_mod

USE departure_pts_mod
USE fields_rhs_mod
USE coriolis_mod
USE gravity_mod
USE eg_parameters_mod

USE print_diag_mod

USE rimtypes
USE lbc_mod
USE eg_conserv_moist_mod
USE eg_conserv_tracers_mod
USE halo_exchange, ONLY: swap_bounds
USE mpp_conf_mod,  ONLY: swap_field_is_vector, swap_field_is_scalar

USE o3intp_mod, ONLY: io3_3dspec, io3_2dspec, io3_2dmasscon, io3_trop_map,  &
                      io3_trop_map_masscon
USE problem_mod, ONLY: standard, dynamical_core, idealised_problem,         &
                       idealised_planet
USE eng_mass_diag_mod, ONLY: eng_mass_diag
USE wet_to_dry_n_calc_mod

USE eg_set_adv_winds_mod

USE physics_tendencies_mod, ONLY:                                 &
    init_stph_tendencies, destroy_phy_tend,                       &
    l_retain_stph_tendencies


!Dynamic tracer schemes
USE dyn_tracers_mod, ONLY: dyn_tr_sav, dyn_tr_slow,               &
                           eg_sl_dyn_tr, dyn_tr_fast,             &
                           dyn_tr_np1, dyn_tr_end

! Stochastic Physics
USE stph_rp2_mod,   ONLY: stph_rp2
USE stph_setup_mod, ONLY: stph_setup
USE stph_seed_mod, ONLY: stph_seed_copy_from_dump, stph_seed_copy_to_dump
USE legendre_poly_comp_stph_mod, ONLY: legendre_poly_comp_stph,   &
                                        create_legendre, destroy_legendre
USE stph_skeb2_mod, ONLY: stph_skeb2
USE spt_main_mod, ONLY: spt_main
USE spt_add_mod, ONLY: spt_add
USE for_pattern_mod, ONLY: for_pattern
USE eg_total_conservation_mod
USE dynamics_input_mod, ONLY:                                     &
          conserve_dry_mass, not_conserved, constant_factor,      &
          linear_factor, linear_factor_IE,                        &
          L_eg_dry_static_adj, tol_sc_fact,                       &
          T_surf_ref,p_surf_ref,NumCycles,L_new_tdisc,            &
          L_LBC_balance,L_qwaterload,                             &
          GCR_use_residual_Tol,GCR_diagnostics,GCR_precon_option, &
          GCR_tol_res,GCR_tol_abs,L_lbc_old,                      &
          eg_vert_damp_coeff,L_RK_dps,L_inc_solver, Ih,           &
          L_conserv_smooth_lap,L_eliminate_rho,L_init_Fnm1,       &
          eg_vert_damp_profile, eta_s,                            &
          innits,l_accel_convergence,                             &
          l_viscosity, horiz_viscosity, vert_viscosity

USE lam_config_inputs_mod, ONLY: n_rims_to_do

USE dynamics_testing_mod, ONLY:                                   &
          L_Physics,L_Run_With_Physics2,                          &
          L_Backwards, L_dry, L_idealised_data,                   &
          L_prognostic_level0, problem_number

USE idealise_run_mod, ONLY:                                       &
          l_fixed_lbcs, l_force_lbc,                              &
          l_shallow, l_const_grav,                                &
          tstep_plot_start, tstep_plot_frequency, L_spec_z0,      &
          t_surface

USE sl_input_mod, ONLY:                                           &
          L_conserve_tracers,halo_lam,halo_phi,look_lam,look_phi, &
          recip_dlam,recip_dphi,L_conserv,L_mono,L_high,          &
          L_Ritchie_high,L_Ritchie_mono,L_2d_sl_geometry,         &
          L_sl_halo_reprod,high_order_scheme,monotone_scheme,     &
          Ritchie_high_order_scheme,Ritchie_monotone_scheme,      &
          Depart_scheme,Depart_order,interp_vertical_search_tol,  &
          Theta_SL,moist_SL,Wind_SL,rho_SL,tracer_SL,             &
          L_priestley_correct_moist, L_priestley_correct_tracers, &
          tr_priestley_opt, L_priestley_correct_thetav

USE diag_ctl_mod, ONLY:                                           &
    rpemax, rpemin, ipesum, rpesum,                               &
    max_w_run, max_wind_run, min_theta1_run,                      &
    dtheta1_run, max_div_run, min_div_run,                        &
    min_lapse_run, max_shear_run, time_max_shear,                 &
    time_div_max, time_div_min, time_lapse_min,                   &
    time_w_max, time_max_wind, time_theta1_min,                   &
    max_KE_run, min_KE_run, max_noise_run,                        &
    time_KE_max, time_KE_min, time_noise_max

! Aerosols
USE aero_ctl_mod_4A,          ONLY: aero_ctl_4A
USE get_sulpc_oxidants_mod,   ONLY: get_sulpc_oxidants
USE set_arcl_clim_mod,        ONLY: set_arcl_clim
USE set_arcl_dimensions_mod,  ONLY: set_arcl_dimensions
USE write_sulpc_oxidants_mod, ONLY: write_sulpc_oxidants

USE dust_parameters_mod, ONLY:  l_dust,                              &
     l_dust_div1,       l_dust_div2,      l_dust_div3,               &
     l_dust_div4,       l_dust_div5,      l_dust_div6,               &
     l_dust_div1_lbc,   l_dust_div2_lbc,  l_dust_div3_lbc,           &
     l_dust_div4_lbc,   l_dust_div5_lbc,  l_dust_div6_lbc,           &
     l_twobin_dust

USE carbon_options_mod, ONLY: l_co2_interactive
USE gen_phys_inputs_mod, ONLY: l_mr_physics
USE lbc_read_data_mod, ONLY: l_int_uvw_lbc, rimweightsa
USE run_aerosol_mod, ONLY: l_ocff_new, l_ocff_agd, l_ocff_cld,  &
                           l_soot_new, l_soot_agd, l_soot_cld,  &
                       l_bmass_new, l_bmass_agd,  l_bmass_cld,  &
         l_sulpc_so2, l_so2_surfem, l_so2_hilem, l_so2_natem,  &
         l_sulpc_ozone, l_sulpc_online_oxidants, l_sulpc_2_way_coupling, &
         l_sulpc_so2_o3_nonbuffered, l_sulpc_nh3, l_nh3, l_nh3_em,  &
         l_soot, l_soot_surem, l_soot_hilem,  &
         l_sulpc_dms, l_dms, l_dms_em, l_dms_em_inter, i_dms_flux,  &
         l_so2, l_so4_aitken, l_so4_accu, l_so4_diss,  &
         l_biomass, l_bmass_surem, l_bmass_hilem,  &
         l_ocff, l_ocff_surem, l_ocff_hilem,  &
         l_tracer1_non_hydro, call_chem_freq,  &
         l_nitr_acc, l_nitr_diss, l_nitrate, l_use_nitrate_sulpc,  &
         l_use_sulphate_sulpc, l_use_seasalt_sulpc, l_use_seasalt_pm,  &
         l_use_ocff_sulpc, l_use_bmass_sulpc,  &
         l_so2_lbc, l_dms_lbc, l_so4_aitken_lbc, l_so4_accu_lbc,  &
         l_so4_diss_lbc, l_nh3_lbc, l_soot_new_lbc, l_soot_agd_lbc,  &
         l_soot_cld_lbc, l_bmass_new_lbc, l_bmass_agd_lbc,  &
         l_bmass_cld_lbc, l_ocff_new_lbc, l_ocff_agd_lbc,  &
         l_ocff_cld_lbc, l_nitr_acc_lbc, l_nitr_diss_lbc

USE easyaerosol_option_mod, ONLY: l_easyaerosol_sw, l_easyaerosol_lw, &
                                  l_easyaerosol_cdnc,                 &
                                  l_easyaerosol_autoconv
USE easyaerosol_read_input_mod, ONLY: easyaerosol_read_input_ctl
USE def_easyaerosol, ONLY: allocate_easyaerosol_rad, &
                           allocate_easyaerosol_cdnc, &
                           deallocate_easyaerosol_rad, &
                           deallocate_easyaerosol_cdnc
USE spec_sw_lw, ONLY: sw_spectrum, lw_spectrum

USE mphys_inputs_mod, ONLY:                                          &
     l_mcr_qcf2,        l_mcr_qcf2_lbc,   l_mcr_qgraup,              &
     l_mcr_qgraup_lbc,  l_mcr_qrain_lbc,  l_mcr_qrain,               &
     l_use_seasalt_autoconv, l_subgrid_qcl_mp, l_casim

USE casim_switches, ONLY: n_casim_progs

USE cloud_inputs_mod, ONLY:                                          &
     i_cld_area, i_cld_vn, l_pc2_reset, l_pc2_lbc, rhcrit
USE pc2_constants_mod, ONLY: acf_off, acf_brooks, i_cld_off,         &
     i_cld_smith, i_cld_pc2
USE river_inputs_mod, ONLY: river_vel, river_mcoef, i_river_vn

USE eng_corr_inputs_mod, ONLY:                                       &
     l_emcorr,          lmass_corr,       lemq_print,                &
     a_energysteps,     lqt_corr,         lenergy

USE nudging_input_mod, ONLY: l_nudging

USE murk_inputs_mod, ONLY: l_murk, l_murk_advect, l_murk_lbc

USE cosp_input_mod, ONLY: l_cosp

USE nlstcall_mod, ONLY: ldump, lexit, ltimer, lcal360, lstashdumptimer

USE science_fixes_mod, ONLY: l_fix_conserv, l_fix_iau_rim_density

USE update_rho_mod
USE update_moisture_fields_mod
USE eg_star_mod
USE init_etadot_mod
USE filter_diag_printing_mod
USE diag_R_mod
USE windmax_mod
USE conservation_diag_mod, ONLY: print_conservation_diag,                  &
                                 initial_timestep, not_initial_timestep
USE init_psiw_mod
USE free_tracers_inputs_mod, ONLY: a_max_trvars, a_tr_active_lbc_index, &
     a_tracer_last, a_tracer_first, l_pv_tracer, l_diab_tracer
USE run_info, ONLY: time_start_atmstep
USE io_configuration_mod, ONLY: print_runtime_info

USE nlsizes_namelist_mod, ONLY:                                            &
    a_len1_coldepc, a_len1_flddepc, a_len1_levdepc, a_len1_rowdepc,        &
    a_len2_coldepc, a_len2_flddepc, a_len2_levdepc, a_len2_lookup,         &
    a_len2_rowdepc, a_len_cfi1, a_len_cfi2, a_len_cfi3, a_len_extcnst,     &
    a_len_inthd, a_len_realhd, aocpl_p_rows, aocpl_row_length, bl_levels,  &
    cloud_levels, global_row_length, global_rows, land_field,              &
    len1_lbc_comp_lookup, len1_lookup,                                     &
    len_dumphist, len_fixhd, len_tot, model_levels, mpp_len1_lookup,       &
    n_cca_lev, n_obj_d1_max, n_rows, ntiles, ozone_levels,                 &
    river_row_length, river_rows, row_length, rows, sm_levels, st_levels,  &
    theta_field_size, theta_off_size, tr_lbc_ukca, tr_lbc_vars, tr_levels, &
    tr_ukca, tr_vars, super_array_size, moisture_array_size,               &
    turb_array_size

USE model_time_mod, ONLY:                                             &
    i_day, i_day_number, i_hour, i_minute, i_month, i_second, i_year, &
    previous_time, stepim
USE errormessagelength_mod, ONLY: errormessagelength

USE model_domain_mod, ONLY: l_regular, model_type, mt_global, mt_lam,      &
                            mt_bi_cyclic_lam

! Subroutines for idealised model stash
USE alloc_ideal_diag_mod,   ONLY: alloc_ideal_diag
USE dealloc_ideal_diag_mod, ONLY: dealloc_ideal_diag
USE diagnostics_ideal_mod,  ONLY: diagnostics_ideal

USE um_stashcode_mod,   ONLY: &
    stashcode_pws_sec,        &
    stashcode_glomap_sec,     &
    stashcode_glomap_clim_sec

USE um_stashcode_mod,     ONLY: stashcode_pws_sec, stashcode_glomap_sec
USE pws_diags_mod,        ONLY: flag_upd_helicity_5k
USE pws_diags_driver_mod, ONLY: pws_diags_driver

USE atm_fields_real_mod

USE atm_boundary_headers_mod, ONLY: rim_stepsa

USE land_soil_dimensions_mod, ONLY: land_points, land_ice_points,       &
     soil_points, land_index, land_ice_index, soil_index
USE disturb_veg_category_mod, ONLY: disturb_veg_pointer

USE river_routing_sizes_mod, ONLY: xpa,xua,xva,                             &
                                   ypa,yua,yva

USE cal_eng_mass_corr_4a_mod, ONLY: cal_eng_mass_corr_4a
USE ls_acf_brooks_mod, ONLY: ls_acf_brooks
USE pc2_pressure_forcing_mod, ONLY: pc2_pressure_forcing
USE qt_bal_cld_mod, ONLY: qt_bal_cld
USE nudging_main1_mod, ONLY: nudging_main1

USE diagnostics_adv_mod, ONLY: diagnostics_adv
USE adv_increments_mod, ONLY: adv_incs_init, adv_incs_dealloc

USE diff_increments_mod, ONLY: diff_incs_init, diff_incs_dealloc

USE diagnostics_solver_mod, ONLY: diagnostics_solver
USE solver_increments_mod, ONLY:                                              &
    solv_incs_init, solv_incs_calc, solv_incs_dealloc

USE diag_adv_correct_mod, ONLY: diagnostics_adv_correct
USE adv_correct_incs_mod, ONLY:                                               &
    adv_correct_incs_init, adv_correct_incs_calc, adv_correct_incs_dealloc

USE eot_increments_mod, ONLY:                                                 &
    eot_incs_init, eot_incs_calc, eot_incs_dealloc, eot_inc_rho

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

!
! Description: Perform a 1-timestep integration of the Atmosphere Model,
!   including assimilation, physics and dynamics processing.
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Declarations:

INTEGER ::                                                        &
  g_p_field                                                       &
                        ! IN : global horiz domain for atmos
, g_r_field             ! IN : global horiz domain for rivers
INTEGER :: obs_flag_len,obs_len
INTEGER :: obs_flag(obs_flag_len)
REAL    :: obs(obs_len)

! ENDGame prognostic variables (not included in the start dump)

REAL :: theta_star(tdims%i_start:tdims%i_end,                           &
                   tdims%j_start:tdims%j_end,                           &
                   tdims%k_start:tdims%k_end),                          &
!       temorarily stored atmos_phys2 initial state for ENDGame
!        source term computation in eg_R_S:
   theta_star_n(tdims%i_start:tdims%i_end,                              &
                tdims%j_start:tdims%j_end,                              &
                tdims%k_start:tdims%k_end)

REAL :: exner_prime_term(pdims_s%i_start:pdims_s%i_end,                 &
                         pdims_s%j_start:pdims_s%j_end,                 &
                         pdims_s%k_start:pdims_s%k_end)

REAL ::                                                                 &
 exner_star(pdims_s%i_start:pdims_s%i_end,                              &
            pdims_s%j_start:pdims_s%j_end,                              &
            pdims_s%k_start:pdims_s%k_end),                             &
 frac_control(land_field,ntype)   !Forcing for land surface (3C)


REAL ::                                                                 &
  rho_n (tdims_s%i_start:tdims_s%i_end,                                 &
         tdims_s%j_start:tdims_s%j_end,                                 &
         tdims_s%k_start:tdims_s%k_end)

REAL :: biogenic (tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end,                            &
                              1:tdims%k_end)

! Subroutine arguments:

! Land mask
LOGICAL :: land   (row_length,rows)
LOGICAL :: cumulus(row_length,rows)
! Top level for turb mixing in any decoupled Sc layer
INTEGER :: ntdsc(row_length,rows)
! Bottom level for turb mixing in any decoupled Sc layer
INTEGER :: nbdsc(row_length,rows)
! Number of model levels in the  turbulently mixed layer
INTEGER :: ntml(row_length,rows)

INTEGER :: ccb(row_length,rows)
INTEGER :: cct(row_length,rows)

! 3-D fields of species to be passed down to radiation
INTEGER, INTENT(IN) :: ngrgas
INTEGER, INTENT(IN) :: grgas_addr(ngrgas)

! Local parameters:
CHARACTER(LEN=*) :: RoutineName
PARAMETER (   RoutineName='ATM_STEP_4A')

INTEGER, SAVE :: dx_frame_cnt
REAL :: exner_res,exner_res_tmp

REAL           :: total_rho, mass_fix_factor,                           &
                  total_gr, total_gr_r, total_rho_r,                    &
                  mass_fix_A, mass_fix_B, total_rho_n
REAL           :: total_aam, total_ke

INTEGER ::                                                              &
 errorstatus      ! Return code : 0 Normal Exit : >0 Error

CHARACTER(LEN=errormessagelength) ::                                    &
 cmessage         ! Error message if return code >0

LOGICAL, PARAMETER :: l_test_tracer=.FALSE.

LOGICAL :: l_call_from_solver

! Used for checking if at end of integration
LOGICAL :: lexitnow

INTEGER, SAVE :: timestep_number_chainsafe

! Local arrays for using the aerosol climatology for NWP

! Internal model switches
LOGICAL :: l_use_arcl(npd_arcl_species)

! Array index of each component
INTEGER :: i_arcl_compnts(npd_arcl_compnts)

INTEGER :: cloud_tol

REAL ::                                                           &
  mag_vector_np (model_levels)                                    &
, dir_vector_np (model_levels)                                    &
, mag_vector_sp (model_levels)                                    &
, dir_vector_sp (model_levels)                                    &
, lambda_a (row_length) ! delta_lambda term for polar wind

! LAM LBC tendency

REAL ::                                                            &
  w_lbc_real_tend(lenrima(fld_type_p,halo_type_extended,           &
                    rima_type_norm),0:model_levels)                &
, exner_lbc_real_tend(lenrima(fld_type_p,halo_type_extended,       &
                    rima_type_norm),model_levels+1)

! Physics arrays needed by dynamics
REAL ::                                                            &
  wet_to_dry_np1(tdims_s%i_start:tdims_s%i_end,                    &
                 tdims_s%j_start:tdims_s%j_end,                    &
                 tdims_s%k_start:tdims_s%k_end)

! arrays holding information to be passed between physics
! routines.

REAL ::                                                            &
  ls_rain    (row_length, rows)                                    &
, ls_snow    (row_length, rows)                                    &
, ls_graup    (row_length, rows)                                   &
, micro_tends(row_length, rows, 2, bl_levels)                      &
!                          ! Tendencies from microphys within BL levels
!                          ! (TL, K/s; QW, kg/kg/s)
, ls_rainfrac(land_points)
                           ! Rain fraction array on land points to be
                           ! passed to atmos_physics2

! height of lcl in a well-mixed BL (types 3 or 4), 0 otherwise
REAL :: zlcl_mixed(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

! Radiation fields 1. SW & common with LW.
REAL ::                                                            &
  photosynth_act_rad(row_length, rows)                             &
                                       ! Net downward
!                                 shortwave radiation in band 1 (w/m2).
, rad_hr(row_length, rows, 2, bl_levels)                           &
!                                 BL (LW,SW) rad heating rates
, dolr(row_length,rows)                                            &
!       local field "dolr" is distinguished from "dolr_field"
!       (in atm_fields_mod)
                             ! TOA - surface upward LW
, sw_tile(land_field,ntiles)                                       &
                             ! Surface net SW on land tiles
, cos_zenith_angle(row_length, rows)

! MPP-related Arrays
INTEGER ::                                                         &
 g_row_length(0:nproc-1)                                           &
                         ! Table of number of points on a row
,g_rows(0:nproc-1)                                                 &
                         ! Table number of rows in theta field
,g_i_pe(1-halo_i:global_row_length+halo_i)                         &
               ! processor on my
!               processor-row holding a given value in i direction
,g_j_pe(1-halo_j:global_rows+halo_j) ! processor on my
!               processor-row holding a given value in j direction

! Useful diagnostic
INTEGER :: istat
REAL :: w_max, u_max, u_at_w, v_at_w

REAL :: stashwork0_dummy(1) ! STASHwork not defined for section 0,
                         !  but required as dummy argument.

! Local arrays for phys1 and phys2 increments for tracers:
REAL :: super_tracer_phys1(tdims_l%i_start:tdims_l%i_end,          &
                           tdims_l%j_start:tdims_l%j_end,          &
                           tdims_l%k_start:tdims_l%k_end,          &
                           super_array_size)

REAL :: super_tracer_phys2(tdims%i_start:tdims%i_end,              &
                           tdims%j_start:tdims%j_end,              &
                           tdims%k_start:tdims%k_end,              &
                           super_array_size)

REAL :: gs1(land_field)

! local variable
REAL ::                                                            &
  tot_dry_mass_final                                               &
                      ! mass at end of energy correction period
, tot_energy_final                                                 &
                      ! energy at end of energy correction period
, tot_moist_final                                                  &
                      ! moist at end of energy correction period
, energy_corr_now
                      ! instanteous energy correction

REAL ::                                                            &
  flux_e(row_length, rows)                                         &
                           ! Surface latent heat flux (W/m^2)
, flux_h(row_length, rows)                                         &
                           ! Surface sensible heat flux (W/m^2)
, ustar_in(row_length, rows)                                       &
                           ! Surface friction velocity (m/s)
, z0m_scm(row_length, rows)                                        &
                           ! SCM specified z0m (m)
, z0h_scm(row_length, rows)! SCM specified z0m (m)

! Dummy scm variables:
INTEGER :: conv_mode
LOGICAL :: l_emcorr_opt

REAL, ALLOCATABLE ::                             &
  q_inc_subs(:,:,:), th_inc_subs(:,:,:)                       &
                            ! subsidence increments
, q_inc_ls(:,:,:), th_inc_ls (:,:,:)                          &
                            ! large scale increments
, u_inc_dmp(:,:,:), q_inc_dmp(:,:,:), th_inc_dmp  (:,:,:)     &
                            !Damping incs
, v_inc_dmp(:,:,:)
! Tolerance for CycleNo >1
REAL ::                                                            &
  gcr_run_tol_abs                                                  &
, gcr_run_tol_res

TYPE(swapable_field_pointer_type) :: fields_to_swap(7)  ! multivariate
                                                        ! swapbounds

INTEGER :: exppxi               ! Function to extract ppxref info

! Local parameters for mixing ratio physics
! Mixing ratios for atmos_physics1 and 2 are defined through namelist
LOGICAL, PARAMETER :: l_mr_qtbalcld=.FALSE. ! Use mr's for qt_bal_cld
LOGICAL, PARAMETER :: l_mr_iau     =.FALSE. ! Use mr's for tfilt_ctl
LOGICAL, PARAMETER :: l_mr_pc2     =.FALSE. ! Use mr's for PC2 routines

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

REAL :: time_end_atmstep

INTEGER :: ierr

LOGICAL :: l_call_from_f1sp_in  ! is the call from f1sp routine?

! Input variables to be used by COSP (row_length,rows,model_levels)
!  Convective rainfall
REAL,POINTER,SAVE :: cosp_crain_3d(:,:,:)
!  Convective snowfall
REAL,POINTER,SAVE :: cosp_csnow_3d(:,:,:)

INTEGER :: lbc_size_new
INTEGER :: fdc_size, spt_size, skeb2_size

TYPE (array_dims) twoddims
TYPE (array_dims) twoddims_no_halo

LOGICAL, SAVE :: subseq_to_dump = .FALSE.
INTEGER, PARAMETER :: number_qs = 6
REAL               :: pseudo_lbflux_mass
REAL               :: pseudo_lbflux_moisture(number_qs)

INTEGER :: IS, ie, js, je

! last item number of UKCA plume scavenging diagnostics
INTEGER :: item_plume_diags_last

! Values of leaf area index (lai_pft) passed into subroutines
! If the RP2b scheme is switched on, then 
! lai_pft_in = lai_mult_rp*lai_pft.
! If the RP2b scheme is switched off, then 
! lai_pft_in = lai_pft
REAL :: lai_pft_in(land_field, npft)

!- End of header

! ----------------------------------------------------------------------
! Section 0.  Initialisation.
! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

IF (ltimer) CALL timer('Atm_Step_4A (AS)',5)

twoddims=tdims_s
twoddims%k_end = twoddims%k_start
twoddims_no_halo=tdims
twoddims_no_halo%k_end = twoddims%k_start

alpha_changed = .FALSE.
reference_profile_changed = .FALSE.

IF ( gcr_use_residual_tol ) THEN
  solver_tolerance = gcr_tol_res*(tol_sc_fact**(numcycles*innits-1))
ELSE
  solver_tolerance = gcr_tol_abs*(tol_sc_fact**(numcycles*innits-1))
END IF

IF (printstatus > prstatus_oper) THEN
  IF (.NOT. couple_app) CALL umPrint('coupling disabled!',src='atm_step_4A')
END IF

errorstatus = 0

!     ENDGame NEEDS the np1 variables
l_new_tdisc = .TRUE.

! DEPENDS ON: Atm_Step_Init
CALL atm_step_init (                                                  &
     lambda_a, g_row_length, g_rows, g_i_pe, g_j_pe, flux_e, flux_h,  &
     ustar_in, z0m_scm, z0h_scm, errorstatus )

check_bottom_levels = MIN(check_bottom_levels, model_levels)
interp_vertical_search_tol = MIN(interp_vertical_search_tol,          &
                                 model_levels/2)

! swapbound the exner from the dump!
IF (first_atmstep_call) THEN
  CALL swap_bounds(exner,                                                &
                   wdims_s%i_len - 2*wdims_s%halo_i,                     &
                   wdims_s%j_len - 2*wdims_s%halo_j,                     &
                   wdims_s%k_len,                                        &
                   wdims_s%halo_i, wdims_s%halo_j,                       &
                   fld_type_p,swap_field_is_scalar)
  CALL swap_bounds(etadot,                                               &
                   tdims_s%i_len - 2*tdims_s%halo_i,                     &
                   tdims_s%j_len - 2*tdims_s%halo_j,                     &
                   tdims_s%k_len,                                        &
                   tdims_s%halo_i, tdims_s%halo_j,                       &
                   fld_type_p,swap_field_is_scalar)
END IF

IF (first_atmstep_call) THEN

   ! allocate storage
  CALL init_fields_rhs(l_skeb2,l_spt)

  !  rho elimination
  alpha_p = alpha_rho
  tau_p   = tau_rho

! We do not deallocate the np1 fields until the very end of the run. To avoid
! multiple allocation we have to if-test with timestep_number == 1
!
! DEPENDS ON: Atm_Step_alloc_4A
  CALL atm_step_alloc_4A(                                              &
    cf_star, cfl_star, cff_star, q_star, qcl_star, qcf_star,           &
    frac_control, r_u, r_v, r_w, errorstatus, 'newtdisc')

!----------------------------------------------------------------------
! All new runs for idealised setups
!----------------------------------------------------------------------
  IF (l_idealised_data) THEN
    ! Allocate arrays to hold initial information for diagnostic output
    CALL alloc_ideal_diag(0)
  END IF
END IF  ! first_atmstep_call

! ---------------------------------------------------------------------
! Section 0.1  Initialisation for idealised test problems
!              For standard runs go to end section 0.1
! ---------------------------------------------------------------------

IF ( first_atmstep_call ) THEN

  CALL eg_idl_set_init(                                                 &
                      row_length, rows, n_rows, halo_i, halo_j,         &
                      offx, offy, model_levels,                         &
                      datastart, global_row_length, global_rows,        &
                      nproc, g_row_length, g_rows,                      &
                      u, v, w, u_adv, v_adv, w_adv, dryrho, exner,      &
                      exner_surf,                                       &
                      m_v_np1, m_cl_np1, m_cf_np1,                      &
                      m_r_np1, m_gr_np1, m_cf2_np1, exner_surf_np1,     &
                      a_realhd(rh_z_top_theta),                         &
                      l_RK_dps, l_dry, alpha_w, ih,                     &
                      etadot, psi_w_surf, psi_w_lid,                    &
                      thetav, m_v, m_cl, m_cf, m_r, m_gr, m_cf2)

  CALL swap_bounds(m_v,                                                  &
                   tdims_s%i_len - 2*tdims_s%halo_i,                     &
                   tdims_s%j_len - 2*tdims_s%halo_j,                     &
                   tdims_s%k_len,                                        &
                   tdims_s%halo_i, tdims_s%halo_j,                       &
                   fld_type_p,swap_field_is_scalar)
  CALL swap_bounds(m_cl,                                                 &
                   tdims_s%i_len - 2*tdims_s%halo_i,                     &
                   tdims_s%j_len - 2*tdims_s%halo_j,                     &
                   tdims_s%k_len,                                        &
                   tdims_s%halo_i, tdims_s%halo_j,                       &
                   fld_type_p,swap_field_is_scalar)
  CALL swap_bounds(m_cf,                                                 &
                   tdims_s%i_len - 2*tdims_s%halo_i,                     &
                   tdims_s%j_len - 2*tdims_s%halo_j,                     &
                   tdims_s%k_len,                                        &
                   tdims_s%halo_i, tdims_s%halo_j,                       &
                   fld_type_p,swap_field_is_scalar)
  CALL swap_bounds(m_cf2,                                                &
                   tdims_s%i_len - 2*tdims_s%halo_i,                     &
                   tdims_s%j_len - 2*tdims_s%halo_j,                     &
                   tdims_s%k_len,                                        &
                   tdims_s%halo_i, tdims_s%halo_j,                       &
                   fld_type_p,swap_field_is_scalar)
  CALL swap_bounds(m_r,                                                  &
                   tdims_s%i_len - 2*tdims_s%halo_i,                     &
                   tdims_s%j_len - 2*tdims_s%halo_j,                     &
                   tdims_s%k_len,                                        &
                   tdims_s%halo_i, tdims_s%halo_j,                       &
                   fld_type_p,swap_field_is_scalar)
  CALL swap_bounds(m_gr,                                                 &
                   tdims_s%i_len - 2*tdims_s%halo_i,                     &
                   tdims_s%j_len - 2*tdims_s%halo_j,                     &
                   tdims_s%k_len,                                        &
                   tdims_s%halo_i, tdims_s%halo_j,                       &
                   fld_type_p,swap_field_is_scalar)
  CALL swap_bounds(exner_surf, tdims%i_len, tdims%j_len,                 &
                   1, tdims_s%halo_i,                                    &
                   tdims_s%halo_j,                                       &
                   fld_type_p,swap_field_is_scalar)
  CALL swap_bounds(exner,                                                &
                   wdims_s%i_len - 2*wdims_s%halo_i,                     &
                   wdims_s%j_len - 2*wdims_s%halo_j,                     &
                   wdims_s%k_len,                                        &
                   wdims_s%halo_i, wdims_s%halo_j,                       &
                   fld_type_p,swap_field_is_scalar)
  CALL swap_bounds(dryrho,                                               &
                   pdims_s%i_len - 2*pdims_s%halo_i,                     &
                   pdims_s%j_len - 2*pdims_s%halo_j,                     &
                   pdims_s%k_len,                                        &
                   pdims_s%halo_i, pdims_s%halo_j,                       &
                   fld_type_p,swap_field_is_scalar)
  CALL swap_bounds(thetav,                                               &
                   tdims_s%i_len - 2*tdims_s%halo_i,                     &
                   tdims_s%j_len - 2*tdims_s%halo_j,                     &
                   tdims_s%k_len,                                        &
                   tdims_s%halo_i, tdims_s%halo_j,                       &
                   fld_type_p,swap_field_is_scalar)
  CALL swap_bounds(wetrho_r_sq_n,                                        &
                   pdims_s%i_len - 2*pdims_s%halo_i,                     &
                   pdims_s%j_len - 2*pdims_s%halo_j,                     &
                   pdims_s%k_len,                                        &
                   pdims_s%halo_i, pdims_s%halo_j,                       &
                   fld_type_p,swap_field_is_scalar)
END IF ! (first_atmstep_call)


!======== LBC Update =========================================================
! ---------------------------------------------------------------------
!   Section 0.2  Update lbcs for LAMs
! ---------------------------------------------------------------------
IF ((ErrorStatus == 0) .AND. (model_type == mt_lam)) THEN
  IF (ltimer) CALL timer('AS LAM_LBCS',5)

  IF ( first_atmstep_call ) THEN
    lbc_size=lenrima(fld_type_u,halo_type_extended, rima_type_norm)
    lbc_size_new = lenrima(fld_type_p,halo_type_extended, rima_type_norm)

    ! Reset u_lbc array to same size at p lbc array
    lenrima(fld_type_u,halo_type_extended, rima_type_norm) =                &
                      lenrima(fld_type_p,halo_type_extended, rima_type_norm)
  END IF

  !--------------------------------------------------------------
  ! Idealised UM LBC forcing
  !  If active, update lateral boundary arrays to contain
  !  idealised namelist profile data interpolated in time.
  !--------------------------------------------------------------
  IF (L_idealised_data .AND. L_force_lbc) THEN
    errorstatus = 1
    cmessage    = 'l_force_bc option still under code development.'
    CALL ereport(routinename, errorstatus, cmessage)
    ! DEPENDS ON: eg_idl_force_lbc
    CALL eg_IDL_Force_LBC (                                           &
                      row_length, rows, offx, offy,                   &
                      halo_i, halo_j,                                 &
                      lenrima(1,1,rima_type_norm),                    &
                      timestep_number,                                &
                      u_lbc, v_lbc,                                   &
                      theta_lbc,q_lbc,                                &
                      u_adv_lbc,v_adv_lbc,                            &
                      exner_lbc,                                      &
                      r_theta_levels, r_rho_levels,                   &
                      eta_theta_levels, eta_rho_levels)
  END IF ! on (L_idealised_data and L_force_lbc)

  IF ( first_atmstep_call ) THEN
    !--------------------------------------------------------------
    !           Update primary fields with LAM LBC data
    !--------------------------------------------------------------
    ! DEPENDS ON: update_lam_lbcs
    CALL update_lam_lbcs(                                              &
         r_rho_levels, r_theta_levels,                                 &
         row_length,rows,n_rows,                                       &
         tr_vars,tr_lbc_vars,tr_levels,                                &
         a_max_trvars,A_tr_active_lbc_index,                           &
         tr_ukca,tr_lbc_ukca,                                          &
         offx,offy,halo_i,halo_j,at_extremity,                         &
         L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,                        &
         L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc, L_pc2_lbc, &
         L_murk, L_murk_lbc,                                           &
         L_LBC_balance, L_int_uvw_lbc,                                 &
         L_dust_div1, L_dust_div1_lbc, L_dust_div2, L_dust_div2_lbc,   &
         L_dust_div3, L_dust_div3_lbc, L_dust_div4, L_dust_div4_lbc,   &
         L_dust_div5, L_dust_div5_lbc, L_dust_div6, L_dust_div6_lbc,   &
         L_so2,L_so2_lbc,L_dms,L_dms_lbc,L_so4_aitken,L_so4_aitken_lbc,&
         L_so4_accu, L_so4_accu_lbc, L_so4_diss, L_so4_diss_lbc,       &
         L_nh3, L_nh3_lbc,                                             &
          L_soot_new_lbc,  L_soot_agd_lbc,       &
          L_soot_cld_lbc, L_bmass_new, L_bmass_new_lbc,     &
         L_bmass_agd, L_bmass_agd_lbc, L_bmass_cld, L_bmass_cld_lbc,   &
          L_ocff_new_lbc,                                   &
          L_ocff_agd_lbc,  L_ocff_cld_lbc,       &
         L_nitr_acc, L_nitr_acc_lbc, L_nitr_diss, L_nitr_diss_lbc,     &
         rimwidtha(rima_type_norm),rimweightsa,                        &
         lenrima(1,1,rima_type_norm),                                  &
         lbc_sizea(1,1,1,rima_type_norm),                              &
         lbc_starta(1,1,1,rima_type_norm),                             &
         theta_lbc,q_lbc,qcl_lbc,                                      &
         qcf_lbc,qcf2_lbc,qrain_lbc,                                   &
         qgraup_lbc, cf_bulk_lbc,cf_liquid_lbc,                        &
         cf_frozen_lbc, rho_lbc,exner_lbc,                             &
         u_lbc,v_lbc,w_lbc,                                            &
         u_adv_lbc,v_adv_lbc,w_adv_lbc,                                &
         murk_lbc,                                                     &
         dust_div1_lbc, dust_div2_lbc, dust_div3_lbc,                  &
         dust_div4_lbc, dust_div5_lbc, dust_div6_lbc,                  &
         so2_lbc, dms_lbc, so4_aitken_lbc,                             &
         so4_accu_lbc,so4_diss_lbc,nh3_lbc,                            &
         soot_new_lbc, soot_agd_lbc, soot_cld_lbc,                     &
         bmass_new_lbc, bmass_agd_lbc, bmass_cld_lbc,                  &
         ocff_new_lbc, ocff_agd_lbc, ocff_cld_lbc,                     &
         nitr_acc_lbc, nitr_diss_lbc,                                  &
         tracer_lbc,tracer_ukca_lbc,                                   &
         thetav,m_v,m_cl,m_cf,                                         &
         m_cf2, m_r, m_gr,                                             &
         cf_bulk,cf_liquid,cf_frozen,                                  &
         dryrho,exner,                                                 &
         u,v,w,                                                        &
         u_adv,v_adv,w_adv,                                            &
         murk,                                                         &
         dust_div1, dust_div2, dust_div3,                              &
         dust_div4, dust_div5, dust_div6,                              &
         so2, dms, so4_aitken,so4_accu,                                &
         so4_diss, nh3,                                                &
         soot_new, soot_agd, soot_cld,                                 &
         bmass_new, bmass_agd, bmass_cld,                              &
         ocff_new, ocff_agd, ocff_cld,                                 &
         nitr_acc, nitr_diss,                                          &
         delta_phi, delta_lambda,                                      &
         base_phi, base_lambda,                                        &
         datastart, lat_rot_NP,                                        &
         global_row_length, global_rows,                               &
         tracer, tracer_ukca )

    ! Fix surface exner
    CALL eg_calc_p_star(model_levels, row_length, rows, exner,    &
                     thetav, m_v, m_cl, m_cf, m_r, m_gr, m_cf2,   &
                     g_theta, exner_surf, psi_w_surf,             &
                     psi_w_lid)

    ! Must now re-calculate the pressure-based variables, namely pressure
    ! on both rho and theta levels, exner on theta levels and pstar so that
    ! they are all consistent with the new LBC-updated values of exner on
    ! rho levels.

    ! *** NOTE ***
    ! This is picked up by the later eg_thetav_theta call, but might need
    ! a call added here as well

  END IF ! first_atmstep_call

  ! time level n fields have been updated so update guess for np1 fields

  IF (ltimer) CALL timer('AS LAM_LBCS',6)
END IF     !   model_type  ==  mt_lam
! ---------------------------------------------------------------------
! End Section 0.2 Update lbcs for LAMs
! ---------------------------------------------------------------------
!======== LBC Update =========================================================

! =====================================================================
! needed because IAU in initial and or LBCs may have changed the fields.
! this was previously at the end of atm_step.

CALL swap_bounds(u,                                                    &
                 udims_s%i_len - 2*udims_s%halo_i,                     &
                 udims_s%j_len - 2*udims_s%halo_j,                     &
                 udims_s%k_len,                                        &
                 udims_s%halo_i, udims_s%halo_j,                       &
                 fld_type_u,swap_field_is_vector)
CALL swap_bounds(v,                                                    &
                 vdims_s%i_len - 2*vdims_s%halo_i,                     &
                 vdims_s%j_len - 2*vdims_s%halo_j,                     &
                 vdims_s%k_len,                                        &
                 vdims_s%halo_i, vdims_s%halo_j,                       &
                 fld_type_v,swap_field_is_vector)
CALL swap_bounds(w,                                                    &
                 wdims_s%i_len - 2*wdims_s%halo_i,                     &
                 wdims_s%j_len - 2*wdims_s%halo_j,                     &
                 wdims_s%k_len,                                        &
                 wdims_s%halo_i, wdims_s%halo_j,                       &
                 fld_type_p,swap_field_is_scalar)

IF (iau_in_initial .OR. (model_type == mt_lam) .OR.                    &
                        (model_type == mt_bi_cyclic_lam)) THEN
  ! etadot bounds updated inside routine.
  CALL init_etadot()

  IF (timestep_number == 1)                                            &
    CALL init_psiw (psi_w_surf, psi_w_lid, l_shallow)
ELSE
  CALL swap_bounds(etadot,                                             &
                   wdims_s%i_len - 2*wdims_s%halo_i,                   &
                   wdims_s%j_len - 2*wdims_s%halo_j,                   &
                   wdims_s%k_len,                                      &
                   wdims_s%halo_i, wdims_s%halo_j,                     &
                   fld_type_p,swap_field_is_scalar)
END IF

IF (iau_in_initial) iau_in_initial=.FALSE.

CALL eg_set_adv_winds(u,v,etadot,                                     &
                      u_adv,v_adv,w_adv,row_length,rows,n_rows,       &
                      model_levels, halo_i, halo_j, l_shallow)

! =====================================================================


! The arrays u/v/w/dryrho do not need to update the bounds at this point.

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE ( i, j, k )
!$OMP DO SCHEDULE(STATIC)
DO j = pdims_s%j_start, pdims_s%j_end
  DO i = pdims_s%i_start, pdims_s%i_end
    exner_surf_np1(i,j) = exner_surf(i,j)
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels
  DO j = udims_s%j_start, udims_s%j_end
    DO i = udims_s%i_start, udims_s%i_end
      u_np1(i,j,k) = u(i,j,k)
    END DO
  END DO
  DO j = vdims_s%j_start, vdims_s%j_end
    DO i = vdims_s%i_start, vdims_s%i_end
      v_np1(i,j,k) = v(i,j,k)
    END DO
  END DO
  DO j = pdims_s%j_start, pdims_s%j_end
    DO i = pdims_s%i_start, pdims_s%i_end
      dryrho_np1(i,j,k) = dryrho(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO k = 0, model_levels
  DO j = tdims_s%j_start, tdims_s%j_end
    DO i = tdims_s%i_start, tdims_s%i_end
      thetav_np1(i,j,k) = thetav(i,j,k)
    END DO
  END DO
  DO j = wdims_s%j_start, wdims_s%j_end
    DO i = wdims_s%i_start, wdims_s%i_end
      w_np1(i,j,k)      = w(i,j,k)
      etadot_np1(i,j,k) = etadot(i,j,k)
    END DO
  END DO
  DO j = tdims_s%j_start, tdims_s%j_end
    DO i = tdims_s%i_start, tdims_s%i_end
      m_v_np1(i,j,k)    = m_v(i,j,k)
      m_cl_np1(i,j,k)   = m_cl(i,j,k)
      m_cf_np1(i,j,k)   = m_cf(i,j,k)
    END DO
  END DO

  DO j = pdims_s%j_start, pdims_s%j_end
    DO i = pdims_s%i_start, pdims_s%i_end
      exner_np1(i,j,k+1) = exner(i,j,k+1)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

IF ( l_mcr_qrain ) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k = 0, model_levels
    DO j = tdims_s%j_start, tdims_s%j_end
      DO i = tdims_s%i_start, tdims_s%i_end
        m_r_np1(i,j,k)   = m_r(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF ( l_mcr_qgraup ) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k = 0, model_levels
    DO j = tdims_s%j_start, tdims_s%j_end
      DO i = tdims_s%i_start, tdims_s%i_end
        m_gr_np1(i,j,k)  = m_gr(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF ( l_mcr_qcf2 ) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k = 0, model_levels
    DO j = tdims_s%j_start, tdims_s%j_end
      DO i = tdims_s%i_start, tdims_s%i_end
        m_cf2_np1(i,j,k) = m_cf2(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO
END IF
!$OMP END PARALLEL

CALL eg_alpha_ramp()

IF (ltimer) CALL timer ('AS Solver',5)

CALL eg_sisl_resetcon(                                                  &
         row_length, rows, n_rows, model_levels, halo_i, halo_j,        &
         offx, offy,  datastart,l_slice, l_test_tracer,                 &
         l_shallow, l_const_grav,                                       &
         a_realhd(rh_z_top_theta),thetav, dryrho, exner_surf, exner,    &
         m_v, m_r, m_gr, m_cl, m_cf, m_cf2,                             &
         f1_comp, f2_comp, f3_comp, ih )

CALL eg_set_helm_lhs(row_length, rows, n_rows, model_levels, ih)

IF (ltimer) CALL timer ('AS Solver',6)

IF ( first_atmstep_call .OR. subseq_to_dump ) THEN
  total_rho_init=eg_total_mass(dryrho, l_exclude_rim=.TRUE.)
  total_aam_init=eg_total_aam(u, dryrho, l_shallow, l_exclude_rim=.TRUE.)
  total_ke_init=eg_total_ke(u, v, w, dryrho, l_exclude_rim=.TRUE.)
END IF
IF (model_type == mt_lam .AND. conserve_dry_mass /= not_conserved) THEN
  total_rho_n = eg_total_mass(dryrho, l_exclude_rim=.TRUE.)
END IF

IF ( subseq_to_dump ) THEN
  CALL swap_bounds(u_adv,                                                &
                   udims_l%i_len - 2*udims_l%halo_i,                     &
                   udims_l%j_len - 2*udims_l%halo_j,                     &
                   udims_l%k_len,                                        &
                   udims_l%halo_i, udims_l%halo_j,                       &
                   fld_type_u,swap_field_is_vector)
  CALL swap_bounds(v_adv,                                                &
                   vdims_l%i_len - 2*vdims_l%halo_i,                     &
                   vdims_l%j_len - 2*vdims_l%halo_j,                     &
                   vdims_l%k_len,                                        &
                   vdims_l%halo_i, vdims_l%halo_j,                       &
                   fld_type_v,swap_field_is_vector)
  CALL swap_bounds(w_adv,                                                &
                   wdims_l%i_len - 2*wdims_l%halo_i,                     &
                   wdims_l%j_len - 2*wdims_l%halo_j,                     &
                   wdims_l%k_len,                                        &
                   wdims_l%halo_i, wdims_l%halo_j,                       &
                   fld_type_p,swap_field_is_scalar)
  subseq_to_dump = .FALSE.
END IF

IF ( first_atmstep_call ) THEN
  IF ( timestep_number == 1 .AND. printstatus >= prstatus_normal) THEN
    IF (ltimer) CALL timer ('AS Diagnostics',5)
      CALL print_windmax_1(u, v, w, conserve_dry_mass)
      IF (nsteps_consv_print > 0) THEN
        CALL print_conservation_diag(total_rho_init, total_aam_init,         &
             total_ke_init, initial_timestep)
      END IF
    IF (ltimer) CALL timer ('AS Diagnostics',6)
  END IF

  dx_frame_cnt = 0
END IF

! Save initial data for plotting if requested
IF ( tstep_plot_start >= 0 .AND.                                            &
    timestep_number >= tstep_plot_start) THEN

  IF ( MOD(timestep_number,tstep_plot_frequency) == 0                     &
       .OR. dx_frame_cnt == 0 ) THEN

    CALL eg_dx_diags(dx_frame_cnt, mype,                                  &
                     row_length, rows, n_rows, model_levels,              &
                     offx, offy, halo_i, halo_j,                          &
                     u, v, w, dryrho, thetav, exner, exner_surf,          &
                     l_dry, m_v, m_cl, m_cf, m_r, m_gr, m_cf2)
  END IF
END IF

!-----------------------------------------------------------------------
! Recompute "starred" Coriolis terms
!-----------------------------------------------------------------------

IF (ltimer) CALL timer ('AS Solver',5)

CALL eg_coriolis_star(dryrho, m_v, m_cl, m_cf, m_r, m_gr, m_cf2)

IF (ltimer) CALL timer ('AS Solver',6)

! ----------------------------------------------------------------------
! Section 1.0  Initialise semi-lagrangian
! ----------------------------------------------------------------------

! Set timelevel n dependent quantities

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE ( i, j, k )
!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels
  DO j = udims_s%j_start, udims_s%j_end
    DO i = udims_s%i_start, udims_s%i_end
      r_u(i,j,k) = 0.0
    END DO
  END DO
  DO j = vdims_s%j_start, vdims_s%j_end
    DO i = vdims_s%i_start, vdims_s%i_end
      r_v(i,j,k) = 0.0
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO k = 0, model_levels
  DO j = wdims%j_start, wdims%j_end
    DO i = wdims%i_start, wdims%i_end
      r_w(i,j,k)     = 0.0
    END DO
  END DO
  DO j = tdims_s%j_start, tdims_s%j_end
    DO i = tdims_s%i_start, tdims_s%i_end
      r_theta(i,j,k) = 0.0
    END DO
  END DO
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      r_m_v(i,j,k)   = 0.0
      r_m_cl(i,j,k)  = 0.0
      r_m_cf(i,j,k)  = 0.0
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

IF ( l_mcr_qrain  ) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k = 0, model_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        r_m_r(i,j,k)   = 0.0
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF ( l_mcr_qgraup ) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k = 0, model_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        r_m_gr(i,j,k)  = 0.0
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF ( l_mcr_qcf2 ) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k = 0, model_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        r_m_cf2(i,j,k) = 0.0
      END DO
    END DO
  END DO
!$OMP END DO
END IF
!$OMP END PARALLEL

l_call_from_solver = .FALSE.

! ***********************************************************************
!  CONVERSION FROM EG to ND prognostics! (start)
IF (ltimer) CALL timer('AS CONVERT',5)

! Compute theta from thetav _AND_ exner_surf _AND_ p_theta_levels
! for use in physics routines
CALL eg_thetav_theta                                                    &
                  (thetav, theta, m_v,                                  &
                   p, pstar, p_theta_levels,                            &
                   exner, exner_surf, exner_theta_levels)

IF (ltimer) CALL timer('AS CONVERT',6)

IF ((.NOT. l_dry) .AND. (.NOT. l_mr_physics)) THEN

  IF (ltimer) CALL timer('AS CONVERT',5)
  CALL eg_mix_to_q                                                      &
                  (tdims_l,tdims_s,                                     &
                   m_v, m_cl, m_cf,                                     &
                   m_cf2, m_r, m_gr,                                    &
                   l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,               &
                   q, qcl, qcf,                                         &
                   qcf2, qrain, qgraup)
  IF (ltimer) CALL timer('AS CONVERT',6)

  !-----------------------------------------------------------------------
  ! Allocate and initialise arrays for total timestep increments
  ! NOTE it is very important this is always done before any process
  !      changes the fields at the beginning of the timestep.
  ! This needs to work on q (specific humidities) as well as mixing
  ! ratios thus requiring ugly use of IF test.
  !-----------------------------------------------------------------------
  IF (ltimer .OR. lstashdumptimer) CALL timer('AS STASH',5)

  IF (sf(0,30) .OR. flag_upd_helicity_5k) CALL eot_incs_init()

  IF (ltimer .OR. lstashdumptimer) CALL timer('AS STASH',6)

ELSE

  ! Ensure using q before overwritten by MRs.
  IF (ltimer .OR. lstashdumptimer) CALL timer('AS STASH',5)

  IF (sf(0,30) .OR. flag_upd_helicity_5k) CALL eot_incs_init()

  IF (ltimer .OR. lstashdumptimer) CALL timer('AS STASH',6)
  IF (ltimer) CALL timer('AS CONVERT',5)
!$OMP PARALLEL DEFAULT(NONE) PRIVATE ( i, j, k )                        &
!$OMP SHARED(model_levels,offy,offx,rows,row_length,q,m_v, qcl, m_cl,   &
!$OMP        qcf,m_cf, l_mcr_qcf2,qcf2,m_cf2,l_mcr_qrain,m_r,qrain,     &
!$OMP        qgraup, m_gr, l_mcr_qgraup)
!$OMP DO SCHEDULE(STATIC)
  DO k = 0, model_levels
    DO j = 1-offy,rows+offy              
      DO i = 1-offx, row_length+offx        
        q(i,j,k)    = m_v(i,j,k)
        qcl(i,j,k)   = m_cl(i,j,k)
        qcf(i,j,k)   = m_cf(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

  IF (l_mcr_qcf2  ) THEN 
!$OMP DO SCHEDULE(STATIC)
    DO k = 0, model_levels
      DO j = 1-offy,rows+offy              
        DO i = 1-offx, row_length+offx        
          qcf2(i,j,k)   = m_cf2(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

  IF (l_mcr_qrain  )   THEN
!$OMP DO SCHEDULE(STATIC)
    DO k = 0, model_levels
      DO j = 1-offy,rows+offy              
        DO i = 1-offx, row_length+offx        
          qrain(i,j,k)   = m_r(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

  IF (l_mcr_qgraup  )     THEN
!$OMP DO SCHEDULE(STATIC)
    DO k = 0, model_levels
      DO j = 1-offy,rows+offy              
        DO i = 1-offx, row_length+offx        
          qgraup(i,j,k)   = m_gr(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO
  END IF

!$OMP END PARALLEL

  IF (ltimer) CALL timer('AS CONVERT',6)

END IF

! ---------------------------------------------------------------------
!  CONVERSION FROM EG to ND prognostics! (end)
!
! ***********************************************************************

! ---------------------------------------------------------------------
!
! compute the wet_to_dry conversion factor at time level n. Needed for
! energy correction and some diagnostics
CALL wet_to_dry_n_calc(q,qcl,qcf,qcf2,qrain,qgraup)

! ---------------------------------------------------------------------

IF (ltimer) CALL timer('AS Stochastic_Phys',5)

!  --------- UM Section 35---- Stochastic Physics Setup -------------
! If SP coefficients were not present in the input dump, allocate
! space in the field of constants to store the coefficients.
IF ( (first_atmstep_call) .AND. (l_rp2 .OR. l_skeb2 .OR. l_spt) ) THEN
  
  ! Check if stochastics physics values were already present in the dump. 
  ! If stochastic physics data is present, we assume all the data in additional
  ! parameters is related to stochastic physics and calculate the positions and
  ! sizes of the relevant component arrays, starting from the beginning of the
  ! additional parameters array.
  ! If the integer header shows that stochastic physics data is not present in
  ! the dump, we need to allocated space for it. If there is already data in 
  ! the additional parameters array, we retain this data and allocated space
  ! after the existing data.
  IF ( a_inthd(ih_stochastic_flag) > stph_seed_present ) THEN
    CALL umPrint('SP data present, assigning to array',src='atm_step_4A')
    CALL assign_a_flddepc_input_values()
  ELSE
    ! if not, allocate space in the field-dependant constants in the header
    CALL umPrint('SP coeffcients not present in dump',src='atm_step_4A')
    CALL umPrint('calling allocate_sp_coefficients to allocate space', &
                 src='atm_step_4A')
    CALL allocate_sp_coefficients()
  END IF

END IF

IF (first_atmstep_call) THEN
  IF (l_rp2 .OR. l_skeb2 .OR. l_spt .OR. i_pert_theta /= off) THEN
    
    !  Called only at the first time-step with no arguments
    CALL stph_setup ( )
      
    ! Read the stochastic physics data from the dump header arrays 
    ! if applicable
    IF (timestep_number > 1 .AND.                                     &
        a_inthd(ih_stochastic_flag) >= stph_seed_present ) THEN
      CALL stph_seed_copy_from_dump ( )
    END IF
  END IF
END IF
!  --------- UM Section 35---- Stochastic Physics Setup END ---------

!  ----- UM Section 35---- Stochastic Physics Random Parameters -----
! Call to the RANDOM PARAMETERS2 (STPH_RP2) subroutine
IF (l_physics .AND. l_rp2) THEN
  CALL stph_rp2(model_levels,                                         &
                rhcrit,rhcrit_max,rhcrit_min,                         &
                m_ci,m_ci_max,m_ci_min,                               &
                charnock)

END IF

! Scale lai by its random parameter if using the RP2b scheme.
!
! The result of this scaling is stored in the local array
! lai_pft_in, which is passed to the subroutines instead of
! lai_pft.  This is done to ensure the scaling is always applied
! to the most recently input ancillary data, held in lai_pft
! (for some model configurations, the LAI ancillary may be updated
! during the forecast).
!
! Note that the scaling should only be used when lai_pft is NOT being updated
! in the subroutines (i.e. not acting as a prognostic).  Since the
! JULES routines TRIFFID and PHENOL both treat lai as a prognostic
! they cannot be used with the RP2b scheme.
!
! If the RP2b scheme is not being used, a straight copy of lai_pft
! is made (lai_pft_in) to be passed to the relevant subroutines.  This
! is then passed back to lai_pft at the end of the subroutines.
! This ensures that for the caes where lai is a prognostic, any changes
! made in the subroutines are preserved between timesteps.

IF ( l_physics .AND. l_rp2 .AND. i_rp_scheme == i_rp2b ) THEN
  ! Scale lai by random parameter
  DO i = 1, npft
    lai_pft_in(:,i) = lai_mult_rp(i)*lai_pft(:,i)
  END DO
ELSE
  lai_pft_in(:,:) = lai_pft(:,:)
END IF

!  ----- UM Section 35---- Stochastic Physics Random Parameters END ---
IF (ltimer) CALL timer('AS Stochastic_Phys',6)

! ---------------------------------------------------------------
!    diagnostic printing at beginning of timestep 1
!    (which is really timestep 0 values so lets record it as such)
! ---------------------------------------------------------------
IF ( L_diag_print .AND. timestep_number ==1 ) THEN
  IF (ltimer) CALL timer ('AS Diagnostics',5)
  CALL Print_diag_4A(u, v, theta, Wetrho_r_sq_n, w, q, qcl, qcf,    &
                 rows, n_rows, row_length, model_levels,         &
                 offx, offy,timestep_number-1,                   &
                 rpemax, rpemin, ipesum, rpesum,                 &
                 max_w_run, max_wind_run, min_theta1_run,        &
                 dtheta1_run, max_div_run, min_div_run,          &
                 min_lapse_run, max_shear_run, time_max_shear,   &
                 time_div_max, time_div_min, time_lapse_min,     &
                 time_w_max, time_max_wind, time_theta1_min,     &
                 max_KE_run, min_KE_run, max_noise_run,          &
                 time_KE_max, time_KE_min, time_noise_max )
  IF ( L_flush6 ) CALL umPrintFlush()
  IF (ltimer) CALL timer ('AS Diagnostics',6)
END IF     !  timestep_number ==1

!=== Polar filter + diffusion section ==================================
IF ( .NOT. first_atmstep_call ) THEN
  ! ----------------------------------------------------------------------
  ! Section 0.4  Filter winds and theta near poles if active
  !              Do horizontal diffusion as a filter if active
  ! ----------------------------------------------------------------------

  IF ( L_filter ) THEN
    IF (ltimer) CALL timer('AS Filter',5)

    CALL filter_diag_printing1()

    CALL eg_NI_filter_Ctl(  thetav, u, v, w, etadot, exner,          &
                         exner_theta_levels,                         &
                         row_length, rows, n_rows, model_levels,     &
                         r_theta_levels, r_rho_levels,               &
                         r_at_u, r_at_v,                             &
                         max_121_rows, u_sweeps, v_sweeps,           &
                         global_u_filter, global_v_filter,           &
                         u_begin, u_end, v_begin, v_end,             &
                         diff_coeff_phi, diff_coeff_u, diff_coeff_v, &
                         first_constant_r_rho_level,                 &
                         first_constant_r_rho_level_m1,              &
                         top_filt_start, top_filt_end,               &
                         up_diff, max_updiff_levels,                 &
                         horizontal_level,                           &
                         offx, offy, halo_i, halo_j,                 &
                         nproc_y, at_extremity,                      &
                         L_pftheta, L_pfuv,                          &
                         L_pfw, L_pfexner, L_diff_exner,             &
                         L_diff_thermo, L_diff_wind, L_diff_w,       &
                         L_pofil_hadgem2,                            &
                         xi1_u, xi1_p, xi2_p, xi2_v,                 &
                         pole_consts, gc_proc_row_group,             &
                         nproc, gc_proc_col_group,                   &
                         global_row_length,                          &
                         csxi2_v, csxi2_p, delta_lambda, delta_phi,  &
                         STASHwork13)

    CALL filter_diag_printing2()

    IF (L_pfexner .AND. L_pofil_new) THEN
      ! DEPENDS ON: consistent_pressure
      CALL     Consistent_Pressure ( exner,                           &
               offx,offy,halo_i,halo_J,                               &
               row_length,rows,model_levels,                          &
               r_theta_levels, r_rho_levels, dryrho,                  &
               p, pstar, p_theta_levels,exner_theta_levels)
    END IF  !  (L_pfexner and l_pofil_new)

    IF (ltimer) CALL timer('AS Filter',6)
  END IF    !  L_filter

END IF ! timestep_number > 1 -> filter
!=== End Polar filter + diffusion section ===================================

! ----------------------------------------------------------------------
! Section 1.0  Call Atmospheric Physics1
! ----------------------------------------------------------------------

IF (ltimer) CALL timer('AS Atmos_Phys1 (AP1)',5)

! Biogenic aerosol climatology, on theta points, but ignores level 0:
IF (l_use_biogenic) THEN
!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                      &
!$OMP& SHARED( tdims, biogenic, arclbiog_bg )                          &
!$OMP& PRIVATE( i, j, k )
  DO k=1, tdims%k_end
    DO j=tdims%j_start, tdims%j_end
      DO i=tdims%i_start, tdims%i_end
        biogenic(i,j,k) = arclbiog_bg(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! Aerosol climatologies - Model switches and climatologies are
! gathered into bigger arrays.

! First set the internal model switches according to the
! CNTLATM settings, and determine how many components we need.
CALL set_arcl_dimensions(l_use_arclbiom, l_use_arclblck,                &
                         l_use_arclsslt, l_use_arclsulp,                &
                         l_use_arcldust, l_use_arclocff,                &
                         l_use_arcldlta, n_arcl_species,                &
                         n_arcl_compnts, l_use_arcl )

! If the aerosol climatology for NWP is used, n_arcl_species
! is larger than 0. In that case, allocate the array gathering
! component mass-mixing ratio and take the values from the
! arrays in atm_fields_mod.

IF (n_arcl_species > 0) THEN
  ! the arcl array is on theta grid, but ignores level 0:
  ALLOCATE(arcl(tdims%i_start:tdims%i_end,                              &
                tdims%j_start:tdims%j_end,                              &
                            1:tdims%k_end,                              &
                            1:n_arcl_compnts) )
  CALL set_arcl_clim(n_arcl_compnts,                                    &
                     ! Internal model switches
                     l_use_arcl,                                        &
                     ! Climatologies from ancillary files
                     arclbiom_fr, arclbiom_ag, arclbiom_ic,             &
                     arclblck_fr, arclblck_ag,                          &
                     arclsslt_fi, arclsslt_jt,                          &
                     arclsulp_ac, arclsulp_ak, arclsulp_di,             &
                     arcldust_b1, arcldust_b2, arcldust_b3,             &
                     arcldust_b4, arcldust_b5, arcldust_b6,             &
                     arclocff_fr, arclocff_ag, arclocff_ic,             &
                     arcldlta_dl,                                       &
                     ! Internal climatology array
                     arcl,                                              &
                     ! Component array indices
                     i_arcl_compnts )

ELSE
  ALLOCATE ( arcl(1,1,1,1) )
END IF

IF (errorstatus == 0) THEN
  IF (l_glomap_mode_clim) THEN
    ! Allocate space for STASH
    IF (.NOT. ALLOCATED(stashwork54))                                          &
      ALLOCATE (stashwork54(stash_maxlen(stashcode_glomap_clim_sec,atmos_im)))
    DO l=1,stash_maxlen(stashcode_glomap_clim_sec,atmos_im)
      stashwork54(l)=0.0
    END DO
  END IF
END IF

!
! UKCA_RADAER: Obtain current UKCA setup and allocate arrays.
!              Most of the work needs only be made for the first
!              Atm_Step_4A() call.
!
IF (errorstatus == 0) THEN

  IF (l_glomap_clim_radaer .OR. l_ukca_radaer) THEN
    
    IF (first_atmstep_call) THEN
      IF (l_glomap_clim_radaer) THEN
        ! CALL mode setup routine to set modes, components...
        SELECT CASE(i_glomap_clim_setup)
        CASE (i_gc_sussocbc_5mode)
          CALL ukca_mode_sussbcoc_5mode
        CASE DEFAULT
          errorstatus = 54000
          WRITE (cmessage,'(A,I0)') 'Unsupported value of ' //                 &
                                    'i_glomap_clim_setup. Value provided : ',  &
                                     i_glomap_clim_setup
          CALL ereport(RoutineName,errorstatus,Cmessage)
        END SELECT
        
        CALL ukca_radaer_init ( errorstatus,                                   &
                                Cmessage,                                      &
                                l_glomap_clim_radaer_sustrat,                  &
                                ukca_radaer )
      ELSE
        CALL ukca_radaer_init ( errorstatus,                                   &
                                Cmessage,                                      &
                                l_ukca_radaer_sustrat,                         &
                                ukca_radaer )
      END IF
      
      IF (errorstatus /= 0) THEN
        CALL Ereport(RoutineName, errorstatus, Cmessage)
      END IF
      
    END IF
    
    !
    ! Allocate those arrays that will receive UKCA output.
    !
    ALLOCATE(ukca_radaer%mix_ratio(row_length, rows,                    &
                       model_levels, ukca_radaer%n_cpnt))
    ALLOCATE(ukca_radaer%comp_vol(row_length, rows,                     &
                       model_levels, ukca_radaer%n_cpnt))

    ALLOCATE(ukca_radaer%dry_diam(row_length, rows,                     &
                             model_levels, ukca_radaer%n_mode))
    ALLOCATE(ukca_radaer%wet_diam(row_length, rows,                     &
                       model_levels, ukca_radaer%n_mode))
    ALLOCATE(ukca_radaer%modal_rho(row_length, rows,                    &
                       model_levels, ukca_radaer%n_mode))
    ALLOCATE(ukca_radaer%modal_wtv(row_length, rows,                    &
                       model_levels, ukca_radaer%n_mode))
    ALLOCATE(ukca_radaer%modal_vol(row_length, rows,                    &
                       model_levels, ukca_radaer%n_mode))
    ALLOCATE(ukca_radaer%modal_nbr(row_length, rows,                    &
                       model_levels, ukca_radaer%n_mode))

    !
    ! Populate ukca_radaer fields
    IF (l_glomap_clim_radaer) THEN
      CALL glomap_clim_radaer_get(errorstatus,Cmessage,ukca_radaer,stashwork54)
    ELSE
      CALL ukca_radaer_get(errorstatus,Cmessage,first_atmstep_call,ukca_radaer)
    END IF
    
    IF (errorstatus /= 0) THEN

      DEALLOCATE(ukca_radaer%mix_ratio)
      DEALLOCATE(ukca_radaer%comp_vol)
      DEALLOCATE(ukca_radaer%dry_diam)
      DEALLOCATE(ukca_radaer%wet_diam)
      DEALLOCATE(ukca_radaer%modal_rho)
      DEALLOCATE(ukca_radaer%modal_wtv)
      DEALLOCATE(ukca_radaer%modal_vol)
      DEALLOCATE(ukca_radaer%modal_nbr)

      CALL Ereport(RoutineName, errorstatus, Cmessage)

    END IF

  ELSE

    ! l_ukca_radaer is false. Allocate minimum space.
    ALLOCATE(ukca_radaer%mix_ratio(1, 1, 1, 1))
    ALLOCATE(ukca_radaer%comp_vol(1, 1, 1, 1))
    ALLOCATE(ukca_radaer%dry_diam(1, 1, 1, 1))
    ALLOCATE(ukca_radaer%wet_diam(1, 1, 1, 1))
    ALLOCATE(ukca_radaer%modal_rho(1, 1, 1, 1))
    ALLOCATE(ukca_radaer%modal_wtv(1, 1, 1, 1))
    ALLOCATE(ukca_radaer%modal_vol(1, 1, 1, 1))
    ALLOCATE(ukca_radaer%modal_nbr(1, 1, 1, 1))

  END IF ! l_ukca_radaer or l_glomap_clim_radaer

END IF ! errorstatus

! UKCA_CDNC: Allocate array and obtain Cloud Droplet No. Concentration from D1
IF (errorstatus == 0) THEN
  IF ( (l_ukca             .AND. ( l_ukca_aie1 .OR. l_ukca_aie2 ) ) .OR.       &
       (l_glomap_mode_clim .AND. ( l_glomap_clim_arg_act .OR.                  &
                                   l_glomap_clim_aie1    .OR.                  &
                                   l_glomap_clim_aie2 ) ) ) THEN
    
    ! ALLOCATE ukca_cdnc%cdnc & %cdnc3 in this CALL to keep atm_step_4A tidy
    IF ( l_glomap_mode_clim ) THEN
      CALL allocate_ukca_cdnc ( ukca_cdnc,                                     &
                                cdnc_dim1,                                     &
                                cdnc_dim2,                                     &
                                cdnc_dim3,                                     &
                                stashwork54 )
    ELSE
      CALL allocate_ukca_cdnc ( ukca_cdnc,                                     &
                                cdnc_dim1,                                     &
                                cdnc_dim2,                                     &
                                cdnc_dim3 )
    END IF
    
  ELSE
    ! Allocate dummy ukca_cdnc%cdnc & %cdnc3 fields to pass to atmos_physics1
    IF (.NOT. ALLOCATED(ukca_cdnc%cdnc))  ALLOCATE(ukca_cdnc%cdnc(1,1,1))
    IF (.NOT. ALLOCATED(ukca_cdnc%cdnc3)) ALLOCATE(ukca_cdnc%cdnc3(1,1,1))
    ukca_cdnc%cdnc(:,:,:)  = 0.0
    ukca_cdnc%cdnc3(:,:,:) = 0.0
    cdnc_dim1 = 1
    cdnc_dim2 = 1
    cdnc_dim3 = 1
  END IF
END IF

IF (errorstatus == 0) THEN
  IF (l_glomap_mode_clim) THEN
    
    ! DEPENDS ON: stash
    CALL stash( atmos_sm, atmos_im, stashcode_glomap_clim_sec, stashwork54,    &
                errorstatus, cmessage)
    
    IF (errorstatus /=  0) THEN
      cmessage=" Section 54: Error in STASH: " // cmessage
      errorstatus = ABS( errorstatus )
      CALL ereport(RoutineName,errorstatus,cmessage)
    END IF
    
    IF (ALLOCATED(stashwork54)) DEALLOCATE(stashwork54)
    
  END IF
END IF

! Search for plume scavenging diagnostic requests and set logical
! if any are found
IF (l_ukca .AND. l_ukca_mode .AND. l_ukca_plume_scav) THEN
  item_plume_diags_last = item1_plume_diags + nmax_plume_diags -1
  IF (ANY(sf(item1_plume_diags:item_plume_diags_last,stashcode_glomap_sec)))  &
    l_ukca_plume_diags = .TRUE.
END IF

! EasyAerosol: Allocate arrays and obtain distributions.

IF (ErrorStatus == 0) THEN
 
  IF (l_easyaerosol_sw   .OR. l_easyaerosol_lw        .OR.  &
      l_easyaerosol_cdnc .OR. l_easyaerosol_autoconv) THEN

    CALL easyaerosol_read_input_ctl(row_length, rows, model_levels, &
                                    global_row_length, global_rows, &
                                    sw_spectrum(1)%basic%n_band, & 
                                    lw_spectrum(1)%basic%n_band, & 
                                    timestep_number, timestep, &
                                    first_atmstep_call, &
                                    easyaerosol_sw, &
                                    easyaerosol_lw, &
                                    easyaerosol_cdnc)

  ELSE
    
    ! EasyAerosol not used: allocate minimum sizes for bound checking.
    CALL allocate_easyaerosol_rad(easyaerosol_sw, 1, 1, 1, 1)
    CALL allocate_easyaerosol_rad(easyaerosol_lw, 1, 1, 1, 1)
    CALL allocate_easyaerosol_cdnc(easyaerosol_cdnc, 1, 1, 1)
    
  END IF

END IF ! ErrorStatus = 0

IF (l_physics .AND. errorstatus == 0) THEN

  ! DEPENDS ON: Atm_Step_phys_init
  CALL atm_step_phys_init(                                              &
     r_v, r_u, theta_star, q_star, qcl_star, qcf_star,                  &
     cf_star, cfl_star, cff_star, exner_lbc_real_tend,                  &
     w_lbc_real_tend, errorstatus, 'tropinit')

  ! DEPENDS ON: o3_to_3d
  CALL o3_to_3d(lexpand_ozone, i_ozone_int,                             &
                rows, row_length, model_levels, ozone_levels,           &
                halo_i, halo_j, offx, offy, at_extremity,               &
                a_realhd(rh_z_top_theta),                               &
                theta(tdims_s%i_start,tdims_s%j_start,1),               &
                exner_theta_levels(tdims_s%i_start,tdims_s%j_start,1),  &
                wetrho_r_sq_n,                                          &
                exner_rho_levels,                                       &
                nd_o3, o3,                                              &
                min_trop_level, max_trop_level,                         &
                l_o3_trop_level,l_o3_trop_height,                       &
                l_t_trop_level,l_t_trop_height,                         &
                o3_trop_level,o3_trop_height,                           &
                t_trop_level,t_trop_height,                             &
                gc_proc_row_group,                                      &
                global_row_length,                                      &
                ozone3d(o3dims2%i_start,o3dims2%j_start,1),             &
                errorstatus, cmessage)

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                      &
!$OMP& SHARED( rows, row_length, ozone3d )                             &
!$OMP& PRIVATE( i, j )
  DO j = 1, rows
    DO i = 1, row_length
      ozone3d(i,j,0)=ozone3d(i,j,1)
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: Atm_Step_phys_init
  CALL atm_step_phys_init(                                              &
     r_v, r_u, theta_star, q_star, qcl_star, qcf_star,                  &
     cf_star, cfl_star, cff_star, exner_lbc_real_tend,                  &
     w_lbc_real_tend, errorstatus, 'ozoninit')

END IF !  L_Physics

! DEPENDS ON: Atm_Step_phys_init
CALL atm_step_phys_init(                                                &
     r_v, r_u, theta_star, q_star, qcl_star, qcf_star,                  &
     cf_star, cfl_star, cff_star, exner_lbc_real_tend,                  &
     w_lbc_real_tend, errorstatus, 'microphy')

! dyn-tracers: Initialize variables
IF (l_pv_tracer .OR. l_diab_tracer) CALL dyn_tr_sav()

IF (l_idealised_data) THEN
  ! Allocate arrays used to hold idealised diagnostics for stash
  ! Note these may not be needed until later but there are two parts to the
  ! idealised forcing (1) before the advection and (2) after the advection
  ! therefore need to allocate now in case required by part (1).
  CALL alloc_ideal_diag(1)
END IF

IF (l_physics .AND. errorstatus == 0) THEN

  ! DEPENDS ON: Atm_Step_diag
  CALL atm_step_diag(1)

  ! only some code for pc2 in atm_step_alloc_4A is needed
  ! Convert to mixing ratios from specific humidities if needed
  ! DEPENDS ON: Atm_Step_alloc_4A
  CALL Atm_Step_alloc_4A(                                               &
            cf_star, cfl_star, cff_star, q_star, qcl_star, qcf_star,    &
            frac_control, r_u, r_v, r_w, errorstatus,'MixRatio')

  ! NB if you are changing the argument list to atmos_physics1, please
  ! do an equivalent change in routine scm_main to keep the single column
  ! model consistent.

  ! NOTE: R_u and  R_v have changed but scm_main has been left untouched

!$OMP  PARALLEL DEFAULT(NONE)                                          &
!$OMP& SHARED( tdims_s, tdims, theta_star, q_star, qcl_star,           &
!$OMP&         qcf_star, qcf2_star, qrain_star, qgraup_star,           &
!$OMP&         l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup )                 &
!$OMP& PRIVATE( i, j, k )
!$OMP DO SCHEDULE(STATIC)
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        theta_star (i,j,k) = 0.0
        q_star  (i,j,k) = 0.0
        qcl_star(i,j,k) = 0.0
        qcf_star(i,j,k) = 0.0
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

  IF (l_mcr_qcf2  )  THEN
!$OMP DO SCHEDULE(STATIC)
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qcf2_star  (i,j,k) = 0.0
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

  IF (l_mcr_qrain )  THEN
!$OMP DO SCHEDULE(STATIC)
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qrain_star (i,j,k) = 0.0
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

  IF (l_mcr_qgraup)  THEN
!$OMP DO SCHEDULE(STATIC)
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qgraup_star(i,j,k) = 0.0
        END DO
      END DO
    END DO
!$OMP END DO
  END IF
!$OMP END PARALLEL

  ! end debugging only

  IF (L_cosp) THEN
    !     This is done until a prognostic is developed
    IF (timestep_number == 1) THEN
      NULLIFY(cosp_crain_3d,cosp_csnow_3d)
    END IF
    IF (.NOT. ASSOCIATED(cosp_crain_3d)) THEN
      IF (PrintStatus >= PrStatus_Oper) THEN
        WRITE(umMessage,'(A,I5)')                                            &
          'COSP in ATM_STEP_4A: allocating cosp_crain_3d in tstep ',         &
          timestep_number
        CALL umPrint(umMessage,src='atm_step_4A')
      END IF
      ALLOCATE(cosp_crain_3d(row_length,rows,model_levels))
      cosp_crain_3d = 0.0
    END IF
    IF (.NOT. ASSOCIATED(cosp_csnow_3d)) THEN
      IF (PrintStatus >= PrStatus_Oper) THEN
        WRITE(umMessage,'(A,I5)')                                            &
          'COSP in ATM_STEP_4A: allocating cosp_csnow_3d in tstep ',         &
          timestep_number
        CALL umPrint(umMessage,src='atm_step_4A')
      END IF
      ALLOCATE(cosp_csnow_3d(row_length,rows,model_levels))
      cosp_csnow_3d = 0.0
    END IF
  ELSE
    IF (.NOT. ASSOCIATED(cosp_crain_3d)) THEN
      ALLOCATE(cosp_crain_3d(1,1,1))
      cosp_crain_3d = 0.0
    END IF
    IF (.NOT. ASSOCIATED(cosp_csnow_3d)) THEN
      ALLOCATE(cosp_csnow_3d(1,1,1))
      cosp_csnow_3d = 0.0
    END IF
  END IF

! Apply idealised forcing
  IF (l_idealised_data) THEN
    CALL eg_idl_forcing(                                                   &
! IN Data Fields.
       u, v, theta, exner_theta_levels, exner, p_theta_levels,             &
! IN/OUT
       theta_star,q, q_star, qcl_star, qcf_star, qcf2_star, qrain_star,    &
       qgraup_star , dryrho,  r_u, r_v, wetrho_r_sq_n,                     &
       t_surface, problem_number,                                          &
! error information
       errorstatus  )
  END IF

  ! DEPENDS ON: atmos_physics1
  CALL Atmos_Physics1(                                             &
  ! Parallel variables
         global_row_length, global_rows, nproc, nproc_x ,nproc_y           &
        ,g_rows, g_row_length                                              &

  ! model dimensions
        , row_length, rows, n_rows, land_field                             &
        , bl_levels, st_levels, sm_levels, Ozone_levels, cloud_levels      &
        , land_ice_points, soil_points, n_cca_lev, ntiles                  &
        , salt_dim1, salt_dim2, salt_dim3, tr_levels, tr_ukca              &
        , cdnc_dim1, cdnc_dim2, cdnc_dim3                                  &
        , co2_dim_len, co2_dim_row, co2_dim_lev                            &
        , n_arcl_species, n_arcl_compnts, i_arcl_compnts                   &

  ! model switches
        ,     l_lbc_old                                                    &
        ,     l_ukca_chem, l_ukca_set_trace_gases                          &
        ,     l_ukca_strat, l_ukca_strattrop                               &
        ,     l_ukca_prescribech4, l_use_arcl                              &

  ! model Parameters
        ,     RHcrit                                                       &
        ,min_trop_level, max_trop_level                                    &

  ! Pass position of greenhouse gases in tracer_ukca array, for chemical coupling
        ,     ngrgas ,grgas_addr                                           &
  ! parameter for stochastic physics random parameters2
        ,     m_ci                                                         &
  ! Vertical coordinate levels.
        , delta_lambda, delta_phi, lat_rot_NP, long_rot_NP                 &
  ! Time stepping information
        , i_year, i_day_number, i_hour, i_minute                           &
        , i_second, previous_time                                          &
  ! diagnostic info
        , STASHwork1,STASHwork2,STASHwork4,STASHwork6,STASHwork14          &
        , STASHwork21                                                      &
  ! Additional variables for SCM diagnostics
        ,     nSCMDpkgs, L_SCMDiags                                        &
  ! Data Fields.
        , theta, q, qcl, qcf, qcf2, qrain, qgraup, wetrho_r_sq_n, u, v, w  &
        , p, pstar, exner, exner_theta_levels, land                        &
        , p_theta_levels, frac_land,frac_control, ukca_cdnc%cdnc           &

  ! ancillary fields and fields needed to be kept from timestep to
  ! timestep
        ,land_index, rgrain_tile ,snsoot, canht_pft, ntml, cumulus         &
        ,ice_fraction,p_ice_fract_rad, p_ice_thick_rad                     &
        ,cca_dp, cca_md, cca_sh                                            &
        ,cca, ccb, cct, cclwp, ccw_rad, lcbase, totalppn                   &
        ,tstar,tstar_land,tstar_sea,p_tstar_sice                           &
        ,sice_alb, land_alb, snodep,p_snodep_sice                          &
        ,ozone3D(o3dims2%i_start,o3dims2%j_start,1)                        &
        , sw_incs, lw_incs, dirpar                                         &
        ,O3_trop_level, O3_trop_height, T_trop_level, T_trop_height        &
        ,zh, orog_sd, orog_grad_xx, orog_grad_xy                           &
        ,orog_grad_yy, cf_area, cf_bulk,cf_liquid,cf_frozen,murk_source    &
        ,arcl, soil_alb, obs_alb_sw, obs_alb_vis, obs_alb_nir, lai_pft_in  &
        ,snodep_tile, frac_typ,tstar_tile,z0_tile,dolr_field,lw_down       &
        ,sw_tile_rts,es_space_interp, rad_mask, cos_zenith_angle           &
        ,easyaerosol_sw, easyaerosol_lw, easyaerosol_cdnc                  &
  ! Variables for COSP
        ,cosp_crain_3d,cosp_csnow_3d                                       &
  ! IN JULES prognostics
        ,snowdepth, lake_h_ice, z0, chloro_sea                             &
  ! Variable to store BL w-variance diagnostic
        ,bl_w_var                                                          &
  ! IN/OUT
        ,theta_star, q_star, qcl_star, qcf_star, qcf2_star, qrain_star     &
        ,qgraup_star, cf_star, cfl_star, cff_star, R_u, R_v                &
        ,a_realhd(rh_energy_corr), net_flux, net_mflux                     &
        ,murk, flash_pot,                                                  &
         dust_div1, dust_div2, dust_div3, dust_div4, dust_div5             &
        ,dust_div6, so2, so4_aitken,so4_accu, so4_diss, nh3, soot_new      &
        ,soot_agd, soot_cld, bmass_new, bmass_agd, bmass_cld, ocff_new,    &
         ocff_agd, ocff_cld, nitr_acc, nitr_diss, co2, tracer, tracer_ukca &
        ,biogenic, a_inthd(23), ukca_radaer                                &
  ! OUT Fields
        ,     ls_rain, ls_rainfrac, ls_snow, ls_graup, micro_tends         &
        ,     dryrho, photosynth_act_rad, rad_hr, dolr, SW_tile            &
  ! error information
        ,     ErrorStatus  )

    ! temporarily store r_u etc in the _p2 arrays.
    ! the _p2 arrays are not used at this point in the code, only
    ! subsequently to eg_f1sp, at which point the stored value is no
    ! longer required.
!$OMP  PARALLEL DEFAULT(NONE)                                          &
!$OMP& SHARED( udims_s, vdims_s, wdims, r_u_p2, r_u, r_v_p2, r_v,      &
!$OMP&         r_w_p2, r_w )                                           &
!$OMP& PRIVATE( i, j, k )
!$OMP DO SCHEDULE(STATIC)
  DO k = udims_s%k_start, udims_s%k_end
    DO j = udims_s%j_start, udims_s%j_end
      DO i = udims_s%i_start, udims_s%i_end
        r_u_p2(i,j,k) = r_u(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO k = vdims_s%k_start, vdims_s%k_end
    DO j = vdims_s%j_start, vdims_s%j_end
      DO i = vdims_s%i_start, vdims_s%i_end
        r_v_p2(i,j,k) = r_v(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO k = wdims%k_start, wdims%k_end
    DO j = wdims%j_start, wdims%j_end
      DO i = wdims%i_start, wdims%i_end
        r_w_p2(i,j,k) = r_w(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL

  ! set zero level star quantities for Endgame
  ! currently set equal to level 1
  CALL set_star_zero_level(                                              &
             theta_star,                                                 &
             q_star,                                                     &
             qcl_star,                                                   &
             qcf_star,                                                   &
             cf_star,                                                    &
             cfl_star,                                                   &
             cff_star,                                                   &
             qcf2_star,                                                  &
             qrain_star,                                                 &
             qgraup_star,                                                &
             L_mcr_qgraup,                                               &
             L_mcr_qrain,                                                &
             L_mcr_qcf2)

  ! ***********************************************************************
  !
  !  CONVERSION FROM ND to EG prognostics! (start)
  IF (ltimer) CALL timer('AS CONVERT',5)

    ! compute mixing ratio star fields from q star fields
  CALL update_m_star()

  !  CONVERSION FROM ND to EG  (end)
  !
  ! ***********************************************************************

  IF (ltimer) CALL timer('AS CONVERT',6)
  IF (ltimer .OR. lstashdumptimer) CALL timer('AS STASH',5)

  ! DEPENDS ON: Atm_Step_stash
  CALL atm_step_stash( &
       errorstatus, 1)

  IF (ltimer .OR. lstashdumptimer) CALL timer('AS STASH',6)

  IF (l_pv_tracer .OR. l_diab_tracer)                                   &
    CALL dyn_tr_slow ( dryrho, exner_theta_levels, dPV_rad, dPV_sw,     &
                       dPV_lw, dPV_mic, dPV_gwd, dPV_ph1, adv_only_PV,  &
                       dtheta_0, dtheta_rad, dtheta_SW, dtheta_LW,      &
                       dtheta_mic, dtheta_slow)


ELSE  ! L_physics =.false.

  ! initialise arrays that hold physics increments to zero
  ! DEPENDS ON: Atm_Step_phys_init
  CALL atm_step_phys_init( &
    r_v, r_u, theta_star, q_star, qcl_star, qcf_star,  &
    cf_star, cfl_star, cff_star, exner_lbc_real_tend,  &
    w_lbc_real_tend, errorstatus, 'zeroincs')

  ! initial NULL change if NOT running any physics code. 
  ! atmos_physics1 creates tendencies (star variables). So if nothing 
  ! happens all tendencies for computation of the source terms are zero.
  ! This is a "unnecessary" memory copy _provided_ that DR's
  ! understanding of the physics increments is correct.
   
!$OMP  PARALLEL DEFAULT(NONE)                                          &
!$OMP& SHARED( tdims_s, tdims, theta_star, q_star, qcl_star,           &
!$OMP&         qcf_star, qcf2_star, qrain_star, qgraup_star,           &
!$OMP&         l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,                  &
!$OMP&     m_star,mcl_star,mcf_star,mcf2_star,mrain_star,mgraup_star)  &
!$OMP& PRIVATE( i, j, k )
!$OMP DO SCHEDULE(STATIC)
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        theta_star (i,j,k) = 0.0
        q_star  (i,j,k)    = 0.0
        qcl_star(i,j,k)    = 0.0
        qcf_star(i,j,k)    = 0.0
        m_star(i,j,k)      = 0.0
        mcl_star(i,j,k)    = 0.0
        mcf_star(i,j,k)    = 0.0
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

  IF (l_mcr_qcf2  )  THEN
!$OMP DO SCHEDULE(STATIC)
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qcf2_star  (i,j,k) = 0.0
          mcf2_star  (i,j,k) = 0.0
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

  IF (l_mcr_qrain )  THEN
!$OMP DO SCHEDULE(STATIC)
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qrain_star (i,j,k) = 0.0
          mrain_star (i,j,k) = 0.0
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

  IF (l_mcr_qgraup)  THEN
!$OMP DO SCHEDULE(STATIC)
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qgraup_star(i,j,k) = 0.0
          mgraup_star(i,j,k) = 0.0
        END DO
      END DO
    END DO
!$OMP END DO
  END IF
!$OMP END PARALLEL

  CALL eg_mix_to_q                                                         &
       (tdims_l,tdims_s,                                                   &
       m_v, m_cl, m_cf,                                                    &
       m_cf2, m_r, m_gr,                                                   &
       l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,                              &
       q, qcl, qcf,                                                        &
       qcf2, qrain, qgraup)

  CALL swap_bounds(q,                                                    &
                   tdims_l%i_len - 2*tdims_l%halo_i,                     &
                   tdims_l%j_len - 2*tdims_l%halo_j,                     &
                   tdims_l%k_len,                                        &
                   tdims_l%halo_i, tdims_l%halo_j,                       &
                   fld_type_w,swap_field_is_scalar)
  CALL swap_bounds(theta,                                                &
                   tdims_s%i_len - 2*tdims_s%halo_i,                     &
                   tdims_s%j_len - 2*tdims_s%halo_j,                     &
                   tdims_s%k_len,                                        &
                   tdims_s%halo_i, tdims_s%halo_j,                       &
                   fld_type_w,swap_field_is_scalar)
  CALL swap_bounds(u,                                                    &
                   udims_s%i_len - 2*udims_s%halo_i,                     &
                   udims_s%j_len - 2*udims_s%halo_j,                     &
                   udims_s%k_len,                                        &
                   udims_s%halo_i, udims_s%halo_j,                       &
                   fld_type_u,swap_field_is_vector)
  CALL swap_bounds(v,                                                    &
                   vdims_s%i_len - 2*vdims_s%halo_i,                     &
                   vdims_s%j_len - 2*vdims_s%halo_j,                     &
                   vdims_s%k_len,                                        &
                   vdims_s%halo_i, vdims_s%halo_j,                       &
                   fld_type_v,swap_field_is_vector)


  IF (l_idealised_data .OR. problem_number /= standard) THEN
  CALL eg_idl_forcing(                                                     &
  ! IN Data Fields.
       u, v, theta, exner_theta_levels, exner, p_theta_levels,             &
  ! IN/OUT
       theta_star,q, q_star, qcl_star, qcf_star, qcf2_star, qrain_star,    &
       qgraup_star , dryrho,  r_u, r_v, wetrho_r_sq_n,                     &
       t_surface, problem_number,                                          &
  ! error information
       errorstatus  )
  END IF

!$OMP  PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(tdims,q_star,q)
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        q_star(i,j,k) = q_star(i,j,k) + q(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  CALL eg_q_to_mix                                                         &
       (tdims,tdims,                                                       &
        q_star(tdims%i_start:tdims%i_end,                                  &
               tdims%j_start:tdims%j_end,                                  &
               tdims%k_start:tdims%k_end),                                 &
        qcl   (tdims%i_start:tdims%i_end,                                  &
               tdims%j_start:tdims%j_end,                                  &
               tdims%k_start:tdims%k_end),                                 &
        qcf   (tdims%i_start:tdims%i_end,                                  &
               tdims%j_start:tdims%j_end,                                  &
               tdims%k_start:tdims%k_end),                                 &
        qcf2  (tdims%i_start:tdims%i_end,                                  &
               tdims%j_start:tdims%j_end,                                  &
               tdims%k_start:tdims%k_end),                                 &
        qrain (tdims%i_start:tdims%i_end,                                  &
               tdims%j_start:tdims%j_end,                                  &
               tdims%k_start:tdims%k_end),                                 &
        qgraup(tdims%i_start:tdims%i_end,                                  &
               tdims%j_start:tdims%j_end,                                  &
               tdims%k_start:tdims%k_end),                                 &
       l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,                              &
       m_star (tdims%i_start:tdims%i_end,                                  &
               tdims%j_start:tdims%j_end,                                  &
               tdims%k_start:tdims%k_end),                                 &
       m_cl   (tdims%i_start:tdims%i_end,                                  &
               tdims%j_start:tdims%j_end,                                  &
               tdims%k_start:tdims%k_end),                                 &
       m_cf   (tdims%i_start:tdims%i_end,                                  &
               tdims%j_start:tdims%j_end,                                  &
               tdims%k_start:tdims%k_end),                                 &
       m_cf2  (tdims%i_start:tdims%i_end,                                  &
               tdims%j_start:tdims%j_end,                                  &
               tdims%k_start:tdims%k_end),                                 &
       m_r    (tdims%i_start:tdims%i_end,                                  &
               tdims%j_start:tdims%j_end,                                  &
               tdims%k_start:tdims%k_end),                                 &
       m_gr   (tdims%i_start:tdims%i_end,                                  &
               tdims%j_start:tdims%j_end,                                  &
               tdims%k_start:tdims%k_end), swap_in=.FALSE. )

! swap bounds need to happen here as we aren't doing them in q_to_mix
CALL swap_bounds(m_cl,                                                     &
                 tdims_s%i_end - tdims_s%i_start+1 - 2*tdims_s%halo_i,     &
                 tdims_s%j_end - tdims_s%j_start+1 - 2*tdims_s%halo_j,     &
                 tdims_s%k_end - tdims_s%k_start + 1,                      &
                 tdims_s%halo_i, tdims_s%halo_j,                           &
                 fld_type_p, swap_field_is_scalar)
CALL swap_bounds(m_cf,                                                     &
                 tdims_s%i_end - tdims_s%i_start+1 - 2*tdims_s%halo_i,     &
                 tdims_s%j_end - tdims_s%j_start+1 - 2*tdims_s%halo_j,     &
                 tdims_s%k_end - tdims_s%k_start + 1,                      &
                 tdims_s%halo_i, tdims_s%halo_j,                           &
                 fld_type_p, swap_field_is_scalar)
IF (l_mcr_qcf2  )                                                          &
  CALL swap_bounds(m_cf2,                                                  &
                   tdims_s%i_end - tdims_s%i_start+1 - 2*tdims_s%halo_i,   &
                   tdims_s%j_end - tdims_s%j_start+1 - 2*tdims_s%halo_j,   &
                   tdims_s%k_end - tdims_s%k_start + 1,                    &
                   tdims_s%halo_i, tdims_s%halo_j,                         &
                   fld_type_p, swap_field_is_scalar)

IF (l_mcr_qrain )                                                          &
  CALL swap_bounds(m_r,                                                    &
                   tdims_s%i_end - tdims_s%i_start+1 - 2*tdims_s%halo_i,   &
                   tdims_s%j_end - tdims_s%j_start+1 - 2*tdims_s%halo_j,   &
                   tdims_s%k_end - tdims_s%k_start + 1,                    &
                   tdims_s%halo_i, tdims_s%halo_j,                         &
                   fld_type_p, swap_field_is_scalar)
IF (l_mcr_qgraup)                                                          &
  CALL swap_bounds(m_gr,                                                   &
                   tdims_s%i_end - tdims_s%i_start+1 - 2*tdims_s%halo_i,   &
                   tdims_s%j_end - tdims_s%j_start+1 - 2*tdims_s%halo_j,   &
                   tdims_s%k_end - tdims_s%k_start + 1,                    &
                   tdims_s%halo_i, tdims_s%halo_j,                         &
                   fld_type_p, swap_field_is_scalar)

 
!$OMP  PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k)                            &
!$OMP  SHARED(tdims,q_star,q,m_star,m_v)
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        q_star(i,j,k) = q_star(i,j,k) - q(i,j,k)
        m_star(i,j,k) = m_star(i,j,k) - m_v(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

END IF ! L_Physics

! Deallocate the arrays of UKCA_RADAER
DEALLOCATE(ukca_radaer%mix_ratio)
DEALLOCATE(ukca_radaer%comp_vol)
DEALLOCATE(ukca_radaer%dry_diam)
DEALLOCATE(ukca_radaer%wet_diam)
DEALLOCATE(ukca_radaer%modal_rho)
DEALLOCATE(ukca_radaer%modal_wtv)
DEALLOCATE(ukca_radaer%modal_vol)
DEALLOCATE(ukca_radaer%modal_nbr)

! Deallocate arrays for UKCA_CDNC
DEALLOCATE(ukca_cdnc%cdnc)
DEALLOCATE(ukca_cdnc%cdnc3)

! Deallocate arrays for EasyAerosol
CALL deallocate_easyaerosol_rad(easyaerosol_sw)
CALL deallocate_easyaerosol_rad(easyaerosol_lw)
CALL deallocate_easyaerosol_cdnc(easyaerosol_cdnc)

! Deallocate the array of the aerosol climatology for NWP
DEALLOCATE(arcl)

IF (.NOT. couple_app) THEN
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE ( i, j, k ) SCHEDULE(STATIC) &
!$OMP SHARED(model_levels,udims_s,r_u,vdims_s,r_v)

  DO k = udims_s%k_start, udims_s%k_end
    DO j = udims_s%j_start, udims_s%j_end
      DO i = udims_s%i_start, udims_s%i_end
        r_u(i,j,k) = 0.0
      END DO
    END DO
    DO j = vdims_s%j_start, vdims_s%j_end
      DO i = vdims_s%i_start, vdims_s%i_end
        r_v(i,j,k) = 0.0
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO 

END IF

IF (printstatus >=  prstatus_diag) CALL diag_R()

! compute the slow physics source terms for theta and moisture
IF (couple_app) THEN
  CALL eg_r( theta_star, m_star,mcl_star,mcf_star,mgraup_star, &
             mrain_star,mcf2_star, m_v, theta )
END IF
l_call_from_f1sp_in = .FALSE.
CALL eg_sisl_init(                                                    &
         row_length,rows,n_rows,model_levels, l_inc_solver,           &
         l_call_from_solver,l_call_from_f1sp_in, ih,g_theta, u, v, w, &
         thetav, dryrho, m_v, m_cl, m_cf, m_r, m_gr, m_cf2,exner,     &
         exner_surf, r_u, r_v, r_w, r_theta, r_rho, r_m_v, r_m_cl,    &
         r_m_cf, r_m_r, r_m_gr, r_m_cf2, etadot,  psi_w_surf,psi_w_lid)

IF (ltimer) CALL timer('AS Atmos_Phys1 (AP1)',6)

IF (l_tracer .AND. errorstatus == 0) THEN

  ! store physics changes
  CALL tr_set_phys_4A(super_array_size, super_tracer_phys1,           &
                   l_co2_interactive, co2,                            &
                   l_murk_advect, murk,                               &
                   l_soot, soot_new, soot_agd, soot_cld,              &
                   l_sulpc_so2, so2, so4_aitken,                      &
                   so4_accu, so4_diss, l_sulpc_nh3, nh3,              &
                   l_sulpc_dms, dms,                                  &
                   l_dust, dust_div1, dust_div2,  dust_div3,          &
                           dust_div4, dust_div5,  dust_div6,          &
                   l_biomass, bmass_new, bmass_agd, bmass_cld,        &
                   l_ocff, ocff_new, ocff_agd, ocff_cld,              &
                   l_nitrate, nitr_acc, nitr_diss,                    &
                   l_use_cariolle, ozone_tracer,                      &
                   tracer, tracer_ukca,                               &
                   row_length, rows,                                  &
                   model_levels, tr_levels, tr_vars, tr_ukca,         &
                   set_halos, tdims_l )

END IF  ! L_tracer and ErrorStatus == 0

! store physics changes for use in Sl_moist_conserve
! only some pc2 code in atm_step_alloc_4A is needed
! DEPENDS ON: Atm_Step_alloc_4A
CALL atm_step_alloc_4A(                                                 &
         cf_star, cfl_star, cff_star, q_star, qcl_star, qcf_star,       &
         frac_control, r_u, r_v, r_w, errorstatus, 'alloc_mr' )

! DEPENDS ON: Atm_Step_phys_reset
CALL atm_step_phys_reset(                                               &
    q_star, qcl_star, qcf_star, r_v, r_u, dolr, theta_star, gs1,        &
    cf_star, cfl_star, cff_star, tstar_tile, snodep_tile, 'cyclrset')

!=== Polar filter + diffusion of increments section ==================================
! ----------------------------------------------------------------------
! Section 0.4  Filter winds and theta near poles if active
!              Do horizontal diffusion as a filter if active
! ----------------------------------------------------------------------

! Call timer for diffusion code

IF ( L_polar_filter_incs .OR. L_filter_incs ) THEN
  IF (ltimer) CALL timer('AS Filter',5)

  CALL filter_diag_printing3()

  CALL eg_ni_filter_incs_ctl(                                       &
                        r_theta, r_u, r_v, r_w, r_rho,              &
                        s_thetav, s_u, s_v, s_w,                    &
                        thetav, u, v, w, dryrho, exner,             &
                        row_length, rows, n_rows, model_levels,     &
                        r_theta_levels, r_rho_levels,               &
                        r_at_u, r_at_v,                             &
                        max_121_rows, u_sweeps, v_sweeps,           &
                        global_u_filter, global_v_filter,           &
                        u_begin, u_end, v_begin, v_end,             &
                        diff_coeff_phi, diff_coeff_u, diff_coeff_v, &
                        first_constant_r_rho_level,                 &
                        first_constant_r_rho_level_m1,              &
                        horizontal_level,                           &
                        offx, offy, halo_i, halo_j,                 &
                        nproc_y, at_extremity,                      &
                        L_polar_filter_incs, L_diff_incs,           &
                        L_pofil_hadgem2,                            &
                        l_eliminate_rho, 0,                         &
                        xi1_u, xi1_p, xi2_p,                        &
                        pole_consts, gc_proc_row_group,             &
                        global_row_length,                          &
                        csxi2_v, csxi2_p)

  CALL filter_diag_printing4()

  IF (ltimer) CALL timer('AS Filter',6)
END IF    ! L_polar_filter_incs .or. L_filter_incs

!=== End Polar filter + diffusion of increments section ===================

! ----------------------------------------------------------------------
! Section 0.5 Calculation of coefficients in turbulence scheme
! ----------------------------------------------------------------------

! DEPENDS ON: Atm_Step_phys_init
CALL atm_step_phys_init( &
  r_v, r_u, theta_star, q_star, qcl_star, qcf_star,  &
  cf_star, cfl_star, cff_star, exner_lbc_real_tend,  &
  w_lbc_real_tend, errorstatus, 'turb_cof')

!------------------------------------------------------------------------
! Set hydrostatic balance in LBC regions
!------------------------------------------------------------------------

IF (model_type == mt_lam) THEN
  IF (ltimer) CALL timer('AS LAM_LBCS',5)

  IF (L_fixed_lbcs) THEN
    L_update_lbcs = .FALSE.
  ELSE
    L_update_lbcs = .TRUE.
  END IF !  L_fixed_lbcs

  L_do_halos=.FALSE.
  L_do_boundaries=.TRUE.

  IF (rim_stepsa  ==  0) THEN
    increment_factor = 0.0
    L_balance = .FALSE.
  ELSE IF ( MOD(Timestep_Number-1,rim_stepsa) == 0) THEN
    L_balance = .TRUE.
    increment_factor = 1.0 / rim_stepsa
  ELSE
    L_balance = .FALSE.
    increment_factor = 1.0 /                                     &
               (rim_stepsa - MOD(Timestep_Number-1,rim_stepsa))
  END IF

  lbc_size=lenrima(fld_type_p,halo_type_extended, rima_type_norm)

  IF ( L_update_lbcs ) THEN
    ! Obtain Exner tendency to pass to the solver to apply lbc

    ! Apply additional balance to EXNER_LBC_TEND, RHO_LBC_TEND
    ! and W_LBC_TEND only on timestep after up_bound has been called.

    IF (L_LBC_balance .AND. L_balance) THEN

      CALL eg_balance_lbc_values(                                &
     exner_lbc_tend, rho_lbc_tend, theta_lbc_tend,               &
     q_lbc_tend, w_lbc_tend, w_adv_lbc_tend,                     &
     u_lbc_tend, v_lbc_tend,                                     &
     qcl_lbc_tend, qcf_lbc_tend, qcf2_lbc_tend,                  &
     qrain_lbc_tend, qgraup_lbc_tend,                            &
     L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,                      &
     r_rho_levels, r_theta_levels,                               &
     row_length, rows, halo_i,halo_j,                            &
     lenrima(fld_type_p,halo_type_extended,rima_type_norm),      &
     lenrima(fld_type_u,halo_type_extended,rima_type_norm),      &
     lenrima(fld_type_v,halo_type_extended,rima_type_norm),      &
     lbc_starta(1,fld_type_p,halo_type_extended,rima_type_norm), &
     lbc_starta(1,fld_type_u,halo_type_extended,rima_type_norm), &
     lbc_starta(1,fld_type_v,halo_type_extended,rima_type_norm), &
     rimwidtha(rima_type_norm), n_rims_to_do,rimweightsa,        &
     at_extremity,                                               &
     delta_phi, delta_lambda,                                    &
     base_phi, base_lambda,                                      &
     datastart,                                                  &
     lat_rot_np, global_row_length, global_rows)
    END IF  !  L_LBC_balance and L_balance

    DO k = 1, model_levels
      DO i=1,lbc_size
        exner_lbc_real_tend(i,k) = increment_factor *              &
                         ( exner_lbc_tend(i,k) - exner_lbc(i,k) )
      END DO
    END DO

  END IF ! L_update_lbcs

  IF (ltimer) CALL timer('AS LAM_LBCS',6)

END IF  ! model_type

!------------------------------------------------------------------------
! End - Set hydrostatic balance in LBC regions
!------------------------------------------------------------------------

! Create Dynamics Advection STASHwork array as it is required later in
! the code as well
IF ( sf(0,12) ) THEN
  ALLOCATE (stashwork12(stash_maxlen(12,atmos_im)))
END IF

! ----------------------------------------------------------------------
! Section 1.0 Start outer loop
! ----------------------------------------------------------------------

IF (total_conv_outer) numcycles = 999

exner_res = 1.0 ! initialise to a "large" value
cycleno=0

outerloop : DO WHILE  (cycleno < numcycles)

  cycleno=cycleno+1

  IF (total_conv_outer .AND. exner_res < 2.0d-3) cycleno = numcycles

  !-----------------------------------------------------------------------
  ! Restore phys1 variables to be used as predictors
  ! DEPENDS ON: Atm_Step_phys_reset
  CALL atm_step_phys_reset(                                             &
    q_star, qcl_star, qcf_star, r_v, r_u, dolr, theta_star, gs1,        &
    cf_star, cfl_star, cff_star, tstar_tile, snodep_tile, 'phy1rest')

  IF ( (L_subfilter_horiz .OR. L_subfilter_vert) .AND.                    &
       (cycleno == 1 .OR. .NOT. l_quick_ap2)    ) THEN
    ! if quick running, only call on 1st cycle

    ! Calculate lambda^2 and S in eg_turb_smagorinsky
    ! DEPENDS ON: eg_turb_smagorinsky
    CALL eg_turb_smagorinsky(                                           &
                             u, v, w,                                   &
                             model_levels,                              &
                             r_theta_levels, r_rho_levels )

  END IF  !   L_subfilter_horiz or L_subfilter_vert and l_quick_ap2

  ! Initialise advection increments in first cycle only
  IF ( cycleno == 1 ) THEN
    IF (sf(0,12)) THEN
      ! here we use the stored value of r_u etc (in the _p2 arrays) which
      ! have not seen the addition of dynamics increments
      CALL adv_incs_init(r_u_p2, r_v_p2, theta_star,                       &
           q_star, qcl_star, qcf_star, qrain_star, qgraup_star, qcf2_star, &
           m_star, mcl_star, mcf_star, mrain_star, mgraup_star, mcf2_star, &
           cf_star, cfl_star, cff_star)
    END IF
  END IF

  ! ----------------------------------------------------------------------
  ! Calculate time level n advective momentum quantities.
  ! ----------------------------------------------------------------------

  ! update LBC's for n+1 variables
  IF (model_type == mt_lam) THEN
    !
    ! Currently forced to match settings for the standard EG SEUKV job.
    ! This will change at a later date when enough LAM tests have been performed
    ! and the best settings have been determined!
    !
    IF ( cycleno > 1) THEN
      ! DEPENDS ON: init_lbc_dynamics
      CALL init_lbc_dynamics(u_np1,v_np1,w_np1, etadot_np1,                 &
                             u_lbc, v_lbc, w_lbc,                           &
                             u_lbc_tend, v_lbc_tend, w_lbc_tend,            &
                             increment_factor, rim_stepsa, timestep_number, &
                             row_length, rows, n_rows, model_levels,        &
                             offx, offy, halo_i, halo_j,                    &
                             lenrima(1,1,rima_type_norm),                   &
                             rimwidtha(rima_type_norm), rimweightsa,        &
                             lbc_sizea(1,1,1,rima_type_norm),               &
                             lbc_startA(1,1,1,rima_type_norm),              &
                             at_extremity,                                  &
                             L_do_boundaries, L_do_halos                    &
                            )
    END IF
  END IF

  IF (ltimer) CALL timer('AS S-L Advect (AA)',5)
  IF (ltimer) CALL timer('AA SL_Full_Wind',5)

  CALL eg_sl_full_wind(g_i_pe,depart_scheme,                            &
               depart_order, high_order_scheme(wind_sl),                &
               monotone_scheme(wind_sl),                                &
               ritchie_high_order_scheme, ritchie_monotone_scheme,      &
               first_constant_r_rho_level,                              &
               interp_vertical_search_tol, check_bottom_levels,         &
               l_shallow,                                               &
               l_rk_dps,                                                &
               l_high(wind_sl),  l_mono(wind_sl),  l_ritchie_high,      &
               l_ritchie_mono, lam_max_cfl,                             &
               etadot, u_np1, v_np1, w_np1, etadot_np1,                 &
               u_adv, v_adv, w_adv, r_u, r_v, r_w,                      &
               r_u_d, r_v_d, r_w_d,                                     &
               errorstatus )

  IF (ltimer) CALL timer('AA SL_Full_Wind',6)

  L_tracer_1_if :  IF ( L_tracer .AND. CycleNo == NumCycles ) THEN  !it1

    IF (ltimer) CALL timer('AA SL_Tracer',5)

    CALL sl_tracer1_4A(super_tracer_phys1,                            &
                     super_array_size, eta_theta_levels,              &
                     row_length, rows, n_rows, model_levels,          &
                     g_i_pe, g_j_pe,                                  &
                     high_order_scheme(tracer_SL),                    &
                     monotone_scheme(tracer_SL),                      &
                     l_high(tracer_SL), l_mono(tracer_SL),            &
                     depart_xi1_w, depart_xi2_w,  depart_xi3_w,       &
                     co2, l_co2_interactive,                          &
                     murk, l_murk_advect,                             &
                     dust_div1,dust_div2,                             &
                     dust_div3,dust_div4,                             &
                     dust_div5,dust_div6, l_dust,                     &
                     soot_new, soot_agd,                              &
                     soot_cld, l_soot,                                &
                     bmass_new, bmass_agd,                            &
                     bmass_cld, l_biomass,                            &
                     ocff_new, ocff_agd, ocff_cld, l_ocff,            &
                     so2, so4_aitken,                                 &
                     so4_accu,                                        &
                     so4_diss, nh3, dms,                              &
                     l_sulpc_so2, l_sulpc_nh3, l_sulpc_dms,           &
                     nitr_acc, nitr_diss, l_nitrate,                  &
                     tracer, tr_vars,                                 &
                     tracer_ukca, tr_ukca,                            &
                     l_use_cariolle, ozone_tracer,                    &
                     dryrho, dryrho_np1, errorstatus)

    IF (ltimer) CALL timer('AA SL_Tracer',6)

  END IF  L_tracer_1_if

  IF (ltimer) CALL timer('AA SL_Thermo',5)

  CALL eg_sl_thermo(g_i_pe,                                             &
                l_inc_solver,high_order_scheme(theta_sl),               &
                monotone_scheme(theta_sl),l_high(theta_sl),             &
                l_mono(theta_sl), r_theta,etadot_np1,                   &
                r_theta_d, dryrho, dryrho_np1, errorstatus )

  IF (ltimer) CALL timer('AA SL_Thermo',6)
  IF (ltimer) CALL timer('AA SL_Rho',5)

  CALL eg_sl_rho(g_i_pe,                                                &
                l_inc_solver, depart_scheme, depart_order,              &
                high_order_scheme(rho_sl), monotone_scheme(rho_sl),     &
                ritchie_high_order_scheme, ritchie_monotone_scheme,     &
                first_constant_r_rho_level,interp_vertical_search_tol,  &
                check_bottom_levels, l_shallow,                         &
                l_rk_dps,                                               &
                l_high(rho_sl),l_mono(rho_sl), l_ritchie_high,          &
                l_ritchie_mono, lam_max_cfl,                            &
                etadot, u_np1, v_np1, w_np1,etadot_np1, u_adv, v_adv,   &
                w_adv, r_rho, r_rho_d, errorstatus)

  IF (ltimer) CALL timer('AA SL_Rho',6)
  IF (ltimer) CALL timer('AA SL_Moisture',5)

  IF (.NOT. l_dry) THEN
    CALL eg_sl_moisture( moisture_array_size,                           &
                row_length, rows, n_rows, model_levels, halo_i,         &
                halo_j, offx, offy, datastart, g_i_pe,                  &
                high_order_scheme(moist_sl),monotone_scheme(moist_sl),  &
                l_high(moist_sl),l_mono(moist_sl),L_mcr_qrain,          &
                l_mcr_qcf2, l_mcr_qgraup,                               &
                r_m_v, r_m_cl, r_m_cf,r_m_r, r_m_gr, r_m_cf2,           &
                cf_bulk, cf_liquid, cf_frozen,exner_theta_levels,       &
                exner_star,r_m_v_d, r_m_cl_d, r_m_cf_d,                 &
                r_m_r_d, r_m_gr_d, r_m_cf2_d,                           &
                cf_star, cfl_star, cff_star,                            &
                conv_prog_1, conv_prog_2, conv_prog_3, conv_prog_precip,&
                dryrho, dryrho_np1, errorstatus)

    IF ( l_casim .AND. n_casim_progs > 0                                &
         .AND. cycleno == numcycles ) THEN

      ! CASIM microphysics active, with additional prognostics not
      ! available to non-CASIM model runs. Advect prognostics on the
      ! last ENDGame cycle only, because these prognostics are not
      ! required in atmos_physics2. Treat as moist prognostic variables
      ! to ensure they are advected in a consistent manner to mass
      ! mixing ratios.

      ! The number of prognostics that are advected by this subroutine
      ! is n_casim_progs, in the module casim_switches.

      CALL eg_sl_casim(                                                  &
              row_length, rows, model_levels, halo_i, halo_j, datastart, &
              g_i_pe, high_order_scheme(moist_sl),                       &
              monotone_scheme(moist_sl), l_high(moist_sl),               &
              l_mono(moist_sl), errorstatus)

    END IF ! l_casim etc

  END IF

  IF (ltimer) CALL timer('AA SL_Moisture',6)

  IF (i_bl_vn == i_bl_vn_1a .AND. l_adv_turb_field) THEN
    CALL eg_sl_turb(                                                    &
                    turb_array_size,                                    &
                    row_length, rows, halo_i, halo_j, datastart, g_i_pe,&
                    e_trb, tsq_trb, qsq_trb, cov_trb,                   &
                    errorstatus)
  END IF

  ! PV-TRACER: Perform interpolation to departure point
  IF ( CycleNo == NumCycles .AND. (l_pv_tracer .OR. l_diab_tracer ))  &
    CALL eg_sl_dyn_tr(                                                &
             row_length, rows, n_rows, model_levels,                  &
             g_i_pe,                                                  &
             high_order_scheme(theta_sl), monotone_scheme(theta_sl),  &
             l_high(theta_sl), l_mono(theta_sl),                      &
             dPV_rad, dPV_sw, dPV_lw,                                 &
             dPV_mic, dPV_gwd, dPV_ph1,                               &
             dPV_conv, dPV_bl, dPV_stph, dPV_cld,                     &
             dPV_iau, dPV_nud, dPV_tot,                               &
             dPV_adv, dPV_sol, dPV_mass, adv_only_PV,                 &
             dtheta_0, dtheta_bl, dtheta_bl_mix, dtheta_bl_LH,        &
             dtheta_conv, dtheta_mic, dtheta_rad,                     &
             dtheta_SW, dtheta_LW, dtheta_slow, dtheta_cld,           &
             errorstatus)

  IF (ltimer) CALL timer('AS S-L Advect (AA)',6)

  ! ----------------------------------------------------------------------
  ! Section 3.0  Call Atmospheric Physics2
  ! ----------------------------------------------------------------------
  IF (ltimer) CALL timer('AS Atmos_Phys2 (AP2)',5)

  IF (.NOT. l_run_with_physics2) THEN
    l_physics_store=l_physics
    l_physics=.FALSE.
  END IF

  IF (l_tracer .AND. cycleno == numcycles ) THEN
    !
    !  super_tracer_phys2 = tracers ( before physic2 )
    !
    CALL tr_set_phys_4A(super_array_size, super_tracer_phys2,   &
                 l_co2_interactive, co2, l_murk_advect, murk,   &
                 l_soot, soot_new, soot_agd,soot_cld,           &
                 l_sulpc_so2, so2,so4_aitken,so4_accu,so4_diss, &
                 l_sulpc_nh3, nh3,l_sulpc_dms, dms,             &
                 l_dust, dust_div1,dust_div2,dust_div3,         &
                         dust_div4,dust_div5,dust_div6,         &
                 l_biomass, bmass_new, bmass_agd, bmass_cld,    &
                 l_ocff, ocff_new, ocff_agd, ocff_cld,          &
                 l_nitrate, nitr_acc, nitr_diss,                &
                 l_use_cariolle, ozone_tracer,                  &
                 tracer, tracer_ukca,                           &
                 row_length, rows,                              &
                 model_levels, tr_levels, tr_vars, tr_ukca,     &
                 do_not_set_halos, tdims )

  END IF  ! L_tracer and CycleNo == NumCycles

  !        F1SP computes the predictors for u,v,w aka R_u, R_v, theta_star
  !        and moisture fields

  IF ( l_inc_solver) THEN
    CALL eg_f1sp_inc(theta_star ,m_star ,mcl_star,mcf_star              &
                      ,mcf2_star,mrain_star,mgraup_star                 &
                      ,u_np1,v_np1,w_np1,u,v,w, exner_np1               &
                      ,exner_surf_np1,dryrho_np1, thetav_np1, etadot_np1&
                      ,m_v_np1, m_cl_np1, m_cf_np1, m_r_np1             &
                      ,m_gr_np1, m_cf2_np1,   ih, l_slice               &
                      ,row_length,rows,n_rows,model_levels              &
                      ,L_mcr_qcf2,L_mcr_qrain,L_mcr_qgraup,             &
                       psi_w_surf,psi_w_lid)

  ELSE
    CALL eg_f1sp(theta_star ,m_star ,mcl_star,mcf_star                  &
                       ,mcf2_star,mrain_star,mgraup_star                &
                       ,u_np1,v_np1,w_np1,u,v,w, exner_np1              &
                       ,exner_surf_np1,dryrho_np1, thetav_np1,etadot_np1&
                       ,m_v_np1, m_cl_np1, m_cf_np1, m_r_np1            &
                       ,m_gr_np1, m_cf2_np1,   ih,l_slice               &
                       ,row_length,rows,n_rows,model_levels             &
                       ,L_mcr_qcf2,L_mcr_qrain,L_mcr_qgraup,            &
                        psi_w_surf,psi_w_lid)
  END IF

  ! ***********************************************************************
  !
  !  CONVERSION FROM EG to ND  (start)

  IF (ltimer) CALL timer('AS CONVERT',5)

  CALL update_q_star()

  IF (ltimer) CALL timer('AS CONVERT',6)
  !  CONVERSION FROM EG to ND  (end)
  !
  ! ***********************************************************************

  ! -------------------------------------------
  ! Section 2.3 Diagnostics at end of advection
  ! -------------------------------------------
  ! Apply diagnostics at final cycle only
  IF ( cycleno == numcycles ) THEN

    ! Section 12: 'dynamics advection' based quantities
    IF ( sf(0,12) .AND. ErrorStatus == 0 ) THEN

      IF (ltimer .OR. lstashdumptimer) CALL timer('AS STASH',5)

      CALL diagnostics_adv(                                             &
                    row_length, rows, n_rows,                           &
      ! primary wind fields:
                    u, v,                                               &
                    theta, q, qcl, qcf, qrain, qgraup, qcf2,            &
                    m_v, m_cl, m_cf, m_r, m_gr, m_cf2,                  &
                    cf_bulk, cf_liquid, cf_frozen,                      &
      ! wind field increments after advection:
                    r_u_p2, r_v_p2, r_w_p2,                             &
                    theta_star, q_star, qcl_star, qcf_star,             &
                    qrain_star, qgraup_star, qcf2_star,                 &
                    m_star, mcl_star, mcf_star,                         &
                    mrain_star, mgraup_star, mcf2_star,                 &
                    cf_star, cfl_star, cff_star,                        &
                    exner_theta_levels,                                 &
      ! departure points for w
                    depart_xi1_rho, depart_xi2_rho, depart_xi3_rho,     &
                    r_theta_levels,                                     &
                    stashwork12)

      CALL adv_incs_dealloc()
 
      IF (ltimer .OR. lstashdumptimer) CALL timer('AS STASH',6)

    END IF !   SF(0,12)

  END IF ! CycleNo == NumCycles

  !        Save star fields to obtain increments after call to
  !        Atmos_Physics2

  !        we do not need to save the moisture, because we have those in
  !        the ENDGame mixing ratio fields anyway!

  !        This section can very likely be replace by some of the logic in
  !        atm_step_local!

!$OMP  PARALLEL DEFAULT(NONE)                                          &
!$OMP& SHARED( udims_s, vdims_s, wdims, r_u_p2_n, r_u_p2, r_v_p2_n,    &
!$OMP&      r_v_p2, r_w_p2_n, r_w_p2, tdims, theta_star, theta_star_n) &
!$OMP& PRIVATE( i, j, k )
!$OMP DO SCHEDULE(STATIC)
  DO k = udims_s%k_start, udims_s%k_end
    DO j = udims_s%j_start, udims_s%j_end
      DO i = udims_s%i_start, udims_s%i_end
        r_u_p2_n(i,j,k) = r_u_p2(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO k = vdims_s%k_start, vdims_s%k_end
    DO j = vdims_s%j_start, vdims_s%j_end
      DO i = vdims_s%i_start, vdims_s%i_end
        r_v_p2_n(i,j,k) = r_v_p2(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO k = wdims%k_start, wdims%k_end
    DO j = wdims%j_start, wdims%j_end
      DO i = wdims%i_start, wdims%i_end
        r_w_p2_n(i,j,k) = r_w_p2(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        theta_star_n(i,j,k) = theta_star(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL

! Calculate idealised surface fluxes
  IF (l_idealised_data) CALL external_force_2(l_spec_z0, flux_e, flux_h,    &
                                              ustar_in,z0m_scm, z0h_scm)

  IF (ltimer) CALL timer('AS Atmos_Phys2 (AP2)',6)

  lphysics :  IF (l_physics) THEN
    IF (ltimer) CALL timer('AS Atmos_Phys2 (AP2)',5)

    ! DEPENDS ON: Atm_Step_diag
    CALL Atm_Step_diag(4)

    ! ensure this is always 0 before ap2 call
    zlcl_mixed = 0.0

    ! NB if you are changing the argument list to atmos_physics2, please
    ! do an equivalent change in routine scm_main to keep the single column
    ! model consistent. Note there are two calls to atmos_physics2 in scm_main.

    ! DEPENDS ON: atmos_physics2
    CALL atmos_physics2(                                           &
    ! Parallel variables
            global_row_length, global_rows, nproc, nproc_x, nproc_y          &
          , g_rows, g_row_length, NumCycles, CycleNo                         &
    ! field dimensions etc.
          ,     row_length, rows, n_rows, land_field                         &
          ,     bl_levels, st_levels, sm_levels, cloud_levels                &
          ,     land_ice_points, soil_points, n_cca_lev, ntiles, tr_levels   &
          ,     first_constant_r_rho_level, dim_cs1, dim_cs2                 &
    ! IN Substepping information and switches
          ,     L_dry, L_lbc_old, lcal360                                    &
    ! Model Parameters
          ,     RHcrit, co2_mmr, tr_vars, tr_ukca                            &
    ! Vertical coordinate levels.
          ,     dryrho                                                       &
          ,     delta_lambda, delta_phi                                      &
    ! The following 2 lines are not used by Endgame
          ,     gdlambda_p, dphi_p, wt_lambda_p, wt_lambda_u                 &
          ,     wt_phi_p, wt_phi_v                                           &
    ! The following line is used by Endgame
          ,     lat_rot_NP, long_rot_NP, f3_at_u         &
    ! Time stepping information
          , i_year, i_day_number, i_hour, i_minute, i_second,                &
    ! River routing
            aocpl_row_length,aocpl_p_rows,xpa,xua,xva,ypa,yua,yva,           &
            g_p_field,g_r_field,a_inthd(16),lasize(1,fld_type_r,             &
            halo_type_no_halo),lasize(2,fld_type_r,halo_type_no_halo),       &
            glsize(1,fld_type_r),glsize(2,fld_type_r),river_vel,             &
            river_mcoef,i_river_vn,riv_direction,riv_sequence,riv_storage,   &
    !  Add inland basin outflow to arguments
            riv_inlandatm,                                                   &
    !  Add lake evaporation:
            acc_lake_evap,                                                   &
    ! Grid-to-grid river routing
            riv_iarea,riv_slope, riv_flowobs1,riv_inext, riv_jnext,riv_land, &
            riv_substore, riv_surfstore,riv_flowin,riv_bflowin,              &

    ! diagnostic info
                STASHwork3,STASHwork5,STASHwork8,STASHwork9,STASHwork19,     &
                STASHwork26,                                                 &
    ! SCM Diagnostics (dummy values in full UM)
            nscmdpkgs, l_scmdiags, conv_mode, l_emcorr_opt,                  &
    ! Data Fields.
           theta, q, qcl, qcf, qrain, qgraup, qcf2,                          &
           Wetrho_r_sq_n, u, v, w, etadot, p ,pstar,                         &
           exner,exner_theta_levels, land, p_theta_levels,                   &
    ! Mixing ratio prognostics; passed in separately to the q fields so
    ! that any convection scheme which works with mixing ratios
    ! can just use them instead of having to convert the q-fields back to
    ! mixing ratios, which would seem perverse!
           m_v, m_cl, m_cf, m_cf2, m_r, m_gr,                                &

    ! ancillary fields and fields needed to be kept from timestep to
    ! timestep
           land_index, land_ice_index, soil_index, canopy_water              &
          ,snodep, therm_cond, therm_cap, vol_smc_crit                       &
          ,vol_smc_wilt, vol_smc_sat, z0m_soil, sthf, sthu                   &
          ,orog_sil,orog_ho2,orog_sd,ice_thickness,ice_fraction              &
          ,u_sea,v_sea,u_0_p,v_0_p,cca_dp,cca_md,cca_sh,cca,ccb              &
          ,cct, cclwp, ccw_rad, lcbase, deep_soil_temp, p_ti, ti, ice_k_cat  &
          ,tstar,z0, p_ice_fract, p_ice_thick                                &
          ,sat_soil_cond,sat_soilw_suction,clapp_horn                        &
          ,smcl, t1_sd, q1_sd, zh, ddmfx, cf_area, cf_bulk, cf_liquid        &
          ,cf_frozen, ls_rain, ls_rainfrac, ls_snow, ls_graup, micro_tends   &
          ,totalppn, photosynth_act_rad, rad_hr, soil_clay                   &
          ,soil_silt,soil_sand, dust_mrel1,dust_mrel2                        &
          ,dust_mrel3,dust_mrel4,dust_mrel5,dust_mrel6                       &
          ,so2_hilem, so2_em, nh3_em, dms_em,soot_hilem, soot_em, ocff_hilem &
          ,ocff_em, co2_emits, co2flux ,deep_flag, past_precip, past_conv_ht &

    ! IN/OUT
          , theta_star,q_star,qcl_star,qcf_star, qrain_star, qgraup_star     &
          , qcf2_star, cf_star,cfl_star,cff_star                             &
          , r_u_p2, r_v_p2, r_w_p2                                           &
          , net_flux, net_mflux, murk,tracer, tracer_ukca                    &
          , dust_div1, dust_div2, dust_div3 ,dust_div4, dust_div5, dust_div6 &
          , so2, dms, so4_aitken, so4_accu, so4_diss, nh3, soot_new          &
          , soot_agd, soot_cld, bmass_new, bmass_agd, bmass_cld, ocff_new    &
          , ocff_agd, ocff_cld, nitr_acc, nitr_diss, co2                     &
    ! IN/OUT River routing
          , tot_surfroff, tot_subroff                                        &
    ! OUT Fields
          , ntml, cumulus, nbdsc, ntdsc                                      &
          , rhcpt, rhc_row_length, rhc_rows                                  &
          , zlcl_mixed                                                       &
    ! Additional variables for MOSES II
          , frac_typ, disturb_veg_pointer, canht_pft, lai_pft_in             &
          , can_water_tile, catch_tile, catch_snow                           &
          , snow_grnd, snodep_tile, z0_tile, z0h_tile                        & 
          , tstar_tile, tsurf_elev_surft                                     &
          , infil_tile, rgrain_tile, cs, gs                                  &
          , co2_dim_row, co2_dim_len, a_inthd(23)                            &
          , stepim(atmos_im), g_lf_pft_acc, g_phlf_pft_acc, npp_pft_acc      &
          , rsp_w_pft_acc, rsa                                               &
          , land_pts_trif, npft_trif, dolr, LW_down, SW_tile                 &
          , frac_land,tstar_land,tstar_sea,p_tstar_sice,tstar_sice           &
          , soil_alb,cos_zenith_angle                                        &
    ! INOUT variables for TKE based turbulence schemes
          , e_trb, tsq_trb, qsq_trb, cov_trb, zhpar_shcu                     &
    ! OUT variable to store BL w-variance diagnostic
          , bl_w_var                                                         &
    ! IN/OUT convection prognostics
          , conv_prog_1, conv_prog_2, conv_prog_3, conv_prog_precip          &
          , ux_ccp, uy_ccp, um_ccp, g_ccp, h_ccp, riso_ccp, rdir_ccp         &
    ! Additional variables required for large-scale hydrology:
          , fexp,gamma_int,ti_mean,ti_sig, fsfc_sat,f_wetland                &
          , water_table,sthzw, a_fsat,c_fsat,a_fwet,c_fwet                   &
    ! IN/OUT JULES 2 prognostics
          , snowdepth, rho_snow_grnd, nsnow                                  &
          , ds, sice, sliq, tsnowlayer, rho_snow, rgrainl                    &
    ! IN/OUT FLake lake scheme prognostics
          , lake_depth, lake_fetch, lake_t_mean, lake_t_mxl                  &
          , lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape                  &
          , lake_g_dt                                                        &
    ! Cariolle ozone
          , ozone_tracer                                                     &

    ! Additional screen-level variables
          , TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans                           &

    ! Variables required for COSP (out)
          ,cosp_crain_3d, cosp_csnow_3d                                      &

    ! SCM and idealised UM surface forcing parameters
          ,flux_e,flux_h,ustar_in,L_spec_z0,z0m_scm,z0h_scm,                 &

    ! error information
                                ErrorStatus )

    IF ( ( .NOT. (l_rp2  .AND. i_rp_scheme == i_rp2b) ) .AND.                &
         cycleno == numcycles ) THEN
      ! LAI may be a prognostic - pass lai_pft_in
      ! back to lai_pft to preserve any changes made in the subroutines
      lai_pft(:,:) = lai_pft_in(:,:)
    END IF


    IF (l_tracer .AND. cycleno == numcycles ) THEN
      !
      ! set super_tracer_phys2 such that it contains only the atmos_physics2 contribution.
      !     super_tracer_phys2 = tracers - super_tracer_phys2(before atmos_physics2)
      !
      CALL tr_reset_4A (super_array_size, super_tracer_phys2,    &
           l_co2_interactive, co2,  l_murk_advect, murk,         &
           l_soot, soot_new, soot_agd, soot_cld,                 &
           l_sulpc_so2, so2, so4_aitken, so4_accu, so4_diss,     &
           l_sulpc_nh3, nh3, l_sulpc_dms, dms, l_dust,           &
           dust_div1, dust_div2, dust_div3, dust_div4,           &
           dust_div5, dust_div6,                                 &
           l_biomass, bmass_new, bmass_agd, bmass_cld,           &
           l_ocff, ocff_new, ocff_agd, ocff_cld,                 &
           l_nitrate, nitr_acc, nitr_diss,                       &
           l_use_cariolle, ozone_tracer, tracer,                 &
           tracer_ukca, row_length, rows,                        &
           model_levels, tr_levels, tr_vars, tr_ukca )
    END IF  ! L_tracer and CycleNo == NumCycles

    ! set zero level star quantities for Endgame
    ! currently set equal to level 1
    IF (i_bl_vn /= i_bl_vn_0) THEN
       CALL set_star_zero_level(                                        &
             theta_star,                                                &
             q_star,                                                    &
             qcl_star,                                                  &
             qcf_star,                                                  &
             cf_star,                                                   &
             cfl_star,                                                  &
             cff_star,                                                  &
             qcf2_star,                                                 &
             qrain_star,                                                &
             qgraup_star,                                               &
             L_mcr_qgraup,                                              &
             L_mcr_qrain,                                               &
             L_mcr_qcf2)
    END IF

    IF (ltimer) CALL timer('AS Atmos_Phys2 (AP2)',6)
    IF (ltimer .OR. lstashdumptimer) CALL timer('AS STASH',5)

    ! DEPENDS ON: Atm_Step_stash
    CALL atm_step_stash( &
         errorstatus, 2)

    IF (ltimer .OR. lstashdumptimer) CALL timer('AS STASH',6)

    !  --------- UM Section 35---- Stochastic Physics SKEB2 -------------
 
    IF (ltimer) CALL timer('AS Stochastic_Phys',5)

    ! Call SKEB2 after the physics have been completed
    IF ( cycleno == 1 ) THEN

      IF (sf(0,35)) THEN
        ! ALLOCATE diagnostic space for STASH
        ALLOCATE (stashwork35(stash_maxlen(35,atmos_im)))
      ELSE
        ! Code to ALLOCATE unity arrays when not used
        ALLOCATE(stashwork35(1))
        stashwork35(1) = 0.0
      END IF
     
     ! Compute associated Legendre Polynomials for SKEB2 and SPT
     ! forcing pattern
      IF (first_atmstep_call) THEN
        IF (l_skeb2 .OR. l_spt) THEN
          CALL legendre_poly_comp_stph(global_rows,delta_phi,create_legendre)
        END IF
      END IF
      
      IF (l_retain_stph_tendencies) CALL init_stph_tendencies()
      
      IF (l_skeb2) THEN
        ! SKEB2 only called once in the outer loop

        IF (ltimer) CALL timer ('SKEB2',3)

        CALL stph_skeb2(                                                &
        ! in
                   row_length, rows, n_rows, model_levels,              &
                   delta_phi, delta_lambda,                             &
                   wetrho_r_sq_n, u, v, w,                              &
        ! inout
                   r_u_skeb, r_v_skeb,                                  &
        ! STASH array sf
                   stashwork35, first_atmstep_call)

        IF (ltimer) CALL timer ('SKEB2',4)

      END IF ! L_SKEB2

      ! Call SPT (inc computing the gaussian forcing pattern
      IF (l_spt) THEN

        CALL for_pattern(                                               &
        ! in
                 row_length, rows, model_levels, first_atmstep_call)

        ! Use dry row if perturbing physics in mixing ratios,
        ! or wetrho if specific humidities
        IF (l_mr_physics) THEN

          CALL SPT_main(                                                &
          !in
                  row_length, rows, n_rows, land_field,                 &
                  exner_theta_levels,p_theta_levels,dryrho,             &
                  first_atmstep_call,                                   &
                  theta_star, q_star, land,orog_sd,                     &
          ! STASH array sf
                  stashwork35)
        ELSE 
          CALL SPT_main(                                                &
          !in
                  row_length, rows, n_rows, land_field,                 &
                  exner_theta_levels,p_theta_levels,wetrho_r_sq_n,      &
                  first_atmstep_call,                                   &
                  theta_star, q_star, land,orog_sd,                     &
          ! STASH array sf
                  stashwork35)
          
        END IF 
      END IF ! (l_spt)

      ! Deallocate array with Legendre Polynomials
      IF (first_atmstep_call) THEN
        IF (l_skeb2 .OR. l_spt) THEN
          CALL legendre_poly_comp_stph(global_rows,delta_phi,destroy_legendre)
        END IF
      END IF

      ! Send diagnostics to STASH
      IF (sf(0,35) .AND. errorstatus == 0) THEN

        ! DEPENDS ON: stash
        CALL stash( atmos_sm, atmos_im, 35, stashwork35,                &
                  errorstatus, cmessage)

      END IF !(sf(0,35))

      DEALLOCATE (stashwork35)

    END IF ! end if over cycle_no

    ! Call Add SPT perturbations for each ENDGame cycle
    IF (l_spt) CALL SPT_add(theta_star,q_star)

    !------ End over stochastic physics loop

    ! PV TRACER:  Compute fast + stph dPV increments
    IF ((l_pv_tracer .OR. l_diab_tracer) .AND. cycleno == numcycles )   &
      CALL dyn_tr_fast( r_u_p2_n, r_v_p2_n, theta_star_n,               &
                        exner_theta_levels,                             &
                        ! PV-tracers
                        dPV_adv, dPV_conv, dPV_bl, dPV_stph,            &
                        ! theta tracers
                        dtheta_bl, dtheta_bl_mix, dtheta_bl_LH,         &
                        dtheta_conv)


    IF (ltimer) CALL timer('AS Stochastic_Phys',6)
    
    !-----------------------------------------------
    ! -------- AC ASSIMILATION SECTION -------------
    !-----------------------------------------------

    IF (LTimer) CALL timer('AS Assimilation',5)

    ! Only call LHN on final outer-loop
    IF (cycleno == numcycles) THEN

      ! DEPENDS ON: Atm_Step_ac_assim
      CALL Atm_Step_ac_assim( &
             obs_flag_len, obs_len, obs_flag, obs,                   &
             q_star, qcl_star, qcf_star, theta_star, ntml, cumulus,  &
              errorstatus)

    END IF
    !-----------------------------------------------

    IF (LTimer) CALL timer('AS Assimilation',6)

  ELSE

    CALL swap_bounds(r_u_p2,                                               &
                     udims_s%i_len - 2*udims_s%halo_i,                     &
                     udims_s%j_len - 2*udims_s%halo_j,                     &
                     udims_s%k_len,                                        &
                     udims_s%halo_i, udims_s%halo_j,                       &
                     fld_type_u,swap_field_is_vector)
    CALL swap_bounds(r_v_p2,                                               &
                     vdims_s%i_len - 2*vdims_s%halo_i,                     &
                     vdims_s%j_len - 2*vdims_s%halo_j,                     &
                     vdims_s%k_len,                                        &
                     vdims_s%halo_i, vdims_s%halo_j,                       &
                     fld_type_v,swap_field_is_vector)
    CALL swap_bounds(u,                                                    &
                     udims_s%i_len - 2*udims_s%halo_i,                     &
                     udims_s%j_len - 2*udims_s%halo_j,                     &
                     udims_s%k_len,                                        &
                     udims_s%halo_i, udims_s%halo_j,                       &
                     fld_type_u,swap_field_is_vector)
    CALL swap_bounds(v,                                                    &
                     vdims_s%i_len - 2*vdims_s%halo_i,                     &
                     vdims_s%j_len - 2*vdims_s%halo_j,                     &
                     vdims_s%k_len,                                        &
                     vdims_s%halo_i, vdims_s%halo_j,                       &
                     fld_type_v,swap_field_is_vector)
    CALL swap_bounds(exner,                                                &
                     pdims_s%i_len - 2*pdims_s%halo_i,                     &
                     pdims_s%j_len - 2*pdims_s%halo_j,                     &
                     pdims_s%k_len,                                        &
                     pdims_s%halo_i, pdims_s%halo_j,                       &
                     fld_type_p,swap_field_is_scalar)
    CALL swap_bounds(theta,                                                &
                     tdims_s%i_len - 2*tdims_s%halo_i,                     &
                     tdims_s%j_len - 2*tdims_s%halo_j,                     &
                     tdims_s%k_len,                                        &
                     tdims_s%halo_i, tdims_s%halo_j,                       &
                     fld_type_p,swap_field_is_scalar)
    CALL swap_bounds(wetrho_r_sq_n,                                        &
                     pdims_s%i_len - 2*pdims_s%halo_i,                     &
                     pdims_s%j_len - 2*pdims_s%halo_j,                     &
                     pdims_s%k_len,                                        &
                     pdims_s%halo_i, pdims_s%halo_j,                       &
                     fld_type_p,swap_field_is_scalar)
    CALL swap_bounds(q,                                                    &
                     tdims_l%i_len - 2*tdims_l%halo_i,                     &
                     tdims_l%j_len - 2*tdims_l%halo_j,                     &
                     tdims_l%k_len,                                        &
                     tdims_l%halo_i, tdims_l%halo_j,                       &
                     fld_type_p,swap_field_is_scalar)

    IF ( cycleno == 1 .AND. l_physics_store ) THEN
      ! RHcrit will be a 1D parametrized array input from user interface
      rhc_row_length = 1
      rhc_rows = 1
      ALLOCATE ( rhcpt(rhc_row_length, rhc_rows, model_levels) )
      DO k=1, model_levels
        rhcpt(1,1,k) = rhcrit(k)
      END DO
    END IF

    IF (l_idealised_data .AND. cycleno == numcycles) THEN
      ! Calculate idealised diagnostics
      CALL eg_idl_forcing2(                                                &
                 u, v,  exner_theta_levels, exner, theta, q, dryrho,       &
                 wetrho_r_sq_n, r_u_p2, r_v_p2, q_star, theta_star,        &
                 qlimit, errorstatus )
    END IF

  END IF lphysics

  ! After all calls to idealised forcing within a timestep output diagnostics
  IF (l_idealised_data .AND. cycleno == numcycles) THEN
    IF ( sf(0,53) ) THEN
      IF (ltimer .OR. lstashdumptimer) CALL timer('AS STASH',5)
    
      ! Allocate stashwork array for idealised diagnostics
      ALLOCATE ( stashwork53(stash_maxlen(53, atmos_im)) )

      ! Transfer diagnostics to stashwork array
      CALL diagnostics_ideal(row_length,rows,                             &
                             stashwork53)

      ! Pass diagnostics to stash
      ! DEPENDS ON: stash
      CALL stash( atmos_sm, atmos_im, 53, stashwork53,                         &
                  errorstatus, cmessage)

      ! Deallocate stashwork array
      DEALLOCATE ( stashwork53 )

      IF (ltimer .OR. lstashdumptimer) CALL timer('AS STASH',6)
    END IF
    ! Deallocate arrays used to hold idealised diagnostics for stash
    CALL dealloc_ideal_diag( )
  END IF

  IF (ltimer) CALL timer('AS Atmos_Phys2 (AP2)',5)

  ! Turn physics back on
  IF (.NOT. l_run_with_physics2 ) THEN
    l_physics=l_physics_store
  END IF

  IF (ltimer) CALL timer('AS CONVERT',5)

  ! Compute fast physics source terms for ENDGame
  ! ***********************************************************************
  !  CONVERSION FROM ND to EG within eg_R_S
  CALL eg_r_s(                                                            &
           theta_star,q_star,qcl_star,qcf_star,qcf2_star,                 &
           qrain_star, qgraup_star,                                       &
           theta_star_n,m_star,mcl_star,mcf_star,mcf2_star,               &
           mrain_star,mgraup_star, couple_app, l_skeb2)

  ! ***********************************************************************

  IF (ltimer) CALL timer('AS CONVERT',6)
  IF (ltimer) CALL timer('AS Atmos_Phys2 (AP2)',6)
  IF (ltimer) CALL timer('AS Diffusion',5)

  IF (errorstatus  ==  0 ) THEN  !it2

    ! -------- DIFFUSION SECTION ------------------

    IF ( L_diff_active .AND. L_subfilter_horiz ) THEN !  it3

      IF ( sf(0,13) ) THEN
        CALL diff_incs_init(s_u, s_v, s_w, s_thetav, s_m_v)
      END IF

      ! DEPENDS ON: eg_diff_ctl
      CALL eg_diff_ctl(                                                   &
                       L_backwards,                                       &
                       pos_timestep, neg_timestep,                        &
                       thetav_np1,m_v_np1,u_np1,v_np1,w_np1,              &
                       m_cl_np1, m_cf_np1, m_r_np1, m_cf2_np1, m_gr_np1,  &
                       r_theta_levels, r_rho_levels,                      &
                       offx, offy, halo_i, halo_j,                        &
                       row_length, rows, n_rows, model_levels,            &
                       first_constant_r_rho_level,                        &
                       s_thetav, s_m_v, s_u, s_v, s_w,                    &
                       s_m_cl, s_m_cf, s_m_r, s_m_cf2, s_m_gr)

      ! ----------------------------------------------------------------------
      ! Section 13.1 Diagnostics at from diffusion and divergence damping
      ! ----------------------------------------------------------------------
      ! Apply diagnostics only at last cycle
      IF ( cycleno == numcycles ) THEN

        ! section 13:
        IF ( sf(0,13) .AND. errorstatus == 0) THEN

          IF (ltimer .OR. lstashdumptimer) CALL timer('AS STASH',5)

          ! Allocate diagnostic space for STASH
          IF (.NOT. l_filter) ALLOCATE (stashwork13(stash_maxlen(13,atmos_im)))

          ! DEPENDS ON: diagnostics_dif
          CALL diagnostics_dif(                                        &
             row_length, rows, n_rows, model_levels, bl_levels,        &
          ! primary  fields:
                   theta, q,                                                 &
          ! wind field increments after  dif :
                   s_u, s_v, s_w,                                            &
          ! Current theta+dtheta and q+dq values
                   s_thetav, s_m_v,                                          &
                   exner_theta_levels,                                       &
                   stashwork13)

          CALL diff_incs_dealloc()

          ! DEPENDS ON: Atm_Step_stash
          CALL atm_step_stash( &
             errorstatus, 3)
          IF (ltimer .OR. lstashdumptimer) CALL timer('AS STASH',6)

        END IF !   SF(0,13)

      END IF ! CycleNo == NumCycles

    END IF      ! L_diff_active  ei3

  END IF       ! ErrorStatus  ==  0  ei2

  IF ( cycleno == numcycles ) THEN

    ! DEPENDS ON: Atm_Step_alloc_4A
    CALL Atm_Step_alloc_4A(                                         &
         cf_star, cfl_star, cff_star, q_star, qcl_star, qcf_star,   &
         frac_control, r_u, r_v, r_w, errorstatus, 'dealsmag' )

  END IF ! cycleno == numcycles

  IF (ltimer) CALL timer('AS Diffusion',6)

  !---------- END DIFFUSION SECTION -------------------------------

  !=== Polar filter + diffusion of fast physics section ==================================
  ! ----------------------------------------------------------------------
  ! Section 0.4  Filter winds and theta near poles if active
  !              Do horizontal diffusion as a filter if active
  ! ----------------------------------------------------------------------
  IF ( L_polar_filter_incs .OR. L_filter_incs ) THEN  !it4
    IF (ltimer) CALL timer ('AS Filter',5)

    CALL filter_diag_printing5()

    CALL eg_NI_filter_incs_Ctl(                                       &
                          r_theta, r_u, r_v, r_w, r_rho,              &
                          s_thetav, s_u, s_v, s_w,                    &
                          thetav, u, v, w, dryrho, exner,             &
                          row_length, rows, n_rows, model_levels,     &
                          r_theta_levels, r_rho_levels,               &
                          r_at_u, r_at_v,                             &
                          max_121_rows, u_sweeps, v_sweeps,           &
                          global_u_filter, global_v_filter,           &
                          u_begin, u_end, v_begin, v_end,             &
                          diff_coeff_phi, diff_coeff_u, diff_coeff_v, &
                          first_constant_r_rho_level,                 &
                          first_constant_r_rho_level_m1,              &
                          horizontal_level,                           &
                          offx, offy, halo_i, halo_j,                 &
                          nproc_y, at_extremity,                      &
                          L_polar_filter_incs, L_diff_incs,           &
                          L_pofil_hadgem2,                            &
                          l_eliminate_rho, cycleno,                   &
                          xi1_u, xi1_p, xi2_p,                        &
                          pole_consts, gc_proc_row_group,             &
                          global_row_length,                          &
                          csxi2_v, csxi2_p)

    CALL filter_diag_printing5()

    IF (ltimer) CALL timer ('AS Filter',6)
  END IF    ! L_polar_filter_incs .or. L_filter_incs  ei4
  !=== End Polar filter + diffusion of increments section ===================

  ! ----------------------------------------------------------------------
  ! Form and solve Helmholtz equation and update variables.
  ! ----------------------------------------------------------------------
  IF (ltimer) CALL timer ('AS Solver',5)

  IF ( cycleno == numcycles ) THEN
    IF (sf(0,10)) THEN
        CALL solv_incs_init(exner_theta_levels, theta_star,            &
                            r_u_p2, r_v_p2, r_w_p2, u, v, w)
    END IF
  END IF

  IF ( l_inc_solver ) THEN
    CALL eg_sl_helmholtz_inc(                                          &
        exner_np1, exner_surf_np1, u_np1, v_np1, w_np1,                &
        etadot_np1, dryrho_np1, thetav_np1, exner_prime_term,          &
        m_v_np1, m_cl_np1, m_cf_np1,m_r_np1, m_gr_np1, m_cf2_np1,      &
        cycleno, ih, innits, gcr_precon_option, gcr_use_residual_tol,  &
        GCR_Diagnostics, r_u_d, r_v_d, r_w_d,r_theta_d, r_rho_d,       &
        r_m_v_d, r_m_cl_d, r_m_cf_d,r_m_r_d, r_m_gr_d, r_m_cf2_d,      &
        solver_tolerance, tol_sc_fact,n_rows, row_length, rows,        &
        model_levels,s_u,s_v,s_w,s_thetav,s_m_v,s_m_cl,s_m_cf,s_m_cf2, &
        s_m_r,s_m_gr,psi_w_surf, psi_w_lid)
  ELSE
    CALL eg_sl_helmholtz (                                             &
        exner_np1, exner_surf_np1, u_np1, v_np1, w_np1,                &
        etadot_np1, dryrho_np1, thetav_np1, exner_prime_term,          &
        m_v_np1, m_cl_np1, m_cf_np1,m_r_np1, m_gr_np1, m_cf2_np1,      &
        cycleno, ih,  innits, gcr_precon_option, gcr_use_residual_tol, &
        GCR_Diagnostics, r_u_d, r_v_d, r_w_d,r_theta_d, r_rho_d,       &
        l_eliminate_rho, r_m_v_d, r_m_cl_d, r_m_cf_d,                  &
        r_m_r_d, r_m_gr_d, r_m_cf2_d,solver_tolerance, tol_sc_fact,    &
        n_rows,row_length, rows, model_levels,                         &
        s_u,s_v,s_w,s_thetav,s_m_v,s_m_cl,s_m_cf,s_m_cf2,              &
        s_m_r,s_m_gr,psi_w_surf, psi_w_lid)
  END IF

  ! ---------------------------------------------------------------
  ! Section 10 Diagnostics at end of helmholtz solver
  ! ---------------------------------------------------------------
  ! Apply diagnostics at final cycle only
  IF ( cycleno == numcycles ) THEN

    ! Section 10: 'dynamics solver' based quantities
    IF ( sf(0,10) .AND. ErrorStatus == 0 ) THEN

      CALL solv_incs_calc(exner_np1, thetav_np1, m_v_np1, u_np1, v_np1, w_np1)

      ! Allocate diagnostic space for STASH
      ALLOCATE (stashwork10(stash_maxlen(10,atmos_im)))

      CALL diagnostics_solver(                                          &
                    row_length, rows, n_rows,                           &
      ! diagnostic array
                    stashwork10)

      ! DEPENDS ON: stash
      CALL stash(atmos_sm,atmos_im,10,stashwork10,errorstatus,cmessage)

      DEALLOCATE (stashwork10)

      CALL solv_incs_dealloc()

    END IF !   SF(0,10)

  END IF ! CycleNo == NumCycles

  IF (ltimer) CALL timer ('AS Solver',6)

  IF (L_tracer .AND. CycleNo == NumCycles) THEN  !it5
    IF (ltimer) CALL timer('AS S-L Advect (AA)',5)

    IF (ltimer) CALL timer('AA SL_Tracer',5)

    !      Cariolle scheme is called to calculate the tracer ozone. All the tracers
    !      are calculated at the end of the timestep. The ozone tracer
    !      calculated here will be used in the radiation scheme on the next timestep.

    IF (l_use_cariolle) THEN

      IF (PrintStatus  >=  PrStatus_Normal .AND.                  &
       first_atmstep_call .AND. mype == 0) THEN
        CALL umPrint('Atm_Step_4A: Calling Cariolle_o3_psc', &
            src='atm_step_4A')
      END IF

      ! DEPENDS ON: cariolle_o3_psc
      CALL cariolle_o3_psc (ozone_tracer,                       &
             o3_prod_loss,   o3_p_l_vmr,                        &
             o3_vmr,         o3_p_l_temp,                       &
             o3_temp,        o3_p_l_colo3,                      &
             o3_colo3,                                          &
             theta,                                             &
             p_theta_levels,                                    &
             offx,offy,theta_off_size,                          &
             rows,row_length,                                   &
             exner_theta_levels,                                &
             model_levels)
      !
      ! Halos updated
      CALL swap_bounds(ozone_tracer,                            &
          pdims_s%i_len - 2*pdims_s%halo_i,                     &
          pdims_s%j_len - 2*pdims_s%halo_j,                     &
          pdims_s%k_len,                                        &
          pdims_s%halo_i, pdims_s%halo_j, fld_type_p,swap_field_is_scalar)

    ELSE
      IF (PrintStatus  >=  PrStatus_Normal .AND.                &
          first_atmstep_call) THEN

        WRITE(umMessage,*) 'Atm_Step_4A: Cariolle scheme not called'
        CALL umPrint(umMessage,src='atm_step_4A')
      END IF
    END IF

    IF (ltimer) CALL timer('AA SL_Tracer',6)
    IF (ltimer) CALL timer('AS S-L Advect (AA)',6)

  END IF  ! L_tracer ei5

  IF (ltimer) CALL timer('AS Diffusion',5)

  !=== Targetted diffusion of m_v_np1 ======================================
  IF ( l_tardiff_q ) THEN
    IF (printstatus >= prstatus_diag) THEN
      max_q_star     = MAXVAL(m_v_np1)
      min_q_star     = MINVAL(m_v_np1)
      CALL gc_rmax(1,nproc,ierr,max_q_star)
      CALL gc_rmin(1,nproc,ierr,min_q_star)
      IF ( mype == 0 ) THEN
        CALL umPrint('===================================================',&
            src='atm_step_4A')
        WRITE(umMessage,FMT='(A)') 'Calling tardiff for m_v_np1'
        CALL umPrint(umMessage,src='atm_step_4A')
        WRITE(umMessage,FMT='(A)') ' Max/Min before polar filter :'
        CALL umPrint(umMessage,src='atm_step_4A')
        WRITE(umMessage,FMT='(A,2E25.10)') 'm_v_np1 = ',             &
            MAX_q_star,MIN_q_star
        CALL umPrint(umMessage,src='atm_step_4A')
      END IF
    END IF

    ! DEPENDS ON: eg_tardiff_q_w
    CALL eg_tardiff_q_w(                                            &
                      m_v_np1, w_np1,                               &
                      r_theta_levels, r_rho_levels,                 &
                      offx, offy, halo_i, halo_j,                   &
                      rows, n_rows, row_length, model_levels,       &
                      w_conv_limit, tardiffq_factor, tardiffq_test, &
                      tardiffq_end,                                 &
                      csxi2_p,csxi2_v, nproc_y, L_pofil_hadgem2,    &
                      u_sweeps, u_begin, u_end, max_121_rows,       &
                      horizontal_level,                             &
                      first_constant_r_rho_level_m1)

    CALL swap_bounds(m_v_np1,                                       &
              tdims_s%i_len - 2*tdims_s%halo_i,                     &
              tdims_s%j_len - 2*tdims_s%halo_j,                     &
              tdims_s%k_len,                                        &
              tdims_s%halo_i, tdims_s%halo_j, fld_type_p,swap_field_is_scalar)

    IF (printstatus >= prstatus_diag) THEN
      max_q_star     = MAXVAL(m_v_np1)
      min_q_star     = MINVAL(m_v_np1)
      CALL gc_rmax(1,nproc,ierr,max_q_star)
      CALL gc_rmin(1,nproc,ierr,min_q_star)
      IF ( mype == 0 ) THEN
        CALL umPrint( ' ',src='atm_step_4A')
        CALL umPrint( ' Max/Min after polar incs filter :', &
            src='atm_step_4A')
        WRITE(umMessage,'(A,2E25.10)') 'm_v_np1 = ',MAX_q_star,MIN_q_star
        CALL umPrint(umMessage,src='atm_step_4A')
        CALL umPrint( &
            '===================================================' , &
            src='atm_step_4A')
      END IF
    END IF
  END IF    ! l_tardiff_q

  IF (ltimer) CALL timer('AS Diffusion',6)
  !=== End Targeted diffusion of m_v_np1 =====================================

  ! If a real data run reset lowest levels to reflect ND behaviour
  IF ( .NOT. L_prognostic_level0 .AND.                              &
       problem_number /= idealised_planet) THEN
!$OMP  PARALLEL DEFAULT(NONE)                                          &
!$OMP& SHARED( pdims_s, l_mcr_qrain, l_mcr_qgraup, l_mcr_qcf2,         &
!$OMP&         thetav_np1, m_v_np1, m_cl_np1, m_cf_np1, m_r_np1,       &
!$OMP&         m_gr_np1, m_cf2_np1 )                                   &
!$OMP& PRIVATE( i, j )
!$OMP DO SCHEDULE(STATIC)
    DO j = pdims_s%j_start, pdims_s%j_end
      DO i = pdims_s%i_start, pdims_s%i_end
        thetav_np1(i,j,0) = thetav_np1(i,j,1)
        m_v_np1(i,j,0)    = m_v_np1(i,j,1)
        m_cl_np1(i,j,0)   = m_cl_np1(i,j,1)
        m_cf_np1(i,j,0)   = m_cf_np1(i,j,1)
      END DO
    END DO
!$OMP END DO NOWAIT

    IF ( l_mcr_qrain ) THEN
!$OMP DO SCHEDULE(STATIC)
      DO j = pdims_s%j_start, pdims_s%j_end
        DO i = pdims_s%i_start, pdims_s%i_end
          m_r_np1(i,j,0)    = m_r_np1(i,j,1)
        END DO
      END DO
!$OMP END DO NOWAIT
    END IF

    IF ( l_mcr_qgraup ) THEN
!$OMP DO SCHEDULE(STATIC)
      DO j = pdims_s%j_start, pdims_s%j_end
        DO i = pdims_s%i_start, pdims_s%i_end
          m_gr_np1(i,j,0)   = m_gr_np1(i,j,1)
        END DO
      END DO
!$OMP END DO NOWAIT
    END IF

    IF ( l_mcr_qcf2 ) THEN
!$OMP DO SCHEDULE(STATIC)
      DO j = pdims_s%j_start, pdims_s%j_end
        DO i = pdims_s%i_start, pdims_s%i_end
          m_cf2_np1(i,j,0)  = m_cf2_np1(i,j,1)
        END DO
      END DO
!$OMP END DO
    END IF
!$OMP END PARALLEL

  END IF

  ! remove any potential instabillity:
  IF (l_eg_dry_static_adj) CALL eg_dry_static_adj( thetav_np1,      &
                                      dryrho_np1,exner_np1,intw_w2rho)

END DO  outerloop ! end iterations for trajectory calc (Outer loop)

! -------------------------------------------------
! Mass Conservation and Corrections
! -------------------------------------------------

IF (sf(0,12)) THEN
  CALL adv_correct_incs_init(exner_np1, thetav_np1,                     &
                             m_v_np1, m_cl_np1, m_cf_np1,               &
                             m_r_np1, m_gr_np1, m_cf2_np1)
END IF

! Pseudo lateral boundary flux for rho if necessary
IF (model_type == mt_lam .AND. conserve_dry_mass /= not_conserved) THEN

  IF (ltimer) CALL timer('AS S-L Advect (AA)',5)
  IF (ltimer) CALL timer('AA SL_RHO_LAM_CONSERV',5)

  CALL eg_rho_pseudo_lbflux(g_i_pe,                                     &
                depart_scheme, depart_order,                            &
                high_order_scheme(rho_sl), monotone_scheme(rho_sl),     &
                ritchie_high_order_scheme, ritchie_monotone_scheme,     &
                first_constant_r_rho_level,interp_vertical_search_tol,  &
                check_bottom_levels, l_shallow,                         &
                l_rk_dps,                                               &
                l_high(rho_sl),l_mono(rho_sl), l_ritchie_high,          &
                l_ritchie_mono, lam_max_cfl,                            &
                etadot, u_np1, v_np1, w_np1,etadot_np1, u_adv, v_adv,   &
                w_adv, r_rho, r_rho_d, errorstatus,                     &
                pseudo_lbflux_mass )

  IF (ltimer) CALL timer('AA SL_RHO_LAM_CONSERV',6)
  IF (ltimer) CALL timer('AS S-L Advect (AA)',6)

ELSE 
  pseudo_lbflux_mass = 0.0
END IF

! Compute total mass at the end of the timestep
CALL eg_total_mass_region(is, ie, js, je, l_exclude_rim=.TRUE.)
total_rho = eg_total_mass(dryrho_np1, l_exclude_rim=.TRUE.)

! Enforce conservation of total mass of dry air
SELECT CASE(conserve_dry_mass)

CASE(not_conserved)
! Print mass-loss error

  IF (PrintStatus > PrStatus_Oper .AND. mype ==0) THEN
    WRITE(umMessage,FMT='(A,E19.12)') 'Dry mass error = ',                     &
                                      (total_rho - total_rho_init) /           &
                                      total_rho_init
    CALL umPrint(umMessage,src='atm_step_4A')
  END IF

CASE(constant_factor)

  IF (model_type == mt_lam) THEN
    mass_fix_factor = (total_rho_n + pseudo_lbflux_mass) / total_rho
  ELSE
    mass_fix_factor = total_rho_init / total_rho
  END IF

  WRITE(umMessage,FMT='(A,E19.12)') 'Mass correction factor = ',               &
                                     mass_fix_factor
  CALL umPrint(umMessage,src='atm_step_4A')
   
  DO k = pdims%k_start, pdims%k_end
    DO j = pdims_s%j_start, pdims_s%j_end
      DO i = pdims_s%i_start, pdims_s%i_end
        dryrho_np1(i,j,k) = dryrho_np1(i,j,k) * mass_fix_factor
      END DO
    END DO
  END DO

CASE(linear_factor, linear_factor_IE)

  total_rho_r = eg_total_mass_r(dryrho_np1)
  total_gr    = eg_total_gr  (dryrho_np1)
  total_gr_r  = eg_total_gr_r(dryrho_np1)

  mass_fix_A  = (total_rho_init*total_gr_r-total_rho_r*total_gr)       &
               /(total_rho*total_gr_r-total_rho_r*total_gr)

  mass_fix_B  = ((total_rho  -total_rho_init)*total_gr)                &
               /(total_rho*total_gr_r-total_rho_r*total_gr)

!$OMP  PARALLEL DEFAULT(NONE)                                          &
!$OMP& SHARED( pdims, pdims_s, mass_fix_A, mass_fix_B, r_rho_levels,   &
!$OMP&         dryrho_np1, r_theta_levels, conserve_dry_mass,          &
!$OMP&         l_prognostic_level0, thetav_np1 )                       &
!$OMP& PRIVATE( i, j, k, mass_fix_factor )
!$OMP DO SCHEDULE(STATIC)
  DO k = pdims%k_start, pdims%k_end
    DO j = pdims_s%j_start, pdims_s%j_end
      DO i = pdims_s%i_start, pdims_s%i_end
        mass_fix_factor   = mass_fix_A + mass_fix_B*r_rho_levels(i,j,k)
        dryrho_np1(i,j,k) = dryrho_np1(i,j,k)*mass_fix_factor
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

! Adjust theta_vd to (approximately) preserve total internal energy
  IF (conserve_dry_mass == linear_factor_IE) THEN
!$OMP DO SCHEDULE(STATIC)
    DO k = pdims%k_start, pdims%k_end
      DO j = pdims_s%j_start, pdims_s%j_end
        DO i = pdims_s%i_start, pdims_s%i_end
          mass_fix_factor   = mass_fix_A + mass_fix_B*r_theta_levels(i,j,k)
          thetav_np1(i,j,k) = thetav_np1(i,j,k) / mass_fix_factor
        END DO
      END DO
    END DO
!$OMP END DO
    IF (.NOT. l_prognostic_level0) THEN
!$OMP DO SCHEDULE(STATIC)
      DO j = pdims_s%j_start, pdims_s%j_end
        DO i = pdims_s%i_start, pdims_s%i_end
          thetav_np1(i,j,0) = thetav_np1(i,j,1)
        END DO
      END DO
!$OMP END DO
    END IF
  END IF
!$OMP END PARALLEL

  IF (PrintStatus > PrStatus_Oper .AND. mype == 0) THEN
    mass_fix_factor = total_rho_init / total_rho
    WRITE(umMessage,FMT='(A,E19.12)') 'Mass Loss corrected for ',              &
                                        mass_fix_factor
    CALL umPrint(umMessage,src='atm_step_4A')
  END IF

END SELECT

! Apply OCF conservation algorithm to thetav
IF (L_priestley_correct_thetav) THEN
  IF (model_type == mt_lam) THEN
    ErrorStatus = 1
    Cmessage = 'l_priestley_correct_thetav = T incompatible with LAM.'
    CALL Ereport(RoutineName, ErrorStatus, Cmessage)
  END IF
  CALL eg_correct_thetav_priestley(g_i_pe, dryrho, dryrho_np1,        &
                                   depart_xi1_w, depart_xi2_w,        &
                                   depart_xi3_w, r_theta, s_thetav,   &
                                   thetav_np1)
END IF

pseudo_lbflux_moisture = 0.0
! Pseudo lateral boundary flux for q or tracers if necessary
IF (model_type == mt_lam) THEN
  IF ( ((.NOT. l_dry) .AND. l_conserv(moist_sl)) ) THEN

    IF (ltimer) CALL timer('AS S-L Advect (AA)',5)
    IF (ltimer) CALL timer('AA SL_UPDATE_DPTS_LAM_CONSERV',5)

    ! Update w departure points for tracers.
    ! Perhaps not necessary, but for consistency as
    !   depart_xi[123]_rho is updated in subroutine eg_sl_rho which is
    !   called from eg_rho_pseudo_lbflux.
    !
    CALL eg_sl_full_wind(g_i_pe,depart_scheme,                          &
               depart_order, high_order_scheme(wind_sl),                &
               monotone_scheme(wind_sl),                                &
               ritchie_high_order_scheme, ritchie_monotone_scheme,      &
               first_constant_r_rho_level,                              &
               interp_vertical_search_tol, check_bottom_levels,         &
               l_shallow,                                               &
               l_rk_dps,                                                &
               l_high(wind_sl),  l_mono(wind_sl),  l_ritchie_high,      &
               l_ritchie_mono, lam_max_cfl,                             &
               etadot, u_np1, v_np1, w_np1, etadot_np1,                 &
               u_adv, v_adv, w_adv, r_u, r_v, r_w,                      &
               r_u_d, r_v_d, r_w_d,                                     &
               errorstatus )

    IF (ltimer) CALL timer('AA SL_UPDATE_DPTS_LAM_CONSERV',6)
    IF (ltimer) CALL timer('AA SL_MOISTURE_LAM_CONSERV',5)

    CALL eg_moisture_pseudo_lbflux( moisture_array_size,                &
                row_length, rows, n_rows, model_levels, halo_i,         &
                halo_j, offx, offy, datastart, g_i_pe,                  &
                high_order_scheme(moist_sl),monotone_scheme(moist_sl),  &
                l_high(moist_sl),l_mono(moist_sl),L_mcr_qrain,          &
                l_mcr_qcf2, l_mcr_qgraup,                               &
                r_m_v, r_m_cl, r_m_cf,r_m_r, r_m_gr, r_m_cf2,           &
                cf_bulk, cf_liquid, cf_frozen,exner_theta_levels,       &
                cf_star, cfl_star, cff_star,                            &
                dryrho, dryrho_np1, l_conserv(moist_sl),                &
                pseudo_lbflux_moisture, number_qs,                      &
                errorstatus)

    IF (ltimer) CALL timer('AA SL_MOISTURE_LAM_CONSERV',6)
    IF (ltimer) CALL timer('AS S-L Advect (AA)',6)

  END IF !( ((.NOT. l_dry) .AND. l_conserv(moist_sl)) )

END IF ! (model_type == mt_lam)
!
! fix moisture mixing ratios to satisfy mass conservation if
! L_conserv(moist_SL) = true; if not then the
! eg_correct_moisture routine will simply clip the fields that
! are below some minimum values [ m >= m_min ]

IF (ltimer) CALL timer('AS S-L Advect (AA)',5)
IF (ltimer) CALL timer('AS CORRECT_MOISTURE',5)

IF ( L_priestley_correct_moist ) THEN

  ! Apply Priestley Scheme to conserve moisture fields
   CALL eg_correct_moisture_priestley(dryrho, dryrho_np1,          &
        r_m_v, r_m_cl, r_m_cf, r_m_r, r_m_gr, r_m_cf2,             &
        m_v_np1, m_cl_np1, m_cf_np1, m_r_np1, m_gr_np1, m_cf2_np1, &
        s_m_v, s_m_cl, s_m_cf, s_m_r, s_m_gr, s_m_cf2,             &
        qlimit, L_mcr_qrain, L_mcr_qgraup, L_mcr_qcf2,             &
        L_conserv(moist_SL),                                       &
        g_i_pe, depart_xi1_w, depart_xi2_w,  depart_xi3_w          )

ELSE     ! Default Correction scheme

  IF ( l_fix_conserv ) THEN
     CALL eg_correct_moisture_fix(dryrho, dryrho_np1,                &
          r_m_v, r_m_cl, r_m_cf, r_m_r, r_m_gr, r_m_cf2,             &
          m_v_np1, m_cl_np1, m_cf_np1, m_r_np1, m_gr_np1, m_cf2_np1, &
          s_m_v, s_m_cl, s_m_cf, s_m_r, s_m_gr, s_m_cf2,             &
          qlimit, L_mcr_qrain, L_mcr_qgraup, L_mcr_qcf2,             &
          L_conserv(moist_SL),  mype, pseudo_lbflux_moisture,        &
          number_qs, L_conserv_smooth_lap)
  ELSE
     CALL eg_correct_moisture(dryrho, dryrho_np1,                    &
          r_m_v, r_m_cl, r_m_cf, r_m_r, r_m_gr, r_m_cf2,             &
          m_v_np1, m_cl_np1, m_cf_np1, m_r_np1, m_gr_np1, m_cf2_np1, &
          s_m_v, s_m_cl, s_m_cf, s_m_r, s_m_gr, s_m_cf2,             &
          qlimit, L_mcr_qrain, L_mcr_qgraup, L_mcr_qcf2,             &
          L_conserv(moist_SL),  mype, pseudo_lbflux_moisture,        &
          number_qs, L_conserv_smooth_lap)
  END IF

END IF         ! Use Priestley scheme ?

IF (sf(0,12)) THEN
  CALL adv_correct_incs_calc(exner_np1, thetav_np1,                     &
                             m_v_np1, m_cl_np1, m_cf_np1,               &
                             m_r_np1, m_gr_np1, m_cf2_np1)
  CALL diagnostics_adv_correct(row_length, rows, stashwork12)
  CALL adv_correct_incs_dealloc()

  ! DEPENDS ON: stash
  CALL stash(atmos_sm,atmos_im,12,stashwork12,errorstatus,cmessage)
  DEALLOCATE (stashwork12)
END IF

IF (ltimer) CALL timer('AS CORRECT_MOISTURE',6)
IF (ltimer) CALL timer('AS S-L Advect (AA)',6)
IF (ltimer) CALL timer('AS CONVERT',5)

! ***********************************************************************
!
!  CONVERSION FROM EG prognostic to ND prognostics! (start)

! compute the new wet rho * r**2
CALL update_wetrho_r_sq(time_level=tl_np1)

!  CONVERSION FROM EG prognostic to ND prognostics! (end)
!
! ***********************************************************************

IF (ltimer) CALL timer('AS CONVERT',6)

!$OMP  PARALLEL DEFAULT(NONE)                                              &
!$OMP& SHARED( sf, pdims, pdims_s, eot_inc_rho, wetrho_r_sq_np1,           &
!$OMP&         wetrho_r_sq_n )                                             &
!$OMP& PRIVATE( i, j, k )
IF ( sf(188,30) ) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k = pdims%k_start, pdims%k_end
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        eot_inc_rho(i,j,k) = wetrho_r_sq_np1(i,j,k) - wetrho_r_sq_n(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

!  ENSURE THAT DENSITY AT TIME LEVEL N IS CONSISTENT AT NEXT TIMESTEP
!  NOTE that we now have a strictly wrong naming of variables, as a variable
!  at time level n+1 is labeled time level n. However, due to systems
!  constraints this was deemed the lesser evil
!$OMP DO SCHEDULE(STATIC)
DO k = pdims_s%k_start, pdims_s%k_end
  DO j = pdims_s%j_start, pdims_s%j_end
    DO i = pdims_s%i_start, pdims_s%i_end
      wetrho_r_sq_n(i,j,k) = wetrho_r_sq_np1(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

IF ( L_tracer ) THEN
  !
  !  fix tracers to satisfy mass conservation if
  !  L_conserve_tracers = true; if not then the
  !  eg_correct_tracers routine will simply clip the fields that
  !  are below some minimum values (t >= t_min (generally zero) )

  IF (ltimer) CALL timer('AS S-L Advect (AA)',5)
  IF (ltimer) CALL timer('AS CORRECT_TRACERS',5)

  IF ( L_priestley_correct_tracers ) THEN

    ! Apply Priestley Scheme to conserve (Classic, murk, idealised) tracers
    CALL eg_correct_tracers_priestley(                            &
                          super_array_size,                       &
                          super_tracer_phys1, super_tracer_phys2, &
                          dryrho, dryrho_np1,                     &
                          co2, L_CO2_interactive,                 &
                          murk, L_murk_advect,                    &
                          soot_new, soot_agd, soot_cld, L_soot,   &
                          bmass_new, bmass_agd, bmass_cld,        &
                          L_biomass,                              &
                          ocff_new, ocff_agd, ocff_cld, l_ocff,   &
                          dust_div1,dust_div2,dust_div3,          &
                          dust_div4,dust_div5,dust_div6,          &
                          l_dust,                                 &
                          so2, so4_aitken, so4_accu,              &
                          so4_diss, nh3, dms,                     &
                          L_sulpc_so2, L_sulpc_nh3, l_sulpc_dms,  &
                          nitr_acc, nitr_diss, L_nitrate,         &
                          l_use_cariolle, ozone_tracer,           &
                          tracer, tr_vars,                        &
                          tr_ukca, tracer_ukca,                   &
                          L_conserve_tracers,                     &
                          g_i_pe, depart_xi1_w,                   &
                          depart_xi2_w,  depart_xi3_w             )
  ELSE

    ! Default conservation scheme

    IF ( l_fix_conserv ) THEN
      CALL eg_correct_tracers_fix(                                &
                          mype, super_array_size,                 &
                          super_tracer_phys1, super_tracer_phys2, &
                          dryrho, dryrho_np1,                     &
                          co2, L_CO2_interactive,                 &
                          murk, L_murk_advect,                    &
                          soot_new, soot_agd, soot_cld, L_soot,   &
                          bmass_new, bmass_agd, bmass_cld,        &
                          L_biomass,                              &
                          ocff_new, ocff_agd, ocff_cld, l_ocff,   &
                          dust_div1,dust_div2,dust_div3,          &
                          dust_div4,dust_div5,dust_div6,          &
                          l_dust,                                 &
                          so2, so4_aitken, so4_accu,              &
                          so4_diss, nh3, dms,                     &
                          L_sulpc_so2, L_sulpc_nh3, l_sulpc_dms,  &
                          nitr_acc, nitr_diss, L_nitrate,         &
                          l_use_cariolle, ozone_tracer,           &
                          tracer, tr_vars,                        &
                          tr_ukca, tracer_ukca,                   &
                          L_conserve_tracers,                     &
                          L_conserv_smooth_lap                    )
    ELSE
      CALL eg_correct_tracers(                                    &
                          mype, super_array_size,                 &
                          super_tracer_phys1, super_tracer_phys2, &
                          dryrho, dryrho_np1,                     &
                          co2, L_CO2_interactive,                 &
                          murk, L_murk_advect,                    &
                          soot_new, soot_agd, soot_cld, L_soot,   &
                          bmass_new, bmass_agd, bmass_cld,        &
                          L_biomass,                              &
                          ocff_new, ocff_agd, ocff_cld, l_ocff,   &
                          dust_div1,dust_div2,dust_div3,          &
                          dust_div4,dust_div5,dust_div6,          &
                          l_dust,                                 &
                          so2, so4_aitken, so4_accu,              &
                          so4_diss, nh3, dms,                     &
                          L_sulpc_so2, L_sulpc_nh3, l_sulpc_dms,  &
                          nitr_acc, nitr_diss, L_nitrate,         &
                          l_use_cariolle, ozone_tracer,           &
                          tracer, tr_vars,                        &
                          tr_ukca, tracer_ukca,                   &
                          L_conserve_tracers,                     &
                          L_conserv_smooth_lap                    )
    END IF
  END IF             ! Apply Priestley scheme to tracers ?

  ! Apply separate (Priestley) conservation to UKCA tracers if requested
  ! and not already done in "eg_correct_tracers_priestley".
  ! The eg_correct_tracers_ukca() routine should not be called with
  ! tr_ukca = 0 as this will result in out of bounds errors when
  ! slicing the super_tracer_phys1 array.
  IF ( (.NOT. l_conserve_ukca_with_tr) .AND. (tr_ukca > 0) ) THEN

    CALL eg_correct_tracers_ukca(                              &
                    row_length, rows, model_levels,            &
                    halo_i, halo_j, offx, offy,                &
                    datastart, g_i_pe,                         &
                    super_array_size,                          &
                    super_tracer_phys1, super_tracer_phys2,    &
                    tracer_ukca, dryrho, dryrho_np1,           &
                    depart_xi1_w, depart_xi2_w,  depart_xi3_w, &
                    tr_ukca                                    )
  END IF

  IF (ltimer) CALL timer('AS CORRECT_TRACERS',6)
  IF (ltimer) CALL timer('AS S-L Advect (AA)',6)

END IF  ! L_tracer

IF ( errorstatus == 0 ) THEN

  ! DEPENDS ON: Atm_Step_alloc_4A
  CALL atm_step_alloc_4A( &
       cf_star, cfl_star, cff_star, q_star, qcl_star, qcf_star, &
       frac_control, r_u, r_v, r_w, errorstatus,'lbc_updt')

END IF        !  ErrorStatus  ==  0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE ( i, j, k )
!$OMP DO SCHEDULE(STATIC)
DO j = tdims_s%j_start, tdims_s%j_end
  DO i = tdims_s%i_start, tdims_s%i_end
    thetav(i,j,0) = thetav_np1(i,j,0)
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO j = wdims_s%j_start, wdims_s%j_end
  DO i = wdims_s%i_start, wdims_s%i_end
    w(i,j,0)   = w_np1(i,j,0)
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels
  DO j = udims_s%j_start, udims_s%j_end
    DO i = udims_s%i_start, udims_s%i_end
      u(i,j,k)      = u_np1(i,j,k)
    END DO
  END DO
  DO j = vdims_s%j_start, vdims_s%j_end
    DO i = vdims_s%i_start, vdims_s%i_end
      v(i,j,k)      = v_np1(i,j,k)
    END DO
  END DO
  DO j = wdims_s%j_start, wdims_s%j_end
    DO i = wdims_s%i_start, wdims_s%i_end
      w(i,j,k)      = w_np1(i,j,k)
      etadot(i,j,k) = etadot_np1(i,j,k)
    END DO
  END DO
  DO j = pdims_s%j_start, pdims_s%j_end
    DO i = pdims_s%i_start, pdims_s%i_end
      dryrho(i,j,k) = dryrho_np1(i,j,k)
      exner(i,j,k)  = exner_np1(i,j,k)
    END DO
  END DO
  DO j = tdims_s%j_start, tdims_s%j_end
    DO i = tdims_s%i_start, tdims_s%i_end
      thetav(i,j,k) = thetav_np1(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

k = model_levels+1
!$OMP DO SCHEDULE(STATIC)
DO j = pdims_s%j_start, pdims_s%j_end
  DO i = pdims_s%i_start, pdims_s%i_end
    exner(i,j,k)  = exner_np1(i,j,k)
    exner_surf(i,j)   = exner_surf_np1(i,j)
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO  k = 0, model_levels
  DO j = tdims_s%j_start, tdims_s%j_end
    DO i = tdims_s%i_start, tdims_s%i_end
      m_v(i,j,k)    = m_v_np1(i,j,k)
      m_cl(i,j,k)   = m_cl_np1(i,j,k)
      m_cf(i,j,k)   = m_cf_np1(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

IF ( l_mcr_qrain ) THEN
!$OMP DO SCHEDULE(STATIC)
  DO  k = 0, model_levels
    DO j = tdims_s%j_start, tdims_s%j_end
      DO i = tdims_s%i_start, tdims_s%i_end
        m_r(i,j,k)    = m_r_np1(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF ( l_mcr_qgraup ) THEN
!$OMP DO SCHEDULE(STATIC)
  DO  k = 0, model_levels
    DO j = tdims_s%j_start, tdims_s%j_end
      DO i = tdims_s%i_start, tdims_s%i_end
        m_gr(i,j,k)   = m_gr_np1(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF ( l_mcr_qcf2 ) THEN
!$OMP DO SCHEDULE(STATIC)
  DO  k = 0, model_levels
    DO j = tdims_s%j_start, tdims_s%j_end
      DO i = tdims_s%i_start, tdims_s%i_end
        m_cf2(i,j,k)  = m_cf2_np1(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (i_cld_vn == i_cld_pc2) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k = tdims%k_start, tdims%k_end
    DO j = 1, rows
      DO i = 1, row_length
        cf_bulk(i,j,k)=cf_star(i,j,k)
        cf_liquid(i,j,k)=cfl_star(i,j,k)
        cf_frozen(i,j,k)=cff_star(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO
END IF
!$OMP END PARALLEL

! ---------------------------------------------------------------------
!   Section 6.1  Update lbcs for LAMs
! ---------------------------------------------------------------------
IF (model_type == mt_lam) THEN
  IF (ltimer) CALL timer('AS LAM_LBCS',5)

  IF (L_Fixed_lbcs) THEN
    L_update_lbcs = .FALSE.
    L_apply_lbcs = .TRUE.
  ELSE IF ( rim_stepsa == 0 ) THEN
    L_update_lbcs = .FALSE.
  ELSE
    L_update_lbcs = .TRUE.
  END IF !  L_Fixed_lbcs

  IF ( L_update_lbcs ) THEN

    ! DEPENDS ON: BOUNDVAL
    CALL boundval(lenrima(1,1,rima_type_norm),                 &
                  L_mcr_qcf2_lbc, L_mcr_qrain_lbc,             &
                  L_mcr_qgraup_lbc, L_pc2_lbc,                 &
                  L_murk_lbc, L_int_uvw_lbc,                   &
                  L_dust_div1_lbc,L_dust_div2_lbc,             &
                  L_dust_div3_lbc,L_dust_div4_lbc,             &
                  L_dust_div5_lbc,L_dust_div6_lbc,             &
                  L_so2_lbc,L_dms_lbc,L_so4_aitken_lbc,        &
                  L_so4_accu_lbc,L_so4_diss_lbc,               &
                  L_nh3_lbc,L_soot_new_lbc,L_soot_agd_lbc,     &
                  L_soot_cld_lbc,L_bmass_new_lbc,              &
                  L_bmass_agd_lbc,L_bmass_cld_lbc,             &
                  L_ocff_new_lbc,                              &
                  L_ocff_agd_lbc,L_ocff_cld_lbc,               &
                  L_nitr_acc_lbc, L_nitr_diss_lbc,             &
                  u_lbc, u_lbc_tend,                           &
                  v_lbc, v_lbc_tend,                           &
                  w_lbc, w_lbc_tend,                           &
                  rho_lbc, rho_lbc_tend,                       &
                  theta_lbc, theta_lbc_tend,                   &
                  q_lbc, q_lbc_tend,                           &
                  qcl_lbc, qcl_lbc_tend,                       &
                  qcf_lbc, qcf_lbc_tend,                       &
                  qcf2_lbc, qcf2_lbc_tend,                     &
                  qrain_lbc, qrain_lbc_tend,                   &
                  qgraup_lbc, qgraup_lbc_tend,                 &
                  cf_bulk_lbc, cf_bulk_lbc_tend,               &
                  cf_liquid_lbc, cf_liquid_lbc_tend,           &
                  cf_frozen_lbc, cf_frozen_lbc_tend,           &
                  exner_lbc, exner_lbc_tend,                   &
                  u_adv_lbc, u_adv_lbc_tend,                   &
                  v_adv_lbc, v_adv_lbc_tend,                   &
                  w_adv_lbc, w_adv_lbc_tend,                   &
                  murk_lbc, murk_lbc_tend,                     &
                  dust_div1_lbc, dust_div1_lbc_tend,           &
                  dust_div2_lbc, dust_div2_lbc_tend,           &
                  dust_div3_lbc, dust_div3_lbc_tend,           &
                  dust_div4_lbc, dust_div4_lbc_tend,           &
                  dust_div5_lbc, dust_div5_lbc_tend,           &
                  dust_div6_lbc, dust_div6_lbc_tend,           &
                  so2_lbc, so2_lbc_tend,                       &
                  dms_lbc, dms_lbc_tend,                       &
                  so4_aitken_lbc, so4_aitken_lbc_tend,         &
                  so4_accu_lbc, so4_accu_lbc_tend,             &
                  so4_diss_lbc, so4_diss_lbc_tend,             &
                  nh3_lbc, nh3_lbc_tend,                       &
                  soot_new_lbc, soot_new_lbc_tend,             &
                  soot_agd_lbc, soot_agd_lbc_tend,             &
                  soot_cld_lbc, soot_cld_lbc_tend,             &
                  bmass_new_lbc, bmass_new_lbc_tend,           &
                  bmass_agd_lbc, bmass_agd_lbc_tend,           &
                  bmass_cld_lbc, bmass_cld_lbc_tend,           &
                  ocff_new_lbc, ocff_new_lbc_tend,             &
                  ocff_agd_lbc, ocff_agd_lbc_tend,             &
                  ocff_cld_lbc, ocff_cld_lbc_tend,             &
                  nitr_acc_lbc, nitr_acc_lbc_tend,             &
                  nitr_diss_lbc, nitr_diss_lbc_tend,           &
                  tracer_lbc, tracer_lbc_tend,                 &
                  tracer_ukca_lbc, tracer_ukca_lbc_tend,       &
                  1, 0, ErrorStatus, cmessage)

  END IF ! L_update_lbcs

  !--------------------------------------------------------------
  !           Update primary fields with LAM LBC data
  !--------------------------------------------------------------
    ! DEPENDS ON: update_lam_lbcs
  CALL update_lam_lbcs(                                                &
         r_rho_levels, r_theta_levels,                                 &
         row_length,rows,n_rows,                                       &
         tr_vars,tr_lbc_vars,tr_levels,                                &
         a_max_trvars,A_tr_active_lbc_index,                           &
         tr_ukca,tr_lbc_ukca,                                          &
         offx,offy,halo_i,halo_j,at_extremity,                         &
         L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,                        &
         L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc, L_pc2_lbc, &
         L_murk, L_murk_lbc,                                           &
         L_LBC_balance, L_int_uvw_lbc,                                 &
         L_dust_div1, L_dust_div1_lbc, L_dust_div2, L_dust_div2_lbc,   &
         L_dust_div3, L_dust_div3_lbc, L_dust_div4, L_dust_div4_lbc,   &
         L_dust_div5, L_dust_div5_lbc, L_dust_div6, L_dust_div6_lbc,   &
         L_so2,L_so2_lbc,L_dms,L_dms_lbc,L_so4_aitken,L_so4_aitken_lbc,&
         L_so4_accu, L_so4_accu_lbc, L_so4_diss, L_so4_diss_lbc,       &
         L_nh3, L_nh3_lbc,                                             &
          L_soot_new_lbc,  L_soot_agd_lbc,       &
          L_soot_cld_lbc, L_bmass_new, L_bmass_new_lbc,     &
         L_bmass_agd, L_bmass_agd_lbc, L_bmass_cld, L_bmass_cld_lbc,   &
          L_ocff_new_lbc,                                   &
          L_ocff_agd_lbc,  L_ocff_cld_lbc,       &
         L_nitr_acc, L_nitr_acc_lbc, L_nitr_diss, L_nitr_diss_lbc,     &
         rimwidtha(rima_type_norm),rimweightsa,                        &
         lenrima(1,1,rima_type_norm),                                  &
         lbc_sizea(1,1,1,rima_type_norm),                              &
         lbc_starta(1,1,1,rima_type_norm),                             &
         theta_lbc, q_lbc, qcl_lbc,                                    &
         qcf_lbc, qcf2_lbc, qrain_lbc,                                 &
         qgraup_lbc, cf_bulk_lbc, cf_liquid_lbc,                       &
         cf_frozen_lbc, rho_lbc,exner_lbc,                             &
         u_lbc, v_lbc, w_lbc,                                          &
         u_adv_lbc, v_adv_lbc, w_adv_lbc,                              &
         murk_lbc,                                                     &
         dust_div1_lbc, dust_div2_lbc, dust_div3_lbc,                  &
         dust_div4_lbc, dust_div5_lbc, dust_div6_lbc,                  &
         so2_lbc, dms_lbc, so4_aitken_lbc,                             &
         so4_accu_lbc,so4_diss_lbc,nh3_lbc,                            &
         soot_new_lbc, soot_agd_lbc, soot_cld_lbc,                     &
         bmass_new_lbc, bmass_agd_lbc, bmass_cld_lbc,                  &
         ocff_new_lbc, ocff_agd_lbc, ocff_cld_lbc,                     &
         nitr_acc_lbc, nitr_diss_lbc,                                  &
         tracer_lbc, tracer_ukca_lbc,                                  &
         thetav, m_v, m_cl, m_cf,                                      &
         m_cf2, m_r, m_gr,                                             &
         cf_bulk, cf_liquid, cf_frozen,                                &
         dryrho, exner,                                                &
         u, v, w, u_adv, v_adv, w_adv, murk,                           &
         dust_div1, dust_div2, dust_div3,                              &
         dust_div4, dust_div5, dust_div6,                              &
         so2, dms, so4_aitken,so4_accu,                                &
         so4_diss, nh3,                                                &
         soot_new, soot_agd, soot_cld,                                 &
         bmass_new, bmass_agd, bmass_cld,                              &
         ocff_new, ocff_agd, ocff_cld,                                 &
         nitr_acc, nitr_diss,                                          &
         delta_phi, delta_lambda,                                      &
         base_phi, base_lambda,                                        &
         datastart, lat_rot_NP,                                        &
         global_row_length, global_rows,                               &
         tracer, tracer_ukca  )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEMP. FIX !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !-----------------------------------------------------------------------!
  ! Re-initialise etadot and exner_star: fixes needed for lateral bcs     !
  !-----------------------------------------------------------------------!

  CALL eg_calc_p_star(model_levels, row_length, rows, exner,           &
                      thetav, m_v, m_cl, m_cf, m_r, m_gr, m_cf2,       &
                      g_theta, exner_surf, psi_w_surf, psi_w_lid)
  CALL init_etadot()

  IF (ltimer) CALL timer('AS LAM_LBCS',6)
END IF     !   model_type  ==  mt_lam
! ---------------------------------------------------------------------
!   End - Section 6.1  Update lbcs for LAMs
! ---------------------------------------------------------------------

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                      &
!$OMP& SHARED( pdims_s, exner_rho_levels, exner )                      &
!$OMP& PRIVATE( i, j, k )
DO k = pdims_s%k_start,pdims_s%k_end + 1
  DO j = pdims_s%j_start,pdims_s%j_end
    DO i = pdims_s%i_start,pdims_s%i_end
      exner_rho_levels(i,j,k) = exner(i,j,k)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

IF (conserve_dry_mass /= not_conserved) THEN
  total_rho = eg_total_mass(dryrho, l_exclude_rim=.TRUE.)
END IF

IF (printstatus >= prstatus_normal) THEN
   CALL print_windmax_2(u, v, w, conserve_dry_mass, mass_fix_factor)

   IF (nsteps_consv_print > 0) THEN
      IF (MOD(timestep_number, nsteps_consv_print) == 0) THEN
         total_rho = eg_total_mass(dryrho, l_exclude_rim=.TRUE.)
         total_aam=eg_total_aam(u, dryrho, l_shallow, L_exclude_rim=.TRUE.)
         total_ke=eg_total_ke(u, v, w, dryrho, L_exclude_rim=.TRUE.)
         CALL print_conservation_diag(total_rho, total_aam, total_ke,&
                                      not_initial_timestep)
      END IF
   END IF
END IF

! DEPENDS ON: exitchek
CALL exitchek(atmos_im, lexitnow)

! ***********************************************************************
!  CONVERSION FROM EG to ND prognostics! (start)
IF (ltimer) CALL timer('AS CONVERT',5)

! update ND's q's from EG's mixing ratios
CALL update_nd_moisture(l_mr_pc2,l_mr_qtbalcld)

! derive ND's theta,p,pstar,p_theta_levels and exner_theta_levels
! from EG prognostics
CALL eg_thetav_theta                                                    &
                  (thetav, theta, m_v,                                  &
                   p, pstar, p_theta_levels,                            &
                   exner, exner_surf, exner_theta_levels)

IF (ltimer) CALL timer('AS CONVERT',6)
!  CONVERSION FROM EG to ND prognostics! (end)
!
! ***********************************************************************

! ------------------------------------------------------------------
! Section 17  Aerosol Modelling - includes Sulphur Cycle, soot, biomass
! and fossil-fuels organic carbon (OCFF)
! ------------------------------------------------------------------
!
IF (ltimer) CALL timer('AS Aerosol Modelling',5)
!
aeroif: IF (l_sulpc_so2 .OR. l_soot .OR. l_biomass .OR.                &
            l_ocff .OR. l_nitrate .OR. l_dust) THEN
  !
  ! Allocate diagnostic space for STASH
  ALLOCATE (stashwork17(stash_maxlen(17,atmos_im)))
  !
  ! Don't call Swap_bounds for fields used in Aero_Ctl_4A
  !
  ALLOCATE(o3_mmr  (tdims%i_start:tdims%i_end,                    &
                    tdims%j_start:tdims%j_end,                    &
                    tdims%k_start:tdims%k_end))
  ALLOCATE(hno3_mmr(tdims%i_start:tdims%i_end,                    &
                    tdims%j_start:tdims%j_end,                    &
                    tdims%k_start:tdims%k_end))
  ALLOCATE(h2o2_mmr(tdims%i_start:tdims%i_end,                    &
                    tdims%j_start:tdims%j_end,                    &
                    tdims%k_start:tdims%k_end))
  ALLOCATE(oh_conc (tdims%i_start:tdims%i_end,                    &
                    tdims%j_start:tdims%j_end,                    &
                    tdims%k_start:tdims%k_end))
  ALLOCATE(ho2_conc(tdims%i_start:tdims%i_end,                    &
                    tdims%j_start:tdims%j_end,                    &
                    tdims%k_start:tdims%k_end))

  ! Get the appropriate oxidant fields for CLASSIC
  CALL get_sulpc_oxidants(                                            &
    l_sulpc_online_oxidants, l_ukca, l_ukca_trop, l_ukca_tropisop,    &
    l_ukca_strattrop, l_ukca_raq, l_ukca_raqaero, l_ukca_offline,     &
    theta, p_theta_levels, exner_theta_levels,                        &
    o3_chem, h2o2_limit, oh, ho2,                                     &
    oh_ukca, ho2_ukca, h2o2_ukca, o3_ukca, hno3_ukca,                 &
    o3_mmr, hno3_mmr, h2o2_mmr, oh_conc, ho2_conc)

  CALL aero_ctl_4A(                                                &
  ! Parallel variables
          halo_i, halo_j, offx, offy, global_row_length, global_rows ,     &
          gc_proc_row_group, gc_proc_col_group,                            &
          at_extremity, nproc, nproc_x, nproc_y  ,                         &
          neighbour, g_rows, g_row_length,  mype ,                         &
  ! model dimensions
          row_length, rows, n_rows, land_field    ,                        &
          bl_levels, n_cca_lev  ,                                          &
          theta_field_size           ,                                     &
          salt_dim1, salt_dim2, salt_dim3 ,                                &
          aero_dim1, aero_dim2, aero_dim3,                                 &
  ! co-ordinate information
          delta_lambda, delta_phi ,                                        &
          lat_rot_np, long_rot_np,                                         &
  ! time stepping information
          i_year, i_month, i_day_number, i_day, i_hour, i_minute,          &
          i_second, timestep_number ,                                      &
          previous_time ,                                                  &
          call_chem_freq  ,                                                &
  ! trig arrays
          sin_theta_longitude, cos_theta_longitude  ,                      &
          fv_cos_theta_latitude ,                                          &
  ! grid-dependent arrays
              f3_at_u, true_longitude, true_latitude,                      &
  !
  ! data fields in
          u, v, tstar, tstar_sea,                                          &
          theta(:,:,1:), q(:,:,1:), qcl(:,:,1:), qcf(:,:,1:),              &
  !
  ! NOTE that this is  wet rho * r**2 at time level (n+1) (see above comment)
          wetrho_r_sq_n,                                                   &
  !
          land, frac_land,                                                 &
          p_theta_levels(:,:,1:),                                          &
          exner_rho_levels,                                                &
          exner_theta_levels(:,:,1:),                                      &
          ice_fraction, snodep  ,                                          &
          cf_bulk(:,:,1:),                                                 &
          oh_conc(:,:,1:), h2o2_mmr(:,:,1:), ho2_conc(:,:,1:),             &
          o3_mmr(:,:,1:), hno3_mmr(:,:,1:),                                &
          so2_em, so2_hilem, so2_natem(:,:,1:),                            &
          dms_em, dms_conc, nh3_em ,                                       &
          soot_em, soot_hilem, bmass_em, bmass_hilem, ocff_em, ocff_hilem, &
          bmass_hilem_h1, bmass_hilem_h2, land_index,                      &
  ! logicals in
          l_sulpc_so2, l_sulpc_dms, l_sulpc_ozone,                         &
          l_sulpc_so2_o3_nonbuffered, l_sulpc_nh3,                         &
          l_sulpc_online_oxidants,                                         &
          l_sulpc_2_way_coupling,                                          &
          l_use_sulphate_sulpc, l_use_seasalt_sulpc, l_soot,               &
          l_biomass, l_use_bmass_sulpc,                  &
          l_ocff, l_use_ocff_sulpc,                                        &
          l_nitrate, l_use_nitrate_sulpc,                                  &
          l_so2_surfem, l_so2_hilem, l_so2_natem, l_dms_em,                &
          l_dms_em_inter, i_dms_flux,                             &
          l_nh3_em, l_ctile,                                               &
          l_soot_surem, l_soot_hilem, l_bmass_surem, l_bmass_hilem,        &
          l_ocff_surem, l_ocff_hilem, l_use_biogenic,                      &
          l_use_seasalt_direct, l_use_seasalt_indirect,                    &
          l_use_seasalt_autoconv, l_use_seasalt_pm, l_dust,                &
  !
  ! data fields in/out
          so2(:,:,1:), dms(:,:,1:),                                        &
          so4_aitken(:,:,1:), so4_accu(:,:,1:), so4_diss(:,:,1:),          &
          h2o2(:,:,1:), nh3(:,:,1:),                                       &
          soot_new(:,:,1:), soot_agd(:,:,1:), soot_cld(:,:,1:),            &
          bmass_new(:,:,1:), bmass_agd(:,:,1:), bmass_cld(:,:,1:),         &
          ocff_new(:,:,1:), ocff_agd(:,:,1:), ocff_cld(:,:,1:),            &
          biogenic(:,:,1:),                                                &
          nitr_acc(:,:,1:), nitr_diss(:,:,1:),                             &
  !
  ! data fields in
          dust_div1(:,:,1:), dust_div2(:,:,1:), dust_div3(:,:,1:),         &
          dust_div4(:,:,1:), dust_div5(:,:,1:), dust_div6(:,:,1:),         &
  !
  ! data fields out
  ! diagnostic info
          stashwork17,                                                     &
  ! error info
          errorstatus )

  IF (l_sulpc_online_oxidants .AND. l_sulpc_2_way_coupling) THEN
    !
    CALL write_sulpc_oxidants(                                             &
        oh_conc, h2o2_mmr, ho2_conc, o3_mmr, hno3_mmr,                     &
        theta, p_theta_levels, exner_theta_levels,                         &
        oh_ukca, ho2_ukca, h2o2_ukca, o3_ukca, hno3_ukca)
    !
  END IF
  !
  ! Don't call Swap_bounds for updated fields
  !
  ! Diagnostics STASHed for Aerosol section 17
  !
  ! DEPENDS ON: stash
  CALL stash(atmos_sm,atmos_im,17,stashwork17,                     &
     errorstatus,cmessage)
  !
  DEALLOCATE (stashwork17)
  !
  IF (ALLOCATED(o3_mmr  )) DEALLOCATE(o3_mmr)
  IF (ALLOCATED(hno3_mmr)) DEALLOCATE(hno3_mmr)
  IF (ALLOCATED(h2o2_mmr)) DEALLOCATE(h2o2_mmr)
  IF (ALLOCATED(oh_conc )) DEALLOCATE(oh_conc)
  IF (ALLOCATED(ho2_conc)) DEALLOCATE(ho2_conc)

  ! set zero level quantities for Endgame
  ! currently set equal to level 1

  ! only copy to level 0 if each aerosol tracer is on
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE ( i, j )
  IF (l_sulpc_so2) THEN
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims_s%j_start, tdims_s%j_end
      DO i = tdims_s%i_start, tdims_s%i_end
        so2(i,j,0)        = so2(i,j,1)
        so4_aitken(i,j,0) = so4_aitken(i,j,1)
        so4_accu(i,j,0)   = so4_accu(i,j,1)
        so4_diss(i,j,0)   = so4_diss(i,j,1)
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF
  IF (l_sulpc_dms) THEN
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims_s%j_start, tdims_s%j_end
      DO i = tdims_s%i_start, tdims_s%i_end
        dms(i,j,0)        = dms(i,j,1)
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF
  IF (l_sulpc_nh3) THEN
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims_s%j_start, tdims_s%j_end
      DO i = tdims_s%i_start, tdims_s%i_end
        nh3(i,j,0)        = nh3(i,j,1)
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF
  IF (l_soot) THEN
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims_s%j_start, tdims_s%j_end
      DO i = tdims_s%i_start, tdims_s%i_end
        soot_new(i,j,0)   = soot_new(i,j,1)
        soot_agd(i,j,0)   = soot_agd(i,j,1)
        soot_cld(i,j,0)   = soot_cld(i,j,1)
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF
  IF (l_biomass) THEN
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims_s%j_start, tdims_s%j_end
      DO i = tdims_s%i_start, tdims_s%i_end
        bmass_new(i,j,0)  = bmass_new(i,j,1)
        bmass_agd(i,j,0)  = bmass_agd(i,j,1)
        bmass_cld(i,j,0)  = bmass_cld(i,j,1)
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF
  IF (l_ocff) THEN
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims_s%j_start, tdims_s%j_end
      DO i = tdims_s%i_start, tdims_s%i_end
        ocff_new(i,j,0)   = ocff_new(i,j,1)
        ocff_agd(i,j,0)   = ocff_agd(i,j,1)
        ocff_cld(i,j,0)   = ocff_cld(i,j,1)
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF
  IF (l_nitrate) THEN
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims_s%j_start, tdims_s%j_end
      DO i = tdims_s%i_start, tdims_s%i_end
        nitr_acc(i,j,0)   = nitr_acc(i,j,1)
        nitr_diss(i,j,0)  = nitr_diss(i,j,1)
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF
  IF (l_dust) THEN
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims_s%j_start, tdims_s%j_end
      DO i = tdims_s%i_start, tdims_s%i_end
        dust_div1(i,j,0)   = dust_div1(i,j,1)
        dust_div2(i,j,0)   = dust_div2(i,j,1)
      END DO
    END DO
!$OMP END DO NOWAIT
    IF ( .NOT. l_twobin_dust ) THEN
!$OMP DO SCHEDULE(STATIC)
      DO j = tdims_s%j_start, tdims_s%j_end
        DO i = tdims_s%i_start, tdims_s%i_end
          dust_div3(i,j,0)   = dust_div3(i,j,1)
          dust_div4(i,j,0)   = dust_div4(i,j,1)
          dust_div5(i,j,0)   = dust_div5(i,j,1)
          dust_div6(i,j,0)   = dust_div6(i,j,1)
        END DO
      END DO
!$OMP END DO
    END IF
  END IF
!$OMP END PARALLEL

END IF aeroif        ! END L_SULPC_SO2.OR.L_SOOT TEST
!
IF (ltimer) CALL timer('AS Aerosol Modelling',6)

IF (l_physics) THEN

  IF ( i_cld_vn > i_cld_off ) THEN

    ! Only call the cloud scheme/qt_bal_cld if i_cld_vn is not off

    ! N. B. Variables t_incr_diagnostic, q_incr_diagnostic and
    ! qcl_incr_diagnostic are allocated within this loop and therefore
    ! must only be deallocated if i_cld_vn is smith/pc2. 

    ! Are we using the PC2 cloud scheme?

    IF (i_cld_vn == i_cld_pc2) THEN

      ! ----------------------------------------------------------------------
      !  PC2: Calculate condensation due to changes in temperature resulting
      !    from adiabatic changes in pressure (mainly from vertical advection)
      !    Also call checking routine from within this subroutine.
      ! ----------------------------------------------------------------------

      CALL atm_step_alloc_4A( &
          cf_star, cfl_star, cff_star, q_star, qcl_star, qcf_star, &
          frac_control, r_u, r_v, r_w, errorstatus, 'allocPC2')

      CALL pc2_pressure_forcing (                                         &
                p,pstar,p_theta_levels(:,:,1:),                           &
                rhc_row_length, rhc_rows, zlcl_mixed,                     &
                rhcpt, theta, cf_bulk,                                    &
                cf_liquid, cf_frozen,                                     &
                q, qcl,qcf,                                               &
                exner_star,                                               &
                exner_theta_levels(:,:,1:),                               &
                ccb, cumulus, rhts, tlts, qtts, ptts,cf_area(:,:,1:),     &
                l_mr_pc2)

      DEALLOCATE(rhts)
      DEALLOCATE(tlts)
      DEALLOCATE(qtts)
      DEALLOCATE(ptts)

    END IF ! i_cld_pc2

    ! ----------------------------------------------------------------------
    ! Add ability to get increments from qt_bal_cld call and output in
    ! section 15
    ! ----------------------------------------------------------------------

    CALL atm_step_diag(8)

    IF (l_run_with_physics2 .AND.                              &
          (i_cld_vn == i_cld_smith .OR. l_pc2_reset )) THEN
      ! ----------------------------------------------------------------------
      ! Call cloud scheme to make cloud consistent with moisture fields
      ! ----------------------------------------------------------------------

      CALL qt_bal_cld( pstar,                                       &
          p_theta_levels(:,:,1:),p,                                 &
          theta,exner_theta_levels,                                 &
          q,qcl,qcf,qcf2,                                           &
          rhcpt, rhc_row_length, rhc_rows, bl_levels,               &
          fv_cos_theta_latitude,                                    &
          l_mcr_qcf2,                                               &
          l_mr_qtbalcld,                                            &
          ntml, cumulus, cf_area, cf_bulk, cf_liquid, cf_frozen,    &
          mype)

      ! calculate changes to T , q etc
      CALL atm_step_diag(9)

    END IF  ! i_cld_smith and L_pc2_reset

  END IF ! i_cld_off

  ! PV-tracer: compute dPV_cld
  IF (l_pv_tracer .OR. l_diab_tracer)                             & 
                   CALL dyn_tr_np1 (dryrho, exner_theta_levels,   &
                                    dPV_sol, dPV_mass, dPV_cld,   &
                                    dtheta_cld)

  DEALLOCATE ( rhcpt )

  ! set zero level quantities for Endgame
  ! currently set equal to level 1
!$OMP  PARALLEL DEFAULT(NONE)                                          &
!$OMP& SHARED( tdims_s, tdims_l, theta, q, qcl, qcf, cf_bulk,          &
!$OMP&         cf_liquid, cf_frozen, l_mcr_qcf2, qcf2 )                &
!$OMP& PRIVATE( i, j )
!$OMP DO SCHEDULE(STATIC)
  DO j = tdims_s%j_start, tdims_s%j_end
    DO i = tdims_s%i_start, tdims_s%i_end
      theta    (i,j,0) = theta    (i,j,1)
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO j = tdims_l%j_start, tdims_l%j_end
    DO i = tdims_l%i_start, tdims_l%i_end
      q        (i,j,0) = q        (i,j,1)
      qcl      (i,j,0) = qcl      (i,j,1)
      qcf      (i,j,0) = qcf      (i,j,1)
      cf_bulk  (i,j,0) = cf_bulk  (i,j,1)
      cf_liquid(i,j,0) = cf_liquid(i,j,1)
      cf_frozen(i,j,0) = cf_frozen(i,j,1)
    END DO
  END DO
!$OMP END DO NOWAIT

  IF (l_mcr_qcf2) THEN
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims_l%j_start, tdims_l%j_end
      DO i = tdims_l%i_start, tdims_l%i_end
        qcf2(i,j,0) = qcf2(i,j,1)
      END DO
    END DO
!$OMP END DO
  END IF
!$OMP END PARALLEL

END IF   ! L_physics

! ----------------------------------------------------------------------
! Section 9.0a Calculate Mass and energy of atmosphere if required
!              using end of timestep values
!             Only done every energy correction period.
! ----------------------------------------------------------------------
IF (ltimer) CALL timer('AS Energy mass   ',5)
IF (l_emcorr) THEN

  ! Set energy correction for use in section 30
  energy_corr_now=a_realhd(rh_energy_corr)

  IF (lenergy) THEN
    ! zero fields to be calculated

    tot_energy_final = 0.0
    tot_dry_mass_final = 0.0
    tot_moist_final = 0.0

    CALL eng_mass_diag (                                           &
                      delta_xi1,delta_xi2,                         &
                      theta , u, v, w,                             &
    ! NOTE that this is  wet rho * r**2 at time level (n+1) (see above comment)
                                wetrho_r_sq_n,                               &
    !
                                q,                                           &
                                qcl, qcf,                                    &
                                wet_to_dry_n,                                &
                                exner_theta_levels,                          &
                                net_mflux,                                   &
                                a_realhd(rh_tot_mass_init),                  &
                                a_realhd(rh_tot_m_init),                     &
    ! logical to indicate mass and moist correction required
                                Lmass_corr,Lqt_corr,Lemq_print,              &
    ! energy correction               
                                a_energysteps,                               &
    ! IN/OUT  results from calculations
                                tot_energy_final,tot_dry_mass_final,         &
                                tot_moist_final)

    CALL cal_eng_mass_corr_4A (                                    &
                      a_energysteps,                               &
                      net_flux,                                    &
                      a_realhd(rh_tot_mass_init),                  &
                      a_realhd(rh_tot_energy_init),                &
                      a_realhd(rh_energy_corr),                    &
                      tot_energy_final )

    ! Swap initial energy and final energy.

    a_realhd(rh_tot_energy_init) = tot_energy_final

    ! Swap initial moisture and final moisture.

    a_realhd(rh_tot_m_init) = tot_moist_final

  END IF   ! LENERGY

ELSE
  ! Set energy correction for use in section 30
  energy_corr_now=0.0
END IF     ! L_EMCORR

IF (ltimer) CALL timer('AS Energy mass   ',6)
IF (ltimer) CALL timer('AS IAU',5)

!----------------------------------------------------------------------
! Incremental Analysis Update (IAU).
!----------------------------------------------------------------------

IF (l_iau .AND. stepim(atmos_im) >= IAU_FirstCallTS .AND. &
                stepim(atmos_im) <= IAU_LastCallTS) THEN

  IF (model_type == mt_lam .AND. l_fix_iau_rim_density) THEN                  
    ! The dryrho rim values were updated during the call to update_lam_lbcs,  
    ! but the wetrho_r_sq_n rim values have not yet been made consistent with 
    ! them. We need wetrho_r_sq_n to be up-to-date before entering the IAU    
    ! routine.                                                                
    CALL update_wetrho_r_sq(time_level=tl_np1_2)                              
  END IF                                                                      

  ! DEPENDS ON: iau
  CALL iau (                                      &
             l_mr_iau,                            & ! in
             frac_typ,                            & ! in
             snowdepth,            nsnow,         & ! in
             vol_smc_wilt,         vol_smc_crit,  & ! in
             vol_smc_sat,                         & ! in
             u,      v,            w,             & ! inout
             theta,                               & ! inout
  ! NOTE that this is  wet rho * r**2 at time level (n+1) (see above comment)
                     wetrho_r_sq_n,                       & ! inout
  !
                     murk,                                & ! inout
                     q,      qCL,          qCF,           & ! inout
                     TStar,  TStar_tile,   Deep_soil_temp,& ! inout
                     smcl,   tsnowlayer,   tstar_anom,    & ! inout
                     Pstar,  p,                           & ! inout
                     p_theta_levels,                      & ! inout
                     exner_rho_levels,                    & ! inout
                     exner_theta_levels,                  & ! inout
                     snodep,                              & ! inout
                     cf_area,                             & ! inout
                     cf_bulk,                             & ! inout
                     cf_liquid,                           & ! inout
                     cf_frozen,                           & ! inout
                     dust_div1, dust_div2, dust_div3,     & ! inout
                     dust_div4, dust_div5, dust_div6,     & ! inout
                     ozone_tracer, SO2, tracer_ukca )       ! inout

  CALL init_etadot()

END IF

IF (ltimer) CALL timer('AS IAU',6)

! Are we using the PC2 cloud scheme to determine area cloud fraction?

IF (i_cld_vn == i_cld_pc2 .AND. .NOT. l_pc2_reset) THEN
  IF (i_cld_area == acf_off) THEN

    ! ----------------------------------------------------------------------
    ! PC2: Set area cloud fraction to the bulk cloud fraction. Use the
    !    D1 arrays directly
    ! ----------------------------------------------------------------------

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                      &
!$OMP& SHARED( model_levels, j_start, j_stop, i_start, i_stop,         &
!$OMP&         cf_area, cf_bulk )                                      &
!$OMP& PRIVATE( i, j, k )
    DO k = 1, model_levels
      DO j = j_start, j_stop
        DO i = i_start, i_stop
          cf_area(i,j,k) = cf_bulk(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

  ELSE IF (i_cld_area == acf_brooks) THEN

    ALLOCATE ( cf_bulk_nohalo(row_length,rows,model_levels) )
    ALLOCATE ( cf_liquid_nohalo(row_length,rows,model_levels) )
    ALLOCATE ( cf_frozen_nohalo(row_length,rows,model_levels) )

    ! Place bulk, liquid and frozen cloud fractions in halo-free arrays
    ! Use indexing over the full row and row_length (including any LAM
    ! boundary rim) since the call to ls_acf_brooks uses this indexing.
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          cf_bulk_nohalo(i,j,k)   = cf_bulk(i,j,k)
          cf_liquid_nohalo(i,j,k) = cf_liquid(i,j,k)
          cf_frozen_nohalo(i,j,k) = cf_frozen(i,j,k)
        END DO
      END DO
    END DO

    CALL ls_acf_brooks (                                      &
       fv_cos_theta_latitude,                                 &
       cf_bulk_nohalo, cf_liquid_nohalo,                      &
       cf_frozen_nohalo, cumulus, cf_area )

    DEALLOCATE ( cf_bulk_nohalo )
    DEALLOCATE ( cf_liquid_nohalo )
    DEALLOCATE ( cf_frozen_nohalo )

  END IF ! i_cld_area

END IF  ! i_cld_pc2 and not L_pc2_reset
!
! ***********************************************************************
!
!  CONVERSION FROM ND to EG prognostics! (start)
IF (ltimer) CALL timer('AS CONVERT',5)
!
! convert moisture back to EG, as PC2 may have changed it:
CALL update_eg_moisture(l_mr_pc2,l_mr_qtbalcld)
!
!  CONVERSION FROM ND to EG prognostics! (end)
!
! ***********************************************************************

! ***********************************************************************
!
!  CONVERSION FROM EG prognostic to ND prognostics! (start)
!  depending on how keen we are to keep the ND quantities derived we will
!  have to update this here now, since EG moisture has changed. It should
!  of course be thermodynamically consistent already! So as an optimisation
!  we may want to skip this conversion

! compute the new wet rho * r**2 (derived from EG density and EG moisture)
CALL update_wetrho_r_sq(time_level=tl_np1_2)

IF (ltimer) CALL timer('AS CONVERT',6)
!  CONVERSION FROM EG prognostic to ND prognostics! (end)
!
! ***********************************************************************

! Nudging with analysis data
nudgingif: IF ( L_nudging ) THEN

  IF (ltimer) CALL timer('AS Nudging',5)
  IF (model_type /= mt_global) THEN
    ErrorStatus = 39
    cmessage = 'Nudging not yet tested for Limited Area Model'

    CALL ereport('ATM_STEP',ErrorStatus,cmessage)
  END IF

  ! Allocate stashwork array
  CALL Atm_Step_diag(39)

  CALL nudging_main1 (                                             &
  ! time stepping information.
         i_year, i_month, i_day, i_hour, i_minute,                   &
         i_second, timestep_number,                                  &

  ! in data fields.
         theta(:,:,1), u, v, p, exner_theta_levels, p_theta_levels,  &
  ! out updated STASH
          STASHwork39)

  theta(:,:,0) = theta(:,:,1)      ! copy lev 1 to lev 0

  ! Swap bounds for modified fields
  CALL swap_bounds(theta,                                                &
                   tdims_s%i_len - 2*tdims_s%halo_i,                     &
                   tdims_s%j_len - 2*tdims_s%halo_j,                     &
                   tdims_s%k_len,                                        &
                   tdims_s%halo_i, tdims_s%halo_j,                       &
                   fld_type_p,swap_field_is_scalar)
  CALL swap_bounds(u,                                                    &
                   udims_s%i_len - 2*udims_s%halo_i,                     &
                   udims_s%j_len - 2*udims_s%halo_j,                     &
                   udims_s%k_len,                                        &
                   udims_s%halo_i, udims_s%halo_j,                       &
                   fld_type_u,swap_field_is_vector)
  CALL swap_bounds(v,                                                    &
                   vdims_s%i_len - 2*vdims_s%halo_i,                     &
                   vdims_s%j_len - 2*vdims_s%halo_j,                     &
                   vdims_s%k_len,                                        &
                   vdims_s%halo_i, vdims_s%halo_j,                       &
                   fld_type_v,swap_field_is_vector)

  IF (ltimer .OR. lstashdumptimer) CALL timer('AS STASH',5)

  ! Copy diagnostics into main D1 array
  ! DEPENDS ON: Atm_Step_stash
  CALL Atm_Step_stash(                                           &
         errorstatus,39)

  IF (ltimer .OR. lstashdumptimer) CALL timer('AS STASH',6)
  IF (ltimer) CALL timer('AS Nudging',6)

END IF nudgingif

! need to recompute thetav for EG
! ***********************************************************************
!
!  CONVERSION FROM ND to EG prognostics! (start)
!
IF (ltimer) CALL timer('AS CONVERT',5)

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                      &
!$OMP& SHARED( tdims, thetav, theta, m_v, recip_epsilon )              &
!$OMP& PRIVATE( i, j, k )
DO k = tdims%k_start, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      thetav(i,j,k) = theta(i,j,k)*(1.0+m_v(i,j,k)*recip_epsilon)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

CALL swap_bounds(thetav,                                               &
                 tdims_s%i_len - 2*tdims_s%halo_i,                     &
                 tdims_s%j_len - 2*tdims_s%halo_j,                     &
                 tdims_s%k_len,                                        &
                 tdims_s%halo_i, tdims_s%halo_j,                       &
                 fld_type_p,swap_field_is_scalar)
                 
IF (ltimer) CALL timer('AS CONVERT',6)                 
!
!  CONVERSION FROM ND to EG prognostics! (end)
!
! ***********************************************************************

! simulate dumping
IF (ldump) THEN
  subseq_to_dump = .TRUE.
  CALL reset_dpt_pts()
END IF

! ----------------------------------------------------------------------
! Section 9.0 Diagnostics at end of timestep
! ----------------------------------------------------------------------

IF (ltimer) CALL timer('AS End TStep Diags',5)

! PV-tracer: Compute dPV at the end of the timestep
IF (l_pv_tracer) CALL dyn_tr_end( dryrho, exner_theta_levels,      &
                                  dPV_iau, dPV_nud, dPV_tot)

! section 15: 'dynamics' based quantities
IF (      sf(0,15) .AND. errorstatus == 0) THEN

  ! DEPENDS ON: st_diag1
  CALL st_diag1(stash_maxlen(15,atmos_im),                         &
 t_incr_diagnostic,q_incr_diagnostic,qcl_incr_diagnostic,          &
 errorstatus,cmessage)

END IF ! SF(0,15)

! Cleanup arrays holding increments

IF ( i_cld_vn > i_cld_off ) THEN

  ! Variables t_incr_diagnostic, q_incr_diagnostic and
  ! qcl_incr_diagnostic are allocated in atm_step_diag
  ! above in the i_cld_vn loop above, with flag
  ! numbers 8 and 9. Therefore these variables must
  ! only be deallocated if i_cld_vn is smith/pc2

  IF (l_physics .AND. errorstatus == 0) THEN
    DEALLOCATE(t_incr_diagnostic)
    DEALLOCATE(q_incr_diagnostic)
    DEALLOCATE(qcl_incr_diagnostic)
  END IF
END IF

! section 16: 'physics' based quantities
IF (      sf(0,16) .AND. errorstatus == 0) THEN

  ! DEPENDS ON: st_diag2
  CALL st_diag2(stash_maxlen(16,atmos_im),                        &
                 errorstatus,cmessage)

END IF !   SF(0,16)

IF (i_cld_vn == i_cld_pc2 .AND. errorstatus == 0 .AND. l_physics) THEN

  CALL atm_step_alloc_4A( &
       cf_star, cfl_star, cff_star, q_star, qcl_star, qcf_star, &
       frac_control, r_u, r_v, r_w, errorstatus, 'dllocPC2')

END IF  ! i_cld_pc2 and ErrorStatus == 0

! section 30: climate diagnostics (also required for PWS updraft helicity)
IF ( (sf(0,30) .OR. flag_upd_helicity_5k) .AND. ErrorStatus == 0) THEN

  ! Calculation of total increments
  ! NOTE - this must be after all processes which alter model prognostics
  CALL eot_incs_calc()

  ! size of diagnostic space
  ALLOCATE (STASHwork30(STASH_maxlen(30,atmos_im)))

  ! DEPENDS ON: st_diag3
  CALL St_diag3(STASHwork30, STASH_maxlen(30,atmos_im),            &
    energy_corr_now, sin_v_latitude, wet_to_dry_n, co2_mmr,        &
    ErrorStatus,CMessage)

  CALL eot_incs_dealloc()

  DEALLOCATE (STASHwork30) ! Clear space

END IF ! SF(0,30)

! Check error condition
IF (errorstatus >  0) THEN
  CALL ereport(routinename,errorstatus,cmessage)
END IF

! PWS diagnostics (migrated from FieldCalc)
! Moved to after st_diag3 in case updraft helicity requested
IF (sf(0,stashcode_pws_sec) .AND. ErrorStatus == 0) THEN

  CALL pws_diags_driver(stash_maxlen(stashcode_pws_sec,atmos_im), &
                          errorstatus,cmessage)

END IF !   SF(0,stashcode_pws_sec)

! Check error condition
IF (errorstatus >  0) THEN
  CALL ereport(routinename,errorstatus,cmessage)
END IF

IF (ltimer) CALL timer('AS End TStep Diags',6)
IF (ltimer .OR. lstashdumptimer) CALL timer('AS STASH',5)

! DEPENDS ON: Atm_Step_stash
CALL atm_step_stash( &
         errorstatus, 4)

IF (ltimer .OR. lstashdumptimer) CALL timer('AS STASH',6)

IF ( l_diag_print ) THEN
  IF (ltimer) CALL timer('AS DIAGNOSTICS',5)
  CALL print_diag_4A(                                            &
                 u, v, theta,                                    &
  ! NOTE that this is  wet rho * r**2 at time level (n+1) (see above comment)
                          wetrho_r_sq_n,                                  &
  !
                          w, q, qcl, qcf,                                 &
                          rows, n_rows, row_length, model_levels,         &
                          offx, offy, timestep_number,                    &
                          rpemax, rpemin, ipesum, rpesum,                 &
                          max_w_run, max_wind_run, min_theta1_run,        &
                          dtheta1_run, max_div_run, min_div_run,          &
                          min_lapse_run, max_shear_run, time_max_shear,   &
                          time_div_max, time_div_min, time_lapse_min,     &
                          time_w_max, time_max_wind, time_theta1_min,     &
                          max_KE_run, min_KE_run, max_noise_run,          &
                          time_KE_max, time_KE_min, time_noise_max )
  IF ( l_flush6 ) CALL umPrintFlush()
  IF (ltimer) CALL timer('AS DIAGNOSTICS',6)
END IF     !  L_diag_print

! ----------------------------------------------------------------------
! Diagnostics at end of timestep
! ----------------------------------------------------------------------

! Logical first_atmstep_call is true on the first call to ATM_STEP_4A
! and is set to false at the end of ATM_STEP_4A (uses the SAVE command)
first_atmstep_call = .FALSE.

! Write random seed to dump headers
IF (ldump) THEN
  IF (l_rp2 .OR. l_skeb2 .OR. l_spt .OR. i_pert_theta /= off ) THEN
    CALL stph_seed_copy_to_dump
  END IF
END IF

IF (l_ukca) THEN

  IF (ltimer) CALL timer('AS UKCA_MAIN1',5)
  CALL ukca_main1(tracer_ukca, q)

  IF (ltimer) CALL timer('AS UKCA_MAIN1',6)

  IF (l_ukca_h2o_feedback) THEN
    ! Update m_v as q has changed
    CALL eg_q_to_mix(tdims_l,tdims_s,                                      &
                     q, qcl, qcf, qcf2,                                    &
                     qrain, qgraup,                                        &
                     l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,                &
                     m_v, m_cl, m_cf,                                      &
                     m_cf2, m_r, m_gr,.TRUE.)

    thetav = theta * (1.0 + m_v*recip_epsilon)

  END IF    ! l_ukca_h2o_feedback

END IF     ! l_ukca

IF (Lexitnow) THEN
  CALL eg_destroy_vert_damp()
  CALL eg_destroy_helmholtz()
  CALL destroy_wet_to_dry_n()
  CALL destroy_track_arr()
  IF (l_skeb2) CALL destroy_r_skeb()
  IF (l_spt)   CALL destroy_r_spt()
  IF (l_spt .OR. l_pv_tracer .OR. l_diab_tracer)                          &
    CALL destroy_phy_tend()
END IF

IF (print_runtime_info) CALL atm_step_info

IF (ltimer) CALL timer('Atm_Step_4A (AS)',6)
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)

END SUBROUTINE atm_step_4A
