! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Boundary layer Implicit solver and Large-scale (Area) Cloud Scheme.


! Code description:
!   language: Fortran 95
!   this code is written to UMDP3 programming standards.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary layer
MODULE ni_imp_ctl_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'NI_IMP_CTL_MOD'
CONTAINS

SUBROUTINE ni_imp_ctl (                                                 &
! IN model dimensions.
  rhc_row_length, rhc_rows, land_points,                                &
  ntiles, bl_levels, cloud_levels, n_cca_levels,                        &
! IN Model switches
  cycleno,                                                              &
! IN model Parameters
  tr_vars,                                                              &
! IN trig arrays
  xx_cos_theta_latitude,                                                &
! IN data fields.
  p_layer_centres, p_layer_boundaries, rho_wet_rsq, rho_wet_tq,         &
  u, v, w,                                                              &
  land_sea_mask, q, qcl, qcf, p_star, theta, exner_theta_levels,        &
! IN ancillary fields and fields needed to be kept from tstep to tstep
  sil_orog_land, ho2r2_orog,                                            &
  ice_fract, di_ncat, ice_fract_ncat, k_sice, u_0, v_0,land_index,      &
  cca, lcbase, ccb0, cct0,                                              &
  ls_rain, ls_snow, conv_rain, conv_snow, l_scrn, l_plsp,               &
! IN variables required from BDY_LAYR in IMP_SOLVER
  alpha1_sea, alpha1_sice, ashtf_sea, ashtf, bq_gb, bt_gb,              &
  dtrdz_charney_grid, rdz_charney_grid, dtrdz_u, dtrdz_v, rdz_u, rdz_v, &
  cdr10m_u, cdr10m_v, cdr10m_n_u, cdr10m_n_v, cd10m_n_u, cd10m_n_v,     &
  z_theta, k_blend_tq, k_blend_uv, uStarGBM, rhokm, rhokm_u, rhokm_v,   &
! IN diagnostics (started or from) BDY_LAYR
  rib_gb, zlcl, zht, zhnl, dzh, qcl_inv_top, zh,                        &
  bl_type_1,bl_type_2,bl_type_3,bl_type_4,bl_type_5,bl_type_6,          &
  bl_type_7, z0m_gb, z0m_eff_gb, z0h_eff_gb,                            &
  ntml, cumulus, l_pc2_diag_sh_pts,                                     &
! IN data required for tracer mixing :
  rho_aresist,aresist,r_b_dust,                                         &
  kent, we_lim, t_frac, zrzi,                                           &
  kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc,                           &
  zhsc,z_half,dust_flux,dust_emiss_frac,                                &
  u_s_t_tile,u_s_t_dry_tile,u_s_std_tile,                               &
! IN additional variables for JULES
  tile_pts,tile_index,tile_frac,canopy,                                 &
  alpha1,fraca,rhokh_tile,smc,chr1p5m,resfs,z0hssi,z0mssi,              &
  canhc_tile,flake,wt_ext_tile,lw_down,lai_ft,canht_ft,                 &
  sw_tile,ashtf_tile,gc,aresist_tile,                                   &
  resft,rhokh_sice, rhokh_sea,z0h_tile,z0m_tile,chr1p5m_sice,           &
  fland,flandg,flandg_u,flandg_v,vshr_land,vshr_ssi,                    &
  emis_tile, t_soil, snow_tile, rib_ssi,                                &
! IN JULES variables for STASH
  Gs,gpp,npp,resp_p,gpp_ft,npp_ft,resp_p_ft,resp_s,                     &
  resp_s_tot,cs,                                                        &
  rib_tile,fsmc,catch,g_leaf,                                           &
  co2_emits, co2flux,                                                   &
!IN additional variables for soil moisture nudging scheme
  wt_ext,                                                               &
! INOUT diagnostic info
 stashwork3, stashwork9,                                                &
! SCM diagnostics (dummy in full UM) & bl_diags
  nSCMDpkgs, L_SCMDiags, BL_diag, sf_diag,                              &
! INOUT (Note ti and ti_gb are IN only if l_sice_multilayers=T)
  TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans,                               &
  cca0,ccw0,cca0_2d, fqw, ftl, taux, tauy, rhokh,                       &
  fqw_ice,ftl_ice,dtstar_tile,dtstar_sea,dtstar_sice, ti,               &
  cf_area, cf_bulk,                                                     &
  t_latest, q_latest, qcl_latest, qcf_latest,                           &
  cf_latest, cfl_latest, cff_latest,                                    &
  r_u, r_v, r_w, cf_liquid, cf_frozen,                                  &
  sum_eng_fluxes,sum_moist_flux,rhcpt,                                  &
! INOUT tracer fields
  aerosol, free_tracers, resist_b,  resist_b_tile,                      &
  dust_div1,dust_div2,dust_div3,dust_div4,dust_div5,dust_div6,          &
  drydep2, so2, dms, so4_aitken, so4_accu, so4_diss, nh3,               &
  soot_new, soot_aged, soot_cld, bmass_new, bmass_agd,                  &
  bmass_cld, ocff_new, ocff_aged, ocff_cld, nitr_acc, nitr_diss,        &
  co2, ozone_tracer,                                                    &
! INOUT additional variables for JULES
  tstar_tile,fqw_tile,epot_tile,ftl_tile,                               &
  radnet_sea,radnet_sice,olr, tstar_sice_cat,tstar_ssi,tstar_sea,       &
  taux_land,taux_ssi,tauy_land,tauy_ssi,error_code,                     &
! OUT fields
  surf_ht_flux_land, zlcl_mixed,                                        &
  theta_star_surf,qv_star_surf,                                         &
! OUT additional variables for JULES
  tstar, ti_gb, ext, snowmelt, tstar_land,tstar_sice, ei_tile,          &
  ecan_tile,melt_tile, surf_htf_tile                                    &
  )

#if !defined(LFRIC)
USE swap_bounds_mv_mod, ONLY: swap_bounds_mv
USE swap_bounds_2d_mv_mod, ONLY: swap_bounds_2d_mv
USE mpp_conf_mod,   ONLY: swap_field_is_scalar
USE diagnostics_bl_mod, ONLY: diagnostics_bl
USE ls_calc_rhcrit_mod, ONLY: ls_calc_rhcrit
#endif
USE atm_step_local, ONLY: dim_cs1, dim_cs2
USE atm_fields_bounds_mod, ONLY:                                        &
    udims, udims_s, vdims, vdims_s, tdims, tdims_s, tdims_l,            &
    pdims, pdims_s, wdims, wdims_s, ScmRowLen, ScmRow
USE bl_option_mod, ONLY:                                                &
    fric_heating, alpha_cd,trweights1, on, off,                         &
    l_bl_mix_qcf, kprof_cu, l_quick_ap2, one_third
USE bl_diags_mod, ONLY: strnewbldiag
USE cv_run_mod, ONLY: l_param_conv,                                     &
    tice, rad_cloud_decay_opt, l_ccrad, l_3d_cca, l_pc2_diag_sh
USE cv_param_mod, ONLY: rad_decay_off
USE cloud_inputs_mod, ONLY:                                             &
    l_fixbug_pc2_mixph, forced_cu, forced_cu_fac,                       &
    i_cld_area, i_rhcpt, i_cld_vn, rhcrit
USE dust_parameters_mod, ONLY: ndiv, ndivh, rhop, drep,                 &
    l_twobin_dust, l_dust
USE dynamics_input_mod, ONLY: l_check_moist_inc, numcycles
USE dynamics_testing_mod, ONLY: l_dry
USE dyn_coriolis_mod, ONLY: f3_at_u
USE eng_corr_inputs_mod, ONLY: l_emcorr
USE free_tracers_inputs_mod, ONLY: l_bl_tracer_mix
USE field_types, ONLY: fld_type_p, fld_type_u, fld_type_v
USE gen_phys_inputs_mod, ONLY: l_mr_physics
USE jules_sea_seaice_mod, ONLY: nice, nice_use
USE jules_soil_mod, ONLY: sm_levels
USE jules_surface_mod, ONLY: IScrnTDiag, IP_ScrnDecpl2, IP_ScrnDecpl3
USE jules_surface_types_mod, ONLY: ntype, npft
USE leonard_incs_mod, ONLY: thetal_inc_leonard, qw_inc_leonard
USE level_heights_mod, ONLY:                                            &
    r_rho_levels, r_theta_levels
USE model_domain_mod, ONLY: model_type, mt_global, mt_single_column
USE mphys_constants_mod, ONLY: mprog_min
USE murk_inputs_mod, ONLY: l_murk
USE pc2_constants_mod,   ONLY: ls_bl0, rhcpt_horiz_var, rhcpt_off,      &
     acf_off, acf_cusack, acf_brooks, cbl_and_cu, cloud_rounding_tol,   &
     i_cld_off, i_cld_smith, i_cld_pc2
USE planet_constants_mod, ONLY: cp, g, lcrcp
USE physics_tendencies_mod,  ONLY:                                      &
    init_bl_tendencies, dt_bl, du_bl, dv_bl, dq_cl_bl,                  &
    l_retain_bl_tendencies, l_retain_q_cl_bl_tendencies
USE pws_diags_mod, ONLY: pws_dustmmr1_em, pws_dustmmr2_em,              &
                         pws_dustmmr3_em, pws_dustmmr4_em,              &
                         pws_dustmmr5_em, pws_dustmmr6_em,              &
                         flag_dustmmr_em
USE run_aerosol_mod, ONLY: l_sulpc_so2,                                 &
    l_soot, l_biomass, l_ocff, l_nitrate
USE rad_pcf, ONLY: ip_cloud_mix_max
USE stochastic_physics_run_mod, ONLY: i_pert_theta, pert_theta_star,    &
    pert_theta_and_moist, pert_theta_mag
USE sf_diags_mod, ONLY: strnewsfdiag
USE s_scmop_mod,  ONLY: default_streams,                                &
    t_avg, t_mult, d_sl, d_tile, scmdiag_surf, scmdiag_land
USE scmoutput_mod,ONLY: scmoutput
USE swapable_field_mod, ONLY:                                           &
    swapable_field_pointer_type
USE stash_array_mod, ONLY: sf
USE timestep_mod, ONLY: timestep
USE turb_diff_mod, ONLY:                                                &
    l_subfilter_vert, l_leonard_term
USE tr_mix_mod, ONLY: tr_mix
USE carbon_options_mod, ONLY: l_co2_interactive
USE um_parparams, ONLY: PNorth, PSouth
USE um_parvars,   ONLY: at_extremity
USE u_to_p_mod, ONLY: u_to_p
USE v_to_p_mod, ONLY: v_to_p
USE water_constants_mod, ONLY: lc
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE flux_diag_mod, ONLY: flux_diag
USE r2_calc_total_cloud_cover_mod, ONLY: r2_calc_total_cloud_cover
USE ls_acf_brooks_mod, ONLY: ls_acf_brooks
USE ls_arcld_mod, ONLY: ls_arcld
USE pc2_delta_hom_turb_mod, ONLY: pc2_delta_hom_turb
USE pc2_hom_arcld_mod, ONLY: pc2_hom_arcld
USE bl_lsp_mod, ONLY: bl_lsp
USE bl_trmix_dd_mod, ONLY: bl_trmix_dd
USE gravsett_mod, ONLY: gravsett

IMPLICIT NONE

! arguments with intent IN. ie: input variables.
! Parallel setup variables

! Model dimensions
INTEGER, INTENT(IN) ::                                                  &
  rhc_row_length,                                                       &
  rhc_rows,                                                             &
  land_points,                                                          &
                  ! IN No.of land points being processed, can be 0.
  ntiles,                                                               &
                  ! IN No. of land-surface tiles ( JULES )
  bl_levels,                                                            &
  cloud_levels,                                                         &
  n_cca_levels
                    ! Number of levels for conv cloud amount :
!                       1 for 2D, nlevs for 3D.

INTEGER, INTENT(IN) ::                                                  &
 cycleno
                 ! Sweep number

INTEGER, INTENT(IN) ::                                                  &
  tr_vars
                           ! IN number of free tracer variables

! Trig arrays
REAL, INTENT(IN) ::                                                     &
  xx_cos_theta_latitude (tdims_s%i_start:tdims_s%i_end,                 &
                         tdims_s%j_start:tdims_s%j_end)

! Data arrays
REAL, INTENT(IN) ::                                                     &
  u(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end,        &
      udims_s%k_start:udims_s%k_end),                                   &
  v(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end,        &
      vdims_s%k_start:vdims_s%k_end),                                   &
  w(wdims_s%i_start:wdims_s%i_end,wdims_s%j_start:wdims_s%j_end,        &
    wdims_s%k_start:wdims_s%k_end),                                     &
  rho_wet_rsq(pdims_s%i_start:pdims_s%i_end,                            &
              pdims_s%j_start:pdims_s%j_end,                            &
              pdims_s%k_start:pdims_s%k_end),                           &
  rho_wet_tq(tdims%i_start:tdims%i_end,                                 &
             tdims%j_start:tdims%j_end,                                 &
             1:tdims%k_end-1),                                          &
  p_layer_centres(tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end,0:tdims%k_end),             &
  p_layer_boundaries(pdims%i_start:pdims%i_end,                         &
                     pdims%j_start:pdims%j_end,0:pdims%k_end),          &
            ! pressure at layer boundaries. Same as p except at
            ! bottom level = pstar, and at top = 0.
  p_star(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
  theta(tdims_s%i_start:tdims_s%i_end,                                  &
        tdims_s%j_start:tdims_s%j_end,                                  &
        tdims_s%k_start:tdims_s%k_end),                                 &
  exner_theta_levels(tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end),                    &
  q(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,        &
    tdims_l%k_start:tdims_l%k_end),                                     &
  qcl(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,      &
      tdims_l%k_start:tdims_l%k_end),                                   &
  qcf(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,      &
      tdims_l%k_start:tdims_l%k_end),                                   &
  ls_rain(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
  ls_snow(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
  conv_rain(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
  conv_snow(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

LOGICAL, INTENT(IN) ::                                                  &
  land_sea_mask(pdims%i_start:pdims%i_end,                              &
                pdims%j_start:pdims%j_end),                             &
  l_scrn,                                                               &
                               ! Logical to control output
                               !    of screen level T,Q,QCL,QCF
  l_plsp                   ! Logical to control output
                               !    of Probability of LS Precip

! ancillary arrays and fields required to be saved from timestep to
! timestep.
INTEGER, INTENT(IN) ::                                                  &
  land_index (land_points)      ! set from land_sea_mask

REAL, INTENT(IN) ::                                                     &
  u_0(udims%i_start:udims%i_end,udims%j_start:udims%j_end),             &
                              ! set to zero
  v_0(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),             &
                              ! set to zero
  sil_orog_land (land_points),                                          &
                                  ! orog/qrparm.orog.as
  ho2r2_orog (land_points),                                             &
                               ! orog/qrparm.orog.h2root2
  ice_fract (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                   ! ice/qrclim.ice.(month)
  ice_fract_ncat (pdims%i_start:pdims%i_end,                            &
                  pdims%j_start:pdims%j_end, nice),                     &
                            !ice fract on categories
  di_ncat(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,          &
          nice),                                                        &
                            ! ice thickness on categories
  k_sice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,nice)
                            ! sea ice effective conductivity in
                            !  sfc layer on categories (W/m2/K)

LOGICAL, INTENT(IN) ::                                                  &
 cumulus (tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),         &
                                 ! *APL bl convection flag
 l_pc2_diag_sh_pts(tdims%i_start:tdims%i_end,                           &
                   tdims%j_start:tdims%j_end)
              ! Carry diagnostic shallow convective information for PC2

INTEGER, INTENT(IN) ::                                                  &
 ntml (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

REAL, INTENT(IN) ::                                                     &
 cca    (tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         n_cca_levels)

INTEGER, INTENT(IN)  ::                                                 &
 lcbase (tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
 ccb0   (tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
 cct0   (tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)

! in variables passed from BDY_LAYR to IMP_SOLVER
REAL, INTENT(IN) ::                                                     &
  alpha1_sea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
  alpha1_sice(pdims%i_start:pdims%i_end,                                &
              pdims%j_start:pdims%j_end,nice_use),                      &
  ashtf_sea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
  ashtf(pdims%i_start:pdims%i_end,                                      &
        pdims%j_start:pdims%j_end,nice_use),                            &
  bq_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,            &
        bl_levels),                                                     &
  bt_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,            &
        bl_levels),                                                     &
  dtrdz_charney_grid(pdims%i_start:pdims%i_end,                         &
                     pdims%j_start:pdims%j_end,                         &
                     bl_levels),                                        &
  rdz_charney_grid(pdims%i_start:pdims%i_end,                           &
                   pdims%j_start:pdims%j_end,                           &
                   bl_levels),                                          &
  dtrdz_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,          &
                          bl_levels),                                   &
  dtrdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,          &
                            bl_levels),                                 &
  rdz_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,            &
                        2:bl_levels),                                   &
  rdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,            &
                          2:bl_levels),                                 &
  cdr10m_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),        &
  cdr10m_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),        &
  cdr10m_n_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),      &
  cdr10m_n_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),      &
  cd10m_n_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),       &
  cd10m_n_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),       &
  z_theta(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,          &
          pdims%k_start:pdims%k_end)

INTEGER, INTENT(IN) ::                                                  &
 k_blend_tq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
 k_blend_uv(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end)

REAL, INTENT(IN) :: uStarGBM(pdims%i_start:pdims%i_end,                 &
                              pdims%j_start:pdims%j_end)
!       ! GBM surface friction velocity for diagnosis of decoupling

REAL, INTENT(IN) ::                                                     &
  rhokm  (pdims_s%i_start:pdims_s%i_end, pdims_s%j_start:pdims_s%j_end, &
                           bl_levels),                                  &
  rhokm_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,          &
                           bl_levels),                                  &
  rhokm_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,          &
                           bl_levels)
  ! NOTE: in passing these rhokm arrays from atmos_physics2 into 
  ! ni_imp_ctl, we are changing the vertical level indexing from 
  ! 0:bl_levels-1 to 1:bl_levels.  This puts them back on the same 
  ! vertical indexing as they are in bdy_layr.

! in diagnostics (started or from) BDY_LAYR
REAL, INTENT(IN) ::                                                     &
  rib_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                  ! Mean bulk Richardson number for
!                                     lowest layer.
    zlcl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
    zht(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                                  ! Max height of turb mixing
    zhnl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                  ! IN non-local BL depth
    dzh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                                  ! IN inversion thickness
    qcl_inv_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),   &
                                  ! IN Parcel water content at inv top
    zh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                                  ! IN max of local and non-local BL depths
    bl_type_1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                                  ! IN Indicator set to 1.0 if stable
!                                 !     b.l. diagnosed, 0.0 otherwise.
    bl_type_2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                                  ! IN Indicator set to 1.0 if Sc over
!                                 !     stable surface layer diagnosed,
!                                 !     0.0 otherwise.
    bl_type_3(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                                  ! IN Indicator set to 1.0 if well
!                                 !     mixed b.l. diagnosed,
!                                 !     0.0 otherwise.
    bl_type_4(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                                  ! IN Indicator set to 1.0 if
!                                 !     decoupled Sc layer (not over
!                                 !     cumulus) diagnosed,
!                                 !     0.0 otherwise.
    bl_type_5(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                                  ! IN Indicator set to 1.0 if
!                                 !     decoupled Sc layer over cumulus
!                                 !     diagnosed, 0.0 otherwise.
    bl_type_6(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                                  ! IN Indicator set to 1.0 if a
!                                 !     cumulus capped b.l. diagnosed,
!                                 !     0.0 otherwise.
    bl_type_7(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                                  ! IN Indicator set to 1.0 if a
!                                 !     shear-dominated b.l. diagnosed,
!                                 !     0.0 otherwise.
    z0m_eff_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
    z0h_eff_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                  ! IN Effective grid-box roughness
                                  !     lengths for diagnostics

REAL, INTENT(IN) ::                                                     &
  rho_aresist(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                              ! IN RHOSTAR*CD_STD*VSHR
                              !     for CLASSIC aerosol scheme
  aresist(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                              ! IN 1/(CD_STD*VSHR)
                              !     for CLASSIC aerosol scheme
  r_b_dust(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,         &
           ndiv),                                                       &
                                      ! IN surf layer res for dust
  we_lim(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),        &
                                  ! IN rho*entrainment rate implied by
!                                   !     placing of subsidence
    zrzi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),        &
                                    ! IN (z-z_base)/(z_i-z_base)
    t_frac(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),      &
                                    ! IN a fraction of the timestep
    we_lim_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),  &
!                                   ! IN rho*entrainment rate implied by
!                                   !     placing of subsidence
    zrzi_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),    &
                                    ! IN (z-z_base)/(z_i-z_base)
    t_frac_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),  &
!                                   ! IN a fraction of the timestep
    z_half(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,         &
           bl_levels),                                                  &
                                    ! IN Z_HALF(*,K) is height of half
!                                   !    level k-1/2.
    zhsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                    ! IN Top of decoupled layer
    wt_ext(land_points,sm_levels) ! IN cumulative fract of trans'n

INTEGER, INTENT(IN) ::                                                  &
  kent(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                                  ! IN grid-level of SML inversion
  kent_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                  ! IN grid-level of DSC inversion

! Mineral dust source flux for tracer mixing
REAL, INTENT(IN) ::                                                     &
  dust_flux(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,        &
            ndiv),                                                      &
  dust_emiss_frac(land_points,ntiles),                                  &
                              ! IN fraction of tile can emit dust
  u_s_t_tile(land_points,ntiles,ndivh),                                 &
                                         !IN threshold frict. vel
  u_s_t_dry_tile(land_points,ntiles,ndivh),                             &
                                             !IN dry soil value
  u_s_std_tile(land_points,ntiles)!IN friction velocity

! IN additional variables for JULES
INTEGER, INTENT(IN) ::                                                  &
 tile_pts(ntype),                                                       &
                               ! IN Number of tile points.
 tile_index(land_points,ntype)
!                                ! IN Index of tile points.

REAL, INTENT(IN) ::                                                     &
 tile_frac(land_points,ntiles),                                         &
                              ! IN fractional coverage for each
                              !    surface tile
 canopy(land_points,ntiles),                                            &
                                 ! IN Surface/canopy water (kg/m2)
 alpha1(land_points,ntiles),                                            &
                                ! IN Mean gradient of saturated
!                                 specific humidity with
!                                 respect to temperature between
!                                 the bottom model layer and the
!                                 tile surfaces.
   fraca(land_points,ntiles),                                           &
                                   ! IN Fraction of surface
                                !            moisture flux with only
                                !            aerodynamic resistance.
   rhokh_tile(land_points,ntiles),                                      &
                                      ! IN
!                                 Tile surface exchange coefficients
!                                 for heat
   smc(land_points),                                                    &
                                ! IN Soil moisture content in root depth
!                                  (kg/m2).
   chr1p5m(land_points,ntiles),                                         &
                                   ! IN Ratio of coefficients reqd for
!                                 calculation of 1.5 m T.
   resfs(land_points,ntiles),                                           &
                                 ! IN Combined soil, stomatal
!                                 and aerodynamicresistance
!                                 factor = PSIS/(1+RS/RA) for
!                                 fraction (1-FRACA)
   z0hssi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
   z0mssi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                                ! IN Roughness lengths over sea (m).
   z0m_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                                ! IN Gridbox mean roughness length
!                               !    for momentum (m).
   canhc_tile(land_points,ntiles),                                      &
!                               ! IN Areal heat capacity of canopy
!                               !    for land tiles (J/K/m2).
   flake(land_points,ntiles),                                           &
                                   ! IN Lake fraction.
   wt_ext_tile(land_points,sm_levels,ntiles),                           &
!                               ! IN Fraction of evapotranspiration
!                               !    which is extracted from each
!                               !    soil layer by each tile.
   lw_down(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                ! IN Surface downward LW radiation
!                               !    (W/m2).
   lai_ft(land_points,npft),                                            &
                                   ! IN LAI on vegetated tiles
   canht_ft(land_points,npft),                                          &
                                   ! IN CANHT on vegetated tiles
   sw_tile(land_points,ntiles),                                         &
                                   ! IN Surface net SW radiation on land
!                               !    tiles (W/m2).
   ashtf_tile(land_points,ntiles),                                      &
!                               ! IN Coefficient to calculate
!                               !    surface heat flux into land
!                               !    tiles.
   resft(land_points,ntiles),                                           &
                                   ! IN Total resistance factor.
!                               !    FRACA+(1-FRACA)*RESFS for
!                               !    snow-free land, 1 for snow.
   rhokh_sice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,      &
                                                             nice_use), &
!                               ! IN Surface exchange coefficients
!                               !    for sea-ice
   rhokh_sea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
!                               ! IN Surface exchange coefficients
!                               !    for sea
   z0h_tile(land_points,ntiles),                                        &
                                ! IN Tile roughness lengths for heat
!                               !    and moisture (m).
   z0m_tile(land_points,ntiles),                                        &
                                ! IN Tile roughness lengths for
!                               !    momentum.
   chr1p5m_sice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),   &
!                               ! IN CHR1P5M for sea and sea-ice
!                               !    (leads ignored).
   fland(land_points),                                                  &
                                ! IN Land fraction on land tiles.
   flandg(pdims_s%i_start:pdims_s%i_end,                                &
          pdims_s%j_start:pdims_s%j_end),                               &
                                ! IN Land fraction on all points.
   flandg_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),       &
!                               ! IN Land frac (on U-grid, with 1st
!                               !    and last rows undefined or, at
!                               !    present, set to "missing data")
   flandg_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),       &
!                               ! IN Land frac (on V-grid, with 1st
!                               !    and last rows undefined or, at
!                               !    present, set to "missing data")
    vshr_land(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
!                               ! IN VSHR over land part of gridbox.
    vshr_ssi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                ! IN VSHR over sea part of gridbox.
    gc(land_points,ntiles),                                             &
                                  !Stomatal conductance to evapn
!                                 !    for land tiles (m/s).
    aresist_tile(land_points,ntiles)
!                                  !1/(CD_STD*VSHR) on land tiles
!                                  !for CLASSIC aerosol scheme

! Additional variables for JULES
REAL, INTENT(IN) ::                                                     &
 emis_tile(land_points,ntiles),                                         &
                                 ! IN  Emissivity for land tiles
  t_soil(land_points,sm_levels),                                        &
                                     ! slt/qrclim.slt_pm(lev).(month)
 snow_tile(land_points,ntiles),                                         &
!                               ! INOUT Snow on tiles (kg/m2).
   rib_ssi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                   ! INOUT Sea mean bulk Richardson
!                                          number for lowest layer.

! IN JULES additional STASH variables
REAL, INTENT(IN) ::                                                     &
 Gs(land_points),                                                       &
                              ! IN "Stomatal" conductance to
!                               !    evaporation (m/s).
   gpp(land_points),                                                    &
                                ! IN Gross primary productivity
!                               !    (kg C/m2/s).
   npp(land_points),                                                    &
                                ! IN Net primary productivity
!                               !    (kg C/m2/s).
   resp_p(land_points),                                                 &
                                ! IN Plant respiration (kg C/m2/s).
   gpp_ft(land_points,npft),                                            &
                                ! IN Gross primary productivity
!                               !    on PFTs (kg C/m2/s).
   npp_ft(land_points,npft),                                            &
                                ! IN Net primary productivity
!                               !    on PFTs (kg C/m2/s).
   resp_p_ft(land_points,npft),                                         &
                                  !IN Plant respiration on PFTs
!                               !     (kg C/m2/s).
   resp_s(land_points,dim_cs1),                                         &
                                   ! IN Soil respiration (kg C/m2/s).
   resp_s_tot(dim_cs2),                                                 &
                                   ! IN Total soil resp'n (kg C/m2/s).
   cs(land_points,dim_cs1),                                             &
                                   ! IN Soil carbon
   rib_tile(land_points,ntiles),                                        &
!                               ! IN RIB for land tiles.
   fsmc(land_points,npft),                                              &
                                ! IN Moisture availability factor.
   catch(land_points,ntiles),                                           &
                                ! IN Surface/canopy water capacity
!                               !    of snow-free land tiles (kg/m2).
   g_leaf(land_points,npft),                                            &
                                ! IN Leaf turnover rate (/360days).
   co2_emits(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                 !IN CO2 Emissions
   co2flux(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                      ! IN ocean CO2 flux

! arguments with intent INOUT. ie: input variables changed on output.
!     Declaration of new BL diagnostics.
TYPE (strnewbldiag), INTENT(INOUT) :: BL_diag
TYPE (strnewsfdiag), INTENT(INOUT) :: sf_diag

! Diagnostics info
REAL, INTENT(INOUT) ::                                                  &
stashwork3(*),                                                          &
                   ! STASH workspace for section 3 (Boundary Layer)
stashwork9(*) ! STASH workspace for section 9 (LS Cloud)

! Additional variables for SCM diagnostics which are dummy in full UM
INTEGER, INTENT(IN) ::                                                  &
  nSCMDpkgs              ! No of SCM diagnostics packages

LOGICAL, INTENT(IN) ::                                                  &
  L_SCMDiags(nSCMDpkgs)  ! Logicals for SCM diagnostics packages

REAL, INTENT(INOUT) :: TScrnDcl_SSI(pdims%i_start:pdims%i_end,          &
                                    pdims%j_start:pdims%j_end)
!                           !    Decoupled screen-level temperature
!                           !    over sea or sea-ice
REAL, INTENT(INOUT) :: TScrnDcl_TILE(land_points,ntiles)
!                           !    Decoupled screen-level temperature
!                           !    over land tiles
REAL, INTENT(INOUT) :: tStbTrans(pdims%i_start:pdims%i_end,             &
                                    pdims%j_start:pdims%j_end)
!                           !    Time since the transition

REAL, INTENT(INOUT) ::                                                  &
 cca0   (tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         n_cca_levels),                                                 &
 ccw0   (tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         tdims%k_end),                                                  &
 cca0_2d(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)

REAL, INTENT(INOUT) ::                                                  &
  resist_b(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                              ! INOUT (1/CH-1/(CD_STD)/VSHR
                              !     for CLASSIC aerosol scheme
    resist_b_tile(land_points,ntiles),                                  &
!                                  !(1/CH-1/CD_STD)/VSHR on land tiles
!                                  !for CLASSIC aerosol scheme
   epot_tile(land_points,ntiles),                                       &
!                               ! INOUT surface tile potential
!                               !       evaporation
  fqw(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,              &
      bl_levels),                                                       &
                                        ! needed as diagnostic ?
  ftl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,              &
      bl_levels),                                                       &
                                        ! needed as diagnostic
  taux(udims%i_start:udims%i_end,udims%j_start:udims%j_end,             &
                         bl_levels),                                    &
                                         ! needed as diagnostic
  tauy(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,             &
                           bl_levels),                                  &
                                           ! needed as diagnostic
  rhokh (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,           &
         bl_levels),                                                    &
  fqw_ice(pdims%i_start:pdims%i_end,                                    &
          pdims%j_start:pdims%j_end,nice_use),                          &
                              ! INOUT Surface FQW for sea-ice
  ftl_ice(pdims%i_start:pdims%i_end,                                    &
          pdims%j_start:pdims%j_end,nice_use),                          &
                              ! INOUT Surface FTL for sea-ice
  dtstar_tile(land_points,ntiles),                                      &
                              ! INOUT Change in TSTAR over timestep
                              !       for land tiles
  dtstar_sea(pdims%i_start:pdims%i_end,                                 &
             pdims%j_start:pdims%j_end),                                &
                              ! INOUT Change is TSTAR over timestep
                              !       for open sea
  dtstar_sice(pdims%i_start:pdims%i_end,                                &
              pdims%j_start:pdims%j_end,nice_use),                      &
                              ! INOUT Change is TSTAR over timestep
                              !       for sea-ice
    ti(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice),       &
                                       ! category sea ice sfc layer temp
                                       ! (IN only if l_sice_multilayers=T)
    r_u(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end,    &
          udims_s%k_start:udims_s%k_end),                               &
    r_v(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end,    &
          vdims_s%k_start:vdims_s%k_end),                               &
    r_w(wdims%i_start:wdims%i_end,wdims%j_start:wdims%j_end,            &
        wdims%k_start:wdims%k_end),                                     &
    t_latest(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             tdims%k_end),                                              &
    q_latest(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             tdims%k_end),                                              &
    qcl_latest(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,     &
               tdims%k_end),                                            &
    qcf_latest(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,     &
               tdims%k_end),                                            &
    cf_latest(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,      &
              tdims%k_end),                                             &
    cfl_latest(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,     &
               tdims%k_end),                                            &
    cff_latest(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,     &
               tdims%k_end),                                            &
    cf_area(tdims%i_start:tdims%i_end,                                  &
            tdims%j_start:tdims%j_end,tdims%k_end),                     &
    cf_bulk(tdims%i_start:tdims%i_end,                                  &
            tdims%j_start:tdims%j_end,tdims%k_end),                     &
    cf_liquid(tdims%i_start:tdims%i_end,                                &
              tdims%j_start:tdims%j_end,tdims%k_end),                   &
    cf_frozen(tdims%i_start:tdims%i_end,                                &
              tdims%j_start:tdims%j_end,tdims%k_end),                   &
    sum_eng_fluxes(pdims%i_start:pdims%i_end,                           &
                   pdims%j_start:pdims%j_end),                          &
    sum_moist_flux(pdims%i_start:pdims%i_end,                           &
                   pdims%j_start:pdims%j_end),                          &
    rhcpt(rhc_row_length, rhc_rows, tdims%k_end)

! Tracer variables
REAL, INTENT(INOUT) ::                                                  &
  aerosol     (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end),                          &
  free_tracers(tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end,                           &
               tr_vars)

REAL, INTENT(INOUT) ::                                                  &
  dust_div1   (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end ),                         &
  dust_div2   (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end ),                         &
  dust_div3   (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end ),                         &
  dust_div4   (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end ),                         &
  dust_div5   (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end ),                         &
  dust_div6   (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end ),                         &
  drydep2(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          ndiv) !dry dep though grav. set.

REAL, INTENT(INOUT) ::                                                  &
  so2         (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end ),                         &
  dms         (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end ),                         &
  so4_aitken  (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end ),                         &
  so4_accu    (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end ),                         &
  so4_diss    (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end ),                         &
  nh3         (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end )

REAL, INTENT(INOUT) ::                                                  &
  soot_new    (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end ),                         &
  soot_aged   (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end ),                         &
  soot_cld    (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end ),                         &
  bmass_new   (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end ),                         &
  bmass_agd   (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end ),                         &
  bmass_cld   (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end ),                         &
  ocff_new    (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end ),                         &
  ocff_aged   (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end ),                         &
  ocff_cld    (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end ),                         &
  nitr_acc    (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end ),                         &
  nitr_diss   (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end ),                         &
  co2         (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end ),                         &
  ozone_tracer(tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end )

! INOUT additional variables for JULES
REAL, INTENT(INOUT) ::                                                  &
 tstar_tile(land_points,ntiles),                                        &
                              ! INOUT Surface tile temperature
 fqw_tile(land_points,ntiles),                                          &
                              ! INOUT surface tile moisture flux
 ftl_tile(land_points,ntiles),                                          &
!                               ! INOUT surface tile heat flux
   radnet_sea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
!                               ! INOUT Open sea surface net radiation.
   radnet_sice(pdims%i_start:pdims%i_end,                               &
               pdims%j_start:pdims%j_end,nice_use),                     &
!                               ! INOUT Sea-ice surface net radiation.
   olr(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                                ! IN    TOA - surface upward LW on
!                               !       last radiation timestep
!                               ! OUT   Corrected TOA outward LW
   tstar_sice_cat(tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end,nice_use),                  &
                                   ! INOUT Sea-ice sfc temperature (K).
   tstar_ssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),      &
                                   ! INOUT Sea mean sfc temperature (K).
   tstar_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),      &
                                   ! INOUT Open sea sfc temperature (K).
   taux_land(udims%i_start:udims%i_end,udims%j_start:udims%j_end),      &
                                   ! INOUT W'ly compt of land sfc wind
!                                  !       stress (N/sq m). (On U-grid
!                                  !       with first and last rows
!                                  !       undefined or, at present,
!                                  !       set to missing data
   taux_ssi(udims%i_start:udims%i_end,udims%j_start:udims%j_end),       &
                                   ! INOUT W'ly compt of sea sfc wind
!                                  !       stress (N/sq m). (On U-grid
!                                  !       with first and last rows
!                                  !       undefined or, at present,
!                                  !       set to missing data
   tauy_land(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),      &
                                   ! INOUT S'ly compt of land sfc wind
!                                  !       stress (N/sq m).  On V-grid;
!                                  !       comments as per TAUX.
   tauy_ssi(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)
                                   ! INOUT S'ly compt of sea sfc wind
!                                  !       stress (N/sq m).  On V-grid;
!                                  !       comments as per TAUX.

INTEGER, INTENT(INOUT) ::                                               &
  error_code

! height of lcl in a well-mixed BL (types 3 or 4), or cumulus (5 or 6)
! 0 otherwise
REAL, INTENT(INOUT) :: zlcl_mixed(pdims%i_start:pdims%i_end,            &
                                pdims%j_start:pdims%j_end)

! OUT surface flux potential temp and moisture scales (K, kg/kg)
! for use in stochastic forcing in the boundary layer
! (applied in stochastic_physics/bl_pert_theta.F90).
! qv_star_surf is only used if i_pert_theta == pert_theta_and_moist
REAL, INTENT(OUT) ::                                                    &
 theta_star_surf(tdims%i_start:tdims%i_end,                             &
                 tdims%j_start:tdims%j_end),                            &
 qv_star_surf(tdims%i_start:tdims%i_end,                                &
                 tdims%j_start:tdims%j_end)

! OUT additional variables for JULES
REAL, INTENT(OUT) ::                                                    &
  tstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),           &
  ti_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),           &
                                             !output seaice temp.
                                             ! sfc layer (ice mean)
  ext(land_points,sm_levels),                                           &
                                  ! Extraction of water from each
!                                    soil layer (kg/m2/s).
    snowmelt(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                        !output from sf_evap.
   tstar_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),     &
                                   ! INOUT Land mean sfc temperature (K)
   tstar_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),     &
                                   ! OUT Sea-ice sfc temperature (K).
                                   !     (Ice mean over categories)
   ei_tile(land_points,ntiles),                                         &
                                   ! OUT EI for land tiles
   ecan_tile(land_points,ntiles),                                       &
                                ! OUT ECAN for land tiles
   melt_tile(land_points,ntiles),                                       &
!                               ! OUT Snowmelt on tiles (kg/m2/s).
   surf_htf_tile(land_points,ntiles),                                   &
!                               ! OUT Net downward surface heat flux
!                               !     on tiles (W/m2)
    surf_ht_flux_land(pdims%i_start:pdims%i_end,                        &
                      pdims%j_start:pdims%j_end)

!---------------------------------------------------------------------
! Local variables
REAL, PARAMETER :: qcl_forced_min = 0.00005
                 ! minimum water content in forced cu clouds

REAL, PARAMETER :: qcl_max_factor = 0.1
                 ! maximum fraction of water vapour forced cu is
                 ! allowed to condense

REAL, PARAMETER :: cf_top = 0.1
                 ! forced cloud fraction at cloud top

REAL :: alpha_tr(bl_levels)  ! Implicit weights for tracers
!                                  !  = alpha_cd from RUN_BL namelist
!                                  ! or = 1 if TRWEIGHTS1 = ON

REAL ::                                                                 &
  denom,                                                                &
                 ! Denominator in PC2 inhomogeneous ice forcing calc.
  q4,                                                                   &
                 ! QCF increment in PC2 inhomog.    ice forcing calc.
  incloudice,                                                           &
                 ! Value of QCF/CFF used for calculating sink of CFF
  deltacff,                                                             &
                 ! Ice cloud fraction increment
  surf_ht_flux_gb(pdims%i_start:pdims%i_end,                            &
                  pdims%j_start:pdims%j_end),                           &
 co2_flux_tot(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                                    !  total CO2 flux
 land_co2(land_points)          !  terrestrial CO2 flux

REAL ::                                                                 &
  e_sea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                               ! needed as diagnostic ?
  h_sea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                               ! needed as diagnostic ?
 radnet_tile(land_points,ntiles),                                       &
!                               ! Tile surface net radiation.
   le_tile(land_points,ntiles),                                         &
                                   ! Surface latent heat flux for
!                               !       land tiles (W/m2).
    e_ssi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                                    ! Evaporation from mean sea
    ei_sice(pdims%i_start:pdims%i_end,                                  &
            pdims%j_start:pdims%j_end,nice_use),                        &
                                    ! Output from sf_evap.
    surf_ht_flux_sice(pdims%i_start:pdims%i_end,                        &
                      pdims%j_start:pdims%j_end,nice),                  &
                                    ! Category sea ice surface heat flux
    ftl_ssi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                                    ! Surface FTL for mean sea
    ecan(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                        !output from sf_evap.
    ei(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                                        !output from sf_evap.
   esoil_tile(land_points,ntiles),                                      &
                                ! Evaporation from bare soil (kg/m2)
   es(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &
                                ! Surface evapotranspiration from
!                               !     soil moisture store (kg/m2/s).
   dust_l1wrk(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                    ! mass of the bottom model level * timestep

!local variables for mineral dust
REAL, ALLOCATABLE :: dust_all(:,:,:,:) !dust mmr

! Arrays for mixing w
REAL, ALLOCATABLE :: w_int(:,:,:)      ! w on interior points
REAL, ALLOCATABLE :: rhokm_mix(:,:,:)  ! rhokm interpolated to rho-levels

REAL :: weight1, weight2, weight3  ! Work variables for interpolation of rhokm

REAL :: cca_at_base(tdims%i_start:tdims%i_end,                          &
                    tdims%j_start:tdims%j_end)

! loop counters
INTEGER ::                                                              &
  i, j, k, iScm, jScm,                                                  &
  kinvert,                                                              &
                              ! vertical index for inverted arrays.
  idiv

INTEGER :: i_field  ! counter for swap_bounds_mv

! Diagnostics controlled by Diagnostic switches
REAL ::                                                                 &
  rho1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                                             ! Density at level 1
  sea_ice_htf(pdims%i_start:pdims%i_end,                                &
              pdims%j_start:pdims%j_end, nice)
                                             !output seaice fcondtop
                                             !(downwards conductive flux
                                             ! used to force ice model)

REAL ::                                                                 &
  rhokh_mix (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,       &
             bl_levels)

! local variables
INTEGER ::                                                              &
  nclds      ! Number of radiation cloud levels ( <=  model levels)

!-------Needed for area_cloud formulation-------------------------------
INTEGER ::                                                              &
  levels_per_level,    & ! 3 is hardwired inside ls_arcld
  large_levels           ! depends on above and model levels

!-----------------------------------------------------------------------

REAL ::                                                                 &
 qcl_forced,                                                            &
 cf_forced,                                                             &
            ! forced cloud water content and fraction
 dqcl,                                                                  &
 dcfl,                                                                  &
            ! forced cloud water content and fraction increments
 qcl_tol,                                                               &
            ! max tolerated forced cloud water content
 cf_base,                                                               &
            ! forced cloud fraction at cloud base
 zc_depth
            ! forced cloud depth

! Local data arrays
REAL ::                                                                 &
  t(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tdims%k_end)

REAL, TARGET :: ext_ice_frac(0:rhc_row_length+1,0:rhc_rows+1)

REAL ::                                                                 &
  u_inc_bl(udims_s%i_start:udims_s%i_end,                               &
            udims_s%j_start:udims_s%j_end,bl_levels),                   &
  v_inc_bl(vdims_s%i_start:vdims_s%i_end,                               &
            vdims_s%j_start:vdims_s%j_end,bl_levels)

REAL ::                                                                 &
  interp(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
  temp3,                                                                &
  ! Workspace variables
  drydep_str(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

REAL, TARGET ::                                                         &
  ext_p_layer_centres(0:rhc_row_length+1,0:rhc_rows+1,                  &
                                         0:tdims%k_end),                &
  ext_tl(0:rhc_row_length+1, 0:rhc_rows+1,tdims%k_end),                 &
  ext_ql(0:rhc_row_length+1, 0:rhc_rows+1,tdims%k_end),                 &
  ext_qcf(0:rhc_row_length+1,0:rhc_rows+1,tdims%k_end)

REAL ::                                                                 &
 work2d_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),         &
                                  ! Single-level work array (cloud)
 plsp(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                  ! Probability of large-scale precip
REAL, ALLOCATABLE :: f3_at_p(:, :)
                                  ! Coriolis parameter at theta-points

REAL :: rho_theta_levs                  ! Rho averaged onto theta levels
REAL :: lwp         (ScmRowLen, ScmRow) ! Liquid water path (kg/m2)
REAL :: iwp         (ScmRowLen, ScmRow) ! Ice water path (kg/m2)
REAL :: bl_alltypes (ScmRowLen, ScmRow) ! Boundary layer all types

CHARACTER (LEN=*), PARAMETER :: RoutineName = 'NI_IMP_CTL'

! Allocatable arrays for diagnostic variables - required to save memory
! when diagnostic not requested
REAL,ALLOCATABLE::                                                      &
 combined_cloud(:,:,:),   & ! Mixed CCA and CF per gridbox
 t_earliest(:,:,:),                                                     &
 q_earliest(:,:,:),                                                     &
 qcl_earliest(:,:,:),                                                   &
 qcf_earliest(:,:,:),                                                   &
 cf_earliest(:,:,:),                                                    &
 cfl_earliest(:,:,:),                                                   &
 cff_earliest(:,:,:),                                                   &
 t_inc_pc2(:,:,:),        & !  temperature     increment due to PC2 homog
 q_inc_pc2(:,:,:),        & !  humidity        increment due to PC2 homog
 qcl_inc_pc2(:,:,:),      & !  qCL             increment due to PC2 homog
 cfl_inc_pc2(:,:,:),      & !  cf_liquid       increment due to PC2 homog
 bcf_inc_pc2(:,:,:)         !  bulk cloud      increment due to PC2 homog

! Variables required for calculating implicit surface stress, ustar_imp
REAL ::                                                                 &
 fb_surf,           & ! surface buoyancy flux
 wm                   ! turbulent velocity scale

REAL :: casqcl
                      ! Temporary qcl used in CASIM.

REAL, PARAMETER :: c_ws=0.25
                      ! constant in turbulent velocity scale
REAL, ALLOCATABLE, TARGET :: taux_halo(:,:)
REAL, ALLOCATABLE, TARGET :: tauy_halo(:,:)
                      ! surface stresses including halo points
REAL, ALLOCATABLE :: tmp_p(:,:)
                      ! temporary in horizontal interpolation
REAL, ALLOCATABLE :: ustar_imp(:,:)
                      ! implicit calculation of ustar

REAL, ALLOCATABLE ::                                                    &
 zeros(:,:,:),             & ! Array of zero values
 tl_force(:,:,:),          & ! Forcing of TL by homogenous processes
 qt_force(:,:,:),          & ! Forcing of QT by homogenous processes
 ccw_cca(:,:,:),           & ! Convective cloud water * frac (i.e. gridbox
                             ! mean ccw)
 cca_3d (:,:,:)              ! 3D array of convective cloud frac

REAL, ALLOCATABLE ::                                                    &
  cca4comb_cld(:,:,:), &! Used to calculate combined cloud
  ccb4comb_cld(:,:),   &! Used to calculate combined cloud
  cct4comb_cld(:,:)     ! Used to calculate combined cloud

! STASHflag switches for increment diagnostics:
LOGICAL ::                                                              &
 l_apply_diag
                           ! flag to determine when to apply
                           ! diagnostics when iterating

! Switches for field calculations to support STASH diagnostics
LOGICAL ::                                                              &
 l_combi_cld,         & ! combined cloud amount
 l_wind_gust,         & ! wind gust diagnostic, 3463
 l_thermal,           & ! BL thermal speed, 3355
 l_theta_star,        & ! theta star, 3487
 l_qv_star,           & ! qv star, 3488
 l_ustar_imp            ! implicit ustar, 3500

LOGICAL :: l_made_rho1 = .FALSE. ! rho1 has been calculated, used for many diags

TYPE(swapable_field_pointer_type) :: fields_to_swap(4) ! mv swapbounds

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!- End of Header
! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF ( error_code  ==  0) THEN
  ! ----------------------------------------------------------------------
  ! Section BL.0 Initialisation of variables.
  ! ----------------------------------------------------------------------
  IF ( trweights1  ==  on ) THEN
    !         ! Set all implict weights used by tracers to one
    !         !  - overweighting not necessary since tracers
    !         !    have no feedback on the diffusion coefficients
    DO k = 1, bl_levels
      alpha_tr(k) = 1.0
    END DO
  ELSE
    !        ! Set implict weights used by tracers to those input
    DO k = 1, bl_levels
      alpha_tr(k) = alpha_cd(k)
    END DO
  END IF

  SELECT CASE (model_type)

  CASE DEFAULT

    ! Apply diags at last cycle only
    l_apply_diag = cycleno == numcycles

    ! Set diagnostic flags required for boundary layer diagnostics from
    ! STASHflags.
    !        !--------------------------------------------------------
    !        ! Note that an equivalent block of code exists in routine
    !        ! ni_bl_ctl, and needs to be kept consistent.
    !        !--------------------------------------------------------
    !        ! Windspeed (227, 230) and u, v at 10m on 'B' or 'C' grid
    sf_diag%su10    = ((sf(209,3) .OR. sf(225,3) .OR. sf(227,3) .OR.    &
               sf(230,3) .OR. sf(463,3)) .AND. l_apply_diag)
    sf_diag%sv10    = ((sf(210,3) .OR. sf(226,3) .OR. sf(227,3) .OR.    &
               sf(230,3) .OR. sf(463,3)) .AND. l_apply_diag)
    sf_diag%slh     =  sf(234,3) .AND. l_apply_diag
    sf_diag%l_lw_surft = (sf(383,3) .OR. sf(384,3)) .AND. l_apply_diag
    sf_diag%sq_t1p5 = (sf(236,3) .OR. sf(237,3) .OR. sf(245,3) .OR.     &
               sf(247,3) .OR. sf(248,3) .OR. sf(250,3) .OR.             &
               l_scrn    .OR. sf(341,3) .OR. sf(342,3) .OR.             &
               sf(253,3) .OR. sf(328,3) .OR. sf(329,3)) .AND.           &
               l_apply_diag
    sf_diag%sq1p5   = sf_diag%sq_t1p5 .AND. l_apply_diag
    sf_diag%st1p5   = ((sf_diag%sq_t1p5 .OR. sf(355,3) .OR. sf(463,3))  &
                           .AND. l_apply_diag)
    sf_diag%l_t10m  = sf(344,3) .AND. l_apply_diag
    sf_diag%l_q10m  = sf(345,3) .AND. l_apply_diag
    sf_diag%suv10m_n = (sf(368,3) .OR. sf(369,3) .OR. sf(370,3) .OR.    &
           sf(371,3) .OR. sf(365,3) .OR. sf(366,3) .OR.                 &
           sf(367,3) ) .AND. l_apply_diag

    ! Sea ice topmelt (single category or multi category)
    sf_diag%simlt   = (sf(235,3) .OR. sf(257,3) ) .AND. l_apply_diag
    sf_diag%smlt    =  sf(258,3) .AND. l_apply_diag
    BL_diag%l_ftl   = (fric_heating /= off) .OR.                        &
              (sf(216,3) .AND. l_apply_diag)
    BL_diag%l_fqw   = (fric_heating /= off) .OR.                        &
              (sf(222,3) .AND. l_apply_diag)
    BL_diag%l_taux  = (fric_heating /= off) .OR.                        &
              (( sf(219,3) .OR. sf(221,3)) .AND.                        &
              l_apply_diag)
    BL_diag%l_tauy  = (fric_heating /= off) .OR.                        &
              (( sf(220,3) .OR. sf(221,3)) .AND.                        &
              l_apply_diag)
    BL_diag%l_u_incr   =  sf(185,3) .AND. l_apply_diag
    BL_diag%l_v_incr   =  sf(186,3) .AND. l_apply_diag
    BL_diag%l_w_incr   =  sf(187,3) .AND. l_apply_diag
    BL_diag%l_t_incr   = (sf(181,9) .OR. sf(181,3))                     &
                         .AND. l_apply_diag
    BL_diag%l_tl_incr  =  sf(189,3) .AND. l_apply_diag
    BL_diag%l_q_incr   = (sf(182,9) .OR. sf(182,3) .OR.                 &
                         l_check_moist_inc) .AND. l_apply_diag
    BL_diag%l_qtl_incr = sf(190,3) .AND. l_apply_diag
    BL_diag%l_qcl_incr =(sf(183,9) .OR. sf(183,3) .OR. sf(170,3)        &
                         .OR. sf(171,3) .OR. l_check_moist_inc)         &
                         .AND. l_apply_diag
    BL_diag%l_qcf_incr = (sf(184,3) .OR. sf(172,3) .OR. sf(173,3)       &
                         .OR. l_check_moist_inc) .AND. l_apply_diag
    BL_diag%l_cf_incr = sf(192,3) .AND. l_apply_diag
    BL_diag%l_cfl_incr = (sf(193,3) .OR. sf(176,3) .OR. sf(177,3))      &
                         .AND. l_apply_diag
    BL_diag%l_cff_incr = (sf(194,3) .OR. sf(178,3) .OR. sf(179,3))      &
                         .AND. l_apply_diag

    ! New sea ice diagnostics for CMIP6
    sf_diag%l_tstar_sice_weighted_cat = (l_apply_diag .AND. sf(534,3))
    sf_diag%l_tstar_sice_weighted     = (l_apply_diag .AND. sf(535,3))
    sf_diag%l_lw_up_sice_weighted_cat = (l_apply_diag .AND. sf(530,3))
    sf_diag%l_lw_up_sice_weighted     = (l_apply_diag .AND. sf(531,3))
    sf_diag%l_ice_present_cat = (l_apply_diag .AND. sf(536,3))
    sf_diag%l_ice_present     = (l_apply_diag .AND. sf(537,3))
    sf_diag%l_ftl_ice_sm      = (l_apply_diag .AND. sf(533,3))

    ! Flag required for pre-calculation of cloud-related fields needed for
    ! cloud and vis diagnostics. Not duplicated in ni_bl_ctl.

    l_combi_cld = ( l_plsp .OR.                                         &
           sf(208,9) .OR. sf(209,9) .OR. sf(210,9) .OR. sf(211,9)       &
      .OR. sf(212,9) .OR. sf(213,9) .OR. sf(214,9) .OR. sf(215,9)       &
      .OR. sf(216,9) .OR. sf(217,9) .OR. sf(223,9) .OR. sf(231,9)       &
      .OR. sf(232,9) .OR. sf(233,9) .OR. sf(234,9)                      &
      ) .AND. l_apply_diag

    ! Diagnostics which depend on implicit ustar calculation
    l_wind_gust  = sf(463,3) .AND. l_apply_diag
    l_thermal    = sf(355,3) .AND. l_apply_diag
    l_theta_star = sf(487,3) .AND. l_apply_diag
    l_qv_star    = sf(488,3) .AND. l_apply_diag
    l_ustar_imp  = sf(500,3) .AND. l_apply_diag

  CASE (mt_single_column)
    ! Always apply diags in SCM if stash flags are on.
    l_apply_diag = .TRUE.

    ! diagnostic switches to be passed in
    ! These could be set using the input diagnostic requests to avoid
    ! un-necessary diagnostic calculations. Currently they are either
    ! on, if the diagnostic can be output, or
    ! off, if it can't.

    sf_diag%su10    = .TRUE.        ! Calculate u10m
    sf_diag%sv10    = .TRUE.        !    "      v10m
    sf_diag%slh     = .TRUE.        !    "      latent_heat
    sf_diag%l_lw_surft = .FALSE.    !    "      lw tile components
    sf_diag%sq1p5   = .TRUE.        !    "      q1p5m
    sf_diag%st1p5   = .TRUE.        !    "      t1p5m
    sf_diag%l_t10m  = .FALSE.
    sf_diag%l_q10m  = .FALSE.
    sf_diag%simlt   = .TRUE.        !    "      sice_mlt_htf
    sf_diag%smlt    = .TRUE.        !    "      snomlt_surf_htf
    sf_diag%sq_t1p5 = .FALSE.       ! This can cause diagnostics_bl to
                            ! change the values of q1p5m and t1p5m.

    !       Neutral winds are switched off since they are a more specialized
    !       diagnostic.
    sf_diag%suv10m_n = .FALSE.         ! Calculate neutral wind diagnostics

    BL_diag%l_ftl   = .TRUE.
    BL_diag%l_fqw   = .TRUE.
    BL_diag%l_taux  = .TRUE.
    BL_diag%l_tauy  = .TRUE.

    ! The following set true or diagnostics_bl will crash
    BL_diag%l_u_incr   = .TRUE.
    BL_diag%l_v_incr   = .TRUE.
    BL_diag%l_w_incr   = .FALSE.
    BL_diag%l_t_incr   = .FALSE.
    BL_diag%l_tl_incr  = .TRUE.
    BL_diag%l_q_incr   = .TRUE.
    BL_diag%l_qtl_incr = .TRUE.
    BL_diag%l_qcl_incr = .TRUE.
    BL_diag%l_qcf_incr = .TRUE.
    BL_diag%l_cf_incr  = .FALSE.
    BL_diag%l_cfl_incr = .FALSE.
    BL_diag%l_cff_incr = .FALSE.
    l_combi_cld        = l_plsp
    l_ustar_imp        = .TRUE.

    ! Diagnostics not currently available in SCM
    l_wind_gust  = .FALSE.
    l_thermal    = .FALSE.
    l_theta_star = .FALSE.
    l_qv_star    = .FALSE.

  END SELECT ! model_type

  IF (sf_diag%su10) THEN
    ALLOCATE(sf_diag%u10m(udims%i_start:udims%i_end,                    &
                          udims%j_start:udims%j_end))
  ELSE
    ALLOCATE(sf_diag%u10m(1,1))
  END IF

  IF (sf_diag%sv10) THEN
    ALLOCATE(sf_diag%v10m(vdims%i_start:vdims%i_end,                    &
                          vdims%j_start:vdims%j_end))
  ELSE
    ALLOCATE(sf_diag%v10m(1,1))
  END IF

  IF (sf_diag%slh) THEN
    ALLOCATE(sf_diag%latent_heat(pdims%i_start:pdims%i_end,             &
                                 pdims%j_start:pdims%j_end))
  ELSE
    ALLOCATE(sf_diag%latent_heat(1,1))
  END IF

  IF (sf_diag%l_lw_surft ) THEN
    ALLOCATE(sf_diag%lw_up_surft(land_points,ntiles))
    ALLOCATE(sf_diag%lw_down_surft(land_points,ntiles))
  ELSE
    ALLOCATE(sf_diag%lw_up_surft(1,1))
    ALLOCATE(sf_diag%lw_down_surft(1,1))
  END IF

  IF (sf_diag%st1p5 .OR. (IScrnTDiag == IP_ScrnDecpl2)                  &
      .OR. (IScrnTDiag == IP_ScrnDecpl3) ) THEN
    ALLOCATE(sf_diag%t1p5m(pdims%i_start:pdims%i_end,                   &
                           pdims%j_start:pdims%j_end))
    ALLOCATE(sf_diag%t1p5m_surft(land_points,ntiles))
  ELSE
    ALLOCATE(sf_diag%t1p5m(1,1))
    ALLOCATE(sf_diag%t1p5m_surft(1,1))
  END IF

  IF (sf_diag%sq1p5 .OR. (IScrnTDiag == IP_ScrnDecpl2)                  &
      .OR. (IScrnTDiag == IP_ScrnDecpl3) ) THEN
    ALLOCATE(sf_diag%q1p5m(pdims%i_start:pdims%i_end,                   &
                           pdims%j_start:pdims%j_end))
    ALLOCATE(sf_diag%q1p5m_surft(land_points,ntiles))
  ELSE
    ALLOCATE(sf_diag%q1p5m(1,1))
    ALLOCATE(sf_diag%q1p5m_surft(1,1))
  END IF
  ! 10m t and q diagnostics over sea/sea-ice
  ! both are required if one is requested due to ls_cld call in diag_bl
  IF (sf_diag%l_t10m .OR. sf_diag%l_q10m) THEN
    ALLOCATE( sf_diag%t10m(tdims%i_start:tdims%i_end,                   &
                           tdims%j_start:tdims%j_end))
    ALLOCATE( sf_diag%q10m(tdims%i_start:tdims%i_end,                   &
                           tdims%j_start:tdims%j_end))
    sf_diag%t10m = 0.0
    sf_diag%q10m = 0.0
  ELSE
    ALLOCATE( sf_diag%t10m(1,1))
    ALLOCATE( sf_diag%q10m(1,1))
  END IF

  IF (sf_diag%simlt) THEN
    ALLOCATE(sf_diag%sice_mlt_htf(pdims%i_start:pdims%i_end,            &
                                  pdims%j_start:pdims%j_end, nice))
  ELSE
    ALLOCATE(sf_diag%sice_mlt_htf(1,1,1))
  END IF

  IF (sf_diag%smlt) THEN
    ALLOCATE(sf_diag%snomlt_surf_htf(pdims%i_start:pdims%i_end,         &
                                     pdims%j_start:pdims%j_end))
  ELSE
    ALLOCATE(sf_diag%snomlt_surf_htf(1,1))
  END IF

  IF (sf_diag%l_tstar_sice_weighted_cat) THEN
    ALLOCATE(sf_diag%tstar_sice_weighted_cat(pdims%i_start:pdims%i_end, &
                                   pdims%j_start:pdims%j_end,nice_use))
  ELSE
    ALLOCATE(sf_diag%tstar_sice_weighted_cat(1,1,1))
  END IF

  IF (sf_diag%l_tstar_sice_weighted) THEN
    ALLOCATE(sf_diag%tstar_sice_weighted(pdims%i_start:pdims%i_end,     &
                                     pdims%j_start:pdims%j_end))
  ELSE
    ALLOCATE(sf_diag%tstar_sice_weighted(1,1))
  END IF

  IF (sf_diag%l_lw_up_sice_weighted_cat) THEN
    ALLOCATE(sf_diag%lw_up_sice_weighted_cat(pdims%i_start:pdims%i_end, &
                                   pdims%j_start:pdims%j_end,nice_use))
  ELSE
    ALLOCATE(sf_diag%lw_up_sice_weighted_cat(1,1,1))
  END IF

  IF (sf_diag%l_lw_up_sice_weighted) THEN
    ALLOCATE(sf_diag%lw_up_sice_weighted(pdims%i_start:pdims%i_end,     &
                                     pdims%j_start:pdims%j_end))
  ELSE
    ALLOCATE(sf_diag%lw_up_sice_weighted(1,1))
  END IF

  IF (sf_diag%l_ice_present_cat) THEN
    ALLOCATE(sf_diag%ice_present_cat(pdims%i_start:pdims%i_end,         &
                                   pdims%j_start:pdims%j_end,nice_use))
  ELSE
    ALLOCATE(sf_diag%ice_present_cat(1,1,1))
  END IF

  IF (sf_diag%l_ice_present) THEN
    ALLOCATE(sf_diag%ice_present(pdims%i_start:pdims%i_end,             &
                                     pdims%j_start:pdims%j_end))
  ELSE
    ALLOCATE(sf_diag%ice_present(1,1))
  END IF

  IF (sf_diag%l_ftl_ice_sm) THEN
    ALLOCATE(sf_diag%ftl_ice_sm(pdims%i_start:pdims%i_end,              &
                                     pdims%j_start:pdims%j_end))
  ELSE
    ALLOCATE(sf_diag%ftl_ice_sm(1,1))
  END IF

    ! Note this will be affected by the anvil scheme if it is
    ! applied, i.e the cca at the anvil base may have been
    ! scaled by the anvil tower_factor.
  IF (l_3d_cca) THEN
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        IF (lcbase(i,j) /= 0) THEN
          cca_at_base(i,j) = cca(i,j,lcbase(i,j))
        ELSE
          cca_at_base(i,j) = 0.0
        END IF
      END DO
    END DO
  ELSE
      ! CCA is only dimensioned with one level
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        IF (lcbase(i,j) /= 0) THEN
          cca_at_base(i,j) = cca(i,j,1)
        ELSE
          cca_at_base(i,j) = 0.0
        END IF
      END DO
    END DO
  END IF ! l_3d_cca

  !---------------------------------------------------------------------
  ! Intercept values of physics increments before being updated by
  ! implicit solver for optional output of bl wind increments
  !---------------------------------------------------------------------
  ! We need to store information about the increments of the temperature
  ! and moisture variables, so copy these to the _earliest variables.
  IF (BL_diag%l_t_incr .OR. BL_diag%l_tl_incr .OR. i_cld_vn==i_cld_pc2  &
       .OR. l_retain_bl_tendencies ) THEN

    ALLOCATE ( t_earliest(tdims%i_start:tdims%i_end,                    &
                          tdims%j_start:tdims%j_end,tdims%k_end) )

    ! Hold initial value of Temperature
!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                       &
!$OMP& SHARED( tdims, t_earliest, t_latest ) PRIVATE( i, j, k)
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          t_earliest(i,j,k) = t_latest(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
!$OMP END PARALLEL DO
  ELSE
    ALLOCATE( t_earliest(1,1,1) )
  END IF                   ! on STASHflags or PC2

  IF (BL_diag%l_q_incr .OR. BL_diag%l_qtl_incr                          &
     .OR. BL_diag%l_qcl_incr .OR. BL_diag%l_tl_incr                     &
     .OR. BL_diag%l_qcf_incr .OR. i_cld_vn==i_cld_pc2                   &
     .OR. l_retain_q_cl_bl_tendencies ) THEN

    ALLOCATE ( q_earliest(tdims%i_start:tdims%i_end,                    &
                          tdims%j_start:tdims%j_end,tdims%k_end) )
    ALLOCATE ( qcl_earliest(tdims%i_start:tdims%i_end,                  &
                            tdims%j_start:tdims%j_end,tdims%k_end) )
    ALLOCATE ( qcf_earliest(tdims%i_start:tdims%i_end,                  &
                            tdims%j_start:tdims%j_end,tdims%k_end) )
    ALLOCATE ( cf_earliest(tdims%i_start:tdims%i_end,                   &
                           tdims%j_start:tdims%j_end,tdims%k_end) )
    ALLOCATE ( cfl_earliest(tdims%i_start:tdims%i_end,                  &
                            tdims%j_start:tdims%j_end,tdims%k_end) )
    ALLOCATE ( cff_earliest(tdims%i_start:tdims%i_end,                  &
                            tdims%j_start:tdims%j_end,tdims%k_end) )

    ! Hold initial values of wet parameters
!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                       &
!$OMP& SHARED( tdims, q_earliest, q_latest, qcl_earliest, qcl_latest,   &
!$OMP&         qcf_earliest, qcf_latest, cf_earliest, cf_latest,        &
!$OMP&         cfl_earliest, cfl_latest, cff_earliest, cff_latest)      &
!$OMP& PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          q_earliest(i,j,k)   = q_latest(i,j,k)
          qcl_earliest(i,j,k) = qcl_latest(i,j,k)
          qcf_earliest(i,j,k) = qcf_latest(i,j,k)
          cf_earliest(i,j,k)  = cf_latest(i,j,k)
          cfl_earliest(i,j,k) = cfl_latest(i,j,k)
          cff_earliest(i,j,k) = cff_latest(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
!$OMP END PARALLEL DO

  ELSE
    ALLOCATE ( q_earliest(1,1,1) )
    ALLOCATE ( qcl_earliest(1,1,1) )
    ALLOCATE ( qcf_earliest(1,1,1) )
    ALLOCATE ( cf_earliest(1,1,1) )
    ALLOCATE ( cfl_earliest(1,1,1) )
    ALLOCATE ( cff_earliest(1,1,1) )

  END IF                  ! on STASHflags or PC2

  ! ----------------------------------------------------------------------
  ! Section BL.1 Calculate T at old time level.
  ! Modified to use latest values to avoid time-level inconsistencies
  ! with cloud data.
  ! ---------------------------------------------------------------------
!$OMP  PARALLEL DEFAULT(NONE)                                          &
!$OMP& SHARED( bl_levels, tdims, t, theta, exner_theta_levels, pdims,  &
!$OMP&         ecan, ei, snowmelt, IScrnTDiag, i_pert_theta,           &
!$OMP&         rho1, rho_wet_rsq, r_rho_levels ) PRIVATE( i, j, k )
!$OMP  DO SCHEDULE(STATIC)
  DO k = 1, bl_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        t(i,j,k) = theta(i,j,k) * exner_theta_levels(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

  !  Initialise output arrays to zero
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      ecan(i,j)=0.0
      ei(i,j) = 0.0
      snowmelt(i,j)=0.0
    END DO
  END DO
!$OMP END DO NOWAIT
  ! ----------------------------------------------------------------------
  ! Section BL.2  Call Implicit solver
  !               Call tracer mixing for qcf.
  !               Call tracer mixing for other tracers.
  ! ----------------------------------------------------------------------

  ! Density of the lowest level is required for the transient screen
  ! diagnostic and stochastic theta perturbations.
  IF ( IScrnTDiag == IP_ScrnDecpl2 .OR.                                 &
       IScrnTDiag == IP_ScrnDecpl3 .OR.                                 &
       i_pert_theta == pert_theta_star) THEN
!$OMP DO SCHEDULE(STATIC)
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        rho1(i,j) = rho_wet_rsq(i,j,1) /                                &
                    (r_rho_levels(i,j,1)*r_rho_levels(i,j,1))
      END DO
    END DO
!$OMP END DO
  END IF
!$OMP END PARALLEL
  ! make a note if we made rho1 (outside the OMP)
  IF ( IScrnTDiag == IP_ScrnDecpl2 .OR.                                 &
       IScrnTDiag == IP_ScrnDecpl3 .OR.                                 &
       i_pert_theta == pert_theta_star) THEN
    l_made_rho1 = .TRUE.
  END IF

  IF (IScrnTDiag == IP_ScrnDecpl2 .OR. IScrnTDiag == IP_ScrnDecpl3) THEN
    !         Allocate the Coriolis parameter on the p-grid and interpolate
    !         from the u-grid.
    ALLOCATE(f3_at_p(pdims%i_start:pdims%i_end,                         &
                     pdims%j_start:pdims%j_end))

    SELECT CASE (model_type)
    CASE DEFAULT

      CALL u_to_p (f3_at_u,                                             &
                    udims_s%i_start,udims_s%i_end,                      &
                    udims_s%j_start,udims_s%j_end,                      &
                    pdims%i_start,pdims%i_end,                          &
                    pdims%j_start,pdims%j_end,                          &
                    1, at_extremity, f3_at_p)

    CASE (mt_single_column)

      ! unsure what will happen with combination of endgame+scm
      f3_at_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)      &
      = f3_at_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end)

    END SELECT ! model_type

  ELSE
    ALLOCATE(f3_at_p(1,1))
  END IF

  !DEPENDS ON: imp_solver
  CALL imp_solver(ntiles, land_points,                                  &
  ! IN values defining model domain, vertical grid of model atmosphere
   bl_levels,alpha_cd,                                                  &
  ! IN U and V momentum fields and increments
   u, v, r_u, r_v,                                                      &
  ! IN soil/vegetation/land surface data :
   land_index, tile_frac,canopy,fland,                                  &
   flandg(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
  ! IN sea/sea-ice data :
   ice_fract,di_ncat,ice_fract_ncat,k_sice, u_0, v_0,                   &
  ! IN cloud data
   q, qcf, qcl, qcf_latest, qcl_latest, t,                              &
  ! IN everything not covered so far :
   rho_wet_tq, p_star,                                                  &
  ! IN variables required in IMP_SOLVER
    alpha1, ashtf, bq_gb, bt_gb,dtrdz_charney_grid, rdz_charney_grid,   &
    dtrdz_u, dtrdz_v, rdz_u, rdz_v, fraca, rhokh_tile, smc, chr1p5m,    &
    resfs,z0hssi,z0mssi,cdr10m_u,cdr10m_v,cdr10m_n_u,cdr10m_n_v,        &
    cd10m_n_u,cd10m_n_v,z_theta,zh,rhokm_u,rhokm_v,                     &
    k_blend_tq, k_blend_uv,                                             &
  ! IN variables for new BL solver
    bl_type_1,bl_type_2,                                                &
  ! IN additional variables for JULES
    tile_pts,tile_index,canhc_tile,flake,wt_ext_tile,lw_down,sw_tile,   &
    alpha1_sea,alpha1_sice,ashtf_sea,ashtf_tile,resft,rhokh_sice,       &
    rhokh_sea,z0h_tile,z0m_tile,chr1p5m_sice,flandg_u,flandg_v,         &
    co2(:, :, 1), rho1, f3_at_p, uStarGBM,                              &
    emis_tile,t_soil,snow_tile,                                         &
  ! INOUT diagnostics :-
    BL_diag, sf_diag,                                                   &
  ! INOUT data :
    dtstar_tile,dtstar_sea,dtstar_sice,ti,tstar_sice_cat,tstar_ssi,     &
    tstar_tile, tstar_sea, t_latest, q_latest,                          &
  ! INOUT  diagnostic started in bdy_layr not requiring stash flags :
    fqw_ice,ftl_ice,fqw,fqw_tile,epot_tile,ftl,ftl_tile,rhokh,          &
    taux,tauy,taux_land,taux_ssi,tauy_land,tauy_ssi,                    &
    TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans,                             &
  ! INOUT additional variables for JULES
    radnet_sice,olr,                                                    &
  ! OUT u and v increments
    u_inc_bl, v_inc_bl,                                                 &
  ! OUT  diagnostic not requiring stash flags :
    rhokh_mix, sea_ice_htf,                                             &
  ! OUT additional variables for JULES
    ti_gb, tstar,tstar_land,tstar_sice,e_sea,h_sea,le_tile,radnet_tile, &
    esoil_tile,surf_ht_flux_gb,surf_ht_flux_land,surf_ht_flux_sice,     &
    surf_htf_tile,ei_tile,ecan_tile,melt_tile,                          &
    e_ssi,ei_sice,ftl_ssi,error_code,                                   &
  ! OUT data required elsewhere in um system :
    ecan, ei, es, ext, snowmelt                                         &
    )

  !       Release space used for the screen diagnostic.
  DEALLOCATE(f3_at_p)

  ! Allocate arrays to store BL increments for physics_tendencies_mod
  IF (l_retain_bl_tendencies) THEN

    CALL init_bl_tendencies()
    ! Store U increment of BL
    DO k = 1, bl_levels
      DO j = udims%j_start, udims%j_end
        DO i = udims%i_start, udims%i_end
          du_bl(i,j,k) = u_inc_bl(i,j,k) - r_u(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
  END IF ! end if over l_retain_bl_tendencies

  ! add boundary layer increment to R_u and R_v
  IF ( BL_diag%l_u_incr ) THEN
    ! add boundary layer increment to R_u and R_v
    ALLOCATE ( BL_diag%u_incr(udims%i_start:udims%i_end,                &
                              udims%j_start:udims%j_end,                &
                              udims%k_start:udims%k_end) )
    BL_diag%u_incr = 0.0
    DO k = 1, bl_levels
      DO j = udims%j_start, udims%j_end
        DO i = udims%i_start, udims%i_end
          BL_diag%u_incr(i,j,k) = u_inc_bl(i,j,k) - r_u(i,j,k)
          r_u(i,j,k) = u_inc_bl(i,j,k)
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE( BL_diag%u_incr(1,1,1) )
!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                       &
!$OMP& SHARED( bl_levels,  udims, r_u, u_inc_bl ) PRIVATE( i, j, k)
    DO k = 1, bl_levels
      DO j = udims%j_start, udims%j_end
        DO i = udims%i_start, udims%i_end
          r_u(i,j,k) = u_inc_bl(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  ! Store V increment of BL
  IF (l_retain_bl_tendencies) THEN
    DO k = 1, bl_levels
      DO j = vdims%j_start, vdims%j_end
        DO i = vdims%i_start, vdims%i_end
          dv_bl(i,j,k) = v_inc_bl(i,j,k) - r_v(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
  END IF ! end if over l_retain_bl_tendencies

  IF ( BL_diag%l_v_incr ) THEN
    ALLOCATE ( BL_diag%v_incr(vdims%i_start:vdims%i_end,                &
                              vdims%j_start:vdims%j_end,                &
                              vdims%k_start:vdims%k_end) )
    BL_diag%v_incr = 0.0
    DO k = 1, bl_levels
      DO j = vdims%j_start, vdims%j_end
        DO i = vdims%i_start, vdims%i_end
          BL_diag%v_incr(i,j,k) = v_inc_bl(i,j,k) - r_v(i,j,k)
          r_v(i,j,k) = v_inc_bl(i,j,k)
        END DO
      END DO
    END DO
  ELSE
    ALLOCATE( BL_diag%v_incr(1,1,1) )
!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                       &
!$OMP& SHARED( bl_levels, vdims, r_v, v_inc_bl ) PRIVATE( i, j, k)
    DO k = 1, bl_levels
      DO j = vdims%j_start, vdims%j_end
        DO i = vdims%i_start, vdims%i_end
          r_v(i,j,k) = v_inc_bl(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  ! Pass in a zero field for source terms.
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      interp(i,j) = 0.0
    END DO
  END DO

  IF (l_bl_mix_qcf .AND. i_cld_vn > i_cld_off ) THEN
    ! Call tr_mix to mix qcf
    ! output qcf_flux in T to save workspace
    ! But only mix qcf in the vertical if we are not running
    ! with the CASIM boundary layer changes

    CALL  tr_mix (                                                      &
    ! IN fields
              bl_levels, alpha_tr, rhokh_mix(:,:,2:), rhokh_mix(:,:,1), &
              dtrdz_charney_grid,interp, interp,                        &
              kent, we_lim, t_frac, zrzi,                               &
              kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc,               &
              zhnl, zhsc, z_half,                                       &
    ! INOUT / OUT fields
              qcf_latest,t,drydep_str                                   &
              )

  END IF ! l_bl_mix_qcf .and. i_cld_vn is on

  !Call tr_mix to mix w in the vertical when the subgrid turbulence
  !scheme is activated

  IF (l_subfilter_vert) THEN

    ALLOCATE (w_int(wdims%i_start:wdims%i_end,                          &
                    wdims%j_start:wdims%j_end,bl_levels) )

    DO k = 1, bl_levels
      DO j = wdims%j_start, wdims%j_end
        DO i = wdims%i_start, wdims%i_end
          w_int(i,j,k) = w(i,j,k)
        END DO
      END DO
    END DO

    ALLOCATE( rhokm_mix(pdims%i_start:pdims%i_end,                      &
                        pdims%j_start:pdims%j_end, bl_levels)  )

    ! Interpolate rhokm from theta-levels to rho-levels, as is done for
    ! rhokh in bdy_expl2.
    ! NOTE: rhokm is defined on theta-levels with the k-indexing offset by
    ! 1 compared to the rest of the UM (k=1 is the surface).
!$OMP  PARALLEL DEFAULT(NONE)                                           &
!$OMP& SHARED(pdims, bl_levels, r_theta_levels, r_rho_levels,           &
!$OMP& rdz_charney_grid, rhokm_mix, rhokm )                             &
!$OMP& PRIVATE(i, j, k, weight1, weight2, weight3)

!$OMP DO SCHEDULE(STATIC)
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        ! Bottom model-level is surface in both arrays, so no interp needed
        ! (for rhokh_mix, this is done in the JULES routine sf_impl2_jls).
        rhokm_mix(i,j,1) = rhokm(i,j,1)
      END DO
    END DO
!$OMP END DO NOWAIT

    DO k = 2, bl_levels-1
!$OMP DO SCHEDULE(STATIC)
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          weight1 = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
          weight2 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
          weight3 = r_rho_levels(i,j,k)   - r_theta_levels(i,j,k-1)          
          rhokm_mix(i,j,k) = (weight3/weight1) * rhokm(i,j,k+1)         &
                           + (weight2/weight1) * rhokm(i,j,k)
          ! Scale exchange coefficients by 1/dz factor, as is done for
          ! rhokh_mix in bdy_impl4
          rhokm_mix(i,j,k) = rhokm_mix(i,j,k) * rdz_charney_grid(i,j,k)
          ! (note this doesn't need to be done for the surface exchange coef)
        END DO
      END DO
!$OMP END DO NOWAIT
    END DO

    k = bl_levels
!$OMP DO SCHEDULE(STATIC)
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        weight1 = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
        weight2 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
        ! assume rhokm(BL_LEVELS+1) is zero
        rhokm_mix(i,j,k) = (weight2/weight1) * rhokm(i,j,k)
        ! Scale exchange coefficients by 1/dz factor, as is done for
        ! rhokh_mix in bdy_impl4
        rhokm_mix(i,j,k) = rhokm_mix(i,j,k) * rdz_charney_grid(i,j,k)
        ! (note this doesn't need to be done for the surface exchange coef)
      END DO
    END DO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

    CALL  tr_mix (                                                      &
    ! In fields
              bl_levels, alpha_cd, rhokm_mix(:,:,2:), rhokm_mix(:,:,1), &
              dtrdz_charney_grid,interp, interp,                        &
              kent, we_lim, t_frac, zrzi,                               &
              kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc,               &
              zhnl, zhsc, z_half,                                       &
    ! INOUT / OUT fields
              w_int,t,drydep_str                                        &
              )

    ! add boundary layer increment to R_w
    IF ( BL_diag%l_w_incr ) THEN
      ALLOCATE ( BL_diag%w_incr(wdims%i_start:wdims%i_end,              &
                                wdims%j_start:wdims%j_end,              &
                                            1:wdims%k_end) )
      BL_diag%w_incr = 0.0
      DO k = 1, bl_levels
        DO j = wdims%j_start, wdims%j_end
          DO i = wdims%i_start, wdims%i_end
            BL_diag%w_incr(i,j,k) = w_int(i,j,k) - w(i,j,k)
          END DO
        END DO
      END DO
    ELSE
      ALLOCATE( BL_diag%w_incr(1,1,1) )
    END IF  ! test on diagnostic requested

    DO k = 1, bl_levels
      DO j = wdims%j_start, wdims%j_end
        DO i = wdims%i_start, wdims%i_end
          r_w(i,j,k) = r_w(i,j,k) + (w_int(i,j,k) - w(i,j,k))
        END DO
      END DO
    END DO

    DEALLOCATE (rhokm_mix)
    DEALLOCATE (w_int)

  END IF     !L_subfilter_vert

  ! apply BL tracer mixing and gravitational settling of dust
  ! on the last cycle only
  IF ( cycleno == numcycles ) THEN

    ! Gravitational settling of mineral dust

    IF (l_dust) THEN

      IF (flag_dustmmr_em) THEN
        ! need to produce a dust mmr at emission, for Visibility diagnostic.
        ! Estimate the dust mmr if the emission were just added,
        ! right now, into the bottom model level:
        k = 1
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            ! the thickness of the lowest dust-containing level is actually
            ! the height of the second rho level (so z_half(2))
            dust_l1wrk(i,j) = timestep / (rho_wet_tq(i,j,1) * z_half(i,j,2))
          END DO
        END DO
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            pws_dustmmr1_em(i,j) = dust_div1(i,j,k) +                      &
                                   (dust_flux(i,j,1) * dust_l1wrk(i,j))
            pws_dustmmr2_em(i,j) = dust_div2(i,j,k) +                      &
                                   (dust_flux(i,j,2) * dust_l1wrk(i,j))
          END DO
        END DO
        IF (.NOT. l_twobin_dust) THEN
          DO j = tdims%j_start, tdims%j_end
            DO i = tdims%i_start, tdims%i_end
              pws_dustmmr3_em(i,j) = dust_div3(i,j,k) +                      &
                                   (dust_flux(i,j,3) * dust_l1wrk(i,j))
              pws_dustmmr4_em(i,j) = dust_div4(i,j,k) +                     &
                                   (dust_flux(i,j,4) * dust_l1wrk(i,j))
              pws_dustmmr5_em(i,j) = dust_div5(i,j,k) +                      &
                                   (dust_flux(i,j,5) * dust_l1wrk(i,j))
              pws_dustmmr6_em(i,j) = dust_div6(i,j,k) +                      &
                               (dust_flux(i,j,6) * dust_l1wrk(i,j))
            END DO
          END DO
        END IF    
      END IF
      ! now proceed with actual dust work:

      ALLOCATE( dust_all(tdims%i_start:tdims%i_end,                     &
                         tdims%j_start:tdims%j_end,                     &
                         1:tdims%k_end,ndiv) )

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                       &
!$OMP& SHARED( tdims, t, theta, exner_theta_levels, dust_all, dust_div1,&
!$OMP&         dust_div2, l_twobin_dust, dust_div3, dust_div4,          &
!$OMP&         dust_div5, dust_div6 ) PRIVATE( i, j, k )
      DO k = 1, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            t(i,j,k)=theta(i,j,k)*exner_theta_levels(i,j,k)
            dust_all(i,j,k,1)=dust_div1(i,j,k)
            dust_all(i,j,k,2)=dust_div2(i,j,k)
            IF (.NOT. l_twobin_dust) THEN
              dust_all(i,j,k,3)=dust_div3(i,j,k)
              dust_all(i,j,k,4)=dust_div4(i,j,k)
              dust_all(i,j,k,5)=dust_div5(i,j,k)
              dust_all(i,j,k,6)=dust_div6(i,j,k)
            END IF
          END DO !i
        END DO !j
      END DO !k
!$OMP END PARALLEL DO

      DO idiv = 1, ndiv

        CALL gravsett(                                                  &
   drep(idiv),rhop,p_layer_centres,p_layer_boundaries,t,                &
   dust_all(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,        &
            pdims%k_start:pdims%k_end,idiv),                            &
   drydep2(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
           idiv))

      END DO

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                       &
!$OMP& SHARED( tdims, dust_div1, dust_all, dust_div2, l_twobin_dust,    &
!$OMP&         dust_div3, dust_div4, dust_div5, dust_div6 )             &
!$OMP& PRIVATE( i, j, k )
      DO k = 1, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            dust_div1(i,j,k)=dust_all(i,j,k,1)
            dust_div2(i,j,k)=dust_all(i,j,k,2)
            IF (.NOT. l_twobin_dust) THEN
              dust_div3(i,j,k)=dust_all(i,j,k,3)
              dust_div4(i,j,k)=dust_all(i,j,k,4)
              dust_div5(i,j,k)=dust_all(i,j,k,5)
              dust_div6(i,j,k)=dust_all(i,j,k,6)
            END IF
          END DO !ROW_LENGTH
        END DO !ROWS
      END DO !MODEL_LEVELS
!$OMP END PARALLEL DO
      DEALLOCATE (dust_all)
    END IF !L_DUST

    ! Mixing for all non-qcf tracers done in subroutine
    IF ( l_bl_tracer_mix .OR. l_sulpc_so2 .OR. l_soot .OR.              &
         l_co2_interactive .OR. l_murk .OR. l_biomass .OR.              &
         l_dust .OR. l_ocff .OR. l_nitrate) THEN

      CALL bl_trmix_dd(                                                 &
      ! IN arguments
                   bl_levels,dtrdz_charney_grid,tr_vars,                &
                   alpha_tr, rhokh_mix, p_star,                         &
      ! IN Emissions fields
                   dust_flux, co2_emits, co2flux, npp, resp_s_tot,      &
      ! IN variables for tr_mix
                   kent, we_lim, t_frac, zrzi,                          &
                   kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc,          &
                   zhnl, zhsc, z_half,                                  &
      ! IN for dry deposition of tracers
                  rho_aresist, aresist,                                 &
                   r_b_dust, tstar, land_points, land_index, ice_fract, &
      ! IN variables for sresfact
                  ntiles, tile_pts, tile_index, tile_frac,              &
                  canopy, catch, snow_tile, gc,                         &
                  aresist_tile, flandg,                                 &
      ! INOUT Fields to mix
                   aerosol, free_tracers, ozone_tracer,                 &
                   drydep_str, resist_b, resist_b_tile,                 &
      ! Mineral dust
                   dust_div1, dust_div2, dust_div3,                     &
                   dust_div4, dust_div5, dust_div6,                     &
      ! Sulphur cycle
                   so2, dms, so4_aitken, so4_accu, so4_diss, nh3,       &
      ! Soot cycle
                   soot_new, soot_aged, soot_cld,                       &
      ! Biomass aerosol
                   bmass_new, bmass_agd, bmass_cld,                     &
      ! Fossil-fuel organic carbon aerosol
                   ocff_new, ocff_aged, ocff_cld,                       &
      ! Ammonium nitrate aerosol
                   nitr_acc, nitr_diss,                                 &
      ! Carbon cycle
                   co2, co2_flux_tot, land_co2,                         &
                  stashwork3                                            &
                  )

    END IF  ! if tracer mixing required

  END IF ! ( CycleNo == NumCycles )

  ! ----------------------------------------------------------------------
  ! Calculate turbulent surface theta scale, theta_star_surf,
  ! turbulent moisture scale, q_star_surf,
  ! also calcualte implicit ustar for diagnostics if required
  ! ----------------------------------------------------------------------
  IF ( (i_pert_theta > pert_theta_mag) .OR. l_qv_star .OR.              &
       l_wind_gust .OR. l_thermal .OR. l_theta_star .OR. l_ustar_imp )  &
       THEN
#if !defined(LFRIC)

    ! x-stress:
    ALLOCATE(taux_halo(udims_s%i_start:udims_s%i_end,                   &
                       udims_s%j_start:udims_s%j_end))
    taux_halo(:,:)  = 0.0
    taux_halo(udims%i_start:udims%i_end,udims%j_start:udims%j_end)      &
           = taux(udims%i_start:udims%i_end,udims%j_start:udims%j_end,1)

    ! Set up for message passing
    i_field = 0
    i_field = i_field + 1
    fields_to_swap(i_field) % field_2d    => taux_halo(:,:)
    fields_to_swap(i_field) % field_type  =  fld_type_u
    fields_to_swap(i_field) % levels      =  1
    fields_to_swap(i_field) % rows        =  udims%j_end
    fields_to_swap(i_field) % vector      =  .TRUE.

    ! y-stress:
    ALLOCATE(tauy_halo(vdims_s%i_start:vdims_s%i_end,                   &
                       vdims_s%j_start:vdims_s%j_end))
    tauy_halo(:,:) = 0.0
    tauy_halo(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)      &
           = tauy(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,1)
    i_field = i_field + 1
    fields_to_swap(i_field) % field_2d    => tauy_halo(:,:)
    fields_to_swap(i_field) % field_type  =  fld_type_v
    fields_to_swap(i_field) % levels      =  1
    fields_to_swap(i_field) % rows        =  vdims%j_end-vdims%j_start+1
    fields_to_swap(i_field) % vector      =  .TRUE.

    CALL swap_bounds_2d_mv( fields_to_swap, i_field,                    &
                        pdims%i_end, pdims_s%halo_i, pdims_s%halo_j)

    ALLOCATE(tmp_p(pdims%i_start:pdims%i_end,                           &
                   pdims%j_start:pdims%j_end))
    ALLOCATE(ustar_imp(pdims%i_start:pdims%i_end,                       &
                       pdims%j_start:pdims%j_end))

    CALL u_to_p(taux_halo,                                              &
                udims_s%i_start,udims_s%i_end,                          &
                udims_s%j_start,udims_s%j_end,                          &
                pdims%i_start,pdims%i_end,                              &
                pdims%j_start,pdims%j_end,                              &
                1,at_extremity,tmp_p)

    ! Begin to form the magnitude of the stress.
    ustar_imp(:,:) = tmp_p(:,:) * tmp_p(:,:)

    CALL v_to_p(tauy_halo,                                              &
                vdims_s%i_start,vdims_s%i_end,                          &
                vdims_s%j_start,vdims_s%j_end,                          &
                pdims%i_start,pdims%i_end,                              &
                pdims%j_start,pdims%j_end,                              &
                1,at_extremity,tmp_p)

    ustar_imp(:,:) = ustar_imp(:,:) + tmp_p(:,:) * tmp_p(:,:)

    DEALLOCATE(taux_halo)
    DEALLOCATE(tauy_halo)
    DEALLOCATE(tmp_p)

    ustar_imp(:,:) = SQRT( SQRT(ustar_imp(:,:)) / rho1(:,:) )

    IF ( (i_pert_theta > pert_theta_mag) .OR.                           &
         l_theta_star .OR. l_qv_star) THEN

      DO j = pdims%j_start,pdims%j_end
        DO i = pdims%i_start,pdims%i_end
          fb_surf = g * ( bt_gb(i,j,1)*(ftl(i,j,1)/cp) +                &
                          bq_gb(i,j,1)*fqw(i,j,1) ) / rho1(i,j)
          IF (ftl(i,j,1) > 0.0 .AND. fb_surf > 0.0) THEN
            wm = ( ustar_imp(i,j)**3.0+c_ws*zh(i,j)*fb_surf )**one_third
            ! Eg: ftl=15 W/m2 and zh=500m gives wm~0.4 and
            !     not interested in very weakly forced convection
            wm = MAX( 0.4, wm )
            theta_star_surf(i,j)=ftl(i,j,1) / (cp*wm*rho1(i,j))
            IF ( (i_pert_theta == pert_theta_and_moist) .OR.            &
                 l_qv_star )                                            &
                 qv_star_surf(i,j)=MAX( 0.0, fqw(i,j,1) ) /             &
                                   (wm*rho1(i,j))
          ELSE
            theta_star_surf(i,j)=0.0
            IF ( (i_pert_theta == pert_theta_and_moist) .OR.            &
                  l_qv_star)                                            &
                 qv_star_surf(i,j)=0.0
          END IF
        END DO
      END DO

    END IF ! i_pert_theta or l_theta_star or l_qv_star
#else
    theta_star_surf=0.0
    qv_star_surf=0.0
#endif
  ELSE
    ALLOCATE(ustar_imp(1,1))
  END IF  ! i_pert_theta or l_wind_gust or l_thermal or l_theta_star
          ! or l_qv_star

  ! ----------------------------------------------------------------------
  ! Section BL.4 Convert and calculate theta and q fields from qT and Tl.
  ! ----------------------------------------------------------------------
  IF (.NOT. l_dry) THEN
    ! If the mixed phase precipitation scheme is used then T and Q are
    ! required to contain T liquid and Q(vapour+liquid) but at this stage
    ! will actually contain T liquid ice and Q(vapour+liquid+ice)

    CALL bl_lsp( bl_levels, qcf_latest, q_latest, t_latest )

  ELSE

    DO k = bl_levels+1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          q_latest(i,j,k) = 0.0
          qcl_latest(i,j,k) = 0.0
          qcf_latest(i,j,k) = 0.0
        END DO
      END DO
    END DO
  END IF

  ! Create Tl and qT outside boundary layer
!$OMP  PARALLEL DEFAULT(NONE)                                           &
!$OMP& SHARED( bl_levels, tdims, t_latest, qcl_latest, cp, q_latest,    &
!$OMP&         qcf_latest ) PRIVATE( i, j, k )
!$OMP DO SCHEDULE(STATIC)
  DO k = bl_levels+1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        t_latest(i,j,k) = t_latest(i,j,k) -                             &
                        (lc * qcl_latest(i,j,k))                        &
                         / cp
        q_latest(i,j,k) = q_latest(i,j,k) + qcl_latest(i,j,k)
        ! Prepare for cloud scheme. Are we using PC2 or not?
        ! zero any negative q_latests
        q_latest(i,j,k) = MAX(q_latest(i,j,k),0.0)
        qcf_latest(i,j,k) = MAX(qcf_latest(i,j,k),0.0)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, bl_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        q_latest(i,j,k) = MAX(q_latest(i,j,k),0.0)
        qcf_latest(i,j,k) = MAX(qcf_latest(i,j,k),0.0)
      END DO
    END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL

  ! Calculate diagnostic RHcrit or read as parameter in from namelist

  ! i_rhcpt_if1:
  IF (i_rhcpt == rhcpt_horiz_var) THEN
#if !defined(LFRIC)
    !       RHCRIT is 3D diagnosed variable
    ! Wet_mlev_do1:
    DO k = 1, tdims%k_end
      ! Rhc_rows_do1:
      DO j = 1, rhc_rows
        ! Rhc_rowlen_do1:
        DO i = 1, rhc_row_length
          ext_p_layer_centres(i,j,k) = p_layer_centres(i,j,k)
          ext_tl(i,j,k) = t_latest(i,j,k)
          ext_ql(i,j,k) = q_latest(i,j,k)
          ext_qcf(i,j,k) = qcf_latest(i,j,k)
        END DO ! Rhc_rowlen_do1
      END DO ! Rhc_rows_do1
    END DO ! Wet_mlev_do1

    ! Rhc_rows_do2:
    DO j = 1, rhc_rows
      ! Rhc_rowlen_do2:
      DO i = 1, rhc_row_length
        ext_p_layer_centres(i,j,0) = p_layer_centres(i,j,0)
        ext_ice_frac(i,j) = ice_fract(i,j)
        ! extended halo land fraction now passed in from AP2
      END DO ! Rhc_rowlen_do2
    END DO ! Rhc_rows_do2

    ! Synchronize haloes.

    i_field = 0
    i_field = i_field + 1
    fields_to_swap(i_field) % field       => ext_p_layer_centres(:,:,:)
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

    CALL swap_bounds_mv( fields_to_swap, i_field,                       &
                         rhc_row_length, 1, 1)

    i_field = 0

    i_field = i_field + 1
    fields_to_swap(i_field) % field_2d    => ext_ice_frac(:,:)
    fields_to_swap(i_field) % field_type  =  fld_type_p
    fields_to_swap(i_field) % levels      =  1
    fields_to_swap(i_field) % rows        =  rhc_rows
    fields_to_swap(i_field) % vector      =  .FALSE.

    CALL swap_bounds_2d_mv( fields_to_swap, i_field,                    &
                         rhc_row_length, 1, 1)

    CALL ls_calc_rhcrit( ext_p_layer_centres,                           &
    !              Array dimensions
             rhc_row_length, rhc_rows,                                  &
    !              Prognostic Fields
             ext_tl, ext_ql, ext_qcf,                                   &
             flandg(0:rhc_row_length,0:rhc_rows), ext_ice_frac,         &
    !              Output
             rhcpt)
#endif
  ELSE IF (i_rhcpt == rhcpt_off) THEN
    !         RHCRIT is 1D Parameter read in from namelist
    DO k = 1, tdims%k_end
      rhcpt(1,1,k) = rhcrit(k)
    END DO
  END IF  ! i_rhcpt_if1

  ! Which cloud scheme are we using?
  ! 3 options are available:
  ! Specific switching off cloud scheme (i_cld_off)
  ! PC2 cloud scheme (i_cld_pc2)
  ! Otherwise use Smith cloud scheme (i_cld_smith)

  IF ( i_cld_vn == i_cld_off ) THEN

    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end

          ! q_latest is total water (qv+qcl) that has been mixed by the
          ! boundary layer scheme. Subtract the original qcl
          ! to find the change in vapour due to mixing.
          ! ( N.B. qcf has been subtracted already in bl_lsp.F90 )

          q_latest(i,j,k) = q_latest(i,j,k) - (qcl_latest(i,j,k))

          ! Ensure that q_latest does no go negative
          IF ( q_latest(i,j,k) < mprog_min ) THEN

            ! Limit qcl to 10% of qt ( 10% is contained in qcl_max_factor )
            casqcl = qcl_max_factor * ( qcl_latest(i,j,k) + q_latest(i,j,k) )

            q_latest(i,j,k) = q_latest(i,j,k) + (qcl_latest(i,j,k) - casqcl)

            qcl_latest(i,j,k) = casqcl

          END IF ! q_latest < 0.0

          ! adjust the temperature to include the original qcl
          t_latest(i,j,k) = t_latest(i,j,k) + ( lcrcp * qcl_latest(i,j,k) )

          ! Set the cloud fractions to all-or-nothing cloud, depending
          ! on whether we have qcl above the minimum value used in the
          ! microphysics or not

          IF ( qcl_latest(i,j,k) > mprog_min ) THEN
            cfl_latest(i,j,k) = 1.0
          ELSE
            cfl_latest(i,j,k) = 0.0
          END IF

          IF ( qcf_latest(i,j,k) > mprog_min ) THEN
            cff_latest(i,j,k) = 1.0
          ELSE
            cff_latest(i,j,k) = 0.0
          END IF

          ! If cff_latest or cff_latest is 1 then bulk cloud fraction to 1
          ! otherwise set it to zero. This can be done avoiding another
          ! expensive If test by using the MAX function and then setting
          ! the other variables to that value.

          cf_bulk(i,j,k)    = MAX( cfl_latest(i,j,k), cff_latest(i,j,k) )

          cf_latest(i,j,k)  = cf_bulk(i,j,k)
          cf_area(i,j,k)    = cf_bulk(i,j,k)
          cf_liquid(i,j,k)  = cfl_latest(i,j,k)
          cf_frozen(i,j,k)  = cff_latest(i,j,k)

        END DO
      END DO
    END DO

  ELSE IF ( i_cld_vn == i_cld_pc2) THEN

    ! ----------------------------------------------------------------------
    ! PC2 cloud scheme
    ! ----------------------------------------------------------------------
    ALLOCATE ( zeros   (tdims%i_start:tdims%i_end,                      &
                        tdims%j_start:tdims%j_end,tdims%k_end) )

    ALLOCATE ( qt_force(tdims%i_start:tdims%i_end,                      &
                        tdims%j_start:tdims%j_end,tdims%k_end) )
    ALLOCATE ( tl_force(tdims%i_start:tdims%i_end,                      &
                        tdims%j_start:tdims%j_end,tdims%k_end) )
    ALLOCATE ( t_inc_pc2(tdims%i_start:tdims%i_end,                     &
                         tdims%j_start:tdims%j_end,tdims%k_end) )
    ALLOCATE ( q_inc_pc2(tdims%i_start:tdims%i_end,                     &
                         tdims%j_start:tdims%j_end,tdims%k_end) )
    ALLOCATE ( qcl_inc_pc2(tdims%i_start:tdims%i_end,                   &
                           tdims%j_start:tdims%j_end,tdims%k_end) )
    ALLOCATE ( cfl_inc_pc2(tdims%i_start:tdims%i_end,                   &
                           tdims%j_start:tdims%j_end,tdims%k_end) )
    ALLOCATE ( bcf_inc_pc2(tdims%i_start:tdims%i_end,                   &
                           tdims%j_start:tdims%j_end,tdims%k_end) )

    ! ----------------------------------------------------------------------
    ! Inhomogenous forcing of ice
    ! ----------------------------------------------------------------------

    ! Calculate in-plume ice content (LS) by assuming that they are equal to
    ! the current values, except when the current value is not defined.
    ! Also calculate the forcing of ice content Q4F

    IF (l_fixbug_pc2_mixph) THEN
!$OMP  PARALLEL DEFAULT(NONE)                                           &
!$OMP& SHARED( tdims, l_fixbug_pc2_mixph, qcf_latest, qcf_earliest,     &
!$OMP&         cff_earliest, cff_latest, cf_latest, zeros, qt_force,    &
!$OMP&         q_latest, q_earliest, qcl_earliest, tl_force, t_latest,  &
!$OMP&         t_earliest, cp )                                         &
!$OMP& PRIVATE( i, j, k, q4, denom, deltacff, incloudice, temp3 )
!$OMP  DO SCHEDULE(DYNAMIC)
      DO k = 1, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end

            ! Calculate Q4. Only perform the calculation if the Q4 is non-zero.
            ! Since ice is mixed by tracer mixing regardless of whether cumulus
            ! is present we still need to provide an ice cloud fraction
            ! increment below cumulus base.

            ! Calculate change in ice cloud fraction differently depending on
            ! whether QCF has increased or decreased.

            q4 = qcf_latest(i,j,k) - qcf_earliest(i,j,k)

            IF (q4 > cloud_rounding_tol) THEN

              ! Source.

              ! Calculate the change in total cloud fraction.
              ! Use a weighted (by ice cloud fraction) average of in-cloud ice
              ! content and a fixed value to specify the plume ice content.
              ! The denominator in the deltaCf calculation is then
              ! ((qcf_earliest/cff_earliest)*cffearliest +
              ! ls_bl0*(1-cff_earliest) - qcf_earliest
              ! and the qcf_earliest terms cancel.

              denom = ls_bl0 * ( 1.0 - cff_earliest(i,j,k) )

              IF ( ABS(denom) > 1.0e-10 ) THEN

                denom = q4 / denom

                deltacff = (1.0-cff_latest(i,j,k)) * denom
                cff_latest(i,j,k) = cff_latest(i,j,k) + deltacff

                ! calculate total cf based on minimum overlap
                cf_latest(i,j,k)  = cf_latest(i,j,k)  + deltacff

              ELSE
                cf_latest(i,j,k)  = 1.0
                cff_latest(i,j,k) = 1.0
              END IF

            ELSE IF (q4 < -cloud_rounding_tol) THEN

              ! Sink.

              ! Given a reduction in qcf, remove some CFF in order to
              ! maintain the same in-cloud IWC.

              incloudice = qcf_latest(i,j,k) /                          &
                MAX(cff_latest(i,j,k),1.0e-6)

              cf_latest(i,j,k) = cf_latest(i,j,k) +                     &
                q4/MAX(incloudice, ls_bl0)

              cff_latest(i,j,k) = cff_latest(i,j,k) +                   &
                q4/MAX(incloudice, ls_bl0)

              IF (cff_latest(i,j,k) > 0.0 .AND.                         &
                !Prevent very high in-cloud values by increasing CFF.
                 (qcf_latest(i,j,k) / cff_latest(i,j,k) > 2.0e-3)) THEN
                temp3=cff_latest(i,j,k)
                cff_latest(i,j,k)=qcf_latest(i,j,k)/2.0e3
                temp3=cff_latest(i,j,k)-temp3
                  ! Update total cloud fraction.
                cf_latest(i,j,k)=cf_latest(i,j,k)+temp3
              END IF

            END IF
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT

      ! ----------------------------------------------------------------------
      ! Homogenous forcing of the liquid cloud
      ! ----------------------------------------------------------------------

      ! Calculate forcing in qT and TL. Currently q_latest contains the vapour
      ! plus liquid content and q_earliest just the initial vapour content

!$OMP DO SCHEDULE(STATIC)
      DO k = 1, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            ! Set zero array
            zeros(i,j,k)        = 0.0
            qt_force(i,j,k) = ( q_latest(i,j,k)                         &
              - (q_earliest(i,j,k) + qcl_earliest(i,j,k)) )
            tl_force(i,j,k) = ( t_latest(i,j,k)                         &
              - (t_earliest(i,j,k)- lc * qcl_earliest(i,j,k) / cp) )
          END DO
        END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL

    ELSE

      DO k = 1, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end

            ! Original code as run in Wilson et al (2008a,b).
            ! Calculation of change in cloud fraction is consitent with
            ! an increase in QCF but not consistent with a decrease.

            q4 = qcf_latest(i,j,k) - qcf_earliest(i,j,k)
            IF (q4  /=  0.0) THEN

              ! Calculate the change in total cloud fraction.
              ! Use a weighted (by ice cloud fraction) average of in-cloud ice
              ! content and a fixed value to specify the plume ice content.
              ! The denominator in the deltaCf calculation is then
              ! (qcf_earliest/cff_earliest + ls_bl0*(1-cff_earliest)
              !   - qcf_earliest
              ! and the qcf_earliest terms cancel.

              denom = ls_bl0 * (1.0 - cff_earliest(i,j,k))

              IF ( ABS(denom)  >   1.0e-10 ) THEN

                denom = q4 / denom
                cf_latest(i,j,k)  = cf_latest(i,j,k)  +                 &
                             (1.0 - cf_latest(i,j,k))  * denom
                cff_latest(i,j,k) = cff_latest(i,j,k) +                 &
                             (1.0 - cff_latest(i,j,k)) * denom

                ! Otherwise cloud fraction will go to one. In theory,
                ! cloud fraction can never be reduced by this process.

              ELSE
                cf_latest(i,j,k)  = 1.0
                cff_latest(i,j,k) = 1.0
              END IF

            END IF !(Q4  /=  0.0)

          END DO
        END DO
      END DO

      ! ----------------------------------------------------------------------
      ! Homogenous forcing of the liquid cloud
      ! ----------------------------------------------------------------------

      ! Calculate forcing in qT and TL. Currently q_latest contains the vapour
      ! plus liquid content and q_earliest just the initial vapour content

      DO k = 1, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            ! Set zero array
            zeros(i,j,k)        = 0.0
            qt_force(i,j,k) = ( q_latest(i,j,k)                         &
              - (q_earliest(i,j,k) + qcl_earliest(i,j,k)) )
            tl_force(i,j,k) = ( t_latest(i,j,k)                         &
              - (t_earliest(i,j,k)- lc * qcl_earliest(i,j,k) / cp) )
          END DO
        END DO
      END DO
    END IF ! l_fixbug_pc2_mixph

    IF ( l_leonard_term ) THEN
      ! Add turbulent increment due to the Leonard terms to the
      ! turbulent forcing of Tl and qw.  These increments were
      ! calculated earlier in atmos_physics2 and get counted as
      ! part of the "non-turbulent" increment in imp_solver, but
      ! should contribute to the turbulent forcing of cloud here.
      DO k = 1, bl_levels
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            qt_force(i,j,k) = qt_force(i,j,k)                           &
                            + qw_inc_leonard(i,j,k)
            tl_force(i,j,k) = tl_force(i,j,k)                           &
                            + thetal_inc_leonard(i,j,k)                 &
                              * exner_theta_levels(i,j,k)
          END DO
        END DO
      END DO
    END IF

    ! Call homogenous forcing routine

    CALL pc2_delta_hom_turb(                                            &
    ! INput variables
            p_layer_centres(1,1,1),                                     &
    ! INput variables
            t_earliest(1,1,1), q_earliest(1,1,1), qcl_earliest(1,1,1),  &
            cf_latest(1,1,1), cfl_latest(1,1,1), cff_latest(1,1,1),     &
            tl_force(1,1,1),qt_force(1,1,1),zeros(1,1,1),zeros(1,1,1),  &
    ! OUTput variables
            t_inc_pc2, q_inc_pc2, qcl_inc_pc2, bcf_inc_pc2, cfl_inc_pc2,&
    ! INput variables (other quantities)
            0.0, 0.0, l_mr_physics)

    IF ( l_leonard_term ) THEN
      ! Leonard term increments were already added in atmos_physics2
      ! and are included in t_earliest / q_earliest here.
      ! We have added them to tl_force / qw_force above to include
      ! the cloud response to Leonard terms.
      ! But later the _latest variables are updated by adding _force
      ! to _earliest, so we'd be double-counting the Leonard terms.
      ! Therefore, subtract Leonard terms off the _force variables here:
      DO k = 1, bl_levels
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            qt_force(i,j,k) = qt_force(i,j,k)                           &
                            - qw_inc_leonard(i,j,k)
            tl_force(i,j,k) = tl_force(i,j,k)                           &
                            - thetal_inc_leonard(i,j,k)                 &
                              * exner_theta_levels(i,j,k)
          END DO
        END DO
      END DO
    END IF

    ! Diagnostic shallow cloud within PC2

    IF ( l_pc2_diag_sh .AND. .NOT. l_ccrad .AND.                        &
         rad_cloud_decay_opt == rad_decay_off ) THEN

          ! The following code is only valid with the follow settings.

          ! i_cld_vn          = i_cld_pc2
          ! l_pc2_diag_sh     = .true.
          ! l_ccrad           = .false.
          ! Rad_cloud_decay_opt = rad_decay_off

      ALLOCATE ( ccw_cca(tdims%i_start:tdims%i_end,                     &
                         tdims%j_start:tdims%j_end,tdims%k_end) )
      ALLOCATE ( cca_3d (tdims%i_start:tdims%i_end,                     &
                         tdims%j_start:tdims%j_end,tdims%k_end) )

      DO k = 1, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end

            ! If below cumulus cloud base or we are using a diagnostic
            ! shallow convective cloud then we simply zero the liquid water
            ! content instead of using the homogenous BL response. The forcings
            ! themselves are still applied but to q and T.
            ! Note that we find, for PC2, that zeroing up to ntml,
            ! rather than ntml+1, is more physically justified and gives
            ! slightly better cloud results.

            IF ( (l_pc2_diag_sh_pts(i,j) .AND. k <= cct0(i,j)) .OR.     &
                 (cumulus(i,j)           .AND. k <= ntml(i,j)) ) THEN
              t_inc_pc2(i,j,k)   = (-lcrcp) * qcl_earliest(i,j,k)
              q_inc_pc2(i,j,k)   = qcl_earliest(i,j,k)
              qcl_inc_pc2(i,j,k) = (-qcl_earliest(i,j,k))
              cfl_inc_pc2(i,j,k) = (-cfl_earliest(i,j,k))
              bcf_inc_pc2(i,j,k) = cff_latest(i,j,k)-cf_earliest(i,j,k)
            END IF

                ! Set convective cloud properties if we are
                ! using a diagnostic shallow convection for PC2
                ! Calculate convective properties

                ! With i_cld_vn==i_cld_pc2, the cca0 which goes
                ! to radiation should only contain only
                ! shallow convection (i.e. no deep or mid-level).

                ! l_pc2_diag_sh and l_ccrad should not be
                ! true at the same time.
            IF (l_3d_cca) THEN

              ccw_cca(i,j,k) = ccw0(i,j,k)*cca0(i,j,k)
              cca_3d(i,j,k)  = cca0(i,j,k)

            ELSE

                  ! Currently with this switch cca is set to be
                  ! a single level field. In which case use cca0_2d
                  ! which should also only contain cca0_2d from
                  ! shallow cloud

              IF (k <= cct0(i,j) .AND. k >= ccb0(i,j)) THEN
                ccw_cca(i,j,k) = ccw0(i,j,k) * cca0_2d(i,j)
                cca_3d(i,j,k)  = cca0_2d(i,j)
              ELSE
                ccw_cca(i,j,k) = 0.0
                cca_3d(i,j,k)  = 0.0
              END IF
            END IF

            IF ( l_pc2_diag_sh_pts(i,j) .AND. k <= cct0(i,j)            &
                          .AND. t_earliest(i,j,k) >= tice) THEN

                  ! Hand over convective attributes to the large scale
              t_inc_pc2(i,j,k)   =  t_inc_pc2(i,j,k)                    &
                                  + lcrcp * ccw_cca(i,j,k)
              q_inc_pc2(i,j,k)   =  q_inc_pc2(i,j,k) -  ccw_cca(i,j,k)
              qcl_inc_pc2(i,j,k) =  qcl_inc_pc2(i,j,k)                  &
                                  + ccw_cca(i,j,k)
              cfl_inc_pc2(i,j,k) =  cfl_inc_pc2(i,j,k)                  &
                                   + cca_3d(i,j,k)
              bcf_inc_pc2(i,j,k) =  bcf_inc_pc2(i,j,k)                  &
                              + cca_3d(i,j,k) *(1.0-cff_latest(i,j,k))

            END IF

                ! Reset the section 0 convective cloud back to zero
                ! now that increments have been included in PC2
                ! increments.
            cca0(i,j,k)  = 0.0
            cca0_2d(i,j) = 0.0
            ccw0(i,j,k)  = 0.0

          END DO  ! i loop
        END DO  ! j
      END DO  ! k

      DEALLOCATE (ccw_cca)
      DEALLOCATE (cca_3d)

    ELSE     ! original code

!$OMP  PARALLEL DEFAULT(NONE)                                           &
!$OMP& SHARED( tdims, cumulus, ntml, forced_cu, bl_type_3, bl_type_4,   &
!$OMP&         z_theta, t_inc_pc2, qcl_earliest, q_inc_pc2,l_param_conv,&
!$OMP&         qcl_inc_pc2, cfl_inc_pc2, cfl_earliest, bcf_inc_pc2,     &
!$OMP&         cff_latest, cf_earliest, zlcl, zlcl_mixed, lcrcp )       &
!$OMP& PRIVATE( i, j, k )
!$OMP DO SCHEDULE(DYNAMIC)
      DO k = 1, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end

            ! If below the LCL in cumulus or unstable (but not shear-dominated)
            ! boundary layers then we simply zero the liquid water content
            ! instead of using the homogenous BL response. The forcings
            ! themselves are still applied but to q and T.
            IF ( ( cumulus(i,j) .AND. k  <=  ntml(i,j) ) .OR.           &
                 ( forced_cu >= on .AND. (bl_type_3(i,j) > 0.5          &
                    .OR. bl_type_4(i,j) > 0.5 )                         &
                    .AND. z_theta(i,j,k)  <  zlcl(i,j) )  ) THEN
              t_inc_pc2(i,j,k)   =  (-lcrcp) * qcl_earliest(i,j,k)
              q_inc_pc2(i,j,k)   =  qcl_earliest(i,j,k)
              qcl_inc_pc2(i,j,k) =  (-qcl_earliest(i,j,k))
              cfl_inc_pc2(i,j,k) =  (-cfl_earliest(i,j,k))
              bcf_inc_pc2(i,j,k) =  cff_latest(i,j,k)                   &
                                      -cf_earliest(i,j,k)
            END IF
          END DO  ! i loop
        END DO  ! j
      END DO  ! k
!$OMP END DO NOWAIT

! To be consistent with the code above, set zlcl_mixed to
! prevent PC2 initiating cloud below this level

! Only needed for cumulus layer if l_param_conv is false - we can't
! initiate in cumulus when convection is active
      IF (.NOT. l_param_conv) THEN
!$OMP DO SCHEDULE(DYNAMIC)
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            IF ( cumulus(i,j) ) zlcl_mixed(i,j) = zlcl(i,j)
          END DO  ! i loop
        END DO  ! j
!$OMP END DO
      END IF

! With forced_cu, also do this for types 3 & 4
      IF (forced_cu >= on) THEN
!$OMP DO SCHEDULE(DYNAMIC)
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            IF (bl_type_3(i,j) > 0.5 .OR. bl_type_4(i,j) > 0.5 )        &
                        zlcl_mixed(i,j) = zlcl(i,j)
          END DO  ! i loop
        END DO  ! j
!$OMP END DO
      END IF
!$OMP END PARALLEL

    END IF    ! test on l_pc2_diag_sh

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                       &
!$OMP& SHARED( tdims, t_latest, t_earliest, tl_force, t_inc_pc2,        &
!$OMP&         q_latest, q_earliest, qt_force, q_inc_pc2, qcl_latest,   &
!$OMP&         qcl_earliest, qcl_inc_pc2, cfl_latest, cfl_earliest,     &
!$OMP&         cfl_inc_pc2, cf_latest, cf_earliest, bcf_inc_pc2 )       &
!$OMP& PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end

          ! Update working version of temperature, moisture and cloud
          ! fields with increments from the PC2 homogenous response.

          t_latest(i,j,k)   = t_earliest(i,j,k) + tl_force(i,j,k)       &
                             + t_inc_pc2(i,j,k)
          q_latest(i,j,k)   = q_earliest(i,j,k) + qt_force(i,j,k)       &
                             + q_inc_pc2(i,j,k)
          qcl_latest(i,j,k) = qcl_earliest(i,j,k)                       &
                             + qcl_inc_pc2(i,j,k)
          cfl_latest(i,j,k) = cfl_earliest(i,j,k)                       &
                             + cfl_inc_pc2(i,j,k)
          cf_latest(i,j,k)  =  cf_earliest(i,j,k)                       &
                             + bcf_inc_pc2(i,j,k)
          !               qcf_latest and cff_latest are not updated

        END DO  ! i loop
      END DO  ! j
    END DO  ! k
!$OMP END PARALLEL DO

    DEALLOCATE (zeros)
    !------------------------------------------------------------
    ! Parametrize "forced cumulus clouds" at top of well-mixed BL
    !------------------------------------------------------------
    IF (forced_cu >= on) THEN

!$OMP  PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(NONE)                      &
!$OMP& SHARED( tdims, zhnl, dzh, zlcl, bl_type_3, z_theta,              &
!$OMP&         cfl_latest, qcl_inv_top, cf_latest, qcl_latest,          &
!$OMP&         q_latest, forced_cu_fac, t_latest, lcrcp )               &
!$OMP& PRIVATE ( i, j, k, zc_depth, cf_base, cf_forced, qcl_forced,     &
!$OMP&           dqcl, qcl_tol, dcfl )
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          zc_depth = zhnl(i,j)+dzh(i,j)-zlcl(i,j)
          IF ( zc_depth > 1.0 .AND. bl_type_3(i,j)>0.5 ) THEN

            DO k = 1, tdims%k_end
              IF ( z_theta(i,j,k) >= zlcl(i,j)                          &
                  .AND. z_theta(i,j,k) <= zhnl(i,j)+dzh(i,j)            &
                  .AND. cfl_latest(i,j,k) < 0.5                         &
                 ) THEN
                  ! Make cloud fraction at cloud base a function of
                  ! the cloud depth:
                cf_base = MAX( cf_top,                                  &
                               0.3 * MIN( 1.0, zc_depth/300.0 ) )
                cf_forced = cf_base - (cf_base-cf_top) *                &
                                   (z_theta(i,j,k)-zlcl(i,j)) / zc_depth
                ! Use diagnostic parcel qcl for in-cloud water content
                qcl_forced = MAX( qcl_forced_min,                       &
                  qcl_inv_top(i,j)*(z_theta(i,j,k)-zlcl(i,j))/zc_depth )
                qcl_forced = forced_cu_fac*qcl_forced ! tuning knob
                qcl_forced = qcl_forced*cf_forced ! GBM water content
                !-----------------------------------------------------
                ! calculate resulting changes in moisture variables
                ! and check for low values
                !-----------------------------------------------------
                dqcl = 0.0
                IF ( qcl_forced > qcl_latest(i,j,k) ) THEN
                  dqcl = qcl_forced - qcl_latest(i,j,k)
                  ! make sure we won't condense more than a fraction of
                  ! the existing water vapour
                  qcl_tol = qcl_max_factor*q_latest(i,j,k)
                  IF ( qcl_forced > qcl_tol ) THEN
                    qcl_forced = qcl_tol
                    dqcl = MAX( 0.0, qcl_forced - qcl_latest(i,j,k) )
                    qcl_forced = qcl_latest(i,j,k) + dqcl
                  END IF
                  IF (q_latest(i,j,k)-dqcl < mprog_min) THEN
                    ! almost no water in this grid-point so leave it!
                    qcl_forced = qcl_latest(i,j,k)
                    dqcl = 0.0
                  END IF
                  qcl_latest(i,j,k) = qcl_forced
                  t_latest(i,j,k)  = t_latest(i,j,k) + lcrcp*dqcl
                  q_latest(i,j,k)  = q_latest(i,j,k) - dqcl
                END IF ! test on qcl_forced

                IF ( cf_forced > cfl_latest(i,j,k) .AND.                &
                     dqcl > qcl_max_factor*mprog_min ) THEN
                  dcfl = cf_forced - cfl_latest(i,j,k)
                  cfl_latest(i,j,k) = cf_forced
                  cf_latest(i,j,k)  = cf_latest(i,j,k) + dcfl
                END IF

              END IF  ! test on z
            END DO  ! loop over k

          END IF  ! test on zc_depth and bl_type3
        END DO
      END DO
!$OMP END PARALLEL DO
    END IF  ! test on forced_cu >= on

    IF ( kprof_cu >= on .AND. forced_cu == cbl_and_cu ) THEN

      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          zc_depth = zhnl(i,j)-zlcl(i,j)
          IF ( zc_depth > 1.0 .AND. bl_type_6(i,j)>0.5 ) THEN
            DO k = 1, tdims%k_end
              !------------------------------------------------------------
              ! check cumulus cloud chunky enough above LCL
              !  - if kprof_cu used then zhnl penetrates above the LCL
              !    and so marks the top of "forced" cumulus
              !------------------------------------------------------------
              IF ( z_theta(i,j,k) >= zlcl(i,j)                          &
                  .AND. z_theta(i,j,k) <= zhnl(i,j)                     &
                  .AND. cfl_latest(i,j,k) < 0.5                         &
                 ) THEN
                ! Make cloud fraction at cloud base a function of
                ! the cloud depth:
                cf_base = MAX( cf_top,                                  &
                               0.3 * MIN( 1.0, zc_depth/300.0 ) )
                cf_forced = cf_base - (cf_base-cf_top) *                &
                                   (z_theta(i,j,k)-zlcl(i,j)) / zc_depth
                ! Use diagnostic parcel qcl for in-cloud water content
                ! For cumulus layers this is taken at zlcl+300m
                qcl_forced = MAX( qcl_forced_min,                       &
                     qcl_inv_top(i,j)*(z_theta(i,j,k)-zlcl(i,j))/300.0 )
                qcl_forced = forced_cu_fac*qcl_forced !tuning knob
                qcl_forced = qcl_forced*cf_forced ! GBM water content
                !-----------------------------------------------------
                ! calculate resulting changes in moisture variables
                ! and check for low values
                !-----------------------------------------------------
                dqcl = 0.0
                IF ( qcl_forced > qcl_latest(i,j,k) ) THEN
                  dqcl = qcl_forced - qcl_latest(i,j,k)
                  ! make sure we won't condense more than a fraction of
                  ! the existing water vapour
                  qcl_tol = qcl_max_factor*q_latest(i,j,k)
                  IF ( qcl_forced > qcl_tol ) THEN
                    qcl_forced = qcl_tol
                    dqcl = MAX( 0.0, qcl_forced - qcl_latest(i,j,k) )
                    qcl_forced = qcl_latest(i,j,k) + dqcl
                  END IF
                  IF (q_latest(i,j,k)-dqcl < mprog_min) THEN
                    ! almost no water in this grid-point so leave it!
                    qcl_forced = qcl_latest(i,j,k)
                    dqcl = 0.0
                  END IF
                  qcl_latest(i,j,k) = qcl_forced
                  t_latest(i,j,k)  = t_latest(i,j,k) + lcrcp*dqcl
                  q_latest(i,j,k)  = q_latest(i,j,k) - dqcl
                END IF ! test on qcl_forced

                IF ( cf_forced > cfl_latest(i,j,k) .AND.                &
                     dqcl > qcl_max_factor*mprog_min  ) THEN
                  dcfl = cf_forced - cfl_latest(i,j,k)
                  cfl_latest(i,j,k) = cf_forced
                  cf_latest(i,j,k)  = cf_latest(i,j,k) + dcfl
                END IF

              END IF  ! test on z

            END DO  ! loop over k
          END IF  ! test on zc_depth and bltype6
        END DO
      END DO

    END IF  ! test on forced_cu eq cbl_and_cu

    DEALLOCATE ( bcf_inc_pc2 )
    DEALLOCATE ( cfl_inc_pc2 )
    DEALLOCATE ( qcl_inc_pc2 )
    DEALLOCATE ( q_inc_pc2 )
    DEALLOCATE ( t_inc_pc2 )
    DEALLOCATE (tl_force)
    DEALLOCATE (qt_force)
    ! ----------------------------------------------------------------------
    ! Copy updated cloud fractions to the in/out variables
    ! ----------------------------------------------------------------------

    IF (i_cld_area == acf_off) THEN

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                       &
!$OMP& SHARED( tdims, cf_area, cf_latest ) PRIVATE( i, j, k )
      DO k = 1, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
                ! For the moment set area cloud fraction
                ! to the bulk cloud fraction
                ! Ensure it has a value between 0.0 and 1.0
            cf_area(i,j,k) = MAX(MIN(cf_latest(i,j,k),1.0),0.0)
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO

    ELSE IF (i_cld_area == acf_brooks) THEN

      CALL ls_acf_brooks (                                              &
           xx_cos_theta_latitude,                                       &
           cf_latest, cfl_latest, cff_latest, cumulus, cf_area )

      ! call sub-level interpolation parameterisation of cloud area
    ELSE IF (i_cld_area == acf_cusack) THEN

      ! Determine number of sublevels for vertical gradient area cloud
      ! Want an odd number of sublevels per level: 3 is hardwired in do loops
      levels_per_level = 3
      large_levels = ((tdims%k_end - 2)*levels_per_level) + 2

      CALL pc2_hom_arcld(p_layer_centres,p_layer_boundaries,            &
       large_levels,levels_per_level,                                   &
       cf_area, t_latest, cf_latest, cfl_latest, cff_latest, q_latest,      &
       qcl_latest, qcf_latest, l_mr_physics)

    END IF ! i_cld_area

    ! ----------------------------------------------------------------------
    ! Provide an estimate of convective cloud fraction for visibility
    ! ----------------------------------------------------------------------

    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        IF (cca_at_base(i,j) == 0.0 .AND.                               &
           (conv_rain(i,j) + conv_snow(i,j)) > 0.0) THEN
               ! Convective precipitation exists but no
               ! estimate for its cloud fraction. Set it
               ! to a constant value of 0.2 as an estimate.
          cca_at_base(i,j) = 0.2
        END IF
      END DO
    END DO

    ! ----------------------------------------------------------------------
    !  PC2: End of cloud section
    ! ----------------------------------------------------------------------

  ELSE IF (i_cld_vn == i_cld_smith) THEN

    ! ----------------------------------------------------------------------
    ! Section BL.4b Call cloud scheme to convert Tl and qT to T, q and qcl
    ! in boundary layer, calculate bulk_cloud fields from qT and qcf
    ! and calculate area_cloud fields.
    ! ----------------------------------------------------------------------

    ! Determine number of sublevels for vertical gradient area cloud
    ! Want an odd number of sublevels per level: 3 is hardwired in do loops
    levels_per_level = 3
    large_levels = ((tdims%k_end - 2)*levels_per_level) + 2

    CALL ls_arcld( p_layer_centres, rhcpt, p_layer_boundaries,          &
                 rhc_row_length, rhc_rows, bl_levels,                   &
                 levels_per_level, large_levels,                        &
                 xx_cos_theta_latitude,                                 &
                 ntml, cumulus, l_mr_physics, qcf_latest,               &
                 t_latest, q_latest, qcl_latest,                        &
                 cf_area, cf_bulk, cf_liquid, cf_frozen,                &
                 error_code)

  END IF   ! i_cld_vn

  ! Store T increment of BL
  IF (l_retain_bl_tendencies) THEN
    ! get temperature
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          dt_bl(i,j,k) = t_latest(i,j,k) - t_earliest(i,j,k)

        END DO ! i
      END DO ! j
    END DO ! k
  END IF

  ! Store q_cl increment of BL
  IF (l_retain_q_cl_bl_tendencies) THEN
    ! get temperature
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          dq_cl_bl(i,j,k) = qcl_latest(i,j,k)  - qcl_earliest(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
  END IF


  IF ( BL_diag%l_t_incr .OR. BL_diag%l_tl_incr ) THEN
    ALLOCATE ( BL_diag%t_incr(tdims%i_start:tdims%i_end,                &
                              tdims%j_start:tdims%j_end,                &
                                          1:tdims%k_end) )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          BL_diag%t_incr(i,j,k) = t_latest(i,j,k) - t_earliest(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
  ELSE
    ALLOCATE ( BL_diag%t_incr(1,1,1) )
  END IF                   ! on STASHflags

  IF ( BL_diag%l_q_incr .OR. BL_diag%l_qtl_incr ) THEN
    ALLOCATE ( BL_diag%q_incr(tdims%i_start:tdims%i_end,                &
                              tdims%j_start:tdims%j_end,                &
                                          1:tdims%k_end) )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          BL_diag%q_incr(i,j,k) = q_latest(i,j,k) - q_earliest(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
  ELSE
    ALLOCATE ( BL_diag%q_incr(1,1,1) )
  END IF                  ! on STASHflags

  IF ( BL_diag%l_qcl_incr .OR. BL_diag%l_tl_incr .OR.                   &
       BL_diag%l_qtl_incr ) THEN
    ALLOCATE ( BL_diag%qcl_incr(tdims%i_start:tdims%i_end,              &
                                tdims%j_start:tdims%j_end,              &
                                            1:tdims%k_end) )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          BL_diag%qcl_incr(i,j,k) = qcl_latest(i,j,k)                   &
                                  - qcl_earliest(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
  ELSE
    ALLOCATE ( BL_diag%qcl_incr(1,1,1) )
  END IF                  ! on STASHflag

  IF ( BL_diag%l_qcf_incr ) THEN
    ALLOCATE (BL_diag%qcf_incr(tdims%i_start:tdims%i_end,               &
                               tdims%j_start:tdims%j_end,               &
                                           1:tdims%k_end))
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          BL_diag%qcf_incr(i,j,k) = qcf_latest(i,j,k)                   &
                                  - qcf_earliest(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
  ELSE
    ALLOCATE ( BL_diag%qcf_incr(1,1,1) )
  END IF                  ! on STASHflag

  IF ( BL_diag%l_cf_incr ) THEN
    ALLOCATE (BL_diag%cf_incr(tdims%i_start:tdims%i_end,                &
                              tdims%j_start:tdims%j_end,                &
                                          1:tdims%k_end))
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          BL_diag%cf_incr(i,j,k) = cf_latest(i,j,k)                     &
                                 - cf_earliest(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
  ELSE
    ALLOCATE ( BL_diag%cf_incr(1,1,1) )
  END IF

  IF ( BL_diag%l_cfl_incr ) THEN
    ALLOCATE (BL_diag%cfl_incr(tdims%i_start:tdims%i_end,               &
                               tdims%j_start:tdims%j_end,               &
                                           1:tdims%k_end))
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          BL_diag%cfl_incr(i,j,k) = cfl_latest(i,j,k)                   &
                                  - cfl_earliest(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
  ELSE
    ALLOCATE ( BL_diag%cfl_incr(1,1,1) )
  END IF

  IF ( BL_diag%l_cff_incr ) THEN
    ALLOCATE (BL_diag%cff_incr(tdims%i_start:tdims%i_end,               &
                               tdims%j_start:tdims%j_end,               &
                                           1:tdims%k_end))
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          BL_diag%cff_incr(i,j,k) = cff_latest(i,j,k)                   &
                                  - cff_earliest(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
  ELSE
    ALLOCATE ( BL_diag%cff_incr(1,1,1) )
  END IF

END IF ! on error code zero

! ----------------------------------------------------------------------
! Section BL.4c Combined cloud field calculation for use by visibility
!               (section 3) and cloud scheme (section 9) diagnostics
!               09208 - 09217 and 09223.
! ----------------------------------------------------------------------

! L_combi_cld_if1:
IF (error_code  <=  0  .AND.  l_combi_cld .AND. l_apply_diag) THEN
      ! Set the combined cloud area fractions in each gridbox.
      ! Convention in Sect 70 (Radiation) is to invert levels, 1 at top.

  ALLOCATE                                                              &
    ( combined_cloud (tdims%i_start:tdims%i_end,                        &
                      tdims%j_start:tdims%j_end,                        &
                      tdims%k_end),                                     &
      cca4comb_cld   (tdims%i_start:tdims%i_end,                        &
                      tdims%j_start:tdims%j_end,                        &
                      n_cca_levels),                                    &
      ccb4comb_cld   (tdims%i_start:tdims%i_end,                        &
                      tdims%j_start:tdims%j_end),                       &
      cct4comb_cld   (tdims%i_start:tdims%i_end,                        &
                      tdims%j_start:tdims%j_end) )

  nclds = MIN(cloud_levels, tdims%k_end)

      ! ***** Code adapted from R2_SET_CLOUD_FIELD. *****

      ! Zero cloud amounts in the upper layers (if necessary).
  ! Nclds_if1:
  IF (tdims%k_end  >   nclds) THEN
    ! Rad_k_do1:
    DO k = 1, tdims%k_end-nclds
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          combined_cloud(i, j, k) = 0.0e+00
        END DO
      END DO
    END DO  ! Rad_k_do1
  END IF  ! Nclds_if1

      ! Use Sec 0 convective cloud
  DO k = 1, n_cca_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        cca4comb_cld(i,j,k) = cca0(i,j,k)
      END DO
    END DO
  END DO
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      ccb4comb_cld(i,j) = ccb0(i,j)
      cct4comb_cld(i,j) = cct0(i,j)
    END DO
  END DO

      ! Calculate combined cloud field
  DO k = tdims%k_end+1-nclds, tdims%k_end

    kinvert = tdims%k_end+1 - k

    IF (l_3d_cca) THEN

      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          combined_cloud(i,j,k) = cca4comb_cld(i,j,kinvert)             &
                                + (1.0 - cca4comb_cld(i,j,kinvert))     &
                                * cf_area(i,j,kinvert)
        END DO
      END DO

    ELSE

      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          IF ( (cct4comb_cld(i,j) >= kinvert+1) .AND.                   &
               (ccb4comb_cld(i,j) <= kinvert) ) THEN
            combined_cloud(i,j,k) = cca4comb_cld(i,j,1)                 &
                                  + (1.0-cca4comb_cld(i,j,1))           &
                                  * cf_area(i,j,kinvert)
          ELSE
            combined_cloud(i,j,k) = cf_area(i,j,kinvert)
          END IF
        END DO
      END DO

    END IF  ! l_3d_cca

  END DO  ! k

ELSE
  ALLOCATE(combined_cloud(1,1,1))
END IF  ! L_combi_cld_if1

! NB: Combined cloud area fractions in each gridbox set up above.
!     Convention in Sect 70 (Radiation) is to invert levels, 1 at top.

! L_plsp_if1:
IF ( error_code <= 0  .AND.  l_plsp .AND. l_apply_diag ) THEN

  CALL r2_calc_total_cloud_cover(                                       &
         pdims%i_end*pdims%j_end, tdims%k_end, nclds,                   &
         ip_cloud_mix_max, combined_cloud(1,1,1), work2d_1,             &
         pdims%i_end*pdims%j_end, tdims%k_end                           &
       )

  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      IF (cca_at_base(i,j) < 1.0) THEN
        plsp(i,j) = MAX( 0.0,                                           &
                         (work2d_1(i,j)-cca_at_base(i,j))               &
                          / (1.0 - cca_at_base(i,j)) )
      ELSE
        plsp(i,j) = 0.0
      END IF
    END DO
  END DO

END IF  ! L_plsp_if1:

IF ( l_apply_diag ) THEN
  ! ----------------------------------------------------------------------
  ! Section BL.4d Cloud scheme (section 9) diagnostics.
  ! ----------------------------------------------------------------------

  ! Check that cloud diagnostics are requested this timestep
  ! DiagSect09_if1:
  IF (model_type /= mt_single_column) THEN
    IF (error_code  ==  0 .AND. sf(0,9)) THEN
#if !defined(LFRIC)
      ! DEPENDS ON: diagnostics_lscld
      CALL diagnostics_lscld(                                           &
                         pdims%i_end,pdims%j_end,                       &
                         rhc_row_length, rhc_rows,                      &
                         bl_levels, cloud_levels,                       &
                    p_layer_centres(1,1,1),p_layer_boundaries(1,1,1),   &
                         BL_diag%t_incr, BL_diag%q_incr,                &
                         BL_diag%qcl_incr,                              &
                         t_latest, q_latest, qcl_latest, qcf_latest,    &
                         cf_area, cf_bulk,                              &
                         cfl_latest, cff_latest,                        &
                         rho_wet_tq,p_star, rhcpt,                      &
                         combined_cloud, aerosol,                       &
                         xx_cos_theta_latitude,                         &
                         stashwork9,                                    &
                         l_combi_cld                                    &
                          )
#endif
    END IF ! DiagSect09_if1
  END IF ! model_type

  ! ----------------------------------------------------------------------
  ! Section BL.5 Energy correction
  ! ----------------------------------------------------------------------
  IF (l_emcorr .AND. error_code  ==  0) THEN

    ! Add surface sensible heat flux into diabatic heating
    ! for use in energy correction procedure.

    CALL flux_diag(ftl, xx_cos_theta_latitude,                          &
              pdims%i_end,pdims%j_end,pdims_s%halo_i,pdims_s%halo_j,    &
                   1.0, sum_eng_fluxes, timestep)
    ! moisture flux level 1 held in fqw
    ! Should be total moisture flux from surface to layer 1 ie evaporation
    CALL flux_diag(fqw, xx_cos_theta_latitude,                          &
              pdims%i_end,pdims%j_end,pdims_s%halo_i,pdims_s%halo_j,    &
                   1.0, sum_moist_flux, timestep)

  END IF   ! L_emcorr
  ! ----------------------------------------------------------------------
  ! Section BL.6 Output Diagnostics requested.
  ! ----------------------------------------------------------------------
  SELECT CASE (model_type)
  CASE DEFAULT

    ! diagnostics requested this timestep
    IF ( sf(0,3) ) THEN
      IF ( error_code  ==  0) THEN

        ! Take rho values on the lowest level for
        ! visibility diagnostics
        ! Don't recalculate if it's already been calculated above
        IF ( .NOT. l_made_rho1 ) THEN
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              rho1(i,j) = rho_wet_rsq(i,j,1) /                          &
                        (r_rho_levels(i,j,1)*r_rho_levels(i,j,1))
            END DO
          END DO
          l_made_rho1 = .TRUE.
        END IF
#if !defined(LFRIC)
        CALL diagnostics_bl(                                            &
        ! IN levels / grids / switches
                             bl_levels, land_points,                    &
                             rhc_row_length, rhc_rows,                  &
                             land_index,ntiles,                         &
        ! IN fields for diagnostics
                             aerosol(1:tdims%i_end,1:tdims%j_end,1),    &
                             plsp, cca_at_base,                         &
                             ls_rain, ls_snow, conv_rain, conv_snow,    &
                             p_star, rhcpt, ntml, cumulus, rho1,        &
                             qcf_latest,                                &
                             zh,                                        &
                             e_sea, h_sea,ei, ustar_imp,                &
                             sea_ice_htf,                               &
                             zht,                                       &
                             bl_type_1,bl_type_2,bl_type_3,bl_type_4,   &
                             bl_type_5,bl_type_6,bl_type_7,             &
                             fqw, ftl, z0m_gb, z0m_eff_gb, z0h_eff_gb,  &
                             rib_gb, taux, tauy,                        &
                             t_soil, surf_ht_flux_gb,                   &
                             surf_ht_flux_land,surf_ht_flux_sice,       &
                             rib_ssi,ftl_ssi,e_ssi,ei_sice,             &
                             ftl_ice,                                   &
                             vshr_land,vshr_ssi,                        &
                             taux_land,taux_ssi,tauy_land,tauy_ssi,     &
                             radnet_sea,radnet_sice,                    &
            flandg(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),&
                             land_sea_mask,                             &
                             sil_orog_land,ho2r2_orog,Gs,gpp,npp,resp_p,&
                             ecan_tile,esoil_tile,gpp_ft,ftl_tile,      &
                             npp_ft,resp_p_ft,resp_s,resp_s_tot,cs,     &
                             rib_tile,es,ecan,fsmc,radnet_tile,sw_tile, &
                             tstar_tile,canopy,catch,z0m_tile,g_leaf,   &
                             le_tile,ei_tile,olr,                       &
                             epot_tile,tile_frac,                       &
                             co2_flux_tot, land_co2, dust_flux,         &
                             dust_emiss_frac,u_s_t_tile,u_s_t_dry_tile, &
                             u_s_std_tile,drydep2, BL_diag, sf_diag,    &
        ! variables required for soil moisture nudging scheme macro
                             rhokh,resfs,chr1p5m,alpha1,wt_ext,         &
                             lai_ft,canht_ft,gc,theta_star_surf,        &
                             qv_star_surf,                              &
        ! MGS extra bl vars for UKCA
              rhokh_mix, rho_aresist, aresist, resist_b, r_b_dust,      &
              dtrdz_charney_grid, kent, we_lim, t_frac, zrzi, kent_dsc, &
              we_lim_dsc, t_frac_dsc, zrzi_dsc, zhsc,                   &
        ! INOUT stash workspace
             stashwork3)
#endif
      END IF   ! error code zero
    END IF      ! sf(0,3)

  CASE (mt_single_column)

    ! Output diagnostics the SCM way.
    ! Initialise arrays
    lwp(:,:)         = 0.0
    iwp(:,:)         = 0.0
    bl_alltypes(:,:) = 0.0

    ! Set single diagnostic of BL type.
    DO j = pdims%j_start, pdims%j_end
      jScm = j - pdims%j_start + 1
      DO i = pdims%i_start, pdims%i_end
        iScm = i - pdims%i_start + 1
        bl_alltypes(iScm,jScm) = bl_type_1(i, j)*1.0                    &
                               + bl_type_2(i, j)*2.0                    &
                               + bl_type_3(i, j)*3.0                    &
                               + bl_type_4(i, j)*4.0                    &
                               + bl_type_5(i, j)*5.0                    &
                               + bl_type_6(i, j)*6.0                    &
                               + bl_type_7(i, j)*7.0
        rho1(i,j) = rho_wet_rsq(i,j,1) /                                &
              (r_rho_levels(i,j,1)*r_rho_levels(i,j,1))
      END DO ! i
    END DO ! j

    ! Calculate density averaged on theta levels, IWP and LWP
    DO k = 1, tdims%k_end-1
      DO j = tdims%j_start, tdims%j_end
        jScm = j - tdims%j_start + 1
        DO i = tdims%i_start, tdims%i_end
          iScm = i - tdims%i_start + 1

          rho_theta_levs = 0.5*                                         &
      ( rho_wet_rsq(i,j,k)/(r_rho_levels(i,j,k)*r_rho_levels(i,j,k))    &
      + rho_wet_rsq(i,j,k+1)                                            &
              /(r_rho_levels(i,j,k+1)*r_rho_levels(i,j,k+1)) )

          lwp(iScm,jScm) = lwp(iScm,jScm)                               &
                     + qcl_latest(i,j,k)*rho_theta_levs*                &
                       (r_rho_levels(i,j,k+1)-r_rho_levels(i,j,k))

          iwp(iScm,jScm) = iwp(iScm,jScm)                               &
                     + qcf_latest(i,j,k)*rho_theta_levs*                &
                       (r_rho_levels(i,j,k+1)-r_rho_levels(i,j,k))
        END DO ! i
      END DO ! j
    END DO ! k

    ! DEPENDS ON: dgnstcs_imp_ctl2
    CALL dgnstcs_imp_ctl2(                                              &
         ! IN
         pdims%i_end,pdims%j_end,rhc_row_length,rhc_rows,               &
         bl_levels, cloud_levels,                                       &
         rhcrit,combined_cloud,nclds,cumulus,ntml,                      &
         plsp,conv_rain,conv_snow,ls_rain,ls_snow,r_u,r_v,              &
         zh,zht,bl_type_1,bl_type_2,bl_type_3,bl_type_4,                &
         bl_type_5,bl_type_6,bl_type_7,bl_alltypes,                     &
         fqw,ftl,t_latest,q_latest,cf_latest,cfl_latest,cff_latest,     &
         t_earliest,q_earliest,qcl_earliest,qcf_earliest,               &
         cf_earliest,cfl_earliest,cff_earliest,                         &
         p_layer_centres(1, 1, 1),lwp,iwp,                              &
         z0m_gb,z0h_eff_gb,z0m_eff_gb,                                  &
         sea_ice_htf,taux,tauy,cf_area,                                 &
         p_star,cca_at_base,                                            &
         rho1,qcl_latest,qcf_latest,                                    &
         surf_ht_flux_sice,rhcpt,surf_ht_flux_gb,ustar_imp,aerosol,     &
         nSCMDpkgs,L_SCMDiags,BL_diag,sf_diag)

    ! SCM diags moved here from atmos_physics2
    IF (L_SCMDiags(scmdiag_land)) THEN

      CALL SCMoutput(ecan,'can_evap',                                   &
           'Canopy evaporation','kg/m2/day',                            &
           t_mult,d_sl,default_streams,'sec_day',RoutineName)

      CALL SCMoutput(le_tile,'lhf_tile',                                &
           'Tile latent heat flux','W/m2',                              &
           t_avg,d_tile,default_streams,'',RoutineName)

      CALL SCMoutput(es,'soil_evap',                                    &
           'Soil evapotranspiration','kg/m2/day',                       &
           t_mult,d_sl,default_streams,'sec_day',RoutineName)

      CALL SCMoutput(sf_diag%t1p5m_surft,'t1p5m_surft',                 &
           'Tile 1.5m temperature','K',                                 &
           t_avg,d_tile,default_streams,'',RoutineName)

    END IF ! L_SCMDiags(SCMDiag_land)

    IF (L_SCMDiags(scmdiag_surf)) THEN

      CALL SCMoutput(ei,'sublim',                                       &
           'Sublimation from lying snow or sea ice','kg/m2/day',        &
            t_mult,d_sl,default_streams,'sec_day',RoutineName)

    END IF ! L_SCMDiags(SCMDiag_surf)

  END SELECT ! model_type


  ! Clear up allocatable arrays
  IF (error_code  <=  0  .AND.  l_combi_cld) THEN
    DEALLOCATE ( cct4comb_cld   )
    DEALLOCATE ( ccb4comb_cld   )
    DEALLOCATE ( cca4comb_cld   )
  END IF

END IF ! L_apply_diag

DEALLOCATE(ustar_imp)
DEALLOCATE ( combined_cloud )
DEALLOCATE (BL_diag%cff_incr)
DEALLOCATE (BL_diag%cfl_incr)
DEALLOCATE (BL_diag%cf_incr)
! Deallocation of BL_diag%q, qcl and qcf _incr moved to atmos_phyiscs2 as
! required if moisture conservation being checked
DEALLOCATE (BL_diag%t_incr)
DEALLOCATE (BL_diag%v_incr)
DEALLOCATE (BL_diag%u_incr)
IF (l_subfilter_vert) DEALLOCATE (BL_diag%w_incr)

DEALLOCATE(sf_diag%u10m)
DEALLOCATE(sf_diag%v10m)
DEALLOCATE(sf_diag%latent_heat)
DEALLOCATE(sf_diag%lw_up_surft)
DEALLOCATE(sf_diag%lw_down_surft)
DEALLOCATE(sf_diag%t1p5m)
DEALLOCATE(sf_diag%q1p5m)
DEALLOCATE(sf_diag%t10m)
DEALLOCATE(sf_diag%q10m)
DEALLOCATE(sf_diag%sice_mlt_htf)
DEALLOCATE(sf_diag%snomlt_surf_htf)
DEALLOCATE(sf_diag%tstar_sice_weighted_cat)
DEALLOCATE(sf_diag%lw_up_sice_weighted_cat)
DEALLOCATE(sf_diag%ice_present_cat)
DEALLOCATE(sf_diag%tstar_sice_weighted)
DEALLOCATE(sf_diag%lw_up_sice_weighted)
DEALLOCATE(sf_diag%ice_present)
DEALLOCATE(sf_diag%ftl_ice_sm)

! Deallocate BL_diags on last substep
IF (l_apply_diag .OR. .NOT. l_quick_ap2) THEN
  DEALLOCATE(sf_diag%ch_ssi)
  DEALLOCATE(sf_diag%cd_ssi)
  DEALLOCATE(sf_diag%fme)
  DEALLOCATE(sf_diag%ra)
  DEALLOCATE(sf_diag%et_stom_ij)
  DEALLOCATE(sf_diag%et_stom_surft)
  DEALLOCATE(sf_diag%resfs_stom)
  DEALLOCATE(sf_diag%fprf)
  DEALLOCATE(sf_diag%fsth)
  DEALLOCATE(sf_diag%ftemp)
  DEALLOCATE(BL_diag%oblen)
  DEALLOCATE(BL_diag%ustar)
  DEALLOCATE(BL_diag%wbsurf)
  DEALLOCATE(BL_diag%gradrich)
  DEALLOCATE(BL_diag%wstar)
  DEALLOCATE(BL_diag%dbdz)
  DEALLOCATE(BL_diag%dvdzm)
  DEALLOCATE(BL_diag%rhokm)
  DEALLOCATE(BL_diag%rhokh)
  DEALLOCATE(BL_diag%tke)
  DEALLOCATE(BL_diag%ostressx)
  DEALLOCATE(BL_diag%ostressy)
  DEALLOCATE(sf_diag%u10m_n)
  DEALLOCATE(sf_diag%v10m_n)
  DEALLOCATE(sf_diag%mu10m_n)
  DEALLOCATE(sf_diag%mv10m_n)
  DEALLOCATE(sf_diag%chr10m)
  DEALLOCATE(BL_diag%smltop)
  DEALLOCATE(BL_diag%dsctop)
  DEALLOCATE(BL_diag%zhlocal)
  DEALLOCATE(BL_diag%zhpar)
  DEALLOCATE(BL_diag%dscbase)
  DEALLOCATE(BL_diag%cldbase)
  DEALLOCATE(BL_diag%weparm)
  DEALLOCATE(BL_diag%weparm_dsc)
  DEALLOCATE(BL_diag%dzh)
  DEALLOCATE(BL_diag%dTfric)
  DEALLOCATE(BL_diag%elm3d)
  DEALLOCATE(BL_diag%elh3d)
  DEALLOCATE(BL_diag%rhokmloc)
  DEALLOCATE(BL_diag%rhokhloc)
  DEALLOCATE(BL_diag%rhokmsurf)
  DEALLOCATE(BL_diag%rhokhsurf)
  DEALLOCATE(BL_diag%rhokmsc)
  DEALLOCATE(BL_diag%rhokhsc)
  DEALLOCATE(BL_diag%weight1d)
  DEALLOCATE(BL_diag%fm)
  DEALLOCATE(BL_diag%fh)
  DEALLOCATE(BL_diag%sgm_trb)
  DEALLOCATE(BL_diag%ql_trb)
  DEALLOCATE(BL_diag%cf_trb)
  DEALLOCATE(BL_diag%wb_ng)
  DEALLOCATE(BL_diag%sh)
  DEALLOCATE(BL_diag%sm)
  DEALLOCATE(BL_diag%tke_dissp)
  DEALLOCATE(BL_diag%tke_boy_prod)
  DEALLOCATE(BL_diag%tke_shr_prod)
  DEALLOCATE(BL_diag%elm)
  DEALLOCATE(BL_diag%rhogamq)
  DEALLOCATE(BL_diag%rhogamt)
  DEALLOCATE(BL_diag%rhogamv)
  DEALLOCATE(BL_diag%rhogamu)
END IF !l_apply_diag

DEALLOCATE (t_earliest)
DEALLOCATE (q_earliest)
DEALLOCATE (qcl_earliest)
DEALLOCATE (qcf_earliest)
DEALLOCATE (cf_earliest)
DEALLOCATE (cfl_earliest)
DEALLOCATE (cff_earliest)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ni_imp_ctl
END MODULE ni_imp_ctl_mod
