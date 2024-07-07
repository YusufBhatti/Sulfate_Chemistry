! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE SET_ATM_POINTERS ---------------------------------------
!
!    Set pointers for primary atmosphere fields
!
!   Programming Standard: Unified Model DP NO. 3, Version 8.3
!
!    Purpose:   Sets integer pointers to atmospheric
!               variables from STASHIN addresses.
!
!    External documentation: UMDP NO. C4
!
!
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Top Level

SUBROUTINE set_atm_pointers(                                      &
                  icode,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ukca_d1_defs, ONLY: ukca_item_sulpc, ukca_sect
USE um_stashcode_mod, ONLY: stashcode_glomap_clim_sec
USE clmchfcg_scenario_mod, ONLY: nsulpat
USE atm_fields_bounds_mod, ONLY: o3dims2, pdims, tdims_s,         &
                                 tdims, udims, vdims, wdims_s
USE run_aerosol_mod, ONLY: l_sulpc_online_oxidants
USE rad_input_mod, ONLY: lexpand_ozone, lexpand_tpps_ozone
USE lookup_addresses, ONLY: lbnpt
USE cv_run_mod, ONLY: l_3d_cca, l_ccrad
USE ukca_option_mod, ONLY: l_ukca
USE dump_headers_mod, ONLY: a_lookup
USE stash_array_mod, ONLY: ppindex, si
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE free_tracers_inputs_mod, ONLY: a_tracer_first, a_tracer_last, &
     a_tr_index, a_tr_stashitem, a_tr_lbc_stashitem,              &
     a_tr_active_lbc_index

USE jules_sea_seaice_mod, ONLY: nice_use

USE nlsizes_namelist_mod, ONLY: &
    a_len1_coldepc, a_len1_rowdepc, a_len2_coldepc, a_len2_levdepc,    &
    a_len2_rowdepc, land_field, len1_lookup,                           &
    len_dumphist, len_fixhd, model_levels,                             &
    mpp_len1_lookup, n_cca_lev, rows, sm_levels, st_levels,            &
    theta_field_size, theta_halo_size, theta_off_size,                 &
    tpps_ozone_levels, tr_lbc_ukca, tr_lbc_vars, tr_ukca, tr_vars,     &
    u_off_size, v_off_size

USE errormessagelength_mod, ONLY: errormessagelength

USE ukca_tracer_stash, ONLY:                                                   &
    a_ukca_last, ukca_tr_active_lbc_index, ukca_tr_lbc_stashitem,              &
    ukca_tr_stashitem, a_ukca_index, a_ukca_first

USE atm_d1_indices_mod, ONLY: jsst, jexnersurf, jdryrho, jetadot, jthetav,     &
    jpsiws, jpsiwl, jmv, jmcl, jmcf, jmcf2, jmrain, jmgraup, ju, jv, jw, jrho, &
    jtheta, jq, jqcl, jqcf, jqcf2, jqrain, jqgraup, jcloudnumber, jrainnumber, &
    jrain3mom, jicenumber, jsnownumber, jsnow3mom, jgraupnumber, jgraup3mom,   &
    jactivesolliquid, jactivesolrain, jactiveinsolice, jactivesolice,          &
    jactiveinsolliquid, jactivesolnumber, jactiveinsolnumber, je_trb, jtsq_trb,&
    jqsq_trb, jcov_trb, jzhpar_shcu, jdPV_rad, jdPV_sw, jdPV_lw, jdPV_mic,     &
    jdPV_gwd, jdPV_ph1, jdPV_conv, jdPV_bl, jdPV_stph, jdPV_cld, jdPV_iau,     &
    jdPV_nud, jdPV_tot, jdPV_adv, jdPV_sol, jdPV_mass, jadv_only_PV,           &
    jdtheta_0, jdtheta_bl,jdtheta_bl_mix, jdtheta_bl_LH, jdtheta_conv,         &
    jdtheta_mic, jdtheta_rad, jdtheta_SW, jdtheta_LW,                          &
    jdtheta_slow, jdtheta_cld, jexner_rho_levels,                              &
    jconv_prog_1, jconv_prog_2, jconv_prog_3, jconv_prog_precip, jtotalppn, jp,&
    jux_ccp, juy_ccp, jum_ccp, jg_ccp, jh_ccp, jriso_ccp, jrdir_ccp,           &
    jp_theta_levels, jexner_theta_levels, jccw_rad, jcca, jcca_dp, jcca_md,    &
    jcca_sh, jcf_area, jcf_bulk, jcf_liquid, jcf_frozen, j_deep_soil_temp,     &
    jsmcl, jsthu, jsthf, jsw_incs, jlw_incs, jdirpar, jozone, jtppsozone,      &
    jtracer, jtr_ukca, jmurk_source, jmurk, jflash_pot, jbl_w_var, jdust_div1, &
    jdust_div2, jdust_div3, jdust_div4, jdust_div5, jdust_div6, jso2, jdms,    &
    jso4_aitken, jso4_accu, jso4_diss, jh2o2, jnh3, jsoot_new, jsoot_agd,      &
    jsoot_cld, jbmass_new, jbmass_agd, jbmass_cld, jocff_new, jocff_agd,       &
    jocff_cld, jso2_natem, joh, jho2, jh2o2_limit, jo3_chem, jco2,             &
    joh_ukca, jho2_ukca, jh2o2_ukca, jo3_ukca, jhno3_ukca, jozone_tracer,      &
    jo3_prod_loss, jo3_p_l_vmr, jo3_vmr, jo3_p_l_temp, jo3_temp, jo3_p_l_colo3,&
    jo3_colo3, jarclbiog_bg, jarclbiom_fr, jarclbiom_ag, jarclbiom_ic,         &
    jarclblck_fr, jarclblck_ag, jarclsslt_fi, jarclsslt_jt, jarclsulp_ac,      &
    jarclsulp_ak, jarclsulp_di, jarcldust_b1, jarcldust_b2, jarcldust_b3,      &
    jarcldust_b4, jarcldust_b5, jarcldust_b6, jarclocff_fr, jarclocff_ag,      &
    jarclocff_ic, jarcldlta_dl, jnitr_acc, jnitr_diss, jgc_nd_nuc_sol,         &
    jgc_nuc_sol_su, jgc_nuc_sol_oc, jgc_nd_ait_sol, jgc_ait_sol_su,            &
    jgc_ait_sol_bc, jgc_ait_sol_oc, jgc_nd_acc_sol, jgc_acc_sol_su,            &
    jgc_acc_sol_bc, jgc_acc_sol_oc, jgc_acc_sol_ss, jgc_nd_cor_sol,            &
    jgc_cor_sol_su, jgc_cor_sol_bc, jgc_cor_sol_oc, jgc_cor_sol_ss,            &
    jgc_nd_ait_ins, jgc_ait_ins_bc, jgc_ait_ins_oc, juser_mult1, juser_mult2,  &
    juser_mult3, juser_mult4, juser_mult5, juser_mult6, juser_mult7,           &
    juser_mult8, juser_mult9, juser_mult10, juser_mult11, juser_mult12,        &
    juser_mult13, juser_mult14, juser_mult15, juser_mult16, juser_mult17,      &
    juser_mult18, juser_mult19, juser_mult20, jorog_lbc, ju_lbc, jv_lbc,       &
    jw_lbc, jrho_lbc, jtheta_lbc, jq_lbc, jqcl_lbc, jqcf_lbc, jqcf2_lbc,       &
    jqrain_lbc, jqgraup_lbc, jcf_bulk_lbc, jcf_liquid_lbc, jcf_frozen_lbc,     &
    jexner_lbc, ju_adv_lbc, jv_adv_lbc, jw_adv_lbc, jmurk_lbc, jtracer_lbc,    &
    jtr_ukca_lbc, jdust_div1_lbc, jdust_div2_lbc, jdust_div3_lbc,              &
    jdust_div4_lbc, jdust_div5_lbc, jdust_div6_lbc, jso2_lbc, jdms_lbc,        &
    jso4_aitken_lbc, jso4_accu_lbc, jso4_diss_lbc, jnh3_lbc, jsoot_new_lbc,    &
    jsoot_agd_lbc, jsoot_cld_lbc, jbmass_new_lbc, jbmass_agd_lbc,              &
    jbmass_cld_lbc, jocff_new_lbc, jocff_agd_lbc, jocff_cld_lbc, jnitr_acc_lbc,&
    jnitr_diss_lbc, ju_lbc_tend, jv_lbc_tend, jw_lbc_tend, jrho_lbc_tend,      &
    jtheta_lbc_tend, jq_lbc_tend, jqcl_lbc_tend, jqcf_lbc_tend, jqcf2_lbc_tend,&
    jqrain_lbc_tend, jqgraup_lbc_tend, jcf_bulk_lbc_tend, jcf_liquid_lbc_tend, &
    jcf_frozen_lbc_tend, jexner_lbc_tend, ju_adv_lbc_tend, jv_adv_lbc_tend,    &
    jw_adv_lbc_tend, jmurk_lbc_tend, jtracer_lbc_tend, jtr_ukca_lbc_tend,      &
    jdust_div1_lbc_tend, jdust_div2_lbc_tend, jdust_div3_lbc_tend,             &
    jdust_div4_lbc_tend, jdust_div5_lbc_tend, jdust_div6_lbc_tend,             &
    jso2_lbc_tend, jdms_lbc_tend, jso4_aitken_lbc_tend, jso4_accu_lbc_tend,    &
    jso4_diss_lbc_tend, jnh3_lbc_tend, jsoot_new_lbc_tend, jsoot_agd_lbc_tend, &
    jsoot_cld_lbc_tend, jbmass_new_lbc_tend, jbmass_agd_lbc_tend,              &
    jbmass_cld_lbc_tend, jocff_new_lbc_tend, jocff_agd_lbc_tend,               &
    jocff_cld_lbc_tend, jnitr_acc_lbc_tend, jnitr_diss_lbc_tend,               &
    jbl_pert_rand_fld, jbl_pert_flag, jtstar, jland,                           &
    jtstar_anom, jfrac_land, jtstar_land, jtstar_sea, jtstar_sice,             &
    jtstar_sice_cat, jsice_alb, jland_alb, jpstar, jlcbase, jccb, jcct, jcclwp,&
    jdeepflag, jpastprecip, jpastconvht, jzh, jddmfx, jt1_sd, jq1_sd,          &
    jtscrndcl_tile, jtscrndcl_ssi, jtstbtrans, jntml, jntdsc, jnbdsc, jcumulus,&
    jsat_soilw_suction, jtherm_cap, jtherm_cond, jvol_smc_crit, jvol_smc_wilt, &
    jvol_smc_sat, jsat_soil_cond, jclapp_horn, jz0m_soil, jcanopy_water, jz0,  &
    jgs, jorog, jorog_sd, jorog_sil,jorog_ho2, jorog_grad_x, jorog_grad_y,     &
    jorog_unfilt, jorog_grad_xx, jorog_grad_xy, jorog_grad_yy, ju_sea, jv_sea, &
    ju_0_p, jv_0_p, jice_fraction, jice_thickness, jti, jice_fract_cat,        &
    jice_thick_cat, jpond_frac_cat, jpond_depth_cat, jtfrz, jti_cat,           &
    jice_k_cat, jchloro_sea, jsnodep, jsnodep_sea, jsnodep_sea_cat,            &
    jcatch_snow, jsnow_grnd, jsnsoot, jsoil_clay, jsoil_silt, jsoil_sand,      &
    jdust_mrel1, jdust_mrel2, jdust_mrel3, jdust_mrel4, jdust_mrel5,           &
    jdust_mrel6, jso2_em, jdms_em, jso2_hilem, jnh3_em, jsoot_em, jsoot_hilem, &
    jbmass_em, jbmass_hilem, jocff_em, jocff_hilem, jdms_conc, jbmass_hilem_h1,&
    jbmass_hilem_h2, juser_anc1, juser_anc2, juser_anc3, juser_anc4,           &
    juser_anc5, juser_anc6, juser_anc7, juser_anc8, juser_anc9, juser_anc10,   &
    juser_anc11, juser_anc12, juser_anc13, juser_anc14, juser_anc15,           &
    juser_anc16, juser_anc17, juser_anc18, juser_anc19, juser_anc20, jnet_flux,&
    jnet_mflux, jfrac_typ, jfrac_con1, jfrac_con2, jfrac_con3, jfrac_con4,     &
    jfrac_con5, jfrac_con6, jfrac_con7, jfrac_con8, jfrac_con9, jlai_pft,      &
    jcanht_pft, jdisturb, jdisturb_prev, jwoodprod_fast, jwoodprod_med,        &
    jwoodprod_slow, jsoil_alb, jobs_alb_sw, jobs_alb_vis, jobs_alb_nir,        &
    jsoil_carb, jsoil_carb1, jsoil_carb2, jsoil_carb3, jsoil_carb4,            &
    jnpp_pft_acc, jg_lf_pft_acc, jg_phlf_pft_acc, jrsp_w_pft_acc, jrsp_s_acc,  &
    jrsp_s_acc1, jrsp_s_acc2, jrsp_s_acc3, jrsp_s_acc4, jcan_water_tile,       &
    jcatch_tile, jinfil_tile, jrgrain_tile, jsnodep_tile, jtstar_tile,         &
    jtsurf_elev_surft,                                                         &
    jz0_tile, jz0h_tile, jdolr, jlw_down, jsw_tile, jurbhgt, jurbhwr, jurbwrr, &
    jurbdisp, jurbztm, jurbalbwl, jurbalbrd, jurbemisw, jurbemisr, j_co2flux,  &
    j_co2_emits, jetatheta, jetarho, jrhcrit, jsoil_thickness,                 &
    jzseak_theta, jck_theta, jzseak_rho, jck_rho, jlambda_input_p,             &
    jlambda_input_u, jphi_input_p, jphi_input_v, jti_mean, jti_sig, jfexp,     &
    jgamma_int, jwater_table, jfsfc_sat, jf_wetland, jsthzw, ja_fsat, jc_fsat, &
    ja_fwet, jc_fwet, jriv_sequence, jriv_direction, jriv_storage,             &
    jtot_surfroff, jtot_subroff, jriv_inlandatm, jacc_lake_evap, jc_solar,     &
    jc_blue, jc_longwave, jc_taux, jc_tauy, jc_w10, jc_sensible, jc_sublim,    &
    jc_evap, jc_fcondtopn, jc_topmeltn, jc_tstar_sicen, jc_lsrain, jc_lssnow,  &
    jc_cvrain, jc_cvsnow, jc_riverout, jc_mslp, jc_calving, jc_surf_co2,       &
    jc_dust_dep, jsnowdepth, jrho_snow_grnd, jnsnow, jds, jsice, jsliq,        &
    jtsnowlayer, jrho_snow, jrgrainl, jlake_depth, jlake_fetch, jlake_t_mean,  &
    jlake_t_mxl, jlake_t_ice, jlake_h_mxl, jlake_h_ice, jlake_shape,           &
    jlake_g_dt, jsoil_nitro1, jsoil_nitro2, jsoil_nitro3, jsoil_nitro4,        &
    jsoil_inorgnit, j_n_deposition,                                            &
    jpasture, jpasture_prev, jagr_crop, jagr_crop_prev, j_triffid_co2


IMPLICIT NONE

!
!   Arguments
!
INTEGER ::                                                        &
    icode                  ! OUT: Error return code
!
CHARACTER(LEN=errormessagelength) ::                              &
    cmessage               ! OUT: Error return message

! local variables

INTEGER ::                                                        &
        ivar,                                                     &
                           ! Loop counts
        jvar,                                                     &
                           ! Loop counts
        lev                                                       &
       ,im_index                                                  &
                      !  Internal Model Index in Stash arrays
       ,Sect_No       !  Stash section number

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SET_ATM_POINTERS'

!     Set to atmosphere internal model
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
im_index  = 1
Sect_No   = 0


! Set pointers for atmospheric primary variables from STASH :-

ju(udims%k_start)                = si(  2,Sect_No,im_index)
jv(vdims%k_start)                = si(  3,Sect_No,im_index)
jtheta(tdims%k_start)            = si(  4,Sect_No,im_index)
jq    (tdims%k_start)            = si( 10,Sect_No,im_index)
jqcf  (tdims%k_start)            = si( 12,Sect_No,im_index)
jtstar                           = si( 24,Sect_No,im_index)
jland                            = si( 30,Sect_No,im_index)
jorog                            = si( 33,Sect_No,im_index)
jw(wdims_s%k_start)              = si(150,Sect_No,im_index)
jrho   (pdims%k_start)           = si(253,Sect_No,im_index)
jqcl   (tdims%k_start)           = si(254,Sect_No,im_index)
jexner_rho_levels(pdims%k_start) = si(255,Sect_No,im_index)
jqcf2  (tdims%k_start)           = si(271,Sect_No,im_index)
jqrain (tdims%k_start)           = si(272,Sect_No,im_index)
jqgraup(tdims%k_start)           = si(273,Sect_No,im_index)

!     ENDGame prognostics (idealised)
jsst                             = si(94,Sect_No,im_index)

!     ENDGame prognostics
jthetav(tdims%k_start)           = si(388,Sect_No,im_index)
jdryrho(pdims%k_start)           = si(389,Sect_No,im_index)
jetadot(tdims%k_start)           = si(387,Sect_No,im_index)
jpsiws                           = si(390,Sect_No,im_index)
jpsiwl                           = si(397,Sect_No,im_index)
jexnersurf                       = si(398,Sect_No,im_index)
jmv(tdims%k_start)               = si(391,Sect_No,im_index)
jmcl(tdims%k_start)              = si(392,Sect_No,im_index)
jmcf(tdims%k_start)              = si(393,Sect_No,im_index)
jmcf2(tdims%k_start)             = si(396,Sect_No,im_index)
jmrain(tdims%k_start)            = si(394,Sect_No,im_index)
jmgraup(tdims%k_start)           = si(395,Sect_No,im_index)

!     for TKE based turbulence scheme
je_trb  (tdims%k_start)          = si(70,Sect_No,im_index)
jtsq_trb(tdims%k_start)          = si(71,Sect_No,im_index)
jqsq_trb(tdims%k_start)          = si(72,Sect_No,im_index)
jcov_trb(tdims%k_start)          = si(73,Sect_No,im_index)
jzhpar_shcu                      = si(74,Sect_No,im_index)

! CASIM microphysics prognostics
jcloudnumber(tdims%k_start)      = si(75,Sect_No,im_index)
jrainnumber (tdims%k_start)      = si(76,Sect_No,im_index)
jrain3mom   (tdims%k_start)      = si(77,Sect_No,im_index)
jicenumber  (tdims%k_start)      = si(78,Sect_No,im_index)
jsnownumber (tdims%k_start)      = si(79,Sect_No,im_index)
jsnow3mom   (tdims%k_start)      = si(80,Sect_No,im_index)
jgraupnumber(tdims%k_start)      = si(81,Sect_No,im_index)
jgraup3mom  (tdims%k_start)      = si(82,Sect_No,im_index)

! CASIM activated aerosol prognostics 
jactivesolliquid  (tdims%k_start) = si( 83,Sect_No,im_index)
jactivesolrain    (tdims%k_start) = si( 84,Sect_No,im_index)
jactiveinsolice   (tdims%k_start) = si( 85,Sect_No,im_index)
jactivesolice     (tdims%k_start) = si( 86,Sect_No,im_index)
jactiveinsolliquid(tdims%k_start) = si( 87,Sect_No,im_index)
jactivesolnumber  (tdims%k_start) = si( 88,Sect_No,im_index)
jActiveinsolnumber(tdims%k_start) = si( 89,Sect_No,im_index)

!     for lightning/electric scheme
jflash_pot(tdims%k_start)        = si(91,Sect_No,im_index)

!     for turbulent qcl production in mixed phase cloud 
jbl_w_var(1)                     = si(99,Sect_No,im_index)

!Stochastic physics prognostics
! used in BL perturbations
jbl_pert_rand_fld = si(595,Sect_No,im_index)
jbl_pert_flag     = si(596,Sect_No,im_index)

! Set extra pointers for coastal tiling.
jfrac_land     = si(505,Sect_No,im_index)
jtstar_land    = si(506,Sect_No,im_index)
jtstar_sea     = si(507,Sect_No,im_index)
jtstar_sice    = si(508,Sect_No,im_index)
jtstar_sice_cat= si(441,Sect_No,im_index)
! Set pointers for seaice and land albedos
jsice_alb      = si(509,Sect_No,im_index)
jland_alb      = si(510,Sect_No,im_index)
!
! Set extra pointers for large-scale hydrology.
jti_mean       = si(274,Sect_No,im_index)
jti_sig        = si(275,Sect_No,im_index)
jfexp          = si(276,Sect_No,im_index)
jgamma_int     = si(277,Sect_No,im_index)
jwater_table   = si(278,Sect_No,im_index)
jfsfc_sat      = si(279,Sect_No,im_index)
jf_wetland     = si(280,Sect_No,im_index)
jsthzw         = si(281,Sect_No,im_index)
ja_fsat        = si(282,Sect_No,im_index)
jc_fsat        = si(283,Sect_No,im_index)
ja_fwet        = si(284,Sect_No,im_index)
jc_fwet        = si(285,Sect_No,im_index)
jacc_lake_evap = si(290,Sect_No,im_index)


DO lev= 1, wdims_s%k_end
  jw(lev)=jw(lev-1)+theta_off_size
END DO

DO lev= udims%k_start+1, udims%k_end
  ju(lev)    =ju(lev-1)     + u_off_size
END DO

DO lev= vdims%k_start+1, vdims%k_end
  jv(lev)    =jv(lev-1)     + v_off_size
END DO

DO lev= pdims%k_start+1, pdims%k_end
  jrho(lev)    =jrho(lev-1)      + theta_off_size
  jdryrho(lev) =jdryrho(lev-1)   + theta_off_size
END DO

DO lev= tdims%k_start+1, tdims%k_end
  jtheta(lev) =jtheta(lev-1)  + theta_off_size
  jthetav(lev)=jthetav(lev-1) + theta_off_size
  jetadot(lev)=jetadot(lev-1) + theta_off_size
END DO

DO lev= pdims%k_start+1, pdims%k_end+1
  jexner_rho_levels(lev)=jexner_rho_levels(lev-1)+theta_off_size
END DO

DO lev= tdims%k_start+1, tdims%k_end
  jq(lev)    =jq(lev-1)     + theta_halo_size
  jqcl(lev)  =jqcl(lev-1)   + theta_halo_size
  jqcf(lev)  =jqcf(lev-1)   + theta_halo_size
  jqcf2(lev)   =jqcf2(lev-1)   + theta_halo_size
  jqrain(lev)  =jqrain(lev-1)  + theta_halo_size
  jqgraup(lev) =jqgraup(lev-1) + theta_halo_size
  jmv(lev)     =jmv(lev-1)     + theta_off_size
  jmcl(lev)    =jmcl(lev-1)    + theta_off_size
  jmcf(lev)    =jmcf(lev-1)    + theta_off_size
  jmcf2(lev)   =jmcf2(lev-1)   + theta_off_size
  jmrain(lev)  =jmrain(lev-1)  + theta_off_size
  jmgraup(lev) =jmgraup(lev-1) + theta_off_size

  jcloudnumber(lev) = jcloudnumber(lev-1) + theta_halo_size
  jrainnumber(lev)  = jrainnumber(lev-1)  + theta_halo_size
  jrain3mom(lev)    = jrain3mom(lev-1)    + theta_halo_size
  jicenumber(lev)   = jicenumber(lev-1)   + theta_halo_size
  jsnownumber(lev)  = jsnownumber(lev-1)  + theta_halo_size
  jsnow3mom(lev)    = jsnow3mom(lev-1)    + theta_halo_size
  jgraupnumber(lev) = jgraupnumber(lev-1) + theta_halo_size
  jgraup3mom(lev)   = jgraup3mom(lev-1)   + theta_halo_size

  jactivesolliquid(lev)   = jactivesolliquid(lev-1)   + theta_halo_size
  jactivesolrain(lev)     = jactivesolrain(lev-1)     + theta_halo_size
  jactiveinsolice(lev)    = jactiveinsolice(lev-1)    + theta_halo_size
  jactivesolice(lev)      = jactivesolice(lev-1)      + theta_halo_size
  jactiveinsolliquid(lev) = jactiveinsolliquid(lev-1) + theta_halo_size
  jactivesolnumber(lev)   = jactivesolnumber(lev-1)   + theta_halo_size
  jactiveinsolnumber(lev) = jactiveinsolnumber(lev-1) + theta_halo_size
END DO

! Lightning prognostics
DO lev = tdims%k_start+1, tdims%k_end
  jflash_pot(lev) = jflash_pot(lev-1) + theta_field_size
END DO

DO lev= tdims%k_start+1, tdims%k_end
  je_trb(lev) = je_trb(lev-1) + theta_field_size
  jtsq_trb(lev) = jtsq_trb(lev-1) + theta_field_size
  jqsq_trb(lev) = jqsq_trb(lev-1) + theta_field_size
  jcov_trb(lev) = jcov_trb(lev-1) + theta_field_size
END DO

! Turbulent production of qcl in mixed phase cloud 
DO lev = 2, tdims%k_end 
  jbl_w_var(lev) = jbl_w_var(lev-1) + theta_field_size
END DO
!
!
!     Convective downdraught mass-flux at cloud base
jddmfx         = si(493,Sect_No,im_index)

! Pointers required to save coupling fields as
! prognostics when employing OASIS as the coupler.
! In non-coupled models these pointers will
! simply end up with a value of 1.
jc_solar = si(171,Sect_No,im_index)
jc_blue =  si(172,Sect_No,im_index)
jc_longwave =  si(174,Sect_No,im_index)
jc_taux =  si(176,Sect_No,im_index)
jc_tauy =  si(177,Sect_No,im_index)
jc_sensible =  si(179,Sect_No,im_index)
IF (nice_use  ==  1) THEN
  jc_sublim =  si(180,Sect_No,im_index)  ! Single cat field
ELSE
  jc_sublim =  si(182,Sect_No,im_index)  ! Multi cat field
END IF
jc_evap =  si(181,Sect_No,im_index)
jc_fcondtopn =  si(184,Sect_No,im_index)
jc_topmeltn =  si(185,Sect_No,im_index)
jc_lsrain =  si(186,Sect_No,im_index)
jc_lssnow =  si(187,Sect_No,im_index)
jc_cvrain =  si(188,Sect_No,im_index)
jc_cvsnow =  si(189,Sect_No,im_index)
jc_calving =  si(190,Sect_No,im_index)
jc_w10 =  si(191,Sect_No,im_index)
jc_riverout =  si(192,Sect_No,im_index)
jc_mslp =  si(193,Sect_No,im_index)
jc_tstar_sicen = si(195,Sect_No,im_index)
jc_surf_CO2 = si(196,Sect_No,im_index)
jc_dust_dep = si(197,Sect_No,im_index)

! Set pointers for optional atmospheric primary variables.
jzh                         = si( 25,Sect_No,im_index)
jntml                       = si(259,Sect_No,im_index)
jnbdsc                      = si(260,Sect_No,im_index)
jntdsc                      = si(261,Sect_No,im_index)
jcumulus                    = si(262,Sect_No,im_index)
jt1_sd                      = si(263,Sect_No,im_index)
jq1_sd                      = si(264,Sect_No,im_index)
jcf_area  (1)               = si(265,Sect_No,im_index)
jcf_bulk  (tdims%k_start)   = si(266,Sect_No,im_index)
jcf_liquid(tdims%k_start)   = si(267,Sect_No,im_index)
jcf_frozen(tdims%k_start)   = si(268,Sect_No,im_index)

IF (l_3d_cca .OR. l_ccrad) THEN
  jcca(1) = si(211,Sect_No,im_index)
  DO lev=2, n_cca_lev      ! n_cca_lev set in dervsize
    jcca(lev)=jcca(lev-1)+theta_field_size
  END DO
ELSE
  jcca(1)      = si( 13,Sect_No,im_index)
END IF

jcca_dp(1) = si(345,Sect_No,im_index)
jcca_md(1) = si(346,Sect_No,im_index)
jcca_sh(1) = si(347,Sect_No,im_index)

DO lev = 2, n_cca_lev
  jcca_dp(lev) = jcca_dp(lev-1)+theta_field_size
  jcca_md(lev) = jcca_md(lev-1)+theta_field_size
  jcca_sh(lev) = jcca_sh(lev-1)+theta_field_size
END DO

jccb           = si( 14,Sect_No,im_index)
jcct           = si( 15,Sect_No,im_index)
jcclwp         = si( 16,Sect_No,im_index)
jdeepflag      = si(342,Sect_No,im_index)
jpastprecip    = si(343,Sect_No,im_index)
jpastconvht    = si(344,Sect_No,im_index)
jlcbase        = si( 21,Sect_No,im_index)
jcanopy_water  = si( 22,Sect_No,im_index)
jccw_rad(1)    = si(212,Sect_No,im_index)

DO lev= 2, tdims%k_end
  jccw_rad(lev) =jccw_rad(lev-1)+theta_field_size
  jcf_area(lev) =jcf_area(lev-1)+theta_field_size
END DO

DO lev= tdims%k_start+1, tdims%k_end
  jcf_bulk(lev)  =jcf_bulk(lev-1)+THETA_halo_SIZE
  jcf_liquid(lev)=jcf_liquid(lev-1)+THETA_halo_SIZE
  jcf_frozen(lev)=jcf_frozen(lev-1)+THETA_halo_SIZE
END DO

jconv_prog_1(tdims%k_start)         = si(136,Sect_No,im_index)
jconv_prog_2(tdims%k_start)         = si(137,Sect_No,im_index)
jconv_prog_3(tdims%k_start)         = si(138,Sect_No,im_index)
jconv_prog_precip(tdims%k_start)    = si(139,Sect_No,im_index)
DO lev= tdims%k_start+1, tdims%k_end
  jconv_prog_1(lev)      = jconv_prog_1(lev-1)      + theta_off_size
  jconv_prog_2(lev)      = jconv_prog_2(lev-1)      + theta_off_size
  jconv_prog_3(lev)      = jconv_prog_3(lev-1)      + theta_off_size
  jconv_prog_precip(lev) = jconv_prog_precip(lev-1) + theta_off_size
END DO

jtotalppn  = si(348,Sect_No,im_index)

!  Set pointers for secondary fields in D1
jexner_theta_levels(tdims%k_start) = si(406,Sect_No,im_index)
jp                 (pdims%k_start) = si(407,Sect_No,im_index)
jp_theta_levels    (tdims%k_start) = si(408,Sect_No,im_index)
jpstar                             = si(409,Sect_No,im_index)
jsw_incs(0)                        = si(410,Sect_No,im_index)
jlw_incs(0)                        = si(411,Sect_No,im_index)

! This code should work for both ND and V-AT-POLES
DO lev= 1, model_levels+1
  jsw_incs(lev)=jsw_incs(lev-1)+theta_field_size
END DO
DO lev= 1, model_levels
  jlw_incs(lev)=jlw_incs(lev-1)+theta_field_size
END DO

! Direct PAR flux
jdirpar               = si(460,Sect_no,im_index)

DO lev= tdims%k_start+1, tdims%k_end
  jexner_theta_levels(lev)=jexner_theta_levels(lev-1)+            &
                             theta_off_size
  jp_theta_levels(lev)    =jp_theta_levels(lev-1)+theta_off_size
END DO

DO lev= pdims%k_start+1, pdims%k_end+1
  jp(lev)                 =jp(lev-1)+theta_off_size
END DO

!  Set pointers for ancillary fields in D1 from STASH
!     Soil fields
jsmcl(1)            = si(  9,Sect_No,im_index)
j_deep_soil_temp(1) = si( 20,Sect_No,im_index)
jvol_smc_wilt       = si( 40,Sect_No,im_index)
jvol_smc_crit       = si( 41,Sect_No,im_index)
jvol_smc_sat        = si( 43,Sect_No,im_index)
jsat_soil_cond      = si( 44,Sect_No,im_index)
jtherm_cap          = si( 46,Sect_No,im_index)
jtherm_cond         = si( 47,Sect_No,im_index)
jsat_soilw_suction  = si( 48,Sect_No,im_index)
jclapp_horn         = si(207,Sect_No,im_index)
jsthu(1)            = si(214,Sect_No,im_index)
jsthf(1)            = si(215,Sect_No,im_index)
jz0m_soil           = si( 97,Sect_No,im_index)

DO lev=2,st_levels
  j_deep_soil_temp(lev)=j_deep_soil_temp(lev-1)+land_field
END DO

DO lev=2,sm_levels
  jsmcl(lev)=jsmcl(lev-1)+land_field
  jsthu(lev)=jsthu(lev-1)+land_field
  jsthf(lev)=jsthf(lev-1)+land_field
END DO

!     Other surface fields
jz0          = si( 26,Sect_No,im_index) ! roughness length (used for sea)
jgs          = si(213,Sect_No,im_index) ! stomatal conductance

! Orography fields
jorog_sil      = si(17,Sect_No,im_index)   ! Silhouette area
jorog_ho2      = si(18,Sect_No,im_index)   ! Peak to trough ht.
jorog_sd       = si(34,Sect_No,im_index)
jorog_grad_x   = si( 5,Sect_No,im_index)
jorog_grad_y   = si( 6,Sect_No,im_index)
jorog_unfilt   = si( 7,Sect_No,im_index)
jorog_grad_xx  = si(35,Sect_No,im_index)
jorog_grad_xy  = si(36,Sect_No,im_index)
jorog_grad_yy  = si(37,Sect_No,im_index)

! Sea/Sea Ice fields
ju_sea         = si( 28,Sect_No,im_index)
jv_sea         = si( 29,Sect_No,im_index)
jice_fraction  = si( 31,Sect_No,im_index)
jice_thickness = si( 32,Sect_No,im_index)
jti            = si( 49,Sect_No,im_index)
jice_fract_cat = si(413,Sect_No,im_index)
jice_thick_cat = si(414,Sect_No,im_index)
jpond_frac_cat = si(428,Sect_No,im_index) 
jpond_depth_cat = si(429,Sect_No,im_index) 
jtfrz          = si(194,Sect_No,im_index) 
jti_cat        = si(415,Sect_No,im_index)
jice_k_cat     = si(440,Sect_No,im_index)
ju_0_p         = si(269,Sect_No,im_index)
jv_0_p         = si(270,Sect_No,im_index)
jchloro_sea    = si( 96,Sect_No,im_index)

! Snow fields
jsnodep        = si( 23,Sect_No,im_index) ! Snow depth over land
jsnodep_sea    = si( 95,Sect_No,im_index) ! Snow depth on sea ice
jsnodep_sea_cat= si(416,Sect_No,im_index) ! Snow depth on ice cats
jsnsoot        = si(221,Sect_No,im_index) ! Snow soot content
jcatch_snow    = si(241,Sect_No,im_index)
jsnow_grnd     = si(242,Sect_No,im_index)

! Decoupled screen temperatures
JTScrnDcl_TILE = si(490,Sect_No,im_index) ! Decoupled screen-level
                                          ! temperature on tiles
JTScrnDcl_SSI  = si(491,Sect_No,im_index) ! Decoupled screen-level
                                          ! temperature on sea/s.ice
JtStbTrans     = si(492,Sect_No,im_index) ! Time since the transition

! Convective cold pools
jux_ccp   = si(569,Sect_No,im_index) ! x cmpt front speed vector sum
juy_ccp   = si(570,Sect_No,im_index) ! y cmpt front speed vector sum
jum_ccp   = si(571,Sect_No,im_index) ! front speed scalar sum
jg_ccp    = si(572,Sect_No,im_index) ! gridbox c.c.p. reduced gravity
jh_ccp    = si(573,Sect_No,im_index) ! gridbox c.c.p. depth
jriso_ccp = si(574,Sect_No,im_index) ! remain counter (isotropic)
jrdir_ccp = si(575,Sect_No,im_index) ! remain counter (directed)

! Ozone
jozone(o3dims2%k_start)     = si(60,Sect_No,im_index)
! Check for zonal ozone and calculate pointers accordingly
lexpand_ozone=.FALSE.
IF (a_lookup(lbnpt,ppindex(60,im_index)) == 1) THEN
  lexpand_ozone = .TRUE.
END IF

DO lev= o3dims2%k_start+1, o3dims2%k_end
  IF (lexpand_ozone) THEN
    !         Ozone held as zonal averages, i.e. one value per row
    jozone(lev)=jozone(lev-1)+rows
  ELSE
    jozone(lev)=jozone(lev-1)+theta_field_size
  END IF
END DO

!! Tropopause-based Ozone
IF (tpps_ozone_levels >  0) THEN
  jtppsozone(1)     = si(341,Sect_No,im_index)

  !Check for zonal tpps_ozone and calculate pointers accordingly
  lexpand_tpps_ozone=.FALSE.
  IF (a_lookup(lbnpt,ppindex(341,im_index)) == 1) THEN
    lexpand_tpps_ozone = .TRUE.
  END IF
END IF

DO lev=2,tpps_ozone_levels
  IF (lexpand_tpps_ozone) THEN
    jtppsozone(lev)=jtppsozone(lev-1)+rows
  ELSE
    jtppsozone(lev)=jtppsozone(lev-1)+theta_field_size
  END IF
END DO

! Add prognostic ozone tracer and cariolle parameters to section 0
!
jozone_tracer = si(480,sect_no,im_index)
jo3_prod_loss = si(481,sect_no,im_index)
jo3_p_l_vmr   = si(482,sect_no,im_index)
jo3_vmr       = si(483,sect_no,im_index)
jo3_p_l_temp  = si(484,sect_no,im_index)
jo3_temp      = si(485,sect_no,im_index)
jo3_p_l_colo3 = si(486,sect_no,im_index)
jo3_colo3     = si(487,sect_no,im_index)

DO  lev = tdims%k_start+1, tdims%k_end
  jozone_tracer(lev) = jozone_tracer(lev-1) + theta_off_size
  jo3_prod_loss(lev) = jo3_prod_loss(lev-1) + rows
  jo3_p_l_vmr(lev)   = jo3_p_l_vmr(lev-1) + rows
  jo3_vmr(lev)       = jo3_vmr(lev-1) + rows
  jo3_p_l_temp(lev)  = jo3_p_l_temp(lev-1) + rows
  jo3_temp(lev)      = jo3_temp(lev-1) + rows
  jo3_p_l_colo3(lev) = jo3_p_l_colo3(lev-1) + rows
  jo3_colo3(lev)     = jo3_colo3(lev-1) + rows
END DO


! STOCHEM fields - removed

! Add sources and aerosol ancillaries

jmurk_source(tdims%k_start)= si(57,Sect_No,im_index) !Murk source
jso2_em                    = si(58,Sect_No,im_index) !Sulphur dioxide emiss.
jdms_em                    = si(59,Sect_No,im_index) !Dimethyl sulphide emiss.
jmurk       (tdims%k_start)= si(90,Sect_No,im_index) !Murk concentration

! Add for Sulphur Cycle
jso2       (tdims%k_start)= si(101,Sect_No,im_index) !Sulphur dioxide gas
jdms       (tdims%k_start)= si(102,Sect_No,im_index) !Dimethyl sulphide gas
jso4_aitken(tdims%k_start)= si(103,Sect_No,im_index) !Aitken mode SO4 aerosol
jso4_accu  (tdims%k_start)= si(104,Sect_No,im_index) !Accumulation mode SO4 aer
jso4_diss  (tdims%k_start)= si(105,Sect_No,im_index) !Dissolved SO4 aerosol
jh2o2      (tdims%k_start)= si(106,Sect_No,im_index) !Hydrogen peroxide mmr
jnh3       (tdims%k_start)= si(107,Sect_No,im_index) !Ammonia gas
jsoot_new  (tdims%k_start)= si(108,Sect_No,im_index) !Fresh soot
jsoot_agd  (tdims%k_start)= si(109,Sect_No,im_index) !Aged soot
jsoot_cld  (tdims%k_start)= si(110,Sect_No,im_index) !Soot in cloud
jbmass_new (tdims%k_start)= si(111,Sect_No,im_index) !Fresh biomass smoke
jbmass_agd (tdims%k_start)= si(112,Sect_No,im_index) !Aged biomass smoke
jbmass_cld (tdims%k_start)= si(113,Sect_No,im_index) !Cloud biomass smoke
jocff_new  (tdims%k_start)= si(114,Sect_No,im_index) !Fresh ocff
jocff_agd  (tdims%k_start)= si(115,Sect_No,im_index) !Aged socff
jocff_cld  (tdims%k_start)= si(116,Sect_No,im_index) !Ocff in cloud
jbmass_hilem_h1 =si(119,Sect_No,im_index)  !Min height for elevated bmass emiss
jbmass_hilem_h2 =si(120,Sect_No,im_index)  !Max height for elevated bmass emiss
jso2_natem (tdims%k_start)= si(121,Sect_No,im_index) !Natural SO2 emissions
joh        (tdims%k_start)= si(122,Sect_No,im_index) !OH 3_D ancillary
jho2       (tdims%k_start)= si(123,Sect_No,im_index) !HO2 3_D ancillary
jh2o2_limit(tdims%k_start)= si(124,Sect_No,im_index) !H2O2 LIMIT 3_D ancillary
jo3_chem   (tdims%k_start)= si(125,Sect_No,im_index) !O3 for chemistry 3_D anc
jso2_hilem    =si(126,Sect_No,im_index)  !High level SO2 emissions
jnh3_em       =si(127,Sect_No,im_index)  !Ammonia surface emiss
jsoot_em      =si(128,Sect_No,im_index)  !Fresh soot surf emiss
jsoot_hilem   =si(129,Sect_No,im_index)  !Fresh soot high emiss
jbmass_em     =si(130,Sect_No,im_index)  !Fresh bmass surf emiss
jbmass_hilem  =si(131,Sect_No,im_index)  !Elevated bmass emiss
jdms_conc     =si(132,Sect_No,im_index)  !DMS conc in seawater
jocff_em      =si(134,Sect_No,im_index)  !Fresh OCFF surf emiss
jocff_hilem   =si(135,Sect_No,im_index)  !Fresh OCFF high emiss

! Aerosol climatologies
jarclbiog_bg  =si(351,Sect_No,im_index)  ! Biogenic aerosol climatology
jarclbiom_fr  =si(352,Sect_No,im_index)  ! Biomass burning (fresh) aerosol clim
jarclbiom_ag  =si(353,Sect_No,im_index)  ! Biomass burning (aged) aerosol clim
jarclbiom_ic  =si(354,Sect_No,im_index)  ! Biomass burning (in-cloud) aerosol clim
jarclblck_fr  =si(355,Sect_No,im_index)  ! Black carbon (fresh) aerosol clim
jarclblck_ag  =si(356,Sect_No,im_index)  ! Black carbon (aged) aerosol clim
jarclsslt_fi  =si(357,Sect_No,im_index)  ! Sea salt (film mode) aerosol clim
jarclsslt_jt  =si(358,Sect_No,im_index)  ! Sea salt (jet mode) aerosol clim
jarclsulp_ac  =si(359,Sect_No,im_index)  ! Sulphate (accumulation mode) aero clim
jarclsulp_ak  =si(360,Sect_No,im_index)  ! Sulphate (Aitken mode) aerosol clim
jarclsulp_di  =si(361,Sect_No,im_index)  ! Sulphate (dissolved) aerosol clim
jarcldust_b1  =si(362,Sect_No,im_index)  ! Dust (bin 1) aerosol climatology
jarcldust_b2  =si(363,Sect_No,im_index)  ! Dust (bin 2) aerosol climatology
jarcldust_b3  =si(364,Sect_No,im_index)  ! Dust (bin 3) aerosol climatology
jarcldust_b4  =si(365,Sect_No,im_index)  ! Dust (bin 4) aerosol climatology
jarcldust_b5  =si(366,Sect_No,im_index)  ! Dust (bin 5) aerosol climatology
jarcldust_b6  =si(367,Sect_No,im_index)  ! Dust (bin 6) aerosol climatology
jarclocff_fr  =si(368,Sect_No,im_index)  ! Org carbon fossil fuel (fresh) aero clim
jarclocff_ag  =si(369,Sect_No,im_index)  ! Org carbon fossil fuel (aged) aero clim
jarclocff_ic  =si(370,Sect_No,im_index)  ! Org carbon fossil fuel (in-cloud) aero clim
jarcldlta_dl  =si(371,Sect_No,im_index)  ! Delta aerosol climatology

DO lev= tdims%k_start+1, tdims%k_end
  jarclbiog_bg(lev) = jarclbiog_bg(lev-1)+theta_field_size
  jarclbiom_fr(lev) = jarclbiom_fr(lev-1)+theta_field_size
  jarclbiom_ag(lev) = jarclbiom_ag(lev-1)+theta_field_size
  jarclbiom_ic(lev) = jarclbiom_ic(lev-1)+theta_field_size
  jarclblck_fr(lev) = jarclblck_fr(lev-1)+theta_field_size
  jarclblck_ag(lev) = jarclblck_ag(lev-1)+theta_field_size
  jarclsslt_fi(lev) = jarclsslt_fi(lev-1)+theta_field_size
  jarclsslt_jt(lev) = jarclsslt_jt(lev-1)+theta_field_size
  jarclsulp_ac(lev) = jarclsulp_ac(lev-1)+theta_field_size
  jarclsulp_ak(lev) = jarclsulp_ak(lev-1)+theta_field_size
  jarclsulp_di(lev) = jarclsulp_di(lev-1)+theta_field_size
  jarcldust_b1(lev) = jarcldust_b1(lev-1)+theta_field_size
  jarcldust_b2(lev) = jarcldust_b2(lev-1)+theta_field_size
  jarcldust_b3(lev) = jarcldust_b3(lev-1)+theta_field_size
  jarcldust_b4(lev) = jarcldust_b4(lev-1)+theta_field_size
  jarcldust_b5(lev) = jarcldust_b5(lev-1)+theta_field_size
  jarcldust_b6(lev) = jarcldust_b6(lev-1)+theta_field_size
  jarclocff_fr(lev) = jarclocff_fr(lev-1)+theta_field_size
  jarclocff_ag(lev) = jarclocff_ag(lev-1)+theta_field_size
  jarclocff_ic(lev) = jarclocff_ic(lev-1)+theta_field_size
  jarcldlta_dl(lev) = jarcldlta_dl(lev-1)+theta_field_size
END DO

! Mineral dust scheme

jsoil_clay   =si(418,Sect_No,im_index)  ! soil clay fraction
jsoil_silt   =si(419,Sect_No,im_index)  ! soil silt fraction
jsoil_sand   =si(420,Sect_No,im_index)  ! soil sand fraction

jdust_mrel1=si(421,Sect_No,im_index) !relative soil mass in div1
jdust_mrel2=si(422,Sect_No,im_index) !relative soil mass in div2
jdust_mrel3=si(423,Sect_No,im_index) !relative soil mass in div3
jdust_mrel4=si(424,Sect_No,im_index) !relative soil mass in div4
jdust_mrel5=si(425,Sect_No,im_index) !relative soil mass in div5
jdust_mrel6=si(426,Sect_No,im_index) !relative soil mass in div6


jdust_div1(tdims%k_start)=si(431,Sect_No,im_index)  ! dust mmr, division 1
jdust_div2(tdims%k_start)=si(432,Sect_No,im_index)  ! dust mmr, division 2
jdust_div3(tdims%k_start)=si(433,Sect_No,im_index)  ! dust mmr, division 3
jdust_div4(tdims%k_start)=si(434,Sect_No,im_index)  ! dust mmr, division 4
jdust_div5(tdims%k_start)=si(435,Sect_No,im_index)  ! dust mmr, division 5
jdust_div6(tdims%k_start)=si(436,Sect_No,im_index)  ! dust mmr, division 6

DO lev = tdims%k_start+1, tdims%k_end
  jdust_div1(lev)=jdust_div1(lev-1)+theta_off_size
  jdust_div2(lev)=jdust_div2(lev-1)+theta_off_size
  jdust_div3(lev)=jdust_div3(lev-1)+theta_off_size
  jdust_div4(lev)=jdust_div4(lev-1)+theta_off_size
  jdust_div5(lev)=jdust_div5(lev-1)+theta_off_size
  jdust_div6(lev)=jdust_div6(lev-1)+theta_off_size
END DO

! Ammonium nitrate scheme
jnitr_acc (tdims%k_start)= si(117,Sect_No,im_index) !Accumulation nitrate MMR
jnitr_diss(tdims%k_start)= si(118,Sect_No,im_index) !Dissolved nitrate MMR
DO lev = tdims%k_start+1, tdims%k_end
  jnitr_acc(lev) = jnitr_acc(lev-1)+theta_off_size
  jnitr_diss(lev) = jnitr_diss(lev-1)+theta_off_size
END DO

! Add for Carbon cycle
j_triffid_co2 = si(249,Sect_No,im_index)
j_co2flux = si(250,Sect_No,im_index)
j_co2_emits  = si(251,Sect_No,im_index)
jco2(tdims%k_start)      = si(252,Sect_No,im_index)
DO lev= tdims%k_start+1, tdims%k_end
  jmurk_source(lev) = jmurk_source(lev-1)+theta_field_size
  jmurk(lev) = jmurk(lev-1)+theta_off_size

  ! For Sulphur Cycle variables
  jso2(lev)=jso2(lev-1)+theta_off_size
  jdms(lev)=jdms(lev-1)+theta_off_size
  jso4_aitken(lev)=jso4_aitken(lev-1)+theta_off_size
  jso4_accu(lev)=jso4_accu(lev-1)+theta_off_size
  jso4_diss(lev)=jso4_diss(lev-1)+theta_off_size
  jh2o2(lev)=jh2o2(lev-1)+theta_off_size
  jso2_natem(lev)=jso2_natem(lev-1)+theta_field_size
  joh(lev) = joh(lev-1)+theta_field_size
  jho2(lev) = jho2(lev-1)+theta_field_size
  jnh3(lev)      = jnh3(lev-1)+theta_off_size
  jsoot_new(lev) = jsoot_new(lev-1)+theta_off_size
  jsoot_agd(lev) = jsoot_agd(lev-1)+theta_off_size
  jsoot_cld(lev) = jsoot_cld(lev-1)+theta_off_size
  jbmass_new(lev) = jbmass_new(lev-1)+theta_off_size
  jbmass_agd(lev) = jbmass_agd(lev-1)+theta_off_size
  jbmass_cld(lev) = jbmass_cld(lev-1)+theta_off_size
  jocff_new(lev) = jocff_new(lev-1)+theta_off_size
  jocff_agd(lev) = jocff_agd(lev-1)+theta_off_size
  jocff_cld(lev) = jocff_cld(lev-1)+theta_off_size
  jh2o2_limit(lev)=jh2o2_limit(lev-1)+theta_field_size
  jo3_chem(lev)=jo3_chem(lev-1)+theta_field_size

  ! For Carbon Cycle variables
  jco2(lev)=jco2(lev-1)+theta_off_size
END DO

! Add user ancillaries

juser_anc1  = si(301,Sect_No,im_index)
juser_anc2  = si(302,Sect_No,im_index)
juser_anc3  = si(303,Sect_No,im_index)
juser_anc4  = si(304,Sect_No,im_index)
juser_anc5  = si(305,Sect_No,im_index)
juser_anc6  = si(306,Sect_No,im_index)
juser_anc7  = si(307,Sect_No,im_index)
juser_anc8  = si(308,Sect_No,im_index)
juser_anc9  = si(309,Sect_No,im_index)
juser_anc10 = si(310,Sect_No,im_index)
juser_anc11 = si(311,Sect_No,im_index)
juser_anc12 = si(312,Sect_No,im_index)
juser_anc13 = si(313,Sect_No,im_index)
juser_anc14 = si(314,Sect_No,im_index)
juser_anc15 = si(315,Sect_No,im_index)
juser_anc16 = si(316,Sect_No,im_index)
juser_anc17 = si(317,Sect_No,im_index)
juser_anc18 = si(318,Sect_No,im_index)
juser_anc19 = si(319,Sect_No,im_index)
juser_anc20 = si(320,Sect_No,im_index)
juser_mult1(1)  = si(321,Sect_No,im_index)
juser_mult2(1)  = si(322,Sect_No,im_index)
juser_mult3(1)  = si(323,Sect_No,im_index)
juser_mult4(1)  = si(324,Sect_No,im_index)
juser_mult5(1)  = si(325,Sect_No,im_index)
juser_mult6(1)  = si(326,Sect_No,im_index)
juser_mult7(1)  = si(327,Sect_No,im_index)
juser_mult8(1)  = si(328,Sect_No,im_index)
juser_mult9(1)  = si(329,Sect_No,im_index)
juser_mult10(1) = si(330,Sect_No,im_index)
juser_mult11(1) = si(331,Sect_No,im_index)
juser_mult12(1) = si(332,Sect_No,im_index)
juser_mult13(1) = si(333,Sect_No,im_index)
juser_mult14(1) = si(334,Sect_No,im_index)
juser_mult15(1) = si(335,Sect_No,im_index)
juser_mult16(1) = si(336,Sect_No,im_index)
juser_mult17(1) = si(337,Sect_No,im_index)
juser_mult18(1) = si(338,Sect_No,im_index)
juser_mult19(1) = si(339,Sect_No,im_index)
juser_mult20(1) = si(340,Sect_No,im_index)

! Set for multi-level user ancillaries
DO lev=2,model_levels
  juser_mult1(lev)  = juser_mult1(lev-1)+theta_field_size
  juser_mult2(lev)  = juser_mult2(lev-1)+theta_field_size
  juser_mult3(lev)  = juser_mult3(lev-1)+theta_field_size
  juser_mult4(lev)  = juser_mult4(lev-1)+theta_field_size
  juser_mult5(lev)  = juser_mult5(lev-1)+theta_field_size
  juser_mult6(lev)  = juser_mult6(lev-1)+theta_field_size
  juser_mult7(lev)  = juser_mult7(lev-1)+theta_field_size
  juser_mult8(lev)  = juser_mult8(lev-1)+theta_field_size
  juser_mult9(lev)  = juser_mult9(lev-1)+theta_field_size
  juser_mult10(lev) = juser_mult10(lev-1)+theta_field_size
  juser_mult11(lev) = juser_mult11(lev-1)+theta_field_size
  juser_mult12(lev) = juser_mult12(lev-1)+theta_field_size
  juser_mult13(lev) = juser_mult13(lev-1)+theta_field_size
  juser_mult14(lev) = juser_mult14(lev-1)+theta_field_size
  juser_mult15(lev) = juser_mult15(lev-1)+theta_field_size
  juser_mult16(lev) = juser_mult16(lev-1)+theta_field_size
  juser_mult17(lev) = juser_mult17(lev-1)+theta_field_size
  juser_mult18(lev) = juser_mult18(lev-1)+theta_field_size
  juser_mult19(lev) = juser_mult19(lev-1)+theta_field_size
  juser_mult20(lev) = juser_mult20(lev-1)+theta_field_size
END DO

! Tiled vegetation and triffid
jfrac_typ     = si(216,Sect_No,im_index) ! surface type fractions
jfrac_con1    = si(442,Sect_No,im_index) ! surface type fractions
jfrac_con2    = si(443,Sect_No,im_index) ! surface type fractions
jfrac_con3    = si(444,Sect_No,im_index) ! surface type fractions
jfrac_con4    = si(445,Sect_No,im_index) ! surface type fractions
jfrac_con5    = si(446,Sect_No,im_index) ! surface type fractions
jfrac_con6    = si(447,Sect_No,im_index) ! surface type fractions
jfrac_con7    = si(448,Sect_No,im_index) ! surface type fractions
jfrac_con8    = si(449,Sect_No,im_index) ! surface type fractions
jfrac_con9    = si(450,Sect_No,im_index) ! surface type fractions
jlai_pft      = si(217,Sect_No,im_index) ! leaf area index of PFTs
jcanht_pft    = si(218,Sect_No,im_index) ! canopy height of PFTs
jdisturb      = si(219,Sect_No,im_index) ! Veg disturbed fraction
jdisturb_prev = si(286,Sect_No,im_index) ! Previous Veg disturbed fraction
jagr_crop     = si(448,Sect_No,im_index) ! Crop disturbed fraction
jagr_crop_prev= si(449,Sect_No,im_index) ! Previous crop disturbed fraction
jpasture      = si(458,Sect_No,im_index) ! Pasture fraction
jpasture_prev = si(459,Sect_No,im_index) ! Previous pasture fraction
jwoodprod_fast = si(287,Sect_No,im_index) ! Wood product pool (fast)      
jwoodprod_med  = si(288,Sect_No,im_index) ! Wood product pool (medium)    
jwoodprod_slow = si(289,Sect_No,im_index) ! Wood product pool (slow)      
jsoil_alb     = si(220,Sect_No,im_index) ! Snow-free soil albedo
jobs_alb_sw   = si(243,Sect_No,im_index) ! Observed snow-free SW albedo
jobs_alb_vis  = si(244,Sect_No,im_index) ! Observed snow-free VIS albedo
jobs_alb_nir  = si(245,Sect_No,im_index) ! Observed snow-free NIR albedo
jsoil_carb    = si(223,Sect_No,im_index) ! Soil carbon content
jsoil_carb1   = si(466,Sect_No,im_index) ! Soil carbon content DPM
jsoil_carb2   = si(467,Sect_No,im_index) ! Soil carbon content RPM
jsoil_carb3   = si(468,Sect_No,im_index) ! Soil carbon content BIO
jsoil_carb4   = si(469,Sect_No,im_index) ! Soil carbon content HUM
jsoil_nitro1   = si(442,Sect_No,im_index) ! Soil nitrogen content DPM
jsoil_nitro2   = si(443,Sect_No,im_index) ! Soil nitrogen content RPM
jsoil_nitro3   = si(444,Sect_No,im_index) ! Soil nitrogen content BIO
jsoil_nitro4   = si(445,Sect_No,im_index) ! Soil nitrogen content HUM
jsoil_inorgnit  = si(446,Sect_No,im_index) ! Soil inorganic nitrogen content 
j_n_deposition  = si(447,Sect_No,im_index) ! Nitrogen deposition on land
jnpp_pft_acc  = si(224,Sect_No,im_index) ! Accumulated NPP on PFTs
jg_lf_pft_acc = si(225,Sect_No,im_index) ! Accumulated leaf
!                                              ! turnover rate on PFTs
jg_phlf_pft_acc=si(226,Sect_No,im_index) ! Accumulat. phenological
!                                              ! leaf turnover rate PFTs
jrsp_w_pft_acc= si(227,Sect_No,im_index) ! Accum. wood resp PFTs
jrsp_s_acc    = si(228,Sect_No,im_index) ! Accumulated soil resp
jrsp_s_acc1 = si(470,Sect_No,im_index)  ! Soil respiration DPM
jrsp_s_acc2 = si(471,Sect_No,im_index)  ! Soil respiration RPM
jrsp_s_acc3 = si(472,Sect_No,im_index)  ! Soil respiration BIO
jrsp_s_acc4 = si(473,Sect_No,im_index)  ! Soil respiration HUM
jcan_water_tile=si(229,Sect_No,im_index) ! Canopy water content
!                                              ! on tiles
jcatch_tile   = si(230,Sect_No,im_index) ! Canopy capacity on
!                                              ! tiles
jrgrain_tile  = si(231,Sect_No,im_index) ! Snow grain size on
!                                              ! tiles
jtstar_tile   = si(233,Sect_No,im_index) ! Tiled surface temp
jtsurf_elev_surft   = si(576,Sect_No,im_index) ! Temperature of elevated 
                                         ! subsurface tiles
jz0_tile      = si(234,Sect_No,im_index) ! Tiled surface roughness
jz0h_tile     = si(246,Sect_No,im_index) ! Tiled surface thermal
!                                              ! roughness
! Stash number for snow tile not finalised yet
jsnodep_tile  = si(240,Sect_No,im_index) ! Tiled snow depth
jinfil_tile   = si(236,Sect_No,im_index) ! Max tile infilt rate
! Stash codes for DORL, LW_DOWN, SW_TILE not finalised yet
jdolr         = si(239,Sect_No,im_index) ! TOA surface up LW
jlw_down      = si(238,Sect_No,im_index) ! Surface down LW
jsw_tile      = si(237,Sect_No,im_index) ! Surface net SW on tiles

! MORUSES urban scheme
jurbhgt             = si(494,Sect_No,im_index) ! Building height
jurbhwr             = si(495,Sect_No,im_index) ! Height to width
jurbwrr             = si(496,Sect_No,im_index) ! Width ratio
jurbdisp            = si(497,Sect_No,im_index) ! Displacement height
jurbztm             = si(498,Sect_No,im_index) ! Effective roughness
                                               ! length for momentum
jurbalbwl           = si(499,Sect_No,im_index) ! Wall albedo
jurbalbrd           = si(500,Sect_No,im_index) ! Road albedo
jurbemisw           = si(501,Sect_No,im_index) ! Wall emmissivity
jurbemisr           = si(502,Sect_No,im_index) ! Road emmissivity

! River routing fields
jriv_sequence  = si(151,Sect_No,im_index) ! River sequence
jriv_direction = si(152,Sect_No,im_index) ! River Direction
jriv_storage   = si(153,Sect_No,im_index) ! River Water Storage
jtot_surfroff  = si(155,Sect_No,im_index) ! Acc. surface runoff
jtot_subroff   = si(156,Sect_No,im_index) ! Acc. sub-surf runoff
! Set pointer for inland basin outflow
jriv_inlandatm    = si(511,Sect_No,im_index)
 !Inland basin outflow

!PV-tracer
jdPV_rad=si(577,sect_No,im_index)
jdPV_sw=si(578,sect_No,im_index)
jdPV_lw=si(579,sect_No,im_index)
jdPV_mic=si(580,sect_No,im_index)
jdPV_gwd=si(581,sect_No,im_index)
jdPV_ph1=si(582,sect_No,im_index)
jdPV_conv=si(583,sect_No,im_index)
jdPV_bl=si(584,sect_No,im_index)
jdPV_stph=si(585,sect_No,im_index)
jdPV_cld=si(586,sect_No,im_index)
jdPV_iau=si(587,sect_No,im_index)
jdPV_nud=si(588,sect_No,im_index)
jdPV_tot=si(589,sect_No,im_index)
jdPV_adv=si(590,sect_No,im_index)
jdPV_sol=si(591,sect_No,im_index)
jdPV_mass=si(592,sect_No,im_index)
jadv_only_PV=si(593,sect_No,im_index)

DO lev= tdims%k_start+1, tdims%k_end
  jdPV_rad(lev) = jdPV_rad(lev-1) + theta_field_size
  jdPV_sw(lev)  = jdPV_sw(lev-1)  + theta_field_size
  jdPV_lw(lev)  = jdPV_lw(lev-1)  + theta_field_size    
  jdPV_mic(lev) = jdPV_mic(lev-1) + theta_field_size
  jdPV_gwd(lev) = jdPV_gwd(lev-1) + theta_field_size
  jdPV_ph1(lev) = jdPV_ph1(lev-1) + theta_field_size
  jdPV_conv(lev)= jdPV_conv(lev-1) + theta_field_size
  jdPV_bl(lev)  = jdPV_bl(lev-1)   + theta_field_size
  jdPV_stph(lev)= jdPV_stph(lev-1) + theta_field_size
  jdPV_cld(lev) = jdPV_cld(lev-1) + theta_field_size
  jdPV_iau(lev) = jdPV_iau(lev-1) + theta_field_size
  jdPV_nud(lev) = jdPV_nud(lev-1) + theta_field_size
  jdPV_tot(lev) = jdPV_tot(lev-1) + theta_field_size
  jdPV_adv(lev) = jdPV_adv(lev-1) + theta_field_size
  jdPV_sol(lev) = jdPV_sol(lev-1) + theta_field_size
  jdPV_mass(lev) = jdPV_mass(lev-1) + theta_field_size
  jadv_only_PV(lev) = jadv_only_PV(lev-1) + theta_field_size
END DO

jdtheta_0       =si(600,sect_No,im_index)
jdtheta_bl      =si(601,sect_No,im_index)
jdtheta_bl_mix  =si(602,sect_No,im_index)
jdtheta_bl_LH   =si(603,sect_No,im_index)
jdtheta_conv    =si(604,sect_No,im_index)
jdtheta_mic     =si(605,sect_No,im_index)
jdtheta_rad     =si(606,sect_No,im_index)
jdtheta_SW      =si(607,sect_No,im_index)
jdtheta_LW      =si(608,sect_No,im_index)
jdtheta_slow    =si(609,sect_No,im_index)
jdtheta_cld     =si(610,sect_No,im_index)

DO lev= tdims%k_start+1, tdims%k_end
  jdtheta_0(lev)    = jdtheta_0(lev-1)    + theta_field_size
  jdtheta_bl(lev)   = jdtheta_bl(lev-1)   + theta_field_size
  jdtheta_bl_mix(lev) = jdtheta_bl_mix(lev-1)    + theta_field_size
  jdtheta_bl_LH(lev)  = jdtheta_bl_LH(lev-1)     + theta_field_size
  jdtheta_conv(lev) = jdtheta_conv(lev-1) + theta_field_size
  jdtheta_mic(lev)  = jdtheta_mic(lev-1)  + theta_field_size
  jdtheta_rad(lev)  = jdtheta_rad(lev-1)  + theta_field_size
  jdtheta_SW(lev)   = jdtheta_SW(lev-1)   + theta_field_size
  jdtheta_LW(lev)   = jdtheta_LW(lev-1)   + theta_field_size
  jdtheta_slow(lev) = jdtheta_slow(lev-1)  + theta_field_size
  jdtheta_cld(lev)  = jdtheta_cld(lev-1)  + theta_field_size    
END DO   

! required for energy correction
jnet_flux=si(222,Sect_No,im_index)    ! store for energy flux
jnet_mflux=si(235,Sect_No,im_index)   ! store for moisture flux

!Fields carried forward from previous version.
jtstar_anom    = si(39,Sect_No,im_index)

! JULES version 2 prognostics
jsnowdepth      = si(376,Sect_No,im_index) ! Snow depth on ground on tiles (m)
jrho_snow_grnd  = si(377,Sect_No,im_index) ! Snowpack bulk density (kg/m3)
jnsnow          = si(380,Sect_No,im_index) ! Number of snow layers on ground on tiles
jds             = si(381,Sect_No,im_index) ! Snow layer thickness (m)
jsice           = si(382,Sect_No,im_index) ! Snow layer ice mass on tiles (Kg/m2)
jsliq           = si(383,Sect_No,im_index) ! Snow layer liquid mass on tiles (Kg/m2)
jtsnowlayer     = si(384,Sect_No,im_index) ! Snow layer temperature (K)
jrho_snow       = si(385,Sect_No,im_index) ! Snow layer densities (kg/m3)
jrgrainl        = si(386,Sect_No,im_index) ! Snow layer grain size on tiles (microns)

! FLake lake scheme prognostics
jlake_depth  = si(291,Sect_No,im_index) ! lake depth (m)
jlake_fetch  = si(292,Sect_No,im_index) ! typical wind fetch (m)
jlake_t_mean = si(293,Sect_No,im_index) ! lake mean temperature (K)
jlake_t_mxl  = si(294,Sect_No,im_index) ! lake mixed-layer temperature (K)
jlake_t_ice  = si(295,Sect_No,im_index) ! temperature at upper boundary of lake ice (K)
jlake_h_mxl  = si(296,Sect_No,im_index) ! lake mixed-layer depth (m)
jlake_h_ice  = si(297,Sect_No,im_index) ! lake ice thickness (m)
jlake_shape  = si(298,Sect_No,im_index) ! thermocline shape factor
jlake_g_dt   = si(299,Sect_No,im_index) ! lake ht.flx / dT (W m-2 K-1)

! Section 54 GLOMAP_CLIM aerosol climatology
sect_no = stashcode_glomap_clim_sec     ! GLOMAP_CLIM section
jgc_nd_nuc_sol (tdims%k_start) = si(101,Sect_No,im_index)
jgc_nuc_sol_su (tdims%k_start) = si(102,Sect_No,im_index)
jgc_nuc_sol_oc (tdims%k_start) = si(126,Sect_No,im_index)
jgc_nd_ait_sol (tdims%k_start) = si(103,Sect_No,im_index)
jgc_ait_sol_su (tdims%k_start) = si(104,Sect_No,im_index)
jgc_ait_sol_bc (tdims%k_start) = si(105,Sect_No,im_index)
jgc_ait_sol_oc (tdims%k_start) = si(106,Sect_No,im_index)
jgc_nd_acc_sol (tdims%k_start) = si(107,Sect_No,im_index)
jgc_acc_sol_su (tdims%k_start) = si(108,Sect_No,im_index)
jgc_acc_sol_bc (tdims%k_start) = si(109,Sect_No,im_index)
jgc_acc_sol_oc (tdims%k_start) = si(110,Sect_No,im_index)
jgc_acc_sol_ss (tdims%k_start) = si(111,Sect_No,im_index)
jgc_nd_cor_sol (tdims%k_start) = si(113,Sect_No,im_index)
jgc_cor_sol_su (tdims%k_start) = si(114,Sect_No,im_index)
jgc_cor_sol_bc (tdims%k_start) = si(115,Sect_No,im_index)
jgc_cor_sol_oc (tdims%k_start) = si(116,Sect_No,im_index)
jgc_cor_sol_ss (tdims%k_start) = si(117,Sect_No,im_index)
jgc_nd_ait_ins (tdims%k_start) = si(119,Sect_No,im_index)
jgc_ait_ins_bc (tdims%k_start) = si(120,Sect_No,im_index)
jgc_ait_ins_oc (tdims%k_start) = si(121,Sect_No,im_index)
! Set for multi-level
DO lev=tdims%k_start+1, tdims%k_end
  ! Note theta_field_size is used for zero halos
  jgc_nd_nuc_sol(lev)  = jgc_nd_nuc_sol(lev-1)+theta_field_size
  jgc_nuc_sol_su(lev)  = jgc_nuc_sol_su(lev-1)+theta_field_size
  jgc_nuc_sol_oc(lev)  = jgc_nuc_sol_oc(lev-1)+theta_field_size
  jgc_nd_ait_sol(lev)  = jgc_nd_ait_sol(lev-1)+theta_field_size
  jgc_ait_sol_su(lev)  = jgc_ait_sol_su(lev-1)+theta_field_size
  jgc_ait_sol_bc(lev)  = jgc_ait_sol_bc(lev-1)+theta_field_size
  jgc_ait_sol_oc(lev)  = jgc_ait_sol_oc(lev-1)+theta_field_size
  jgc_nd_acc_sol(lev)  = jgc_nd_acc_sol(lev-1)+theta_field_size
  jgc_acc_sol_su(lev)  = jgc_acc_sol_su(lev-1)+theta_field_size
  jgc_acc_sol_bc(lev)  = jgc_acc_sol_bc(lev-1)+theta_field_size
  jgc_acc_sol_oc(lev)  = jgc_acc_sol_oc(lev-1)+theta_field_size
  jgc_acc_sol_ss(lev)  = jgc_acc_sol_ss(lev-1)+theta_field_size
  jgc_nd_cor_sol(lev)  = jgc_nd_cor_sol(lev-1)+theta_field_size
  jgc_cor_sol_su(lev)  = jgc_cor_sol_su(lev-1)+theta_field_size
  jgc_cor_sol_bc(lev)  = jgc_cor_sol_bc(lev-1)+theta_field_size
  jgc_cor_sol_oc(lev)  = jgc_cor_sol_oc(lev-1)+theta_field_size
  jgc_cor_sol_ss(lev)  = jgc_cor_sol_ss(lev-1)+theta_field_size
  jgc_nd_ait_ins(lev)  = jgc_nd_ait_ins(lev-1)+theta_field_size
  jgc_ait_ins_bc(lev)  = jgc_ait_ins_bc(lev-1)+theta_field_size
  jgc_ait_ins_oc(lev)  = jgc_ait_ins_oc(lev-1)+theta_field_size
END DO

! Set all pointers referencing Lateral Boundary Conditions

Sect_No=31  ! LBC section

jorog_lbc     = si(1,Sect_No,im_index)
ju_lbc        = si(2,Sect_No,im_index)
jv_lbc        = si(3,Sect_No,im_index)
jw_lbc        = si(4,Sect_No,im_index)
jrho_lbc      = si(5,Sect_No,im_index)
jtheta_lbc    = si(6,Sect_No,im_index)
jq_lbc        = si(7,Sect_No,im_index)
jqcl_lbc      = si(8,Sect_No,im_index)
jqcf_lbc      = si(9,Sect_No,im_index)
jexner_lbc    = si(10,Sect_No,im_index)
ju_adv_lbc    = si(11,Sect_No,im_index)
jv_adv_lbc    = si(12,Sect_No,im_index)
jw_adv_lbc    = si(13,Sect_No,im_index)
jqcf2_lbc     = si(14,Sect_No,im_index)
jqrain_lbc    = si(15,Sect_No,im_index)
jqgraup_lbc   = si(16,Sect_No,im_index)
jcf_bulk_lbc  = si(17,Sect_No,im_index)
jcf_liquid_lbc= si(18,Sect_No,im_index)
jcf_frozen_lbc= si(19,Sect_No,im_index)
jmurk_lbc     = si(20,Sect_No,im_index)

jdust_div1_lbc = si(23,Sect_No,im_index)
jdust_div2_lbc = si(24,Sect_No,im_index)
jdust_div3_lbc = si(25,Sect_No,im_index)
jdust_div4_lbc = si(26,Sect_No,im_index)
jdust_div5_lbc = si(27,Sect_No,im_index)
jdust_div6_lbc = si(28,Sect_No,im_index)
jso2_lbc       = si(29,Sect_No,im_index)
jdms_lbc       = si(30,Sect_No,im_index)
jso4_aitken_lbc= si(31,Sect_No,im_index)
jso4_accu_lbc  = si(32,Sect_No,im_index)
jso4_diss_lbc  = si(33,Sect_No,im_index)
jnh3_lbc       = si(35,Sect_No,im_index)
jsoot_new_lbc  = si(36,Sect_No,im_index)
jsoot_agd_lbc  = si(37,Sect_No,im_index)
jsoot_cld_lbc  = si(38,Sect_No,im_index)
jbmass_new_lbc = si(39,Sect_No,im_index)
jbmass_agd_lbc = si(40,Sect_No,im_index)
jbmass_cld_lbc = si(41,Sect_No,im_index)
jocff_new_lbc  = si(42,Sect_No,im_index)
jocff_agd_lbc  = si(43,Sect_No,im_index)
jocff_cld_lbc  = si(44,Sect_No,im_index)
jnitr_acc_lbc  = si(45,Sect_No,im_index)
jnitr_diss_lbc = si(46,Sect_No,im_index)

ju_lbc_tend     = si(257,Sect_No,im_index)
jv_lbc_tend     = si(258,Sect_No,im_index)
jw_lbc_tend     = si(259,Sect_No,im_index)
jrho_lbc_tend   = si(260,Sect_No,im_index)
jtheta_lbc_tend = si(261,Sect_No,im_index)
jq_lbc_tend     = si(262,Sect_No,im_index)
jqcl_lbc_tend   = si(263,Sect_No,im_index)
jqcf_lbc_tend   = si(264,Sect_No,im_index)
jexner_lbc_tend = si(265,Sect_No,im_index)
ju_adv_lbc_tend = si(266,Sect_No,im_index)
jv_adv_lbc_tend = si(267,Sect_No,im_index)
jw_adv_lbc_tend = si(268,Sect_No,im_index)
jqcf2_lbc_tend   = si(269,Sect_No,im_index)
jqrain_lbc_tend  = si(270,Sect_No,im_index)
jqgraup_lbc_tend = si(271,Sect_No,im_index)
jcf_bulk_lbc_tend  = si(272,Sect_No,im_index)
jcf_liquid_lbc_tend= si(273,Sect_No,im_index)
jcf_frozen_lbc_tend= si(274,Sect_No,im_index)
jmurk_lbc_tend   = si(275,Sect_No,im_index)

jdust_div1_lbc_tend = si(276,Sect_No,im_index)
jdust_div2_lbc_tend = si(277,Sect_No,im_index)
jdust_div3_lbc_tend = si(278,Sect_No,im_index)
jdust_div4_lbc_tend = si(279,Sect_No,im_index)
jdust_div5_lbc_tend = si(280,Sect_No,im_index)
jdust_div6_lbc_tend = si(281,Sect_No,im_index)
jso2_lbc_tend       = si(282,Sect_No,im_index)
jdms_lbc_tend       = si(283,Sect_No,im_index)
jso4_aitken_lbc_tend= si(284,Sect_No,im_index)
jso4_accu_lbc_tend  = si(285,Sect_No,im_index)
jso4_diss_lbc_tend  = si(286,Sect_No,im_index)
jnh3_lbc_tend       = si(288,Sect_No,im_index)
jsoot_new_lbc_tend  = si(289,Sect_No,im_index)
jsoot_agd_lbc_tend  = si(290,Sect_No,im_index)
jsoot_cld_lbc_tend  = si(291,Sect_No,im_index)
jbmass_new_lbc_tend = si(292,Sect_No,im_index)
jbmass_agd_lbc_tend = si(293,Sect_No,im_index)
jbmass_cld_lbc_tend = si(294,Sect_No,im_index)
jocff_new_lbc_tend  = si(295,Sect_No,im_index)
jocff_agd_lbc_tend  = si(296,Sect_No,im_index)
jocff_cld_lbc_tend  = si(297,Sect_No,im_index)
jnitr_acc_lbc_tend  = si(298,Sect_No,im_index)
jnitr_diss_lbc_tend = si(299,Sect_No,im_index)

! Set pointers for Tracer prognostics
! Tracer prognostics are now in section 33, not section 0
sect_no = 33   ! tracers section
jvar=0         ! JVAR+1 is the current tracer to be found
IF (tr_vars >  0) THEN
  DO ivar=a_tracer_first,a_tracer_last
    IF (si(ivar,sect_no,im_index) /= 1) THEN ! tracer in use
      jvar=jvar+1
      jtracer(tdims_s%k_start,jvar) = si(ivar,sect_no,im_index)
      ! Set up array containing stash item codes for active
      ! tracer prognostics (1:TR_VARS)
      A_TR_StashItem(jvar) = ivar
      DO lev = tdims_s%k_start+1, tdims_s%k_end
        jtracer(lev,jvar)=jtracer(lev-1,jvar)+theta_off_size
      END DO
      a_tr_index(ivar-a_tracer_first+1)=jvar
    ELSE
      ! If tracer not active, set value to -1
      a_tr_index(ivar-a_tracer_first+1) = -1
    END IF
  END DO
ELSE
  ! Ensure a sensible address even if no tracers
  jtracer(tdims_s%k_start,1)=1
END IF
IF (jvar /= tr_vars) THEN
  WRITE(umMessage,*) 'STATMPT: TR_VARS and SI are inconsistent'
  CALL umPrint(umMessage,src='set_atm_pointers')
  WRITE(umMessage,*) 'TR_VARS=',tr_vars,' .     But, SI implies :',jvar
  CALL umPrint(umMessage,src='set_atm_pointers')
  cmessage=  'STATMPT: TR_VARS and SI  inconsistent, see output'
  icode=100
  GO TO 9999 ! error return
END IF

! UKCA tracer prognostics are in section 34.
sect_no = ukca_sect   ! UKCA tracers section
jvar=0         ! JVAR+1 is the current tracer to be found
IF (tr_ukca >  0) THEN
  DO ivar=a_ukca_first,a_ukca_last
    IF (si(ivar,sect_no,im_index) /= 1) THEN ! tracer in use
      jvar=jvar+1
      jtr_ukca(tdims_s%k_start,jvar) = si(ivar,sect_no,im_index)
      ! Set up array containing stash item codes for active
      ! tracer prognostics (1:TR_UKCA)
      UKCA_TR_StashItem(jvar) = ivar
      DO lev = tdims_s%k_start+1, tdims_s%k_end
        jtr_ukca(lev,jvar)=jtr_ukca(lev-1,jvar)+theta_off_size
      END DO
      a_ukca_index(ivar-a_ukca_first+1)=jvar
    ELSE
      ! If tracer not active, set value to -1
      a_ukca_index(ivar-a_ukca_first+1) = -1
    END IF
  END DO
ELSE
  ! Ensure a sensible address when no tracers
  jtr_ukca(tdims_s%k_start,1)=1
END IF

IF (jvar /= tr_ukca) THEN
  WRITE(umMessage,*) 'STATMPT: TR_UKCA and SI are inconsistent'
  CALL umPrint(umMessage,src='set_atm_pointers')
  WRITE(umMessage,*) 'TR_UKCA=',tr_ukca,' .     But, SI implies :',jvar
  CALL umPrint(umMessage,src='set_atm_pointers')
  cmessage=  'STATMPT: TR_UKCA and SI  inconsistent, see output'
  icode=100
  GO TO 9999 ! error return
END IF

IF (l_sulpc_online_oxidants .AND. l_ukca) THEN
  jo3_ukca(tdims%k_start)   =si(ukca_item_sulpc(1),sect_no,im_index)
  jhno3_ukca(tdims%k_start) =si(ukca_item_sulpc(2),sect_no,im_index)
  jh2o2_ukca(tdims%k_start) =si(ukca_item_sulpc(3),sect_no,im_index)
  joh_ukca(tdims%k_start)   =si(ukca_item_sulpc(4),sect_no,im_index)
  jho2_ukca(tdims%k_start)  =si(ukca_item_sulpc(5),sect_no,im_index)
END IF

! Set pointers for Tracer lateral boundary data

      ! Free tracer lbc data is in section 36
      ! The tracer lbc stash item codes in A_TR_LBC_StashItem
      ! are set up in INBOUNDA

sect_no = 36   ! tracer lbcs section
jvar = 0
IF (tr_lbc_vars > 0) THEN
  DO ivar= a_tracer_first,a_tracer_last
    IF (si(ivar,sect_no,im_index) /= 1) THEN
      jvar = jvar + 1
      jtracer_lbc(jvar)        = si(ivar,sect_no,im_index)
      jtracer_lbc_tend(jvar)   = si(256+ivar,sect_no,im_index)
      A_TR_LBC_StashItem(jvar) = ivar
    END IF
  END DO
ELSE
  ! Ensure a sensible address even if no tracer lbcs
  jtracer_lbc(1)=1
  jtracer_lbc_tend(1)=1
END IF ! IF (TR_LBC_VARS > 0)

! Set up array (1:TR_VARS) pointing to location in TRACER_LBC
! array for each prognostic tracer (A_tr_active_lbc_index)
! This allows tracer prognostics to be active even if
! there are no lbcs, and lbc data to be ignored if the
! tracer prognostic is not active. If there is no lbc data
! for a particular tracer, the array is set to -1.

DO ivar= 1,tr_vars
  A_tr_active_lbc_index(ivar) = -1
  DO jvar= 1,tr_lbc_vars
    IF (A_TR_LBC_StashItem(jvar) == A_TR_StashItem(ivar)) THEN
      A_tr_active_lbc_index(ivar) = jvar
      WRITE(umMessage,*) 'A_tr_active_lbc_index:',ivar,jvar,             &
            A_TR_LBC_StashItem(jvar),A_TR_StashItem(ivar),        &
            A_tr_active_lbc_index(ivar)
      CALL umPrint(umMessage,src='set_atm_pointers')
    END IF
  END DO
END DO

! UKCA tracer lbc data is in section 37
! The tracer lbc stash item codes in UKCA_TR_LBC_StashItem
! are set up in item_bounda_mod

sect_no = 37   ! UKCA tracer lbcs section
jvar = 0
IF (tr_lbc_ukca > 0) THEN
  DO ivar= a_ukca_first,a_ukca_last
    IF (si(ivar,sect_no,im_index) /= 1) THEN
      jvar = jvar + 1
      jtr_ukca_lbc(jvar)        = si(ivar,sect_no,im_index)
      jtr_ukca_lbc_tend(jvar)   = si(256+ivar,sect_no,im_index)
      UKCA_TR_LBC_StashItem(jvar) = ivar
    END IF
  END DO
ELSE
  ! Ensure a sensible address even if no tracer lbcs
  jtr_ukca_lbc(1)=1
  jtr_ukca_lbc_tend(1)=1
END IF ! IF (TR_LBC_UKCA > 0)

! Set up array (1:TR_UKCA) pointing to location in TRACER_LBC_UKCA
! array for each prognostic tracer (UKCA_tr_active_lbc_index)
! This allows tracer prognostics to be active even if
! there are no lbcs, and lbc data to be ignored if the
! tracer prognostic is not active. If there is no lbc data
! for a particular tracer, the array is set to -1.

DO ivar= 1,tr_ukca
  UKCA_tr_active_lbc_index(ivar) = -1
  DO jvar= 1,tr_lbc_ukca
    IF (UKCA_TR_LBC_StashItem(jvar) == UKCA_TR_StashItem(ivar))   &
    THEN
      UKCA_tr_active_lbc_index(ivar) = jvar
      WRITE(umMessage,*) 'UKCA_tr_active_lbc_index:',ivar,jvar,          &
            UKCA_TR_LBC_StashItem(jvar),UKCA_TR_StashItem(ivar),  &
            UKCA_tr_active_lbc_index(ivar)
      CALL umPrint(umMessage,src='set_atm_pointers')
    END IF
  END DO
END DO

! Set pointers to level dependent constants for atmosphere.
jetatheta      =1              ! eta_theta_levels(0:model_levels)
jetarho        =jetatheta+model_levels+1 ! eta_rho_levels
jrhcrit        =jetarho+model_levels+1   ! rhcrit
jsoil_thickness=jrhcrit+model_levels+1   ! soil level depths
! For height definition z(i,j,k) = zsea(k) + C(k)*zorog(i,j)
Jzseak_theta   =jsoil_thickness+model_levels+1 ! zsea (theta levs)
JCk_theta      =Jzseak_theta   +model_levels+1 ! C    (theta levs)
Jzseak_rho     =JCk_theta      +model_levels+1 ! zsea (rho levs)
JCk_rho        =Jzseak_rho     +model_levels+1 ! C    (rho levs)

IF (a_len2_rowdepc > 0 .AND. a_len2_coldepc > 0) THEN
  ! Set pointers to Row dependent constants for atmosphere.
  jphi_input_p       =1
  jphi_input_v       =jphi_input_p + a_len1_rowdepc
  ! Set pointers to Col dependent constants for atmosphere.
  jlambda_input_p    =1
  jlambda_input_u    =jlambda_input_p + a_len1_coldepc
ELSE
  jphi_input_p       =1
  jphi_input_v       =1
  jlambda_input_p    =1
  jlambda_input_u    =1
END IF

! The following if block should only be required for test purposes
! during the transition of vn5.2. This ensures old values for
! setting blev/bulev/brlev on old style dumps (before 'smooth'
! algorithm introduced and new definitions introduced).
IF (a_len2_levdepc  <=  4) THEN ! ie before 5.2 change
  JCk_rho=jetarho
  Jck_theta=jetatheta
  Jzseak_rho=jetarho
  Jzseak_theta=jetatheta
  WRITE(umMessage,*) 'SETCONA_CTL: WARNING. This dump has not been '//  &
  'reconfigured with new 5.2 set of level dependent constants:'
  CALL umPrint(umMessage,src='set_atm_pointers')
  WRITE(umMessage,*) 'PP headers will revert to 5.0/5.1 style '//       &
  'blev/bulev/brlev definitions'
  CALL umPrint(umMessage,src='set_atm_pointers')
END IF

9999 CONTINUE ! ERROR GOTO point.
 
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE set_atm_pointers
