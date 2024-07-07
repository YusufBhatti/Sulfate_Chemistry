! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Points atmosphere fields to the appropriate sections of D1
!
! Subroutine Interface:
SUBROUTINE Set_Atm_Fields ( &
      d1)

USE UM_ParParams,  ONLY: halo_type_extended, halo_type_single
USE Field_Types,   ONLY: fld_type_p, fld_type_u, fld_type_v

USE atm_fields_mod   ! atmosphere fields
USE field_length_mod ! field_length function
USE atm_fields_bounds_mod, ONLY : o3dims2, pdims, pdims_l, pdims_s,        &
                                  rdims2, tdims, tdims_l, tdims_s,         &
                                  udims, udims_l, udims_s, vdims, vdims_l, &
                                  vdims_s, wdims, wdims_l, wdims_s
USE jules_snow_mod, ONLY: nsmax
USE jules_surface_types_mod, ONLY: npft,ntype
USE atm_step_local, ONLY: land_pts_trif, npft_trif, dim_cs1

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE rimtypes, ONLY: rima_type_norm
USE lbc_mod, ONLY: lenrima
USE run_aerosol_mod, ONLY: l_sulpc_online_oxidants
USE rad_input_mod, ONLY: lexpand_ozone, lexpand_tpps_ozone
USE ukca_option_mod, ONLY: l_ukca

USE submodel_mod, ONLY: atmos_im
USE stash_array_mod, ONLY: si
USE free_tracers_inputs_mod, ONLY: a_tracer_first, a_tracer_last
USE jules_sea_seaice_mod, ONLY: nice, nice_use, l_sice_meltponds_cice, &
                                l_saldep_freeze
USE jules_vegetation_mod, ONLY: l_triffid

USE nlsizes_namelist_mod, ONLY:                                        &
    len_tot, model_levels, n_cca_lev, ntiles, rows, row_length, n_rows,&
    sm_levels, st_levels, tpps_ozone_levels, tr_lbc_ukca, tr_lbc_vars, &
    land_field

USE casim_switches, ONLY : n_casim_progs
USE casim_prognostics, ONLY :                                          &
! prognostics for casim cloud and ice microphysics
     cloudnumber, rainnumber, rain3mom, icenumber, snownumber,         &
     Snow3mom, graupnumber, graup3mom,                                 &
! prognostics for activated aerosol
     activesolliquid, activesolrain, activeinsolice, activesolice,     &
     activeinsolliquid, activesolnumber, activeinsolnumber

USE ukca_tracer_stash, ONLY: a_ukca_last, a_ukca_first

! Switches for stochastic physics
USE stochastic_physics_run_mod,  ONLY: i_pert_theta, i_pert_theta_type, &
                                       pert_theta_correl_seq
USE bl_option_mod, ONLY: off

USE atm_d1_indices_mod, ONLY: jsst, jexnersurf, jdryrho, jetadot, jthetav,     &
    jpsiws, jpsiwl, jmv, jmcl, jmcf, jmcf2, jmrain, jmgraup, ju, jv, jw, jrho, &
    jtheta, jq, jqcl, jqcf, jqcf2, jqrain, jqgraup, jcloudnumber, jrainnumber, &
    jrain3mom, jicenumber, jsnownumber, jsnow3mom, jgraupnumber, jgraup3mom,   &
    jactivesolliquid, jactivesolrain, jactiveinsolice, jactivesolice,          &
    jactiveinsolliquid, jactivesolnumber, jactiveinsolnumber, je_trb, jtsq_trb,&
    jqsq_trb, jcov_trb, jzhpar_shcu, jdPV_rad, jdPV_sw, jdPV_lw, jdPV_mic,     &
    jdPV_gwd, jdPV_ph1, jdPV_conv, jdPV_bl, jdPV_stph, jdPV_cld, jdPV_iau,     &
    jdPV_nud, jdPV_tot, jdPV_adv, jdPV_sol, jdPV_mass, jadv_only_PV,           &
    jdtheta_0, jdtheta_bl, jdtheta_bl_mix, jdtheta_bl_LH, jdtheta_conv,        &
    jdtheta_mic, jdtheta_rad,jdtheta_SW, jdtheta_LW,                           &
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
    j_co2_emits, jsoil_thickness,                                              &
    jzseak_theta, jck_theta, jzseak_rho, jck_rho,                              &
    jti_mean, jti_sig, jfexp,                                                  &
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
! Description:
!   Routine to point atmosphere fields to the appropriate sections of D1.
!   After calling this subroutine, the fields can be used directly without
!   referring to D1
!
! Method:
!   Assuming SET_ATM_POINTERS has been called beforehand, this subroutine
!   points each field to an area of D1 starting at the corresponding
!   "jpointer" (at the first level) and ending at the "jpointer" (at the
!   last level) plus the size of a single level of that field.
!
!   Tracers are dealt with differently:   First the number of active tracers
!   is computed so that the correct sections of the corresponding tracer
!   jpointers can be used in pointing the tracer fields to D1.  If no tracers
!   are active then the fields are pointed to a dummy array
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
! Code description:
!   Language:  Fortran 90.
!   This code is written to UM programming standards UMDP3 vn 8.3.
!

! Subroutine arguments

REAL,    TARGET, INTENT(IN) :: d1(len_tot)

! Local variables

INTEGER :: nTracer ! loop counter over available tracers
INTEGER :: nActiveTracers ! number of tracers actually being used

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SET_ATM_FIELDS'

! End of header

! 1.0 Start of subroutine code; point fields to D1

!    atmospheric primary variables
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

sst  (tdims_s%i_start:tdims_s%i_end, tdims_s%j_start:tdims_s%j_end)       &
      => d1(jsst : jsst + field_length(theta_points,single_halo,1) -1)

Exner(pdims_s%i_start:pdims_s%i_end,                                      &
      pdims_s%j_start:pdims_s%j_end,                                      &
      pdims_s%k_start:pdims_s%k_end+1)                                    &
           => d1(jexner_rho_levels(pdims%k_start) :                       &
                 jexner_rho_levels(pdims%k_start) +                       &
                 field_length(theta_points,single_halo,pdims%k_len+1) -1)

Exner_Surf(pdims_s%i_start:pdims_s%i_end,                                 &
           pdims_s%j_start:pdims_s%j_end)                                 &
           => d1(jexnersurf : jexnersurf +                                &
                 field_length(theta_points,single_halo,1) -1)

dryrho(pdims_s%i_start:pdims_s%i_end,                                     &
       pdims_s%j_start:pdims_s%j_end,                                     &
       pdims_s%k_start:pdims_s%k_end)                                     &
           => d1(jdryrho(pdims%k_start) : jdryrho(pdims%k_start) +        &
                 field_length(theta_points,single_halo,                   &
                              pdims%k_len) -1)

etadot(tdims_s%i_start:tdims_s%i_end,                                     &
      tdims_s%j_start:tdims_s%j_end,                                      &
      tdims_s%k_start:tdims_s%k_end)                                      &
           => d1(jetadot(wdims%k_start) : jetadot(wdims%k_start) +        &
                 field_length(theta_points,single_halo,                   &
                              wdims%k_len) -1)

thetav(tdims_s%i_start:tdims_s%i_end,                                     &
       tdims_s%j_start:tdims_s%j_end,                                     &
       tdims_s%k_start:tdims_s%k_end)                                     &
           => d1(jthetav(tdims%k_start) : jthetav(tdims%k_start) +        &
                 field_length(theta_points,single_halo,                   &
                              tdims%k_len) -1)

psi_w_surf(tdims%i_start:tdims%i_end,                                     &
           tdims%j_start:tdims%j_end)                                     &
           => d1(jpsiws : jpsiws + field_length(theta_points,no_halo,1) -1)

psi_w_lid (tdims%i_start:tdims%i_end,                                     &
           tdims%j_start:tdims%j_end)                                     &
           => d1(jpsiwl : jpsiwl + field_length(theta_points,no_halo,1) -1)

m_v (tdims_s%i_start:tdims_s%i_end,                                       &
     tdims_s%j_start:tdims_s%j_end,                                       &
     tdims_s%k_start:tdims_s%k_end)                                       &
           => d1(jmv(tdims%k_start) : jmv(tdims%k_start) +                &
                 field_length(theta_points,single_halo,                   &
                         tdims%k_len) -1)

m_cl (tdims_s%i_start:tdims_s%i_end,                                      &
      tdims_s%j_start:tdims_s%j_end,                                      &
      tdims_s%k_start:tdims_s%k_end)                                      &
           => d1(jmcl(tdims%k_start) : jmcl(tdims%k_start) +              &
                 field_length(theta_points,single_halo,                   &
                              tdims%k_len) -1)

m_cf (tdims_s%i_start:tdims_s%i_end,                                      &
      tdims_s%j_start:tdims_s%j_end,                                      &
      tdims_s%k_start:tdims_s%k_end)                                      &
           => d1(jmcf(tdims%k_start) : jmcf(tdims%k_start) +              &
                 field_length(theta_points,single_halo,                   &
                              tdims%k_len) -1)

m_cf2 (tdims_s%i_start:tdims_s%i_end,                                     &
       tdims_s%j_start:tdims_s%j_end,                                     &
       tdims_s%k_start:tdims_s%k_end)                                     &
           => d1(jmcf2(tdims%k_start) : jmcf2(tdims%k_start) +            &
                 field_length(theta_points,single_halo,                   &
                              tdims%k_len) -1)

m_gr (tdims_s%i_start:tdims_s%i_end,                                      &
      tdims_s%j_start:tdims_s%j_end,                                      &
      tdims_s%k_start:tdims_s%k_end)                                      &
           => d1(jmgraup(tdims%k_start) : jmgraup(tdims%k_start) +        &
                 field_length(theta_points,single_halo,                   &
                         tdims%k_len) -1)

m_r  (tdims_s%i_start:tdims_s%i_end,                                      &
      tdims_s%j_start:tdims_s%j_end,                                      &
      tdims_s%k_start:tdims_s%k_end)                                      &
           => d1(jmrain(tdims%k_start) : jmrain(tdims%k_start) +          &
                 field_length(theta_points,single_halo,                   &
                              tdims%k_len) -1)


u   (udims_s%i_start:udims_s%i_end,                                       &
     udims_s%j_start:udims_s%j_end,                                       &
     udims_s%k_start:udims_s%k_end)                                       &
                           => d1(ju(udims%k_start) : ju(udims%k_start) +  &
        field_length(u_points,single_halo,udims%k_len) -1)

v   (vdims_s%i_start:vdims_s%i_end,                                       &
     vdims_s%j_start:vdims_s%j_end,                                       &
     vdims_s%k_start:vdims_s%k_end)                                       &
                           => d1(jv(vdims%k_start) : jv(vdims%k_start) +  &
        field_length(v_points,single_halo,vdims%k_len) -1)

theta(tdims_s%i_start:tdims_s%i_end,                                      &
      tdims_s%j_start:tdims_s%j_end,                                      &
      tdims_s%k_start:tdims_s%k_end)                                      &
                   => d1(jtheta(tdims%k_start) : jtheta(tdims%k_start) +  &
        field_length(theta_points,single_halo,tdims%k_len) -1)

q   (tdims_l%i_start:tdims_l%i_end,                                       &
     tdims_l%j_start:tdims_l%j_end,                                       &
     tdims_l%k_start:tdims_l%k_end)                                       &
                           => d1(jq(tdims%k_start) : jq(tdims%k_start) +  &
     field_length(theta_points,extended_halo,tdims%k_len) -1)

qcf (tdims_l%i_start:tdims_l%i_end,                                       &
     tdims_l%j_start:tdims_l%j_end,                                       &
     tdims_l%k_start:tdims_l%k_end)                                       &
                       => d1(jqcf(tdims%k_start) : jqcf(tdims%k_start) +  &
     field_length(theta_points,extended_halo,tdims%k_len) -1)

tstar(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)               &
          => d1(jtstar : jtstar+field_length(theta_points,no_halo,1) -1)

land (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)               &
          => d1(jland  : jland +field_length(theta_points,no_halo,1) -1)

orography(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)           &
          => d1(jorog  : jorog +field_length(theta_points,no_halo,1) -1)

w   (wdims_s%i_start:wdims_s%i_end,                                       &
     wdims_s%j_start:wdims_s%j_end,                                       &
     wdims_s%k_start:wdims_s%k_end)                                       &
                        => d1(jw(wdims_s%k_start) : jw(wdims_s%k_start) + &
       field_length(theta_points,single_halo,wdims%k_len) -1)

wetrho_r_sq_n(pdims_s%i_start:pdims_s%i_end,                              &
              pdims_s%j_start:pdims_s%j_end,                              &
              pdims_s%k_start:pdims_s%k_end)                              &
              => d1(jrho(pdims%k_start) : jrho(pdims%k_start)+            &
       field_length(theta_points,single_halo,pdims%k_len) -1)

rho  (pdims_s%i_start:pdims_s%i_end,                                      &
      pdims_s%j_start:pdims_s%j_end,                                      &
      pdims_s%k_start:pdims_s%k_end)                                      &
           => d1(jrho(pdims%k_start) : jrho(pdims%k_start)+               &
       field_length(theta_points,single_halo,pdims%k_len) -1)

qcl (tdims_l%i_start:tdims_l%i_end,                                       &
     tdims_l%j_start:tdims_l%j_end,                                       &
     tdims_l%k_start:tdims_l%k_end)                                       &
           => d1(jqcl(tdims%k_start) : jqcl(tdims%k_start)+               &
     field_length(theta_points,extended_halo,tdims%k_len) -1)

qcf2(tdims_l%i_start:tdims_l%i_end,                                       &
     tdims_l%j_start:tdims_l%j_end,                                       &
     tdims_l%k_start:tdims_l%k_end)                                       &
          => d1(jqcf2(tdims%k_start) : jqcf2(tdims%k_start)+              &
     field_length(theta_points,extended_halo,tdims%k_len) -1)

qrain(tdims_l%i_start:tdims_l%i_end,                                      &
      tdims_l%j_start:tdims_l%j_end,                                      &
      tdims_l%k_start:tdims_l%k_end)                                      &
          => d1(jqrain(tdims%k_start) : jqrain(tdims%k_start)+            &
     field_length(theta_points,extended_halo,tdims%k_len) -1)

qgraup(tdims_l%i_start:tdims_l%i_end,                                     &
       tdims_l%j_start:tdims_l%j_end,                                     &
       tdims_l%k_start:tdims_l%k_end)                                     &
          => d1(jqgraup(tdims%k_start) : jqgraup(tdims%k_start)+          &
     field_length(theta_points,extended_halo,tdims%k_len) -1)

IF ( n_casim_progs > 0 ) THEN

  ! CASIM Microphysics
  CloudNumber( tdims_l%i_start:tdims_l%i_end,                                &
               tdims_l%j_start:tdims_l%j_end,                                &
               tdims_l%k_start:tdims_l%k_end)                                &
  => D1(JCloudNumber(tdims%k_start) : JCloudNumber(tdims%k_start) +          &
        field_length(theta_points,extended_halo,tdims%k_len) -1)

  RainNumber( tdims_l%i_start:tdims_l%i_end,                                 &
              tdims_l%j_start:tdims_l%j_end,                                 &
              tdims_l%k_start:tdims_l%k_end)                                 &
  => D1(JRainNumber(tdims%k_start) : JRainNumber(tdims%k_start) +            &
        field_length(theta_points,extended_halo,tdims%k_len) -1)

  Rain3mom( tdims_l%i_start:tdims_l%i_end,                                   &
            tdims_l%j_start:tdims_l%j_end,                                   &
            tdims_l%k_start:tdims_l%k_end)                                   &
  => D1(JRain3mom(tdims%k_start) : JRain3mom(tdims%k_start)+                 &
        field_length(theta_points,extended_halo,tdims%k_len) -1)

  IceNumber( tdims_l%i_start:tdims_l%i_end,                                  &
             tdims_l%j_start:tdims_l%j_end,                                  &
             tdims_l%k_start:tdims_l%k_end)                                  &
  => D1( JIceNumber(tdims%k_start) : JIceNumber(tdims%k_start)+              &
         field_length(theta_points,extended_halo,tdims%k_len) -1)

  SnowNumber( tdims_l%i_start:tdims_l%i_end,                                 &
              tdims_l%j_start:tdims_l%j_end,                                 &
              tdims_l%k_start:tdims_l%k_end)                                 &
  => D1( JSnowNumber(tdims%k_start) : JSnowNumber(tdims%k_start)+            &
         field_length(theta_points,extended_halo,tdims%k_len) -1)

  Snow3mom( tdims_l%i_start:tdims_l%i_end,                                   &
            tdims_l%j_start:tdims_l%j_end,                                   &
            tdims_l%k_start:tdims_l%k_end)                                   &
  => D1( JSnow3mom(tdims%k_start) : JSnow3mom(tdims%k_start)+                &
         field_length(theta_points,extended_halo,tdims%k_len) -1)

  GraupNumber( tdims_l%i_start:tdims_l%i_end,                                &
               tdims_l%j_start:tdims_l%j_end,                                &
               tdims_l%k_start:tdims_l%k_end)                                &
  => D1( JGraupNumber(tdims%k_start) : JGraupNumber(tdims%k_start)+          &
         field_length(theta_points,extended_halo,tdims%k_len) -1)

  Graup3mom( tdims_l%i_start:tdims_l%i_end,                                  &
             tdims_l%j_start:tdims_l%j_end,                                  &
             tdims_l%k_start:tdims_l%k_end)                                  &
  => D1( JGraup3mom(tdims%k_start) : JGraup3mom(tdims%k_start)+              &
         field_length(theta_points,extended_halo,tdims%k_len) -1)

        ! CASIM activated aerosol prognostics
  ActiveSolLiquid( tdims_l%i_start:tdims_l%i_end,                            &
                   tdims_l%j_start:tdims_l%j_end,                            &
                   tdims_l%k_start:tdims_l%k_end)                            &
  => D1( JActiveSolLiquid(tdims%k_start) : JActiveSolLiquid(tdims%k_start)+  &
         field_length(theta_points,extended_halo,tdims%k_len) -1)

  ActiveSolRain( tdims_l%i_start:tdims_l%i_end,                              &
                 tdims_l%j_start:tdims_l%j_end,                              &
                 tdims_l%k_start:tdims_l%k_end)                              &
  => D1( JActiveSolRain(tdims%k_start) : JActiveSolRain(tdims%k_start)+      &
         field_length(theta_points,extended_halo,tdims%k_len) -1)

  ActiveInsolIce( tdims_l%i_start:tdims_l%i_end,                             &
                  tdims_l%j_start:tdims_l%j_end,                             &
                  tdims_l%k_start:tdims_l%k_end)                             &
  => D1( JActiveInSolIce(tdims%k_start) : JActiveInSolIce(tdims%k_start)+    &
           field_length(theta_points,extended_halo,tdims%k_len) -1)

  ActiveSolIce( tdims_l%i_start:tdims_l%i_end,                               &
                tdims_l%j_start:tdims_l%j_end,                               &
                tdims_l%k_start:tdims_l%k_end)                               &
  => D1( JActiveSolIce(tdims%k_start) : JActiveSolIce(tdims%k_start)+        &
         field_length(theta_points,extended_halo,tdims%k_len) -1)

  ActiveInsolLiquid( tdims_l%i_start:tdims_l%i_end,                           &
                     tdims_l%j_start:tdims_l%j_end,                           &
                     tdims_l%k_start:tdims_l%k_end)                           &
  => D1( JActiveInSolLiquid(tdims%k_start):JActiveInSolLiquid(tdims%k_start)+ &
         field_length(theta_points,extended_halo,tdims%k_len) -1)

  ActiveSolNumber( tdims_l%i_start:tdims_l%i_end,                             &
                   tdims_l%j_start:tdims_l%j_end,                             &
                   tdims_l%k_start:tdims_l%k_end)                             &
  => D1( JActiveSolNumber(tdims%k_start) : JActiveSolNumber(tdims%k_start)+   &
         field_length(theta_points,extended_halo,tdims%k_len) -1)

  ActiveInSolNumber( tdims_l%i_start:tdims_l%i_end,                           &
                     tdims_l%j_start:tdims_l%j_end,                           &
                     tdims_l%k_start:tdims_l%k_end)                           &
  => D1( JActiveInSolNumber(tdims%k_start):JActiveInSolNumber(tdims%k_start)+ &
         field_length(theta_points,extended_halo,tdims%k_len) -1)

END IF


exner_rho_levels(pdims_s%i_start:pdims_s%i_end,                           &
                 pdims_s%j_start:pdims_s%j_end,                           &
                 pdims_s%k_start:pdims_s%k_end+1)                         &
          => d1(jexner_rho_levels(pdims%k_start) :                        &
                jexner_rho_levels(pdims%k_start) +                        &
                 field_length(theta_points,single_halo,pdims%k_len+1) -1)

e_trb (tdims%i_start:tdims%i_end,                                         &
       tdims%j_start:tdims%j_end,                                         &
       tdims%k_start:tdims%k_end)                                         &
          => d1(je_trb(tdims%k_start) : je_trb(tdims%k_start)+            &
 field_length(theta_points,no_halo,tdims%k_len) -1)

tsq_trb(tdims%i_start:tdims%i_end,                                        &
        tdims%j_start:tdims%j_end,                                        &
        tdims%k_start:tdims%k_end)                                        &
          => d1(jtsq_trb(tdims%k_start) : jtsq_trb(tdims%k_start)+        &
                  field_length(theta_points,no_halo,tdims%k_len) -1)

qsq_trb(tdims%i_start:tdims%i_end,                                        &
        tdims%j_start:tdims%j_end,                                        &
        tdims%k_start:tdims%k_end)                                        &
          => d1(jqsq_trb(tdims%k_start) : jqsq_trb(tdims%k_start)+        &
                  field_length(theta_points,no_halo,tdims%k_len) -1)

cov_trb(tdims%i_start:tdims%i_end,                                        &
        tdims%j_start:tdims%j_end,                                        &
        tdims%k_start:tdims%k_end)                                        &
          => d1(jcov_trb(tdims%k_start) : jcov_trb(tdims%k_start)+        &
                  field_length(theta_points,no_halo,tdims%k_len) -1)

zhpar_shcu(tdims%i_start:tdims%i_end,                                     &
           tdims%j_start:tdims%j_end)                                     &
  => d1(jzhpar_shcu:jzhpar_shcu + field_length(theta_points,no_halo,1) -1)

dPV_rad(tdims%i_start:tdims%i_end,                                        &
        tdims%j_start:tdims%j_end,                                        &
        tdims%k_start:tdims%k_end)                                        &
        => d1(jdPV_rad(tdims%k_start) : jdPV_rad(tdims%k_start)  +        &
                  field_length(theta_points,no_halo,tdims%k_len) -1)

dPV_sw(tdims%i_start:tdims%i_end,                                         &
       tdims%j_start:tdims%j_end,                                         &
       tdims%k_start:tdims%k_end)                                         &
       => d1(jdPV_sw(tdims%k_start) : jdPV_sw(tdims%k_start)  +           &
                 field_length(theta_points,no_halo,tdims%k_len) -1)

dPV_lw(tdims%i_start:tdims%i_end,                                         &
       tdims%j_start:tdims%j_end,                                         &
       tdims%k_start:tdims%k_end)                                         &
       => d1(jdPV_lw(tdims%k_start) : jdPV_lw(tdims%k_start)  +           &
                 field_length(theta_points,no_halo,tdims%k_len) -1)

dPV_mic(tdims%i_start:tdims%i_end,                                        &
        tdims%j_start:tdims%j_end,                                        &
        tdims%k_start:tdims%k_end)                                        &
        => d1(jdPV_mic(tdims%k_start) : jdPV_mic(tdims%k_start)  +        &
                  field_length(theta_points,no_halo,tdims%k_len) -1)

dPV_gwd(tdims%i_start:tdims%i_end,                                        &
        tdims%j_start:tdims%j_end,                                        &
        tdims%k_start:tdims%k_end)                                        &
        => d1(jdPV_gwd(tdims%k_start) : jdPV_gwd(tdims%k_start)  +        &
                  field_length(theta_points,no_halo,tdims%k_len) -1)

dPV_ph1(tdims%i_start:tdims%i_end,                                        &
        tdims%j_start:tdims%j_end,                                        &
        tdims%k_start:tdims%k_end)                                        &
        =>  d1(jdPV_ph1(tdims%k_start) : jdPV_ph1(tdims%k_start)  +       &
                  field_length(theta_points,no_halo,tdims%k_len) -1)

dPV_conv(tdims%i_start:tdims%i_end,                                       &
         tdims%j_start:tdims%j_end,                                       &
         tdims%k_start:tdims%k_end)                                       &
         => d1(jdPV_conv(tdims%k_start) : jdPV_conv(tdims%k_start) +      &
                  field_length(theta_points,no_halo,tdims%k_len) -1)

dPV_bl(tdims%i_start:tdims%i_end,                                         &
       tdims%j_start:tdims%j_end,                                         &
       tdims%k_start:tdims%k_end)                                         &
        => d1(jdPV_bl(tdims%k_start) : jdPV_bl(tdims%k_start)  +          &
                  field_length(theta_points,no_halo,tdims%k_len) -1)

dPV_stph(tdims%i_start:tdims%i_end,                                       &
        tdims%j_start:tdims%j_end,                                        &
        tdims%k_start:tdims%k_end)                                        &
        => d1(jdPV_stph(tdims%k_start) : jdPV_stph(tdims%k_start)  +      &
                  field_length(theta_points,no_halo,tdims%k_len) -1)

dPV_cld(tdims%i_start:tdims%i_end,                                        &
        tdims%j_start:tdims%j_end,                                        &
        tdims%k_start:tdims%k_end)                                        &
        => d1(jdPV_cld(tdims%k_start) : jdPV_cld(tdims%k_start)  +        &
                  field_length(theta_points,no_halo,tdims%k_len) -1)

dPV_iau(tdims%i_start:tdims%i_end,                                        &
        tdims%j_start:tdims%j_end,                                        &
        tdims%k_start:tdims%k_end)                                        &
        => d1(jdPV_iau(tdims%k_start) : jdPV_iau(tdims%k_start)  +        &
                  field_length(theta_points,no_halo,tdims%k_len) -1)

dPV_nud(tdims%i_start:tdims%i_end,                                        &
        tdims%j_start:tdims%j_end,                                        &
        tdims%k_start:tdims%k_end)                                        &
        => d1(jdPV_nud(tdims%k_start) : jdPV_nud(tdims%k_start)  +        &
                  field_length(theta_points,no_halo,tdims%k_len) -1)

dPV_tot(tdims%i_start:tdims%i_end,                                        &
        tdims%j_start:tdims%j_end,                                        &
        tdims%k_start:tdims%k_end)                                        &
        => d1(jdPV_tot(tdims%k_start) : jdPV_tot(tdims%k_start)  +        &
                  field_length(theta_points,no_halo,tdims%k_len) -1)

dPV_adv(tdims%i_start:tdims%i_end,                                        &
        tdims%j_start:tdims%j_end,                                        &
        tdims%k_start:tdims%k_end)                                        &
        => d1(jdPV_adv(tdims%k_start) : jdPV_adv(tdims%k_start)  +        &
                  field_length(theta_points,no_halo,tdims%k_len) -1)

dPV_sol(tdims%i_start:tdims%i_end,                                        &
        tdims%j_start:tdims%j_end,                                        &
        tdims%k_start:tdims%k_end)                                        &
        => d1(jdPV_sol(tdims%k_start) : jdPV_sol(tdims%k_start)  +        &
                  field_length(theta_points,no_halo,tdims%k_len) -1)

dPV_mass(tdims%i_start:tdims%i_end,                                       &
        tdims%j_start:tdims%j_end,                                        &
        tdims%k_start:tdims%k_end)                                        &
        => d1(jdPV_mass(tdims%k_start) : jdPV_mass(tdims%k_start)  +      &
                  field_length(theta_points,no_halo,tdims%k_len) -1)

adv_only_PV(tdims%i_start:tdims%i_end,                                    &
            tdims%j_start:tdims%j_end,                                    &
            tdims%k_start:tdims%k_end)                                    &
        => d1(jadv_only_PV(tdims%k_start) : jadv_only_PV(tdims%k_start) + &
                  field_length(theta_points,no_halo,tdims%k_len) -1)

dtheta_0(tdims%i_start:tdims%i_end,                                       &
         tdims%j_start:tdims%j_end,                                       &
         tdims%k_start:tdims%k_end)                                       &
        => d1(jdtheta_0(tdims%k_start) : jdtheta_0(tdims%k_start)  +      &
                  field_length(theta_points,no_halo,tdims%k_len) -1)

dtheta_bl(tdims%i_start:tdims%i_end,                                      &
         tdims%j_start:tdims%j_end,                                       &
         tdims%k_start:tdims%k_end)                                       &
        => d1(jdtheta_bl(tdims%k_start) : jdtheta_bl(tdims%k_start)  +    &
                  field_length(theta_points,no_halo,tdims%k_len) -1)

dtheta_bl_mix(tdims%i_start:tdims%i_end,                                  &
              tdims%j_start:tdims%j_end,                                  &
              tdims%k_start:tdims%k_end)                                  &
        => d1(jdtheta_bl_mix(tdims%k_start):jdtheta_bl_mix(tdims%k_start)+&
                  field_length(theta_points,no_halo,tdims%k_len) -1)

dtheta_bl_LH(tdims%i_start:tdims%i_end,                                   &
             tdims%j_start:tdims%j_end,                                   &
             tdims%k_start:tdims%k_end)                                   &
        => d1(jdtheta_bl_LH(tdims%k_start):jdtheta_bl_LH(tdims%k_start) + &
                  field_length(theta_points,no_halo,tdims%k_len) -1)

dtheta_conv(tdims%i_start:tdims%i_end,                                    &
         tdims%j_start:tdims%j_end,                                       &
         tdims%k_start:tdims%k_end)                                       &
        => d1(jdtheta_conv(tdims%k_start) : jdtheta_conv(tdims%k_start) + &
                  field_length(theta_points,no_halo,tdims%k_len) -1)

dtheta_mic(tdims%i_start:tdims%i_end,                                     &
         tdims%j_start:tdims%j_end,                                       &
         tdims%k_start:tdims%k_end)                                       &
        => d1(jdtheta_mic(tdims%k_start) : jdtheta_mic(tdims%k_start)  +  &
                  field_length(theta_points,no_halo,tdims%k_len) -1)

dtheta_rad(tdims%i_start:tdims%i_end,                                     &
         tdims%j_start:tdims%j_end,                                       &
         tdims%k_start:tdims%k_end)                                       &
        => d1(jdtheta_rad(tdims%k_start) : jdtheta_rad(tdims%k_start)  +  &
                  field_length(theta_points,no_halo,tdims%k_len) -1)

dtheta_SW(tdims%i_start:tdims%i_end,                                      &
         tdims%j_start:tdims%j_end,                                       &
         tdims%k_start:tdims%k_end)                                       &
        => d1(jdtheta_SW(tdims%k_start) : jdtheta_SW(tdims%k_start)  +    &
                  field_length(theta_points,no_halo,tdims%k_len) -1)

dtheta_LW(tdims%i_start:tdims%i_end,                                      &
         tdims%j_start:tdims%j_end,                                       &
         tdims%k_start:tdims%k_end)                                       &
        => d1(jdtheta_LW(tdims%k_start) : jdtheta_LW(tdims%k_start)  +    &
                  field_length(theta_points,no_halo,tdims%k_len) -1)

dtheta_slow(tdims%i_start:tdims%i_end,                                    &
         tdims%j_start:tdims%j_end,                                       &
         tdims%k_start:tdims%k_end)                                       &
        => d1(jdtheta_slow(tdims%k_start) : jdtheta_slow(tdims%k_start) + &
                  field_length(theta_points,no_halo,tdims%k_len) -1)

dtheta_cld(tdims%i_start:tdims%i_end,                                     &
         tdims%j_start:tdims%j_end,                                       &
         tdims%k_start:tdims%k_end)                                       &
        => d1(jdtheta_cld(tdims%k_start) : jdtheta_cld(tdims%k_start)  +  &
                  field_length(theta_points,no_halo,tdims%k_len) -1)

flash_pot(tdims%i_start:tdims%i_end,                                      &
          tdims%j_start:tdims%j_end,                                      &
          tdims%k_start:tdims%k_end)                                      &
           => d1(jflash_pot(tdims%k_start) : jflash_pot(tdims%k_start) +  &
                      field_length(theta_points, no_halo, tdims%k_len) -1)

bl_w_var (tdims%i_start:tdims%i_end,                                      &
          tdims%j_start:tdims%j_end,                                      &
          tdims%k_start:tdims%k_end)                                      &
           => d1(jbl_w_var(1) : jbl_w_var(1) + &
 field_length(theta_points, no_halo, tdims%k_len) -1)

!    Coastal Tiling
frac_land  => d1(jfrac_land:jfrac_land  + &
                                    field_length(land_points,no_halo,1) -1)
tstar_land(tdims%i_start:tdims%i_end,                                     &
           tdims%j_start:tdims%j_end)                                     &
     => d1(jtstar_land:jtstar_land+ field_length(theta_points,no_halo,1)-1)

tstar_sea(tdims%i_start:tdims%i_end,                                      &
          tdims%j_start:tdims%j_end)                                      &
     => d1(jtstar_sea:jtstar_sea  + field_length(theta_points,no_halo,1)-1)

tstar_sice(tdims%i_start:tdims%i_end,                                     &
           tdims%j_start:tdims%j_end,1:1)                                 &
     => d1(jtstar_sice:jtstar_sice+ field_length(theta_points,no_halo,1)-1)

tstar_sice_cat(tdims%i_start:tdims%i_end,                                 &
               tdims%j_start:tdims%j_end,1:nice_use)                      &
     => d1(jtstar_sice_cat:jtstar_sice_cat+ &
                             field_length(theta_points,no_halo,nice_use)-1)

!    SeaIce and Land albedos
sice_alb (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)           &
     => d1(jsice_alb : jsice_alb + field_length(theta_points,no_halo,1)-1)
land_alb (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)           &
     => d1(jland_alb : jland_alb + field_length(theta_points,no_halo,1)-1)

!    Large-Scale hydrology
ti_mean   => d1(jti_mean:jti_mean+field_length(land_points,no_halo,1) -1)
ti_sig    => d1(jti_sig:jti_sig  +field_length(land_points,no_halo,1) -1)
fexp      => d1(jfexp:jfexp      +field_length(land_points,no_halo,1) -1)
gamma_int => d1(jgamma_int:jgamma_int + &
                                  field_length(land_points,no_halo,1) -1)
fsfc_sat  => d1(jfsfc_sat:jfsfc_sat   + &
                                  field_length(land_points,no_halo,1) -1)
f_wetland => d1(jf_wetland:jf_wetland + &
                                  field_length(land_points,no_halo,1) -1)
water_table => d1(jwater_table:jwater_table+ &
                                  field_length(land_points,no_halo,1) -1)

sthzw   => d1(jsthzw  : jsthzw  + field_length(land_points,no_halo,1) -1)
a_fsat  => d1(ja_fsat : ja_fsat + field_length(land_points,no_halo,1) -1)
c_fsat  => d1(jc_fsat : jc_fsat + field_length(land_points,no_halo,1) -1)
a_fwet  => d1(ja_fwet : ja_fwet + field_length(land_points,no_halo,1) -1)
c_fwet  => d1(jc_fwet : jc_fwet + field_length(land_points,no_halo,1) -1)

!    Optional atmospheric primary variables

zh    (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)              &
           => d1(jzh    : jzh    +field_length(theta_points,no_halo,1) -1)

ddmfx (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)              &
           => d1(jddmfx : jddmfx +field_length(theta_points,no_halo,1) -1)

  ALLOCATE(u_adv_nodump( (udims_l%i_len)*   &
                         (udims_l%j_len)*   &
                         (udims_l%k_len) ) )
  ALLOCATE(v_adv_nodump( (vdims_l%i_len)*   &
                         (vdims_l%j_len)*   &
                         (vdims_l%k_len) ) )
  ALLOCATE(w_adv_nodump( (wdims_l%i_len)*   &
                         (wdims_l%j_len)*   &
                         (wdims_l%k_len) ) )

  u_adv(udims_l%i_start:udims_l%i_end,                                      &
        udims_l%j_start:udims_l%j_end,                                      &
        udims_l%k_start:udims_l%k_end)  => u_adv_nodump
  v_adv(vdims_l%i_start:vdims_l%i_end,                                      &
        vdims_l%j_start:vdims_l%j_end,                                      &
        vdims_l%k_start:vdims_l%k_end)  => v_adv_nodump
  w_adv(wdims_l%i_start:wdims_l%i_end,                                      &
        wdims_l%j_start:wdims_l%j_end,                                      &
        wdims_l%k_start:wdims_l%k_end)  => w_adv_nodump

ntml   (pdims%i_start:pdims%i_end, pdims%j_start:pdims%j_end)               &
         => d1(jntml  : jntml +field_length(theta_points,no_halo,1) -1)
nbdsc  (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)               &
         => d1(jnbdsc : jnbdsc+field_length(theta_points,no_halo,1) -1)
ntdsc  (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)               &
         => d1(jntdsc : jntdsc+field_length(theta_points,no_halo,1) -1)
cumulus(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)               &
         => d1(jcumulus : jcumulus+field_length(theta_points,no_halo,1)-1)

t1_sd(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)                 &
          => d1(jt1_sd : jt1_sd+field_length(theta_points,no_halo,1) -1)

q1_sd(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)                 &
          => d1(jq1_sd : jq1_sd+field_length(theta_points,no_halo,1) -1)

! Note: no zeroth k level
cf_area(tdims%i_start:tdims%i_end,                                          &
        tdims%j_start:tdims%j_end,                                          &
                    1:tdims%k_end)                                          &
   => d1(jcf_area(1) : jcf_area(1) +                                        &
                    field_length(theta_points,no_halo      ,tdims%k_end) -1)

cf_bulk(tdims_l%i_start:tdims_l%i_end,                                      &
        tdims_l%j_start:tdims_l%j_end,                                      &
        tdims_l%k_start:tdims_l%k_end)                                      &
   => d1(jcf_bulk(tdims%k_start) : jcf_bulk(tdims%k_start) +                &
                    field_length(theta_points,extended_halo,tdims%k_len) -1)

cf_liquid(tdims_l%i_start:tdims_l%i_end,                                    &
          tdims_l%j_start:tdims_l%j_end,                                    &
          tdims_l%k_start:tdims_l%k_end)                                    &
   => d1(jcf_liquid(tdims%k_start) : jcf_liquid(tdims%k_start) +            &
                    field_length(theta_points,extended_halo,tdims%k_len) -1)

cf_frozen(tdims_l%i_start:tdims_l%i_end,                                    &
          tdims_l%j_start:tdims_l%j_end,                                    &
          tdims_l%k_start:tdims_l%k_end)                                    &
   => d1(jcf_frozen(tdims%k_start) : jcf_frozen(tdims%k_start) +            &
                    field_length(theta_points,extended_halo,tdims%k_len) -1)

! size of cca varies according to L_3D_CCA - n_cca_lev set in dervsize
cca (tdims%i_start:tdims%i_end,                                             &
     tdims%j_start:tdims%j_end, 1:n_cca_lev) =>                             &
    d1(jcca(1):jcca(1) + field_length(theta_points,no_halo,n_cca_lev) -1)

cca_dp (tdims%i_start:tdims%i_end,                                          &
        tdims%j_start:tdims%j_end, 1:n_cca_lev)                             &
           => d1(jcca_dp(1):jcca_dp(1) +                                    &
                         field_length(theta_points,no_halo,n_cca_lev) -1)

cca_md (tdims%i_start:tdims%i_end,                                          &
        tdims%j_start:tdims%j_end, 1:n_cca_lev)                             &
             => d1(jcca_md(1):jcca_md(1) +                                  &
                         field_length(theta_points,no_halo,n_cca_lev) -1)

cca_sh (tdims%i_start:tdims%i_end,                                          &
        tdims%j_start:tdims%j_end, 1:n_cca_lev)                             &
             => d1(jcca_sh(1):jcca_sh(1) +                                  &
                         field_length(theta_points,no_halo,n_cca_lev) -1)

ccb (tdims%i_start:tdims%i_end,                                             &
     tdims%j_start:tdims%j_end)                                             &
             => d1(jccb   : jccb  +field_length(theta_points,no_halo,1) -1)

cct  (tdims%i_start:tdims%i_end,                                            &
      tdims%j_start:tdims%j_end)                                            &
             => d1(jcct   : jcct  +field_length(theta_points,no_halo,1) -1)

cclwp  (tdims%i_start:tdims%i_end,                                          &
        tdims%j_start:tdims%j_end)                                          &
             => d1(jcclwp : jcclwp+field_length(theta_points,no_halo,1) -1)

deep_flag (tdims%i_start:tdims%i_end,                                       &
           tdims%j_start:tdims%j_end)                                       &
             => d1(jdeepflag : jdeepflag +                                  &
                               field_length(theta_points,no_halo,1) -1)
Past_precip (tdims%i_start:tdims%i_end,                                     &
             tdims%j_start:tdims%j_end)                                     &
             => d1(jpastprecip : jpastprecip +                              &
                               field_length(theta_points,no_halo,1) -1)
Past_conv_ht (tdims%i_start:tdims%i_end,                                    &
              tdims%j_start:tdims%j_end)                                    &
             => d1(jpastconvht : jpastconvht +                              &
                               field_length(theta_points,no_halo,1) -1)

canopy_water => d1(jcanopy_water : jcanopy_water+                           &
                                field_length(land_points,no_halo,1) -1)

lcbase (tdims%i_start:tdims%i_end,                                          &
        tdims%j_start:tdims%j_end)                                          &
             => d1(jlcbase : jlcbase + field_length(theta_points,no_halo,1) -1)

! Note: no zeroth k level
ccw_rad (tdims%i_start:tdims%i_end,                                         &
         tdims%j_start:tdims%j_end,                                         &
                     1:tdims%k_end)                                         &
             => d1(jccw_rad(1) : jccw_rad(1) +                              &
                         field_length(theta_points,no_halo,tdims%k_end) -1)

conv_prog_1(tdims_s%i_start:tdims_s%i_end,          &
            tdims_s%j_start:tdims_s%j_end,          &
            tdims_s%k_start:tdims_s%k_end)          &
             => d1(jconv_prog_1(tdims%k_start) : jconv_prog_1(tdims%k_start) + &
                       field_length(theta_points,single_halo,tdims%k_len) -1)
conv_prog_2(tdims_s%i_start:tdims_s%i_end,          &
            tdims_s%j_start:tdims_s%j_end,          &
            tdims_s%k_start:tdims_s%k_end)          &
             => d1(jconv_prog_2(tdims%k_start) : jconv_prog_2(tdims%k_start) + &
                       field_length(theta_points,single_halo,tdims%k_len) -1)
conv_prog_3(tdims_s%i_start:tdims_s%i_end,          &
            tdims_s%j_start:tdims_s%j_end,          &
            tdims_s%k_start:tdims_s%k_end)          &
             => d1(jconv_prog_3(tdims%k_start) : jconv_prog_3(tdims%k_start) + &
                       field_length(theta_points,single_halo,tdims%k_len) -1)
conv_prog_precip(tdims_s%i_start:tdims_s%i_end,          &
                 tdims_s%j_start:tdims_s%j_end,          &
                 tdims_s%k_start:tdims_s%k_end)          &
             => d1(jconv_prog_precip(tdims%k_start) :                    &
                   jconv_prog_precip(tdims%k_start) +                    &
                       field_length(theta_points,single_halo,tdims%k_len) -1)

totalppn(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)           &
             => d1(jtotalppn : jtotalppn +                               &
                       field_length(theta_points,no_halo,1) -1)

! Stochastic physics prognostics for BL perturbations
IF (i_pert_theta /= off .AND.                                           &
    i_pert_theta_type == pert_theta_correl_seq) THEN
  bl_pert_rand_fld(tdims%i_start:tdims%i_end,                           &
             tdims%j_start:tdims%j_end)                                 &
   => d1(jbl_pert_rand_fld : jbl_pert_rand_fld +                        &
                               field_length(theta_points,no_halo,1) -1)
  bl_pert_flag(tdims%i_start:tdims%i_end,                               &
             tdims%j_start:tdims%j_end)                                 &
   => d1(jbl_pert_flag : jbl_pert_flag +                                &
                               field_length(theta_points,no_halo,1) -1)
END IF

!    Secondary Fields in D1

exner_theta_levels(tdims_s%i_start:tdims_s%i_end,          &
                   tdims_s%j_start:tdims_s%j_end,          &
                   tdims_s%k_start:tdims_s%k_end)          &
 => d1(jexner_theta_levels(tdims%k_start):jexner_theta_levels(tdims%k_start)+ &
                         field_length(theta_points,single_halo,tdims%k_len) -1)

p (pdims_s%i_start:pdims_s%i_end,                                       &
   pdims_s%j_start:pdims_s%j_end,                                       &
   pdims_s%k_start:pdims_s%k_end) =>                                    &
            d1(jp(pdims%k_start):jp(pdims%k_start) +                    &
                  field_length(theta_points,single_halo, pdims%k_len) -1)

p_theta_levels(tdims_s%i_start:tdims_s%i_end,                                 &
               tdims_s%j_start:tdims_s%j_end,                                 &
               tdims_s%k_start:tdims_s%k_end) =>                              &
      d1(jp_theta_levels(tdims%k_start):jp_theta_levels(tdims%k_start)+       &
                     field_length(theta_points,single_halo,tdims%k_len) -1)

pstar(pdims%i_start:pdims%i_end, pdims%j_start:pdims%j_end) =>                &
      d1(jpstar : jpstar +field_length(theta_points,no_halo,1) -1)


sw_incs(rdims2%i_start:rdims2%i_end,                                          &
        rdims2%j_start:rdims2%j_end,                                          &
        rdims2%k_start:rdims2%k_end)                                          &
   => d1(jsw_incs(0) : jsw_incs(0) + &
                  field_length(theta_points, no_halo, rdims2%k_len)-1)

lw_incs(tdims%i_start:tdims%i_end,                                            &
        tdims%j_start:tdims%j_end,                                            &
        tdims%k_start:tdims%k_end)                                            &
   => d1(jlw_incs(0) : jlw_incs(0) + &
                  field_length(theta_points, no_halo, tdims%k_len)-1)

!    Direct PAR flux for STOCHEM
dirpar(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)                  &
        => d1(jdirpar : jdirpar+field_length(theta_points,no_halo,1) -1)

!    Soil Fields
smcl(1:land_field, 1:sm_levels)                                               &
     => d1(jsmcl(1):jsmcl(1)+ field_length(land_points,no_halo,sm_levels)-1)

deep_soil_temp(1:land_field, 1:sm_levels)                                     &
     => d1(j_deep_soil_temp(1):j_deep_soil_temp(1)+                           &
                         field_length(land_points,no_halo,sm_levels)-1)

vol_smc_wilt   => d1(jvol_smc_wilt:jvol_smc_wilt+ &
                                 field_length(land_points,no_halo,1)-1)
vol_smc_crit   => d1(jvol_smc_crit:jvol_smc_crit+ &
                                 field_length(land_points,no_halo,1)-1)
vol_smc_sat    => d1(jvol_smc_sat:jvol_smc_sat+ &
                                 field_length(land_points,no_halo,1)-1)
sat_soil_cond  => d1(jsat_soil_cond:jsat_soil_cond+ &
                                 field_length(land_points,no_halo,1)-1)
therm_cap  =>d1(jtherm_cap:jtherm_cap + &
                                field_length(land_points,no_halo,1) -1)
therm_cond =>d1(jtherm_cond:jtherm_cond+ &
                                 field_length(land_points,no_halo,1)-1)
clapp_horn =>d1(jclapp_horn:jclapp_horn+ &
                                 field_length(land_points,no_halo,1)-1)
z0m_soil => d1(jz0m_soil:jz0m_soil+ &
                                 field_length(land_points,no_halo,1)-1)
sat_soilw_suction => d1(jsat_soilw_suction:jsat_soilw_suction+ &
                                 field_length(land_points,no_halo,1)-1)

sthu(1:land_field, 1:sm_levels)                                               &
     => d1(jsthu(1):jsthu(1)+ field_length(land_points,no_halo,sm_levels)-1)

sthf(1:land_field, 1:sm_levels)                                               &
     => d1(jsthf(1):jsthf(1)+ field_length(land_points,no_halo,sm_levels)-1)

!    Roughness lenght of sea points
z0(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)                      &
     => d1(jz0:jz0           +field_length(theta_points,no_halo,1)-1)

gs   => d1(jgs:jgs           +field_length(land_points,no_halo,1)-1)

!    Orography Fields
orog_sil &
    => d1(jorog_sil : jorog_sil + field_length(land_points,no_halo,1)-1)
orog_ho2 &
    => d1(jorog_ho2 : jorog_ho2 + field_length(land_points,no_halo,1)-1)
orog_sd &
     => d1(jorog_sd  : jorog_sd + field_length(land_points,no_halo,1)-1)
orog_grad_x &
     => d1(jorog_grad_x : jorog_grad_x + field_length(land_points,no_halo,1)-1)
orog_grad_y  &
    => d1(jorog_grad_y : jorog_grad_y + field_length(land_points,no_halo,1)-1)
orog_unfilt &
     => d1(jorog_unfilt : jorog_unfilt + field_length(land_points,no_halo,1)-1)
orog_grad_xx &
    => d1(jorog_grad_xx : jorog_grad_xx + field_length(land_points,no_halo,1)-1)
orog_grad_xy &
    => d1(jorog_grad_xy : jorog_grad_xy + field_length(land_points,no_halo,1)-1)
orog_grad_yy &
    => d1(jorog_grad_yy : jorog_grad_yy + field_length(land_points,no_halo,1)-1)

!    Sea/Sea Ice Fields
u_sea(udims%i_start:udims%i_end, udims%j_start:udims%j_end)                   &
    => d1(ju_sea : ju_sea    +field_length(u_points,no_halo,1) -1)

v_sea(vdims%i_start:vdims%i_end, vdims%j_start:vdims%j_end)                   &
    => d1(jv_sea : jv_sea    +field_length(v_points,no_halo,1) -1)

ice_fraction(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end,1:1)        &
    => d1(jice_fraction  : jice_fraction +                                    &
                   field_length(theta_points_sea_only,no_halo,1) -1)

ice_thickness(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end,1:1)       &
    => d1(jice_thickness : jice_thickness +                                   &
                   field_length(theta_points_sea_only,no_halo,1) -1)

ti (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end,1:1)                 &
    => d1(jti : jti+field_length(theta_points_sea_only,no_halo,1) -1)

ice_fract_cat(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, 1:nice)   &
   => d1(jice_fract_cat : jice_fract_cat +                                    &
                         field_length(theta_points,no_halo,nice) -1)

ice_thick_cat(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, 1:nice)   &
   => d1(jice_thick_cat : jice_thick_cat +                                    &
                         field_length(theta_points,no_halo,nice) -1)

IF (l_sice_meltponds_cice) THEN
  pond_frac_cat(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, 1:nice) &
       => d1(jpond_frac_cat : jpond_frac_cat+                                 &
                        field_length(theta_points,no_halo,nice) -1)

  pond_depth_cat(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, 1:nice)&
       => d1(jpond_depth_cat : jpond_depth_cat+                               &
                        field_length(theta_points,no_halo,nice) -1)
ELSE
  pond_frac_cat => dummy_fld3
  pond_depth_cat => dummy_fld3
END IF

IF (l_saldep_freeze) THEN
  !    Sea surface freezing temperature
  sstfrz(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)             &
       => d1(jtfrz:jtfrz  + field_length(theta_points,no_halo,1)-1)
ELSE
  sstfrz => dummy_fld2
END IF

ti_cat (tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, 1:nice)          &
   => d1(jti_cat : jti_cat + field_length(theta_points,no_halo,nice) -1)

ice_k_cat (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, 1:nice)      &
   => d1(jice_k_cat : jice_k_cat + field_length(theta_points,no_halo,nice) -1)

u_0_p(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)                   &
      => d1(ju_0_p : ju_0_p+field_length(theta_points,no_halo,1) -1)

v_0_p(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)                   &
      => d1(jv_0_p : jv_0_p+field_length(theta_points,no_halo,1) -1)

chloro_sea(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)              &
      => d1(jchloro_sea : jchloro_sea + field_length(theta_points,no_halo,1)-1)

!    Snow Fields
snodep(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)                  &
      => d1(jsnodep : jsnodep+field_length(theta_points,no_halo,1) -1)

snodep_sea (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end,1:1)         &
      => d1(jsnodep_sea : jsnodep_sea+                                        &
                      field_length(theta_points_sea_only,no_halo,1) -1)

snodep_sea_cat(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, 1:nice)   &
      => d1(jsnodep_sea_cat : jsnodep_sea_cat+                                &
                            field_length(theta_points,no_halo,nice) -1)

! SNSOOT may not be used as of vn6.6
snsoot(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)                  &
              => d1(jsnsoot : jsnsoot+field_length(theta_points,no_halo,1) -1)

catch_snow(1:land_field, 1:ntiles) => d1(jcatch_snow : jcatch_snow +          &
                                 field_length(land_points,no_halo,ntiles) -1)

snow_grnd(1:land_field, 1:ntiles) => d1(jsnow_grnd : jsnow_grnd +             &
                           field_length(land_points,no_halo,ntiles) -1)


!    Decoupled screen temperatures
TScrnDcl_TILE(1:land_field,1:ntiles)                                         &
    => d1(JTScrnDcl_TILE : JTScrnDcl_TILE +                                  &
                                   field_length(land_points,no_halo,ntiles) -1)

TScrnDcl_SSI(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)           &
    => d1(JTScrnDcl_SSI : JTScrnDcl_SSI +                                    &
                                   field_length(theta_points,no_halo,1) -1)

tStbTrans   (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)           &
    => d1(JtStbTrans : JtStbTrans+ field_length(theta_points,no_halo,1) -1)

! convective cold pools
ux_ccp(pdims%i_start:pdims%i_end,          &
       pdims%j_start:pdims%j_end)          &
   => d1(jux_ccp : jux_ccp + field_length(theta_points,no_halo,1) -1)
uy_ccp(pdims%i_start:pdims%i_end,          &
       pdims%j_start:pdims%j_end)          &
   => d1(juy_ccp : juy_ccp + field_length(theta_points,no_halo,1) -1)
um_ccp(pdims%i_start:pdims%i_end,          &
       pdims%j_start:pdims%j_end)          &
   => d1(jum_ccp : jum_ccp + field_length(theta_points,no_halo,1) -1)
g_ccp(pdims%i_start:pdims%i_end,           &
       pdims%j_start:pdims%j_end)          &
   => d1(jg_ccp : jg_ccp + field_length(theta_points,no_halo,1) -1)
h_ccp(pdims%i_start:pdims%i_end,           &
       pdims%j_start:pdims%j_end)          &
   => d1(jh_ccp : jh_ccp + field_length(theta_points,no_halo,1) -1)
riso_ccp(pdims%i_start:pdims%i_end,        &
       pdims%j_start:pdims%j_end)          &
   => d1(jriso_ccp : jriso_ccp + field_length(theta_points,no_halo,1) -1)
rdir_ccp(pdims%i_start:pdims%i_end,        &
       pdims%j_start:pdims%j_end)          &
   => d1(jrdir_ccp : jrdir_ccp + field_length(theta_points,no_halo,1) -1)

!    OZONE (has extra surface level for V-AT-POLES)
IF (lexpand_ozone) THEN
  !      Ozone held as zonal averages, i.e. one value per row
  ! O3dims2 = tdims

  o3(1: rows *o3dims2%k_len)                                    &
     => d1(jozone(o3dims2%k_start):jozone(o3dims2%k_start)+     &
                                      rows * (o3dims2%k_len) -1)
ELSE
! o3dims%i_end = row_length*rows*ozone_levels

  o3(o3dims%i_start:o3dims%i_end)                           &
     => d1(jozone(o3dims2%k_start):jozone(o3dims2%k_start)+ &
                      field_length(ozone_points,no_halo,o3dims2%k_len) -1)
END IF

!    Tropopause-based Ozone
IF (tpps_ozone_levels > 0) THEN
  IF (lexpand_tpps_ozone) THEN
    tppsozone(1:rows*tpps_ozone_levels) => &
       d1(jtppsozone(o3dims2%k_start): jtppsozone(o3dims2%k_start) + &
                                                  rows*tpps_ozone_levels -1)
  ELSE
    tppsozone(1:row_length*rows*tpps_ozone_levels) => &
           d1(jtppsozone(o3dims2%k_start): jtppsozone(o3dims2%k_start)+ &
                    field_length(ozone_points,no_halo,tpps_ozone_levels) -1)
  END IF
ELSE
  tppsozone => dummy_field
END IF

!    Ozone tracer field and cariolle parameters
ozone_tracer(tdims_s%i_start:tdims_s%i_end,                             &
             tdims_s%j_start:tdims_s%j_end,                             &
             tdims_s%k_start:tdims_s%k_end) =>                          &
          d1(jozone_tracer(tdims%k_start):jozone_tracer(tdims%k_start)+ &
  field_length(theta_points,single_halo,tdims%k_len) -1)

o3_prod_loss(tdims%j_start:tdims%j_end, tdims%k_start:tdims%k_end) =>   &
          d1(jo3_prod_loss(tdims%k_start):jo3_prod_loss(tdims%k_start)+ &
                        (rows*(tdims%k_len) ) -1)

o3_p_l_vmr(tdims%j_start:tdims%j_end, tdims%k_start:tdims%k_end)   =>   &
          d1(jo3_p_l_vmr(tdims%k_start)  :jo3_p_l_vmr(tdims%k_start)+   &
                        (rows*(tdims%k_len) ) -1)

o3_vmr (tdims%j_start:tdims%j_end, tdims%k_start:tdims%k_end)      =>   &
          d1(jo3_vmr (tdims%k_start)     :jo3_vmr(tdims%k_start)+       &
                        (rows*(tdims%k_len) ) -1)

o3_p_l_temp (tdims%j_start:tdims%j_end, tdims%k_start:tdims%k_end) =>   &
          d1(jo3_p_l_temp(tdims%k_start) :jo3_p_l_temp(tdims%k_start)+  &
                        (rows*(tdims%k_len) ) -1)

o3_temp (tdims%j_start:tdims%j_end, tdims%k_start:tdims%k_end)     =>   &
          d1(jo3_temp(tdims%k_start)     :jo3_temp(tdims%k_start)+      &
                        (rows*(tdims%k_len) ) -1)

o3_p_l_colo3(tdims%j_start:tdims%j_end, tdims%k_start:tdims%k_end) =>   &
          d1(jo3_p_l_colo3(tdims%k_start):jo3_p_l_colo3(tdims%k_start)+ &
                        (rows*(tdims%k_len) ) -1)

o3_colo3 (tdims%j_start:tdims%j_end, tdims%k_start:tdims%k_end)    =>   &
          d1(jo3_colo3(tdims%k_start)    :jo3_colo3(tdims%k_start)+     &
                        (rows*(tdims%k_len) ) -1)

!    Sources and Aerosol Ancillaries
murk_source(tdims%i_start:tdims%i_end,                                  &
            tdims%j_start:tdims%j_end,                                  &
            tdims%k_start:tdims%k_end)                                  &
     => d1(jmurk_source(tdims%k_start) : jmurk_source(tdims%k_start) +  &
                      field_length(theta_points,no_halo,tdims%k_len) -1)

so2_em(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)            &
     => d1(jso2_em : jso2_em + field_length(theta_points,no_halo,1) -1)

dms_em(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)            &
     => d1(jdms_em : jdms_em + field_length(theta_points,no_halo,1) -1)

murk   (tdims_s%i_start:tdims_s%i_end,                                  &
        tdims_s%j_start:tdims_s%j_end,                                  &
        tdims_s%k_start:tdims_s%k_end)                                  &
     => d1(jmurk(tdims%k_start) : jmurk(tdims%k_start) +                &
                  field_length(theta_points,single_halo,tdims%k_len) -1)

!    Sulphur cycle
so2(tdims_s%i_start:tdims_s%i_end,                                      &
    tdims_s%j_start:tdims_s%j_end,                                      &
    tdims_s%k_start:tdims_s%k_end)                                      &
     => d1(jso2(tdims%k_start) : jso2(tdims%k_start) +                  &
                  field_length(theta_points,single_halo,tdims%k_len) -1)

dms(tdims_s%i_start:tdims_s%i_end,                                      &
    tdims_s%j_start:tdims_s%j_end,                                      &
    tdims_s%k_start:tdims_s%k_end)                                      &
     => d1(jdms(tdims%k_start) : jdms(tdims%k_start) +                  &
                  field_length(theta_points,single_halo,tdims%k_len) -1)

so4_aitken(tdims_s%i_start:tdims_s%i_end,                               &
           tdims_s%j_start:tdims_s%j_end,                               &
           tdims_s%k_start:tdims_s%k_end)                               &
     => d1(jso4_aitken(tdims%k_start) : jso4_aitken(tdims%k_start) +    &
                  field_length(theta_points,single_halo,tdims%k_len) -1)

so4_accu  (tdims_s%i_start:tdims_s%i_end,                               &
           tdims_s%j_start:tdims_s%j_end,                               &
           tdims_s%k_start:tdims_s%k_end)                               &
    => d1(jso4_accu(tdims%k_start) : jso4_accu(tdims%k_start) +         &
                  field_length(theta_points,single_halo,tdims%k_len) -1)

so4_diss  (tdims_s%i_start:tdims_s%i_end,                               &
           tdims_s%j_start:tdims_s%j_end,                               &
           tdims_s%k_start:tdims_s%k_end)                               &
    => d1(jso4_diss(tdims%k_start) : jso4_diss(tdims%k_start) +         &
                  field_length(theta_points,single_halo,tdims%k_len) -1)

h2o2      (tdims_s%i_start:tdims_s%i_end,                               &
           tdims_s%j_start:tdims_s%j_end,                               &
           tdims_s%k_start:tdims_s%k_end)                               &
        => d1(jh2o2(tdims%k_start) : jh2o2(tdims%k_start) +             &
                  field_length(theta_points,single_halo,tdims%k_len) -1)

nh3       (tdims_s%i_start:tdims_s%i_end,                               &
           tdims_s%j_start:tdims_s%j_end,                               &
           tdims_s%k_start:tdims_s%k_end)                               &
         => d1(jnh3(tdims%k_start) : jnh3(tdims%k_start) +              &
                  field_length(theta_points,single_halo,tdims%k_len) -1)

soot_new  (tdims_s%i_start:tdims_s%i_end,                               &
           tdims_s%j_start:tdims_s%j_end,                               &
          tdims_s%k_start:tdims_s%k_end)                                &
    => d1(jsoot_new(tdims%k_start) : jsoot_new(tdims%k_start) +         &
                  field_length(theta_points,single_halo,tdims%k_len) -1)

soot_agd  (tdims_s%i_start:tdims_s%i_end,                               &
           tdims_s%j_start:tdims_s%j_end,                               &
           tdims_s%k_start:tdims_s%k_end)                               &
    => d1(jsoot_agd(tdims%k_start) : jsoot_agd(tdims%k_start) +         &
                  field_length(theta_points,single_halo,tdims%k_len) -1)

soot_cld  (tdims_s%i_start:tdims_s%i_end,                               &
           tdims_s%j_start:tdims_s%j_end,                               &
           tdims_s%k_start:tdims_s%k_end)                               &
    => d1(jsoot_cld(tdims%k_start) : jsoot_cld(tdims%k_start) +         &
                  field_length(theta_points,single_halo,tdims%k_len) -1)

bmass_new (tdims_s%i_start:tdims_s%i_end,                               &
           tdims_s%j_start:tdims_s%j_end,                               &
           tdims_s%k_start:tdims_s%k_end)                               &
   => d1(jbmass_new(tdims%k_start) :jbmass_new(tdims%k_start) +         &
                  field_length(theta_points,single_halo,tdims%k_len) -1)

bmass_agd (tdims_s%i_start:tdims_s%i_end,                               &
           tdims_s%j_start:tdims_s%j_end,                               &
           tdims_s%k_start:tdims_s%k_end)                               &
   => d1(jbmass_agd(tdims%k_start) : jbmass_agd(tdims%k_start) +        &
                  field_length(theta_points,single_halo,tdims%k_len) -1)

bmass_cld (tdims_s%i_start:tdims_s%i_end,                               &
           tdims_s%j_start:tdims_s%j_end,                               &
           tdims_s%k_start:tdims_s%k_end)                               &
   => d1(jbmass_cld(tdims%k_start) : jbmass_cld(tdims%k_start) +        &
                  field_length(theta_points,single_halo,tdims%k_len) -1)

ocff_new  (tdims_s%i_start:tdims_s%i_end,                               &
           tdims_s%j_start:tdims_s%j_end,                               &
           tdims_s%k_start:tdims_s%k_end)                               &
    => d1(jocff_new(tdims%k_start) : jocff_new(tdims%k_start) +         &
                  field_length(theta_points,single_halo,tdims%k_len) -1)

ocff_agd  (tdims_s%i_start:tdims_s%i_end,                               &
           tdims_s%j_start:tdims_s%j_end,                               &
           tdims_s%k_start:tdims_s%k_end)                               &
    => d1(jocff_agd(tdims%k_start) : jocff_agd(tdims%k_start) +         &
                  field_length(theta_points,single_halo,tdims%k_len) -1)

ocff_cld  (tdims_s%i_start:tdims_s%i_end,                               &
           tdims_s%j_start:tdims_s%j_end,                               &
           tdims_s%k_start:tdims_s%k_end)                               &
    => d1(jocff_cld(tdims%k_start) : jocff_cld(tdims%k_start) +         &
                  field_length(theta_points,single_halo,tdims%k_len) -1)

so2_natem(tdims%i_start:tdims%i_end,                                    &
          tdims%j_start:tdims%j_end,                                    &
          tdims%k_start:tdims%k_end)                                    &
    => d1(jso2_natem(tdims%k_start) : jso2_natem(tdims%k_start) +       &
                  field_length(theta_points,no_halo,tdims%k_len) -1)

oh       (tdims%i_start:tdims%i_end,                                    &
          tdims%j_start:tdims%j_end,                                    &
          tdims%k_start:tdims%k_end)                                    &
    => d1(joh(tdims%k_start) : joh(tdims%k_start) +                     &
                  field_length(theta_points,no_halo,tdims%k_len) -1)

ho2      (tdims%i_start:tdims%i_end,                                    &
          tdims%j_start:tdims%j_end,                                    &
          tdims%k_start:tdims%k_end)                                    &
    => d1(jho2(tdims%k_start) : jho2(tdims%k_start) +                   &
                  field_length(theta_points,no_halo,tdims%k_len) -1)

h2o2_limit(tdims%i_start:tdims%i_end,                                   &
           tdims%j_start:tdims%j_end,                                   &
           tdims%k_start:tdims%k_end)                                   &
    => d1(jh2o2_limit(tdims%k_start) :jh2o2_limit(tdims%k_start) +      &
                  field_length(theta_points,no_halo,tdims%k_len) -1)
o3_chem   (tdims%i_start:tdims%i_end,                                   &
           tdims%j_start:tdims%j_end,                                   &
           tdims%k_start:tdims%k_end)                                   &
    => d1(jo3_chem(tdims%k_start) : jo3_chem(tdims%k_start) +           &
                  field_length(theta_points,no_halo,tdims%k_len) -1)

so2_hilem(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)         &
   => d1(jso2_hilem: jso2_hilem+ field_length(theta_points,no_halo,1) -1)

nh3_em(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)            &
   => d1(jnh3_em: jnh3_em + field_length(theta_points,no_halo,1) -1)

soot_em(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)           &
   => d1(jsoot_em: jsoot_em + field_length(theta_points,no_halo,1) -1)

soot_hilem(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)        &
   => d1(jsoot_hilem: jsoot_hilem + field_length(theta_points,no_halo,1) -1)

bmass_em(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)          &
   => d1(jbmass_em: jbmass_em + field_length(theta_points,no_halo,1) -1)

bmass_hilem(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)       &
   => d1(jbmass_hilem: jbmass_hilem + field_length(theta_points,no_halo,1) -1)

bmass_hilem_h1(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)    &
   => d1(jbmass_hilem_h1:                                               &
         jbmass_hilem_h1 + field_length(theta_points,no_halo,1) -1)

bmass_hilem_h2(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)    &
   => d1(jbmass_hilem_h2:                                               &
         jbmass_hilem_h2 + field_length(theta_points,no_halo,1) -1)

ocff_em(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)           &
   => d1(jocff_em: jocff_em + field_length(theta_points,no_halo,1) -1)

ocff_hilem(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)        &
   => d1(jocff_hilem: jocff_hilem + field_length(theta_points,no_halo,1) -1)

dms_conc(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)          &
   => d1(jdms_conc: jdms_conc + field_length(theta_points,no_halo,1) -1)

!
! Ammonium nitrate scheme:
nitr_acc (tdims_s%i_start:tdims_s%i_end,                                &
          tdims_s%j_start:tdims_s%j_end,                                &
          tdims_s%k_start:tdims_s%k_end)                                &
  => d1(jnitr_acc(tdims%k_start) : jnitr_acc(tdims%k_start) +           &
                   field_length(theta_points,single_halo,tdims%k_len) -1)
nitr_diss(tdims_s%i_start:tdims_s%i_end,                                &
          tdims_s%j_start:tdims_s%j_end,                                &
          tdims_s%k_start:tdims_s%k_end)                                &
 => d1(jnitr_diss(tdims%k_start) : jnitr_diss(tdims%k_start) +          &
                   field_length(theta_points,single_halo,tdims%k_len) -1)

!
! Aerosol climatologies
arclbiog_bg(tdims%i_start:tdims%i_end,                     &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end) =>                  &
          d1(jarclbiog_bg(tdims%k_start) : jarclbiog_bg(tdims%k_start) + &
                       field_length(theta_points,no_halo,tdims%k_len) -1)

arclbiom_fr(tdims%i_start:tdims%i_end,                     &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end) =>                  &
          d1(jarclbiom_fr(tdims%k_start) : jarclbiom_fr(tdims%k_start) + &
                       field_length(theta_points,no_halo,tdims%k_len) -1)

arclbiom_ag(tdims%i_start:tdims%i_end,                     &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end) =>                  &
          d1(jarclbiom_ag(tdims%k_start) : jarclbiom_ag(tdims%k_start) + &
                       field_length(theta_points,no_halo,tdims%k_len) -1)

arclbiom_ic(tdims%i_start:tdims%i_end,                     &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end) =>                  &
          d1(jarclbiom_ic(tdims%k_start) : jarclbiom_ic(tdims%k_start) + &
                       field_length(theta_points,no_halo,tdims%k_len) -1)

arclblck_fr(tdims%i_start:tdims%i_end,                     &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end) =>                  &
          d1(jarclblck_fr(tdims%k_start) : jarclblck_fr(tdims%k_start) + &
                       field_length(theta_points,no_halo,tdims%k_len) -1)

arclblck_ag(tdims%i_start:tdims%i_end,                     &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end) =>                  &
          d1(jarclblck_ag(tdims%k_start) : jarclblck_ag(tdims%k_start) + &
                       field_length(theta_points,no_halo,tdims%k_len) -1)

arclsslt_fi(tdims%i_start:tdims%i_end,                     &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end) =>                  &
          d1(jarclsslt_fi(tdims%k_start) : jarclsslt_fi(tdims%k_start) + &
                       field_length(theta_points,no_halo,tdims%k_len) -1)

arclsslt_jt(tdims%i_start:tdims%i_end,                     &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end) =>                  &
          d1(jarclsslt_jt(tdims%k_start) : jarclsslt_jt(tdims%k_start) + &
                       field_length(theta_points,no_halo,tdims%k_len) -1)

arclsulp_ac(tdims%i_start:tdims%i_end,                     &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end) =>                  &
          d1(jarclsulp_ac(tdims%k_start) : jarclsulp_ac(tdims%k_start) + &
                       field_length(theta_points,no_halo,tdims%k_len) -1)

arclsulp_ak(tdims%i_start:tdims%i_end,                     &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end) =>                  &
          d1(jarclsulp_ak(tdims%k_start) : jarclsulp_ak(tdims%k_start) + &
                       field_length(theta_points,no_halo,tdims%k_len) -1)

arclsulp_di(tdims%i_start:tdims%i_end,                     &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end) =>                  &
          d1(jarclsulp_di(tdims%k_start) : jarclsulp_di(tdims%k_start) + &
                       field_length(theta_points,no_halo,tdims%k_len) -1)

arcldust_b1(tdims%i_start:tdims%i_end,                     &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end) =>                  &
          d1(jarcldust_b1(tdims%k_start) : jarcldust_b1(tdims%k_start) + &
                       field_length(theta_points,no_halo,tdims%k_len) -1)

arcldust_b2(tdims%i_start:tdims%i_end,                     &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end) =>                  &
          d1(jarcldust_b2(tdims%k_start) : jarcldust_b2(tdims%k_start) + &
                       field_length(theta_points,no_halo,tdims%k_len) -1)

arcldust_b3(tdims%i_start:tdims%i_end,                     &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end) =>                  &
          d1(jarcldust_b3(tdims%k_start) : jarcldust_b3(tdims%k_start) + &
                       field_length(theta_points,no_halo,tdims%k_len) -1)

arcldust_b4(tdims%i_start:tdims%i_end,                     &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end) =>                  &
          d1(jarcldust_b4(tdims%k_start) : jarcldust_b4(tdims%k_start) + &
                       field_length(theta_points,no_halo,tdims%k_len) -1)

arcldust_b5(tdims%i_start:tdims%i_end,                     &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end) =>                  &
          d1(jarcldust_b5(tdims%k_start) : jarcldust_b5(tdims%k_start) + &
                       field_length(theta_points,no_halo,tdims%k_len) -1)

arcldust_b6(tdims%i_start:tdims%i_end,                     &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end) =>                  &
          d1(jarcldust_b6(tdims%k_start) : jarcldust_b6(tdims%k_start) + &
                       field_length(theta_points,no_halo,tdims%k_len) -1)

arclocff_fr(tdims%i_start:tdims%i_end,                     &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end) =>                  &
          d1(jarclocff_fr(tdims%k_start) : jarclocff_fr(tdims%k_start) + &
                       field_length(theta_points,no_halo,tdims%k_len) -1)

arclocff_ag(tdims%i_start:tdims%i_end,                     &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end) =>                  &
          d1(jarclocff_ag(tdims%k_start) : jarclocff_ag(tdims%k_start) + &
                       field_length(theta_points,no_halo,tdims%k_len) -1)

arclocff_ic(tdims%i_start:tdims%i_end,                     &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end) =>                  &
          d1(jarclocff_ic(tdims%k_start) : jarclocff_ic(tdims%k_start) + &
                       field_length(theta_points,no_halo,tdims%k_len) -1)

arcldlta_dl(tdims%i_start:tdims%i_end,                     &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end) =>                  &
          d1(jarcldlta_dl(tdims%k_start) : jarcldlta_dl(tdims%k_start) + &
                       field_length(theta_points,no_halo,tdims%k_len) -1)



! Section 54 GLOMAP_CLIM Aerosol climatologies

gc_nd_nuc_sol(tdims%i_start:tdims%i_end,                   &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end)  =>                 &
          d1(jgc_nd_nuc_sol(tdims%k_start) : jgc_nd_nuc_sol(tdims%k_start) + &
          field_length(theta_points,no_halo,tdims%k_len) -1)

gc_nuc_sol_su(tdims%i_start:tdims%i_end,                   &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end)  =>                 &
          d1(jgc_nuc_sol_su(tdims%k_start) : jgc_nuc_sol_su(tdims%k_start) + &
          field_length(theta_points,no_halo,tdims%k_len) -1)

gc_nuc_sol_oc(tdims%i_start:tdims%i_end,                   &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end)  =>                 &
          d1(jgc_nuc_sol_oc(tdims%k_start) : jgc_nuc_sol_oc(tdims%k_start) + &
          field_length(theta_points,no_halo,tdims%k_len) -1)

gc_nd_ait_sol(tdims%i_start:tdims%i_end,                   &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end)  =>                 &
          d1(jgc_nd_ait_sol(tdims%k_start) : jgc_nd_ait_sol(tdims%k_start) + &
          field_length(theta_points,no_halo,tdims%k_len) -1)

gc_ait_sol_su(tdims%i_start:tdims%i_end,                   &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end)  =>                 &
          d1(jgc_ait_sol_su(tdims%k_start) : jgc_ait_sol_su(tdims%k_start) + &
          field_length(theta_points,no_halo,tdims%k_len) -1)

gc_ait_sol_bc(tdims%i_start:tdims%i_end,                   &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end)  =>                 &
          d1(jgc_ait_sol_bc(tdims%k_start) : jgc_ait_sol_bc(tdims%k_start) + &
          field_length(theta_points,no_halo,tdims%k_len) -1)

gc_ait_sol_oc(tdims%i_start:tdims%i_end,                   &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end)  =>                 &
          d1(jgc_ait_sol_oc(tdims%k_start) : jgc_ait_sol_oc(tdims%k_start) + &
          field_length(theta_points,no_halo,tdims%k_len) -1)

gc_nd_acc_sol(tdims%i_start:tdims%i_end,                   &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end)  =>                 &
          d1(jgc_nd_acc_sol(tdims%k_start) : jgc_nd_acc_sol(tdims%k_start) + &
          field_length(theta_points,no_halo,tdims%k_len) -1)

gc_acc_sol_su(tdims%i_start:tdims%i_end,                   &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end)  =>                 &
          d1(jgc_acc_sol_su(tdims%k_start) : jgc_acc_sol_su(tdims%k_start) + &
          field_length(theta_points,no_halo,tdims%k_len) -1)

gc_acc_sol_bc(tdims%i_start:tdims%i_end,                   &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end)  =>                 &
          d1(jgc_acc_sol_bc(tdims%k_start) : jgc_acc_sol_bc(tdims%k_start) + &
          field_length(theta_points,no_halo,tdims%k_len) -1)

gc_acc_sol_oc(tdims%i_start:tdims%i_end,                   &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end)  =>                 &
          d1(jgc_acc_sol_oc(tdims%k_start) : jgc_acc_sol_oc(tdims%k_start) + &
          field_length(theta_points,no_halo,tdims%k_len) -1)

gc_acc_sol_ss(tdims%i_start:tdims%i_end,                   &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end)  =>                 &
          d1(jgc_acc_sol_ss(tdims%k_start) : jgc_acc_sol_ss(tdims%k_start) + &
          field_length(theta_points,no_halo,tdims%k_len) -1)

gc_nd_cor_sol(tdims%i_start:tdims%i_end,                   &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end)  =>                 &
          d1(jgc_nd_cor_sol(tdims%k_start) : jgc_nd_cor_sol(tdims%k_start) + &
          field_length(theta_points,no_halo,tdims%k_len) -1)

gc_cor_sol_su(tdims%i_start:tdims%i_end,                   &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end)  =>                 &
          d1(jgc_cor_sol_su(tdims%k_start) : jgc_cor_sol_su(tdims%k_start) + &
          field_length(theta_points,no_halo,tdims%k_len) -1)

gc_cor_sol_bc(tdims%i_start:tdims%i_end,                   &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end)  =>                 &
          d1(jgc_cor_sol_bc(tdims%k_start) : jgc_cor_sol_bc(tdims%k_start) + &
          field_length(theta_points,no_halo,tdims%k_len) -1)

gc_cor_sol_oc(tdims%i_start:tdims%i_end,                   &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end)  =>                 &
          d1(jgc_cor_sol_oc(tdims%k_start) : jgc_cor_sol_oc(tdims%k_start) + &
          field_length(theta_points,no_halo,tdims%k_len) -1)

gc_cor_sol_ss(tdims%i_start:tdims%i_end,                   &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end)  =>                 &
          d1(jgc_cor_sol_ss(tdims%k_start) : jgc_cor_sol_ss(tdims%k_start) + &
          field_length(theta_points,no_halo,tdims%k_len) -1)

gc_nd_ait_ins(tdims%i_start:tdims%i_end,                   &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end)  =>                 &
          d1(jgc_nd_ait_ins(tdims%k_start) : jgc_nd_ait_ins(tdims%k_start) + &
          field_length(theta_points,no_halo,tdims%k_len) -1)

gc_ait_ins_bc(tdims%i_start:tdims%i_end,                   &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end)  =>                 &
          d1(jgc_ait_ins_bc(tdims%k_start) : jgc_ait_ins_bc(tdims%k_start) + &
          field_length(theta_points,no_halo,tdims%k_len) -1)

gc_ait_ins_oc(tdims%i_start:tdims%i_end,                   &
            tdims%j_start:tdims%j_end,                     &
            tdims%k_start:tdims%k_end)  =>                 &
          d1(jgc_ait_ins_oc(tdims%k_start) : jgc_ait_ins_oc(tdims%k_start) + &
          field_length(theta_points,no_halo,tdims%k_len) -1)

!    Mineral Dust Scheme
soil_clay(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)            &
    => d1(jsoil_clay : jsoil_clay+ field_length(theta_points,no_halo,1) -1)
soil_silt(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)            &
    => d1(jsoil_silt : jsoil_silt+ field_length(theta_points,no_halo,1) -1)
soil_sand(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)            &
    => d1(jsoil_sand : jsoil_sand+ field_length(theta_points,no_halo,1) -1)
dust_mrel1(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)            &
   => d1(jdust_mrel1 : jdust_mrel1+ field_length(theta_points,no_halo,1) -1)
dust_mrel2(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)            &
   => d1(jdust_mrel2 : jdust_mrel2+ field_length(theta_points,no_halo,1) -1)
dust_mrel3(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)            &
   => d1(jdust_mrel3 : jdust_mrel3+ field_length(theta_points,no_halo,1) -1)
dust_mrel4(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)            &
   => d1(jdust_mrel4 : jdust_mrel4+ field_length(theta_points,no_halo,1) -1)
dust_mrel5(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)            &
   => d1(jdust_mrel5 : jdust_mrel5+ field_length(theta_points,no_halo,1) -1)
dust_mrel6(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)            &
   => d1(jdust_mrel6 : jdust_mrel6+ field_length(theta_points,no_halo,1) -1)

dust_div1(tdims_s%i_start:tdims_s%i_end,                                    &
          tdims_s%j_start:tdims_s%j_end,                                    &
          tdims_s%k_start:tdims_s%k_end)                                    &
     => d1(jdust_div1(tdims%k_start) : jdust_div1(tdims%k_start)+           &
                      field_length(theta_points,single_halo,tdims%k_len) -1)

dust_div2(tdims_s%i_start:tdims_s%i_end,                                    &
          tdims_s%j_start:tdims_s%j_end,                                    &
          tdims_s%k_start:tdims_s%k_end)                                    &
     => d1(jdust_div2(tdims%k_start) : jdust_div2(tdims%k_start)+           &
                      field_length(theta_points,single_halo,tdims%k_len) -1)

dust_div3(tdims_s%i_start:tdims_s%i_end,                                    &
          tdims_s%j_start:tdims_s%j_end,                                    &
          tdims_s%k_start:tdims_s%k_end)                                    &
     => d1(jdust_div3(tdims%k_start) : jdust_div3(tdims%k_start)+           &
                      field_length(theta_points,single_halo,tdims%k_len) -1)

dust_div4(tdims_s%i_start:tdims_s%i_end,                                    &
          tdims_s%j_start:tdims_s%j_end,                                    &
          tdims_s%k_start:tdims_s%k_end)                                    &
      => d1(jdust_div4(tdims%k_start) : jdust_div4(tdims%k_start)+          &
                      field_length(theta_points,single_halo,tdims%k_len) -1)

dust_div5(tdims_s%i_start:tdims_s%i_end,                                    &
          tdims_s%j_start:tdims_s%j_end,                                    &
          tdims_s%k_start:tdims_s%k_end)                                    &
      => d1(jdust_div5(tdims%k_start) : jdust_div5(tdims%k_start)+          &
                      field_length(theta_points,single_halo,tdims%k_len) -1)

dust_div6(tdims_s%i_start:tdims_s%i_end,                                    &
          tdims_s%j_start:tdims_s%j_end,                                    &
          tdims_s%k_start:tdims_s%k_end)                                    &
      => d1(jdust_div6(tdims%k_start) : jdust_div6(tdims%k_start)+          &
                      field_length(theta_points,single_halo,tdims%k_len) -1)

!    Carbon Cycle
triffid_co2_d1(1:land_field) =>  d1(j_triffid_co2:j_triffid_co2 +           &
                                  field_length(land_points,no_halo,1)-1)

co2flux (tdims%i_start:tdims%i_end,                                        &
         tdims%j_start:tdims%j_end)                                        &
    => d1(j_co2flux:j_co2flux + field_length(theta_points,no_halo,1) -1)

co2_emits(tdims%i_start:tdims%i_end,                                       &
          tdims%j_start:tdims%j_end)                                       &
    => d1(j_co2_emits:j_co2_emits + field_length(theta_points,no_halo,1) -1)

co2  (tdims_s%i_start:tdims_s%i_end,                                       &
      tdims_s%j_start:tdims_s%j_end,                                       &
      tdims_s%k_start:tdims_s%k_end)                                       &
    => d1(jco2(tdims%k_start):jco2(tdims%k_start) +                        &
                 field_length(theta_points,single_halo,tdims%k_len) -1)

!    level dependent constants
zseak_theta => d1(jzseak_theta : jzseak_theta+(model_levels+1) -1)
Ck_theta    => d1(jck_theta    : jck_theta   +(model_levels+1) -1)
zseak_rho   => d1(jzseak_rho   : jzseak_rho  +(model_levels+1) -1)
Ck_rho      => d1(jck_rho      : jck_rho     +(model_levels+1) -1)
soil_thickness => d1(jsoil_thickness : jsoil_thickness + st_levels)

!    User ancillaries (all 1d, including user_mult)
user_anc1   => d1(juser_anc1  : juser_anc1 + &
 field_length(theta_points,no_halo,1) -1)
user_anc2   => d1(juser_anc2  : juser_anc2 + &
 field_length(theta_points,no_halo,1) -1)
user_anc3   => d1(juser_anc3  : juser_anc3 + &
 field_length(theta_points,no_halo,1) -1)
user_anc4   => d1(juser_anc4  : juser_anc4 + &
 field_length(theta_points,no_halo,1) -1)
user_anc5   => d1(juser_anc5  : juser_anc5 + &
 field_length(theta_points,no_halo,1) -1)
user_anc6   => d1(juser_anc6  : juser_anc6 + &
 field_length(theta_points,no_halo,1) -1)
user_anc7   => d1(juser_anc7  : juser_anc7 + &
 field_length(theta_points,no_halo,1) -1)
user_anc8   => d1(juser_anc8  : juser_anc8 + &
 field_length(theta_points,no_halo,1) -1)
user_anc9   => d1(juser_anc9  : juser_anc9 + &
 field_length(theta_points,no_halo,1) -1)
user_anc10  => d1(juser_anc10 : juser_anc10+ &
 field_length(theta_points,no_halo,1) -1)
user_anc11  => d1(juser_anc11 : juser_anc11+ &
 field_length(theta_points,no_halo,1) -1)
user_anc12  => d1(juser_anc12 : juser_anc12+ &
 field_length(theta_points,no_halo,1) -1)
user_anc13  => d1(juser_anc13 : juser_anc13+ &
 field_length(theta_points,no_halo,1) -1)
user_anc14  => d1(juser_anc14 : juser_anc14+ &
 field_length(theta_points,no_halo,1) -1)
user_anc15  => d1(juser_anc15 : juser_anc15+ &
 field_length(theta_points,no_halo,1) -1)
user_anc16  => d1(juser_anc16 : juser_anc16+ &
 field_length(theta_points,no_halo,1) -1)
user_anc17  => d1(juser_anc17 : juser_anc17+ &
 field_length(theta_points,no_halo,1) -1)
user_anc18  => d1(juser_anc18 : juser_anc18+ &
 field_length(theta_points,no_halo,1) -1)
user_anc19  => d1(juser_anc19 : juser_anc19+ &
 field_length(theta_points,no_halo,1) -1)
user_anc20  => d1(juser_anc20 : juser_anc20+ &
 field_length(theta_points,no_halo,1) -1)
user_mult1  => d1(juser_mult1(1)  : juser_mult1(1) + &
 field_length(theta_points,no_halo,model_levels) -1)
user_mult2  => d1(juser_mult2(1)  : juser_mult2(1) + &
 field_length(theta_points,no_halo,model_levels) -1)
user_mult3  => d1(juser_mult3(1)  : juser_mult3(1) + &
 field_length(theta_points,no_halo,model_levels) -1)
user_mult4  => d1(juser_mult4(1)  : juser_mult4(1) + &
 field_length(theta_points,no_halo,model_levels) -1)
user_mult5  => d1(juser_mult5(1)  : juser_mult5(1) + &
 field_length(theta_points,no_halo,model_levels) -1)
user_mult6  => d1(juser_mult6(1)  : juser_mult6(1) + &
 field_length(theta_points,no_halo,model_levels) -1)
user_mult7  => d1(juser_mult7(1)  : juser_mult7(1) + &
 field_length(theta_points,no_halo,model_levels) -1)
user_mult8  => d1(juser_mult8(1)  : juser_mult8(1) + &
 field_length(theta_points,no_halo,model_levels) -1)
user_mult9  => d1(juser_mult9(1)  : juser_mult9(1) + &
 field_length(theta_points,no_halo,model_levels) -1)
user_mult10 => d1(juser_mult10(1) : juser_mult10(1)+ &
 field_length(theta_points,no_halo,model_levels) -1)
user_mult11 => d1(juser_mult11(1) : juser_mult11(1)+ &
 field_length(theta_points,no_halo,model_levels) -1)
user_mult12 => d1(juser_mult12(1) : juser_mult12(1)+ &
 field_length(theta_points,no_halo,model_levels) -1)
user_mult13 => d1(juser_mult13(1) : juser_mult13(1)+ &
 field_length(theta_points,no_halo,model_levels) -1)
user_mult14 => d1(juser_mult14(1) : juser_mult14(1)+ &
 field_length(theta_points,no_halo,model_levels) -1)
user_mult15 => d1(juser_mult15(1) : juser_mult15(1)+ &
 field_length(theta_points,no_halo,model_levels) -1)
user_mult16 => d1(juser_mult16(1) : juser_mult16(1)+ &
 field_length(theta_points,no_halo,model_levels) -1)
user_mult17 => d1(juser_mult17(1) : juser_mult17(1)+ &
 field_length(theta_points,no_halo,model_levels) -1)
user_mult18 => d1(juser_mult18(1) : juser_mult18(1)+ &
 field_length(theta_points,no_halo,model_levels) -1)
user_mult19 => d1(juser_mult19(1) : juser_mult19(1)+ &
 field_length(theta_points,no_halo,model_levels) -1)
user_mult20 => d1(juser_mult20(1) : juser_mult20(1)+ &
 field_length(theta_points,no_halo,model_levels) -1)

!    Tiled vegetation and triffid
frac_typ(1:land_field,1:ntype) =>d1(jfrac_typ:jfrac_typ  + &
                                  field_length(land_points,no_halo,ntype)-1)

frac_con1=>d1(jfrac_con1:jfrac_con1+ &
                                  field_length(land_points,no_halo,1)-1)
frac_con2=>d1(jfrac_con2:jfrac_con2+ &
                                  field_length(land_points,no_halo,1)-1)
frac_con3=>d1(jfrac_con3:jfrac_con3+ &
                                  field_length(land_points,no_halo,1)-1)
frac_con4=>d1(jfrac_con4:jfrac_con4+ &
                                  field_length(land_points,no_halo,1)-1)
frac_con5=>d1(jfrac_con5:jfrac_con5+ &
                                  field_length(land_points,no_halo,1)-1)
frac_con6=>d1(jfrac_con6:jfrac_con6+ &
                                  field_length(land_points,no_halo,1)-1)
frac_con7=>d1(jfrac_con7:jfrac_con7+ &
                                  field_length(land_points,no_halo,1)-1)
frac_con8=>d1(jfrac_con8:jfrac_con8+ &
                                  field_length(land_points,no_halo,1)-1)
frac_con9=>d1(jfrac_con9:jfrac_con9+ &
                                  field_length(land_points,no_halo,1)-1)

lai_pft(1:land_field,1:npft)  =>d1(jlai_pft:jlai_pft    + &
                                  field_length(land_points,no_halo,npft) -1)
canht_pft(1:land_field,1:npft) =>d1(jcanht_pft:jcanht_pft+ &
                                  field_length(land_points,no_halo,npft)-1)
disturb_veg =>  d1(jdisturb:jdisturb+ &
                                  field_length(land_points,no_halo,1)-1)
disturb_veg_prev =>  d1(jdisturb_prev:jdisturb_prev+ &  
                                  field_length(land_points,no_halo,1)-1)
pasture_frac_d1 =>  d1(jpasture:jpasture+ &
                                  field_length(land_points,no_halo,1)-1)
pasture_frac_prev_d1 =>  d1(jpasture_prev:jpasture_prev+ &  
                                  field_length(land_points,no_halo,1)-1)
agr_crop_frac_d1 =>  d1(jagr_crop:jagr_crop+ &
                                  field_length(land_points,no_halo,1)-1)
agr_crop_frac_prev_d1 =>  d1(jagr_crop_prev:jagr_crop_prev+ &  
                                  field_length(land_points,no_halo,1)-1)
wood_prod_fast_d1 =>  d1(jwoodprod_fast:jwoodprod_fast+ &
                                  field_length(land_points,no_halo,1)-1)
wood_prod_med_d1 =>  d1(jwoodprod_med:jwoodprod_med+ &   
                                  field_length(land_points,no_halo,1)-1)
wood_prod_slow_d1 =>  d1(jwoodprod_slow:jwoodprod_slow+ &
                                  field_length(land_points,no_halo,1)-1)

soil_alb  =>  d1(jsoil_alb:jsoil_alb+ &
                                  field_length(land_points,no_halo,1)-1)
obs_alb_sw  =>  d1(jobs_alb_sw:jobs_alb_sw+ &
                                  field_length(land_points,no_halo,1)-1)
obs_alb_vis  =>  d1(jobs_alb_vis:jobs_alb_vis+ &
                                  field_length(land_points,no_halo,1)-1)
obs_alb_nir  =>  d1(jobs_alb_nir:jobs_alb_nir+ &
                                  field_length(land_points,no_halo,1)-1)

soil_carb(1:land_field,1:dim_cs1) =>d1(jsoil_carb:jsoil_carb + &
                                  field_length(land_points,no_halo,dim_cs1)-1)

soil_carb1(1:land_field,1:dim_cs1) => d1(jsoil_carb1:jsoil_carb1+ &
                                  field_length(land_points,no_halo,dim_cs1) -1)

soil_carb2(1:land_field,1:dim_cs1) => d1(jsoil_carb2:jsoil_carb2+ &
                                  field_length(land_points,no_halo,1) -1)
soil_carb3(1:land_field,1:dim_cs1) => d1(jsoil_carb3:jsoil_carb3+ &
                                  field_length(land_points,no_halo,1) -1)
soil_carb4(1:land_field,1:dim_cs1) => d1(jsoil_carb4:jsoil_carb4+ &
                                  field_length(land_points,no_halo,1) -1)

soil_nitro1(1:land_field) => d1(jsoil_nitro1:jsoil_nitro1+ &
                                  field_length(land_points,no_halo,1) -1)
soil_nitro2(1:land_field) => d1(jsoil_nitro2:jsoil_nitro2+ &
                                  field_length(land_points,no_halo,1) -1)
soil_nitro3(1:land_field) => d1(jsoil_nitro3:jsoil_nitro3+ &
                                  field_length(land_points,no_halo,1) -1)
soil_nitro4(1:land_field) => d1(jsoil_nitro4:jsoil_nitro4+ &
                                  field_length(land_points,no_halo,1) -1)
soil_inorgnit(1:land_field) => d1(jsoil_inorgnit:jsoil_inorgnit+ &
                                  field_length(land_points,no_halo,1) -1)
nitrogen_deposition_d1(1:land_field) => d1(j_n_deposition:j_n_deposition+ &
                                  field_length(land_points,no_halo,1) -1)

g_lf_pft_acc(1:land_field,1:npft)   => d1(jg_lf_pft_acc:jg_lf_pft_acc + &
                               field_length(land_points,no_halo,npft) -1)

g_phlf_pft_acc(1:land_field,1:npft) => d1(jg_phlf_pft_acc:jg_phlf_pft_acc + &
                               field_length(land_points,no_halo,npft) -1)

IF (l_triffid) THEN

  npp_pft_acc(1:land_pts_trif,1:npft_trif)   =>   &
                   d1(jnpp_pft_acc:jnpp_pft_acc + &
                               field_length(land_points,no_halo,npft_trif)-1)

  rsp_w_pft_acc(1:land_pts_trif,1:npft_trif) =>       &
                   d1(jrsp_w_pft_acc:jrsp_w_pft_acc + &
                               field_length(land_points,no_halo,npft_trif)-1)

ELSE

  npp_pft_acc  => dummy_fld2

  rsp_w_pft_acc => dummy_fld2

END IF

rsp_s_acc(1:land_field,1:dim_cs1)  => d1(jrsp_s_acc:jrsp_s_acc + &
                                  field_length(land_points,no_halo,dim_cs1)-1)

rsp_s_acc1(1:land_field,1:dim_cs1) => d1(jrsp_s_acc1:jrsp_s_acc1+ &
                                  field_length(land_points,no_halo,dim_cs1)-1)

rsp_s_acc2(1:land_field,1:dim_cs1) => d1(jrsp_s_acc2:jrsp_s_acc2+ &
                                  field_length(land_points,no_halo,1) -1)

rsp_s_acc3(1:land_field,1:dim_cs1) => d1(jrsp_s_acc3:jrsp_s_acc3+ &
                                  field_length(land_points,no_halo,1) -1)

rsp_s_acc4(1:land_field,1:dim_cs1) => d1(jrsp_s_acc4:jrsp_s_acc4+ &
                                  field_length(land_points,no_halo,1) -1)

can_water_tile(1:land_field,1:ntiles)                                     &
        => d1(jcan_water_tile:jcan_water_tile +                            &
                                  field_length(land_points,no_halo,ntiles) -1)

catch_tile(1:land_field,1:ntiles)                                         &
        => d1(jcatch_tile:jcatch_tile +                                    &
                                  field_length(land_points,no_halo,ntiles) -1)
rgrain_tile(1:land_field,1:ntiles)                                        &
        => d1(jrgrain_tile:jrgrain_tile +                                  &
                                  field_length(land_points,no_halo,ntiles) -1)

tstar_tile(1:land_field,1:ntiles)                                         &
        => d1(jtstar_tile:jtstar_tile +                                    &
                                  field_length(land_points,no_halo,ntiles) -1)
tsurf_elev_surft(1:land_field,1:ntiles)                                   &
        => d1(jtsurf_elev_surft:jtsurf_elev_surft +                        &
                                  field_length(land_points,no_halo,ntiles) -1)

z0_tile(1:land_field,1:ntiles)                                            &
        => d1(jz0_tile:jz0_tile+field_length(land_points,no_halo,ntiles) -1)

z0h_tile(1:land_field,1:ntiles)                                           &
        => d1(jz0h_tile:jz0h_tile+field_length(land_points,no_halo,ntiles) -1)

snodep_tile(1:land_field,1:ntiles)                                        &
        => d1(jsnodep_tile:jsnodep_tile +                                  &
                                  field_length(land_points,no_halo,ntiles) -1)

infil_tile(1:land_field,1:ntiles)                                         &
        => d1(jinfil_tile:jinfil_tile +                                    &
                                  field_length(land_points,no_halo,ntiles) -1)

dolr_field(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)           &
        => d1(jdolr:jdolr + field_length(theta_points,no_halo,1) -1)

lw_down(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)              &
        => d1(jlw_down:jlw_down + field_length(theta_points,no_halo,1) -1)

sw_tile_rts(1:land_field,1:ntiles)                                        &
        => d1(jsw_tile:jsw_tile+field_length(land_points,no_halo,ntiles) -1)

! MORUSES - new urban two-tile scheme
hgt   => d1(jurbhgt  :jurbhgt   +field_length(land_points,no_halo,1) -1)
! building height
hwr   => d1(jurbhwr  :jurbhwr   +field_length(land_points,no_halo,1) -1)
! height to width
wrr   => d1(jurbwrr  :jurbwrr   +field_length(land_points,no_halo,1) -1)
! width ratio
disp  => d1(jurbdisp :jurbdisp  +field_length(land_points,no_halo,1) -1)
! displacement height
ztm   => d1(jurbztm  :jurbztm   +field_length(land_points,no_halo,1) -1)
!
albwl => d1(jurbalbwl:jurbalbwl +field_length(land_points,no_halo,1) -1)
! wall albedo
albrd => d1(jurbalbrd:jurbalbrd +field_length(land_points,no_halo,1) -1)
! road albedo
emisw => d1(jurbemisw:jurbemisw +field_length(land_points,no_halo,1) -1)
! wall emissivity
emisr => d1(jurbemisr:jurbemisr +field_length(land_points,no_halo,1) -1)
! road emissivity

!    River routing fields
riv_sequence  => d1(jriv_sequence : jriv_sequence+ &
                                field_length(river_points,no_halo,1) -1)
riv_direction => d1(jriv_direction : jriv_direction+ &
                                field_length(river_points,no_halo,1) -1)
riv_storage   => d1(jriv_storage : jriv_storage+ &
                                field_length(river_points,no_halo,1) -1)
tot_surfroff  => d1(jtot_surfroff : jtot_surfroff+ &
                                field_length(land_points,no_halo,1) -1)
tot_subroff   => d1(jtot_subroff : jtot_subroff+ &
                                field_length(land_points,no_halo,1) -1)
riv_inlandatm => d1(jriv_inlandatm : jriv_inlandatm+ &
                                field_length(land_points,no_halo,1)  -1)
! these are uninitialised upon entering ATM_STEP
riv_iarea     => dummy_field !D1(1:1+row_length*rows)
riv_slope     => dummy_field !D1(1:1+row_length*rows)
riv_flowobs1  => dummy_field !D1(1:1+row_length*rows)
riv_inext     => dummy_field !D1(1:1+row_length*rows)
riv_jnext     => dummy_field !D1(1:1+row_length*rows)
riv_land      => dummy_field !D1(1:1+row_length*rows)
riv_substore  => dummy_field !D1(1:1+row_length*rows)
riv_surfstore => dummy_field !D1(1:1+row_length*rows)
riv_flowin    => dummy_field !D1(1:1+row_length*rows)
riv_bflowin   => dummy_field !D1(1:1+row_length*rows)

!    Required for water conservation correction due to lake evaporation
acc_lake_evap => d1(jacc_lake_evap:jacc_lake_evap                  &
              +field_length(theta_points,no_halo,1) -1)

!    Fields to be retained in dumps for coupled models using OASIS
c_solar(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end) &
   => d1(jc_solar : jc_solar + &
               field_length(theta_points,no_halo,1) -1)

c_blue(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end) &
   =>  d1(jc_blue : jc_blue + &
               field_length(theta_points,no_halo,1) -1)

c_longwave(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end) &
   => d1(jc_longwave : jc_longwave + &
               field_length(theta_points,no_halo,1) -1)

c_taux(1:row_length, udims%j_start:udims%j_end) &
   => d1(jc_taux : jc_taux + &
               field_length(u_points,no_halo,1) -1)

c_tauy(vdims%i_start:vdims%i_end, 1:n_rows) &
   => d1(jc_tauy : jc_tauy + &
               field_length(v_points,no_halo,1) -1)

c_w10(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end) &
   => d1(jc_w10 : jc_w10 + &
               field_length(theta_points,no_halo,1) -1)

c_sensible(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end) &
   => d1(jc_sensible : jc_sensible + &
               field_length(theta_points,no_halo,1) -1)

c_sublim(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, 1:nice_use) &
   =>  d1(jc_sublim : jc_sublim + &
               field_length(theta_points,no_halo,nice_use) -1)

c_evap(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end) &
   =>  d1(jc_evap : jc_evap + &
               field_length(theta_points,no_halo,1) -1)

c_fcondtopn(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, 1:nice) &
   => d1(jc_fcondtopn : jc_fcondtopn + &
               field_length(theta_points,no_halo,nice) -1)

c_topmeltn(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, 1:nice) &
   => d1(jc_topmeltn : jc_topmeltn + &
               field_length(theta_points,no_halo,nice) -1)

c_tstar_sicen(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, 1:nice) &
   => d1(jc_tstar_sicen:jc_tstar_sicen + &
                    field_length(theta_points,no_halo,nice) -1)

c_lsrain(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end) &
   =>  d1(jc_lsrain : jc_lsrain + &
               field_length(theta_points,no_halo,1) -1)

c_lssnow(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end) &
   =>  d1(jc_lssnow : jc_lssnow + &
               field_length(theta_points,no_halo,1) -1)

c_cvrain(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end) &
   =>  d1(jc_cvrain : jc_cvrain + &
               field_length(theta_points,no_halo,1) -1)

c_cvsnow(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end) &
   =>  d1(jc_cvsnow : jc_cvsnow + &
               field_length(theta_points,no_halo,1) -1)

c_riverout(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end) &
   => d1(jc_riverout : jc_riverout + &
               field_length(theta_points,no_halo,1) -1)

c_calving(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end) &
   => d1(jc_calving : jc_calving + &
               field_length(theta_points,no_halo,1) -1)

c_mslp(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end) &
   => d1(jc_mslp : jc_mslp + &
               field_length(theta_points,no_halo,1) -1)

c_surf_CO2(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end) &
   => d1(jc_surf_CO2 : jc_surf_CO2 + &
               field_length(theta_points,no_halo,1) -1)

c_dust_dep(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end) &
   => d1(jc_dust_dep : jc_dust_dep + &
               field_length(theta_points,no_halo,1) -1)

! JULES 2 prognostics
snowdepth(1:land_field, 1:ntiles)                                           &
  =>  d1(jsnowdepth :jsnowdepth + field_length(land_points,no_halo,ntiles) -1)

rho_snow_grnd(1:land_field, 1:ntiles)                                       &
  =>  d1(jrho_snow_grnd :jrho_snow_grnd +                                    &
                                  field_length(land_points,no_halo,ntiles) -1)

nsnow (1:land_field, 1:ntiles)                                              &
  =>  d1(jnsnow :jnsnow +         field_length(land_points,no_halo,ntiles) -1)

ds (1:land_field, 1:ntiles, 1:nsmax)                                        &
  =>  d1(jds :jds         + field_length(land_points,no_halo,ntiles*nsmax) -1)

sice (1:land_field, 1:ntiles, 1:nsmax)                                      &
  =>  d1(jsice :jsice     + field_length(land_points,no_halo,ntiles*nsmax) -1)

sliq (1:land_field, 1:ntiles, 1:nsmax)                                      &
  =>  d1(jsliq :jsliq     + field_length(land_points,no_halo,ntiles*nsmax) -1)

tsnowlayer (1:land_field, 1:ntiles, 1:nsmax)                                &
  =>  d1(jtsnowlayer :jtsnowlayer  +                                         &
                            field_length(land_points,no_halo,ntiles*nsmax) -1)

rho_snow (1:land_field, 1:ntiles, 1:nsmax)                                  &
  =>  d1(jrho_snow :jrho_snow +                                              &
                            field_length(land_points,no_halo,ntiles*nsmax) -1)

rgrainl  (1:land_field, 1:ntiles, 1:nsmax)                                  &
  =>  d1(jrgrainl:jrgrainl +field_length(land_points,no_halo,ntiles*nsmax) -1)

! FLake lake scheme prognostics
lake_depth     =>  d1(jlake_depth     :jlake_depth     + &
  field_length(land_points,no_halo,1) -1)
lake_fetch     =>  d1(jlake_fetch     :jlake_fetch     + &
  field_length(land_points,no_halo,1) -1)
lake_t_mean    =>  d1(jlake_t_mean    :jlake_t_mean    + &
  field_length(land_points,no_halo,1) -1)
lake_t_mxl     =>  d1(jlake_t_mxl     :jlake_t_mxl     + &
  field_length(land_points,no_halo,1) -1)
lake_t_ice     =>  d1(jlake_t_ice     :jlake_t_ice     + &
  field_length(land_points,no_halo,1) -1)
lake_h_mxl     =>  d1(jlake_h_mxl     :jlake_h_mxl     + &
  field_length(land_points,no_halo,1) -1)
lake_h_ice     =>  d1(jlake_h_ice     :jlake_h_ice     + &
  field_length(land_points,no_halo,1) -1)
lake_shape     =>  d1(jlake_shape     :jlake_shape     + &
  field_length(land_points,no_halo,1))
lake_g_dt      =>  d1(jlake_g_dt      :jlake_g_dt      + &
  field_length(land_points,no_halo,1))

!    Required for energy correction
net_flux(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)               &
       => d1(jnet_flux:jnet_flux  + field_length(theta_points,no_halo,1) -1)
net_mflux(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)              &
       => d1(jnet_mflux:jnet_mflux + field_length(theta_points,no_halo,1) -1)

!    Fields carried forward from previous version
tstar_anom => d1(jtstar_anom : jtstar_anom + &
                                field_length(theta_points,no_halo,1) -1)

!    lateral boundary conditions

orog_lbc  => d1(jorog_lbc : jorog_lbc + &
  lenrima(fld_type_p,halo_type_extended,1)*1  -1)

u_lbc(1:lenrima(fld_type_u,halo_type_extended,1),udims_s%k_start:udims_s%k_end)&
     => d1(ju_lbc : ju_lbc + lenrima(fld_type_u,halo_type_extended,1)* &
                                  (udims_l%k_len) -1)

v_lbc(1:lenrima(fld_type_v,halo_type_extended,1),vdims_s%k_start:vdims_s%k_end)&
     => d1(jv_lbc : jv_lbc +  lenrima(fld_type_v,halo_type_extended,1)* &
                                  (vdims_l%k_len) -1)

! model_levels+1
w_lbc(1:lenrima(fld_type_p,halo_type_extended,1),wdims_s%k_start:wdims_s%k_end)&
 => d1(jw_lbc :jw_lbc +lenrima(fld_type_p,halo_type_extended,rima_type_norm)* &
                                  (wdims_l%k_len) -1)

rho_lbc(1:lenrima(fld_type_p,halo_type_extended,1), &
                                               pdims_s%k_start:pdims_s%k_end) &
 => d1(jrho_lbc : jrho_lbc + lenrima(fld_type_p,halo_type_extended,1)* &
                                  (pdims_l%k_len) -1)

theta_lbc(1:lenrima(fld_type_p,halo_type_extended,1), &
                                            tdims_s%k_start:tdims_s%k_end) &
 => d1(jtheta_lbc : jtheta_lbc + lenrima(fld_type_p,halo_type_extended,1)* &
                                  (tdims_l%k_len) -1)

q_lbc (1:lenrima(fld_type_p,halo_type_extended,1), &
                                                tdims_l%k_start:tdims_l%k_end) &
 => d1(jq_lbc : jq_lbc + lenrima(fld_type_p,halo_type_extended,1)* &
                                  (tdims_l%k_len) -1)

qcl_lbc(1:lenrima(fld_type_p,halo_type_extended,1), &
                                                 tdims_l%k_start:tdims_l%k_end)&
 => d1(jqcl_lbc : jqcl_lbc + lenrima(fld_type_p,halo_type_extended,1)* &
                                  (tdims_l%k_len) -1)

qcf_lbc(1:lenrima(fld_type_p,halo_type_extended,1), &
                                                 tdims_l%k_start:tdims_l%k_end)&
  => d1(jqcf_lbc : jqcf_lbc + lenrima(fld_type_p,halo_type_extended,1)* &
                                  (tdims_l%k_len) -1)

! model_levels+1
exner_lbc(1:lenrima(fld_type_p,halo_type_extended,1), &
                                           pdims_s%k_start:pdims_s%k_end+1) &
  => d1(jexner_lbc : jexner_lbc + lenrima(fld_type_p,halo_type_extended,1)* &
                                  (pdims_l%k_len+1) -1)

u_adv_lbc(1:lenrima(fld_type_u,halo_type_extended,1), &
                                             udims_s%k_start:udims_s%k_end) &
  => d1(ju_adv_lbc : ju_adv_lbc + lenrima(fld_type_u,halo_type_extended,1)* &
                                   (udims_l%k_len) -1)

v_adv_lbc(1:lenrima(fld_type_v,halo_type_extended,1), &
                                             vdims_s%k_start:vdims_s%k_end) &
  => d1(jv_adv_lbc : jv_adv_lbc + lenrima(fld_type_v,halo_type_extended,1)* &
                                   (vdims_l%k_len) -1)

! model_levels+1
w_adv_lbc(1:lenrima(fld_type_p,halo_type_extended,1), &
                                             wdims_s%k_start:wdims_s%k_end) &
  => d1(jw_adv_lbc : jw_adv_lbc + lenrima(fld_type_p,halo_type_extended,1)* &
                                   (wdims_l%k_len) -1)

qcf2_lbc(1:lenrima(fld_type_p,halo_type_extended,1), &
                                             tdims_l%k_start:tdims_l%k_end) &
  => d1(jqcf2_lbc : jqcf2_lbc + lenrima(fld_type_p,halo_type_extended,1)*   &
                                  (tdims_l%k_len) -1)

qrain_lbc(1:lenrima(fld_type_p,halo_type_extended,1), &
                                             tdims_l%k_start:tdims_l%k_end) &
  => d1(jqrain_lbc : jqrain_lbc + lenrima(fld_type_p,halo_type_extended,1)* &
                                  (tdims_l%k_len) -1)

qgraup_lbc(1:lenrima(fld_type_p,halo_type_extended,1), &
                                             tdims_l%k_start:tdims_l%k_end)  &
   => d1(jqgraup_lbc :jqgraup_lbc +lenrima(fld_type_p,halo_type_extended,1)* &
                                  (tdims_l%k_len) -1)

cf_bulk_lbc(1:lenrima(fld_type_p,halo_type_extended,1), &
                                             tdims_l%k_start:tdims_l%k_end)    &
   => d1(jcf_bulk_lbc :jcf_bulk_lbc +lenrima(fld_type_p,halo_type_extended,1)* &
                                  (tdims_l%k_len) -1)

cf_liquid_lbc(1:lenrima(fld_type_p,halo_type_extended,1), &
                                             tdims_l%k_start:tdims_l%k_end) &
   => d1(jcf_liquid_lbc : jcf_liquid_lbc +                                  &
                                  lenrima(fld_type_p,halo_type_extended,1)* &
                                  (tdims_l%k_len) -1)

cf_frozen_lbc(1:lenrima(fld_type_p,halo_type_extended,1), &
                                             tdims_l%k_start:tdims_l%k_end) &
   => d1(jcf_frozen_lbc : jcf_frozen_lbc +                                  &
                                  lenrima(fld_type_p,halo_type_extended,1)* &
                                  (tdims_l%k_len) -1)

murk_lbc(1:lenrima(fld_type_p,halo_type_single,1), &
                                                tdims_s%k_start:tdims_s%k_end) &
 => d1(jmurk_lbc : jmurk_lbc + lenrima(fld_type_p,halo_type_single,1)* &
                                  (tdims_s%k_len) -1)

dust_div1_lbc(1:lenrima(fld_type_p,halo_type_single,1), &
                                                tdims_s%k_start:tdims_s%k_end) &
 => d1(jdust_div1_lbc :jdust_div1_lbc +lenrima(fld_type_p,halo_type_single,1)* &
                                  (tdims_s%k_len) -1)

dust_div2_lbc(1:lenrima(fld_type_p,halo_type_single,1), &
                                                tdims_s%k_start:tdims_s%k_end) &
 => d1(jdust_div2_lbc :jdust_div2_lbc +lenrima(fld_type_p,halo_type_single,1)* &
                                  (tdims_s%k_len) -1)

dust_div3_lbc(1:lenrima(fld_type_p,halo_type_single,1), &
                                                tdims_s%k_start:tdims_s%k_end) &
 => d1(jdust_div3_lbc :jdust_div3_lbc +lenrima(fld_type_p,halo_type_single,1)* &
                                  (tdims_s%k_len) -1)

dust_div4_lbc(1:lenrima(fld_type_p,halo_type_single,1), &
                                              tdims_s%k_start:tdims_s%k_end) &
 => d1(jdust_div4_lbc :jdust_div4_lbc + &
  lenrima(fld_type_p,halo_type_single,1)* &
                                  (tdims_s%k_len) -1)

dust_div5_lbc(1:lenrima(fld_type_p,halo_type_single,1), &
                                              tdims_s%k_start:tdims_s%k_end) &
 => d1(jdust_div5_lbc :jdust_div5_lbc +lenrima(fld_type_p,halo_type_single,1)* &
                                  (tdims_s%k_len) -1)

dust_div6_lbc(1:lenrima(fld_type_p,halo_type_single,1), &
                                              tdims_s%k_start:tdims_s%k_end) &
 => d1(jdust_div6_lbc :jdust_div6_lbc +lenrima(fld_type_p,halo_type_single,1)* &
                                  (tdims_s%k_len) -1)

so2_lbc(1:lenrima(fld_type_p,halo_type_single,1), &
                                              tdims_s%k_start:tdims_s%k_end) &
 => d1(jso2_lbc :jso2_lbc + lenrima(fld_type_p,halo_type_single,1)* &
                                  (tdims_s%k_len) -1)

dms_lbc(1:lenrima(fld_type_p,halo_type_single,1), &
                                              tdims_s%k_start:tdims_s%k_end) &
 => d1(jdms_lbc :jdms_lbc + lenrima(fld_type_p,halo_type_single,1)* &
                                  (tdims_s%k_len) -1)

so4_aitken_lbc(1:lenrima(fld_type_p,halo_type_single,1), &
                                              tdims_s%k_start:tdims_s%k_end) &
 => d1(jso4_aitken_lbc :jso4_aitken_lbc + &
                                     lenrima(fld_type_p,halo_type_single,1)* &
                                  (tdims_s%k_len) -1)

so4_accu_lbc(1:lenrima(fld_type_p,halo_type_single,1), &
                                               tdims_s%k_start:tdims_s%k_end) &
 => d1(jso4_accu_lbc :jso4_accu_lbc + lenrima(fld_type_p,halo_type_single,1)* &
                                  (tdims_s%k_len) -1)

so4_diss_lbc(1:lenrima(fld_type_p,halo_type_single,1), &
                                              tdims_s%k_start:tdims_s%k_end) &
 => d1(jso4_diss_lbc :jso4_diss_lbc + lenrima(fld_type_p,halo_type_single,1)* &
                                  (tdims_s%k_len) -1)

nh3_lbc(1:lenrima(fld_type_p,halo_type_single,1), &
                                              tdims_s%k_start:tdims_s%k_end) &
 => d1(jnh3_lbc :jnh3_lbc + lenrima(fld_type_p,halo_type_single,1)* &
                                  (tdims_s%k_len) -1)

soot_new_lbc(1:lenrima(fld_type_p,halo_type_single,1), &
                                              tdims_s%k_start:tdims_s%k_end) &
 => d1(jsoot_new_lbc :jsoot_new_lbc + lenrima(fld_type_p,halo_type_single,1)* &
                                  (tdims_s%k_len) -1)

soot_agd_lbc(1:lenrima(fld_type_p,halo_type_single,1), &
                                              tdims_s%k_start:tdims_s%k_end) &
 => d1(jsoot_agd_lbc :jsoot_agd_lbc + lenrima(fld_type_p,halo_type_single,1)* &
                                  (tdims_s%k_len) -1)

soot_cld_lbc(1:lenrima(fld_type_p,halo_type_single,1), &
                                              tdims_s%k_start:tdims_s%k_end) &
 => d1(jsoot_cld_lbc :jsoot_cld_lbc + lenrima(fld_type_p,halo_type_single,1)* &
                                  (tdims_s%k_len) -1)

bmass_new_lbc(1:lenrima(fld_type_p,halo_type_single,1), &
                                              tdims_s%k_start:tdims_s%k_end) &
 => d1(jbmass_new_lbc :jbmass_new_lbc +lenrima(fld_type_p,halo_type_single,1)* &
                                  (tdims_s%k_len) -1)

bmass_agd_lbc(1:lenrima(fld_type_p,halo_type_single,1), &
                                              tdims_s%k_start:tdims_s%k_end) &
 => d1(jbmass_agd_lbc :jbmass_agd_lbc +lenrima(fld_type_p,halo_type_single,1)* &
                                  (tdims_s%k_len) -1)

bmass_cld_lbc(1:lenrima(fld_type_p,halo_type_single,1), &
                                              tdims_s%k_start:tdims_s%k_end) &
 => d1(jbmass_cld_lbc :jbmass_cld_lbc +lenrima(fld_type_p,halo_type_single,1)* &
                                  (tdims_s%k_len) -1)

ocff_new_lbc(1:lenrima(fld_type_p,halo_type_single,1), &
                                              tdims_s%k_start:tdims_s%k_end) &
 => d1(jocff_new_lbc :jocff_new_lbc +lenrima(fld_type_p,halo_type_single,1)* &
                                  (tdims_s%k_len) -1)

ocff_agd_lbc(1:lenrima(fld_type_p,halo_type_single,1), &
                                              tdims_s%k_start:tdims_s%k_end) &
 => d1(jocff_agd_lbc :jocff_agd_lbc +lenrima(fld_type_p,halo_type_single,1)* &
                                  (tdims_s%k_len) -1)

ocff_cld_lbc(1:lenrima(fld_type_p,halo_type_single,1), &
                                              tdims_s%k_start:tdims_s%k_end) &
 => d1(jocff_cld_lbc :jocff_cld_lbc +lenrima(fld_type_p,halo_type_single,1)* &
                                  (tdims_s%k_len) -1)

nitr_acc_lbc(1:lenrima(fld_type_p,halo_type_single,1), &
                                              tdims_s%k_start:tdims_s%k_end) &
 => d1(jnitr_acc_lbc :jnitr_acc_lbc +lenrima(fld_type_p,halo_type_single,1)* &
                                  (tdims_s%k_len) -1)

nitr_diss_lbc(1:lenrima(fld_type_p,halo_type_single,1), &
                                                tdims_s%k_start:tdims_s%k_end) &
 => d1(jnitr_diss_lbc :jnitr_diss_lbc +lenrima(fld_type_p,halo_type_single,1)* &
                                  (tdims_s%k_len) -1)

u_lbc_tend(1:lenrima(fld_type_u,halo_type_extended,1), &
                                               udims_s%k_start:udims_s%k_end)&
 => d1(ju_lbc_tend : ju_lbc_tend + lenrima(fld_type_u,halo_type_extended,1)* &
                                  (udims_l%k_len) -1)

v_lbc_tend(1:lenrima(fld_type_v,halo_type_extended,1), &
                                               vdims_s%k_start:vdims_s%k_end)&
 => d1(jv_lbc_tend : jv_lbc_tend + lenrima(fld_type_v,halo_type_extended,1)* &
                                  (vdims_l%k_len) -1)

! model_levels+1
w_lbc_tend(1:lenrima(fld_type_p,halo_type_extended,1),&
                                               wdims_s%k_start:wdims_s%k_end)&
 => d1(jw_lbc_tend : jw_lbc_tend + &
                      lenrima(fld_type_p,halo_type_extended,rima_type_norm)* &
                                  (wdims_l%k_len) -1)

rho_lbc_tend(1:lenrima(fld_type_p,halo_type_extended,1), &
                                              pdims_s%k_start:pdims_s%k_end) &
 => d1(jrho_lbc_tend : jrho_lbc_tend + &
  lenrima(fld_type_p,halo_type_extended,1)* &
                                  (pdims_l%k_len) -1)

theta_lbc_tend(1:lenrima(fld_type_p,halo_type_extended,1), &
                                              tdims_s%k_start:tdims_s%k_end) &
 => d1(jtheta_lbc_tend : jtheta_lbc_tend + &
                                   lenrima(fld_type_p,halo_type_extended,1)* &
                                  (tdims_l%k_len) -1)

q_lbc_tend(1:lenrima(fld_type_p,halo_type_extended,1), &
                                              tdims_l%k_start:tdims_l%k_end) &
  => d1(jq_lbc_tend : jq_lbc_tend + &
                                   lenrima(fld_type_p,halo_type_extended,1)* &
                                  (tdims_l%k_len) -1)

qcl_lbc_tend(1:lenrima(fld_type_p,halo_type_extended,1), &
                                               tdims_l%k_start:tdims_l%k_end)&
  => d1(jqcl_lbc_tend : jqcl_lbc_tend + &
                                   lenrima(fld_type_p,halo_type_extended,1)* &
                                  (tdims_l%k_len) -1)

qcf_lbc_tend(1:lenrima(fld_type_p,halo_type_extended,1), &
                                               tdims_l%k_start:tdims_l%k_end)&
  => d1(jqcf_lbc_tend : jqcf_lbc_tend + &
                                   lenrima(fld_type_p,halo_type_extended,1)* &
                                  (tdims_l%k_len) -1)

! model_levels+1
exner_lbc_tend(1:lenrima(fld_type_p,halo_type_extended,1), &
                                            pdims_s%k_start:pdims_s%k_end+1) &
 => d1(jexner_lbc_tend : jexner_lbc_tend + &
  lenrima(fld_type_p,halo_type_extended,1)*(pdims_l%k_len+1) -1)

u_adv_lbc_tend(1:lenrima(fld_type_u,halo_type_extended,1), &
                                             udims_s%k_start:udims_s%k_end) &
 => d1(ju_adv_lbc_tend : ju_adv_lbc_tend + &
                                   lenrima(fld_type_u,halo_type_extended,1)* &
                                   (udims_l%k_len) -1)

v_adv_lbc_tend(1:lenrima(fld_type_u,halo_type_extended,1), &
                                              udims_s%k_start:udims_s%k_end) &
 => d1(jv_adv_lbc_tend : jv_adv_lbc_tend + &
                                   lenrima(fld_type_v,halo_type_extended,1)* &
                                   (vdims_l%k_len) -1)

! model_levels+1
w_adv_lbc_tend(1:lenrima(fld_type_p,halo_type_extended,1), &
                                             wdims_s%k_start:wdims_s%k_end) &
 => d1(jw_adv_lbc_tend : jw_adv_lbc_tend + &
                                  lenrima(fld_type_p,halo_type_extended,1)* &
                                   (wdims_l%k_len) -1)

qcf2_lbc_tend(1:lenrima(fld_type_p,halo_type_extended,1), &
                                             tdims_l%k_start:tdims_l%k_end) &
 => d1(jqcf2_lbc_tend : jqcf2_lbc_tend + &
                                  lenrima(fld_type_p,halo_type_extended,1)* &
                                   (tdims_l%k_len) -1)

qrain_lbc_tend(1:lenrima(fld_type_p,halo_type_extended,1), &
                                             tdims_l%k_start:tdims_l%k_end) &
 => d1(jqrain_lbc_tend : jqrain_lbc_tend + &
                                  lenrima(fld_type_p,halo_type_extended,1)* &
                                   (tdims_l%k_len) -1)

qgraup_lbc_tend(1:lenrima(fld_type_p,halo_type_extended,1), &
                                            tdims_l%k_start:tdims_l%k_end)  &
 => d1(jqgraup_lbc_tend : jqgraup_lbc_tend + &
                                  lenrima(fld_type_p,halo_type_extended,1)* &
                                   (tdims_l%k_len) -1)

cf_bulk_lbc_tend(1:lenrima(fld_type_p,halo_type_extended,1), &
                                             tdims_l%k_start:tdims_l%k_end) &
 => d1(jcf_bulk_lbc_tend : jcf_bulk_lbc_tend + &
                                  lenrima(fld_type_p,halo_type_extended,1)* &
                                   (tdims_l%k_len) -1)

cf_liquid_lbc_tend(1:lenrima(fld_type_p,halo_type_extended,1), &
                                             tdims_l%k_start:tdims_l%k_end) &
 => d1(jcf_liquid_lbc_tend : jcf_liquid_lbc_tend + &
                                  lenrima(fld_type_p,halo_type_extended,1)* &
                                   (tdims_l%k_len) -1)

cf_frozen_lbc_tend(1:lenrima(fld_type_p,halo_type_extended,1), &
                                             tdims_l%k_start:tdims_l%k_end) &
 => d1(jcf_frozen_lbc_tend : jcf_frozen_lbc_tend + &
                                  lenrima(fld_type_p,halo_type_extended,1)* &
                                   (tdims_l%k_len) -1)

murk_lbc_tend(1:lenrima(fld_type_p,halo_type_single,1), &
                                             tdims_s%k_start:tdims_s%k_end) &
 => d1(jmurk_lbc_tend : jmurk_lbc_tend + &
                                    lenrima(fld_type_p,halo_type_single,1)* &
                                   (tdims_s%k_len) -1)

dust_div1_lbc_tend(1:lenrima(fld_type_p,halo_type_single,1), &
                                             tdims_s%k_start:tdims_s%k_end) &
 => d1(jdust_div1_lbc_tend : jdust_div1_lbc_tend + &
                                    lenrima(fld_type_p,halo_type_single,1)* &
                                   (tdims_s%k_len) -1)

dust_div2_lbc_tend(1:lenrima(fld_type_p,halo_type_single,1), &
                                             tdims_s%k_start:tdims_s%k_end) &
 => d1(jdust_div2_lbc_tend : jdust_div2_lbc_tend + &
                                    lenrima(fld_type_p,halo_type_single,1)* &
                                   (tdims_s%k_len) -1)

dust_div3_lbc_tend(1:lenrima(fld_type_p,halo_type_single,1), &
                                             tdims_s%k_start:tdims_s%k_end) &
 => d1(jdust_div3_lbc_tend : jdust_div3_lbc_tend + &
                                    lenrima(fld_type_p,halo_type_single,1)* &
                                   (tdims_s%k_len) -1)

dust_div4_lbc_tend(1:lenrima(fld_type_p,halo_type_single,1), &
                                             tdims_s%k_start:tdims_s%k_end) &
 => d1(jdust_div4_lbc_tend : jdust_div4_lbc_tend + &
                                    lenrima(fld_type_p,halo_type_single,1)* &
                                   (tdims_s%k_len) -1)

dust_div5_lbc_tend(1:lenrima(fld_type_p,halo_type_single,1), &
                                             tdims_s%k_start:tdims_s%k_end) &
 => d1(jdust_div5_lbc_tend : jdust_div5_lbc_tend + &
                                    lenrima(fld_type_p,halo_type_single,1)* &
                                   (tdims_s%k_len) -1)

dust_div6_lbc_tend(1:lenrima(fld_type_p,halo_type_single,1), &
                                             tdims_s%k_start:tdims_s%k_end) &
 => d1(jdust_div6_lbc_tend : jdust_div6_lbc_tend + &
                                    lenrima(fld_type_p,halo_type_single,1)* &
                                   (tdims_s%k_len) -1)

so2_lbc_tend(1:lenrima(fld_type_p,halo_type_single,1), &
                                             tdims_s%k_start:tdims_s%k_end) &
  => d1(jso2_lbc_tend : jso2_lbc_tend + &
                                    lenrima(fld_type_p,halo_type_single,1)* &
                                   (tdims_s%k_len) -1)

dms_lbc_tend(1:lenrima(fld_type_p,halo_type_single,1), &
                                             tdims_s%k_start:tdims_s%k_end) &
 => d1(jdms_lbc_tend : jdms_lbc_tend + &
                                    lenrima(fld_type_p,halo_type_single,1)* &
                                   (tdims_s%k_len) -1)

so4_aitken_lbc_tend(1:lenrima(fld_type_p,halo_type_single,1), &
                                             tdims_s%k_start:tdims_s%k_end) &
 => d1(jso4_aitken_lbc_tend : jso4_aitken_lbc_tend + &
                                    lenrima(fld_type_p,halo_type_single,1)* &
                                   (tdims_s%k_len) -1)

so4_accu_lbc_tend(1:lenrima(fld_type_p,halo_type_single,1), &
                                             tdims_s%k_start:tdims_s%k_end) &
 => d1(jso4_accu_lbc_tend : jso4_accu_lbc_tend + &
                                    lenrima(fld_type_p,halo_type_single,1)* &
                                   (tdims_s%k_len) -1)

so4_diss_lbc_tend(1:lenrima(fld_type_p,halo_type_single,1), &
                                             tdims_s%k_start:tdims_s%k_end) &
 => d1(jso4_diss_lbc_tend : jso4_diss_lbc_tend + &
                                    lenrima(fld_type_p,halo_type_single,1)* &
                                   (tdims_s%k_len) -1)

nh3_lbc_tend(1:lenrima(fld_type_p,halo_type_single,1), &
                                             tdims_s%k_start:tdims_s%k_end) &
 => d1(jnh3_lbc_tend : jnh3_lbc_tend + &
                                    lenrima(fld_type_p,halo_type_single,1)* &
                                   (tdims_s%k_len) -1)

soot_new_lbc_tend(1:lenrima(fld_type_p,halo_type_single,1), &
                                             tdims_s%k_start:tdims_s%k_end) &
 => d1(jsoot_new_lbc_tend : jsoot_new_lbc_tend + &
                                    lenrima(fld_type_p,halo_type_single,1)* &
                                   (tdims_s%k_len) -1)

soot_agd_lbc_tend(1:lenrima(fld_type_p,halo_type_single,1), &
                                             tdims_s%k_start:tdims_s%k_end) &
 => d1(jsoot_agd_lbc_tend : jsoot_agd_lbc_tend + &
                                    lenrima(fld_type_p,halo_type_single,1)* &
                                   (tdims_s%k_len) -1)

soot_cld_lbc_tend(1:lenrima(fld_type_p,halo_type_single,1), &
                                             tdims_s%k_start:tdims_s%k_end) &
 => d1(jsoot_cld_lbc_tend : jsoot_cld_lbc_tend + &
                                    lenrima(fld_type_p,halo_type_single,1)* &
                                   (tdims_s%k_len) -1)

bmass_new_lbc_tend(1:lenrima(fld_type_p,halo_type_single,1), &
                                             tdims_s%k_start:tdims_s%k_end) &
 => d1(jbmass_new_lbc_tend : jbmass_new_lbc_tend + &
                                    lenrima(fld_type_p,halo_type_single,1)* &
                                   (tdims_s%k_len) -1)

bmass_agd_lbc_tend(1:lenrima(fld_type_p,halo_type_single,1), &
                                             tdims_s%k_start:tdims_s%k_end) &
 => d1(jbmass_agd_lbc_tend : jbmass_agd_lbc_tend + &
                                    lenrima(fld_type_p,halo_type_single,1)* &
                                   (tdims_s%k_len) -1)

bmass_cld_lbc_tend(1:lenrima(fld_type_p,halo_type_single,1), &
                                             tdims_s%k_start:tdims_s%k_end) &
 => d1(jbmass_cld_lbc_tend : jbmass_cld_lbc_tend + &
                                    lenrima(fld_type_p,halo_type_single,1)* &
                                   (tdims_s%k_len) -1)

ocff_new_lbc_tend(1:lenrima(fld_type_p,halo_type_single,1), &
                                             tdims_s%k_start:tdims_s%k_end) &
 => d1(jocff_new_lbc_tend : jocff_new_lbc_tend + &
                                    lenrima(fld_type_p,halo_type_single,1)* &
                                   (tdims_s%k_len) -1)

ocff_agd_lbc_tend(1:lenrima(fld_type_p,halo_type_single,1), &
                                             tdims_s%k_start:tdims_s%k_end) &
 => d1(jocff_agd_lbc_tend : jocff_agd_lbc_tend + &
                                    lenrima(fld_type_p,halo_type_single,1)* &
                                   (tdims_s%k_len) -1)

ocff_cld_lbc_tend(1:lenrima(fld_type_p,halo_type_single,1), &
                                             tdims_s%k_start:tdims_s%k_end) &
 => d1(jocff_cld_lbc_tend : jocff_cld_lbc_tend + &
                                    lenrima(fld_type_p,halo_type_single,1)* &
                                   (tdims_s%k_len) -1)

nitr_acc_lbc_tend(1:lenrima(fld_type_p,halo_type_single,1), &
                                             tdims_s%k_start:tdims_s%k_end) &
 => d1(jnitr_acc_lbc_tend : jnitr_acc_lbc_tend + &
                                    lenrima(fld_type_p,halo_type_single,1)* &
                                   (tdims_s%k_len) -1)

nitr_diss_lbc_tend(1:lenrima(fld_type_p,halo_type_single,1), &
                                             tdims_s%k_start:tdims_s%k_end) &
 => d1(jnitr_diss_lbc_tend : jnitr_diss_lbc_tend + &
                                    lenrima(fld_type_p,halo_type_single,1)* &
                                   (tdims_s%k_len) -1)


! Oxidant concentrations from UKCA for use in HadGEM sulphur
! cycle and ammonium nitrate scheme (these are in Section 33):
IF (l_sulpc_online_oxidants .AND. l_ukca) THEN
  oh_ukca(tdims%i_start:tdims%i_end,                                    &
          tdims%j_start:tdims%j_end,                                    &
          tdims%k_start:tdims%k_end)   =>                               &
               d1(joh_ukca(tdims%k_start) : joh_ukca(tdims%k_start)+    &
                field_length(theta_points,no_halo,tdims%k_len) -1)

  h2o2_ukca(tdims_s%i_start:tdims_s%i_end,                              &
            tdims_s%j_start:tdims_s%j_end,                              &
            tdims_s%k_start:tdims_s%k_end)   =>                         &
            d1(jh2o2_ukca(tdims%k_start): jh2o2_ukca(tdims%k_start)+    &
                field_length(theta_points,single_halo,tdims%k_len) -1)

  ho2_ukca(tdims%i_start:tdims%i_end,                                   &
           tdims%j_start:tdims%j_end,                                   &
           tdims%k_start:tdims%k_end)  =>                               &
                d1(jho2_ukca(tdims%k_start): jho2_ukca(tdims%k_start)+  &
                field_length(theta_points,no_halo,tdims%k_len) -1)

  o3_ukca (tdims_s%i_start:tdims_s%i_end,                               &
           tdims_s%j_start:tdims_s%j_end,                               &
           tdims_s%k_start:tdims_s%k_end)   =>                          &
                d1(jo3_ukca(tdims%k_start) : jo3_ukca(tdims%k_start)+   &
                field_length(theta_points,single_halo,tdims%k_len) -1)

  hno3_ukca(tdims_s%i_start:tdims_s%i_end,                              &
            tdims_s%j_start:tdims_s%j_end,                              &
            tdims_s%k_start:tdims_s%k_end)  =>                          &
                d1(jhno3_ukca(tdims%k_start):jhno3_ukca(tdims%k_start)+ &
                field_length(theta_points,single_halo,tdims%k_len) -1)
ELSE
  oh_ukca   => dummy_fld3
  h2o2_ukca => dummy_fld3
  ho2_ukca  => dummy_fld3
  o3_ukca   => dummy_fld3
  hno3_ukca => dummy_fld3
END IF

! 1.1 point tracer fields to D1

      ! find out how many tracers are active
nActiveTracers=0
DO nTracer=a_tracer_first,a_tracer_last
  IF (si(nTracer,33,atmos_im) /= 1) THEN
    nActiveTracers = nActiveTracers+1
  END IF
END DO ! nTracer

IF (nActiveTracers /= 0) THEN

  ! set the pointer to the appropriate section of D1
  tracer(tdims_s%i_start:tdims_s%i_end,                     &
         tdims_s%j_start:tdims_s%j_end,                     &
         tdims_s%k_start:tdims_s%k_end,                     &
         1:nActiveTracers)                                  &
    => d1( jtracer(tdims_s%k_start,a_tracer_first) :        &
           jtracer(tdims_s%k_start,a_tracer_first) +        &
    field_length(theta_points,single_halo,tdims%k_len)*nActiveTracers -1)

ELSE
  ! or set it to something non-null if there are no active tracers
  tracer => dummy_fld4
END IF

! do the same for section 34 (UKCA) tracers
nActiveTracers=0
DO nTracer=a_ukca_first,a_ukca_last
  IF (si(nTracer,34,atmos_im) /= 1) THEN
    nActiveTracers = nActiveTracers+1
  END IF
END DO ! nTracer

IF (nActiveTracers /= 0) THEN
  ! sec 34 ukca tracers (up to a_ukca_last) havs single point halos
  tracer_ukca(tdims_s%i_start:tdims_s%i_end,                &
              tdims_s%j_start:tdims_s%j_end,                &
              tdims_s%k_start:tdims_s%k_end,                &
              1:nActiveTracers)                             &
    => d1( jtr_ukca(tdims_s%k_start,a_ukca_first) :         &
           jtr_ukca(tdims_s%k_start,a_ukca_first) +         &
    field_length(theta_points,single_halo,tdims%k_len)*nActiveTracers -1)

ELSE
  tracer_ukca => dummy_fld4
END IF

! find out how many free tracer LBCs are active
nActiveTracers=0
DO nTracer=1,tr_lbc_vars
  IF (si(nTracer,36,atmos_im) /= 1) THEN
    nActiveTracers = nActiveTracers+1
  END IF
END DO ! nTracer

IF (nActiveTracers /= 0) THEN
  ! set the pointer to the appropriate section of D1
  ! sec 36 tracers have extended halos
  tracer_lbc(1:lenrima(fld_type_p,halo_type_extended,1),           &
             tdims_l%k_start:tdims_l%k_end,                        &
             1:nActiveTracers)                                     &
     => d1(jtracer_lbc(1) :                                        &
   jtracer_lbc(nActiveTracers) +                                   &
   (lenrima(fld_type_p,halo_type_extended,1)*tdims_l%k_len ) )

  tracer_lbc_tend(1:lenrima(fld_type_p,halo_type_extended,1),      &
                  tdims_l%k_start:tdims_l%k_end,                   &
                  1:nActiveTracers)                                &
     => d1(jtracer_lbc_tend(1) :                                   &
   jtracer_lbc_tend(nActiveTracers) +                              &
   (lenrima(fld_type_p,halo_type_extended,1)*tdims_l%k_len ) )
ELSE
  ! or set it to something non-null if there are no active tracer LBCs
  tracer_lbc => dummy_fld3
  tracer_lbc_tend => dummy_fld3
END IF

! find out how many UKCA tracer LBCs are active
nActiveTracers=0
DO nTracer=1,tr_lbc_ukca
  IF (si(nTracer,37,atmos_im) /= 1) THEN
    nActiveTracers = nActiveTracers+1
  END IF
END DO ! nTracer

IF (nActiveTracers /= 0) THEN
  ! set the pointer to the appropriate section of D1
  ! sec 37 tracers have extended halos
  tracer_ukca_lbc(1:lenrima(fld_type_p,halo_type_extended,1),      &
                  tdims_l%k_start:tdims_l%k_end,                   &
                  1:nActiveTracers)                                &
      => d1(jtr_ukca_lbc(1) :                                      &
   jtr_ukca_lbc(nActiveTracers) +                                  &
   (lenrima(fld_type_p,halo_type_extended,1)*tdims_l%k_len) )

  tracer_ukca_lbc_tend(1:lenrima(fld_type_p,halo_type_extended,1), &
                       tdims_l%k_start:tdims_l%k_end,              &
                       1:nActiveTracers)                           &
      => d1(jtr_ukca_lbc_tend(1) :                                 &
   jtr_ukca_lbc_tend(nActiveTracers) +                             &
   (lenrima(fld_type_p,halo_type_extended,1)*tdims_l%k_len) )

ELSE
  ! set to something non-null if there are no active ukca tracer LBCs
  tracer_ukca_lbc => dummy_fld3
  tracer_ukca_lbc_tend => dummy_fld3

END IF
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE Set_Atm_Fields
