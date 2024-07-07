! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Update the Lateral Boundary Conditions (LBCs) of LAM fields
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: LBC Input

SUBROUTINE update_lam_lbcs(                                       &
  r_rho_levels, r_theta_levels,                                   &
  row_length,rows,n_rows,                                         &
  tr_vars,tr_lbc_vars,tr_levels,                                  &
  a_max_trvars,A_tr_active_lbc_index,                             &
  tr_ukca,tr_lbc_ukca,                                            &
  offx,offy,halo_i,halo_j,at_extremity,                           &
  L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,                          &
  L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc, L_pc2_lbc,   &
  L_murk, L_murk_lbc,                                             &
  L_LBC_balance, L_int_uvw_lbc,                                   &
  L_dust_div1, L_dust_div1_lbc, L_dust_div2, L_dust_div2_lbc,     &
  L_dust_div3, L_dust_div3_lbc, L_dust_div4, L_dust_div4_lbc,     &
  L_dust_div5, L_dust_div5_lbc, L_dust_div6, L_dust_div6_lbc,     &
  L_so2,L_so2_lbc,L_dms,L_dms_lbc,L_so4_aitken,L_so4_aitken_lbc,  &
  L_so4_accu, L_so4_accu_lbc, L_so4_diss, L_so4_diss_lbc,         &
  L_nh3, L_nh3_lbc,                                               &
  L_soot_new_lbc,  L_soot_agd_lbc, L_soot_cld_lbc, L_bmass_new,   &
  L_bmass_new_lbc,                                                &
  L_bmass_agd, L_bmass_agd_lbc, L_bmass_cld, L_bmass_cld_lbc,     &
  L_ocff_new_lbc, L_ocff_agd_lbc,  L_ocff_cld_lbc,                &
  L_nitr_acc, L_nitr_acc_lbc, L_nitr_diss, L_nitr_diss_lbc,       &
  rimwidth,rimweights,lenrim,lbc_size,lbc_start,                  &
  theta_lbc, q_lbc, qcl_lbc, qcf_lbc,                             &
  qcf2_lbc, qrain_lbc, qgraup_lbc,                                &
  cf_bulk_lbc, cf_liquid_lbc, cf_frozen_lbc,                      &
  rho_lbc, exner_lbc,                                             &
  u_lbc, v_lbc, w_lbc, u_adv_lbc, v_adv_lbc, w_adv_lbc,           &
  murk_lbc,                                                       &
  dust_div1_lbc, dust_div2_lbc, dust_div3_lbc,                    &
  dust_div4_lbc, dust_div5_lbc, dust_div6_lbc,                    &
  so2_lbc, dms_lbc, so4_aitken_lbc, so4_accu_lbc, so4_diss_lbc,   &
  nh3_lbc, soot_new_lbc, soot_agd_lbc, soot_cld_lbc,              &
  bmass_new_lbc, bmass_agd_lbc, bmass_cld_lbc,                    &
  ocff_new_lbc, ocff_agd_lbc, ocff_cld_lbc,                       &
  nitr_acc_lbc, nitr_diss_lbc,                                    &
  tracer_lbc, tracer_ukca_lbc,                                    &
  theta, q, qcl, qcf, qcf2, qrain, qgraup,                        &
  cf_bulk, cf_liquid, cf_frozen,                                  &
  rho, exner, u, v, w, u_adv, v_adv, w_adv, murk,                 &
  dust_div1, dust_div2, dust_div3,                                &
  dust_div4, dust_div5, dust_div6,                                &
  so2, dms, so4_aitken, so4_accu, so4_diss,                       &
  nh3, soot_new, soot_agd, soot_cld,                              &
  bmass_new, bmass_agd, bmass_cld,                                &
  ocff_new, ocff_agd, ocff_cld,                                   &
  nitr_acc, nitr_diss,                                            &
  delta_phi, delta_lambda,                                        &
  base_phi, base_lambda,                                          &
  datastart, lat_rot_NP,                                          &
  global_row_length, global_rows,                                 &
  free_tracers, tracer_ukca                                       &
    )

USE yomhook, ONLY: lhook, dr_hook

USE parkind1, ONLY: jprb, jpim

USE atm_fields_bounds_mod, ONLY : array_dims,                     &
                                  pdims, pdims_s,                 &
                                  tdims, tdims_l, tdims_s,        &
                                  udims, udims_l, udims_s,        &
                                  vdims, vdims_l, vdims_s,        &
                                  wdims, wdims_l, wdims_s
USE cloud_inputs_mod, ONLY: i_cld_vn
USE pc2_constants_mod, ONLY: i_cld_pc2
USE nlsizes_namelist_mod, ONLY: model_levels

USE run_aerosol_mod, ONLY: L_ocff_new, L_ocff_agd, L_ocff_cld,  &
                         ! true if fossil fuel aerosol active
                           L_soot_new, L_soot_agd, L_soot_cld
                         ! true if soot active

USE umPrintMgr

USE UM_ParParams

USE eg_balance_lbc_values_mod

USE ukca_tracer_stash, ONLY: ukca_tr_active_lbc_index

IMPLICIT NONE
! Updates all the input fields with their LBCs
!

! Parameters required for dimensioning some of the arguments
! Arguments

INTEGER ::                                                        &
  row_length                                                      &
                    ! IN : Length of a model row
, rows                                                            &
                    ! IN : Number of rows for theta,u fields
, n_rows                                                          &
                    ! IN : Number of rows for v fields
, offx                                                            &
                    ! IN : Size of "single" halo (EW direction)
, offy                                                            &
                    ! IN : Size of "single" halo (NS direction)
, halo_i                                                          &
                    ! IN : Size of extended halo (EW direction)
, halo_j                                                          &
                    ! IN : Size of extended halo (NS direction)
, rimwidth                                                        &
                    ! IN : Size of boundary region
, lenrim(Nfld_max,NHalo_max)                                      &
                    ! IN : Size of single level of LBC
, lbc_size(4,Nfld_max,NHalo_max)                                  &
                    ! IN : Size of a side of a LBC
, lbc_start(4,Nfld_max,NHalo_max)                                 &
                    ! IN : Start of a side in a LBC
, datastart(3)                                                    &
                    ! IN : position of first data point for PE
, global_row_length                                               &
                    ! IN : no. of columns in LAM
, global_rows                                                     &
                    ! IN : no. or rows in LAM
, tr_vars                                                         &
                    ! Number of tracer prognostics
, tr_ukca                                                         &
                    ! Number of UKCA tracer prognostics
, tr_lbc_vars                                                     &
                    ! Number of tracer lbcs
, tr_lbc_ukca                                                     &
                    ! Number of UKCA tracer lbcs
, tr_levels                                                       &
                    ! Number of tracer levels
, a_max_trvars                                                    &
                    ! Max number of tracer prognostics
, A_tr_active_lbc_index(a_max_trvars)
                    ! >0 if lbc active

LOGICAL, INTENT (IN) ::                                           &
  at_extremity(4)                                                 &
                ! IN : At an edge?
, L_mcr_qcf2                                                      &
                ! true if prognostic 2nd cloud ice active
, L_mcr_qrain                                                     &
                ! true if prognostic rain active
, L_mcr_qgraup                                                    &
                ! true if prognostic graupel active
, L_mcr_qcf2_lbc                                                  &
                    ! true if prognostic 2nd cloud ice in lbcs
, L_mcr_qrain_lbc                                                 &
                    ! true if prognostic rain in lbcs
, L_mcr_qgraup_lbc                                                &
                    ! true if prognostic graupel in lbcs
, L_pc2_lbc                                                       &
                    ! true if prognostic cloud fracs in lbcs
, L_murk                                                          &
                    ! true if murk aerosol active
, L_murk_lbc                                                      &
                    ! true if murk aerosol lbcs active
, L_int_uvw_lbc                                                   &
                    ! true if using interpolated advecting winds
                    !  in lateral boundaries
, L_LBC_balance                                                   &
                    ! true if imposing balance in vertical
                    !      momentum equation in LBC regions
, L_dust_div1                                                     &
, L_dust_div2                                                     &
, L_dust_div3                                                     &
, L_dust_div4                                                     &
, L_dust_div5                                                     &
, L_dust_div6                                                     &
                    ! true if DUST active
, L_dust_div1_lbc                                                 &
, L_dust_div2_lbc                                                 &
, L_dust_div3_lbc                                                 &
, L_dust_div4_lbc                                                 &
, L_dust_div5_lbc                                                 &
, L_dust_div6_lbc                                                 &
                    ! true if DUST lbcs active
, L_so2                                                           &
                    ! true if SO2 active
, L_so2_lbc                                                       &
                    ! true if SO2 lbcs active
, L_dms                                                           &
                    ! true if DMS active
, L_dms_lbc                                                       &
                    ! true if DMS lbcs aactive
, L_so4_aitken                                                    &
                    ! true if SO4_AITKEN active
, L_so4_aitken_lbc                                                &
                    ! true if SO4_AITKEN lbcs active
, L_so4_accu                                                      &
                    ! true if SO4_ACCU active
, L_so4_accu_lbc                                                  &
                    ! true if SO4_ACCU lbcs active
, L_so4_diss                                                      &
                    ! true if SO4_DISS active
, L_so4_diss_lbc                                                  &
                    ! true if SO4_DISS lbcs active
, L_nh3                                                           &
                    ! true if NH3 active
, L_nh3_lbc                                                       &
                    ! true if NH3 lbcs active
, L_soot_new_lbc                                                  &
, L_soot_agd_lbc                                                  &
, L_soot_cld_lbc                                                  &
                    ! true if soot lbcs active
, L_bmass_new                                                     &
, L_bmass_agd                                                     &
, L_bmass_cld                                                     &
                    ! true if biomass active
, L_bmass_new_lbc                                                 &
, L_bmass_agd_lbc                                                 &
, L_bmass_cld_lbc                                                 &
                    ! true if biomass lbcs active
, L_ocff_new_lbc                                                  &
, L_ocff_agd_lbc                                                  &
, L_ocff_cld_lbc                                                  &
                    ! true if fossil fuel aerosol lbcs active
, L_nitr_acc                                                      &
, L_nitr_diss                                                     &
                    ! true if nitrate aerosol active
, L_nitr_acc_lbc                                                  &
, L_nitr_diss_lbc   ! true if nitrate aerosol lbcs active


REAL, INTENT(IN) ::                                               &
     ! vertical co-ordinate information
  r_theta_levels(1-halo_i:row_length+halo_i,                      &
                 1-halo_j:rows+halo_j,0:model_levels)             &
, r_rho_levels(1-halo_i:row_length+halo_i,                        &
               1-halo_j:rows+halo_j, model_levels)

REAL, INTENT (IN) ::                                              &
  delta_phi                                                       &
              ! grid spacing (latitude)
, delta_lambda                                                    &
              ! grid spacing (longitude)
, base_phi                                                        &
              ! first latitude
, base_lambda                                                     &
              ! first longitude
, lat_rot_NP                                                      &
              ! lat of rotated pole
, rimweights(rimwidth)  ! weight to apply to LBC

REAL, INTENT (INOUT) ::                                           &
  theta_lbc(lenrim(fld_type_p,halo_type_extended),                &
            tdims_s%k_start:tdims_s%k_end)                        &
, q_lbc(lenrim(fld_type_p,halo_type_extended),                    &
        tdims_l%k_start:tdims_l%k_end)                            &
, qcl_lbc(lenrim(fld_type_p,halo_type_extended),                  &
          tdims_l%k_start:tdims_l%k_end)                          &
, qcf_lbc(lenrim(fld_type_p,halo_type_extended),                  &
          tdims_l%k_start:tdims_l%k_end)                          &
, qcf2_lbc(lenrim(fld_type_p,halo_type_extended),                 &
          tdims_l%k_start:tdims_l%k_end)                          &
, qrain_lbc(lenrim(fld_type_p,halo_type_extended),                &
          tdims_l%k_start:tdims_l%k_end)                          &
, qgraup_lbc(lenrim(fld_type_p,halo_type_extended),               &
          tdims_l%k_start:tdims_l%k_end)                          &
, cf_bulk_lbc(lenrim(fld_type_p,halo_type_extended),              &
          tdims_l%k_start:tdims_l%k_end)                          &
, cf_liquid_lbc(lenrim(fld_type_p,halo_type_extended),            &
          tdims_l%k_start:tdims_l%k_end)                          &
, cf_frozen_lbc(lenrim(fld_type_p,halo_type_extended),            &
          tdims_l%k_start:tdims_l%k_end)                          &
, rho_lbc(lenrim(fld_type_p,halo_type_extended),                  &
          pdims_s%k_start:pdims_s%k_end)                          &
, exner_lbc(lenrim(fld_type_p,halo_type_extended),                &
            pdims_s%k_start:pdims_s%k_end+1)                      &
, u_lbc(lenrim(fld_type_u,halo_type_extended),                    &
        udims_s%k_start:udims_s%k_end)                            &
, v_lbc(lenrim(fld_type_v,halo_type_extended),                    &
        vdims_s%k_start:vdims_s%k_end)                            &
, w_lbc(lenrim(fld_type_p,halo_type_extended),                    &
        wdims_s%k_start:wdims_s%k_end)                            &
, u_adv_lbc(lenrim(fld_type_u,halo_type_extended),                &
            udims_s%k_start:udims_s%k_end)                        &
, v_adv_lbc(lenrim(fld_type_v,halo_type_extended),                &
            vdims_s%k_start:vdims_s%k_end)                        &
, w_adv_lbc(lenrim(fld_type_p,halo_type_extended),                &
            wdims_s%k_start:wdims_s%k_end)                        &
, murk_lbc(lenrim(fld_type_p,halo_type_single),                   &
            tdims_s%k_start:tdims_s%k_end)                        &
, dust_div1_lbc(lenrim(fld_type_p,halo_type_single),              &
            tdims_s%k_start:tdims_s%k_end)                        &
, dust_div2_lbc(lenrim(fld_type_p,halo_type_single),              &
            tdims_s%k_start:tdims_s%k_end)                        &
, dust_div3_lbc(lenrim(fld_type_p,halo_type_single),              &
            tdims_s%k_start:tdims_s%k_end)                        &
, dust_div4_lbc(lenrim(fld_type_p,halo_type_single),              &
            tdims_s%k_start:tdims_s%k_end)                        &
, dust_div5_lbc(lenrim(fld_type_p,halo_type_single),              &
            tdims_s%k_start:tdims_s%k_end)                        &
, dust_div6_lbc(lenrim(fld_type_p,halo_type_single),              &
            tdims_s%k_start:tdims_s%k_end)                        &
, so2_lbc(lenrim(fld_type_p,halo_type_single),                    &
            tdims_s%k_start:tdims_s%k_end)                        &
, dms_lbc(lenrim(fld_type_p,halo_type_single),                    &
            tdims_s%k_start:tdims_s%k_end)                        &
, so4_aitken_lbc(lenrim(fld_type_p,halo_type_single),             &
            tdims_s%k_start:tdims_s%k_end)                        &
, so4_accu_lbc(lenrim(fld_type_p,halo_type_single),               &
            tdims_s%k_start:tdims_s%k_end)                        &
, so4_diss_lbc(lenrim(fld_type_p,halo_type_single),               &
            tdims_s%k_start:tdims_s%k_end)                        &
, nh3_lbc(lenrim(fld_type_p,halo_type_single),                    &
            tdims_s%k_start:tdims_s%k_end)                        &
, soot_new_lbc(lenrim(fld_type_p,halo_type_single),               &
            tdims_s%k_start:tdims_s%k_end)                        &
, soot_agd_lbc(lenrim(fld_type_p,halo_type_single),               &
            tdims_s%k_start:tdims_s%k_end)                        &
, soot_cld_lbc(lenrim(fld_type_p,halo_type_single),               &
            tdims_s%k_start:tdims_s%k_end)                        &
, bmass_new_lbc(lenrim(fld_type_p,halo_type_single),              &
            tdims_s%k_start:tdims_s%k_end)                        &
, bmass_agd_lbc(lenrim(fld_type_p,halo_type_single),              &
            tdims_s%k_start:tdims_s%k_end)                        &
, bmass_cld_lbc(lenrim(fld_type_p,halo_type_single),              &
            tdims_s%k_start:tdims_s%k_end)                        &
, ocff_new_lbc(lenrim(fld_type_p,halo_type_single),               &
            tdims_s%k_start:tdims_s%k_end)                        &
, ocff_agd_lbc(lenrim(fld_type_p,halo_type_single),               &
            tdims_s%k_start:tdims_s%k_end)                        &
, ocff_cld_lbc(lenrim(fld_type_p,halo_type_single),               &
            tdims_s%k_start:tdims_s%k_end)                        &
, nitr_acc_lbc(lenrim(fld_type_p,halo_type_single),               &
            tdims_s%k_start:tdims_s%k_end)                        &
, nitr_diss_lbc(lenrim(fld_type_p,halo_type_single),              &
            tdims_s%k_start:tdims_s%k_end)                        &
, tracer_lbc(lenrim(fld_type_p,halo_type_extended),               &
            tdims_l%k_start:tdims_l%k_end,tr_lbc_vars)            &
, tracer_ukca_lbc(lenrim(fld_type_p,halo_type_extended),          &
            tdims_l%k_start:tdims_l%k_end,tr_lbc_ukca)

REAL, INTENT (INOUT) ::                                           &
  theta(tdims_s%i_start:tdims_s%i_end,                            &
        tdims_s%j_start:tdims_s%j_end,                            &
        tdims_s%k_start:tdims_s%k_end)                            &
, q(tdims_s%i_start:tdims_s%i_end,                                &
    tdims_s%j_start:tdims_s%j_end,                                &
    tdims_s%k_start:tdims_s%k_end)                                &
, qcl(tdims_s%i_start:tdims_s%i_end,                              &
      tdims_s%j_start:tdims_s%j_end,                              &
      tdims_s%k_start:tdims_s%k_end)                              &
, qcf(tdims_s%i_start:tdims_s%i_end,                              &
      tdims_s%j_start:tdims_s%j_end,                              &
      tdims_s%k_start:tdims_s%k_end)                              &
, qcf2(tdims_s%i_start:tdims_s%i_end,                             &
       tdims_s%j_start:tdims_s%j_end,                             &
       tdims_s%k_start:tdims_s%k_end)                             &
, qrain(tdims_s%i_start:tdims_s%i_end,                            &
        tdims_s%j_start:tdims_s%j_end,                            &
        tdims_s%k_start:tdims_s%k_end)                            &
, qgraup(tdims_s%i_start:tdims_s%i_end,                           &
         tdims_s%j_start:tdims_s%j_end,                           &
         tdims_s%k_start:tdims_s%k_end)                           &
, cf_bulk(tdims_l%i_start:tdims_l%i_end,                          &
        tdims_l%j_start:tdims_l%j_end,                            &
        tdims_l%k_start:tdims_l%k_end)                            &
, cf_liquid(tdims_l%i_start:tdims_l%i_end,                        &
        tdims_l%j_start:tdims_l%j_end,                            &
        tdims_l%k_start:tdims_l%k_end)                            &
, cf_frozen(tdims_l%i_start:tdims_l%i_end,                        &
        tdims_l%j_start:tdims_l%j_end,                            &
        tdims_l%k_start:tdims_l%k_end)                            &
, rho(pdims_s%i_start:pdims_s%i_end,                              &
    pdims_s%j_start:pdims_s%j_end,                                &
    pdims_s%k_start:pdims_s%k_end)                                &
, exner(pdims_s%i_start:pdims_s%i_end,                            &
    pdims_s%j_start:pdims_s%j_end,                                &
    pdims_s%k_start:pdims_s%k_end+1)                              &
, u(udims_s%i_start:udims_s%i_end,                                &
        udims_s%j_start:udims_s%j_end,                            &
        udims_s%k_start:udims_s%k_end)                            &
, v(vdims_s%i_start:vdims_s%i_end,                                &
        vdims_s%j_start:vdims_s%j_end,                            &
        vdims_s%k_start:vdims_s%k_end)                            &
, w(wdims_s%i_start:wdims_s%i_end,                                &
    wdims_s%j_start:wdims_s%j_end,                                &
    wdims_s%k_start:wdims_s%k_end)                                &
, u_adv(udims_l%i_start:udims_l%i_end,                            &
        udims_l%j_start:udims_l%j_end,                            &
        udims_l%k_start:udims_l%k_end)                            &
, v_adv(vdims_l%i_start:vdims_l%i_end,                            &
        vdims_l%j_start:vdims_l%j_end,                            &
        vdims_l%k_start:vdims_l%k_end)                            &
, w_adv(wdims_l%i_start:wdims_l%i_end,                            &
    wdims_l%j_start:wdims_l%j_end,                                &
    wdims_l%k_start:wdims_l%k_end)                                &
, murk(tdims_s%i_start:tdims_s%i_end,                             &
    tdims_s%j_start:tdims_s%j_end,                                &
    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT (INOUT) ::                                           &
  dust_div1(tdims_s%i_start:tdims_s%i_end,                        &
    tdims_s%j_start:tdims_s%j_end,                                &
    tdims_s%k_start:tdims_s%k_end)                                &
, dust_div2(tdims_s%i_start:tdims_s%i_end,                        &
    tdims_s%j_start:tdims_s%j_end,                                &
    tdims_s%k_start:tdims_s%k_end)                                &
, dust_div3(tdims_s%i_start:tdims_s%i_end,                        &
    tdims_s%j_start:tdims_s%j_end,                                &
    tdims_s%k_start:tdims_s%k_end)                                &
, dust_div4(tdims_s%i_start:tdims_s%i_end,                        &
    tdims_s%j_start:tdims_s%j_end,                                &
    tdims_s%k_start:tdims_s%k_end)                                &
, dust_div5(tdims_s%i_start:tdims_s%i_end,                        &
    tdims_s%j_start:tdims_s%j_end,                                &
    tdims_s%k_start:tdims_s%k_end)                                &
, dust_div6(tdims_s%i_start:tdims_s%i_end,                        &
    tdims_s%j_start:tdims_s%j_end,                                &
    tdims_s%k_start:tdims_s%k_end)                                &
, so2(tdims_s%i_start:tdims_s%i_end,                              &
    tdims_s%j_start:tdims_s%j_end,                                &
    tdims_s%k_start:tdims_s%k_end)                                &
, dms(tdims_s%i_start:tdims_s%i_end,                              &
    tdims_s%j_start:tdims_s%j_end,                                &
   tdims_s%k_start:tdims_s%k_end)                                 &
, so4_aitken(tdims_s%i_start:tdims_s%i_end,                       &
    tdims_s%j_start:tdims_s%j_end,                                &
    tdims_s%k_start:tdims_s%k_end)                                &
, so4_accu(tdims_s%i_start:tdims_s%i_end,                         &
    tdims_s%j_start:tdims_s%j_end,                                &
    tdims_s%k_start:tdims_s%k_end)                                &
, so4_diss(tdims_s%i_start:tdims_s%i_end,                         &
    tdims_s%j_start:tdims_s%j_end,                                &
    tdims_s%k_start:tdims_s%k_end)                                &
, nh3(tdims_s%i_start:tdims_s%i_end,                              &
    tdims_s%j_start:tdims_s%j_end,                                &
    tdims_s%k_start:tdims_s%k_end)                                &
, soot_new(tdims_s%i_start:tdims_s%i_end,                         &
    tdims_s%j_start:tdims_s%j_end,                                &
    tdims_s%k_start:tdims_s%k_end)                                &
, soot_agd(tdims_s%i_start:tdims_s%i_end,                         &
    tdims_s%j_start:tdims_s%j_end,                                &
    tdims_s%k_start:tdims_s%k_end)                                &
, soot_cld(tdims_s%i_start:tdims_s%i_end,                         &
    tdims_s%j_start:tdims_s%j_end,                                &
    tdims_s%k_start:tdims_s%k_end)                                &
, bmass_new(tdims_s%i_start:tdims_s%i_end,                        &
    tdims_s%j_start:tdims_s%j_end,                                &
    tdims_s%k_start:tdims_s%k_end)                                &
, bmass_agd(tdims_s%i_start:tdims_s%i_end,                        &
    tdims_s%j_start:tdims_s%j_end,                                &
    tdims_s%k_start:tdims_s%k_end)                                &
, bmass_cld(tdims_s%i_start:tdims_s%i_end,                        &
    tdims_s%j_start:tdims_s%j_end,                                &
    tdims_s%k_start:tdims_s%k_end)                                &
, ocff_new(tdims_s%i_start:tdims_s%i_end,                         &
    tdims_s%j_start:tdims_s%j_end,                                &
    tdims_s%k_start:tdims_s%k_end)                                &
, ocff_agd(tdims_s%i_start:tdims_s%i_end,                         &
    tdims_s%j_start:tdims_s%j_end,                                &
    tdims_s%k_start:tdims_s%k_end)                                &
, ocff_cld(tdims_s%i_start:tdims_s%i_end,                         &
    tdims_s%j_start:tdims_s%j_end,                                &
    tdims_s%k_start:tdims_s%k_end)                                &
, nitr_acc(tdims_s%i_start:tdims_s%i_end,                         &
    tdims_s%j_start:tdims_s%j_end,                                &
    tdims_s%k_start:tdims_s%k_end)                                &
, nitr_diss(tdims_s%i_start:tdims_s%i_end,                        &
    tdims_s%j_start:tdims_s%j_end,                                &
    tdims_s%k_start:tdims_s%k_end)                                &
, free_tracers(tdims_s%i_start:tdims_s%i_end,                     &
               tdims_s%j_start:tdims_s%j_end,                     &
               tdims_s%k_start:tdims_s%k_end, tr_vars)            &
, tracer_ukca(tdims_s%i_start:tdims_s%i_end,                      &
              tdims_s%j_start:tdims_s%j_end,                      &
              tdims_s%k_start:tdims_s%k_end,tr_ukca)

! Local variables

LOGICAL ::                                                        &
  L_Do_Boundaries                                                 &
, L_Do_Halos
INTEGER ::                                                        &
  n_rims_to_do                                                    &
, lbc_index                                                       &
, tracer

! This pointer is used to distinguish between EG and ND
! moisture field sizes (avoids changing the argument lists)


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UPDATE_LAM_LBCS'

! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

L_Do_Boundaries = .TRUE.
L_Do_Halos = .TRUE.
n_rims_to_do = rimwidth

IF (L_LBC_balance) THEN

  ! Reset LBC array values for Exner, rho and w

  CALL eg_balance_lbc_values(                                    &
   exner_lbc, rho_lbc, theta_lbc,q_lbc, w_lbc, w_adv_lbc,        &
   u_lbc, v_lbc,                                                 &
   qcl_lbc, qcf_lbc, qcf2_lbc, qrain_lbc, qgraup_lbc,            &
   L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,                        &
   r_rho_levels, r_theta_levels,                                 &
   row_length, rows, halo_i, halo_j,                             &
   lenrim(fld_type_p,halo_type_extended),                        &
   lenrim(fld_type_u,halo_type_extended),                        &
   lenrim(fld_type_v,halo_type_extended),                        &
   lbc_start(1,fld_type_p,halo_type_extended),                   &
   lbc_start(1,fld_type_u,halo_type_extended),                   &
   lbc_start(1,fld_type_v,halo_type_extended),                   &
   rimwidth, n_rims_to_do, rimweights, at_extremity,             &
   delta_phi, delta_lambda,                                      &
   base_phi, base_lambda,                                        &
   datastart, lat_rot_np,                                        &
   global_row_length, global_rows)

END IF

! Theta
! DEPENDS ON: set_lateral_boundaries
IF (rimwidth > 0) THEN
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_s%halo_i,                                               &
   tdims_s%halo_j,                                               &
   tdims_s%k_len,                                                &
   fld_type_p,theta,                                             &
   lenrim(fld_type_p,halo_type_extended),                        &
   lbc_size(1,fld_type_p,halo_type_extended),                    &
   lbc_start(1,fld_type_p,halo_type_extended),                   &
   halo_i,halo_j,theta_lbc,                                      &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)

! Q
! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_s%halo_i,                                               &
   tdims_s%halo_j,                                               &
   tdims_s%k_len,                                                &
   fld_type_p,q,                                                 &
   lenrim(fld_type_p,halo_type_extended),                        &
   lbc_size(1,fld_type_p,halo_type_extended),                    &
   lbc_start(1,fld_type_p,halo_type_extended),                   &
   halo_i,halo_j,q_lbc,                                          &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)

! QCL
! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_s%halo_i,                                               &
   tdims_s%halo_j,                                               &
   tdims_s%k_len,                                                &
   fld_type_p,qcl,                                               &
   lenrim(fld_type_p,halo_type_extended),                        &
   lbc_size(1,fld_type_p,halo_type_extended),                    &
   lbc_start(1,fld_type_p,halo_type_extended),                   &
   halo_i,halo_j,qcl_lbc,                                        &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)

! QCF
! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_s%halo_i,                                               &
   tdims_s%halo_j,                                               &
   tdims_s%k_len,                                                &
   fld_type_p,qcf,                                               &
   lenrim(fld_type_p,halo_type_extended),                        &
   lbc_size(1,fld_type_p,halo_type_extended),                    &
   lbc_start(1,fld_type_p,halo_type_extended),                   &
   halo_i,halo_j,qcf_lbc,                                        &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)
END IF


! QCF2 - Update if prognostic is active and data is in lbcs
IF (L_mcr_qcf2 .AND. L_mcr_qcf2_lbc .AND. (rimwidth > 0) ) THEN
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_s%halo_i,                                               &
   tdims_s%halo_j,                                               &
   tdims_s%k_len,                                                &
   fld_type_p,qcf2,                                              &
   lenrim(fld_type_p,halo_type_extended),                        &
   lbc_size(1,fld_type_p,halo_type_extended),                    &
   lbc_start(1,fld_type_p,halo_type_extended),                   &
   halo_i,halo_j,qcf2_lbc,                                       &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)
END IF

! QRAIN - Update if prognostic is active and data is in lbcs
IF (L_mcr_qrain .AND. L_mcr_qrain_lbc .AND. (rimwidth > 0) ) THEN
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_s%halo_i,                                               &
   tdims_s%halo_j,                                               &
   tdims_s%k_len,                                                &
   fld_type_p,qrain,                                             &
   lenrim(fld_type_p,halo_type_extended),                        &
   lbc_size(1,fld_type_p,halo_type_extended),                    &
   lbc_start(1,fld_type_p,halo_type_extended),                   &
   halo_i,halo_j,qrain_lbc,                                      &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)
END IF

! QGRAUP - Update if prognostic is active and data is in lbcs
IF (L_mcr_qgraup .AND. L_mcr_qgraup_lbc .AND. rimwidth > 0) THEN
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_s%halo_i,                                               &
   tdims_s%halo_j,                                               &
   tdims_s%k_len,                                                &
   fld_type_p,qgraup,                                            &
   lenrim(fld_type_p,halo_type_extended),                        &
   lbc_size(1,fld_type_p,halo_type_extended),                    &
   lbc_start(1,fld_type_p,halo_type_extended),                   &
   halo_i,halo_j,qgraup_lbc,                                     &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)
END IF

! Cloud fractions - Update if prognostics are active and data in lbcs
IF (i_cld_vn == i_cld_pc2 .AND. L_pc2_lbc .AND. rimwidth>0) THEN

  ! CF_BULK
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_l%halo_i,                                               &
   tdims_l%halo_j,                                               &
   tdims_l%k_len,                                                &
   fld_type_p,cf_bulk,                                           &
   lenrim(fld_type_p,halo_type_extended),                        &
   lbc_size(1,fld_type_p,halo_type_extended),                    &
   lbc_start(1,fld_type_p,halo_type_extended),                   &
   halo_i,halo_j,cf_bulk_lbc,                                    &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)

  ! CF_LIQUID
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_l%halo_i,                                               &
   tdims_l%halo_j,                                               &
   tdims_l%k_len,                                                &
   fld_type_p,cf_liquid,                                         &
   lenrim(fld_type_p,halo_type_extended),                        &
   lbc_size(1,fld_type_p,halo_type_extended),                    &
   lbc_start(1,fld_type_p,halo_type_extended),                   &
   halo_i,halo_j,cf_liquid_lbc,                                  &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)

  ! CF_FROZEN
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_l%halo_i,                                               &
   tdims_l%halo_j,                                               &
   tdims_l%k_len,                                                &
   fld_type_p,cf_frozen,                                         &
   lenrim(fld_type_p,halo_type_extended),                        &
   lbc_size(1,fld_type_p,halo_type_extended),                    &
   lbc_start(1,fld_type_p,halo_type_extended),                   &
   halo_i,halo_j,cf_frozen_lbc,                                  &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)

END IF ! i_cld_pc2 .AND. L_pc2_lbc

! RHO
! DEPENDS ON: set_lateral_boundaries
IF ( rimwidth > 0 ) THEN
  CALL set_lateral_boundaries(                                   &
   pdims%i_len,                                                  &
   pdims%j_len,                                                  &
   pdims_s%halo_i,                                               &
   pdims_s%halo_j,                                               &
   pdims_s%k_len,                                                &
   fld_type_p,rho,                                               &
   lenrim(fld_type_p,halo_type_extended),                        &
   lbc_size(1,fld_type_p,halo_type_extended),                    &
   lbc_start(1,fld_type_p,halo_type_extended),                   &
   halo_i,halo_j,rho_lbc,                                        &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)

! EXNER
! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   pdims%i_len,                                                  &
   pdims%j_len,                                                  &
   pdims_s%halo_i,                                               &
   pdims_s%halo_j,                                               &
   pdims_s%k_end - pdims_s%k_start + 2,                          &
   fld_type_p,exner,                                             &
   lenrim(fld_type_p,halo_type_extended),                        &
   lbc_size(1,fld_type_p,halo_type_extended),                    &
   lbc_start(1,fld_type_p,halo_type_extended),                   &
   halo_i,halo_j,exner_lbc,                                      &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)

! U
! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   udims%i_len,                                                  &
   udims%j_len,                                                  &
   udims_s%halo_i,                                               &
   udims_s%halo_j,                                               &
   udims_s%k_len,                                                &
   fld_type_u,u,                                                 &
   lenrim(fld_type_u,halo_type_extended),                        &
   lbc_size(1,fld_type_u,halo_type_extended),                    &
   lbc_start(1,fld_type_u,halo_type_extended),                   &
   halo_i, halo_j,u_lbc,                                         &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)

! V
! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   vdims%i_len,                                                  &
   vdims%j_len,                                                  &
   vdims_s%halo_i,                                               &
   vdims_s%halo_j,                                               &
   vdims_s%k_len,                                                &
   fld_type_v,v,                                                 &
   lenrim(fld_type_v,halo_type_extended),                        &
   lbc_size(1,fld_type_v,halo_type_extended),                    &
   lbc_start(1,fld_type_v,halo_type_extended),                   &
   halo_i, halo_j,v_lbc,                                         &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)

! W
! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   wdims%i_len,                                                  &
   wdims%j_len,                                                  &
   wdims_s%halo_i,                                               &
   wdims_s%halo_j,                                               &
   wdims_s%k_len,                                                &
   fld_type_p,w,                                                 &
   lenrim(fld_type_p,halo_type_extended),                        &
   lbc_size(1,fld_type_p,halo_type_extended),                    &
   lbc_start(1,fld_type_p,halo_type_extended),                   &
   halo_i, halo_j,w_lbc,                                         &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)
END IF

IF ( .NOT. L_int_uvw_lbc .AND. rimwidth > 0 ) THEN

  ! U_ADV
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   udims%i_len,                                                  &
   udims%j_len,                                                  &
   udims_l%halo_i,                                               &
   udims_l%halo_j,                                               &
   udims_l%k_len,                                                &
   fld_type_u,                                                   &
   u_adv,lenrim(fld_type_u,halo_type_extended),                  &
   lbc_size(1,fld_type_u,halo_type_extended),                    &
   lbc_start(1,fld_type_u,halo_type_extended),                   &
   halo_i,halo_j,u_adv_lbc,                                      &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)

  ! V_ADV
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   vdims%i_len,                                                  &
   vdims%j_len,                                                  &
   vdims_l%halo_i,                                               &
   vdims_l%halo_j,                                               &
   vdims_l%k_len,                                                &
   fld_type_v,                                                   &
   v_adv,lenrim(fld_type_v,halo_type_extended),                  &
   lbc_size(1,fld_type_v,halo_type_extended),                    &
   lbc_start(1,fld_type_v,halo_type_extended),                   &
   halo_i,halo_j,v_adv_lbc,                                      &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)

  ! W_ADV
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   wdims%i_len,                                                  &
   wdims%j_len,                                                  &
   wdims_l%halo_i,                                               &
   wdims_l%halo_j,                                               &
   wdims_l%k_len,                                                &
   fld_type_p,                                                   &
   w_adv,lenrim(fld_type_p,halo_type_extended),                  &
   lbc_size(1,fld_type_p,halo_type_extended),                    &
   lbc_start(1,fld_type_p,halo_type_extended),                   &
   halo_i,halo_j,w_adv_lbc,                                      &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)

END IF ! .NOT. L_int_uvw_lbc

! MURK Aerosol
      ! If murk aerosol is in use and murk lbcs are active,
IF (L_murk .AND. L_murk_lbc .AND. rimwidth > 0 ) THEN
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_s%halo_i,                                               &
   tdims_s%halo_j,                                               &
   tdims_s%k_len,                                                &
   fld_type_p,murk,                                              &
   lenrim(fld_type_p,halo_type_single),                          &
   lbc_size(1,fld_type_p,halo_type_single),                      &
   lbc_start(1,fld_type_p,halo_type_single),                     &
   offx,offy,murk_lbc,                                           &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)
END IF

! DUST_DIV1
      ! If dust is in use and dust lbcs are active,
IF (L_dust_div1 .AND. L_dust_div1_lbc .AND. rimwidth > 0 ) THEN
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_s%halo_i,                                               &
   tdims_s%halo_j,                                               &
   tdims_s%k_len,                                                &
   fld_type_p,dust_div1,                                         &
   lenrim(fld_type_p,halo_type_single),                          &
   lbc_size(1,fld_type_p,halo_type_single),                      &
   lbc_start(1,fld_type_p,halo_type_single),                     &
   offx,offy,dust_div1_lbc,                                      &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)
END IF

! DUST_DIV2
      ! If dust is in use and dust lbcs are active,
IF (L_dust_div2 .AND. L_dust_div2_lbc .AND. rimwidth > 0 ) THEN
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_s%halo_i,                                               &
   tdims_s%halo_j,                                               &
   tdims_s%k_len,                                                &
   fld_type_p,dust_div2,                                         &
   lenrim(fld_type_p,halo_type_single),                          &
   lbc_size(1,fld_type_p,halo_type_single),                      &
   lbc_start(1,fld_type_p,halo_type_single),                     &
   offx,offy,dust_div2_lbc,                                      &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)
END IF

! DUST_DIV3
      ! If dust is in use and dust lbcs are active,
IF (L_dust_div3 .AND. L_dust_div3_lbc .AND. rimwidth > 0 ) THEN
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_s%halo_i,                                               &
   tdims_s%halo_j,                                               &
   tdims_s%k_len,                                                &
   fld_type_p,dust_div3,                                         &
   lenrim(fld_type_p,halo_type_single),                          &
   lbc_size(1,fld_type_p,halo_type_single),                      &
   lbc_start(1,fld_type_p,halo_type_single),                     &
   offx,offy,dust_div3_lbc,                                      &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)
END IF

! DUST_DIV4
      ! If dust is in use and dust lbcs are active,
IF (L_dust_div4 .AND. L_dust_div4_lbc .AND. rimwidth > 0 ) THEN
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_s%halo_i,                                               &
   tdims_s%halo_j,                                               &
   tdims_s%k_len,                                                &
   fld_type_p,dust_div4,                                         &
   lenrim(fld_type_p,halo_type_single),                          &
   lbc_size(1,fld_type_p,halo_type_single),                      &
   lbc_start(1,fld_type_p,halo_type_single),                     &
   offx,offy,dust_div4_lbc,                                      &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)
END IF

! DUST_DIV5
      ! If dust is in use and dust lbcs are active,
IF (L_dust_div5 .AND. L_dust_div5_lbc .AND. rimwidth > 0 ) THEN
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_s%halo_i,                                               &
   tdims_s%halo_j,                                               &
   tdims_s%k_len,                                                &
   fld_type_p,dust_div5,                                         &
   lenrim(fld_type_p,halo_type_single),                          &
   lbc_size(1,fld_type_p,halo_type_single),                      &
   lbc_start(1,fld_type_p,halo_type_single),                     &
   offx,offy,dust_div5_lbc,                                      &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)
END IF

! DUST_DIV6
      ! If dust is in use and dust lbcs are active,
IF (L_dust_div6 .AND. L_dust_div6_lbc .AND. rimwidth > 0 ) THEN
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_s%halo_i,                                               &
   tdims_s%halo_j,                                               &
   tdims_s%k_len,                                                &
   fld_type_p,dust_div6,                                         &
   lenrim(fld_type_p,halo_type_single),                          &
   lbc_size(1,fld_type_p,halo_type_single),                      &
   lbc_start(1,fld_type_p,halo_type_single),                     &
   offx,offy,dust_div6_lbc,                                      &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)
END IF

! SO2
      ! If so2 is in use and so2 lbcs are active,
IF (L_so2 .AND. L_so2_lbc .AND. rimwidth > 0 ) THEN
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_s%halo_i,                                               &
   tdims_s%halo_j,                                               &
   tdims_s%k_len,                                                &
   fld_type_p,so2,                                               &
   lenrim(fld_type_p,halo_type_single),                          &
   lbc_size(1,fld_type_p,halo_type_single),                      &
   lbc_start(1,fld_type_p,halo_type_single),                     &
   offx,offy,so2_lbc,                                            &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)
END IF

! DMS
      ! If dms is in use and dms lbcs are active,
IF (L_dms .AND. L_dms_lbc .AND. rimwidth > 0 ) THEN
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_s%halo_i,                                               &
   tdims_s%halo_j,                                               &
   tdims_s%k_len,                                                &
   fld_type_p,dms,                                               &
   lenrim(fld_type_p,halo_type_single),                          &
   lbc_size(1,fld_type_p,halo_type_single),                      &
   lbc_start(1,fld_type_p,halo_type_single),                     &
   offx,offy,dms_lbc,                                            &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)
END IF

! SO4_AITKEN
      ! If so4_aitken is in use and so4_aitken lbcs are active,
IF (L_so4_aitken .AND. L_so4_aitken_lbc .AND. rimwidth > 0) THEN
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_s%halo_i,                                               &
   tdims_s%halo_j,                                               &
   tdims_s%k_len,                                                &
   fld_type_p,so4_aitken,                                        &
   lenrim(fld_type_p,halo_type_single),                          &
   lbc_size(1,fld_type_p,halo_type_single),                      &
   lbc_start(1,fld_type_p,halo_type_single),                     &
   offx,offy,so4_aitken_lbc,                                     &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)
END IF

! SO4_ACCU
      ! If so4_accu is in use and so4_accu lbcs are active,
IF (L_so4_accu .AND. L_so4_accu_lbc .AND. rimwidth > 0 ) THEN
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_s%halo_i,                                               &
   tdims_s%halo_j,                                               &
   tdims_s%k_len,                                                &
   fld_type_p,so4_accu,                                          &
   lenrim(fld_type_p,halo_type_single),                          &
   lbc_size(1,fld_type_p,halo_type_single),                      &
   lbc_start(1,fld_type_p,halo_type_single),                     &
   offx,offy,so4_accu_lbc,                                       &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)
END IF

! SO4_DISS
      ! If so4_diss is in use and so4_diss lbcs are active,
IF (L_so4_diss .AND. L_so4_diss_lbc .AND. rimwidth > 0 ) THEN
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_s%halo_i,                                               &
   tdims_s%halo_j,                                               &
   tdims_s%k_len,                                                &
   fld_type_p,so4_diss,                                          &
   lenrim(fld_type_p,halo_type_single),                          &
   lbc_size(1,fld_type_p,halo_type_single),                      &
   lbc_start(1,fld_type_p,halo_type_single),                     &
   offx,offy,so4_diss_lbc,                                       &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)
END IF

! NH3
      ! If nh3 is in use and nh3 lbcs are active,
IF (L_nh3 .AND. L_nh3_lbc .AND. rimwidth > 0 ) THEN
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_s%halo_i,                                               &
   tdims_s%halo_j,                                               &
   tdims_s%k_len,                                                &
   fld_type_p,nh3,                                               &
   lenrim(fld_type_p,halo_type_single),                          &
   lbc_size(1,fld_type_p,halo_type_single),                      &
   lbc_start(1,fld_type_p,halo_type_single),                     &
   offx,offy,nh3_lbc,                                            &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)
END IF

! SOOT_NEW
      ! If soot_new is in use and soot_new lbcs are active,
IF (L_soot_new .AND. L_soot_new_lbc .AND. rimwidth > 0 ) THEN
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_s%halo_i,                                               &
   tdims_s%halo_j,                                               &
   tdims_s%k_len,                                                &
   fld_type_p,soot_new,                                          &
   lenrim(fld_type_p,halo_type_single),                          &
   lbc_size(1,fld_type_p,halo_type_single),                      &
   lbc_start(1,fld_type_p,halo_type_single),                     &
   offx,offy,soot_new_lbc,                                       &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)
END IF

! SOOT_AGD
      ! If soot_agd is in use and soot_agd lbcs are active,
IF (L_soot_agd .AND. L_soot_agd_lbc .AND. rimwidth > 0 ) THEN
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_s%halo_i,                                               &
   tdims_s%halo_j,                                               &
   tdims_s%k_len,                                                &
   fld_type_p,soot_agd,                                          &
   lenrim(fld_type_p,halo_type_single),                          &
   lbc_size(1,fld_type_p,halo_type_single),                      &
   lbc_start(1,fld_type_p,halo_type_single),                     &
   offx,offy,soot_agd_lbc,                                       &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)
END IF

! SOOT_CLD
      ! If soot_cld is in use and soot_cld lbcs are active,
IF (L_soot_cld .AND. L_soot_cld_lbc .AND. rimwidth > 0 ) THEN
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_s%halo_i,                                               &
   tdims_s%halo_j,                                               &
   tdims_s%k_len,                                                &
   fld_type_p,soot_cld,                                          &
   lenrim(fld_type_p,halo_type_single),                          &
   lbc_size(1,fld_type_p,halo_type_single),                      &
   lbc_start(1,fld_type_p,halo_type_single),                     &
   offx,offy,soot_cld_lbc,                                       &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)
END IF

! BMASS_NEW
      ! If bmass_new is in use and bmass_new lbcs are active,
IF (L_bmass_new .AND. L_bmass_new_lbc .AND. rimwidth > 0 ) THEN
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_s%halo_i,                                               &
   tdims_s%halo_j,                                               &
   tdims_s%k_len,                                                &
   fld_type_p,bmass_new,                                         &
   lenrim(fld_type_p,halo_type_single),                          &
   lbc_size(1,fld_type_p,halo_type_single),                      &
   lbc_start(1,fld_type_p,halo_type_single),                     &
   offx,offy,bmass_new_lbc,                                      &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)
END IF

! BMASS_AGD
      ! If bmass_agd is in use and bmass_agd lbcs are active,
IF (L_bmass_agd .AND. L_bmass_agd_lbc .AND. rimwidth > 0 ) THEN
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_s%halo_i,                                               &
   tdims_s%halo_j,                                               &
   tdims_s%k_len,                                                &
   fld_type_p,bmass_agd,                                         &
   lenrim(fld_type_p,halo_type_single),                          &
   lbc_size(1,fld_type_p,halo_type_single),                      &
   lbc_start(1,fld_type_p,halo_type_single),                     &
   offx,offy,bmass_agd_lbc,                                      &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)
END IF

! BMASS_CLD
      ! If bmass_cld is in use and bmass_cld lbcs are active,
IF (L_bmass_cld .AND. L_bmass_cld_lbc .AND. rimwidth > 0 ) THEN
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_s%halo_i,                                               &
   tdims_s%halo_j,                                               &
   tdims_s%k_len,                                                &
   fld_type_p,bmass_cld,                                         &
   lenrim(fld_type_p,halo_type_single),                          &
   lbc_size(1,fld_type_p,halo_type_single),                      &
   lbc_start(1,fld_type_p,halo_type_single),                     &
   offx,offy,bmass_cld_lbc,                                      &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)
END IF

! OCFF_NEW
      ! If ocff_new is in use and ocff_new lbcs are active,
IF (L_ocff_new .AND. L_ocff_new_lbc .AND. rimwidth > 0 ) THEN
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_s%halo_i,                                               &
   tdims_s%halo_j,                                               &
   tdims_s%k_len,                                                &
   fld_type_p,ocff_new,                                          &
   lenrim(fld_type_p,halo_type_single),                          &
   lbc_size(1,fld_type_p,halo_type_single),                      &
   lbc_start(1,fld_type_p,halo_type_single),                     &
   offx,offy,ocff_new_lbc,                                       &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)
END IF

! OCFF_AGD
      ! If ocff_agd is in use and ocff_agd lbcs are active,
IF (L_ocff_agd .AND. L_ocff_agd_lbc .AND. rimwidth > 0 ) THEN
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_s%halo_i,                                               &
   tdims_s%halo_j,                                               &
   tdims_s%k_len,                                                &
   fld_type_p,ocff_agd,                                          &
   lenrim(fld_type_p,halo_type_single),                          &
   lbc_size(1,fld_type_p,halo_type_single),                      &
   lbc_start(1,fld_type_p,halo_type_single),                     &
   offx,offy,ocff_agd_lbc,                                       &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)
END IF

! OCFF_CLD
      ! If ocff_cld is in use and ocff_cld lbcs are active,
IF (L_ocff_cld .AND. L_ocff_cld_lbc .AND. rimwidth > 0 ) THEN
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_s%halo_i,                                               &
   tdims_s%halo_j,                                               &
   tdims_s%k_len,                                                &
   fld_type_p,ocff_cld,                                          &
   lenrim(fld_type_p,halo_type_single),                          &
   lbc_size(1,fld_type_p,halo_type_single),                      &
   lbc_start(1,fld_type_p,halo_type_single),                     &
   offx,offy,ocff_cld_lbc,                                       &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)
END IF

! NITR_ACC
      ! If nitr_acc is in use and nitr_acc lbcs are active,
IF (L_nitr_acc .AND. L_nitr_acc_lbc .AND. rimwidth > 0 ) THEN
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_s%halo_i,                                               &
   tdims_s%halo_j,                                               &
   tdims_s%k_len,                                                &
   fld_type_p,nitr_acc,                                          &
   lenrim(fld_type_p,halo_type_single),                          &
   lbc_size(1,fld_type_p,halo_type_single),                      &
   lbc_start(1,fld_type_p,halo_type_single),                     &
   offx,offy,nitr_acc_lbc,                                       &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)
END IF

! NITR_DISS
      ! If nitr_diss is in use and nitr_diss lbcs are active,
IF (L_nitr_diss .AND. L_nitr_diss_lbc .AND. rimwidth > 0 ) THEN
  ! DEPENDS ON: set_lateral_boundaries
  CALL set_lateral_boundaries(                                   &
   tdims%i_len,                                                  &
   tdims%j_len,                                                  &
   tdims_s%halo_i,                                               &
   tdims_s%halo_j,                                               &
   tdims_s%k_len,                                                &
   fld_type_p,nitr_diss,                                         &
   lenrim(fld_type_p,halo_type_single),                          &
   lbc_size(1,fld_type_p,halo_type_single),                      &
   lbc_start(1,fld_type_p,halo_type_single),                     &
   offx,offy,nitr_diss_lbc,                                      &
   rimwidth,n_rims_to_do,rimweights,at_extremity,                &
   L_Do_Boundaries,L_Do_Halos)
END IF

! TRACER LBCs
! If tracer lbcs are present, update relevant tracer field
IF (tr_vars > 0 .AND. tr_lbc_vars > 0) THEN

  ! Loop over active tracer prognostics
  DO tracer = 1,tr_vars

    lbc_index = A_tr_active_lbc_index(tracer)

    ! If tracer lbc data is present, update the lbcs
    IF (lbc_index /= -1) THEN

      IF (printstatus >= prstatus_diag) THEN
        WRITE(umMessage,*)                &
            'UPDATE_LAM_LBCS:',tracer,lbc_index
        CALL umPrint(umMessage,src='update_lam_lbcs-updlbc1a')
      END IF

      ! Update tracer field boundaries
      ! DEPENDS ON: set_lateral_boundaries
      IF ( rimwidth > 0 ) THEN
        CALL set_lateral_boundaries(                             &
         tdims%i_len,                                            &
         tdims%j_len,                                            &
         tdims_s%halo_i,                                         &
         tdims_s%halo_j,                                         &
         tdims_s%k_len   ,                                       &
         fld_type_p,                                             &
         free_tracers(tdims_s%i_start:tdims_s%i_end,             &
                      tdims_s%j_start:tdims_s%j_end,             &
                      tdims_s%k_start:tdims_s%k_end,             &
                      tracer),                                   &
         lenrim(fld_type_p,halo_type_extended),                  &
         lbc_size(1,fld_type_p,halo_type_extended),              &
         lbc_start(1,fld_type_p,halo_type_extended),             &
         halo_i,halo_j,tracer_lbc(1,1,lbc_index),                &
         rimwidth,n_rims_to_do,rimweights,at_extremity,          &
         L_Do_Boundaries,L_Do_Halos)
      END IF

    END IF

  END DO

END IF  ! on (TR_VARS > 0 .AND. TR_LBC_VARS > 0)

! UKCA TRACER LBCs
! If ukca tracer lbcs are present, update relevant tracer field
IF (tr_ukca > 0 .AND. tr_lbc_ukca > 0) THEN

  ! Loop over active tracer prognostics
  DO tracer = 1,tr_ukca

    lbc_index = ukca_tr_active_lbc_index(tracer)

    ! If tracer lbc data is present, update the lbcs
    IF (lbc_index /= -1) THEN

      IF (printstatus >= prstatus_diag) THEN
        WRITE(umMessage,*)                &
        'UPDATE_LAM_LBCS (UKCA):',tracer,lbc_index
        CALL umPrint(umMessage,src='update_lam_lbcs-updlbc1a')
      END IF

      ! Update tracer field boundaries
      ! DEPENDS ON: set_lateral_boundaries
      IF ( rimwidth > 0 ) THEN
        CALL set_lateral_boundaries(                             &
         tdims%i_len,                                            &
         tdims%j_len,                                            &
         tdims_s%halo_i,                                         &
         tdims_s%halo_j,                                         &
         tdims_s%k_len   ,                                       &
         fld_type_p,                                             &
         tracer_ukca(tdims_s%i_start:tdims_s%i_end,              &
                     tdims_s%j_start:tdims_s%j_end,              &
                     tdims_s%k_start:tdims_s%k_end,              &
                     tracer),                                    &
         lenrim(fld_type_p,halo_type_extended),                  &
         lbc_size(1,fld_type_p,halo_type_extended),              &
         lbc_start(1,fld_type_p,halo_type_extended),             &
         halo_i,halo_j,tracer_ukca_lbc(1,1,lbc_index),           &
         rimwidth,n_rims_to_do,rimweights,at_extremity,          &
         L_Do_Boundaries,L_Do_Halos)
       END IF

    END IF

  END DO

END IF  ! on (TR_UKCA > 0 .AND. TR_LBC_UKCA > 0)

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE UPDATE_LAM_LBCs
