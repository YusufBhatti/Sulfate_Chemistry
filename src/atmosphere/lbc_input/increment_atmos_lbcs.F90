! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine INCREMENT_ATMOS_LBCS
!
! Purpose : Increments atmosphere LBCs
!
! ---------------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: LBC Input

SUBROUTINE increment_atmos_lbcs(                                  &
  steps_to_next_update,                                           &
  lenrim,                                                         &
  tr_lbc_vars,tr_levels,tr_lbc_ukca,                              &
  L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc, l_pc2_lbc,   &
  L_murk_lbc, L_int_uvw_lbc,                                      &
  L_dust_div1_lbc,L_dust_div2_lbc,                                &
  L_dust_div3_lbc,L_dust_div4_lbc,                                &
  L_dust_div5_lbc,L_dust_div6_lbc,                                &
  L_so2_lbc,L_dms_lbc,L_so4_aitken_lbc,                           &
  L_so4_accu_lbc,L_so4_diss_lbc,                                  &
  L_nh3_lbc,L_soot_new_lbc,L_soot_agd_lbc,                        &
  L_soot_cld_lbc,L_bmass_new_lbc,                                 &
  L_bmass_agd_lbc,L_bmass_cld_lbc,                                &
  L_ocff_new_lbc,L_ocff_agd_lbc,L_ocff_cld_lbc,                   &
  L_nitr_acc_lbc,L_nitr_diss_lbc,                                 &
  u_lbc, u_lbc_tend,                                              &
  v_lbc, v_lbc_tend,                                              &
  w_lbc, w_lbc_tend,                                              &
  rho_lbc, rho_lbc_tend,                                          &
  theta_lbc, theta_lbc_tend,                                      &
  q_lbc, q_lbc_tend,                                              &
  qcl_lbc, qcl_lbc_tend,                                          &
  qcf_lbc, qcf_lbc_tend,                                          &
  qcf2_lbc, qcf2_lbc_tend,                                        &
  qrain_lbc, qrain_lbc_tend,                                      &
  qgraup_lbc, qgraup_lbc_tend,                                    &
  cf_bulk_lbc, cf_bulk_lbc_tend,                                  &
  cf_liquid_lbc, cf_liquid_lbc_tend,                              &
  cf_frozen_lbc, cf_frozen_lbc_tend,                              &
  exner_lbc, exner_lbc_tend,                                      &
  u_adv_lbc, u_adv_lbc_tend,                                      &
  v_adv_lbc, v_adv_lbc_tend,                                      &
  w_adv_lbc, w_adv_lbc_tend,                                      &
  murk_lbc, murk_lbc_tend,                                        &
  dust_div1_lbc,dust_div1_lbc_tend,                               &
  dust_div2_lbc,dust_div2_lbc_tend,                               &
  dust_div3_lbc,dust_div3_lbc_tend,                               &
  dust_div4_lbc,dust_div4_lbc_tend,                               &
  dust_div5_lbc,dust_div5_lbc_tend,                               &
  dust_div6_lbc,dust_div6_lbc_tend,                               &
  so2_lbc,so2_lbc_tend,                                           &
  dms_lbc,dms_lbc_tend,                                           &
  so4_aitken_lbc,so4_aitken_lbc_tend,                             &
  so4_accu_lbc,so4_accu_lbc_tend,                                 &
  so4_diss_lbc,so4_diss_lbc_tend,                                 &
  nh3_lbc,nh3_lbc_tend,                                           &
  soot_new_lbc,soot_new_lbc_tend,                                 &
  soot_agd_lbc,soot_agd_lbc_tend,                                 &
  soot_cld_lbc,soot_cld_lbc_tend,                                 &
  bmass_new_lbc,bmass_new_lbc_tend,                               &
  bmass_agd_lbc,bmass_agd_lbc_tend,                               &
  bmass_cld_lbc,bmass_cld_lbc_tend,                               &
  ocff_new_lbc,ocff_new_lbc_tend,                                 &
  ocff_agd_lbc,ocff_agd_lbc_tend,                                 &
  ocff_cld_lbc,ocff_cld_lbc_tend,                                 &
  nitr_acc_lbc,nitr_acc_lbc_tend,                                 &
  nitr_diss_lbc,nitr_diss_lbc_tend,                               &
  tracer_lbc, tracer_lbc_tend,                                    &
  tracer_ukca_lbc, tracer_ukca_lbc_tend,                          &
  rim_stepsa,icode)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE atm_fields_bounds_mod, ONLY: udims_s, vdims_s, wdims_s,       &
                                 pdims_s, tdims_s, tdims_l
USE nlsizes_namelist_mod, ONLY: model_levels
USE ereport_mod, ONLY: ereport
USE UM_ParParams, ONLY: NHalo_max, halo_type_single, nfld_max,    &
                        halo_type_extended, fld_type_u,           &
                        fld_type_v, fld_type_w, fld_type_p
USE umPrintMgr, ONLY: umPrint, umMessage
IMPLICIT NONE


! Parameters required for argument declarations

! Arguments:

INTEGER ::                                                        &
  steps_to_next_update                                            &
                              ! IN : Number of timesteps before
                              !      next LBC update
, lenrim(Nfld_max,NHalo_Max)                                      &
                              ! IN : Size of a level of LBC
, tr_lbc_vars                                                     &
                              ! IN : Number of tracer LBCs
, tr_lbc_ukca                                                     &
                              ! IN : Number of UKCA tracer LBCs
, tr_levels                                                       &
                              ! IN : Number of tracer levels
, rim_stepsa                  ! IN : Number of timesteps in LBC
                              !      period

LOGICAL, INTENT (IN) ::                                           &
  L_mcr_qcf2_lbc                                                  &
                   ! true if using second cloud ice lbcs
, L_mcr_qrain_lbc                                                 &
                   ! true if using rain lbcs
, L_mcr_qgraup_lbc                                                &
                   ! true if using graupel lbcs
, L_pc2_lbc                                                       &
                   ! true if using cloud fraction lbcs
, L_int_uvw_lbc                                                   &
                    ! true if using interpolated advecting winds
                    !  in lateral boundaries
, L_murk_lbc                                                      &
                    ! true if using murk aerosol lbcs
, L_dust_div1_lbc                                                 &
, L_dust_div2_lbc                                                 &
, L_dust_div3_lbc                                                 &
, L_dust_div4_lbc                                                 &
, L_dust_div5_lbc                                                 &
, L_dust_div6_lbc                                                 &
                   ! true if DUST in lbcs
, L_so2_lbc                                                       &
                   ! true if SO2 in lbcs
, L_dms_lbc                                                       &
                   ! true if DMS in lbcs
, L_so4_aitken_lbc                                                &
                   ! true if SO4_AITKEN in lbcs
, L_so4_accu_lbc                                                  &
                   ! true if SO4_ACCU in lbcs
, L_so4_diss_lbc                                                  &
                   ! true if SO4_DISS in lbcs
, L_nh3_lbc                                                       &
                   ! true if NH3 in lbcs
, L_soot_new_lbc                                                  &
, L_soot_agd_lbc                                                  &
, L_soot_cld_lbc                                                  &
                   ! true if soot in lbcs
, L_bmass_new_lbc                                                 &
, L_bmass_agd_lbc                                                 &
, L_bmass_cld_lbc                                                 &
                   ! true if biomass in lbcs
, L_ocff_new_lbc                                                  &
, L_ocff_agd_lbc                                                  &
, L_ocff_cld_lbc                                                  &
                   ! true if fossil fuel aerosol in lbcs
, L_nitr_acc_lbc                                                  &
, L_nitr_diss_lbc   ! true if nitrate aerosol in lbcs

! Note: LBCs are the current value of the LBC (and will be updated)
!       LBC_TENDs are the value of the LBC at the end of the LBC
!       period, towards which the LBC will tend.

REAL, INTENT (INOUT) ::                                           &
  u_lbc(lenrim(fld_type_u,halo_type_extended),                    &
        udims_s%k_start:udims_s%k_end)                            &
                              ! IN/OUT : U LBC
, u_lbc_tend(lenrim(fld_type_u,halo_type_extended),               &
                    udims_s%k_start:udims_s%k_end)                &
                              ! IN : U LBC tendency
, v_lbc(lenrim(fld_type_v,halo_type_extended),                    &
        vdims_s%k_start:vdims_s%k_end)                            &
                              ! IN/OUT : V LBC
, v_lbc_tend(lenrim(fld_type_v,halo_type_extended),               &
             vdims_s%k_start:vdims_s%k_end)                       &
                              ! IN : V LBC tendency
, w_lbc(lenrim(fld_type_p,halo_type_extended),                    &
        wdims_s%k_start:wdims_s%k_end)                            &
                              ! IN/OUT : V LBC
, w_lbc_tend(lenrim(fld_type_p,halo_type_extended),               &
             wdims_s%k_start:wdims_s%k_end)                       &
                              ! IN : V LBC tendency
, rho_lbc(lenrim(fld_type_p,halo_type_extended),                  &
          pdims_s%k_start:pdims_s%k_end)                          &
                              ! IN/OUT : rho LBC
, rho_lbc_tend(lenrim(fld_type_p,halo_type_extended),             &
               pdims_s%k_start:pdims_s%k_end)                     &
                              ! IN : rho LBC tendency
, theta_lbc(lenrim(fld_type_p,halo_type_extended),                &
            tdims_s%k_start:tdims_s%k_end)                        &
                              ! IN/OUT : theta LBC
, theta_lbc_tend(lenrim(fld_type_p,halo_type_extended),           &
                 tdims_s%k_start:tdims_s%k_end)                   &
                              ! IN : theta LBC tendency
, q_lbc(lenrim(fld_type_p,halo_type_extended),                    &
        tdims_l%k_start:tdims_l%k_end)                            &
                              ! IN/OUT : Q LBC
, q_lbc_tend(lenrim(fld_type_p,halo_type_extended),               &
             tdims_l%k_start:tdims_l%k_end)                       &
                              ! IN : Q LBC tendency
, qcl_lbc(lenrim(fld_type_p,halo_type_extended),                  &
          tdims_l%k_start:tdims_l%k_end)                          &
                              ! IN/OUT : QCL LBC
, qcl_lbc_tend(lenrim(fld_type_p,halo_type_extended),             &
               tdims_l%k_start:tdims_l%k_end)                     &
                              ! IN : QCL LBC tendency
, qcf_lbc(lenrim(fld_type_p,halo_type_extended),                  &
          tdims_l%k_start:tdims_l%k_end)                          &
                              ! IN/OUT : QCL LBC
, qcf_lbc_tend(lenrim(fld_type_p,halo_type_extended),             &
               tdims_l%k_start:tdims_l%k_end)                     &
                              ! IN : QCL LBC tendency
, qcf2_lbc(lenrim(fld_type_p,halo_type_extended),                 &
           tdims_l%k_start:tdims_l%k_end)                         &
                              ! IN/OUT : QCF2 LBC
, qcf2_lbc_tend(lenrim(fld_type_p,halo_type_extended),            &
                tdims_l%k_start:tdims_l%k_end)                    &
                              ! IN : QCF2 LBC tendency
, qrain_lbc(lenrim(fld_type_p,halo_type_extended),                &
            tdims_l%k_start:tdims_l%k_end)                        &
                              ! IN/OUT : QRAIN LBC
, qrain_lbc_tend(lenrim(fld_type_p,halo_type_extended),           &
                 tdims_l%k_start:tdims_l%k_end)                   &
                              ! IN : QRAIN LBC tendency
, qgraup_lbc(lenrim(fld_type_p,halo_type_extended),               &
             tdims_l%k_start:tdims_l%k_end)                       &
                              ! IN/OUT : QGRAUP LBC
, qgraup_lbc_tend(lenrim(fld_type_p,halo_type_extended),          &
                  tdims_l%k_start:tdims_l%k_end)                  &
                              ! IN : QGRAUP LBC tendency
, cf_bulk_lbc(lenrim(fld_type_p,halo_type_extended),              &
              tdims_l%k_start:tdims_l%k_end)                      &
                              ! IN/OUT : CF_BULK LBC
, cf_bulk_lbc_tend(lenrim(fld_type_p,halo_type_extended),         &
                   tdims_l%k_start:tdims_l%k_end)                 &
                              ! IN : CF_BULK LBC tendency
, cf_liquid_lbc(lenrim(fld_type_p,halo_type_extended),            &
                tdims_l%k_start:tdims_l%k_end)                    &
                              ! IN/OUT : CF_LIQUID LBC
, cf_liquid_lbc_tend(lenrim(fld_type_p,halo_type_extended),       &
                            tdims_l%k_start:tdims_l%k_end)        &
                              ! IN : CF_LIQUID LBC tendency
, cf_frozen_lbc(lenrim(fld_type_p,halo_type_extended),            &
                tdims_l%k_start:tdims_l%k_end)                    &
                              ! IN/OUT : CF_FROZEN LBC
, cf_frozen_lbc_tend(lenrim(fld_type_p,halo_type_extended),       &
                     tdims_l%k_start:tdims_l%k_end)               &
                              ! IN : CF_FROZEN LBC tendency
, exner_lbc(lenrim(fld_type_p,halo_type_extended),                &
            pdims_s%k_start:pdims_s%k_end+1)                      &
                              ! IN/OUT : Exner LBC
, exner_lbc_tend(lenrim(fld_type_p,halo_type_extended),           &
                 pdims_s%k_start:pdims_s%k_end+1)                 &
                                   ! IN : Exner LBC tendency
, u_adv_lbc(lenrim(fld_type_u,halo_type_extended),                &
            udims_s%k_start:udims_s%k_end)                        &
                              ! IN/OUT : u_adv LBC
, u_adv_lbc_tend(lenrim(fld_type_u,halo_type_extended),           &
                 udims_s%k_start:udims_s%k_end)                   &
                              ! IN : u_adv LBC tendency
, v_adv_lbc(lenrim(fld_type_v,halo_type_extended),                &
            vdims_s%k_start:vdims_s%k_end)                        &
                              ! IN/OUT : v_adv LBC
, v_adv_lbc_tend(lenrim(fld_type_v,halo_type_extended),           &
                 vdims_s%k_start:vdims_s%k_end)                   &
                              ! IN : v_adv LBC tendency
, w_adv_lbc(lenrim(fld_type_p,halo_type_extended),                &
                   wdims_s%k_start:wdims_s%k_end)                 &
                              ! IN/OUT : W LBC
, w_adv_lbc_tend(lenrim(fld_type_p,halo_type_extended),           &
                 wdims_s%k_start:wdims_s%k_end)                   &
                                 ! IN : W LBC tendency
, murk_lbc(lenrim(fld_type_p,halo_type_single),                   &
           tdims_s%k_start:tdims_s%k_end)                         &
                              ! IN/OUT : MURK LBC
, murk_lbc_tend(lenrim(fld_type_p,halo_type_single),              &
                tdims_s%k_start:tdims_s%k_end) ! IN : MURK LBC tendency

REAL, INTENT (INOUT) ::                                           &
  dust_div1_lbc(lenrim(fld_type_p,halo_type_single),              &
                tdims_s%k_start:tdims_s%k_end)                    &
                              ! IN/OUT : dust_div1 LBC
, dust_div1_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                     tdims_s%k_start:tdims_s%k_end)               &
                              ! IN : dust_div1 LBC tendency
, dust_div2_lbc(lenrim(fld_type_p,halo_type_single),              &
                tdims_s%k_start:tdims_s%k_end)                    &
                              ! IN/OUT : dust_div2 LBC
, dust_div2_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                     tdims_s%k_start:tdims_s%k_end)               &
                              ! IN : dust_div2 LBC tendency
, dust_div3_lbc(lenrim(fld_type_p,halo_type_single),              &
                tdims_s%k_start:tdims_s%k_end)                    &
                              ! IN/OUT : dust_div3 LBC
, dust_div3_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                     tdims_s%k_start:tdims_s%k_end)               &
                              ! IN : dust_div3 LBC tendency
, dust_div4_lbc(lenrim(fld_type_p,halo_type_single),              &
                tdims_s%k_start:tdims_s%k_end)                    &
                              ! IN/OUT : dust_div4 LBC
, dust_div4_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                     tdims_s%k_start:tdims_s%k_end)               &
                              ! IN : dust_div4 LBC tendency
, dust_div5_lbc(lenrim(fld_type_p,halo_type_single),              &
                tdims_s%k_start:tdims_s%k_end)                    &
                              ! IN/OUT : dust_div5 LBC
, dust_div5_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                     tdims_s%k_start:tdims_s%k_end)               &
                              ! IN : dust_div5 LBC tendency
, dust_div6_lbc(lenrim(fld_type_p,halo_type_single),              &
                tdims_s%k_start:tdims_s%k_end)                    &
                              ! IN/OUT : dust_div6 LBC
, dust_div6_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                     tdims_s%k_start:tdims_s%k_end)               &
                              ! IN : dust_div6 LBC tendency
, SO2_lbc(lenrim(fld_type_p,halo_type_single),                    &
          tdims_s%k_start:tdims_s%k_end)                          &
                              ! IN/OUT : SO2 LBC
, SO2_lbc_tend(lenrim(fld_type_p,halo_type_single),               &
               tdims_s%k_start:tdims_s%k_end)                     &
                              ! IN : SO2 LBC tendency
, dms_lbc(lenrim(fld_type_p,halo_type_single),                    &
          tdims_s%k_start:tdims_s%k_end)                          &
                              ! IN/OUT : DMS LBC
, dms_lbc_tend(lenrim(fld_type_p,halo_type_single),               &
               tdims_s%k_start:tdims_s%k_end)                     &
                              ! IN : DMS LBC tendency
, SO4_aitken_lbc(lenrim(fld_type_p,halo_type_single),             &
                 tdims_s%k_start:tdims_s%k_end)                   &
                              ! IN/OUT : SO4_AITKEN LBC
, SO4_aitken_lbc_tend(lenrim(fld_type_p,halo_type_single),        &
                      tdims_s%k_start:tdims_s%k_end)              &
                              ! IN : SO4_AITKEN LBC tendency
, SO4_accu_lbc(lenrim(fld_type_p,halo_type_single),               &
               tdims_s%k_start:tdims_s%k_end)                     &
                              ! IN/OUT : SO4_accU LBC
, SO4_accu_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                    tdims_s%k_start:tdims_s%k_end)                &
                              ! IN : SO4_accU LBC tendency
, SO4_diss_lbc(lenrim(fld_type_p,halo_type_single),               &
               tdims_s%k_start:tdims_s%k_end)                     &
                              ! IN/OUT : SO4_DISS LBC
, SO4_diss_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                    tdims_s%k_start:tdims_s%k_end)                &
                              ! IN : SO4_DISS LBC tendency
, NH3_lbc(lenrim(fld_type_p,halo_type_single),                    &
          tdims_s%k_start:tdims_s%k_end)                          &
                              ! IN/OUT : NH3 LBC
, NH3_lbc_tend(lenrim(fld_type_p,halo_type_single),               &
               tdims_s%k_start:tdims_s%k_end)                     &
                              ! IN : NH3 LBC tendency
, soot_new_lbc(lenrim(fld_type_p,halo_type_single),               &
               tdims_s%k_start:tdims_s%k_end)                     &
                              ! IN/OUT : soot_NEW LBC
, soot_new_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                    tdims_s%k_start:tdims_s%k_end)                &
                              ! IN : soot_NEW LBC tendency
, soot_agd_lbc(lenrim(fld_type_p,halo_type_single),               &
               tdims_s%k_start:tdims_s%k_end)                     &
                              ! IN/OUT : soot_agd LBC
, soot_agd_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                    tdims_s%k_start:tdims_s%k_end)                &
                              ! IN : soot_agd LBC tendency
, soot_cld_lbc(lenrim(fld_type_p,halo_type_single),               &
               tdims_s%k_start:tdims_s%k_end)                     &
                              ! IN/OUT : soot_CLD LBC
, soot_cld_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                    tdims_s%k_start:tdims_s%k_end)                &
                              ! IN : soot_CLD LBC tendency
, bmass_new_lbc(lenrim(fld_type_p,halo_type_single),              &
                tdims_s%k_start:tdims_s%k_end)                    &
                              ! IN/OUT : bmass_NEW LBC
, bmass_new_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                     tdims_s%k_start:tdims_s%k_end)               &
                              ! IN : bmass_NEW LBC tendency
, bmass_agd_lbc(lenrim(fld_type_p,halo_type_single),              &
                tdims_s%k_start:tdims_s%k_end)                    &
                              ! IN/OUT : bmass_agd LBC
, bmass_agd_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                     tdims_s%k_start:tdims_s%k_end)               &
                              ! IN : bmass_agd LBC tendency
, bmass_cld_lbc(lenrim(fld_type_p,halo_type_single),              &
                tdims_s%k_start:tdims_s%k_end)                    &
                              ! IN/OUT : bmass_CLD LBC
, bmass_cld_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                     tdims_s%k_start:tdims_s%k_end)               &
                              ! IN : bmass_CLD LBC tendency
, ocff_new_lbc(lenrim(fld_type_p,halo_type_single),               &
               tdims_s%k_start:tdims_s%k_end)                     &
                              ! IN/OUT : ocff_NEW LBC
, ocff_new_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                    tdims_s%k_start:tdims_s%k_end)                &
                              ! IN : ocff_NEW LBC tendency
, ocff_agd_lbc(lenrim(fld_type_p,halo_type_single),               &
               tdims_s%k_start:tdims_s%k_end)                     &
                              ! IN/OUT : ocff_agd LBC
, ocff_agd_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                    tdims_s%k_start:tdims_s%k_end)                &
                              ! IN : ocff_agd LBC tendency
, ocff_cld_lbc(lenrim(fld_type_p,halo_type_single),               &
               tdims_s%k_start:tdims_s%k_end)                     &
                              ! IN/OUT : ocff_CLD LBC
, ocff_cld_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                    tdims_s%k_start:tdims_s%k_end)                &
                              ! IN : ocff_CLD LBC tendency
, nitr_acc_lbc(lenrim(fld_type_p,halo_type_single),               &
               tdims_s%k_start:tdims_s%k_end)                     &
                              ! IN/OUT : nitr_acc LBC
, nitr_acc_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                    tdims_s%k_start:tdims_s%k_end)                &
                              ! IN : nitr_acc LBC tendency
, nitr_diss_lbc(lenrim(fld_type_p,halo_type_single),              &
                tdims_s%k_start:tdims_s%k_end)                    &
                              ! IN/OUT : nitr_DISS LBC
, nitr_diss_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                     tdims_s%k_start:tdims_s%k_end)               &
                              ! IN : nitr_DISS LBC tendency
, tracer_lbc(lenrim(fld_type_p,halo_type_extended),               &
                    tdims_l%k_start:tdims_l%k_end,tr_lbc_vars)    &
                                     ! IN/OUT : Tracer LBCs
, tracer_lbc_tend(lenrim(fld_type_p,halo_type_extended),          &
                  tdims_l%k_start:tdims_l%k_end,tr_lbc_vars)      &
                                     ! IN : Tracer LBCs tendency
, tracer_ukca_lbc(lenrim(fld_type_p,halo_type_extended),          &
                  tdims_l%k_start:tdims_l%k_end,tr_lbc_ukca)      &
                                     ! IN/OUT : UKCA Tracer LBCs
, tracer_ukca_lbc_tend(lenrim(fld_type_p,halo_type_extended),     &
                       tdims_l%k_start:tdims_l%k_end,tr_lbc_ukca)
                                     ! IN : UKCA Tracer LBCs tend


INTEGER ::                                                        &
  icode                           ! Error code

! Local variables

REAL ::                                                           &
  increment_factor   ! Factor to increment LBC by

INTEGER ::                                                        &
  tracer                                                          &
            ! loop counter over tracer variables
, k                                                               &
            ! loop counter over levels
, i         ! loop counter over horizontal dimension

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INCREMENT_ATMOS_LBCS'

! ---------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! Calculate increment_factor

IF (steps_to_next_update  >   rim_stepsa) THEN
  WRITE(umMessage,*) 'STEPS_TO_NEXT_UPDATE must be < ',rim_stepsa
  CALL umPrint(umMessage,src='increment_atmos_lbcs')
  WRITE(umMessage,*) 'Value supplied was ',steps_to_next_update
  CALL umPrint(umMessage,src='increment_atmos_lbcs')
  icode=1

  CALL Ereport("INCREMENT_ATMOS_LBCS ", icode,                    &
               "STEPS_TO_NEXT_UPDATE>BCSTEPS" )
END IF

increment_factor=1.0/steps_to_next_update

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k,i,tracer)

!$OMP DO SCHEDULE(STATIC)
DO k=udims_s%k_start,udims_s%k_end
  DO i=1,lenrim(fld_type_u,halo_type_extended)
    u_lbc(i,k)=u_lbc(i,k)+increment_factor*                       &
                          (u_lbc_tend(i,k)-u_lbc(i,k))
  END DO ! i
END DO !k=udims_s%k_start,udims_s%k_end
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO k=vdims_s%k_start,vdims_s%k_end
  DO i=1,lenrim(fld_type_v,halo_type_extended)
    v_lbc(i,k)=v_lbc(i,k)+increment_factor*                       &
                          (v_lbc_tend(i,k)-v_lbc(i,k))
  END DO ! i
END DO !k=vdims_s%k_start,vdims_s%k_end
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO k=pdims_s%k_start,pdims_s%k_end
  DO i=1,lenrim(fld_type_p,halo_type_extended)
    rho_lbc(i,k)=rho_lbc(i,k)+increment_factor*                   &
                              (rho_lbc_tend(i,k)-rho_lbc(i,k))
  END DO ! i
END DO ! k=pdims_s%k_start,pdims_s%k_end
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO k=tdims_s%k_start,tdims_s%k_end
  DO i=1,lenrim(fld_type_p,halo_type_extended)
    theta_lbc(i,k)=theta_lbc(i,k)+increment_factor*               &
                                  (theta_lbc_tend(i,k)-           &
                                   theta_lbc(i,k))
  END DO ! i
END DO ! k=tdims_s%k_start,tdims_s%k_end
!$OMP END DO NOWAIT


!$OMP DO SCHEDULE(STATIC)
DO k=tdims_l%k_start,tdims_l%k_end
  DO i=1,lenrim(fld_type_p,halo_type_extended)
    q_lbc(i,k)=q_lbc(i,k)+increment_factor*                       &
                          (q_lbc_tend(i,k)-q_lbc(i,k))
    qcl_lbc(i,k)=qcl_lbc(i,k)+increment_factor*                   &
                              (qcl_lbc_tend(i,k)-qcl_lbc(i,k))
    qcf_lbc(i,k)=qcf_lbc(i,k)+increment_factor*                   &
                              (qcf_lbc_tend(i,k)-qcf_lbc(i,k))
  END DO ! i
END DO ! k
!$OMP END DO NOWAIT

IF (L_mcr_qcf2_lbc) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k=tdims_l%k_start,tdims_l%k_end
    DO i=1,lenrim(fld_type_p,halo_type_extended)
      qcf2_lbc(i,k) = qcf2_lbc(i,k) + increment_factor *          &
                       (qcf2_lbc_tend(i,k)-qcf2_lbc(i,k))
    END DO ! i
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (L_mcr_qrain_lbc) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k=tdims_l%k_start,tdims_l%k_end
    DO i=1,lenrim(fld_type_p,halo_type_extended)
      qrain_lbc(i,k) = qrain_lbc(i,k) + increment_factor *        &
                       (qrain_lbc_tend(i,k)-qrain_lbc(i,k))
    END DO ! i
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (L_mcr_qgraup_lbc) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k=tdims_l%k_start,tdims_l%k_end
    DO i=1,lenrim(fld_type_p,halo_type_extended)
      qgraup_lbc(i,k) = qgraup_lbc(i,k) + increment_factor *      &
                         (qgraup_lbc_tend(i,k)-qgraup_lbc(i,k))
    END DO ! i
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (L_pc2_lbc) THEN  ! Cloud fractions are in lbcs
!$OMP DO SCHEDULE(STATIC)
  DO k=tdims_l%k_start,tdims_l%k_end
    DO i=1,lenrim(fld_type_p,halo_type_extended)
      cf_bulk_lbc(i,k) = cf_bulk_lbc(i,k) + increment_factor *    &
                     (cf_bulk_lbc_tend(i,k)-cf_bulk_lbc(i,k))
      cf_liquid_lbc(i,k) = cf_liquid_lbc(i,k)+increment_factor*   &
                     (cf_liquid_lbc_tend(i,k)-cf_liquid_lbc(i,k))
      cf_frozen_lbc(i,k) = cf_frozen_lbc(i,k)+increment_factor*   &
                     (cf_frozen_lbc_tend(i,k)-cf_frozen_lbc(i,k))
    END DO ! i
  END DO ! k
!$OMP END DO NOWAIT
END IF ! L_pc2_lbc

!$OMP DO SCHEDULE(STATIC)
DO k=wdims_s%k_start,wdims_s%k_end
  DO i=1,lenrim(fld_type_p,halo_type_extended)
    w_lbc(i,k)=w_lbc(i,k)+increment_factor*                       &
                          (w_lbc_tend(i,k)-w_lbc(i,k))
  END DO ! i
END DO ! k
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO k=pdims_s%k_start,pdims_s%k_end+1
  DO i=1,lenrim(fld_type_p,halo_type_extended)
    exner_lbc(i,k)=exner_lbc(i,k)+increment_factor*              &
                                 (exner_lbc_tend(i,k)-           &
                                  exner_lbc(i,k))
  END DO ! i
END DO ! k
!$OMP END DO NOWAIT

IF ( .NOT. L_int_uvw_lbc) THEN

!$OMP DO SCHEDULE(STATIC)
  DO k = udims_s%k_start, udims_s%k_end
    DO i = 1, lenrim(fld_type_u,halo_type_extended)
      u_adv_lbc(i,k)=u_adv_lbc(i,k)+increment_factor *            &
                                  (u_adv_lbc_tend(i,k) -          &
                                   u_adv_lbc(i,k))
    END DO ! i
  END DO !k = udims_s%k_start, udims_s%k_end
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO k = vdims_s%k_start, vdims_s%k_end
    DO i = 1, lenrim(fld_type_v,halo_type_extended)
      v_adv_lbc(i,k)=v_adv_lbc(i,k)+increment_factor *            &
                                  (v_adv_lbc_tend(i,k) -          &
                                   v_adv_lbc(i,k))
    END DO ! i
  END DO !  k = vdims_s%k_start, vdims_s%k_end
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO k = wdims_s%k_start,wdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_extended)
      w_adv_lbc(i,k)=w_adv_lbc(i,k)+increment_factor *            &
                                  (w_adv_lbc_tend(i,k) -          &
                                   w_adv_lbc(i,k))
    END DO ! i
  END DO ! k = wdims_s%k_start,wdims_s%k_end
!$OMP END DO NOWAIT

END IF !  .NOT. L_int_uvw_lbc

IF (L_murk_lbc) THEN
  ! murk aerosol turned on and murk lbcs active
!$OMP DO SCHEDULE(STATIC)
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      murk_lbc(i,k) = murk_lbc(i,k) + increment_factor *          &
                       (murk_lbc_tend(i,k)-murk_lbc(i,k))
    END DO ! i
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (L_dust_div1_lbc) THEN
  ! dust lbcs active
!$OMP DO SCHEDULE(STATIC)
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      dust_div1_lbc(i,k) = dust_div1_lbc(i,k) + increment_factor * &
                  (dust_div1_lbc_tend(i,k)-dust_div1_lbc(i,k))
    END DO ! i
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (L_dust_div2_lbc) THEN
  ! dust lbcs active
!$OMP DO SCHEDULE(STATIC)
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      dust_div2_lbc(i,k) = dust_div2_lbc(i,k) + increment_factor * &
                  (dust_div2_lbc_tend(i,k)-dust_div2_lbc(i,k))
    END DO ! i
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (L_dust_div3_lbc) THEN
  ! dust lbcs active
!$OMP DO SCHEDULE(STATIC)
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      dust_div3_lbc(i,k) = dust_div3_lbc(i,k) + increment_factor * &
                  (dust_div3_lbc_tend(i,k)-dust_div3_lbc(i,k))
    END DO ! i
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (L_dust_div4_lbc) THEN
  ! dust lbcs active
!$OMP DO SCHEDULE(STATIC)
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      dust_div4_lbc(i,k) = dust_div4_lbc(i,k) + increment_factor * &
                  (dust_div4_lbc_tend(i,k)-dust_div4_lbc(i,k))
    END DO ! i
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (L_dust_div5_lbc) THEN
  ! dust lbcs active
!$OMP DO SCHEDULE(STATIC)
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      dust_div5_lbc(i,k) = dust_div5_lbc(i,k) + increment_factor * &
                  (dust_div5_lbc_tend(i,k)-dust_div5_lbc(i,k))
    END DO ! i
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (L_dust_div6_lbc) THEN
  ! dust lbcs active
!$OMP DO SCHEDULE(STATIC)
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      dust_div6_lbc(i,k) = dust_div6_lbc(i,k) + increment_factor * &
                  (dust_div6_lbc_tend(i,k)-dust_div6_lbc(i,k))
    END DO ! i
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (L_so2_lbc) THEN
  ! so2 lbcs active
!$OMP DO SCHEDULE(STATIC)
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      so2_lbc(i,k) = so2_lbc(i,k) + increment_factor *            &
                             (so2_lbc_tend(i,k)-so2_lbc(i,k))
    END DO ! i
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (L_dms_lbc) THEN
  ! dms lbcs active
!$OMP DO SCHEDULE(STATIC)
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      dms_lbc(i,k) = dms_lbc(i,k) + increment_factor *            &
                             (dms_lbc_tend(i,k)-dms_lbc(i,k))
    END DO ! i
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (L_so4_aitken_lbc) THEN
  ! so4_aitken lbcs active
!$OMP DO SCHEDULE(STATIC)
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      so4_aitken_lbc(i,k) = so4_aitken_lbc(i,k) + increment_factor * &
                  (so4_aitken_lbc_tend(i,k)-so4_aitken_lbc(i,k))
    END DO ! i
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (L_so4_accu_lbc) THEN
  ! so4_accu lbcs active
!$OMP DO SCHEDULE(STATIC)
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      so4_accu_lbc(i,k) = so4_accu_lbc(i,k) + increment_factor *  &
               (so4_accu_lbc_tend(i,k)-so4_accu_lbc(i,k))
    END DO ! i
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (L_so4_diss_lbc) THEN
  ! so4_diss lbcs active
!$OMP DO SCHEDULE(STATIC)
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      so4_diss_lbc(i,k) = so4_diss_lbc(i,k) + increment_factor *  &
               (so4_diss_lbc_tend(i,k)-so4_diss_lbc(i,k))
    END DO ! i
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (L_nh3_lbc) THEN
  ! nh3 lbcs active
!$OMP DO SCHEDULE(STATIC)
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      nh3_lbc(i,k) = nh3_lbc(i,k) + increment_factor *            &
               (nh3_lbc_tend(i,k)-nh3_lbc(i,k))
    END DO ! i
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (L_soot_new_lbc) THEN
  ! soot_new lbcs active
!$OMP DO SCHEDULE(STATIC)
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      soot_new_lbc(i,k) = soot_new_lbc(i,k) + increment_factor *  &
               (soot_new_lbc_tend(i,k)-soot_new_lbc(i,k))
    END DO ! i
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (L_soot_agd_lbc) THEN
  ! soot_agd lbcs active
!$OMP DO SCHEDULE(STATIC)
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      soot_agd_lbc(i,k) = soot_agd_lbc(i,k) + increment_factor *  &
               (soot_agd_lbc_tend(i,k)-soot_agd_lbc(i,k))
    END DO ! i
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (L_soot_cld_lbc) THEN
  ! soot_cld lbcs active
!$OMP DO SCHEDULE(STATIC)
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      soot_cld_lbc(i,k) = soot_cld_lbc(i,k) + increment_factor *  &
               (soot_cld_lbc_tend(i,k)-soot_cld_lbc(i,k))
    END DO ! i
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (L_bmass_new_lbc) THEN
  ! bmass_new lbcs active
!$OMP DO SCHEDULE(STATIC)
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      bmass_new_lbc(i,k) = bmass_new_lbc(i,k) + increment_factor * &
                (bmass_new_lbc_tend(i,k)-bmass_new_lbc(i,k))
    END DO ! i
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (L_bmass_agd_lbc) THEN
  ! bmass_agd lbcs active
!$OMP DO SCHEDULE(STATIC)
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      bmass_agd_lbc(i,k) = bmass_agd_lbc(i,k) + increment_factor * &
                (bmass_agd_lbc_tend(i,k)-bmass_agd_lbc(i,k))
    END DO ! i
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (L_bmass_cld_lbc) THEN
  ! bmass_cld lbcs active
!$OMP DO SCHEDULE(STATIC)
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      bmass_cld_lbc(i,k) = bmass_cld_lbc(i,k) + increment_factor * &
                (bmass_cld_lbc_tend(i,k)-bmass_cld_lbc(i,k))
    END DO ! i
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (L_ocff_new_lbc) THEN
  ! ocff_new lbcs active
!$OMP DO SCHEDULE(STATIC)
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      ocff_new_lbc(i,k) = ocff_new_lbc(i,k) + increment_factor *  &
               (ocff_new_lbc_tend(i,k)-ocff_new_lbc(i,k))
    END DO ! i
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (L_ocff_agd_lbc) THEN
  ! ocff_agd lbcs active
!$OMP DO SCHEDULE(STATIC)
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      ocff_agd_lbc(i,k) = ocff_agd_lbc(i,k) + increment_factor *  &
               (ocff_agd_lbc_tend(i,k)-ocff_agd_lbc(i,k))
    END DO ! i
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (L_ocff_cld_lbc) THEN
  ! ocff_cld lbcs active
!$OMP DO SCHEDULE(STATIC)
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      ocff_cld_lbc(i,k) = ocff_cld_lbc(i,k) + increment_factor *  &
               (ocff_cld_lbc_tend(i,k)-ocff_cld_lbc(i,k))
    END DO ! i
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (L_nitr_acc_lbc) THEN
  ! nitr_new lbcs active
!$OMP DO SCHEDULE(STATIC)
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      nitr_acc_lbc(i,k) = nitr_acc_lbc(i,k) + increment_factor *  &
               (nitr_acc_lbc_tend(i,k)-nitr_acc_lbc(i,k))
    END DO ! i
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (L_nitr_diss_lbc) THEN
  ! nitr_diss lbcs active
!$OMP DO SCHEDULE(STATIC)
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      nitr_diss_lbc(i,k) = nitr_diss_lbc(i,k) + increment_factor *  &
               (nitr_diss_lbc_tend(i,k)-nitr_diss_lbc(i,k))
    END DO ! i
  END DO ! k
!$OMP END DO NOWAIT
END IF

!$OMP DO SCHEDULE(STATIC)
DO tracer=1,tr_lbc_vars
  DO k=tdims_l%k_start,tdims_l%k_end
    DO i=1,lenrim(fld_type_p,halo_type_extended)
      tracer_lbc(i,k,tracer)=tracer_lbc(i,k,tracer)+              &
        increment_factor*                                         &
        (tracer_lbc_tend(i,k,tracer)-tracer_lbc(i,k,tracer))
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO tracer=1,tr_lbc_ukca
  DO k=tdims_l%k_start,tdims_l%k_end
    DO i=1,lenrim(fld_type_p,halo_type_extended)
      tracer_ukca_lbc(i,k,tracer)=tracer_ukca_lbc(i,k,tracer)+    &
        increment_factor*                                         &
  (tracer_ukca_lbc_tend(i,k,tracer)-tracer_ukca_lbc(i,k,tracer))
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE increment_atmos_lbcs
