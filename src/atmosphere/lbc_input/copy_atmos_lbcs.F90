! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine COPY_ATMOS_LBCS
!
! Purpose : Copies LBC_TEND array to LBC array, ready for new
!           tendencies to be read into LBC_TEND
!
! ---------------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: LBC Input

SUBROUTINE copy_atmos_lbcs(                                       &
  lenrim,                                                         &
  tr_lbc_vars,tr_levels,tr_lbc_ukca,                              &
  L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc, l_pc2_lbc,   &
  L_murk_lbc,                                                     &
  L_dust_div1_lbc,L_dust_div2_lbc,                                &
  L_dust_div3_lbc,L_dust_div4_lbc,                                &
  L_dust_div5_lbc,L_dust_div6_lbc,                                &
  L_so2_lbc,L_dms_lbc,L_so4_aitken_lbc,                           &
  L_so4_accu_lbc,L_so4_diss_lbc,                                  &
  L_nh3_lbc,L_soot_new_lbc,L_soot_agd_lbc,                        &
  L_soot_cld_lbc,L_bmass_new_lbc,                                 &
  L_bmass_agd_lbc,L_bmass_cld_lbc,                                &
  L_ocff_new_lbc,L_ocff_agd_lbc,L_ocff_cld_lbc,                   &
  L_nitr_acc_lbc, L_nitr_diss_lbc,                                &
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
  dust_div1_lbc, dust_div1_lbc_tend,                              &
  dust_div2_lbc, dust_div2_lbc_tend,                              &
  dust_div3_lbc, dust_div3_lbc_tend,                              &
  dust_div4_lbc, dust_div4_lbc_tend,                              &
  dust_div5_lbc, dust_div5_lbc_tend,                              &
  dust_div6_lbc, dust_div6_lbc_tend,                              &
  so2_lbc, so2_lbc_tend,                                          &
  dms_lbc, dms_lbc_tend,                                          &
  so4_aitken_lbc, so4_aitken_lbc_tend,                            &
  so4_accu_lbc, so4_accu_lbc_tend,                                &
  so4_diss_lbc, so4_diss_lbc_tend,                                &
  nh3_lbc, nh3_lbc_tend,                                          &
  soot_new_lbc, soot_new_lbc_tend,                                &
  soot_agd_lbc, soot_agd_lbc_tend,                                &
  soot_cld_lbc, soot_cld_lbc_tend,                                &
  bmass_new_lbc, bmass_new_lbc_tend,                              &
  bmass_agd_lbc, bmass_agd_lbc_tend,                              &
  bmass_cld_lbc, bmass_cld_lbc_tend,                              &
  ocff_new_lbc, ocff_new_lbc_tend,                                &
  ocff_agd_lbc, ocff_agd_lbc_tend,                                &
  ocff_cld_lbc, ocff_cld_lbc_tend,                                &
  nitr_acc_lbc, nitr_acc_lbc_tend,                                &
  nitr_diss_lbc, nitr_diss_lbc_tend,                              &
  tracer_lbc, tracer_lbc_tend,                                    &
  tracer_ukca_lbc, tracer_ukca_lbc_tend                           &
  )

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE atm_fields_bounds_mod, ONLY: udims_s, vdims_s, wdims_s,       &
                                 pdims_s, tdims_s, tdims_l
USE UM_ParParams
IMPLICIT NONE

! Parameters required for argument declarations

! Arguments:

INTEGER ::                                                        &
  lenrim(Nfld_max,NHalo_Max)                                      &
                              ! IN : Size of a level of LBC
, tr_lbc_vars                                                     &
                              ! IN : Number of tracer LBCs
, tr_lbc_ukca                                                     &
                              ! IN : Number of UKCA tracer LBCs
, tr_levels                   ! IN : Number of tracer levels

LOGICAL, INTENT (IN) ::                                           &
  L_mcr_qcf2_lbc                                                  &
                    ! true if second cloud ice lbcs active
, L_mcr_qrain_lbc                                                 &
                    ! true if rain lbcs active
, L_mcr_qgraup_lbc                                                &
                    ! true if graupel lbcs active
, L_pc2_lbc                                                       &
                    ! true if cloud fractions in lbcs
, L_murk_lbc                                                      &
                    ! true if murk aerosol in lbcs
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
, L_nitr_diss_lbc
                    ! true if nitrate aerosol in lbcs

! Note: LBCs are the current value of the LBC (and will be updated)
!       LBC_TENDs are the value of the LBC at the end of the LBC
!       period, towards which the LBC will tend.

REAL ::                                                           &
  u_lbc(lenrim(fld_type_u,halo_type_extended),                    &
        udims_s%k_start:udims_s%k_end)                            &
                              !  U LBC
, u_lbc_tend(lenrim(fld_type_u,halo_type_extended),               &
                    udims_s%k_start:udims_s%k_end)                &
                              !  U LBC tendency
, v_lbc(lenrim(fld_type_v,halo_type_extended),                    &
        vdims_s%k_start:vdims_s%k_end)                            &
                              !  V LBC
, v_lbc_tend(lenrim(fld_type_v,halo_type_extended),               &
             vdims_s%k_start:vdims_s%k_end)                       &
                              !  V LBC tendency
, w_lbc(lenrim(fld_type_p,halo_type_extended),                    &
        wdims_s%k_start:wdims_s%k_end)                            &
                              !  V LBC
, w_lbc_tend(lenrim(fld_type_p,halo_type_extended),               &
             wdims_s%k_start:wdims_s%k_end)                       &
                              !  V LBC tendency
, rho_lbc(lenrim(fld_type_p,halo_type_extended),                  &
          pdims_s%k_start:pdims_s%k_end)                          &
                              !  rho LBC
, rho_lbc_tend(lenrim(fld_type_p,halo_type_extended),             &
               pdims_s%k_start:pdims_s%k_end)                     &
                              !  rho LBC tendency
, theta_lbc(lenrim(fld_type_p,halo_type_extended),                &
            tdims_s%k_start:tdims_s%k_end)                        &
                              !  theta LBC
, theta_lbc_tend(lenrim(fld_type_p,halo_type_extended),           &
                 tdims_s%k_start:tdims_s%k_end)                   &
                              !  theta LBC tendency
, q_lbc(lenrim(fld_type_p,halo_type_extended),                    &
        tdims_l%k_start:tdims_l%k_end)                            &
                              !  Q LBC
, q_lbc_tend(lenrim(fld_type_p,halo_type_extended),               &
             tdims_l%k_start:tdims_l%k_end)                       &
                              !  Q LBC tendency
, qcl_lbc(lenrim(fld_type_p,halo_type_extended),                  &
          tdims_l%k_start:tdims_l%k_end)                          &
                              !  QCL LBC
, qcl_lbc_tend(lenrim(fld_type_p,halo_type_extended),             &
               tdims_l%k_start:tdims_l%k_end)                     &
                              !  QCL LBC tendency
, qcf_lbc(lenrim(fld_type_p,halo_type_extended),                  &
          tdims_l%k_start:tdims_l%k_end)                          &
                              !  QCL LBC
, qcf_lbc_tend(lenrim(fld_type_p,halo_type_extended),             &
               tdims_l%k_start:tdims_l%k_end)                     &
                              !  QCL LBC tendency
, qcf2_lbc(lenrim(fld_type_p,halo_type_extended),                 &
           tdims_l%k_start:tdims_l%k_end)                         &
                              !  QCF2 LBC
, qcf2_lbc_tend(lenrim(fld_type_p,halo_type_extended),            &
                tdims_l%k_start:tdims_l%k_end)                    &
                              !  QCF2 LBC tendency
, qrain_lbc(lenrim(fld_type_p,halo_type_extended),                &
            tdims_l%k_start:tdims_l%k_end)                        &
                              !  QRAIN LBC
, qrain_lbc_tend(lenrim(fld_type_p,halo_type_extended),           &
                 tdims_l%k_start:tdims_l%k_end)                   &
                              !  QRAIN LBC tendency
, qgraup_lbc(lenrim(fld_type_p,halo_type_extended),               &
             tdims_l%k_start:tdims_l%k_end)                       &
                              !  QGRAUP LBC
, qgraup_lbc_tend(lenrim(fld_type_p,halo_type_extended),          &
                  tdims_l%k_start:tdims_l%k_end)                  &
                              !  QGRAUP LBC tendency
, cf_bulk_lbc(lenrim(fld_type_p,halo_type_extended),              &
              tdims_l%k_start:tdims_l%k_end)                      &
                              !  CF_BULK LBC
, cf_bulk_lbc_tend(lenrim(fld_type_p,halo_type_extended),         &
                   tdims_l%k_start:tdims_l%k_end)                 &
                              !  CF_BULK LBC tendency
, cf_liquid_lbc(lenrim(fld_type_p,halo_type_extended),            &
                tdims_l%k_start:tdims_l%k_end)                    &
                              !  CF_LIQUID LBC
, cf_liquid_lbc_tend(lenrim(fld_type_p,halo_type_extended),       &
                            tdims_l%k_start:tdims_l%k_end)        &
                              !  CF_LIQUID LBC tendency
, cf_frozen_lbc(lenrim(fld_type_p,halo_type_extended),            &
                tdims_l%k_start:tdims_l%k_end)                    &
                              !  CF_FROZEN LBC
, cf_frozen_lbc_tend(lenrim(fld_type_p,halo_type_extended),       &
                     tdims_l%k_start:tdims_l%k_end)               &
                              !  CF_FROZEN LBC tendency
, exner_lbc(lenrim(fld_type_p,halo_type_extended),                &
            pdims_s%k_start:pdims_s%k_end+1)                      &
                              !  Exner LBC
, exner_lbc_tend(lenrim(fld_type_p,halo_type_extended),           &
                 pdims_s%k_start:pdims_s%k_end+1)                 &
                                   !  Exner LBC tendency
, u_adv_lbc(lenrim(fld_type_u,halo_type_extended),                &
            udims_s%k_start:udims_s%k_end)                        &
                              !  u_adv LBC
, u_adv_lbc_tend(lenrim(fld_type_u,halo_type_extended),           &
                 udims_s%k_start:udims_s%k_end)                   &
                              !  u_adv LBC tendency
, v_adv_lbc(lenrim(fld_type_v,halo_type_extended),                &
            vdims_s%k_start:vdims_s%k_end)                        &
                              !  v_adv LBC
, v_adv_lbc_tend(lenrim(fld_type_v,halo_type_extended),           &
                 vdims_s%k_start:vdims_s%k_end)                   &
                              !  v_adv LBC tendency
, w_adv_lbc(lenrim(fld_type_p,halo_type_extended),                &
                   wdims_s%k_start:wdims_s%k_end)                 &
                              !  W LBC
, w_adv_lbc_tend(lenrim(fld_type_p,halo_type_extended),           &
                 wdims_s%k_start:wdims_s%k_end)                   &
                                 !  W LBC tendency
, murk_lbc(lenrim(fld_type_p,halo_type_single),                   &
           tdims_s%k_start:tdims_s%k_end)                         &
                              !  MURK LBC
, murk_lbc_tend(lenrim(fld_type_p,halo_type_single),              &
                tdims_s%k_start:tdims_s%k_end)                    &
                              !  MURK LBC tendency
, dust_div1_lbc(lenrim(fld_type_p,halo_type_single),              &
                tdims_s%k_start:tdims_s%k_end)                    &
                              !  dust_div1 LBC
, dust_div1_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                     tdims_s%k_start:tdims_s%k_end)               &
                              !  dust_div1 LBC tendency
, dust_div2_lbc(lenrim(fld_type_p,halo_type_single),              &
                tdims_s%k_start:tdims_s%k_end)                    &
                              !  dust_div2 LBC
, dust_div2_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                     tdims_s%k_start:tdims_s%k_end)               &
                              !  dust_div2 LBC tendency
, dust_div3_lbc(lenrim(fld_type_p,halo_type_single),              &
                tdims_s%k_start:tdims_s%k_end)                    &
                              !  dust_div3 LBC
, dust_div3_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                     tdims_s%k_start:tdims_s%k_end)               &
                              !  dust_div3 LBC tendency
, dust_div4_lbc(lenrim(fld_type_p,halo_type_single),              &
                tdims_s%k_start:tdims_s%k_end)                    &
                              !  dust_div4 LBC
, dust_div4_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                     tdims_s%k_start:tdims_s%k_end)               &
                              !  dust_div4 LBC tendency
, dust_div5_lbc(lenrim(fld_type_p,halo_type_single),              &
                tdims_s%k_start:tdims_s%k_end)                    &
                              !  dust_div5 LBC
, dust_div5_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                     tdims_s%k_start:tdims_s%k_end)               &
                              !  dust_div5 LBC tendency
, dust_div6_lbc(lenrim(fld_type_p,halo_type_single),              &
                tdims_s%k_start:tdims_s%k_end)                    &
                              !  dust_div6 LBC
, dust_div6_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                     tdims_s%k_start:tdims_s%k_end)               &
                              !  dust_div6 LBC tendency
, SO2_lbc(lenrim(fld_type_p,halo_type_single),                    &
          tdims_s%k_start:tdims_s%k_end)                          &
                              !  SO2 LBC
, SO2_lbc_tend(lenrim(fld_type_p,halo_type_single),               &
               tdims_s%k_start:tdims_s%k_end)                     &
                              !  SO2 LBC tendency
, dms_lbc(lenrim(fld_type_p,halo_type_single),                    &
          tdims_s%k_start:tdims_s%k_end)                          &
                              !  DMS LBC
, dms_lbc_tend(lenrim(fld_type_p,halo_type_single),               &
               tdims_s%k_start:tdims_s%k_end)                     &
                              !  DMS LBC tendency
, SO4_aitken_lbc(lenrim(fld_type_p,halo_type_single),             &
                 tdims_s%k_start:tdims_s%k_end)                   &
                              !  SO4_AITKEN LBC
, SO4_aitken_lbc_tend(lenrim(fld_type_p,halo_type_single),        &
                      tdims_s%k_start:tdims_s%k_end)              &
                              !  SO4_AITKEN LBC tendency
, SO4_accu_lbc(lenrim(fld_type_p,halo_type_single),               &
               tdims_s%k_start:tdims_s%k_end)                     &
                              !  SO4_accU LBC
, SO4_accu_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                    tdims_s%k_start:tdims_s%k_end)                &
                              !  SO4_accU LBC tendency
, SO4_diss_lbc(lenrim(fld_type_p,halo_type_single),               &
               tdims_s%k_start:tdims_s%k_end)                     &
                              !  SO4_DISS LBC
, SO4_diss_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                    tdims_s%k_start:tdims_s%k_end)                &
                              !  SO4_DISS LBC tendency
, NH3_lbc(lenrim(fld_type_p,halo_type_single),                    &
          tdims_s%k_start:tdims_s%k_end)                          &
                              !  NH3 LBC
, NH3_lbc_tend(lenrim(fld_type_p,halo_type_single),               &
               tdims_s%k_start:tdims_s%k_end)                     &
                              !  NH3 LBC tendency
, soot_new_lbc(lenrim(fld_type_p,halo_type_single),               &
               tdims_s%k_start:tdims_s%k_end)                     &
                              !  soot_NEW LBC
, soot_new_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                    tdims_s%k_start:tdims_s%k_end)                &
                              !  soot_NEW LBC tendency
, soot_agd_lbc(lenrim(fld_type_p,halo_type_single),               &
               tdims_s%k_start:tdims_s%k_end)                     &
                              !  soot_agd LBC
, soot_agd_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                    tdims_s%k_start:tdims_s%k_end)                &
                              !  soot_agd LBC tendency
, soot_cld_lbc(lenrim(fld_type_p,halo_type_single),               &
               tdims_s%k_start:tdims_s%k_end)                     &
                              !  soot_CLD LBC
, soot_cld_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                    tdims_s%k_start:tdims_s%k_end)                &
                              !  soot_CLD LBC tendency
, bmass_new_lbc(lenrim(fld_type_p,halo_type_single),              &
                tdims_s%k_start:tdims_s%k_end)                    &
                              !  bmass_NEW LBC
, bmass_new_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                     tdims_s%k_start:tdims_s%k_end)               &
                              !  bmass_NEW LBC tendency
, bmass_agd_lbc(lenrim(fld_type_p,halo_type_single),              &
                tdims_s%k_start:tdims_s%k_end)                    &
                              !  bmass_agd LBC
, bmass_agd_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                     tdims_s%k_start:tdims_s%k_end)               &
                              !  bmass_agd LBC tendency
, bmass_cld_lbc(lenrim(fld_type_p,halo_type_single),              &
                tdims_s%k_start:tdims_s%k_end)                    &
                              !  bmass_CLD LBC
, bmass_cld_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                     tdims_s%k_start:tdims_s%k_end)               &
                              !  bmass_CLD LBC tendency
, ocff_new_lbc(lenrim(fld_type_p,halo_type_single),               &
               tdims_s%k_start:tdims_s%k_end)                     &
                              !  ocff_NEW LBC
, ocff_new_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                    tdims_s%k_start:tdims_s%k_end)                &
                              !  ocff_NEW LBC tendency
, ocff_agd_lbc(lenrim(fld_type_p,halo_type_single),               &
               tdims_s%k_start:tdims_s%k_end)                     &
                              !  ocff_agd LBC
, ocff_agd_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                    tdims_s%k_start:tdims_s%k_end)                &
                              !  ocff_agd LBC tendency
, ocff_cld_lbc(lenrim(fld_type_p,halo_type_single),               &
               tdims_s%k_start:tdims_s%k_end)                     &
                              !  ocff_CLD LBC
, ocff_cld_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                    tdims_s%k_start:tdims_s%k_end)                &
                              !  ocff_CLD LBC tendency
, nitr_acc_lbc(lenrim(fld_type_p,halo_type_single),               &
               tdims_s%k_start:tdims_s%k_end)                     &
                              !  nitr_acc LBC
, nitr_acc_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                    tdims_s%k_start:tdims_s%k_end)                &
                              !  nitr_acc LBC tendency
, nitr_diss_lbc(lenrim(fld_type_p,halo_type_single),              &
                tdims_s%k_start:tdims_s%k_end)                    &
                              !  nitr_DISS LBC
, nitr_diss_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                     tdims_s%k_start:tdims_s%k_end)               &
                              !  nitr_DISS LBC tendency
, tracer_lbc(lenrim(fld_type_p,halo_type_extended),               &
                    tdims_l%k_start:tdims_l%k_end,tr_lbc_vars)    &
                                     !  Tracer LBCs
, tracer_lbc_tend(lenrim(fld_type_p,halo_type_extended),          &
                  tdims_l%k_start:tdims_l%k_end,tr_lbc_vars)      &
                                     !  Tracer LBCs tendency
, tracer_ukca_lbc(lenrim(fld_type_p,halo_type_extended),          &
                  tdims_l%k_start:tdims_l%k_end,tr_lbc_ukca)      &
                                     !  UKCA Tracer LBCs
, tracer_ukca_lbc_tend(lenrim(fld_type_p,halo_type_extended),     &
                       tdims_l%k_start:tdims_l%k_end,tr_lbc_ukca)
                                     !  UKCA Tracer LBCs tend

! Local variables

INTEGER ::                                                        &
  tracer                                                          &
            ! loop counter over tracer variables
, k                                                               &
            ! loop counter over levels
, i         ! loop counter over horizontal dimension

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_ATMOS_LBCS'

! ---------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

DO k=udims_s%k_start,udims_s%k_end
  DO i=1,lenrim(fld_type_u,halo_type_extended)
    u_lbc(i,k)=u_lbc_tend(i,k)
    u_adv_lbc(i,k)=u_adv_lbc_tend(i,k)
  END DO ! i
END DO !k=udims_s%k_start,udims_s%k_end

DO k=vdims_s%k_start,vdims_s%k_end
  DO i=1,lenrim(fld_type_v,halo_type_extended)
    v_lbc(i,k)=v_lbc_tend(i,k)
    v_adv_lbc(i,k)=v_adv_lbc_tend(i,k)
  END DO ! i
END DO ! k

DO k=pdims_s%k_start,pdims_s%k_end
  DO i=1,lenrim(fld_type_p,halo_type_extended)
    rho_lbc(i,k)=rho_lbc_tend(i,k)
  END DO ! i
END DO ! k

DO k=tdims_s%k_start,tdims_s%k_end
  DO i=1,lenrim(fld_type_p,halo_type_extended)
    theta_lbc(i,k)=theta_lbc_tend(i,k)
  END DO ! i
END DO ! k

DO k=tdims_l%k_start,tdims_l%k_end
  DO i=1,lenrim(fld_type_p,halo_type_extended)
    q_lbc(i,k)=q_lbc_tend(i,k)
    qcl_lbc(i,k)=qcl_lbc_tend(i,k)
    qcf_lbc(i,k)=qcf_lbc_tend(i,k)
  END DO ! i
END DO ! k

IF (L_mcr_qcf2_lbc) THEN  ! qcf2 lbcs active
  DO k=tdims_l%k_start,tdims_l%k_end
    DO i=1,lenrim(fld_type_p,halo_type_extended)
      qcf2_lbc(i,k) = qcf2_lbc_tend(i,k)
    END DO ! i
  END DO ! k
END IF

IF (L_mcr_qrain_lbc) THEN  ! qrain lbcs active
  DO k=tdims_l%k_start,tdims_l%k_end
    DO i=1,lenrim(fld_type_p,halo_type_extended)
      qrain_lbc(i,k) = qrain_lbc_tend(i,k)
    END DO ! i
  END DO ! k
END IF

IF (L_mcr_qgraup_lbc) THEN  ! qgraup lbcs active
  DO k=tdims_l%k_start,tdims_l%k_end
    DO i=1,lenrim(fld_type_p,halo_type_extended)
      qgraup_lbc(i,k) = qgraup_lbc_tend(i,k)
    END DO ! i
  END DO ! k
END IF

IF (L_pc2_lbc) THEN  ! cloud fraction lbc's active
  DO k=tdims_l%k_start,tdims_l%k_end
    DO i=1,lenrim(fld_type_p,halo_type_extended)
      cf_bulk_lbc(i,k)   = cf_bulk_lbc_tend(i,k)
      cf_liquid_lbc(i,k) = cf_liquid_lbc_tend(i,k)
      cf_frozen_lbc(i,k) = cf_frozen_lbc_tend(i,k)
    END DO ! i
  END DO ! k
END IF

DO k=wdims_s%k_start,wdims_s%k_end
  DO i=1,lenrim(fld_type_p,halo_type_extended)
    w_lbc(i,k)=w_lbc_tend(i,k)
    w_adv_lbc(i,k)=w_adv_lbc_tend(i,k)
  END DO ! i
END DO ! k

DO k=pdims_s%k_start,pdims_s%k_end+1
  DO i=1,lenrim(fld_type_p,halo_type_extended)
    exner_lbc(i,k)=exner_lbc_tend(i,k)
  END DO ! i
END DO ! k

IF (L_murk_lbc) THEN
  ! murk aerosol turned on and murk lbcs active
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      murk_lbc(i,k) = murk_lbc_tend(i,k)
    END DO ! i
  END DO ! k
END IF

IF (L_dust_div1_lbc) THEN
  ! dust turned on and dust lbcs active
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      dust_div1_lbc(i,k) = dust_div1_lbc_tend(i,k)
    END DO ! i
  END DO ! k
END IF

IF (L_dust_div2_lbc) THEN
  ! dust turned on and dust lbcs active
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      dust_div2_lbc(i,k) = dust_div2_lbc_tend(i,k)
    END DO ! i
  END DO ! k
END IF

IF (L_dust_div3_lbc) THEN
  ! dust turned on and dust lbcs active
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      dust_div3_lbc(i,k) = dust_div3_lbc_tend(i,k)
    END DO ! i
  END DO ! k
END IF

IF (L_dust_div4_lbc) THEN
  ! dust turned on and dust lbcs active
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      dust_div4_lbc(i,k) = dust_div4_lbc_tend(i,k)
    END DO ! i
  END DO ! k
END IF

IF (L_dust_div5_lbc) THEN
  ! dust turned on and dust lbcs active
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      dust_div5_lbc(i,k) = dust_div5_lbc_tend(i,k)
    END DO ! i
  END DO ! k
END IF

IF (L_dust_div6_lbc) THEN
  ! dust turned on and dust lbcs active
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      dust_div6_lbc(i,k) = dust_div6_lbc_tend(i,k)
    END DO ! i
  END DO ! k
END IF

IF (L_so2_lbc) THEN
  ! so2 turned on and so2 lbcs active
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      so2_lbc(i,k) = so2_lbc_tend(i,k)
    END DO ! i
  END DO ! k
END IF

IF (L_dms_lbc) THEN
  ! dms turned on and dms lbcs active
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      dms_lbc(i,k) = dms_lbc_tend(i,k)
    END DO ! i
  END DO ! k
END IF

IF (L_so4_aitken_lbc) THEN
  ! so4_aitken turned on and so4_aitken lbcs active
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      so4_aitken_lbc(i,k) = so4_aitken_lbc_tend(i,k)
    END DO ! i
  END DO ! k
END IF

IF (L_so4_accu_lbc) THEN
  ! so4_accu turned on and so4_accu lbcs active
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      so4_accu_lbc(i,k) = so4_accu_lbc_tend(i,k)
    END DO ! i
  END DO ! k
END IF

IF (L_so4_diss_lbc) THEN
  ! so4_diss turned on and so4_diss lbcs active
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      so4_diss_lbc(i,k) = so4_diss_lbc_tend(i,k)
    END DO ! i
  END DO ! k
END IF

IF (L_nh3_lbc) THEN
  ! nh3 turned on and nh3 lbcs active
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      nh3_lbc(i,k) = nh3_lbc_tend(i,k)
    END DO ! i
  END DO ! k
END IF

IF (L_soot_new_lbc) THEN
  ! soot_new turned on and soot_new lbcs active
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      soot_new_lbc(i,k) = soot_new_lbc_tend(i,k)
    END DO ! i
  END DO ! k
END IF

IF (L_soot_agd_lbc) THEN
  ! soot_agd turned on and soot_agd lbcs active
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      soot_agd_lbc(i,k) = soot_agd_lbc_tend(i,k)
    END DO ! i
  END DO ! k
END IF

IF (L_soot_cld_lbc) THEN
  ! soot_cld turned on and soot_cld lbcs active
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      soot_cld_lbc(i,k) = soot_cld_lbc_tend(i,k)
    END DO ! i
  END DO ! k
END IF

IF (L_bmass_new_lbc) THEN
  ! bmass_new turned on and bmass_new lbcs active
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      bmass_new_lbc(i,k) = bmass_new_lbc_tend(i,k)
    END DO ! i
  END DO ! k
END IF

IF (L_bmass_agd_lbc) THEN
  ! bmass_agd turned on and bmass_agd lbcs active
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      bmass_agd_lbc(i,k) = bmass_agd_lbc_tend(i,k)
    END DO ! i
  END DO ! k
END IF

IF (L_bmass_cld_lbc) THEN
  ! bmass_cld turned on and bmass_cld lbcs active
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      bmass_cld_lbc(i,k) = bmass_cld_lbc_tend(i,k)
    END DO ! i
  END DO ! k
END IF

IF (L_ocff_new_lbc) THEN
  ! ocff_new turned on and ocff_new lbcs active
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      ocff_new_lbc(i,k) = ocff_new_lbc_tend(i,k)
    END DO ! i
  END DO ! k
END IF

IF (L_ocff_agd_lbc) THEN
  ! ocff_agd turned on and ocff_agd lbcs active
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      ocff_agd_lbc(i,k) = ocff_agd_lbc_tend(i,k)
    END DO ! i
  END DO ! k
END IF

IF (L_ocff_cld_lbc) THEN
  ! ocff_cld turned on and ocff_cld lbcs active
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      ocff_cld_lbc(i,k) = ocff_cld_lbc_tend(i,k)
    END DO ! i
  END DO ! k
END IF

IF (L_nitr_acc_lbc) THEN
  ! nitr_acc turned on and nitr_acc lbcs active
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      nitr_acc_lbc(i,k) = nitr_acc_lbc_tend(i,k)
    END DO ! i
  END DO ! k
END IF

IF (L_nitr_diss_lbc) THEN
  ! nitr_diss turned on and nitr_diss lbcs active
  DO k=tdims_s%k_start,tdims_s%k_end
    DO i=1,lenrim(fld_type_p,halo_type_single)
      nitr_diss_lbc(i,k) = nitr_diss_lbc_tend(i,k)
    END DO ! i
  END DO ! k
END IF

DO tracer=1,tr_lbc_vars
  DO k=tdims_l%k_start,tdims_l%k_end
    DO i=1,lenrim(fld_type_p,halo_type_extended)
      tracer_lbc(i,k,tracer)=tracer_lbc_tend(i,k,tracer)
    END DO
  END DO
END DO

DO tracer=1,tr_lbc_ukca
  DO k=tdims_l%k_start,tdims_l%k_end
    DO i=1,lenrim(fld_type_p,halo_type_extended)
      tracer_ukca_lbc(i,k,tracer)=tracer_ukca_lbc_tend(i,k,tracer)
    END DO
  END DO
END DO

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE copy_atmos_lbcs
