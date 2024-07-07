! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine BOUNDVAL
!
! Purpose : Checks whether a boundary incrementing step and increments
!           boundary values if required.
!
!
! ---------------------------------------------------------
!
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: LBC Input
SUBROUTINE boundval(                                              &
  lenrim                                                          &
, l_mcr_qcf2_lbc, l_mcr_qrain_lbc, l_mcr_qgraup_lbc, l_pc2_lbc    &
, l_murk_lbc, l_int_uvw_lbc                                       &
, l_dust_div1_lbc,l_dust_div2_lbc                                 &
, l_dust_div3_lbc,l_dust_div4_lbc                                 &
, l_dust_div5_lbc,l_dust_div6_lbc                                 &
, l_SO2_lbc,l_dms_lbc,l_SO4_aitken_lbc                            &
, l_SO4_accu_lbc,l_SO4_diss_lbc                                   &
, l_nh3_lbc,l_soot_new_lbc,l_soot_agd_lbc                         &
, l_soot_cld_lbc,l_bmass_new_lbc                                  &
, l_bmass_agd_lbc,l_bmass_cld_lbc                                 &
, l_ocff_new_lbc,l_ocff_agd_lbc,l_ocff_cld_lbc                    &
, l_nitr_acc_lbc, l_nitr_diss_lbc                                 &
, u_lbc, u_lbc_tend                                               &
, v_lbc, v_lbc_tend                                               &
, w_lbc, w_lbc_tend                                               &
, rho_lbc, rho_lbc_tend                                           &
, theta_lbc, theta_lbc_tend                                       &
, q_lbc, q_lbc_tend                                               &
, qcl_lbc, qcl_lbc_tend                                           &
, qcf_lbc, qcf_lbc_tend                                           &
, qcf2_lbc, qcf2_lbc_tend                                         &
, qrain_lbc, qrain_lbc_tend                                       &
, qgraup_lbc, qgraup_lbc_tend                                     &
, cf_bulk_lbc, cf_bulk_lbc_tend                                   &
, cf_liquid_lbc, cf_liquid_lbc_tend                               &
, cf_frozen_lbc, cf_frozen_lbc_tend                               &
, exner_lbc, exner_lbc_tend                                       &
, u_adv_lbc, u_adv_lbc_tend                                       &
, v_adv_lbc, v_adv_lbc_tend                                       &
, w_adv_lbc, w_adv_lbc_tend                                       &
, murk_lbc, murk_lbc_tend                                         &
, dust_div1_lbc, dust_div1_lbc_tend                               &
, dust_div2_lbc, dust_div2_lbc_tend                               &
, dust_div3_lbc, dust_div3_lbc_tend                               &
, dust_div4_lbc, dust_div4_lbc_tend                               &
, dust_div5_lbc, dust_div5_lbc_tend                               &
, dust_div6_lbc, dust_div6_lbc_tend                               &
, SO2_lbc, SO2_lbc_tend                                           &
, dms_lbc, dms_lbc_tend                                           &
, SO4_aitken_lbc, SO4_aitken_lbc_tend                             &
, SO4_accu_lbc, SO4_accu_lbc_tend                                 &
, SO4_diss_lbc, SO4_diss_lbc_tend                                 &
, NH3_lbc, NH3_lbc_tend                                           &
, soot_new_lbc, soot_new_lbc_tend                                 &
, soot_agd_lbc, soot_agd_lbc_tend                                 &
, soot_cld_lbc, soot_cld_lbc_tend                                 &
, bmass_new_lbc, bmass_new_lbc_tend                               &
, bmass_agd_lbc, bmass_agd_lbc_tend                               &
, bmass_cld_lbc, bmass_cld_lbc_tend                               &
, ocff_new_lbc, ocff_new_lbc_tend                                 &
, ocff_agd_lbc, ocff_agd_lbc_tend                                 &
, ocff_cld_lbc, ocff_cld_lbc_tend                                 &
, nitr_acc_lbc, nitr_acc_lbc_tend                                 &
, nitr_diss_lbc, nitr_diss_lbc_tend                               &
, tracer_lbc, tracer_lbc_tend                                     &
, tracer_ukca_lbc, tracer_ukca_lbc_tend                           &
, io1, io2                                                        &
, icode, cmessage)

USE submodel_mod, ONLY: atmos_im
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE lbc_mod
USE ereport_mod, ONLY: ereport
USE UM_ParParams
USE Control_Max_Sizes
USE atm_fields_bounds_mod, ONLY: udims_s, vdims_s, wdims_s,       &
                                 pdims_s, tdims_s, tdims_l
USE lbc_read_data_mod, ONLY: current_lbc_step
USE lookup_addresses

USE nlsizes_namelist_mod, ONLY:                                   &
    a_len_inthd, a_len_realhd, len1_lbc_comp_lookup, len1_lookup, &
    len_fixhd, model_levels, tr_lbc_ukca, tr_lbc_vars, tr_levels

USE model_time_mod, ONLY: &
    bndary_offsetim, stepim
USE errormessagelength_mod, ONLY: errormessagelength

USE atm_boundary_headers_mod, ONLY: rim_stepsa

IMPLICIT NONE

! Arguments:

INTEGER ::                                                        &
  lenrim(Nfld_max,NHalo_Max)                                      &
                            ! IN : Size of a level of LBC
, io1, io2                  ! Offsets to allow for old or new lbcs

LOGICAL, INTENT (IN) ::                                           &
  l_mcr_qcf2_lbc                                                  &
                   ! true if using second cloud ice lbcs
, l_mcr_qrain_lbc                                                 &
                   ! true if using rain lbcs
, l_mcr_qgraup_lbc                                                &
                   ! true if using graupel lbcs
, l_pc2_lbc                                                       &
                   ! true if cloud fractions in lbcs
, l_int_uvw_lbc                                                   &
                    ! true if using interpolated advecting winds
                    !  in lateral boundaries
, l_murk_lbc                                                      &
                    ! true if murk aerosol in lbcs
, l_dust_div1_lbc                                                 &
, l_dust_div2_lbc                                                 &
, l_dust_div3_lbc                                                 &
, l_dust_div4_lbc                                                 &
, l_dust_div5_lbc                                                 &
, l_dust_div6_lbc                                                 &
                   ! true if DUST in lbcs
, l_SO2_lbc                                                       &
                   ! true if SO2 in lbcs
, l_dms_lbc                                                       &
                   ! true if DMS in lbcs
, l_SO4_aitken_lbc                                                &
                   ! true if SO4_AITKEN in lbcs
, l_SO4_accu_lbc                                                  &
                   ! true if SO4_accU in lbcs
, l_SO4_diss_lbc                                                  &
                   ! true if SO4_DISS in lbcs
, l_nh3_lbc                                                       &
                   ! true if NH3 in lbcs
, l_soot_new_lbc                                                  &
, l_soot_agd_lbc                                                  &
, l_soot_cld_lbc                                                  &
                   ! true if soot in lbcs
, l_bmass_new_lbc                                                 &
, l_bmass_agd_lbc                                                 &
, l_bmass_cld_lbc                                                 &
                   ! true if biomass in lbcs
, l_ocff_new_lbc                                                  &
, l_ocff_agd_lbc                                                  &
, l_ocff_cld_lbc                                                  &
                   ! true if fossil fuel aerosol in lbcs
, l_nitr_acc_lbc                                                  &
, l_nitr_diss_lbc   ! true if fossil fuel aerosol in lbcs

! Note: LBCs are the current value of the LBC (and will be updated)
!       LBC_tends are the value of the LBC at the end of the LBC
!       period, towards which the LBC will tend.

REAL, INTENT (INOUT) ::                                           &
  u_lbc(lenrim(fld_type_u,halo_type_extended),                    &
        udims_s%k_start:udims_s%k_end)                            &
                              ! IN/OUT : U LBC
, u_lbc_tend(lenrim(fld_type_u,halo_type_extended),               &
                    udims_s%k_start:udims_s%k_end)                &
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

CHARACTER(LEN=errormessagelength) ::                              &
  cmessage                        ! Error message

! Local variables

INTEGER       :: steps_to_next_update

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BOUNDVAL'

! --------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

cmessage=' '

! Check that the lbcs haven't lost synchronisation and then call
! increment_lbcs to increment them

IF (current_lbc_step  /=  stepim(atmos_im) + io1) THEN

  IF ((current_lbc_STEP /= stepim(atmos_im) - io2) .AND.          &
         (MOD(bndary_offsetim(atmos_im)+stepim(atmos_im) + io1,   &
              rim_stepsa) /= 2) .AND.                             &
         (MOD(bndary_offsetim(atmos_im)+stepim(atmos_im) + io1,   &
              rim_stepsa) /= 0)) THEN

    icode=10

    CALL Ereport("BOUNDVA1 ", icode,                              &
                 "LBC values have lost synchronisation" )
  END IF

  steps_to_next_update = rim_stepsa -                             &
               MOD(bndary_offsetim(atmos_im)+stepim(atmos_im)-1,  &
                   rim_stepsa) + io2

  ! DEPENDS ON: increment_atmos_lbcs
  CALL increment_atmos_lbcs(                                      &
        steps_to_next_update,                                     &
        lenrim, tr_lbc_vars,tr_levels, tr_lbc_ukca,               &
        l_mcr_qcf2_lbc, l_mcr_qrain_lbc, l_mcr_qgraup_lbc,        &
        l_pc2_lbc, l_murk_lbc, l_int_uvw_lbc,                     &
        l_dust_div1_lbc,l_dust_div2_lbc,                          &
        l_dust_div3_lbc,l_dust_div4_lbc,                          &
        l_dust_div5_lbc,l_dust_div6_lbc,                          &
        l_SO2_lbc,l_dms_lbc,l_SO4_aitken_lbc,                     &
        l_SO4_accu_lbc,l_SO4_diss_lbc,                            &
        l_nh3_lbc,l_soot_new_lbc,l_soot_agd_lbc,                  &
        l_soot_cld_lbc,l_bmass_new_lbc,                           &
        l_bmass_agd_lbc,l_bmass_cld_lbc,                          &
        l_ocff_new_lbc,l_ocff_agd_lbc,l_ocff_cld_lbc,             &
        l_nitr_acc_lbc,l_nitr_diss_lbc,                           &
        u_lbc,u_lbc_tend,                                         &
        v_lbc,v_lbc_tend,                                         &
        w_lbc,w_lbc_tend,                                         &
        rho_lbc,rho_lbc_tend,                                     &
        theta_lbc,theta_lbc_tend,                                 &
        q_lbc,q_lbc_tend,                                         &
        qcl_lbc,qcl_lbc_tend,                                     &
        qcf_lbc,qcf_lbc_tend,                                     &
        qcf2_lbc,qcf2_lbc_tend,                                   &
        qrain_lbc,qrain_lbc_tend,                                 &
        qgraup_lbc,qgraup_lbc_tend,                               &
        cf_bulk_lbc,cf_bulk_lbc_tend,                             &
        cf_liquid_lbc,cf_liquid_lbc_tend,                         &
        cf_frozen_lbc,cf_frozen_lbc_tend,                         &
        exner_lbc,exner_lbc_tend,                                 &
        u_adv_lbc,u_adv_lbc_tend,                                 &
        v_adv_lbc,v_adv_lbc_tend,                                 &
        w_adv_lbc,w_adv_lbc_tend,                                 &
        murk_lbc,murk_lbc_tend,                                   &
        dust_div1_lbc,dust_div1_lbc_tend,                         &
        dust_div2_lbc,dust_div2_lbc_tend,                         &
        dust_div3_lbc,dust_div3_lbc_tend,                         &
        dust_div4_lbc,dust_div4_lbc_tend,                         &
        dust_div5_lbc,dust_div5_lbc_tend,                         &
        dust_div6_lbc,dust_div6_lbc_tend,                         &
        SO2_lbc,SO2_lbc_tend,                                     &
        dms_lbc,dms_lbc_tend,                                     &
        SO4_aitken_lbc,SO4_aitken_lbc_tend,                       &
        SO4_accu_lbc,SO4_accu_lbc_tend,                           &
        SO4_diss_lbc,SO4_diss_lbc_tend,                           &
        NH3_lbc,NH3_lbc_tend,                                     &
        soot_new_lbc,soot_new_lbc_tend,                           &
        soot_agd_lbc,soot_agd_lbc_tend,                           &
        soot_cld_lbc,soot_cld_lbc_tend,                           &
        bmass_new_lbc,bmass_new_lbc_tend,                         &
        bmass_agd_lbc,bmass_agd_lbc_tend,                         &
        bmass_cld_lbc,bmass_cld_lbc_tend,                         &
        ocff_new_lbc,ocff_new_lbc_tend,                           &
        ocff_agd_lbc,ocff_agd_lbc_tend,                           &
        ocff_cld_lbc,ocff_cld_lbc_tend,                           &
        nitr_acc_lbc, nitr_acc_lbc_tend,                          &
        nitr_diss_lbc, nitr_diss_lbc_tend,                        &
        tracer_lbc,tracer_lbc_tend,                               &
        tracer_ukca_lbc,tracer_ukca_lbc_tend,                     &
        rim_stepsa,icode)

  !       IF (ICODE  >   0) THEN
  !             WRITE(6,*) 'Failure in INCREMENT_ATMOS_lbcS while ',      &
  !    &                   'attempting to update LBCS for timestep ',     &
  !    &                   stepim(atmos_im)
  !         GOTO 9999
  !       ENDIF

  current_lbc_STEP = stepim(atmos_im) + io1

END IF ! IF (current_lbc_STEP  /=  stepim(atmos_im))

9999 CONTINUE
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE boundval
