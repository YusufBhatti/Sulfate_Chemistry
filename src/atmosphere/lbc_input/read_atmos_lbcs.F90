! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine READ_ATMOS_LBCS
!
! Purpose : Reads in all the Atmos LBCs
!
! ---------------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: LBC Input

SUBROUTINE read_atmos_lbcs(                                       &
  lenrim, global_LENRIM,                                          &
  tr_lbc_vars,tr_levels,tr_lbc_ukca,                              &
  L_int_uvw_lbc, L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,           &
  L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc, l_pc2_lbc,   &
  L_murk, L_murk_lbc,                                             &
  L_dust_div1, L_dust_div1_lbc, L_dust_div2, L_dust_div2_lbc,     &
  L_dust_div3, L_dust_div3_lbc, L_dust_div4, L_dust_div4_lbc,     &
  L_dust_div5, L_dust_div5_lbc, L_dust_div6, L_dust_div6_lbc,     &
  L_so2,L_so2_lbc,L_dms,L_dms_lbc,L_so4_aitken,L_so4_aitken_lbc,  &
  L_so4_accu, L_so4_accu_lbc, L_so4_diss, L_so4_diss_lbc,         &
  L_nh3, L_nh3_lbc, L_soot_new, L_soot_new_lbc, L_soot_agd,       &
  L_soot_agd_lbc, L_soot_cld, L_soot_cld_lbc, L_bmass_new,        &
  L_bmass_new_lbc, L_bmass_agd, L_bmass_agd_lbc, L_bmass_cld,     &
  L_bmass_cld_lbc, L_ocff_new, L_ocff_new_lbc, L_ocff_agd,        &
  L_ocff_agd_lbc, L_ocff_cld, L_ocff_cld_lbc,                     &
  L_nitr_acc, L_nitr_acc_lbc, L_nitr_diss, L_nitr_diss_lbc,       &
  len1_lookup,len_fixhd,nftin, n_lookups,lookup, fixhd,           &
  u_lbc,v_lbc,w_lbc,rho_lbc,theta_lbc,q_lbc,qcl_lbc,qcf_lbc,      &
  qcf2_lbc, qrain_lbc, qgraup_lbc,                                &
  cf_bulk_lbc, cf_liquid_lbc, cf_frozen_lbc,                      &
  exner_lbc, u_adv_lbc, v_adv_lbc, w_adv_lbc,                     &
  murk_lbc,                                                       &
  dust_div1_lbc, dust_div2_lbc, dust_div3_lbc,                    &
  dust_div4_lbc, dust_div5_lbc, dust_div6_lbc,                    &
  so2_lbc, dms_lbc, so4_aitken_lbc, so4_accu_lbc, so4_diss_lbc,   &
  nh3_lbc, soot_new_lbc, soot_agd_lbc, soot_cld_lbc,              &
  bmass_new_lbc, bmass_agd_lbc, bmass_cld_lbc,                    &
  ocff_new_lbc, ocff_agd_lbc, ocff_cld_lbc,                       &
  nitr_acc_lbc, nitr_diss_lbc,                                    &
  tracer_lbc, tracer_ukca_lbc,                                    &
  icode,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE atm_fields_bounds_mod, ONLY: udims_s, vdims_s, wdims_s,       &
                                 pdims_s, tdims_s, tdims_l
USE cloud_inputs_mod, ONLY: i_cld_vn
USE pc2_constants_mod, ONLY: i_cld_pc2
USE UM_ParParams
USE lookup_addresses

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

USE errormessagelength_mod, ONLY: errormessagelength

USE item_bounda_mod

USE readflds_mod, ONLY: readflds

USE d1_array_mod, ONLY: other

IMPLICIT NONE


! Parameters required for argument declarations

! Arguments:

INTEGER ::                                                        &
  lenrim(Nfld_max,NHalo_Max)                                      &
                              ! IN : Size of a level of LBC
, global_LENRIM(Nfld_max,NHalo_Max)                               &
                              ! IN : Full (disk) size of a
                              !      level of LBC
, tr_lbc_vars                                                     &
                              ! IN : Number of tracer LBCs
, tr_lbc_ukca                                                     &
                              ! IN : Number of UKCA tracer LBCs
, tr_levels                                                       &
                              ! IN : Number of tracer levels
, len1_lookup                                                     &
                              ! IN : Size of LOOKUP header
, len_fixhd                                                       &
                              ! IN : Size of fixed length header
, nftin                                                           &
                              ! IN : Unit for LBC file
, n_lookups                                                       &
                              ! IN  : Number of lookup headers
, lookup(len1_lookup,n_lookups)                                   &
                              ! IN : LOOKUP headers
, fixhd(len_fixhd)            ! IN : Fixed header

LOGICAL, INTENT (IN) ::                                           &
  L_int_uvw_lbc                                                   &
                    ! if TRUE, do not need advected winds
, L_mcr_qcf2                                                      &
                    ! true if second cloud ice active
, L_mcr_qrain                                                     &
                    ! true if rain active
, L_mcr_qgraup                                                    &
                    ! true if graupel active
, L_mcr_qcf2_lbc                                                  &
                    ! true if second cloud ice lbcs active
, L_mcr_qrain_lbc                                                 &
                    ! true if rain lbcs active
, L_mcr_qgraup_lbc                                                &
                    ! true if graupel lbcs active
, L_pc2_lbc                                                       &
                    ! true if cloud fractions in lbcs
, L_murk                                                          &
                    ! true if prognostic murk aerosol active
, L_murk_lbc                                                      &
                    ! true if murk aerosol in lbcs
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
, L_soot_new                                                      &
, L_soot_agd                                                      &
, L_soot_cld                                                      &
                    ! true if soot active
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
, L_ocff_new                                                      &
, L_ocff_agd                                                      &
, L_ocff_cld                                                      &
                    ! true if fossil fuel aerosol active
, L_ocff_new_lbc                                                  &
, L_ocff_agd_lbc                                                  &
, L_ocff_cld_lbc                                                  &
                    ! true if fossil fuel aerosol lbcs active
, L_nitr_acc                                                      &
, L_nitr_diss                                                     &
                    ! true if nitrate aerosol active
, L_nitr_acc_lbc                                                  &
, L_nitr_diss_lbc   ! true if nitrate aerosol lbcs active

REAL ::                                                           &
  u_lbc(lenrim(fld_type_u,halo_type_extended),                    &
        udims_s%k_start:udims_s%k_end)                            &
                              ! OUT : U LBC
, v_lbc(lenrim(fld_type_v,halo_type_extended),                    &
        vdims_s%k_start:vdims_s%k_end)                            &
                              ! OUT : V LBC
, w_lbc(lenrim(fld_type_p,halo_type_extended),                    &
        wdims_s%k_start:wdims_s%k_end)                            &
                              ! OUT : V LBC
, rho_lbc(lenrim(fld_type_p,halo_type_extended),                  &
          pdims_s%k_start:pdims_s%k_end)                          &
                              ! OUT : Rho LBC
, theta_lbc(lenrim(fld_type_p,halo_type_extended),                &
            tdims_s%k_start:tdims_s%k_end)                        &
                              ! OUT : Theta LBC
, q_lbc(lenrim(fld_type_p,halo_type_extended),                    &
        tdims_l%k_start:tdims_l%k_end)                            &
                              ! OUT : Q LBC
, qcl_lbc(lenrim(fld_type_p,halo_type_extended),                  &
          tdims_l%k_start:tdims_l%k_end)                          &
                              ! OUT : QCL LBC
, qcf_lbc(lenrim(fld_type_p,halo_type_extended),                  &
          tdims_l%k_start:tdims_l%k_end)                          &
                              ! OUT : QCL LBC
, qcf2_lbc(lenrim(fld_type_p,halo_type_extended),                 &
           tdims_l%k_start:tdims_l%k_end)                         &
                              ! OUT : QCF2 LBC
, qrain_lbc(lenrim(fld_type_p,halo_type_extended),                &
            tdims_l%k_start:tdims_l%k_end)                        &
                              ! OUT : QRAIN LBC
, qgraup_lbc(lenrim(fld_type_p,halo_type_extended),               &
             tdims_l%k_start:tdims_l%k_end)                       &
                              ! OUT : QGRAUP LBC
, cf_bulk_lbc(lenrim(fld_type_p,halo_type_extended),              &
              tdims_l%k_start:tdims_l%k_end)                      &
                              ! OUT : CF_BULK LBC
, cf_liquid_lbc(lenrim(fld_type_p,halo_type_extended),            &
                tdims_l%k_start:tdims_l%k_end)                    &
                              ! OUT : CF_LIQUID LBC
, cf_frozen_lbc(lenrim(fld_type_p,halo_type_extended),            &
                tdims_l%k_start:tdims_l%k_end)                    &
                              ! OUT : CF_FROZEN LBC
, exner_lbc(lenrim(fld_type_p,halo_type_extended),                &
            pdims_s%k_start:pdims_s%k_end+1)                      &
                              ! OUT : Exner LBC
, u_adv_lbc(lenrim(fld_type_u,halo_type_extended),                &
            udims_s%k_start:udims_s%k_end)                        &
                              ! OUT : U_ADV LBC
, v_adv_lbc(lenrim(fld_type_v,halo_type_extended),                &
            vdims_s%k_start:vdims_s%k_end)                        &
                              ! OUT : V_ADV LBC
, w_adv_lbc(lenrim(fld_type_p,halo_type_extended),                &
                   wdims_s%k_start:wdims_s%k_end)                 &
                              ! OUT : W LBC
, murk_lbc(lenrim(fld_type_p,halo_type_single),                   &
           tdims_s%k_start:tdims_s%k_end)                         &
                                 ! OUT : MURK LBC
, dust_div1_lbc(lenrim(fld_type_p,halo_type_single),              &
                tdims_s%k_start:tdims_s%k_end)                    &
                              ! OUT : DUST_DIV1 LBC
, dust_div2_lbc(lenrim(fld_type_p,halo_type_single),              &
                tdims_s%k_start:tdims_s%k_end)                    &
                              ! OUT : DUST_DIV2 LBC
, dust_div3_lbc(lenrim(fld_type_p,halo_type_single),              &
                tdims_s%k_start:tdims_s%k_end)                    &
                              ! OUT : DUST_DIV3 LBC
, dust_div4_lbc(lenrim(fld_type_p,halo_type_single),              &
                tdims_s%k_start:tdims_s%k_end)                    &
                              ! OUT : DUST_DIV4 LBC
, dust_div5_lbc(lenrim(fld_type_p,halo_type_single),              &
                tdims_s%k_start:tdims_s%k_end)                    &
                              ! OUT : DUST_DIV5 LBC
, dust_div6_lbc(lenrim(fld_type_p,halo_type_single),              &
                tdims_s%k_start:tdims_s%k_end)                    &
                              ! OUT : DUST_DIV6 LBC
, so2_lbc(lenrim(fld_type_p,halo_type_single),                    &
                 tdims_s%k_start:tdims_s%k_end)                   &
                              ! OUT : SO2 LBC
, dms_lbc(lenrim(fld_type_p,halo_type_single),                    &
                 tdims_s%k_start:tdims_s%k_end)                   &
                              ! OUT : DMS LBC
, so4_aitken_lbc(lenrim(fld_type_p,halo_type_single),             &
                 tdims_s%k_start:tdims_s%k_end)                   &
                              ! OUT : SO4_AITKEN LBC
, so4_accu_lbc(lenrim(fld_type_p,halo_type_single),               &
                      tdims_s%k_start:tdims_s%k_end)              &
                              ! OUT : SO4_ACCU LBC
, so4_diss_lbc(lenrim(fld_type_p,halo_type_single),               &
               tdims_s%k_start:tdims_s%k_end)                     &
                              ! OUT : SO4_DISS LBC
, nh3_lbc(lenrim(fld_type_p,halo_type_single),                    &
                 tdims_s%k_start:tdims_s%k_end)                   &
                              ! OUT : NH3 LBC
, soot_new_lbc(lenrim(fld_type_p,halo_type_single),               &
                      tdims_s%k_start:tdims_s%k_end)              &
                              ! OUT : SOOT_NEW LBC
, soot_agd_lbc(lenrim(fld_type_p,halo_type_single),               &
                      tdims_s%k_start:tdims_s%k_end)              &
                              ! OUT : SOOT_AGD LBC
, soot_cld_lbc(lenrim(fld_type_p,halo_type_single),               &
                      tdims_s%k_start:tdims_s%k_end)              &
                              ! OUT : SOOT_CLD LBC
, bmass_new_lbc(lenrim(fld_type_p,halo_type_single),              &
                       tdims_s%k_start:tdims_s%k_end)             &
                              ! OUT : BMASS_NEW LBC
, bmass_agd_lbc(lenrim(fld_type_p,halo_type_single),              &
                tdims_s%k_start:tdims_s%k_end)                    &
                              ! OUT : BMASS_AGD LBC
, bmass_cld_lbc(lenrim(fld_type_p,halo_type_single),              &
                tdims_s%k_start:tdims_s%k_end)                    &
                              ! OUT : BMASS_CLD LBC
, ocff_new_lbc(lenrim(fld_type_p,halo_type_single),               &
                      tdims_s%k_start:tdims_s%k_end)              &
                              ! OUT : OCFF_NEW LBC
, ocff_agd_lbc(lenrim(fld_type_p,halo_type_single),               &
               tdims_s%k_start:tdims_s%k_end)                     &
                              ! OUT : OCFF_AGD LBC
, ocff_cld_lbc(lenrim(fld_type_p,halo_type_single),               &
               tdims_s%k_start:tdims_s%k_end)                     &
                              ! OUT : OCFF_CLD LBC
, nitr_acc_lbc(lenrim(fld_type_p,halo_type_single),               &
               tdims_s%k_start:tdims_s%k_end)                     &
                              ! OUT : NITR_ACC LBC
, nitr_diss_lbc(lenrim(fld_type_p,halo_type_single),              &
                tdims_s%k_start:tdims_s%k_end)                    &
                              ! OUT : NITR_DISS LBC
, tracer_lbc(lenrim(fld_type_p,halo_type_extended),               &
             tdims_l%k_start:tdims_l%k_end,tr_lbc_vars)           &
                              ! OUT : Tracer LBCs
, tracer_ukca_lbc(lenrim(fld_type_p,halo_type_extended),          &
                  tdims_l%k_start:tdims_l%k_end,tr_lbc_ukca)
                              ! OUT : UKCA Tracer LBCs

INTEGER ::                                                        &
  icode                           ! Error code

CHARACTER(LEN=errormessagelength) ::                              &
  cmessage                        ! Error message

! Local variables

INTEGER ::                                                        &
  tracer                                                          &
           ! loop counter over tracer variables
, lbc_num  ! number of active lbcs



INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_ATMOS_LBCS'

INTEGER :: typemap(1)

! ---------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! Check the field about to be read is U_LBC - if not, something
! has gone wrong.
IF ( lookup(item_code,lookup_subscript_u)  /=  31002) THEN
  WRITE(umMessage,*) 'Expected to find U_LBC (31002) in LBC file '
  CALL umPrint(umMessage,src='read_atmos_lbcs')
  WRITE(umMessage,*) 'But found ',lookup(item_code,1),' instead.'
  CALL umPrint(umMessage,src='read_atmos_lbcs')
  icode=1
  cmessage='READ_ATMOS_LBCS : Found wrong field in LBC file'
  GO TO 9999
END IF

! Now read in all the fields

typemap(1)=other

CALL readflds(nftin,                                                           &
              1,                                                               &
              1,                                                               &
              lookup(:,lookup_subscript_u:),                                   &
              u_lbc,                                                           &
              fixhd,                                                           &
              1,                                                               &
              icode,                                                           &
              cmessage,                                                        &
              typemap=typemap)

IF (icode  >   0) THEN
  WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading U_LBC'
  CALL umPrint(umMessage,src='read_atmos_lbcs')
  GO TO 9999
END IF

CALL readflds(nftin,                                                           &
              1,                                                               &
              1,                                                               &
              lookup(:,lookup_subscript_v:),                                   &
              v_lbc,                                                           &
              fixhd,                                                           &
              1,                                                               &
              icode,                                                           &
              cmessage,                                                        &
              typemap=typemap)

IF (icode  >   0) THEN
  WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading V_LBC'
  CALL umPrint(umMessage,src='read_atmos_lbcs')
  GO TO 9999
END IF

CALL readflds(nftin,                                                           &
              1,                                                               &
              1,                                                               &
              lookup(:,lookup_subscript_w:),                                   &
              w_lbc,                                                           &
              fixhd,                                                           &
              1,                                                               &
              icode,                                                           &
              cmessage,                                                        &
              typemap=typemap)

IF (icode  >   0) THEN
  WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading W_LBC'
  CALL umPrint(umMessage,src='read_atmos_lbcs')
  GO TO 9999
END IF

CALL readflds(nftin,                                                           &
              1,                                                               &
              1,                                                               &
              lookup(:,lookup_subscript_rho:),                                 &
              rho_lbc,fixhd,                                                   &
              1,                                                               &
              icode,                                                           &
              cmessage,                                                        &
              typemap=typemap)

IF (icode  >   0) THEN
  WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading RHO_LBC'
  CALL umPrint(umMessage,src='read_atmos_lbcs')
  GO TO 9999
END IF

CALL readflds(nftin,                                                           &
              1,                                                               &
              1,                                                               &
              lookup(:,lookup_subscript_theta:),                               &
              theta_lbc,                                                       &
              fixhd,                                                           &
              1,                                                               &
              icode,                                                           &
              cmessage,                                                        &
              typemap=typemap)

IF (icode  >   0) THEN
  WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading THETA_LBC'
  CALL umPrint(umMessage,src='read_atmos_lbcs')
  GO TO 9999
END IF

CALL readflds(nftin,                                                           &
              1,                                                               &
              1,                                                               &
              lookup(:,lookup_subscript_q:),                                   &
              q_lbc,                                                           &
              fixhd,                                                           &
              1,                                                               &
              icode,                                                           &
              cmessage,                                                        &
              typemap=typemap)

IF (icode  >   0) THEN
  WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading Q_LBC'
  CALL umPrint(umMessage,src='read_atmos_lbcs')
  GO TO 9999
END IF

CALL readflds(nftin,                                                           &
              1,                                                               &
              1,                                                               &
              lookup(:,lookup_subscript_qcl:),                                 &
              qcl_lbc,                                                         &
              fixhd,                                                           &
              1,                                                               &
              icode,                                                           &
              cmessage,                                                        &
              typemap=typemap)

IF (icode  >   0) THEN
  WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading QCL_LBC'
  CALL umPrint(umMessage,src='read_atmos_lbcs')
  GO TO 9999
END IF

CALL readflds(nftin,                                                           &
              1,                                                               &
              1,                                                               &
              lookup(:,lookup_subscript_qcf:),                                 &
              qcf_lbc,                                                         &
              fixhd,                                                           &
              1,                                                               &
              icode,                                                           &
              cmessage,                                                        &
              typemap=typemap)

IF (icode  >   0) THEN
  WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading QCF_LBC'
  CALL umPrint(umMessage,src='read_atmos_lbcs')
  GO TO 9999
END IF


CALL readflds(nftin,                                                           &
              1,                                                               &
              1,                                                               &
              lookup(:,lookup_subscript_exner:),                               &
              exner_lbc,                                                       &
              fixhd,                                                           &
              1,                                                               &
              icode,                                                           &
              cmessage,                                                        &
              typemap=typemap)

IF (icode  >   0) THEN
  WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading EXNER_LBC'
  CALL umPrint(umMessage,src='read_atmos_lbcs')
  GO TO 9999
END IF

! treat advected winds as optional LBCs
! using the interpolated winds logical switch
IF (.NOT. L_int_uvw_lbc) THEN

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lookup_subscript_u_adv:),                             &
                u_adv_lbc,                                                     &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading U_ADV_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lookup_subscript_v_adv:),                             &
                v_adv_lbc,                                                     &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading V_ADV_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lookup_subscript_w_adv:),                             &
                w_adv_lbc,fixhd,                                               &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading W_ADV_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

END IF  ! on L_int_uvw_lbc

! ----------------------------------------------------------------
! Read ice crystal (qcf2) lbcs if present in lbc file
! ----------------------------------------------------------------
!recover lbc_num variable from item_bounda_mod module
!need one less than value stored as no orography field here
lbc_num = lbc_num_stored - 1


IF (L_mcr_qcf2_lbc) THEN

  ! qcf2 is in input lbc file
  ! Increment number of lbcs
  lbc_num = lbc_num + 1

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lbc_num:),                                            &
                qcf2_lbc,                                                      &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading QCF2_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

ELSE IF (L_mcr_qcf2) THEN

  ! qcf2 prognostic is active but qcf2 lbcs are not in input file
  ! set lbc array to zero
  qcf2_lbc(:,:) = 0.0

END IF  ! on L_mcr_qcf2, L_mcr_qcf2_lbc

! ----------------------------------------------------------------
! Read rain lbcs if present in lbc file
! ----------------------------------------------------------------

IF (L_mcr_qrain_lbc) THEN

  ! qrain is in input lbc file
  ! Increment number of lbcs
  lbc_num = lbc_num + 1

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lbc_num:),                                            &
                qrain_lbc,                                                     &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading QRAIN_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

ELSE IF (L_mcr_qrain) THEN

  ! qrain is active but qrain lbcs are not in input file
  ! set lbc array to zero
  qrain_lbc(:,:) = 0.0

END IF  ! on L_mcr_qrain, L_mcr_qrain_lbc

! ----------------------------------------------------------------
! Read graupel lbcs if present in lbc file
! ----------------------------------------------------------------

IF (L_mcr_qgraup_lbc) THEN

  ! qgraup is in input lbc file
  ! Increment number of lbcs
  lbc_num = lbc_num + 1

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lbc_num:),                                            &
                qgraup_lbc,                                                    &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading QGRAUP_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

ELSE IF (L_mcr_qgraup) THEN

  ! qgraup is active but qgraup lbcs are not in input file
  ! set lbc array to zero
  qgraup_lbc(:,:) = 0.0

END IF  ! on L_mcr_qgraup, L_mcr_qgraup_lbc

! ----------------------------------------------------------------
! Read cloud fraction lbcs present in lbc file
! ----------------------------------------------------------------

! Set cloud fraction lbcs
IF (L_pc2_lbc) THEN

  ! Cloud fractions variables are in input lbc file
  ! Increment number of lbcs
  lbc_num = lbc_num + 1

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lbc_num:),                                            &
                cf_bulk_lbc,                                                   &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading CF_BULK_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

  lbc_num = lbc_num + 1

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lbc_num:),                                            &
                cf_liquid_lbc,                                                 &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading CF_LIQUID_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

  lbc_num = lbc_num + 1

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lbc_num:),                                            &
                cf_frozen_lbc,                                                 &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading CF_FROZEN_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

ELSE IF (i_cld_vn == i_cld_pc2) THEN

  ! Cloud fraction prognostics are active but cloud fraction lbcs
  ! are not in input file so set the lbc arrays to zero.
  cf_bulk_lbc(:,:)   = 0.0
  cf_liquid_lbc(:,:) = 0.0
  cf_frozen_lbc(:,:) = 0.0

END IF  ! on i_cld_pc2, L_pc2_lbc

! ----------------------------------------------------------------
! Read murk aerosol lbcs if expected in lbc file
! ----------------------------------------------------------------

IF (L_murk_lbc) THEN   ! murk lbcs are in lbc file

  ! Increment number of lbcs
  lbc_num = lbc_num + 1

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lbc_num:),                                            &
                murk_lbc,                                                      &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading MURK_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

ELSE IF (L_murk) THEN

  ! murk prognostic is active but murk lbcs are not in input file
  ! set lbc array to zero
  murk_lbc(:,:) = 0.0

END IF ! on L_murk_lbc, L_murk

! ----------------------------------------------------------------
! Read dust lbcs if expected in lbc file
! ----------------------------------------------------------------

IF (L_dust_div1_lbc) THEN   ! dust lbcs are in lbc file

  ! Increment number of lbcs
  lbc_num = lbc_num + 1

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lbc_num:),                                            &
                dust_div1_lbc,                                                 &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading DUST_DIV1_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

ELSE IF (l_dust_div1) THEN

  ! dust prognostic is active but dust lbcs are not in input file
  ! set lbc array to zero
  dust_div1_lbc(:,:) = 0.0

END IF ! on L_dust_div1_lbc, L_dust_div1

!-----------------------------------------------------------------

IF (L_dust_div2_lbc) THEN   ! dust lbcs are in lbc file

  ! Increment number of lbcs
  lbc_num = lbc_num + 1

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lbc_num:),                                            &
                dust_div2_lbc,                                                 &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading DUST_DIV2_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

ELSE IF (l_dust_div2) THEN

  ! dust prognostic is active but dust lbcs are not in input file
  ! set lbc array to zero
  dust_div2_lbc(:,:) = 0.0

END IF ! on L_dust_div2_lbc, L_dust_div2

!-----------------------------------------------------------------

IF (L_dust_div3_lbc) THEN   ! dust lbcs are in lbc file

  ! Increment number of lbcs
  lbc_num = lbc_num + 1

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lbc_num:),                                            &
                dust_div3_lbc,                                                 &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading DUST_DIV3_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

ELSE IF (l_dust_div3) THEN

  ! dust prognostic is active but dust lbcs are not in input file
  ! set lbc array to zero
  dust_div3_lbc(:,:) = 0.0

END IF ! on L_dust_div3_lbc, L_dust_div3

!-----------------------------------------------------------------

IF (L_dust_div4_lbc) THEN   ! dust lbcs are in lbc file

  ! Increment number of lbcs
  lbc_num = lbc_num + 1

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lbc_num:),                                            &
                dust_div4_lbc,                                                 &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading DUST_DIV4_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

ELSE IF (l_dust_div4) THEN

  ! dust prognostic is active but dust lbcs are not in input file
  ! set lbc array to zero
  dust_div4_lbc(:,:) = 0.0

END IF ! on L_dust_div4_lbc, L_dust_div4

!-----------------------------------------------------------------

IF (L_dust_div5_lbc) THEN   ! dust lbcs are in lbc file

  ! Increment number of lbcs
  lbc_num = lbc_num + 1

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lbc_num:),                                            &
                dust_div5_lbc,                                                 &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading DUST_DIV5_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

ELSE IF (l_dust_div5) THEN

  ! dust prognostic is active but dust lbcs are not in input file
  ! set lbc array to zero
  dust_div5_lbc(:,:) = 0.0

END IF ! on L_dust_div5_lbc, L_dust_div5

!-----------------------------------------------------------------

IF (L_dust_div6_lbc) THEN   ! dust lbcs are in lbc file

  ! Increment number of lbcs
  lbc_num = lbc_num + 1

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lbc_num:),                                            &
                dust_div6_lbc,                                                 &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading DUST_DIV6_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

ELSE IF (l_dust_div6) THEN

  ! dust prognostic is active but dust lbcs are not in input file
  ! set lbc array to zero
  dust_div6_lbc(:,:) = 0.0

END IF ! on L_dust_div6_lbc, L_dust_div6

! ----------------------------------------------------------------
! Read SO2 lbcs if expected in lbc file
! ----------------------------------------------------------------

IF (L_so2_lbc) THEN   ! so2 lbcs are in lbc file

  ! Increment number of lbcs
  lbc_num = lbc_num + 1

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lbc_num:),                                            &
                so2_lbc,                                                       &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading SO2_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

ELSE IF (l_so2) THEN

  ! so2 prognostic is active but so2 lbcs are not in input file
  ! set lbc array to zero
  so2_lbc(:,:) = 0.0

END IF ! on L_so2_lbc, L_so2

! ----------------------------------------------------------------
! Read DMS lbcs if expected in lbc file
! ----------------------------------------------------------------

IF (L_dms_lbc) THEN   ! dms lbcs are in lbc file

  ! Increment number of lbcs
  lbc_num = lbc_num + 1

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lbc_num:),                                            &
                dms_lbc,                                                       &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading DMS_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

ELSE IF (l_dms) THEN

  ! dms prognostic is active but dms lbcs are not in input file
  ! set lbc array to zero
  dms_lbc(:,:) = 0.0

END IF ! on L_dms_lbc, L_dms

! ----------------------------------------------------------------
! Read SO4_AITKEN lbcs if expected in lbc file
! ----------------------------------------------------------------

IF (L_so4_aitken_lbc) THEN   ! so4_aitken lbcs are in lbc file

  ! Increment number of lbcs
  lbc_num = lbc_num + 1

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lbc_num:),                                            &
                so4_aitken_lbc,                                                &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading SO4_AITKEN_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

ELSE IF (l_so4_aitken) THEN

  ! so4_aitken prognostic is active but so4_aitken lbcs
  ! are not in input file.
  ! Set lbc array to zero
  so4_aitken_lbc(:,:) = 0.0

END IF ! on L_so4_aitken_lbc, L_so4_aitken

! ----------------------------------------------------------------
! Read SO4_ACCU lbcs if expected in lbc file
! ----------------------------------------------------------------

IF (L_so4_accu_lbc) THEN   ! so4_accu lbcs are in lbc file

  ! Increment number of lbcs
  lbc_num = lbc_num + 1

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lbc_num:),                                            &
                so4_accu_lbc,                                                  &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading SO4_ACCU_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

ELSE IF (l_so4_accu) THEN

  ! so4_accu prognostic is active but so4_accu lbcs
  ! are not in input file.
  ! Set lbc array to zero
  so4_accu_lbc(:,:) = 0.0

END IF ! on L_so4_accu_lbc, L_so4_accu

! ----------------------------------------------------------------
! Read SO4_DISS lbcs if expected in lbc file
! ----------------------------------------------------------------

IF (L_so4_diss_lbc) THEN   ! so4_diss lbcs are in lbc file

  ! Increment number of lbcs
  lbc_num = lbc_num + 1

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lbc_num:),                                            &
                so4_diss_lbc,                                                  &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading SO4_DISS_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

ELSE IF (l_so4_diss) THEN

  ! so4_diss prognostic is active but so4_diss lbcs
  ! are not in input file.
  ! Set lbc array to zero
  so4_diss_lbc(:,:) = 0.0

END IF ! on L_so4_diss_lbc, L_so4_diss

! ----------------------------------------------------------------
! Read NH3 lbcs if expected in lbc file
! ----------------------------------------------------------------

IF (L_nh3_lbc) THEN   ! nh3 lbcs are in lbc file

  ! Increment number of lbcs
  lbc_num = lbc_num + 1

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lbc_num:),                                            &
                nh3_lbc,                                                       &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading NH3_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

ELSE IF (l_nh3) THEN

  ! nh3 prognostic is active but nh3 lbcs
  ! are not in input file.
  ! Set lbc array to zero
  nh3_lbc(:,:) = 0.0

END IF ! on L_nh3_lbc, L_nh3

! ----------------------------------------------------------------
! Read soot lbcs if expected in lbc file
! ----------------------------------------------------------------

IF (L_soot_new_lbc) THEN   ! soot_new lbcs are in lbc file

  ! Increment number of lbcs
  lbc_num = lbc_num + 1

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lbc_num:),                                            &
                soot_new_lbc,                                                  &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading SOOT_NEW_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

ELSE IF (l_soot_new) THEN

  ! soot_new prognostic is active but soot_new lbcs
  ! are not in input file.
  ! Set lbc array to zero
  soot_new_lbc(:,:) = 0.0

END IF ! on L_soot_new, L_soot_new

!-----------------------------------------------------------------

IF (L_soot_agd_lbc) THEN   ! soot_agd lbcs are in lbc file

  ! Increment number of lbcs
  lbc_num = lbc_num + 1

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lbc_num:),                                            &
                soot_agd_lbc,                                                  &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading SOOT_AGD_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

ELSE IF (l_soot_agd) THEN

  ! soot_agd prognostic is active but soot_agd lbcs
  ! are not in input file.
  ! Set lbc array to zero
  soot_agd_lbc(:,:) = 0.0

END IF ! on L_soot_agd, L_soot_agd

!-----------------------------------------------------------------

IF (L_soot_cld_lbc) THEN   ! soot_cld lbcs are in lbc file

  ! Increment number of lbcs
  lbc_num = lbc_num + 1

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lbc_num:),                                            &
                soot_cld_lbc,                                                  &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading SOOT_CLD_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

ELSE IF (l_soot_cld) THEN

  ! soot_cld prognostic is active but soot_cld lbcs
  ! are not in input file.
  ! Set lbc array to zero
  soot_cld_lbc(:,:) = 0.0

END IF ! on L_soot_cld, L_soot_cld

! ----------------------------------------------------------------
! Read biomass lbcs if expected in lbc file
! ----------------------------------------------------------------

IF (L_bmass_new_lbc) THEN   ! bmass_new lbcs are in lbc file

  ! Increment number of lbcs
  lbc_num = lbc_num + 1

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lbc_num:),                                            &
                bmass_new_lbc,                                                 &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading BMASS_NEW_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

ELSE IF (l_bmass_new) THEN

  ! bmass_new prognostic is active but bmass_new lbcs
  ! are not in input file.
  ! Set lbc array to zero
  bmass_new_lbc(:,:) = 0.0

END IF ! on L_bmass_new, L_bmass_new

!-----------------------------------------------------------------

IF (L_bmass_agd_lbc) THEN   ! bmass_agd lbcs are in lbc file

  ! Increment number of lbcs
  lbc_num = lbc_num + 1

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lbc_num:),                                            &
                bmass_agd_lbc,                                                 &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading BMASS_AGD_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

ELSE IF (l_bmass_agd) THEN

  ! bmass_agd prognostic is active but bmass_agd lbcs
  ! are not in input file.
  ! Set lbc array to zero
  bmass_agd_lbc(:,:) = 0.0

END IF ! on L_bmass_agd, L_bmass_agd

!-----------------------------------------------------------------

IF (L_bmass_cld_lbc) THEN   ! bmass_cld lbcs are in lbc file

  ! Increment number of lbcs
  lbc_num = lbc_num + 1

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lbc_num:),                                            &
                bmass_cld_lbc,                                                 &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading BMASS_CLD_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

ELSE IF (l_bmass_cld) THEN

  ! bmass_cld prognostic is active but bmass_cld lbcs
  ! are not in input file.
  ! Set lbc array to zero
  bmass_cld_lbc(:,:) = 0.0

END IF ! on L_bmass_cld, L_bmass_cld

! ----------------------------------------------------------------
! Read fossil fuel aerosol lbcs if expected in lbc file
! ----------------------------------------------------------------

IF (L_ocff_new_lbc) THEN   ! ocff_new lbcs are in lbc file

  ! Increment number of lbcs
  lbc_num = lbc_num + 1

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lbc_num:),                                            &
                ocff_new_lbc,                                                  &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading OCFF_NEW_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

ELSE IF (l_ocff_new) THEN

  ! ocff_new prognostic is active but ocff_new lbcs
  ! are not in input file.
  ! Set lbc array to zero
  ocff_new_lbc(:,:) = 0.0

END IF ! on L_ocff_new, L_ocff_new

!-----------------------------------------------------------------

IF (L_ocff_agd_lbc) THEN   ! ocff_agd lbcs are in lbc file

  ! Increment number of lbcs
  lbc_num = lbc_num + 1

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lbc_num:),                                            &
                ocff_agd_lbc,                                                  &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading OCFF_AGD_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

ELSE IF (l_ocff_agd) THEN

  ! ocff_agd prognostic is active but ocff_agd lbcs
  ! are not in input file.
  ! Set lbc array to zero
  ocff_agd_lbc(:,:) = 0.0

END IF ! on L_ocff_agd, L_ocff_agd

!-----------------------------------------------------------------

IF (L_ocff_cld_lbc) THEN   ! ocff_cld lbcs are in lbc file

  ! Increment number of lbcs
  lbc_num = lbc_num + 1

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lbc_num:),                                            &
                ocff_cld_lbc,                                                  &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading OCFF_CLD_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

ELSE IF (l_ocff_cld) THEN

  ! ocff_cld prognostic is active but ocff_cld lbcs
  ! are not in input file.
  ! Set lbc array to zero
  ocff_cld_lbc(:,:) = 0.0

END IF ! on L_ocff_cld, L_ocff_cld

! ----------------------------------------------------------------
! Read nitrate aerosol lbcs if expected in lbc file
! ----------------------------------------------------------------

IF (L_nitr_acc_lbc) THEN   ! nitr_acc lbcs are in lbc file

  ! Increment number of lbcs
  lbc_num = lbc_num + 1

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lbc_num:),                                            &
                nitr_acc_lbc,                                                  &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading NITR_ACC_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

ELSE IF (l_nitr_acc) THEN

  ! nitr_acc prognostic is active but nitr_acc lbcs
  ! are not in input file.
  ! Set lbc array to zero
  nitr_acc_lbc(:,:) = 0.0

END IF ! on L_nitr_acc_lbc, L_nitr_acc

!-----------------------------------------------------------------

IF (L_nitr_diss_lbc) THEN   ! nitr_diss lbcs are in lbc file

  ! Increment number of lbcs
  lbc_num = lbc_num + 1

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lbc_num:),                                            &
                nitr_diss_lbc,                                                 &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading NITR_DISS_LBC'
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

ELSE IF (l_nitr_diss) THEN

  ! nitr_diss prognostic is active but nitr_diss lbcs
  ! are not in input file.
  ! Set lbc array to zero
  nitr_diss_lbc(:,:) = 0.0

END IF ! on L_nitr_diss_lbc, L_nitr_diss

! ----------------------------------------------------------------
! Read tracer lbcs if expected in lbc file
! ----------------------------------------------------------------

DO tracer=1,tr_lbc_vars

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lbc_num+tracer:),                                     &
                tracer_lbc(1:,1:,tracer:),                                     &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'READ_ATMOS_LBCS : Problem reading Tracer LBC ', &
        tracer
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

END DO  ! tracer

!Increment number of lbcs
lbc_num = lbc_num + tr_lbc_vars

! ----------------------------------------------------------------
! Read UKCA tracer lbcs if expected in lbc file
! ----------------------------------------------------------------

DO tracer=1,tr_lbc_ukca

  CALL readflds(nftin,                                                         &
                1,                                                             &
                1,                                                             &
                lookup(:,lbc_num+tracer:),                                     &
                tracer_ukca_lbc(1:,1:,tracer:),                                &
                fixhd,                                                         &
                1,                                                             &
                icode,                                                         &
                cmessage,                                                      &
                typemap=typemap)

  IF (icode  >   0) THEN
    WRITE(umMessage,*) &
        'READ_ATMOS_LBCS : Problem reading UKCA Tracer LBC ',tracer
    CALL umPrint(umMessage,src='read_atmos_lbcs')
    GO TO 9999
  END IF

END DO  ! tracer


9999 CONTINUE

! DEPENDS ON: convert_lbcs
CALL convert_lbcs(theta_lbc, q_lbc, qcl_lbc, qcf_lbc,             &
                   qcf2_lbc, qrain_lbc, qgraup_lbc,               &
                   l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,         &
                   lenrim(fld_type_p,halo_type_extended) )

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE read_atmos_lbcs
