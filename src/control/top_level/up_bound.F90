! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine UP_BOUND
!
! Purpose: Updates the boundary conditions in D1 with data from disk.
!          Atmosphere: On the first call, the orography LBC is read
!                      in, together with the LBCs at the start
!                      (FIELD_LBC) and end (FIELD_LBC_TEND) of the
!                      update period. (The only exception to this is
!                      when the LBC_TEND data is to come from a
!                      different boundary file. In this case, the data
!                      from the first file that would normally be read
!                      into LBC is read into LBC_TEND instead. An extra
!                      call is then made to this routine in which the
!                      LBC_TEND data is copied into LBC, and data from
!                      the second boundary file is read into LBC_TEND.)
!                        At subsequent boundary updating steps, the
!                      value that was in FIELD_LBC_TEND) is copied
!                      to FIELD_LBC (this is now the LBC at the start
!                      of the new period) and the next record is read
!                      from disk into FIELD_LBC_TEND.
!          Other     : At the first step, two records are read from
!                      disk. The first record is stored in the dump,
!                      as is the tendency (the difference between the
!                      first and second records).
!                      At subsequent boundary updating sets, the next
!                      record is read, and the difference between this
!                      and the current LBC value is stored in the dump
!                      as the new tendency.
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Top Level

SUBROUTINE up_bound(i_ao,                                         &
                  icode,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE atm_fields_bounds_mod
USE atm_fields_mod,     ONLY: orog_lbc, u_lbc, v_lbc, w_lbc, rho_lbc,       &
                              theta_lbc, q_lbc, qcl_lbc, qcf_lbc,           &
                              qcf2_lbc, qrain_lbc, qgraup_lbc,              &
                              cf_bulk_lbc, cf_liquid_lbc, cf_frozen_lbc,    &
                              exner_lbc, u_adv_lbc, v_adv_lbc, w_adv_lbc,   &
                              murk_lbc, tracer_lbc, tracer_ukca_lbc,        &
                              dust_div1_lbc, dust_div2_lbc, dust_div3_lbc,  &
                              dust_div4_lbc, dust_div5_lbc, dust_div6_lbc,  &
                              so2_lbc, dms_lbc, so4_aitken_lbc,             &
                              so4_accu_lbc, so4_diss_lbc, nh3_lbc,          &
                              soot_new_lbc, soot_agd_lbc, soot_cld_lbc,     &
                              bmass_new_lbc, bmass_agd_lbc, bmass_cld_lbc,  &
                              ocff_new_lbc, ocff_agd_lbc, ocff_cld_lbc,     &
                              nitr_acc_lbc, nitr_diss_lbc, u_lbc_tend,      &
                              v_lbc_tend, w_lbc_tend, rho_lbc_tend,         &
                              theta_lbc_tend, q_lbc_tend, qcl_lbc_tend,     &
                              qcf_lbc_tend, qcf2_lbc_tend, qrain_lbc_tend,  &
                              qgraup_lbc_tend, cf_bulk_lbc_tend,            &
                              cf_liquid_lbc_tend, cf_frozen_lbc_tend,       &
                              exner_lbc_tend, u_adv_lbc_tend,               &
                              v_adv_lbc_tend, w_adv_lbc_tend,               &
                              murk_lbc_tend, tracer_lbc_tend,               &
                              tracer_ukca_lbc_tend, dust_div1_lbc_tend,     &
                              dust_div2_lbc_tend, dust_div3_lbc_tend,       &
                              dust_div4_lbc_tend, dust_div5_lbc_tend,       &
                              dust_div6_lbc_tend, so2_lbc_tend,             &
                              dms_lbc_tend, so4_aitken_lbc_tend,            &
                              so4_accu_lbc_tend, so4_diss_lbc_tend,         &
                              nh3_lbc_tend, soot_new_lbc_tend,              &
                              soot_agd_lbc_tend, soot_cld_lbc_tend,         &
                              bmass_new_lbc_tend, bmass_agd_lbc_tend,       &
                              bmass_cld_lbc_tend, ocff_new_lbc_tend,        &
                              ocff_agd_lbc_tend, ocff_cld_lbc_tend,         &
                              nitr_acc_lbc_tend, nitr_diss_lbc_tend
USE umPrintMgr
USE UM_ParParams
USE Control_Max_Sizes
USE rimtypes
USE lbc_mod
USE lookup_addresses
USE dust_parameters_mod, ONLY: l_dust,                            &
     l_dust_div1,       l_dust_div2,        l_dust_div3,          &
     l_dust_div4,       l_dust_div5,        l_dust_div6,          &
     l_dust_div1_lbc,   l_dust_div2_lbc,    l_dust_div3_lbc,      &
     l_dust_div4_lbc,   l_dust_div5_lbc,    l_dust_div6_lbc


USE run_aerosol_mod, ONLY:  &
     l_so2, l_dms, l_so4_aitken, l_so4_accu, l_so4_diss,     &
     l_soot_new, l_soot_agd, l_soot_cld, l_bmass_new,        &
     l_bmass_agd, l_bmass_cld, l_ocff_new, l_ocff_agd, l_ocff_cld,&
     l_nh3, l_nitr_acc, l_nitr_diss, l_so2_lbc, l_dms_lbc,   &
     l_so4_aitken_lbc, l_so4_accu_lbc, l_so4_diss_lbc, l_nh3_lbc, &
     l_soot_new_lbc,    l_soot_agd_lbc,     l_soot_cld_lbc,  &
     l_bmass_new_lbc,   l_bmass_agd_lbc,    l_bmass_cld_lbc, &
     l_ocff_new_lbc,    l_ocff_agd_lbc,     l_ocff_cld_lbc,  &
     l_nitr_acc_lbc,    l_nitr_diss_lbc

USE mphys_inputs_mod, ONLY:                                       &
     l_mcr_qcf2,        l_mcr_qrain,        l_mcr_qgraup,         &
     l_mcr_qcf2_lbc,    l_mcr_qrain_lbc,    l_mcr_qgraup_lbc
USE cloud_inputs_mod, ONLY: l_pc2_lbc
USE murk_inputs_mod, ONLY: l_murk, l_murk_lbc
USE lbc_read_data_mod, ONLY: albc_num, albc_swapstep, l_int_uvw_lbc
USE submodel_mod, ONLY: atmos_im
USE nlstcall_mod, ONLY: Num_ALBCs

USE nlsizes_namelist_mod, ONLY:                                        &
    a_len1_coldepc, a_len1_flddepc, a_len1_levdepc, a_len1_rowdepc,    &
    a_len2_coldepc, a_len2_flddepc, a_len2_levdepc, a_len2_lookup,     &
    a_len2_rowdepc, a_len_cfi1, a_len_cfi2, a_len_cfi3, a_len_extcnst, &
    a_len_inthd, a_len_realhd, len1_lbc_comp_lookup, len1_lookup,      &
    len_dumphist, len_fixhd, len_tot,                                  &
    model_levels, mpp_len1_lookup, n_cca_lev, n_obj_d1_max, sm_levels, &
    st_levels, tpps_ozone_levels, tr_lbc_ukca, tr_lbc_vars, tr_levels, &
    tr_ukca, tr_vars

USE file_manager, ONLY: get_file_unit_by_id

USE model_time_mod, ONLY: &
    bndary_offsetim, stepim
USE errormessagelength_mod, ONLY: errormessagelength

USE model_domain_mod, ONLY: model_type, mt_global

USE atm_boundary_headers_mod, ONLY: lookup_bounda, rim_stepsa,         & 
                                    lookup_comp_bounda, fixhd_bounda,  &
                                    nbound_lookup

USE readflds_mod, ONLY: readflds

USE d1_array_mod, ONLY: other

IMPLICIT NONE



INTEGER ::                                                        &
       i_ao,                                                      &
                     !  atmosphere/Ocean indicator
       icode         ! Error code = 0 Normal Exit
!                          !            > 0 Error Condition

CHARACTER(LEN=errormessagelength) ::                              &
       cmessage      ! Error message

!

! Local variables

INTEGER ::                                                        &
       i,j                                                        &
,      pretend_ts                                                 &
,      first_ts                                                   &
                            ! first timestep
,      last_ts                                                    &
                            ! last  timestep
,      steps_to_next_update                                       &
,      len_buf                                                    &
                            ! length of buffer for readflds
,      item_in_file         ! stash item code
LOGICAL ::                                                        &
       periodic,                                                  &
                          ! True if periodic lateral boundary data
       Between_ALBC_files ! True if about to swap atmos bndy files

INTEGER           :: steps_from_bdi_start ! Timesteps between
                                          ! start of bdi and
                                          ! start of run
LOGICAL, SAVE     :: first_atm_call = .TRUE.
CHARACTER (LEN=4) :: ch_lbc_time
INTEGER           :: lbc_unit

INTEGER :: typemap(1)

! --------------------------------------------------------------------

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UP_BOUND'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
icode=0
cmessage=' '

IF (model_type /= mt_global) THEN

  ! 1.1   See whether the LBC and LBC_TEND data is to come from separate
  !       boundary files, a situation that could arise if we're using two
  !       boundary files. If this is the case, the LBC data from the first
  !       file will initially be read into LBC_TEND. An extra call to this
  !       routine will then be made after the second boundary file has
  !       been opened. On this extra call, the LBC_TEND data will be
  !       copied to LBC, and the LBC_TEND data will be replaced with data
  !       from the second boundary file.

          ! Initialise:
  Between_ALBC_files = .FALSE.
  Steps_from_bdi_start = bndary_offsetim(atmos_im)

  IF (Num_ALBCS == 2)                                             &
    Between_ALBC_files = ALBC_num == 1 .AND.                      &
                         stepim(atmos_im) >= ALBC_SwapStep

  ! 1.2 Read atmosphere lateral boundary field, first step.

  IF (first_atm_call .AND. i_ao == 1) THEN

    lbc_unit = get_file_unit_by_id("lbc_input", handler="portio")

    ! 1.2.1 Read Orography data

    typemap(1)=other

    CALL readflds(lbc_unit,                                                    &
                  1,                                                           &
                  1,                                                           &
                  lookup_bounda,                                               &
                  orog_lbc(:),                                                 &
                  RESHAPE(fixhd_bounda(:,:),[SIZE(fixhd_bounda)]),             &
                  1,                                                           &
                  icode,                                                       &
                  cmessage,                                                    &
                  typemap=typemap)

    IF (icode  >   0) THEN
      WRITE(umMessage,*) 'UP_BOUND : Problem in READFLDS reading ',        &
                 'atmosphere orography'
      CALL umPrint(umMessage,src='up_bound')
      WRITE(umMessage,*) 'ICODE: ',icode
      CALL umPrint(umMessage,src='up_bound')
      WRITE(umMessage,*) 'CMESSAGE ',cmessage
      CALL umPrint(umMessage,src='up_bound')
      GO TO 9999
    END IF

    ! 1.2.2 Update LOOKUP_BOUNDA with the correct information for the
    !       current set of LOOKUP headers

    DO i=2,rim_lookupsa     ! Loop over lookup headers for the
                            ! first record (ignoring orog.)

      j=nbound_lookup(1)+i-2  ! The "real" lookup header number
                              ! from the LBC file

      lookup_bounda(lbyr,i)=lookup_comp_bounda(lbcc_lbyr,j)
      lookup_bounda(lbmon,i)=lookup_comp_bounda(lbcc_lbmon,j)
      lookup_bounda(lbdat,i)=lookup_comp_bounda(lbcc_lbdat,j)
      lookup_bounda(lbhr,i)=lookup_comp_bounda(lbcc_lbhr,j)
      lookup_bounda(lbmin,i)=lookup_comp_bounda(lbcc_lbmin,j)
      lookup_bounda(lbsec,i)=lookup_comp_bounda(lbcc_lbsec,j)
      lookup_bounda(lbegin,i)=lookup_comp_bounda(lbcc_lbegin,j)
      lookup_bounda(naddr,i)=lookup_comp_bounda(lbcc_naddr,j)

    END DO ! i

    IF (.NOT. Between_ALBC_files) THEN

      ! 1.2.3 Read in the LBCs for timestep 0

      ! DEPENDS ON: read_atmos_lbcs
      CALL read_atmos_lbcs(                                         &
        lenrima(1,1,rima_type_norm),                                &
        global_LENRIMA(1,1,rima_type_norm),                         &
        tr_lbc_vars,tr_levels,tr_lbc_ukca,                          &
        L_int_uvw_lbc,L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,        &
        L_mcr_qcf2_lbc, L_mcr_qrain_lbc,                            &
        L_mcr_qgraup_lbc, L_pc2_lbc, L_murk, L_murk_lbc,            &
        L_dust_div1, L_dust_div1_lbc, L_dust_div2, L_dust_div2_lbc, &
        L_dust_div3, L_dust_div3_lbc, L_dust_div4, L_dust_div4_lbc, &
        L_dust_div5, L_dust_div5_lbc, L_dust_div6, L_dust_div6_lbc, &
        L_so2, L_so2_lbc, L_dms, L_dms_lbc,                         &
        L_so4_aitken, L_so4_aitken_lbc, L_so4_accu, L_so4_accu_lbc, &
        L_so4_diss, L_so4_diss_lbc, L_nh3, L_nh3_lbc,               &
        L_soot_new, L_soot_new_lbc, L_soot_agd, L_soot_agd_lbc,     &
        L_soot_cld, L_soot_cld_lbc, L_bmass_new, L_bmass_new_lbc,   &
        L_bmass_agd, L_bmass_agd_lbc, L_bmass_cld, L_bmass_cld_lbc, &
        L_ocff_new, L_ocff_new_lbc,                                 &
        L_ocff_agd, L_ocff_agd_lbc, L_ocff_cld, L_ocff_cld_lbc,     &
        L_nitr_acc, L_nitr_acc_lbc, L_nitr_diss, L_nitr_diss_lbc,   &
        len1_lookup, len_fixhd, lbc_unit,                           &
        rim_lookupsa-1, lookup_bounda(1,2), fixhd_bounda,           &
        u_lbc, v_lbc, w_lbc, rho_lbc,                               &
        theta_lbc, q_lbc, qcl_lbc, qcf_lbc,                         &
        qcf2_lbc, qrain_lbc, qgraup_lbc,                            &
        cf_bulk_lbc, cf_liquid_lbc, cf_frozen_lbc,                  &
        exner_lbc, u_adv_lbc, v_adv_lbc,                            &
        w_adv_lbc, murk_lbc,                                        &
        dust_div1_lbc, dust_div2_lbc, dust_div3_lbc,                &
        dust_div4_lbc, dust_div5_lbc, dust_div6_lbc,                &
        so2_lbc, dms_lbc, so4_aitken_lbc,                           &
        so4_accu_lbc,so4_diss_lbc, nh3_lbc,                         &
        soot_new_lbc, soot_agd_lbc, soot_cld_lbc,                   &
        bmass_new_lbc, bmass_agd_lbc, bmass_cld_lbc,                &
        ocff_new_lbc, ocff_agd_lbc, ocff_cld_lbc,                   &
        nitr_acc_lbc, nitr_diss_lbc,                                &
        tracer_lbc,  tracer_ukca_lbc,                               &
        icode,cmessage)

      IF (icode  >   0) THEN
        WRITE(umMessage,*) 'Problem with READ_ATMOS_LBCS reading initial ', &
                   'boundary data for timestep 0'
        CALL umPrint(umMessage,src='up_bound')
        WRITE(umMessage,*) 'ICODE= ',icode
        CALL umPrint(umMessage,src='up_bound')
        WRITE(umMessage,*) 'CMESSAGE= ',cmessage
        CALL umPrint(umMessage,src='up_bound')
        GO TO 9999
      END IF

      IF (PrintStatus >= PrStatus_Normal) THEN
        WRITE (ch_lbc_time(1:2),'(I2.2)') lookup_bounda(lbhr,2)
        WRITE (ch_lbc_time(3:4),'(I2.2)') lookup_bounda(lbmin,2)
        WRITE(umMessage,*)                                                 &
        'Up_Bound: Timestep ',stepim(atmos_im),' : LBCs read in for ', &
        ch_lbc_time,'Z ',lookup_bounda(lbdat,2),'/',                &
        lookup_bounda(lbmon,2),'/',lookup_bounda(lbyr,2)
        CALL umPrint(umMessage,src='up_bound')
      END IF

      ! 1.2.4  Increment lookup header to start of next update period

      nbound_lookup(1)=nbound_lookup(1)+(rim_lookupsa-1)

      DO i=2,rim_lookupsa     ! Loop over lookup headers for the
                              ! first record (ignoring orog.)

        j=nbound_lookup(1)+i-2  ! The "real" lookup header number
                                ! from the LBC file

        lookup_bounda(lbyr,i)=lookup_comp_bounda(lbcc_lbyr,j)
        lookup_bounda(lbmon,i)=lookup_comp_bounda(lbcc_lbmon,j)
        lookup_bounda(lbdat,i)=lookup_comp_bounda(lbcc_lbdat,j)
        lookup_bounda(lbhr,i)=lookup_comp_bounda(lbcc_lbhr,j)
        lookup_bounda(lbmin,i)=lookup_comp_bounda(lbcc_lbmin,j)
        lookup_bounda(lbsec,i)=lookup_comp_bounda(lbcc_lbsec,j)
        lookup_bounda(lbegin,i)=lookup_comp_bounda(lbcc_lbegin,j)
        lookup_bounda(naddr,i)=lookup_comp_bounda(lbcc_naddr,j)

      END DO ! i

    END IF ! (.NOT.Between_ALBC_files)

    ! 1.2.5 Read in the LBC_TEND for the end of the update period

    ! DEPENDS ON: read_atmos_lbcs
    CALL read_atmos_lbcs(                                           &
      lenrima(1,1,rima_type_norm),                                  &
      global_LENRIMA(1,1,rima_type_norm),                           &
      tr_lbc_vars,tr_levels,tr_lbc_ukca,                            &
      L_int_uvw_lbc,L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,          &
      L_mcr_qcf2_lbc, L_mcr_qrain_lbc,                              &
      L_mcr_qgraup_lbc, L_pc2_lbc, L_murk, L_murk_lbc,              &
      L_dust_div1, L_dust_div1_lbc, L_dust_div2, L_dust_div2_lbc,   &
      L_dust_div3, L_dust_div3_lbc, L_dust_div4, L_dust_div4_lbc,   &
      L_dust_div5, L_dust_div5_lbc, L_dust_div6, L_dust_div6_lbc,   &
      L_so2,L_so2_lbc,L_dms,L_dms_lbc,L_so4_aitken,L_so4_aitken_lbc,&
      L_so4_accu, L_so4_accu_lbc, L_so4_diss, L_so4_diss_lbc,       &
      L_nh3, L_nh3_lbc,                                             &
      L_soot_new, L_soot_new_lbc, L_soot_agd, L_soot_agd_lbc,       &
      L_soot_cld, L_soot_cld_lbc, L_bmass_new, L_bmass_new_lbc,     &
      L_bmass_agd, L_bmass_agd_lbc, L_bmass_cld, L_bmass_cld_lbc,   &
      L_ocff_new, L_ocff_new_lbc,                                   &
      L_ocff_agd, L_ocff_agd_lbc, L_ocff_cld, L_ocff_cld_lbc,       &
      L_nitr_acc, L_nitr_acc_lbc, L_nitr_diss, L_nitr_diss_lbc,     &
      len1_lookup, len_fixhd, lbc_unit,                             &
      rim_lookupsa-1, lookup_bounda(1,2), fixhd_bounda,             &
      u_lbc_tend, v_lbc_tend, w_lbc_tend,                           &
      rho_lbc_tend, theta_lbc_tend, q_lbc_tend,                     &
      qcl_lbc_tend, qcf_lbc_tend,                                   &
      qcf2_lbc_tend, qrain_lbc_tend,                                &
      qgraup_lbc_tend, cf_bulk_lbc_tend,                            &
      cf_liquid_lbc_tend, cf_frozen_lbc_tend,                       &
      exner_lbc_tend,                                               &
      u_adv_lbc_tend, v_adv_lbc_tend,                               &
      w_adv_lbc_tend, murk_lbc_tend,                                &
      dust_div1_lbc_tend, dust_div2_lbc_tend,                       &
      dust_div3_lbc_tend, dust_div4_lbc_tend,                       &
      dust_div5_lbc_tend, dust_div6_lbc_tend,                       &
      so2_lbc_tend, dms_lbc_tend,                                   &
      so4_aitken_lbc_tend,                                          &
      so4_accu_lbc_tend, so4_diss_lbc_tend,                         &
      nh3_lbc_tend, soot_new_lbc_tend,                              &
      soot_agd_lbc_tend, soot_cld_lbc_tend,                         &
      bmass_new_lbc_tend, bmass_agd_lbc_tend,                       &
      bmass_cld_lbc_tend, ocff_new_lbc_tend,                        &
      ocff_agd_lbc_tend, ocff_cld_lbc_tend,                         &
      nitr_acc_lbc_tend, nitr_diss_lbc_tend,                        &
      tracer_lbc_tend, tracer_ukca_lbc_tend,                        &
      icode,cmessage)

    IF (icode  >   0) THEN
      WRITE(umMessage,*) 'Problem with READ_ATMOS_LBCS reading tendency ',  &
                 'boundary data for timestep 0'
      CALL umPrint(umMessage,src='up_bound')
      WRITE(umMessage,*) 'ICODE= ',icode
      CALL umPrint(umMessage,src='up_bound')
      WRITE(umMessage,*) 'CMESSAGE= ',cmessage
      CALL umPrint(umMessage,src='up_bound')
      GO TO 9999
    END IF

    IF (PrintStatus >= PrStatus_Normal) THEN
      WRITE (ch_lbc_time(1:2),'(I2.2)') lookup_bounda(lbhr,2)
      WRITE (ch_lbc_time(3:4),'(I2.2)') lookup_bounda(lbmin,2)
      WRITE(umMessage,*)                                                   &
      'Up_Bound: Timestep ',stepim(atmos_im),' : LBCs read in for ', &
      ch_lbc_time,'Z ',lookup_bounda(lbdat,2),'/',                  &
      lookup_bounda(lbmon,2),'/',lookup_bounda(lbyr,2)
      CALL umPrint(umMessage,src='up_bound')
    END IF

    ! Increment the lookup header, ready for the next data

    nbound_lookup(1)=nbound_lookup(1)+(rim_lookupsa-1)

    ! Check to see if we need to increment the LBCs to the correct timestep
    ! (If between boundary files, this step will be done in an additional
    ! call to this routine.)

    IF (MOD(Steps_from_bdi_start+stepim(atmos_im),rim_stepsa) /= 0  &
        .AND. .NOT. Between_ALBC_files) THEN

      IF ((stepim(atmos_im) < rim_stepsa) .AND.                     &
          (stepim(atmos_im) /= 0)) THEN  

        ! First timestep between two LBC update intervals
        ! We need to increment the LBC to the correct timestep

        first_ts = 0
        last_ts  = stepim(atmos_im)-1

      ELSE

        ! The current timestep falls between two LBC update
        ! intervals (this is probably a continuation run).
        ! We need to increment the LBC to the correct timestep

        first_ts = 0
        last_ts  = MOD(Steps_from_bdi_start + stepim(atmos_im),     &
                       rim_stepsa) - 1

      END IF

      DO pretend_ts = first_ts, last_ts

        steps_to_next_update=rim_stepsa-pretend_ts

        ! DEPENDS ON: increment_atmos_lbcs
        CALL increment_atmos_lbcs(                                  &
          steps_to_next_update,                                     &
          lenrima(1,1,rima_type_norm),                              &
          tr_lbc_vars,tr_levels, tr_lbc_ukca,                       &
          L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc,        &
          L_pc2_lbc, L_murk_lbc, L_int_uvw_lbc,                     &
          L_dust_div1_lbc,L_dust_div2_lbc,                          &
          L_dust_div3_lbc,L_dust_div4_lbc,                          &
          L_dust_div5_lbc,L_dust_div6_lbc,                          &
          L_so2_lbc,L_dms_lbc,L_so4_aitken_lbc,                     &
          L_so4_accu_lbc,L_so4_diss_lbc,                            &
          L_nh3_lbc,L_soot_new_lbc,L_soot_agd_lbc,                  &
          L_soot_cld_lbc,L_bmass_new_lbc,                           &
          L_bmass_agd_lbc,L_bmass_cld_lbc,                          &
          L_ocff_new_lbc,L_ocff_agd_lbc,L_ocff_cld_lbc,             &
          L_nitr_acc_lbc, L_nitr_diss_lbc,                          &
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
          cf_bulk_lbc  ,cf_bulk_lbc_tend  ,                         &
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
          so2_lbc,so2_lbc_tend,                                     &
          dms_lbc,dms_lbc_tend,                                     &
          so4_aitken_lbc,so4_aitken_lbc_tend,                       &
          so4_accu_lbc,so4_accu_lbc_tend,                           &
          so4_diss_lbc,so4_diss_lbc_tend,                           &
          nh3_lbc,nh3_lbc_tend,                                     &
          soot_new_lbc,soot_new_lbc_tend,                           &
          soot_agd_lbc,soot_agd_lbc_tend,                           &
          soot_cld_lbc,soot_cld_lbc_tend,                           &
          bmass_new_lbc,bmass_new_lbc_tend,                         &
          bmass_agd_lbc,bmass_agd_lbc_tend,                         &
          bmass_cld_lbc,bmass_cld_lbc_tend,                         &
          ocff_new_lbc,ocff_new_lbc_tend,                           &
          ocff_agd_lbc,ocff_agd_lbc_tend,                           &
          ocff_cld_lbc,ocff_cld_lbc_tend,                           &
          nitr_acc_lbc,nitr_acc_lbc_tend,                           &
          nitr_diss_lbc,nitr_diss_lbc_tend,                         &
          tracer_lbc,tracer_lbc_tend,                               &
          tracer_ukca_lbc,tracer_ukca_lbc_tend,                     &
          rim_stepsa,                                               &
          icode)

        IF (icode  >   0) THEN
          WRITE(umMessage,*) 'Failure in INCREMENT_ATMOS_LBCS while ',      &
                     'attempting to set LBCS for first timestep'
          CALL umPrint(umMessage,src='up_bound')
          GO TO 9999
        END IF

      END DO  ! pretend_ts

    END IF  ! IF (MOD(stepim(atmos_im),RIM_STEPSA)  /=  0 .AND.
           !     .NOT.Between_ALBC_files)


  END IF  ! (first_atm_call .AND. I_AO == 1)


  ! 2.1 Read atmosphere lateral boundary fields, general update step


  IF (.NOT. first_atm_call .AND. i_ao == 1) THEN

    lbc_unit = get_file_unit_by_id("lbc_input", handler="portio")

    IF (MOD(steps_from_bdi_start + stepim(atmos_im),                &
                                                 ! LBC update step
            rim_stepsa) == 0 ) THEN

      IF (nbound_lookup(1)  >=  fixhd_bounda(152,1)) THEN
        icode=11
        cmessage='UP_BOUND : Reached end of atmosphere LBC file'
        GO TO 9999
      END IF

      ! 2.1.1 Copy data from LBC_TEND to LBC - this is the new LBC data
      !       starting point for the new LBC update period

      ! DEPENDS ON: copy_atmos_lbcs
      CALL copy_atmos_lbcs(                                         &
        lenrima(1,1,rima_type_norm),                                &
        tr_lbc_vars,tr_levels,                                      &
        tr_lbc_ukca,                                                &
        L_mcr_qcf2_lbc, L_mcr_qrain_lbc,                            &
        L_mcr_qgraup_lbc, L_pc2_lbc, L_murk_lbc,                    &
        L_dust_div1_lbc,L_dust_div2_lbc,                            &
        L_dust_div3_lbc,L_dust_div4_lbc,                            &
        L_dust_div5_lbc,L_dust_div6_lbc,                            &
        L_so2_lbc,L_dms_lbc,L_so4_aitken_lbc,                       &
        L_so4_accu_lbc,L_so4_diss_lbc,                              &
        L_nh3_lbc,L_soot_new_lbc,L_soot_agd_lbc,                    &
        L_soot_cld_lbc,L_bmass_new_lbc,                             &
        L_bmass_agd_lbc,L_bmass_cld_lbc,                            &
        L_ocff_new_lbc,L_ocff_agd_lbc,L_ocff_cld_lbc,               &
        L_nitr_acc_lbc, L_nitr_diss_lbc,                            &
        u_lbc, u_lbc_tend,                                          &
        v_lbc, v_lbc_tend,                                          &
        w_lbc, w_lbc_tend,                                          &
        rho_lbc, rho_lbc_tend,                                      &
        theta_lbc, theta_lbc_tend,                                  &
        q_lbc, q_lbc_tend,                                          &
        qcl_lbc, qcl_lbc_tend,                                      &
        qcf_lbc, qcf_lbc_tend,                                      &
        qcf2_lbc, qcf2_lbc_tend,                                    &
        qrain_lbc, qrain_lbc_tend,                                  &
        qgraup_lbc, qgraup_lbc_tend,                                &
        cf_bulk_lbc  , cf_bulk_lbc_tend  ,                          &
        cf_liquid_lbc, cf_liquid_lbc_tend,                          &
        cf_frozen_lbc, cf_frozen_lbc_tend,                          &
        exner_lbc, exner_lbc_tend,                                  &
        u_adv_lbc, u_adv_lbc_tend,                                  &
        v_adv_lbc, v_adv_lbc_tend,                                  &
        w_adv_lbc, w_adv_lbc_tend,                                  &
        murk_lbc, murk_lbc_tend,                                    &
        dust_div1_lbc,dust_div1_lbc_tend,                           &
        dust_div2_lbc,dust_div2_lbc_tend,                           &
        dust_div3_lbc,dust_div3_lbc_tend,                           &
        dust_div4_lbc,dust_div4_lbc_tend,                           &
        dust_div5_lbc,dust_div5_lbc_tend,                           &
        dust_div6_lbc,dust_div6_lbc_tend,                           &
        so2_lbc,so2_lbc_tend,                                       &
        dms_lbc,dms_lbc_tend,                                       &
        so4_aitken_lbc,so4_aitken_lbc_tend,                         &
        so4_accu_lbc,so4_accu_lbc_tend,                             &
        so4_diss_lbc,so4_diss_lbc_tend,                             &
        nh3_lbc,nh3_lbc_tend,                                       &
        soot_new_lbc,soot_new_lbc_tend,                             &
        soot_agd_lbc,soot_agd_lbc_tend,                             &
        soot_cld_lbc,soot_cld_lbc_tend,                             &
        bmass_new_lbc,bmass_new_lbc_tend,                           &
        bmass_agd_lbc,bmass_agd_lbc_tend,                           &
        bmass_cld_lbc,bmass_cld_lbc_tend,                           &
        ocff_new_lbc,ocff_new_lbc_tend,                             &
        ocff_agd_lbc,ocff_agd_lbc_tend,                             &
        ocff_cld_lbc,ocff_cld_lbc_tend,                             &
        nitr_acc_lbc, nitr_acc_lbc_tend,                            &
        nitr_diss_lbc, nitr_diss_lbc_tend,                          &
        tracer_lbc, tracer_lbc_tend,                                &
        tracer_ukca_lbc, tracer_ukca_lbc_tend                       &
        )

      ! 2.1.2 Update LOOKUP_BOUNDA with the correct information for the
      !       current set of LOOKUP headers

      DO i=2,rim_lookupsa     ! Loop over lookup headers for the
                              ! first record (ignoring orog.)

        j=nbound_lookup(1)+i-2  ! The "real" lookup header number
                                ! from the LBC file

        lookup_bounda(lbyr,i)=lookup_comp_bounda(lbcc_lbyr,j)
        lookup_bounda(lbmon,i)=lookup_comp_bounda(lbcc_lbmon,j)
        lookup_bounda(lbdat,i)=lookup_comp_bounda(lbcc_lbdat,j)
        lookup_bounda(lbhr,i)=lookup_comp_bounda(lbcc_lbhr,j)
        lookup_bounda(lbmin,i)=lookup_comp_bounda(lbcc_lbmin,j)
        lookup_bounda(lbsec,i)=lookup_comp_bounda(lbcc_lbsec,j)
        lookup_bounda(lbegin,i)=lookup_comp_bounda(lbcc_lbegin,j)
        lookup_bounda(naddr,i)=lookup_comp_bounda(lbcc_naddr,j)

      END DO ! i

      ! 2.1.3 Read in the new set of LBC_TENDs

      ! DEPENDS ON: read_atmos_lbcs
      CALL read_atmos_lbcs(                                         &
        lenrima(1,1,rima_type_norm),                                &
        global_LENRIMA(1,1,rima_type_norm),                         &
        tr_lbc_vars,tr_levels, tr_lbc_ukca,                         &
        L_int_uvw_lbc,L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,        &
        L_mcr_qcf2_lbc, L_mcr_qrain_lbc,                            &
        L_mcr_qgraup_lbc, L_pc2_lbc, L_murk, L_murk_lbc,            &
        L_dust_div1, L_dust_div1_lbc, L_dust_div2, L_dust_div2_lbc, &
        L_dust_div3, L_dust_div3_lbc, L_dust_div4, L_dust_div4_lbc, &
        L_dust_div5, L_dust_div5_lbc, L_dust_div6, L_dust_div6_lbc, &
        L_so2, L_so2_lbc, L_dms, L_dms_lbc,                         &
        L_so4_aitken, L_so4_aitken_lbc, L_so4_accu, L_so4_accu_lbc, &
        L_so4_diss, L_so4_diss_lbc, L_nh3, L_nh3_lbc,               &
        L_soot_new, L_soot_new_lbc, L_soot_agd, L_soot_agd_lbc,     &
        L_soot_cld, L_soot_cld_lbc, L_bmass_new, L_bmass_new_lbc,   &
        L_bmass_agd, L_bmass_agd_lbc, L_bmass_cld, L_bmass_cld_lbc, &
        L_ocff_new, L_ocff_new_lbc,                                 &
        L_ocff_agd, L_ocff_agd_lbc, L_ocff_cld, L_ocff_cld_lbc,     &
        L_nitr_acc, L_nitr_acc_lbc, L_nitr_diss, L_nitr_diss_lbc,   &
        len1_lookup, len_fixhd, lbc_unit,                           &
        rim_lookupsa-1, lookup_bounda(1,2), fixhd_bounda,           &
        u_lbc_tend, v_lbc_tend, w_lbc_tend,                         &
        rho_lbc_tend, theta_lbc_tend, q_lbc_tend,                   &
        qcl_lbc_tend, qcf_lbc_tend,                                 &
        qcf2_lbc_tend, qrain_lbc_tend,                              &
        qgraup_lbc_tend,cf_bulk_lbc_tend,                           &
        cf_liquid_lbc_tend, cf_frozen_lbc_tend,                     &
        exner_lbc_tend,                                             &
        u_adv_lbc_tend, v_adv_lbc_tend,                             &
        w_adv_lbc_tend, murk_lbc_tend,                              &
        dust_div1_lbc_tend, dust_div2_lbc_tend,                     &
        dust_div3_lbc_tend, dust_div4_lbc_tend,                     &
        dust_div5_lbc_tend, dust_div6_lbc_tend,                     &
        so2_lbc_tend, dms_lbc_tend,                                 &
        so4_aitken_lbc_tend,                                        &
        so4_accu_lbc_tend, so4_diss_lbc_tend,                       &
        nh3_lbc_tend,soot_new_lbc_tend,                             &
        soot_agd_lbc_tend, soot_cld_lbc_tend,                       &
        bmass_new_lbc_tend, bmass_agd_lbc_tend,                     &
        bmass_cld_lbc_tend, ocff_new_lbc_tend,                      &
        ocff_agd_lbc_tend, ocff_cld_lbc_tend,                       &
        nitr_acc_lbc_tend, nitr_diss_lbc_tend,                      &
        tracer_lbc_tend,tracer_ukca_lbc_tend,                       &
        icode,cmessage)

      IF (icode  >   0) THEN
        WRITE(umMessage,*) 'Problem with READ_ATMOS_LBCS reading ',         &
                   'LBC_TENDs for general update step.'
        CALL umPrint(umMessage,src='up_bound')
        WRITE(umMessage,*) 'ICODE= ',icode
        CALL umPrint(umMessage,src='up_bound')
        WRITE(umMessage,*) 'CMESSAGE= ',cmessage
        CALL umPrint(umMessage,src='up_bound')
        GO TO 9999
      END IF

      IF (PrintStatus >= PrStatus_Normal) THEN
        WRITE (ch_lbc_time(1:2),'(I2.2)') lookup_bounda(lbhr,2)
        WRITE (ch_lbc_time(3:4),'(I2.2)') lookup_bounda(lbmin,2)
        WRITE(umMessage,*)                                                 &
        'Up_Bound: Timestep ',stepim(atmos_im),' : LBCs read in for ', &
        ch_lbc_time,'Z ',lookup_bounda(lbdat,2),'/',                &
        lookup_bounda(lbmon,2),'/',lookup_bounda(lbyr,2)
        CALL umPrint(umMessage,src='up_bound')
      END IF

      IF (MOD(stepim(atmos_im),rim_stepsa)  /=  0 .AND.             &
          Between_ALBC_files) THEN

        ! The current timestep falls between two LBC update
        ! intervals (this is probably a continuation run).
        ! We need to increment the LBC to the correct timestep

        DO pretend_ts=0,MOD(stepim(atmos_im),rim_stepsa)-1

          steps_to_next_update=rim_stepsa-pretend_ts

          ! DEPENDS ON: increment_atmos_lbcs
          CALL increment_atmos_lbcs(                                &
            steps_to_next_update,                                   &
            lenrima(1,1,rima_type_norm),                            &
            tr_lbc_vars,tr_levels, tr_lbc_ukca,                     &
            L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc,      &
            L_pc2_lbc, L_murk_lbc, L_int_uvw_lbc,                   &
            L_dust_div1_lbc,L_dust_div2_lbc,                        &
            L_dust_div3_lbc,L_dust_div4_lbc,                        &
            L_dust_div5_lbc,L_dust_div6_lbc,                        &
            L_so2_lbc,L_dms_lbc,L_so4_aitken_lbc,                   &
            L_so4_accu_lbc,L_so4_diss_lbc,                          &
            L_nh3_lbc,L_soot_new_lbc,L_soot_agd_lbc,                &
            L_soot_cld_lbc,L_bmass_new_lbc,                         &
            L_bmass_agd_lbc,L_bmass_cld_lbc,                        &
            L_ocff_new_lbc,L_ocff_agd_lbc,L_ocff_cld_lbc,           &
            L_nitr_acc_lbc, L_nitr_diss_lbc,                        &
            u_lbc,u_lbc_tend,                                       &
            v_lbc,v_lbc_tend,                                       &
            w_lbc,w_lbc_tend,                                       &
            rho_lbc,rho_lbc_tend,                                   &
            theta_lbc,theta_lbc_tend,                               &
            q_lbc,q_lbc_tend,                                       &
            qcl_lbc,qcl_lbc_tend,                                   &
            qcf_lbc,qcf_lbc_tend,                                   &
            qcf2_lbc,qcf2_lbc_tend,                                 &
            qrain_lbc,qrain_lbc_tend,                               &
            qgraup_lbc,qgraup_lbc_tend,                             &
            cf_bulk_lbc  ,cf_bulk_lbc_tend  ,                       &
            cf_liquid_lbc,cf_liquid_lbc_tend,                       &
            cf_frozen_lbc,cf_frozen_lbc_tend,                       &
            exner_lbc,exner_lbc_tend,                               &
            u_adv_lbc,u_adv_lbc_tend,                               &
            v_adv_lbc,v_adv_lbc_tend,                               &
            w_adv_lbc,w_adv_lbc_tend,                               &
            murk_lbc, murk_lbc_tend,                                &
            dust_div1_lbc,dust_div1_lbc_tend,                       &
            dust_div2_lbc,dust_div2_lbc_tend,                       &
            dust_div3_lbc,dust_div3_lbc_tend,                       &
            dust_div4_lbc,dust_div4_lbc_tend,                       &
            dust_div5_lbc,dust_div5_lbc_tend,                       &
            dust_div6_lbc,dust_div6_lbc_tend,                       &
            so2_lbc,so2_lbc_tend,                                   &
            dms_lbc,dms_lbc_tend,                                   &
            so4_aitken_lbc,so4_aitken_lbc_tend,                     &
            so4_accu_lbc,so4_accu_lbc_tend,                         &
            so4_diss_lbc,so4_diss_lbc_tend,                         &
            nh3_lbc,nh3_lbc_tend,                                   &
            soot_new_lbc,soot_new_lbc_tend,                         &
            soot_agd_lbc,soot_agd_lbc_tend,                         &
            soot_cld_lbc,soot_cld_lbc_tend,                         &
            bmass_new_lbc,bmass_new_lbc_tend,                       &
            bmass_agd_lbc,bmass_agd_lbc_tend,                       &
            bmass_cld_lbc,bmass_cld_lbc_tend,                       &
            ocff_new_lbc,ocff_new_lbc_tend,                         &
            ocff_agd_lbc,ocff_agd_lbc_tend,                         &
            ocff_cld_lbc,ocff_cld_lbc_tend,                         &
            nitr_acc_lbc,nitr_acc_lbc_tend,                         &
            nitr_diss_lbc,nitr_diss_lbc_tend,                       &
            tracer_lbc,tracer_lbc_tend,                             &
            tracer_ukca_lbc,tracer_ukca_lbc_tend,                   &
            rim_stepsa,                                             &
            icode)

          IF (icode  >   0) THEN
            WRITE(umMessage,*) 'Failure in INCREMENT_ATMOS_LBCS while ',    &
                       'attempting to set LBCS for first timestep'
            CALL umPrint(umMessage,src='up_bound')
            GO TO 9999
          END IF

        END DO  ! pretend_ts

      END IF  ! IF (MOD(stepim(atmos_im),RIM_STEPSA)  /=  0 .AND.
             !     Between_ALBC_files)

      ! 2.1.4  Increment lookup header to start of next update period

      nbound_lookup(1)=nbound_lookup(1)+(rim_lookupsa-1)

    END IF ! IF it's an LBC update step

  ELSE

    first_atm_call = .FALSE.

  END IF !  (.NOT.first_atm_call .AND. I_AO == 1)

END IF  ! .NOT. GLOBAL

9999 CONTINUE
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE up_bound
