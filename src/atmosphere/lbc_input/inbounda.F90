! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Subroutine INBOUNDA
!
!   Purpose : Takes as input,the code defining whether updates of
!    boundary data are required.  The physical files required are
!    identified, and the headers lookup tables are read in.
!    Reads the update intervals from the boundary datasets.
!    Where the update interval is in months or years, the check will be
!    made daily.
!
!   Programming standard; Unified Model Documentation Paper No. 3
!   version no. 1, dated 15/01/90
!
!   Logical components covered : C720
!
!   System task : C7
!
!   Documentation : Unified Model Documentation Paper No C7
!
!
!
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: LBC Input

SUBROUTINE inbounda(                                              &
   a_len1_levdepcda,a_len2_levdepcda,                             &
   a_len1_rowdepcda,a_len2_rowdepcda,                             &
   a_len1_coldepcda,a_len2_coldepcda)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE io, ONLY: setpos
USE io_constants, ONLY: ioNoDelete, ioOpenReadOnly
USE atm_fields_bounds_mod
USE atm_d1_indices_mod, ONLY: jetatheta, jetarho
USE check_iostat_mod
USE dump_headers_mod, ONLY: a_fixhd, a_inthd, a_realhd, a_levdepc,&
                            a_rowdepc, a_coldepc
USE submodel_mod, ONLY: internal_id_max, atmos_im
USE nlstgen_mod, ONLY: steps_per_periodim, secs_per_periodim
USE near_equal_real_mod, ONLY: near_equal_real
USE ereport_mod, ONLY: ereport
USE umPrintMgr
USE UM_ParVars
USE Field_Types, ONLY: fld_type_p
USE Control_Max_Sizes
USE lbc_mod
USE lookup_addresses
USE dust_parameters_mod, ONLY:             l_dust,                &
     l_dust_div1,       l_dust_div2,       l_dust_div3,           &
     l_dust_div4,       l_dust_div5,       l_dust_div6,           &
     l_dust_div1_lbc,   l_dust_div2_lbc,   l_dust_div3_lbc,       &
     l_dust_div4_lbc,   l_dust_div5_lbc,   l_dust_div6_lbc

USE run_aerosol_mod,  ONLY:                                  &
     l_so2_lbc,        l_dms_lbc,       l_so4_aitken_lbc,    &
     l_so4_accu_lbc,   l_so4_diss_lbc,  l_nh3_lbc,           &
     l_soot_new_lbc,   l_soot_agd_lbc,  l_soot_cld_lbc,      &
     l_bmass_new_lbc,  l_bmass_agd_lbc, l_bmass_cld_lbc,     &
     l_ocff_new_lbc,   l_ocff_agd_lbc,  l_ocff_cld_lbc,      &
     l_nitr_acc_lbc,   l_nitr_diss_lbc

USE model_file, ONLY: model_file_open, model_file_close
USE idealise_run_mod, ONLY: l_force_lbc, l_fixed_lbcs

USE lbc_read_data_mod, ONLY: albc_num, albc_swapstep,             &
                             albc2_starttime_steps,               &
                             current_lbc_step, rimweightsa
USE mphys_inputs_mod, ONLY:                                       &
     l_mcr_qcf2,       l_mcr_qcf2_lbc,    l_mcr_qgraup,           &
     l_mcr_qgraup_lbc, l_mcr_qrain_lbc,   l_mcr_qrain

USE murk_inputs_mod,  ONLY: l_murk, l_murk_lbc

USE nlstcall_mod, ONLY: Num_ALBCs, lcal360 

USE nlcfiles_namelist_mod, ONLY: &
    lbc_file_1 => alabcin1, lbc_file_2 => alabcin2
USE filenamelength_mod, ONLY: filenamelength
USE file_manager, ONLY: assign_file_unit, release_file_unit, &
                        get_file_unit_by_id

USE item_bounda_mod

USE nlsizes_namelist_mod, ONLY:                                        &
    a_len1_coldepc, a_len1_flddepc, a_len1_levdepc, a_len1_rowdepc,    &
    a_len2_coldepc, a_len2_flddepc, a_len2_levdepc, a_len2_lookup,     &
    a_len2_rowdepc, a_len_cfi1, a_len_cfi2, a_len_cfi3, a_len_extcnst, &
    a_len_inthd, a_len_realhd, len1_lbc_comp_lookup, len1_lookup,      &
    len_dumphist, len_fixhd, model_levels,                             &
    mpp_len1_lookup, n_cca_lev, sm_levels, st_levels,                  &
    tpps_ozone_levels, tr_lbc_ukca, tr_lbc_vars, tr_ukca, tr_vars

USE conversions_mod, ONLY: isec_per_day

USE model_time_mod, ONLY:                                                &
    basis_time_days, basis_time_secs, bndary_offsetim, boundary_stepsim, &
    i_day, i_hour, i_minute, i_month, i_second, i_year, stepim
USE errormessagelength_mod, ONLY: errormessagelength

USE model_domain_mod, ONLY: model_type, mt_cyclic_lam, mt_bi_cyclic_lam

USE atm_boundary_headers_mod, ONLY: fixhd_bounda, inthd_bounda,        &
                                    lookup_bounda, realhd_bounda,      &
                                    lookup_comp_bounda, rim_stepsa,    &
                                    nbound_lookup


IMPLICIT NONE


INTEGER :: a_len1_levdepcda   ! IN : copy of A_LEN1_LEVDEPC
INTEGER :: a_len2_levdepcda   ! IN : copy of A_LEN2_LEVDEPC
INTEGER :: a_len1_rowdepcda   ! IN : copy of A_LEN1_ROWDEPC
INTEGER :: a_len2_rowdepcda   ! IN : copy of A_LEN2_ROWDEPC
INTEGER :: a_len1_coldepcda   ! IN : copy of A_LEN1_COLDEPC
INTEGER :: a_len2_coldepcda   ! IN : copy of A_LEN2_COLDEPC

!    Local variables

INTEGER ::                                                        &
        i,                                                        &
        j,                                                        &
        j1,                                                       &
        lbc_num,                                                  &
        start_block,                                              &
        im_index,                                                 &
                             ! Internal model index
        elapsed_days,                                             &
                             ! Days since basis time
        elapsed_secs,                                             &
                             ! Secs since basis time
        current_time_days,                                        &
                             ! No. of days to current time
        current_time_secs,                                        &
                             ! No. of secs-in-day to current time
        days_to_data_start,                                       &
                             ! Days  to start of boundary data
        secs_to_data_start,                                       &
                             ! Secs  to start of boundary data
        days_to_data_end,                                         &
                             ! Days  to end   of boundary data
        secs_to_data_end,                                         &
                             ! Secs  to end   of boundary data
        steps_to_data_start,                                      &
                             ! Steps to start of boundary data
        steps_to_data_end,                                        &
                             ! Steps to end   of boundary data
        rim_stepsa_OLD,                                           &
                             ! Data interval for last bndy file
        steps_to_bdi_start,                                       &
                             ! Steps to start/end of
        steps_to_bdi_end,                                         &
                             ! current boundary data interval
        basis_to_data_start_steps,                                &
                                   ! Steps from basis time to
                                   ! start of boundary data
        item_bounda(rim_lookupsa)  ! Boundary updatable item list

LOGICAL :: this_albc_for_bdi_end ! True IF same boundary file to
                                 ! be used for end of current
                                 ! boundary data interval

REAL ::                                                              &
      a_levdepc_bo(a_len1_levdepcda,a_len2_levdepcda),               &
      a_rowdepc_bo(MAX(a_len1_rowdepcda,1),MAX(a_len2_rowdepcda,1)), &
      a_coldepc_bo(MAX(a_len1_coldepcda,1),MAX(a_len2_coldepcda,1))

INTEGER :: full_lookup_bounda(len1_lookup,bound_lookupsa)

INTEGER, PARAMETER :: dummy =1


INTEGER             :: ErrorStatus      ! Return code
CHARACTER (LEN=errormessagelength) :: cmessage         ! Error message
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'INBOUNDA'

INTEGER                        :: lbc_unit
CHARACTER (LEN=filenamelength) :: filename

LOGICAL, SAVE :: L_FirstCall = .TRUE.

INTEGER :: k

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

!       Internal Structure

ErrorStatus=0
cmessage=' '

IF (L_FirstCall) THEN

  !   1.0 Initialise variables, lest undefined for some
  !       choices of boundary updating, eg. -DEF,GLOBAL but not DEF,FLOOR,
  !       as used in section 2 to set boundary_stepsim for atmos
  rim_stepsa = 0
  !  Initialise bndary_offsetim & boundary_stepsim in model_time_mod
  DO i=1,internal_id_max
    bndary_offsetim(i) = 0
    boundary_stepsim(i) = 0
  END DO

END IF ! (L_FirstCall)

!   1.1 Update interval for lateral boundaries for atmosphere
!       Read headers and test whether boundary updating required

IF (l_fixed_lbcs .OR. l_force_lbc .OR.                          &
    model_type == mt_bi_cyclic_lam .OR.                       &
    model_type == mt_cyclic_lam) THEN
  rim_stepsa=0
ELSE

  !       Open input boundary file and read headers

  nbound_lookup(1)=1

  IF (ALBC_num == 2) THEN

    lbc_unit = get_file_unit_by_id("lbc_input", handler="portio")
    CALL model_file_close(lbc_unit, lbc_file_1,               &
        delete=ioNoDelete, error=errorstatus)
    IF (ErrorStatus /= 0) THEN
      WRITE (cmessage,'(A)') 'Failure closing first boundary file.'
      CALL Ereport(RoutineName,ErrorStatus,cmessage)
    END IF
    CALL release_file_unit(lbc_unit, handler="portio")
    filename = lbc_file_2
  ELSE
    filename = lbc_file_1
  END IF

  CALL assign_file_unit(filename, lbc_unit,                         &
                        handler="portio", id="lbc_input")

  CALL model_file_open(lbc_unit, filename,                          &
         read_write=ioOpenReadOnly, error=ErrorStatus)

  !       Read in fixed header to get array dimensions
  ! DEPENDS ON: read_flh
  CALL read_flh(lbc_unit,fixhd_bounda(1,1),                          &
                       len_fixhd,ErrorStatus,cmessage)
  IF (ErrorStatus >  0) THEN
    WRITE (cmessage,'(A)')                                        &
          'INBOUNDA : Error in READ_FLH for BOUNDA(1,1)'

    CALL Ereport(RoutineName,ErrorStatus,cmessage)
  END IF

  !       Check for negative dimensions
  IF (fixhd_bounda(101,1) <= 0) fixhd_bounda(101,1)=1
  IF (fixhd_bounda(106,1) <= 0) fixhd_bounda(106,1)=1
  IF (fixhd_bounda(111,1) <= 0) fixhd_bounda(111,1)=1
  IF (fixhd_bounda(112,1) <= 0) fixhd_bounda(112,1)=1
  IF (fixhd_bounda(116,1) <= 0) fixhd_bounda(116,1)=1
  IF (fixhd_bounda(117,1) <= 0) fixhd_bounda(117,1)=1
  IF (fixhd_bounda(121,1) <= 0) fixhd_bounda(121,1)=1
  IF (fixhd_bounda(122,1) <= 0) fixhd_bounda(122,1)=1
  IF (fixhd_bounda(151,1) <= 0) fixhd_bounda(151,1)=1
  IF (fixhd_bounda(152,1) <= 0) fixhd_bounda(152,1)=1
  IF (fixhd_bounda(161,1) <= 0) fixhd_bounda(161,1)=1

  !       Check IF sufficient space allocated for LOOKUP table
  IF (fixhd_bounda(152,1) >  bound_lookupsa) THEN
    CALL umPrint(' INBOUNDA; not enough space for LBC lookup headers.')
    CALL umPrint('           try increasing value specified in namelist')
    CALL umPrint(' nrim_timesa  in gui  LBC Related Options')
    WRITE(cmessage,'(A)')                                         &
        'INBOUNDA: Insufficient space for Lookup Table'
    Errorstatus = 2

    CALL Ereport(RoutineName,ErrorStatus,cmessage)
  END IF


  CALL setpos (lbc_unit,0,ErrorStatus)
  IF (ErrorStatus >  0) THEN
    CALL umPrint( 'INBOUNDA: Problem with SETPOS for BOUNDA(1,1)')
    WRITE(umMessage,'(A,I8,A,I8)')                                       &
        'ErrorStatus ',ErrorStatus,' lbc_unit ',lbc_unit
    CALL umPrint(umMessage)
    WRITE(cmessage,'(A)') 'Problem with Setpos for Bounda(1,1)'

    CALL Ereport(RoutineName,Errorstatus,cmessage)
  END IF

  ! DEPENDS ON: readhead
  CALL readhead(lbc_unit,                                         &
                fixhd_bounda(1,1),len_fixhd,                      &
                inthd_bounda(1,1),fixhd_bounda(101,1),            &
                realhd_bounda(1,1),fixhd_bounda(106,1),           &
                a_levdepc_bo(1,1),                                &
                fixhd_bounda(111,1),fixhd_bounda(112,1),          &
                a_rowdepc_bo(1,1),                                &
                fixhd_bounda(116,1),fixhd_bounda(117,1),          &
                a_coldepc_bo(1,1),                                &
                fixhd_bounda(121,1),fixhd_bounda(122,1),          &
                dummy,dummy,dummy,                                &
                dummy,dummy,                                      &
                dummy,dummy,                                      &
                dummy,dummy,                                      &
                dummy,dummy,                                      &
                dummy,dummy,                                      &
                full_lookup_bounda,                               &
                fixhd_bounda(151,1),fixhd_bounda(152,1),          &
                fixhd_bounda(161,1),                              &
                start_block,ErrorStatus,cmessage)

  IF (ErrorStatus >  0) THEN
    CALL umPrint( 'INBOUNDA: Problem with READHEAD for BOUNDA(1,1)')
    WRITE(umMessage,'(A,I8,A,A)')'ErrorStatus ',ErrorStatus,' cmessage ',&
        cmessage
    CALL umPrint(umMessage)
    WRITE(cmessage,'(A)') 'Problem with READHEAD for BOUNDA(1,1)'

    CALL Ereport(RoutineName,errorStatus,cmessage)
  END IF

  ! Copy the first set of headers into LOOKUP_BOUNDA
  DO i=1,rim_lookupsa
    DO j=1,len1_lookup
      lookup_bounda(j,i)=full_lookup_bounda(j,i)
    END DO ! j
  END DO ! i

  ! Copy all the varying items from the LOOKUP into COMP_LOOKUP_BOUNDA
  DO i=1,bound_lookupsa

    lookup_comp_bounda(lbcc_lbyr,i)=full_lookup_bounda(lbyr,i)
    lookup_comp_bounda(lbcc_lbmon,i)=full_lookup_bounda(lbmon,i)
    lookup_comp_bounda(lbcc_lbdat,i)=full_lookup_bounda(lbdat,i)
    lookup_comp_bounda(lbcc_lbhr,i)=full_lookup_bounda(lbhr,i)
    lookup_comp_bounda(lbcc_lbmin,i)=full_lookup_bounda(lbmin,i)
    lookup_comp_bounda(lbcc_lbsec,i)=full_lookup_bounda(lbsec,i)
    lookup_comp_bounda(lbcc_lbegin,i)=full_lookup_bounda(lbegin,i)
    lookup_comp_bounda(lbcc_naddr,i)=full_lookup_bounda(naddr,i)

  END DO ! i

  ! Check validity of headers

  ! Integer headers
  IF (fixhd_bounda(100,1)  >   0) THEN

    IF (inthd_bounda(6,1) /= glsize(1,fld_type_p)) THEN
      CALL umPrint( 'LBC Integer Header Mismatch:')
      WRITE(umMessage,'(A,I8)')                                          &
          'ROW_LENGTH from INTHD: ',inthd_bounda(6,1)
      CALL umPrint(umMessage)
      WRITE(umMessage,'(A,I8)')                                          &
          'Model ROW_LENGTH: ',glsize(1,fld_type_p)
      CALL umPrint(umMessage)

      ErrorStatus=3
      WRITE(cmessage,'(A)')'Integer header (ROW_LENGTH) mismatch'

      CALL Ereport(RoutineName,ErrorStatus,cmessage)
    END IF

    IF (inthd_bounda(7,1) /= glsize(2,fld_type_p)) THEN
      CALL umPrint( 'LBC Integer Header Mismatch:')
      WRITE(umMessage,'(A,I8)') 'Number of rows from INTHD: ',           &
          inthd_bounda(7,1)
      CALL umPrint(umMessage)
      WRITE(umMessage,'(A,I8)') 'Model number of rows: ',                &
          glsize(2,fld_type_p)
      CALL umPrint(umMessage)

      ErrorStatus=4
      WRITE(cmessage,'(A)')'Integer header (N_ROWS) mismatch'

      CALL Ereport(Routinename,Errorstatus,cmessage)
    END IF

    IF ( inthd_bounda(17,1) /= a_inthd(17) ) THEN
      CALL umPrint( 'LBC Integer Header Mismatch:')
      WRITE(umMessage,'(A,I8)') 'LBC   : Height Generator Method : ',            &
                  inthd_bounda(17,1)
      CALL umPrint(umMessage)
      WRITE(umMessage,'(A,I8)') 'Model : Height Generator Method : ',            &
                  a_inthd(17)
      CALL umPrint(umMessage)

      ErrorStatus=5
      WRITE (cmessage,'(A)') 'INBOUNDA : Mis-match in height '//         &
                         'generator method.'

      CALL Ereport(RoutineName,ErrorStatus,cmessage)
    END IF

    IF ( inthd_bounda(24,1) /= a_inthd(24) ) THEN
      CALL umPrint( 'LBC Integer Header Mismatch:')
      WRITE(umMessage,'(A,I8)')                                          &
          'LBC  : First rho level with constant height: ',               &
          inthd_bounda(24,1)
      CALL umPrint(umMessage)
      WRITE(umMessage,'(A,I8)')                                          &
          'Model: First rho level with constant height: ',a_inthd(24)
      CALL umPrint(umMessage)

      ErrorStatus=6
      WRITE (cmessage,'(A)') 'INBOUNDA : Mis-match in height '//         &
                         'generator method.'

      CALL Ereport(RoutineName,ErrorStatus,cmessage)
    END IF

  END IF ! IF (fixhd_bounda(100,1)  >   0)

  ! Real constants

  IF (fixhd_bounda(105,1)  >   0) THEN

    ! Check real headers in LBC file against those in dump to
    ! 32 bit accuracy
    DO j=1,6
      IF (near_equal_real(realhd_bounda(j,1),a_realhd(j))) THEN
        WRITE(umMessage,'(A,I8,A)')                                      &
            'LBC Real Header mismatch at position ',j,':'
        CALL umPrint(umMessage)
        WRITE(umMessage,'(A,E16.8)')                                     &
            'Value from LBC file is ',realhd_bounda(j,1)
        CALL umPrint(umMessage)
        WRITE(umMessage,'(A,E16.8)')                                     &
            'Value from model dump is ',a_realhd(j)
        CALL umPrint(umMessage)

        ErrorStatus=7
        WRITE(cmessage,'(A)') 'INBOUNDA : Real header mismatch'

        CALL Ereport(RoutineName,ErrorStatus,cmessage)
      END IF
    END DO ! j

    IF (near_equal_real(realhd_bounda(16,1),a_realhd(16))) THEN
      CALL umPrint( 'LBC Real Header mismatch at position 16 :')
      WRITE(umMessage,'(A,E16.8)')'LBC file : Height at top of model ',  &
          realhd_bounda(16,1)
      CALL umPrint(umMessage)
      WRITE(umMessage,'(A,E16.8)') 'Model    : Height at top of model ', &
          a_realhd(16)
      CALL umPrint(umMessage)

      ErrorStatus=8
      WRITE(cmessage,'(A)') 'INBOUNDA : Model Height mismatch'

      CALL Ereport(RoutineName,ErrorStatus,cmessage)
    END IF

  END IF ! IF (fixhd_bounda(105,1)  >   0)

  ! Level dependent constants

  IF (fixhd_bounda(110,1) > 0) THEN

    ! ---------------------------------
    ! Check eta values for theta levels
    ! ---------------------------------

    DO j=1, model_levels+1
      IF ( near_equal_real( a_levdepc_bo(j,1),                           &
                 a_levdepc(jetatheta+j-1) ) ) THEN

        CALL umPrint( 'LBC Level Dependent Constants Mismatch')
        WRITE(umMessage,'(A,A,I8)')                                      &
            'Eta values for theta levels mismatch for ','level ',j
        CALL umPrint(umMessage)
        WRITE(umMessage,'(A,E16.8)') 'Value from LBC file   : ',         &
            a_levdepc_bo(j,1)
        CALL umPrint(umMessage)
        WRITE(umMessage,'(A,E16.8)') 'Value from model dump : ',         &
            a_levdepc(jetatheta+j-1)
        CALL umPrint(umMessage)
        ErrorStatus=9
        WRITE(cmessage,'(A)')                                            &
            'INBOUNDA : Level dependent constant mismatch'

        CALL Ereport(RoutineName,ErrorStatus,cmessage)

      END IF
    END DO ! j

    ! -------------------------------
    ! Check eta values for rho levels
    ! -------------------------------
    DO j=1, model_levels
      IF ( near_equal_real( a_levdepc_bo(j,2),                           &
          a_levdepc(jetarho+j-1) ) ) THEN

        CALL umPrint( 'LBC Level Dependent Constants Mismatch')
        WRITE(umMessage,'(A,A,I8)')                                      &
            'Eta values for rho levels mismatch for ','level ',j
        CALL umPrint(umMessage)
        WRITE(umMessage,'(A,E16.8)') 'Value from LBC file   : ',         &
            a_levdepc_bo(j,2)
        CALL umPrint(umMessage)
        WRITE(umMessage,'(A,E16.8)') 'Value from model dump : ',         &
            a_levdepc(jetarho+j-1)
        CALL umPrint(umMessage)
        ErrorStatus=10
        WRITE(cmessage,'(A)')                                            &
            'INBOUNDA : Level dependent constant mismatch'

        CALL Ereport(RoutineName,ErrorStatus,cmessage)

      END IF
    END DO ! j

  END IF ! IF (fixhd_bounda(110,1) > 0)

  ! Check that for variable resolution input, the LBC file has the same
  ! grid as the input dump
  IF ( (a_fixhd(115) > 0 .AND. a_fixhd(120) > 0) .OR.              &
      (fixhd_bounda(115,1) > 0 .AND. fixhd_bounda(120,1) > 0) ) THEN
    ! Input dump is variable resolution
    IF (fixhd_bounda(116,1) /= a_fixhd(116) .OR.                  &
        fixhd_bounda(117,1) /= a_fixhd(117) .OR.                  &
        fixhd_bounda(121,1) /= a_fixhd(121) .OR.                  &
        fixhd_bounda(122,1) /= a_fixhd(122) ) THEN

      WRITE(cmessage,'(A)') "Variable resolution grid " //        &
        "dimensions in LBC file do not match input dump"
      errorstatus = -10
      CALL ereport(RoutineName, errorstatus, cmessage)
    END IF

    DO j = 1, a_fixhd(117)
      DO i = 1, a_fixhd(116)
        k = (j-1) * a_fixhd(116) + i
        IF (near_equal_real(a_rowdepc(k), a_rowdepc_bo(i,j))) THEN
          WRITE(cmessage,'(A)') "Row dependent constants " //     &
        "in LBC file do not match input dump"
          errorstatus = -11
          CALL ereport(RoutineName, errorstatus, cmessage)
        END IF
      END DO
    END DO

    DO j = 1, a_fixhd(122)
      DO i = 1, a_fixhd(121)
        k = (j-1) * a_fixhd(121) + i
        IF (near_equal_real(a_coldepc(k), a_coldepc_bo(i,j))) THEN
          WRITE(cmessage,'(A)') "Column dependent constants " //  &
        "in LBC file do not match input dump"
          errorstatus = -12
          CALL ereport(RoutineName, errorstatus, cmessage)
        END IF
      END DO
    END DO
  END IF ! IF variable resolution input dump




  !       Set update interval
  !       IF update interval includes months or years, a 360 day
  !       calender assumed.

          ! Save previous value of rim_stepsa:
  IF (.NOT. L_FirstCall) rim_stepsa_OLD = rim_stepsa

  rim_stepsa=((fixhd_bounda(35,1)*8640+fixhd_bounda(36,1)*720     &
   +fixhd_bounda(37,1)*24+fixhd_bounda(38,1))*3600                &
   +fixhd_bounda(39,1)*60+fixhd_bounda(40,1))                     &
   *steps_per_periodim(atmos_im)/secs_per_periodim(atmos_im)

  ! Check that rim_stepsa has not changed:
  IF (.NOT. L_FirstCall) THEN
    IF (rim_stepsa /= rim_stepsa_OLD) THEN
      ErrorStatus = 1
      WRITE (cmessage,'(A,I8,A,I8,A)')                            &
      'Boundary updating period (rim_stepsa) has changed from ',  &
      rim_stepsa_OLD, ' to ', rim_stepsa, '. Not allowed!'

      CALL EReport (RoutineName, ErrorStatus, cmessage)
    END IF
  END IF

  ! Initialise Current_LBC_Step

  IF (L_FirstCall) THEN
    Current_LBC_Step = 1 + stepim(atmos_im)
    IF (PrintStatus >= PrStatus_Normal) THEN
      WRITE(umMessage,'(A,I8,A,I8)') ' INBOUNDA : Timestep ',            &
            stepim(atmos_im), ' Current_LBC_Step ',Current_LBC_Step
      CALL umPrint(umMessage)
    END IF
  END IF

  IF (Num_ALBCs == 2 .AND. L_FirstCall) THEN

    ! Calculate step on which to swap boundary files:
    ALBC_SwapStep = ALBC2_StartTime_steps - rim_stepsa

    ! Check that ALBC_SwapStep is non-negative:
    IF (ALBC_SwapStep < 0) THEN
      ErrorStatus = 1
      WRITE (cmessage,'(A,I8,A)')                                 &
        'Step on which to swap boundary files (ALBC_SwapStep = ', &
        ALBC_SwapStep, ') must be non-negative'

      CALL EReport (RoutineName, ErrorStatus, cmessage)
    END IF

    ! Check that ALBC_SwapStep is a multiple of rim_stepsa:
    IF (MOD(ALBC_SwapStep, rim_stepsa) /= 0) THEN
      ErrorStatus = 1
      WRITE (cmessage,'(A,I8,A,I8)')                              &
        'Step on which to swap boundary files (ALBC_SwapStep = ', &
        ALBC_SwapStep, ') must be a a multiple of the boundary '//&
        'updating period (rim_stepsa = ', rim_stepsa, ')'
      ErrorStatus=101

      CALL EReport (RoutineName, ErrorStatus, cmessage)
    END IF

  END IF

END IF

!   2   Set interval for setting any boundary field

boundary_stepsim(atmos_im) = rim_stepsa

!   3    Check LOOKUP Table

j1=0
IF (.NOT. (model_type == mt_bi_cyclic_lam .OR.                  &
    model_type == mt_cyclic_lam    .OR.                         &
    l_force_lbc .OR. l_fixed_lbcs) ) THEN
  j1=fixhd_bounda(152,1)

  IF (fixhd_bounda(150,1) >  0) THEN

    ! Set up list of variables expected to be boundary updated.

    CALL assign_item_bounda(item_bounda)

    ! DEPENDS ON: chk_look_bounda
    CALL chk_look_bounda(                                         &
      item_bounda,full_lookup_bounda,                             &
                         ErrorStatus,cmessage)

    IF (ErrorStatus  /=  0) THEN
      ! Use message returned by chk_look_bounda in CALL to ereport

      CALL ereport(RoutineName,ErrorStatus,cmessage)
    END IF

    !  Find start position in lookup tables

              ! Get days/seconds since basis time:
    ! DEPENDS ON: time2sec
    CALL time2sec (i_year,             i_month,                   &
                   i_day,              i_hour,                    &
                   i_minute,           i_second,                  &
                   basis_time_days,    basis_time_secs,           &
                   elapsed_days,       elapsed_secs,              &
                   lcal360)

    ! Get current model time in days/seconds-in-day:
    current_time_days = basis_time_days + elapsed_days +          &
                     (basis_time_secs + elapsed_secs)/isec_per_day

    current_time_secs = MOD(basis_time_secs+elapsed_secs,         &
                                            isec_per_day)

    ! Get days/seconds to start of boundary data:
    ! DEPENDS ON: time2sec
    CALL time2sec (fixhd_bounda(21,1), fixhd_bounda(22,1),        &
                   fixhd_bounda(23,1), fixhd_bounda(24,1),        &
                   fixhd_bounda(25,1), fixhd_bounda(26,1),        &
                   current_time_days,  current_time_secs,         &
                   days_to_data_start, secs_to_data_start,        &
                   lcal360)

    ! Get days/seconds to end of boundary data:
    ! DEPENDS ON: time2sec
    CALL time2sec (fixhd_bounda(28,1), fixhd_bounda(29,1),        &
                   fixhd_bounda(30,1), fixhd_bounda(31,1),        &
                   fixhd_bounda(32,1), fixhd_bounda(33,1),        &
                   current_time_days,  current_time_secs,         &
                   days_to_data_end,   secs_to_data_end,          &
                   lcal360)

    ! Get steps to start of boundary data:
    ! DEPENDS ON: tim2step
    CALL tim2step (days_to_data_start,                            &
                   secs_to_data_start,                            &
                   steps_per_periodim(atmos_im),                  &
                   secs_per_periodim(atmos_im),                   &
                   steps_to_data_start)

    ! Get steps to end of boundary data:
    ! DEPENDS ON: tim2step
    CALL tim2step (days_to_data_end,                              &
                   secs_to_data_end,                              &
                   steps_per_periodim(atmos_im),                  &
                   secs_per_periodim(atmos_im),                   &
                   steps_to_data_end)

    ! Get steps from basis time to start of boundary data:
    basis_to_data_start_steps = stepim(atmos_im) + steps_to_data_start

    ! Check that the above is a multiple of rim_stepsa only
    ! when two boundary files exist
    IF ((Num_ALBCs == 2) .AND.                                    &
        (MOD(basis_to_data_start_steps, rim_stepsa) /= 0)) THEN
      ErrorStatus = 1
      WRITE (cmessage,'(A,I8,A,I8,A)')                            &
        'Steps from basis time to start of boundary data (',      &
        basis_to_data_start_steps, ') must be a multiple of '//   &
        'the boundary updating period (', rim_stepsa, ')'

      CALL EReport (RoutineName, ErrorStatus, cmessage)
    END IF

    ! There are two situations we need to cater for here:
    !
    !   1. There is no boundary data in memory.
    !   2. The boundary data valid at the start of the current
    !      boundary data interval is already in memory, but
    !      data for the end of the interval needs to be read in.
    !
    ! The first case applies IF and only IF we are on the first
    ! CALL to this routine.
    !
    ! In the first case, we need to make sure that there is
    ! boundary data valid at the start of the current boundary
    ! data interval, and THEN point to it. We also need to make
    ! sure that data is available for the end of the interval. The
    ! only exception to the latter is the case where the data for
    ! the end of the interval is to come from a different boundary
    ! file, which may be the case IF two boundary files are being
    ! used.
    !
    ! In the second case, we need to make sure that there is
    ! boundary data valid at the end of the current boundary data
    ! interval, and THEN point to it.

    ! Steps to start of current boundary data interval:
    steps_to_bdi_start = MOD(basis_to_data_start_steps,           &
                             rim_stepsa)

    ! Steps to end   of current boundary data interval:
    steps_to_bdi_end   = steps_to_bdi_start + rim_stepsa

    ! Initialise bndary_offset for use in model_time_mod
    bndary_offsetim(atmos_im) = MOD(-basis_to_data_start_steps,   &
                                 rim_stepsa)

    IF (L_FirstCall) THEN ! Case 1

      IF (steps_to_data_start <= steps_to_bdi_start .AND.         &
          steps_to_data_end   >= steps_to_bdi_start) THEN
        nbound_lookup(1) = (-steps_to_data_start / rim_stepsa)    &
                         * (rim_lookupsa-1) + 2
      ELSE IF (steps_to_data_start > steps_to_bdi_start) THEN
        ErrorStatus = 101
        WRITE (cmessage,'(A)')                                    &
          'Boundary data starts after start of current '//        &
          'boundary data interval'

        CALL EReport (RoutineName, ErrorStatus, cmessage)
      ELSE
        ErrorStatus = 102
        WRITE (cmessage,'(A)')                                    &
          'Boundary data ends before start of current '//         &
          'boundary data interval'

        CALL EReport (RoutineName, ErrorStatus, cmessage)
      END IF

      ! Is the data for the end of the current boundary data
      ! interval to come from the same boundary file?
      this_albc_for_bdi_end =.TRUE.
      IF (NUM_ALBCs == 2) THEN
        IF (stepim(atmos_im) <= ALBC_SwapStep) THEN
          this_albc_for_bdi_end =.FALSE.
        END IF
      END IF

      ! IF so, check that the data exists:
      IF (this_albc_for_bdi_end) THEN
        IF (steps_to_data_start > steps_to_bdi_end) THEN
          ErrorStatus = 101
          WRITE (cmessage,'(A)')                                  &
            'Boundary data starts after end of current '//        &
            'boundary data interval'

          CALL EReport (RoutineName, ErrorStatus, cmessage)
        ELSE IF (steps_to_data_end < steps_to_bdi_end) THEN
          ErrorStatus = 102
          WRITE (cmessage,'(A)')                                  &
            'Boundary data ends before end of current '//         &
            'boundary data interval'

          CALL EReport (RoutineName, ErrorStatus, cmessage)
        END IF
      END IF

    ELSE ! Not first CALL (Case 2)

      IF (steps_to_data_start <= steps_to_bdi_end .AND.           &
          steps_to_data_end   >= steps_to_bdi_end) THEN
        nbound_lookup(1) = (1 - steps_to_data_start / rim_stepsa) &
                         * (rim_lookupsa-1) + 2
      ELSE IF (steps_to_data_start > steps_to_bdi_end) THEN
        ErrorStatus = 101
        WRITE (cmessage,'(A)')                                    &
          'Boundary data starts after end of current '//          &
          'boundary data interval'

        CALL EReport (RoutineName, ErrorStatus, cmessage)
      ELSE
        ErrorStatus = 102
        WRITE (cmessage,'(A)')                                    &
          'Boundary data ends before end of current '//           &
          'boundary data interval'

        CALL EReport (RoutineName, ErrorStatus, cmessage)
      END IF

    END IF ! (L_FirstCall)

  END IF

END IF ! lateral boundary

L_FirstCALL = .FALSE.

!   4   End of routine

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE inbounda

