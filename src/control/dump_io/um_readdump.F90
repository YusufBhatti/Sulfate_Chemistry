! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    SUBROUTINE UM_READDUMP---------------------------------------
!
!    Purpose: Reads in model dump on unit NFTIN and checks model
!             and dump dimensions for consistency.
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Dump I/O

! Subroutine Interface

MODULE um_readdump_mod

IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UM_READDUMP_MOD'

CONTAINS

SUBROUTINE um_readdump(nftin,fixhd,len_fixhd                      &
 ,inthd,len_inthd                                                 &
 ,realhd,len_realhd                                               &
 ,levdepc,len1_levdepc,len2_levdepc                               &
 ,rowdepc,len1_rowdepc,len2_rowdepc                               &
 ,coldepc,len1_coldepc,len2_coldepc                               &
 ,flddepc,len1_flddepc,len2_flddepc                               &
 ,extcnst,len_extcnst                                             &
 ,dumphist,len_dumphist                                           &
 ,cfi1,len_cfi1                                                   &
 ,cfi2,len_cfi2                                                   &
 ,cfi3,len_cfi3                                                   &
 ,lookup,len1_lookup,len2_lookup                                  &
 ,mpp_lookup,mpp_len1_lookup                                      &
 ,submodel_id,n_objs_d1,d1_addr                                   &
 ,len_data,d1                                                     &
 ,read_header                                                     &
  )

USE stextend_mod, ONLY: in_s
USE application_description, ONLY: isSmallExec
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE io
USE ios_common, ONLY: IOS_use_helpers
USE ereport_mod, ONLY: ereport
USE umPrintMgr
USE UM_ParVars
USE decomp_params, ONLY: decomp_standard_atmos
USE Decomp_DB
USE lookup_addresses
USE rimtypes, ONLY: rima_type_orog, rima_type_norm
USE lbc_mod, ONLY: lenrima
USE submodel_mod, ONLY: atmos_sm

USE d1_array_mod, ONLY: d1_list_len, d1_object_type, d1_section, d1_item,      &
                        d1_no_levels, diagnostic, d1_grid_type, d1_halo_type,  &
                        d1_length, d1_gridpoint_code, d1_north_code,           &
                        d1_east_code, d1_west_code, d1_south_code,             &
                        d1_proc_no_code, d1_imodl

USE cppxref_mod, ONLY: ppx_pt_code, ppx_pl_code, ppx_atm_tzonal, ppx_atm_ozone

USE pr_look_mod, ONLY: pr_look
! JULES
USE land_tile_ids,  ONLY: surface_type_ids, ml_snow_type_ids
USE max_dimensions, ONLY: ntype_max, snow_layers_max
USE jules_surface_mod, ONLY: l_aggregate
USE ppxlook_mod, ONLY: exppxi
!$ USE omp_lib
USE errormessagelength_mod, ONLY: errormessagelength
USE missing_data_mod, ONLY: imdi
USE model_domain_mod, ONLY: FH_GridStagger_Endgame

USE read_serial_mod, ONLY: read_serial

USE readflds_mod, ONLY: readflds

USE stparam_mod, ONLY: st_accum_code, st_time_mean_code, st_time_series_mean,  &
                       st_min_code, st_max_code, vert_mean_base,               &
                       vert_mean_top, global_mean_base, global_mean_top,       &
                       stash_nmdi_mask_code, zonal_mean_base, st_input_code,   &
                       st_output_addr

USE stash_array_mod, ONLY:                                                     &
    stlist, stindex

USE um_types, ONLY: integer64

USE fort2c_portio_interfaces, ONLY: portiodetachallhelpers, portioaddhelper

IMPLICIT NONE

INTEGER, INTENT(IN) ::                                            &
  nftin                                                           &
                 !IN Unit number of dump
, len_fixhd                                                       &
                 !IN Length of fixed length header
, len_inthd                                                       &
                 !IN Length of integer header
, len_realhd                                                      &
                 !IN Length of real header
, len1_levdepc                                                    &
                 !IN 1st dim of level dep consts
, len2_levdepc                                                    &
                 !IN 2nd dim of level dep consts
, len1_rowdepc                                                    &
                 !IN 1st dim of row dep consts
, len2_rowdepc                                                    &
                 !IN 2nd dim of row dep consts
, len1_coldepc                                                    &
                 !IN 1st dim of column dep consts
, len2_coldepc                                                    &
                 !IN 2nd dim of column dep consts
, len1_flddepc                                                    &
                 !IN 1st dim of field dep consts
, len2_flddepc                                                    &
                 !IN 2nd dim of field dep consts
, len_extcnst                                                     &
                 !IN Length of extra constants
, len_dumphist                                                    &
                 !IN Length of history block
, len_cfi1                                                        &
                 !IN Length of comp field index 1
, len_cfi2                                                        &
                 !IN Length of comp field index 2
, len_cfi3                                                        &
                 !IN Length of comp field index 3
, len1_lookup                                                     &
                 !IN 1st dim of lookup
, len2_lookup                                                     &
                 !IN 2nd dim of lookup
, mpp_len1_lookup                                                 &
                 !IN 1st dim of MPP lookup
, submodel_id                                                     &
                 !IN submodel of dump
, n_objs_d1                                                       &
                 !IN number of objects (3D fields) in D1
, len_data        !IN length of model data

INTEGER ::                                                        &
  fixhd(len_fixhd)                                                &
                     !IN Fixed length header
, inthd(len_inthd)                                                &
                     !IN Integer header
, lookup(len1_lookup,len2_lookup)                                 &
                     !IN PP lookup tables
, cfi1(len_cfi1+1)                                                &
                     !IN Compressed field index no 1
, cfi2(len_cfi2+1)                                                &
                     !IN Compressed field index no 2
, cfi3(len_cfi3+1)                                                &
                     !IN Compressed field index no 3

, mpp_lookup(mpp_len1_lookup,len2_lookup)
                     !OUT Local processor lookup

REAL ::                                                           &
  realhd(len_realhd)                                              &
                     !IN Real header
, levdepc(1+len1_levdepc*len2_levdepc)                            &
                                       !IN Lev dep consts
, rowdepc(1+len1_rowdepc*len2_rowdepc)                            &
                                       !IN Row dep consts
, coldepc(1+len1_coldepc*len2_coldepc)                            &
                                       !IN Col dep consts
, flddepc(1+len1_flddepc*len2_flddepc)                            &
                                       !IN Field dep consts
, extcnst(len_extcnst+1)                                          &
                           !IN Extra constants
, dumphist(len_dumphist+1)                                        &
                           !IN History block

, d1(len_data)       !OUT Local subdomain of dump

LOGICAL ::                                                        &
 read_header         !IN  True if header is to be read in

INTEGER ::                                                        &
  d1_addr(d1_list_len,n_objs_d1)
                     ! IN D1 addressing info.

! Local variables

INTEGER ::                                                        &
  start_block                                                     &
                ! first word of field data in dump (unused)
, object_index                                                    &
                ! pointer to entry in D1_ADDR
, deriv_object_index                                              &
                ! pointer to entry in D1_ADDR
, input_obj                                                       &
                ! stash object index for deriv_object_index
, level                                                           &
                ! level number of multi-level field
, d1_item_code                                                    &
                ! sec/item in d1_addr converted into single code
, number_of_fields                                                &
                ! total number of fields to read in
, field_start                                                     &
                ! start address of a field in the file
, data_size                                                       &
                ! number of words of data on disk for a field
, data_read_size                                                  &
                ! total number of words to read for a field
, data_full_size                                                  &
                ! total number of words after any unpacking
, len_io                                                          &
                ! number of words of data successfully read
, k, j, i                                                         &
                ! loop counter over fields
, orig_decomp                                                     &
                ! current decomposition on entry
, address                                                         &
                ! start address of field in D1
, d1_start_read                                                   &
                ! Same as start address, but may be modified in
                ! case of zero-length fields at the end of the
                ! d1-array
, local_len                                                       &
                ! number of words of data put into D1 on this
                ! processor
, isfc                                                            &
                ! Loop counter for surface types
, ipseudo                                                         &
                ! pseudo type code
, ipseudl                                                         &
                ! pseudo last level code
, im_ident                                                        &
                ! model number
, isec                                                            &
                ! section number
, iitm, iitm_pre, iitm_next                                       &
                ! item number
, ps_lev_count                                                    &
                ! counter
, n_plevels                                                       &
                ! Number of psuedo levels based on ipseudl
, expected_lbplev(ntype_max*snow_layers_max)                      &
! Expected header given JULES surface type config
, idummy                                                          &
, d1_proc_no                                                      &
, check_proccode

! Error reporting
INTEGER ::    icode       ! =0 normal exit; >0 error exit
CHARACTER(LEN=errormessagelength) :: cmessage    ! Error message
CHARACTER(LEN=*), PARAMETER :: RoutineName='UM_READDUMP'
CHARACTER(LEN=20) :: umFormat

!$ INTEGER :: thread
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER :: first_field

INTEGER, ALLOCATABLE :: obj_type_map(:), north(:), south(:),&
                        east(:), west(:), gridpoint(:), proccode(:)

INTEGER :: expand

INTEGER, PARAMETER :: unset = -1

INTEGER, PARAMETER :: temporal_proccode = 128 + 4096 + 8192

LOGICAL :: repeat_stash ! logical to test if a field in the dump has the same
                        ! STASH code as the previous field. Relevant for
                        ! partially processed fields.

LOGICAL :: partial_processed_from_diagnostic ! logical to test if a field in
                                             ! the dump is partially processed
                                             ! from a diganostic field.

INTEGER(KIND=integer64) :: checksum_shft
INTEGER :: checksum_int, check_len
INTEGER(KIND=integer64), ALLOCATABLE :: d1_field(:)
CHARACTER(LEN=128) :: checksum_str

!--------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
icode=0
cmessage=''

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(thread)
!$ thread=omp_get_thread_num()

!$ IF (thread==0) THEN

IF (mype  ==  0) THEN
  CALL umPrint('',src='um_readdump')
  WRITE(umMessage,'('' READING UNIFIED MODEL DUMP ON UNIT'',I3)')nftin
  CALL umPrint(umMessage,src='um_readdump')
  CALL umPrint(' #####################################',src='um_readdump')
  CALL umPrint('',src='um_readdump')
END IF

! Change to the relevant decomposition type for this dump

orig_decomp=current_decomp_type

IF (submodel_id  ==  atmos_sm) THEN
  IF (current_decomp_type  /=  decomp_standard_atmos)                          &
    CALL change_decomposition(decomp_standard_atmos,icode)

ELSE  ! unsupported decomposition type
  WRITE(umMessage,'(A,A,I0)')                                                  &
    'UM_READDUMP : Could not change to decomposition required ',               &
    'for submodel type ',submodel_id
  CALL umPrint(umMessage,src='um_readdump')
  icode=1
  cmessage='Unsupported submodel for MPP code'

  CALL ereport ( routinename, icode, cmessage )

END IF

IF (icode  /=  0) THEN
  icode=2
  CALL umPrint('UM_READDUMP : Error - Could not set decomposition '//          &
      'for selected submodel.',src='um_readdump')
  cmessage='Unsupported decomposition selected for MPP code'

  CALL ereport(routinename,icode,cmessage)
END IF

! Read in the header records, and do a consistency check

IF (read_header) THEN

  ! DEPENDS ON: readhead
  CALL readhead(nftin,fixhd,len_fixhd,                            &
                inthd,len_inthd,                                  &
                realhd,len_realhd,                                &
                levdepc,len1_levdepc,len2_levdepc,                &
                rowdepc,len1_rowdepc,len2_rowdepc,                &
                coldepc,len1_coldepc,len2_coldepc,                &
                flddepc,len1_flddepc,len2_flddepc,                &
                extcnst,len_extcnst,                              &
                dumphist,len_dumphist,                            &
                cfi1,len_cfi1,                                    &
                cfi2,len_cfi2,                                    &
                cfi3,len_cfi3,                                    &
                lookup,len1_lookup,len2_lookup,                   &
                len_data,                                         &
                start_block,icode,cmessage)

  IF (icode  /=  0) THEN
    WRITE(umMessage,'(A,A,I0)') 'UM_READDUMP : Error reading dump header ',    &
        'on unit ',nftin
    CALL umPrint(umMessage,src='um_readdump')
    WRITE(umMessage,'(A,I0,A,A)') 'Return code from READHEAD was ',icode,      &
        ' and error message was ',cmessage
    CALL umPrint(umMessage,src='um_readdump')
    icode=3
    cmessage='Error reading dump header'

    CALL ereport ( routinename, icode, cmessage )
  END IF

END IF  ! IF (READ_HEADER)


IF (fixhd(9) /= FH_GridStagger_Endgame) THEN
  WRITE(umMessage,'(A)')                                                       &
      'UM_READDUMP : The UM only supports reading of ENDGame dumps.'
  CALL umPrint(umMessage,src='um_readdump')
  icode=4
  cmessage='Dump not identified as an ENDGame dump. fixhd(9) /= 6'
END IF


IF (fixhd(160)  >   0) THEN ! If there is data to read

  ! Loop over fields and read into D1

  number_of_fields=fixhd(152)

  ALLOCATE(obj_type_map(number_of_fields))
  ALLOCATE(north(number_of_fields))
  ALLOCATE(south(number_of_fields))
  ALLOCATE(east(number_of_fields))
  ALLOCATE(west(number_of_fields))
  ALLOCATE(gridpoint(number_of_fields))
  ALLOCATE(proccode(number_of_fields))

  address=1
  object_index=1
  level=1

  iitm = imdi
  ps_lev_count = 0

  DO k = 1,number_of_fields  ! loop over fields to read in

    obj_type_map(k) = d1_addr(d1_object_type,object_index)
    north(k) = d1_addr(d1_north_code,object_index)
    south(k) = d1_addr(d1_south_code,object_index)
    east(k) = d1_addr(d1_east_code,object_index)
    west(k) = d1_addr(d1_west_code,object_index)
    gridpoint(k) = d1_addr(d1_gridpoint_code,object_index)
    d1_proc_no = d1_addr(d1_proc_no_code,object_index)
    proccode(k) = d1_proc_no
    check_proccode = 0

    ! If zonal mean, add 2**6=64 to proccode
    IF (gridpoint(k) >= zonal_mean_base) THEN
      check_proccode = check_proccode + 64
    END IF

    ! If accumulation or time mean, add 2**7=128
    IF ( (d1_proc_no == st_accum_code) .OR. (d1_proc_no == st_time_mean_code)  &
         .OR. (d1_proc_no == st_time_series_mean) ) THEN
      check_proccode = check_proccode + 128
    END IF

    ! If minimum value, add 2**12=4096
    IF (d1_proc_no == st_min_code) THEN
      check_proccode = check_proccode + 4096
    END IF

    ! If maximum value, add 2**13=8192
    IF (d1_proc_no == st_max_code) THEN
      check_proccode = check_proccode + 8192
    END IF

    ! If mean over vertical layer(s), add 2**11=2048
    IF ( (gridpoint(k)>=vert_mean_base .AND. gridpoint(k)<vert_mean_top) .OR. &
         (gridpoint(k)>=global_mean_base .AND. gridpoint(k)<global_mean_top) )&
    THEN
      check_proccode = check_proccode + 2048
    END IF

    ! If minimal mdi masking has been used, add 2**18=262144
    IF (gridpoint(k) == stash_nmdi_mask_code) THEN
      check_proccode = check_proccode + 262144
    END IF

    IF (lookup(lblrec,k)  >   0) THEN ! If there's data in
                                      ! the field

      ! Check that DATA_TYPE is valid no: +/-1 to +/-3
      IF (( ABS(lookup(data_type,k))  >=  1) .AND.                             &
          ( ABS(lookup(data_type,k))  <=  3)) THEN

        ! Check that the diagnostic in the dump matches that expected
        ! from D1_ADDR

        IF (obj_type_map(k)  ==  diagnostic) THEN

          d1_item_code= (d1_addr(d1_section,object_index)*1000) +              &
                         d1_addr(d1_item,object_index)

          IF ( lookup(item_code,k) /= d1_item_code ) THEN

            WRITE(umMessage,'(A,I0,A,I0,A,I0)')                                &
                'UM_READDUMP : Dump field ',k,                                 &
                ' does not match STASH request for section ',                  &
                d1_addr(d1_section,object_index),                              &
                ' item ', d1_addr(d1_item,object_index)

            CALL umPrint(umMessage,src='um_readdump')
            WRITE(umMessage,'(A,I0)') 'Expected item code ',d1_item_code
            CALL umPrint(umMessage,src='um_readdump')
            WRITE(umMessage,'(A,I0)') 'Found item code ',lookup(item_code,k)
            CALL umPrint(umMessage,src='um_readdump')
            cmessage='UM_READDUMP Dump does not match STASH list'
            icode=5

            CALL ereport(routinename,icode,cmessage)

          END IF ! IF (LOOKUP(ITEM_CODE,K)  /=  d1_item_code)

          ! check the prcessing code we have from the dump file matches what we
          ! expect from our calculated processing code.

          IF ( lookup(lbproc,k) /= check_proccode ) THEN

            ! Dumps may contain a series of one or more partially
            ! temporally-processed diagnostic fields.
            !
            ! Where these diagnostics are calculated from a prognostic field,
            ! they immediatly follow the primary (unprocessed prognostic)
            ! version of the field.
            !
            ! If however they are calculated from a diagnostic field, there
            ! is no primary version in the dump.
            !
            ! Unfortunately, these partially processed fields do not have
            ! complete lookup headers (most elements are set to imdi/rmdi).
            !
            ! We therefore need to check that if the processing codes didn't
            ! match, this is because either:
            !  i)  this is a partially processed diagnostic, which is derived
            !      from a prognostic, and so has the same stash code as the
            !      previous field.
            !  ii) this is a partially processed diagnostic, but it is derived
            !      from a diagnostic field, and so there is no primary field.

            repeat_stash = .FALSE.

            IF ( lookup(lbproc,k) == imdi .AND. k>1 ) THEN
              IF (lookup(item_code,k)==lookup(item_code,k-1)) THEN
                IF (IAND(check_proccode, temporal_proccode) > 0) THEN
                  ! this field is a repeat of the same STASH code with partial
                  ! processing applied
                  repeat_stash = .TRUE.
                END IF
              END IF
            END IF

            IF (.NOT.repeat_stash) THEN
              ! This is not a repeat of the same STASH code with partial
              ! processing applied, so was is derived from a diagnostic?

              partial_processed_from_diagnostic = .FALSE.

              ! look up the field this diagnostic is dervived from
              input_obj = stlist(st_input_code,                                &
                                 stindex(1,                                    &
                                         d1_addr(d1_item,object_index),        &
                                         d1_addr(d1_section,object_index),     &
                                         d1_addr(d1_imodl,object_index)))

              ! Check if that field was a diagnostic
              IF ( input_obj < 0 ) THEN

                deriv_object_index = stlist(st_output_addr,-input_obj)

                IF ( d1_addr(d1_object_type,deriv_object_index)                &
                     ==  diagnostic ) THEN
                  IF (IAND(check_proccode, temporal_proccode) > 0) THEN
                    ! this field is a field with partial processing applied,
                    ! and dervived from a diagnostic field.
                    partial_processed_from_diagnostic = .TRUE.
                  END IF
                END IF

              ELSE IF ( input_obj == 1 ) THEN

                IF (IAND(check_proccode, temporal_proccode) > 0) THEN
                  ! this field is a field with partial processing applied,
                  ! and dervived from a diagnostic in STASHwork space.
                  partial_processed_from_diagnostic = .TRUE.
                END IF

              END IF

              ! If this is neither a repeated STASH code, nor a partially
              ! processed diagnostic derived from another diagnostic, then fail.
              IF (.NOT.partial_processed_from_diagnostic) THEN
                WRITE(umMessage,'(A,I0,A,I0,A,I0)')                            &
                      'UM_READDUMP : Dump field ',k,                           &
                      ' does not match STASH processing code for section ',    &
                      d1_addr(d1_section,object_index),                        &
                      ' item ', d1_addr(d1_item,object_index)
                CALL umPrint(umMessage,src='um_readdump')
                WRITE(umMessage,'(A,I0)') 'Expected processing code ',         &
                                          check_proccode
                CALL umPrint(umMessage,src='um_readdump')
                WRITE(umMessage,'(A,I0)') 'Found processing code ',            &
                                          lookup(lbproc,k)
                CALL umPrint(umMessage,src='um_readdump')

                cmessage='UM_READDUMP Dump does not match STASH list:' //      &
                         ' Wrong processing code.'
                icode=6

                CALL ereport(routinename,icode,cmessage)
              END IF

            END IF

          END IF ! lookup(lbproc,k) /= check_proccode

        END IF ! IF (D1_ADDR(d1_object_type,object_index)  ==  diagnostic)

        ! For atmosphere zonal ozone fields - set to zonal grid type
        IF ( (d1_addr(d1_grid_type,object_index)  ==  ppx_atm_ozone)           &
            .AND.  (lookup(lbnpt,k)  ==  1) ) THEN

          d1_addr(d1_grid_type,object_index) = ppx_atm_tzonal
          WRITE(umMessage,'(A,I0,A)')  'UM_READDUMP : Dump field ', k,         &
                              ' switching from ppx_atm_ozone to ppx_atm_tzonal '
          CALL umPrint(umMessage,src='um_readdump')

        END IF

      ELSE ! Error in LOOKUP(DATA_TYPE,K)

        IF (( fixhd(5)  <   6) .OR.                                            &
            ( fixhd(5)  >   8)) THEN ! Not AC, Var or Cx

          CALL pr_look(lookup,k)

        END IF

        WRITE(umMessage,'(A,I0,A,I0)') 'um_readdump : failure for field ',k,   &
            ' of ',number_of_fields
        CALL umPrint(umMessage,src='um_readdump')
        WRITE(umMessage,'(A,I0)') 'LOOKUP(DATA_TYPE,K)= ',lookup(data_type,k)
        CALL umPrint(umMessage,src='um_readdump')
        icode=7
        cmessage='Invalid data type ( LOOKUP(DATA_TYPE,K) )'

        CALL ereport(routinename, icode, cmessage)

      END IF

    END IF ! If there was data in the field

#if !defined(UTILIO)

    ! JULES flexible tiles facility
    ! Check surface tile types in the input dump tiles prognostics against
    !   the tile types specified in the UI job set up.
    ! We require that every surface tiles variable in the dump matches the
    !   tile types specified in the job set up; if any does not, we abort.
    ! Note that pseudl=9 denotes a tiles variable defined on all tile types,
    !   pseudl=8 denotes a tiles variable defined on all vegetation tile types
    ! Cannot do this check in the aggregated tile case

    IF ( .NOT. l_aggregate ) THEN
      im_ident = lookup(model_code,k)
      isec     = 0
      iitm_pre = iitm
      iitm     = lookup(item_code,k)

      IF (iitm < 1000) THEN  ! Section 0 item

        ipseudo = exppxi(im_ident, isec, iitm, ppx_pt_code, icode, cmessage)
        ipseudl = exppxi(im_ident, isec, iitm, ppx_pl_code, icode, cmessage)

        IF ( icode == 0 ) THEN

          IF ( ipseudo == 9 ) THEN      ! Tiled variable

            ! DEPENDS ON: pslevcod
            CALL pslevcod(ipseudl,n_plevels,'L',idummy,cmessage)
            IF ( ipseudl == 11 ) THEN
              expected_lbplev(1:n_plevels) = ml_snow_type_ids(1:n_plevels)
            ELSE
              expected_lbplev(1:n_plevels) = surface_type_ids(1:n_plevels)
            END IF

            ps_lev_count = ps_lev_count + 1

            ! Check to make sure STASH number has changed when expected
            IF ( ps_lev_count == 1 .AND. iitm_pre == iitm ) THEN
              icode = 8
              WRITE(cmessage,'(A, I5, A)')                                     &
                    'Supposed to be 1st level of new STASH item', iitm,        &
                    ' however STASH number has not changed from previous'
              CALL ereport(routinename, icode, cmessage)
            ELSE IF ( ( ps_lev_count > 1 .AND. ps_lev_count <= n_plevels )     &
                     .AND. iitm_pre /= iitm) THEN
              icode = 8
              WRITE(cmessage,'(A, I5)')                                        &
                    'STASH number has changed in the midst of STASH item',     &
                    iitm
              CALL ereport(routinename, icode, cmessage)
            END IF

            ! Check surface type matches for this tiled variable
            ! Surface types & levels should be in the right order otherwise
            ! recon needs to be run
            DO isfc = 1, n_plevels
              IF ( expected_lbplev(isfc) == lookup(lbplev,k) ) THEN
                IF ( isfc /= ps_lev_count ) THEN
                  WRITE(umMessage,'(A, I5, A)')                                &
                        'UM_READDUMP : surface stash item ',                   &
                        iitm, ' does not match surface type configuration'
                  CALL umPrint(umMessage,src='um_readdump')
                  WRITE(umFormat,'(A,I3,A)') '(A,',n_plevels,'(I7))'
                  WRITE(umMessage,umFormat) 'Dump:    ',                       &
                        lookup(lbplev,k-ps_lev_count+1:k-ps_lev_count+n_plevels)
                  CALL umPrint(umMessage,src='um_readdump')
                  WRITE(umMessage,umFormat) 'Expected:',                       &
                        expected_lbplev(1:n_plevels)
                  CALL umPrint(umMessage,src='um_readdump')
                  icode = 8
                  cmessage = 'Dump does not match surface type ' //            &
                             'configuration; please check then run recon'
                  CALL ereport(routinename, icode, cmessage)
                END IF
                CYCLE
              END IF
            END DO

            ! Reached the end of the pt_code 9 field so reset counter
            IF ( ps_lev_count == n_plevels ) THEN
              IF ( k /= number_of_fields ) THEN
                ! Double check the stash number of the next field is different
                iitm_next = lookup(item_code,k+1)
                IF ( iitm_next /= iitm ) THEN
                  ps_lev_count = 0
                ELSE
                  icode = 8
                  WRITE(cmessage,'(A, I5, A, I3)')                             &
                        'Number of levels for STASH item', iitm,               &
                        ' is not what was expected:', n_plevels
                  CALL ereport(routinename, icode, cmessage)
                END IF
              END IF
            END IF
          END IF ! ipseudo == 9

        END IF ! icode == 0

      END IF  ! JULES iitm < 1000

    END IF ! .NOT. l_aggregate

#endif

    level=level+1
    IF (level  >   d1_addr(d1_no_levels,object_index)) THEN
      level=1
      object_index=object_index+1
    END IF

  END DO ! K : loop over fields to read in

  expand = 1
  first_field = 1

  CALL readflds(nftin,                                                         &
                number_of_fields,                                              &
                first_field,                                                   &
                lookup,                                                        &
                d1,                                                            &
                fixhd,                                                         &
                expand,                                                        &
                icode,                                                         &
                cmessage,                                                      &
                typemap = obj_type_map,                                        &
                north = north,                                                 &
                south = south,                                                 &
                east = east,                                                   &
                west = west,                                                   &
                gridpoint = gridpoint,                                         &
                proccode = proccode,                                           &
                mpp_lookup = mpp_lookup)

  DEALLOCATE(proccode)
  DEALLOCATE(gridpoint)
  DEALLOCATE(west)
  DEALLOCATE(east)
  DEALLOCATE(south)
  DEALLOCATE(north)
  DEALLOCATE(obj_type_map)

  IF (printstatus  >=  prstatus_normal) THEN
    IF (mype  ==  0) THEN
      CALL umPrint('Data successfully read',src='um_readdump')
      WRITE(umMessage,'(I0,A,I0)') fixhd(161),' words read from unit ',nftin
      CALL umPrint(umMessage,src='um_readdump')
      IF ((fixhd(5)  >=  6) .AND. & ! AC/Var
          (fixhd(5)  <=  8)) THEN   ! Obs/Cx
        CALL umPrint('(Observational data)',src='um_readdump')
      ELSE
        CALL umPrint('(Model data)',src='um_readdump')
      END IF
    END IF ! IF (mype  ==  0)
  END IF ! IF (PrintStatus  >=  PrStatus_Normal)

END IF ! IF (FIXHD(160)  >   0)

! Reset to original decomposition type
CALL change_decomposition(orig_decomp,icode)

!!! Note - we are still in an openmp region in a thread==0 clause
!!!        thread zero now detatches helpers
!$   IF (IOS_use_helpers) THEN
!$     CALL portioDetachAllHelpers()
!$     CALL umPrint('um_readddump: read completed, detatch request completed', &
!$          src='um_readdump')
!$   END IF

!!! the top thread will be the helper
!$ ELSE IF (thread==omp_get_num_threads()-1) THEN
!$   IF (IOS_use_helpers) THEN
!$     CALL umPrint('um_readdump: Reader thread helper attaching', &
!$          src='um_readdump')
!$     CALL portioAddHelper(1)
!$   END IF
!!! other threads do nothing (yet)
!$ ELSE
!$   WRITE(umMessage,'(A,I3,A)')'um_readdump: thread',thread,' is inactive'
!$   CALL umPrint(umMessage,src='um_readdump')
!$ END IF

!$OMP END PARALLEL

IF (printstatus  >=  prstatus_diag) THEN

  ! Print a debug checksum string for each field.
  ! This will provide a value to allow for quick validation that the data was
  ! read in correctly, by comparison to a known good value for the fields in the
  ! dump.

  DO k = 1,number_of_fields  ! loop over fields to read in

    checksum_str = ''

    i = mpp_lookup(p_naddr,k)
    check_len = mpp_lookup(p_lblrec,k)

    IF (check_len>0 .AND. i>0) THEN

      ALLOCATE(d1_field(check_len))

      ! construct an array of 64-bit integers containing the field data bits.
      d1_field = TRANSFER(d1(i:i+check_len-1),[1_integer64])

      DO j = 0,63

        ! Create a 64-bit mask with only the j-th bit set.
        checksum_shft = ISHFT(INT(b'1',KIND=integer64),j)

        ! Create the checksum value by counting the number of times the j-th bit
        ! is set (has a value of one) in the field data.
        ! This is done by filtering the data using the IAND function to select
        ! only the j-th bits from the array.
        checksum_int = COUNT( IAND(d1_field,checksum_shft) == checksum_shft )

        ! Ensure the checksum value printed is in the range 1-99.
        ! This will prevent it overflowing 2 digits, but leave it
        ! distinguishable from the empty field case.
        WRITE(checksum_str,'(A,I2.2)') TRIM(checksum_str),1+MOD(checksum_int,99)

      END DO

      DEALLOCATE(d1_field)

    ELSE

      ! Empty field case - set all the checksum values to 0.
      checksum_int = 0
      WRITE(checksum_str,'(64("00"))')

    END IF

    WRITE(umMessage,'(A,I0,A,A)') 'Field ',k,' checksum = ',TRIM(checksum_str)
    CALL umPrint(umMessage,src='um_readdump')

  END DO ! K : loop over fields to read in
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE um_readdump
END MODULE um_readdump_mod

