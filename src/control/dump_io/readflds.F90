! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Reads in a number of fields from UM format file
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Dump I/O

MODULE readflds_mod

USE pr_look_mod, ONLY: pr_look

IMPLICIT NONE

PRIVATE
PUBLIC :: readflds

INTERFACE readflds
MODULE PROCEDURE                                                               &
  readflds1d,                                                                  &
  readflds2d,                                                                  &
  readflds3d
END INTERFACE

CHARACTER(LEN=*), PARAMETER :: ModuleName='READFLDS_MOD'

CONTAINS

!------------------------------------------------------------------------------!

! 3D array Subroutine Interface:
SUBROUTINE readflds3d(nftin, number_of_fields, first_field, lookup, d1, fixhd, &
                      expand, icode, cmessage, typemap, north, south, east,    &
                      west, gridpoint, proccode, mpp_lookup)

USE errormessagelength_mod, ONLY: errormessagelength

USE iso_c_binding, ONLY: C_LOC, C_F_POINTER

IMPLICIT NONE

! Subroutine Arguments:

INTEGER, INTENT(IN) ::                                                         &
  nftin,                & ! IN: unit number to read data from
  number_of_fields,     & ! IN: number of fields to read in
  first_field,          & ! IN: first field to read in
  fixhd(*),             & ! IN: fixed length header
  expand                  ! IN: (=1 if WGDOS or RLE packed data
                          !      is to be expanded)
                          ! Only used for small execs etc

INTEGER, INTENT(INOUT) ::                                                      &
  lookup(:,:)             ! IN: lookup table starting at field 1

INTEGER, INTENT(OUT) ::                                                        &
  icode                   ! OUT: return code

REAL, INTENT(INOUT), TARGET ::                                                 &
  d1(:,:,:)               ! INOUT: array to return the data in

CHARACTER(LEN=errormessagelength), INTENT(OUT) ::                              &
  cmessage                ! OUT: Error message if ICODE <> 0

INTEGER, INTENT(IN), OPTIONAL ::                                               &
  typemap(:),           & ! IN: defines whether field is prognostic, diagnostic
                          !     or other
  north(:),             & ! IN: defines the field's north code
  south(:),             & ! IN: defines the field's south code
  east(:),              & ! IN: defines the field's east code
  west(:),              & ! IN: defines the field's west code
  gridpoint(:),         & ! IN: defines the field's gridpoint code
  proccode(:)             ! IN: defines the field's processing code

INTEGER, INTENT(INOUT), OPTIONAL ::                                            &
  mpp_lookup(:,:)         ! INOUT: local PE lookup table

! Local variables
REAL, POINTER :: d1_remap(:)

  ! C_F_POINTER pointer voodoo (d1 must be contiguous)
  CALL C_F_POINTER(C_LOC(d1(1,1,1)), d1_remap, [SIZE(d1)])

  CALL readflds_core(nftin, number_of_fields, first_field, lookup, d1_remap,   &
                     fixhd, expand, icode, cmessage, typemap, north, south,    &
                     east, west, gridpoint, proccode, mpp_lookup)

END SUBROUTINE readflds3d

!------------------------------------------------------------------------------!

! 2D array Subroutine Interface:
SUBROUTINE readflds2d(nftin, number_of_fields, first_field, lookup, d1, fixhd, &
                      expand, icode, cmessage, typemap, north, south, east,    &
                      west, gridpoint, proccode, mpp_lookup)

USE errormessagelength_mod, ONLY: errormessagelength

USE iso_c_binding, ONLY: C_LOC, C_F_POINTER

IMPLICIT NONE

! Subroutine Arguments:

INTEGER, INTENT(IN) ::                                                         &
  nftin,                & ! IN: unit number to read data from
  number_of_fields,     & ! IN: number of fields to read in
  first_field,          & ! IN: first field to read in
  fixhd(*),             & ! IN: fixed length header
  expand                  ! IN: (=1 if WGDOS or RLE packed data
                          !      is to be expanded)
                          ! Only used for small execs etc

INTEGER, INTENT(INOUT) ::                                                      &
  lookup(:,:)             ! IN: lookup table starting at field 1

INTEGER, INTENT(OUT) ::                                                        &
  icode                   ! OUT: return code

REAL, INTENT(INOUT), TARGET ::                                                 &
  d1(:,:)                 ! INOUT: array to return the data in

CHARACTER(LEN=errormessagelength), INTENT(OUT) ::                              &
  cmessage                ! OUT: Error message if ICODE <> 0

INTEGER, INTENT(IN), OPTIONAL ::                                               &
  typemap(:),           & ! IN: defines whether field is prognostic, diagnostic
                          !     or other
  north(:),             & ! IN: defines the field's north code
  south(:),             & ! IN: defines the field's south code
  east(:),              & ! IN: defines the field's east code
  west(:),              & ! IN: defines the field's west code
  gridpoint(:),         & ! IN: defines the field's gridpoint code
  proccode(:)             ! IN: defines the field's processing code

INTEGER, INTENT(INOUT), OPTIONAL ::                                            &
  mpp_lookup(:,:)         ! INOUT: local PE lookup table

! Local variables
REAL, POINTER :: d1_remap(:)

  ! C_F_POINTER pointer voodoo (d1 must be contiguous)
  CALL C_F_POINTER(C_LOC(d1(1,1)), d1_remap, [SIZE(d1)])

  CALL readflds_core(nftin, number_of_fields, first_field, lookup, d1_remap,   &
                     fixhd, expand, icode, cmessage, typemap, north, south,    &
                     east, west, gridpoint,proccode, mpp_lookup)

END SUBROUTINE readflds2d

!------------------------------------------------------------------------------!

! 1D array Subroutine Interface:
SUBROUTINE readflds1d(nftin, number_of_fields, first_field, lookup, d1, fixhd, &
                      expand, icode, cmessage, typemap, north, south, east,    &
                      west, gridpoint,proccode, mpp_lookup)

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Subroutine Arguments:

INTEGER, INTENT(IN) ::                                                         &
  nftin,                & ! IN: unit number to read data from
  number_of_fields,     & ! IN: number of fields to read in
  first_field,          & ! IN: first field to read in
  fixhd(*),             & ! IN: fixed length header
  expand                  ! IN: (=1 if WGDOS or RLE packed data
                          !      is to be expanded)
                          ! Only used for small execs etc

INTEGER, INTENT(INOUT) ::                                                      &
  lookup(:,:)             ! IN: lookup table starting at field 1

INTEGER, INTENT(OUT) ::                                                        &
  icode                   ! OUT: return code

REAL, INTENT(INOUT), TARGET ::                                                 &
  d1(:)                   ! INOUT: array to return the data in

CHARACTER(LEN=errormessagelength), INTENT(OUT) ::                              &
  cmessage                ! OUT: Error message if ICODE <> 0

INTEGER, INTENT(IN), OPTIONAL ::                                               &
  typemap(:),           & ! IN: defines whether field is prognostic, diagnostic
                          !     or other
  north(:),             & ! IN: defines the field's north code
  south(:),             & ! IN: defines the field's south code
  east(:),              & ! IN: defines the field's east code
  west(:),              & ! IN: defines the field's west code
  gridpoint(:),         & ! IN: defines the field's gridpoint code
  proccode(:)             ! IN: defines the field's processing code

INTEGER, INTENT(INOUT), OPTIONAL ::                                            &
  mpp_lookup(:,:)         ! INOUT: local PE lookup table

  CALL readflds_core(nftin, number_of_fields, first_field, lookup, d1, fixhd,  &
                     expand, icode, cmessage, typemap, north, south, east,     &
                     west, gridpoint, proccode,mpp_lookup)

END SUBROUTINE readflds1d

!------------------------------------------------------------------------------!

! Core Subroutine Interface:
SUBROUTINE readflds_core(nftin, number_of_fields, first_field, lookup, d1,     &
                         fixhd, expand, icode, cmessage, typemap, north,       &
                         south, east, west, gridpoint,proccode, mpp_lookup)

USE stextend_mod, ONLY:                                                        &
    in_s

USE application_description, ONLY:                                             &
    isSmallExec

USE yomhook, ONLY:                                                             &
    lhook, dr_hook

USE parkind1, ONLY:                                                            &
    jprb, jpim

USE io, ONLY:                                                                  &
    isAllLocal, broadcast_read, setpos

USE io_configuration_mod, ONLY:                                                &
    io_alltoall_readflds

USE ereport_mod, ONLY:                                                         &
    ereport

USE umPrintMgr, ONLY:       &
    umPrint,                &
    umMessage,              &
    newline,                &
    PrintStatus,            &
    PrMin,                  &
    PrStatus_Oper,          &
    PrStatus_Diag

USE UM_ParVars

USE UM_ParCore, ONLY: mype, nproc

USE lookup_addresses, ONLY:                                                    &
    p_lblrec, lbpack, lbnpt, p_naddr, lbegin, lbrow, lbext, item_code,         &
    model_code, lbhem, lbyr, data_type, lblrec, lbnrec, lbproc

USE submodel_mod, ONLY:                                                        &
    atmos_im, atmos_sm

USE d1_array_mod, ONLY:                                                        &
    d1_object_type, d1_imodl, d1_section, d1_item,  d1_no_levels, prognostic,  &
    d1_addr

USE errormessagelength_mod, ONLY:                                              &
    errormessagelength

USE lbc_mod, ONLY:                                                             &
    g_lenrima

USE read_serial_mod, ONLY:                                                     &
    read_serial

USE d1_array_mod, ONLY:                                                        &
    diagnostic

USE cppxref_mod, ONLY:                                                         &
    ppx_atm_rim,                                                               &
    ppx_atm_compressed,                                                        &
    ppx_atm_cuall,                                                             &
    ppx_atm_cvall,                                                             &
    ppx_atm_tzonal,                                                            &
    ppx_atm_ozone,                                                             &
    ppx_atm_uzonal,                                                            &
    ppx_atm_umerid,                                                            &
    ppx_atm_tmerid,                                                            &
    ppx_atm_lbc_u,                                                             &
    ppx_atm_lbc_v,                                                             &
    ppx_atm_lbc_theta,                                                         &
    ppx_atm_lbc_orog,                                                          &
    ppx_atm_tall

USE atm_land_sea_mask, ONLY:                                                   &
    atmos_landmask_local,                                                      &
    atmos_landmask

#if !defined(UTILIO)

USE ppxlook_mod, ONLY:                                                         &
    exppxi

USE read_multi_mod, ONLY:                                                      &
    read_multi

USE mpl, ONLY:                                                                 &
    mpl_real8

#endif

USE nlsizes_namelist_mod, ONLY:                                                &
    n_obj_d1_max

USE rimtypes

USE mask_compression, ONLY :                                                   &
    compress_to_mask,                                                          &
    expand_from_mask

USE lbc_calc_size_mod, ONLY: lbc_calc_size

USE field_types, ONLY:                                                         &
    fld_type_p,                                                                &
    fld_type_u,                                                                &
    fld_type_v

USE um_parparams, ONLY:                                                        &
    halo_type_no_halo

USE um_parparams, ONLY:                                                        &
    pnorth, peast

USE atm_land_sea_mask, ONLY:                                                   &
    atmos_number_of_landpts_proc, atmos_number_of_landpts

USE Field_Types,  ONLY: Nfld_Max

USE sterr_mod, ONLY: st_no_data

USE stparam_mod, ONLY: st_time_series_code, st_time_series_mean

USE packing_codes_mod, ONLY:                                                   &
    PC_WGDOS_packing,                                                          &
    PC_RunLength_Packing,                                                      &
    PC_Cray32_packing,                                                         &
    PC_GRIB_Packing

IMPLICIT NONE

! Description:

!  Reads in NUMBER_OF_FIELDS fields from file on unit NFTIN,
!  starting at field number FIRST_FIELD. The data is returned
!  in the D1 array.

! Subroutine Arguments:

INTEGER, INTENT(IN) ::                                                         &
  nftin,                & ! IN: unit number to read data from
  number_of_fields,     & ! IN: number of fields to read in
  first_field,          & ! IN: first field to read in
  fixhd(*),             & ! IN: fixed length header
  expand                  ! IN: (=1 if WGDOS or RLE packed data
                          !      is to be expanded)
                          ! Only used for small execs etc

INTEGER, INTENT(INOUT) ::                                                      &
  lookup(:,:)             ! IN: lookup table starting at field 1

INTEGER, INTENT(OUT) ::                                                        &
  icode                   ! OUT: return code

REAL, INTENT(INOUT), TARGET ::                                                 &
  d1(:)                   ! INOUT: array to return the data in

CHARACTER(LEN=errormessagelength), INTENT(OUT) ::                              &
  cmessage                ! OUT: Error message if ICODE <> 0

INTEGER, INTENT(IN), OPTIONAL ::                                               &
  typemap(:),           & ! IN: defines whether field is prognostic, diagnostic
                          !     or other
  north(:),             & ! IN: defines the field's north code
  south(:),             & ! IN: defines the field's south code
  east(:),              & ! IN: defines the field's east code
  west(:),              & ! IN: defines the field's west code
  gridpoint(:),         & ! IN: defines the field's gridpoint code
  proccode(:)             ! IN: defines the field's processing code

INTEGER, INTENT(INOUT), OPTIONAL ::                                            &
  mpp_lookup(:,:)         ! INOUT: local PE lookup table

INTEGER, PARAMETER :: fh_lookupsize2 = 152

! Local variables

INTEGER ::                                                                     &
  k,                    & ! loop over fields to read in
  d1_off,               & ! local offset into D1 for this field
  pack_code,            & ! packing code for field
  field_start,          & ! location of field on disk
  data_size,            & ! number of words of data on disk
                          ! (including padding for WFIO)
  data_read_size,       & ! number of words to read from disk
  data_full_size,       & ! number of words after any unpacking
  len_io,               & ! number of words read from disk
  field_item,           & ! Item number of field
  field_sect,           & ! Section number of field
  field_model,          & ! Model ID of field
  grid_type,            & ! grid type code
  fld_type,             & ! field type (P,U or V)
  halo_type,            & ! halo type code
  i,                    & ! loop index
  local_len,            & ! size of field section on this PE
  ipseudo,              & ! pseudo type code
  im_ident,             & ! model number
  num_levels,           & ! number of levels to be read by read_multi
  field_length,         & ! Record length of field
  field_type,           &
  d1_start_read

REAL    ::                                                                     &
  a_io                 ! Return code from BUFFIN

REAL, ALLOCATABLE ::                                                           &
  send_buf(:), recv_buf(:)

LOGICAL ::                                                                     &
  serial,   & ! switch between serial and parallel reads:
              ! .TRUE. if read_serial() will be used instead of read_multi() or
              ! all-to-all reads. The default is .FALSE., unless UTILIO is
              ! defined, or there is only one MPI rank.
  alltoall, & ! switch between using the all-to-all method or not for parallel
              ! reads:
              ! .TRUE. if read_alltoall() will be used instead of read_multi()
              ! Default is .FALSE.
  subdomain, &
  abort_alltoall

! Per processor arrays
INTEGER ::                                                                     &
  send_size(0:nproc-1),                                                        &
  send_offset(0:nproc-1),                                                      &
  recv_size(0:nproc-1),                                                        &
  recv_offset(0:nproc-1),                                                      &
  start_field(0:nproc),                                                        &
  fields_per_proc(0:nproc-1),                                                  &
  address(0:nproc-1),                                                          &
  remote_send_size(0:nproc-1)

REAL, ALLOCATABLE ::                                                           &
  read_field(:),                                                               &
  unpacked_field(:)

INTEGER ::                                                                     &
  igpos, igpos1, ipos, j, maskpoints, xproc, yproc, my_comm, info, north_loc,  &
  south_loc, east_loc, west_loc, gridpoint_loc, proccode_loc,  typemap_loc,    &
  iside, full_lbc_row_len, full_lbc_nrows, rim_type, decomp_lbc_row_len,       &
  full_lbc_start_pt, decomp_lbc_nrows_no_halo, decomp_lbc_start_pt_no_halo,    &
  decomp_lbc_nrows, full_lbc_nrows_no_halo, decomp_lbc_row_len_no_halo, lbck,  &
  full_lbc_row_len_no_halo, igposl, full_lbc_start_pt_no_halo, igposg,         &
  decomp_lbc_start_pt, iproc, grid_type_proc, local_start_row, local_end_row,  &
  local_start_col, local_end_col

! Parameters
INTEGER, PARAMETER ::                                                          &
  unset=-1             ! unset values

CHARACTER(LEN=*), PARAMETER ::                                                 &
  RoutineName='READFLDS_CORE'

! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Functions

INTEGER :: get_fld_type

!--------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

abort_alltoall = .FALSE.

#if defined(UTILIO)

! Always read serially when UTILIO defined
serial = .TRUE.
alltoall = .FALSE.

#else

! By default, serial reading and all-to-all are off
serial = .FALSE.
alltoall = .FALSE.

! Use all to all if that option is selected, provided that this fuctionality
! is possible (which requires the file to be opened in ioAllLocal mode).
IF( io_alltoall_readflds .AND. isAllLocal(nftin) ) THEN

  alltoall = .TRUE.

ELSE

  ! If we only have one rank, or the file must be read in fully on every rank,
  ! then broadcast_read(nftin) will be .FALSE.
  ! In this case, revert to serial
  IF( .NOT.broadcast_read(nftin) ) THEN
    serial = .TRUE.
  END IF

END IF

IF ( PrintStatus>=PrStatus_Oper ) THEN

  IF( alltoall ) THEN
    WRITE(umMessage,'(A,I0,A)') 'READFLDS: reading ', number_of_fields,        &
                                ' fields with all-to-all algorithm'
    CALL umPrint(umMessage,src='readflds')
  ELSE IF( serial ) THEN
    WRITE(umMessage,'(A,I0,A)') 'READFLDS: reading ', number_of_fields,        &
                                ' fields with serial algorithm'
    CALL umPrint(umMessage,src='readflds')
  ELSE
    WRITE(umMessage,'(A,I0,A)') 'READFLDS: reading ', number_of_fields,        &
                                ' fields with multi algorithm'
    CALL umPrint(umMessage,src='readflds')
  END IF

END IF

#endif

! If an ancillary this may be generated by the new ancil program.
! For ancils the version has nothing to do with UM version.

IF (fixhd(12) <  403 .AND. fixhd(5) /= 4 ) THEN

  WRITE(umMessage,'(A,I0)') 'READFLDS: file created by UM version ', fixhd(12)
  CALL umPrint(umMessage,src='readflds')
  icode=1
  cmessage='READFLDS: Cannot read fields files from before UM vn4.3'

  CALL ereport ( routinename, icode, cmessage )

END IF

IF (alltoall) THEN

#if !defined(UTILIO)

  ! Distribute fields across all the constituent processors
  fields_per_proc(:) = number_of_fields / nproc
  DO iproc=1,MOD(number_of_fields, nproc)
    fields_per_proc(iproc-1) = fields_per_proc(iproc-1) + 1
  END DO

  IF(SUM(fields_per_proc) /= number_of_fields) THEN
    CALL umPrint("READFLDS_ALLTOALL:"//                                        &
                 " distributed fields does not match fields",                  &
                 src='readflds')
    icode = 11
    cmessage = "READFLDS_ALLTOALL: Field distribution has gone wrong"
    CALL ereport ( routinename, icode, cmessage )
  END IF

  ! Calculate the position of each PE's fields
  ! NB: start_field is one element larger than nproc to include the end point
  start_field(0) = first_field
  DO iproc=1,nproc
     start_field(iproc) = start_field(iproc-1) + fields_per_proc(iproc-1)
  END DO

  ! Calculate how much data to send to each processor
  ! and how much data is coming and going from/to each rank
  send_size(:) = 0
  recv_size(:) = 0
  i = 0

  ! calculate start d1_off on this rank
  d1_off=LBOUND(d1,dim=1)-1

  WRITE(umMessage,'(A,I0,A,I0,A,I0,A)') 'readflds[',mype,'] : reading fields ',&
                                        start_field(mype), ' to ',             &
                                        start_field(mype+1)-1, ' on this PE.'
  CALL umPrint(umMessage,src='readflds')

  ! Process fields read in by lower ranks than us

  DO k = first_field,start_field(mype)-1

    ! Reset local_len
    local_len = unset

    field_item=MOD(lookup(item_code,k),1000)
    field_sect=(lookup(item_code,k)-field_item)/1000
    field_model=lookup(model_code,k)

    CALL get_grid(fixhd(5),                                                    &
                  field_model,                                                 &
                  field_sect,                                                  &
                  field_item,                                                  &
                  lookup(lbhem,k),                                             &
                  lookup(lbpack,k),                                            &
                  lookup(lbnpt,k),                                             &
                  grid_type,                                                   &
                  halo_type,                                                   &
                  fld_type,                                                    &
                  icode)

    IF ( icode /= 0 ) THEN
      cmessage='READFLDS : get_grid call failed'
      GO TO 9999
    END IF

    ! Set the number of levels in the field

    IF (  grid_type == ppx_atm_lbc_theta                                       &
         .OR. grid_type == ppx_atm_lbc_u                                       &
         .OR.  grid_type == ppx_atm_lbc_v) THEN

      CALL get_lbc_levels(field_model,                                         &
                          field_sect,                                          &
                          field_item,                                          &
                          num_levels,                                          &
                          icode)

    ELSE

       num_levels = 1

    END IF

    CALL readflds_local_field_length(field_length,                             &
                                     grid_type,                                &
                                     fld_type,                                 &
                                     halo_type,                                &
                                     lookup(lbproc,k),                         &
                                     num_levels,                               &
                                     k,                                        &
                                     mype,                                     &
                                     .TRUE.,                                   &
                                     typemap,                                  &
                                     north,                                    &
                                     south,                                    &
                                     east,                                     &
                                     west,                                     &
                                     gridpoint,                                &
                                     st_no_data_test_fix=.TRUE.)


    IF (lookup(lblrec,k) > 0) THEN ! If there's data in the field

      ! Check that DATA_TYPE is valid no: +/-1 to +/-3
      IF (( ABS(lookup(data_type,k))  >=  1) .AND.                             &
          ( ABS(lookup(data_type,k))  <=  3)) THEN

        !   ------------------------
        !   Perform sanity checks
        !   ------------------------

        d1_start_read = d1_off + 1

        ! Check address is within array bounds
        IF (d1_start_read > UBOUND(d1,dim=1)) THEN

          ! for LBCs we need some special logic...
          IF ( (lookup(lbhem,k)  >=  99) .AND.                                 &
             (lookup(lbhem,k)  <   1000)) THEN ! This is a LBC field

            IF (field_model  /=  atmos_im) THEN
              icode=2
              WRITE(umMessage,'(A,I0)') 'READFLDS: Cannot process LBC for ' // &
                                      'model type ', field_model
              CALL umPrint(umMessage,src='readflds')
              cmessage='READFLDS : Cannot read LBCS for this model type'
              GO TO 9999
            END IF

          END IF

          IF ( field_length > 0) THEN
            icode=5
            WRITE(umMessage,'(4(A,I0))') "Readflds Address error for field ",  &
                                         k, ":" // newline //                  &
                                         "Calculated (address)" //             &
                                         "/(field_length): (", d1_start_read,  &
                                         ')/(', field_length, ')' // newline //&
                                         "Maximum D1 address: ",               &
                                         UBOUND(d1,dim=1)
            CALL ereport(routinename, icode, umMessage)
          ELSE

            ! A zero sized field at the end of the d1-array
            ! set the d1_start_address to equal the length
            ! of the d1 array to prevent out of bounds access
            d1_start_read = UBOUND(d1,dim=1) - 1

          END IF

          local_len = unset

        END IF

        ! Check field is included in in_s array (not applicable to small execs)
        IF (.NOT. isSmallExec()) THEN
          IF (in_s(1,field_model,field_sect,field_item) == 0) THEN
            WRITE(umMessage,'(A)') 'ERROR: Updated field not used in this ' // &
                                   'model run.'
            CALL umPrint(umMessage,src='readflds',level=PrMin)
            WRITE(umMessage,'(A)') 'Check your ancillary files are correct ' //&
                                   'as it may'
            CALL umPrint(umMessage,src='readflds',level=PrMin)
            WRITE(umMessage,'(A)') 'be that the updating has exceeded the ' // &
                                   'end of an'
            CALL umPrint(umMessage,src='readflds',level=PrMin)
            WRITE(umMessage,'(A)') 'ancillary file and it has gone onto ' //   &
                                   'the next'
            CALL umPrint(umMessage,src='readflds',level=PrMin)
            WRITE(umMessage,'(A)') 'field (which is this unacceptable one).' //&
                                   ' Look'
            CALL umPrint(umMessage,src='readflds',level=PrMin)
            WRITE(umMessage,'(A)') 'at the ancillary file associated with ' // &
                                   'the stash'
            CALL umPrint(umMessage,src='readflds',level=PrMin)
            WRITE(umMessage,'(A)') 'code before the following one:'
            CALL umPrint(umMessage,src='readflds',level=PrMin)
            WRITE(umMessage,'(A,I0)') 'Model ID : ',field_model
            CALL umPrint(umMessage,src='readflds',level=PrMin)
            WRITE(umMessage,'(A,I0)') 'Section  : ',field_sect
            CALL umPrint(umMessage,src='readflds',level=PrMin)
            WRITE(umMessage,'(A,I0)') 'Item     : ',field_item
            CALL umPrint(umMessage,src='readflds',level=PrMin)
            icode=101
            cmessage='Updated field not used in this model run. See output.'
            GO TO 9999
          END IF
        END IF

      ELSE ! Error in LOOKUP(DATA_TYPE,K)

        IF (( fixhd(5)  <   6) .OR.                                            &
            ( fixhd(5)  >   8)) THEN ! Not AC, Var or Cx

          CALL pr_look(lookup,k)

        END IF

        WRITE(umMessage,'(A,I0,A,I0)') 'readflds : failure for field ', k,     &
                                       ' of ', number_of_fields
        CALL umPrint(umMessage,src='readflds')
        WRITE(umMessage,'(A,I0)') 'LOOKUP(DATA_TYPE,K)= ',lookup(data_type,k)
        CALL umPrint(umMessage,src='readflds')
        icode=7
        cmessage='Invalid data type ( LOOKUP(DATA_TYPE,K) )'

        CALL ereport(routinename, icode, cmessage)

      END IF

    END IF ! If there was data in the field

    d1_off=d1_off+field_length

    ! ppx_atm_compressed fields will be unpacked before scattering
    IF (atmos_number_of_landpts_proc(mype) > 0) THEN
      IF (grid_type == ppx_atm_compressed) THEN
        SELECT CASE( fld_type )
          CASE ( fld_type_p )
            grid_type = ppx_atm_tall

          CASE ( fld_type_u )
            grid_type = ppx_atm_cuall

          CASE ( fld_type_v )
            grid_type = ppx_atm_cvall
        END SELECT
      END IF
    END IF

    ! recv_size includeds no halos, no temp_fix
    CALL readflds_local_field_length(field_length,                             &
                                     grid_type,                                &
                                     fld_type,                                 &
                                     halo_type,                                &
                                     lookup(lbproc,k),                         &
                                     num_levels,                               &
                                     k,                                        &
                                     mype,                                     &
                                     .FALSE.,                                  &
                                     typemap,                                  &
                                     north,                                    &
                                     south,                                    &
                                     east,                                     &
                                     west,                                     &
                                     gridpoint,                                &
                                     st_no_data_test_fix=.FALSE.)

    recv_size(i) = recv_size(i) + field_length

    ! Increment receiving processor when we reach its section of fields
    IF (start_field(i+1)-1 == k) THEN

      i = i + 1

    END IF

    IF (PRESENT(north)) THEN
      north_loc = north(k)
    ELSE
      north_loc = glsize(2,fld_type)
    END IF

    IF (PRESENT(south)) THEN
      south_loc = south(k)
    ELSE
      south_loc = 1
    END IF

    IF (PRESENT(east)) THEN
      east_loc = east(k)
    ELSE
      east_loc = glsize(1,fld_type)
    END IF

    IF (PRESENT(west)) THEN
      west_loc = west(k)
    ELSE
      west_loc = 1
    END IF

    IF (PRESENT(gridpoint)) THEN
      gridpoint_loc = gridpoint(k)
    ELSE
      gridpoint_loc = unset
    END IF

    IF (PRESENT(typemap)) THEN
      typemap_loc = typemap(k)
    ELSE
      typemap_loc = unset
    END IF

    subdomain = field_issubdomain(fld_type,                                    &
                                  lookup(lbproc,k),                            &
                                  typemap_loc,                                 &
                                  north_loc,                                   &
                                  south_loc,                                   &
                                  east_loc,                                    &
                                  west_loc,                                    &
                                  gridpoint_loc)

    IF (subdomain) THEN
      abort_alltoall = .TRUE.
    END IF

    IF (PRESENT(proccode)) THEN
      proccode_loc=proccode(k)
    ELSE
      proccode_loc=unset
    END IF

    IF ((typemap_loc  ==  diagnostic) .AND.                                    &
        ((proccode_loc  ==  st_time_series_code) .OR.                          &
         (proccode_loc  ==  st_time_series_mean))) THEN

      abort_alltoall = .TRUE.

    END IF

  END DO

  ! Process fields read in by us

  DO k = start_field(mype),start_field(mype+1)-1

    field_item=MOD(lookup(item_code,k),1000)
    field_sect=(lookup(item_code,k)-field_item)/1000
    field_model=lookup(model_code,k)

    CALL get_grid(fixhd(5),                                                    &
                  field_model,                                                 &
                  field_sect,                                                  &
                  field_item,                                                  &
                  lookup(lbhem,k),                                             &
                  lookup(lbpack,k),                                            &
                  lookup(lbnpt,k),                                             &
                  grid_type,                                                   &
                  halo_type,                                                   &
                  fld_type,                                                    &
                  icode)

    IF ( icode /= 0 ) THEN
      cmessage='READFLDS : get_grid call failed'
      GO TO 9999
    END IF

    ! Check if LBC, and if so, if its type is valid
    IF ( (lookup(lbhem,k)  >=  99) .AND.                                       &
         (lookup(lbhem,k)  <   1000)) THEN ! This is a LBC field

      IF (field_model /= atmos_im) THEN
        icode=2
        WRITE(umMessage,'(A,I0)') 'READFLDS: Cannot process LBC for model' //  &
                                  ' type ',field_model
        CALL umPrint(umMessage,src='readflds')
        cmessage='READFLDS : Cannot read LBCS for this model type'
        GO TO 9999
      END IF

    END IF

    ! Set the number of levels in the field

    IF ( grid_type == ppx_atm_lbc_theta                                        &
         .OR. grid_type == ppx_atm_lbc_u                                       &
         .OR.  grid_type == ppx_atm_lbc_v) THEN

      CALL get_lbc_levels(field_model,                                         &
                          field_sect,                                          &
                          field_item,                                          &
                          num_levels,                                          &
                          icode)

    ELSE

      num_levels = 1

    END IF

    IF(lookup(lblrec,k) > 0) THEN

      DO iproc = 0, nproc-1

       ! ppx_atm_compressed fields will be unpacked before scattering
        IF (grid_type == ppx_atm_compressed .AND.                              &
            atmos_number_of_landpts_proc(iproc) /= 0) THEN
          SELECT CASE( fld_type )
            CASE ( fld_type_p )
              grid_type_proc = ppx_atm_tall

            CASE ( fld_type_u )
              grid_type_proc = ppx_atm_cuall

            CASE ( fld_type_v )
              grid_type_proc = ppx_atm_cvall
          END SELECT
        ELSE
          grid_type_proc = grid_type
        END IF

        CALL readflds_local_field_length(field_length,                         &
                                         grid_type_proc,                       &
                                         fld_type,                             &
                                         halo_type,                            &
                                         lookup(lbproc,k),                     &
                                         num_levels,                           &
                                         k,                                    &
                                         iproc,                                &
                                         .FALSE.,                              &
                                         typemap,                              &
                                         north,                                &
                                         south,                                &
                                         east,                                 &
                                         west,                                 &
                                         gridpoint,                            &
                                         st_no_data_test_fix=.FALSE.)

        send_size(iproc) = send_size(iproc) + field_length

      END DO

      ! ppx_atm_compressed fields will be unpacked before scattering
      IF (atmos_number_of_landpts_proc(mype) > 0) THEN
        IF (grid_type == ppx_atm_compressed) THEN
          SELECT CASE( fld_type )
            CASE ( fld_type_p )
              grid_type = ppx_atm_tall

            CASE ( fld_type_u )
              grid_type = ppx_atm_cuall

            CASE ( fld_type_v )
              grid_type = ppx_atm_cvall
          END SELECT
        END IF
      END IF

      CALL readflds_local_field_length(field_length,                           &
                                       grid_type,                              &
                                       fld_type,                               &
                                       halo_type,                              &
                                       lookup(lbproc,k),                       &
                                       num_levels,                             &
                                       k,                                      &
                                       mype,                                   &
                                       .FALSE.,                                &
                                       typemap,                                &
                                       north,                                  &
                                       south,                                  &
                                       east,                                   &
                                       west,                                   &
                                       gridpoint,                              &
                                       st_no_data_test_fix=.FALSE.)

      recv_size(i) = recv_size(i) + field_length

    END IF

    ! Increment receiving processor when we reach its section of fields
    IF (start_field(i+1)-1 == k) THEN
      i = i + 1
    END IF

    IF (PRESENT(north)) THEN
      north_loc = north(k)
    ELSE
      north_loc = glsize(2,fld_type)
    END IF

    IF (PRESENT(south)) THEN
      south_loc = south(k)
    ELSE
      south_loc = 1
    END IF

    IF (PRESENT(east)) THEN
      east_loc = east(k)
    ELSE
      east_loc = glsize(1,fld_type)
    END IF

    IF (PRESENT(west)) THEN
      west_loc = west(k)
    ELSE
      west_loc = 1
    END IF

    IF (PRESENT(gridpoint)) THEN
      gridpoint_loc = gridpoint(k)
    ELSE
      gridpoint_loc = unset
    END IF

    IF (PRESENT(typemap)) THEN
      typemap_loc = typemap(k)
    ELSE
      typemap_loc = unset
    END IF

    subdomain = field_issubdomain(fld_type,                                    &
                                  lookup(lbproc,k),                            &
                                  typemap_loc,                                 &
                                  north_loc,                                   &
                                  south_loc,                                   &
                                  east_loc,                                    &
                                  west_loc,                                    &
                                  gridpoint_loc)

    IF (subdomain) THEN
      abort_alltoall = .TRUE.
    END IF

    IF (PRESENT(proccode)) THEN
      proccode_loc=proccode(k)
    ELSE
      proccode_loc=unset
    END IF

    IF ((typemap_loc  ==  diagnostic) .AND.                                    &
        ((proccode_loc  ==  st_time_series_code) .OR.                          &
         (proccode_loc  ==  st_time_series_mean))) THEN

      abort_alltoall = .TRUE.

    END IF

  END DO

  ! Process fields read in by higher ranks than us

  DO k = start_field(mype+1), first_field+number_of_fields-1

    field_item=MOD(lookup(item_code,k),1000)
    field_sect=(lookup(item_code,k)-field_item)/1000
    field_model=lookup(model_code,k)

    CALL get_grid(fixhd(5),                                                    &
                  field_model,                                                 &
                  field_sect,                                                  &
                  field_item,                                                  &
                  lookup(lbhem,k),                                             &
                  lookup(lbpack,k),                                            &
                  lookup(lbnpt,k),                                             &
                  grid_type,                                                   &
                  halo_type,                                                   &
                  fld_type,                                                    &
                  icode)

    IF ( icode /= 0 ) THEN
      cmessage='READFLDS : get_grid call failed'
      GO TO 9999
    END IF

    ! ppx_atm_compressed fields will be unpacked before scattering
    IF (atmos_number_of_landpts_proc(mype) > 0) THEN
      IF (grid_type == ppx_atm_compressed) THEN
        SELECT CASE( fld_type )
          CASE ( fld_type_p )
            grid_type = ppx_atm_tall

          CASE ( fld_type_u )
            grid_type = ppx_atm_cuall

          CASE ( fld_type_v )
            grid_type = ppx_atm_cvall
        END SELECT
      END IF
    END IF

    ! Set the number of levels in the field

    IF (  grid_type == ppx_atm_lbc_theta                                       &
         .OR. grid_type == ppx_atm_lbc_u                                       &
         .OR.  grid_type == ppx_atm_lbc_v) THEN

      CALL get_lbc_levels(field_model,                                         &
                          field_sect,                                          &
                          field_item,                                          &
                          num_levels,                                          &
                          icode)

    ELSE

       num_levels = 1

    END IF

    CALL readflds_local_field_length(field_length,                             &
                                     grid_type,                                &
                                     fld_type,                                 &
                                     halo_type,                                &
                                     lookup(lbproc,k),                         &
                                     num_levels,                               &
                                     k,                                        &
                                     mype,                                     &
                                     .FALSE.,                                  &
                                     typemap,                                  &
                                     north,                                    &
                                     south,                                    &
                                     east,                                     &
                                     west,                                     &
                                     gridpoint,                                &
                                     st_no_data_test_fix=.FALSE.)

    recv_size(i) = recv_size(i) + field_length

    ! Increment receiving processor when we reach its section of fields
    IF (start_field(i+1)-1 == k) THEN
      i = i + 1
    END IF

    IF (PRESENT(north)) THEN
      north_loc = north(k)
    ELSE
      north_loc = glsize(2,fld_type)
    END IF

    IF (PRESENT(south)) THEN
      south_loc = south(k)
    ELSE
      south_loc = 1
    END IF

    IF (PRESENT(east)) THEN
      east_loc = east(k)
    ELSE
      east_loc = glsize(1,fld_type)
    END IF

    IF (PRESENT(west)) THEN
      west_loc = west(k)
    ELSE
      west_loc = 1
    END IF

    IF (PRESENT(gridpoint)) THEN
      gridpoint_loc = gridpoint(k)
    ELSE
      gridpoint_loc = unset
    END IF

    IF (PRESENT(typemap)) THEN
      typemap_loc = typemap(k)
    ELSE
      typemap_loc = unset
    END IF

    subdomain = field_issubdomain(fld_type,                                    &
                                  lookup(lbproc,k),                            &
                                  typemap_loc,                                 &
                                  north_loc,                                   &
                                  south_loc,                                   &
                                  east_loc,                                    &
                                  west_loc,                                    &
                                  gridpoint_loc)

    IF (subdomain) THEN
      abort_alltoall = .TRUE.
    END IF

    IF (PRESENT(proccode)) THEN
      proccode_loc=proccode(k)
    ELSE
      proccode_loc=unset
    END IF

    IF ((typemap_loc  ==  diagnostic) .AND.                                    &
        ((proccode_loc  ==  st_time_series_code) .OR.                          &
         (proccode_loc  ==  st_time_series_mean))) THEN

         abort_alltoall = .TRUE.

    END IF

  END DO

  ! Calculate the offsets based on the data sizes
  send_offset(0) = 0
  recv_offset(0) = 0
  ! Use the address array to keep track of where data should be copied
  ! to in the send array
  address(0) = 1

  ! Calculate my send and receive size to each rank, starting with PE 0

  IF (mype == 0) THEN
    remote_send_size(0:nproc-1) = send_size(0:nproc-1)
  END IF

  ! confirm we agree on the sizes with other ranks

  CALL gc_ibcast(1,nproc,0,nproc,info,remote_send_size)

  IF (recv_size(0) /= remote_send_size(mype)) THEN
    WRITE(umMessage,'(A,I0,A,I0,A,I0,A)') 'readflds[', mype,                   &
                                          '] : [ERROR]: Send size [',          &
                                          remote_send_size(mype),              &
                                          '] /= receive size [',               &
                                          recv_size(0), '] from PE 0'
    CALL umPrint(umMessage,src='readflds')
    icode = 3011
    cmessage = 'READFLDS: Send/receive sizes do not match'
    CALL ereport ( routinename, icode, cmessage )
  END IF

  ! Now calculate my send and receive size for PEs > 0

  DO iproc = 1,nproc-1

    recv_offset(iproc) = recv_offset(iproc-1) + recv_size(iproc-1)
    send_offset(iproc) = send_offset(iproc-1) + send_size(iproc-1)
    address(iproc) = send_offset(iproc)+1

    IF (mype == iproc) THEN

      remote_send_size(0:nproc-1) = send_size(0:nproc-1)

    END IF

    ! confirm we agree on the sizes with other ranks

    CALL gc_ibcast(1+iproc,nproc,iproc,nproc,info,remote_send_size)

    IF (recv_size(iproc) /= remote_send_size(mype) ) THEN
      WRITE(umMessage,'(A,I0,A,I0,A,I0,A,I0)') 'readflds[', mype,              &
                                               '] : [ERROR]: Send size [',     &
                                               remote_send_size(mype),         &
                                               '] /= recieve size [',          &
                                               recv_size(iproc), '] from PE ', &
                                               iproc
      CALL umPrint(umMessage,src='readflds')
      icode = 3011
      cmessage = 'READFLDS: Send/receive sizes do not match'
      CALL ereport ( routinename, icode, cmessage )
    END IF

  END DO

  ! Allocate the send and receive buffers
  ALLOCATE(send_buf(SUM(send_size(:))))
  ALLOCATE(recv_buf(SUM(recv_size(:))))

#endif

END IF

IF (abort_alltoall) THEN

  ! Deallocate the temporary buffers created
  DEALLOCATE(recv_buf)
  DEALLOCATE(send_buf)

  ! Reset the method
  alltoall = .FALSE.
  IF( .NOT.broadcast_read(nftin) ) THEN
    serial = .TRUE.
  END IF

  ! Warn we have changed method.
  icode = -23
  cmessage = 'READFLDS: alltoall was selected, but a field being read is ' //  &
             'currently un-handled. Resetting read method.'
  CALL ereport ( routinename, icode, cmessage )

END IF

IF (.NOT.alltoall) THEN

  ! All ranks read every field, so...

  ! Set the start_field of mype to the first field
  start_field(mype) = first_field

  ! Set the end point (the start of the next pe) to be the total fields
  start_field(mype+1) = first_field+number_of_fields


  d1_off=LBOUND(d1,dim=1)-1

END IF

DO k = start_field(mype),start_field(mype+1)-1

  ! Reset local_len
  local_len = unset

  IF (PRESENT(mpp_lookup) .AND. .NOT.alltoall) THEN
    mpp_lookup(p_lblrec,k) = 0
    mpp_lookup(p_naddr,k)  = d1_off + 1
  END IF

  IF (lookup(lblrec,k)  >   0) THEN ! If there's data in
                                    ! the field

    ! Check that DATA_TYPE is valid no: +/-1 to +/-3
    IF (( ABS(lookup(data_type,k))  >=  1) .AND.                               &
        ( ABS(lookup(data_type,k))  <=  3)) THEN

      ! Get some information about this field

      field_item=MOD(lookup(item_code,k),1000)
      field_sect=(lookup(item_code,k)-field_item)/1000
      field_model=lookup(model_code,k)
      IF (PRESENT(typemap)) THEN
        field_type=typemap(k)
      ELSE
        field_type=prognostic
      END IF

      pack_code=MOD((lookup(lbpack,k)),10)

      !   ---------------------------------------
      !   Determine location of the field on disk
      !   ---------------------------------------

      ! Set up the location of the field on disk and how much data
      ! needs to be read in
      field_start=lookup(lbegin,k) ! position of field in file

      IF (field_start <= 0) THEN
        WRITE(umMessage,'(A,I0)') 'READFLDS: start address =',field_start
        CALL umPrint(umMessage,src='readflds')
        icode = 20
        cmessage = 'READFLDS: start address of field not given'
        CALL ereport ( routinename, icode, cmessage )
      END IF

      !   --------------------
      !   Determine data sizes
      !   --------------------

      ! DATA_SIZE : contains the number of words of data used to store the
      ! field on disk (needs to be halved if 32 bit packing has been used)

      IF (pack_code == PC_Cray32_packing) THEN
        data_size=(lookup(lblrec,k)+1)/2
      ELSE
        data_size=lookup(lblrec,k)
      END IF

      ! DATA_FULL_SIZE : is the number of words required to store the field
      ! in memory after any unpacking is done.

      ! This is to give buf the correct size in RDUNPCK, as
      ! buf will be the final expanded size of whole field
      ! including extra data
      ! warning lbext - may be -32768 missing value !

      IF ((pack_code == PC_RunLength_packing) .AND. (lookup(lbext, k) > 0)) THEN
        data_full_size=MAX(lookup(lbrow, k)*lookup(lbnpt, k)+lookup(lbext, k), &
                           lookup(lblrec,k))
      ELSE
        data_full_size=MAX(lookup(lbrow, k)*lookup(lbnpt, k),                  &
                           lookup(lblrec,k))
      END IF

      IF ((lookup(lbrow,k) <  0) .OR. (lookup(lbnpt,k) <  0)) THEN
        data_full_size=lookup(lblrec,k)
      END IF

      ! data_read_size contains the number of words to data that need to
      ! be read in for a field. Each field has extra words of dummy data
      ! added at the end to ensure each field starts on a disk sector
      ! boundary.
      data_read_size = lookup(lbnrec,k)

      ! The last field on a dump does not have these extra words
      ! added. So check against number of lookups and if next lookup is -99
      ! assume we are at the end.  Otherwise assume we are the last field.
      IF (k < fixhd(fh_lookupsize2)) THEN
        IF (lookup(lbyr,k+1) == -99) THEN
          data_Read_Size = data_size
        END IF
      ELSE
        data_read_size = data_size
      END IF

      IF (data_read_size < 0) THEN
        WRITE(umMessage,'(A,I0)') 'READFLDS: number of words to read =',       &
                                  data_read_size
        CALL umPrint(umMessage,src='readflds')
        icode = 30
        cmessage = 'READFLDS: number of words to read not given'
        CALL ereport ( routinename, icode, cmessage )
      END IF

      ! data_full_size needs to be at least as big as data_read_size since
      ! it is used to dimension the BUF array in READ_MULTI.

      data_full_size = MAX(data_full_size, data_read_size)

      !   -------------------------------------------
      !   Move file pointer to the start of the field
      !   -------------------------------------------

      CALL setpos(nftin,field_start,icode)

      IF (icode  /=  0) THEN
        WRITE(umMessage,'(A,I0,A,I0)') 'READFLDS - SETPOS failed to move ' //  &
                                       'file pointer to ', field_start,        &
                                       ' on unit ', nftin
        CALL umPrint(umMessage,src='readflds')
        WRITE(umMessage,'(A,I0)') 'SETPOS returned error code ',icode
        CALL umPrint(umMessage,src='readflds')
        icode=5
        cmessage='SETPOS failed while reading dump. See output.'

        CALL ereport ( routinename, icode, cmessage )
      END IF

      d1_start_read = d1_off + 1

      !   ------------------------
      !   Perform sanity checks
      !   ------------------------

      ! Check address is within array bounds
      IF (d1_start_read > UBOUND(d1,dim=1)) THEN

        CALL get_grid(fixhd(5),                                                &
                      field_model,                                             &
                      field_sect,                                              &
                      field_item,                                              &
                      lookup(lbhem,k),                                         &
                      lookup(lbpack,k),                                        &
                      lookup(lbnpt,k),                                         &
                      grid_type,                                               &
                      halo_type,                                               &
                      fld_type,                                                &
                      icode)


        IF ( icode /= 0 ) THEN
            cmessage='READFLDS : get_grid call failed'
            GO TO 9999
        END IF

        num_levels=1

        ! for LBCs we need some special logic...
        IF ( (lookup(lbhem,k)  >=  99) .AND.                                   &
           (lookup(lbhem,k)  <   1000)) THEN ! This is a LBC field

          IF (field_model  ==  atmos_im) THEN
            IF (lookup(lbhem,k)  /=  99) THEN
              ! New style LBCs with different field types
              ! See UMDP C04

              CALL get_lbc_levels(field_model,                                 &
                                  field_sect,                                  &
                                  field_item,                                  &
                                  num_levels,                                  &
                                  icode)

            END IF
          ELSE
            icode=2
            WRITE(umMessage,'(A,I0)') 'READFLDS: Cannot process LBC for ' //   &
                                      'model type ', field_model
            CALL umPrint(umMessage,src='readflds')
            cmessage='READFLDS : Cannot read LBCS for this model type'
            GO TO 9999
          END IF

        END IF

        CALL readflds_local_field_length(local_len,                            &
                                         grid_type,                            &
                                         fld_type,                             &
                                         halo_type,                            &
                                         lookup(lbproc,k),                     &
                                         num_levels,                           &
                                         k,                                    &
                                         mype,                                 &
                                         .TRUE.,                               &
                                         typemap,                              &
                                         north,                                &
                                         south,                                &
                                         east,                                 &
                                         west,                                 &
                                         gridpoint,                            &
                                         st_no_data_test_fix=.TRUE.)

        IF ( local_len > 0) THEN
          icode=5
          WRITE(umMessage,'(4(A,I0))') "Readflds Address error for field ", k, &
                                       ":"//newline//"Calculated (address)"//  &
                                       "/(field_length): (", d1_start_read,    &
                                       ')/(', local_len, ')'//newline//        &
                                       "Maximum D1 address: ", UBOUND(d1,dim=1)
          CALL ereport(routinename, icode, umMessage)
        ELSE
          ! A zero sized field at the end of the d1-array
          ! set the d1_start_address to equal the length
          ! of the d1 array to prevent out of bounds access
          d1_start_read = UBOUND(d1,dim=1) - 1

        END IF

        local_len = unset

      END IF

      ! Check field is included in in_s array (not applicable to small execs)
      IF (.NOT. isSmallExec()) THEN
        IF (in_s(1,field_model,field_sect,field_item) == 0) THEN
          WRITE(umMessage,'(A)') 'ERROR: Updated field not used in this ' //   &
                                 'model run.'
          CALL umPrint(umMessage,src='readflds',level=PrMin)
          WRITE(umMessage,'(A)') 'Check your ancillary files are correct ' //  &
                                 'as it may'
          CALL umPrint(umMessage,src='readflds',level=PrMin)
          WRITE(umMessage,'(A)') 'be that the updating has exceeded the ' //   &
                                 'end of an'
          CALL umPrint(umMessage,src='readflds',level=PrMin)
          WRITE(umMessage,'(A)') 'ancillary file and it has gone onto the next'
          CALL umPrint(umMessage,src='readflds',level=PrMin)
          WRITE(umMessage,'(A)') 'field (which is this unacceptable one). Look'
          CALL umPrint(umMessage,src='readflds',level=PrMin)
          WRITE(umMessage,'(A)') 'at the ancillary file associated with ' //   &
                                 'the stash'
          CALL umPrint(umMessage,src='readflds',level=PrMin)
          WRITE(umMessage,'(A)') 'code before the following one:'
          CALL umPrint(umMessage,src='readflds',level=PrMin)
          WRITE(umMessage,'(A,I0)') 'Model ID : ',field_model
          CALL umPrint(umMessage,src='readflds',level=PrMin)
          WRITE(umMessage,'(A,I0)') 'Section  : ',field_sect
          CALL umPrint(umMessage,src='readflds',level=PrMin)
          WRITE(umMessage,'(A,I0)') 'Item     : ',field_item
          CALL umPrint(umMessage,src='readflds',level=PrMin)
          icode=101
          cmessage='Updated field not used in this model run. See output.'
          GO TO 9999
        END IF
      END IF

      !   ------------------------
      !   Read the field from disk
      !   ------------------------

      IF (serial) THEN

        ! Each rank reads all fields serially

        CALL read_serial(nftin, d1(d1_off+1:), data_read_size, data_full_size, &
                         len_io, lookup(:,k), fixhd(12), field_model,          &
                         field_sect, field_item, icode, cmessage, expand)

#if !defined(UTILIO)

      ELSE IF (alltoall) THEN

        ! Read a different fraction of the fields serially on each rank, then
        ! distribute to each rank using all-to-all later on

        ALLOCATE(read_field(data_full_size))

        CALL read_serial(nftin, read_field(:), data_read_size, data_full_size, &
                         len_io, lookup(:,k), fixhd(12), field_model,          &
                         field_sect, field_item, icode, cmessage, expand)

        IF ( icode /= 0 ) THEN
          cmessage='READFLDS : all-to-all read_serial call failed'
          GO TO 9999
        END IF

        CALL get_grid(fixhd(5),                                                &
                      field_model,                                             &
                      field_sect,                                              &
                      field_item,                                              &
                      lookup(lbhem,k),                                         &
                      lookup(lbpack,k),                                        &
                      lookup(lbnpt,k),                                         &
                      grid_type,                                               &
                      halo_type,                                               &
                      fld_type,                                                &
                      icode)

        IF ( icode /= 0 ) THEN
          cmessage='READFLDS : get_grid call failed'
          GO TO 9999
        END IF

        ! Set the number of levels in the field

        IF (  grid_type == ppx_atm_lbc_theta                                   &
            .OR. grid_type == ppx_atm_lbc_u                                    &
            .OR.  grid_type == ppx_atm_lbc_v) THEN

          CALL get_lbc_levels(field_model,                                     &
                              field_sect,                                      &
                              field_item,                                      &
                              num_levels,                                      &
                              icode)

        ELSE

          num_levels = 1

        END IF

        CALL readflds_local_field_length(field_length,                         &
                                         grid_type,                            &
                                         fld_type,                             &
                                         halo_type,                            &
                                         lookup(lbproc,k),                     &
                                         num_levels,                           &
                                         k,                                    &
                                         mype,                                 &
                                         .TRUE.,                               &
                                         typemap,                              &
                                         north,                                &
                                         south,                                &
                                         east,                                 &
                                         west,                                 &
                                         gridpoint,                            &
                                         st_no_data_test_fix=.TRUE. )

        d1_off=d1_off+field_length

        IF (PRESENT(north)) THEN
          north_loc = north(k)
        ELSE
          north_loc = glsize(2,fld_type)
        END IF

        IF (PRESENT(south)) THEN
          south_loc = south(k)
        ELSE
          south_loc = 1
        END IF

        IF (PRESENT(east)) THEN
          east_loc = east(k)
        ELSE
          east_loc = glsize(1,fld_type)
        END IF

        IF (PRESENT(west)) THEN
          west_loc = west(k)
        ELSE
          west_loc = 1
        END IF

        IF (PRESENT(gridpoint)) THEN
          gridpoint_loc = gridpoint(k)
        ELSE
          gridpoint_loc = unset
        END IF

        IF (PRESENT(typemap)) THEN
          typemap_loc = typemap(k)
        ELSE
          typemap_loc = unset
        END IF

        subdomain = field_issubdomain(fld_type,                                &
                                      lookup(lbproc,k),                        &
                                      typemap_loc,                             &
                                      north_loc,                               &
                                      south_loc,                               &
                                      east_loc,                                &
                                      west_loc,                                &
                                      gridpoint_loc)


        ! pack the field into the send buffer - this does not include halos.

        IF (subdomain) THEN

          ! We should not be able to get here yet.
          icode=1
          cmessage = 'READFLDS: ' //                                           &
                     'Reached a subdomained field with all-to-all enabled.'
          CALL ereport(routinename, icode, cmessage)

        ELSE

          SELECT CASE( grid_type )
            CASE ( ppx_atm_tzonal, ppx_atm_uzonal )

              ! All zonal fields to scatter should have no halos. Throw an error
              ! if this isn't the case
              IF (halo_type /= halo_type_no_halo) THEN
                icode = 1
                WRITE(cmessage, '(a)') 'Cannot scatter zonal fields with a halo'
                CALL ereport(routinename, icode, cmessage)
              END IF

              igpos = 0
              DO yproc = 0, nproc_y-1
                DO xproc = 0, nproc_x-1
                  iproc = (yproc * nproc_x) + xproc
                  ipos = address(iproc)
                  DO i=1,g_blsize(2, fld_type, iproc)
                    send_buf(ipos) = read_field(igpos+i)
                    ipos = ipos + 1
                  END DO
                  address(iproc) = ipos
                END DO
                igpos = igpos + g_blsize(2, fld_type, yproc*nproc_x)
              END DO

            CASE ( ppx_atm_tmerid, ppx_atm_umerid )

              ! All meridional fields to scatter should have no halos. Throw an
              ! error if this isn't the case
              IF (halo_type /= halo_type_no_halo) THEN
                icode = 1
                WRITE(cmessage, '(a)') 'Cannot scatter meridional fields'//    &
                                     ' with a halo'
                CALL ereport(routinename, icode, cmessage)
              END IF

              igpos = 0
              DO xproc = 0, nproc_x-1
                DO yproc = 0, nproc_y-1
                  iproc = (yproc * nproc_x) + xproc
                  ipos = address(iproc)
                  DO i=1,g_blsize(1, fld_type, iproc)
                    send_buf(ipos) = read_field(igpos+i)
                    ipos = ipos + 1
                  END DO
                  address(iproc) = ipos
                END DO
                igpos = igpos + g_blsize(1, fld_type, xproc)
              END DO

            CASE ( ppx_atm_rim, ppx_atm_lbc_theta, ppx_atm_lbc_u,              &
                   ppx_atm_lbc_v, ppx_atm_lbc_orog )

              IF (grid_type  ==  ppx_atm_lbc_orog) THEN
                rim_type=rima_type_orog
              ELSE
                rim_type=rima_type_norm
              END IF

              DO iproc=0, nproc-1

                ipos = address(iproc) ! start of array

                DO iside=1,4
                  IF (g_at_extremity(iside,iproc)) THEN

                    ! This processor is at edge type iside and so needs LBC data
                    ! In order to copy data from full_lbc array to data_lbc
                    ! (that contains data to be scattered) the following
                    ! variables need to be calculated

                    ! full_lbc_row_len: East-West dimension of the full lbc side
                    ! full_lbc_nrows: North-South dimension of the full lbc side
                    ! decomp_lbc_row_len: East-West dimension of decomp lbc side
                    ! decomp_lbc_nrows: North-South dimension of decomp lbc side
                    ! full_lbc_start_pt: First point of the decomposed lbc side
                    !                     inside the 1d full_lbc array
                    ! decomp_lbc_start_pt: First point of the decomposed lbc
                    !                      side inside the 1d decomposed lbc
                    !                      array
                    ! data_lbc_start: First point inside the data lbc array
                    !                 (corresponds to the first point of the
                    !                  distributed lbc_comp)


                    CALL lbc_calc_size(iside,                                  &
                                       full_lbc_row_len,                       &
                                       full_lbc_nrows,                         &
                                       decomp_lbc_row_len,                     &
                                       decomp_lbc_nrows,                       &
                                       full_lbc_start_pt,                      &
                                       decomp_lbc_start_pt,                    &
                                       fld_type,                               &
                                       halo_type,                              &
                                       rim_type,                               &
                                       iproc)

                    CALL lbc_calc_size(iside,                                  &
                                       full_lbc_row_len_no_halo,               &
                                       full_lbc_nrows_no_halo,                 &
                                       decomp_lbc_row_len_no_halo,             &
                                       decomp_lbc_nrows_no_halo,               &
                                       full_lbc_start_pt_no_halo,              &
                                       decomp_lbc_start_pt_no_halo,            &
                                       fld_type,                               &
                                       halo_type_no_halo,                      &
                                       rim_type,                               &
                                       iproc)

                    DO lbck=0,num_levels-1

                      ! Skip the already done levels
                      igposl = lbck * decomp_lbc_nrows * decomp_lbc_row_len

                      DO j = 0, decomp_lbc_nrows_no_halo-1

                        ! global start point () is given by:
                        ! the global start point of this PE (full_lbc_start_pt)
                        ! + the points aready done from the previous level
                        !   (igposl)
                        ! + the points from rows we have already done at this
                        !   level (j * decomp_lbc_row_len)
                        igpos1 = full_lbc_start_pt + igposl                    &
                                 + j * decomp_lbc_row_len

                        ! if we are at pnorth, we also need to skip the upper
                        ! external halo
                        IF (g_at_extremity(pnorth,iproc)) THEN
                          igpos1 = igpos1 + decomp_lbc_row_len                 &
                                 * (decomp_lbc_nrows - decomp_lbc_nrows_no_halo)
                        END IF

                        ! if we are at peast, we also need to skip the leftmost
                        ! external halo
                        IF (g_at_extremity(peast,iproc)) THEN
                          igpos1 = igpos1 + (decomp_lbc_row_len -              &
                                                     decomp_lbc_row_len_no_halo)
                        END IF

                        DO i = 1, decomp_lbc_row_len_no_halo
                          igpos =  i + igpos1
                          send_buf(ipos) = read_field(igpos)
                          ipos = ipos + 1
                        END DO
                      END DO

                    END DO

                  END IF ! At extremity
                END DO ! iside loop

                address(iproc) = ipos

              END DO ! iproc loop

            CASE (ppx_atm_compressed )
              ! ppx_atm_compressed, which will be unpacked prior to comms.

              IF (atmos_number_of_landpts > 0) THEN

                ALLOCATE(unpacked_field(glsize(1,fld_type)*glsize(2,fld_type)))

                CALL expand_from_mask(unpacked_field,                          &
                                      read_field,                              &
                                      atmos_landmask,                          &
                                      glsize(1,fld_type)*glsize(2,fld_type),   &
                                      maskpoints)

                DEALLOCATE(read_field)
                CALL MOVE_ALLOC(unpacked_field,read_field)


                DO iproc = 0, nproc-1

                  IF (atmos_number_of_landpts_proc(iproc) > 0) THEN

                    ipos = address(iproc) ! start of array

                    DO j = 0, g_blsize(2, fld_type, iproc)-1
                      igpos1 = ((g_datastart_f(2,fld_type,iproc) + j - 1) *    &
                                glsize(1,fld_type)) +                          &
                                (g_datastart_f(1,fld_type,iproc) - 1)

                      DO i = 1, g_blsize(1, fld_type, iproc)
                        igpos =  i + igpos1
                        send_buf(ipos) = read_field(igpos)
                        ipos = ipos + 1
                      END DO

                    END DO

                    address(iproc) = ipos

                  END IF

                END DO

              END IF

            CASE DEFAULT ! all except zonal, meridional, ppx_atm_compressed,
                         ! and LBC fields.

              DO iproc = 0, nproc-1
                ipos = address(iproc) ! start of array
                DO j = 0, g_blsize(2, fld_type, iproc)-1
                  igpos1 = ((g_datastart_f(2,fld_type,iproc) + j - 1) *        &
                            glsize(1,fld_type)) +                              &
                            (g_datastart_f(1,fld_type,iproc) - 1)

                  DO i = 1, g_blsize(1, fld_type, iproc)
                    igpos =  i + igpos1
                    send_buf(ipos) = read_field(igpos)
                    ipos = ipos + 1
                  END DO
                END DO

                address(iproc) = ipos
              END DO

          END SELECT

        END IF

        DEALLOCATE(read_field)

      ELSE ! read multi

        ! Rank 0 reads in serially, then distributes fields to all other ranks

        CALL get_grid(fixhd(5),                                                &
                      field_model,                                             &
                      field_sect,                                              &
                      field_item,                                              &
                      lookup(lbhem,k),                                         &
                      lookup(lbpack,k),                                        &
                      lookup(lbnpt,k),                                         &
                      grid_type,                                               &
                      halo_type,                                               &
                      fld_type,                                                &
                      icode)

        IF ( icode /= 0 ) THEN
          cmessage='READFLDS : get_grid call failed'
          GO TO 9999
        END IF

        ! Check if LBC, and if so, if its type is valid
        IF ((lookup(lbhem,k)  >=  99) .AND.                                    &
            (lookup(lbhem,k)  <   1000)) THEN ! This is a LBC field

          IF (field_model  /=  atmos_im) THEN
            icode=2
            WRITE(umMessage,'(A,I0)') 'READFLDS: Cannot process LBC for ' //   &
                                      'model type ', field_model
            CALL umPrint(umMessage,src='readflds')
            cmessage='READFLDS : Cannot read LBCS for this model type'
            GO TO 9999
          END IF

        END IF

        ! Set the number of levels in the field

        IF (  grid_type == ppx_atm_lbc_theta                                   &
            .OR. grid_type == ppx_atm_lbc_u                                    &
            .OR.  grid_type == ppx_atm_lbc_v) THEN

          CALL get_lbc_levels(field_model,                                     &
                              field_sect,                                      &
                              field_item,                                      &
                              num_levels,                                      &
                              icode)

        ELSE

          num_levels = 1

        END IF

        CALL readflds_local_field_length(field_length,                         &
                                         grid_type,                            &
                                         fld_type,                             &
                                         halo_type,                            &
                                         lookup(lbproc,k),                     &
                                         num_levels,                           &
                                         k,                                    &
                                         mype,                                 &
                                         .TRUE.,                               &
                                         typemap,                              &
                                         north,                                &
                                         south,                                &
                                         east,                                 &
                                         west,                                 &
                                         gridpoint,                            &
                                         st_no_data_test_fix=.TRUE.)

        local_len=field_length

        IF (PRESENT(north)) THEN
          north_loc = north(k)
        ELSE
          north_loc = glsize(2,fld_type)
        END IF

        IF (PRESENT(south)) THEN
          south_loc = south(k)
        ELSE
          south_loc = 1
        END IF

        IF (PRESENT(east)) THEN
          east_loc = east(k)
        ELSE
          east_loc = glsize(1,fld_type)
        END IF

        IF (PRESENT(west)) THEN
          west_loc = west(k)
        ELSE
          west_loc = 1
        END IF

        IF (PRESENT(gridpoint)) THEN
          gridpoint_loc = gridpoint(k)
        ELSE
          gridpoint_loc = unset
        END IF

        IF (PRESENT(proccode)) THEN
          proccode_loc=proccode(k)
        ELSE
          proccode_loc=unset
        END IF

        ! Read the field from disk and distribute it over the processors

        CALL read_multi(nftin,                                                 &
                        d1(d1_start_read:),                                    &
                        data_read_size,                                        &
                        data_full_size,                                        &
                        len_io,                                                &
                        local_len,                                             &
                        lookup(1:64,k),                                        &
                        fixhd(12),                                             &
                        grid_type,                                             &
                        halo_type,                                             &
                        field_type,                                            &
                        proccode_loc,                                          &
                        field_length,                                          &
                        north_loc,                                             &
                        south_loc,                                             &
                        east_loc,                                              &
                        west_loc,                                              &
                        gridpoint_loc,                                         &
                        field_model,                                           &
                        field_sect,                                            &
                        field_item,                                            &
                        num_levels,                                            &
                        icode,                                                 &
                        cmessage,                                              &
                        expand)

#endif

      END IF

      IF (PRESENT(mpp_lookup) .AND. .NOT.alltoall) THEN
        mpp_lookup(p_lblrec,k)=local_len
      END IF

      !   --------------------------------
      !   Check that data been read in OK?
      !   --------------------------------

      IF ((icode /= 0) .OR. (len_io /= data_read_size)) THEN

        WRITE(umMessage,'(A,I0,A,I0,A,I0)') 'readflds - Error while ' //       &
                                            'attempting to read field ', k,    &
                                            ' of ', number_of_fields,          &
                                            ' from unit ', nftin
        CALL umPrint(umMessage,src='readflds')

        IF (icode /= 0) THEN
          WRITE(umMessage,'(A,I0,A,A)') 'Return code from READ_MULTI was ',    &
                                        icode, ' and error message was ',      &
                                        TRIM(cmessage)
          CALL umPrint(umMessage,src='readflds')
        END IF

        CALL umPrint('Field Information: ',src='readflds')

        WRITE(umMessage,'(A,I0)') 'Number of words requested : ', data_read_size
        CALL umPrint(umMessage,src='readflds')

        WRITE(umMessage,'(A,I0)') 'Number of words returned : ',len_io
        CALL umPrint(umMessage,src='readflds')

        WRITE(umMessage,'(A,I0)') 'MODEL ID ',field_model
        CALL umPrint(umMessage,src='readflds')

        WRITE(umMessage,'(A,I0)') 'SECTION ',field_sect
        CALL umPrint(umMessage,src='readflds')

        WRITE(umMessage,'(A,I0)') 'ITEM ',field_item
        CALL umPrint(umMessage,src='readflds')

        WRITE(umMessage,'(A,I0)') 'Disk address : ',field_start
        CALL umPrint(umMessage,src='readflds')

        WRITE(umMessage,'(A,I0)') 'D1 address : ',d1_start_read
        CALL umPrint(umMessage,src='readflds')

        IF (fixhd(5) <  6 .OR. fixhd(5) >  10) THEN

          CALL pr_look(lookup,k)

        END IF

        cmessage = 'READFLDS: I/O error'
        icode=6
        CALL ereport(routinename, icode, cmessage)

      END IF ! If an error was detected reading the field

    ELSE ! Error in LOOKUP(DATA_TYPE,K)

      IF (( fixhd(5)  <   6) .OR.                                              &
          ( fixhd(5)  >   8)) THEN ! Not AC, Var or Cx

        CALL pr_look(lookup,k)

      END IF

      WRITE(umMessage,'(A,I0,A,I0)') 'readflds : failure for field ', k,       &
                                     ' of ', number_of_fields
      CALL umPrint(umMessage,src='readflds')
      WRITE(umMessage,'(A,I0)') 'LOOKUP(DATA_TYPE,K)= ',lookup(data_type,k)
      CALL umPrint(umMessage,src='readflds')
      icode=7
      cmessage='Invalid data type ( LOOKUP(DATA_TYPE,K) )'

      CALL ereport(routinename, icode, cmessage)

    END IF

#if defined(UTILIO) && !defined(VOMEXT)
    IF (printstatus >= prstatus_diag) THEN
      ! Write out information about the field just read in
      IF (mype  ==  0) THEN
        IF (fixhd(5) <  6 .OR. fixhd(5) >  8) THEN

          ! Print out header and summary of data field

          CALL pr_look(lookup,k)

          IF ( (pack_code  ==  PC_WGDOS_packing) .AND. expand /= 1)  THEN
            WRITE(umMessage,'(A,A)') 'WGDOS packing not supported .',          &
                                     'Field summary omitted.'
            CALL umPrint(umMessage,src='readflds')
          ELSE IF (pack_code  ==  PC_GRIB_Packing) THEN
            WRITE(umMessage,'(A,A)') 'GRIB compression not supported .',       &
                                     'Field summary omitted.'
            CALL umPrint(umMessage,src='readflds')
          ELSE IF ((pack_code == PC_RunLength_Packing) .AND. (expand /= 1)) THEN
            WRITE(umMessage,'(A,A)') 'RLE packing not supported .',            &
                                     'Field summary omitted.'
            CALL umPrint(umMessage,src='readflds')
          ELSE IF ((pack_code == PC_Cray32_packing) .AND. (expand /= 1)) THEN
            WRITE(umMessage,'(A,A)') 'Cray32 packing not supported .',         &
                                     'Field summary omitted.'
            CALL umPrint(umMessage,src='readflds')
          ELSE IF (fixhd(5)  ==  5) THEN
            WRITE(umMessage,'(A,A)') 'Boundary dataset .',                     &
                                     'Field summary omitted.'
            CALL umPrint(umMessage,src='readflds')
          ELSE IF (lookup(data_type,k)  <   0) THEN
            WRITE(umMessage,'(A,A)') 'Time series field. ',                    &
                                     'Field summary omitted.'
            CALL umPrint(umMessage,src='readflds')
          ELSE
            ! Write out summary of field

            IF ( (fixhd(2)  ==  1) .AND. (lookup(item_code,k) == 30) ) THEN
              ! Land-sea mask

              ! DEPENDS ON: pr_lfld
              CALL pr_lfld(lookup,lookup,SIZE(lookup,DIM=1),d1(d1_off+1),k)

            ELSE IF (lookup(data_type,k) == 1) THEN  !  Real
              ! DEPENDS ON: pr_rfld
              CALL pr_rfld(lookup,lookup,d1(d1_off+1),k)

            ELSE IF (lookup(data_type,k) == 2) THEN  !  Integer
              ! DEPENDS ON: pr_ifld
              CALL pr_ifld(lookup,lookup,d1(d1_off+1),k)

            ELSE IF (lookup(data_type,k) == 3) THEN  !  Logical
              ! DEPENDS ON: pr_lfld
              CALL pr_lfld(lookup,lookup,SIZE(lookup,DIM=1),d1(d1_off+1),k)

            END IF ! type of field

          END IF ! if this field can have a summary

        END IF ! if this file is suitable for summaries

      END IF ! IF (mype  ==  0)
    END IF ! PrintStatus >= PrStatus_Diag
#endif

  END IF ! If there was data in the field

  ! Move to start address of next field to read in,
  ! i.e. the start of the current field + plus the
  ! length of the field we have just read in

  IF (.NOT.alltoall) THEN
    IF (local_len == unset) THEN
      d1_off=d1_off+lookup(lblrec,k)
    ELSE
      d1_off=d1_off+local_len
    END IF
  END IF

END DO ! K : loop over fields to read in

IF (alltoall) THEN

  ! Do checks for remainder of fields
  DO k = start_field(mype+1), first_field+number_of_fields-1

    ! Reset local_len
    local_len = unset

    IF (lookup(lblrec,k)  >   0) THEN ! If there's data in
                                      ! the field

      ! Check that DATA_TYPE is valid no: +/-1 to +/-3
      IF (( ABS(lookup(data_type,k))  >=  1) .AND.                             &
          ( ABS(lookup(data_type,k))  <=  3)) THEN

        ! Get some information about this field

        field_item=MOD(lookup(item_code,k),1000)
        field_sect=(lookup(item_code,k)-field_item)/1000
        field_model=lookup(model_code,k)

        !   ------------------------
        !   Perform sanity checks
        !   ------------------------

        d1_start_read = d1_off + 1

        ! Check address is within array bounds
        IF (d1_start_read > UBOUND(d1,dim=1)) THEN

          CALL get_grid(fixhd(5),                                              &
                        field_model,                                           &
                        field_sect,                                            &
                        field_item,                                            &
                        lookup(lbhem,k),                                       &
                        lookup(lbpack,k),                                      &
                        lookup(lbnpt,k),                                       &
                        grid_type,                                             &
                        halo_type,                                             &
                        fld_type,                                              &
                        icode)


          IF ( icode /= 0 ) THEN
            cmessage='READFLDS : get_grid call failed'
            GO TO 9999
          END IF

          num_levels=1

          ! for LBCs we need some special logic...
          IF ( (lookup(lbhem,k)  >=  99) .AND.                                 &
               (lookup(lbhem,k)  <   1000)) THEN ! This is a LBC field

            IF (field_model  ==  atmos_im) THEN
              IF (lookup(lbhem,k)  /=  99) THEN
                ! New style LBCs with different field types
                ! See UMDP C04

                CALL get_lbc_levels(field_model,                               &
                                    field_sect,                                &
                                    field_item,                                &
                                    num_levels,                                &
                                    icode)

              END IF
            ELSE
              icode=2
              WRITE(umMessage,'(A,I0)') 'READFLDS: Cannot process LBC for ' // &
                                        'model type ', field_model
              CALL umPrint(umMessage,src='readflds')
              cmessage='READFLDS : Cannot read LBCS for this model type'
              GO TO 9999
            END IF

          END IF

          CALL readflds_local_field_length(local_len,                          &
                                           grid_type,                          &
                                           fld_type,                           &
                                           halo_type,                          &
                                           lookup(lbproc,k),                   &
                                           num_levels,                         &
                                           k,                                  &
                                           mype,                               &
                                           .TRUE.,                             &
                                           typemap,                            &
                                           north,                              &
                                           south,                              &
                                           east,                               &
                                           west,                               &
                                           gridpoint,                          &
                                           st_no_data_test_fix=.TRUE.)

          IF ( local_len > 0) THEN
            icode=5
            WRITE(umMessage,'(4(A,I0))') "Readflds Address error for field ",  &
                                         k, ":"  //  newline  //               &
                                         "Calculated (address)"//              &
                                         "/(field_length): (", d1_start_read,  &
                                         ')/(', local_len, ')' // newline //   &
                                         "Maximum D1 address: ",               &
                                         UBOUND(d1,dim=1)
            CALL ereport(routinename, icode, umMessage)
          ELSE
            ! A zero sized field at the end of the d1-array
            ! set the d1_start_address to equal the length
            ! of the d1 array to prevent out of bounds access
            d1_start_read = UBOUND(d1,dim=1) - 1

          END IF

          local_len = unset

        END IF

        ! Check field is included in in_s array (not applicable to small execs)
        IF (.NOT. isSmallExec()) THEN
          IF (in_s(1,field_model,field_sect,field_item) == 0) THEN
            WRITE(umMessage,'(A)') 'ERROR: Updated field not used in this ' // &
                                   'model run.'
            CALL umPrint(umMessage,src='readflds',level=PrMin)
            WRITE(umMessage,'(A)') 'Check your ancillary files are correct ' //&
                                   'as it may'
            CALL umPrint(umMessage,src='readflds',level=PrMin)
            WRITE(umMessage,'(A)') 'be that the updating has exceeded the ' // &
                                   'end of an'
            CALL umPrint(umMessage,src='readflds',level=PrMin)
            WRITE(umMessage,'(A)') 'ancillary file and it has gone onto the' //&
                                   ' next'
            CALL umPrint(umMessage,src='readflds',level=PrMin)
            WRITE(umMessage,'(A)') 'field (which is this unacceptable one).' //&
                                   ' Look'
            CALL umPrint(umMessage,src='readflds',level=PrMin)
            WRITE(umMessage,'(A)') 'at the ancillary file associated with ' // &
                                   'the stash'
            CALL umPrint(umMessage,src='readflds',level=PrMin)
            WRITE(umMessage,'(A)') 'code before the following one:'
            CALL umPrint(umMessage,src='readflds',level=PrMin)
            WRITE(umMessage,'(A,I0)') 'Model ID : ',field_model
            CALL umPrint(umMessage,src='readflds',level=PrMin)
            WRITE(umMessage,'(A,I0)') 'Section  : ',field_sect
            CALL umPrint(umMessage,src='readflds',level=PrMin)
            WRITE(umMessage,'(A,I0)') 'Item     : ',field_item
            CALL umPrint(umMessage,src='readflds',level=PrMin)
            icode=101
            cmessage='Updated field not used in this model run. See output.'
            GO TO 9999
          END IF
        END IF

        !   ------------------------
        !   calculate offset for next field
        !   ------------------------

#if !defined(UTILIO)


        CALL get_grid(fixhd(5),                                                &
                      field_model,                                             &
                      field_sect,                                              &
                      field_item,                                              &
                      lookup(lbhem,k),                                         &
                      lookup(lbpack,k),                                        &
                      lookup(lbnpt,k),                                         &
                      grid_type,                                               &
                      halo_type,                                               &
                      fld_type,                                                &
                      icode)

        IF ( icode /= 0 ) THEN
          cmessage='READFLDS : get_grid call failed'
          GO TO 9999
        END IF

        ! Set the number of levels in the field

        IF (  grid_type == ppx_atm_lbc_theta                                   &
            .OR. grid_type == ppx_atm_lbc_u                                    &
            .OR.  grid_type == ppx_atm_lbc_v) THEN

          CALL get_lbc_levels(field_model,                                     &
                              field_sect,                                      &
                              field_item,                                      &
                              num_levels,                                      &
                              icode)

        ELSE

          num_levels = 1

        END IF

        CALL readflds_local_field_length(field_length,                         &
                                         grid_type,                            &
                                         fld_type,                             &
                                         halo_type,                            &
                                         lookup(lbproc,k),                     &
                                         num_levels,                           &
                                         k,                                    &
                                         mype,                                 &
                                         .TRUE.,                               &
                                         typemap,                              &
                                         north,                                &
                                         south,                                &
                                         east,                                 &
                                         west,                                 &
                                         gridpoint,                            &
                                         st_no_data_test_fix=.TRUE. )

        d1_off=d1_off+field_length

#endif

      ELSE ! Error in LOOKUP(DATA_TYPE,K)

        IF (( fixhd(5)  <   6) .OR.                                            &
            ( fixhd(5)  >   8)) THEN ! Not AC, Var or Cx

          CALL pr_look(lookup,k)

        END IF

        WRITE(umMessage,'(A,I0,A,I0)') 'readflds : failure for field ', k,     &
                                       ' of ', number_of_fields
        CALL umPrint(umMessage,src='readflds')
        WRITE(umMessage,'(A,I0)') 'LOOKUP(DATA_TYPE,K)= ',lookup(data_type,k)
        CALL umPrint(umMessage,src='readflds')
        icode=7
        cmessage='Invalid data type ( LOOKUP(DATA_TYPE,K) )'

        CALL ereport(routinename, icode, cmessage)

      END IF

    END IF ! If there was data in the field

  END DO ! K : loop over fields to read in

#if !defined(UTILIO)
  ! Do the data swap between all ranks. Replaces scatterv
  CALL gc_get_communicator(my_comm, info)
  CALL MPL_Alltoallv(send_buf, send_size, send_offset, mpl_real8, recv_buf,    &
                     recv_size, recv_offset, mpl_real8, my_comm, info)

  ! Unpack the data into the d1 array. Do manual halo placement
  igpos = 1
  ipos = 1

  DO K = 1,number_of_fields

    ! The current offset in the d1 is the igpos array
    IF (PRESENT(mpp_lookup)) THEN
      mpp_lookup(p_naddr,k) = igpos
    END IF

    IF (lookup(lblrec,k) > 0) THEN

      field_item=MOD(lookup(item_code,k),1000)
      field_sect=(lookup(item_code,k)-field_item)/1000
      field_model=lookup(model_code,k)

      CALL get_grid(fixhd(5),                                                  &
                    field_model,                                               &
                    field_sect,                                                &
                    field_item,                                                &
                    lookup(lbhem,k),                                           &
                    lookup(lbpack,k),                                          &
                    lookup(lbnpt,k),                                           &
                    grid_type,                                                 &
                    halo_type,                                                 &
                    fld_type,                                                  &
                    icode)

      IF (PRESENT(north)) THEN
        north_loc = north(k)
      ELSE
        north_loc = glsize(2,fld_type)
      END IF

      IF (PRESENT(south)) THEN
        south_loc = south(k)
      ELSE
        south_loc = 1
      END IF

      IF (PRESENT(east)) THEN
        east_loc = east(k)
      ELSE
        east_loc = glsize(1,fld_type)
      END IF

      IF (PRESENT(west)) THEN
        west_loc = west(k)
      ELSE
        west_loc = 1
      END IF

      IF (PRESENT(gridpoint)) THEN
        gridpoint_loc = gridpoint(k)
      ELSE
        gridpoint_loc = unset
      END IF

      IF (PRESENT(typemap)) THEN
        typemap_loc = typemap(k)
      ELSE
        typemap_loc = unset
      END IF

      subdomain = field_issubdomain(fld_type,                                  &
                                    lookup(lbproc,k),                          &
                                    typemap_loc,                               &
                                    north_loc,                                 &
                                    south_loc,                                 &
                                    east_loc,                                  &
                                    west_loc,                                  &
                                    gridpoint_loc)


      ! pack the field into the send buffer - this does not include halos.

      IF (subdomain) THEN

        ! We should not be able to get here yet.
        icode=1
        cmessage = 'READFLDS: ' //                                             &
                   'Reached a subdomained field with all-to-all enabled.'
        CALL ereport(routinename, icode, cmessage)

      END IF

      SELECT CASE( grid_type )

        ! If data is zonal then copy across the E-W dimension
        ! (we are only supporting the no-halo case for zonal fields)
        CASE (ppx_atm_tzonal,ppx_atm_uzonal)

          IF (subdomain) THEN

            ! We should not be able to get here yet.
            icode=1
            cmessage = 'READFLDS: ' //                                         &
                       'Reached a subdomained field with all-to-all enabled.'
            CALL ereport(routinename, icode, cmessage)

          ELSE

            DO j = 1,g_blsize(2,fld_type,mype)
              d1(igpos) = recv_buf(ipos)
              ipos = ipos + 1
              igpos = igpos + 1
            END DO

          END IF

        ! If data is zonal then copy across the N-S dimension
        ! (we are only supporting the no-halo case for meridional fields)
        CASE (ppx_atm_tmerid, ppx_atm_umerid)

          IF (subdomain) THEN

            ! We should not be able to get here yet.
            icode=1
            cmessage = 'READFLDS: ' //                                         &
                       'Reached a subdomained field with all-to-all enabled.'
            CALL ereport(routinename, icode, cmessage)

          ELSE

            DO j = 1,g_blsize(1,fld_type,mype)
              d1(igpos) = recv_buf(ipos)
              ipos = ipos + 1
              igpos = igpos + 1
            END DO

          END IF

        ! If the data is a compressed field, it has been scattered uncompressed.
        ! Therefore, re-compress it.
        CASE (ppx_atm_compressed)

          IF (subdomain) THEN

            ! We should not be able to get here yet.
            icode=1
            cmessage = 'READFLDS: ' //                                         &
                       'Reached a subdomained field with all-to-all enabled.'
            CALL ereport(routinename, icode, cmessage)

          ELSE

            IF (atmos_number_of_landpts_proc(mype) > 0) THEN
              CALL compress_to_mask(recv_buf(ipos:), d1(igpos:),               &
                                    atmos_landmask_local,                      &
                                    g_lasize(1,fld_type,halo_type,mype)*       &
                                    g_lasize(2,fld_type,halo_type,mype),       &
                                    maskpoints)

              igpos = igpos + maskpoints
              ipos = ipos + g_lasize(1,fld_type,halo_type,mype) *              &
                            g_lasize(2,fld_type,halo_type,mype)
            END IF

          END IF

        CASE ( ppx_atm_rim, ppx_atm_lbc_theta, ppx_atm_lbc_u, ppx_atm_lbc_v,   &
               ppx_atm_lbc_orog )

          ! LBC fields cannot be diagnostics, and therefore cannot be
          ! subdomained.

          IF (grid_type  ==  ppx_atm_lbc_orog) THEN
            rim_type=rima_type_orog
          ELSE
            rim_type=rima_type_norm
          END IF

          DO iside=1,4
            IF (g_at_extremity(iside,iproc)) THEN

              igposg = igpos

              ! This processor is at edge type iside and so needs LBC data
              ! In order to copy data from full_lbc array to data_lbc ( that
              ! contains data to be scattered ) the following variables need to
              ! be calculated

              ! full_lbc_row_len : East-West dimension of the full lbc side
              ! full_lbc_nrows : North-South dimension of the full lbc side
              ! decomp_lbc_row_len : East-West dimension of decomp lbc side
              ! decomp_lbc_nrows : North-South dimension of decomp lbc side
              ! full_lbc_start_pt : First point of the decomposed lbc side
              !                     inside the 1d full_lbc array
              ! decomp_lbc_start_pt : First point of the decomposed lbc side
              !                       inside the 1d decomposed lbc array
              ! data_lbc_start : First point inside the data lbc array
              !                  ( corresponds to the first point of the
              !                    distributed lbc_comp )


              CALL lbc_calc_size(iside,                                        &
                                 full_lbc_row_len,                             &
                                 full_lbc_nrows,                               &
                                 decomp_lbc_row_len,                           &
                                 decomp_lbc_nrows,                             &
                                 full_lbc_start_pt,                            &
                                 decomp_lbc_start_pt,                          &
                                 fld_type,                                     &
                                 halo_type,                                    &
                                 rim_type,                                     &
                                 iproc)

              CALL lbc_calc_size(iside,                                        &
                                 full_lbc_row_len_no_halo,                     &
                                 full_lbc_nrows_no_halo,                       &
                                 decomp_lbc_row_len_no_halo,                   &
                                 decomp_lbc_nrows_no_halo,                     &
                                 full_lbc_start_pt_no_halo,                    &
                                 decomp_lbc_start_pt_no_halo,                  &
                                 fld_type,                                     &
                                 halo_type_no_halo,                            &
                                 rim_type,                                     &
                                 iproc)

              DO lbck=0,num_levels-1

                ! Skip the already done levels
                igposl = lbck * decomp_lbc_nrows * decomp_lbc_row_len

                DO j = 0, decomp_lbc_nrows_no_halo-1

                  ! global start point () is given by:
                  ! the global start point of this PE (full_lbc_start_pt)
                  ! + the points aready done from the previous level (igposl)
                  ! + the points from rows we have already done at this level
                  !   (j * decomp_lbc_row_len)
                  igpos1 = full_lbc_start_pt + igposl + j * decomp_lbc_row_len

                  ! if we are at pnorth, we also need to skip the upper
                  ! external halo
                  IF (g_at_extremity(pnorth,iproc)) THEN
                    igpos1 = igpos1 + decomp_lbc_row_len *                     &
                                   (decomp_lbc_nrows - decomp_lbc_nrows_no_halo)
                  END IF

                  ! if we are at peast, we also need to skip the leftmost
                  ! external halo
                  IF (g_at_extremity(peast,iproc)) THEN
                    igpos1 = igpos1 + (decomp_lbc_row_len -                    &
                                                     decomp_lbc_row_len_no_halo)
                  END IF

                  DO i = 1, decomp_lbc_row_len_no_halo
                    igpos = igposg + i + igpos1
                    d1(igpos) = recv_buf(ipos)
                    send_buf(ipos) = read_field(igpos)
                    ipos = ipos + 1
                  END DO

                END DO

              END DO

              igpos = igposg + decomp_lbc_nrows * decomp_lbc_row_len

            END IF ! At extremity
          END DO ! iside loop

        ! All other fields, copy from receive buffer to d1 skipping halos
        CASE DEFAULT

          IF (subdomain) THEN

            ! We should not be able to get here yet.
            icode=1
            cmessage = 'READFLDS: ' //                                         &
                       'Reached a subdomained field with all-to-all enabled.'
            CALL ereport(routinename, icode, cmessage)

          ELSE

            ! The first halo rows
            igpos = igpos + (halosize(2,halo_type) *                           &
                           g_lasize(1,fld_type,halo_type,mype))
            DO j = 1,g_blsize(2,fld_type,mype)
              ! The first halo columns
              igpos = igpos + halosize(1,halo_type)
              ! Copy the raw data
              DO i=1,g_blsize(1,fld_type,mype)
                d1(igpos) = recv_buf(ipos)
                ipos = ipos + 1
                igpos = igpos + 1
              END DO
              ! The final halo columns
              igpos = igpos + halosize(1,halo_type)
            END DO
            ! The final halo rows
            igpos = igpos + (halosize(2,halo_type) *                           &
                           g_lasize(1,fld_type,halo_type,mype))

          END IF

      END SELECT

    END IF

  END DO

  ! Deallocate the temporary buffers created
  DEALLOCATE(recv_buf)
  DEALLOCATE(send_buf)

  IF (PRESENT(mpp_lookup)) THEN

    DO k = 1,number_of_fields - 1

      mpp_lookup(p_lblrec,k) = mpp_lookup(p_naddr,k+1) - mpp_lookup(p_naddr,k)

    END DO

    k = number_of_fields

    mpp_lookup(p_lblrec,k) = igpos - mpp_lookup(p_naddr,k)

  END IF

#endif
END IF

9999 CONTINUE

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE readflds_core

!------------------------------------------------------------------------------!

SUBROUTINE get_grid(fixhd_5, field_model, field_sect, field_item, lbhem,       &
                    lbpack, lbnpt, grid_type, halo_type, fld_type, icode)

! get_grid() calculates the grid type, halo type and field type of a field from
! its STASH code and lookup headers.

USE field_types, ONLY:                                                         &
    fld_type_p,                                                                &
    fld_type_u,                                                                &
    fld_type_v

USE cppxref_mod, ONLY:                                                         &
    ppx_grid_type,                                                             &
    ppx_halo_type,                                                             &
    ppx_atm_rim,                                                               &
    ppx_atm_compressed,                                                        &
    ppx_atm_tall,                                                              &
    ppx_atm_cuall,                                                             &
    ppx_atm_cvall,                                                             &
    ppx_atm_ozone,                                                             &
    ppx_atm_tzonal

USE umPrintMgr, ONLY:                                                          &
    umPrint,                                                                   &
    umMessage

USE um_parparams, ONLY:                                                        &
    halo_type_no_halo

USE errormessagelength_mod, ONLY:                                              &
    errormessagelength

USE submodel_mod, ONLY:                                                        &
    atmos_im

USE ppxlook_mod, ONLY:                                                         &
    exppxi

USE packing_codes_mod, ONLY:                                                   &
    PC_No_CompressType

IMPLICIT NONE

! Arguments

INTEGER, INTENT(IN) ::                                                         &
  fixhd_5,                                                                     &
  field_model,                                                                 &
  field_sect,                                                                  &
  field_item,                                                                  &
  lbhem,                                                                       &
  lbpack,                                                                      &
  lbnpt

INTEGER, INTENT(OUT) ::                                                        &
  grid_type,                                                                   &
  halo_type,                                                                   &
  fld_type,                                                                    &
  icode

! Local Variables

CHARACTER(LEN=errormessagelength) ::                                           &
  cmessage

! Functions

INTEGER, EXTERNAL ::                                                           &
  get_fld_type

icode = 0

IF (fixhd_5  >=  6 .AND. fixhd_5  <=  9) THEN
  ! Set grid_type and halo_type for ACobs and VARobs, Cx and CovStats
  ! (FIXHD(5)=6, 7, 8 and 9, respectively).
  grid_type = ppx_atm_tall
  halo_type = halo_type_no_halo
ELSE
  grid_type=exppxi(field_model, field_sect, field_item, ppx_grid_type, icode,  &
                   cmessage)
  halo_type=exppxi(field_model, field_sect, field_item, ppx_halo_type, icode,  &
                   cmessage)
  IF (icode  /=  0) THEN
    WRITE(umMessage,'(A)') 'READFLDS - EXPPXI failed to get PPXREF ' //        &
                           'information for field: '
    CALL umPrint(umMessage,src='readflds')
    WRITE(umMessage,'(A,I0)') 'Model ID : ',field_model
    CALL umPrint(umMessage,src='readflds')
    WRITE(umMessage,'(A,I0)') 'Section  : ',field_sect
    CALL umPrint(umMessage,src='readflds')
    WRITE(umMessage,'(A,I0)') 'Item     : ',field_item
    CALL umPrint(umMessage,src='readflds')
    WRITE(umMessage,'(A,I0)') 'Error code was ',icode
    CALL umPrint(umMessage,src='readflds')
    WRITE(umMessage,'(A,A)') 'Error message was ',cmessage
    CALL umPrint(umMessage,src='readflds')
    icode=1
  END IF
END IF

IF (icode  ==  0) THEN
  ! DEPENDS ON: get_fld_type
  fld_type=get_fld_type(grid_type) ! field type P,U or V

  ! for LBCs & Ancillaries we need some special logic...
  IF ( (lbhem  >=  99) .AND.                                                   &
       (lbhem  <   1000)) THEN ! This is a LBC field

    IF (field_model  ==  atmos_im) THEN
      IF (lbhem  ==  99) THEN ! Old style LBCs
        grid_type=ppx_atm_rim
      END IF
    ELSE
      icode=2
      WRITE(umMessage,'(A,I0)') 'READFLDS: Cannot process LBC for model type ',&
                                field_model
      CALL umPrint(umMessage,src='readflds')
    END IF

  ELSE IF (fixhd_5 == 4) THEN ! This is an ancillary File

    IF ( (MOD( lbpack/10, 10 ) == PC_No_CompressType )                         &
          .AND. grid_type == ppx_atm_compressed ) THEN

      ! Compressed in stashmaster but uncompressed in header
      SELECT CASE( fld_type )
        CASE ( fld_type_p )
          grid_type = ppx_atm_tall

        CASE ( fld_type_u )
          grid_type = ppx_atm_cuall

        CASE ( fld_type_v )
          grid_type = ppx_atm_cvall
      END SELECT

    END IF

  END IF

  ! For atmosphere zonal ozone fields - set to zonal grid type
  IF (( grid_type  ==  ppx_atm_ozone) .AND.  ( lbnpt  ==  1)) THEN
    grid_type=ppx_atm_tzonal
  END IF

END IF

END SUBROUTINE get_grid

!------------------------------------------------------------------------------!

SUBROUTINE get_lbc_levels(field_model,field_sect,field_item,lbc_levs,icode)

! LBCs can contain multiple levels in each field. This routine calculates the
! number of levels (lbc_levs) contained in each field for a given STASH code

USE errormessagelength_mod, ONLY:                                              &
    errormessagelength

USE ppxlook_mod, ONLY:                                                         &
    exppxi

USE cppxref_mod, ONLY:                                                         &
    ppx_lb_code, ppx_lt_code

USE nlsizes_namelist_mod, ONLY:                                                &
    model_levels

USE ereport_mod, ONLY:                                                         &
    ereport

IMPLICIT NONE

! Arguments

INTEGER, INTENT(IN) ::                                                         &
  field_model,field_sect,field_item

INTEGER, INTENT(OUT) ::                                                        &
  lbc_levs, icode

! Local variables

INTEGER ::                                                                     &
  lt_code, lb_code

CHARACTER(LEN=*), PARAMETER ::                                                 &
  RoutineName='get_lbc_levels'

CHARACTER(LEN=errormessagelength) ::                                           &
  cmessage

lt_code = exppxi(field_model, field_sect, field_item, ppx_lt_code, icode,      &
                 cmessage)

SELECT CASE( lt_code )
  CASE ( -1, 1 )
    ! Top level is level one
    lbc_levs = 1

  CASE ( 2, 3 )
    ! Top level is level model_levels
    lbc_levs = model_levels

  CASE ( 19 )
    ! Top level is level model_levels + 1
    lbc_levs = model_levels+1

  CASE DEFAULT
    ! Unknown
    icode=7
    WRITE(cmessage,'(A,I0,A)') "Readflds (LBC) levels error for ppx_lt_code ", &
                                lt_code, " - this code is not handled"
    CALL ereport(routinename, icode, cmessage)

END SELECT

lb_code = exppxi(field_model, field_sect, field_item, ppx_lb_code, icode,      &
                 cmessage)

SELECT CASE( lb_code )
  CASE ( -1, 1 )
    ! Bottom level is level one (do nothing)

  CASE ( 38, 40 )
    ! Bottom level is level zero
    lbc_levs = lbc_levs + 1

  CASE DEFAULT
    ! Unknown
    icode=8
    WRITE(cmessage,'(A,I0,A)') "Readflds (LBC) levels error for ppx_lb_code ", &
                               lb_code, " - this code is not handled"
    CALL ereport(routinename, icode, cmessage)

END SELECT

END SUBROUTINE get_lbc_levels

!------------------------------------------------------------------------------!

SUBROUTINE readflds_local_field_length(local_field_len, grid_type, fld_type,   &
                                       halo_type, proccode, num_levels, lev,   &
                                       proc, halos, typemap, north, south,     &
                                       east, west, gridpoint,                  &
                                       st_no_data_test_fix)

! This routine is used to calculate the length (local_field_len) of the
! decomposed section of a global field, which is held on on rank "proc"

USE um_parvars

USE um_parparams, ONLY:                                                        &
    halo_type_no_halo

USE atm_land_sea_mask, ONLY:                                                   &
    atmos_number_of_landpts_proc

USE lbc_mod, ONLY: g_lenrima

USE cppxref_mod, ONLY:                                                         &
    ppx_atm_rim,                                                               &
    ppx_atm_compressed,                                                        &
    ppx_atm_tall,                                                              &
    ppx_atm_cuall,                                                             &
    ppx_atm_cvall,                                                             &
    ppx_atm_tzonal,                                                            &
    ppx_atm_ozone,                                                             &
    ppx_atm_uzonal,                                                            &
    ppx_atm_umerid,                                                            &
    ppx_atm_tmerid,                                                            &
    ppx_atm_lbc_u,                                                             &
    ppx_atm_lbc_v,                                                             &
    ppx_atm_lbc_theta,                                                         &
    ppx_atm_lbc_orog,                                                          &
    ppx_atm_tsea,                                                              &
    ppx_atm_river,                                                             &
    ppx_atm_uall,                                                              &
    ppx_atm_tland

USE d1_array_mod, ONLY: diagnostic

USE ereport_mod, ONLY: ereport

USE errormessagelength_mod, ONLY: errormessagelength

USE rimtypes

USE sterr_mod

IMPLICIT NONE

! SUBROUTINE arguments

INTEGER, INTENT(OUT) ::          &
  local_field_len                  ! OUT: the field length on this field on proc

INTEGER, INTENT(IN) ::           &
  num_levels,                    &
  proc,                          &
  grid_type,                     &
  fld_type,                      &
  halo_type,                     &
  proccode,                      & ! IN: defines the field's processing code
  lev

LOGICAL, INTENT(IN) ::           &
  halos

INTEGER, INTENT(IN), OPTIONAL :: &
  typemap(:),                    & ! IN: defines whether field is prognostic,
                                   !      diagnostic, or other
  north(:),                      & ! IN: defines the field's north code
  south(:),                      & ! IN: defines the field's south code
  east(:),                       & ! IN: defines the field's east code
  west(:),                       & ! IN: defines the field's west code
  gridpoint(:)                     ! IN: defines the   field's gridpoint code

LOGICAL, INTENT(IN), OPTIONAL :: &
  st_no_data_test_fix              ! IN: temporary fix logical (see below)

! Local Variables

LOGICAL ::                                                                     &
  subdomain, &
  loc_st_no_data_test_fix

INTEGER ::                                                                     &
  north_loc, south_loc, east_loc, west_loc, gridpoint_loc, typemap_loc,        &
  local_start_row, local_end_col, local_end_row, local_start_col

INTEGER ::                                                                     &
  mean_type, icode

CHARACTER(LEN=errormessagelength) ::                                           &
  cmessage

! Parameters
INTEGER, PARAMETER ::                                                          &
  unset=-1             ! unset values

CHARACTER(LEN=*), PARAMETER ::                                                 &
  RoutineName='readflds_local_field_length'

IF (PRESENT(st_no_data_test_fix)) THEN
  loc_st_no_data_test_fix = st_no_data_test_fix
ELSE
  loc_st_no_data_test_fix = .TRUE.
END IF

IF (PRESENT(north)) THEN
  north_loc = north(lev)
ELSE
  north_loc = glsize(2,fld_type)
END IF

IF (PRESENT(south)) THEN
  south_loc = south(lev)
ELSE
  south_loc = 1
END IF

IF (PRESENT(east)) THEN
  east_loc = east(lev)
ELSE
  east_loc = glsize(1,fld_type)
END IF

IF (PRESENT(west)) THEN
  west_loc = west(lev)
ELSE
  west_loc = 1
END IF

IF (PRESENT(gridpoint)) THEN
  gridpoint_loc = gridpoint(lev)
ELSE
  gridpoint_loc = unset
END IF

IF (PRESENT(typemap)) THEN
  typemap_loc = typemap(lev)
ELSE
  typemap_loc = unset
END IF

subdomain = field_issubdomain(fld_type,                                        &
                              proccode,                                        &
                              typemap_loc,                                     &
                              north_loc,                                       &
                              south_loc,                                       &
                              east_loc,                                        &
                              west_loc,                                        &
                              gridpoint_loc)

IF (subdomain) THEN

  ! DEPENDS ON: global_to_local_subdomain
  CALL global_to_local_subdomain(halos,                                        &
                                 halos,                                        &
                                 grid_type,                                    &
                                 halo_type,                                    &
                                 proc,                                         &
                                 south_loc,                                    &
                                 east_loc,                                     &
                                 north_loc,                                    &
                                 west_loc,                                     &
                                 local_start_row,                              &
                                 local_end_col,                                &
                                 local_end_row,                                &
                                 local_start_col)

  ! This section should really test against st_no_data. But errors elsewhere in
  ! the code mean this would break the d1_offset (because the D1 sizes are
  ! incorrectly allocated). We therefore protect against this with the
  ! st_no_data_test_fix logical. For more information see: ticket #2894

  ! Set st_no_data_test_fix = .TRUE. to maintain backwards compatible behaviour

  IF ((local_start_row==st_no_data .OR. local_start_col == st_no_data)         &
      .AND. .NOT.loc_st_no_data_test_fix )                                     &
  THEN
    local_field_len = 0
  ELSE

    IF (local_end_col < local_start_col) THEN
       local_end_col=glsize(1,fld_type)-(local_start_col-local_end_col+1)
       local_start_col=1
    END IF

    IF (local_end_row < local_start_row) THEN
       local_end_row=glsize(2,fld_type)-(local_start_row-local_end_row+1)
       local_start_row=1
    END IF

    local_field_len = (local_end_row - local_start_row + 1) *                  &
                      (local_end_col - local_start_col + 1)
  END IF

ELSE

  ! not a diagnostic field - must be full domain

  SELECT CASE( grid_type )
    CASE ( ppx_atm_compressed )
      ! This is a field compressed to land points
      local_field_len = atmos_number_of_landpts_proc(proc)

    CASE ( ppx_atm_tzonal, ppx_atm_uzonal )
      ! This is a zonal field
      IF (halos) THEN
        local_field_len = g_lasize(2,fld_type,halo_type, proc)
      ELSE
        local_field_len = g_blsize(2, fld_type, proc)
      END IF

    CASE ( ppx_atm_tmerid, ppx_atm_umerid )
      ! This is a meridional field
      IF (halos) THEN
        local_field_len = g_lasize(1,fld_type,halo_type,proc)
      ELSE
        local_field_len = g_blsize(1, fld_type, proc)
      END IF

    CASE ( ppx_atm_cuall, ppx_atm_cvall, ppx_atm_tall, ppx_atm_tsea,           &
           ppx_atm_tland, ppx_atm_ozone, ppx_atm_uall, ppx_atm_river )
      IF (halos) THEN
        local_field_len = g_lasize(1,fld_type,halo_type,proc) *                &
                          g_lasize(2,fld_type,halo_type,proc)
      ELSE
        local_field_len = g_blsize(1, fld_type, proc) * &
                                g_blsize(2, fld_type, proc)
      END IF

    CASE ( ppx_atm_rim, ppx_atm_lbc_theta, ppx_atm_lbc_u, ppx_atm_lbc_v )
      IF (halos) THEN
        local_field_len = num_levels * g_lenrima(fld_type, halo_type,          &
                                                 rima_type_norm, proc)
      ELSE
        local_field_len = num_levels * g_lenrima(fld_type, halo_type_no_halo,  &
                                                 rima_type_norm, proc)
      END IF

    CASE ( ppx_atm_lbc_orog )
      IF (halos) THEN
        local_field_len = g_lenrima(fld_type, halo_type, rima_type_orog, proc)
      ELSE
        local_field_len = g_lenrima(fld_type, halo_type_no_halo,               &
                                    rima_type_orog, proc)
      END IF

    CASE DEFAULT
      WRITE(cmessage,'(A,A,I0)') 'UM_READDUMP: Dump grid type not handled - ', &
                                 grid_type
      icode=5
      CALL ereport(routinename,icode,cmessage)

  END SELECT
END IF

END SUBROUTINE readflds_local_field_length

!------------------------------------------------------------------------------!

FUNCTION field_issubdomain(fld_type, proccode, typemap, north, south, east,    &
                           west, gridpoint)

USE d1_array_mod, ONLY: diagnostic
USE um_parvars, ONLY: glsize

IMPLICIT NONE

LOGICAL ::                                                                     &
  field_issubdomain

! SUBROUTINE arguments

INTEGER, INTENT(INOUT) ::          &
  north,                      & ! IN: defines the field's north code
  south,                      & ! IN: defines the field's south code
  east,                       & ! IN: defines the field's east code
  west                          ! IN: defines the field's west code

INTEGER, INTENT(IN) ::           &
  fld_type,                      &
  proccode                         ! IN: defines the field's processing code

INTEGER, INTENT(IN) :: &
  typemap,                    & ! IN: defines whether field is prognostic,
                                   !      diagnostic, or other

  gridpoint                     ! IN: defines the   field's gridpoint code

! Local Variables

INTEGER ::                                                                     &
  mean_type

CHARACTER(LEN=*), PARAMETER ::                                                 &
  RoutineName='field_issubdomain'

field_issubdomain = .FALSE.

! diagnostics may be subdomained
IF (typemap == diagnostic) THEN

  ! Actually a subdomain
  IF (north /= glsize(2,fld_type)                                              &
      .OR. south /= 1                                                          &
      .OR. east /= glsize(1,fld_type)                                          &
      .OR. west /= 1                                                           &
     ) THEN
    field_issubdomain = .TRUE.
  END IF

  ! Spatial means are effectively subdomains
  IF (IAND(proccode,64) == 64) THEN
    field_issubdomain = .TRUE.
    mean_type=gridpoint/10
    IF (mean_type  ==  2) THEN ! zonal mean
      west = 1
      east = west
    ELSE IF (mean_type  ==  3) THEN ! meridional mean
      south = 1
      north = south
    ELSE IF (mean_type  >=  4) THEN ! field/global mean
      west = 1
      south = 1
      east = west
      north = south
    END IF
  END IF

END IF

END FUNCTION field_issubdomain

!------------------------------------------------------------------------------!

END MODULE readflds_mod
