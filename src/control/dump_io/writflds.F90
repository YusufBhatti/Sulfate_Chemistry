! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Write out fields to a UM format file.

SUBROUTINE writflds ( UNIT,                                       &
                                    ! in
                      numfields,                                  &
                                    ! in
                      POSITION,                                   &
                                    ! in
                      lookup,                                     &
                                    ! in
                      len1lookup,                                 &
                                    ! in
                      d1,                                         &
                                    ! in
                      buflen,                                     &
                                    ! in
                      fixhd,                                      &
                                    ! in
                      icode,                                      &
                                    ! out
                      cmessage )    ! out

! Description:

!   Buffers out NumFields fields from D1 to a UM format file on unit
!   Unit, starting at field number Position.


! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dump I/O

! Declarations:

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE io
USE UM_ParVars
USE io_configuration_mod, ONLY: io_field_padding
USE lookup_addresses
USE cppxref_mod, ONLY: &
    ppx_atm_rim,       &
    ppx_halo_type,     &
    ppx_grid_type
USE submodel_mod, ONLY: atmos_im
USE ppxlook_mod, ONLY: exppxi
USE d1_array_mod, ONLY: d1_list_len, d1_object_type, d1_imodl,   &
                        d1_section, d1_item, d1_length,          &
                        d1_grid_type, d1_no_levels, d1_halo_type,&
                        prognostic
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE missing_data_mod, ONLY: imdi

USE errormessagelength_mod, ONLY: errormessagelength

USE Packing_Codes_Mod, ONLY: PC_Cray32_Packing

IMPLICIT NONE

! Subroutine arguments:

INTEGER,       INTENT(IN) ::                                      &
  len1lookup            ! Lookup dimension.

INTEGER,       INTENT(IN) ::                                      &

  UNIT,                                                           &
                        ! Unit number of UM file.
  numfields,                                                      &
                        ! Number of fields to write out.
  POSITION,                                                       &
                        ! Field number from which to begin I/O.
  lookup(len1lookup,*)  ! Lookup table for UM file.

REAL,          INTENT(IN) ::                                      &

  d1(*)                 ! Array containing fields to be written.

INTEGER,       INTENT(IN) ::                                      &

  buflen,                                                         &
                        ! Length of I/O buffer.
  fixhd(*)              ! Fixed-length header for UM file.

INTEGER,       INTENT(OUT) ::                                     &

  icode                 ! Return code. >0 => error.

CHARACTER(LEN=errormessagelength), INTENT(OUT) ::                 &

  cmessage              ! Error message.

! Local variables:

INTEGER :: i, j, k,                                               &
           lenio,                                                 &
           ipts,                                                  &
                           ! No. of values to be written to disk.
           wordaddress,                                           &
                           ! Address from which to begin I/O.
           l_ipts,                                                &
                           ! Record length during index search.
           um_sector_ipts,                                        &
                           ! No. of words to write, rounded up.
           ipts_write      ! No. of words written to disk.

REAL :: buf( ( (buflen+io_field_padding-1)/io_field_padding ) *       &
             io_field_padding )

!dir$ cache_align buf

INTEGER :: address,                                               &
                           ! Start address of D1 field.
           locallen,                                              &
           item,                                                  &
           section,                                               &
           model,                                                 &
           jcode,                                                 &
                           ! Return code.
           fld_type,                                              &
           grid_type,                                             &
           halo_type,                                             &
           fake_d1_addr(d1_list_len) ! Fake D1_addr record to
                                     ! be fed to write_multi.

! Integer functions called:

INTEGER  :: get_fld_type

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='WRITFLDS'

!- End of header ------------------------------------------------------

!----------------------------------------------------------------------
! [1]: Initialize.
!----------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
address = 1

icode    = 0
j        = 0
cmessage = ' '

!----------------------------------------------------------------------
! [2]: Write fields.
!----------------------------------------------------------------------

DO k = POSITION, POSITION + numfields - 1

  ! See whether data is stored in 32-bit format on disk.
  IF (MOD(lookup(lbpack, k), 10) == PC_Cray32_Packing) THEN
    ipts = ( lookup(lblrec, k) + 1 ) / 2
  ELSE
    ipts = ( lookup(lblrec, k) + 1 )
  END IF

  IF ( lookup(lbnrec, k) == 0    .OR.                             &
                                      ! Old format dumps.
       lookup(lbnrec, k) == imdi .OR.                             &
                                      ! \ Ocean ACOBS
       lookup(lbegin, k) == imdi .OR.                             &
                                      ! / files?
       ( lookup(lbnrec, k) == imdi                                &
                                      ! \ Prognostic lookup
         .AND. fixhd(12) < 301 )                                  &
                                      ! / in pre-vn3.2 dump.
     ) THEN

    wordaddress = 1
    DO i = 2, k
      IF (MOD(lookup(lbpack, i-1), 10) == PC_Cray32_Packing) THEN
        l_ipts = ( lookup(lblrec, i-1) + 1 ) / 2
      ELSE
        l_ipts =   lookup(lblrec, i-1)
      END IF
      wordaddress = wordaddress + l_ipts
    END DO
    wordaddress = wordaddress + fixhd(160) - 2
    um_sector_ipts = ipts

  ELSE ! PP type files and post vn4.3 dumps.

    wordaddress = lookup(lbegin, k)

    ! Use the stored rounded-up value:
    um_sector_ipts = lookup(lbnrec, k)

  END IF

  ipts_write = um_sector_ipts

  ! Position file pointer:

  CALL setpos (UNIT, wordaddress, icode)

  !----------------------------------------------------------------------
  ! [2.2]: MPP write.
  !----------------------------------------------------------------------

  ! Set up fake D1_addr record:

  DO i = 1, d1_list_len
    fake_d1_addr(i) = 0
  END DO

  item    = MOD(lookup(item_code, k), 1000)
  section = ( lookup(item_code, k) - item ) / 1000
  model   = lookup(model_code, k)

  halo_type = exppxi ( model, section, item, ppx_halo_type,       &
                       jcode, cmessage )

  grid_type = exppxi ( model, section, item, ppx_grid_type,       &
                       jcode, cmessage )

#if defined(FLDOP) || defined(MERGE) || defined(PPTOANC)
  grid_type=1
  halo_type=3
#endif
  IF (jcode /= 0) THEN
    WRITE(umMessage,*) ''
    CALL umPrint(umMessage,src='writflds')
    WRITE(umMessage,*) 'WRITFLDS: Failed to get PPXREF info.'
    CALL umPrint(umMessage,src='writflds')
    WRITE(umMessage,*) ''
    CALL umPrint(umMessage,src='writflds')
    WRITE(umMessage,*) '  Model ID:      ', model
    CALL umPrint(umMessage,src='writflds')
    WRITE(umMessage,*) '  Section:       ', section
    CALL umPrint(umMessage,src='writflds')
    WRITE(umMessage,*) '  Item:          ', item
    CALL umPrint(umMessage,src='writflds')
    WRITE(umMessage,*) '  Error code:    ', jcode
    CALL umPrint(umMessage,src='writflds')
    WRITE(umMessage,*) '  Error message: ', cmessage
    CALL umPrint(umMessage,src='writflds')
    WRITE(umMessage,*) ''
    CALL umPrint(umMessage,src='writflds')
    icode    = 1
    cmessage = 'WRITFLDS: Failed to get PPXREF info.'
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF

  fake_d1_addr(d1_object_type) = prognostic
  fake_d1_addr(d1_imodl)       = model
  fake_d1_addr(d1_section)     = section
  fake_d1_addr(d1_item)        = item
  fake_d1_addr(d1_halo_type)   = halo_type

  ! Grid type: for LBCs we need some special logic...
  IF (lookup(lbhem, k) == 99) THEN

    IF ( lookup(model_code, k) == atmos_im ) THEN
      fake_d1_addr(d1_grid_type) = ppx_atm_rim
    ELSE
      icode = 2
      WRITE(umMessage,*) ''
      CALL umPrint(umMessage,src='writflds')
      WRITE(umMessage,*) 'WRITFLDS: Cannot process LBC for model type',  &
                             lookup(model_code, k)
      CALL umPrint(umMessage,src='writflds')
      WRITE(umMessage,*) ''
      CALL umPrint(umMessage,src='writflds')
      cmessage = 'Cannot write LBCs for this model type.'
      IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
      RETURN
    END IF

  ELSE ! Not an LBC.

    fake_d1_addr(d1_grid_type) = grid_type

  END IF

  ! DEPENDS ON: get_fld_type
  fld_type = get_fld_type(grid_type)

  fake_d1_addr(d1_length)    = lasize(1, fld_type, halo_type) *   &
                               lasize(2, fld_type, halo_type)
  fake_d1_addr(d1_no_levels) = 1

  ! Write field:

  ! DEPENDS ON: write_multi
  CALL write_multi (                                              &
    UNIT,         d1(address), um_sector_ipts,                    &
                                                 ! in
    lenio,        locallen,                                       &
                                                 ! out
    lookup(1,k),  fixhd(12),   buf,                               &
                                                 ! in
    fake_d1_addr,                                                 &
                                                 ! in
    jcode,        cmessage )                     ! out

  address = address + locallen
  IF (locallen == 0) THEN
    address=address+lookup(lblrec,k)
  END IF

  ! Check for errors:

  IF (lenio /= um_sector_ipts) THEN

    WRITE(umMessage,*) ''
    CALL umPrint(umMessage,src='writflds')
    WRITE(umMessage,*) 'WRITFLDS: Error writing field number ', k,       &
                         ' on unit ', UNIT
    CALL umPrint(umMessage,src='writflds')
    WRITE(umMessage,*) ''
    CALL umPrint(umMessage,src='writflds')
    WRITE(umMessage,*) '  write_multi error code:    ', jcode
    CALL umPrint(umMessage,src='writflds')
    WRITE(umMessage,*) '  write_multi error message: ', cmessage
    CALL umPrint(umMessage,src='writflds')
    WRITE(umMessage,*) ''
    CALL umPrint(umMessage,src='writflds')
    icode    = jcode + 1
    cmessage = 'WRITFLDS: Error from write_multi.'
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN

  END IF

END DO ! k


IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE writflds

