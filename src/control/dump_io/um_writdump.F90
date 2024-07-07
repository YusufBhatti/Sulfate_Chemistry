! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    SUBROUTINE UM_WRITDUMP------------------------------------
!
!    Programming standard: Unified Model Documentation Paper No 3
!
!    Purpose: Writes out model dump on unit NFTOUT and checks model
!             and dump dimensions for consistency.
!
!    Documentation: Unified Model Documentation Paper No F3
!
!    Code Owner: Please refer to the UM file CodeOwners.txt
!    This file belongs in section: Dump I/O

SUBROUTINE um_writdump(nftout,fixhd,len_fixhd,    &
    inthd,len_inthd,                              &
    realhd,len_realhd,                            &
    levdepc,len1_levdepc,len2_levdepc,            &
    rowdepc,len1_rowdepc,len2_rowdepc,            &
    coldepc,len1_coldepc,len2_coldepc,            &
    flddepc,len1_flddepc,len2_flddepc,            &
    extcnst,len_extcnst,                          &
    dumphist,len_dumphist,                        &
    cfi1,len_cfi1,                                &
    cfi2,len_cfi2,                                &
    cfi3,len_cfi3,                                &
    lookup,len1_lookup,len2_lookup,               &
    mpp_lookup,mpp_len1_lookup,                   &
    buflen,                                       &
    submodel_id,                                  &
    n_objs_d1,d1_addr,                            &
    len_data,d1)

USE io_configuration_mod, ONLY: io_field_padding
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE io
USE IOS_Stash_Common
USE IOS_Constants
USE IOS_dump, ONLY:                            &
    IOS_dumpSetFieldParams,                     &
    IOS_dump_write_data
USE IOS_stash, ONLY:                           &
    IOS_stash_next_buffer,                      &
    IOS_stash_dispatch
USE ereport_mod, ONLY: ereport
USE umPrintMgr, ONLY: umPrint, umMessage, newline                 &
                    , printstatus, prstatus_diag, prstatus_normal

USE UM_ParVars
USE um_types,      ONLY: umFortranIntegerSize
USE decomp_params, ONLY: decomp_standard_atmos
USE Decomp_DB
USE writhead_mod
USE lookup_addresses
USE submodel_mod, ONLY: atmos_im
USE d1_array_mod, ONLY: d1_list_len, d1_section, d1_item,         &
                        d1_length, d1_grid_type, d1_no_levels
USE errormessagelength_mod, ONLY: errormessagelength
USE pr_look_mod, ONLY: pr_look
USE Packing_Codes_Mod, ONLY: PC_Cray32_Packing
IMPLICIT NONE

INTEGER ::                                                        &
 nftout                                                           &
               !IN Unit no of dump
,len_fixhd                                                        &
               !IN Length of fixed length header
,len_inthd                                                        &
               !IN Length of integer header
,len_realhd                                                       &
               !IN Length of real header
,len1_levdepc                                                     &
               !IN 1st dim of level dep consts
,len2_levdepc                                                     &
               !IN 2ndt dim of level dep consts
,len1_rowdepc                                                     &
               !IN 1st dim of row dep consts
,len2_rowdepc                                                     &
               !IN 2nd dim of row dep consts
,len1_coldepc                                                     &
               !IN 1st dim of column dep consts
,len2_coldepc                                                     &
               !IN 2nd dim of column dep consts
,len1_flddepc                                                     &
               !IN 1st dim of field dep consts
,len2_flddepc                                                     &
               !IN 2nd dim of field dep consts
,len_extcnst                                                      &
               !IN Length of extra constants
,len_dumphist                                                     &
               !IN Length of history block
,len_cfi1                                                         &
               !IN Length of comp field index 1
,len_cfi2                                                         &
               !IN Length of comp field index 2
,len_cfi3                                                         &
               !IN Length of comp field index 3
,len1_lookup                                                      &
               !IN 1st dim of lookup
,len2_lookup                                                      &
               !IN 2nd dim of lookup
,mpp_len1_lookup                                                  &
               !IN 1st dim of MPP lookup

,submodel_id                                                      &
               !IN submodel of dump
,n_objs_d1     !IN number of objects (3D fields) in D1

INTEGER ::                                                        &
 fixhd(len_fixhd)                                                 &
                                 !IN Fixed length header
,inthd(len_inthd)                                                 &
                                 !IN Integer header
,lookup(len1_lookup,len2_lookup)                                  &
                                 !IN PP lookup tables
,mpp_lookup(mpp_len1_lookup,len2_lookup)                          &
                                 !OUT PP lookup tables
,cfi1(len_cfi1+1)                                                 &
                  !IN Compressed field index no 1
,cfi2(len_cfi2+1)                                                 &
                  !IN Compressed field index no 2
,cfi3(len_cfi3+1) !IN Compressed field index no 3

INTEGER ::                                                        &
 len_data       !IN Length of real data
REAL ::                                                           &
 realhd(len_realhd)                                               &
                                      !IN Real header
,levdepc(1+len1_levdepc*len2_levdepc)                             &
                                      !IN Lev dep consts
,rowdepc(1+len1_rowdepc*len2_rowdepc)                             &
                                      !IN Row dep consts
,coldepc(1+len1_coldepc*len2_coldepc)                             &
                                      !IN Col dep consts
,flddepc(1+len1_flddepc*len2_flddepc)                             &
                                      !IN Field dep consts
,extcnst(len_extcnst+1)                                           &
                                      !IN Extra constants
,dumphist(len_dumphist+1)                                         &
                                      !IN History block
,d1(len_data)     !IN Real equivalence of data block

INTEGER ::                                                        &
  d1_addr(d1_list_len,n_objs_d1)  ! IN D1 addressing info

INTEGER ::                                                        &
 buflen         !IN Maximum length of single field in dump

! Local variables
INTEGER ::                                                        &
  start_block                                                     &
                     ! first word of field data in dump (unused)
, object_index                                                    &
                     ! pointer to entry in D1_ADDR
, level                                                           &
                     ! level number of multi-level field
, number_of_fields                                                &
                     ! total number of fields to write out
, field_start                                                     &
                     ! start address of a field in the file
, data_size                                                       &
                     ! number of words of data on disk for a field
, data_write_size                                                 &
                     ! total number of words to read for a field
, len_io                                                          &
                     ! number of words actually written
, orig_decomp                                                     &
                     ! original decomposition type
, local_len                                                       &
                     ! length of local field from buffout
, address                                                         &
                     ! address of field in local D1 array
, d1_start_read                                                   &
                     ! Same as start address, but may be modified
                     ! in case of zero-length fields at the end
                     ! of the d1-array
, k,i                                                             &
                     ! loop indicies
, number_of_data_words_in_memory                                  &
                                 ! unused return argument
, number_of_data_words_on_disk                                    &
                                 ! unused return argument
, disk_address       ! unused return argument

INTEGER :: buffer_start     ! start field in output buffer
INTEGER :: buffer_pos       ! current position in output buffer
INTEGER :: buffered_fields  ! count of buffered fields

LOGICAL ::                                                        &
  packed_field       ! TRUE if a field has been packed to 32 bits

REAL    :: IOSTAT    ! return from buffout

! Error reporting
INTEGER ::    icode       ! =0 normal exit; >0 error exit
CHARACTER(LEN=errormessagelength) :: cmessage    ! Error message
CHARACTER(LEN=*) :: RoutineName
PARAMETER (   RoutineName='UM_WRITDUMP')


INTEGER:: level_blocking

! Local arrays
REAL, ALLOCATABLE :: blocked_out_buf(:)

REAL :: temp_buf(buflen)   ! For 32 bit packing

INTEGER :: IOS_Q_Slot ! Opaque handle for referencing async operations
INTEGER :: IOS_packing_flag
INTEGER :: IOS_packing_type
INTEGER :: IOS_fullfield_flag
INTEGER :: Get_Fld_Type
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!--------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)


level_blocking = 1
IF (isUsingAsyncDumps()) THEN
  level_blocking=IOS_AsyncMaxFieldsInPack
  IF (mype == 0 .AND. printstatus >= prstatus_diag) THEN
    WRITE(umMessage,'(A,I3)')'UM_WritDump: Output blocking factor:', &
        level_blocking
    CALL umPrint(umMessage,src='um_writdump')
  END IF
END IF

! blocked_out_buf is only used by pe0 when the IOS is not active
IF (.NOT. isUsingAsyncDumps() .AND. mype == 0) THEN
  ALLOCATE(blocked_out_buf((((buflen+io_field_padding)/io_field_padding)       &
      * io_field_padding) * level_blocking))
ELSE
  ALLOCATE(blocked_out_buf(0))
END IF

icode=0
cmessage=''

IF (mype  ==  0) THEN
  CALL umPrint('',src='um_writdump')
  WRITE(umMessage,'('' WRITING UNIFIED MODEL DUMP ON UNIT'',I3)')nftout
  CALL umPrint(umMessage,src='um_writdump')
  WRITE(umMessage,'('' #####################################'')')
  CALL umPrint(umMessage,src='um_writdump')
  CALL umPrint('',src='um_writdump')
END IF

! Select the relevant decomposition type for this dump

orig_decomp=current_decomp_type

IF (submodel_id  ==  atmos_im) THEN
  IF (current_decomp_type  /=  decomp_standard_atmos)             &
  CALL change_decomposition(decomp_standard_atmos,icode)

ELSE  ! unsupported decomposition type
  WRITE(umMessage,'(A,A,I3)')                                             &
    'UM_WRITEDUMP : Could not change to decomposition required ', &
    'for submodel type ',submodel_id
  CALL umPrint(umMessage,src='um_writdump')
  icode=1
  cmessage='Unsupported submodel for MPP code'

  CALL ereport(routinename,icode,cmessage)
END IF

IF (icode  /=  0) THEN
  icode=2
  WRITE(umMessage,'(A,A)')                                                &
    'UM_WRITDUMP : Error - Could not set decomposition ',         &
    'for selected submodel.'
  CALL umPrint(umMessage,src='um_writdump')
  cmessage='Unsupported decomposition selected for MPP code'

  CALL ereport(routinename,icode,cmessage)
END IF

! Reset the disk addresses and lengths for well-formed I/O
! DEPENDS ON: set_dumpfile_address
CALL set_dumpfile_address(fixhd, len_fixhd,                       &
                          lookup, len1_lookup,                    &
                          len2_lookup,                            &
                          number_of_data_words_in_memory,         &
                          number_of_data_words_on_disk,           &
                          disk_address)

! Write out the header records
CALL writhead(nftout,fixhd,len_fixhd,                             &
              inthd,len_inthd,                                    &
              realhd,len_realhd,                                  &
              levdepc,len1_levdepc,len2_levdepc,                  &
              rowdepc,len1_rowdepc,len2_rowdepc,                  &
              coldepc,len1_coldepc,len2_coldepc,                  &
              flddepc,len1_flddepc,len2_flddepc,                  &
              extcnst,len_extcnst,                                &
              dumphist,len_dumphist,                              &
              cfi1,len_cfi1,                                      &
              cfi2,len_cfi2,                                      &
              cfi3,len_cfi3,                                      &
              lookup,len1_lookup,len2_lookup,len_data,            &
              umFortranIntegerSize()*8,                           &
              start_block,icode,cmessage)

IF (icode  /=  0) THEN
  WRITE(umMessage,'(A,A,I3)') 'UM_WRITDUMP : Error writing dump header ', &
             'on unit ',nftout
  CALL umPrint(umMessage,src='um_writdump')
  WRITE(umMessage,'(A,I4,A,A)') 'Return code from WRITHEAD was ',icode,   &
             ' and error message was ',cmessage
  CALL umPrint(umMessage,src='um_writdump')
  icode=3
  cmessage='Error writing dump header'

  CALL ereport(routinename,icode,cmessage)
END IF


IF (fixhd(160)  >   0) THEN ! If there is data to write

  ! Loop over fields and write from D1

  number_of_fields= fixhd(152)
  address         = 1
  object_index    = 1
  level           = 1
  buffered_fields = 0
  buffer_start    = 0
  buffer_pos      = 1


  ! If we are using the async stash approach then we prepare a new
  ! dispatch buffer. We will fill this with all output arising from
  ! the subsequent levels loop.

  IOS_Q_Slot=-1 ! default flag for normal general_gather_field
  IF (isUsingAsyncDumps()) THEN
    IF (printstatus >= prstatus_diag) THEN
      CALL umPrint('UM_WritDump: Getting 1st buffer for async dump')
    END IF
    CALL IOS_stash_next_buffer( nftout, IOS_Q_Slot, &
        bufType=IOS_Action_StashWriteDumpData)
    IF (IOS_Q_Slot<=0) THEN

      CALL ereport(routinename,icode,"Cannot lcoate an async slot")
    END IF
    IF (printstatus >= prstatus_diag) THEN
      WRITE(umMessage,'(A,I4)')'UM_WritDump: New buffer slot=',IOS_Q_Slot
      CALL umPrint(umMessage,src='um_writdump')
    END IF
  END IF


  DO k=1,number_of_fields  ! loop over fields to write out

    mpp_lookup(p_lblrec,k)=0
    mpp_lookup(p_naddr,k)=address

    IF (lookup(lblrec,k)  >   0) THEN  ! if there's data in
                                       ! the field

      ! Check that DATA_TYPE is valid no: +/-1 to +/-3
      IF (( ABS(lookup(data_type,k))  >=  1) .AND.                &
          ( ABS(lookup(data_type,k))  <=  3)) THEN

        ! OK - we're gonna write this one out. Setup buffering parameters
        buffered_fields = buffered_fields + 1

        IF (buffer_start == 0) THEN
          buffer_start = k
        END IF

        ! Check that the "blocked_out_buf" array is big enough for this field

        IF (lookup(lblrec,k)  >   buflen) THEN
          WRITE(umMessage,*)                                              &
            'UM_WRITDUMP : Field length longer than buffer'
          CALL umPrint(umMessage,src='um_writdump')
          WRITE(umMessage,*) 'Field :',k
          CALL umPrint(umMessage,src='um_writdump')
          WRITE(umMessage,*) 'Field length on disk ',lookup(lblrec,k)
          CALL umPrint(umMessage,src='um_writdump')
          WRITE(umMessage,*) 'Buffer length ',buflen
          CALL umPrint(umMessage,src='um_writdump')
          icode=4
          cmessage='Field length longer than buffer'

          CALL ereport(routinename,icode,cmessage)
        END IF

        ! Set "packed_field" to .TRUE. if 32bit packing has been used
        packed_field=(MOD((lookup(lbpack,k)),10) == PC_Cray32_Packing)

        ! Set up the location of the field on disk and how much data
        ! needs to be written
        field_start=lookup(lbegin,buffer_start)

        ! data_size contains the number of words of data used to store
        ! the field on disk
        IF (packed_field) THEN
          data_size=(lookup(lblrec,k)+1)/2
        ELSE
          data_size=lookup(lblrec,k)
        END IF

        ! data_write_size contains the number of words of data that need to
        ! be written out for a field. Each field has extra words of dummy data
        ! added at the end to ensure each field starts on a disk sector
        ! boundary.
        data_write_size=lookup(lbnrec,k)

        IF (isUsingAsyncDumps()) THEN
          IOS_packing_flag=IOS_no_packing
          IOS_packing_type=IOS_no_packing
          IF (packed_field) THEN
            IF (lookup(data_type,k)  ==  1) THEN
              IOS_packing_flag=IOS_packing
              IOS_packing_type=IOS_packing_type_pack21
            END IF
          END IF
          CALL IOS_dumpSetFieldParams( &
               IOS_packing_flag,       &
               IOS_packing_type,       &
               data_write_size,        &
               field_start)
        END IF

        d1_start_read = address

        ! Trap case of where address may attempt to access beyond length
        ! of d1 array. This could happen in the valid case of a zero length
        ! land-packed field on the end of the d1-array.
        IF (d1_start_read > len_data) THEN
          IF (d1_addr(d1_length, object_index) > 0) THEN
            icode = 5
            WRITE(umMessage,'(5(A,I0))')                                       &
              'Model is about to attempt to access an index [',                &
              d1_start_read, ']' //                                   newline//&
              'which is beyond the length of the d1 array [', len_data, '].' //&
              newline//                                                        &
              "Field object: ", k,                                    newline//&
              "Stash number: ", d1_addr(d1_section,object_index), ', ',        &
              d1_addr(d1_item,object_index)
            CALL ereport(routinename,icode,cmessage)

          ELSE
            ! A zero sized field at the end of the d1-array
            ! set the d1_start_address to equal the length
            ! of the d1 array to prevent out of bounds access
            ! This can happen if a land-packed field is on the
            ! end of the d1 array
            d1_start_read = len_data

          END IF
        END IF

        ! Gather the field (into tempbuf) onto PE 0 for packing and output
        ! DEPENDS ON: general_gather_field
        CALL general_gather_field(                                       &
            d1(d1_start_read),temp_buf,local_len,lookup(lblrec,k),1,     &
            d1_addr(1,object_index),0,                                   &
            IOS_Q_Slot,icode,cmessage)

        IF (icode  /=  0) THEN
          WRITE(umMessage,*) 'WRITEDUMP: Call to GENERAL_GATHER_FIELD ', &
                     'failed'
          CALL umPrint(umMessage,src='um_writdump')
          WRITE(umMessage,*) 'Return code was ',icode
          CALL umPrint(umMessage,src='um_writdump')
          WRITE(umMessage,*) 'Error message was ',cmessage
          CALL umPrint(umMessage,src='um_writdump')
          WRITE(umMessage,*) 'Field number ',lookup(item_code,k)
          CALL umPrint(umMessage,src='um_writdump')
          WRITE(umMessage,*) 'Dimensions ',lookup(lbnpt,k),              &
                                 ' x ',lookup(lbrow,k)
          CALL umPrint(umMessage,src='um_writdump')
          WRITE(umMessage,*) 'Grid type ',                               &
                      d1_addr(d1_grid_type,object_index)
          CALL umPrint(umMessage,src='um_writdump')
          WRITE(umMessage,*) 'Field was not written out'
          CALL umPrint(umMessage,src='um_writdump')

          icode=300
          cmessage='Failure to gather field'

          CALL ereport(routinename,icode,cmessage)
        END IF

        IF (.NOT. isUsingAsyncDumps()) THEN

          !  Do the packing of the data

          IF (mype  ==  0) THEN
            ! Does this field need to be compressed?
            IF ( packed_field ) THEN
              IF (lookup(data_type,k)  ==  1) THEN
                ! DEPENDS ON: pack21
                CALL pack21(lookup(lblrec,k),temp_buf,              &
                    blocked_out_buf(buffer_pos))
              END IF
            ELSE ! no compression required - just do a copy
              DO i=1,lookup(lblrec,k)
                blocked_out_buf(buffer_pos + i - 1) = temp_buf(i)
              END DO
            END IF
          END IF ! am I PE 0 ?

        END IF ! async stash

        ! If the field was compressed for writing on disk, we need to compress
        ! and expand the field in memory. This ensures the same field exists in
        ! memory that would exist if this dump was read back in.

        IF ( packed_field ) THEN
          IF (lookup(data_type,k)  ==  1) THEN
            ! DEPENDS ON: pack21
            CALL pack21(local_len,d1(d1_start_read),temp_buf)
            ! DEPENDS ON: expand21
            CALL expand21(local_len,temp_buf,d1(d1_start_read))
          END IF
        END IF

        ! Increment position in buffer
        buffer_pos = buffer_pos + data_write_size

        ! Now write out the buffered fields if we've reached the limit
        IF (buffered_fields == level_blocking .OR.                &
            k == number_of_fields) THEN

          ! Set the position in the dump

          CALL setpos(nftout,field_start,icode)
          IF (icode  /=  0) THEN
            WRITE(umMessage,*)                                            &
            'UM_WRITDUMP - SETPOS failed to move file pointer ',  &
            ' to ',field_start,' on unit ',nftout
            CALL umPrint(umMessage,src='um_writdump')
            WRITE(umMessage,*) 'SETPOS returned error code ',icode
            CALL umPrint(umMessage,src='um_writdump')
            icode=6
            cmessage='SETPOS failed while writing dump.'

            CALL ereport(routinename,icode,cmessage)
          END IF

          IF (isUsingAsyncDumps()) THEN
            ! Now we signal that the packed buffer is available for
            ! transmit to the IO server. As the trasmission is asyncronous
            ! We do not test for completion, the buffer will not be recycled
            ! until message delivery.

            ! Release data buffers to the IO server

            CALL IOS_stash_dispatch( IOS_Q_Slot )

            ! Send the control message to the IO Server.
            CALL IOS_dump_write_data( IOS_Q_Slot )

            IF (k /= number_of_fields) THEN
              ! Prepare a new buffer if there are still more fields
              IF (printstatus >= prstatus_diag) THEN
                WRITE(umMessage,'(A,I4,A,I4)')                      &
                    'UM_WritDump: New buffer for field=', &
                    k,' of ',number_of_fields
                CALL umPrint(umMessage,src='um_writdump')
              END IF
              CALL IOS_stash_next_buffer( nftout, IOS_Q_Slot, &
                  bufType=IOS_Action_StashWriteDumpData)
              IF (IOS_Q_Slot<=0) THEN

                CALL ereport(routinename,icode,"Cannot lcoate an async slot")
              END IF
              IF (printstatus >= prstatus_diag) THEN
                WRITE(umMessage,'(A,I4)')             &
                    'UM_WritDump: New buffer slot=',IOS_Q_Slot
                CALL umPrint(umMessage,src='um_writdump')
              END IF
            END IF
          ELSE

            ! Write out the data on PE 0
            IF (mype == 0) THEN

              IF ( packed_field ) THEN
                IF (lookup(data_type,k)  ==  1) THEN
                  ! Data is packed using CRAY 32 bit method - note that we need to write
                  ! out 2*ISIZE 32 bit words using BUFFO32_f77 (because the array is 64 bit)
                  ! DEPENDS ON: buffout32_f77
                  CALL buffout32_f77(nftout,blocked_out_buf,2*(buffer_pos - 1),&
                                      len_io,IOSTAT)
                  ! And then halve LEN_IO to satisfy tests against ISIZE
                  len_io = len_io/2
                END IF
              ELSE
                ! For non-packed data
                CALL buffout(nftout,blocked_out_buf,buffer_pos - 1,      &
                    len_io,IOSTAT)
              END IF

              IF ((IOSTAT  /=  -1.0) .OR.                           &
                  (len_io  /=  (buffer_pos - 1 ))) THEN
                WRITE(umMessage,*) 'WRITEDUMP: Error in call to ',          &
                    'BUFFOUT'
                CALL umPrint(umMessage,src='um_writdump')
                WRITE(umMessage,*) 'Field : ', k
                CALL umPrint(umMessage,src='um_writdump')
                WRITE(umMessage,*) 'LEN_IO : ',len_io
                CALL umPrint(umMessage,src='um_writdump')
                WRITE(umMessage,*) 'IOSTAT : ',IOSTAT
                CALL umPrint(umMessage,src='um_writdump')

                icode=400
                cmessage='Failure writing out field'

                CALL ereport(routinename,icode,cmessage)
              END IF
            END IF  ! mype == 0
          END IF
          ! Reset data for indexing in buffer
          buffer_pos = 1
          buffer_start = 0
          buffered_fields = 0
        END IF ! Writing out data


        mpp_lookup(p_lblrec,k)=local_len
        address=address+local_len

        IF (icode  /=  0) THEN
          WRITE(umMessage,*)                                              &
            'UM_WRITDUMP - Error while attempting to ',           &
            'write field ',k,' of ',number_of_fields,             &
            ' from unit ',nftout
          CALL umPrint(umMessage,src='um_writdump')
          WRITE(umMessage,*) 'Return code from WRITE_MULTI was ',icode,   &
            'and error message was ',cmessage
          CALL umPrint(umMessage,src='um_writdump')
          WRITE(umMessage,*) 'Field Information: '
          CALL umPrint(umMessage,src='um_writdump')
          WRITE(umMessage,*) 'Section ',d1_addr(d1_section,object_index), &
                     ' Item ',d1_addr(d1_item,object_index)
          CALL umPrint(umMessage,src='um_writdump')
          WRITE(umMessage,*) 'Disk address : ',field_start
          CALL umPrint(umMessage,src='um_writdump')
          WRITE(umMessage,*) 'D1 address : ',d1_start_read
          CALL umPrint(umMessage,src='um_writdump')
          WRITE(umMessage,*) 'Number of words requested : ',              &
                     data_write_size
          CALL umPrint(umMessage,src='um_writdump')
          WRITE(umMessage,*) 'Number of words returned : ',len_io
          CALL umPrint(umMessage,src='um_writdump')

          icode=7
          cmessage='Error writing field to dump'

          CALL ereport(routinename,icode,cmessage)

        END IF ! If an error was detected writing the field


      ELSE ! invalid data type

        IF (( fixhd(5)  <   6) .OR.                               &
                                     ! Not AC, Var or Cx
            ( fixhd(5)  >   8)) THEN

          CALL pr_look(lookup,k)
        END IF

        WRITE(umMessage,*) 'UM_WRITDUMP : Failure for field ',k,' of ',   &
                   number_of_fields
        CALL umPrint(umMessage,src='um_writdump')
        WRITE(umMessage,*) 'LOOKUP(DATA_TYPE,K)= ',lookup(data_type,k)
        CALL umPrint(umMessage,src='um_writdump')
        icode=8
        cmessage='Invalid data type ( LOOKUP(DATA_TYPE,K) )'

        CALL ereport(routinename,icode,cmessage)

      END IF ! Check the LOOKUP(DATA_TYPE,K) is valid

#if defined (UTILIO)
      IF (printstatus >= prstatus_diag) THEN
        IF ((fixhd(5)  <   6) .OR.                                &
                                    ! Not AC/Var
            (fixhd(5)  >   8)) THEN !     Obs/Cx

          ! Print out header and summary of data field

          CALL pr_look(lookup,k)
          IF (fixhd(5)  /= 5 ) THEN   !  Skip if boundary dataset
            IF (lookup(data_type,k) == 1) THEN  !  Real
              ! DEPENDS ON: pr_rfld
              CALL pr_rfld(lookup,lookup,d1(lookup(naddr,k)),k)
            ELSE IF (lookup(data_type,k) == 2) THEN  !  Integer
              ! DEPENDS ON: pr_ifld
              CALL pr_ifld(lookup,lookup,d1(lookup(naddr,k)),k)
            ELSE IF (lookup(data_type,k) == 3) THEN  !  Logical
              ! DEPENDS ON: pr_lfld
              CALL pr_lfld(lookup,lookup,len1_lookup,             &
                   d1(lookup(naddr,k)),k)
            END IF ! Test for different data types
          END IF ! test for boundary dataset
        END IF ! test for field type
      END IF ! PrintStatus >= PrStatus_Diag
#endif
    END IF ! If there was data in the field

    level=level+1
    IF (level  >   d1_addr(d1_no_levels,object_index)) THEN
      level=1
      object_index=object_index+1
    END IF

  END DO ! K : loop over fields to write out

  IF (printstatus  >=  prstatus_normal) THEN
    IF (mype  ==  0) THEN
      WRITE(umMessage,*) 'Data successfully written'
      CALL umPrint(umMessage,src='um_writdump')
      WRITE(umMessage,*) fixhd(161),' words written to unit ',nftout
      CALL umPrint(umMessage,src='um_writdump')
      IF ((fixhd(5)  >=  6) .AND.                                 &
                                  ! AC/Var
          (fixhd(5)  <=  8)) THEN ! Obs/Cx
        WRITE(umMessage,*) '(Observational data)'
        CALL umPrint(umMessage,src='um_writdump')
      ELSE
        WRITE(umMessage,*) '(Model data)'
        CALL umPrint(umMessage,src='um_writdump')
      END IF
    END IF ! IF (mype  ==  0)
  END IF ! IF (PrintStatus  >=  PrStatus_Normal)

END IF ! IF (FIXHD(160)  >   0)

! Reset to original decomposition type
CALL change_decomposition(orig_decomp,icode)
DEALLOCATE(blocked_out_buf)

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE um_writdump
