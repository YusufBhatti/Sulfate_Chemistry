! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
!      Subroutine:
!      ACUMPS
!
!      Purpose:
!      To accumulate partial sums of climate mean tagged diagnostics
!      and create dumps containing them. Also to overwrite the D1
!      diagnostic with the partial sum for use by MEANPS. This
!      saves MEANPS having to reread the partial sum dump.
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

SUBROUTINE acumps(                                                &
  n_objs_d1,d1_addr                                               &
  ,len_data,d1                                                    &
  ,maxsize,means_total                                            &
  ,flag,nftin,nftout,lclimrealyr,meanlev                          &
  ,i_month,i_year                                                 &
  ,head_out,head_len,head_size,                                   &
  timestep,cmitems,fixhd12,                                       &
  icode,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE io
USE UM_ParVars
USE UM_ParCore,   ONLY: mype
USE lookup_addresses
USE d1_array_mod, ONLY: d1_list_len, d1_imodl, d1_section,        &
                        d1_item, d1_address, d1_length,           &
                        d1_no_levels, d1_south_code,              &
                        d1_north_code, d1_east_code,              &
                        d1_west_code, d1_gridpoint_code,          &
                        d1_proc_no_code, d1_object_type,          &
                        d1_halo_type, d1_grid_type
USE stash_array_mod, ONLY: totitems, stlist
USE stparam_mod, ONLY: st_macrotag, st_d1pos, s_modl, &
                       st_dump_level_output_length

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE missing_data_mod, ONLY: rmdi

USE errormessagelength_mod, ONLY: errormessagelength

USE setperlen_mod, ONLY: setperlen

IMPLICIT NONE

INTEGER ::                                                        &
  n_objs_d1               !IN No objects in D1 array

INTEGER ::                                                        &
  d1_addr(d1_list_len,n_objs_d1) !IN Addressing of D1 array

INTEGER ::                                                        &
  maxsize,                                                        &
                          ! IN dimension of largest data field
  len_data,                                                       &
                          ! IN Length of model data
  flag,                                                           &
                          ! IN Flag for reading partial sum dump
  nftin,                                                          &
                          ! IN Unit no for reading partial sums
  nftout,                                                         &
                          ! IN Unit no for writing partial sums
  icode,                                                          &
                          ! OUT Return code; successful=0
                          !                  error>0
  meanlev,                                                        &
                          ! IN level of climate meaning
  means_total,                                                    &
                          ! IN Indicates a meaning period
  i_month,                                                        &
                          ! IN Current model time (months)
  i_year,                                                         &
                          ! IN Current model time (years)
  fixhd12,                                                        &
                          ! IN Version of model
  cmitems,                                                        &
                          ! IN Number of items being meaned
  timestep                ! IN Submodel timestep
!
CHARACTER(LEN=errormessagelength) ::                              &
       cmessage             ! OUT Error message if ICODE>0
!
REAL ::                                                            &
   d1(len_data)            ! IN/OUT Real equivalence of data block
!
LOGICAL ::                                                         &
   lclimrealyr             ! IN Real-period climate meaning


INTEGER ::                                                        &
  head_len                                                        &
  ,head_size

INTEGER ::                                                        &
  head_out(head_len,totitems)                                     &
                              ! IN Header contains packing
                          !    info for output ps file
  ,head_buf(head_size)

! Header formatted as follows:
! HEAD_OUT(1,*): No of words per level in field
! HEAD_OUT(2,*): 2 for packed, 1 for unpacked
! HEAD_OUT(3,*): No of words per level on disk

!
!      Local variables
!
INTEGER ::                                                        &
  i,j,k                                                           &
                          ! Loop indices
  ,len_io                                                         &
                          ! Actual IO length
  ,citems                                                         &
                          ! Count variable
  ,periodlen                                                      &
                          ! Current meaning period in days
  ,tag                                                            &
                          ! Stash tag
  ,ptd1                                                           &
                          ! Pointer to D1_ADDR information
  ,address                                                        &
                          ! Address in local D1
  ,levels                                                         &
                          ! Number of levels per diagnostic
  ,length                                                         &
                          ! Length of each level in local D1
  ,global_length                                                  &
                          ! Length of global field
  ,offset                 ! Indexing offset for WORK array
!
INTEGER ::                                                        &
  header(2)                                                       &
                          ! Initial header
  ,head_in(head_len,cmitems) ! Packing info for input ps file
                          ! Will differ from HEAD_OUT if packing
                          ! codes have changed mid-run

REAL ::                                                           &
  IOSTAT,                                                         &
                          ! IO error code
  realperiodlen,                                                  &
                          ! explicitly real equivalent
                          ! of PERIODLEN
  cwork,                                                          &
                          ! Accumulative sum for WORK array
  ckwork,                                                         &
                          ! Checksum for WORK array
  ckworko                 ! Packed CKWORK

REAL, PARAMETER :: small_number = TINY(1.0)
                          ! Small number for comparing REALS

!
!      Local arrays
!
REAL ::                                                           &
  d1_data(maxsize)        ! Work area for fields
REAL ::                                                           &
  work(maxsize+4)        ! Work area and IO buffer

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ACUMPS'

!
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
IF (icode /= 0) GO TO 9999

! Align for well-formed io

! Arrays sent to BUFFIN/BUFFOUT need to be cache aligned for well-formed
! io to work. The following adds an offset to the index of
! WORK which resolves the problem.
! This is an offset from 0, so an offset of 1 really means no offset!
offset=1

!   Set up variables needed for weighting accumulations if real-period
!   climate meaning is selected. Partial sums are normalised elsewhere.

IF (lclimrealyr) THEN
  CALL setperlen(meanlev,i_month,i_year,periodlen)
  realperiodlen=REAL(periodlen)
END IF

! STEP 1: Read in headers of previous partial sum and write out
!         header of new.

IF (flag /= 1) THEN       ! PS data exist on disk
  ! Read headers for input partial sum file

  CALL buffin(nftin,head_buf,head_size,len_io,IOSTAT)
  IF (IOSTAT /= -1.0 .OR. len_io /= head_size) THEN
    WRITE(umMessage,*)'ACUMPS: Error reading header: IO code ',           &
      IOSTAT,' on unit ',nftin
    CALL umPrint(umMessage,src='acumps')
    WRITE(umMessage,*)'Words requested ',head_size,                       &
      ' Words read ',len_io
    CALL umPrint(umMessage,src='acumps')
    icode=1
    cmessage='ACUMPS: BUFFIN error - see output'
    GO TO 9999
  END IF
  ! Transfer header information from buffer to header arrays
  header(1)=head_buf(1) ! Timestep of creation
  header(2)=head_buf(2) ! Number of records
  k=3
  DO i=1,cmitems
    DO j=1,head_len
      head_in(j,i)=head_buf(k)
      k=k+1
    END DO
  END DO

  IF (header(1) >= timestep .OR. header(2) /= cmitems) THEN
    WRITE(umMessage,*)'ACUMPS1: Partial sum file inconsistent'
    CALL umPrint(umMessage,src='acumps')
    WRITE(umMessage,*)'PS file holds ', &
        header(2),' items and written at STEP ',header(1)
    CALL umPrint(umMessage,src='acumps')
    WRITE(umMessage,*)'Expected timestep should be < ',timestep
    CALL umPrint(umMessage,src='acumps')
    WRITE(umMessage,*)'Expected number of items ',cmitems
    CALL umPrint(umMessage,src='acumps')
    cmessage='ACUMPS1: Partial sum file inconsistent. See Output'
    icode=2
    GO TO 9999
  END IF
ELSE
  ! No input sum, so initialise header array
  DO i=1,head_size
    head_buf(i)=0
  END DO
END IF

! Write headers for new partial sum file
! Transfer information to io buffer
head_buf(1)=timestep
head_buf(2)=cmitems
k=3
DO i=1,cmitems
  DO j=1,head_len
    head_buf(k)=head_out(j,i)
    k=k+1
  END DO
END DO

CALL buffout(nftout,head_buf,head_size,len_io,IOSTAT)
IF (IOSTAT /= -1.0 .OR. len_io /= head_size) THEN
  WRITE(umMessage,*)'ACUMPS: Error writing header: IO code ',             &
    IOSTAT,' on unit ',nftout
  CALL umPrint(umMessage,src='acumps')
  WRITE(umMessage,*)'Words requested ',head_size,                         &
    ' Words written ',len_io
  CALL umPrint(umMessage,src='acumps')
  icode=4
  cmessage='ACUMPS: BUFFOUT error - see output'
  GO TO 9999
END IF

! STEP 2 : Loop over all STASH items. For each tagged item, gather
!          current data to D1_DATA array, read partial sum into WORK
!          array (if there is a partial sum), sum the two and write
!          out to new partial sum file.
!           Also, if this is a meaning period, overwrite the field
!          in D1 with the complete sum, to be picked up by MEANPS.

!     Start of loop over STASH items
citems=0
DO k=1,totitems
  tag=stlist(st_macrotag,k)/1000
  ptd1=stlist(st_d1pos,k)
  IF (tag/=0) THEN
    IF (stlist(s_modl,k)==d1_addr(d1_imodl,ptd1)) THEN
      ! Object tagged for climate meaning and in relevant internal model
      address=d1_addr(d1_address,ptd1)
      levels=d1_addr(d1_no_levels,ptd1)
      length=d1_addr(d1_length,ptd1)/levels
      global_length=stlist(st_dump_level_output_length,k)
      citems=citems+1
      DO j=1,levels
        ! Copy current field from D1 to D1_DATA
        ! by gathering full field to pe0
        ! DEPENDS ON: general_gather_field
        CALL general_gather_field(                                  &
          d1(address),d1_data,length,                               &
          global_length,1,                                          &
          d1_addr(1,ptd1),0,                                        &
          -1,icode,cmessage)
        IF (icode /= 0) GO TO 9999
        DO i=global_length+1,maxsize
          d1_data(i)=0.0
        END DO
        ! Set initial value for CKWORKO and CWORK
        cwork=0.0
        ckworko=0.0
        ! If partial sum exists on disk, read it in and add to current field
        IF (flag /= 1) THEN ! PS data exist on disk
          ! Read in one level of partial sum field
          IF (head_in(2,citems)  ==  2) THEN

            ! Data is packed using CRAY 32 bit method - note that we need to read
            ! in 2*HEAD_IN(3,CITEMS) 32 bit words using BUFFIN32_f77 (because the
            ! array is 64 bit)

            ! DEPENDS ON: buffin32_f77
            CALL BUFFIN32_f77(nftin,work(offset:), &
                 2*head_in(3,citems),len_io,IOSTAT)

            ! And then halve LEN_IO to satisfy tests against HEAD_IN(3,CITEMS)
            len_io = len_io/2
          ELSE ! For non-packed data
            CALL buffin(nftin,work(offset:),head_in(3,citems)        &
             ,len_io,IOSTAT)
          END IF
          IF (IOSTAT /= -1.0 .OR. len_io /= head_in(3,citems)) THEN
            WRITE(umMessage,*)'ACUMPS: Error reading partial sum IO code ', &
              IOSTAT,' on unit ',nftin
            CALL umPrint(umMessage,src='acumps')
            WRITE(umMessage,*)'Words requested ',head_in(3,citems),         &
              ' Words read ',len_io
            CALL umPrint(umMessage,src='acumps')
            icode=6
            cmessage='ACUMPS: BUFFIN error - see output'
            GO TO 9999
          END IF
          IF (mype == 0) THEN
            ! Valid data exists on pe0 only
            ! Unpack if data on disk was packed
            IF (head_in(2,citems) == 2) THEN
              ! DEPENDS ON: expand32b
              CALL expand32b(global_length+1,work(offset),fixhd12)
              ! Calculate a checksum, this is the mean of the field
              DO i=1,global_length
                ckwork=REAL(i-1)*(ckwork - work(i+offset-1))/REAL(i) +   &
                       work(i+offset-1)
              END DO
              ! Pack and umpack checksum to force it losing precision in order to do
              ! the comparison
              ! DEPENDS ON: pack21
              CALL pack21(1,ckwork,ckworko)
              ! DEPENDS ON: expand32b
              CALL expand32b(1,ckworko,fixhd12)
            ELSE
              ! Calculate a checksum, this is the mean of the field
              DO i=1,global_length
                ckworko=REAL(i-1)*(ckworko - work(i+offset-1))/REAL(i) + &
                        work(i+offset-1)
              END DO
            END IF
            IF (ABS(ckworko - work(global_length+offset)) > small_number ) THEN
              WRITE(umMessage,*)'ERROR: checksum failure in climate mean'
              CALL umPrint(umMessage,src='acumps')
              WRITE(umMessage,*)'Section ',d1_addr(d1_section,ptd1), &
                 ' item ',d1_addr(d1_item,ptd1)
              CALL umPrint(umMessage,src='acumps')
              WRITE(umMessage,*) 'This can be due to invalid values in field, ', &
                 'or corruption of partial sum file'
              CALL umPrint(umMessage,src='acumps')
              WRITE(umMessage,*)'Remove or fix diagnostic, and rerun'
              CALL umPrint(umMessage,src='acumps')
              icode=4
              cmessage='ACUMPS: Diagnostic error. See output for item no.'
              GO TO 9999
            END IF
            ! Sum with field in D1 - Scale data if 365 day calendar
            IF (lclimrealyr) THEN
              DO i=1,global_length
                IF (work(i+offset-1) == rmdi) THEN
                  d1_data(i)=rmdi
                ELSE
                  d1_data(i)=work(i+offset-1)+                      &
                    (realperiodlen*d1_data(i))
                END IF
              END DO
            ELSE
              ! 360 day calendar
              DO i=1,global_length
                IF (work(i+offset-1) == rmdi) THEN
                  d1_data(i)=rmdi
                ELSE
                  d1_data(i)=work(i+offset-1)+d1_data(i)
                END IF
              END DO
            END IF
          END IF
        ELSE
          ! First data for this period - no partial sum to add
          IF (lclimrealyr) THEN
            ! Scale initial data if 365 day calendar
            IF (mype == 0) THEN
              DO i=1,global_length
                IF (d1_data(i) /= rmdi) THEN
                  d1_data(i)=realperiodlen*d1_data(i)
                END IF
              END DO
            END IF
          END IF
        END IF ! End of adding PS data

        !         Write out sum to PS file
        IF (mype == 0) THEN
          ! Copy data to WORK array, packing if necessary
          IF (head_out(2,citems) == 2) THEN
            DO i=head_out(1,citems),maxsize
              work(i+offset-1)=0.0
            END DO
            ! DEPENDS ON: pack21
            CALL pack21(global_length+1,d1_data,                    &
                        work(offset))
            ! DEPENDS ON: expand32b
            CALL expand32b(global_length+1,work(offset),fixhd12)
            ! Calculate a checksum, this is the mean of the field
            DO i=1,global_length
              cwork=REAL(i-1)*(cwork - work(i+offset-1))/REAL(i)+ &
                     work(i+offset-1)
            END DO
            work(global_length+offset)=cwork
            ! DEPENDS ON: pack21
            CALL pack21(global_length+1,work(offset),               &
                        work(offset) )

            ! If data not packed, calculate checksum straight away
          ELSE
            DO i=1,global_length
              work(i+offset-1)=d1_data(i)
              cwork=REAL(i-1)*(cwork - work(i+offset-1))/REAL(i)+ &
                     work(i+offset-1)
            END DO
            work(global_length+offset)=cwork
          END IF
        END IF


        ! Output partial sum to file
        IF (head_out(2,citems)  ==  2) THEN

          ! Data is packed using CRAY 32 bit method - note that we need to write
          ! out 2*HEAD_OUT(3,CITEMS) 32 bit words using BUFFOUT32_F77 (because
          ! we have a supplied 64 bit array)

          ! DEPENDS ON: buffout32_f77
          CALL BUFFOUT32_f77(nftout,work(offset:),                  &
               2*head_out(3,citems),len_io,IOSTAT)
          ! And then halve LEN_IO to satisfy tests against HEAD_OUT(3,CITEMS)
          len_io = len_io/2
        ELSE
          ! For non-packed data

          CALL buffout(nftout,work(offset:),                        &
            head_out(3,citems),len_io,IOSTAT)
        END IF
        IF (IOSTAT /= -1.0 .OR. len_io /= head_out(3,citems)) THEN
          WRITE(umMessage,*)'ACUMPS: Error writing partial sum. Code ',     &
            IOSTAT,' on unit ',nftout
          CALL umPrint(umMessage,src='acumps')
          WRITE(umMessage,*)'Words requested ',head_out(3,citems),          &
            ' Words written ',len_io
          CALL umPrint(umMessage,src='acumps')
          icode=7
          cmessage='ACUMPS: BUFFOUT error - see output'
          GO TO 9999
        END IF
        IF (means_total /= 0) THEN
          ! Overwrite field in D1 with partial sum for use by MEANPS
          IF (mype == 0) THEN
            ! Pack and unpack for bit comparison with old system
            IF (head_out(2,citems) == 2) THEN
              DO i=1,head_out(1,citems)
                d1_data(i)=work(i+offset-1)
              END DO
              ! DEPENDS ON: expand32b
              CALL expand32b(global_length,d1_data,fixhd12)
            END IF
          END IF
          ! DEPENDS ON: general_scatter_field
          CALL general_scatter_field(                               &
           d1(address),d1_data,length,global_length,1,              &
           d1_addr(d1_grid_type,ptd1),                              &
           d1_addr(d1_halo_type,ptd1),                              &
           d1_addr(d1_object_type,ptd1),                            &
           d1_addr(d1_proc_no_code,ptd1),                           &
           d1_addr(d1_length,ptd1),                                 &
           d1_addr(d1_no_levels,ptd1),                              &
           d1_addr(d1_north_code,ptd1),                             &
           d1_addr(d1_south_code,ptd1),                             &
           d1_addr(d1_east_code,ptd1),                              &
           d1_addr(d1_west_code,ptd1),                              &
           d1_addr(d1_gridpoint_code,ptd1),                         &
           0,icode,cmessage)
          IF (icode /= 0) GO TO 9999
        END IF
        address=address + length ! Point to next level
      END DO                 ! End loop over levels
    END IF
  END IF                   ! End tagged for meaning
END DO                    ! End of loop over STASH list

9999 CONTINUE
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE acumps
