! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    Routine: CRAY_IBM
!
!    Purpose: To read a direct access PP file  and convert it to a
!    sequential file read to be passed across to the IBM
!
!    Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!
!    Logical components covered: C41
!
!    Project task: C4
!
!    External documentation:
!
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Small execs

SUBROUTINE cray_ibm(npts,num_values,ppunit,                       &
             len1_lookup,pp_len2_lookup,pp_fixhd,lookup,          &
             rookup,entry_no,data_add,model_flag,                 &
             cos_ppunit,iextra,iextraw,last,oper,                 &
             icode,cmessage,lcal360)
USE pp_header_manips, ONLY: header_manip
USE lookup_addresses
USE errormessagelength_mod, ONLY: errormessagelength

USE fort2c_data_conv_interfaces, ONLY:                                         &
  ieee2ibm,                                                                    &
  integer_type,                                                                &
  real_type,                                                                   &
  logical_type

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE um_types, ONLY: real64, integer64

USE packing_codes_mod, ONLY: PC_No_Packing, PC_RunLength_Packing

IMPLICIT NONE

!     Arguments
CHARACTER ::                                                      &
     cmessage*(errormessagelength)           !OUT error messages
LOGICAL ::                                                        &
     last                                                         &
                            !IN indicates last record process
    ,oper                                                         &
                            !IN indicates whether operational
    ,model_flag                                                   &
                            !IN True => dumps, False => fieldsfile
    ,lcal360
INTEGER ::                                                        &
     ppunit                                                       &
                            !IN unit no of required fieldsfile
    ,cos_ppunit                                                   &
                            !IN unit no of COS output file
    ,num_values                                                   &
                            !IN No of data points NROWS*NCOLS
    ,npts                                                         &
                            !IN NUM_VALUES rounded to an even no
!                                 !  used to dimension The output array
          ,data_add                                                     &
                                  !IN The word address of the data.
          ,len1_lookup                                                  &
                                  !IN First dimension of the lookup
          ,pp_len2_lookup                                               &
                                  !IN Size of the LOOKUP on the file
          ,iextra(10)                                                   &
                                  !IN Used within READFF
          ,iextraw                                                      &
                                  !IN no of words of extra data.
          ,entry_no                                                     &
                                  !IN Lookup entry no of the Field.
          ,pp_fixhd(*)                                                  &
                                  !IN PPfile fixed header
          ,lookup(len1_lookup,pp_len2_lookup)                           &
                                               !IN integer lookup
          ,icode                  !OUT error code
REAL ::                                                           &
     rookup(len1_lookup,pp_len2_lookup)    !IN Real lookup
! ---------------------------------------------------------------------
!     arguments for called routines
INTEGER ::                                                        &
     max_len_ilabel                                               &
                       ! maximum length of INT part of pp header
    ,max_len_rlabel    ! maximum length of REAL part of pp header
PARAMETER (max_len_ilabel=45,max_len_rlabel=32)
INTEGER ::                                                        &
     end_year                                                     &
                     ! )
    ,end_month                                                    &
                     ! )
    ,end_day                                                      &
                     ! )  arguments
    ,end_hour                                                     &
                     ! )
    ,end_minute                                                   &
                     ! )     for
    ,end_second                                                   &
                     ! )
    ,end_day_number                                               &
                     ! )
    ,end_time_days                                                &
                     ! )
    ,end_time_secs                                                &
                     ! )  date/time
    ,start_time_secs                                              &
                     ! )
    ,start_time_days                                              &
                     ! )
    ,data_year                                                    &
                     ! )  conversion
    ,data_month                                                   &
                     ! )
    ,data_day                                                     &
                     ! )     when
    ,data_hour                                                    &
                     ! )
    ,data_minute                                                  &
                     ! )  OPER is TRUE
    ,data_second                                                  &
                     ! )
    ,data_day_number                                              &
                     ! )
    ,addr                                                         &
                     ! address in fld, used to process extra data
    ,ibm_addr                                                     &
                     ! address in ibm fld where extra data going.
    ,bit_off                                                      &
                     ! what bit offset are we using
!                            (32 for odd, 0 for even values of addr)
          ,ier                                                          &
                           ! error RETURN CODE from conversion
          ,iv                                                           &
                           ! value of integer code for vectors
          ,len_ilabel                                                   &
                           ! number of values in ILABEL
          ,len_rlabel                                                   &
                           ! number of values in RLABEL
          ,data_values     ! number of values in real extra data


INTEGER(KIND=integer64) ::                                                     &
  ilabel(max_len_ilabel),      & ! holds integer part of LOOKUP
  ibm_label((len1_lookup+1)/2)   ! holds IBM conversion of LABEL

REAL(KIND=real64) ::                                                           &
  field(npts),           & ! array holding data
  ibm_field(npts/2),     & ! array holding IBM data
  rlabel(max_len_rlabel)   ! holds real part of LOOKUP

! ---------------------------------------------------------------------
!    LOCAL VARIABLES
INTEGER ::                                                        &
     i                                                            &
                    ! local counter
    ,pack_type                                                    &
                    ! packing type N1 of LBPACK
    ,number_format                                                &
                    ! number format
    ,fcst_prd

LOGICAL :: packed      ! indicates whether the data is packed



DO i=1,npts        ! make sure FIELD is initialised. An odd
  field(i)=0.0       ! number of points might upset conversion
END DO

!     Initialise output field holding IBM data
DO i=1,npts/2
  ibm_field(i)=0.0
END DO
packed=.FALSE.

!  access the Fields File.
! DEPENDS ON: readff_convpp
CALL READFF_convpp(ppunit,field,npts,entry_no,                    &
ilabel,rlabel,iextra,pp_len2_lookup,len1_lookup,                  &
pp_fixhd,lookup,rookup,data_add,                                  &
model_flag,max_len_ilabel,max_len_rlabel,                         &
len_ilabel,len_rlabel,                                            &
icode,cmessage)
!
IF (icode /= 0) RETURN

!-----------------------------------------------------------------

! The data has now been read in and has 1) Been read in packed
! and left packed or 2) read in as packed and then un-packed or
! 3) The data was never packed at all.  If packed FIELD will have
! LBLREC/2 values if a DUMP and LBLREC values if a PP_FILE. If
! the data is not packed FIELD will have the no of data points
! length LBROW*LBNPT+LBEXT if a pp_file and LBLREC if a dump file.
!
! For a dump LBLREC will hold original no of data points.  For a
! pp_file LBLREC will hold the no of CRAY words needed to hold
! the data (if un-packed also no of data points)
!
! The value returned in ILABEL(LBLREC) may have to change because
! the IBM only has a 32 bit word length compared to the CRAY's 64
! bit word length. On the IBM ILABEL(LBLREC) will be no of IBM
! words needed to hold the data . If the data is not packed (or
! it has been un-packed) then this will be the no of data points.
! If the data is left packed the value of ILABEL(LBLREC) on the
! IBM will have to be doubled as the no of IBM words needed to
! hold the data will twice that on the CRAY.

! On output the data will either have been converted to IBM
! numbers and stored in IBM_FIELD or left packed in FIELD.  If packed
! then LBLREC/2  words of FIELD are written as LBLREC is now
! the no of IBM words. If un-packed IBM_FIELD which has size
! npts/2 (or NUM_VALUES/2) is written as it is.

!-----------------------------------------------------------------
!     decode LBPACK
pack_type = MOD(ilabel(lbpack),10)
number_format = ilabel(lbpack)/1000

!     Run length encoded data needs to be converted into IBM
!     numbers unlike wgdos packed data.
IF (pack_type > PC_No_Packing                                                  &
    .AND. pack_type /= PC_RunLength_Packing) packed=.TRUE.

IF (packed) THEN              ! Data left in packed form. Number of
  ilabel(lblrec)=ilabel(lblrec)*2 ! IBM words needed is 2*CRAY
END IF
! verify that don't have extra data and packing at once
IF (iextraw >  0 .AND. packed) THEN
  cmessage='CONVPP: Extra data with packing not supported'
  icode=1
  RETURN
END IF

CALL header_manip(ilabel)

!     now native format for front-end
ilabel(lbpack) = ilabel(lbpack) -number_format*1000
!
! Note in the header this is in the IBM number format
!ilabel(lbpack) = ilabel(lbpack) -number_format*1000 + 1000

!  Convert ILABEL to IBM(Hitachi) Integers
bit_off = 0
ibm_addr=1

ier = ieee2ibm(integer_type, len_ilabel,                                       &
               ibm_label(ibm_addr:ibm_addr+(len_ilabel+1)/2-1), bit_off,       &
               ilabel(1:len_ilabel), 1, 64, 32)

IF (ier /= 0) THEN
  icode=1
  cmessage=' CRAY_IBM error converting INT for IBM_LABEL'
  RETURN
END IF
!  Convert RLABEL to IBM(Hitachi) Real.
ibm_addr=len_ilabel/2
IF (ibm_addr*2 /= len_ilabel) bit_off=32
ibm_addr=ibm_addr+1

ier = ieee2ibm(real_type, len_rlabel,                                          &
               ibm_label(ibm_addr:ibm_addr+(len_rlabel+1)/2-1), bit_off,       &
               rlabel(1:len_rlabel), 1, 64, 32)

IF (ier /= 0) THEN
  icode=1
  cmessage=' CRAY_IBM error converting REAL for IBM_LABEL'
  RETURN
END IF
bit_off = 0
IF (.NOT. packed) THEN
  !  Convert Real DATA to IBM(Hitachi) Real if not packed.
  IF (ilabel(data_type) == 1) THEN        !Data Type Real
    IF (ilabel(32) == 74) THEN

      ! The output land sea mask has been converted to a real field
      ! in READFF, earlier. (LSM is usually a logical field.)
      ! Add notifcation to CONVPP output file for conversion.

      WRITE(umMessage,'(A,A)') 'Convert type 74 (=landsea mask) from logical', &
                 ' to real. Datatype already labelled as real.'
      CALL umPrint(umMessage,src='cray_ibm')
    END IF
    IF (pack_type == PC_RunLength_Packing) THEN !  Run Length encoding

      ier = ieee2ibm(real_type, ilabel(lblrec),                                &
                     ibm_field(1:(ilabel(lblrec)+1)/2),bit_off,                &
                     field(1:ilabel(lblrec)), 1, 64, 32)

      IF (ier /= 0) THEN
        icode=1
        cmessage='CRAY_IBM error converting real for IBM_FIELD'
        RETURN
      END IF
    ELSE

      ier = ieee2ibm(real_type, num_values-iextraw,                            &
                     ibm_field(1:(num_values-iextraw+1)/2), bit_off,           &
                     field(1:num_values-iextraw), 1, 64, 32)

      IF (ier /= 0) THEN
        icode=1
        cmessage='CRAY_IBM error converting real for IBM_FIELD'
        RETURN
      END IF
    END IF
    !  Convert Integer data to IBM(Hitachi) Integer.
  ELSE IF (ilabel(data_type) == 2) THEN      !Data Type Integer

    ier = ieee2ibm(integer_type, num_values-iextraw,                           &
                   ibm_field(1:(num_values-iextraw+1)/2), bit_off,             &
                   field(1:num_values-iextraw), 1, 64, 32)

    IF (ier /= 0) THEN
      icode=1
      cmessage='CRAY_IBM error converting int for IBM_FIELD'
      RETURN
    END IF
  ELSE IF (ilabel(data_type) == 3) THEN      !Data Type Logical

    ier = ieee2ibm(logical_type, num_values-iextraw,                           &
                   ibm_field(1:(num_values-iextraw+1)/2), bit_off,             &
                   field(1:num_values-iextraw), 1, 64, 32)

    IF (ier /= 0) THEN
      icode=1
      cmessage='CRAY_IBM error converting logical for IBM_FIELD'
      RETURN
    END IF
  END IF
END IF

!  process extra data

              ! About BIT OFFSET
              !
              ! 1         2         3         4         5 (addr)
              ! |---------|---------|---------|---------|  FIELD
              ! .        .         .
              ! .      .       .
              ! .    .    .
              ! |----|----|----|----|----|----|----|----|  IBM_FIELD
              ! 1         2         3         4         5 (ibm_addr)
              !                     |    |
              ! <--------->         |    |
              !  a "word"           | bit_off=32
              !                 bit_off=0
              ! Example:
              !  if ADDR=2, IBM_ADDR=3/2=1
              !  IBM_ADDR*2 eq 1;  so BIT_OFF=32
              !

IF (iextraw >  0) THEN ! process extra data as got some
  !  init values for while loop
  addr=num_values-iextraw+1 ! start address in field for extra dat
  ibm_addr=(addr+1)/2
  IF (ibm_addr*2 == addr) THEN
    bit_off=32
  ELSE
    bit_off=0
  END IF

  DO WHILE (addr <  num_values)
    ! IV is integer header for extra data vector which contains
    ! encoded info for extra data - vector length & data type
    ! Decode IV: data_values will be vector length
    ! Details about extra data, see UM documentation Paper F3
    ! NB. integer header for extra data vector is converted to its
    !     real EQUIVALENCE during model run.  Hence, TRANSFER
    !     serves to convert it back to INTEGER
    iv = TRANSFER(field(addr), iv)

    ! DEPENDS ON: check_extra
    CALL check_extra(iv,data_values,icode,cmessage)
    IF (icode /= 0) THEN
      RETURN
    END IF

    ier = ieee2ibm(integer_type, 1, ibm_field(ibm_addr), bit_off,              &
                   field(addr), 1, 64, 32)

    !         convert the integer from cray format to ibm format
    IF (ier /= 0) THEN
      icode=1
      cmessage='CRAY_IBM: failed in integer conv of extra data'
      RETURN
    END IF

    !          update bit_off, addr and ibm_addr
    IF (bit_off == 0) THEN
      bit_off=32
    ELSE
      bit_off=0
      ibm_addr=ibm_addr+1 ! GONE ON ANOTHER WORD..
    END IF
    addr=addr+1           ! INCREMENT ADDRESS

    ! now to convert REAL vector to IBM format.
    ier = ieee2ibm(real_type, data_values,                                     &
                   ibm_field(ibm_addr:ibm_addr+(data_values+1)/2-1),           &
                   bit_off, field(addr:addr+data_values-1), 1, 64, 32)

    !         convert the real data values
    IF (ier /= 0) THEN
      icode=1
      cmessage='CRAY_IBM: FAILED IN REAL CONV OF EXTRA DATA'
      RETURN
    END IF
    !  update loop variables.
    addr=addr+data_values
    ibm_addr=ibm_addr+data_values/2
    IF ((data_values/2)*2 /= data_values) THEN ! ODD NO. OF VALUES
      IF (bit_off == 0) THEN
        bit_off=32
      ELSE
        bit_off=0
        ibm_addr=ibm_addr+1 ! GONE ON ANOTHER WORD..
      END IF
    END IF
  END DO                   ! continue until run out of data....
  !  Verify addr and ibm_addr have correct values at end of whileloop
  !  first check that addr is ok
  IF (addr /= num_values+1) THEN
    WRITE(cmessage,'(A,I5,A,I5)') 'CRAY_IBM: addr',addr, &
    ' <> num_values+1',num_values+1
    icode=1
    RETURN
  END IF
  !  and so is ibm_addr
  IF (bit_off == 0) ibm_addr=ibm_addr-1
  IF (ibm_addr /= (num_values+1)/2) THEN
    WRITE(cmessage,'(A,I5,A,I5)') 'CRAY_IBM: ibm_addr ', &
    ibm_addr,' <> (num_values+1)/2',(num_values+1)/2
    icode=1
    RETURN
  END IF
END IF ! end processing of extra data

WRITE(cos_ppunit) ibm_label
IF (packed) THEN
  WRITE(cos_ppunit) (field(i),i=1,ilabel(lblrec)/2)
ELSE
  WRITE(cos_ppunit) ibm_field
END IF
!
!   The last field has been processed. An extra field is now written
!   to act as a delimiter for the M08 software. This extra field is
!   a duplicate, but with a PP field code of -99 .
IF (last) THEN
  CALL umPrint( '  WRITING LAST RECORD IN THE COS FILE ',src='cray_ibm')
  ilabel(lbfc)=-99
  !  Convert ILABEL to IBM(Hitachi) Integers
  bit_off = 0
  ibm_addr=1

  ier = ieee2ibm(integer_type, len_ilabel,                                     &
                 ibm_label(ibm_addr:ibm_addr+(len_ilabel+1)/2-1),              &
                 bit_off, ilabel(1:len_ilabel), 1, 64, 32)

  IF (ier /= 0) THEN
    icode=1
    cmessage=' CRAY_IBM error converting INT for IBM_LABEL'
    RETURN
  END IF
  !  Convert RLABEL to IBM(Hitachi) Real.
  ibm_addr=len_ilabel/2
  IF (ibm_addr*2 /= len_ilabel) bit_off=32
  ibm_addr=ibm_addr+1

  ier = ieee2ibm(real_type, len_rlabel,                                        &
                 ibm_label(ibm_addr:ibm_addr+(len_rlabel+1)/2-1), bit_off,     &
                 rlabel(1:len_rlabel), 1, 64, 32)

  IF (ier /= 0) THEN
    icode=1
    cmessage=' CRAY_IBM error converting REAL for IBM_LABEL'
    RETURN
  END IF
  WRITE(cos_ppunit) ibm_label
  IF (packed) THEN
    WRITE(cos_ppunit) (field(i),i=1,ilabel(lblrec)/2)
  ELSE
    WRITE(cos_ppunit) ibm_field
  END IF
END IF
9999 CONTINUE
RETURN
END SUBROUTINE cray_ibm
