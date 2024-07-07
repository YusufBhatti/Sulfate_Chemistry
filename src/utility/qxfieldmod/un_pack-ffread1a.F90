! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Interface and arguments: ------------------------------------------
!  Routine: UN_PACK  -------------------------------------------------
!
!  Purpose: To unpack data from the input array FIELD and return
!  the data in FIELD.
!
!  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!
!  Logical components covered: ...
!
!  Project task: ...
!
!  External documentation:
!
!  -------------------------------------------------------------------
!  Interface and arguments: ------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Small execs
SUBROUTINE UN_PACK_ffread1a(pack_type,full_size,field             &
                   ,num_cray_words,ilabel,amdi,pp_fixhd,icode     &
                   ,cmessage)
USE wgdos_packing_mod, ONLY: wgdos_expand_field
USE ereport_mod, ONLY: ereport
USE lookup_addresses
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE errormessagelength_mod, ONLY: errormessagelength
USE Packing_Codes_Mod, ONLY: PC_WGDOS_Packing, PC_RunLength_Packing
IMPLICIT NONE

INTEGER ::                                                        &
     pack_type                                                    &
                          !IN  The type of packing used
    ,full_size                                                    &
                          !IN  The full unpacked size of a field
    ,ilabel(45)                                                   &
                          !OUT holds integer part of LOOKUP
    ,icode                                                        &
                          !OUT Non zero for any error
    ,pp_fixhd(*)          !IN  PPfile fixed length header
REAL ::                                                           &
     field(full_size)                                             &
                          !INOUT On Input contains data.On output
    ,amdi                 !IN  Missing data indicator.
!                               ! contains the un-packed data.
CHARACTER(LEN=errormessagelength) :: cmessage   
                                !OUT Will contain any error mesages.
!
!     LOCAL  VARIABLES
REAL ::                                                           &
     work_array(full_size)  !WORK array used for un_packing
INTEGER ::                                                        &
     int_field(full_size)                                         &
                            ! Integer representation of field array,
                            ! used by wgdos_expand_field
                            ! The length of a FULL_WORD
    ,ixx                                                          &
                            ! Returned X dimension from wgdos_expand_field
    ,iyy                                                          &
                            ! Returned Y dimension from wgdos_expand_field
    ,idum                                                         &
                            ! Dummy variable
    ,num_cray_words                                               &
                            ! IN no of values in an input field
    ,i                                                            &
                            ! Loop counter
    ,num_unpack_values      ! Number of numbers originally packed
!
!
IF (pack_type == PC_WGDOS_Packing) THEN     ! WGDOS packing
  int_field = TRANSFER(field, int_field)
  CALL wgdos_expand_field(work_array,full_size,int_field,         &
                          num_cray_words,idum,ixx,iyy,idum,       &
                          amdi,ilabel(item_code),icode)
  num_unpack_values=ixx*iyy

ELSE IF (pack_type == PC_RunLength_Packing) THEN ! Run length encoded data
  num_unpack_values = ilabel(lbnpt) * ilabel(lbrow)
  ! DEPENDS ON: runlen_decode
  CALL runlen_decode(work_array,full_size,field,num_cray_words,   &
                     amdi,icode,cmessage)
ELSE
  icode=6
  cmessage=' UN_PACK_FFREAD1A - packing type not supported'
END IF
DO i=1,num_unpack_values
  field(i)=work_array(i)
END DO
ilabel(data_type)=1  ! The data type must now be real
ilabel(lbpack)=ilabel(lbpack)-pack_type ! data no longer packed

IF (icode  /=  0) CALL ereport("UN_PACK_FFREAD1A", icode, cmessage)
RETURN
END SUBROUTINE UN_PACK_ffread1a
