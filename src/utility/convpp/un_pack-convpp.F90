! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    Suboutine: UN_PACK
!
!    Purpose: To unpack data from the input array FIELD and return
!    the data in FIELD.
!
!    Tested under compiler:   cft77
!    Tested under OS version: UNICOS 5.1
!
!    Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!
!    Logical components covered:
!
!    Project task:
!
!    External documentation:
!
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Small execs
SUBROUTINE UN_PACK_convpp(pack_type,full_size,field,num_cray_words,&
          ilabel,len_ilabel,amdi,pp_fixhd,icode,cmessage)

USE wgdos_packing_mod, ONLY: wgdos_expand_field
USE ereport_mod, ONLY: ereport
USE lookup_addresses

USE errormessagelength_mod, ONLY: errormessagelength
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE packing_codes_mod, ONLY: PC_No_Packing, PC_WGDOS_Packing,                  &
    PC_Cray32_Packing, PC_RunLength_Packing

IMPLICIT NONE

!     arguments
CHARACTER ::                                                      &
     cmessage*(errormessagelength)         !OUT error mesages.
INTEGER ::                                                        &
     pack_type                                                    &
                          !INOUT Type of packing used
    ,full_size                                                    &
                          !IN    full unpacked size of a field
    ,pp_fixhd(*)                                                  &
                          !IN    PPfile fixed length header
    ,num_cray_words                                               &
                          !IN    length of input field
    ,len_ilabel                                                   &
                          !IN    length of ilabel array
    ,ilabel(len_ilabel)                                           &
                          !INOUT holds integer part of LOOKUP
    ,icode                !OUT   Non zero for any error
REAL ::                                                           &
     field(full_size)                                             &
                          !INOUT Input contains packed data.
!                               !      Output contains un-packed data.
          ,amdi                 !IN    Missing data indicator.
! ---------------------------------------------------------------------
!     arguments for called routines
INTEGER ::                                                        &
     ixx                                                          &
                            ! Returned X dimension from wgdos_expand_field
    ,iyy                                                          &
                            ! Returned Y dimension from wgdos_expand_field
    ,idum                   ! Dummy variable
REAL ::                                                           &
     work_array(full_size)  !WORK array used for un_packing

!     LOCAL  VARIABLES
INTEGER ::                                                        &
     int_field(full_size)  & ! Integer representation of field array,
                             ! used by wgdos_expand_field
    ,num_unpack_values                                            &
                            ! Number of numbers originally packed
    ,i                      ! loop counter
!
!
IF (pack_type == PC_WGDOS_Packing) THEN     ! WGDOS packing
  int_field = TRANSFER(field, int_field)
  CALL wgdos_expand_field(work_array,full_size,int_field,num_cray_words, &
                          idum,ixx,iyy,idum,amdi,ilabel(item_code),icode)
  num_unpack_values=ixx*iyy
  ilabel(lblrec)=ilabel(lbrow)*ilabel(lbnpt)+ilabel(lbext)
ELSE IF (pack_type == PC_Cray32_Packing) THEN !  32 Bit CRAY packing
  num_cray_words=num_cray_words*2
  ! DEPENDS ON: expand21
  CALL expand21(num_cray_words,field,work_array)
  num_unpack_values=num_cray_words

ELSE IF (pack_type == PC_RunLength_Packing) THEN !  Run Length encoding
  num_unpack_values = ilabel(lbnpt) * ilabel(lbrow)
  ! DEPENDS ON: runlen_decode
  CALL runlen_decode(work_array,full_size,field,num_cray_words,   &
                     amdi,icode,cmessage)
ELSE
  icode=6
  cmessage=' UNPACK_CONVPP - packing type not supported'
END IF
DO i=1,num_unpack_values
  field(i)=work_array(i)
END DO
ilabel(data_type)=1  ! data must now be real
ilabel(lbpack)=ilabel(lbpack)-pack_type ! data no longer packed
pack_type=PC_No_Packing          ! data now not packed
RETURN
END SUBROUTINE un_pack_convpp
