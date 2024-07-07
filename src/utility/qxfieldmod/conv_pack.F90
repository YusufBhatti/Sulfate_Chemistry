! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Small execs

SUBROUTINE conv_pack(ilabel,rlabel,pack_code,                     &
                     input_pack_type,output_pack_type,            &
                     field,IDIM,len_field,                        &
                     pp_fixhd,icode,cmessage)
USE lookup_addresses
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE errormessagelength_mod, ONLY: errormessagelength

USE packing_codes_mod, ONLY:                                                   &
    PC_No_Packing,                                                             &
    PC_WGDOS_packing,                                                          &
    PC_RunLength_Packing,                                                      &
    PC_GRIB_Packing,                                                           &
    PC_Cray32_Packing

IMPLICIT NONE


INTEGER ::                                                        &
     ilabel(50)                                                   &
    ,pack_code                                                    &
    ,IDIM                                                         &
    ,pp_fixhd(*)                                                  &
    ,len_field                                                    &
    ,icode
REAL ::                                                           &
     rlabel(19)                                                   &
    ,field(IDIM)
CHARACTER ::                                                      &
     input_pack_type*6                                            &
    ,output_pack_type*6                                           &
    ,cmessage*(errormessagelength)
REAL ::                                                           &
     amdi

amdi=rlabel(18)

! DEPENDS ON: un_pack_ffread1a
CALL UN_PACK_ffread1a(pack_code,IDIM,field,len_field,             &
             ilabel,amdi,pp_fixhd,icode,cmessage)

IF (icode /= 0) THEN
  WRITE(7,*) icode
  WRITE(7,*) cmessage
END IF

len_field = ilabel(lbrow) * ilabel(lbnpt)
WRITE(umMessage,*) input_pack_type,' NOW UNPACKED'
CALL umPrint(umMessage,src='conv_pack')

IF (output_pack_type == 'NONE  ') THEN
  pack_code=PC_No_Packing   ! no repacking needed
ELSE IF (output_pack_type == 'WGDOS ') THEN
  pack_code=PC_WGDOS_packing   ! repack using WGDOS
  ! DEPENDS ON: re_pack
  CALL re_pack(pack_code,IDIM,field,len_field,                    &
             ilabel,rlabel,pp_fixhd,icode,cmessage)
ELSE IF (output_pack_type == 'CRAY32') THEN
  pack_code=PC_Cray32_Packing   ! repack using cray 32
  WRITE(umMessage,*) 'packing not supported'
  CALL umPrint(umMessage,src='conv_pack')

ELSE IF (output_pack_type == 'GRIB  ') THEN
  pack_code=PC_GRIB_Packing   ! repack using grib
  ! DEPENDS ON: re_pack
  CALL re_pack(pack_code,IDIM,field,len_field,                    &
             ilabel,rlabel,pp_fixhd,icode,cmessage)
ELSE IF (output_pack_type == 'RUNLEN') THEN
  pack_code=PC_RunLength_Packing   ! repack using run length encoding
  ! DEPENDS ON: re_pack
  CALL re_pack(pack_code,IDIM,field,len_field,                    &
             ilabel,rlabel,pp_fixhd,icode,cmessage)
END IF

IF (icode /= 0) THEN
  WRITE(7,*) icode
  WRITE(7,*) cmessage
END IF
ilabel(lblrec) = len_field

WRITE(umMessage,*) 'NOW PACKED INTO ',output_pack_type
CALL umPrint(umMessage,src='conv_pack')
!     WRITE(6,*)'len_field=',len_field

RETURN
END SUBROUTINE conv_pack
