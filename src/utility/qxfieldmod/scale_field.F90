! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
SUBROUTINE scale_field(pdata,rdata,npoints,scale_factor,          &
                       lrec,pdata_len,pack_code,amdi,             &
                      plen, icode, cmessage)
!    subroutine to unpack a field, multiply by a scale factor,
!    then repack data.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Small execs
USE wgdos_packing_mod, ONLY: wgdos_expand_field, wgdos_compress_field
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE errormessagelength_mod, ONLY: errormessagelength

USE packing_codes_mod, ONLY:                                                   &
    PC_No_Packing,                                                             &
    PC_WGDOS_packing,                                                          &
    PC_RunLength_Packing

IMPLICIT NONE

INTEGER :: npoints,pdata_len,pack_code
REAL :: field(npoints),rdata(npoints),scale_factor,amdi
INTEGER :: pdata(npoints),nrow,ncol,isc
INTEGER :: lrec, plen
INTEGER :: icode
INTEGER :: fake_stash_code
INTEGER :: i
CHARACTER(LEN=errormessagelength) :: cmessage
icode = 0
fake_stash_code = 0

!     Initialise FIELD variable
DO i=1,npoints
  field(i) = 0.0
END DO
IF (pack_code == PC_WGDOS_packing) THEN
  CALL wgdos_expand_field(field,npoints,pdata,npoints,pdata_len,  &
                          nrow,ncol,isc,amdi,fake_stash_code,     &
                          icode)

  DO i=1,ncol*nrow
    IF (field(i) /= amdi) THEN
      field(i) = field(i) * scale_factor
    END IF
  END DO

  CALL wgdos_compress_field(field,npoints,pdata,npoints,nrow,     &
                            pdata_len,isc,amdi,fake_stash_code,   &
                            icode)

  plen = (pdata_len + 1) /2

ELSE IF (pack_code == PC_RunLength_Packing) THEN

  ! DEPENDS ON: runlen_decode
  CALL runlen_decode(field,npoints,pdata,lrec,                    &
                     amdi,icode,cmessage )
  DO i=1,npoints
    IF (field(i) /= amdi) THEN
      field(i) = field(i) * scale_factor
    END IF
  END DO

  ! DEPENDS ON: runlen_encode
  CALL runlen_encode(field,npoints,pdata,npoints,                 &
                     plen,amdi,icode,cmessage)
ELSE IF (pack_code == PC_No_Packing) THEN
  plen = lrec
  DO i=1,pdata_len
    IF (rdata(i) /= amdi) THEN
      rdata(i) = rdata(i) * scale_factor
    END IF
  END DO
ELSE
  WRITE(umMessage,*)pack_code,' not yet coded'
  CALL umPrint(umMessage,src='scale_field')
END IF

RETURN
END SUBROUTINE scale_field
