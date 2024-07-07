! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

SUBROUTINE thin_field(pdata,rdata,pdata_len,ixx,iyy,              &
                      ixxstep,iyystep,full_size,pack_code,amdi,   &
                      lrec, icode, cmessage)
!
!    Subroutine to unpack a field, thin, then repack data.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Small execs
!
!
USE wgdos_packing_mod, ONLY: wgdos_expand_field, wgdos_compress_field
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE errormessagelength_mod, ONLY: errormessagelength
USE packing_codes_mod, ONLY: PC_No_Packing, PC_WGDOS_Packing,                  &
    PC_RunLength_Packing

IMPLICIT NONE

INTEGER :: full_size,pdata_len,pack_code,ixxstep,iyystep
INTEGER :: lrec,icode,fake_stash_code
INTEGER :: pdata(full_size),ixx,iyy,isc
INTEGER :: i,j,k,kk,ix1,iy1
INTEGER :: countx,county
REAL :: rdata(full_size),field(full_size),amdi
CHARACTER(LEN=errormessagelength) :: cmessage

icode = 0
fake_stash_code = 0

IF (pack_code == PC_WGDOS_Packing) THEN
  CALL wgdos_expand_field(field,full_size,pdata,full_size,        &
                          pdata_len,ixx,iyy,isc,amdi,             &
                          fake_stash_code,icode)

  ! If IXX and IYY are not decreased by 1 then GRDSET ( a PP routine)
  ! will fail and give the message 'BAD GRID DEFINITION'.
  ! Unfortunately the same failure occurs if IXX and IYY are decreased
  ! when a step size of 1 is specified so IXX and IYY will only be
  ! decreased for step sizes > 1 (in case anyone uses a step size of 1
  ! instead of using SELECT in the namelist)

  IF (ixxstep >  1) THEN
    ix1 = ixx - 1
  ELSE
    ix1 = ixx
  END IF
  IF (iyystep >  1) THEN
    iy1 = iyy - 1
  ELSE
    iy1 = iyy
  END IF

  county = 0

  k = 1
  DO j=1,iy1,iyystep
    countx = 0
    DO i=1,ix1,ixxstep
      kk = (j-1) * ixx + i
      field(k) = field(kk)
      k = k + 1
      countx = countx + 1
    END DO
    county = county + 1
  END DO

  ixx = (ix1 + ixxstep - 1) / ixxstep
  iyy = (iy1 + iyystep - 1) / iyystep

  CALL wgdos_compress_field(field,full_size,pdata,full_size,      &
                            ixx,pdata_len,isc,amdi,               &
                            fake_stash_code,icode)
  IF (icode /= 0) THEN
    WRITE(umMessage,'(I0,A)') icode,                              &
                         'THIN_FIELD FAILED TRYING TO PACK FIELD:'
    CALL umPrint(umMessage,src='thin_field')
  END IF

ELSE IF (pack_code == PC_RunLength_Packing) THEN

  ! DEPENDS ON: runlen_decode
  CALL runlen_decode(field,ixx*iyy,pdata,lrec,                    &
                     amdi,icode,cmessage )

  IF (ixxstep >  1) THEN
    ix1 = ixx - 1
  ELSE
    ix1 = ixx
  END IF
  IF (iyystep >  1) THEN
    iy1 = iyy - 1
  ELSE
    iy1 = iyy
  END IF

  county = 0

  k = 1
  DO j=1,iy1,iyystep
    countx = 0
    DO i=1,ix1,ixxstep
      kk = (j-1) * ixx + i
      field(k) = field(kk)
      k = k + 1
      countx = countx + 1
    END DO
    county = county + 1
  END DO

  ixx = (ix1 + ixxstep - 1) / ixxstep
  iyy = (iy1 + iyystep - 1) / iyystep

  ! DEPENDS ON: runlen_encode
  CALL runlen_encode(field,ixx*iyy,pdata,ixx*iyy,                 &
                     pdata_len,amdi,icode,cmessage)
ELSE IF (pack_code == PC_No_Packing) THEN

  IF (ixxstep >  1) THEN
    ix1 = ixx - 1
  ELSE
    ix1 = ixx
  END IF
  IF (iyystep >  1) THEN
    iy1 = iyy - 1
  ELSE
    iy1 = iyy
  END IF

  county = 0

  k = 1
  DO j=1,iy1,iyystep
    countx = 0
    DO i=1,ix1,ixxstep
      kk = (j-1) * ixx + i
      rdata(k) = rdata(kk)
      k = k + 1
      countx = countx + 1
    END DO
    county = county + 1
  END DO

  ixx = (ix1 + ixxstep - 1) / ixxstep
  iyy = (iy1 + iyystep - 1) / iyystep
  pdata_len = ixx * iyy

ELSE
  WRITE(umMessage,*)pack_code,' not yet coded'
  CALL umPrint(umMessage,src='thin_field')
END IF

RETURN
END SUBROUTINE thin_field
