! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
SUBROUTINE wind_10m_fix(pdata,rdata,pdata_len,                    &
                        fct,itype,level,iproj,ppunit1,            &
                        wind_10m_scale,wind_10m_orog,             &
                        model_orog,ilabel_orog,rlabel_orog,       &
                        full_size,pack_code,amdi)
!
!    subroutine to unpack a 10m winds and replace if posible by
!    the level 1 wind scaled using wind_10m_scale
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Small execs
!
USE wgdos_packing_mod, ONLY: wgdos_expand_field,wgdos_compress_field
USE lookup_addresses
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE errormessagelength_mod, ONLY: errormessagelength

USE packing_codes_mod, ONLY:                                                   &
    PC_No_Packing,                                                             &
    PC_WGDOS_packing,                                                          &
    PC_RunLength_Packing

IMPLICIT NONE

INTEGER :: full_size,pdata_len,pack_code
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: cmessage
REAL :: rdata(full_size),field(full_size),field1(full_size),amdi
REAL :: model_orog(full_size),rlabel_orog(19),rlabel(19)
REAL :: wind_10m_orog,wind_10m_scale
INTEGER :: pdata(full_size),nrow,ncol,isc
INTEGER :: ilabel_orog(45),ilabel(45)
INTEGER :: fct,itype,itype1,level,level1,iproj,ppunit1
INTEGER :: iextra(10)
INTEGER :: i
INTEGER :: ixx,iyy

DO i=1,10
  iextra(i)=0
END DO

WRITE(umMessage,*) ' read level1 winds'
CALL umPrint(umMessage,src='wind_10m_fix')
itype1 = 6
IF (itype == 75) itype1 = 5
level1 = 1
! DEPENDS ON: ffread
CALL ffread(iproj,fct,itype1,level1,ppunit1,field1,full_size,     &
                ilabel,rlabel,iextra,icode,cmessage)

WRITE(umMessage,*) 'icode=',icode
CALL umPrint(umMessage,src='wind_10m_fix')
IF (icode == 0) THEN
  IF (pack_code == PC_WGDOS_packing) THEN
    WRITE(umMessage,'(A)') 'call wgdos_expand_field'
    CALL umPrint(umMessage,src='wind_10m_fix')
    CALL wgdos_expand_field(field,full_size,pdata,full_size,      &
                            pdata_len,nrow,ncol,isc,amdi,         &
                            ilabel(item_code),icode)
          
    WRITE(umMessage,*)'loop field'
    CALL umPrint(umMessage,src='wind_10m_fix')
    DO i=1,ncol*nrow
      IF (field(i) /= amdi) THEN
        IF (model_orog(i) >= wind_10m_orog) THEN
          WRITE(umMessage,*)i,model_orog(i),field(i),field1(i),field1(i)*.8
          CALL umPrint(umMessage,src='wind_10m_fix')
          field(i) = field1(i) * wind_10m_scale
        END IF
      END IF
    END DO

    CALL wgdos_compress_field(field,full_size,pdata,full_size,    &
                              nrow,pdata_len,isc,amdi,            &
                              ilabel(item_code),icode)
  ELSE IF (pack_code == PC_RunLength_Packing) THEN

    ! DEPENDS ON: runlen_decode
    CALL runlen_decode(field,full_size,pdata,full_size,           &
                       amdi,icode,cmessage )

    DO i=1,pdata_len
      IF (field(i) /= amdi) THEN
        IF (model_orog(i) >= wind_10m_orog) THEN

          field(i) = field1(i) * wind_10m_scale
        END IF
      END IF
    END DO
    ixx=ilabel(18)
    iyy=ilabel(19)
    ! DEPENDS ON: runlen_encode
    CALL runlen_encode(field,ixx*iyy,pdata,ixx*iyy,               &
                       pdata_len,amdi,icode,cmessage)
  ELSE IF (pack_code == PC_No_Packing) THEN
    DO i=1,pdata_len
      IF (rdata(i) /= amdi) THEN
        IF (model_orog(i) >= wind_10m_orog) THEN
          rdata(i) = field1(i) * wind_10m_scale
        END IF
      END IF
    END DO
  ELSE
    WRITE(umMessage,*)pack_code,' not yet coded'
    CALL umPrint(umMessage,src='wind_10m_fix')
  END IF
ELSE
  WRITE(7,*) icode
  WRITE(7,*) '10M WIND FAILED WITH THE FOLLOWING REASON:'
  WRITE(7,*) cmessage
END IF

RETURN
END SUBROUTINE wind_10m_fix
