! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine interface:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP
SUBROUTINE read_unpack(buf,isize,unpack_size,lookup,fixhd12,      &
expand,                                                          &
icode,cmessage)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE umPrintMgr
USE lookup_addresses
USE missing_data_mod, ONLY: rmdi
USE errormessagelength_mod, ONLY: errormessagelength
USE wgdos_packing_mod, ONLY: wgdos_expand_field
USE um_types, ONLY: real64, integer64

USE packing_codes_mod, ONLY:                                                   &
    PC_No_Packing,                                                             &
    PC_WGDOS_packing,                                                          &
    PC_RunLength_Packing,                                                      &
    PC_Cray32_packing,                                                         &
    PC_GRIB_Packing

IMPLICIT NONE

! Description: Unpacking codes from READFL1A & RDMULT1A is
!              combined into this new deck.  This enables
!              readin data to be unpacked for both UM and SX


! Subroutine Arguments:

INTEGER :: pack_code
INTEGER, INTENT(IN) :: fixhd12      !  IN : 12th element of fixed length header
INTEGER, INTENT(IN) :: isize        !  IN :
INTEGER, INTENT(IN) :: unpack_size  !  IN : Size of data when unpacked

INTEGER, INTENT(OUT) :: icode       !  OUT: return code

INTEGER, INTENT(IN) :: expand

REAL(KIND=real64), INTENT(INOUT) ::  buf(unpack_size)
                                         ! INOUT: Holds field to be unpacked
INTEGER, INTENT(INOUT) :: lookup(64)     ! INOUT: LOOKUP header from dump

CHARACTER(LEN=errormessagelength) ::    cmessage

!local variables
INTEGER ::  i
INTEGER ::  num_unpack_values   ! used to detect error
INTEGER ::  process_size

INTEGER :: dimx,dimy,idum

REAL(KIND=real64) ::  tempbuf(unpack_size) ! Holds field while it's unpacked
                                           ! from buf
REAL ::  probuf(isize),extbuf(isize)
INTEGER(KIND=integer64) ::  int_buf(unpack_size) ! Holds field to be unpacked
                                                 ! as Integer array

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_UNPACK'

! ----------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
pack_code=MOD((lookup(lbpack)),10)

IF (pack_code == PC_Cray32_packing) THEN ! 32 Bit CRAY packing
  IF (lookup(data_type)  ==  1) THEN ! if it's REAL
    IF (expand == 1) THEN
    ! DEPENDS ON: expand32b
    CALL expand32b( lookup(lblrec) , buf, fixhd12 )
    END IF
  ELSE ! not a real field, so not allowed to be compressed
    WRITE(umMessage,*) 'READ_UNPACK : Error trying to uncompress field'
    CALL umPrint(umMessage,src='read_unpack')
    WRITE(umMessage,*) 'This field cannot be compressed/uncompressed ',   &
               'at it is type ',lookup(data_type)
    CALL umPrint(umMessage,src='read_unpack')
    icode=200
    cmessage='Failure uncompressing field'
    GO TO 9999
  END IF  ! if it's REAL

  num_unpack_values=lookup(lblrec)

ELSE IF (pack_code == PC_WGDOS_packing) THEN ! WGDOS packing
  IF (expand == 1) THEN  ! Field is to be expanded
    idum=0

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                  &
!$OMP&            SHARED(isize,buf,int_buf)                 &
!$OMP&            PRIVATE(i)
    DO i=1,isize
      int_buf(i) = TRANSFER(buf(i), int_buf(i))
    END DO
!$OMP END PARALLEL DO

    CALL wgdos_expand_field(buf, unpack_size, int_buf, isize, idum, dimx, &
                            dimy, idum, rmdi, lookup(item_code), icode)

    num_unpack_values=dimx*dimy
  END IF

ELSE IF (pack_code == PC_RunLength_Packing) THEN

  IF (expand == 1) THEN
    num_unpack_values=lookup(lbrow) * lookup(lbnpt)
    ! Enabling expansion of pack_code 4 files with extra data
    IF (lookup(lbext) >  0) THEN

      ! process size not lbrow*lbnpt as data is packed
      process_size=lookup(lblrec)-lookup(lbext)

      ! decompose BUF into grid data array & extra data array
      DO i=1,process_size
        probuf(i)=buf(i)
      END DO
      DO i=1,lookup(lbext)
        extbuf(i)=buf(process_size+i)
      END DO

      ! DEPENDS ON: runlen_decode
      CALL runlen_decode(tempbuf,num_unpack_values,probuf,        &
                         process_size,rmdi,icode,cmessage)

      ! combine expanded grid data with extra data
      DO i=1,num_unpack_values
        buf(i)=tempbuf(i)
      END DO


      DO i=1,lookup(lbext)
        buf(num_unpack_values+i)=extbuf(i)
      END DO

    ELSE   ! if no extra data

      ! DEPENDS ON: runlen_decode
      CALL runlen_decode(tempbuf,num_unpack_values,buf,             &
                          lookup(lblrec),rmdi,icode,cmessage)
    END IF
  END IF


ELSE IF (pack_code /= PC_No_Packing) THEN
  icode=6
  cmessage='READ_UNPACK: packing type not supported'
  IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

IF ( (pack_code /= PC_No_Packing) .AND. (expand == 1) ) THEN
  IF (unpack_size <  num_unpack_values) THEN
    icode=7
    WRITE(umMessage,*)'READ_UNPACK: UNPACK_SIZE = ',unpack_size,          &
      ' but NUM_UNPACK_VALUES = ',num_unpack_values
    CALL umPrint(umMessage,src='read_unpack')
    cmessage='READ_UNPACK: workspace too small for Unpacked Data'
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END IF

! Adjust the value of the data record length if data has been unpacked
IF (((pack_code == PC_WGDOS_packing .OR. pack_code == PC_RunLength_Packing)    &
     .AND. expand == 1) .OR. pack_code == PC_GRIB_Packing) THEN
  IF (lookup(lblrec) /= num_unpack_values) THEN

    IF (lookup(lbext) >  0) THEN
      num_unpack_values=num_unpack_values+lookup(lbext)
    END IF
    IF (printstatus >= prstatus_diag) THEN
      WRITE(umMessage,*)'Record Length Changed from ',                  &
        lookup(lblrec), ' to ',num_unpack_values
      CALL umPrint(umMessage,src='read_unpack')
    END IF

    lookup(lblrec)=num_unpack_values
  END IF
END IF


9999 CONTINUE

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE read_unpack
