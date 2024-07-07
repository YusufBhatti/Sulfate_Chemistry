! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Serial UM interface to BUFFIN

MODULE read_serial_mod

IMPLICIT NONE

PRIVATE
PUBLIC :: read_serial

CHARACTER(LEN=*), PARAMETER :: ModuleName='READ_SERIAL_MOD'

CONTAINS

! Subroutine Interface:
SUBROUTINE read_serial(nft, d1, isize, unpack_size, len_io, lookup, fixhd12,   &
                       imodl, section, item, icode, cmessage, expand)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE io

USE lookup_addresses, ONLY: lbpack

USE umPrintMgr, ONLY:                                                          &
    umPrint,                                                                   &
    umMessage

USE errormessagelength_mod, ONLY: errormessagelength

USE Packing_Codes_Mod, ONLY: PC_Cray32_Packing

IMPLICIT NONE

! Description:
!  This routine provides a serial interface to BUFFIN for the
!  Unified Model. It is used where only one PE must read in
!  the whole of a field.

! Method:
!  Reads in the field from file (via BUFFIN), then unpacks

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP

! Subroutine Arguments:

INTEGER, INTENT(IN) ::  nft         !  FORTRAN unit number
INTEGER, INTENT(IN) ::  isize       !  no. of words to be read in
                                    !  (global field)
INTEGER, INTENT(IN) ::  unpack_size !  no. of words after any unpacking (global)
INTEGER, INTENT(OUT) :: len_io      !  no. of words read in (global field)

INTEGER, INTENT(IN) :: expand      ! expansion for small exec

INTEGER, INTENT(IN) :: fixhd12     ! 12th element of fixed length header
INTEGER, INTENT(OUT) :: icode      ! Return code


INTEGER, INTENT(IN) :: imodl        ! model code of the field to be read in
INTEGER, INTENT(IN) :: section      ! section code of the field to be read in
INTEGER, INTENT(IN) :: item         ! item code of the field to be read in

INTEGER, INTENT(INOUT) :: lookup(64)               !   LOOKUP header from dump

CHARACTER(LEN=errormessagelength), INTENT(OUT) ::  cmessage !  Error message

REAL, INTENT(OUT) ::  d1(*)          !  Array to read data in to

! Local variables
INTEGER :: i             ! loop counter

REAL    :: buf_icode    ! return code from BUFFIN

LOGICAL :: broadcast_read_flag ! original file broadcast behaviour

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_SERIAL'

! ------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! IO may be set to broadcast data read in for this unit to other PEs.
! If we are reading serially, we need to temporarily suspend this behaviour

! Check current broadcast behaviour for this unit
broadcast_read_flag = broadcast_read(nft)

! Switch off IO broadcasting
IF (broadcast_read_flag) CALL set_unit_bcast_flag(nft)

IF (MOD(lookup(lbpack),10)  ==  PC_Cray32_Packing) THEN

  ! Data is packed using CRAY 32 bit method - note that we need to read
  ! in 2*ISIZE 32 bit words using BUFFIN32_F77 (as the array is 64 bit)
  ! DEPENDS ON: buffin32_f77
  CALL buffin32_f77(nft,d1,2*isize,len_io,buf_icode)
  ! And then halve LEN_IO to satisfy tests against ISIZE

  len_io = len_io/2

ELSE
  ! Data is not 32-bit packed
  CALL buffin(nft,d1(1:isize),isize,len_io,buf_icode)
END IF

! Restore original broadcast behaviour
IF (broadcast_read_flag) CALL clear_unit_bcast_flag(nft)

! Has the data been read in OK?
IF ((buf_icode  /=  -1.0) .OR. (len_io  /=  isize)) THEN
  WRITE(umMessage,'(A)') 'READ_SERIAL : Error in call to BUFFIN'
  CALL umPrint(umMessage,src='read_serial')
  WRITE(umMessage,'(A,I0)') 'LEN_IO : ',len_io
  CALL umPrint(umMessage,src='read_serial')
  WRITE(umMessage,'(A,E15.7)') 'IOSTAT : ',buf_icode
  CALL umPrint(umMessage,src='read_serial')
  WRITE(umMessage,'(A,I0,1X,I0,1X,I0)')                                        &
        'Attempting to read field (Model,Section,Item) ',                      &
        imodl, section, item
  CALL umPrint(umMessage,src='read_serial')
  icode=100
  cmessage='Failure reading in field'
  GO TO 9999
END IF

! ----------------------------------------------
! If it's a compressed REAL field, expand it out
! ----------------------------------------------

! DEPENDS ON: read_unpack
CALL read_unpack(d1, isize, unpack_size, lookup, fixhd12, expand, icode,       &
                 cmessage)

IF (icode /= 0) THEN
  cmessage='Failure unpacking field'
  GO TO 9999
END IF

9999  CONTINUE

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE read_serial

END MODULE read_serial_mod
