! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Reduce storage requirements by removing land points from
!          fields within ocean fieldsfiles. This is done using run
!          length encoding to compress the sequences of missing data
!          values that represent the land points. Using this method,
!          a sequence of missing data values are mapped into a single
!          missing data value followed by a value indicating the
!          length of the run of missing data values.
!
!----------------------------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: STASH

SUBROUTINE runlen_encode(unpacked,unpacked_size,packed,           &
                         packed_size,packed_end,bmdi,             &
                         icode,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


INTEGER :: i,j
INTEGER :: nmdi
INTEGER :: tmdi
INTEGER :: other

INTEGER :: unpacked_size    ! IN Size of input field used to hold
                         !    unpacked data.
INTEGER :: packed_size      ! IN Size of field used to hold packed
                         !    data.
INTEGER :: packed_end       ! OUT Contains the size of the encoded
                         !     field on return.
INTEGER :: icode, ocode     ! IN/OUT Error codes >0 indicates error
                         !                    =0 success.

REAL :: unpacked(unpacked_size) ! IN Input field representing
                         !    unpacked data.
REAL :: packed(packed_size) ! OUT output field representing run
                         !     length decoded fields.
REAL :: bmdi                ! IN Real missing data values

CHARACTER(LEN=errormessagelength) :: cmessage  
                         ! OUT Returned error message if icode >0.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RUNLEN_ENCODE'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
j       = 1
packed_end = 1
nmdi    = 0
tmdi    = 0
other   = 0

!    Check if an error has already been encountered, and get out
!    if it has.
ocode = 0
IF (icode  >   0) THEN
  GO TO 999
ELSE IF (icode  <   0) THEN
  ocode = icode
  icode = 0
END IF

DO i=1,unpacked_size

  IF (unpacked(i)  ==  bmdi) THEN
    nmdi = nmdi + 1
    tmdi = tmdi + 1
  ELSE
    IF (nmdi  >   0) THEN
      IF ((tmdi + other)  <=  unpacked_size) THEN
        packed(j) = bmdi
        packed(j+1) = nmdi
        j = j + 2
        packed_end = packed_end + 2
        nmdi = 0
      ELSE
        WRITE(umMessage,*)                                               &
        'RLENCODE : ERROR : Ocean Run length encoding failed'
        CALL umPrint(umMessage,src='runlen_encode')
        cmessage = 'Ocean Run length encoding failed'
        icode = 1
        GO TO 999
      END IF
    END IF
    packed(j) = unpacked(i)
    packed_end    = packed_end + 1
    other      = other + 1
    j          = j + 1
  END IF
END DO
IF (nmdi  >   0) THEN
  IF ((tmdi + other)  <=  unpacked_size) THEN
    packed(j) = bmdi
    packed(j+1) = nmdi
    j = j + 2
    packed_end = packed_end + 2
  ELSE
    WRITE(umMessage,*)                                                   &
      'RLENCODE : ERROR : Ocean Run length encoding failed'
    CALL umPrint(umMessage,src='runlen_encode')
    cmessage = 'Ocean Run length encoding failed'
    icode = 1
    GO TO 999
  END IF
END IF

999  CONTINUE
!     If we have found an error, leave it in icode.  If no error
!     occurred then check if the original input value of icode was]
!     non-zero (a previous untrapped error/warning), and copy this
!     back into ICODE before eaving the routine.
IF (icode  ==  0 .AND. ocode  /=  0) THEN
  icode = ocode
END IF

packed_end = packed_end - 1

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE runlen_encode



!----------------------------------------------------------------------
