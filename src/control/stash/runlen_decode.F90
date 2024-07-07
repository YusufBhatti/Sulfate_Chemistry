! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------------
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


SUBROUTINE runlen_decode(unpacked,unpacked_size,packed,           &
                         packed_size,bmdi,icode,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE


INTEGER :: i,j,k
INTEGER :: nmdi

INTEGER :: unpacked_size    ! IN  Size of input field used to hold
                         !     unpacked data
INTEGER :: packed_size      ! IN  Size of field used to hold packed
                         !     data.
INTEGER :: icode, ocode     ! IN/OUT Error codes >0 indicates error
                         !                    =0 success.

REAL :: unpacked(unpacked_size) ! OUT Input field representing
                         !     unpacked data.
REAL :: packed(packed_size) ! IN  output field representing run
                         !     length decoded fields.
REAL :: bmdi                ! IN  Real missing data values

CHARACTER(LEN=errormessagelength) :: cmessage  
                            ! OUT Returned error message if icode >0

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RUNLEN_DECODE'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
i    = 0
k    = 1
nmdi = 0

!    Check if an error has already been encountered, and get out
!    if it has.
ocode = 0
IF (icode  >   0) THEN
  GO TO 999
ELSE IF (icode  <   0) THEN
  ocode = icode
  icode = 0
END IF

DO WHILE(k <= packed_size)

  IF (packed(k)  ==  bmdi) THEN
    nmdi = packed(k+1)

    DO j=1, nmdi
      i = i + 1
      unpacked(i) = bmdi

      IF (i  >   unpacked_size) THEN
        WRITE(umMessage,*)                                               &
          'RLENCODE : ERROR : Ocean Run length decoding failed'
        CALL umPrint(umMessage,src='runlen_decode')
        cmessage = 'Ocean Run length decoding failed'
        icode = 1
        GO TO 999
      END IF
    END DO

    k = k + 2
  ELSE
    i = i + 1
    unpacked(i) = packed(k)
    k = k + 1

    IF (i  >   unpacked_size) THEN
      WRITE(umMessage,*)                                                 &
        'RLENCODE : ERROR : Ocean Run length decoding failed'
      CALL umPrint(umMessage,src='runlen_decode')
      cmessage = 'Ocean Run length decoding failed'
      icode = 1
      GO TO 999
    END IF

  END IF

END DO

999  CONTINUE

!     If we have found an error, leave it in icode.  If no error
!     occurred then check if the original input value of icode was]
!     non-zero (a previous untrapped error/warning), and copy this
!     back into ICODE before eaving the routine.
IF (icode  ==  0 .AND. ocode  /=  0) THEN
  icode = ocode
END IF

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE runlen_decode


!----------------------------------------------------------------------
