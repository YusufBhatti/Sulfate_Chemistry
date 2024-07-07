! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    Routine: CHECK_EXTRA ----------------------------------------------
!
!    Purpose: To check that code is correct for vector
!
!    Programming standard: UMDP3
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Small execs

SUBROUTINE check_extra(code,data_values,icode,cmessage)
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE

!     arguments
INTEGER, INTENT(IN) :: code
         !IN Code to be checked
INTEGER, INTENT(OUT) :: data_values
         !OUT Number of data values in vector
INTEGER, INTENT(OUT) :: icode 
         !OUT Error code
CHARACTER(LEN=errormessagelength), INTENT(OUT) ::  cmessage
         !OUT Error message if icode > 0 

! Local parameters

! Define valid extra data vector types for use in pp fields.
INTEGER, PARAMETER :: no_extra_vectors=14
INTEGER :: extra_vector(no_extra_vectors) =                       &
          (/ 1,2,3,4,5,6,7,8,9,11,12,13,14,15 /)
! Vector type 10 is not currently supported

! Local variables
LOGICAL :: valid_type ! Flag to indicate valid vector type
INTEGER :: i ! Loop counter
INTEGER :: this_type

cmessage=''
icode=0

data_values=code/1000
this_type=code-data_values*1000
valid_type = .FALSE.
DO i = 1, no_extra_vectors
  IF (this_type == extra_vector(i)) valid_type = .TRUE.
END DO

IF (.NOT. valid_type) THEN
  icode = 1
  cmessage='CHECK_EXTRA: Unsupported extra data vector code'
  WRITE(umMessage,"(A,I0)") "Unsupported Extra Data code", this_type
  CALL umPrint(umMessage,src='check_extra')
END IF

RETURN
END SUBROUTINE check_extra



