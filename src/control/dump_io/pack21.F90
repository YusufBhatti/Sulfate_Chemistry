! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  SUBROUTINE PACK21
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Dump I/O

!  Purpose: Packs IEEE 64-bit data into IEEE 32-bit data.


SUBROUTINE pack21(n, IN, OUT)
!--compresses input array 'in' from 64-bit into 'out' in 32-bit

!  n       the number of floating point words to convert
!  in      the input array of 64-bit numbers
!  out     the output array of 32-bit numbers

!--

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE um_types
IMPLICIT NONE


! Argument variables
INTEGER (KIND=integer64), INTENT(IN) :: n
REAL (KIND=real64), INTENT(IN)       :: IN(1:n)
REAL (KIND=real32), INTENT(OUT)      :: OUT(1:n)

! Local real parameter
REAL, PARAMETER :: tiny32=TINY(OUT(1))
! Local variables
INTEGER (KIND=integer64) :: i

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PACK21'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
DO i = 1,n
  IF (ABS(IN(i))  <   tiny32) THEN
    ! Prevent 'denormalized' numbers
    OUT(i) = 0.0
  ELSE
    OUT(i) = IN(i)
  END IF
END DO

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE pack21
