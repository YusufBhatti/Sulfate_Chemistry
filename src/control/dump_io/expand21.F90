! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  SUBROUTINE EXPAND21
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Dump I/O

!  Purpose: Unpacks IEEE 32-bit data into IEEE 64-bit data.


SUBROUTINE expand21(n, in_array, out_array)
!  expands input array 'in_array' from packed 32-bit into 'out_array' in 64-bit
!
!  n          the number of floating point words to convert
!  in_array   the input array of packed 32-bit numbers
!  out_array  the output array of 64-bit numbers
!

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE um_types, ONLY: integer64, real64, real32

IMPLICIT NONE

! Argument variables
INTEGER (KIND=integer64), INTENT(IN)  :: n
REAL       (KIND=real64), INTENT(IN)  :: in_array ( 1 : (n+1)/2 )
REAL       (KIND=real64), INTENT(OUT) :: out_array( 1 :       n )

! DrHook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EXPAND21'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

!Transfer packed 32-bit values into expanded 64-bit array
out_array( 1 : n ) = TRANSFER( in_array( 1 : (n+1)/2 ),                      &
                               REAL(1.0,KIND=real32), n )

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)

END SUBROUTINE expand21
