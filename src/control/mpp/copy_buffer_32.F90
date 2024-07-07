! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  A copy from an input 32-bit buffer to an output 32-bit buffer

SUBROUTINE copy_buffer_32( IN, OUT, COUNT, off_out, off_in )


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE um_types
USE fort2c_pio_byteswap_interfaces, ONLY: get_machine_endianism, littleEndian

IMPLICIT NONE
!
! Description:
!   Copies input buffer (32 bit data) to output buffer (32 bit data)
!   possibly with offsets in either input or output.
!
! Method:
!   Very simple copy of input to output with offsets.
!   Only required as a seperate deck to ensure 32 bit sizes are
!   respected (on NEC).
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP
!
! Code Description:
!   Language: FORTRAN 90 + common extensions.
!   This code is written to UMDP3 v8 programming standards.
!
! Declarations:
!
! Global variables (#include statements etc):

! Arguments
INTEGER (KIND=integer64), INTENT(IN)    ::                &
 COUNT                                                    &
         ! size of data to copy
,off_in                                                   &
         ! input offset
,off_out
         ! output offset

INTEGER (KIND=integer32), INTENT(IN) ::                   &
 IN(*)
         ! input data
INTEGER (KIND=integer32), INTENT(INOUT) ::                &
 OUT(*)
         ! output data

! Local variables ----------------------------------------
INTEGER     ::  i        ! Looper

INTEGER     ::                                            &
 offset_in                                                &
            ! looper offset
,offset_out
            ! looper offset

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_BUFFER_32'

! End of header
! --------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

IF (get_machine_endianism()==littleEndian) THEN

! Determine location of first relevant byte in output buffer (offset = +/-1)
offset_out = ( -2 * MOD(off_out, 2) ) + 1

! Determine location of first relevant byte in input buffer (offset = +/-1)
offset_in = ( -2 * MOD(off_in, 2) ) + 1

! Copy alternate 32-bit bytes while shifting +/-1 to reverse byte order
DO i = 1, COUNT - off_in, 2
  OUT(off_out + i + offset_out) = IN(off_in + i + offset_in)
END DO

! Copy remaining 32-bit bytes while shifting -/+1 to reverse byte order
DO i = 2, COUNT - off_in , 2
  OUT(off_out + i - offset_out) = IN(off_in + i - offset_in)
END DO

ELSE

DO i = 1, COUNT - off_in
  OUT(off_out + i) = IN(i + off_in)
END DO

END IF

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE copy_buffer_32
