! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   -----------------------------------------------------------
!
!   Purpose: Generate extra data for the timeseries.
!   This extra data provides information about what processing was done
!   to produce the timeseries. This information will hopefully be of some
!   use to users doing further processing of the timeseries data.
!   EXTRA_TS_INFO : which generates the codes and sets up the space
!                   for the extra data.
!
!   To some extent this routine has much in common with the
!   multi_spatial routine, but as it has a different function
!   viz generate info on timeseries rather than generating a single time
!   for the timeseries, it is coded separately.
!   When modifying multi_spatial be sure also to modify this routine and
!   vice versa
!
!   Programming Standard: UM DOC Paper3, Verion 4 (05/02/92)
!
!   System Task:C4
!

!   Interface and arguments ------------------------------------
!
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: STASH

SUBROUTINE extra_ts_info(extra_data,extra_data_len,no_records)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

INTEGER :: no_records ! IN how many timeseries records are there ?
INTEGER :: extra_data_len ! IN  size of extra data required
REAL :: extra_data(extra_data_len) ! OUT the extra data array
!   LOCAL PARAMETERS
!   --------------------------------------------------------------------
INTEGER :: no_extra_blocks ! how many blocks of extra data we got ?
PARAMETER(no_extra_blocks=6) ! 6 words to describe
!   Local variables
!   -------------------------------------------------
INTEGER :: record_len ! size of block for extra data
INTEGER :: hdr(no_extra_blocks) ! the headers for each block
! order is lat, long, 2nd lat, 2nd long, first level, 2nd level
DATA hdr/3,4,5,6,7,8/ ! codes for above
INTEGER :: addr ! address in array for writting/reading data
INTEGER :: i ! loop count

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EXTRA_TS_INFO'

! -------------------------------------------------------------------
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
record_len=no_records+1 ! how much info in a block
addr=1
DO i=1,no_extra_blocks ! put the headers into the extra data

  ! We cannot do extra_data(addr)=1000*no_records+hdr(i) as this will
  ! put in a floating point conversion of the integer. We actually want
  ! to save the binary representation of the integer which is done
  ! by STUFF_INT
  ! DEPENDS ON: stuff_int
  CALL stuff_int(extra_data(addr),                                &
    1000*no_records+hdr(i))
  addr=addr+record_len
END DO
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE extra_ts_info
