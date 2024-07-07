! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Programming Standard: UM DOC Paper3, Verion 4 (05/02/92)
!
!   System Task:C4
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: STASH

!   Interface and arguments: -----------------------------

!      Subroutine: STUFF_INT-----------------------------------------
!
!      Purpose: To put the binary representation of an integer into a
!      real variable through hidden equivalencing via argument passing,
!      when generating extra data for the timeseries.
!

SUBROUTINE stuff_int(array_out,data_in)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE
REAL :: array_out    ! OUT: Real array element in calling routine
REAL :: data_in      ! IN: Integer data in calling routine

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='STUFF_INT'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
array_out=data_in
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE stuff_int
