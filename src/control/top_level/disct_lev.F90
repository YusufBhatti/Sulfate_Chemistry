! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Test whether level type is discrete (model) or continuous (non-model)
! Function Interface:
LOGICAL FUNCTION disct_lev(lev_code,ErrorStatus,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!  System task:               Sub-Models Project
!

! Function arguments:
!   Scalar arguments with intent(in):
INTEGER :: lev_code !Level code from STASHmaster

! ErrorStatus
INTEGER :: ErrorStatus
CHARACTER(LEN=errormessagelength) :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DISCT_LEV'

!- End of Header ----------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
IF (lev_code == 0 .OR. lev_code == 1 .OR. lev_code == 2 .OR.      &
    lev_code == 6 .OR. lev_code == 10) THEN
  disct_lev=.TRUE.
ELSE IF (lev_code  >=  0 .AND. lev_code  <=  10) THEN
  disct_lev=.FALSE.
ELSE
  disct_lev=.FALSE.
  ErrorStatus=1
  cmessage='DISCT_LEV : Invalid level type in STASHmaster'
END IF
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END FUNCTION disct_lev
!- End of Function code --------------------------------------------
