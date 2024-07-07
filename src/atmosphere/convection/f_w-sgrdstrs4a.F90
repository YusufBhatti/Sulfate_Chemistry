! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Returns an ensemble vertical velocity profile for shallow cumulus
!
#if !defined(NEC)
REAL FUNCTION f_w(x)   ! ENSEMBLE VERTICAL VELOCITY PROFILE

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE
!
! Description:
!   Calculates non-dimensional vertical velocity profile for shallow cumulus
!
! Method:
!   Derived from large-eddy simulations.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.
!-----------------------------------------------------------------------------

REAL, INTENT(IN) :: x             ! non-dimesional height in cumulus layer.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='F_W'


!-----------------------------------------------------------------------------
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

IF (x <= 1.0) THEN
  f_w=6.0*x
ELSE
  f_w=6.0
END IF

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END FUNCTION f_w
#endif
