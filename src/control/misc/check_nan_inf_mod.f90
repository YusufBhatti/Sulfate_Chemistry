! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved. 
! For further details please refer to the file COPYRIGHT.txt 
! which you should have received as part of this distribution. 
! *****************************COPYRIGHT*******************************

MODULE check_nan_inf_mod

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook
IMPLICIT NONE

! Description:
! Subroutine to fail if array contains NaNs or Infinites

! Method:
! Check contents of array using um_has_nan and um_has_inf and if either
! gives a positive result the fail using EReport

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Control/Misc

! Language: Fortran 2003.
! This code is written to UMDP3 standards.


PRIVATE
PUBLIC :: fail_if_nan_inf
INTERFACE fail_if_nan_inf
  MODULE PROCEDURE fail_if_nan_inf64,    fail_if_nan_inf32, &
                   fail_if_nan_inf64_2d, fail_if_nan_inf32_2d, &
                   fail_if_nan_inf64_3d, fail_if_nan_inf32_3d, &
                   fail_if_nan_inf64_4d, fail_if_nan_inf32_4d, &
                   fail_if_nan_inf64_5d, fail_if_nan_inf32_5d
END INTERFACE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
CHARACTER(LEN=*), PARAMETER :: ModuleName='CHECK_NAN_INF_MOD'

CONTAINS

  FUNCTION nan_inf_error_msg(varname, has_nan, has_inf)

  USE errormessagelength_mod, ONLY: errormessagelength
  USE umPrintMgr, ONLY: newline

  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN) :: varname
  LOGICAL, INTENT(IN):: has_nan, has_inf

  CHARACTER(LEN=errormessagelength) :: nan_inf_error_msg   ! Function result
  CHARACTER(LEN=errormessagelength) :: nanmsg, infmsg

  IF (has_nan) THEN
    WRITE(nanmsg, '(A,A)') 'Found NaN in ', varname
  ELSE
    WRITE(nanmsg, '(A,A)') 'No NaN in ', varname
  END IF

  IF (has_inf) THEN
    WRITE(infmsg, '(A,A)') 'Found Inf in ', varname
  ELSE
    WRITE(infmsg, '(A,A)') 'No Inf in ', varname
  END IF

  WRITE(nan_inf_error_msg, '(A)') TRIM(nanmsg) // newline // TRIM(infmsg)

END FUNCTION nan_inf_error_msg

!***************************************************************************
! 1D Array 64-bit version
!***************************************************************************
SUBROUTINE fail_if_nan_inf64(x, varname, warn)
USE um_is_nan_mod, ONLY: um_has_nan
USE um_is_inf_mod, ONLY: um_has_inf
USE um_types, ONLY : real64
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Function argument
REAL (KIND=real64), INTENT(IN) :: x(:)
CHARACTER(LEN=*), INTENT(IN) :: varname
LOGICAL, OPTIONAL :: warn

! Local data
LOGICAL                           :: has_nan, has_inf
REAL(KIND=jprb)                   :: zhook_handle
CHARACTER(LEN=*), PARAMETER       :: RoutineName='FAIL_IF_NAN_INF64'
CHARACTER(LEN=errormessagelength) :: Cmessage      ! used for EReport
INTEGER                           :: ErrorStatus   ! used for EReport

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

has_nan = um_has_nan(x)
has_inf = um_has_inf(x)
IF (has_nan .OR. has_inf) THEN
  Cmessage = nan_inf_error_msg(varname, has_nan, has_inf)
  ErrorStatus = 64
  IF (PRESENT(warn)) THEN
    IF (warn) ErrorStatus = -64
  END IF
  CALL EReport( RoutineName, ErrorStatus, Cmessage )
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE fail_if_nan_inf64

!***************************************************************************
! 1D Array 32-bit version
!***************************************************************************
SUBROUTINE fail_if_nan_inf32(x, varname, warn)
USE um_is_nan_mod, ONLY: um_has_nan
USE um_is_inf_mod, ONLY: um_has_inf
USE um_types, ONLY : real32
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Function argument
REAL (KIND=real32), INTENT(IN) :: x(:)
CHARACTER(LEN=*), INTENT(IN) :: varname
LOGICAL, OPTIONAL :: warn

! Local data
LOGICAL                           :: has_nan, has_inf
REAL(KIND=jprb)                   :: zhook_handle
CHARACTER(LEN=*), PARAMETER       :: RoutineName='FAIL_IF_NAN_INF32'
CHARACTER(LEN=errormessagelength) :: Cmessage      ! used for EReport
INTEGER                           :: ErrorStatus   ! used for EReport

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

has_nan = um_has_nan(x)
has_inf = um_has_inf(x)
IF (has_nan .OR. has_inf) THEN
  Cmessage = nan_inf_error_msg(varname, has_nan, has_inf)
  ErrorStatus = 32
  IF (PRESENT(warn)) THEN
    IF (warn) ErrorStatus = -32
  END IF
  CALL EReport( RoutineName, ErrorStatus, Cmessage )
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE fail_if_nan_inf32

!***************************************************************************
! 2D Array 64-bit version
!***************************************************************************
SUBROUTINE fail_if_nan_inf64_2d(x, varname, warn)
USE um_types, ONLY : real64

IMPLICIT NONE

! Function argument
REAL (KIND=real64), INTENT(IN) :: x(:,:)
CHARACTER(LEN=*), INTENT(IN) :: varname
LOGICAL, OPTIONAL :: warn

! Local data
LOGICAL                     :: lwarn
REAL(KIND=jprb)             :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='FAIL_IF_NAN_INF64_2D'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Sort out optional input
IF (PRESENT(warn)) THEN
  lwarn = warn
ELSE
  lwarn = .FALSE.
END IF

! Reshape array and pass through to 1d array version
CALL fail_if_nan_inf64(RESHAPE(x, (/SIZE(x)/)), varname, lwarn)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE fail_if_nan_inf64_2d

!***************************************************************************
! 2D Array 32-bit version
!***************************************************************************
SUBROUTINE fail_if_nan_inf32_2d(x, varname, warn)
USE um_types, ONLY : real32

IMPLICIT NONE

! Function argument
REAL (KIND=real32), INTENT(IN) :: x(:,:)
CHARACTER(LEN=*), INTENT(IN) :: varname
LOGICAL, OPTIONAL :: warn

! Local data
LOGICAL                     :: lwarn
REAL(KIND=jprb)             :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='FAIL_IF_NAN_INF32_2D'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Sort out optional input
IF (PRESENT(warn)) THEN
  lwarn = warn
ELSE
  lwarn = .FALSE.
END IF

! Reshape array and pass through to 1d array version
CALL fail_if_nan_inf32(RESHAPE(x, (/SIZE(x)/)), varname, lwarn)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE fail_if_nan_inf32_2d

!***************************************************************************
! 3D Array 64-bit version
!***************************************************************************
SUBROUTINE fail_if_nan_inf64_3d(x, varname, warn)
USE um_types, ONLY : real64

IMPLICIT NONE

! Function argument
REAL (KIND=real64), INTENT(IN) :: x(:,:,:)
CHARACTER(LEN=*), INTENT(IN) :: varname
LOGICAL, OPTIONAL :: warn

! Local data
LOGICAL                     :: lwarn
REAL(KIND=jprb)             :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='FAIL_IF_NAN_INF64_3D'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Sort out optional input
IF (PRESENT(warn)) THEN
  lwarn = warn
ELSE
  lwarn = .FALSE.
END IF

! Reshape array and pass through to 1d array version
CALL fail_if_nan_inf64(RESHAPE(x, (/SIZE(x)/)), varname, lwarn)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE fail_if_nan_inf64_3d

!***************************************************************************
! 3D Array 32-bit version
!***************************************************************************
SUBROUTINE fail_if_nan_inf32_3d(x, varname, warn)
USE um_types, ONLY : real32

IMPLICIT NONE

! Function argument
REAL (KIND=real32), INTENT(IN) :: x(:,:,:)
CHARACTER(LEN=*), INTENT(IN) :: varname
LOGICAL, OPTIONAL :: warn

! Local data
LOGICAL                     :: lwarn
REAL(KIND=jprb)             :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='FAIL_IF_NAN_INF32_3D'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Sort out optional input
IF (PRESENT(warn)) THEN
  lwarn = warn
ELSE
  lwarn = .FALSE.
END IF

! Reshape array and pass through to 1d array version
CALL fail_if_nan_inf32(RESHAPE(x, (/SIZE(x)/)), varname, lwarn)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE fail_if_nan_inf32_3d

!***************************************************************************
! 4D Array 64-bit version
!***************************************************************************
SUBROUTINE fail_if_nan_inf64_4d(x, varname, warn)
USE um_types, ONLY : real64

IMPLICIT NONE

! Function argument
REAL (KIND=real64), INTENT(IN) :: x(:,:,:,:)
CHARACTER(LEN=*), INTENT(IN) :: varname
LOGICAL, OPTIONAL :: warn

! Local data
LOGICAL                     :: lwarn
REAL(KIND=jprb)             :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='FAIL_IF_NAN_INF64_4D'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Sort out optional input
IF (PRESENT(warn)) THEN
  lwarn = warn
ELSE
  lwarn = .FALSE.
END IF

! Reshape array and pass through to 1d array version
CALL fail_if_nan_inf64(RESHAPE(x, (/SIZE(x)/)), varname, lwarn)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE fail_if_nan_inf64_4d

!***************************************************************************
! 4D Array 32-bit version
!***************************************************************************
SUBROUTINE fail_if_nan_inf32_4d(x, varname, warn)
USE um_types, ONLY : real32

IMPLICIT NONE

! Function argument
REAL (KIND=real32), INTENT(IN) :: x(:,:,:,:)
CHARACTER(LEN=*), INTENT(IN) :: varname
LOGICAL, OPTIONAL :: warn

! Local data
LOGICAL                     :: lwarn
REAL(KIND=jprb)             :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='FAIL_IF_NAN_INF32_4D'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Sort out optional input
IF (PRESENT(warn)) THEN
  lwarn = warn
ELSE
  lwarn = .FALSE.
END IF

! Reshape array and pass through to 1d array version
CALL fail_if_nan_inf32(RESHAPE(x, (/SIZE(x)/)), varname, lwarn)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE fail_if_nan_inf32_4d

!***************************************************************************
! 5D Array 64-bit version
!***************************************************************************
SUBROUTINE fail_if_nan_inf64_5d(x, varname, warn)
USE um_types, ONLY : real64

IMPLICIT NONE

! Function argument
REAL (KIND=real64), INTENT(IN) :: x(:,:,:,:,:)
CHARACTER(LEN=*), INTENT(IN) :: varname
LOGICAL, OPTIONAL :: warn

! Local data
LOGICAL                     :: lwarn
REAL(KIND=jprb)             :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='FAIL_IF_NAN_INF64_5D'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Sort out optional input
IF (PRESENT(warn)) THEN
  lwarn = warn
ELSE
  lwarn = .FALSE.
END IF

! Reshape array and pass through to 1d array version
CALL fail_if_nan_inf64(RESHAPE(x, (/SIZE(x)/)), varname, lwarn)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE fail_if_nan_inf64_5d

!***************************************************************************
! 5D Array 32-bit version
!***************************************************************************
SUBROUTINE fail_if_nan_inf32_5d(x, varname, warn)
USE um_types, ONLY : real32

IMPLICIT NONE

! Function argument
REAL (KIND=real32), INTENT(IN) :: x(:,:,:,:,:)
CHARACTER(LEN=*), INTENT(IN) :: varname
LOGICAL, OPTIONAL :: warn

! Local data
LOGICAL                     :: lwarn
REAL(KIND=jprb)             :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='FAIL_IF_NAN_INF32_5D'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Sort out optional input
IF (PRESENT(warn)) THEN
  lwarn = warn
ELSE
  lwarn = .FALSE.
END IF

! Reshape array and pass through to 1d array version
CALL fail_if_nan_inf32(RESHAPE(x, (/SIZE(x)/)), varname, lwarn)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE fail_if_nan_inf32_5d

END MODULE check_nan_inf_mod
