! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved. 
! For further details please refer to the file COPYRIGHT.txt 
! which you should have received as part of this distribution. 
! *****************************COPYRIGHT*******************************

MODULE um_is_nan_mod

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook
IMPLICIT NONE

! Description:
! A simple module to wrap, or to provide equivalent functionality for,
! the ieee_is_nan function. This allows a fully portable
! NaN test to be performed in both 32 and 64 bit reals.

! Method:
! Use the ieee_is_nan function if available, otherwise test value for
! equality against itself (as NaNs can't equal anything)

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Control/Misc

! Language: Fortran 2003.
! This code is written to UMDP3 standards.


PRIVATE
PUBLIC :: um_is_nan, um_has_nan
INTERFACE um_is_nan
  MODULE PROCEDURE um_is_nan64, um_is_nan32
END INTERFACE
INTERFACE um_has_nan
  MODULE PROCEDURE um_has_nan64,    um_has_nan32, &
                   um_has_nan64_2d, um_has_nan32_2d, &
                   um_has_nan64_3d, um_has_nan32_3d, &
                   um_has_nan64_4d, um_has_nan32_4d, &
                   um_has_nan64_5d, um_has_nan32_5d
END INTERFACE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
CHARACTER(LEN=*), PARAMETER :: ModuleName='UM_IS_NAN_MOD'

CONTAINS

!***************************************************************************
! Scalar 64-bit version
!***************************************************************************
LOGICAL FUNCTION um_is_nan64(x)
#if !(defined(GNU_FORTRAN) && GNU_FORTRAN<4010000)
! gfortran only supports ieee_aritmetic from 4.10
USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_nan, ieee_support_nan
#endif
USE um_types, ONLY : real64

IMPLICIT NONE

! Function argument
REAL (KIND=real64), INTENT(IN) :: x

! Local data
REAL(KIND=jprb)             :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='UM_IS_NAN64'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

#if !(defined(GNU_FORTRAN) && GNU_FORTRAN<4010000)
! Use the ieee version if supported
IF (ieee_support_nan(x)) THEN
  um_is_nan64 = ieee_is_nan(x)
ELSE
#endif
  ! Otherwise we cook up our own version
  IF (x /= x) THEN
    um_is_nan64 = .TRUE.
  ELSE
    um_is_nan64 = .FALSE.
  END IF
#if !(defined(GNU_FORTRAN) && GNU_FORTRAN<4010000)
END IF
#endif

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION um_is_nan64

!***************************************************************************
! Scalar 32-bit version
!***************************************************************************
LOGICAL FUNCTION um_is_nan32(x)
#if !(defined(GNU_FORTRAN) && GNU_FORTRAN<4010000)
! gfortran only supports ieee_aritmetic from 4.10
USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_nan, ieee_support_nan
#endif
USE um_types, ONLY : real32

IMPLICIT NONE

! Function argument
REAL (KIND=real32), INTENT(IN) :: x

! Local data
REAL(KIND=jprb)             :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='UM_IS_NAN32'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

#if !(defined(GNU_FORTRAN) && GNU_FORTRAN<4010000)
! Use the ieee version if supported
IF (ieee_support_nan(x)) THEN
  um_is_nan32 = ieee_is_nan(x)
ELSE
#endif
  ! Otherwise we cook up our own version
  IF (x /= x) THEN
    um_is_nan32 = .TRUE.
  ELSE
    um_is_nan32 = .FALSE.
  END IF
#if !(defined(GNU_FORTRAN) && GNU_FORTRAN<4010000)
END IF
#endif

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION um_is_nan32

!***************************************************************************
! 1D Array 64-bit version
!***************************************************************************
LOGICAL FUNCTION um_has_nan64(x)
USE um_types, ONLY : real64

IMPLICIT NONE

! Function argument
REAL (KIND=real64), INTENT(IN) :: x(:)

! Local data
INTEGER                     :: ix
REAL(KIND=jprb)             :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='UM_HAS_NAN64'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Loop over elements of x and determine if any are NaNs
! Exit immediately if any are found
DO ix=1,SIZE(x)
  um_has_nan64 = um_is_nan64(x(ix))
  IF (um_has_nan64) EXIT
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION um_has_nan64

!***************************************************************************
! 1D Array 32-bit version
!***************************************************************************
LOGICAL FUNCTION um_has_nan32(x)
USE um_types, ONLY : real32

IMPLICIT NONE

! Function argument
REAL (KIND=real32), INTENT(IN) :: x(:)

! Local data
INTEGER                     :: ix
REAL(KIND=jprb)             :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='UM_HAS_NAN32'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Loop over elements of x and determine if any are NaNs
! Exit immediately if any are found
DO ix=1,SIZE(x)
  um_has_nan32 = um_is_nan32(x(ix))
  IF (um_has_nan32) EXIT
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION um_has_nan32

! To use for multi-dimensional arrays you can call um_has_nan with the array
! reshaped, e.g. um_has_nan(RESHAPE(x, (/SIZE(x)/)))

!***************************************************************************
! 2D Array 64-bit version
!***************************************************************************
LOGICAL FUNCTION um_has_nan64_2d(x)
USE um_types, ONLY : real64

IMPLICIT NONE

! Function argument
REAL (KIND=real64), INTENT(IN) :: x(:,:)

! Local data
REAL(KIND=jprb)             :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='UM_HAS_NAN64_2D'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Reshape array and pass through 1d array version
um_has_nan64_2d = um_has_nan64(RESHAPE(x, (/SIZE(x)/)))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION um_has_nan64_2d

!***************************************************************************
! 2D Array 32-bit version
!***************************************************************************
LOGICAL FUNCTION um_has_nan32_2d(x)
USE um_types, ONLY : real32

IMPLICIT NONE

! Function argument
REAL (KIND=real32), INTENT(IN) :: x(:,:)

! Local data
REAL(KIND=jprb)             :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='UM_HAS_NAN32_2D'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Reshape array and pass through 1d array version
um_has_nan32_2d = um_has_nan32(RESHAPE(x, (/SIZE(x)/)))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION um_has_nan32_2d

!***************************************************************************
! 3D Array 64-bit version
!***************************************************************************
LOGICAL FUNCTION um_has_nan64_3d(x)
USE um_types, ONLY : real64

IMPLICIT NONE

! Function argument
REAL (KIND=real64), INTENT(IN) :: x(:,:,:)

! Local data
REAL(KIND=jprb)             :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='UM_HAS_NAN64_3D'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Reshape array and pass through 1d array version
um_has_nan64_3d = um_has_nan64(RESHAPE(x, (/SIZE(x)/)))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION um_has_nan64_3d

!***************************************************************************
! 3D Array 32-bit version
!***************************************************************************
LOGICAL FUNCTION um_has_nan32_3d(x)
USE um_types, ONLY : real32

IMPLICIT NONE

! Function argument
REAL (KIND=real32), INTENT(IN) :: x(:,:,:)

! Local data
REAL(KIND=jprb)             :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='UM_HAS_NAN32_3D'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Reshape array and pass through 1d array version
um_has_nan32_3d = um_has_nan32(RESHAPE(x, (/SIZE(x)/)))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION um_has_nan32_3d

!***************************************************************************
! 4D Array 64-bit version
!***************************************************************************
LOGICAL FUNCTION um_has_nan64_4d(x)
USE um_types, ONLY : real64

IMPLICIT NONE

! Function argument
REAL (KIND=real64), INTENT(IN) :: x(:,:,:,:)

! Local data
REAL(KIND=jprb)             :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='UM_HAS_NAN64_4D'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Reshape array and pass through 1d array version
um_has_nan64_4d = um_has_nan64(RESHAPE(x, (/SIZE(x)/)))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION um_has_nan64_4d

!***************************************************************************
! 4D Array 32-bit version
!***************************************************************************
LOGICAL FUNCTION um_has_nan32_4d(x)
USE um_types, ONLY : real32

IMPLICIT NONE

! Function argument
REAL (KIND=real32), INTENT(IN) :: x(:,:,:,:)

! Local data
REAL(KIND=jprb)             :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='UM_HAS_NAN32_4D'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Reshape array and pass through 1d array version
um_has_nan32_4d = um_has_nan32(RESHAPE(x, (/SIZE(x)/)))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION um_has_nan32_4d

!***************************************************************************
! 5D Array 64-bit version
!***************************************************************************
LOGICAL FUNCTION um_has_nan64_5d(x)
USE um_types, ONLY : real64

IMPLICIT NONE

! Function argument
REAL (KIND=real64), INTENT(IN) :: x(:,:,:,:,:)

! Local data
REAL(KIND=jprb)             :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='UM_HAS_NAN64_5D'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Reshape array and pass through 1d array version
um_has_nan64_5d = um_has_nan64(RESHAPE(x, (/SIZE(x)/)))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION um_has_nan64_5d

!***************************************************************************
! 5D Array 32-bit version
!***************************************************************************
LOGICAL FUNCTION um_has_nan32_5d(x)
USE um_types, ONLY : real32

IMPLICIT NONE

! Function argument
REAL (KIND=real32), INTENT(IN) :: x(:,:,:,:,:)

! Local data
REAL(KIND=jprb)             :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='UM_HAS_NAN32_5D'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Reshape array and pass through 1d array version
um_has_nan32_5d = um_has_nan32(RESHAPE(x, (/SIZE(x)/)))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION um_has_nan32_5d

END MODULE um_is_nan_mod
