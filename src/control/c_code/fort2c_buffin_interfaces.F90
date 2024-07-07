! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! This module contains the interfaces to call c code within fortran.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: C Code
!

MODULE fort2c_buffin_interfaces

! DEPENDS ON: portio2a.o
! DEPENDS ON: portio2b.o

USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_LOC, C_PTR
USE um_types, ONLY: real32, real64, integer32, integer64, logical64, logical32,&
                    integer16

IMPLICIT NONE

PRIVATE

PUBLIC ::                                                                      &

  ! Interfaces
  buffin8_single,                                                              &
  buffin32_single,                                                             &
  buffin16_single,                                                             &
  buffin64_single,                                                             &
  buffin_single,                                                               &

  ! Error Codes
  Buffin_Arguments_Error

!------------------------------------------------------------------------------!
! Interfaces to wrap buffin8/32/64 to be backwards compatible                  !
!------------------------------------------------------------------------------!

INTERFACE buffin_single
MODULE PROCEDURE                                                               &
  buffin32_single_r32,                                                         &
  buffin32_single_i,                                                           &
  buffin32_single_i_scalar,                                                    &
  buffin32_single_i2D,                                                         &
  buffin32_single_r2D,                                                         &
  buffin32_single_r3D,                                                         &
  buffin32_single_r4D,                                                         &
  buffin32_single_r5D,                                                         &
  buffin32_single_l,                                                           &
  buffin64_single_r,                                                           &
  buffin64_single_i,                                                           &
  buffin64_single_l2D,                                                         &
  buffin64_single_r2D,                                                         &
  buffin64_single_i2D,                                                         &
  buffin64_single_l,                                                           &
  buffin8_single_char
END INTERFACE

!------------------------------------------------------------------------------!

INTERFACE buffin8_single

SUBROUTINE c_buffin8_single(unit, array, maxlen, length, status)               &
           BIND(c,NAME="buffin8_single")

IMPORT :: C_PTR

IMPLICIT NONE

INTEGER :: unit, maxlen, length
REAL :: status
TYPE(C_PTR), VALUE :: array

END SUBROUTINE

MODULE PROCEDURE                                                               &
  buffin8_single_char

END INTERFACE

!------------------------------------------------------------------------------!

INTERFACE buffin16_single

SUBROUTINE c_buffin16_single(unit, array, maxlen, length, status)              &
           BIND(c,NAME="buffin16_single")

IMPORT :: C_PTR

IMPLICIT NONE

INTEGER :: unit, maxlen, length
REAL :: status
TYPE(C_PTR), VALUE :: array

END SUBROUTINE

MODULE PROCEDURE                                                               &
  buffin16_single_i

END INTERFACE


!------------------------------------------------------------------------------!

INTERFACE buffin32_single

SUBROUTINE c_buffin32_single(unit, array, maxlen, length, status)              &
           BIND(c,NAME="buffin32_single")

IMPORT :: C_PTR

IMPLICIT NONE

INTEGER :: unit, maxlen, length
REAL :: status
TYPE(C_PTR), VALUE :: array

END SUBROUTINE

MODULE PROCEDURE                                                               &
  buffin32_single_r32,                                                         &
  buffin32_single_r64,                                                         &
  buffin32_single_i,                                                           &
  buffin32_single_i_scalar,                                                    &
  buffin32_single_i2D,                                                         &
  buffin32_single_r2D,                                                         &
  buffin32_single_r3D,                                                         &
  buffin32_single_r4D,                                                         &
  buffin32_single_r5D,                                                         &
  buffin32_single_l

END INTERFACE

!------------------------------------------------------------------------------!

INTERFACE buffin64_single

SUBROUTINE c_buffin64_single(unit, array, maxlen, length, status)              &
           BIND(c,NAME="buffin64_single")

IMPORT :: C_PTR

IMPLICIT NONE

INTEGER :: unit, maxlen, length
REAL :: status
TYPE(C_PTR), VALUE :: array

END SUBROUTINE

MODULE PROCEDURE                                                               &
  buffin64_single_r,                                                           &
  buffin64_single_i,                                                           &
  buffin64_single_l2D,                                                         &
  buffin64_single_r2D,                                                         &
  buffin64_single_i2D,                                                         &
  buffin64_single_l

END INTERFACE

!------------------------------------------------------------------------------!
! Error Number Parameters                                                      !
!------------------------------------------------------------------------------!

REAL, PARAMETER ::                                                             &
  Buffin_Arguments_Error = 2.0

!------------------------------------------------------------------------------!

CONTAINS

!------------------------------------------------------------------------------!

SUBROUTINE buffin8_single_char(unit, array, maxlen, length, status)

IMPLICIT NONE

INTEGER :: unit, maxlen, length
CHARACTER, TARGET :: array(maxlen)
REAL :: status

CALL c_buffin8_single(unit, C_LOC(array), maxlen, length, status)

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE buffin32_single_l(unit, array, maxlen, length, status)

IMPLICIT NONE

INTEGER :: unit, maxlen, length
LOGICAL(KIND=logical32), TARGET :: array(maxlen)
REAL :: status

! C_LOC does not work with logical arrays (no corresponding c type)
! [ c_bool is too small; F2K3 standard requires arrays to be of
!   interoperable types                                         ]

! Therefore use address of first element of the array

! This is OK, because array will be contigous in memory, given even
! if the original argument is not, the compiler will create & use an array
! tempory which is.

CALL c_buffin32_single(unit, C_LOC(array(1)), maxlen, length, status)

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE buffin64_single_l(unit, array, maxlen, length, status)

IMPLICIT NONE

INTEGER :: unit, maxlen, length
LOGICAL(KIND=logical64), TARGET :: array(maxlen)
REAL :: status

! C_LOC does not work with logical arrays (no corresponding c type)
! [ c_bool is too small; F2K3 standard requires arrays to be of
!   interoperable types                                         ]

! Therefore use address of first element of the array

! This is OK, because array will be contigous in memory, given even
! if the original argument is not, the compiler will create & use an array
! tempory which is.

CALL c_buffin64_single(unit, C_LOC(array(1)), maxlen, length, status)

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE buffin16_single_i(unit, array, maxlen, length, status)

IMPLICIT NONE

INTEGER :: unit, maxlen, length
INTEGER(KIND=integer16), TARGET :: array(maxlen)
REAL :: status

CALL c_buffin16_single(unit, C_LOC(array), maxlen, length, status)

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE buffin32_single_i(unit, array, maxlen, length, status)

IMPLICIT NONE

INTEGER :: unit, maxlen, length
INTEGER(KIND=integer32), TARGET :: array(maxlen)
REAL :: status

CALL c_buffin32_single(unit, C_LOC(array), maxlen, length, status)

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE buffin32_single_i_scalar(unit, scalar, maxlen, length, status)

IMPLICIT NONE

INTEGER :: unit, maxlen, length
INTEGER(KIND=integer32), TARGET :: scalar
REAL :: status

IF (maxlen/=1) THEN
  status = Buffin_Arguments_Error
ELSE
  CALL c_buffin32_single(unit, C_LOC(scalar), 1, length, status)
END IF

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE buffin64_single_i(unit, array, maxlen, length, status)

IMPLICIT NONE

INTEGER :: unit, maxlen, length
INTEGER(KIND=integer64), TARGET :: array(maxlen)
REAL :: status

CALL c_buffin64_single(unit, C_LOC(array), maxlen, length, status)

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE buffin32_single_i2D(unit, array, maxlen, length, status)

IMPLICIT NONE

INTEGER :: unit, maxlen, length
INTEGER(KIND=integer32), TARGET :: array(:,:)
REAL :: status

! C_LOC does not work with assumed shape arrays (may be an array slice
! with non-unitary stride, or stored discontiguously in memory)

! Therefore use address of (1,1)th element of the array

! NB: this means the user must ensure the array is contiguous in memory
!     and stride one

CALL c_buffin32_single(unit, C_LOC(array(1,1)), maxlen, length, status)

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE buffin32_single_r2D(unit, array, maxlen, length, status)

IMPLICIT NONE

INTEGER :: unit, maxlen, length
REAL(KIND=real32), TARGET :: array(:,:)
REAL :: status

! C_LOC does not work with assumed shape arrays (may be an array slice
! with non-unitary stride, or stored discontiguously in memory)

! Therefore use address of (1,1)th element of the array

! NB: this means the user must ensure the array is contiguous in memory
!     and stride one

CALL c_buffin32_single(unit, C_LOC(array(1,1)), maxlen, length, status)

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE buffin32_single_r3D(unit, array, maxlen, length, status)

IMPLICIT NONE

INTEGER :: unit, maxlen, length
REAL(KIND=real32), TARGET :: array(:,:,:)
REAL :: status

! C_LOC does not work with assumed shape arrays (may be an array slice
! with non-unitary stride, or stored discontiguously in memory)

! Therefore use address of (1,1,1)th element of the array

! NB: this means the user must ensure the array is contiguous in memory
!     and stride one

CALL c_buffin32_single(unit, C_LOC(array(1,1,1)), maxlen, length, status)

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE buffin32_single_r4D(unit, array, maxlen, length, status)

IMPLICIT NONE

INTEGER :: unit, maxlen, length
REAL(KIND=real32), TARGET :: array(:,:,:,:)
REAL :: status

! C_LOC does not work with assumed shape arrays (may be an array slice
! with non-unitary stride, or stored discontiguously in memory)

! Therefore use address of (1,1,1,1)th element of the array

! NB: this means the user must ensure the array is contiguous in memory
!     and stride one

CALL c_buffin32_single(unit, C_LOC(array(1,1,1,1)), maxlen, length, status)

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE buffin32_single_r5D(unit, array, maxlen, length, status)

IMPLICIT NONE

INTEGER :: unit, maxlen, length
REAL(KIND=real32), TARGET :: array(:,:,:,:,:)
REAL :: status

! C_LOC does not work with assumed shape arrays (may be an array slice
! with non-unitary stride, or stored discontiguously in memory)

! Therefore use address of (1,1,1,1,1)th element of the array

! NB: this means the user must ensure the array is contiguous in memory
!     and stride one

CALL c_buffin32_single(unit, C_LOC(array(1,1,1,1,1)), maxlen, length, status)

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE buffin32_single_r64(unit, array, maxlen, length, status)

IMPLICIT NONE

INTEGER :: unit, maxlen, length
REAL(KIND=real64), TARGET :: array((maxlen+1)/2)
REAL :: status

CALL c_buffin32_single(unit, C_LOC(array), maxlen, length, status)

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE buffin32_single_r32(unit, array, maxlen, length, status)

IMPLICIT NONE

INTEGER :: unit, maxlen, length
REAL(KIND=real32), TARGET :: array(maxlen)
REAL :: status

CALL c_buffin32_single(unit, C_LOC(array), maxlen, length, status)

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE buffin64_single_r(unit, array, maxlen, length, status)

IMPLICIT NONE

INTEGER :: unit, maxlen, length
REAL(KIND=real64), TARGET :: array(maxlen)
REAL :: status

CALL c_buffin64_single(unit, C_LOC(array), maxlen, length, status)

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE buffin64_single_r2D(unit, array, maxlen, length, status)

IMPLICIT NONE

INTEGER :: unit, maxlen, length
REAL(KIND=real64), TARGET :: array(:,:)
REAL :: status

! C_LOC does not work with assumed shape arrays (may be an array slice
! with non-unitary stride, or stored discontiguously in memory)

! Therefore use address of (1,1)th element of the array

! NB: this means the user must ensure the array is contiguous in memory
!     and stride one

CALL c_buffin64_single(unit, C_LOC(array(1,1)), maxlen, length, status)

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE buffin64_single_i2D(unit, array, maxlen, length, status)

IMPLICIT NONE

INTEGER :: unit, maxlen, length
INTEGER(KIND=integer64), TARGET :: array(:,:)
REAL :: status

! C_LOC does not work with assumed shape arrays (may be an array slice
! with non-unitary stride, or stored discontiguously in memory)

! Therefore use address of (1,1)th element of the array

! NB: this means the user must ensure the array is contiguous in memory
!     and stride one

CALL c_buffin64_single(unit, C_LOC(array(1,1)), maxlen, length, status)

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE buffin64_single_l2D(unit, array, maxlen, length, status)

IMPLICIT NONE

INTEGER :: unit, maxlen, length
LOGICAL(KIND=logical64), TARGET :: array(:,:)
REAL :: status

! C_LOC does not work with assumed shape arrays (may be an array slice
! with non-unitary stride, or stored discontiguously in memory)

! Therefore use address of (1,1)th element of the array

! NB: this means the user must ensure the array is contiguous in memory
!     and stride one

CALL c_buffin64_single(unit, C_LOC(array(1,1)), maxlen, length, status)

END SUBROUTINE

!------------------------------------------------------------------------------!

END MODULE

