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

MODULE fort2c_buffo_interfaces

! DEPENDS ON: portio2a.o
! DEPENDS ON: portio2b.o

USE, INTRINSIC :: ISO_C_BINDING, ONLY:                                         &
  C_LOC,                                                                       &
  C_PTR

USE um_types, ONLY:                                                            &
  real32,                                                                      &
  real64,                                                                      &
  integer8,                                                                    &
  integer32,                                                                   &
  integer64,                                                                   &
  logical64

IMPLICIT NONE

PRIVATE

PUBLIC ::                                                                      &

  ! Interfaces
  buffout_single,                                                              &
  buffout8_single,                                                             &
  buffout32_single,                                                            &
  buffout64_single,                                                            &

  ! Error Codes
  Buffo_Arguments_Error

!------------------------------------------------------------------------------!
! Interfaces to wrap buffout32/64 to be backwards compatible                   !
!------------------------------------------------------------------------------!

INTERFACE buffout_single
MODULE PROCEDURE                                                               &
  buffout64_single_r,                                                          &
  buffout64_single_i,                                                          &
  buffout64_single_l,                                                          &
  buffout64_single_i2D,                                                        &
  buffout64_single_r2D,                                                        &
  buffout64_single_l2D,                                                        &
  buffout32_single_r32,                                                        &
  buffout32_single_i,                                                          &
  buffout32_single_i_scalar,                                                   &
  buffout32_single_r2D,                                                        &
  buffout32_single_r3D,                                                        &
  buffout32_single_r4D,                                                        &
  buffout32_single_r5D,                                                        &
  buffout32_single_i2D
END INTERFACE

!------------------------------------------------------------------------------!

INTERFACE buffout64_single

SUBROUTINE c_buffout64_single(unit, array, maxlen, length, status)             &
                 BIND(c,NAME="buffout64_single")

IMPORT :: C_PTR

IMPLICIT NONE

INTEGER :: unit, maxlen, length
TYPE(C_PTR), VALUE :: array
REAL :: status

END SUBROUTINE

MODULE PROCEDURE                                                               &
  buffout64_single_r,                                                          &
  buffout64_single_i,                                                          &
  buffout64_single_l,                                                          &
  buffout64_single_i2D,                                                        &
  buffout64_single_r2D,                                                        &
  buffout64_single_l2D

END INTERFACE

!------------------------------------------------------------------------------!

INTERFACE buffout32_single

SUBROUTINE c_buffout32_single(unit, array, maxlen, length, status)             &
           BIND(c,NAME="buffout32_single")

IMPORT :: C_PTR

IMPLICIT NONE

INTEGER :: unit, maxlen, length
REAL :: status
TYPE(C_PTR), VALUE :: array
END SUBROUTINE

MODULE PROCEDURE                                                               &
  buffout32_single_r32,                                                        &
  buffout32_single_r64,                                                        &
  buffout32_single_i,                                                          &
  buffout32_single_i_scalar,                                                   &
  buffout32_single_r2D,                                                        &
  buffout32_single_r3D,                                                        &
  buffout32_single_r4D,                                                        &
  buffout32_single_r5D,                                                        &
  buffout32_single_i2D

END INTERFACE

!------------------------------------------------------------------------------!

INTERFACE buffout8_single

SUBROUTINE c_buffout8_single(unit, array, maxlen, length, status)              &
           BIND(c,NAME="buffout8_single")

IMPORT :: C_PTR

IMPLICIT NONE

INTEGER :: unit, maxlen, length
REAL :: status
TYPE(C_PTR), VALUE :: array
END SUBROUTINE

MODULE PROCEDURE                                                               &
  buffout8_single_i

END INTERFACE

!------------------------------------------------------------------------------!
! Error Number Parameters                                                      !
!------------------------------------------------------------------------------!

REAL, PARAMETER ::                                                             &
  Buffo_Arguments_Error = 2.0

!------------------------------------------------------------------------------!

CONTAINS

!------------------------------------------------------------------------------!

SUBROUTINE buffout64_single_r(unit, array, maxlen, length, status)

IMPLICIT NONE
INTEGER :: unit, maxlen, length
REAL(KIND=real64), TARGET :: array(maxlen)
REAL :: status

CALL c_buffout64_single(unit, C_LOC(array), maxlen, length, status)

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE buffout64_single_i(unit, array, maxlen, length, status)

IMPLICIT NONE

INTEGER :: unit, maxlen, length
INTEGER(KIND=integer64), TARGET :: array(maxlen)
REAL :: status

CALL c_buffout64_single(unit, C_LOC(array), maxlen, length, status)

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE buffout64_single_l(unit, array, maxlen, length, status)

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

CALL c_buffout64_single(unit, C_LOC(array(1)), maxlen, length, status)

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE buffout64_single_r2D(unit, array, maxlen, length, status)

IMPLICIT NONE
INTEGER :: unit, maxlen, length
REAL(KIND=real64), TARGET :: array(:,:)
REAL :: status

! C_LOC does not work with assumed shape arrays (may be an array slice
! with non-unitary stride, or stored discontiguously in memory)

! Therefore use address of (1,1)th element of the array

! NB: this means the user must ensure the array is contiguous in memory
!     and stride one

CALL c_buffout64_single(unit, C_LOC(array(1,1)), maxlen, length, status)

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE buffout64_single_i2D(unit, array, maxlen, length, status)

IMPLICIT NONE

INTEGER :: unit, maxlen, length
INTEGER(KIND=integer64), TARGET :: array(:,:)
REAL :: status

! C_LOC does not work with assumed shape arrays (may be an array slice
! with non-unitary stride, or stored discontiguously in memory)

! Therefore use address of (1,1)th element of the array

! NB: this means the user must ensure the array is contiguous in memory
!     and stride one

CALL c_buffout64_single(unit, C_LOC(array(1,1)), maxlen, length, status)

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE buffout64_single_l2D(unit, array, maxlen, length, status)

IMPLICIT NONE

INTEGER :: unit, maxlen, length
LOGICAL(KIND=logical64), TARGET :: array(:,:)
REAL :: status

! C_LOC does not work with assumed shape arrays (may be an array slice
! with non-unitary stride, or stored discontiguously in memory)

! Therefore use address of (1,1)th element of the array

! NB: this means the user must ensure the array is contiguous in memory
!     and stride one

CALL c_buffout64_single(unit, C_LOC(array(1,1)), maxlen, length, status)

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE buffout32_single_r64(unit, array, maxlen, length, status)

IMPLICIT NONE

INTEGER :: unit, maxlen, length
REAL(KIND=real64), TARGET :: array((maxlen+1)/2)
REAL :: status

CALL c_buffout32_single(unit, C_LOC(array), maxlen, length, status)

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE buffout32_single_r32(unit, array, maxlen, length, status)

IMPLICIT NONE

INTEGER :: unit, maxlen, length
REAL(KIND=real32), TARGET :: array(maxlen)
REAL :: status

CALL c_buffout32_single(unit, C_LOC(array), maxlen, length, status)

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE buffout32_single_i(unit, array, maxlen, length, status)

IMPLICIT NONE
INTEGER :: unit, maxlen, length
INTEGER(KIND=integer32), TARGET :: array(maxlen)
REAL :: status

CALL c_buffout32_single(unit, C_LOC(array), maxlen, length, status)

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE buffout32_single_i_scalar(unit, scalar, maxlen, length, status)

IMPLICIT NONE
INTEGER :: unit, maxlen, length
INTEGER(KIND=integer32), TARGET :: scalar
REAL :: status

IF (maxlen/=1) THEN
  status = Buffo_Arguments_Error
ELSE
  CALL c_buffout32_single(unit, C_LOC(scalar), 1, length, status)
END IF

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE buffout32_single_i2D(unit, array, maxlen, length, status)

IMPLICIT NONE
INTEGER :: unit, maxlen, length
INTEGER(KIND=integer32), TARGET :: array(:,:)
REAL :: status

! C_LOC does not work with assumed shape arrays (may be an array slice
! with non-unitary stride, or stored discontiguously in memory)

! Therefore use address of (1,1)th element of the array

! NB: this means the user must ensure the array is contiguous in memory
!     and stride one

CALL c_buffout32_single(unit, C_LOC(array(1,1)), maxlen, length, status)

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE buffout32_single_r2D(unit, array, maxlen, length, status)

IMPLICIT NONE
INTEGER :: unit, maxlen, length
REAL(KIND=real32), TARGET :: array(:,:)
REAL :: status

! C_LOC does not work with assumed shape arrays (may be an array slice
! with non-unitary stride, or stored discontiguously in memory)

! Therefore use address of (1,1)th element of the array

! NB: this means the user must ensure the array is contiguous in memory
!     and stride one

CALL c_buffout32_single(unit, C_LOC(array(1,1)), maxlen, length, status)

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE buffout32_single_r3D(unit, array, maxlen, length, status)

IMPLICIT NONE
INTEGER :: unit, maxlen, length
REAL(KIND=real32), TARGET :: array(:,:,:)
REAL :: status

! C_LOC does not work with assumed shape arrays (may be an array slice
! with non-unitary stride, or stored discontiguously in memory)

! Therefore use address of (1,1,1)th element of the array

! NB: this means the user must ensure the array is contiguous in memory
!     and stride one

CALL c_buffout32_single(unit, C_LOC(array(1,1,1)), maxlen, length, status)

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE buffout32_single_r4D(unit, array, maxlen, length, status)

IMPLICIT NONE
INTEGER :: unit, maxlen, length
REAL(KIND=real32), TARGET :: array(:,:,:,:)
REAL :: status

! C_LOC does not work with assumed shape arrays (may be an array slice
! with non-unitary stride, or stored discontiguously in memory)

! Therefore use address of (1,1,1,1)th element of the array

! NB: this means the user must ensure the array is contiguous in memory
!     and stride one

CALL c_buffout32_single(unit, C_LOC(array(1,1,1,1)), maxlen, length, status)

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE buffout32_single_r5D(unit, array, maxlen, length, status)

IMPLICIT NONE
INTEGER :: unit, maxlen, length
REAL(KIND=real32), TARGET :: array(:,:,:,:,:)
REAL :: status

! C_LOC does not work with assumed shape arrays (may be an array slice
! with non-unitary stride, or stored discontiguously in memory)

! Therefore use address of (1,1,1,1,1)th element of the array

! NB: this means the user must ensure the array is contiguous in memory
!     and stride one

CALL c_buffout32_single(unit, C_LOC(array(1,1,1,1,1)), maxlen, length, status)

END SUBROUTINE

!------------------------------------------------------------------------------!

SUBROUTINE buffout8_single_i(unit, array, maxlen, length, status)

IMPLICIT NONE
INTEGER :: unit, maxlen, length
INTEGER(KIND=integer8), TARGET :: array(maxlen)
REAL :: status

CALL c_buffout8_single(unit, C_LOC(array), maxlen, length, status)

END SUBROUTINE

!------------------------------------------------------------------------------!

END MODULE