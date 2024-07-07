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

MODULE fort2c_memcpy_interfaces

! DEPENDS ON: portutils.o

USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_LOC, C_PTR, C_INT64_T
USE um_types, ONLY: real32, real64, integer32, integer64, logical64, logical32
USE ios_constants, ONLY: IOS_tukind
USE IOS_types, ONLY: IOS_status

IMPLICIT NONE

PRIVATE

PUBLIC :: um_memcpy64, um_memcpy32

! -----------------------------------------------------------------------------!

INTERFACE um_memcpy64
SUBROUTINE c_um_memcpy64(dest,src,n) BIND(c,NAME="um_memcpy64")

import :: C_PTR, C_INT64_T

IMPLICIT NONE

TYPE(C_PTR), VALUE :: dest
TYPE(C_PTR), VALUE :: src
INTEGER(KIND=C_INT64_T), VALUE :: n

END SUBROUTINE

MODULE PROCEDURE um_memcpy64_IOS_Status, um_memcpy64_r, um_memcpy64_l,     &
                 um_memcpy64_i, um_memcpy64_r2D, um_memcpy64_l2D,          &
                 um_memcpy64_i2D, um_memcpy64_IOS_tukind_r2D,              &
                 um_memcpy64_IOS_tukind_r, um_memcpy64_IOS_tukind_l,       &
                 um_memcpy64_IOS_tukind_i, um_memcpy64_IOS_tukind_i2D,     &
                 um_memcpy64_IOS_tukind_l2D, um_memcpy64_i32,              &
                 um_memcpy64_r32_r2D, um_memcpy64_IOS_tukind_i2D32,        &
                 um_memcpy64_IOS_tukind_r32, um_memcpy64_IOS_tukind_l32,   &
                 um_memcpy64_r_r32, um_memcpy64_r32_r
END INTERFACE

INTERFACE um_memcpy32
SUBROUTINE c_um_memcpy32(dest,src,n) BIND(c,NAME="um_memcpy32")

import :: C_PTR, C_INT64_T

IMPLICIT NONE

TYPE(C_PTR), VALUE :: dest
TYPE(C_PTR), VALUE :: src
INTEGER(KIND=C_INT64_T), VALUE :: n

END SUBROUTINE

MODULE PROCEDURE um_memcpy32_r, um_memcpy32_i, um_memcpy32_i2D,            &
                 um_memcpy32_l, um_memcpy32_IOS_tukind_i2D,                &
                 um_memcpy32_IOS_tukind_r, um_memcpy32_r_r,                &
                 um_memcpy32_r64_r
END INTERFACE

! -----------------------------------------------------------------------------!

CONTAINS

! -----------------------------------------------------------------------------!

SUBROUTINE um_memcpy64_IOS_Status(dest,src,n)

IMPLICIT NONE

TYPE(IOS_Status), TARGET         :: dest
INTEGER(KIND=IOS_tukind), TARGET :: src(*)
INTEGER :: n

CALL c_um_memcpy64(C_LOC(dest),C_LOC(src),INT(n,KIND=C_INT64_T))

END SUBROUTINE

SUBROUTINE um_memcpy32_r(dest,src,n)

IMPLICIT NONE

REAL(KIND=real32), TARGET        :: dest(*)
INTEGER(KIND=IOS_tukind), TARGET :: src(*)
INTEGER :: n

CALL c_um_memcpy32(C_LOC(dest),C_LOC(src),INT(n,KIND=C_INT64_T))

END SUBROUTINE

SUBROUTINE um_memcpy32_r_r(dest,src,n)

IMPLICIT NONE

REAL(KIND=real32), TARGET :: dest(*)
REAL(KIND=real32), TARGET :: src(*)
INTEGER :: n

CALL c_um_memcpy32(C_LOC(dest),C_LOC(src),INT(n,KIND=C_INT64_T))

END SUBROUTINE

SUBROUTINE um_memcpy32_r64_r(dest,src,n)

IMPLICIT NONE

REAL(KIND=real64), TARGET :: dest(*)
REAL(KIND=real32), TARGET :: src(*)
INTEGER :: n

CALL c_um_memcpy32(C_LOC(dest),C_LOC(src),INT(n,KIND=C_INT64_T))

END SUBROUTINE

SUBROUTINE um_memcpy32_IOS_tukind_r(dest,src,n)

IMPLICIT NONE

REAL(KIND=real32), TARGET        :: src(*)
INTEGER(KIND=IOS_tukind), TARGET :: dest(*)
INTEGER :: n

CALL c_um_memcpy32(C_LOC(dest),C_LOC(src),INT(n,KIND=C_INT64_T))

END SUBROUTINE

SUBROUTINE um_memcpy64_IOS_tukind_r(dest,src,n)

IMPLICIT NONE

REAL(KIND=real64), TARGET        :: src(*)
INTEGER(KIND=IOS_tukind), TARGET :: dest(*)
INTEGER :: n

CALL c_um_memcpy64(C_LOC(dest),C_LOC(src),INT(n,KIND=C_INT64_T))

END SUBROUTINE

SUBROUTINE um_memcpy64_IOS_tukind_r32(dest,src,n)

IMPLICIT NONE

REAL(KIND=real32), TARGET        :: src(*)
INTEGER(KIND=IOS_tukind), TARGET :: dest(*)
INTEGER :: n

CALL c_um_memcpy64(C_LOC(dest),C_LOC(src),INT(n,KIND=C_INT64_T))

END SUBROUTINE

SUBROUTINE um_memcpy32_i(dest,src,n)

IMPLICIT NONE

INTEGER(KIND=integer32), TARGET  :: dest(*)
INTEGER(KIND=IOS_tukind), TARGET :: src(*)
INTEGER :: n

CALL c_um_memcpy32(C_LOC(dest),C_LOC(src),INT(n,KIND=C_INT64_T))

END SUBROUTINE

SUBROUTINE um_memcpy32_i2D(dest,src,n)

IMPLICIT NONE

INTEGER(KIND=integer32), TARGET  :: dest(:,:)
INTEGER(KIND=IOS_tukind), TARGET :: src(*)
INTEGER :: n

! C_LOC does not work with assumed shape arrays (may be an array slice
! with non-unitary stride, or stored discontiguously in memory)

! Thererefore use address of (1,1)th element of the array

!NB: this means the user must ensure the array is contiguos in memory
!    and stride one

CALL c_um_memcpy32(C_LOC(dest(1,1)),C_LOC(src),INT(n,KIND=C_INT64_T))

END SUBROUTINE

SUBROUTINE um_memcpy32_IOS_tukind_i2D(dest,src,n)

IMPLICIT NONE

INTEGER(KIND=integer32), TARGET  :: src(:,:)
INTEGER(KIND=IOS_tukind), TARGET :: dest(*)
INTEGER :: n

! C_LOC does not work with assumed shape arrays (may be an array slice
! with non-unitary stride, or stored discontiguously in memory)

! Thererefore use address of (1,1)th element of the array

!NB: this means the user must ensure the array is contiguos in memory
!    and stride one

CALL c_um_memcpy32(C_LOC(dest),C_LOC(src(1,1)),INT(n,KIND=C_INT64_T))

END SUBROUTINE


SUBROUTINE um_memcpy64_r2D(dest,src,n)

IMPLICIT NONE

REAL(KIND=real64), TARGET        :: dest(:,:)
INTEGER(KIND=IOS_tukind), TARGET :: src(*)
INTEGER :: n

! C_LOC does not work with assumed shape arrays (may be an array slice
! with non-unitary stride, or stored discontiguously in memory)

! Thererefore use address of (1,1)th element of the array

!NB: this means the user must ensure the array is contiguos in memory
!    and stride one

CALL c_um_memcpy64(C_LOC(dest(1,1)),C_LOC(src),INT(n,KIND=C_INT64_T))

END SUBROUTINE

SUBROUTINE um_memcpy64_IOS_tukind_r2D(dest,src,n)

IMPLICIT NONE

REAL(KIND=real64), TARGET        :: src(:,:)
INTEGER(KIND=IOS_tukind), TARGET :: dest(*)
INTEGER :: n

! C_LOC does not work with assumed shape arrays (may be an array slice
! with non-unitary stride, or stored discontiguously in memory)

! Thererefore use address of (1,1)th element of the array

!NB: this means the user must ensure the array is contiguos in memory
!    and stride one

CALL c_um_memcpy64(C_LOC(dest),C_LOC(src(1,1)),INT(n,KIND=C_INT64_T))

END SUBROUTINE

SUBROUTINE um_memcpy64_i2D(dest,src,n)

IMPLICIT NONE

INTEGER(KIND=integer64), TARGET  :: dest(:,:)
INTEGER(KIND=IOS_tukind), TARGET :: src(*)
INTEGER :: n

! C_LOC does not work with assumed shape arrays (may be an array slice
! with non-unitary stride, or stored discontiguously in memory)

! Thererefore use address of (1,1)th element of the array

!NB: this means the user must ensure the array is contiguos in memory
!    and stride one

CALL c_um_memcpy64(C_LOC(dest(1,1)),C_LOC(src),INT(n,KIND=C_INT64_T))

END SUBROUTINE

SUBROUTINE um_memcpy64_IOS_tukind_i2D(dest,src,n)

IMPLICIT NONE

INTEGER(KIND=integer64), TARGET  :: src(:,:)
INTEGER(KIND=IOS_tukind), TARGET :: dest(*)
INTEGER :: n

! C_LOC does not work with assumed shape arrays (may be an array slice
! with non-unitary stride, or stored discontiguously in memory)

! Thererefore use address of (1,1)th element of the array

!NB: this means the user must ensure the array is contiguos in memory
!    and stride one

CALL c_um_memcpy64(C_LOC(dest),C_LOC(src(1,1)),INT(n,KIND=C_INT64_T))

END SUBROUTINE

SUBROUTINE um_memcpy64_IOS_tukind_i2D32(dest,src,n)

IMPLICIT NONE

INTEGER(KIND=integer32), TARGET  :: src(:,:)
INTEGER(KIND=IOS_tukind), TARGET :: dest(*)
INTEGER :: n

! C_LOC does not work with assumed shape arrays (may be an array slice
! with non-unitary stride, or stored discontiguously in memory)

! Thererefore use address of (1,1)th element of the array

!NB: this means the user must ensure the array is contiguos in memory
!    and stride one

CALL c_um_memcpy64(C_LOC(dest),C_LOC(src(1,1)),INT(n,KIND=C_INT64_T))

END SUBROUTINE


SUBROUTINE um_memcpy64_l2D(dest,src,n)

IMPLICIT NONE

LOGICAL(KIND=logical64), TARGET  :: dest(:,:)
INTEGER(KIND=IOS_tukind), TARGET :: src(*)
INTEGER :: n

! C_LOC does not work with assumed shape arrays (may be an array slice
! with non-unitary stride, or stored discontiguously in memory)

! Thererefore use address of (1,1)th element of the array

!NB: this means the user must ensure the array is contiguos in memory
!    and stride one

CALL c_um_memcpy64(C_LOC(dest(1,1)),C_LOC(src),INT(n,KIND=C_INT64_T))

END SUBROUTINE

SUBROUTINE um_memcpy64_IOS_tukind_l2D(dest,src,n)

IMPLICIT NONE

LOGICAL(KIND=logical64), TARGET  :: src(:,:)
INTEGER(KIND=IOS_tukind), TARGET :: dest(*)
INTEGER :: n

! C_LOC does not work with assumed shape arrays (may be an array slice
! with non-unitary stride, or stored discontiguously in memory)

! Thererefore use address of (1,1)th element of the array

!NB: this means the user must ensure the array is contiguos in memory
!    and stride one

CALL c_um_memcpy64(C_LOC(dest),C_LOC(src(1,1)),INT(n,KIND=C_INT64_T))

END SUBROUTINE

SUBROUTINE um_memcpy64_r(dest,src,n)

IMPLICIT NONE

REAL(KIND=real64), TARGET        :: dest(*)
INTEGER(KIND=IOS_tukind), TARGET :: src(*)
INTEGER :: n

CALL c_um_memcpy64(C_LOC(dest),C_LOC(src),INT(n,KIND=C_INT64_T))

END SUBROUTINE

SUBROUTINE um_memcpy64_r_r32(dest,src,n)

IMPLICIT NONE

REAL(KIND=real64), TARGET :: dest(*)
REAL(KIND=real32), TARGET :: src(*)
INTEGER :: n

CALL c_um_memcpy64(C_LOC(dest),C_LOC(src),INT(n,KIND=C_INT64_T))

END SUBROUTINE

SUBROUTINE um_memcpy64_r32_r(dest,src,n)

IMPLICIT NONE

REAL(KIND=real32), TARGET :: dest(*)
REAL(KIND=real64), TARGET :: src(*)
INTEGER :: n

CALL c_um_memcpy64(C_LOC(dest),C_LOC(src),INT(n,KIND=C_INT64_T))

END SUBROUTINE

SUBROUTINE um_memcpy64_r32_r2D(dest,src,n)

IMPLICIT NONE

REAL(KIND=real32), TARGET :: dest(*)
REAL(KIND=real64), TARGET :: src(:,:)
INTEGER :: n

! C_LOC does not work with assumed shape arrays (may be an array slice
! with non-unitary stride, or stored discontiguously in memory)

! Thererefore use address of (1,1)th element of the array

!NB: this means the user must ensure the array is contiguos in memory
!    and stride one

CALL c_um_memcpy64(C_LOC(dest),C_LOC(src(1,1)),INT(n,KIND=C_INT64_T))

END SUBROUTINE

SUBROUTINE um_memcpy64_l(dest,src,n)

IMPLICIT NONE

LOGICAL(KIND=logical64), TARGET  :: dest(*)
INTEGER(KIND=IOS_tukind), TARGET :: src(*)
INTEGER :: n

! C_LOC does not work with logical arrays (no corresponding c type)
! [ c_bool is too small; F2K3 standard requires arrays to be of
!   interoperable types                                         ]

! Thererefore use address of first element of the array

! This is OK, because array will be contigous in memory, given even
! if the original argument is not, the compiler will create & use an array
! tempory which is.

CALL c_um_memcpy64(C_LOC(dest(1)),C_LOC(src),INT(n,KIND=C_INT64_T))

END SUBROUTINE

SUBROUTINE um_memcpy64_IOS_tukind_l(dest,src,n)

IMPLICIT NONE

LOGICAL(KIND=logical64), TARGET  :: src(*)
INTEGER(KIND=IOS_tukind), TARGET :: dest(*)
INTEGER :: n

! C_LOC does not work with logical arrays (no corresponding c type)
! [ c_bool is too small; F2K3 standard requires arrays to be of
!   interoperable types                                         ]

! Thererefore use address of first element of the array

! This is OK, because array will be contigous in memory, given even
! if the original argument is not, the compiler will create & use an array
! tempory which is.

CALL c_um_memcpy64(C_LOC(dest),C_LOC(src(1)),INT(n,KIND=C_INT64_T))

END SUBROUTINE

SUBROUTINE um_memcpy64_IOS_tukind_l32(dest,src,n)

IMPLICIT NONE

LOGICAL(KIND=logical32), TARGET  :: src(*)
INTEGER(KIND=IOS_tukind), TARGET :: dest(*)
INTEGER :: n

! C_LOC does not work with logical arrays (no corresponding c type)
! [ c_bool is too small; F2K3 standard requires arrays to be of
!   interoperable types                                         ]

! Thererefore use address of first element of the array

! This is OK, because array will be contigous in memory, given even
! if the original argument is not, the compiler will create & use an array
! tempory which is.

CALL c_um_memcpy64(C_LOC(dest),C_LOC(src(1)),INT(n,KIND=C_INT64_T))

END SUBROUTINE

SUBROUTINE um_memcpy64_i(dest,src,n)

IMPLICIT NONE

INTEGER(KIND=integer64), TARGET  :: dest(*)
INTEGER(KIND=IOS_tukind), TARGET :: src(*)
INTEGER :: n

CALL c_um_memcpy64(C_LOC(dest),C_LOC(src),INT(n,KIND=C_INT64_T))

END SUBROUTINE

SUBROUTINE um_memcpy64_i32(dest,src,n)

IMPLICIT NONE

INTEGER(KIND=integer32), TARGET  :: dest(*)
INTEGER(KIND=IOS_tukind), TARGET :: src(*)
INTEGER :: n

CALL c_um_memcpy64(C_LOC(dest),C_LOC(src),INT(n,KIND=C_INT64_T))

END SUBROUTINE

SUBROUTINE um_memcpy64_IOS_tukind_i(dest,src,n)

IMPLICIT NONE

INTEGER(KIND=integer64), TARGET  :: src(*)
INTEGER(KIND=IOS_tukind), TARGET :: dest(*)
INTEGER :: n

CALL c_um_memcpy64(C_LOC(dest),C_LOC(src),INT(n,KIND=C_INT64_T))

END SUBROUTINE

SUBROUTINE um_memcpy32_l(dest,src,n)

IMPLICIT NONE

LOGICAL(KIND=logical32), TARGET  :: dest(*)
INTEGER(KIND=IOS_tukind), TARGET :: src(*)
INTEGER :: n

! C_LOC does not work with logical arrays (no corresponding c type)
! [ c_bool is too small; F2K3 standard requires arrays to be of
!   interoperable types                                         ]

! Thererefore use address of first element of the array

! This is OK, because array will be contigous in memory, given even
! if the original argument is not, the compiler will create & use an array
! tempory which is.

CALL c_um_memcpy32(C_LOC(dest(1)),C_LOC(src),INT(n,KIND=C_INT64_T))

END SUBROUTINE

END MODULE

