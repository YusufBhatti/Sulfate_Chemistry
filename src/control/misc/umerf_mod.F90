! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   -------------------- MODULE umErf_mod -----------------------------
!
!   Purpose: Provides a portable wrapper function for the erf intrinsic
!
!   -------------------------------------------------------------------
!
!   Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Misc

MODULE umErf_mod

USE, INTRINSIC :: iso_c_binding, ONLY: C_DOUBLE

IMPLICIT NONE

PRIVATE
PUBLIC :: umErf

#if defined(UMERF_USE_LIBC)

INTERFACE
FUNCTION c_erf(x) BIND(c,NAME="erf")

IMPORT :: C_DOUBLE

IMPLICIT NONE

REAL(KIND=C_DOUBLE)             :: c_erf
REAL(KIND=C_DOUBLE), INTENT(IN) :: x

END FUNCTION
END INTERFACE

! Use the libc standard C99 version

INTERFACE umErf
MODULE PROCEDURE umErf_scalar, umErf_vector
END INTERFACE

#endif

!------------------------------------------------------------------------------!
CONTAINS
!------------------------------------------------------------------------------!

#if !defined(UMERF_USE_LIBC)

ELEMENTAL FUNCTION umErf(x)

  IMPLICIT NONE

  REAL             :: umErf
  REAL, INTENT(IN) :: x

  ! Use the non-F2k3 standard version - which is ELEMENTAL
  umErf = ERF(x)

END FUNCTION umErf

!------------------------------------------------------------------------------!

#else

FUNCTION umErf_scalar(x)

IMPLICIT NONE

REAL, INTENT(IN)            :: x

REAL                        :: umErf_scalar

umErf_scalar = c_erf(REAL(x,KIND=c_double))

END FUNCTION umErf_scalar

!------------------------------------------------------------------------------!

FUNCTION umErf_vector(x)

IMPLICIT NONE

REAL, INTENT(IN)            :: x(:)

REAL                        :: umErf_vector(1:SIZE(x))

INTEGER                     :: i

DO i=1,SIZE(x)
  umErf_vector(i)=c_erf(REAL(x(i),KIND=c_double))
END DO

END FUNCTION umErf_vector

!------------------------------------------------------------------------------!

#endif

END MODULE umErf_mod
