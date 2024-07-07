! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Misc

! Parameters for 32 and 64 bit kinds

MODULE um_types

#if !defined(LFRIC)
USE fort2c_addr_mod, ONLY: um_addr_diff, um_addr_size
#endif

USE, INTRINSIC :: ISO_C_BINDING, ONLY:                                         &
#if !defined(LFRIC)
  C_LOC, C_INTPTR_T,                                                           &
#endif
  C_INT8_T, C_INT16_T

!$ USE omp_lib

IMPLICIT NONE
! Precision and range for 64 bit real
INTEGER, PARAMETER :: prec64  = 15
INTEGER, PARAMETER :: range64 = 307

! Precision and range for 32 bit real
INTEGER, PARAMETER :: prec32  = 6
INTEGER, PARAMETER :: range32 = 37

! Range for integers
INTEGER, PARAMETER :: irange64=15
INTEGER, PARAMETER :: irange32=9

! Range for small logicals
INTEGER, PARAMETER :: lrange1=1

! Kind for 64 bit real
INTEGER, PARAMETER :: real64  = SELECTED_REAL_KIND(prec64,range64)
! Kind for 32 bit real
INTEGER, PARAMETER :: real32  = SELECTED_REAL_KIND(prec32,range32)
! Kind for 64 bit integer
INTEGER, PARAMETER :: integer64 = SELECTED_INT_KIND(irange64)
! Kind for 32 bit integer
INTEGER, PARAMETER :: integer32 = SELECTED_INT_KIND(irange32)
! Kind for 16 bit integer
INTEGER, PARAMETER :: integer16 = C_INT16_T
! Kind for 8 bit integer
INTEGER, PARAMETER :: integer8 = C_INT8_T

!Scheme-specific precisions

!Large scale precipitation scheme
#if defined(LSPREC_32B)
INTEGER, PARAMETER :: real_lsprec = real32
#else
INTEGER, PARAMETER :: real_lsprec = real64
#endif

! Kind for use with OpenMP functions (from omp_lib).
! Equal to KIND(0) if no OpenMP, else is equal to KIND(openmp_version)
INTEGER, PARAMETER :: integer_omp = KIND(0)                                    &
!$                                  * 0 + KIND(openmp_version)                 &
                                    + 0

! Kinds for 64 and 32 bit logicals. Note that there is no
! "selected_logical_kind", but using the equivalent integer kind is a
! workaround that works on every platform we have tested.
INTEGER, PARAMETER :: logical64 = integer64
INTEGER, PARAMETER :: logical32 = integer32
! Kind for small logicals
INTEGER, PARAMETER :: log_small = SELECTED_INT_KIND(lrange1)

#if !defined(LFRIC)
CONTAINS

! Discover the size of a fortran default real
FUNCTION umFortranRealSize() RESULT (r)

IMPLICIT NONE
INTEGER :: r
INTEGER(KIND=C_INTPTR_T) :: rintptr
REAL, TARGET    :: a(2)

CALL um_addr_diff(C_LOC(a(1)),C_LOC(a(2)),rintptr)

r = INT(rintptr)

END FUNCTION umFortranRealSize

! Discover the size of a fortran default integer
FUNCTION umFortranIntegerSize() RESULT (r)

IMPLICIT NONE
INTEGER :: r
INTEGER, TARGET :: a(2)
INTEGER(KIND=C_INTPTR_T) :: rintptr

CALL um_addr_diff(C_LOC(a(1)),C_LOC(a(2)),rintptr)

r = INT(rintptr)

END FUNCTION umFortranIntegerSize

! Discover the size of a fortran default logical
FUNCTION umFortranLogicalSize() RESULT (r)

IMPLICIT NONE
INTEGER :: r
LOGICAL, TARGET :: a(2)
INTEGER(KIND=C_INTPTR_T) :: rintptr

CALL um_addr_diff(C_LOC(a(1)),C_LOC(a(2)),rintptr)

r = INT(rintptr)

END FUNCTION umFortranLogicalSize

! Discover the size of a pointer (we assume c and fortran are the same)
FUNCTION umFortranPointerSize() RESULT (r)
IMPLICIT NONE
INTEGER :: r
r=um_addr_size()
END FUNCTION umFortranPointerSize
#endif

END MODULE um_types

