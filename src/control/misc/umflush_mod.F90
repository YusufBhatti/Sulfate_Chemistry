! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   -------------------- SUBROUTINE umFlush ---------------------
!
!   Purpose: A wrapper script for the flush intrinsic
!
!   -------------------------------------------------------------------
!
!   Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Misc


MODULE umFlush_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UMFLUSH_MOD'

CONTAINS
SUBROUTINE umFlush(lunit,icode)

#if defined(LINUX_NAG_COMPILER)
  ! Required if using the NAG compiler
USE f90_unix_io, ONLY: flush
#elif defined(INTEL_FORTRAN) && (INTEL_FORTRAN < 12000000)
! Required for versions of ifort which don't implement the standard
! flush intrinsic. The version number is a guess.
USE ifcore
#endif

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE um_types

IMPLICIT NONE

!     The subroutine's arguments, whatever the compiler
INTEGER, INTENT(IN)  :: lunit
INTEGER, INTENT(OUT) :: icode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UMFLUSH'

!32-bit arguments needed for some platforms.
#if defined(LINUX_NAG_COMPILER) || defined(_X1) || defined(XD1) \
|| defined(XT3)
INTEGER(KIND=integer32) :: icode1
INTEGER(KIND=integer32) :: lunit1
#endif

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!     If on NAG or X1, require two 32 bit arguments to flush.
#if defined(LINUX_NAG_COMPILER) || defined(XT3)

lunit1 = lunit
CALL flush(lunit1,icode1)
icode = icode1

!     All others use one 64 bit argument

#elif defined(INTEL_FORTRAN) && (INTEL_FORTRAN < 12000000)
    ! Finding information on when Intel implemented the standard instrinsic
    ! has proved hard. I have assumed that anything before our current version
    ! didn't have it.
commitqq(lunit)
icode = 0
#elif defined(GNU_FORTRAN) && (GNU_FORTRAN < 4000000)
! Version 4.0.0 of GNU Fortran brought support for the standard approach
! so this non-standard method is only used for earlier versions.
CALL flush(lunit)
icode = 0
#else
! This is the Fortran 2003 standard for flushing a unit.
flush(lunit)
icode = 0
#endif
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE umFlush
END MODULE umFlush_mod
