! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!

MODULE vectlib_mod

! Description:
! This routine acts as an interface to vector versions
! of intrinsics functions on a platform.
!
! Supported libraries:
! IBM's VMASS (compile using VMASS def)
! Documentation: http://www-01.ibm.com/support/docview.wss?uid=swg27005473
!
! Intel's MKL (compile with MKL def)
! Documentation:
! https://software.intel.com/en-us/mkl-developer-reference-fortran
!
! Default compiles to equivalent do loop over array
!
! INTERFACEs for 32/64 bit versions created as required. Implemented so far:
! -powr_v
! -exp_v
! -oneover_v
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Misc

IMPLICIT NONE

PRIVATE :: exp_v_32b, powr_v_32b, oneover_v_32b

INTERFACE exp_v_interface
  MODULE PROCEDURE exp_v, exp_v_32b
END INTERFACE

INTERFACE powr_v_interface
  MODULE PROCEDURE powr_v, powr_v_32b
END INTERFACE

INTERFACE oneover_v_interface
  MODULE PROCEDURE oneover_v, oneover_v_32b
END INTERFACE

CONTAINS

SUBROUTINE exp_v(n,x,y)
USE um_types
IMPLICIT NONE
! Sets y(i) to the exponential function of x(i), for i=1,..,n
INTEGER :: n
REAL (KIND=real64) :: y(n), x(n)
INTEGER (KIND=integer32) :: l_n
INTEGER :: i

#if defined(MKL_EXP_V)
  ! Interfaces for MKL
  INCLUDE 'mkl_vml.f90'
#endif

#if defined(VMASS)
  l_n=n
  CALL vexp (y, x, l_n)
#elif defined(MKL_EXP_V)
  CALL vdexp(n, x, y)
#else
  DO i=1, n
    y(i) = EXP(x(i))
  END DO
#endif
RETURN
END SUBROUTINE exp_v

SUBROUTINE exp_v_32b(n,x,y)
USE um_types
IMPLICIT NONE
! Sets y(i) to the exponential function of x(i), for i=1,..,n
INTEGER :: n
REAL (KIND=real32) :: y(n), x(n)
INTEGER (KIND=integer32) :: l_n
INTEGER :: i

#if defined(MKL_EXP_V)
  ! Interfaces for MKL
  INCLUDE 'mkl_vml.f90'
#endif

#if defined(VMASS)
  l_n=n
  CALL vsexp (y, x, l_n)
#elif defined(MKL_EXP_V)
  CALL vsexp(n, x, y)
#else
  DO i=1, n
    y(i) = EXP(x(i))
  END DO
#endif
RETURN
END SUBROUTINE exp_v_32b

!-----------------------------------------------------------

SUBROUTINE powr_v(n, x, power, z)
USE um_types
IMPLICIT NONE
! Sets z(i) to x(i) raised to the power power, for i=1,..,n
INTEGER :: n, i
REAL (KIND=real64) :: z(n), x(n), y(n), power
INTEGER (KIND=integer32) :: l_n

#if defined(MKL_POWR_V)
  ! Interfaces for MKL
  INCLUDE 'mkl_vml.f90'
#endif

#if defined(VMASS)
  l_n=n
  DO i=1, n
    y(i)=power
  END DO
  CALL vpow (z, x, y, l_n)
#elif defined(MKL_POWR_V)
  CALL vdpowx(n, x, power, z)
#else
  DO i=1, n
    z(i) = x(i)**power
  END DO
#endif
RETURN
END SUBROUTINE powr_v

SUBROUTINE powr_v_32b(n, x, power, z)
USE um_types
IMPLICIT NONE
! Sets z(i) to x(i) raised to the power power, for i=1,..,n
INTEGER :: n, i
REAL (KIND=real32) :: z(n), x(n), y(n), power
INTEGER (KIND=integer32) :: l_n

#if defined(MKL_POWR_V)
  ! Interfaces for MKL
  INCLUDE 'mkl_vml.f90'
#endif

#if defined(VMASS)
  l_n=n
  DO i=1, n
    y(i)=power
  END DO
  CALL vspow (z, x, y, l_n)
#elif defined(MKL_POWR_V)
  CALL vspowx(n, x, power, z)
#else
!Allowing this loop to vectorise with the Cray compiler version before
!8.4.0 and 32bits, breaks PROC comparability, even at safe
#if defined (CRAY_FORTRAN) && (CRAY_FORTRAN <8004000)
!DIR$ NOVECTOR
#endif
  DO i=1, n
    z(i) = x(i)**power
  END DO
#endif
RETURN
END SUBROUTINE powr_v_32b

!-----------------------------------------------------------

SUBROUTINE rtor_v(n, x, y, z)
USE um_types
IMPLICIT NONE

! Sets z(i) to x(i) raised to the power y(i), for i=1,..,n

INTEGER :: n
REAL (KIND=real64) :: z(n), x(n), y(n)
INTEGER (KIND=integer32) :: l_n
INTEGER :: i

#if defined(MKL_RTOR_V)
! Interfaces for MKL
INCLUDE 'mkl_vml.f90'
#endif

l_n=n

#if defined(VMASS)
CALL vpow (z, x, y, l_n)

#elif defined(MKL_RTOR_V)
CALL vdpow(n, x, y, z)

#else
DO i=1, n
  z(i) = x(i)**y(i)
END DO
#endif

RETURN
END SUBROUTINE rtor_v

!-----------------------------------------------------------

SUBROUTINE sqrt_v(n, x, y)
USE um_types
IMPLICIT NONE

! Sets y(i) to the square root of x(i), for i=1,..,n

INTEGER :: n
REAL (KIND=real64) :: x(n), y(n)
INTEGER (KIND=integer32) :: l_n
INTEGER :: i

#if defined(MKL_SQRT_V)
! Interfaces for MKL
INCLUDE 'mkl_vml.f90'
#endif

l_n=n


#if defined(VMASS)
CALL vsqrt (y, x, l_n)

#elif defined(MKL_SQRT_V)
CALL vdsqrt(n, x, y)

#else
DO i=1, n
  y(i) = SQRT(x(i))
END DO
#endif

RETURN
END SUBROUTINE sqrt_v

!-----------------------------------------------------------

SUBROUTINE oneover_v(n, x, y)
USE um_types
IMPLICIT NONE
! Sets y(i) to the reciprocal of x(i), for i=1,..,n
INTEGER :: n
REAL (KIND=real64) :: x(n), y(n)
INTEGER (KIND=integer32) :: l_n
INTEGER :: i

#if defined(VMASS)
  l_n=n
  CALL vrec (y, x, l_n)
#else
  DO i=1, n
    y(i) = 1/x(i)
  END DO
#endif
RETURN
END SUBROUTINE oneover_v

SUBROUTINE oneover_v_32b(n, x, y)
USE um_types
IMPLICIT NONE
! Sets y(i) to the reciprocal of x(i), for i=1,..,n
INTEGER :: n
REAL (KIND=real32) :: x(n), y(n)
INTEGER (KIND=integer32) :: l_n
INTEGER :: i

#if defined(VMASS)
  l_n=n
  CALL vsrec (y, x, l_n)
#else
  DO i=1, n
    y(i) = 1/x(i)
  END DO
#endif
RETURN
END SUBROUTINE oneover_v_32b

!-----------------------------------------------------------

SUBROUTINE log_v (n, x, y)
USE um_types
IMPLICIT NONE

! Sets y(i) to the natural logarithm of x(i), for i=1,..,n

INTEGER :: n
REAL (KIND=real64) :: x(n), y(n)
INTEGER (KIND=integer32) :: l_n
INTEGER :: i

#if defined(MKL_LOG_V)
! Interfaces for MKL
INCLUDE 'mkl_vml.f90'
#endif

l_n=n

#if defined(VMASS)
CALL vlog (y, x, l_n)

#elif defined(MKL_LOG_V)
CALL vdln( n, x, y )

#else
DO i=1, n
  y(i) = LOG(x(i))
END DO
#endif

RETURN
END SUBROUTINE log_v

!-----------------------------------------------------------

SUBROUTINE sin_v(n,x,y)
USE um_types
IMPLICIT NONE

! Sets y(i) to the sin function of x(i), for i=1,..,n

INTEGER :: n
REAL (KIND=real64) :: y(n), x(n)
INTEGER (KIND=integer32) :: l_n
INTEGER :: i

#if defined(MKL_SIN_V)
! Interfaces for MKL
INCLUDE 'mkl_vml.f90'
#endif

l_n=n

#if defined(VMASS)
CALL vsin (y, x, l_n)

#elif defined(MKL_SIN_V)
CALL vdsin(n, x, y)

#else
DO i=1, n
  y(i) = SIN(x(i))
END DO
#endif

RETURN
END SUBROUTINE sin_v

!-----------------------------------------------------------

SUBROUTINE cos_v(n,x,y)
USE um_types
IMPLICIT NONE

! Sets y(i) to the cos function of x(i), for i=1,..,n

INTEGER :: n
REAL (KIND=real64) :: y(n), x(n)
INTEGER (KIND=integer32) :: l_n
INTEGER :: i

#if defined(MKL_COS_V)
! Interfaces for MKL
INCLUDE 'mkl_vml.f90'
#endif

l_n=n

#if defined(VMASS)
CALL vcos (y, x, l_n)

#elif defined(MKL_COS_V)
CALL vdcos(n, x, y)

#else
DO i=1, n
  y(i) = COS(x(i))
END DO
#endif

RETURN
END SUBROUTINE cos_v

!-----------------------------------------------------------
SUBROUTINE acos_v(n,x,y)
USE um_types
IMPLICIT NONE

! Sets y(i) to the cos function of x(i), for i=1,..,n

INTEGER :: n
REAL (KIND=real64) :: y(n), x(n)
INTEGER (KIND=integer32) :: l_n
INTEGER :: i

#if defined(MKL_ACOS_V)
! Interfaces for MKL
INCLUDE 'mkl_vml.f90'
#endif

l_n=n

#if defined(VMASS)
CALL vacos (y, x, l_n)

#elif defined(MKL_ACOS_V)
CALL vdacos(n, x, y)

#else
DO i=1, n
  y(i) = ACOS(x(i))
END DO
#endif

RETURN
END SUBROUTINE acos_v

!-----------------------------------------------------------

SUBROUTINE asin_v(n,x,y)
USE um_types
IMPLICIT NONE

! Sets y(i) to the asin function of x(i), for i=1,..,n

INTEGER :: n
REAL (KIND=real64) :: y(n), x(n)
INTEGER (KIND=integer32) :: l_n
INTEGER :: i

#if defined(MKL_ASIN_V)
! Interfaces for MKL
INCLUDE 'mkl_vml.f90'
#endif

l_n=n

#if defined(VMASS)
CALL vasin (y, x, l_n)

#elif defined(MKL_ASIN_V)
CALL vdasin(n, x, y)

#else
DO i=1, n
  y(i) = ASIN(x(i))
END DO
#endif

RETURN
END SUBROUTINE asin_v

!-----------------------------------------------------------

SUBROUTINE atan2_v(n,a,b,y)
USE um_types
IMPLICIT NONE

! Sets y(i) to the atan2 function of a(i),b(i), for i=1,..,n

INTEGER :: n
REAL (KIND=real64) :: y(n), a(n),b(n)
INTEGER (KIND=integer32) :: l_n
INTEGER :: i

#if defined(MKL_ATAN2_V)
! Interfaces for MKL
INCLUDE 'mkl_vml.f90'
#endif

l_n=n

#if defined(VMASS)
CALL vatan2 (y, a, b, l_n)

#elif defined(MKL_ATAN2_V)
CALL vdatan2(n, a, b, y)

#else
DO i=1, n
  y(i) = ATAN2(a(i),b(i))
END DO
#endif

RETURN
END SUBROUTINE atan2_v

!-----------------------------------------------------------

SUBROUTINE cubrt_v(n, x, y)
USE um_types
IMPLICIT NONE

! Sets y(i) to the cube root of x(i), for i=1,..,n

INTEGER, INTENT(IN)  :: n
REAL (KIND=real64), INTENT(IN) :: x(n)
REAL (KIND=real64), INTENT(OUT) :: y(n)

! local variables
#if defined(VMASS)
INTEGER (KIND=integer32) :: l_n
#endif
INTEGER :: i

#if defined(MKL_CUBRT_V)
! Interfaces for MKL
INCLUDE 'mkl_vml.f90'
#endif


#if defined(VMASS)
!! use IBM's vector maths library if available
l_n=n
CALL vcbrt ( y, x, l_n)
#elif defined(MKL_CUBRT_V)
!! use INTEL's vector maths library if available
CALL vdcbrt(n, x, y)
#else
DO i=1, n
  y(i) = x(i)**(1.0/3.0)
END DO
#endif

RETURN
END SUBROUTINE cubrt_v

!-----------------------------------------------------------

END MODULE vectlib_mod
