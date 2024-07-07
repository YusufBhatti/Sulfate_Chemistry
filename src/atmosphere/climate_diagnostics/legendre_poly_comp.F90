! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! This routine computes the Associate Legendre Polynomials and populates
! a matrix Ymn with their values for different latitudes, meridional
! wavenumber 'l' and zonal wavenumber 'm'.
! It employs recurrence relationships from Y(0,0,lat) to fill the
! matrix up to Ntop (maximum wavenumber)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Climate Diagnostics
MODULE legendre_poly_comp_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LEGENDRE_POLY_COMP_MOD'

CONTAINS

SUBROUTINE legendre_poly_comp(global_rows,delta_phi)

USE track_mod,          ONLY: ntop_850,ntop_tc,Ymn
USE conversions_mod,    ONLY: pi
USE yomhook,            ONLY: lhook, dr_hook
USE parkind1,           ONLY: jprb, jpim

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: global_rows   ! Number of latitude points
REAL,    INTENT(IN) :: delta_phi     ! Grid latitude spacing in radians

! Local variables
INTEGER :: ntop    ! Truncation (n1~lower and n2~upper) limit.
INTEGER :: l,m,j   ! Indexes for DO loops
INTEGER :: lp1,lm1 ! l+1 and l-1 for indexing arrays

REAL ::  sin_v(global_rows+1) ! array for sine to build up Fm matrix
REAL ::  cos_v(global_rows+1) ! array for cosine to build up Fm matrix
REAL ::  SH_fac               ! Factor for SH before j (lat loop)
REAL ::  SH_lm1               ! Factor for the Y(l-1,m) term

CHARACTER(LEN=*), PARAMETER  :: RoutineName='LEGENDRE_POLY_COMP'

! DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! End of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Create sin_v: Use Cos (lat-Pi/2) to remap the latitude
DO j=1,global_rows+1
  sin_v(j)=SIN(delta_phi*REAL(j-1))
END DO

! Create sin_v -> mind the '-' in the sin(phi -Pi/2) = - cos(phi)
DO j=1,global_rows+1
  cos_v(j)=COS( delta_phi*REAL(j-1))
END DO

! Declaration of Ymn matrix. Set maximun wavenumber and allocate Ymn
ntop = MAX(ntop_850,ntop_tc)
IF (.NOT. ALLOCATED(Ymn)) ALLOCATE(Ymn( global_rows+1,0:ntop,0:ntop) )

! Populate matrix at m=0, n=0
! Compute P(0,0) = 1. / sqrt (4 Pi) for all latitudes
!
! Example of values computed for Ymn in this loop:
! * denotes values computed in this loop and x values of Ymn
! | *         |
! | x x       |
! | x x x     |
! | x x x x   |
! | x x x x x |

!Amplitude factor
SH_fac = 0.5 / SQRT(pi)

! Loop over all latitudes
DO j=1,global_rows+1
  Ymn(j,0,0) = SH_fac
END DO

! Populate diagonal (l=m) using recurrence relationship
! where Y(l+1,l+1) depends on Y(l,l). For l=0, ... , ntop-1
!
! Example of values computed for Ymn in this loop:
! * denotes values computed in this loop and x values of Ymn
! | x         |
! | x *       |
! | x x *     |
! | x x x *   |
! | x x x x * |

DO l=0,ntop-1
  lp1=l+1
  ! Compute amplitude factor
  SH_fac = (-1.0)*SQRT( REAL(2*l + 3)/REAL(2*l+2) )
  ! Do loop over latitude for Y(l+1,l+1)
  DO j=1,global_rows+1
    Ymn(j,lp1,lp1) = SH_fac * sin_v(j) * Ymn(j,l,l)
  END DO
END DO

! Populate lower bi-diagonal using recurrence relationship
! where Y(l+1,l) depends on Y(l,l).
! For l=1,..., ntop and m=0, ... , ntop-1.
!
! Example of values computed for Ymn in this loop:
! * denotes values computed in this loop and x values of Ymn
! | x         |
! | * x       |
! | x * x     |
! | x x * x   |
! | x x x * x |

DO l=0,ntop-1
  lp1=l+1
  ! Factor for Y(l,l)
  SH_fac = SQRT( REAL(2*l + 3))
  ! Compute all latitudes
  DO j=1,global_rows+1
    Ymn(j,lp1,l) = SH_fac * cos_v(j) * Ymn(j,l,l)
  END DO
END DO

! Populate firts column (where l=0) using Legendre Polynomials
! recurrence relationship where Y(l+1,0) depends on Y(l,0) and Y(l-1,0).
! For n=1, ..., ntop-1 and m=0
!
! Example of values computed for Ymn in this loop:
! * denotes values computed in this loop and x values of Ymn
! | x         |
! | x x       |
! | * x x     |
! | * x x x   |
! | * x x x x |

DO l=1,ntop-1
  lp1=l+1
  lm1=l-1
  ! Compute amplitude factor for Y(l,0)
  SH_fac = 1.0 / REAL(l+1) * SQRT( REAL((2*l + 3)*(2*l+1)) )
  ! Compute amplitude factor for Y(l-1,0)
  SH_lm1 = (-1.0)*REAL(l) / REAL(l+1) * SQRT( REAL(2*l+3) / REAL(2*l - 1) )

  ! Compute Y(l+1,0) for all latitudes
  DO j=1,global_rows+1
    Ymn(j,lp1,0) = SH_fac * cos_v(j)*Ymn(j,l,0) + SH_lm1 * Ymn(j,lm1,0)
  END DO
END DO

! Populate matrix for columns using recurrence relationship where
! Y(l+1,m) depends on Y(l,m) and Y(l-1,m). For each m (from 1 to ntop-2)
! a loop is performend from l+1 to ntop-1
!
!
! Example of values computed for Ymn in this loop:
! * denotes values computed in this loop and x values of Ymn
! | x         |
! | x x       |
! | x x x     |
! | x * x x   |
! | x * * x x |

DO m=1,ntop-2

  ! Populate values of each column from l=m+1 to ntop-1 (mind the l+1
  ! of the lhs)
  DO l=m+1,ntop-1
    lp1= l + 1
    lm1= l - 1

    ! Amplitude factor for Y(l,m)
    SH_fac = SQRT( (4.0* REAL((l+1)*(l+1)) -1.0)  /                         &
                    (REAL((l+1)*(l+1)) - REAL(m*m) ) )

    ! Amplitude factor for Y(l-1,m)
    SH_lm1 = SH_fac *                                                     &
             SQRT( (REAL(l*l) - REAL(m*m)) / (4.0* REAL(l*l) -1.0) )

    ! Compute Y(l+1,m) for each latitude
    DO j=1,global_rows+1
      Ymn(j,lp1,m) = SH_fac* cos_v(j) * Ymn(j,l,m) -                      &
                               SH_lm1 * Ymn(j,lm1,m)
    END DO

  END DO
END DO

! Set Ymn to zero for m > n above diagonal.
!
! Example of values computed for Ymn in this loop:
! * denotes values computed in this loop and x values of Ymn
! | x * * * * |
! | x x * * * |
! | x x x * * |
! | x x x x * |
! | x x x x x |

DO m=1,ntop
  DO l=0,m-1
    DO j=1,global_rows+1
      Ymn(j,l,m) = 0.0
    END DO
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE legendre_poly_comp

END MODULE legendre_poly_comp_mod
