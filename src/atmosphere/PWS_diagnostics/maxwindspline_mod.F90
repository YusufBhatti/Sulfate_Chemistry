! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Routine to calculate Max Wind value and height

MODULE MaxWindSpline_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE MaxWindSpline (NumLevs, i, j,    & 
                          UFields, VFields, PFields, MaxWLev,  &
                          Uinc, Vinc, Pinc )

! Description: Routine to calculate Max Wind value and height
!
! Method: as was done in fieldcalc.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS_diagnostics
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE atm_fields_bounds_mod

IMPLICIT NONE

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumLevs
INTEGER, INTENT(IN) :: i, j       ! Row, Col on grid


REAL, INTENT(IN) :: UFields(udims%i_start:udims%i_end,    &
                            vdims%j_start:vdims%j_end,    &
                            NumLevs) ! U-wind on B-grid

REAL, INTENT(IN) :: VFields(udims%i_start:udims%i_end,    &
                            vdims%j_start:vdims%j_end,    &
                            NumLevs) ! v-wind on B-grid

REAL, INTENT(IN) :: PFields(udims%i_start:udims%i_end,    &
                            vdims%j_start:vdims%j_end,    &
                            NumLevs) ! P on wind B-grid

INTEGER, INTENT(IN) :: MaxWLev(udims%i_start:udims%i_end, &
                               vdims%j_start:vdims%j_end)

INTEGER, PARAMETER :: ninc=16           ! number of increments used

REAL, INTENT(OUT) :: Uinc(2*ninc)  ! u at increment points
REAL, INTENT(OUT) :: Vinc(2*ninc)  ! v at increment points
REAL, INTENT(OUT) :: Pinc(2*ninc)  ! p at increment points

! Local constants

INTEGER, PARAMETER :: kmax=5            ! total levels needed for spline
INTEGER, PARAMETER :: khalf=(kmax+1)/2  ! number of levels each side

! Local variables:
INTEGER :: k
REAL    :: P_Lwr, P_Mid, P_Upr

! Variables involving spline interpolation
INTEGER :: spl(2*ninc)
REAL    :: xsp(kmax), ysp(kmax)
REAL    :: bsp(kmax), csp(kmax), dsp(kmax)
REAL    :: dx  (2*ninc)


! Use levels above and below to create splines of U and V
! Use the splines to evaluate U and V at intervals between
! level-1 and level+1      OLD : level-1/2 and level+1/2

! P_Lwr : P on Rho Level below maximum
! P_Mid : P on Rho Level at    maximum
! P_Upr : P on Rho Level above maximum


P_Lwr = PFields(i,j,(MaxWLev(i,j)-1))
P_Mid = PFields(i,j,MaxWLev(i,j)  )
P_Upr = PFields(i,j,(MaxWLev(i,j)+1))

DO k = 1, ninc                               ! Calc p increments
  Pinc(k)      = P_Lwr + (P_Mid-P_Lwr)*REAL(k)/REAL(ninc)
  Pinc(k+ninc) = P_Mid + (P_Upr-P_Mid)*REAL(k)/REAL(ninc)
END DO

!-- Pressure is x -- going down through levels so it is increasing
DO k = 1,kmax                      ! P on surrounding levels
  xsp(k) = PFields(i,j,(MaxWLev(i,j)+khalf-k) )
END DO

spl (1     :  ninc) = khalf
spl (1+ninc:2*ninc) = khalf-1
dx(:) = Pinc(:) - xsp(spl(:))

!-- Get U values --
DO k = 1,kmax                      ! U on surrounding levels
  ysp(k) = UFields(i,j,(MaxWLev(i,j)+khalf-k))
END DO

CALL SplineSetup( kmax, xsp, ysp, bsp, csp, dsp )
Uinc(:) = ysp(spl) + dx*( bsp(spl) + dx*( csp(spl) + dx*dsp(spl)))

!-- Get V values --
DO k = 1,kmax                      ! V on surrounding levels
  ysp(k) = VFields(i,j,(MaxWLev(i,j)+khalf-k))
END DO

CALL SplineSetup( kmax, xsp, ysp, bsp, csp, dsp )
Vinc(:) = ysp(spl) + dx*( bsp(spl) + dx*( csp(spl) + dx*dsp(spl)))

END SUBROUTINE MaxWindSpline


! ====================================================================


SUBROUTINE SplineSetup (n, x, y, b, c, d) 

IMPLICIT NONE

! Subroutine Arguments:
INTEGER, INTENT(IN)  :: n
REAL,    INTENT(IN)  :: x(n), y(n)
REAL,    INTENT(OUT) :: b(n), c(n), d(n)

! Local constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "SplineSetup"

! Local variables:
INTEGER :: i
REAL :: t
REAL :: deltax(n-1), deltay(n-1), dydx(n-1)
!-----------------------------------------------------------------------
!  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
!  for a cubic interpolating spline
!    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!    for  x(i) < x < x(i+1)
!  input..
!    n = the number of data points or knots (n >= 2)
!    x = the abscissas of the knots in strictly increasing order
!    y = the ordinates of the knots
!  output..
!    b, c, d  = arrays of spline coefficients as defined above.
!  using  p  to denote differentiation,
!    y(i) = s(x(i))
!    b(i) = sp(x(i))
!    c(i) = spp(x(i))/2
!    d(i) = sppp(x(i))/6  (derivative from the right)
!  the accompanying function subprogram SplineEval can be used
!  to evaluate the spline.
!-----------------------------------------------------------------------

IF ( n >= 2 ) THEN

  deltax(1:n-1) = x(2:n) - x(1:n-1)
  deltay(1:n-1) = y(2:n) - y(1:n-1)
  dydx  (1:n-1) = deltay(1:n-1)/deltax(1:n-1)

  IF ( n == 2 ) THEN
    b(1:2) = dydx(1)
    c(1:2) = 0.0
    d(1:2) = 0.0
  ELSE

    !--------------------------------------------------------------------
    ! Tridiagonal system:  b = diagonal, d = offdiagonal, c = RH side.
    ! Third derivs at x(1) and x(n) obtained from divided differences
    !--------------------------------------------------------------------
    b(1)     = -deltax(1)
    b(2:n-1) =  deltax(1:n-2) + deltax(2:n-1)
    b(n)     = -deltax(n-1)

    c(2:n-1) =  dydx(2:n-1)-dydx(1:n-2)
    IF ( n /= 3 ) THEN
      c(1) = c(  3)/b(  3) - c(  2)/b(  2)
      c(n) = c(n-2)/b(n-2) - c(n-1)/b(n-1)
      c(1) = c(1)*b(1)**2/(b(  3)-b(1))
      c(n) = c(n)*b(n)**2/(b(n-2)-b(n))
    ELSE
      c(1) = 0.0
      c(n) = 0.0
    END IF
    b(2:n-1) = 2.0*b(2:n-1)
    !----------------------------------------------  forward elimination
    DO i = 2, n
      t=deltax(i-1)/b(i-1)
      b(i) = b(i) - t*deltax(i-1)
      c(i) = c(i) - t*c(i-1)
    END DO
    !------------------------------------------------  back substitution
    c(n) = c(n)/b(n)
    DO i = n-1, 1, -1
      c(i) = (c(i) - deltax(i)*c(i+1))/b(i) !2nd derivative/ 6
    END DO
    !------------------------------ c(i) is now the sigma(i) of the text
    !---------------------------------------- compute polynomial coeff.s
    !-----------------------------------------------------------------------
    d(1:n-1) =(c(2:n)-c(1:n-1))/deltax(1:n-1) !3rd derivative/ 6
    d(n)     = d(n-1)                         !

    b(1:n-1) = dydx(1:n-1) - deltax(1:n-1)*(c(2:n) + 2.0*c(1:n-1))
    b(n)     = dydx(n-1)   + deltax(n-1)  *(c(n-1) + 2.0*c(n)    )

    c(1:n)   = 3.0*c(1:n)  ! 2nd derivative/2

  END IF
END IF

END SUBROUTINE SplineSetup

!=======================================================================

END MODULE MaxWindSpline_mod
