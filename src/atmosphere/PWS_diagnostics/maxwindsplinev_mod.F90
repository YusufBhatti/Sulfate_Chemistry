! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE MaxWindSplineV_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE MaxWindSplineV(NumCols, NumRows, NumLevs,           &
                          UFields, VFields, PFields, MaxWLev,  &
                          Uinc, Vinc, Pinc )
! Description: Routine to calculate Max Wind value and height
!
! Method: as was used in fieldcalc.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS diagnostics
!
! Code Description:
!   Language:           Fortran 90


IMPLICIT NONE

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumCols, NumRows, NumLevs

REAL, INTENT(IN) :: UFields(NumCols,NumRows,NumLevs) ! U-wind on B-grid
REAL, INTENT(IN) :: VFields(NumCols,NumRows,NumLevs) ! V-wind on B-grid
REAL, INTENT(IN) :: PFields(NumCols,NumRows,NumLevs) ! P on wind B-grid

INTEGER, INTENT(IN) :: MaxWLev(NumCols,NumRows)

INTEGER, PARAMETER :: ninc=16           ! number of increments used

REAL, INTENT(OUT) :: Uinc(NumCols,NumRows,2*ninc)  ! u at increment pts
REAL, INTENT(OUT) :: Vinc(NumCols,NumRows,2*ninc)  ! v at increment pts
REAL, INTENT(OUT) :: Pinc(NumCols,NumRows,2*ninc)  ! p at increment pts

! Local constants

INTEGER, PARAMETER :: kmax=5            ! total levels needed for spline
INTEGER, PARAMETER :: khalf=(kmax+1)/2  ! number of levels each side

! Local variables:
INTEGER :: i, j, k, kint
REAL :: P_Lwr(NumCols,NumRows)
REAL :: P_Mid(NumCols,NumRows) 
REAL :: P_Upr(NumCols,NumRows) 

#if defined(VECTOR)
REAL :: Arr_3D(NumCols,NumRows,NumLevs)
#endif

! Variables involving spline interpolation
INTEGER :: spl(2*ninc)
REAL :: xsp(NumCols,NumRows,kmax)
REAL :: ysp(NumCols,NumRows,kmax)
REAL :: bsp(NumCols,NumRows,kmax)
REAL :: csp(NumCols,NumRows,kmax)
REAL :: dsp(NumCols,NumRows,kmax)
REAL :: dx(NumCols,NumRows,2*ninc)


! Use levels above and below to create splines of U and V
! Use the splines to evaluate U and V at intervals between
! level-1 and level+1      OLD : level-1/2 and level+1/2

! Since MaxWLev may be 0 on places where U or V is missing
! we have to define an array which has save values there


#if defined(VECTOR)
! Store PFields into a conventional 3D array for vectorization
DO k=1,NumLevs
  Arr_3D(:,:,k) = PFields(:,:,k)
END DO
#endif

! P_Lwr : P on Rho Level below maximum
! P_Mid : P on Rho Level at    maximum
! P_Upr : P on Rho Level above maximum

DO j=1,NumRows
  DO i=1,NumCols
    IF (maxwlev(i,j) /= 0) THEN
#if defined(VECTOR)
      P_Lwr(i,j) = Arr_3D(i,j,MaxWLev(i,j)-1)
      P_Mid(i,j) = Arr_3D(i,j,MaxWLev(i,j)  )
      P_Upr(i,j) = Arr_3D(i,j,MaxWLev(i,j)+1)
#else
      P_Lwr(i,j) = PFields(i,j,(MaxWLev(i,j)-1))
      P_Mid(i,j) = PFields(i,j,MaxWLev(i,j)) 
      P_Upr(i,j) = PFields(i,j,(MaxWLev(i,j)+1)) 
#endif
    END IF
  END DO
END DO

DO k = 1, ninc                               ! Calc p increments
  DO j=1,NumRows
    DO i=1,NumCols
      IF (maxwlev(i,j) /= 0) THEN
        Pinc(i,j,k)      = P_Lwr(i,j) + &
                    (P_Mid(i,j)-P_Lwr(i,j))*REAL(k)/REAL(ninc)
        Pinc(i,j,k+ninc) = P_Mid(i,j) + &
                    (P_Upr(i,j)-P_Mid(i,j))*REAL(k)/REAL(ninc)
      END IF
    END DO
  END DO
END DO

!-- Pressure is x -- going down through levels so it is increasing
DO k = 1,kmax                      ! P on surrounding levels
  DO j=1,NumRows
    DO i=1,NumCols
      IF (maxwlev(i,j) /= 0) THEN
#if defined(VECTOR)
        xsp(i,j,k) = Arr_3D(i,j,MaxWLev(i,j)+khalf-k)
#else
        xsp(i,j,k) = PFields(i,j,(MaxWLev(i,j)+khalf-k)) 
#endif
      END IF
    END DO
  END DO
END DO


spl (1     :  ninc) = khalf
spl (1+ninc:2*ninc) = khalf-1


DO k=1,2*ninc
  kint = spl(k)
  DO j=1,NumRows
    DO i=1,NumCols
      IF (maxwlev(i,j) /= 0) THEN
        dx(i,j,k) = Pinc(i,j,k) - xsp(i,j,kint)
      END IF
    END DO
  END DO
END DO


#if defined(VECTOR)
!-- Store UFields into a conventional 3D array for vectorization
DO k=1,NumLevs
  Arr_3D(:,:,k) = UFields(:,:,k) 
END DO
#endif

!-- Get U values --
DO k = 1,kmax                      ! U on surrounding levels
  DO j=1,NumRows
    DO i=1,NumCols
      IF (maxwlev(i,j) /= 0) THEN
#if defined(VECTOR)
        ysp(i,j,k) = Arr_3D(i,j,MaxWLev(i,j)+khalf-k)
#else
        ysp(i,j,k) = UFields(i,j,(MaxWLev(i,j)+khalf-k))
#endif
      END IF
    END DO
  END DO
END DO

CALL SplineSetupV(NumCols, NumRows, kmax, xsp, ysp, bsp, csp, &
                   dsp, MaxWLev )

DO k=1,2*ninc
  kint = spl(k)
  DO j=1,NumRows
    DO i=1,NumCols
      IF (maxwlev(i,j) /= 0) THEN
        Uinc(i,j,k) = ysp(i,j,kint) + dx(i,j,k)*( bsp(i,j,kint) + &
              dx(i,j,k)*(csp(i,j,kint) + dx(i,j,k)*dsp(i,j,kint)))
      END IF
    END DO
  END DO
END DO

#if defined(VECTOR)
!-- Store VFields into a conventional 3D array for vectorization
DO k=1,NumLevs
  Arr_3D(:,:,k) = VFields(:,:,k) 
END DO
#endif

!-- Get V values --
DO k = 1,kmax                      ! V on surrounding levels
  DO j=1,NumRows
    DO i=1,NumCols
      IF (maxwlev(i,j) /= 0) THEN
#if defined(VECTOR)
        ysp(i,j,k) = Arr_3D(i,j,MaxWLev(i,j)+khalf-k)
#else
        ysp(i,j,k) = VFields(i,j,(MaxWLev(i,j)+khalf-k)) 
#endif
      END IF
    END DO
  END DO
END DO


CALL SplineSetupV(NumCols, NumRows, kmax, xsp, ysp, bsp, csp, &
                  dsp, MaxWLev )

DO k=1,2*ninc
  kint = spl(k)
  DO j=1,NumRows
    DO i=1,NumCols
      IF (maxwlev(i,j) /= 0) THEN
        Vinc(i,j,k) = ysp(i,j,kint) + dx(i,j,k)*( bsp(i,j,kint) + &
        dx(i,j,k)*( csp(i,j,kint) + dx(i,j,k)*dsp(i,j,kint)))
      END IF
    END DO
  END DO
END DO

END SUBROUTINE MaxWindSplineV


SUBROUTINE SplineSetupV(m, n, kmax, x, y, b, c, d, MaxWLev) 

IMPLICIT NONE

! Subroutine Arguments:
INTEGER, INTENT(IN) :: m, n,kmax
REAL,    INTENT(IN) :: x(m,n,kmax), y(m,n,kmax)
REAL,    INTENT(OUT) :: b(m,n,kmax), c(m,n,kmax), d(m,n,kmax)

INTEGER, INTENT(IN) :: MaxWLev(m,n)

! Local constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "SplineSetupV"

! Local variables:
INTEGER :: i, j, k
REAL :: t(m,n)
REAL :: deltax(m,n,kmax-1), deltay(m,n,kmax-1), dydx(m,n,kmax-1)

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

IF ( kmax >= 2 ) THEN

  DO k=1, kmax-1
    DO j=1, n
      DO i=1, m
        IF (maxwlev(i,j) /=0 ) THEN
          deltax(i,j,k)=x(i,j,k+1)-x(i,j,k)
          deltay(i,j,k)=y(i,j,k+1)-y(i,j,k)
          dydx(i,j,k)=deltay(i,j,k)/deltax(i,j,k)
        END IF
      END DO
    END DO
  END DO

  IF ( kmax == 2 ) THEN
    DO j=1, n
      DO i=1, m
        IF (maxwlev(i,j) /=0 ) THEN
          b(i,j,1) = dydx(i,j,1)
          b(i,j,2) = dydx(i,j,1)
          c(i,j,1:2) = 0.0
          d(i,j,1:2) = 0.0
        END IF
      END DO
    END DO

  ELSE
    !--------------------------------------------------------------------
    ! Tridiagonal system:  b = diagonal, d = offdiagonal, c = RH side.
    ! Third derivs at x(1) and x(n) obtained from divided differences
    !--------------------------------------------------------------------
    DO j=1, n
      DO i=1, m
        IF (maxwlev(i,j) /=0 ) THEN
          b(i,j,1)     = -deltax(i,j,1)
          b(i,j,kmax)  = -deltax(i,j,kmax-1)
        END IF
      END DO
    END DO

    DO k=2, kmax-1
      DO j=1, n
        DO i=1, m
          IF (maxwlev(i,j) /=0 ) THEN
            b(i,j,k) =  deltax(i,j,k-1) + deltax(i,j,k)
            c(i,j,k) =  dydx(i,j,k)-dydx(i,j,k-1)
          END IF
        END DO
      END DO
    END DO

    DO j=1, n
      DO i=1, m
        IF (maxwlev(i,j) /=0 ) THEN
          IF ( kmax /= 3 ) THEN
            c(i,j,1) = c(i,j,3)/b(i,j,3) - c(i,j,2)/b(i,j,2)
            c(i,j,kmax) = c(i,j,kmax-2)/b(i,j,kmax-2) - &
                          c(i,j,kmax-1)/b(i,j,kmax-1)
            c(i,j,1) = c(i,j,1)*b(i,j,1)**2/(b(i,j,3)-b(i,j,1))
            c(i,j,kmax) = c(i,j,kmax)*b(i,j,kmax)**2 &
                          /(b(i,j,kmax-2)-b(i,j,kmax))
          ELSE
            c(i,j,1) = 0.0
            c(i,j,n) = 0.0
          END IF
        END IF
      END DO
    END DO

    DO k=2, kmax-1
      DO j=1, n
        DO i=1, m
          IF (maxwlev(i,j) /=0 ) THEN
            b(i,j,k) = 2.0*b(i,j,k)
          END IF
        END DO
      END DO
    END DO

    !----------------------------------------------  forward elimination
    DO k = 2, kmax
      DO j=1, n
        DO i=1, m
          IF (maxwlev(i,j) /=0 ) THEN
            t(i,j)=deltax(i,j,k-1)/b(i,j,k-1)
            b(i,j,k) = b(i,j,k) - t(i,j)*deltax(i,j,k-1)
            c(i,j,k) = c(i,j,k) - t(i,j)*c(i,j,k-1)
          END IF
        END DO
      END DO
    END DO

    !------------------------------------------------  back substitution

    DO j=1, n
      DO i=1, m
        IF (maxwlev(i,j) /=0 ) THEN
          c(i,j,kmax) = c(i,j,kmax)/b(i,j,kmax)
        END IF
      END DO
    END DO

    DO k = kmax-1, 1, -1
      DO j=1, n
        DO i=1, m
          IF (maxwlev(i,j) /=0 ) THEN
            c(i,j,k) = (c(i,j,k) - deltax(i,j,k)*c(i,j,k+1)) &
                        /b(i,j,k) !2nd derivative/ 6
          END IF
        END DO
      END DO
    END DO

    !------------------------------ c(:,:,i) is now the sigma(i) of the text
    !---------------------------------------- compute polynomial coeff.s
    !-----------------------------------------------------------------------
    DO k=1, kmax-1
      DO j=1, n
        DO i=1, m
          IF (maxwlev(i,j) /=0 ) THEN
            d(i,j,k) =(c(i,j,k+1)-c(i,j,k)) &
                       /deltax(i,j,k) !3rd derivative/ 6
          END IF
        END DO
      END DO
    END DO

    DO j=1, n
      DO i=1, m
        IF (maxwlev(i,j) /=0 ) THEN
          d(i,j,kmax)     = d(i,j,kmax-1)                         !
        END IF
      END DO
    END DO

    DO k=1, kmax-1
      DO j=1, n
        DO i=1, m
          IF (maxwlev(i,j) /=0 ) THEN
            b(i,j,k)=dydx(i,j,k)-deltax(i,j,k)* &
                      (c(i,j,k+1)+2.0*c(i,j,k))
          END IF
        END DO
      END DO
    END DO

    DO j=1, n
      DO i=1, m
        IF (maxwlev(i,j) /=0 ) THEN
          b(i,j,kmax)=dydx(i,j,kmax-1)  +deltax(i,j,kmax-1) &
            *(c(i,j,kmax-1)+2.0*c(i,j,kmax))
        END IF
      END DO
    END DO

    DO k=1,kmax
      DO j=1, n
        DO i=1, m
          IF (maxwlev(i,j) /=0 ) THEN
            c(i,j,k)   = 3.0*c(i,j,k)  ! 2nd derivative/2
          END IF
        END DO
      END DO
    END DO

  END IF
END IF

END SUBROUTINE SplineSetupV

!=======================================================================


END MODULE MaxWindSplineV_mod
