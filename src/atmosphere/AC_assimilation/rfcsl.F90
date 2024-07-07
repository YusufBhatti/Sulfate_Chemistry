! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation

!     SUBROUTINE RFCSL
!     RECURSIVE FILTER - CONSTANT COEFF - FOR A SCALAR LTD.AREA FIELD
!
!     THIS VERSION FOR A VARIABLE RESOLUTION GRID AND DLAT AND DLONG
!     ARE VECTORS OF THE GRID SPACING
!
!     ARGUMENTS:
MODULE rfcsl_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RFCSL_MOD'

CONTAINS

SUBROUTINE rfcsl(f,n_rows,row_length,m_grid,f_boundary,           &
 coslat,dlat,dlong,SCALE,npass                                    &
  )

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
!$    USE omp_lib
IMPLICIT NONE
!---- INPUT
INTEGER :: n_rows      ! NUMBER OF ROWS IN GRID
INTEGER :: row_length  ! NUMBER OF POINTS IN EACH ROW
INTEGER :: m_grid      ! 0=USE F_BOUNDARY, -1=NO SPECIAL BOUNDARY
REAL :: f_boundary     ! F ASSUMED OUTSIDE LTD.AREA
REAL :: coslat(n_rows) ! COSINE(LATITUDE)
REAL :: dlat(n_rows)             ! ROW SPACING (RADIANS)
REAL :: dlong(row_length)        ! POINT SPACING (RADIANS)
REAL :: SCALE                    ! FILTER SCALE (RADIANS)
INTEGER :: npass       ! NUMBER OF PASSES IN EACH DIRECTION
!             NPASS=2 WILL GIVE EQUIVALENT OF SOAR FUNCTION (USUAL)
!             NPASS=LARGE WILL GIVE EQUIVALENT OF GAUSSIAN FUNCTION
!---- IN-OUT
REAL :: f(row_length,n_rows)     ! FIELD TO BE FILTERED
!
!---- WORK SPACE
REAL :: a_ns(n_rows)                    ! FILTER COEFFICIENTS N-S-N
REAL :: a_we(row_length,n_rows)         ! W-E-W FILTER COEFFICIENTS
REAL :: e     ! INTERMEDIATE VALUE IN FORMULA FOR FILTER COEFFICIENT
REAL :: z     ! WORK VARIABLE FOR CALCULATING COEFFS
INTEGER :: jpass,jrow,jpt     ! LOOP COUNTERS FOR PASS,ROW,PT IN ROW
INTEGER :: jj, omp_block      ! FOR OMP LOOP BLOCKING

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RFCSL'

!
! ----------------------------------------------------------------------
!           PROGRAMMING NOTES
!     THE ORDER OF NESTING OF DOUBLE LOOPS ON JPT & JROW IS OPTIONAL.
!     SUCH LOOPS ARE INDICATED BY !)) ON THE ENDDO STATEMENT.
!     THE LOOP IN THE DIRECTION OF FILTERING IS RECURSIVE, AND SO CANNOT
!     BE VECTORIZED; HERE IT IS MADE THE OUTER LOOP.
!     A SCALAR COMPUTER MAY BE MORE EFFICIENT WITH A RECURSIVE INNER
!-----------------------------------------------------------------------
! **        BRIEF DESCRIPTION OF METHOD
!     THE INPUT FIELD IS F
!     FIRST PASS IN ONE DIRECTION GIVES G (OVERWRITING ARRAY F):
!      G(I)=A(I)*G(I-1)+B(I)*F(I)  :  I=1,N
!     PASS IN OPPOSITE DIRECTION GIVES H (OVERWRITING ARRAY F):
!      H(I)=A(I)*H(I+1)+B(I)*G(I)  :  I=N,1
!     FOR FILTER ALONG ROW (CALLED E-W-E IN COMMENTS) THE GRID IS
!      REGULAR, SO B=1-A, AND A IS PRECALCULATED.
!     FOR FILTER OF LIMITED AREA THE GRID IS
!      ASSUMED REGULAR, SO B=1-A, AND A IS PRECALCULATED.
!      BECAUSE OF THIS APPROX, THE GRID SHOULD NOT GO NEAR POLE.
!-----------------------------------------------------------------------
!     NO EXTERNAL SUBROUTINE CALLS
!-----------------------------------------------------------------------
!  ** PRECALCULATE COEFFICIENTS
!
!---- CALCULATE COEFFS FOR E->W SWEEP (ALSO USED FOR W->E)
!       USING 6.5.13 & 6.5.9
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP  PARALLEL DEFAULT(NONE) SHARED(row_length, n_rows, coslat, dlong,  &
!$OMP& a_we, a_ns, dlat, npass, scale) PRIVATE(z, e, jrow, jpt)

!$OMP DO SCHEDULE(STATIC)
DO jpt=1,row_length
  DO jrow=1,n_rows
    z=npass*(dlong(jpt)*coslat(jrow))**2*0.25
    e=z/SCALE          **2
    a_we(jpt,jrow)= 1.0+e-SQRT(e*(e+2.0))
  END DO ! JROW
END DO   ! JPT
!$OMP END DO

!
!---- CALCULATE COEFFS FOR N->S SWEEP (ALSO USED FOR S->N)
!       USING 6.5.13 & 6.5.9

!$OMP DO SCHEDULE(STATIC)
DO jrow=1,n_rows
  z=npass*dlat(jrow)**2*0.25
  e=z/SCALE          **2
  a_ns(jrow)      = 1.0+e-SQRT(e*(e+2.0))
END DO ! JROW
!$OMP END DO

!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(jrow, jpt, jpass, omp_block, jj)

!-----------------------------------------------------------------------
!
!  ** START LOOP OVER PASSES
DO jpass=1,npass     ! ===========================================
  !
  !  **  FILTER W->E              -----------------------------------
  !



  IF (m_grid == 0) THEN! BCS SPECIFIED TO ALLOW FOR CONSTANT OUTSIDE
    !      DO W BC
    !      THE BOUNDARY CONDITIONS DEPEND ON THE NUMBER OF PREVIOUS PASSES
    !      IN THE OPPOSITE DIRECTION.
    !      IF F_BOUNDARY=0 THE FORMULAE FOR 0,1,2 PREVIOUS PASSES ARE:-
    !        G1=(1-A)F1                          (6.5.26)
    !        G1=(1/(1+A))F1                      (6.5.27)
    !        G1=(1/(1+A)(1-A**2))(F1-A**3*F2)    (6.5.28)
    !      IF F_BOUNDARY /= 0 IT MUST BE SUBTRACTED BEFORE USING FORMULAE.
    IF (jpass == 1) THEN    !     NO PREVIOUS E->W PASS

!$OMP DO SCHEDULE(STATIC)
      DO      jrow=1,n_rows
        f(1,jrow)=f_boundary+(f(1,jrow)-f_boundary)                   &
        *(1.0-a_we(1,jrow))
      END DO ! JROW
!$OMP END DO

    ELSE IF (jpass == 2) THEN    !    ONE PREVIOUS E->W PASS

!$OMP DO SCHEDULE(STATIC)
      DO      jrow=1,n_rows
        f(1,jrow)=f_boundary+(f(1,jrow)-f_boundary)                   &
                      /(1.0+a_we(1,jrow))
      END DO ! JROW
!$OMP END DO
    ELSE IF (jpass >= 3) THEN    !    TWO PREVIOUS E->W PASSES
      !                                       BC FOR 2 IS ALSO USED FOR >2

!$OMP DO SCHEDULE(STATIC)
      DO      jrow=1,n_rows
        f(1,jrow)=f_boundary+(f(1,jrow)-f_boundary-                   &
             a_we(1,jrow)**3*(f(2,jrow)-f_boundary))                  &
         /((1.0+a_we(1,jrow))*(1.0-a_we(1,jrow)**2))
      END DO ! JROW
!$OMP END DO

    END IF ! JPASS
    !      ELSE  ! M_GRID /= 0  NO SPECIAL TREATMENT OF BOUNDARY
  END IF ! M_GRID
  !
  !          FILTER W->E PASS      EQN 6.5.1

  omp_block = n_rows
  !      ENSURE AS LARGE A STRIDE AS POSSIBLE FOR EACH OMP THREAD
!$     OMP_BLOCK = CEILING(REAL(N_ROWS)/OMP_GET_NUM_THREADS())

!$OMP DO SCHEDULE(STATIC)
  DO        jj=1, n_rows, omp_block
    DO        jpt=2,row_length
      DO        jrow=jj,MIN(jj+omp_block-1,n_rows)
        f(jpt,jrow)=f(jpt,jrow)+                                       &
            (f(jpt-1,jrow)-f(jpt,jrow))*a_we(jpt,jrow)
      END DO !)) JROW
    END DO !)) JPT
  END DO !)) JJ
!$OMP END DO

  !
  !  **  FILTER E->W              -----------------------------------
  !
  IF (m_grid == 0) THEN! BCS SPECIFIED TO ALLOW FOR CONSTANT OUTSIDE
    !      DO E BC
    IF (jpass == 1) THEN    !    ONE PREVIOUS W->E PASS

!$OMP DO SCHEDULE(STATIC)
      DO      jrow=1,n_rows
        f(row_length,jrow)=f_boundary+(f(row_length,jrow)-f_boundary) &
                               /(1.0+a_we(row_length,jrow))
      END DO ! JROW
!$OMP END DO
    ELSE IF (jpass >= 2) THEN    !    TWO PREVIOUS W->E PASSES
      !                                       BC FOR 2 IS ALSO USED FOR >2
!$OMP DO SCHEDULE(STATIC)
      DO      jrow=1,n_rows
        f(row_length,jrow)=f_boundary+(f(row_length,jrow)-f_boundary- &
         a_we(row_length,jrow)**3*(f(row_length-1,jrow)-f_boundary))  &
         /((1.0+a_we(row_length,jrow))*(1.0-a_we(row_length,jrow)**2))
      END DO ! JROW
!$OMP END DO

    END IF ! JPASS
    !      ELSE  ! M_GRID /= 0  NO SPECIAL TREATMENT OF BOUNDARY
  END IF ! M_GRID
  !
  !          FILTER E->W PASS

  omp_block = n_rows
  !      ENSURE AS LARGE A STRIDE AS POSSIBLE FOR EACH OMP THREAD
!$     OMP_BLOCK = CEILING(REAL(N_ROWS)/OMP_GET_NUM_THREADS())

!$OMP DO SCHEDULE(STATIC)
  DO        jj=1,n_rows,omp_block
    DO        jpt=row_length-1,1,-1
      DO        jrow=jj,MIN(jj+omp_block-1, n_rows)
        f(jpt,jrow)=f(jpt,jrow)+                                       &
            (f(jpt+1,jrow)-f(jpt,jrow))*a_we(jpt,jrow)
      END DO !)) JROW
    END DO !)) JPT
  END DO !)) JJ
!$OMP END DO

  !
  !  **   FILTER N->S             -----------------------------------
  !
  IF (m_grid == 0) THEN! BCS SPECIFIED TO ALLOW FOR CONSTANT OUTSIDE
    !      DO N BC
    IF (jpass == 1) THEN    !     NO PREVIOUS S->N PASS

!$OMP DO SCHEDULE(STATIC)
      DO      jpt=1,row_length
        f(jpt,1)=f_boundary+(f(jpt,1)-f_boundary)                     &
         *(1.0-a_ns(1))
      END DO ! JPT
!$OMP END DO

    ELSE IF (jpass == 2) THEN    !    ONE PREVIOUS S->N PASS

!$OMP DO SCHEDULE(STATIC)
      DO      jpt=1,row_length
        f(jpt,1)=f_boundary+(f(jpt,1)-f_boundary)                     &
                     /(1.0+a_ns(1))
      END DO ! JPT
!$OMP END DO
    ELSE IF (jpass >= 3) THEN    !    TWO PREVIOUS S->N PASSES
      !                                       BC FOR 2 IS ALSO USED FOR >2
!$OMP DO SCHEDULE(STATIC)
      DO      jpt=1,row_length
        f(jpt,1)=f_boundary+(f(jpt,1)-f_boundary-                     &
         a_ns(1)    **3*(f(jpt,2)-f_boundary))                        &
         /((1.0+a_ns(1))*(1.0-a_ns(1)**2))
      END DO ! JPT
!$OMP END DO
    END IF ! JPASS
    !      ELSE  ! M_GRID /= 0  NO SPECIAL TREATMENT OF BOUNDARY
  END IF ! M_GRID
  !

  omp_block = row_length
  !      ENSURE AS LARGE A STRIDE AS POSSIBLE FOR EACH OMP THREAD
!$     OMP_BLOCK = CEILING(REAL(ROW_LENGTH)/OMP_GET_NUM_THREADS())

  !          FILTER N->S PASS
!$OMP DO SCHEDULE(STATIC)
  DO jj=1, row_length, omp_block
    DO        jrow=2,n_rows
      DO        jpt=jj,MIN(jj+omp_block-1,row_length)
        f(jpt,jrow)=f(jpt,jrow)+                                       &
            (f(jpt,jrow-1)-f(jpt,jrow))*a_ns(jrow)
      END DO !)) JPT
    END DO !)) JROW
  END DO !)) JJ
!$OMP END DO

  !
  !  **   FILTER S->N             -----------------------------------
  !
  IF (m_grid == 0) THEN! BCS SPECIFIED TO ALLOW FOR CONSTANT OUTSIDE
    !      DO S BC
    IF (jpass == 1) THEN    !    ONE PREVIOUS N->S PASS
      !        BC FOR AFTER A SINGLE SWEEP IS USED FOR ALL LATER SWEEPS
!$OMP DO SCHEDULE(STATIC)
      DO      jpt=1,row_length
        f(jpt,n_rows)=f_boundary+(f(jpt,n_rows)-f_boundary)           &
                          /(1.0+a_ns(n_rows)    )
      END DO ! JPT
!$OMP END DO
    ELSE IF (jpass >= 2) THEN    !    TWO PREVIOUS S->N PASSES
      !                                       BC FOR 2 IS ALSO USED FOR >2
!$OMP DO SCHEDULE(STATIC)
      DO      jpt=1,row_length
        f(jpt,n_rows)=f_boundary+(f(jpt,n_rows)-f_boundary-           &
         a_ns(n_rows)    **3*(f(jpt,n_rows-1)-f_boundary))            &
         /((1.0+a_ns(n_rows)    )*(1.0-a_ns(n_rows)    **2))
      END DO ! JPT
!$OMP END DO
    END IF ! JPASS
    !      ELSE  ! M_GRID /= 0  NO SPECIAL TREATMENT OF BOUNDARY
  END IF ! M_GRID
  !
  !          FILTER N->S PASS

  omp_block = row_length
  !      ENSURE AS LARGE A STRIDE AS POSSIBLE FOR EACH OMP THREAD
!$     OMP_BLOCK = CEILING(REAL(ROW_LENGTH)/OMP_GET_NUM_THREADS())

!$OMP DO SCHEDULE(STATIC)
  DO        jj=1,row_length,omp_block
    DO        jrow=n_rows-1,1,-1
      DO        jpt=jj,MIN(jj+omp_block-1,row_length)
        f(jpt,jrow)=f(jpt,jrow)+                                       &
            (f(jpt,jrow+1)-f(jpt,jrow))*a_ns(jrow)
      END DO !)) JPT
    END DO !)) JROW
  END DO !)) JJ
!$OMP END DO


  !
END DO ! JPASS ====================================================

!$OMP END PARALLEL
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE rfcsl

END MODULE rfcsl_mod
