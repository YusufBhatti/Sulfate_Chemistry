! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!     SUBROUTINE RFCSLR
!     RECURSIVE FILTER - CONSTANT COEFF - FOR A SCALAR LTD.AREA FIELD
!
!     THIS VERSION IS FOR A REGULAR GRID AND DLAT AND DLONG ARE
!     SCALARS.
!
!     DOCUMENTATION:
!     THIS IS EXACTLY LIKE RFVSL (SEE BELOW),
!     EXCEPT THAT SCALE IS A SINGLE VALUE, APPLICABLE EVERYWHERE.
!
!     ARGUMENTS:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation

MODULE rfcslr_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RFCSLR_MOD'

CONTAINS

SUBROUTINE rfcslr(f,n_rows,row_length,m_grid,f_boundary,          &
 coslat,dlat,dlong,SCALE,npass                                    &
  )
!FPP$ NOCONCUR R

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE
!---- INPUT
INTEGER :: n_rows      ! NUMBER OF ROWS IN GRID
INTEGER :: row_length  ! NUMBER OF POINTS IN EACH ROW
INTEGER :: m_grid      ! 0=USE F_BOUNDARY, -1=NO SPECIAL BOUNDARY
REAL :: f_boundary     ! F ASSUMED OUTSIDE LTD.AREA
REAL :: coslat(n_rows) ! COSINE(LATITUDE)
REAL :: dlat                     ! ROW SPACING (RADIANS)
REAL :: dlong                    ! POINT SPACING (RADIANS)
REAL :: SCALE                    ! FILTER SCALE (RADIANS)
INTEGER :: npass       ! NUMBER OF PASSES IN EACH DIRECTION
!             NPASS=2 WILL GIVE EQUIVALENT OF SOAR FUNCTION (USUAL)
!             NPASS=LARGE WILL GIVE EQUIVALENT OF GAUSSIAN FUNCTION
!---- IN-OUT
REAL :: f(row_length,n_rows)     ! FIELD TO BE FILTERED
!
!---- WORK SPACE
REAL :: a_ns                            ! FILTER COEFFICIENTS N-S-N
REAL :: a_we(           n_rows)         ! W-E-W FILTER COEFFICIENTS
REAL :: e     ! INTERMEDIATE VALUE IN FORMULA FOR FILTER COEFFICIENT
REAL :: z     ! WORK VARIABLE FOR CALCULATING COEFFS
INTEGER :: jpass,jrow,jpt     ! LOOP COUNTERS FOR PASS,ROW,PT IN ROW

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RFCSLR'
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
DO      jrow=1,n_rows
  z=npass*(dlong*coslat(jrow))**2*0.25
  e=z/SCALE          **2
  a_we(    jrow)= 1.0+e-SQRT(e*(e+2.0))
END DO ! JROW
!
!---- CALCULATE COEFFS FOR N->S SWEEP (ALSO USED FOR S->N)
!       USING 6.5.13 & 6.5.9
z=npass*dlat**2*0.25
e=z/SCALE          **2
a_ns          = 1.0+e-SQRT(e*(e+2.0))

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
      DO      jrow=1,n_rows
        f(1,jrow)=f_boundary+(f(1,jrow)-f_boundary)                   &
        *(1.0-a_we(  jrow))
      END DO ! JROW
    ELSE IF (jpass == 2) THEN    !    ONE PREVIOUS E->W PASS
      DO      jrow=1,n_rows
        f(1,jrow)=f_boundary+(f(1,jrow)-f_boundary)                   &
                      /(1.0+a_we(  jrow))
      END DO ! JROW
    ELSE IF (jpass >= 3) THEN    !    TWO PREVIOUS E->W PASSES
      !                                       BC FOR 2 IS ALSO USED FOR >2
      DO      jrow=1,n_rows
        f(1,jrow)=f_boundary+(f(1,jrow)-f_boundary-                   &
             a_we(  jrow)**3*(f(2,jrow)-f_boundary))                  &
         /((1.0+a_we(  jrow))*(1.0-a_we(  jrow)**2))
      END DO ! JROW
    END IF ! JPASS
    !      ELSE  ! M_GRID /= 0  NO SPECIAL TREATMENT OF BOUNDARY
  END IF ! M_GRID
  !
  !          FILTER W->E PASS      EQN 6.5.1
  DO        jpt=2,row_length
    DO        jrow=1,n_rows
      f(jpt,jrow)=f(jpt,jrow)+                                       &
          (f(jpt-1,jrow)-f(jpt,jrow))*a_we(    jrow)
    END DO !)) JROW
  END DO !)) JPT
  !
  !  **  FILTER E->W              -----------------------------------
  !
  IF (m_grid == 0) THEN! BCS SPECIFIED TO ALLOW FOR CONSTANT OUTSIDE
    !      DO E BC
    IF (jpass == 1) THEN    !    ONE PREVIOUS W->E PASS
      DO      jrow=1,n_rows
        f(row_length,jrow)=f_boundary+(f(row_length,jrow)-f_boundary) &
                               /(1.0+a_we(           jrow))
      END DO ! JROW
    ELSE IF (jpass >= 2) THEN    !    TWO PREVIOUS W->E PASSES
      !                                       BC FOR 2 IS ALSO USED FOR >2
      DO      jrow=1,n_rows
        f(row_length,jrow)=f_boundary+(f(row_length,jrow)-f_boundary- &
         a_we(           jrow)**3*(f(row_length-1,jrow)-f_boundary))  &
         /((1.0+a_we(           jrow))*(1.0-a_we(           jrow)**2))
      END DO ! JROW
    END IF ! JPASS
    !      ELSE  ! M_GRID /= 0  NO SPECIAL TREATMENT OF BOUNDARY
  END IF ! M_GRID
  !
  !          FILTER E->W PASS
  DO        jpt=row_length-1,1,-1
    DO        jrow=1,n_rows
      f(jpt,jrow)=f(jpt,jrow)+                                       &
          (f(jpt+1,jrow)-f(jpt,jrow))*a_we(    jrow)
    END DO !)) JROW
  END DO !)) JPT
  !
  !  **   FILTER N->S             -----------------------------------
  !
  IF (m_grid == 0) THEN! BCS SPECIFIED TO ALLOW FOR CONSTANT OUTSIDE
    !      DO N BC
    IF (jpass == 1) THEN    !     NO PREVIOUS S->N PASS
      DO      jpt=1,row_length
        f(jpt,1)=f_boundary+(f(jpt,1)-f_boundary)                     &
         *(1.0-a_ns       )
      END DO ! JPT
    ELSE IF (jpass == 2) THEN    !    ONE PREVIOUS S->N PASS
      DO      jpt=1,row_length
        f(jpt,1)=f_boundary+(f(jpt,1)-f_boundary)                     &
                     /(1.0+a_ns       )
      END DO ! JPT
    ELSE IF (jpass >= 3) THEN    !    TWO PREVIOUS S->N PASSES
      !                                       BC FOR 2 IS ALSO USED FOR >2
      DO      jpt=1,row_length
        f(jpt,1)=f_boundary+(f(jpt,1)-f_boundary-                     &
         a_ns       **3*(f(jpt,2)-f_boundary))                        &
         /((1.0+a_ns       )*(1.0-a_ns       **2))
      END DO ! JPT
    END IF ! JPASS
    !      ELSE  ! M_GRID /= 0  NO SPECIAL TREATMENT OF BOUNDARY
  END IF ! M_GRID
  !
  !          FILTER N->S PASS
  DO        jrow=2,n_rows
    DO        jpt=1,row_length
      f(jpt,jrow)=f(jpt,jrow)+                                       &
          (f(jpt,jrow-1)-f(jpt,jrow))*a_ns
    END DO !)) JPT
  END DO !)) JROW
  !
  !  **   FILTER S->N             -----------------------------------
  !
  IF (m_grid == 0) THEN! BCS SPECIFIED TO ALLOW FOR CONSTANT OUTSIDE
    !      DO S BC
    IF (jpass == 1) THEN    !    ONE PREVIOUS N->S PASS
      !        BC FOR AFTER A SINGLE SWEEP IS USED FOR ALL LATER SWEEPS
      DO      jpt=1,row_length
        f(jpt,n_rows)=f_boundary+(f(jpt,n_rows)-f_boundary)           &
                          /(1.0+a_ns            )
      END DO ! JPT
    ELSE IF (jpass >= 2) THEN    !    TWO PREVIOUS S->N PASSES
      !                                       BC FOR 2 IS ALSO USED FOR >2
      DO      jpt=1,row_length
        f(jpt,n_rows)=f_boundary+(f(jpt,n_rows)-f_boundary-           &
         a_ns            **3*(f(jpt,n_rows-1)-f_boundary))            &
         /((1.0+a_ns            )*(1.0-a_ns            **2))
      END DO ! JPT
    END IF ! JPASS
    !      ELSE  ! M_GRID /= 0  NO SPECIAL TREATMENT OF BOUNDARY
  END IF ! M_GRID
  !
  !          FILTER N->S PASS
  DO        jrow=n_rows-1,1,-1
    DO        jpt=1,row_length
      f(jpt,jrow)=f(jpt,jrow)+                                       &
          (f(jpt,jrow+1)-f(jpt,jrow))*a_ns
    END DO !)) JPT
  END DO !)) JROW
  !
END DO ! JPASS ====================================================
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE rfcslr

END MODULE rfcslr_mod
