! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!      The U.M. version of this deck was made from the PP module PPRFILT
!      To keep the two codes equivalent,  changes should be made to
!      PPRFILT, and the U.M. deck re-created from it.
!      Instructions for doing this are in PPRFILT.
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation

!     SUBROUTINE RFVVL
!     RECURSIVE FILTER - VARIABLE COEFF - FOR A VECTOR LTD.AREA FIELD

MODULE rfvvl_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RFVVL_MOD'

CONTAINS

SUBROUTINE rfvvl(u,v,n_rows,row_length,m_grid,                    &
 coslat,firstrow,dlat,dlong,SCALE,npass                           &
 )

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE
!---- INPUT
INTEGER :: n_rows      ! NUMBER OF ROWS IN GRID
INTEGER :: row_length  ! NUMBER OF POINTS IN EACH ROW
INTEGER :: m_grid      ! 0=ASSUME 0 OUTSIDE. -1=NO SPECIAL BOUNDARY
REAL :: coslat(n_rows) ! COSINE(LATITUDE)
!      EXCEPT FOR ROWS AT POLE, WHICH HAVE COSLAT(1)=COSLAT(2)/8
!      & COSLAT(NNS)=COSLAT(NNS-1)/8 IF M_GRID=1.
!      TO ALLOW FOR THE NON-ZERO AREA OF GRID TRIANGLE SEGMENTS AT POLE
REAL :: firstrow                 ! CO-LATITUDE OF FIRST ROW (RADIANS)
REAL :: dlat                     ! ROW SPACING (RADIANS)
REAL :: dlong                    ! POINT SPACING (RADIANS)
REAL :: SCALE(row_length,n_rows) ! FILTER SCALE (RADIANS)
INTEGER :: npass       ! NUMBER OF PASSES IN EACH DIRECTION
!             NPASS=2 WILL GIVE EQUIVALENT OF SOAR FUNCTION (USUAL)
!             NPASS=LARGE WILL GIVE EQUIVALENT OF GAUSSIAN FUNCTION
!---- IN-OUT
REAL :: u(row_length,n_rows)     ! U FIELD TO BE FILTERED
REAL :: v(row_length,n_rows)     ! V FIELD TO BE FILTERED
!
!---- WORK SPACE
!
REAL :: a_ns(row_length,n_rows)         ! FILTER COEFFICIENTS N-S-N
REAL :: a_we(row_length,n_rows)         ! W-E-W FILTER COEFFICIENTS
REAL :: zu(n_rows)   ! RECURSIVELY CALCULATED U VALUE IN W-E-W SWEEPS
REAL :: zv(n_rows)   ! RECURSIVELY CALCULATED V VALUE IN W-E-W SWEEPS
REAL :: c_rot(n_rows)  ! COS(ROTATION) : ROTATION=SIN(LAT)*DLONG
REAL :: s_rot(n_rows)  ! SIN(ROTATION) : ROTATION=SIN(LAT)*DLONG

REAL :: zzu,zzv     ! TEMP U,V IN W-E-W FILTER
REAL :: e     ! INTERMEDIATE VALUE IN FORMULA FOR FILTER COEFFICIENT
REAL :: z     ! WORK VARIABLE FOR CALCULATING COEFFS
INTEGER :: jpass,jrow,jpt     ! LOOP COUNTERS FOR PASS,ROW,PT IN ROW

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RFVVL'

! ----------------------------------------------------------------------
!           PROGRAMMING NOTES
!     THE ORDER OF NESTING OF DOUBLE LOOPS ON JPT & JROW IS OPTIONAL.
!     SUCH LOOPS ARE INDICATED BY !)) ON THE ENDDO STATEMENT.
!     THE LOOP IN THE DIRECTION OF FILTERING IS RECURSIVE, AND SO CANNOT
!     BE VECTORIZED; HERE IT IS MADE THE OUTER LOOP.
!     A SCALAR COMPUTER MAY BE MORE EFFICIENT WITH A RECURSIVE INNER
!-----------------------------------------------------------------------
! **        BRIEF DESCRIPTION OF METHOD
!     THE INPUT FIELD IS VECTOR F=(U,V)
!     FIRST PASS IN ONE DIRECTION GIVES VECTOR G (OVERWRITING ARRAY F):
!      G(I)=A(I)*G(I-1)+B(I)*F(I)  :  I=1,N
!     PASS IN OPPOSITE DIRECTION GIVES VECTOR H (OVERWRITING ARRAY F):
!      H(I)=A(I)*H(I+1)+B(I)*G(I)  :  I=N,1
!     FOR FILTER ALONG ROW (CALLED E-W-E IN COMMENTS) THE GRID IS
!      REGULAR, SO B=1-A, AND A IS PRECALCULATED.
!     FOR FILTER ALONG COLUMNS (CALLED N-S-N IN COMMENTS) THE GRID IS
!      ASSUMED REGULAR, SO B=1-A, AND A IS PRECALCULATED.
!      BECAUSE OF THIS APPROX, THE GRID SHOULD NOT GO NEAR POLE.
!-----------------------------------------------------------------------
!     NO EXTERNAL SUBROUTINE CALLS
!-----------------------------------------------------------------------
!  ** PRECALCULATE COEFFICIENTS
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
z=firstrow       ! N.B. THIS IS CO-LATITUDE
DO jrow=1,n_rows ! CALCULATE ROTATION COEFFICIENTS
  e=COS(z)*dlong
  c_rot(jrow)=COS(e)
  s_rot(jrow)=SIN(e)
  z=z+dlat
END DO ! JROW
!
!---- CALCULATE COEFFS FOR E->W SWEEP (ALSO USED FOR W->E)
DO      jrow=1,n_rows
  z=npass*(dlong*coslat(jrow))**2*0.25
  DO      jpt=1,row_length
    !      STRICTLY, SCALE SHOULD BE AVERAGE JPT & JPT+(-)1 FOR W->E (E->W)
    !       NOT DOING SO MAKES THE CODE SIMPLER.
    !       SINCE SCALE IS SMOOTH, THE DIFFERENCE IS NEGLIGABLE.
    e=z/SCALE(jpt,jrow)**2
    a_we(jpt,jrow)= 1.0+e-SQRT(e*(e+2.0))
  END DO ! JPT
END DO ! JROW
!
!---- CALCULATE COEFFS FOR N->S SWEEP (ALSO USED FOR S->N)
z=npass*dlat**2*0.25
DO      jrow=1,n_rows
  DO      jpt=1,row_length
    !      STRICTLY, SCALE SHOULD BE AVERAGE JROW & JROW+(-)1
    !       NOT DOING SO MAKES THE CODE SIMPLER.
    !       SINCE SCALE IS SMOOTH, THE DIFFERENCE IS NEGLIGABLE.
    e=z/SCALE(jpt,jrow)**2
    a_ns(jpt,jrow)= 1.0+e-SQRT(e*(e+2.0))
  END DO ! JPT
END DO ! JROW
!-----------------------------------------------------------------------
!
!  ** START LOOP OVER PASSES
DO jpass=1,npass     ! ===========================================
  !
  !  **  FILTER W->E              -----------------------------------
  !
  !
  IF (m_grid == 0) THEN! BCS SPECIFIED TO ALLOW FOR CONSTANT OUTSIDE
    !      DO W BC
    !      THE BOUNDARY CONDITIONS DEPEND ON THE NUMBER OF PREVIOUS PASSES
    !      IN THE OPPOSITE DIRECTION.
    !      THE FORMULAE FOR 0,1,2 PREVIOUS PASSES ARE:-
    !        G1=(1-A)F1
    !        G1=(1/(1+A))F1
    !        G1=(1/(1+A)(1-A**2))(F1-A**3*F2)
    !      SEE METO11 WORKING PAPER 91 BY R J PURSER (1987) FOR DETAILS
    IF (jpass == 1) THEN    !     NO PREVIOUS E->W PASS
      DO      jrow=1,n_rows
        u(1,jrow)=u(1,jrow)*(1.0-a_we(1,jrow))
        v(1,jrow)=v(1,jrow)*(1.0-a_we(1,jrow))
        zu(jrow)=u(1,jrow)
        zv(jrow)=v(1,jrow)
      END DO ! JROW
    ELSE IF (jpass == 2) THEN    !    ONE PREVIOUS E->W PASS
      DO      jrow=1,n_rows
        u(1,jrow)=u(1,jrow)/(1.0+a_we(1,jrow))
        v(1,jrow)=v(1,jrow)/(1.0+a_we(1,jrow))
        zu(jrow)=u(1,jrow)
        zv(jrow)=v(1,jrow)
      END DO ! JROW
    ELSE IF (jpass >= 3) THEN    !    TWO PREVIOUS E->W PASSES
      !                                       BC FOR 2 IS ALSO USED FOR >2
      DO      jrow=1,n_rows     !    EQN 6.5.28
        u(1,jrow)=(u(1,jrow)-a_we(1,jrow)**3*u(2,jrow))               &
         /((1.0+a_we(1,jrow))*(1.0-a_we(1,jrow)**2))
        v(1,jrow)=(v(1,jrow)-a_we(1,jrow)**3*v(2,jrow))               &
         /((1.0+a_we(1,jrow))*(1.0-a_we(1,jrow)**2))
        zu(jrow)=u(1,jrow)
        zv(jrow)=v(1,jrow)
      END DO ! JROW
    END IF ! JPASS
    !      ELSE  ! M_GRID /= 0  NO SPECIAL TREATMENT OF BOUNDARY
  END IF ! M_GRID
  !
  !          FILTER W->E PASS
  DO        jpt=2,row_length
    DO        jrow=1,n_rows
      u(jpt,jrow)=u(jpt,jrow)+a_we(jpt,jrow)*                        &
       (c_rot(jrow)*zu(jrow)+s_rot(jrow)*zv(jrow)-u(jpt,jrow))
      v(jpt,jrow)=v(jpt,jrow)+a_we(jpt,jrow)*                        &
       (c_rot(jrow)*zv(jrow)-s_rot(jrow)*zu(jrow)-v(jpt,jrow))
      zu(jrow)=u(jpt,jrow)
      zv(jrow)=v(jpt,jrow)
    END DO !)) JROW
  END DO !)) JPT
  !
  !  **  FILTER E->W              -----------------------------------
  !
  IF (m_grid == 0) THEN! BCS SPECIFIED TO ALLOW FOR CONSTANT OUTSIDE
    !      DO E BC
    IF (jpass == 1) THEN    !    ONE PREVIOUS W->E PASS
      DO      jrow=1,n_rows
        u(row_length,jrow)=u(row_length,jrow)                         &
                               /(1.0+a_we(row_length,jrow))
        v(row_length,jrow)=v(row_length,jrow)                         &
                               /(1.0+a_we(row_length,jrow))
        zu(jrow)=u(row_length,jrow)
        zv(jrow)=v(row_length,jrow)
      END DO ! JROW
    ELSE IF (jpass >= 2) THEN    !    TWO PREVIOUS W->E PASSES
      !                                       BC FOR 2 IS ALSO USED FOR >2
      DO      jrow=1,n_rows     !    EQN 6.5.28
        u(row_length,jrow)=(u(row_length,jrow)-                       &
         a_we(row_length,jrow)**3*u(row_length-1,jrow))               &
         /((1.0+a_we(row_length,jrow))*(1.0-a_we(row_length,jrow)**2))
        v(row_length,jrow)=(v(row_length,jrow)-                       &
         a_we(row_length,jrow)**3*v(row_length-1,jrow))               &
         /((1.0+a_we(row_length,jrow))*(1.0-a_we(row_length,jrow)**2))
        zu(jrow)=u(row_length,jrow)
        zv(jrow)=v(row_length,jrow)
      END DO ! JROW
    END IF ! JPASS
    !      ELSE  ! M_GRID /= 0  NO SPECIAL TREATMENT OF BOUNDARY
  END IF ! M_GRID
  !
  !          FILTER E->W PASS
  DO        jpt=row_length-1,1,-1
    DO        jrow=1,n_rows
      u(jpt,jrow)=u(jpt,jrow)+a_we(jpt,jrow)*                        &
       (c_rot(jrow)*zu(jrow)-s_rot(jrow)*zv(jrow)-u(jpt,jrow))
      v(jpt,jrow)=v(jpt,jrow)+a_we(jpt,jrow)*                        &
       (c_rot(jrow)*zv(jrow)+s_rot(jrow)*zu(jrow)-v(jpt,jrow))
      zu(jrow)=u(jpt,jrow)
      zv(jrow)=v(jpt,jrow)
    END DO !)) JROW
  END DO !)) JPT
  !
  !  **   FILTER N->S             -----------------------------------
  !
  IF (m_grid == 0) THEN! BCS SPECIFIED TO ALLOW FOR CONSTANT OUTSIDE
    !      DO N BC
    IF (jpass == 1) THEN    !     NO PREVIOUS S->N PASS
      DO      jpt=1,row_length
        u(jpt,1)=u(jpt,1)*(1.0-a_ns(jpt,1))
        v(jpt,1)=v(jpt,1)*(1.0-a_ns(jpt,1))
      END DO ! JPT
    ELSE IF (jpass == 2) THEN    !    ONE PREVIOUS S->N PASS
      DO      jpt=1,row_length
        u(jpt,1)=u(jpt,1)/(1.0+a_ns(jpt,1))
        v(jpt,1)=v(jpt,1)/(1.0+a_ns(jpt,1))
      END DO ! JPT
    ELSE IF (jpass >= 3) THEN    !    TWO PREVIOUS S->N PASSES
      !                                       BC FOR 2 IS ALSO USED FOR >2
      DO      jpt=1,row_length  !    EQN 6.5.28
        u(jpt,1)=(u(jpt,1)-a_ns(jpt,1)**3*u(jpt,2))                   &
         /((1.0+a_ns(jpt,1))*(1.0-a_ns(jpt,1)**2))
        v(jpt,1)=(v(jpt,1)-a_ns(jpt,1)**3*v(jpt,2))                   &
         /((1.0+a_ns(jpt,1))*(1.0-a_ns(jpt,1)**2))
      END DO ! JPT
    END IF ! JPASS
    !      ELSE  ! M_GRID /= 0  NO SPECIAL TREATMENT OF BOUNDARY
  END IF ! M_GRID
  !
  !          FILTER N->S PASS
  DO        jrow=2,n_rows
    DO        jpt=1,row_length
      u(jpt,jrow)=u(jpt,jrow)+                                       &
                 (u(jpt,jrow-1)-u(jpt,jrow))*a_ns(jpt,jrow)
      v(jpt,jrow)=v(jpt,jrow)+                                       &
                 (v(jpt,jrow-1)-v(jpt,jrow))*a_ns(jpt,jrow)
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
        u(jpt,n_rows)=u(jpt,n_rows)/(1.0+a_ns(jpt,n_rows))
        v(jpt,n_rows)=v(jpt,n_rows)/(1.0+a_ns(jpt,n_rows))
      END DO ! JPT
    ELSE IF (jpass >= 2) THEN    !    TWO PREVIOUS S->N PASSES
      !                                       BC FOR 2 IS ALSO USED FOR >2
      DO      jpt=1,row_length  !    EQN 6.5.28
        u(jpt,n_rows)=                                                &
          (u(jpt,n_rows)-a_ns(jpt,n_rows)**3*u(jpt,n_rows-1))         &
         /((1.0+a_ns(jpt,n_rows))*(1.0-a_ns(jpt,n_rows)**2))
        v(jpt,n_rows)=                                                &
          (v(jpt,n_rows)-a_ns(jpt,n_rows)**3*v(jpt,n_rows-1))         &
         /((1.0+a_ns(jpt,n_rows))*(1.0-a_ns(jpt,n_rows)**2))
      END DO ! JPT
    END IF ! JPASS
    !      ELSE  ! M_GRID /= 0  NO SPECIAL TREATMENT OF BOUNDARY
  END IF ! M_GRID
  !
  !          FILTER N->S PASS
  DO        jrow=n_rows-1,1,-1
    DO        jpt=1,row_length
      u(jpt,jrow)=u(jpt,jrow)+                                       &
          (u(jpt,jrow+1)-u(jpt,jrow))*a_ns(jpt,jrow)
      v(jpt,jrow)=v(jpt,jrow)+                                       &
          (v(jpt,jrow+1)-v(jpt,jrow))*a_ns(jpt,jrow)
    END DO !)) JPT
  END DO !)) JROW
  !
END DO ! JPASS ====================================================
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE rfvvl
END MODULE rfvvl_mod
