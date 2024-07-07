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
!     SUBROUTINE RFVSG
!     RECURSIVE FILTER - VARIABLE COEFF - FOR A SCALAR GLOBAL FIELD
!     DOCUMENTATION:
!       U.M. DOCUMENTATION PAPER 30
!       "THE ANALYSIS CORRECTION DATA ASSIMILATION SCHEME" SECTION 6.5
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation
MODULE rfvsg_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RFVSG_MOD'

CONTAINS

SUBROUTINE rfvsg(f,n_rows,row_length,m_grid,                      &
 coslat,dlat,dlong,SCALE,npass                                    &
  )

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE
!---- INPUT
INTEGER :: n_rows      ! NUMBER OF ROWS IN GRID
INTEGER :: row_length  ! NUMBER OF POINTS IN EACH ROW
INTEGER :: m_grid      ! 1=(PTS AT POLE), 2=(STAG'D FROM POLE)
REAL :: coslat(n_rows) ! COSINE(LATITUDE)
!      EXCEPT FOR ROWS AT POLE, WHICH HAVE COSLAT(1)=COSLAT(2)/8
!      & COSLAT(NNS)=COSLAT(NNS-1)/8 IF M_GRID=1.
!      TO ALLOW FOR THE NON-ZERO AREA OF GRID TRIANGLE SEGMENTS AT POLE
REAL :: dlat                     ! ROW SPACING (RADIANS)
REAL :: dlong                    ! POINT SPACING (RADIANS)
REAL :: SCALE(row_length,n_rows) ! FILTER SCALE (RADIANS)
INTEGER :: npass       ! NUMBER OF PASSES IN EACH DIRECTION
!             NPASS=2 WILL GIVE EQUIVALENT OF SOAR FUNCTION (USUAL)
!             NPASS=LARGE WILL GIVE EQUIVALENT OF GAUSSIAN FUNCTION
!---- IN-OUT
REAL :: f(row_length,n_rows)     ! FIELD TO BE FILTERED
!
!---- WORK SPACE
REAL :: a(row_length,n_rows)         ! FILTER COEFFICIENTS N-S-N
REAL :: b(row_length,n_rows)         ! FILTER COEFFICIENTS N-S-N
REAL :: a_we(row_length,n_rows)         ! W-E-W FILTER COEFFICIENTS
REAL :: a_product(n_rows)     ! PRODUCT OF A_WE ROUND LATITUDE CIRCLE
REAL :: zf(n_rows)     ! RECURSIVELY CALCULATED VALUE IN W-E-W SWEEPS
REAL :: ones(row_length) ! RECURSIVE RESULT OF FILTERING ONES   N_S_N
!
REAL :: scale_n     ! FILTER SCALE AT N POLE
REAL :: scale_s     ! FILTER SCALE AT S POLE
REAL :: e     ! INTERMEDIATE VALUE IN FORMULA FOR FILTER COEFFICIENT
REAL :: z     ! WORK VARIABLE FOR CALCULATING COEFFS
REAL :: areg     ! A COEFF IF THE GRIDBOXES DID NOT CHANGE SIZE N-S
REAL :: f_pole     ! FILTERED VALUE OF FIELD AT POLE
REAL :: rw1,rw2     ! RATIOS OF WEIGHTS BETWEEN ADJACENT ROWS
REAL :: cotlat     ! COTAN(GRID LATITUDE)
REAL :: rscale     ! 1/SCALE FOR EXTRA SWEEP E-W-E
INTEGER :: irf,irg,irh,irl ! LIMITS TO ROWS DONE IN W-E-W SWEEPS
INTEGER :: jpass,jrow,jpt     ! LOOP COUNTERS FOR PASS,ROW,PT IN ROW

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RFVSG'

! ----------------------------------------------------------------------
!           PROGRAMMING NOTES
!     THE ORDER OF NESTING OF DOUBLE LOOPS ON JPT & JROW IS OPTIONAL.
!     SUCH LOOPS ARE INDICATED BY !)) ON THE ENDDO STATEMENT.
!     THE LOOP IN THE DIRECTION OF FILTERING IS RECURSIVE, AND SO CANNOT
!     BE VECTORIZED; HERE IT IS MADE THE OUTER LOOP, AND ZF IS A
!     TEMPORARY VECTOR, SO THE INNER LOOP CAN BE VECTORIZED.
!     A SCALAR COMPUTER MAY BE MORE EFFICIENT WITH A RECURSIVE INNER
!     LOOP, WITH ZF IN A TEMPORARY SCALAR.
!-----------------------------------------------------------------------
! **        BRIEF DESCRIPTION OF METHOD
!     THE INPUT FIELD IS F
!     FIRST PASS IN ONE DIRECTION GIVES G (OVERWRITING ARRAY F):
!      G(I)=A(I)*G(I-1)+B(I)*F(I)  :  I=1,N
!     PASS IN OPPOSITE DIRECTION GIVES H (OVERWRITING ARRAY F):
!      H(I)=C(I)*H(I+1)+D(I)*G(I)  :  I=N,1
!     FOR FILTER ALONG ROW (CALLED E-W-E IN COMMENTS) THE GRID IS
!      REGULAR, SO D=B=1-C=1-A, AND A IS PRECALCULATED.
!     FOR FILTER ALONG COLUMNS (CALLED N-S-N IN COMMENTS) THE GRID BOX
!      SIZE CHANGES, SO A B C D ARE CALCULATED TO ALLOW FOR THIS.
!      A B C D ARE ALPHA BETA GAMA DELTA IN DOCUMENTATION
!      C IS STORED IN A, D IS STORED IN B.
!-----------------------------------------------------------------------
!     NO EXTERNAL SUBROUTINE CALLS
!-----------------------------------------------------------------------
!  ** PRECALCULATE COEFFICIENTS
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
scale_n=0.0
scale_s=0.0
DO      jpt=1,row_length
  scale_n=scale_n+SCALE(jpt,1  )
  scale_s=scale_s+SCALE(jpt,n_rows)
END DO ! JPT
scale_n=scale_n/row_length
scale_s=scale_s/row_length
!
IF (m_grid == 1) THEN! ROWS 1&N_ROWS ARE POLES: NO W->E FILTER
  irf=2
  irl=n_rows-1
ELSE                ! FIRST & LAST ROWS ARE 1/2 G-L FROM POLES
  irf=1
  irl=n_rows
END IF ! M_GRID == 1
!
!---- CALCULATE COEFFS FOR E->W SWEEP (ALSO USED FOR W->E)
!               EQNS. 6.5.13 & 6.5.9
DO      jrow=irf,irl
  z=npass*(dlong*coslat(jrow))**2*0.25
  DO      jpt=1,row_length
    !      STRICTLY, SCALE SHOULD BE AVERAGE JPT & JPT+(-)1 FOR W->E (E->W)
    !       NOT DOING SO MAKES THE CODE SIMPLER.
    !       SINCE SCALE IS SMOOTH, THE DIFFERENCE IS NEGLIGABLE.
    e=z/SCALE(jpt,jrow)**2
    a_we(jpt,jrow)= 1.0+e-SQRT(e*(e+2.0))
  END DO ! JPT
  a_product(jrow)=a_we(1 ,jrow)
END DO ! JROW
!              FOR USE IN 6.5.23
DO        jpt=2,row_length
  DO        jrow=irf,irl
    a_product(jrow)=a_product(jrow)*a_we(jpt,jrow)
  END DO !)) JROW
END DO !)) JPT
!-----------------------------------------------------------------------
!
!  ** START LOOP OVER PASSES
DO jpass=1,npass     ! ===========================================
  !
  !  **  FILTER W->E              -----------------------------------
  !
  IF (m_grid == 1) THEN! REPLACE ROWS 1 & N_ROWS BY AVERAGE
    f_pole=0.0
    DO      jpt=1,row_length
      f_pole=f_pole+f(jpt,1  )
    END DO ! JPT
    f_pole=f_pole/row_length
    DO      jpt=1,row_length
      f(jpt,1  )=f_pole
    END DO ! JPT
    f_pole=0.0
    DO      jpt=1,row_length
      f_pole=f_pole+f(jpt,n_rows)
    END DO ! JPT
    f_pole=f_pole/row_length
    DO      jpt=1,row_length
      f(jpt,n_rows)=f_pole
    END DO ! JPT
  END IF ! M_GRID == 1
  !
  !      SET UP W BC FOR WRAP-ROUND BY FILTERING ACROSS GRID
  DO      jrow=irf,irl
    zf(jrow)=(1.0-a_we(1 ,jrow))*f(1 ,jrow)
  END DO ! JROW
  DO        jpt=2,row_length
    DO        jrow=irf,irl
      zf(jrow)=f(jpt,jrow)+(zf(jrow)-f(jpt,jrow))*a_we(jpt,jrow)
    END DO !)) JROW
  END DO !)) JPT
  !      SOLVE IMPLICIT EQ. 6.5.23 FOR WRAP-ROUND OF LATITUDE CIRCLE
  DO      jrow=irf,irl
    zf(jrow)=zf(jrow)/(1.0-a_product(jrow))
  END DO ! JROW
  !
  !          FILTER W->E PASS  (5.6.1A)
  DO        jpt=1,row_length
    DO        jrow=irf,irl
      zf(jrow)=f(jpt,jrow)+(zf(jrow)-f(jpt,jrow))*a_we(jpt,jrow)
      f(jpt,jrow)=zf(jrow)
    END DO !)) JROW
  END DO !)) JPT
  !
  !  **  FILTER E->W              -----------------------------------
  !
  !      SET UP E BC FOR WRAP-ROUND BY FILTERING ACROSS GRID
  DO      jrow=irf,irl
    zf(jrow)=(1.0-a_we(row_length,jrow))*f(row_length,jrow)
  END DO ! JROW
  DO        jpt=row_length-1,1,-1
    DO        jrow=irf,irl
      zf(jrow)=f(jpt,jrow)+(zf(jrow)-f(jpt,jrow))*a_we(jpt,jrow)
    END DO !)) JROW
  END DO !)) JPT
  !      SOLVE IMPLICIT EQ. 6.5.23 FOR WRAP-ROUND OF LATITUDE CIRCLE
  DO      jrow=irf,irl
    zf(jrow)=zf(jrow)/(1.0-a_product(jrow))
  END DO ! JROW
  !
  !          FILTER E->W PASS (5.6.1B)
  DO        jpt=row_length,1,-1
    DO        jrow=irf,irl
      zf(jrow)=f(jpt,jrow)+(zf(jrow)-f(jpt,jrow))*a_we(jpt,jrow)
      f(jpt,jrow)=zf(jrow)
    END DO !)) JROW
  END DO !)) JPT
  !
  !      N-S-N & S-N-S SWEEPS ARE ALTERNATED TO EQUALIZE FILTER AT POLES
  IF (MOD(jpass,2) == 1) THEN! SWEEP N->S->N
    !  **   CALCULATE COEFFICIENTS FOR A N->S SWEEP
    !         B FOR S.POLE ROW
    DO jpt=1,row_length
      b(jpt,n_rows)= 1.0
    END DO ! JPT
    !       CALCULATE COEFFS USING A RECURSIVE RELATIONSHIP TO CONSERVE
    !       THE (AREA WEIGHTED) OUTPUT FROM EACH GRIDPOINT.
    !       A FROM 6.5.22, B FROM 6.5.18
    DO      jrow=n_rows,2,-1
      rw1=coslat(jrow)/coslat(jrow-1)
      rw2=MAX(rw1,1.0/rw1)
      DO      jpt=1,row_length
        e=npass*dlat**2/(SCALE(jpt,jrow-1)+SCALE(jpt,jrow))**2
        areg= 1.0+e-SQRT(e*(e+2.0))
        a(jpt,jrow)=rw2*areg*(b(jpt,jrow)/(1.0-areg))
        b(jpt,jrow-1)=1.0/(1.0+rw1*a(jpt,jrow)/b(jpt,jrow))
      END DO ! JPT
    END DO ! JROW
    !
    !  **   FILTER N->S             -----------------------------------
    !
    !       AT THE SAME TIME CALCULATE COEFFICIENTS FOR A S->N SWEEP
    !       SPECIAL TREATMENT FOR ROW 1
    DO      jpt=1,row_length
      f(jpt,1)=f(jpt,1)*b(jpt,1)
      ones(jpt)=b(jpt,1)
      b(jpt,1)=1.0             ! D FOR NEXT (S->N) SWEEP 6.5.19
      a(jpt,1)=1.0-ones(jpt)   ! C FOR NEXT (S->N) SWEEP 6.5.20
    END DO ! JPT
    !       LOOP OVER ROWS 2,N_ROWS
    DO        jrow=2,n_rows
      DO        jpt=1,row_length
        f(jpt,jrow)=f(jpt,jrow)*b(jpt,jrow)+f(jpt,jrow-1)*a(jpt,jrow)
        ones(jpt)=            b(jpt,jrow)+ones(jpt)  *a(jpt,jrow)
        !          D(S->N) IS CALCULATED TO CONSERVE EACH INPUT GRID VALUE
        b(jpt,jrow)=1.0/ (1.0+(coslat(jrow-1)*a(jpt,jrow-1))/          &
                              (coslat(jrow)*b(jpt,jrow-1)))
        !          C(S->N) IS CALCULATED SUCH THAT A FIELD OF ONES IS UNCHANGED
        a(jpt,jrow)=1.0-ones(jpt)*b(jpt,jrow)
      END DO !)) JPT
    END DO !)) JROW
    !
    !  **   FILTER S->N             -----------------------------------
    !
    !       SPECIAL TREATMENT FOR ROW N_ROWS
    DO      jpt=1,row_length
      f(jpt,n_rows)=f(jpt,n_rows)/ones(jpt)
    END DO ! JPT
    DO        jrow=n_rows-1,1,-1
      DO        jpt=1,row_length
        f(jpt,jrow)=f(jpt,jrow)*b(jpt,jrow)+f(jpt,jrow+1)*a(jpt,jrow)
      END DO !)) JPT
    END DO !)) JROW
    !
    !      N-S-N & S-N-S SWEEPS ARE ALTERNATED TO EQUALIZE FILTER AT POLES
  ELSE  ! SWEEP S->N->S ON EVEN SWEEPS TO EQUALIZE FILTER OF POLES
    !
    !  **   CALCULATE COEFFICIENTS FOR A S->N SWEEP
    DO      jpt=1,row_length
      b(jpt,1  )= 1.0
    END DO ! JPT
    !       CALCULATE COEFFS USING A RECURSIVE RELATIONSHIP TO CONSERVE
    !       THE (AREA WEIGHTED) OUTPUT FROM EACH GRIDPOINT.
    !       A FROM 6.5.22, B FROM 6.5.18
    DO      jrow=1,n_rows-1
      rw1=coslat(jrow)/coslat(jrow+1)
      rw2=MAX(rw1,1.0/rw1)
      DO      jpt=1,row_length
        e=npass*dlat**2/(SCALE(jpt,jrow+1)+SCALE(jpt,jrow))**2
        areg= 1.0+e-SQRT(e*(e+2.0))
        a(jpt,jrow)=rw2*areg*(b(jpt,jrow)/(1.0-areg))
        b(jpt,jrow+1)=1.0/(1.0+rw1*a(jpt,jrow)/b(jpt,jrow))
      END DO ! JPT
    END DO ! JROW
    !
    !  **   FILTER S->N             -----------------------------------
    !
    !       AT THE SAME TIME CALCULATE COEFFICIENTS FOR A N->S SWEEP
    !       SPECIAL TREATMENT FOR ROW N_ROWS
    DO      jpt=1,row_length
      f(jpt,n_rows)=f(jpt,n_rows)*b(jpt,n_rows)
      ones(jpt)=b(jpt,n_rows)
      b(jpt,n_rows)=1.0
      a(jpt,n_rows)=1.0-ones(jpt)
    END DO ! JPT
    !       LOOP OVER ROWS N_ROWS-1,1
    DO        jrow=n_rows-1,1,-1
      DO        jpt=1,row_length
        f(jpt,jrow)=f(jpt,jrow)*b(jpt,jrow)+f(jpt,jrow+1)*a(jpt,jrow)
        ones(jpt)=            b(jpt,jrow)+ones(jpt)  *a(jpt,jrow)
        b(jpt,jrow)=1.0/ (1.0+(coslat(jrow+1)*a(jpt,jrow+1))/          &
                              (coslat(jrow)*b(jpt,jrow+1)))
        a(jpt,jrow)=1.0-ones(jpt)*b(jpt,jrow)
      END DO !)) JPT
    END DO !)) JROW
    !
    !  **   FILTER N->S             -----------------------------------
    !
    !       SPECIAL TREATMENT FOR ROW 1
    DO      jpt=1,row_length
      f(jpt,1  )=f(jpt,1  )/ones(jpt)
    END DO ! JPT
    DO        jrow=2,n_rows
      DO        jpt=1,row_length
        f(jpt,jrow)=f(jpt,jrow)*b(jpt,jrow)+f(jpt,jrow-1)*a(jpt,jrow)
      END DO !)) JPT
    END DO !)) JROW
    !      N-S-N & S-N-S SWEEPS ARE ALTERNATED TO EQUALIZE FILTER AT POLES
  END IF ! END OF (ODD) N-S-N OR (EVEN) S-N-S SWEEP
  !
END DO ! JPASS ====================================================
!
!  ** EXTRA W-E-W PASS NEAR POLES
!  **  FILTER W->E              -----------------------------------
IF (m_grid == 1) THEN!REPLACE  1&N_ROWS BY AVERAGE
  f_pole=0.0
  DO      jpt=1,row_length
    f_pole=f_pole+f(jpt,1  )
  END DO ! JPT
  f_pole=f_pole/row_length
  DO      jpt=1,row_length
    f(jpt,1  )=f_pole
  END DO ! JPT
  f_pole=0.0
  DO      jpt=1,row_length
    f_pole=f_pole+f(jpt,n_rows)
  END DO ! JPT
  f_pole=f_pole/row_length
  DO      jpt=1,row_length
    f(jpt,n_rows)=f_pole
  END DO ! JPT
END IF ! M_GRID == 1
irg=n_rows/2
irh=irg+1
DO jrow=irf,irl
  !       CALCULATE COEFF FOR THIS EXTRA SWEEP USING 6.5.24 & 6.5.25
  cotlat=coslat(jrow)/MAX(.000001,1.0-coslat(jrow)**2)
  IF (jrow <  n_rows/2) THEN! NEAR N.POLE
    rscale=MAX(cotlat/scale_n,1.0)/scale_n
  ELSE ! NEAR S.POLE
    rscale=MAX(cotlat/scale_s,1.0)/scale_s
  END IF ! NEAR N. OR S.POLE
  !       A_WE IS ASSUMED CONSTANT ALONE ROW FOR THIS EXTRA SWEEP
  !       SO ONLY THE FIRST POINT IN ROW IS USED.
  a_we(1,jrow)=MAX(0.0,1.0-2.000*dlong*coslat(jrow)*rscale)
  IF (a_we(1,jrow) >  0.0) THEN! NEAR POLE
    a_product(jrow)=a_we(1,jrow)**row_length
    !         SET UP W BC FOR WRAP-ROUND BY FILTERING ACROSS GRID
    zf(jrow)=(1.0-a_we(1,jrow))*f(1 ,jrow)
  ELSE ! DETERMINE LIMITS WHERE FILTER IS DOING NOTHING
    irg=MIN(irg,jrow-1)
    irh=MAX(irh,jrow+1)
  END IF ! NEAR POLE
END DO ! JROW
DO        jpt=2,row_length
  DO        jrow=irf,irg
    zf(jrow)=f(jpt,jrow)+(zf(jrow)-f(jpt,jrow))*a_we(1,jrow)
  END DO !)) JROW
END DO !)) JPT
DO        jpt=2,row_length
  DO        jrow=irh,irl
    zf(jrow)=f(jpt,jrow)+(zf(jrow)-f(jpt,jrow))*a_we(1,jrow)
  END DO !)) JROW
END DO !)) JPT
!      SOLVE IMPLICIT EQUATION FOR WRAP-ROUND FILTER OF LATITUDE CIRCLE
DO      jrow=irf,irg
  zf(jrow)=zf(jrow)/(1.0-a_product(jrow))
END DO ! JROW
DO      jrow=irh,irl
  zf(jrow)=zf(jrow)/(1.0-a_product(jrow))
END DO ! JROW
!          FILTER W->E PASS
DO        jpt=1,row_length
  DO        jrow=irf,irg
    zf(jrow)=f(jpt,jrow)+(zf(jrow)-f(jpt,jrow))*a_we(1,jrow)
    f(jpt,jrow)=zf(jrow)
  END DO !)) JROW
END DO !)) JPT
DO        jpt=1,row_length
  DO        jrow=irh,irl
    zf(jrow)=f(jpt,jrow)+(zf(jrow)-f(jpt,jrow))*a_we(1,jrow)
    f(jpt,jrow)=zf(jrow)
  END DO !)) JROW
END DO !)) JPT
!  **  FILTER E->W              -----------------------------------
!      SET UP E BC FOR WRAP-ROUND BY FILTERING ACROSS GRID
DO      jrow=irf,irg
  zf(jrow)=(1.0-a_we(1,jrow))*f(row_length,jrow)
END DO ! JROW
DO      jrow=irh,irl
  zf(jrow)=(1.0-a_we(1,jrow))*f(row_length,jrow)
END DO ! JROW
DO        jpt=row_length-1,1,-1
  DO        jrow=irf,irg
    zf(jrow)=f(jpt,jrow)+(zf(jrow)-f(jpt,jrow))*a_we(1,jrow)
  END DO !)) JROW
END DO !)) JPT
DO        jpt=row_length-1,1,-1
  DO        jrow=irh,irl
    zf(jrow)=f(jpt,jrow)+(zf(jrow)-f(jpt,jrow))*a_we(1,jrow)
  END DO !)) JROW
END DO !)) JPT
!      SOLVE IMPLICIT EQUATION FOR WRAP-ROUND FILTER OF LATITUDE CIRCLE
DO      jrow=irf,irg
  zf(jrow)=zf(jrow)/(1.0-a_product(jrow))
END DO ! JROW
DO      jrow=irh,irl
  zf(jrow)=zf(jrow)/(1.0-a_product(jrow))
END DO ! JROW
!          FILTER E->W PASS
DO        jpt=row_length,1,-1
  DO        jrow=irf,irg
    zf(jrow)=f(jpt,jrow)+(zf(jrow)-f(jpt,jrow))*a_we(1,jrow)
    f(jpt,jrow)=zf(jrow)
  END DO !)) JROW
END DO !)) JPT
DO        jpt=row_length,1,-1
  DO        jrow=irh,irl
    zf(jrow)=f(jpt,jrow)+(zf(jrow)-f(jpt,jrow))*a_we(1,jrow)
    f(jpt,jrow)=zf(jrow)
  END DO !)) JROW
END DO !)) JPT
! *** END OF EXTRA W-E-W PASS NEAR POLES
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE rfvsg

END MODULE rfvsg_mod
