! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: SSTpert library

MODULE sstpert_mod

IMPLICIT NONE 

CONTAINS

SUBROUTINE sstpert(factor, dt, num_rows, num_cols, fieldclim, fieldout) &
                   BIND(c, NAME="sstpert")

! Description:  Generates a SST perturbation field based on dimensions
!               of input field
!
! Method:       Shares code with Section35 (Stochastic Physics) to
!               create Gaussian Field with prescribed power spectrum.
!               Modulates perturbation field by interpolating input
!               monthly climatology to analysis date.

USE, INTRINSIC :: iso_c_binding, ONLY: C_INT64_T, C_DOUBLE, C_F_POINTER, C_LOC
USE sst_genpatt_mod, ONLY: sst_genpatt
USE stph_seed_mod, ONLY: stph_seed_gen_from_date

IMPLICIT NONE

REAL(KIND=C_DOUBLE), INTENT(IN) :: factor
                       ! SST perturbation inflation factor (alpha)

INTEGER(KIND=C_INT64_T), INTENT(IN) :: dt(8)
                       ! Array of date information (used in first stage of
                       ! random seed generation).  Should contain:
                       !  * (1) Year of input T* field
                       !  * (2) Month of input T* field
                       !  * (3) Day of input T* field
                       !  * (4) UTC Shift of input T* field
                       !  * (5) Hour of input T* field
                       !  * (6) Minute of input T* field
                       !  * (7) Ensemble member number
                       !  * (8) Ensemble member numer + 100

INTEGER(KIND=C_INT64_T), INTENT(IN) :: num_rows
                       ! Number of rows in field data

INTEGER(KIND=C_INT64_T), INTENT(IN) :: num_cols
                       ! Number of columns in field data

REAL(KIND=C_DOUBLE), INTENT(IN) :: fieldclim(12, num_cols, num_rows)
                       ! Monthly DeltaSST climatology arrays 

REAL(KIND=C_DOUBLE), INTENT(OUT), TARGET :: fieldout(num_cols, num_rows)
                       ! Perturbed output field

REAL(KIND=C_DOUBLE), POINTER :: field2d(:,:)

INTEGER :: i_month1, i_month2
                       ! months adjacent to current date
INTEGER :: i,j
REAL :: r_month1, r_month2
                       ! months adjacent to current date
REAL, ALLOCATABLE ::                                                    &
        rnum(:)
REAL, ALLOCATABLE ::                                                    &
        psif(:,:)

CALL C_F_POINTER (C_LOC(fieldout), field2d, [num_cols, num_rows])

! Setup the random seed using the provided target date
CALL stph_seed_gen_from_date(dt)

! Find the two adjacent months nearest AINITIAL date
!   and calculate ratio of each month
IF ( dt(3) < 15 ) THEN
  i_month1 = dt(2) - 1
  r_month1 = (16 - dt(3))/30.0
  r_month2 = 1.0 - r_month1
ELSE
  i_month1 = dt(2)
  r_month2 = (dt(3) - 15)/30.0
  r_month1 = 1.0 - r_month2
END IF
i_month2 = i_month1 + 1
IF ( i_month2 > 12 ) i_month2 = i_month2 - 12
IF ( i_month1 < 1 ) i_month1 = i_month1 + 12

! Generate Random field with specified power spectrum
ALLOCATE(psif(num_cols, num_rows))

CALL sst_genpatt(num_cols, num_rows, psif)

! Find absolute sqrt of forcing pattern (power is in units K^2)
DO j = 1, num_rows
  DO i = 1, num_cols
    psif(i,j) = SIGN( SQRT( ABS( psif(i,j) )), psif(i,j))
  END DO
END DO

!  multiply by interpolated DSST field and add to surface temp
DO j = 1, num_rows
  DO i = 1, num_cols
    field2d(i,j) = psif(i,j) * factor *                               &
                   (r_month1 * fieldclim(i_month1,i,j) +               &
                    r_month2 * fieldclim(i_month2,i,j))
  END DO
END DO

NULLIFY(field2d)
DEALLOCATE(psif)

END SUBROUTINE sstpert

END MODULE sstpert_mod
