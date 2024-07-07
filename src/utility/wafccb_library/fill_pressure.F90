! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE fill_pressure_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE fill_pressure( Input_fld,    &
                          Input_mask,   &
                          Output_fld,   &
                          rmdi)      

! Description:
!   Interpolate field values within Cb mask area.
!
! Method:
! Performs linear interpolation horizontally E-W and vertically N-S
! to fill in the field values (pressure or heights of Cb bases and tops
! or horizontal extent of Cb) at points within the Cb mask that have
! missing values.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: WAFC CB Library
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6
!

IMPLICIT NONE


REAL, INTENT(INOUT) :: Input_fld(:, :)  ! Input/output field index
REAL, INTENT(IN)    :: Input_mask(:, :) ! Mask field index
REAL, INTENT(INOUT) :: Output_fld(:, :) ! Dummy field index (can
                                        ! also be used for output)
REAL, INTENT(IN)    :: rmdi

INTEGER :: i, j, k                                  ! loop counters
REAL    :: a                                        ! holds adjusted field value

!  Variables p & q contain the input field values in a 3x3 box around the point
!  of interest (p5/q5)
!
!  p1,p2,p3,
!  p4,p5,p6,
!  p7,p8,p9,
!
!  q1,q2,q3,
!  q4,q5,q6,
!  q7,q8,q9
!
! Only p5, q3, q5, q6, q9 are used in this routine
!
REAL    :: p5, q3, q5, q6, q9

INTEGER :: rows
INTEGER :: cols

! End of header --------------------------------------------------------

cols   = SIZE(input_mask, 1)
rows   = SIZE(input_mask, 2)


! Copy input field to dummy field ready for checking the field and making
! adjustments.
Output_fld(:,:) = Input_fld(:,:)


!---- Check through the field and adjust according to neighbouring points ----
!
! If input field has a missing value at point of interest but mask field
! indicates Cb, check if there is an input field value at the southerly,
! northery or easterly neighbouring gridpoints. If so, set value at point
! of interest to be this value. Repeat check through field a further 9 times.
!
! Optimisation comments:
! kk Exchanging I and J loop avoids bank conflicts AND allows vectorization.

DO k = 1, 10
  DO j = 1, rows
    DO i = 1, cols

      p5 = Input_mask(i,j)

      q3 = rmdi
      IF ((j /= rows) .AND. &
          (i /= cols)) THEN
        q3 = Output_fld(i+1,j+1)
      END IF

      q5 = Output_fld(i,j)

      q6 = rmdi
      IF (i /= cols) THEN
        q6 = Output_fld(i+1,j)
      END IF

      q9 = rmdi
      IF ((j /= 1) .AND. (i /= cols)) THEN
        q9 = Output_fld(i+1,j-1)
      END IF

      IF ( (q5 == rmdi ) .AND. ( p5 /= rmdi)) THEN
        a=rmdi
        IF (q9 /= rmdi ) THEN  
          a=q9
        END IF
        IF (q3 /= rmdi ) THEN
          a=q3
        END IF
        IF (q6 /= rmdi ) THEN
          a=q6
        END IF
        Output_fld(i,j)=a
      END IF
    END DO
  END DO
END DO


! Copy adjusted field back to Pfields(Inf1) ready for output
Input_fld(:,:) = Output_fld(:,:)

END SUBROUTINE fill_pressure

END MODULE fill_pressure_mod
