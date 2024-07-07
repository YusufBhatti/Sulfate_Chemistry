! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE NEAR_PT-------------------------------------------------
!
!    Purpose:  To produce gather indices which map each point on the
!              target grid onto the nearest point on the source grid.
!              This allows interpolation by choosing the value of the
!              nearest neighbour. The code uses coefficients and gather
!              indices calculated by subroutine H_INT_CO.
!
!    Programming standard:
!             Unified Model Documentation Paper No 3
!             Version No 1 15/1/90
!
!    System component:
!
!    System task: S123
!
!    Documentation: The interpolation formulae are described in
!                   unified model on-line documentation paper S1.
!
!    -------------------------------------------------------------------
!    Arguments:---------------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration

SUBROUTINE near_pt                                                &
(index_b_l,index_b_r,weight_t_r,weight_b_r,weight_t_l,weight_b_l  &
,points,points_lambda_srce,index_nearest)
USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim
IMPLICIT NONE

INTEGER ::                                                        &
 points_lambda_srce                                               &
                      !IN Number of lambda points on source grid
,points                                                           &
                      !IN Total number of points on target grid
,index_b_l(points)                                                &
                      !IN  Index of bottom lefthand corner
                      !    of source gridbox
,index_b_r(points)                                                &
                      !IN  Index of bottom righthand corner
                      !    of source gridbox
,index_nearest(points)!OUT Index of nearest source point to
                      ! each target point

REAL ::                                                           &
 weight_t_r(points)                                               &
                    !IN  Weight applied to value at top right
                    !    hand corner of source gridbox
,weight_b_l(points)                                               &
                    !IN  Weight applied to value at bottom left
                    !    hand corner of source gridbox
,weight_b_r(points)                                               &
                    !IN  Weight applied to value at bottom right
                    !    hand corner of source gridbox
,weight_t_l(points) !IN  Weight applied to value at top left
                    !    hand corner of source gridbox

! Local arrays:---------------------------------------------------------
INTEGER ::                                                        &
 index_temp(points,4)   ! Index of 4 sourrounding source points
                        ! ordered by distance

REAL ::                                                           &
 max_weight(points,4)   ! Linear interpolation weights ordered by
                        ! distance

!   External subroutines called:----------------------------------------
! None
! ----------------------------------------------------------------------
! Local variables:------------------------------------------------------
REAL :: temp
INTEGER :: i,j,k,itemp
CHARACTER (LEN=*), PARAMETER  :: RoutineName='NEAR_PT'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
! ----------------------------------------------------------------------
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! 1.  Accumulate source weights and indices associated with
!     each coastal point on target grid.

DO i=1,points

  max_weight(i,1)=weight_b_l(i)
  max_weight(i,2)=weight_b_r(i)
  max_weight(i,3)=weight_t_l(i)
  max_weight(i,4)=weight_t_r(i)
  index_temp(i,1)=index_b_l(i)
  index_temp(i,2)=index_b_r(i)
  index_temp(i,3)=index_b_l(i)                                      &
                  +points_lambda_srce
  index_temp(i,4)=index_b_r(i)                                      &
                  +points_lambda_srce
END DO

! 2.  Sort gather indices of the 4 surrounding source
!     gridpoints according to distance from target gridpoint;
!     arranged so that nearest point comes first in list (ie K=1).

DO k=1,3
  DO j=k+1,4
    DO i=1,points
      IF (max_weight(i,k) <  max_weight(i,j)) THEN
        temp=max_weight(i,k)
        max_weight(i,k)=max_weight(i,j)
        max_weight(i,j)=temp
        itemp=index_temp(i,k)
        index_temp(i,k)=index_temp(i,j)
        index_temp(i,j)=itemp
      END IF
    END DO
  END DO
END DO

! 3. Assign index of nearest source point to output array

DO i=1,points
  index_nearest(i)=index_temp(i,1)
END DO

IF (lhook) CALL dr_hook(RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE near_pt
