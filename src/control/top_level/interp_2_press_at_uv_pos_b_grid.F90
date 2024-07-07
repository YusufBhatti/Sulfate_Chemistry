! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Interpolate vertically from model levels onto pressure levels and
! horizontally from atmospheric u points on the c grid to uv positions
! on the B grid.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

SUBROUTINE Interp_2_Press_At_UV_Pos_B_Grid(                       &
           im_index, item, section,                               &
           from, row_length, rows, n_rows, num_p_levels,          &
           model_levels, off_x,off_y, exner_theta_levels, interpolated)

USE stash_array_mod, ONLY: stindex, stlist, stash_levels
USE planet_constants_mod, ONLY: kappa, p_zero

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE interpor_mod, ONLY: interp_order_linear
USE uc_to_ub_mod, ONLY: uc_to_ub
USE vert_interp2_mod, ONLY: vert_interp2

IMPLICIT NONE

INTEGER :: i,j,k
INTEGER :: im_index
INTEGER :: item, section
INTEGER :: row_length, rows, n_rows
INTEGER :: num_p_levels, model_levels
INTEGER :: off_x, off_y

REAL :: exner_theta_levels                                      &
   (1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels)

REAL :: from(row_length,rows,model_levels)
REAL :: exner_at_theta(row_length,rows,model_levels)
REAL :: interpolated(row_length, n_rows, num_p_levels)
REAL :: work(row_length, rows, num_p_levels)

INTEGER :: idx

REAL :: exner

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INTERP_2_PRESS_AT_UV_POS_B_GRID'


IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

!  Calculate exner at theta points. Store in exner_at_theta

DO k = 1, model_levels
  DO j = 1, rows
    DO i = 1, row_length
      exner_at_theta(i,j,k) = exner_theta_levels(i,j,k)
    END DO ! i
  END DO ! j
END DO ! k

idx = -stlist(10,stindex(1,item,section,im_index))
DO k = 1, num_p_levels
  exner = (((stash_levels(k+1,idx)/1000.0) * 100.0)/p_zero)**kappa

  CALL vert_interp2 (from                                       &
    ,row_length, rows, model_levels                             &
    ,exner                                                      &
    ,0,0,0,0                                                    &
    ,exner_at_theta, interp_order_linear                        &
    ,work(:,:,k))
END DO

! Perform simple horizontal interpolation from 'C' to 'B' grid
CALL  uC_to_uB(work, row_length, rows, n_rows, num_p_levels,    &
              off_x, off_y, interpolated)

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE Interp_2_Press_At_UV_Pos_B_Grid
