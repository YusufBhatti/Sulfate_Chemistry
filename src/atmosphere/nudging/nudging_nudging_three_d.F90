! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------
! Description:
!  Nudges one variable towards another

!  Part of the Nudged model (see nudging_main.F90)

!  Called from NUDGING_NUDGING_CONTROL

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Nudging

! Code description:
!   Language: FORTRAN 90

!------------------------------------------------------------------
MODULE nudging_nudging_three_d_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER,             &
                 PRIVATE :: ModuleName = 'NUDGING_NUDGING_THREE_D_MOD'
CONTAINS

SUBROUTINE nudging_nudging_three_d(     &
  model_levels,                         &  ! GLOBAL No. of levels
  proc_row_length_min,                  &  ! Mimimum column
  proc_row_length_max,                  &  ! Maximum column
  proc_rows_min,                        &  ! Minimum row
  proc_rows_max,                        &  ! Maximum row
  variable_model,                       &  ! Model (both input and out)
  variable_data,                        &  ! input data
  relaxation_parameter)                    ! Relaxation parameter

USE nudging_control

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook
USE umPrintMgr

IMPLICIT NONE

INTEGER, INTENT(IN) :: model_levels         ! model levels
INTEGER, INTENT(IN) :: proc_row_length_min  ! min column
INTEGER, INTENT(IN) :: proc_row_length_max  ! max column
INTEGER, INTENT(IN) :: proc_rows_min        ! min row
INTEGER, INTENT(IN) :: proc_rows_max        ! max row

! model variable
REAL, INTENT(INOUT) :: variable_model(                                &
 proc_row_length_min:proc_row_length_max,                             &
 proc_rows_min:proc_rows_max,1:model_levels )

! data variable interpolated onto model levels
REAL, INTENT(IN)    :: variable_data (                                &
 proc_row_length_min:proc_row_length_max,                             &
 proc_rows_min:proc_rows_max,1:model_levels )

! Relaxation Parameter
REAL, INTENT(IN)    :: relaxation_parameter (                         &
 proc_row_length_min:proc_row_length_max,                             &
 proc_rows_min:proc_rows_max,1:model_levels)

REAL                :: delta_variable   ! Data Model Difference
INTEGER             :: i, j, k          ! Looping variables

INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='NUDGING_NUDGING_THREE_D'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ***************************************************************************
! End of Header

! initialise the diff
delta_variable = 0

! Loop round relevant rows, columns and levels performing the nudging
DO i = proc_rows_min, proc_rows_max
  DO j = proc_row_length_min , proc_row_length_max
    DO k = ndg_lev_bottom, ndg_lev_top

      delta_variable =                                               &
        variable_data(j,i,k) - variable_model(j,i,k)
      delta_variable =                                               &
        delta_variable*relaxation_parameter(j,i,k)
      variable_model(j,i,k) =                                        &
        variable_model(j,i,k) + delta_variable

    END DO  ! levels
  END DO ! columns
END DO ! 3D nudging (rows)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE nudging_nudging_three_d

END MODULE nudging_nudging_three_d_mod
