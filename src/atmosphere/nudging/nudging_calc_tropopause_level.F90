! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------
! Description:
!  Converts tropopause pressure to model level
!
!  Part of the Nudged model (see nudging_main.F90)
!
!  Called from NUDGING_NUDGING_CONTROL
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Nudging
!
! Code description:
!   Language: FORTRAN 90 (formatted)
!
!------------------------------------------------------------------
MODULE nudging_calc_tropopause_level_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER,                &
              PRIVATE :: ModuleName = 'NUDGING_CALC_TROPOPAUSE_LEVEL_MOD'
CONTAINS

SUBROUTINE nudging_calc_tropopause_level(  &
  proc_row_length_min,                     & ! Min column
  proc_row_Length_max,                     & ! Max column
  proc_rows_min,                           & ! Min row
  proc_rows_max,                           & ! Max row
  model_levels,                            & ! No. levels
  model_pressure_model_levels,             & ! Pressure on model levels
  model_trop_pressure,                     & ! Tropopause pressure
  tropopause_level)                          ! Output tropopause levl

USE nudging_control

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook
USE umPrintMgr

IMPLICIT NONE

! ******************************************************************

INTEGER, INTENT(IN)    :: proc_row_length_min
INTEGER, INTENT(IN)    :: proc_row_length_max
INTEGER, INTENT(IN)    :: proc_rows_min
INTEGER, INTENT(IN)    :: proc_rows_max
INTEGER, INTENT(IN)    :: model_levels
REAL, INTENT(IN)       :: model_pressure_model_levels( &
     proc_row_length_min:proc_row_length_max,          &
     proc_rows_min:proc_rows_max, 1:model_levels)
REAL, INTENT(IN)       :: model_trop_pressure(         &
       proc_row_length_min:proc_row_length_max,        &
       proc_rows_min:proc_rows_max)
INTEGER, INTENT(OUT)   :: tropopause_level(            &
       proc_row_length_min:proc_row_length_max,        &
       proc_rows_min:proc_rows_max)

INTEGER                :: i, j, k

INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='NUDGING_CALC_TROPOPAUSE_LEVEL'

!-----------------------------------
!EOH

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO i=proc_row_length_min, proc_row_length_max
  DO j=proc_rows_min, proc_rows_max

    k=1
    DO WHILE(model_trop_pressure(i,j) <                     &
         model_pressure_model_levels(i,j,k)                 &
         .AND. k < model_levels)
      k=k+1
    END DO

    tropopause_level(i,j) = k

  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nudging_calc_tropopause_level

END MODULE nudging_calc_tropopause_level_mod
