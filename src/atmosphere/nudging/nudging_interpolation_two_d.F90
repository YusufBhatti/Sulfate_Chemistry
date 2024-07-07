! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------
! Description:
!   Inetrpolates a 3D variable (ie 1D + 2D)

!  Part of the Nudged model (see nudging_main.F90)

!  Called from NUDGING_MAIN.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Nudging

! Code description:
!   Language: FORTRAN 90

!------------------------------------------------------------------
MODULE nudging_interpolation_two_d_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER,                   &
               PRIVATE :: ModuleName = 'NUDGING_INTERPOLATION_TWO_D_MOD'
CONTAINS

SUBROUTINE nudging_interpolation_two_d(      &
  input_dim1,                                & ! extra dimension 1
  input_dim2,                                & ! extra dimension 2
  input_limit1,                              & ! starting variable
  input_limit2,                              & ! finishing variable
  step,                                      & ! Fraction between start and end
  variable,                                  & ! Interpolated variable
  interp_type)                                 ! Interp type

USE nudging_control

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook
USE umPrintMgr

IMPLICIT NONE

! *****************************************************************************

INTEGER, INTENT(IN) :: input_dim1       ! extra dimension 1
INTEGER, INTENT(IN) :: input_dim2       ! extra dimension 2

! Variable at start at end of interpolation
REAL, INTENT(IN)    :: input_limit1(1:input_dim1,1:input_dim2)
REAL, INTENT(IN)    :: input_limit2(1:input_dim1,1:input_dim2)
! Interpolated variable
REAL, INTENT(OUT)   :: variable (1:input_dim1, 1:input_dim2)

REAL, INTENT(IN)    :: step             ! Fraction between 1 & 2
INTEGER, INTENT(IN) :: interp_type      ! Type of Interpolation

! Slope for linear interpolation
REAL                :: delta_variable(1:input_dim1,1:input_dim2)
INTEGER             :: i, j, k          ! Looping variables

INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='NUDGING_INTERPOLATION_TWO_D'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ****************************************************************************
! End of Header

! Linear Interpolation  (At present only have linear interpolation)

! Loop over the arrays performing the interpolation
DO i=1,input_dim1
  DO j=1,input_dim2
    delta_variable(i,j) = input_limit2(i,j) - input_limit1(i,j)
    variable(i,j) = input_limit1(i,j) + delta_variable(i,j)*step
  END DO
END DO ! 3D Linear Interpolation

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE nudging_interpolation_two_d

END MODULE nudging_interpolation_two_d_mod
