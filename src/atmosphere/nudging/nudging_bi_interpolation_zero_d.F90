! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------
! Description:
! One dimensional bi-linear interpolation (i.e 0 extra dimensions)

!  Part of the Nudged model (see nudging_main.F90)

!  Called from various

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Nudging

! Code description:
!   Language: FORTRAN 90

!------------------------------------------------------------------
MODULE nudging_bi_interpolation_zero_d_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER,                  &
              PRIVATE :: ModuleName = 'NUDGING_BI_INTERPOLATION_ZERO_D_MOD'
CONTAINS

SUBROUTINE nudging_bi_interpolation_zero_d(  &
 input_limit11,                              &  ! Variable at x_0, y_0
 input_limit12,                              &  ! Variable at x_0, y_1
 input_limit21,                              &  ! Variable at x_1, y_0
 input_limit22,                              &  ! Variable at x_1, y_1
 step1,                                      &  ! fraction between x_0 and x_1
 step2,                                      &  ! fraction between y_0 and y_1
 variable,                                   &  ! intrepolated variable
 interp_type)                                   ! Type of interpolation

USE nudging_control

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook
USE umPrintMgr


IMPLICIT NONE

! ***************************************************************************

REAL, INTENT(IN)    :: input_limit11
REAL, INTENT(IN)    :: input_limit12
REAL, INTENT(IN)    :: input_limit21
REAL, INTENT(IN)    :: input_limit22
REAL, INTENT(OUT)   :: variable
REAL, INTENT(IN)    :: step1
REAL, INTENT(IN)    :: step2
INTEGER, INTENT(IN) :: interp_type


! ****************************************
! Variables for Linear Interpolation
REAL :: variable_11
REAL :: variable_21
REAL :: variable_12
REAL :: variable_22

INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='NUDGING_BI_INTERPOLATION_ZERO_D'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! **************************************************************************
! End of Header

! ***********************************
! Linear Interpolation  (At present only have linear interpolation)

! Performing the interpolation
variable_11 =  input_limit11 *(1-step1)*(1-step2)
variable_12 =  input_limit12 *(1-step1)*(step2)
variable_21 =  input_limit21 *(step1)*(1-step2)
variable_22 =  input_limit22 *(step1)*(step2)

variable       = variable_11 + variable_12                            &
               + variable_21 + variable_22

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE nudging_bi_interpolation_zero_d

END MODULE nudging_bi_interpolation_zero_d_mod
