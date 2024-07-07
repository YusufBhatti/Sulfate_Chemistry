! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Updates time-correlated perturbation
!
MODULE update_pert_mod

IMPLICIT NONE
SAVE

! Description:
!  This routine updates a given perturbation using a simple
!  AR1 process.  This results in an updated perturbation
!  correlated in time with the original perturbation.
!
!  The perturbation is assumed to have zero mean with max and min
!  values specified at input.  Any perturbation calculated outside of 
!  the max or min is reflected back into the specified range.
!  The initial value of the perturbation must be set outside of this 
!  module.
!
!  The algorithm to update the perturbation includes a random shock
!  term with a specified amplitude.  For the case when a time-correlated
!  random number sequence is required, the shock amplitude should be set
!  to 1.
!
!  The auto-correlation coefficient, mu, should be calculated outside of
!  this module and typically takes the form EXP(-delta_t/delta_tau)
!  where delta_t is the timestep between perturbation updates, and 
!  delta_tau is the decorrelation timescale.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: stochastic_physics
!
! Code description:
!  Language: Fortran 95.
!  This code is written to UMDP3 standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UPDATE_PERT_MOD'

CONTAINS

! subroutine interface
SUBROUTINE update_pert(pert, mu, rand_no, shock_amp, pert_min, pert_max)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Arguments with INTENT IN. ie: Input variables.
REAL, INTENT(IN) :: mu             ! Auto-correlation coefficient in
                                   ! AR1 process
REAL, INTENT(IN) :: rand_no        ! Random number between 0 and 1
REAL, INTENT(IN) :: shock_amp      ! Noise amplitude in AR1 process
REAL, INTENT(IN) :: pert_min       ! Maximum value of perturbation
REAL, INTENT(IN) :: pert_max       ! Minimum value of perturbation

! Arguments with INTENT INOUT
REAL, INTENT(INOUT) :: pert        ! perturbation to update

! Local Variables
REAL :: rand_mult  ! Random number used in the AR1 process

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UPDATE_PERT'
!----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Adjust random number to vary between -1 and 1
rand_mult=2.0*rand_no - 1.0

! Evolve perturbation with an AR1 process
pert = mu*pert + SQRT(1.0-mu*mu)*rand_mult*shock_amp

! Reflect perturbation at specified maximum and minimum
IF (pert > pert_max) pert = 2.0*pert_max - pert
IF (pert < pert_min) pert = 2.0*pert_min - pert

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE update_pert

END MODULE update_pert_mod
