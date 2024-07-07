! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Stochastic Physics (sect35) Random Parameters Ver. 2
MODULE stph_rp_pert_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='STPH_RP_PERT_MOD'

CONTAINS


SUBROUTINE stph_rp_pert(l_firstcall, rp_rand0, rp_rand, rp_0, rp_max,   &
                        rp_min, rp_pert)

USE stochastic_physics_run_mod, ONLY:                                   &
    ran_max, ran_count, rp2_callfreq, rp2_decorr_ts,                    & 
    i_rp_scheme, i_rp2, i_rp2b

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

! Use switch for single amplitude of shock term in AR1 process
USE science_fixes_mod,    ONLY: l_fix_rp_shock_amp

IMPLICIT NONE
!
! Description:
!   Performs calculation of new perturbed value of parameter passed
!   from stph_rp2
!
! Method:
!   Increment perturbation using AR1 process
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Stochastic Physics
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

!
! IN variables
!
REAL, INTENT(IN)    :: rp_max            ! Max value of random parameter
REAL, INTENT(IN)    :: rp_min            ! Min value of random parameter
REAL, INTENT(IN)    :: rp_rand0(ran_max) ! Random number for firstcall
REAL, INTENT(IN)    :: rp_rand(ran_max)  ! Random number
LOGICAL, INTENT(IN) :: l_firstcall       ! 1st call to RP?
!
! IN/OUT variables
!
REAL, INTENT(INOUT) :: rp_pert           ! Perturbed parameter value
REAL, INTENT(INOUT) :: rp_0              ! Default (deterministic value)
!
! Local variables
!
REAL                :: shock             ! shock term in perturbation
REAL                :: maxran, minran    ! max and min values of paramter
!
! Local variables RP2b only
!
REAL                :: rand_mult         ! multiplier based on random number
REAL                :: phi_pert          ! perturbed parameter value in
                                         ! phi space
!
! Variables associated with autoregression model
!
REAL                :: corr              ! Correlation value for the
                                         ! first-order autoregression

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='STPH_RP_PERT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Perturb parameters using AR1 process.
! Different methods are used depending on which RP
! set-up is called (i_rp2 or i_rp2b).
!----------------------------------------------------------------------
IF (i_rp_scheme == i_rp2) THEN
  !--------------------------------------------------------------------
  ! On first call set gui/namelist default field and randomise initial
  ! pert value within range (uniform distribution)
  !--------------------------------------------------------------------
  IF (l_firstcall) THEN
    rp_0 = rp_pert
    rp_pert = rp_rand0(ran_count) * (rp_max - rp_min) + rp_min
  END IF

  !--------------------------------------------------------------------
  ! Perturb parameter using AR1 process
  !--------------------------------------------------------------------
  maxran = (rp_max - rp_min)/3.0
  minran = SIGN(maxran, -1.0)
  shock  = rp_rand(ran_count) * (maxran - minran) + minran

  corr = 1.0 - (rp2_callfreq / rp2_decorr_ts)
  rp_pert = rp_0 + ((rp_pert - rp_0)*corr) + shock
  rp_pert = MAX( rp_pert, rp_min )
  rp_pert = MIN( rp_pert, rp_max )

ELSE ! i_rp_scheme == i_rp2b
  !--------------------------------------------------------------------
  !  On first call set the default field and randomise initial pert 
  !  value and correlation value
  !
  !  For all other calls, transfer rp_pert to phi space to perform 
  !  update rp_min --> phi_pert = -1, rp_0 --> phi_pert = 0, 
  !  rp_max --> phi_pert = 1
  !
  !  For the case rp_min = rp_0 = rp_max, phi_pert is set to zero
  !--------------------------------------------------------------------
  IF (l_firstcall) THEN
    rand_mult = 2.0*rp_rand0(ran_count) - 1.0
    rp_0 = rp_pert
    phi_pert = rand_mult
  ELSE
    IF (rp_pert < rp_0) THEN
       phi_pert = (rp_pert - rp_0)/(rp_0 - rp_min)
    ELSE IF (rp_pert > rp_0) THEN
       phi_pert = (rp_pert - rp_0)/(rp_max - rp_0)
    ELSE
       phi_pert = 0.0
    END IF
  END IF

  !--------------------------------------------------------------------
  !  Perturb parameter in phi space using AR1 process
  !--------------------------------------------------------------------
  corr = 1.0 - (rp2_callfreq / rp2_decorr_ts)
  rand_mult = 2.0*rp_rand(ran_count) - 1.0
  IF (l_fix_rp_shock_amp) THEN
    shock = rand_mult * SQRT(1.0 - corr*corr)
  ELSE
    shock = 2.0 * rand_mult * SQRT(1.0 - corr*corr)
  END IF
  phi_pert = phi_pert * corr + shock

  !--------------------------------------------------------------------
  !  Reflect parameter value at boundaries
  !--------------------------------------------------------------------
  IF (phi_pert > 1.0) phi_pert = 2.0 - phi_pert
  IF (phi_pert < -1.0) phi_pert = - 2.0 - phi_pert

  !--------------------------------------------------------------------
  !  Transfer back to parameter space
  !--------------------------------------------------------------------
  IF (phi_pert < 0.0) THEN
    rp_pert = rp_0 + (rp_0 - rp_min) * phi_pert
  ELSE
    rp_pert = rp_0 + (rp_max - rp_0) * phi_pert
  END IF

END IF

! Increment random array counter (limit by maximum size)
ran_count = ran_count + 1
ran_count = MIN(ran_count, ran_max)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE stph_rp_pert
END MODULE stph_rp_pert_mod
