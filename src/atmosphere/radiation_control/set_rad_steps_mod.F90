! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Description:
!   Module contains subroutines to set which timesteps to call the SW and LW
!   radiation schemes on, in prog or diag mode.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
! Method:
!   Subroutine set_a_radstep is called in readlsta (called from um_shell to
!     read the UM namelists), to set the number of radiation timesteps per 
!     model timestep.  This is needed in stash_proc, which is called from 
!     um_shell.  set_a_radstep also calculates the radiation timestep length
!     used in atm_step.
!   Subroutine set_l_rad_step is called in settsctl, to set the flags for
!     prog / diag radiation timesteps, used in u_model / u_model_4A.
!     They are needed here, inside the timestep loop but outside of
!     atm_step, because UKCA uses them too.
!
!   The calling structure of the control routines that calculate the 
!   radiation timestep variables, and those that use them, is shown below:
!
!   UM: um_shell
!           CALL readlsta
!                   CALL read_nml_run_radiation
!                   CALL check_run_radiation
!                   CALL set_a_radstep
!           CALL stash_proc
!           CALL u_model(_4A)
!                   [start time-step loop]
!                   CALL settsctl
!                           CALL set_l_rad_step
!                   CALL atm_step(_4A)
!                           CALL atmos_physics1
!                           CALL atmos_physics2
!                   CALL ukca_main1
!                   [end time-step loop]
!
!
!   In the SCM, set_a_radstep must be called lower down the calling tree,
!   because it requires the model time-step length, which isn't read until
!   the call to read_scm_nml in scm_main.
!   The calling structure is as follows:
!
!   SCM: scm_shell
!           READ(10,RUN_Radiation)
!           CALL check_run_radiation
!           CALL scm_main
!                   CALL read_scm_nml
!                   CALL set_a_radstep
!                   [start time-step loop]
!                   CALL pre_physics
!                           CALL set_l_rad_step
!                   CALL atmos_physics1
!                   CALL atmos_physics2
!                   [end time-step loop]


! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3

MODULE set_rad_steps_mod

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

! Flags for whether each step is a full or diagnostic radiation timestep
LOGICAL, SAVE :: l_rad_step_prog = .FALSE.
LOGICAL, SAVE :: l_rad_step_diag = .FALSE.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SET_RAD_STEPS_MOD'

CONTAINS



!----------------------------------------------------------------------------
! Subroutine to calculate the number of model timesteps per radiation 
! timestep, for prog and diag radiation calls.  These  are valid for
! the whole run.  Also calculates the radiation timestep length.
!----------------------------------------------------------------------------

SUBROUTINE set_a_radstep()

! Various settings and parameters that are needed here:
USE model_domain_mod, ONLY: model_type, mt_single_column
USE timestep_mod,     ONLY: timestep
USE nlstgen_mod,      ONLY: secs_per_periodim, steps_per_periodim
USE submodel_mod,     ONLY: atmos_im
USE conversions_mod,  ONLY: isec_per_day
USE rad_input_mod,    ONLY: i_rad_extra_call, ip_increment_call,         &
                            i_sw_radstep_perday_prog,                    &
                            i_lw_radstep_perday_prog,                    &
                            i_sw_radstep_perday_diag,                    &
                            i_lw_radstep_perday_diag

! a_radstep variables that are set by this routine:
USE rad_input_mod,    ONLY: a_sw_radstep_prog, a_lw_radstep_prog,        &
                            a_sw_radstep_diag, a_lw_radstep_diag

! Radiation time-step lengths, calculated here:
USE timestep_mod,     ONLY: radiation_tstep_prog, radiation_tstep_diag

IMPLICIT NONE

INTEGER :: secs_per_timestep
INTEGER :: timesteps_per_day

REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SET_A_RADSTEP'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


! Find the timestep length:
IF (model_type == mt_single_column) THEN
  ! In the SCM, this subroutine is called after the timestep has been set,
  ! so we can just use the value of timestep
  secs_per_timestep = INT(timestep)

ELSE
  ! In the full UM, the timestep hasn't been set yet, so we must calculate
  ! it from the nlstgen namelist inputs:
  secs_per_timestep = secs_per_periodim(atmos_im) / steps_per_periodim(atmos_im)

END IF

! Calculate number of timesteps per day
timesteps_per_day = isec_per_day / secs_per_timestep


! Number of timesteps per main SW and LW radiation step:
a_sw_radstep_prog = timesteps_per_day / i_sw_radstep_perday_prog
a_lw_radstep_prog = timesteps_per_day / i_lw_radstep_perday_prog

! Number of timesteps per "diag" SW and LW radiation step:
IF (i_rad_extra_call == ip_increment_call) THEN
  ! If doing extra radiation steps for the incremental 
  ! timestepping scheme, calculate number of timesteps per extra rad step:
  a_sw_radstep_diag = timesteps_per_day / i_sw_radstep_perday_diag
  a_lw_radstep_diag = timesteps_per_day / i_lw_radstep_perday_diag

ELSE
  ! If not doing extra radiation steps, set the number of timesteps per diag
  ! step the same as for the full radiation steps (for calculating diagnostics).
  a_sw_radstep_diag = a_sw_radstep_prog
  a_lw_radstep_diag = a_lw_radstep_prog

END IF


! Calculate the radiation timestep length (this is actually the SW radiation 
! timestep) for full radiation (_prog) and fast / diagnostic radiation (_diag)
! radiation calls.
radiation_tstep_prog = REAL(secs_per_timestep) * a_sw_radstep_prog
radiation_tstep_diag = REAL(secs_per_timestep) * a_sw_radstep_diag

! Note: radiation_tstep_diag is only used in rad_ctl if using the incremental
! radiative timestepping scheme (i_rad_extra_call==ip_increment_call).


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE set_a_radstep



!----------------------------------------------------------------------------
! Subroutine to be called each timestep, to set whether each type of 
! radiation call (prog or diag) should be done this timestep.
!----------------------------------------------------------------------------

SUBROUTINE set_l_rad_step(stepcount)

! Various settings and parameters that are needed here:
USE rad_input_mod,    ONLY: a_sw_radstep_prog, a_sw_radstep_diag,     &
                            a_lw_radstep_prog, a_lw_radstep_diag,     &
                            it_rad1

IMPLICIT NONE

INTEGER, INTENT(IN) :: stepcount  ! Current timestep number

REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SET_L_RAD_STEP'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


! Set whether this is a radiation timestep, using the a_radstep variables
! (number of model timesteps per radiation timestep) set earlier.
l_rad_step_prog = (MOD(stepcount-it_rad1,a_sw_radstep_prog) == 0)     &
             .OR. (MOD(stepcount-it_rad1,a_lw_radstep_prog) == 0)
l_rad_step_diag = (MOD(stepcount-it_rad1,a_sw_radstep_diag) == 0)     &
             .OR. (MOD(stepcount-it_rad1,a_lw_radstep_diag) == 0)

! Note: unless i_rad_extra_call==ip_increment_call, a_sw_radstep_diag and
! a_lw_radstep_diag should both be the same as their _prog counterparts, so
! that radiation diagnostics are calculated on (and only on) the same 
! timesteps as the full radiation calls.


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE set_l_rad_step



END MODULE set_rad_steps_mod
