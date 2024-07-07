! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Time step initialisation section for Atm_Step
!
! Subroutine Interface:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

SUBROUTINE atm_step_timestep_init

USE timestep_mod
USE submodel_mod, ONLY: atmos_im
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE umPrintMgr
USE UM_ParVars
USE UM_ParCore, ONLY: mype
USE Control_Max_Sizes
USE dynamics_testing_mod, ONLY:  l_backwards
USE turb_diff_mod, ONLY: l_print_pe
USE run_info, ONLY: start_time, time_start_atmstep
USE io_configuration_mod, ONLY: print_runtime_info
USE IOS_print_mgr, ONLY: IOS_print_start_time

USE model_time_mod, ONLY:                                                &
    i_day, i_hour, i_minute, i_month, i_second, i_year, secs_per_stepim, &
    stepim

USE get_wallclock_time_mod, ONLY: get_wallclock_time
 
IMPLICIT NONE

! --------------------------------------------------------------------------


! timestep information

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ATM_STEP_TIMESTEP_INIT'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
timestep             = SECS_PER_STEPim(atmos_im)      ! timestep in seconds
timestep_number      = STEPim(atmos_im)               ! no. of steps since basis
recip_timestep       = 1.0 / timestep                 ! recip model timestep

pos_timestep =  timestep
neg_timestep = -timestep

IF (L_Backwards) timestep = neg_timestep

! Extra whitespace and line of asterisks to better delineate timesteps
IF (.NOT. (IOS_print_start_time .OR. print_runtime_info)) THEN
  CALL umPrint( '',src='atm_step_timestep_init')
  CALL umPrint( '****************************************'                  &
               //'****************************************',                &
               src='atm_step_timestep_init')
  CALL umPrint( '',src='atm_step_timestep_init')
END IF

IF (print_runtime_info) THEN
  IF ( L_print_pe .OR. mype == 0 ) THEN
    time_start_atmstep = get_wallclock_time()
    WRITE(umMessage,'(A,I8,A,F10.3,A)')                         &
      'Atm_Step: Info: Starting timestep ', timestep_number,    &
      ' at time=',time_start_atmstep - Start_time,' seconds'
    CALL umPrint(umMessage,src='atm_step_timestep_init')
    CALL umPrint( '',src='atm_step_timestep_init')
  END IF ! L_print_pe .or. mype == 0
END IF  ! print_runtime_info

IF (PrintStatus  >=  PrStatus_Normal) THEN

  IF (timestep_number  ==  1 .AND. L_Backwards) THEN
    CALL umPrint( '',src='atm_step_timestep_init')
    CALL umPrint( '          *************************', &
        src='atm_step_timestep_init')
    CALL umPrint( 'Atm_Step: * INTEGRATING BACKWARDS *', &
        src='atm_step_timestep_init')
    CALL umPrint( '          *************************', &
        src='atm_step_timestep_init')
    CALL umPrint( '',src='atm_step_timestep_init')
  END IF

END IF     ! PrintStatus

IF (PrintStatus  >=  PrStatus_Min) THEN

  IF ( L_print_pe .OR. mype == 0 ) THEN
    WRITE(umMessage,'(A,I8,A,I6,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)')         &
          'Atm_Step: Timestep ', timestep_number, '   Model time: ', i_year, &
          '-',i_month,'-', i_day, ' ', i_hour, ':', i_minute, ':',i_second
    CALL umPrint(umMessage,src='atm_step_timestep_init')
  END IF ! L_print_pe .or. mype == 0

END IF     ! PrintStatus

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE atm_step_timestep_init
