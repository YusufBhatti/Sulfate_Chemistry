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

MODULE atm_step_timestep_mod

USE get_wallclock_time_mod, ONLY: get_wallclock_time
 
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ATM_STEP_TIMESTEP_MOD'

CONTAINS
SUBROUTINE atm_step_time_init

USE timestep_mod, ONLY: timestep, timestep_number, recip_timestep,      &
                        pos_timestep, neg_timestep
USE submodel_mod, ONLY: atmos_im
USE umPrintMgr, ONLY: umPrint, umMessage, PrintStatus, PrStatus_Normal, &
                      PrStatus_Min

USE um_parcore, ONLY: mype
USE dynamics_testing_mod, ONLY:  l_backwards
USE turb_diff_mod, ONLY: l_print_pe
USE run_info, ONLY: start_time, time_start_atmstep

USE io_configuration_mod, ONLY: print_runtime_info
USE IOS_print_mgr, ONLY: IOS_print_start_time
USE model_time_mod, ONLY:                                                &
    i_day, i_hour, i_minute, i_month, i_second, i_year, secs_per_stepim, &
    stepim

USE parkind1, ONLY: jpim, jprb
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

! -------------------------------------------------------------------
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


CHARACTER(LEN=*), PARAMETER :: RoutineName='ATM_STEP_TIME_INIT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

timestep        = secs_per_stepim(atmos_im)  ! timestep in seconds
timestep_number = stepim(atmos_im)           ! no. of steps since basis

recip_timestep  = 1.0 / timestep                ! recip model timestep

pos_timestep =  timestep
neg_timestep = -timestep

IF (L_Backwards) timestep = neg_timestep

! Extra whitespace and line of asterisks to better delineate timesteps
IF (.NOT. (IOS_print_start_time .OR. print_runtime_info)) THEN
  CALL umPrint( '',src='atm_step_time_init')
  CALL umPrint( '****************************************'                  &
               //'****************************************',                &
               src='atm_step_time_init')
  CALL umPrint( '',src='atm_step_time_init')
END IF

IF (print_runtime_info) THEN
  IF ( L_print_pe .OR. mype == 0 ) THEN
    time_start_atmstep = get_wallclock_time()
    WRITE(umMessage,'(A,I8,A,F10.3,A)')                         &
      'Atm_Step: Info: Starting timestep ', timestep_number,    &
      ' at time=',time_start_atmstep - Start_time,' seconds'
    CALL umPrint(umMessage,src='atm_step_time_init')
    CALL umPrint( '',src='atm_step_time_init')
  END IF ! L_print_pe .or. mype == 0
END IF  ! print_runtime_info

IF (PrintStatus  >=  PrStatus_Normal) THEN

  IF ( L_Backwards ) THEN
    CALL umPrint( '',src='atm_step_time_init')
    CALL umPrint( '          *************************', &
        src='atm_step_time_init')
    CALL umPrint( 'Atm_Step: * INTEGRATING BACKWARDS *', &
        src='atm_step_time_init')
    CALL umPrint( '          *************************', &
        src='atm_step_time_init')
    CALL umPrint( '',src='atm_step_time_init')
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

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE atm_step_time_init

SUBROUTINE atm_step_timestep

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr, ONLY: umPrint, umMessage, PrintStatus, PrStatus_Min
USE um_parcore, ONLY: mype
USE timestep_mod, ONLY: timestep_number
USE submodel_mod, ONLY: atmos_im
USE turb_diff_mod, ONLY: l_print_pe
USE run_info, ONLY: start_time, time_start_atmstep
USE io_configuration_mod, ONLY: print_runtime_info
USE IOS_print_mgr, ONLY: IOS_print_start_time
USE model_time_mod, ONLY:                                                &
    i_day, i_hour, i_minute, i_month, i_second, i_year, stepim

IMPLICIT NONE

! ----------------------------------------------------------------------


! timestep information

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ATM_STEP_TIMESTEP'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

timestep_number  = stepim(atmos_im)        ! no. of steps since basis

! Extra whitespace and line of asterisks to better delineate timesteps
IF (.NOT. (IOS_print_start_time .OR. print_runtime_info)) THEN
  CALL umPrint( '',src='atm_step_timestep')
  CALL umPrint( '****************************************'              &
               //'****************************************',            &
               src='atm_step_timestep')
  CALL umPrint( '',src='atm_step_timestep')
END IF

IF (print_runtime_info) THEN
  IF ( L_print_pe .OR. mype == 0 ) THEN
    time_start_atmstep = get_wallclock_time()
    WRITE(umMessage,'(A,I8,A,F10.3,A)')                         &
      'Atm_Step: Info: Starting timestep ', timestep_number,    &
      ' at time=',time_start_atmstep - Start_time,' seconds'
    CALL umPrint(umMessage,src='atm_step_timestep')
    CALL umPrint( '',src='atm_step_timestep')
  END IF ! L_print_pe .or. mype == 0
END IF  ! print_runtime_info

IF (PrintStatus  >=  PrStatus_Min) THEN

  IF ( L_print_pe .OR. mype == 0 ) THEN
    WRITE(umMessage,'(A,I8,A,I6,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)')       &
          'Atm_Step: Timestep ', timestep_number, '   Model time: ', i_year,&
          '-',i_month,'-', i_day, ' ', i_hour, ':', i_minute, ':',i_second
    CALL umPrint(umMessage,src='atm_step_timestep')
  END IF ! L_print_pe .or. mype == 0

END IF     ! PrintStatus

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE atm_step_timestep

SUBROUTINE atm_step_info

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr, ONLY: umPrint, umMessage
USE um_parcore, ONLY: mype
USE timestep_mod, ONLY: timestep_number
USE submodel_mod, ONLY: atmos_im
USE turb_diff_mod, ONLY: l_print_pe
USE run_info, ONLY: time_start_atmstep

IMPLICIT NONE

! ----------------------------------------------------------------------

REAL :: time_end_atmstep

! timestep information

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ATM_STEP_INFO'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

time_end_atmstep = get_wallclock_time()

IF ( L_print_pe .OR. mype ==0 ) THEN
  CALL umPrint( '',src='atm_step_info')
  CALL umPrint( '****************************************'            &
         //'****************************************',                &
         src='atm_step_info')
  WRITE(umMessage,'(A,I8,A,F10.3,A)')                                 &
    'Atm_Step: Info: timestep ', timestep_number,                     &
    ' took ', time_end_atmstep - time_start_atmstep, ' seconds'
  CALL umPrint(umMessage,src='atm_step_info')
  CALL umPrint( '****************************************'            &
         //'****************************************',                &
         src='atm_step_info')
END IF ! L_print_pe .or. mype ==0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE atm_step_info

END MODULE atm_step_timestep_mod
