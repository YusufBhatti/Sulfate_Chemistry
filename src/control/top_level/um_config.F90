! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Purpose: Provide top level interface to model/application initialisation
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Top Level

MODULE UM_Config

! DEPENDS ON: exceptions

! for dependencies
USE io_dependencies
USE thread_utils
USE um_types
USE io, ONLY: ioInit, ioShutDown
USE application_description, ONLY:                                           &
    exe_UM, exe_RCF, exe_scm,                                                &
    exe_combine, exe_merge,                                                  &
    exe_hreset,                                                              &
    exe_convpp,                                                              &
    exe_setup, exe_fldmod, exe_hprint,                                       &
    exe_pptoanc, exe_pickup,                                                 &
    exe_vomext,exe_crmstyle_coarse_grid,                                     &
    setApplicationDesc,                                                      &
    getExeType,                                                              &
    isParallel,                                                              &
    isSmallExec
USE app_banner, ONLY:                                                       &
    reportApplication
USE umPrintMgr
USE UM_ParCore, ONLY:                                                       &
    mype,                                                                    &
    nproc_max
USE file_manager, ONLY: init_file_manager
USE fort2c_exceptions_interfaces, ONLY:                                      &
    signal_set_verbose, signal_rank, signal_command,                         &
    signal_trap, signal_add_handler, signal_traceback, signal_off
IMPLICIT NONE

CONTAINS

SUBROUTINE appInit(exe)
IMPLICIT NONE
INTEGER, INTENT(IN) :: exe

! Note that the UM gc initialisation is complex due to
! (a) threading
! (b) oasis
! (c) flume
!
! ... So we will only init gcom for small execs, the UM must do
! its own thing before calling appInit.

IF (exe/=exe_UM) THEN
  CALL gc_init(' ',mype, nproc_max)
END IF

! Init output management, parallel executables (except SCM) will redirect
! output, other execs will stick with stdout
IF (exe==exe_UM  .OR.                                                      &
    exe==exe_RCF .OR.                                                      &
    exe==exe_crmstyle_coarse_grid) THEN
  CALL umPrintSetTarget()
END IF

! Init file manager
CALL init_file_manager()

! Init app registry
CALL setApplicationDesc(exe)
CALL umPrintLoadOptions()
CALL umPrintSetLevel()

! Init exception handling
CALL umSetApplicationExceptions(exe)
CALL reportApplication()

! Initialise IO
CALL ioInit()

END SUBROUTINE appInit

SUBROUTINE appTerminate()
USE ereport_mod, ONLY: ereport_finalise
USE fort2c_exceptions_interfaces, ONLY: signal_unregister_callbacks
IMPLICIT NONE

CALL ioShutdown()

CALL ereport_finalise()

CALL umPrintFinalise()

IF (getExeType()/=exe_um) THEN
  CALL gc_exit()
END IF

! free memory structures for callbacks added with signal_add_handler()
CALL signal_unregister_callbacks()

END SUBROUTINE appTerminate

! Abort handler for um_abort_mod to call
SUBROUTINE gcom_signal_abort(errcode)

USE um_types, ONLY: integer64
USE UM_ParCore, ONLY: mype, nproc_max
USE fort2c_exceptions_interfaces, ONLY: signal_controlled_exit

IMPLICIT NONE 

INTEGER, INTENT(IN) :: errcode

! DEPENDS ON: exceptions
CALL signal_controlled_exit(INT(errcode, integer64))

CALL GC_Abort(mype, nproc_max, "um_abort called")

END SUBROUTINE gcom_signal_abort

! turn on signal handling accordingly
SUBROUTINE umSetApplicationExceptions(exe)

USE um_abort_mod, ONLY: set_abort_handler

IMPLICIT NONE
INTEGER, INTENT(IN)     :: exe
INTEGER(KIND=integer64) :: trapping_option
CHARACTER(LEN=256)      :: command

! All of the signal_* C interfaces use int64_t for integer arguments.

IF ( exe==exe_UM .OR. exe==exe_RCF ) THEN
  trapping_option=signal_traceback ! no core, but traceback
ELSE
  trapping_option=signal_off       ! off, completely
END IF

IF (PrintStatus>=PrDiag) CALL signal_set_verbose()

! Tell the signal handler the mpi rank, so that it can tag output.
CALL signal_rank(INT(mype,integer64))

! Tell the signal handler the program name. Not all fortran can do this
! however it is not essential.
CALL get_command_argument(0,command)
CALL signal_command(TRIM(command),INT(LEN_TRIM(command),integer64))

! Set up OS exception handling
! Must be in a serial region because some c calls used are not reentrant.
CALL signal_trap(trapping_option)

    ! Note that registered handlers will be called both on applciation
    ! exceptions (if specified by signal_trap), and on 'clean' shutdowns
    ! via ereport (or umPrintError if from umPrintMgr) - all exit should go
    ! through one of these.
CALL signal_add_handler(umPrintExceptionHandler)

! Set up abort function
CALL set_abort_handler(gcom_signal_abort)

END SUBROUTINE umSetApplicationExceptions
END MODULE UM_Config
