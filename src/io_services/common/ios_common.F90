! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: C96

! Module defining interfaces that might be called by client or server tasks

MODULE IOS_Common
USE IOS_types, ONLY:     &
    IOS_metadata_type, IOS_status
USE IOS_constants, ONLY: &
   IOS_Action_strings, IOS_Action_strlen,                                    &
   IOS_Action_Finish, IOS_Action_Process, IOS_Action_StashWritePPData,       &
   IOS_Action_StashWriteDumpData, IOS_Action_DumpInitModel, IOS_Num_Actions, &
   IOS_Action_Assign_Unit, IOS_Action_StashInitModel, IOS_Action_Config,     &
   IOS_Action_Open, IOS_Action_Close, IOS_Action_Sync,                       &
   IOS_Action_Sync_Barrier, IOS_BytesPerWord64, IOS_BytesPerWord32
USE IOS_communicators
USE missing_data_mod, ONLY: imdi
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE IOS_print_mgr, ONLY: &
    IOS_print, IOS_message

USE um_types, ONLY:                                                            &
  integer_omp

IMPLICIT NONE


! params/vars  for dr_hook
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

INTEGER :: file_op_pseudo_unit

INTEGER :: IOS_Server_Groups=0
INTEGER :: IOS_Tasks_Per_Server=imdi  ! namelist item

! Sequencing: These items track the order of operations to ensure
! communications are always processed in order. There is a different
! sequence for each IO server, that the client must track.

INTEGER, PARAMETER :: IOS_Sequence_Max = 100000
INTEGER, POINTER   :: IOS_Sequence_ID(:,:)

! Unit ranges allowed (MUST conform to the equivalent variables in io module
! and portio)
INTEGER, PARAMETER :: minUnit = 1
INTEGER, PARAMETER :: maxUnit = 300

! Internally defined I/O Processor for a unit
INTEGER            :: io_server_for_unit_lookup                              &
    (minUnit:maxUnit)

! Which global ranks in the application do IO
INTEGER, POINTER   :: io_servers(:,:) => NULL()

! Requests for MPI help for an operation.
INTEGER            :: IOS_ReadCompletionRequested
INTEGER            :: IOS_EnqCompletionRequested

! If I am an io server
LOGICAL            :: l_io_server = .FALSE.
LOGICAL            :: l_io_leader = .FALSE.

! do we need to get lock?
LOGICAL, PARAMETER :: needing_lock   = .TRUE.
LOGICAL, PARAMETER :: no_lock_needed = .FALSE.

! File status
TYPE(ios_status)   :: IOS_unit_status

! Start time
REAL               :: IOS_start_time

! WordLengths
INTEGER            :: IOS_BytesPerReal
INTEGER            :: IOS_BytesPerInteger

!Control the behaviour of IOS
!  Non-namelist items
INTEGER            :: threading_model
LOGICAL            :: use_blocking_recvs       = .FALSE. ! for debugging
LOGICAL            :: serialize_all_ops        = .FALSE. ! for debugging
!  Namelist items
INTEGER            :: IOS_backoff_interval     = imdi
INTEGER            :: IOS_Unit_Alloc_Policy    = imdi
INTEGER            :: IOS_concurrency          = imdi
INTEGER            :: IOS_concurrency_max_mem  = imdi
INTEGER            :: IOS_Timeout              = imdi
!            How long we sit and wait before giving up and aborting the model

INTEGER(KIND=integer_omp) ::                                                   &
  IOS_num_threads = imdi ! if different from model

LOGICAL            :: IOS_use_helpers

LOGICAL            :: IOS_local_ro_files       = .FALSE.
LOGICAL            :: IOS_acquire_model_prsts  = .FALSE.
LOGICAL            :: IOS_serialise_mpi_calls  = .FALSE.
LOGICAL            :: IOS_thread_0_calls_mpi   = .FALSE.
LOGICAL            :: IOS_Lock_Meter           = .FALSE.
LOGICAL            :: IOS_Enable_mpiio         = .FALSE.
LOGICAL            :: IOS_RelayToSlaves        = .TRUE.
INTEGER            :: IOS_Decomp_Model         = 0

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='IOS_COMMON'

CONTAINS

LOGICAL FUNCTION assert_client()
IMPLICIT NONE

INTEGER                         :: ErrorCode = 99
CHARACTER (LEN=*), PARAMETER    :: RoutineName =                           &
    'IOS_Common:ASSERT_CLIENT'

assert_client=.TRUE.
IF (l_io_server) THEN
  assert_client=.FALSE.
  CALL IOS_ereport( RoutineName, ErrorCode,                                &
      'ASSERT_CLIENT FAILURE: ROUTINE CALLED BY IO SERVER PROCESS' )
END IF
END FUNCTION assert_client

LOGICAL FUNCTION assert_server()
IMPLICIT NONE

INTEGER                         :: ErrorCode = 99
CHARACTER (LEN=*), PARAMETER    ::                                         &
    RoutineName = 'IOS_COMMON:ASSERT_SERVER'
assert_server=.TRUE.
IF (.NOT. l_io_server) THEN
  assert_server=.FALSE.
  CALL IOS_ereport( RoutineName, ErrorCode,                                &
      'ASSERT_SERVER FAILURE: ROUTINE CALLED BY IO CLIENT PROCESS')
END IF
END FUNCTION assert_server

SUBROUTINE IOS_ereport(r,code,m,md,UNIT)
USE ereport_mod, ONLY: ereport
IMPLICIT NONE

INTEGER, INTENT(IN)               :: code
TYPE(IOS_metadata_type), OPTIONAL :: md
INTEGER, OPTIONAL                 :: UNIT
CHARACTER(LEN=*), INTENT(IN)      :: r
CHARACTER(LEN=*), INTENT(IN)      :: m
INTEGER                           :: lcode
CHARACTER (LEN=*), PARAMETER      ::                                       &
    RoutineName = 'IOS_COMMON:IOS_EREPORT'

lcode=code
WRITE(IOS_message,*)'-------- IOS ERROR REPORT ---------------'
CALL IOS_print(IOS_message,src='ios_common')
IF (PRESENT(UNIT)) THEN
!$OMP CRITICAL(internal_write)
  WRITE(IOS_message,'(A,I8)')'Problem with unit ',UNIT
!$OMP END CRITICAL(internal_write)
  CALL IOS_print(IOS_message,src='ios_common')
ELSE IF (PRESENT(md)) THEN
  IF (md%UNIT > 0) THEN
!$OMP CRITICAL(internal_write)
    WRITE(IOS_message,'(A,I8)')'Problem with unit ',md%UNIT
!$OMP END CRITICAL(internal_write)
    CALL IOS_print(IOS_message,src='ios_common')
  END IF
END IF

IF (PRESENT(md))CALL ios_report_md(md)
CALL ereport(r,lcode,m)

END SUBROUTINE IOS_ereport

FUNCTION IOS_ActionName(i) RESULT(r)
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE
INTEGER, INTENT(IN)              :: i
INTEGER, PARAMETER               :: warn_code = -1
CHARACTER(LEN=IOS_Action_Strlen) :: r
CHARACTER(LEN=errormessagelength):: error_message
CHARACTER(LEN=*), PARAMETER      :: routinename="IOS_Common :: IOS_ActionName"


IF (i<0 .OR. i>IOS_Num_Actions) THEN
! This is an action we don't recongnise. We will throw a warning and return
! an action name that tells us of a problem, in case we can continue
!$OMP CRITICAL(internal_write)
  WRITE(error_message,'(A,I0)')'IOS Action is outside valid range: ',i
  WRITE(r,'(A)')'BAD_ACTION'
!$OMP END CRITICAL(internal_write)
  CALL IOS_Ereport(routinename, warn_code, error_message)
ELSE
!$OMP CRITICAL(internal_write)
  WRITE(r,'(A)')IOS_Action_strings((i-1)*                                  &
      IOS_Action_Strlen+1:i*IOS_Action_Strlen)
!$OMP END CRITICAL(internal_write)
END IF
END FUNCTION IOS_ActionName

SUBROUTINE IOS_Report_MD(md)
IMPLICIT NONE

TYPE(IOS_metadata_type), INTENT(IN)  :: md
CALL IOS_print(   '------------ IOS METADATA REPORT -----------',          &
    src='ios_common')
!$OMP CRITICAL(internal_write)
WRITE(IOS_message,'(A,I0)')'   md%action=',md%ACTION
!$OMP END CRITICAL(internal_write)
CALL IOS_print(IOS_message,src='ios_common')
IF (md%ACTION >= 1 .AND. md%ACTION <= IOS_Num_Actions) THEN
  WRITE(IOS_message,'(A,A)')'   md%action is ',                            &
      TRIM(IOS_ActionName(md%ACTION))
  CALL IOS_print(IOS_message,src='ios_common')
ELSE
  CALL IOS_print('   md%action is INVALID',src='ios_common')
END IF
!$OMP CRITICAL(internal_write)
WRITE(IOS_message,'(A,I0)')'   md%unit=',md%UNIT
!$OMP END CRITICAL(internal_write)
CALL IOS_print(IOS_message,src='ios_common')
!$OMP CRITICAL(internal_write)
WRITE(IOS_message,'(A,I0)')'   md%name_length=',md%name_length
!$OMP END CRITICAL(internal_write)
CALL IOS_print(IOS_message,src='ios_common')
!$OMP CRITICAL(internal_write)
WRITE(IOS_message,'(A,I0)')'   md%subtype=',md%subtype
!$OMP END CRITICAL(internal_write)
CALL IOS_print(IOS_message,src='ios_common')
!$OMP CRITICAL(internal_write)
WRITE(IOS_message,'(A,I0)')'   md%delete=',md%delete
!$OMP END CRITICAL(internal_write)
CALL IOS_print(IOS_message,src='ios_common')
!$OMP CRITICAL(internal_write)
WRITE(IOS_message,'(A,I0)')'   md%data_size=',md%data_size
!$OMP END CRITICAL(internal_write)
CALL IOS_print(IOS_message,src='ios_common')
!$OMP CRITICAL(internal_write)
WRITE(IOS_message,'(A,I0)')'   md%address=',md%address
!$OMP END CRITICAL(internal_write)
CALL IOS_print(IOS_message,src='ios_common')
!$OMP CRITICAL(internal_write)
WRITE(IOS_message,'(A,I0)')'   md%client=',md%client
!$OMP END CRITICAL(internal_write)
CALL IOS_print(IOS_message,src='ios_common')
!$OMP CRITICAL(internal_write)
WRITE(IOS_message,'(A,I0)')'   md%handle=',md%handle
!$OMP END CRITICAL(internal_write)
CALL IOS_print(IOS_message,src='ios_common')
!$OMP CRITICAL(internal_write)
WRITE(IOS_message,'(A,I0)')'   md%Originating_Slot=',md%Originating_Slot
!$OMP END CRITICAL(internal_write)
CALL IOS_print(IOS_message,src='ios_common')
WRITE(IOS_message,'(A,A)') '   md%string=',md%string
CALL IOS_print(IOS_message,src='ios_common')
CALL IOS_print(   '------------ IOS METADATA REPORT -----------',          &
    src='ios_common')
END SUBROUTINE IOS_Report_MD

LOGICAL FUNCTION L_IOS_active()
IMPLICIT NONE

L_IOS_active=ASSOCIATED(io_servers)
END FUNCTION L_IOS_active

INTEGER FUNCTION io_server_for_unit(UNIT)
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE
INTEGER, INTENT(IN)             :: UNIT
INTEGER                         :: errorCode
CHARACTER (LEN=errormessagelength)              :: message
CHARACTER (LEN=*), PARAMETER    ::                                         &
    RoutineName = 'IOS_COMMON:IO_SERVER_FOR_UNIT'

errorCode=99

IF (.NOT. l_IOS_Active()) THEN
  CALL IOS_ereport(RoutineName,errorCode,'IOS appears deactivated')
END IF

IF (UNIT == file_op_pseudo_unit) THEN
  ! This file unit must go to the lead IO server
  io_server_for_unit = io_servers(1,1)
ELSE

  IF (UNIT < minUnit .OR. UNIT > maxUnit) THEN
!$OMP CRITICAL(internal_write)
    WRITE(message,'(A,I8,A,I3,A,I3,A)')                                    &
        'Supplied unit (',UNIT,') outside allowed range [',                &
        minUnit,' - ',maxUnit,']'
!$OMP END CRITICAL(internal_write)
    CALL IOS_ereport(RoutineName,errorCode,message)
  END IF

  io_server_for_unit = io_server_for_unit_lookup(UNIT)

END IF

END FUNCTION io_server_for_unit


! Returns the server rank ID (i.e. 1..numServers) from
! a global rank

FUNCTION ioServerNo(globRank) RESULT(server)
IMPLICIT NONE

INTEGER, INTENT(IN)              :: globRank
INTEGER                          :: server
INTEGER                          :: i
INTEGER                          :: j
REAL(KIND=jprb)                  :: zhook_handle
CHARACTER (LEN=*), PARAMETER    ::                                         &
    RoutineName = 'IOSERVERNO'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

server=-1

oloop: DO j=1,IOS_tasks_per_server
  iloop: DO i=1,IOS_Server_Groups
    IF (globRank==io_servers(i,j)) THEN
      server=i
      EXIT oloop
    END IF
  END DO iloop
END DO oloop

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION ioServerNo

INTEGER FUNCTION ioServerRank(globRank) RESULT(theRank)
IMPLICIT NONE

INTEGER, INTENT(IN)              :: globRank
INTEGER                          :: i
INTEGER                          :: j
REAL(KIND=jprb)                  :: zhook_handle
CHARACTER (LEN=*), PARAMETER     ::                                        &
    RoutineName = 'IOSERVERRANK'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

theRank=-1

oloop: DO j=1,IOS_tasks_per_server
  iloop: DO i=1,IOS_Server_Groups
    IF (globRank==io_servers(i,j)) THEN
      theRank=j-1
      EXIT oloop
    END IF
  END DO iloop
END DO oloop

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION ioServerRank

! Returns true if the given action is collective over the whole
! parallel IO Server. Simple convenient enumeration - must match
! other code.
LOGICAL FUNCTION IOS_FullTeamNeeded(ACTION)
IMPLICIT NONE
INTEGER, INTENT(IN)              :: ACTION

IOS_FullTeamNeeded=.FALSE.
IF (                                                                       &
    ACTION == IOS_Action_Finish .OR.                                       &
    ACTION == IOS_Action_Process .OR.                                      &
    ACTION == IOS_Action_StashWritePPData .OR.                             &
    ACTION == IOS_Action_StashWriteDumpData .OR.                           &
    ACTION == IOS_Action_DumpInitModel .OR.                                &
    ACTION == IOS_Action_Assign_Unit .OR.                                  &
    ACTION == IOS_Action_StashInitModel .OR.                               &
    ACTION == IOS_Action_Config .OR.                                       &
    IOS_enable_mpiio .AND.                                                 &
    ( ACTION == IOS_Action_Open .OR.                                       &
      ACTION == IOS_Action_Close .OR.                                      &
      ACTION == IOS_Action_Sync .OR.                                       &
      ACTION == IOS_Action_Sync_Barrier )                                  &
    ) THEN
  IOS_FullTeamNeeded=.TRUE.
END IF

END FUNCTION IOS_FullTeamNeeded

INTEGER FUNCTION IOS_Metadata_Receivers()
IMPLICIT NONE

IF (IOS_RelayToSlaves) THEN
  IOS_Metadata_Receivers=1
ELSE
  IOS_Metadata_Receivers=IOS_tasks_per_server
END IF

END FUNCTION IOS_Metadata_Receivers

END MODULE IOS_Common
