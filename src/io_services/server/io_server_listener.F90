! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: C96

MODULE IO_Server_Listener

USE IOS_Common
USE IOS_Constants
USE IOS_types,              ONLY:                                           &
    IOS_metadata_type,                                                       &
    IOS_State_Size
USE IOS_Queue_Mod
USE IOS_Stash_common
USE IOS_Stash_Server,       ONLY:                                           &
    IOS_Stash_Server_Init_Recvs,                                             &
    levels_in_pack
USE IOS_Model_Geometry,     ONLY:                                           &
    IOS_Server_Geometry_Init,                                                &
    getMaxFieldDomain,                                                       &
    atm_numprocs
USE IOS_Server_Coupler, ONLY:                                               &
    IOS_Server_Coupler_Init,                                                 &
    procs
USE mpl
USE yomhook,                ONLY:                                           &
    lhook,                                                                   &
    dr_hook
USE parkind1,               ONLY:                                           &
    jprb,                                                                    &
    jpim
USE IOS_print_mgr, ONLY:                                                    &
    IOS_print,                                                               &
    IOS_print_flush,                                                         &
    IOS_Verbosity,                                                           &
    IOS_message,                                                             &
    IOS_PrStatus_Debug,                                                      &
    IOS_PrStatus_Diag,                                                       &
    IOS_PrStatus_Oper

USE errormessagelength_mod, ONLY: errormessagelength

USE get_wallclock_time_mod, ONLY: get_wallclock_time
 
IMPLICIT NONE


CHARACTER (LEN=errormessagelength), PRIVATE :: IOS_listener_message
INTEGER                      :: timestep
! params/vars  for dr_hook
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='IO_SERVER_LISTENER'

CONTAINS

LOGICAL FUNCTION assert_permitted_client(metadata)
IMPLICIT NONE

TYPE(IOS_metadata_type),INTENT(IN)      :: metadata
IF ( metadata%client /= 0 )                                                &
    CALL IOS_ereport( 'IOS_listener',999,                                  &
    'client rank !=0 not permitted',md=metadata)
assert_permitted_client=.TRUE.
END FUNCTION assert_permitted_client

SUBROUTINE IOS_SlaveListener()
USE IOS_Comms, ONLY:                                                      &
    IOS_Bcast
USE fort2c_interfaces, ONLY: um_sleep
IMPLICIT NONE

CHARACTER (LEN=*), PARAMETER :: RoutineName = 'IOS_SLAVELISTENER'
TYPE(IOS_metadata_type)      :: metadata
REAL(KIND=jprb)              :: zhook_handle1
REAL(KIND=jprb)              :: zhook_handle2

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle1)
CALL IOS_print( 'Info: SlaveListener: Process started',                    &
    src='io_Server_listener')
timestep=0
metadata%ACTION = IOS_Action_Unset
! Keep listening until we get told to stop
DO WHILE (metadata%ACTION /= IOS_Action_Finish)
  CALL IOS_Bcast(metadata)
  CALL IOS_process_normal(metadata,0)
END DO

CALL IOS_print( 'Info: SlaveListener: Exited listen loop',                 &
    src='io_Server_listener')

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle2)

IF ( IOS_thread_0_calls_mpi ) THEN
  DO WHILE (IOS_getQueueItems() > 0)
    CALL IOS_Aux_Thread1_Assist()
    CALL um_sleep(IOS_backoff_interval)
  END DO
END IF

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle2)

CALL IOS_print( 'Info: SlaveListener: Process closing',                    &
    src='io_Server_listener')
IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle1)

END SUBROUTINE IOS_SlaveListener

SUBROUTINE IOS_Listener()
USE IOS_Comms, ONLY:                                                      &
    IOS_Bcast,                                                             &
    acquire_lock_mpi,                                                      &
    release_lock_mpi
USE fort2c_interfaces, ONLY: um_sleep
USE IOS_Constants, ONLY: IOS_Action_Strlen
IMPLICIT NONE


CHARACTER (LEN=*), PARAMETER :: RoutineName = 'IOS_LISTENER'
CHARACTER (LEN=IOS_Action_Strlen) :: action_name
INTEGER                      :: ierror
INTEGER                      :: STATUS(mpl_status_size)
INTEGER                      :: command_request !an MPL request obj
INTEGER                      :: express_command_request !an MPL request obj
INTEGER                      :: reserveSpace
TYPE(IOS_metadata_type)      :: metadata
TYPE(IOS_metadata_type)      :: express_metadata
TYPE(IOS_node_type), POINTER :: node
INTEGER, ALLOCATABLE         :: idata(:)
INTEGER                      :: mdTag
INTEGER                      :: exTag
INTEGER                      :: strTag
INTEGER                      :: payloadTag
INTEGER                      :: tagOffset
LOGICAL                      :: have_message
LOGICAL                      :: ok
LOGICAL                      :: process_express_item
LOGICAL                      :: process_normal_item
INTEGER                      :: sequenceID
REAL                         :: t
REAL(KIND=jprb):: zhook_handle1
REAL(KIND=jprb):: zhook_handle2


IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle1)

CALL IOS_print( 'Info: Listener: Process started',src='io_server_listener')
timestep=0
sequenceID=1;

metadata%ACTION = IOS_Action_Unset

! Start these flags off as true in order to post their initial receives
process_express_item=.TRUE.
process_normal_item=.TRUE.

! Keep listening until we get told to stop
DO WHILE (metadata%ACTION /= IOS_Action_Finish)
  NULLIFY(node)

  ! Calculate tags for this loop trip
  tagOffset=MOD(sequenceID,IOS_Request_Tag_Gap)
  mdTag      =IOS_Request_Tag_Base    +tagOffset
  exTag      =IOS_Request_Tag_Express
  strTag     =IOS_Request_Tag_Str     +tagOffset
  payloadTag =IOS_Request_Tag_Payload +tagOffset

  !----------------------------------------------------------------
  ! Receive metadata. This is complicated because this thread may
  ! not block. If he does, he may not notice that thread 1 is
  ! asking for assistance in the case that a single threaded MPI
  ! is in use, which would result in a deadlock.
  !----------------------------------------------------------------
  CALL acquire_lock_mpi()

  ! If we have not got an outstanding receive for an express item post one
  IF (process_express_item) THEN
    CALL MPL_iRecv(express_metadata , IOS_md_len, mpl_integer,             &
        mpl_any_source, exTag, Global_Comm,                                &
        express_command_request, ierror)
    process_express_item=.FALSE.
  END IF

  ! If we have not got an outstanding receive for an regular item post one
  IF (process_normal_item) THEN
    CALL MPL_iRecv(metadata , IOS_md_len, mpl_integer,                     &
        mpl_any_source, mdTag, Global_Comm,                                &
        command_request, ierror)
    process_normal_item=.FALSE.
  END IF

  CALL release_lock_mpi()

  ! Now wait for something to arrive
  have_message=.FALSE.
  IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
!$OMP CRITICAL(internal_write)
    WRITE(IOS_message,'(A,I8)')'Info: Listener: Waiting for operation ',   &
        sequenceID
!$OMP END CRITICAL(internal_write)
    CALL IOS_print(IOS_message,src='io_server_listener')
  END IF
  DO WHILE (.NOT. have_message)
    CALL acquire_lock_mpi()
    ! First check the express tag (note we don't release MPI - its express!)
    CALL MPL_Test(express_command_request, have_message, STATUS, ierror)

    ! If we have an express message then we need to set the flag to post another
    IF (have_message) THEN
      process_express_item=.TRUE.

      ! Otherwise check on normal messages (and release MPI)
    ELSE
      CALL MPL_Test(command_request, have_message, STATUS, ierror)
      CALL release_lock_mpi()

      ! If we have a normal message then we need to set the flag to post another
      ! and update the sequence ID
      IF (have_message) THEN
        process_normal_item=.TRUE.
        sequenceID=sequenceID+1
      ELSE
        IF ( IOS_thread_0_calls_mpi ) THEN
          !Check to see if thread 1 needs attention
          CALL IOS_Aux_Thread1_Assist()
        END IF
      END IF
    END IF
    IF ( .NOT. have_message ) CALL um_sleep(IOS_backoff_interval)
  END DO

  IF (process_express_item) THEN
    IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
      action_name=TRIM(IOS_ActionName(express_metadata%action))
!$OMP CRITICAL(internal_write)
      WRITE(IOS_message,'(A,A)')                                           &
          'Info: Listener: Received High Priority Action: ',               &
          action_name
!$OMP END CRITICAL(internal_write)
      CALL IOS_print(IOS_message,src='io_server_listener')
      CALL IOS_Print_Flush()
    END IF
  ELSE
    IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
      action_name=TRIM(IOS_ActionName(metadata%action))
!$OMP CRITICAL(internal_write)
      WRITE(IOS_message,'(A,A)')                                           &
          'Info: Listener: Received a transaction: ',                      &
          action_name
!$OMP END CRITICAL(internal_write)
      CALL IOS_print(IOS_message,src='io_server_listener')
      CALL IOS_Print_Flush()
    END IF
  END IF

  ok=assert_permitted_client(metadata)
  ! The client filled in the detail with his "model_rank" id.
  ! We may need to send back data, so we need to switch
  ! that to his global ID after we checked it.
  metadata%client = STATUS(mpl_source)

  ! For debugging we may want to pause here until
  ! the writer isn't busy this will allow everything
  ! to be processed in natural "1 in - 1 out" sequence
  ! with little activity overlap between threads.
  IF ( serialize_all_ops ) THEN
    CALL IOS_WaitForQueueDrain()
  END IF

  ! Now we have received something give it to the process routine
  ! to decide how to respond.
  IF (process_express_item) THEN
    SELECT CASE (express_metadata%ACTION)

    CASE (IOS_Action_LoadStatus)

      !--------------------------------------------------------------
      ! Process the request for my loading
      !--------------------------------------------------------------

      ALLOCATE( idata ( loadBalanceDataSize ) )
      idata(loadBalanceQueueLen)  = IOS_getQueueItems()
      idata(loadBalanceQueueData) = IOS_getQueuePayload()
      CALL MPL_Send(idata,                                                 &
          loadBalanceDataSize*IOS_tuperword64,                             &
          IOS_tutype, express_metadata%client,                             &
          exTag, Global_Comm, ierror)
      DEALLOCATE( idata  )

    CASE DEFAULT
!$OMP CRITICAL(internal_write)
      WRITE(IOS_listener_message,'(A,I3,A)') 'Express Action ',            &
          express_metadata%ACTION,' not recognised!'
!$OMP END CRITICAL(internal_write)
      CALL IOS_ereport( RoutineName, 60, IOS_listener_message,             &
          md=express_metadata )

    END SELECT

    ! Release the MPI lock from when we checked for the express message
    CALL release_lock_mpi()
  ELSE
    IF (IOS_FullTeamNeeded(metadata%ACTION) .AND. IOS_RelayToSlaves) THEN
      IF (IOS_Verbosity>=IOS_PrStatus_Debug) THEN
        CALL IOS_print('Info: Listener: Sending '//                        &
            TRIM(IOS_ActionName(metadata%ACTION)) //' to slaves',          &
            src='io_Server_listener')
      END IF
      CALL IOS_Bcast(metadata)
    END IF
    CALL IOS_process_normal(metadata,tagOffset)
  END IF

END DO

! Cancel the preemption receive request, so that we don't have any
! pending receives at shutdown.
CALL MPL_Cancel(express_command_request, ierror)

CALL IOS_print( 'Info: Listener: Exited listen loop',                      &
    src='io_Server_listener')

! We cannot just return at this point. In MPI models where only
! thread zero can call MPI, the writer thread may need to proxy MPI
! commands via this thread. We need to hang around and help him
! clear the queue.
IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle2)

IF ( IOS_thread_0_calls_mpi ) THEN
  DO WHILE (IOS_getQueueItems() > 0)
    CALL IOS_Aux_Thread1_Assist()
    CALL um_sleep(IOS_backoff_interval)
  END DO
END IF

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle2)

! Final output
CALL IOS_print( 'Info: Listener: Process closing',                         &
    src='io_Server_listener')
IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle1)

RETURN
END SUBROUTINE IOS_Listener



SUBROUTINE IOS_process_normal(metadata, tagOffset)
USE IOS_Comms, ONLY:                                                      &
    IOS_Bcast,                                                             &
    IOS_SoftSync,                                                          &
    acquire_lock_mpi,                                                      &
    release_lock_mpi
USE io_constants, ONLY: ioFileTypeMPIIO
USE ios_common, ONLY: needing_lock
IMPLICIT NONE
TYPE(IOS_metadata_type)      :: metadata
INTEGER, INTENT(IN)          :: tagOffset
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'IOS_PROCESS_NORMAL'
INTEGER, ALLOCATABLE         :: idata(:)
TYPE(IOS_node_type), POINTER :: node
INTEGER                      :: payload_request
INTEGER                      :: reserveSpace
INTEGER                      :: iError
INTEGER                      :: strTag
INTEGER                      :: payloadTag
INTEGER                      :: STATUS(mpl_status_size)
REAL                         :: theTime
REAL(KIND=jprb)              :: zhook_handle1

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle1)

strTag     =IOS_Request_Tag_Str     +tagOffset
payloadTag =IOS_Request_Tag_Payload +tagOffset

SELECT CASE (metadata%ACTION)

CASE (IOS_Action_Assign_Unit)

  !--------------------------------------------------------------
  ! Process the assignment of a unit
  !--------------------------------------------------------------

  node => make_new_node()
  node%metadata = metadata
  CALL IOS_Put_Node_In_Queue(node)

  !--------------------------------------------------------------
  ! Process Writes (64 bit/stash)
  !--------------------------------------------------------------

CASE (IOS_Action_Write64,                                                  &
    IOS_Action_Write_PP_Prepacked)

  CALL IOS_WaitForFreeSpace(metadata%data_size*                            &
      IOS_BytesPerWord64)

  node => make_new_node()
  node%metadata = metadata

  CALL IOS_Increment_Queue_Length                                          &
      (metadata%data_size*IOS_BytesPerWord64,needing_lock)
  ALLOCATE( node%real_data( metadata%data_size) )

  CALL acquire_lock_mpi()
  CALL MPL_Recv(node%real_data, metadata%data_size, mpl_real,              &
      metadata%client, payloadTag , Global_Comm,                           &
      STATUS, ierror)
  CALL release_lock_mpi()
  ! Put data on Queue
  CALL IOS_Put_Node_In_Queue(node)

  !--------------------------------------------------------------
  ! Process Writes (32 bit)
  !--------------------------------------------------------------

CASE (IOS_Action_Write32)

  CALL IOS_WaitForFreeSpace(metadata%data_size*                            &
      IOS_BytesPerWord32)

  node => make_new_node()
  node%metadata = metadata

  CALL IOS_Increment_Queue_Length                                          &
      (metadata%data_size*IOS_BytesPerWord32,needing_lock)
  ALLOCATE( node%real32_data( metadata%data_size) )

  CALL acquire_lock_mpi()
  CALL MPL_Recv(node%real32_data, metadata%data_size, mpl_real,            &
      metadata%client, PayloadTag, Global_Comm,                            &
      STATUS, ierror)
  CALL release_lock_mpi()
  ! Put data on Queue
  CALL IOS_Put_Node_In_Queue(node)

  !--------------------------------------------------------------
  ! Process Dump initialisation
  !--------------------------------------------------------------

CASE ( IOS_Action_DumpInitModel )

  CALL IOS_WaitForFreeSpace(metadata%data_size*                            &
      IOS_BytesPerWord64)

  node => make_new_node()
  node%metadata = metadata

  CALL IOS_Increment_Queue_Length(                                         &
      metadata%data_size*IOS_BytesPerInteger,needing_lock)
  ALLOCATE( node%integer_data( metadata%data_size) )

  IF (model_rank==0 .OR. .NOT. IOS_RelayToSlaves) THEN
    CALL acquire_lock_mpi()
    CALL MPL_Recv(node%integer_data, metadata%data_size, mpl_integer,      &
        metadata%client, PayloadTag, Global_Comm,                          &
        STATUS, ierror)
    CALL release_lock_mpi()
  END IF

  IF (IOS_RelayToSlaves) CALL IOS_Bcast(node%integer_data)

  ! Put data on Queue
  CALL IOS_Put_Node_In_Queue(node)

  !--------------------------------------------------------------
  ! Process Reads
  !--------------------------------------------------------------

CASE (IOS_Action_Read64,                                                   &
    IOS_Action_Read32,                                                     &
    IOS_Action_Read32_Integer,                                             &
    IOS_Action_Read64_Integer)
  node => make_new_node()
  node%metadata   = metadata
  node%payloadTag = payloadTag
  CALL IOS_Put_Node_In_Queue(node)

  !--------------------------------------------------------------
  ! Process Enquiry
  !--------------------------------------------------------------

CASE (IOS_Action_Enquire)
  node => make_new_node()
  node%metadata   = metadata
  node%payloadTag = payloadTag
  CALL IOS_Put_Node_In_Queue(node)

  !--------------------------------------------------------------
  ! Process A unit sync OP
  !--------------------------------------------------------------

CASE (IOS_Action_Sync)
  CALL IOS_Put_Metadata_In_Queue(metadata)

CASE (IOS_Action_Sync_Barrier)
  CALL IOS_Put_Metadata_In_Queue(metadata)
  IF ( IOS_thread_0_calls_mpi ) THEN
    CALL IOS_WaitForQueueDrain()
    ! Send 1 integer to the client to indicate completion
    ! In the MPI-IO case, the client is only expecting data from
    ! the lead IO task. In the POSIX case, only the lead IO task 
    ! recieves this action. 
    IF ( model_rank == 0 )                                                 &
         CALL MPL_Send(metadata,                                           &
              1, mpl_integer, metadata%client,                             &
              node%payloadTag, Global_Comm,                                &
              ierror)
  END IF

  !--------------------------------------------------------------
  ! Process Open/Close
  !--------------------------------------------------------------

CASE (IOS_Action_Open,IOS_Action_Close)
  ! Get extra metadata and put action on queue
   IF (model_rank==0 .OR. & 
        ( IOS_enable_mpiio.AND..NOT. IOS_RelayToSlaves) ) THEN
      CALL acquire_lock_mpi()
      CALL MPL_Recv(metadata%string, IOS_string_max, mpl_character,        &
           metadata%client, strTag, Global_Comm,                           &
           STATUS, ierror)
      CALL release_lock_mpi()
      CALL IOS_Put_Metadata_In_Queue(metadata)
   END IF
   ! If we are opening an MPI-IO file, broadcast the name to the other IO
   ! server tasks and ensure the metadata is queued on all tasks.
   IF ( IOS_enable_mpiio .AND. metadata%address == ioFileTypeMPIIO         &
        .AND. IOS_RelayToSlaves ) THEN 
      CALL IOS_bcast(metadata%string)
      IF ( model_rank /= 0 ) CALL IOS_Put_Metadata_In_Queue(metadata)
   END IF

CASE (IOS_Action_FileOp)
  CALL acquire_lock_mpi()
  CALL MPL_Recv(metadata%string, IOS_string_max, mpl_character,            &
      metadata%client, strTag, Global_Comm,                                &
      STATUS, ierror)
  CALL release_lock_mpi()

  IF ( IOS_thread_0_calls_mpi ) THEN
    IF ( IOS_Verbosity >= IOS_PrStatus_Oper ) THEN
!$OMP CRITICAL(internal_write)
      WRITE(IOS_message,'(A,A)')'Info: Listener: ',                        &
          'Received Release or FileOp, stalling listener task'
!$OMP END CRITICAL(internal_write)
      CALL IOS_print(IOS_message,src='io_Server_listener')
    END IF
    CALL IOS_WaitForQueueDrain()
    IF ( IOS_Verbosity >= IOS_PrStatus_Oper ) THEN
      CALL IOS_print('Info: Listener: Queue drained',                      &
          src='io_Server_listener')
    END IF

  END IF
  CALL IOS_Put_metadata_in_queue(metadata)

  !--------------------------------------------------------------
  ! Process Fence requests
  !--------------------------------------------------------------

CASE (IOS_Action_Fence)
  IF ( IOS_thread_0_calls_mpi ) THEN
    CALL IOS_WaitForQueueDrain()
    CALL IOS_SoftSync(metadata%address)
  END IF
  CALL IOS_Put_metadata_in_queue(metadata)

  !--------------------------------------------------------------
  ! Process Config
  !--------------------------------------------------------------

CASE (IOS_Action_Config)
  CALL IOS_Put_metadata_in_queue(metadata)

  !--------------------------------------------------------------
  ! Process Flush/Setpos
  !--------------------------------------------------------------

CASE (IOS_Action_Flush,                                                    &
    IOS_Action_Setpos,                                                     &
    IOS_Action_StashSetPos )
  CALL IOS_Put_Metadata_In_Queue(metadata)

  !--------------------------------------------------------------
  ! Process Finish
  !--------------------------------------------------------------

CASE (IOS_Action_Finish)
  CALL IOS_Put_Metadata_In_Queue(metadata)

  theTime=get_wallclock_time()
  CALL IOS_print('*******************************************************',&
      src='io_server_listener')
!$OMP CRITICAL(internal_write)
  WRITE(IOS_message,'(A,F10.3,A)')                                         &
      '* IO SERVER:LISTENER RECEIVED FINISH CMD AT ',                      &
      theTime-IOS_Start_Time,' *'
!$OMP END CRITICAL(internal_write)
  CALL IOS_print(IOS_message,src='io_server_listener')
  CALL IOS_print('*******************************************************',&
      src='io_server_listener')

  !--------------------------------------------------------------
  ! Process Start of New Timestep
  !--------------------------------------------------------------

CASE (IOS_Action_Process)
  timestep=timestep+1
  theTime=get_wallclock_time()
  IF ( IOS_Verbosity >= IOS_PrStatus_Oper ) THEN
!$OMP CRITICAL(internal_write)
    WRITE(IOS_message,'(A,I4,A,F10.3)')                                    &
        'Info: Listener: Entering Timestep ',timestep,                     &
        ' at time=',theTime-IOS_start_time
!$OMP END CRITICAL(internal_write)
    CALL IOS_print(IOS_message,src='io_Server_listener')
  END IF
  CALL IOS_Put_Metadata_In_Queue(metadata)

  !-------------------------------------------------------------------
  ! Process Initialising a fixed Header
  !-------------------------------------------------------------------

CASE (IOS_Action_StashInitHeader)
  node => make_new_node()
  node % metadata = metadata
  CALL IOS_Put_Node_In_Queue(node)

  !-------------------------------------------------------------------
  ! Process Setting the content of a fixed Header
  !-------------------------------------------------------------------

CASE (IOS_Action_StashSetHeader)

  CALL IOS_WaitForFreeSpace(metadata % data_size*                          &
      IOS_BytesPerWord64)
  node => make_new_node()
  node % metadata = metadata
  CALL IOS_Increment_Queue_Length(                                         &
      metadata % data_size * IOS_BytesPerInteger,needing_lock)
  ALLOCATE( node % integer_data( metadata % data_size) )

  CALL acquire_lock_mpi()
  CALL MPL_Recv(node % integer_data, metadata % data_size,                 &
      mpl_integer,metadata%client, PayloadTag, Global_Comm,                &
      STATUS, ierror)
  CALL release_lock_mpi()

  CALL IOS_Put_Node_In_Queue(node)

  !-------------------------------------------------------------------
  ! Process Writes of PP lookup data
  !-------------------------------------------------------------------

CASE (IOS_ACTION_StashWritePPLookup,                                       &
    IOS_action_mergepplookup)

  CALL IOS_WaitForFreeSpace(metadata % data_size*                          &
      IOS_BytesPerWord64)
  node => make_new_node()
  node % metadata = metadata
  CALL IOS_Increment_Queue_Length(                                         &
      metadata % data_size * IOS_BytesPerInteger,needing_lock)
  ALLOCATE( node % integer_data( metadata % data_size) )

  CALL acquire_lock_mpi()
  CALL MPL_Recv(node % integer_data, metadata % data_size,                 &
      mpl_integer, metadata%client, PayloadTag,                            &
      Global_Comm, STATUS, ierror)
  CALL release_lock_mpi()

  ! Put data on Queue
  CALL IOS_Put_Node_In_Queue(node)

  !-------------------------------------------------------------------
  ! Process The initialisation of PP lookup buffers
  !-------------------------------------------------------------------

CASE (IOS_Action_StashInitPPLookup)
  CALL IOS_WaitForFreeSpace(metadata % data_size*                          &
      IOS_BytesPerWord64)
  node => make_new_node()
  node % metadata = metadata
  CALL IOS_Increment_Queue_Length(                                         &
      metadata % data_size * IOS_BytesPerInteger,needing_lock)
  ALLOCATE( node % integer_data( metadata % data_size) )

  CALL acquire_lock_mpi()
  CALL MPL_Recv(node % integer_data, metadata % data_size,                 &
      mpl_integer, metadata%client, payloadTag,                            &
      Global_Comm,                                                         &
      STATUS, ierror)
  CALL release_lock_mpi()

  ! Put data on Queue
  CALL IOS_Put_Node_In_Queue(node)

  !-------------------------------------------------------------------
  ! Process The initialisation of model geometry
  !-------------------------------------------------------------------

CASE (IOS_Action_StashInitModel)

  ALLOCATE( idata ( metadata % data_size) )

  IF (model_rank==0 .OR. .NOT. IOS_RelayToSlaves) THEN
    CALL acquire_lock_mpi()
    CALL MPL_Recv(idata, metadata % data_size, mpl_integer,                &
        metadata%client, payloadTag,                                       &
        Global_Comm,                                                       &
        STATUS, ierror)
    CALL release_lock_mpi()
  END IF
  IF (IOS_RelayToSlaves) CALL IOS_Bcast(idata)

  ! Do not put data on Queue, just call the stash server with the details
  ! directly but let us drain the queue first out of sheer paranoia
  ! Just in case the other thread is in the middle of doing something
  ! that we havn't thought carefully about....
  CALL IOS_WaitForQueueDrain()

  CALL acquire_lock_mpi()
  CALL IOS_Server_Geometry_Init(idata)
  CALL release_lock_mpi()
  CALL IOS_Server_Coupler_Init()
  DEALLOCATE(idata)

  !-------------------------------------------------------------------
  ! Process Write PP file data
  !-------------------------------------------------------------------

CASE (IOS_Action_StashWritePPData,                                         &
    IOS_Action_StashWriteDumpData)
  CALL IOS_WaitForFreeSpace(metadata % data_size*                          &
      IOS_BytesPerWord64)

  ! we need to allow for the queue to accommodate the receive
  ! buffers even if we do not allocate them here
  ! otherwise we may see a deadlock at the actual allocation time

  node => make_new_node()
  node % metadata = metadata
  ALLOCATE( node % integer_data( metadata % data_size) )

  IF (model_rank==0 .OR. .NOT. IOS_RelayToSlaves) THEN
    CALL acquire_lock_mpi()
    CALL MPL_Recv(node % integer_data, metadata % data_size,               &
        mpl_integer, metadata%client, payloadTag, Global_Comm,             &
        STATUS, ierror)
    CALL release_lock_mpi()
  END IF

  IF (IOS_RelayToSlaves) CALL IOS_bcast(node % integer_data)

  ! Cheap sanity check that the buffer we received looks like a stash
  ! descriptor block.....
  IF (node%integer_data(1) /= IOS_stash_record_start) THEN
!$OMP CRITICAL(internal_write)
    WRITE(IOS_listener_message,'(A,I0)')                                   &
        'IOS check failure, 1st word of control buffer is not ',           &
        IOS_stash_record_start
!$OMP END CRITICAL(internal_write)
    CALL IOS_ereport( RoutineName, 999, IOS_listener_message,              &
        md=metadata )
  END IF

  reserveSpace=                                                            &
      getMaxFieldDomain()*                                                 &
      levels_in_pack(node)*                                                &
      procs
  CALL IOS_WaitForFreeSpace(reserveSpace*                                  &
      IOS_BytesPerWord64)
  CALL IOS_Increment_Queue_Length(reserveSpace * IOS_BytesPerReal,         &
                                needing_lock)

  CALL IOS_Stash_server_init_recvs(node)
  CALL IOS_Put_Node_In_Queue(node)
  CALL IOS_Increment_Queue_Length(                                         &
      metadata % data_size * IOS_BytesPerInteger,needing_lock)


  !--------------------------------------------------------------
  ! Process any other unrecognised action
  !--------------------------------------------------------------

CASE DEFAULT
!$OMP CRITICAL(internal_write)
  WRITE(IOS_listener_message,'(A,I3,A)') 'Action ',                        &
      metadata%ACTION,' not recognised!'
!$OMP END CRITICAL(internal_write)
  CALL IOS_ereport( RoutineName, 60, IOS_listener_message,                 &
      md=metadata )

END SELECT

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle1)

END SUBROUTINE IOS_process_normal

SUBROUTINE construct_status(UNIT,STATUS)
  ! Note we use the 1=true convention for the status items
  ! Status only knows about open/closed thus far
USE io, ONLY: is_unit_open
IMPLICIT NONE

INTEGER, INTENT(IN) :: UNIT
INTEGER, INTENT(OUT):: STATUS(1)
INTEGER             :: opn

! Is the file open?
STATUS(1)=0
IF (is_unit_open(UNIT))STATUS(1)=1

END SUBROUTINE construct_status

SUBROUTINE IOS_Aux_Thread1_Assist()
USE IOS_Stash_Server,       ONLY:                                         &
    IOS_Stash_Server_Finish_Recvs
IMPLICIT NONE


TYPE(IOS_node_type), POINTER :: q_head
LOGICAL                      :: test_status

NULLIFY(q_head)

!$OMP FLUSH
IF (IOS_ReadCompletionRequested==1) THEN
  !Thread 1 wants me to CALL MPI for him
  CALL IOS_Get_Node_From_Queue(q_head)
  CALL IOS_Aux_Read_Assist(q_head)
  NULLIFY (q_head)
  IOS_ReadCompletionRequested=0
!$OMP FLUSH
ELSE IF (IOS_AsyncCompletionRequested==1) THEN
  !Thread 1 wants me to CALL MPI for him
  CALL IOS_Get_Node_From_Queue(q_head)
  test_status=IOS_Stash_Server_Finish_Recvs(q_head)
  NULLIFY (q_head)
  IF (test_status)IOS_AsyncCompletionRequested=0
!$OMP FLUSH
ELSE IF (IOS_EnqCompletionRequested==1) THEN
  !Thread 1 wants me to CALL MPI for him
  CALL IOS_Get_Node_From_Queue(q_head)
  CALL IOS_Aux_Enq_Assist(q_head)
  NULLIFY (q_head)
  IOS_EnqCompletionRequested=0
!$OMP FLUSH
END IF
END SUBROUTINE IOS_Aux_Thread1_Assist


SUBROUTINE IOS_Aux_Read_Assist(node)
USE IOS_Comms, ONLY:                                                      &
    acquire_lock_mpi,                                                      &
    release_lock_mpi
IMPLICIT NONE

TYPE(IOS_node_type),POINTER  :: node
INTEGER                      :: ierror
REAL(KIND=jprb)              :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='IOS_AUX_READ_ASSIST'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

SELECT CASE (node%metadata%ACTION)

CASE (IOS_Action_Read64)
  CALL acquire_lock_mpi()
  IF (node%metadata%subtype==IOS_Read_Broadcast) THEN
    CALL mpl_bcast(                                                        &
        node%real_data,                                                    &
        node%metadata%data_size*IOS_tuperword64,                           &
        IOS_tutype,                                                        &
        bcast_procs-1,                                                     &
        IOS_BCast_Server_Comm,                                             &
        ierror)
  ELSE
    CALL MPL_Send(node%real_data,                                          &
        node%metadata%data_size*IOS_tuperword64,                           &
        IOS_tutype, node%metadata%client,                                  &
        node%payloadTag, Global_Comm, ierror)
  END IF
  CALL release_lock_mpi()
CASE (IOS_Action_Read32 )
  CALL acquire_lock_mpi()
  IF (node%metadata%subtype==IOS_Read_Broadcast) THEN
    CALL mpl_bcast(                                                        &
        node%real32_data,                                                  &
        node%metadata%data_size*IOS_tuperword32,                           &
        IOS_tutype,                                                        &
        bcast_procs-1,                                                     &
        IOS_BCast_Server_Comm,                                             &
        ierror)
  ELSE
    CALL MPL_Send(node%real32_data,                                        &
        node%metadata%data_size*IOS_tuperword32,                           &
        IOS_tutype, node%metadata%client,                                  &
        node%payloadTag, Global_Comm, ierror)
  END IF
  CALL release_lock_mpi()
CASE (IOS_Action_Read64_Integer)
  CALL acquire_lock_mpi()
  IF (node%metadata%subtype==IOS_Read_Broadcast) THEN
    CALL mpl_bcast(                                                        &
        node%integer_data,                                                 &
        node%metadata%data_size*IOS_tuperword64,                           &
        IOS_tutype,                                                        &
        bcast_procs-1,                                                     &
        IOS_BCast_Server_Comm,                                             &
        ierror)
  ELSE
    CALL MPL_Send(node%integer_data,                                       &
        node%metadata%data_size*IOS_tuperword64,                           &
        IOS_tutype, node%metadata%client,                                  &
        node%payloadTag, Global_Comm, ierror)
  END IF
  CALL release_lock_mpi()
CASE DEFAULT
  CALL IOS_Ereport('IOS_Aux_Read_Assist', 10,                              &
      'Received a non-read operation',                                     &
      md=node%metadata)

END SELECT

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE IOS_Aux_Read_Assist

SUBROUTINE IOS_Aux_Enq_Assist(node)
USE IOS_Comms, ONLY:                                                      &
    acquire_lock_mpi,                                                      &
    release_lock_mpi
IMPLICIT NONE

TYPE(IOS_node_type),POINTER  :: node
INTEGER                      :: ierror
REAL(KIND=jprb)              :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='IOS_AUX_ENQ_ASSIST'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL acquire_lock_mpi()
CALL MPL_Send(IOS_Unit_status,                                             &
    IOS_State_Size, mpl_integer,                                           &
    node%metadata%client, node%payloadTag,                                 &
    Global_Comm, ierror)
CALL release_lock_mpi()

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE IOS_Aux_Enq_Assist

SUBROUTINE IOS_WaitForQueueDrain()
IMPLICIT NONE

LOGICAL                     :: done
REAL(KIND=jprb)             :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='IOS_WAITFORQUEUEDRAIN'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

done=.FALSE.
DO WHILE (.NOT. done)
  IF (.NOT. waitforDrain()) THEN
    !it returned before space was available so offer help
    !Check to see if thread 1 needs attention
    CALL IOS_Aux_Thread1_Assist()
  ELSE
    done=.TRUE.
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE IOS_WaitForQueueDrain

SUBROUTINE IOS_WaitForFreeSpace(sizeNeeded)
IMPLICIT NONE
INTEGER, INTENT(IN)         :: sizeNeeded
LOGICAL                     :: done
REAL(KIND=jprb)             :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='IOS_WAITFORFREESPACE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

done=.FALSE.
DO WHILE (.NOT. done)
  IF (.NOT. HasFreeSpace(sizeNeeded)) THEN
    !it returned before space was available so offer help
    !Check to see if thread 1 needs attention
    CALL IOS_Aux_Thread1_Assist()
  ELSE
    done=.TRUE.
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE IOS_WaitForFreeSpace

END MODULE IO_Server_Listener
