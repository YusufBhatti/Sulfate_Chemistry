! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: C96

! Implements a ring buffer for client side IOS operations
! and manages dispatch and data receipt into elements of that queue

MODULE IOS_Client_Queue

USE IOS_Common, ONLY:                                                        &
    IOS_Sequence_ID, io_servers, L_IOS_Active, Model_Comm, Model_Procs,      &
    model_rank, assert_client, IOS_Unit_Alloc_Policy, IOS_Timeout,           &
    IOS_start_time, io_server_for_unit_lookup, io_server_for_unit,           &
    IOS_concurrency, IOS_concurrency_max_mem, IOS_Server_Groups,             &
    IOS_Tasks_Per_Server, IOS_ereport, IOS_Metadata_Receivers,               &
    IOS_ActionName, ioServerNo, ioServerRank, IOS_RelayToSlaves,             &
    IOS_Report_MD, file_op_pseudo_unit

USE IOS_communicators, ONLY:                                                 &
    global_comm

USE IOS_Constants, ONLY:                                                     &
    IOS_Unit_Alloc_Static, IOS_queue_slot_unused, IOS_Request_Tag_Express,   &
    loadBalanceQueueData, IOS_Action_Sync_Barrier, IOS_err_strlen,           &
    IOS_Unit_Alloc_AtFirstUse, IOS_Unit_Alloc_Dynamic_Rotate, IOS_No_Server, &
    IOS_Action_Assign_Unit, loadBalanceDataSize, IOS_Unit_Alloc_Dynamic_LB,  &
    IOS_tuPerWord32, IOS_Action_StashInitModel, IOS_Action_Write32,          &
    IOS_tukind, IOS_No_Location, IOS_Action_LoadStatus, IOS_BytesPerTU,      &
    IOS_Action_Read32_Integer, IOS_Action_Enquire, IOS_BytesPerWord32,       &
    IOS_Action_Read32, ios_queue_slot_dispatched, IOS_Action_DumpInitModel,  &
    IOS_Action_Read64, IOS_Action_Write_PP_Prepacked, IOS_md_len,            &
    IOS_Action_StashWritePPLookup, IOS_Read_Broadcast, IOS_Action_Write64,   &
    IOS_Action_StashSetHeader, IOS_Action_StashWritePPData, IOS_String_Max,  &
    IOS_Action_Read64_Integer, IOS_Request_Tag_Gap, IOS_Request_Tag_Payload, &
    IOS_queue_slot_initialized, IOS_Action_StashInitPPLookup, IOS_tutype,    &
    IOS_queue_slot_partfilled, loadBalanceQueueLen, IOS_BytesPerWord64,      &
    IOS_Action_Finish, IOS_Request_Tag_Str, IOS_Action_StashWriteDumpData,   &
    IOS_tuperword64, IOS_Request_Tag_Base, IOS_Action_MergePPLookup,         &
    IOS_Action_StashInitHeader

USE IOS_types, ONLY:                                                         &
    IOS_metadata_type, IOS_async_object,                                     &
    cl_as_requests, cl_as_request_md, cl_as_request_str, cl_as_request_payl

USE yomhook, ONLY: lhook, dr_hook

USE parkind1, ONLY: jprb, jpim

USE IOS_print_mgr, ONLY:                                                     &
    IOS_print,                                                               &
    IOS_Verbosity,                                                           &
    IOS_message,                                                             &
    IOS_PrStatus_Oper,                                                       &
    IOS_PrStatus_Normal,                                                     &
    IOS_PrStatus_Diag,                                                       &
    IOS_PrStatus_Debug

USE fort2c_memcpy_interfaces, ONLY: um_memcpy64

USE get_wallclock_time_mod, ONLY: get_wallclock_time
 
IMPLICIT NONE

INTEGER, POINTER                :: ioGlobalNo(:,:)=>NULL()
INTEGER                         :: IOS_Client_Queue_Slot
INTEGER                         :: IOS_Client_Queue_Size=0
REAL                            :: IOS_cl_queue_stall=0.0
REAL                            :: IOS_lb_stall=0.0   ! stall time from 
                                                      ! getting load balance
#if defined(CPP_IOS_STATIC_DATA)
TYPE(IOS_async_object), TARGET  :: IOS_Dispatch_Queue(cpp_ios_concurrency)
#else
TYPE(IOS_async_object), POINTER :: IOS_Dispatch_Queue(:)
#endif
CHARACTER (LEN=IOS_err_strlen),                                              &
    PRIVATE                     :: IOS_clq_message

! params/vars  for dr_hook
INTEGER(KIND=jpim),                                                          &
    PARAMETER, PRIVATE          :: zhook_in  = 0
INTEGER(KIND=jpim),                                                          &
    PARAMETER, PRIVATE          :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='IOS_CLIENT_QUEUE'

CONTAINS


SUBROUTINE IOS_Client_Init()
USE mpl, ONLY: mpl_request_null
IMPLICIT NONE
INTEGER             :: i
REAL(KIND=jprb)     :: zhook_handle
CHARACTER (LEN=*),                                                         &
    PARAMETER :: RoutineName = 'IOS_CLIENT_INIT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

#if !defined(CPP_IOS_STATIC_DATA)
ALLOCATE(IOS_Dispatch_Queue(IOS_Concurrency))
! Never deallocated.
#endif
ALLOCATE(IOS_Sequence_id(IOS_Server_Groups,IOS_Tasks_per_server))
! Never deallocated.
ALLOCATE(ioGlobalNo     (IOS_Server_Groups,IOS_Tasks_per_server))
! Never deallocated.
IOS_Sequence_id(:,:)=0
IOS_Client_Queue_Slot=0
DO i=1,IOS_Concurrency
  IOS_Dispatch_Queue(i)%state       =IOS_queue_slot_unused
  IOS_Dispatch_Queue(i)%request(:)  =mpl_request_null
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE IOS_Client_Init


SUBROUTINE IOS_assign_server_for_unit(UNIT)
USE mpl, ONLY: mpl_integer
IMPLICIT NONE

INTEGER, INTENT(IN) :: UNIT
INTEGER, PARAMETER  :: loadBalalanceBcastRoot  = 0
INTEGER, PARAMETER  :: loadBalalanceBcastItems = 1
INTEGER             :: server
INTEGER             :: subtask
INTEGER             :: qHandle
INTEGER,SAVE        :: nextserver=0
INTEGER             :: targetObject
INTEGER             :: loads(loadBalanceDataSize)
INTEGER             :: minLoad
INTEGER             :: minSize
INTEGER             :: curLoad
INTEGER             :: curSize
INTEGER             :: minLoad_server
INTEGER             :: minSize_server
INTEGER             :: errorCode
INTEGER(IOS_tukind),                                                       &
    POINTER         :: recvBuffer(:)
LOGICAL             :: reallocation
REAL                :: t_start   ! Timer
REAL                :: t_end     ! Timer
TYPE(IOS_metadata_type),                                                   &
    POINTER         :: metadata

CHARACTER (LEN=*),                                                         &
    PARAMETER :: RoutineName = 'IOS_CLIENT_QUEUE:IOS_ASSIGN_SERVER_FOR_UNIT'

ErrorCode=-10
targetObject=-1
reallocation=.FALSE.

IF (L_IOS_Active()) THEN

  SELECT CASE(IOS_Unit_Alloc_Policy)
  CASE (IOS_Unit_Alloc_Static)

    ! assignment was made at init time,
    ! so this is a null op or an error

    IF (io_server_for_unit(UNIT) == IOS_No_Server )                        &
        CALL IOS_Ereport('IOS_Assign_server_for_unit', ErrorCode,          &
        'Unallocated unit was presented in a static scheme' )

  CASE (IOS_Unit_Alloc_AtFirstUse)

    ! assignment was not made at init time,
    ! but we may not reassign existing mappings

    IF (io_server_for_unit(UNIT) == IOS_No_Server ) THEN

      nextserver=nextserver+1
      IF (nextserver > IOS_Server_Groups) nextserver=1
      targetObject=io_servers(nextServer,1)

    END IF

  CASE (IOS_Unit_Alloc_Dynamic_Rotate)

    nextserver=nextserver+1
    IF (nextserver > IOS_Server_Groups) nextserver=1
    targetObject=io_servers(nextServer,1)

    IF (io_server_for_unit(UNIT) /= IOS_No_Server )                        &
        reallocation=.TRUE.

  CASE (IOS_Unit_Alloc_Dynamic_LB)

    IF ( io_server_for_unit(UNIT) == IOS_No_Server ) THEN

      nextserver=nextserver+1
      IF (nextserver > IOS_Server_Groups) nextserver=1
      targetObject=io_servers(nextServer,1)

    ELSE

      t_start= get_wallclock_time()
      IF (model_rank == 0 ) THEN
        minLoad=HUGE(minload)
        minSize=HUGE(minSize)
        minLoad_server=-1
        minSize_server=-1

        DO server=1,IOS_Server_Groups
          qHandle=IOS_init_md(-1*io_servers(server,1),                     &
              IOS_No_Location, IOS_Action_LoadStatus,                      &
              datasize=loadBalanceDataSize)
          CALL IOS_Send(qHandle)
          recvBuffer => IOS_attach_recvBuffer(qHandle)

          loads(:)=-99
          CALL um_memcpy64(loads,recvBuffer,loadBalanceDataSize)

          IF (io_server_for_unit(UNIT)==io_servers(server,1)) THEN
            IF (IOS_Verbosity>=IOS_prstatus_Oper) THEN
!$OMP CRITICAL(internal_write)
              WRITE(IOS_message,'(A,F10.3,A,I4,A,I4)')                     &
                  'Info: Current queue size=',                             &
                  loads(loadBalanceQueueData)/1024.0/1024.0,               &
                  'MB items=',                                             &
                  loads(loadBalanceQueueLen),                              &
                  ' on server ',                                           &
                  io_servers(server,1)
!$OMP END CRITICAL(internal_write)
              CALL IOS_print(IOS_message,src='ios_client_queue')
            END IF
          END IF

          IF (loads(loadBalanceQueueLen) < minLoad) THEN
            minload        = loads(loadBalanceQueueLen)
            minLoad_server = server
          END IF
          IF (loads(loadBalanceQueueData) < minSize) THEN
            minSize        = loads(loadBalanceQueueData)
            minSize_server = server
          END IF
        END DO

        IF (IOS_Verbosity>=IOS_prstatus_Oper) THEN
!$OMP CRITICAL(internal_write)
          WRITE(IOS_message,'(A,F8.2,A,I4,A,I4,A,I4,A)')                   &
              'Info: Min queue size=',                                     &
              minSize/1024.0/1024.0,                                       &
              ' MB (',                                                     &
              io_servers(minSize_Server,1),                                &
              ') Min queue load=',                                         &
              minLoad,                                                     &
              ' items (',                                                  &
              io_servers(minLoad_Server,1),')'
!$OMP END CRITICAL(internal_write)
          CALL IOS_print(IOS_message,src='ios_client_queue')
        END IF

      END IF

      IF (Model_Procs > 1) THEN
        CALL MPL_Bcast(minSize_Server,                                     &
            loadBalalanceBcastItems, mpl_integer,                          &
            loadBalalanceBcastRoot, Model_Comm,                            &
            errorCode)
      END IF
      targetObject=io_servers(minSize_server,1)

      IF (io_server_for_unit(UNIT) == targetObject) THEN

        ! Do not waste the servers energy, if we are remapping
        ! to the same place

        IF (IOS_Verbosity>=IOS_prstatus_Oper) THEN
!$OMP CRITICAL(internal_write)
          WRITE(IOS_message,'(A,I3,A,I4)')                                 &
              'Info: Unit ',UNIT,                                          &
              ' is not reassigned, and remains on ',                       &
              io_server_for_unit(UNIT)
!$OMP END CRITICAL(internal_write)
          CALL IOS_print(IOS_message,src='ios_client_queue')
        END IF
        reallocation=.FALSE.
        targetObject=-1
      ELSE
        reallocation=.TRUE.
      END IF

      t_end = get_wallclock_time()
      IOS_lb_stall=IOS_lb_stall+(t_end-t_start)
      IF (IOS_Verbosity>=IOS_prstatus_oper) THEN
        WRITE(IOS_message,'(A,F10.3,A,F8.3)')                              &
            'Info: Stall getting load balance at ',                        &
            t_start-IOS_Start_time,' of ',t_end-t_start
        CALL IOS_print(IOS_message,src='ios_client_queue')
      END IF

    END IF

  CASE DEFAULT

    ErrorCode=-10
    CALL IOS_Ereport('IOS_Assign_server_for_unit', ErrorCode,              &
        'Unknown server allocation policy' )

  END SELECT


  IF (targetObject >= 0) THEN
    io_server_for_unit_lookup(UNIT)=targetObject

    IF (IOS_Verbosity>=IOS_prstatus_Oper) THEN
      IF (reallocation) THEN
!$OMP CRITICAL(internal_write)
        WRITE(IOS_message,'(A,I3,A,I4)')                                   &
            'Info: Unit ',UNIT,                                            &
            ' is reassigned to IO server ',                                &
            targetObject
!$OMP END CRITICAL(internal_write)
        CALL IOS_print(IOS_message,src='ios_client_queue')
      ELSE
!$OMP CRITICAL(internal_write)
        WRITE(IOS_message,'(A,I3,A,I4)')                                   &
            'Info: Unit ',UNIT,                                            &
            ' is assigned to IO server ',                                  &
            targetObject
!$OMP END CRITICAL(internal_write)
        CALL IOS_print(IOS_message,src='ios_client_queue')
      END IF
    END IF

    IF (model_rank == 0 ) THEN

      DO server=1,IOS_Server_Groups
        DO subtask=1,IOS_Metadata_Receivers()
          qHandle=IOS_init_md(-1*io_servers(server,subtask),               &
              IOS_No_Location, IOS_Action_Assign_Unit)
          metadata => IOS_Attach_Metadata(qHandle)
          metadata % subtype = UNIT
          metadata % address = targetObject
          CALL IOS_Send(qHandle)
        END DO
      END DO
    END IF
  END IF

END IF

END SUBROUTINE IOS_assign_server_for_unit


FUNCTION IOS_getDestPe(handle) RESULT(pe)
IMPLICIT NONE

INTEGER, INTENT(IN) :: handle
INTEGER             :: pe
pe=IOS_Dispatch_Queue(handle)%pe
END FUNCTION IOS_getDestPe


INTEGER FUNCTION IOS_init_md                                                 &
    (targetObject, location, ACTION, dataSize, bcast,                        &
    targetRank) RESULT(handle)
IMPLICIT NONE

INTEGER, INTENT(IN)                  :: targetObject
! positive values are a unit
! negative values (or 0) are a processor.
INTEGER, INTENT(IN)                  :: ACTION
INTEGER, INTENT(IN)                  :: location
INTEGER, OPTIONAL                    :: dataSize
INTEGER, OPTIONAL                    :: bcast
INTEGER, OPTIONAL                    :: targetRank
! targetRank is used to specify a processor, when the targetObject is a
! unit, and the command must reach all team members.
INTEGER                              :: nextSlot
LOGICAL                              :: ok
LOGICAL                              :: flag
LOGICAL                              :: l_bcast
LOGICAL                              :: waitforbuffer
INTEGER                              :: errorFlag
REAL(KIND=jprb)                      :: zhook_handle
CHARACTER (LEN=*),                                                         &
    PARAMETER :: RoutineName = 'IOS_INIT_MD'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

ok=assert_client()

IF (IOS_Verbosity>=IOS_prstatus_debug) THEN
  WRITE(IOS_message,'(A,A)')                                               &
      'Info: Getting protocol queue slot for ',                            &
      TRIM(IOS_ActionName(ACTION))
  CALL IOS_print(IOS_message,src='ios_client_queue')
END IF
NextSlot=IOS_Client_Queue_Slot+1
IF (NextSlot>IOS_Concurrency)NextSlot=1

IF (IOS_Dispatch_Queue(NextSlot)%state==                                   &
    IOS_queue_slot_initialized) THEN

  WRITE(IOS_clq_message,'(A)')                                             &
      'The next queue slot is already initialzed'
  errorFlag=99
  CALL IOS_Ereport( RoutineName, errorFlag , IOS_clq_message )

ELSE IF (IOS_Dispatch_Queue(NextSlot)%state==                               &
    ios_queue_slot_dispatched) THEN
  waitforbuffer = .TRUE.
  flag=IOS_QueryBuffer(NextSlot,'IOS_init_md: Wating for buffer ',waitforbuffer)

ELSE IF (IOS_Dispatch_Queue(NextSlot)%state/=IOS_queue_slot_unused) THEN

!$OMP CRITICAL(internal_write)
  WRITE(IOS_clq_message,'(A,I6)')'IOS_Client_Queue: UNKNOWN STATE=',       &
      IOS_Dispatch_Queue(NextSlot)%state
!$OMP END CRITICAL(internal_write)
  errorFlag=99
  CALL IOS_Ereport( RoutineName, errorFlag , IOS_clq_message )

END IF

! We can now safely reuse our target slot
IOS_Client_Queue_Slot=nextSlot

! The handle for this request is the slot id
handle=IOS_Client_Queue_Slot

! Set defaults from arguments
IOS_Dispatch_Queue(handle)%md%ACTION  = ACTION
IOS_Dispatch_Queue(handle)%md%UNIT    = targetObject
IOS_Dispatch_Queue(handle)%md%address = location

! Blank other fields !
IOS_Dispatch_Queue(handle)%md%name_length = 0
IOS_Dispatch_Queue(handle)%md%subtype     = 0
IOS_Dispatch_Queue(handle)%md%delete      = 0
IOS_Dispatch_Queue(handle)%md%data_size   = 0
IOS_Dispatch_Queue(handle)%md%client      = model_rank
IOS_Dispatch_Queue(handle)%md%handle      = 0
IOS_Dispatch_Queue(handle)%md%Originating_Slot = -99

NULLIFY(IOS_Dispatch_Queue(handle)%payload)
WRITE(IOS_Dispatch_Queue(handle)%md%string,'(A)')'NotSet'

! Figure out the destination pe
IF (targetObject <= 0) THEN
  IOS_Dispatch_Queue(handle)%pe=-1*targetObject
ELSE IF (targetObject == file_op_pseudo_unit) THEN
!$OMP CRITICAL(internal_write)
  WRITE(IOS_clq_message,'(A,I5,A)') 'A target object of ',                 &
      targetObject,' is not allowed'
!$OMP END CRITICAL(internal_write)
  errorFlag=61
  CALL IOS_ereport( RoutineName, errorFlag, IOS_clq_message )
ELSE
  IOS_Dispatch_Queue(handle)%pe=io_server_for_unit(targetObject)
  IF (IOS_Dispatch_Queue(handle)%pe == IOS_No_Server) THEN
!$OMP CRITICAL(internal_write)
    WRITE(IOS_clq_message,'(A,I3,A)') 'No server for unit ',               &
        targetObject,' has been assigned'
!$OMP END CRITICAL(internal_write)
    errorFlag=61
    CALL IOS_Ereport( RoutineName, errorFlag, IOS_clq_message )
  END IF

  IF (PRESENT(targetRank)) THEN
    IOS_Dispatch_Queue(handle)%pe=                                         &
        io_servers(ioServerNo(IOS_Dispatch_Queue(handle)%pe),targetRank)
  END IF
END IF

! Do we need to set the broadcast flag?
l_bcast=.FALSE.
IF (PRESENT(bcast)) THEN
  IF (bcast==IOS_Read_Broadcast) THEN
    l_bcast=.TRUE.
    IOS_Dispatch_Queue(handle)%md%subtype=bcast
  END IF
END IF

! Do we need a data buffer?
IF (PRESENT(dataSize)) THEN
  IOS_Dispatch_Queue(handle)%md%data_size=dataSize

  IF (.NOT. l_bcast) THEN

    SELECT CASE(ACTION)

    CASE (IOS_Action_Write32,                                              &
        IOS_Action_Read32_Integer,                                         &
        IOS_Action_Read32                                                  &
        )
      CALL IOS_consume_client_mem(handle,dataSize*IOS_BytesPerWord32)
      ALLOCATE (IOS_Dispatch_Queue(handle)%payload                         &
          (dataSize*ios_tuperword32))

    CASE (                                                                 &
        IOS_Action_StashInitModel,                                         &
        IOS_Action_StashWritePPData,                                       &
        IOS_Action_StashWriteDumpData,                                     &
        IOS_Action_DumpInitModel,                                          &
        IOS_Action_Enquire,                                                &
        IOS_Action_LoadStatus,                                             &
        IOS_Action_Write64,                                                &
        IOS_Action_Read64,                                                 &
        IOS_Action_Read64_Integer,                                         &
        IOS_Action_Sync_Barrier,                                           &
        IOS_Action_StashInitPPLookup,                                      &
        IOS_Action_Write_PP_Prepacked,                                     &
        IOS_Action_StashWritePPLookup,                                     &
        IOS_Action_StashInitHeader,                                        &
        IOS_Action_StashSetHeader,                                         &
        IOS_Action_MergePPLookup)
      CALL IOS_consume_client_mem(handle,dataSize*IOS_BytesPerWord64)
      ALLOCATE (IOS_Dispatch_Queue(handle)%payload                         &
          (dataSize*IOS_tuperword64))
    CASE DEFAULT
!$OMP CRITICAL(internal_write)
      WRITE(IOS_clq_message,'(A,I3,A)') 'Action ',                         &
          ACTION,' not recognised!'
!$OMP END CRITICAL(internal_write)
      errorFlag=60
      CALL IOS_Ereport( RoutineName, errorFlag, IOS_clq_message )

    END SELECT
  END IF
END IF

IF (IOS_Verbosity>=IOS_prstatus_debug) THEN
!$OMP CRITICAL(internal_write)
  WRITE(IOS_message,'(A,I4)')'Info: New queue slot @ ',handle
!$OMP END CRITICAL(internal_write)
  CALL IOS_print(IOS_message,src='ios_client_queue')
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END FUNCTION IOS_init_md


FUNCTION IOS_Attach_Metadata(handle) RESULT(md)
IMPLICIT NONE

INTEGER, INTENT(IN)              :: handle
TYPE(IOS_Metadata_Type), POINTER :: md
REAL(KIND=jprb)                  :: zhook_handle
CHARACTER (LEN=*),                                                         &
    PARAMETER :: RoutineName = 'IOS_ATTACH_METADATA'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

md => IOS_Dispatch_Queue(handle)%md

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END FUNCTION IOS_Attach_Metadata


FUNCTION IOS_Attach_SendBuffer(handle) RESULT(sb)
IMPLICIT NONE

INTEGER, INTENT(IN)               :: handle
INTEGER(KIND=IOS_tukind), POINTER :: sb(:)
REAL(KIND=jprb)                   :: zhook_handle
INTEGER                           :: errorFlag
CHARACTER (LEN=*),                                                         &
    PARAMETER :: RoutineName = 'IOS_ATTACH_SENDBUFFER'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

errorFlag=99

IF (.NOT. ASSOCIATED(IOS_Dispatch_Queue(handle)%payload)) THEN
  WRITE(IOS_clq_message,'(A)')                                             &
      'Attach_SendBuffer Failed, no buffer allocated'
  CALL IOS_Ereport( RoutineName, errorFlag , IOS_clq_message )
END IF

sb => IOS_Dispatch_Queue(handle)%payload(:)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END FUNCTION IOS_Attach_SendBuffer


FUNCTION IOS_Attach_RecvBuffer(handle) RESULT(rb)
IMPLICIT NONE

INTEGER, INTENT(IN)               :: handle
INTEGER(KIND=IOS_tukind), POINTER :: rb(:)
REAL(KIND=jprb)                   :: zhook_handle
INTEGER                           :: errorFlag
CHARACTER (LEN=*),                                                         &
    PARAMETER :: RoutineName = 'IOS_ATTACH_RECVBUFFER'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF (.NOT. ASSOCIATED(IOS_Dispatch_Queue(handle)%payload)) THEN
  WRITE(IOS_clq_message,'(A)')                                             &
      'Attach_RecvBuffer Failed, no buffer allocated'
  errorFlag=99
  CALL IOS_Ereport( RoutineName,errorFlag , IOS_clq_message )
END IF

rb => IOS_Dispatch_Queue(handle)%payload

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END FUNCTION IOS_Attach_RecvBuffer


SUBROUTINE IOS_Send(handle,hasString)
USE mpl, ONLY:                                                            &
    mpl_status_size,                                                       &
    mpl_integer,                                                           &
    mpl_character

IMPLICIT NONE

INTEGER, INTENT(IN)   :: handle
LOGICAL, OPTIONAL     :: hasString
INTEGER               :: errorFlag
INTEGER               :: ioServer
INTEGER               :: ioServerR
INTEGER               :: tagOffset
INTEGER               :: mdtag
INTEGER               :: pltag
INTEGER               :: strtag
INTEGER               :: STATUS(mpl_status_size)
REAL(KIND=jprb)       :: zhook_handle

CHARACTER (LEN=*),                                                         &
    PARAMETER :: RoutineName = 'IOS_SEND'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

ioServer  = ioServerNo   (IOS_Dispatch_Queue(handle)%pe)
ioServerR = ioServerRank (IOS_Dispatch_Queue(handle)%pe)

IF ( IOS_Dispatch_Queue(handle)%md%ACTION ==                               &
    IOS_Action_LoadStatus) THEN
  mdtag  = IOS_Request_Tag_Express
  pltag  = IOS_Request_Tag_Express
  strtag = IOS_Request_Tag_Express
ELSE
  IOS_Sequence_ID(ioServer,ioServerR+1)=                                   &
      IOS_Sequence_ID(ioServer,ioServerR+1)+1
  tagOffset=MOD(IOS_Sequence_ID(ioServer,ioServerR+1),                     &
      IOS_Request_Tag_Gap)

  mdtag  = IOS_Request_Tag_Base+tagOffset
  pltag  = IOS_Request_Tag_Payload+tagOffset
  strtag = IOS_Request_Tag_Str+tagOffset
END IF

! Leave in handy for debugging.
IF (.FALSE.) THEN
  CALL IOS_print('IOS SEND: -------------------------------------------',  &
      src='ios_client_queue')
!$OMP CRITICAL(internal_write)
  WRITE(IOS_message,'(A,I0,A,I0)')                                         &
      'IOS SEND: targ pe=',IOS_Dispatch_Queue(handle)%pe,                  &
      '  seq=',IOS_Sequence_ID(ioServer,ioServerR+1)
!$OMP END CRITICAL(internal_write)
  CALL IOS_print(IOS_message,src='ios_client_queue')
!$OMP CRITICAL(internal_write)
  WRITE(IOS_message,'(A,I0)')'IOS SEND: server=',ioServer
!$OMP END CRITICAL(internal_write)
  CALL IOS_print(IOS_message,src='ios_client_queue')
!$OMP CRITICAL(internal_write)
  WRITE(IOS_message,'(A,I0)')'IOS SEND: rank=',ioServerR
!$OMP END CRITICAL(internal_write)
  CALL IOS_print(IOS_message,src='ios_client_queue')
!$OMP CRITICAL(internal_write)
  WRITE(IOS_message,'(A,I0)')'IOS SEND: handle=',handle
!$OMP END CRITICAL(internal_write)
  CALL IOS_print(IOS_message,src='ios_client_queue')
!$OMP CRITICAL(internal_write)
  WRITE(IOS_message,'(A,4I0)')'IOS SEND: tagOffset=',                      &
      tagOffset,mdtag,pltag,strtag
!$OMP END CRITICAL(internal_write)
  CALL IOS_print(IOS_message,src='ios_client_queue')
!$OMP CRITICAL(internal_write)
  WRITE(IOS_message,'(A,I0)')'IOS SEND: action=',                          &
      IOS_Dispatch_Queue(handle)%md%ACTION
!$OMP END CRITICAL(internal_write)
  CALL IOS_print(IOS_message,src='ios_client_queue')
  WRITE(IOS_message,'(A,A)')'IOS SEND: action=',                           &
      TRIM(IOS_ActionName(IOS_Dispatch_Queue(handle)%md%ACTION))
  CALL IOS_print(IOS_message,src='ios_client_queue')
!$OMP CRITICAL(internal_write)
  WRITE(IOS_message,'(A,I0)')'IOS SEND: subtype=',                         &
      IOS_Dispatch_Queue(handle)%md%subtype
!$OMP END CRITICAL(internal_write)
  CALL IOS_print(IOS_message,src='ios_client_queue')
  CALL IOS_print('IOS SEND: -------------------------------------------',  &
      src='ios_client_queue')
END IF
IF (ioServerR == -1) THEN
!$OMP CRITICAL(internal_write)
  WRITE(IOS_clq_message,'(A,I5,A)')                                        &
      'IOS_Send: Not permitted to address pe ',                            &
      IOS_Dispatch_Queue(handle)%pe,                                       &
      ' as it is not reporting an IO server rank'
!$OMP END CRITICAL(internal_write)
  errorFlag=99
  CALL IOS_Ereport( RoutineName, errorFlag , IOS_clq_message )
END IF

IF (IOS_RelayToSlaves .AND. ioServerR /= 0) THEN
!$OMP CRITICAL(internal_write)
  WRITE(IOS_clq_message,'(A,I5,A,I5)')                                     &
      'IOS_Send: Not permitted to address pe ',                            &
      IOS_Dispatch_Queue(handle)%pe,                                       &
      ' with IOS_RelayToSlaves rank= ',                                    &
      ioServerR
!$OMP END CRITICAL(internal_write)
  errorFlag=99
  CALL IOS_Ereport( RoutineName, errorFlag , IOS_clq_message )
END IF

CALL MPL_iSend(                                                            &
    IOS_Dispatch_Queue(handle)%md,                                         &
    IOS_md_len,                                                            &
    mpl_integer,                                                           &
    IOS_Dispatch_Queue(handle)%pe,                                         &
    mdtag,                                                                 &
    Global_Comm,                                                           &
    IOS_Dispatch_Queue(handle)%request(cl_as_request_md),                  &
    errorFlag)

IF (PRESENT(hasString)) THEN
  IF (hasString) THEN
    CALL MPL_iSend(                                                        &
        IOS_Dispatch_Queue(handle)%md%string,                              &
        IOS_String_Max,                                                    &
        mpl_character,                                                     &
        IOS_Dispatch_Queue(handle)%pe,                                     &
        strtag,                                                            &
        Global_Comm,                                                       &
        IOS_Dispatch_Queue(handle)%request(cl_as_request_str),             &
        errorFlag)
  END IF
END IF

IF (ASSOCIATED(IOS_Dispatch_Queue(handle)%payload)) THEN
  SELECT CASE(IOS_Dispatch_Queue(handle)%md%ACTION)
  CASE (                                                                   &
      IOS_Action_Enquire,                                                  &
      IOS_Action_LoadStatus,                                               &
      IOS_Action_Read32_Integer,                                           &
      IOS_Action_Read32,                                                   &
      IOS_Action_Read64_Integer,                                           &
      IOS_Action_Read64                                                    &
      )
    IF (IOS_Dispatch_Queue(handle) % md%subtype /= IOS_Read_Broadcast) THEN
      CALL MPL_Recv(                                                       &
          IOS_Dispatch_Queue(handle)%payload,                              &
          SIZE(IOS_Dispatch_Queue(handle)%payload),                        &
          IOS_tutype,                                                      &
          IOS_Dispatch_Queue(handle)%pe,                                   &
          pltag,                                                           &
          Global_Comm,                                                     &
          STATUS,                                                          &
          errorFlag)
    END IF
  CASE DEFAULT
    CALL MPL_iSend(                                                        &
        IOS_Dispatch_Queue(handle)%payload,                                &
        SIZE(IOS_Dispatch_Queue(handle)%payload),                          &
        IOS_tutype,                                                        &
        IOS_Dispatch_Queue(handle)%pe,                                     &
        pltag,                                                             &
        Global_Comm,                                                       &
        IOS_Dispatch_Queue(handle)%request(cl_as_request_payl),            &
        errorFlag)
  END SELECT
END IF
IOS_Dispatch_Queue(handle)%state=IOS_queue_slot_dispatched

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE IOS_Send

!---------------------------------------------------------------------
! Finalisation - terminate remote IO Servers
!---------------------------------------------------------------------
SUBROUTINE IOS_Finalise()
IMPLICIT NONE

INTEGER                            :: server
INTEGER                            :: subtask
TYPE(IOS_metadata_type),POINTER    :: metadata
INTEGER                            :: qHandle
CHARACTER (LEN=*),                                                         &
    PARAMETER :: RoutineName = 'IOS_CLIENT_QUEUE:IOS_FINALISE'

IF (model_rank == 0 ) THEN

  DO server=1,IOS_Server_Groups
    DO subtask=1,IOS_Metadata_Receivers()
      qHandle=IOS_init_md(-1*io_servers(server,subtask),                   &
          -1,IOS_Action_Finish)
      CALL IOS_Send(qHandle)
    END DO
  END DO
END IF

RETURN
END SUBROUTINE IOS_Finalise

!---------------------------------------------------------------------
! Shutdown, Close down servers, and tidy up pending messages
!---------------------------------------------------------------------
SUBROUTINE IOS_Shutdown()

IMPLICIT NONE

INTEGER                              :: LastUsedSlot
INTEGER                              :: NextSlot
INTEGER                              :: errorFlag
LOGICAL                              :: flag
LOGICAL                              :: waitforbuffer
REAL                                 :: t1
REAL                                 :: t2
REAL(KIND=jprb)                      :: zhook_handle
CHARACTER (LEN=*),                                                         &
    PARAMETER :: RoutineName = 'IOS_SHUTDOWN'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF (model_rank == 0) THEN

  CALL IOS_Finalise()

  ! Wait completion of queue
  lastUsedSlot=IOS_Client_Queue_Slot
  NextSlot    =-1

  IF (NextSlot>IOS_Concurrency)NextSlot=1

  t1=get_wallclock_time()

  DO WHILE (NextSlot /= LastUsedSlot)
    IF (NextSlot == -1) THEN
      NextSlot=lastUsedSlot+1
    ELSE
      NextSlot=NextSlot+1
    END IF

    IF (NextSlot>IOS_Concurrency)NextSlot=1

    IF (IOS_Verbosity>=IOS_prstatus_diag) THEN
!$OMP CRITICAL(internal_write)
      WRITE(ios_message,'(A,I3)')                                          &
          'Shutdown: Waiting completion of Q slot: ',NextSlot
!$OMP END CRITICAL(internal_write)
    END IF
    SELECT CASE(IOS_Dispatch_Queue(NextSlot)%state)

    CASE (IOS_queue_slot_initialized)
      IF (IOS_Verbosity>=IOS_prstatus_diag) THEN
        CALL IOS_print( 'Shutdown:  --> IOS_QUEUE_SLOT_INITIALIZED',       &
            src='ios_client_queue')
      END IF
    CASE (IOS_queue_slot_partfilled)
      IF (IOS_Verbosity>=IOS_prstatus_diag) THEN
        CALL IOS_print( 'Shutdown:  --> IOS_QUEUE_SLOT_PARTFILLED',        &
            src='ios_client_queue')
      END IF
    CASE (IOS_queue_slot_dispatched)
      IF (IOS_Verbosity>=IOS_prstatus_diag) THEN
        CALL IOS_print( 'Shutdown:  --> IOS_QUEUE_SLOT_DISPATCHED',        &
            src='ios_client_queue')
      END IF
      waitforbuffer= .TRUE.
      flag=IOS_QueryBuffer                                                 &
          (NextSlot,'Shutdown:  --> Waiting for buffer ',waitforbuffer)

    CASE (IOS_queue_slot_unused)
      IF (IOS_Verbosity>=IOS_prstatus_diag) THEN
        CALL IOS_print('Shutdown:  --> IOS_QUEUE_SLOT_UNUSED',             &
            src='ios_client_queue')
      END IF
    CASE DEFAULT
!$OMP CRITICAL(internal_write)
      WRITE(IOS_clq_message,'(A,I0)')'IOS_Client_Queue: UNKNOWN STATE=',   &
          IOS_Dispatch_Queue(NextSlot)%state
!$OMP END CRITICAL(internal_write)
      errorFlag=99
      CALL IOS_Ereport( RoutineName, errorFlag , IOS_clq_message )

    END SELECT

  END DO
  t2=get_wallclock_time()

  IF (IOS_Verbosity>=IOS_PrStatus_Normal) THEN
    WRITE(IOS_message,'(A,F8.3)')                                          &
        'Info: total stall time in protocol queue = ',                     &
        IOS_cl_queue_stall
    CALL IOS_print(IOS_message,src='ios_client_queue')

    ! Only write out Stall time due to getting load balance info if this 
    ! policy is selected
    IF (IOS_Unit_Alloc_Policy == IOS_Unit_Alloc_Dynamic_LB) THEN
      WRITE(IOS_message,'(A,F8.3)')                                        &
          'Info: total stall time obtaining load balance = ',              &
          IOS_lb_stall
      CALL IOS_print(IOS_message,src='ios_client_queue')
    END IF

    WRITE(IOS_message,'(A,F8.3)')                                          &
        'Info: total time waiting for shutdown    = ',                     &
        (t2-t1)
    CALL IOS_print(IOS_message,src='ios_client_queue')
  END IF

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE IOS_Shutdown

! Enquire about the status of a protocol buffer.
!  If wait is true, the we will spin until the operation completes.
!  We return truew if the op is completed, and false otherwise.
!  If an op is completed, we will reset the state accordingly.
LOGICAL FUNCTION IOS_QueryBuffer(slot,messg,waitForBuffer)
USE mpl, ONLY:                                                            &
    mpl_request_null,                                                      &
    mpl_status_size
IMPLICIT NONE
INTEGER, INTENT(IN)    :: slot
LOGICAL, INTENT(IN)    :: waitForBuffer
CHARACTER(LEN=*)       :: messg
INTEGER                :: errorFlag
INTEGER                :: mpi_statuses(mpl_status_size,3)
REAL                   :: t_start
REAL                   :: t_end
LOGICAL                :: completionFlag
CHARACTER (LEN=*),                                                         &
    PARAMETER :: RoutineName = 'IOS_QUERYBUFFER'
REAL(KIND=jprb):: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Assume incomplete
IOS_QueryBuffer=.FALSE.

IF (IOS_Dispatch_Queue(Slot)%state/=ios_queue_slot_dispatched) THEN
  ! We don't need to wait for it!
  IOS_QueryBuffer=.TRUE.
ELSE
  completionFlag=.FALSE.
  CALL mpl_testall                                                         &
      (cl_as_requests,IOS_Dispatch_Queue(Slot)%request,completionFlag,     &
      mpi_statuses,errorFlag)

  IF (waitForBuffer .AND. .NOT. completionFlag) THEN

    IF (IOS_Verbosity>=IOS_prstatus_oper) THEN
!$OMP CRITICAL(internal_write)
      WRITE(IOS_message,'(A,I4)') messg, slot
!$OMP END CRITICAL(internal_write)
      CALL IOS_print(IOS_message,src='ios_client_queue')
    END IF
    t_start=get_wallclock_time()
    DO WHILE (.NOT. completionFlag)
      CALL MPL_Testall(cl_as_requests,IOS_Dispatch_Queue(Slot)%request,    &
          completionFlag,mpi_statuses,errorFlag)
      t_end=get_wallclock_time()
      IF (t_end-t_start>IOS_Timeout) THEN
        CALL IOS_print_protocol_slot(slot)
!$OMP CRITICAL(internal_write)
        WRITE(IOS_clq_message,'(A,I5,A,I5,A)')                             &
            'Time out waiting for protocol buffer ',slot,                  &
            ' addressing pe ',IOS_Dispatch_Queue(Slot)%pe
!$OMP END CRITICAL(internal_write)
        errorFlag=99
        CALL IOS_Ereport( RoutineName, errorFlag , IOS_clq_message )
      END IF
    END DO ! waiting for completion

    IOS_cl_queue_stall=IOS_cl_queue_stall+(t_end-t_start)
    IF (IOS_Verbosity>=IOS_prstatus_oper) THEN
      WRITE(IOS_message,'(A,F10.3,A,F8.3)')                                &
          'Info: Stall getting protocol queue slot at ',                   &
          t_start-IOS_Start_time,' of ',t_end-t_start
      CALL IOS_print(IOS_message,src='ios_client_queue')
    END IF

  END IF ! wait for completion

  IF (completionFlag) THEN

    IOS_QueryBuffer=.TRUE.

    IF (IOS_Verbosity>=IOS_prstatus_debug) THEN
!$OMP CRITICAL(internal_write)
      WRITE(IOS_message,'(A,I5)') 'Info: Completed slot: ',slot
!$OMP END CRITICAL(internal_write)
      CALL IOS_print(IOS_message,src='ios_client_queue')
    END IF

    IF (ASSOCIATED(IOS_Dispatch_Queue(Slot)%payload)) THEN
      IOS_Client_Queue_Size=IOS_Client_Queue_Size-                         &
          SIZE(IOS_Dispatch_Queue(Slot)%payload)*IOS_BytesPerTU

      IF (IOS_Verbosity>=IOS_prstatus_debug) THEN
!$OMP CRITICAL(internal_write)
        WRITE(IOS_message,'(A,I20,A,F8.2,A)')                              &
            'Info: Freeing ',                                              &
            SIZE(IOS_Dispatch_Queue(Slot)%payload)*IOS_BytesPerTU,         &
            ' bytes. Consumed payload memory now ',                        &
            IOS_Client_Queue_Size/(1024.0*1024.0),' MB'
!$OMP END CRITICAL(internal_write)
        CALL IOS_print(IOS_message,src='ios_client_queue')
      END IF
      DEALLOCATE(IOS_Dispatch_Queue(Slot)%payload)
      IF (IOS_Client_Queue_Size<0) THEN
        WRITE(IOS_clq_message,'(A)')                                       &
            'IOS Client thinks memory use is negative'
        errorFlag=99
        CALL IOS_Ereport( RoutineName, errorFlag , IOS_clq_message )
      END IF
    END IF

    IOS_Dispatch_Queue(Slot)%request(:)  =mpl_request_null
    IOS_Dispatch_Queue(Slot)%state       =ios_queue_slot_unused

  ELSE
    IF (IOS_Verbosity>=IOS_prstatus_debug) THEN
!$OMP CRITICAL(internal_write)
      WRITE(IOS_message,'(A,I20,A)')                                       &
          'Info: Slot ',slot,' still busy'
!$OMP END CRITICAL(internal_write)
      CALL IOS_print(IOS_message,src='ios_client_queue')
    END IF
  END IF ! actions on completion

END IF ! whether it was initially dispatched

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END FUNCTION IOS_QueryBuffer

! Wait for sufficient memory (mem bytes) to use for slot
SUBROUTINE IOS_consume_client_mem(slot,mem)
IMPLICIT NONE
INTEGER, INTENT(IN) :: slot
INTEGER, INTENT(IN) :: mem
INTEGER             :: i
INTEGER             :: errorFlag
LOGICAL             :: flag
LOGICAL             :: waitforbuffer
REAL                :: t_start
REAL                :: t_end
CHARACTER (LEN=*),                                                         &
    PARAMETER :: RoutineName = 'IOS_CONSUME_CLIENT_MEM'
REAL(KIND=jprb):: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF (mem > IOS_concurrency_max_mem) THEN
!$OMP CRITICAL(internal_write)
  WRITE(IOS_clq_message,'(A,I12,A,I12,A)')                                 &
      'IOS Client memory cap (',IOS_concurrency_max_mem,                   &
      ') is too small for payload (',mem,')'
!$OMP END CRITICAL(internal_write)
  errorFlag=99
  CALL IOS_Ereport( RoutineName, errorFlag , IOS_clq_message )
END IF


i=slot+1
IF (i>IOS_Concurrency)i=1

t_start=get_wallclock_time()
IF (IOS_Client_Queue_Size + mem > IOS_concurrency_max_mem) THEN
  IF (IOS_Verbosity>=IOS_prstatus_oper) THEN
!$OMP CRITICAL(internal_write)
    WRITE(IOS_message,'(A,I20,A)')'Info: Waiting for ',mem,                &
        ' bytes of protocol space'
!$OMP END CRITICAL(internal_write)
    CALL IOS_print(IOS_message,src='ios_client_queue')
  END IF
  DO WHILE (IOS_Client_Queue_Size + mem > IOS_concurrency_max_mem)
    IF (ASSOCIATED(IOS_Dispatch_Queue(i)%payload) .AND. i/=slot) THEN
      waitforbuffer = .FALSE.
      flag=IOS_QueryBuffer                                                 &
          (i,'IOS_consume_client_mem: query buffer ',waitforbuffer)
    END IF
    i=i+1
    IF (i>IOS_Concurrency)i=1

    t_end=get_wallclock_time()
    IF (t_end-t_start>IOS_Timeout) THEN
      DO i=1,IOS_Concurrency
        CALL IOS_print_protocol_slot(i)
      END DO
      WRITE(IOS_clq_message,'(A)')                                         &
          'Time out waiting for memory to become available '
      errorFlag=99
      CALL IOS_Ereport( RoutineName, errorFlag , IOS_clq_message )
    END IF


  END DO
  t_end=get_wallclock_time()

  IF (IOS_Verbosity>=IOS_prstatus_oper) THEN
    WRITE(IOS_message,'(A,F10.3,A,F10.3)')                                  &
        'Info: Stall waiting for protocol queue memory at ',               &
        t_start-IOS_Start_time,' of ',t_end-t_start
    CALL IOS_print(IOS_message,src='ios_client_queue')
  END IF

END IF

IOS_Client_Queue_Size = IOS_Client_Queue_Size + mem
IF (IOS_Verbosity>=IOS_prstatus_debug) THEN
  WRITE(IOS_message,'(A,F8.2,A,F8.2,A)')                                   &
      'Info: Consumed payload memory now ',                                &
      IOS_Client_Queue_Size/(1024.0*1024.0),                               &
      ' MB, Maximum is ',                                                  &
      IOS_Concurrency_max_mem/(1024.0*1024.0),' MB'
  CALL IOS_print(IOS_message,src='ios_client_queue')
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE IOS_consume_client_mem

! Diagnostic reporting
SUBROUTINE IOS_print_protocol_slot(slot)
USE mpl, ONLY: mpl_request_null
IMPLICIT NONE
INTEGER, INTENT(IN) :: slot

!$OMP CRITICAL(internal_write)
WRITE(IOS_message,'(A,I5)')  'State of slot ',slot
!$OMP END CRITICAL(internal_write)
CALL IOS_print(IOS_message,src='ios_client_queue')
!$OMP CRITICAL(internal_write)
WRITE(IOS_message,'(A,3I20)')'Requests: ',                                 &
    IOS_Dispatch_Queue(slot)%request(:)
!$OMP END CRITICAL(internal_write)
CALL IOS_print(IOS_message,src='ios_client_queue')
!$OMP CRITICAL(internal_write)
WRITE(IOS_message,'(A,I20)') '   Note: REQ_NULL is : ',mpl_request_null
!$OMP END CRITICAL(internal_write)
CALL IOS_print(IOS_message,src='ios_client_queue')
!$OMP CRITICAL(internal_write)
WRITE(IOS_message,'(A,I5)')  'State: ',IOS_Dispatch_Queue(slot)%state
!$OMP END CRITICAL(internal_write)
CALL IOS_print(IOS_message,src='ios_client_queue')
!$OMP CRITICAL(internal_write)
WRITE(IOS_message,'(A,I5)')  'Dest PE: ',IOS_Dispatch_Queue(slot)%pe
!$OMP END CRITICAL(internal_write)
CALL IOS_print(IOS_message,src='ios_client_queue')

IF (ASSOCIATED(IOS_Dispatch_Queue(slot)%payload)) THEN
!$OMP CRITICAL(internal_write)
  WRITE(IOS_message,'(A,I20)')'Payload Size: ',                            &
      SIZE(IOS_Dispatch_Queue(slot)%payload)
!$OMP END CRITICAL(internal_write)
  CALL IOS_print(IOS_message,src='ios_client_queue')
ELSE
  CALL IOS_print('No Payload',src='ios_client_queue')
END IF

CALL IOS_print('Metadata attached to slot is:',src='ios_client_queue')
CALL IOS_Report_MD(IOS_Dispatch_Queue(slot)%md)

END SUBROUTINE IOS_print_protocol_slot

END MODULE IOS_Client_Queue

