! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: C96

!
! Handle asynchronous stash requests on an io server
!

MODULE IOS_Stash_Server
USE io, ONLY:                                                               &
    setpos,                                                                  &
    setpos8,                                                                 &
    buffout,                                                                 &
    mpiio_fh
USE mask_compression, ONLY:                                                 &
    compress_to_mask
USE unite_output_files_mod, ONLY:                                           &
    unite_wgdos_files
USE IOS_Queue_Mod
USE IOS_Stash_Wgdos
USE IOS_Stash_Common
USE IOS_Model_Geometry, ONLY: getmaxfielddomain, atm_numprocs,               &
    atm_global_points, size_map, offset_map, land_mask
USE UM_types
USE mpl, ONLY: mpl_real, mpl_integer, mpl_integer4,                          &
    mpl_comm_world, mpl_status_size, mpl_integer8,                           &
    mpl_sum, mpl_request_null, mpl_real4
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE IOS_print_mgr, ONLY:                                                    &
    IOS_print,                                                               &
    IOS_Verbosity,                                                           &
    IOS_message,                                                             &
    IOS_PrStatus_Debug,                                                      &
    IOS_PrStatus_Diag,                                                       &
    IOS_PrStatus_Oper

USE errormessagelength_mod, ONLY: errormessagelength
USE get_wallclock_time_mod, ONLY: get_wallclock_time
USE Packing_Codes_Mod, ONLY: PC_No_Packing, PC_WGDOS_Packing

IMPLICIT NONE

INTERFACE IOS_dump_pack
MODULE PROCEDURE                                                           &
    IOS_dump_pack_1D,                                                      &
    IOS_dump_pack_2D
END INTERFACE

INTEGER, PARAMETER     :: stashLogUnit=10
INTEGER, POINTER       :: atm_to_global(:)
CHARACTER (LEN=errormessagelength),                                        &
    PRIVATE            :: iosStashServerMessage

! Debugging type parameters
LOGICAL                :: disable_packing
LOGICAL                :: disable_writes
LOGICAL                :: disable_subdomaining

! Statistics for MPI comms
INTEGER, PARAMETER     :: async_recv_calls  =1
INTEGER, PARAMETER     :: async_recv_bytes  =2
INTEGER, PARAMETER     :: async_recv_bytes2 =3
INTEGER, POINTER       :: async_bytes(:,:)

INTEGER,                                                                     &
    PARAMETER, PRIVATE :: rootPE=0

! params/vars  for dr_hook
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='IOS_STASH_SERVER'

CONTAINS

LOGICAL FUNCTION isIOS(x)
IMPLICIT NONE
INTEGER, INTENT(IN) :: x
INTEGER             :: i
INTEGER             :: j

isIOS=.FALSE.
DO i=1,IOS_Server_Groups
  DO j=1,model_procs
    IF (io_servers(i,j)==x)isIOS=.TRUE.
  END DO
END DO
END FUNCTION isIOS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialise the stash async subsystem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE IOS_stash_server_init(                                            &
    stash_active,                                                            &
    dumps_active,                                                            &
    disable_packing_in,                                                      &
    disable_writes_in,                                                       &
    disable_subdomaining_in)

IMPLICIT NONE

LOGICAL, INTENT(IN) :: disable_writes_in
LOGICAL, INTENT(IN) :: stash_active
LOGICAL, INTENT(IN) :: dumps_active
LOGICAL, INTENT(IN) :: disable_packing_in
LOGICAL, INTENT(IN) :: disable_subdomaining_in
INTEGER             :: f
INTEGER             :: l
INTEGER             :: ierr
INTEGER             :: myRank
INTEGER             :: i
INTEGER             :: k
REAL(KIND=jprb)     :: zhook_handle
CHARACTER (LEN=120) :: stashLogName
CHARACTER (LEN=*),                                                         &
    PARAMETER :: RoutineName = 'IOS_STASH_SERVER_INIT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

disable_packing=disable_packing_in
disable_writes=disable_writes_in
disable_subdomaining=disable_subdomaining_in
CALL useAsyncStash(stash_active)
CALL useAsyncDump (dumps_active)

IF (stash_active .OR. dumps_active) THEN

  atm_numprocs=global_procs-SIZE(io_servers)

  ALLOCATE(atm_to_global(0:atm_numprocs-1))!never deallocated

  IF (IOS_asyncDoStats) THEN
    ALLOCATE(async_bytes(3,0:atm_numprocs-1))!never deallocated
    async_bytes(:,:)=0
  END IF

  k=0
  DO i=0,global_procs-1
    IF (.NOT. isIOS(i)) THEN
      atm_to_global(k)=i
      k=k+1
    END IF
  END DO

!$OMP CRITICAL(internal_write)
  WRITE(IOS_message,'(A,I0,A,I0,A)')                                       &
      'Info: Stash Server: Initialised: There are ',                       &
      atm_numprocs,                                                        &
      ' atm procs and ',                                                   &
      SIZE(io_servers),' io procs'
!$OMP END CRITICAL(internal_write)
  CALL IOS_print(IOS_message,src='ios_stash_server')
  CALL mpl_comm_rank(mpl_comm_world,myRank,ierr)

  IF (model_rank == 0) THEN
!$OMP CRITICAL(internal_write)
    WRITE(stashLogName,'(A,I5.5)')'IOServer_stash_log.',myRank
    WRITE(IOS_message,'(A,A)')'Info: Stash Server: Logging to ',           &
        TRIM(stashLogName)
!$OMP END CRITICAL(internal_write)
    CALL IOS_print(IOS_message,src='ios_stash_server')
    OPEN(UNIT=stashLogUnit,FILE=TRIM(stashLogName))

    WRITE(stashLogUnit,'(A,A)')                                            &
        '   time hnd unt    position  datasize ',                          &
        ' disksize  blk fld full   S    N    W    E pack'
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE IOS_stash_server_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shut down the stash async subsystem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE IOS_stash_server_fini()

USE IOS_Server_Coupler, ONLY:                                             &
    proc_start,                                                            &
    proc_end

IMPLICIT NONE

INTEGER             :: proc
REAL                :: ave
REAL                :: var
REAL(KIND=jprb)     :: zhook_handle
CHARACTER (LEN=*),                                                         &
    PARAMETER :: RoutineName = 'IOS_STASH_SERVER_FINI'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (isUsingAsyncStash() .OR. isUsingAsyncDumps()) THEN

  IF (model_rank == 0)                                                     &
      CLOSE(stashLogUnit)

  IF (IOS_asyncDoStats) THEN

    CALL IOS_print(' ')
    CALL IOS_print('MPI counts for asynchronous operations')
    CALL IOS_print(' ')
    WRITE(IOS_message,'(A4,A10,A12,A12)')                                  &
        '----','----------','------------',                                &
        '------------'
    CALL IOS_print(IOS_message,src='ios_stash_server')
    WRITE(IOS_message,'(A4,A10,A12,A12)')                                  &
        ' CPU','  Receives','  Ave. Bytes',                                &
        '        S.D.'
    CALL IOS_print(IOS_message,src='ios_stash_server')
    WRITE(IOS_message,'(A4,A10,A12,A12)')                                  &
        '----','----------','------------',                                &
        '------------'
    CALL IOS_print(IOS_message,src='ios_stash_server')
    DO proc=proc_start,proc_end
      IF (async_bytes(async_recv_calls,proc)>0) THEN
        ave=1.0*async_bytes(async_recv_bytes,proc)/                        &
            async_bytes(async_recv_calls,proc)
        var=1.0*async_bytes(async_recv_bytes2,proc)/                       &
            async_bytes(async_recv_calls,proc)-ave*ave
      ELSE
        ave=0.0
        var=0.0
      END IF
!$OMP CRITICAL(internal_write)
      WRITE(IOS_message,'(I4,I0,F12.2,F12.2)')                             &
          proc,                                                            &
          async_bytes(async_recv_calls,proc),                              &
          ave,                                                             &
          SQRT(var)
!$OMP END CRITICAL(internal_write)
      CALL IOS_print(IOS_message,src='ios_stash_server')
    END DO
    WRITE(IOS_message,'(A4,A10,A12,A12)')                                  &
        '----','----------','------------',                                &
        '------------'
    CALL IOS_print(IOS_message,src='ios_stash_server')
    DEALLOCATE(async_bytes)
    NULLIFY(async_bytes)
  END IF

  DEALLOCATE(atm_to_global)
  NULLIFY(atm_to_global)

END IF
CALL IOS_print('Info: Async service terminated.')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE IOS_stash_server_fini

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Post receive calls for inbound components of field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE IOS_stash_server_init_recvs(node)!called by thr 0
USE mpl, ONLY:                                                            &
    mpl_request_null
USE IOS_Server_Coupler, ONLY:                                             &
    proc_start,                                                            &
    proc_end
USE IOS_Comms, ONLY:                                                      &
    acquire_lock_mpi,                                                      &
    release_lock_mpi
IMPLICIT NONE

TYPE(IOS_node_type),                                                       &
    INTENT(INOUT)  :: node
INTEGER            :: receive_buffer_size
INTEGER            :: proc
INTEGER            :: FieldType
INTEGER            :: j,i
INTEGER            :: control_block_len
INTEGER            :: num_records
INTEGER            :: num_fields
INTEGER            :: i_error
REAL(KIND=jprb)    :: zhook_handle
CHARACTER (LEN=*), PARAMETER                                               &
    :: RoutineName ='IOS_STASH_SERVER_INIT_RECVS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

receive_buffer_size=getMaxFieldDomain()*levels_in_pack(node)

! we need to allow for the queue to accommodate the receive
! buffers even if we dont allocate them here
! otherwise we may see a deadlock at the actual allocation time

IF (IOS_Verbosity>=IOS_PrStatus_Diag                                       &
    .AND. model_rank ==0) THEN
!$OMP CRITICAL(internal_write)
    WRITE(stashLogUnit,'(F8.2,2I4,A)')                                     &
    get_wallclock_time()-IOS_Start_Time,                                   &
    node%metadata%handle,node%metadata%UNIT,                               &
    ' Timing: Receive transaction'
!$OMP END CRITICAL(internal_write)
END IF


NULLIFY (node%receive_tags)
NULLIFY (node%receive_data_len)
NULLIFY (node%receive_requests)
NULLIFY (node%distributed_data)

IF (.NOT. use_blocking_recvs) THEN
  ! If we are using non-blocking then we need to track the
  ! request objects
  ALLOCATE( node%receive_requests(proc_start:proc_end))
  !deallocate in fini_recvs
  ALLOCATE( node%distributed_data                                          &
      (receive_buffer_size,proc_start:proc_end) )
  !deallocate in stash_server_deallocate
END IF

ALLOCATE (node%receive_tags(proc_start:proc_end))
ALLOCATE (node%receive_data_len(proc_start:proc_end))

! Pre calculation of tags/errors etc
DO proc=proc_start,proc_end
  node%receive_data_len(proc) = payload_size(node,proc)
  node%receive_tags(proc)     = receive_tag(node,proc)

  IF (node%receive_data_len(proc)<0) THEN
!$OMP CRITICAL(internal_write)
    WRITE(IOSStashServerMessage,'(A,I0,A)')                                &
        'recv ',                                                           &
        node%receive_data_len(proc),                                       &
        ' is negative '
!$OMP END CRITICAL(internal_write)
    CALL IOS_Ereport( RoutineName, 99, IOSStashServerMessage,              &
        md=node%metadata,UNIT=node%metadata%UNIT)
  END IF

  IF (node%receive_data_len(proc)>receive_buffer_size) THEN
!$OMP CRITICAL(internal_write)
    WRITE(IOSStashServerMessage,'(A,I0,A,I0,A,I0,A,I0)')                   &
        'recv ',                                                           &
        node%receive_data_len(proc),                                       &
        ' from ',                                                          &
        proc,'/',atm_numprocs,                                             &
        ' exceeds local buffer size of ',receive_buffer_size
!$OMP END CRITICAL(internal_write)
    CALL IOS_Ereport( RoutineName, 99, IOSStashServerMessage,              &
        md=node%metadata,UNIT=node%metadata%UNIT)
  END IF

  IF (IOS_asyncDoStats) THEN
    IF (node%receive_data_len(proc)>0) THEN

      async_bytes(async_recv_calls,proc)=                                  &
          async_bytes(async_recv_calls,proc)+1
      async_bytes(async_recv_bytes,proc)=                                  &
          async_bytes(async_recv_bytes,proc)+                              &
          (node%receive_data_len(proc)*8)
      async_bytes(async_recv_bytes2,proc)=                                 &
          async_bytes(async_recv_bytes2,proc)+                             &
          (node%receive_data_len(proc)*8)*                                 &
          (node%receive_data_len(proc)*8)

    END IF
  END IF

END DO

IF (.NOT. use_blocking_recvs) THEN
  CALL acquire_lock_mpi()

  DO proc=proc_start,proc_end

    IF (node%receive_data_len(proc)>0) THEN
      CALL mpl_irecv( node%distributed_data(1,proc),                       &
          node%receive_data_len(proc),                                     &
          mpl_real,                                                        &
          atm_to_global(proc),                                             &
          node%receive_tags(proc),                                         &
          IOS_Async_Comm, node%receive_requests(proc),                     &
          i_error)
    ELSE
      node%receive_requests(proc)=mpl_request_null
    END IF

  END DO

  CALL release_lock_mpi()

  DEALLOCATE (node%receive_tags)
  DEALLOCATE (node%receive_data_len)
  NULLIFY(node%receive_data_len)
  NULLIFY(node%receive_tags)

END IF

IF (lhook)CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE IOS_stash_server_init_recvs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Wait for data to arrive for the given node
!
! This function to be called by thread 1. The function returns when
! receives associated with the async op have completed. This function
! does not actually perform the receives.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE IOS_stash_server_wait_for_data(node)
USE fort2c_interfaces, ONLY: um_sleep
IMPLICIT NONE
TYPE(IOS_node_type),INTENT(INOUT) :: node
LOGICAL                           :: flag
REAL(KIND=jprb)                   :: zhook_handle
CHARACTER (LEN=*), PARAMETER                                               &
    :: RoutineName ='IOS_STASH_SERVER_WAIT_FOR_DATA'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (IOS_Verbosity>=IOS_PrStatus_Diag                                       &
    .AND. model_rank ==0) THEN
!$OMP CRITICAL(internal_write)
    WRITE(stashLogUnit,'(F8.2,2I4,A)')                                     &
    get_wallclock_time()-IOS_Start_Time,                                   &
    node%metadata%handle,node%metadata%UNIT,                               &
    ' Timing: Transaction at q head'
!$OMP END CRITICAL(internal_write)
END IF

IF (IOS_thread_0_calls_mpi) THEN
  ! We need to tell thread 0 that we need their help
  IOS_AsyncCompletionRequested=1
!$OMP FLUSH
  CALL IOS_print(                                                          &
      'Waiting for thread 0 to complete a request for stash',              &
      src='ios_stash_server')
  DO WHILE  (IOS_AsyncCompletionRequested==1)
    CALL um_sleep(IOS_backoff_interval)
!$OMP FLUSH
  END DO
ELSE
  ! This is called from thread 1, so flag can be ignored
  ! completion occurs on return
  flag=IOS_Stash_Server_Finish_Recvs(node)
END IF
IF (IOS_Verbosity>=IOS_PrStatus_Diag .AND. model_rank ==0) THEN
!$OMP CRITICAL(internal_write)
  WRITE(stashLogUnit,'(F8.2,2I4,A)')                                       &
      get_wallclock_time()-IOS_Start_Time,                                 &
      node%metadata%handle,node%metadata%UNIT,                             &
      ' Timing: Data transmission complete'
!$OMP END CRITICAL(internal_write)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE IOS_stash_server_wait_for_data


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Complete receiving inbound components of field
!
! If we are called from thread 0 we may return .FALSE.
! if the operation is not yet complete and use_blocking_receives
! is not set. This is so that the listener thread can continue
! listening for other work whilst the comms are ongoing. If we are
! called from thread N!=0 we will block on completion of recvs and
! the function will return always return .TRUE. or deadlock waiting.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
LOGICAL FUNCTION IOS_stash_server_finish_recvs(node)

USE IOS_Server_Coupler, ONLY:                                             &
    proc_start,                                                            &
    proc_end
USE IOS_Comms, ONLY:                                                      &
    acquire_lock_mpi,                                                      &
    release_lock_mpi,                                                      &
    IOS_WaitAll

IMPLICIT NONE

TYPE(IOS_node_type),INTENT(INOUT):: node

INTEGER,POINTER    :: STATUS(:,:)
INTEGER            :: proc
REAL               :: t1
REAL               :: t2
REAL               :: data_len
INTEGER            :: thread
LOGICAL            :: flag
INTEGER            :: mpiStat(mpl_status_size)
INTEGER            :: j
INTEGER            :: mymax
INTEGER            :: receive_buffer_size
INTEGER            :: i_error
REAL(KIND=jprb)    :: zhook_handle
CHARACTER (LEN=*), PARAMETER                                               &
    :: RoutineName ='IOS_STASH_SERVER_FINISH_RECVS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Set the thread to the main thread first, then overwrite it with the actual
! thread number if OpenMP is active (to allow this routine to work for non
! OpenMP builds)
thread = 0
!$    thread=omp_get_thread_num()
IOS_Stash_Server_Finish_Recvs=.FALSE.

IF (use_blocking_recvs) THEN
  receive_buffer_size=getMaxFieldDomain()*levels_in_pack(node)

  ALLOCATE(node%distributed_data                                           &
      (receive_buffer_size,proc_start:proc_end))
  !deallocate in stash_server_deallocate

  IF (.NOT. ASSOCIATED(node%receive_tags)) THEN
    WRITE(IOSStashServerMessage,'(A)')                                     &
        'MPI tag buffer not associated, exiting'
    CALL IOS_Ereport( RoutineName, 99, IOSStashServerMessage,              &
        md=node%metadata,UNIT=node%metadata%UNIT)
  END IF
  IF (.NOT. ASSOCIATED(node%receive_data_len)) THEN
    WRITE(IOSStashServerMessage,'(A)')                                     &
        'MPI length buffer not associated, exiting'
    CALL IOS_Ereport( RoutineName, 99, IOSStashServerMessage,              &
        md=node%metadata,UNIT=node%metadata%UNIT)

  END IF
  CALL acquire_lock_mpi()

  t1=get_wallclock_time()
  data_len=0
  DO proc=proc_start,proc_end
    IF (node%receive_data_len(proc)>0) THEN
      CALL mpl_recv( node%distributed_data(1,proc),                        &
          node%receive_data_len(proc),                                     &
          mpl_real, atm_to_global(proc),                                   &
          node%receive_tags(proc), IOS_async_comm,                         &
          mpiStat, i_error)

    END IF
    data_len=data_len+node%receive_data_len(proc)
  END DO
  t2=get_wallclock_time()

  data_len=data_len*8.0/1024.0/1024.0
  IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
    WRITE(IOS_message,'(A,F6.2,A,F8.2,A)')                                 &
        'Info: Recv: ',data_len,'MB @',data_len/(t2-t1),'MB/s'
    CALL IOS_print(IOS_message,src='ios_stash_server')
  END IF
  CALL release_lock_mpi()

  DEALLOCATE ( node%receive_data_len)!allocated in init_recvs
  DEALLOCATE ( node%receive_tags )!allocated in init_recvs
  NULLIFY(node%receive_data_len)
  NULLIFY(node%receive_tags)
  IOS_Stash_Server_Finish_Recvs=.TRUE.
ELSE
  ALLOCATE (STATUS(mpl_status_size,proc_start:proc_end))
                                   ! deallocate at end of block

  IF (.NOT. ASSOCIATED(node%receive_requests)) THEN
    WRITE(IOSStashServerMessage,'(A)')                                     &
        'MPI request buffer not associated, exiting'
    CALL IOS_Ereport( RoutineName, 99, IOSStashServerMessage,              &
        md=node%metadata,UNIT=node%metadata%UNIT)
  END IF
  IF (.NOT. ASSOCIATED(node%distributed_data)) THEN
    WRITE(IOSStashServerMessage,'(A)')                                     &
        'MPI receive buffer not associated, exiting'
    CALL IOS_Ereport( RoutineName, 99, IOSStashServerMessage,              &
        md=node%metadata,UNIT=node%metadata%UNIT)
  END IF
  IF (thread > 0) THEN
    CALL IOS_WaitAll(node%receive_requests, SIZE(node%receive_requests))
    IOS_Stash_Server_Finish_Recvs=.TRUE.
  ELSE! a CALL from thread 0 shouldnt block, nor do we need to lock
    CALL MPL_Testall(proc_end-proc_start+1,node%receive_requests,          &
        flag,STATUS,i_error)
    IOS_Stash_Server_Finish_Recvs=flag
  END IF
  IF (IOS_Stash_Server_Finish_Recvs) THEN
    ! If the operation completed deallocate the
    ! receive request objects
    DEALLOCATE(node%receive_requests)!allocated in init_recvs
    NULLIFY (node%receive_requests)
  END IF

  DEALLOCATE(STATUS)! allocated at start of block
  NULLIFY (STATUS)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END FUNCTION IOS_stash_server_finish_recvs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Process a stash object from the IO queue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE IOS_stash_server_process(node)

USE model_file, ONLY:                                                     &
    MF_Data_Address,                                                       &
    Get_File_Address,                                                      &
    attachLookups
USE IOS_Server_Coupler, ONLY:                                             &
    proc_start,                                                            &
    proc_end,                                                              &
    grid_row_start,                                                        &
    grid_row_end,                                                          &
    grid_point_start,                                                      &
    grid_point_end
USE IOS_Constants
USE IOS_geometry_utils
USE IOS_Comms, ONLY:                                                      &
    IOS_Gather,                                                           &
    acquire_lock_mpi,                                                     &
    release_lock_mpi

USE lookup_addresses

USE fort2c_pio_byteswap_interfaces, ONLY: get_machine_endianism,           &
     littleEndian, pio_byteswap
USE fort2c_data_conv_interfaces, ONLY: ibm2ieee, ieee2ibm, integer_type
USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_LOC

IMPLICIT NONE

TYPE(IOS_node_type)  :: node
INTEGER              :: FieldType
INTEGER,                                                                   &
    ALLOCATABLE      :: dataOffset(:)

INTEGER, POINTER     :: lookupTable(:,:)
INTEGER              :: lenOut32
INTEGER              :: lenOut64
INTEGER              :: ControlBlockLength
INTEGER              :: x_range
INTEGER              :: y_range
INTEGER              :: recordNumber
INTEGER              :: diskLocation
INTEGER              :: totalDiskSize
INTEGER              :: totalDataSize
INTEGER              :: currentTotalDataSize
INTEGER              :: numLandPoints
INTEGER, ALLOCATABLE :: ControlRecord(:)
INTEGER, ALLOCATABLE :: AutoControlRecord(:)
INTEGER, ALLOCATABLE :: partialLengths(:)
INTEGER, ALLOCATABLE :: displacements(:)
REAL, ALLOCATABLE    :: finalData(:)
REAL, ALLOCATABLE    :: gatheredData(:)
INTEGER              :: repeatCount
LOGICAL              :: doWrite
LOGICAL              :: firstRecordInPack
REAL,                                                                      &
    POINTER          :: field(:,:) !Stores the global field
REAL(KIND=real32), TARGET,                                                 &
    ALLOCATABLE      :: packedPartialData32(:)
REAL, ALLOCATABLE    :: landCompressedField(:)
TYPE(box)            :: domainBox

! Parameters from the input record
INTEGER              :: preprocess
INTEGER              :: packing
INTEGER              :: subdomainType
INTEGER              :: packingType
INTEGER              :: compressionAccuracy
INTEGER              :: landmaskCompression
INTEGER              :: fullField
INTEGER              :: n,s,e,w !inclusive bounds for subfields
INTEGER              :: DiskBlockSize  !write data in lumps of this
INTEGER              :: DiskBlockStart !write address
REAL                 :: missingDataIndicator

! identification of stash fields
INTEGER              :: section
INTEGER              :: code
INTEGER              :: level
INTEGER              :: dummy

! Various local counters/loop vars etc
INTEGER              :: j
INTEGER              :: iy
INTEGER              :: xlow
INTEGER              :: xhi
INTEGER              :: ylow
INTEGER              :: yhi
INTEGER              :: proc

! MPI-IO
INTEGER              :: offset_array(model_procs)
INTEGER              :: offset_reqs(model_procs)
INTEGER              :: mpl_status(MPL_STATUS_SIZE)
INTEGER              :: offset_statuses(MPL_STATUS_SIZE,model_procs)
INTEGER              :: offset_tag

INTEGER              :: start_of_write_offset
INTEGER              :: write_offset
INTEGER              :: end_of_write_offset

INTEGER              :: headerSendBuf(2)
INTEGER              :: headerRecvBuf(2)
INTEGER              :: headerSendBufSize

INTEGER              :: idxStashWrite
INTEGER              :: lenStashWrite
INTEGER              :: stashArrayExchangeReqs(2)
INTEGER              :: stashArrayExchangeTags(2)
INTEGER              :: stashArrayExchangeStatuses(MPL_STATUS_SIZE,2)

INTEGER, PARAMETER   :: idxSend = 1
INTEGER, PARAMETER   :: idxRecv = 2

INTEGER, PARAMETER   :: offsetTagBase = IOS_Request_Tag_Express +          &
     IOS_Request_Tag_Gap
INTEGER, PARAMETER   :: arrayExchangeTagBase = offsetTagBase + IOS_maxServers

INTEGER              :: dummyerr
INTEGER              :: ix

! Error codes, placeholders, and fluff
INTEGER              :: len_io
INTEGER              :: ErrorStatus
REAL                 :: IOSTAT
REAL(KIND=jprb)      :: zhook_handle
CHARACTER (LEN=*), PARAMETER                                               &
    :: RoutineName ='IOS_STASH_SERVER_PROCESS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

NULLIFY(lookupTable)
NULLIFY(field)

ALLOCATE(dataOffset(0:atm_numprocs-1))
ALLOCATE(AutoControlRecord(1:IOS_stash_control_auto_len))
ALLOCATE(ControlRecord(1:IOS_async_control_sz_max))
ALLOCATE(partialLengths(model_procs))
ALLOCATE(displacements(model_procs))

dataOffset(0:atm_numprocs-1)=1
firstRecordInPack=.TRUE.
repeatCount=0


! Just keep moving through the control block, we don't know how long it is
j=1
DO WHILE (j<node%metadata%data_size)
  AutoControlRecord(1:IOS_stash_control_auto_len) =                        &
      node%integer_data(j:j+IOS_stash_control_auto_len-1)
  ! IF this word is a start marker....
  IF (AutoControlRecord(loc_record_type)                                   &
      ==IOS_stash_record_start) THEN
    ControlBlockLength=                                                    &
        -1*AutoControlRecord(loc_record_len_control)
    j=j+IOS_stash_control_auto_len!1st element provided by user
    IF (node%integer_data(j)==IOS_Repeat_Record) THEN
      repeatCount=repeatCount+1

    ELSE IF (node%integer_data(j)==                                        &
        IOS_stash_distributed_field) THEN

      repeatCount=0
      ControlRecord(1:IOS_async_control_sz_max) =                          &
          node%integer_data(j:j+IOS_async_control_sz_max-1)
    ELSE
!$OMP CRITICAL(internal_write)
      WRITE(IOS_message,'(A,I0,I0)')'Bad Marker: ',                        &
          node%integer_data(j),ControlRecord(loc_record_type)
!$OMP END CRITICAL(internal_write)
      CALL IOS_print(IOS_message,src='ios_stash_server')
    END IF
    IF (ControlRecord(loc_record_type)==                                   &
        IOS_stash_distributed_field) THEN
      ! This looks like a field and we so we process it
      ! Lets get the contents into more meaningful variables
      FieldType           = ControlRecord(loc_fld_type)
      preprocess          = ControlRecord(loc_preprocess_flag)
      fullField           = ControlRecord(loc_subdomain_flag)
      packing             = ControlRecord(loc_packing_flag)
      packingType         = ControlRecord(loc_pack_type)
      compressionAccuracy = ControlRecord(loc_comp_accry)
      DiskBlockSize       = ControlRecord(loc_disk_block)
      DiskBlockStart      = ControlRecord(loc_disk_block_start)
      recordNumber        = ControlRecord(loc_seek_target)
      landmaskCompression = ControlRecord(loc_landmaskcompress)
      CALL IOS_unpack4(controlRecord(loc_boundary),n,s,e,w)
      ! recover that real value from the integer array...
      CALL um_memcpy_f                                                     &
          (MissingDataIndicator,ControlRecord(loc_dmi),1)

      ! Do not use the control record below here as we will modify values

      ! Because we may have a repeat record, update recordNumber
      recordNumber=recordNumber+repeatCount

      IF (FieldType < 1 .OR. FieldType > nfld_types ) THEN
!$OMP CRITICAL(internal_write)
        WRITE(IOSStashServerMessage,'(A,I0)')                              &
            'ERROR INVALID FLD TYPE ',FieldType
!$OMP END CRITICAL(internal_write)
        CALL IOS_Ereport( RoutineName, 99, IOSStashServerMessage,          &
            md=node%metadata,UNIT=node%metadata%UNIT)
      END IF

      ALLOCATE(field(atm_global_points(1,FieldType),                       &
          atm_global_points(2,FieldType)))

      ! If the disable subdomaining flag is set, then
      ! zero the field so that parts we don't have data for
      ! (ie aren't in the requested domain) aren't garbage.

      IF (disable_subdomaining)                                            &
          field(:,:)=missingDataIndicator


      ! ****************************************************
      ! Extract into field data from everyone who sent us
      ! some component of the field
      ! ****************************************************

      !Assemble the field from the multitude of receive buffers
      domainBox%n = n
      domainBox%s = s
      domainBox%e = e
      domainBox%w = w

      DO proc=proc_start,proc_end

        IF (fullField == IOS_partial_field) THEN
          subdomainType=classify_subdomain(proc,FieldType,domainBox)
        ELSE
          subdomainType=complete_intersection
        END IF

        IF (subdomainType /= no_intersection                               &
            .OR. IOS_AsyncSendNull) THEN

          xlow=offset_map(1,FieldType,proc)
          xhi=offset_map(1,FieldType,proc)+                                &
              size_map(1,FieldType,proc)-1
          ylow=offset_map(2,FieldType,proc)
          yhi=offset_map(2,FieldType,proc)+                                &
              size_map(2,FieldType,proc)-1


          DO iy=ylow,yhi
            field(xlow:xhi,iy)=                                            &
                node%distributed_data(dataOffset(proc):                    &
                dataOffset(proc)+xhi-xlow,proc)
            dataOffset(proc)=dataOffset(proc)+xhi-xlow+1
          END DO
        END IF
      END DO

      IF (preprocess==IOS_stash_preprocess) THEN
        ! Some kind of data processing required

        IF (node%metadata%ACTION==IOS_Action_StashWritePPData              &
            .AND. disable_packing) THEN
          IF (packing==ios_packing) THEN
            CALL IOS_print(                                                &
                'Warning: DEBUG: Disabling packing option active!',        &
                src='ios_stash_server')
            packing=ios_no_packing
            packingType=0
          END IF
        END IF

        doWrite=.TRUE.
        IF (disable_writes) THEN
!$OMP CRITICAL(internal_write)
          WRITE(IOS_message,'(A,I0)')                                      &
              'Warning: DEBUG: Disabling data writes for unit ',           &
              node%metadata%UNIT
!$OMP END CRITICAL(internal_write)
          CALL IOS_print(IOS_message,src='ios_stash_server')
          doWrite=.FALSE.
        END IF

        IF (node%metadata%ACTION==IOS_Action_StashWritePPData              &
            .AND. disable_subdomaining) THEN
          IF (fullField==ios_partial_field) THEN
            CALL IOS_print(                                                &
                'Warning: DEBUG: Disabling '//                             &
                'subdomaining option active!',                             &
                src='ios_stash_server')
            fullField=ios_full_field
          END IF
        END IF

        ! ****************************************************
        ! From our assembled field construct a correctly sized
        ! array and copy in the data
        ! subdomain = INTERSECTION(req. domain, this IOS task)
        ! ****************************************************

        ! Ensure N,S,E,W are sane for full fields
        IF (fullField==IOS_full_field) THEN
          s=1
          n=atm_global_points(2,FieldType)
          w=1
          e=atm_global_points(1,FieldType)
        END IF

        CALL IOS_stash_server_subdomain(                                   &
            field,        & ! on input a pointer to the global field array
            fullField,    & ! whether or not this is a full field request
            FieldType,    & !
            MIN(n,grid_row_end  (fieldtype)),                              &
            MAX(s,grid_row_start(fieldtype)),                              &
            e,                                                             &
            w,                                                             &
            x_range,                                                       &
            y_range)


        ! If there was a subdomain 'field' now points to the
        ! smaller region.

        IF (node%metadata%ACTION==IOS_Action_StashWritePPData) THEN

          ! ****************************************************
          ! Operations for stash data
          ! ****************************************************


          ! Always call this to make section and code available later
            CALL IOS_unpack4(controlRecord(loc_id),                        &
                section,                                                   &
                code,                                                      &
                level,                                                     &
                dummy)
          IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
!$OMP CRITICAL(internal_write)
            WRITE(IOS_message,'(A,I3,A,I0)')                               &
                'Info: Stash Server: unit=',                               &
                node % metadata % UNIT,                                    &
                ' field number=',recordNumber
!$OMP END CRITICAL(internal_write)
            CALL IOS_print(IOS_message,src='ios_stash_server')
!$OMP CRITICAL(internal_write)
            WRITE(IOS_message,'(A,I0,A,I0)')                               &
                'Info: Stash Server: Section= ',                           &
                section,' Code= ',code
!$OMP END CRITICAL(internal_write)
            CALL IOS_print(IOS_message,src='ios_stash_server')
          END IF

          ! Pack as needed

          ALLOCATE (packedPartialData32                                    &
              (2*(x_range+packedheaderwords32)*y_range+packedheaderwords32+1))

          CALL IOS_stash_pack(                                             &
              field,                                                       &
              packedPartialData32,                                         &
              packing,                                                     &
              packingType,                                                 &
              compressionAccuracy,                                         &
              MissingDataIndicator,                                        &
              lenOut32,                                                    &
              lenOut64,                                                    &
              FieldType,                                                   &
              x_range,                                                     &
              y_range,                                                     &
              1000*section+code,                                           &
              DiskBlockSize)

          IF ( IOS_enable_mpiio ) THEN
             IF ( packing == IOS_packing ) THEN
                ! Reduce data in local WGDOS headers on to rank 0 in order to
                ! generate global WGDOS header on rank 0, and to inform all 
                ! IO tasks of the size of the record.
                dummyerr=ibm2ieee(integer_type,1,packedPartialData32(1), 0,  &
                         headerSendBuf(1),1,BIT_SIZE(headerSendBuf(1)),32)
                dummyerr=ibm2ieee(integer_type,1,packedPartialData32(3), 16, &
                         headerSendBuf(2),1,BIT_SIZE(headerSendBuf(2)),16)
                ! Remove header size on all but rank 0
                IF ( model_rank > 0 )                                        &
                     headerSendBuf(1) = headerSendBuf(1) - packedheaderwords32
                headerSendBufSize = 2
             ELSE
                ! If packing is not used, populate the headerSendBuf
                ! with lenOut32
                headerSendBuf(1) = lenOut32
                headerSendBufSize = 1
             END IF
             ! When MPI3 is approved, replace with MPL_Iallreduce
             IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
             CALL MPL_Allreduce(headerSendBuf, headerRecvBuf,              &
                  headerSendBufSize, MPL_INTEGER, MPL_SUM, model_comm,     &
                  errorStatus )
             IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()
             ! At the start of each pack, there is no guarantee that the 
             ! file pointer is where we left it at the end of the last pack.
             ! In the case of CRUNs, there may not have been a last pack.
             IF ( firstRecordInPack ) THEN
                IF ( model_rank == 0 ) diskLocation =                       &
                     IOS_stash_server_seek(node%metadata%UNIT, recordNumber)
                IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
                CALL MPL_Bcast( diskLocation, 1, MPL_INTEGER, 0, model_comm,&
                     errorStatus )
                IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()
                IF ( model_rank /= 0 )                                      &
                     CALL setpos( node%metadata%UNIT, diskLocation,         &
                     errorStatus )
                firstRecordInPack = .FALSE.
             END IF
             ! Recover position in bytes (even if we just set it)
             IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
             CALL MPL_File_get_position( mpiio_fh(node%metadata%UNIT),     &
                  start_of_write_offset, ErrorStatus )
             ! Set disk location for the IOS stash server log file
             diskLocation=start_of_write_offset/IOS_BytesPerWord64
             ! Distribute file offsets
             DO proc=1,model_procs
                IF ( proc == model_rank+1 ) THEN
                   offset_array(proc) = 0
                   offset_reqs(proc)  = MPL_REQUEST_NULL  
                ELSE IF ( proc < model_rank+1 ) THEN
                   offset_tag = offsetTagBase + proc
                   CALL MPL_Irecv(offset_array(proc),1,MPL_INTEGER,proc-1,  &
                        offset_tag,model_comm,offset_reqs(proc),            &
                        ErrorStatus )
                ELSE
                   offset_array(proc) = headerSendBuf(1) ! = lenOut32 on 
                                                         ! rank 0, lenOut32
                                                         ! - 3 on other ranks
                   offset_tag = offsetTagBase + model_rank + 1
                   CALL MPL_Isend(offset_array(proc),1,MPL_INTEGER,proc-1,   &
                        offset_tag,model_comm,offset_reqs(proc),ErrorStatus )
                END IF
             END DO
             IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()
                ! Ensure all data written will be multiples of 64 bits, 
                ! rearrange arrays if necessary to ensure correct endianness
             IF ( packing == IOS_packing ) THEN
                stashArrayExchangeReqs = MPL_REQUEST_NULL
                ! Need offset values for proper endian conversion, so can't
                ! hide communication cost with the endian conversion.
                IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
                CALL MPL_Waitall( model_procs, offset_reqs,              &
                     offset_statuses, ErrorStatus )
                IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()     
                ! Calculate per-rank write offsets in 32-bits for now, will
                ! convert to byte offsets after any array exchanges are complete
                write_offset = SUM(offset_array(1:model_rank))
                IF ( get_machine_endianism()==littleEndian ) THEN
                   ! Rank 0 always writes entire array, and potentially an
                   ! extra value if its data length is not a multiple of 64 bits
                   IF ( model_rank == 0 ) THEN
                      idxStashWrite = 1
                      ! Need the initial element from rank 1
                      IF ( MODULO(lenOut32,2) == 1 ) THEN
                         lenStashWrite=(lenOut32+1)/2
                         IF ( model_procs > 1 ) THEN
                            stashArrayExchangeTags(idxRecv) =                  &
                                 arrayExchangeTagBase + model_rank + 1
                         ! Only receive if there is more than 1 IO task
                            IF ( IOS_serialise_mpi_calls )                     &
                                 CALL acquire_lock_mpi()
                            ! Receive to what will be the last element of
                            ! packedPartialData32 after endian change
                            CALL MPL_Irecv( packedPartialData32(lenOut32), 1,  &
                                 MPL_REAL4,model_rank + 1,                     &
                                 stashArrayExchangeTags(idxRecv),              &
                                 model_comm,stashArrayExchangeReqs(idxRecv),   &
                                 errorStatus )
                            IF ( IOS_serialise_mpi_calls )                     &
                                 CALL release_lock_mpi()
                         END IF
                      ELSE
                         lenStashWrite = lenOut32/2
                      END IF
                      ! Construct global WGDOS header
                      dummyerr=ieee2ibm(integer_type,1,packedPartialData32(1), &
                           0,headerRecvBuf(1),1,BIT_SIZE(headerSendBuf(1)),32)
                      dummyerr=ieee2ibm(integer_type,1,packedPartialData32(3), &
                           16,headerRecvBuf(2),1,BIT_SIZE(headerSendBuf(2)),16)
                   ELSE
                      ! For the other ranks, determine if the write starts at
                      ! a multiple of 64 bits or not.
                      IF( MODULO(write_offset,2)==0 ) THEN
                         ! This rank is writing on a 64-bit boundary to disk,
                         ! however, it will not be writing the 96 bit header, 
                         ! so will not be writing on a 64-bit boundary inside
                         ! packedPartialData32. Reordering needs to take place
                         ! to ensure data is correct after endian conversion
                         DO ix=2,headerSendBuf(1),2
                            packedPartialData32(ix) = packedPartialData32(ix+4)
                         END DO
                         IF ( MODULO(headerSendBuf(1),2) == 1 ) THEN
                            IF( model_rank == model_procs - 1 ) THEN
                            ! Special case for last proc in IO server - 
                            ! Last (after endian swap) element contains
                            ! garbage, so zero it.
                               packedPartialData32(lenOut32-2) = 0
                            ELSE
                            ! All other ranks receive from rank+1
                               stashArrayExchangeTags(idxRecv) =              &
                                    arrayExchangeTagBase + model_rank + 1
                               IF ( IOS_serialise_mpi_calls )                 &
                                    CALL acquire_lock_mpi()
                               ! Receive to what will be the last written 
                               ! element of the reordered packedPartialData32 
                               ! after endian change.
                               CALL MPL_Irecv(                                &
                                    packedPartialData32(lenOut32-2),          &
                                    1, MPL_REAL4,model_rank+1,                &
                                    stashArrayExchangeTags(idxRecv),          &
                                    model_comm,                               &
                                    stashArrayExchangeReqs(idxRecv),          &
                                    errorStatus )
                               IF ( IOS_serialise_mpi_calls )                 &
                                    CALL release_lock_mpi()
                            END IF
                            lenStashWrite=(lenOut32-(packedheaderwords32-1))/2
                         ELSE
                            ! No action taken on packedPartialData32
                            ! when writing to 64-bit boundary
                            lenStashWrite=(lenOut32-packedheaderwords32)/2
                         END IF
                         ! Since packedPartialData32 has been reordered, write
                         ! starts from what was inside the WGDOS header
                         idxStashWrite=2
                      ELSE
                         ! write_offset%2 != 0
                         ! Need to send first element of actual data to 
                         ! model_rank - 1
                         stashArrayExchangeTags(idxSend) =                    &
                              arrayExchangeTagBase + model_rank
                         IF ( IOS_serialise_mpi_calls )                       &
                              CALL acquire_lock_mpi()
                         ! WGDOS header info is first 3 elements of 
                         ! packedPartialData32 after endian change, so 
                         ! first data value is located at index 3 before 
                         ! endian change
                         CALL MPL_Isend(packedPartialData32(3), 1, MPL_REAL4, &
                              model_rank-1, stashArrayExchangeTags(idxSend),  &
                              model_comm, stashArrayExchangeReqs(idxSend),    &
                              errorStatus)
                         IF ( IOS_serialise_mpi_calls )                       &
                              CALL release_lock_mpi()
                         ! If this rank holds a multiple of 64 bits in 
                         ! packedPartialData32, it also needs to receive
                         IF ( MODULO(headerSendBuf(1),2) == 0 ) THEN
                            ! Unless it is the last IO task in its server
                            IF ( model_rank /= model_procs -1 ) THEN
                               stashArrayExchangeTags(idxRecv) =              &
                                    arrayExchangeTagBase + model_rank + 1
                               IF ( IOS_serialise_mpi_calls )                 &
                                    CALL acquire_lock_mpi()
                               ! Receive to what will be the last element of
                               ! packedPartialData32 after endian change
                               CALL MPL_Irecv(                                &
                                    packedPartialData32(lenOut32),            &
                                    1, MPL_REAL4,model_rank + 1,              &
                                    stashArrayExchangeTags(idxRecv),          &
                                    model_comm,                               &
                                    stashArrayExchangeReqs(idxRecv),          &
                                    errorStatus )
                               IF ( IOS_serialise_mpi_calls )                 &
                                    CALL release_lock_mpi()
                            END IF
                            lenStashWrite = (lenOut32-packedheaderwords32)/2
                         ELSE
                            ! Don't need to receive, but since model_rank - 1
                            ! is writing part of this ranks array, make
                            ! sure it doesn't over write 
                            IF ( model_rank == model_procs - 1 ) THEN
                               lenStashWrite = (lenOut32-packedheaderwords32)/2
                            ELSE
                               lenStashWrite =                                &
                                    (lenOut32-(packedheaderwords32+1))/2
                            END IF
                         END IF
                         ! Set array index for write start, and correct offset
                         ! to account for shifted write
                         idxStashWrite = packedheaderwords32 + 2
                         write_offset = write_offset + 1
                      END IF
                   END IF
                   ! Wait for outstanding send/recv requests (if any)
                   IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
                   CALL MPL_Waitall( 2, stashArrayExchangeReqs,                &
                        stashArrayExchangeStatuses, errorStatus )
                   IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()
                   ! Change packedPartialData32 endian
                   ErrorStatus = pio_byteswap(                                 &
                        C_LOC(packedPartialData32(idxStashWrite)),             &
                        lenStashWrite,IOS_BytesPerWord64 )
                ELSE
                   IF ( model_rank == 0 ) THEN
                      idxStashWrite = 1
                      lenStashWrite = (lenOut32+1)/2
                   ELSE
                      idxStashWrite = packedheaderwords32+1
                      lenStashWrite = (lenOut32-packedheaderwords32+1)/2
                   END IF
                END IF
                ! Convert write_offset to bytes from start of file.
                write_offset = write_offset * IOS_BytesPerWord32 +             &
                     start_of_write_offset
             ELSE
                IF ( get_machine_endianism()==littleEndian )                   &
                     ErrorStatus = pio_byteswap(                               &
                     C_LOC(packedPartialData32),lenOut64,IOS_BytesPerWord64 )
                idxStashWrite = 1
                lenStashWrite = lenOut64
                IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
                CALL MPL_Waitall( model_procs, offset_reqs,                    &
                     offset_statuses, ErrorStatus )
                IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()
                ! Calculate per-rank write offsets
                write_offset = start_of_write_offset +                         &
                     SUM(offset_array(1:model_rank)) * IOS_BytesPerWord32
             END IF
             IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
             CALL MPL_File_write_at_all(mpiio_fh(node%metadata%unit),          &
                  write_offset,packedPartialData32(idxStashWrite),             &
                  lenStashWrite,MPL_Integer8,mpl_status,errorStatus)
             IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()
             ! Set end of write offset to nearest multiple of DiskBlockSize
             ! headerRecvBuf(1) is in 32-bits, DiskBlockSize is in 64-bits
             end_of_write_offset = ( (headerRecvBuf(1) + 2*DiskBlockSize - 1) /&
                  (2*DiskBlockSize) ) * 2*DiskBlockSize
             ! Add on to start_of_write_offset
             end_of_write_offset = end_of_write_offset * IOS_BytesPerWord32 +  &
                  start_of_write_offset
             ! Move local pointer to new location
             CALL setpos8( node%metadata%UNIT, end_of_write_offset, errorStatus)
             ! Logging/header calculations
             IF ( model_rank == 0 ) THEN
                totalDataSize = (headerRecvBuf(1)+1) / 2
                totalDiskSize = (end_of_write_offset - start_of_write_offset)  &
                     / IOS_BytesPerWord64
             END IF
          ELSE

             CALL IOS_gather(partialLengths,lenOut64)

             IF (model_rank==0) THEN
                displacements(1)=0
                DO proc=2,model_procs
                   displacements(proc)=displacements(proc-1)+                  &
                        partialLengths(proc-1)
                END DO

                totalDataSize=SUM(partialLengths(:))

                ALLOCATE(finalData((totalDataSize)+DiskBlockSize))
                ALLOCATE(gatheredData((totalDataSize)+DiskBlockSize))
             ELSE
                ALLOCATE(finalData(1))
                ALLOCATE(gatheredData(1))
             END IF

             CALL  IOS_SS_Stash_Gather(                                       &
                  packedPartialData32,                                         &
                  lenOut64,                                                    &
                  gatheredData,                                                &
                  partialLengths)

             IF (model_rank==0) THEN

            ! Move the file pointer to the right place
                diskLocation=                                                  &
                     IOS_stash_server_seek                                     &
                     (node%metadata%UNIT, recordNumber)

                IF ( packing == IOS_packing ) THEN

                   DO proc=1,model_procs
                      CALL unite_wgdos_files(                                 &
                           gatheredData(displacements(proc)+1),              &
                           finaldata,                                        &
                           totalDataSize,                                    &
                           proc-1)
                   END DO

                   ! convert totalDataSize to whole 64 bit words
                   totalDataSize=(totalDataSize+1)/2

                   ! Write it onto disk
                   CALL IOS_stash_server_write(                              &
                        node,                                                &
                        finalData,                                           &
                        totalDataSize,                                       &
                        DiskBlockSize,                                       &
                        totalDiskSize,                                       &
                  doWrite)
                ELSE
                   CALL IOS_stash_server_write(                              &
                        node,                                                &
                        gatheredData,                                        &
                        totalDataSize,                                       &
                        DiskBlockSize,                                       &
                        totalDiskSize,                                       &
                        doWrite)
                   
                END IF

             END IF
             DEALLOCATE(finalData)
             DEALLOCATE(gatheredData)

          END IF
            ! ****************************************************
            ! For stash, update the lookup tables
            ! ****************************************************

            ! now: output_size    = data length (words)
            !      totalDiskSize  = amount written to disk (words)
            !      diskLocation   = location in file (words)
          IF ( model_rank == 0 ) THEN
            lookupTable=>attachLookups(node % metadata % UNIT)

            lookupTable(lbegin , recordNumber)=diskLocation
            lookupTable(lblrec , recordNumber)=totalDataSize
            lookupTable(lbnrec , recordNumber)=totalDiskSize

            ! Note for fields files NADDR is the same as LBEGIN
            lookupTable(naddr  , recordNumber)=diskLocation
            lookupTable(lbpack , recordNumber)=PC_No_Packing
            IF (packing==IOS_packing) THEN
              lookupTable(lbpack , recordNumber)=PC_WGDOS_Packing
            END IF

            ! If we disabled subdomains then we should also
            ! correct entries in the lookup which will be wrong
            ! compared to what the atmos model will later send.
            IF (disable_subdomaining) THEN
              lookupTable(lbrow , recordNumber)=                           &
                  atm_global_points(2,FieldType)
              lookupTable(lbnpt , recordNumber)=                           &
                  atm_global_points(1,FieldType)
              lookupTable(bzx , recordNumber)= 0
              lookupTable(bzy , recordNumber)= 0
              lookupTable(lbhem , recordNumber)= 0
              CALL IOS_print(                                              &
                  'Warning: Exact details of grid are unknown',            &
                  src='ios_stash_server')
              WRITE(IOS_message,'(A,A)')                                   &
                  '         Field will be written',                        &
                  ' as global 0,0 origined'
              CALL IOS_print(IOS_message,src='ios_stash_server')

              ! We should set LBHEM, but we don't have enough data to do so.
            END IF

            NULLIFY(lookupTable)

            ! Make a note of what happened in the log.
            CALL IOS_stash_server_log(node,FieldType,                      &
                preprocess, fullField,packing,s,n,w,e,                     &
                MissingDataIndicator,DiskBlockSize,                        &
                diskLocation,totalDiskSize*8,(totalDataSize+1)/2)
            
          END IF
         
          DEALLOCATE(packedPartialData32)
          
        ELSE IF (node%metadata%ACTION ==                                   &
            IOS_Action_StashWriteDumpData) THEN

          ! ****************************************************
          ! Operations for dump data
          ! ****************************************************

          IF (landmaskCompression==IOS_packing_type_landmask) THEN
            ! DiskBlockSize is the size of the record on disk
            ! not the sector size as it is for stash
            ALLOCATE(landCompressedField(x_range*y_range))
            CALL compress_to_mask(                                         &
                field,                                                     &
                landCompressedField,                                       &
                land_mask(                                                 &
                grid_point_start(fieldType):                               &
                grid_point_end(fieldType)),                                &
                x_range*y_range,                                           &
                numLandPoints)
            ALLOCATE (packedPartialData32(2*numLandPoints))
            CALL IOS_dump_pack(                                            &
                landCompressedField,                                       &
                packedPartialData32,      &! field in/out
                packing,                                                   &
                packingType,              &! packing args
                lenOut32,                                                  &
                lenOut64,                 &! data lengths out
                FieldType,                                                 &
                numLandPoints,                                             &
                DiskBlockSize)
            DEALLOCATE(landCompressedField)
          ELSE
            ALLOCATE (packedPartialData32(2*x_range*y_range))
            CALL IOS_dump_pack(                                            &
                field,                                                     &
                packedPartialData32,      & ! field in/out
                packing,                  & ! packing args
                packingType,                                               &
                lenOut32,                 & ! data lengths out
                lenOut64,                                                  &
                FieldType,                                                 &
                x_range,                                                   &
                y_range,                                                   &
                DiskBlockSize)

          END IF

          IF (model_rank==0 .OR. IOS_Enable_mpiio) THEN
            IF (firstRecordInPack) THEN
              ! we need to setpos
              IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
!$OMP CRITICAL(internal_write)
                WRITE(IOS_message,'(A,A,I0)')                              &
                    'Info: Stash Server: ',                                &
                    'Seek position for dump output=',                      &
                    DiskBlockStart
!$OMP END CRITICAL(internal_write)
                CALL IOS_print(IOS_message,src='ios_stash_server')
              END IF
              CALL setpos(node%metadata%UNIT,DiskBlockStart,               &
                  errorStatus)
              firstRecordInPack=.FALSE.
            END IF
          END IF

          IF ( IOS_enable_mpiio ) THEN
             ! Send data lengths for offset calculations
             CALL acquire_lock_mpi()
             DO proc=1,model_procs
                IF ( proc == model_rank+1 ) THEN
                   offset_array(proc) = 0
                   offset_reqs(proc)  = MPL_REQUEST_NULL  
                ELSE IF ( proc < model_rank+1 ) THEN
                   offset_tag = offsetTagBase + proc
                   CALL MPL_Irecv(offset_array(proc),1,MPL_INTEGER,proc-1,  &
                        offset_tag,model_comm,offset_reqs(proc),ErrorStatus )
                ELSE
                   offset_array(proc) = lenOut32
                   offset_tag = offsetTagBase + model_rank + 1
                   CALL MPL_Isend(offset_array(proc),1,MPL_INTEGER,proc-1,  &
                        offset_tag,model_comm,offset_reqs(proc),ErrorStatus )
                END IF
             END DO             
             ! Recover the local pointer location and calculate the 
             ! offset for the end of the record in bytes
             CALL MPL_File_get_position( mpiio_fh(node%metadata%UNIT),   &
                  start_of_write_offset, ErrorStatus )
             CALL release_lock_mpi()
             end_of_write_offset = start_of_write_offset +               &
                     DiskBlockSize*IOS_BytesPerWord64
             IF (packingType==ios_packing_type_pack21) THEN
                ! Convert endian if necessary (32 bit)
                IF ( get_machine_endianism()==littleEndian )             &
                     ErrorStatus = pio_byteswap(                         &
                     C_LOC(packedPartialData32),lenOut32,IOS_BytesPerWord32 )
                CALL acquire_lock_mpi()
                CALL MPL_Waitall( model_procs, offset_reqs,              &
                     offset_statuses, ErrorStatus )
                ! Calculate per-rank write offsets
                write_offset = start_of_write_offset +                   &
                     SUM(offset_array(1:model_rank)) * IOS_BytesPerWord32
                ! Write data (32-bit)
                CALL MPL_File_write_at_all( mpiio_fh(node%metadata%UNIT),&
                     write_offset, PackedPartialData32, lenOut32,        &
                     MPL_Integer4, mpl_Status, errorStatus )
                CALL release_lock_mpi()
             ELSE
                ! Convert endian if necessary (64 bit)
                IF ( get_machine_endianism()==littleEndian )             &
                     ErrorStatus = pio_byteswap(                         &
                     C_LOC(packedPartialData32),lenOut64,IOS_BytesPerWord64 )
                CALL acquire_lock_mpi()
                CALL MPL_Waitall( model_procs, offset_reqs,              &
                     offset_statuses, ErrorStatus )
                ! Calculate per-rank write offsets
                write_offset = start_of_write_offset +                   &
                     SUM(offset_array(1:model_rank)) * IOS_BytesPerWord32
                ! Write data (64-bit)
                CALL MPL_File_write_at_all( mpiio_fh(node%metadata%UNIT),&
                     write_offset, packedPartialData32, lenOut64,        &
                     MPL_Integer8, mpl_Status, errorStatus )
                CALL release_lock_mpi()
             END IF
             ! Update the local pointer to the end of the write
             CALL setpos8( node%metadata%UNIT, end_of_write_offset,      &
                      errorStatus)
          ELSE

             CALL IOS_gather(partialLengths,lenOut32)

             IF (model_rank==0) THEN
                displacements(1)=0
                DO proc=2,model_procs
                   displacements(proc)=                                    &
                        displacements(proc-1)+partialLengths(proc-1)
                END DO
                totalDataSize=SUM(partialLengths(:))
                ALLOCATE(finalData(DiskBlockSize))
             ELSE
                ALLOCATE(finalData(1))
             END IF
             
             CALL  IOS_SS_Dump_Gather(                                        &
                  packedPartialData32,                                        &
                  lenOut32,                                                   &
                  finalData,                                                  &
                  partialLengths)

             IF (packingType==IOS_packing_type_pack21) THEN
                
                IF (model_rank==0) THEN
                   ! DEPENDS ON: buffout32_f77
                   CALL buffout32_f77(node%metadata%UNIT,                    &
                        finalData,                                           &
                        DiskBlockSize*2,                                     &
                        len_io,IOSTAT)
                END IF
                
             ELSE
                IF (model_rank==0) THEN
                   CALL buffout(node%metadata%UNIT,                          &
                        finalData,                                           &
                        DiskBlockSize,                                       &
                        len_io,IOSTAT)
                END IF
            
             END IF

             DEALLOCATE(finalData)
             
          END IF

          DEALLOCATE(packedPartialData32)
          

        ELSE ! Wrong action somehow
          CALL IOS_Ereport(RoutineName,99,                                 &
              'Wrong action in Stash_Server',                              &
              md=node%metadata,UNIT=node%metadata%UNIT)

        END IF

      ELSE

!$OMP CRITICAL(internal_write)
        WRITE(IOSStashServerMessage,'(A,I0)')                              &
            'unknown flag passed for preprocessing: ',preprocess
!$OMP END CRITICAL(internal_write)
        CALL IOS_Ereport( RoutineName, 99, IOSStashServerMessage,          &
            md=node%metadata,UNIT=node%metadata%UNIT)
      END IF
      DEALLOCATE(field)
      NULLIFY (field)
    ELSE
!$OMP CRITICAL(internal_write)
      WRITE(IOSStashServerMessage,'(A,I0)')                                &
          'Unknown data object in control record at ',j
!$OMP END CRITICAL(internal_write)
      CALL IOS_Ereport( RoutineName, 99, IOSStashServerMessage,            &
          md=node%metadata,UNIT=node%metadata%UNIT)
    END IF ! is a distributed field
    ! Advance to the next record
    j=j+ControlBlockLength
  ELSE
    !We didnt get a start of record marker where expected :-(
!$OMP CRITICAL(internal_write)
    WRITE(IOSStashServerMessage,'(A,I0)')                                  &
        'process: unexpected lack of record ',j
!$OMP END CRITICAL(internal_write)
    CALL IOS_Ereport( RoutineName, 99, IOSStashServerMessage,              &
        md=node%metadata,UNIT=node%metadata%UNIT)
  END IF
  IF (IOS_Verbosity>=IOS_PrStatus_Debug) THEN
!$OMP CRITICAL(internal_write)
    WRITE(IOS_message,'(A,I0,A)')'Info: Stash Server: unit=',              &
        node % metadata % UNIT,                                            &
        ' done record'
!$OMP END CRITICAL(internal_write)
    CALL IOS_print(IOS_message,src='ios_stash_server')
  END IF
END DO

IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
!$OMP CRITICAL(internal_write)
  WRITE(IOS_message,'(A,I0,A)')'Info: Stash Server: unit=',                &
      node % metadata % UNIT,                                              &
      ' done pack'
!$OMP END CRITICAL(internal_write)
  CALL IOS_print(IOS_message,src='ios_stash_server')
END IF

DEALLOCATE(dataOffset)
DEALLOCATE(AutoControlRecord)
DEALLOCATE(ControlRecord)
DEALLOCATE(partialLengths)
DEALLOCATE(displacements)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE IOS_stash_server_process

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set the file position to the location for a record
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION IOS_stash_server_seek(UNIT, recordNumber)                           &
    RESULT(seekAddress)

USE model_file, ONLY:                                                     &
    get_file_address,                                                      &
    MF_data_address,                                                       &
    attachLookups
USE lookup_addresses

IMPLICIT NONE

INTEGER, INTENT(IN)         :: UNIT
INTEGER, INTENT(IN)         :: recordNumber
INTEGER, POINTER            :: fixed_header(:)
INTEGER, POINTER            :: ipplook(:,:)
INTEGER                     :: seekAddress
INTEGER                     :: dataAddress
INTEGER                     :: ierror
REAL(KIND=jprb)             :: zhook_handle
CHARACTER (LEN=*), PARAMETER                                               &
    :: RoutineName ='IOS_STASH_SERVER_SEEK'

NULLIFY(fixed_header)
NULLIFY(ipplook)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

ipplook=>attachLookups(UNIT)
dataAddress=get_file_address(UNIT,MF_data_address)
IF (recordNumber < 1) THEN
  CALL IOS_Ereport('Stash_seek',99,'invalid record number',                &
      UNIT=unit)

ELSE IF (recordNumber==1) THEN
  seekAddress = dataAddress
ELSE
  ! previous record location + prev record size
  seekAddress =                                                            &
      ipplook(lbegin, recordNumber-1)+                                     &
      ipplook(lbnrec, recordNumber-1)
END IF

IF (seekAddress < MAX(0,dataAddress)) THEN
!$OMP CRITICAL(internal_write)
  WRITE(IOSStashServerMessage,'(A,I0,A,I0,A)')                             &
      'Seek value (',seekAddress,                                          &
      ') lower than data start address (',dataAddress,                     &
      ') from fixed header'
!$OMP END CRITICAL(internal_write)
  CALL IOS_Ereport('stash_server_seek',99,                                 &
      IOSStashServerMessage,UNIT=unit)
END IF

IF (IOS_Verbosity>=IOS_PrStatus_Debug) THEN
!$OMP CRITICAL(internal_write)
  WRITE(IOS_message,'(A,I0,A,I0)')                                         &
      'Info: Stash Server: decoded record ',recordNumber,                  &
      ' to disk address ',seekAddress
!$OMP END CRITICAL(internal_write)
  CALL IOS_print(IOS_message,src='ios_stash_server')
END IF
CALL setpos(UNIT,seekAddress, ierror)
NULLIFY(ipplook)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END FUNCTION IOS_stash_server_seek

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Perform compression on a stash object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE IOS_stash_pack(                                                   &
    input,                                                                   &
    output,                                                                  &
    packing,                                                                 &
    packingType,                                                             &
    compressionAccuracy,                                                     &
    MissingDataIndicator,                                                    &
    lenOut32,                                                                &
    lenOut64,                                                                &
    FieldType,                                                               &
    x_range,                                                                 &
    y_range,                                                                 &
    StashItem,                                                               &
    DiskBlockSize)

USE IOS_stash_common
USE IOS_stash_wgdos
USE fort2c_memcpy_interfaces, ONLY: um_memcpy64

IMPLICIT NONE
INTEGER, INTENT(IN)          :: DiskBlockSize
INTEGER, INTENT(IN)          :: compressionAccuracy
INTEGER, INTENT(IN)          :: packingType
INTEGER, INTENT(IN)          :: FieldType
INTEGER, INTENT(IN)          :: x_range,y_range
INTEGER, INTENT(IN)          :: StashItem
INTEGER, INTENT(INOUT)       :: packing
INTEGER, INTENT(OUT)         :: lenOut32
INTEGER, INTENT(OUT)         :: lenOut64
REAL, INTENT(IN)             :: MissingDataIndicator
REAL, POINTER                :: input(:,:)
REAL(KIND=real32),                                                         &
    INTENT(OUT)              :: output(:)
REAL(KIND=jprb)              :: zhook_handle
CHARACTER (LEN=*), PARAMETER                                               &
    :: RoutineName ='IOS_STASH_PACK'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF (packing==IOS_packing) THEN
  IF (packingType==IOS_packing_type_wgdos) THEN
    IF (compressionAccuracy>-99) THEN
      CALL IOS_stash_pack_wgdos(input,output,                              &
          lenOut32,compressionAccuracy,                                    &
          StashItem, MissingDataIndicator)
      IF (lenout32==0) THEN
        WRITE(IOS_message,'(A,A)')'Info: Stsh: Packing Failed...',         &
            ' writing field unpacked.'
        CALL IOS_print(IOS_message,src='ios_stash_server')
      END IF
    ELSE
      IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
!$OMP CRITICAL(internal_write)
        WRITE(IOS_message,'(A,I0)')                                        &
            'Info: Stash Server: No packing: accuracy=',                   &
            compressionAccuracy
!$OMP END CRITICAL(internal_write)
        CALL IOS_print(IOS_message,src='ios_stash_server')
      END IF
      packing=IOS_no_packing
    END IF
  ELSE IF (packingType==ios_packing_type_pack21) THEN
    WRITE(IOS_message,'(A,A)')'Info: Stsh: PACK2D: pac21 packing...',      &
        'not implemented'
    CALL IOS_print(IOS_message,src='ios_stash_server')
    packing=IOS_no_packing
  ELSE
    CALL IOS_print('Info: Stsh: Unknown packing code',                     &
        src='ios_stash_server')
    packing=ios_no_packing
  END IF
  lenOut64=(lenOut32+1)/2
END IF

! We may enter this because WGDOS packing failed.
IF (packing==IOS_no_packing) THEN
  ! Literal copy of real data into integer output storage.
  CALL um_memcpy64(output,input,x_range*y_range)
  lenOut64=x_range*y_range
  lenOut32=lenOut64*2
END IF

IF (packing/=IOS_packing .AND. packing/=IOS_no_packing) THEN
!$OMP CRITICAL(internal_write)
  WRITE(IOSStashServerMessage,'(A,I0)')                                    &
      'Packing flag incorrectly set ',packing
!$OMP END CRITICAL(internal_write)
  CALL IOS_Ereport( RoutineName, 99, IOSStashServerMessage )
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE IOS_stash_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Perform compression on a dump object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE IOS_dump_pack_2D(                                                 &
    input,                                                                   &
    output32,                                                                &
    packing,                                                                 &
    packingType,                                                             &
    lenOut32,                                                                &
    lenOut64,                                                                &
    FieldType,                                                               &
    x_range,                                                                 &
    y_range,                                                                 &
    DiskBlockSize)

USE IOS_stash_common
USE IOS_stash_wgdos
USE fort2c_memcpy_interfaces, ONLY: um_memcpy64

IMPLICIT NONE
REAL, INTENT(IN)             :: input(:,:)
REAL(KIND=real32)            :: output32(:)
INTEGER, INTENT(OUT)         :: lenOut32
INTEGER, INTENT(OUT)         :: lenOut64
INTEGER, INTENT(INOUT)       :: packing
INTEGER, INTENT(IN)          :: packingType
INTEGER, INTENT(IN)          :: FieldType
INTEGER, INTENT(IN)          :: x_range
INTEGER, INTENT(IN)          :: y_range
INTEGER                      :: DiskBlockSize
REAL(KIND=jprb)              :: zhook_handle
CHARACTER (LEN=*), PARAMETER                                               &
    :: RoutineName ='IOS_DUMP_PACK_2D'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF (packing==IOS_packing) THEN
  IF (packingType==IOS_packing_type_pack21) THEN
    CALL pack21(x_range*y_range,input,output32)
  ELSE
    CALL IOS_print('unknown packing code',src='ios_stash_server')
    packing=ios_no_packing
  END IF
  lenOut32=x_range*y_range
  lenOut64=(lenOut32+1)/2
END IF

! We may enter this because of an unknown packing code
IF (packing==IOS_no_packing) THEN
  ! Literal copy of real data into integer output storage.
  CALL um_memcpy64(output32,input,x_range*y_range)
  lenOut64=x_range*y_range
  lenOut32=lenOut64*2
END IF

IF (packing/=IOS_packing .AND. packing/=IOS_no_packing) THEN
!$OMP CRITICAL(internal_write)
  WRITE(IOSStashServerMessage,'(A,I0)')                                    &
      'Packing flag incorrectly set ',packing
!$OMP END CRITICAL(internal_write)
  CALL IOS_Ereport( RoutineName, 99, IOSStashServerMessage )
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE IOS_dump_pack_2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Perform compression on an object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE IOS_dump_pack_1D(input,output32,                                  &
    packing,packingType,                                                     &
    lenOut32,lenOut64,                                                       &
    FieldType,num_points,DiskBlockSize)
USE IOS_stash_common
USE IOS_stash_wgdos
USE fort2c_memcpy_interfaces, ONLY: um_memcpy64

IMPLICIT NONE
INTEGER, INTENT(IN)          :: packingType
INTEGER, INTENT(IN)          :: FieldType
INTEGER, INTENT(IN)          :: num_points
INTEGER, INTENT(IN)          :: DiskBlockSize
INTEGER, INTENT(OUT)         :: lenOut32
INTEGER, INTENT(OUT)         :: lenOut64
INTEGER, INTENT(INOUT)       :: packing
REAL, INTENT(IN)             :: input(:)
REAL(KIND=real32)            :: output32(:)
REAL(KIND=jprb)              :: zhook_handle
CHARACTER (LEN=*), PARAMETER                                               &
    :: RoutineName ='IOS_DUMP_PACK_1D'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF (packing==IOS_packing) THEN
  IF (packingType==IOS_packing_type_pack21) THEN
    CALL pack21(num_points,input,output32)
  ELSE
    CALL IOS_print('Unknown packing code (pack 1D)',                       &
        src='ios_stash_server')
    packing=ios_no_packing
  END IF
  lenOut32 = num_points
  lenOut64=(lenOut32+1)/2
END IF

! We may enter this because WGDOS packing failed.
IF (packing==IOS_no_packing) THEN
  ! Literal copy of real data into integer output storage.

  CALL um_memcpy64(output32,input,num_points)
  lenOut64=num_points
  lenOut32=lenOut64*2

END IF

IF (packing/=IOS_packing .AND. packing/=IOS_no_packing) THEN
!$OMP CRITICAL(internal_write)
  WRITE(IOSStashServerMessage,'(A,I0)')                                    &
      'Packing flag incorrectly set ',packing
!$OMP END CRITICAL(internal_write)
  CALL IOS_Ereport( RoutineName, 99, IOSStashServerMessage )
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE IOS_dump_pack_1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write a stash record to disk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE IOS_stash_server_write(node,outData,data_length,                  &
    block_size,written,do_write)

IMPLICIT NONE

TYPE(IOS_node_type),                                                       &
    INTENT(INOUT)        :: node
REAL, INTENT(IN)         :: outData(:)
INTEGER, INTENT(IN)      :: data_length
INTEGER, INTENT(IN)      :: block_size
LOGICAL, INTENT(IN)      :: do_write
INTEGER, INTENT(OUT)     :: written
INTEGER                  :: io_len
REAL                     :: rstat
REAL(KIND=jprb)          :: zhook_handle
CHARACTER (LEN=*), PARAMETER                                               &
    :: RoutineName ='IOS_STASH_SERVER_WRITE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

io_len =((data_length+block_size-1)/block_size)*block_size

IF (do_write) THEN
  CALL buffout( node%metadata%UNIT, outData,                               &
      io_len, written, rstat )
  IF (written /= io_len .OR. rstat /= -1.0) THEN
    IF (written /= io_len) THEN
!$OMP CRITICAL(internal_write)
      WRITE(IOSStashServerMessage,'(A,A,I0,I0)')                           &
          'IOS:server:stash: ',                                            &
          'Mismatch between bytes written/requested',                      &
          written,io_len
!$OMP END CRITICAL(internal_write)
      CALL IOS_Ereport( RoutineName, 99, IOSStashServerMessage )
    END IF
    IF (rstat /= -1.0) THEN
      WRITE(IOSStashServerMessage,'(A,F8.4,F8.4)')                         &
          'IOS:server:stash: Bad Status:',rstat,-1.0
      CALL IOS_Ereport( RoutineName, 99, IOSStashServerMessage )
    END IF
  END IF
ELSE
  CALL IOS_print('DEBUG: Skipping write per debug setting ',               &
      src='ios_stash_server')
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE IOS_stash_server_write

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Deallocate memory associated with a stash object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE IOS_stash_server_deallocate(node)

USE IOS_Server_Coupler, ONLY:                                             &
    procs
USE ios_common, ONLY: needing_lock
IMPLICIT NONE

TYPE(IOS_node_type),INTENT(INOUT) :: node
INTEGER                           :: receive_buffer_size
REAL(KIND=jprb)                   :: zhook_handle
CHARACTER (LEN=*), PARAMETER                                               &
    :: RoutineName ='IOS_STASH_SERVER_DEALLOCATE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF (ASSOCIATED(node%distributed_data)) THEN

  receive_buffer_size=                                                     &
      getMaxFieldDomain()*                                                 &
      levels_in_pack(node)*                                                &
      procs

  DEALLOCATE(node%distributed_data)
  NULLIFY(node%distributed_data)

  ! Release space associated with receive buffers
  CALL IOS_Increment_Queue_Length                                          &
      (-1*receive_buffer_size*IOS_BytesPerReal,needing_lock)
ELSE
  WRITE(IOSStashServerMessage,'(A)')                                       &
      'Deallocation requested, but no allocation was present'
  CALL IOS_Ereport( RoutineName, -99, IOSStashServerMessage )
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE IOS_stash_server_deallocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute the receive tag to use for a given request
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER FUNCTION receive_tag(node,processor)

IMPLICIT NONE
TYPE(IOS_node_type), INTENT(IN) :: node
INTEGER, INTENT(IN)             :: processor

receive_tag=node%metadata%handle

END FUNCTION receive_tag


INTEGER FUNCTION levels_in_pack(node)
IMPLICIT NONE
TYPE(IOS_node_type),                                                       &
    INTENT(IN)      :: node
INTEGER             :: num_records
INTEGER             :: num_fields
INTEGER             :: fieldtype
INTEGER             :: i
INTEGER             :: j
INTEGER             :: ControlBlockLength
REAL(KIND=jprb)     :: zhook_handle
CHARACTER (LEN=*), PARAMETER                                               &
    :: RoutineName ='LEVELS_IN_PACK'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

j=1
levels_in_pack=0
num_records=0
num_fields=0
DO WHILE(j<node%metadata%data_size)
  IF (node%integer_data(j)==IOS_stash_record_start) THEN!new record
    num_records=num_records+1
    ControlBlockLength=-1*node%integer_data(j+loc_record_len_control-1)

    !Advance to next record
    j=j+IOS_stash_control_auto_len+ControlBlockLength
    num_fields=num_fields+1
  ELSE
!$OMP CRITICAL(internal_write)
    WRITE(IOSStashServerMessage,'(A,I0)')                                  &
        'levels_in_pack: Control Block FAULT: j=',j
!$OMP END CRITICAL(internal_write)
    CALL IOS_Ereport( RoutineName, 99, IOSStashServerMessage,              &
        md=node%metadata,UNIT=node%metadata%UNIT)
  END IF
END DO

levels_in_pack=num_records

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END FUNCTION levels_in_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Convenience function to return the total payload accross all
! processors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER FUNCTION total_payload_size(node)

IMPLICIT NONE
TYPE(IOS_node_type),                                                       &
    INTENT(IN)      :: node
INTEGER             :: processor
total_payload_size=0
DO processor=0,atm_numprocs-1
  total_payload_size=total_payload_size+payload_size(node,processor)
END DO
END FUNCTION total_payload_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Figure out how many items there are in a package from a
! given cpu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER FUNCTION payload_size(node,processor)

USE IOS_geometry_utils

IMPLICIT NONE

TYPE(IOS_node_type),                                                       &
    INTENT(IN)      :: node
INTEGER, INTENT(IN) :: processor
INTEGER             :: j,num_records,num_fields
INTEGER             :: FieldType
INTEGER             :: i
INTEGER             :: n,s,e,w
INTEGER             :: FieldDomain
INTEGER             :: subdomainType
INTEGER             :: ControlBlockLength
INTEGER             :: last
TYPE(box)           :: domainBox
REAL(KIND=jprb)     :: zhook_handle
CHARACTER (LEN=*), PARAMETER                                               &
    :: RoutineName ='PAYLOAD_SIZE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

j=1
payload_size=0
num_records=0
num_fields=0
last=0
DO WHILE (j < node%metadata%data_size)
  IF (node%integer_data(j)==IOS_stash_record_start) THEN!new record
    num_records=num_records+1
    ControlBlockLength=-1*node%integer_data(j+loc_record_len_control-1)

    !Advance to 1st element provided by user
    j=j+IOS_stash_control_auto_len

    IF (node%integer_data(j)==IOS_stash_distributed_field) THEN

      FieldType   = node%integer_data(j+loc_fld_type-1)
      FieldDomain = node%integer_data(j+loc_subdomain_flag-1)
      CALL IOS_unpack4(node%integer_data(j+loc_boundary-1),                &
          domainBox%n,                                                     &
          domainBox%s,                                                     &
          domainBox%e,                                                     &
          domainBox%w)

      IF (fieldDomain == IOS_partial_field) THEN
        subdomainType=                                                     &
            classify_subdomain(processor,FieldType,domainBox)
      ELSE
        subdomainType=complete_intersection
      END IF

      IF (FieldType < 1 .OR. FieldType > nfld_types) THEN
!$OMP CRITICAL(internal_write)
        WRITE(IOSStashServerMessage,'(A,I0)')                              &
            'ERROR INVALID FIELD TYPE ',FieldType
!$OMP END CRITICAL(internal_write)
        DO i=1,SIZE(node%integer_data)
!$OMP CRITICAL(internal_write)
          WRITE(IOS_message,'(I8,A,I0,A,I0)')i,'/',                        &
              SIZE(node%integer_data),' value: ',                          &
              node%integer_data(i)
!$OMP END CRITICAL(internal_write)
          CALL IOS_print(IOS_message,src='ios_stash_server')
        END DO
        CALL IOS_Ereport( RoutineName, 99, IOSStashServerMessage,          &
            md=node%metadata,UNIT=node%metadata%UNIT)
      ELSE
        num_fields=num_fields+1
        IF (subdomainType /= no_intersection                               &
            .OR. IOS_AsyncSendNull) THEN
          payload_size=payload_size+                                       &
              size_map(1,FieldType,processor)*                             &
              size_map(2,FieldType,processor)
          last=                                                            &
              size_map(1,FieldType,processor)*                             &
              size_map(2,FieldType,processor)
        ELSE
          last=0
        END IF
      END IF
    ELSE IF (node%integer_data(j)==IOS_repeat_record) THEN
      payload_size=payload_size+last
    ELSE
!$OMP CRITICAL(internal_write)
      WRITE(IOSStashServerMessage,'(A,I0)')                                &
          'Unknown data object in control record at ',j
!$OMP END CRITICAL(internal_write)
      DO i=1,SIZE(node%integer_data)
!$OMP CRITICAL(internal_write)
        WRITE(IOS_message,'(I0,A,I0,A,I0)')                                &
            i,'/',SIZE(node%integer_data),' value: ',                      &
            node%integer_data(i)
!$OMP END CRITICAL(internal_write)
        CALL IOS_print(IOS_message,src='ios_stash_server')
      END DO
      CALL IOS_Ereport( RoutineName, 99, IOSStashServerMessage,            &
          md=node%metadata,UNIT=node%metadata%UNIT)
    END IF
    j=j+ControlBlockLength
  ELSE
!$OMP CRITICAL(internal_write)
    WRITE(IOSStashServerMessage,'(A,I0)')                                  &
        'payload_size: Control Block FAULT: j=',j
!$OMP END CRITICAL(internal_write)
    DO i=1,SIZE(node%integer_data)
!$OMP CRITICAL(internal_write)
      WRITE(IOS_message,'(I0,A,I0,A,I0)')i,'/',                            &
          SIZE(node%integer_data),' value: ',                              &
          node%integer_data(i)
!$OMP END CRITICAL(internal_write)
      CALL IOS_print(IOS_message,src='ios_stash_server')
    END DO
    CALL IOS_Ereport( RoutineName, 99, IOSStashServerMessage,              &
        md=node%metadata,UNIT=node%metadata%UNIT)
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END FUNCTION payload_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Log to disk what we just did
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE IOS_stash_server_log(node,FieldType,preprocess,                   &
    fullField,packing,s,n,w,e,                                               &
    MissingDataIndicator,DiskBlockSize,location,sz,                          &
    data_len)

IMPLICIT NONE

TYPE(IOS_node_type) :: node
INTEGER,INTENT(IN)  :: DiskBlockSize!write data in lumps of this
INTEGER,INTENT(IN)  :: n,s,e,w
INTEGER,INTENT(IN)  :: fieldtype
INTEGER,INTENT(IN)  :: packing
INTEGER,INTENT(IN)  :: fullField
INTEGER,INTENT(IN)  :: preprocess
INTEGER,INTENT(IN)  :: location
INTEGER,INTENT(IN)  :: sz
INTEGER,INTENT(IN)  :: data_len
REAL   ,INTENT(IN)  :: MissingDataIndicator
REAL                :: t
REAL(KIND=jprb)     :: zhook_handle
CHARACTER (LEN=*), PARAMETER                                               &
    :: RoutineName ='IOS_STASH_SERVER:IOS_STASH_SERVER_LOG'

t=get_wallclock_time()-IOS_Start_Time
IF (IOS_Verbosity>=IOS_PrStatus_Oper .AND. model_rank ==0) THEN
!$OMP CRITICAL(internal_write)
  WRITE(stashLogUnit,'(F8.2,I5,I4,A,I0,A,I0,A,I0,A,2I4,6I5)')              &
      t,                                                                   &
      node%metadata%handle,                                                &
      node%metadata%UNIT,                                                  &
      ' ',location,                                                        &
      ' ',data_len*8,                                                      &
      ' ',sz,                                                              &
      ' ',DiskBlockSize,fieldtype,                                         &
      fullField,s,n,w,e,packing
!$OMP END CRITICAL(internal_write)
END IF

END SUBROUTINE IOS_stash_server_log

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Pull out the part of the field wanted for subdomaining
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE  IOS_stash_server_subdomain(global,fullField,                     &
    FieldType,n,s,e,w,x_range,y_range)

IMPLICIT NONE

REAL, POINTER         :: working(:,:)
REAL, POINTER         :: global (:,:)
INTEGER,INTENT(IN)    :: n,s,w
INTEGER,INTENT(INOUT) :: e
INTEGER,INTENT(IN)    :: fullField
INTEGER,INTENT(IN)    :: FieldType
INTEGER,INTENT(OUT)   :: x_range,y_range
INTEGER               :: first_slice
INTEGER               :: second_slice
REAL(KIND=jprb)       :: zhook_handle
CHARACTER (LEN=*), PARAMETER                                               &
    :: RoutineName ='IOS_STASH_SERVER_SUBDOMAIN'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

NULLIFY(working)
x_range=0
y_range=0

IF (fullField==IOS_full_field .AND. model_procs==1) THEN
  ! the subdomain is the same as the global field.
  x_range=atm_global_points(1,FieldType)
  y_range=atm_global_points(2,FieldType)

ELSE ! For other cases we need to create a new object
     ! representing the subdomain

  ! Normalise our E/W ordering.
  IF (w > e) THEN
    e=e+atm_global_points(1,FieldType)
  END IF

  IF (s > n) THEN ! The subdomain is empty

    y_range=0
    x_range=0
    ALLOCATE(working(1:0,1:0)) ! a zero lengthed array

  ELSE IF (e > atm_global_points(1,FieldType)) THEN

    x_range=e-w+1
    y_range=n-s+1
    first_slice  = atm_global_points(1,FieldType)-w+1
    second_slice = e-atm_global_points(1,FieldType)
    ALLOCATE(working(1:x_range,1:y_range))
    working(1:first_slice,:)=                                              &
        global(w:atm_global_points(1,FieldType),s:n)
    working(first_slice+1:x_range,:)=                                      &
        global(1:e-atm_global_points(1,FieldType),s:n)

  ELSE

    x_range=e-w+1
    y_range=n-s+1
    ALLOCATE(working(1:x_range,1:y_range))
    working(:,:)=global(w:e,s:n)

  END IF

  ! As the subdomain is now allocated and populated, retarget
  ! the pointer at the correctly sized object.
  DEALLOCATE (global)
  global => working
  NULLIFY (working)

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE IOS_stash_server_subdomain


SUBROUTINE IOS_SS_Stash_Gather(myData,myLength,gatheredData,                 &
    partialLengths)
USE IOS_Comms, ONLY:                                                      &
    acquire_lock_mpi,                                                      &
    release_lock_mpi,                                                      &
    IOS_WaitAll
USE fort2c_memcpy_interfaces, ONLY: um_memcpy64

IMPLICIT NONE
INTEGER, INTENT(IN)  :: myLength
REAL(KIND=real32),                                                         &
    INTENT(IN)       :: myData(:)
INTEGER, INTENT(IN)  :: partialLengths(:)
REAL, INTENT(OUT)    :: gatheredData(:)
INTEGER              :: errorStatus
INTEGER              :: tag
INTEGER              :: recv_offset
INTEGER              :: proc
INTEGER              :: num_requests
INTEGER, ALLOCATABLE :: requests(:)
CHARACTER (LEN=*), PARAMETER                                               &
    :: RoutineName ='IOS_STASH_SERVER:IOS_SS_STASH_GATHER'

tag=2
IF (model_rank==0) THEN

  CALL um_memcpy64(gatheredData,myData,myLength)
  recv_offset=myLength

  num_requests = model_procs - 1
  ALLOCATE (requests(num_requests))

  IF (model_procs>1) THEN
    CALL acquire_lock_mpi()
    DO proc=1,model_procs-1
      CALL MPL_iRecv(                                                      &
          gatheredData(1+recv_offset),                                     &
          partialLengths(proc+1),                                          &
          mpl_integer,                                                     &
          proc,                                                            &
          tag,                                                             &
          model_comm,                                                      &
          requests(proc),                                                  &
          errorStatus)
      recv_offset=recv_offset+partialLengths(proc+1)
    END DO
    CALL release_lock_mpi()
  END IF
ELSE

  num_requests = 1
  ALLOCATE (requests(num_requests))

  CALL acquire_lock_mpi()
  CALL MPL_Isend(                                                          &
      myData,                                                              &
      myLength,                                                            &
      mpl_integer,                                                         &
      rootPE,                                                              &
      tag,                                                                 &
      model_comm,                                                          &
      requests(1),                                                         &
      errorStatus)
  CALL release_lock_mpi()
  gatheredData(:)=-1
END IF

CALL IOS_WaitAll(requests, num_requests)

DEALLOCATE (requests)

END SUBROUTINE IOS_SS_Stash_Gather

SUBROUTINE IOS_SS_Dump_Gather(myData,myLength,gatheredData,partialLengths)
USE IOS_Comms, ONLY:                                                      &
    acquire_lock_mpi,                                                      &
    release_lock_mpi,                                                      &
    IOS_WaitAll
USE fort2c_memcpy_interfaces, ONLY: um_memcpy32
IMPLICIT NONE
INTEGER, INTENT(IN)  :: myLength          ! In 32 bit words
REAL(KIND=real32),                                                         &
    INTENT(IN)       :: myData(:)         ! Note this is 64 bit buffer
INTEGER, INTENT(IN)  :: partialLengths(:) ! In 32 bit words
REAL, INTENT(OUT)    :: gatheredData(:)   ! Note this is 64 bit buffer

REAL(KIND=real32),                                                         &
    ALLOCATABLE      :: tempData(:)       ! 32 bit temp buffer
INTEGER              :: errorStatus
INTEGER              :: tag
INTEGER              :: recv_offset
INTEGER              :: proc
INTEGER              :: num_requests
INTEGER, ALLOCATABLE :: requests(:)
CHARACTER (LEN=*), PARAMETER                                               &
    :: RoutineName ='IOS_STASH_SERVER:IOS_SS_DUMP_GATHER'

tag=3
IF (model_rank==0) THEN
  ALLOCATE(tempData(SIZE(gatheredData)*2))
  CALL um_memcpy32(tempData,myData,myLength)
  recv_offset=myLength

  num_requests = model_procs-1
  ALLOCATE (requests(model_procs-1))

  IF (model_procs>1) THEN
    CALL acquire_lock_mpi()
    DO proc=1,model_procs-1
      CALL MPL_iRecv(                                                      &
          tempData(1+recv_offset),                                         &
          partialLengths(proc+1),                                          &
          mpl_integer4,                                                    &
          proc,                                                            &
          tag,                                                             &
          model_comm,                                                      &
          requests(proc),                                                  &
          errorStatus)
      recv_offset=recv_offset+partialLengths(proc+1)
    END DO
    CALL release_lock_mpi()
  END IF
ELSE

  num_requests = 1
  ALLOCATE (requests(num_requests))

  CALL acquire_lock_mpi()
  CALL MPL_Isend(                                                          &
      myData,                                                              &
      myLength,                                                            &
      mpl_integer4,                                                        &
      rootPE,                                                              &
      tag,                                                                 &
      model_comm,                                                          &
      requests(1),                                                         &
      errorStatus)
  CALL release_lock_mpi()
  gatheredData(:)=-1
END IF

CALL IOS_WaitAll(requests, num_requests)

IF (model_rank==0) THEN
  CALL um_memcpy32(gatheredData,tempData,SUM(partialLengths))
  DEALLOCATE(tempData)
END IF

DEALLOCATE (requests)

END SUBROUTINE IOS_SS_Dump_Gather

END MODULE IOS_Stash_Server


