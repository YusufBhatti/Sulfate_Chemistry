! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: C96

!---------------------------------------------------------------------
! Initialisation of IOS, called whether or not IOS to be used
!---------------------------------------------------------------------

MODULE IOS_Init
USE UM_Types
USE IOS_Common
!$ USE omp_lib
USE IOS_Comms, ONLY:                                                         &
    IOS_Comms_init,                                                          &
    IOS_Comms_fini
USE IOS_Queue_Mod,      ONLY:                                                &
    IOS_Queue_Initialise,                                                    &
    IOS_Queue_Finalise
USE IO_Server_listener, ONLY:                                                &
    IOS_Listener,                                                            &
    IOS_SlaveListener
USE IO_Server_writer,   ONLY:                                                &
    IOS_Writer
USE IOS_Client_Queue
USE IOS_Stash,          ONLY:                                                &
    IOS_Stash_Client_Init
USE IOS_Stash_Server,   ONLY:                                                &
    IOS_Stash_Server_Init,                                                   &
    IOS_Stash_Server_Fini
USE IOS_Stash_Common,   ONLY:                                                &
    IOS_Async_Comm,                                                          &
    IOS_asyncMaxFieldsInPack,                                                &
    IOS_asyncSendNull,                                                       &
    IOS_asyncDoStats,                                                        &
    IOS_asyncnumslots
USE IOS_MPI_error_handlers
USE missing_data_mod, ONLY: imdi
USE ereport_mod, ONLY: ereport
USE umprintmgr, ONLY: newline
USE IOS_print_mgr, ONLY:                                                     &
    IOS_PrintSetLevel,                                                       &
    IOS_print,                                                               &
    IOS_Verbosity,                                                           &
    IOS_print_start_time,                                                    &
    IOS_message
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength

USE get_wallclock_time_mod, ONLY: get_wallclock_time
 
IMPLICIT NONE

! DrHook-related parameters.
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

INTEGER, PRIVATE                :: IOS_spacing
INTEGER, PRIVATE                :: IOS_buffer_size
INTEGER, PRIVATE                :: IOS_Offset
INTEGER, PRIVATE                :: IOS_as_concurrency
INTEGER, PRIVATE                :: IOS_async_levs_per_pack
INTEGER, PRIVATE                :: IOS_force_threading_mode
! where:
!  0 = single
!  1 = funneled
!  2 = serialized
!  3 = multiple
LOGICAL, PRIVATE                :: IOS_Interleave
LOGICAL, PRIVATE                :: IOS_use_async_stash
LOGICAL, PRIVATE                :: IOS_use_async_dump
LOGICAL, PRIVATE                :: IOS_debug_no_write
LOGICAL, PRIVATE                :: IOS_debug_no_packing
LOGICAL, PRIVATE                :: IOS_debug_no_subdomaining
LOGICAL, PRIVATE                :: IOS_async_send_null
LOGICAL, PRIVATE                :: IOS_async_stats

INTEGER, PRIVATE                :: listener_thread=-1
INTEGER, PRIVATE                :: writer_thread=-1
INTEGER, PRIVATE                :: reader_thread=-1

PRIVATE ioscntl
NAMELIST / ioscntl /                                                         &
    IOS_tasks_per_server,                                                    &
    IOS_Spacing,                                                             &
    IOS_Interleave,                                                          &
    IOS_Offset,                                                              &
    IOS_buffer_size,                                                         &
    IOS_use_async_stash,                                                     &
    IOS_use_async_dump,                                                      &
    IOS_serialise_mpi_calls,                                                 &
    IOS_thread_0_calls_mpi,                                                  &
    IOS_force_threading_mode,                                                &
    IOS_num_threads,                                                         &
    IOS_debug_no_packing,                                                    &
    IOS_debug_no_write,                                                      &
    IOS_debug_no_subdomaining,                                               &
    IOS_Verbosity,                                                           &
    IOS_backoff_interval,                                                    &
    IOS_timeout,                                                             &
    IOS_acquire_model_prsts,                                                 &
    IOS_local_ro_files,                                                      &
    IOS_concurrency,                                                         &
    IOS_concurrency_max_mem,                                                 &
    IOS_as_concurrency,                                                      &
    IOS_Unit_Alloc_Policy,                                                   &
    IOS_async_levs_per_pack,                                                 &
    IOS_async_send_null,                                                     &
    IOS_async_stats,                                                         &
    IOS_use_helpers,                                                         &
    IOS_RelayToSlaves,                                                       &
    IOS_Decomp_Model,                                                        &
    IOS_print_start_time,                                                    &
    IOS_Lock_Meter,                                                          &
    IOS_Enable_mpiio

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='IOS_INIT'

CONTAINS

LOGICAL FUNCTION IOS_Setup( numIOServers )
USE mpl, ONLY:                                                            &
    mpl_comm_null,                                                         &
    mpl_thread_multiple,                                                   &
    mpl_thread_serialized,                                                 &
    mpl_thread_funneled,                                                   &
    mpl_thread_single,                                                     &
    mpl_integer
USE FilenameLength_mod, ONLY:                                             &
    filenamelength

USE IOS_Decompose, ONLY:                                                  &
    IOS_Decompose_Plan,                                                    &
    IOS_Decompose_Max_Plan

USE umPrintMgr, ONLY:                                                      &
    printerActive,                                                         &
    prnt_writers,                                                          &
    outputZeroIOS

USE um_parcore, ONLY: mype

USE get_env_var_mod, ONLY: get_env_var

IMPLICIT NONE

! Argument
INTEGER, INTENT(IN)    :: numIOServers

INTEGER, PARAMETER     :: atm_coloured=0
INTEGER, PARAMETER     :: io_coloured=1
INTEGER, PARAMETER     :: IOS_Namelist_Unit=127

! Note, formatting of output copes with 1,000,000 mpi ranks.
INTEGER, PARAMETER     :: IOS_Max_Reported=8

INTEGER                :: io_server_counter
INTEGER                :: io_server_rank

INTEGER                :: handler
INTEGER                :: ierror
INTEGER                :: num_threads
LOGICAL                :: isIOSCapable

! Various local vars for arithmetic
INTEGER                :: atm_next_key
INTEGER                :: num_units
INTEGER                :: remaining_ios
INTEGER                :: allocated_ios
INTEGER                :: colour
INTEGER                :: key
INTEGER                :: subkey
INTEGER                :: subcolour
INTEGER                :: i,j
INTEGER                :: error_status
INTEGER                :: length
LOGICAL                :: Flag

REAL                   :: t1

! Strings for IO
CHARACTER (LEN=80)             :: IOS_ini_message
CHARACTER (LEN=80)             :: IOS_capable_string
CHARACTER (LEN=errormessagelength)  :: iomessage
CHARACTER (LEN=FileNameLength) :: NamelistFile
CHARACTER (LEN=*), PARAMETER   :: RoutineName = 'IOS_INIT:IOS_SETUP'


INTEGER                        :: icode

! Get current communicator from GCOM
CALL GC_Get_Communicator(global_comm, ierror)
CALL MPL_Comm_Rank(global_comm,global_rank  ,ierror)
CALL MPL_Comm_Size(global_comm,global_procs ,ierror)

! Set some sensible defaults and initialise variables
! All namelist items not listed here are initialised in IOS_common
IOS_Setup                   = .FALSE.
IOS_force_threading_mode    = imdi
IOS_ReadCompletionRequested = 0
IOS_EnqCompletionRequested  = 0
IOS_spacing                 = imdi
IOS_Offset                  = imdi
IOS_Interleave              = .FALSE.
IOS_use_async_stash         = .FALSE.
IOS_use_async_dump          = .FALSE.
l_io_server                 = .FALSE.
l_io_leader                 = .FALSE.
IOS_debug_no_packing        = .FALSE.
IOS_debug_no_write          = .FALSE.
IOS_debug_no_subdomaining   = .FALSE.
IOS_async_send_null         = .FALSE.
IOS_async_stats             = .FALSE.
IOS_buffer_size             = imdi
IOS_concurrency             = imdi
IOS_concurrency_max_mem     = imdi
IOS_as_concurrency          = imdi
IOS_async_levs_per_pack     = imdi

! Read Namelist for control
CALL get_env_var("IOSCNTL",NameListFile,allow_missing=.TRUE.,length=length)
IF ( length < 0) THEN
  NameListFile='IOSCNTL.namelist'
END IF
WRITE(IOS_message,'(A,A)')'Info: Control file:',TRIM(NameListFile)
CALL IOS_print(IOS_message,src='ios_init')
IF (mype == 0) OPEN(IOS_Namelist_Unit,FILE=NameListFile,FORM='FORMATTED',  &
    STATUS='OLD', ACTION='READ', IOSTAT=ierror, IOMSG=iomessage)
CALL mpl_bcast(ierror,1,mpl_integer,0,global_comm,icode)
IF ( ierror == 0 ) THEN
  IF (IOS_Verbosity > IOS_PrStatus_Oper .AND. mype == 0) THEN
    CALL IOS_print('Info: Reading IOS control file')
  END IF
  CALL read_nml_ioscntl(IOS_Namelist_Unit, ierror)
  IF ( ierror /= 0 .AND. mype == 0) THEN
    CALL IOS_Ereport(RoutineName,123,                                      &
        'Error reading IOS control file (please check setup).')
  END IF

  IF (mype == 0) CLOSE(IOS_Namelist_Unit)
ELSE
  CALL IOS_Ereport(RoutineName,124,                                         &
      'Failed to open IOS control file: '// TRIM(NameListFile) //  newline//&
      'IoMsg: '//TRIM(iomessage))
END IF

! Some checks and resets:
IF (numIOServers <= 0) THEN

  ! No I/O servers
  IOS_spacing = global_procs

ELSE

  IF (IOS_Spacing == 0) THEN
    ! User has requested evenly distributed IO servers
    IOS_Spacing = global_procs/numIOServers
    IOS_Offset  = 0
  ELSE IF (IOS_Spacing < 0) THEN
    ! Spacing not set or invalid input
    WRITE(IOS_ini_message,'(A)')                                           &
        'Spacing of IO servers is negative or undefined'
    CALL IOS_Ereport(RoutineName,10,IOS_ini_message)
  END IF

END IF

#if defined(CPP_IOS_STATIC_DATA)
CALL IOS_print('Overriding sizes with compile time values',src='ios_init')
IOS_concurrency             = cpp_ios_concurrency
IOS_as_concurrency          = cpp_ios_as_concurrency
IOS_async_levs_per_pack     = cpp_ios_as_levels
#endif

! Set namelist flags in other modules
!$  IOS_AsyncNumSlots        = IOS_as_concurrency
!$  IOS_AsyncMaxFieldsInPack = IOS_async_levs_per_pack
!$  IOS_AsyncSendNull        = IOS_async_send_null
!$  IOS_AsyncDoStats         = IOS_async_stats

! Set default word lengths for types
IOS_BytesPerInteger      = umFortranIntegerSize()
IOS_BytesPerReal         = umFortranRealSize()

! Determine capability for IOS of runtime environment
isIOSCapable=.FALSE.
WRITE(IOS_capable_string,'(A)')'Not OpenMP compiled'
!$  isIOSCapable=.TRUE.  ! If openMP compiled
!$  WRITE(IOS_capable_string,'(A)')'IOS is possible'

! Override settings according to MPI's capabilities
CALL mpl_query_thread (threading_model,ierror)

IF (ios_force_threading_mode /= imdi) THEN
  error_status = -10
  CALL Ereport("IOS_INIT", error_status,                          &
       "MPI's advertised threading capability is being overriden.")
END IF

! Check for a namelist override
IF (ios_force_threading_mode==0) THEN
  threading_model=mpl_thread_single
ELSE IF (ios_force_threading_mode==1) THEN
  threading_model=mpl_thread_funneled
ELSE IF (ios_force_threading_mode==2) THEN
  threading_model=mpl_thread_serialized
ELSE IF (ios_force_threading_mode==3) THEN
  threading_model=mpl_thread_multiple
ELSE IF (ios_force_threading_mode/=imdi) THEN
  error_status = 10
  CALL Ereport("IOS_INIT", error_status,                          &
       "Incorrect value set for ios_force_threading_mode.")
END IF ! Any other value retains the result of MPL_QUERY_THREAD()

IF ( threading_model == mpl_thread_multiple ) THEN
  CALL IOS_print('Info: Full Multithreading available.',src='ios_init')
ELSE IF ( threading_model == mpl_thread_serialized ) THEN
  IOS_serialise_mpi_calls=.TRUE.
  CALL IOS_print('Info: Serialized threading available.',src='ios_init')
ELSE IF ( threading_model == mpl_thread_funneled ) THEN
  IOS_thread_0_calls_mpi=.TRUE.
  CALL IOS_print('Info: Funneled threading available.',src='ios_init')
ELSE IF ( threading_model == mpl_thread_single ) THEN
  isIOSCapable=.FALSE.
  WRITE(IOS_capable_string,'(A,A)')                                        &
      'MPI implementation not capable ',                                   &
      'of running threaded jobs.'
  WRITE(IOS_message,'(A,A)')                                               &
      'If you are using OpenMPI or think your MPI is thread capable, ',    &
      'try setting ios_force_threading_mode appropriately in the UI'
  CALL IOS_print(IOS_message,src='ios_init')
ELSE
  isIOSCapable=.FALSE.
  WRITE(IOS_capable_string,'(A)')'Broken MPI cannot describe itself'
  WRITE(IOS_ini_message,'(A,I8)')                                          &
      'MPI reported an unknown threading model:',threading_model
  CALL IOS_Ereport(RoutineName,-99,IOS_ini_message)
END IF

! Saneness for single threaded MPI
IF (IOS_thread_0_calls_mpi  .AND.                                          &
    (IOS_use_async_stash .OR. IOS_use_async_dump)) THEN
  IOS_use_async_stash=.FALSE.
  IOS_use_async_dump=.FALSE.
  WRITE(IOS_message,'(A,A)')'Info: IOS Accelerated dump/STASH disabled:',  &
      ' Require at least MPL_THREAD_SERIALIZED'
  CALL IOS_print(IOS_message,src='ios_init')
  CALL IOS_print('Info: On the IBM parallel environment use' &
                 // ' MP_SINGLE_THREAD=NO to correct this',  &
                 src='ios_init')
END IF

IF ( IOS_thread_0_calls_mpi .AND. IOS_Enable_mpiio ) THEN
   IOS_Enable_mpiio=.FALSE.
   WRITE(IOS_message,'(A,A)')'Info: MPI-IO disabled: requires at least ',  &
        'MPL_THREAD_SERIALIZED'
   CALL IOS_print(IOS_message,src='ios_init')
END IF


! Print out parameters being used
IF ( global_rank == 0 .AND.                                                &
    IOS_Verbosity >= IOS_PrStatus_Oper) THEN
  WRITE(IOS_message,'(A,I5,A)')'Info: Task spacing               = ',      &
      IOS_Spacing, ' tasks'
  CALL IOS_print(IOS_message,src='ios_init')
  WRITE(IOS_message,'(A,I5,A)')'Info: Task offset                = ',      &
      IOS_Offset, ' tasks'
  CALL IOS_print(IOS_message,src='ios_init')
  WRITE(IOS_message,'(A,L1,A)')'Info: Interleaved servers        = ',      &
      IOS_Interleave
  CALL IOS_print(IOS_message,src='ios_init')
  WRITE(IOS_message,'(A,I5,A)')'Info: Server size                = ',      &
      IOS_tasks_per_server, ' tasks'
  CALL IOS_print(IOS_message,src='ios_init')
  WRITE(IOS_message,'(A,I5,A)')'Info: Buffer size                = ',      &
      IOS_buffer_size,' MB'
  CALL IOS_print(IOS_message,src='ios_init')
  WRITE(IOS_message,'(A,I1)')'Info: Verbosity                  = ',        &
      IOS_Verbosity
  CALL IOS_print(IOS_message,src='ios_init')
  WRITE(IOS_message,'(A,L1)')'Info: Asynchronous stash         = ',        &
      IOS_Use_Async_Stash
  CALL IOS_print(IOS_message,src='ios_init')
  WRITE(IOS_message,'(A,L1)')'Info: Asynchronous dumps         = ',        &
      IOS_Use_Async_Dump
  CALL IOS_print(IOS_message,src='ios_init')
  WRITE(IOS_message,'(A,L1)')'Info: Serialise all mpi calls    = ',        &
      IOS_serialise_mpi_calls
  CALL IOS_print(IOS_message,src='ios_init')
  WRITE(IOS_message,'(A,L1)')'Info: Only thread zero calls MPI = ',        &
      IOS_thread_0_calls_mpi
  CALL IOS_print(IOS_message,src='ios_init')
  WRITE(IOS_message,'(A,L1)')'Info: Read only files stay local = ',        &
      IOS_local_ro_files
  CALL IOS_print(IOS_message,src='ios_init')
  WRITE(IOS_message,'(A,I0)')'Info: Unit allocation policy     = ',        &
      IOS_Unit_Alloc_Policy
  CALL IOS_print(IOS_message,src='ios_init')
  WRITE(IOS_message,'(A,I0)')'Info: Polling Interval           = ',        &
      IOS_backoff_interval
  CALL IOS_print(IOS_message,src='ios_init')
  WRITE(IOS_message,'(A,I0)')'Info: Timeout Interval           = ',        &
      IOS_timeout
  CALL IOS_print(IOS_message,src='ios_init')
  WRITE(IOS_message,'(A,L1)')'Info: Lock Metering              = ',        &
      IOS_lock_meter
  CALL IOS_print(IOS_message,src='ios_init')
  WRITE(IOS_message,'(A,L1)')'Info: Relay protocol to slaves   = ',        &
      IOS_RelayToSlaves
  CALL IOS_print(IOS_message,src='ios_init')
  WRITE(IOS_message,'(A,I0)')'Info: Server Decomposition Model = ',        &
      IOS_Decomp_Model
  CALL IOS_print(IOS_message,src='ios_init')
  WRITE(IOS_message,'(A,L1)')'Info: Acquire model output level = ',        &
      IOS_acquire_model_prsts
  WRITE(IOS_message,'(A,L1)')'Info: Info: Use helper threads   = ',        &
      IOS_use_helpers
  CALL IOS_print(IOS_message,src='ios_init')
  WRITE(IOS_message,'(A,L1)')'Info: Info: Enable MPI-IO        = ',        &
      IOS_enable_mpiio
  CALL IOS_print(IOS_message,src='ios_init')
  WRITE(IOS_message,'(A,I0)')'Info: Async Dispatch slots       = ',        &
      IOS_Concurrency
  CALL IOS_print(IOS_message,src='ios_init')
  WRITE(IOS_message,'(A,I0,A)')'Info: Dispatch slots mem cap     = ',      &
      IOS_Concurrency_max_mem,' MB'
  CALL IOS_print(IOS_message,src='ios_init')
!$    WRITE(IOS_message,'(A,I3)')'Info: Async Stash Dispatch slots = ',        &
!$        IOS_AsyncNumSlots
!$    CALL IOS_print(IOS_message,src='ios_init')
!$    WRITE(IOS_message,'(A,I3)')'Info: Async fields/levs per pack = ',        &
!$        IOS_AsyncMaxFieldsInPack
!$    CALL IOS_print(IOS_message,src='ios_init')
!$    WRITE(IOS_message,'(A,L1)')'Info: Async send empty tiles     = ',        &
!$        IOS_AsyncSendNull
!$    CALL IOS_print(IOS_message,src='ios_init')
!$    WRITE(IOS_message,'(A,L1)')'Info: Async stats profiling      = ',        &
!$        IOS_AsyncDoStats
!$    CALL IOS_print(IOS_message,src='ios_init')
  WRITE(IOS_message,'(A,L1)')'Info: IOS print start time       = ',        &
      IOS_print_start_time
  CALL IOS_print(IOS_message,src='ios_init')
  IF (IOS_debug_no_packing)                                                &
      CALL IOS_print('Info: Debug Option: No stash packing.',src='ios_init')
  IF (IOS_debug_no_subdomaining)                                           &
      CALL IOS_print                                                       &
      ('Info: Debug Option: No stash subdomaining.',src='ios_init')
  IF (IOS_debug_no_write)                                                  &
      CALL IOS_print('Info: Debug Option: No data writes.',                &
      src='ios_init')
END IF

! Convert this to bytes from human friendly MB provided in the namelist
IOS_Concurrency_max_mem=IOS_Concurrency_max_mem*1024*1024

! Set the IOS team decomposition strategy
IF (IOS_Decomp_Model > IOS_Decompose_Max_Plan .OR.                         &
    IOS_Decomp_Model < 0) THEN
  CALL IOS_Ereport(RoutineName,-99,                                        &
      "Invalid decomposition: IOS_Decomp_Model, using default")
  IOS_Decomp_Model=0
END IF

IOS_Decompose_Plan=IOS_Decomp_Model

IF (IOS_debug_no_packing .OR. IOS_debug_no_subdomaining .OR.               &
    IOS_debug_no_write) THEN
  CALL IOS_Ereport(RoutineName,-99,                                        &
      "One or more debug options are set, output may be invalid")
END IF

IF (IOS_tasks_per_server /= 0) THEN
  IOS_server_groups=numIOServers/IOS_tasks_per_server
ELSE
  IF ( numIOServers<=0 ) THEN
    WRITE(ios_ini_message,'(A)')                                           &
      'No IO servers but resetting IOS_tasks_per_server to be greater '//  & 
      'than zero.'
    CALL IOS_Ereport(RoutineName,-98,IOS_ini_message)
    IOS_tasks_per_server=1
    IOS_server_groups=numIOServers/IOS_tasks_per_server
  ELSE
    WRITE(IOS_message,'(A)')                                   &
      'IF you wish to turn off IO servers set FLUME_IOS_NPROC to zero '//  &
      'and reset IOS_tasks_per_server to be greater than zero.'
    CALL IOS_print(IOS_message,src='ios_init')
    WRITE(ios_ini_message,'(A)')                                           &
      'The number of IOS_tasks_per_server must be greater than zero.'
    CALL IOS_Ereport(RoutineName,98,IOS_ini_message)
  END IF
END IF

IF (IOS_server_groups*IOS_tasks_per_server                                 &
    /= numIOServers) THEN
  WRITE(ios_ini_message,'(A,I0,A,I0,A)')                                   &
      'The number of IO tasks (',numIOServers,                             &
      ') must be a multiple of server size (',IOS_tasks_per_server,')'
  CALL IOS_Ereport(RoutineName,99,IOS_ini_message)
END IF

! Allocate Storage
ALLOCATE(io_servers                                                        &
    (IOS_server_groups,IOS_tasks_per_Server))

! Allocate global processor numbers to IOS or ATM models
atm_next_key      = 0
io_server_counter = 0
io_server_rank    = 0
allocated_ios     = 0

! Loop over cpus and allocate them to a role
DO i=0,global_procs-1

  remaining_ios=numIOServers-allocated_ios
  IF (( i-IOS_Offset >= 0                            .AND.                 &
      MOD(i-IOS_Offset,IOS_Spacing) == IOS_Spacing-1 .AND.                 &
      allocated_ios    < numIOServers )              .OR.                  &
      global_procs-i-1 < remaining_ios ) THEN

    IF ( global_rank == i ) THEN
      colour = io_coloured+io_server_counter
      key    = io_server_rank
    END IF

    io_servers(io_server_counter+1,io_server_rank+1)=i
    allocated_ios = allocated_ios+1

    IF (IOS_Interleave) THEN
      io_server_counter=io_server_counter+1
      IF (io_server_counter == IOS_Server_groups ) THEN
        io_server_counter = 0
        io_server_rank    = io_server_rank+1
      END IF
    ELSE
      io_server_rank = io_server_rank+1
      IF (io_server_rank == IOS_tasks_per_server) THEN
        io_server_rank    = 0
        io_server_counter = io_server_counter+1
      END IF
    END IF

  ELSE ! not an io server
    IF ( global_rank == i ) THEN
      colour     = atm_coloured
      key        = atm_next_key
    END IF
    atm_next_key = atm_next_key+1
  END IF
END DO

IF ( numIOServers > 0 ) THEN
  CALL IOS_print(' ',src='ios_init')
  WRITE(ios_ini_message,*)                                                 &
      '(A,I5,A',(',I6',j=1,MIN(IOS_Max_Reported,IOS_tasks_per_server)),')'
  DO i=1,IOS_Server_groups
    ! This string is 30 + 6* ios_tasks_per_server
    ! so it can only support (30 is the number of literal chars)
    ! (80-30)/6 = 8 = IOS_Max_Reported servers
    WRITE(IOS_message,ios_ini_message)'Info: IO Server ',i,' is ',         &
        (io_servers(i,j),j=1,MIN(IOS_Max_Reported,IOS_tasks_per_server))
    CALL IOS_print(IOS_message,src='ios_init')
    IF (IOS_tasks_per_server > IOS_Max_Reported) THEN
      WRITE(IOS_message,'(A,I5,A)')'         :    and an additional',      &
          IOS_tasks_per_server-IOS_Max_Reported,' other processes '
      CALL IOS_print(IOS_message,src='ios_init')
    END IF
    IF (prnt_writers==outputZeroIOS) THEN
      IF (mype==io_servers(i,1)) printerActive=.TRUE.
    END IF
  END DO
ELSE
  CALL IOS_print('Info: IO servers are not configured',src='ios_init')
END IF

!
! Set up all the communicators we want....
!
!
! Split the communicator into IO and ATM tasks
subcolour=atm_coloured
subkey=key
IF (colour>atm_coloured) THEN
  subcolour = io_coloured
  subkey    = key+colour*IOS_tasks_per_server
END IF
CALL MPL_Comm_Split(global_comm, subcolour, subkey, io_comm, ierror)
CALL MPL_Comm_Rank (io_comm, io_rank ,ierror)
CALL MPL_Comm_Size (io_comm, io_procs ,ierror)

! IO communicator is useless on atmos ranks so zap it.
IF (subcolour == atm_coloured) THEN
  io_comm  = mpl_comm_null
  io_procs = numIOServers
  io_rank  = -1
END IF

! Split the communicator into leaders and non-leaders
IF (key==0) THEN
  subcolour=0
  subkey=colour
ELSE
  subcolour=1
  subkey=0
END IF
CALL MPL_Comm_Split(global_comm, subcolour, subkey, leader_comm,           &
    ierror)
CALL MPL_Comm_Rank (leader_comm, leader_rank ,ierror)
CALL MPL_Comm_Size (leader_comm, leader_procs ,ierror)
! leader communicator is useless on non-leader ranks so zap it.
IF (subcolour == 1) THEN
  leader_comm  = mpl_comm_null
  leader_procs = IOS_Server_groups+1
  leader_rank  = -1
END IF

! Split the communicator into the main IO Groups and ATM tasks
CALL MPL_Comm_Split(global_comm, colour, key, model_comm, ierror)
CALL MPL_Comm_Rank (model_comm , model_rank  ,ierror)
CALL MPL_Comm_Size (model_comm , model_procs ,ierror)

! Lets keep funky communications in a funky communicator....
!$  CALL MPL_Comm_dup(global_comm,IOS_async_comm,ierror)

! Set up communicators for broadcasts....
IF ( colour > atm_coloured ) THEN
  l_io_server = .TRUE.
  IF (leader_comm /= mpl_comm_null) THEN
    l_io_leader = .TRUE.
  END IF
END IF

ALLOCATE(IOS_BCast_Comm(numIOServers))
! Generate an MPI group suitable for broadcasting from an IO server

! set the server comm to something sane in case of inadvertant use.
IOS_BCast_Server_Comm=mpl_comm_null

! There is a broadcast comm for each io server group.
DO i=1,IOS_Server_groups
  IF (l_io_leader .AND. colour==i) THEN
    subcolour=1
    subkey=global_procs-io_procs
  ELSE IF (l_io_server) THEN
    subcolour=2
    subkey=model_rank
  ELSE
    subcolour=1
    subkey=model_rank
  END IF

  CALL MPL_Comm_Split(global_comm, subcolour,                              &
      subkey, IOS_bcast_comm(i), ierror)
  IF (subcolour /= 2) THEN
    CALL MPL_Comm_Size( IOS_bcast_comm(i), bcast_procs ,ierror)
    CALL MPL_Comm_Rank( IOS_bcast_comm(i), bcast_rank  ,ierror)
  END IF

  IF (l_io_leader .AND. colour == i) THEN
    IOS_BCast_Server_Comm=IOS_bcast_comm(i)
  END IF

END DO






! Everyone synchronise their watches...

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(t1)
!$OMP BARRIER
  t1 = get_wallclock_time()
!$OMP CRITICAL(internal_write)
  WRITE(IOS_message,'(A,F10.3,A )') 'IO servers starting at time= ',  &
        t1 - IOS_Start_Time, ' seconds'
!$OMP END CRITICAL(internal_write)

 CALL IOS_print(IOS_message,src='ios_init')
!$OMP END PARALLEL 

! Set the default communicator for GCOM
CALL GC_Set_Communicator                                                   &
    (model_comm, model_rank, model_procs, ierror)

IF (IOS_Verbosity >= IOS_PrStatus_Oper) THEN

  WRITE(IOS_message,'(A,I7,A,I7,A,I7)')                                    &
      'Info: Original Rank=',global_rank,' Global_Rank=',                  &
      global_rank,' Model_Rank=',model_rank
  CALL IOS_print(IOS_message,src='ios_init')
  WRITE(IOS_message,'(A,I7,A,I7)')                                         &
      'Info: Total Size=',global_procs,' Submodel size=',                  &
      model_procs
  CALL IOS_print(IOS_message,src='ios_init')

  CALL IOS_print(' ',src='ios_init')
  CALL IOS_print('Communicators....',src='ios_init')
  CALL IOS_print(' ',src='ios_init')
  WRITE(IOS_message,'(A12,A12,A12,A12)')'Function',                        &
      'Comm ID','Processors','My Rank'
  CALL IOS_print(IOS_message,src='ios_init')
  WRITE(IOS_message,'(A12,A12,A12,A12)')'--------',                        &
      '-------','----------','-------'
  CALL IOS_print(IOS_message,src='ios_init')
  WRITE(IOS_message,'(A8,I12,I12,I12)')'Global  ',                         &
      global_comm,global_procs,global_rank
  CALL IOS_print(IOS_message,src='ios_init')
  WRITE(IOS_message,'(A8,I12,I12,I12)')'Model   ',                         &
      model_comm,model_procs,model_rank
  CALL IOS_print(IOS_message,src='ios_init')
  WRITE(IOS_message,'(A8,I12,I12,I12)')'IO      ',                         &
      io_comm,io_procs,io_rank
  CALL IOS_print(IOS_message,src='ios_init')
  WRITE(IOS_message,'(A8,I12,I12,I12)')'Leader  ',                         &
      leader_comm,leader_procs,leader_rank
  CALL IOS_print(IOS_message,src='ios_init')
  CALL IOS_print(' ',src='ios_init')

END IF

! The primary return of the function
IOS_Setup = l_io_server

! Set printing according to settings
CALL IOS_PrintSetLevel(    &
    level=IOS_Verbosity,   &                ! provide our current level
    setParent=IOS_Setup,   &                ! IO servers set UM print level
    setFromParent=IOS_acquire_model_prsts)  ! unless told to use the UM val.

! register error handlers:

CALL MPL_Comm_Create_ErrHandler(model_mpi_error_handler,                   &
    handler,ierror)
CALL MPL_Comm_Set_ErrHandler(model_comm,handler,ierror)

!$  CALL MPL_Comm_Create_ErrHandler                                            &
!$      (async_mpi_error_handler,handler,ierror)
!$  CALL MPL_Comm_Set_ErrHandler(IOS_async_comm,handler,ierror)

CALL MPL_Comm_Create_ErrHandler(global_mpi_error_handler,                  &
    handler,ierror)
CALL MPL_Comm_Set_ErrHandler(global_comm,handler,ierror)

IF ( colour > atm_coloured ) THEN
  CALL MPL_Comm_Create_ErrHandler(io_mpi_error_handler,                    &
      handler,ierror)
  CALL MPL_Comm_Set_ErrHandler(io_comm,handler,ierror)
END IF

IF ( leader_rank >= 0 ) THEN
  CALL MPL_Comm_Create_ErrHandler(leader_mpi_error_handler,                &
      handler,ierror)
  CALL MPL_Comm_Set_ErrHandler(leader_comm,handler,ierror)
END IF

! Check for no servers defined
IF ( numIOServers <= 0 ) THEN
  isIOSCapable=.FALSE.
  WRITE(IOS_capable_string,'(A)')'No IO server processes assigned'
END IF
! Deactivate if insufficient threads
num_threads=2 ! allow for the 'not openmp message' to propogate

! Set number of threads from the namelist if appropriate
! Note only do this on I/O servers
!$ IF (ios_num_threads > 1 .AND. l_io_server)                           &
!$   CALL omp_set_num_threads(ios_num_threads)

!$OMP PARALLEL DEFAULT(NONE) SHARED(num_threads)
!$OMP MASTER
!$  num_threads=omp_get_num_threads()
!$OMP END MASTER
!$OMP END PARALLEL

! If we've changed the number of threads for the server, write
! out the number of threads available for the server
IF (global_rank == 0 .AND. ios_num_threads > 1 .AND.                    &
    IOS_Verbosity >= IOS_PrStatus_Oper) THEN
  WRITE(IOS_message,'(A,I3)') 'Info: Number of threads for I/O server = ',  &
          ios_num_threads
  CALL IOS_print(IOS_message,src='ios_init')
END IF

! If either normal threading count or the ios changed threading count
! are too low, then we aren't capable of running I/O server
IF ( num_threads < 2 .AND. ios_num_threads < 2 .AND. isIOSCapable) THEN
  isIOSCapable=.FALSE.
  WRITE(IOS_capable_string,'(A,I3)')'Not enough OpenMP Threads:',          &
      num_threads
END IF

! At present only 3 threads do anything, so we'll space them out
! over available threads...
Listener_thread =0
Writer_thread   =num_threads/2

IF (num_threads>2 .AND. IOS_use_helpers)                                   &
    reader_thread = num_threads-1

IF ( isIOSCapable ) THEN
  ! The following code depends on modules only available in OpenMP
  ! compilation so it should not compile on non-threaded builds (!$)
  !
  ! --- Define the unit mapping - which IOS pes deal with which units ---
  !
!$    IF (IOS_Unit_Alloc_Policy /= IOS_Unit_Alloc_Static) THEN
!$      io_server_for_unit_lookup(:) = IOS_No_Server  
!$    ELSE
!$      DO i = minUnit, maxUnit
!$        IF      (i < minUnit ) THEN
!$          io_server_for_unit_lookup(i) = io_servers(1,1)
!$        ELSE IF ( i > maxUnit ) THEN
!$          io_server_for_unit_lookup(i) = io_servers(IOS_Server_groups,1)
!$        ELSE
!$          io_server_for_unit_lookup(i) =                                     &
!$              io_servers(1+MOD(i - minUnit, IOS_Server_groups),1)
!$        END IF
!$      END DO
!$    END IF
!$
!$    CALL IOS_Comms_init()
!$
!$    IF ( l_io_server ) THEN
!$      WRITE(IOS_message,'(A,I7,A)')                                          &
!$          'Info: PE ', global_rank, ' is an I/O server'
!$      CALL IOS_print(IOS_message,src='ios_init')
!$      num_units=0
!$      IF ( l_io_leader ) THEN
!$        DO i = minUnit, maxUnit
!$          IF ( io_server_for_unit(i) == global_rank ) THEN
!$            num_units=num_units+1
!$          END IF
!$        END DO
!$      
!$        WRITE(IOS_message,'(A,I3,A)')                                        &
!$            'Info: I am responsible for ',num_units,' units'
!$        CALL IOS_print(IOS_message,src='ios_init')
!$      END IF
  !
  ! --- Initialise IOS queue structure ---
  !
!$      CALL IOS_Queue_Initialise(IOS_buffer_size)
  !
  ! --- Initialise IOS Stash support (server side) ---
  ! ---  this is incomplete initialisation, geometry will
  ! ---  need initialising later after the model is loaded
  !
!$      CALL IOS_stash_server_init(                                            &
!$          IOS_use_async_stash,                                               &
!$          IOS_use_async_dump,                                                &
!$          IOS_debug_no_packing,                                              &
!$          IOS_debug_no_write,                                                &
!$          IOS_debug_no_subdomaining)
!$    ELSE
  !
  ! --- Initialise IOS Stash support (client side) ---
  !
!$      CALL IOS_client_init()
!$      CALL IOS_stash_client_init(IOS_use_async_stash,                        &
!$          IOS_use_async_dump,IOS_debug_no_subdomaining)
!$    END IF

ELSE
  WRITE(IOS_message,'(A,A)')'Info: Configuration not IOS capable,',        &
      ' deactivating IOS'
  CALL IOS_print(IOS_message,src='ios_init')
  WRITE(IOS_message,'(A,A)')'Info: Reason:',TRIM(IOS_capable_string)
  CALL IOS_print(IOS_message,src='ios_init')

  DEALLOCATE(io_servers)
  NULLIFY(io_servers)
  DEALLOCATE(IOS_BCast_Comm)
  NULLIFY(IOS_BCast_Comm)
END IF

RETURN

END FUNCTION IOS_Setup

SUBROUTINE ios_run()
USE errormessagelength_mod, ONLY: errormessagelength

USE fort2c_portio_interfaces, ONLY: portiodetachallhelpers, portioaddhelper

IMPLICIT NONE
CHARACTER (LEN=*), PARAMETER   :: RoutineName = 'IOS_INIT:IOS_RUN'
CHARACTER (LEN=errormessagelength)             :: message
INTEGER :: warnCode
INTEGER :: this_thread

IF (.NOT. L_IO_Server) THEN
  WRITE(message,'(A)')'Task is not an IO Server, so should not call me'
  warnCode=-1
  CALL IOS_Ereport(RoutineName,99,message)
ELSE IF (.NOT. L_IOS_Active() ) THEN
  WRITE(message,'(A)')'IO Server is not active'
  warnCode=-2
  CALL IOS_Ereport(RoutineName,99,message)
ELSE

  ! --- Allocate subtasks to server threads ---
  !
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(this_thread)
!$    this_thread=omp_get_thread_num()
!$    IF      ( this_thread == listener_thread ) THEN
!$      IF (IOS_RelayToSlaves.AND.model_rank>0) THEN
!$        CALL IOS_SlaveListener()
!$      ELSE
!$        CALL IOS_listener()
!$      END IF
!$    ELSE IF ( this_thread == writer_thread ) THEN
!$      CALL IOS_writer()
!$      CALL IOS_print('Info: Writer completed, detatching any helpers',       &
!$          src='ios_init')
!$      CALL portioDetachAllHelpers()
!$      CALL IOS_print('Info: Writer completed, detatch request completed',    &
!$          src='ios_init')
!$    ELSE IF ( this_thread == reader_thread ) THEN
!$      CALL IOS_print('Info: Reader thread helper attaching',src='ios_init')
!$      CALL portioAddHelper(1)
!$      CALL IOS_print('Info: Reader thread has detatched',src='ios_init')
!$    ELSE
!$      WRITE(IOS_message,'(A,I2,A)')'Info: Thread ',this_thread,              &
!$          ' is inactive in IOS'
!$      CALL IOS_print(IOS_message,src='ios_init')
!$    END IF
!$OMP END PARALLEL
  !
  ! --- Tidy up IOS Stash support (server side) ---
  !
!$    CALL IOS_stash_server_fini()
!$    CALL IOS_Comms_fini()
  !
  ! --- Close IOS queue structure ---
  !
!$    CALL IOS_Queue_Finalise()


END IF

END SUBROUTINE ios_run

SUBROUTINE read_nml_ioscntl(unit_in, ErrorStatus)

USE um_parcore, ONLY: mype

USE setup_namelist, ONLY: setup_nml_type

IMPLICIT NONE

INTEGER, INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER, INTENT(OUT) :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_IOSCNTL'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 2
INTEGER, PARAMETER :: n_int = 15
INTEGER, PARAMETER :: n_log = 17

TYPE my_namelist
  SEQUENCE
  INTEGER :: IOS_tasks_per_server
  INTEGER ::  IOS_Spacing
  INTEGER ::  IOS_Offset
  INTEGER ::  IOS_buffer_size
  INTEGER ::  IOS_force_threading_mode
  INTEGER ::  IOS_Verbosity
  INTEGER ::  IOS_backoff_interval
  INTEGER ::  IOS_timeout
  INTEGER ::  IOS_concurrency
  INTEGER ::  IOS_concurrency_max_mem
  INTEGER ::  IOS_as_concurrency
  INTEGER ::  IOS_Unit_Alloc_Policy
  INTEGER ::  IOS_async_levs_per_pack
  INTEGER ::  IOS_Decomp_Model
  INTEGER ::  IOS_num_threads
  LOGICAL ::  IOS_Interleave
  LOGICAL ::  IOS_use_async_stash
  LOGICAL ::  IOS_use_async_dump
  LOGICAL ::  IOS_serialise_mpi_calls
  LOGICAL ::  IOS_thread_0_calls_mpi
  LOGICAL ::  IOS_debug_no_packing
  LOGICAL ::  IOS_debug_no_write
  LOGICAL ::  IOS_debug_no_subdomaining
  LOGICAL ::  IOS_acquire_model_prsts
  LOGICAL ::  IOS_local_ro_files
  LOGICAL ::  IOS_async_send_null
  LOGICAL ::  IOS_async_stats
  LOGICAL ::  IOS_use_helpers
  LOGICAL ::  IOS_RelayToSlaves
  LOGICAL ::  IOS_print_start_time
  LOGICAL ::  IOS_Lock_Meter
  LOGICAL ::  IOS_Enable_mpiio
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,   &
                    n_log_in=n_log)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=ioscntl, IOSTAT=ErrorStatus, IOMSG=iomessage)

  my_nml % IOS_tasks_per_server     = IOS_tasks_per_server
  my_nml % IOS_Spacing              = IOS_Spacing
  my_nml % IOS_Offset               = IOS_Offset
  my_nml % IOS_buffer_size          = IOS_buffer_size
  my_nml % IOS_force_threading_mode = IOS_force_threading_mode
  my_nml % IOS_Verbosity            = IOS_Verbosity
  my_nml % IOS_backoff_interval     = IOS_backoff_interval
  my_nml % IOS_timeout              = IOS_timeout
  my_nml % IOS_concurrency          = IOS_concurrency
  my_nml % IOS_concurrency_max_mem  = IOS_concurrency_max_mem
  my_nml % IOS_as_concurrency       = IOS_as_concurrency
  my_nml % IOS_Unit_Alloc_Policy    = IOS_Unit_Alloc_Policy
  my_nml % IOS_async_levs_per_pack  = IOS_async_levs_per_pack
  my_nml % IOS_Decomp_Model         = IOS_Decomp_Model
  my_nml % IOS_num_threads          = IOS_num_threads
  ! end of integers
  my_nml % IOS_Interleave            = IOS_Interleave
  my_nml % IOS_use_async_stash       = IOS_use_async_stash
  my_nml % IOS_use_async_dump        = IOS_use_async_dump
  my_nml % IOS_serialise_mpi_calls   = IOS_serialise_mpi_calls
  my_nml % IOS_thread_0_calls_mpi    = IOS_thread_0_calls_mpi
  my_nml % IOS_debug_no_packing      = IOS_debug_no_packing
  my_nml % IOS_debug_no_write        = IOS_debug_no_write
  my_nml % IOS_debug_no_subdomaining = IOS_debug_no_subdomaining
  my_nml % IOS_acquire_model_prsts   = IOS_acquire_model_prsts
  my_nml % IOS_local_ro_files        = IOS_local_ro_files
  my_nml % IOS_async_send_null       = IOS_async_send_null
  my_nml % IOS_async_stats           = IOS_async_stats
  my_nml % IOS_use_helpers           = IOS_use_helpers
  my_nml % IOS_RelayToSlaves         = IOS_RelayToSlaves
  my_nml % IOS_print_start_time      = IOS_print_start_time
  my_nml % IOS_Lock_Meter            = IOS_Lock_Meter
  my_nml % IOS_Enable_mpiio          = IOS_Enable_mpiio

END IF

CALL mpl_bcast(my_nml,1,MPL_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  IOS_tasks_per_server     = my_nml % IOS_tasks_per_server
  IOS_Spacing              = my_nml % IOS_Spacing
  IOS_Offset               = my_nml % IOS_Offset
  IOS_buffer_size          = my_nml % IOS_buffer_size
  IOS_force_threading_mode = my_nml % IOS_force_threading_mode
  IOS_Verbosity            = my_nml % IOS_Verbosity
  IOS_backoff_interval     = my_nml % IOS_backoff_interval
  IOS_timeout              = my_nml % IOS_timeout
  IOS_concurrency          = my_nml % IOS_concurrency
  IOS_concurrency_max_mem  = my_nml % IOS_concurrency_max_mem
  IOS_as_concurrency       = my_nml % IOS_as_concurrency
  IOS_Unit_Alloc_Policy    = my_nml % IOS_Unit_Alloc_Policy
  IOS_async_levs_per_pack  = my_nml % IOS_async_levs_per_pack
  IOS_Decomp_Model         = my_nml % IOS_Decomp_Model
  IOS_num_threads          = my_nml % IOS_num_threads
  ! end of integers
  IOS_Interleave            = my_nml % IOS_Interleave
  IOS_use_async_stash       = my_nml % IOS_use_async_stash
  IOS_use_async_dump        = my_nml % IOS_use_async_dump
  IOS_serialise_mpi_calls   = my_nml % IOS_serialise_mpi_calls
  IOS_thread_0_calls_mpi    = my_nml % IOS_thread_0_calls_mpi
  IOS_debug_no_packing      = my_nml % IOS_debug_no_packing
  IOS_debug_no_write        = my_nml % IOS_debug_no_write
  IOS_debug_no_subdomaining = my_nml % IOS_debug_no_subdomaining
  IOS_acquire_model_prsts   = my_nml % IOS_acquire_model_prsts
  IOS_local_ro_files        = my_nml % IOS_local_ro_files
  IOS_async_send_null       = my_nml % IOS_async_send_null
  IOS_async_stats           = my_nml % IOS_async_stats
  IOS_use_helpers           = my_nml % IOS_use_helpers
  IOS_RelayToSlaves         = my_nml % IOS_RelayToSlaves
  IOS_print_start_time      = my_nml % IOS_print_start_time
  IOS_Lock_Meter            = my_nml % IOS_Lock_Meter
  IOS_Enable_mpiio          = my_nml % IOS_Enable_mpiio

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_ioscntl

END MODULE IOS_Init
