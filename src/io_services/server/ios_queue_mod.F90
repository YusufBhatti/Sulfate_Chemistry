! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: C96

MODULE IOS_Queue_Mod
USE um_types
USE IOS_Constants, ONLY:                                                     &
    ios_md_len, IOS_BytesPerWord32, IOS_BytesPerWord64, IOS_Action_Unset
USE IOS_types, ONLY:                                                         &
    ios_metadata_type
USE IOS_Common, ONLY:                                                        &
    serialize_all_ops, IOS_thread_0_calls_mpi, IOS_ActionName,               &
    IOS_start_time, IOS_backoff_interval, IOS_BytesPerInteger, IOS_ereport,  &
    IOS_Lock_meter, no_lock_needed
USE IOS_communicators, ONLY:                                                 &
    global_rank
!$ USE omp_lib
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
    IOS_PrStatus_Oper, IOS_PrStatus_Diag

USE errormessagelength_mod, ONLY: errormessagelength

USE get_wallclock_time_mod, ONLY: get_wallclock_time
 
IMPLICIT NONE


!-----------------------------------------------------------------------
! TYPE definitions
!-----------------------------------------------------------------------
TYPE IOS_node_type
   ! Queue related
  INTEGER                           :: transaction_number
  INTEGER                           :: payloadTag
  TYPE(IOS_metadata_type)           :: metadata
  TYPE(IOS_node_type), POINTER      :: next
  TYPE(IOS_node_type), POINTER      :: prev

   ! Payload related:
  INTEGER, POINTER                  :: integer_data(:)
  INTEGER(KIND=integer32), POINTER  :: integer32_data(:)
  INTEGER, POINTER                  :: receive_requests(:)
  INTEGER, POINTER                  :: receive_tags(:)
  INTEGER, POINTER                  :: receive_data_len(:)
  REAL, POINTER                     :: distributed_data(:,:)
  REAL(KIND=real32), POINTER        :: real32_data(:)
  REAL, POINTER                     :: real_data(:)
END TYPE IOS_node_type

!-----------------------------------------------------------------------
! Common variables
!-----------------------------------------------------------------------

INTEGER, PRIVATE, PARAMETER   :: logunit                    = 9
INTEGER, PRIVATE              :: IOS_transaction_number     = 0
INTEGER, PRIVATE              :: IOS_transactions_completed = 0
INTEGER, PRIVATE              :: IOS_buffer_max
INTEGER, PRIVATE, volatile    :: Q_size                     = 0
INTEGER, PRIVATE, volatile    :: Q_items                    = 0
INTEGER                       :: IOS_AsyncCompletionRequested
! Only needed by OpenMP builds
!$  INTEGER (KIND=omp_lock_kind)  :: Q_lock_var
TYPE(IOS_node_type), POINTER  :: Q_first
TYPE(IOS_node_type), POINTER  :: Q_last
CHARACTER (LEN=errormessagelength)           :: qmessage
CHARACTER (LEN=*),                                                           &
    PARAMETER, PRIVATE        :: RoutineName = 'IOS_QUEUE'

! Lock metering
INTEGER, PARAMETER, PRIVATE   :: words_per_cacheline        = 64
! Only needed by OpenMP builds
!$  INTEGER, POINTER, PRIVATE     :: LockTries(:,:)
!$  REAL, POINTER, PRIVATE        :: LockTime (:,:)


INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in          = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out         = 1

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='IOS_QUEUE_MOD'

CONTAINS

!-----------------------------------------------------------------------
! SUBROUTINE to put node into Queue (node allocated outside)
!-----------------------------------------------------------------------
SUBROUTINE IOS_Put_Node_In_Queue(node)

USE IOS_Constants, ONLY: IOS_Action_Strlen
IMPLICIT NONE

CHARACTER (LEN=IOS_Action_Strlen) :: action_name
TYPE(IOS_node_type), POINTER :: node
INTEGER                      :: errCode
LOGICAL                      :: dummy

REAL(KIND=jprb):: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='IOS_PUT_NODE_IN_QUEUE'
IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( .NOT. ASSOCIATED(node) ) THEN
  WRITE(qmessage,'(A)')'Tried to put an unassociated POINTER on the Q'
  CALL IOS_Ereport( RoutineName,99, Qmessage )
END IF

node%transaction_number=IOS_transaction_number
IOS_transaction_number=IOS_transaction_number+1

IF ( serialize_all_ops ) THEN
  IF ( IOS_thread_0_calls_mpi ) THEN
    WRITE(IOS_message,'(A,A)')                                             &
        'IOS_queue: Warning: Cant serialize operations',                   &
        ' in this MPI mode'
    CALL IOS_print(IOS_message,src='ios_queue_mod')
  ELSE
    dummy=WaitForDrain()
  END IF
END IF

CALL acquire_lock()
CALL IOS_Increment_Queue_Length(ios_md_len*umFortranIntegerSize(),          &
                                no_lock_needed)
!$OMP FLUSH
IF ( .NOT. ASSOCIATED( Q_first) ) THEN
  Q_first => node
!$OMP FLUSH
END IF

IF ( ASSOCIATED( Q_last ) ) THEN
  Q_last%next => node
  node%prev => Q_last
  NULLIFY (node%next)
END IF

Q_last => node
!$OMP FLUSH
Q_items=Q_items+1
!$OMP FLUSH

IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
  action_name=TRIM(IOS_ActionName(node%metadata%ACTION))
!$OMP CRITICAL(internal_write)
  WRITE(qmessage,'(A,A,A,I0,A,I0,A)')                                      &
      'Info: Queue: Added action ',                                        &
      action_name,                                                         &
      ' trns_no: ',node%transaction_number ,                               &
      ' to queue, now ',Q_items,' items'
!$OMP END CRITICAL(internal_write)
  CALL IOS_print(TRIM(qmessage),src='ios_queue_mod')
  CALL IOS_print_flush()
END IF

CALL release_lock()

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE IOS_Put_Node_In_Queue

!-----------------------------------------------------------------------
! SUBROUTINE to put metadata into Queue
!-----------------------------------------------------------------------
SUBROUTINE IOS_Put_Metadata_In_Queue(metadata)
IMPLICIT NONE

TYPE(IOS_metadata_type), INTENT(IN) :: metadata
TYPE(IOS_node_type), POINTER        :: node
CHARACTER(LEN=80) :: mes
node => make_new_node()
node%metadata = metadata
CALL IOS_Put_Node_In_Queue(node)
END SUBROUTINE IOS_Put_Metadata_In_Queue

!-----------------------------------------------------------------------
! SUBROUTINE to remove last request and data from Queue
!-----------------------------------------------------------------------
SUBROUTINE IOS_Remove_Data_From_Queue()
USE IOS_Constants, ONLY: IOS_Action_Strlen
IMPLICIT NONE

TYPE(IOS_node_type), POINTER         :: node
INTEGER                              :: errCode
REAL(KIND=jprb)                      :: zhook_handle

CHARACTER (LEN=IOS_Action_Strlen) :: action_name
CHARACTER(LEN=*), PARAMETER :: RoutineName='IOS_REMOVE_DATA_FROM_QUEUE'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL acquire_lock()
!$OMP FLUSH
IF ( ASSOCIATED (Q_first) ) THEN
  node => Q_first
  IF ( ASSOCIATED (Q_first%next ) ) THEN
    Q_first => Q_first%next
  ELSE! There is only 1 thing in the queue and we are
      ! about to delete it
    NULLIFY(Q_first)
    NULLIFY(Q_last)
  END IF
!$OMP FLUSH
  Q_items=Q_items-1
!$OMP FLUSH

  IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
    action_name=TRIM(IOS_ActionName(node%metadata%ACTION))
    IF ( Q_items == 0 ) THEN
!$OMP CRITICAL(internal_write)
      WRITE(qmessage,'(A,A,A,I0,A)')                                       &
          'Info: Queue: Removing ',                                        &
          action_name,                                                     &
          ' Trns no ',node%transaction_number,                             &
          ' from queue, queue is empty'
!$OMP END CRITICAL(internal_write)
    ELSE
!$OMP CRITICAL(internal_write)
      WRITE(qmessage,'(A,A,A,I0,A,I0,A)')                                  &
          'Info: Queue: Removing ',                                        &
          action_name,                                                     &
          ' Trns no. ',node%transaction_number,                            &
          ' from queue, now ',Q_items,' items'
!$OMP END CRITICAL(internal_write)
    END IF
    CALL IOS_print(TRIM(qmessage),src='ios_queue_mod')
    CALL IOS_print_flush()
  END IF
  CALL destroy_node( node )
ELSE
  CALL IOS_Ereport( RoutineName,99,                                        &
      'Remove data from queue called on empty queue')
END IF
CALL IOS_Increment_Queue_Length(                                           &
    -1*ios_md_len*umFortranIntegerSize(),no_lock_needed)
CALL release_lock()

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE IOS_Remove_Data_From_Queue

!-----------------------------------------------------------------------
! SUBROUTINE to Get Node from Queue (from the front ;-) )
!-----------------------------------------------------------------------
SUBROUTINE IOS_Get_Node_From_Queue(node)
IMPLICIT NONE

TYPE(IOS_node_type), POINTER  :: node

CALL acquire_lock()
!$OMP FLUSH
IF ( ASSOCIATED (Q_first) ) THEN
  node => Q_first
ELSE
  CALL IOS_Ereport( RoutineName,99,                                        &
      'Tried to get a node from an empty queue')
END IF
CALL release_lock()
END SUBROUTINE IOS_Get_Node_From_Queue

!-----------------------------------------------------------------------
! FUNCTION to size of the queue (this is the payload size in words)
!-----------------------------------------------------------------------
INTEGER FUNCTION IOS_getQueuePayload()
IMPLICIT NONE

CALL acquire_lock()
!$OMP FLUSH
IOS_getQueuePayload=Q_size
CALL release_lock()
RETURN
END FUNCTION IOS_getQueuePayload

!-----------------------------------------------------------------------
! Subroutine to modify the size of the queue in bytes
!-----------------------------------------------------------------------
SUBROUTINE IOS_Increment_Queue_Length(amount,need_lock)
IMPLICIT NONE

INTEGER, INTENT(IN)           :: amount
LOGICAL, INTENT(IN)           :: need_lock
INTEGER                       :: errCode
REAL(KIND=real64)             :: timeStamp
REAL(KIND=jprb):: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='IOS_INCREMENT_QUEUE_LENGTH'
IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( need_lock ) CALL acquire_lock()
!$OMP FLUSH
Q_size=Q_size+amount
!$OMP FLUSH
IF ( Q_size > IOS_buffer_max ) THEN
!$OMP CRITICAL(internal_write)
  WRITE(qmessage,'(A,I0,A,I0,A)')'Adding ',amount,                         &
      ' to the queue exceeded IOS_buffer_max (',                           &
      IOS_buffer_max,')'
!$OMP END CRITICAL(internal_write)
  CALL IOS_Ereport( RoutineName,-99, Qmessage )
ELSE IF ( Q_size < 0 ) THEN
!$OMP CRITICAL(internal_write)
  WRITE(qmessage,'(A,I0,A,I0,A)')                                          &
      ' Decrementing the Q length by ',-1*amount,                          &
      'resulted in a negative length (',Q_Size,')'
!$OMP END CRITICAL(internal_write)
  CALL IOS_Ereport( RoutineName,-99, Qmessage )
END IF
timeStamp=get_wallclock_time()
!$OMP FLUSH

! Intel Fortran seems to have some form of non-thread safety for any
! file output...
!$OMP CRITICAL(OUT_UNIT)
WRITE(logunit,'(F10.3,A,I12,A,I12)')                                          &
    timeStamp-IOS_start_time,' ',Q_size,' ',amount
!$OMP END CRITICAL(OUT_UNIT)
!$OMP FLUSH
IF ((1.0*Q_Size)/(1.0*IOS_buffer_max) > 0.9) THEN
  IF (IOS_Verbosity >= IOS_PrStatus_Oper) THEN
    WRITE(IOS_message,'(A,F8.3,A)')'WARNING: Queue:',                      &
        (100.0*Q_Size)/(1.0*IOS_buffer_max),'% capacity'
    CALL IOS_print(IOS_message,src='ios_queue_mod')
    CALL IOS_print_flush()
  END IF
END IF


IF ( need_lock )CALL release_lock()
IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE IOS_Increment_Queue_Length


!-----------------------------------------------------------------------
! Subroutines to guard against concurrent updates to the queue
! NOTE: these routines are completely nullified for non-OpenMP builds
!       using the OMP sentinel !$ but must still exist to satisfy
!       dependencies
!-----------------------------------------------------------------------
SUBROUTINE acquire_lock()
IMPLICIT NONE

!$    REAL :: t1,t2
!$    IF (IOS_Lock_meter) t1=get_wallclock_time()
!$    CALL   omp_set_lock(Q_lock_var)
!$    IF (IOS_Lock_meter) THEN
!$      t2=get_wallclock_time()
!$      LockTries(1,omp_get_thread_num())=                                       &
!$          LockTries(1,omp_get_thread_num())+1
!$      LockTime(1,omp_get_thread_num())=                                        &
!$          LockTime(1,omp_get_thread_num())+(t2-t1)
!$    END IF
END SUBROUTINE acquire_lock

SUBROUTINE release_lock()
IMPLICIT NONE

!$    CALL  omp_unset_lock(Q_lock_var)
END SUBROUTINE release_lock



!-----------------------------------------------------------------------
! Subroutine to initialize the queue
!-----------------------------------------------------------------------
SUBROUTINE IOS_Queue_initialise(buffer)
IMPLICIT NONE

INTEGER, INTENT(IN)           :: buffer
INTEGER                       :: errCode
CHARACTER (LEN=20)            :: logname
CHARACTER (LEN=80)            :: buffer_env
REAL(KIND=jprb):: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='IOS_QUEUE_INITIALISE'
IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise variables
logname = "dummy filename"
NULLIFY(Q_First)
NULLIFY(Q_Last )
Q_size = 0
Q_items= 0
IOS_AsyncCompletionRequested=0
! Lock the queue if this is an OpenMP run (not relevant otherwise)
!$    CALL omp_init_lock( Q_lock_var)

! Set the maxumum payload size of the queue
! Unit as given is in MB - convert to bytes
IOS_buffer_max = buffer  * 1024 * 1024

! Open a file to log usage to
!$OMP CRITICAL(internal_write)
WRITE(logname,'(A,I5.5)')'ioserver_log.',global_rank
!$OMP END CRITICAL(internal_write)
WRITE(IOS_message,'(A,A)')'IOS_queue: Logging to: ',TRIM(logname)
CALL IOS_print(IOS_message,src='ios_queue_mod')

! Intel fortran needs the critical region
!$OMP CRITICAL(OUT_UNIT)
OPEN(UNIT=logunit,FILE=TRIM(logname))
!$OMP END CRITICAL(OUT_UNIT)

! Setup the lock arrays if this is an OpenMP build
!$    IF (IOS_Lock_meter) THEN
!$      ALLOCATE(LockTime                                                      &
!$          (words_per_cacheline,0:omp_get_max_threads()-1))
!$      ALLOCATE(LockTries                                                     &
!$          (words_per_cacheline,0:omp_get_max_threads()-1))
!$      LockTime (1,0:omp_get_max_threads()-1)=0.0
!$      LockTries(1,0:omp_get_max_threads()-1)=0
!$    END IF

! Record the start time
IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE IOS_Queue_initialise

!-----------------------------------------------------------------------
! Subroutine to finalise the queue
! NOTE: this routine is completely nullified for non-OpenMP builds using
!       the OMP sentinel !$ but must still exist to satisfy dependencies
!-----------------------------------------------------------------------
SUBROUTINE IOS_Queue_Finalise()
IMPLICIT NONE

!$    INTEGER :: thread

!$    CALL omp_destroy_lock( Q_lock_var)

!$    IF (IOS_Lock_meter) THEN
!$      CALL IOS_print('')
!$      CALL IOS_print('IOS Queue Access:')
!$      WRITE(IOS_message,'(A6,A12,A12)')'Thread ','Locks','Time'
!$      CALL IOS_print(IOS_message,src='ios_queue_mod')
!$      WRITE(IOS_message,'(A30)')'------------------------------'
!$      CALL IOS_print(IOS_message,src='ios_queue_mod')
!$      DO thread = 0,omp_get_max_threads()-1
!$OMP CRITICAL(internal_write)
!$        WRITE(IOS_message,'(I6,I12,F12.2)')thread,                             &
!$            LockTries(1,thread),LockTime(1,thread)
!$OMP END CRITICAL(internal_write)
!$        CALL IOS_print(IOS_message,src='ios_queue_mod')
!$      END DO
!$      CALL IOS_print('')
!$    END IF

!$    CALL IOS_print('Info: Queue service terminated.')

END SUBROUTINE IOS_Queue_Finalise

!-----------------------------------------------------------------------
! FUNCTION that reports whether the queue has space to accommodate another
! len words. It may wait until the condition is met, and it may not,
! hence the return code.
!-----------------------------------------------------------------------
LOGICAL FUNCTION HasFreeSpace(amountNeeded)
USE fort2c_interfaces, ONLY: um_sleep
IMPLICIT NONE

INTEGER, INTENT(IN) :: amountNeeded
INTEGER             :: Q_size_copy
REAL                :: t1,t2
REAL(KIND=jprb):: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='HASFREESPACE'
IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

HasFreeSpace=.TRUE.

!$OMP FLUSH
IF ( amountNeeded  > IOS_buffer_max ) THEN
!$OMP CRITICAL(internal_write)
  WRITE(Qmessage,'(A,I0,A,I0)')                                            &
      'The request for ',amountNeeded,                                     &
      ' words of space on the Q cannot be satisfied, maximum set to ',     &
      IOS_buffer_max
!$OMP END CRITICAL(internal_write)
  CALL IOS_Ereport( RoutineName,99, Qmessage )
END IF

! Q_size could be modified by the other thread, so need to do tests on
! it in a thread safe manner. May introduce a little inefficiency, but
! safer.
CALL acquire_lock()
!$OMP FLUSH
Q_size_copy  = Q_size
CALL release_lock()

IF ( Q_size_copy + amountNeeded  > IOS_buffer_max ) THEN
  t1=get_wallclock_time()
!$OMP FLUSH
  IF (IOS_Verbosity>=IOS_PrStatus_Oper) THEN
!$OMP CRITICAL(internal_write)
    WRITE(IOS_Message,'(A,I0,A)')'Info: Queue: Waiting for ',              &
        amountNeeded,' units of free space'
!$OMP END CRITICAL(internal_write)
    CALL IOS_print(IOS_message,src='ios_queue_mod')
  END IF

  DO WHILE(Q_size_copy + amountNeeded  > IOS_buffer_max .AND.              &
      IOS_AsyncCompletionRequested==0)
    CALL um_sleep(IOS_backoff_interval)
    ! Update the value of Q_size we hold
    CALL acquire_lock()
!$OMP FLUSH
    Q_size_copy  = Q_size
    CALL release_lock()
  END DO
  t2=get_wallclock_time()
  IF (IOS_Verbosity>=IOS_PrStatus_Oper) THEN
    WRITE(IOS_Message,'(A,F8.2,A)')'Info: Queue: Wait done: ',             &
        t2-t1,'s'
    CALL IOS_print(IOS_Message,src='ios_queue_mod')

  END IF
  IF ( IOS_AsyncCompletionRequested /= 0 ) HasFreeSpace=.FALSE.
END IF
IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION HasFreeSpace

!-----------------------------------------------------------------------
! FUNCTION that reports whether the queue is fully drained
! It may wait until the condition is met, and it may not.
!-----------------------------------------------------------------------
LOGICAL FUNCTION WaitForDrain()
USE fort2c_interfaces, ONLY: um_sleep
IMPLICIT NONE

INTEGER            :: counter
INTEGER, PARAMETER :: counter_targ = 300
REAL               :: t1
REAL               :: t2
REAL(KIND=jprb)    :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='WAITFORDRAIN'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
WaitForDrain=.TRUE.
counter = 0
t1=get_wallclock_time()
!$OMP FLUSH
IF ( Q_items > 0 ) THEN
  IF (IOS_Verbosity>=IOS_PrStatus_Oper)                                    &
      CALL IOS_print('Info: Queue: Waiting for drain',src='ios_queue_mod')
  DO WHILE(Q_items > 0 .AND. IOS_AsyncCompletionRequested==0)
    CALL um_sleep(IOS_backoff_interval)
    counter=counter+1
    IF (counter==counter_targ) THEN
      counter = 0
      t2=get_wallclock_time()
      IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
        WRITE(IOS_Message,'(A,F8.2,A)')'Info: Queue: Waiting: ',           &
            t2-t1,'s'
        CALL IOS_print(IOS_message,src='ios_queue_mod')
      END IF
    END IF
!$OMP FLUSH
  END DO
  IF ( IOS_AsyncCompletionRequested /= 0 ) WaitForDrain=.FALSE.
  IF (IOS_Verbosity>=IOS_PrStatus_Oper)                                    &
      CALL IOS_print('Info: Queue: done drain',src='ios_queue_mod')
END IF
IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION WaitForDrain

!-----------------------------------------------------------------------
! FUNCTION that reports the number of distinct items in the queue
!-----------------------------------------------------------------------
INTEGER FUNCTION IOS_getQueueItems()
IMPLICIT NONE

CALL acquire_lock()
!$OMP FLUSH
IOS_getQueueItems=Q_items
CALL release_lock()
RETURN
END FUNCTION IOS_getQueueItems

!-----------------------------------------------------------------------
! FUNCTION: A factory that makes nodes to put in the queue
!-----------------------------------------------------------------------
FUNCTION make_new_node() RESULT(n)

IMPLICIT NONE

TYPE(IOS_node_type),POINTER :: n

ALLOCATE(n)

NULLIFY(n%next)
NULLIFY(n%prev)
NULLIFY(n%real_data)
NULLIFY(n%real32_data)
NULLIFY(n%distributed_data)
NULLIFY(n%integer_data)
NULLIFY(n%integer32_data)
NULLIFY(n%receive_requests)
NULLIFY(n%receive_tags)
NULLIFY(n%receive_data_len)

n%metadata%ACTION=IOS_Action_unset
n%transaction_number=-1

END FUNCTION make_new_node

!-----------------------------------------------------------------------
! SUBROUTINE: A scrapyard for recyling old nodes safely
!-----------------------------------------------------------------------
SUBROUTINE destroy_node(node)
IMPLICIT NONE

TYPE(IOS_node_type), POINTER :: node
INTEGER                      :: wordsRemoved

REAL(KIND=jprb):: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DESTROY_NODE'
IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
wordsRemoved=0

IF ( ASSOCIATED(node%real_data) ) THEN
  DEALLOCATE  (node%real_data)
  NULLIFY     (node%real_data)
  wordsRemoved=wordsRemoved-node%metadata%data_size*                       &
      IOS_BytesPerWord64
END IF

IF ( ASSOCIATED(node%integer_data) ) THEN
  DEALLOCATE  (node%integer_data)
  NULLIFY     (node%integer_data)
  wordsRemoved=wordsRemoved-node%metadata%data_size*                       &
      IOS_BytesPerInteger
END IF

IF ( ASSOCIATED(node%real32_data) ) THEN
  DEALLOCATE  (node%real32_data)
  NULLIFY     (node%real32_data)
  wordsRemoved=wordsRemoved-node%metadata%data_size*                       &
      IOS_BytesPerWord32
END IF

IF ( ASSOCIATED(node%integer32_data) ) THEN
  DEALLOCATE  (node%integer32_data)
  NULLIFY     (node%integer32_data)
  wordsRemoved=wordsRemoved-node%metadata%data_size*                       &
      IOS_BytesPerWord32
END IF

IF ( ASSOCIATED(node%receive_requests) ) THEN
  DEALLOCATE  (node%receive_requests)
  NULLIFY     (node%receive_requests)
END IF

IF ( ASSOCIATED(node%receive_tags) ) THEN
  DEALLOCATE (node%receive_tags)
  NULLIFY    (node%receive_tags)
END IF

IF ( ASSOCIATED(node%receive_data_len) ) THEN
  DEALLOCATE  (node%receive_data_len)
  NULLIFY     (node%receive_data_len)
END IF

IF ( ASSOCIATED(node%distributed_data) ) THEN
  DEALLOCATE  (node%distributed_data)
  NULLIFY     (node%distributed_data)
END IF

!$OMP FLUSH
CALL IOS_Increment_Queue_Length(wordsRemoved,no_lock_needed)
!$OMP FLUSH

DEALLOCATE(node)
NULLIFY(node)
IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE destroy_node

!-----------------------------------------------------------------------
! SUBROUTINE: A means of dumping the state of the queue for debug
!-----------------------------------------------------------------------
SUBROUTINE IOS_Queue_Report
IMPLICIT NONE


CALL acquire_lock()
!$OMP FLUSH
!$OMP CRITICAL(internal_write)
WRITE(IOS_message,'(A,I6,A,I6,A,F5.2,A)')                                  &
    'IOS Q: ',                                                             &
    Q_size*8/1024.0/1024.0,' MB in ',                                      &
    Q_items,                                                               &
    ' items (',                                                            &
    (Q_size*100.0)/IOS_buffer_max,'%)'
!$OMP END CRITICAL(internal_write)
CALL IOS_print(IOS_message,src='ios_queue_mod')
IF ( Q_items == 0 .AND. Q_size /= 0 ) THEN
  WRITE(qmessage,'(A)')'Q_items and Q_size mismatch'
  CALL IOS_Ereport( RoutineName,-99, Qmessage )
END IF
IF ( Q_items == 0 ) THEN
  CALL IOS_print('Q contains no transactions ',src='ios_queue_mod')
ELSE IF ( Q_items == 1 ) THEN
!$OMP CRITICAL(internal_write)
  WRITE(IOS_message,'(A,I0)')'Q contains transaction ',                    &
      Q_first%transaction_number
!$OMP END CRITICAL(internal_write)
  CALL IOS_print(IOS_message,src='ios_queue_mod')
ELSE IF ( Q_items >= 1 ) THEN
!$OMP CRITICAL(internal_write)
  WRITE(IOS_message,'(A,I0,A,I0)')'Q contains transactions ',              &
      Q_first%transaction_number,                                          &
      ' to ',Q_last%transaction_number
!$OMP END CRITICAL(internal_write)
  CALL IOS_print(IOS_message,src='ios_queue_mod')
END IF
CALL release_lock()
END SUBROUTINE IOS_Queue_Report

END MODULE IOS_Queue_Mod
