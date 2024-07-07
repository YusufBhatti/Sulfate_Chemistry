! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: C96

MODULE IOS_comms
!$ USE omp_lib
USE IOS_Common
USE IOS_Constants, ONLY:  &
    ios_md_len,           &
    IOS_String_max
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE IOS_print_mgr

USE get_wallclock_time_mod, ONLY: get_wallclock_time
 
IMPLICIT NONE

INTERFACE ios_bcast
MODULE PROCEDURE                                                           &
    IOS_BCast_i,                                                           &
    ios_bcast_md,                                                          &
    ios_bcast_str
END INTERFACE

PRIVATE IOS_BCast_i
PRIVATE ios_BCast_md
PRIVATE IOS_BCast_str

! Only needed by OpenMP builds
!$  INTEGER(KIND=omp_lock_kind),                                               &
!$      PRIVATE                   :: IOS_comms_access_lock

  ! Lock metering
INTEGER, PARAMETER, PRIVATE   :: cachepad = 16
INTEGER, POINTER, PRIVATE     :: LockTries(:,:)
REAL, POINTER, PRIVATE        :: LockTime (:,:)

! params/vars  for dr_hook
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='IOS_COMMS'

CONTAINS

SUBROUTINE IOS_Comms_Init()
IMPLICIT NONE
! NOTE: this routine is completely nullified for non-OpenMP builds using
!       the OMP sentinel !$ but must still exist to satisfy dependencies
!$    CALL omp_init_lock( IOS_comms_access_lock)

!$    IF (IOS_Lock_meter) THEN
!$      ALLOCATE(LockTime(cachepad,0:omp_get_max_threads()-1))
!$      ALLOCATE(LockTries(cachepad,0:omp_get_max_threads()-1))
!$      LockTime (1,0:omp_get_max_threads()-1)=0.0
!$      LockTries(1,0:omp_get_max_threads()-1)=0
!$    END IF

END SUBROUTINE IOS_Comms_Init


SUBROUTINE IOS_Comms_Fini()
IMPLICIT NONE
! NOTE: this routine is completely nullified for non-OpenMP builds using
!       the OMP sentinel !$ but must still exist to satisfy dependencies
!$    INTEGER :: thread

!$    CALL omp_destroy_lock( IOS_comms_access_lock)

!$    IF (IOS_Lock_meter) THEN
!$      CALL IOS_print(' ',src='ios_comms')
!$      CALL IOS_print('MPI Access:',src='ios_comms')
!$      WRITE(IOS_message,'(A6,A12,A12)')'Thread ','Locks','Time'
!$      CALL IOS_print(IOS_message,src='ios_comms')
!$      WRITE(IOS_message,'(A30)')'------------------------------'
!$      CALL IOS_print(IOS_message,src='ios_comms')
!$      DO thread=0,omp_get_max_threads()-1
!$OMP CRITICAL(internal_write)
!$        WRITE(IOS_message,'(I6,I12,F12.2)')                                 &
!$            thread,LockTries(1,thread),LockTime(1,thread)
!$OMP END CRITICAL(internal_write)
!$        CALL IOS_print(IOS_message,src='ios_comms')
!$      END DO
!$      DEALLOCATE(LockTime)
!$      DEALLOCATE(LockTries)
!$    END IF

END SUBROUTINE IOS_Comms_Fini

!----------------------------------------------------------------------
! Subroutines to guard against concurrent calls to MPI
! (which some MPI doesn't like)
!----------------------------------------------------------------------
SUBROUTINE acquire_lock_MPI()
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE
! NOTE: this routine is completely nullified for non-OpenMP builds using
!       the OMP sentinel !$ but must still exist to satisfy dependencies
!$    REAL :: t1,t2
!$    CHARACTER (LEN=errormessagelength)           :: message
!$    CHARACTER (LEN=*),                                                         &
!$        PARAMETER                 :: RoutineName = 'ACQUIRE_LOCK_MPI'

!$    IF ( IOS_thread_0_calls_mpi .AND. omp_get_thread_num() /= 0 ) THEN
!$OMP CRITICAL(internal_write)
!$      WRITE(message,'(A,I2,A)')                                             &
!$          'Thread ',omp_get_thread_num(),                                   &
!$          ' tried to call mpi, but only thread 0 is allowed'
!$OMP END CRITICAL(internal_write)
!$      CALL IOS_Ereport( RoutineName,99, message )
!$    END IF

!$    IF ( IOS_serialise_mpi_calls ) THEN
!$      IF (IOS_Lock_meter) t1=get_wallclock_time()
!$      CALL omp_set_lock(IOS_comms_access_lock)
!$      IF (IOS_Lock_meter) THEN
!$        t2=get_wallclock_time()
!$        LockTries (1,omp_get_thread_num())=                                    &
!$            LockTries (1,omp_get_thread_num())+1
!$        LockTime  (1,omp_get_thread_num())=                                    &
!$            LockTime  (1,omp_get_thread_num())+(t2-t1)
!$      END IF
!$    END IF

END SUBROUTINE acquire_lock_MPI

SUBROUTINE release_lock_MPI()
IMPLICIT NONE
! NOTE: this routine is completely nullified for non-OpenMP builds using
!       the OMP sentinel !$ but must still exist to satisfy dependencies
!$    IF ( IOS_serialise_mpi_calls ) CALL omp_unset_lock(IOS_comms_access_lock)

END SUBROUTINE release_lock_MPI

! Wait for a thing without thread hogging MPI.
SUBROUTINE IOS_Wait(request)
IMPLICIT NONE
INTEGER, INTENT(INOUT)  :: request
INTEGER                 :: req(1)
req(1)=request
CALL IOS_WaitAll(req,SIZE(req))
request=req(1)
END SUBROUTINE IOS_Wait

! Wait for some things without thread hogging MPI.
SUBROUTINE ios_waitall(requests,items)
USE mpl, ONLY: mpl_status_size
IMPLICIT NONE
INTEGER, INTENT(IN   )  :: items
INTEGER, INTENT(INOUT)  :: requests(items)
INTEGER, ALLOCATABLE    :: STATUS(:,:)
INTEGER                 :: errorStatus
LOGICAL                 :: flag
CHARACTER (LEN=*),                                                         &
    PARAMETER                   :: RoutineName = 'IOS_WAITALL'


IF (items>0) THEN
  ALLOCATE(STATUS(mpl_status_size,items))
  flag=.FALSE.

  DO WHILE (.NOT. flag)
    CALL acquire_lock_mpi()
    CALL MPL_Testall(                                                      &
        items,                                                             &
        requests,                                                          &
        flag,                                                              &
        STATUS,                                                            &
        errorStatus)
    CALL release_lock_mpi()
  END DO
  DEALLOCATE(STATUS)
END IF

END SUBROUTINE ios_waitall

! Gather integer without thread hogging MPI.
SUBROUTINE IOS_gather(gatheredData,myDatum)
USE mpl, ONLY: mpl_integer
IMPLICIT NONE
INTEGER, INTENT(IN)  :: myDatum
INTEGER, INTENT(OUT) :: gatheredData(:)
INTEGER              :: errorStatus
INTEGER              :: tag
INTEGER              :: proc
INTEGER              :: num_requests
INTEGER, ALLOCATABLE :: requests(:)
INTEGER, PARAMETER   :: rootPE=0
CHARACTER (LEN=*),                                                         &
    PARAMETER                   :: RoutineName = 'IOS_GATHER'

tag=1

IF (model_rank==0) THEN

  gatheredData(1)=myDatum
  num_requests = model_procs-1
  ALLOCATE (requests(num_requests))

  IF (model_procs>1) THEN
    CALL acquire_lock_mpi()
    DO proc=1,model_procs-1
      CALL MPL_iRecv(                                                      &
          gatheredData(proc+1),                                            &
          1,                                                               &
          mpl_integer,                                                     &
          proc,                                                            &
          tag,                                                             &
          model_comm,                                                      &
          requests(proc),                                                  &
          errorStatus)
    END DO
    CALL release_lock_mpi()
  END IF

ELSE

  num_requests = 1
  ALLOCATE (requests(num_requests))
  CALL acquire_lock_mpi()
  CALL MPL_Isend(                                                          &
      myDatum,                                                             &
      1,                                                                   &
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

END SUBROUTINE IOS_gather


SUBROUTINE IOS_Bcast_md(theData)
USE mpl, ONLY: mpl_integer
IMPLICIT NONE

TYPE(IOS_metadata_type),                                                   &
    INTENT(INOUT)      :: theData
TYPE(IOS_metadata_type),                                                   &
    ALLOCATABLE        :: theData_copy(:)   ! Copies of theData for sending
INTEGER                :: errorStatus
INTEGER                :: tag
INTEGER                :: proc
INTEGER                :: num_requests
INTEGER, ALLOCATABLE   :: requests(:)
INTEGER, PARAMETER     :: rootPE=0
REAL(KIND=jprb)        :: zhook_handle
CHARACTER (LEN=*),                                                         &
    PARAMETER          :: RoutineName = 'IOS_BCAST_MD'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
tag=1

IF (model_rank==0) THEN

  num_requests = model_procs-1
  ALLOCATE (requests(num_requests))
  ALLOCATE (theData_copy(num_requests))

  IF (model_procs>1) THEN
    DO proc=1,model_procs-1
      ! Need to copy the data to send to avoid multiple open sends on
      ! the same data - as this is outlawed by the MPI standard
      theData_Copy(proc) = theData

      CALL acquire_lock_mpi()
      CALL MPL_iSend(                                                      &
          theData_copy(proc),                                              &
          ios_md_len,                                                      &
          mpl_integer,                                                     &
          proc,                                                            &
          tag,                                                             &
          model_comm,                                                      &
          requests(proc),                                                  &
          errorStatus)
      CALL release_lock_mpi()

    END DO
  END IF

ELSE

  num_requests = 1
  ALLOCATE (requests(num_requests))

  CALL acquire_lock_mpi()
  CALL MPL_iRecv(                                                          &
      theData,                                                             &
      ios_md_len,                                                          &
      mpl_integer,                                                         &
      rootPE,                                                              &
      tag,                                                                 &
      model_comm,                                                          &
      requests(1),                                                         &
      errorStatus)
  CALL release_lock_mpi()

END IF

CALL IOS_WaitAll(requests, num_requests)

IF (model_rank == 0) DEALLOCATE(theData_copy)
DEALLOCATE (requests)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE IOS_Bcast_md

SUBROUTINE IOS_Bcast_i(theData)
USE mpl, ONLY: mpl_integer
IMPLICIT NONE

INTEGER, INTENT(INOUT) :: theData(:)
INTEGER                :: errorStatus
INTEGER                :: tag
INTEGER                :: proc
INTEGER                :: num_requests
INTEGER, ALLOCATABLE   :: requests(:)
INTEGER, ALLOCATABLE   :: thedata_copy(:,:)  ! copies of thedata for sending
INTEGER, PARAMETER     :: rootPE=0
REAL(KIND=jprb)        :: zhook_handle
CHARACTER (LEN=*),                                                         &
    PARAMETER          :: RoutineName = 'IOS_BCAST_I'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
tag=1

IF (model_rank==0) THEN

  num_requests = model_procs-1
  ALLOCATE (requests(num_requests))
  ALLOCATE (theData_copy(SIZE(theData),num_requests))

  IF (model_procs>1) THEN
    DO proc=1,model_procs-1
      ! Need to copy the data to send to avoid multiple open sends on
      ! the same data - as this is outlawed by the MPI standard
      theData_copy(:,proc) = theData(:)

      CALL acquire_lock_mpi()
      CALL MPL_iSend(                                                      &
          theData_copy(:,proc),                                            &
          SIZE(theData),                                                   &
          mpl_integer,                                                     &
          proc,                                                            &
          tag,                                                             &
          model_comm,                                                      &
          requests(proc),                                                  &
          errorStatus)
      CALL release_lock_mpi()

    END DO
  END IF

ELSE

  num_requests = 1
  ALLOCATE (requests(num_requests))

  CALL acquire_lock_mpi()
  CALL MPL_iRecv(                                                          &
      theData,                                                             &
      SIZE(theData),                                                       &
      mpl_integer,                                                         &
      rootPE,                                                              &
      tag,                                                                 &
      model_comm,                                                          &
      requests(1),                                                         &
      errorStatus)
  CALL release_lock_mpi()

END IF

CALL IOS_WaitAll(requests, num_requests)

IF (model_rank == 0) DEALLOCATE(theData_copy)
DEALLOCATE (requests)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE IOS_Bcast_i

SUBROUTINE IOS_Bcast_str(theData)
USE mpl, ONLY: mpl_character
IMPLICIT NONE

CHARACTER(len=IOS_string_max),                                             &
    INTENT(INOUT)      :: theData
CHARACTER(len=IOS_string_max),                                             &
    ALLOCATABLE        :: theData_copy(:)   ! Copies of theData for sending
INTEGER                :: errorStatus
INTEGER                :: tag
INTEGER                :: proc
INTEGER                :: num_requests
INTEGER, ALLOCATABLE   :: requests(:)
INTEGER, PARAMETER     :: rootPE=0
REAL(KIND=jprb)        :: zhook_handle
CHARACTER (LEN=*),                                                         &
    PARAMETER          :: RoutineName = 'IOS_BCAST_STR'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
tag=1

IF (model_rank==0) THEN

  num_requests = model_procs-1
  ALLOCATE (requests(num_requests))
  ALLOCATE (theData_copy(num_requests))

  IF (model_procs>1) THEN
    DO proc=1,model_procs-1
      ! Need to copy the data to send to avoid multiple open sends on
      ! the same data - as this is outlawed by the MPI standard
      theData_Copy(proc) = theData

      CALL acquire_lock_mpi()
      CALL MPL_iSend(                                                      &
          theData_copy(proc),                                              &
          ios_string_max,                                                  &
          mpl_character,                                                   &
          proc,                                                            &
          tag,                                                             &
          model_comm,                                                      &
          requests(proc),                                                  &
          errorStatus)
      CALL release_lock_mpi()

    END DO
  END IF

ELSE

  num_requests = 1
  ALLOCATE (requests(num_requests))
 
  CALL acquire_lock_mpi()
  CALL MPL_iRecv(                                                          &
      theData,                                                             &
      ios_string_max,                                                      &
      mpl_character,                                                       &
      rootPE,                                                              &
      tag,                                                                 &
      model_comm,                                                          &
      requests(1),                                                         &
      errorStatus)
  CALL release_lock_mpi()

END IF

CALL IOS_WaitAll(requests, num_requests)

IF (model_rank == 0) DEALLOCATE(theData_copy)
DEALLOCATE (requests)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE IOS_Bcast_str

! Limited synchronisation primitive
!
! Collective over rank 0 of all IO servers and the process identified
! by "root" on the global communicator.
!
! Purpose: Delay root until other involved processes make the call.
SUBROUTINE IOS_SoftSync(root)
USE mpl, ONLY: mpl_integer, mpl_request_null, mpl_status_size
USE IOS_Constants, ONLY:IOS_ControlSync_Tag
IMPLICIT NONE

INTEGER, INTENT(IN)  :: root ! The global rank owning the epoch
INTEGER              :: server_pe
INTEGER              :: i
INTEGER              :: ierror
INTEGER              :: flag
INTEGER              :: num_requests
INTEGER, SAVE        :: epoch_count=0
INTEGER, ALLOCATABLE :: requests(:)
REAL(KIND=jprb)      :: zhook_handle
CHARACTER (LEN=*),                                                         &
    PARAMETER        :: RoutineName = 'IOS_SOFTSYNC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

epoch_count=epoch_count+1

IF (root<0 .OR. root>= global_procs) THEN
  CALL IOS_ereport(RoutineName, 99, 'Invalid root process provided' )
END IF

flag=1
IF (global_rank==root) THEN
  num_requests = IOS_Server_Groups
  ALLOCATE(requests(num_requests))
  CALL acquire_lock_mpi()
  DO i=1,IOS_Server_Groups
    server_pe=io_servers(i,1)
    IF (server_pe/=global_rank) THEN
      CALL MPL_iRecv(flag,1,mpl_integer,                                   &
          server_pe,IOS_ControlSync_Tag,                                   &
          global_comm,requests(i),ierror)
    ELSE
      requests(i)=mpl_request_null
    END IF
  END DO
  CALL release_lock_mpi()
  IF (IOS_Verbosity>=IOS_PrStatus_Oper) THEN
!$OMP CRITICAL(internal_write)
    WRITE(IOS_message,'(A,I0)')'Info: SoftSync: Waiting for '//            &
        'other IOS to pass epoch ',epoch_count
!$OMP END CRITICAL(internal_write)
    CALL IOS_print(IOS_message,src='ios_comms')
  END IF
ELSE
  num_requests = 1
  ALLOCATE(requests(num_requests))
  CALL acquire_lock_mpi()
  CALL MPL_iSend(flag,1,mpl_integer,root,                                  &
      IOS_ControlSync_Tag,global_comm,requests(1),ierror)
  CALL release_lock_mpi()
END IF

CALL IOS_WaitAll(requests, num_requests)
DEALLOCATE(requests)

IF (IOS_Verbosity>=IOS_PrStatus_Oper) THEN
  IF (global_rank==root) THEN
!$OMP CRITICAL(internal_write)
    WRITE(IOS_message,'(A,I0,A)')                                          &
        'Info: Epoch ',epoch_count,' complete.'
!$OMP END CRITICAL(internal_write)
    CALL IOS_print(IOS_message,src='ios_comms')
  ELSE
!$OMP CRITICAL(internal_write)
    WRITE(IOS_message,'(A,I0,A,I0,A)')'Info: Signalled to',root,           &
        ' epoch ',epoch_count,' complete.'
!$OMP END CRITICAL(internal_write)
    CALL IOS_print(IOS_message,src='ios_comms')
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE IOS_SoftSync

END MODULE IOS_comms
