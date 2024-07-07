! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: C96

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module defining interfaces for clients of an IO server
! usually called from io.F90 with exceptions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE ios

USE IOS_Common
USE IOS_Client_Queue
USE IOS_types, ONLY:                                                         &
     IOS_state_size
USE mpl,ONLY:                                                               &
    MPL_Integer,                                                             &
    MPL_Integer8,                                                            &
    MPL_Logical4,                                                            &
    MPL_Logical8,                                                            &
    mpl_real4,                                                               &
    mpl_real8,                                                               &
    mpl_real,                                                                &
    mpl_character,                                                           &
    mpl_status_size
USE um_types, ONLY:                                                          &
    real32, real64, logical32, logical64, integer64, integer32
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE IOS_print_mgr, ONLY:                                                    &
    IOS_print,                                                               &
    IOS_Verbosity,                                                           &
    IOS_message,                                                             &
    IOS_print_start_time

USE fort2c_memcpy_interfaces, ONLY: um_memcpy32, um_memcpy64

USE io_constants, ONLY: ioFileTypeMPIIO

USE IOS_Constants, ONLY:                                                       &
    IOS_Read_NoBroadcast,                                                      &
    IOS_Action_Fence,                                                          &
    IOS_Action_Setpos,                                                         &
    IOS_Action_StashSetPos

USE get_wallclock_time_mod, ONLY: get_wallclock_time
 
IMPLICIT NONE

INTERFACE IOS_write32
MODULE PROCEDURE                                                           &
    IOS_Write32_r,                                                         &
    IOS_Write32_i,                                                         &
    IOS_write32_i2d
END INTERFACE

INTERFACE IOS_write64
MODULE PROCEDURE                                                           &
    IOS_Write64_r,                                                         &
    IOS_Write64_i,                                                         &
    IOS_Write64_l,                                                         &
    IOS_write64_r2d,                                                       &
    IOS_write64_i2d,                                                       &
    IOS_write64_l2d
END INTERFACE

INTERFACE IOS_read32
MODULE PROCEDURE                                                           &
    IOS_Read32_r,                                                          &
    IOS_Read32_i,                                                          &
    IOS_read32_i2d,                                                        &
    IOS_Read32_l
END INTERFACE

INTERFACE IOS_read64
MODULE PROCEDURE                                                           &
    IOS_Read64_r,                                                          &
    IOS_Read64_i,                                                          &
    IOS_Read64_l,                                                          &
    IOS_read64_r2d,                                                        &
    IOS_read64_i2d,                                                        &
    IOS_read64_l2d
END INTERFACE

! params/vars  for dr_hook
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

INTEGER, PRIVATE :: ierror ! generic error code used everywhere

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='IOS'

CONTAINS

SUBROUTINE IOS_fileState(UNIT,location,state)

IMPLICIT NONE
INTEGER, INTENT(IN)          :: UNIT
INTEGER, INTENT(IN)          :: location
INTEGER                      :: qHandle
TYPE(IOS_Status)             :: state
INTEGER(KIND=IOS_tukind),                                                  &
    POINTER                  :: recvBuffer(:)
REAL(KIND=jprb)              :: zhook_handle
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'IOS_FILESTATE'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

qHandle=IOS_init_md(UNIT,location,IOS_Action_Enquire,                      &
    datasize=IOS_State_Size)
CALL IOS_Send(qHandle)
recvBuffer => IOS_attach_recvBuffer(qHandle)
CALL um_memcpy64(state,recvBuffer,IOS_State_Size)

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE IOS_fileState


LOGICAL FUNCTION IOS_Is_Unit_Open(UNIT)

IMPLICIT NONE
INTEGER, INTENT(IN)          :: UNIT
TYPE(IOS_Status)             :: state
REAL(KIND=jprb)              :: zhook_handle
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'IOS_IS_UNIT_OPEN'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

CALL IOS_fileState(UNIT,IOS_No_Location,state)
IOS_Is_unit_open=state%isOpen

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END FUNCTION IOS_Is_Unit_Open


SUBROUTINE IOS_Close(UNIT,NAME,delete,fileType)

IMPLICIT NONE

INTEGER, INTENT(IN)          :: UNIT   ! unit number
INTEGER, INTENT(IN)          :: delete ! =0 read only, otherwise r/w
CHARACTER (LEN=*), INTENT(IN):: NAME    ! nameironment var/filename
INTEGER, OPTIONAL, INTENT(IN):: fileType ! Determine if we need to 
                                         ! distribute the close
INTEGER                      :: targetRank
INTEGER                      :: targetRankMax
INTEGER                      :: errorCode=99
INTEGER                      :: qHandle
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'IOS_CLOSE'
TYPE(IOS_metadata_type),                                                   &
    POINTER                  :: metadata
REAL(KIND=jprb):: zhook_handle

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Check name is small enough to tranmit.
IF ( LEN_TRIM(NAME) > IOS_String_Max ) THEN
  IOS_message    = 'String is too big to transmit'
  CALL IOS_Ereport(RoutineName, ErrorCode, IOS_message)
END IF

IF ( PRESENT(fileType)) THEN
   IF ( IOS_Enable_mpiio .AND. fileType == ioFileTypeMPIIO ) THEN 
      targetRankMax = IOS_Metadata_Receivers()
   ELSE
      targetRankMax = 1
   END IF
ELSE
   targetRankMax = 1
END IF

DO targetRank=1,targetRankMax
   ! Determine the IO server responsible
   qHandle=IOS_init_md(UNIT,-1,IOS_Action_Close,targetRank=targetRank)
   metadata => IOS_Attach_Metadata(qHandle)
   metadata % name_length              = LEN_TRIM(NAME)
   metadata % delete                   = delete
   metadata % string(1:metadata % name_length)                                &
        = NAME(1:metadata % name_length)
   IF ( PRESENT(fileType) ) metadata % address = fileType
   CALL IOS_Send(qHandle,hasString=.TRUE.)
END DO

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE IOS_Close


SUBROUTINE IOS_Fence(UNIT,root)

USE IOS_Comms, ONLY: IOS_SoftSync

IMPLICIT NONE

INTEGER, INTENT(IN)          :: UNIT
INTEGER, INTENT(INOUT)       :: root
INTEGER                      :: qHandle
INTEGER                      :: pe
REAL(KIND=jprb)              :: zhook_handle
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'IOS_FENCE'
TYPE(IOS_metadata_type)      :: metadata

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF ( L_IOS_active() .AND. model_rank == 0 ) THEN

  !The root is the global task id owning the unit, if its owned by
  !atmosphere then its set to zero.
  IF ( root == 0 ) THEN
    root=global_rank ! because atm rank 0 is not COMM_WORLD 0
  ELSE
    root=io_server_for_unit(UNIT)
  END IF

!$OMP CRITICAL(internal_write)
  WRITE(IOS_message,'(A,I3,A,I5)')'Fencing unit ',UNIT,' on rank ',root
!$OMP END CRITICAL(internal_write)
  CALL IOS_print(IOS_message,src='ios')

  DO pe=1,IOS_Server_Groups

    !recycle address to encode 'owner' in this case
    qHandle=IOS_init_md(-1*io_servers(pe,1),root,IOS_Action_Fence)
    CALL IOS_Send(qHandle)

  END DO

  ! If, I, atmos rank 0, own the file I have to participate in the
  ! operation. Of course I have no queue, so I don't need to wait

  IF ( root == global_rank ) THEN
    CALL IOS_print('Atmosphere waiting for fence',src='ios')
!$      CALL IOS_SoftSync(root)
  END IF

END IF

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE IOS_Fence


SUBROUTINE IOS_DiskSync(UNIT,fileType,wait)

IMPLICIT NONE

INTEGER, INTENT(IN)          :: UNIT! unit number
INTEGER, INTENT(IN)          :: fileType
LOGICAL, OPTIONAL            :: wait! length of env
INTEGER                      :: errorCode=99
INTEGER                      :: qHandle
INTEGER                      :: targetRank
INTEGER                      :: targetRankMax
LOGICAL                      :: myWaitFlag
REAL(KIND=jprb)              :: zhook_handle
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'IOS_DISKSYNC'
TYPE(IOS_metadata_type),                                                   &
    POINTER                  :: metadata

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

myWaitFlag=.TRUE. ! Default is to block (ie slow but safe)
IF ( PRESENT(wait) ) THEN
  IF ( .NOT. wait ) THEN
    myWaitFlag=.FALSE.
  END IF
END IF

IF ( IOS_Enable_mpiio .AND. fileType == ioFileTypeMPIIO ) THEN
   targetRankMax=IOS_Metadata_Receivers()
ELSE
   targetRankMax=1
END IF

DO targetRank=targetRankMax,1,-1
   ! To avoid deadlocks when waiting for the response from 
   ! IOS_Action_Sync_Barrier, only IO rank 0 communicates back to the 
   ! client. Therefore, rank 0 should recieve its metadata last.
   IF ( myWaitFlag ) THEN
      qHandle=IOS_init_md(UNIT,-1,IOS_Action_Sync_Barrier,datasize=1,      &
           targetRank=targetRank)
   ELSE
      qHandle=IOS_init_md(UNIT,-1,IOS_Action_Sync,targetRank=targetRank)
   END IF
   ! Need the filetype to determine whether to use sync_single or
   ! MPL_File_sync
   metadata => IOS_Attach_Metadata(qHandle)
   metadata % address = fileType
   CALL IOS_Send(qHandle)
END DO

! Return code - assume OK, errors handled by IO Server

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE IOS_DiskSync


FUNCTION IOS_Open(UNIT,NAME,read_write,filetype,allowRemap) RESULT(r)

IMPLICIT NONE

INTEGER, INTENT(IN)          :: UNIT      ! unit number
INTEGER, INTENT(IN)          :: read_write! ==0 read only,
                                          ! otherwise r/w
INTEGER, INTENT(IN)          :: fileType
LOGICAL, OPTIONAL            :: allowRemap
INTEGER                      :: r
INTEGER                      :: qHandle
INTEGER                      :: targetRank
INTEGER                      :: targetRankMax
INTEGER                      :: errorCode=99
REAL(KIND=jprb)              :: zhook_handle
CHARACTER (LEN=*), INTENT(IN):: NAME ! environment var/filename
TYPE(IOS_metadata_type),                                                   &
    POINTER                  :: metadata
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'IOS_OPEN'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Check name is small enough to transmit.
IF (LEN_TRIM(NAME) > IOS_String_Max) THEN
!$OMP CRITICAL(internal_write)
  WRITE(IOS_Message,'(A,I6,A,I6,A)')                                       &
      'Filename name too big (',                                           &
      LEN_TRIM(NAME),' > maximum of ', IOS_String_Max,')'
!$OMP END CRITICAL(internal_write)
  CALL IOS_Ereport(RoutineName, ErrorCode, IOS_message)
END IF

IF (PRESENT(allowRemap)) THEN
  IF (allowRemap) THEN
    CALL IOS_assign_server_for_unit(UNIT)
  END IF
END IF

! If we are using a dynamic allocation, then we need to
! check and assign a server for the unit. All ranks must
! call this
IF (io_server_for_unit(UNIT) == IOS_No_Server)                             &
    CALL IOS_assign_Server_for_unit(UNIT)

r=-1
IF ( IOS_Enable_mpiio .AND. fileType == ioFileTypeMPIIO ) THEN  
   targetRankMax=IOS_Metadata_Receivers()
ELSE
   targetRankMax=1
END IF

IF (model_rank==0) THEN
   DO targetRank=1,targetRankMax
      qHandle  = IOS_init_md(UNIT,-1,IOS_Action_Open,targetRank=targetRank)
      metadata => IOS_Attach_Metadata(qHandle)
      metadata % name_length              = LEN_TRIM(NAME)
      metadata % subtype                  = read_write
      metadata % address                  = fileType
      metadata % string(1:metadata % name_length)                              &
           = NAME(1:metadata % name_length)
         ! Only return the rank of the lead IO task
      IF( targetRank == 1 ) r = IOS_GetDestPe(qHandle)
      CALL IOS_Send(qHandle,hasString=.TRUE.)
   END DO   
END IF

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END FUNCTION IOS_Open


SUBROUTINE IOS_Process(timestep)

USE io_configuration_mod, ONLY: print_runtime_info

IMPLICIT NONE
INTEGER, INTENT(IN)          :: timestep! current timestep
INTEGER                      :: server
INTEGER                      :: subtask
INTEGER                      :: qHandle
REAL                         :: time_now
REAL(KIND=jprb)              :: zhook_handle
TYPE(IOS_metadata_type),                                                   &
    POINTER                  :: metadata
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'IOS_PROCESS'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF (model_rank == 0 ) THEN

  DO server=1,IOS_Server_Groups
    DO subtask=1,IOS_Metadata_Receivers()
      qHandle=IOS_init_md(-1*io_servers(server,subtask),                   &
          -1,IOS_Action_Process)
      CALL IOS_Send(qHandle)
    END DO
  END DO
  IF (IOS_print_start_time) THEN
    time_now = get_wallclock_time()
    IF ( .NOT. print_runtime_info ) THEN
      CALL IOS_print('',src='ios')
      CALL IOS_Print( '****************************************'           &
           //'****************************************',src='ios')
      CALL IOS_print('',src='ios')
    END IF
!$OMP CRITICAL(internal_write)
    WRITE(IOS_message,'(A,I8,A,F10.3,A)')                                  &
        'IO server: Info: Entering timestep ',timestep,                    &
        ' at time=',time_now-IOS_Start_time,' seconds'
!$OMP END CRITICAL(internal_write)
    CALL IOS_print(IOS_message,src='ios')
    CALL IOS_print('',src='ios')
  END IF  ! IOS_print_start_time
END IF

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE IOS_Process


SUBROUTINE IOS_Read32_r(UNIT, location, array, numWords, RdBroadCast)

IMPLICIT NONE
INTEGER, INTENT(IN)          :: UNIT! unit
INTEGER, INTENT(IN)          :: location
INTEGER, INTENT(IN)          :: numWords
INTEGER, INTENT(IN)          :: RdBroadCast
REAL(KIND=real32),                                                         &
    INTENT(OUT)              :: array(:)! data read
INTEGER                      :: qHandle
INTEGER                      :: errorCode
REAL(KIND=jprb)              :: zhook_handle
INTEGER(KIND=IOS_tukind),                                                  &
    POINTER                  :: recvBuffer(:)
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'IOS_READ32_R'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF (model_rank==0) THEN
  qHandle=IOS_init_md(UNIT,location,IOS_Action_Read32,                     &
      dataSize=numWords,bcast=RdBroadCast)
  CALL IOS_Send(qHandle)
END IF

IF (RdBroadCast == IOS_Read_NoBroadcast) THEN
  recvBuffer => IOS_attach_recvBuffer(qHandle)
  CALL um_memcpy32(array,recvBuffer,numWords)
ELSE
  CALL mpl_bcast(array,                                                    &
      numWords*IOS_tuperword32,                                            &
      IOS_tutype,                                                          &
      bcast_procs-1,                                                       &
      IOS_BCast_Comm(IOServerNo(io_server_for_unit(UNIT))),                &
      errorCode)
END IF

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE IOS_Read32_r


SUBROUTINE IOS_Read32_i(UNIT, location, array, numWords, RdBroadCast)

IMPLICIT NONE

INTEGER, INTENT(IN)        :: UNIT! unit
INTEGER, INTENT(IN)        :: location
INTEGER, INTENT(IN)        :: numWords
INTEGER, INTENT(IN)        :: RdBroadCast
INTEGER                    :: errorCode
INTEGER(KIND=integer32),                                                   &
    INTENT(OUT)            :: array(:)! data read
INTEGER                    :: qHandle
REAL(KIND=jprb)            :: zhook_handle
INTEGER(KIND=IOS_tukind),                                                  &
    POINTER                :: recvBuffer(:)
CHARACTER (LEN=*),                                                         &
    PARAMETER              :: RoutineName = 'IOS_READ32_I'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF (model_rank==0) THEN
  qHandle=IOS_init_md(UNIT,location,IOS_Action_Read32,                     &
      dataSize=numWords,bcast=RdBroadCast)
  CALL IOS_Send(qHandle)
END IF

IF (RdBroadCast == IOS_Read_NoBroadcast) THEN
  recvBuffer => IOS_attach_recvBuffer(qHandle)
  CALL um_memcpy32(array,recvBuffer,numWords)
ELSE
  CALL mpl_bcast(array,                                                    &
      numWords*IOS_tuperword32,                                            &
      IOS_tutype,                                                          &
      bcast_procs-1,                                                       &
      IOS_BCast_Comm(IOServerNo(io_server_for_unit(UNIT))),                &
      errorCode)
END IF

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE IOS_Read32_i


SUBROUTINE IOS_read32_i2d(UNIT, location, array, numWords, RdBroadCast)

IMPLICIT NONE

INTEGER, INTENT(IN)        :: UNIT! unit
INTEGER, INTENT(IN)        :: location
INTEGER, INTENT(IN)        :: numWords
INTEGER, INTENT(IN)        :: RdBroadCast
INTEGER(KIND=integer32),                                                   &
    INTENT(OUT)            :: array(:,:)! data read
INTEGER                    :: qHandle
INTEGER                    :: errorCode
INTEGER                    :: STATUS(mpl_status_size)
REAL(KIND=jprb)            :: zhook_handle
INTEGER(KIND=IOS_tukind),                                                  &
    POINTER                :: recvBuffer(:)
CHARACTER (LEN=*),                                                         &
    PARAMETER              :: RoutineName = 'IOS_READ32_I2D'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF (model_rank==0) THEN
  qHandle=IOS_init_md(UNIT,location,IOS_Action_Read32,                     &
      dataSize=numWords,bcast=RdBroadCast)
  CALL IOS_Send(qHandle)
END IF

IF (RdBroadCast == IOS_Read_NoBroadcast) THEN
  recvBuffer => IOS_attach_recvBuffer(qHandle)
  CALL um_memcpy32(array,recvBuffer,numWords)
ELSE
  CALL mpl_bcast(array,                                                    &
      numWords*IOS_tuperword32,                                            &
      IOS_tutype,                                                          &
      bcast_procs-1,                                                       &
      IOS_BCast_Comm(IOServerNo(io_server_for_unit(UNIT))),                &
      errorCode)
END IF

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE IOS_read32_i2d


SUBROUTINE IOS_Read32_l(UNIT, location, array, numWords, RdBroadCast)

IMPLICIT NONE

INTEGER, INTENT(IN)        :: UNIT! unit
INTEGER, INTENT(IN)        :: location
INTEGER, INTENT(IN)        :: numWords
INTEGER, INTENT(IN)        :: RdBroadCast
LOGICAL(KIND=logical32),                                                   &
    INTENT(OUT)            :: array(:)! data read
INTEGER                    :: qHandle
INTEGER                    :: errorCode
REAL(KIND=jprb)            :: zhook_handle
INTEGER(KIND=IOS_tukind),                                                  &
    POINTER                :: recvBuffer(:)
CHARACTER (LEN=*),                                                         &
    PARAMETER              :: RoutineName = 'IOS_READ32_L'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF (model_rank==0) THEN
  qHandle=IOS_init_md(UNIT,location,IOS_Action_Read32,                     &
      dataSize=numWords,bcast=RdBroadCast)
  CALL IOS_Send(qHandle)
END IF

IF (RdBroadCast == IOS_Read_NoBroadcast) THEN
  recvBuffer => IOS_attach_recvBuffer(qHandle)
  CALL um_memcpy32(array,recvBuffer,numWords)
ELSE
  CALL mpl_bcast(array,                                                    &
      numWords*IOS_tuperword32,                                            &
      IOS_tutype,                                                          &
      bcast_procs-1,                                                       &
      IOS_BCast_Comm(IOServerNo(io_server_for_unit(UNIT))),                &
      errorCode)
END IF

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE IOS_Read32_l


SUBROUTINE IOS_Read64_r(UNIT, location, array, numWords, RdBroadCast)

IMPLICIT NONE

INTEGER, INTENT(IN)        :: UNIT! unit
INTEGER, INTENT(IN)        :: location
INTEGER, INTENT(IN)        :: numWords
INTEGER, INTENT(IN)        :: RdBroadCast
REAL(KIND=real64),                                                         &
    INTENT(OUT)            :: array(:)! data read
INTEGER                    :: qHandle
INTEGER                    :: errorCode
REAL(KIND=jprb)            :: zhook_handle
INTEGER(KIND=IOS_tukind),                                                  &
    POINTER                :: recvBuffer(:)
CHARACTER (LEN=*),                                                         &
    PARAMETER              :: RoutineName = 'IOS_READ64_R'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF (model_rank==0) THEN
  qHandle=IOS_init_md(UNIT,location,IOS_Action_Read64,                     &
      dataSize=numWords,bcast=RdBroadCast)
  CALL IOS_Send(qHandle)
END IF

IF (RdBroadCast == IOS_Read_NoBroadcast) THEN
  recvBuffer => IOS_attach_recvBuffer(qHandle)
  CALL um_memcpy64(array,recvBuffer,numWords)
ELSE
  CALL mpl_bcast(array,                                                    &
      numWords*IOS_tuperword64,                                            &
      IOS_tutype,                                                          &
      bcast_procs-1,                                                       &
      IOS_BCast_Comm(IOServerNo(io_server_for_unit(UNIT))),                &
      errorCode)
END IF

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE IOS_Read64_r


SUBROUTINE IOS_Read64_l(UNIT, location, array, numWords, RdBroadCast)

IMPLICIT NONE

INTEGER, INTENT(IN)        :: UNIT! unit
INTEGER, INTENT(IN)        :: location
INTEGER, INTENT(IN)        :: numWords
INTEGER, INTENT(IN)        :: RdBroadCast
LOGICAL(KIND=logical64),                                                   &
    INTENT(OUT)            :: array(:)! data read
INTEGER                    :: qHandle
INTEGER                    :: errorCode
REAL(KIND=jprb)            :: zhook_handle
INTEGER(KIND=IOS_tukind),                                                  &
    POINTER                :: recvBuffer(:)
CHARACTER (LEN=*),                                                         &
    PARAMETER              :: RoutineName = 'IOS_READ64_L'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF (model_rank==0) THEN
  qHandle=IOS_init_md(UNIT,location,IOS_Action_Read64,                     &
      dataSize=numWords,bcast=RdBroadCast)
  CALL IOS_Send(qHandle)
END IF

IF (RdBroadCast == IOS_Read_NoBroadcast) THEN
  recvBuffer => IOS_attach_recvBuffer(qHandle)
  CALL um_memcpy64(array,recvBuffer,numWords)
ELSE
  CALL mpl_bcast(array,                                                    &
      numWords*IOS_tuperword64,                                            &
      IOS_tutype,                                                          &
      bcast_procs-1,                                                       &
      IOS_BCast_Comm(IOServerNo(io_server_for_unit(UNIT))),                &
      errorCode)
END IF

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE IOS_Read64_l


SUBROUTINE IOS_Read64_i(UNIT, location, array, numWords, RdBroadCast)

IMPLICIT NONE

INTEGER, INTENT(IN)        :: UNIT! unit
INTEGER, INTENT(IN)        :: location
INTEGER, INTENT(IN)        :: numWords
INTEGER, INTENT(IN)        :: RdBroadCast
INTEGER(KIND=integer64),                                                   &
    INTENT(OUT)            :: array(:)! data read
INTEGER                    :: qHandle
INTEGER                    :: errorCode
REAL(KIND=jprb)            :: zhook_handle
INTEGER(KIND=IOS_tukind),                                                  &
    POINTER                :: recvBuffer(:)
TYPE(IOS_metadata_type),                                                   &
    POINTER                :: metadata
CHARACTER (LEN=*),                                                         &
    PARAMETER              :: RoutineName = 'IOS_READ64_I'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF (model_rank==0) THEN
  qHandle=IOS_init_md(UNIT,location,IOS_Action_Read64,                     &
      dataSize=numWords,bcast=RdBroadCast)
  CALL IOS_Send(qHandle)
END IF

IF (RdBroadCast == IOS_Read_NoBroadcast) THEN
  recvBuffer => IOS_attach_recvBuffer(qHandle)
  CALL um_memcpy64(array,recvBuffer,numWords)
ELSE
  CALL mpl_bcast(array,                                                    &
      numWords*IOS_tuperword64,                                            &
      IOS_tutype,                                                          &
      bcast_procs-1,                                                       &
      IOS_BCast_Comm(IOServerNo(io_server_for_unit(UNIT))),                &
      errorCode)
END IF

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE IOS_Read64_i


SUBROUTINE IOS_read64_r2d(UNIT, location, array, numWords, RdBroadCast)

IMPLICIT NONE

INTEGER, INTENT(IN)        :: UNIT! unit
INTEGER, INTENT(IN)        :: location
INTEGER, INTENT(IN)        :: numWords
INTEGER, INTENT(IN)        :: RdBroadCast
REAL(KIND=real64),                                                         &
    INTENT(OUT)            :: array(:,:)! data read
INTEGER                    :: qHandle
INTEGER                    :: errorCode
REAL(KIND=jprb)            :: zhook_handle
INTEGER(KIND=IOS_tukind),                                                  &
    POINTER                :: recvBuffer(:)
CHARACTER (LEN=*),                                                         &
    PARAMETER              :: RoutineName = 'IOS_READ64_R2D'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF (model_rank==0) THEN
  qHandle=IOS_init_md(UNIT,location,IOS_Action_Read64,                     &
      dataSize=numWords,bcast=RdBroadCast)
  CALL IOS_Send(qHandle)
END IF

IF (RdBroadCast == IOS_Read_NoBroadcast) THEN
  recvBuffer => IOS_attach_recvBuffer(qHandle)
  CALL um_memcpy64(array,recvBuffer,numWords)
ELSE
  CALL mpl_bcast(array,                                                    &
      numWords*IOS_tuperword64,                                            &
      IOS_tutype,                                                          &
      bcast_procs-1,                                                       &
      IOS_BCast_Comm(IOServerNo(io_server_for_unit(UNIT))),                &
      errorCode)
END IF

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE IOS_read64_r2d


SUBROUTINE IOS_read64_l2d(UNIT, location, array, numWords, RdBroadCast)

IMPLICIT NONE

INTEGER, INTENT(IN)        :: UNIT! unit
INTEGER, INTENT(IN)        :: location
INTEGER, INTENT(IN)        :: numWords
INTEGER, INTENT(IN)        :: RdBroadCast
LOGICAL(KIND=logical64),                                                   &
    INTENT(OUT)            :: array(:,:)! data read
INTEGER                    :: qHandle
INTEGER                    :: errorCode
REAL(KIND=jprb)            :: zhook_handle
INTEGER(KIND=IOS_tukind),                                                  &
    POINTER                :: recvBuffer(:)
CHARACTER (LEN=*),                                                         &
    PARAMETER              :: RoutineName = 'IOS_READ64_L2D'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF (model_rank==0) THEN
  qHandle=IOS_init_md(UNIT,location,IOS_Action_Read64,                     &
      dataSize=numWords,bcast=RdBroadCast)
  CALL IOS_Send(qHandle)
END IF

IF (RdBroadCast == IOS_Read_NoBroadcast) THEN
  recvBuffer => IOS_attach_recvBuffer(qHandle)
  CALL um_memcpy64(array,recvBuffer,numWords)
ELSE
  CALL mpl_bcast(array,                                                    &
      numWords*IOS_tuperword64,                                            &
      IOS_tutype,                                                          &
      bcast_procs-1,                                                       &
      IOS_BCast_Comm(IOServerNo(io_server_for_unit(UNIT))),                &
      errorCode)
END IF

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE IOS_read64_l2d


SUBROUTINE IOS_read64_i2d(UNIT, location, array, numWords, RdBroadCast)

IMPLICIT NONE

INTEGER, INTENT(IN)        :: UNIT! unit
INTEGER, INTENT(IN)        :: location
INTEGER, INTENT(IN)        :: numWords
INTEGER, INTENT(IN)        :: RdBroadCast
INTEGER(KIND=integer64),                                                   &
    INTENT(OUT)            :: array(:,:)! data read
INTEGER                    :: qHandle
INTEGER                    :: errorCode
REAL(KIND=jprb)            :: zhook_handle
INTEGER(KIND=IOS_tukind),                                                  &
    POINTER                :: recvBuffer(:)
CHARACTER (LEN=*),                                                         &
    PARAMETER              :: RoutineName = 'IOS_READ64_I2D'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF (model_rank==0) THEN
  qHandle=IOS_init_md(UNIT,location,IOS_Action_Read64,                     &
      dataSize=numWords,bcast=RdBroadCast)
  CALL IOS_Send(qHandle)
END IF

IF (RdBroadCast == IOS_Read_NoBroadcast) THEN
  recvBuffer => IOS_attach_recvBuffer(qHandle)
  CALL um_memcpy64(array,recvBuffer,numWords)
ELSE
  CALL mpl_bcast(array,                                                    &
      numWords*IOS_tuperword64,                                            &
      IOS_tutype,                                                          &
      bcast_procs-1,                                                       &
      IOS_BCast_Comm(IOServerNo(io_server_for_unit(UNIT))),                &
      errorCode)
END IF

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE IOS_read64_i2d


SUBROUTINE IOS_Setpos(UNIT, pos, icode)
! Deprecated? completely obsolete?
IMPLICIT NONE

INTEGER, INTENT(IN)        :: UNIT! unit
INTEGER, INTENT(IN)        :: pos! file seek position
INTEGER, INTENT(OUT)       :: icode! error flag
INTEGER                    :: qHandle
REAL(KIND=jprb)            :: zhook_handle
TYPE(IOS_metadata_type)    :: metadata
CHARACTER (LEN=*),                                                         &
    PARAMETER              :: RoutineName = 'IOS_SETPOS'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Set metadata
qHandle=IOS_init_md(UNIT,pos,IOS_Action_Setpos)
CALL IOS_Send(qHandle)
icode = 0

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE IOS_Setpos


SUBROUTINE IOS_Write32_r(UNIT, location, array, numWords)

IMPLICIT NONE

INTEGER, INTENT(IN)        :: UNIT! unit
INTEGER, INTENT(IN)        :: location
INTEGER, INTENT(IN)        :: numWords
REAL(KIND=real32),                                                         &
    INTENT(IN)             :: array(:)! data to write
INTEGER                    :: qHandle
REAL(KIND=jprb)            :: zhook_handle
INTEGER(KIND=IOS_tukind),                                                  &
    POINTER                :: sendBuffer(:)
CHARACTER (LEN=*),                                                         &
    PARAMETER              :: RoutineName = 'IOS_WRITE32_R'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

qHandle=IOS_init_md(UNIT,location,                                         &
    IOS_Action_Write32, dataSize=numWords)
sendBuffer => IOS_attach_SendBuffer(qHandle)
CALL um_memcpy32(sendBuffer,array,numWords)
CALL IOS_Send(qHandle)

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE IOS_Write32_r


SUBROUTINE IOS_Write32_i(UNIT, location, array, numWords)

IMPLICIT NONE

INTEGER, INTENT(IN)        :: UNIT! unit
INTEGER, INTENT(IN)        :: location
INTEGER, INTENT(IN)        :: numWords
INTEGER(KIND=integer32),                                                   &
    INTENT(IN)             :: array(:)! data to write
INTEGER                    :: qHandle
REAL(KIND=jprb)            :: zhook_handle
INTEGER(KIND=IOS_tukind),                                                  &
    POINTER                :: sendBuffer(:)
CHARACTER (LEN=*),                                                         &
    PARAMETER              :: RoutineName = 'IOS_WRITE32_I'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

qHandle=IOS_init_md(UNIT,location,                                         &
    IOS_Action_Write32, dataSize=numWords)
sendBuffer => IOS_attach_SendBuffer(qHandle)
CALL um_memcpy32(sendBuffer,array,numWords)
CALL IOS_Send(qHandle)

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE IOS_Write32_i


SUBROUTINE IOS_write32_i2d(UNIT, location, array, numWords)

IMPLICIT NONE

INTEGER, INTENT(IN)        :: UNIT! unit
INTEGER, INTENT(IN)        :: location
INTEGER, INTENT(IN)        :: numWords
INTEGER(KIND=integer32),                                                   &
    INTENT(IN)           :: array(:,:)! data to write
INTEGER                    :: qHandle
INTEGER(KIND=IOS_tukind),                                                  &
    POINTER                :: sendBuffer(:)
REAL(KIND=jprb)            :: zhook_handle
CHARACTER (LEN=*),                                                         &
    PARAMETER              :: RoutineName = 'IOS_WRITE32_I2D'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

qHandle=IOS_init_md(UNIT,location,                                         &
    IOS_Action_Write32, dataSize=numWords)
sendBuffer => IOS_attach_SendBuffer(qHandle)
CALL um_memcpy32(sendBuffer,array,numWords)
CALL IOS_Send(qHandle)

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE IOS_write32_i2d


SUBROUTINE IOS_Write64_r(UNIT, location, array, numWords)

IMPLICIT NONE

INTEGER, INTENT(IN)        :: UNIT! unit
INTEGER, INTENT(IN)        :: location
INTEGER, INTENT(IN)        :: numWords
REAL(KIND=real64),                                                         &
    INTENT(IN)             :: array(:)! data to write
INTEGER                    :: qHandle
REAL(KIND=jprb)            :: zhook_handle
INTEGER(KIND=IOS_tukind),                                                  &
    POINTER                :: sendBuffer(:)
CHARACTER (LEN=*),                                                         &
    PARAMETER              :: RoutineName = 'IOS_WRITE64_R'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

qHandle=IOS_init_md(UNIT,location,IOS_Action_Write64,dataSize=numWords)
sendBuffer => IOS_attach_SendBuffer(qHandle)
CALL um_memcpy64(sendBuffer,array,numWords)
CALL IOS_Send(qHandle)

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE IOS_Write64_r


SUBROUTINE IOS_Write64_i(UNIT, location, array, numWords)

IMPLICIT NONE

INTEGER, INTENT(IN)        :: UNIT! unit
INTEGER, INTENT(IN)        :: location
INTEGER, INTENT(IN)        :: numWords
INTEGER(KIND=integer64),                                                   &
    INTENT(IN)             :: array(:)! data to write
INTEGER                    :: qHandle
REAL(KIND=jprb)            :: zhook_handle
INTEGER(KIND=IOS_tukind),                                                  &
    POINTER                :: sendBuffer(:)
CHARACTER (LEN=*),                                                         &
    PARAMETER              :: RoutineName = 'IOS_WRITE64_I'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

qHandle=IOS_init_md(UNIT,location,IOS_Action_Write64,dataSize=numWords)
sendBuffer => IOS_attach_SendBuffer(qHandle)
CALL um_memcpy64(sendBuffer,array,numWords)
CALL IOS_Send(qHandle)

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE IOS_Write64_i


SUBROUTINE IOS_Write64_l(UNIT, location, array, numWords)

IMPLICIT NONE

INTEGER, INTENT(IN)        :: UNIT! unit
INTEGER, INTENT(IN)        :: location
INTEGER, INTENT(IN)        :: numWords
LOGICAL(KIND=logical64),                                                   &
    INTENT(IN)             :: array(:)! data to write
INTEGER                    :: qHandle
REAL(KIND=jprb)            :: zhook_handle
INTEGER(KIND=IOS_tukind),                                                  &
    POINTER                :: sendBuffer(:)
CHARACTER (LEN=*),                                                         &
    PARAMETER              :: RoutineName = 'IOS_WRITE64_L'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

qHandle=IOS_init_md(UNIT,location,IOS_Action_Write64,dataSize=numWords)
sendBuffer => IOS_attach_SendBuffer(qHandle)
CALL um_memcpy64(sendBuffer,array,numWords)
CALL IOS_Send(qHandle)

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE IOS_Write64_l


!2D...............
SUBROUTINE IOS_write64_r2d(UNIT, location, array, numWords)

IMPLICIT NONE

INTEGER, INTENT(IN)        :: UNIT! unit
INTEGER, INTENT(IN)        :: location
INTEGER, INTENT(IN)        :: numWords
REAL(KIND=real64),                                                         &
    INTENT(IN)             :: array(:,:)! data to write
INTEGER                    :: qHandle
REAL(KIND=jprb)            :: zhook_handle
INTEGER(KIND=IOS_tukind),                                                  &
    POINTER                :: sendBuffer(:)
CHARACTER (LEN=*),                                                         &
    PARAMETER              :: RoutineName = 'IOS_WRITE64_R2D'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

qHandle=IOS_init_md(UNIT,location,IOS_Action_Write64,dataSize=numWords)
sendBuffer => IOS_attach_SendBuffer(qHandle)
CALL um_memcpy64(sendBuffer,array,numWords)
CALL IOS_Send(qHandle)

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE IOS_write64_r2d


SUBROUTINE IOS_write64_i2d(UNIT, location, array, numWords)

IMPLICIT NONE

INTEGER, INTENT(IN)        :: UNIT! unit
INTEGER, INTENT(IN)        :: location
INTEGER, INTENT(IN)        :: numWords
INTEGER(KIND=integer64),                                                   &
    INTENT(IN)             :: array(:,:)! data to write
INTEGER                    :: qHandle
REAL(KIND=jprb)            :: zhook_handle
INTEGER(KIND=IOS_tukind),                                                  &
    POINTER                :: sendBuffer(:)
CHARACTER (LEN=*),                                                         &
    PARAMETER              :: RoutineName = 'IOS_WRITE64_I2D'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

qHandle=IOS_init_md(UNIT,location,IOS_Action_Write64,dataSize=numWords)
sendBuffer => IOS_attach_SendBuffer(qHandle)
CALL um_memcpy64(sendBuffer,array,numWords)
CALL IOS_Send(qHandle)

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE IOS_write64_i2d


SUBROUTINE IOS_write64_l2d(UNIT, location, array, numWords)

IMPLICIT NONE

INTEGER, INTENT(IN)        :: UNIT! unit
INTEGER, INTENT(IN)        :: location
INTEGER, INTENT(IN)        :: numWords
LOGICAL(KIND=logical64),                                                   &
    INTENT(IN)             :: array(:,:)! data to write
INTEGER                    :: qHandle
REAL(KIND=jprb)            :: zhook_handle
INTEGER(KIND=IOS_tukind),                                                  &
    POINTER                :: sendBuffer(:)
CHARACTER (LEN=*),                                                         &
    PARAMETER              :: RoutineName = 'IOS_WRITE64_L2D'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

qHandle=IOS_init_md(UNIT,location,IOS_Action_Write64,dataSize=numWords)
sendBuffer => IOS_attach_SendBuffer(qHandle)
CALL um_memcpy64(sendBuffer,array,numWords)
CALL IOS_Send(qHandle)

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE IOS_write64_l2d


SUBROUTINE IOS_Init_Header(UNIT,LEN)

USE IOS_stash_common, ONLY: isUsingAsyncStash

IMPLICIT NONE

INTEGER, INTENT(IN)        :: UNIT
INTEGER, INTENT(IN)        :: LEN
INTEGER                    :: qHandle
REAL(KIND=jprb)            :: zhook_handle
CHARACTER (LEN=*),                                                         &
    PARAMETER              :: RoutineName = 'IOS_INIT_HEADER'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF (model_rank==0 .AND. L_IOS_Active() .AND. isUsingAsyncStash()) THEN
  qHandle=IOS_init_md(UNIT,LEN,IOS_Action_StashInitHeader)
  CALL IOS_Send(qHandle)
END IF

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE IOS_Init_Header


SUBROUTINE IOS_Init_PP_Lookup(args)

USE IOS_stash_common, ONLY: isUsingAsyncStash

IMPLICIT NONE

INTEGER, PARAMETER         :: num_args=5
INTEGER, INTENT(IN)        :: args(num_args)
INTEGER                    :: qHandle
REAL(KIND=jprb)            :: zhook_handle
INTEGER(KIND=IOS_tukind),                                                  &
    POINTER                :: sendBuffer(:)
CHARACTER (LEN=*),                                                         &
    PARAMETER              :: RoutineName = 'IOS_INIT_PP_LOOKUP'


IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF ( model_rank == 0 .AND.                                                 &
    L_IOS_Active()  .AND.                                                  &
    isUsingAsyncStash()) THEN

  qHandle=IOS_init_md(args(1),-1,IOS_Action_StashInitPPLookup,             &
      dataSize=num_args)
  sendBuffer => IOS_attach_SendBuffer(qHandle)
  CALL um_memcpy64(sendBuffer,args,num_args)
  CALL IOS_Send(qHandle)

END IF

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE IOS_Init_PP_Lookup


SUBROUTINE IOS_Set_Header(UNIT,header_data)

USE IOS_stash_common, ONLY: isUsingAsyncStash

IMPLICIT NONE

INTEGER, INTENT(IN)          :: UNIT
INTEGER, POINTER             :: header_data(:)
INTEGER                      :: qHandle
REAL(KIND=jprb)              :: zhook_handle
INTEGER(KIND=IOS_tukind),                                                  &
    POINTER                  :: sendBuffer(:)
CHARACTER (LEN=*),                                                         &
    PARAMETER                :: RoutineName = 'IOS_SET_HEADER'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF ( model_rank == 0 .AND.                                                 &
    L_IOS_Active()  .AND.                                                  &
    (.NOT. l_io_server) .AND.                                              &
    isUsingAsyncStash()) THEN

  qHandle=IOS_init_md(UNIT,-1,IOS_Action_StashSetHeader,                   &
      dataSize=SIZE(header_data))
  sendBuffer => IOS_attach_SendBuffer(qHandle)
  CALL um_memcpy64(sendBuffer,header_data,SIZE(header_data))
  CALL IOS_Send(qHandle)

END IF

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE IOS_Set_Header


SUBROUTINE IOS_Setpos_Stash(UNIT, theRecord, icode)

USE IOS_stash_common, ONLY: isUsingAsyncStash

IMPLICIT NONE

INTEGER, INTENT(IN)        :: UNIT! unit
INTEGER, INTENT(IN)        :: theRecord! file seek position
INTEGER, INTENT(OUT)       :: icode! error flag
INTEGER                    :: qHandle
REAL(KIND=jprb)            :: zhook_handle
CHARACTER (LEN=*),                                                         &
    PARAMETER              :: RoutineName = 'IOS_SETPOS_STASH'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF ( model_rank == 0 .AND.                                                 &
    L_IOS_Active()  .AND.                                                  &
    isUsingAsyncStash()) THEN

  qHandle=IOS_init_md(UNIT,theRecord,IOS_Action_StashSetPos)
  CALL IOS_Send(qHandle)

END IF

icode = 0

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE IOS_Setpos_Stash


SUBROUTINE IOS_MergeLookup(UNIT, array,isize)

USE IOS_stash_common, ONLY: isUsingAsyncStash

IMPLICIT NONE

INTEGER, INTENT(IN)        :: UNIT! unit
INTEGER, INTENT(IN)        :: isize! data size
INTEGER, INTENT(IN)        :: array(:,:)! data to write
INTEGER                    :: io_slave_pe
INTEGER                    :: qHandle
REAL(KIND=jprb)            :: zhook_handle
TYPE(IOS_metadata_type)    :: metadata
INTEGER(KIND=IOS_tukind),                                                  &
    POINTER                :: sendBuffer(:)
CHARACTER (LEN=*),                                                         &
    PARAMETER              :: RoutineName = 'IOS_MERGELOOKUP'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF ( model_rank == 0 .AND.                                                 &
    L_IOS_Active()  .AND.                                                  &
    isUsingAsyncStash()) THEN

  qHandle=IOS_init_md(UNIT,IOS_No_Location,IOS_action_mergepplookup,       &
      dataSize=isize)
  sendBuffer => IOS_attach_SendBuffer(qHandle)
  CALL um_memcpy64(sendBuffer,array,isize)
  CALL IOS_Send(qHandle)

END IF

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE IOS_MergeLookup


SUBROUTINE IOS_Write_Lookup(UNIT, location, array)

USE IOS_stash_common, ONLY: isUsingAsyncStash

IMPLICIT NONE

INTEGER, INTENT(IN)        :: UNIT! unit
INTEGER, INTENT(IN)        :: location
INTEGER, INTENT(IN)        :: array(:,:)! data to write
INTEGER                    :: qHandle
REAL(KIND=jprb)            :: zhook_handle
INTEGER(KIND=IOS_tukind),                                                  &
    POINTER                :: sendBuffer(:)
CHARACTER (LEN=*),                                                         &
    PARAMETER              :: RoutineName = 'IOS_WRITE_LOOKUP'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF ( model_rank == 0 .AND.                                                 &
    L_IOS_Active()  .AND.                                                  &
    isUsingAsyncStash()) THEN

  qHandle=IOS_init_md(UNIT,location,IOS_ACTION_StashWritePPLookup,         &
      dataSize=SIZE(array))
  sendBuffer => IOS_attach_SendBuffer(qHandle)
  CALL um_memcpy64(sendBuffer,array,SIZE(array))
  CALL IOS_Send(qHandle)

END IF

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE IOS_Write_Lookup


SUBROUTINE IOS_Write_Stash_PrePacked_Data                                    &
    (UNIT ,theRecord, array)

USE IOS_stash_common, ONLY: isUsingAsyncStash

IMPLICIT NONE

INTEGER, INTENT(IN)        :: UNIT! unit
INTEGER, INTENT(IN)        :: theRecord! the stash record
REAL, INTENT(IN)           :: array(:)! data to write
INTEGER                    :: qHandle
REAL(KIND=jprb)            :: zhook_handle
INTEGER(KIND=IOS_tukind),                                                  &
    POINTER                :: sendBuffer(:)
TYPE(IOS_metadata_type),                                                   &
    POINTER                :: metadata
CHARACTER (LEN=*),                                                         &
    PARAMETER              :: RoutineName = 'IOS_WRITE_STASH_PREPACKED_DATA'

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF ( model_rank == 0 .AND.                                                 &
    L_IOS_Active()  .AND.                                                  &
    isUsingAsyncStash()) THEN

  ! Set metadata (note we use location to encode record)
  qHandle=IOS_init_md(UNIT,theRecord,                                      &
      IOS_Action_Write_PP_Prepacked,                                       &
      datasize = SIZE(array))
  sendBuffer => IOS_attach_SendBuffer(qHandle)
  CALL um_memcpy64(sendBuffer,array,SIZE(array))
  CALL IOS_Send(qHandle)

END IF

IF ( lhook ) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE IOS_Write_Stash_PrePacked_Data

END MODULE ios

