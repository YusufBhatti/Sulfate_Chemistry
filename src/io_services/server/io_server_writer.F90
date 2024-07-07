! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: C96

! Discharge events from an IOS queue
MODULE IO_Server_Writer

USE io, ONLY: setpos, fileState, is_unit_open, File_Open, buffin, buffout,   &
    buffin32, buffout32, iofilestate, file_close, mpiio_fh
USE io_constants, ONLY: &
    ioFileTypeUM, ioFileTypeMPIIO, ioAllLocal, ioNameProvided
USE IOS_Common
USE IOS_Constants
USE IOS_types, ONLY:                                                        &
    IOS_State_Size
USE IOS_Queue_Mod
USE IOS_Stash_Server, ONLY:                                                 &
    IOS_Stash_Server_Init_Recvs,                                             &
    IOS_Stash_Server_Wait_For_Data,                                          &
    IOS_Stash_server_Process,                                                &
    IOS_Stash_server_Deallocate
USE IOS_Model_Geometry, ONLY:                                               &
    IOS_DumpInitModel
USE MPPIO_file_utils, ONLY:                                                 &
    file_action
USE mpl, ONLY:                                                              &
    mpl_integer,                                                             &
    mpl_character
USE model_file,     ONLY:                                                   &
    model_file_open,                                                         &
    model_file_close,                                                        &
    setRecordDiskLength,                                                     &
    setRecordDiskStart,                                                      &
    attachLookups,                                                           &
    attachHeader,                                                            &
    initLookups,                                                             &
    setLookups,                                                              &
    initHeader,                                                              &
    FixedHeader,                                                             &
    setHeader
USE file_manager, ONLY: get_file_by_unit, um_file_type, assign_file_unit,    &
                        release_file_unit
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE IOS_print_mgr
USE lookup_addresses

USE errormessagelength_mod, ONLY: errormessagelength

USE get_wallclock_time_mod, ONLY: get_wallclock_time
 
IMPLICIT NONE

! Due to the asynchronous nature of the server the listener and
! writer timestep independently. Whilst we do not make use of
! timestepping in this version it is an important concept in many
! coupling frameworks.
INTEGER :: timestep

! params/vars  for dr_hook
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='IO_SERVER_WRITER'

CONTAINS
SUBROUTINE IOS_writer()
USE IOS_Comms, ONLY:                                                      &
    IOS_SoftSync,                                                          &
    acquire_lock_mpi,                                                      &
    release_lock_mpi
USE IOS_Constants, ONLY:  IOS_Action_Strlen
USE ios_common, ONLY: needing_lock
USE fort2c_interfaces, ONLY: um_sleep
USE fort2c_portio_interfaces, ONLY: sync_single

IMPLICIT NONE

TYPE(IOS_node_type), POINTER      :: node
INTEGER, POINTER                  :: ipplook(:,:)
INTEGER, POINTER                  :: fixed_header(:)
INTEGER                           :: ierror
INTEGER                           :: len_io
INTEGER                           :: record_id
LOGICAL                           :: done
TYPE(fileState)                   :: fState
REAL                              :: io_stat
REAL                              :: t
REAL                              :: t1
REAL                              :: t2
REAL(KIND=jprb)                   :: zhook_handle
CHARACTER (LEN=*), PARAMETER      :: RoutineName = 'IOS_WRITER'
CHARACTER (LEN=errormessagelength)               :: IOS_Wrt_Message
CHARACTER (LEN=IOS_Action_Strlen) :: action_name
TYPE(um_file_type), POINTER       :: um_file

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
timestep=0
NULLIFY(ipplook)
NULLIFY(fixed_header)
done=.FALSE.

DO WHILE(.NOT. done)

  IF (IOS_getQueueItems() == 0) THEN
    CALL um_sleep(IOS_backoff_interval)
  ELSE

    CALL IOS_Get_Node_From_Queue( node )
    IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
      t1=get_wallclock_time()
      action_name=TRIM(IOS_ActionName(node%metadata%ACTION))
!$OMP CRITICAL(internal_write)
      WRITE(IOS_message,'(A,I8,A,I3,A,A)')                                 &
          'Info: Writer: transaction: ',                                   &
          node%transaction_number,' is for unit ',                         &
          node%metadata%UNIT,' :',                                         &
          action_name
!$OMP END CRITICAL(internal_write)
      CALL IOS_print(IOS_message,src='io_server_writer')
      CALL IOS_print_flush()
    END IF

    IF (node%metadata%UNIT>0) THEN
      IF (node%metadata%ACTION /= IOS_Action_StashWritePPData              &
          .AND.                                                            &
          node%metadata%ACTION /= IOS_Action_StashWriteDumpData            &
          .AND. (.NOT. IOS_enable_mpiio .OR.                               &
          ( node%metadata%ACTION /= IOS_Action_Open                        &
          .AND.                                                            &
            node%metadata%ACTION /= IOS_Action_Sync                        &
          .AND.                                                            &
            node%metadata%ACTION /= IOS_Action_Sync_Barrier                &
          .AND.                                                            &
            node%metadata%ACTION /= IOS_Action_Close ) )                   &
          ) THEN
        IF (global_rank/=io_server_for_unit(node%metadata%UNIT)) THEN
!$OMP CRITICAL(internal_write)
          WRITE(IOS_wrt_message,'(A,I3,A,I6)')                             &
              'Transaction for unit ',                                     &
              node%metadata%UNIT,                                          &
              ' came to me, but I think it should be for PE:',             &
              io_server_for_unit(node%metadata%UNIT)
!$OMP END CRITICAL(internal_write)
          CALL IOS_Ereport( RoutineName, 10,                         &
              IOS_wrt_message, md=node%metadata )
        END IF
      END IF
    END IF

    SELECT CASE (node%metadata%ACTION)

    CASE (IOS_Action_assign_unit)
      io_server_for_unit_lookup(node%metadata%subtype)=                    &
          node%metadata%address
!$OMP CRITICAL(internal_write)
      WRITE(IOS_message,'(A,I3,A,I4)')                                     &
          'Info: Writer: unit ',node%metadata%subtype,                     &
          ' is assigned to IO server rank ',                               &
          node%metadata%address
!$OMP END CRITICAL(internal_write)
      CALL IOS_print(IOS_message,src='io_server_writer')
      CALL IOS_print_flush()
!$OMP FLUSH

                !--------------------------------------------------------------
                ! Process
                !--------------------------------------------------------------

    CASE (IOS_Action_Process)
      timestep=timestep+1
      t=get_wallclock_time()
      IF ( IOS_Verbosity >= IOS_PrStatus_Oper ) THEN
!$OMP CRITICAL(internal_write)
        WRITE(IOS_message,'(A,I4,A,F10.3)')                                &
            'Info: Writer: Entering Timestep ',timestep,                   &
            ' at time=',t-IOS_start_time
!$OMP END CRITICAL(internal_write)
        CALL IOS_print(IOS_message,src='io_server_writer')
      END IF
      !--------------------------------------------------------------
      ! Write out data (64 bit)
      !--------------------------------------------------------------

    CASE (IOS_Action_Write64,                                              &
        IOS_Action_Write_PP_Prepacked)

      IF (node%metadata%address >= 0 ) THEN
        ! Metadata contains embedded setpos address
        IF (node%metadata%ACTION == IOS_Action_Write_PP_Prepacked) THEN

          ! We are given a stash record in the address field.
          ! So we should decode and record its location in the lookup table.

          record_id=node%metadata%address

          IF (record_id==1) THEN
            fixed_header=>attachHeader(node%metadata%UNIT,FixedHeader)
            node%metadata%address=Fixed_Header(160)-1
            NULLIFY(fixed_header)
          ELSE
            ipplook=>attachLookups(node%metadata%UNIT)
            node%metadata%address=                                         &
                ipplook(lbegin, record_id-1)+                              &
                ipplook(lbnrec, record_id-1)
            NULLIFY(ipplook)
          END IF

          IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
!$OMP CRITICAL(internal_write)
            WRITE(IOS_message,'(A,I8)')                                    &
                'Info: Writer: prepacked data record=',                    &
                record_id
!$OMP END CRITICAL(internal_write)
            CALL IOS_print(IOS_message,src='io_server_writer')
!$OMP CRITICAL(internal_write)
            WRITE(IOS_message,'(A,I8)')                                    &
                'Info: Writer:  decoded to disk address ',                 &
                node%metadata%address
!$OMP END CRITICAL(internal_write)
            CALL IOS_print(IOS_message,src='io_server_writer')
          END IF

          CALL setRecordDiskStart(                                         &
              node%metadata%UNIT,                                          &
              record_id,                                                   &
              node%metadata % address)
        END IF ! Whether this was stash
        IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
!$OMP CRITICAL(internal_write)
          WRITE(IOS_message,'(A,I12)')                                     &
              'Info: Writer: Process embedded setpos: ',                   &
              node%metadata%address
!$OMP END CRITICAL(internal_write)
          CALL IOS_print(IOS_message,src='io_server_writer')
        END IF
        CALL Setpos( node%metadata%UNIT,                                   &
            node%metadata%address, ierror)
      END IF ! Whether the address field was positive

      IF (.NOT. ASSOCIATED(node%real_data)) THEN
        WRITE(IOS_wrt_message,'(A)')                                       &
            'FAULT: there is no data on this (Action_Write) node'
        CALL IOS_Ereport( RoutineName, 10, IOS_Wrt_Message,          &
            md=node%metadata )
      END IF

      IO_Stat = -1.0
      Len_IO  = node%metadata%data_size

      CALL Buffout(node%metadata%UNIT, node%real_data,                     &
          node%metadata%data_size, Len_IO, IO_Stat )

      IF (Len_IO /= node%metadata%data_size .OR. IO_Stat/=-1.0) THEN
!$OMP CRITICAL(internal_write)
        WRITE(IOS_Wrt_Message,'(A,I4,I10,I10,F6.2)')                       &
            'Error in Writing to file on unit ',                           &
            node%metadata%UNIT, Len_IO,                                    &
            node%metadata%data_size, IO_Stat
!$OMP END CRITICAL(internal_write)
        CALL IOS_Ereport( RoutineName, 11, IOS_Wrt_Message,          &
            md=node%metadata )
      END IF

      ! For stash ops. We were given a record number, so we need to
      ! update the correct record with the length of the write
      IF (node%metadata%ACTION==IOS_Action_Write_PP_Prepacked)             &
          CALL setRecordDiskLength(                                        &
          node%metadata%UNIT,                                              &
          record_id,                                                       &
          node%metadata % data_size)

      !--------------------------------------------------------------
      ! Write out data (32 bit)
      !--------------------------------------------------------------

    CASE (IOS_Action_Write32)

      IF (node%metadata%address >= 0 ) THEN
        ! Metadata contains embedded setpos address

        IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
!$OMP CRITICAL(internal_write)
          WRITE(IOS_message,'(A,I12)')                                     &
              'Info: Writer: Process embedded setpos: ',                   &
              node%metadata%address
!$OMP END CRITICAL(internal_write)
          CALL IOS_print(IOS_message,src='io_server_writer')
        END IF
        CALL Setpos( node%metadata%UNIT,                                   &
            node%metadata%address, ierror)
      END IF

      IF (.NOT. ASSOCIATED(node%real32_data)) THEN
        WRITE(IOS_wrt_message,'(A)')                                       &
            'FAULT: there is no data on this (Action_Write32) node'
        CALL IOS_Ereport( RoutineName, 10, IOS_Wrt_Message,          &
            md=node%metadata )
      END IF
      IO_Stat = -1.0
      Len_IO  = node%metadata%data_size

      CALL Buffout32(node%metadata%UNIT, node%real32_data,                 &
          node%metadata%data_size, Len_IO, IO_Stat )

      IF (Len_IO /= node%metadata%data_size .OR. IO_Stat/=-1.0) THEN
!$OMP CRITICAL(internal_write)
        WRITE(IOS_Wrt_Message,'(A,I4,I10,I10,F6.2)')                       &
            'Error in Writing to file on unit ',                           &
            node%metadata%UNIT, Len_IO,                                    &
            node%metadata%data_size, IO_Stat
!$OMP END CRITICAL(internal_write)
        CALL IOS_Ereport( RoutineName, 12, IOS_Wrt_Message,          &
            md=node%metadata )
      END IF

      !--------------------------------------------------------------
      ! Read Data
      !--------------------------------------------------------------

    CASE (IOS_Action_Read64)

      IF (node%metadata%address >= 0 ) THEN
        ! Metadata contains embedded setpos address
        IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
!$OMP CRITICAL(internal_write)
          WRITE(IOS_message,'(A,I12)')                                     &
              'Info: Writer: Process embedded setpos: ',                   &
              node%metadata%address
!$OMP END CRITICAL(internal_write)
          CALL IOS_print(IOS_message,src='io_server_writer')
        END IF
        CALL Setpos( node%metadata%UNIT,                                   &
            node%metadata%address, ierror)
      END IF

      CALL IOS_Increment_Queue_Length                                      &
          (node%metadata%data_size*IOS_BytesPerWord64,needing_lock)
      ALLOCATE(node%real_data( node%metadata%data_size) )
      CALL Buffin( node%metadata%UNIT,                                     &
          node%real_data,                                                  &
          node%metadata%data_size,                                         &
          len_IO, IO_stat )
      IF (.NOT. IOS_thread_0_calls_mpi) THEN
        CALL acquire_lock_mpi()
        IF (node%metadata%subtype==IOS_Read_Broadcast) THEN
          CALL mpl_bcast(                                                  &
              node%real_data,                                              &
              node%metadata%data_size*IOS_tuperword64,                     &
              IOS_tutype,                                                  &
              bcast_procs-1,                                               &
              IOS_BCast_Server_Comm,                                       &
              ierror)
        ELSE
          CALL MPL_Send(                                                   &
              node%real_data,                                              &
              node%metadata%data_size*IOS_tuperword64,                     &
              IOS_tutype, node%metadata%client,                            &
              Node%payloadTag, Global_Comm,                                &
              ierror)
        END IF
        CALL release_lock_mpi()
      ELSE
        IOS_ReadCompletionRequested=1
!$OMP FLUSH
        IF (IOS_Verbosity>=IOS_PrStatus_Oper) THEN
          CALL IOS_print(                                                  &
              'Waiting for thread 0 to complete a request for read',       &
              src='io_server_writer')
        END IF
        DO WHILE (IOS_ReadCompletionRequested==1)
          CALL um_sleep(IOS_backoff_interval)
!$OMP FLUSH
        END DO
      END IF

    CASE (IOS_Action_Read32 )

      IF (node%metadata%address >= 0 ) THEN
        ! Metadata contains embedded setpos address
        IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
!$OMP CRITICAL(internal_write)
          WRITE(IOS_message,'(A,I12)')                                     &
              'Info: Writer: Process embedded setpos: ',                   &
              node%metadata%address
!$OMP END CRITICAL(internal_write)
          CALL IOS_print(IOS_message,src='io_server_writer')
        END IF
        CALL Setpos( node%metadata%UNIT,                                   &
            node%metadata%address, ierror)
      END IF

      CALL IOS_Increment_Queue_Length                                      &
          (node%metadata%data_size*IOS_BytesPerWord32,needing_lock)
      ALLOCATE(node%real32_data( node %metadata%data_size) )
      CALL Buffin32( node%metadata%UNIT,                                   &
          node%real32_data, node%metadata%data_size,                       &
          len_IO, IO_stat)
      IF (.NOT. IOS_thread_0_calls_mpi) THEN
        CALL acquire_lock_mpi()
        IF (node%metadata%subtype==IOS_Read_Broadcast) THEN
          CALL mpl_bcast(                                                  &
              node%real32_data,                                            &
              node%metadata%data_size*IOS_tuperword32,                     &
              IOS_tutype,                                                  &
              bcast_procs-1,                                               &
              IOS_BCast_Server_Comm,                                       &
              ierror)
        ELSE
          CALL MPL_Send(node%real32_data,                                  &
              node%metadata%data_size*IOS_tuperword32,                     &
              IOS_tutype,                                                  &
              node%metadata%client,                                        &
              Node%payloadTag, Global_Comm, ierror)
        END IF
        CALL release_lock_mpi()
      ELSE
        IOS_ReadCompletionRequested=1
!$OMP FLUSH
        IF (IOS_Verbosity>=IOS_PrStatus_Oper) THEN
          CALL IOS_print(                                                  &
              'Waiting for thread 0 to complete a request for read',       &
              src='io_server_writer')
        END IF
        DO WHILE (IOS_ReadCompletionRequested==1)
          CALL um_sleep(IOS_backoff_interval)
!$OMP FLUSH
        END DO
      END IF

    CASE (IOS_Action_Read64_Integer)

      IF (node%metadata%address >= 0 ) THEN
        ! Metadata contains embedded setpos address
        IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
!$OMP CRITICAL(internal_write)
          WRITE(IOS_message,'(A,I12)')                                     &
              'Info: Writer: Process embedded setpos: ',                   &
              node%metadata%address
!$OMP END CRITICAL(internal_write)
          CALL IOS_print(IOS_message,src='io_server_writer')
        END IF

        CALL Setpos( node%metadata%UNIT,                                   &
            node%metadata%address, ierror)
      END IF

      CALL IOS_Increment_Queue_Length                                      &
          (node%metadata%data_size*IOS_BytesPerWord64,needing_lock)
      ALLOCATE(node%integer_data( node%metadata%data_size) )
      CALL Buffin( node%metadata%UNIT,                                     &
          node%integer_data,                                               &
          node%metadata%data_size,  len_IO, IO_stat )
      IF (.NOT. IOS_thread_0_calls_mpi) THEN
        CALL acquire_lock_mpi()
        IF (node%metadata%subtype==IOS_Read_Broadcast) THEN
          CALL mpl_bcast(                                                  &
              node%integer_data,                                           &
              node%metadata%data_size*IOS_tuperword64,                     &
              IOS_tutype,                                                  &
              bcast_procs-1,                                               &
              IOS_BCast_Server_Comm,                                       &
              ierror)
        ELSE
          CALL MPL_Send(node%integer_data,                                 &
              node%metadata%data_size*IOS_tuperword64,                     &
              IOS_tutype,                                                  &
              node%metadata%client, Node%payloadTag,                       &
              Global_Comm, ierror)
        END IF
        CALL release_lock_mpi()
      ELSE
        IOS_ReadCompletionRequested=1
        IF (IOS_Verbosity>=IOS_PrStatus_Oper) THEN
          CALL IOS_print(                                                  &
              'Waiting for thread 0 to complete a request for read',       &
              src='io_server_writer')
        END IF
!$OMP FLUSH
        DO WHILE (IOS_ReadCompletionRequested==1)
          CALL um_sleep(IOS_backoff_interval)
!$OMP FLUSH
        END DO
      END IF

      !--------------------------------------------------------------
      ! Enquiry
      !--------------------------------------------------------------

    CASE (IOS_Action_Enquire)

      IF (node%metadata%address >= 0 ) THEN
        ! Metadata contains embedded setpos address
        IF (IOS_Verbosity>=IOS_PrStatus_diag) THEN
!$OMP CRITICAL(internal_write)
          WRITE(IOS_message,'(A,I12)')                                     &
              'Info: Writer: Process embedded setpos: ',                   &
              node%metadata%address
!$OMP END CRITICAL(internal_write)
          CALL IOS_print(IOS_message,src='io_server_writer')
        END IF
        CALL Setpos( node%metadata%UNIT,                                   &
            node%metadata%address, ierror)
      END IF

      CALL ioFileState(node%metadata%UNIT,fState)
      IOS_unit_status%extent   = fstate%fileExtent
      IOS_unit_status%POSITION = fstate%filePosition
      IF (.NOT. IOS_thread_0_calls_mpi) THEN
        CALL acquire_lock_mpi()
        CALL MPL_Send(IOS_unit_status,                                     &
            IOS_State_Size*IOS_tuperword64,                                &
            IOS_tutype,                                                    &
            node%metadata%client, Node%payloadTag,                         &
            Global_Comm, ierror)
        CALL release_lock_mpi()
      ELSE
        IOS_EnqCompletionRequested=1
!$OMP FLUSH
        DO WHILE (IOS_EnqCompletionRequested==1)
          CALL um_sleep(IOS_backoff_interval)
!$OMP FLUSH
        END DO
      END IF


      !--------------------------------------------------------------
      ! Open File
      !--------------------------------------------------------------

    CASE (IOS_Action_Open)

      ! Here we add a file object to the IO PEs to correspond to the object
      ! added by the model PEs, and use the FORCE keyword to ensure it is
      ! assigned the correct unit number
      CALL assign_file_unit(node%metadata%string,                          &
          node%metadata%UNIT,                                              &
          force=node%metadata%UNIT, handler="portio")

      IF (node%metadata%address == ioFileTypeUM .OR.                       &
           node%metadata%address == ioFileTypeMPIIO ) THEN

        CALL Model_File_Open(node%metadata%UNIT,                           &
            node%metadata%string,                                          &
            node%metadata%name_length,                                     &
            node%metadata%subtype,                                         &
            ioNameProvided,                                                &
            ierror,                                                        &
            ioLocality=ioAllLocal,                                         &
            fileType=node%metadata%address)

      ELSE

        CALL File_Open(node%metadata%UNIT,                                 &
            node%metadata%string,                                          &
            node%metadata%name_length,                                     &
            node%metadata%subtype,                                         &
            ioNameProvided,                                                &
            ierror,                                                        &
            ioLocality=ioAllLocal,                                         &
            fileType=node%metadata%address)

      END IF

      IF (ierror /= 0) THEN
!$OMP CRITICAL(internal_write)
        WRITE(IOS_Wrt_Message,'(A,I3)')                                    &
            'Error in Opening file on unit ',                              &
            node%metadata%UNIT
!$OMP END CRITICAL(internal_write)
        CALL IOS_Ereport( RoutineName, 10, IOS_Wrt_Message,          &
            md=node%metadata )
      END IF

      !--------------------------------------------------------------
      ! Close File
      !--------------------------------------------------------------

    CASE (IOS_Action_Close)

      IF (is_unit_open(node%metadata%UNIT)) THEN! Unit is Open
        NULLIFY(um_file)
        um_file => get_file_by_unit(node%metadata%UNIT, handler="portio")
        IF (um_file % pp_meta % managed) THEN
          CALL Model_File_Close(node%metadata%UNIT,                        &
              node%metadata%string,                                        &
              node%metadata%name_length,                                   &
              1,                                                           &
              node%metadata%delete,                                        &
              ierror )
        ELSE
          CALL File_Close(node%metadata%UNIT,                              &
              node%metadata%string,                                        &
              node%metadata%name_length,                                   &
              1,                                                           &
              node%metadata%delete,                                        &
              ierror )
        END IF

        IF (ierror /= 0) THEN
!$OMP CRITICAL(internal_write)
          WRITE(IOS_Wrt_Message,'(A,I3)')                                  &
              'Error in Closing File on unit ',                            &
              node%metadata%UNIT
!$OMP END CRITICAL(internal_write)
          CALL IOS_Ereport( RoutineName, 10, IOS_Wrt_Message,        &
              md=node%metadata )
        END IF

        CALL release_file_unit(node%metadata%UNIT, handler="portio")

      ELSE
        IF (IOS_Verbosity>=IOS_PrStatus_Normal) THEN
!$OMP CRITICAL(internal_write)
          WRITE(IOS_message,'(A,I3,A)')                                    &
              'Tried to close unit (',node%metadata%UNIT,                  &
              ') which is not open'
!$OMP END CRITICAL(internal_write)
          CALL IOS_print(IOS_message,src='io_server_writer')
        END IF
      END IF

      !--------------------------------------------------------------
      ! Setpos
      !--------------------------------------------------------------

    CASE (IOS_Action_Setpos)
      IF (Is_unit_open(node%metadata%UNIT)) THEN! Unit is Open
        CALL Setpos( node%metadata%UNIT,                                   &
            node%metadata%address, ierror)
        IF (ierror /= 0) THEN
!$OMP CRITICAL(internal_write)
          WRITE(IOS_Wrt_Message,'(A,I3)')                                  &
              'Error in setpos on unit ',                                  &
              node%metadata%UNIT
!$OMP END CRITICAL(internal_write)
          CALL IOS_Ereport( RoutineName, 10, IOS_Wrt_Message,        &
              md=node%metadata )
        END IF
      ELSE
!$OMP CRITICAL(internal_write)
        WRITE(IOS_Wrt_Message,'(A,I3)')                                    &
            'Setpos called on non-open unit:',                             &
            node%metadata%UNIT
!$OMP END CRITICAL(internal_write)
        CALL IOS_Ereport( RoutineName, -999, IOS_Wrt_Message,        &
            md=node%metadata )
      END IF

      !-------------------------------------------------------------------
      ! sync Ops
      !-------------------------------------------------------------------

    CASE (IOS_Action_Sync)
       IF ( IOS_enable_mpiio .AND.                                         &
            node%metadata%address == ioFileTypeMPIIO) THEN
          CALL acquire_lock_mpi()
          CALL MPL_File_sync( mpiio_fh(node%metadata%UNIT), ierror )
          CALL release_lock_mpi()
       ELSE
          CALL Sync_single(node%metadata%UNIT,ierror)
       END IF

    CASE (IOS_Action_Sync_Barrier)
       IF ( IOS_enable_mpiio .AND.                                         &
            node%metadata%address == ioFileTypeMPIIO ) THEN
          CALL acquire_lock_mpi()
          CALL MPL_File_sync( mpiio_fh(node%metadata%UNIT), ierror )
          CALL release_lock_mpi()
       ELSE
          CALL Sync_single(node%metadata%UNIT,ierror)
       END IF
       ! Check on model_rank necessary to prevent deadlock when
       ! IOS_enable_mpiio is true
       IF (.NOT. IOS_thread_0_calls_mpi .AND. model_rank == 0 ) THEN
        ! Send 1 integer to the client to indicate completion
          CALL MPL_Send(node%metadata,                                     &
            1 , mpl_integer, node%metadata%client,                         &
            Node%payloadTag, Global_comm,                                  &
            ierror)
       END IF

      !--------------------------------------------------------------
      ! File operations
      !--------------------------------------------------------------

    CASE (IOS_Action_FileOp)
      t=get_wallclock_time()
      IF (IOS_Verbosity>=IOS_PrStatus_Oper) THEN
!$OMP CRITICAL(internal_write)
        WRITE(IOS_message,'(A,F10.3)')                                     &
            'Info: Writer: FileOp, syncing at time=',                      &
            t-IOS_Start_time
!$OMP END CRITICAL(internal_write)
        CALL IOS_print(IOS_message,src='io_server_writer')
      END IF
      IF (.NOT. IOS_thread_0_calls_mpi) THEN
        CALL IOS_SoftSync(io_server_for_unit(file_op_pseudo_unit))
      END IF
      t=get_wallclock_time()
      IF (IOS_Verbosity>=IOS_PrStatus_Oper) THEN
!$OMP CRITICAL(internal_write)
        WRITE(IOS_message,'(A,F10.3)')                                     &
            'Info: Writer: FileOp, syncing done at time=',                 &
            t-IOS_Start_time
!$OMP END CRITICAL(internal_write)
        CALL IOS_print(IOS_message,src='io_server_writer')
      END IF
      IF (global_rank==io_server_for_unit(file_op_pseudo_unit)) THEN
        CALL file_action(TRIM(node%metadata%string))
      END IF

      !--------------------------------------------------------------
      ! Fence Ops
      !--------------------------------------------------------------

    CASE (IOS_Action_Fence)

      IF (.NOT. IOS_thread_0_calls_mpi) THEN

        t=get_wallclock_time()
        IF (IOS_Verbosity>=IOS_PrStatus_Oper) THEN
!$OMP CRITICAL(internal_write)
          WRITE(IOS_message,'(A,I4,A,F10.3)')                              &
              'Info: Writer: Fence: server=',                              &
              node%metadata%address,' at time=',                           &
              t-IOS_Start_time
!$OMP END CRITICAL(internal_write)
          CALL IOS_print(IOS_message,src='io_server_writer')
        END IF
        CALL IOS_SoftSync(node%metadata%address)

        t=get_wallclock_time()
        IF (IOS_Verbosity>=IOS_PrStatus_Oper) THEN
!$OMP CRITICAL(internal_write)
          WRITE(IOS_message,'(A,F10.3)')                                   &
              'Info: Writer: Fence completed at time=',                    &
              t-IOS_Start_time
!$OMP END CRITICAL(internal_write)
          CALL IOS_print(IOS_message,src='io_server_writer')
        END IF
      ELSE
        IF (IOS_Verbosity>=IOS_PrStatus_Normal) THEN
!$OMP CRITICAL(internal_write)
          WRITE(IOS_message,'(A,A)')'Info: Writer:',                       &
              ' Fence occured previously due to MPI semantics'
!$OMP END CRITICAL(internal_write)
          CALL IOS_print(IOS_message,src='io_server_writer')
        END IF
      END IF

      !--------------------------------------------------------------
      ! Finish - close down server
      !--------------------------------------------------------------

    CASE (IOS_Action_Finish)

      t=get_wallclock_time()
      CALL IOS_print                                                       &
          ('********************************************************',     &
          src='io_server_writer')
!$OMP CRITICAL(internal_write)
      WRITE(IOS_message,'(A,F10.3,A)')                                     &
          '* IO SERVER:WRITER   RECEIVED FINISH CMD AT TIME=',             &
          t-IOS_Start_Time,' *'
!$OMP END CRITICAL(internal_write)
      CALL IOS_print(IOS_message,src='io_server_writer')
      CALL IOS_print                                                       &
          ('********************************************************',     &
          src='io_server_writer')
      done=.TRUE.

      !-------------------------------------------------------------------
      ! Process a setpos. For Stash setpos we are given a
      ! unit/record from which to compute a seek address.
      !-------------------------------------------------------------------

    CASE (IOS_Action_StashSetPos)

      IF ( Is_unit_open( node%metadata % UNIT ) ) THEN
        fixed_header=>attachHeader(node%metadata%UNIT,FixedHeader)
        ipplook=>attachLookups(node%metadata%UNIT)
        IF (node%metadata%address==1) THEN
          node%metadata%address=Fixed_Header(160)-1
        ELSE
          node%metadata%address=                                           &
              ipplook(lbegin, node%metadata%address-1)+                    &
              ipplook(lbnrec, node%metadata%address-1)
        END IF

        CALL Setpos( node%metadata % UNIT,                                 &
            node%metadata % address, ierror)
        IF (ierror /= 0) THEN
!$OMP CRITICAL(internal_write)
          WRITE(IOS_Wrt_Message,'(A,I8)') 'Error in setpos on unit ',      &
              node%metadata % UNIT
!$OMP END CRITICAL(internal_write)
          CALL IOS_Ereport( RoutineName, 50, IOS_Wrt_Message,        &
              md=node%metadata )
        END IF
      ELSE
!$OMP CRITICAL(internal_write)
        WRITE(IOS_Wrt_Message,'(A,I8)') 'Setpos called on non-open unit:', &
            node%metadata % UNIT
!$OMP END CRITICAL(internal_write)
        CALL IOS_Ereport( RoutineName, -999, IOS_Wrt_Message,        &
            md=node%metadata )
      END IF
      NULLIFY(ipplook)
      NULLIFY(fixed_header)

      !-------------------------------------------------------------------
      ! Process The merging of PP lookup tables
      !-------------------------------------------------------------------

    CASE (IOS_action_mergepplookup)

      IF (.NOT. ASSOCIATED(node%integer_data)) THEN
        WRITE(IOS_wrt_message,'(A)')                                       &
            'FAULT: there is no data on this (Action_WritePPLook) node'
        CALL IOS_Ereport( RoutineName, 10, IOS_Wrt_Message,          &
            md=node%metadata )
      END IF

      CALL setLookups(node%metadata%UNIT, node%integer_data)
      IF (IOS_Verbosity>=IOS_PrStatus_Oper) THEN
!$OMP CRITICAL(internal_write)
        WRITE(IOS_Message,'(A,A,I3)')'Writer: ',                           &
            'Merged lookup records for unit ',                             &
            node%metadata % UNIT
!$OMP END CRITICAL(internal_write)
        CALL IOS_print(IOS_message,src='io_server_writer')
      END IF

      !-------------------------------------------------------------------
      ! Process The writing of PP lookup tables
      !-------------------------------------------------------------------

    CASE (IOS_ACTION_StashWritePPLookup)
      IF (.NOT. ASSOCIATED(node%integer_data)) THEN
        WRITE(IOS_wrt_message,'(A)')                                       &
            'FAULT: there is no data on this (Action_WritePPLook) node'
        CALL IOS_Ereport( RoutineName, 10, IOS_Wrt_Message,          &
            md=node%metadata )
      END IF

      IF (node%metadata%address >= 0 ) THEN
        ! Metadata contains a setpos address
        IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
!$OMP CRITICAL(internal_write)
          WRITE(IOS_message,'(A,I12)')                                     &
              'Info: Writer: Process embedded setpos: ',                   &
              node%metadata%address
!$OMP END CRITICAL(internal_write)
          CALL IOS_print(IOS_message,src='io_server_writer')
        END IF

        CALL Setpos( node%metadata%UNIT,                                   &
            node%metadata%address, ierror)
      END IF

      IO_Stat = -1.0
      Len_IO  = node%metadata % data_size

      CALL setLookups(node%metadata % UNIT,node%integer_data)
      IF (IOS_Verbosity>=IOS_PrStatus_Oper) THEN
!$OMP CRITICAL(internal_write)
        WRITE(ios_message,'(A,A,I3)')'Writer: ',                           &
            'Merged lookup records for unit ',                             &
            node%metadata % UNIT
!$OMP END CRITICAL(internal_write)
        CALL ios_print(ios_message,src='io_server_writer')
      END IF
      ipplook=>attachLookups(node%metadata % UNIT)

      CALL Buffout(node%metadata % UNIT,ipplook,                           &
          node%metadata % data_size, Len_IO, IO_Stat )

      IF (Len_IO /= node%metadata % data_size .OR.                         &
          IO_Stat /= -1.0) THEN
!$OMP CRITICAL(internal_write)
        WRITE(IOS_Wrt_Message,'(A,I4,I10,I10,F6.2)')                       &
            'Error in Writing to file on unit ',                           &
            node%metadata % UNIT, Len_IO,                                  &
            node%metadata % data_size,                                     &
            IO_Stat
!$OMP END CRITICAL(internal_write)
        CALL IOS_Ereport( RoutineName, 10, IOS_Wrt_Message,          &
            md=node%metadata )
      END IF

      NULLIFY(ipplook)

      !-------------------------------------------------------------------
      ! Process The initialisation of PP lookup buffers
      !-------------------------------------------------------------------

    CASE (IOS_Action_StashInitPPLookup)
      IF (.NOT. ASSOCIATED(node%integer_data)) THEN
        WRITE(IOS_wrt_message,'(A)')                                       &
            'FAULT: no data in Action_StashInitPPLookup node'
        CALL IOS_Ereport( RoutineName, 10, IOS_Wrt_Message,          &
            md=node%metadata )
      END IF

      CALL InitLookups(ipplook,                                            &
          node%integer_data(1),                                            &
          node%integer_data(2),                                            &
          node%integer_data(3),                                            &
          node%integer_data(4),                                            &
          node%integer_data(5))
      NULLIFY(ipplook)

      !-------------------------------------------------------------------
      ! Process The initialisation of a fixed header
      !-------------------------------------------------------------------

    CASE (IOS_Action_StashInitHeader)
      CALL initHeader(node%metadata % UNIT, FixedHeader,                   &
          node%metadata%address)

      !-------------------------------------------------------------------
      ! Process The setting of a fixed header
      !-------------------------------------------------------------------

    CASE (IOS_Action_StashSetHeader)
      IF (.NOT. ASSOCIATED(node%integer_data)) THEN
        WRITE(IOS_wrt_message,'(A)')                                       &
            'FAULT: no data in Action_StashSetFXH node'
        CALL IOS_Ereport( RoutineName, 10, IOS_Wrt_Message,          &
            md=node%metadata )
      END IF
      CALL setHeader(node%metadata % UNIT, FixedHeader, node%integer_data)

      !--------------------------------------------------------------
      ! Dump initialisation
      !--------------------------------------------------------------

    CASE (IOS_Action_DumpInitModel)
      CALL IOS_DumpInitModel(                                              &
          node%integer_data,                                               &
          node%metadata%data_size )

      !--------------------------------------------------------------
      ! Write out pp data
      !--------------------------------------------------------------

    CASE (IOS_Action_StashWritePPData)
      CALL IOS_Stash_Server_Wait_For_Data(node)
      CALL IOS_Stash_server_process(node)
      CALL IOS_Stash_server_deallocate(node)

      !--------------------------------------------------------------
      ! Write out dump data
      !--------------------------------------------------------------

    CASE (IOS_Action_StashWriteDumpData)
      CALL IOS_Stash_Server_Wait_For_Data(node)
      CALL IOS_Stash_server_process(node)
      CALL IOS_Stash_server_deallocate(node)

      !--------------------------------------------------------------
      ! Modify IOS behaviour
      !--------------------------------------------------------------

    CASE (IOS_Action_Config)

      CALL IOS_print('IOS_Action_Config: not currently used',              &
          src='io_server_writer')

      !--------------------------------------------------------------
      ! Unrecognised action
      !--------------------------------------------------------------

    CASE DEFAULT
!$OMP CRITICAL(internal_write)
      WRITE(IOS_Wrt_Message,'(A,I8)')                                      &
          'Action not recognised:',                                        &
          node%metadata%ACTION
!$OMP END CRITICAL(internal_write)
      CALL IOS_Ereport( RoutineName, 10, IOS_Wrt_Message,            &
          md=node%metadata )
    END SELECT

    IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
      t2=get_wallclock_time()
!$OMP CRITICAL(internal_write)
      WRITE(IOS_message,'(A,I8,A,F8.3)')                                   &
          'Info: Writer: transaction: ',                                   &
          node%transaction_number,' is completed in ',                     &
          t2-t1
!$OMP END CRITICAL(internal_write)
      CALL IOS_print(IOS_message,src='io_server_writer')
      CALL IOS_print_flush()
    END IF
    CALL IOS_Remove_Data_From_Queue()
  END IF
END DO

CALL IOS_print( 'Info: Writer: Process closing',src='io_server_writer')
CALL IOS_print_flush()
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE IOS_Writer
END MODULE IO_Server_Writer
