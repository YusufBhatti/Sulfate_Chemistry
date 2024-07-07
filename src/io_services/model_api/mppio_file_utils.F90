! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: C96

MODULE MPPIO_file_utils

  ! Implement basic file operations. File ops originated through this
  ! module will execute immediately if no IO servers are present or will be
  ! queue on an IO server with fence semantics if IO servers are active. Fence
  ! semantics mandate that the operation will take place in-order with respect
  ! to any other IO server queued operations and the operation will be delayed
  ! until pending IO ops are completed.
  !
  ! Commands should normally be issued via high level interfaces;
  !    file_copy(<src>,<dest>)
  !    file_delete(<file>)
  !    create_directory(<name>)
  !
  ! which will be translated into a cannonical form and passed to
  ! the generic function;
  !    file_action()
  ! which in turn calls an implementation in c (portio).
  !    <error code>=file_op(<cannonical string>)
  !
  ! The optional force_local argument specifies that (if .TRUE.) the
  ! operation will be perfomed locally and immediately if the calling rank is
  ! rank 0 and the op can be assumed to have been completed on rank 0 upon
  ! return. Other ranks can not make this assumption.
  !
  ! An explicit initialisation on all ranks is required before use.

USE io
USE IOS_Client_Queue
USE ios_common, ONLY: file_op_pseudo_unit, L_IO_Server, L_IOS_active

USE yomhook, ONLY:                                                          &
    lhook,                                                                   &
    dr_hook
USE parkind1, ONLY:                                                         &
    jprb,                                                                    &
    jpim

USE um_parvars
USE ereport_mod
USE umPrintMgr

USE IOS_Constants, ONLY: IOS_Action_FileOp

IMPLICIT NONE

! Profiling
INTEGER(KIND=jpim),PRIVATE, PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim),PRIVATE, PARAMETER :: zhook_out = 1

! Seperator character in cannonical command strings
CHARACTER(LEN=*), PARAMETER             :: sep=":"

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='MPPIO_FILE_UTILS'

CONTAINS

SUBROUTINE mppio_file_utils_init()
IMPLICIT NONE
REAL(KIND=jprb):: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MPPIO_FILE_UTILS_INIT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Ensure a mapping exists for file ops.
! Although these do not use a unit, we still need to define
! which IOS will handle them.
IF (L_IOS_active() .AND. .NOT. L_IO_Server) THEN
  IF (io_server_for_unit(file_op_pseudo_unit)==IOS_No_Server)              &
      CALL IOS_assign_server_for_unit(file_op_pseudo_unit)
END IF

IF (printstatus>=prstatus_normal) THEN
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I3)')'MPPIO_File_Utils: Initialised file utils using unit ', &
      file_op_pseudo_unit
!$OMP END CRITICAL(internal_write)
  CALL umPrint(umMessage,src='mppio_file_utils')
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE mppio_file_utils_init

SUBROUTINE file_touch(theFile,force_local)
IMPLICIT NONE
CHARACTER (LEN=*), INTENT(IN)  :: theFile  ! to read
LOGICAL, OPTIONAL              :: force_local
CHARACTER (LEN=*), PARAMETER   :: cmd="touch"
CHARACTER (LEN=IOS_String_Max) :: lineBuffer
INTEGER                        :: err_code
REAL(KIND=jprb)                :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FILE_TOUCH'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Check string is small enough to transmit.
IF (LEN_TRIM(theFile) +LEN_TRIM(cmd)+ 3*LEN_TRIM(sep)                      &
    > IOS_String_Max) THEN
  WRITE(umMessage,'(A,A)')'String too long theFile =',theFile
  CALL umPrint(umMessage,src='mppio_file_utils')
  err_code = 99
  CALL ereport('file_utils:file_delete',err_code,                          &
      'Arguments too long')
END IF

WRITE(lineBuffer,'(A,A,A,A,A)')                                            &
    TRIM(cmd),                                                             &
    TRIM(sep),                                                             &
    TRIM(sep),                                                             &
    TRIM(theFile),                                                         &
    TRIM(sep)

CALL file_action(TRIM(lineBuffer),force_local)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE file_touch

SUBROUTINE file_copy(src,dest,force_local)
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN)    :: src  ! to read
CHARACTER(LEN=*),INTENT(IN)    :: dest ! to write
LOGICAL, OPTIONAL            :: force_local
INTEGER                      :: err_code
CHARACTER(LEN=*), PARAMETER    :: cmd="cp"
CHARACTER (LEN=IOS_String_Max)  :: lineBuffer
REAL(KIND=jprb)              :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FILE_COPY'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Check string is small enough to transmit.
IF (LEN_TRIM(src)+LEN_TRIM(dest)+LEN_TRIM(cmd)+3*LEN_TRIM(sep)             &
    > IOS_String_Max) THEN
  WRITE(umMessage,'(A,A)')'Str too long src =',src
  CALL umPrint(umMessage,src='mppio_file_utils')
  WRITE(umMessage,'(A,A)')'             dest=',dest
  CALL umPrint(umMessage,src='mppio_file_utils')
  err_code = 99
  CALL ereport('file_utils:file_copy',err_code,                            &
      'Arguments too long')
END IF

WRITE(lineBuffer,'(A,A,A,A,A,A)')                                          &
    TRIM(cmd),                                                             &
    TRIM(sep),                                                             &
    TRIM(sep),                                                             &
    TRIM(src),                                                             &
    TRIM(sep),                                                             &
    TRIM(dest)

CALL file_action(TRIM(lineBuffer),force_local)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE file_copy

SUBROUTINE dir_create(theDirectory,force_local)
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN)    :: theDirectory
LOGICAL, OPTIONAL            :: force_local
INTEGER                      :: err_code
CHARACTER(LEN=*), PARAMETER    :: cmd="mkdir"
CHARACTER (LEN=IOS_String_Max)  :: lineBuffer
REAL(KIND=jprb)              :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DIR_CREATE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Check string is small enough to transmit.
IF (LEN_TRIM(theDirectory) +LEN_TRIM(cmd)+ 3*LEN_TRIM(sep)                 &
    > IOS_String_Max) THEN
  WRITE(umMessage,'(A,A)')'Str too long theDirectory =',theDirectory
  CALL umPrint(umMessage,src='mppio_file_utils')
  err_code = 99
  CALL ereport('file_utils:file_delete',err_code,                          &
      'Arguments too long')
END IF

WRITE(lineBuffer,'(A,A,A,A,A)')                                            &
    TRIM(cmd),                                                             &
    TRIM(sep),                                                             &
    TRIM(sep),                                                             &
    TRIM(theDirectory),                                                    &
    TRIM(sep)

CALL file_action(TRIM(lineBuffer),force_local)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE dir_create

SUBROUTINE file_delete(theFile,force_local,silent)
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN)    :: theFile  ! to delete
LOGICAL, OPTIONAL              :: force_local
LOGICAL, OPTIONAL              :: silent
LOGICAL                        :: lsilent
INTEGER                        :: err_code
CHARACTER(LEN=*), PARAMETER    :: cmd="rm"
CHARACTER(LEN=3)               :: option
CHARACTER (LEN=IOS_String_Max) :: lineBuffer
REAL(KIND=jprb)                :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FILE_DELETE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Check string is small enough to transmit.
IF (LEN_TRIM(theFile) + LEN_TRIM(cmd)+ 3*LEN_TRIM(sep)                     &
    > IOS_String_Max) THEN
  WRITE(umMessage,'(A,A)')'Str too long theFile =',theFile
  CALL umPrint(umMessage,src='mppio_file_utils')
  err_code = 99
  CALL ereport('file_utils:file_delete',err_code,                          &
      'Arguments too long')
END IF

lsilent=.FALSE.
IF (PRESENT(silent)) lsilent=silent

!$OMP CRITICAL(internal_write)
WRITE(option,'(A)')''
IF (lsilent) WRITE(option,'(A)')'-f'

WRITE(lineBuffer,'(A,A,A,A,A,A)')                                          &
    TRIM(cmd),                                                             &
    TRIM(sep),                                                             &
    TRIM(option),                                                          &
    TRIM(sep),                                                             &
    TRIM(theFile),                                                         &
    TRIM(sep)
!$OMP END CRITICAL(internal_write)

CALL file_action(TRIM(lineBuffer),force_local)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE file_delete

SUBROUTINE file_action(command,force_local)
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN)   :: command
LOGICAL, OPTIONAL            :: force_local
LOGICAL                      :: l_force_local ! localised copy of above
INTEGER                      :: pe,i          ! loop counters
INTEGER                      :: error         ! return code from c
INTEGER                      :: qHandle       ! IOS operation handle
INTEGER                      :: err_code
INTEGER, EXTERNAL            :: file_op       ! c method to do the work
TYPE(IOS_metadata_type),                                                   &
    POINTER                  :: metadata
REAL(KIND=jprb)              :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FILE_ACTION'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

l_force_local=.FALSE.
IF (PRESENT(force_local)) l_force_local=force_local

IF (L_IOS_active()    .AND.                                                &
    .NOT. L_IO_Server .AND.                                                &
    .NOT. l_force_local) THEN

  IF (io_server_for_unit(file_op_pseudo_unit) ==                           &
      IOS_No_Server) THEN
    err_code = 99
    CALL ereport('file_utils:file_action',err_code,                       &
        'Module used before initialisation')
  END IF

  IF (model_rank == 0) THEN

    DO pe=1,IOS_Server_Groups
      qHandle=ios_init_md(-1*io_servers(pe,1),                             &
          IOS_No_location,                                                 &
          IOS_Action_FileOp)
      metadata => IOS_Attach_Metadata(qHandle)

      ! Command
      DO i=1,IOS_String_MAX
        metadata % string(i:i)=' '
      END DO
      WRITE(metadata % string,'(A)')command
      metadata % name_length=LEN_TRIM(command)
      metadata % subtype = -1
      CALL IOS_Send(qHandle,hasString=.TRUE.)
    END DO
  END IF

ELSE
  IF (model_rank == 0) THEN
    WRITE(umMessage,'(A,A)')'MPPIO: file op: ',TRIM(command)
    CALL umPrint(umMessage,src='mppio_file_utils')
    error=file_op(command,LEN_TRIM(command))
    IF (error /= 0) THEN
      CALL ereport('file_utils:file_action',error,                         &
          'A file action failed in the C level implementation,'//          &
          ' see prior messages')
    ELSE
      CALL umPrint('MPPIO: file op completed',src='mppio_file_utils')
    END IF
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE file_action

END MODULE MPPIO_file_utils
