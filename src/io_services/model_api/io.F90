! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: C96

MODULE io
! A module to contain Fortran interfaces to file operations

! Description:
!  This module provides BUFFOUT,BUFFIN,OPEN,CLOSE,SETPOS routines
!  for the Parallel Unified Model, Reconfiguration, and Small Executables.
!
! Complete list of interfaces
!
!             buffin,buffin64,buffin32
!             buffout,buffout64,buffout32
!
!  Note that routines are selected by array type, and hence the generic
!  (non-lengthed) interface is preferred. The array type supplied MUST
!  reflect its content in order to facilitate correct byteswapping.
!
! Complete list of routines
!
! ioInit()
! ioShutdown()
! checkUnit(unit)
! file_open
!  (unit,env,env_len,read_write,name_in_environ,error,ioLocality)
! file_close(unit,env,env_len,name_in_environ,delete,err)
! buffout64_r(unit,array,numWords,wordsWritten,errorCode)
! buffout64_i(unit,array,numWords,wordsWritten,errorCode)
! buffout64_l(unit,array,numWords,wordsWritten,errorCode)
! buffout64_r2D(unit,array,numWords,wordsWritten,errorCode)
! buffout64_i2D(unit,array,numWords,wordsWritten,errorCode)
! buffout64_l2D(unit,array,numWords,wordsWritten,errorCode)
! buffin64_r(unit,array,numWords,wordsRead,errorCode)
! buffin64_i(unit,array,numWords,wordsRead,errorCode)
! buffin64_l(unit,array,numWords,wordsRead,errorCode)
! buffin64_i2D(unit,array,numWords,wordsRead,errorCode)
! buffin64_l2D(unit,array,numWords,wordsRead,errorCode)
! buffin64_r2D(unit,array,numWords,wordsRead,errorCode)
! buffin32_r(unit,array,numWords,wordsRead,errorCode)
! buffin32_i(unit,array,numWords,wordsRead,errorCode)
! buffin32_i2D(unit,array,numWords,wordsRead,errorCode)
! buffin32_l(unit,array,numWords,wordsRead,errorCode)
! buffin_character(unit,array,numWords,csize,wordsRead,errorCode)
! setpos(unit,ipos,icode)
! setpos8(unit,ipos,icode)
! setpos32(unit,ipos,icode)
! ioFileState(unit,state)
! readbleWords(unit,wordLength)
! set_unit_bcast_flag(unit)
! clear_unit_bcast_flag(unit)
! logical broadcast_read(unit)
!
! PRIVATE:
!       setAttribute(var,mask)
!       clearAttribute(var,mask)
!       LOGICAL testAttribute(var,mask)
!       performOperation(unit)
!       io_ereport
!
!  Note that setpos expresses units of the default integer/real type
!  and there is currently no explicit setpos64, and that 'getpos'
!  functionality is not exposed at the fortran layer - use ioFileState
!  instead to obtain the current file position in bytes.
!
! Method:
!  The C interfaces *_single are called only by the IO module
!  appropriately on tasks designated to perform IO. Calls to C
!  versions should not be made elsewhere. The designation of who
!  performs IO is specified at file open time via the ioLocality
!  optional argument and is persistant until file close (with the
!  exception of buffin, for which broadcasting can be switched on and off)
!
!  Failing to provide the ioLocality argument results in a default
!  behaviour in which IO Servers perform IO if available for writeable
!  files, otherwise rank zero. Rank zero will perform IO for non-writeable
!  files (unless overriden via the IO Server configuration). By default
!  buffin operations are collective and broadcast read data to all tasks.
!
!  Specifying ioLocality=ioAllLocal results in each task performing IO
!  completely independently (take care to write to different files).
!
!  When using any option other than ioAllLocal, calls should be regarded
!  as collective and the call made on all tasks, with correct arguments
!  (e.g. unit numbers, modes etc) even if you think the task
!  will not actually participate. This allows for future optimisations
!  to perform in parallel and/or distribute IO load.
!
!  Some exceptions are made to the above. They are there to support the
!  existing UM code and such use is depricated. Modification to existing UM IO
!  code should attempt to make it collective.

USE um_types
USE filenamelength_mod, ONLY:                                               &
    maxFileNameLength => fileNameLength
USE ios
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr
USE UM_ParCore, ONLY: mype, nproc
USE IO_dependencies ! Lower layer dependencies
USE fort2c_buffo_interfaces, ONLY: buffout64_single, buffout32_single
USE fort2c_buffin_interfaces, ONLY:                                         &
  buffin64_single, buffin32_single, buffin8_single
USE fort2c_portio_interfaces, ONLY: open_single, close_single, sync_single,    &
                                    getextent, portioinit, portioShutdown

USE fort2c_getpos_interfaces, ONLY:                                         &
  getpos8

USE io_constants, ONLY: ioFileTypeMPIIO, ioOpenWriteOnly, ioOpenReadOnly,   &
                        ioNameinEnv, ioNameProvided, ioReadReplicate,       &
                        ioNoFileName, ioPolicyDefault, ioPolicyLax,         &
                        ioNotKnown, ioDelete, ioPolicySloppy, ioNotSet,     &
                        ioPolicyStrict, ioWriteOnly, ioNoEnvName,           &
                        ioNoLocation, ioFileTypeUnknown, ioReadOnly,        &
                        ioRemote, ioDefault, ioOpen, ioLocal, ioAllLocal,   &
                        ioOpenReadWrite, ioNoDelete

! MPI-IO
USE fort2c_pio_byteswap_interfaces, ONLY: get_machine_endianism,            &
     littleEndian, pio_byteswap
USE mpl, ONLY: MPL_STATUS_SIZE, MPL_INFO_NULL, MPL_SUCCESS, MPL_FILE_NULL,  &
     MPL_OFFSET_KIND, MPL_SEEK_SET, MPL_MODE_WRONLY, MPL_MODE_RDONLY,       &
     MPL_MODE_RDWR, MPL_MODE_CREATE, MPL_ERRORS_ARE_FATAL, MPL_BYTE,        &
     MPL_INTEGER4
USE IOS_Comms, ONLY: acquire_lock_mpi, release_lock_mpi
USE lustre_control_mod, ONLY: apply_file_striping_mpiio
USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_LOC

IMPLICIT NONE
PRIVATE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MPI-IO file handler and info object storage
INTEGER, PUBLIC :: mpiio_fh(minUnit:maxUnit)   = MPL_FILE_NULL
INTEGER, PUBLIC :: mpiio_info(minUnit:maxUnit) = MPL_INFO_NULL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Private storage for the attributes.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Type for containing file/unit state
TYPE unitAttributes
  INTEGER :: state
  INTEGER :: policy
  INTEGER :: nextLocation
  INTEGER :: fileType
  LOGICAL :: from_environment
  CHARACTER (LEN=maxFileNameLength) :: fileName
  CHARACTER (LEN=maxFileNameLength) :: env
END TYPE unitAttributes

! Type for reporting the state of a file, e.g. for enquiry functions
TYPE, PUBLIC :: fileState
  SEQUENCE
  INTEGER :: UNIT
  INTEGER :: fileExtent
  INTEGER :: filePosition
END TYPE fileState

! Length of the above type for transmitting it via e.g. MPI
INTEGER, PARAMETER :: fileState_len=3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Intefaces and access
!
! Nomenclature for implementations:
!
! buff[in|out]<wordsize>_<array_type><rank>

!read numWords words (of some type) from disk
INTERFACE buffin
MODULE PROCEDURE                                                           &
    buffin64_r,                                                            &
    buffin64_i,                                                            &
    buffin64_l,                                                            &
    buffin64_r2D,                                                          &
    buffin64_i2D,                                                          &
    buffin64_l2D,                                                          &
    buffin_character,                                                      &
    buffin32_r,                                                            &
    buffin32_i,                                                            &
    buffin32_i2D,                                                          &
    buffin32_l
END INTERFACE

!read numWords 64 bit words from disk
INTERFACE buffin64
MODULE PROCEDURE                                                           &
    buffin64_r,                                                            &
    buffin64_i,                                                            &
    buffin64_l,                                                            &
    buffin64_r2D,                                                          &
    buffin64_i2D,                                                          &
    buffin64_l2D,                                                          &
    buffin_character
END INTERFACE

!read numWords 32 bit words from disk
INTERFACE buffin32
MODULE PROCEDURE                                                           &
    buffin32_r,                                                            &
    buffin32_i,                                                            &
    buffin32_i2D,                                                          &
    buffin32_l
END INTERFACE

!write numWords words (of some type) to disk
INTERFACE buffout
MODULE PROCEDURE                                                           &
    buffout64_r,                                                           &
    buffout64_i,                                                           &
    buffout64_l,                                                           &
    buffout64_r2D,                                                         &
    buffout64_i2D,                                                         &
    buffout64_l2D,                                                         &
    buffout32_i,                                                           &
    buffout32_r,                                                           &
    buffout32_i2D
END INTERFACE

!write numWords 64 bit words to disk
INTERFACE buffout64
MODULE PROCEDURE                                                           &
    buffout64_r,                                                           &
    buffout64_i,                                                           &
    buffout64_l,                                                           &
    buffout64_r2D,                                                         &
    buffout64_i2D,                                                         &
    buffout64_l2D
END INTERFACE

!write numWords 32 bit words to disk
INTERFACE buffout32
MODULE PROCEDURE                                                           &
    buffout32_i,                                                           &
    buffout32_r,                                                           &
    buffout32_i2D
END INTERFACE

! Access control
PRIVATE performOperation
PRIVATE setAttribute
PRIVATE clearAttribute
PRIVATE testAttribute
PUBLIC  ioInit
PUBLIC  ioShutdown
PUBLIC  io_ereport
PUBLIC  ioDiskSynchronise
PUBLIC  ioFence
PUBLIC  setpos
PUBLIC  setpos8
PUBLIC  setpos32
PUBLIC  set_unit_bcast_flag
PUBLIC  clear_unit_bcast_flag
PUBLIC  buffin
PUBLIC  buffout
PUBLIC  buffout32
PUBLIC  buffin32
PUBLIC  buffout64
PUBLIC  buffin64
PUBLIC  is_unit_open
PUBLIC  isRemote
PUBLIC  file_open
PUBLIC  file_close
PUBLIC  io_timestep
PUBLIC  ioFileState
PUBLIC  readableWords
PUBLIC  ReturnFileName
PUBLIC  isReadOnly
PUBLIC  broadcast_read
PUBLIC  isAllLocal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Global variables (access is not thread safe)

! File States
TYPE(unitAttributes)               :: fileAttributes(minUnit:maxUnit)

! General purpose message
CHARACTER (LEN=2*maxFileNameLength):: io_message

! params/vars  for dr_hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='IO'

CONTAINS


!  Top level initialisation,
!     load config
!     initialise C
!     initialise this module
!
SUBROUTINE ioInit()
USE io_configuration_mod
USE ios_common, ONLY: file_op_pseudo_unit
USE file_manager, ONLY: assign_file_unit
USE fort2c_interfaces, ONLY: logical_to_bool

IMPLICIT NONE
LOGICAL, SAVE  :: initted=.FALSE.
INTEGER        :: UNIT
REAL(KIND=jprb):: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='IOINIT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (.NOT. initted) THEN

  CALL assign_file_unit("pseudo-file for UNIX operations", &
      file_op_pseudo_unit, handler="portio", id="io_reserved_unit")

  CALL ioLoadConfig()
  CALL portioInit(io_wbuffer_size,     &
                  io_rbuffer_size,     &
                  io_rbuffer_count,    &
                  logical_to_bool(io_rbuffer_update==1),   &
                  io_rbuffer_prefetch, &
                  logical_to_bool(io_timing==1))

#if defined(C95_2A)
  ! Initialise IO Timer
  IF (io_timing==io_timing_on) THEN
    CALL io_timing_init()
  END IF
#endif

  DO UNIT=minUnit,maxUnit
    CALL resetFileAttributes(UNIT)
  END DO
  IF (printstatus>=prstatus_normal)                                        &
      CALL umPrint('IO: Initialised IO',src='io')
  initted=.TRUE.

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ioInit

!  Close down io activity, tidy up and maybe print timing data etc.
!
SUBROUTINE ioShutdown()
USE io_configuration_mod, ONLY:                                                &
    io_timing,                                                                 &
    io_timing_on,                                                              &
    release_io_filesystem_profile
USE file_manager, ONLY: release_file_unit
#if defined(C95_2A)
USE io_timing_interfaces_mod, ONLY: io_total_timings
#endif
IMPLICIT NONE
REAL(KIND=jprb):: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='IOSHUTDOWN'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

#if defined(C95_2A)
! Write out PE 0 IO timings if required, for 2B this is handled
! in the portioshutdown call sequence.
IF (( mype == 0)  .OR. L_IO_Server ) THEN
  IF (io_timing==io_timing_on) THEN
!$OMP CRITICAL(internal_write)
    WRITE(umMessage,'(A,I5)') 'IO timings for PE ', mype
!$OMP END CRITICAL(internal_write)
    CALL umPrint(umMessage,src='io')
    CALL io_total_timings()
  END IF
END IF
#endif

CALL release_io_filesystem_profile()

CALL portioShutdown()

CALL release_file_unit("io_reserved_unit", handler="portio")

! Report and flush the Portio shutdown
!$OMP CRITICAL(internal_write)
WRITE(umMessage,'(A)') 'IO: Portio' //                                         &
                       ' has been succesfully finalised and shut down.'
!$OMP END CRITICAL(internal_write)
CALL umPrint(umMessage,src='io')
CALL umPrintFlush()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ioShutdown

!  Internal subroutine to reset the unit state description to defaults
SUBROUTINE resetFileAttributes(UNIT)
IMPLICIT NONE
INTEGER, INTENT(IN) :: UNIT ! Fortran unit
REAL(KIND=jprb)     :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RESETFILEATTRIBUTES'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
fileAttributes(UNIT)%fileName(:)     =' '
fileAttributes(UNIT)%env(:)          =' '
fileAttributes(UNIT)%fileName        =ioNoFileName
fileAttributes(UNIT)%env             =ioNoEnvName
fileAttributes(UNIT)%from_environment=.FALSE.
fileAttributes(UNIT)%policy          =ioNotSet
fileAttributes(UNIT)%state           =ioNotSet
fileAttributes(UNIT)%nextLocation    =ioNoLocation
fileAttributes(UNIT)%fileType        =ioFileTypeUnknown
!$OMP FLUSH(fileAttributes)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE resetFileAttributes


!  Internal subroutine to wrap ereport and provide more io specific
!  system state
SUBROUTINE io_ereport(r,code,m,UNIT)
USE ereport_mod, ONLY: ereport
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN) :: r
CHARACTER(LEN=*),INTENT(IN) :: m
INTEGER, INTENT(IN)         :: code ! code to pass through to ereport
INTEGER, INTENT(IN),                                                       &
    OPTIONAL                :: UNIT ! unit responsible for the error
INTEGER                     :: lcode
INTEGER                     :: i
REAL(KIND=jprb)             :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='IO_EREPORT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
lcode=code ! fix ereport_mod/ereport API mismatch
WRITE(umMessage,'(A,A)')                                                   &
    '****************** IO Error Report ',                                 &
    '***********************************'
CALL umPrint(umMessage,src='io')
IF (PRESENT(UNIT)) THEN
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I5)')'Unit Generating error=',UNIT
!$OMP END CRITICAL(internal_write)
  CALL umPrint(umMessage,src='io')
END IF
IF (printstatus>=prstatus_diag .OR. code>0) THEN
  CALL umPrint('---File States --------------------------')
  DO i=minUnit,maxUnit
    CALL ioStatusPrint(i)
  END DO
  CALL umPrint('---End File States ----------------------')
ELSE
  WRITE(umMessage,'(A,A)')                                                 &
      '*** File states can be reported by ',                               &
      'setting diagnostic output levels **'
  CALL umPrint(umMessage,src='io')
END IF
CALL ereport(r,lcode,m)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE io_ereport


!  Internal subroutine to report the state of a unit to stdout
SUBROUTINE ioStatusPrint(UNIT)

IMPLICIT NONE
INTEGER, INTENT(IN) :: UNIT ! unit to report about
INTEGER             :: state
REAL(KIND=jprb)     :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='IOSTATUSPRINT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

state=fileAttributes(UNIT)%state
IF (testAttribute(state,ioOpen)) THEN
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,i3,A,A)')'Unit ',UNIT,' open on filename ',          &
      TRIM(fileAttributes(UNIT)%filename)
!$OMP END CRITICAL(internal_write)
  CALL umPrint(umMessage,src='io')

  IF (fileAttributes(UNIT)%from_environment) THEN
    WRITE(umMessage,'(A,A)')'  --> Opened from environment variable:',     &
        TRIM(fileAttributes(UNIT)%env)
    CALL umPrint(umMessage,src='io')
  END IF

!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I3,A,L1,A,L1)')                                      &
      '  --> File Type: ',fileAttributes(UNIT)%fileType,                   &
      ' , Read Only: ',testAttribute(state,ioReadOnly),                    &
      ' , Write Only: ',testAttribute(state,ioWriteOnly)
!$OMP END CRITICAL(internal_write)
  CALL umPrint(umMessage,src='io')
  WRITE(umMessage,'(A,L1,A,L1,A,L1,A,L1)')                                         &
      '  --> Local: ',testAttribute(state,ioLocal),                        &
      ' AllLocal: ',testAttribute(state,ioAllLocal),                       &
      ' Remote: ',testAttribute(state,ioRemote),                           &
      ' Broadcast: ',testAttribute(state,ioReadReplicate)
  CALL umPrint(umMessage,src='io')
  CALL umPrint(umMessage,src='io')
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ioStatusPrint


!  Internal Subroutine to check if the calling process has a valid
!  unit number.
SUBROUTINE checkUnit(UNIT)

IMPLICIT NONE
INTEGER, INTENT(IN)   :: UNIT
INTEGER               :: abortCode=10

IF ( UNIT > maxUnit .OR. UNIT < minUnit ) THEN
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I10,A,i3,A,i3)')                                     &
      'Fortran unit out of allowed range, unit = ',UNIT,                   &
      ' RANGE = ',minUnit,' - ',maxUnit
!$OMP END CRITICAL(internal_write)
  CALL umPrint(umMessage,src='io')
  CALL io_ereport ('io:checkUnit',abortCode,                               &
      'Fortran Unit provided is out of range',UNIT=unit)
END IF

END SUBROUTINE checkUnit


!  Internal function to check if the calling process
!  will participate in the operation. This enables us
!  to be able to check parameters sanely
LOGICAL FUNCTION performOperation(UNIT,reading,writing)
IMPLICIT NONE
INTEGER, INTENT(IN) :: UNIT    !Fortran unit
LOGICAL, INTENT(IN),                                                       &
    OPTIONAL        :: reading !Are we reading from the file
LOGICAL, INTENT(IN),                                                       &
    OPTIONAL        :: writing !Are we writing to the file
LOGICAL             :: Lreading
LOGICAL             :: Lwriting
INTEGER             :: abortCode=11

Lreading=.FALSE.
Lwriting=.FALSE.
IF (PRESENT(reading))Lreading=reading
IF (PRESENT(writing))Lwriting=writing

CALL CheckUnit(UNIT)

performOperation=.FALSE.
IF (testAttribute(fileAttributes(UNIT)%state,ioAllLocal)) THEN
  performOperation=.TRUE.
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioLocal)                 &
    .AND. mype==0) THEN
  performOperation=.TRUE.
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioRemote)                &
    .AND. mype==0) THEN
  performOperation=.TRUE.
END IF

! If we are reading we may participate in a broadcast op.
IF (Lreading .AND. broadcast_read(UNIT)) THEN
  performOperation=.TRUE.
END IF

!If we are going to be involved in an IO in any way,
!do some saneness checks
IF (performOperation) THEN
  !the unit should be recognisably open
  IF (.NOT. testAttribute(fileAttributes(UNIT)%state,ioOpen))               &
      CALL io_ereport('io:performOperation',                               &
      abortCode,'Operation attempted to a closed unit',UNIT=unit)
  !We shouldn't  write to a 'read only' defined unit
  IF (Lwriting .AND.                                                        &
      testAttribute(fileAttributes(UNIT)%state,ioReadOnly))                &
      CALL io_ereport('io:performOperation',                               &
      abortCode,'Operation to write to a ReadOnly unit',UNIT=unit)
  !We cant read from a 'write only' defined unit
  IF (Lreading .AND.                                                        &
      testAttribute(fileAttributes(UNIT)%state,ioWriteOnly))               &
      CALL io_ereport('io:performOperation',                               &
      abortCode,'Operation to read from a WriteOnly unit',UNIT=unit)
END IF
END FUNCTION performOperation

!  Report is a unit is remote or not
LOGICAL FUNCTION isRemote(UNIT)

IMPLICIT NONE
INTEGER, INTENT(IN) :: UNIT          !File Unit number
INTEGER             :: STATUS
REAL(KIND=jprb)     :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ISREMOTE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL checkUnit(UNIT)

isRemote=testAttribute(fileAttributes(UNIT)%state,ioRemote)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END FUNCTION isRemote

!  Report is a unit is open or not
LOGICAL FUNCTION is_unit_open(UNIT)

IMPLICIT NONE
INTEGER, INTENT(IN) :: UNIT          !File Unit number
INTEGER             :: STATUS
REAL(KIND=jprb)     :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='IS_UNIT_OPEN'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL ioInit()
CALL checkUnit(UNIT)

is_unit_open=.FALSE.
IF (testAttribute(fileAttributes(UNIT)%state,ioOpen)) THEN
  is_unit_open=.TRUE.
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END FUNCTION is_unit_open


!  Open a file on the given unit
SUBROUTINE file_open (                                                       &
    UNIT,env,env_len,read_write,name_in_environ,error,                       &
    ioLocality,ioPolicy,allowRemap,fileType)
USE IOS_Common, ONLY: IOS_local_ro_files
USE io_configuration_mod, ONLY:                                              &
     io_filesystem_profile,                                                  &
     lustre_filesystem_profile
IMPLICIT NONE
INTEGER, INTENT(IN)               :: UNIT
                                   !Unit number for I/O
CHARACTER(LEN=*),INTENT(IN)       :: env
                                   !Environment name or
                                   !explicit file name
INTEGER, INTENT(IN), OPTIONAL     :: env_len
                                   !Length of env
INTEGER, INTENT(IN), OPTIONAL     :: read_write
                                   ! == 0 read only access,
                                   ! /= 0 read+write
INTEGER, INTENT(IN), OPTIONAL     :: name_in_environ
                                   ! == 0 file name stored in
                                   !environment variable, otherwise
                                   !name specified explicitly
INTEGER, INTENT(IN), OPTIONAL     :: ioPolicy
                                   !Degree of inconsistency allowed
                                   !between file ops before
                                   !generating errors
INTEGER, INTENT(IN), OPTIONAL     :: ioLocality
                                   !Where IO on this unit should
                                   !take place (e.g. all processors,
                                   ! remote services)
LOGICAL, INTENT(IN), OPTIONAL     :: allowRemap
                                   !Specifies that this file has a new
                                   !filename and hence can be remapped
                                   !to a different IO server.
INTEGER, INTENT(OUT), OPTIONAL    :: error
                                   ! == 0 file OPENED
                                   ! /= 0 file NOT OPENED due to error
INTEGER, INTENT(IN), OPTIONAL     :: fileType

! Local copies of optional args
INTEGER                           :: Local_fileType
INTEGER                           :: Local_read_write
INTEGER                           :: Local_nameinenv
INTEGER                           :: Local_envlen
INTEGER                           :: Local_Error

! Need to pass equivalent of Local_read_write to MPI_File_open
! while retaining the read/write only status in the file 
! attributes
INTEGER                           :: mpiio_mode

INTEGER                           :: length          ! length of returned string
INTEGER                           :: info
INTEGER                           :: abortCode=12
INTEGER                           :: warnCode=-12
INTEGER                           :: fileNameLength
INTEGER                           :: newPolicy
INTEGER                           :: newState
INTEGER                           :: mymin
INTEGER                           :: mymax
INTEGER                           :: mymodel_com
INTEGER                           :: server
CHARACTER (LEN=maxFileNameLength) :: fileName
CHARACTER (LEN=30)                :: serveString
REAL(KIND=jprb)                   :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FILE_OPEN'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL ioInit()
CALL checkUnit(UNIT)

Local_fileType=ioFileTypeUnknown
IF (PRESENT(fileType)) THEN
  Local_fileType=fileType
END IF


Local_nameinenv=ioNameProvided
IF (PRESENT(name_in_environ)) THEN
  IF (name_in_environ==ioNameinEnv)                                        &
      Local_nameinenv=ioNameInEnv
END IF

Local_read_write = ioOpenReadWrite
IF (PRESENT(read_write)) THEN
   IF ( IOS_enable_mpiio .AND. Local_filetype == ioFileTypeMPIIO ) THEN
      SELECT CASE( read_write )
      CASE (ioOpenWriteOnly)
         mpiio_mode = IOR(MPL_MODE_WRONLY,MPL_MODE_CREATE)
      CASE (ioOpenReadOnly)
         mpiio_mode = MPL_MODE_RDONLY
      CASE DEFAULT
         mpiio_mode = IOR(MPL_MODE_RDWR,MPL_MODE_CREATE)
      END SELECT
   END IF
   IF (read_write/= ioOpenReadOnly .AND. read_write /= ioOpenWriteOnly) THEN
      Local_read_write = ioOpenReadWrite
   ELSE
      Local_read_write = read_write
   END IF
END IF

IF (PRESENT(env_len)) THEN
  local_envlen=env_len
ELSE
  local_envlen=LEN_TRIM(env)
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First, we perform sanity checks
! - Set params for optional arguments
! - adapt arguments to the IO Server configuration
! - trap inconsistent and illegal requests
! - enforce collective semantics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Check arguments and set default policy
IF (.NOT. PRESENT(ioPolicy)) THEN
  newPolicy=ioPolicyDefault
ELSE
  IF (ioPolicy==ioPolicyStrict .OR.                                          &
      ioPolicy==ioPolicySloppy .OR.                                         &
      ioPolicy==ioPolicyLax) THEN
    newPolicy=ioPolicy
  ELSE
!$OMP CRITICAL(internal_write)
    WRITE(umMessage,'(A,i5,A)')'Policy value of ',ioPolicy,' is invalid'
!$OMP END CRITICAL(internal_write)
    CALL umPrint(umMessage,src='io')
    CALL io_ereport('io:file_open',abortCode,                              &
        'Invalid Optional Policy Argument',UNIT=unit)
  END IF
END IF

! Check arguments and set default Locality
IF (.NOT. PRESENT(ioLocality)) THEN
  newState=ioDefault
ELSE
  newState=ioLocality
END IF

! If there is no IO server and remote is set,
! then set the local flag instead
IF ( testAttribute(newState,ioRemote) .AND. .NOT. l_ios_active()) THEN
  CALL setAttribute(newState,ioLocal)
  CALL clearAttribute(newState,ioRemote)
  CALL umPrint('IO: Switching file mode to local because '//               &
      'there is no IO server',src='io',level=prstatus_diag)
END IF

! Check that user doesn't ask for read broadcast and per task IO
IF (testAttribute(newState,ioAllLocal) .AND.                               &
    testAttribute(newState,ioReadReplicate)) THEN
  CALL io_ereport('io:file_open',abortCode,                                &
      'Requested per task IO and broadcast :-o',UNIT=unit)
END IF

!Route read only files locally if requested.
IF ( IOS_local_ro_files ) THEN
  IF (local_read_write==0) THEN
    IF ( testAttribute(newState,ioRemote)) THEN
      CALL setAttribute(newState,ioLocal)
      CALL clearAttribute(newState,ioRemote)
      IF (printstatus>=prstatus_diag)                                      &
          CALL umPrint('IO: Routing read only file to local',src='io')
    END IF
  END IF
END IF

! Clear the broadcast bit if there is only one processor, so that we
! avoid calling gcom in serial applications
IF (testAttribute(newState,ioReadReplicate) .AND.                            &
    model_procs==1) THEN
  CALL clearAttribute(newState,ioReadReplicate)
  CALL umPrint('IO: Clearing the broadcast read flag '//                   &
      'because there is only one CPU',src='io',level=prstatus_diag)
END IF

IF (printstatus>=prstatus_diag) THEN
  IF (testAttribute(newState,ioReadReplicate)) THEN
!$OMP CRITICAL(internal_write)
    WRITE(umMessage,'(A,I3,A)')'IO: Opening unit ',UNIT,                   &
        ' with collective(broadcast) semantics'
!$OMP END CRITICAL(internal_write)
    CALL umPrint(umMessage,src='io')
  ELSE
!$OMP CRITICAL(internal_write)
    WRITE(umMessage,'(A,I3,A)')'IO: Opening unit ',UNIT,                   &
        ' without collective(broadcast) semantics'
!$OMP END CRITICAL(internal_write)
    CALL umPrint(umMessage,src='io')
  END IF
  IF (local_read_write==0)                                                 &
      CALL umPrint('IO: Read Only mode',src='io')
END IF

Local_Error=0

! Check if the unit is already attached
IF (is_unit_open(UNIT)) THEN
  IF (newPolicy==ioPolicyStrict .OR.                                       &
      fileAttributes(UNIT)%policy==ioPolicyStrict) THEN
    CALL io_ereport('io:file_open',abortCode,                              &
        'Unit already open',UNIT)
  ELSE
    IF (printstatus>=prstatus_diag) THEN
      CALL io_ereport('io:file_open',warnCode,                             &
          'Unit already open - attempting close',UNIT)
    END IF
    CALL file_close(UNIT,env,env_len,name_in_environ,ioNoDelete,Local_Error)
  END IF
END IF


! For small jobs (where its cheap) check for consistency
! We also avoid this for small execs as currently the IBM
! 32 bit builds link with a gcom 64 bit library, so we will skip it
! for the nproc=1 case too.

IF (.NOT. testAttribute(newstate,ioAllLocal) .AND.                           &
    nproc<=256 .AND. nproc>1) THEN
  IF (printstatus>=prstatus_diag) THEN
    CALL umPrint('IO: Checking consistency of unit open request...')
    CALL umPrintFlush()
  END IF
  CALL gc_get_communicator(mymodel_com,Local_Error)

  mymin=UNIT
  mymax=UNIT
  CALL GC_imin(1,nproc,Local_Error,mymin)
  CALL GC_imax(1,nproc,Local_Error,mymax)
  IF (mymin/=mymax) CALL io_ereport('file_open',abortCode,                 &
      'unit no. mismatch',UNIT)

  mymin=local_read_write
  mymax=local_read_write
  CALL GC_imin(1,nproc,Local_Error,mymin)
  CALL GC_imax(1,nproc,Local_Error,mymax)
  IF (mymin/=mymax) CALL io_ereport('file_open',abortCode,                 &
      'read_write flag mismatch',UNIT)

  mymin=Local_nameinenv
  mymax=Local_nameinenv
  CALL GC_imin(1,nproc,Local_Error,mymin)
  CALL GC_imax(1,nproc,Local_Error,mymax)
  IF (mymin/=mymax) CALL io_ereport('file_open',abortCode,                 &
      'name_in_environ mismatch',UNIT)

  mymin=newstate
  mymax=newstate
  CALL GC_imin(1,nproc,Local_Error,mymin)
  CALL GC_imax(1,nproc,Local_Error,mymax)
  IF (mymin/=mymax) CALL io_ereport('file_open',abortCode,                 &
      'locale mismatch',UNIT)

  IF (printstatus>=prstatus_diag) THEN
    CALL umPrint('IO: Valid request')
    CALL umPrintFlush()
  END IF

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Second, we recover the filename from the environment if needed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF (Local_nameinenv == ioNameInEnv) THEN
! Call get_environment_variable directly so we can pass errors to io_ereport.
  CALL GET_ENVIRONMENT_VARIABLE(env,filename,length,info)
  fileNameLength=LEN_TRIM(filename)
  IF (info /= 0 .OR. length == 0 ) THEN
    WRITE(io_message,'(A,A)')'Problem reading environment variable ', TRIM(env)
    CALL io_ereport('io:file_open',abortCode,io_message,                   &
        UNIT=unit)
  END IF
ELSE
  IF (local_envlen>maxFileNameLength) THEN
!$OMP CRITICAL(internal_write)
    WRITE(io_message,'(A,I10,A,I3)')'File name too long (',                &
        local_envlen,') must be less than ',maxFileNameLength
!$OMP END CRITICAL(internal_write)
    CALL io_ereport('io:file_open',abortCode,                              &
        io_message,UNIT=unit)
  END IF
  fileNameLength=local_envlen
  filename  =env(1:fileNameLength)
END IF
!Blank the end of the filename string so we can use trim()
! with impunity
filename(fileNameLength+1:)=' '

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Third, we actually open the file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Local_Error=0
serveString=""
IF (testAttribute(newState,ioAllLocal)) THEN
   IF ( IOS_enable_mpiio .AND. Local_filetype == ioFileTypeMPIIO ) THEN
      IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
      CALL MPL_Info_create(mpiio_info(unit),Local_error)
      IF ( io_filesystem_profile == lustre_filesystem_profile )            &
           CALL apply_file_striping_mpiio( TRIM(filename), mpiio_info(unit) )
      CALL MPL_File_open(model_comm, TRIM(filename), mpiio_mode,           &
           mpiio_info(unit), mpiio_fh(unit), Local_error)
      IF (Local_Error/=0) THEN
         IF ( fileAttributes(UNIT)%policy==ioPolicyStrict) THEN
            CALL io_ereport('io:file_open',abortCode,                      &
                 'An error occurred opening a file',UNIT=unit)
         ELSE
            CALL io_ereport('io:file_open',warnCode,                       &
                 'An error occurred opening a file',UNIT=unit)
         END IF
      END IF
      CALL MPL_File_set_errhandler(mpiio_fh(unit),MPL_ERRORS_ARE_FATAL,    &
           Local_error)
      IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()
   ELSE
      CALL open_single(UNIT,                                               &
           TRIM(filename),LEN_TRIM(filename),local_read_write,1,Local_Error)
   END IF
ELSE IF (testAttribute(newState,ioLocal) .AND. mype==0) THEN
  CALL open_single(UNIT,                                                   &
      TRIM(filename),LEN_TRIM(filename),local_read_write,1,Local_Error)
ELSE IF (testAttribute(newState,ioRemote)) THEN
  server=IOS_Open(UNIT,                                                    &
      TRIM(filename),local_read_write,Local_filetype,allowRemap)
  IF (server /= -1) THEN
!$OMP CRITICAL(internal_write)
    WRITE(serveString,'(A,I4,A)')' (Server = ',server,')'
!$OMP END CRITICAL(internal_write)
  ELSE
!$OMP CRITICAL(internal_write)
    WRITE(serveString,'(A,I4,A)')' (ServerUnknown)'
!$OMP END CRITICAL(internal_write)
  END IF
ELSE IF (mype==0) THEN
  CALL umPrint('IO: FILE_OPEN: Bad Request.')
  Local_Error=-999
ELSE
  Local_Error=0
END IF

IF (Local_Error/=0) THEN
  IF ( fileAttributes(UNIT)%policy==ioPolicyStrict) THEN
    CALL io_ereport('io:file_open',abortCode,                              &
        'An error occurred opening a file',UNIT=unit)
  ELSE
    CALL io_ereport('io:file_open',warnCode,                               &
        'An error occurred opening a file',UNIT=unit)
  END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fourth, we set the state of the file according to what
! actually happened, or die if we failed.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF (Local_Error==0) THEN
  ! File was opened so set any other attributes
  ! even if its not open local to this rank
  IF (Local_nameinenv == ioNameInEnv) THEN
    fileAttributes(UNIT)%env=env
    fileAttributes(UNIT)%env(local_envlen+1:)=' '
    fileAttributes(UNIT)%filename=filename
    fileAttributes(UNIT)%filename(fileNameLength+1:)=' '
    fileAttributes(UNIT)%from_environment=.TRUE.
  ELSE
    fileAttributes(UNIT)%fileName=env
    fileAttributes(UNIT)%fileName(local_envlen+1:)=' '
    fileAttributes(UNIT)%env=ioNoEnvName
    fileAttributes(UNIT)%from_environment=.FALSE.
  END IF
  fileAttributes(UNIT)%state=newState
  fileAttributes(UNIT)%policy=newPolicy
  fileAttributes(UNIT)%fileType=Local_fileType
  CALL setAttribute(fileAttributes(UNIT)%state,ioOpen)
  IF (local_read_write==ioOpenReadOnly)                                    &
      CALL setAttribute(fileAttributes(UNIT)%state,ioReadOnly)
  IF (local_read_write==ioOpenWriteOnly)                                   &
      CALL setAttribute(fileAttributes(UNIT)%state,ioWriteOnly)


!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,A,A,I3,A)')                                          &
      'IO: Open: ',TRIM(fileAttributes(UNIT)%fileName),                    &
      ' on unit ',UNIT,serveString
!$OMP END CRITICAL(internal_write)
  CALL umPrint(umMessage,src='io')
  IF (Local_nameinenv == ioNameInEnv) THEN
    WRITE(umMessage,'(A,A)')'IO: from environment variable ',              &
        TRIM(fileAttributes(UNIT)%env)
    CALL umPrint(umMessage,src='io')
  END IF
!$OMP FLUSH(fileAttributes)
ELSE
  ! We failed.
  IF (Local_nameinenv == ioNameInEnv) THEN
    WRITE(umMessage,'(A,A,A)')                                             &
        'IO: Failure to open file specified by environment ',              &
        'variable: ',                                                      &
        TRIM(env)
    CALL umPrint(umMessage,src='io')
  END IF
  WRITE(umMessage,'(A,A)') 'IO: Failure to open file: ',TRIM(filename)
  CALL umPrint(umMessage,src='io')
  ! Reset the file attributes
  CALL io_ereport('io:file_open',Local_Error,'Failed to open file ' //     &
      TRIM(filename), UNIT=unit)
  CALL resetFileAttributes(UNIT)
END IF

IF ( PRESENT(error) ) error = Local_Error

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE file_open


!  Close the file on the given unit
SUBROUTINE file_close(UNIT,env,env_len,name_in_environ,delete,error)

IMPLICIT NONE
INTEGER, INTENT(IN):: UNIT           !Unit number for I/O
CHARACTER(LEN=*), INTENT(IN)                                               &
              :: env                 !Environment name or filename
INTEGER, INTENT(IN),                                                       &
    OPTIONAL  :: env_len             !Length of env
INTEGER, INTENT(IN), &               ! == 0, do not delete file
    OPTIONAL       :: delete         ! /= 0, delete file
INTEGER, INTENT(IN), &               ! == 0 file name stored in env var
    OPTIONAL  :: name_in_environ     ! /= 0 file name given
INTEGER, INTENT(OUT),                                                      &
    OPTIONAL  :: error               !return code (0=no error)

! Local copies of optional args
INTEGER            :: Local_Error
INTEGER            :: Local_nameinenv
INTEGER            :: Local_delete
INTEGER            :: Local_envlen

INTEGER            :: length         ! length of returned string
INTEGER            :: abortCode=33   !for ereport
INTEGER            :: warnCode=-34   !for ereport
REAL(KIND=jprb)    :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FILE_CLOSE'
INTEGER                           :: fileNameLength
CHARACTER (LEN=maxFileNameLength) :: fileName
INTEGER                           :: info

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

Local_nameinenv=ioNameProvided
IF (PRESENT(name_in_environ)) THEN
  IF (name_in_environ==ioNameInEnv)                                        &
      Local_nameinenv=ioNameInEnv
END IF

Local_delete=ioNoDelete
IF (PRESENT(delete)) THEN
  IF (delete/=ioNoDelete)                                                  &
      Local_delete=ioDelete
END IF

IF (PRESENT(env_len)) THEN
  local_envlen=env_len
ELSE
  local_envlen=LEN_TRIM(env)
END IF



Local_Error = 0
CALL checkUnit(UNIT)

IF (Local_nameinenv == ioNameInEnv) THEN
! Call get_environment_variable directly so we can pass errors to io_ereport.
  CALL GET_ENVIRONMENT_VARIABLE(env,filename,length,info)
  fileNameLength=LEN_TRIM(filename)
  IF (info /= 0 .OR. length == 0 ) THEN
    WRITE(io_message,'(A,A)')'Problem reading environment variable ', TRIM(env)
    CALL io_ereport('io:file_close',abortCode,io_message,                  &
        UNIT=unit)
  END IF
ELSE
  IF (local_envlen>maxFileNameLength) THEN
!$OMP CRITICAL(internal_write)
    WRITE(io_message,'(A,I10,A,I3)')'File name too long (',                &
        local_envlen,') must be less than ',maxFileNameLength
!$OMP END CRITICAL(internal_write)
    CALL io_ereport('io:file_close',abortCode,                             &
        io_message,UNIT=unit)
  END IF
  fileNameLength=local_envlen
  filename      =env(1:fileNameLength)
END IF

!Blank the end of the filename string so we can use trim()
! with impunity
filename(fileNameLength+1:)=' '


IF (is_unit_open(UNIT)) THEN
  IF (testAttribute(fileAttributes(UNIT)%state,ioAllLocal)) THEN
     IF ( IOS_enable_mpiio .AND.                                           &
          testAttribute(fileAttributes(UNIT)%fileType,ioFileTypeMPIIO) ) THEN
        IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
        CALL MPL_File_close( mpiio_fh(unit),Local_error )
        CALL MPL_Info_free( mpiio_info(unit),Local_error ) 
        IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()
     ELSE
        CALL close_single(UNIT,TRIM(filename),LEN_TRIM(filename),1,        &
             Local_delete,Local_Error)
     END IF
  ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioLocal)               &
      .AND. mype==0) THEN
    CALL close_single(UNIT,TRIM(filename),LEN_TRIM(filename),1,            &
        Local_delete,Local_Error)
  ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioRemote)              &
      .AND. mype==0) THEN
     CALL IOS_Close(UNIT,TRIM(filename),Local_delete,                      &
          fileType=fileAttributes(UNIT)%fileType)
  END IF

  IF (printstatus>=prstatus_oper) THEN
!$OMP CRITICAL(internal_write)
    WRITE(umMessage,'(A,A,A,I3)')                                          &
        'IO: Close: ',TRIM(fileAttributes(UNIT)%fileName),                 &
        ' on unit ',UNIT
!$OMP END CRITICAL(internal_write)
    CALL umPrint(umMessage,src='io')
    IF (fileAttributes(UNIT)%from_environment) THEN
      WRITE(umMessage,'(A,A)')'IO: Originally from environment variable ', &
          TRIM(fileAttributes(UNIT)%env)
      CALL umPrint(umMessage,src='io')
    END IF
  END IF

  IF (TRIM(fileName) /= TRIM(fileAttributes(UNIT)%fileName)) THEN
    IF (printstatus>=prstatus_oper .AND. Local_delete /= 0) THEN
      CALL io_ereport('io:file_close ',warncode,                           &
          'Filename to delete does not match with open',                   &
          UNIT=unit)
    END IF
  END IF
ELSE
  IF (printstatus>=prstatus_diag) THEN
    CALL io_ereport('io:file_close',warncode,                              &
        'Tried to close a unit which was not open',UNIT=unit)
  END IF
END IF ! IS Unit Open?



CALL resetFileAttributes(UNIT)

IF ( PRESENT(error) ) error = Local_Error

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE file_close


!  Parallel UM version of buffout
SUBROUTINE buffout64_r                                                       &
    (UNIT,array,numWords,wordsWritten,errorCode)

IMPLICIT NONE
INTEGER, INTENT(IN) :: UNIT
INTEGER, INTENT(OUT),                                                      &
    OPTIONAL        :: wordsWritten    ! Words written
INTEGER, INTENT(IN),                                                       &
    OPTIONAL        :: numWords        ! No. of words to write out
REAL, INTENT(OUT),                                                         &
    OPTIONAL        :: errorCode
REAL(KIND=real64), TARGET,                                                     &
    INTENT(IN)      :: array(:)        ! Array to write out

! Local copies of opt. args
INTEGER             :: Local_WordsWritten    ! Words written
INTEGER             :: Local_NumWords
REAL                :: Local_Error      ! Exit Code
INTEGER             :: mpl_write_err
INTEGER             :: mpl_count_err
INTEGER             :: byteswap_err
INTEGER             :: mpl_write_stat(MPL_STATUS_SIZE)

INTEGER             :: info
INTEGER             :: abortCode=13
REAL(KIND=jprb)     :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BUFFOUT64_R'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (PRESENT(numWords)) THEN
  Local_NumWords=numWords
ELSE
  Local_NumWords=SIZE(array)
END IF
Local_Error=-1.0
Local_WordsWritten=Local_NumWords

IF (performOperation(UNIT,writing=.TRUE.) .AND.                             &
    Local_NumWords>SIZE(array)) THEN
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I0,A,I0)')'Request for ',Local_NumWords,           &
      ' words, is not supported by a buffer of size ',SIZE(array)
!$OMP END CRITICAL(internal_write)
  CALL umPrint(umMessage,src='io')
  CALL io_ereport('io:buffout',abortCode,                                  &
      'Supplied buffer too small',UNIT=unit)
END IF

IF (testAttribute(fileAttributes(UNIT)%state,ioAllLocal)) THEN
   IF (testAttribute(fileAttributes(UNIT)%filetype,ioFileTypeMPIIO) .AND. &
        IOS_Enable_mpiio ) THEN
      IF (get_machine_endianism()==littleEndian)                           &
           byteswap_err = pio_byteswap(C_LOC(array),Local_NumWords,        &
           IOS_BytesPerWord64)
      IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
      CALL MPL_File_write( mpiio_fh(UNIT),array,Local_numwords,MPL_REAL8,  &
           mpl_write_stat,mpl_write_err )
! Errors will be handled by the MPI error handler, but set local_error anyway
      IF ( mpl_write_err/=MPL_SUCCESS ) Local_Error=REAL(mpl_write_err)
      CALL MPL_Get_count( mpl_write_stat,MPL_REAL8,Local_WordsWritten,     &
           mpl_count_err )
      IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()
      IF (get_machine_endianism()==littleEndian)                           &
           byteswap_err = pio_byteswap(C_LOC(array),Local_NumWords,        &
           IOS_BytesPerWord64)
   ELSE
      CALL buffout64_single                                                &
      (UNIT,array,Local_NumWords,Local_WordsWritten,Local_Error)
   END IF
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioLocal)                 &
    .AND. mype==0) THEN
  CALL buffout64_single                                                    &
      (UNIT,array,Local_NumWords,Local_WordsWritten,Local_Error)
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioRemote)                &
    .AND. mype==0) THEN
  CALL ios_write64     (UNIT,fileAttributes(UNIT)%nextLocation,            &
      array(1:Local_NumWords),Local_NumWords)
  fileAttributes(UNIT)%nextLocation=ioNoLocation
END IF

IF ( PRESENT(wordsWritten) ) wordsWritten = Local_WordsWritten
IF ( PRESENT(errorCode) ) errorCode = Local_Error

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE buffout64_r


!  Output 64 bit data to the given unit
SUBROUTINE buffout64_i                                                       &
    (UNIT,array,numWords,wordsWritten,errorCode)

IMPLICIT NONE
INTEGER, INTENT(IN) :: UNIT
INTEGER, INTENT(OUT),                                                      &
    OPTIONAL        :: wordsWritten    ! Words written
INTEGER, INTENT(IN),                                                       &
    OPTIONAL        :: numWords        ! no. of words to write out
REAL, INTENT(OUT),                                                         &
    OPTIONAL        :: errorCode
INTEGER(KIND=integer64), TARGET,                                               &
    INTENT(IN)      :: array(:)

! Local copies of opt. args
INTEGER             :: Local_NumWords
REAL                :: Local_Error       ! Exit Code       ! Exit Code
INTEGER             :: Local_WordsWritten    ! Words written
INTEGER             :: mpl_write_err
INTEGER             :: mpl_count_err
INTEGER             :: byteswap_err
INTEGER             :: mpl_write_stat(MPL_STATUS_SIZE)

INTEGER             :: abortCode=14
INTEGER             :: info
REAL(KIND=jprb)     :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BUFFOUT64_I'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (PRESENT(numWords)) THEN
  Local_NumWords=numWords
ELSE
  Local_NumWords=SIZE(array)
END IF
Local_Error=-1.0
Local_WordsWritten=Local_NumWords

IF (performOperation(UNIT,writing=.TRUE.) .AND.                             &
    Local_NumWords>SIZE(array)) THEN
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I0,A,I0)')'Request for ',Local_NumWords,           &
      ' words, is not supported by a buffer of size ',SIZE(array)
!$OMP END CRITICAL(internal_write)
  CALL umPrint(umMessage,src='io')
  CALL io_ereport('io:buffout64',abortCode,                                &
      'Supplied buffer too small',UNIT=unit)
END IF

IF (testAttribute(fileAttributes(UNIT)%state,ioAllLocal)) THEN
   IF (testAttribute(fileAttributes(UNIT)%filetype,ioFileTypeMPIIO) .AND. &
        IOS_Enable_mpiio ) THEN
      IF (get_machine_endianism()==littleEndian)                           &
           byteswap_err = pio_byteswap(C_LOC(array),Local_NumWords,        &
           IOS_BytesPerWord64)
      IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
      CALL MPL_File_write( mpiio_fh(UNIT),array,Local_numwords,MPL_INTEGER8, &
           mpl_write_stat,mpl_write_err )
! Errors will be handled by the MPI error handler, but set local_error anyway
      IF ( mpl_write_err/=MPL_SUCCESS ) Local_Error=REAL(mpl_write_err)
      CALL MPL_Get_count( mpl_write_stat,MPL_INTEGER8,Local_WordsWritten,  &
           mpl_count_err )
      IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()
      IF (get_machine_endianism()==littleEndian)                           &
           byteswap_err = pio_byteswap(C_LOC(array),Local_NumWords,        &
           IOS_BytesPerWord64)
   ELSE
      CALL buffout64_single                                                &
           (UNIT,array,Local_NumWords,Local_WordsWritten,Local_Error)
   END IF
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioLocal) .AND.            &
    mype==0) THEN
  CALL buffout64_single                                                    &
      (UNIT,array,Local_NumWords,Local_WordsWritten,Local_Error)
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioRemote) .AND.           &
    mype==0) THEN
  CALL ios_write64     (UNIT,fileAttributes(UNIT)%nextLocation,            &
      array(1:Local_NumWords),Local_NumWords)
  fileAttributes(UNIT)%nextLocation=ioNoLocation
END IF

IF ( PRESENT(wordsWritten) ) wordsWritten = Local_WordsWritten
IF ( PRESENT(errorCode) ) errorCode = Local_Error

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE buffout64_i


!  Output 64 bit data to the given unit
SUBROUTINE buffout64_l                                                       &
    (UNIT,array,numWords,wordsWritten,errorCode)

IMPLICIT NONE
INTEGER, INTENT(IN) :: UNIT
INTEGER, INTENT(OUT),                                                      &
    OPTIONAL        :: wordsWritten ! Words written
INTEGER, INTENT(IN),                                                       &
    OPTIONAL        :: numWords     ! no. of words to write out
REAL, INTENT(OUT),                                                         &
    OPTIONAL        :: errorCode
LOGICAL(KIND=logical64), TARGET,                                               &
    INTENT(IN)      :: array(:)     ! Array to write out - WARNING: This array
                                    ! must be contiguous for the call to C_LOC
                                    ! to execute correctly

! Local copies of opt. args
INTEGER             :: Local_WordsWritten
INTEGER             :: Local_NumWords
REAL                :: Local_Error
INTEGER             :: mpl_write_err
INTEGER             :: mpl_count_err
INTEGER             :: byteswap_err
INTEGER             :: mpl_write_stat(MPL_STATUS_SIZE)

INTEGER             :: info
INTEGER             :: abortCode=15
REAL(KIND=jprb)     :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BUFFOUT64_L'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (PRESENT(numWords)) THEN
  Local_NumWords=numWords
ELSE
  Local_NumWords=SIZE(array)
END IF
Local_Error=-1.0
Local_WordsWritten=Local_NumWords

IF (performOperation(UNIT,writing=.TRUE.) .AND.                             &
    Local_NumWords>SIZE(array)) THEN
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I0,A,I0)')'Request for ',Local_NumWords,           &
      ' words, is not supported by a buffer of size ',SIZE(array)
!$OMP END CRITICAL(internal_write)
  CALL umPrint(umMessage,src='io')
  CALL io_ereport('io:buffout64',abortCode,                                &
      'Supplied buffer too small',UNIT=unit)
END IF

IF (testAttribute(fileAttributes(UNIT)%state,ioAllLocal)) THEN
   IF (testAttribute(fileAttributes(UNIT)%filetype,ioFileTypeMPIIO) .AND. &
        IOS_Enable_mpiio ) THEN
      IF (get_machine_endianism()==littleEndian)                           &
           byteswap_err = pio_byteswap(C_LOC(array(1)),Local_NumWords,        &
           IOS_BytesPerWord64)
      IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
      CALL MPL_File_write( mpiio_fh(UNIT),array,Local_numwords,MPL_INTEGER8, &
           mpl_write_stat,mpl_write_err )
! Errors will be handled by the MPI error handler, but set local_error anyway
      IF ( mpl_write_err/=MPL_SUCCESS ) Local_Error=REAL(mpl_write_err)
      CALL MPL_Get_count( mpl_write_stat,MPL_INTEGER8,Local_WordsWritten,  &
           mpl_count_err )
      IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()
      IF (get_machine_endianism()==littleEndian)                           &
           byteswap_err = pio_byteswap(C_LOC(array(1)),Local_NumWords,        &
           IOS_BytesPerWord64)
   ELSE
      CALL buffout64_single                                                 &
           (UNIT,array,Local_NumWords,Local_WordsWritten,Local_Error)
   END IF
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioLocal) .AND.            &
    mype==0) THEN
  CALL buffout64_single                                                    &
      (UNIT,array,Local_NumWords,Local_WordsWritten,Local_Error)
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioRemote) .AND.           &
    mype==0) THEN
  CALL ios_write64     (UNIT,fileAttributes(UNIT)%nextLocation,            &
      array(1:Local_NumWords),Local_NumWords)
  fileAttributes(UNIT)%nextLocation=ioNoLocation
END IF

IF ( PRESENT(wordsWritten) ) wordsWritten = Local_WordsWritten
IF ( PRESENT(errorCode) ) errorCode = Local_Error

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE buffout64_l


!  Output 64 bit data to the given unit
SUBROUTINE buffout64_r2D                                                     &
    (UNIT,array,numWords,wordsWritten,errorCode)

IMPLICIT NONE
INTEGER, INTENT(IN) :: UNIT
INTEGER, INTENT(IN),                                                       &
    OPTIONAL        :: numWords     ! no. of words to write out
INTEGER, INTENT(OUT),                                                      &
    OPTIONAL        :: wordsWritten ! Words written
REAL, INTENT(OUT),                                                         &
    OPTIONAL        :: errorCode
REAL(KIND=real64), TARGET,                                                     &
    INTENT(IN)      :: array(:,:)   ! Array to write data from

! Local copies of opt. args
INTEGER             :: Local_NumWords
REAL                :: Local_Error
INTEGER             :: Local_WordsWritten
INTEGER             :: mpl_write_err
INTEGER             :: mpl_count_err
INTEGER             :: byteswap_err
INTEGER             :: mpl_write_stat(MPL_STATUS_SIZE)

INTEGER             :: info
INTEGER             :: abortCode=16
REAL(KIND=jprb)     :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BUFFOUT64_R2D'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (PRESENT(numWords)) THEN
  Local_NumWords=numWords
ELSE
  Local_NumWords=SIZE(array)
END IF
Local_Error=-1.0
Local_WordsWritten=Local_NumWords

IF (performOperation(UNIT,writing=.TRUE.) .AND.                             &
    Local_NumWords>SIZE(array)) THEN
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I0,A,I0)')'Request for ',Local_NumWords,           &
      ' words, is not supported by a buffer of size ',SIZE(array)
!$OMP END CRITICAL(internal_write)
  CALL umPrint(umMessage,src='io')
  CALL io_ereport('io:buffout64',abortCode,                                &
      'Supplied buffer too small',UNIT=unit)
END IF

IF (testAttribute(fileAttributes(UNIT)%state,ioAllLocal)) THEN
   IF (testAttribute(fileAttributes(UNIT)%filetype,ioFileTypeMPIIO) .AND. &
        IOS_Enable_mpiio ) THEN
      IF (get_machine_endianism()==littleEndian)                           &
           byteswap_err = pio_byteswap(C_LOC(array),Local_NumWords,        &
           IOS_BytesPerWord64)
      IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
      CALL MPL_File_write( mpiio_fh(UNIT),array,Local_numwords,MPL_REAL8,  &
           mpl_write_stat,mpl_write_err )
! Errors will be handled by the MPI error handler, but set local_error anyway
      IF ( mpl_write_err/=MPL_SUCCESS ) Local_Error=REAL(mpl_write_err)
      CALL MPL_Get_count( mpl_write_stat,MPL_REAL8,Local_WordsWritten,     &
           mpl_count_err )
      IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()
      IF (get_machine_endianism()==littleEndian)                           &
           byteswap_err = pio_byteswap(C_LOC(array),Local_NumWords,        &
           IOS_BytesPerWord64)
   ELSE
      CALL buffout64_single                                                 &
           (UNIT,array,Local_NumWords,Local_WordsWritten,Local_Error)
   END IF
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioLocal) .AND.            &
    mype==0) THEN
  CALL buffout64_single                                                    &
      (UNIT,array,Local_NumWords,Local_WordsWritten,Local_Error)
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioRemote) .AND.          &
    mype==0) THEN
  CALL ios_write64     (UNIT,fileAttributes(UNIT)%nextLocation,            &
      array(:,:),Local_NumWords)
  fileAttributes(UNIT)%nextLocation=ioNoLocation
END IF

IF ( PRESENT(wordsWritten) ) wordsWritten = Local_WordsWritten
IF ( PRESENT(errorCode) ) errorCode = Local_Error

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE buffout64_r2D


!  Output 64 bit data to the given unit
SUBROUTINE buffout64_i2D                                                     &
    (UNIT,array,numWords,wordsWritten,errorCode)

IMPLICIT NONE
INTEGER, INTENT(IN) :: UNIT
INTEGER(KIND=integer64), TARGET,                                               &
    INTENT(IN)      :: array(:,:)   ! Array to write data from
INTEGER, INTENT(IN),                                                       &
    OPTIONAL        :: numWords     ! no. of words to write out
INTEGER, INTENT(OUT),                                                      &
    OPTIONAL        :: wordsWritten ! Words written
REAL, INTENT(OUT),                                                         &
    OPTIONAL        :: errorCode

! Local copies of opt. args
INTEGER             :: Local_NumWords
REAL                :: Local_Error
INTEGER             :: Local_WordsWritten
INTEGER             :: mpl_write_err
INTEGER             :: mpl_count_err
INTEGER             :: byteswap_err
INTEGER             :: mpl_write_stat(MPL_STATUS_SIZE)

INTEGER             :: info
INTEGER             :: abortCode=17
REAL(KIND=jprb)     :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BUFFOUT64_I2D'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (PRESENT(numWords)) THEN
  Local_NumWords=numWords
ELSE
  Local_NumWords=SIZE(array)
END IF
Local_Error=-1.0
Local_WordsWritten=Local_NumWords

IF (performOperation(UNIT,writing=.TRUE.) .AND.                             &
    Local_NumWords>SIZE(array)) THEN
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I0,A,I0)')'Request for ',Local_NumWords,           &
      ' words, is not supported by a buffer of size ',SIZE(array)
!$OMP END CRITICAL(internal_write)
  CALL umPrint(umMessage,src='io')
  CALL io_ereport('io:buffout64',abortCode,                                &
      'Supplied buffer too small',UNIT=unit)
END IF

IF (testAttribute(fileAttributes(UNIT)%state,ioAllLocal)) THEN
   IF (testAttribute(fileAttributes(UNIT)%filetype,ioFileTypeMPIIO) .AND. &
        IOS_Enable_mpiio ) THEN
      IF (get_machine_endianism()==littleEndian)                           &
           byteswap_err = pio_byteswap(C_LOC(array),Local_NumWords,        &
           IOS_BytesPerWord64)
      IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
      CALL MPL_File_write( mpiio_fh(UNIT),array,Local_numwords,MPL_INTEGER8, &
           mpl_write_stat,mpl_write_err )
! Errors will be handled by the MPI error handler, but set local_error anyway
      IF ( mpl_write_err/=MPL_SUCCESS ) Local_Error=REAL(mpl_write_err)
      CALL MPL_Get_count( mpl_write_stat,MPL_INTEGER8,Local_WordsWritten,  &
           mpl_count_err )
      IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()
      IF (get_machine_endianism()==littleEndian)                           &
           byteswap_err = pio_byteswap(C_LOC(array),Local_NumWords,        &
           IOS_BytesPerWord64)
   ELSE
      CALL buffout64_single                                                 &
           (UNIT,array,Local_NumWords,Local_WordsWritten,Local_Error)
   END IF
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioLocal) .AND.            &
    mype==0) THEN
  CALL buffout64_single                                                    &
      (UNIT,array,Local_NumWords,Local_WordsWritten,Local_Error)
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioRemote) .AND.           &
    mype==0) THEN
  CALL ios_write64     (UNIT,fileAttributes(UNIT)%nextLocation,            &
      array(:,:),Local_NumWords)
  fileAttributes(UNIT)%nextLocation=ioNoLocation
END IF

IF ( PRESENT(wordsWritten) ) wordsWritten = Local_WordsWritten
IF ( PRESENT(errorCode) ) errorCode = Local_Error

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE buffout64_i2D


!  Output 64 bit data to the given unit
SUBROUTINE buffout64_l2D                                                     &
    (UNIT,array,numWords,wordsWritten,errorCode)

IMPLICIT NONE
INTEGER, INTENT(IN) :: UNIT
LOGICAL(KIND=logical64), TARGET,                                               &
    INTENT(IN)      :: array(:,:)   ! Array to write data from - WARNING: This
                                    ! array must be contiguous for the call to
                                    ! C_LOC to execute correctly
INTEGER, INTENT(IN),                                                       &
    OPTIONAL        :: numWords     ! no. of words to write out
INTEGER, INTENT(OUT),                                                      &
    OPTIONAL        :: wordsWritten ! Words written
REAL, INTENT(OUT),                                                         &
    OPTIONAL        :: errorCode

! Local copies of opt. args
INTEGER             :: Local_NumWords
REAL                :: Local_Error
INTEGER             :: Local_WordsWritten
INTEGER             :: mpl_write_err
INTEGER             :: mpl_count_err
INTEGER             :: byteswap_err
INTEGER             :: mpl_write_stat(MPL_STATUS_SIZE)

INTEGER             :: info
INTEGER             :: abortCode=18
REAL(KIND=jprb)     :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BUFFOUT64_L2D'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (PRESENT(numWords)) THEN
  Local_NumWords=numWords
ELSE
  Local_NumWords=SIZE(array)
END IF
Local_Error=-1.0
Local_WordsWritten=Local_NumWords

IF (performOperation(UNIT,writing=.TRUE.) .AND.                             &
    Local_NumWords>SIZE(array)) THEN
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I0,A,I0)')'Request for ',Local_NumWords,           &
      ' words, is not supported by a buffer of size ',SIZE(array)
!$OMP END CRITICAL(internal_write)
  CALL umPrint(umMessage,src='io')
  CALL io_ereport('io:buffout64',abortCode,                                &
      'Supplied buffer too small',UNIT=unit)
END IF

IF (testAttribute(fileAttributes(UNIT)%state,ioAllLocal)) THEN
   IF (testAttribute(fileAttributes(UNIT)%filetype,ioFileTypeMPIIO) .AND. &
        IOS_Enable_mpiio ) THEN
      IF (get_machine_endianism()==littleEndian)                           &
           byteswap_err = pio_byteswap(C_LOC(array(1,1)),Local_NumWords,       &
           IOS_BytesPerWord64)
      IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
      CALL MPL_File_write( mpiio_fh(UNIT),array,Local_numwords,MPL_INTEGER8, &
           mpl_write_stat,mpl_write_err )
! Errors will be handled by the MPI error handler, but set local_error anyway
      IF ( mpl_write_err/=MPL_SUCCESS ) Local_Error=REAL(mpl_write_err)
      CALL MPL_Get_count( mpl_write_stat,MPL_INTEGER8,Local_WordsWritten,  &
           mpl_count_err )
      IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()
      IF (get_machine_endianism()==littleEndian)                           &
           byteswap_err = pio_byteswap(C_LOC(array(1,1)),Local_NumWords,       &
           IOS_BytesPerWord64)
   ELSE
      CALL buffout64_single                                                 &
           (UNIT,array,Local_NumWords,Local_WordsWritten,Local_Error)
   END IF
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioLocal) .AND.            &
    mype==0) THEN
  CALL buffout64_single                                                    &
      (UNIT,array,Local_NumWords,Local_WordsWritten,Local_Error)
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioRemote) .AND.          &
    mype==0) THEN
  CALL ios_write64     (UNIT,fileAttributes(UNIT)%nextLocation,            &
      array(:,:),Local_NumWords)
  fileAttributes(UNIT)%nextLocation=ioNoLocation
END IF

IF ( PRESENT(wordsWritten) ) wordsWritten = Local_WordsWritten
IF ( PRESENT(errorCode) ) errorCode = Local_Error

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE buffout64_l2D


!  Output 32 bit data to the given unit
SUBROUTINE buffout32_r                                                       &
    (UNIT,array,numWords,wordsWritten,errorCode)

IMPLICIT NONE
INTEGER, INTENT(IN)  :: UNIT         ! Fortran unit number
INTEGER, INTENT(IN),                                                       &
    OPTIONAL         :: numWords     ! no. of words to write out
REAL(KIND=real32), TARGET,                                                     &
    INTENT(IN)       :: array(:)     ! Array to write out
REAL, INTENT(OUT),                                                         &
    OPTIONAL         :: errorCode
INTEGER, INTENT(OUT),                                                      &
    OPTIONAL        :: wordsWritten ! Words written

! Local copies of opt. args
REAL                 :: Local_Error
INTEGER              :: Local_WordsWritten
INTEGER              :: Local_NumWords
INTEGER              :: mpl_write_err
INTEGER              :: mpl_count_err
INTEGER              :: byteswap_err
INTEGER              :: mpl_write_stat(MPL_STATUS_SIZE)

INTEGER              :: info         ! Broadcast information code
INTEGER              :: abortCode=19
REAL(KIND=jprb)      :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BUFFOUT32_R'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (PRESENT(numWords)) THEN
  Local_NumWords=numWords
ELSE
  Local_NumWords=SIZE(array)
END IF
Local_Error=-1.0
Local_WordsWritten=Local_NumWords

IF (performOperation(UNIT,writing=.TRUE.) .AND.                             &
    Local_NumWords>2*SIZE(array)) THEN
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I0,A)')'Request for ',Local_NumWords,' 32 bit words'
!$OMP END CRITICAL(internal_write)
  CALL umPrint(umMessage,src='io')
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,A,I0)')'is not supported by a 64 bit ',             &
      'real buffer of size ',                                              &
      SIZE(array)
!$OMP END CRITICAL(internal_write)
  CALL umPrint(umMessage,src='io')
  CALL io_ereport('io:buffout32',abortCode,                                &
      'Supplied buffer too small',UNIT=unit)
END IF

IF (testAttribute(fileAttributes(UNIT)%state,ioAllLocal)) THEN
   IF (testAttribute(fileAttributes(UNIT)%filetype,ioFileTypeMPIIO) .AND. &
        IOS_Enable_mpiio ) THEN
      IF (get_machine_endianism()==littleEndian)                           &
           byteswap_err = pio_byteswap(C_LOC(array),Local_NumWords,        &
           IOS_BytesPerWord32)
      IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
      CALL MPL_File_write( mpiio_fh(UNIT),array,Local_numwords,MPL_REAL4,  &
           mpl_write_stat,mpl_write_err )
! Errors will be handled by the MPI error handler, but set local_error anyway
      IF ( mpl_write_err/=MPL_SUCCESS ) Local_Error=REAL(mpl_write_err)
      CALL MPL_Get_count( mpl_write_stat,MPL_REAL4,Local_WordsWritten,     &
           mpl_count_err )
      IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()
      IF (get_machine_endianism()==littleEndian)                           &
           byteswap_err = pio_byteswap(C_LOC(array),Local_NumWords,        &
           IOS_BytesPerWord32)
   ELSE
      CALL buffout32_single                                                &
           (UNIT,array,Local_NumWords,Local_WordsWritten,Local_Error)
   END IF
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioLocal) .AND.            &
    mype==0) THEN
  CALL buffout32_single                                                    &
      (UNIT,array,Local_NumWords,Local_WordsWritten,Local_Error)
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioRemote) .AND.           &
    mype==0) THEN
  CALL ios_write32   (UNIT,fileAttributes(UNIT)%nextLocation,              &
      array(1:Local_NumWords),Local_NumWords)
  fileAttributes(UNIT)%nextLocation=ioNoLocation
END IF

IF ( PRESENT(wordsWritten) ) wordsWritten = Local_WordsWritten
IF ( PRESENT(errorCode) ) errorCode = Local_Error

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE buffout32_r


!  Output 32 bit data to the given unit
SUBROUTINE buffout32_i                                                       &
    (UNIT,array,numWords,wordsWritten,errorCode)

IMPLICIT NONE
INTEGER, INTENT(IN)  :: UNIT         ! Fortran unit number
INTEGER, INTENT(IN),                                                       &
    OPTIONAL         :: numWords     ! no. of words to write out
INTEGER(KIND=integer32), TARGET,                                               &
    INTENT(IN)       :: array(:)     ! Array to write out
REAL, INTENT(OUT),                                                         &
    OPTIONAL         :: errorCode
INTEGER, INTENT(OUT),                                                      &
    OPTIONAL        :: wordsWritten ! Words written

! Local copies of opt. args
INTEGER              :: Local_WordsWritten
INTEGER              :: Local_NumWords
REAL                 :: Local_Error
INTEGER              :: mpl_write_err
INTEGER              :: mpl_count_err
INTEGER              :: byteswap_err
INTEGER              :: mpl_write_stat(MPL_STATUS_SIZE)

INTEGER              :: info         ! Broadcast information code
INTEGER              :: abortCode=20
REAL(KIND=jprb)      :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BUFFOUT32_I'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (PRESENT(numWords)) THEN
  Local_NumWords=numWords
ELSE
  Local_NumWords=SIZE(array)
END IF
Local_Error=-1.0
Local_WordsWritten=Local_NumWords

IF (performOperation(UNIT,writing=.TRUE.) .AND.                             &
    Local_NumWords>2*SIZE(array)) THEN
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I0,A)')'Request for ',Local_NumWords,' 32 bit words'
!$OMP END CRITICAL(internal_write)
  CALL umPrint(umMessage,src='io')
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,A,I0)')'is not supported by a 64 bit int ',         &
      'buffer of size ',                                                   &
      SIZE(array)
!$OMP END CRITICAL(internal_write)
  CALL umPrint(umMessage,src='io')
  CALL io_ereport('io:buffout32',abortCode,                                &
      'Supplied buffer too small',UNIT=unit)
END IF

IF (testAttribute(fileAttributes(UNIT)%state,ioAllLocal)) THEN
   IF (testAttribute(fileAttributes(UNIT)%filetype,ioFileTypeMPIIO) .AND. &
        IOS_Enable_mpiio ) THEN
      IF (get_machine_endianism()==littleEndian)                           &
           byteswap_err = pio_byteswap(C_LOC(array),Local_NumWords,        &
           IOS_BytesPerWord32)
      IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
      CALL MPL_File_write( mpiio_fh(UNIT),array,Local_numwords,MPL_INTEGER4, &
           mpl_write_stat,mpl_write_err )
! Errors will be handled by the MPI error handler, but set local_error anyway
      IF ( mpl_write_err/=MPL_SUCCESS ) Local_Error=REAL(mpl_write_err)
      CALL MPL_Get_count( mpl_write_stat,MPL_INTEGER4,Local_WordsWritten,  &
           mpl_count_err )
      IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()
      IF (get_machine_endianism()==littleEndian)                           &
           byteswap_err = pio_byteswap(C_LOC(array),Local_NumWords,        &
           IOS_BytesPerWord32)
   ELSE
      CALL buffout32_single                                                &
           (UNIT,array,Local_NumWords,Local_WordsWritten,Local_Error)
   END IF
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioLocal) .AND.            &
    mype==0) THEN
  CALL buffout32_single                                                    &
      (UNIT,array,Local_NumWords,Local_WordsWritten,Local_Error)
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioRemote) .AND.          &
    mype==0) THEN
  CALL ios_write32   (UNIT,fileAttributes(UNIT)%nextLocation,              &
      array(1:Local_NumWords),Local_NumWords)
  fileAttributes(UNIT)%nextLocation=ioNoLocation
END IF

IF ( PRESENT(wordsWritten) ) wordsWritten = Local_WordsWritten
IF ( PRESENT(errorCode) ) errorCode = Local_Error

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE buffout32_i


!  Output 32 bit data to the given unit
SUBROUTINE buffout32_i2D                                                     &
    (UNIT,array,numWords,wordsWritten,errorCode)

IMPLICIT NONE
INTEGER, INTENT(IN)  :: UNIT         ! Fortran unit number
INTEGER, INTENT(IN),                                                       &
    OPTIONAL         :: numWords     ! no. of words to write out
INTEGER(KIND=integer32), TARGET,                                               &
    INTENT(IN)       :: array(:,:)   ! Array to write out
REAL, INTENT(OUT),                                                         &
    OPTIONAL         :: errorCode
INTEGER, INTENT(OUT),                                                      &
    OPTIONAL        :: wordsWritten ! Words written

! Local copies of opt. args
INTEGER              :: Local_WordsWritten
INTEGER              :: Local_NumWords
REAL                 :: Local_Error
INTEGER              :: mpl_write_err
INTEGER              :: mpl_count_err
INTEGER              :: byteswap_err
INTEGER              :: mpl_write_stat(MPL_STATUS_SIZE)

INTEGER              :: info         ! Broadcast information code
INTEGER              :: abortCode=21
REAL(KIND=jprb)      :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BUFFOUT32_I2D'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (PRESENT(numWords)) THEN
  Local_NumWords=numWords
ELSE
  Local_NumWords=SIZE(array)
END IF
Local_Error=-1.0
Local_WordsWritten=Local_NumWords

IF (performOperation(UNIT,writing=.TRUE.) .AND.                             &
    Local_NumWords>2*SIZE(array)) THEN
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I0,A)')'Request for ',Local_NumWords,' 32 bit words'
!$OMP END CRITICAL(internal_write)
  CALL umPrint(umMessage,src='io')
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I0)')                                               &
      'is not supported by a 64 bit 2D int buffer of size ',SIZE(array)
!$OMP END CRITICAL(internal_write)
  CALL umPrint(umMessage,src='io')
  CALL io_ereport('io:buffout32',abortCode,                                &
      'Supplied buffer too small',UNIT=unit)
END IF

IF (testAttribute(fileAttributes(UNIT)%state,ioAllLocal)) THEN
   IF (testAttribute(fileAttributes(UNIT)%filetype,ioFileTypeMPIIO) .AND. &
        IOS_Enable_mpiio ) THEN
      IF (get_machine_endianism()==littleEndian)                           &
           byteswap_err = pio_byteswap(C_LOC(array),Local_NumWords,        &
           IOS_BytesPerWord32)
      IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
      CALL MPL_File_write( mpiio_fh(UNIT),array,Local_numwords,MPL_INTEGER4, &
           mpl_write_stat,mpl_write_err )
! Errors will be handled by the MPI error handler, but set local_error anyway
      IF ( mpl_write_err/=MPL_SUCCESS ) Local_Error=REAL(mpl_write_err)
      CALL MPL_Get_count( mpl_write_stat,MPL_INTEGER4,Local_WordsWritten,  &
           mpl_count_err )
      IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()
      IF (get_machine_endianism()==littleEndian)                           &
           byteswap_err = pio_byteswap(C_LOC(array),Local_NumWords,        &
           IOS_BytesPerWord32)
   ELSE
      CALL buffout32_single                                                 &
           (UNIT,array,Local_NumWords,Local_WordsWritten,Local_Error)
   END IF
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioLocal) .AND.            &
    mype==0) THEN
  CALL buffout32_single                                                    &
      (UNIT,array,Local_NumWords,Local_WordsWritten,Local_Error)
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioRemote) .AND.           &
    mype==0) THEN
  CALL ios_write32   (UNIT,fileAttributes(UNIT)%nextLocation,              &
      array(:,:),Local_NumWords)
  fileAttributes(UNIT)%nextLocation=ioNoLocation
END IF

IF ( PRESENT(wordsWritten) ) wordsWritten = Local_WordsWritten
IF ( PRESENT(errorCode) ) errorCode = Local_Error

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE buffout32_i2D


!  Read 64 bit data from the given unit
SUBROUTINE buffin64_r(UNIT,array,numWords,wordsRead,errorCode)

IMPLICIT NONE
INTEGER, INTENT(IN)  :: UNIT         ! Fortran unit number
INTEGER, INTENT(IN),                                                       &
    OPTIONAL         :: numWords     ! no. of words to read in
REAL(KIND=real64), TARGET,                                                     &
    INTENT(OUT)      :: array(:)     ! Array to read in
REAL, INTENT(OUT),                                                         &
    OPTIONAL         :: errorCode
INTEGER, INTENT(OUT),                                                      &
    OPTIONAL         :: wordsRead    ! No. of words read

! Local copies of opt. args
INTEGER              :: Local_NumWords
INTEGER              :: Local_WordsRead
REAL                 :: Local_Error
INTEGER              :: mpl_read_err
INTEGER              :: mpl_count_err
INTEGER              :: byteswap_err
INTEGER              :: mpl_read_stat(MPL_STATUS_SIZE)

INTEGER              :: abortCode=22
INTEGER              :: info         ! Broadcast information code
REAL(KIND=jprb)      :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BUFFIN64_R'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (PRESENT(numWords)) THEN
  Local_NumWords=numWords
ELSE
  Local_NumWords=SIZE(array)
END IF
Local_Error=-1.0
Local_WordsRead=Local_NumWords

IF (performOperation(UNIT,reading=.TRUE.) .AND.                             &
    Local_NumWords>SIZE(array)) THEN
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I0,A,I0)')'Request for ',Local_NumWords,           &
      ' words, is not supported by a buffer of size ',SIZE(array)
!$OMP END CRITICAL(internal_write)
  CALL umPrint(umMessage,src='io')
  CALL io_ereport('io:buffin',abortCode,                                   &
      'Supplied buffer too small',UNIT=unit)
END IF

IF (testAttribute(fileAttributes(UNIT)%state,ioAllLocal)) THEN
   IF (testAttribute(fileAttributes(UNIT)%filetype,ioFileTypeMPIIO) .AND.  &
        IOS_Enable_mpiio ) THEN
      IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
      CALL MPL_File_read( mpiio_fh(UNIT),array,Local_numwords,MPL_REAL8,   &
           mpl_read_stat,mpl_read_err )
! Errors will be handled by the MPI error handler, but set local_error anyway
      IF ( mpl_read_err/=MPL_SUCCESS ) Local_Error=REAL(mpl_read_err)
      CALL MPL_Get_count( mpl_read_stat,MPL_REAL8,Local_WordsRead,         &
           mpl_count_err )
      IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()
      IF (get_machine_endianism()==littleEndian)                           &
           byteswap_err = pio_byteswap(C_LOC(array),Local_NumWords,        &
           IOS_BytesPerWord64)
   ELSE
      CALL buffin64_single                                                 &
           (UNIT,array,Local_NumWords,Local_WordsRead,Local_Error)
   END IF
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioLocal)) THEN
  IF (mype==0)                                                             &
      CALL buffin64_single                                                 &
      (UNIT,array,Local_NumWords,Local_WordsRead,Local_Error)
  IF (broadcast_read(UNIT))                                                 &
      CALL gc_rbcast(1,Local_NumWords,0,nproc,info,array)
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioRemote)) THEN
  IF (broadcast_read(UNIT)) THEN
    CALL ios_read64     (UNIT,fileAttributes(UNIT)%nextLocation,           &
        array(1:Local_NumWords),Local_NumWords,IOS_Read_Broadcast)
  ELSE
    CALL ios_read64     (UNIT,fileAttributes(UNIT)%nextLocation,           &
        array(1:Local_NumWords),Local_NumWords,IOS_Read_NoBroadcast)
  END IF
  fileAttributes(UNIT)%nextLocation=ioNoLocation
END IF

IF (Local_Error /= -1 .OR. Local_WordsRead /= Local_NumWords ) THEN
!$OMP CRITICAL(internal_write)
  WRITE(io_message,'(A,F5.2,A,I0,A,I0)')                                 &
      'Error in buffin errorCode=',Local_Error,                            &
      ' len=',Local_WordsRead,'/',Local_NumWords
!$OMP END CRITICAL(internal_write)
  CALL io_ereport('io:buffin',abortCode,io_message,                        &
      UNIT=unit)
END IF

IF ( PRESENT(wordsRead) ) wordsRead = Local_WordsRead
IF ( PRESENT(errorCode) ) errorCode = Local_Error

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE buffin64_r


!  Read 64 bit data from the given unit
SUBROUTINE buffin64_i(UNIT,array,numWords,wordsRead,errorCode)

IMPLICIT NONE
INTEGER, INTENT(IN)  :: UNIT         ! Fortran unit number
INTEGER, INTENT(IN),                                                       &
    OPTIONAL         :: numWords     ! no. of words to read in
INTEGER(KIND=integer64), TARGET,                                               &
    INTENT(OUT)      :: array(:)     ! Array to read in
INTEGER, INTENT(OUT),                                                      &
    OPTIONAL         :: wordsRead    ! No. of words read
REAL, INTENT(OUT),                                                         &
    OPTIONAL         :: errorCode

! Local copies of opt. args
INTEGER              :: Local_NumWords
INTEGER              :: Local_WordsRead
REAL                 :: Local_Error
INTEGER              :: mpl_read_err
INTEGER              :: mpl_count_err
INTEGER              :: byteswap_err
INTEGER              :: mpl_read_stat(MPL_STATUS_SIZE)

INTEGER              :: info         ! Broadcast information code
INTEGER              :: abortCode=24
REAL(KIND=jprb)      :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BUFFIN64_I'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (PRESENT(numWords)) THEN
  Local_NumWords=numWords
ELSE
  Local_NumWords=SIZE(array)
END IF
Local_Error=-1.0
Local_WordsRead=Local_NumWords

IF (performOperation(UNIT,reading=.TRUE.) .AND.                             &
    Local_NumWords>SIZE(array)) THEN
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I0,A,I0)')'Request for ',Local_NumWords,           &
      ' words, is not supported by a buffer of size ',SIZE(array)
!$OMP END CRITICAL(internal_write)
  CALL umPrint(umMessage,src='io')
  CALL io_ereport('io:buffin',abortCode,                                   &
      'Supplied buffer too small',UNIT=unit)
END IF

IF (testAttribute(fileAttributes(UNIT)%state,ioAllLocal)) THEN
   IF (testAttribute(fileAttributes(UNIT)%filetype,ioFileTypeMPIIO) .AND.  &
        IOS_Enable_mpiio ) THEN
      IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
      CALL MPL_File_read( mpiio_fh(UNIT),array,Local_numwords,MPL_INTEGER8,&
           mpl_read_stat,mpl_read_err )
! Errors will be handled by the MPI error handler, but set local_error anyway
      IF ( mpl_read_err/=MPL_SUCCESS ) Local_Error=REAL(mpl_read_err)
      CALL MPL_Get_count( mpl_read_stat,MPL_INTEGER8,Local_WordsRead,      &
           mpl_count_err )
      IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()
      IF (get_machine_endianism()==littleEndian)                           &
           byteswap_err = pio_byteswap(C_LOC(array),Local_NumWords,        &
           IOS_BytesPerWord64)
   ELSE
      CALL buffin64_single                                                 &
           (UNIT,array,Local_NumWords,Local_WordsRead,Local_Error)
   END IF
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioLocal)) THEN
  IF (mype==0)                                                             &
      CALL buffin64_single                                                 &
      (UNIT,array,Local_NumWords,Local_WordsRead,Local_Error)
  IF (broadcast_read(UNIT))                                                 &
      CALL gc_ibcast(1,Local_NumWords,0,nproc,info,array)
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioRemote)) THEN
  IF (broadcast_read(UNIT)) THEN
    CALL ios_read64     (UNIT,fileAttributes(UNIT)%nextLocation,           &
        array(1:Local_NumWords),Local_NumWords,IOS_Read_Broadcast)
  ELSE
    CALL ios_read64     (UNIT,fileAttributes(UNIT)%nextLocation,           &
        array(1:Local_NumWords),Local_NumWords,IOS_Read_NoBroadcast)
  END IF
  fileAttributes(UNIT)%nextLocation=ioNoLocation
END IF

IF (Local_Error /= -1 .OR. Local_WordsRead /= Local_NumWords ) THEN
!$OMP CRITICAL(internal_write)
  WRITE(io_message,'(A,F5.2,A,I0,A,I0)')                                 &
      'Error in buffin errorCode=',Local_Error,                            &
      ' len=',Local_WordsRead,'/',Local_NumWords
!$OMP END CRITICAL(internal_write)
  CALL io_ereport('io:buffin',abortCode,io_message,                        &
      UNIT=unit)
END IF

IF ( PRESENT(wordsRead) ) wordsRead = Local_WordsRead
IF ( PRESENT(errorCode) ) errorCode = Local_Error

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE buffin64_i


!  Read 64 bit data from the given unit
SUBROUTINE buffin64_l(UNIT,array,numWords,wordsRead,errorCode)

IMPLICIT NONE
INTEGER, INTENT(IN)  :: UNIT         ! Fortran unit number
INTEGER, INTENT(IN),                                                       &
    OPTIONAL         :: numWords     ! no. of words to read in
LOGICAL(KIND=logical64), TARGET,                                               &
    INTENT(OUT)      :: array(:)     ! Array to read into - WARNING: This
                                     ! array must be contiguous for the call to
                                     ! C_LOC to execute correctly
REAL, INTENT(OUT),                                                         &
    OPTIONAL         :: errorCode
INTEGER, INTENT(OUT),                                                      &
    OPTIONAL         :: wordsRead    ! No. of words read

! Local copies of opt. args
INTEGER              :: Local_WordsRead
INTEGER              :: Local_NumWords
REAL                 :: Local_Error
INTEGER              :: mpl_read_err
INTEGER              :: mpl_count_err
INTEGER              :: byteswap_err
INTEGER              :: mpl_read_stat(MPL_STATUS_SIZE)

INTEGER              :: abortCode=25
INTEGER              :: info         ! Broadcast information code
REAL(KIND=jprb)      :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BUFFIN64_L'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (PRESENT(numWords)) THEN
  Local_NumWords=numWords
ELSE
  Local_NumWords=SIZE(array)
END IF
Local_Error=-1.0
Local_WordsRead=Local_NumWords

IF (performOperation(UNIT,reading=.TRUE.) .AND.                             &
    Local_NumWords>SIZE(array)) THEN
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I0,A,I0)')'Request for ',Local_NumWords,           &
      ' words, is not supported by a buffer of size ',SIZE(array)
!$OMP END CRITICAL(internal_write)
  CALL umPrint(umMessage,src='io')
  CALL io_ereport('io:buffin',abortCode,                                   &
      'Supplied buffer too small',UNIT=unit)
END IF

IF (testAttribute(fileAttributes(UNIT)%state,ioAllLocal)) THEN
   IF (testAttribute(fileAttributes(UNIT)%filetype,ioFileTypeMPIIO) .AND.  &
        IOS_Enable_mpiio ) THEN
      IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
      CALL MPL_File_read( mpiio_fh(UNIT),array,Local_numwords,MPL_INTEGER8,&
           mpl_read_stat,mpl_read_err )
! Errors will be handled by the MPI error handler, but set local_error anyway
      IF ( mpl_read_err/=MPL_SUCCESS ) Local_Error=REAL(mpl_read_err)
      CALL MPL_Get_count( mpl_read_stat,MPL_INTEGER8,Local_WordsRead,      &
           mpl_count_err )
      IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()
      IF (get_machine_endianism()==littleEndian)                           &
           byteswap_err = pio_byteswap(C_LOC(array(1)),Local_NumWords,        &
           IOS_BytesPerWord64)
   ELSE
      CALL buffin64_single                                                 &
           (UNIT,array,Local_NumWords,Local_WordsRead,Local_Error)
   END IF
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioLocal)) THEN
  IF (mype==0)                                                             &
      CALL buffin64_single                                                 &
      (UNIT,array,Local_NumWords,Local_WordsRead,Local_Error)
  IF (broadcast_read(UNIT))                                                 &
      CALL gc_ibcast(1,Local_NumWords,0,nproc,info,array)
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioRemote)) THEN
  IF (broadcast_read(UNIT)) THEN
    CALL ios_read64     (UNIT,fileAttributes(UNIT)%nextLocation,           &
        array(1:Local_NumWords),Local_NumWords,IOS_Read_Broadcast)
  ELSE
    CALL ios_read64     (UNIT,fileAttributes(UNIT)%nextLocation,           &
        array(1:Local_NumWords),Local_NumWords,IOS_Read_NoBroadcast)
  END IF
  fileAttributes(UNIT)%nextLocation=ioNoLocation
END IF

IF (Local_Error /= -1 .OR. Local_WordsRead /= Local_NumWords ) THEN
!$OMP CRITICAL(internal_write)
  WRITE(io_message,'(A,F5.2,A,I0,A,I0)')                                 &
      'Error in buffin errorCode=',Local_Error,                            &
      ' len=',Local_WordsRead,'/',Local_NumWords
!$OMP END CRITICAL(internal_write)
  CALL io_ereport('io:buffin',abortCode,io_message,                        &
      UNIT=unit)
END IF

IF ( PRESENT(wordsRead) ) wordsRead = Local_WordsRead
IF ( PRESENT(errorCode) ) errorCode = Local_Error

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE buffin64_l


!  Read 64 bit data from the given unit
SUBROUTINE buffin64_i2D                                                      &
    (UNIT,array,numWords,wordsRead,errorCode)

IMPLICIT NONE
INTEGER, INTENT(IN)  :: UNIT         ! Fortran unit number
INTEGER, INTENT(IN),                                                       &
    OPTIONAL         :: numWords     ! no. of words to read in
INTEGER(KIND=integer64), TARGET,                                               &
    INTENT(OUT)      :: array(:,:)   ! Array to read into
REAL, INTENT(OUT),                                                         &
    OPTIONAL         :: errorCode
INTEGER, INTENT(OUT),                                                      &
    OPTIONAL         :: wordsRead    ! No. of words read

! Local copies of opt. args
INTEGER              :: Local_WordsRead
INTEGER              :: Local_NumWords
REAL                 :: Local_Error
INTEGER              :: mpl_read_err
INTEGER              :: mpl_count_err
INTEGER              :: byteswap_err
INTEGER              :: mpl_read_stat(MPL_STATUS_SIZE)

INTEGER              :: info         ! Broadcast information code
INTEGER              :: abortCode=26
REAL(KIND=jprb)      :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BUFFIN64_I2D'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (PRESENT(numWords)) THEN
  Local_NumWords=numWords
ELSE
  Local_NumWords=SIZE(array)
END IF
Local_Error=-1.0
Local_WordsRead=Local_NumWords

IF (performOperation(UNIT,reading=.TRUE.) .AND.                             &
    Local_NumWords>SIZE(array)) THEN
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I0,A,I0)')'Request for ',Local_NumWords,           &
      ' words, is not supported by a buffer of size ',SIZE(array)
!$OMP END CRITICAL(internal_write)
  CALL umPrint(umMessage,src='io')
  CALL io_ereport('io:buffin',abortCode,                                   &
      'Supplied buffer too small',UNIT=unit)
END IF

IF (testAttribute(fileAttributes(UNIT)%state,ioAllLocal)) THEN
   IF (testAttribute(fileAttributes(UNIT)%filetype,ioFileTypeMPIIO) .AND.  &
        IOS_Enable_mpiio ) THEN
      IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
      CALL MPL_File_read( mpiio_fh(UNIT),array,Local_numwords,MPL_INTEGER8,&
           mpl_read_stat,mpl_read_err )
! Errors will be handled by the MPI error handler, but set local_error anyway
      IF ( mpl_read_err/=MPL_SUCCESS ) Local_Error=REAL(mpl_read_err)
      CALL MPL_Get_count( mpl_read_stat,MPL_INTEGER8,Local_WordsRead,      &
           mpl_count_err )
      IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()
      IF (get_machine_endianism()==littleEndian)                           &
           byteswap_err = pio_byteswap(C_LOC(array),Local_NumWords,        &
           IOS_BytesPerWord64)
   ELSE
      CALL buffin64_single                                                 &
           (UNIT,array,Local_NumWords,Local_WordsRead,Local_Error)
   END IF
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioLocal)) THEN
  IF (mype==0)                                                             &
      CALL buffin64_single                                                 &
      (UNIT,array,Local_NumWords,Local_WordsRead,Local_Error)
  IF (broadcast_read(UNIT))                                                 &
      CALL gc_ibcast(1,Local_NumWords,0,nproc,info,array)
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioRemote)) THEN
  IF (broadcast_read(UNIT)) THEN
    CALL ios_read64     (UNIT,fileAttributes(UNIT)%nextLocation,           &
        array,Local_NumWords,IOS_Read_Broadcast)
  ELSE
    CALL ios_read64     (UNIT,fileAttributes(UNIT)%nextLocation,           &
        array,Local_NumWords,IOS_Read_NoBroadcast)
  END IF
  fileAttributes(UNIT)%nextLocation=ioNoLocation
END IF

IF (Local_Error /= -1 .OR. Local_WordsRead /= Local_NumWords ) THEN
!$OMP CRITICAL(internal_write)
  WRITE(io_message,'(A,F5.2,A,I0,A,I0)')                                 &
      'Error in buffin errorCode=',Local_Error,                            &
      ' len=',Local_WordsRead,'/',Local_NumWords
!$OMP END CRITICAL(internal_write)
  CALL io_ereport('io:buffin',abortCode,io_message,                        &
      UNIT=unit)
END IF

IF ( PRESENT(wordsRead) ) wordsRead = Local_WordsRead
IF ( PRESENT(errorCode) ) errorCode = Local_Error

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE buffin64_i2D


!  Read 64 bit data from the given unit
SUBROUTINE buffin64_l2D                                                      &
    (UNIT,array,numWords,wordsRead,errorCode)

IMPLICIT NONE
INTEGER, INTENT(IN)  :: UNIT         ! Fortran unit number
INTEGER, INTENT(IN),                                                       &
    OPTIONAL         :: numWords     ! no. of words to read in
LOGICAL(KIND=logical64), TARGET,                                               &
    INTENT(OUT)      :: array(:,:)   ! Array to read into - WARNING: This array
                                     ! must be contiguous for the call to C_LOC
                                     ! to execute correctly
REAL, INTENT(OUT),                                                         &
    OPTIONAL         :: errorCode
INTEGER, INTENT(OUT),                                                      &
    OPTIONAL         :: wordsRead    ! No. of words read

! Local copies of opt. args
INTEGER              :: Local_WordsRead
INTEGER              :: Local_NumWords
REAL                 :: Local_Error
INTEGER              :: mpl_read_err
INTEGER              :: mpl_count_err
INTEGER              :: byteswap_err
INTEGER              :: mpl_read_stat(MPL_STATUS_SIZE)

INTEGER              :: info         ! Broadcast information code
INTEGER              :: abortCode=27
REAL(KIND=jprb)      :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BUFFIN64_L2D'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (PRESENT(numWords)) THEN
  Local_NumWords=numWords
ELSE
  Local_NumWords=SIZE(array)
END IF
Local_Error=-1.0
Local_WordsRead=Local_NumWords

IF (performOperation(UNIT,reading=.TRUE.) .AND.                             &
    Local_NumWords>SIZE(array)) THEN
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I0,A,I0)')'Request for ',Local_NumWords,           &
      ' words, is not supported by a buffer of size ',SIZE(array)
!$OMP END CRITICAL(internal_write)
  CALL umPrint(umMessage,src='io')
  CALL io_ereport('io:buffin',abortCode,                                   &
      'Supplied buffer too small',UNIT=unit)
END IF

IF (testAttribute(fileAttributes(UNIT)%state,ioAllLocal)) THEN
      IF (testAttribute(fileAttributes(UNIT)%filetype,ioFileTypeMPIIO) .AND.  &
        IOS_Enable_mpiio ) THEN
      IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
      CALL MPL_File_read( mpiio_fh(UNIT),array,Local_numwords,MPL_INTEGER8,&
           mpl_read_stat,mpl_read_err )
! Errors will be handled by the MPI error handler, but set local_error anyway
      IF ( mpl_read_err/=MPL_SUCCESS ) Local_Error=REAL(mpl_read_err)
      CALL MPL_Get_count( mpl_read_stat,MPL_INTEGER8,Local_WordsRead,      &
           mpl_count_err )
      IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()
      IF (get_machine_endianism()==littleEndian)                           &
           byteswap_err = pio_byteswap(C_LOC(array(1,1)),Local_NumWords,       &
           IOS_BytesPerWord64)
   ELSE
      CALL buffin64_single                                                 &
           (UNIT,array,Local_NumWords,Local_WordsRead,Local_Error)
   END IF
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioLocal)) THEN
  IF (mype==0)                                                             &
      CALL buffin64_single                                                 &
      (UNIT,array,Local_NumWords,Local_WordsRead,Local_Error)
  IF (broadcast_read(UNIT))                                                 &
      CALL gc_ibcast(1,Local_NumWords,0,nproc,info,array)
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioRemote)) THEN
  IF (broadcast_read(UNIT)) THEN
    CALL ios_read64     (UNIT,fileAttributes(UNIT)%nextLocation,           &
        array,Local_NumWords,IOS_Read_Broadcast)
  ELSE
    CALL ios_read64     (UNIT,fileAttributes(UNIT)%nextLocation,           &
        array,Local_NumWords,IOS_Read_NoBroadcast)
  END IF
  fileAttributes(UNIT)%nextLocation=ioNoLocation
END IF

IF (Local_Error /= -1 .OR. Local_WordsRead /= Local_NumWords ) THEN
!$OMP CRITICAL(internal_write)
  WRITE(io_message,'(A,F5.2,A,I0,A,I0)')                                 &
      'Error in buffin errorCode=',Local_Error,                            &
      ' len=',Local_WordsRead,'/',Local_NumWords
!$OMP END CRITICAL(internal_write)
  CALL io_ereport('io:buffin',abortCode,io_message,                        &
      UNIT=unit)
END IF

IF ( PRESENT(wordsRead) ) wordsRead = Local_WordsRead
IF ( PRESENT(errorCode) ) errorCode = Local_Error

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE buffin64_l2D


!  Read 64 bit data from the given unit
SUBROUTINE buffin64_r2D                                                      &
    (UNIT,array,numWords,wordsRead,errorCode)

IMPLICIT NONE
INTEGER, INTENT(IN)  :: UNIT         ! Fortran unit number
INTEGER, INTENT(IN),                                                       &
    OPTIONAL         :: numWords     ! no. of words to read in
REAL(KIND=real64), TARGET,                                                     &
    INTENT(OUT)      :: array(:,:)   ! Array to read into
REAL, INTENT(OUT),                                                         &
    OPTIONAL         :: errorCode
INTEGER, INTENT(OUT),                                                      &
    OPTIONAL         :: wordsRead    ! No. of words read

! Local copies of opt. args
INTEGER              :: Local_WordsRead
INTEGER              :: Local_NumWords
REAL                 :: Local_Error
INTEGER              :: mpl_read_err
INTEGER              :: mpl_count_err
INTEGER              :: byteswap_err
INTEGER              :: mpl_read_stat(MPL_STATUS_SIZE)

INTEGER              :: info         ! Broadcast information code
INTEGER              :: abortCode=28
REAL(KIND=jprb)      :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BUFFIN64_R2D'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (PRESENT(numWords)) THEN
  Local_NumWords=numWords
ELSE
  Local_NumWords=SIZE(array)
END IF
Local_Error=-1.0
Local_WordsRead=Local_NumWords

IF (performOperation(UNIT,reading=.TRUE.) .AND.                             &
    Local_NumWords>SIZE(array)) THEN
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I0,A,I0)')'Request for ',Local_NumWords,           &
      ' words, is not supported by a buffer of size ',SIZE(array)
!$OMP END CRITICAL(internal_write)
  CALL umPrint(umMessage,src='io')
  CALL io_ereport('io:buffin',abortCode,                                   &
      'Supplied buffer too small',UNIT=unit)
END IF

IF (testAttribute(fileAttributes(UNIT)%state,ioAllLocal)) THEN
   IF (testAttribute(fileAttributes(UNIT)%filetype,ioFileTypeMPIIO) .AND.  &
        IOS_Enable_mpiio ) THEN
      IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
      CALL MPL_File_read( mpiio_fh(UNIT),array,Local_numwords,MPL_REAL8,   &
           mpl_read_stat,mpl_read_err )
! Errors will be handled by the MPI error handler, but set local_error anyway
      IF ( mpl_read_err/=MPL_SUCCESS ) Local_Error=REAL(mpl_read_err)
      CALL MPL_Get_count( mpl_read_stat,MPL_REAL8,Local_WordsRead,         &
           mpl_count_err )
      IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()
      IF (get_machine_endianism()==littleEndian)                           &
           byteswap_err = pio_byteswap(C_LOC(array),Local_NumWords,        &
           IOS_BytesPerWord64)
   ELSE
      CALL buffin64_single                                                 &
           (UNIT,array,Local_NumWords,Local_WordsRead,Local_Error)
   END IF
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioLocal)) THEN
  IF (mype==0)                                                             &
      CALL buffin64_single                                                 &
      (UNIT,array,Local_NumWords,Local_WordsRead,Local_Error)
  IF (broadcast_read(UNIT))                                                 &
      CALL gc_rbcast(1,Local_NumWords,0,nproc,info,array)
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioRemote)) THEN
  IF (broadcast_read(UNIT)) THEN
    CALL ios_read64     (UNIT,fileAttributes(UNIT)%nextLocation,           &
        array,Local_NumWords,IOS_Read_Broadcast)
  ELSE
    CALL ios_read64     (UNIT,fileAttributes(UNIT)%nextLocation,           &
        array,Local_NumWords,IOS_Read_NoBroadcast)
  END IF
  fileAttributes(UNIT)%nextLocation=ioNoLocation
END IF

IF (Local_Error /= -1 .OR. Local_WordsRead /= Local_NumWords ) THEN
!$OMP CRITICAL(internal_write)
  WRITE(io_message,'(A,F5.2,A,I0,A,I0)')                                 &
      'Error in buffin errorCode=',Local_Error,                            &
      ' len=',Local_WordsRead,'/',Local_NumWords
!$OMP END CRITICAL(internal_write)
  CALL io_ereport('io:buffin',abortCode,io_message,                        &
      UNIT=unit)
END IF

IF ( PRESENT(wordsRead) ) wordsRead = Local_WordsRead
IF ( PRESENT(errorCode) ) errorCode = Local_Error

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE buffin64_r2D


!  Read 32 bit data from the given unit
!
SUBROUTINE buffin32_r(UNIT,array,numWords,wordsRead,errorCode)

IMPLICIT NONE
INTEGER, INTENT(IN)  :: UNIT         ! Fortran unit number
INTEGER, INTENT(IN),                                                       &
    OPTIONAL         :: numWords     ! no. of words to read in
REAL(KIND=real32), TARGET,                                                     &
    INTENT(OUT)      :: array(:)     ! Array to read into
REAL, INTENT(OUT),                                                         &
    OPTIONAL         :: errorCode
INTEGER, INTENT(OUT),                                                      &
    OPTIONAL         :: wordsRead    ! No. of words read

! Local copies of opt. args
INTEGER              :: Local_WordsRead
INTEGER              :: Local_NumWords
REAL                 :: Local_Error
INTEGER              :: mpl_read_err
INTEGER              :: mpl_count_err
INTEGER              :: byteswap_err
INTEGER              :: mpl_read_stat(MPL_STATUS_SIZE)

INTEGER              :: info         ! Broadcast information code
INTEGER              :: abortCode=29
REAL(KIND=jprb)      :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BUFFIN32_R'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (PRESENT(numWords)) THEN
  Local_NumWords=numWords
ELSE
  Local_NumWords=SIZE(array)
END IF
Local_Error=-1.0
Local_WordsRead=Local_NumWords

IF (performOperation(UNIT,reading=.TRUE.) .AND.                             &
    Local_NumWords>SIZE(array)) THEN
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I0,A,I0)')'Request for ',Local_NumWords,           &
      ' words, is not supported by a buffer of size ',SIZE(array)
!$OMP END CRITICAL(internal_write)
  CALL umPrint(umMessage,src='io')
  CALL io_ereport('io:buffin32',abortCode,                                 &
      'Supplied buffer too small',UNIT=unit)
END IF

IF (testAttribute(fileAttributes(UNIT)%state,ioAllLocal)) THEN
   IF (testAttribute(fileAttributes(UNIT)%filetype,ioFileTypeMPIIO) .AND.  &
        IOS_Enable_mpiio ) THEN
      IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
      CALL MPL_File_read( mpiio_fh(UNIT),array,Local_numwords,MPL_REAL4,   &
           mpl_read_stat,mpl_read_err )
! Errors will be handled by the MPI error handler, but set local_error anyway
      IF ( mpl_read_err/=MPL_SUCCESS ) Local_Error=REAL(mpl_read_err)
      CALL MPL_Get_count( mpl_read_stat,MPL_REAL4,Local_WordsRead,         &
           mpl_count_err )
      IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()
      IF (get_machine_endianism()==littleEndian)                           &
           byteswap_err = pio_byteswap(C_LOC(array),Local_NumWords,        &
           IOS_BytesPerWord32)
   ELSE
      CALL buffin32_single                                                 &
           (UNIT,array,Local_NumWords,Local_WordsRead,Local_Error)
   END IF
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioLocal)) THEN
  IF (mype==0)                                                             &
      CALL buffin32_single                                                 &
      (UNIT,array,Local_NumWords,Local_WordsRead,Local_Error)
  IF (broadcast_read(UNIT))                                                 &
      CALL gc_rbcast(1,Local_NumWords/2,0,nproc,info,array)
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioRemote)) THEN
  IF (broadcast_read(UNIT)) THEN
    CALL ios_read32     (UNIT,fileAttributes(UNIT)%nextLocation,           &
        array(1:Local_NumWords),Local_NumWords,IOS_Read_Broadcast)
  ELSE
    CALL ios_read32     (UNIT,fileAttributes(UNIT)%nextLocation,           &
        array(1:Local_NumWords),Local_NumWords,IOS_Read_NoBroadcast)
  END IF
  fileAttributes(UNIT)%nextLocation=ioNoLocation
END IF

IF (Local_Error /= -1 .OR. Local_WordsRead /= Local_NumWords ) THEN
!$OMP CRITICAL(internal_write)
  WRITE(io_message,'(A,F5.2,A,I0,A,I0)')                                 &
      'Error in buffin errorCode=',Local_Error,                            &
      ' len=',Local_WordsRead,'/',Local_NumWords
!$OMP END CRITICAL(internal_write)
  CALL io_ereport('io:buffin32',abortCode,io_message,                      &
      UNIT=unit)
END IF

IF ( PRESENT(wordsRead) ) wordsRead = Local_WordsRead
IF ( PRESENT(errorCode) ) errorCode = Local_Error

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE buffin32_r


!  Read 32 bit data from the given unit
SUBROUTINE buffin32_i(UNIT,array,numWords,wordsRead,errorCode)
IMPLICIT NONE
INTEGER, INTENT(IN)  :: UNIT         ! Fortran unit number
INTEGER, INTENT(IN),                                                       &
    OPTIONAL         :: numWords     ! no. of words to read in
INTEGER(KIND=integer32), TARGET,                                               &
    INTENT(OUT)      :: array(:)     ! Array to read into
REAL, INTENT(OUT),                                                         &
    OPTIONAL         :: errorCode
INTEGER, INTENT(OUT),                                                      &
    OPTIONAL         :: wordsRead    ! No. of words read

! Local copies of opt. args
INTEGER              :: Local_WordsRead
INTEGER              :: Local_NumWords
REAL                 :: Local_Error
INTEGER              :: mpl_read_err
INTEGER              :: mpl_count_err
INTEGER              :: byteswap_err
INTEGER              :: mpl_read_stat(MPL_STATUS_SIZE)

INTEGER              :: info         ! Broadcast information code
INTEGER              :: abortCode=29
REAL(KIND=jprb)      :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BUFFIN32_I'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (PRESENT(numWords)) THEN
  Local_NumWords=numWords
ELSE
  Local_NumWords=SIZE(array)
END IF
Local_Error=-1.0
Local_WordsRead=Local_NumWords

IF (performOperation(UNIT,reading=.TRUE.) .AND.                             &
    Local_NumWords>SIZE(array)) THEN
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I0,A,I0)')'Request for ',Local_NumWords,           &
      ' words, is not supported by a buffer of size ',SIZE(array)
!$OMP END CRITICAL(internal_write)
  CALL umPrint(umMessage,src='io')
  CALL io_ereport('io:buffin32',abortCode,                                 &
      'Supplied buffer too small',UNIT=unit)
END IF

IF (testAttribute(fileAttributes(UNIT)%state,ioAllLocal)) THEN
   IF (testAttribute(fileAttributes(UNIT)%filetype,ioFileTypeMPIIO) .AND.  &
        IOS_Enable_mpiio ) THEN
      IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
      CALL MPL_File_read( mpiio_fh(UNIT),array,Local_numwords,MPL_INTEGER4,&
           mpl_read_stat,mpl_read_err )
! Errors will be handled by the MPI error handler, but set local_error anyway
      IF ( mpl_read_err/=MPL_SUCCESS ) Local_Error=REAL(mpl_read_err)
      CALL MPL_Get_count( mpl_read_stat,MPL_INTEGER4,Local_WordsRead,      &
           mpl_count_err )
      IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()
      IF (get_machine_endianism()==littleEndian)                           &
           byteswap_err = pio_byteswap(C_LOC(array),Local_NumWords,        &
           IOS_BytesPerWord32)
   ELSE
      CALL buffin32_single                                                 &
           (UNIT,array,Local_NumWords,Local_WordsRead,Local_Error)
   END IF
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioLocal)) THEN
  IF (mype==0)                                                             &
      CALL buffin32_single                                                 &
      (UNIT,array,Local_NumWords,Local_WordsRead,Local_Error)
  IF (broadcast_read(UNIT))                                                 &
      CALL gc_ibcast(1,Local_NumWords/2,0,nproc,info,array)
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioRemote)) THEN
  IF (broadcast_read(UNIT)) THEN
    CALL ios_read32     (UNIT,fileAttributes(UNIT)%nextLocation,           &
        array(1:Local_NumWords),Local_NumWords,IOS_Read_Broadcast)
  ELSE
    CALL ios_read32     (UNIT,fileAttributes(UNIT)%nextLocation,           &
        array(1:Local_NumWords),Local_NumWords,IOS_Read_NoBroadcast)
  END IF
  fileAttributes(UNIT)%nextLocation=ioNoLocation
ELSE
!$OMP CRITICAL(internal_write)
  WRITE(io_message,'(A,I0)')                                              &
      'Error in buffin - unable to resolve filestate ',                    &
      fileAttributes(UNIT)%state
!$OMP END CRITICAL(internal_write)
  CALL io_ereport('io:buffin32',abortCode,io_message,                      &
      UNIT=unit)
END IF

IF (Local_Error /= -1 .OR. Local_WordsRead /= Local_NumWords ) THEN
!$OMP CRITICAL(internal_write)
  WRITE(io_message,'(A,F5.2,A,I0,A,I0)')                                 &
      'Error in buffin errorCode=',Local_Error,                            &
      ' len=',Local_WordsRead,'/',Local_NumWords
!$OMP END CRITICAL(internal_write)
  CALL io_ereport('io:buffin32',abortCode,io_message,                      &
      UNIT=unit)
END IF

IF ( PRESENT(wordsRead) ) wordsRead = Local_WordsRead
IF ( PRESENT(errorCode) ) errorCode = Local_Error

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE buffin32_i


!  Read 32 bit data from the given unit
SUBROUTINE buffin32_i2D(UNIT,array,numWords,wordsRead,errorCode)
IMPLICIT NONE
INTEGER, INTENT(IN)  :: UNIT         ! Fortran unit number
INTEGER, INTENT(IN),                                                       &
    OPTIONAL         :: numWords     ! no. of words to read in
INTEGER(KIND=integer32), TARGET,                                               &
    INTENT(OUT)      :: array(:,:)   ! Array to read into
REAL, INTENT(OUT),                                                         &
    OPTIONAL         :: errorCode
INTEGER, INTENT(OUT),                                                      &
    OPTIONAL         :: wordsRead    ! No. of words read

! Local copies of opt. args
INTEGER              :: Local_WordsRead
INTEGER              :: Local_NumWords
REAL                 :: Local_Error
INTEGER              :: mpl_read_err
INTEGER              :: mpl_count_err
INTEGER              :: byteswap_err
INTEGER              :: mpl_read_stat(MPL_STATUS_SIZE)

INTEGER              :: info         ! Broadcast information code
INTEGER              :: abortCode=29
REAL(KIND=jprb)      :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BUFFIN32_I2D'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (PRESENT(numWords)) THEN
  Local_NumWords=numWords
ELSE
  Local_NumWords=SIZE(array)
END IF
Local_Error=-1.0
Local_WordsRead=Local_NumWords

IF (performOperation(UNIT,reading=.TRUE.) .AND.                             &
    Local_NumWords>SIZE(array)) THEN
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I0,A,I0)')'Request for ',Local_NumWords,           &
      ' words, is not supported by a buffer of size ',SIZE(array)
!$OMP END CRITICAL(internal_write)
  CALL umPrint(umMessage,src='io')
  CALL io_ereport('io:buffin32',abortCode,                                 &
      'Supplied buffer too small',UNIT=unit)
END IF

IF (testAttribute(fileAttributes(UNIT)%state,ioAllLocal)) THEN
   IF (testAttribute(fileAttributes(UNIT)%filetype,ioFileTypeMPIIO) .AND.  &
        IOS_Enable_mpiio ) THEN
      IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
      CALL MPL_File_read( mpiio_fh(UNIT),array,Local_numwords,MPL_INTEGER4,&
           mpl_read_stat,mpl_read_err )
! Errors will be handled by the MPI error handler, but set local_error anyway
      IF ( mpl_read_err/=MPL_SUCCESS ) Local_Error=REAL(mpl_read_err)
      CALL MPL_Get_count( mpl_read_stat,MPL_INTEGER4,Local_WordsRead,      &
           mpl_count_err )
      IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()
      IF (get_machine_endianism()==littleEndian)                           &
           byteswap_err = pio_byteswap(C_LOC(array),Local_NumWords,        &
           IOS_BytesPerWord32)
   ELSE
      CALL buffin32_single                                                 &
           (UNIT,array,Local_NumWords,Local_WordsRead,Local_Error)
   END IF
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioLocal)) THEN
  IF (mype==0)                                                             &
      CALL buffin32_single                                                 &
      (UNIT,array,Local_NumWords,Local_WordsRead,Local_Error)
  IF (broadcast_read(UNIT))                                                 &
      CALL gc_ibcast(1,Local_NumWords/2,0,nproc,info,array)
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioRemote)) THEN
  IF (broadcast_read(UNIT)) THEN
    CALL ios_read32     (UNIT,fileAttributes(UNIT)%nextLocation,           &
        array,Local_NumWords,IOS_Read_Broadcast)
  ELSE
    CALL ios_read32     (UNIT,fileAttributes(UNIT)%nextLocation,           &
        array,Local_NumWords,IOS_Read_NoBroadcast)
  END IF
  fileAttributes(UNIT)%nextLocation=ioNoLocation
END IF

IF (Local_Error /= -1 .OR. Local_WordsRead /= Local_NumWords ) THEN
!$OMP CRITICAL(internal_write)
  WRITE(io_message,'(A,F5.2,A,I0,A,I0)')                                 &
      'Error in buffin errorCode=',Local_Error,                            &
      ' len=',Local_WordsRead,'/',Local_NumWords
!$OMP END CRITICAL(internal_write)
  CALL io_ereport('io:buffin32',abortCode,io_message,                      &
      UNIT=unit)
END IF

IF ( PRESENT(wordsRead) ) wordsRead = Local_WordsRead
IF ( PRESENT(errorCode) ) errorCode = Local_Error

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE buffin32_i2D


!  Read 32 bit data from the given unit
SUBROUTINE buffin32_l(UNIT,array,numWords,wordsRead,errorCode)
IMPLICIT NONE
INTEGER, INTENT(IN)  :: UNIT         ! Fortran unit number
INTEGER, INTENT(IN),                                                       &
    OPTIONAL         :: numWords     ! no. of words to read in
LOGICAL(KIND=logical32), TARGET,                                               &
    INTENT(OUT)      :: array(:)     ! Array to read into - WARNING: This
                                     ! array must be contiguous for the call to
                                     ! C_LOC to execute correctly
REAL, INTENT(OUT),                                                         &
    OPTIONAL         :: errorCode
INTEGER, INTENT(OUT),                                                      &
    OPTIONAL         :: wordsRead    ! No. of words read

! Local copies of opt. args
INTEGER              :: Local_WordsRead
INTEGER              :: Local_NumWords
REAL                 :: Local_Error
INTEGER              :: mpl_read_err
INTEGER              :: mpl_count_err
INTEGER              :: byteswap_err
INTEGER              :: mpl_read_stat(MPL_STATUS_SIZE)


INTEGER              :: info         ! Broadcast information code
INTEGER              :: abortCode=29
REAL(KIND=jprb)      :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BUFFIN32_L'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (PRESENT(numWords)) THEN
  Local_NumWords=numWords
ELSE
  Local_NumWords=SIZE(array)
END IF
Local_Error=-1.0
Local_WordsRead=Local_NumWords

IF (performOperation(UNIT,reading=.TRUE.) .AND.                             &
    Local_NumWords>SIZE(array)) THEN
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I0,A,I0)')'Request for ',Local_NumWords,           &
      ' words, is not supported by a buffer of size ',SIZE(array)
!$OMP END CRITICAL(internal_write)
  CALL umPrint(umMessage,src='io')
  CALL io_ereport('io:buffin32',abortCode,                                 &
      'Supplied buffer too small',UNIT=unit)
END IF

IF (testAttribute(fileAttributes(UNIT)%state,ioAllLocal)) THEN
   IF (testAttribute(fileAttributes(UNIT)%filetype,ioFileTypeMPIIO) .AND.  &
        IOS_Enable_mpiio ) THEN
      IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
      CALL MPL_File_read( mpiio_fh(UNIT),array,Local_numwords,MPL_INTEGER4,&
           mpl_read_stat,mpl_read_err )
! Errors will be handled by the MPI error handler, but set local_error anyway
      IF ( mpl_read_err/=MPL_SUCCESS ) Local_Error=REAL(mpl_read_err)
      CALL MPL_Get_count( mpl_read_stat,MPL_INTEGER4,Local_WordsRead,      &
           mpl_count_err )
      IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()
      IF (get_machine_endianism()==littleEndian)                           &
           byteswap_err = pio_byteswap(C_LOC(array(1)),Local_NumWords,        &
           IOS_BytesPerWord32)
      CALL buffin32_single                                                 &
           (UNIT,array,Local_NumWords,Local_WordsRead,Local_Error)
   END IF
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioLocal)) THEN
  IF (mype==0)                                                             &
      CALL buffin32_single                                                 &
      (UNIT,array,Local_NumWords,Local_WordsRead,Local_Error)
  IF (broadcast_read(UNIT))                                                 &
      CALL gc_ibcast(1,Local_NumWords/2,0,nproc,info,array)
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioRemote)) THEN
  IF (broadcast_read(UNIT)) THEN
    CALL ios_read32     (UNIT,fileAttributes(UNIT)%nextLocation,           &
        array(1:Local_NumWords),Local_NumWords,IOS_Read_Broadcast)
  ELSE
    CALL ios_read32     (UNIT,fileAttributes(UNIT)%nextLocation,           &
        array(1:Local_NumWords),Local_NumWords,IOS_Read_NoBroadcast)
  END IF
  fileAttributes(UNIT)%nextLocation=ioNoLocation
END IF

IF (Local_Error /= -1 .OR. Local_WordsRead /= Local_NumWords ) THEN
!$OMP CRITICAL(internal_write)
  WRITE(io_message,'(A,F5.2,A,I0,A,I0)')                                 &
      'Error in buffin errorCode=',Local_Error,                            &
      ' len=',Local_WordsRead,'/',Local_NumWords
!$OMP END CRITICAL(internal_write)
  CALL io_ereport('io:buffin32',abortCode,io_message,                      &
      UNIT=unit)
END IF

IF ( PRESENT(wordsRead) ) wordsRead = Local_WordsRead
IF ( PRESENT(errorCode) ) errorCode = Local_Error

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE buffin32_l


!  Read character data from the given unit
SUBROUTINE buffin_character                                                  &
    (UNIT,array,numWords,csize,wordsRead,errorCode)

IMPLICIT NONE
INTEGER,  INTENT(IN) :: UNIT            ! Fortran unit number
INTEGER,  INTENT(IN) :: numWords        ! no. of elements to read
INTEGER,  INTENT(IN) :: csize           ! Size of an element
INTEGER, INTENT(OUT),                                                      &
    OPTIONAL         :: wordsRead    ! No. of words read
REAL, INTENT(OUT),                                                         &
    OPTIONAL         :: errorCode
CHARACTER(LEN=csize),                                                      &
    INTENT(OUT)      :: array(numWords) ! Array to read in

! Local copies of opt. args
INTEGER              :: Local_WordsRead ! local copy
REAL                 :: Local_Error     ! local copy

INTEGER              :: mpl_read_err
INTEGER              :: mpl_count_err
INTEGER              :: mpl_read_stat(MPL_STATUS_SIZE)

INTEGER              :: info            ! Broadcast information code
INTEGER              :: abortCode=23
REAL(KIND=jprb)      :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BUFFIN_CHARACTER'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

Local_Error=-1.0
Local_WordsRead=numWords*csize

IF (testAttribute(fileAttributes(UNIT)%state,ioAllLocal)) THEN
   IF (testAttribute(fileAttributes(UNIT)%filetype,ioFileTypeMPIIO) .AND.  &
        IOS_Enable_mpiio ) THEN
      IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
      CALL MPL_File_read( mpiio_fh(UNIT),array,numWords*csize,MPL_BYTE,    &
           mpl_read_stat,mpl_read_err )
! Errors will be handled by the MPI error handler, but set local_error anyway
      IF ( mpl_read_err/=MPL_SUCCESS ) Local_Error=REAL(mpl_read_err)
      CALL MPL_Get_count( mpl_read_stat,MPL_BYTE,Local_WordsRead,          &
           mpl_count_err )
      IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()
   ELSE
      CALL buffin8_single(UNIT,array,numWords*                             &
           csize,Local_WordsRead,Local_Error)
   END IF
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioLocal)) THEN
  IF (mype==0)                                                             &
      CALL buffin8_single(UNIT,array,numWords*                             &
      csize,Local_WordsRead,Local_Error)
  IF (broadcast_read(UNIT))                                                 &
      CALL gc_cbcast(1,numWords*csize,0,nproc,info,array)
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioRemote) .AND.           &
    mype==0) THEN
  CALL io_ereport('io:buffin_character',abortCode,                         &
      'Remote access via buffin8 not supported',UNIT=unit)
END IF

IF ( PRESENT(wordsRead) ) wordsRead = Local_WordsRead
IF ( PRESENT(errorCode) ) errorCode = Local_Error

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE buffin_character


!  Set the file pointer to a position in the file
!  (a 64 bit word address)
SUBROUTINE setpos(UNIT,ipos,error)
USE fort2c_setpos_interfaces, ONLY: setpos_single
IMPLICIT NONE
INTEGER, INTENT(IN)     :: UNIT   ! Fortran unit number
INTEGER, INTENT(IN)     :: ipos   ! Position in file (word aligned)
INTEGER, INTENT(OUT),                                                      &
    OPTIONAL            :: error
INTEGER                 :: Local_Error  ! Return code
INTEGER                 :: abortCode=30
REAL(KIND=jprb)         :: zhook_handle

INTEGER(KIND=MPL_OFFSET_KIND) :: ipos_byte

CHARACTER(LEN=*), PARAMETER :: RoutineName='SETPOS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

Local_Error = 0
CALL checkUnit(UNIT)
IF (testAttribute(fileAttributes(UNIT)%state,ioAllLocal)) THEN
   IF (testAttribute(fileAttributes(UNIT)%filetype,ioFileTypeMPIIO) .AND. &
        IOS_Enable_mpiio ) THEN
      ipos_byte=ipos*IOS_BytesPerInteger
      IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
      CALL MPL_File_seek( mpiio_fh(UNIT),ipos_byte,MPL_SEEK_SET,Local_error )
      IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()
   ELSE
      CALL setpos_single(UNIT,ipos,Local_Error)
   END IF
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioLocal) .AND.            &
    mype==0) THEN
  CALL setpos_single(UNIT,ipos,Local_Error)
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioRemote) .AND.           &
    mype==0) THEN
  fileAttributes(UNIT)%nextLocation=ipos
  !      CALL IOS_Setpos   (unit,ipos,Local_Error)
END IF

IF ( PRESENT(error) ) error = Local_Error

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE setpos


!  Set the file pointer to a position in the file
!  (an 8 bit word address)
SUBROUTINE setpos8(UNIT,ipos,error)
USE fort2c_setpos_interfaces, ONLY: setpos8_single
IMPLICIT NONE
INTEGER, INTENT(IN)     :: UNIT   ! Fortran unit number
INTEGER, INTENT(IN)     :: ipos   ! Position in file (word aligned)
INTEGER, INTENT(OUT),                                                      &
    OPTIONAL            :: error
INTEGER                 :: Local_Error  ! Return code
INTEGER                 :: abortCode=30
REAL(KIND=jprb)         :: zhook_handle

INTEGER(KIND=MPL_OFFSET_KIND) :: ipos_byte

CHARACTER(LEN=*), PARAMETER :: RoutineName='SETPOS8'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

Local_Error = 0
CALL checkUnit(UNIT)
IF (testAttribute(fileAttributes(UNIT)%state,ioAllLocal)) THEN
   IF (testAttribute(fileAttributes(UNIT)%filetype,ioFileTypeMPIIO) .AND. &
        IOS_Enable_mpiio ) THEN
      ipos_byte=ipos
      IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
      CALL MPL_File_seek( mpiio_fh(UNIT),ipos_byte,MPL_SEEK_SET,Local_error )
      IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()
   ELSE
      CALL setpos8_single(UNIT,ipos,Local_Error)
   END IF
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioLocal) .AND.            &
    mype==0) THEN
  CALL setpos8_single(UNIT,ipos,Local_Error)
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioRemote) .AND.           &
    mype==0) THEN
  CALL io_ereport('io:setpos8',abortcode,                                  &
      'setpos8 not implemented for IO Servers!')
  !      fileAttributes(unit)%nextLocation=ipos
  !      CALL IOS_Setpos   (unit,ipos,Local_Error)
END IF

IF ( PRESENT(error) ) error = Local_Error

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE setpos8


!  Set the file pointer to a position in the file
!  (a 32 bit word address)
SUBROUTINE setpos32(UNIT,ipos,error)
USE fort2c_setpos_interfaces, ONLY: setpos32_single
IMPLICIT NONE
INTEGER, INTENT(IN)     :: UNIT   ! Fortran unit number
INTEGER, INTENT(IN)     :: ipos   ! Position in file (word aligned)
INTEGER, INTENT(OUT),                                                      &
    OPTIONAL            :: error
INTEGER                 :: abortCode=30
INTEGER                 :: Local_Error
REAL(KIND=jprb)         :: zhook_handle

INTEGER(KIND=MPL_OFFSET_KIND) :: ipos_byte

CHARACTER(LEN=*), PARAMETER :: RoutineName='SETPOS32'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

Local_Error = 0
CALL checkUnit(UNIT)
IF (testAttribute(fileAttributes(UNIT)%state,ioAllLocal)) THEN
   IF (testAttribute(fileAttributes(UNIT)%filetype,ioFileTypeMPIIO) .AND. &
        IOS_Enable_mpiio ) THEN
      ipos_byte=ipos*IOS_BytesPerWord32
      IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
      CALL MPL_File_seek( mpiio_fh(UNIT),ipos_byte,MPL_SEEK_SET,Local_error)
      IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()
   ELSE
      CALL setpos32_single(UNIT,ipos,Local_Error)
   END IF
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioLocal) .AND.            &
    mype==0) THEN
  CALL setpos32_single(UNIT,ipos,Local_Error)
ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioRemote) .AND.           &
    mype==0) THEN
  CALL io_ereport('io:setpos32',abortcode,                                 &
      'setpos32 not implemented for IO Servers!')
  !      fileAttributes(unit)%nextLocation=ipos
  !      CALL IOS_Setpos   (unit,ipos,Local_Error)
END IF

IF ( PRESENT(error) ) error = Local_Error

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE setpos32


!  Set a bit in a vaiable
SUBROUTINE setAttribute(var,mask)

IMPLICIT NONE
INTEGER, INTENT(INOUT) :: var  ! Variable
INTEGER, INTENT(IN)    :: mask ! Bitmask to set
var=IOR(var,mask)
END SUBROUTINE setAttribute


!  Clear a bit in a variable
SUBROUTINE clearAttribute(var,mask)

IMPLICIT NONE
INTEGER, INTENT(INOUT) :: var  ! Variable
INTEGER, INTENT(IN)    :: mask ! Bitmask to clear
INTEGER, PARAMETER     :: minus_one=-1
INTEGER                :: tmp
tmp=IEOR(minus_one,mask)
var=IAND(var,tmp)
END SUBROUTINE clearAttribute


!  Test whether a bit is set in a variable
LOGICAL FUNCTION testAttribute(var,mask)

IMPLICIT NONE
INTEGER, INTENT(IN)    :: var  ! Variable
INTEGER, INTENT(IN)    :: mask ! Bitmask to compare
testAttribute=.FALSE.
IF (IAND(var,mask)==mask) testAttribute=.TRUE.
RETURN
END FUNCTION testAttribute


!  Cause a unit to NOT broadcast read data to all tasks
SUBROUTINE set_unit_bcast_flag(UNIT)

IMPLICIT NONE
INTEGER, INTENT(IN) :: UNIT
CALL checkUnit(UNIT)
CALL clearAttribute(fileAttributes(UNIT)%state,ioReadReplicate)
END SUBROUTINE set_unit_bcast_flag


!  Cause a unit to broadcast read data to all tasks
SUBROUTINE clear_unit_bcast_flag(UNIT)

IMPLICIT NONE
INTEGER, INTENT(IN) :: UNIT
CALL checkUnit(UNIT)
CALL setAttribute(fileAttributes(UNIT)%state,ioReadReplicate)
END SUBROUTINE clear_unit_bcast_flag


!  Enquire as to whether a unit will broadcast read data
LOGICAL FUNCTION broadcast_read(UNIT)

IMPLICIT NONE
INTEGER, INTENT(IN) :: UNIT
CALL checkUnit(UNIT)
broadcast_read=                                                            &
    testAttribute(fileAttributes(UNIT)%state,ioReadReplicate)
RETURN
END FUNCTION broadcast_read


!  Pass a timestep event to the subsystem
!   (important for coupled models)
SUBROUTINE io_timestep(ts)
  ! Timestep notification

IMPLICIT NONE
INTEGER, INTENT(IN):: ts
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='IO_TIMESTEP'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (l_ios_active() .AND. mype==0) CALL IOS_Process(ts)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE io_timestep


!  Commit data to the disk subsystem with optional wait on completion
!   (important for coupled models)
SUBROUTINE ioDiskSynchronise(UNIT,wait)

IMPLICIT NONE
INTEGER, INTENT(IN) :: UNIT      ! Unit to sync to disk
LOGICAL, OPTIONAL   :: wait      ! If true, the caller should block until
                                 ! the opertaion comletes. If the unit
                                 ! is open locally then this will have
                                 ! little effect, but remote IO units
                                 ! setting this flag will see a significant

INTEGER             :: icode     ! Local error code
LOGICAL             :: L_wait    ! local copy of optional argument
REAL(KIND=jprb)     :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='IODISKSYNCHRONISE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

L_Wait=.FALSE.
IF (PRESENT(wait)) THEN
  IF (wait)L_wait=wait
END IF

! As this can be expensive, lets have some documentation that it happened
IF (printstatus>=prstatus_oper) THEN
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I3,A)')'IO: Synchronising unit: ',UNIT,' with disk.'
!$OMP END CRITICAL(internal_write)
  CALL umPrint(umMessage,src='io')
  IF (l_wait) THEN
!$OMP CRITICAL(internal_write)
    WRITE(umMessage,'(A,I3)')'IO: Disk sync blocking until completion...'
!$OMP END CRITICAL(internal_write)
    CALL umPrint(umMessage,src='io')
  END IF
END IF

IF (testAttribute(fileAttributes(UNIT)%state,ioOpen)) THEN
  IF (testAttribute(fileAttributes(UNIT)%state,ioAllLocal)) THEN
    CALL sync_single(UNIT,icode)             ! Flush unix layers
  ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioLocal)               &
      .AND. mype==0) THEN
    CALL sync_single(UNIT,icode)
  ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioRemote)              &
      .AND. mype==0) THEN
     CALL IOS_DiskSync(UNIT,fileAttributes(UNIT)%fileType,l_wait)
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ioDiskSynchronise


!  Generate a fence operation. Define a model epoch such as to make guarantees
!  about IO operation ordering
!
!  Calling this is equivilent to the statement "Nothing further will happen
!  regarding unit 'unit', until all other units instructions issued before
!  point this have been completed". The epoch is described to have passed
!  when this condition is satisfied.
!
!  Implementation notes:
!
!   Operations issued for other units after the epoch may still complete
!   before the epoch completes.
!
!   Atmosphere does not know that the epoch has completed.
!
!   If 'unit' is open locally then Atmosphere will stall waiting for the
!   epoch to pass. Therefore it is not recommended to use this method
!   on a local file if IO servers are active and processing other IO channels.
!
!   If 'unit' is remote then the IO server responsible will stall. This may
!   cause performance problems (ie a stall in Atmosphere) if his queue becomes
!   full, and cannot receive new work from the model.

SUBROUTINE ioFence(UNIT)
IMPLICIT NONE
INTEGER, INTENT(IN) :: UNIT
INTEGER             :: root
INTEGER             :: errorcode=41
REAL(KIND=jprb)     :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='IOFENCE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (testAttribute(fileAttributes(UNIT)%state,ioOpen)) THEN
  IF (testAttribute(fileAttributes(UNIT)%state,ioAllLocal)) THEN
    CALL io_ereport('io:Fence',errorcode,                                  &
        'Cannot perform fence operation on a unit with many writers')
  ELSE IF (mype==0) THEN
    IF (testAttribute(fileAttributes(UNIT)%state,ioLocal)) THEN
      root=mype
    ELSE
      root=-1
    END IF
    CALL IOS_Fence(UNIT,root)
  END IF

END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE ioFence

INTEGER FUNCTION readableWords(UNIT,wordLength)
IMPLICIT NONE
INTEGER, INTENT(IN) :: UNIT
INTEGER, OPTIONAL   :: wordLength
INTEGER             :: local_wordLength
TYPE(fileState)     :: state
REAL(KIND=jprb)     :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READABLEWORDS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! If not specified we'll assume they want default type lengths
local_wordLength=umFortranRealSize()

IF ( PRESENT(wordLength) ) local_wordLength=wordLength

CALL ioFileState(UNIT,state)
readableWords=(state%fileExtent-state%filePosition)/                       &
    local_wordLength

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION readableWords

SUBROUTINE ioFileState(UNIT,state)
IMPLICIT NONE
INTEGER, INTENT(IN) :: UNIT
INTEGER             :: info ! generic error code
TYPE(fileState)     :: state
TYPE(IOS_Status)    :: rstate
REAL(KIND=jprb)     :: zhook_handle

INTEGER(KIND=MPL_OFFSET_KIND) :: filePosition
INTEGER(KIND=MPL_OFFSET_KIND) :: fileExtent

CHARACTER(LEN=*), PARAMETER :: RoutineName='IOFILESTATE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

state%UNIT         = UNIT
state%fileExtent   = ioNotKnown
state%filePosition = ioNotKnown
info               = 0

IF (testAttribute(fileAttributes(UNIT)%state,ioOpen)) THEN
  IF (testAttribute(fileAttributes(UNIT)%state,ioAllLocal)) THEN
     IF (testAttribute(fileAttributes(UNIT)%filetype,ioFileTypeMPIIO) .AND. &
          IOS_Enable_mpiio ) THEN
        IF ( IOS_serialise_mpi_calls ) CALL acquire_lock_mpi()
        CALL MPL_File_get_position(mpiio_fh(UNIT),filePosition,info)
        CALL MPL_File_get_size(mpiio_fh(UNIT),fileExtent,info)
        IF ( IOS_serialise_mpi_calls ) CALL release_lock_mpi()
        state%fileExtent = fileExtent
        state%filePosition = filePosition
     ELSE
        CALL getpos8  (UNIT,state%filePosition)
        CALL getextent(UNIT,state%fileExtent,info)
     END IF
  ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioLocal)               &
      .AND. mype==0) THEN
    CALL getpos8  (UNIT,state%filePosition)
    CALL getextent(UNIT,state%fileExtent,info)
  ELSE IF (testAttribute(fileAttributes(UNIT)%state,ioRemote)              &
      .AND. mype==0) THEN
    CALL IOS_fileState(UNIT,fileAttributes(UNIT)%nextLocation,rstate)
    state%filePosition =rstate%POSITION
    state%fileExtent   =rstate%extent
  END IF
END IF

IF (.NOT. testAttribute(fileAttributes(UNIT)%state,ioAllLocal)) THEN
  IF (broadcast_read(UNIT)) THEN
    CALL gc_ibcast(1,filestate_len,0,nproc,info,state)
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ioFileState

CHARACTER(LEN=maxFileNameLength) FUNCTION ReturnFileName(UNIT)
IMPLICIT NONE
INTEGER, INTENT(IN) :: UNIT

ReturnFileName=fileAttributes(UNIT)%filename

END FUNCTION ReturnFileName

LOGICAL FUNCTION isReadOnly(UNIT)
IMPLICIT NONE
INTEGER, INTENT(IN) :: UNIT

isReadOnly=testAttribute(fileAttributes(UNIT)%state,ioReadOnly)

END FUNCTION isReadOnly

LOGICAL FUNCTION isAllLocal(unit)
IMPLICIT NONE
INTEGER, INTENT(IN) :: unit

isAllLocal=testAttribute(fileAttributes(unit)%state,ioAllLocal)

END FUNCTION isAllLocal

END MODULE io
