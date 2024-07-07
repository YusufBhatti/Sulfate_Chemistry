! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: C96

! Module defining buffering constants and namelists

MODULE io_configuration_mod

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength
USE missing_data_mod, ONLY: imdi

IMPLICIT NONE

INTEGER, PARAMETER :: io_namelist_unit = 127
INTEGER, PARAMETER :: io_timing_on     = 1
INTEGER, PARAMETER :: io_timing_off    = 0

INTEGER :: io_wbuffer_size    = 524288    ! In words. Size of buffer when
                                          ! buffered I/O is active.
                                          ! Default is 4MB.

INTEGER :: io_rbuffer_size    = 524288    ! In words. Size of buffer when
                                          ! buffered I/O is active.
                                          ! Default is 4MB.

INTEGER :: io_rbuffer_count    = 4        ! Number of read buffers/stream

INTEGER :: io_rbuffer_update   = 1        !
INTEGER :: io_rbuffer_prefetch = 1        !

! In words. Alignment of data block in
INTEGER :: io_data_alignment = 524288
                                          ! dumps/fieldsfiles. Defaults is 4MB.

INTEGER :: io_field_padding  = 512        ! Previously um_sector_size. The
                                          ! padding between  fields. Now mostly
                                          ! defunct due to buffered I/O.

LOGICAL :: io_omp_passive_dump_reading = .FALSE. 
                                          ! Use OMP passive waiting
                                          ! policy for dump reading.

! ---------------------------------------------------------------------------- !
! Profile numbers for filesystem specific settings.
! Can have the following values: -
!
! 0 / no_filesystem_profile     = no profile used (no custom settings required)
! 1 / lustre_filesystem_profile = LUSTRE file system (tunings for file striping)
!
! other values currently have no meaning

INTEGER, PARAMETER :: no_filesystem_profile     = 0
INTEGER, PARAMETER :: lustre_filesystem_profile = 1

INTEGER :: io_filesystem_profile = imdi

INTEGER, PARAMETER :: io_filesystem_profile_max = 1

! ---------------------------------------------------------------------------- !

INTEGER :: io_timing  = io_timing_off

LOGICAL :: io_external_control = .FALSE.

LOGICAL :: io_use_unique_units = .FALSE.

LOGICAL :: io_alltoall_readflds = .FALSE.

LOGICAL :: print_runtime_info = .FALSE.

LOGICAL :: l_postp = .FALSE.              ! Control of ".arch" files

! A namelist will allow us to control these via the gui.
NAMELIST /io_control/                                                          &
    io_wbuffer_size,                                                           &
    io_rbuffer_size,                                                           &
    io_rbuffer_count,                                                          &
    io_rbuffer_update,                                                         &
    io_rbuffer_prefetch,                                                       &
    io_data_alignment,                                                         &
    io_field_padding,                                                          &
    io_timing,                                                                 &
    io_filesystem_profile,                                                     &
    io_use_unique_units,                                                       &
    io_alltoall_readflds,                                                      &
    print_runtime_info,                                                        &
    io_external_control,                                                       &
    l_postp,                                                                   &
    io_omp_passive_dump_reading

! DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

! ---------------------------------------------------------------------------- !

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='IO_CONFIGURATION_MOD'

CONTAINS

! ---------------------------------------------------------------------------- !

SUBROUTINE ioLoadConfig()
USE FilenameLength_mod, ONLY:                                             &
    filenamelength
USE ereport_mod
USE application_description, ONLY: isSmallExec
USE um_parcore, ONLY: mype
USE mpl, ONLY: mpl_integer
USE file_manager, ONLY: portio_unique_units
USE get_env_var_mod, ONLY: get_env_var
USE umPrintMgr, ONLY: PrintStatus, PrStatus_Normal

IMPLICIT NONE
INTEGER :: length  ! Length of env var contents
INTEGER :: error
INTEGER :: errCode
CHARACTER (LEN=FileNameLength) :: NamelistFile
CHARACTER (LEN=errormessagelength) :: message
CHARACTER (LEN=errormessagelength) :: iomessage
INTEGER :: icode
INTEGER :: my_comm


NamelistFile = "dummy file"
errCode=-10

CALL gc_get_communicator(my_comm, icode)

CALL get_env_var("IOSCNTL",NameListFile, allow_missing=.TRUE., length=length)

IF (length > 0) THEN
  IF (mype == 0) OPEN(IO_Namelist_Unit,FILE=NameListFile,FORM='FORMATTED', &
      STATUS='OLD', ACTION='READ',IOSTAT=error, IOMSG=iomessage)
  CALL mpl_bcast(error,1,mpl_integer,0,my_comm,icode)
  IF ( error == 0 ) THEN
    ! Read the I/O buffering control namelist too
    CALL read_nml_io_control(IO_Namelist_Unit)
    IF (mype == 0) CLOSE(IO_Namelist_Unit)

    IF (PrintStatus >= PrStatus_Normal .AND. mype == 0) THEN
      CALL print_nlist_io_control()
    END IF

    ! Set the io control option in file manager
    portio_unique_units = io_use_unique_units

    CALL parse_io_filesystem_profile()
  ELSE
    ! Stifle warnings from small execs which caused user concern.
    ! Small execs don't need to find the file as default values work
    ! perfectly well.
    IF (mype == 0) THEN
      IF (.NOT. isSmallExec()) THEN
        WRITE(message,'(A)')                                           &
          'Failed to open IO control file ('                         //&
          TRIM(NameListFile) // ')' //                        newline//&
          'IoMsg: '//TRIM(iomessage)
        CALL ereport('ioLoadConfig',errCode,message)
      END IF
    END IF
  END IF
END IF

END SUBROUTINE ioLoadConfig

! ---------------------------------------------------------------------------- !

SUBROUTINE read_nml_io_control(unit_in)

USE um_parcore, ONLY: mype

USE check_iostat_mod, ONLY: check_iostat

USE setup_namelist, ONLY: setup_nml_type

IMPLICIT NONE

INTEGER, INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_IO_CONTROL'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 2
INTEGER, PARAMETER :: n_int = 9
INTEGER, PARAMETER :: n_log = 6

TYPE my_namelist
  SEQUENCE
  INTEGER :: io_wbuffer_size
  INTEGER :: io_rbuffer_size
  INTEGER :: io_rbuffer_count
  INTEGER :: io_rbuffer_update
  INTEGER :: io_rbuffer_prefetch
  INTEGER :: io_data_alignment
  INTEGER :: io_field_padding
  INTEGER :: io_timing
  INTEGER :: io_filesystem_profile
  LOGICAL :: io_use_unique_units
  LOGICAL :: print_runtime_info
  LOGICAL :: io_external_control
  LOGICAL :: l_postp
  LOGICAL :: io_omp_passive_dump_reading
  LOGICAL :: io_alltoall_readflds
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,   &
                    n_log_in=n_log)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=io_control, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "NAMELIST io_control", iomessage)

  my_nml % io_wbuffer_size       = io_wbuffer_size
  my_nml % io_rbuffer_size       = io_rbuffer_size
  my_nml % io_rbuffer_count      = io_rbuffer_count
  my_nml % io_rbuffer_update     = io_rbuffer_update
  my_nml % io_rbuffer_prefetch   = io_rbuffer_prefetch
  my_nml % io_data_alignment     = io_data_alignment
  my_nml % io_field_padding      = io_field_padding
  my_nml % io_timing             = io_timing
  my_nml % io_filesystem_profile = io_filesystem_profile
  my_nml % io_use_unique_units   = io_use_unique_units
  my_nml % print_runtime_info    = print_runtime_info
  my_nml % io_external_control   = io_external_control
  my_nml % l_postp               = l_postp
  my_nml % io_omp_passive_dump_reading = io_omp_passive_dump_reading
  my_nml % io_alltoall_readflds  = io_alltoall_readflds

END IF

CALL mpl_bcast(my_nml,1,MPL_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  io_wbuffer_size       = my_nml % io_wbuffer_size
  io_rbuffer_size       = my_nml % io_rbuffer_size
  io_rbuffer_count      = my_nml % io_rbuffer_count
  io_rbuffer_update     = my_nml % io_rbuffer_update
  io_rbuffer_prefetch   = my_nml % io_rbuffer_prefetch
  io_data_alignment     = my_nml % io_data_alignment
  io_field_padding      = my_nml % io_field_padding
  io_timing             = my_nml % io_timing
  io_filesystem_profile = my_nml % io_filesystem_profile
  io_use_unique_units   = my_nml % io_use_unique_units
  print_runtime_info    = my_nml % print_runtime_info
  io_external_control   = my_nml % io_external_control
  l_postp               = my_nml % l_postp
  io_omp_passive_dump_reading = my_nml % io_omp_passive_dump_reading
  io_alltoall_readflds  = my_nml % io_alltoall_readflds

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_io_control

! ---------------------------------------------------------------------------- !

SUBROUTINE parse_io_filesystem_profile()

USE umPrintMgr, ONLY:                                                          &
  umPrint,                                                                     &
  PrintStatus,                                                                 &
  umMessage,                                                                   &
  PrStatus_Normal

USE ereport_mod, ONLY:                                                         &
  ereport

USE lustre_control_mod, ONLY:                                                  &
  Load_Lustre_Config,                                                          &
  Print_Lustre_Config

IMPLICIT NONE

INTEGER :: errCode
CHARACTER (LEN=errormessagelength) :: message
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PARSE_IO_FILESYSTEM_PROFILE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Load the profile specific settings
SELECT CASE (io_filesystem_profile)

   CASE (no_filesystem_profile)
      IF (PrintStatus>PrStatus_Normal) THEN
        WRITE(umMessage,'(A)') '[INFO] no IO Filesystem Profile in use'
        CALL umPrint(umMessage)
      END IF

   CASE (lustre_filesystem_profile)
      CALL Load_Lustre_Config()
      IF (PrintStatus>PrStatus_Normal) THEN
        WRITE(umMessage,'(A)') '[INFO] Using Lustre IO Filesystem Profile'
        CALL umPrint(umMessage)
        CALL Print_Lustre_Config()
      END IF

   CASE DEFAULT
     errCode=10
!$OMP CRITICAL(internal_write)
     WRITE(message,'(A,I0,A,I0,A,I0)')                                         &
                                  'Unrecognised io_filesystem_profile value: ',&
                                  io_filesystem_profile,                       &
                                  ' is not in the range ',                     &
                                  no_filesystem_profile, ' to ',               &
                                  io_filesystem_profile_max
!$OMP END CRITICAL(internal_write)

     CALL ereport('parse_io_filesystem_profile',errCode,message)

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE

! ---------------------------------------------------------------------------- !

SUBROUTINE release_io_filesystem_profile()

USE lustre_control_mod, ONLY:                                                  &
  Release_Lustre_Config

IMPLICIT NONE

! unload the profile specific settings
SELECT CASE (io_filesystem_profile)

   CASE (no_filesystem_profile)
     ! Nothing to do for no profile

   CASE (lustre_filesystem_profile)
      CALL Release_Lustre_Config()

END SELECT

END SUBROUTINE

SUBROUTINE print_nlist_io_control()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_IO_CONTROL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist io_control', &
    src='io_configuration_mod')

WRITE(lineBuffer,'(A,I0)') ' io_wbuffer_size             = ',io_wbuffer_size
CALL umPrint(lineBuffer,src='io_configuration_mod')
WRITE(lineBuffer,'(A,I0)') ' io_rbuffer_size             = ',io_rbuffer_size
CALL umPrint(lineBuffer,src='io_configuration_mod')
WRITE(lineBuffer,'(A,I0)') ' io_rbuffer_count            = ',io_rbuffer_count
CALL umPrint(lineBuffer,src='io_configuration_mod')
WRITE(lineBuffer,'(A,I0)') ' io_rbuffer_update           = ',io_rbuffer_update
CALL umPrint(lineBuffer,src='io_configuration_mod')
WRITE(lineBuffer,'(A,I0)') ' io_rbuffer_prefetch         = ',io_rbuffer_prefetch
CALL umPrint(lineBuffer,src='io_configuration_mod')
WRITE(lineBuffer,'(A,I0)') ' io_data_alignment           = ',io_data_alignment
CALL umPrint(lineBuffer,src='io_configuration_mod')
WRITE(lineBuffer,'(A,I0)') ' io_field_padding            = ',io_field_padding
CALL umPrint(lineBuffer,src='io_configuration_mod')
WRITE(lineBuffer,'(A,I0)') ' io_timing                   = ',io_timing
CALL umPrint(lineBuffer,src='io_configuration_mod')
WRITE(lineBuffer,'(A,I0)') ' io_filesystem_profile       = ',                  &
                           io_filesystem_profile
CALL umPrint(lineBuffer,src='io_configuration_mod')
WRITE(lineBuffer,'(A,L1)')' io_use_unique_units         = ',io_use_unique_units
CALL umPrint(lineBuffer,src='io_configuration_mod')
WRITE(lineBuffer,'(A,L1)')' print_runtime_info          = ',print_runtime_info
CALL umPrint(lineBuffer,src='io_configuration_mod')
WRITE(lineBuffer,'(A,L1)')' io_external_control         = ',io_external_control
CALL umPrint(lineBuffer,src='io_configuration_mod')
WRITE(lineBuffer,'(A,L1)')' l_postp                     = ',l_postp
CALL umPrint(lineBuffer,src='io_configuration_mod')
WRITE(lineBuffer,'(A,L1)') ' io_omp_passive_dump_reading = ',                  &
                           io_omp_passive_dump_reading
CALL umPrint(lineBuffer,src='io_configuration_mod')
WRITE(lineBuffer,'(A,L1)') ' io_alltoall_readflds        = ',                  &
                           io_alltoall_readflds
CALL umPrint(lineBuffer,src='io_configuration_mod')

CALL umPrint('- - - - - - end of namelist - - - - - -', &
    src='io_configuration_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_io_control

END MODULE io_configuration_mod
