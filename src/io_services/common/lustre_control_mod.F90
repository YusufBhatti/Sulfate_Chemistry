! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: io_services

! Module defining Lustre control settings

MODULE lustre_control_mod

USE yomhook,  ONLY:                                                            &
  lhook,                                                                       &
  dr_hook

USE parkind1, ONLY:                                                            &
  jprb,                                                                        &
  jpim

USE filenamelength_mod, ONLY:                                                  &
  filenamelength

USE errormessagelength_mod, ONLY:                                              &
  errormessagelength

USE file_manager, ONLY:                                                        &
  assign_file_unit,                                                            &
  release_file_unit

USE missing_data_mod, ONLY:                                                    &
  imdi

USE IOS_Constants, ONLY:                                                       &
  IOS_info_strlen

IMPLICIT NONE
PRIVATE

INTEGER :: default_stripe_size    = imdi
INTEGER :: default_stripe_offset  = imdi
INTEGER :: default_stripe_count   = imdi
INTEGER :: default_stripe_pattern = imdi

INTEGER :: number_custom_files = imdi

CHARACTER(LEN=filenamelength), ALLOCATABLE :: striped_file(:)
INTEGER, ALLOCATABLE                       :: stripe_size(:)
INTEGER, ALLOCATABLE                       :: stripe_offset(:)
INTEGER, ALLOCATABLE                       :: stripe_count(:)
INTEGER, ALLOCATABLE                       :: stripe_pattern(:)


! A namelist will allow us to control these via the gui.
NAMELIST /lustre_control/                                                      &
    default_stripe_size,                                                       &
    default_stripe_offset,                                                     &
    default_stripe_count,                                                      &
    default_stripe_pattern,                                                    &
    number_custom_files

NAMELIST /lustre_control_custom_files/                                         &
    stripe_size,                                                               &
    stripe_offset,                                                             &
    stripe_count,                                                              &
    stripe_pattern,                                                            &
    striped_file

! DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

PUBLIC ::                                                                      &
  Load_Lustre_Config,                                                          &
  Release_Lustre_Config,                                                       &
  Print_Lustre_Config,                                                         &
  apply_file_striping_mpiio

! ---------------------------------------------------------------------------- !

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LUSTRE_CONTROL_MOD'

CONTAINS

! ---------------------------------------------------------------------------- !

SUBROUTINE Load_Lustre_Config()

USE ereport_mod

USE get_env_var_mod, ONLY: get_env_var

USE um_parcore, ONLY:                                                          &
  mype

USE mpl, ONLY:                                                                 &
  mpl_integer

IMPLICIT NONE

INTEGER :: length   ! Length of env var contents
INTEGER :: error
INTEGER :: errCode
CHARACTER (LEN=FileNameLength) :: NamelistFile
CHARACTER (LEN=errormessagelength) :: message
CHARACTER (LEN=errormessagelength) :: iomessage
INTEGER :: icode
INTEGER :: my_comm
INTEGER :: Lustre_Namelist_Unit
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LOAD_LUSTRE_CONFIG'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

NamelistFile = "dummy file"
errCode=10

CALL gc_get_communicator(my_comm, icode)

CALL get_env_var("IOSCNTL",NameListFile, allow_missing=.TRUE., length=length)

IF (length > 0) THEN

  IF (mype == 0) THEN
    CALL assign_file_unit(NameListFile, Lustre_Namelist_Unit, HANDLER="fortran")
    OPEN(Lustre_Namelist_Unit, FILE=NameListFile, FORM='FORMATTED',            &
         STATUS='OLD', IOSTAT=error, IOMSG=iomessage, ACTION='READ')
  END IF

  CALL mpl_bcast(error,1,mpl_integer,0,my_comm,icode)

  IF ( error == 0 ) THEN

    CALL read_nml_lustre_control(Lustre_Namelist_Unit)

    CALL read_nml_lustre_control_custom_files(Lustre_Namelist_Unit)

    IF (mype == 0) THEN
      CLOSE(Lustre_Namelist_Unit)
      CALL release_file_unit(Lustre_Namelist_Unit, HANDLER="fortran")
    END IF

  ELSE

    IF (mype == 0) THEN
        WRITE(message,'(A,A,A,A)')'Failed to open Lustre control file (',      &
          TRIM(NameListFile),')\nError: ', TRIM(iomessage)
        CALL ereport('Load_Lustre_Config',errCode,message)
    END IF
  END IF

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE Load_Lustre_Config

! ---------------------------------------------------------------------------- !

SUBROUTINE read_nml_lustre_control(unit_in)

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

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_LUSTRE_CONTROL'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: num_lustre_cntl_nml_types = 1
INTEGER, PARAMETER :: n_int = 5

TYPE lustre_cntl_nml_type
  SEQUENCE
  INTEGER :: default_stripe_size
  INTEGER :: default_stripe_offset
  INTEGER :: default_stripe_count
  INTEGER :: default_stripe_pattern
  INTEGER :: number_custom_files
END TYPE lustre_cntl_nml_type

TYPE (lustre_cntl_nml_type) :: lustre_cntl_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

! Read in default settings

CALL setup_nml_type(num_lustre_cntl_nml_types, mpl_nml_type, n_int_in=n_int)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=lustre_control, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "NAMELIST lustre_control", iomessage)

  lustre_cntl_nml % default_stripe_size    = default_stripe_size
  lustre_cntl_nml % default_stripe_offset  = default_stripe_offset
  lustre_cntl_nml % default_stripe_count   = default_stripe_count
  lustre_cntl_nml % default_stripe_pattern = default_stripe_pattern
  lustre_cntl_nml % number_custom_files    = number_custom_files

END IF

CALL mpl_bcast(lustre_cntl_nml,1,MPL_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  default_stripe_size    = lustre_cntl_nml % default_stripe_size
  default_stripe_offset  = lustre_cntl_nml % default_stripe_offset
  default_stripe_count   = lustre_cntl_nml % default_stripe_count
  default_stripe_pattern = lustre_cntl_nml % default_stripe_pattern
  number_custom_files    = lustre_cntl_nml % number_custom_files

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_lustre_control

! ---------------------------------------------------------------------------- !

SUBROUTINE read_nml_lustre_control_custom_files(unit_in)

USE um_parcore, ONLY: mype
USE check_iostat_mod, ONLY: check_iostat
USE mpl, ONLY: mpl_packed, mpl_integer, mpl_character

IMPLICIT NONE

INTEGER, INTENT(IN) :: unit_in

INTEGER :: my_comm
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: &
  RoutineName='READ_NML_LUSTRE_CONTROL_CUSTOM_FILES'

INTEGER :: n_int
INTEGER :: n_chars
INTEGER :: packsize_int, packsize_char
INTEGER :: position
CHARACTER, ALLOCATABLE :: buffer(:)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (number_custom_files>0) THEN

  n_int = (number_custom_files*4)
  n_chars = (number_custom_files*filenamelength)

  CALL gc_get_communicator(my_comm, icode)

  ! Allocate arrays for file specific settings

  ALLOCATE (striped_file(number_custom_files))
  ALLOCATE (stripe_size(number_custom_files))
  ALLOCATE (stripe_offset(number_custom_files))
  ALLOCATE (stripe_count(number_custom_files))
  ALLOCATE (stripe_pattern(number_custom_files))

  CALL mpl_pack_size(n_int, mpl_integer, my_comm, packsize_int, icode)
  CALL mpl_pack_size(n_chars, mpl_character, my_comm, packsize_char, icode)

  ! Allocate an extra 10% to buffer for safety
  ALLOCATE(buffer(INT(1.1*(packsize_int+packsize_char))))

  IF (mype == 0) THEN

    ! Now read in file specific settings
    READ(UNIT=unit_in, NML=lustre_control_custom_files, IOSTAT=ErrorStatus,    &
         IOMSG=iomessage)
    CALL check_iostat(errorstatus, "NAMELIST lustre_control_custom_files",     &
                      iomessage)

    position = 0

    CALL mpl_pack(stripe_size, number_custom_files, mpl_integer, buffer,       &
                  SIZE(buffer), position, my_comm, icode)
    CALL mpl_pack(stripe_offset, number_custom_files, mpl_integer, buffer,     &
                  SIZE(buffer), position, my_comm, icode)
    CALL mpl_pack(stripe_count, number_custom_files, mpl_integer, buffer,      &
                  SIZE(buffer), position, my_comm, icode)
    CALL mpl_pack(stripe_pattern, number_custom_files, mpl_integer, buffer,    &
                  SIZE(buffer), position, my_comm, icode)
    CALL mpl_pack(striped_file, n_chars, mpl_character, buffer,                &
                  SIZE(buffer), position, my_comm, icode)

  END IF

  CALL mpl_bcast(buffer, SIZE(buffer), mpl_packed, 0, my_comm, icode)

  IF (mype /= 0) THEN

    position = 0

    CALL mpl_unpack(buffer, SIZE(buffer), position, stripe_size,               &
                    number_custom_files, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, SIZE(buffer), position, stripe_offset,             &
                    number_custom_files, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, SIZE(buffer), position, stripe_count,              &
                    number_custom_files, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, SIZE(buffer), position, stripe_pattern,            &
                    number_custom_files, mpl_integer, my_comm, icode)
    CALL mpl_unpack(buffer, SIZE(buffer), position, striped_file,              &
                    n_chars, mpl_character, my_comm, icode)

  END IF

  DEALLOCATE (buffer)

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_lustre_control_custom_files

! ---------------------------------------------------------------------------- !

SUBROUTINE Release_Lustre_Config()

IMPLICIT NONE

IF( ALLOCATED(striped_file) ) THEN

  DEALLOCATE (stripe_pattern)
  DEALLOCATE (stripe_count)
  DEALLOCATE (stripe_offset)
  DEALLOCATE (stripe_size)
  DEALLOCATE (striped_file)

END IF

END SUBROUTINE Release_Lustre_Config

! ---------------------------------------------------------------------------- !

SUBROUTINE Print_Lustre_Config()

USE umprintmgr, ONLY:                                                          &
  umprint,                                                                     &
  umMessage

IMPLICIT NONE

INTEGER :: i
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_LUSTRE_CONFIG'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP CRITICAL(internal_write)
WRITE(umMessage,'(A)') ''
!$OMP END CRITICAL(internal_write)
CALL umprint(umMessage)
!$OMP CRITICAL(internal_write)
WRITE(umMessage,'(A)') '[INFO] *** LUSTRE CONFIG PRINT ***'
!$OMP END CRITICAL(internal_write)
CALL umprint(umMessage)
!$OMP CRITICAL(internal_write)
WRITE(umMessage,'(A)') ''
!$OMP END CRITICAL(internal_write)
CALL umprint(umMessage)
!$OMP CRITICAL(internal_write)
WRITE(umMessage,'(A,I0)') 'default_stripe_size    = ', default_stripe_size
!$OMP END CRITICAL(internal_write)
CALL umprint(umMessage)
!$OMP CRITICAL(internal_write)
WRITE(umMessage,'(A,I0)') 'default_stripe_offset  = ', default_stripe_offset
!$OMP END CRITICAL(internal_write)
CALL umprint(umMessage)
!$OMP CRITICAL(internal_write)
WRITE(umMessage,'(A,I0)') 'default_stripe_count   = ', default_stripe_count
!$OMP END CRITICAL(internal_write)
CALL umprint(umMessage)
!$OMP CRITICAL(internal_write)
WRITE(umMessage,'(A,I0)') 'default_stripe_pattern = ', default_stripe_pattern
!$OMP END CRITICAL(internal_write)
CALL umprint(umMessage)
!$OMP CRITICAL(internal_write)
WRITE(umMessage,'(A,I0)') 'number_custom_files    = ', number_custom_files
!$OMP END CRITICAL(internal_write)
CALL umprint(umMessage)
!$OMP CRITICAL(internal_write)
WRITE(umMessage,'(A)') ''
!$OMP END CRITICAL(internal_write)
CALL umprint(umMessage)

DO i = 1 ,number_custom_files

!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I0,A,A)') 'File ', i, ': ', TRIM(striped_file(i))
!$OMP END CRITICAL(internal_write)
  CALL umprint(umMessage)
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I0)') 'stripe_size            = ', stripe_size(i)
!$OMP END CRITICAL(internal_write)
  CALL umprint(umMessage)
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I0)') 'stripe_offset          = ', stripe_offset(i)
!$OMP END CRITICAL(internal_write)
  CALL umprint(umMessage)
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I0)') 'stripe_count           = ', stripe_count(i)
!$OMP END CRITICAL(internal_write)
  CALL umprint(umMessage)
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I0)') 'stripe_pattern         = ', stripe_pattern(i)
!$OMP END CRITICAL(internal_write)
  CALL umprint(umMessage)
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A)') ''
!$OMP END CRITICAL(internal_write)
  CALL umprint(umMessage)

END DO

!$OMP CRITICAL(internal_write)
WRITE(umMessage,'(A)') '*** END LUSTRE CONFIG PRINT ***'
!$OMP END CRITICAL(internal_write)
CALL umprint(umMessage)
!$OMP CRITICAL(internal_write)
WRITE(umMessage,'(A)') ''
!$OMP END CRITICAL(internal_write)
CALL umprint(umMessage)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE Print_Lustre_Config

! ---------------------------------------------------------------------------- !

SUBROUTINE c_get_file_striping(c_stripe_count, c_stripe_size, c_stripe_offset, &
                               c_stripe_pattern, c_filename)                   &
                               BIND(c, NAME="f_get_file_striping")

USE, INTRINSIC :: iso_c_binding, ONLY:                                         &
  c_int64_t,                                                                   &
  C_CHAR

USE f_shum_string_conv_mod, ONLY:                                              &
  f_shum_c2f_string

IMPLICIT NONE

CHARACTER(KIND=C_CHAR,LEN=1), INTENT(IN) ::                                    &
  c_filename(*)

INTEGER(KIND=c_int64_t), INTENT(OUT) ::                                        &
  c_stripe_count,                                                              &
  c_stripe_size,                                                               &
  c_stripe_offset,                                                             &
  c_stripe_pattern

CHARACTER(LEN=filenamelength) ::                                               &
  lookup_filename


INTEGER ::                                                                     &
  lookup_stripe_count,                                                         &
  lookup_stripe_size,                                                          &
  lookup_stripe_offset,                                                        &
  lookup_stripe_pattern

  lookup_filename = f_shum_c2f_string(c_filename)

  CALL get_file_striping( lookup_stripe_count, lookup_stripe_size,           &
                          lookup_stripe_offset, lookup_stripe_pattern,       &
                          lookup_filename )

  c_stripe_count = lookup_stripe_count
  c_stripe_size = lookup_stripe_size
  c_stripe_offset = lookup_stripe_offset
  c_stripe_pattern = lookup_stripe_pattern

END SUBROUTINE c_get_file_striping

SUBROUTINE apply_file_striping_mpiio( lookup_filename, mpiio_info )

USE umprintmgr, ONLY:                                                          &
  umprint,                                                                     &
  umMessage,                                                                   &
  PrintStatus,                                                                 &
  PrStatus_Normal

IMPLICIT NONE

  CHARACTER(len=filenamelength), INTENT(IN) ::                                 &
       lookup_filename

  INTEGER, INTENT(IN)                       ::                                 &
       mpiio_info

  INTEGER ::                                                                   &
       lookup_stripe_count,                                                    &
       lookup_stripe_size,                                                     &
       lookup_stripe_offset,                                                   &
       lookup_stripe_pattern

  CHARACTER(len=IOS_info_strlen) ::                                            &
       info_val_str

  INTEGER ::                                                                   &
       icode,                                                                  &
       str_len

  LOGICAL ::                                                                   &
       defined

  CALL get_file_striping( lookup_stripe_count, lookup_stripe_size,             &
                          lookup_stripe_offset, lookup_stripe_pattern,         &
                          lookup_filename )

  WRITE(info_val_str,'(i0)') lookup_stripe_count
  CALL MPL_Info_set( mpiio_info, 'striping_factor', TRIM(info_val_str), icode )
  WRITE(info_val_str,'(i0)') lookup_stripe_size
  CALL MPL_Info_set( mpiio_info, 'striping_unit', TRIM(info_val_str), icode )
  WRITE(info_val_str,'(i0)') lookup_stripe_offset
  CALL MPL_Info_set( mpiio_info, 'start_iodevice', TRIM(info_val_str), icode )
  ! stripe_pattern not supported.

  IF (PrintStatus>PrStatus_Normal) THEN
!$OMP CRITICAL(internal_write)
     WRITE(umMessage,'(A,A)') 'Lustre settings for MPI-IO file ',       &
          TRIM(lookup_filename)
!$OMP END CRITICAL(internal_write)
     CALL umprint(umMessage)
     CALL MPL_Info_get( mpiio_info, 'striping_factor', str_len, info_val_str, &
                        defined, icode )
!$OMP CRITICAL(internal_write)
     WRITE(umMessage,'(A)') 'stripe_count: '//TRIM(info_val_str)
!$OMP END CRITICAL(internal_write)
     CALL umprint(umMessage)
     CALL MPL_Info_get( mpiio_info, 'striping_unit', str_len, info_val_str,   &
                        defined, icode )
!$OMP CRITICAL(internal_write)
     WRITE(umMessage,'(A,I0)') 'stripe_size: '//TRIM(info_val_str)
!$OMP END CRITICAL(internal_write)
     CALL umprint(umMessage)
     CALL MPL_Info_get( mpiio_info, 'start_iodevice', str_len, info_val_str,  &
                        defined, icode )
!$OMP CRITICAL(internal_write)
     WRITE(umMessage,'(A,I0)') 'stripe_offset: '//TRIM(info_val_str)
!$OMP END CRITICAL(internal_write)
     CALL umprint(umMessage)
  END IF

END SUBROUTINE apply_file_striping_mpiio


SUBROUTINE get_file_striping( lookup_stripe_count, lookup_stripe_size,         &
                              lookup_stripe_offset, lookup_stripe_pattern,     &
                              lookup_filename )

USE umPrintMgr, ONLY:                                                          &
  umPrint,                                                                     &
  PrintStatus,                                                                 &
  umMessage,                                                                   &
  PrStatus_Normal

IMPLICIT NONE

  CHARACTER(LEN=filenamelength), INTENT(IN) ::                                 &
       lookup_filename

  INTEGER, INTENT(OUT) ::                                                      &
       lookup_stripe_count,                                                    &
       lookup_stripe_size,                                                     &
       lookup_stripe_offset,                                                   &
       lookup_stripe_pattern

  INTEGER :: i

  ! If we have no custom files, the results will obviously be the defaults
  IF ( number_custom_files == 0 ) THEN

    lookup_stripe_count = default_stripe_count
    lookup_stripe_size = default_stripe_size
    lookup_stripe_offset = default_stripe_offset
    lookup_stripe_pattern = default_stripe_pattern

    RETURN

 END IF

 IF (PrintStatus>PrStatus_Normal) THEN
!$OMP CRITICAL(internal_write)
    WRITE(umMessage,'(A,A)') '[info] querying lustre settings for file: ',     &
         TRIM(lookup_filename)
!$OMP END CRITICAL(internal_write)
    CALL umprint(umMessage)
 END IF

  ! Try to find if we have a match for this file
 DO i=1,number_custom_files

!$OMP CRITICAL(internal_write)
    WRITE(umMessage,'(A,A)') '[info] Trying file: ',     &
         TRIM(striped_file(i))
!$OMP END CRITICAL(internal_write)
    CALL umprint(umMessage)

    IF(INDEX(lookup_filename,TRIM(striped_file(i)))==0) CYCLE

    lookup_stripe_count = stripe_count(i)
    lookup_stripe_size = stripe_size(i)
    lookup_stripe_offset = stripe_offset(i)
    lookup_stripe_pattern = stripe_pattern(i)

    IF (PrintStatus>PrStatus_Normal) THEN
!$OMP CRITICAL(internal_write)
      WRITE(umMessage,'(A,A)') 'Found custom Lustre settings for file ',       &
          TRIM(lookup_filename)
!$OMP END CRITICAL(internal_write)
      CALL umprint(umMessage)

!$OMP CRITICAL(internal_write)
      WRITE(umMessage,'(A,I0)') 'stripe_count: ', lookup_stripe_count
!$OMP END CRITICAL(internal_write)
      CALL umprint(umMessage)
!$OMP CRITICAL(internal_write)
      WRITE(umMessage,'(A,I0)') 'stripe_size: ', lookup_stripe_size
!$OMP END CRITICAL(internal_write)
      CALL umprint(umMessage)
!$OMP CRITICAL(internal_write)
      WRITE(umMessage,'(A,I0)') 'stripe_offset: ', lookup_stripe_offset
!$OMP END CRITICAL(internal_write)
      CALL umprint(umMessage)
!$OMP CRITICAL(internal_write)
      WRITE(umMessage,'(A,I0)') 'stripe_pattern: ', lookup_stripe_pattern
!$OMP END CRITICAL(internal_write)
      CALL umprint(umMessage)
    END IF

    RETURN

  END DO

  ! If we get this far, there is no custom striping for this file
  ! ...so just use the defaults.
  lookup_stripe_count = default_stripe_count
  lookup_stripe_size = default_stripe_size
  lookup_stripe_offset = default_stripe_offset
  lookup_stripe_pattern = default_stripe_pattern
  
  IF (PrintStatus>PrStatus_Normal) THEN
!$OMP CRITICAL(internal_write)
    WRITE(umMessage,'(A,A,A)') 'Cannot find custom Lustre settings for file ', &
                               TRIM(lookup_filename), ' so using defaults.'
!$OMP END CRITICAL(internal_write)
    CALL umprint(umMessage)
 END IF

END SUBROUTINE get_file_striping

END MODULE lustre_control_mod
