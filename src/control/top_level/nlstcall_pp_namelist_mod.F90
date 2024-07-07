! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description: Module based interface to the nlstcall_pp namelist, 
!              which configures output streams and files
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!             This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 90

MODULE nlstcall_pp_namelist_mod

USE filenamelength_mod, ONLY: filenamelength
USE errormessagelength_mod, ONLY: errormessagelength
USE iso_fortran_env, ONLY: iostat_end
USE file_manager, ONLY: &
    assign_file_unit,   &
    release_file_unit,  &
    um_file_type,       &
    get_file_by_unit,   &
    init_file_loop,     &
    get_file_by_id
USE nlstcall_subs_mod, ONLY: totime, real_months, timesteps
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE profilename_length_mod, ONLY: fileid_length

IMPLICIT NONE 

! Read variables for nlstcall_pp
CHARACTER(LEN=filenamelength) :: filename
CHARACTER(LEN=fileid_length)  :: file_id
CHARACTER(LEN=filenamelength) :: filename_base

LOGICAL :: l_reinit         ! Re-initialisation
INTEGER :: reinit_unit      ! Re-init time unit
INTEGER :: reinit_end       ! Re-init end time
INTEGER :: reinit_step      ! Re-init frequency
INTEGER :: reinit_start     ! Re-init start time
INTEGER :: reserved_headers ! Reserved headers
INTEGER :: packing          ! Packing code

! Used only when read from HISTORY file
LOGICAL :: partial    ! Partially written
INTEGER :: last_field ! Last field written



! DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

NAMELIST / nlstcall_pp / filename,         &
                         file_id,          &
                         l_reinit,         &
                         reinit_end,       &
                         reinit_step,      &
                         reinit_start,     &
                         reinit_unit,      &
                         filename_base,    &
                         reserved_headers, &
                         packing

NAMELIST / nlstcall_pp_hist / filename,    &
                              file_id,     &
                              last_field

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='NLSTCALL_PP_NAMELIST_MOD'

CONTAINS

! This version of the routine is intended to read the history file
!-------------------------------------------------------------------------------
SUBROUTINE read_nml_nlstcall_pp_hist(unit_in)

USE um_parcore, ONLY: mype
USE check_iostat_mod, ONLY: check_iostat
USE setup_namelist, ONLY: setup_nml_type
USE umprintmgr, ONLY: printstatus, prstatus_oper

IMPLICIT NONE 

INTEGER, INTENT(IN)  :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
INTEGER :: pp_unit
INTEGER :: outtime

CHARACTER(LEN=errormessagelength) :: iomessage

REAL(KIND=jprb)             :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_NLSTCALL_PP_HIST'
TYPE(um_file_type), POINTER :: um_file

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 2
INTEGER, PARAMETER :: n_int = 2
INTEGER, PARAMETER :: n_chars = filenamelength + fileid_length

TYPE my_namelist
  SEQUENCE
  INTEGER :: last_field
  INTEGER :: error
  CHARACTER (LEN=filenamelength) :: filename
  CHARACTER (LEN=fileid_length)  :: file_id
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,        &
                    n_chars_in=n_chars)

! Start of loop over namelists
DO

  ! Setup defaults 
  last_field = 0
  filename   = ""
  file_id    = ""

  IF (mype == 0) THEN

    READ(UNIT=unit_in, NML=nlstcall_pp_hist, IOSTAT=ErrorStatus,      &
         IOMSG=iomessage)
    ! Don't check for negative error status (EOF indicator)
    IF (ErrorStatus /= iostat_end) THEN
      CALL check_iostat(errorstatus, "namelist NLSTCALL_PP_HIST", iomessage)
    END IF

    my_nml % last_field  = last_field
    my_nml % filename    = filename
    my_nml % file_id     = file_id
    my_nml % error       = ErrorStatus

  END IF

  CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

  IF (mype /= 0) THEN

    last_field       = my_nml % last_field
    filename         = my_nml % filename
    file_id          = my_nml % file_id
    ErrorStatus      = my_nml % error

  END IF

  ! If this attempt to read the namelist hit the end of the file, stop
  ! (note ErrorStatus is broadcast WITH the nml variables to ensure all 
  ! PEs know when further reads are required)
  IF (ErrorStatus == iostat_end) EXIT

  IF (printstatus >= prstatus_oper .AND. mype == 0) THEN
    CALL print_nlist_nlstcall_pp_hist()
  END IF

  ! Create a new file object for this stream
  NULLIFY(um_file)
  CALL assign_file_unit(filename, pp_unit, &
                        handler="portio", id=file_id)
  um_file => get_file_by_unit(pp_unit, handler="portio")

  ! Set the flag to indicate it is partially written
  um_file % meta % partially_written = .TRUE.
  um_file % pp_meta % last_written_field = last_field

END DO

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_nlstcall_pp_hist

!-------------------------------------------------------------------------------
SUBROUTINE print_nlist_nlstcall_pp_hist()

USE umprintMgr, ONLY: umprint

IMPLICIT NONE

CHARACTER(LEN=50000) :: linebuffer
REAL(KIND=jprb)      :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_NLSTCALL_PP_HIST'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umprint('Contents of namelist nlstcall_pp_hist', src=modulename)

WRITE(linebuffer,"(A,I0)")' last_field = ', last_field
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,A)")' filename = ', TRIM(filename)
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,A)")' file_id = ', TRIM(file_id)
CALL umprint(linebuffer,src=modulename)

CALL umprint('- - - - - - end of namelist - - - - - -', src=modulename)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_nlstcall_pp_hist

! This routine will read as many nlstcall_pp namelists as it finds in the
! file attached to the provided unit
!-------------------------------------------------------------------------------
SUBROUTINE read_nml_nlstcall_pp(unit_in)

USE um_parcore, ONLY: mype
USE check_iostat_mod, ONLY: check_iostat
USE setup_namelist, ONLY: setup_nml_type
USE umprintmgr, ONLY: printstatus, prstatus_oper

IMPLICIT NONE

INTEGER, INTENT(IN)  :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
INTEGER :: pp_unit
INTEGER :: outtime

CHARACTER(LEN=errormessagelength) :: iomessage

REAL(KIND=jprb)             :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_NLSTCALL_PP'
TYPE(um_file_type), POINTER :: um_file

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int = 7
INTEGER, PARAMETER :: n_log = 1
INTEGER, PARAMETER :: n_chars = filenamelength*2 + fileid_length

TYPE my_namelist
  SEQUENCE
  INTEGER :: reinit_end
  INTEGER :: reinit_step
  INTEGER :: reinit_start
  INTEGER :: reinit_unit
  INTEGER :: reserved_headers
  INTEGER :: packing
  INTEGER :: error
  LOGICAL :: l_reinit
  CHARACTER (LEN=filenamelength) :: filename_base
  CHARACTER (LEN=filenamelength) :: filename
  CHARACTER (LEN=fileid_length)  :: file_id
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,        &
                    n_log_in=n_log, n_chars_in=n_chars)

! Start of loop over namelists
DO

  ! Setup defaults 
  reinit_end       = 0
  reinit_step      = 0
  reinit_start     = 0
  reinit_unit      = 1
  filename_base    = ""
  reserved_headers = 0
  packing          = 0
  l_reinit         = .FALSE.
  filename         = ""
  file_id          = ""

  IF (mype == 0) THEN

    READ(UNIT=unit_in, NML=nlstcall_pp, IOSTAT=ErrorStatus, IOMSG=iomessage)
    ! Don't check for negative error status (EOF indicator)
    IF (ErrorStatus /= iostat_end) THEN
      CALL check_iostat(errorstatus, "namelist NLSTCALL_PP", iomessage)
    END IF

    my_nml % reinit_end       = reinit_end
    my_nml % reinit_step      = reinit_step
    my_nml % reinit_start     = reinit_start
    my_nml % reinit_unit      = reinit_unit
    my_nml % filename_base    = filename_base
    my_nml % reserved_headers = reserved_headers
    my_nml % packing          = packing
    my_nml % l_reinit         = l_reinit
    my_nml % filename         = filename
    my_nml % file_id          = file_id
    my_nml % error            = ErrorStatus

  END IF

  CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

  IF (mype /= 0) THEN

    reinit_end       = my_nml % reinit_end
    reinit_step      = my_nml % reinit_step
    reinit_start     = my_nml % reinit_start
    reinit_unit      = my_nml % reinit_unit
    filename_base    = my_nml % filename_base
    reserved_headers = my_nml % reserved_headers
    packing          = my_nml % packing
    l_reinit         = my_nml % l_reinit
    filename         = my_nml % filename
    file_id          = my_nml % file_id
    ErrorStatus      = my_nml % error

  END IF

  ! If this attempt to read the namelist hit the end of the file, stop
  ! (note ErrorStatus is broadcast WITH the nml variables to ensure all 
  ! PEs know when further reads are required)
  IF (ErrorStatus == iostat_end) EXIT

  IF (printstatus >= prstatus_oper .AND. mype == 0) THEN
    CALL print_nlist_nlstcall_pp()
  END IF
  CALL check_nlstcall_pp()

  ! Otherwise all PEs should now initialise an object to hold the stream
  ! It's possible that an earlier read of the history namelist will have
  ! already initialised any streams which are going to continue as part
  ! of a CRUN, in these cases we will find the object already exists
  NULLIFY(um_file)
  um_file => get_file_by_id(file_id, handler="portio", &
                            ignore_missing=.TRUE.)
  IF (.NOT. ASSOCIATED(um_file)) THEN
    ! If the file object doesn't already exist the file will be new
    ! to this run, so initialise the object here instead.

    ! If the file is tagged for re-initialisation it won't be
    ! provided with a filename - so provide a placeholder for
    ! any logging messages in the meantime
    IF (l_reinit) THEN
      filename = "Reserved unit for re-initialised stream ("//&
                  TRIM(filename_base)//")"
    END IF
    ! Create a new file object for this stream
    CALL assign_file_unit(filename, pp_unit, &
                          handler="portio", id=file_id)
    um_file => get_file_by_unit(pp_unit, handler="portio")
  END IF

  ! Simple settings and flags from the namelist
  um_file % meta % packing_code = packing
  um_file % pp_meta % reserved_headers = reserved_headers
  um_file % pp_meta % init_file_type = "p"
  um_file % meta % is_output_file = .TRUE.

  ! Re-initialisation settings (if used)
  IF (l_reinit) THEN
    um_file % meta % filename_base = filename_base
    IF ((reinit_unit) == real_months) THEN
      um_file % meta % init_steps      = -reinit_step
      um_file % meta % init_start_step =  reinit_start
      um_file % meta % init_end_step   =  reinit_end
    ELSE
      CALL totime(reinit_step, reinit_unit, outtime)
      um_file % meta % init_steps = outtime
      CALL totime(reinit_start, reinit_unit, outtime)
      um_file % meta % init_start_step = outtime
      CALL totime(reinit_end, reinit_unit, outtime)
      um_file % meta % init_end_step = outtime
    END IF
  END IF

  ! Special settings for PPVAR stream - although the streams are designed to
  ! be generic, the VAR system refuses to change to accomodate the UM output
  ! so here we inject a bit of an awful hack to detect a stream specifically 
  ! named as the stream used by VAR and force a few of its settings to 
  ! non-standard values to preserve the old behaviour
  !
  ! The result of these changes are that the var stream's first file in the 
  ! sequence will be shifted by one timestep, causing it to initialise in the
  ! previous re-initialisation period.  It's type letter will also change from
  ! "p" to "c" which will appear in the filenames and also force the first file
  ! to be named as if no shift had taken place
  IF (TRIM(file_id) == "ppvar" .OR. TRIM(file_id) == "ppmbc") THEN
    um_file % meta % init_start_step = &
        um_file % meta % init_start_step - 1
    um_file % meta % init_end_step = &
        um_file % meta % init_end_step - 1
    um_file % pp_meta % init_file_type = "c"
  END IF

END DO

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_nlstcall_pp

!-------------------------------------------------------------------------------
SUBROUTINE print_nlist_nlstcall_pp()

USE umprintMgr, ONLY: umprint

IMPLICIT NONE

CHARACTER(LEN=50000) :: linebuffer
REAL(KIND=jprb)      :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_NLSTCALL_PP'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umprint('Contents of namelist nlstcall_pp', src=modulename)

WRITE(linebuffer,"(A,I0)")' reinit_end   = ', reinit_end
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,I0)")' reinit_step  = ', reinit_step
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,I0)")' reinit_start = ', reinit_start
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,I0)")' reinit_unit  = ', reinit_unit
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,A)") ' filename_base    = ', TRIM(filename_base)
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,I0)")' reserved_headers = ', reserved_headers
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,I0)")' packing  = ', packing
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,L1)")' l_reinit = ', l_reinit
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,A)") ' filename = ', TRIM(filename)
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,A)") ' file_id  = ', TRIM(file_id)
CALL umprint(linebuffer,src=modulename)

CALL umprint('- - - - - - end of namelist - - - - - -', src=modulename)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_nlstcall_pp


!-------------------------------------------------------------------------------
SUBROUTINE check_nlstcall_pp()
! Description:
!   Subroutine to apply logic controls and set control variables based on the
!   options selected in the nlstcall_pp  namelist.
!
! Dr Hook Modules
USE yomhook,           ONLY: lhook, dr_hook
USE parkind1,          ONLY: jprb, jpim

USE chk_opts_mod,      ONLY: chk_var, def_src

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CHECK_NLSTCALL_PP'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
def_src = RoutineName

CALL chk_var( packing,            'packing',                 '[0,1,2,4,5,7]' )
CALL chk_var( reinit_unit,        'reinit_unit',             '[0,1,2,3,4]' )

def_src = ""
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_nlstcall_pp

!-------------------------------------------------------------------------------
! The inverse of the above routine - write current output streams to
! a file, for use in history
!-------------------------------------------------------------------------------
SUBROUTINE write_nml_nlstcall_pp_hist(unit_in)

USE um_parcore, ONLY: mype
USE ereport_mod, ONLY: ereport

IMPLICIT NONE 

INTEGER, INTENT(IN)  :: unit_in

INTEGER :: icode
REAL(KIND=jprb)             :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='WRITE_NML_NLSTCALL_PP_HIST'
TYPE(um_file_type), POINTER :: um_file

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! This should never be called by a PE other than PE0 but protect it 
! anyway just to be sure
IF (mype == 0) THEN
  
  ! Loop over all binary/non-portable files
  NULLIFY(um_file)
  um_file => init_file_loop(handler="portio")
  DO WHILE (ASSOCIATED(um_file))

    IF (um_file % meta % is_output_file) THEN
      ! Partially written files must be indicated in the namelist to
      ! allow a CRUN to pick up from the final written position
      IF (um_file % meta % partially_written) THEN
        filename    = um_file % filename
        file_id     = um_file % id
        last_field  = um_file % pp_meta % last_written_field

        ! Namelist is now populated and ready for output
        WRITE(unit_in, nlstcall_pp_hist, IOSTAT=icode)
        CALL ereport('write_nml_nlstcall_pp_hist',icode,                       &
                     'Write ERROR for NLSTCALL_PP_HIST namelist')
      END IF
    END IF
    ! Increment file pointer for next pass of loop
    um_file => um_file % next
  END DO
ELSE
  ! This shouldn't be called by other PEs, but warn instead of crashing since
  ! it isn't a critical mistake
  icode = -10
  CALL ereport('write_nml_nlstcall_pp_hist', icode,                            &
      'Should not be called by PEs other than PE0')
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE write_nml_nlstcall_pp_hist

END MODULE nlstcall_pp_namelist_mod

