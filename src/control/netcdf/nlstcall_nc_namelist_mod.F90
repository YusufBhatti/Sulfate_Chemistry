! *****************************COPYRIGHT*******************************
! 
! Copyright 2017-2018 University of Reading
! 
! Redistribution and use in source and binary forms, with or without 
! modification, are permitted provided that the following conditions are met:
! 
! 1. Redistributions of source code must retain the above copyright notice, this
! list of conditions and the following disclaimer.
! 
! 2. Redistributions in binary form must reproduce the above copyright notice, 
! this list of conditions and the following disclaimer in the documentation 
! and/or other materials provided with the distribution.
! 
! 3. Neither the name of the copyright holder nor the names of its contributors 
! may be used to endorse or promote products derived from this software without 
! specific prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! 
! *****************************COPYRIGHT*******************************
!
! Description: Module based interface to the nlstcall_nc namelist, 
!              which configures netCDF format output streams and files
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!             This file belongs in section: NetCDF output
!
! Code Description:
!   Language: FORTRAN 90

MODULE nlstcall_nc_namelist_mod

USE filenamelength_mod, ONLY: filenamelength
USE errormessagelength_mod, ONLY: errormessagelength
USE missing_data_mod, ONLY: imdi
USE file_manager, ONLY: &
    assign_file_unit,   &
    um_file_type,       &
    get_file_by_unit,   &
    init_file_loop,     &
    get_file_by_id
USE nlstcall_subs_mod, ONLY: totime,real_months
USE yomhook,  ONLY: lhook,dr_hook
USE parkind1, ONLY: jprb,jpim

USE profilename_length_mod, ONLY: fileid_length

USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: IOSTAT_END


IMPLICIT NONE 

! Read variables for nlstcall_nc
CHARACTER(LEN=filenamelength) :: filename
CHARACTER(LEN=fileid_length)  :: file_id
CHARACTER(LEN=filenamelength) :: filename_base

INTEGER :: reinit_unit      ! Re-init time unit
INTEGER :: reinit_end       ! Re-init end time
INTEGER :: reinit_step      ! Re-init frequency
INTEGER :: reinit_start     ! Re-init start time
INTEGER :: packing          ! Packing code
INTEGER :: nccomp_level     ! Deflate level, value between 1 and 9

LOGICAL :: l_compress       ! NetCDF 4 compression switch
LOGICAL :: l_shuffle        ! Switch to turn on netCDF4 shuffle filter

! Read variables for nlstcall_nc_options
INTEGER :: ncformat         ! NetCDF output format
INTEGER :: ncdim_prec       ! Floating point precision of dimensions
INTEGER :: ncvar_prec       ! Floating point precision of variables

LOGICAL :: l_netcdf         ! Main switch to turn on netCDF output
LOGICAL :: l_checksum       ! Switch to turn on netCDF4 fletcher32
                            ! checksum filter

! Used only when reading from HISTORY file
LOGICAL :: partial          ! Partially written

! DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

NAMELIST / nlstcall_nc / file_id,            &
                         reinit_end,         &
                         reinit_step,        &
                         reinit_start,       &
                         reinit_unit,        &
                         filename_base,      &
                         packing,            &
                         l_compress,         &
                         nccomp_level,       &
                         l_shuffle

NAMELIST / nlstcall_nc_options / l_netcdf,   &
                                 ncformat,   &
                                 ncdim_prec, &
                                 ncvar_prec, &
                                 l_checksum

NAMELIST / nlstcall_nc_hist / filename,      &
                              file_id

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='NLSTCALL_NC_NAMELIST_MOD'

CONTAINS

! This version of the routine is intended to read the history file
!-------------------------------------------------------------------------------
SUBROUTINE read_nml_nlstcall_nc_hist(unit_in)

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
INTEGER :: nc_unit
INTEGER :: outtime

CHARACTER(LEN=errormessagelength) :: iomessage

REAL(KIND=jprb)             :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_NLSTCALL_NC_HIST'
TYPE(um_file_type), POINTER :: um_file

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 2
INTEGER, PARAMETER :: n_int = 1
INTEGER, PARAMETER :: n_chars = filenamelength + fileid_length

TYPE my_namelist
  SEQUENCE
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
  filename   = ""
  file_id    = ""

  IF (mype == 0) THEN

    READ(UNIT=unit_in, NML=nlstcall_nc_hist, IOSTAT=ErrorStatus,      &
         IOMSG=iomessage)
    ! Don't check for negative error status (EOF indicator)
    IF (ErrorStatus /= IOSTAT_END) THEN
      CALL check_iostat(errorstatus, "namelist NLSTCALL_NC_HIST", iomessage)
    END IF

    my_nml % filename    = filename
    my_nml % file_id     = file_id
    my_nml % error       = ErrorStatus

  END IF

  CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

  IF (mype /= 0) THEN

    filename         = my_nml % filename
    file_id          = my_nml % file_id
    ErrorStatus      = my_nml % error

  END IF

  ! If this attempt to read the namelist hit the end of the file, stop
  ! (note ErrorStatus is broadcast WITH the nml variables to ensure all 
  ! PEs know when further reads are required)
  IF (ErrorStatus == IOSTAT_END) EXIT

  IF (printstatus >= prstatus_oper .AND. mype == 0) THEN
    CALL print_nlist_nlstcall_nc_hist()
  END IF

  ! Create a new file object for this stream
  NULLIFY(um_file)
  CALL assign_file_unit(filename, nc_unit, &
                        handler="netcdf", id=file_id)
  um_file => get_file_by_unit(nc_unit, handler="netcdf")

  ! Set the flag to indicate it is partially written
  um_file % meta % partially_written = .TRUE.

END DO

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_nlstcall_nc_hist

!-------------------------------------------------------------------------------
SUBROUTINE print_nlist_nlstcall_nc_hist()

USE umprintMgr, ONLY: umprint

IMPLICIT NONE

CHARACTER(LEN=50000) :: linebuffer
REAL(KIND=jprb)      :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_NLSTCALL_NC_HIST'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umprint('Contents of namelist nlstcall_nc_hist', src=modulename)

WRITE(linebuffer,"(A,A)")' filename = ', TRIM(filename)
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,A)")' file_id = ', TRIM(file_id)
CALL umprint(linebuffer,src=modulename)

CALL umprint('- - - - - - end of namelist - - - - - -', src=modulename)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_nlstcall_nc_hist

! This routine will read as many nlstcall_nc namelists as it finds in the
! file attached to the provided unit
!-------------------------------------------------------------------------------
SUBROUTINE read_nml_nlstcall_nc(unit_in)

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
INTEGER :: nc_unit
INTEGER :: outtime

CHARACTER(LEN=errormessagelength) :: iomessage

REAL(KIND=jprb)             :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_NLSTCALL_NC'
TYPE(um_file_type), POINTER :: um_file

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int = 7
INTEGER, PARAMETER :: n_log = 2
INTEGER, PARAMETER :: n_chars = filenamelength + fileid_length

TYPE my_namelist
  SEQUENCE
  INTEGER :: reinit_end
  INTEGER :: reinit_step
  INTEGER :: reinit_start
  INTEGER :: reinit_unit
  INTEGER :: packing
  INTEGER :: error
  INTEGER :: nccomp_level
  LOGICAL :: l_compress
  LOGICAL :: l_shuffle
  CHARACTER (LEN=filenamelength) :: filename_base
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
  packing          = 0
  filename_base    = ""
  file_id          = ""
  l_compress       = .FALSE.
  nccomp_level     = imdi
  l_shuffle        = .FALSE.

  IF (mype == 0) THEN

    READ(UNIT=unit_in, NML=nlstcall_nc, IOSTAT=ErrorStatus, IOMSG=iomessage)
    ! Don't check for negative error status (EOF indicator)
    IF (ErrorStatus /= IOSTAT_END) THEN
      CALL check_iostat(errorstatus, "namelist NLSTCALL_NC", iomessage)
    END IF

    my_nml % reinit_end       = reinit_end
    my_nml % reinit_step      = reinit_step
    my_nml % reinit_start     = reinit_start
    my_nml % reinit_unit      = reinit_unit
    my_nml % packing          = packing
    my_nml % filename_base    = filename_base
    my_nml % file_id          = file_id
    my_nml % l_compress       = l_compress
    my_nml % nccomp_level     = nccomp_level
    my_nml % l_shuffle        = l_shuffle
    my_nml % error            = ErrorStatus

  END IF

  CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

  IF (mype /= 0) THEN

    reinit_end       = my_nml % reinit_end
    reinit_step      = my_nml % reinit_step
    reinit_start     = my_nml % reinit_start
    reinit_unit      = my_nml % reinit_unit
    packing          = my_nml % packing
    filename_base    = my_nml % filename_base
    file_id          = my_nml % file_id
    l_compress       = my_nml % l_compress
    nccomp_level     = my_nml % nccomp_level
    l_shuffle        = my_nml % l_shuffle
    ErrorStatus      = my_nml % error

  END IF

  ! If this attempt to read the namelist hit the end of the file, stop
  ! (note ErrorStatus is broadcast WITH the nml variables to ensure all 
  ! PEs know when further reads are required)
  IF (ErrorStatus == IOSTAT_END) EXIT

  IF (printstatus >= prstatus_oper .AND. mype == 0) THEN
    CALL print_nlist_nlstcall_nc()
  END IF

  ! Otherwise all PEs should now initialise an object to hold the stream
  ! It's possible that an earlier read of the history namelist will have
  ! already initialised any streams which are going to continue as part
  ! of a CRUN, in these cases we will find the object already exists
  NULLIFY(um_file)
  um_file => get_file_by_id(file_id, handler="netcdf", &
                            ignore_missing=.TRUE.)
  IF (.NOT. ASSOCIATED(um_file)) THEN
    ! If the file object doesn't already exist the file will be new
    ! to this run, so initialise the object here instead.

    ! If the file is tagged for re-initialisation it won't be
    ! provided with a filename - so provide a placeholder for
    ! any logging messages in the meantime
    filename = "Reserved unit for re-initialised stream ("//&
                TRIM(filename_base)//")"
    ! Create a new file object for this stream
    CALL assign_file_unit(filename, nc_unit, &
                          handler="netcdf", id=file_id)
    um_file => get_file_by_unit(nc_unit, handler="netcdf")
  END IF

  ! Simple settings and flags from the namelist
  um_file % meta % packing_code = packing
  um_file % meta % is_output_file = .TRUE.

  ! Re-initialisation settings
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

  um_file % nc_meta % l_compress = l_compress
  um_file % nc_meta % nccomp_level = nccomp_level
  um_file % nc_meta % l_shuffle = l_shuffle

END DO

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_nlstcall_nc

!-------------------------------------------------------------------------------
SUBROUTINE print_nlist_nlstcall_nc()

USE umprintMgr, ONLY: umprint

IMPLICIT NONE

CHARACTER(LEN=50000) :: linebuffer
REAL(KIND=jprb)      :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_NLSTCALL_NC'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umprint('Contents of namelist nlstcall_nc', src=modulename)

WRITE(linebuffer,"(A,I0)")' reinit_end   = ', reinit_end
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,I0)")' reinit_step  = ', reinit_step
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,I0)")' reinit_start = ', reinit_start
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,I0)")' reinit_unit  = ', reinit_unit
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,I0)")' packing  = ', packing
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,A)") ' filename_base    = ', TRIM(filename_base)
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,A)") ' file_id  = ', TRIM(file_id)
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,L1)")' l_compress = ', l_compress
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,I0)")' nccomp_level = ', nccomp_level
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,L1)")' l_shuffle = ', l_shuffle
CALL umprint(linebuffer,src=modulename)

CALL umprint('- - - - - - end of namelist - - - - - -', src=modulename)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_nlstcall_nc

! This routine will read as many nlstcall_nc_options namelists as it finds
! in the file attached to the provided unit
!-------------------------------------------------------------------------------
SUBROUTINE read_nml_nlstcall_nc_options(unit_in)

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
INTEGER :: nc_unit
INTEGER :: outtime

CHARACTER(LEN=errormessagelength) :: iomessage

REAL(KIND=jprb)             :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_NLSTCALL_NC_OPTIONS'
TYPE(um_file_type), POINTER :: um_file

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 2
INTEGER, PARAMETER :: n_int = 4
INTEGER, PARAMETER :: n_log = 2

TYPE my_namelist
  SEQUENCE
  INTEGER :: ncformat
  INTEGER :: ncdim_prec
  INTEGER :: ncvar_prec
  INTEGER :: error
  LOGICAL :: l_netcdf
  LOGICAL :: l_checksum
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int, n_log_in=n_log)

! Setup defaults 
l_netcdf         = .FALSE.
ncformat         = imdi
ncdim_prec       = imdi
ncvar_prec       = imdi
l_checksum       = .FALSE.

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=nlstcall_nc_options, &
       IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist NLSTCALL_NC_OPTIONS", iomessage)

  my_nml % l_netcdf         = l_netcdf
  my_nml % ncformat         = ncformat
  my_nml % ncdim_prec       = ncdim_prec
  my_nml % ncvar_prec       = ncvar_prec
  my_nml % l_checksum       = l_checksum
  my_nml % error            = ErrorStatus

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  l_netcdf         = my_nml % l_netcdf
  ncformat         = my_nml % ncformat
  ncdim_prec       = my_nml % ncdim_prec
  ncvar_prec       = my_nml % ncvar_prec
  l_checksum       = my_nml % l_checksum
  ErrorStatus      = my_nml % error

END IF

IF (printstatus >= prstatus_oper .AND. mype == 0) THEN
  CALL print_nlist_nlstcall_nc_options()
END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_nlstcall_nc_options

!-------------------------------------------------------------------------------
SUBROUTINE print_nlist_nlstcall_nc_options()

USE umprintMgr, ONLY: umprint

IMPLICIT NONE

CHARACTER(LEN=50000) :: linebuffer
REAL(KIND=jprb)      :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_NLSTCALL_NC_OPTIONS'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umprint('Contents of namelist nlstcall_nc_options', src=modulename)

WRITE(linebuffer,"(A,L1)")' l_netcdf = ', l_netcdf
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,I0)")' ncformat = ', ncformat
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,I0)")' ncdim_prec = ', ncdim_prec
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,I0)")' ncvar_prec = ', ncvar_prec
CALL umprint(linebuffer,src=modulename)
WRITE(linebuffer,"(A,L1)")' l_checksum = ', l_checksum
CALL umprint(linebuffer,src=modulename)

CALL umprint('- - - - - - end of namelist - - - - - -', src=modulename)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_nlstcall_nc_options

! The inverse of the above routine - write current output streams to
! a file, for use in history
!-------------------------------------------------------------------------------
SUBROUTINE write_nml_nlstcall_nc_hist(unit_in)

USE um_parcore, ONLY: mype
USE ereport_mod, ONLY: ereport

IMPLICIT NONE 

INTEGER, INTENT(IN)  :: unit_in

INTEGER :: icode
REAL(KIND=jprb)             :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='WRITE_NML_NLSTCALL_NC_HIST'
TYPE(um_file_type), POINTER :: um_file

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! This should never be called by a PE other than PE0 but protect it 
! anyway just to be sure
IF (mype == 0) THEN
  
  ! Loop over all binary/non-portable files
  NULLIFY(um_file)
  um_file => init_file_loop(handler="netcdf")
  DO WHILE (ASSOCIATED(um_file))

    IF (um_file % meta % is_output_file) THEN
      ! Partially written files must be indicated in the namelist to
      ! allow a CRUN to pick up from the final written position
      IF (um_file % meta % partially_written) THEN
        filename    = um_file % filename
        file_id     = um_file % id

        ! Namelist is now populated and ready for output
        WRITE(unit_in, nlstcall_nc_hist, IOSTAT=icode)
        CALL ereport('write_nml_nlstcall_nc_hist',icode,                       &
                     'Write ERROR for NLSTCALL_NC_HIST namelist')
      END IF
    END IF
    ! Increment file pointer for next pass of loop
    um_file => um_file % next
  END DO
ELSE
  ! This shouldn't be called by other PEs, but warn instead of crashing since
  ! it isn't a critical mistake
  icode = -10
  CALL ereport('write_nml_nlstcall_nc_hist', icode,                            &
      'Should not be called by PEs other than PE0')
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE write_nml_nlstcall_nc_hist

END MODULE nlstcall_nc_namelist_mod

