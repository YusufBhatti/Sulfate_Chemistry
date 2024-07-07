! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Module to hold data and procedures for handling NetCDF emission files.
!
!  Method:
!    Contains a number of subroutines which call some Fortran NetCDF
!    procedures to get information/values of dimensions/variables
!    present in NetCDF emission files. It is generally assumed that dimensions
!    of variables are ordered x, y, z, t as recommended by CF conventions
!    (though any dimension can be omitted).
!    For further information see the header of each subroutine
!    contained here.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE emiss_io_mod

USE um_types,               ONLY: integer32, real32  ! For NetCDF arguments
USE netcdf
USE um_parcore,             ONLY: mype, nproc
USE field_types,            ONLY: fld_type_p
USE UM_Parvars,             ONLY: gc_all_proc_group
USE UM_ParParams,           ONLY: halo_type_no_halo
USE mpl,                    ONLY: mpl_character, mpl_integer, mpl_real
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE umPrintMgr,             ONLY: umPrint, umMessage, PrStatus_Diag, PrintStatus
USE conversions_mod,        ONLY: isec_per_day, rday_per_hour
USE setup_namelist,         ONLY: setup_nml_type
USE parkind1,               ONLY: jpim, jprb         ! DrHook
USE yomhook,                ONLY: lhook, dr_hook     ! DrHook
USE nlstcall_mod,           ONLY: lcal360

IMPLICIT NONE

! Default private
PRIVATE

! Local variables
!
INTEGER (KIND=integer32), SAVE, PUBLIC :: nfid     ! ID of NetCDF file
INTEGER                        :: ecode            ! return code for ereport
INTEGER, PARAMETER             :: max_dims = 5  ! Max nr of dimens per variable

! Expected max No. of values in any attribute, just in case some attributes
! may be arrays (could be changed to higher value if needed)
INTEGER, PARAMETER             :: max_att_len = 10

! Expected max. length of a variable name in a NetCDF file
! (this is equal to the length of 'var_name' within the
! data type 'ukca_em_struct' declared in UKCA_EMISS_MOD).
INTEGER, PARAMETER             :: varname_len = 80

! Expected max. length of the units in the time variable
INTEGER, PARAMETER             :: t_units_len = 40

! Lenght for the strings indicating variable type (e.g. 'CHAR', 'REAL', etc.)
INTEGER, PARAMETER             :: vartype_len = 4

! Declare integer parameters indicating the possible ways of updating
! emissions, following the same conventions as for ancillary files
INTEGER, PARAMETER, PUBLIC :: iupd_single = 0  ! single time (not tested yet)
INTEGER, PARAMETER, PUBLIC :: iupd_serial = 1  ! time series
INTEGER, PARAMETER, PUBLIC :: iupd_cyclic = 2  ! perodic time series

! Arbitrary year for cyclic files
INTEGER, PARAMETER, PUBLIC :: cyclic_year = 1

! Variables to hold information on the record ('time') dimension.
! These will be populated/allocated (unless not requested) during em_fopen 
! and reset/deallocated (if populated) with the em_fclose call
CHARACTER(LEN=varname_len), SAVE :: tim_name = 'empty'
CHARACTER(LEN=t_units_len), SAVE :: ct_units = 'empty'
INTEGER,                    SAVE :: varid_t  = -99
INTEGER,                    SAVE :: dimid_t  = -99
LOGICAL,                    SAVE :: lcal360file = .FALSE.
REAL, ALLOCATABLE                :: time_axis(:)
 
! Subroutines and overloaded interfaces available outside this module
PUBLIC :: em_fopen, em_fclose, nd_error
PUBLIC :: em_var_check_dims, em_get_time_rec
PUBLIC :: em_get_file_info, em_get_var_info
PUBLIC :: em_get_data
PUBLIC :: em_get_var_att

! ------------------------------------------------------------------------
! Interface blocks used for overloading. 
!
! A generic call to EM_GET_DATA will select the appropriate subroutine 
! automatically, depending on the type of the INTENT(INOUT) argument
! (variable values).
!
! A generic call to EM_GET_VAR_ATT will select the appropriate subroutine
! automatically, depending on the type of the INTENT(OUT) argument 
! (attribute value).
!
! Note that both EM_GET_DATA and EM_GET_VAR_ATT are public. Therefore the
! generic interfaces can be called from outside this module, while the
! type specific subroutines can only be used inside this module.
! ------------------------------------------------------------------------

INTERFACE em_get_data
MODULE PROCEDURE            &
  em_get_data_char,         &
  em_get_data_int,          &
  em_get_data_real,         &
  em_get_data_real2D,       &
  em_get_data_real2D_zonal, &
  em_get_data_real3D,       &
  em_get_data_real3d_zonal, &
  em_get_data_real4D
END INTERFACE

INTERFACE em_get_var_att
MODULE PROCEDURE            &
  em_get_var_att_char,      &
  em_get_var_att_int,       &
  em_get_var_att_real
END INTERFACE

INTERFACE em_var_check_dims
MODULE PROCEDURE                  &
  em_var_check_dims_2d_zonal,     &
  em_var_check_dims_3d_zonal,     &
  em_var_check_dims_3d,           &
  em_var_check_dims_4d
END INTERFACE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EMISS_IO_MOD'

CONTAINS

! ------------------------------------------------------------------------
! Description:
!   Check for existence of NetCDF emission file and open it.
!
! Method:
!   Call NF90_OPEN to open a NetCDF file and return nstatus.
!   Call ND_ERROR, which will stop the model if nstatus reports errors.
!   By default, check that a record ('time') dimension exists in the 
!   file, through a call to EM_GET_TIME_INFO. This can be overriden by
!   optional argument ignore_time = True
! ------------------------------------------------------------------------
SUBROUTINE em_fopen (filename, fileid, ignore_time)

IMPLICIT NONE

! Subroutine arguments
CHARACTER (LEN=*), INTENT(IN)  :: filename

INTEGER,           INTENT(OUT) :: fileid

! Optional flag that will be set to True only if it is needed to prevent
! the default attempt to detect the time dimension in a file.
LOGICAL, OPTIONAL, INTENT(IN)  :: ignore_time

! Local variables
INTEGER (KIND=integer32)       :: nstatus       ! netCDF return code

! Local version of "ignore_time" 
LOGICAL :: l_ignore_time 

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EM_FOPEN'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF ( mype == 0 ) THEN
  nstatus = nf90_open (filename, nf90_nowrite, nfid)
  CALL nd_error (nstatus, 'EM_FOPEN', TRIM(filename))
  fileid = nfid
ELSE
  ! Set file-id to dummy on other PEs, as they
  ! would not access the files anyway
  fileid = 99  
END IF

! Copy the value of argument "ignore_time" to the local variable
! else - if not specified - set to False so that information on 
!  'time' variable/dimension is searched by default in all files
IF (PRESENT(ignore_time)) THEN
  l_ignore_time = ignore_time
ELSE
  l_ignore_time = .FALSE.
END IF

! By default EM_GET_TIME_INFO is called when opening a file to look
! for the time dimension (1-D, 'days since' or 'hours since') and
! stores some values in global variables. However make sure that
! l_ignore_time is not set to True just in case there are situations
! when the call is not needed.
IF ( .NOT. l_ignore_time ) THEN
  CALL em_get_time_info (fileid)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE em_fopen

! --------------------------------------------------------------------------
! Description:
!   Obtain the following information on a NetCDF file whose id has been
!   provided: number of variables, number of dimensions, number of (global)
!   attributes and record/time dimension
!
! Method:
!    Call NF90_INQUIRE and return only the information required. The data
!    is read on PEO and broadcast to all PEs. 
! --------------------------------------------------------------------------
SUBROUTINE em_get_file_info (fileid,  num_vars, num_dims, num_att, unlim_id)

IMPLICIT NONE

! Subroutine arguments

INTEGER,           INTENT(IN)    :: fileid          ! File id

INTEGER, OPTIONAL, INTENT(OUT)   :: num_vars        ! Nr of variables
INTEGER, OPTIONAL, INTENT(OUT)   :: num_dims        ! Nr of dimensions
INTEGER, OPTIONAL, INTENT(OUT)   :: num_att         ! Nr of Attributes
INTEGER, OPTIONAL, INTENT(OUT)   :: unlim_id        ! ID of  the record
                                                    ! or 'time' dimension

! Local variables
INTEGER (KIND=integer32) :: nvars32, ndims32, natts32, undimid32
INTEGER (KIND=integer32) :: nstatus ! netCDF return code

! For Parallel I/O
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: icode

INTEGER, PARAMETER :: no_of_types = 1 ! integer
INTEGER, PARAMETER :: n_int = 4       ! num_vars & num_dims & num_att & unlim_id

! Derived type used to hold information about a NetCDF file
! (read on PEO and broadcast to all PEs).
TYPE nc_fileinfo
  SEQUENCE
  INTEGER :: num_vars
  INTEGER :: num_dims
  INTEGER :: num_att
  INTEGER :: unlim_id
END TYPE nc_fileinfo

TYPE (nc_fileinfo) :: file_info

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EM_GET_FILE_INFO'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Obtain the communication group for this run
CALL gc_get_communicator (my_comm, icode)

! Obtain an object representing the multi-datatype structure which is
! going to be broadcast. This is necessary as the broadcast operation
! needs to know exactly how many bytes of data need to be transferred.
! SETUP_NML_TYPE uses the number of variables and their size to calculate
! the resulting number of bytes as an MPL datatype.
CALL setup_nml_type (no_of_types, mpl_nml_type, n_int_in = n_int)

! PE0 reads data from the file and adds it to the derived type file_info
IF ( mype == 0 ) THEN
  nfid = fileid

  ! Get all file information anyway
  nstatus = nf90_inquire ( nfid, ndims32, nvars32, natts32, undimid32 ) 
  CALL nd_error (nstatus, 'EM_GET_FILE_INFO','INQUIRE')

  file_info % num_vars = nvars32
  file_info % num_dims = ndims32
  file_info % num_att  = natts32
  file_info % unlim_id = undimid32
END IF

! Data structure broadcast to all processors
icode = 0
CALL mpl_bcast (file_info, 1, mpl_nml_type, 0, my_comm, icode)
IF ( icode /= 0 ) THEN
  ecode = fileid
  CALL ereport ('EM_GET_FILE_INFO', ecode, 'NML Broadcast Error')
END IF

! Local variables have been filled on PE0. Now copy
! data structure (from the derived type file_info)
! to the local variables for the rest of PEs.
IF ( mype /= 0 ) THEN
  nvars32   = file_info % num_vars
  ndims32   = file_info % num_dims
  natts32   = file_info % num_att
  undimid32 = file_info % unlim_id
END IF

! Check what information is required and transfer data from
! the local variables to the output arguments. Operation
! done on all PEs.
IF ( PRESENT(num_vars) ) num_vars = nvars32
IF ( PRESENT(num_dims) ) num_dims = ndims32
IF ( PRESENT(num_att ) ) num_att  = natts32
IF ( PRESENT(unlim_id) ) unlim_id = undimid32

! Reset the multi-datatype structure allocated for broadcasting
CALL mpl_type_free (mpl_nml_type, icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE em_get_file_info

! --------------------------------------------------------------------------
! Description:
!   Obtain the following information on a variable (in a NetCDF file)
!   whose name or id has been provided: variable id or name, number of
!   dimensions and corresponding sizes, whether there is record/
!   unlimited dimension, and variable type.
!
! Method:
! 1) Check whether name or ID of variable in NetCDF file is supplied
!    (if name is supplied then get the id). Store id as varid_local and
!    pass it to NF90_INQUIRE_VARIABLE to get info on the variable.
! 2) Loop through all dimensions in the variable and call
!    NF90_INQUIRE_DIMENSION for each one to get their size.
! 3) If the optional argument lrec is present then call NF90_INQUIRE
!    and get the time record dimension, assuming it is the last dims
!    in the file.
! 4) If the optional argument var_type is present then get it as string.
! Note that these operations are done by PE0 and the relevant information
! is broadcast to all PEs if needed.
! --------------------------------------------------------------------------
SUBROUTINE em_get_var_info (ncid,    varid, varnam,                   &
                            ndims, vdimsize, lrec, var_type)

IMPLICIT NONE

! Subroutine arguments

INTEGER,           INTENT(IN)    :: ncid                ! File id

INTEGER,           INTENT(INOUT) :: varid               ! Variable id

CHARACTER (LEN=*), INTENT(INOUT) :: varnam              ! Variable name

INTEGER, OPTIONAL, INTENT(OUT)   :: ndims               ! Nr of dims
INTEGER, OPTIONAL, INTENT(OUT)   :: vdimsize(max_dims)  ! Size of dims

! Flag set to True if variable contains record/time dimension
LOGICAL, OPTIONAL, INTENT(OUT)   :: lrec

CHARACTER(LEN=vartype_len), OPTIONAL, INTENT(OUT) :: var_type ! Variable type


! Local variables
INTEGER (KIND=integer32)    :: nstatus        ! netCDF return code
INTEGER (KIND=integer32)    :: varid_local, itype, idumm, n, undid
INTEGER (KIND=integer32)    :: ndim_local
INTEGER (KIND=integer32)    :: dimid       (max_dims)
INTEGER (KIND=integer32)    :: dsize_local (max_dims)

LOGICAL                     :: lrec_local
CHARACTER (LEN=vartype_len) :: var_type_local
CHARACTER (LEN=varname_len) :: varnam_local

! For Parallel I/O
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: icode

INTEGER, PARAMETER :: no_of_types = 3           ! int + logical + char
INTEGER, PARAMETER :: n_int   = 2 + max_dims    ! varid & ndims & max_dims
INTEGER, PARAMETER :: n_log   = 1               ! lrec

! Nr. of characters for varname and var_type (see below)
INTEGER, PARAMETER :: n_chars = varname_len + vartype_len

! Derived type used to hold information about a variable in a NetCDF file
! (read on PEO and broadcast to all PEs).
TYPE nc_fileinfo
  SEQUENCE
  INTEGER                    :: varid
  INTEGER                    :: ndims
  INTEGER                    :: vdimsize(max_dims)
  LOGICAL                    :: lrec
  CHARACTER(LEN=varname_len) :: varname
  CHARACTER(LEN=vartype_len) :: var_type
END TYPE nc_fileinfo

TYPE (nc_fileinfo) :: file_info

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EM_GET_VAR_INFO'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Obtain the communication group for this run
CALL gc_get_communicator(my_comm, icode)

! Obtain an object representing the multi-datatype structure which is
! going to be broadcast (see further details above)
CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in = n_int,    &
                    n_log_in = n_log, n_chars_in = n_chars )

! Check whether Name or ID of variable in NetCDF file is supplied.
! PE0 reads the relevant data and adds it to the derived type file_info.
IF ( mype == 0 ) THEN
  nfid = ncid

  ! The following lines include NetCDF calls to populate local variables
  ! (e.g. varid_local, ndim_local, dsize_local(:), varnam_local,
  ! var_type_local) which will be transferred to the broadcast type.
  ! There is no need to initialise most of those local variables here,
  ! because if they were not correctly populated the NetCDF calls would
  ! report errors anyway.

  ! If varnam provided then get varid
  IF ( TRIM(varnam) /= "" ) THEN
    ! Note: 'varnam' is intent(in) in nf90_inq_varid.
    ! See note below about 'varnam_local'.
    nstatus = nf90_inq_varid (nfid, TRIM(varnam), varid_local)
    CALL nd_error ( nstatus, 'EM_GET_VAR_INFO', 'INQ_VARID ' // varnam )
    varid = varid_local

  ! If varnam not provided then make sure that varid has
  ! been provided; if that is not the case then report error.
  ELSE
    IF ( varid <= 0 ) THEN
      ecode = 1
      CALL ereport ('EM_GET_VAR_INFO', ecode,           &
        'Either variable name or ID should be provided')
    END IF
    varid_local = varid

  END IF            ! check if name/id provided

  ! Get type information, dimension ids, etc.
  ! 
  ! Here we get the name of the variable, now 'varnam_local',
  ! which unlike 'varnam' above needs to be intent(out) in
  ! nf90_inquire_variable. That is the reason why two
  ! 'varnams' are used. It is necessary to get 'varnam_local'
  ! because this routine can be called with 'varnam' empty.
  !
  ! Note also that if the user did pass in 'varnam' it will be
  ! overwritten by 'varnam_local', but that is OK because 
  ! the code above has ensured that the variable name 
  ! corresponds to the id being inquired about.
  !
  nstatus = nf90_inquire_variable (nfid,  varid_local, varnam_local, itype, &
                                   ndim_local, dimid, idumm)

  CALL nd_error (nstatus, 'EM_GET_VAR_INFO',                                & 
                'INQ_VARIABLE ' // TRIM(varnam_local) )
  varnam = varnam_local

  ! Get size of each dimension
  IF ( PRESENT(vdimsize) ) THEN
    DO n = 1, ndim_local
      nstatus = nf90_inquire_dimension (nfid, dimid(n), LEN = dsize_local(n))
      CALL nd_error (nstatus, 'EM_GET_VAR_INFO', 'INQ_DIMSIZE ' // TRIM(varnam))
    END DO
    vdimsize (1:ndim_local) = dsize_local (1:ndim_local)
  END IF

  ! Check if the variable contains time dimension (mandatory for files). 
  ! The dimension ID for time should have been stored during file-open.
  ! Compare this ID with that of the last dimension in the variable,
  ! as we assume that time record is always the last dim.
  lrec_local = .FALSE.
  IF ( PRESENT(lrec) ) THEN
    IF ( dimid_t < 0 ) THEN
      WRITE (umMessage,'(A,A)') 'EM_GET_VAR_INFO: ' //                   &
       'Time dimension info not found, check if em_fopen was called ' // &
       'with option l_ignore_time = .TRUE.'
      CALL umPrint(umMessage,src='em_get_var_info')
      nstatus = nf90_enotvar
      CALL nd_error (nstatus,  'EM_GET_VAR_INFO',                        &
                    'No information on time dimension for this file')      
    ELSE
      IF (dimid(ndim_local) == dimid_t ) lrec_local = .TRUE.
    END IF
  END IF

  ! Convert type into string (only the following common types for now)
  IF ( PRESENT(var_type) ) THEN
    SELECT CASE(itype)
    CASE (nf90_char)
      var_type_local = 'CHAR'
    CASE (nf90_int)
      var_type_local = 'INT'
    CASE (nf90_float)
      var_type_local = 'REAL'
    CASE (nf90_double)
      var_type_local = 'DBLE'
    CASE DEFAULT
      WRITE (umMessage,'(A,1X,I2)') 'EM_GET_VAR_INFO: ' //               &
                       'Unknown variable type defined ', itype
      CALL umPrint (umMessage, src='em_get_var_info')
      WRITE (umMessage,'(A)')  'Should be one of: CHAR / INT / REAL / DBLE'
      CALL umPrint (umMessage, src='em_get_var_info')
      nstatus = nf90_ebadname
      CALL nd_error (nstatus, 'EM_GET_VAR_INFO', 'VAR_TYPE ' // TRIM(varnam))
    END SELECT
  END IF 

  ! Fill in derived type variable to later do broadcast
  file_info % varid       = varid_local
  file_info % ndims       = ndim_local
  file_info % vdimsize(:) = dsize_local(:)
  file_info % varname     = varnam_local
  file_info % lrec        = lrec_local
  file_info % var_type    = var_type_local

END IF

! Data structure broadcast to all processors
icode = 0
CALL mpl_bcast (file_info, 1, mpl_nml_type, 0, my_comm, icode)
IF ( icode /= 0 ) THEN
  ecode = ncid
  CALL ereport ('EM_GET_VAR_INFO', ecode, 'NML Broadcast Error')
END IF

! Local variables have been filled on PE0. Now transfer
! data from the broadcast type into the local variables
! for the rest of PEs.
IF ( mype /= 0 ) THEN
  varid_local    = file_info % varid
  varnam_local   = file_info % varname
  ndim_local     = file_info % ndims
  dsize_local(:) = file_info % vdimsize(:)
  lrec_local     = file_info % lrec
  var_type_local = file_info % var_type
END IF

! Copy the local variables back onto the output variables
varid  = varid_local
varnam = varnam_local
IF ( PRESENT(ndims) )    ndims       = ndim_local
IF ( PRESENT(vdimsize) ) vdimsize(:) = dsize_local(:)
IF ( PRESENT(lrec) )     lrec        = lrec_local
IF ( PRESENT(var_type) ) var_type    = var_type_local

! Reset the multi-datatype structure allocated for broadcasting
CALL mpl_type_free (mpl_nml_type, icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE em_get_var_info

! --------------------------------------------------------------------------
! Description:
!   Obtain the Variable name, Units, Variable ID and Dimension ID for the 
!   Record or 'time' dimension in a NetCDF file.
!   Assumes the 'time' variable is uni-dimensional and has an 'units' attribute
!   starting with 'days since' or 'hours since'
!
! Method:
! 1) Loop through all variables, identify the uni-dimensional ones and search 
!    for above strings in 'units'. This will give the variable name and id.
! 2) Obtain the dimension id of that variable and verify that it exists as the
!    last dimension of the first 4-D variable encountered.
! 3) The name and ids for 'time' are stored in variables global to this module.
! 4) Uses actual netcdf calls and not wrappers from this module as most of the 
!    work is done by PE0 only and there is no need to broadcast any
!    intermediate data.
! --------------------------------------------------------------------------
SUBROUTINE em_get_time_info (ncid)

IMPLICIT NONE

! Subroutine arguments

INTEGER,           INTENT(IN)    :: ncid                ! File id

! Local variables
INTEGER (KIND=integer32) :: nstatus      ! netCDF return code
INTEGER (KIND=integer32) :: ndim_local
INTEGER (KIND=integer32) :: dimid       (max_dims)
INTEGER (KIND=integer32) :: dsize_local (max_dims)
INTEGER (KIND=integer32) :: itype
INTEGER (KIND=integer32) :: varid32
INTEGER (KIND=integer32) :: nvars32, ndims32, natts32, undimid32, icalflex32

! Use nf90_max_name to allocate the maximum possible length for these two
! local variables, but they will be truncated anyway when copied back to
! the global varibles tim_name and ct_units.
CHARACTER(LEN=nf90_max_name) :: varnam_local  ! Name read from file
CHARACTER(LEN=nf90_max_name) :: cunits        ! For 'units' attribute
CHARACTER(LEN=nf90_max_name) :: ccalendar     ! For 'calendar' attribute

INTEGER :: v               ! loop counter for variables
LOGICAL :: l_verify_time   ! Additional verification

! For Parallel I/O
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: icode

INTEGER, PARAMETER :: no_of_types = 3                     ! int + char + logical
INTEGER, PARAMETER :: n_int       = 2                     ! varid_t & dimid_t
INTEGER, PARAMETER :: n_chars = varname_len + t_units_len ! varname_t & cunits_t
INTEGER, PARAMETER :: n_log       = 1                     ! lcal360file

! Derived type used to hold time information from a NetCDF file
! (read on PEO and broadcast to all PEs).
TYPE nc_fileinfo
  SEQUENCE
  INTEGER                    :: varid_t
  INTEGER                    :: dimid_t
  CHARACTER(LEN=varname_len) :: varname_t
  CHARACTER(LEN=t_units_len) :: cunits_t
  LOGICAL                    :: lcal360file
END TYPE nc_fileinfo

TYPE (nc_fileinfo) :: file_info

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EM_GET_TIME_INFO'
CHARACTER(LEN=errormessagelength) :: cmessage  ! error message

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Obtain the communication group for this run
CALL gc_get_communicator (my_comm, icode)

! Obtain an object representing the multi-datatype structure which is
! going to be broadcast (see further details above) 
CALL setup_nml_type (no_of_types, mpl_nml_type, n_int_in = n_int,    &
                     n_chars_in = n_chars, n_log_in = n_log )

! File access by PE0, which reads the data needed and adds it
! to the derived type file_info
IF ( mype == 0 ) THEN
  nfid = ncid

  ! Get information on file 
  nstatus = nf90_inquire ( nfid, ndims32, nvars32, natts32, undimid32 ) 
  CALL nd_error (nstatus, 'EM_GET_TIME_INFO','INQUIRE_FILE')

  ! First loop over variables to search for any unidimensional ones. If
  ! found then search for the time dimension (by checking attribute units).
  DO v = 1, nvars32
    
    varid32 = v
    nstatus = nf90_inquire_variable (nfid,  varid32, varnam_local, itype, &
                                     ndim_local, dimid)
    CALL nd_error ( nstatus, 'EM_GET_TIME_INFO', 'INQ_VARIABLE ' )
    
    IF ( ndim_local == 1 ) THEN 
      
      ! Get 'units' attribute if present
      cunits = ' '
      nstatus = nf90_get_att (nfid, varid32, 'units', cunits) 
      IF ( nstatus == nf90_noerr) THEN          
       ! Attribute exists - check for string
        IF ( cunits(1:10) == 'days since'     .OR.     &
             cunits(1:11) == 'hours since' )  THEN  
          ! Found time. Populate global variables holding information
          ! on time dimension and  exit loop.
          varid_t   = varid32
          dimid_t   = dimid(1)
          tim_name  = TRIM(varnam_local)
          ct_units  = TRIM(cunits)
          EXIT   ! Exit because found
        ELSE
          CYCLE  ! Go to next variable
        END IF
      ELSE
        CYCLE    ! Attribute not found or error
      END IF

    END IF      ! 1-dimensional ?

  END DO

  ! If matching variable not found - abort
  IF ( TRIM(tim_name) == 'empty' .OR. varid_t == -99 )   THEN
    nstatus = nf90_enotvar   ! Status = Variable not found
    CALL nd_error (nstatus, 'EM_GET_TIME_INFO','Time variable not found')
  END IF

  ! Time coordinate variable found: now check 'calendar' attribute
  ccalendar = ' '
  varid32 = varid_t
  nstatus = nf90_get_att (nfid, varid32, 'calendar', ccalendar) 
  IF ( nstatus == nf90_noerr) THEN
    ! Set file calendar logical
    IF (TRIM(ccalendar) == '360_day' .OR. &
        TRIM(ccalendar) == '360') THEN
      ! 360-day model calendar. Older versions of CF allowed '360' instead of
      ! '360_day', so check for both.
      lcal360file = .TRUE.
    ELSE IF (TRIM(ccalendar) == 'gregorian' .OR. &
             TRIM(ccalendar) == 'standard') THEN
      lcal360file = .FALSE.
    ELSE
      WRITE (cmessage,'(A,A,A,A)')                               &
        'Time coordinate: ', TRIM(tim_name), &
        ' has unrecognised calendar: ', TRIM(ccalendar)
      ecode = nfid
      CALL ereport (RoutineName, ecode, cmessage)
    END IF
  ELSE
    ! Attribute not found - all time coordinates must have a calendar att
    WRITE (cmessage,'(A,A)')                                    &
      'No calendar attribute for time coordinate variable: ', TRIM(tim_name)
    CALL nd_error (nstatus, 'EM_GET_TIME_INFO', cmessage)
  END IF

  ! Also get calendar_flexible attribute, and use it to set icalflex.
  ! NetCDF3 doesn't have boolean type, so using integer = 0/1 instead
  nstatus = nf90_get_att (nfid, varid32, 'calendar_flexible', icalflex32)
  IF (nstatus /= nf90_noerr) THEN
    ! Attribute doesn't exist, default is non-flexible
    icalflex32 = 0
  END IF

  ! Check compatability of calendar
  IF (lcal360 .NEQV. lcal360file) THEN
    ! File calendar doesn't match model calendar. Raise an error unless time
    ! coordinate has attribute calendar_flexible=1.
    IF (icalflex32 /= 1) THEN
      WRITE (cmessage,'(A,A,A)')                               &
          'NetCDF calendar ', TRIM(ccalendar) , ' does not match model' // &
          ' calendar, and attribute calendar_flexible /= 1 or not set'
      ecode = nfid
      CALL ereport (RoutineName, ecode, cmessage)
    END IF
  END IF

  ! Second loop over variables to check that the time dimension is
  ! the last dimension in 4D or 5D variables.
  l_verify_time = .FALSE.
  DO v = 1, nvars32
    
    varid32 = v
    nstatus = nf90_inquire_variable (nfid,  varid32, varnam_local, itype, &
                                     ndim_local, dimid)
    CALL nd_error ( nstatus, 'EM_GET_TIME_INFO', 'INQ_VARIABLE 4D' )
    
    IF ( (ndim_local == 4 .AND. dimid(4) == dimid_t) .OR. &
         (ndim_local == 5 .AND. dimid(5) == dimid_t) ) THEN      
      l_verify_time = .TRUE.
      EXIT   ! Exit because found
    ELSE
      CYCLE  ! Go to next variable
    END IF

  END DO

  ! Time dimension does not seem to exist in any 4-D variable
  IF ( .NOT. l_verify_time ) THEN
    nstatus = nf90_enotvar   ! Status = Variable not found
    WRITE (umMessage,'(A,A)')                                    &
     'Time dimension exists, but missing from data variables: ', &
      TRIM(tim_name)
    CALL nd_error (nstatus, 'EM_GET_TIME_INFO', umMessage)
  END IF

  ! Copy variables into broadcast type
  file_info % varid_t   = varid_t
  file_info % dimid_t   = dimid_t
  file_info % varname_t = tim_name
  file_info % cunits_t  = ct_units
  file_info % lcal360file = lcal360file

END IF   ! PE0 work

! Sync PEs and broadcast time dimen-related values
CALL gc_gsync(nproc,icode)
icode = 0

! Data structure broadcast to all processors
CALL mpl_bcast (file_info, 1, mpl_nml_type, 0, my_comm, icode)
IF ( icode /= 0 ) THEN
  ecode = ncid
  CALL ereport ('EM_GET_TIME_INFO', ecode, 'Broadcast Error')
END IF

! Global variables holding time information have been filled
! on PE0. Now copy data structure to those variables
! for the rest of PEs.
IF ( mype /= 0 ) THEN
  varid_t  = file_info % varid_t
  dimid_t  = file_info % dimid_t
  tim_name = file_info % varname_t
  ct_units = file_info % cunits_t
  lcal360file = file_info % lcal360file
END IF

! Reset the multi-datatype structure allocated for broadcasting
CALL mpl_type_free (mpl_nml_type, icode)

! Finally read the time axis and store in a global array
CALL em_get_time_axis (ncid)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE em_get_time_info

! --------------------------------------------------------------------------
! Description:
!   Check that dimensions of a variable in a NetCDF file match global sizes,
!   in the case of a two-dimensional, zonally-averaged distribution.
!
! Method:
! 1) Call EM_GET_VAR_INFO to get the dimensions and check they match.
! 2) If wrong dimensions or if there is no time record then stop with
!    call to EREPORT
! --------------------------------------------------------------------------
SUBROUTINE em_var_check_dims_2d_zonal (global_rows, model_levels, &
                           l_zonal_mean, fileid, varname, dimsize)

IMPLICIT NONE

! Subroutine arguments
INTEGER,           INTENT(IN)    :: global_rows
INTEGER,           INTENT(IN)    :: model_levels
LOGICAL,           INTENT(IN)    :: l_zonal_mean
INTEGER,           INTENT(IN)    :: fileid

CHARACTER (LEN=*), INTENT(INOUT) :: varname

INTEGER,           INTENT(OUT)   :: dimsize(max_dims)

! Local variables
INTEGER                        :: vid    ! variable id
LOGICAL                        :: lreco  ! true if contains record dim

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EM_VAR_CHECK_DIMS_2D_ZONAL'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! First just check if this routine should have been called
IF (.NOT. l_zonal_mean) THEN
  WRITE (umMessage,'(A,A,A)') 'em_var_check_dims_2d_zonal called for ' // &
                             TRIM(varname) // &
                            ' but with l_zonal_mean as .false. '
  CALL umPrint(umMessage,src='em_var_check_dims_2d_zonal')
  ecode = fileid
  CALL ereport ('EM_VAR_CHECK_DIMS_2D_ZONAL', ecode,                      &
      'l_zonal_mean was .false. so em_var_check_dims_2d_zonal: ' //       &
      ' should not have been called for ' // TRIM(varname))
END IF

! Emission files should be defined as (y,lev,t)
!
! Get information on variable. One of the two INOUT arguments
! of EM_GET_VAR_INFO is supplied (varname); assign an arbitrary
! negative value to the other one since it is not known yet (vid).
vid = -1
CALL em_get_var_info (fileid, vid, varname, vdimsize=dimsize, lrec=lreco)

! Check that dimensions (y, lev) match those of the model
IF ( dimsize(1) /= global_rows .OR.  &
    (dimsize(2) > 1 .AND. dimsize(2) /= model_levels)               .OR.  &
     dimsize(2) < 1 ) THEN 

  WRITE (umMessage,'(A,A)') 'ERROR: Variable dimensions do not match for ',&
                   TRIM(varname)
  CALL umPrint(umMessage,src='em_var_check_dims')
  WRITE (umMessage,'(A,2(1x,I0))') 'Expected:',  &
                    global_rows, model_levels
  CALL umPrint(umMessage,src='em_var_check_dims')
  WRITE (umMessage,'(A,2(1x,I0))') 'Found:', dimsize(1:2)
  CALL umPrint(umMessage,src='em_var_check_dims')
  ecode = fileid
  CALL ereport ('EM_VAR_CHECK_DIMS', ecode,                               &
      'Dimension mismatch for: ' // TRIM(varname))
END IF

! Check that there is time record
IF (.NOT. lreco ) THEN 
  WRITE (umMessage,'(A,A)') 'No record dimension found for ', TRIM(varname)
  CALL umPrint(umMessage,src='em_var_check_dims_2d_zonal')
  ecode = fileid
  CALL ereport ('EM_VAR_CHECK_DIMS_2D_ZONAL', ecode,                       &
      'No record dimension found for ' // TRIM(varname))
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE em_var_check_dims_2d_zonal

! --------------------------------------------------------------------------
! Description:
!   Check that dimensions of a variable in a NetCDF file match global sizes,
!   in the case of a full-atmosphere three-dimensional distribution.
!
! Method:
! 1) Call EM_GET_VAR_INFO to get the dimensions and check they match.
! 2) If wrong dimensions or if there is no time record then stop with
!    call to EREPORT
! --------------------------------------------------------------------------
SUBROUTINE em_var_check_dims_3d (global_row_length, global_rows, model_levels, &
                                 fileid, varname, dimsize)

IMPLICIT NONE

! Subroutine arguments
INTEGER,           INTENT(IN)    :: global_row_length
INTEGER,           INTENT(IN)    :: global_rows
INTEGER,           INTENT(IN)    :: model_levels
INTEGER,           INTENT(IN)    :: fileid

CHARACTER (LEN=*), INTENT(INOUT) :: varname

INTEGER,           INTENT(OUT)   :: dimsize(max_dims)

! Local variables
INTEGER                        :: vid    ! variable id
LOGICAL                        :: lreco  ! true if contains record dim

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EM_VAR_CHECK_DIMS_3D'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Emission files should be defined as (x,y,lev,t)
!
! Get information on variable. One of the two INOUT arguments
! of EM_GET_VAR_INFO is supplied (varname); assign an arbitrary
! negative value to the other one since it is not known yet (vid).
vid = -1
CALL em_get_var_info (fileid, vid, varname, vdimsize=dimsize, lrec=lreco)

! Check that dimensions (x, y, lev) match those of the model
IF ( dimsize(1) /= global_row_length .OR. dimsize(2) /= global_rows .OR.  &
    (dimsize(3) > 1 .AND. dimsize(3) /= model_levels)               .OR.  &
     dimsize(3) < 1 ) THEN  

  WRITE (umMessage,'(A,A)') 'ERROR: Variable dimensions do not match for ',&
                   TRIM(varname)
  CALL umPrint(umMessage,src='em_var_check_dims')
  WRITE (umMessage,'(A,3(1x,I5))') 'Expected:', global_row_length,  &
                    global_rows, model_levels
  CALL umPrint(umMessage,src='em_var_check_dims')
  WRITE (umMessage,'(A,3(1x,I5))') 'Found:', dimsize(1:3)
  CALL umPrint(umMessage,src='em_var_check_dims')
  ecode = fileid
  CALL ereport ('EM_VAR_CHECK_DIMS_3D', ecode,                            &
      'Dimension mismatch for: ' // TRIM(varname))
END IF

! Check that there is time record
IF (.NOT. lreco ) THEN  
  WRITE (umMessage,'(A,A)') 'No record dimension found for ', TRIM(varname)
  CALL umPrint(umMessage,src='em_var_check_dims')
  ecode = fileid
  CALL ereport ('EM_VAR_CHECK_DIMS_3D', ecode,                            &
      'Dimension mismatch for: ' // TRIM(varname))
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE em_var_check_dims_3d

! --------------------------------------------------------------------------
! Description:
!   Check that dimensions of a variable in a NetCDF file match global sizes,
!   in the case of a time-dependent, zonally-averaged distribution.
!
! Method:
! 1) Call EM_GET_VAR_INFO to get the dimensions and check they match.
! 2) If wrong dimensions or if there is no time record then stop with
!    call to EREPORT
! --------------------------------------------------------------------------
SUBROUTINE em_var_check_dims_3d_zonal (global_rows, model_levels, &
                    dimension_3, l_zonal_mean, fileid, varname, dimsize)

IMPLICIT NONE

! Subroutine arguments
INTEGER,           INTENT(IN)    :: global_rows
INTEGER,           INTENT(IN)    :: model_levels
INTEGER,           INTENT(IN)    :: dimension_3
LOGICAL,           INTENT(IN)    :: l_zonal_mean
INTEGER,           INTENT(IN)    :: fileid

CHARACTER (LEN=*), INTENT(INOUT) :: varname

INTEGER,           INTENT(OUT)   :: dimsize(max_dims)

! Local variables
INTEGER                        :: vid    ! variable id
LOGICAL                        :: lreco  ! true if contains record dim

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EM_VAR_CHECK_DIMS_3D_ZONAL'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! First just check if this routine should have been called
IF (.NOT. l_zonal_mean) THEN
  WRITE (umMessage,'(A,A,A)') 'em_var_check_dims_3d_zonal called for ' // &
                             TRIM(varname) // &
                            ' but with l_zonal_mean as .false. '
  CALL umPrint(umMessage,src='em_var_check_dims_3d_zonal')
  ecode = fileid
  CALL ereport ('EM_VAR_CHECK_DIMS_3D_ZONAL', ecode,                      &
      'l_zonal_mean was .false. so em_var_check_dims_3d_zonal: ' //       &
      ' should not have been called for ' // TRIM(varname))
END IF

! Emission files should be defined as (y,lev,dimension_3,t)
!
! Get information on variable. One of the two INOUT arguments
! of EM_GET_VAR_INFO is supplied (varname); assign an arbitrary
! negative value to the other one since it is not known yet (vid).
vid = -1
CALL em_get_var_info (fileid, vid, varname, vdimsize=dimsize, lrec=lreco)

! Check that dimensions (x, y, lev) match those of the model
IF ( dimsize(1) /= global_rows .OR.  &
    (dimsize(2) > 1 .AND. dimsize(2) /= model_levels)               .OR.  &
     dimsize(2) < 1                                                 .OR.  &
     dimsize(3) /= dimension_3) THEN 

  WRITE (umMessage,'(A,A)') 'ERROR: Variable dimensions do not match for ',&
                   TRIM(varname)
  CALL umPrint(umMessage,src='em_var_check_dims_3d_zonal')
  WRITE (umMessage,'(A,3(1x,I0))') 'Expected:',  &
                    global_rows, model_levels, dimension_3
  CALL umPrint(umMessage,src='em_var_check_dims_3d_zonal')
  WRITE (umMessage,'(A,3(1x,I0))') 'Found:', dimsize(1:3)
  CALL umPrint(umMessage,src='em_var_check_dims')
  ecode = fileid
  CALL ereport ('EM_VAR_CHECK_DIMS_3D_ZONAL', ecode,                       &
      'Dimension mismatch for: ' // TRIM(varname))
END IF

! Check that there is time record
IF (.NOT. lreco ) THEN 
  WRITE (umMessage,'(A,A)') 'No record dimension found for ', TRIM(varname)
  CALL umPrint(umMessage,src='em_var_check_dims')
  ecode = fileid
  CALL ereport ('EM_VAR_CHECK_DIMS_3D_ZONAL', ecode,                       &
      'No record dimension found for ' // TRIM(varname))
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE em_var_check_dims_3d_zonal

! --------------------------------------------------------------------------
! Description:
!   Check that dimensions of a variable in a NetCDF file match global sizes,
!   in the case of a full-atmosphere, time-dependent 3D distribution.
!
! Method:
! 1) Call EM_GET_VAR_INFO to get the dimensions and check they match.
! 2) If wrong dimensions or if there is no time record then stop with
!    call to EREPORT
! --------------------------------------------------------------------------
SUBROUTINE em_var_check_dims_4d (global_row_length, global_rows, model_levels, &
                                 dimension_4, fileid, varname, dimsize)

IMPLICIT NONE

! Subroutine arguments
INTEGER,           INTENT(IN)    :: global_row_length
INTEGER,           INTENT(IN)    :: global_rows
INTEGER,           INTENT(IN)    :: model_levels
INTEGER,           INTENT(IN)    :: dimension_4
INTEGER,           INTENT(IN)    :: fileid

CHARACTER (LEN=*), INTENT(INOUT) :: varname

INTEGER,           INTENT(OUT)   :: dimsize(max_dims)

! Local variables
INTEGER                        :: vid    ! variable id
LOGICAL                        :: lreco  ! true if contains record dim

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EM_VAR_CHECK_DIMS_4D'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! 4D emission files should be defined as (x,y,lev,dimension_4,t)
!
! Get information on variable. One of the two INOUT arguments
! of EM_GET_VAR_INFO is supplied (varname); assign an arbitrary
! negative value to the other one since it is not known yet (vid).
vid = -1
CALL em_get_var_info (fileid, vid, varname, vdimsize=dimsize, lrec=lreco)

! Check that dimensions (x, y, lev, dimension_4) match those of the model
IF ( dimsize(1) /= global_row_length .OR. dimsize(2) /= global_rows .OR.  &
    (dimsize(3) > 1 .AND. dimsize(3) /= model_levels)               .OR.  &
     dimsize(3) < 1                                                 .OR.  &
     dimsize(4) /= dimension_4 ) THEN

  WRITE (umMessage,'(A,A)') 'ERROR: Variable dimensions do not match for ',&
                   TRIM(varname)
  CALL umPrint(umMessage,src='em_var_check_dims')
  WRITE (umMessage,'(A,4(1x,I0))') 'Expected:', global_row_length,  &
                    global_rows, model_levels, dimension_4
  CALL umPrint(umMessage,src='em_var_check_dims')
  WRITE (umMessage,'(A,4(1x,I0))') 'Found:', dimsize(1:4)
  CALL umPrint(umMessage,src='em_var_check_dims')
  ecode = fileid
  CALL ereport ('EM_VAR_CHECK_DIMS_4D', ecode,                             &
      'Dimension mismatch for: ' // TRIM(varname))
END IF

! Check that there is time record
IF (.NOT. lreco ) THEN
  WRITE (umMessage,'(A,A)') 'No record dimension found for ', TRIM(varname)
  CALL umPrint(umMessage,src='em_var_check_dims')
  ecode = fileid
  CALL ereport ('EM_VAR_CHECK_DIMS_4D', ecode,                             &
      'No record dimension found for ' // TRIM(varname))
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE em_var_check_dims_4d

! --------------------------------------------------------------------------
! Description:
!   Get time information (in days and seconds) for specified record.
!
! Method:
!   The name and units attribute of the time variable have been read
!   previously during an EM_FOPEN call. They are stored in the
!   global variables 'tim_name' and 'ct_units', respectively. 
!   This routine uses that information to:
! 1) Get absolute date for record time (calculations done in file's calendar)
!    1a. Get time origin from ct_units in day/sec since 0000/01/01
!    1b. Add time offset to this to get record time in day/sec since 0000/01/01
!    1c. Convert this record time to absolute date
! 2) Change year if required for periodic ancils
! 3) Convert to tday/tsec since 0000/01/01 (in model's calendar) 
! --------------------------------------------------------------------------
SUBROUTINE em_get_time_rec (row_length, rows, fileid, irec, iupd_type, &
                            L_first, L_new_year, tday, tsec, L_notime)

IMPLICIT NONE

! Subroutine arguments

INTEGER, INTENT(IN)  :: row_length  ! horiz. model dimensions
INTEGER, INTENT(IN)  :: rows
INTEGER, INTENT(IN)  :: fileid      ! File id
INTEGER, INTENT(IN)  :: irec        ! Which record to check
INTEGER, INTENT(IN)  :: iupd_type   ! update_type: iupd_serial = 1
                                    !              iupd_cyclic = 2

LOGICAL, INTENT(IN)  :: L_first     ! First time
LOGICAL, INTENT(IN)  :: L_new_year  ! True if we need to get time
                                    ! for a new year

INTEGER, INTENT(OUT) :: tday        ! time (in days) found for irec
INTEGER, INTENT(OUT) :: tsec        ! time (in secs) found for irec

LOGICAL, INTENT(OUT) :: L_notime    ! True if no time record in file

! Local variables
INTEGER          :: nvars, undid
INTEGER          :: vid, ndim, dim_size(max_dims)
INTEGER          :: n

REAL             :: ftime       ! Fractional time read from the file

REAL             :: ftfrac      ! used to separate days & sec in ftime

! Dummy characters while decoding the time attribute "units"
CHARACTER(LEN=1) :: cd1         ! for '-', ' ' or ':'
CHARACTER(LEN=4) :: cdumm4      ! for "days"
CHARACTER(LEN=5) :: cdumm5      ! for "hours"
INTEGER          :: file_origt(6)  ! File origin time, also read from units
                                   ! 6 posit to hold: yyyy, mm, dd, hh, mm, ss
INTEGER          :: rec_time(7)    ! Record validity time
                                   ! 7 posit to hold: y,m,d,h,m,s,daynum

INTEGER          :: dayfilest, secfilest ! file orig time in days/sec
INTEGER          :: rec_yr    ! year of file record (imposed if cyclic)

LOGICAL          :: lreco       ! true if single record dim is found
LOGICAL          :: L_unit_dy   ! true if time is in days, otherwise in hours

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0  ! for DrHook
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EM_GET_TIME_REC'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Decode what is held in units for time variable, which should be either 
!    "days since yyyy-01-01 00:00:00"
!   "hours since yyyy-01-01 00:00:00"
! Only units allowed are 'days' or 'hours'.
IF ( ct_units(1:4) /= 'days' .AND. ct_units(1:4) /= 'hour' ) THEN
  ecode = fileid
  CALL ereport('EM_GET_TIME_REC',ecode,                      &
          'Time units should be in days or hours only')
END IF

L_unit_dy = .FALSE.                        ! units are hours
L_unit_dy = ( ct_units(1:4) == 'days' )    ! units are days

! Continue decoding the units attribute. Now get file_orig (1:6),
! which holds (year, month, day, hour, min, sec)
IF ( L_unit_dy ) THEN     ! assume start with the 4-letter word 'days'
  READ (ct_units,                                             &
      '(A4,1x,A5,1x,I4,A1,I2,A1,I2,1x,I2,A1,I2,A1,I2)')       &
       cdumm4, cdumm5,                                        & ! 'days since'
       file_origt(1), cd1, file_origt(2), cd1, file_origt(3), & ! 'yyyy-mm-dd'
       file_origt(4), cd1, file_origt(5), cd1, file_origt(6)    ! 'hh:mm:ss'

ELSE                       ! assume start with the 5-letter word 'hours'
  READ (ct_units,                                            &
     '(A5,1x,A5,1x,I4,A1,I2,A1,I2,1x,I2,A1,I2,A1,I2)')       &
      cdumm5, cdumm5,                                        & ! 'hours since'
      file_origt(1), cd1, file_origt(2), cd1, file_origt(3), & ! 'yyyy-mm-dd'
      file_origt(4), cd1, file_origt(5), cd1, file_origt(6)    ! 'hh:mm:ss'
END IF

IF ((file_origt(1) < 1) .AND. (iupd_type == iupd_cyclic) ) THEN
  file_origt(1) = 1
  IF (L_first .AND. irec == 1) THEN
    IF (PrintStatus >= PrStatus_Diag) THEN
      umMessage = "File basis year for periodic file increased to 0001 " // &
                  "because time2sec doesn't handle smaller years."
      CALL umPrint(umMessage,src='em_get_time_rec')
    END IF
  END IF
END IF

! Get file data start time in days/secs, e.g. from 0000-01-01.
! Pass file_orig (1:6) to time2sec and get the starting time as days and 
! secs (in the file's calendar) since basis time: dayfilest & secfilest
! DEPENDS ON: time2sec
CALL time2sec (file_origt(1), file_origt(2), file_origt(3),  &
               file_origt(4), file_origt(5), file_origt(6),  &
               0, 0,                                         &
               dayfilest, secfilest, lcal360file )

! Print out time origin in NetCDF file if first time and extra diags.
IF (L_first .AND. irec == 1) THEN
  IF (PrintStatus >= PrStatus_Diag) THEN
    WRITE (umMessage,'(A,I4,5(1x,I2),2(1x,I8))') 'EM_GET_TIME_REC: '// &
           'File starts ', file_origt(1:6), dayfilest, secfilest
    CALL umPrint(umMessage,src='em_get_time_rec')
  END IF
END IF

! ----------------------------------------------------------------------------
! Get record for requested time (i.e. record irec from variable
! with name tim_name) and store it in ftime
IF (ALLOCATED(time_axis)) THEN
  ! check irec is within the range of time_axis dimension
  IF (irec < 1 .OR. irec > SIZE(time_axis)) L_notime = .TRUE.
ELSE
  ! time_axis not allocated - so file doesn't have time coordinate
  L_notime = .TRUE.
END IF
IF ( L_notime ) THEN     ! no time record found in file
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
  RETURN
END IF

ftime = time_axis(irec)

! If time is in hours then convert to days
IF (.NOT. L_unit_dy) ftime = ftime * rday_per_hour

! Now calculate time in full days & seconds: tday & tsec.
! For that check if ftime value is full days or decimal.
! This is done in a crude way, which can lead to very
! low decimal diffs at REAL8. Seconds may be rounded off,
! but this is not expected to have any impact.
ftfrac = REAL  ( INT(ftime) )                ! Fractional day
tsec   = INT   ( (ftime-ftfrac) * isec_per_day )
                                             ! Fractional part of ftime in secs
tday   = INT   (ftime)                       ! Full days

! Finally, get the final time in days and seconds by adding the starting
! days and seconds that we read from the attribute 'units'.
tday = tday + dayfilest
tsec = tsec + secfilest

! Convert tday,tsec (in file's calendar) to absolute year, month, day, etc 
! to allow year change for periodic files, or to allow conversion from file's
! calendar to model calendar if required.
! DEPENDS ON: sec2time
CALL sec2time(tday, tsec, &
              0, 0, &
              rec_time(1), rec_time(2), rec_time(3), &
              rec_time(4), rec_time(5), rec_time(6), &
              rec_time(7), lcal360file)

! If file is periodic then ignore year. Note that the year is set to
! cyclic_year in order to be consistent with get_emfile_rec.
IF (iupd_type == iupd_cyclic) THEN
  rec_yr = cyclic_year

  ! Add one year if indicated by GET_EMFILE_REC (in UKCA_EMISS_MOD)
  IF (L_new_year) THEN
    rec_yr = rec_yr + 1
  END IF

  rec_time(1) = rec_yr

END IF

! Finally convert absolute date to tday, tsec since date 0000-01-01 in model's
! calendar.
! DEPENDS ON: time2sec
CALL time2sec(rec_time(1), rec_time(2), rec_time(3),  &
              rec_time(4), rec_time(5), rec_time(6),  &
              0, 0,                                   &
              tday, tsec, lcal360)

! tday and tsec can now be passed to other emission routines, where
! they will be compared with the current model time in order to
! perform interpolations of the emission fields.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE em_get_time_rec

!-------------------------------------------------------------------------
! Description:
!   Read data values of character variable in NetCDF file.
!   Part of the interface block EM_GET_DATA used for overloading; this
!   subroutine is called when the variable is of type character.
!
! Method:
!   Get variable ID (varid), number of dimensions in variable (ndim), 
!   and vector of dimension sizes (dim_size) by calling EM_GET_VAR_INFO.
!   Then get domain/indices to extract and call NF90_GET_VAR to get the 
!   variable values: cvalue. Values are read by PE0 and broadcast
!   to all PEs.
!   If the netCDF variable has a time dimension, then only the data for
!   time index = trec are read. In this case the netCDF variable must be a 1D 
!   array of strings, and a single string is read. It is assumed
!   that time is the last dimension for this variable: (str_len, time).
!-------------------------------------------------------------------------
SUBROUTINE em_get_data_char (row_length, rows, fileid, trec,  &
                             var_name, cvalue, L_norec)

IMPLICIT NONE

!  Subroutine arguments
INTEGER,          INTENT(IN)    :: row_length ! model horiz dimensions
INTEGER,          INTENT(IN)    :: rows
INTEGER,          INTENT(IN)    :: fileid    ! file id
INTEGER,          INTENT(IN)    :: trec      ! time record wanted

CHARACTER(LEN=*), INTENT(IN)    :: var_name  ! variable name

CHARACTER(LEN=*), INTENT(INOUT) :: cvalue(*) ! value out

LOGICAL,          INTENT(OUT)   :: L_norec   ! T if record not in file (EOF?)

!  Local variables

INTEGER(KIND=integer32) :: nstatus ! netCDF return code

!  Vector specifying the index in the variable from which
!  the first data value will be read along each dimension
INTEGER(KIND=integer32) :: start   (max_dims)

!  Vector specifying the number of counters to read along each dimension
INTEGER(KIND=integer32) :: counter (max_dims)

INTEGER(KIND=integer32) :: varid32
INTEGER                 :: varid, ndim, dim_size(max_dims), n, ilen, i_get_dims
LOGICAL                 :: lreco

CHARACTER (LEN=varname_len) :: varname_local

! For Parallel I/O
INTEGER :: my_comm
INTEGER :: icode

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EM_GET_DATA_CHAR'

!  End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Obtain the communication group for this run
CALL gc_get_communicator (my_comm, icode)

! Get information on variable. One of the two INOUT arguments of
! EM_GET_VAR_INFO is supplied (varname_local); assign an arbitrary
! negative value to the other one since it is not known yet (varid).
varname_local = var_name    ! need to pass INOUT argument
varid         = -1
CALL em_get_var_info (fileid,      varid,             varname_local, &
                      ndims=ndim,  vdimsize=dim_size, lrec=lreco)

! Compute ilen: total size needed for upper bound of assumed-size
! array used in nf90_get_var. The total array size is the product
! of all dimension sizes (excluding last dimension if lreco=true).
IF (lreco) THEN
  i_get_dims = ndim - 1
ELSE
  i_get_dims = ndim
END IF
ilen = 1
DO n =1, i_get_dims
  ilen = ilen * dim_size(n)
END DO
! Check that this length is no longer than the string it  will be read into
IF (ilen > LEN(cvalue)) THEN
  ecode = ilen
  CALL ereport ('EM_GET_DATA_CHAR', ecode,                       &
      'NetCDF character variable too long in: ' // TRIM(var_name))
END IF

! Only PE0 reads data values from the file
IF ( mype == 0 ) THEN
  ! Default indices/domain to extract.
  ! Default is to get complete data array from NetCDF file so set
  ! all dimensions equal to the size of the variable in the file.
  start   (1:ndim) = 1
  counter (1:ndim) = dim_size (1:ndim)

  IF (lreco) THEN             ! If var contains time dimension
    L_norec = .TRUE.

    ! Record out of bounds and return that no record found
    IF ( trec < 1 .OR. trec > dim_size(ndim) ) THEN
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,       &
                              zhook_handle)
      RETURN
    END IF
                                       
    L_norec        = .FALSE.  ! Indicate to main routine that time rec found

    ! Read a single time from the file
    start   (ndim) = trec     ! last index corresponds to record number
    counter (ndim) = 1        ! read only one record
  END IF

  nfid    = fileid
  varid32 = varid    ! 32-bit integer needed
  nstatus = nf90_get_var (nfid, varid32, cvalue(1:ilen),   &
                         start(1:ndim), counter(1:ndim))
  CALL nd_error (nstatus, 'EM_GET_DATA_CHAR', var_name)
END IF   ! I/O PE

! Distribute data to all other PEs
icode = 0
CALL mpl_bcast (cvalue, LEN(cvalue), mpl_character, 0, my_comm, icode)
IF ( icode /= 0 ) THEN
  ecode= fileid
  CALL ereport ('EM_GET_DATA_CHAR', ecode, 'Broadcast Error')
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE em_get_data_char

!-------------------------------------------------------------------------
! Description:
!   Read data values of integer variable in NetCDF file.
!   Part of the interface block EM_GET_DATA used for overloading; this
!   subroutine is called when the variable is of type integer.
!
! Method:
!   Get variable ID (varid), number of dimensions in variable (ndim), 
!   and vector of dimension sizes (dim_size) by calling EM_GET_VAR_INFO.
!   Then get domain/indices to extract and call NF90_GET_VAR to get the 
!   variable values: ivalue. Values are read by PE0 and broadcast
!   to all PEs.
!   If the netCDF variable has a time dimension, then only the data for
!   time index = trec are read. In this case the netCDF variable is 1D and
!   a single value is read, or is 2D and a 1D slice is read. It is assumed
!   that time is the last dimension for this variable.
!-------------------------------------------------------------------------
SUBROUTINE em_get_data_int (row_length, rows, fileid, trec,  &
                            var_name, ivalue, L_norec)

IMPLICIT NONE

! Subroutine arguments
INTEGER,           INTENT (IN)    :: row_length ! model horiz dimensions
INTEGER,           INTENT (IN)    :: rows
INTEGER,           INTENT (IN)    :: fileid   ! file id
INTEGER,           INTENT (IN)    :: trec     ! time record wanted

CHARACTER (LEN=*), INTENT (IN)    :: var_name ! variable name

INTEGER,           INTENT (INOUT) :: ivalue (:) ! variable values to read

LOGICAL,           INTENT (OUT)   :: L_norec  ! T if rec not in file, EOF

! Local variables

INTEGER (KIND=integer32) :: nstatus ! netCDF return code

!  Vector specifying the index in the variable from which
!  the first data value will be read along each dimension
INTEGER (KIND=integer32) :: start   (max_dims)

!  Vector specifying the number of counters to read along each dimension
INTEGER (KIND=integer32) :: counter (max_dims)

INTEGER (KIND=integer32) :: varid32
INTEGER                  :: varid, ndim, dim_size(max_dims), n, ilen, i_get_dims
LOGICAL                  :: lreco    ! true if contains time record dim

CHARACTER (LEN=varname_len) :: varname_local

! For Parallel I/O
INTEGER :: my_comm
INTEGER :: icode

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EM_GET_DATA_INT'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Obtain the communication group for this run
CALL gc_get_communicator (my_comm, icode)

! Get information on variable. One of the two INOUT arguments of
! EM_GET_VAR_INFO is supplied (varname_local); assign an arbitrary
! negative value to the other one since it is not known yet (varid).
varname_local = var_name    ! need to pass INOUT argument
varid         = -1
CALL em_get_var_info (fileid,     varid,             varname_local, &
                      ndims=ndim, vdimsize=dim_size, lrec=lreco)

! Compute ilen: total size needed for upper bound of assumed-size
! array used in nf90_get_var. The total array size is the product
! of all dimension sizes (excluding last dimension if lreco=true).
IF (lreco) THEN
  i_get_dims = ndim - 1
ELSE
  i_get_dims = ndim
END IF
ilen = 1
DO n =1, i_get_dims
  ilen = ilen * dim_size(n)
END DO

! Only PE0 reads data values from the file
IF ( mype == 0 ) THEN
  ! Default indices/domain to extract.
  ! Default is to get complete data array from NetCDF file so set
  ! all dimensions equal to the size of the variable in the file.
  start   (1:ndim) = 1
  counter (1:ndim) = dim_size(1:ndim)

  ! If var contains time dimension then select record to read
  IF (lreco) THEN
    L_norec = .TRUE.

    ! Record out of bounds and return that no record found
    IF ( trec < 1 .OR. trec > dim_size(ndim) ) THEN
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,       &
                              zhook_handle)
      RETURN
    END IF

    L_norec        = .FALSE.  ! Indicate to main routine that time rec found

    ! Read a single time from the file
    start   (ndim) = trec     ! last index corresponds to record number
    counter (ndim) = 1        ! read only one record
  END IF

  ! Check that array with data values to read (ivalue) has the right dimensions
  IF (SIZE (ivalue,1) /= counter (1)) THEN
    ecode = 1
    CALL ereport ('EM_GET_DATA_INT', ecode,           &
        'Dimension mismatch for: ' // TRIM(var_name))
  END IF

  !  Get values of the variable as ivalue, specifying what to read
  !  with the optional start and counter arguments
  nfid    = fileid
  varid32 = varid    ! 32-bit integer needed
  nstatus = nf90_get_var (nfid, varid32, ivalue(1:ilen),   &
                          start(1:ndim), counter (1:ndim))
  CALL nd_error (nstatus, 'EM_GET_DATA_INT', var_name)

END IF

! Check that the array to be broadcast has the same size as the buffer length
! given to mpl_bcast
IF (SIZE(ivalue,1) /= ilen) THEN
  ecode = ilen
  CALL ereport ('EM_GET_DATA_INT', ecode,                       &
      'Dimension mismatch for: ' // TRIM(var_name))
END IF

! Distribute data to all other PEs
icode = 0
CALL mpl_bcast (ivalue, ilen, mpl_integer, 0, my_comm, icode)
IF ( icode /= 0 ) THEN
  ecode = fileid
  CALL ereport ('EM_GET_DATA_INT', ecode, 'Broadcast Error')
END IF


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE em_get_data_int

!-------------------------------------------------------------------------
! Description:
!   Read data values of one-dimensional real variable in NetCDF file.
!   Part of the interface block EM_GET_DATA used for overloading; this
!   subroutine is called when the variable is of type real 1-dim.
!
! Method:
!   Get variable ID (varid), number of dimensions in variable (ndim), 
!   and vector of dimension sizes (dim_size) by calling EM_GET_VAR_INFO.
!   Then get domain/indices to extract and call NF90_GET_VAR to get the 
!   variable values: rvalue. Values are read by PE0 and broadcast
!   to all PEs.
!   If the netCDF variable has a time dimension, then only the data for
!   time index = trec are read. In this case the netCDF variable is 1D and
!   a single value is read, or is 2D and a 1D slice is read. It is assumed
!   that time is the last dimension for this variable.
!-------------------------------------------------------------------------
SUBROUTINE em_get_data_real (row_length, rows, fileid, trec, &
                             var_name, rvalue, L_norec)

IMPLICIT NONE

! Subroutine arguments
INTEGER,           INTENT (IN)    :: row_length ! model horiz dimensions
INTEGER,           INTENT (IN)    :: rows
INTEGER,           INTENT (IN)    :: fileid     ! file id
INTEGER,           INTENT (IN)    :: trec       ! time record wanted

CHARACTER (LEN=*), INTENT (IN)    :: var_name   ! variable name

REAL,              INTENT (INOUT) :: rvalue (:) ! variable values to read

LOGICAL,           INTENT (OUT)   :: L_norec    ! T if rec not in file, EOF

! Local variables

INTEGER (KIND=integer32) :: nstatus ! netCDF return code

!  Vector specifying the index in the variable from which
!  the first data value will be read along each dimension
INTEGER (KIND=integer32) :: start   (max_dims)

!  Vector specifying the number of counters to read along each dimension
INTEGER (KIND=integer32) :: counter (max_dims)

INTEGER (KIND=integer32) :: varid32
INTEGER                  :: varid, ndim, dim_size(max_dims), n, ilen, i_get_dims
LOGICAL                  :: lreco    ! true if contains time record dim

CHARACTER (LEN=varname_len) :: varname_local

! For Parallel I/O
INTEGER :: my_comm
INTEGER :: icode

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EM_GET_DATA_REAL'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Obtain the communication group for this run
CALL gc_get_communicator (my_comm, icode)

! Get information on variable. One of the two INOUT arguments of
! EM_GET_VAR_INFO is supplied (varname_local); assign an arbitrary
! negative value to the other one since it is not known yet (varid).
varname_local = var_name    ! need to pass INOUT argument
varid         = -1
CALL em_get_var_info (fileid,     varid,             varname_local,  &
                      ndims=ndim, vdimsize=dim_size, lrec=lreco)

! Compute ilen: total size needed for upper bound of assumed-size
! array used in nf90_get_var. The total array size is the product
! of all dimension sizes (excluding last dimension if lreco=true).
IF (lreco) THEN
  i_get_dims = ndim - 1
ELSE
  i_get_dims = ndim
END IF
ilen = 1
DO n =1, i_get_dims
  ilen = ilen * dim_size(n)
END DO

! Only PE0 reads data values from the file
IF ( mype == 0 ) THEN
  ! Default indices/domain to extract.
  ! Default is to get complete data array from NetCDF file so set
  ! all dimensions equal to the size of the variable in the file.
  start   (1:ndim) = 1
  counter (1:ndim) = dim_size (1:ndim)

  ! If var contains time dimension then select record to read
  IF (lreco) THEN
    L_norec = .TRUE.

    ! Record out of bounds and return that no record found
    IF ( trec < 1 .OR. trec > dim_size(ndim) ) THEN
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,       &
                              zhook_handle)
      RETURN
    END IF

    L_norec        = .FALSE.  ! Indicate to main routine that time rec found

    ! Read a single time from the file
    start   (ndim) = trec     ! last index corresponds to record number
    counter (ndim) = 1        ! read only one record
  END IF

  ! Check that array with data values to read (rvalue) has the right dimensions
  IF (SIZE (rvalue,1) /= counter (1)) THEN
    ecode = 1
    CALL ereport ('EM_GET_DATA_REAL', ecode,                       &
        'Dimension mismatch for: ' // TRIM(var_name))
  END IF

  !  Get values of the variable as rvalue, specifying what to read
  !  with the optional start and counter arguments
  nfid    = fileid
  varid32 = varid    ! 32-bit integer needed
  nstatus = nf90_get_var (nfid, varid32, rvalue, start   (1:ndim), &
                                                 counter (1:ndim))
  CALL nd_error (nstatus, 'EM_GET_DATA_REAL', var_name)
END IF   ! I/O PE

! Check that the array to be broadcast has the same size as the buffer length
! given to mpl_bcast
IF (SIZE(rvalue,1) /= ilen) THEN
  ecode = ilen
  CALL ereport ('EM_GET_DATA_REAL', ecode,                       &
      'Dimension mismatch for: ' // TRIM(var_name))
END IF

! Distribute data to all other PEs
icode = 0
CALL mpl_bcast (rvalue, ilen, mpl_real, 0, my_comm, icode)
IF ( icode /= 0 ) THEN
  ecode = fileid
  CALL ereport ('EM_GET_DATA_REAL', ecode, 'Broadcast Error')
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE em_get_data_real

!-------------------------------------------------------------------------
! Description:
!   Read data values of two-dimensional real variable in NetCDF file.
!   Part of the interface block EM_GET_DATA used for overloading; this
!   subroutine is called when the variable is of type real 2-dim.
!
! Method:
!   Get variable ID (varid), number of dimensions in variable (ndim), 
!   and vector of dimension sizes (dim_size) by calling EM_GET_VAR_INFO.
!   Then get domain/indices to extract and call NF90_GET_VAR to get the 
!   variable values, which are read by PEO as global_var. Then the values
!   are broadcast to all PEs and returned as rvalue. 
!   Assumes (lon, lat).
!-------------------------------------------------------------------------
SUBROUTINE em_get_data_real2D (row_length, rows, fileid, trec, &
                               var_name, rvalue, L_norec)

USE nlsizes_namelist_mod,  ONLY:global_row_length, global_rows

IMPLICIT NONE

! Subroutine arguments
INTEGER,           INTENT (IN)    :: row_length ! model horiz dimensions
INTEGER,           INTENT (IN)    :: rows
INTEGER,           INTENT (IN)    :: fileid     ! file id
INTEGER,           INTENT (IN)    :: trec       ! time record wanted

CHARACTER (LEN=*), INTENT (IN)    :: var_name   ! variable name

REAL,              INTENT (INOUT) :: rvalue (:,:) ! variable values to read

LOGICAL,           INTENT (OUT)   :: L_norec    ! T if rec not in file, EOF

! Local variables

INTEGER (KIND=integer32) :: nstatus ! netCDF return code

!  Vector specifying the index in the variable from which
!  the first data value will be read along each dimension
INTEGER (KIND=integer32) :: start   (max_dims)

!  Vector specifying the number of counters to read along each dimension
INTEGER (KIND=integer32) :: counter (max_dims)

INTEGER (KIND=integer32) :: varid32
INTEGER                  :: varid, ndim, dim_size(max_dims)
LOGICAL                  :: lreco    ! true if contains time record dim

CHARACTER (LEN=varname_len) :: varname_local

! Variable to hold file data (lon, lat)
REAL                     :: global_var(global_row_length,global_rows)

INTEGER:: icode
CHARACTER (LEN=errormessagelength) :: cmessage  ! error message

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EM_GET_DATA_REAL2D'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Get information on variable. One of the two INOUT arguments of
! EM_GET_VAR_INFO is supplied (varname_local); assign an arbitrary
! negative value to the other one since it is not known yet (varid).
varname_local = var_name    ! need to pass INOUT argument
varid         = -1
CALL em_get_var_info (fileid,      varid,             varname_local, &
                      ndims=ndim,  vdimsize=dim_size, lrec=lreco)

! Get values of variable from PE0 and store in global_var
IF ( mype == 0 ) THEN
  ! Default indices/domain to extract.
  ! Default is to get complete data array from NetCDF file so set
  ! all dimensions equal to the size of the variable in the file.
  start   (1:ndim) = 1
  counter (1:ndim) = dim_size(1:ndim)

  ! If var contains time dimension then select record to read
  IF (lreco) THEN
    L_norec = .TRUE.

    ! Record out of bounds and return that no record found
    IF ( trec > dim_size(ndim) ) THEN
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,       &
                              zhook_handle)
      RETURN
    END IF

    L_norec        = .FALSE.  ! Indicate to main routine that time rec found

    ! Read a single time from the file
    start   (ndim) = trec     ! last index corresponds to record number
    counter (ndim) = 1        ! read only one record
  END IF

  ! Check that data in NetCDF file matches Global Dimensions
  IF (counter (1) /= global_row_length .OR.                               &
      counter (2) /= global_rows ) THEN
    ecode = 1
    CALL ereport ('EM_GET_DATA_REAL2D', ecode,                            &
         'Data does not match global dimensions for: ' // TRIM(var_name))
  END IF

  !  Get values of the variable as global_var, specifying what to read
  !  with the optional start and counter arguments
  nfid    = fileid
  varid32 = varid    ! 32-bit integer needed
  nstatus = nf90_get_var (nfid, varid32, global_var, start (1:ndim),      &
                          counter (1:ndim))
  CALL nd_error (nstatus, 'EM_GET_DATA_REAL2D', var_name)
END IF    ! I/O PE

icode = 0
CALL gc_gsync(nproc,icode)   ! Synchronise PEs before scatter

! Distribute data to all other PEs
! DEPENDS ON: scatter_field
CALL scatter_field (rvalue, global_var, row_length, rows,   &
                    global_row_length, global_rows,         &
                    fld_type_p, halo_type_no_halo,          &
                    0, gc_all_proc_group)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE em_get_data_real2D

!-------------------------------------------------------------------------
! Description:
!   Read data values of two-dimensional real variable in NetCDF file,
!   (zonal mean data with dimensions: lat, lev), then expand these into
!   3D data field applied across all longitudes (1:global_row_length).   
!   Part of the interface block EM_GET_DATA used for overloading; this
!   subroutine is called when the variable is of type real 3-dim
!   but logical l_zonal_mean is supplied, and set to true.
!
! Method:
!   Get variable ID (varid), number of dimensions in variable (ndim),
!   and vector of dimension sizes (dim_size) by calling EM_GET_VAR_INFO.
!   Then get domain/indices to extract and call NF90_GET_VAR to get the
!   variable values, which are read by PEO as global_var_2d. Then the values
!   are expanded in global_var_3d, then broadcast to all PEs and returned as
!   rvalue. Assumes (lon, lat).
!-------------------------------------------------------------------------
SUBROUTINE em_get_data_real2D_zonal (row_length, rows, fileid, trec, &
                               l_zonal_mean, var_name, rvalue, L_norec)

USE nlsizes_namelist_mod,  ONLY: global_row_length, global_rows

IMPLICIT NONE

! Subroutine arguments
INTEGER,           INTENT (IN)    :: row_length
INTEGER,           INTENT (IN)    :: rows
INTEGER,           INTENT (IN)    :: fileid     ! file id
INTEGER,           INTENT (IN)    :: trec       ! time record wanted
LOGICAL,           INTENT (IN)    :: l_zonal_mean ! file contains zonal means

CHARACTER (LEN=*), INTENT (IN)    :: var_name   ! variable name

REAL,              INTENT (INOUT) :: rvalue (:,:,:) ! variable values to read

LOGICAL,           INTENT (OUT)   :: L_norec    ! T if rec not in file, EOF

! Local variables

INTEGER (KIND=integer32) :: nstatus ! netCDF return code

!  Vector specifying the index in the variable from which
!  the first data value will be read along each dimension
INTEGER (KIND=integer32) :: start   (max_dims)

!  Vector specifying the number of counters to read along each dimension
INTEGER (KIND=integer32) :: counter (max_dims)

INTEGER (KIND=integer32) :: varid32
INTEGER                  :: varid, ndim, dim_size(max_dims)
LOGICAL                  :: lreco    ! true if contains time record dim

CHARACTER (LEN=varname_len) :: varname_local

! Variable to hold file data (lat, lev)
REAL, ALLOCATABLE        :: global_var_2d(:,:)

! Variable to hold expanded data (lon, lat, lev)
REAL, ALLOCATABLE        :: global_var_3d(:,:,:)


INTEGER :: i, k      ! loop counters
INTEGER:: icode
CHARACTER (LEN=errormessagelength) :: cmessage  ! error message

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EM_GET_DATA_REAL2D_ZONAL'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! First just check if this routine should have been called
IF (.NOT. l_zonal_mean) THEN
  WRITE (umMessage,'(A,A,A)') 'EM_GET_DATA_REAL2D_ZONAL called for ' // &
                             TRIM(var_name) // &
                            ' but with l_zonal_mean as .false. '
  CALL umPrint(umMessage,src='em_get_data_real2D_zonal')
  ecode = fileid
  CALL ereport ('EM_GET_DATA_REAL2D_ZONAL', ecode,                      &
      'l_zonal_mean was .false. so EM_GET_DATA_REAL2D_ZONAL: ' //       &
      ' should not have been called for ' // TRIM(var_name))
END IF

! Get information on variable. One of the two INOUT arguments of
! EM_GET_VAR_INFO is supplied (varname_local); assign an arbitrary
! negative value to the other one since it is not known yet (varid).
varname_local = var_name    ! need to pass INOUT argument
varid         = -1
CALL em_get_var_info (fileid,      varid,             varname_local, &
                      ndims=ndim,  vdimsize=dim_size, lrec=lreco)

! Allocate the variables that will hold the data.
! The zonal mean data will be copied from global_var_2d into global_var_4d.
! PE0 will fill the array and broadcast to other PEs.
ALLOCATE ( global_var_2d (dim_size(1), dim_size(2) ))
ALLOCATE ( global_var_3d (global_row_length, dim_size(1), dim_size(2) ))

! Get values of variable from PE0 and store in global_var_2d
IF ( mype == 0 ) THEN
  ! Default indices/domain to extract.
  ! Default is to get complete data array from NetCDF file so set
  ! all dimensions equal to the size of the variable in the file.
  start   (1:ndim) = 1
  counter (1:ndim) = dim_size(1:ndim)

  ! If var contains time dimension then select record to read
  IF (lreco) THEN
    L_norec = .TRUE.

    ! Record out of bounds and return that no record found
    IF ( trec > dim_size(ndim) ) THEN
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, &
                              zhook_out, zhook_handle)
      RETURN
    END IF
    L_norec        = .FALSE.  ! Indicate to main routine that time rec found

    ! Read a single time from the file
    start   (ndim) = trec     ! last index corresponds to record number
    counter (ndim) = 1        ! read only one record
  END IF

  ! Check that data in NetCDF file matches Global Dimensions
  IF (dim_size (1)   /= global_rows  .OR.                               &
      SIZE(rvalue,3) /= dim_size(2)) THEN
    ecode = 1
    CALL ereport ('EM_GET_DATA_REAL2D_ZONAL', ecode,                      &
         'Data does not match global dimensions for: ' // TRIM(var_name))
  END IF

  !  Get values of the variable as global_var_2d, specifying what to read
  !  with the optional start and counter arguments
  nfid    = fileid
  varid32 = varid    ! 32-bit integer needed
  nstatus = nf90_get_var (nfid, varid32, global_var_2d, start (1:ndim),      &
                          counter (1:ndim))
  CALL nd_error (nstatus, 'EM_GET_DATA_REAL2D_ZONAL', var_name)
END IF    ! I/O PE

! Expand zonal mean out to all longitudes
DO i = 1, global_row_length
   global_var_3d(i,:,:) = global_var_2d(:,:)
END DO

icode = 0
CALL gc_gsync(nproc,icode)   ! Synchronise PEs before scatter

! Distribute data to all other PEs
DO k = 1, dim_size(2)
  ! DEPENDS ON: scatter_field
  CALL scatter_field (rvalue(:,:,k), global_var_3d(:,:,k),    &
                      row_length, rows,                       &
                      global_row_length, global_rows,         &
                      fld_type_p, halo_type_no_halo,          &
                      0, gc_all_proc_group)
END DO

DEALLOCATE(global_var_3d)
DEALLOCATE(global_var_2d)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE em_get_data_real2D_zonal

!-------------------------------------------------------------------------
! Description:
!   Read data values of three-dimensional real variable in NetCDF file.
!   Part of the interface block EM_GET_DATA used for overloading; this
!   subroutine is called when the variable is of type real 3-dim.
!
! Method:
!   Get variable ID (varid), number of dimensions in variable (ndim), 
!   and vector of dimension sizes (dim_size) by calling EM_GET_VAR_INFO.
!   Then get domain/indices to extract and call NF90_GET_VAR to get the 
!   variable values, which are read by PEO as global_var. The values
!   are broadcast to all PEs, looping over levels, and returned as rvalue.
!   Assumes (lon, lat, lev).
!-------------------------------------------------------------------------
SUBROUTINE em_get_data_real3D (row_length, rows, fileid, trec,   &
                               var_name, rvalue, L_norec)

USE nlsizes_namelist_mod,  ONLY:global_row_length, global_rows

IMPLICIT NONE

! Subroutine arguments
INTEGER,           INTENT (IN)    :: row_length ! model horiz dimensions
INTEGER,           INTENT (IN)    :: rows
INTEGER,           INTENT (IN)    :: fileid    ! file id
INTEGER,           INTENT (IN)    :: trec      ! time record wanted

CHARACTER (LEN=*), INTENT (IN)    :: var_name  ! variable name

! variable values to read (lon, lat, lev)
REAL,              INTENT (INOUT) :: rvalue (:,:,:)

LOGICAL,           INTENT (OUT)   :: L_norec     ! T if rec not in file, EOF

! Local variables

INTEGER (KIND=integer32) :: nstatus ! netCDF return code

! Vector specifying the index in the variable from which
! the first data value will be read along each dimension.
INTEGER (KIND=integer32) :: start   (max_dims)

! Vector specifying the number of counters to read along each dimension
INTEGER (KIND=integer32) :: counter (max_dims)

INTEGER (KIND=integer32) :: varid32
INTEGER                  :: varid, ndim, dim_size(max_dims)
LOGICAL                  :: lreco    ! true if contains time record dim

CHARACTER (LEN=varname_len) :: varname_local

! Variable to hold file data
REAL, ALLOCATABLE        :: global_var(:,:,:)

INTEGER:: icode, k
CHARACTER (LEN=errormessagelength) :: cmessage  ! error message

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EM_GET_DATA_REAL3D'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Get information on variable (need to do this from all PEs so that
! global_var can be allocated with the right 3rd dimension)
!
! One of the two INOUT arguments of EM_GET_VAR_INFO is supplied
! (varname_local); assign an arbitrary negative value to the other 
! one since it is not known yet (varid).
varname_local = var_name    ! need to pass INOUT argument
varid         = -1
CALL em_get_var_info (fileid,      varid,             varname_local,   &
                      ndims=ndim,  vdimsize=dim_size, lrec=lreco)

! Check that data in NetCDF file matches Global Dimensions and 
! output variable rvalue has correct space for levels
IF (dim_size (1)   /= global_row_length .OR.            &
    dim_size (2)   /= global_rows       .OR.            &
    SIZE(rvalue,3) /= dim_size(3))      THEN
    ecode = 1
    CALL ereport ('EM_GET_DATA_REAL3D', ecode,          &
       'Dimension mismatch for: ' // TRIM(var_name))
END IF

! Allocate the variable that will hold the data. It will be filled
! by PE0 and broadcast to other PEs.
ALLOCATE ( global_var (dim_size(1), dim_size(2), dim_size(3)) )

! Get values of variable from PE0 and store in global_var
IF ( mype == 0 ) THEN
  ! Default indices for the first data value and number of counters
  ! along two dimensions: Global lon and lat. Note that counter(3) will
  ! often be different from the number of model vertical levels, as
  ! it happens in the case of surface emissions.
  start   (1:3) = 1
  counter (1)   = dim_size (1)
  counter (2)   = dim_size (2)
  counter (3)   = dim_size (3)

  ! If var contains time dimension then select record to read
  IF (lreco) THEN
    L_norec = .TRUE.

    ! Record out of bounds and return that no record found
    IF ( trec > dim_size(ndim) ) THEN
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,       &
                              zhook_handle)
      RETURN
    END IF

    L_norec        = .FALSE.  ! Indicate to main routine that time rec found

    ! Read a single time from the file
    start   (ndim) = trec     ! last index corresponds to record number
    counter (ndim) = 1        ! read only one record
  END IF

  ! Get values of the variable as global_var, specifying what to read
  ! with the optional start and counter arguments
  nfid    = fileid
  varid32 = varid    ! 32-bit integer needed
  nstatus = nf90_get_var (nfid, varid32, global_var, start  (1:ndim),  &
                                                     counter(1:ndim))
  CALL nd_error (nstatus, 'EM_GET_DATA_REAL3D', var_name)
END IF    ! I/O PE

icode = 0
CALL gc_gsync(nproc,icode)   ! Synchronise PEs before scatter

! Distribute data to all PEs, looping over levels
DO k = 1, dim_size(3)
  ! DEPENDS ON: scatter_field
  CALL scatter_field(rvalue(:,:,k), global_var(:,:,k),       &
                   row_length, rows,                         &
                   global_row_length, global_rows,           &
                   fld_type_p, halo_type_no_halo,            &
                   0, gc_all_proc_group)
END DO

DEALLOCATE(global_var)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE em_get_data_real3D

!-------------------------------------------------------------------------
! Description:
!   Read data values of three-dimensional real variable in NetCDF file,
!   (zonal mean data with dimensions: lat, lev, 3rd dimen), then expand these
!   into 4D data field applied across all longitudes (1:global_row_length).   
!   Part of the interface block EM_GET_DATA used for overloading; this
!   subroutine is called when the variable is of type real 4-dim
!   but logical l_zonal_mean is supplied, and set to true.
!
! Method:
!   Get variable ID (varid), number of dimensions in variable (ndim),
!   and vector of dimension sizes (dim_size) by calling EM_GET_VAR_INFO.
!   Then get domain/indices to extract and call NF90_GET_VAR to get the
!   variable values, which are read by PEO as global_var_3d. Then the values
!   are expanded in global_var_4d, then broadcast to all PEs and returned as
!   rvalue. Assumes (lon, lat).
!-------------------------------------------------------------------------
SUBROUTINE em_get_data_real3D_zonal (row_length, rows, fileid, trec, &
                               l_zonal_mean, var_name, rvalue, L_norec)

USE nlsizes_namelist_mod,  ONLY: global_row_length, global_rows

IMPLICIT NONE

! Subroutine arguments
INTEGER,           INTENT (IN)    :: row_length
INTEGER,           INTENT (IN)    :: rows
INTEGER,           INTENT (IN)    :: fileid     ! file id
INTEGER,           INTENT (IN)    :: trec       ! time record wanted
LOGICAL,           INTENT (IN)    :: l_zonal_mean ! file contains zonal means

CHARACTER (LEN=*), INTENT (IN)    :: var_name   ! variable name

REAL,              INTENT (INOUT) :: rvalue (:,:,:,:) ! variable values to read

LOGICAL,           INTENT (OUT)   :: L_norec    ! T if rec not in file, EOF

! Local variables

INTEGER (KIND=integer32) :: nstatus ! netCDF return code

!  Vector specifying the index in the variable from which
!  the first data value will be read along each dimension
INTEGER (KIND=integer32) :: start   (max_dims)

!  Vector specifying the number of counters to read along each dimension
INTEGER (KIND=integer32) :: counter (max_dims)

INTEGER (KIND=integer32) :: varid32
INTEGER                  :: varid, ndim, dim_size(max_dims)
LOGICAL                  :: lreco    ! true if contains time record dim

CHARACTER (LEN=varname_len) :: varname_local

! Variable to hold file data(lat, lev, waveband)
REAL, ALLOCATABLE        :: global_var_3d(:,:,:)

! Variable to hold expanded data (lon, lat, lev, waveband)
REAL, ALLOCATABLE        :: global_var_4d(:,:,:,:)


INTEGER :: i, k, l      ! loop counters
INTEGER:: icode
CHARACTER (LEN=errormessagelength) :: cmessage  ! error message

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EM_GET_DATA_REAL3D_ZONAL'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! First just check if this routine should have been called
IF (.NOT. l_zonal_mean) THEN
  WRITE (umMessage,'(A,A,A)') 'EM_GET_DATA_REAL3D_ZONAL called for ' // &
                             TRIM(var_name) // &
                            ' but with l_zonal_mean as .false. '
  CALL umPrint(umMessage,src='em_get_data_real3d_zonal')
  ecode = fileid
  CALL ereport ('EM_GET_DATA_REAL3D_ZONAL', ecode,                      &
      'l_zonal_mean was .false. so EM_GET_DATA_REAL3D_ZONAL: ' //       &
      ' should not have been called for ' // TRIM(var_name))
END IF

! Get information on variable. One of the two INOUT arguments of
! EM_GET_VAR_INFO is supplied (varname_local); assign an arbitrary
! negative value to the other one since it is not known yet (varid).
varname_local = var_name    ! need to pass INOUT argument
varid         = -1
CALL em_get_var_info (fileid,      varid,             varname_local, &
                      ndims=ndim,  vdimsize=dim_size, lrec=lreco)


! Allocate the variables that will hold the data.
! The zonal mean data will be copied from global_var_3d into global_var_4d.
! PE0 will fill the array and broadcast to other PEs.
ALLOCATE ( global_var_3d (dim_size(1), dim_size(2), dim_size(3)) )
ALLOCATE ( global_var_4d (global_row_length, &
     dim_size(1), dim_size(2), dim_size(3)) )

! Get values of variable from PE0 and store in global_var_3d
IF ( mype == 0 ) THEN
  ! Default indices/domain to extract.
  ! Default is to get complete data array from NetCDF file so set
  ! all dimensions equal to the size of the variable in the file.
  start   (1:ndim) = 1
  counter (1:ndim) = dim_size(1:ndim)

  ! If var contains time dimension then select record to read
  IF (lreco) THEN
    L_norec = .TRUE.

    ! Record out of bounds and return that no record found
    IF ( trec > dim_size(ndim) ) THEN
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, &
                              zhook_out, zhook_handle)
      RETURN
    END IF
    L_norec        = .FALSE.  ! Indicate to main routine that time rec found

    ! Read a single time from the file
    start   (ndim) = trec     ! last index corresponds to record number
    counter (ndim) = 1        ! read only one record
  END IF

  ! Check that data in NetCDF file matches Global Dimensions
  IF (dim_size (1)   /= global_rows  .OR.                               &
      SIZE(rvalue,3) /= dim_size(2)  .OR.                               &
      SIZE(rvalue,4) /= dim_size(3)) THEN

    ecode = 1
    CALL ereport ('EM_GET_DATA_REAL3D_ZONAL', ecode,                      &
         'Data does not match global dimensions for: ' // TRIM(var_name))
  END IF

  !  Get values of the variable as global_var_3d, specifying what to read
  !  with the optional start and counter arguments
  nfid    = fileid
  varid32 = varid    ! 32-bit integer needed
  nstatus = nf90_get_var (nfid, varid32, global_var_3d, start (1:ndim),      &
                          counter (1:ndim))
  CALL nd_error (nstatus, 'EM_GET_DATA_REAL3D_ZONAL', var_name)
END IF    ! I/O PE

! Expand zonal mean out to all longitudes
DO i = 1, global_row_length
   global_var_4d(i,:,:,:) = global_var_3d(:,:,:)
END DO

icode = 0
CALL gc_gsync(nproc,icode)   ! Synchronise PEs before scatter

! Distribute data to all other PEs
DO l = 1, dim_size(3)
  DO k = 1, dim_size(2)
    ! DEPENDS ON: scatter_field
    CALL scatter_field (rvalue(:,:,k,l), global_var_4d(:,:,k,l), &
                        row_length, rows,                        &
                        global_row_length, global_rows,          &
                        fld_type_p, halo_type_no_halo,           &
                        0, gc_all_proc_group)
  END DO  ! k
END DO  ! l

DEALLOCATE(global_var_4d)
DEALLOCATE(global_var_3d)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE em_get_data_real3D_zonal

!-------------------------------------------------------------------------
! Description:
!   Read data values of four-dimensional real variable in NetCDF file.
!   Part of the interface block EM_GET_DATA used for overloading; this
!   subroutine is called when the variable is of type real 4-dim.
!
! Method:
!   Get variable ID (varid), number of dimensions in variable (ndim),
!   and vector of dimension sizes (dim_size) by calling EM_GET_VAR_INFO.
!   Then get domain/indices to extract and call NF90_GET_VAR to get the
!   variable values, which are read by PEO as global_var. The values
!   are broadcast to all PEs, looping over levels, and returned as rvalue.
!   Assumes (lon, lat, lev, dimension_4).
!-------------------------------------------------------------------------
SUBROUTINE em_get_data_real4D (row_length, rows, fileid, trec,   &
                               var_name, rvalue, L_norec)

USE nlsizes_namelist_mod,  ONLY:global_row_length, global_rows

IMPLICIT NONE

! Subroutine arguments
INTEGER,           INTENT (IN)    :: row_length ! model horiz dimensions
INTEGER,           INTENT (IN)    :: rows
INTEGER,           INTENT (IN)    :: fileid    ! file id
INTEGER,           INTENT (IN)    :: trec      ! time record wanted

CHARACTER (LEN=*), INTENT (IN)    :: var_name  ! variable name

! variable values to read (lon, lat, lev, dimension_4)
REAL,              INTENT (INOUT) :: rvalue (:,:,:,:)

LOGICAL,           INTENT (OUT)   :: L_norec     ! T if rec not in file, EOF

! Local variables

INTEGER (KIND=integer32) :: nstatus ! netCDF return code

! Vector specifying the index in the variable from which
! the first data value will be read along each dimension.
INTEGER (KIND=integer32) :: start   (max_dims)

! Vector specifying the number of counters to read along each dimension
INTEGER (KIND=integer32) :: counter (max_dims)

INTEGER (KIND=integer32) :: varid32
INTEGER                  :: varid, ndim, dim_size(max_dims)
LOGICAL                  :: lreco    ! true if contains time record dim

CHARACTER (LEN=varname_len) :: varname_local

! Variable to hold file data
REAL, ALLOCATABLE        :: global_var(:,:,:,:)

INTEGER:: icode, k, i
CHARACTER (LEN=errormessagelength) :: cmessage  ! error message

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EM_GET_DATA_REAL4D'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Get information on variable (need to do this from all PEs so that
! global_var can be allocated with the right 4th dimension)
!
! One of the two INOUT arguments of EM_GET_VAR_INFO is supplied
! (varname_local); assign an arbitrary negative value to the other
! one since it is not known yet (varid).
varname_local = var_name    ! need to pass INOUT argument
varid         = -1
CALL em_get_var_info (fileid,      varid,             varname_local,   &
                      ndims=ndim,  vdimsize=dim_size, lrec=lreco)

! Check that data in NetCDF file matches Global Dimensions and
! output variable rvalue has correct space for levels
IF (dim_size (1)   /= global_row_length .OR.            &
    dim_size (2)   /= global_rows       .OR.            &
    SIZE(rvalue,3) /= dim_size(3)       .OR.            &
    SIZE(rvalue,4) /= dim_size(4))      THEN
    ecode = 1
    CALL ereport ('EM_GET_DATA_REAL3D', ecode,          &
       'Dimension mismatch for: ' // TRIM(var_name))
END IF

! Allocate the variable that will hold the data. It will be filled
! by PE0 and broadcast to other PEs.
ALLOCATE ( global_var (dim_size(1), dim_size(2),        &
                       dim_size(3), dim_size(4)) )

! Get values of variable from PE0 and store in global_var
IF ( mype == 0 ) THEN
  ! Default indices for the first data value and number of counters
  ! along two dimensions: Global lon and lat. Note that counter(3) will
  ! often be different from the number of model vertical levels, as
  ! it happens in the case of surface emissions.
  start   (1:4) = 1
  counter (1)   = dim_size (1)
  counter (2)   = dim_size (2)
  counter (3)   = dim_size (3)
  counter (4)   = dim_size (4)

  ! If var contains time dimension then select record to read
  IF (lreco) THEN
    L_norec = .TRUE.

    ! Record out of bounds and return that no record found
    IF ( trec > dim_size(ndim) ) THEN
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,  &
                              zhook_out, zhook_handle)
      RETURN
    END IF
    L_norec        = .FALSE.  ! Indicate to main routine that time rec found

    ! Read a single time from the file
    start   (ndim) = trec     ! last index corresponds to record number
    counter (ndim) = 1        ! read only one record
  END IF

  ! Get values of the variable as global_var, specifying what to read
  ! with the optional start and counter arguments
  nfid    = fileid
  varid32 = varid    ! 32-bit integer needed
  nstatus = nf90_get_var (nfid, varid32, global_var, start  (1:ndim),  &
                                                     counter(1:ndim))
  CALL nd_error (nstatus, 'EM_GET_DATA_REAL4D', var_name)
END IF    ! I/O PE

icode = 0
CALL gc_gsync(nproc,icode)   ! Synchronise PEs before scatter

! Distribute data to all PEs, looping over levels
DO i = 1, dim_size(4)
  DO k = 1, dim_size(3)
    ! DEPENDS ON: scatter_field
    CALL scatter_field(rvalue(:,:,k,i), global_var(:,:,k,i),   &
                     row_length, rows,                         &
                     global_row_length, global_rows,           &
                     fld_type_p, halo_type_no_halo,            &
                     0, gc_all_proc_group)
  END DO ! k
END DO ! i

DEALLOCATE(global_var)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE em_get_data_real4D

!-------------------------------------------------------------------------
! Description:
!   Read data values of one-dimensional time axis variable in NetCDF file.
!   Similar to em_get_data_real, except it only works for varid_t, and reads the
!   whole axis rather than a single time record.
!
! Method:
!   First get size of time dimension and allocate axis variable on all PEs
!   Then read data into array and broadcast
!   Result is stored in global variable time_axis
!-------------------------------------------------------------------------
SUBROUTINE em_get_time_axis (fileid)

IMPLICIT NONE

! Subroutine arguments
INTEGER,           INTENT (IN)    :: fileid     ! file id

! Local variables

INTEGER (KIND=integer32) :: nstatus ! netCDF return code

INTEGER (KIND=integer32) :: varid32
INTEGER                  :: ndim, dim_size(max_dims), ilen
CHARACTER(LEN=varname_len) :: var_name

! For Parallel I/O
INTEGER :: my_comm
INTEGER :: icode

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EM_GET_TIME_AXIS'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Obtain the communication group for this run
CALL gc_get_communicator (my_comm, icode)

! Get information on variable, using the global varid_t to specify it.
var_name = ""
CALL em_get_var_info (fileid, varid_t, var_name, ndims=ndim, vdimsize=dim_size)

! Check that the variable has only a single dimension and allocate time_axis to
! length of this dimension. time_axis is deallocated in em_fclose, but check
! this, in case em_fclose wasn't called for some reason.
IF (ndim /= 1) THEN
  ecode = ndim
  CALL ereport (RoutineName, ecode,   &
      'Time axis variable ' //TRIM(var_name)// ' has more than one dimension')
END IF
ilen = dim_size(1)
IF (ALLOCATED(time_axis)) DEALLOCATE(time_axis)
ALLOCATE(time_axis(ilen))

! Only PE0 reads data values from the file
IF ( mype == 0 ) THEN
  !  Get values of the variable as rvalue, specifying what to read
  !  with the optional start and counter arguments
  nfid    = fileid
  varid32 = varid_t    ! 32-bit integer needed
  nstatus = nf90_get_var (nfid, varid32, time_axis)
  CALL nd_error (nstatus, RoutineName, var_name)
END IF   ! I/O PE

! Distribute data to all other PEs
icode = 0
CALL mpl_bcast (time_axis, ilen, mpl_real, 0, my_comm, icode)
IF ( icode /= 0 ) THEN
  ecode = fileid
  CALL ereport (RoutineName, ecode, 'Broadcast Error')
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE em_get_time_axis

! --------------------------------------------------------------------------
! Description:
!   Read attributes of type character from a NetCDF file, given the
!   variable and attribute names.
!   Part of the interface block EM_GET_VAR_ATT used for overloading; this
!   subroutine is called when the required attribute is of type character.
!
! Method:
!   If a global attribute is required, variable name = 'global'.
!   If variable name is passed call NF90_INQ_VARID to get variable-id.
!   Call NF90_GET_ATT to get attribute value: cvalue.
!   Attribute values are read on PE0 and broadcast to all PEs.
! --------------------------------------------------------------------------
SUBROUTINE em_get_var_att_char (fileid, ilen, var_name, att_name, cvalue , &
                                att_exists)

IMPLICIT NONE

!  Subroutine arguments
INTEGER,          INTENT(IN)    :: fileid    ! file id
INTEGER,          INTENT(IN)    :: ilen      ! num items in value
                                             ! (not used for chars)

CHARACTER(LEN=*), INTENT(IN)    :: var_name  ! variable name
CHARACTER(LEN=*), INTENT(IN)    :: att_name  ! attribute name

CHARACTER(LEN=*), INTENT(OUT)   :: cvalue    ! value out

! Optional argument used to check if an attribute is present. If
! attribute found then l_exist will be set to True, otherwise
! l_exist will be False but the model will not stop with error.
LOGICAL, OPTIONAL, INTENT(INOUT) :: att_exists

!  Local variables
INTEGER   (KIND=integer32) :: varid32
INTEGER   (KIND=integer32) :: nstatus          ! netCDF return code
INTEGER,  PARAMETER        :: charlength = 256
CHARACTER (LEN=charlength) :: cval             ! Dummy char to hold value
LOGICAL                    :: l_exist

! For Parallel I/O
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: icode

INTEGER, PARAMETER :: no_of_types = 2          ! logical + char
INTEGER, PARAMETER :: n_log   = 1              ! l_exist
INTEGER, PARAMETER :: n_chars = 1 * charlength ! att_value

! Derived type used to hold attribute values from a NetCDF file
! (read on PEO and broadcast to all PEs).
TYPE nc_fileinfo
  SEQUENCE
  LOGICAL                   :: l_exist
  CHARACTER(LEN=charlength) :: att_value
END TYPE nc_fileinfo

TYPE (nc_fileinfo) :: file_info

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EM_GET_VAR_ATT_CHAR'

!  End of header
   
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Obtain the communication group for this run
CALL gc_get_communicator(my_comm, icode)

! Obtain an object representing the multi-datatype structure which is
! going to be broadcast (see further details above)
CALL setup_nml_type(no_of_types, mpl_nml_type, n_log_in = n_log, &
                    n_chars_in = n_chars)

! PE0 reads attribute values and adds them to the derived type file_info
IF ( mype == 0 ) THEN
  nfid = fileid
  ! Get variable id
  IF ( TRIM(var_name) == 'GLOBAL' .OR. TRIM(var_name) == 'global' .OR. &
       TRIM(var_name) == 'Global' ) THEN
    varid32 = NF90_GLOBAL
  ELSE
    nstatus = nf90_inq_varid (nfid, var_name, varid32)
    CALL nd_error (nstatus, 'EM_GET_VAR_ATT_CHAR: INQ Varid', var_name)
  END IF

  ! Get value of attribute
  cval    = ' '
  nstatus = nf90_get_att (nfid, varid32, att_name, cval)

  ! If argument att_exists is present that means that the attribute is not
  ! compulsory. Therefore, if attribute found in the previous call to 
  ! nf90_get_att then set the local variable l_exist to True; otherwise 
  ! leave it as False but there is no need to flag error.
  l_exist = .FALSE.
  IF ( PRESENT(att_exists) ) THEN
    IF (nstatus == nf90_noerr) l_exist = .TRUE. ! attrib exists

  ! If argument att_exists not present then the attribute is needed. Therefore
  ! call nd_error, which will stop the model only if the previous call to
  ! nf90_get_att reported error (as given by nstatus).
  ELSE
    CALL nd_error (nstatus, 'EM_GET_VAR_ATT_CHAR: Get Attribute for ' // &
                   TRIM(var_name), att_name)    ! Check for error
  END IF

  ! Fill in output argument.
  cvalue = TRIM(cval)

  ! Copy into broadcast type
  file_info % att_value = cvalue
  file_info % l_exist   = l_exist

END IF

! Distribute data to all other PEs
icode = 0
CALL mpl_bcast (file_info, 1, mpl_nml_type, 0, my_comm, icode)
IF ( icode /= 0 ) THEN
  ecode = fileid
  CALL ereport ('EM_GET_VAR_ATT_CHAR', ecode, 'NML Broadcast Error')
END IF

! Operation on all PEs different from PE0: Transfer data from broadcast type
! into output argument cvalue and local variable l_exist.
IF ( mype /= 0 ) THEN
  cvalue  = file_info % att_value
  l_exist = file_info % l_exist
END IF

! Fill in optional argument (if present) with value from local variable.
IF ( PRESENT(att_exists) ) att_exists = l_exist

! Reset the multi-datatype structure allocated for broadcasting
CALL mpl_type_free (mpl_nml_type, icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE em_get_var_att_char

! --------------------------------------------------------------------------
! Description:
!   Read attributes of type integer from a NetCDF file, given the
!   variable and attribute names.
!   Part of the interface block EM_GET_VAR_ATT used for overloading; this
!   subroutine is called when the required attribute is of type integer.
!   It can also convert character attributes to integers for compatibility 
!   with older files which contain integer attributes in strings, e.g. "1".
!
! Method:
!   If a global attribute is required, variable name = 'global'.
!   If variable name is passed call NF90_INQ_VARID to get variable-id.
!   Call NF90_GET_ATT to get attribute value: ivalue.
!   Attribute values are read on PE0 and broadcast to all PEs.
! --------------------------------------------------------------------------
SUBROUTINE em_get_var_att_int (fileid, ilen, var_name, att_name, ivalue, &
                               att_exists )

IMPLICIT NONE

!  Subroutine arguments
INTEGER,          INTENT(IN)    :: fileid       ! file id
INTEGER,          INTENT(IN)    :: ilen         ! num items in value (just
                                                ! in case it is an array)

CHARACTER(LEN=*), INTENT(IN)    :: var_name     ! variable name
CHARACTER(LEN=*), INTENT(IN)    :: att_name     ! attribute name

INTEGER,          INTENT(OUT)   :: ivalue(ilen) ! value out

! Optional argument used to check if an attribute is present. If
! attribute found then l_exist will be set to True, otherwise
! l_exist will be False but the model will not stop with error.
LOGICAL, OPTIONAL, INTENT(INOUT) :: att_exists

!  Local variables
LOGICAL :: l_exist
INTEGER,  PARAMETER        :: charlength = 12
CHARACTER (LEN=charlength) :: cval             ! Dummy char to hold value

! 32-bit integers to actually pass to NC routines
INTEGER(KIND=integer32) :: varid32, ival32(ilen)
INTEGER(KIND=integer32) :: nstatus               ! netCDF return code
INTEGER(KIND=integer32) :: iattlen               ! length of attribute in file
INTEGER(KIND=integer32) :: xtype                 ! netcdf data type

! For Parallel I/O
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: icode

INTEGER, PARAMETER :: no_of_types = 2     ! integer + logical
INTEGER, PARAMETER :: n_int = max_att_len ! max possible size
INTEGER, PARAMETER :: n_log = 1           ! l_exist

! Derived type used to hold attribute values from a NetCDF file
! (read on PEO and broadcast to all PEs).
TYPE nc_fileinfo
  SEQUENCE
  INTEGER :: att_value(max_att_len)
  LOGICAL :: l_exist
END TYPE nc_fileinfo

TYPE (nc_fileinfo) :: file_info

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EM_GET_VAR_ATT_INT'

!  End of header
   
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Obtain the communication group for this run
CALL gc_get_communicator(my_comm, icode)

! Obtain an object representing the multi-datatype structure which is
! going to be broadcast (see further details above)
CALL setup_nml_type (no_of_types, mpl_nml_type, n_int_in = n_int, &
                     n_log_in = n_log)

! PE0 reads attribute values and adds them to the derived type file_info
IF ( mype == 0 ) THEN
  ! Get variable id
  nfid = fileid
  IF ( TRIM(var_name) == 'GLOBAL'   .OR. TRIM(var_name) == 'global' .OR. &
       TRIM(var_name) == 'Global' ) THEN
    varid32 = NF90_GLOBAL
  ELSE
    nstatus = nf90_inq_varid (nfid, var_name, varid32)
    CALL nd_error (nstatus, 'EM_GET_VAR_ATT_INT: INQ Varid', var_name)
  END IF

  ! Check data type of attribute and whether it exists
  nstatus = nf90_inquire_attribute(nfid, varid32, att_name, xtype=xtype,  &
                                   len=iattlen)

  ! Set l_exist according to whether this call returned an error (nstatus).
  ! If optional argument att_exists is present, then attribute is allowed to be
  ! missing and the value of l_exist is returned. If not, call nd_error which
  ! will stop the model since nstatus /= nf90_noerr.
  IF (nstatus == nf90_noerr) THEN
    l_exist = .TRUE.  ! attrib exists
  ELSE
    l_exist = .FALSE.
    IF (.NOT. PRESENT(att_exists) ) THEN
      CALL nd_error (nstatus, 'EM_GET_VAR_ATT_INT: Inquire Attribute for ' // &
                     TRIM(var_name), att_name)
    END IF
  END IF

  ! If attribute does exist, read it. Convert data type here if required.
  ival32(:) = 0
  IF (l_exist) THEN
    SELECT CASE(xtype)
    CASE(NF90_CHAR)
      ! Treat string as a single integer, iattlen is length of string.
      ! '(I12)' format statement must agree with charlength variable. 12 digits
      ! is enough to hold any 32-bit integer.
      cval = ' '
      nstatus = nf90_get_att (nfid, varid32, att_name, cval(1:iattlen))
      READ (cval, '(I12)') ival32(1)
      
    CASE DEFAULT
      ! All other types are numbers so netcdf will convert to int on the fly
      nstatus   = nf90_get_att (nfid, varid32, att_name, ival32(1:iattlen))

    END SELECT
    CALL nd_error (nstatus, 'EM_GET_VAR_ATT_INT: Get Attribute for ' // &
                   TRIM(var_name), att_name)
  END IF

  ! Store local variables on broadcast type
  file_info % l_exist           = l_exist
  file_info % att_value(1:ilen) = ival32(1:ilen) 

END IF

! Distribute data to all other PEs
icode = 0
CALL mpl_bcast (file_info, 1, mpl_nml_type, 0, my_comm, icode)
IF ( icode /= 0 ) THEN
  ecode = fileid
  CALL ereport ('EM_GET_ATT_INT', ecode, 'NML Broadcast Error')
END IF
  
! Operation on all PEs different from PE0: Transfer data from
! broadcast type into local variables.
IF ( mype /= 0 ) THEN
  ival32(1:ilen) = file_info % att_value(1:ilen)
  l_exist        = file_info % l_exist
END IF

! Fill in the output argument 'ivalue' and (if present) the optional
! argument 'att_exists' with the values taken from the local variables.
! Operation done on all processors.
ivalue = ival32
IF ( PRESENT(att_exists) ) att_exists = l_exist

! Reset the multi-datatype structure allocated for broadcasting
CALL mpl_type_free (mpl_nml_type, icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE em_get_var_att_int

! --------------------------------------------------------------------------
! Description:
!   Read attributes of type real from a NetCDF file, given the
!   variable and attribute names.
!   Part of the interface block EM_GET_VAR_ATT used for overloading; this
!   subroutine is called when the required attribute is of type real.
!
! Method:
!   If a global attribute is required, variable name = 'global'.
!   If variable name is passed call NF90_INQ_VARID to get variable-id.
!   Call NF90_GET_ATT to get attribute value: rvalue.
!   Attribute values are read on PE0 and broadcast to all PEs.
! --------------------------------------------------------------------------
SUBROUTINE em_get_var_att_real (fileid, ilen, var_name, att_name, rvalue, &
                                att_exists)

IMPLICIT NONE

!  Subroutine arguments
INTEGER,           INTENT(IN)    :: fileid       ! file id
INTEGER,           INTENT(IN)    :: ilen         ! num items in value (just in
                                                 ! case it is an array)

CHARACTER(LEN=*),  INTENT(IN)    :: var_name     ! variable name
CHARACTER(LEN=*),  INTENT(IN)    :: att_name     ! attribute name

REAL,              INTENT(OUT)   :: rvalue(ilen) ! value out

! Optional argument used to check if an attribute is present. If
! attribute found then l_exist will be set to True, otherwise
! l_exist will be False but the model will not stop with error.
LOGICAL, OPTIONAL, INTENT(INOUT) :: att_exists

!  Local variables
INTEGER(KIND=integer32) :: varid32
INTEGER(KIND=integer32) :: nstatus       ! netCDF return code
REAL   (KIND=real32)    :: rval32(ilen)
LOGICAL                 :: l_exist 

! For Parallel I/O
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: icode

INTEGER, PARAMETER :: no_of_types = 2        ! real + logical
INTEGER, PARAMETER :: n_real = max_att_len   ! max possible size
INTEGER, PARAMETER :: n_log = 1              ! l_exist

! Derived type used to hold attribute values from a NetCDF file
! (read on PEO and broadcast to all PEs).
TYPE nc_fileinfo
  SEQUENCE
  REAL    :: att_value(max_att_len)
  LOGICAL :: l_exist
END TYPE nc_fileinfo

TYPE (nc_fileinfo) :: file_info

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EM_GET_VAR_ATT_REAL'

!  End of header
   
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Obtain the communication group for this run
CALL gc_get_communicator(my_comm, icode)

! Obtain an object representing the multi-datatype structure which is
! going to be broadcast (see futher details above)
CALL setup_nml_type(no_of_types, mpl_nml_type, n_real_in = n_real, &
                    n_log_in = n_log)

! PE0 reads attribute values and adds them to the derived type file_info
IF ( mype == 0 ) THEN
  ! Get variable id
  nfid = fileid
  IF ( TRIM(var_name) == 'GLOBAL'   .OR.  TRIM(var_name) == 'global' .OR. &
       TRIM(var_name) == 'Global' ) THEN
    varid32 = NF90_GLOBAL
  ELSE
    nstatus = nf90_inq_varid (nfid, var_name, varid32)
    CALL nd_error (nstatus, 'EM_GET_VAR_ATT_REAL: INQ Varid', var_name)
  END IF

  ! Get values of attribute
  rval32(:) = 0.0
  nstatus   = nf90_get_att (nfid, varid32, att_name, rval32)

  ! If argument att_exists is present that means that the attribute is not
  ! compulsory. Therefore, if attribute found in the previous call to 
  ! nf90_get_att then change the local variable l_exist to True; otherwise
  ! leave it as False but there is no need to flag error.
  l_exist = .FALSE.
  IF ( PRESENT(att_exists)  ) THEN
    IF (nstatus == nf90_noerr) l_exist = .TRUE.  ! attrib exists

  ! If argument att_exists not present then the attribute is needed. Therefore
  ! call nd_error, which will stop the model only if the previous call to
  ! nf90_get_att reported error (as given by nstatus).
  ELSE
    CALL nd_error (nstatus, 'EM_GET_VAR_ATT_REAL: Get Attribute for ' //  &
                   TRIM(var_name), att_name)
  END IF

  ! Copy local variables into broadcast type
  file_info % l_exist           = l_exist
  file_info % att_value(1:ilen) = rval32(1:ilen)

END IF

! Distribute data to all other PEs
icode = 0
CALL mpl_bcast (file_info, 1, mpl_nml_type, 0, my_comm, icode)
IF ( icode /= 0 ) THEN
  ecode = fileid
  CALL ereport ('EM_GET_VAR_ATT_REAL', ecode, 'NML Broadcast Error')
END IF

! Operation on all PEs different from PE0: Transfer data from
! broadcast type to local variables.
IF ( mype /= 0 ) THEN
  l_exist        = file_info % l_exist
  rval32(1:ilen) = file_info % att_value(1:ilen)
END IF

! Operation on all PEs: Transfer data from local variables to
! OUT and (if present) optional INOUT arguments.
rvalue = rval32
IF ( PRESENT(att_exists) ) att_exists = l_exist

! Reset the multi-datatype structure allocated for broadcasting
CALL mpl_type_free (mpl_nml_type, icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE em_get_var_att_real


! ------------------------------------------------------------------------
! Description:
!   Close NetCDF emission file.
!
! Method:
!   Call NF90_CLOSE to close a NetCDF file and return nstatus. Also call
!   ND_ERROR, which will stop the model if nstatus reports errors. Note
!   that this is done only by PE0.
!   If needed reset time dimension related values since they are specific
!   to each file.
! ------------------------------------------------------------------------
SUBROUTINE em_fclose (fileid)

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: fileid

! Local variables
INTEGER (KIND=integer32)       :: nstatus        ! netCDF return code

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EM_FCLOSE'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Close file only by PE0
IF ( mype == 0 ) THEN
  nfid    = fileid
  nstatus = nf90_close(nfid)
  CALL nd_error (nstatus, 'EM_FCLOSE', ' ')
END IF

! Reset time dimension related values (specific to each file)
! if they are populated
IF ( varid_t > 0 .OR. dimid_t > 0 ) THEN   ! Default -99
  varid_t  = -99
  dimid_t  = -99
  tim_name = 'empty'
  ct_units = 'empty'
END IF

IF (ALLOCATED(time_axis)) DEALLOCATE(time_axis)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE em_fclose

! -------------------------------------------------------------------------
! Description:
!   Checks return status of netcdf call and stop if needed.
!
! Method:
!   If the return status (nstatus) indicates error then:
!    * Close the NetCDF file by using NF90_CLOSE.
!    * Obtain error message from NF90_STRERROR and concatenate it with
!      the string cmessage, which contains more detailed information passed
!      from the routine where the error was found. Then print out
!      those messages and call EREPORT to stop the model.
! -------------------------------------------------------------------------
SUBROUTINE nd_error (nstatus, routine_name, cmessage)

IMPLICIT NONE

! Subroutine arguments
INTEGER (KIND=integer32), INTENT(IN) :: nstatus  ! netCDF return code

CHARACTER (LEN=*), INTENT(IN) :: routine_name ! routine where error found
CHARACTER (LEN=*), INTENT(IN) :: cmessage     ! error messg from that routine

! Local variables
INTEGER (KIND=integer32)       :: ncode

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ND_ERROR'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF ( nstatus /= nf90_noerr ) THEN
  ! Close file only by PE0
  IF ( mype == 0 ) THEN
    ncode = nf90_close(nfid)
  END IF
  ecode = 1

  ! Printing error messages first to avoid passing
  ! a very long string argument to ereport
  WRITE (umMessage,'(a)') 'Error in ' // routine_name    //  ": "  // &
                  TRIM (nf90_strerror(nstatus))  //  " - " // &
                  TRIM (cmessage)
  CALL umPrint(umMessage,src='nd_error')
  CALL ereport (routine_name, ecode, 'NetCDF error')
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN

END SUBROUTINE nd_error
! -------------------------------------------------------

END MODULE emiss_io_mod
