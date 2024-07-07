! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Module to hold data and procedures for handling NetCDF climatology files.
!
!  Method:
!    Contains a number of subroutines which call some Fortran NetCDF
!    procedures to get information/values of dimensions/variables
!    present in NetCDF climatology files. It is generally assumed that
!    dimensions of variables are ordered x, y, z, t as recommended by CF
!    conventions (though any dimension can be omitted).
!    For further information see the header of each subroutine
!    contained here.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: GLOMAP_CLIM
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE glomap_clim_netcdf_io_mod

USE conversions_mod,          ONLY: isec_per_day
USE ereport_mod,              ONLY: ereport
USE errormessagelength_mod,   ONLY: errormessagelength

USE glomap_clim_netcdf_parameter_mod, ONLY: &
    max_att_len,                            &
    max_dims_per_var,                       &
    t_calendar_len,                         &
    t_units_len,                            &
    vartype_len

USE items_nml_mod,            ONLY: netcdf_varname_len

USE mpl,                      ONLY: &
    mpl_character,                  &
    mpl_integer,                    &
    mpl_real,                       &
    mpl_logical

USE netcdf,                   ONLY: &
    nf90_char,                      &
    nf90_close,                     &
    nf90_double,                    &
    nf90_ebadname,                  &
    nf90_enotvar,                   &
    nf90_float,                     &
    nf90_get_att,                   &
    nf90_get_var,                   &
    nf90_global,                    &
    nf90_inq_varid,                 &
    nf90_inquire,                   &
    nf90_inquire_dimension,         &
    nf90_inquire_variable,          &
    nf90_int,                       &
    nf90_noerr,                     &
    nf90_nowrite,                   &
    nf90_open,                      &
    nf90_strerror

USE nlstcall_mod,             ONLY: lcal360

USE parkind1,                 ONLY: &
    jpim,                           &
    jprb

USE setup_namelist,           ONLY: setup_nml_type
 
USE UM_Parvars,               ONLY: &
    gc_all_proc_group,              &
    current_decomp_type,            &
    change_decomposition

USE UM_ParCore,               ONLY: &
    mype,                           &
    nproc

USE um_types,                 ONLY: &
    integer32,                      &
    real32

USE umPrintMgr,               ONLY: &
    umPrint,                        &
    umMessage,                      &
    PrStatus_Diag,                  &
    PrintStatus

USE yomhook,                  ONLY: &
    lhook,                          &
    dr_hook

IMPLICIT NONE

! Default private
PRIVATE

! Local variables
INTEGER (KIND=integer32), SAVE :: ncdfid ! ID of NetCDF file
INTEGER                        :: ecode  ! return code for ereport

! Length of character strings used for padding
INTEGER, PARAMETER         :: padding_a_len  = 36
INTEGER, PARAMETER         :: padding_b_len  = 188

! Numbers zero and one
INTEGER, PARAMETER         :: one  = 1
INTEGER, PARAMETER         :: zero = 0

! Declare integer parameters indicating the possible ways of updating
! climatology, following the same conventions as for ancillary files
INTEGER, PARAMETER, PUBLIC :: iupdn_single = 0  ! single time (not tested yet)
INTEGER, PARAMETER, PUBLIC :: iupdn_serial = 1  ! time series
INTEGER, PARAMETER, PUBLIC :: iupdn_cyclic = 2  ! perodic time series

! Arbitrary year for cyclic files
INTEGER, PARAMETER, PUBLIC :: cyclic_year_netcdf = 1

! Variables to hold information on the record ('time') dimension.
! These will be populated (unless not requested) during netcdf_fopen
! and reset (if populated) with the netcdf_fclose call
CHARACTER(LEN=netcdf_varname_len), SAVE :: tim_name    = 'unset'
CHARACTER(LEN=t_units_len),        SAVE :: ct_units    = 'unset'
CHARACTER(LEN=t_calendar_len),     SAVE :: ct_calendar = 'unset'
INTEGER,                           SAVE :: varid_t     = -99
INTEGER,                           SAVE :: dimid_t     = -99

! Subroutines and overloaded interfaces available outside this module
PUBLIC :: netcdf_fopen
PUBLIC :: netcdf_get_file_info
PUBLIC :: netcdf_get_var_info
!         netcdf_get_time_info  is only called within this module
!         netcdf_get_time_var   is only called within this module
PUBLIC :: netcdf_var_check_dims
PUBLIC :: netcdf_get_time_rec
PUBLIC :: netcdf_get_data
PUBLIC :: netcdf_get_var_att
PUBLIC :: netcdf_fclose
!         netcdf_error          is only called within this module

! ------------------------------------------------------------------------
! Interface blocks used for overloading. 
!
! A generic call to NETCDF_GET_DATA will select the appropriate subroutine 
! automatically, depending on the type of the INTENT(INOUT) argument
! (variable values).
!
! A generic call to NETCDF_GET_VAR_ATT will select the appropriate subroutine
! automatically, depending on the type of the INTENT(OUT) argument 
! (attribute value).
!
! Note that both NETCDF_GET_DATA and NETCDF_GET_VAR_ATT are public. Therefore
! the generic interfaces can be called from outside this module, while the
! type specific subroutines can only be used inside this module.
! ------------------------------------------------------------------------

INTERFACE netcdf_get_data
MODULE PROCEDURE             &
  netcdf_get_data_real_1D,   &
  netcdf_get_data_real_field
END INTERFACE

INTERFACE netcdf_get_var_att
MODULE PROCEDURE             &
  netcdf_get_var_att_char,   &
  netcdf_get_var_att_int
END INTERFACE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='GLOMAP_CLIM_NETCDF_IO_MOD'

CONTAINS

! ------------------------------------------------------------------------
! Description:
!   Check for existence of NetCDF climatology file and open it.
!
! Method:
!   Call NF90_OPEN to open a NetCDF file and return nstatus.
!   Call NETCDF_ERROR, which will stop the model if nstatus reports errors.
!   By default, check that a record ('time') dimension exists in the 
!   file, through a call to NETCDF_GET_TIME_INFO. This can be overriden by
!   optional argument ignore_time = True
! ------------------------------------------------------------------------
SUBROUTINE netcdf_fopen (filename, fileid, ignore_time)

IMPLICIT NONE

! Subroutine arguments
CHARACTER (LEN=*), INTENT(IN)  :: filename

INTEGER,           INTENT(OUT) :: fileid

! Optional flag that will be set to True only if it is needed to prevent
! the default attempt to detect the time dimension in a file.
LOGICAL, OPTIONAL, INTENT(IN)  :: ignore_time

! Local variables
INTEGER (KIND=integer32)       :: nstatus       ! netCDF return code

! Local version of "ignore_time" with default setting
LOGICAL :: l_ignore_time = .FALSE.

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='NETCDF_FOPEN'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF ( mype == 0 ) THEN
  nstatus = nf90_open (filename, nf90_nowrite, ncdfid)
  CALL netcdf_error (nstatus, RoutineName, TRIM(filename))
  fileid = ncdfid
ELSE
  !Set file-id to dummy on other PEs, as they would not access the files anyway
  fileid = 99  
END IF

! Copy the value of argument "ignore_time" to the local variable
IF (PRESENT(ignore_time)) l_ignore_time = ignore_time

! By default NETCDF_GET_TIME_INFO is called when opening a file to look
! for the time dimension (1-D, 'days since' or 'hours since') and
! stores some values in global variables. However make sure that
! l_ignore_time is not set to True just in case there are situations
! when the call is not needed.
IF ( .NOT. l_ignore_time ) THEN
  CALL netcdf_get_time_info (filename, fileid)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE netcdf_fopen

! ------------------------------------------------------------------------
! Description:
!   The NF90_INQUIRE subroutine returns information about an open netCDF
!   dataset, given its netCDF ID. The subroutine can be called from either
!   define mode or data mode, and returns values for any or all of the
!   following:
!   the number of dimensions,
!   the number of variables,
!   the number of global attributes,
!   the dimension ID of the dimension defined with unlimited length, if any.
!
! Method:
!   Call NF90_INQUIRE and return only the information required. The data
!    is read on PEO and broadcast to all PEs. 
!
! ------------------------------------------------------------------------
SUBROUTINE netcdf_get_file_info (fileid, nDimensions, nVariables,             &
                                 nAttributes, unlimitedDimId)

IMPLICIT NONE

! Subroutine arguments

INTEGER, INTENT(IN)           :: fileid         ! from previous NF90_OPEN call
INTEGER, OPTIONAL,INTENT(OUT) :: nDimensions    ! no. dimensions of dataset
INTEGER, OPTIONAL,INTENT(OUT) :: nVariables     ! no. variables of dataset
INTEGER, OPTIONAL,INTENT(OUT) :: nAttributes    ! no. attributes of dataset
INTEGER, OPTIONAL,INTENT(OUT) :: unlimitedDimId ! returned ID of the 
                                                ! unlimited dimension
                                                ! unlimitedDimId returns -1
                                                ! if there
                                                ! are no unlimited dimensions

! Local variables
INTEGER (KIND=integer32) :: nvars32, ndims32, natts32, undimid32
INTEGER (KIND=integer32) :: nstatus       ! netCDF return code

! For Parallel I/O
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: icode

INTEGER, PARAMETER :: no_of_types = 1 ! integer
INTEGER, PARAMETER :: n_int = 4

! Derived type used to hold information about a NetCDF file
! (read on PEO and broadcast to all PEs).
TYPE nc_fileinfo
  SEQUENCE
  INTEGER :: nVariables
  INTEGER :: nDimensions
  INTEGER :: nAttributes
  INTEGER :: unlimitedDimId
END TYPE nc_fileinfo

TYPE (nc_fileinfo) :: file_info

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='NETCDF_GET_FILE_INFO'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Obtain the communication group for this run
CALL gc_get_communicator (my_comm, icode)

! Obtain an object representing the multi-datatype structure which is
! going to be broadcast.
CALL setup_nml_type (no_of_types, mpl_nml_type, n_int_in = n_int)

IF ( mype == 0 ) THEN
  ncdfid = fileid
  ! Get all file information 
  nstatus = nf90_inquire(ncdfid, ndims32, nvars32, natts32, undimid32 )
  CALL netcdf_error (nstatus, RoutineName, 'INQUIRE')
  file_info % nVariables = nvars32
  file_info % nDimensions = ndims32
  file_info % nAttributes  = natts32
  file_info % unlimitedDimId = undimid32
END IF

! Data structure broadcast to all processors
icode = 0
CALL mpl_bcast (file_info, one, mpl_nml_type, zero, my_comm, icode)
IF ( icode /= 0 ) THEN 
  ecode = fileid
  CALL ereport (RoutineName, ecode,'NML Broadcast Error')
END IF

! Copy data broadcast from PE0 to PE local variables.
IF ( mype /= 0 ) THEN
  nvars32   = file_info % nVariables
  ndims32   = file_info % nDimensions
  natts32   = file_info % nAttributes
  undimid32 = file_info % unlimitedDimId
END IF

! Check what information is required and transfer data from
! the local variables to the output arguments. Operation done on all PEs.
IF ( PRESENT(nVariables    ) ) nVariables     = nvars32
IF ( PRESENT(nDimensions   ) ) nDimensions    = ndims32
IF ( PRESENT(nAttributes   ) ) nAttributes    = natts32
IF ( PRESENT(unlimitedDimId) ) unlimitedDimId = undimid32

! Reset the multi-datatype structure allocated for broadcasting
CALL mpl_type_free (mpl_nml_type, icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE netcdf_get_file_info

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
SUBROUTINE netcdf_get_var_info (fileid, varid, varname,                       &
                                ndims, vdimsize, lrec, var_type)

IMPLICIT NONE

! Subroutine arguments
INTEGER,           INTENT(IN)                     :: fileid   ! File id
INTEGER,           INTENT(INOUT)                  :: varid    ! Variable id
CHARACTER (LEN=netcdf_varname_len), INTENT(INOUT) :: varname  ! Variable name
INTEGER, OPTIONAL, INTENT(OUT)                    :: ndims    ! Nr of dims
INTEGER, OPTIONAL, INTENT(OUT)  :: vdimsize(max_dims_per_var) ! Size of dims
LOGICAL, OPTIONAL, INTENT(OUT)  :: lrec ! Flag set to True if variable 
                                        ! contains record/time dimension
CHARACTER (LEN=vartype_len), OPTIONAL, INTENT(OUT) :: var_type ! Variable type

! Local variables
INTEGER (KIND=integer32)           :: nstatus
INTEGER (KIND=integer32)           :: varid_local, itype, idumm, n
INTEGER (KIND=integer32)           :: ndim_local
INTEGER (KIND=integer32)           :: dimid       (max_dims_per_var)
INTEGER (KIND=integer32)           :: dsize_local (max_dims_per_var)
LOGICAL                            :: lrec_local
CHARACTER (LEN=vartype_len)        :: var_type_local
CHARACTER (LEN=netcdf_varname_len) :: varname_local
CHARACTER (LEN=padding_a_len)      :: padding_local

! For Parallel I/O
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: icode

INTEGER, PARAMETER :: no_of_types = 3                ! int + logical + char
INTEGER, PARAMETER :: n_int   = 2 + max_dims_per_var 
INTEGER, PARAMETER :: n_log   = 1

! Nr. of characters for varname and var_type (see below)
INTEGER, PARAMETER :: n_chars = netcdf_varname_len + vartype_len +            &
                                padding_a_len

TYPE nc_fileinfo
  SEQUENCE
  INTEGER                           :: varid
  INTEGER                           :: ndims
  INTEGER                           :: vdimsize(max_dims_per_var)
  LOGICAL                           :: lrec
  CHARACTER(LEN=netcdf_varname_len) :: varname
  CHARACTER(LEN=vartype_len)        :: var_type
  CHARACTER(LEN=padding_a_len)      :: padding
END TYPE nc_fileinfo

! Padding is used to ensure total length is a mulitple of longest length

TYPE (nc_fileinfo) :: file_info

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER    :: RoutineName='NETCDF_GET_VAR_INFO'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Obtain the communication group for this run
CALL gc_get_communicator(my_comm, icode)

! Obtain an object representing the multi-datatype structure which is
! going to be broadcast (see further details above)
CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in = n_int,              &
                    n_log_in = n_log, n_chars_in = n_chars )

! Check whether Name or ID of variable in NetCDF file is supplied.
! PE0 reads the relevant data and adds it to the derived type file_info.
IF ( mype == 0 ) THEN
  ncdfid = fileid
  
  ! If varname provided then get varid
  IF ( TRIM(varname) /= "" ) THEN
    ! Note: 'varname' is intent(in) in nf90_inq_varid.
    ! See note below about 'varname_local'.
    nstatus = nf90_inq_varid (ncdfid, TRIM(varname), varid_local)
    CALL netcdf_error ( nstatus, RoutineName,                                 & 
                                 'INQ_VARID ' // TRIM(varname) )
    varid = varid_local
  ELSE
    ! If varname not provided then make sure that varid has
    ! been provided; if that is not the case then report error.
    IF ( varid <= 0 ) THEN
      ecode = 54901
      CALL ereport (RoutineName, ecode,                                       &
                   'Either variable name or ID should be provided')
      END IF
      varid_local = varid
  END IF            ! check if name/id provided
  
  ! Get type information, dimension ids, etc.
  ! 
  ! Here we get the name of the variable, now 'varname_local',
  ! which unlike 'varname' above needs to be intent(out) in
  ! nf90_inquire_variable. That is the reason why two
  ! 'varnames' are used. It is necessary to get 'varname_local'
  ! because this routine can be called with 'varname' unset.
  !
  ! Note also that if the user did pass in 'varname' it will be
  ! overwritten by 'varname_local', but that is OK because 
  ! the code above has ensured that the variable name 
  ! corresponds to the id being inquired about.
  !
  nstatus = nf90_inquire_variable (ncdfid, varid_local, name=varname_local,   &
                                   xtype=itype, ndims=ndim_local,             &
                                   dimids=dimid, nAtts=idumm)
  
  CALL netcdf_error (nstatus, RoutineName,                                    & 
                    'INQ_VARIABLE ' // TRIM(varname_local) )
  varname = varname_local
  
  ! Get size of each dimension
  IF ( PRESENT(vdimsize) ) THEN
    DO n = 1, ndim_local
      nstatus = nf90_inquire_dimension (ncdfid, dimid(n),LEN=dsize_local(n))
      CALL netcdf_error (nstatus, RoutineName, 'INQ_DIMSIZE ' // TRIM(varname))
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
      WRITE (umMessage,'(A)') TRIM(RoutineName) //                            &
              ' : Time dimension info not found, check if netcdf_fopen ' //   &
              'was called with option l_ignore_time = .TRUE.'
      CALL umPrint(umMessage,src=RoutineName)
      nstatus = nf90_enotvar
      CALL netcdf_error (nstatus, RoutineName,                                &
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
      var_type_local = 'INTR'
    CASE (nf90_float)
      var_type_local = 'REAL'
    CASE (nf90_double)
      var_type_local = 'DBLE'
    CASE DEFAULT
      WRITE (umMessage,'(A,1X,I2)') TRIM(RoutineName) //                      &
                              ': Unknown variable type defined ', itype
      CALL umPrint (umMessage, src=RoutineName)
      
      WRITE (umMessage,'(A)') 'Should be one of: CHAR / INTR / REAL / DBLE'
      CALL umPrint (umMessage, src=RoutineName)
      
      nstatus = nf90_ebadname
      CALL netcdf_error ( nstatus, RoutineName, 'VAR_TYPE ' // TRIM(varname) )
    END SELECT
  END IF
  
  padding_local = 'ABCDEFGH'
      
  ! Fill in derived type variable to later do broadcast
  file_info % varid       = varid_local
  file_info % ndims       = ndim_local
  file_info % vdimsize(:) = dsize_local(:)
  file_info % varname     = varname_local
  file_info % lrec        = lrec_local
  file_info % var_type    = var_type_local
  file_info % padding     = padding_local
END IF

! Data structure broadcast to all processors
icode = 0
CALL mpl_bcast (file_info, one, mpl_nml_type, zero, my_comm, icode)
IF ( icode /= 0 ) THEN 
  ecode = fileid
  CALL ereport (RoutineName, ecode, 'NML Broadcast Error')
END IF

! Copy data broadcast from PE0 to PE local variables.
IF ( mype /= 0 ) THEN
  varid_local    = file_info % varid
  varname_local  = file_info % varname
  ndim_local     = file_info % ndims
  dsize_local(:) = file_info % vdimsize(:)
  lrec_local     = file_info % lrec
  var_type_local = file_info % var_type
  padding_local  = file_info % padding
END IF

! Copy the local variables back onto the output variables
varid  = varid_local
varname = varname_local
IF ( PRESENT(ndims)    ) ndims       = ndim_local
IF ( PRESENT(vdimsize) ) vdimsize(:) = dsize_local(:)
IF ( PRESENT(lrec)     ) lrec        = lrec_local
IF ( PRESENT(var_type) ) var_type    = var_type_local

! Reset the multi-datatype structure allocated for broadcasting
CALL mpl_type_free (mpl_nml_type, icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE netcdf_get_var_info

! --------------------------------------------------------------------------
! Description:
!   Check that there is a time dimension for every 4D variable.
!   Assumes that there is a 1D 'time' variable is uni-dimensional which has 
!   a 'units' attribute starting with 'days since' or 'hours since'.
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
SUBROUTINE netcdf_get_time_info (filename, fileid)

IMPLICIT NONE

! Subroutine arguments
CHARACTER (LEN=*), INTENT(IN)  :: filename
INTEGER,           INTENT(IN)  :: fileid

! Local variables
INTEGER                  :: ierror       ! Error code
INTEGER                  :: time_counter ! No. of 1D time variables present
INTEGER                  :: v            ! loop counter for variables
INTEGER (KIND=integer32) :: nstatus      ! netCDF return code
INTEGER (KIND=integer32) :: ndim_local
INTEGER (KIND=integer32) :: dimid (max_dims_per_var)
INTEGER (KIND=integer32) :: itype
INTEGER (KIND=integer32) :: varid32
INTEGER (KIND=integer32) :: nvars32, ndims32, natts32, undimid32
LOGICAL                  :: l_verify_time   ! Additional verification
CHARACTER (LEN=1)        :: caxis
CHARACTER (LEN=errormessagelength) :: cmessage
CHARACTER (LEN=netcdf_varname_len) :: varnam_local  ! Name read from file
CHARACTER (LEN=t_units_len)        :: cunits        ! For 'units' attribute
CHARACTER (LEN=padding_b_len)      :: padding_local

! For Parallel I/O
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: icode

INTEGER, PARAMETER :: no_of_types = 2                     ! int + char
INTEGER, PARAMETER :: n_int       = 2                     ! varid_t + dimid_t
INTEGER, PARAMETER :: n_chars = netcdf_varname_len + t_units_len +            &
                                t_calendar_len + padding_b_len

! Derived type used to hold time information from a NetCDF file
! (read on PEO and broadcast to all PEs).
TYPE nc_fileinfo
  SEQUENCE
  INTEGER                           :: varid_t
  INTEGER                           :: dimid_t
  CHARACTER(LEN=netcdf_varname_len) :: varname_t
  CHARACTER(LEN=t_units_len)        :: cunits_t
  CHARACTER(LEN=t_calendar_len)     :: ccalendar_t
  CHARACTER(LEN=padding_b_len)      :: padding
END TYPE nc_fileinfo

! Padding is used to ensure total length is a mulitple of longest length

TYPE (nc_fileinfo) :: file_info

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='NETCDF_GET_TIME_INFO'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Obtain the communication group for this run
CALL gc_get_communicator (my_comm, icode)

! Obtain an object representing the multi-datatype structure which is
! going to be broadcast (see further details above) 
CALL setup_nml_type (no_of_types, mpl_nml_type, n_int_in = n_int,             &
                     n_chars_in = n_chars )

! File access by PE0, which reads the data needed and adds it
! to the derived type file_info
IF ( mype == 0 ) THEN
  ncdfid = fileid

  ! Get information on file 
  nstatus = nf90_inquire ( ncdfid, ndims32, nvars32, natts32, undimid32 ) 
  CALL netcdf_error (nstatus, RoutineName,'INQUIRE_FILE')

  ! Count number of variables with a time dimension
  time_counter = 0
  
  ! First loop over variables to search for any unidimensional ones.
  ! Then search for the time dimension (first by checking attribute axis).
  var_loop: DO v = 1, nvars32
     
    varid32 = v
    nstatus = nf90_inquire_variable (ncdfid,  varid32, varnam_local, itype,   &
              ndim_local, dimid)
    CALL netcdf_error ( nstatus, RoutineName, 'INQ_VARIABLE ' )
    
    IF ( .NOT. ndim_local == 1 ) CYCLE var_loop

    ! Get 'units' attribute if present
    cunits = ' '
    nstatus = nf90_get_att (ncdfid, varid32, 'units', cunits)
    IF ( nstatus /= nf90_noerr) THEN
      WRITE (umMessage,'(A)') TRIM(RoutineName) // 'File non cf compilent' // &
           ' - all variables must contain "unit" attribute'
      CALL umPrint(umMessage,src=RoutineName)
    END IF
    CALL netcdf_error (nstatus, RoutineName, 'GET_ATT variable ' //           &
                       TRIM(varnam_local) // ' in file: ' // TRIM(filename) )
    
    ! Get 'axis' attribute if present
    caxis = ' '
    nstatus = nf90_get_att (ncdfid, varid32, 'axis', caxis)
    IF ( nstatus == nf90_noerr) THEN
      
      ! Attribute exists - check for string
      IF ( .NOT. caxis == 'T' ) CYCLE var_loop
      
      IF ( .NOT. ( (cunits(1:10) == 'days since')  .OR.                       &
                   (cunits(1:11) == 'hours since') ) ) THEN
        ! The NetCDF reading routines expect time information
        ! in the form "days since" or "hours since"
        ierror   = 54902
        cmessage = 'Units of time varible ' // TRIM(varnam_local) //          &
                   ' with axis==T does not match expected form' //            &
                   ' "days since" or "hours since".'
        CALL ereport(RoutineName, ierror, cmessage)
      END IF
      
      CALL netcdf_get_time_var (filename, fileid, varid32)
      ! Found time. Force time_counter to 1
      time_counter = 1
      EXIT var_loop  ! Exit because found
      
    ELSE ! 'axis' not an attribute
      IF ( .NOT. ( (cunits(1:10) == 'days since')  .OR.                       &
                   (cunits(1:11) == 'hours since') ) ) THEN
        ! Cycle to next variable
        CYCLE var_loop
      END IF
      
      CALL netcdf_get_time_var (filename, fileid, varid32)
      time_counter = time_counter + 1
      
    END IF ! 'axis' not an attribute
  END DO var_loop
  
  ! Populate global variables holding information on time dimension.
  varid_t      = varid32
  dimid_t      = dimid(1)
  tim_name     = TRIM(varnam_local)
  ct_units     = TRIM(cunits)
  
  ! If matching variable not found - abort
  IF ( TRIM(tim_name) == 'unset' .OR. varid_t == -99 )   THEN
     nstatus = nf90_enotvar   ! Status = Variable not found
     CALL netcdf_error (nstatus,RoutineName,'Time variable missing')
  END IF

  ! If matching variable not found - abort
  IF ( time_counter > 1 )   THEN
     ierror   = 54903
     cmessage = 'There are more than one 1D variables describing the time' // &
          ' dimension. Include "axis==T" attribute in the relevant variable'
     CALL ereport(RoutineName, ierror, cmessage)
  END IF

  ! Second loop over variables to check that the time dimension is
  ! the 4th one in any 4-D variable.
  ! Note that this condition is enough because all NetCDF files should
  ! contain 4 dimensions (lon, lat, lev, time). Lev is the 3rd dim 
  ! (always present but with only one element in same cases such as
  ! for 2-D fields) and time is the 4th one.
  l_verify_time = .FALSE.
  DO v = 1, nvars32
    
    varid32 = v
    nstatus = nf90_inquire_variable (ncdfid,  varid32, varnam_local, itype,   &
                                     ndim_local, dimid)
    CALL netcdf_error ( nstatus, RoutineName, 'INQ_VARIABLE 4D' )
    
    IF ( ndim_local == 4 .AND. dimid(4) == dimid_t ) THEN
      l_verify_time = .TRUE.
      EXIT   ! Exit because found
    ELSE
      CYCLE  ! Go to next variable
    END IF
     
  END DO
  
  ! Time dimension does not seem to exist in any 4-D variable
  IF ( .NOT. l_verify_time ) THEN
    nstatus = nf90_enotvar   ! Status = Variable not found
    WRITE (umMessage,'(A,A)')                                                 &
         'Time dimension exists, but missing from data variables: ',          &
          TRIM(tim_name)
    CALL netcdf_error (nstatus, RoutineName, umMessage)
  END IF

  padding_local = 'ABCDEFGH'

  ! Copy variables into broadcast type
  file_info % varid_t     = varid_t
  file_info % dimid_t     = dimid_t
  file_info % varname_t   = tim_name
  file_info % cunits_t    = ct_units
  file_info % ccalendar_t = ct_calendar
  file_info % padding     = padding_local
END IF   ! PE0 work

! Sync PEs and broadcast time dimen-related values
CALL gc_gsync(nproc,icode)
icode = 0

! Data structure broadcast to all processors
CALL mpl_bcast (file_info, one, mpl_nml_type, zero, my_comm, icode)

IF ( icode /= 0 ) THEN
  ecode = ncdfid
  CALL ereport (RoutineName, ecode, 'Broadcast Error')
END IF

! Copy data broadcast from PE0 to PE local variables.
IF ( mype /= 0 ) THEN
  varid_t       = file_info % varid_t
  dimid_t       = file_info % dimid_t
  tim_name      = file_info % varname_t
  ct_units      = file_info % cunits_t
  ct_calendar   = file_info % ccalendar_t
  padding_local = file_info % padding
END IF

! Reset the multi-datatype structure allocated for broadcasting
CALL mpl_type_free (mpl_nml_type, icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE netcdf_get_time_info

! --------------------------------------------------------------------------
! Description:
!   Works out which variable contains the time dimension
!
! --------------------------------------------------------------------------

SUBROUTINE netcdf_get_time_var (filename, fileid, varid32)

IMPLICIT NONE

! Subroutine arguments
CHARACTER (LEN=*),        INTENT(IN)    :: filename
INTEGER,                  INTENT(IN)    :: fileid
INTEGER (KIND=integer32), INTENT(INOUT) :: varid32

! Local variables
INTEGER (KIND=integer32)         :: nstatus      ! netCDF return code
CHARACTER (LEN=t_calendar_len)   :: ccalendar    ! For 'calendar' attribute

INTEGER (KIND=jpim), PARAMETER   :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER   :: zhook_out = 1
REAL (KIND=jprb)                 :: zhook_handle

CHARACTER (LEN=*), PARAMETER     :: RoutineName='NETCDF_GET_TIME_VAR'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF ( mype == 0 ) THEN
  ncdfid = fileid
  ! Now check 'calendar' attribute if present
  ccalendar = ' '
  nstatus = nf90_get_att (ncdfid, varid32, 'calendar', ccalendar)
  
  ! Return if variable does not contain calendar attribute
  IF (nstatus /= nf90_noerr) THEN
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,          &
                            zhook_handle)
    RETURN
  END IF
  
  IF ( (                          ccalendar == 'none'       )   .OR.          &
       (          lcal360  .AND. (ccalendar == '360_day')   )   .OR.          &
       (   (.NOT. lcal360) .AND. (ccalendar == 'standard')  )   .OR.          &
       (   (.NOT. lcal360) .AND. (ccalendar == 'gregorian') ) ) THEN
    ! ct_calendar gets saved
    ct_calendar  = TRIM(ccalendar)                  
    WRITE (umMessage,'(A,A)') 'Calendar = ',TRIM(ccalendar)
    CALL umPrint(umMessage,src=RoutineName)
  ELSE IF ( (     lcal360  .AND. (ccalendar == 'standard')  )   .OR.          &
            (     lcal360  .AND. (ccalendar == 'gregorian') ) ) THEN
    ecode = fileid
    CALL ereport (RoutineName, ecode,                                        &
         'Expecting 360 calendar.  Instead file: ' // TRIM(filename) //       &
         ' contains calendar format: ' // TRIM(ccalendar) )
  ELSE IF ((.NOT. lcal360) .AND. (ccalendar == '360_day')     ) THEN
    ecode =fileid
    CALL ereport (RoutineName, ecode,                                        &
         'Not expecting 360 calendar.  Instead file: ' // TRIM(filename) //   &
         ' contains calendar format: ' // TRIM(ccalendar) )
  ELSE IF ( (ccalendar == 'proleptic_gregorian')              .OR.            &
            (ccalendar == 'julian')                           .OR.            &
            (ccalendar == 'noleap')                           .OR.            &
            (ccalendar == '365_day')                          .OR.            &
            (ccalendar == 'all_leap')                         .OR.            &
            (ccalendar == '366_day') )                        THEN
    ecode = fileid
    CALL ereport (RoutineName, ecode,                                        &
         'Calendar ' // TRIM(ccalendar) // ' in file : ' // TRIM(filename) // &
         ' not currently supported. Please use "360_day" or "standard" ' //   &
         ' or "gregorian".' )
   ELSE
     ecode = fileid
     CALL ereport (RoutineName, ecode, 'Calendar ' // TRIM(ccalendar) //     &
          ' in file : ' // TRIM(filename) // ' non cf compilent. ' //         &
          'Please use "360_day" or "standard" or "gregorian" calendar.')
   END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN

END SUBROUTINE netcdf_get_time_var

! --------------------------------------------------------------------------
! Description:
!   Check that dimensions of a variable in a NetCDF file match global sizes.
!
! Method:
! 1) Call NETCDF_GET_VAR_INFO to get the dimensions and check they match.
! 2) If wrong dimensions or if there is no time record then stop with
!    call to EREPORT
! --------------------------------------------------------------------------
SUBROUTINE netcdf_var_check_dims (row_length_in, rows_in,             &
                                  levels_in, fileid, varname, dimsize)

IMPLICIT NONE

! Subroutine arguments
INTEGER,                            INTENT(IN)    :: row_length_in
INTEGER,                            INTENT(IN)    :: rows_in
INTEGER,                            INTENT(IN)    :: levels_in
INTEGER,                            INTENT(IN)    :: fileid
INTEGER,                            INTENT(OUT)   :: dimsize(max_dims_per_var)
CHARACTER (LEN=netcdf_varname_len), INTENT(INOUT) :: varname

! Local variables
INTEGER                        :: vid    ! variable id
LOGICAL                        :: lreco  ! true if contains record dim

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='NETCDF_VAR_CHECK_DIMS'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Climatology files should be defined as (x,y,lev,t)
!
! Get information on variable. One of the two INOUT arguments
! of NETCDF_GET_VAR_INFO is supplied (varname); assign an arbitrary
! negative value to the other one since it is not known yet (vid).
vid = -1
CALL netcdf_get_var_info (fileid, vid, varname, vdimsize=dimsize, lrec=lreco)

! Check that dimensions (x, y, lev) match those of the model
IF ( dimsize(1) /= row_length_in .OR.                                         &
     dimsize(2) /= rows_in .OR.                                               &
    (dimsize(3) /= levels_in .AND. dimsize(3) /= 1) ) THEN
  
  WRITE (umMessage,'(A,A)') 'ERROR: Variable dimensions do not match for ',   &
                             TRIM(varname)
  CALL umPrint(umMessage,src=RoutineName)
  
  WRITE (umMessage,'(A,3(1x,I5))') 'Expected (multi-layer):',                 &
                                    row_length_in,                            &
                                    rows_in,                                  &
                                    levels_in
  CALL umPrint(umMessage,src=RoutineName)
  
  WRITE (umMessage,'(A,3(1x,I5))') 'Expected (single layer):',                &
                                    row_length_in,                            &
                                    rows_in,                                  &
                                    1
  CALL umPrint(umMessage,src=RoutineName)
  
  WRITE (umMessage,'(A,3(1x,I5))') 'Found dimsize(1:3) : ', dimsize(1:3)
  CALL umPrint(umMessage,src=RoutineName)
  
  ecode = 54908
  CALL ereport (RoutineName, ecode,                                           &
       'Dimension mismatch for: ' // TRIM(varname))
END IF

! Check that there is time record
IF (.NOT. lreco ) THEN
  WRITE (umMessage,'(A,A)') 'No record dimension found for ', TRIM(varname)
  CALL umPrint(umMessage,src=RoutineName)
  ecode = fileid
  CALL ereport (RoutineName, ecode,                                           &
               'Dimension mismatch for: ' // TRIM(varname))
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE netcdf_var_check_dims

! --------------------------------------------------------------------------
! Description:
!   Get time information (in days and seconds) for specified record.
!
! Method:
!   The name and units attribute of the time variable have been read
!   previously during an NETCDF_FOPEN call. They are stored in the
!   global variables 'tim_name' and 'ct_units', respectively. 
!   This routine uses that information to:
! 1) Check that attribute units of time variable (i.e. 'ct_units') is
!    'days' or 'hours'.
! 2) Get time origin, given by 'ct_units',
!    (expected to be something like "days since yyyy-01-01 00:00:00"  or
!                                  "hours since yyyy-01-01 00:00:00")
!    and decode it. Once decoded pass it to TIME2SEC to get that starting
!    time in days and seconds: dayfilest, secfilest.
! 3) Call NETCDF_GET_DATA (passing 'tim_name' as argument) to get the values
!    of time in hours (ftime). Decode ftime to get that time in 
!    days and seconds: tday, tsec.
! 4) Update (tday, tsec) by adding the starting times (dayfilest, secfilest)
!    obtained in (2).
! --------------------------------------------------------------------------
SUBROUTINE netcdf_get_time_rec ( fileid, irec, iupdn_type,                    &
                                 l_first, l_new_year, tday, tsec, l_notime)

IMPLICIT NONE

! Subroutine arguments

INTEGER, INTENT(IN)  :: fileid     ! File id
INTEGER, INTENT(IN)  :: irec       ! Which record to check
INTEGER, INTENT(IN)  :: iupdn_type ! update_type: iupdn_serial = 1
                                   !              iupdn_cyclic = 2

LOGICAL, INTENT(IN)  :: l_first    ! First time
LOGICAL, INTENT(IN)  :: l_new_year ! True if we need to get time
                                   ! for a new year

INTEGER, INTENT(OUT) :: tday       ! time (in days) found for irec
INTEGER, INTENT(OUT) :: tsec       ! time (in secs) found for irec

LOGICAL, INTENT(OUT) :: l_notime   ! True if no time record in file

! Local variables

REAL             :: ftime(1)       ! Fractional time read from the file
                                   ! Decl as (1) to match nc_get() args

REAL             :: ftfrac         ! used to separate days + sec in ftime

! Dummy characters while decoding the time attribute "units"
CHARACTER(LEN=1) :: cd1            ! for '-', ' ' or ':'
CHARACTER(LEN=4) :: cdumm4         ! for "days"
CHARACTER(LEN=5) :: cdumm5         ! for "hours"
INTEGER          :: file_origt(6)  ! File origin time, also read from units
                                   ! 6 posit to hold: yyyy, mm, dd, hh, mm, ss
INTEGER          :: rec_time(7)    ! Record validity time
                                   ! 7 posit to hold: y,m,d,h,m,s,daynum

INTEGER          :: dayfilest, secfilest ! file orig time in days/sec
INTEGER          :: rec_yr         ! year of file record (imposed if cyclic)


LOGICAL          :: l_unit_dy      ! true if time in days, otherwise in hours

REAL,  PARAMETER :: hour_to_day = 1.0/24.0  ! conv factor, hours to days

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0  ! for DrHook
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER    :: RoutineName='NETCDF_GET_TIME_REC'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Decode what is held in units for time variable, which should be either 
!    "days since yyyy-01-01 00:00:00"
!   "hours since yyyy-01-01 00:00:00"
! Only units allowed are 'days' or 'hours'.
IF ( ct_units(1:4) /= 'days' .AND. ct_units(1:4) /= 'hour' ) THEN
     ecode = fileid
     CALL ereport(RoutineName,ecode,                                          &
                 'Time units should be in days or hours only')
END IF

l_unit_dy = ( ct_units(1:4) == 'days' )    ! true if units are days

! Continue decoding the units attribute. Now get file_orig (1:6),
! which holds (year, month, day, hour, min, sec)
IF ( l_unit_dy ) THEN     ! assume start with the 4-letter word 'days'
  READ (ct_units,                                                             &
      '(A4,1x,A5,1x,I4,A1,I2,A1,I2,1x,I2,A1,I2,A1,I2)')                       &
       cdumm4, cdumm5,                                                        & 
       file_origt(1), cd1, file_origt(2), cd1, file_origt(3),                 &
       file_origt(4), cd1, file_origt(5), cd1, file_origt(6)

ELSE                       ! assume start with the 5-letter word 'hours'
  READ (ct_units,                                                             &
     '(A5,1x,A5,1x,I4,A1,I2,A1,I2,1x,I2,A1,I2,A1,I2)')                        &
      cdumm5, cdumm5,                                                         &
      file_origt(1), cd1, file_origt(2), cd1, file_origt(3),                  &
      file_origt(4), cd1, file_origt(5), cd1, file_origt(6)
END IF

IF ((file_origt(1) < 1) .AND. (iupdn_type == iupdn_cyclic) ) THEN
  file_origt(1) = 1
  IF (l_first .AND. irec == 1) THEN
    IF (PrintStatus >= PrStatus_Diag) THEN
      umMessage = "File basis year for periodic file increased to 0001 " //   &
                  "because time2sec doesn't handle smaller years."
      CALL umPrint(umMessage,src=RoutineName)
    END IF
  END IF
END IF

! Get file data start time in days/secs, e.g. from 0000-01-01.
! Pass file_orig (1:6) to time2sec and get the starting time as
! days and secs since basis time: dayfilest + secfilest
! DEPENDS ON: time2sec
CALL time2sec (file_origt(1), file_origt(2), file_origt(3),                   &
               file_origt(4), file_origt(5), file_origt(6),                   &
               0, 0, dayfilest, secfilest, lcal360 )

! Print out time origin in NetCDF file if first time and extra diags.
IF (l_first .AND. irec == 1) THEN
  IF (PrintStatus >= PrStatus_Diag) THEN
    WRITE (umMessage,'(A,I4,5(1x,I2),2(1x,I8))') TRIM(RoutineName) //         &
           ' : File starts ', file_origt(1:6), dayfilest, secfilest
    CALL umPrint(umMessage,src=RoutineName)
  END IF
END IF

! ----------------------------------------------------------------------------
! Get record for requested time (i.e. record irec from variable
! with name tim_name) and store it in ftime
CALL netcdf_get_data (fileid, irec, tim_name, ftime, l_notime)

IF ( L_notime ) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
  RETURN   ! no time record found in file
END IF

! If time is in hours then convert to days
IF (.NOT. l_unit_dy) ftime(:) = ftime(:) * hour_to_day

! Now calculate time in full days + seconds: tday + tsec.
! For that check if ftime value is full days or decimal.
! This is done in a crude way, which can lead to very
! low decimal diffs at REAL8. Seconds may be rounded off,
! but this is not expected to have any impact.
ftfrac = REAL  ( INT(ftime(1)) )                    ! Fractional day
tsec   = INT   ( (ftime(1)-ftfrac) * isec_per_day ) ! Fractional part
                                                    ! of ftime in secs
tday   = INT   (ftime(1))                           ! Full days

! Finally, get the final time in days and seconds by adding the starting
! days and seconds that we read from the attribute 'units'.
tday = tday + dayfilest
tsec = tsec + secfilest

! If file is periodic then ignore year. Note that year is set to 
! cyclic_year_netcdf in 
! order to be consistent with glomap_clim_get_netcdffile_rec subroutine.
! To change the year, it converts day/sec to a date, adjusts the year,
! and converts back to day/sec.
IF (iupdn_type == iupdn_cyclic) THEN
  rec_yr = cyclic_year_netcdf

  ! Add one year if indicated by l_new_year
  IF (l_new_year) THEN
    rec_yr = rec_yr + 1
  END IF

  ! DEPENDS ON: sec2time
  CALL sec2time(tday, tsec, 0, 0,                                             &
                rec_time(1), rec_time(2), rec_time(3),                        &
                rec_time(4), rec_time(5), rec_time(6),                        &
                rec_time(7), lcal360)

  rec_time(1) = rec_yr

  ! DEPENDS ON: time2sec
  CALL time2sec(rec_time(1), rec_time(2), rec_time(3),                        &
                rec_time(4), rec_time(5), rec_time(6),                        &
                0, 0, tday, tsec, lcal360 )


END IF

! tday and tsec can now be passed to other climatology routines, where
! they will be compared with the current model time in order to
! perform interpolations of the climatology fields.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE netcdf_get_time_rec

!-------------------------------------------------------------------------
! Description:
!   Read data values of one-dimensional real variable in NetCDF file.
!   Part of the interface block NETCDF_GET_DATA used for overloading; this
!   subroutine is called when the variable is of type real 1-dim.
!
! Method:
!   Get variable ID (fileid), number of dimensions in variable (ndim), 
!   and vector of dimension sizes (dim_size) by calling NETCDF_GET_VAR_INFO.
!   Then get domain/indices to extract and call NF90_GET_VAR to get the 
!   variable values: rvalue. Values are read by PE0 and broadcast
!   to all PEs.
!   If the netCDF variable has a time dimension, then only the data for
!   time index = trec are read. In this case the netCDF variable is 1D and
!   a single value is read, or is 2D and a 1D slice is read. It is assumed
!   that time is the last dimension for this variable.
!-------------------------------------------------------------------------
SUBROUTINE netcdf_get_data_real_1D (fileid, trec, varname, rvalue, l_norec)

IMPLICIT NONE

! Subroutine arguments
INTEGER,           INTENT (IN)    :: fileid     ! file id
INTEGER,           INTENT (IN)    :: trec       ! time record wanted
CHARACTER (LEN=*), INTENT (IN)    :: varname    ! variable name
REAL,              INTENT (INOUT) :: rvalue (:) ! variable values to read
LOGICAL,           INTENT (OUT)   :: l_norec    ! T if rec not in file, EOF

! Local variables


INTEGER (KIND=integer32) :: nstatus ! netCDF return code
!  Vector specifying the index in the variable from which
!  the first data value will be read along each dimension
INTEGER (KIND=integer32) :: start   (max_dims_per_var)
!  Vector specifying the number of counters to read along each dimension
INTEGER (KIND=integer32) :: counter (max_dims_per_var)
INTEGER (KIND=integer32) :: varid32
INTEGER                  :: varid, ndim, dim_size(max_dims_per_var), n
INTEGER                  :: i_len, i_get_dims
LOGICAL                  :: lreco    ! true if contains time record dim

CHARACTER (LEN=netcdf_varname_len) :: varname_local

! For Parallel I/O
INTEGER :: my_comm
INTEGER :: icode

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER    :: RoutineName='NETCDF_GET_DATA_REAL_1D'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Obtain the communication group for this run
CALL gc_get_communicator (my_comm, icode)

l_norec = .FALSE. ! Initialised to false

! Get information on variable. One of the two INOUT arguments of
! NETCDF_GET_VAR_INFO is supplied (varname_local); assign an arbitrary
! negative value to the other one since it is not known yet (varid).
varname_local = varname    ! need to pass INOUT argument
varid         = -1
CALL netcdf_get_var_info (fileid, varid, varname_local,                       &
                          ndims=ndim, vdimsize=dim_size, lrec=lreco)

! Compute i_len: total size needed for upper bound of assumed-size
! array used in nf90_get_var. The total array size is the product
! of all dimension sizes (excluding last dimension if lreco=true).
IF (lreco) THEN
  i_get_dims = ndim - 1
ELSE
  i_get_dims = ndim
END IF
i_len = 1
DO n =1, i_get_dims
  i_len = i_len * dim_size(n)
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
    l_norec = .TRUE.
    
    ! Record out of bounds and return that no record found
    IF ( trec < 1 .OR. trec > dim_size(ndim) ) THEN
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,        &
                              zhook_handle)
      RETURN
    END IF
    
    l_norec        = .FALSE.  ! Indicate to main routine that time rec found
    
    ! Read a single time from the file
    start   (ndim) = trec     ! last index corresponds to record number
    counter (ndim) = 1        ! read only one record
  END IF
  
  ! Check that array with data values to read (rvalue) has the right dimensions
  IF (SIZE (rvalue,1) /= counter (1)) THEN
    ecode = 54905
    CALL ereport (RoutineName, ecode,                                         &
                 'Dimension mismatch for: ' // TRIM(varname))
  END IF
  
  !  Get values of the variable as rvalue, specifying what to read
  !  with the optional start and counter arguments
  ncdfid  = fileid
  varid32 = varid    ! 32-bit integer needed
  nstatus = nf90_get_var (ncdfid, varid32, rvalue, start (1:ndim),            &
                          counter (1:ndim))
  CALL netcdf_error (nstatus, RoutineName, varname)
END IF   ! I/O PE

! Check that the array to be broadcast has the same size as the buffer length
! given to mpl_bcast
IF (SIZE(rvalue,1) /= i_len) THEN
  ecode = 54000 + i_len
  CALL ereport (RoutineName, ecode,                                           &
                'Dimension mismatch for: ' // TRIM(varname))
END IF

! Distribute data to all other PEs
icode = 0
CALL mpl_bcast (rvalue,  i_len, mpl_real,    zero, my_comm, icode)
IF ( icode /= 0 ) THEN
  ecode = fileid
  CALL ereport (RoutineName, ecode, 'Broadcast Error')
END IF

CALL mpl_bcast (l_norec, one,   mpl_logical, zero, my_comm, icode)
IF ( icode /= 0 ) THEN
  ecode = fileid 
  CALL ereport (RoutineName, ecode, 'Broadcast Error')
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE netcdf_get_data_real_1D

!-------------------------------------------------------------------------
! Description:
!   Read data values of three-dimensional real variable in NetCDF file.
!   Part of the interface block NETCDF_GET_DATA used for overloading; this
!   subroutine is called when the variable is of type real 3-dim.
!
! Method:
!   Get variable ID (varid), number of dimensions in variable (ndim), 
!   and vector of dimension sizes (dim_size) by calling NETCDF_GET_VAR_INFO.
!   Then get domain/indices to extract and call NF90_GET_VAR to get the 
!   variable values, which are read by PEO as global_var. The values
!   are broadcast to all PEs, looping over levels, and returned as rvalue.
!   Assumes (lon, lat, lev).
!-------------------------------------------------------------------------
SUBROUTINE netcdf_get_data_real_field ( glob_row_length_in, glob_rows_in,     &
                                        loc_row_length_in, loc_rows_in,       &
                                        levels_in, fileid,                    &
                                        output_decomposition,                 &
                                        trec, varname, rvalue, l_norec )

USE field_types,  ONLY: &
    fld_type_p

USE UM_ParParams, ONLY: &
    halo_type_no_halo

IMPLICIT NONE

! Subroutine arguments
INTEGER,           INTENT (IN)    :: glob_row_length_in ! global horizontal dim
INTEGER,           INTENT (IN)    :: glob_rows_in       ! global horizontal dim
INTEGER,           INTENT (IN)    :: loc_row_length_in  ! local horizontal dim
INTEGER,           INTENT (IN)    :: loc_rows_in        ! local horizontal dim
INTEGER,           INTENT (IN)    :: levels_in          ! model_levels or 1
INTEGER,           INTENT (IN)    :: fileid             ! file id
INTEGER,           INTENT (IN)    :: output_decomposition ! Setup output
                                                          ! decomposition for
                                                          ! writing
INTEGER,           INTENT (IN)    :: trec               ! time record wanted
LOGICAL,           INTENT (OUT)   :: l_norec            ! T if rec not in file
                                                        !  or reached EOF
CHARACTER (LEN=*), INTENT (IN)    :: varname            ! variable name
REAL,              INTENT (INOUT) :: rvalue (:,:)

! Local variables

INTEGER (KIND=integer32)           :: nstatus ! netCDF return code
! Vector specifying the index in the variable from which
! the first data value will be read along each dimension.
INTEGER (KIND=integer32)           :: start   (max_dims_per_var)
! Vector specifying the number of counters to read along each dimension
INTEGER (KIND=integer32)           :: counter (max_dims_per_var)
INTEGER (KIND=integer32)           :: varid32
INTEGER                            :: varid, ndim, dim_size(max_dims_per_var)
INTEGER                            :: orig_decomp
LOGICAL                            :: lreco ! true if contains time record dim
CHARACTER (LEN=netcdf_varname_len) :: varname_local

REAL, ALLOCATABLE           :: global_var(:,:,:) ! Variable to hold file data

INTEGER :: k

! For Parallel I/O
INTEGER :: my_comm, icode

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='NETCDF_GET_DATA_REAL_FIELD'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Obtain the communication group for this run
CALL gc_get_communicator (my_comm, icode)

l_norec = .FALSE. ! Initialised to false

! Get information on variable (need to do this from all PEs so that
! global_var can be allocated with the right 3rd dimension)
!
! One of the two INOUT arguments of NETCDF_GET_VAR_INFO is supplied
! (varname_local); assign an arbitrary negative value to the other 
! one since it is not known yet (varid).
varname_local = varname    ! need to pass INOUT argument
varid         = -1
CALL netcdf_get_var_info ( fileid, varid, varname_local,                      &
                           ndims=ndim, vdimsize=dim_size, lrec=lreco )

! Check that data in NetCDF file matches Global Dimensions and 
! output variable rvalue has correct space for levels_in
CALL netcdf_var_check_dims ( glob_row_length_in, glob_rows_in,                &
                             levels_in, fileid, varname_local, dim_size )

! Allocate the variable that will hold the data.
! It will be filled by PE0 and broadcast to other PEs.
ALLOCATE ( global_var (dim_size(1), dim_size(2), dim_size(3)) )

! Get values of variable from PE0 and store in global_var
IF ( mype == 0 ) THEN
  ! Default indices for the first data value and number of counters
  ! along two dimensions: Global lon and lat. Note that counter(3) can
  ! be different from the number of model vertical levels.
  start   (1:3) = 1
  counter (1)   = dim_size (1)
  counter (2)   = dim_size (2)
  counter (3)   = dim_size (3)
  
  ! If var contains time dimension then select record to read
  IF (lreco) THEN
    l_norec = .TRUE.
    
    ! Record out of bounds and return that no record found
    IF ( trec > dim_size(ndim) ) THEN
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,        &
                              zhook_handle)
      RETURN
    END IF
    
    l_norec        = .FALSE.  ! Indicate to main routine that time rec found
    
    ! Read a single time from the file
    start   (ndim) = trec     ! last index corresponds to record number
    counter (ndim) = 1        ! read only one record
  END IF
  
  ! Get values of the variable as global_var, specifying what to read
  ! with the optional start and counter arguments
  ncdfid  = fileid
  varid32 = varid    ! 32-bit integer needed
  nstatus = nf90_get_var ( ncdfid, varid32, global_var, start(1:ndim),        &
                           counter(1:ndim) )
  CALL netcdf_error (nstatus, RoutineName, varname)
END IF    ! I/O PE

! Distribute data to all other PEs
icode = 0
CALL mpl_bcast (l_norec, one,   mpl_logical, zero, my_comm, icode)
IF ( icode /= 0 ) THEN
  ecode = fileid  
  CALL ereport (RoutineName, ecode, 'Broadcast Error')
END IF

CALL gc_gsync(nproc,icode)   ! Synchronise PEs before scatter

! ===============================================
! Setup output decomposition for writing.
orig_decomp = current_decomp_type
IF (orig_decomp /= output_decomposition) THEN
  CALL Change_Decomposition( output_decomposition )
END IF

! Scatter data, looping over levels
DO k = 1, dim_size(3)
  !   Scatter_Field_Real expects in this order:
  !   returned_local_field(:), global_field(:), local_row_len, local_rows,    &
  !   global_row_len, global_rows, fld_type, halo_type, scatter_pe, proc_group
  
  !DEPENDS ON: scatter_field
  CALL scatter_field ( rvalue(:,k),                                           &
          RESHAPE(global_var(:,:,k),(/ glob_row_length_in * glob_rows_in /) ),&
                       loc_row_length_in,                                     &
                       loc_rows_in,                                           &
                       glob_row_length_in,                                    &
                       glob_rows_in,                                          &
                       fld_type_p,                                            &
                       halo_type_no_halo,                                     &
                       zero,                                                  &
                       gc_all_proc_group )
END DO

DEALLOCATE(global_var)

! ===============================================
! Change back to original if different.
IF (current_decomp_type /= orig_decomp) THEN
  CALL Change_Decomposition( orig_decomp )
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE netcdf_get_data_real_field

! --------------------------------------------------------------------------
! Description:
!   Read attributes of type character from a NetCDF file, given the
!   variable and attribute names.
!   Part of the interface block NETCDF_GET_VAR_ATT used for overloading; this
!   subroutine is called when the required attribute is of type character.
!
! Method:
!   If a global attribute is required, variable name = 'global'.
!   If variable name is passed call NF90_INQ_VARID to get variable-id.
!   Call NF90_GET_ATT to get attribute value: cvalue.
!   Attribute values are read on PE0 and broadcast to all PEs.
! --------------------------------------------------------------------------
SUBROUTINE netcdf_get_var_att_char (fileid, i_len, varname, att_name, cvalue, &
                                    att_exists)

IMPLICIT NONE

!  Subroutine arguments
INTEGER,          INTENT(IN)     :: fileid    ! file id
INTEGER,          INTENT(IN)     :: i_len     ! num items in value
                                              ! (not used for chars)
CHARACTER(LEN=*), INTENT(IN)     :: varname   ! variable name
CHARACTER(LEN=*), INTENT(IN)     :: att_name  ! attribute name
CHARACTER(LEN=*), INTENT(OUT)    :: cvalue    ! value out

! Optional argument used to check if an attribute is present. If
! attribute found then l_exist will be set to True, otherwise
! l_exist will be False but the model will not stop with error.
LOGICAL, OPTIONAL, INTENT(INOUT) :: att_exists

!  Local variables
INTEGER   (KIND=integer32) :: varid32
INTEGER   (KIND=integer32) :: nstatus          ! netCDF return code
INTEGER,  PARAMETER        :: charlength = 256
INTEGER                    :: idummy1          ! only used to compile code
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

CHARACTER(LEN=*), PARAMETER :: RoutineName='NETCDF_GET_VAR_ATT_CHAR'

!  End of header
   
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! This dummy line will get rid of compilation error
idummy1 = i_len

! Obtain the communication group for this run
CALL gc_get_communicator(my_comm, icode)

! Obtain an object representing the multi-datatype structure which is
! going to be broadcast (see further details above)
CALL setup_nml_type(no_of_types, mpl_nml_type, n_log_in = n_log,              &
                    n_chars_in = n_chars)

! PE0 reads attribute values and adds them to the derived type file_info
IF ( mype == 0 ) THEN
  ncdfid = fileid
  ! Get variable id
  IF ( TRIM(varname) == 'GLOBAL' .OR.                                         &
       TRIM(varname) == 'global' .OR.                                         &
       TRIM(varname) == 'Global' ) THEN
    varid32 = NF90_GLOBAL
  ELSE
    nstatus = nf90_inq_varid (ncdfid, varname, varid32)
    CALL netcdf_error (nstatus, TRIM(RoutineName) // ' : INQ Varid', varname)
  END IF

  ! Get value of attribute
  cval    = ' '
  nstatus = nf90_get_att (ncdfid, varid32, att_name, cval)

  ! If argument att_exists is present that means that the attribute is not
  ! compulsory. Therefore, if attribute found in the previous call to 
  ! nf90_get_att then set the local variable l_exist to True; otherwise 
  ! leave it as False but there is no need to flag error.
  l_exist = .FALSE.
  IF ( PRESENT(att_exists) ) THEN
    IF (nstatus == nf90_noerr) l_exist = .TRUE. ! attrib exists

  ! If argument att_exists not present then the attribute is needed. Therefore
  ! call netcdf_error, which will stop the model only if the previous call to
  ! nf90_get_att reported error (as given by nstatus).
  ELSE
    CALL netcdf_error (nstatus, TRIM(RoutineName) // ' : Get Attribute for '  &
                             // TRIM(varname), att_name) ! Check for error
  END IF

  ! Fill in output argument.
  cvalue = TRIM(cval)

  ! Copy into broadcast type
  file_info % att_value = cvalue
  file_info % l_exist   = l_exist
END IF

! Distribute data to all other PEs
icode = 0
CALL mpl_bcast (file_info, one, mpl_nml_type, zero, my_comm, icode)
IF ( icode /= 0 )  THEN
  ecode = fileid
  CALL ereport (RoutineName, ecode, 'NML Broadcast Error')
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
END SUBROUTINE netcdf_get_var_att_char

! --------------------------------------------------------------------------
! Description:
!   Read attributes of type integer from a NetCDF file, given the
!   variable and attribute names.
!   Part of the interface block NETCDF_GET_VAR_ATT used for overloading; this
!   subroutine is called when the required attribute is of type integer.
!
! Method:
!   If a global attribute is required, variable name = 'global'.
!   If variable name is passed call NF90_INQ_VARID to get variable-id.
!   Call NF90_GET_ATT to get attribute value: ivalue.
!   Attribute values are read on PE0 and broadcast to all PEs.
! --------------------------------------------------------------------------
SUBROUTINE netcdf_get_var_att_int (fileid, i_len, varname, att_name, ivalue,  &
                                   att_exists )

IMPLICIT NONE

!  Subroutine arguments
INTEGER,          INTENT(IN)     :: fileid        ! file id
INTEGER,          INTENT(IN)     :: i_len         ! num items in value (just
                                                  ! in case it is an array)
CHARACTER(LEN=*), INTENT(IN)     :: varname       ! variable name
CHARACTER(LEN=*), INTENT(IN)     :: att_name      ! attribute name
INTEGER,          INTENT(OUT)    :: ivalue(i_len) ! value out

! Optional argument used to check if an attribute is present. If
! attribute found then l_exist will be set to True, otherwise
! l_exist will be False but the model will not stop with error.
LOGICAL, OPTIONAL, INTENT(INOUT) :: att_exists

!  Local variables
LOGICAL :: l_exist

! 32-bit integers to actually pass to NC routines
INTEGER(KIND=integer32) :: varid32, ival32(i_len)
INTEGER(KIND=integer32) :: nstatus               ! netCDF return code

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

CHARACTER(LEN=*), PARAMETER :: RoutineName='NETCDF_GET_VAR_ATT_INT'

!  End of header
   
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Obtain the communication group for this run
CALL gc_get_communicator(my_comm, icode)

! Obtain an object representing the multi-datatype structure which is
! going to be broadcast (see further details above)
CALL setup_nml_type (no_of_types, mpl_nml_type, n_int_in = n_int,             &
                     n_log_in = n_log)

! PE0 reads attribute values and adds them to the derived type file_info
IF ( mype == 0 ) THEN
  ! Get variable id
  ncdfid = fileid
  IF ( TRIM(varname) == 'GLOBAL' .OR.                                         &
       TRIM(varname) == 'global' .OR.                                         &
       TRIM(varname) == 'Global' ) THEN
    varid32 = NF90_GLOBAL
  ELSE
    nstatus = nf90_inq_varid (ncdfid, varname, varid32)
    CALL netcdf_error (nstatus, TRIM(RoutineName) // ' : INQ Varid', varname)
  END IF

  ! Get values of attribute
  ival32(:) = 0
  nstatus   = nf90_get_att (ncdfid, varid32, att_name, ival32)

  ! If argument att_exists is present that means that the attribute is not
  ! compulsory. Therefore, if attribute found in the previous call to 
  ! nf90_get_att then set the local variable l_exist to True; otherwise 
  ! leave it as False but there is no need to flag error.
  l_exist = .FALSE.
  IF ( PRESENT(att_exists) ) THEN
    IF (nstatus == nf90_noerr) l_exist = .TRUE.  ! attrib exists

  ! If argument att_exists not present then the attribute is needed. Therefore
  ! call netcdf_error, which will stop the model only if the previous call to
  ! nf90_get_att reported error (as given by nstatus).
  ELSE
    CALL netcdf_error (nstatus,TRIM(RoutineName) // ' : Get Attribute for '   &
                            // TRIM(varname), att_name)
  END IF

  ! Store local variables on broadcast type
  file_info % l_exist           = l_exist
  file_info % att_value(1:i_len) = ival32(1:i_len)
END IF

! Distribute data to all other PEs
icode = 0
CALL mpl_bcast (file_info, one, mpl_nml_type, zero, my_comm, icode)
IF ( icode /= 0 ) THEN
  ecode = fileid
  CALL ereport (RoutineName, ecode, 'NML Broadcast Error')
END IF

! Operation on all PEs different from PE0: Transfer data from
! broadcast type into local variables.
IF ( mype /= 0 ) THEN
  ival32(1:i_len) = file_info % att_value(1:i_len)
  l_exist         = file_info % l_exist
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
END SUBROUTINE netcdf_get_var_att_int

! ------------------------------------------------------------------------
! Description:
!   Close NetCDF climatology file.
!
! Method:
!   Call NF90_CLOSE to close a NetCDF file and return nstatus. Also call
!   NETCDF_ERROR, which will stop the model if nstatus reports errors. Note
!   that this is done only by PE0.
!   If needed reset time dimension related values since they are specific
!   to each file.
! ------------------------------------------------------------------------
SUBROUTINE netcdf_fclose (fileid)

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN)            :: fileid

! Local variables
INTEGER (KIND=integer32)       :: nstatus        ! netCDF return code

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='NETCDF_FCLOSE'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Close file only by PE0
IF ( mype == 0 ) THEN
  ncdfid    = fileid
  nstatus = nf90_close(ncdfid)
  CALL netcdf_error (nstatus, RoutineName, ' ')
END IF

! Reset time dimension related values (specific to each file)
! if they are populated
IF ( varid_t > 0 .OR. dimid_t > 0 ) THEN   ! Default -99
  varid_t  = -99
  dimid_t  = -99
  tim_name = 'unset'
  ct_units = 'unset'
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE netcdf_fclose

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
SUBROUTINE netcdf_error (nstatus, routine_name, cmessage)

IMPLICIT NONE

! Subroutine arguments
INTEGER (KIND=integer32), INTENT(IN) :: nstatus      ! netCDF return code
CHARACTER (LEN=*), INTENT(IN)        :: routine_name ! routine passing error
CHARACTER (LEN=*), INTENT(IN)        :: cmessage     ! error message

! Local variables
INTEGER (KIND=integer32)             :: ncode
INTEGER (KIND=jpim), PARAMETER       :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER       :: zhook_out = 1
REAL    (KIND=jprb)                  :: zhook_handle
CHARACTER(LEN=*), PARAMETER          :: RoutineName='NETCDF_ERROR'
CHARACTER(LEN=errormessagelength)    :: cmessageX

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF ( nstatus /= nf90_noerr ) THEN
  ! Close file only by PE0
  IF ( mype == 0 ) ncode = nf90_close(ncdfid)
  ecode = 54909
  
  ! Printing error messages first to avoid passing
  ! a very long string argument to ereport
  WRITE (umMessage,'(a)') 'Error in ' // routine_name //  ": "  //            &
                           TRIM (nf90_strerror(nstatus)) //  " - " //         &
                           TRIM (cmessage)
  CALL umPrint(umMessage,src=RoutineName)
  cmessageX = 'Error with NetCDF io'
  CALL ereport (routine_name, ecode, cmessageX)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN

END SUBROUTINE netcdf_error

END MODULE glomap_clim_netcdf_io_mod
