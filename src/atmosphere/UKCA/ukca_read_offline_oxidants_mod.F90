! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Module to contain code to enable the offline oxidant files in
!  netCDF format to be read in and updated as required.
!
! Method:
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds, University of Oxford, and
!  The Met Office. See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.
!
! Contained subroutines in this module:
!   ukca_read_offline_oxidants_ctl     (Public)
!   ukca_offline_oxidants_diags        (Public)
!   ukca_offline_oxidants_init         (Internal)
!   ukca_offline_struct_init           (Internal)
!   ukca_offline_oxidants_update       (Internal)
!   ukca_fill_offline_oxidants         (Internal)
!
! --------------------------------------------------------------------------

MODULE ukca_read_offline_oxidants_mod

USE ukca_option_mod,  ONLY: ukca_offline_dir, ukca_offline_files,              &
  l_ukca_offline, max_offline_files
USE missing_data_mod, ONLY: imdi
USE nlstcall_mod,     ONLY: lcal360, ancil_reftime
USE ukca_emiss_mod,   ONLY: ukca_em_struct
USE parkind1,         ONLY: jpim, jprb      ! DrHook
USE yomhook,          ONLY: lhook, dr_hook  ! DrHook
USE errormessagelength_mod, ONLY: errormessagelength

USE submodel_mod, ONLY: atmos_im
USE stash_array_mod, ONLY:                                                     &
    len_stlist, stindex, stlist, num_stash_levels, stash_levels, si, sf

IMPLICIT NONE

! Default private
PRIVATE

INTEGER, SAVE            :: offline_update_freq(max_offline_files)
! update freq. (hours)
INTEGER, SAVE            :: offline_update_type(max_offline_files)
! 1 = serial, 2 = periodic
CHARACTER(LEN=256), SAVE :: offline_filename(max_offline_files)
! names of source files

INTEGER, SAVE            :: num_file_namel            ! No. of files in namelist
INTEGER, SAVE            :: num_cdf_offline_flds      ! No. of NetCDF fields

! Super array of offline oxidants
TYPE(ukca_em_struct), SAVE, ALLOCATABLE :: offline_ox(:)

PUBLIC :: ukca_read_offline_oxidants_ctl
PUBLIC :: ukca_offline_oxidants_diags

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='UKCA_READ_OFFLINE_OXIDANTS_MOD'

CONTAINS

! ##############################################################################

SUBROUTINE ukca_read_offline_oxidants_ctl(row_length, rows, model_levels,      &
  global_row_length, global_rows,                                              &
  timestep_number, timestep, l_first)
!---------------------------------------------------------------------------
! Description:
!   Calls routines to intialise the offline oxidants structure, update the
!    data fields when necessary, and fill the aapropriate arrays.
!---------------------------------------------------------------------------

USE ereport_mod,        ONLY: ereport
USE umPrintMgr,         ONLY: umPrint, umMessage, PrintStatus, PrStatus_Diag, &
  PrStatus_Normal

IMPLICIT NONE

! Input arguments with info on model dimensions
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels
INTEGER, INTENT(IN) :: global_row_length
INTEGER, INTENT(IN) :: global_rows
INTEGER, INTENT(IN) :: timestep_number     ! No. timesteps since start

REAL,    INTENT(IN) :: timestep            ! Timestep length (sec)

LOGICAL, INTENT(IN) :: l_first             ! T if first call to UKCA

INTEGER            :: errcode    ! Error code for ereport
CHARACTER (LEN=errormessagelength) :: cmessage   ! Error message

REAL (KIND=jprb)   :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_READ_OFFLINE_OXIDANTS_CTL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in,zhook_handle)

IF (l_first) THEN

  ! Read oxidant NetCDF files and look for the oxidants fields in them
  !  to allocate all variables in the oxidants structure.

  CALL ukca_offline_oxidants_init(row_length, rows, model_levels,              &
    global_row_length, global_rows)
END IF   ! l_first

! Check if it is time to update the fields that are read from NetCDF files
! (depending on time step and update frequency). If needed update the files.

CALL ukca_offline_oxidants_update(row_length, rows, model_levels,              &
  global_row_length, global_rows, timestep_number,                             &
  timestep, l_first)

! Fill the offline oxidant arrays
CALL ukca_fill_offline_oxidants(row_length,rows,model_levels)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_read_offline_oxidants_ctl

! ##############################################################################

SUBROUTINE ukca_offline_oxidants_init(row_length, rows, model_levels,          &
  global_row_length, global_rows)
!---------------------------------------------------------------------------
! Description:
!   Allocates the offline_ox structure and reads attributes from the list
!    of files supplied. Checks that there are no missing fields, and that all
!    fields relate to species defined in the chemistry scheme.
!---------------------------------------------------------------------------

USE asad_mod,         ONLY: speci, ctype, jpcf, jpsp, jpspec
USE um_types,         ONLY: integer32
USE ereport_mod,      ONLY: ereport
USE umPrintMgr,       ONLY: umPrint, umMessage, PrStatus_Normal, PrintStatus
USE emiss_io_mod,     ONLY: em_fclose, em_fopen, em_var_check_dims,           &
                            em_get_file_info, em_get_var_info, em_get_var_att

IMPLICIT NONE

! Subroutine arguments
! model dimensions
INTEGER,      INTENT(IN) :: row_length
INTEGER,      INTENT(IN) :: rows
INTEGER,      INTENT(IN) :: model_levels
INTEGER,      INTENT(IN) :: global_row_length
INTEGER,      INTENT(IN) :: global_rows

! Local
INTEGER, PARAMETER       :: max_dims = 5     ! Max no of dims per var
INTEGER                  :: nvars, varid
INTEGER                  :: dimen_size(max_dims)
INTEGER                  :: ierror           ! Error code
INTEGER                  :: l                ! variable id (for use in loop)
INTEGER                  :: n                ! loop counter for files
INTEGER                  :: ispec            ! loop counter for species
INTEGER                  :: ecount           ! elements counter
INTEGER                  :: fileid           ! File id

! Flag to check presence of an attrib
LOGICAL                  :: l_exist  
! Flags for presence of standard_name & long_name attributes
LOGICAL                  :: l_exist_sn, l_exist_ln 

! both fields are read from global attrib.
CHARACTER (LEN=errormessagelength) :: cmessage              ! error message
CHARACTER (LEN=256) :: variable_name         ! Names to contain file attributes
CHARACTER (LEN=256) :: standard_name
CHARACTER (LEN=256) :: long_name
CHARACTER (LEN=10)  :: offline_tracer
CHARACTER (LEN=30)  :: offline_units
CHARACTER (LEN=256) :: att_name

LOGICAL             :: l_found               ! Indicates if species found

REAL (KIND=jprb)    :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_OFFLINE_OXIDANTS_INIT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)
ierror   = 0

! Initialise variables to default values before
! reading them from namelist
offline_filename    (:) = ' '       ! File names
offline_update_freq (:) = imdi      ! Update frequencies
offline_update_type (:) = imdi      ! Update type: serial, periodic, ...

! Go through all possible files in the namelist and compose dir+filename.
!  Exit on encountering first empty file.
num_file_namel = 0
DO n = 1, max_offline_files
  ! Exit when first empty file found (note that the array ukca_offline_files
  !  was initialised as 'ukca_offline_files is unset' in UKCA_OPTION_MOD)
  IF (ukca_offline_files(n) == ' '                       .OR.                  &
    ukca_offline_files(n) == 'ukca_offline_files is unset') EXIT

  offline_filename(n) = TRIM(ukca_offline_dir)//'/'//                          &
    TRIM(ukca_offline_files(n))
  num_file_namel    = num_file_namel + 1
END DO

! Stop if directory is not set
IF (ukca_offline_dir == 'ukca_offline_dir is unset') THEN
  ierror   = 1
  cmessage = ' Directory for offline oxidant netCDF files is unset'
  CALL ereport('UKCA_OFFLINE_OXIDANTS_INIT', ierror, cmessage)
END IF

! Stop with error message if no files found
IF (num_file_namel == 0) THEN
  ierror   = 1
  cmessage = 'No offline oxidant netCDF files found'
  CALL ereport('UKCA_OFFLINE_OXIDANTS_INIT', ierror, cmessage)
END IF

IF (PrintStatus >= PrStatus_Normal) THEN
  WRITE (umMessage,'(I3,1X,A)') num_file_namel,'files found in offline '//     &
    'namelist'
  CALL umPrint(umMessage,src='ukca_offline_oxidants_init')
END IF

! 1: --------------------------------------------------------------------------
! First loop through all files to get their update
!  frequencies (stored as global attribute) as well as the
!  total number of NetCDF oxidant fields (num_cdf_offline_flds).
num_cdf_offline_flds  = 0

DO n = 1, num_file_namel
  CALL em_fopen(TRIM(offline_filename(n)),fileid)
  
  ! Get the update frequency of the file (in hours). Needs to be read
  !  as character before storing it in a numeric array.
  CALL em_get_var_att (fileid, 1, 'GLOBAL', 'update_freq_in_hours',            &
                       offline_update_freq(n:n))

  ! Make sure that update frequency is filled with an integer greater
  ! or equal to 1
  IF (offline_update_freq(n) < 1) THEN
    ierror   = 1
    cmessage = 'Update frequency needs to be multiple of one hour for ' //     &
      TRIM(offline_filename(n))
    CALL ereport('UKCA_EMISS_INIT', ierror, cmessage)
  END IF

  ! Similarly read the indicator for times at which data is provided.
  !  We adopt the same conventions as for ancillary files:
  !    0 Single time           (not tested and therefore not allowed yet)
  !    1 Time series
  !    2 Periodic time series
  CALL em_get_var_att (fileid, 1, 'GLOBAL', 'update_type',                     &
                       offline_update_type(n:n))

  ! For the moment only allow the two options that have been tested:
  !   time series - update_type = 1
  !   periodic    - update_type = 2
  IF (offline_update_type(n) /= 1 .AND. offline_update_type(n) /= 2) THEN
    ierror   = 1
    cmessage = 'Time series (1) and periodic (2) are the only ' //             &
      'options allowed for ' // TRIM(offline_filename(n))
    CALL ereport('UKCA_OFFLINE_OXIDANTS_INIT', ierror, cmessage)
  END IF

  ! Get number of variables in an offline file
  CALL em_get_file_info(fileid, num_vars = nvars)

  ! Go through all variables in each offline file to calculate the
  ! total number of offline fields (considering all files)
  DO varid = 1, nvars
    ! Set variable_name to ' ' so that EM_GET_VAR_INFO returns
    ! the variable name associated with the variable ID.
    variable_name = ' '
    l = varid     ! cannot directly pass loop counter as arg declared INOUT
    CALL em_get_var_info (fileid, l, variable_name)

    ! Some variables (e.g dimension related) will not have all the standard 
    ! attributes required for offline oxidants. To prevent NetCDF errors for 
    ! missing attributes, first check  whether the variable has attribute 
    ! 'tracer_name', else continue to the next variable.
    ! This is done by setting the flag l_exist to .FALSE; which will cause
    ! em_get_var_att to return .TRUE. if attribute exists, otherwise .FALSE.
    l_exist = .FALSE.
    CALL em_get_var_att (fileid, 1, variable_name, 'tracer_name',              &
                         offline_tracer, l_exist)

    IF ( .NOT. l_exist ) CYCLE ! 'tracer_name' not found, continue to next varid

    ! Select only variables which are consistent with the definition of the
    ! chemistry scheme. For that the attribute 'tracer_name' needs to
    ! be present in the NetCDF field and the variable where that attribute is
    ! stored (i.e. offline_tracer) needs to be equal to one of the elements in
    ! the species name array (i.e. array speci), and that the species is one to
    ! be externally defined (constant field). H2O2 is treated as a special
    ! case as it is a prognostic with an offline field provided as a limiter.
    l_found = .FALSE.
    DO ispec=1,jpspec
      IF ((TRIM(speci(ispec)) == TRIM(offline_tracer) .AND.                    &
        ctype(ispec) == jpcf) .OR.                                             &
        (TRIM(speci(ispec)) == TRIM(offline_tracer) .AND.                      &
        TRIM(offline_tracer) == 'H2O2' .AND.                                   &
        ctype(ispec) == jpsp)) THEN
        num_cdf_offline_flds = num_cdf_offline_flds + 1
        IF (PrintStatus >= PrStatus_Normal) THEN
          WRITE (umMessage,'(A,A,A)') offline_tracer, ' found in ',            &
            TRIM(offline_filename(n))
          CALL umPrint(umMessage,src='ukca_offline_oxidants_init')
        END IF
        l_found = .TRUE.
        EXIT
      END IF
    END DO

    IF (.NOT. l_found) THEN
      ! Report a warning and ignore the field if it is not present in the
      !   current chemistry scheme, or if it is not defined as a constant
      !   field.  Note that with the condition below
      !   the warning message will not appear for variables which are not
      !   offline fields (e.g. dimension variables, variables related to the
      !   "grid_mapping" of the NetCDF files, etc.), because they will not
      !   contain the attribute 'tracer_name'.
      IF (TRIM(offline_tracer) /= '') THEN
        ierror   = -1
        cmessage = 'offline field ' // TRIM(offline_tracer) //                 &
          ' not present in the chemistry scheme and therefore ignored'
        CALL ereport('UKCA_OFFLINE_OXIDANTS_INIT', ierror, cmessage)
      END IF
    END IF       ! not l_found

  END DO         ! varid loop through variables in an emission file

  CALL em_fclose(fileid)
END DO    ! loop through all NetCDF offline files

IF (PrintStatus >= PrStatus_Normal) THEN
  WRITE (umMessage,'(A,1X,I3)') 'Number of NetCDF offline fields:',            &
    num_cdf_offline_flds
  CALL umPrint(umMessage,src='ukca_offline_oxidants_init')
END IF

! Allocate the offline_ox structure and initialise some
!  of its elements to default values
ALLOCATE (offline_ox(num_cdf_offline_flds))
CALL ukca_offline_struct_init()


! 2: --------------------------------------------------------------------------
! Now repeat loops through emission files and fields, to fill in (or allocate)
!  some items in the offline_ox structure

ecount = 0   ! initialise counter for elements in emissions structure
DO n = 1, num_file_namel
  CALL em_fopen(TRIM(offline_filename(n)),fileid)
 
  ! Get number of variables in each emission file
  CALL em_get_file_info(fileid, num_vars = nvars)

  ! Go through all variables in each offline file to calculate the
  ! total number of offline fields (considering all files)
  DO varid = 1, nvars
    ! Set variable_name to ' ' so that EM_GET_VAR_INFO returns
    ! the variable name associated with the variable ID.
    variable_name = ' '
    l = varid     ! cannot directly pass loop counter as arg declared INOUT
    CALL em_get_var_info (fileid, l, variable_name)

    ! As before, do not look for attributes for dimension related variables
    ! (this will be the case if attribute 'tracer_name' is absent)
    l_exist = .FALSE.
    CALL em_get_var_att (fileid, 1, variable_name, 'tracer_name',             &
                         offline_tracer, l_exist)    

    IF ( .NOT. l_exist ) CYCLE ! 'tracer_name' not found, continue to next varid

    ! If an offline field is found then do some updates
    ! to the offline_ox structure
    DO ispec=1,jpspec
      IF ((TRIM(speci(ispec)) == TRIM(offline_tracer) .AND.                   &
        ctype(ispec) == jpcf) .OR.                                            &
        (TRIM(offline_tracer) == 'H2O2' .AND. ctype(ispec) == jpsp)) THEN
        ecount = ecount + 1

        ! Get values of some attributes for the offline field

        ! First read mandatory attributes from this variable
        CALL em_get_var_att (fileid, 1, variable_name, 'units', offline_units)

        ! Read some optional attributes from the variable. This is done by
        ! passing the INOUT argument l_exist_x, which has input value = .FALSE.
        ! The output value will be changed to .TRUE. only if att exists.
        ! Atts which are missing will remain with the initialised value of ' '.
        ! However at least one of two atts (standard_name and long_name) 
        ! should be present to comply with CF conventions.
        standard_name = ' '
        long_name     = ' '

        l_exist_sn = .FALSE.
        CALL em_get_var_att (fileid, 1, variable_name, 'standard_name',       &
                                        standard_name, l_exist_sn)

        l_exist_ln = .FALSE.
        CALL em_get_var_att (fileid, 1, variable_name, 'long_name',           &
                                        long_name,      l_exist_ln)

        IF ( (.NOT. l_exist_sn) .AND. (.NOT. l_exist_ln) ) THEN
          ierror   = 1
          cmessage = 'standard_name and long_name missing in '      //    &
                     TRIM (variable_name)                           //    &
                     ': at least one of them needs to be present'
          CALL ereport('UKCA_EM_STRUCT_SET', ierror, cmessage)
        END IF

        ! Fill in file_name and update_freq with information read from the namelist
        !  -these two variables are the same for all fields within the same file
        offline_ox(ecount)%file_name   = offline_filename(n)
        offline_ox(ecount)%update_freq = offline_update_freq(n) ! hours
        offline_ox(ecount)%update_type = offline_update_type(n)

        ! The other elements of the structure can be read from the variables/attributes
        !  in the files
        offline_ox(ecount)%var_name    = variable_name
        offline_ox(ecount)%tracer_name = offline_tracer
        offline_ox(ecount)%std_name    = standard_name
        offline_ox(ecount)%lng_name    = long_name
        offline_ox(ecount)%units       = offline_units

        ! Check number of dimensions of the field
        CALL em_var_check_dims(global_row_length, global_rows, model_levels,   &
          fileid, variable_name, dimen_size)

        IF (dimen_size(3) == model_levels) THEN   ! 3D field
          ALLOCATE (offline_ox(ecount)%values(row_length, rows, model_levels))
          ALLOCATE (offline_ox(ecount)%diags (row_length, rows, model_levels))
          ! Update the lowest and highest levels where emiss are injected
          !  (just for consistency) and indicate that it is a 3D emiss field
          offline_ox(ecount)%lowest_lev  = 1
          offline_ox(ecount)%highest_lev = model_levels
          offline_ox(ecount)%three_dim   = .TRUE.
        ELSE
          ierror   = fileid
          cmessage = ' Dimension size (3) not equal to model_levels'
          CALL ereport('UKCA_OFFLINE_OXIDANTS_INIT', ierror, cmessage)
        END IF

        offline_ox(ecount)%values(:,:,:) = 0.0  ! will be filled when
        offline_ox(ecount)%diags (:,:,:) = 0.0  ! emissions are updated

        EXIT
      END IF
    END DO       ! ispec=1,jpspec
  END DO       ! varid = 1, nvars

  CALL em_fclose(fileid)
END DO       ! loop through all NetCDF emission fields


! 3: --------------------------------------------------------------------------
! Make sure that there are no missing fields in the offline oxidant files,
!  i.e. all names in speci identified as 'CF' are present in the structure
!  Water vapour is excluded as this can be updated internally
DO ispec = 1, jpspec
  IF (ctype(ispec) == jpcf) THEN
    IF (.NOT. ANY (offline_ox(:)%tracer_name == speci(ispec)) .AND.            &
      speci(ispec) /= 'H2O       ') THEN
      ierror   = ispec
      cmessage = TRIM(speci(ispec)) // ' missing in the emission files'
      CALL ereport ('UKCA_OFFLINE_OXIDANTS_INIT', ierror, cmessage)
    END IF
  END IF
END DO


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_offline_oxidants_init

! ##############################################################################

SUBROUTINE ukca_offline_struct_init
!---------------------------------------------------------------------------
! Description:
!   Initialisation of offline_ox structure to default values.
! Method:
!   Initialise elements of the offline_ox structure, which has been
!   previously allocated with upper bound equal to 'num_cdf_offline_flds'.
!   Note that as indicated below some elements cannot be
!   allocated yet within this subroutine.
!---------------------------------------------------------------------------

USE ereport_mod,        ONLY: ereport
IMPLICIT NONE

INTEGER                        :: ierror         ! Error code
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_OFFLINE_STRUCT_INIT'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

ierror = 0

IF ( .NOT. ALLOCATED(offline_ox)) THEN
  ierror = 1
  CALL ereport('UKCA_OFFLINE_STRUCT_INIT', ierror,                             &
    ' offline_ox super array should be allocated before this routine')
END IF

offline_ox(:)%file_name   = ' '
offline_ox(:)%var_name    = ' '
offline_ox(:)%tracer_name = ' '
offline_ox(:)%std_name    = ' '
offline_ox(:)%lng_name    = ' '
offline_ox(:)%units       = ' '

offline_ox(:)%update_freq = 0
offline_ox(:)%update_type = 0
offline_ox(:)%l_update    = .FALSE.
offline_ox(:)%three_dim   = .FALSE.

offline_ox(:)%base_fact   = 1.0 ! by default do not apply any conversion factor

offline_ox(:)%hourly_fact = ' '
offline_ox(:)%daily_fact  = ' '
offline_ox(:)%vert_fact   = ' '
offline_ox(:)%lowest_lev  = 1   ! Initially assume surface (changed later
offline_ox(:)%highest_lev = 1   ! in the code for both 3D oxidants)

! Note that the pointers offline_ox(:)%values and offline_ox(:)%diags
! cannot be initialised yet. Before that they need to be allocated
! independently for each element of offline_ox(:). This is done in
! ukca_offline_oxidants_init after calling this routine.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_offline_struct_init

! #############################################################################

SUBROUTINE ukca_offline_oxidants_update(row_length, rows, model_levels,        &
  global_row_length, global_rows, timestep_number,                             &
  timestep, l_first)
! ---------------------------------------------------------------------
! Description:
!   Check if it is time to update the values in the offline_ox structure,
!   and do it when needed.
! Method:
! 1) Check first if full hour. If so go through all oxidant fields and
!    decide to update their values only when the hours from start of the
!    run are multiple of their update frequency (in hours).
! 2) When it is time to update then take the current model time
!    (i_year, i_month, i_day, i_hour, i_minute, i_second) from model_time_mod,
!    and pass it to GET_EMFILE_REC to get the previous and next register.
!    Then call EM_GET_DATA to get the emission values for those two registers,
!    and interpolate the values if none of the registers coincides exactly
!    with the model time.
! ---------------------------------------------------------------------

USE ukca_emiss_mod,  ONLY: get_emfile_rec
USE emiss_io_mod,    ONLY: em_fclose, em_fopen, em_var_check_dims, &
  em_get_data, iupd_single
USE umPrintMgr,      ONLY: umPrint, umMessage, PrintStatus, PrStatus_Diag,    &
  Prstatus_Oper
USE ereport_mod,     ONLY: ereport
USE t_int_mod,       ONLY: t_int

USE model_time_mod, ONLY: &
    i_day, i_hour, i_minute, i_month, i_second, i_year

IMPLICIT NONE

! Subroutine arguments
! model dimensions
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels
INTEGER, INTENT(IN) :: global_row_length
INTEGER, INTENT(IN) :: global_rows
INTEGER, INTENT(IN) :: timestep_number  ! No of timesteps since start

REAL,    INTENT(IN) :: timestep         ! Timestep length (sec)

LOGICAL, INTENT(IN) :: l_first          ! True if first call

! Local variables
INTEGER :: ancil_ref_days, ancil_ref_secs  ! Ancil reference time in days and
                                           ! seconds since 0,0 basis time

INTEGER :: i_year_trg, i_month_trg, i_day_trg      ! Target interpolation time:
INTEGER :: i_hour_trg, i_minute_trg, i_second_trg  ! year, month, etc

INTEGER :: day, sec          ! Current time in days/seconds since ancil ref time
REAL    :: hours_from_ancref ! Hours from ancil ref time to current time
REAL    :: fhr               ! Hours from ancil ref time to mid-point of update
                             ! interval
INTEGER :: fdaynum           ! Day number returned by sec2time (not used)
INTEGER :: anc_steps         ! Number of ancil update intervals between ancil
                             ! ref time and current time

INTEGER             :: fid              ! NetCDF file id

INTEGER             :: ireco(2)      ! NetCDF time record indices: prev & next
INTEGER             :: creco(2,6)    ! Time records: 1st dim is for prev & next
! 2nd dim is for yr, mon, day, hr, min, sec

REAL                :: ftreco(2)     ! Time of prev/next time records (fract
! hours since 0:00 h on 1st Jan)

REAL                :: fhournow      ! Current time (fractional hours since
! 0:00 h on 1st Jan)

REAL, PARAMETER     :: fsec_to_hour = 1.0/3600.0  ! Conv factor, seconds to hour
INTEGER, PARAMETER  :: max_dims = 5            ! Max no of dims per var
INTEGER             :: dimsize(max_dims)       ! array of sizes
REAL, ALLOCATABLE   :: tmpdat1 (:,:,:)         ! Temporal data arrays
REAL, ALLOCATABLE   :: tmpdat2 (:,:,:)
INTEGER             :: l                       ! Index
INTEGER             :: ierror                  ! Error code
LOGICAL             :: l_noreco                ! True if record not found (EOF)
CHARACTER (LEN=errormessagelength) :: cmessage                ! error message

REAL (KIND=jprb)    :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_OFFLINE_OXIDANTS_UPDATE'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Reset l_update to false for all fields.
!  Note that this is changed to TRUE below separately for each
!  NetCDF field if it is updated in the current timestep.
offline_ox(:)%l_update = .FALSE.

! Get ancil reference time as days/secs since 0/0 basis time and use this to
! calculate hours from ancil reference time to now.
CALL time2sec(ancil_reftime(1), ancil_reftime(2), ancil_reftime(3), &
              ancil_reftime(4), ancil_reftime(5), ancil_reftime(6), &
              0,0,ancil_ref_days,ancil_ref_secs,lcal360)

CALL time2sec(i_year,i_month,i_day,i_hour,i_minute,i_second,         &
              ancil_ref_days,ancil_ref_secs,day,sec ,lcal360)
hours_from_ancref = day*24.0 + sec*fsec_to_hour

! Loop over NetCDF emission fields
DO l = 1, num_cdf_offline_flds

  ! Decide if updating, and set target time for time-interpolation
  IF (offline_ox(l)%update_type == iupd_single) THEN
    ! Single-time field read on first timestep and no other
    ! Set target time to now for get_emfile_rec (not actually used)
    IF (l_first) THEN
      offline_ox(l)%l_update = .TRUE.
      i_year_trg = i_year
      i_month_trg = i_month
      i_day_trg = i_day
      i_hour_trg = i_hour
      i_minute_trg = i_minute
      i_second_trg = i_second
    END IF

  ELSE

    ! Period or timeseries offline_ox.
    ! Compute integer number of ancil steps from ref time to the begining of the
    ! present update period (=anc_steps). Use FLOOR to round down (even if 
    ! hours_from_ancref is negative) to get time which is at or before now.
    ! Force update when anc_steps increases from that stored at last update.
    ! Always update on first timestep.
    ! Set target time for time-interp to be mid-point of update interval:
    ! start of this update interval plus half the update frequency
    anc_steps = FLOOR(hours_from_ancref / offline_ox(l)%update_freq)

    IF (l_first .OR. anc_steps /= offline_ox(l)%last_update) THEN
      offline_ox(l)%l_update = .TRUE.
      offline_ox(l)%last_update = anc_steps

      fhr = (REAL(anc_steps) + 0.5) * REAL(offline_ox(l)%update_freq)
      day = INT(fhr/ 24.0)
      sec = INT((fhr - day*24) * 3600.0)
      CALL sec2time (day, sec, ancil_ref_days, ancil_ref_secs, &
        i_year_trg, i_month_trg, i_day_trg, &
        i_hour_trg, i_minute_trg, i_second_trg, &
        fdaynum, lcal360 )
    END IF
  END IF


  IF (offline_ox(l)%l_update) THEN

    ! Debug update times only for first field
    IF (l == 1 .AND. PrintStatus >= PrStatus_Diag) THEN
      WRITE (umMessage, '(A,6(1x,I4))')                                        &
       'UKCA_OFFLINE_OXIDANTS_UPDATE: Update first data field for target time',&
       i_year_trg, i_month_trg, i_day_trg,                                     &
       i_hour_trg, i_minute_trg, i_second_trg
      CALL umPrint(umMessage,src='ukca_offline_oxidants_update')
    END IF

    ! Open file to read an oxidant field
    CALL em_fopen(offline_ox(l)%file_name, fid)

    ! Pass current model time (i_year, ..., i_second) and other data to
    !  GET_EMFILE_REC in order to get:
    !  * current hour since basis time: fhournow,
    !  * indices of the matching NetCDF records used for
    !     interpolation: ireco (1:2)
    !  * yr, mon, day, hr, min, sec for those two records: creco (1:2,1:6)
    !  * hour since basis time for those records: ftreco (1:2)

    CALL get_emfile_rec(row_length, rows, fid, l,                            &
      offline_ox(l)%update_type,                                             &
      (/i_year_trg, i_month_trg, i_day_trg,                                  &
        i_hour_trg, i_minute_trg, i_second_trg /),                           &
      l_first, ireco, creco, fhournow, ftreco)

    ! Later we will call EM_GET_DATA to get the offline_ox values for registers
    !  ireco(1) & ireco(2), at hours ftreco(1) and ftreco(2), so that we can
    !  interpolate them to the current model time in hours (i.e. to hournow).

    ! Check that matching time records are found
    IF ( ANY(ireco(:) == -1) ) THEN  
      WRITE (umMessage,'(A,6(1x,I4),1x,F15.5)') 'At time:', i_year,          &
        i_month, i_day, i_hour, i_minute, i_second, fhournow
      CALL umPrint(umMessage,src='ukca_emiss_init')
      CALL em_fclose(fid)
      ierror = l
      cmessage =                                                             &
        'Required time record not in file ' // offline_ox(l)%file_name
      CALL ereport('UKCA_OFFLINE_OXIDANTS_UPDATE', ierror, cmessage)
    END IF

    ! Go ahead with reading values and interpolations

    ! First get dimensions, in particular the vertical levels, so
    !  that we can allocate the temporary data arrays tmpdat1 and tmpdat2.
    CALL em_var_check_dims(global_row_length, global_rows, model_levels,     &
      fid, offline_ox(l)%var_name, dimsize)

    l_noreco = .FALSE.  ! initially assume that records are found

    IF ( ireco(1) == ireco(2) ) THEN
      ! No interpolation if exact time match. Get only emiss values for
      !  register ireco(1).
      CALL em_get_data (row_length, rows, fid, ireco(1),                     &
        offline_ox(l)%var_name,                                              &
        offline_ox(l)%values(:, :, 1:dimsize(3)), l_noreco)
      IF (l_noreco) THEN  
        CALL em_fclose(fid)
        ierror = l
        cmessage =  'Error reading data1: ' // offline_ox(l)%var_name
        CALL ereport('UKCA_OFFLINE_OXIDANTS_UPDATE', ierror, cmessage)
      END IF

    ELSE
      ! Get two values tmpdat1 and tmpdat2, for registers
      !  ireco(1) and ireco(2), to later interpolate.
      ALLOCATE(tmpdat1(row_length, rows, dimsize(3)) )
      ALLOCATE(tmpdat2(row_length, rows, dimsize(3)) )

      ! Read in data for previous matching time. Stop if record not found.
      CALL em_get_data(row_length, rows, fid, ireco(1),                      &
        offline_ox(l)%var_name,                                              &
        tmpdat1 (:,:,:), l_noreco)
      IF (l_noreco) THEN  
        CALL em_fclose(fid)
        ierror = l
        cmessage = 'Error reading data1: ' // offline_ox(l)%var_name
        CALL ereport('UKCA_OFFLINE_OXIDANTS_UPDATE', ierror, cmessage)
      END IF

      ! Read in data for next matching time.
      CALL em_get_data (row_length, rows, fid, ireco(2),                     &
        offline_ox(l)%var_name,                                              &
        tmpdat2 (:,:,:), l_noreco)
      IF (l_noreco) THEN  
        CALL em_fclose(fid)
        ierror = l
        cmessage = 'Error reading data2: ' // offline_ox(l)%var_name
        CALL ereport('UKCA_OFFLINE_OXIDANTS_UPDATE', ierror, cmessage)
      END IF

      ! Interpolate using same UM routine as REPLANCA
      CALL t_int (                                                           &
        tmpdat1 (:,:,:), ftreco(1),               &  ! Data1 and time1
        tmpdat2 (:,:,:), ftreco(2),               &  ! Data2 and time2
        offline_ox(l)%values(:, :, 1:dimsize(3)), &  ! Data interpolated
        fhournow,                                 &  ! to current time
        row_length*rows*dimsize(3))                  ! No of points

      DEALLOCATE(tmpdat2)
      DEALLOCATE(tmpdat1)
    END IF      ! Interpolation required?

    ! Close emission file after reading an emission field
    CALL em_fclose(fid)

    ! Write diagnostic results
    IF (PrintStatus >= PrStatus_Oper) THEN
      WRITE(umMessage,'(A26,6(1x,I4),F15.5)') 'Updated Oxidant at time: ',   &
        i_year, i_month, i_day, i_hour, i_minute, i_second, fhournow
      CALL umPrint(umMessage,src='ukca_offline_oxidants_update')
      WRITE(umMessage,'(2A15)') 'Oxidant name: ', offline_ox(l)%tracer_name
      CALL umPrint(umMessage,src='ukca_offline_oxidants_update')
      WRITE(umMessage,'(A35,3E12.3)') 'Oxidant level 1: (max, min, mean) ',  &
        MAXVAL(offline_ox(l)%values(:,:,1)),                                 &
        MINVAL(offline_ox(l)%values(:,:,1)),                                 &
        SUM(offline_ox(l)%values(:,:,1))/SIZE(offline_ox(l)%values(:,:,1))
      CALL umPrint(umMessage,src='ukca_offline_oxidants_update')
    END IF

  END IF        ! Update time?
END DO          ! loop over oxidant species


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_offline_oxidants_update

! #############################################################################

SUBROUTINE ukca_fill_offline_oxidants(row_length,rows,model_levels)
! ---------------------------------------------------------------------
! Description:
!  To allocate the offline oxidant arrays, and fill with data from the values
!  read into the offline_ox structure from the netCDF input fields.
! ---------------------------------------------------------------------

USE ukca_chem_offline,    ONLY: o3_offline, oh_offline, no3_offline,           &
  ho2_offline, h2o2_offline
IMPLICIT NONE

! Input arguments with info on model dimensions
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels

INTEGER :: l        ! counter

REAL (KIND=jprb)    :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_FILL_OFFLINE_OXIDANTS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF (.NOT. ALLOCATED(o3_offline))                                               &
  ALLOCATE(o3_offline(row_length,rows,model_levels))
IF (.NOT. ALLOCATED(oh_offline))                                               &
  ALLOCATE(oh_offline(row_length,rows,model_levels))
IF (.NOT. ALLOCATED(no3_offline))                                              &
  ALLOCATE(no3_offline(row_length,rows,model_levels))
IF (.NOT. ALLOCATED(ho2_offline))                                              &
  ALLOCATE(ho2_offline(row_length,rows,model_levels))
IF (.NOT. ALLOCATED(h2o2_offline))                                             &
  ALLOCATE(h2o2_offline(row_length,rows,model_levels))

DO l = 1,num_cdf_offline_flds
  IF (TRIM(offline_ox(l)%tracer_name) == 'O3' .AND.                            &
    offline_ox(l)%l_update) THEN
    o3_offline(:,:,:) = offline_ox(l)%values(:,:,:)
  END IF
  IF (TRIM(offline_ox(l)%tracer_name) == 'OH' .AND.                            &
    offline_ox(l)%l_update) THEN
    oh_offline(:,:,:) = offline_ox(l)%values(:,:,:)
  END IF
  IF (TRIM(offline_ox(l)%tracer_name) == 'NO3' .AND.                           &
    offline_ox(l)%l_update) THEN
    no3_offline(:,:,:) = offline_ox(l)%values(:,:,:)
  END IF
  IF (TRIM(offline_ox(l)%tracer_name) == 'HO2' .AND.                           &
    offline_ox(l)%l_update) THEN
    ho2_offline(:,:,:) = offline_ox(l)%values(:,:,:)
  END IF
  IF (TRIM(offline_ox(l)%tracer_name) == 'H2O2' .AND.                          &
    offline_ox(l)%l_update) THEN
    h2o2_offline(:,:,:) = offline_ox(l)%values(:,:,:)
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_fill_offline_oxidants

! #############################################################################

SUBROUTINE ukca_offline_oxidants_diags(row_length, rows, model_levels,         &
  totnodens,                                                                   &
  len_stashwork,STASHwork50)
! ---------------------------------------------------------------------
! Description:
!  To place offline oxidants fields into the section 50 STASHwork array
! Method:
!  For the oxidants which are constant fields, the diagnostic fields are
!  allocated and filled in routine UKCA_SET_DIURNAL_OX.  For the H2O2
!  limiting field, there is no requirement for a seperate diagnostic field
!  as it is not given a diurnal cycle.
! ---------------------------------------------------------------------

USE ukca_chem_offline,   ONLY: o3_offline_diag, oh_offline_diag,               &
  no3_offline_diag, ho2_offline_diag
USE ukca_chem_offline,   ONLY: h2o2_offline
USE ukca_constants,      ONLY: c_o3, c_oh, c_no3, c_ho2
USE ukca_d1_defs,        ONLY: UKCA_diag_sect
USE ereport_mod,         ONLY: ereport
USE um_parvars,          ONLY: at_extremity
IMPLICIT NONE


! Input arguments
INTEGER, INTENT(IN) :: row_length         ! No. of rows
INTEGER, INTENT(IN) :: rows               ! No. of columns
INTEGER, INTENT(IN) :: model_levels       ! No. of levels
INTEGER, INTENT(IN) :: len_stashwork      ! Length of stashwork array

REAL, INTENT(IN)    :: totnodens(row_length,rows,model_levels)
! Density in molecules/m^3
REAL, INTENT(INOUT) :: STASHwork50(len_stashwork)  ! Work array

! Local variables

INTEGER :: section                        ! Stash section
INTEGER :: item                           ! Stash item
INTEGER :: im_index                       ! internal model index
INTEGER :: icode = 0                      ! error code
INTEGER :: k                              ! loop index
INTEGER :: errcode                        ! Error code for ereport

REAL :: field3d(row_length,rows,model_levels)  ! To contain 3d diagnostic

CHARACTER(LEN=errormessagelength) :: cmessage            ! Error return message

REAL (KIND=jprb)   :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_OFFLINE_OXIDANTS_DIAGS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

im_index = 1

! Offline oxidant diagnostics

section = ukca_diag_sect
DO item = 206,210
  IF (sf(item,section) .AND. item == 206) THEN
    DO k=1,model_levels
      field3d(:,:,k) = RESHAPE(o3_offline_diag(:,k),(/row_length,rows/))*      &
        c_o3/(totnodens(:,:,k)/1.0e6)
    END DO
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork50(si(item,section,im_index)),                  &
      field3d,                                                                 &
      row_length,rows,model_levels,0,0,0,0, at_extremity,                      &
      stlist(1,stindex(1,item,section,im_index)),len_stlist,                   &
      stash_levels,num_stash_levels+1,                                         &
      atmos_im,section,item,icode,cmessage)
    DEALLOCATE(o3_offline_diag)
  ELSE IF (sf(item,section) .AND. item == 207) THEN
    DO k=1,model_levels
      field3d(:,:,k) = RESHAPE(oh_offline_diag(:,k),(/row_length,rows/))*      &
        c_oh/(totnodens(:,:,k)/1.0e6)
    END DO
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork50(si(item,section,im_index)),                  &
      field3d,                                                                 &
      row_length,rows,model_levels,0,0,0,0, at_extremity,                      &
      stlist(1,stindex(1,item,section,im_index)),len_stlist,                   &
      stash_levels,num_stash_levels+1,                                         &
      atmos_im,section,item,icode,cmessage)
    DEALLOCATE(oh_offline_diag)
  ELSE IF (sf(item,section) .AND. item == 208) THEN
    DO k=1,model_levels
      field3d(:,:,k) = RESHAPE(no3_offline_diag(:,k),(/row_length,rows/))*     &
        c_no3/(totnodens(:,:,k)/1.0e6)
    END DO
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork50(si(item,section,im_index)),                  &
      field3d,                                                                 &
      row_length,rows,model_levels,0,0,0,0, at_extremity,                      &
      stlist(1,stindex(1,item,section,im_index)),len_stlist,                   &
      stash_levels,num_stash_levels+1,                                         &
      atmos_im,section,item,icode,cmessage)
    DEALLOCATE(no3_offline_diag)
  ELSE IF (sf(item,section) .AND. item == 209) THEN
    DO k=1,model_levels
      field3d(:,:,k) = RESHAPE(ho2_offline_diag(:,k),(/row_length,rows/))*     &
        c_ho2/(totnodens(:,:,k)/1.0e6)
    END DO
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork50(si(item,section,im_index)),                  &
      field3d,                                                                 &
      row_length,rows,model_levels,0,0,0,0, at_extremity,                      &
      stlist(1,stindex(1,item,section,im_index)),len_stlist,                   &
      stash_levels,num_stash_levels+1,                                         &
      atmos_im,section,item,icode,cmessage)
    DEALLOCATE(ho2_offline_diag)
  ELSE IF (sf(item,section) .AND. item == 210) THEN
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork50(si(item,section,im_index)),                  &
      h2o2_offline(:,:,:),                                                     &
      row_length,rows,model_levels,0,0,0,0, at_extremity,                      &
      stlist(1,stindex(1,item,section,im_index)),len_stlist,                   &
      stash_levels,num_stash_levels+1,                                         &
      atmos_im,section,item,icode,cmessage)
  END IF

  IF (icode >  0) THEN 
    errcode = section*1000+item
    CALL ereport('UKCA_OFFLINE_OXIDANTS_DIAGS',errcode,cmessage)
  END IF
END DO       ! item

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE ukca_offline_oxidants_diags

END MODULE ukca_read_offline_oxidants_mod
