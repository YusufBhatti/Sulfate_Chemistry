! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Module gathering code that reads in and updates the EasyAerosol
!  distributions
!
! Method:
!
!  EasyAerosol files are in netCDF format.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
! Code description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.
!
! Contained subroutines in this module:
!   easyaerosol_read_input_ctl         (Public)
!   easyaerosol_init_input             (Private)
!   easyaerosol_update                 (Private)
!   easyaerosol_fill                   (Private)
!
! --------------------------------------------------------------------------

MODULE easyaerosol_read_input_mod

USE easyaerosol_option_mod,  ONLY: l_easyaerosol_sw, &
                                   l_easyaerosol_lw, &
                                   l_easyaerosol_cdnc, &
                                   l_easyaerosol_autoconv, &
                                   l_easyaerosol_zonal, &
                                   easyaerosol_dir,   &
                                   easyaerosol_files, &
                                   n_easyaerosol_files, &
                                   n_easy
USE missing_data_mod, ONLY: imdi
USE easyaerosol_mod,         ONLY: easyaerosol_struct
USE parkind1,         ONLY: jpim, jprb      ! DrHook
USE yomhook,          ONLY: lhook, dr_hook  ! DrHook
USE errormessagelength_mod, ONLY: errormessagelength

USE nlstcall_mod,          ONLY: ancil_reftime, lcal360
USE conversions_mod,       ONLY: rhour_per_sec

USE filenamelength_mod, ONLY: &
    filenamelength

IMPLICIT NONE

! Default private
PRIVATE

INTEGER, SAVE            :: easyaerosol_update_freq(n_easyaerosol_files)
! update freq. (hours)
INTEGER, SAVE            :: easyaerosol_update_type(n_easyaerosol_files)
! 1 = serial, 2 = periodic
CHARACTER(LEN=filenamelength), SAVE :: easyaerosol_filename(n_easyaerosol_files)
! names of source files

! No. of EasyAerosol distributions
INTEGER, SAVE            :: n_easyaerosol_dist

! Arrays of EasyAerosol distributions in the shortwave, longwave, and CDNC.
TYPE(easyaerosol_struct), SAVE, ALLOCATABLE :: easyaerosol_dist(:)

PUBLIC :: easyaerosol_read_input_ctl

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EASYAEROSOL_READ_INPUT_MOD'

CONTAINS

! ##############################################################################

SUBROUTINE easyaerosol_read_input_ctl(row_length, rows, model_levels,    &
  global_row_length, global_rows, h_swbands, h_lwbands,                  &
  timestep_number, timestep, l_first,                                    &
  easyaerosol_sw, easyaerosol_lw, easyaerosol_cdnc)

!---------------------------------------------------------------------------
! Description:
!   Calls routines to intialise the EasyAerosol structure, update the
!    data fields when necessary, and fill the appropriate arrays.
!---------------------------------------------------------------------------

USE def_easyaerosol,    ONLY: t_easyaerosol_rad, t_easyaerosol_cdnc

IMPLICIT NONE

! Input arguments with info on model dimensions
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels
INTEGER, INTENT(IN) :: global_row_length
INTEGER, INTENT(IN) :: global_rows
INTEGER, INTENT(IN) :: h_swbands
INTEGER, INTENT(IN) :: h_lwbands
INTEGER, INTENT(IN) :: timestep_number     ! No. timesteps since start

REAL,    INTENT(IN) :: timestep            ! Timestep length (sec)

LOGICAL, INTENT(IN) :: l_first             ! T if first call to Atm_Step

TYPE (t_easyaerosol_rad), INTENT(INOUT) :: easyaerosol_sw
TYPE (t_easyaerosol_rad), INTENT(INOUT) :: easyaerosol_lw
TYPE (t_easyaerosol_cdnc), INTENT(INOUT) :: easyaerosol_cdnc

! Local variables
REAL (KIND=jprb)   :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EASYAEROSOL_READ_INPUT_CTL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in,zhook_handle)

IF (l_first) THEN

  ! Read EasyAerosol NetCDF files and look for the expected fields in them
  !  to allocate all variables in the EasyAerosol structure.

  CALL easyaerosol_init_input(row_length, rows, model_levels,                  &
    global_row_length, global_rows, h_swbands, h_lwbands)
END IF   ! l_first

! Check if it is time to update the fields that are read from NetCDF files
! (depending on time step and update frequency). If needed update the files.

IF (l_easyaerosol_zonal) THEN

! For files with zonal mean climatology (no longitude dimension)
  CALL easyaerosol_update_zonal(row_length, rows, model_levels,           &
    global_rows, timestep_number,                                         &
    timestep, l_first)

ELSE

! For files with full 3D climatology (default)
  CALL easyaerosol_update(row_length, rows, model_levels,                 &
    global_row_length, global_rows, timestep_number,                      &
    timestep, l_first)

END IF

! Fill the EasyAerosol arrays
CALL easyaerosol_fill(row_length, rows, model_levels, h_swbands, h_lwbands, &
                      easyaerosol_sw, easyaerosol_lw, easyaerosol_cdnc)


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,zhook_handle)
RETURN
END SUBROUTINE easyaerosol_read_input_ctl

! ##############################################################################

SUBROUTINE easyaerosol_init_input(row_length, rows, model_levels,             &
  global_row_length, global_rows, h_swbands, h_lwbands)
!---------------------------------------------------------------------------
! Description:
!   Allocates the EasyAerosol structure and reads attributes from the list
!    of files supplied. Checks that there are no missing fields.
!---------------------------------------------------------------------------

USE ereport_mod,      ONLY: ereport
USE umPrintMgr,       ONLY: umPrint, umMessage, PrStatus_Normal, PrintStatus
USE emiss_io_mod,     ONLY: em_fclose, em_fopen, em_var_check_dims,           &
                            em_get_var_info, em_get_var_att

IMPLICIT NONE

! Subroutine arguments
! model dimensions
INTEGER,      INTENT(IN) :: row_length
INTEGER,      INTENT(IN) :: rows
INTEGER,      INTENT(IN) :: model_levels
INTEGER,      INTENT(IN) :: global_row_length
INTEGER,      INTENT(IN) :: global_rows
INTEGER,      INTENT(IN) :: h_swbands
INTEGER,      INTENT(IN) :: h_lwbands

! Local
INTEGER, PARAMETER       :: max_dims = 5     ! Max no of dims per var
INTEGER                  :: nvars, varid
INTEGER                  :: dimen_size(max_dims)
INTEGER                  :: ierror           ! Error code
INTEGER                  :: l                ! variable id (for use in loop)
INTEGER                  :: n                ! loop counter for files
INTEGER                  :: ivar             ! loop counter for variables
INTEGER                  :: ecount           ! elements counter
INTEGER                  :: fileid           ! File id

! Names of variables expected in the input files.
! EasyAerosol files contain only one distribution each.
CHARACTER (LEN=18) :: easyaerosol_varname(n_easyaerosol_files) = &
  (/ 'easy_extinction_sw', 'easy_absorption_sw', 'easy_asymmetry_sw ', &
     'easy_extinction_lw', 'easy_absorption_lw', 'easy_asymmetry_lw ', &
     'easy_cdnc         ' /)

! Arrays of flags indicating whether the file/variable, and whether than
! variable has a fourth dimension (typically, number of spectral bands).
! EasyAerosol netCDF files contain one variable each (not including dimension
! and bounds variables).
LOGICAL :: l_easyaerosol_required(n_easyaerosol_files)
INTEGER :: n_easyaerosol_specbands(n_easyaerosol_files)

INTEGER i_easyaerosol_dist

CHARACTER (LEN=4)   :: char_update_freq ! Emiss update freq (in hours),
CHARACTER (LEN=1)   :: char_update_type ! serial, periodic, single point
! both fields are read from global attrib.
CHARACTER (LEN=errormessagelength) :: cmessage              ! error message
CHARACTER (LEN=256) :: variable_name         ! Names to contain file attributes
CHARACTER (LEN=256) :: standard_name
CHARACTER (LEN=256) :: long_name
CHARACTER (LEN=30)  :: easy_units
CHARACTER (LEN=256) :: att_name

LOGICAL             :: l_found               ! Indicates if file/variable found
LOGICAL             :: l_exist_sn
LOGICAL             :: l_exist_ln

REAL (KIND=jprb)    :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EASYAEROSOL_INIT_INPUT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)
ierror   = 0

! Initialise variables to default values before updating them from namelist
easyaerosol_filename    (:) = 'unset'   ! File names
easyaerosol_update_freq (:) = imdi      ! Update frequencies
easyaerosol_update_type (:) = imdi      ! Update type: serial, periodic, ...

! Determine which distributions are required depending on EasyAerosol requests.
n_easyaerosol_dist = 0
l_easyaerosol_required(:) = .FALSE.

IF (l_easyaerosol_sw) THEN
  
  ! SW variables (indices 1 to n_easy in 
  ! easyaerosol_varname) are required.
  ! 4th dimension is the number of shortwave wavebands.
  n_easyaerosol_dist = n_easyaerosol_dist + n_easy
  l_easyaerosol_required(1:n_easy) = .TRUE.
  n_easyaerosol_specbands(1:n_easy) = h_swbands

END IF

IF (l_easyaerosol_lw) THEN
  
  ! LW variables (indices n_easy+1 to 2*n_easy in 
  ! easyaerosol_varname) are required.
  ! 4th dimension is the number of longwave wavebands.
  n_easyaerosol_dist = n_easyaerosol_dist + n_easy
  l_easyaerosol_required(n_easy+1:2*n_easy) = .TRUE.
  n_easyaerosol_specbands(n_easy+1:2*n_easy) = h_lwbands
  

END IF

IF (l_easyaerosol_cdnc .OR. l_easyaerosol_autoconv) THEN
  
  ! CDNC variable (last index in 
  ! easyaerosol_varname) is required.
  ! There is no 4th dimension for that distribution.
  n_easyaerosol_dist = n_easyaerosol_dist + 1
  l_easyaerosol_required(2*n_easy+1) = .TRUE.
  n_easyaerosol_specbands(2*n_easy+1) = 0

END IF

IF (n_easyaerosol_dist == 0) THEN
  ierror   = 1
  cmessage = 'EasyAerosol being initialised but not used.'
  CALL ereport('EASYAEROSOL_INIT_INPUT', ierror, cmessage)
END IF

! Allocate and initiliase EasyAerosol structure
ALLOCATE(easyaerosol_dist(n_easyaerosol_dist))
CALL easyaerosol_struct_init()

! Stop if directory is empty or not set
IF (LEN_TRIM(easyaerosol_dir) == 0 .OR. easyaerosol_dir == 'unset') THEN
  ierror   = 1
  cmessage = ' Directory for EasyAerosol netCDF files is not set'
  CALL ereport('EASYAEROSOL_INIT_INPUT', ierror, cmessage)
END IF

! Go through all required files in the namelist and compose dir+filename.
! Missing at least one required file name is an error.
i_easyaerosol_dist = 1
DO n = 1, n_easyaerosol_files

  IF (l_easyaerosol_required(n)) THEN
  
    IF (LEN_TRIM(easyaerosol_files(n)) == 0 .OR. &
        easyaerosol_files(n) == 'unset') THEN

      ierror   = 1
      WRITE(cmessage, '(A,I0,A)') 'Required EasyAerosol netCDF file number ', &
                                  n, ' is not set'
      CALL ereport('EASYAEROSOL_INIT_INPUT', ierror, cmessage)    
    
    ELSE
    
      easyaerosol_filename(n) = &
        TRIM(easyaerosol_dir) // '/' // TRIM(easyaerosol_files(n))
      easyaerosol_dist(i_easyaerosol_dist)%n_specbands = &
        n_easyaerosol_specbands(n)
  
      ! Get the file's update frequency (stored as global attribute) and 
      ! check that it contains the distribution that we require (EasyAerosol 
      ! files contain only one distribution each).

      CALL em_fopen(easyaerosol_filename(n), fileid)
  
      ! Get the update frequency of the file (in hours). Needs to be read
      !  as character before storing it in a numeric array.
      CALL em_get_var_att(fileid, 1, 'GLOBAL', 'update_freq_in_hours',         &
                          char_update_freq)
      READ(char_update_freq, '(I4)') easyaerosol_update_freq(n)

      ! Make sure that update frequency is filled with an integer greater
      ! or equal to 1
      IF (easyaerosol_update_freq(n) < 1) THEN
        ierror   = 1
        cmessage = 'Update frequency needs to be multiple of one hour for ' // &
                   TRIM(easyaerosol_filename(n))
        CALL ereport('EASYAEROSOL_INIT_INPUT', ierror, cmessage)
      END IF

      ! Similarly read the indicator for times at which data is provided.
      !  We adopt the same conventions as for ancillary files:
      !    0 Single time           (not tested and therefore not allowed yet)
      !    1 Time series
      !    2 Periodic time series
      CALL em_get_var_att(fileid, 1, 'GLOBAL', 'update_type',                  &
                          char_update_type)
      READ(char_update_type, '(I1)') easyaerosol_update_type(n)

      ! For the moment only allow the two options that have been tested:
      !   time series - update_type = 1
      !   periodic    - update_type = 2
      IF (easyaerosol_update_type(n) /= 1 .AND. &
          easyaerosol_update_type(n) /= 2) THEN
        ierror   = 1
        cmessage = 'Time series (1) and periodic (2) are the only ' //         &
                   'options allowed for ' // &
                   TRIM(easyaerosol_filename(n))
        CALL ereport('EASYAEROSOL_INIT_INPUT', ierror, cmessage)
      END IF

      ! Get the variable ID corresponding to the EasyAerosol variable
      ! expected in that file, and store it.
      CALL em_get_var_info(fileid, l, easyaerosol_varname(n))
      easyaerosol_dist(i_easyaerosol_dist)%var_name = easyaerosol_varname(n)
      easyaerosol_dist(i_easyaerosol_dist)%varid = l
      
      IF (PrintStatus >= PrStatus_Normal) THEN
        WRITE (umMessage,'(A,A,A,A)') 'EasyAerosol variable ', &
            TRIM(easyaerosol_varname(n)), ' found in ',        &
            TRIM(easyaerosol_filename(n))
        CALL umPrint(umMessage,src='easyaerosol_init_input')
      END IF

      ! Get values of some attributes for the EasyAerosol field

      ! First read mandatory attributes from this variable
      CALL em_get_var_att(fileid, 1, easyaerosol_varname(n), &
                          'units', easy_units)

      ! Read some optional attributes from the variable. This is done by
      ! passing the INOUT argument l_exist_x, which has input value = .FALSE.
      ! The output value will be changed to .TRUE. only if att exists.
      ! Atts which are missing will remain with the initialised value of ' '.
      ! However at least one of two atts (standard_name and long_name) 
      ! should be present to comply with CF conventions.
      standard_name = ' '
      long_name     = ' '

      l_exist_sn = .FALSE.
      CALL em_get_var_att(fileid, 1, easyaerosol_varname(n), 'standard_name', &
                                     standard_name, l_exist_sn)

      l_exist_ln = .FALSE.
      CALL em_get_var_att(fileid, 1, easyaerosol_varname(n), 'long_name', &
                                     long_name, l_exist_ln)

      IF ( (.NOT. l_exist_sn) .AND. (.NOT. l_exist_ln) ) THEN
        ierror   = 1
        cmessage = 'standard_name and long_name missing in '      //    &
                   TRIM(easyaerosol_varname(n))                   //    &
                   ': at least one of them needs to be present'
        CALL ereport('EASYAEROSOL_INIT_INPUT', ierror, cmessage)
      END IF

      ! Update elements of the structure that were read from the 
      ! variables/attributes in the file.
      easyaerosol_dist(i_easyaerosol_dist)%std_name    = standard_name
      easyaerosol_dist(i_easyaerosol_dist)%long_name   = long_name
      easyaerosol_dist(i_easyaerosol_dist)%units       = easy_units
      
      ! Assuming that data files contain full 3D spatial fields
      ! (default) rather than zonal means
      ! 
      IF (.NOT. l_easyaerosol_zonal) THEN
        ! Some distributions have additional dimension for spectral bands
        ! 
        IF (easyaerosol_dist(i_easyaerosol_dist)%n_specbands /= 0) THEN
      
          ! Check number of dimensions of the field
          CALL em_var_check_dims(global_row_length, global_rows, model_levels, &
                                 n_easyaerosol_specbands(n),                   &
                                 fileid, easyaerosol_varname(n), dimen_size)
 
          ALLOCATE(easyaerosol_dist(i_easyaerosol_dist)%values_4d(row_length, &
                   rows, model_levels, n_easyaerosol_specbands(n)))
          easyaerosol_dist(i_easyaerosol_dist)%values_4d(:,:,:,:) = 0.0
        
          easyaerosol_dist(i_easyaerosol_dist)%ndims = 4

        ELSE
        
          ! Check number of dimensions of the field
          CALL em_var_check_dims(global_row_length, global_rows, model_levels, &
                                 fileid, easyaerosol_varname(n), dimen_size)

          ALLOCATE(easyaerosol_dist(i_easyaerosol_dist)%values_3d(row_length, &
                   rows, model_levels))
          easyaerosol_dist(i_easyaerosol_dist)%values_3d(:,:,:) = 0.0

          easyaerosol_dist(i_easyaerosol_dist)%ndims = 3
 
        END IF

      ! Now consider the case that the data file contains zonal mean fields
      ! so there is no longitude dimension
      ELSE ! (l_easyaerosol_zonal = .true.)

        ! Some distributions have additional dimension for spectral bands
        ! 
        IF (easyaerosol_dist(i_easyaerosol_dist)%n_specbands /= 0) THEN

          ! Check number of dimensions of the field
          CALL em_var_check_dims(global_rows, model_levels,  &
                                 n_easyaerosol_specbands(n), &
                                 l_easyaerosol_zonal,        &
                                 fileid, easyaerosol_varname(n), dimen_size)
 
          ALLOCATE(easyaerosol_dist(i_easyaerosol_dist)%values_4d(&
              row_length, rows, model_levels, n_easyaerosol_specbands(n)))
          easyaerosol_dist(i_easyaerosol_dist)%values_4d(:,:,:,:) = 0.0
        
          easyaerosol_dist(i_easyaerosol_dist)%ndims = 4

        ELSE
        
          ! Check number of dimensions of the field
          CALL em_var_check_dims(global_rows, model_levels,  &
                                 l_easyaerosol_zonal,        &
                               fileid, easyaerosol_varname(n), dimen_size)

          ALLOCATE(easyaerosol_dist(i_easyaerosol_dist)%values_3d( &
               row_length, rows, model_levels))
          easyaerosol_dist(i_easyaerosol_dist)%values_3d(:,:,:) = 0.0

          easyaerosol_dist(i_easyaerosol_dist)%ndims = 3

        END IF
        
      END IF      

      CALL em_fclose(fileid)
      
      !
      ! Fill in file_name, update frequency and type with information
      ! from the namelist.
      !
      easyaerosol_dist(i_easyaerosol_dist)%file_name   = &
        easyaerosol_filename(n)
      easyaerosol_dist(i_easyaerosol_dist)%update_freq = &
        easyaerosol_update_freq(n)
      easyaerosol_dist(i_easyaerosol_dist)%update_type = &
        easyaerosol_update_type(n)

      IF (PrintStatus >= PrStatus_Normal) THEN

        WRITE (umMessage,'(A,I0,A)') 'EasyAerosol distribution number ', &
            i_easyaerosol_dist, ' initialised.'
        CALL umPrint(umMessage,src='easyaerosol_init_input')

        WRITE (umMessage,'(A,A)') 'Source file ', &
            TRIM(easyaerosol_dist(i_easyaerosol_dist)%file_name)
        CALL umPrint(umMessage,src='easyaerosol_init_input')

        WRITE (umMessage,'(A,A)') 'Source variable name ', &
            TRIM(easyaerosol_dist(i_easyaerosol_dist)%var_name)
        CALL umPrint(umMessage,src='easyaerosol_init_input')

        WRITE (umMessage,'(A,A)') 'Source standard name ', &
            TRIM(easyaerosol_dist(i_easyaerosol_dist)%std_name)
        CALL umPrint(umMessage,src='easyaerosol_init_input')

        WRITE (umMessage,'(A,A)') 'Source long name ', &
            TRIM(easyaerosol_dist(i_easyaerosol_dist)%long_name)
        CALL umPrint(umMessage,src='easyaerosol_init_input')

        WRITE (umMessage,'(A,A)') 'Source variable units ', &
            TRIM(easyaerosol_dist(i_easyaerosol_dist)%units)
        CALL umPrint(umMessage,src='easyaerosol_init_input')

        WRITE (umMessage,'(A,3(1X,I0))') &
            'Variable ID, ndims, n_specbands', &
            easyaerosol_dist(i_easyaerosol_dist)%varid, &
            easyaerosol_dist(i_easyaerosol_dist)%ndims, &
            easyaerosol_dist(i_easyaerosol_dist)%n_specbands
        CALL umPrint(umMessage,src='easyaerosol_init_input')

        WRITE (umMessage,'(A,2(1X,I0))') &
            'Update frequency and type', &
            easyaerosol_dist(i_easyaerosol_dist)%update_freq, &
            easyaerosol_dist(i_easyaerosol_dist)%update_type
        CALL umPrint(umMessage,src='easyaerosol_init_input')

      END IF
      
      i_easyaerosol_dist = i_easyaerosol_dist + 1
      
    END IF ! if file name is set
  
  END IF ! if file is required
  
END DO ! n

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE easyaerosol_init_input

! ##############################################################################

SUBROUTINE easyaerosol_struct_init
!---------------------------------------------------------------------------
! Description:
!   Initialisation of EasyAerosol structures to default values.
! Method:
!   Initialise elements of EasyAerosols structures.
!   Note that as indicated below some elements cannot be
!   allocated yet within this subroutine.
!---------------------------------------------------------------------------

USE ereport_mod,        ONLY: ereport
IMPLICIT NONE

INTEGER                        :: ierror         ! Error code
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EASYAEROSOL_STRUCT_INIT'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)


ierror = 0

IF ( .NOT. ALLOCATED(easyaerosol_dist)) THEN
  ierror = 1
  CALL ereport('EASYAEROSOL_STRUCT_INIT', ierror,                             &
    ' EasyAerosol super array should have been allocated before this routine')
END IF

easyaerosol_dist(:)%file_name   = ' '
easyaerosol_dist(:)%var_name    = ' '
easyaerosol_dist(:)%varid       = -1
easyaerosol_dist(:)%std_name    = ' '
easyaerosol_dist(:)%long_name   = ' '
easyaerosol_dist(:)%units       = ' '

easyaerosol_dist(:)%update_freq = 0
easyaerosol_dist(:)%update_type = 0
easyaerosol_dist(:)%l_update    = .FALSE.
easyaerosol_dist(:)%ndims       = 3
easyaerosol_dist(:)%n_specbands = 0

! Note that the pointers easyaerosol_*(:)%values_3d or values_4d
! cannot be initialised yet. Before that they need to be allocated
! independently for each element of easyaerosol_*(:). This is done in
! easyaerosol_init_input after calling this routine.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE easyaerosol_struct_init

! #############################################################################

SUBROUTINE easyaerosol_update(row_length, rows, model_levels,                  &
  global_row_length, global_rows, timestep_number,                             &
  timestep, l_first)
! ---------------------------------------------------------------------
! Description:
!   Check if it is time to update the values in the EasyAerosol structure,
!   and do it when needed.
! Method:
! 1) Check first if full hour. If so go through all EasyAerosol fields and
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
USE emiss_io_mod,    ONLY: em_fclose, em_fopen, em_var_check_dims, em_get_data,&
     iupd_single
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

REAL                :: fhournow      ! Target time (fractional hours since
! 0:00 h on 1st Jan)

INTEGER, PARAMETER  :: max_dims = 5            ! Max no of dims per var
INTEGER             :: dimsize(max_dims)       ! array of sizes
REAL, ALLOCATABLE   :: tmpdat1 (:,:,:)         ! Temporal data arrays (3D)
REAL, ALLOCATABLE   :: tmpdat2 (:,:,:)
REAL, ALLOCATABLE   :: tmpdat1_4d(:,:,:,:)     ! Temporal data arrays (4D)
REAL, ALLOCATABLE   :: tmpdat2_4d(:,:,:,:)
INTEGER             :: l                       ! Index
INTEGER             :: k                       ! Loop counter
INTEGER             :: ierror                  ! Error code
LOGICAL             :: l_noreco                ! True if record not found (EOF)
CHARACTER (LEN=errormessagelength) :: cmessage                ! error message

REAL (KIND=jprb)    :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EASYAEROSOL_UPDATE'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Reset l_update to false for all fields.
!  Note that this is changed to TRUE below separately for each
!  NetCDF field if it is updated in the current timestep.
easyaerosol_dist(:)%l_update = .FALSE.

! Get ancil reference time as days/secs since 0/0 basis time and use this to
! calculate hours from ancil reference time to now.
CALL time2sec(ancil_reftime(1), ancil_reftime(2), ancil_reftime(3), &
              ancil_reftime(4), ancil_reftime(5), ancil_reftime(6), &
              0,0,ancil_ref_days,ancil_ref_secs,lcal360)

CALL time2sec(i_year,i_month,i_day,i_hour,i_minute,i_second,         &
              ancil_ref_days,ancil_ref_secs,day,sec ,lcal360)
hours_from_ancref = day*24.0 + sec*rhour_per_sec


! Loop over EasyAerosol fields
DO l = 1, n_easyaerosol_dist

  ! Decide if updating, and set target time for time-interpolation
  IF (easyaerosol_dist(l)%update_type == iupd_single) THEN
  ! Single-time easyaerosol_dist read on first timestep and no other
  ! Set target time to now for get_emfile_rec (not actually used)
    IF (l_first) THEN
      easyaerosol_dist(l)%l_update = .TRUE.
      i_year_trg = i_year
      i_month_trg = i_month
      i_day_trg = i_day
      i_hour_trg = i_hour
      i_minute_trg = i_minute
      i_second_trg = i_second
    END IF
    
  ELSE

    ! Period or timeseries easyaerosol_dist.
    ! Compute integer number of ancil steps from ref time to the begining of the
    ! present update period (=anc_steps). Use FLOOR to round down (even if 
    ! hours_from_ancref is negative) to get time which is at or before now.
    ! Force update when anc_steps increases from that stored at last update.
    ! Always update on first timestep.
    ! Set target time for time-interp to be mid-point of update interval:
    ! start of this update interval plus half the update frequency
    anc_steps = FLOOR(hours_from_ancref / easyaerosol_dist(l)%update_freq)

    IF (l_first .OR. anc_steps /= easyaerosol_dist(l)%last_update) THEN
      easyaerosol_dist(l)%l_update = .TRUE.
      easyaerosol_dist(l)%last_update = anc_steps

      fhr = (REAL(anc_steps) + 0.5) * REAL(easyaerosol_dist(l)%update_freq)
      day = INT(fhr/ 24.0)
      sec = INT((fhr - day*24) * 3600.0)
      CALL sec2time (day, sec, ancil_ref_days, ancil_ref_secs, &
        i_year_trg, i_month_trg, i_day_trg, &
        i_hour_trg, i_minute_trg, i_second_trg, &
        fdaynum, lcal360 )
    END IF
  END IF


  IF (easyaerosol_dist(l)%l_update) THEN
  
      ! Debug update times only for first field
      IF (l == 1 .AND. PrintStatus >= PrStatus_Diag) THEN
        WRITE (umMessage, '(A,6(1x,I0))')                                      &
          'EASYAEROSOL_UPDATE_ZONAL: Update first data field for target time', &
           i_year_trg, i_month_trg, i_day_trg,                                 &
           i_hour_trg, i_minute_trg, i_second_trg
        CALL umPrint(umMessage,src='easyaerosol_update')
      END IF

      ! Open file to read an EasyAerosol field
      CALL em_fopen(easyaerosol_dist(l)%file_name, fid)

      ! Pass current model time (i_year, ..., i_second) and other data to
      !  GET_EMFILE_REC in order to get:
      !  * current hour since basis time: fhournow,
      !  * indices of the matching NetCDF records used for
      !     interpolation: ireco (1:2)
      !  * yr, mon, day, hr, min, sec for those two records: creco (1:2,1:6)
      !  * hour since basis time for those records: ftreco (1:2)

      CALL get_emfile_rec(row_length, rows, fid, l,                            &
        easyaerosol_dist(l)%update_type,                                       &
        (/i_year_trg, i_month_trg, i_day_trg,                                  &
          i_hour_trg, i_minute_trg, i_second_trg /),                           &
        l_first, ireco, creco, fhournow, ftreco)

      ! Later we will call EM_GET_DATA to get the values for registers
      !  ireco(1) & ireco(2), at hours ftreco(1) and ftreco(2), so that we can
      !  interpolate them to the current model time in hours (i.e. to hournow).

      ! Check that matching time records are found
      IF ( ANY(ireco(:) == -1) ) THEN  
        WRITE (umMessage,'(A,6(1x,I0),1x,F15.5)') 'At time:', i_year,          &
          i_month, i_day, i_hour, i_minute, i_second, fhournow
        CALL umPrint(umMessage,src='easyaerosol_update')
        CALL em_fclose(fid)
        ierror = l
        cmessage =                                                             &
          'Required time record not in file ' // easyaerosol_dist(l)%file_name
        CALL ereport('EASYAEROSOL_UPDATE', ierror, cmessage)
      END IF

      ! Go ahead with reading values and interpolations

      ! First get dimensions, in particular the vertical levels, so
      !  that we can allocate the temporary data arrays tmpdat1 and tmpdat2.
      IF (easyaerosol_dist(l)%ndims == 3) THEN
        
        ! 3 dimensions + time
        CALL em_var_check_dims(global_row_length, global_rows, model_levels,   &
                               fid, easyaerosol_dist(l)%var_name, dimsize)
      ELSE
      
        ! 4 dimensions + time
        CALL em_var_check_dims(global_row_length, global_rows, model_levels,   &
                               easyaerosol_dist(l)%n_specbands,                &
                               fid, easyaerosol_dist(l)%var_name, dimsize)

      END IF
      
      l_noreco = .FALSE.  ! initially assume that records are found

      IF ( ireco(1) == ireco(2) ) THEN
        ! No interpolation if exact time match. Get only emiss values for
        !  register ireco(1).
        IF (easyaerosol_dist(l)%ndims == 3) THEN
  
          ! 3 dimensions + time
          CALL em_get_data (row_length, rows, fid, ireco(1),    &
            easyaerosol_dist(l)%var_name,                       &
            easyaerosol_dist(l)%values_3d(:, :, 1:dimsize(3)),  &
            l_noreco)

        ELSE

          ! 4 dimensions + time
          CALL em_get_data (row_length, rows, fid, ireco(1),                 &
            easyaerosol_dist(l)%var_name,                                    &
            easyaerosol_dist(l)%values_4d(:, :, 1:dimsize(3), 1:dimsize(4)), &
            l_noreco)

        END IF

        IF (l_noreco) THEN  
          CALL em_fclose(fid)
          ierror = l
          cmessage =  'Error reading data1: ' // easyaerosol_dist(l)%var_name
          CALL ereport('EASYAEROSOL_UPDATE', ierror, cmessage)
        END IF

      ELSE
        ! Get two values tmpdat1 and tmpdat2, for registers
        !  ireco(1) and ireco(2), to later interpolate.

        IF (easyaerosol_dist(l)%ndims == 3) THEN

          ! 3 dimensions + time

          ALLOCATE(tmpdat1(row_length, rows, dimsize(3)) )
          ALLOCATE(tmpdat2(row_length, rows, dimsize(3)) )

          ! Read in data for previous matching time. Stop if record not found.
          CALL em_get_data(row_length, rows, fid, ireco(1),                    &
                           easyaerosol_dist(l)%var_name,                       &
                           tmpdat1 (:,:,:), l_noreco)
          IF (l_noreco) THEN  
            CALL em_fclose(fid)
            ierror = l
            cmessage = 'Error reading data1: ' // easyaerosol_dist(l)%var_name
            CALL ereport('EASYAEROSOL_UPDATE', ierror, cmessage)
          END IF

          ! Read in data for next matching time.
          CALL em_get_data (row_length, rows, fid, ireco(2),                   &
                            easyaerosol_dist(l)%var_name,                      &
                            tmpdat2 (:,:,:), l_noreco)
          IF (l_noreco) THEN  
            CALL em_fclose(fid)
            ierror = l
            cmessage = 'Error reading data2: ' // easyaerosol_dist(l)%var_name
            CALL ereport('EASYAEROSOL_UPDATE', ierror, cmessage)
          END IF

          ! Interpolate using same UM routine as REPLANCA
          CALL t_int (                                                         &
            tmpdat1 (:,:,:), ftreco(1),               &  ! Data1 and time1
            tmpdat2 (:,:,:), ftreco(2),               &  ! Data2 and time2
            easyaerosol_dist(l)%values_3d(:, :, 1:dimsize(3)), &  ! Data interp
            fhournow,                                 &  ! to current time
            row_length*rows*dimsize(3))                  ! No of points

          DEALLOCATE(tmpdat2)
          DEALLOCATE(tmpdat1)

        ELSE
  
          ! 4 dimensions + time

          ALLOCATE(tmpdat1_4d(row_length, rows, dimsize(3), dimsize(4)) )
          ALLOCATE(tmpdat2_4d(row_length, rows, dimsize(3), dimsize(4)) )

          ! Read in data for previous matching time. Stop if record not found.
          CALL em_get_data(row_length, rows, fid, ireco(1),                    &
                           easyaerosol_dist(l)%var_name,                       &
                           tmpdat1_4d(:,:,:,:), l_noreco)
          IF (l_noreco) THEN  
            CALL em_fclose(fid)
            ierror = l
            cmessage = 'Error reading data1: ' // easyaerosol_dist(l)%var_name
            CALL ereport('EASYAEROSOL_UPDATE', ierror, cmessage)
          END IF

          ! Read in data for next matching time.
          CALL em_get_data (row_length, rows, fid, ireco(2),                   &
                            easyaerosol_dist(l)%var_name,                      &
                            tmpdat2_4d(:,:,:,:), l_noreco)
          IF (l_noreco) THEN  
            CALL em_fclose(fid)
            ierror = l
            cmessage = 'Error reading data2: ' // easyaerosol_dist(l)%var_name
            CALL ereport('EASYAEROSOL_UPDATE', ierror, cmessage)
          END IF
          
          DO k = 1, dimsize(4)
  
            ! Interpolate using same UM routine as REPLANCA
            CALL t_int (                                                       &
              tmpdat1_4d(:,:,:,k), ftreco(1),               & ! Data1 and time1
              tmpdat2_4d(:,:,:,k), ftreco(2),               & ! Data2 and time2
              easyaerosol_dist(l)%values_4d(:, :, &
                                         1:dimsize(3), k), &  ! Data interp
              fhournow,                                 &  ! to current time
              row_length*rows*dimsize(3))                  ! No of points
          
          END DO ! k
  
          DEALLOCATE(tmpdat2_4d)
          DEALLOCATE(tmpdat1_4d)  

        END IF

      END IF      ! Interpolation required?

      ! Close emission file after reading an emission field
      CALL em_fclose(fid)

      ! Write diagnostic results
      IF (PrintStatus >= PrStatus_Oper) THEN
        WRITE(umMessage,'(A,6(1x,I0),F15.5)') &
          'Updated EasyAerosol at time: ',   &
          i_year, i_month, i_day, i_hour, i_minute, i_second, fhournow
        CALL umPrint(umMessage,src='easyaerosol_update')
        
        IF (easyaerosol_dist(l)%ndims == 3) THEN
  
          ! 3 dimensions
  
          WRITE(umMessage,'(A,3E18.6)') &
            'EasyAerosol level 1: (max, min, mean) ',       &
            MAXVAL(easyaerosol_dist(l)%values_3d(:,:,1)),   &
            MINVAL(easyaerosol_dist(l)%values_3d(:,:,1)),   &
            SUM(easyaerosol_dist(l)%values_3d(:,:,1))/      &
            SIZE(easyaerosol_dist(l)%values_3d(:,:,1))
          CALL umPrint(umMessage,src='easyaerosol_update')
        
        ELSE
  
          ! 4 dimensions
          DO k = 1, dimsize(4)
  
            WRITE(umMessage,'(A,I0,1X,3E18.6)') &
              'EasyAerosol level 1: (max, min, mean) Band: ',  &
              k,                                               &
              MAXVAL(easyaerosol_dist(l)%values_4d(:,:,1,k)),  &
              MINVAL(easyaerosol_dist(l)%values_4d(:,:,1,k)),  &
              SUM(easyaerosol_dist(l)%values_4d(:,:,1,k))/     &
             SIZE(easyaerosol_dist(l)%values_4d(:,:,1,k))
          CALL umPrint(umMessage,src='easyaerosol_update')

          END DO ! k  
  
        END IF

      END IF

  END IF        ! Update time?
END DO          ! loop over aerosol dist elements



IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE easyaerosol_update

! #############################################################################

SUBROUTINE easyaerosol_fill(row_length,rows,model_levels, &
  h_swbands, h_lwbands, easyaerosol_sw, easyaerosol_lw, easyaerosol_cdnc)
! ---------------------------------------------------------------------
! Description:
!  To allocate the EasyAerosol arrays, and fill with data from the values
!  read into the easyaerosol_dist structure from the netCDF input fields.
! ---------------------------------------------------------------------

USE def_easyaerosol,    ONLY: t_easyaerosol_rad, t_easyaerosol_cdnc, &
                              allocate_easyaerosol_rad, &
                              allocate_easyaerosol_cdnc

IMPLICIT NONE

! Input arguments with info on model dimensions
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels
INTEGER, INTENT(IN) :: h_swbands
INTEGER, INTENT(IN) :: h_lwbands

! EasyAerosol distributions
TYPE (t_easyaerosol_rad), INTENT(INOUT) :: easyaerosol_sw
TYPE (t_easyaerosol_rad), INTENT(INOUT) :: easyaerosol_lw
TYPE (t_easyaerosol_cdnc), INTENT(INOUT) :: easyaerosol_cdnc

INTEGER :: i_ext, i_abs, i_asy, i_cdn ! array indices

REAL (KIND=jprb)    :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EASYAEROSOL_FILL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF (l_easyaerosol_sw) THEN
  
  CALL allocate_easyaerosol_rad(easyaerosol_sw, row_length, rows, &
                                model_levels, h_swbands)
  
  ! SW distributions are the first 3 distributions of easyaerosol_dist in
  ! all cases.
  i_ext = 1
  i_abs = 2
  i_asy = 3
  
  easyaerosol_sw%extinction(:,:,:,:)= easyaerosol_dist(i_ext)%values_4d(:,:,:,:)
  easyaerosol_sw%absorption(:,:,:,:)= easyaerosol_dist(i_abs)%values_4d(:,:,:,:)
  easyaerosol_sw%asymmetry(:,:,:,:) = easyaerosol_dist(i_asy)%values_4d(:,:,:,:)

ELSE

  ! Allocate smallest space for unused distributions
  CALL allocate_easyaerosol_rad(easyaerosol_sw, 1, 1, 1, 1)
  
END IF

IF (l_easyaerosol_lw) THEN
  
  CALL allocate_easyaerosol_rad(easyaerosol_lw, row_length, rows, &
                                model_levels, h_lwbands)
  
  ! LW distributions are the first 3 distributions of easyaerosol_dist if
  ! there are no SW distributions, and distributions 4 to 6 otherwise.
  IF (.NOT. l_easyaerosol_sw) THEN
  
    i_ext = 1
    i_abs = 2
    i_asy = 3
  
  ELSE
  
    i_ext = 4
    i_abs = 5
    i_asy = 6
  
  END IF
  
  easyaerosol_lw%extinction(:,:,:,:)= easyaerosol_dist(i_ext)%values_4d(:,:,:,:)
  easyaerosol_lw%absorption(:,:,:,:)= easyaerosol_dist(i_abs)%values_4d(:,:,:,:)
  easyaerosol_lw%asymmetry(:,:,:,:) = easyaerosol_dist(i_asy)%values_4d(:,:,:,:)
  
ELSE

  ! Allocate smallest space for unused distributions
  CALL allocate_easyaerosol_rad(easyaerosol_lw, 1, 1, 1, 1)
  
END IF

IF (l_easyaerosol_cdnc .OR. l_easyaerosol_autoconv) THEN

  CALL allocate_easyaerosol_cdnc(easyaerosol_cdnc, row_length, rows, &
                                 model_levels)

  ! CDNC distribution is the first distribution of easyaerosol_dist, unless
  ! -- SW or LW distributions take the first 3 slots
  ! -- SW and LW distributions take the first 6 slots.
  i_cdn = 1
  
  IF (l_easyaerosol_sw) THEN
  
    i_cdn = i_cdn + 3
  
  END IF
  
  IF (l_easyaerosol_lw) THEN
  
    i_cdn = i_cdn + 3
  
  END IF
  
  easyaerosol_cdnc%cdnc(:,:,:) = easyaerosol_dist(i_cdn)%values_3d(:,:,:)
  
  ! Apply a minimum value of 5 cm-3 to cloud droplet number, following
  ! what is done in number_droplet().
  WHERE (easyaerosol_cdnc%cdnc < 5.0e06) easyaerosol_cdnc%cdnc = 5.0e06
  
ELSE

  ! Allocate smallest space for unused distributions
  CALL allocate_easyaerosol_cdnc(easyaerosol_cdnc, 1, 1, 1)

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE easyaerosol_fill

! #############################################################################

SUBROUTINE easyaerosol_update_zonal(row_length, rows, model_levels,            &
                     global_rows, timestep_number,                             &
                     timestep, l_first)
! ---------------------------------------------------------------------
! Description:
!   Check if it is time to update the values in the EasyAerosol structure,
!   and do it when needed.
! Method:
! 1) Check first if full hour. If so go through all EasyAerosol fields and
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
USE emiss_io_mod,    ONLY: em_fclose, em_fopen, em_var_check_dims, em_get_data,&
     iupd_single
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

REAL                :: fhournow      ! Target time (fractional hours since
! 0:00 h on 1st Jan)

INTEGER, PARAMETER  :: max_dims = 5            ! Max no of dims per var
INTEGER             :: dimsize(max_dims)       ! array of sizes
REAL, ALLOCATABLE   :: tmpdat1 (:,:,:)         ! Temporal data arrays (3D)
REAL, ALLOCATABLE   :: tmpdat2 (:,:,:)
REAL, ALLOCATABLE   :: tmpdat1_4d(:,:,:,:)     ! Temporal data arrays (4D)
REAL, ALLOCATABLE   :: tmpdat2_4d(:,:,:,:)
INTEGER             :: l                       ! Index
INTEGER             :: k                       ! Loop counter
INTEGER             :: ierror                  ! Error code
LOGICAL             :: l_noreco                ! True if record not found (EOF)
CHARACTER (LEN=errormessagelength) :: cmessage                ! error message

REAL (KIND=jprb)    :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EASYAEROSOL_UPDATE_ZONAL'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Reset l_update to false for all fields.
!  Note that this is changed to TRUE below separately for each
!  NetCDF field if it is updated in the current timestep.
easyaerosol_dist(:)%l_update = .FALSE.

! Get ancil reference time as days/secs since 0/0 basis time and use this to
! calculate hours from ancil reference time to now.
CALL time2sec(ancil_reftime(1), ancil_reftime(2), ancil_reftime(3), &
              ancil_reftime(4), ancil_reftime(5), ancil_reftime(6), &
              0,0,ancil_ref_days,ancil_ref_secs,lcal360)

CALL time2sec(i_year,i_month,i_day,i_hour,i_minute,i_second,         &
              ancil_ref_days,ancil_ref_secs,day,sec ,lcal360)
hours_from_ancref = day*24.0 + sec*rhour_per_sec


! Loop over EasyAerosol fields
DO l = 1, n_easyaerosol_dist

  ! Decide if updating, and set target time for time-interpolation
  IF (easyaerosol_dist(l)%update_type == iupd_single) THEN
  ! Single-time easyaerosol_dist read on first timestep and no other
  ! Set target time to now for get_emfile_rec (not actually used)
    IF (l_first) THEN
      easyaerosol_dist(l)%l_update = .TRUE.
      i_year_trg = i_year
      i_month_trg = i_month
      i_day_trg = i_day
      i_hour_trg = i_hour
      i_minute_trg = i_minute
      i_second_trg = i_second
    END IF

  ELSE

    ! Period or timeseries easyaerosol_dist.
    ! Compute integer number of ancil steps from ref time to the begining of the
    ! present update period (=anc_steps). Use FLOOR to round down (even if 
    ! hours_from_ancref is negative) to get time which is at or before now.
    ! Force update when anc_steps increases from that stored at last update.
    ! Always update on first timestep.
    ! Set target time for time-interp to be mid-point of update interval:
    ! start of this update interval plus half the update frequency
    anc_steps = FLOOR(hours_from_ancref / easyaerosol_dist(l)%update_freq)

    IF (l_first .OR. anc_steps /= easyaerosol_dist(l)%last_update) THEN
      easyaerosol_dist(l)%l_update = .TRUE.
      easyaerosol_dist(l)%last_update = anc_steps

      fhr = (REAL(anc_steps) + 0.5) * REAL(easyaerosol_dist(l)%update_freq)
      day = INT(fhr/ 24.0)
      sec = INT((fhr - day*24) * 3600.0)
      CALL sec2time (day, sec, ancil_ref_days, ancil_ref_secs, &
        i_year_trg, i_month_trg, i_day_trg, &
        i_hour_trg, i_minute_trg, i_second_trg, &
        fdaynum, lcal360 )
    END IF
  END IF


  IF (easyaerosol_dist(l)%l_update) THEN
  
    ! Debug update times only for first field
    IF (l == 1 .AND. PrintStatus >= PrStatus_Diag) THEN
      WRITE (umMessage, '(A,6(1x,I0))')                                      &
        'EASYAEROSOL_UPDATE_ZONAL: Update first data field for target time', &
         i_year_trg, i_month_trg, i_day_trg,                                 &
         i_hour_trg, i_minute_trg, i_second_trg
      CALL umPrint(umMessage,src='easyaerosol_update_zonal')
    END IF

    ! Open file to read an EasyAerosol field
    CALL em_fopen(easyaerosol_dist(l)%file_name, fid)

    ! Pass current model time (i_year, ..., i_second) and other data to
    !  GET_EMFILE_REC in order to get:
    !  * current hour since basis time: fhournow,
    !  * indices of the matching NetCDF records used for
    !     interpolation: ireco (1:2)
    !  * yr, mon, day, hr, min, sec for those two records: creco (1:2,1:6)
    !  * hour since basis time for those records: ftreco (1:2)

    CALL get_emfile_rec(1, rows, fid, l,                                     &
      easyaerosol_dist(l)%update_type,                                       &
      (/i_year_trg, i_month_trg, i_day_trg,                                  &
        i_hour_trg, i_minute_trg, i_second_trg /),                           &
      l_first, ireco, creco, fhournow, ftreco)

    ! Later we will call EM_GET_DATA to get the values for registers
    !  ireco(1) & ireco(2), at hours ftreco(1) and ftreco(2), so that we can
    !  interpolate them to the current model time in hours (i.e. to hournow).

    ! Check that matching time records are found
    IF ( ANY(ireco(:) == -1) ) THEN  
      WRITE (umMessage,'(A,6(1x,I0),1x,F15.5)') 'At time:', i_year,          &
        i_month, i_day, i_hour, i_minute, i_second, fhournow
      CALL umPrint(umMessage,src='easyaerosol_update_zonal')
      CALL em_fclose(fid)
      ierror = l
      cmessage =                                                             &
        'Required time record not in file ' // easyaerosol_dist(l)%file_name
      CALL ereport('EASYAEROSOL_UPDATE_ZONAL', ierror, cmessage)
    END IF

    ! Go ahead with reading values and interpolations

    ! First get dimensions, in particular the vertical levels, so
    !  that we can allocate the temporary data arrays tmpdat1 and tmpdat2.
    IF (easyaerosol_dist(l)%ndims == 3) THEN
        
      ! File should have 3 dimensions (lat, lev, time)
      CALL em_var_check_dims(global_rows, model_levels,    &
                             l_easyaerosol_zonal,          &
                             fid, easyaerosol_dist(l)%var_name, dimsize)
    ELSE
      
      ! File should have 4 dimensions (lat, lev, waveband, time)
      CALL em_var_check_dims(global_rows, model_levels,       &
                             easyaerosol_dist(l)%n_specbands, &
                             l_easyaerosol_zonal,             &
                             fid, easyaerosol_dist(l)%var_name, dimsize)
    END IF
      
    l_noreco = .FALSE.  ! initially assume that records are found

    IF ( ireco(1) == ireco(2) ) THEN
      ! No interpolation if exact time match. Get only emiss values for
      !  register ireco(1).
      IF (easyaerosol_dist(l)%ndims == 3) THEN
  
        ! File has 2 spatial dimensions (lat, lev) but these will be expanded
        ! in em_get_data_real2d_zonal to 3 spatial dimensions
        ! so data fields in values_3d are (lon, lat, lev)
        !
        CALL em_get_data(row_length, rows, fid, ireco(1),               &
          l_easyaerosol_zonal, easyaerosol_dist(l)%var_name,            &
          easyaerosol_dist(l)%values_3d(:,:, 1:dimsize(2)),             &
          l_noreco)

      ELSE

        ! File has 3 dimensions (lat, lev, spectral bands)
        ! but the data will be expanded in em_get_data_3d_zonal 
        ! to 3 spatial dimensions so data fields in values_4d
        ! are (lon, lat, lev, spectral bands)
        !    
        CALL em_get_data(row_length, rows, fid, ireco(1),                 &
          l_easyaerosol_zonal, easyaerosol_dist(l)%var_name,              &
          easyaerosol_dist(l)%values_4d(:,:, 1:dimsize(2), 1:dimsize(3)), &
          l_noreco)

      END IF

      IF (l_noreco) THEN  
        CALL em_fclose(fid)
        ierror = l
        cmessage =  'Error reading data1: ' // easyaerosol_dist(l)%var_name
        CALL ereport('EASYAEROSOL_UPDATE_ZONAL', ierror, cmessage)
      END IF

    ELSE
      ! Get two values tmpdat1 and tmpdat2, for registers
      !  ireco(1) and ireco(2), to later interpolate.

      IF (easyaerosol_dist(l)%ndims == 3) THEN

        ! Data field has 3 dimensions (lon, lat, lev) + time

        ALLOCATE(tmpdat1(row_length, rows, dimsize(2)) )
        ALLOCATE(tmpdat2(row_length, rows, dimsize(2)) )

        ! Read in data for previous matching time. Stop if record not found.
        CALL em_get_data(row_length, rows, fid, ireco(1),                 &
          l_easyaerosol_zonal, easyaerosol_dist(l)%var_name,              &
          tmpdat1 (:,:,:), l_noreco)

        IF (l_noreco) THEN  
          CALL em_fclose(fid)
          ierror = l
          cmessage = 'Error reading data1: ' // easyaerosol_dist(l)%var_name
          CALL ereport('EASYAEROSOL_UPDATE_ZONAL', ierror, cmessage)
        END IF

        ! Read in data for next matching time.
        CALL em_get_data (row_length, rows, fid, ireco(2),                   &
          l_easyaerosol_zonal, easyaerosol_dist(l)%var_name,                 &
          tmpdat2 (:,:,:), l_noreco)

        IF (l_noreco) THEN  
          CALL em_fclose(fid)
          ierror = l
          cmessage = 'Error reading data2: ' // easyaerosol_dist(l)%var_name
          CALL ereport('EASYAEROSOL_UPDATE_ZONAL', ierror, cmessage)
        END IF

        ! Interpolate using same UM routine as REPLANCA
        CALL t_int (                                                         &
          tmpdat1 (:,:,:), ftreco(1),               &  ! Data1 and time1
          tmpdat2 (:,:,:), ftreco(2),               &  ! Data2 and time2
          easyaerosol_dist(l)%values_3d(:, :, 1:dimsize(2)), &  ! Data interp
          fhournow,                                 &  ! to current time
          row_length*rows*dimsize(2))                  ! No of points

        DEALLOCATE(tmpdat2)
        DEALLOCATE(tmpdat1)

      ELSE
  
        ! Data fields have 4 dimensions (lon, lat, lev, waveband) + time

        ALLOCATE(tmpdat1_4d(row_length, rows, dimsize(2), dimsize(3)) )
        ALLOCATE(tmpdat2_4d(row_length, rows, dimsize(2), dimsize(3)) )

        ! Read in data for previous matching time. Stop if record not found.
        CALL em_get_data(row_length, rows, fid, ireco(1),                  &
          l_easyaerosol_zonal, easyaerosol_dist(l)%var_name,               &
          tmpdat1_4d(:,:,:,:), l_noreco)

        IF (l_noreco) THEN  
          CALL em_fclose(fid)
          ierror = l
          cmessage = 'Error reading data1: ' // easyaerosol_dist(l)%var_name
          CALL ereport('EASYAEROSOL_UPDATE_ZONAL', ierror, cmessage)
        END IF

        ! Read in data for next matching time.
        CALL em_get_data(row_length, rows, fid, ireco(2),                   &
          l_easyaerosol_zonal, easyaerosol_dist(l)%var_name,                &
          tmpdat2_4d(:,:,:,:), l_noreco)

        IF (l_noreco) THEN  
          CALL em_fclose(fid)
          ierror = l
          cmessage = 'Error reading data2: ' // easyaerosol_dist(l)%var_name
          CALL ereport('EASYAEROSOL_UPDATE_ZONAL', ierror, cmessage)
        END IF
          
        DO k = 1, dimsize(3)
  
          ! Interpolate using same UM routine as REPLANCA
          CALL t_int (                                                       &
            tmpdat1_4d(:,:,:,k), ftreco(1),               & ! Data1 and time1
            tmpdat2_4d(:,:,:,k), ftreco(2),               & ! Data2 and time2
            easyaerosol_dist(l)%values_4d(:, :, &
                                        1:dimsize(2), k), &  ! Data interp
            fhournow,                                 &  ! to current time
            row_length*rows*dimsize(2))                  ! No of points
          
        END DO ! k
  
        DEALLOCATE(tmpdat2_4d)
        DEALLOCATE(tmpdat1_4d)  

      END IF

    END IF      ! Interpolation required?

    ! Close emission file after reading an emission field
    CALL em_fclose(fid)

    ! Write diagnostic results
    IF (PrintStatus >= PrStatus_Oper) THEN
      WRITE(umMessage,'(A,6(1x,I0),F15.5)') &
        'Updated EasyAerosol at time: ',   &
        i_year, i_month, i_day, i_hour, i_minute, i_second, fhournow
      CALL umPrint(umMessage,src='easyaerosol_update_zonal')
        
      IF (easyaerosol_dist(l)%ndims == 3) THEN
  
        ! Data field has 3 dimensions 
  
        WRITE(umMessage,'(A,3E18.6)') &
          'EasyAerosol level 1: (max, min, mean) ',       &
          MAXVAL(easyaerosol_dist(l)%values_3d(:,:,1)),   &
          MINVAL(easyaerosol_dist(l)%values_3d(:,:,1)),   &
          SUM(easyaerosol_dist(l)%values_3d(:,:,1))/      &
          SIZE(easyaerosol_dist(l)%values_3d(:,:,1))
        CALL umPrint(umMessage,src='easyaerosol_update_zonal')
        
      ELSE
  
        ! Data field has 4 dimensions
        DO k = 1, dimsize(4)
  
          WRITE(umMessage,'(A,I0,1X,3E18.6)') &
            'EasyAerosol level 1: (max, min, mean) Band: ',  &
            k,                                               &
            MAXVAL(easyaerosol_dist(l)%values_4d(:,:,1,k)),  &
            MINVAL(easyaerosol_dist(l)%values_4d(:,:,1,k)),  &
            SUM(easyaerosol_dist(l)%values_4d(:,:,1,k))/     &
            SIZE(easyaerosol_dist(l)%values_4d(:,:,1,k))
          CALL umPrint(umMessage,src='easyaerosol_update')

        END DO ! k  

      END IF

    END IF

  END IF        ! Update time?
END DO          ! loop over aerosol dist elements



IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE easyaerosol_update_zonal

END MODULE easyaerosol_read_input_mod
