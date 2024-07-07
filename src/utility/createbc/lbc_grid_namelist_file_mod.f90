! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE lbc_grid_namelist_file_mod

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.
!
! Contains the namelist lbc_grid which controls the operation of CreateBC,
! and the procedure which reads the namelist and loads the contents into the
! lbc_output_control_type from where it is used by the rest of the code.

USE lbc_output_control_mod, ONLY: lbc_output_control_type,        &
                                  max_stash_items,                &
                                  max_input_files, create_frame,  &
                                  packing_unchanged

USE Packing_Codes_Mod, ONLY: PC_WGDOS_Packing

USE file_mod,           ONLY: file_type
USE filenamelength_mod, ONLY: filenamelength
USE missing_data_mod,   ONLY: imdi, rmdi
! DrHook Modules
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

REAL    :: dlat                         = rmdi ! Latitudinal spacing
REAL    :: dlong                        = rmdi ! Longitudinal spacing
INTEGER :: end_time(6)                  = imdi ! End time for LBC generation
INTEGER :: horizontal_interpolation_method = imdi ! Horiz interp method
INTEGER :: nlat                         = imdi ! Number of latitudinal rows
INTEGER :: nlong                        = imdi ! Number of longitudinal columns
INTEGER :: num_levels                   = imdi ! Number of vertical levels
INTEGER :: output_grid_stagger          = imdi ! Output grid stagger:
                                               !         3 = ND, 6 = EG
REAL    :: pole_lat                     = rmdi ! Latitude of rotated pole
REAL    :: pole_long                    = rmdi ! Longitude of rotated pole
REAL    :: start_lat                    = rmdi ! Latitude of first p-point
REAL    :: start_long                   = rmdi ! Longitude of first p-point
INTEGER :: start_time(6)                = imdi ! Start time for LBC generation
INTEGER :: stash_codes(max_stash_items) = imdi ! Array of stash items 
LOGICAL :: variable_resolution = .FALSE.       !  var-res output flag
INTEGER :: vertical_interpolation_method = imdi ! Vertical interp method
INTEGER :: halo_lat                     = imdi  ! Halo size in y direction
INTEGER :: halo_long                    = imdi ! Halo size in x direction
INTEGER :: rim_width                    = imdi ! Rim width
REAL    :: q_min                        = rmdi ! Reset specific humidity to 
                                               !      minimum if value below
                                               !      minimum
INTEGER :: num_reserved_headers         = imdi ! Number of lookup headers in 
                                               !      output file
INTEGER :: frames_cutout_adjust_north   = imdi
INTEGER :: frames_cutout_adjust_south   = imdi
INTEGER :: frames_cutout_adjust_east    = imdi
INTEGER :: frames_cutout_adjust_west    = imdi
INTEGER :: frames_packing_option        = imdi

INTEGER :: time_interval                = imdi ! Time interval between LBCs
                                               ! in seconds

LOGICAL :: write_header_only_once       = .FALSE. ! Write the header once (true) 
                                                  ! or after every field (false)

! Filenames of input files
CHARACTER(LEN=filenamelength) :: input_data(max_input_files) = 'unset'

! Filename for variable res horizontal namelists
CHARACTER(LEN=filenamelength) :: horizontal_grid_file = 'unset'  

! Output LBC/Frame filename
CHARACTER(LEN=filenamelength) :: output_filename = 'unset'

! Vertical levels namelist file
CHARACTER(LEN=filenamelength) :: vertical_levels_file = 'unset'

! Path of directory containing STASHmaster_A file
CHARACTER(LEN=filenamelength) :: stashmaster_path = 'unset'

NAMELIST /lbc_grid/ dlat,                                 &
                    dlong,                                &
                    nlat,                                 &
                    nlong,                                &
                    num_levels,                           &
                    output_grid_stagger,                  &
                    pole_lat,                             &
                    pole_long,                            &
                    start_lat,                            & 
                    start_long,                           &
                    stash_codes,                          &
                    variable_resolution,                  &
                    vertical_interpolation_method,        &
                    end_time,                             &
                    start_time,                           &
                    input_data,                           &
                    output_filename,                      &
                    vertical_levels_file,                 &
                    horizontal_grid_file,                 & 
                    horizontal_interpolation_method,      &
                    halo_lat,                             &
                    halo_long,                            &
                    rim_width,                            &
                    q_min,                                &
                    stashmaster_path,                     &
                    num_reserved_headers,                 &
                    frames_cutout_adjust_north,           &
                    frames_cutout_adjust_south,           &
                    frames_cutout_adjust_east,            &
                    frames_cutout_adjust_west,            &
                    time_interval,                        &
                    frames_packing_option,                &
                    write_header_only_once

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                        &
  ModuleName = 'LBC_GRID_NAMELIST_FILE_MOD'

TYPE, EXTENDS(file_type) :: lbc_grid_namelist_file_type

  CONTAINS
  PROCEDURE, PASS :: read_namelist

END TYPE lbc_grid_namelist_file_type


CONTAINS


FUNCTION read_namelist(this) RESULT(namelist_in)
! Read the namelist and populate an lbc_output_control_type object with the 
! contents

USE um_stashcode_mod, ONLY: stashcode_dust1_mmr, stashcode_dust2_mmr,   &
      stashcode_dust3_mmr, stashcode_dust4_mmr, stashcode_dust5_mmr, stashcode_dust6_mmr
USE lbc_stashcode_mapping_mod, ONLY: lbc_stashcode_mapping

USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE check_iostat_mod, ONLY: check_iostat
USE fieldsfile_constants_mod, ONLY: new_dynamics, endgame

IMPLICIT NONE

TYPE(lbc_output_control_type) :: namelist_in
CLASS(lbc_grid_namelist_file_type), INTENT(INOUT) :: this
INTEGER :: i, l, m, stash
INTEGER :: standard_stash_codes(max_stash_items)
INTEGER :: dust_stash_codes(max_stash_items)

INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER :: routinename = 'READ_NAMELIST'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Read the namelist from file
READ(this%unit_num, NML=lbc_grid, IOSTAT=icode, IOMSG=cmessage)

! Check the namelist read successfully
CALL check_iostat(icode, "namelist lbc_grid", cmessage)

! do defensive checking of inputs
! If we call ereport from this function, it should always abort with a 
! fatal error
icode = 1

! Vertical grid checks
IF (num_levels < 1) THEN
  cmessage = 'Output grid must have a positive number of levels'
  CALL ereport(routinename, icode, cmessage)
END IF
IF (TRIM(vertical_levels_file) == 'unset') THEN
  cmessage = 'Vertical levels file must be set in namelist'
  CALL ereport(routinename, icode, cmessage)
END IF

! Horizontal grid checks
IF (.NOT. (output_grid_stagger == new_dynamics .OR.                         &
           output_grid_stagger == endgame)) THEN
  cmessage = 'Output grid staggering not implemented'
  CALL ereport(routinename, icode, cmessage)
END IF

IF (pole_lat == rmdi .OR. pole_long == rmdi) THEN
  cmessage = 'Pole co-ordinates not specified'
  CALL ereport(routinename, icode, cmessage)
END IF

IF (variable_resolution) THEN
  IF (TRIM(horizontal_grid_file) == 'unset') THEN
    cmessage = 'Variable resolution flag set but no horizontal grid file'//  &
               ' specified'
    CALL ereport(routinename, icode, cmessage)
  END IF
ELSE
  IF (start_lat == rmdi .OR. start_long == rmdi) THEN
    cmessage = 'Starting latitude and/or longitude not specified '//         &
               '(start_lat and start_long)'
    CALL ereport(routinename, icode, cmessage)
  END IF
  IF (dlat == rmdi .OR. dlong == rmdi) THEN
    cmessage = 'Latitudinal and/or longitudinal spacing not specified '//    &
               '(dlat and dlong)'
    CALL ereport(routinename, icode, cmessage)
  END IF
  IF (nlat == imdi .OR. nlong == imdi) THEN
    cmessage = 'Number of rows and/or columns not specified (nlat and nlong)'
    CALL ereport(routinename, icode, cmessage)
  END IF
END IF

! Rim and haloes
IF (rim_width < 0) THEN
  cmessage = 'Rim width cannot be negative'
  CALL ereport(routinename, icode, cmessage)
END IF
IF (halo_long < 0) THEN
  cmessage = 'Longitudinal halo cannot be negative'
  CALL ereport(routinename, icode, cmessage)
END IF

IF (halo_lat < 0) THEN
  cmessage = 'Latitudinal halo cannot be negative'
  CALL ereport(routinename, icode, cmessage)
END IF

! Interpolation methods
IF (horizontal_interpolation_method < 0 .OR.                                 &
    horizontal_interpolation_method > 2 .OR.                                 &
    horizontal_interpolation_method == 1) THEN
    cmessage = 'Invalid horizontal interpolation method'
    CALL ereport(routinename, icode, cmessage)
END IF
IF (vertical_interpolation_method < 0 .OR.                                   &
    vertical_interpolation_method > 5 .OR.                                   &
    vertical_interpolation_method == 4) THEN
  cmessage = 'Invalid vertical interpolation method'
  CALL ereport(routinename, icode, cmessage)
END IF

! Input/output/STASHmaster files
IF (TRIM(input_data(1)) == 'unset') THEN
  cmessage = 'No input files specified'
  CALL ereport(routinename, icode, cmessage)
END IF
IF (TRIM(output_filename) == 'unset') THEN
  cmessage = 'No output file specified'
  CALL ereport(routinename, icode, cmessage)
END IF
IF (TRIM(stashmaster_path) == 'unset') THEN
  cmessage = 'No STASHmaster path specified'
  CALL ereport(routinename, icode, cmessage)
END IF

! Science settings
IF (stash_codes(1) == imdi) THEN
  cmessage = 'No STASH codes to interpolate specified'
  CALL ereport(routinename, icode, cmessage)
END IF

IF ((frames_packing_option < packing_unchanged                                 &
     .OR. frames_packing_option > PC_WGDOS_Packing)                            &
    .AND. horizontal_interpolation_method == create_frame)  THEN
  cmessage = 'Unrecognised packing option for Frames'
  icode = 99
  CALL ereport(routinename, icode, cmessage)
END IF
    
! Copy the contents of the namelist into an lbc_output_control_type object
namelist_in%dlat = dlat
namelist_in%dlong = dlong
namelist_in%nlat = nlat
namelist_in%nlong = nlong
namelist_in%num_levels = num_levels
namelist_in%output_grid_stagger = output_grid_stagger
namelist_in%pole_lat = pole_lat
namelist_in%pole_long = pole_long
namelist_in%start_lat = start_lat
namelist_in%start_long = start_long
namelist_in%variable_resolution = variable_resolution
namelist_in%vertical_interpolation_method = vertical_interpolation_method
namelist_in%end_time = end_time
namelist_in%start_time = start_time
namelist_in%input_data = input_data
namelist_in%output_filename = output_filename
namelist_in%vertical_levels_file = vertical_levels_file
namelist_in%horizontal_grid_file = horizontal_grid_file
namelist_in%horizontal_interpolation_method = horizontal_interpolation_method
namelist_in%halo_lat = halo_lat
namelist_in%halo_long = halo_long
namelist_in%rim_width = rim_width
namelist_in%q_min = q_min
namelist_in%stashmaster_path = stashmaster_path
namelist_in%num_reserved_headers = num_reserved_headers
namelist_in%frames_cutout_adjust_north = frames_cutout_adjust_north
namelist_in%frames_cutout_adjust_south = frames_cutout_adjust_south
namelist_in%frames_cutout_adjust_east = frames_cutout_adjust_east
namelist_in%frames_cutout_adjust_west = frames_cutout_adjust_west
namelist_in%time_interval = time_interval
namelist_in%frames_packing_option = frames_packing_option
namelist_in%write_header_only_once = write_header_only_once

! STASH-specific handling
l = 0
m = 0
dust_stash_codes = imdi
standard_stash_codes = imdi

! Separate out the dush stash codes, they need special handling to allow
! for dust bin conversion
DO i = 1, max_stash_items
  SELECT CASE (stash_codes(i))
  CASE (stashcode_dust1_mmr, stashcode_dust2_mmr, stashcode_dust3_mmr,       &
       stashcode_dust4_mmr, stashcode_dust5_mmr, stashcode_dust6_mmr)
    l = l + 1
    dust_stash_codes(l) = stash_codes(i)
  CASE DEFAULT
    m = m + 1
    standard_stash_codes(m) = stash_codes(i)
  END SELECT
END DO

namelist_in%stash_codes = standard_stash_codes
namelist_in%dust_stash_codes = dust_stash_codes

! Keep a list of the original order of STASH codes, converting to section 31
! codes in the process. This is used to sort the lookup.
DO i = 1, max_stash_items
  IF (stash_codes(i) == imdi) THEN
    namelist_in%target_stash_codes(i) = imdi
  ELSE 
    CALL lbc_stashcode_mapping(stash_codes(i), stash)
    namelist_in%target_stash_codes(i) = stash
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION read_namelist


END MODULE lbc_grid_namelist_file_mod
