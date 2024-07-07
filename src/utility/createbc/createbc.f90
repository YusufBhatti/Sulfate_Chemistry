! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
PROGRAM createbc

! Classes
USE datafile_mod,                 ONLY: datafile_type
USE fieldsfile_mod,               ONLY: fieldsfile_type
USE field_mod,                    ONLY: field_type, ASSIGNMENT(=)
USE three_dimensional_grid_mod,   ONLY: three_dimensional_grid_type
USE lbcfile_mod,                  ONLY: lbcfile_type
USE wind_rotation_coeff_mod,      ONLY: wind_rotation_coeff_type

! CreateBC Modules
USE lbc_output_control_mod,       ONLY: lbc_output_control_type, create_frame,           &
                                        max_input_files, max_fields_per_file,            &
                                        packing_unchanged
USE lbc_grid_namelist_file_mod,   ONLY: lbc_grid_namelist_file_type
USE stashmaster_utils_mod,        ONLY: load_stashmaster, number_of_dust_fields
USE stashmaster_constants_mod,    ONLY: p_points, u_points, v_points
USE dust_field_conversion_mod,    ONLY: process_dust
USE copy_file_header_mod,         ONLY: copy_file_header
USE interp_control_mod,           ONLY: interp_control
USE interp_weights_mod,           ONLY: weights_index, u_to_enlarged_p, v_to_enlarged_p, &
                                        p_to_enlarged_p, interp_weights_type,            &
                                        num_halo_sizes, num_grid_types
USE check_pole_rotation_mod,      ONLY: check_pole_rotation
USE update_file_header_mod,       ONLY: update_file_header
USE check_vertical_interp_mod,    ONLY: check_vertical_interp
USE process_orography_mod,        ONLY: process_orography, p_to_u_points, p_to_v_points, &
                                        u_points_orog, v_points_orog, p_points_orog,     &
                                        process_orography_enlarged_grid
USE post_interp_transform_mod,    ONLY: post_interp_transform
USE file_manager,                 ONLY: assign_file_unit, release_file_unit
USE autodetect_file_mod,          ONLY: detect_file_type, filetype_fieldsfile,           &
                                        filetype_grib2
USE process_winds_mod,            ONLY: process_winds
USE calc_wind_rotation_coeff_mod, ONLY: calc_wind_rotation_coeff
USE calc_frame_grids_mod,         ONLY: calc_frame_grids
USE update_frame_field_grid_mod,  ONLY: update_frame_field_grid
USE find_fields_to_interp_mod,    ONLY: find_fields_to_interp,               &
                                        found_fields_to_interp,              &
                                        fields_to_interp,                    &
                                        clear_fields_to_interp

! UM modules
USE um_parparams,            ONLY: halo_type_single, halo_type_extended, halo_type_no_halo
USE hostname_mod,            ONLY: get_hostname
USE filenamelength_mod,      ONLY: filenamelength
USE ereport_mod,             ONLY: ereport
USE errormessagelength_mod,  ONLY: errormessagelength
USE um_stashcode_mod,        ONLY: stashcode_orog, stashcode_dust1_mmr, stashcode_dust2_mmr, &
     stashcode_dust3_mmr, stashcode_dust4_mmr, stashcode_dust5_mmr, stashcode_dust6_mmr,     &
     stashcode_u, stashcode_v, stashcode_u_adv, stashcode_v_adv
USE umPrintMgr,              ONLY: umPrint, umMessage
USE UM_Config,               ONLY: appInit, appTerminate
USE application_description, ONLY: exe_createbc
USE planet_constants_mod,    ONLY: i_planet, ip_earth, set_planet_constants

! DrHook Modules
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE drhook_control_mod, ONLY: drhook_control_enable, drhook_control_disable

!$ USE omp_lib, ONLY: omp_get_num_threads, openmp_version

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

! Description:
!   Top level control of createbc
!
! Method:
! 1. Read namelist(s)
! 2. Create output file header
! 3. Loop over files 
! 3a. Read file header
! 3b. Look for fields to interpolate
! 3c. Loop over fields in that file
! 4. Finalise
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.

INTEGER :: j, k, p, halo_code, input_file_num  ! Loop counters
CHARACTER(LEN=filenamelength)           :: namelist_filename
CLASS(fieldsfile_type),    ALLOCATABLE  :: input_file(:)
CLASS(datafile_type),      ALLOCATABLE  :: output_file
TYPE(interp_weights_type), ALLOCATABLE  :: interp_weights_array(:,:)
TYPE(interp_weights_type), ALLOCATABLE  :: rotation_interp_weights_array(:)
TYPE(wind_rotation_coeff_type) :: wind_rotation_coeff
TYPE(lbc_output_control_type)          :: lbc_output_control
TYPE(lbc_grid_namelist_file_type)       :: nlfile

INTEGER :: new_field_num ! Store reference number of field just added
TYPE(field_type) :: output_field

CHARACTER(LEN=19) :: datetime_string
LOGICAL :: copied_header

! Used to check grid is consistent between all input files
TYPE(three_dimensional_grid_type) :: previous_file_p_grid

! Orography 
TYPE(field_type) :: orography_lbc_fields(3,3)
TYPE(field_type) :: orography_enlarged_grid_lbc_field
LOGICAL :: processed_orography

! Vertical interpolation control
LOGICAL :: l_vert_interp = .FALSE.

! Error reporting
INTEGER :: icode, arglen
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER :: routinename = 'CREATEBC'

! Delete existing output file variable
INTEGER :: tmp_unit
LOGICAL :: output_exists = .FALSE.

! Input file format
INTEGER :: input_file_format

! Dust 
LOGICAL :: convert_dust_bins

! Detect if we've written a field to disk
LOGICAL :: fields_written = .FALSE.

! Remember if we've warned about failures to change packing
LOGICAL :: packing_warned = .FALSE.

! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! End of variable declarations, main program starts here

!   Permit calls to DrHook, then call the first top-level caliper.
CALL drhook_control_enable()

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! 0. Need a call to appInit  - Now sets up output PEs
CALL appInit(exe_createbc)

CALL timer( 'createbc', 1)

!$OMP PARALLEL DEFAULT(NONE)
!$OMP MASTER
!$  WRITE(umMessage,'(A,I0,A)') '[INFO] I am running with ',                 &
!$      omp_get_num_threads(),' OpenMP thread(s).'
!$  CALL umPrint(umMessage,src='createbc')
!$  WRITE(umMessage,'(A,I0)') '[INFO] OpenMP Specification: ',openmp_version
!$  CALL umPrint(umMessage,src='createbc')
!$OMP END MASTER
!$OMP END PARALLEL

! Planet constants: set planet as Earth
i_planet = ip_earth
CALL set_planet_constants()

! 1. Read namelist(s)

! Get the filename from command line
CALL get_command_argument(1, namelist_filename, arglen, icode)
SELECT CASE(icode)
CASE (0)
  CONTINUE
CASE (1) 
  cmessage = 'No namelist filename provided on command line'
  CALL ereport(routinename, icode, cmessage)
CASE (-1)
  cmessage = 'The filename and path to the namelist file is too long ' //    &
             'for currently compiled string declaration. Please '     //     &
             'recompile with a larger filenamelength parameter or '   //     &
             'reconsider location of the provided namelist.'
  icode = 1 ! Force fatal error
  CALL ereport(routinename, icode, cmessage)
CASE DEFAULT
  cmessage = 'Unknown error reading namelist file from command line.'
  CALL ereport(routinename, icode, cmessage)
END SELECT    

WRITE(umMessage,'(A,A)') '[INFO] Namelist file is ',namelist_filename
CALL umPrint(umMessage, src='createbc')

! Write out hostname
WRITE(umMessage,'(A,A)') '[INFO] Host is ',TRIM(get_hostname())
CALL umPrint(umMessage, src='createbc')

! Read and process the namelists, generating the LBC grids
CALL nlfile%open_file(namelist_filename)
lbc_output_control = nlfile%read_namelist()
CALL lbc_output_control%process()
CALL nlfile%close_file()

! Load the STASHmaster
CALL load_stashmaster(lbc_output_control%stashmaster_path)

! 2. Create output file header
IF (lbc_output_control%horizontal_interpolation_method == create_frame) THEN
  ! We're creating a Frame, therefore the output file is a fieldsfile
  ALLOCATE(fieldsfile_type :: output_file)
  WRITE(umMessage,'(A)') '[INFO] Output file type is FieldsFile'
  CALL umPrint(umMessage, src='createbc')
  WRITE(umMessage,'(A)') '[INFO] FRAMES mode selected'
  CALL umPrint(umMessage, src='createbc')
ELSE
  ! We're creating an LBC file
  ALLOCATE(lbcfile_type :: output_file)
  WRITE(umMessage,'(A)') '[INFO] Output file type is LBC File'
  CALL umPrint(umMessage, src='createbc')

  ! If the output file is an LBC, configure that appropriately
  ! These values go in the file header, so the file needs to know about them
  SELECT TYPE(output_file)
  CLASS IS (lbcfile_type)
    CALL output_file%set_rim_width(lbc_output_control%rim_width)
    WRITE(umMessage,'(A,I8)') '[INFO] Longitudinal Halo: ',              &
                                           lbc_output_control%halo_long
    CALL umPrint(umMessage, src='createbc')
    WRITE(umMessage,'(A,I8)') '[INFO] Latitudinal Halo: ',               &
                                           lbc_output_control%halo_lat
    CALL umPrint(umMessage, src='createbc')
    WRITE(umMessage,'(A,I8)') '[INFO] Rim width: ',                      &
                                           lbc_output_control%rim_width
    CALL umPrint(umMessage, src='createbc')
  END SELECT
END IF

! Check output file exists, print name, and open file for writing
WRITE(umMessage, '(A,A)') '[INFO] Output file is ',           &
                                    TRIM(lbc_output_control%output_filename)
CALL umPrint(umMessage, src='createbc')

INQUIRE(FILE=lbc_output_control%output_filename, EXIST=output_exists)
IF (output_exists) THEN
  CALL assign_file_unit(lbc_output_control%output_filename, tmp_unit, handler="fortran")
  OPEN(UNIT=tmp_unit, FILE=lbc_output_control%output_filename)
  WRITE(umMessage, '(A)') '[WARN] Output file exists. Deleting old version.'
  CALL umPrint(umMessage, src='createbc')
  CLOSE(UNIT=tmp_unit, STATUS='delete')
  CALL release_file_unit(tmp_unit, handler="fortran")
END IF

CALL output_file%open_file(lbc_output_control%output_filename)


! Initialise variables prior to loop over files

! Orography should be the first field in the first file, we track if we've
! found and processed it
processed_orography = .FALSE.

convert_dust_bins = .FALSE.

! Track to see if we've initialised the output file header using a valid
! input file as a template
copied_header = .FALSE.

input_file_format = detect_file_type(lbc_output_control%input_data(1))
SELECT CASE(input_file_format)
CASE (filetype_fieldsfile)
  WRITE(umMessage,'(A)') '[INFO] Input file format: fieldsfile'
  CALL umPrint(umMessage, src='createbc')
  ALLOCATE(fieldsfile_type :: input_file(lbc_output_control%num_input_files))
CASE (filetype_grib2)
  WRITE(umMessage,'(A)') '[INFO] Input file format: GRIB2'
  CALL umPrint(umMessage, src='createbc')
  icode = 11
  cmessage = 'GRIB2 input files not supported'
  CALL ereport(routinename, icode, cmessage)
CASE DEFAULT
  icode = 10
  cmessage = 'Unrecognised input file type'
  CALL ereport(routinename, icode, cmessage)
END SELECT

! ####################################
! 3. Loop over files 
! ####################################
DO input_file_num = 1, lbc_output_control%num_input_files

  IF (input_file_format /= detect_file_type(lbc_output_control%input_data(input_file_num))) THEN
    icode = 13
    cmessage = 'Inconsistent input file formats. All input files must have the same format.'
    CALL ereport(routinename, icode, cmessage)
  END IF

  ! 3a. Read file header
  CALL input_file(input_file_num)%allow_read_only(.TRUE.)
  CALL input_file(input_file_num)%open_file(lbc_output_control%input_data(input_file_num))
  CALL input_file(input_file_num)%read_header()

  IF (input_file_num == 1) THEN
    ! Determine if source and target are on the same rotation and set logical flags
    CALL check_pole_rotation(input_file(input_file_num)%pole_lat,                       &
                             input_file(input_file_num)%pole_long,                      &
                             lbc_output_control%pole_lat, lbc_output_control%pole_long, &
                             lbc_output_control%l_same_rotation,                        &
                             input_file(input_file_num)%l_source_rotated,               &
                             lbc_output_control%l_target_rotated)

    ! Do we require vertical interpolation?
    IF (lbc_output_control%horizontal_interpolation_method /= create_frame) THEN
      CALL check_vertical_interp(input_file(input_file_num), lbc_output_control, l_vert_interp)
    ELSE ! Will not perform vert interp if generating a frame
      l_vert_interp = .FALSE.
    END IF
    previous_file_p_grid = input_file(input_file_num)%p_grid
  ELSE
    ! Check that current file grid matches the previous file, otherwise interpolation
    ! weights will not be valid
    CALL input_file(input_file_num)%p_grid%compare_grids(previous_file_p_grid)
  END IF
  
  ! Create interpolation weights, only need to do this for first input file
  IF (input_file_num == 1) THEN 
    ! Allocate the interp_weights array for the five grid transformations:
    ! u->u, v->v, p->p, p->u and p->v and the three halo sizes: zero, single point
    ! and extended halos.
    ALLOCATE(interp_weights_type :: interp_weights_array(num_grid_types, num_halo_sizes))
    DO halo_code = 1,3 ! Loop over three halo types
      ! Call the calc_lbc_interp_weights method for each array element to calculate the 
      ! weights for each of the three supported horizontal grids.  The weight_index function takes
      ! the grid_code as an argument and returns the appropiate interp_weights_array index.
      CALL interp_weights_array(weights_index(u_points),halo_code)%calc_lbc_interp_weights(       &
                                                                   input_file(input_file_num),    &
                                                                   lbc_output_control, halo_code, &
                                                                   u_points, u_points)
      CALL interp_weights_array(weights_index(v_points),halo_code)%calc_lbc_interp_weights(       &
                                                                   input_file(input_file_num),    &
                                                                   lbc_output_control,            &
                                                                   halo_code, v_points, v_points)
      CALL interp_weights_array(weights_index(p_points),halo_code)%calc_lbc_interp_weights(       &
                                                                   input_file(input_file_num),    &
                                                                   lbc_output_control, halo_code, &
                                                                   p_points, p_points)
      ! Also calculate weights for interpolating the orography from P grid to U and V grid
      ! This is needed when vertically interpolating U and V fields
      CALL interp_weights_array(p_to_u_points, halo_code)%calc_lbc_interp_weights(       &
                                                          input_file(input_file_num),    &
                                                          lbc_output_control, halo_code, &
                                                          p_points, u_points)
      CALL interp_weights_array(p_to_v_points, halo_code)%calc_lbc_interp_weights(       &
                                                          input_file(input_file_num),    &
                                                          lbc_output_control, halo_code, &
                                                          p_points, v_points)
    END DO
    IF (.NOT. lbc_output_control%l_same_rotation) THEN
      ! If grids are not on the same rotation we will need to perform a rotation
      ! of the vector field quantities (i.e. U and V winds). Here we calculate the interpolation
      ! coefficients for the enlarged grids needed when rotating. Only needed for extended halos.
      ALLOCATE(interp_weights_type :: rotation_interp_weights_array(3))
      CALL rotation_interp_weights_array(u_to_enlarged_p)%calc_lbc_interp_weights(                &
                                                          input_file(input_file_num),             &
                                                          lbc_output_control, halo_type_extended, &
                                                          u_points, p_points,                     &
                                                          l_enlarged_p_grid = .TRUE.)
      CALL rotation_interp_weights_array(v_to_enlarged_p)%calc_lbc_interp_weights(                &
                                                          input_file(input_file_num),             &
                                                          lbc_output_control, halo_type_extended, &
                                                          v_points, p_points,                     &
                                                          l_enlarged_p_grid = .TRUE.)
      CALL rotation_interp_weights_array(p_to_enlarged_p)%calc_lbc_interp_weights(                &
                                                          input_file(input_file_num),             &
                                                          lbc_output_control, halo_type_extended, &
                                                          p_points, p_points,                     &
                                                          l_enlarged_p_grid = .TRUE.)
      ! Calculate the wind coefficients used by later routines to rotate the fields
      CALL calc_wind_rotation_coeff(input_file(input_file_num), lbc_output_control, wind_rotation_coeff)
      ! If we are creating a frame we want to determine the extent of the new LAM
      ! frame grid.
      IF (lbc_output_control%horizontal_interpolation_method == create_frame) THEN
        ! Pass in the rotation_interp_weights array to the calc_frame_grids routine if grids are not 
        ! on the same rotation.
        CALL calc_frame_grids(input_file(input_file_num), output_file, lbc_output_control, &
             interp_weights_array, rotation_interp_weights_array)
      END IF
    ELSE
      IF (lbc_output_control%horizontal_interpolation_method == create_frame) THEN
        ! If grids are on the same rotation then only need to pass the standard interp_weights array to
        ! the calc_frame_grids routine.
        CALL calc_frame_grids(input_file(input_file_num), output_file, lbc_output_control, &
             interp_weights_array)
      END IF
    END IF
  END IF ! input_file_num == 1

  ! 3b. Look for fields to interpolate
  CALL find_fields_to_interp(input_file(input_file_num), lbc_output_control)

  IF (found_fields_to_interp()) THEN
    IF (.NOT. copied_header) THEN
      ! Copy first input file header to LBC file as a starting point.
      CALL copy_file_header(input_file(input_file_num), output_file)
      copied_header = .TRUE.
      CALL update_file_header(output_file, lbc_output_control)
      ! Calculate size of header and set data_start so that we know where to write
      ! the first field on disk later
      CALL output_file%set_header_sizes()
    END IF

    WRITE(umMessage, '(A,I8,A,A)') '[INFO] Found ',                      &
           SIZE(fields_to_interp) , ' fields to interpolate in file ',   &
           TRIM(input_file(input_file_num)%filename)
    CALL umPrint(umMessage, src='createbc')

    ! ########################################
    ! 3c. Loop over fields in the current file
    ! ########################################

    DO j = 1, SIZE(fields_to_interp)

      ! The field we want to interpolate is k-th field in current input file
      k = fields_to_interp(j)
      
      IF (input_file(input_file_num)%fields(k)%has_stashcode(stashcode_orog) .AND. processed_orography) THEN
        CYCLE
      END IF
      
      ! Check grid for this field is consistent for FF
      SELECT TYPE(input_file)
        CLASS IS (fieldsfile_type)
          CALL input_file(input_file_num)%validate_grid(k)
      END SELECT
  
      WRITE(umMessage, '(A,I5,A,I4,A)') '[INFO] Retrieving field with STASH ',            &
            input_file(input_file_num)%fields(k)%quantity_ident, ' (',                    &
            input_file(input_file_num)%fields(k)%get_num_levels() ,' levels)'
      CALL umPrint(umMessage, src='createbc')

      IF (lbc_output_control%num_dust_fields > 0) THEN
        ! Check if the current field is a dust field in a fieldsfile. This is the 
        ! only supported way of converting between two- and six-bin dust.
        SELECT TYPE(input_file)
        CLASS IS (fieldsfile_type)
          IF (input_file(input_file_num)%fields(k)%has_stashcode(stashcode_dust1_mmr)) THEN
            ! If this is the first dust bin and the number of dust bins required is equal to the number
            ! of dust bins in the input file, no special handling is required, so convert_dust_bins remains
            ! FALSE.
            IF (input_file(input_file_num)%number_of_dust_bins_in_file() == lbc_output_control%num_dust_fields) THEN
              WRITE(umMessage, '(A)') '[INFO] No dust bin conversion required'
              CALL umPrint(umMessage, src='createbc')
            ELSE
              ! If this is the first dust bin and conversion is required, set the flag to TRUE.
              convert_dust_bins = .TRUE.
            END IF
          ELSE IF (input_file(input_file_num)%fields(k)%has_stashcode(stashcode_dust2_mmr) .AND. convert_dust_bins) THEN
            ! If we're doing dust bin conversion, this field is dust bin 2, and there are only two dust bins in the 
            ! file, the next field is non-dust and does not require special handling, so reset convert_dust_bins to FALSE
            IF (input_file(input_file_num)%number_of_dust_bins_in_file() == 2) THEN
              convert_dust_bins = .FALSE.
            END IF
            ! If we get to here we've already dealt with bin 2 as part of the dust bin conversion.
            CYCLE
            ! If we're doing dust bin conversion and this is a dust field in bins 3-6, we've already converted and 
            ! written those fields so we just ignore the field.
          ELSE IF (input_file(input_file_num)%fields(k)%has_stashcode(stashcode_dust3_mmr) .AND. convert_dust_bins) THEN
            CYCLE
          ELSE IF (input_file(input_file_num)%fields(k)%has_stashcode(stashcode_dust4_mmr) .AND. convert_dust_bins) THEN
            CYCLE
          ELSE IF (input_file(input_file_num)%fields(k)%has_stashcode(stashcode_dust5_mmr) .AND. convert_dust_bins) THEN
            CYCLE
          ELSE IF (input_file(input_file_num)%fields(k)%has_stashcode(stashcode_dust6_mmr) .AND. convert_dust_bins) THEN
            ! If this is bin 6, we ignore the field, but the next field is not dust, so we reset the flag to FALSE.
            convert_dust_bins = .FALSE.
            CYCLE
          END IF
        CLASS DEFAULT
          ! No dust for anything other than fieldsfiles
          convert_dust_bins = .FALSE.
        END SELECT
      END IF

      IF (lbc_output_control% horizontal_interpolation_method == create_frame) THEN
        ! Do not want to modify input fields if simply creating a frame
        convert_dust_bins = .FALSE.
      END IF
      
      IF (convert_dust_bins) THEN
        ! If we get this far as convert_dust_bins is TRUE, the current field is dust bin 1 and the next
        ! one or five fields in the input file are the remaining dust bins. The process_dust subroutine
        ! handles reading, converting, and writing the dust fields, so we call that. No further dust processing
        ! is required for this time.
        SELECT TYPE(input_file)
        CLASS IS (fieldsfile_type)

          CALL process_dust(j, lbc_output_control,                                      &
                            input_file(input_file_num), output_file, l_vert_interp,     &
                            interp_weights_array, orography_lbc_fields)

        END SELECT ! Check for fieldsfile input for dust
      ELSE
        ! Non-dust field, or dust not requiring bin conversion
        
        ! If rotations of the wind fields are necessary, call wind processing code.  This will 
        ! perform all rotations and interpolations to the LBC grid and write the output to file.
        !
        ! It is necessary to process both U and V components of the wind at the same time, as
        ! both components are needed when rotating.  We are currently in a loop over all fields
        ! in the current file that are to be interpolated, these fields have been sorted to 
        ! match the order of the fields specified in the input namelist.  It is a requirement that
        ! the v (advected v) field follows immediately after the u (advected u) field in the list
        ! of stash codes in the namelist. When the loop encounters a u (advected u) field and the
        ! input and output grids are not on the same rotation the process_winds routine will be called.
        ! This routine will also process the next field, and check it is the equivalent v component.
        ! On the next interation of the loop the v (advected v) we can cycle the loop as the v
        ! component has already been processed.
        IF (.NOT. lbc_output_control%l_same_rotation) THEN
          IF (input_file(input_file_num)%fields(k)%has_stashcode(stashcode_u) .OR.                     &
              input_file(input_file_num)%fields(k)%has_stashcode(stashcode_u_adv)) THEN
            CALL process_winds(j, lbc_output_control, input_file(input_file_num),                      &
                 output_file, orography_enlarged_grid_lbc_field, orography_lbc_fields,                 &
                 wind_rotation_coeff, interp_weights_array, rotation_interp_weights_array, l_vert_interp)
            CYCLE
          END IF
          IF (input_file(input_file_num)%fields(k)%has_stashcode(stashcode_v) .OR.                     &
              input_file(input_file_num)%fields(k)%has_stashcode(stashcode_v_adv)) THEN
            CYCLE
          END IF
        END IF
        
        ! 3c i. Read field
        CALL input_file(input_file_num)%read_field(k)
        ! Check that number of levels is consistent with first and last level numbers
        CALL input_file(input_file_num)%fields(k)%grid%vert_grid%check_num_levels()

        ! 3c.ii Process orography - orography must appear first in the namelist of stash items
        IF (.NOT. processed_orography .AND. input_file(input_file_num)%fields(k)%quantity_ident ==     &
                stashcode_orog) THEN
          CALL process_orography(input_file(input_file_num), lbc_output_control, orography_lbc_fields, &
                                 interp_weights_array, k)
          IF (.NOT. lbc_output_control%l_same_rotation .AND.                                           &
               lbc_output_control%l_target_rotated) THEN
            ! If wind rotations are applied to LBC U and V fields, will need enlarged P grid orography
            CALL process_orography_enlarged_grid(input_file(input_file_num), lbc_output_control,       &
                                  orography_enlarged_grid_lbc_field, rotation_interp_weights_array, k)
          END IF
          processed_orography = .TRUE.
          WRITE(umMessage,'(A)') '[INFO] Saving orography fields'
          CALL umPrint(umMessage, src='createbc')
        END IF

        ! ##################################################
        ! 3c iii. Perform the interpolation / frames cutout
        ! ##################################################

        ! Interpolate the kth input field and add the output LBC field object returned 
        ! by interpolation control to the output file object.  The add_field function
        ! then returns the new_field_num integer which is used to reference the new
        ! field in the array of fields contained in the output_file.
        
        output_field = interp_control(lbc_output_control,                         &
                                      input_file(input_file_num),                 &
                                      interp_weights_array, orography_lbc_fields, &
                                      l_vert_interp, field_number=k)
        new_field_num = output_file%add_field(output_field)

        WRITE(umMessage,'(A,I5,A,A)') '[INFO] Generating field for STASH ',                            & 
             output_file%fields(new_field_num)%quantity_ident, ' at time ',                            &
             output_file%fields(new_field_num)%return_validity_time_string()
        CALL umPrint(umMessage, src='createbc')

        !   Once we've got the interpolated field, we no longer need the original 
        !   in memory.
        CALL input_file(input_file_num)%fields(k)%unload_data()

        IF (lbc_output_control% horizontal_interpolation_method == create_frame) THEN
          ! If creating a frame need to move the field onto the smaller frame grid
          CALL update_frame_field_grid(output_file, new_field_num)
        ELSE
          ! Perform post interpolation transformation
          CALL post_interp_transform(output_file%fields(new_field_num),                      &
               lbc_output_control%q_min)
        END IF

        IF (lbc_output_control%horizontal_interpolation_method == create_frame &
             .AND. lbc_output_control%frames_packing_option /= packing_unchanged ) THEN

          SELECT TYPE(output_file)
            
            TYPE IS(fieldsfile_type)
              CALL output_file%modify_packing_method(new_field_num,       &
                                      lbc_output_control%frames_packing_option) 
              IF (.NOT. packing_warned) THEN
                WRITE(umMessage,'(A,I0)') '[INFO] Changing packing of Frames to packing option code ', &
                                           lbc_output_control%frames_packing_option
                CALL umPrint(umMessage, src='createbc')
                packing_warned = .TRUE.
              END IF
            CLASS DEFAULT
              IF (.NOT. packing_warned) THEN
                WRITE(umMessage,'(A)') '[WARN] Changing packing of Frames not supported for this file type'
                CALL umPrint(umMessage, src='createbc')
                packing_warned = .TRUE.
              END IF
          END SELECT
        END IF

        ! 3c. iv  Write field
        CALL output_file%write_field(new_field_num)
        fields_written = .TRUE.

        ! 3c. v. Tidy memory
        CALL output_file%fields(new_field_num)%unload_data()

        ! 3c. vi. Write header - if we do it here, for the cost of a little IO 
        !                        we get a vaguely useful looking file if 
        !                        something causes createbc to abort.
        IF (.NOT. lbc_output_control%write_header_only_once) THEN
          CALL output_file%write_header()
        END IF

      END IF ! Doing a dust field requiring bin conversion or not
    END DO ! Loop over fields
    ! Clear array now we've done all interpolation required
    CALL clear_fields_to_interp()

  ELSE ! no fields to interpolate
    WRITE(umMessage,'(A)') '[WARN] No fields to interpolate found in file'
    CALL umPrint(umMessage, src='createbc')
  END IF ! Any fields to interpolate in this file?

  DEALLOCATE (input_file(input_file_num)%fields)
  CALL input_file(input_file_num)%close_file()

END DO ! Loop over files

! 4. Finalise
!  Write header of output file to disk.
IF (fields_written) THEN
  CALL output_file%write_header()
ELSE
  cmessage = 'No fields written to disk.'
  icode = 999
  CALL ereport(routinename, icode, cmessage)
END IF

!  Tidy up output file and deallocate
CALL output_file%close_file()
DEALLOCATE(input_file)

WRITE(umMessage,'(A)') '[DONE] CreateBC completed.'
CALL umPrint(umMessage, src='createbc')

CALL timer( 'createbc', 2)

CALL appTerminate()

! Final top-level DrHook caliper, then prevent any further calls to DrHook.
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
CALL drhook_control_disable()

END PROGRAM createbc
