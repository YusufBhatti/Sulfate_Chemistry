! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE lbc_output_control_mod

USE filenamelength_mod,         ONLY: filenamelength
USE missing_data_mod,           ONLY: rmdi, imdi
USE three_dimensional_grid_mod, ONLY: three_dimensional_grid_type
USE file_mod,                   ONLY: file_type
USE fieldsfile_constants_mod,   ONLY: new_dynamics, endgame
! DrHook Modules
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

! Description:
!    Module to contain data and routines related to the lbc_grid namelist.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.

! Fixed parameters
INTEGER, PARAMETER :: max_stash_items = 100
INTEGER, PARAMETER :: max_input_files = 1000
INTEGER, PARAMETER :: max_fields_per_file = 10000

! Parameters for horizontal_interpolation_method
INTEGER, PARAMETER :: create_frame = 0
INTEGER, PARAMETER :: bilinear_interp = 2

! Parameters for packing
INTEGER, PARAMETER :: packing_unchanged = -1 ! As input field
REAL,    PARAMETER :: packing_unpacked_accuracy = -99.0

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'LBC_OUTPUT_CONTROL_MOD'

TYPE :: lbc_output_control_type

  INTEGER, PUBLIC :: num_input_files = 0
  INTEGER, PUBLIC :: num_stash_items = 0
  INTEGER, PUBLIC :: num_dust_fields = 0

  ! Namelist variables
  REAL    :: dlat                 = rmdi    ! Latitudinal spacing
  REAL    :: dlong                = rmdi    ! Longitudinal spacing
  INTEGER :: end_time(6)          = imdi    ! End time for LBC generation
  INTEGER :: horizontal_interpolation_method = imdi ! Horiz interp method
  INTEGER :: nlat                 = imdi    ! Number of latitudinal rows
  INTEGER :: nlong                = imdi    ! Number of longitudinal columns
  INTEGER :: num_levels           = imdi    ! Number of vertical levels
  INTEGER :: output_grid_stagger  = imdi    ! Output grid staggering, 
                                            !        3 = ND, 6 = EG
  REAL    :: pole_lat             = rmdi    ! Latitude of rotated pole
  REAL    :: pole_long            = rmdi    ! Longitude of rotated pole
  REAL    :: start_lat            = rmdi    ! Latitude of first p-point
  REAL    :: start_long           = rmdi    ! Longitude of first p-point
  INTEGER :: start_time(6)        = imdi    ! Start time for LBC generation
  INTEGER :: stash_codes(max_stash_items) = imdi ! Array of stash items 
  LOGICAL :: variable_resolution  = .FALSE. !  var-res output flag
  INTEGER :: vertical_interpolation_method = imdi ! Vertical interp method
  INTEGER :: halo_lat             = imdi    ! Halo size in y direction
  INTEGER :: halo_long            = imdi    ! Halo size in x direction
  INTEGER :: rim_width            = imdi    ! Rim width
  REAL    :: q_min                = rmdi    ! Reset to q_min if specific 
                                            !    humidity value is below q_min
  INTEGER :: num_reserved_headers = imdi    ! Number of lookup headers in 
                                            !    output file
  LOGICAL :: l_same_rotation      = .FALSE. ! Input and output on same rotation
  LOGICAL :: l_target_rotated     = .FALSE. ! Target is on rotated grid
  LOGICAL :: write_header_only_once = .FALSE. ! Write the header once (true) 
                                              ! or after every field (false)

  INTEGER :: frames_cutout_adjust_north = imdi ! Adjust the mdi region surrounding
  INTEGER :: frames_cutout_adjust_south = imdi ! the cutout frame
  INTEGER :: frames_cutout_adjust_east  = imdi
  INTEGER :: frames_cutout_adjust_west  = imdi
  INTEGER :: frames_packing_option = imdi ! New packing code for a Frame
  
  INTEGER :: time_interval        = imdi    ! Time interval between LBCs in 
                                            ! seconds

  ! Array of dust stash items - stored separately to handle dust-bin conversion
  INTEGER :: dust_stash_codes(max_stash_items) = imdi 

  ! Original order of STASH codes for things which are sensitive to 
  ! order in output file. This is corrected to turn into LBC section 31 STASH.
  INTEGER :: target_stash_codes(max_stash_items) = imdi 

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

  ! Grid objects, array of length three for the three halo codes
  TYPE(three_dimensional_grid_type) :: p_grid(3), v_grid(3), u_grid(3), p_grid_enlarged

  CONTAINS
  
  PROCEDURE, PASS :: process
  PROCEDURE, PASS :: count_input_files
  PROCEDURE, PASS :: check_limits
  PROCEDURE, PASS :: generate_grids

END TYPE lbc_output_control_type


CONTAINS

!-------------------------------------------------------------------------------

SUBROUTINE process(this)
! Using the contents of the namelist, derive other auxiliary quantities

USE transform_stash_order_mod, ONLY: initialise_master_stash
IMPLICIT NONE
CLASS(lbc_output_control_type), INTENT(INOUT) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'PROCESS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL this%check_limits()
CALL this%generate_grids()
CALL this%count_input_files()
CALL initialise_master_stash()
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE process

!-------------------------------------------------------------------------------

SUBROUTINE count_input_files(this)
! Count the number of input files, stash codes, and dust stash codes in the
! namelist

USE errormessagelength_mod, ONLY: errormessagelength
USE ereport_mod, ONLY: ereport
IMPLICIT NONE
CLASS(lbc_output_control_type), INTENT(INOUT) :: this
INTEGER :: i
CHARACTER(LEN=errormessagelength) :: cmessage
INTEGER :: icode
CHARACTER(LEN=*), PARAMETER :: routinename='COUNT_INPUT_FILES'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Count the number of input files
this%num_input_files = max_input_files
DO i = 1, max_input_files
  IF (this%input_data(i) == 'unset' .OR. TRIM(this%input_data(i)) == '') THEN
    this%num_input_files = i - 1
    EXIT
  END IF
END DO

! Abort if no input files are specified
IF (this%num_input_files == 0 ) THEN
  icode = 1
  cmessage = 'No input files specified.'
  CALL ereport(routinename, icode, cmessage)
END IF

! Count the number of non-dust STASH items

! Initialise to max_stash_items for the (unlikely) case where the number of
! stash codes equals the maximum, so the IF condition will never be true and
! num_stash_items would otherwise be unset.
this%num_stash_items = max_stash_items
DO i = 1, max_stash_items
  IF (this%stash_codes(i) == imdi) THEN
    this%num_stash_items = i - 1
    EXIT
  END IF
END DO

! Count the number of dust STASH items

! Initialise to max_stash_items for the same reason as above.
this%num_dust_fields = max_stash_items
DO i = 1, max_stash_items
  IF (this%dust_stash_codes(i) == imdi) THEN
    this%num_dust_fields = i - 1
    EXIT
  END IF
END DO
IF (this%num_dust_fields == max_stash_items) THEN
  this%num_dust_fields = 0
END IF

! Abort if no stash codes are specified
IF (this%num_dust_fields + this%num_stash_items == 0 ) THEN
  icode = 1
  cmessage = 'No STASH codes specified.'
  CALL ereport(routinename, icode, cmessage)
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE count_input_files

!-------------------------------------------------------------------------------

SUBROUTINE check_limits(this)
! If halo_lat, halo_long and/or rim_width are inappropriately set, abort.
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE

CLASS(lbc_output_control_type), INTENT(INOUT) :: this

INTEGER :: icode
CHARACTER(LEN=*), PARAMETER :: routinename='CHECK_LIMITS'
CHARACTER(LEN=errormessagelength) :: cmessage

! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (this%halo_long < 0 .OR. this%halo_lat < 0 .OR. this%rim_width < 0) THEN
  cmessage = 'Haloes and rim width must be zero or positive.'
  icode = 140
  CALL ereport(routinename, icode, cmessage)
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE check_limits

!-------------------------------------------------------------------------------

SUBROUTINE generate_grids(this)

USE vertnamelist_mod,          ONLY: vertlevs, eta_theta, eta_rho, z_top_of_model,   &
                                     first_constant_r_rho_level
USE vrhoriz_grid_mod,          ONLY: horizgrid, lambda_input_p, lambda_input_u,      &
                                     phi_input_p, phi_input_v
USE file_manager,              ONLY: assign_file_unit, release_file_unit
USE um_parparams,              ONLY: halo_type_single, halo_type_extended, halo_type_no_halo
USE stashmaster_constants_mod, ONLY: p_points, u_points, v_points
USE ereport_mod,               ONLY: ereport
USE errormessagelength_mod,    ONLY: errormessagelength

IMPLICIT NONE
CLASS(lbc_output_control_type), INTENT(INOUT) :: this

INTEGER :: vert_unit_num, var_res_unit_num ! unit numbers for namelist IO
INTEGER :: i,j ! loop counter
INTEGER :: halo_code ! Halo code loop
REAL :: lat_p, long_p, lat_u, long_u, lat_v, long_v
INTEGER :: nlat_halo, nlong_halo, adjust_v_nlat, adjust_u_nlong
REAL, ALLOCATABLE :: lambda_halo_p(:), lambda_halo_u(:)
REAL, ALLOCATABLE :: phi_halo_v(:), phi_halo_p(:)
REAL, ALLOCATABLE :: temp_latitude(:), temp_longitude(:) 
REAL :: edge_spacing_lat, edge_spacing_long
REAL :: start_lat_halo, start_long_halo
INTEGER :: local_halo_long, local_halo_lat

CHARACTER(LEN=errormessagelength) :: cmessage
INTEGER :: icode
CHARACTER(LEN=*), PARAMETER :: routinename = 'GENERATE_GRIDS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Will calculate the lbc grids three times, once for each halo type
! This is seen as preferable to creating just the extended halo grid
! and then using offsets when working with fields on single point or
! zero point halos
DO halo_code = 1,3 ! Loop over halo_type_single, halo_type_extended, and no halo
  IF (halo_code == halo_type_single) THEN 
    local_halo_long = 1 
    local_halo_lat  = 1 
  ELSE IF (halo_code == halo_type_extended) THEN 
    local_halo_long = this%halo_long
    local_halo_lat  = this%halo_lat
  ELSE IF (halo_code == halo_type_no_halo) THEN 
    local_halo_long = 0
    local_halo_lat  = 0
  END IF

  IF (this%variable_resolution) THEN
    CALL assign_file_unit(this%horizontal_grid_file, var_res_unit_num, handler="fortran")
    lambda_input_p = rmdi
    phi_input_p = rmdi
    OPEN(UNIT=var_res_unit_num, FILE=this%horizontal_grid_file, &
         ACTION='READ')
    READ(var_res_unit_num, horizgrid)
    CLOSE(var_res_unit_num)
    CALL release_file_unit(var_res_unit_num, handler="fortran")

    DO i = 1, SIZE(lambda_input_p)
      IF (lambda_input_p(i) == rmdi) THEN
        this%nlong = i - 1
        EXIT
      END IF
    END DO
    
    DO i = 1, SIZE(phi_input_p)
      IF (phi_input_p(i) == rmdi) THEN
        this%nlat = i - 1
        EXIT
      END IF
    END DO
  END IF
  
  nlong_halo = this%nlong + local_halo_long * 2 
  nlat_halo = this%nlat + local_halo_lat * 2 
  
  ! Start lat and long for first row and column of grid - regular grid including halo
  ! note that for ND this is a P row and P column but for EG this is a V row and a 
  ! U column.
  start_lat_halo = this%start_lat - this%dlat * local_halo_lat
  start_long_halo = this%start_long - this%dlong * local_halo_long

  IF (this%output_grid_stagger == new_dynamics) THEN
    ! For ND grid using current lbc behaviour the U grid has one less column
    ! and V grid has one less row compared to P grid
    adjust_u_nlong  = -1
    adjust_v_nlat   = -1
    ! Calculate offset from start lat long to first P, U and V points.
    lat_p  = start_lat_halo
    long_p = start_long_halo
    lat_u  = start_lat_halo
    long_u = start_long_halo + this%dlong/2.0
    lat_v  = start_lat_halo + this%dlat/2.0
    long_v = start_long_halo
  ELSE IF (this%output_grid_stagger == endgame) THEN
    ! For EG grid, V has an extra row compared to P grid
    adjust_u_nlong  =  0
    adjust_v_nlat   =  1
    ! Calculate offset from start lat long to first P, U and V points.
    lat_p  = start_lat_halo + this%dlat/2.0
    long_p = start_long_halo + this%dlong/2.0
    lat_u  = start_lat_halo + this%dlat/2.0
    long_u = start_long_halo
    lat_v  = start_lat_halo
    long_v = start_long_halo + this%dlong/2.0
  ELSE
    WRITE(cmessage, '(A,I8,A)') 'Horizontal grid type ', this%output_grid_stagger, &
         'not known'
    icode = 2
    CALL ereport(routinename, icode, cmessage)
  END IF

  IF (this%variable_resolution) THEN

    ! Allocate to maximum possible size 
    ALLOCATE(lambda_halo_p(nlong_halo))
    ALLOCATE(lambda_halo_u(nlong_halo + adjust_u_nlong))
    ALLOCATE(phi_halo_p(nlat_halo))
    ALLOCATE(phi_halo_v(nlat_halo + adjust_v_nlat))

    lambda_halo_p = rmdi
    lambda_halo_u = rmdi
    phi_halo_p = rmdi
    phi_halo_v = rmdi

    ! Copy central region excluding halos
    lambda_halo_p(local_halo_long + 1:local_halo_long + this%nlong) =  &
         lambda_input_p(1:this%nlong)
    ! U columns may need adjusting compared to P grid
    lambda_halo_u(local_halo_long + 1:local_halo_long + this%nlong + adjust_u_nlong) =  &
         lambda_input_u(1:this%nlong + adjust_u_nlong)
    phi_halo_p(local_halo_lat + 1:local_halo_lat + this%nlat) =        &
         phi_input_p(1:this%nlat)
    ! V rows need adjusting compared to P grid
    phi_halo_v(local_halo_lat + 1:local_halo_lat + this%nlat + adjust_v_nlat) = &
         phi_input_v(1:this%nlat + adjust_v_nlat)

    ! Fill edges by extrapolation
    edge_spacing_long = lambda_input_p(2) - lambda_input_p(1)
    edge_spacing_lat = phi_input_p(2) - phi_input_p(1)

    IF (local_halo_long > 0) THEN
      ! Western edge
      DO i = local_halo_long, 1, -1
        lambda_halo_p(i) = lambda_halo_p(i+1) - edge_spacing_long
        lambda_halo_u(i) = lambda_halo_u(i+1) - edge_spacing_long
      END DO

      ! Eastern edge
      DO i = local_halo_long + this%nlong + 1, nlong_halo
        lambda_halo_p(i) = lambda_halo_p(i-1) + edge_spacing_long
      END DO
      ! Adjust to account for number of U columns
      DO i = local_halo_long + this%nlong + adjust_u_nlong + 1, nlong_halo + adjust_u_nlong
        lambda_halo_u(i) = lambda_halo_u(i-1) + edge_spacing_long
      END DO
    END IF

    IF (local_halo_lat > 0) THEN
      ! Southern edge
      DO i = local_halo_lat, 1, -1
        phi_halo_p(i) = phi_halo_p(i+1) - edge_spacing_lat
        phi_halo_v(i) = phi_halo_v(i+1) - edge_spacing_lat
      END DO

      ! Northern edge
      DO i = local_halo_lat + this%nlat + 1, nlat_halo
        phi_halo_p(i) = phi_halo_p(i-1) + edge_spacing_lat
      END DO
      ! Adjust to account for number of V rows
      DO i = local_halo_lat + this%nlat + adjust_v_nlat + 1, nlat_halo + adjust_v_nlat
        phi_halo_v(i) = phi_halo_v(i-1) + edge_spacing_lat
      END DO
    END IF

    CALL this%p_grid(halo_code)%horiz_grid%set_defined_horizontal_grid( &
         nx      =  nlong_halo,                                         &
         ny      =  nlat_halo,                                          &
         list_x  =  lambda_halo_p(1:nlong_halo),                        &
         list_y  =  phi_halo_p(1:nlat_halo))
    CALL this%u_grid(halo_code)%horiz_grid%set_defined_horizontal_grid( &
         nx      =  nlong_halo + adjust_u_nlong,                        &
         ny      =  nlat_halo,                                          &
         list_x  =  lambda_halo_u(1:nlong_halo),                        &
         list_y  =  phi_halo_p(1:nlat_halo))
    CALL this%v_grid(halo_code)%horiz_grid%set_defined_horizontal_grid( &
         nx      =  nlong_halo,                                         &
         ny      =  nlat_halo + adjust_v_nlat,                          &
         list_x  =  lambda_halo_p(1:nlong_halo),                        &
         list_y  =  phi_halo_v(1:nlat_halo))
  
    DEALLOCATE(phi_halo_v)
    DEALLOCATE(phi_halo_p)
    DEALLOCATE(lambda_halo_u)
    DEALLOCATE(lambda_halo_p)

  ELSE ! Regular grid
    CALL this%p_grid(halo_code)%horiz_grid%set_regular_horizontal_grid( &
         nx     =  nlong_halo,                                          &
         ny     =  nlat_halo,                                           &
         startx =  long_p,                                              &
         starty =  lat_p,                                               &
         dx     =  this%dlong,                                          & 
         dy     =  this%dlat)
    CALL this%u_grid(halo_code)%horiz_grid%set_regular_horizontal_grid( &
         nx      =  nlong_halo + adjust_u_nlong,                        &
         ny      =  nlat_halo,                                          &
         startx  =  long_u,                                             &
         starty  =  lat_u,                                              &
         dx      =  this%dlong,                                         &
         dy=this%dlat)
    CALL this%v_grid(halo_code)%horiz_grid%set_regular_horizontal_grid( &
         nx      =  nlong_halo,                                         &
         ny      =  nlat_halo + adjust_v_nlat,                          &
         startx  =  long_v,                                             &
         starty  =  lat_v,                                              &
         dx      =  this%dlong,                                         &
         dy      =  this%dlat)
  END IF


  CALL assign_file_unit(this%vertical_levels_file, vert_unit_num, handler="fortran")
  OPEN(UNIT=vert_unit_num, FILE=this%vertical_levels_file, ACTION='READ')
  READ(vert_unit_num, vertlevs)
  ! Read in eta_theta, eta_rho, first_constant_r_rho_level, z_top_of_model
  CLOSE(vert_unit_num)
  CALL release_file_unit(vert_unit_num, handler="fortran")

  CALL this%p_grid(halo_code)%vert_grid%set_defined_vertical_grid(      &
       this%num_levels,                                                 &
       num_model_levels              =  this%num_levels,                &
       first_rho_of_constant_height  =  first_constant_r_rho_level,     &
       height_at_top_theta_level     =  z_top_of_model,                 &
       list_eta_theta                =  eta_theta(1:this%num_levels+1), &
       list_eta_rho                  =  eta_rho(1:this%num_levels))
  CALL this%u_grid(halo_code)%vert_grid%set_defined_vertical_grid(      &
       this%num_levels,                                                 &
       num_model_levels              =  this%num_levels,                &
       first_rho_of_constant_height  =  first_constant_r_rho_level,     &
       height_at_top_theta_level     =  z_top_of_model,                 &
       list_eta_theta                =  eta_theta(1:this%num_levels+1), &
       list_eta_rho                  =  eta_rho(1:this%num_levels))
  CALL this%v_grid(halo_code)%vert_grid%set_defined_vertical_grid(      &
       this%num_levels,                                                 &
       num_model_levels              =  this%num_levels,                &
       first_rho_of_constant_height  =  first_constant_r_rho_level,     &
       height_at_top_theta_level     =  z_top_of_model,                 &
       list_eta_theta                =  eta_theta(1:this%num_levels+1), &
       list_eta_rho                  =  eta_rho(1:this%num_levels))

  ! Set grid type
  CALL this%p_grid(halo_code)%set_horiz_grid_type(p_points)
  CALL this%u_grid(halo_code)%set_horiz_grid_type(u_points)
  CALL this%v_grid(halo_code)%set_horiz_grid_type(v_points)
  ! Set pole lat and long
  CALL this%p_grid(halo_code)%set_pole(this%pole_lat, this%pole_long)
  CALL this%u_grid(halo_code)%set_pole(this%pole_lat, this%pole_long)
  CALL this%v_grid(halo_code)%set_pole(this%pole_lat, this%pole_long)
  ! Set the halo size in the horizontal grids
  CALL this%p_grid(halo_code)%set_halo_ns(local_halo_lat)
  CALL this%u_grid(halo_code)%set_halo_ns(local_halo_lat)
  CALL this%v_grid(halo_code)%set_halo_ns(local_halo_lat)
  CALL this%p_grid(halo_code)%set_halo_ew(local_halo_long)
  CALL this%u_grid(halo_code)%set_halo_ew(local_halo_long)
  CALL this%v_grid(halo_code)%set_halo_ew(local_halo_long)

END DO ! End loop over halo codes

! Enlarged P grid is needed when rotating winds
IF (this%output_grid_stagger == new_dynamics) THEN
  ! New Dynamics P grid surrounds U and V points
  this%p_grid_enlarged = this%p_grid(halo_type_extended)
ELSE IF (this%output_grid_stagger == endgame) THEN
  ! Need an extra row on north and south and extra column on west
  ! First take the standard p grid as a starting point
  this%p_grid_enlarged = this%p_grid(halo_type_extended)
  ALLOCATE(temp_latitude(this%p_grid_enlarged%get_num_rows() + 2))
  ALLOCATE(temp_longitude(this%p_grid_enlarged%get_num_cols() + 1))
  ! Fill in known value of points from original arrays
  DO i = 2, this%p_grid_enlarged%get_num_cols() + 1
    temp_longitude(i) = this%p_grid_enlarged%horiz_grid%longitudes(i-1)
  END DO
  DO j = 2, this%p_grid_enlarged%get_num_rows() + 1
    temp_latitude(j) = this%p_grid_enlarged%horiz_grid%latitudes(j-1)
  END DO
  ! Calculate extra points
  temp_longitude(1) = temp_longitude(2) - ( temp_longitude(3) - temp_longitude(2))
  temp_latitude(1) = temp_latitude(2) - ( temp_latitude(3) - temp_latitude(2))
  temp_latitude(this%p_grid_enlarged%get_num_rows() + 2) =            &
       temp_latitude(this%p_grid_enlarged%get_num_rows()+1)           &
       + (  temp_latitude(this%p_grid_enlarged%get_num_rows()+1) -    &
       temp_latitude(this%p_grid_enlarged%get_num_rows())         )
  ! Now update the grid object lat/long arrays
  CALL  this%p_grid_enlarged%set_num_rows( this%p_grid_enlarged%get_num_rows() + 2 )
  CALL  this%p_grid_enlarged%set_num_cols( this%p_grid_enlarged%get_num_cols() + 1 )
  CALL this%p_grid_enlarged%set_latitudes(temp_latitude)
  CALL this%p_grid_enlarged%set_longitudes(temp_longitude)
  DEALLOCATE(temp_longitude)
  DEALLOCATE(temp_latitude)
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE generate_grids

!-------------------------------------------------------------------------------

END MODULE lbc_output_control_mod
