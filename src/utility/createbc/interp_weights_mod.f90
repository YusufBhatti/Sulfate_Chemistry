! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE interp_weights_mod
! DrHook Modules
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

! Description:
!   A base class to hold interpolation weights and frames mask
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.

INTEGER, PUBLIC, PARAMETER :: num_halo_sizes = 3
INTEGER, PUBLIC, PARAMETER :: num_grid_types = 5

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'INTERP_WEIGHTS_MOD'

TYPE, PUBLIC :: interp_weights_type

  PRIVATE
  
  ! Indices of the input points which contribute to the new
  ! interpolated output value.  Only calculate bottom left/right
  ! as can access the top left/right by adding a row
  INTEGER,  PUBLIC, ALLOCATABLE :: bilinear_index_b_l(:)
  INTEGER,  PUBLIC, ALLOCATABLE :: bilinear_index_b_r(:)
  INTEGER,  PUBLIC, ALLOCATABLE :: bilinear_index_t_l(:)
  INTEGER,  PUBLIC, ALLOCATABLE :: bilinear_index_t_r(:)
  ! Weights of the four input points which contribute to the new interpolated
  ! value
  REAL,     PUBLIC, ALLOCATABLE :: weight_t_r(:)
  REAL,     PUBLIC, ALLOCATABLE :: weight_b_r(:)
  REAL,     PUBLIC, ALLOCATABLE :: weight_t_l(:)
  REAL,     PUBLIC, ALLOCATABLE :: weight_b_l(:)
  ! Mask array used when generating a frame. Should be the same size as the 
  ! single level input field.  TRUE if the input point is used by the interpolation.
  LOGICAL,  PUBLIC, ALLOCATABLE :: frame_mask(:,:)
  INTEGER,  PUBLIC              :: frame_start_row
  INTEGER,  PUBLIC              :: frame_start_col
  INTEGER,  PUBLIC              :: frame_num_rows
  INTEGER,  PUBLIC              :: frame_num_cols

  CONTAINS
  PROCEDURE, PUBLIC, PASS :: set_frames_domain
  PROCEDURE, PUBLIC, PASS :: get_frames_domain
  PROCEDURE, PUBLIC, PASS :: calc_lbc_interp_weights

END TYPE interp_weights_type

INTEGER, PUBLIC, PARAMETER :: u_to_enlarged_p = 1 ! Parameters used to index the interpolation coeffs
INTEGER, PUBLIC, PARAMETER :: v_to_enlarged_p = 2 ! used when working with enlarged grids for rotations
INTEGER, PUBLIC, PARAMETER :: p_to_enlarged_p = 3 

CONTAINS

!-------------------------------------------------------------------------------

SUBROUTINE get_frames_domain(this, start_row, start_col, num_rows, num_cols)
IMPLICIT NONE
CLASS(interp_weights_type),   INTENT(IN) :: this
INTEGER, INTENT(OUT) :: start_row
INTEGER, INTENT(OUT) :: start_col
INTEGER, INTENT(OUT) :: num_rows
INTEGER, INTENT(OUT) :: num_cols
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_FRAMES_DOMAIN'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

start_row = this%frame_start_row
start_col = this%frame_start_col
num_cols  = this%frame_num_cols
num_rows  = this%frame_num_rows

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE get_frames_domain

!-------------------------------------------------------------------------------

SUBROUTINE set_frames_domain(this, start_row, start_col, num_rows, num_cols)
IMPLICIT NONE
CLASS(interp_weights_type),   INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: start_row
INTEGER, INTENT(IN) :: start_col
INTEGER, INTENT(IN) :: num_rows
INTEGER, INTENT(IN) :: num_cols
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_FRAMES_DOMAIN'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

this%frame_start_row = start_row     
this%frame_start_col = start_col 
this%frame_num_cols  = num_cols  
this%frame_num_rows  = num_rows

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_frames_domain

!-------------------------------------------------------------------------------

INTEGER FUNCTION weights_index(grid_type)
USE stashmaster_constants_mod, ONLY: u_points, v_points, p_points
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE
INTEGER, INTENT(IN) :: grid_type
INTEGER            :: icode = 0
CHARACTER(LEN=errormessagelength) :: cmessage 
CHARACTER(LEN=*), PARAMETER  :: routinename = 'WEIGHTS_INDEX'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! An array of interpolation weight types is generated in the main createbc
! program.  Each index refers to a set of weights for a grid type. This
! function will map the grid type to the index so that the grid type is all
! that is needed to access the correct weights.
IF (grid_type == u_points) THEN
  weights_index = 1
ELSE IF (grid_type == v_points) THEN
  weights_index = 2
ELSE IF (grid_type == p_points) THEN
  weights_index = 3
ELSE
  icode = 55
  WRITE(cmessage, '(A,I8)') "Cannot ascertain index of weights"// &
       " array for grid type: ", grid_type
  CALL ereport(routinename, icode, cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION weights_index

!-------------------------------------------------------------------------------

SUBROUTINE calc_lbc_interp_weights(this, input_file, lbc_output_control, halo_code, input_grid_type, &
                                   output_grid_type, l_enlarged_p_grid)

USE lbc_output_control_mod, ONLY: lbc_output_control_type, create_frame
USE three_dimensional_grid_mod, ONLY: three_dimensional_grid_type
USE datafile_mod, ONLY: datafile_type
USE stashmaster_constants_mod,  ONLY: u_points, v_points, p_points
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE calc_lbc_coords_mod, ONLY: calc_lbc_coords
USE lam_inclusion_mod, ONLY: lam_inclusion

USE fort2c_interfaces, ONLY: logical_to_bool
USE latlon_eq_rotation_mod, ONLY: rotate_eq_to_latlon, rotate_latlon_to_eq
USE f_shum_horizontal_field_interp_mod, ONLY :                                 &
                               f_shum_horizontal_field_bi_lin_interp_get_coeffs
USE umPrintMgr,   ONLY: umPrint, umMessage

IMPLICIT NONE

CLASS(interp_weights_type),     INTENT(INOUT) :: this
CLASS(datafile_type),           TARGET, INTENT(IN) :: input_file
TYPE(lbc_output_control_type), TARGET, INTENT(IN) :: lbc_output_control
INTEGER,           INTENT(IN)  :: halo_code
INTEGER,           INTENT(IN)  :: input_grid_type
INTEGER,           INTENT(IN)  :: output_grid_type
LOGICAL, OPTIONAL, INTENT(IN)  :: l_enlarged_p_grid

TYPE(three_dimensional_grid_type), POINTER :: lbc_grid
TYPE(three_dimensional_grid_type), POINTER :: input_grid
REAL, ALLOCATABLE :: latitude_lbc(:)    ! Latitude points for LBC grid
REAL, ALLOCATABLE :: longitude_lbc(:)   ! Longitude points for LBC grid
REAL, ALLOCATABLE :: latitude_temp(:)   ! Latitude points from source grid temp array
REAL, ALLOCATABLE :: longitude_temp(:)  ! Longitude points from source grid temp array
INTEGER           :: lbc_size
LOGICAL           :: enlarged_p_grid

INTEGER            :: icode = 0
CHARACTER(LEN=errormessagelength) :: cmessage 
CHARACTER(LEN=*), PARAMETER  :: routinename = 'CALC_LBC_INTERP_WEIGHTS'
INTEGER, PARAMETER :: global = 0
INTEGER :: i, corner_index
INTEGER :: num_rows, num_cols, col, row, l 
INTEGER :: index_bl, index_br, index_tl, index_tr
INTEGER :: index_array(4)
REAL :: lbc_phi_min
REAL :: lbc_phi_max
REAL :: grid_tolerance
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Determine which grid the interpolation weights are being calculated for
! Point at LBC grid based on grid code and halo code
! Point at input grid based on grid tyoe
IF ( input_grid_type == u_points ) THEN
  input_grid => input_file%u_grid
ELSE IF ( input_grid_type == v_points ) THEN
  input_grid => input_file%v_grid
ELSE IF ( input_grid_type == p_points ) THEN
  input_grid => input_file%p_grid
ELSE
  icode = 10
  WRITE(cmessage, '(A,I8)') "Unsupported input grid type: ", input_grid_type
  CALL ereport(routinename, icode, cmessage)
END IF

IF ( output_grid_type == u_points ) THEN
  lbc_grid   => lbc_output_control%u_grid(halo_code)
ELSE IF ( output_grid_type == v_points ) THEN
  lbc_grid   => lbc_output_control%v_grid(halo_code)
ELSE IF ( output_grid_type == p_points ) THEN
  lbc_grid   => lbc_output_control%p_grid(halo_code)
ELSE
  icode = 20
  WRITE(cmessage, '(A,I8)') "Unsupported output grid type: ", output_grid_type
  CALL ereport(routinename, icode, cmessage)
END IF

! Check if the interpolation weights are being calculated for the enlarged p grid
IF ( PRESENT(l_enlarged_p_grid) ) THEN
  enlarged_p_grid = l_enlarged_p_grid
ELSE
  enlarged_p_grid = .FALSE.
END IF

IF (enlarged_p_grid) THEN
  ! If rotating the U and V winds need to interpolate to the enlarged P grid. This is 
  ! to ensure that there is enough data to linearly interpolate from the enlarged P grid
  ! safely back to the U and V grid. Will also need the weights to be generated for moving 
  ! from standard P grid to enlarged P grid in order to generate orography used in vertical 
  !interpolation of the enlarged field grids.
  lbc_grid   => lbc_output_control%p_grid_enlarged
END IF

! Get the lbc target grid information and calculate lat and long for each LBC point
CALL calc_lbc_coords(enlarged_p_grid, lbc_output_control, lbc_grid, &
                     latitude_lbc, longitude_lbc, lbc_size)

! Need to perform rotations if source and lbc grid are not
! on same rotation. Will use latitude_temp and longitude_temp to avoid
! array aliasing in argument list
IF (.NOT. lbc_output_control%l_same_rotation) THEN
  ALLOCATE( longitude_temp (lbc_size) )
  ALLOCATE( latitude_temp    (lbc_size) )
  IF (lbc_output_control%l_target_rotated) THEN
    ! Rotate the LBC points onto the standard poles 90, 0
    CALL rotate_eq_to_latlon(latitude_lbc, longitude_lbc,   &
                             latitude_temp, longitude_temp, &
         lbc_output_control%pole_lat, lbc_output_control%pole_long, lbc_size)
    ! Copy back the rotated values from rotate_eq_to_latlon
    latitude_lbc    = latitude_temp
    longitude_lbc   = longitude_temp
  END IF
  IF (input_file%l_source_rotated) THEN
    ! If the source grid is not on a standard pole rotation will 
    ! need to rotate the LBC points to this rotation
    CALL rotate_latlon_to_eq(latitude_lbc, longitude_lbc,   &
                             latitude_temp, longitude_temp, &
         input_file%pole_lat, input_file%pole_long, lbc_size)
    latitude_lbc  = latitude_temp
    longitude_lbc = longitude_temp
  END IF
  DEALLOCATE( latitude_temp  )
  DEALLOCATE( longitude_temp )
END IF


! Perform spatial consistency check to validate that the target domain is
! entirely with the source domain. 
IF (.NOT. input_grid%get_hemisphere_indicator() == global) THEN
  lbc_phi_min = MINVAL(latitude_lbc)
  lbc_phi_max = MAXVAL(latitude_lbc)
  grid_tolerance = SQRT(EPSILON(1.0))
  ! The halo size is included in the num_rows and num_cols, so pass in halo_y/x as zero
  CALL lam_inclusion(src_rows=input_grid%get_num_rows(),          src_row_len=input_grid%get_num_cols(),       &
                     src_halo_x=0,                                src_halo_y=0,                                &
                     src_lambda=input_grid%get_longitudes(),      src_phi=input_grid%get_latitudes(),          &
                     src_pole_lat=input_grid%get_pole_latitude(), src_pole_lon=input_grid%get_pole_longitude(),& 
                     targ_rows=lbc_grid%get_num_rows(),           targ_row_len=lbc_grid%get_num_cols(),        &
                     targ_halo_x=0,                               targ_halo_y=0,                               &
                     targ_lambda=lbc_grid%get_longitudes(),       targ_phi=lbc_grid%get_latitudes(),           &
                     targ_phi_min=lbc_phi_min,                    targ_phi_max=lbc_phi_max,                    &
                     targ_pole_lat=lbc_grid%get_pole_latitude(),  targ_pole_lon=lbc_grid%get_pole_longitude(), &  
                     l_same_rotation=lbc_output_control%l_same_rotation,                                       &
                     l_src_rotated=input_file%l_source_rotated,   grid_tol_arg=grid_tolerance)
END IF
! Allocate the weight and index arrays of the interp_weights_type
ALLOCATE( this%bilinear_index_b_l(lbc_size) ) 
ALLOCATE( this%bilinear_index_b_r(lbc_size) ) 
ALLOCATE( this%bilinear_index_t_l(lbc_size) ) 
ALLOCATE( this%bilinear_index_t_r(lbc_size) ) 
ALLOCATE( this%weight_t_r(lbc_size) ) 
ALLOCATE( this%weight_b_r(lbc_size) ) 
ALLOCATE( this%weight_t_l(lbc_size) ) 
ALLOCATE( this%weight_b_l(lbc_size) )


! Call general utility routine to calculate weights
CALL f_shum_horizontal_field_bi_lin_interp_get_coeffs                      &
              ( this%bilinear_index_b_l,        this%bilinear_index_b_r,   &
                this%bilinear_index_t_l,        this%bilinear_index_t_r,   &
                this%weight_t_r,                this%weight_b_r,           &
                this%weight_t_l,                this%weight_b_l,           &
                input_grid%get_longitudes(),    input_grid%get_latitudes(),&
                longitude_lbc,                  latitude_lbc,              &
                input_grid%get_num_cols(),      input_grid%get_num_rows(), &
                lbc_size,                       logical_to_bool(           &
                                                  input_file%l_source_cyclic))

IF (lbc_output_control%horizontal_interpolation_method == create_frame) THEN
  num_cols = input_grid%get_num_cols()
  num_rows = input_grid%get_num_rows()
  ALLOCATE(this%frame_mask(num_cols, num_rows))
  this%frame_mask = .FALSE.
  DO l = 1, SIZE(this%bilinear_index_b_l)
    ! Set mask array to TRUE for index of input data points which would have been 
    ! used by bilinear interpolation 
    ! The indices stored in bilinear_index_b_l/r assume that the field is a 1D array. These
    ! need to be converted to their 2D counterparts.        
    index_bl = this%bilinear_index_b_l(l)             
    index_br = this%bilinear_index_b_r(l)             
    index_tl = this%bilinear_index_t_l(l)
    index_tr = this%bilinear_index_t_r(l)
    ! For each point which contributes to the interpolation we should also set the frames
    ! mask to TRUE for the four surrounding points. These may be needed in the case of 
    ! input wind rotation
    index_array = [index_bl, index_br, index_tl, index_tr]
    DO corner_index = 1, 4
      row = CEILING(REAL(index_array(corner_index))/REAL(num_cols))
      col = index_array(corner_index) - ((row - 1) * num_cols)
      this%frame_mask(col,row-1:row+1)   = .TRUE.
      IF (col == num_cols) THEN ! Wrap round to first column
        IF (input_file%l_source_cyclic) THEN
          this%frame_mask(1,row-1:row+1) = .TRUE.
        ELSE
          icode = 30
          WRITE(cmessage, '(A)') "Target domain too close to edge of source domain"
          CALL ereport(routinename, icode, cmessage)
        END IF
      ELSE
        this%frame_mask(col+1,row-1:row+1) = .TRUE.
      END IF
      IF (col == 1) THEN ! Wrap round to last column
        IF (input_file%l_source_cyclic) THEN
          this%frame_mask(num_cols,row-1:row+1) = .TRUE.
        ELSE
          icode = 40
          WRITE(cmessage, '(A)') "Target domain too close to edge of source domain"
          CALL ereport(routinename, icode, cmessage)
        END IF
      ELSE
        this%frame_mask(col-1,row-1:row+1) = .TRUE.
      END IF
    END DO
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE calc_lbc_interp_weights

!-------------------------------------------------------------------------------

END MODULE interp_weights_mod
