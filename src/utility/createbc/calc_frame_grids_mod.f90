! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE calc_frame_grids_mod
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.

! DrHook Modules
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

PRIVATE
PUBLIC :: calc_frame_grids

CHARACTER(LEN=*), PARAMETER :: ModuleName = 'CALC_FRAME_GRIDS_MOD'

CONTAINS

SUBROUTINE calc_frame_grids(input_file, output_file, lbc_output_control, interp_weights, &
                            rotation_interp_weights)
! Description:
!  Calculate the size of the reduced input domain.  This will contain only the data needed 
!  to produce an LBC for the given output domain.
!
USE datafile_mod,              ONLY: datafile_type, endgame, new_dynamics
USE interp_weights_mod,        ONLY: interp_weights_type, p_to_enlarged_p, &
                                     weights_index, u_to_enlarged_p, v_to_enlarged_p, &
                                     num_halo_sizes, num_grid_types
USE lbc_output_control_mod,    ONLY: lbc_output_control_type
USE stashmaster_constants_mod, ONLY: p_points, u_points, v_points
USE process_orography_mod,     ONLY: p_to_u_points, p_to_v_points
USE um_parparams,              ONLY: halo_type_extended
USE ereport_mod,               ONLY: ereport
USE errormessagelength_mod,    ONLY: errormessagelength
USE missing_data_mod,          ONLY: imdi

IMPLICIT NONE
CLASS(datafile_type),           INTENT(INOUT) :: input_file
CLASS(datafile_type),           INTENT(INOUT) :: output_file
TYPE(lbc_output_control_type), INTENT(INOUT) :: lbc_output_control
CLASS(interp_weights_type),     INTENT(INOUT) :: interp_weights(num_grid_types, num_halo_sizes)
CLASS(interp_weights_type), OPTIONAL, INTENT(INOUT) :: rotation_interp_weights(3)
! Error reporting variables
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER  :: routinename = 'CALC_FRAME_GRIDS'
! Local variables
INTEGER :: source_east_col, source_west_col, source_north_row, source_south_row
INTEGER :: i,j,m, halo_code
LOGICAL :: rotation_weights_present = .FALSE.
INTEGER :: lam_num_rows_p, lam_num_cols_p
INTEGER :: lam_num_rows_u, lam_num_cols_u
INTEGER :: lam_num_rows_v, lam_num_cols_v
REAL, ALLOCATABLE :: lam_p_longitudes(:)
REAL, ALLOCATABLE :: lam_p_latitudes(:)
REAL, ALLOCATABLE :: lam_u_longitudes(:)
REAL, ALLOCATABLE :: lam_u_latitudes(:)
REAL, ALLOCATABLE :: lam_v_longitudes(:)
REAL, ALLOCATABLE :: lam_v_latitudes(:)
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Need to determine the extent of the largest frames mask. i.e the most southerly/northerly
! westerly/easterly points. Rather than work with latitudes and longitudes we can just
! use the row and column numbers.  Once we have determined these points the P, U and V
! grids can be calculated for the new cutout LAM domain.
IF (PRESENT(rotation_interp_weights)) rotation_weights_present = .TRUE.


! Use the P grid to establish the size of the cutout
IF (rotation_weights_present) THEN
  ! Use enlarged p grid if wind rotation is being performed
  CALL calc_domain_extent(rotation_interp_weights(p_to_enlarged_p)%frame_mask, source_east_col,  &
       source_west_col, source_north_row, source_south_row)
ELSE
  ! Use standard P grid
  CALL calc_domain_extent(interp_weights(weights_index(p_points),halo_type_extended)%frame_mask, &
       source_east_col, source_west_col, source_north_row, source_south_row)
END IF

! Increase the size of the cutout for safety allow overide from namelist variables if user needs to 
! extend or reduce the cutout area. Please note this will not change the size of the frames mask, it will
! simply increase or reduce the size of the missing data area surrounding the mask.

! If the variables have not been set in the namelist then extend the cutout by 1 row/col
IF (lbc_output_control%frames_cutout_adjust_east == imdi) lbc_output_control%frames_cutout_adjust_east = 1
IF (lbc_output_control%frames_cutout_adjust_west == imdi) lbc_output_control%frames_cutout_adjust_west = 1
IF (lbc_output_control%frames_cutout_adjust_north == imdi) lbc_output_control%frames_cutout_adjust_north = 1
IF (lbc_output_control%frames_cutout_adjust_south == imdi) lbc_output_control%frames_cutout_adjust_south = 1

source_east_col  = source_east_col  + lbc_output_control%frames_cutout_adjust_east
source_west_col  = source_west_col  - lbc_output_control%frames_cutout_adjust_west
source_north_row = source_north_row + lbc_output_control%frames_cutout_adjust_north
source_south_row = source_south_row - lbc_output_control%frames_cutout_adjust_south

! Get P grid number of columns and rows
lam_num_rows_p = source_north_row - source_south_row + 1
! Calculation of the number of columns will depend on whether the domain crosses the meridian
! if this is the case then the eastern extent of the domain will have a lower column number
! than the western extent
IF (source_west_col > source_east_col) THEN
  lam_num_cols_p = source_east_col + (input_file%p_grid%get_num_cols() - source_west_col + 1)
ELSE
  lam_num_cols_p = source_east_col - source_west_col  + 1
END IF
! U grid
lam_num_cols_u = lam_num_cols_p
lam_num_rows_u = lam_num_rows_p
! V grid
lam_num_cols_v = lam_num_cols_p
IF (input_file%grid_staggering == endgame) THEN
  lam_num_rows_v = lam_num_rows_p + 1
ELSE IF (input_file%grid_staggering == new_dynamics) THEN
  lam_num_rows_v = lam_num_rows_p - 1
ELSE
  icode = 10
  cmessage = "Unsupported grid staggering"
  CALL ereport(routinename, icode, cmessage)
END IF

ALLOCATE(lam_p_latitudes(lam_num_rows_p))
ALLOCATE(lam_p_longitudes(lam_num_cols_p))
ALLOCATE(lam_u_latitudes(lam_num_rows_u))
ALLOCATE(lam_u_longitudes(lam_num_cols_u))
ALLOCATE(lam_v_latitudes(lam_num_rows_v))
ALLOCATE(lam_v_longitudes(lam_num_cols_v))

! Update the output file p grid
output_file%p_grid = input_file%p_grid
j = source_south_row
DO i = 1, lam_num_rows_p
  lam_p_latitudes(i) = output_file%p_grid%horiz_grid%latitudes(j)
  j = j + 1
END DO
j = source_west_col
DO i = 1, lam_num_cols_p
  IF (j > output_file%p_grid%get_num_cols()) j = 1 ! Crosses meridian
  lam_p_longitudes(i) = output_file%p_grid%horiz_grid%longitudes(j)
  j = j + 1
END DO
! Use MOVE_ALLOC to reallocate the file grid lat/long and populate with the newly
! calculated values.
CALL MOVE_ALLOC(lam_p_latitudes, output_file%p_grid%horiz_grid%latitudes)
CALL MOVE_ALLOC(lam_p_longitudes, output_file%p_grid%horiz_grid%longitudes)
CALL output_file%p_grid%set_num_rows(lam_num_rows_p)
CALL output_file%p_grid%set_num_cols(lam_num_cols_p)
! Non wrapping LAM
CALL output_file%p_grid%set_hemisphere_indicator(3)

! Update the output file u grid
output_file%u_grid = input_file%u_grid
j = source_south_row
DO i = 1, lam_num_rows_u
  lam_u_latitudes(i) = output_file%u_grid%horiz_grid%latitudes(j)
  j = j + 1
END DO
j = source_west_col
DO i = 1, lam_num_cols_u
  IF (j > output_file%u_grid%get_num_cols()) j = 1 ! Crosses meridian
  lam_u_longitudes(i) = output_file%u_grid%horiz_grid%longitudes(j)
  j = j + 1
END DO
! Use MOVE_ALLOC to reallocate the file grid lat/long and populate with the newly
! calculated values.
CALL MOVE_ALLOC(lam_u_latitudes, output_file%u_grid%horiz_grid%latitudes)
CALL MOVE_ALLOC(lam_u_longitudes, output_file%u_grid%horiz_grid%longitudes)
CALL output_file%u_grid%set_num_rows(lam_num_rows_u)
CALL output_file%u_grid%set_num_cols(lam_num_cols_u)
! Non wrapping LAM
CALL output_file%u_grid%set_hemisphere_indicator(3)

! Update the output file v grid
output_file%v_grid = input_file%v_grid
j = source_south_row
DO i = 1, lam_num_rows_v
  lam_v_latitudes(i) = output_file%v_grid%horiz_grid%latitudes(j)
  j = j + 1
END DO
j = source_west_col
DO i = 1, lam_num_cols_v
  IF (j > output_file%v_grid%get_num_cols()) j = 1 ! Crosses meridian
  lam_v_longitudes(i) = output_file%v_grid%horiz_grid%longitudes(j)
  j = j + 1
END DO
! Use MOVE_ALLOC to reallocate the file grid lat/long and populate with the newly
! calculated values.
CALL MOVE_ALLOC(lam_v_latitudes, output_file%v_grid%horiz_grid%latitudes)
CALL MOVE_ALLOC(lam_v_longitudes, output_file%v_grid%horiz_grid%longitudes)
CALL output_file%v_grid%set_num_rows(lam_num_rows_v)
CALL output_file%v_grid%set_num_cols(lam_num_cols_v)
! Non wrapping LAM
CALL output_file%v_grid%set_hemisphere_indicator(3)

! Need to add the start row/col and number of rows and columns to each 
! interp_weights_type.  When the frames masks is being applied to a field
! these values can be used to cutout the lam domain.
DO halo_code = 1, 3
  CALL interp_weights(weights_index(p_points),halo_code)%set_frames_domain( &
       source_south_row, source_west_col, lam_num_rows_p, lam_num_cols_p)
  CALL interp_weights(weights_index(u_points),halo_code)%set_frames_domain( &
       source_south_row, source_west_col, lam_num_rows_u, lam_num_cols_u)
  CALL interp_weights(weights_index(v_points),halo_code)%set_frames_domain( &
       source_south_row, source_west_col, lam_num_rows_v, lam_num_cols_v)
  CALL interp_weights(p_to_u_points,halo_code)%set_frames_domain(           &
       source_south_row, source_west_col, lam_num_rows_p, lam_num_cols_p)
  CALL interp_weights(p_to_v_points,halo_code)%set_frames_domain(           &
       source_south_row, source_west_col, lam_num_rows_p, lam_num_cols_p)
END DO
IF (rotation_weights_present) THEN
  CALL rotation_interp_weights(u_to_enlarged_p)%set_frames_domain(          &
       source_south_row, source_west_col, lam_num_rows_u, lam_num_cols_u)
  CALL rotation_interp_weights(v_to_enlarged_p)%set_frames_domain(          &
       source_south_row, source_west_col, lam_num_rows_v, lam_num_cols_v)
  CALL rotation_interp_weights(p_to_enlarged_p)%set_frames_domain(          &
       source_south_row, source_west_col, lam_num_rows_p, lam_num_cols_p)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE calc_frame_grids

!-------------------------------------------------------------------------------

SUBROUTINE calc_domain_extent(frame_mask, east_col, west_col, north_row, south_row)
! Determine the min and max column which contains data
IMPLICIT NONE
LOGICAL, INTENT(IN) :: frame_mask(:,:)
INTEGER, INTENT(OUT) :: east_col
INTEGER, INTENT(OUT) :: west_col
INTEGER, INTENT(OUT) :: north_row
INTEGER, INTENT(OUT) :: south_row
! Local vars
INTEGER :: num_cols 
INTEGER :: num_rows 
INTEGER :: i
CHARACTER(LEN=*), PARAMETER :: routinename = 'CALC_DOMAIN_EXTENT'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

num_cols = SIZE(frame_mask, 1)
num_rows = SIZE(frame_mask, 2)

IF (ANY(frame_mask(num_cols,:)) .AND.  ANY(frame_mask(1,:)) ) THEN
  ! In this case the frames domain crosses the east/west boundary of input domain. 
  ! Most easterly point on LAM domain  will have lower column number than the most 
  ! westerly point.
  !
  ! Examine each column starting from the western edge of the grid and exit the loop
  ! when the column does not contain any TRUE values this means the previous column 
  ! will be the easterly extent of the LAM domain
  east_col = 1
  DO i = 1, num_cols
    IF (.NOT. ANY(frame_mask(i,:))) THEN
      EXIT
    END IF
  END DO
  east_col = i - 1
  ! Examine each column starting from the eastern edge of the grid and exit the loop
  ! when the column does not contain any TRUE values this means the previous column 
  ! will be the westernly extent of the LAM domain
  west_col = num_cols
  DO i = num_cols, 1, -1
    IF (.NOT. ANY(frame_mask(i,:))) THEN
      EXIT
    END IF
  END DO
  west_col = i + 1
ELSE
  ! In this case the most easterly point on LAM domain will have higher column number than 
  ! the most westerly point.
  !
  ! Examine each column starting from the western edge of the grid and exit the loop
  ! when the column does not contain any FALSE values this means that this column 
  ! will be the westernly extent of the LAM domain
  west_col = 1
  DO i = 1, num_cols
    IF (ANY(frame_mask(i,:))) THEN
      EXIT
    END IF
  END DO
  west_col = i
  ! Examine each column starting from the eastern edge of the grid and exit the loop
  ! when the column does not contain any FALSE values this means that this column 
  ! will be the westernly extent of the LAM domain
  east_col = num_cols
  DO i = num_cols, 1, -1
    IF (ANY(frame_mask(i,:))) THEN
      EXIT
    END IF
  END DO
  east_col = i
END IF

! Examine each row starting from the southern edge of the domain. Once the row contains 
! any true values exit the loop this means the row will be the southerly extent of the
! LAM domain
south_row = 1
DO i = 1, num_rows
  IF (ANY(frame_mask(:,i))) THEN
    EXIT
  END IF
END DO
south_row = i
! Examine each row starting from the northern edge of the domain. Once the row contains 
! any true values exit the loop this means the row will be the northernly extent of the
! LAM domain
north_row = num_rows
DO i = num_rows, 1, -1
  IF (ANY(frame_mask(:,i))) THEN
    EXIT
  END IF
END DO
north_row = i

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE calc_domain_extent

END MODULE calc_frame_grids_mod
