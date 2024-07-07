! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE three_dimensional_grid_mod

USE horizontal_lat_long_grid_mod, ONLY: horizontal_lat_long_grid_type
USE vertical_grid_mod,            ONLY: vertical_grid_type
! DrHook Modules
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

! Description:
!   A class to describe a three-dimensional grid, made up of a 2D horiz grid
!   and a 1D vertical grid. Further extension to this class could allow 
!   different shaped grid types, but for the time being we restrict this to
!   that case.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                        &
  ModuleName = 'THREE_DIMENSIONAL_GRID_MOD'

TYPE, PUBLIC :: three_dimensional_grid_type
  PRIVATE
  TYPE(horizontal_lat_long_grid_type), PUBLIC :: horiz_grid
  TYPE(vertical_grid_type), PUBLIC :: vert_grid
  CONTAINS
  PROCEDURE, PASS :: compare_grids
  PROCEDURE, PASS :: get_grid_dimensions
  PROCEDURE, PASS :: get_horizontal_coordinates
  PROCEDURE, PASS :: get_horiz_grid_type
  PROCEDURE, PASS :: get_latitudes
  PROCEDURE, PASS :: get_longitudes
  PROCEDURE, PASS :: set_latitudes
  PROCEDURE, PASS :: set_longitudes
  PROCEDURE, PASS :: get_regular_horizontal_grid
  PROCEDURE, PASS :: get_vertical_levels
  PROCEDURE, PASS :: get_vertical_level_values
  PROCEDURE, PASS :: get_grid_code
  PROCEDURE, PASS :: get_hemisphere_indicator
  PROCEDURE, PASS :: get_halo_code
  PROCEDURE, PASS :: get_halo_ns
  PROCEDURE, PASS :: get_halo_ew
  PROCEDURE, PASS :: get_num_rows
  PROCEDURE, PASS :: get_num_cols
  PROCEDURE, PASS :: get_num_levels
  PROCEDURE, PASS :: get_num_model_levels
  PROCEDURE, PASS :: get_pole
  PROCEDURE, PASS :: set_pole
  PROCEDURE, PASS :: get_pole_latitude
  PROCEDURE, PASS :: set_pole_latitude
  PROCEDURE, PASS :: get_pole_longitude
  PROCEDURE, PASS :: set_pole_longitude
  PROCEDURE, PASS :: set_regular_horizontal_grid
  PROCEDURE, PASS :: set_defined_horizontal_grid
  PROCEDURE, PASS :: set_horiz_grid_type
  PROCEDURE, PASS :: set_num_levels
  PROCEDURE, PASS :: set_num_model_levels
  PROCEDURE, PASS :: set_grid_code
  PROCEDURE, PASS :: set_hemisphere_indicator
  PROCEDURE, PASS :: set_halo_code
  PROCEDURE, PASS :: set_halo_ns
  PROCEDURE, PASS :: set_halo_ew
  PROCEDURE, PASS :: set_num_rows
  PROCEDURE, PASS :: set_num_cols

END TYPE three_dimensional_grid_type

CONTAINS

!-------------------------------------------------------------------------------

SUBROUTINE compare_grids(this, compare_grid)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(IN) :: this
CLASS(three_dimensional_grid_type), INTENT(IN) :: compare_grid
CHARACTER(LEN=*), PARAMETER :: routinename = 'COMPARE_GRIDS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL this%vert_grid%compare_grids(compare_grid%vert_grid)
CALL this%horiz_grid%compare_grids(compare_grid%horiz_grid)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE compare_grids

!-------------------------------------------------------------------------------

SUBROUTINE get_grid_dimensions(this, nx, ny, nz)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(IN) :: this
INTEGER, INTENT(OUT) :: nx, ny, nz
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_GRID_DIMENSIONS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
nz = this%get_num_levels()
ny = this%get_num_rows()
nx = this%get_num_cols()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE get_grid_dimensions

!-------------------------------------------------------------------------------

INTEGER FUNCTION get_horiz_grid_type(this)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(INOUT) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_HORIZ_GRID_TYPE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

get_horiz_grid_type = this%horiz_grid%get_grid_type()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION get_horiz_grid_type

!-------------------------------------------------------------------------------

SUBROUTINE get_horizontal_coordinates(this, latitudes, longitudes)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(IN) :: this
REAL, ALLOCATABLE, INTENT(OUT) :: latitudes(:), longitudes(:)
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_HORIZONTAL_COORDINATES'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL this%horiz_grid%get_horizontal_coordinates(latitudes, longitudes)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE get_horizontal_coordinates

!-------------------------------------------------------------------------------

FUNCTION get_latitudes(this) RESULT(latitudes)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(IN) :: this
REAL, ALLOCATABLE :: local_latitudes(:), latitudes(:)
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_LATITUDES'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

local_latitudes = this%horiz_grid%get_latitudes()
CALL MOVE_ALLOC(local_latitudes, latitudes)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION get_latitudes

!-------------------------------------------------------------------------------

FUNCTION get_longitudes(this) RESULT(longitudes)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(IN) :: this
REAL, ALLOCATABLE  :: local_longitudes(:), longitudes(:)
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_LONGITUDES'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

local_longitudes = this%horiz_grid%get_longitudes()
CALL MOVE_ALLOC(local_longitudes, longitudes)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION get_longitudes

!-------------------------------------------------------------------------------

SUBROUTINE set_longitudes(this, longitudes)
IMPLICIT NONE 
CLASS(three_dimensional_grid_type), INTENT(INOUT) :: this
REAL, INTENT(IN)    :: longitudes(:)
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_LONGITUDES'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL this%horiz_grid%set_longitudes(longitudes)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_longitudes

!-------------------------------------------------------------------------------

SUBROUTINE set_latitudes(this, latitudes)
IMPLICIT NONE 
CLASS(three_dimensional_grid_type), INTENT(INOUT) :: this
REAL, INTENT(IN)    :: latitudes(:)
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_LATITUDES'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL this%horiz_grid%set_latitudes(latitudes)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_latitudes

!-------------------------------------------------------------------------------

SUBROUTINE get_regular_horizontal_grid(this, nx, ny, startx, starty, dx, dy)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(INOUT) :: this
INTEGER, INTENT(OUT) :: nx, ny
REAL, INTENT(OUT) :: startx, starty, dx, dy
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_REGULAR_HORIZONTAL_GRID'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL this%horiz_grid%get_regular_horizontal_grid(nx, ny, startx, starty, dx, dy)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE get_regular_horizontal_grid

!-------------------------------------------------------------------------------

INTEGER FUNCTION get_grid_code(this)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(IN) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_GRID_CODE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

get_grid_code = this%horiz_grid%get_grid_code()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION get_grid_code

!-------------------------------------------------------------------------------

INTEGER FUNCTION get_hemisphere_indicator(this)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(IN) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_HEMISPHERE_INDICATOR'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

get_hemisphere_indicator = this%horiz_grid%get_hemisphere_indicator()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION get_hemisphere_indicator

!-------------------------------------------------------------------------------

INTEGER FUNCTION get_halo_code(this)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(IN) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_HALO_CODE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

get_halo_code = this%horiz_grid%get_halo_code()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION get_halo_code

!-------------------------------------------------------------------------------

INTEGER FUNCTION get_halo_ns(this)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(IN) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_HALO_NS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

get_halo_ns = this%horiz_grid%get_halo_ns()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION get_halo_ns

!-------------------------------------------------------------------------------

INTEGER FUNCTION get_halo_ew(this)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(IN) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_HALO_EW'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

get_halo_ew = this%horiz_grid%get_halo_ew()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION get_halo_ew

!-------------------------------------------------------------------------------

INTEGER FUNCTION get_num_rows(this)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(IN) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_NUM_ROWS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

get_num_rows = this%horiz_grid%get_num_rows()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION get_num_rows

!-------------------------------------------------------------------------------

INTEGER FUNCTION get_num_cols(this)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(IN) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_NUM_COLS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

get_num_cols = this%horiz_grid%get_num_cols()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION get_num_cols

!-------------------------------------------------------------------------------

INTEGER FUNCTION get_num_levels(this)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(IN) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_NUM_LEVELS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

get_num_levels = this%vert_grid%get_num_levels()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION get_num_levels

!-------------------------------------------------------------------------------

INTEGER FUNCTION get_num_model_levels(this)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(IN) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_NUM_MODEL_LEVELS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

get_num_model_levels = this%vert_grid%get_num_model_levels()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION get_num_model_levels

!-------------------------------------------------------------------------------

SUBROUTINE get_vertical_levels(this, levels)
USE missing_data_mod, ONLY: imdi
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(IN) :: this
INTEGER, INTENT(INOUT) :: levels(:)
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_VERTICAL_LEVELS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (ALLOCATED(this%vert_grid%levels)) THEN
  levels = this%vert_grid%levels
ELSE
  levels = imdi
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE get_vertical_levels

!-------------------------------------------------------------------------------

SUBROUTINE get_vertical_level_values(this, level_values)
USE missing_data_mod, ONLY: rmdi
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(IN) :: this
REAL, INTENT(INOUT) :: level_values(:)
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_VERTICAL_LEVEL_VALUES'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (ALLOCATED(this%vert_grid%levels)) THEN
  level_values = this%vert_grid%level_values
ELSE
  level_values = rmdi
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE get_vertical_level_values

!-------------------------------------------------------------------------------

SUBROUTINE get_pole(this, pole_latitude, pole_longitude)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(IN) :: this
REAL, INTENT(OUT) :: pole_latitude, pole_longitude
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_POLE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL this%horiz_grid%get_pole(pole_latitude, pole_longitude)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE get_pole

!-------------------------------------------------------------------------------

FUNCTION get_pole_latitude(this)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(IN) :: this
REAL :: get_pole_latitude
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_POLE_LATITUDE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

get_pole_latitude = this%horiz_grid%get_pole_latitude()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION get_pole_latitude

!-------------------------------------------------------------------------------

FUNCTION get_pole_longitude(this)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(IN) :: this
REAL :: get_pole_longitude
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_POLE_LONGITUDE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

get_pole_longitude = this%horiz_grid%get_pole_longitude()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION get_pole_longitude

!-------------------------------------------------------------------------------

SUBROUTINE set_horiz_grid_type(this, grid_type)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: grid_type
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_HORIZ_GRID_TYPE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL this%horiz_grid%set_grid_type(grid_type)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_horiz_grid_type

!-------------------------------------------------------------------------------

SUBROUTINE set_num_levels(this, num_levels)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: num_levels
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_NUM_LEVELS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL this%vert_grid%set_num_levels(num_levels)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_num_levels

!-------------------------------------------------------------------------------

SUBROUTINE set_num_model_levels(this, num_model_levels)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: num_model_levels
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_NUM_MODEL_LEVELS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL this%vert_grid%set_num_model_levels(num_model_levels)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_num_model_levels

!-------------------------------------------------------------------------------

SUBROUTINE set_grid_code(this, grid_code)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: grid_code
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_GRID_CODE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL this%horiz_grid%set_grid_code(grid_code)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_grid_code

!-------------------------------------------------------------------------------

SUBROUTINE set_hemisphere_indicator(this, hemisphere_indicator)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: hemisphere_indicator
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_HEMISPHERE_INDICATOR'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL this%horiz_grid%set_hemisphere_indicator(hemisphere_indicator)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_hemisphere_indicator

!-------------------------------------------------------------------------------

SUBROUTINE set_halo_code(this, halo_code)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: halo_code
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_HALO_CODE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL this%horiz_grid%set_halo_code(halo_code)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_halo_code

!-------------------------------------------------------------------------------

SUBROUTINE set_halo_ns(this, halo_ns)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: halo_ns
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_HALO_NS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL this%horiz_grid%set_halo_ns(halo_ns)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_halo_ns

!-------------------------------------------------------------------------------

SUBROUTINE set_halo_ew(this, halo_ew)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: halo_ew
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_HALO_EW'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL this%horiz_grid%set_halo_ew(halo_ew)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_halo_ew

!-------------------------------------------------------------------------------

SUBROUTINE set_num_rows(this, num_rows)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: num_rows
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_NUM_ROWS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL this%horiz_grid%set_num_rows(num_rows)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_num_rows

!-------------------------------------------------------------------------------

SUBROUTINE set_num_cols(this, num_cols)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: num_cols
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_NUM_COLS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL this%horiz_grid%set_num_cols(num_cols)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_num_cols

!-------------------------------------------------------------------------------

SUBROUTINE set_pole(this, pole_latitude, pole_longitude)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(INOUT) :: this
REAL, INTENT(IN) :: pole_latitude, pole_longitude
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_POLE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL this%horiz_grid%set_pole(pole_latitude, pole_longitude)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_pole

!-------------------------------------------------------------------------------

SUBROUTINE set_pole_latitude(this, pole_latitude)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(INOUT) :: this
REAL, INTENT(IN) :: pole_latitude
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_POLE_LATITUDE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL this%horiz_grid%set_pole_latitude(pole_latitude)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_pole_latitude

!-------------------------------------------------------------------------------

SUBROUTINE set_pole_longitude(this, pole_longitude)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(INOUT) :: this
REAL, INTENT(IN) :: pole_longitude
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_POLE_LONGITUDE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL this%horiz_grid%set_pole_longitude(pole_longitude)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_pole_longitude

!-------------------------------------------------------------------------------

SUBROUTINE set_regular_horizontal_grid(this, nx, ny, startx, starty, dx, dy)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: nx, ny
REAL, INTENT(IN) :: startx, starty, dx, dy
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_REGULAR_HORIZONTAL_GRID'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL this%horiz_grid%set_regular_horizontal_grid(nx, ny, startx, starty, dx, dy)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_regular_horizontal_grid

!-------------------------------------------------------------------------------

SUBROUTINE set_defined_horizontal_grid(this, nx, ny, list_x, list_y)
IMPLICIT NONE
CLASS(three_dimensional_grid_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: nx, ny
REAL, INTENT(IN) :: list_x(nx), list_y(ny)
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_DEFINED_HORIZONTAL_GRID'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL this%horiz_grid%set_defined_horizontal_grid(nx, ny, list_x, list_y)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_defined_horizontal_grid

!-------------------------------------------------------------------------------

END MODULE three_dimensional_grid_mod
  
  
