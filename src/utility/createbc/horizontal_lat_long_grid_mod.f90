! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE horizontal_lat_long_grid_mod

! DrHook Modules
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

! Description:
!   A class to describe a horizontal lat-long grid. 
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                        &
  ModuleName = 'HORIZONTAL_LAT_LONG_GRID_MOD'

TYPE, PUBLIC :: horizontal_lat_long_grid_type
  PRIVATE
  REAL, ALLOCATABLE, PUBLIC :: latitudes(:), longitudes(:)
  REAL, PRIVATE :: pole_latitude, pole_longitude

  ! Grid code from pp lookup 16 i.e. lat long rotated
  INTEGER, PRIVATE :: grid_code
  ! Type of grid points i.e. u points on c grid
  INTEGER, PRIVATE :: grid_type 
  ! Domain of grid i.e. global or LAM
  INTEGER, PRIVATE :: hemisphere_indicator
  ! Size of halos (if included in grid)
  INTEGER, PRIVATE :: halo_code
  ! Size of halo region
  INTEGER, PRIVATE :: halo_ns = 0
  INTEGER, PRIVATE :: halo_ew = 0

  INTEGER, PRIVATE :: num_rows, num_cols
  CONTAINS
  PROCEDURE, PASS :: compare_grids
  PROCEDURE, PASS :: set_regular_horizontal_grid
  PROCEDURE, PASS :: set_defined_horizontal_grid
  PROCEDURE, PASS :: get_regular_horizontal_grid
  PROCEDURE, PASS :: get_defined_horizontal_grid
  PROCEDURE, PASS :: get_horizontal_coordinates
  PROCEDURE, PASS :: get_grid_type
  PROCEDURE, PASS :: get_latitudes
  PROCEDURE, PASS :: get_longitudes
  PROCEDURE, PASS :: set_latitudes
  PROCEDURE, PASS :: set_longitudes
  PROCEDURE, PASS :: get_num_rows
  PROCEDURE, PASS :: get_halo_ns
  PROCEDURE, PASS :: get_halo_ew
  PROCEDURE, PASS :: get_num_cols
  PROCEDURE, PASS :: get_grid_code
  PROCEDURE, PASS :: get_hemisphere_indicator
  PROCEDURE, PASS :: get_halo_code
  PROCEDURE, PASS :: set_pole
  PROCEDURE, PASS :: get_pole
  PROCEDURE, PASS :: get_pole_latitude
  PROCEDURE, PASS :: set_pole_latitude
  PROCEDURE, PASS :: get_pole_longitude
  PROCEDURE, PASS :: set_pole_longitude
  PROCEDURE, PASS :: set_grid_code_rotation
  PROCEDURE, PASS :: set_grid_type
  PROCEDURE, PASS :: set_num_rows
  PROCEDURE, PASS :: set_grid_code
  PROCEDURE, PASS :: set_hemisphere_indicator
  PROCEDURE, PASS :: set_halo_code
  PROCEDURE, PASS :: set_halo_ns
  PROCEDURE, PASS :: set_halo_ew
  PROCEDURE, PASS :: set_num_cols
  
END TYPE horizontal_lat_long_grid_type

CONTAINS

!-------------------------------------------------------------------------------

SUBROUTINE compare_grids(this, compare_grid)

USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

CLASS(horizontal_lat_long_grid_type), INTENT(IN) :: this
CLASS(horizontal_lat_long_grid_type), INTENT(IN) :: compare_grid

! Error reporting variables
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER :: routinename = 'COMPARE_GRIDS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (this%num_rows /= compare_grid%num_rows) THEN
  WRITE(cmessage, '(A,I0,A,I0)') "Number of rows not consistent "//                &
       "between the two grid objects. Grid 1: ",this%num_rows,                     &
       " Grid 2: ",compare_grid%num_rows
  icode = 10
  CALL ereport(routinename, icode, cmessage)
END IF

IF (this%num_cols /= compare_grid%num_cols) THEN
  WRITE(cmessage, '(A,I0,A,I0)') "Number of columns not consistent "//             &
       "between the two grid objects. Grid 1: ",this%num_cols,                     &
       " Grid 2: ",compare_grid%num_cols
  icode = 20
  CALL ereport(routinename, icode, cmessage)
END IF

IF (ABS(this%pole_latitude - compare_grid%pole_latitude)                           &
     > EPSILON(1.0)) THEN
  WRITE(cmessage, '(A,F20.10,A,F20.10)') "Pole latitude not consistent "//         &
       "between the two grid objects. Grid 1: ",this%pole_latitude,                &
       " Grid 2: ",compare_grid%pole_latitude
  icode = 30
  CALL ereport(routinename, icode, cmessage)
END IF

IF (ABS(this%pole_longitude - compare_grid%pole_longitude)                         &
     > EPSILON(1.0)) THEN
  WRITE(cmessage, '(A,F20.10,A,F20.10)') "Pole longitude not consistent "//        &
       "between the two grid objects. Grid 1: ",this%pole_longitude,               &
       " Grid 2: ",compare_grid%pole_longitude
  icode = 40
  CALL ereport(routinename, icode, cmessage)
END IF

IF (ABS(this%latitudes(1) - compare_grid%latitudes(1))                             &
     > EPSILON(1.0)) THEN
  WRITE(cmessage, '(A,F20.10,A,F20.10)') "First latitude point not consistent "//  &
       "between the two grid objects. Grid 1: ",this%latitudes(1),                 &
       " Grid 2: ",compare_grid%latitudes(1)
  icode = 50
  CALL ereport(routinename, icode, cmessage)
END IF

IF (ABS(this%longitudes(1) - compare_grid%longitudes(1))                           &
     > EPSILON(1.0)) THEN
  WRITE(cmessage, '(A,F20.10,A,F20.10)') "First longitude point not consistent "// &
       "between the two grid objects. Grid 1: ",this%longitudes(1),                &
       " Grid 2: ",compare_grid%longitudes(1)
  icode = 60
  CALL ereport(routinename, icode, cmessage)
END IF

IF (this%grid_code /= compare_grid%grid_code) THEN
  WRITE(cmessage, '(A,I0,A,I0)') "Grid code not consistent "//                     &
       "between the two grid objects. Grid 1: ",this%grid_code,                    &
       " Grid 2: ",compare_grid%grid_code
  icode = 70
  CALL ereport(routinename, icode, cmessage)
END IF

IF (this%grid_type /= compare_grid%grid_type) THEN
  WRITE(cmessage, '(A,I0,A,I0)') "Grid type not consistent "//                     &
       "between the two grid objects. Grid 1: ",this%grid_type,                    &
       " Grid 2: ",compare_grid%grid_type
  icode = 80
  CALL ereport(routinename, icode, cmessage)
END IF

IF (this%hemisphere_indicator /= compare_grid%hemisphere_indicator) THEN
  WRITE(cmessage, '(A,I0,A,I0)') "hemisphere_indicator not consistent "//          &
       "between the two grid objects. Grid 1: ",this%hemisphere_indicator,         &
       " Grid 2: ",compare_grid%hemisphere_indicator
  icode = 90
  CALL ereport(routinename, icode, cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE compare_grids

!-------------------------------------------------------------------------------

SUBROUTINE get_horizontal_coordinates(this, latitudes, longitudes)
IMPLICIT NONE
CLASS(horizontal_lat_long_grid_type), INTENT(IN) :: this
REAL, ALLOCATABLE, INTENT(OUT) :: latitudes(:), longitudes(:)
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_HORIZONTAL_COORDINATES'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
latitudes = this%get_latitudes()
longitudes =  this%get_longitudes()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE get_horizontal_coordinates

!-------------------------------------------------------------------------------

SUBROUTINE set_longitudes(this, longitudes)
IMPLICIT NONE 
CLASS(horizontal_lat_long_grid_type), INTENT(INOUT) :: this
REAL, INTENT(IN)    :: longitudes(:)
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_LONGITUDES'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (ALLOCATED(this%longitudes)) DEALLOCATE(this%longitudes)
this%longitudes = longitudes

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_longitudes

!-------------------------------------------------------------------------------

SUBROUTINE set_latitudes(this, latitudes)
IMPLICIT NONE 
CLASS(horizontal_lat_long_grid_type), INTENT(INOUT) :: this
REAL, INTENT(IN)    :: latitudes(:)
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_LATITUDES'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (ALLOCATED(this%latitudes)) DEALLOCATE(this%latitudes)
this%latitudes = latitudes

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_latitudes

!-------------------------------------------------------------------------------

FUNCTION get_longitudes(this) RESULT(longitudes)
IMPLICIT NONE 
CLASS(horizontal_lat_long_grid_type), INTENT(IN) :: this
REAL, ALLOCATABLE :: longitudes(:)
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_LONGITUDES'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE(longitudes(SIZE(this%longitudes)))
longitudes(:) = this%longitudes(:)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION get_longitudes

!-------------------------------------------------------------------------------

FUNCTION get_latitudes(this) RESULT(latitudes)
IMPLICIT NONE 
CLASS(horizontal_lat_long_grid_type), INTENT(IN) :: this
REAL, ALLOCATABLE :: latitudes(:)
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_LATITUDES'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE(latitudes(SIZE(this%latitudes)))
latitudes(:) = this%latitudes(:)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION get_latitudes

!-------------------------------------------------------------------------------

SUBROUTINE set_regular_horizontal_grid(this, nx, ny, startx, starty, dx, dy)
! Set up a regular grid given the number of points, starting lat/long and the
! spacing
IMPLICIT NONE
CLASS(horizontal_lat_long_grid_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: nx, ny
REAL, INTENT(IN) :: startx, starty, dx, dy
INTEGER :: col, row
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_REGULAR_HORIZONTAL_GRID'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
this%num_rows = ny
this%num_cols = nx
ALLOCATE(this%latitudes(this%num_rows))
ALLOCATE(this%longitudes(this%num_cols))
DO row = 1, ny
  this%latitudes(row) = starty + (row - 1) * dy
END DO
DO col = 1, nx
  this%longitudes(col) = startx + (col - 1) * dx
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_regular_horizontal_grid

!-------------------------------------------------------------------------------

SUBROUTINE get_regular_horizontal_grid(this, nx, ny, startx, starty, dx, dy)
! Return the number of points, starting position, and spacing for lat and long
! for a given grid object. 
IMPLICIT NONE
CLASS(horizontal_lat_long_grid_type), INTENT(INOUT) :: this
INTEGER, INTENT(OUT) :: nx, ny
REAL, INTENT(OUT) :: startx, starty, dx, dy
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_REGULAR_HORIZONTAL_GRID'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
dx = this%longitudes(2) - this%longitudes(1)
dy = this%latitudes(2) - this%latitudes(1)
nx = this%get_num_cols()
ny = this%get_num_rows()
startx = this%longitudes(1)
starty = this%latitudes(1)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE get_regular_horizontal_grid

!-------------------------------------------------------------------------------

SUBROUTINE set_defined_horizontal_grid(this, nx, ny, list_x, list_y)
! Set up a grid exactly given the number of points in lat/long together with
! an array containing the actual lat/long values.
IMPLICIT NONE
CLASS(horizontal_lat_long_grid_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: nx, ny
REAL, INTENT(IN) :: list_x(nx), list_y(ny)
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_DEFINED_HORIZONTAL_GRID'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
this%num_rows = ny
this%num_cols = nx
ALLOCATE(this%latitudes(this%num_rows))
ALLOCATE(this%longitudes(this%num_cols))
this%longitudes = list_x
this%latitudes = list_y

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_defined_horizontal_grid

!-------------------------------------------------------------------------------

SUBROUTINE get_defined_horizontal_grid(this, nx, ny, list_x, list_y)
! Return the number of lat/long points and arrays containing their lat/longs.
IMPLICIT NONE
CLASS(horizontal_lat_long_grid_type), INTENT(IN) :: this
INTEGER, INTENT(OUT) :: nx, ny
REAL, INTENT(OUT) :: list_x(this%num_cols), list_y(this%num_rows)
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_DEFINED_HORIZONTAL_GRID'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
nx = this%num_cols
ny = this%num_rows
list_x = this%longitudes
list_y = this%latitudes

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE get_defined_horizontal_grid

!-------------------------------------------------------------------------------

INTEGER FUNCTION get_grid_code(this)
IMPLICIT NONE
CLASS(horizontal_lat_long_grid_type), INTENT(IN) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_GRID_CODE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
get_grid_code = this%grid_code

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION get_grid_code

!-------------------------------------------------------------------------------

INTEGER FUNCTION get_hemisphere_indicator(this)
IMPLICIT NONE
CLASS(horizontal_lat_long_grid_type), INTENT(IN) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_HEMISPHERE_INDICATOR'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
get_hemisphere_indicator = this%hemisphere_indicator

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION get_hemisphere_indicator

!-------------------------------------------------------------------------------

INTEGER FUNCTION get_halo_code(this)
IMPLICIT NONE
CLASS(horizontal_lat_long_grid_type), INTENT(IN) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_HALO_CODE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
get_halo_code = this%halo_code

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION get_halo_code

!-------------------------------------------------------------------------------

INTEGER FUNCTION get_halo_ns(this)
IMPLICIT NONE
CLASS(horizontal_lat_long_grid_type), INTENT(IN) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_HALO_NS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
get_halo_ns = this%halo_ns

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION get_halo_ns

!-------------------------------------------------------------------------------

INTEGER FUNCTION get_halo_ew(this)
IMPLICIT NONE
CLASS(horizontal_lat_long_grid_type), INTENT(IN) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_HALO_EW'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
get_halo_ew = this%halo_ew

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION get_halo_ew

!-------------------------------------------------------------------------------

INTEGER FUNCTION get_num_rows(this)
IMPLICIT NONE
CLASS(horizontal_lat_long_grid_type), INTENT(IN) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_NUM_ROWS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
get_num_rows = this%num_rows

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION get_num_rows

!-------------------------------------------------------------------------------


INTEGER FUNCTION get_num_cols(this)
IMPLICIT NONE
CLASS(horizontal_lat_long_grid_type), INTENT(IN) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_NUM_COLS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
get_num_cols = this%num_cols

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION get_num_cols

!-------------------------------------------------------------------------------

SUBROUTINE set_pole(this, pole_latitude, pole_longitude)
IMPLICIT NONE
CLASS(horizontal_lat_long_grid_type), INTENT(INOUT) :: this
REAL, INTENT(IN) :: pole_latitude, pole_longitude
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_POLE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
this%pole_latitude = pole_latitude
this%pole_longitude = pole_longitude
CALL this%set_grid_code_rotation()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_pole

!-------------------------------------------------------------------------------

SUBROUTINE get_pole(this, pole_latitude, pole_longitude)
IMPLICIT NONE
CLASS(horizontal_lat_long_grid_type), INTENT(IN) :: this
REAL, INTENT(OUT) :: pole_latitude, pole_longitude
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_POLE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
pole_latitude = this%pole_latitude
pole_longitude = this%pole_longitude

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE get_pole

!-------------------------------------------------------------------------------

SUBROUTINE set_pole_latitude(this, pole_latitude)
IMPLICIT NONE
CLASS(horizontal_lat_long_grid_type), INTENT(INOUT) :: this
REAL, INTENT(IN) :: pole_latitude
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_POLE_LATITUDE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
this%pole_latitude = pole_latitude
CALL this%set_grid_code_rotation()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_pole_latitude

!-------------------------------------------------------------------------------

FUNCTION get_pole_latitude(this)
IMPLICIT NONE
CLASS(horizontal_lat_long_grid_type), INTENT(IN) :: this
REAL ::  get_pole_latitude
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_POLE_LATITUDE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
get_pole_latitude = this%pole_latitude

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION get_pole_latitude

!-------------------------------------------------------------------------------

SUBROUTINE set_pole_longitude(this, pole_longitude)
IMPLICIT NONE
CLASS(horizontal_lat_long_grid_type), INTENT(INOUT) :: this
REAL, INTENT(IN) :: pole_longitude
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_POLE_LONGITUDE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
this%pole_longitude = pole_longitude
CALL this%set_grid_code_rotation()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_pole_longitude

!-------------------------------------------------------------------------------

FUNCTION get_pole_longitude(this)
IMPLICIT NONE
CLASS(horizontal_lat_long_grid_type), INTENT(IN) :: this
REAL :: get_pole_longitude
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_POLE_LONGITUDE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
get_pole_longitude = this%pole_longitude

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION get_pole_longitude

!-------------------------------------------------------------------------------

INTEGER FUNCTION get_grid_type(this)
IMPLICIT NONE
CLASS(horizontal_lat_long_grid_type), INTENT(INOUT) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_GRID_TYPE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
get_grid_type = this%grid_type

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION get_grid_type

!-------------------------------------------------------------------------------

SUBROUTINE set_grid_type(this, grid_type)
IMPLICIT NONE
CLASS(horizontal_lat_long_grid_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: grid_type
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_GRID_TYPE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
this%grid_type = grid_type

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_grid_type

!-------------------------------------------------------------------------------

SUBROUTINE set_grid_code_rotation(this)
IMPLICIT NONE
CLASS(horizontal_lat_long_grid_type), INTENT(INOUT) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_GRID_CODE_ROTATION'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! Check if the poles of the grid are on a non-standard rotation
IF (this%pole_latitude /= 90.0 .OR. &
    this%pole_longitude /= 0.0) THEN
  this%grid_code = 101
ELSE
  this%grid_code = 1
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_grid_code_rotation

!-------------------------------------------------------------------------------

SUBROUTINE set_halo_ew(this, halo_ew)
IMPLICIT NONE
CLASS(horizontal_lat_long_grid_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: halo_ew
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_HALO_EW'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
this%halo_ew = halo_ew

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_halo_ew

!-------------------------------------------------------------------------------

SUBROUTINE set_halo_ns(this, halo_ns)
IMPLICIT NONE
CLASS(horizontal_lat_long_grid_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: halo_ns
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_HALO_NS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
this%halo_ns = halo_ns

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_halo_ns

!-------------------------------------------------------------------------------

SUBROUTINE set_grid_code(this, grid_code)
IMPLICIT NONE
CLASS(horizontal_lat_long_grid_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: grid_code
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_GRID_CODE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
this%grid_code = grid_code

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_grid_code

!-------------------------------------------------------------------------------

SUBROUTINE set_hemisphere_indicator(this, hemisphere_indicator)
IMPLICIT NONE
CLASS(horizontal_lat_long_grid_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: hemisphere_indicator
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_HEMISPHERE_INDICATOR'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
this%hemisphere_indicator = hemisphere_indicator

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_hemisphere_indicator

!-------------------------------------------------------------------------------

SUBROUTINE set_halo_code(this, halo_code)
IMPLICIT NONE
CLASS(horizontal_lat_long_grid_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: halo_code
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_HALO_CODE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
this%halo_code = halo_code

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_halo_code

!-------------------------------------------------------------------------------

SUBROUTINE set_num_rows(this, num_rows)
IMPLICIT NONE
CLASS(horizontal_lat_long_grid_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: num_rows
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_NUM_ROWS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
this%num_rows = num_rows

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_num_rows

!-------------------------------------------------------------------------------


SUBROUTINE set_num_cols(this, num_cols)
IMPLICIT NONE
CLASS(horizontal_lat_long_grid_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: num_cols
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_NUM_COLS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
this%num_cols = num_cols

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_num_cols

!-------------------------------------------------------------------------------

END MODULE horizontal_lat_long_grid_mod
