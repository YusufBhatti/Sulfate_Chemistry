! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE field_mod

USE three_dimensional_grid_mod, ONLY: three_dimensional_grid_type
! DrHook Modules
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE
PRIVATE

! Description:
!   A class to describe a field. 
!   This contains the data itself together with any metadata to describe it.
!   It specifically does NOT contain any information about how the data is to
!   be stored on disk (see the data_location class for that).
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.

PUBLIC :: ASSIGNMENT(=), field_type

CHARACTER(LEN=*), PARAMETER :: ModuleName = 'FIELD_MOD'

TYPE :: field_type

  PRIVATE
  ! Data and validity times (year, month, date, hour, min, sec)
  INTEGER, PUBLIC :: validity_time(6)    ! Fieldsfile lookup items 1-6
  INTEGER, PUBLIC :: data_time(6)        ! Fieldsfile lookup items 7-12

  ! Calendar type
  INTEGER, PUBLIC :: calendar            ! fieldsfile fixed header item 8

  ! Time indicator
  INTEGER, PUBLIC :: time_indicator      ! Fieldsfile lookup item 13

  ! The forecast period
  INTEGER, PUBLIC :: fctime              ! Fieldsfile lookup item 14

  ! The packing method
  INTEGER, PUBLIC :: packing_method      ! Fieldsfile lookup item 21

  ! Two variables to contain a number describing the contents scientifically.
  ! For the UM, quantity_ident is the STASH code - lookup item 42
  ! This is expected to be the stashcode in CreateBC, so other file formats will
  ! need to convert as appropriate when reading the header.
  INTEGER, PUBLIC :: quantity_ident

  ! Store processing has been done on this field
  ! For the UM, this is LBPROC, lookup item 25
  INTEGER, PUBLIC :: processing_code   

  ! A user-defined reference number, fieldsfile lookup item 28
  INTEGER, PUBLIC :: user_reference

  ! Variables to store information about the source of this field
  ! These correspond to fieldsfile fixed header items 6 and 7
  INTEGER, PUBLIC :: run_ident
  INTEGER, PUBLIC :: expt_ident

  ! Type of data: 1 = real, 2 = integer, 3 = logical
  ! Corresponds to fieldsfile lookup item 39
  INTEGER, PUBLIC :: data_type 

  ! Meaning of vertical dimension
  ! For the UM, this is 0 for normal levels, or 1 for pseudo-level or a single
  ! level field (8888)
  INTEGER, PUBLIC :: meaning_of_k_dimension

  ! Source model identifier
  INTEGER, PUBLIC :: source_model_identifier ! Fieldsfile lookup item 38

  ! Datum constant, subtracted from each field value
  ! Fieldsfile lookup item 50
  REAL, PUBLIC :: datum_constant

  ! Accuracy of packing
  ! Fieldsfile lookup item 51
  REAL, PUBLIC :: packing_accuracy
  
  ! Missing data indicator - fieldsfile lookup item 63
  REAL, PUBLIC :: mdi  

  ! The data itself
  REAL, ALLOCATABLE, PUBLIC :: rdata(:,:,:)
  INTEGER, ALLOCATABLE, PUBLIC :: idata(:,:,:)
  LOGICAL, ALLOCATABLE, PUBLIC :: ldata(:,:,:)
  ! LBC data 
  REAL, ALLOCATABLE, PUBLIC :: lbc_rdata(:,:) ! lbc_points, levels

  ! The grid the data is stored on
  TYPE(three_dimensional_grid_type), PUBLIC :: grid

  CONTAINS
  PROCEDURE, PASS :: unload_data
  PROCEDURE, PASS :: get_grid_dimensions
  PROCEDURE, PASS :: get_horiz_grid_type
  PROCEDURE, PASS :: get_horizontal_coordinates
  PROCEDURE, PASS :: get_vertical_levels
  PROCEDURE, PASS :: get_vertical_level_values
  PROCEDURE, PASS :: get_first_level
  PROCEDURE, PASS :: get_last_level
  PROCEDURE, PASS :: get_latitudes
  PROCEDURE, PASS :: get_longitudes
  PROCEDURE, PASS :: set_latitudes
  PROCEDURE, PASS :: set_longitudes
  PROCEDURE, PASS :: get_num_levels
  PROCEDURE, PASS :: get_num_model_levels
  PROCEDURE, PASS :: get_grid_code
  PROCEDURE, PASS :: get_hemisphere_indicator
  PROCEDURE, PASS :: get_halo_code
  PROCEDURE, PASS :: get_halo_ns
  PROCEDURE, PASS :: get_halo_ew
  PROCEDURE, PASS :: get_num_rows
  PROCEDURE, PASS :: get_num_cols
  PROCEDURE, PASS :: maximum
  PROCEDURE, PASS :: minimum
  PROCEDURE, PASS :: surface_mean
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
  PROCEDURE, PASS :: get_pole
  PROCEDURE, PASS :: set_pole
  PROCEDURE, PASS :: get_pole_latitude
  PROCEDURE, PASS :: set_pole_latitude
  PROCEDURE, PASS :: get_pole_longitude
  PROCEDURE, PASS :: set_pole_longitude
  PROCEDURE, PASS :: has_stashcode
  PROCEDURE, PASS :: return_data_time_string
  PROCEDURE, PASS :: return_validity_time_string
  PROCEDURE, PASS :: update_um_vn
  PROCEDURE, PASS :: calc_header_release
  PROCEDURE, PASS :: copy_metadata_from_field
  PROCEDURE, PASS :: modify_packing_method

END TYPE field_type

INTERFACE ASSIGNMENT(=)
  MODULE PROCEDURE assign_field_type
END INTERFACE ASSIGNMENT(=)

CONTAINS

!-------------------------------------------------------------------------------

SUBROUTINE unload_data(this)
IMPLICIT NONE
CLASS(field_type), INTENT(INOUT) :: this

CHARACTER(LEN=*), PARAMETER :: routinename = 'UNLOAD_DATA'

! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (ALLOCATED(this%rdata)) THEN
  DEALLOCATE(this%rdata)
END IF
IF (ALLOCATED(this%lbc_rdata)) THEN
  DEALLOCATE(this%lbc_rdata)
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE unload_data

!-------------------------------------------------------------------------------

SUBROUTINE get_grid_dimensions(this, nx, ny, nz)
IMPLICIT NONE
CLASS(field_type), INTENT(INOUT) :: this
INTEGER, INTENT(OUT) :: nx, ny, nz
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_GRID_DIMENSIONS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL this%grid%get_grid_dimensions(nx, ny, nz)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE get_grid_dimensions

!-------------------------------------------------------------------------------

INTEGER FUNCTION get_horiz_grid_type(this)
IMPLICIT NONE
CLASS(field_type), INTENT(INOUT) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_HORIZ_GRID_TYPE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
get_horiz_grid_type = this%grid%get_horiz_grid_type()
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION get_horiz_grid_type

!-------------------------------------------------------------------------------

SUBROUTINE get_horizontal_coordinates(this, latitudes, longitudes)
IMPLICIT NONE
CLASS(field_type), INTENT(INOUT) :: this
REAL, ALLOCATABLE, INTENT(INOUT) :: latitudes(:), longitudes(:)
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_HORIZONTAL_COORDINATES'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL this%grid%get_horizontal_coordinates(latitudes, longitudes)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE get_horizontal_coordinates

!-------------------------------------------------------------------------------

SUBROUTINE get_vertical_levels(this, levels)
IMPLICIT NONE
CLASS(field_type), INTENT(INOUT) :: this
INTEGER, INTENT(INOUT) :: levels(:)
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_VERTICAL_LEVELS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL this%grid%get_vertical_levels(levels)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE get_vertical_levels

!-------------------------------------------------------------------------------

FUNCTION get_latitudes(this) RESULT(latitudes)
IMPLICIT NONE
CLASS(field_type), INTENT(IN) :: this
REAL, ALLOCATABLE :: local_latitudes(:), latitudes(:)
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_LATITUDES'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
local_latitudes = this%grid%get_latitudes()
CALL MOVE_ALLOC(local_latitudes, latitudes)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION get_latitudes

!-------------------------------------------------------------------------------

FUNCTION get_longitudes(this) RESULT(longitudes)
IMPLICIT NONE
CLASS(field_type), INTENT(IN) :: this
REAL, ALLOCATABLE :: local_longitudes(:), longitudes(:)
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_LONGITUDES'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
local_longitudes = this%grid%get_longitudes()
CALL MOVE_ALLOC(local_longitudes, longitudes)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION get_longitudes

!-------------------------------------------------------------------------------

SUBROUTINE set_longitudes(this, longitudes)
IMPLICIT NONE 
CLASS(field_type), INTENT(INOUT) :: this
REAL, INTENT(IN)    :: longitudes(:)
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_LONGITUDES'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL this%grid%set_longitudes(longitudes)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_longitudes

!-------------------------------------------------------------------------------

SUBROUTINE set_latitudes(this, latitudes)
IMPLICIT NONE 
CLASS(field_type), INTENT(INOUT) :: this
REAL, INTENT(IN)    :: latitudes(:)
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_LATITUDES'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL this%grid%set_latitudes(latitudes)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_latitudes

!-------------------------------------------------------------------------------

SUBROUTINE get_pole(this, pole_latitude, pole_longitude)
IMPLICIT NONE
CLASS(field_type), INTENT(IN) :: this
REAL, INTENT(OUT) :: pole_latitude, pole_longitude
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_POLE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL this%grid%get_pole(pole_latitude, pole_longitude)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE get_pole

!-------------------------------------------------------------------------------

FUNCTION get_pole_latitude(this)
IMPLICIT NONE
CLASS(field_type), INTENT(IN) :: this
REAL :: get_pole_latitude
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_POLE_LATITUDE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
get_pole_latitude = this%grid%get_pole_latitude()
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION get_pole_latitude

!-------------------------------------------------------------------------------

FUNCTION get_pole_longitude(this)
IMPLICIT NONE
CLASS(field_type), INTENT(IN) :: this
REAL :: get_pole_longitude
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_POLE_LONGITUDE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
get_pole_longitude = this%grid%get_pole_longitude()
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION get_pole_longitude

!-------------------------------------------------------------------------------

SUBROUTINE set_pole(this, pole_latitude, pole_longitude)
IMPLICIT NONE
CLASS(field_type), INTENT(INOUT) :: this
REAL, INTENT(IN) :: pole_latitude, pole_longitude
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_POLE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL this%grid%set_pole(pole_latitude, pole_longitude)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_pole

!-------------------------------------------------------------------------------

SUBROUTINE set_pole_latitude(this, pole_latitude)
IMPLICIT NONE
CLASS(field_type), INTENT(INOUT) :: this
REAL, INTENT(IN) :: pole_latitude
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_POLE_LATITUDE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL this%grid%set_pole_latitude(pole_latitude)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_pole_latitude

!-------------------------------------------------------------------------------

SUBROUTINE set_pole_longitude(this, pole_longitude)
IMPLICIT NONE
CLASS(field_type), INTENT(INOUT) :: this
REAL, INTENT(IN) :: pole_longitude
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_POLE_LONGITUDE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL this%grid%set_pole_longitude(pole_longitude)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_pole_longitude

!-------------------------------------------------------------------------------

SUBROUTINE get_vertical_level_values(this, level_values)
IMPLICIT NONE
CLASS(field_type), INTENT(INOUT) :: this
REAL, INTENT(INOUT) :: level_values(:)
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_VERTICAL_LEVEL_VALUES'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL this%grid%get_vertical_level_values(level_values)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE get_vertical_level_values

!-------------------------------------------------------------------------------

INTEGER FUNCTION get_first_level(this)
IMPLICIT NONE
CLASS(field_type), INTENT(IN) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_FIRST_LEVEL'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
get_first_level = this%grid%vert_grid%first_level
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION  get_first_level

!-------------------------------------------------------------------------------

INTEGER FUNCTION get_last_level(this)
IMPLICIT NONE
CLASS(field_type), INTENT(IN) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_LAST_LEVEL'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
get_last_level = this%grid%vert_grid%last_level
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION get_last_level

!-------------------------------------------------------------------------------

INTEGER FUNCTION get_num_levels(this)
IMPLICIT NONE
CLASS(field_type), INTENT(IN) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_NUM_LEVELS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
get_num_levels = this%grid%get_num_levels()
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION get_num_levels

!-------------------------------------------------------------------------------

INTEGER FUNCTION get_num_model_levels(this)
IMPLICIT NONE
CLASS(field_type), INTENT(IN) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_NUM_MODEL_LEVELS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
get_num_model_levels = this%grid%get_num_model_levels()
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION get_num_model_levels

!-------------------------------------------------------------------------------

SUBROUTINE set_grid_code(this, grid_code)
IMPLICIT NONE
CLASS(field_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: grid_code
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_GRID_CODE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL this%grid%set_grid_code(grid_code)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_grid_code

!-------------------------------------------------------------------------------

SUBROUTINE set_hemisphere_indicator(this, hemisphere_indicator)
IMPLICIT NONE
CLASS(field_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: hemisphere_indicator
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_HEMISPHERE_INDICATOR'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL this%grid%set_hemisphere_indicator(hemisphere_indicator)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_hemisphere_indicator

!-------------------------------------------------------------------------------

SUBROUTINE set_halo_code(this, halo_code)
IMPLICIT NONE
CLASS(field_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: halo_code
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_HALO_CODE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL this%grid%set_halo_code(halo_code)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_halo_code

!-------------------------------------------------------------------------------

SUBROUTINE set_halo_ns(this, halo_ns)
IMPLICIT NONE
CLASS(field_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: halo_ns
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_HALO_NS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL this%grid%set_halo_ns(halo_ns)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_halo_ns

!-------------------------------------------------------------------------------

SUBROUTINE set_halo_ew(this, halo_ew)
IMPLICIT NONE
CLASS(field_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: halo_ew
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_HALO_EW'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL this%grid%set_halo_ew(halo_ew)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_halo_ew

!-------------------------------------------------------------------------------

INTEGER FUNCTION get_grid_code(this)
IMPLICIT NONE
CLASS(field_type), INTENT(IN) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_GRID_CODE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
get_grid_code = this%grid%get_grid_code()
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION get_grid_code

!-------------------------------------------------------------------------------

INTEGER FUNCTION get_hemisphere_indicator(this)
IMPLICIT NONE
CLASS(field_type), INTENT(IN) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_HEMISPHERE_INDICATOR'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
get_hemisphere_indicator = this%grid%get_hemisphere_indicator()
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION get_hemisphere_indicator

!-------------------------------------------------------------------------------

INTEGER FUNCTION get_halo_code(this)
IMPLICIT NONE
CLASS(field_type), INTENT(IN) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_HALO_CODE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
get_halo_code = this%grid%get_halo_code()
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION get_halo_code

!-------------------------------------------------------------------------------

INTEGER FUNCTION get_halo_ns(this)
IMPLICIT NONE
CLASS(field_type), INTENT(IN) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_HALO_NS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
get_halo_ns = this%grid%get_halo_ns()
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION get_halo_ns

!-------------------------------------------------------------------------------

INTEGER FUNCTION get_halo_ew(this)
IMPLICIT NONE
CLASS(field_type), INTENT(IN) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_HALO_EW'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
get_halo_ew = this%grid%get_halo_ew()
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION get_halo_ew

!-------------------------------------------------------------------------------

INTEGER FUNCTION get_num_rows(this)
IMPLICIT NONE
CLASS(field_type), INTENT(IN) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_NUM_ROWS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
get_num_rows = this%grid%get_num_rows()
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION get_num_rows

!-------------------------------------------------------------------------------

INTEGER FUNCTION get_num_cols(this)
IMPLICIT NONE
CLASS(field_type), INTENT(IN) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_NUM_COLS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
get_num_cols = this%grid%get_num_cols()
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION get_num_cols

!-------------------------------------------------------------------------------

LOGICAL FUNCTION has_stashcode(this, stashcode)
CLASS(field_type), INTENT(IN) :: this
INTEGER :: stashcode
CHARACTER(LEN=*), PARAMETER :: routinename = 'HAS_STASHCODE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (this%quantity_ident == stashcode) THEN
  has_stashcode = .TRUE.
ELSE
  has_stashcode = .FALSE.
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION has_stashcode

!-------------------------------------------------------------------------------

REAL FUNCTION maximum(this) RESULT(max_value)
IMPLICIT NONE
CLASS(field_type), INTENT(IN) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'MAXIMUM'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (ALLOCATED(this%rdata)) THEN
  max_value = MAXVAL(this%rdata)
ELSE
  max_value = this%mdi
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION maximum

!-------------------------------------------------------------------------------

REAL FUNCTION minimum(this) RESULT(min_value)
IMPLICIT NONE
CLASS(field_type), INTENT(IN) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'MINIMUM'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (ALLOCATED(this%rdata)) THEN
  min_value = MINVAL(this%rdata)
ELSE
  min_value = this%mdi
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION minimum

!-------------------------------------------------------------------------------

SUBROUTINE set_horiz_grid_type(this, grid_type)
IMPLICIT NONE
CLASS(field_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: grid_type
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_HORIZ_GRID_TYPE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL this%grid%set_horiz_grid_type(grid_type)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_horiz_grid_type

!-------------------------------------------------------------------------------

SUBROUTINE set_num_levels(this, num_levels)
IMPLICIT NONE
CLASS(field_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: num_levels
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_NUM_LEVELS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL this%grid%set_num_levels(num_levels)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_num_levels

!-------------------------------------------------------------------------------

SUBROUTINE set_num_model_levels(this, num_model_levels)
IMPLICIT NONE
CLASS(field_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: num_model_levels
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_NUM_MODEL_LEVELS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL this%grid%set_num_model_levels(num_model_levels)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_num_model_levels

!-------------------------------------------------------------------------------

SUBROUTINE set_num_rows(this, num_rows)
IMPLICIT NONE
CLASS(field_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: num_rows
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_NUM_ROWS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL this%grid%set_num_rows(num_rows)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_num_rows

!-------------------------------------------------------------------------------

SUBROUTINE set_num_cols(this, num_cols)
IMPLICIT NONE
CLASS(field_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: num_cols
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_NUM_COLS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL this%grid%set_num_cols(num_cols)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_num_cols

!-------------------------------------------------------------------------------

REAL FUNCTION surface_mean(this)
IMPLICIT NONE
INTEGER :: i,j, nx,ny,nz, num
REAL :: total
CLASS(field_type) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'SURFACE_MEAN'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL this%get_grid_dimensions(nx, ny, nz)
  
total = 0.0
num = 0
IF (ALLOCATED(this%rdata)) THEN
  DO j = 1, ny
    DO i = 1, nx
      IF (this%rdata(i,j,1) /= this%mdi) THEN
        total = total + this%rdata(i,j,1)
        num = num + 1
      END IF
    END DO
  END DO
  surface_mean = total / REAL(num)
ELSE
  surface_mean = this%mdi
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION surface_mean

!-------------------------------------------------------------------------------

CHARACTER(LEN=19) FUNCTION return_data_time_string(this)
IMPLICIT NONE
CLASS(field_type), INTENT(IN) :: this
CHARACTER(LEN=19) :: dstring
CHARACTER(LEN=*), PARAMETER :: routinename = 'RETURN_DATA_TIME_STRING'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
WRITE(dstring,'(I4,A,I2.2,A,I2,A,I2.2,A,I2.2,A,I2.2)')                   &
    this%data_time(1),'-',this%data_time(2),'-', this%data_time(3), ' ', &
    this%data_time(4),':',this%data_time(5),'/', this%data_time(6)
  
return_data_time_string = dstring
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION return_data_time_string

!-------------------------------------------------------------------------------

CHARACTER(LEN=19) FUNCTION return_validity_time_string(this)
IMPLICIT NONE
CLASS(field_type), INTENT(IN) :: this
CHARACTER(LEN=19) :: dstring
CHARACTER(LEN=*), PARAMETER :: routinename = 'RETURN_VALIDITY_TIME_STRING'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
WRITE(dstring,'(I4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)')                    &
    this%validity_time(1), '-', this%validity_time(2), '-',                 &
    this%validity_time(3), ' ', this%validity_time(4), ':',                 &
    this%validity_time(5), '/', this%validity_time(6)
  
return_validity_time_string = dstring
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION return_validity_time_string

!-------------------------------------------------------------------------------

SUBROUTINE update_um_vn(this)
USE um_version_mod, ONLY: um_version_int
IMPLICIT NONE
CLASS(field_type), INTENT(INOUT) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'UPDATE_UM_VN'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Update the model identifier based on the code version being run
! 1111 is the code for unified model, see UMDP F3

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
this%source_model_identifier = (um_version_int * 10000) + 1111
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE update_um_vn

!-------------------------------------------------------------------------------

INTEGER FUNCTION calc_header_release(this)
IMPLICIT NONE
CLASS(field_type), INTENT(INOUT) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'CALC_HEADER_RELEASE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! UM version number is less than 8.1 then LBREL (UMDP F3) is
! 2. If it is greater then LBREL is 3.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (this%source_model_identifier < 8010000) THEN
  calc_header_release = 2
ELSE 
  calc_header_release = 3
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION calc_header_release

!-------------------------------------------------------------------------------

SUBROUTINE copy_metadata_from_field(this, source_field)
IMPLICIT NONE
CLASS(field_type), INTENT(INOUT) :: this
TYPE(field_type), INTENT(IN) :: source_field
CHARACTER(LEN=*), PARAMETER :: routinename = 'COPY_METADATA_FROM_FIELD'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Copies field metadata from another field object

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
this%validity_time            = source_field%validity_time            
this%data_time                = source_field%data_time                
this%calendar                 = source_field%calendar                 
this%time_indicator           = source_field%time_indicator           
this%fctime                   = source_field%fctime                   
this%packing_method           = source_field%packing_method           
this%quantity_ident           = source_field%quantity_ident           
this%processing_code          = source_field%processing_code          
this%user_reference           = source_field%user_reference           
this%run_ident                = source_field%run_ident                
this%expt_ident               = source_field%expt_ident               
this%data_type                = source_field%data_type                
this%meaning_of_k_dimension   = source_field%meaning_of_k_dimension   
this%source_model_identifier  = source_field%source_model_identifier  
this%datum_constant           = source_field%datum_constant           
this%packing_accuracy         = source_field%packing_accuracy         
this%mdi                      = source_field%mdi        
this%grid                     = source_field%grid
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE copy_metadata_from_field

!-------------------------------------------------------------------------------
SUBROUTINE modify_packing_method(this, new_packing_method, new_accuracy)
IMPLICIT NONE
CLASS(field_type), INTENT(INOUT) :: this
! Modify the packing method of a field
INTEGER, INTENT(IN) :: new_packing_method
REAL,    INTENT(IN) :: new_accuracy
CHARACTER(LEN=*), PARAMETER :: routinename = 'MODIFY_PACKING_METHOD'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
this%packing_method   = new_packing_method
this%packing_accuracy = new_accuracy
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE modify_packing_method
!-------------------------------------------------------------------------------

SUBROUTINE assign_field_type(a_field,b_field)

! Subroutine to copy the contents of one field to another.
! This subroutine is used by the ASSIGNEMENT(=) operator
! i.e. it should perform an optimised equivalent of a_field = b_field

IMPLICIT NONE

CLASS(field_type), INTENT(INOUT) :: a_field
CLASS(field_type), INTENT(IN) :: b_field

CHARACTER(LEN=*), PARAMETER :: routinename = 'ASSIGN_FIELD_TYPE'

! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER :: i,j,k
INTEGER :: i_size,j_size,k_size

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! The field contains two parts which must be copied:
!  * the metadata, which can be coupied with the copy_metadata_from_field()
!    method.
!  * the actual field data, which is contained in one of four arrays - rdata,
!    idata, ldatd, or lbc_rdata

! Start with the metadata

CALL a_field%copy_metadata_from_field(b_field)

! Now the actual field data. The data arrays are allocatable, so for each array
! we must check that if the field we are copying from is allocated, that the
! field we are copying to is also allocated, and with the correct size.
! If the field we are copying from is not allocated, we must ensure the field
! we are copying to is also deallocated (if its not already).

! Deal with each of the four arrays in turn:

! 1 - field in rdata:

! is the field we are copying from (b_field) allocated?
IF (ALLOCATED(b_field%rdata)) THEN

  ! b_field is allocated - check if a_field is allocated with the correct size

  i_size = SIZE(b_field%rdata,DIM=1)
  j_size = SIZE(b_field%rdata,DIM=2)
  k_size = SIZE(b_field%rdata,DIM=3)

  IF (ALLOCATED(a_field%rdata)) THEN
    IF (      SIZE(a_field%rdata,DIM=1) /= i_size                              &
         .OR. SIZE(a_field%rdata,DIM=2) /= j_size                              &
         .OR. SIZE(a_field%rdata,DIM=3) /= k_size                              &
       ) THEN

      ! a_field is the wrong size. Deallocate it, so that it is re-allocated
      ! later with the correct size
      DEALLOCATE(a_field%rdata)

    END IF
  END IF

  ! If a_field has no allocation, allocate it now.
  IF (.NOT. ALLOCATED(a_field%rdata)) THEN
    ALLOCATE(a_field%rdata(i_size,j_size,k_size))
  END IF

  ! The allocations are now the same size - do the actual data copy.

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& PRIVATE(i,j,k)                                                          &
!$OMP& SHARED(a_field,b_field,i_size,j_size,k_size)
  DO k = 1, k_size
    DO j = 1, j_size
      DO i = 1, i_size
        a_field%rdata(i,j,k)=b_field%rdata(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

ELSE IF (ALLOCATED(a_field%rdata)) THEN

  ! b_field is not allocated, but a_field is - deallocate it

  DEALLOCATE (a_field%rdata)

END IF

! 2 - field in idata:

! is the field we are copying from (b_field) allocated?
IF (ALLOCATED(b_field%idata)) THEN

  ! b_field is allocated - check if a_field is allocated with the correct size

  i_size = SIZE(b_field%idata,DIM=1)
  j_size = SIZE(b_field%idata,DIM=2)
  k_size = SIZE(b_field%idata,DIM=3)

  IF (ALLOCATED(a_field%idata)) THEN
    IF (      SIZE(a_field%idata,DIM=1) /= i_size                              &
         .OR. SIZE(a_field%idata,DIM=2) /= j_size                              &
         .OR. SIZE(a_field%idata,DIM=3) /= k_size                              &
       ) THEN

      ! a_field is the wrong size. Deallocate it, so that it is re-allocated
      ! later with the correct size
      DEALLOCATE(a_field%idata)

    END IF
  END IF

  ! If a_field has no allocation, allocate it now.
  IF (.NOT. ALLOCATED(a_field%idata)) THEN
    ALLOCATE(a_field%idata(i_size,j_size,k_size))
  END IF

  ! The allocations are now the same size - do the actual data copy.

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& PRIVATE(i,j,k)                                                          &
!$OMP& SHARED(a_field,b_field,i_size,j_size,k_size)
  DO k = 1, k_size
    DO j = 1, j_size
      DO i = 1, i_size
        a_field%idata(i,j,k)=b_field%idata(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

ELSE IF (ALLOCATED(a_field%idata)) THEN

  ! b_field is not allocated, but a_field is - deallocate it
  DEALLOCATE (a_field%idata)

END IF

! 3 - field in ldata:

! is the field we are copying from (b_field) allocated?
IF (ALLOCATED(b_field%ldata)) THEN

  ! b_field is allocated - check if a_field is allocated with the correct size

  i_size = SIZE(b_field%ldata,DIM=1)
  j_size = SIZE(b_field%ldata,DIM=2)
  k_size = SIZE(b_field%ldata,DIM=3)

  IF (ALLOCATED(a_field%ldata)) THEN
    IF (      SIZE(a_field%ldata,DIM=1) /= i_size                              &
         .OR. SIZE(a_field%ldata,DIM=2) /= j_size                              &
         .OR. SIZE(a_field%ldata,DIM=3) /= k_size                              &
       ) THEN

      ! a_field is the wrong size. Deallocate it, so that it is re-allocated
      ! later with the correct size
      DEALLOCATE(a_field%ldata)

    END IF
  END IF

  ! If a_field has no allocation, allocate it now.
  IF (.NOT. ALLOCATED(a_field%ldata)) THEN
    ALLOCATE(a_field%ldata(i_size,j_size,k_size))
  END IF

  ! The allocations are now the same size - do the actual data copy.

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& PRIVATE(i,j,k)                                                          &
!$OMP& SHARED(a_field,b_field,i_size,j_size,k_size)
  DO k = 1, k_size
    DO j = 1, j_size
      DO i = 1, i_size
        a_field%ldata(i,j,k)=b_field%ldata(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

ELSE IF (ALLOCATED(a_field%ldata)) THEN

  ! b_field is not allocated, but a_field is - deallocate it
  DEALLOCATE (a_field%ldata)

END IF

! 4 - field in lbc_rdata:

! is the field we are copying from (b_field) allocated?
IF (ALLOCATED(b_field%lbc_rdata)) THEN

  ! b_field is allocated - check if a_field is allocated with the correct size

  i_size = SIZE(b_field%lbc_rdata,DIM=1)
  k_size = SIZE(b_field%lbc_rdata,DIM=2)

  IF (ALLOCATED(a_field%lbc_rdata)) THEN
    IF (      SIZE(a_field%lbc_rdata,DIM=1) /= i_size                          &
         .OR. SIZE(a_field%lbc_rdata,DIM=2) /= k_size                          &
       ) THEN

      ! a_field is the wrong size. Deallocate it, so that it is re-allocated
      ! later with the correct size
      DEALLOCATE(a_field%lbc_rdata)

    END IF
  END IF

  ! If a_field has no allocation, allocate it now.
  IF (.NOT. ALLOCATED(a_field%lbc_rdata)) THEN
    ALLOCATE(a_field%lbc_rdata(i_size,k_size))
  END IF

  ! The allocations are now the same size - do the actual data copy.

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& PRIVATE(i,k)                                                            &
!$OMP& SHARED(a_field,b_field,i_size,k_size)
  DO k = 1, k_size
    DO i = 1, i_size
      a_field%lbc_rdata(i,k)=b_field%lbc_rdata(i,k)
    END DO
  END DO
!$OMP END PARALLEL DO

ELSE IF (ALLOCATED(a_field%lbc_rdata)) THEN

  ! b_field is not allocated, but a_field is - deallocate it
  DEALLOCATE (a_field%lbc_rdata)

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE assign_field_type

END MODULE field_mod
