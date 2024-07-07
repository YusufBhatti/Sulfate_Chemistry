! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE calc_wind_rotation_coeff_mod

! Description:
!  Routine to calculate the coefficients used to rotate the wind fields
!  for both the input and output grid
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

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                        &
  ModuleName = 'CALC_WIND_ROTATION_COEFF_MOD'

CONTAINS

SUBROUTINE calc_wind_rotation_coeff(input_file, lbc_output_control, wind_rotation_coeff)

USE datafile_mod,               ONLY: datafile_type
USE lbc_output_control_mod,     ONLY: lbc_output_control_type
USE wind_rotation_coeff_mod,    ONLY: wind_rotation_coeff_type
USE three_dimensional_grid_mod, ONLY: three_dimensional_grid_type
USE latlon_eq_rotation_mod,     ONLY: rotate_eq_to_latlon,  &
                                      eq_latlon_vector_coeffs
USE calc_lbc_coords_mod,        ONLY: calc_lbc_coords

IMPLICIT NONE

CLASS(datafile_type), INTENT(IN) :: input_file
TYPE(lbc_output_control_type), INTENT(IN), TARGET :: lbc_output_control
TYPE(wind_rotation_coeff_type), INTENT(INOUT) :: wind_rotation_coeff

! Local variables
! Input grid coords
REAL, ALLOCATABLE :: longitude_input_coords(:)  ! Single row of longitude values
REAL, ALLOCATABLE :: latitude_input_coords(:)   ! Single column of latitude values
REAL, ALLOCATABLE :: longitude_input_standard(:)! Standard pole - longitude values for all input points
REAL, ALLOCATABLE :: latitude_input_standard(:) ! Standard pole - latitude values for all input points
REAL, ALLOCATABLE :: longitude_input_rotated(:) ! Rotated pole  - longitude values for all input points
REAL, ALLOCATABLE :: latitude_input_rotated(:)  ! Rotated pole  - latitude values for all input points
! LBC grid coords
REAL, ALLOCATABLE :: latitude_lbc_standard(:)   ! Standard pole - latitude points for LBC grid
REAL, ALLOCATABLE :: longitude_lbc_standard(:)  ! Standard pole - longitude points for LBC grid
REAL, ALLOCATABLE :: latitude_lbc_rotated(:)    ! Rotated pole  - latitude points for LBC grid
REAL, ALLOCATABLE :: longitude_lbc_rotated(:)   ! Rotated pole  - longitude points for LBC grid
REAL :: input_pole_latitude, input_pole_longitude, lbc_pole_latitude, lbc_pole_longitude
INTEGER :: input_num_points
INTEGER :: input_num_cols
INTEGER :: input_num_rows
INTEGER :: i,j,counter
INTEGER :: lbc_size
LOGICAL :: enlarged_p_grid
CHARACTER(LEN=*), PARAMETER :: routinename = 'CALC_WIND_ROTATION_COEFF'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (input_file%l_source_rotated) THEN
  ! Get coordinate information from the input grid
  longitude_input_coords = input_file%p_grid%get_longitudes()
  latitude_input_coords =  input_file%p_grid%get_latitudes()
  input_num_cols   = input_file%p_grid%get_num_cols()
  input_num_rows   = input_file%p_grid%get_num_rows()
  input_num_points = input_num_cols * input_num_rows
  ALLOCATE(longitude_input_rotated(input_num_points))
  ALLOCATE(latitude_input_rotated(input_num_points))
  counter = 0
  ! Fill in longitude values for all input points
  DO j = 1, input_num_rows 
    DO i = 1, input_num_cols
      counter = counter + 1
      longitude_input_rotated(counter) = longitude_input_coords(i)
      latitude_input_rotated(counter)  = latitude_input_coords(j)
    END DO
  END DO
  ! Calculate unrotated coordinates
  ALLOCATE(longitude_input_standard(input_num_points))
  ALLOCATE(latitude_input_standard(input_num_points))
  CALL input_file%p_grid%get_pole(input_pole_latitude, input_pole_longitude)
  CALL rotate_eq_to_latlon(latitude_input_rotated, longitude_input_rotated,  &
                           latitude_input_standard, longitude_input_standard,&
                           input_pole_latitude, input_pole_longitude,        &
                           input_num_points)
  ! Calculate input wind rotation coefficients used to move from a rotated grid to the
  ! standard un-rotated pole
  ALLOCATE(wind_rotation_coeff%input_wind_coeff1(input_num_points))
  ALLOCATE(wind_rotation_coeff%input_wind_coeff2(input_num_points))
  CALL eq_latlon_vector_coeffs(wind_rotation_coeff%input_wind_coeff1,  &
                               wind_rotation_coeff%input_wind_coeff2,  &
                 longitude_input_standard, longitude_input_rotated,    &
                 input_pole_latitude, input_pole_longitude, input_num_points)
  DEALLOCATE(latitude_input_standard)
  DEALLOCATE(longitude_input_standard)
  DEALLOCATE(latitude_input_rotated)
  DEALLOCATE(longitude_input_rotated)
  DEALLOCATE(latitude_input_coords)
  DEALLOCATE(longitude_input_coords)
END IF ! input_file%l_source_rotated

IF (lbc_output_control%l_target_rotated) THEN
  ! LBC wind is rotated whilst on an enlarged p grid rather than u/v
  enlarged_p_grid = .TRUE.
  ! Get the latitude/longitude values and lbc size of the target grid
  ! These will be returned in the format of the four lbc sectors
  CALL calc_lbc_coords(enlarged_p_grid, lbc_output_control,       &
       lbc_output_control%p_grid_enlarged,                        &
       latitude_lbc_rotated, longitude_lbc_rotated, lbc_size)
  ! Get standard non-rotated latitude and longitude of lbc grid
  ALLOCATE(latitude_lbc_standard(lbc_size))
  ALLOCATE(longitude_lbc_standard(lbc_size))
  CALL rotate_eq_to_latlon(latitude_lbc_rotated, longitude_lbc_rotated,   &
                           latitude_lbc_standard, longitude_lbc_standard, &
               lbc_output_control%pole_lat, lbc_output_control%pole_long, &
               lbc_size)
  ALLOCATE(wind_rotation_coeff%output_wind_coeff1(lbc_size))
  ALLOCATE(wind_rotation_coeff%output_wind_coeff2(lbc_size))
  CALL eq_latlon_vector_coeffs(wind_rotation_coeff%output_wind_coeff1, &
                               wind_rotation_coeff%output_wind_coeff2, &
            longitude_lbc_standard, longitude_lbc_rotated,             &
            lbc_output_control%pole_lat, lbc_output_control%pole_long, &
            lbc_size)
  DEALLOCATE(longitude_lbc_standard)
  DEALLOCATE(latitude_lbc_standard)
END IF ! lbc_output_control%l_target_rotated
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE calc_wind_rotation_coeff

END MODULE calc_wind_rotation_coeff_mod
