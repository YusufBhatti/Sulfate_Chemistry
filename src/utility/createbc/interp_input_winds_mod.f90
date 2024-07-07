! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE interp_input_winds_mod

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.

! DrHook Modules
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'INTERP_INPUT_WINDS_MOD'

CONTAINS

!-------------------------------------------------------------------------------

SUBROUTINE interp_input_wind_p_to_u(input_u_field, input_file)
! Description:
!  Perform a linear interpolation of the U wind from the P grid 
!  back to the U grid.
USE datafile_mod,              ONLY: datafile_type, endgame, new_dynamics
USE field_mod,                 ONLY: field_type, ASSIGNMENT(=)
USE stashmaster_constants_mod, ONLY: p_points, u_points
USE ereport_mod,               ONLY: ereport
USE errormessagelength_mod,    ONLY: errormessagelength

IMPLICIT NONE
! Arguments
TYPE(field_type), INTENT(INOUT) :: input_u_field
CLASS(datafile_type), INTENT(INOUT) :: input_file
! Local variables
INTEGER :: column, row, level, icode
REAL, ALLOCATABLE :: temp_field_data(:,:,:)
CHARACTER(LEN=*), PARAMETER :: routinename = 'INTERP_INPUT_WIND_P_TO_U'
CHARACTER(LEN=errormessagelength) :: cmessage
REAL, ALLOCATABLE :: weight1(:), weight2(:) ! Linear interpolation weights
REAL, ALLOCATABLE :: p_grid_long(:) ! Longitude arrays
REAL, ALLOCATABLE :: u_grid_long(:)
INTEGER :: p_grid_cols, p_grid_rows
INTEGER :: u_grid_cols, u_grid_rows
INTEGER :: num_levels
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Check that U wind field is on P grid
IF (input_u_field%get_horiz_grid_type() /= p_points) THEN
  icode = 10
  WRITE(cmessage, '(A,I0,A,I0)') "Expected grid type of U field to be p_points: ", &
       p_points, " but field has grid type code: ", input_u_field%get_horiz_grid_type()
  CALL ereport(routinename, icode, cmessage)
END IF

! Get grid dimensions
p_grid_long = input_u_field%get_longitudes()
u_grid_long = input_file%u_grid%get_longitudes()
p_grid_cols = input_u_field%get_num_cols()
p_grid_rows = input_u_field%get_num_rows()
u_grid_rows = input_file%u_grid%get_num_rows()
u_grid_cols = input_file%u_grid%get_num_cols()
num_levels  = input_u_field%get_num_levels()

ALLOCATE(weight1(u_grid_cols))
ALLOCATE(weight2(u_grid_cols))
ALLOCATE(temp_field_data(u_grid_cols, u_grid_rows, num_levels))
IF (input_file%grid_staggering == new_dynamics) THEN
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(level, weight2, weight1, column, row)     &
!$OMP SHARED(num_levels, u_grid_cols, u_grid_rows, p_grid_long,                   &
!$OMP u_grid_long, temp_field_data, input_u_field)
  DO level = 1, num_levels
    DO row = 1, u_grid_rows
      DO column = 1, u_grid_cols-1
        ! Only calc longitude interp weights for first row
        IF (row == 1) THEN
          weight2(column) = (u_grid_long(column) - p_grid_long(column)) / &
               (p_grid_long(column+1) - p_grid_long(column))
          weight1(column) = 1 - weight2(column)
        END IF
        temp_field_data(column, row, level) =                             &
             weight1(column) * input_u_field%rdata(column, row, level)    &
           + weight2(column) * input_u_field%rdata(column+1, row, level)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  ! Copy last U column values from to last P column
  temp_field_data(u_grid_cols, :, :) = input_u_field%rdata(u_grid_cols, :, :)
ELSE IF (input_file%grid_staggering == endgame) THEN
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(level, weight2, weight1, column, row) &
!$OMP SHARED(num_levels, u_grid_cols, u_grid_rows, p_grid_long,               &
!$OMP u_grid_long, temp_field_data, input_u_field)
  DO level = 1, num_levels
    DO row = 1, u_grid_rows
      DO column = 2, u_grid_cols
        ! Only calc longitude interp weights for first row
        IF (row == 1) THEN
          weight2(column) = (u_grid_long(column) - p_grid_long(column-1)) / &
               (p_grid_long(column) - p_grid_long(column-1))
          weight1(column) = 1 - weight2(column)
        END IF
        temp_field_data(column, row, level) =                               &
               weight1(column) * input_u_field%rdata(column-1, row, level)  &
             + weight2(column) * input_u_field%rdata(column, row, level)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  ! Copy first U column values from to first P column
  temp_field_data(1, :, :) = input_u_field%rdata(1, :, :)
ELSE
  icode = 20
  WRITE(cmessage, '(A,I0)') "Unsupported grid staggering value: ", input_file%grid_staggering
  CALL ereport(routinename, icode, cmessage)
END IF

! Replace P grid data with newly interpolated U grid values 
CALL MOVE_ALLOC(temp_field_data, input_u_field%rdata)
! Set field latitudes and longitudes to U grid values
CALL input_u_field%set_latitudes(input_file%u_grid%get_latitudes())
CALL input_u_field%set_longitudes(input_file%u_grid%get_longitudes())
CALL input_u_field%set_num_rows(input_file%u_grid%get_num_rows())
CALL input_u_field%set_num_cols(input_file%u_grid%get_num_cols())
CALL input_u_field%set_horiz_grid_type(u_points)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE interp_input_wind_p_to_u

!-------------------------------------------------------------------------------

SUBROUTINE interp_input_wind_p_to_v(input_v_field, input_file)
! Description:
!  Perform a linear interpolation of the V wind from the P grid 
!  back to the V grid.
USE datafile_mod,              ONLY: datafile_type, endgame, new_dynamics
USE field_mod,                 ONLY: field_type, ASSIGNMENT(=)
USE stashmaster_constants_mod, ONLY: p_points, v_points
USE ereport_mod,               ONLY: ereport
USE errormessagelength_mod,    ONLY: errormessagelength

IMPLICIT NONE
! Arguments
TYPE(field_type), INTENT(INOUT) :: input_v_field
CLASS(datafile_type), INTENT(INOUT) :: input_file
! Local variables
INTEGER :: column, row, level, icode
REAL, ALLOCATABLE :: temp_field_data(:,:,:)
CHARACTER(LEN=*), PARAMETER :: routinename = 'INTERP_INPUT_WIND_P_TO_V'
CHARACTER(LEN=errormessagelength) :: cmessage
REAL, ALLOCATABLE :: weight1(:), weight2(:) ! Linear interpolation weights
REAL, ALLOCATABLE :: p_grid_lat(:) ! Latitude arrays
REAL, ALLOCATABLE :: v_grid_lat(:)
INTEGER :: p_grid_cols, p_grid_rows
INTEGER :: v_grid_cols, v_grid_rows
INTEGER :: num_levels
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Check that U wind field is on P grid
IF (input_v_field%get_horiz_grid_type() /= p_points) THEN
  icode = 20
  WRITE(cmessage, '(A,I0,A,I0)') "Expected grid type of V field to be p_points: ", &
       p_points, " but field has grid type code: ", input_v_field%get_horiz_grid_type()
  CALL ereport(routinename, icode, cmessage)
END IF

! Get grid dimensions
p_grid_lat  = input_v_field%get_latitudes()
v_grid_lat  = input_file%v_grid%get_latitudes()
p_grid_cols = input_v_field%get_num_cols()
p_grid_rows = input_v_field%get_num_rows()
v_grid_rows = input_file%v_grid%get_num_rows()
v_grid_cols = input_file%v_grid%get_num_cols()
num_levels  = input_v_field%get_num_levels()

ALLOCATE(weight1(v_grid_rows))
ALLOCATE(weight2(v_grid_rows))
ALLOCATE(temp_field_data(v_grid_cols, v_grid_rows, num_levels))
IF (input_file%grid_staggering == new_dynamics) THEN
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(level, weight2, weight1, column, row)     &
!$OMP SHARED(num_levels, v_grid_cols, v_grid_rows, p_grid_lat,                    &
!$OMP v_grid_lat, temp_field_data, input_v_field)
  DO level = 1, num_levels
    DO column = 1, v_grid_cols
      DO row = 1, v_grid_rows
        ! Only calc latitude interp weights for first column
        IF (column == 1) THEN
          weight2(row) = (v_grid_lat(row) - p_grid_lat(row)) /                &
               (p_grid_lat(row+1) - p_grid_lat(row))
          weight1(row) = 1 - weight2(row)
        END IF
        temp_field_data(column, row, level) =                                 &
               weight1(row) * input_v_field%rdata(column, row, level)         &
             + weight2(row) * input_v_field%rdata(column, row+1, level)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE IF (input_file%grid_staggering == endgame) THEN
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(level, weight2, weight1, column, row) &
!$OMP SHARED(num_levels, v_grid_cols, v_grid_rows, p_grid_lat,                &
!$OMP v_grid_lat, temp_field_data, input_v_field)
  DO level = 1, num_levels
    DO column = 1, v_grid_cols
      DO row = 2, v_grid_rows-1
         ! Only calc latitude interp weights for first column
        IF (column == 1) THEN
          weight2(row) = (v_grid_lat(row) - p_grid_lat(row-1)) /              &
               (p_grid_lat(row) - p_grid_lat(row-1))
          weight1(row) = 1 - weight2(row)
        END IF
        temp_field_data(column, row, level) =                                 &
               weight1(row) * input_v_field%rdata(column, row-1, level)       &
             + weight2(row) * input_v_field%rdata(column, row, level)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  ! Copy for start and end points
  temp_field_data(:, 1, :) = input_v_field%rdata(:, 1, :)
  temp_field_data(:, v_grid_rows, :) = input_v_field%rdata(:, p_grid_rows, :)
ELSE
  icode = 20
  WRITE(cmessage, '(A,I0)') "Unsupported grid staggering value: ", input_file%grid_staggering
  CALL ereport(routinename, icode, cmessage)
END IF

! Replace P grid data with newly interpolated V grid values
CALL MOVE_ALLOC(temp_field_data, input_v_field%rdata)
! Set field latitudes and longitudes to V grid values
CALL input_v_field%set_latitudes(input_file%v_grid%get_latitudes())
CALL input_v_field%set_longitudes(input_file%v_grid%get_longitudes())
CALL input_v_field%set_num_rows(input_file%v_grid%get_num_rows())
CALL input_v_field%set_num_cols(input_file%v_grid%get_num_cols())
CALL input_v_field%set_horiz_grid_type(v_points)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE interp_input_wind_p_to_v

SUBROUTINE interp_input_wind_u_to_p(input_u_field, input_file)
! Description:
!  Perform a linear interpolation of the U wind to the P grid. Any 
!  P column which is not surrounded by U points will have its values
!  taken from the neighbouring column.  This is sufficient as when 
!  generating LBC areas one expects the LBC area to be not be at the
!  edge of the source domain.
USE datafile_mod,              ONLY: datafile_type, endgame, new_dynamics
USE field_mod,                 ONLY: field_type, ASSIGNMENT(=)
USE stashmaster_constants_mod, ONLY: u_points, p_points
USE ereport_mod,               ONLY: ereport
USE errormessagelength_mod,    ONLY: errormessagelength

IMPLICIT NONE

! Arguments
TYPE(field_type), INTENT(INOUT) :: input_u_field
CLASS(datafile_type), INTENT(INOUT) :: input_file
! Local variables
INTEGER :: column, row, level, icode
REAL, ALLOCATABLE :: temp_field_data(:,:,:)
CHARACTER(LEN=*), PARAMETER :: routinename = 'INTERP_INPUT_WIND_U_TO_P'
CHARACTER(LEN=errormessagelength) :: cmessage
REAL, ALLOCATABLE :: weight1(:), weight2(:) ! Linear interpolation weights
REAL, ALLOCATABLE :: p_grid_long(:) ! Longitude arrays
REAL, ALLOCATABLE :: u_grid_long(:)
INTEGER :: p_grid_cols, p_grid_rows
INTEGER :: u_grid_cols, u_grid_rows
INTEGER :: num_levels
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL timer( 'u_wind_to_p_grid', 5)
! Check that field is on U grid
IF (input_u_field%get_horiz_grid_type() /= u_points) THEN
  icode = 10
  WRITE(cmessage, '(A,I0,A,I0)') "Expected grid type of U field to be u_points: ", &
       u_points, " but field has grid type code: ", input_u_field%get_horiz_grid_type()
  CALL ereport(routinename, icode, cmessage)
END IF

! Get grid dimensions
u_grid_long = input_u_field%get_longitudes()
p_grid_long = input_file%p_grid%get_longitudes()
u_grid_cols = input_u_field%get_num_cols()
u_grid_rows = input_u_field%get_num_rows()
p_grid_rows = input_file%p_grid%get_num_rows()
p_grid_cols = input_file%p_grid%get_num_cols()
num_levels  = input_u_field%get_num_levels()

ALLOCATE(weight1(p_grid_cols))
ALLOCATE(weight2(p_grid_cols))
ALLOCATE(temp_field_data(p_grid_cols, p_grid_rows, num_levels))

IF (input_file%grid_staggering == new_dynamics) THEN
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(level, weight2, weight1, column, row)&
!$OMP SHARED(num_levels, u_grid_cols, u_grid_rows, p_grid_long,              &
!$OMP u_grid_long, temp_field_data, input_u_field)
  DO level = 1, num_levels
    DO row = 1, u_grid_rows
      DO column = 2, u_grid_cols
        ! Only calc longitude interp weights for first row
        IF (row == 1) THEN
          weight2(column) = (p_grid_long(column) - u_grid_long(column-1)) / &
               (u_grid_long(column) - u_grid_long(column-1))
          weight1(column) = 1 - weight2(column)
        END IF
        temp_field_data(column, row, level) =                               &
               weight1(column) * input_u_field%rdata(column-1, row, level)  &
             + weight2(column) * input_u_field%rdata(column, row, level)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  ! Direct copy for first P column values from first U column
  temp_field_data(1, :, :) = input_u_field%rdata(1, :, :)
      
ELSE IF (input_file%grid_staggering == endgame) THEN
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(level, weight2, weight1, column, row) &
!$OMP SHARED(num_levels, u_grid_cols, u_grid_rows, p_grid_long,               &
!$OMP u_grid_long, temp_field_data, input_u_field)
  DO level = 1, num_levels
    DO row = 1, u_grid_rows
      DO column = 1, u_grid_cols-1
        ! Only calc longitude interp weights for first row
        IF (row == 1) THEN
          weight2(column) = (p_grid_long(column) - u_grid_long(column)) /     &
               (u_grid_long(column+1) - u_grid_long(column))
          weight1(column) = 1 - weight2(column)
        END IF
        temp_field_data(column, row, level) =                                 &
               weight1(column) * input_u_field%rdata(column, row, level)      &
             + weight2(column) * input_u_field%rdata(column+1, row, level)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  ! Copy last P column values from to last U column
  temp_field_data(p_grid_cols, :, :) = input_u_field%rdata(u_grid_cols, :, :)
ELSE
  icode = 20
  WRITE(cmessage, '(A,I0)') "Unsupported grid staggering value: ", input_file%grid_staggering
  CALL ereport(routinename, icode, cmessage)
END IF

! Replace original U field data with P grid values 
CALL MOVE_ALLOC(temp_field_data, input_u_field%rdata)
CALL input_u_field%set_latitudes(input_file%p_grid%get_latitudes())
CALL input_u_field%set_longitudes(input_file%p_grid%get_longitudes())
CALL input_u_field%set_num_rows(input_file%p_grid%get_num_rows())
CALL input_u_field%set_num_cols(input_file%p_grid%get_num_cols())
CALL input_u_field%set_horiz_grid_type(p_points)

CALL timer( 'u_wind_to_p_grid', 6)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE interp_input_wind_u_to_p

!-------------------------------------------------------------------------------

SUBROUTINE interp_input_wind_v_to_p(input_v_field, input_file)
! Description:
!  Perform a linear interpolation of the U wind to the P grid. Any 
!  P column which is not surrounded by U points will have its values
!  taken from the neighbouring column.  This is sufficient as when 
!  generating LBC areas one expects the LBC area to be not be at the
!  edge of the source domain.
USE datafile_mod,              ONLY: datafile_type, endgame, new_dynamics
USE field_mod,                 ONLY: field_type, ASSIGNMENT(=)
USE stashmaster_constants_mod, ONLY: v_points, p_points
USE ereport_mod,               ONLY: ereport
USE errormessagelength_mod,    ONLY: errormessagelength

IMPLICIT NONE

! Arguments
TYPE(field_type), INTENT(INOUT) :: input_v_field
CLASS(datafile_type), INTENT(INOUT) :: input_file
! Local variables
INTEGER :: column, row, level, icode
REAL, ALLOCATABLE :: temp_field_data(:,:,:)
CHARACTER(LEN=*), PARAMETER :: routinename = 'INTERP_INPUT_WIND_V_TO_P'
CHARACTER(LEN=errormessagelength) :: cmessage
REAL, ALLOCATABLE :: weight1(:), weight2(:) ! Linear interpolation weights
REAL, ALLOCATABLE :: p_grid_lat(:) ! Latitude arrays
REAL, ALLOCATABLE :: v_grid_lat(:)
INTEGER :: p_grid_cols, p_grid_rows
INTEGER :: v_grid_cols, v_grid_rows
INTEGER :: num_levels
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL timer( 'v_wind_to_p_grid', 5)

! Check that field is on U grid
IF (input_v_field%get_horiz_grid_type() /= v_points) THEN
  icode = 10
  WRITE(cmessage, '(A,I0,A,I0)') "Expected grid type of V field to be v_points: ", &
       v_points, " but field has grid type code: ", input_v_field%get_horiz_grid_type()
  CALL ereport(routinename, icode, cmessage)
END IF

! Get grid dimensions
v_grid_lat  = input_v_field%get_latitudes()
p_grid_lat  = input_file%p_grid%get_latitudes()
v_grid_cols = input_v_field%get_num_cols()
v_grid_rows = input_v_field%get_num_rows()
p_grid_cols = input_file%p_grid%get_num_cols()
p_grid_rows = input_file%p_grid%get_num_rows()
num_levels  = input_v_field%get_num_levels()

ALLOCATE(weight1(p_grid_rows))
ALLOCATE(weight2(p_grid_rows))
ALLOCATE(temp_field_data(p_grid_cols, p_grid_rows, num_levels))

IF (input_file%grid_staggering == new_dynamics) THEN
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(level, weight2, weight1, column, row)     &
!$OMP SHARED(num_levels, v_grid_cols, v_grid_rows, p_grid_lat,     &
!$OMP v_grid_lat, temp_field_data, input_v_field)
  DO level = 1, num_levels
    DO column = 1, v_grid_cols
      DO row = 2, v_grid_rows
        ! Only calc latitude interp weights for first column
        IF (column == 1) THEN
          weight2(row) = (p_grid_lat(row) - v_grid_lat(row-1)) / &
               (v_grid_lat(row) - v_grid_lat(row-1))
          weight1(row) = 1 - weight2(row)
        END IF
        temp_field_data(column, row, level) =                              &
               weight1(row) * input_v_field%rdata(column, row-1, level)    &
             + weight2(row) * input_v_field%rdata(column, row, level)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  ! Copy for start and end points
  temp_field_data(:, 1, :) = input_v_field%rdata(:, 1, :)
  temp_field_data(:, p_grid_rows, :) = input_v_field%rdata(:, v_grid_rows, :)
ELSE IF (input_file%grid_staggering == endgame) THEN
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(level, weight2, weight1, column, row)  &
!$OMP SHARED(num_levels, v_grid_cols, v_grid_rows, p_grid_lat,     &
!$OMP v_grid_lat, temp_field_data, input_v_field)
  DO level = 1, num_levels
    DO column = 1, v_grid_cols
      DO row = 1, v_grid_rows-1
         ! Only calc latitude interp weights for first column
        IF (column == 1) THEN
          weight2(row) = (p_grid_lat(row) - v_grid_lat(row)) / &
               (v_grid_lat(row+1) - v_grid_lat(row))
          weight1(row) = 1 - weight2(row)
        END IF
        temp_field_data(column, row, level) =                              &
               weight1(row) * input_v_field%rdata(column, row, level)      &
             + weight2(row) * input_v_field%rdata(column, row+1, level)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  icode = 20
  WRITE(cmessage, '(A,I0)') "Unsupported grid staggering value: ", input_file%grid_staggering
  CALL ereport(routinename, icode, cmessage)
END IF

! Replace original V field data with P grid values
CALL MOVE_ALLOC(temp_field_data, input_v_field%rdata)
! Deallocate field V grid latitudes and latitudes
CALL input_v_field%set_latitudes(input_file%p_grid%get_latitudes())
CALL input_v_field%set_longitudes(input_file%p_grid%get_longitudes())
CALL input_v_field%set_num_rows(input_file%p_grid%get_num_rows())
CALL input_v_field%set_num_cols(input_file%p_grid%get_num_cols())
CALL input_v_field%set_horiz_grid_type(p_points)

CALL timer( 'v_wind_to_p_grid', 6)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE interp_input_wind_v_to_p

!-------------------------------------------------------------------------------

END MODULE interp_input_winds_mod
