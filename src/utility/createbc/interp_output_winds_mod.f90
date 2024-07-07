! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE interp_output_winds_mod
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.
USE interp_weights_mod,        ONLY: interp_weights_type, p_to_enlarged_p, u_to_enlarged_p, &
                                     v_to_enlarged_p, weights_index
USE field_mod,                 ONLY: field_type, ASSIGNMENT(=)
USE lbc_output_control_mod,    ONLY: lbc_output_control_type, new_dynamics, endgame
USE stashmaster_constants_mod, ONLY: p_points, u_points, v_points
USE um_parparams,              ONLY: halo_type_extended
USE missing_data_mod,          ONLY: rmdi
! DrHook Modules
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE
PRIVATE
PUBLIC :: interp_output_wind_p_to_u, interp_output_wind_p_to_v

CHARACTER(LEN=*), PARAMETER :: ModuleName = 'INTERP_OUTPUT_WINDS_MOD'

CONTAINS

SUBROUTINE interp_output_wind_p_to_u(output_u_field, lbc_output_control, interp_weights)
! Description:
!  Interpolate the output U wind field from the enlarged P grid to the 
!  standard LBC U grid.
!
IMPLICIT NONE
! Arguments
TYPE(lbc_output_control_type),   TARGET, INTENT(IN) :: lbc_output_control
TYPE(field_type), INTENT(INOUT) :: output_u_field
CLASS(interp_weights_type), INTENT(IN) :: interp_weights(:,:)
! Local variables
REAL, ALLOCATABLE :: p_grid_long(:) ! Longitude arrays
REAL, ALLOCATABLE :: u_grid_long(:)
REAL, ALLOCATABLE :: weight1(:)
REAL, ALLOCATABLE :: weight2(:)
INTEGER :: p_grid_cols, p_grid_rows
INTEGER :: u_grid_cols, u_grid_rows
INTEGER :: halo_ns, halo_ew, rim_width
INTEGER :: num_levels
INTEGER :: lbc_u_level_size
INTEGER :: counter_u_lbc, counter_p_lbc, k, i, j
REAL, ALLOCATABLE :: temp_field_data(:,:)
CHARACTER(LEN=*), PARAMETER :: routinename = 'INTERP_OUTPUT_WIND_P_TO_U'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Get longitudes of U and P grids
! output_u_field is on the P grid
p_grid_long = output_u_field%get_longitudes()
u_grid_long = lbc_output_control%u_grid(halo_type_extended)%get_longitudes()
p_grid_cols = output_u_field%get_num_cols()
p_grid_rows = output_u_field%get_num_rows()
u_grid_cols = lbc_output_control%u_grid(halo_type_extended)%get_num_cols()
u_grid_rows = lbc_output_control%u_grid(halo_type_extended)%get_num_rows()
num_levels  = output_u_field%get_num_levels()
halo_ns     = lbc_output_control%u_grid(halo_type_extended)%get_halo_ns()
halo_ew     = lbc_output_control%u_grid(halo_type_extended)%get_halo_ew()
rim_width   = lbc_output_control%rim_width

! Get LBC level size of U field grid with extended haloes
lbc_u_level_size = SIZE(interp_weights(weights_index(u_points), &
                                       halo_type_extended)%bilinear_index_b_l)
ALLOCATE(temp_field_data(lbc_u_level_size, num_levels))
temp_field_data = rmdi

! Calculate longitude linear weights
ALLOCATE(weight1(u_grid_cols))
ALLOCATE(weight2(u_grid_cols))
DO i = 1, u_grid_cols
  weight2(i) = (u_grid_long(i) - p_grid_long(i)) / &
       (p_grid_long(i+1) - p_grid_long(i))
  weight1(i) = 1 - weight2(i)
END DO

! Each level of data is contained in a 1D array, indexed following the LBC structure
! described in C71. As the enlarged P grid and the u grid have a different LBC level
! size, care must be taken to make sure the correct LBC data points are accessed
! when performing the linear interpolation. Each LBC section will have more P points
! than U points.

! Linearly interpolate each LBC section (N,E,S,W) seperately.

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(k, i, j, counter_u_lbc, counter_p_lbc)      &
!$OMP SHARED(num_levels, u_grid_cols, u_grid_rows, rim_width, halo_ns, halo_ew,     &
!$OMP weight1, weight2, temp_field_data, p_grid_cols, output_u_field,               &
!$OMP lbc_output_control)
IF (lbc_output_control%output_grid_stagger == new_dynamics) THEN
!$OMP DO
  DO k = 1, num_levels
    counter_u_lbc = 0
    ! Linearly interpolate each LBC section (N,E,S,W) seperately.

    ! North LBC section:
    ! The enlarged P grid has extra row of points on the bottom of the north section
    ! these are for the benefit of interpolating the V field and are not needed by 
    ! the U field. We need to skip this row and start at the beginning of the next 
    ! row.
    counter_p_lbc = p_grid_cols
    !  The enlarged P grid has extra row of points on the bottom of the north section.
    DO j = 1, rim_width + halo_ns
      DO i = 1, u_grid_cols
        counter_u_lbc = counter_u_lbc + 1
        counter_p_lbc = counter_p_lbc + 1
        temp_field_data(counter_u_lbc, k) =                              &
             (weight1(i) * output_u_field%lbc_rdata(counter_p_lbc, k)) + &
             (weight2(i) * output_u_field%lbc_rdata(counter_p_lbc+1, k))
      END DO
      ! One more P point per row
      counter_p_lbc = counter_p_lbc + 1
    END DO
    ! East LBC section:
    !  The enlarged P grid has extra column of points on the left of the east section.
    DO j = 1, u_grid_rows - (2*rim_width) - (2*halo_ns)
      DO i = 1, rim_width + halo_ew
        counter_u_lbc = counter_u_lbc + 1
        counter_p_lbc = counter_p_lbc + 1
        temp_field_data(counter_u_lbc, k) =                              &
             (weight1(u_grid_cols-rim_width-halo_ew+i) *                 &
             output_u_field%lbc_rdata(counter_p_lbc, k)) +               &
             (weight2(u_grid_cols-rim_width-halo_ew+i) *                 &
             output_u_field%lbc_rdata(counter_p_lbc+1, k))
      END DO
      ! One more P point per row
      counter_p_lbc = counter_p_lbc + 1
    END DO
    ! South LBC section:
    !  The enlarged P grid has extra row of points on the top of the south section.
    DO j = 1, rim_width + halo_ns
      DO i = 1, u_grid_cols
        counter_u_lbc = counter_u_lbc + 1
        counter_p_lbc = counter_p_lbc + 1
        temp_field_data(counter_u_lbc, k) =                              &
             (weight1(i) * output_u_field%lbc_rdata(counter_p_lbc, k)) + & 
             (weight2(i) * output_u_field%lbc_rdata(counter_p_lbc+1, k))
      END DO
      ! One more P point per row
      counter_p_lbc = counter_p_lbc + 1
    END DO
    ! The enlarged P grid has extra row of points on the top of the south section
    ! these are for the benefit of interpolating the V field and are not needed by 
    ! the U field. We need to skip this row before moving onto the west section.
    counter_p_lbc = counter_p_lbc + p_grid_cols
    ! West LBC section:
    !  The enlarged P grid has extra column of points on the left of the west section.
    DO j = 1, u_grid_rows - (2*rim_width) - (2*halo_ns)
      DO i = 1, rim_width + halo_ew
        counter_u_lbc = counter_u_lbc + 1
        counter_p_lbc = counter_p_lbc + 1
        temp_field_data(counter_u_lbc, k) =                              &
             (weight1(i) * output_u_field%lbc_rdata(counter_p_lbc, k)) + &
             (weight2(i) * output_u_field%lbc_rdata(counter_p_lbc+1, k))
      END DO
      ! One more P point per row
      counter_p_lbc = counter_p_lbc + 1
    END DO
  END DO
!$OMP END DO
ELSE IF (lbc_output_control%output_grid_stagger == endgame) THEN
  ! The EG enlarged P grid object an extra row on north and south and 
  ! extra column on west as such the num_rows is increased by 2 and num
  ! cols increased by 1 compared to the standard P grid.
!$OMP DO
  DO k = 1, num_levels
    counter_u_lbc = 0
    counter_p_lbc = 0
    ! North LBC section:
    !  The enlarged P grid has extra row of points on the top of the north section and
    !  an extra column on left of north section.
    DO j = 1, rim_width + halo_ns
      DO i = 1, u_grid_cols
        counter_u_lbc = counter_u_lbc + 1
        counter_p_lbc = counter_p_lbc + 1
        temp_field_data(counter_u_lbc, k) =                              &
             (weight1(i) * output_u_field%lbc_rdata(counter_p_lbc, k)) + &
             (weight2(i) * output_u_field%lbc_rdata(counter_p_lbc+1, k))
      END DO
      ! One more P point per row
      counter_p_lbc = counter_p_lbc + 1
    END DO
    ! Need to skip the extra row of P points before we move onto the east section
    counter_p_lbc = counter_p_lbc + p_grid_cols
    ! East LBC section:
    !  The enlarged LBC grid east section has extra row at top and bottom and 
    !  and extra column on the left.
    ! First skip the extra row at bottom as this is only needed for V grid interpolation
    counter_p_lbc = counter_p_lbc + rim_width + halo_ew + 1
    DO j = 1, u_grid_rows - (2*rim_width) - (2*halo_ns)
      DO i = 1, rim_width + halo_ew
        counter_u_lbc = counter_u_lbc + 1
        counter_p_lbc = counter_p_lbc + 1
        temp_field_data(counter_u_lbc, k) =                              &
             (weight1(u_grid_cols-rim_width-halo_ew+i) *                 &
             output_u_field%lbc_rdata(counter_p_lbc, k)) +               &
             (weight2(u_grid_cols-rim_width-halo_ew+i) *                 &
             output_u_field%lbc_rdata(counter_p_lbc+1, k))
      END DO
      ! One more P point per row
      counter_p_lbc = counter_p_lbc + 1
    END DO
    ! Need to skip the extra row of P points before we move onto the south section
    counter_p_lbc = counter_p_lbc + rim_width + halo_ew + 1
    ! South LBC section:
    !  The enlarged P grid south section has an extra row of points on bottom edge and extra
    !  column on left side
    ! Need to skip the extra row on the bottom as these are only needed for V interpolation
    counter_p_lbc =  counter_p_lbc + p_grid_cols
    DO j = 1, rim_width + halo_ns
      DO i = 1, u_grid_cols
        counter_u_lbc = counter_u_lbc + 1
        counter_p_lbc = counter_p_lbc + 1
        temp_field_data(counter_u_lbc, k) =                              &
             (weight1(i) * output_u_field%lbc_rdata(counter_p_lbc, k)) + &
             (weight2(i) * output_u_field%lbc_rdata(counter_p_lbc+1, k))
      END DO
      ! One more P point per row
      counter_p_lbc = counter_p_lbc + 1
    END DO
    ! West LBC section:
    !  The enlarged LBC grid west section has extra row at top and bottom and 
    !  and extra column on the left.
    ! First skip the extra row at bottom as this is only needed for V grid interpolation
    counter_p_lbc = counter_p_lbc + rim_width + halo_ew + 1
    DO j = 1, u_grid_rows - (2*rim_width) - (2*halo_ns)
      DO i = 1, rim_width + halo_ew
        counter_u_lbc = counter_u_lbc + 1
        counter_p_lbc = counter_p_lbc + 1
        temp_field_data(counter_u_lbc, k) =                              &
             (weight1(i) * output_u_field%lbc_rdata(counter_p_lbc, k)) + &
             (weight2(i) * output_u_field%lbc_rdata(counter_p_lbc+1, k))
      END DO
      ! One more P point per row
      counter_p_lbc = counter_p_lbc + 1
    END DO
  END DO
!$OMP END DO
END IF
!$OMP END PARALLEL

DEALLOCATE(weight1)
DEALLOCATE(weight2)

! Replace P grid data with newly interpolated U grid values
CALL MOVE_ALLOC(temp_field_data, output_u_field%lbc_rdata)
! Set field latitudes and longitudes to U grid values
CALL output_u_field%set_latitudes(lbc_output_control%u_grid(halo_type_extended)%get_latitudes())
CALL output_u_field%set_longitudes(lbc_output_control%u_grid(halo_type_extended)%get_longitudes())
CALL output_u_field%set_num_rows(u_grid_rows)
CALL output_u_field%set_num_cols(u_grid_cols)
CALL output_u_field%set_horiz_grid_type(u_points)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE interp_output_wind_p_to_u


SUBROUTINE interp_output_wind_p_to_v(output_v_field, lbc_output_control, interp_weights)
! Description:
!  Interpolate the output V wind field from the enlarged P grid to the 
!  standard LBC V grid.
!
IMPLICIT NONE
TYPE(lbc_output_control_type),   TARGET, INTENT(IN) :: lbc_output_control
TYPE(field_type), INTENT(INOUT) :: output_v_field
CLASS(interp_weights_type), INTENT(IN) :: interp_weights(:,:)
! Local variables
REAL, ALLOCATABLE :: p_grid_lat(:) ! Latitude arrays
REAL, ALLOCATABLE :: v_grid_lat(:)
REAL, ALLOCATABLE :: weight1(:)
REAL, ALLOCATABLE :: weight2(:)
INTEGER :: p_grid_cols, p_grid_rows
INTEGER :: v_grid_cols, v_grid_rows
INTEGER :: halo_ns, halo_ew, rim_width
INTEGER :: num_levels
INTEGER :: lbc_v_level_size
INTEGER :: north_offset, east_offset ! Offsets needed to determine LBC index of the P point directly
INTEGER :: south_offset, west_offset ! north of the current V point
INTEGER :: counter_v_lbc, counter_p_lbc, k, i, j
REAL, ALLOCATABLE :: temp_field_data(:,:)
CHARACTER(LEN=*), PARAMETER :: routinename = 'INTERP_OUTPUT_WIND_P_TO_V'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Get latitudes of V and P grids
! output_v_field is on the P grid
p_grid_lat  = output_v_field%get_latitudes()
v_grid_lat  = lbc_output_control%v_grid(halo_type_extended)%get_latitudes()
p_grid_cols = output_v_field%get_num_cols()
p_grid_rows = output_v_field%get_num_rows()
v_grid_cols = lbc_output_control%v_grid(halo_type_extended)%get_num_cols()
v_grid_rows = lbc_output_control%v_grid(halo_type_extended)%get_num_rows()
num_levels  = output_v_field%get_num_levels()
halo_ns     = lbc_output_control%v_grid(halo_type_extended)%get_halo_ns()
halo_ew     = lbc_output_control%v_grid(halo_type_extended)%get_halo_ew()
rim_width   = lbc_output_control%rim_width

! Get LBC level size of V field grid with extended haloes
lbc_v_level_size = SIZE(interp_weights(weights_index(v_points), &
                                       halo_type_extended)%bilinear_index_b_l)
ALLOCATE(temp_field_data(lbc_v_level_size, num_levels))
temp_field_data = rmdi
! Calculate latitude linear weights
ALLOCATE(weight1(v_grid_rows))
ALLOCATE(weight2(v_grid_rows))
DO i = 1, v_grid_rows
  weight2(i) = (v_grid_lat(i) - p_grid_lat(i)) / &
       (p_grid_lat(i+1) - p_grid_lat(i))
  weight1(i) = 1 - weight2(i)
END DO
! Each level of data is contained in a 1D array, indexed following the LBC structure
! described in C71. As the enlarged P grid and the V grid have a different LBC level
! size, care must be taken to make sure the correct LBC data points are accessed
! when performing the linear interpolation. Each LBC section will have more P points
! than V points.

! Linearly interpolate each LBC section (N,E,S,W) seperately.

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(k, i, j, counter_v_lbc, counter_p_lbc,         &
!$OMP north_offset, east_offset, south_offset, west_offset)                         &
!$OMP SHARED(num_levels, v_grid_cols, v_grid_rows, rim_width, halo_ns, halo_ew,     &
!$OMP weight1, weight2, temp_field_data, p_grid_cols, output_v_field,               &
!$OMP lbc_output_control)
IF (lbc_output_control%output_grid_stagger == new_dynamics) THEN
!$OMP DO
  DO k = 1, num_levels
    counter_v_lbc = 0
    counter_p_lbc = 0
    ! North LBC section:
    !  The enlarged P grid has extra row of points on the bottom of the north section.
    !  This is to ensure that the V points in the north section are surrounded by P
    !  points
    !
    !  Offset to determine the lbc index of the P point directly north of the V point
    !  we wish to interpolate to. For the north section this is the number of P points in a
    !  row
    north_offset = p_grid_cols
    DO j = 1, rim_width + halo_ns
      DO i = 1, v_grid_cols
        counter_v_lbc = counter_v_lbc + 1
        counter_p_lbc = counter_p_lbc + 1
        temp_field_data(counter_v_lbc, k) = (weight1(v_grid_rows-rim_width-halo_ns+j) *      &
                                             output_v_field%lbc_rdata(counter_p_lbc, k)) +   &
                                            (weight2(v_grid_rows-rim_width-halo_ns+j) *      &
                                             output_v_field%lbc_rdata(counter_p_lbc+north_offset, k))
      END DO
    END DO
    ! Extra row of P points compared to V points
    counter_p_lbc = counter_p_lbc + p_grid_cols

    ! East LBC section:
    !  The enlarged P grid has extra column of points on the left of the east section.
    !
    !  Offset to determine the lbc index of the P point directly north of the V point
    !  we wish to interpolate to. For the east section this is the east-west halo size + rim width
    !  + 1 to account for the extra column of P points
    east_offset = rim_width + halo_ew + 1
    DO j = 1, v_grid_rows - (2*rim_width) - (2*halo_ns)
      ! East section has extra column of P points on its left side, need to skip these as they are used
      ! when interpolating U field and are not needed for the V field
      counter_p_lbc = counter_p_lbc + 1
      DO i = 1, rim_width + halo_ew
        counter_v_lbc = counter_v_lbc + 1
        counter_p_lbc = counter_p_lbc + 1
        temp_field_data(counter_v_lbc, k) = (weight1(rim_width+halo_ns-1+j) *                &
             output_v_field%lbc_rdata(counter_p_lbc, k)) + (weight2(rim_width+halo_ns-1+j) * &
             output_v_field%lbc_rdata(counter_p_lbc+east_offset, k))
      END DO
    END DO
    ! The east section P grid has an extra row compared to the V grid. Need to skip these
    ! before we can move onto the south section
    counter_p_lbc = counter_p_lbc + rim_width + halo_ew + 1
    ! South LBC section:
    !  The enlarged P grid has extra row of points on the top of the south section.
    !
    !  Offset to determine the lbc index of the P point directly north of the V point
    !  we wish to interpolate to. For the south section this is the the number of P points in a
    !  row
    south_offset = p_grid_cols
    DO j = 1, rim_width + halo_ns
      DO i = 1, v_grid_cols
        counter_v_lbc = counter_v_lbc + 1
        counter_p_lbc = counter_p_lbc + 1
        temp_field_data(counter_v_lbc, k) = (weight1(j) *                                    &
                                             output_v_field%lbc_rdata(counter_p_lbc, k)) +   & 
                                            (weight2(j) *                                    &
                                             output_v_field%lbc_rdata(counter_p_lbc+south_offset, k))
      END DO
    END DO
    ! The enlarged P grid has extra row of points on the top of the south section.
    ! We need to skip this row before moving onto the west section.
    counter_p_lbc = counter_p_lbc + p_grid_cols
    ! West LBC section:
    !  The enlarged P grid has extra column of points on the left of the west section.
    !
    !  Offset to determine the lbc index of the P point directly north of the V point
    !  we wish to interpolate to. For the west section this is rim width + east west halo
    !  size + 1 to account for the extra column of P points.
    west_offset = rim_width + halo_ew + 1
    DO j = 1, v_grid_rows - (2*rim_width) - (2*halo_ns)
      DO i = 1, rim_width + halo_ew
        counter_v_lbc = counter_v_lbc + 1
        counter_p_lbc = counter_p_lbc + 1
        temp_field_data(counter_v_lbc, k) = (weight1(rim_width+halo_ns-1+j) *                &
                                             output_v_field%lbc_rdata(counter_p_lbc, k)) +   &
                                            (weight2(rim_width+halo_ns-1+j) *                &
                                            output_v_field%lbc_rdata(counter_p_lbc+west_offset, k))
      END DO
      ! West section has extra column of P points on its right side, need to skip these as they are used
      ! when interpolating U field and are not needed for the V field
      counter_p_lbc = counter_p_lbc + 1
    END DO
  END DO
!$OMP END DO
ELSE IF (lbc_output_control%output_grid_stagger == endgame) THEN
  ! The EG enlarged P grid object an extra row on north and south and 
  ! extra column on west as such the num_rows is increased by 2 and num
  ! cols increased by 1 compared to the standard P grid.
!$OMP DO    
  DO k = 1, num_levels
    counter_v_lbc = 0
    counter_p_lbc = 0
    ! North LBC section:
    !  The enlarged P grid has extra row of points on the top of the north section and
    !  an extra column on left of north section.
    !  Offset to determine the lbc index of the P point directly north of the V point
    !  we wish to interpolate to. For the north section this is the number of P points in a
    !  row
    north_offset = p_grid_cols
    DO j = 1, rim_width + halo_ns
      ! North section has extra column of P points on its left side, need to skip these.
      counter_p_lbc = counter_p_lbc + 1
      DO i = 1, v_grid_cols
        counter_v_lbc = counter_v_lbc + 1
        counter_p_lbc = counter_p_lbc + 1
        temp_field_data(counter_v_lbc, k) = (weight1(v_grid_rows-rim_width-halo_ns+j) *      &
                                             output_v_field%lbc_rdata(counter_p_lbc, k)) +   &
                                            (weight2(v_grid_rows-rim_width-halo_ns+j) *      &
                                             output_v_field%lbc_rdata(counter_p_lbc+north_offset, k))
      END DO
    END DO
    ! Extra row of P points compared to V points
    counter_p_lbc = counter_p_lbc + p_grid_cols
    ! East LBC section:
    !  The enlarged LBC grid east section has extra row at top and bottom and 
    !  and extra column on the left.
    !  Offset to determine the lbc index of the P point directly north of the V point
    !  we wish to interpolate to. For the east section this is the east-west halo size + rim width
    !  + 1 to account for the extra column of P points
    east_offset = rim_width + halo_ew + 1
    DO j = 1, v_grid_rows - (2*rim_width) - (2*halo_ns)
      ! East section has extra column of P points on its left side, need to skip these as they are used
      ! when interpolating U field and are not needed for the V field
      counter_p_lbc = counter_p_lbc + 1
      DO i = 1, rim_width + halo_ew
        counter_v_lbc = counter_v_lbc + 1
        counter_p_lbc = counter_p_lbc + 1
        temp_field_data(counter_v_lbc, k) = (weight1(rim_width+halo_ns-1+j) *              &
                                             output_v_field%lbc_rdata(counter_p_lbc, k)) + &
                                             (weight2(rim_width+halo_ns-1+j) *             &
                                             output_v_field%lbc_rdata(counter_p_lbc+east_offset, k))
      END DO
    END DO
    ! The east section P grid has an extra row compared to the V grid. Need to skip these
    ! before we can move onto the south section
    counter_p_lbc = counter_p_lbc + rim_width + halo_ew + 1
    ! South LBC section:
    !  The enlarged P grid south section has an extra row of points on bottom edge and extra
    !  column on left side
    !  Offset to determine the lbc index of the P point directly north of the V point
    !  we wish to interpolate to. For the south section this is the the number of P points in a
    !  row
    south_offset = p_grid_cols
    DO j = 1, rim_width + halo_ns
      ! South section has extra column of P points on its left side, need to skip these as they are used
      ! when interpolating U field and are not needed for the V field
      counter_p_lbc = counter_p_lbc + 1
      DO i = 1, v_grid_cols
        counter_v_lbc = counter_v_lbc + 1
        counter_p_lbc = counter_p_lbc + 1
        temp_field_data(counter_v_lbc, k) = (weight1(j) *                                  &
                                             output_v_field%lbc_rdata(counter_p_lbc, k)) + & 
                                            (weight2(j) *                                  &
                                             output_v_field%lbc_rdata(counter_p_lbc+south_offset, k))
      END DO
    END DO
     ! The east section P grid has an extra row compared to the V grid. Need to skip these
    ! before we can move onto the south section
    counter_p_lbc = counter_p_lbc + p_grid_cols
    ! West LBC section:
    !  The enlarged LBC grid west section has extra row at top and bottom and 
    !  and extra column on the left.
    !  Offset to determine the lbc index of the P point directly north of the V point
    !  we wish to interpolate to. For the west section this is rim width + east west halo
    !  size + 1 to account for the extra column of P points.
    west_offset = rim_width + halo_ew + 1
    DO j = 1, v_grid_rows - (2*rim_width) - (2*halo_ns)
      ! West section has extra column of P points on its left side, need to skip these as they are used
      ! when interpolating U field and are not needed for the V field
      counter_p_lbc = counter_p_lbc + 1
      DO i = 1, rim_width + halo_ew
        counter_v_lbc = counter_v_lbc + 1
        counter_p_lbc = counter_p_lbc + 1
        temp_field_data(counter_v_lbc, k) = (weight1(rim_width+halo_ns-1+j) *              &
                                             output_v_field%lbc_rdata(counter_p_lbc, k)) + &
                                            (weight2(rim_width+halo_ns-1+j) *              &
                                             output_v_field%lbc_rdata(counter_p_lbc+west_offset, k))
      END DO
    END DO
  END DO ! End level loop
!$OMP END DO
END IF
!$OMP END PARALLEL
DEALLOCATE(weight1)
DEALLOCATE(weight2)

! Replace P grid data with newly interpolated V grid values
CALL MOVE_ALLOC(temp_field_data, output_v_field%lbc_rdata)
! Set field latitudes and longitudes to V grid values
CALL output_v_field%set_latitudes(lbc_output_control%v_grid(halo_type_extended)%get_latitudes())
CALL output_v_field%set_longitudes(lbc_output_control%v_grid(halo_type_extended)%get_longitudes())
CALL output_v_field%set_num_rows(v_grid_rows)
CALL output_v_field%set_num_cols(v_grid_cols)
CALL output_v_field%set_horiz_grid_type(v_points)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE interp_output_wind_p_to_v

END MODULE interp_output_winds_mod
