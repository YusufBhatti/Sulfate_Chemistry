! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE calc_lbc_coords_mod

! DrHook Modules
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

! Description:
!  Routine to calculate the size of an LBC field and calculate latitude and longitude
!  of each lbc point.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'CALC_LBC_COORDS_MOD'

CONTAINS

SUBROUTINE calc_lbc_coords(enlarged_p_grid, lbc_output_control, lbc_grid,    &
                           latitude_lbc, longitude_lbc, lbc_size)

USE lbc_output_control_mod,     ONLY: lbc_output_control_type, new_dynamics, &
                                      endgame
USE three_dimensional_grid_mod, ONLY: three_dimensional_grid_type
USE ereport_mod,                ONLY: ereport
USE errormessagelength_mod,     ONLY: errormessagelength

IMPLICIT NONE

! Arguments
LOGICAL, INTENT(IN) :: enlarged_p_grid
TYPE(lbc_output_control_type), INTENT(IN)    :: lbc_output_control
TYPE(three_dimensional_grid_type), INTENT(IN) :: lbc_grid
REAL, ALLOCATABLE, INTENT(OUT) :: latitude_lbc(:)    ! Latitude points for LBC grid
REAL, ALLOCATABLE, INTENT(OUT) :: longitude_lbc(:)   ! Longitude points for LBC grid
INTEGER, INTENT(OUT) :: lbc_size ! Size of single level of lbc field

! Local variables
REAL, ALLOCATABLE :: latitude_target(:)  ! Lat points from LBC grid definition namelist
REAL, ALLOCATABLE :: longitude_target(:) ! Longitude points from LBC grid definition namelist
INTEGER           :: lbc_row_len         ! Row length of LBC grid 
INTEGER           :: lbc_rows            ! Rows of LBC grid       
INTEGER           :: row, col, pass_number
INTEGER           :: lbc_index
! Offsets needed for the enlarged P grid
INTEGER :: extra_nt = 0
INTEGER :: extra_nb = 0
INTEGER :: extra_nl = 0
INTEGER :: extra_et = 0
INTEGER :: extra_eb = 0
INTEGER :: extra_el = 0
INTEGER :: extra_sl = 0
INTEGER :: extra_sr = 0
INTEGER :: extra_st = 0
INTEGER :: extra_sb = 0
INTEGER :: extra_wl = 0
INTEGER :: extra_wr = 0
INTEGER :: extra_wt = 0
INTEGER :: extra_wb = 0
! Loop bounds
INTEGER :: north_row_start_bound 
INTEGER :: north_row_end_bound   
INTEGER :: north_col_start_bound 
INTEGER :: north_col_end_bound   
INTEGER :: east_row_start_bound  
INTEGER :: east_row_end_bound    
INTEGER :: east_col_start_bound  
INTEGER :: east_col_end_bound    
INTEGER :: south_row_start_bound 
INTEGER :: south_row_end_bound   
INTEGER :: south_col_start_bound 
INTEGER :: south_col_end_bound   
INTEGER :: west_row_start_bound  
INTEGER :: west_row_end_bound    
INTEGER :: west_col_start_bound  
INTEGER :: west_col_end_bound    


! Error reporting
INTEGER :: icode
CHARACTER(LEN=*), PARAMETER :: routinename = 'CALC_LBC_COORDS'
CHARACTER(LEN=errormessagelength) :: cmessage
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! EG and ND need different enlarged grids, due to their different grid staggering.  The U and V
! points must be surrounded by P points.  Each LBC section (north, south, east, west) will be
! separately enlarged as needed, even if this results in some duplication (i.e. if you increase 
! the top of the south section it will overlap with the bottom of the east and west sections).
! This duplication will have a trivial computational cost but will make the linear interpolation
! from the enlarged P grid back to the standard U and V grid less complicated, as this will be done 
! on each section separately
IF (enlarged_p_grid) THEN
  ! After interpolation both U and V need to be on the same enlarged P grid. The 
  ! enlarged P grid must therefore surround both U and V grids.
  IF (lbc_output_control%output_grid_stagger == endgame) THEN
    ! ENDGame grid. In order to surround:
    ! i)  U points - need extra point on left side of all four sections
    ! ii) V points - need extra point on top and bottom of east and westion sections,
    !                and extra point on top of north section and bottom of south
    extra_nt = 1
    extra_nb = 0
    extra_nl = 1
    extra_et = 1
    extra_eb = 1
    extra_el = 1
    extra_st = 0
    extra_sb = 1
    extra_sl = 1
    extra_sr = 0
    extra_wt = 1
    extra_wb = 1
    extra_wl = 1
    extra_wr = 0
  ELSE
    ! New Dynamics grid. In order to surround:
    ! i)  U points - need extra point on right side of west section and on
    !                the left of east section.
    ! ii) V points - need extra point on top of south section and bottom of
    !                north section.
    extra_nt = 0
    extra_nb = 1
    extra_nl = 0
    extra_et = 0
    extra_eb = 0
    extra_el = 1
    extra_sl = 0
    extra_sr = 0
    extra_st = 1
    extra_sb = 0
    extra_wl = 0
    extra_wr = 1
    extra_wt = 0
    extra_wb = 0
  END IF
END IF

! Get the target domain grid information from the lbc grid object
! lbc_rows and row_len are for the full LBC grid inc. halos
lbc_row_len = lbc_grid%get_num_cols()
lbc_rows    = lbc_grid%get_num_rows()
CALL lbc_grid%get_horizontal_coordinates(latitude_target, longitude_target)

! Calculate the LBC grid lat and longitude values from the input namelist.
! This will be stored in a 1D array and filled in the order expected by
! the LBC files, ie. the four boundary areas: north, east, west and south.
! See UMDP C71 for more details on the LBC file definition.
!
! Note that the LBC grid lat and long arrays are indexed from 1 to a num_rows/num_cols
! that includes the halos.
!
! The enlarged LBC P grid for EG has an extra row on north and south and and extra column 
! on the west.  The LBC P grid for ND is just the standard P grid as P points surround
! all U and V points.

! Loop bounds
! North
north_row_start_bound = lbc_rows - lbc_grid%get_halo_ns() -            &
                        lbc_output_control%rim_width - extra_nb - extra_nt + 1 
north_row_end_bound   = lbc_rows
north_col_start_bound = 1
north_col_end_bound   = lbc_row_len
! East
east_row_start_bound  = 1 + lbc_grid%get_halo_ns() +                   &
                        lbc_output_control%rim_width
east_row_end_bound    = lbc_rows - lbc_grid%get_halo_ns() -            &
                        lbc_output_control%rim_width
east_col_start_bound  = lbc_row_len - lbc_grid%get_halo_ew() -         &
                        lbc_output_control%rim_width - extra_el + 1
east_col_end_bound    = lbc_row_len
! South
south_row_start_bound = 1
south_row_end_bound   = lbc_grid%get_halo_ns() +                       &
                        lbc_output_control%rim_width + extra_st + extra_sb 
south_col_start_bound = 1
south_col_end_bound   = lbc_row_len
! West
west_row_start_bound  = 1 + lbc_grid%get_halo_ns() +                   &
                        lbc_output_control%rim_width
west_row_end_bound    = lbc_rows - lbc_grid%get_halo_ns() -            &
                        lbc_output_control%rim_width
west_col_start_bound  = 1
west_col_end_bound    = lbc_grid%get_halo_ew() +                       &
                        lbc_output_control%rim_width + extra_wr + extra_wl

! First work through each LBC section to determine the LBC size
lbc_size = 0
lbc_size = lbc_size + ( north_row_end_bound - north_row_start_bound + 1 ) * &
                      ( north_col_end_bound - north_col_start_bound + 1 )
lbc_size = lbc_size + ( east_row_end_bound  - east_row_start_bound  + 1 ) * &
                      ( east_col_end_bound  - east_col_start_bound  + 1 )
lbc_size = lbc_size + ( south_row_end_bound - south_row_start_bound + 1 ) * &
                      ( south_col_end_bound - south_col_start_bound + 1 )
lbc_size = lbc_size + ( west_row_end_bound  - west_row_start_bound  + 1 ) * &
                      ( west_col_end_bound  - west_col_start_bound  + 1 )

ALLOCATE (longitude_lbc(lbc_size))
ALLOCATE (latitude_lbc(lbc_size))
lbc_index = 0

! Now that the longitude and latitude arrays have been allocated they can be populated

!  Northern Boundary 
DO row = north_row_start_bound, north_row_end_bound
  DO col = north_col_start_bound, north_col_end_bound
    lbc_index = lbc_index + 1
    longitude_lbc(lbc_index) = longitude_target(col)
    latitude_lbc(lbc_index)  = latitude_target(row)
  END DO
END DO

!  Eastern Boundary
DO row = east_row_start_bound, east_row_end_bound
  DO col = east_col_start_bound, east_col_end_bound
    lbc_index = lbc_index + 1
    longitude_lbc(lbc_index) = longitude_target(col)
    latitude_lbc(lbc_index)  = latitude_target(row)
  END DO
END DO

!  Southern Boundary 
DO row = south_row_start_bound, south_row_end_bound
  DO col = south_col_start_bound, south_col_end_bound
    lbc_index = lbc_index + 1
    longitude_lbc(lbc_index) = longitude_target(col)
    latitude_lbc(lbc_index)  = latitude_target(row)
  END DO
END DO

!  Western Boundary
DO row = west_row_start_bound, west_row_end_bound
  DO col = west_col_start_bound, west_col_end_bound
    lbc_index = lbc_index + 1
    longitude_lbc(lbc_index) = longitude_target(col)
    latitude_lbc(lbc_index)  = latitude_target(row)
  END DO
END DO

IF ( lbc_index /= lbc_size ) THEN
  WRITE(cmessage, '(A,I0,A,I0)') "Error calculating size of lbc field. lbc_index = ", &
       lbc_index, " lbc_size = ", lbc_size
  icode = 10
  CALL ereport(routinename, icode, cmessage)
END IF
  
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE calc_lbc_coords

END MODULE calc_lbc_coords_mod

