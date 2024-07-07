! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE update_frame_field_grid_mod
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
  ModuleName = 'UPDATE_FRAME_FIELD_GRID_MOD'

CONTAINS

SUBROUTINE update_frame_field_grid(output_file, field_num)
!
!  Description: Update the vertical and horizontal grids
!   contained in the output field object with the appropiate frame 
!   grids from the output file
!   
USE datafile_mod,                 ONLY: datafile_type
USE horizontal_lat_long_grid_mod, ONLY: horizontal_lat_long_grid_type
USE stashmaster_constants_mod,    ONLY: u_points, v_points, p_points
USE ereport_mod,                  ONLY: ereport
USE errormessagelength_mod,       ONLY: errormessagelength

IMPLICIT NONE
! Args
CLASS(datafile_type), TARGET, INTENT(INOUT) :: output_file
INTEGER,              INTENT(IN)    :: field_num
! Local variables
TYPE(horizontal_lat_long_grid_type), POINTER :: field_horiz_grid
TYPE(horizontal_lat_long_grid_type), POINTER :: frames_horiz_grid
INTEGER :: num_rows, num_cols
CHARACTER(LEN=*), PARAMETER :: routinename = 'UPDATE_FRAME_FIELD_GRID'
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: cmessage

! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

field_horiz_grid => output_file%fields(field_num)%grid%horiz_grid

SELECT CASE (field_horiz_grid%get_grid_type())
CASE (u_points)
  frames_horiz_grid => output_file%u_grid%horiz_grid
CASE (v_points)
  frames_horiz_grid => output_file%v_grid%horiz_grid
CASE (p_points)
  frames_horiz_grid => output_file%p_grid%horiz_grid
CASE DEFAULT
  icode = 100
  WRITE(cmessage, '(A,I8)') "Unsupported grid type: ",field_horiz_grid%get_grid_type()
  CALL ereport(routinename, icode, cmessage)
END SELECT

! Update grid code
CALL field_horiz_grid%set_grid_code(frames_horiz_grid%get_grid_code())
! Update hemisphere indicator
CALL field_horiz_grid%set_hemisphere_indicator(frames_horiz_grid%get_hemisphere_indicator())
! Update num rows and cols
CALL field_horiz_grid%set_num_rows(frames_horiz_grid%get_num_rows())
CALL field_horiz_grid%set_num_cols(frames_horiz_grid%get_num_cols())
! Update lat and long arrays
CALL field_horiz_grid%set_latitudes(frames_horiz_grid%get_latitudes())
CALL field_horiz_grid%set_longitudes(frames_horiz_grid%get_longitudes())

NULLIFY(field_horiz_grid)
NULLIFY(frames_horiz_grid)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE update_frame_field_grid

END MODULE update_frame_field_grid_mod
