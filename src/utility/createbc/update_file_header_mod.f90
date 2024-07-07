! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE update_file_header_mod
! Description:
!   A routine to set header values based on the LBC definition
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

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UPDATE_FILE_HEADER_MOD'

CONTAINS

SUBROUTINE update_file_header(output_file, lbc_output_control)

USE fieldsfile_mod,          ONLY: fieldsfile_type
USE lbcfile_mod,             ONLY: lbcfile_type
USE datafile_mod,            ONLY: datafile_type
USE lbc_output_control_mod,  ONLY: lbc_output_control_type
USE um_version_mod,          ONLY: um_version_int
USE um_parparams,            ONLY: halo_type_no_halo
USE missing_data_mod,        ONLY: imdi

USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
  
IMPLICIT NONE

CLASS(datafile_type), INTENT(INOUT) :: output_file
CLASS(lbc_output_control_type), INTENT(IN) :: lbc_output_control

CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER :: routinename='UPDATE_FILE_HEADER'
INTEGER :: icode
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

SELECT TYPE(output_file)
CLASS IS (fieldsfile_type)
  output_file%source_model_version = um_version_int
  IF (output_file%p_grid%get_grid_code() == 101) THEN
    output_file%horizontal_grid_indicator = 103
  ELSE IF (output_file%p_grid%get_grid_code() == 1) THEN
    output_file%horizontal_grid_indicator = 3
  ELSE
    output_file%horizontal_grid_indicator = imdi
  END IF
  IF (lbc_output_control%num_reserved_headers /= imdi) THEN
    output_file%num_reserved_headers = lbc_output_control%num_reserved_headers
  END IF

CLASS IS (lbcfile_type)
  output_file%source_model_version = um_version_int
  ! Copy the vertical grid with the new eta theta and eta rho values
  output_file%file_theta_rho_levels = lbc_output_control%p_grid(1)%vert_grid
  ! CreateBC only supports the "Smooth" method of height generation
  output_file%algorithm_to_generate_height_fields = 2
  ! Check grid staggering
  output_file%grid_staggering = lbc_output_control%output_grid_stagger
  ! Copy the U, V and P grids to the output file
  output_file%p_grid = lbc_output_control%p_grid(halo_type_no_halo)
  output_file%u_grid = lbc_output_control%u_grid(halo_type_no_halo)
  output_file%v_grid = lbc_output_control%v_grid(halo_type_no_halo)
  ! Update variable resolution flag
  output_file%variable_resolution = lbc_output_control%variable_resolution
  ! Set number of levels
  output_file%num_model_levels = lbc_output_control%num_levels
  ! Pole lat / long
  output_file%pole_lat = lbc_output_control%pole_lat
  output_file%pole_long = lbc_output_control%pole_long
  ! LBC files only have 4 sets of level dependant constants
  output_file%len2_lev_dep_constants = 4
  ! Overide number of reserved headers with value from namelist
  IF (lbc_output_control%num_reserved_headers /= imdi) THEN
    output_file%num_reserved_headers = lbc_output_control%num_reserved_headers
  END IF
  IF (output_file%p_grid%get_grid_code() == 101) THEN
    output_file%horizontal_grid_indicator = 103
  ELSE IF (output_file%p_grid%get_grid_code() == 1) THEN
    output_file%horizontal_grid_indicator = 3
  ELSE
    output_file%horizontal_grid_indicator = imdi
  END IF

CLASS DEFAULT
  cmessage = 'Unknown type of output file, unable to update header'
  icode = 10
  CALL ereport(routinename, icode, cmessage)
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE update_file_header

END MODULE update_file_header_mod
