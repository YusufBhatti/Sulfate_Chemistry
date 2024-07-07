! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE process_orography_mod
!
! Description:
!   Interpolate the input orography from the P grid onto the LBC
!   P, U and V grids and save in an array.  Do this for each size of LBC halos.  
!   These LBC orography fields can then be used when vertically interpolating
!   multi-level fields.
!   
!   The orography LBC field which will be output to the LBC file will
!   not be taken from this array of fields, it will just be processed
!   along with the other fields.  Interpolating a single level field such
!   as orography is not computationally expensive and this is preferred
!   rather than introduce code to copy a field from the orography array.
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

! LBC orography field array indices
INTEGER, PARAMETER :: u_points_orog = 1
INTEGER, PARAMETER :: v_points_orog = 2
INTEGER, PARAMETER :: p_points_orog = 3

! Weight arrays indices
INTEGER, PARAMETER :: p_to_u_points = 4
INTEGER, PARAMETER :: p_to_v_points = 5

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'PROCESS_OROGRAPHY_MOD'

CONTAINS

!-------------------------------------------------------------------------------
                                             
SUBROUTINE process_orography(input_file, lbc_output_control,     &
                             orography_lbc_fields,               &
                             interp_weights,                     &
                             input_orog_field_number)


USE datafile_mod,              ONLY: datafile_type
USE lbc_output_control_mod,    ONLY: lbc_output_control_type
USE field_mod,                 ONLY: field_type, ASSIGNMENT(=)
USE three_dimensional_grid_mod,ONLY: three_dimensional_grid_type
USE interp_weights_mod,        ONLY: interp_weights_type, weights_index
USE stashmaster_constants_mod, ONLY: u_points, v_points, p_points
USE interp_lbc_mod,       ONLY: interp_lbc

IMPLICIT NONE

CLASS(datafile_type), INTENT(INOUT)        :: input_file
TYPE(lbc_output_control_type), INTENT(IN) :: lbc_output_control
TYPE(field_type), INTENT(INOUT)            :: orography_lbc_fields(3,3)
CLASS(interp_weights_type), INTENT(IN)     :: interp_weights(:,:)
INTEGER, INTENT(IN)                        :: input_orog_field_number

INTEGER :: halo_code
LOGICAL :: field_vert_interp = .FALSE. ! Vertical interpolation not needed for orography
CHARACTER(LEN=*), PARAMETER :: routinename = 'PROCESS_OROGRAPHY'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO halo_code=1,3
  ! Interpolate orography onto U grid
  orography_lbc_fields(u_points_orog, halo_code) =                                                 &
                                      ! Call interp_lbc function to return interpolated field
                                      interp_lbc(input_file%fields(input_orog_field_number),       &
                                      lbc_output_control%u_grid(halo_code),                        &
                                      lbc_output_control%horizontal_interpolation_method,          &
                                      lbc_output_control%vertical_interpolation_method,            &
                                      interp_weights(p_to_u_points, halo_code), field_vert_interp, &
                                      lbc_output_control%output_grid_stagger)
  ! Interpolate orography onto V grid
  orography_lbc_fields(v_points_orog, halo_code) =                                                 &
                                      interp_lbc(input_file%fields(input_orog_field_number),       &
                                      lbc_output_control%v_grid(halo_code),                        &
                                      lbc_output_control%horizontal_interpolation_method,          &
                                      lbc_output_control%vertical_interpolation_method,            &
                                      interp_weights(p_to_v_points, halo_code), field_vert_interp, &
                                      lbc_output_control%output_grid_stagger)
  ! Interpolate orography onto P grid 
  orography_lbc_fields(p_points_orog, halo_code) =                                                 &
                                      interp_lbc(input_file%fields(input_orog_field_number),       &
                                      lbc_output_control%p_grid(halo_code),                        &
                                      lbc_output_control%horizontal_interpolation_method,          &
                                      lbc_output_control%vertical_interpolation_method,            &
                                      interp_weights(weights_index(p_points), halo_code),          &
                                      field_vert_interp,                                           &
                                      lbc_output_control%output_grid_stagger)
END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE process_orography

!-------------------------------------------------------------------------------

SUBROUTINE process_orography_enlarged_grid(input_file, lbc_output_control, &
                             orography_enlarged_grid_lbc_field,            &
                             rotation_interp_weights,                      &
                             input_orog_field_number)

USE datafile_mod,              ONLY: datafile_type
USE lbc_output_control_mod,    ONLY: lbc_output_control_type
USE field_mod,                 ONLY: field_type, ASSIGNMENT(=)
USE interp_weights_mod,        ONLY: interp_weights_type, p_to_enlarged_p
USE umPrintMgr,                ONLY: umPrint, umMessage
USE interp_lbc_mod,       ONLY: interp_lbc

IMPLICIT NONE

CLASS(datafile_type), INTENT(INOUT)        :: input_file
TYPE(lbc_output_control_type), INTENT(IN) :: lbc_output_control
TYPE(field_type), INTENT(INOUT)            :: orography_enlarged_grid_lbc_field
CLASS(interp_weights_type), INTENT(IN)     :: rotation_interp_weights(:)
INTEGER, INTENT(IN)                        :: input_orog_field_number

LOGICAL :: field_vert_interp = .FALSE. ! Vertical interpolation not needed for orography
CHARACTER(LEN=*), PARAMETER :: routinename = 'PROCESS_OROGRAPHY_ENLARGED_GRID'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! If wind rotations are being performed on the output lbc fields then the U and V winds will
! be on an enlarged P grid therefore need to calculate orography on the enlarged P grid.
IF (.NOT. lbc_output_control%l_same_rotation  &
     .AND. lbc_output_control%l_target_rotated) THEN
  ! Calc enlarged P grid orography
  WRITE(umMessage, '(A)') "[INFO] Calculating orography for the enlarged LBC p grid"
  CALL umPrint(umMessage, src='process_orography_enlarged_grid')
  orography_enlarged_grid_lbc_field =                          &
       interp_lbc( input_file%fields(input_orog_field_number), &
       lbc_output_control%p_grid_enlarged,                     &
       lbc_output_control%horizontal_interpolation_method,     &
       lbc_output_control%vertical_interpolation_method,       &
       rotation_interp_weights(p_to_enlarged_p),               &
       field_vert_interp, lbc_output_control%output_grid_stagger)
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE process_orography_enlarged_grid

!-------------------------------------------------------------------------------

END MODULE process_orography_mod
