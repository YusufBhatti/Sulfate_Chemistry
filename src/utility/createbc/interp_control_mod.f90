! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE interp_control_mod
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.
USE field_mod, ONLY: field_type, ASSIGNMENT(=)
! DrHook Modules
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE
PRIVATE
PUBLIC :: interp_control

CHARACTER(LEN=*), PARAMETER :: ModuleName = 'INTERP_CONTROL_MOD'

CONTAINS

TYPE(field_type) FUNCTION interp_control(lbc_output_control,                 &
                 input_file, interp_weights, orography_lbc_fields,           &
                 l_vert_interp, field_number, field_in)    RESULT(output_field)
! Description: A wrapper routine to interp_lbc to determine which  
!  input/output grids and orography are appropiate for the current field. These
!  can then be passed into the argument list of interp_lbc
!
USE lbc_output_control_mod,    ONLY: lbc_output_control_type, create_frame
USE datafile_mod,              ONLY: datafile_type
USE three_dimensional_grid_mod,ONLY: three_dimensional_grid_type
USE stashmaster_constants_mod, ONLY: u_points, v_points, p_points,           &
                                     rho_levels, theta_levels, single_level, &
                                     halo_type_code
USE fieldsfile_mod,            ONLY: fieldsfile_type
USE interp_weights_mod,        ONLY: interp_weights_type, weights_index
USE ereport_mod,               ONLY: ereport
USE errormessagelength_mod,    ONLY: errormessagelength
USE stashmaster_utils_mod,     ONLY: query_stashmaster
USE process_orography_mod,     ONLY: u_points_orog, v_points_orog, p_points_orog
USE interp_lbc_mod,            ONLY: interp_lbc
USE lbc_stashcode_mapping_mod, ONLY: lbc_stashcode_mapping

IMPLICIT NONE

CLASS(datafile_type), TARGET, INTENT(INOUT) :: input_file
CLASS(lbc_output_control_type), INTENT(IN) :: lbc_output_control
CLASS(interp_weights_type), INTENT(IN)      :: interp_weights(:,:)
CLASS(field_type), TARGET, INTENT(IN)  :: orography_lbc_fields(:,:)
LOGICAL, INTENT(IN) :: l_vert_interp
INTEGER, OPTIONAL, INTENT(IN)               :: field_number 
CLASS(field_type), OPTIONAL, TARGET, INTENT(INOUT)  :: field_in

TYPE(three_dimensional_grid_type), SAVE :: source_grid, target_grid
CLASS(field_type), POINTER :: input_field
CLASS(field_type), POINTER :: interp_orog
LOGICAL :: field_vert_interp

CHARACTER(LEN=*), PARAMETER :: routinename = 'INTERP_CONTROL'
INTEGER :: icode
INTEGER :: lbc_halo_code! Halo code of output field
INTEGER :: stash_out    ! Temp variable to hold output STASHcode
CHARACTER(LEN=errormessagelength) :: cmessage
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Get the input field which is going to be interpolated. This can either be
! taken directly from a field object passed in from the argument list, or if
! the field object is not present the field number will be used to get the
! input field from the input file.
IF (PRESENT(field_in)) THEN
  input_field => field_in
ELSE IF (PRESENT(field_number)) THEN
  input_field => input_file%fields(field_number)
ELSE
  icode = 60
  WRITE(cmessage, '(A)') "Either field number of field object must be passed to routine"
  CALL ereport(routinename, icode, cmessage)
END IF

! Set the field vertical interpolation flag to the global value determined by the
! check_vert_interp routine
field_vert_interp = l_vert_interp

! Use input field STASHcode to determine output LBC stashcode
CALL lbc_stashcode_mapping(input_field%quantity_ident, stash_out)
! Set halo code based on new STASHcode
lbc_halo_code = query_stashmaster(stash_out, halo_type_code)

! LBCs have 9 possible target grids, three halo types x three grid types

SELECT CASE (input_field%get_horiz_grid_type())
CASE (u_points)
  ! Set interpolated orography to the appropiate grid type
  interp_orog => orography_lbc_fields(u_points_orog, lbc_halo_code)
  target_grid = lbc_output_control%u_grid(lbc_halo_code)    
  IF (input_field%grid%vert_grid%level_type == single_level) THEN
    field_vert_interp = .FALSE.
  END IF
CASE (v_points)
  ! Set interpolated orography to the appropiate grid type
  interp_orog => orography_lbc_fields(v_points_orog, lbc_halo_code)
  target_grid = lbc_output_control%v_grid(lbc_halo_code)
  IF (input_field%grid%vert_grid%level_type == single_level) THEN
    field_vert_interp = .FALSE.
  END IF
CASE (p_points)
  ! Set interpolated orography to the appropiate grid type
  interp_orog => orography_lbc_fields(p_points_orog, lbc_halo_code)
  target_grid = lbc_output_control%p_grid(lbc_halo_code)
  IF (input_field%grid%vert_grid%level_type == single_level) THEN
    field_vert_interp = .FALSE.
  END IF
CASE DEFAULT
  icode = 100
  WRITE(cmessage, '(A,I8)') "Unsupported grid type: ",input_field%get_horiz_grid_type()
  CALL ereport(routinename, icode, cmessage)

END SELECT

CALL timer( 'interp_lbc', 5)

output_field = interp_lbc(                                                &
     input_field, target_grid,                                                 &
     lbc_output_control%horizontal_interpolation_method,                      &
     lbc_output_control%vertical_interpolation_method,                        &
     ! Use grid type and halo code to select the correct interpolation weights
     interp_weights(                                                           &
     weights_index(input_field%get_horiz_grid_type()), lbc_halo_code),     &
     field_vert_interp,                                                        &
     lbc_output_control%output_grid_stagger,                                  &
     interp_orog=interp_orog)

CALL timer( 'interp_lbc', 6)

IF (lbc_output_control%horizontal_interpolation_method /= create_frame) THEN
  ! Now we have an output field which has the LBC data on the LBC grid
  ! Set the LBC stashcode
  output_field%quantity_ident = stash_out
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION interp_control

END MODULE interp_control_mod

