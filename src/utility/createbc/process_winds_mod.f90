! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE process_winds_mod
!
! Description:
!  Process wind fields if rotations are necessary.  If the source and the target
!  grids are on the same pole rotation then this routine will not be called and 
!  the winds will be processed in the same manner as the standard scalar fields.
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

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'PROCESS_WINDS_MOD'

CONTAINS

SUBROUTINE process_winds(field_number, lbc_output_control,                           &
                         input_file, output_file, orography_enlarged_grid_lbc_field, &
                         orography_lbc_fields, wind_rotation_coeff, interp_weights,  &
                         rotation_interp_weights, l_vert_interp)

USE datafile_mod,              ONLY: datafile_type
USE lbc_output_control_mod,    ONLY: lbc_output_control_type, create_frame
USE lbcfile_mod,               ONLY: lbcfile_type
USE field_mod,                 ONLY: field_type, ASSIGNMENT(=)
USE ereport_mod,               ONLY: ereport
USE errormessagelength_mod,    ONLY: errormessagelength
USE um_stashcode_mod,          ONLY: stashcode_u, stashcode_v, stashcode_u_adv, stashcode_v_adv
USE missing_data_mod,          ONLY: imdi
USE wind_rotation_coeff_mod,   ONLY: wind_rotation_coeff_type
USE unrotate_input_winds_mod,  ONLY: unrotate_input_winds
USE umPrintMgr,                ONLY: umPrint, umMessage
USE post_interp_transform_mod, ONLY: post_interp_transform
USE process_orography_mod,     ONLY: u_points_orog, v_points_orog, p_points_orog
USE interp_weights_mod,        ONLY: interp_weights_type, p_to_enlarged_p, u_to_enlarged_p, &
                                     v_to_enlarged_p, weights_index
USE three_dimensional_grid_mod,ONLY: three_dimensional_grid_type
USE stashmaster_constants_mod, ONLY: u_points, v_points
USE interp_lbc_mod,            ONLY: interp_lbc
USE um_parparams,              ONLY: halo_type_extended
USE rotate_output_winds_mod,   ONLY: rotate_output_winds
USE find_fields_to_interp_mod, ONLY: fields_to_interp

USE interp_input_winds_mod,   ONLY: interp_input_wind_v_to_p, interp_input_wind_u_to_p, &
                                    interp_input_wind_p_to_v, interp_input_wind_p_to_u
USE interp_output_winds_mod,  ONLY: interp_output_wind_p_to_v, interp_output_wind_p_to_u
USE update_frame_field_grid_mod,  ONLY: update_frame_field_grid
IMPLICIT NONE

INTEGER, INTENT(IN)    :: field_number
TYPE(lbc_output_control_type),   TARGET, INTENT(IN) :: lbc_output_control
TYPE(field_type),                 TARGET, INTENT(IN) :: orography_enlarged_grid_lbc_field
TYPE(field_type),                 TARGET, INTENT(IN) :: orography_lbc_fields(3,3)
TYPE(wind_rotation_coeff_type), INTENT(IN) :: wind_rotation_coeff
CLASS(interp_weights_type), TARGET, INTENT(IN) :: interp_weights(:,:)
CLASS(interp_weights_type), TARGET, INTENT(IN) :: rotation_interp_weights(:)
LOGICAL, INTENT(IN) :: l_vert_interp ! Is vertical interpolation required?
CLASS(datafile_type), TARGET, INTENT(INOUT) :: input_file
CLASS(datafile_type), INTENT(INOUT) :: output_file

! Local variables
INTEGER :: icode = 0
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER :: routinename = 'PROCESS_WINDS'
! Output fields
TYPE(field_type), SAVE :: output_u_field
TYPE(field_type), SAVE :: output_v_field
! Input fields, locations in file
INTEGER :: field_location
INTEGER :: next_field_location
! Pointers which will be associated with different types depending on 
! the particular rotations of the input and output grids
TYPE(field_type), POINTER :: input_u_field
TYPE(field_type), POINTER :: input_v_field
TYPE(field_type), POINTER :: u_lbc_orography
TYPE(field_type), POINTER :: v_lbc_orography
TYPE(interp_weights_type), POINTER :: u_interp_weights
TYPE(interp_weights_type), POINTER :: v_interp_weights
TYPE(three_dimensional_grid_type), POINTER :: u_target_grid
TYPE(three_dimensional_grid_type), POINTER :: v_target_grid
! Output field numbers
INTEGER :: u_field_num, v_field_num
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Get the field location of the current field
field_location = fields_to_interp(field_number)
! Get the field location of the next field due to be processed
next_field_location = fields_to_interp(field_number + 1)

! Check that the field number which has been passed into routine was a U (adv U)
! field
IF (input_file%fields(field_location)%quantity_ident /= stashcode_u .AND.                   &
    input_file%fields(field_location)%quantity_ident /= stashcode_u_adv) THEN
  WRITE(cmessage, '(A,I0)') "Process wind routine requires U or Advected U field " //       &
        "to be first wind field, please check stash_codes array in namelist. " //           &
        "STASHcode of field passed to routine was: ",                                       &
       input_file%fields(field_location)%quantity_ident
  icode = 10
  CALL ereport(routinename, icode, cmessage)
END IF

IF (next_field_location <=  1) THEN
  ! The fields_to_interp array is defaulted to -99, if the field number passed
  ! into this routine is the last item in the stash_code namelist array then the next
  ! element will be -99. Catch this and ereport.
  WRITE(cmessage, '(A,I0)') "Cannot determine next field. Process wind routine " //         &
       "requires V or adv V field STASHcode to be the next STASHcode in the " //            &
       "stash_codes namelist array. STASHcode of current field was: ",                      &
       input_file%fields(field_location)%quantity_ident
  icode = 11
  CALL ereport(routinename, icode, cmessage)
END IF

! Check that the next field in the fields_to_interp object is the corresponding
! V (adv V) field
IF (input_file%fields(field_location)%quantity_ident == stashcode_u) THEN
  IF (input_file%fields(next_field_location)%quantity_ident /= stashcode_v) THEN
    WRITE(cmessage, '(A,I0)') "If rotating winds CreateBC requires V field STASHcode " //   &
         "to immediatly follow the U field stashcode in the list of STASHcodes read " //    &
         "in via namelist.  However process_winds found the next STASHcode to be: ",        &
         input_file%fields(next_field_location)%quantity_ident
    icode = 20
    CALL ereport(routinename, icode, cmessage)
  END IF
END IF
IF (input_file%fields(field_location)%quantity_ident == stashcode_u_adv) THEN
  IF (input_file%fields(next_field_location)%quantity_ident /= stashcode_v_adv) THEN
    WRITE(cmessage, '(A,I0)') "If rotating winds CreateBC requires advected V field " //    &
         "STASHcode to immediatly follow the advected U field stashcode in the list of " // &
         "STASHcodes read in via namelist.  However, found the next STASHcode " //          &
         "to be: ", input_file%fields(next_field_location)%quantity_ident
    icode = 30
    CALL ereport(routinename, icode, cmessage)
  END IF
END IF

! Read in U wind and V winds
CALL input_file%read_field(field_location)
input_u_field => input_file%fields(field_location)
CALL input_file%read_field(next_field_location)
input_v_field => input_file%fields(next_field_location)

IF (lbc_output_control% horizontal_interpolation_method /= create_frame) THEN
  ! Do not want to modify input data if simply creating frame
  IF (input_file%l_source_rotated) THEN
    ! Interpolate U and V winds to the P grid 
    CALL interp_input_wind_u_to_p(input_u_field, input_file)
    CALL interp_input_wind_v_to_p(input_v_field, input_file)
    ! Apply the rotation coefficients to rotate the wind vectors
    CALL unrotate_input_winds(input_u_field, input_v_field, wind_rotation_coeff)
    ! Interpolate U and V winds from the P grid back to the original U and V grids
    CALL interp_input_wind_p_to_u(input_u_field, input_file)
    CALL interp_input_wind_p_to_v(input_v_field, input_file) 
  END IF ! input_file%l_source_rotated
END IF

! Set pointers to interpolation weight, orography field and target grids
IF (lbc_output_control%l_target_rotated) THEN
  ! Input grid is on U and V points and need to interpolate to the LBC enlarged P grid
  u_interp_weights => rotation_interp_weights(u_to_enlarged_p)
  v_interp_weights => rotation_interp_weights(v_to_enlarged_p)
  ! Use the enlarged LBC orography on p points
  u_lbc_orography => orography_enlarged_grid_lbc_field
  v_lbc_orography => orography_enlarged_grid_lbc_field
  ! Use enlarged P target grid
  u_target_grid => lbc_output_control%p_grid_enlarged
  v_target_grid => lbc_output_control%p_grid_enlarged
ELSE 
  ! Input grid is on U and V points and need to interpolate to the LBC U and V grids
  u_interp_weights => interp_weights(weights_index(u_points), halo_type_extended)
  v_interp_weights => interp_weights(weights_index(v_points), halo_type_extended)
  ! Use the standard LBC orography on u and v points
  u_lbc_orography => orography_lbc_fields(u_points_orog, halo_type_extended)
  v_lbc_orography => orography_lbc_fields(v_points_orog, halo_type_extended)
  ! Interpolate to the standard U and V target grids
  u_target_grid => lbc_output_control%u_grid(halo_type_extended)
  v_target_grid => lbc_output_control%v_grid(halo_type_extended)
END IF

! Check arrays are allocated
IF (.NOT. ALLOCATED(u_interp_weights%bilinear_index_b_l)) THEN
  icode = 40
  cmessage = "u_interp_weights is not allocated, unable to perform interpolation"
ELSE IF (.NOT. ALLOCATED(v_interp_weights%bilinear_index_b_l)) THEN
  icode = 41
  cmessage = "v_interp_weights is not allocated, unable to perform interpolation"
END IF
IF (lbc_output_control% horizontal_interpolation_method /= create_frame) THEN
  ! Only need to check orography fields if not creating frame
  IF (.NOT. ALLOCATED(u_lbc_orography%lbc_rdata)) THEN
    icode = 42
    cmessage = "u lbc orography field is not allocated, unable to perform interpolation"
  ELSE IF (.NOT. ALLOCATED(v_lbc_orography%lbc_rdata)) THEN
    icode = 43
    cmessage = "v lbc orography field is not allocated, unable to perform interpolation"
  END IF
END IF
IF (icode > 0) CALL ereport(routinename, icode, cmessage)

! Perform the interpolation
output_u_field = interp_lbc(input_u_field, u_target_grid,            &
                 lbc_output_control%horizontal_interpolation_method, &
                 lbc_output_control%vertical_interpolation_method,   &
                 u_interp_weights, l_vert_interp,                    &
                 lbc_output_control%output_grid_stagger,             &
                 interp_orog=u_lbc_orography)
output_v_field = interp_lbc(input_v_field, v_target_grid,            &
                 lbc_output_control%horizontal_interpolation_method, &
                 lbc_output_control%vertical_interpolation_method,   &
                 v_interp_weights, l_vert_interp,                    &
                 lbc_output_control%output_grid_stagger,             &
                 interp_orog=v_lbc_orography)
! No longer need input data from memory
CALL input_u_field%unload_data()
CALL input_v_field%unload_data()

WRITE(umMessage,'(A,I5,A,A)') '[INFO] Generating field for STASH ',  & 
     output_u_field%quantity_ident,' at time ',output_u_field%return_validity_time_string()
CALL umPrint(umMessage, src='createbc')
WRITE(umMessage,'(A,I5,A,A)') '[INFO] Generating field for STASH ',  & 
     output_v_field%quantity_ident,' at time ',output_v_field%return_validity_time_string()
CALL umPrint(umMessage, src='createbc')

IF (lbc_output_control% horizontal_interpolation_method /= create_frame) THEN
  ! Do not want to modify output field if producing frame
  ! Rotate winds to pole rotation of the target grid
  IF (lbc_output_control%l_target_rotated) THEN
    ! Rotate the wind U and V vectors
    CALL rotate_output_winds(wind_rotation_coeff, output_u_field, output_v_field)
    ! Move data from enlarged P grid to the U and V grids
    CALL interp_output_wind_p_to_u(output_u_field, lbc_output_control, interp_weights)
    CALL interp_output_wind_p_to_v(output_v_field, lbc_output_control, interp_weights) 
  END IF ! lbc_output_control%l_target_rotated
  ! Post interp transformations and corrections
  CALL post_interp_transform(output_u_field, lbc_output_control%q_min)
  CALL post_interp_transform(output_v_field, lbc_output_control%q_min)
END IF

! Add output fields to the output file object and write field to disk
u_field_num = output_file%add_field(output_u_field)
IF (lbc_output_control% horizontal_interpolation_method == create_frame) THEN
  CALL update_frame_field_grid(output_file, u_field_num)
END IF
CALL output_file%write_field(u_field_num)
CALL output_file%fields(u_field_num)%unload_data()
v_field_num = output_file%add_field(output_v_field)
IF (lbc_output_control% horizontal_interpolation_method == create_frame) THEN
  CALL update_frame_field_grid(output_file, v_field_num)
END IF
CALL output_file%write_field(v_field_num)
CALL output_file%fields(v_field_num)%unload_data()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE process_winds

END MODULE process_winds_mod