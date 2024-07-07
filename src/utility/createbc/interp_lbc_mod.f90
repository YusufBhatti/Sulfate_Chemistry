! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE interp_lbc_mod
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

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'INTERP_LBC_MOD'

CONTAINS

TYPE(field_type) FUNCTION interp_lbc(input_field,                           &
                       target_grid,                                         &
                       horizontal_interpolation_method,                     &
                       vertical_interpolation_method,                       &
                       interp_weights, field_vert_interp,                   &
                       output_grid_stagger, interp_orog                     &
                       ) RESULT(output_field)

USE lbc_output_control_mod,     ONLY: lbc_output_control_type, create_frame, bilinear_interp
USE interp_weights_mod,         ONLY: interp_weights_type
USE f_shum_horizontal_field_interp_mod, ONLY :                                 &
                                     f_shum_horizontal_field_bi_lin_interp_calc
USE generate_heights_mod,       ONLY: generate_heights
USE ereport_mod,                ONLY: ereport
USE errormessagelength_mod,     ONLY: errormessagelength
USE stashmaster_constants_mod,  ONLY: first_level_code, last_level_code, &
                                      vert_level_type, theta_levels, rho_levels
USE stashmaster_utils_mod,      ONLY: query_stashmaster
USE vert_interp_mod,            ONLY: vert_interp
USE missing_data_mod,           ONLY: imdi
USE generate_frame_mod,         ONLY: generate_frame
USE field_mod,                  ONLY: field_type, ASSIGNMENT(=)
USE three_dimensional_grid_mod, ONLY: three_dimensional_grid_type
USE lbc_stashcode_mapping_mod,  ONLY: lbc_stashcode_mapping

USE packing_codes_mod, ONLY: PC_Cray32_Packing

! Method:
!  This function takes the input_field object from the argument list and 
!  returns the output_field object.  In order to facilitate the horizontal
!  and vertical interpolation an intermeditate field object is used, named
!  middle_field.
!
!  If no interpolation is being performed and the program is generating a
!  frame file rather than an LBC then the input_field is simply copied 
!  across to the output_field and the frames calculations performed on the
!  output field object.
!
!  If horizontal interpolation is being performed then input_field is 
!  bilinearly interpolated into the new middle_field object.  This field
!  is still on the same vertical level set as the original field.  If 
!  interpolation is necessary the middle_field is then interpolated
!  into the output_field object which is then on the output horizontal
!  and vertical grid.

IMPLICIT NONE

CLASS(field_type), INTENT(INOUT) :: input_field    ! The input_field

CLASS(three_dimensional_grid_type), INTENT(IN) :: target_grid
INTEGER, INTENT(IN) :: horizontal_interpolation_method
INTEGER, INTENT(IN) :: vertical_interpolation_method
CLASS(interp_weights_type), INTENT(IN) :: interp_weights
LOGICAL, INTENT(IN) :: field_vert_interp
INTEGER, INTENT(IN) :: output_grid_stagger
CLASS(field_type), OPTIONAL, INTENT(IN) :: interp_orog

! 2D arrays - 1D slice containing heights for each LBC level
REAL, ALLOCATABLE :: input_rho_heights(:,:)
REAL, ALLOCATABLE :: input_theta_heights(:,:)
REAL, ALLOCATABLE :: output_rho_heights(:,:)
REAL, ALLOCATABLE :: output_theta_heights(:,:)

TYPE(field_type) :: middle_field

INTEGER :: lbc_level_size
INTEGER :: ii, level, j, i ! Looper
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER :: routinename = 'INTERP_LBC'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (input_field%data_type /= 1) THEN
  ! Check if data type is real
  icode = 31
  WRITE(cmessage, '(A,I8)') "Currently only REAL data type supported, " // &
  "data type attempted was: ",input_field%data_type 
  CALL ereport(routinename, icode, cmessage)
END IF
! Set size of single level - this has already been determined when
! calculating the interpolation weights
lbc_level_size = SIZE(interp_weights%bilinear_index_b_l)

IF (horizontal_interpolation_method == bilinear_interp) THEN
  ! Copy across field metadata from input field as starting point
  CALL middle_field%copy_metadata_from_field(input_field)

  ! Cray-32 pack LBC fields
  middle_field%packing_method = PC_Cray32_Packing

  ! Allocate the middle field
  ALLOCATE( middle_field%lbc_rdata(SIZE(interp_weights%bilinear_index_b_l), &
       input_field%get_num_levels()))
  ! Perform the horizontal interpolation
  DO level=1, input_field%get_num_levels()
    CALL f_shum_horizontal_field_bi_lin_interp_calc(        &
         rows_in       = input_field%get_num_rows(),        &
         row_length_in = input_field%get_num_cols(),        &
         len_field_out = lbc_level_size,                    &
         index_b_l     = interp_weights%bilinear_index_b_l, &
         index_b_r     = interp_weights%bilinear_index_b_r, &
         index_t_l     = interp_weights%bilinear_index_t_l, &
         index_t_r     = interp_weights%bilinear_index_t_r, &
         data_in       = input_field%rdata(1,1,level),      &
         weight_b_l    = interp_weights%weight_b_l,         &
         weight_b_r    = interp_weights%weight_b_r,         &
         weight_t_l    = interp_weights%weight_t_l,         &
         weight_t_r    = interp_weights%weight_t_r,         &
         data_out      = middle_field%lbc_rdata(1,level))
  END DO
  ! Assign LBC grid to the middle field now that it has been interpolated
  middle_field%grid%horiz_grid = target_grid%horiz_grid
  ! Determine new stash code of the LBC field
  CALL lbc_stashcode_mapping(input_field%quantity_ident, middle_field%quantity_ident)

  IF (field_vert_interp) THEN
    ! Theta is dimensioned by model levels+1 due to zeroth level
    ! Rho is on model levels, but due to interpolation also need
    ! height at one level above top of model, so also dimension by model levels +1
    ALLOCATE(input_theta_heights(lbc_level_size,0:middle_field%get_num_model_levels()))
    ALLOCATE(input_rho_heights(lbc_level_size,1:middle_field%get_num_model_levels()+1))
    ALLOCATE(output_theta_heights(lbc_level_size,0:target_grid%get_num_model_levels()))
    ALLOCATE(output_rho_heights(lbc_level_size,1:target_grid%get_num_model_levels()+1))
    ! Generate heights for input vertical grid
    CALL generate_heights(input_field%grid%vert_grid, lbc_level_size, &
                          interp_orog%lbc_rdata(:,1),                 &
                          input_theta_heights, input_rho_heights)
    ! Generate heights for LBC output vertical grid
    CALL generate_heights(target_grid%vert_grid, lbc_level_size,      &
                          interp_orog%lbc_rdata(:,1),                 &
                          output_theta_heights, output_rho_heights)
    ! Copy field metadata and grid across to begin with
    CALL output_field%copy_metadata_from_field(middle_field)
    ! Assign LBC vertical grid to the horizontally interpolated field
    output_field%grid%vert_grid = target_grid%vert_grid
    ! Assign vertical level type i.e. theta or rho
    output_field%grid%vert_grid%level_type = query_stashmaster(       &
         output_field%quantity_ident, vert_level_type)
    ! Assign first and last level codes based on LBC stashcode
    output_field%grid%vert_grid%first_level_code = query_stashmaster( &
       output_field%quantity_ident, first_level_code)
    output_field%grid%vert_grid%last_level_code = query_stashmaster(  &
       output_field%quantity_ident, last_level_code)
    ! Assign the first and last level numbers of the vertical grid level codes
    CALL output_field%grid%vert_grid%calc_first_level(output_grid_stagger)
    CALL output_field%grid%vert_grid%calc_last_level(output_grid_stagger)
    ! Check that levels values returned are useful numbers
    IF (output_field%grid%vert_grid%first_level == imdi .OR.          &
         output_field%grid%vert_grid%last_level == imdi) THEN
      icode = 40
      WRITE(cmessage, '(A,I8,A)') "Level code for STASH item: ", output_field%quantity_ident, &
           " not supported"
      CALL ereport(routinename, icode, cmessage)
    END IF
    ! Now update the number of levels 
    CALL output_field%grid%vert_grid%recalc_num_levels()

    ! Allocate the space needed for the vertically interpolated field
    ALLOCATE(output_field%lbc_rdata(lbc_level_size, output_field%get_num_levels()))

    SELECT CASE (output_field%grid%vert_grid%level_type)
    CASE (theta_levels)
      j = 1
      DO level = output_field%grid%vert_grid%first_level, output_field%grid%vert_grid%last_level
        CALL vert_interp(                                                                     &
             data_in      = middle_field%lbc_rdata,                                           &
             data_points  = lbc_level_size,                                                   &
             data_levels  = middle_field%get_num_levels(),                                    &
             desired_r    = output_theta_heights(1:,level),                                   &
             r_at_data    = input_theta_heights(1:,middle_field%grid%vert_grid%first_level:), &
             interp_order = vertical_interpolation_method,                                    &
             data_out     = output_field%lbc_rdata(:,j))
        j = j + 1
      END DO
    CASE (rho_levels)
      j = 1
      DO level = output_field%grid%vert_grid%first_level, output_field%grid%vert_grid%last_level
        CALL vert_interp(                                                                     &
             data_in      = middle_field%lbc_rdata,                                           &
             data_points  = lbc_level_size,                                                   &
             data_levels  = middle_field%get_num_levels(),                                    &
             desired_r    = output_rho_heights(1:,level),                                     &
             r_at_data    = input_rho_heights(1:,middle_field%grid%vert_grid%first_level:),   &
             interp_order = vertical_interpolation_method,                                    &
             data_out     = output_field%lbc_rdata(:,j))
        j = j + 1
      END DO
    END SELECT
    ! Now the output field contains the target horizontal and vertical grids as well
    ! as the vertically interpolated data

    DEALLOCATE(output_rho_heights)
    DEALLOCATE(output_theta_heights)
    DEALLOCATE(input_rho_heights)
    DEALLOCATE(input_theta_heights)
  ELSE
    ! If not doing any vertical interpolation take the input vertical grid
    middle_field%grid%vert_grid = input_field%grid%vert_grid
    ! Copy everything from middle to output field, apart from the data
    ! The data may need adjustment to add and remove the zeroth level
    CALL output_field%copy_metadata_from_field(middle_field)

    ! Need to update first and last level numbers even if not doing vertical
    ! interpolation.  Some fields have different level numbers depending 
    ! on if they are a New Dynamics or ENDGame field.

    ! Update first and last level codes based on LBC stashcode
    output_field%grid%vert_grid%first_level_code = query_stashmaster(                         &
         output_field%quantity_ident, first_level_code)
    output_field%grid%vert_grid%last_level_code = query_stashmaster(                          &
         output_field%quantity_ident, last_level_code)
    ! Assign the first and last level numbers of the vertical grid level codes
    CALL output_field%grid%vert_grid%calc_first_level(output_grid_stagger)
    CALL output_field%grid%vert_grid%calc_last_level(output_grid_stagger)
    ! Check that levels values returned are useful numbers
    IF (output_field%grid%vert_grid%first_level == imdi .OR.                                  &
         output_field%grid%vert_grid%last_level == imdi) THEN
      icode = 40
      WRITE(cmessage, '(A,I8,A)') "Level code for STASH item: ", output_field%quantity_ident, &
           " not supported"
      CALL ereport(routinename, icode, cmessage)
    END IF
    ! Now update the number of levels 
    CALL output_field%grid%vert_grid%recalc_num_levels()
    ! Allocate space in the output field
    ALLOCATE(output_field%lbc_rdata(lbc_level_size,                                           &
                                    1:output_field%get_num_levels()))
    IF ( input_field%grid%vert_grid%first_level == 0 .AND.                                    &
         output_field%grid%vert_grid%first_level == 1 ) THEN
      ! ENDGame -> New Dynamics
      ! Need to remove the zeroth level by copying across from the second
      ! element of the data array onwards
      DO level = output_field%grid%vert_grid%first_level, output_field%grid%vert_grid%last_level
        output_field%lbc_rdata(:,level) = middle_field%lbc_rdata(:,level+1)
      END DO

    ELSE IF ( input_field%grid%vert_grid%first_level == 1 .AND.                               &
         output_field%grid%vert_grid%first_level == 0 ) THEN
      ! New Dynamics -> ENDGame
      ! Copy across data from first atmos level to last level, and leave
      ! room in the array for the zeroth level to be added in
      DO level = 1, output_field%grid%vert_grid%last_level
        output_field%lbc_rdata(:,level+1) = middle_field%lbc_rdata(:,level)
      END DO
      ! Need to add the zeroth level by copying across data from first atmos level
      output_field%lbc_rdata(:,1) = middle_field%lbc_rdata(:,1)
    ELSE
      ! New Dynamics -> New Dynamics
      ! ENDGame -> ENDGame
      DO level = 1, output_field%get_num_levels()
        output_field%lbc_rdata(:,level) = middle_field%lbc_rdata(:,level)
      END DO
    END IF
      
  END IF
ELSE IF (horizontal_interpolation_method == create_frame) THEN
  ! For now if you do not select horizontal interpolation will just copy across 
  ! full input field
  output_field = input_field
  ! Generate mask array using the input data array indexes contained in the interpolation
  ! weights type. If any of the input data points are not due to be used by the interpolation 
  ! routine we can set their value to missing data.
  CALL generate_frame(interp_weights, output_field)
ELSE
  icode = 50
  WRITE(cmessage, '(A,I0)') "Unsupported horizontal interpolation method: ", &
       horizontal_interpolation_method
  CALL ereport(routinename, icode, cmessage)
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION interp_lbc

END MODULE interp_lbc_mod
