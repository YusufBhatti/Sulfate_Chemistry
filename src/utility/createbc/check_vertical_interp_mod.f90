! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE check_vertical_interp_mod
!
! Description:
!   Determine if vertical interpolation is required
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

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'CHECK_VERTICAL_INTERP_MOD'

CONTAINS

SUBROUTINE check_vertical_interp(input_file, lbc_output_control, l_vert_interp)

USE datafile_mod,            ONLY: datafile_type
USE lbc_output_control_mod,  ONLY: lbc_output_control_type
USE interpor_mod,            ONLY: interp_order_linear, interp_order_linear_noex,     &
                                   interp_order_cubic, interp_order_quintic
USE ereport_mod,             ONLY: ereport
USE errormessagelength_mod,  ONLY: errormessagelength
USE umPrintMgr,              ONLY: umMessage, umPrint

IMPLICIT NONE

CLASS(datafile_type), INTENT(INOUT) :: input_file
CLASS(lbc_output_control_type), INTENT(IN) :: lbc_output_control
LOGICAL, INTENT(INOUT) :: l_vert_interp

CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER :: routinename='CHECK_VERTICAL_INTERP'
INTEGER :: icode
INTEGER :: level
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Check vertical interpolation is valid
IF (lbc_output_control%vertical_interpolation_method /= interp_order_linear        .AND. &
     lbc_output_control%vertical_interpolation_method /= interp_order_linear_noex  .AND. &
     lbc_output_control%vertical_interpolation_method /= interp_order_cubic        .AND. &
     lbc_output_control%vertical_interpolation_method /= interp_order_quintic) THEN
  icode = 10
  cmessage    = "[INFO] Vertical interpolation method not recognised."
  CALL ereport(routinename, icode, cmessage)
END IF

! Check if number of model levels has changed
IF (input_file%num_model_levels /= lbc_output_control%num_levels) THEN
  WRITE(umMessage,'(A)') '[INFO] Difference in no of model levels, turning on '//&
       'vertical interpolation'
  CALL umPrint(umMessage, src='check_vertical_interp')
  l_vert_interp = .TRUE.
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Check first rho level at which height is constant
IF (input_file%file_theta_rho_levels%first_rho_of_constant_height /=                       &
     lbc_output_control%p_grid(1)%vert_grid%first_rho_of_constant_height) THEN
  WRITE(umMessage,'(A)') '[INFO] Difference in first constant rho level, turning on '//    &
       'vertical interpolation'
  CALL umPrint(umMessage, src='check_vertical_interp')
  l_vert_interp = .TRUE.
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Check height of model
IF (ABS(input_file%file_theta_rho_levels%height_at_top_theta_level -                       &
     lbc_output_control%p_grid(1)%vert_grid%height_at_top_theta_level) > EPSILON(1.0)) THEN
  WRITE(umMessage,'(A)') '[INFO] Difference in height at top of model, turning on '//      &
       'vertical interpolation'
  CALL umPrint(umMessage, src='check_vertical_interp')
  l_vert_interp = .TRUE.
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Check if height values of the theta levels are the same
DO level=0, input_file%num_model_levels
  IF (ABS(input_file%file_theta_rho_levels%eta_theta(level) -                              &
       lbc_output_control%p_grid(1)%vert_grid%eta_theta(level))                            &
       > EPSILON(1.0)) THEN
    l_vert_interp = .TRUE.
    WRITE(umMessage,'(A)') '[INFO] Difference in eta values for theta levels, ' //         &
         'turning on vertical interpolation'
    CALL umPrint(umMessage, src='check_vertical_interp')
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END DO

! Check if height values of the rho levels are the same
DO level=1, input_file%num_model_levels
  IF (ABS(input_file%file_theta_rho_levels%eta_rho(level) -                                &
       lbc_output_control%p_grid(1)%vert_grid%eta_rho(level))                              &
       > EPSILON(1.0)) THEN
    l_vert_interp = .TRUE.
    WRITE(umMessage,'(A)') '[INFO] Difference in eta values for rho levels, turning on '// &
         'vertical interpolation'
    CALL umPrint(umMessage, src='check_vertical_interp')
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE check_vertical_interp

END MODULE check_vertical_interp_mod
