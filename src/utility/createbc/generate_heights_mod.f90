! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE generate_heights_mod
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

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'GENERATE_HEIGHTS_MOD'

CONTAINS

SUBROUTINE generate_heights(vert_grid, lbc_level_size, interp_orog, theta_heights, rho_heights)
! Description:
!   Calculate the theta and rho heights from a given set of eta values,
!   the height at the top of the model, the first constant rho level and
!   the interpolated orography on the LBC grid
!
!   Uses orography which has been interpolated to the LBC points on the
!   appropiate grid for the field (i.e. may need orography on u or v
!   points).
!
!   CreateBC only supports the "Smooth" method of height generation
!
USE planet_constants_mod,   ONLY: planet_radius
USE ereport_mod,            ONLY: ereport
USE vertical_grid_mod,      ONLY: vertical_grid_type
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Arguments
TYPE(vertical_grid_type), INTENT(IN) :: vert_grid
REAL, INTENT(IN)  :: interp_orog(:)
INTEGER, INTENT(IN) :: lbc_level_size
REAL, INTENT(OUT) :: theta_heights(lbc_level_size, 0:vert_grid%num_model_levels) 
REAL, INTENT(OUT) :: rho_heights(lbc_level_size, 1:vert_grid%num_model_levels+1)  

! Local variables
REAL :: ref_theta_heights(0:vert_grid%num_model_levels)
REAL :: ref_rho_heights(1:vert_grid%num_model_levels)

INTEGER :: i, k, icode

! Error reporting
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER :: routinename = 'GENERATE_HEIGHTS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Check that the LBC size matches orography array size
IF (lbc_level_size /= SIZE(interp_orog)) THEN
  icode = 10
  WRITE(cmessage, '(A,I10,A,I10)') "Mismatch between level size of field being interpolated =  ", &
      lbc_level_size, " and orography level size =  ",  SIZE(interp_orog)
  CALL ereport(routinename, icode, cmessage)
END IF

! Generate reference heights
DO k = 1, vert_grid%get_num_model_levels()
  ref_theta_heights(k) = vert_grid%eta_theta(k) * vert_grid%height_at_top_theta_level
  ref_rho_heights(k) = vert_grid%eta_rho(k) * vert_grid%height_at_top_theta_level
END DO

! Set bottom level equal to orography
DO i = 1, lbc_level_size
  theta_heights(i, 0) = interp_orog(i) + planet_radius
END DO

! All levels above and including the first constant rho level are just calculated 
! with reference to the top of the model
DO k = vert_grid%first_rho_of_constant_height, vert_grid%get_num_model_levels()
  DO i = 1, lbc_level_size
    theta_heights(i,k) = planet_radius + ref_theta_heights(k)
    rho_heights(i,k) = planet_radius + ref_rho_heights(k)
  END DO
END DO

! Levels below first constant rho level are set using a quadratic scheme
DO k = 1, vert_grid%first_rho_of_constant_height - 1
  DO i = 1, lbc_level_size
    theta_heights(i,k) = ref_theta_heights(k) + planet_radius + &
         interp_orog(i) * (1.0 - vert_grid%eta_theta(k) /       &
         vert_grid%eta_rho(vert_grid%first_rho_of_constant_height))**2

    rho_heights(i,k) = ref_rho_heights(k) + planet_radius +     &
         interp_orog(i) * (1.0 - vert_grid%eta_rho(k) /         &
         vert_grid%eta_rho(vert_grid%first_rho_of_constant_height))**2
  END DO
END DO

! Derive height for the extra rho level above model top    
DO i = 1, lbc_level_size
  rho_heights(i, vert_grid%get_num_model_levels()+1) = theta_heights(i, vert_grid%get_num_model_levels()) +   &
       (theta_heights(i, vert_grid%get_num_model_levels()) - rho_heights(i, vert_grid%get_num_model_levels()))
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE generate_heights

END MODULE generate_heights_mod
