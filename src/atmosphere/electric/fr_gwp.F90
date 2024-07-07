! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE fr_gwp_mod

! Purpose: Calculates flash rate using a simple graupel water path
!          relationship

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Electric

USE atm_fields_bounds_mod,  ONLY: tdims
USE electric_inputs_mod,    ONLY: g1, g2
USE electric_constants_mod, ONLY: i, j
USE cderived_mod,           ONLY: delta_lambda, delta_phi
USE level_heights_mod,      ONLY: r_theta_levels
USE trignometric_mod,       ONLY: fv_cos_theta_latitude, cos_theta_latitude

! Dr Hook modules
USE yomhook,                ONLY: lhook, dr_hook
USE parkind1,               ONLY: jprb, jpim

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='FR_GWP_MOD'

CONTAINS

SUBROUTINE fr_gwp(storm_field, gwp, flash)

IMPLICIT NONE

REAL, INTENT(IN)  :: gwp(   tdims%i_start : tdims%i_end,                       &
                            tdims%j_start : tdims%j_end )
! Graupel water path (kg m-2)

REAL, INTENT(INOUT) :: flash( tdims%i_start : tdims%i_end,                     &
                              tdims%j_start : tdims%j_end )
! Flash rate (s-1)

LOGICAL, INTENT(IN) :: storm_field(tdims%i_start : tdims%i_end,                &
                                   tdims%j_start : tdims%j_end )

! Local variables
REAL :: gbs  ! grid box size (square metres)

REAL, POINTER :: xx_cos_theta_latitude(:,:)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FR_GWP'

!==================================================================
! Start the subroutine
!==================================================================

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Grid box size variable dependent on latitude

! Calculate the trignometric constant dependent
! on dynamical core

xx_cos_theta_latitude => cos_theta_latitude

! Perform calculate of grid box size

DO j = tdims%j_start, tdims%j_end

  DO i = tdims%i_start, tdims%i_end

    ! Determine grid box size first
    gbs = r_theta_levels(i,j,1) * delta_lambda                          &
        * r_theta_levels(i,j,1) * delta_phi                             &
        * xx_cos_theta_latitude(i,j)

    IF ( storm_field(i,j) ) flash(i,j) = ( gbs * g1 * gwp(i,j) ) - g2

  END DO

END DO


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE fr_gwp

END MODULE fr_gwp_mod
