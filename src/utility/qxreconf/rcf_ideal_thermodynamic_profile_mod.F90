! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Idealised thermodynamic profile setup

MODULE rcf_ideal_thermodynamic_profile_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName = 'RCF_IDEAL_THERMODYNAMIC_PROFILE_MOD'

CONTAINS

! Description:
!   Sets up 1D profiles of thermodynamic variables for idealised problems.
!
! Method:
!   Given 1D surface temperature and pressure, and the temperature
!   profile option, 1D profiles of exner and potential temperature
!   are returned. These match the specified temperature profile
!   at theta-levels and are hydrostatically balanced.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: FORTRAN 2003
!   This code is written to UMDP3 programming standards.


SUBROUTINE rcf_ideal_thermodynamic_profile(model_levels, p_surface, z_theta,   &
                                           z_rho, g_theta, thetavd, exner)

USE prof_interp_mod,        ONLY: interp_linear, extrap_constant,              &
                                  prof_interp
USE profiles_mod,           ONLY: theta_dry, temperature, dtheta_dz,           &
                                  theta_init, Brunt_Vaisala
USE rcf_ideal_hydrostatic_from_temp_mod, ONLY: rcf_ideal_hydrostatic_from_temp
USE rcf_brunt_vaisala_prof_mod, ONLY: rcf_brunt_vaisala_prof
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

USE parkind1,               ONLY: jpim, jprb       !DrHook
USE yomhook,                ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE

! Input
INTEGER, INTENT(IN) :: model_levels
REAL, INTENT(IN)    :: p_surface
REAL, INTENT(IN)    :: z_theta(0:model_levels)
REAL, INTENT(IN)    :: z_rho(1:model_levels)
REAL, INTENT(IN)    :: g_theta(0:model_levels)

! Output
REAL, INTENT(OUT)   :: thetavd(0:model_levels)
REAL, INTENT(OUT)   :: exner(0:model_levels+1)

! Local
INTEGER             :: ierr, k
REAL                :: dz
REAL                :: temp(0:model_levels)

! DR HOOK
INTEGER(KIND=jpim), PARAMETER     :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER     :: zhook_out = 1
REAL(KIND=jprb)                   :: zhook_handle

CHARACTER(LEN=*), PARAMETER       :: & 
  RoutineName = 'RCF_IDEAL_THERMODYNAMIC_PROFILE'

! Error reporting
CHARACTER(LEN=errormessagelength) :: Cmessage

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (theta_init % field_type == dtheta_dz) THEN
! Construct theta profile: first element of input data is theta, subsequent
! elements are d(theta)/dz.
  DO k = 2, theta_init % n_heights
    dz                    = theta_init % height(k) - theta_init % height(k-1)
    theta_init % vprof(k,1) = theta_init % vprof(k-1,1) +                      &
                              theta_init % vprof(k,1) * dz
  END DO
END IF

! Interpolate input profile onto model grid
IF (theta_init % field_type /= Brunt_Vaisala) THEN
  CALL prof_interp(interp_linear, extrap_constant,                             &
                 theta_init % n_heights, theta_init % height,                  &
                 theta_init % vprof, model_levels+1, z_theta, temp)

END IF

SELECT CASE(theta_init % field_type)
CASE(theta_dry, dtheta_dz)
  thetavd(:) = temp(:)
CASE(temperature)
  CALL rcf_ideal_hydrostatic_from_temp(model_levels, temp, p_surface, z_theta, &
                                       z_rho, g_theta, thetavd, exner)
CASE(Brunt_Vaisala)
  CALL rcf_brunt_vaisala_prof(model_levels, theta_init % n_heights,            &
                              theta_init % height, theta_init % vprof,         &
                              z_theta, temp)
  thetavd(:) = temp(:)


CASE DEFAULT
  ierr = 1
  WRITE(Cmessage, '(''Unknown theta_init_field_type = '', I0)')                &
                  theta_init % field_type
  CALL ereport(RoutineName, ierr, Cmessage)
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE rcf_ideal_thermodynamic_profile
END MODULE rcf_ideal_thermodynamic_profile_mod
