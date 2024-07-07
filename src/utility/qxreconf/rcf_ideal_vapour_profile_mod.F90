! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Vapour mixing ratio profile setup
MODULE rcf_ideal_vapour_profile_mod

! Description:
!   Sets up 1D temp/pressure/mixing ratio for idealised problems.
!
! Method:
!   Given surface temperature and pressure, and the temperature
!   profile option, profiles of Exner, potential temperature and
!   vapour mixing ratio are returned. These match the specified
!   temperature profile at theta-levels and are hydrostatically balanced.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName = 'RCF_IDEAL_VAPOUR_PROFILE_MOD'

CONTAINS

SUBROUTINE rcf_ideal_vapour_profile(model_levels, z_theta, theta, exner, mv)

USE rcf_interp_weights_mod, ONLY: &
    intw_rho2w

USE prof_interp_mod,        ONLY: interp_linear, extrap_constant,              &
                                  prof_interp
USE profiles_mod,           ONLY: rel_humidity, mv_init

USE planet_constants_mod,   ONLY: pref, recip_kappa

USE parkind1,               ONLY: jpim, jprb       !DrHook
USE yomhook,                ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE

! Input
INTEGER, INTENT(IN) :: model_levels
REAL, INTENT(IN)    :: z_theta(0:model_levels)

! Output
REAL, INTENT(INOUT)   :: theta(0:model_levels)
REAL, INTENT(INOUT)   :: exner(0:model_levels+1)
REAL, INTENT(OUT)     :: mv(0:model_levels)

! Local
INTEGER             :: k
REAL                :: exner_theta(0:model_levels)
REAL                :: T(0:model_levels) ! Absolute temperature (K)
REAL                :: p_theta(0:model_levels), m_sat(0:model_levels)

! DR HOOK
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER   :: RoutineName = 'RCF_IDEAL_VAPOUR_PROFILE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL prof_interp(interp_linear, extrap_constant,                               &
                 mv_init % n_heights, mv_init % height, mv_init % vprof,       &
                 model_levels+1, z_theta, mv)

IF (mv_init % field_type == rel_humidity) THEN

  p_theta(0) = pref*exner(0)**recip_kappa
  exner_theta(0) = exner(0)
  DO k = 1, model_levels-1
    exner_theta(k) = intw_rho2w(k,1)*exner(k+1)+intw_rho2w(k,2)*exner(k)
    p_theta(k)     = pref*exner_theta(k)**recip_kappa
  END DO
  exner_theta(model_levels) = 0.5*(exner(model_levels)+exner(model_levels+1))
  p_theta(model_levels)     = pref*exner_theta(model_levels)**recip_kappa

  DO k = 0, model_levels
    T(k) = theta(k)*exner_theta(k)
  END DO

  ! DEPENDS ON: qsat_mix
  CALL qsat_mix(m_sat, T, p_theta, model_levels+1, .TRUE.)
  DO k = 0, model_levels
    mv(k) = mv(k)*m_sat(k)
  END DO

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE rcf_ideal_vapour_profile
END MODULE rcf_ideal_vapour_profile_mod
