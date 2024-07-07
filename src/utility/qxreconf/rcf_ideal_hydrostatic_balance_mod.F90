! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! hydrostatic balance setup
MODULE rcf_ideal_hydrostatic_balance_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName = 'RCF_IDEAL_HYDROSTATIC_BALANCE_MOD'

CONTAINS

! Description:
!   Derives dry potential temperature (thetavd) and in turn
!   calculates hydrostatically balanced Exner and dry rho.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: FORTRAN 2003
!   This code is written to UMDP3 programming standards.
!
SUBROUTINE rcf_ideal_hydrostatic_balance(model_levels, p_surface, z_theta, &
                                         z_rho, g_theta, theta, mv,        &
                                         theta_vd, exner, dry_rho)

USE rcf_interp_weights_mod, ONLY: &
    intw_w2rho

USE planet_constants_mod, ONLY: cp, kappa, pref, R, recip_epsilon

USE parkind1,             ONLY: jpim, jprb       !DrHook
USE yomhook,              ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE

! Input
INTEGER, INTENT(IN) :: model_levels
REAL, INTENT(IN)    :: p_surface
REAL, INTENT(IN)    :: theta(0:model_levels)
REAL, INTENT(IN)    :: mv(0:model_levels)
REAL, INTENT(IN)    :: z_theta(0:model_levels)
REAL, INTENT(IN)    :: z_rho(1:model_levels)
REAL, INTENT(IN)    :: g_theta(0:model_levels)

! Output
REAL, INTENT(OUT)   :: theta_vd(0:model_levels)
REAL, INTENT(OUT)   :: exner(0:model_levels+1)
REAL, INTENT(OUT)   :: dry_rho(model_levels)

! Local
INTEGER             :: k
REAL                :: exner_lid, scale_ht_lid
! Virtual potential temperature
REAL                :: theta_v(0:model_levels)

! DR HOOK
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER   :: RoutineName = 'RCF_IDEAL_HYDROSTATIC_BALANCE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Virtual and virtual dry potential temperature
DO k = 0, model_levels
  theta_vd(k) = theta(k) * (1.0 + recip_epsilon * mv(k))
  theta_v(k)  = theta_vd(k) / (1.0 + mv(k))
END DO

! Exner at the surface
exner(0) = (p_surface/pref)**kappa - g_theta(0)*z_theta(0) / (cp*theta_v(0))

! Exner at first rho-level
exner(1 )= exner(0) - g_theta(0)*(z_rho(1)-z_theta(0)) / (cp*theta_v(0))

DO k = 1, model_levels-1
  exner(k+1) = exner(k) - g_theta(k)*(z_rho(k+1)-z_rho(k)) / (cp*theta_v(k))
END DO

! Exner at model_levels+1 is such that, when averaged with the value at
! model_levels, it provides Exner at the upper boundary that is
! obtained by linearly extrapolating log(exner)
scale_ht_lid = (LOG(exner(model_levels-1)) - LOG(exner(model_levels))) /       &
               (z_rho(model_levels)        - z_rho(model_levels-1))

exner_lid = exner(model_levels)*EXP(-scale_ht_lid*(z_theta(model_levels) -     &
                                                   z_rho(model_levels)))
exner(model_levels+1) = 2*exner_lid-exner(model_levels)

! Use equation of state to diagnose density
DO k = 1, model_levels
  dry_rho(k) = (pref*exner(k)**(1.0/kappa-1.0)) /                              &
           (R*(intw_w2rho(k,2)*theta_vd(k-1) + intw_w2rho(k,1)*theta_vd(k)))
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE rcf_ideal_hydrostatic_balance
END MODULE rcf_ideal_hydrostatic_balance_mod
