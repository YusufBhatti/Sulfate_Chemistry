! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Hydrostatic profile for a given temp.
MODULE rcf_ideal_hydrostatic_from_temp_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName = 'RCF_IDEAL_HYDROSTATIC_FROM_TEMP_MOD'

CONTAINS

! Description:
!   Calculates hydrostatically balanced Exner and dry potential 
!   temperature (theta) from a given absolute temperature profile (T)
!   and surface pressure. The product of theta and 
!   vertically-averaged Exner matches T.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

SUBROUTINE rcf_ideal_hydrostatic_from_temp(model_levels, T, p_surface,         &
                                           z_theta, z_rho, g_theta,            &
                                           theta, exner)

USE rcf_interp_weights_mod, ONLY: &
    intw_rho2w

USE planet_constants_mod,   ONLY: cp, kappa, pref

USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

USE parkind1,               ONLY: jpim, jprb       !DrHook
USE yomhook,                ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE

! Input
INTEGER, INTENT(IN) :: model_levels
REAL, INTENT(IN)    :: p_surface
REAL, INTENT(IN)    :: T(0:model_levels)
REAL, INTENT(IN)    :: z_theta(0:model_levels)
REAL, INTENT(IN)    :: z_rho(1:model_levels)
REAL, INTENT(IN)    :: g_theta(0:model_levels)

! Output
REAL, INTENT(OUT)   :: theta(0:model_levels)
REAL, INTENT(OUT)   :: exner(0:model_levels+1)

! Local
INTEGER, PARAMETER  :: itn_max=1000
REAL, PARAMETER     :: error_tol=1.0e-11

INTEGER             :: itn, k
INTEGER             :: ErrorStatus
REAL                :: scale_ht
REAL                :: exner_surface
REAL                :: err_theta, err_exner, error
REAL                :: f_theta(0:model_levels)
REAL                :: f_exner(model_levels)
REAL                :: del_theta(0:model_levels)
REAL                :: del_exner(model_levels)
REAL                :: diag_theta
REAL                :: off_diag_theta
REAL                :: diag_exner
REAL                :: off_diag_exner

REAL                :: exner_lid, scale_ht_lid

! DR HOOK
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=errormessagelength) :: Cmessage
CHARACTER(LEN=*), PARAMETER       :: &
        RoutineName = 'RCF_IDEAL_HYDROSTATIC_FROM_TEMP'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! First guess
exner_surface = (p_surface/pref)**kappa
scale_ht      = g_theta(0) / (cp*T(0))
theta(0)      = T(0) / exner_surface

DO k = 1, model_levels
  exner(k) = exner_surface*EXP(-scale_ht*z_rho(k))
  theta(k) = theta(0)*EXP(scale_ht*z_theta(k))
END DO

del_theta(:) = 0.0
del_exner(:) = 0.0

DO itn = 1, itn_max

  f_theta(0) = theta(0)*exner_surface - T(0)
  f_exner(1) = g_theta(0) + cp*theta(0)*(exner(1)-exner_surface) /             &
                            (z_rho(1)-z_theta(0))

  DO k = 1, model_levels-1

    f_theta(k) = theta(k)*(intw_rho2w(k,2)*exner(k) +                          &
                           intw_rho2w(k,1)*exner(k+1)) - T(k)

    f_exner(k+1) = g_theta(k) + cp*theta(k)*(exner(k+1) - exner(k)) /          &
                                            (z_rho(k+1) - z_rho(k))

  END DO

  err_theta = MAX(MAXVAL(ABS(del_theta)), MAXVAL(ABS(f_theta)))
  err_exner = MAX(MAXVAL(ABS(del_exner)), MAXVAL(ABS(f_exner)))
  error     = MAX(err_theta, err_exner)
  IF (error < error_tol) EXIT

  diag_theta     = exner_surface
  off_diag_theta = 0.0
  del_theta(0)   = -(f_theta(0) + off_diag_theta) / diag_theta

  diag_exner     = cp*theta(0) / (z_rho(1) - z_theta(0))
  off_diag_exner = cp*del_theta(0)*(exner(1) - exner_surface) /                &
                                   (z_rho(1) - z_theta(0))
  del_exner(1)   = -(f_exner(1) + off_diag_exner) / diag_exner

  DO k = 1, model_levels-1

    diag_theta     = intw_rho2w(k,2)*exner(k) + intw_rho2w(k,1)*exner(k+1)
    off_diag_theta = theta(k)*(intw_rho2w(k,2)*del_exner(k) +                  &
                               intw_rho2w(k,1)*del_exner(k+1))
    del_theta(k)   = -(f_theta(k) + off_diag_theta) / diag_theta

    diag_exner     = cp*theta(k) / (z_rho(k+1) - z_rho(k))
    off_diag_exner = cp*del_theta(k)*(exner(k+1) - exner(k)) /                 &
                                     (z_rho(k+1) - z_rho(k))                   &
                     - cp*theta(k)*del_exner(k) / (z_rho(k+1) - z_rho(k))

    del_exner(k+1) = -(f_exner(k+1) + off_diag_exner) / diag_exner

  END DO

  theta(:)              = theta(:) + del_theta(:)
  exner(1:model_levels) = exner(1:model_levels) + del_exner(:)

END DO

IF (itn > itn_max) THEN
  ErrorStatus = 1
  WRITE(Cmessage, '(''Iteration failed: err_theta = '',' //                    &
                  'e12.4,'', err_exner = '', e12.4)') err_theta, err_exner
  CALL Ereport(RoutineName, ErrorStatus, Cmessage)
END IF

! Exner at model_levels+1 is such that, when averaged with the value at
! model_levels, it provides Exner at the upper boundary that is
! obtained by linearly extrapolating log(exner)
scale_ht_lid = (LOG(exner(model_levels-1)) - LOG(exner(model_levels))) /       &
               (z_rho(model_levels)        - z_rho(model_levels-1))

exner_lid    = exner(model_levels)*EXP(-scale_ht_lid*(z_theta(model_levels) -  &
                                                      z_rho(model_levels)))

exner(model_levels+1) = 2*exner_lid - exner(model_levels)
theta(model_levels)   = T(model_levels) / exner_lid

exner(0) = exner_surface

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE rcf_ideal_hydrostatic_from_temp
END MODULE rcf_ideal_hydrostatic_from_temp_mod
