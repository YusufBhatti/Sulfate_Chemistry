! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine p_profile
MODULE p_profile_mod

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

USE missing_data_mod, ONLY: imdi

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'P_PROFILE_MOD'

CONTAINS

SUBROUTINE p_profile(p_surface, rho, p_prof_theta, p_prof_rho)
! Purpose:
!          Generate reference pressure profile.
!
! Method:
!          Calculate hydrostatically balanced pressure from horizontally
!          averaged density.
!
USE atm_fields_bounds_mod,  ONLY: tdims, pdims, pdims_s
USE nlsizes_namelist_mod,   ONLY: global_row_length, global_rows
USE global_2d_sums_mod,     ONLY: global_2d_sums
USE gravity_mod,            ONLY: g_ref_rho
USE level_heights_mod,      ONLY: z_ref_theta
USE idealise_run_mod,       ONLY: zero_orography

USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

USE parkind1,               ONLY: jpim, jprb       !DrHook
USE yomhook,                ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE

! Input
REAL, INTENT(IN)    :: p_surface
REAL, INTENT(IN)    :: rho(pdims_s%i_start:pdims_s%i_end,                      &
                           pdims_s%j_start:pdims_s%j_end,                      &
                           pdims_s%k_start:pdims_s%k_end)

! Output
REAL, INTENT(OUT)   :: p_prof_theta(tdims%k_start:tdims%k_end)
REAL, INTENT(OUT)   :: p_prof_rho(pdims%k_start:pdims%k_end)

! Local
INTEGER             :: k
REAL                :: rho_ave(pdims%k_start:pdims%k_end)

! Error reporting
CHARACTER(LEN=errormessagelength)  :: Cmessage
INTEGER                            :: ErrorStatus

! DR HOOK
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER   :: RoutineName = 'P_PROFILE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Check that orography is flat.
IF (.NOT. zero_orography) THEN
  ErrorStatus = 1
  Cmessage = 'Profiles specified against pressure levels - zero orography only'
  CALL ereport(RoutineName, ErrorStatus, Cmessage)
END IF

! Calculate horizontal average of density
CALL global_2d_sums(rho, pdims%i_len,    pdims%j_len,                          &
                         pdims_s%halo_i, pdims_s%halo_j,                       &
                         pdims%k_len,    rho_ave)

rho_ave(:) = rho_ave(:) / (global_row_length*global_rows)

p_prof_theta(0) = p_surface
DO k = 1, tdims%k_end
  p_prof_theta(k) = p_prof_theta(k-1) - rho_ave(k) * g_ref_rho(k) *            &
                                        (z_ref_theta(k) - z_ref_theta(k-1))
END DO

! Average in log pressure to rho-points - assumed half-way between theta-points
DO k = pdims%k_start, pdims%k_end
  p_prof_rho(k) = SQRT(p_prof_theta(k-1) * p_prof_theta(k))
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE p_profile

END MODULE p_profile_mod
