! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE geostrophic_forcing_mod

IMPLICIT NONE
  ! Description: Apply geostrophic forcing
  !  
  !
  ! Method: Linearly interpolates, in time and height, time-varying
  !         profiles of geostrophic wind components. These are
  !         applied to the wind as forcing increments.
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
  ! This file belongs in section: IDEALISED
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !   This code is written to UMDP3 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='GEOSTROPHIC_FORCING_MOD'

CONTAINS

SUBROUTINE geostrophic_forcing(r_u, r_v)

USE atm_fields_bounds_mod,     ONLY: udims_s, vdims_s, udims, vdims
USE timestep_mod,              ONLY: timestep
USE profiles_mod,              ONLY: u_geostrophic, v_geostrophic
USE dyn_coriolis_mod,          ONLY: f3_at_u, f3_at_v
USE field_types,               ONLY: fld_type_u, fld_type_v
USE profile_uv_geo_mod,        ONLY: profile_uv_geo

USE umPrintMgr
USE parkind1,                  ONLY: jpim, jprb       !DrHook
USE yomhook,                   ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE

! Increments to horizontal wind components
REAL, INTENT (INOUT) :: r_u(udims_s%i_start:udims_s%i_end,                     &
                            udims_s%j_start:udims_s%j_end,                     &
                            udims_s%k_start:udims_s%k_end)

REAL, INTENT (INOUT) :: r_v(vdims_s%i_start:vdims_s%i_end,                     &
                            vdims_s%j_start:vdims_s%j_end,                     &
                            vdims_s%k_start:vdims_s%k_end)
 
! Local
INTEGER :: i, j, k
REAL    :: u_g(udims%k_start:udims%k_end)
REAL    :: v_g(vdims%k_start:vdims%k_end)

! DR HOOK
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GEOSTROPHIC_FORCING'

! 1.0 Start of subroutine code: perform the calculation.
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-------------- For time- and vertical-varying geostrophic forcing  ------------

CALL profile_uv_geo(u_geostrophic, udims, fld_type_u, u_g)

CALL profile_uv_geo(v_geostrophic, vdims, fld_type_v, v_g)


DO k = udims%k_start, udims%k_end
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      r_u(i,j,k) = r_u(i,j,k) - f3_at_u(i,j) * v_g(k) * timestep
    END DO
  END DO
END DO

DO k = vdims%k_start, vdims%k_end
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      r_v(i,j,k) = r_v(i,j,k) + f3_at_v(i,j) * u_g(k) * timestep
    END DO
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE geostrophic_forcing

END MODULE geostrophic_forcing_mod
