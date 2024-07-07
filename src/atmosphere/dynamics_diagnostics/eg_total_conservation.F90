! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate sum of global mass, axial angular momentum and kinetic energy
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Atmosphere Dynamics Diagnostics

MODULE eg_total_conservation_mod
IMPLICIT NONE

REAL, SAVE :: total_rho_init, total_aam_init, total_ke_init

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_TOTAL_CONSERVATION_MOD'

CONTAINS

FUNCTION eg_total_gr(rho)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE eg_helmholtz_mod,      ONLY: ec_vol

USE atm_fields_bounds_mod, ONLY: pdims, pdims_s
USE global_2d_sums_mod, ONLY: global_2d_sums
USE level_heights_mod,  ONLY: xi3_at_rho=>r_rho_levels
USE gravity_mod,        ONLY: g_rho
USE planet_constants_mod,  ONLY: planet_radius

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_TOTAL_GR'

! Input
REAL    :: rho  (pdims_s%i_start:pdims_s%i_end,                  &
                 pdims_s%j_start:pdims_s%j_end,                  &
                 pdims_s%k_start:pdims_s%k_end)
! Output
REAL    :: eg_total_gr

! Local

INTEGER :: i, j, k
REAL    :: pe_mass(1)
REAL    :: k_sum(pdims%i_start:pdims%i_end,                      &
                 pdims%j_start:pdims%j_end)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

k_sum = 0.0
DO k = pdims%k_start,pdims%k_end
  DO j = pdims%j_start,pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      k_sum(i,j) = k_sum(i,j) + g_rho(i,j,k)*rho(i,j,k)       &
                  *(xi3_at_rho(i,j,k)-planet_radius)          &
                  *ec_vol(i,j,k)
    END DO
  END DO
END DO

CALL global_2d_sums(k_sum,pdims%i_end,pdims%j_end,0,0,1,pe_mass)

eg_total_gr = pe_mass(1)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION eg_total_gr


FUNCTION eg_total_mass(rho, l_exclude_rim)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE eg_helmholtz_mod,      ONLY: ec_vol

USE atm_fields_bounds_mod, ONLY: pdims, pdims_s
USE global_2d_sums_mod, ONLY: global_2d_sums
USE eg_total_mass_region_mod, ONLY: eg_total_mass_region

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_TOTAL_MASS'

! Input
REAL    :: rho  (pdims_s%i_start:pdims_s%i_end,                  &
                 pdims_s%j_start:pdims_s%j_end,                  &
                 pdims_s%k_start:pdims_s%k_end)
LOGICAL, INTENT(IN) :: l_exclude_rim

! Output
REAL    :: eg_total_mass

! Local

INTEGER :: i, j, k
REAL    :: pe_mass(1)
REAL    :: k_sum(pdims%i_start:pdims%i_end,                      &
                 pdims%j_start:pdims%j_end)
INTEGER :: IS, ie, js, je

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL eg_total_mass_region(IS, ie, js, je, l_exclude_rim)

k_sum = 0.0
DO k = pdims%k_start,pdims%k_end
  DO j = js, je
    DO i = IS, ie
      k_sum(i,j) = k_sum(i,j) + rho(i,j,k)*ec_vol(i,j,k)
    END DO
  END DO
END DO

CALL global_2d_sums(k_sum, pdims%i_end, pdims%j_end, 0, 0, 1, pe_mass)

eg_total_mass = pe_mass(1)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)


END FUNCTION eg_total_mass


FUNCTION eg_total_gr_r(rho)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE atm_fields_bounds_mod, ONLY: pdims, pdims_s
USE eg_helmholtz_mod,      ONLY: ec_vol

USE global_2d_sums_mod, ONLY: global_2d_sums
USE level_heights_mod,  ONLY: xi3_at_rho=>r_rho_levels
USE gravity_mod,        ONLY: g_rho
USE planet_constants_mod,  ONLY: planet_radius

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_TOTAL_GR_R'

! Input
REAL    :: rho  (pdims_s%i_start:pdims_s%i_end,                  &
                 pdims_s%j_start:pdims_s%j_end,                  &
                 pdims_s%k_start:pdims_s%k_end)
! Output
REAL    :: eg_total_gr_r

! Local

INTEGER :: i, j, k
REAL    :: pe_mass(1)
REAL    :: k_sum(pdims%i_start:pdims%i_end,                      &
                 pdims%j_start:pdims%j_end)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

k_sum = 0.0
DO k = pdims%k_start,pdims%k_end
  DO j = pdims%j_start,pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      k_sum(i,j) = k_sum(i,j) + g_rho(i,j,k)*rho(i,j,k)       &
                  *(xi3_at_rho(i,j,k)-planet_radius)          &
                  *xi3_at_rho(i,j,k)*ec_vol(i,j,k)
    END DO
  END DO
END DO

CALL global_2d_sums(k_sum,pdims%i_end,pdims%j_end,0,0,1,pe_mass)

eg_total_gr_r = pe_mass(1)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION eg_total_gr_r


FUNCTION eg_total_mass_r(rho)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE atm_fields_bounds_mod, ONLY: pdims, pdims_s
USE eg_helmholtz_mod,      ONLY: ec_vol

USE global_2d_sums_mod, ONLY: global_2d_sums
USE level_heights_mod,  ONLY: xi3_at_rho=>r_rho_levels

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_TOTAL_MASS_R'

! Input
REAL    :: rho  (pdims_s%i_start:pdims_s%i_end,                  &
                 pdims_s%j_start:pdims_s%j_end,                  &
                 pdims_s%k_start:pdims_s%k_end)
! Output
REAL    :: eg_total_mass_r

! Local

INTEGER :: i, j, k
REAL    :: pe_mass(1)
REAL    :: k_sum(pdims%i_start:pdims%i_end,                      &
                 pdims%j_start:pdims%j_end)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,k)                      &
!$OMP          SHARED(k_sum,pdims,rho,xi3_at_rho,ec_vol)
!$OMP DO SCHEDULE(STATIC)
DO j = pdims%j_start,pdims%j_end
  DO i = pdims%i_start,pdims%i_end
    k_sum(i,j)= 0.0
  END DO
END DO
!$OMP END DO
DO k = pdims%k_start,pdims%k_end
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start,pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      k_sum(i,j) = k_sum(i,j)                                 &
                  + rho(i,j,k)*xi3_at_rho(i,j,k)*ec_vol(i,j,k)
    END DO
  END DO
!$OMP END DO NOWAIT
END DO
!$OMP END PARALLEL

CALL global_2d_sums(k_sum, pdims%i_end, pdims%j_end, 0, 0, 1, pe_mass)

eg_total_mass_r = pe_mass(1)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION eg_total_mass_r

FUNCTION eg_total_ke(u, v, w, rho, l_exclude_rim)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE atm_fields_bounds_mod, ONLY: pdims, pdims_s, udims_s,        &
     vdims_s, wdims_s
USE horiz_grid_mod,        ONLY: intw_u2p, intw_v2p,intw_w2rho
USE eg_helmholtz_mod,      ONLY: ec_vol

USE global_2d_sums_mod, ONLY: global_2d_sums
USE eg_total_mass_region_mod, ONLY: eg_total_mass_region


IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_TOTAL_KE'

! Input
REAL    :: rho  (pdims_s%i_start:pdims_s%i_end,                  &
                 pdims_s%j_start:pdims_s%j_end,                  &
                 pdims_s%k_start:pdims_s%k_end)
REAL    :: u    (udims_s%i_start:udims_s%i_end,                  &
                 udims_s%j_start:udims_s%j_end,                  &
                 udims_s%k_start:udims_s%k_end)
REAL    :: v    (vdims_s%i_start:vdims_s%i_end,                  &
                 vdims_s%j_start:vdims_s%j_end,                  &
                 vdims_s%k_start:vdims_s%k_end)
REAL    :: w    (wdims_s%i_start:wdims_s%i_end,                  &
                 wdims_s%j_start:wdims_s%j_end,                  &
                 wdims_s%k_start:wdims_s%k_end)

LOGICAL, INTENT(IN) :: l_exclude_rim

! Output
REAL    :: eg_total_ke

! Local

REAL    :: u_ave, v_ave, w_ave

INTEGER :: i, j, k
REAL    :: pe_ke(1)
REAL    :: k_sum(pdims%i_start:pdims%i_end,                      &
                 pdims%j_start:pdims%j_end)
INTEGER :: IS, ie, js, je

! KE integral=Integral_Volume(0.5*M**(u^2+v^2+w^2))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL eg_total_mass_region(IS, ie, js, je, l_exclude_rim)

k_sum = 0.0
DO k = pdims%k_start,pdims%k_end
  DO j = js, je
    DO i = IS, ie
      u_ave=intw_u2p(i,1)*u(i-1,j,k)+intw_u2p(i,2)*u(i,j,k)
      v_ave=intw_v2p(j,1)*v(i,j-1,k)+intw_v2p(j,2)*v(i,j,k)
      w_ave=intw_w2rho(k,1)*w(i,j,k)+intw_w2rho(k,2)*w(i,j,k-1)
      k_sum(i,j) = k_sum(i,j) + 0.5*rho(i,j,k)*ec_vol(i,j,k)*         &
                                  (u_ave**2+v_ave**2+w_ave**2)
    END DO
  END DO
END DO

CALL global_2d_sums(k_sum, pdims%i_end, pdims%j_end, 0, 0, 1, pe_ke)

eg_total_ke = pe_ke(1)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION eg_total_ke

FUNCTION eg_total_aam(u, rho, L_shallow, l_exclude_rim)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE atm_fields_bounds_mod, ONLY: pdims, pdims_s, udims_s
USE horiz_grid_mod,        ONLY: intw_u2p, Csxi2_p
USE planet_constants_mod,  ONLY: planet_radius, omega
USE eg_helmholtz_mod,      ONLY: ec_vol

USE global_2d_sums_mod, ONLY: global_2d_sums
USE level_heights_mod,  ONLY: xi3_at_rho=>r_rho_levels
USE eg_total_mass_region_mod, ONLY: eg_total_mass_region

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_TOTAL_AAM'

! Input
REAL    :: rho  (pdims_s%i_start:pdims_s%i_end,                  &
                 pdims_s%j_start:pdims_s%j_end,                  &
                 pdims_s%k_start:pdims_s%k_end)
REAL    :: u    (udims_s%i_start:udims_s%i_end,                  &
                 udims_s%j_start:udims_s%j_end,                  &
                 udims_s%k_start:udims_s%k_end)

LOGICAL, INTENT(IN) :: L_shallow
LOGICAL, INTENT(IN) :: l_exclude_rim

! Output
REAL    :: eg_total_aam

! Local

REAL    :: r_prime
REAL    :: u_ave

INTEGER :: i, j, k
REAL    :: pe_aam(1)
REAL    :: k_sum(pdims%i_start:pdims%i_end,                      &
                 pdims%j_start:pdims%j_end)
INTEGER :: IS, ie, js, je

! aam=Integral_Volume(M*Omega*r*cos(phi)+u)*r*cos(phi)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL eg_total_mass_region(IS, ie, js, je, l_exclude_rim)

k_sum = 0.0
DO k = pdims%k_start,pdims%k_end
  DO j = js, je
     IF (L_shallow) r_prime=planet_radius*Csxi2_p(j)

    DO i = IS, ie
       IF (.NOT. L_shallow) r_prime=xi3_at_rho(i,j,k)*Csxi2_p(j)

       u_ave=intw_u2p(i,1)*u(i-1,j,k)+intw_u2p(i,2)*u(i,j,k)
       k_sum(i,j) = k_sum(i,j) + rho(i,j,k)*ec_vol(i,j,k)*          &
            r_prime*(omega*r_prime+u_ave)
    END DO
  END DO
END DO

CALL global_2d_sums(k_sum, pdims%i_end, pdims%j_end, 0, 0, 1, pe_aam)

eg_total_aam = pe_aam(1)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION eg_total_aam

END MODULE
