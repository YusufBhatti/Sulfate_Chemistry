! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_sisl_consts_mod

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

IMPLICIT NONE

! Used in eg_sisl_init routines
REAL :: del_rho     
REAL :: rho_switch

! Used in eg_sisl_init_uvw routines
REAL, ALLOCATABLE  ::   dxi1_p(:), dxi2_p(:), deta_rho(:),         &
                        dxi1_u(:), dxi2_v(:), deta_w(:)
! Used in eg_helm_rhs_star routine
REAL, ALLOCATABLE  ::   rdxi1p(:), rdxi2p(:),                   &
                        rdxiu(:), rdxiv(:), rdxiw(:), rdrho(:)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_SISL_CONSTS_MOD'

CONTAINS

SUBROUTINE eg_sisl_consts( )

USE parkind1,          ONLY: jpim, jprb       !DrHook
USE yomhook,           ONLY: lhook, dr_hook   !DrHook

USE atm_fields_bounds_mod, ONLY: udims, vdims, pdims
USE nlsizes_namelist_mod,  ONLY: model_levels
USE level_heights_mod,     ONLY: eta_theta_levels, eta_rho_levels

USE planet_constants_mod, ONLY: kappa
USE eg_parameters_mod,    ONLY: l_rho_av_zz
USE horiz_grid_mod,       ONLY: xi1_u, xi1_p, xi2_v, xi2_p
USE eg_helmholtz_mod,     ONLY: hm_b, hm_pp, hm_rhoz, hm_theta

IMPLICIT NONE

! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_SISL_CONSTS'

! Local variables

INTEGER :: i, j, k

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ----------------------------------------------------
! 1.0 These are constant for whole run so set here
! ----------------------------------------------------

hm_pp = ( 1.0 - kappa ) / kappa

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j)                           &
!$OMP& SHARED( pdims, model_levels, hm_theta, hm_rhoz )

!$OMP DO SCHEDULE(STATIC)
DO j=pdims%j_start, pdims%j_end
  DO i=pdims%i_start, pdims%i_end
    hm_theta(i,j,0)            = 0.0
    hm_theta(i,j,model_levels) = 0.0
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC) 
DO j=pdims%j_start, pdims%j_end
  DO i=pdims%i_start, pdims%i_end
    hm_rhoz(i,j,0) = 0.0
    hm_rhoz(i,j,model_levels) = 0.0
  END DO
END DO
!$OMP END DO

!$OMP END PARALLEL

DO j=pdims%j_start, pdims%j_end
  DO i=pdims%i_start, pdims%i_end
    hm_b(i,j,0) = 0.0            !   hm_b(i,j,1)
    hm_b(i,j,model_levels) = 0.0 !   hm_b(i,j,model_levels-1)
  END DO
END DO

ALLOCATE ( dxi1_u(udims%i_start:udims%i_end) )
ALLOCATE ( dxi2_v(vdims%j_start:vdims%j_end) )
ALLOCATE ( dxi1_p(pdims%i_start:pdims%i_end) )
ALLOCATE ( dxi2_p(pdims%j_start:pdims%j_end) )
ALLOCATE ( deta_w(0:model_levels) )
ALLOCATE ( deta_rho(model_levels) )
ALLOCATE ( rdxiu(udims%i_start:udims%i_end) )
ALLOCATE ( rdxiv(vdims%j_start:vdims%j_end) )
ALLOCATE ( rdxi1p(pdims%i_start:pdims%i_end) )
ALLOCATE ( rdxi2p(pdims%j_start:pdims%j_end) )
ALLOCATE ( rdxiw(model_levels) )
ALLOCATE ( rdrho(model_levels) )

del_rho = 1.0

rho_switch = 1.0
IF ( l_rho_av_zz ) rho_switch = 0.0

DO i = pdims%i_start, pdims%i_end
  dxi1_p(i) = xi1_u(i)-xi1_u(i-1)
  rdxi1p(i) = 1.0 / dxi1_p(i)
END DO

DO i = udims%i_start, udims%i_end
  dxi1_u(i) = xi1_p(i+1)-xi1_p(i)
  rdxiu(i) = 1.0 / dxi1_u(i)
END DO

DO j = vdims%j_start, vdims%j_end
  dxi2_v(j) = xi2_p(j+1)-xi2_p(j)
  rdxiv(j) = 1.0 / dxi2_v(j)
END DO

DO j = pdims%j_start, pdims%j_end
  dxi2_p(j) = xi2_v(j)-xi2_v(j-1)
  rdxi2p(j) = 1.0 / dxi2_p(j)
END DO

k = 0
deta_w(k) = eta_rho_levels(k+1) - eta_theta_levels(k)

DO k = 1, model_levels - 1
  deta_w(k) = eta_rho_levels(k+1) - eta_rho_levels(k)
  deta_rho(k)   = eta_theta_levels(k) - eta_theta_levels(k-1)
  rdxiw(k) = 1.0 / deta_w(k)
  rdrho(k) = 1.0 / deta_rho(k)
END DO

k = model_levels
deta_rho(k) = eta_theta_levels(k) - eta_theta_levels(k-1)
rdrho(k) = 1.0 / deta_rho(k)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_sisl_consts
END MODULE eg_sisl_consts_mod
