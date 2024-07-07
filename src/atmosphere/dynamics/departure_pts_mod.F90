MODULE departure_pts_mod
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Description: Contains the ENDGame departure points
!
!
! Method:
!
! Documentation: ENDGame formulation version 1.01
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Atmosphere Dynamics
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

IMPLICIT NONE


REAL,SAVE, ALLOCATABLE ::                                             &
  depart_xi1_u  (:,:,:), depart_xi2_u  (:,:,:), depart_xi3_u  (:,:,:),&
  depart_xi1_v  (:,:,:), depart_xi2_v  (:,:,:), depart_xi3_v  (:,:,:),&
  depart_xi1_w  (:,:,:), depart_xi2_w  (:,:,:), depart_xi3_w  (:,:,:),&
  depart_xi1_rho(:,:,:), depart_xi2_rho(:,:,:), depart_xi3_rho(:,:,:)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DEPARTURE_PTS_MOD'

CONTAINS

SUBROUTINE init_depart_pts()

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE ereport_mod
USE atm_fields_bounds_mod

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_DEPART_PTS'

INTEGER :: ierr

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE (depart_xi1_u(udims%i_start:udims%i_end,                     &
                       udims%j_start:udims%j_end,                     &
                       udims%k_start:udims%k_end),                    &
          depart_xi2_u(udims%i_start:udims%i_end,                     &
                       udims%j_start:udims%j_end,                     &
                       udims%k_start:udims%k_end),                    &
          depart_xi3_u(udims%i_start:udims%i_end,                     &
                       udims%j_start:udims%j_end,                     &
                       udims%k_start:udims%k_end),                    &
          depart_xi1_v(vdims%i_start:vdims%i_end,                     &
                       vdims%j_start:vdims%j_end,                     &
                       vdims%k_start:vdims%k_end),                    &
          depart_xi2_v(vdims%i_start:vdims%i_end,                     &
                       vdims%j_start:vdims%j_end,                     &
                       vdims%k_start:vdims%k_end),                    &
          depart_xi3_v(vdims%i_start:vdims%i_end,                     &
                       vdims%j_start:vdims%j_end,                     &
                       vdims%k_start:vdims%k_end),                    &
          depart_xi1_w(wdims%i_start:wdims%i_end,                     &
                       wdims%j_start:wdims%j_end,                     &
                       wdims%k_start:wdims%k_end),                    &
          depart_xi2_w(wdims%i_start:wdims%i_end,                     &
                       wdims%j_start:wdims%j_end,                     &
                       wdims%k_start:wdims%k_end),                    &
          depart_xi3_w(wdims%i_start:wdims%i_end,                     &
                       wdims%j_start:wdims%j_end,                     &
                       wdims%k_start:wdims%k_end),                    &
        depart_xi1_rho(pdims%i_start:pdims%i_end,                     &
                       pdims%j_start:pdims%j_end,                     &
                       pdims%k_start:pdims%k_end),                    &
        depart_xi2_rho(pdims%i_start:pdims%i_end,                     &
                       pdims%j_start:pdims%j_end,                     &
                       pdims%k_start:pdims%k_end),                    &
        depart_xi3_rho(pdims%i_start:pdims%i_end,                     &
                       pdims%j_start:pdims%j_end,                     &
                       pdims%k_start:pdims%k_end),                    &
        stat=ierr)

IF (ierr/=0) CALL Ereport("init_depart_pts",ierr, "Unable to allocate.")

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE init_depart_pts

SUBROUTINE reset_dpt_pts()

USE horiz_grid_mod
USE atm_fields_bounds_mod, ONLY: udims,vdims,pdims
USE level_heights_mod,     ONLY: eta_rho_levels,eta_theta_levels
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER :: i,j,k
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RESET_DPT_PTS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,k)                               &
!$OMP          SHARED(pdims,udims,depart_xi1_u,depart_xi2_u,depart_xi3_u, &
!$OMP                 xi1_u,xi2_p,eta_rho_levels,vdims,depart_xi1_v,      &
!$OMP                 depart_xi2_v,depart_xi3_v,xi1_p,xi2_v,depart_xi1_w, &
!$OMP                 depart_xi2_w,depart_xi3_w,depart_xi1_rho,           &
!$OMP                 depart_xi2_rho,depart_xi3_rho,eta_theta_levels)
!$OMP DO SCHEDULE(STATIC)
DO k = 1, pdims%k_end
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      depart_xi1_u(i,j,k) = xi1_u(i)
      depart_xi2_u(i,j,k) = xi2_p(j)
      depart_xi3_u(i,j,k) = eta_rho_levels(k)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
DO k = 1, pdims%k_end
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      depart_xi1_v(i,j,k) = xi1_p(i)
      depart_xi2_v(i,j,k) = xi2_v(j)
      depart_xi3_v(i,j,k) = eta_rho_levels(k)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
DO k = 0, pdims%k_end
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      depart_xi1_w(i,j,k) = xi1_p(i)
      depart_xi2_w(i,j,k) = xi2_p(j)
      depart_xi3_w(i,j,k) = eta_theta_levels(k)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
DO k = 1, pdims%k_end
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      depart_xi1_rho(i,j,k) = xi1_p(i)
      depart_xi2_rho(i,j,k) = xi2_p(j)
      depart_xi3_rho(i,j,k) = eta_rho_levels(k)
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE


END MODULE departure_pts_mod
