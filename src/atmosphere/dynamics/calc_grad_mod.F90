! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE calc_grad_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CALC_GRAD_MOD'

CONTAINS
SUBROUTINE calc_grad(field,field_surf,grad_u,grad_v,grad_w,l_compute_into_halo)

USE parkind1,              ONLY: jpim, jprb       !DrHook
USE yomhook,               ONLY: lhook, dr_hook   !DrHook

USE atm_fields_bounds_mod, ONLY: udims, vdims, wdims,                          &
                                 udims_s, vdims_s, wdims_s,                    &
                                 pdims, pdims_s

USE horiz_grid_mod,        ONLY: xi1_u, xi2_v, xi1_p, xi2_p,                   &
                                 intw_rho2w, intw_p2u, intw_p2v

USE level_heights_mod,     ONLY: eta_theta_levels, eta_rho_levels,             &
                                 xi3_at_theta => r_theta_levels,               &
                                 xi3_at_rho   => r_rho_levels,                 &
                                 xi3_at_u     => r_at_u,                       &
                                 xi3_at_v     => r_at_v

USE metric_terms_mod,      ONLY: deta_xi3_u, deta_xi3_v, deta_xi3,             &
                                 deta_xi3_theta,                               &
                                 h1_xi1_u, h2_xi2_v, h3_p_eta

IMPLICIT NONE
!
! Description:
!   Calculates gradient of a scalar field.
!
! Method: ENDGame formulation version 1.01,
!         section 7.2.
!  
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: DYNAMICS
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

! Subroutine arguments

! Input
REAL, INTENT(IN) :: field(pdims_s%i_start:pdims_s%i_end,                       &
                          pdims_s%j_start:pdims_s%j_end,                       &
                          pdims_s%k_start:pdims_s%k_end)

REAL, INTENT(IN) :: field_surf(pdims_s%i_start:pdims_s%i_end,                  &
                               pdims_s%j_start:pdims_s%j_end)

! Output
REAL, INTENT(OUT) :: grad_u(udims%i_start:udims%i_end,                         &
                            udims%j_start:udims%j_end,                         &
                            udims%k_start:udims%k_end)

REAL, INTENT(OUT) :: grad_v(vdims%i_start:vdims%i_end,                         &
                            vdims%j_start:vdims%j_end,                         &
                            vdims%k_start:vdims%k_end)

REAL, INTENT(OUT) :: grad_w(wdims%i_start:wdims%i_end,                         &
                            wdims%j_start:wdims%j_end,                         &
                            wdims%k_start:wdims%k_end)

LOGICAL, OPTIONAL :: l_compute_into_halo

! Local
INTEGER :: i, j, k                             ! loop counters
REAL    :: metric_coeff
REAL    :: dxi1_u(udims%i_start:udims%i_end),                                  &
           dxi1_p(pdims%i_start:pdims%i_end),                                  &
           dxi2_v(vdims%j_start:vdims%j_end),                                  &
           dxi2_p(pdims%j_start:pdims%j_end),                                  &
           deta_w(wdims%k_start:wdims%k_end),                                  &
           deta_rho(pdims%k_start:pdims%k_end)

REAL    :: work(wdims_s%i_start:wdims_s%i_end,                                 &
                wdims_s%j_start:wdims_s%j_end,                                 &
                wdims_s%k_start:wdims_s%k_end)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_GRAD'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Grid spacing
DO i=pdims%i_start, pdims%i_end
  dxi1_p(i) = xi1_u(i)-xi1_u(i-1)
END DO

DO i=udims%i_start, udims%i_end
  dxi1_u(i) = xi1_p(i+1)-xi1_p(i)
END DO

DO j=vdims%j_start, vdims%j_end
  dxi2_v(j) = xi2_p(j+1)-xi2_p(j)
END DO

DO j=pdims%j_start, pdims%j_end
  dxi2_p(j) = xi2_v(j)-xi2_v(j-1)
END DO

deta_w(0) = eta_rho_levels(1) - eta_theta_levels(0)
DO k=1, wdims%k_end-1
  deta_w(k) = eta_rho_levels(k+1) - eta_rho_levels(k)
END DO
deta_w(wdims%k_end) = eta_theta_levels(wdims%k_end) -                       &
     eta_rho_levels(wdims%k_end)

DO k=pdims%k_start,pdims%k_end
  deta_rho(k) = eta_theta_levels(k) - eta_theta_levels(k-1)
END DO

! ----------------------------------------------------------------------------
! U-component
! ----------------------------------------------------------------------------

DO j=udims%j_start, udims%j_end
  DO i=udims%i_start, udims%i_end
    work(i,j,0) = (intw_p2u(i,1)*field_surf(i,j) +                             &
                   intw_p2u(i,2)*field_surf(i+1,j)) *                          &
                  (xi3_at_theta(i+1,j,0)-xi3_at_theta(i,j,0)) / dxi1_u(i)
  END DO
END DO

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& PRIVATE(i,j,k)                                                          &
!$OMP& SHARED(udims,work,intw_rho2w,intw_p2u,field,xi3_at_theta,dxi1_u)
DO k=udims%k_start, udims%k_end-1
  DO j=udims%j_start, udims%j_end
    DO i=udims%i_start, udims%i_end
      work(i,j,k)=(intw_rho2w(k,1)*(intw_p2u(i,1)*field(i,j,k+1) +             &
                                    intw_p2u(i,2)*field(i+1,j,k+1)) +          &
                   intw_rho2w(k,2)*(intw_p2u(i,1)*field(i,j,k)+                &
                                    intw_p2u(i,2)*field(i+1,j,k))) *           &
                  (xi3_at_theta(i+1,j,k)-xi3_at_theta(i,j,k))/dxi1_u(i)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

! Assume that d(xi3)/d(xi1)=0.0 on upper boundary
DO j=udims%j_start, udims%j_end
  DO i=udims%i_start, udims%i_end
    work(i,j,udims%k_end)=0.0
  END DO
END DO

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& PRIVATE(i,j,k,metric_coeff)                                             &
!$OMP& SHARED(udims,grad_u,field,work,deta_xi3,dxi1_u,deta_rho,h1_xi1_u,       &
!$OMP& deta_xi3_u)
DO k=udims%k_start, udims%k_end
  DO j=udims%j_start, udims%j_end
    DO i=udims%i_start, udims%i_end
      metric_coeff=1.0/(h1_xi1_u(i,j,k)*deta_xi3_u(i,j,k))
      grad_u(i,j,k)=metric_coeff*((field(i+1,j,k)*deta_xi3(i+1,j,k) -          &
                                   field(i,j,k)*deta_xi3(i,j,k))/dxi1_u(i) -   &
                                  (work(i,j,k)-work(i,j,k-1))/deta_rho(k))
    END DO
  END DO
END DO
!$OMP END PARALLEL DO


! ----------------------------------------------------------------------------
! V-component
! ----------------------------------------------------------------------------
DO j=vdims%j_start, vdims%j_end
  DO i=vdims%i_start, vdims%i_end
    work(i,j,0) = (intw_p2v(j,1)*field_surf(i,j) +                             &
                   intw_p2v(j,2)*field_surf(i,j+1)) *                          &
                  (xi3_at_theta(i,j+1,0)-xi3_at_theta(i,j,0)) / dxi2_v(j)
  END DO
END DO

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& PRIVATE(i,j,k)                                                          &
!$OMP& SHARED(vdims,work,intw_rho2w,intw_p2v,field,xi3_at_theta,dxi2_v)
DO k=vdims%k_start, vdims%k_end-1
  DO j=vdims%j_start, vdims%j_end
    DO i=vdims%i_start, vdims%i_end
      work(i,j,k)=(intw_rho2w(k,1)*(intw_p2v(j,1)*field(i,j,k+1) +             &
                                    intw_p2v(j,2)*field(i,j+1,k+1)) +          &
                   intw_rho2w(k,2)*(intw_p2v(j,1)*field(i,j,k)+                &
                                    intw_p2v(j,2)*field(i,j+1,k))) *           &
                  (xi3_at_theta(i,j+1,k)-xi3_at_theta(i,j,k))/dxi2_v(j)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

! Assume that d(xi3)/d(xi2)=0.0 on upper boundary
DO j=vdims%j_start, vdims%j_end
  DO i=vdims%i_start, vdims%i_end
    work(i,j,vdims%k_end)=0.0
  END DO
END DO

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& PRIVATE(i,j,k,metric_coeff)                                             &
!$OMP& SHARED(vdims,grad_v,field,work,deta_xi3,dxi2_v,deta_rho,h2_xi2_v,       &
!$OMP&        deta_xi3_v)
DO k=vdims%k_start, vdims%k_end
  DO j=vdims%j_start, vdims%j_end
    DO i=vdims%i_start, vdims%i_end
      metric_coeff=1.0/(h2_xi2_v(i,j,k)*deta_xi3_v(i,j,k))
      grad_v(i,j,k)=metric_coeff*((field(i,j+1,k)*deta_xi3(i,j+1,k) -          &
                                   field(i,j,k)*deta_xi3(i,j,k))/dxi2_v(j) -   &
                                  (work(i,j,k)-work(i,j,k-1))/deta_rho(k))
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------------
! W-component
! ----------------------------------------------------------------------------

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& PRIVATE(i,j,k,metric_coeff)                                             &
!$OMP& SHARED(wdims,grad_w,field,deta_w,h3_p_eta,deta_xi3_theta)
DO k=wdims%k_start+1, wdims%k_end-1
  DO j=wdims%j_start, wdims%j_end
    DO i=wdims%i_start, wdims%i_end
      metric_coeff=1.0/(h3_p_eta(i,j,k)*deta_xi3_theta(i,j,k))
      grad_w(i,j,k)=metric_coeff*(field(i,j,k+1)-field(i,j,k))/deta_w(k)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

! Default option is zero vertical gradient at top and bottom boundaries
DO j=wdims%j_start, wdims%j_end
  DO i=wdims%i_start, wdims%i_end
    grad_w(i,j,wdims%k_start)=0.0
    grad_w(i,j,wdims%k_end)=0.0
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE calc_grad
END MODULE calc_grad_mod
