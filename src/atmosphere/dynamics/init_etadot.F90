! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE init_etadot_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='INIT_ETADOT_MOD'

CONTAINS

SUBROUTINE init_etadot()

USE atm_fields_bounds_mod
USE horiz_grid_mod, ONLY: intw_u2p,intw_v2p, intw_rho2w
USE atm_fields_mod, ONLY: u,v,w,etadot
USE metric_terms_mod,    ONLY:  h2_p_eta,h3_p_eta,h1_p_eta,&
                                deta_xi3_theta,dxi2_xi3,dxi1_xi3
USE field_types
USE halo_exchange, ONLY: swap_bounds
USE mpp_conf_mod,  ONLY: swap_field_is_scalar
USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE
!
! Description:
!
!
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Solver
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

INTEGER :: i,j,k
REAL :: u_at_w,v_at_w
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_ETADOT'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in, zhook_handle)

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)           &
!$OMP& PRIVATE(i,j,k,u_at_w,v_at_w) SHARED(pdims,intw_rho2w, &
!$OMP& intw_u2p,u,v,intw_v2p,etadot,w,h3_p_eta,dxi1_xi3,h1_p_eta, &
!$OMP& dxi2_xi3,h2_p_eta,deta_xi3_theta)
DO k = pdims%k_start, pdims%k_end-1
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

      u_at_w = intw_rho2w(k,1)*( intw_u2p(i,1)*u(i-1,j,k+1) +   &
                                 intw_u2p(i,2)*u(i,j,k+1) ) +   &
               intw_rho2w(k,2)*( intw_u2p(i,1)*u(i-1,j,k) +     &
                                 intw_u2p(i,2)*u(i,j,k) )

      v_at_w = intw_rho2w(k,1)*( intw_v2p(j,1)*v(i,j-1,k+1) +   &
                                 intw_v2p(j,2)*v(i,j,k+1) ) +   &
               intw_rho2w(k,2)*( intw_v2p(j,1)*v(i,j-1,k) +     &
                                 intw_v2p(j,2)*v(i,j,k) )

      etadot(i,j,k) = ( w(i,j,k)/h3_p_eta(i,j,k) -                  &
                         u_at_w*dxi1_xi3(i,j,k)/                    &
                                       h1_p_eta(i,j,k) -            &
                         v_at_w*dxi2_xi3(i,j,k)/                    &
                                       h2_p_eta(i,j,k) ) /          &
                                         deta_xi3_theta(i,j,k)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

etadot(:,:,0) = 0.0
etadot(:,:,pdims%k_end) = 0.0

CALL swap_bounds(etadot,                                               &
                 wdims_s%i_len - 2*wdims_s%halo_i,                     &
                 wdims_s%j_len - 2*wdims_s%halo_j,                     &
                 wdims_s%k_len,                                        &
                 wdims_s%halo_i, wdims_s%halo_j,                       &
                 fld_type_p,swap_field_is_scalar)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out, zhook_handle)

END SUBROUTINE
END MODULE
