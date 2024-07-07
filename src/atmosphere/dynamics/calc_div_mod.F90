! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Description: Calculate divergence (see subroutine description)
!
! Method: See subroutine description
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: DYNAMICS ADVECTION
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

MODULE calc_div_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CALC_DIV_MOD'

CONTAINS
SUBROUTINE calc_div(u,v,etadot,divdims,div,l_compute_into_halo)

USE parkind1,              ONLY: jpim, jprb       !DrHook
USE yomhook,               ONLY: lhook, dr_hook   !DrHook

USE atm_fields_bounds_mod, ONLY: udims, vdims, wdims,                          &
                                 udims_s, vdims_s, wdims_s,                    &
                                 pdims, pdims_s, array_dims

USE horiz_grid_mod,        ONLY: xi1_u, xi2_v

USE level_heights_mod,     ONLY: eta_theta_levels, eta_rho_levels,             &
                                 xi3_at_theta => r_theta_levels,               &
                                 xi3_at_rho   => r_rho_levels,                 &
                                 xi3_at_u     => r_at_u,                       &
                                 xi3_at_v     => r_at_v

USE metric_terms_mod,      ONLY: deta_xi3_u, deta_xi3_v, deta_xi3,             &
                                 deta_xi3_theta,                               &
                                 h1_p, h1_xi2_v, h1_p_eta,                     &
                                 h2_p, h2_xi1_u, h2_p_eta,                     &
                                 h3_p, h3_xi1_u, h3_xi2_v, h3_p_eta

USE halo_exchange,         ONLY: swap_bounds
USE mpp_conf_mod,          ONLY: swap_field_is_scalar
USE field_types,           ONLY: fld_type_p

IMPLICIT NONE
!
! Description:
!   Calculates divergence of a vector field, whose components are
!   distributed on a C-grid precisely as the wind field.
!
! Method: ENDGame formulation version 1.01,
!         section 7.2.
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

! Input
REAL, INTENT(IN) :: u(udims_s%i_start:udims_s%i_end,                           &
                      udims_s%j_start:udims_s%j_end,                           &
                      udims_s%k_start:udims_s%k_end)

REAL, INTENT(IN) :: v(vdims_s%i_start:vdims_s%i_end,                           &
                      vdims_s%j_start:vdims_s%j_end,                           &
                      vdims_s%k_start:vdims_s%k_end)

REAL, INTENT(IN) :: etadot(wdims_s%i_start:wdims_s%i_end,                      &
                           wdims_s%j_start:wdims_s%j_end,                      &
                           wdims_s%k_start:wdims_s%k_end)

TYPE(array_dims), INTENT(IN) :: divdims

LOGICAL, OPTIONAL :: l_compute_into_halo

! Output
REAL, INTENT(OUT) :: div(divdims%i_start:divdims%i_end,                        &
                         divdims%j_start:divdims%j_end,                        &
                         divdims%k_start:divdims%k_end)

! Local
INTEGER :: i, j, k                             ! loop counters
INTEGER :: i_start, i_end
INTEGER :: j_start, j_end
INTEGER :: k_start, k_end

REAL    :: rdxi1, rdxi2, rdxi3                 ! reciprocals of grid spacing
REAL    :: d_xi1_term, d_xi2_term, deta_term   ! three terms of divergence
REAL    :: Jacobian

LOGICAL :: l_swap_bounds

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_DIV'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set loop limits
i_start=divdims%i_start+divdims%halo_i
i_end=divdims%i_end-divdims%halo_i

j_start=divdims%j_start+divdims%halo_j
j_end=divdims%j_end-divdims%halo_j

k_start=divdims%k_start
k_end=divdims%k_end

! Do swap bounds if divergence array has halos
l_swap_bounds=(divdims%halo_i /= 0 .OR. divdims%halo_j /= 0)

! Modify loop limits if computing into the halos
IF (PRESENT(l_compute_into_halo)) THEN
  IF (l_compute_into_halo .AND. divdims%halo_i * divdims%halo_j > 0) THEN
    i_start=divdims%i_start+divdims%halo_i
    i_end=divdims%i_end
    j_start=divdims%j_start+divdims%halo_j
    j_end=divdims%j_end
    l_swap_bounds=.FALSE.
  END IF
END IF

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& PRIVATE(i,j,k,rdxi1,rdxi2,rdxi3,d_xi1_term,d_xi2_term,deta_term,        &
!$OMP&         Jacobian)                                                       &
!$OMP& SHARED(i_start,i_end,j_start,j_end,k_start,k_end,                       &
!$OMP&        eta_theta_levels,xi1_u,xi2_v,h2_xi1_u,h3_xi1_u,                  &
!$OMP&        deta_xi3_u,h1_xi2_v,h3_xi2_v,deta_xi3_v,h1_p_eta,                &
!$OMP&        h2_p_eta,h3_p_eta,deta_xi3_theta,h1_p,h2_p,h3_p,deta_xi3,        &
!$OMP&        u,v,etadot,div)
DO k=k_start+1, k_end-1
  rdxi3 = 1.0/( eta_theta_levels(k) - eta_theta_levels(k-1) )
  DO j=j_start, j_end
    rdxi2 = 1.0/( xi2_v(j) - xi2_v(j-1) )
    DO i=i_start, i_end
      rdxi1 = 1.0/( xi1_u(i) - xi1_u(i-1) )

      d_xi1_term = ( h2_xi1_u(i,j,k)*h3_xi1_u(i,j,k)*deta_xi3_u(i,j,k)*        &
                     u(i,j,k) -                                                &
                     h2_xi1_u(i-1,j,k)*h3_xi1_u(i-1,j,k)*deta_xi3_u(i-1,j,k)*  &
                     u(i-1,j,k) ) * rdxi1

      d_xi2_term = ( h1_xi2_v(i,j,k)*h3_xi2_v(i,j,k)*deta_xi3_v(i,j,k)*        &
                     v(i,j,k) -                                                &
                     h1_xi2_v(i,j-1,k)*h3_xi2_v(i,j-1,k)*deta_xi3_v(i,j-1,k)*  &
                     v(i,j-1,k) ) * rdxi2

      deta_term = ( h1_p_eta(i,j,k)*h2_p_eta(i,j,k)*                           &
                    h3_p_eta(i,j,k)*deta_xi3_theta(i,j,k)*                     &
                    etadot(i,j,k) -                                            &
                    h1_p_eta(i,j,k-1)*h2_p_eta(i,j,k-1)*                       &
                    h3_p_eta(i,j,k-1)*deta_xi3_theta(i,j,k-1)*                 &
                    etadot(i,j,k-1)) * rdxi3

      Jacobian  = h1_p(i,j,k)*h2_p(i,j,k)*h3_p(i,j,k)*deta_xi3(i,j,k)

      div(i,j,k) = ( d_xi1_term + d_xi2_term + deta_term ) / Jacobian

    END DO
  END DO
END DO
!$OMP END PARALLEL DO

DO j=j_start, j_end
  rdxi2 = 1.0/( xi2_v(j) - xi2_v(j-1) )
  DO i=i_start, i_end
    rdxi1 = 1.0/( xi1_u(i) - xi1_u(i-1) )

! Lower boundary
    k=k_start
    rdxi3 = 1.0/( eta_theta_levels(k) - eta_theta_levels(k-1) )
    d_xi1_term = ( h2_xi1_u(i,j,k)*h3_xi1_u(i,j,k)*deta_xi3_u(i,j,k)*          &
                   u(i,j,k) -                                                  &
                   h2_xi1_u(i-1,j,k)*h3_xi1_u(i-1,j,k)*deta_xi3_u(i-1,j,k)*    &
                   u(i-1,j,k) ) * rdxi1

    d_xi2_term = ( h1_xi2_v(i,j,k)*h3_xi2_v(i,j,k)*deta_xi3_v(i,j,k)*          &
                   v(i,j,k) -                                                  &
                   h1_xi2_v(i,j-1,k)*h3_xi2_v(i,j-1,k)*deta_xi3_v(i,j-1,k)*    &
                   v(i,j-1,k) ) * rdxi2

    deta_term = ( h1_p_eta(i,j,k)*h2_p_eta(i,j,k)*                             &
                  h3_p_eta(i,j,k)*deta_xi3_theta(i,j,k)*                       &
                  etadot(i,j,k) ) * rdxi3

    Jacobian  = h1_p(i,j,k)*h2_p(i,j,k)*h3_p(i,j,k)*deta_xi3(i,j,k)

    div(i,j,k) = ( d_xi1_term + d_xi2_term + deta_term ) / Jacobian

! Upper boundary
    k=k_end
    rdxi3 = 1.0/( eta_theta_levels(k) - eta_theta_levels(k-1) )

    d_xi1_term = ( h2_xi1_u(i,j,k)*h3_xi1_u(i,j,k)*deta_xi3_u(i,j,k)*          &
                   u(i,j,k) -                                                  &
                   h2_xi1_u(i-1,j,k)*h3_xi1_u(i-1,j,k)*deta_xi3_u(i-1,j,k)*    &
                   u(i-1,j,k) ) * rdxi1

    d_xi2_term = ( h1_xi2_v(i,j,k)*h3_xi2_v(i,j,k)*deta_xi3_v(i,j,k)*          &
                   v(i,j,k) -                                                  &
                   h1_xi2_v(i,j-1,k)*h3_xi2_v(i,j-1,k)*deta_xi3_v(i,j-1,k)*    &
                   v(i,j-1,k) ) * rdxi2

    deta_term = ( -h1_p_eta(i,j,k-1)*h2_p_eta(i,j,k-1)*                        &
                  h3_p_eta(i,j,k-1)*deta_xi3_theta(i,j,k-1)*                   &
                  etadot(i,j,k-1)) * rdxi3

    Jacobian  = h1_p(i,j,k)*h2_p(i,j,k)*h3_p(i,j,k)*deta_xi3(i,j,k)

    div(i,j,k) = ( d_xi1_term + d_xi2_term + deta_term ) / Jacobian

  END DO
END DO

IF (l_swap_bounds) CALL swap_bounds(div, divdims%i_len - 2*divdims%halo_i,     &
                                         divdims%j_len - 2*divdims%halo_j,     &
                                         divdims%k_len,                        &
                                         divdims%halo_i, divdims%halo_j,       &
                                         fld_type_p, swap_field_is_scalar)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE calc_div
END MODULE calc_div_mod
