! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose:
! This module contains two subroutines
! 
! calc_pv_new_or_old: Chooses between the default method or "new method",
! the former calls calc_pv_at_theta and latter calc_PV_at_theta_opt described
! below.
! 
! Calc_PV_at_theta_opt: Calculates the components of
! vorticity and grad(theta) separately and combines them.
! The horizontal components of vorticity and grad(theta) are
! evaluated on theta levels and the vertical components are evaluated on rho
! levels due to the staggering of grid points. So the vertical components
! are combined then interpolated in the vertical before adding to the
! horizontal components. This routine is employed instead of the 
! calc_pv_at_theta as it has been optimized and includes further components
! (e.g. vertical velocity, full coriolis vector)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics diagnostics
!
! Code Description: Chooses between the new and default PV calculations. 
!                   Contains the routine that computes PV under the new method
!                   which includes the vertical velocity component and
!                   the 3 components from coriolis vector.

MODULE calc_PV_full_mod

USE UM_ParParams,          ONLY: PNorth, PSouth
USE UM_ParVars,            ONLY: at_extremity

! Grid parameters
USE level_heights_mod,     ONLY: r_at_u, r_at_v, r_theta_levels,        &
                                 r_rho_levels


USE cderived_mod,          ONLY: delta_lambda, delta_phi

USE atm_fields_bounds_mod, ONLY: array_dims, tdims, pdims,              &
                                 udims_s, vdims_s, wdims_s, tdims_s,    &
                                 pdims_s, udims_l, vdims_l, tdims_l

! Use h_i factors rather than (r cos phi, r,1) as these are designed for
! derivates and got the same haloes as u,v,w
USE metric_terms_mod,      ONLY: h1_p_eta,h2_p_eta,                     &
                                 h1_p, h2_p, h1_xi1_u, h1_xi2_v,        &
                                 h2_xi1_u, h2_xi2_v

! Use horiz_grid_mod to weight vertical levels
USE horiz_grid_mod, ONLY: intw_rho2w

! Run type
USE model_domain_mod,      ONLY: model_type, mt_global

! Dr Hook
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE
CHARACTER(LEN=*), PRIVATE, PARAMETER :: ModuleName = 'CALC_PV_FULL_MOD'

PRIVATE
PUBLIC :: calc_pv_new_or_old,calc_pv_at_theta_opt

CONTAINS
  !----------------------------------------------------------------------
SUBROUTINE calc_pv_new_or_old(l_coriolis, u, v, theta, rho, pv)

! Chooses between the old and new method to compute PV

USE atm_fields_mod,             ONLY: w
USE dyn_coriolis_mod,           ONLY: f1_at_v, f2_at_u, f3_at_v

USE trignometric_mod,           ONLY: sec_v_latitude, tan_v_latitude,   &
                                      sec_theta_latitude


USE free_tracers_inputs_mod,    ONLY: l_calc_pv_full

IMPLICIT NONE
!-Input-Variables
LOGICAL, INTENT(IN) ::                                                  &
  l_coriolis
  ! Set coriolis terms to its value if True,
  ! or to zero if False in case we are computing the friction term 
  ! from wind increments


REAL, INTENT(IN) ::                                                     &
  ! Wind fields
  u(udims_s%i_start:udims_s%i_end, udims_s%j_start:udims_s%j_end,       &
    udims_s%k_start:udims_s%k_end)                                      &
, v(vdims_s%i_start:vdims_s%i_end, vdims_s%j_start:vdims_s%j_end,       &
    vdims_s%k_start:vdims_s%k_end)                                      &

  ! Potential temperature
, theta(tdims_s%i_start:tdims_s%i_end, tdims_s%j_start:tdims_s%j_end,   &
        tdims_s%k_start:tdims_s%k_end)                                  &

  ! Density
, rho(pdims_s%i_start:pdims_s%i_end, pdims_s%j_start:pdims_s%j_end,     &
      pdims_s%k_start:pdims_s%k_end)

REAL, INTENT(OUT) ::                                                    &
! PV at theta points
  pv (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end,             &
      1:tdims%k_end)

REAL ::                                                                 &
  zero_f_v (vdims_s%i_start:vdims_s%i_end,                              &
            vdims_s%j_start:vdims_s%j_end)                              &

, zero_f_u (udims_s%i_start:udims_s%i_end,                              &
            udims_s%j_start:udims_s%j_end)

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL   (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_PV_NEW_OR_OLD'

!-End-Header---------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (l_coriolis) THEN

  ! Call any of the two routines available to compute PV
  IF (l_calc_pv_full) THEN

    CALL Calc_PV_at_theta_opt (u, v, w, theta, rho,                     &
                               f1_at_v, f2_at_u, f3_at_v,               &
                               pv)
  ELSE
    CALL calc_pv_at_theta(u, v, theta, rho,                             &
                          r_theta_levels, r_rho_levels,                 &
                          r_at_u, r_at_v,                               &
                          sec_v_latitude, tan_v_latitude,               &
                          sec_theta_latitude, f3_at_v,                  &
                          delta_lambda, delta_phi,                      &
                          pv)
  END IF
ELSE

  zero_f_u = 0.0
  zero_f_v = 0.0

  ! Call any of the two routines available to compute PV
  IF (l_calc_pv_full) THEN
    CALL Calc_PV_at_theta_opt (u, v, w, theta, rho,                     &
                               zero_f_v, zero_f_u, zero_f_v,            &
                               pv)
  ELSE
    CALL calc_pv_at_theta(u, v, theta, rho,                             &
                          r_theta_levels, r_rho_levels,                 &
                          r_at_u, r_at_v,                               &
                          sec_v_latitude, tan_v_latitude,               &
                          sec_theta_latitude, zero_f_v,                 &
                          delta_lambda, delta_phi,                      &
                          pv)
  END IF

END IF ! end if over coriolis

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE calc_pv_new_or_old

  !----------------------------------------------------------------------
SUBROUTINE Calc_PV_at_theta_opt (u, v, w, theta, rho,                   &
                                 f1_at_v, f2_at_u, f3_at_v,             &
                                 PV)

! Horizontal Staggering
! -Vertical derivatives are calculated at midpoints so appear on alternate
! -levels. Horizontal derivatives use centred differences between -1/+1
! -indices so appear at the same position as the variable
!
! Rho Levels
! _______________ _______________ _______________
!|               |               |               |
!|       v       |               |               |
!|_______________|_______________|_______________|
!|               |               |               |
!|               |               |               |
!|_______________|_______________|_______________|
!|               |               |               |
!|rho, dtheta/dz |               |       u       |
!|_______________|_______________|_______________|

! Theta Levels
! _______________ _______________ _______________
!|               |               |               |
!|     dv/dz     |               |               |
!|_______________|_______________|_______________|
!|               |               |               |
!|               |               |               |
!|_______________|_______________|_______________|
!|               |               |               |
!|   theta, w    |               |     du/dz     |
!|_______________|_______________|_______________|

IMPLICIT NONE

!-Input-Variables
REAL, INTENT(IN) ::                                                     &
  ! Wind fields
  u(udims_s%i_start:udims_s%i_end, udims_s%j_start:udims_s%j_end,       &
    udims_s%k_start:udims_s%k_end)                                      &
, v(vdims_s%i_start:vdims_s%i_end, vdims_s%j_start:vdims_s%j_end,       &
    vdims_s%k_start:vdims_s%k_end)                                      &
, w(wdims_s%i_start:wdims_s%i_end, wdims_s%j_start:wdims_s%j_end,       &
    wdims_s%k_start:wdims_s%k_end)                                      &

  ! Potential temperature
, theta(tdims_s%i_start:tdims_s%i_end, tdims_s%j_start:tdims_s%j_end,   &
        tdims_s%k_start:tdims_s%k_end)                                  &

  ! Density
, rho(pdims_s%i_start:pdims_s%i_end, pdims_s%j_start:pdims_s%j_end,     &
      pdims_s%k_start:pdims_s%k_end)                                    &

  ! components of coriolis force at theta (rho) points
, f1_at_v(vdims_s%i_start:vdims_s%i_end, vdims_s%j_start:vdims_s%j_end) &
, f2_at_u(udims_s%i_start:udims_s%i_end, udims_s%j_start:udims_s%j_end) &
, f3_at_v(vdims_s%i_start:vdims_s%i_end, vdims_s%j_start:vdims_s%j_end)


REAL, INTENT(OUT) ::                                                    &
! PV at theta points
  PV (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end,             &
      1:tdims%k_end)

!-Local-Variables
! Loop counters
INTEGER :: i, j, k

REAL ::                                                                 &
  ! Components of vorticity
  lambda_term(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end,     &
              1:tdims%k_end)                                            &
, phi_term   (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end,     &
              1:tdims%k_end)                                            &
, r_term     (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end,     &
              pdims%k_start:pdims%k_end)                                &

  ! Grad(theta) components
, dtheta_dlambda(tdims_s%i_start+1:tdims_s%i_end-1,                     &
                 tdims_s%j_start:tdims_s%j_end,                         &
                 tdims_s%k_start:tdims_s%k_end)                         &
, dtheta_dphi   (tdims_s%i_start:tdims_s%i_end,                         &
                 tdims_s%j_start+1:tdims_s%j_end-1,                     &
                 tdims_s%k_start:tdims_s%k_end)                         &
, dtheta_dr     (tdims_s%i_start:tdims_s%i_end,                         &
                 tdims_s%j_start:tdims_s%j_end,                         &
                 tdims_s%k_start:tdims_s%k_end)

REAL ::                                                                 &
  ! Components of PV
  PV_r(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end,            &
       pdims%k_start:pdims%k_end)                                       &
, PV_lambda, PV_phi                                                     &
  ! Z-component of PV averaged to theta levels
, PV_r_ave                                                              &
  ! Density averaged to theta levels
, rho_theta

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL   (KIND=jprb)            :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_PV_AT_THETA_OPT'

!-End-Header---------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Calculate vorticity
! dw/dphi - d(h2v)/dr (at theta points)
lambda_term = curl_lambda(v, w)
! d(h1u)/dr - dw/dlambda (at theta points)
phi_term    = curl_phi(u, w)
! d(h2v)/dlambda - d(h1u)/dphi (at rho points)
r_term      = curl_r(u, v)

! Calculate grad(theta) components
! dtheta/dlambda (at theta points)
dtheta_dlambda = d_dlambda(theta, tdims_s, h1_p_eta, tdims_s, tdims_s)
! dtheta/dphi (at theta points)
dtheta_dphi    = d_dphi(theta, tdims_s, h2_p_eta, tdims_s, tdims_s)
! dtheta/dr (at rho points)
dtheta_dr      = d_dr(theta, r_theta_levels, tdims_s, tdims_l, 1)

! Calculate vertical component of PV since it is still at rho-levels
! and needs to be averaged
DO k = pdims%k_start, pdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      PV_r(i,j,k) = (r_term(i,j,k) + 0.5 * (f3_at_v(i,j) +              &
                                            f3_at_v(i,j-1))) *          &
                    dtheta_dr(i,j,k)
    END DO
  END DO
END DO

DO k = 1, tdims%k_end-1
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      ! Calculate horizontal PV components
      PV_lambda = (lambda_term(i,j,k) +                                   &
              0.5 * (f1_at_v(i,j) + f1_at_v(i,j-1))) * dtheta_dlambda(i,j,k)
      PV_phi = (phi_term(i,j,k) +                                         &
              0.5 * (f2_at_u(i,j) + f2_at_u(i-1,j))) * dtheta_dphi(i,j,k)

      ! Average vertical PV component to theta levels
      PV_r_ave = PV_r(i,j,k) * intw_rho2w(k,2) +                        &
                 PV_r(i,j,k+1) *intw_rho2w(k,1)

      ! Average density to theta points
      rho_theta = rho(i,j,k) * intw_rho2w(k,2) +                        &
                  rho(i,j,k+1) *intw_rho2w(k,1)

      ! Combine to calculate full PV
      PV(i,j,k) = (PV_lambda + PV_phi + PV_r_ave) / rho_theta
    END DO
  END DO
END DO

! Set PV at the top level equal to the level below since it is undefined
DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end
    PV(i,j,tdims%k_end) = PV(i,j,tdims%k_end-1)
  END DO
END DO

! PV is unresolved at the poles so set to PV at the next equatorward level
IF (model_type == mt_global) THEN
  IF (at_extremity(PSouth)) THEN
    DO k = 1, tdims%k_end
      DO i = tdims%i_start, tdims%i_end
        PV(i,tdims%j_start,k) = PV(i,tdims%j_start+1,k)
      END DO
    END DO
  END IF

  IF (at_extremity(PNorth)) THEN
    DO k = 1, tdims%k_end
      DO i = tdims%i_start, tdims%i_end
        PV(i,tdims%j_end,k) = PV(i,tdims%j_end-1,k)
      END DO
    END DO
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE calc_pv_at_theta_opt
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
FUNCTION curl_lambda(v, w)
! Calculate the lambda-component of curl at theta points assuming standard
! grid staggering for v and w

IMPLICIT NONE

!-Input-Variables
REAL, INTENT(IN) ::                                                     &
  v(vdims_s%i_start:vdims_s%i_end, vdims_s%j_start:vdims_s%j_end,       &
    vdims_s%k_start:vdims_s%k_end)                                      &
, w(wdims_s%i_start:wdims_s%i_end, wdims_s%j_start:wdims_s%j_end,       &
    wdims_s%k_start:wdims_s%k_end)

! Function output
REAL ::                                                                 &
  curl_lambda(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end,     &
              1:tdims%k_end)
!-Local-Variables
INTEGER :: i, j, k

REAL ::                                                                 &
  h2v  (vdims_s%i_start:vdims_s%i_end, vdims_s%j_start:vdims_s%j_end,   &
        vdims_s%k_start:vdims_s%k_end)                                  &
, dv_dr(vdims_s%i_start:vdims_s%i_end, vdims_s%j_start:vdims_s%j_end,   &
        vdims_s%k_start:vdims_s%k_end)                                  &
, dw_dphi(tdims_s%i_start:tdims_s%i_end,                                &
          tdims_s%j_start+1:tdims_s%j_end-1,                            &
          tdims_s%k_start:tdims_s%k_end)                                &
, dv_dr_ave
!-End-Header------------------------------------------------------------

! Compute "rv" for dv_dr
DO k = vdims_s%k_start, vdims_s%k_end
  DO j = vdims_s%j_start, vdims_s%j_end
    DO i = vdims_s%i_start, vdims_s%i_end
      h2v(i,j,k)=h2_xi2_v(i,j,k)*v(i,j,k)
    END DO
  END DO
END DO

! Calculate d(h2v)_dr at v-points on theta levels
dv_dr = d_dr(h2v, r_at_v, vdims_s, vdims_l, 0)
! Calculate dw_dphi at theta
dw_dphi = d_dphi(w, wdims_s, h2_p_eta, tdims_s, tdims_s)

! Calculate curl_lambda
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      ! Average to theta points
      dv_dr_ave = 0.5 * (dv_dr(i,j,k) + dv_dr(i,j-1,k))/h2_p_eta(i,j,k)
      ! Calculate output
      curl_lambda(i,j,k) = dw_dphi(i,j,k) - dv_dr_ave
    END DO
  END DO
END DO
END FUNCTION curl_lambda
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
FUNCTION curl_phi(u, w)
! Calculate the phi-component of curl at theta points assuming standard grid
! staggering for u and w

IMPLICIT NONE

!-Input-Variables
REAL, INTENT(IN) ::                                                     &
  u(udims_s%i_start:udims_s%i_end, udims_s%j_start:udims_s%j_end,       &
    udims_s%k_start:udims_s%k_end)                                      &
, w(wdims_s%i_start:wdims_s%i_end, wdims_s%j_start:wdims_s%j_end,       &
    wdims_s%k_start:wdims_s%k_end)
! Function output
REAL ::                                                                 &
  curl_phi(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end,        &
           1:tdims%k_end)
!-Local-Variables
INTEGER :: i, j, k
REAL ::                                                                 &
  h1u  (udims_s%i_start:udims_s%i_end, udims_s%j_start:udims_s%j_end,   &
        udims_s%k_start:udims_s%k_end)                                  &
, du_dr(udims_s%i_start:udims_s%i_end, udims_s%j_start:udims_s%j_end,   &
        udims_s%k_start:udims_s%k_end)                                  &
, dw_dlambda(tdims_s%i_start+1:tdims_s%i_end-1,                         &
             tdims_s%j_start:tdims_s%j_end,                             &
             tdims_s%k_start:tdims_s%k_end)                             &
, du_dr_ave
!-End-Header------------------------------------------------------------

! Compute "h1u" for du_dr
DO k = udims_s%k_start, udims_s%k_end
  DO j = udims_s%j_start, udims_s%j_end
    DO i = udims_s%i_start, udims_s%i_end
      h1u(i,j,k)=h1_xi1_u(i,j,k)*u(i,j,k)
    END DO
  END DO
END DO

! Calculate du_dr at u-points on theta levels
du_dr = d_dr(h1u, r_at_u, udims_s, udims_l, 0)
! Calculate dw_dlambda at theta points
dw_dlambda = d_dlambda(w, wdims_s, h1_p_eta, tdims_s, tdims_s)

! Calculate curl_phi
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      ! Average to theta points
      du_dr_ave = 0.5 * (du_dr(i,j,k) + du_dr(i-1,j,k))/h1_p_eta(i,j,k)
      ! Calculate output
      curl_phi(i,j,k) = du_dr_ave - dw_dlambda(i,j,k)
    END DO
  END DO
END DO
END FUNCTION curl_phi
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
FUNCTION curl_r(u, v)
! Calculate the vertical component of curl at rho points

IMPLICIT NONE

!-Input-Variables
REAL, INTENT(IN) ::                                                     &
  u(udims_s%i_start:udims_s%i_end, udims_s%j_start:udims_s%j_end,       &
    udims_s%k_start:udims_s%k_end)                                      &
, v(vdims_s%i_start:vdims_s%i_end, vdims_s%j_start:vdims_s%j_end,       &
    vdims_s%k_start:vdims_s%k_end)

! Function output
REAL ::                                                                 &
  curl_r(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end,          &
         pdims%k_start:pdims%k_end)
!-Local-Variables
INTEGER :: i, j, k

REAL :: du_dphi_ave, dv_dlambda_ave
REAL ::                                                                 &
  h1u  (udims_s%i_start:udims_s%i_end, udims_s%j_start:udims_s%j_end,   &
        udims_s%k_start:udims_s%k_end)                                  &
, h2v  (vdims_s%i_start:vdims_s%i_end, vdims_s%j_start:vdims_s%j_end,   &
        vdims_s%k_start:vdims_s%k_end)                                  &
, du_dphi(udims_s%i_start:udims_s%i_end,                                &
        udims_s%j_start+1:udims_s%j_end-1,                              &
        udims_s%k_start:udims_s%k_end)                                  &
, dv_dlambda(vdims_s%i_start+1:vdims_s%i_end-1,                         &
             vdims_s%j_start:vdims_s%j_end,                             &
             vdims_s%k_start:vdims_s%k_end)

!-End-Header------------------------------------------------------------

! Compute "h1u" for du_dphi
DO k = udims_s%k_start, udims_s%k_end
  DO j = udims_s%j_start, udims_s%j_end
    DO i = udims_s%i_start, udims_s%i_end
      h1u(i,j,k)=h1_xi1_u(i,j,k)*u(i,j,k)
    END DO
  END DO
END DO

! Compute "rv" for dv_dlambda
DO k = vdims_s%k_start, vdims_s%k_end
  DO j = vdims_s%j_start, vdims_s%j_end
    DO i = vdims_s%i_start, vdims_s%i_end
      h2v(i,j,k)=h2_xi2_v(i,j,k)*v(i,j,k)
    END DO
  END DO
END DO

! Calculate du_dphi at u-points
du_dphi = d_dphi(h1u, udims_s, h2_xi1_u, udims_s, udims_s)
! Calculate dv_dlambda at v-points
dv_dlambda = d_dlambda(h2v, vdims_s, h1_xi2_v, vdims_s, vdims_s)

! Calculate curl_r
DO k = pdims%k_start, pdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      ! Average to rho points
      du_dphi_ave = 0.5 * (du_dphi(i,j,k) + du_dphi(i-1,j,k))/          &
                           h1_p(i,j,k)
      dv_dlambda_ave = 0.5 * (dv_dlambda(i,j,k) + dv_dlambda(i,j-1,k))/ &
                           h2_p(i,j,k)
      ! Calculate output
      curl_r(i,j,k) = dv_dlambda_ave - du_dphi_ave
    END DO
  END DO
END DO
END FUNCTION curl_r
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
FUNCTION d_dlambda(array, dims_a, h1, dims_h1, loop)
! Calculates the derivative in the lambda-direction using centred 
! differencing. Therefore the output array will have a lambda-dimension 
! reduced by two

IMPLICIT NONE

!-Input-Variables
! Array dimensions for inputs and indices to loop over
TYPE(array_dims), INTENT(IN) :: dims_a, dims_h1, loop

REAL, INTENT(IN) ::                                                     &
  array(dims_a%i_start:dims_a%i_end, dims_a%j_start:dims_a%j_end,       &
        dims_a%k_start:dims_a%k_end)                                    &
, h1   (dims_h1%i_start:dims_h1%i_end, dims_h1%j_start:dims_h1%j_end,   &
        dims_h1%k_start:dims_h1%k_end)
! Function output
REAL ::                                                                 &
  d_dlambda(loop%i_start+1:loop%i_end-1, loop%j_start:loop%j_end,       &
            loop%k_start:loop%k_end)
!-Local-Variables
! Loop counters
INTEGER :: i, j, k
!-End-Header------------------------------------------------------------

DO k = loop%k_start, loop%k_end
  DO j = loop%j_start, loop%j_end
    DO i = loop%i_start+1, loop%i_end-1
      d_dlambda(i,j,k) = (array(i+1,j,k) - array(i-1,j,k)) /            &
                          (h1(i,j,k) * 2 * delta_lambda)
    END DO
  END DO
END DO
END FUNCTION d_dlambda
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
FUNCTION d_dphi(array, dims_a, h2, dims_h2, loop)
! Calculates the derivative in the phi-direction using centred differencing
! Therefore the output array will have a phi-dimension reduced by two

IMPLICIT NONE

!-Input-Variables
! Array dimensions with small and large haloes
TYPE(array_dims), INTENT(IN) :: dims_a, dims_h2, loop

REAL, INTENT(IN) ::                                                     &
  array(dims_a%i_start:dims_a%i_end, dims_a%j_start:dims_a%j_end,       &
        dims_a%k_start:dims_a%k_end)                                    &
, h2   (dims_h2%i_start:dims_h2%i_end, dims_h2%j_start:dims_h2%j_end,   &
        dims_h2%k_start:dims_h2%k_end)
! Function output
REAL ::                                                                 &
  d_dphi(loop%i_start:loop%i_end, loop%j_start+1:loop%j_end-1,          &
         loop%k_start:loop%k_end)
!-Local-Variables
! Loop counters
INTEGER :: i, j, k
!-End-Header------------------------------------------------------------
DO k = loop%k_start, loop%k_end
  DO j = loop%j_start+1, loop%j_end-1
    DO i = loop%i_start, loop%i_end
      d_dphi(i,j,k) = (array(i,j+1,k) - array(i,j-1,k)) /               &
                      (h2(i,j,k) * 2 * delta_phi)
    END DO
  END DO
END DO
END FUNCTION d_dphi
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
FUNCTION d_dr(array, r, dims_s, dims_l, offset)
! Calculate vertical derivatives at midpoints in the vertical

IMPLICIT NONE

! For differences on theta levels: dlambda_d(rho)_i = lambda_k+1 - lambda_k
! For differences on rho levels: dlambda_d(theta)_i = lambda_k - lambda_k-1
! This results in either the bottom or top level being unspecified

!-Input-Variables
! Array dimensions with small and large haloes
TYPE(array_dims), INTENT(IN) :: dims_s, dims_l

! Offset on array indices
! Output on theta levels (input on rho levels): offset=0
! Output on rho levels (input on theta levels): offset=1
INTEGER, INTENT(IN) :: offset

REAL, INTENT(IN) ::                                                     &
  array(dims_s%i_start:dims_s%i_end, dims_s%j_start:dims_s%j_end,       &
        dims_s%k_start:dims_s%k_end)                                    &
, r(dims_l%i_start:dims_l%i_end, dims_l%j_start:dims_l%j_end,           &
    dims_l%k_start:dims_l%k_end)
! Function output
REAL ::                                                                 &
  d_dr(dims_s%i_start:dims_s%i_end, dims_s%j_start:dims_s%j_end,        &
       dims_s%k_start:dims_s%k_end)
!-Local-Variables
! Loop counters
INTEGER :: i, j, k
!-End-Header------------------------------------------------------------

DO k = 1, dims_s%k_end-1
  DO j = dims_s%j_start, dims_s%j_end
    DO i = dims_s%i_start, dims_s%i_end
      d_dr(i,j,k+offset) = (array(i,j,k+1) - array(i,j,k)) /            &
                           (r(i,j,k+1) - r(i,j,k) )
    END DO
  END DO
END DO

IF (offset==0) THEN
  DO j = dims_s%j_start, dims_s%j_end
    DO i = dims_s%i_start, dims_s%i_end
      ! Output is on theta levels so top level is missing
      ! Assume continuous gradient at the top boundary
      d_dr(i,j,dims_s%k_end) = d_dr(i,j,dims_s%k_end-1)
    END DO
  END DO
ELSE IF (offset==1) THEN
  ! Output is on rho levels so bottom level is missing
  ! Assume continuous gradient at the bottom boundary
  DO j = dims_s%j_start, dims_s%j_end
    DO i = dims_s%i_start, dims_s%i_end
      d_dr(i,j,1) = d_dr(i,j,2)
    END DO
  END DO
END IF
END FUNCTION d_dr
!-----------------------------------------------------------------------
END MODULE calc_pv_full_mod
