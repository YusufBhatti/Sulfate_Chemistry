! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Calculate PV at all theta points.

SUBROUTINE Calc_PV_at_theta ( u, v, theta, rho,                   &
                                                            ! in
                              r_theta_levels, r_rho_levels,       &
                                                            ! in
                              r_at_u, r_at_v,                     &
                                                            ! in
                              sec_v_latitude,                     &
                                                            ! in
                              tan_v_latitude,                     &
                                                            ! in
                              sec_theta_latitude,                 &
                                                            ! in
                              f3_at_v,                            &
                                                            ! in
                              delta_lambda, delta_phi,            &
                                                            ! in
                              pv_at_theta )                 ! out

! Description:
!
!   Calculate PV at all theta points.
!
! Method:
!
!   1. Call Calc_PV to obtain PV midway E-W between v points on rho
!      levels.
!   2. Add haloes to resulting field.
!   3. Interpolate horizontally and vertically to theta points.
!      (PV at top theta level is set to PV at top rho level.)
!   4. Reset polar rows to their mean values.
!
! Note:
!
!   This routine scales poorly at high resolutions. An alternative module
!   (calc_pv_full_mod) contains Calc_PV_at_theta_opt, which has been 
!   optimized and contains more components to compute PV 
!   (e.g. vertical velocity, full Coriolis vector).
!   It is recommended to adapt the new routine (calc_pv_at_theta_opt) 
!   and deprecate this one (calc_at_theta) where possible


! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Diagnostics
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! Declarations:

USE atm_fields_bounds_mod, ONLY:                                  &
    udims, vdims, udims_s, vdims_s, tdims_s, pdims_s,             &
    udims_l, vdims_l

USE global_2d_sums_mod, ONLY: global_2d_sums
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE Field_Types, ONLY: fld_type_v
USE nlsizes_namelist_mod, ONLY: &
    global_row_length, model_levels, n_rows, row_length, rows

USE mpp_conf_mod,         ONLY: swap_field_is_scalar

USE model_domain_mod, ONLY: model_type, mt_global

IMPLICIT NONE

! Subroutine arguments:

REAL, INTENT(IN) ::                                               &

  u(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end,  &
    udims_s%k_start:udims_s%k_end),                               &

  v(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end,  &
    vdims_s%k_start:vdims_s%k_end),                               &

  theta(tdims_s%i_start:tdims_s%i_end,                            &
         tdims_s%j_start:tdims_s%j_end,                           &
         tdims_s%k_start:tdims_s%k_end),                          &

  rho(pdims_s%i_start:pdims_s%i_end,                              &
       pdims_s%j_start:pdims_s%j_end,                             &
       pdims_s%k_start:pdims_s%k_end),                            &

  r_theta_levels     ( 1 - halo_i : row_length + halo_i,          &
                       1 - halo_j : rows       + halo_j,          &
                       0          : model_levels ),               &

  r_rho_levels       ( 1 - halo_i : row_length + halo_i,          &
                       1 - halo_j : rows       + halo_j,          &
                       1          : model_levels ),               &

  r_at_u (udims_l%i_start:udims_l%i_end,                          &
          udims_l%j_start:udims_l%j_end,                          &
          udims_l%k_start:udims_l%k_end),                         &

  r_at_v (vdims_l%i_start:vdims_l%i_end,                          &
          vdims_l%j_start:vdims_l%j_end,                          &
          vdims_l%k_start:vdims_l%k_end),                         &

  sec_v_latitude (vdims_s%i_start:vdims_s%i_end,                  &
                  vdims_s%j_start:vdims_s%j_end),                 &

  tan_v_latitude(vdims%i_start:vdims%i_end,                       &
                 vdims%j_start:vdims%j_end),                      &

  sec_theta_latitude ( 1 - offx   : row_length + offx,            &
                       1 - offy   : rows       + offy ),          &

  f3_at_v (vdims_s%i_start:vdims_s%i_end,                         &
           vdims_s%j_start:vdims_s%j_end),                        &

  delta_lambda,                                                   &
  delta_phi

REAL, INTENT(OUT) ::                                              &

  pv_at_theta (row_length, rows, model_levels)

! Local variables:

INTEGER :: i, j, k
INTEGER :: ICode

REAL :: pv (udims%i_start:udims%i_end,                            &
            vdims%j_start:vdims%j_end, model_levels)

REAL :: pv_plus_haloes ( udims_s%i_start:udims_s%i_end,           &
                         vdims_s%j_start:vdims_s%j_end,           &
                         1          : model_levels )

REAL :: polar_sums (model_levels)
REAL :: polar_means(model_levels)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_PV_AT_THETA'

!- End of header ------------------------------------------------------

!----------------------------------------------------------------------
! [1]: Calculate PV midway E-W between v points on rho levels.
!----------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
! DEPENDS ON: calc_pv
CALL Calc_PV  ( u, v, theta, rho,                                 &
                                                ! in
                r_theta_levels, r_rho_levels,                     &
                                                ! in
                r_at_u, r_at_v,                                   &
                                                ! in
                sec_v_latitude, tan_v_latitude,                   &
                                                ! in
                sec_theta_latitude, f3_at_v,                      &
                                                ! in
                delta_lambda, delta_phi,                          &
                                                ! in
                row_length, rows, n_rows,                         &
                                                ! in
                model_levels,                                     &
                                                ! in
                offx, offy, halo_i, halo_j,                       &
                                                ! in
                at_extremity,                                     &
                                                ! in
                pv )                            ! out

!----------------------------------------------------------------------
! [2]: Add haloes.
!----------------------------------------------------------------------

pv_plus_haloes(:,:,:) = 0.0

DO k = 1, model_levels
  DO j = vdims%j_start,vdims%j_end
    DO i = udims%i_start,udims%i_end
      pv_plus_haloes(i,j,k) = pv(i,j,k)
    END DO
  END DO
END DO

! DEPENDS ON: swap_bounds
CALL Swap_Bounds ( pv_plus_haloes,                                &
                                                       ! inout
                   row_length, n_rows, model_levels,              &
                                                       ! in
                   offx, offy, fld_type_v, swap_field_is_scalar)   ! in

!----------------------------------------------------------------------
! [3]: Interpolate to theta points.
!----------------------------------------------------------------------

DO k = 1, model_levels - 1
  DO j = 1, rows
    DO i = 1, row_length
      pv_at_theta(i,j,k) = 0.25                                   &
                           * ( ( pv_plus_haloes(i,  j,  k)        &
                               + pv_plus_haloes(i,  j-1,k)        &
                               + pv_plus_haloes(i-1,j,  k)        &
                               + pv_plus_haloes(i-1,j-1,k)        &
                               )                                  &
                             * ( r_rho_levels  (i,j,k+1)          &
                               - r_theta_levels(i,j,k)            &
                               )                                  &
                             + ( pv_plus_haloes(i,  j,  k+1)      &
                               + pv_plus_haloes(i,  j-1,k+1)      &
                               + pv_plus_haloes(i-1,j,  k+1)      &
                               + pv_plus_haloes(i-1,j-1,k+1)      &
                               )                                  &
                             * ( r_theta_levels(i,j,k)            &
                               - r_rho_levels  (i,j,k)            &
                             ) )                                  &
                           / ( r_rho_levels(i,j,k+1)              &
                             - r_rho_levels(i,j,k)                &
                             )
    END DO
  END DO
END DO

! Set PV at top theta level equal to PV at top rho level.
k = model_levels
DO j = 1, rows
  DO i = 1, row_length
    pv_at_theta(i,j,k) = 0.25                                     &
                       * ( pv_plus_haloes(i,  j,  k)              &
                         + pv_plus_haloes(i,  j-1,k)              &
                         + pv_plus_haloes(i-1,j,  k)              &
                         + pv_plus_haloes(i-1,j-1,k)              &
                         )
  END DO
END DO


IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE Calc_PV_at_theta
