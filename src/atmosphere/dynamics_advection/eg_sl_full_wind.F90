! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_sl_full_wind_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_SL_FULL_WIND_MOD'

CONTAINS
SUBROUTINE eg_sl_full_wind(g_i_pe,depart_scheme,                               &
                depart_order, high_order_scheme, monotone_scheme,              &
                depart_high_order_scheme, depart_monotone_scheme,              &
                first_constant_r_rho_level,                                    &
                interp_vertical_search_tol, check_bottom_levels,               &
                l_shallow, l_rk_dps,                                           &
                l_high, l_mono,                                                &
                l_depart_high, l_depart_mono, lam_max_cfl,                     &
                etadot, u_np1, v_np1, w_np1, etadot_np1,                       &
                u, v, w, r_u, r_v, r_w, r_u_d, r_v_d, r_w_d,                   &
                error_code )


USE parkind1,          ONLY: jpim, jprb       !DrHook
USE yomhook,           ONLY: lhook, dr_hook   !DrHook

USE nlsizes_namelist_mod,  ONLY: global_row_length,                            &
                                 row_length, rows, n_rows, model_levels

USE um_parvars,            ONLY: at_extremity,                                 &
                                 halo_i, halo_j, offx, offy

USE timestep_mod,      ONLY: timestep
USE level_heights_mod, ONLY: eta_theta_levels, eta_rho_levels,                 &
                             xi3_at_theta=>r_theta_levels,                     &
                             xi3_at_rho=>r_rho_levels,                         &
                             xi3_at_u=>r_at_u, xi3_at_v=>r_at_v,               &
                             xi3_at_u_w=>r_at_u_w,                             &
                             xi3_at_v_w=>r_at_v_w
USE eg_v_at_poles_mod
USE eg_sl_wind_u_mod
USE eg_sl_wind_v_mod
USE eg_sl_wind_w_mod
USE atm_fields_bounds_mod

USE UM_ParParams
USE departure_pts_mod
USE halo_exchange, ONLY: swap_bounds
USE mpp_conf_mod,  ONLY: swap_field_is_scalar,swap_field_is_vector

USE model_domain_mod, ONLY: model_type, mt_global

USE atm_step_local, ONLY: L_print_L2norms, cycleno
USE sl_wind_norm_mod, ONLY: sl_wind_norm
USE umPrintMgr, ONLY: umMessage, umPrint
USE turb_diff_mod, ONLY: norm_lev_start, norm_lev_end
USE lam_config_inputs_mod, ONLY: n_rims_to_do

IMPLICIT NONE
!
! Description:
!   Find departure points on u,v,w-grid
!   and apply rotation on timelevel n
!   quantities R_u, R_v, R_w.
!
!
! Method: ENDGame formulation version 1.01,
!         section 7.2.
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_SL_FULL_WIND'

INTEGER, INTENT(IN) :: first_constant_r_rho_level

! MPP options

INTEGER, INTENT(IN) ::                                                &
  g_i_pe(1-halo_i:global_row_length+halo_i),                          &
                     ! processor on my processor-row
                     ! holding a given value in i direction
  lam_max_cfl(2)     ! Max CFL for a LAM allowed near the
                     ! boundaries

! Loop index bounds for arrays defined on p, u, v points respectively

! Integer parameters for advection

INTEGER, INTENT(IN) ::                                                &
  high_order_scheme,                                                  &
                     ! a code saying which high order scheme to
                     ! use. 1 = tensor tri-cubic lagrange order
                     ! (j,i,k) no other options available at
                     ! present.
  monotone_scheme,                                                    &
                     ! a code saying which monotone scheme to use.
                     ! 1 = tri-linear order (j,i,k)
                     ! no other options available at present.
  depart_scheme,                                                      &
                     ! code saying which departure point scheme to
                     ! use.
  depart_order,                                                       &
                     ! for the chosen departure point scheme how
                     ! many iterations/terms to use.
  depart_high_order_scheme,                                           &
                     ! code choosing high order
                     ! interpolation scheme used in
                     ! Departure_point routine
  depart_monotone_scheme,                                             &
                     ! code choosing monotone
                     ! interpolation scheme used in
                     ! Departure_point routine
  interp_vertical_search_tol,                                         &
                      ! used in interpolation code.
  check_bottom_levels ! used in interpolation code, and is
                      ! the number of levels to check to see
                      ! if the departure point lies inside the
                      ! orography.

! Logical switches for advection

LOGICAL, INTENT(IN) ::                                                &
  l_high,                                                             &
                   ! True, if high order interpolation required.
  l_mono,                                                             &
                   ! True, if interpolation required to be monotone.
  l_depart_high,                                                      &
                   ! True if high order interpolation scheme to be
                   ! used in Departure scheme
  l_depart_mono
                   ! True if monotone interpolation scheme to be
                   ! used in Departure scheme

LOGICAL, INTENT(IN) :: l_shallow, l_rk_dps

REAL, INTENT(INOUT) ::                                                &
  etadot(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),     &
  u(-halo_i:row_length-1+halo_i,1-halo_j:rows+halo_j,                 &
         model_levels),                                               &
  v(1-halo_i:row_length+halo_i,-halo_j:n_rows-1+halo_j,               &
        model_levels),                                                &
  w(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,                 &
         0:model_levels),                                             &
  u_np1(-offx:row_length-1+offx,1-offy:rows+offy, model_levels),      &
  v_np1(1-offx:row_length+offx,-offy:n_rows-1+offy, model_levels),    &
  w_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
  etadot_np1(1-offx:row_length+offx,1-offy:rows+offy,                 &
             0:model_levels)

! Timelevel n arrival point quantities

REAL, INTENT(INOUT) ::                                                &
  r_u(-offx:row_length-1+offx,1-offy:rows+offy,model_levels),         &
  r_v(1-offx:row_length+offx,-offy:n_rows-1+offy,model_levels),       &
  r_w(row_length,rows,0:model_levels)


INTEGER :: error_code   ! Non-zero on exit if error detected.

! Timelevel n departure point quantities

REAL, INTENT(OUT) ::                                                  &
  r_u_d(-offx:row_length-1+offx,1-offy:rows+offy,model_levels),       &
  r_v_d(1-offx:row_length+offx,-offy:n_rows-1+offy,model_levels),     &
  r_w_d(row_length,rows,0:model_levels)

! Halo-ed copies of R_u, R_v, R_w for interpolation

REAL ::                                                               &
  rwork_u(-halo_i:row_length-1+halo_i,1-halo_j:rows+halo_j,           &
           model_levels),                                             &
  rwork_v(1-halo_i:row_length+halo_i,-halo_j:n_rows-1+halo_j,         &
           model_levels),                                             &
  rwork_w(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,            &
            0:model_levels)

! Local variables

INTEGER :: i,j,k

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (model_type == mt_global) THEN
  IF ( at_extremity(psouth) ) THEN

    CALL eg_v_at_poles(r_u,r_v, 1.0, udims%j_start, vdims%j_start,&
                       udims_s,vdims_s)

  END IF
  IF ( at_extremity(pnorth) ) THEN

    CALL eg_v_at_poles(r_u,r_v, -1.0, udims%j_end, vdims%j_end,&
                       udims_s,vdims_s)

  END IF
END IF

! set workspace

k = 0
DO j=pdims%j_start, pdims%j_end
  DO i=pdims%i_start, pdims%i_end
    rwork_w(i,j,k) = r_w(i,j,k)
  END DO
END DO

!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                &
!$OMP& SHARED(model_levels,udims,vdims,pdims,r_u,r_v,r_w,     &
!$OMP& rwork_u,rwork_v,rwork_w) SCHEDULE(STATIC)
DO k=1, model_levels

  DO j=udims%j_start, udims%j_end
    DO i=udims%i_start, udims%i_end
      rwork_u(i,j,k) = r_u(i,j,k)
    END DO
  END DO

  DO j=vdims%j_start, vdims%j_end
    DO i=vdims%i_start, vdims%i_end
      rwork_v(i,j,k) = r_v(i,j,k)
    END DO
  END DO

  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end
      rwork_w(i,j,k) = r_w(i,j,k)
    END DO
  END DO

END DO
!$OMP END PARALLEL DO

CALL swap_bounds(rwork_u,                                              &
                 udims_l%i_len - 2*udims_l%halo_i,                     &
                 udims_l%j_len - 2*udims_l%halo_j,                     &
                 udims_l%k_len,                                        &
                 udims_l%halo_i, udims_l%halo_j,                       &
                 fld_type_u,swap_field_is_vector)
CALL swap_bounds(rwork_v,                                              &
                 vdims_l%i_len - 2*vdims_l%halo_i,                     &
                 vdims_l%j_len - 2*vdims_l%halo_j,                     &
                 vdims_l%k_len,                                        &
                 vdims_l%halo_i, vdims_l%halo_j,                       &
                 fld_type_v,swap_field_is_vector)
CALL swap_bounds(rwork_w,                                              &
                 wdims_l%i_len - 2*wdims_l%halo_i,                     &
                 wdims_l%j_len - 2*wdims_l%halo_j,                     &
                 wdims_l%k_len,                                        &
                 wdims_l%halo_i, wdims_l%halo_j,                       &
                 fld_type_p, swap_field_is_scalar)

CALL eg_sl_wind_w(g_i_pe, depart_scheme, l_rk_dps,                    &
                depart_order, high_order_scheme, monotone_scheme,     &
                depart_high_order_scheme, depart_monotone_scheme,     &
                first_constant_r_rho_level,                           &
                interp_vertical_search_tol, check_bottom_levels,      &
                l_high, l_mono,                                       &
                l_depart_high, l_depart_mono,                         &
                l_shallow, lam_max_cfl,                               &
                etadot, u_np1, v_np1, w_np1, etadot_np1,              &
                u, v, w, rwork_u, rwork_v, rwork_w, r_w_d,            &
                error_code )

CALL eg_sl_wind_u(g_i_pe, depart_scheme, l_rk_dps,                    &
                depart_order, high_order_scheme, monotone_scheme,     &
                depart_high_order_scheme,depart_monotone_scheme,      &
                first_constant_r_rho_level,                           &
                interp_vertical_search_tol, check_bottom_levels,      &
                l_high, l_mono,                                       &
                l_depart_high, l_depart_mono,                         &
                l_shallow, lam_max_cfl,                               &
                etadot, u_np1, v_np1, w_np1, etadot_np1,              &
                u, v, w, rwork_u, rwork_v, rwork_w, r_u_d,            &
                error_code )

CALL eg_sl_wind_v(g_i_pe,depart_scheme, l_rk_dps,                     &
                depart_order, high_order_scheme, monotone_scheme,     &
                depart_high_order_scheme,depart_monotone_scheme,      &
                first_constant_r_rho_level,                           &
                interp_vertical_search_tol, check_bottom_levels,      &
                l_high, l_mono,                                       &
                l_depart_high, l_depart_mono,                         &
                l_shallow, lam_max_cfl,                               &
                etadot, u_np1, v_np1, w_np1, etadot_np1,              &
                u, v, w, rwork_u, rwork_v, rwork_w, r_v_d,error_code )


CALL swap_bounds(r_u_d,                                                &
                 udims_s%i_len - 2*udims_s%halo_i,                     &
                 udims_s%j_len - 2*udims_s%halo_j,                     &
                 udims_s%k_len,                                        &
                 udims_s%halo_i, udims_s%halo_j,                       &
                 fld_type_u,swap_field_is_vector, do_east_arg=.TRUE.)
CALL swap_bounds(r_v_d,                                                &
                 vdims_s%i_len - 2*vdims_s%halo_i,                     &
                 vdims_s%j_len - 2*vdims_s%halo_j,                     &
                 vdims_s%k_len,                                        &
                 vdims_s%halo_i, vdims_s%halo_j,                       &
                 fld_type_v,swap_field_is_vector, do_north_arg=.TRUE.)

 IF ( L_print_L2norms ) THEN
    WRITE(umMessage,'(A, I2, A)') ' ** cycleno = ',  cycleno             &
                            , ' **  L2 norms after eg_sl_full_wind  **'
    CALL umPrint(umMessage,src='EG_SL_FULL_WIND')
    CALL sl_wind_norm(                                                   &
                      norm_lev_start, norm_lev_end, n_rims_to_do,        &
                      etadot, u_np1, v_np1, w_np1, etadot_np1,           &
                      u, v, w, r_u, r_v, r_w, r_u_d, r_v_d, r_w_d,       &
                      .TRUE., .FALSE., .FALSE. )
  END IF !  L_print_L2norms

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_sl_full_wind
END MODULE eg_sl_full_wind_mod
