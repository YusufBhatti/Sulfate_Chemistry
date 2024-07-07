! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_rho_pseudo_lbflux_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_RHO_PSEUDO_LBFLUX_MOD'

CONTAINS
SUBROUTINE eg_rho_pseudo_lbflux(g_i_pe,                               &
                depart_scheme, depart_order,                          &
                high_order_scheme, monotone_scheme,                   &
                depart_high_order_scheme,depart_monotone_scheme,      &
                first_constant_r_rho_level,                           &
                interp_vertical_search_tol, check_bottom_levels,      &
                l_shallow, l_rk_dps,                                  &
                l_high, l_mono,                                       &
                l_depart_high, l_depart_mono, lam_max_cfl,            &
                etadot, u_np1, v_np1, w_np1,                          &
                etadot_np1, u, v, w,r_rho, r_rho_d,                   &
                error_code, pseudo_lbflux )

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE nlsizes_namelist_mod,  ONLY: global_row_length,                   &
                                 row_length, rows, n_rows, model_levels
USE um_parvars,            ONLY:  halo_i, halo_j, offx, offy, datastart

USE timestep_mod,      ONLY: timestep
USE level_heights_mod, ONLY: eta_theta_levels, eta_rho_levels,        &
                              xi3_at_theta=>r_theta_levels,           &
                              xi3_at_rho=>r_rho_levels,               &
                              xi3_at_u=>r_at_u, xi3_at_v=>r_at_v

USE eg_interpolation_eta_pmf_mod, ONLY: eg_interpolation_eta_pmf
USE horiz_grid_mod
USE ref_pro_mod
USE atm_fields_bounds_mod
USE departure_pts_mod
USE metric_terms_mod
USE ereport_mod, ONLY: ereport
USE Field_Types

USE eg_sl_rho_mod,            ONLY: eg_sl_rho
USE eg_helmholtz_mod,         ONLY: ec_vol
USE eg_lam_domain_kind_mod,   ONLY: lam_domain_kind
USE global_2d_sums_mod,       ONLY: global_2d_sums
USE eg_total_conservation_mod,ONLY: eg_total_mass
USE eg_total_mass_region_mod, ONLY: eg_total_mass_region
USE eg_alpha_mod,             ONLY: alpha_rho

USE model_domain_mod, ONLY: model_type, mt_lam

IMPLICIT NONE
! Description:
!             This routine gives the pseudo lateral boundary flux (PLF)
!             of dry density.
!
! Method: ENDGame formulation version 4.xx
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: dynamics advection
!
! Code description:
! Language: Fortran 95.
! This code is written to UMDP3 standards.

! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_RHO_PSEUDO_LBFLUX'


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
                     ! interpolation scheme used in Depart routine
  depart_monotone_scheme,                                             &
                     ! code choosing monotone
                     ! interpolation scheme used in Depart routine
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
                   ! True if high order interpolation scheme to
                   ! be used in Depart scheme
  l_depart_mono
                   ! True if monotone interpolation scheme to
                   ! be used in Depart scheme

LOGICAL, INTENT(IN) :: l_shallow, l_rk_dps

! wind components

REAL, INTENT(IN) ::                                                   &
  u_np1(-offx:row_length-1+offx,1-offy:rows+offy, model_levels),      &
  v_np1(1-offx:row_length+offx,-offy:n_rows-1+offy, model_levels),    &
  w_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
  u(-halo_i:row_length-1+halo_i,1-halo_j:rows+halo_j,                 &
         model_levels),                                               &
  v(1-halo_i:row_length+halo_i,-halo_j:n_rows-1+halo_j,               &
        model_levels),                                                &
  w(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,                 &
        0:model_levels)

! Timelevel n arrival point quantities

REAL, INTENT(IN) ::                                                   &
  r_rho(1-offx:row_length+offx,1-offy:rows+offy,model_levels),        &
  etadot(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),     &
  etadot_np1(1-offx:row_length+offx,1-offy:rows+offy,                 &
             0:model_levels)

INTEGER :: error_code   ! Non-zero on exit if error detected.

! Timelevel n (rotated) w-departure point quantity

REAL, INTENT(INOUT) ::                                                &
  r_rho_d(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

REAL, INTENT(OUT) :: pseudo_lbflux

! Local variables

INTEGER :: i, j, k, kk

REAL :: total_mass_before
REAL :: total_mass_after
REAL :: mass_factor

INTEGER :: IS, ie, js, je
REAL, ALLOCATABLE :: work_in(:,:,:,:)
REAL, ALLOCATABLE :: work_out(:,:,:,:)
REAL, ALLOCATABLE :: local_sums(:,:)
REAL              :: global_sums(1)
LOGICAL           :: L_mono_lbflux
REAL :: rdxi3(1:model_levels)
REAL :: rdxi2(pdims%j_start:pdims%j_end)
REAL :: rdxi1(pdims%i_start:pdims%i_end)
REAL :: vol
REAL :: one_plus_divu_alpha_np1
REAL :: r_one_plus_divu_alpha_np1
REAL :: deta_term
REAL :: d_xi2_term
REAL :: d_xi1_term

LOGICAL, PARAMETER :: inc_solver_on = .TRUE.

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE( local_sums(tdims%i_start:tdims%i_end, &
                       tdims%j_start:tdims%j_end) )
ALLOCATE( work_in(1-offx:row_length+offx,   &
                    1-offy:rows+offy,       &
                    model_levels,2)         )
ALLOCATE( work_out(1-offx:row_length+offx,   &
                     1-offy:rows+offy,       &
                     model_levels,2)         )
CALL eg_total_mass_region(IS, ie, js, je, L_exclude_rim=.TRUE.)

local_sums = 0.0
work_in    = 0.0
work_out   = 0.0
!$OMP PARALLEL DO DEFAULT(NONE)                           &
!$OMP& PRIVATE(k, j, i)                                   &
!$OMP& SHARED(model_levels, offy, rows, offx, row_length, &
!$OMP&        lam_domain_kind, work_in, r_rho)
DO k = 1, model_levels
  DO j = 1-offy, rows+offy
    DO i = 1-offx ,row_length+offx
      IF (lam_domain_kind(i,j) == 1) THEN
        ! forecast domain close to the boundaries
        work_in(i,j,k,1) = r_rho(i,j,k)
        work_in(i,j,k,2) = 0.0
      ELSE IF (lam_domain_kind(i,j) == 0) THEN
        ! rim
        work_in(i,j,k,1) = 0.0
        work_in(i,j,k,2) = r_rho(i,j,k)
      ELSE IF (lam_domain_kind(i,j) == 2) THEN
        ! forecast domain far enough from the boundaries
        work_in(i,j,k,1) = 0.0
        work_in(i,j,k,2) = 0.0
      END IF
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

! compute initial masses including the rim
total_mass_before = eg_total_mass(work_in(:,:,:,1),L_exclude_rim=.FALSE.)
 
l_mono_lbflux = .FALSE.
DO kk = 1, 2
  CALL eg_sl_rho( g_i_pe,                                               &
                  inc_solver_on, depart_scheme, depart_order,           &
                  high_order_scheme, monotone_scheme,                   &
                  depart_high_order_scheme,depart_monotone_scheme,      &
                  first_constant_r_rho_level,                           &
                  interp_vertical_search_tol, check_bottom_levels,      &
                  l_shallow, l_rk_dps,                                  &
                  l_high, l_mono_lbflux,                                &
                  l_depart_high, l_depart_mono, lam_max_cfl,            &
                  etadot, u_np1, v_np1, w_np1,                          &
                  etadot_np1, u, v, w,                                  &
                  work_in(:,:,:,kk), work_out(:,:,:,kk),                &
                  error_code)
END DO

DO k = 1, model_levels
  rdxi3(k) = 1.0/( eta_theta_levels(k) - eta_theta_levels(k-1) )
END DO
DO j=pdims%j_start, pdims%j_end
  rdxi2(j) = 1.0/( xi2_v(j) - xi2_v(j-1) )
END DO
DO i=pdims%i_start, pdims%i_end
  rdxi1(i) = 1.0/( xi1_u(i) - xi1_u(i-1) )
END DO

!$OMP PARALLEL DO DEFAULT(NONE)                                    &
!$OMP& PRIVATE(k, j, i, d_xi1_term, d_xi2_term, deta_term, vol,    &
!$OMP&         one_plus_divu_alpha_np1, r_one_plus_divu_alpha_np1) &
!$OMP& SHARED(pdims, model_levels,                                 &
!$OMP&        h2_xi1_u, h3_xi1_u, deta_xi3_u, u_np1, rdxi1,        &
!$OMP&        h1_xi2_v, h3_xi2_v, deta_xi3_v, v_np1, rdxi2,        &
!$OMP&        h1_p_eta, h2_p_eta, h3_p_eta, deta_xi3_theta,        &
!$OMP&        etadot_np1, rdxi3,                                   &
!$OMP&        h1_p, h2_p, h3_p, deta_xi3,                          &
!$OMP&        alpha_rho, timestep,                                 &
!$OMP&        work_out)
DO j=pdims%j_start, pdims%j_end
  k = 1
  DO i=pdims%i_start, pdims%i_end
    d_xi1_term = ( h2_xi1_u(i,j,k)*h3_xi1_u(i,j,k)*               &
                   deta_xi3_u(i,j,k)*u_np1(i,j,k) -               &
                   h2_xi1_u(i-1,j,k)*h3_xi1_u(i-1,j,k)*           &
                   deta_xi3_u(i-1,j,k)*u_np1(i-1,j,k) )*rdxi1(i)

    d_xi2_term = ( h1_xi2_v(i,j,k)*h3_xi2_v(i,j,k)*               &
                   deta_xi3_v(i,j,k)*v_np1(i,j,k) -               &
                   h1_xi2_v(i,j-1,k)*h3_xi2_v(i,j-1,k)*           &
                   deta_xi3_v(i,j-1,k)*v_np1(i,j-1,k) )*rdxi2(j)

    deta_term = ( h1_p_eta(i,j,k)*h2_p_eta(i,j,k)*                &
                    h3_p_eta(i,j,k)*deta_xi3_theta(i,j,k)*        &
                    etadot_np1(i,j,k) )*rdxi3(k)

    vol=h1_p(i,j,k)*h2_p(i,j,k)*h3_p(i,j,k)*deta_xi3(i,j,k)

    one_plus_divu_alpha_np1 =                             &
          1.0                                             &
          + alpha_rho * timestep                          &
            * (d_xi1_term + d_xi2_term + deta_term) / vol
    r_one_plus_divu_alpha_np1 = 1.0 / one_plus_divu_alpha_np1
    work_out(i,j,k,1) = work_out(i,j,k,1) * r_one_plus_divu_alpha_np1
    work_out(i,j,k,2) = work_out(i,j,k,2) * r_one_plus_divu_alpha_np1
  END DO

  DO k = 2, model_levels-1
    DO i=pdims%i_start, pdims%i_end
      d_xi1_term = ( h2_xi1_u(i,j,k)*h3_xi1_u(i,j,k)*               &
                     deta_xi3_u(i,j,k)*u_np1(i,j,k) -               &
                     h2_xi1_u(i-1,j,k)*h3_xi1_u(i-1,j,k)*           &
                     deta_xi3_u(i-1,j,k)*u_np1(i-1,j,k) )*rdxi1(i)

      d_xi2_term = ( h1_xi2_v(i,j,k)*h3_xi2_v(i,j,k)*               &
                     deta_xi3_v(i,j,k)*v_np1(i,j,k) -               &
                     h1_xi2_v(i,j-1,k)*h3_xi2_v(i,j-1,k)*           &
                     deta_xi3_v(i,j-1,k)*v_np1(i,j-1,k) )*rdxi2(j)

      deta_term = ( h1_p_eta(i,j,k)*h2_p_eta(i,j,k)*                &
                  h3_p_eta(i,j,k)*deta_xi3_theta(i,j,k)*            &
                  etadot_np1(i,j,k) -                               &
                     h1_p_eta(i,j,k-1)*h2_p_eta(i,j,k-1)*           &
                     h3_p_eta(i,j,k-1)*deta_xi3_theta(i,j,k-1)*     &
                     etadot_np1(i,j,k-1) )*rdxi3(k)

      vol=h1_p(i,j,k)*h2_p(i,j,k)*h3_p(i,j,k)*deta_xi3(i,j,k)

      one_plus_divu_alpha_np1 =                             &
            1.0                                             &
            + alpha_rho * timestep                          &
              * (d_xi1_term + d_xi2_term + deta_term) / vol
      r_one_plus_divu_alpha_np1 = 1.0 / one_plus_divu_alpha_np1
      work_out(i,j,k,1) = work_out(i,j,k,1) * r_one_plus_divu_alpha_np1
      work_out(i,j,k,2) = work_out(i,j,k,2) * r_one_plus_divu_alpha_np1
    END DO
  END DO

  k = model_levels
  DO i=pdims%i_start, pdims%i_end

    d_xi1_term = ( h2_xi1_u(i,j,k)*h3_xi1_u(i,j,k)*               &
                   deta_xi3_u(i,j,k)*u_np1(i,j,k) -               &
                   h2_xi1_u(i-1,j,k)*h3_xi1_u(i-1,j,k)*           &
                   deta_xi3_u(i-1,j,k)*u_np1(i-1,j,k) )*rdxi1(i)

    d_xi2_term = ( h1_xi2_v(i,j,k)*h3_xi2_v(i,j,k)*               &
                   deta_xi3_v(i,j,k)*v_np1(i,j,k) -               &
                   h1_xi2_v(i,j-1,k)*h3_xi2_v(i,j-1,k)*           &
                   deta_xi3_v(i,j-1,k)*v_np1(i,j-1,k) )*rdxi2(j)

    deta_term = ( - h1_p_eta(i,j,k-1)*h2_p_eta(i,j,k-1)*          &
                  h3_p_eta(i,j,k-1)*deta_xi3_theta(i,j,k-1)*      &
                  etadot_np1(i,j,k-1) )*rdxi3(k)

    vol=h1_p(i,j,k)*h2_p(i,j,k)*h3_p(i,j,k)*deta_xi3(i,j,k)

    one_plus_divu_alpha_np1 =                             &
          1.0                                             &
          + alpha_rho * timestep                          &
            * (d_xi1_term + d_xi2_term + deta_term) / vol
    r_one_plus_divu_alpha_np1 = 1.0 / one_plus_divu_alpha_np1
    work_out(i,j,k,1) = work_out(i,j,k,1) * r_one_plus_divu_alpha_np1
    work_out(i,j,k,2) = work_out(i,j,k,2) * r_one_plus_divu_alpha_np1
  END DO
END DO
!$OMP END PARALLEL DO

! rescale work_out(:,:,:,1) by a factor to conserve the mass at start
 
total_mass_after = eg_total_mass(work_out(:,:,:,1), L_exclude_rim=.FALSE.)  
mass_factor = total_mass_before / total_mass_after

!$OMP PARALLEL DO DEFAULT(NONE)                           &
!$OMP& PRIVATE(k, j, i)                                   &
!$OMP& SHARED(model_levels, offy, rows, offx, row_length, &
!$OMP&        work_out, mass_factor)
DO k = 1, model_levels
  DO j = 1-offy, rows+offy
    DO i = 1-offx ,row_length+offx
      work_out(i,j,k,1) = work_out(i,j,k,1) * mass_factor
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(NONE)                           &
!$OMP& PRIVATE(k, j, i)                                   &
!$OMP& SHARED(js, je, model_levels, is, ie,               &
!$OMP&        local_sums, ec_vol, work_out, work_in)
DO j = js, je
  DO k = 1, model_levels
    DO i = IS, ie
      local_sums(i,j) = local_sums(i,j)                   &
                      + ec_vol(i,j,k) * (                 &
                                 work_out(i,j,k,1)        &
                               + work_out(i,j,k,2)        &
                                - work_in(i,j,k,1)    )
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

CALL global_2d_sums(local_sums, tdims%i_len,   &
                                  tdims%j_len,   &
                                  0, 0, 1, global_sums)

pseudo_lbflux = global_sums(1)

DEALLOCATE( local_sums  )
DEALLOCATE( work_in     )
DEALLOCATE( work_out    )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_rho_pseudo_lbflux
END MODULE eg_rho_pseudo_lbflux_mod
