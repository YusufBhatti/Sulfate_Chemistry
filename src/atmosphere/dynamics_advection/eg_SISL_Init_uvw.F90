! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_sisl_init_uvw_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_SISL_INIT_UVW_MOD'

CONTAINS
SUBROUTINE eg_sisl_init_uvw(                                          &
         row_length, rows, n_rows, model_levels, ih, g_theta,         &
         u, v, w, etadot, thetav, rho, m_v, m_cl, m_cf, m_r,          &
         m_gr, m_cf2, exner, exner_star, r_u, r_v, r_w,               &
         l_call_from_solver, l_call_from_f1sp,psi_w_surf,psi_w_lid)

USE level_heights_mod, ONLY: eta_theta_levels, eta_rho_levels,        &
                              xi3_at_theta=>r_theta_levels
USE um_parvars,        ONLY: offx, offy, halo_i, halo_j, at_extremity
USE um_parparams,      ONLY: Pnorth, Psouth
USE eg_alpha_mod,      ONLY: alpha_u, alpha_v, alpha_w
USE timestep_mod,      ONLY: timestep

USE metric_terms_mod,  ONLY: h1_xi1_u, h1_xi2_v, h1_p_eta, h2_xi1_u, &
                              h2_xi2_v, h2_p_eta, h3_xi1_u, h3_xi2_v, &
                              h3_p_eta, deta_xi3, deta_xi3_theta,     &
                              deta_xi3_u,  deta_xi3_v
USE coriolis_mod,      ONLY: f1_star, f2_star, f3_star,              &
                              f1_comp, f2_comp, f3_comp
USE parkind1,          ONLY: jpim, jprb       !DrHook
USE yomhook,           ONLY: lhook, dr_hook   !DrHook
USE dynamics_input_mod, ONLY: l_simple_coriolis, l_viscosity,        &
                              horiz_viscosity, vert_viscosity
USE planet_constants_mod, ONLY: cp
USE horiz_grid_mod
USE atm_fields_bounds_mod
USE Field_Types
USE non_blocking_halo_exchange, ONLY: &
  begin_swap_bounds, end_swap_bounds
USE mpp_conf_mod,  ONLY: swap_field_is_vector
USE calc_vector_Laplacian_mod, ONLY: calc_vector_Laplacian

USE model_domain_mod, ONLY: model_type, mt_cyclic_lam, mt_lam
USE dynamics_testing_mod, ONLY: l_dry

IMPLICIT NONE
!
! Description: computes time level n arrival point quantities:
!              Ru, Rv, Rw, Rtheta, Rrho, Rm_v, Rm_cl, Rm_cf
!
!
! Method: ENDGame formulation version 1.01,
!         section 11 (solution procedure), paragraph 1.

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

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_SISL_INIT_UVW'

INTEGER, INTENT(IN) :: row_length, rows, n_rows, model_levels

REAL ::                                                                 &
  psi_w_surf(row_length,rows),                                          &
  psi_w_lid (row_length,rows)


! Loop index bounds for arrays defined on p, u, v points respectively

! SI time weights     & hydrostatic switch

REAL, INTENT(IN) :: ih

! Gravity arrays

REAL, INTENT(IN) ::                                                   &
  g_theta (1-offx:row_length+offx, 1-offy:rows+offy,                  &
           0:model_levels)

REAL, INTENT(IN) ::                                                   &
  u(-offx:row_length-1+offx,1-offy:rows+offy,model_levels),           &
  v(1-offx:row_length+offx,-offy:n_rows-1+offy,model_levels),         &
  w(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),          &
  etadot(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),     &
  exner(1-offx:row_length+offx,1-offy:rows+offy,model_levels),        &
  exner_star(1-offx:row_length+offx,1-offy:rows+offy),                &
  rho(1-offx:row_length+offx,1-offy:rows+offy,model_levels ),         &
  m_v(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels ),       &
  m_cl(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels ),      &
  m_cf(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels ),      &
  m_r(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels ),       &
  m_gr(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels ),      &
  m_cf2(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels )


LOGICAL, INTENT(IN) :: l_call_from_solver, l_call_from_f1sp

! Timelevel n arrival point quantities

REAL, INTENT(INOUT) ::                                                &
  r_u(-offx:row_length-1+offx,1-offy:rows+offy,model_levels),         &
  r_v(1-offx:row_length+offx,-offy:n_rows-1+offy,model_levels),       &
  r_w(row_length,rows,0:model_levels)

REAL, INTENT(IN)    ::                                                &
  thetav(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)

! Local variables

INTEGER :: i,j,k, kp1

REAL :: beta_u_dt, beta_v_dt, beta_w_dt
REAL :: p_grad_coeff, thetav_ave, rho_bar

! Allocate temps to aid loop vectorization and readability

REAL ::                                                               &
  rho_wet(1-offx:row_length+offx,1-offy:rows+offy,model_levels),      &
  theta_wet(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),  &
  mix_fact(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
  ustar(-offx:row_length-1+offx,1-offy:rows+offy,model_levels),       &
  vstar(1-offx:row_length+offx,-offy:n_rows-1+offy,model_levels),     &
  wstar(row_length,rows,0:model_levels),                              &
  work1(1-offx:row_length+offx,1-offy:rows+offy,model_levels),        &
  worku(1-offx:row_length+offx,1-offy:rows+offy,model_levels),        &
  workv(1-offx:row_length+offx,1-offy:rows+offy,model_levels),        &
  workw(1-offx:row_length+offx,1-offy:rows+offy,model_levels),        &
  work2(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
  dxi1_u(udims%i_start:udims%i_end),                                  &
  dxi1_p(pdims%i_start:pdims%i_end),                                  &
  dxi2_v(vdims%j_start:vdims%j_end),                                  &
  dxi2_p(pdims%j_start:pdims%j_end),                                  &
  deta_w(0:model_levels),                                             &
  deta_rho(model_levels)

REAL :: visco_u(udims%i_start:udims%i_end,                            &
                udims%j_start:udims%j_end,                            &
                udims%k_start:udims%k_end)

REAL :: visco_v(vdims%i_start:vdims%i_end,                            &
                vdims%j_start:vdims%j_end,                            &
                vdims%k_start:vdims%k_end)

REAL :: visco_w(wdims%i_start:wdims%i_end,                            &
                wdims%j_start:wdims%j_end,                            &
                wdims%k_start:wdims%k_end)


! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


! Calculate total rho and virtual potential temperature at all
! grid points including haloes - saves the fill!

IF (l_dry) THEN
  DO k=tdims_s%k_start, tdims_s%k_end
    DO j=tdims_s%j_start, tdims_s%j_end
      DO i=tdims_s%i_start, tdims_s%i_end
        theta_wet(i,j,k) = thetav(i,j,k)
      END DO
    END DO
  END DO
  DO k=pdims_s%k_start, pdims_s%k_end
    DO j=pdims_s%j_start, pdims_s%j_end
      DO i=pdims_s%i_start, pdims_s%i_end
        rho_wet(i,j,k) =rho(i,j,k)
      END DO
    END DO
  END DO
ELSE
!$OMP  PARALLEL DEFAULT(NONE)                                          &
!$OMP& SHARED( tdims_s, mix_fact, m_v, m_cf, m_cl, m_r, m_gr, m_cf2,   &
!$OMP&         theta_wet, thetav, pdims_s, rho_wet, rho, intw_w2rho,   &
!$OMP&         model_levels ) PRIVATE( i, j, k )
!$OMP  DO SCHEDULE(STATIC)
  DO k = 0, model_levels
    DO j = tdims_s%j_start, tdims_s%j_end
      DO i = tdims_s%i_start, tdims_s%i_end
        mix_fact(i,j,k) = 1.0 + m_v(i,j,k) + m_cf(i,j,k) + m_cl(i,j,k) &
                            + m_r(i,j,k) + m_gr(i,j,k) + m_cf2(i,j,k)
        theta_wet(i,j,k) = thetav(i,j,k)/mix_fact(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
  DO k = pdims_s%k_start, pdims_s%k_end
    DO j = pdims_s%j_start, pdims_s%j_end
      DO i = pdims_s%i_start, pdims_s%i_end
        rho_wet(i,j,k) = rho(i,j,k)*( intw_w2rho(k,1)*mix_fact(i,j,k)  &
                                   +intw_w2rho(k,2)*mix_fact(i,j,k-1))
      END DO
    END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL

END IF

IF (model_type == mt_lam) THEN
  CALL fill_external_halos(rho_wet, pdims%i_len, pdims%j_len,          &
                   model_levels, pdims_s%halo_i, pdims_s%halo_j)
  CALL fill_external_halos(theta_wet, tdims%i_len, tdims%j_len,        &
                   model_levels+1, tdims_s%halo_i, tdims_s%halo_j)
END IF

IF (l_call_from_solver .OR. l_call_from_f1sp) THEN
  beta_u_dt = alpha_u*timestep
  beta_v_dt = alpha_v*timestep
  beta_w_dt = alpha_w*timestep
ELSE
  beta_u_dt = (1.0-alpha_u)*timestep
  beta_v_dt = (1.0-alpha_v)*timestep
  beta_w_dt = (1.0-alpha_w)*timestep
END IF

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

k = 0
deta_w(k) = eta_rho_levels(k+1) - eta_theta_levels(k)

DO k=1, model_levels-1
  deta_w(k) = eta_rho_levels(k+1) - eta_rho_levels(k)
  deta_rho(k)   = eta_theta_levels(k) - eta_theta_levels(k-1)
END DO

k = model_levels
deta_rho(k) = eta_theta_levels(k) - eta_theta_levels(k-1)

IF ( .NOT. l_simple_coriolis ) THEN
  !----------------------------------------------------------------------
  ! Compute u-Coriolis acceleration term (Omega x U)u eqns:
  !   (7.10) ENDGame formulation vn1.01.
  !----------------------------------------------------------------------
!$OMP  PARALLEL DEFAULT(NONE)                                         &
!$OMP& SHARED( model_levels, vdims, vstar, dxi1_p, deta_rho,          &
!$OMP&         h1_xi2_v, h3_xi2_v, deta_xi3_v, intw_p2v, rho_wet, v,  &
!$OMP&         udims, ustar, dxi2_p, h2_xi1_u, h3_xi1_u, deta_xi3_u,  &
!$OMP&         intw_p2u, u )  PRIVATE(i,j,k)
!$OMP  DO SCHEDULE(STATIC)
  DO k=1, model_levels
    DO j=vdims%j_start, vdims%j_end
      DO i=vdims%i_start, vdims%i_end
        vstar(i,j,k) = dxi1_p(i)*deta_rho(k)*h1_xi2_v(i,j,k)*           &
                       h3_xi2_v(i,j,k)*deta_xi3_v(i,j,k)*               &
                       ( intw_p2v(j,1)*rho_wet(i,j,k) +                 &
                       intw_p2v(j,2)*rho_wet(i,j+1,k) )*v(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

  !----------------------------------------------------------------------
  ! Compute v-Coriolis acceleration term (Omega x U)v eqns:
  !   (7.11) ENDGame formulation vn1.01.
  !----------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
  DO k=1, model_levels
    DO j=udims%j_start, udims%j_end
      DO i=udims%i_start, udims%i_end
        ustar(i,j,k) = dxi2_p(j)*deta_rho(k)*h2_xi1_u(i,j,k)*           &
                       h3_xi1_u(i,j,k)*deta_xi3_u(i,j,k)*               &
                       ( intw_p2u(i,1)*rho_wet(i,j,k) +                 &
                         intw_p2u(i,2)*rho_wet(i+1,j,k) )*u(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL

  CALL begin_swap_bounds(vstar,                                          &
                   vdims_s%i_len - 2*vdims_s%halo_i,                     &
                   vdims_s%j_len - 2*vdims_s%halo_j,                     &
                   vdims_s%k_len,                                        &
                   vdims_s%halo_i, vdims_s%halo_j,                       &
                   fld_type_v,swap_field_is_vector,                      &
                   do_east_arg=.FALSE.,  do_west_arg=.FALSE.,            &
                   do_south_arg=.FALSE., do_north_arg=.TRUE.,            &
                   do_corners_arg=.FALSE.)


  CALL begin_swap_bounds(ustar,                                          &
                   udims_s%i_len - 2*udims_s%halo_i,                     &
                   udims_s%j_len - 2*udims_s%halo_j,                     &
                   udims_s%k_len,                                        &
                   udims_s%halo_i, udims_s%halo_j,                       &
                   fld_type_u,swap_field_is_vector,                      &
                   do_east_arg=.TRUE.,  do_west_arg=.FALSE.,             &
                   do_south_arg=.FALSE., do_north_arg=.FALSE.,           &
                   do_corners_arg=.FALSE.)


  ! Bottom BC: rho(k=0) = rho(k=1/2)

!$OMP  PARALLEL DEFAULT(NONE)                                          &
!$OMP& SHARED( rows, row_length, pdims, model_levels, dxi1_p, dxi2_p,  &
!$OMP&         wstar, h1_p_eta, h2_p_eta, intw_rho2w, rho_wet, w )     &
!$OMP& PRIVATE(i,j,k)
  k = 0
!$OMP  DO SCHEDULE(STATIC)
  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end
      wstar(i,j,k) = dxi1_p(i)*dxi2_p(j)*h1_p_eta(i,j,k)*             &
                     h2_p_eta(i,j,k)*rho_wet(i,j,k+1)*w(i,j,k)
    END DO
  END DO
!$OMP END DO NOWAIT

  ! Top BC:
  k = model_levels
!$OMP DO SCHEDULE(STATIC)
  DO j=1, rows
    DO i=1, row_length
      wstar(i,j,k) = 0.0
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO k=1, model_levels-1
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        wstar(i,j,k) = dxi1_p(i)*dxi2_p(j)*h1_p_eta(i,j,k)*             &
                 h2_p_eta(i,j,k)*( intw_rho2w(k,1)*rho_wet(i,j,k+1) +   &
                 intw_rho2w(k,2)*rho_wet(i,j,k) )*w(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL
END IF

! Compute worku = (<vstar>^xi2)*f3_star - (<wstar>^eta)*f2_star

IF ( l_simple_coriolis ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)       &
!$OMP& SHARED(model_levels,pdims,v,f3_comp,w,f2_comp,worku,           &
!$OMP&        intw_v2p,intw_w2rho)
  DO k=1, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        worku(i,j,k) = f3_comp(i,j)*                                 &
                         ( intw_v2p(j,1)*v(i,j-1,k)                  &
                          +intw_v2p(j,2)*v(i,j,k) )                  &
                      -f2_comp(i,j)*                                 &
                         ( intw_w2rho(k,1)*w(i,j,k)                  &
                          +intw_w2rho(k,2)*w(i,j,k-1) )
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE

  CALL end_swap_bounds(vstar)

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)       &
!$OMP& SHARED(model_levels,pdims,vstar,f3_star,wstar,f2_star,worku)
  DO k=1, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        worku(i,j,k) = 0.5*(vstar(i,j,k)+vstar(i,j-1,k))*               &
                f3_star(i,j,k) - 0.5*(wstar(i,j,k)+wstar(i,j,k-1))*     &
                f2_star(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

CALL begin_swap_bounds(worku,                                          &
                 pdims_s%i_len - 2*pdims_s%halo_i,                     &
                 pdims_s%j_len - 2*pdims_s%halo_j,                     &
                 pdims_s%k_len,                                        &
                 pdims_s%halo_i, pdims_s%halo_j,                       &
                 fld_type_p,swap_field_is_vector,                      &
                 do_east_arg=.FALSE.,  do_west_arg=.TRUE.,             &
                 do_south_arg=.FALSE., do_north_arg=.FALSE.,           &
                 do_corners_arg=.FALSE.)

! Compute workv = (<wstar>^eta)*f1_star - (<ustar>^xi1)*f3_star

IF ( l_simple_coriolis ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)       &
!$OMP& SHARED(model_levels,pdims,w,f1_comp,u,f3_comp,workv,           &
!$OMP&        intw_w2rho,intw_u2p)
  DO k=1, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        workv(i,j,k) = f1_comp(i,j)*                                   &
                          ( intw_w2rho(k,1)*w(i,j,k)                   &
                           +intw_w2rho(k,2)*w(i,j,k-1) )               &
                      -f3_comp(i,j)*                                   &
                          ( intw_u2p(i,1)*u(i-1,j,k)                   &
                           +intw_u2p(i,2)*u(i,j,k) )
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE

  CALL end_swap_bounds(ustar)

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)       &
!$OMP& SHARED(model_levels,pdims,wstar,f1_star,ustar,f3_star,workv)
  DO k=1, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        workv(i,j,k) = 0.5*(wstar(i,j,k)+wstar(i,j,k-1))*               &
                       f1_star(i,j,k) -                                 &
                       0.5*(ustar(i,j,k)+ustar(i-1,j,k))*               &
                       f3_star(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

CALL begin_swap_bounds(workv,                                          &
                 pdims_s%i_len - 2*pdims_s%halo_i,                     &
                 pdims_s%j_len - 2*pdims_s%halo_j,                     &
                 pdims_s%k_len,                                        &
                 pdims_s%halo_i, pdims_s%halo_j,                       &
                 fld_type_p,swap_field_is_vector,                      &
                 do_east_arg=.FALSE.,  do_west_arg=.FALSE.,            &
                 do_south_arg=.TRUE., do_north_arg=.FALSE.,            &
                 do_corners_arg=.FALSE.)


CALL end_swap_bounds(worku)

IF ( l_simple_coriolis ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)         &
!$OMP& SHARED(model_levels,udims,u,r_u,beta_u_dt,worku,intw_p2u)
  DO k=1, model_levels
    DO j=udims%j_start, udims%j_end
      DO i=udims%i_start, udims%i_end
        r_u(i,j,k) = u(i,j,k) + r_u(i,j,k) + beta_u_dt *                &
                     ( intw_p2u(i,1)*worku(i,j,k)                       &
                      +intw_p2u(i,2)*worku(i+1,j,k) )
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)       &
!$OMP& SHARED(model_levels,udims,u,r_u,beta_u_dt,worku,h1_xi1_u,dxi1_u)
  DO k=1, model_levels
    DO j=udims%j_start, udims%j_end
      DO i=udims%i_start, udims%i_end
        r_u(i,j,k) = u(i,j,k) + r_u(i,j,k) + beta_u_dt *                &
                     0.5*(worku(i,j,k)+worku(i+1,j,k)) /                &
                     ( h1_xi1_u(i,j,k)*dxi1_u(i) )
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF


CALL end_swap_bounds(workv)

IF ( l_simple_coriolis ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)       &
!$OMP& SHARED(model_levels,vdims,r_v,beta_v_dt,workv,v,intw_p2v)
  DO k=1, model_levels
    DO j=vdims%j_start, vdims%j_end
      DO i=vdims%i_start, vdims%i_end
        r_v(i,j,k) = v(i,j,k) + r_v(i,j,k) + beta_v_dt *                &
                     ( intw_p2v(j,1)*workv(i,j,k)                       &
                      +intw_p2v(j,2)*workv(i,j+1,k) )
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)       &
!$OMP& SHARED(model_levels,vdims,r_v,beta_v_dt,h2_xi2_v,dxi2_v,workv,v)
  DO k=1, model_levels
    DO j=vdims%j_start, vdims%j_end
      DO i=vdims%i_start, vdims%i_end
        r_v(i,j,k) = v(i,j,k) + r_v(i,j,k) + beta_v_dt *                &
                     0.5*(workv(i,j,k)+workv(i,j+1,k)) /                &
                     ( h2_xi2_v(i,j,k)*dxi2_v(j) )
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

!----------------------------------------------------------------------
! Compute w-Coriolis acceleration term (Omega x U)w eqns:
!   (7.12) ENDGame formulation vn1.01.
!----------------------------------------------------------------------


! workw = (<ustar>^xi1)*f2_star - (<vstar>^xi2)*f1_star

! approximate ustar, vstar on surface with ustar(1/2), vstar(1/2).
! Coriolis terms are also approximated with Coriolis(1/2).
! This is consistent with the approximation done for rho at the
! surface: rho(surf) ~ rho(1/2)

IF ( l_simple_coriolis ) THEN
  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end
      workw(i,j,1)= f2_comp(i,j)*                                       &
                       ( intw_u2p(i,1)*u(i-1,j,1)                       &
                        +intw_u2p(i,2)*u(i,j,1) )                       &
                   -f1_comp(i,j)*                                       &
                       ( intw_v2p(j,1)*v(i,j-1,1)                       &
                        +intw_v2p(j,2)*v(i,j,1) )
    END DO
  END DO

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k,kp1)     &
!$OMP& SHARED(model_levels,pdims,r_w,beta_w_dt,workw,u,v,w,ih,          &
!$OMP&        f1_comp,f2_comp,intw_u2p,intw_v2p,intw_rho2w)
  DO k=1, model_levels-1
    kp1 = k+1
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        workw(i,j,kp1) = f2_comp(i,j)*                                  &
                           ( intw_u2p(i,1)*u(i-1,j,kp1)                 &
                            +intw_u2p(i,2)*u(i,j,kp1) )                 &
                        -f1_comp(i,j)*                                  &
                           ( intw_v2p(j,1)*v(i,j-1,kp1)                 &
                            +intw_v2p(j,2)*v(i,j,kp1) )
        ! workw(k) already computed
        r_w(i,j,k) = r_w(i,j,k) + ih*w(i,j,k) +                         &
                     beta_w_dt*( intw_rho2w(k,1)*workw(i,j,kp1)         &
                                +intw_rho2w(k,2)*workw(i,j,k))
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
!$OMP  PARALLEL DEFAULT(NONE)                                           &
!$OMP& SHARED( pdims, workw, ustar, f2_star, vstar, f1_star,            &
!$OMP&         model_levels, r_w, ih, w, beta_w_dt, h3_p_eta,           &
!$OMP&         deta_xi3_theta, deta_w ) PRIVATE( i, j, k, kp1 )
!$OMP  DO SCHEDULE(STATIC) 
  DO k=1, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        workw(i,j,k) = 0.5*(ustar(i,j,k)+ustar(i-1,j,k))*              &
                           f2_star(i,j,k) -                            &
                       0.5*(vstar(i,j,k)+vstar(i,j-1,k))*              &
                           f1_star(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
  DO k=1, model_levels-1
    kp1 = k+1
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        ! workw(k) already computed
        r_w(i,j,k) = r_w(i,j,k) + ih*w(i,j,k) +                         &
                     beta_w_dt*0.5*(workw(i,j,kp1) +                    &
                     workw(i,j,k))/(h3_p_eta(i,j,k)*                    &
                     deta_xi3_theta(i,j,k)*deta_w(k))
      END DO
    END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL
END IF

!----------------------------------------------------------------------
! Compute components of pressure gradient in Exner-thetav form:
!   (7.2) ENDGame formulation vn2.02.
!----------------------------------------------------------------------


!----------------------------------------------------------------------
! U Component of Pressure Gradient in R_u Eqn(7.7) of EG2.02
!----------------------------------------------------------------------

! work1: avg(thetav)^(eta) in pressure grad coeff of eqn 7.7.
! work2: avg(exner)^(xi1 eta)*(d(xi3)/d(xi1)) (last term of eqn 7.7)

! *** NB. (1) c_pd=CP=const is assumed and taken outside averaging.
!         (2) It is assumed that horizontal and vertical averages
!             commute, ie. avg(X)^(xi eta) = avg(X)^(eta xi)

! Calculate averaged thetav
!$OMP  PARALLEL DEFAULT(NONE)                                         &
!$OMP& SHARED( model_levels, pdims_s, work1, intw_w2rho, theta_wet,   &
!$OMP&         udims, work2, intw_p2u, exner_star, xi3_at_theta,      &
!$OMP&         dxi1_u, intw_rho2w, exner, cp, h1_xi1_u, deta_xi3_u,   &
!$OMP&         r_u, beta_u_dt, deta_xi3, deta_rho, vdims, intw_p2v,   &
!$OMP&         dxi2_v, r_v, beta_v_dt, r_w, ih, beta_w_dt, psi_w_surf,&
!$OMP&         psi_w_lid, pdims, h3_p_eta, deta_w, g_theta, h2_xi2_v, &
!$OMP&         deta_xi3_v, w, deta_xi3_theta )                        &
!$OMP&         PRIVATE( i, j, k, p_grad_coeff, thetav_ave )
!$OMP  DO SCHEDULE(STATIC)
DO k=1, model_levels
  DO j=pdims_s%j_start, pdims_s%j_end
    DO i=pdims_s%i_start, pdims_s%i_end
      work1(i,j,k) = intw_w2rho(k,1)*theta_wet(i,j,k) +               &
                     intw_w2rho(k,2)*theta_wet(i,j,k-1)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

! Calculate vertical Exner gradient

!$OMP DO SCHEDULE(STATIC)
!CDIR NOUNROLL
DO j=udims%j_start, udims%j_end
  DO i=udims%i_start, udims%i_end
    work2(i,j,0) = (intw_p2u(i,1)*exner_star(i,j) +                   &
                    intw_p2u(i,2)*exner_star(i+1,j)) *                &
                   (xi3_at_theta(i+1,j,0) -                           &
                    xi3_at_theta(i,j,0)) / dxi1_u(i)
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO k=1, model_levels-1
  DO j=udims%j_start, udims%j_end
    DO i=udims%i_start, udims%i_end
      work2(i,j,k) = (intw_rho2w(k,1) *                               &
                      (intw_p2u(i,1)*exner(i,j,k+1) +                 &
                       intw_p2u(i,2)*exner(i+1,j,k+1)) +              &
                      intw_rho2w(k,2) *                               &
                      (intw_p2u(i,1)*exner(i,j,k)+                    &
                       intw_p2u(i,2)*exner(i+1,j,k))) *               &
                     (xi3_at_theta(i+1,j,k) -                         &
                      xi3_at_theta(i,j,k)) / dxi1_u(i)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

k = model_levels
!$OMP DO SCHEDULE(STATIC)
DO j=udims%j_start, udims%j_end
  DO i=udims%i_start, udims%i_end
    work2(i,j,k) = 0.0
  END DO
END DO
!$OMP END DO

! Add pressure gradient term to R_u

!$OMP  DO SCHEDULE(STATIC)
DO k=1, model_levels
  DO j=udims%j_start, udims%j_end
    DO i=udims%i_start, udims%i_end

      p_grad_coeff = cp / ( h1_xi1_u(i,j,k)*deta_xi3_u(i,j,k) )

      thetav_ave = intw_p2u(i,1)*work1(i,j,k) +                       &
                   intw_p2u(i,2)*work1(i+1,j,k)

      r_u(i,j,k) = r_u(i,j,k) -                                       &
                   beta_u_dt * p_grad_coeff * thetav_ave *            &
                   ( (exner(i+1,j,k)*deta_xi3(i+1,j,k) -              &
                      exner(i,j,k)*deta_xi3(i,j,k)) / dxi1_u(i) -     &
                     (work2(i,j,k)-work2(i,j,k-1)) / deta_rho(k) )


    END DO
  END DO
END DO
!$OMP END DO

!----------------------------------------------------------------------
! V Component of Pressure Gradient in R_v Eqn(7.8) of EG2.02
!----------------------------------------------------------------------

! work2: avg(exner)^(xi2 eta)*(d(xi3)/d(xi2)) (last term of eqn 7.8)

! Calculate vertical Exner gradient

k = 0
!$OMP DO SCHEDULE(STATIC)
!CDIR NOUNROLL
DO j=vdims%j_start, vdims%j_end
  DO i=vdims%i_start, vdims%i_end
    work2(i,j,k) = (intw_p2v(j,1)*exner_star(i,j) +                 &
                    intw_p2v(j,2)*exner_star(i,j+1)) *              &
                   (xi3_at_theta(i,j+1,k) -                         &
                    xi3_at_theta(i,j,k)) / dxi2_v(j)
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO k=1, model_levels-1
  DO j=vdims%j_start, vdims%j_end
    DO i=vdims%i_start, vdims%i_end
      work2(i,j,k) = (intw_rho2w(k,1) *                               &
                      (intw_p2v(j,1)*exner(i,j,k+1) +                 &
                       intw_p2v(j,2)*exner(i,j+1,k+1)) +              &
                      intw_rho2w(k,2) *                               &
                      (intw_p2v(j,1)*exner(i,j,k) +                   &
                       intw_p2v(j,2)*exner(i,j+1,k))) *               &
                     (xi3_at_theta(i,j+1,k) -                         &
                      xi3_at_theta(i,j,k)) / dxi2_v(j)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

k = model_levels
!$OMP DO SCHEDULE(STATIC)
DO j=vdims%j_start, vdims%j_end
  DO i=vdims%i_start, vdims%i_end
    work2(i,j,k) = 0.0
  END DO
END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
DO k=1, model_levels
  DO j=vdims%j_start, vdims%j_end
    DO i=vdims%i_start, vdims%i_end

      p_grad_coeff = cp / ( h2_xi2_v(i,j,k)*deta_xi3_v(i,j,k) )

      thetav_ave = intw_p2v(j,1)*work1(i,j,k) +                       &
                   intw_p2v(j,2)*work1(i,j+1,k)

      r_v(i,j,k) = r_v(i,j,k) -                                       &
                   beta_v_dt * p_grad_coeff * thetav_ave *            &
                   ( (exner(i,j+1,k)*deta_xi3(i,j+1,k) -              &
                      exner(i,j,k)*deta_xi3(i,j,k)) / dxi2_v(j) -     &
                     (work2(i,j,k)-work2(i,j,k-1)) / deta_rho(k) )

    END DO
  END DO
END DO
!$OMP END DO NOWAIT

!----------------------------------------------------------------------
! W Component of Pressure Gradient in R_w Eqn(7.9) of EG2.02
!----------------------------------------------------------------------

! vpg term at bottom, k=0

!$OMP DO SCHEDULE(STATIC)
DO j=pdims%j_start, pdims%j_end
  DO i=pdims%i_start, pdims%i_end

    k = 0
    r_w(i,j,k) = r_w(i,j,k) + ih*w(i,j,k) + beta_w_dt*psi_w_surf(i,j)

    k = model_levels
    r_w(i,j,k) = r_w(i,j,k) + ih*w(i,j,k) + beta_w_dt*psi_w_lid(i,j)

  END DO
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO k=1, model_levels-1
  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end

      p_grad_coeff = cp/(h3_p_eta(i,j,k)*deta_xi3_theta(i,j,k))

      r_w(i,j,k) = r_w(i,j,k) -                                       &
                   beta_w_dt * (p_grad_coeff * theta_wet(i,j,k)*      &
                   (exner(i,j,k+1)-exner(i,j,k)) / deta_w(k) +        &
                   g_theta(i,j,k))

    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

! ---------------------------------------------------------------------
!   Compute Components of Viscosity
! ---------------------------------------------------------------------
IF (l_viscosity .AND.                                                   &
    .NOT.(l_call_from_solver .OR. l_call_from_f1sp)) THEN

  beta_u_dt=timestep
  beta_v_dt=timestep
  beta_w_dt=timestep

  CALL calc_vector_Laplacian(u, v, w, etadot, visco_u, visco_v, visco_w)

 IF (horiz_viscosity > 0.0) THEN
!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP& PRIVATE(i,j,k,rho_bar)                                           &
!$OMP& SHARED(udims,rho_wet,R_u,beta_u_dt,horiz_viscosity,visco_u)
    DO k=udims%k_start, udims%k_end
       DO j=udims%j_start, udims%j_end
          DO i=udims%i_start, udims%i_end
             rho_bar=0.5*(rho_wet(i,j,k)+rho_wet(i+1,j,k))
             R_u(i,j,k)=R_u(i,j,k)+beta_u_dt*(horiz_viscosity/rho_bar) *   &
                  visco_u(i,j,k)
          END DO
       END DO
    END DO
!$OMP END PARALLEL DO

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP& PRIVATE(i,j,k,rho_bar)                                           &
!$OMP& SHARED(vdims,rho_wet,R_v,beta_v_dt,horiz_viscosity,visco_v)
    DO k=vdims%k_start, vdims%k_end
       DO j=vdims%j_start, vdims%j_end
          DO i=vdims%i_start, vdims%i_end
             rho_bar=0.5*(rho_wet(i,j,k)+rho_wet(i,j+1,k))
             R_v(i,j,k)=R_v(i,j,k)+beta_v_dt*(horiz_viscosity/rho_bar) *   &
                  visco_v(i,j,k)
          END DO
       END DO
    END DO
!$OMP END PARALLEL DO
 END IF

  IF (vert_viscosity > 0.0) THEN
!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP& PRIVATE(i,j,k,rho_bar)                                           &
!$OMP& SHARED(wdims,rho_wet,R_w,beta_w_dt,vert_viscosity,visco_w)
    DO k=wdims%k_start+1, wdims%k_end-1
      DO j=wdims%j_start, wdims%j_end
        DO i=wdims%i_start, wdims%i_end
          rho_bar=0.5*(rho_wet(i,j,k)+rho_wet(i,j,k+1))
          R_w(i,j,k)=R_w(i,j,k)+beta_w_dt*(vert_viscosity/rho_bar) *   &
                              visco_w(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

END IF

! Zero R_v on bloundary for channel flows

IF (model_type == mt_cyclic_lam) THEN
  IF ( at_extremity(Psouth) ) THEN
    R_v(:,vdims%j_start,:) = 0.0
  END IF
  IF ( at_extremity(Pnorth) ) THEN
    R_v(:,vdims%j_end,:) = 0.0
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_sisl_init_uvw
END MODULE eg_sisl_init_uvw_mod
