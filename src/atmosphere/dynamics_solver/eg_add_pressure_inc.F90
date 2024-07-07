! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE eg_add_pressure_inc_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_ADD_PRESSURE_INC_MOD'

CONTAINS
SUBROUTINE eg_add_pressure_inc(exner_prime, exner_np1, exner_star_np1,   &
           R_v_South, R_v_North,                                         &
           u_np1, v_np1, w_np1, thetav_np1, etadot_np1, rho_np1,         &
           m_v_np1, m_cl_np1, m_cf_np1, m_r_np1, m_gr_np1,               &
           m_cf2_np1, R_w_d, R_u, R_v, R_w, R_theta, R_etadot, R_p,      &
           Ih,n_rows, row_length, rows, model_levels, w_surf, w_lid)

USE um_parvars,          ONLY: halo_i, halo_j, at_extremity
USE level_heights_mod,   ONLY: eta_theta_levels, eta_rho_levels,        &
                               xi3_at_theta=>r_theta_levels,            &
                               xi3_at_rho=>r_rho_levels
USE timestep_mod,        ONLY: timestep
USE eg_alpha_mod,        ONLY: alpha_w
USE planet_constants_mod, ONLY: p_zero, r
USE eg_vert_damp_mod,    ONLY: mu_w
USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim
USE eg_helmholtz_mod
USE eg_coriolis_star_mod
USE eg_v_at_poles_mod
USE eg_calc_p_star_mod

USE horiz_grid_mod
USE ref_pro_mod
USE atm_fields_bounds_mod
USE UM_ParParams
USE metric_terms_mod
USE helmholtz_const_matrix_mod
USE coriolis_mod
USE gravity_mod
USE eg_parameters_mod,   ONLY: pole_consts
USE halo_exchange, ONLY: swap_bounds
USE mpp_conf_mod,  ONLY: swap_field_is_vector, swap_field_is_scalar
USE pressure_grad_mod

USE model_domain_mod, ONLY: model_type, mt_global, mt_cyclic_lam

IMPLICIT NONE

!
! Description:
!         Update values to obtain the timelevel n+1 fields
!
!
! Method: ENDGame formulation version 3.02
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Solver
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!

! Array dimensions

INTEGER,       INTENT(IN)    :: row_length
INTEGER,       INTENT(IN)    :: rows
INTEGER,       INTENT(IN)    :: model_levels
INTEGER,       INTENT(IN)    :: n_rows
REAL,          INTENT(IN)    :: Ih

! Fields at timelevel n+1

REAL, INTENT(INOUT) :: u_np1(udims_s%i_start:udims_s%i_end,     &
                           udims_s%j_start:udims_s%j_end,       &
                           udims_s%k_start:udims_s%k_end)

REAL, INTENT(INOUT) :: v_np1(vdims_s%i_start:vdims_s%i_end,     &
                           vdims_s%j_start:vdims_s%j_end,       &
                           vdims_s%k_start:vdims_s%k_end)

REAL, INTENT(INOUT) :: w_np1(wdims_s%i_start:wdims_s%i_end,     &
                           wdims_s%j_start:wdims_s%j_end,       &
                           wdims_s%k_start:wdims_s%k_end)

REAL, INTENT(INOUT) :: thetav_np1(tdims_s%i_start:tdims_s%i_end,&
                                tdims_s%j_start:tdims_s%j_end,  &
                                tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(INOUT) :: rho_np1(pdims_s%i_start:pdims_s%i_end,   &
                             pdims_s%j_start:pdims_s%j_end,     &
                             pdims_s%k_start:pdims_s%k_end)

REAL, INTENT(INOUT) :: etadot_np1(wdims_s%i_start:wdims_s%i_end,&
                                wdims_s%j_start:wdims_s%j_end,  &
                                wdims_s%k_start:wdims_s%k_end)

REAL, INTENT(INOUT) :: exner_np1(pdims_s%i_start:pdims_s%i_end, &
                               pdims_s%j_start:pdims_s%j_end,   &
                               pdims_s%k_start:pdims_s%k_end+1)

REAL, INTENT(IN) ::     m_v_np1(tdims_s%i_start:tdims_s%i_end,  &
                                tdims_s%j_start:tdims_s%j_end,  &
                                tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN) ::    m_cl_np1(tdims_s%i_start:tdims_s%i_end,  &
                                tdims_s%j_start:tdims_s%j_end,  &
                                tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN) ::    m_cf_np1(tdims_s%i_start:tdims_s%i_end,  &
                                tdims_s%j_start:tdims_s%j_end,  &
                                tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN) ::     m_r_np1(tdims_s%i_start:tdims_s%i_end,  &
                                tdims_s%j_start:tdims_s%j_end,  &
                                tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN) ::    m_gr_np1(tdims_s%i_start:tdims_s%i_end,  &
                                tdims_s%j_start:tdims_s%j_end,  &
                                tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN) ::   m_cf2_np1(tdims_s%i_start:tdims_s%i_end,  &
                                tdims_s%j_start:tdims_s%j_end,  &
                                tdims_s%k_start:tdims_s%k_end)

REAL ::          exner_star_np1(pdims_s%i_start:pdims_s%i_end,  &
                                pdims_s%j_start:pdims_s%j_end)

! RHS terms from time level n

REAL, INTENT(IN)    :: R_w_d(wdims%i_start:wdims%i_end,         &
                             wdims%j_start:wdims%j_end,         &
                             wdims%k_start:wdims%k_end)

REAL,  INTENT(IN)   :: R_v_North(vdims_s%i_start:vdims_s%i_end, &
                                 vdims_s%k_start:vdims_s%k_end)

REAL,  INTENT(IN)   :: R_v_South(vdims_s%i_start:vdims_s%i_end, &
                                 vdims_s%k_start:vdims_s%k_end)

! estimates of RHS terms at time level n+1

REAL, INTENT(INOUT) :: R_u(udims_s%i_start:udims_s%i_end,       &
                           udims_s%j_start:udims_s%j_end,       &
                           udims_s%k_start:udims_s%k_end)

REAL, INTENT(INOUT) :: R_v(vdims_s%i_start:vdims_s%i_end,       &
                           vdims_s%j_start:vdims_s%j_end,       &
                           vdims_s%k_start:vdims_s%k_end)

REAL, INTENT(INOUT) ::   R_w(wdims%i_start:wdims%i_end,         &
                             wdims%j_start:wdims%j_end,         &
                             wdims%k_start:wdims%k_end)

REAL, INTENT(INOUT) ::  R_theta(tdims_s%i_start:wdims_s%i_end,  &
                                tdims_s%j_start:tdims_s%j_end,  &
                                tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(INOUT) ::  R_p(pdims_s%i_start:wdims_s%i_end,      &
                            pdims_s%j_start:pdims_s%j_end,      &
                            pdims_s%k_start:pdims_s%k_end)

REAL, INTENT(INOUT) :: R_etadot(wdims%i_start:wdims%i_end,      &
                                wdims%j_start:wdims%j_end,      &
                                wdims%k_start:wdims%k_end)

! Current estimate of the pressure perturbation

REAL :: exner_prime(pdims_s%i_start:wdims_s%i_end,              &
                    pdims_s%j_start:pdims_s%j_end,              &
                    pdims_s%k_start:pdims_s%k_end)

REAL :: w_surf(wdims%i_start:wdims%i_end,                       &
               wdims%j_start:wdims%j_end)

REAL :: w_lid(wdims%i_start:wdims%i_end,                        &
              wdims%j_start:wdims%j_end)

! Local variables

REAL    :: d
REAL    :: rdxi1, rdxi2, rdxi3
REAL    :: t2
INTEGER :: i, j, k

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_ADD_PRESSURE_INC'

! Arrays holding the pressure gradient

REAL                :: Pw(wdims%i_start:wdims%i_end,            &
                          wdims%j_start:wdims%j_end,            &
                          wdims%k_start:wdims%k_end)
REAL                :: Pu(udims_s%i_start:udims_s%i_end,        &
                          udims_s%j_start:udims_s%j_end,        &
                          udims_s%k_start:udims_s%k_end)

REAL                :: Pv(vdims_s%i_start:vdims_s%i_end,        &
                          vdims_s%j_start:vdims_s%j_end,        &
                          vdims_s%k_start:vdims_s%k_end)

REAL :: field_inc, u_at_w, v_at_w

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL pressure_grad(exner_prime,Pu, Pv, Pw)

! Back Subs to get u and v
!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(i,j,k,field_inc)                         &
!$OMP&             SHARED(model_levels, udims, vdims, pdims,                   &
!$OMP&                    r_u, r_v, R_w, R_theta, HM_theta,                    &
!$OMP&                    u_np1, v_np1, w_np1, thetav_np1,  Pu, Pv, Pw)
!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      field_inc    = R_u(i,j,k) - Pu(i,j,k)
      u_np1(i,j,k) = u_np1(i,j,k) + field_inc
    END DO
  END DO

  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      field_inc    = R_v(i,j,k) - Pv(i,j,k)
      v_np1(i,j,k) = v_np1(i,j,k) + field_inc
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

! Back Subs to get etadot and theta
!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels-1
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

      field_inc         = R_w(i,j,k) - Pw(i,j,k)

      w_np1(i,j,k)      = w_np1(i,j,k) + field_inc

! Store theta' in R_theta

      R_theta(i,j,k)    = R_theta(i,j,k) + Hm_theta(i,j,k)*Pw(i,j,k)

      thetav_np1(i,j,k) = thetav_np1(i,j,k) + R_theta(i,j,k)

    END DO
  END DO
END DO
!$OMP END DO NOWAIT

k = model_levels
!$OMP DO SCHEDULE(STATIC)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    thetav_np1(i,j,0) = R_theta(i,j,0) + thetav_np1(i,j,0)
    thetav_np1(i,j,k) = R_theta(i,j,k) + thetav_np1(i,j,k)
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

IF (model_type == mt_global) THEN
  IF ( at_extremity(PSouth) ) THEN
    CALL eg_v_at_poles(u_np1,v_np1, 1.0, udims%j_start,         &
                       vdims%j_start,                           &
                       udims_s,vdims_s)
  END IF

  IF ( at_extremity(PNorth) ) THEN
    CALL eg_v_at_poles(u_np1,v_np1,-1.0, udims%j_end,           &
                       vdims%j_end,                             &
                       udims_s,vdims_s)
  END IF
END IF

CALL swap_bounds(u_np1, udims%i_len, udims%j_len, udims%k_len,     &
                 udims_s%halo_i, udims_s%halo_j,                   &
                 fld_type_u,swap_field_is_vector)
CALL swap_bounds(v_np1, vdims%i_len, vdims%j_len, vdims%k_len,     &
                 vdims_s%halo_i, vdims_s%halo_j,                   &
                 fld_type_v,swap_field_is_vector)


! Surface value used for BCs
etadot_np1(:,:,0)            = 0.0
etadot_np1(:,:,model_levels) = 0.0
w_np1(:,:,model_levels)      = 0.0

d = 1.0/(alpha_w*timestep)
! Back Subs to get u and v
!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(i,j,k,field_inc,u_at_w,v_at_w)           &
!$OMP&             SHARED(pdims,model_levels,d,Ih,                             &
!$OMP&                    intw_u2p,intw_rho2w, intw_v2p,                       &
!$OMP&                    R_w_d, dxi1_xi3,dxi2_xi3,deta_xi3_theta,             &
!$OMP&                    h3_p_eta,h1_p_eta,h2_p_eta,                          &
!$OMP&                    w_surf,w_lid,u_np1, v_np1, w_np1, etadot_np1)
!$OMP DO SCHEDULE(STATIC)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
      u_at_w = intw_u2p(i,1)*u_np1(i-1,j,1) +                       &
               intw_u2p(i,2)*u_np1(i,j,1)

      v_at_w = intw_v2p(j,1)*v_np1(i,j-1,1) +                       &
               intw_v2p(j,2)*v_np1(i,j,1)

    field_inc = ( u_at_w*dxi1_xi3(i,j,0)/h1_p_eta(i,j,0)            &
                 +v_at_w*dxi2_xi3(i,j,0)/h2_p_eta(i,j,0) )          &
                 *h3_p_eta(i,j,0)


! w back subs not enforcing R_etadot = 0

    w_np1(i,j,0) = field_inc

    w_surf(i,j)   = d*( R_w_d(i,j,0) - Ih*w_np1(i,j,0) )
    w_lid(i,j)    = d*R_w_d(i,j,model_levels)
  END DO
END DO
!$OMP END DO NOWAIT

! update etadot with full equation
!$OMP DO SCHEDULE(STATIC)
DO k = pdims%k_start, pdims%k_end-1
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

      u_at_w = intw_rho2w(k,1)*( intw_u2p(i,1)*u_np1(i-1,j,k+1) +      &
                                 intw_u2p(i,2)*u_np1(i,j,k+1) ) +      &
               intw_rho2w(k,2)*( intw_u2p(i,1)*u_np1(i-1,j,k) +        &
                                 intw_u2p(i,2)*u_np1(i,j,k) )

      v_at_w = intw_rho2w(k,1)*( intw_v2p(j,1)*v_np1(i,j-1,k+1) +      &
                                 intw_v2p(j,2)*v_np1(i,j,k+1) ) +      &
               intw_rho2w(k,2)*( intw_v2p(j,1)*v_np1(i,j-1,k) +        &
                                 intw_v2p(j,2)*v_np1(i,j,k) )

      etadot_np1(i,j,k) = ( w_np1(i,j,k)/h3_p_eta(i,j,k) -             &
                            u_at_w*dxi1_xi3(i,j,k)/h1_p_eta(i,j,k) -   &
                            v_at_w*dxi2_xi3(i,j,k)/h2_p_eta(i,j,k) )   &
                          /deta_xi3_theta(i,j,k)

    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

etadot_np1(:,:,0)           = 0.0
etadot_np1(:,:,pdims%k_end) = 0.0

!------------------------------------------------------------------------------

! Update density using eqn(9.32) of EG2.02 and, finally, pressure.

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k,T2)    &
!$OMP&            SHARED(model_levels, pdims, intw_w2rho, r_theta,    &
!$OMP&                   thetav_ref_pro, rho_np1, rho_ref_pro, hm_pp, &
!$OMP&                   exner_prime, exner_ref_pro, r_p, exner_np1)
DO k = 1, model_levels
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

! R_theta now contains theta'
      t2 =  intw_w2rho(k,1)*thetav_ref_pro(i,j,k)                      &
           +intw_w2rho(k,2)*thetav_ref_pro(i,j,k-1)

      t2 = (intw_w2rho(k,1)*R_theta(i,j,k)                             &
           +intw_w2rho(k,2)*R_theta(i,j,k-1) )/t2

      rho_np1(i,j,k) = rho_np1(i,j,k)                                  &
                           +rho_ref_pro(i,j,k)*(                       &
                   Hm_pp*exner_prime(i,j,k)/exner_ref_pro(i,j,k)       &
                            -t2) - R_p(i,j,k)

! Use exner_prime to update exner and eqn of state to update rho
      exner_np1(i,j,k) = exner_np1(i,j,k) + exner_prime(i,j,k)

    END DO
  END DO
END DO
!$OMP END PARALLEL DO

!----------------------------------------------------------------------------
! Compute surface pressure
!----------------------------------------------------------------------------

CALL EG_Calc_P_star(                                                   &
                  model_levels,  row_length, rows, exner_np1,          &
                  thetav_np1,                                          &
                  m_v_np1, m_cl_np1, m_cf_np1, m_r_np1, m_gr_np1,      &
                  m_cf2_np1, g_theta, exner_star_np1, w_surf, w_lid)

CALL swap_bounds(w_np1, wdims%i_len, wdims%j_len, wdims%k_len,         &
                 wdims_s%halo_i, wdims_s%halo_j,                       &
                 fld_type_p, swap_field_is_scalar)

CALL swap_bounds(etadot_np1, wdims%i_len, wdims%j_len, wdims%k_len,    &
                 wdims_s%halo_i, wdims_s%halo_j,                       &
                 fld_type_p, swap_field_is_scalar)

CALL swap_bounds(thetav_np1,tdims%i_len, tdims%j_len, tdims%k_len,     &
                 tdims_s%halo_i, tdims_s%halo_j,                       &
                 fld_type_p, swap_field_is_scalar)

CALL swap_bounds(rho_np1, pdims%i_len, pdims%j_len, pdims%k_len,       &
                 pdims_s%halo_i, pdims_s%halo_j,                       &
                 fld_type_p, swap_field_is_scalar)
CALL swap_bounds(exner_np1, wdims%i_len, wdims%j_len, wdims%k_len,     &
                 wdims_s%halo_i, wdims_s%halo_j,                       &
                 fld_type_p, swap_field_is_scalar)

! IF cyclic LAM fix boundary pressure to balance v-forcing terms

IF (model_type == mt_cyclic_lam) THEN
  DO k = 1, model_levels
    DO i = pdims%i_start, pdims%i_end

      ! North boundary
      j = vdims%j_end
      d = HM_v(i,j,k)*deta_xi3(i,j,k)/(xi2_p(j+1) - xi2_p(j))
      exner_np1(i,j+1,k) = exner_np1(i,j,k) + R_v_North(i,k)/d

      ! South boundary
      j = vdims%j_start
      d = HM_v(i,j,k)*deta_xi3(i,j,k)/(xi2_p(j+1) - xi2_p(j))
      exner_np1(i,j,k) = exner_np1(i,j+1,k) - R_v_South(i,k)/d

    END DO
  END DO
END IF

!-------------------------------------------------------------------------
! Compute "starred" Coriolis terms at new timelevel
!-------------------------------------------------------------------------

w_surf  = -w_surf
w_lid   = -w_lid
CALL eg_coriolis_star(rho_np1, m_v_np1, m_cl_np1, m_cf_np1,     &
                      m_r_np1, m_gr_np1, m_cf2_np1)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_add_pressure_inc
END MODULE eg_add_pressure_inc_mod
