! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE eg_helm_rhs_inc_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_HELM_RHS_INC_MOD'

CONTAINS
SUBROUTINE eg_helm_rhs_inc(R_u, R_v, R_w, R_p, R_theta, R_rho,           &
                            R_etadot, rho_divu_out,                      &
                            R_m_v_d, R_m_cl_d, R_m_cf_d,                 &
                            R_m_r_d, R_m_gr_d, R_m_cf2_d,                &
                            u, v, w, exner, exner_star, rho,             &
                            thetav, etadot,m_v, m_cl, m_cf, m_r, m_gr,   &
                            m_cf2, R_u_d, R_v_d, R_w_d, R_theta_d,       &
                            R_rho_d, Ih, row_length, rows, n_rows,       &
                            model_levels, inner_it_val,                  &
                            S_u,S_v,S_w,S_thetav,S_m_v,S_m_cl,S_m_cf,    &
                            S_m_cf2,S_m_r,S_m_gr, rhs,psi_w_surf,        &
                            psi_w_lid)

USE timestep_mod,      ONLY: timestep
USE level_heights_mod, ONLY: eta_theta_levels
USE um_parcore,        ONLY: mype, nproc
USE um_parvars,        ONLY: halo_i, halo_j, offx, offy,                 &
                             at_extremity
USE eg_alpha_mod,      ONLY: alpha_w
USE yomhook,           ONLY: lhook, dr_hook
USE parkind1,          ONLY: jprb, jpim
USE eg_vert_damp_mod,  ONLY: mu_w
USE eg_helmholtz_mod
USE eg_sisl_init_mod
USE eg_v_at_poles_mod
USE horiz_grid_mod
USE ref_pro_mod
USE atm_fields_bounds_mod
USE UM_ParParams
USE metric_terms_mod
USE update_moisture_mod
USE eg_dxout_mod
USE helmholtz_const_matrix_mod
USE div_pu_mod
USE coriolis_mod
USE gravity_mod
USE eg_parameters_mod,    ONLY: pole_consts
USE planet_constants_mod, ONLY: p_zero, r
USE nlsizes_namelist_mod, ONLY: global_row_length, global_rows
USE umPrintMgr,           ONLY: umPrint, umMessage, PrintStatus,       &
                                PrStatus_Diag
USE halo_exchange, ONLY: swap_bounds
USE mpp_conf_mod,  ONLY: swap_field_is_vector

USE model_domain_mod, ONLY: model_type, mt_global, mt_lam

IMPLICIT NONE

!
! Description: Code to calculate the star variables for use
!              in the Helmholtz problem and updating of the
!              time level (n+1) fields.
!
! Method: ENDGame formulation version 1.01
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Solver
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!

! In

INTEGER, INTENT(IN) ::  row_length, rows, n_rows, model_levels,         &
                        inner_it_val

REAL ::                                                                 &
psi_w_surf(row_length,rows),                                            &
psi_w_lid (row_length,rows)

! Time step, hydrostatic flag  and off-centre weights

REAL, INTENT(IN) :: Ih

! sources at departure points

REAL, INTENT(IN) ::    R_u_d(udims_s%i_start:udims_s%i_end,     &
                             udims_s%j_start:udims_s%j_end,     &
                             udims_s%k_start:udims_s%k_end)

REAL, INTENT(IN) ::    R_v_d(vdims_s%i_start:vdims_s%i_end,     &
                             vdims_s%j_start:vdims_s%j_end,     &
                             vdims_s%k_start:vdims_s%k_end)

REAL, INTENT(IN) ::      R_w_d(wdims%i_start:wdims%i_end,       &
                               wdims%j_start:wdims%j_end,       &
                               wdims%k_start:wdims%k_end)

REAL, INTENT(IN) ::     R_theta_d(tdims_s%i_start:tdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN) ::  R_rho_d(pdims_s%i_start:tdims_s%i_end,     &
                             pdims_s%j_start:pdims_s%j_end,     &
                             pdims_s%k_start:pdims_s%k_end)


REAL, INTENT(IN) ::   R_m_v_d(tdims%i_start:tdims%i_end,        &
                              tdims%j_start:tdims%j_end,        &
                              tdims%k_start:tdims%k_end)

REAL, INTENT(IN) ::  R_m_cl_d(tdims%i_start:tdims%i_end,        &
                              tdims%j_start:tdims%j_end,        &
                              tdims%k_start:tdims%k_end)

REAL, INTENT(IN) ::  R_m_cf_d(tdims%i_start:tdims%i_end,        &
                              tdims%j_start:tdims%j_end,        &
                              tdims%k_start:tdims%k_end)

REAL, INTENT(IN) ::   R_m_r_d(tdims%i_start:tdims%i_end,        &
                              tdims%j_start:tdims%j_end,        &
                              tdims%k_start:tdims%k_end)

REAL, INTENT(IN) ::  R_m_gr_d(tdims%i_start:tdims%i_end,        &
                              tdims%j_start:tdims%j_end,        &
                              tdims%k_start:tdims%k_end)

REAL, INTENT(IN) :: R_m_cf2_d(tdims%i_start:tdims%i_end,        &
                              tdims%j_start:tdims%j_end,        &
                              tdims%k_start:tdims%k_end)

! Output

REAL :: rhs(pdims_s%i_start:pdims_s%i_end,                     &
            pdims_s%j_start:pdims_s%j_end,                     &
            pdims_s%k_start:pdims_s%k_end)


REAL, INTENT(OUT)   :: R_u(udims_s%i_start:udims_s%i_end,       &
                           udims_s%j_start:udims_s%j_end,       &
                           udims_s%k_start:udims_s%k_end)

REAL, INTENT(OUT)   :: R_v(vdims_s%i_start:vdims_s%i_end,       &
                           vdims_s%j_start:vdims_s%j_end,       &
                           vdims_s%k_start:vdims_s%k_end)

REAL, INTENT(OUT)   ::   R_w(wdims%i_start:wdims%i_end,         &
                             wdims%j_start:wdims%j_end,         &
                             wdims%k_start:wdims%k_end)

REAL, INTENT(OUT)   ::  R_theta(tdims_s%i_start:wdims_s%i_end,  &
                                tdims_s%j_start:tdims_s%j_end,  &
                                tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(OUT)   :: R_rho(pdims_s%i_start:wdims_s%i_end,     &
                             pdims_s%j_start:pdims_s%j_end,     &
                             pdims_s%k_start:pdims_s%k_end)


REAL, INTENT(OUT)   ::  R_p(pdims_s%i_start:wdims_s%i_end,      &
                            pdims_s%j_start:pdims_s%j_end,      &
                            pdims_s%k_start:pdims_s%k_end)

REAL, INTENT(OUT)   :: R_etadot(wdims%i_start:wdims%i_end,      &
                                wdims%j_start:wdims%j_end,      &
                                wdims%k_start:wdims%k_end)

REAL, INTENT(INOUT) ::                                          &
               rho_divu_out(pdims_s%i_start:wdims_s%i_end,      &
                            pdims_s%j_start:pdims_s%j_end,      &
                            pdims_s%k_start:pdims_s%k_end)


! time level (n+1) fields

REAL, INTENT(INOUT) ::     u(udims_s%i_start:udims_s%i_end,       &
                           udims_s%j_start:udims_s%j_end,         &
                           udims_s%k_start:udims_s%k_end)

REAL, INTENT(INOUT) ::     v(vdims_s%i_start:vdims_s%i_end,       &
                           vdims_s%j_start:vdims_s%j_end,         &
                           vdims_s%k_start:vdims_s%k_end)

REAL, INTENT(INOUT) ::     w(wdims_s%i_start:wdims_s%i_end,       &
                           wdims_s%j_start:wdims_s%j_end,         &
                           wdims_s%k_start:wdims_s%k_end)

REAL, INTENT(INOUT) ::     exner(pdims_s%i_start:pdims_s%i_end,   &
                               pdims_s%j_start:pdims_s%j_end,     &
                               pdims_s%k_start:pdims_s%k_end+1)

REAL, INTENT(INOUT) :: exner_star(pdims_s%i_start:pdims_s%i_end,  &
                                pdims_s%j_start:pdims_s%j_end)

REAL, INTENT(INOUT) ::     rho(pdims_s%i_start:pdims_s%i_end,     &
                             pdims_s%j_start:pdims_s%j_end,       &
                             pdims_s%k_start:pdims_s%k_end)

REAL, INTENT(INOUT) ::     thetav(tdims_s%i_start:tdims_s%i_end,  &
                                tdims_s%j_start:tdims_s%j_end,    &
                                tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(INOUT) ::      m_v(tdims_s%i_start:tdims_s%i_end,  &
                                tdims_s%j_start:tdims_s%j_end,  &
                                tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(INOUT) ::     m_cl(tdims_s%i_start:tdims_s%i_end,  &
                                tdims_s%j_start:tdims_s%j_end,  &
                                tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(INOUT) ::     m_cf(tdims_s%i_start:tdims_s%i_end,  &
                                tdims_s%j_start:tdims_s%j_end,  &
                                tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(INOUT) ::      m_r(tdims_s%i_start:tdims_s%i_end,  &
                                tdims_s%j_start:tdims_s%j_end,  &
                                tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(INOUT) ::     m_gr(tdims_s%i_start:tdims_s%i_end,  &
                                tdims_s%j_start:tdims_s%j_end,  &
                                tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(INOUT) ::    m_cf2(tdims_s%i_start:tdims_s%i_end,  &
                                tdims_s%j_start:tdims_s%j_end,  &
                                tdims_s%k_start:tdims_s%k_end)


REAL, INTENT(INOUT) ::   etadot(wdims_s%i_start:wdims_s%i_end,  &
                                wdims_s%j_start:wdims_s%j_end,  &
                                wdims_s%k_start:wdims_s%k_end)

!     Fast Physics source terms
REAL                :: S_u(udims_s%i_start:udims_s%i_end,       &
                           udims_s%j_start:udims_s%j_end,       &
                           udims_s%k_start:udims_s%k_end)

REAL                :: S_v(vdims_s%i_start:vdims_s%i_end,       &
                           vdims_s%j_start:vdims_s%j_end,       &
                           vdims_s%k_start:vdims_s%k_end)

REAL                ::   S_w(wdims%i_start:wdims%i_end,         &
                             wdims%j_start:wdims%j_end,         &
                             wdims%k_start:wdims%k_end)

REAL                :: S_thetav(tdims_s%i_start:wdims_s%i_end,  &
                                tdims_s%j_start:tdims_s%j_end,  &
                                tdims_s%k_start:tdims_s%k_end)

REAL              ::    S_m_v(tdims%i_start:tdims%i_end,        &
                              tdims%j_start:tdims%j_end,        &
                              tdims%k_start:tdims%k_end)

REAL              ::   S_m_cl(tdims%i_start:tdims%i_end,        &
                              tdims%j_start:tdims%j_end,        &
                              tdims%k_start:tdims%k_end)

REAL              ::   S_m_cf(tdims%i_start:tdims%i_end,        &
                              tdims%j_start:tdims%j_end,        &
                              tdims%k_start:tdims%k_end)

REAL              ::    S_m_r(tdims%i_start:tdims%i_end,        &
                              tdims%j_start:tdims%j_end,        &
                              tdims%k_start:tdims%k_end)

REAL              ::   S_m_gr(tdims%i_start:tdims%i_end,        &
                              tdims%j_start:tdims%j_end,        &
                              tdims%k_start:tdims%k_end)

REAL              ::  S_m_cf2(tdims%i_start:tdims%i_end,        &
                              tdims%j_start:tdims%j_end,        &
                              tdims%k_start:tdims%k_end)

! Local variables

INTEGER :: i, j, k
LOGICAL :: L_call_from_solver
LOGICAL :: l_inc_solver_in
LOGICAL :: l_call_from_f1sp
REAL    :: t1, t3, t4, d1, a_p

REAL    :: rdxi1, rdxi2, rdxi3
REAL    :: u_at_w, v_at_w

REAL    :: R_p_err,R_rho_err,R_etadot_err,R_w_err,R_theta_err,  &
           R_u_err, R_v_err

INTEGER :: istat


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_HELM_RHS_INC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!
! Add the fast physics source terms
!

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)  &
!$OMP&            SHARED(model_levels,udims,R_u, S_u,R_u_d,u,    &
!$OMP&            vdims,R_v,S_v,R_v_d,v,pdims,R_w,S_w,R_w_d,     &
!$OMP&            Ih,mu_w,w)
DO k = 1, model_levels
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      R_u(i,j,k) = S_u(i,j,k) + R_u_d(i,j,k) - 2.0*u(i,j,k)
    END DO
  END DO

  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      R_v(i,j,k) = S_v(i,j,k) + R_v_d(i,j,k) - 2.0*v(i,j,k)
    END DO
  END DO

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      R_w(i,j,k) = S_w(i,j,k) + R_w_d(i,j,k)                    &
                          - (2.0*Ih + mu_w(i,j,k))*w(i,j,k)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

IF ( inner_it_val == 1 ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP&            SHARED( model_levels, pdims,R_theta,          &
!$OMP&            S_thetav,R_theta_d,thetav)
  DO k = 0, model_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        R_theta (i,j,k) = S_thetav(i,j,k) + R_theta_d(i,j,k)    &
                           - thetav(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(k) &
!$OMP&            SHARED(R_theta,model_levels)
  DO k = 0, model_levels
    R_theta (:,:,k) = 0.0
  END DO
!$OMP END PARALLEL DO
END IF

IF ( inner_it_val == 1 ) THEN

  CALL update_moisture(m_v, m_cl, m_cf, m_r,m_gr, m_cf2,        &
                       R_m_v_d, R_m_cl_d, R_m_cf_d,             &
                       R_m_r_d, R_m_gr_d, R_m_cf2_d,            &
                       S_m_v, S_m_cl, S_m_cf,                   &
                       S_m_r, S_m_gr, S_m_cf2)
END IF

! Adjust arrival point values to account for extra terms

!------------------------------------------------------------------------------
! Compute Psi_w_surf to be used in eg_sisl_init to compute bottom R_w
!------------------------------------------------------------------------------
a_p = 1.0/(alpha_w*timestep)
k = 0
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    psi_w_surf(i,j) = ( Ih*(                                           &
                          (h3_p_eta(i,j,0)/h1_p_eta(i,j,0))*           &
                      dxi1_xi3(i,j,0)*(intw_u2p(i,1)*u(i-1,j,1)        &
                                      +intw_u2p(i,2)*u(i,j,1)) +       &
                          (h3_p_eta(i,j,0)/h2_p_eta(i,j,0))*           &
                      dxi2_xi3(i,j,0)*(intw_v2p(j,1)*v(i,j-1,1)        &
                                      +intw_v2p(j,2)*v(i,j,1))) -      &
                                R_w_d(i,j,0) )*a_p

    psi_w_lid(i,j) = -R_w_d(i,j,model_levels)*a_p

    ! This does not need setting but it is accessed so set to zero
    R_w(i,j,0)     = 0.0
  END DO
END DO

L_call_from_solver = .TRUE.
l_inc_solver_in = .TRUE.
l_call_from_f1sp = .FALSE.
CALL eg_sisl_init(                                              &
         row_length, rows, n_rows, model_levels,                &
         l_inc_solver_in,L_call_from_solver,l_call_from_f1sp,   &
         Ih,g_theta,u, v, w,  thetav, rho,                      &
         m_v, m_cl, m_cf, m_r, m_gr, m_cf2,                     &
         exner, exner_star, R_u, R_v, R_w, R_theta, R_rho,      &
         S_m_v, S_m_cl, S_m_cf, S_m_r, S_m_gr, S_m_cf2, etadot, &
         psi_w_surf, psi_w_lid)

! Compute R_pi

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k,T1)  &
!$OMP&            SHARED(model_levels,pdims,p_zero,r,               &
!$OMP&            intw_w2rho,thetav,R_p,rho,exner,Hm_pp,rho_ref_pro)
DO k = 1, model_levels
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      t1 = intw_w2rho(k,1)*thetav(i,j,k)                               &
         + intw_w2rho(k,2)*thetav(i,j,k-1)
      R_p(i,j,k) =-p_zero/(r*t1)*exner(i,j,k)**Hm_pp + rho(i,j,k)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

! Compute R_rho
IF ( inner_it_val == 1 ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)    &
!$OMP&            SHARED(model_levels,pdims,                       &
!$OMP&            rho_divu_out, rho, R_rho, rho_ref_pro, R_rho_d)
  DO k = 1, model_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
! Compute (rho-rho_ref)*div(u) (R_rho contains rho - rho*div(u)
        rho_divu_out(i,j,k) = (rho(i,j,k) - R_rho(i,j,k))       &
                        * (1.0 - rho_ref_pro(i,j,k)/rho(i,j,k))
! first inner loop use R_rho as usual
        R_rho(i,j,k) = R_rho(i,j,k) + R_rho_d(i,j,k)            &
                       - 2.0*rho(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  DO k = 1, model_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        R_rho(i,j,k) = (rho(i,j,k) - R_rho(i,j,k))              &
                     * (1.0 - rho_ref_pro(i,j,k)/rho(i,j,k))    &
                     - rho_divu_out(i,j,k)
      END DO
    END DO
  END DO
END IF

DO k = 1, model_levels
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

! R_w after elimination of theta'

      R_w(i,j,k) = Hlm_Ck(i,j,k)                                       &
                   *(R_w(i,j,k) - HM_b(i,j,k)*R_theta(i,j,k))

      R_theta(i,j,k) = R_theta(i,j,k) - HM_theta(i,j,k)*R_w(i,j,k)

    END DO
  END DO
END DO

! Fix R_v at the poles

IF ( model_type == mt_global ) THEN
  IF ( at_extremity(PSouth) ) THEN
    CALL eg_v_at_poles(R_u,R_v, 1.0, udims%j_start,             &
                       vdims%j_start,                           &
                       udims_s,vdims_s)
  END IF

  IF ( at_extremity(PNorth) ) THEN
    CALL eg_v_at_poles(R_u,R_v,-1.0, udims%j_end, vdims%j_end,  &
                       udims_s,vdims_s)
  END IF
END IF

CALL swap_bounds(R_u, udims%i_len, udims%j_len, udims%k_len,           &
                 udims_s%halo_i, udims_s%halo_j,                       &
                 fld_type_u,swap_field_is_vector)
CALL swap_bounds(R_v, vdims%i_len, vdims%j_len, vdims%k_len,           &
                 vdims_s%halo_i, vdims_s%halo_j,                       &
                 fld_type_v,swap_field_is_vector)

DO k = pdims%k_start, pdims%k_end-1
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

      u_at_w = intw_rho2w(k,1)*( intw_u2p(i,1)*R_u(i-1,j,k+1) +   &
                                 intw_u2p(i,2)*R_u(i,j,k+1) ) +   &
               intw_rho2w(k,2)*( intw_u2p(i,1)*R_u(i-1,j,k)   +   &
                                 intw_u2p(i,2)*R_u(i,j,k) )

      v_at_w = intw_rho2w(k,1)*( intw_v2p(j,1)*R_v(i,j-1,k+1) +   &
                                 intw_v2p(j,2)*R_v(i,j,k+1) ) +   &
               intw_rho2w(k,2)*( intw_v2p(j,1)*R_v(i,j-1,k) +     &
                                 intw_v2p(j,2)*R_v(i,j,k) )

      R_etadot(i,j,k) = ( R_w(i,j,k)/h3_p_eta(i,j,k) -              &
                         u_at_w*dxi1_xi3(i,j,k)/                    &
                                       h1_p_eta(i,j,k) -            &
                         v_at_w*dxi2_xi3(i,j,k)/                    &
                                       h2_p_eta(i,j,k) ) /          &
                                         deta_xi3_theta(i,j,k)

    END DO
  END DO
END DO

R_etadot(:,:,pdims%k_end) = 0.0
k = 0
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
     R_w(i,j,0) = -w(i,j,0)

     u_at_w = intw_u2p(i,1)*R_u(i-1,j,k+1) + intw_u2p(i,2)*R_u(i,j,k+1)

     v_at_w = intw_v2p(j,1)*R_v(i,j-1,k+1) + intw_v2p(j,2)*R_v(i,j,k+1)

     R_etadot(i,j,k) = ( R_w(i,j,k)/h3_p_eta(i,j,k) -              &
                         u_at_w*dxi1_xi3(i,j,k)/                    &
                                       h1_p_eta(i,j,k) -            &
                         v_at_w*dxi2_xi3(i,j,k)/                    &
                                       h2_p_eta(i,j,k) ) /          &
                                         deta_xi3_theta(i,j,k)

  END DO
END DO
rhs = 0.0
CALL div_Pu(RHS, R_u,R_v, R_etadot)

! Now build RHS
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP&            PRIVATE(i,j,k,D1)                                    &
!$OMP&            SHARED(model_levels, eta_theta_levels, pdims,        &
!$OMP&                   rhs, R_rho, R_theta, R_p, Hlm_Lp,             &
!$OMP&                   rho_ref_pro, intw_w2rho, thetav_ref_pro)
DO k = 1, model_levels
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

      r_rho(i,j,k) = R_rho(i,j,k) - rhs(i,j,k)

      d1 = rho_ref_pro(i,j,k)                                          &
          /( intw_w2rho(k,1)*thetav_ref_pro(i,j,k)                     &
            +intw_w2rho(k,2)*thetav_ref_pro(i,j,k-1) )
       d1 = d1*( intw_w2rho(k,1)*R_theta(i,j,k)                        &
                +intw_w2rho(k,2)*R_theta(i,j,k-1) )

      rhs(i,j,k) = r_rho(i,j,k) + R_p(i,j,k) + d1

      rhs(i,j,k) = rhs(i,j,k)*Hlm_Lp(i,j,k)

    END DO
  END DO
END DO
!$OMP END PARALLEL DO

! Exner fixed on boundary (exner prime =0)

IF ( model_type == mt_lam ) THEN
  i = pdims%i_end - 1         ! one less p-point in LAM's
  IF ( at_extremity(PWest) ) THEN
    rhs(1,:,:) = 0.0
  END IF
  IF ( at_extremity(PEast) ) THEN
    rhs(i,:,:)   = 0.0
    rhs(i+1,:,:) = 0.0
  END IF

  j = pdims%j_end
  IF ( at_extremity(PSouth) ) THEN
    rhs(:,1,:) = 0.0
  END IF
  IF ( at_extremity(PNorth) ) THEN
    rhs(:,j,:) = 0.0
  END IF
END IF

! Output initial errors
IF ( PrintStatus >= PrStatus_Diag) THEN
  R_p_err      = 0.0
  R_rho_err    = 0.0
  R_etadot_err = 0.0
  R_w_err      = 0.0
  R_theta_err  = 0.0
  R_u_err      = 0.0
  R_v_err      = 0.0

  DO k = 1, model_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        R_p_err      = R_p_err      + R_p(i,j,k)*R_p(i,j,k)
        R_rho_err    = R_rho_err    + R_rho(i,j,k)*R_rho(i,j,k)
        R_etadot_err = R_etadot_err + R_etadot(i,j,k)*R_etadot(i,j,k)
        R_w_err      = R_w_err      + R_w(i,j,k)*R_w(i,j,k)
        R_theta_err  = R_theta_err  + R_theta(i,j,k)*R_theta(i,j,k)
      END DO
    END DO
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        R_u_err = R_u_err + R_u(i,j,k)*R_u(i,j,k)
      END DO
    END DO

    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        R_v_err = R_v_err + R_v(i,j,k)*R_v(i,j,k)
      END DO
    END DO
  END DO
  k = 0
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      R_etadot_err = R_etadot_err + R_etadot(i,j,k)*R_etadot(i,j,k)
      R_theta_err  = R_theta_err  + R_theta(i,j,k)*R_theta(i,j,k)
    END DO
  END DO

  CALL gc_rmax(1,nproc,istat,R_p_err)
  CALL gc_rmax(1,nproc,istat,R_rho_err)
  CALL gc_rmax(1,nproc,istat,R_etadot_err)
  CALL gc_rmax(1,nproc,istat,R_w_err)
  CALL gc_rmax(1,nproc,istat,R_theta_err)
  CALL gc_rmax(1,nproc,istat,R_u_err)
  CALL gc_rmax(1,nproc,istat,R_v_err)


  T1 = 1.0/REAL(global_row_length*global_rows*model_levels)
  R_p_err      = SQRT(R_p_err*T1)
  R_rho_err    = SQRT(R_rho_err*T1)
  R_u_err      = SQRT(R_u_err*T1)
  R_v_err      = SQRT(R_v_err*T1)

  T1 = 1.0/REAL(global_row_length*global_rows*(model_levels+1))
  R_etadot_err = SQRT(R_etadot_err*T1)
  R_w_err      = SQRT(R_w_err*T1)
  R_theta_err  = SQRT(R_theta_err*T1)

  IF ( mype == 0 ) THEN
    CALL umPrint(                                          &
          '_________________________________________________')
    WRITE(umMessage,'(A,E25.5)') '    Error in R_u      ',          &
                             R_u_err
    CALL umPrint(umMessage,src='eg_helm_rhs_inc')
    WRITE(umMessage,'(A,E25.5)') '    Error in R_v      ',          &
                             R_v_err
    CALL umPrint(umMessage,src='eg_helm_rhs_inc')
    WRITE(umMessage,'(A,E25.5)') '    Error in R_w      ',          &
                             R_w_err
    CALL umPrint(umMessage,src='eg_helm_rhs_inc')
    WRITE(umMessage,'(A,E25.5)') '    Error in R_theta  ',          &
                             R_theta_err
    CALL umPrint(umMessage,src='eg_helm_rhs_inc')
    WRITE(umMessage,'(A,E25.5)') '    Error in R_etadot ',          &
                             R_etadot_err
    CALL umPrint(umMessage,src='eg_helm_rhs_inc')
    WRITE(umMessage,'(A,E25.5)') '    Error in R_rho    ',          &
                             R_rho_err
    CALL umPrint(umMessage,src='eg_helm_rhs_inc')
    WRITE(umMessage,'(A,E25.5)') '    Error in R_p      ',          &
                             R_p_err
    CALL umPrint(umMessage,src='eg_helm_rhs_inc')
    CALL umPrint(                                                   &
          '_________________________________________________')

  END IF
END IF


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_helm_rhs_inc
END MODULE eg_helm_rhs_inc_mod
